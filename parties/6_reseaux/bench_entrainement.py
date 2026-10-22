"""Partie 6.3 : entraîner un neurone quaternionique, Broyden quaternionique « à droite ».

Pour une couche y = w·x (w ∈ ℍⁿ ligne de poids, x ∈ ℍⁿ), la dérivée par rapport aux poids est δw ↦ δw·x : une multiplication À
DROITE. Sur N exemples, les résidus rₜ(w) = σ(Σⱼ wⱼ·xₜⱼ + c) − yₜ ont une jacobienne à blocs R(xₜⱼ) (multipliés à gauche par
Dσ). On normalise donc le quotient à droite, ρₜ = Fₜ(w + h)·Fₜ(w)⁻¹, ce qui garde la structure (sonde_hn.py, modèle "D").
Problème carré (N = n exemples, n poids quaternioniques, dimension 4n) avec un professeur w* : le zéro existe (interpolation). Le
cas moindres carrés à N > n (minimum non nul) n'est pas traité : la partie 2 montre que la couleur n'y gagne que si |f| est connu.

Activations : "lin" (σ = identité, problème affine), "radial-faible" (tanh radial, entrées 3 fois plus petites : σ presque
linéaire, défaut à droite 0,05–0,09), "radial" (tanh(‖z‖)·z/‖z‖, défaut 0,25–0,42), "split" (tanh par composante).
Méthodes : sonde-…-quaternionique-à-droite (D), sonde-couleur-mémoire-quaternionique (H, le mauvais côté), sonde-couleur-mémoire
(full), Broyden-sans-couleur, Broyden-mêmes-mesures. 24 bits par composante, budget 600 (2000 pour n = 16).
Usage : python bench_entrainement.py [--smoke] [--reprise]"""
import os, sys, json, time
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "src")); sys.path.insert(0, HERE)
from multiprocessing import Pool
import banc_hn as B
from quat import qmul, Lmat, Rmat, E4, rand_unit
from rnn import act, dact_block

LB = [Lmat(E4[k]) for k in range(4)]; RB = [Rmat(E4[k]) for k in range(4)]


def defaut_cote(J, side):
    """Distance relative de J aux matrices à blocs L(·) (side='L') ou R(·) (side='R')."""
    base = LB if side == "L" else RB; n = J.shape[0] // 4; P = np.zeros_like(J)
    for i in range(n):
        for j in range(n):
            Bk = J[4 * i:4 * i + 4, 4 * j:4 * j + 4]
            P[4 * i:4 * i + 4, 4 * j:4 * j + 4] = sum(np.sum(Bk * L) / 4 * L for L in base)
    return float(np.linalg.norm(J - P) / np.linalg.norm(J))


class Neurone:
    def __init__(self, n, kind="radial", seed=0):
        rng = np.random.default_rng(63000 + 101 * seed + 7 * n)
        sc = 0.3 if kind.endswith("-faible") else 1.0          # « faible » : entrées 3 fois plus petites, σ presque linéaire
        kind = kind.replace("-faible", "")
        self.n, self.kind = n, kind
        self.ws = np.array([rand_unit(rng) for _ in range(n)]) * rng.uniform(0.6, 1.2, (n, 1))      # professeur
        self.X = sc * rng.standard_normal((n, n, 4)) / np.sqrt(2 * n)                               # n exemples × n entrées
        self.c = sc * 0.3 * rng.standard_normal((n, 4))
        self.Y = np.array([self._sig(self._z(self.ws, t)) for t in range(n)])
        self.qs = self.ws.reshape(-1)

    def _z(self, w, t):
        return sum(qmul(w[j], self.X[t, j]) for j in range(self.n)) + self.c[t]

    def _sig(self, z):
        return z if self.kind == "lin" else act(z, self.kind)

    def f(self, q):
        w = np.asarray(q, float).reshape(self.n, 4)
        return np.concatenate([self._sig(self._z(w, t)) - self.Y[t] for t in range(self.n)])

    def J(self, q):
        w = np.asarray(q, float).reshape(self.n, 4); n = self.n; Jm = np.zeros((4 * n, 4 * n))
        for t in range(n):
            D = np.eye(4) if self.kind == "lin" else dact_block(self._z(w, t), self.kind)
            for j in range(n):
                Jm[4 * t:4 * t + 4, 4 * j:4 * j + 4] = D @ Rmat(self.X[t, j])
        return Jm


for _k in ("lin", "radial-faible", "radial", "split"):
    B.PROBLEMES["neurone-" + _k] = (lambda kind: (lambda n, seed: (lambda m: (m.f, m.qs, {}))(Neurone(n, kind, seed))))(_k)
METH = {"sonde-D": ("sonde", dict(modele="D")), "sonde-H": ("sonde", dict(modele="H")), "sonde-full": ("sonde", dict(modele="full")),
        "broyden": ("broyden", {}), "broyden-deg": ("broyden-deg", {})}
RAD = 0.5

if __name__ == "__main__":
    smoke = "--smoke" in sys.argv
    NS = [2] if smoke else [2, 4, 8, 16]
    KINDS = ("lin", "radial-faible", "radial", "split")
    if "--kinds" in sys.argv: KINDS = tuple(sys.argv[sys.argv.index("--kinds") + 1].split(","))
    if "--ns" in sys.argv: NS = [int(x) for x in sys.argv[sys.argv.index("--ns") + 1].split(",")]
    seeds, nst = ([0], 4) if smoke else ([0, 1, 2], 20)
    SIG = [0.0] if smoke else [0.0, 0.02]
    out = os.path.join(HERE, "bench_entrainement_smoke.json" if smoke else "bench_entrainement_results.json")
    res, defauts = [], {}
    if os.path.exists(out) and "--reprise" in sys.argv:          # reprise : les cases déjà faites sont gardées
        J = json.load(open(out)); res, defauts = J["results"], J["defauts"]
    have = {(r["kind"], r["n"], r["sigma"]) for r in res}
    t0 = time.time()
    with Pool(2) as pool:
        for kind in KINDS:
            prob = "neurone-" + kind
            for n in NS:
                defauts[f"{kind}/{n}"] = {s: [defaut_cote(Neurone(n, kind, k).J(Neurone(n, kind, k).qs), s) for k in seeds] for s in "LR"}
                budget = 2000 if n > 8 else 600
                for sg in SIG:
                    if (kind, n, sg) in have:
                        continue
                    tasks = [(prob, n, s, i, m, mo, sg, budget, RAD) for m, mo in METH.items() for s in seeds for i in range(nst)
                             if not (m == "sonde-full" and n > 8)]
                    rr = pool.map(B.run, tasks, chunksize=2)
                    for r in rr: r["kind"] = kind
                    res += rr
                    d = defauts[f"{kind}/{n}"]
                    line = "  ".join(f"{m}:{B.cellule(B.resume(rr, m))}" for m in METH if B.resume(rr, m))
                    print(f"{kind:6s} n={n:2d} σ={sg} déf-droite={np.mean(d['R']):.2f} déf-gauche={np.mean(d['L']):.2f} | {line} [{time.time()-t0:.0f}s]", flush=True)
                json.dump(dict(results=res, defauts=defauts), open(out, "w"))
