"""Partie 6.4 : entraîner une COUCHE entière (m sorties), suite de la partie 6.3 (un seul neurone).

Couche y = σ(W·x + c), W ∈ ℍ^{m×p} (le poids est à GAUCHE de l'entrée), p exemples (problème carré par sortie, professeur W*).
Résidus r_{i,t}(W) = σ(Σⱼ W_ij·x_tj + c_i) − y_{t,i} : n = m·p quaternions d'inconnues et de résidus, dimension 4mp.
Deux structures s'empilent :
  1. à droite (partie 6.3) : ∂r_{i,t}/∂W_ij · h = Dσ·(h·x_tj), des blocs R(x_tj) ;
  2. diagonale par sortie : r_{i,t} ne dépend que de la ligne i de W, donc la jacobienne a m blocs diagonaux p × p.
Avec la couleur, la structure 2 se lit gratuitement : une sonde sur le poids j de TOUTES les sorties à la fois donne la colonne j
de tous les blocs (les sorties ne se mélangent pas). Lecture : p sondes au lieu de m·p (structure à droite seule) ou 4mp (rien).

Méthodes (couleur-complète-angles 24 bits par résidu) :
  sonde-Dbloc   sonde-couleur-mémoire-quaternionique-à-droite-par-sortie (les deux structures ; src/sonde_hn.py, modèle "Dbloc")
  sonde-D       la sonde de la partie 6.3 (à droite seulement, m·p sondes)
  sonde-full    sans structure (4mp sondes ; seulement si 4mp ≤ 64)
  broyden       Broyden-sans-couleur (valeur exacte, différences finies 4mp lectures)
  broyden-blocs Broyden-sans-couleur-par-sortie : valeur exacte ET structure diagonale connue (4p lectures de départ, mises à jour
                bloc par bloc) : la référence classique qui sait la même chose que la sonde
  broyden-deg   Broyden-mêmes-mesures
Usage : python bench_couche.py [--smoke] [--reprise] [--kinds lin,radial] [--cases 1x4,4x4] [--sigmas 0,0.02] [--procs N]   (cas = m x p)"""
import os, sys, json, time
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "src")); sys.path.insert(0, HERE)
from multiprocessing import Pool
import banc_hn as B
import sonde_hn as S
from quat import qmul, Rmat, rand_unit
from rnn import act, dact_block
from bench_entrainement import defaut_cote


class Couche:
    def __init__(self, m, p, kind="radial-faible", seed=0):
        rng = np.random.default_rng(64000 + 101 * seed + 7 * p + 1009 * m)
        sc = 0.3 if kind.endswith("-faible") else 1.0
        self.kind = kind.replace("-faible", ""); self.m, self.p = m, p
        self.Ws = np.array([[rand_unit(rng) for _ in range(p)] for _ in range(m)]) * rng.uniform(0.6, 1.2, (m, p, 1))
        self.X = sc * rng.standard_normal((p, p, 4)) / np.sqrt(2 * p)          # p exemples × p entrées (les mêmes pour toutes les sorties)
        self.c = sc * 0.3 * rng.standard_normal((m, 4))
        self.Y = self._out(self.Ws)                                             # (m, p, 4) : sortie i, exemple t
        self.qs = self.Ws.reshape(-1)

    def _z(self, W):
        return np.array([[sum(qmul(W[i, j], self.X[t, j]) for j in range(self.p)) + self.c[i] for t in range(self.p)] for i in range(self.m)])

    def _sig(self, z):
        if self.kind == "lin": return z
        return np.array([[act(z[i, t], self.kind) for t in range(self.p)] for i in range(self.m)])

    def _out(self, W):
        return self._sig(self._z(W))

    def f(self, q):
        W = np.asarray(q, float).reshape(self.m, self.p, 4)
        return (self._out(W) - self.Y).reshape(-1)

    def J(self, q):
        m, p = self.m, self.p; W = np.asarray(q, float).reshape(m, p, 4); Z = self._z(W)
        Jm = np.zeros((4 * m * p, 4 * m * p))
        for i in range(m):
            for t in range(p):
                Dm = np.eye(4) if self.kind == "lin" else dact_block(Z[i, t], self.kind)
                for j in range(p):
                    r, c = i * p + t, i * p + j
                    Jm[4 * r:4 * r + 4, 4 * c:4 * c + 4] = Dm @ Rmat(self.X[t, j])
        return Jm


def broyden_blocs(M, q0, m=1, bt_max=6, hfd=1e-6):
    """Broyden-sans-couleur qui connaît la structure diagonale par sortie : départ par différences finies en 4p lectures (la
    composante k du poids j perturbée pour toutes les sorties à la fois), puis Broyden « bon » sur chaque bloc (Schubert)."""
    q = np.array(q0, float); F = M.read(q); n4 = len(q); p4 = n4 // m
    def lire():
        Jb = np.zeros((m, p4, p4))
        for c in range(p4):
            e = np.zeros(n4); e[c::p4] = hfd
            dF = (M.read(q + e) - F) / hfd
            for i in range(m):
                Jb[i, :, c] = dF[i * p4:(i + 1) * p4]
        return Jb
    Jb = lire(); nrm = np.linalg.norm
    while True:
        d = np.concatenate([-np.linalg.lstsq(Jb[i], F[i * p4:(i + 1) * p4], rcond=1e-10)[0] for i in range(m)])
        t, Fn = 1.0, None
        for _ in range(bt_max + 1):
            Fn = M.read(q + t * d)
            if nrm(Fn) < nrm(F): break
            t /= 2; Fn = None
        if Fn is None:
            Jb = lire(); continue
        s, y = t * d, Fn - F
        for i in range(m):
            si, yi = s[i * p4:(i + 1) * p4], y[i * p4:(i + 1) * p4]
            if si @ si > 0:
                Jb[i] += np.outer(yi - Jb[i] @ si, si) / (si @ si)
        q, F = q + s, Fn


def p_couche(kind):
    def build(n, seed, m):
        L = Couche(m, n // m, kind, seed); return L.f, L.qs, {}
    return build


RAD = 0.5


def run(a):
    kind, m, p, seed, i, meth, sigma, budget = a
    L = Couche(m, p, kind, seed); n = m * p
    q0 = B.depart(L.qs, n, i, seed, RAD)
    if meth == "broyden-blocs":
        M = S.MeterLog(L.f, n, "exact", None, budget)
        try:
            broyden_blocs(M, q0, m); ok = False
        except S.Hit:
            ok = True
        except S.ERRS:
            ok = False
        c = M.cnt if ok else budget
    else:
        genre, opts = {"sonde-Dbloc": ("sonde", dict(modele="Dbloc", m=m)), "sonde-D": ("sonde", dict(modele="D")),
                       "sonde-full": ("sonde", dict(modele="full")), "broyden": ("broyden", dict(hfd=1e-6)),
                       "broyden-deg": ("broyden-deg", dict(hfd=1e-2))}[meth]
        if genre == "sonde":
            opts = {**B.P0, **opts}; opts["cap"] = opts["cap"] * max(1.0, np.sqrt(n / 4))
        ok, c, *_ = S.lancer(L.f, q0, n, genre, budget=budget, sigma=sigma, i=i, **opts)
    return dict(kind=kind, m=m, p=p, n=n, seed=seed, i=i, meth=meth, sigma=sigma, ok=bool(ok), cost=int(c))


METH = ["sonde-Dbloc", "sonde-D", "sonde-full", "broyden", "broyden-blocs", "broyden-deg"]

if __name__ == "__main__":
    smoke = "--smoke" in sys.argv
    CASES = [(1, 4), (2, 4), (4, 4)] if smoke else [(1, 4), (2, 4), (4, 4), (8, 4), (2, 8), (4, 8), (8, 8)]   # (m, p)
    KINDS = ("radial-faible",) if smoke else ("lin", "radial-faible", "radial")
    seeds, nst = ([0], 3) if smoke else ([0, 1, 2], 10)
    SIG = [0.0] if smoke else [0.0, 0.02]
    arg = lambda k: sys.argv[sys.argv.index(k) + 1] if k in sys.argv else None
    if arg("--kinds"): KINDS = tuple(arg("--kinds").split(","))
    if arg("--cases"): CASES = [tuple(int(v) for v in c.split("x")) for c in arg("--cases").split(",")]
    if arg("--sigmas"): SIG = [float(v) for v in arg("--sigmas").split(",")]
    out = os.path.join(HERE, "bench_couche_smoke.json" if smoke else "bench_couche_results.json")
    res, defauts = [], {}
    if os.path.exists(out) and "--reprise" in sys.argv:
        J = json.load(open(out)); res, defauts = J["results"], J["defauts"]
    have = {(r["kind"], r["m"], r["p"], r["sigma"]) for r in res}
    t0 = time.time()
    with Pool(int(arg("--procs") or 2)) as pool:
        for kind in KINDS:
            for m, p in CASES:
                n = m * p
                L0 = [Couche(m, p, kind, s) for s in seeds]
                defauts[f"{kind}/{m}/{p}"] = dict(R=[defaut_cote(L.J(L.qs), "R") for L in L0])
                budget = 600 if n <= 8 else 2000
                for sg in SIG:
                    if (kind, m, p, sg) in have: continue
                    meths = [x for x in METH if not (x == "sonde-full" and n > 16) and not (x == "sonde-D" and n > 32)]
                    tasks = [(kind, m, p, s, i, x, sg, budget) for x in meths for s in seeds for i in range(nst)]
                    rr = pool.map(run, tasks, chunksize=1)
                    res += rr
                    line = "  ".join(f"{x}:{B.cellule(B.resume(rr, x))}" for x in meths)
                    print(f"{kind:13s} m={m} p={p} (dim {4*n:3d}) σ={sg} déf-droite={np.mean(defauts[f'{kind}/{m}/{p}']['R']):.2f} | {line} [{time.time()-t0:.0f}s]", flush=True)
                    json.dump(dict(results=res, defauts=defauts), open(out, "w"))
