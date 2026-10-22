"""Partie 1.5 (M5) : lire le degré topologique dans les couleurs.

5a  récupération du degré (+1, -1, +2, -2, +3 et 0) sur la sphère de rayon 1, en fonction des bits, du bruit, de la profondeur
    de triangulation (32 ou 192 lectures) et de l'encodage de l'axe (RGB ou 2 angles). Fonctions perturbées au hasard
    (f = base + 0.15 A q, ‖A‖ = 1 : le degré est préservé).
5b  détection d'un zéro : f(q) = base(q - xs), sphère de rayon 1 centrée en 0 ; un degré non nul prouve que |xs| < 1.
    Taux de bonne classification (dedans / dehors) en fonction de |xs|.
5c  la descente de la sonde jusqu'à un zéro de degré -2 (zéro double : la jacobienne s'annule), puis certificat autour du point trouvé.
Le capteur S1 (sans axe) ne permet pas de lire la direction de f, donc pas de degré : il n'apparaît pas ici.
Usage : python bench_degre.py [--smoke]
"""
import os, sys, json
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "..", "src"))
from quat import qmul, qconj, rand_unit
from sensors import Sensor
from degree import certificate, n_reads
from probes import probe_fd, probe_secant, Meter, Hit, OutOfBudget

BASES = {
    "q": (lambda q: np.asarray(q, float), 1),
    "q̄": (lambda q: qconj(q), -1),
    "q²": (lambda q: qmul(q, q), 2),
    "q̄²": (lambda q: qmul(qconj(q), qconj(q)), -2),
    "q³": (lambda q: qmul(q, qmul(q, q)), 3),
    "q+3 (sans zéro)": (lambda q: np.asarray(q, float) + np.array([3.0, 0, 0, 0]), 0),
}


def perturbed(base, rng, amp=0.15):
    A = rng.standard_normal((4, 4)); A /= np.linalg.norm(A, 2)
    return lambda q: base(q) + amp * (A @ np.asarray(q, float))


def sensor(enc, B, sigma, seed):
    return Sensor("s2", B, enc, sigma=(sigma, sigma, sigma), seed=seed)


def run_5a(smoke):
    npert, nseed = (1, 1) if smoke else (4, 3)
    cfgs = []
    for enc in ("rgb", "sph"):
        for B in (12, 16, 24, 36, 48):
            for depth in (1, 2):
                cfgs.append(dict(enc=enc, B=B, sigma=0.0, depth=depth))
    for enc in ("rgb", "sph"):
        for sg in (0.01, 0.03, 0.1, 0.2):
            for depth in (1, 2):
                cfgs.append(dict(enc=enc, B=48, sigma=sg, depth=depth))
    out = []
    for cfg in cfgs:
        ok, tot, agree, per_base = 0, 0, [], {k: [0, 0] for k in BASES}
        for bi, (bn, (base, d)) in enumerate(BASES.items()):
            for p in range(npert):
                f = perturbed(base, np.random.default_rng(100 * bi + p))
                for s in range(nseed):
                    sen = sensor(cfg["enc"], cfg["B"], cfg["sigma"], seed=s + 7 * p)
                    got, ag, nr = certificate(f, sen, np.zeros(4), 1.0, depth=cfg["depth"], seed=s)
                    good = int(got == d)
                    ok += good; tot += 1; agree.append(ag); per_base[bn][0] += good; per_base[bn][1] += 1
        out.append(dict(cfg, reads=n_reads(cfg["depth"]), rate=ok / tot, n=tot, agree=float(np.mean(agree)),
                        per_base={k: v[0] / v[1] for k, v in per_base.items()}))
        print(f"5a {cfg['enc']:3s} B={cfg['B']:2d} σ={cfg['sigma']:.2f} depth={cfg['depth']} ({out[-1]['reads']:3d} lectures) : "
              f"degré exact {100*out[-1]['rate']:.0f} %", flush=True)
    return out


def run_5b(smoke):
    ds = [0.2, 0.5, 0.8, 0.9, 1.1, 1.2, 1.5, 1.8]
    ndir, nseed = (2, 1) if smoke else (10, 3)
    cfgs = [dict(enc="sph", B=24, sigma=0.0, depth=2), dict(enc="sph", B=24, sigma=0.03, depth=2),
            dict(enc="sph", B=36, sigma=0.1, depth=2), dict(enc="sph", B=24, sigma=0.0, depth=1),
            dict(enc="rgb", B=24, sigma=0.0, depth=2)]
    out = []
    for cfg in cfgs:
        res = {}
        for bn in ("q", "q̄²"):
            base, deg = BASES[bn]
            acc, degok = [], []
            for d in ds:
                a, g = [], []
                for k in range(ndir):
                    xs = d * rand_unit(np.random.default_rng(1000 * k + int(100 * d)))
                    f = lambda q, xs=xs, base=base: base(np.asarray(q, float) - xs)
                    for s in range(nseed):
                        sen = sensor(cfg["enc"], cfg["B"], cfg["sigma"], seed=s + 13 * k)
                        got, ag, nr = certificate(f, sen, np.zeros(4), 1.0, depth=cfg["depth"], seed=s)
                        inside = d < 1.0
                        a.append((got != 0) == inside)
                        if inside:
                            g.append(got == deg)
                acc.append(float(np.mean(a))); degok.append(float(np.mean(g)) if g else None)
            res[bn] = dict(acc=acc, degree_ok=degok)
        out.append(dict(cfg, reads=n_reads(cfg["depth"]), d=ds, res=res))
        print(f"5b {cfg} classification q: {[round(x, 2) for x in res['q']['acc']]}", flush=True)
    return out


def run_5c(smoke):
    """La sonde converge vers le zéro DOUBLE de q̄² (degré -2) ; certificat à 36 bits autour du point trouvé."""
    f0 = BASES["q̄²"][0]
    rng = np.random.default_rng(5)
    rows = []
    for i in range(3 if smoke else 20):
        xs = rng.standard_normal(4) * 0.5
        f = lambda q, xs=xs: f0(np.asarray(q, float) - xs)
        q0 = xs + 1.2 * rand_unit(rng)
        sen = Sensor("s2", 36, "sph", seed=i)
        M = Meter(f, sen, budget=800, thr=1e-3)
        try:
            probe_secant(M, q0, model="rank1", c=0.3, step=0.5, cap=1.0)
            continue
        except Hit:
            q_end, cost = M.last, M.n
        except OutOfBudget:
            continue
        sen2 = Sensor("s2", 36, "sph", seed=100 + i)
        got, ag, nr = certificate(f, sen2, q_end, 0.5, depth=2)
        rows.append(dict(cost=cost, dist=float(np.linalg.norm(q_end - xs)), degree=got, agree=ag, reads=nr))
    print(f"5c {len(rows)} descentes réussies ; degrés lus : {sorted(set(r['degree'] for r in rows))}", flush=True)
    return rows


if __name__ == "__main__":
    smoke = "--smoke" in sys.argv
    out = dict(a=run_5a(smoke), b=run_5b(smoke), c=run_5c(smoke))
    json.dump(out, open(os.path.join(HERE, "bench_degre_smoke.json" if smoke else "bench_degre_results.json"), "w"), default=lambda o: o.item() if hasattr(o, "item") else str(o))
