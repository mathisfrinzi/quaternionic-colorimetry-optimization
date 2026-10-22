"""Partie 2 (version réduite, voir README) : minimum NON NUL. f : H -> H sans zéro ; on minimise l = |f|^2, minimum |f| = 1.
Succès : une lecture avec |f| < 1.05 (perte à 10 % du minimum). La sonde à cible (repo complexe, vendored) ne lit que la
COULEUR (luminosité) du scalaire l - m0 ; le gradient vient de différences de quotients, la hessienne de BFGS.
Comparaison à lectures égales : bfgs-FD, lm-FD, nelder-mead, cma-es (lectures exactes de f).
DÉVIATION du plan : on n'exploite pas le canal d'axe de S2 pour la partie 4 (l est un scalaire) ; l'ablation axe/teinte face au
bruit n'est pas faite ici. Familles : nozero (conditionnement 3), nozero_ill (100). 3 instances × 30 départs, budget 600."""
import os, sys, json
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "src")); sys.path.insert(0, os.path.join(HERE, "..", "..", "vendor", "complex_colorimetry"))
import bench, classic
from sonde_cible import sonde_cible, Arret
from probes import OutOfBudget
pick = lambda rng, xs: xs[int(rng.integers(len(xs)))]


def cible(M, q0, bits=24, c=0.1, L=1.0):
    def m(x):
        F = M.read(x)
        return float(F @ F)
    try:
        sonde_cible(m, np.array(q0, float), c=c, bits=bits, L=L)
    except Arret:
        raise OutOfBudget


for B in (24, 36):
    bench.register(f"cible-{B}", cible, lambda rng: dict(c=pick(rng, [0.05, 0.1, 0.3]), L=pick(rng, [0.5, 1.0, 2.0])), ncfg=6, fixed=dict(bits=B))
bench.register("bfgs-FD", classic.bfgs_fd, lambda rng: dict(hfd=pick(rng, [1e-6, 1e-4])), ncfg=2)
bench.register("lm-FD", classic.lm_fd, lambda rng: dict(diff_step=pick(rng, [1e-6, 1e-4])), ncfg=2)
bench.register("nelder-mead", classic.nelder_mead, lambda rng: dict(simplex=pick(rng, [0.3, 1.0]), loss="sq"), ncfg=2)
bench.register("cma-es", classic.cma_es, lambda rng: dict(sigma=pick(rng, [0.3, 0.6, 1.0]), pop=8), ncfg=3)

if __name__ == "__main__":
    smoke = "--smoke" in sys.argv
    methods = list(bench.REGISTRY)
    fam = [("gen", k, s) for k in ("nozero", "nozero_ill") for s in (0, 1, 2)]
    res = bench.run_bench(methods, [("gen", "nozero", 10), ("gen", "nozero_ill", 10)], fam, ntune=2 if smoke else 6, ntest=3 if smoke else 30, budget=600, ncfg=1 if smoke else None)
    json.dump(res, open(os.path.join(HERE, "bench_smoke.json" if smoke else "bench_minimum_results.json"), "w"))
    for m in methods:
        for k in ("nozero", "nozero_ill"):
            keys = [x for x in res["data"][m] if k + '"' in x]
            s = bench.summarize(res, m, keys); print(f"{m:14s} {k:11s} {100*s['rate']:5.1f}%  coût médian {s['cost']:.0f}")
