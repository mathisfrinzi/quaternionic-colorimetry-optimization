"""Partie 1.2 (M2) : la sonde 4 directions sur fonctions génériques de R^4 dans R^4.

Familles : zéro isolé (lin, quad, cub, ill), courbe de zéros, surface de zéros, aucun zéro. 3 instances aléatoires par famille
au test (graines 0-2), 2 autres au réglage (graines 10-11). Budget commun : 600 lectures.
Les sondes lisent un capteur quantifié (24 ou 48 bits par lecture) ; les méthodes classiques LISENT LA VALEUR EXACTE
(comparaison favorable à la méthode classique ; la partie 6 égalise les mesures).
Usage : python bench_generique.py [--smoke]
"""
import os, sys, json
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", "..", ".."))
sys.path.insert(0, os.path.join(ROOT, "src"))
import bench, classic
from probes import probe_fd, probe_secant, probe_s1
from funcs import GENERIC_KINDS

pick = lambda rng, xs: xs[int(rng.integers(len(xs)))]


def fd_sample(rng):
    return dict(c=pick(rng, [0.15, 0.3, 0.45]), step=pick(rng, [0.5, 0.65, 0.8]), cap=pick(rng, [0.5, 1.0, 2.0]))


def register_all():
    R = bench.register
    # sondes à lectures de couleur
    R("sonde-S0-24b", probe_fd, fd_sample, ncfg=4, sensor=dict(kind="s0", bits=24), fixed=dict(model="full"))
    R("sonde-S0-48b", probe_fd, fd_sample, ncfg=6, sensor=dict(kind="s0", bits=48), fixed=dict(model="full"))
    R("sonde-S1-24b", probe_s1, fd_sample, ncfg=6, sensor=dict(kind="s1", bits=24))
    R("sonde-S1-48b", probe_s1, fd_sample, ncfg=6, sensor=dict(kind="s1", bits=48))
    R("sonde-S2rgb-24b", probe_fd, fd_sample, ncfg=6, sensor=dict(kind="s2", bits=24, enc="rgb"), fixed=dict(model="full"))
    R("sonde-S2rgb-48b", probe_fd, fd_sample, ncfg=6, sensor=dict(kind="s2", bits=48, enc="rgb"), fixed=dict(model="full"))
    R("sonde-S2sph-24b", probe_fd, fd_sample, ncfg=6, sensor=dict(kind="s2", bits=24, enc="sph"), fixed=dict(model="full"))
    R("sonde-exacte", probe_fd, fd_sample, ncfg=6, sensor=dict(kind="s2", bits=24, exact=True), fixed=dict(model="full"))
    R("sonde-S2rgb-24b-Broyden", probe_secant, fd_sample, ncfg=6, sensor=dict(kind="s2", bits=24, enc="rgb"), fixed=dict(model="rank1"))
    R("sonde-S2rgb-48b-Broyden", probe_secant, fd_sample, ncfg=6, sensor=dict(kind="s2", bits=48, enc="rgb"), fixed=dict(model="rank1"))
    # classiques (lecture exacte)
    R("newton-FD", classic.newton_fd, lambda rng: dict(hfd=pick(rng, [1e-6, 1e-4])), ncfg=2)
    R("newton-Jexacte", classic.newton_exactJ, lambda rng: {}, ncfg=1)
    R("broyden", classic.broyden, lambda rng: dict(hfd=pick(rng, [1e-6, 1e-4])), ncfg=2)
    R("bfgs-FD", classic.bfgs_fd, lambda rng: {}, ncfg=1)
    R("lm-FD", classic.lm_fd, lambda rng: {}, ncfg=1)
    R("nelder-mead", classic.nelder_mead, lambda rng: dict(simplex=pick(rng, [0.3, 1.0]), loss=pick(rng, ["log", "sq"])), ncfg=3)
    R("powell", classic.powell, lambda rng: dict(loss=pick(rng, ["log", "sq"])), ncfg=2)
    R("cma-es", classic.cma_es, lambda rng: dict(sigma=pick(rng, [0.3, 1.0]), pop=pick(rng, [6, 10]), loss=pick(rng, ["log", "sq"])), ncfg=3)


if __name__ == "__main__":
    smoke = "--smoke" in sys.argv
    register_all()
    methods = list(bench.REGISTRY)
    tune = [("gen", k, s) for k in GENERIC_KINDS for s in (10, 11)]
    test = [("gen", k, s) for k in GENERIC_KINDS for s in (0, 1, 2)]
    res = bench.run_bench(methods, tune, test, ntune=2 if smoke else 5, ntest=3 if smoke else 40, budget=600, ncfg=1 if smoke else None)
    json.dump(res, open(os.path.join(HERE, "bench_generique_smoke.json" if smoke else "bench_generique_results.json"), "w"))
    for m in methods:
        s = bench.summarize(res, m)
        print(f"{m:26s} {100*s['rate']:5.1f}% ({100*s['lo']:.0f}-{100*s['hi']:.0f})  coût médian {s['cost']:.0f}")
