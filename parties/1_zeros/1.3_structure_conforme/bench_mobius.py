"""Partie 1.3 (M3) : le pari structurel. f_ε = Möbius + ε · perturbation lisse aléatoire (même zéro q*).

ε = 0 : f exactement conforme (jacobienne = L(α)R(β), 7 paramètres). ε grand : générique (défaut de conformité ≈ 0.6).
Question : à BUDGET DE LECTURES ÉGAL, un modèle local conforme (3 sondes au lieu de 4 ; ou 4 sondes + projection ; ou ajustement
(α, β) sur une fenêtre de sécantes) bat-il le même algorithme sans structure, et jusqu'à quelle distance de la conformité ?
Capteur S2-rgb à B ∈ {16, 24, 36} bits. Les méthodes classiques lisent la valeur exacte (favorable aux classiques).
Appariement : toutes les méthodes reçoivent exactement les mêmes départs (comparaison par départ possible).
Usage : python bench_mobius.py [--smoke]
"""
import os, sys, json
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", "..", ".."))
sys.path.insert(0, os.path.join(ROOT, "src"))
import bench, classic
from probes import probe_fd, probe_secant

pick = lambda rng, xs: xs[int(rng.integers(len(xs)))]
EPS = [0.0, 0.03, 0.1, 0.3, 1.0, 3.0]
BITS = [16, 24, 36]
SEEDS_TEST = [0, 1, 2, 3, 4]


def fd_sample(rng):
    return dict(c=pick(rng, [0.15, 0.3, 0.45]), step=pick(rng, [0.5, 0.65, 0.8]), cap=pick(rng, [0.5, 1.0, 2.0]))


def sec_sample(rng):
    return dict(**fd_sample(rng), K=pick(rng, [2, 3]), ridge=pick(rng, [1e-3, 1e-2]))


def rank1_sample(rng):
    return fd_sample(rng)


def register_all():
    R = bench.register
    for B in BITS:
        sen = dict(kind="s2", bits=B, enc="rgb")
        R(f"fd-full@{B}", probe_fd, fd_sample, ncfg=6, sensor=sen, fixed=dict(model="full"))
        R(f"fd-conf3@{B}", probe_fd, fd_sample, ncfg=6, sensor=sen, fixed=dict(model="conf3"))
        R(f"fd-confproj@{B}", probe_fd, fd_sample, ncfg=6, sensor=sen, fixed=dict(model="confproj"))
        R(f"sec-rank1@{B}", probe_secant, rank1_sample, ncfg=6, sensor=sen, fixed=dict(model="rank1"))
        R(f"sec-multi@{B}", probe_secant, sec_sample, ncfg=8, sensor=sen, fixed=dict(model="multi"))
        R(f"sec-conf@{B}", probe_secant, sec_sample, ncfg=8, sensor=sen, fixed=dict(model="conf"))
    R("broyden", classic.broyden, lambda rng: dict(hfd=pick(rng, [1e-6, 1e-4])), ncfg=2)
    R("newton-FD", classic.newton_fd, lambda rng: dict(hfd=pick(rng, [1e-6, 1e-4])), ncfg=2)
    R("newton-Jexacte", classic.newton_exactJ, lambda rng: {}, ncfg=1)
    R("bfgs-FD", classic.bfgs_fd, lambda rng: {}, ncfg=1)
    R("nelder-mead", classic.nelder_mead, lambda rng: dict(simplex=pick(rng, [0.3, 1.0]), loss=pick(rng, ["log", "sq"])), ncfg=3)


if __name__ == "__main__":
    smoke = "--smoke" in sys.argv
    register_all()
    methods = list(bench.REGISTRY)
    tune = [("mobius", s, e) for e in (0.0, 1.0) for s in (10, 11)]
    test = [("mobius", s, e) for e in EPS for s in SEEDS_TEST]
    res = bench.run_bench(methods, tune, test, ntune=2 if smoke else 4, ntest=3 if smoke else 30, budget=600, ncfg=1 if smoke else None)
    res["eps"] = EPS
    json.dump(res, open(os.path.join(HERE, "bench_mobius_smoke.json" if smoke else "bench_mobius_results.json"), "w"))
    for m in methods:
        s = bench.summarize(res, m)
        print(f"{m:20s} {100*s['rate']:5.1f}%  coût médian {s['cost']:.0f}")
