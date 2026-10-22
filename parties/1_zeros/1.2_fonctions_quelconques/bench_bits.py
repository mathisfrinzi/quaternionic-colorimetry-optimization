"""Partie 1.2 (suite) : « le canal axe paie-t-il ses bits ? » Succès et coût en fonction du budget TOTAL de bits par lecture B.

Paramètres de sonde fixés (c = 0.3, pas = 0.65, cap = 1) pour ne pas mélanger l'effet du capteur et celui du réglage.
Familles à zéro isolé (lin, quad, cub, ill), 3 instances × 40 départs. Capteurs : S0 (4 composantes log-flottant), S1 (ℓ + θ, 2
équations par pas), S2-rgb (axe sur 3 canaux), S2-sph (axe en 2 angles), à B bits par lecture.
Usage : python bench_bits.py
"""
import os, sys, json
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "..", "src"))
import bench
from probes import probe_fd, probe_s1

BS = [16, 24, 32, 36, 48, 64]
KINDS = ["lin", "quad", "cub", "ill"]
FIX = dict(c=0.3, step=0.65, cap=1.0)


def register_all():
    R = bench.register
    for B in BS:
        R(f"S0@{B}", probe_fd, lambda rng: {}, ncfg=1, sensor=dict(kind="s0", bits=B), fixed={**FIX, "model": "full"})
        R(f"S1@{B}", probe_s1, lambda rng: {}, ncfg=1, sensor=dict(kind="s1", bits=B), fixed=dict(c=0.3, step=0.65, cap=1.0))
        R(f"S2rgb@{B}", probe_fd, lambda rng: {}, ncfg=1, sensor=dict(kind="s2", bits=B, enc="rgb"), fixed={**FIX, "model": "full"})
        R(f"S2sph@{B}", probe_fd, lambda rng: {}, ncfg=1, sensor=dict(kind="s2", bits=B, enc="sph"), fixed={**FIX, "model": "full"})
    R("exacte", probe_fd, lambda rng: {}, ncfg=1, sensor=dict(kind="s2", bits=24, exact=True), fixed={**FIX, "model": "full"})


if __name__ == "__main__":
    register_all()
    methods = list(bench.REGISTRY)
    test = [("gen", k, s) for k in KINDS for s in (0, 1, 2)]
    res = bench.run_bench(methods, test[:1], test, ntune=1, ntest=40, budget=600, chosen={m: {} for m in methods})
    out = dict(bits=BS, kinds=KINDS, budget=600)
    for fam in ("S0", "S1", "S2rgb", "S2sph"):
        out[fam] = [bench.summarize(res, f"{fam}@{B}") for B in BS]
    out["exacte"] = bench.summarize(res, "exacte")
    json.dump(out, open(os.path.join(HERE, "bench_bits_results.json"), "w"))
    for fam in ("S0", "S1", "S2rgb", "S2sph"):
        print(f"{fam:6s}", "  ".join(f"B={B}:{100*s['rate']:.0f}%/{s['cost']:.0f}" for B, s in zip(BS, out[fam])))
    print("exacte", out["exacte"])
