"""Partie 4.2 : points fixes de réservoirs quaternioniques, n = 2 à 32 (dimension 8 à 128), et défaut ℍ-linéaire contrôlé.

A. banc principal : n ∈ {2, 8, 16, 32}, activation ∈ {split, radial, global}, rayon spectral ∈ {0,8 ; 0,95}, entrée β = 1.
B. défaut contrôlé (n = 8, rayon 0,95) : force de l'entrée β ∈ {0,25 … 4} (le point de fonctionnement de l'activation), et une pièce
   réelle non quaternionique ajoutée à W, ‖P‖ = t·‖W‖, t ∈ {0,1 ; 0,3 ; 0,6} (ce que deviendrait une pièce apprise sans contrainte).

Sondes de couleur (S2 sphérique 24 bits par composante) : sec-H, sec-H2full, sec-full (partie 8), et sec-Hid : même a priori que
broyden-id (jacobienne ≈ diagonale par blocs), lu avec une seule sonde simultanée (2 lectures), puis Broyden quaternionique.
Références, lectures EXACTES de F : broyden (jacobienne initiale par différences finies, 4n lectures), broyden-id (jacobienne
initiale = identité, 0 lecture : le bon choix pour un point fixe, car J = I − Dσ·W), picard (h ← σ(Wh + b), une lecture par pas).
Référence dégradée : broyden-deg. 3 réservoirs × 10 départs (h* + perturbation de norme 0,5·√n), budget 2000.
Usage : python bench_reservoirs.py [--smoke]"""
import os, sys, json
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "..", "src")); sys.path.insert(0, os.path.join(HERE, "..", "..", "3_dimension_4n", "3.1_jacobienne_H_lineaire"))
from multiprocessing import Pool
from sensors import Sensor
import classic
from probes import Hit, OutOfBudget
from bench_hn import MeterN, probe_hn
from reservoir import Reservoir

_CACHE = {}


def res(key):
    if key not in _CACHE:
        n, rho, beta, kind, t, seed = key
        _CACHE[key] = Reservoir(n, rho, beta, kind, seed, t)
    return _CACHE[key]


def broyden_id(M, q0, bt_max=6):
    """Broyden « bon » avec jacobienne initiale = identité (aucune lecture d'initialisation)."""
    q = np.array(q0, float); F = M.read(q); H = np.eye(len(q))
    while True:
        d = -H @ F
        t, Fn = classic._line(M, q, F, d, bt_max)
        if t is None:
            H = np.eye(len(q)); continue
        s, y = t * d, Fn - F; Hy = H @ y; den = s @ Hy
        if abs(den) > 1e-14:
            H = H + np.outer(s - Hy, s @ H) / den
        q, F = q + s, Fn


def picard(M, q0):
    q = np.array(q0, float)
    while True:
        q = q - M.read(q)                 # h ← h − F(h) = σ(W h + b)


def starts(R, k, seed):
    rng = np.random.default_rng(4242 + seed)
    out = []
    for _ in range(k):
        v = rng.standard_normal(4 * R.n); v /= np.linalg.norm(v)
        out.append(R.hs + 0.5 * np.sqrt(R.n) * v)
    return out


def run(args):
    name, params, key, i, q0, sigma, budget = args
    R = res(key); n = R.n
    mode = "color" if name.startswith("sec-") else ("degraded" if name.endswith("-deg") else "exact")
    sens = [Sensor("s2", 24, "sph", sigma=(sigma,) * 3, seed=1000 * i + k) for k in range(n)]
    M = MeterN(R.F, n, mode, sens, budget)
    try:
        if name.startswith("sec-"): probe_hn(M, q0, name, n, **dict(params, cap=params["cap"] * max(1.0, np.sqrt(n / 4))))
        elif name == "broyden-id": broyden_id(M, q0)
        elif name == "picard": picard(M, q0)
        else: classic.broyden(M, q0, **params)
    except Hit:
        return True, M.cnt
    except (OutOfBudget, np.linalg.LinAlgError, FloatingPointError, ValueError, ZeroDivisionError, OverflowError):
        return False, budget
    return False, budget


PC = [dict(c=0.3, step=0.65, cap=1.5), dict(c=0.3, step=0.8, cap=3.0), dict(c=0.3, step=0.9, cap=3.0)]
METHODS = {"sec-full": PC, "sec-H": PC, "sec-H2full": PC, "sec-Hid": PC, "broyden": [dict(hfd=1e-6), dict(hfd=1e-4)],
           "broyden-id": [{}], "picard": [{}], "broyden-deg": [dict(hfd=1e-2), dict(hfd=1e-1)]}


def evaluate(name, p, cfg, seeds, k, sigma, budget, pool):
    tasks = []
    for seed in seeds:
        key = (*cfg, seed); R = res(key)
        for i, q0 in enumerate(starts(R, k, seed)):
            tasks.append((name, p, key, i, q0, sigma, budget))
    out = pool.map(run, tasks, chunksize=2)
    return np.array([o[0] for o in out]), np.array([o[1] for o in out], float)


def bench(cfgs, methods, seeds, k, budget, pool, tag, out, smoke):
    for cfg in cfgs:                                   # cfg = (n, rho, beta, kind, t_real)
        n, rho, beta, kind, t = cfg
        defects = [res((*cfg, s)).defect() for s in seeds]
        for name in methods:
            best = None
            for p in METHODS[name]:
                ok, c = evaluate(name, p, cfg, [100], 3 if smoke else 5, 0.0, budget, pool)
                key = (ok.mean(), -np.median(np.where(ok, c, budget)))
                if best is None or key > best[0]: best = (key, p)
            ok, c = evaluate(name, best[1], cfg, seeds, k, 0.0, budget, pool)
            out.append(dict(part=tag, n=n, rho=rho, beta=beta, kind=kind, t_real=t, defect=float(np.mean(defects)), method=name,
                            params=best[1], rate=float(ok.mean()), cost=float(np.median(c[ok])) if ok.any() else None,
                            ok=ok.tolist(), costs=c.tolist()))
            print(tag, cfg, f"défaut {np.mean(defects):.3f}", f"{name:12s} {100*ok.mean():4.0f}%  {out[-1]['cost']}", flush=True)
        json.dump(out, open(os.path.join(HERE, "bench_smoke.json" if smoke else "bench_reservoirs_results.json"), "w"))


if __name__ == "__main__":
    smoke = "--smoke" in sys.argv
    seeds, k, budget = ([0], 3, 2000) if smoke else ([0, 1, 2], 10, 2000)
    out = []
    with Pool(2) as pool:
        NS = [2, 8] if smoke else [2, 8, 16, 32]
        cfgA = [(n, rho, 1.0, kind, 0.0) for n in NS for rho in ([0.95] if smoke else [0.8, 0.95]) for kind in ("split", "radial", "global")]
        bench(cfgA, list(METHODS), seeds, k, budget, pool, "A", out, smoke)
        if not smoke:
            cfgB = [(8, 0.95, beta, kind, 0.0) for kind in ("split", "radial", "global") for beta in (0.25, 0.5, 2.0, 4.0)]
            cfgB += [(8, 0.95, 1.0, "split", t) for t in (0.1, 0.3, 0.6)]
            bench(cfgB, ["sec-H", "sec-Hid", "sec-full", "broyden", "broyden-id", "picard"], seeds, k, budget, pool, "B", out, smoke)
