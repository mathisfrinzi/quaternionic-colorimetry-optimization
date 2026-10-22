"""Partie 3.2 : la même famille que bench_hn.py, en grande dimension (n = 8, 16, 32, soit 32 à 128 réels).
Version vectorisée : A est stockée comme la matrice réelle 4n×4n de h ↦ A·h (blocs L(aᵢⱼ)). Mêmes sondes, budget 2000.
On garde les sondes qui passent à l'échelle (sec-full, sec-H, sec-H2full, fd-H) ; fd-full (4n + 1 lectures par pas) et BFGS n'ont plus
de sens à n = 32. Broyden exact s'initialise par différences finies (4n lectures) : c'est la boîte noire honnête.
Usage : python bench_hn_grand.py [--smoke]"""
import os, sys, json
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "..", "src")); sys.path.insert(0, HERE); sys.path.insert(0, os.path.join(HERE, "..", "3.1_jacobienne_H_lineaire"))
from multiprocessing import Pool
from quat import rand_unit
from sensors import Sensor
import classic
from probes import Hit, OutOfBudget
from bench_hn import MeterN, probe_hn, Lblock


def problem(n, eps, seed):
    """Même définition que bench_hn.problem (mêmes graines), calculée par matrices."""
    rng = np.random.default_rng(5000 + 97 * seed + 13 * n)
    A = np.zeros((n, n, 4))
    for i in range(n):
        for j in range(n):
            A[i, j] = 0.35 / np.sqrt(n) * rng.standard_normal(4)
        A[i, i] += rand_unit(rng) * 1.0
    G = rng.standard_normal((4 * n, 4 * n)) / np.sqrt(4 * n)
    qs = 0.5 * rng.standard_normal(4 * n)
    LA = Lblock(A)
    def f(q):
        d = (np.asarray(q, float) - qs).reshape(n, 4)
        s = 1.0 + 0.5 * np.sum(d * d, axis=1)
        return LA @ (d * s[:, None]).reshape(-1) + eps * (G @ d.reshape(-1))
    return f, qs


def starts(n, qs, k, seed):
    rng = np.random.default_rng(777 + seed)
    out = []
    for _ in range(k):
        v = rng.standard_normal(4 * n); v /= np.linalg.norm(v)
        out.append(qs + 1.5 * np.sqrt(n) * v)
    return out


def run(args):
    name, params, n, eps, seed, i, q0, bits, sigma, budget = args
    f, qs = problem(n, eps, seed)
    mode = "color" if name.startswith(("fd-", "sec-")) else ("degraded" if name.endswith("-deg") else "exact")
    sens = [Sensor("s2", bits, "sph", sigma=(sigma,) * 3, seed=1000 * i + k) for k in range(n)]
    M = MeterN(f, n, mode, sens, budget)
    try:
        if name.startswith(("fd-", "sec-")):
            # le plafond de pas suit la distance typique au zéro, qui croît en √n (inchangé pour n ≤ 4, comme en 8a)
            probe_hn(M, q0, name, n, **dict(params, cap=params["cap"] * max(1.0, np.sqrt(n / 4))))
        else:
            classic.broyden(M, q0, **params)
    except Hit:
        return True, M.cnt
    except (OutOfBudget, np.linalg.LinAlgError, FloatingPointError, ValueError, ZeroDivisionError, OverflowError):
        return False, budget
    return False, budget


PC = [dict(c=0.3, step=0.65, cap=1.5), dict(c=0.3, step=0.8, cap=3.0), dict(c=0.3, step=0.9, cap=3.0)]
METHODS = {"fd-H": PC, "sec-full": PC, "sec-H": PC, "sec-H2full": PC,
           "broyden": [dict(hfd=1e-6), dict(hfd=1e-4)], "broyden-deg": [dict(hfd=1e-2), dict(hfd=1e-1)]}


def evaluate(name, p, n, eps, sg, seeds, k, budget, pool):
    tasks = []
    for seed in seeds:
        f, qs = problem(n, eps, seed)
        for i, q0 in enumerate(starts(n, qs, k, seed)):
            tasks.append((name, p, n, eps, seed, i, q0, 24, sg, budget))
    res = pool.map(run, tasks, chunksize=2)
    return np.array([r[0] for r in res]), np.array([r[1] for r in res], float)


if __name__ == "__main__":
    smoke = "--smoke" in sys.argv
    NS = [8] if smoke else [8, 16, 32]
    EPSS = [0.0] if smoke else [0.0, 0.3]
    SIGS = [0.0] if smoke else [0.0, 0.02]
    seeds, k, budget = ([0], 3, 2000) if smoke else ([0, 1, 2], 10, 2000)
    out = []
    with Pool(2) as pool:
        for n in NS:
            for eps in EPSS:
                for sg in SIGS:
                    for name, cfgs in METHODS.items():
                        best = None
                        for p in cfgs:
                            ok, c = evaluate(name, p, n, eps, sg, [100], 3 if smoke else 5, budget, pool)
                            key = (ok.mean(), -np.median(np.where(ok, c, budget)))
                            if best is None or key > best[0]: best = (key, p)
                        ok, c = evaluate(name, best[1], n, eps, sg, seeds, k, budget, pool)
                        out.append(dict(n=n, eps=eps, bits=24, sigma=sg, method=name, params=best[1], rate=float(ok.mean()),
                                        cost=float(np.median(c[ok])) if ok.any() else None, ok=ok.tolist(), costs=c.tolist()))
                        print(n, eps, sg, f"{name:12s} {100*ok.mean():4.0f}%  {out[-1]['cost']}", flush=True)
                        json.dump(out, open(os.path.join(HERE, "bench_hn_grand_smoke.json" if smoke else "bench_hn_grand_results.json"), "w"))
