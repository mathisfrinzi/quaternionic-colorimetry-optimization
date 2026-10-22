"""Banc commun : protocole à coût d'évaluations égal, repris du repo complexe (bench2d.py).

  * Unité de coût : une lecture (= une évaluation de f). Budget commun pour toutes les méthodes.
  * Succès : une lecture tombe sur un point où |f| < seuil (vraie valeur, que l'algorithme ne voit pas) avant épuisement du budget.
  * Réglage : pour chaque méthode, `ncfg` configurations tirées au hasard, essayées sur des départs de RÉGLAGE (graine différente
    de celle du test) ; on garde celle qui a le meilleur succès puis le plus faible coût médian. Test sur d'autres départs.
  * Statistiques : taux de succès avec intervalle de Wilson 95 %, coût médian parmi les succès.

Usage typique (voir parties/*/bench_*.py) :
    bench.register("nom", kind="color"|"classic", fn=..., sample=lambda rng: {...}, ncfg=6, sensor=dict(kind="s2", bits=24))
    res = bench.run_bench(methods, problems, ...)
"""
import os, sys, json, math, time, functools
import numpy as np
from multiprocessing import Pool

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from funcs import build
from probes import Meter, Hit, OutOfBudget
from sensors import Sensor

np.seterr(all="ignore")
REGISTRY = {}


def register(name, fn, sample, ncfg=6, sensor=None, absolute=False, starts="ball", fixed=None):
    """fn(M, q0, **params). sensor : dict d'arguments de Sensor (None = lecture exacte). absolute : vue à échelle absolue
    (méthodes classiques de la partie 6). starts : 'ball' (boule de R^4) ou 'slice' (tranche C_i).
    fixed : paramètres imposés (non réglés), fusionnés avec ceux tirés par `sample`."""
    REGISTRY[name] = dict(fn=fn, sample=sample, ncfg=ncfg, sensor=sensor, absolute=absolute, starts=starts, fixed=fixed or {})


@functools.lru_cache(maxsize=None)
def _problem(spec):
    return build(spec)


def run_one(name, params, pspec, q0, budget, seed=0):
    R = REGISTRY[name]
    P = _problem(pspec)
    sen = None if R["sensor"] is None else Sensor(seed=seed, **R["sensor"])
    M = Meter(P.f, sen, budget=budget, thr=P.thr, absolute=R["absolute"])
    try:
        R["fn"](M, q0, **{**R["fixed"], **params})
    except Hit:
        return True, M.n
    except (OutOfBudget, np.linalg.LinAlgError, FloatingPointError, ValueError, ZeroDivisionError, OverflowError):
        return False, budget
    return False, budget


def _starts(P, mode, n, seed):
    return P.slice_starts(n, seed) if mode == "slice" else P.starts(n, seed)


def configs(name):
    R = REGISTRY[name]
    rng = np.random.default_rng(sum(map(ord, name)))
    return [R["sample"](rng) for _ in range(R["ncfg"])]


def tune_task(args):
    name, ci, params, pspecs, ntune, budget, seed = args
    R = REGISTRY[name]
    ok, cost = [], []
    for k, ps in enumerate(pspecs):
        P = _problem(ps)
        for i, q0 in enumerate(_starts(P, R["starts"], ntune, seed + 31 * k)):
            o, c = run_one(name, params, ps, q0, budget, seed=i)
            ok.append(o); cost.append(c if o else budget)
    return name, ci, float(np.mean(ok)), float(np.median(cost))


def test_task(args):
    name, params, ps, ntest, budget, seed = args
    R = REGISTRY[name]
    P = _problem(ps)
    res = [run_one(name, params, ps, q0, budget, seed=1000 + i) for i, q0 in enumerate(_starts(P, R["starts"], ntest, seed))]
    return name, ps, [bool(o) for o, _ in res], [int(c) for _, c in res]


def wilson(k, n, z=1.96):
    if n == 0:
        return float("nan"), float("nan")
    ph = k / n; d = 1 + z * z / n
    c = (ph + z * z / (2 * n)) / d
    h = z * math.sqrt(ph * (1 - ph) / n + z * z / (4 * n * n)) / d
    return c - h, c + h


def run_bench(methods, tune_specs, test_specs, ntune=6, ntest=60, budget=600, seed_tune=2, seed_test=1, ncfg=None,
              procs=2, out=None, verbose=True, chosen=None):
    """methods : liste de noms enregistrés. tune_specs / test_specs : listes de specs de problèmes (cf. funcs.build).
    chosen : dict {méthode: params} pour sauter le réglage. Retourne un dict sérialisable."""
    t0 = time.time()
    if ncfg is not None:
        for m in methods:
            REGISTRY[m]["ncfg"] = min(REGISTRY[m]["ncfg"], ncfg)
    with Pool(procs) as pool:
        if chosen is None:
            tasks = [(m, i, p, tune_specs, ntune, budget, seed_tune) for m in methods for i, p in enumerate(configs(m))]
            res = pool.map(tune_task, tasks, chunksize=1)
            best = {}
            for m, i, sr, md in res:
                key = (sr, -md)
                if m not in best or key > best[m][0]:
                    best[m] = (key, i)
            chosen = {m: configs(m)[best[m][1]] for m in methods}
            if verbose:
                print(f"  réglage fini en {time.time() - t0:.0f}s", flush=True)
        tasks = [(m, chosen[m], ps, ntest, budget, seed_test) for m in methods for ps in test_specs]
        res = pool.map(test_task, tasks, chunksize=1)
    if verbose:
        print(f"  test fini en {time.time() - t0:.0f}s", flush=True)
    data = {m: {} for m in methods}
    for m, ps, ok, cost in res:
        data[m][json.dumps(ps)] = dict(ok=ok, cost=cost)
    return dict(budget=budget, ntune=ntune, ntest=ntest, chosen=chosen, data=data, seconds=time.time() - t0)


def summarize(res, method, keys=None):
    """(taux, borne basse, borne haute, coût médian parmi les succès, n) agrégés sur les problèmes `keys` (liste de clés json)."""
    d = res["data"][method]
    keys = list(d) if keys is None else keys
    ok = np.concatenate([d[k]["ok"] for k in keys]); cost = np.concatenate([d[k]["cost"] for k in keys])
    n, s = len(ok), int(ok.sum())
    lo, hi = wilson(s, n)
    med = float(np.median(cost[ok])) if s else float("nan")
    return dict(rate=s / n, lo=lo, hi=hi, cost=med, n=n)
