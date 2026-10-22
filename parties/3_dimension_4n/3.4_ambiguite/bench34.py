"""Partie 3.4 : banc comparatif des remèdes à l'ambiguïté de la luminosité.

Même famille de fonctions que les parties 3.2 et 3.3 (ε = 0), couleur-complète-angles 24 bits par composante, budget 2000.
Mêmes départs pour toutes les méthodes (comparaisons appariées).

Usage :
  python bench34.py --smoke                 quelques départs, n = 16 et 32, pour vérifier et chronométrer
  python bench34.py --calib                 choix de k et τ sur des graines de CALIBRATION (300–302), jamais réutilisées
  python bench34.py --eval [--ns 16,32] [--sigmas 0,0.02] [--meth a,b] [--nseeds 10] [--k 10]
                                            banc final, graines 100–109 × 10 départs = 100 départs par case
Chaque case est ajoutée à son JSON dès qu'elle est finie ; une case déjà présente est sautée (reprise possible)."""
import os, sys, json, time
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
D33 = os.path.join(HERE, "..", "3.3_decrochage"); D31 = os.path.join(HERE, "..", "3.1_jacobienne_H_lineaire"); D32 = os.path.join(HERE, "..", "3.2_grande_dimension")
for p in (HERE, os.path.join(HERE, "..", "..", "..", "src"), D31, D32):
    sys.path.insert(0, p)
from multiprocessing import Pool
from sensors import Sensor
from probes import Hit, OutOfBudget
import classic
from bench_hn import MeterN
import bench_hn_grand as G
from sonde34 import probe34

P0 = dict(c=0.3, step=0.65, cap=1.5)
ERRS = (OutOfBudget, np.linalg.LinAlgError, FloatingPointError, ValueError, ZeroDivisionError, OverflowError)

# noms du rapport → options de probe34 ; "broyden" = Broyden-sans-couleur (lit F exactement, hfd = 1e-6 comme en 3.2)
METHODES = {
    "sonde-garde-ambiguite":  dict(),
    "sonde-depliage-vers-zero": dict(lum="prior0"),
    "sonde-oracle-decade":    dict(lum="desamb"),
    "broyden-sans-couleur":   "broyden",
    # réglage figé après calibration (bench34_calib.json, graines 300–304) : les configurations à double dépliage sont à égalité
    # dans le bruit de 30 départs ; on garde la plus simple (k = 2, τ = 0,25, sans repli ni mise à jour supplémentaire)
    "R1-relecture-naive":     dict(verif_k=4, tau=0.25),
    "R1-relecture":           dict(verif_k=2, tau=0.25, double=True),
    "R2-par-composante":      dict(partiel=True),
    "R1+R2":                  dict(verif_k=2, tau=0.25, double=True, partiel=True),
}


def run(args):
    meth, opts, n, sigma, seed, i, budget = args
    f, qs = G.problem(n, 0.0, seed)
    q0 = G.starts(n, qs, i + 1, seed)[i]
    t0 = time.time()
    st = {}
    if opts == "broyden":
        M = MeterN(f, n, "exact", None, budget)
        try:
            classic.broyden(M, q0, hfd=1e-6); ok = False
        except Hit:
            ok = True
        except ERRS:
            ok = False
    else:
        sens = [Sensor("s2", 24, "sph", sigma=(sigma,) * 3, seed=1000 * i + k) for k in range(n)]
        M = MeterN(f, n, "color", sens, budget)
        try:
            probe34(M, q0, n, **dict(P0, cap=P0["cap"] * max(1.0, np.sqrt(n / 4))), **opts, f_vrai=f, stats=st)
            ok = False
        except Hit:
            ok = True
        except ERRS:
            ok = False
    return dict(meth=meth, n=n, sigma=sigma, seed=seed, i=i, ok=ok, cost=M.cnt if ok else budget,
                reads=min(M.cnt, budget), sec=time.time() - t0, **st)


def case(pool, meth, opts, n, sigma, seeds, k, budget):
    tasks = [(meth, opts, n, sigma, s, i, budget) for s in seeds for i in range(k)]
    return pool.map(run, tasks, chunksize=1)


def summary(res):
    ok = np.array([r["ok"] for r in res]); c = np.array([r["cost"] for r in res], float)
    pas = sum(r.get("pas", 0) for r in res) or 1
    return (f"succès {100*ok.mean():3.0f}%  lect méd {np.median(c[ok]) if ok.any() else float('nan'):6.1f}"
            f"  rej amb {sum(r.get('rej_amb',0) for r in res)/pas:.2f} desc {sum(r.get('rej_desc',0) for r in res)/pas:.2f}"
            f"  verif/pas {sum(r.get('verif',0) for r in res)/pas:.2f}  [{sum(r['sec'] for r in res):.0f} s cpu]")


def arg(name, default, cast=str):
    return cast(sys.argv[sys.argv.index(name) + 1]) if name in sys.argv else default


if __name__ == "__main__":
    budget = 2000
    if "--smoke" in sys.argv:
        NS, SIGS, seeds, k, out = [16, 32], [0.0, 0.02], [100], 2, "bench34_smoke.json"
        meths = {m: METHODES[m] for m in arg("--meth", ",".join(METHODES)).split(",")}
    elif "--calib" in sys.argv:
        NS, SIGS, seeds, k, out = [32], [0.0, 0.02], [300, 301, 302, 303, 304], 6, "bench34_calib.json"
        meths = {}
        # le smoke test et les 4 premières configurations (double = 0, gardées dans le JSON) montrent que sans double prior
        # les décades fausses muettes passent : on ne balaie que double = 1
        for kk in (2, 4):
            for tau in (0.25, 0.4):
                for repli in (False, True):
                    meths[f"R1 k={kk} tau={tau} double=1 repli={int(repli)} maj=0"] = dict(verif_k=kk, tau=tau, double=True, repli=repli)
        meths["R1 k=2 tau=0.4 double=1 repli=1 maj=1"] = dict(verif_k=2, tau=0.4, double=True, repli=True, verif_maj=True)
    else:
        NS = [int(x) for x in arg("--ns", "16,32").split(",")]
        SIGS = [float(x) for x in arg("--sigmas", "0,0.02").split(",")]
        seeds = list(range(100, 100 + arg("--nseeds", 10, int))); k = arg("--k", 10, int)
        out = arg("--out", "bench34_results.json")
        meths = {m: METHODES[m] for m in arg("--meth", ",".join(METHODES)).split(",")}
    path = os.path.join(HERE, out)
    done = json.load(open(path)) if os.path.exists(path) else []
    have = {(r["meth"], r["n"], r["sigma"]) for r in done}
    t0 = time.time()
    with Pool(2) as pool:
        for n in NS:
            for sg in SIGS:
                for meth, opts in meths.items():
                    if (meth, n, sg) in have:
                        continue
                    if opts == "broyden" and sg > 0 and (meth, n, 0.0) in have:
                        # Broyden-sans-couleur lit F exactement : σ ne le concerne pas, on recopie la case σ = 0
                        res = [dict(r, sigma=sg) for r in done if r["meth"] == meth and r["n"] == n and r["sigma"] == 0.0]
                    else:
                        res = case(pool, meth, opts, n, sg, seeds, k, budget)
                    done += res
                    json.dump(done, open(path, "w"))
                    print(f"n={n:2d} σ={sg:<4} {meth:26s} {summary(res)}  [{time.time()-t0:.0f} s]", flush=True)
