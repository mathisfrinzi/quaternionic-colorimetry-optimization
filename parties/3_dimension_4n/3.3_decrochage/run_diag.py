"""Partie 3.3 : balayage des diagnostics. Usage : python run_diag.py [--smoke]"""
import os, sys, json, time
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE); sys.path.insert(0, os.path.join(HERE, "..", "..", "..", "src")); sys.path.insert(0, os.path.join(HERE, "..", "3.1_jacobienne_H_lineaire")); sys.path.insert(0, os.path.join(HERE, "..", "3.2_grande_dimension"))
from multiprocessing import Pool
from sensors import Sensor
from probes import Hit, OutOfBudget
from bench_hn import MeterN
import bench_hn_grand as G
from diag11 import probe_diag

# chaque variante : options de probe_diag + capmult (facteur sur le plafond de pas) + bits du capteur + sigma
VARIANTES = {
    "base":        dict(),
    "lum_exacte":  dict(lum="exact"),
    "W_vraie":     dict(W_oracle=True),
    "pas/2":       dict(capmult=0.5),
    "pas/4":       dict(capmult=0.25),
    "36_bits":     dict(bits=36),
    "48_bits":     dict(bits=48),
    "desamb":      dict(lum="desamb"),
    "amorti0.3":   dict(amorti=0.3),
    "amorti1":     dict(amorti=1.0),
    "amorti3":     dict(amorti=3.0),
    "prior0":      dict(lum="prior0"),
    "noflag":      dict(lum="noflag"),
    "desamb+am1":  dict(lum="desamb", amorti=1.0),
}
P0 = dict(c=0.3, step=0.65, cap=1.5)


def run(args):
    var, n, eps, sigma, bits, seed, i, budget = args
    opts = dict(VARIANTES[var]); capmult = opts.pop("capmult", 1.0); bits = opts.pop("bits", bits)
    f, qs = G.problem(n, eps, seed)
    q0 = G.starts(n, qs, i + 1, seed)[i]
    sens = [Sensor("s2", bits, "sph", sigma=(sigma,) * 3, seed=1000 * i + k) for k in range(n)]
    M = MeterN(f, n, "color", sens, budget)
    log = []
    ok = False
    try:
        p = dict(P0, cap=capmult * P0["cap"] * max(1.0, np.sqrt(n / 4)))
        probe_diag(M, q0, n, **p, **opts, log=log, f_vrai=f)
    except Hit:
        ok = True
    except (OutOfBudget, np.linalg.LinAlgError, FloatingPointError, ValueError, ZeroDivisionError, OverflowError):
        pass
    pas = [e for e in log if e["kind"] == "pas"]; son = [e for e in log if e["kind"] == "sonde"]
    nc = max(len(pas), 1)
    w = [e["werr"] for e in pas if e["werr"] is not None]
    return dict(ok=ok, cost=M.cnt if ok else budget, n_pas=len(pas), n_sonde=len(son),
                f_amb=sum(e["cause"] == "ambigu" for e in pas) / nc,
                f_desc=sum(e["cause"] == "descente" for e in pas) / nc,
                f_acc=sum(e["cause"] == "accepte" for e in pas) / nc,
                comp_amb=float(np.mean([e["n_amb"] / n for e in pas])) if pas else 0.0,
                nh_med=float(np.median([e["nh"] for e in pas])) if pas else 0.0,
                perr=float(np.mean([e["perr"] for e in pas])) if pas else 0.0,
                werr_med=float(np.median(w)) if w else None,
                sonde_rejet=sum(e["ok"] is False for e in son))


if __name__ == "__main__":
    smoke = "--smoke" in sys.argv
    SIG = float(sys.argv[sys.argv.index("--sigma") + 1]) if "--sigma" in sys.argv else 0.0
    NS = [8] if smoke else ([int(x) for x in sys.argv[sys.argv.index("--ns") + 1].split(",")] if "--ns" in sys.argv else [8, 16, 32])
    ONLY = sys.argv[sys.argv.index("--var") + 1].split(",") if "--var" in sys.argv else list(VARIANTES)
    seeds, k = ([100], 2) if smoke else ([100, 101, 102], 4)
    budget = 2000
    out = []
    t0 = time.time()
    with Pool(2) as pool:
        for n in NS:
            for var in ONLY:
                tasks = [(var, n, 0.0, SIG, 24, s, i, budget) for s in seeds for i in range(k)]
                res = pool.map(run, tasks, chunksize=1)
                ok = np.array([r["ok"] for r in res])
                cost = np.array([r["cost"] for r in res], float)
                agg = {key: float(np.mean([r[key] for r in res])) for key in ("f_amb", "f_desc", "f_acc", "comp_amb", "nh_med", "perr", "sonde_rejet")}
                wl = [r["werr_med"] for r in res if r["werr_med"] is not None]
                agg["werr_med"] = float(np.median(wl)) if wl else None
                row = dict(n=n, variante=var, bits=24, sigma=SIG, rate=float(ok.mean()),
                           cost=float(np.median(cost[ok])) if ok.any() else None, **agg)
                out.append(row)
                print(f"n={n:2d} {var:11s} succès {100*ok.mean():3.0f}%  lect {row['cost']}  rejets: amb {agg['f_amb']:.2f} desc {agg['f_desc']:.2f}"
                      f"  comp_amb {agg['comp_amb']:.3f}  |h| {agg['nh_med']:.2f}  werr {agg['werr_med']}  [{time.time()-t0:.0f}s]", flush=True)
                json.dump(out, open(os.path.join(HERE, "diag_smoke.json" if smoke else f"diag_{'+'.join(ONLY)}_sigma{SIG}.json" if "--var" in sys.argv else f"diag_results_sigma{SIG}.json"), "w"), indent=1)
