"""Partie 5.7 : un capteur dont l'AMPLITUDE est peu fiable (gain multiplicatif aléatoire à chaque lecture), angle et axe précis.

Un interféromètre ou un polarimètre mesure bien une phase ou une direction et mal une amplitude. On modélise un gain aléatoire par
lecture : ln r reçoit un bruit gaussien d'écart-type σℓ (en népers ; σℓ = 1 ≈ facteur 2,7), alors que l'angle et l'axe ne sont pas
bruités. Les quotients de la sonde éliminent le module d'un gain COMMUN aux deux lectures, mais pas d'un gain indépendant à chaque
lecture ; la question est donc de savoir laquelle des deux familles, couleur ou valeur absolue, souffre le moins.

Même problème ℍ-linéaire que les parties 3.2 / 3.5 / 5.5 (graines 100–102 × 10 départs, départs à 1,5·√n), 24 bits.
Méthodes : sonde-couleur-mémoire-quaternionique, sonde-deux-échelles + test géométrique, Broyden-quaternionique-mêmes-mesures
(luminosité absolue bruitée ; pas de différence finie 0,5 et 1, le meilleur par case) ; bornes sans bruit : Broyden-quaternionique-sans-couleur
(valeur exacte) et Broyden-sans-couleur.
Usage : python bench_gain_amplitude.py [--smoke] [--procs N]   (reprise automatique)"""
import os, sys, json, time
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.join(HERE, "..", "..")
for d in ("src", "parties/3_dimension_4n/3.2_grande_dimension", "parties/3_dimension_4n/3.1_jacobienne_H_lineaire"):
    sys.path.insert(0, os.path.join(ROOT, d))
from multiprocessing import Pool
import sonde_hn as S
import banc_hn as B
import bench_hn_grand as G

P0 = dict(c=0.3, step=0.65, cap=1.5)


def run(a):
    n, sl, seed, i, meth, budget = a
    f, qs = G.problem(n, 0.0, seed); q0 = G.starts(n, qs, i + 1, seed)[i]
    sg = (sl, 0.0, 0.0); cap = P0["cap"] * max(1.0, np.sqrt(n / 4))
    if meth == "sonde":              genre, kw = "sonde", dict(P0, cap=cap, modele="H")
    elif meth == "sonde-2e-geo":     genre, kw = "sonde", dict(P0, cap=cap, modele="H", deux_echelles=True, test="geo")
    elif meth.startswith("bq-deg-"): genre, kw = "broyden-q-deg", dict(side="L", hfd=float(meth.split("-")[-1]))
    elif meth == "bq":               genre, kw = "broyden-q", dict(side="L")
    else:                            genre, kw = "broyden", dict(hfd=1e-6)
    ok, c, *_ = S.lancer(f, q0, n, genre, budget=budget, sigma=sg, i=i, **kw)
    return dict(n=n, sl=sl, seed=seed, i=i, meth=meth, ok=bool(ok), cost=int(c))


METH = ["sonde", "sonde-2e-geo", "bq-deg-0.5", "bq-deg-1", "bq", "broyden"]

if __name__ == "__main__":
    smoke = "--smoke" in sys.argv
    arg = lambda key: sys.argv[sys.argv.index(key) + 1] if key in sys.argv else None
    NS = [4] if smoke else [4, 8, 16]
    SL = [0.0, 0.5] if smoke else [0.0, 0.05, 0.1, 0.2, 0.5, 1.0]
    seeds, nst = ([100], 3) if smoke else ([100, 101, 102], 10)
    out = os.path.join(HERE, "bench_gain_amplitude_smoke.json" if smoke else "bench_gain_amplitude_results.json")
    res = json.load(open(out)) if os.path.exists(out) else []
    have = {(r["n"], r["sl"], r["meth"]) for r in res}
    t0 = time.time()
    with Pool(int(arg("--procs") or 2)) as pool:
        for n in NS:
            for sl in SL:
                todo = [m for m in METH if (n, sl, m) not in have]
                if not todo: continue
                rr = pool.map(run, [(n, sl, s, i, m, 2000) for m in todo for s in seeds for i in range(nst)], chunksize=1)
                res += rr; json.dump(res, open(out, "w"))
                print(f"n={n:2d} σℓ={sl} | " + "  ".join(f"{m}:{B.cellule(B.resume(rr, m))}" for m in todo) + f" [{time.time() - t0:.0f}s]", flush=True)
