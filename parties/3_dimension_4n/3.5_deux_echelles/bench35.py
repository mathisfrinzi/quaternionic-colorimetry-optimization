"""Partie 3.5, suite de 3.4 : les pistes restantes contre l'ambiguïté de la luminosité, puis le test de descente bruité.

Mêmes problèmes, mêmes départs (graines 100–109 × 10 = 100 départs par case), mêmes capteurs que la partie 3.4 : les méthodes de
la 3.4 (garde-ambiguïté, dépliage-vers-zéro, relecture, oracle-décade, Broyden-sans-couleur) ne sont PAS relancées, leurs
résultats sont relus dans ../3.4_ambiguite/bench34_results.json (src/sonde_hn.py reproduit exactement probe34 : tests/test_sonde_hn.py).

Nouvelles méthodes (src/sonde_hn.py) :
  sonde-relecture-économe     relit seulement si le double dépliage est en désaccord (pas sur le seul écart à la prédiction)
  sonde-vote                  décalage commun des composantes sûres (médiane) appliqué aux composantes douteuses
  sonde-deux-échelles         capteur à deux échelles À BITS ÉGAUX (24 bits : luminosité fine 5, angle 5, axe 5 + 5,
                              luminosité grossière 4 bits sur 10 décades) : la décade est lue, l'ambiguïté disparaît à la source
  …+ test géométrique         test de descente sur la moyenne des ln‖ρᵢ‖ au lieu de la moyenne des ‖ρᵢ‖²
  …+ test pondéré             test sur Σ mᵢ²‖ρᵢ‖² / Σ mᵢ² avec mᵢ = ‖Fᵢ‖ suivi par les luminosités (la vraie perte ‖F‖²)
  sonde-oracle-décade + test géométrique / pondéré   (TRICHE, borne) : ce que vaut le test seul, ambiguïté levée
Diagnostic (vraie fonction, sans effet sur l'algorithme) : part des pas refusés par le test de descente alors que la vraie
moyenne des ‖ρᵢ‖² passait (« refus à tort »).
Usage : python bench35.py [--smoke] [--ns 16,32] [--sigmas 0,0.02] [--meth a,b]   (reprise possible : les cases faites sont sautées)"""
import os, sys, json, time
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
R = os.path.join(HERE, "..", "..", "..")
sys.path.insert(0, os.path.join(R, "src")); sys.path.insert(0, os.path.join(HERE, "..", "3.2_grande_dimension")); sys.path.insert(0, os.path.join(HERE, "..", "3.1_jacobienne_H_lineaire"))
from multiprocessing import Pool
import sonde_hn as S
import bench_hn_grand as G

P0 = dict(c=0.3, step=0.65, cap=1.5)
METHODES = {
    "sonde-relecture-econome": dict(verif_k=2, econome=True),
    "sonde-vote": dict(vote=True),
    "sonde-deux-echelles": dict(deux_echelles=True),
    "sonde-deux-echelles-geo": dict(deux_echelles=True, test="geo"),
    "sonde-deux-echelles-pond": dict(deux_echelles=True, test="pond"),
    "sonde-oracle-decade": dict(lum="desamb"),          # relancée pour le diagnostic « refus à tort » (et contrôle : mêmes chiffres qu'en 3.4)
    "sonde-oracle-decade-geo": dict(lum="desamb", test="geo"),
    "sonde-oracle-decade-pond": dict(lum="desamb", test="pond"),
}


def run(a):
    meth, opts, n, sigma, seed, i, budget = a
    f, qs = G.problem(n, 0.0, seed); q0 = G.starts(n, qs, i + 1, seed)[i]
    t0 = time.time()
    ok, c, _, _, st = S.lancer(f, q0, n, "sonde", budget=budget, sigma=sigma, i=i, f_vrai=f,
                               **dict(P0, cap=P0["cap"] * max(1.0, np.sqrt(n / 4))), **opts)
    return dict(meth=meth, n=n, sigma=sigma, seed=seed, i=i, ok=bool(ok), cost=int(c), sec=time.time() - t0, **st)


def arg(name, default, cast=str):
    return cast(sys.argv[sys.argv.index(name) + 1]) if name in sys.argv else default


if __name__ == "__main__":
    smoke = "--smoke" in sys.argv
    NS = [int(x) for x in arg("--ns", "16,32").split(",")]
    SIGS = [float(x) for x in arg("--sigmas", "0,0.02").split(",")]
    meths = {m: METHODES[m] for m in arg("--meth", ",".join(METHODES)).split(",")}
    seeds, k = ([100], 2) if smoke else (list(range(100, 110)), 10)
    path = os.path.join(HERE, "bench35_smoke.json" if smoke else "bench35_results.json")
    done = json.load(open(path)) if os.path.exists(path) else []
    have = {(r["meth"], r["n"], r["sigma"]) for r in done}
    t0 = time.time()
    with Pool(2) as pool:
        for n in NS:
            for sg in SIGS:
                for meth, opts in meths.items():
                    if (meth, n, sg) in have:
                        continue
                    res = pool.map(run, [(meth, opts, n, sg, s, i, 2000) for s in seeds for i in range(k)], chunksize=1)
                    done += res
                    json.dump(done, open(path, "w"))
                    ok = np.array([r["ok"] for r in res]); c = np.array([r["cost"] for r in res])
                    pas = sum(r["pas"] for r in res) or 1
                    print(f"n={n:2d} σ={sg:<4} {meth:26s} succès {100*ok.mean():3.0f}%  lect méd {np.median(c[ok]) if ok.any() else float('nan'):6.1f}"
                          f"  refus ambiguïté {sum(r['rej_amb'] for r in res)/pas:.2f} descente {sum(r['rej_desc'] for r in res)/pas:.2f}"
                          f" (dont à tort {sum(r['rej_tort'] for r in res)/max(sum(r['rej_desc'] for r in res),1):.2f})  [{time.time()-t0:.0f} s]", flush=True)
