"""Partie 6.2 : inverser un réseau quaternionique en changeant l'activation.

f(x) = V·σ(U·x + b) − y*, même tirage que la famille « inversion réseau » de la partie 4.1, avec trois activations (partie 4.2) :
tanh par composante (celle de la partie 4.1, où toutes les sondes échouent), tanh radial, saturation globale.
Méthodes : sonde-couleur-mémoire-quaternionique (H), sonde-couleur-mémoire (full, n ≤ 4), Broyden-sans-couleur, Broyden-mêmes-mesures.
Couleur-complète-angles 24 bits par composante, budget 600 (n ≤ 4) ou 2000 (n = 8), départs à 1,0·√n du zéro (comme en 4.1).
Usage : python bench_inversion.py [--smoke]"""
import os, sys, json, time
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "src")); sys.path.insert(0, HERE)
from multiprocessing import Pool
import banc_hn as B
from rnn import QNN, hdefect

for _k in ("split", "radial", "global"):
    B.PROBLEMES["qnn-" + _k] = (lambda kind: (lambda n, seed: (lambda net: (net.f, net.qs, {}))(QNN(n, kind, seed))))(_k)
METH = {"sonde-H": ("sonde", dict(modele="H")), "sonde-full": ("sonde", dict(modele="full")),
        "broyden": ("broyden", {}), "broyden-deg": ("broyden-deg", {})}

if __name__ == "__main__":
    smoke = "--smoke" in sys.argv
    NS = [2] if smoke else [2, 3, 4, 8]
    seeds, nst = ([0], 4) if smoke else ([0, 1, 2], 20)
    SIG = [0.0] if smoke else [0.0, 0.02]
    out = os.path.join(HERE, "bench_inversion_smoke.json" if smoke else "bench_inversion_results.json")
    res, defauts = [], {}
    t0 = time.time()
    with Pool(2) as pool:
        for kind in ("split", "radial", "global"):
            prob = "qnn-" + kind
            for n in NS:
                defauts[f"{kind}/{n}"] = [hdefect(QNN(n, kind, s).J(QNN(n, kind, s).qs)) for s in seeds]
                budget = 2000 if n > 4 else 600
                for sg in SIG:
                    tasks = [(prob, n, s, i, m, mo, sg, budget, 1.0) for m, mo in METH.items() for s in seeds for i in range(nst)
                             if not (m == "sonde-full" and n > 4)]
                    rr = pool.map(B.run, tasks, chunksize=2)
                    for r in rr: r["kind"] = kind
                    res += rr
                    line = "  ".join(f"{m}:{B.cellule(B.resume(rr, m))}" for m in METH if B.resume(rr, m))
                    print(f"{kind:6s} n={n} σ={sg} déf={np.mean(defauts[f'{kind}/{n}']):.2f} | {line} [{time.time()-t0:.0f}s]", flush=True)
                json.dump(dict(results=res, defauts=defauts), open(out, "w"))
