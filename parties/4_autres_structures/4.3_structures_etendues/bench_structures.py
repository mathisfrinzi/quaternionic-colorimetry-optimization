"""Partie 4.3 : étendre la structure quand elle manque.

La partie 4.1 a montré que « projeter puis affiner » ne rattrape pas une structure fausse. On teste trois modèles intermédiaires
entre 4n² (ℍ-linéaire) et 16n² paramètres (réel quelconque), chacun sur la famille pour laquelle il est fait (sonde_hn.py) :

  modèle (nom du rapport)                        famille (zéro q*, d = q − q*)                                         paramètres
  sonde-…-conforme-par-blocs   (blocconf)        bilatère : fᵢ = Σⱼ aᵢⱼ·dⱼ·bᵢⱼ · (1 + ‖d‖²/2)  (partie 4.1, défaut ≈ 0,8–0,9)      7n²
  sonde-…-quaternionique+rang-k (Hrang, k = 2)   ℍ-linéaire + rang 2 : A·d·(1 + ‖dⱼ‖²/2) + t·P·Qᵀ·d, P, Q réelles 4n×2      4n² + 16n
  sonde-…-boîte-grise          (grise)           inversion d'un réseau, tanh par composante (partie 4.1) : J = L(V)·Dσ·L(U),  4n² + 3n
                                                  U connue, Dσ diagonale réelle inconnue

Pour chaque famille : la sonde ℍ-linéaire (sonde-couleur-mémoire-quaternionique), le modèle étendu, la sonde sans structure
(sonde-couleur-mémoire), Broyden-sans-couleur, Broyden-mêmes-mesures. Couleur-complète-angles 24 bits par composante, budget 600.
Usage : python bench_structures.py [--smoke]"""
import os, sys, json, time
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
R = os.path.join(HERE, "..", "..", "..")
for p in (os.path.join(R, "src"), os.path.join(HERE, "..", "4.1_hors_structure"), os.path.join(R, "parties", "6_reseaux")):
    sys.path.insert(0, p)
from multiprocessing import Pool
import banc_hn as B
from bench_hors_structure import problem as prob41, hdefect
from rnn import QNN, Lblock

TP = 0.6       # taille relative de la pièce de rang 2


def p_bilat(n, seed):
    f, qs = prob41("bilat", n, seed); return f, qs, {}


def p_hrang(n, seed):
    f0, qs = prob41("hlin", n, seed)
    rng = np.random.default_rng(4300 + 31 * seed + n)
    P, Q = rng.standard_normal((2, 4 * n, 2)) / np.sqrt(4 * n)
    LA = np.linalg.norm(np.column_stack([(f0(qs + 1e-6 * e) - f0(qs - 1e-6 * e)) / 2e-6 for e in np.eye(4 * n)]))
    PQ = P @ Q.T; PQ *= TP * LA / np.linalg.norm(PQ)
    f = lambda q: f0(q) + PQ @ (np.asarray(q, float) - qs)
    return f, qs, {}


def p_qnn(n, seed):
    net = QNN(n, "split", seed)
    return net.f, net.qs, dict(U=net.U)


B.PROBLEMES.update(bilat=p_bilat, hrang=p_hrang, qnn=p_qnn)
RAD = dict(bilat=1.5, hrang=1.5, qnn=1.0)
COMMUNES = {"sonde-H": ("sonde", dict(modele="H")), "sonde-full": ("sonde", dict(modele="full")),
            "broyden": ("broyden", {}), "broyden-deg": ("broyden-deg", {})}
ETENDUE = {"bilat": {"sonde-blocconf": ("sonde", dict(modele="blocconf"))},
           "hrang": {"sonde-Hrang": ("sonde", dict(modele="Hrang", k=2))},
           "qnn": {"sonde-grise": ("sonde", dict(modele="grise"))}}


def taches(prob, n, seeds, nst, sigma, budget):
    meths = {**ETENDUE[prob], **COMMUNES}
    out = []
    for s in seeds:
        extra = B.PROBLEMES[prob](n, s)[2]
        for m, (g, o) in meths.items():
            o = dict(o, **(extra if o.get("modele") == "grise" else {}))
            for i in range(nst):
                out.append((prob, n, s, i, m, (g, o), sigma, budget, RAD[prob]))
    return out


if __name__ == "__main__":
    smoke = "--smoke" in sys.argv
    NS = [2] if smoke else [2, 3, 4]
    seeds, nst = ([0], 4) if smoke else ([0, 1, 2], 20)
    SIG = [0.0] if smoke else [0.0, 0.02]
    out = os.path.join(HERE, "bench_structures_smoke.json" if smoke else "bench_structures_results.json")
    res, defauts = [], {}
    t0 = time.time()
    with Pool(2) as pool:
        for prob in ("bilat", "hrang", "qnn"):
            for n in NS:
                defauts[f"{prob}/{n}"] = [hdefect(*B.PROBLEMES[prob](n, s)[:2]) for s in seeds]
                for sg in SIG:
                    rr = pool.map(B.run, taches(prob, n, seeds, nst, sg, 600), chunksize=2)
                    res += rr
                    line = "  ".join(f"{m}:{B.cellule(B.resume(rr, m))}" for m in dict.fromkeys(r["meth"] for r in rr))
                    print(f"{prob} n={n} σ={sg} déf={np.mean(defauts[f'{prob}/{n}']):.2f} | {line} [{time.time()-t0:.0f}s]", flush=True)
                json.dump(dict(results=res, defauts=defauts), open(out, "w"))
