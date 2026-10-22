"""Partie 6.5 : pourquoi la sonde s'effondre en régime saturé, alors que le Broyden quaternionique tient.

Constat (parties 5.5 et 6.3) : sur un neurone `tanh` radial (défaut à droite 0,27), la sonde-oracle-décade (qui a EXACTEMENT l'information de la
valeur) réussit 7 % des départs à n = 8, et le Broyden-quaternionique-sans-couleur 100 %. Le défaut de structure est le même pour les deux
modèles (blocs R(b) sur la jacobienne brute), et l'ambiguïté de décade est écartée par l'oracle. Il reste trois candidats dans
l'algorithme de la sonde ; on les retire un à un (la décade est donnée par l'oracle partout) :
  A. le test de descente : « moyenne des ‖ρᵢ‖² < 0,92² » pèse chaque résidu de la même façon, quelle que soit sa taille, alors que la vraie perte
     est ‖F‖² (test pondéré, test="pond") ; le test géométrique (moyenne des ln ‖ρᵢ‖) est l'intermédiaire ;
  B. l'amortissement du pas (step = 0,65) ;
  C. le seuil d'acceptation (acc = 0,92, soit 8 % de baisse exigée) ; acc = 1 : toute baisse est acceptée.
Variantes : base, pond, geo, pas 1,0, acc 1,0, pond + pas 1,0, pond + pas 1,0 + acc 1,0, et pour comparaison le Broyden-quaternionique-sans-couleur.
Neurone `tanh` radial (partie 6.3), n = 4 et 8, 3 instances × 20 départs, budget 2000, σ = 0.
Usage : python bench_diag_sature.py [--smoke] [--procs N] [--seulement v1,v2]"""
import os, sys, json, time
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "src")); sys.path.insert(0, HERE)
from multiprocessing import Pool
import sonde_hn as S
import banc_hn as B
from bench_entrainement import Neurone

P0 = dict(c=0.3, step=0.65, cap=1.5)
VARIANTES = {
    "base":               dict(),
    "pond":               dict(test="pond"),
    "geo":                dict(test="geo"),
    "pas-1":              dict(step=1.0),
    "acc-1":              dict(acc=1.0),
    "pond+pas-1":         dict(test="pond", step=1.0),
    "pond+pas-1+acc-1":   dict(test="pond", step=1.0, acc=1.0),
    # hypothèses écartées (voir README 6.5) : coupure des petites valeurs singulières, taille des sondes après un échec
    "rcond-1e-10":        dict(rcond=1e-10),
    "eta-fixe":           dict(eta_fixe=True),
    "pond+pas-1+rcond+eta-fixe": dict(test="pond", step=1.0, rcond=1e-10, eta_fixe=True),
    # test discriminant : AUCUN test de descente (acc = 10 : tout pas est accepté), pas complet ou amorti
    "aucun-test":         dict(step=1.0, acc=10.0),
    "aucun-test-pas-065": dict(step=0.65, acc=10.0),
}


def run(a):
    n, seed, i, nom, kind = a
    N = Neurone(n, kind, seed); q0 = B.depart(N.qs, n, i, seed, 0.5)
    if nom == "bq":
        ok, c, *_ = S.lancer(N.f, q0, n, "broyden-q", budget=2000, side="R", hfd=1e-6)
    elif nom == "bq-sans-rl":          # Broyden-quaternionique sans recherche linéaire (même absence de test que « aucun-test »)
        ok, c, *_ = S.lancer(N.f, q0, n, "broyden-q", budget=2000, side="R", hfd=1e-6, suivi=True)
    else:
        kw = {**P0, **VARIANTES[nom]}; kw["cap"] = kw["cap"] * max(1.0, np.sqrt(n / 4))
        ok, c, *_ = S.lancer(N.f, q0, n, "sonde", budget=2000, sigma=0.0, i=i, modele="D", lum="desamb", f_vrai=N.f, **kw)
    return dict(n=n, seed=seed, i=i, meth=nom, kind=kind, ok=bool(ok), cost=int(c))


if __name__ == "__main__":
    smoke = "--smoke" in sys.argv
    arg = lambda key: sys.argv[sys.argv.index(key) + 1] if key in sys.argv else None
    NS = [4] if smoke else [4, 8]
    seeds, nst = ([0], 4) if smoke else ([0, 1, 2], 20)
    out = os.path.join(HERE, "bench_diag_sature_smoke.json" if smoke else "bench_diag_sature_results.json")
    seul = arg("--seulement").split(",") if arg("--seulement") else None          # ne lancer que ces variantes (ajoutées au JSON existant)
    res, t0 = ([json.load(open(out))] if (seul and os.path.exists(out)) else [[]])[0], time.time()
    with Pool(int(arg("--procs") or 2)) as pool:
        for kind in ("radial", "radial-faible"):
            for n in NS:
                deja = {(r["kind"], r["n"], r["meth"]) for r in res}          # reprise : ne relance pas ce qui est déjà dans le JSON
                todo = [m for m in (seul or [*VARIANTES, "bq"]) if (kind, n, m) not in deja]
                if not todo: continue
                rr = pool.map(run, [(n, s, i, m, kind) for m in todo for s in seeds for i in range(nst)], chunksize=1)
                res += rr; json.dump(res, open(out, "w"))
                print(f"{kind} n={n} | " + "  ".join(f"{m}:{B.cellule(B.resume(rr, m))}" for m in todo) + f" [{time.time() - t0:.0f}s]", flush=True)
