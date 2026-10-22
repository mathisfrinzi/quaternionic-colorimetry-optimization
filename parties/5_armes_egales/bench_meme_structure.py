"""Partie 5.5 : MÊME STRUCTURE pour tous. La couleur contre la valeur, à structure égale.

Dans les parties 3, 4 et 6, la sonde exploite une structure (blocs L(a) ℍ-linéaires, blocs R(b) à droite, diagonale par sortie)
que les méthodes classiques ignoraient : Broyden-sans-couleur lit 4n directions. Mais la structure ne demande pas la couleur :
si un bloc vaut L(a) ou R(b), une différence finie le long de la direction réelle 1 donne a ou b, avec la valeur exacte aussi.
La référence qui manquait est donc un Broyden QUATERNIONIQUE classique (src/sonde_hn.broyden_q) : même modèle de jacobienne que
la sonde, mêmes mises à jour qui gardent la structure, mais il lit la valeur de f.

Méthodes nouvelles (lancées ici) :
  bq           Broyden-quaternionique-sans-couleur : valeur exacte, même structure que la sonde
  bq-deg-0.5   Broyden-quaternionique-mêmes-mesures, pas de différence finie 0,5 (mesure dégradée de la partie 5 : mêmes bits,
  bq-deg-1     même bruit, luminosité avec sa décade) ; deux pas, réglés sur d'autres graines ; on garde le meilleur par case
  sonde-oracle sonde avec la décade donnée (lum="desamb") : la sonde qui a EXACTEMENT l'information de la mesure dégradée
  (pour ℍ-linéaire n = 4, 8 : aussi la sonde de base, la sonde à deux échelles + test géométrique et Broyden-sans-couleur)
Méthodes relues (mêmes départs, comparaisons appariées) : parties 3.4 / 3.5 (ℍ-linéaire n = 16, 32), 6.3 (neurone), 6.4 (couche).

Familles : "H" (partie 3.2, n = 4, 8, 16, 32 ; graines 100–109 × 10 départs), "neurone" (6.3 ; n = 4, 8, 16 ; 3 × 20), "couche"
(6.4 ; m × p = 4×4, 8×4, 4×8, 8×8 ; 3 × 10). σ = 0 et 0,02, 24 bits.
Usage : python bench_meme_structure.py [--smoke] [--fams H,neurone,couche] [--procs N]   (reprise automatique)"""
import os, sys, json, time
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.join(HERE, "..", "..")
for d in ("src", "parties/6_reseaux", "parties/3_dimension_4n/3.2_grande_dimension", "parties/3_dimension_4n/3.1_jacobienne_H_lineaire"):
    sys.path.insert(0, os.path.join(ROOT, d))
from multiprocessing import Pool
import sonde_hn as S
import banc_hn as B
import bench_hn_grand as G
from bench_entrainement import Neurone
from bench_couche import Couche

P0 = dict(c=0.3, step=0.65, cap=1.5)


def probleme(fam, kind, n, m, seed, i):
    """(f, q0, côté, m) avec exactement les départs des parties d'origine."""
    if fam == "H":
        f, qs = G.problem(n, 0.0, seed); return f, G.starts(n, qs, i + 1, seed)[i], "L", 1
    if fam == "neurone":
        N = Neurone(n, kind, seed); return N.f, B.depart(N.qs, n, i, seed, 0.5), "R", 1
    L = Couche(m, n // m, kind, seed); return L.f, B.depart(L.qs, n, i, seed, 0.5), "R", m


def run(a):
    fam, kind, n, m, seed, i, meth, sigma, budget = a
    f, q0, side, mm = probleme(fam, kind, n, m, seed, i)
    if meth.startswith("bq"):
        genre = "broyden-q" if meth == "bq" else "broyden-q-deg"
        kw = dict(side=side, m=mm) if meth == "bq" else dict(side=side, m=mm, hfd=float(meth.split("-")[-1]))
    elif meth == "broyden":
        genre, kw = "broyden", dict(hfd=1e-6)
    else:
        modele = "H" if fam == "H" else ("D" if fam == "neurone" else "Dbloc")
        extra = dict(sonde=dict(), **{"sonde-oracle": dict(lum="desamb"), "sonde-2e-geo": dict(deux_echelles=True, test="geo")})[meth]
        genre = "sonde"
        kw = dict(P0, cap=P0["cap"] * max(1.0, np.sqrt(n / 4)), modele=modele, **extra)
        if modele == "Dbloc": kw["m"] = mm
        kw["f_vrai"] = f          # utilisé seulement par l'oracle-décade (et les diagnostics)
    ok, c, *_ = S.lancer(f, q0, n, genre, budget=budget, sigma=sigma, i=i, **kw)
    return dict(fam=fam, kind=kind, n=n, m=m, seed=seed, i=i, meth=meth, sigma=sigma, ok=bool(ok), cost=int(c))


NOUV = ["bq", "bq-deg-0.5", "bq-deg-1", "sonde-oracle"]


def cases(smoke):
    out = []   # (fam, kind, n, m, seeds, nst, budget, méthodes)
    if smoke:
        return [("H", "", 4, 1, [100], 3, 2000, NOUV), ("neurone", "lin", 4, 1, [0], 3, 600, NOUV),
                ("couche", "lin", 16, 4, [0], 3, 2000, NOUV)]
    for n in (4, 8, 16, 32):
        extra = ["sonde", "sonde-2e-geo", "broyden"] if n <= 8 else []
        out.append(("H", "", n, 1, list(range(100, 110)), 10, 2000, NOUV + extra))
    for kind in ("lin", "radial-faible", "radial"):
        for n in (4, 8, 16):
            out.append(("neurone", kind, n, 1, [0, 1, 2], 20, 2000 if n > 8 else 600, NOUV))
        for m, p in ((4, 4), (8, 4), (4, 8), (8, 8)):
            # couche en tanh radial : la sonde-oracle est déjà à 0-47 % pour un neurone (n ≤ 16) et ne finirait pas en temps raisonnable
            out.append(("couche", kind, m * p, m, [0, 1, 2], 10, 2000, [x for x in NOUV if not (kind == "radial" and x == "sonde-oracle")]))
    return out


if __name__ == "__main__":
    smoke = "--smoke" in sys.argv
    arg = lambda k: sys.argv[sys.argv.index(k) + 1] if k in sys.argv else None
    fams = arg("--fams").split(",") if arg("--fams") else ["H", "neurone", "couche"]
    out = os.path.join(HERE, "bench_meme_structure_smoke.json" if smoke else "bench_meme_structure_results.json")
    res = json.load(open(out)) if os.path.exists(out) else []
    have = {(r["fam"], r["kind"], r["n"], r["m"], r["sigma"], r["meth"]) for r in res}
    t0 = time.time()
    with Pool(int(arg("--procs") or 2)) as pool:
        for fam, kind, n, m, seeds, nst, budget, meths in cases(smoke):
            if fam not in fams: continue
            for sg in ([0.0] if smoke else [0.0, 0.02]):
                todo = [x for x in meths if (fam, kind, n, m, sg, x) not in have]
                if not todo: continue
                rr = pool.map(run, [(fam, kind, n, m, s, i, x, sg, budget) for x in todo for s in seeds for i in range(nst)], chunksize=1)
                res += rr
                json.dump(res, open(out, "w"))
                line = "  ".join(f"{x}:{B.cellule(B.resume(rr, x))}" for x in todo)
                print(f"{fam:7s} {kind:13s} n={n:2d} m={m} σ={sg} | {line} [{time.time() - t0:.0f}s]", flush=True)
