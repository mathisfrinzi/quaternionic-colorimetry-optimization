"""Partie 5.6 : une boîte noire ℍ-linéaire INSTABLE (mal conditionnée), à bits égaux.

Question : la dent de scie de la luminosité (précision relative constante à toute échelle, plage illimitée par continuité) donne-t-elle
un avantage quand F traverse beaucoup de décades ? La référence « mêmes mesures » de la partie 5.5 lisait une luminosité absolue de
plage illimitée sans payer de bits : c'était un cadeau. Ici la mesure absolue est PAYÉE : c'est le capteur à deux échelles (24 bits :
luminosité fine 5 + angle 5 + axe 5 + 5 + luminosité grossière 4 sur une fenêtre de 10 décades), lu comme une valeur. Hors de la
fenêtre, la luminosité grossière est écrêtée.

Deux familles, toutes deux avec A ∈ ℍ^{n×n}, A = U·diag(10^{k·sᵢ})·V†, U, V quaternioniques unitaires (Gram–Schmidt), sᵢ sur [0, 1] :
  « cond »     f(q) = A·(q − q*), valeurs singulières de 1 à 10^k (conditionnement 10^k). Défaut ℍ-linéaire nul pour tout k.
  « explosif » k = 0 (A unitaire) et f(q) = A·(dⱼ·(1 + ½‖dⱼ‖²)^p)ⱼ avec d = q − q* : bien conditionné au zéro (jacobienne A), mais ‖f‖ croît
               comme ‖d‖^{2p+1} : plusieurs décades parcourues pendant la descente, comme une dynamique qui explose.
Succès : ‖f‖ < 10⁻³ (donc ‖q − q*‖ < 10⁻³ puisque la plus petite valeur singulière vaut 1).
Fenêtre (−4, 6) de la mesure absolue à bits égaux : elle contient le seuil de succès ; au-dessus de 10⁶, la luminosité grossière est écrêtée.

Méthodes : sonde-couleur-mémoire-quaternionique (luminosité en dent de scie, plage illimitée), sonde-oracle-décade (la sonde avec la décade
donnée : borne), sonde-deux-échelles + test géométrique (même capteur à bits égaux que la mesure absolue, fenêtre finie),
Broyden-quaternionique-sans-couleur (valeur exacte : borne), Broyden-quaternionique-mêmes-mesures (luminosité absolue de plage
illimitée, 24 bits, comme en 5.5) et Broyden-quaternionique-mesure-absolue-à-bits-égaux (capteur à deux échelles, fenêtre finie),
Broyden-sans-couleur (sans structure). Pas de différence finie de Broyden-quaternionique-mêmes-mesures : 0,5 et 1 (le meilleur par case).
Usage : python bench_boite_instable.py [--smoke] [--procs N]   (reprise automatique)"""
import os, sys, json, time
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.join(HERE, "..", "..")
sys.path.insert(0, os.path.join(ROOT, "src"))
from multiprocessing import Pool
from quat import qmul, qconj
import sonde_hn as S
import banc_hn as B

P0 = dict(c=0.3, step=0.65, cap=1.5)
FEN = (-4.0, 6.0)


def qunitaire(n, rng):
    """Matrice quaternionique unitaire n×n par Gram–Schmidt : colonnes u_j, <u,v> = Σ conj(u_k)·v_k, projection u·<u,v> (scalaire à droite)."""
    V = rng.standard_normal((n, n, 4)); U = np.zeros_like(V)
    for j in range(n):
        v = V[:, j].copy()
        for i in range(j):
            c = sum(qmul(qconj(U[k, i]), v[k]) for k in range(n))
            for k in range(n):
                v[k] = v[k] - qmul(U[k, i], c)
        U[:, j] = v / np.sqrt(np.sum(v * v))
    return U


def matrice_A(n, k, seed):
    rng = np.random.default_rng(8800 + 31 * seed + n)
    U, V = qunitaire(n, rng), qunitaire(n, rng)
    d = 10.0 ** (k * np.linspace(0.0, 1.0, n)) if n > 1 else np.ones(1)
    # A = U·diag(d)·V†  : (U d V†)_{ij} = Σ_l U_il d_l conj(V_jl)
    A = np.zeros((n, n, 4))
    for i in range(n):
        for j in range(n):
            A[i, j] = sum(d[l] * qmul(U[i, l], qconj(V[j, l])) for l in range(n))
    return A


def probleme(n, k, p, seed):
    A = matrice_A(n, k, seed); LA = S.Lblock(A)
    rng = np.random.default_rng(9900 + 17 * seed + n)
    qs = 0.5 * rng.standard_normal(4 * n)

    def f(q):
        d = (np.asarray(q, float) - qs).reshape(n, 4)
        if p:
            d = d * ((1.0 + 0.5 * np.sum(d * d, axis=1)) ** p)[:, None]
        return LA @ d.reshape(-1)
    return f, qs, LA


def depart(qs, n, i, seed):
    rng = np.random.default_rng(777 + seed + 13 * i)
    v = rng.standard_normal(4 * n); v /= np.linalg.norm(v)
    return qs + 1.5 * np.sqrt(n) * v


def run(a):
    n, k, p, seed, i, meth, sigma, budget = a
    f, qs, _ = probleme(n, k, p, seed); q0 = depart(qs, n, i, seed)
    cap = P0["cap"] * max(1.0, np.sqrt(n / 4))
    if meth == "sonde":                genre, kw = "sonde", dict(P0, cap=cap, modele="H")
    elif meth == "sonde-oracle":       genre, kw = "sonde", dict(P0, cap=cap, modele="H", lum="desamb", f_vrai=f)
    elif meth == "sonde-2e-geo":       genre, kw = "sonde", dict(P0, cap=cap, modele="H", deux_echelles=True, test="geo", fenetre=FEN)
    elif meth == "bq":                 genre, kw = "broyden-q", dict(side="L")
    elif meth.startswith("bq-deg2"):   genre, kw = "broyden-q-deg2", dict(side="L", hfd=float(meth.split("-")[-1]), fenetre=FEN)
    elif meth.startswith("bq-deg-"):   genre, kw = "broyden-q-deg", dict(side="L", hfd=float(meth.split("-")[-1]))
    else:                              genre, kw = "broyden", dict(hfd=1e-6)
    ok, c, *_ = S.lancer(f, q0, n, genre, budget=budget, sigma=sigma, i=i, **kw)
    return dict(n=n, k=k, p=p, seed=seed, i=i, meth=meth, sigma=sigma, ok=bool(ok), cost=int(c))


METH = ["sonde", "sonde-oracle", "sonde-2e-geo", "bq", "bq-deg-0.5", "bq-deg-1", "bq-deg2-0.5", "bq-deg2-1", "broyden"]

if __name__ == "__main__":
    smoke = "--smoke" in sys.argv
    arg = lambda key: sys.argv[sys.argv.index(key) + 1] if key in sys.argv else None
    NS = [4] if smoke else [4, 8, 16]
    CAS = [(0, 0), (6, 0), (0, 4)] if smoke else [(0, 0), (1, 0), (2, 0), (3, 0), (0, 2), (0, 4), (0, 8)]     # (k, p)
    seeds, nst = ([100], 3) if smoke else ([100, 101, 102], 10)
    SIG = [0.0] if smoke else [0.0, 0.02]
    out = os.path.join(HERE, "bench_boite_instable_smoke.json" if smoke else "bench_boite_instable_results.json")
    res = json.load(open(out)) if os.path.exists(out) else []
    have = {(r["n"], r["k"], r["p"], r["sigma"], r["meth"]) for r in res}
    t0 = time.time()
    with Pool(int(arg("--procs") or 2)) as pool:
        for n in NS:
            for sg in SIG:
                if n == 16 and sg > 0: continue          # temps de calcul : à n = 16, sans bruit seulement
                for k, p in CAS:
                    todo = [m for m in METH if (n, k, p, sg, m) not in have]
                    if not todo: continue
                    rr = pool.map(run, [(n, k, p, s, i, m, sg, 3000) for m in todo for s in seeds for i in range(nst)], chunksize=1)
                    res += rr; json.dump(res, open(out, "w"))
                    line = "  ".join(f"{m}:{B.cellule(B.resume(rr, m))}" for m in todo)
                    print(f"n={n:2d} k={k:2d} p={p} σ={sg} | {line} [{time.time() - t0:.0f}s]", flush=True)
