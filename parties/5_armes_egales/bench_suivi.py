"""Partie 5.8 : suivre un zéro qui bouge (optimisation en temps réel), à nombre de lectures par pas de temps donné.

Le monde avance pendant qu'on lit : la fonction vue à la k-ième lecture est f(·, t) avec t = k / r, où r est le nombre de lectures
par pas de temps (r = 1, 2, 4, 8). Aucune méthode n'est modifiée : elle voit simplement une fonction qui change (aucun arrêt au succès).
  F(q, t) = A(t)·(q − q*(t)),  A(t) = A₀ + 0,3·sin(0,15 t)·B   (A₀, B ∈ ℍ^{n×n} comme en 3.2, donc ℍ-linéaire à tout instant),
  q*(t) = q*(0) + v·t·e   (dérive rectiligne à vitesse v par pas de temps, e unitaire).
Départ à 1,5·√n de q*(0). 80 pas de temps. Mesure : à chaque pas, la meilleure lecture (distance ‖q − q*(t)‖ du point lu au zéro du
moment), puis la médiane de ces distances sur les 40 derniers pas ; on la divise par la vitesse v : c'est un retard, en pas de temps
(un suivi parfait avec une lecture de retard donnerait 1/r).
Méthodes, toutes en variante « suivi » (pas de Newton complet toujours accepté, aucun test de descente ni recherche linéaire : ils
comparent des valeurs prises à des instants différents et refusent presque tous les pas) : sonde-couleur-mémoire-quaternionique,
sonde-deux-échelles (capteur à bits égaux, la décade se lit), sonde-oracle-décade (borne), Broyden-quaternionique-sans-couleur (valeur
exacte), Broyden-quaternionique-mêmes-mesures (luminosité absolue de plage illimitée) et à bits égaux (capteur à deux échelles lu comme une
valeur), tous avec un pas de différence finie de 0,5. Contrôle : « bq-sans-suivi », le Broyden-quaternionique-sans-couleur de la partie 5.5
inchangé (avec recherche linéaire).
n = 4 et 8, σ = 0 et 0,02, v = 0,02 ; 0,1 ; 0,3 ; 3 matrices × 5 départs.
Usage : python bench_suivi.py [--smoke] [--procs N]"""
import os, sys, json, time
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.join(HERE, "..", "..")
sys.path.insert(0, os.path.join(ROOT, "src"))
from multiprocessing import Pool
import sonde_hn as S
import banc_hn as B

P0 = dict(c=0.3, step=0.65, cap=1.5)
T = 80


class Mobile:
    def __init__(self, n, seed, v, r):
        rng = np.random.default_rng(4100 + 53 * seed + n)
        self.n, self.v, self.r, self.k = n, v, r, 0
        def mat():
            A = np.zeros((n, n, 4))
            for i in range(n):
                for j in range(n):
                    A[i, j] = 0.35 / np.sqrt(n) * rng.standard_normal(4)
                A[i, i, 0] += 1.0
            return S.Lblock(A)
        self.A0, self.B = mat(), mat() - mat()
        self.q0 = 0.5 * rng.standard_normal(4 * n)
        e = rng.standard_normal(4 * n); self.e = e / np.linalg.norm(e)

    def cible(self, t):
        return self.q0 + self.v * t * self.e

    def pur(self, q):
        t = self.k / self.r
        return (self.A0 + 0.3 * np.sin(0.15 * t) * self.B) @ (np.asarray(q, float) - self.cible(t))

    def f(self, q):
        F = self.pur(q); self.k += 1
        return F


def run(a):
    n, v, r, seed, i, meth, sigma = a
    W = Mobile(n, seed, v, r); budget = T * r
    rng = np.random.default_rng(555 + 17 * i + seed)
    d = rng.standard_normal(4 * n); q0 = W.cible(0) + 1.5 * np.sqrt(n) * d / np.linalg.norm(d)
    cap = P0["cap"] * max(1.0, np.sqrt(n / 4))
    exact = meth in ("bq", "bq-sans-suivi")
    if exact:
        mode, sens = "exact", None
    elif meth == "sonde-2e":
        mode, sens = "color", S.capteurs(n, sigma, i, True)
    elif meth == "bq-deg2":
        mode, sens = "degraded2", S.capteurs(n, sigma, i, True)
    elif meth.startswith("sonde"):
        mode, sens = "color", S.capteurs(n, sigma, i)
    else:
        mode, sens = "degraded", S.capteurs(n, sigma, i)
    M = S.MeterLog(W.f, n, mode, sens, budget, thr=-1.0, log=True)
    try:
        # variantes « suivi » : aucun test de descente (la fonction change entre deux lectures), pas de Newton complet.
        # Pas de différence finie de Broyden-quaternionique : 0,5 pour toutes (un pas de 10⁻⁶ serait noyé par le mouvement du zéro).
        suivi = dict(P0, cap=cap, step=1.0, acc=10.0)
        if meth == "sonde":
            S.sonde(M, q0, n, modele="H", **suivi)
        elif meth == "sonde-2e":
            S.sonde(M, q0, n, modele="H", **suivi)
        elif meth == "sonde-oracle":
            S.sonde(M, q0, n, modele="H", lum="desamb", f_vrai=W.pur, **suivi)
        elif meth in ("bq", "bq-deg", "bq-deg2"):
            S.broyden_q(M, q0, n, side="L", hfd=0.5, suivi=True)
        else:                                         # sans suivi : la méthode de la partie 5.5 telle quelle
            S.broyden_q(M, q0, n, side="L", hfd=1e-6)
    except S.ERRS:
        pass
    # distance au zéro du moment, à chaque lecture (la lecture numéro c a été faite à t = (c − 1) / r)
    dist = [(c, float(np.linalg.norm(q - W.cible((c - 1) / r)))) for c, nf, q in M.log]
    par_pas = {}
    for c, dd in dist:
        s_ = (c - 1) // r; par_pas[s_] = min(par_pas.get(s_, np.inf), dd)
    fin = [par_pas.get(s_, np.inf) for s_ in range(T // 2, T)]
    err = float(np.median(fin))
    return dict(n=n, v=v, r=r, seed=seed, i=i, meth=meth, sigma=sigma, err=err if np.isfinite(err) else None,
                retard=None if not np.isfinite(err) else err / v, lus=len(M.log))


METH = ["sonde", "sonde-2e", "sonde-oracle", "bq", "bq-deg", "bq-deg2", "bq-sans-suivi"]

if __name__ == "__main__":
    smoke = "--smoke" in sys.argv
    arg = lambda key: sys.argv[sys.argv.index(key) + 1] if key in sys.argv else None
    NS = [4] if smoke else [4, 8]
    VS = [0.1] if smoke else [0.02, 0.1, 0.3]
    RS = [1, 4] if smoke else [1, 2, 4, 8]
    SIG = [0.0] if smoke else [0.0, 0.02]
    seeds, nst = ([0], 2) if smoke else ([0, 1, 2], 5)
    out = os.path.join(HERE, "bench_suivi_smoke.json" if smoke else "bench_suivi_results.json")
    res = json.load(open(out)) if os.path.exists(out) else []
    have = {(r["n"], r["v"], r["r"], r["sigma"]) for r in res}
    t0 = time.time()
    with Pool(int(arg("--procs") or 2)) as pool:
        for n in NS:
            for sg in SIG:
                for v in VS:
                    for r in RS:
                        if (n, v, r, sg) in have: continue
                        rr = pool.map(run, [(n, v, r, s, i, m, sg) for m in METH for s in seeds for i in range(nst)], chunksize=1)
                        res += rr; json.dump(res, open(out, "w"))
                        line = "  ".join(f"{m}:{np.median([x['retard'] for x in rr if x['meth'] == m and x['retard'] is not None] or [np.nan]):.2f}" for m in METH)
                        print(f"n={n} σ={sg} v={v} r={r} | retard médian (pas) {line} [{time.time() - t0:.0f}s]", flush=True)
