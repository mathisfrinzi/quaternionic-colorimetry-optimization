"""Partie 6.1 : points fixes de RNN quaternioniques.

A. Pièce apprise. RNN h ← σ(W·h + U·x + b) avec entrée constante x, W = W₀ (gelée, ℍ-linéaire) + P (apprise, rnn.py) ; on cherche
   le point fixe de chaque entrée apprise depuis des départs perturbés (norme 0,5·√n, comme en partie 4.2). Le défaut ℍ-linéaire
   de J au point fixe est mesuré pour chaque façon d'apprendre P (aucune, ℍ-linéaire, libre, pénalisée).
B. Plusieurs points fixes. Réseau multistable (n = 1 et 2, W ≈ 2,5·rotation faible + bruit, tanh par composante) : depuis 40
   départs, combien de points fixes DISTINCTS chaque méthode trouve-t-elle, et pour quel coût total ? Pour n = 1, chaque point
   trouvé par la sonde est certifié par son degré topologique lu dans les couleurs (partie 1.5, 192 lectures), comparé au signe
   de det J (le degré d'un zéro non dégénéré).
Usage : python bench_points_fixes.py [--smoke] [--A] [--B]"""
import os, sys, json, time
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "src")); sys.path.insert(0, HERE)
from multiprocessing import Pool
from quat import qexp
from sensors import Sensor
import sonde_hn as S
from rnn import RNN, Lblock, hdefect
import degree

P0 = dict(c=0.3, step=0.65, cap=1.5)
METH_A = {"sonde-H": ("sonde", dict(modele="H")), "sonde-Hid": ("sonde", dict(modele="Hid")),
          "sonde-H-relecture": ("sonde", dict(modele="H", verif_k=2, tau=0.25, double=True)),
          "sonde-full": ("sonde", dict(modele="full")),
          "broyden": ("broyden", dict(hfd=1e-6)), "broyden-id": ("broyden-id", {})}
PIECES = {"aucune": dict(piece="aucune"), "H": dict(piece="H"), "libre": dict(piece="libre"),
          "penalise-0.1": dict(piece="penalise", lam=0.1), "penalise-1": dict(piece="penalise", lam=1.0)}


def reseau_A(n, kind, piece, seed):
    return RNN(n, kind, rho0=0.6, K=2, seed=seed, cible=0.3, **PIECES[piece])


def run_A(a):
    meth, n, kind, piece, seed, i, sigma, budget = a
    r = reseau_A(n, kind, piece, seed)
    k = i % 2
    hs = r.point_fixe(k, r.H[k])
    if hs is None:
        return None
    rng = np.random.default_rng(4242 + 31 * i + seed)
    v = rng.standard_normal(4 * n); v /= np.linalg.norm(v)
    q0 = hs + 0.5 * np.sqrt(n) * v
    f = lambda h: r.F(h, k)
    kind_m, opts = METH_A[meth]
    opts = dict(opts)
    if kind_m == "sonde":
        opts.update(P0, cap=P0["cap"] * max(1.0, np.sqrt(n / 4)))
    ok, c, *_ = S.lancer(f, q0, n, kind_m, budget=budget, sigma=sigma, i=i, **opts)
    return dict(meth=meth, n=n, kind=kind, piece=piece, sigma=sigma, seed=seed, i=i, ok=ok, cost=c)


# ------------------------------------------------------------------------------------------------ B : multistable
def reseau_B(n, seed, g=2.5, phi=0.2):
    r = RNN(n, "split", piece="aucune", K=1, seed=seed, beta=0.3)
    rng = np.random.default_rng(seed + 50)
    w = np.zeros((n, n, 4))
    for i in range(n):
        v = rng.standard_normal(3); v /= np.linalg.norm(v)
        w[i, i] = g * qexp(np.concatenate([[0], phi * v]))
    w += 0.3 * rng.standard_normal((n, n, 4)) / np.sqrt(4 * n)
    r.W = Lblock(w)
    return r


def tous_points_fixes(r, n, ntry=400):
    rng = np.random.default_rng(0); pts = []
    for _ in range(ntry):
        h = r.point_fixe(0, rng.uniform(-1.2, 1.2, 4 * n))
        if h is not None and not any(np.linalg.norm(h - p) < 1e-4 for p in pts):
            pts.append(h)
    return pts


def run_B(a):
    meth, n, seed, i, budget = a
    r = reseau_B(n, seed)
    rng = np.random.default_rng(900 + 17 * i + seed)
    q0 = rng.uniform(-1.2, 1.2, 4 * n)
    f = lambda h: r.F(h, 0)
    kind_m, opts = METH_A[meth]
    opts = dict(opts)
    if kind_m == "sonde":
        opts.update(P0)
    ok, c, log, *_ = S.lancer(f, q0, n, kind_m, budget=budget, i=i, log=True, **opts)
    arr = log[-1][2].tolist() if ok else None
    return dict(meth=meth, n=n, seed=seed, i=i, ok=ok, cost=c, arr=arr)


def certifier(r, pts_trouves, pts_vrais, seed):
    """Degré lu dans les couleurs (24 bits, angles, profondeur 2) autour de chaque point fixe trouvé (n = 1)."""
    out = []
    sens = Sensor("s2", 24, "sph", seed=seed)
    for p in pts_trouves:
        p = np.asarray(p)
        # le vrai point fixe le plus proche (le succès s'arrête à ‖F‖ < 10⁻³, pas exactement dessus)
        dv = [np.linalg.norm(p - v) for v in pts_vrais]; j = int(np.argmin(dv)); v = pts_vrais[j]
        others = [np.linalg.norm(v - w) for k, w in enumerate(pts_vrais) if k != j]
        rad = min(0.3, 0.45 * min(others)) if others else 0.3
        d, agree, nl = degree.certificate(lambda h: r.F(h, 0), sens, p, rad, depth=2, seed=seed)
        out.append(dict(deg_lu=int(d), deg_vrai=int(np.sign(np.linalg.det(r.J(v)))), accord=float(agree), lectures=nl,
                        rayon=float(rad), dist=float(dv[j])))
    return out


if __name__ == "__main__":
    smoke = "--smoke" in sys.argv
    doA = "--A" in sys.argv or "--B" not in sys.argv
    doB = "--B" in sys.argv or "--A" not in sys.argv
    out = os.path.join(HERE, "bench_points_fixes_smoke.json" if smoke else "bench_points_fixes_results.json")
    res = json.load(open(out)) if os.path.exists(out) else {}
    t0 = time.time()
    with Pool(2) as pool:
        if doA:
            NS = [4] if smoke else [4, 8]
            seeds, nst = ([0], 4) if smoke else ([0, 1, 2], 10)
            A = []
            for n in NS:
                for kind in ("split", "radial", "global"):
                    for piece in PIECES:
                        dj, dP, em = [], [], []
                        for s in seeds:
                            r = reseau_A(n, kind, piece, s); a, b = r.defauts(); dj.append(a); dP.append(b); em.append(r.erreur_memoire())
                        for sg in ((0.0,) if smoke else (0.0, 0.02)):
                            tasks = [(m, n, kind, piece, s, i, sg, 2000) for m in METH_A for s in seeds for i in range(nst)
                                     if not (m == "sonde-full" and n > 4)]
                            rr = [x for x in pool.map(run_A, tasks, chunksize=2) if x is not None]
                            for x in rr:
                                x.update(defaut_J=float(np.nanmean(dj)), defaut_P=float(np.mean(dP)), err_memoire=float(np.nanmean(em)))
                            A += rr
                            line = "  ".join(f"{m}:{100*np.mean([x['ok'] for x in rr if x['meth']==m]):.0f}%/"
                                             f"{np.median([x['cost'] for x in rr if x['meth']==m and x['ok']] or [np.nan]):.0f}" for m in METH_A
                                             if any(x['meth'] == m for x in rr))
                            print(f"A n={n} {kind:6s} {piece:12s} σ={sg} défJ={np.nanmean(dj):.3f} défP={np.mean(dP):.2f} err={np.nanmean(em):.3f} | {line} [{time.time()-t0:.0f}s]", flush=True)
            res["A"] = A
            json.dump(res, open(out, "w"))
        if doB:
            B = []
            cfg = [(1, s) for s in (0, 1)] + [(2, s) for s in (0, 1)]
            if smoke: cfg = cfg[:1]
            nst = 6 if smoke else 40
            for n, seed in cfg:
                r = reseau_B(n, seed)
                vrais = tous_points_fixes(r, n)
                dfs = [hdefect(r.J(p)) for p in vrais]
                tasks = [(m, n, seed, i, 600) for m in ("sonde-H", "sonde-full", "broyden") for i in range(nst)]
                rr = pool.map(run_B, tasks, chunksize=2)
                for m in ("sonde-H", "sonde-full", "broyden"):
                    xs = [x for x in rr if x["meth"] == m]
                    arr = [np.array(x["arr"]) for x in xs if x["ok"]]
                    # point fixe atteint = le vrai le plus proche
                    ids = sorted({int(np.argmin([np.linalg.norm(a - v) for v in vrais])) for a in arr}) if vrais else []
                    tot = int(sum(x["cost"] for x in xs))
                    cert = []
                    if n == 1 and m == "sonde-H" and arr:
                        uniq = []
                        for a in arr:
                            j = int(np.argmin([np.linalg.norm(a - v) for v in vrais]))
                            if j not in [u[0] for u in uniq]: uniq.append((j, a))
                        cert = certifier(r, [u[1] for u in uniq], vrais, seed)
                    B.append(dict(meth=m, n=n, seed=seed, n_vrais=len(vrais), defaut_moyen=float(np.mean(dfs)),
                                  succes=float(np.mean([x["ok"] for x in xs])), distincts=len(ids), lectures_total=tot,
                                  lectures_med=float(np.median([x["cost"] for x in xs if x["ok"]] or [np.nan])), certificats=cert))
                    cs = f" certificats {sum(c['deg_lu']==c['deg_vrai'] for c in cert)}/{len(cert)}" if cert else ""
                    print(f"B n={n} seed={seed} vrais={len(vrais)} déf={np.mean(dfs):.2f} {m:10s} succès {100*B[-1]['succes']:.0f}% "
                          f"distincts {len(ids)} total {tot}{cs} [{time.time()-t0:.0f}s]", flush=True)
            res["B"] = B
            json.dump(res, open(out, "w"))
