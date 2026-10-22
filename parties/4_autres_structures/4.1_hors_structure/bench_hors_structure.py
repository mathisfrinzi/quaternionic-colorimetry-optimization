"""Partie 4.1 : ℍⁿ → ℍⁿ HORS de la structure construite exprès. Le gain de la partie 8 survit-il sur des fonctions qui ne sont
pas fabriquées pour avoir une jacobienne ℍ-linéaire ? Et sinon, « projeter sur les ℍ-linéaires puis affiner » le rattrape-t-il ?

Familles (zéro connu en q*, d = q - q*, coefficients quaternioniques aléatoires) :
  hlin     f_i = Σ_j a_ij d_j (1 + |d_j|²/2)                  témoin : la famille de la partie 8 (ε = 0)
  bilat    f_i = Σ_j a_ij d_j b_ij                           produits des DEUX côtés : chaque bloc est conforme, la somme ne l'est pas
  poly     f_i = Σ_j a_ij d_j + Σ_jk c_ijk d_j d_k            polynôme non commutatif : ℍ-linéaire AU zéro, de moins en moins en s'éloignant
  qnn      f(x) = V tanh(U x + b) - y*,  y* = même réseau en x*   inversion d'un petit réseau quaternionique (tanh composante par composante)
  qres     f(h) = h - tanh(W h + b)                            point fixe d'un réservoir quaternionique (W ℍ-linéaire, rayon spectral 0,8)
  gen      f = G d (1 + |d|²/2), G réelle                      témoin générique
Pour chaque instance, le DÉFAUT ℍ-linéaire de la jacobienne au zéro, ‖J − P(J)‖/‖J‖ (P : projection orthogonale sur les matrices
dont chaque bloc 4×4 est une multiplication à gauche), est mesuré : 0 = ℍ-linéaire, ≈ 0,87 pour une matrice réelle quelconque.
Sondes (S2 sphérique 24 bits par composante) : sec-full, sec-H, sec-Hinit (lecture ℍ-linéaire puis Broyden réel), sec-H2full
(Broyden quaternionique, puis réel au premier échec). Références : broyden (lectures exactes), broyden-deg (mêmes mesures).
Usage : python bench_hors_structure.py [--smoke]"""
import os, sys, json
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "..", "src")); sys.path.insert(0, os.path.join(HERE, "..", "..", "3_dimension_4n", "3.1_jacobienne_H_lineaire"))
from multiprocessing import Pool
from quat import qmul, rand_unit, Lmat, BASIS
from sensors import Sensor
import classic
from probes import Hit, OutOfBudget
from bench_hn import MeterN, probe_hn

FAMS = ["hlin", "bilat", "poly", "qnn", "qres", "gen"]


def qmat(rng, n, diag=1.0, off=0.35):
    A = off / np.sqrt(n) * rng.standard_normal((n, n, 4))
    for i in range(n):
        A[i, i] += diag * rand_unit(rng)
    return A


def apply(A, x):
    """(A x)_i = Σ_j A_ij x_j (produit de Hamilton, coefficient à gauche)."""
    n = A.shape[0]; out = np.zeros((n, 4))
    for i in range(n):
        for j in range(n):
            out[i] += qmul(A[i, j], x[j])
    return out


def problem(fam, n, seed):
    rng = np.random.default_rng(9000 + 131 * seed + 17 * n + 7 * FAMS.index(fam))
    qs = 0.5 * rng.standard_normal(4 * n)
    if fam == "hlin":
        A = qmat(rng, n)
        def g(d):
            return apply(A, d * (1 + 0.5 * np.sum(d * d, axis=1))[:, None])
    elif fam == "bilat":
        A = qmat(rng, n); Bq = np.array([[rand_unit(rng) for _ in range(n)] for _ in range(n)])
        def g(d):
            out = np.zeros((n, 4))
            for i in range(n):
                for j in range(n):
                    out[i] += qmul(qmul(A[i, j], d[j]), Bq[i, j])
            return out * (1 + 0.5 * np.sum(d * d))
    elif fam == "poly":
        A = qmat(rng, n); C = 0.5 / n * rng.standard_normal((n, n, n, 4))
        def g(d):
            out = apply(A, d)
            for i in range(n):
                for j in range(n):
                    for k in range(n):
                        out[i] += qmul(qmul(C[i, j, k], d[j]), d[k])
            return out
    elif fam == "qnn":
        U = qmat(rng, n, 1.0, 0.6); V = qmat(rng, n, 1.0, 0.6); b = 0.8 * rng.standard_normal((n, 4))
        xs = qs.reshape(n, 4)
        net = lambda x: apply(V, np.tanh(apply(U, x) + b))
        ys = net(xs)
        def g(d):
            return net(xs + d) - ys
    elif fam == "qres":
        Wr = qmat(rng, n, 0.0, 1.0)
        Wr *= 0.8 / max(abs(np.linalg.eigvals(_L(Wr))))
        b = 0.6 * rng.standard_normal((n, 4))
        hs = np.zeros((n, 4))                                    # point fixe par itération (contraction)
        for _ in range(2000):
            hs = np.tanh(apply(Wr, hs) + b)
        qs = hs.reshape(-1)
        def g(d):
            h = hs + d
            return h - np.tanh(apply(Wr, h) + b)
    elif fam == "gen":
        G = rng.standard_normal((4 * n, 4 * n)) / np.sqrt(4 * n) + np.eye(4 * n)
        def g(d):
            return (G @ d.reshape(-1)).reshape(n, 4) * (1 + 0.5 * np.sum(d * d))
    else:
        raise ValueError(fam)
    f = lambda q: g((np.asarray(q, float) - qs).reshape(n, 4)).reshape(-1)
    return f, qs


def _L(A):
    n = A.shape[0]; W = np.zeros((4 * n, 4 * n))
    for i in range(n):
        for j in range(n):
            W[4 * i:4 * i + 4, 4 * j:4 * j + 4] = Lmat(A[i, j])
    return W


def hdefect(f, qs, h=1e-6):
    """Défaut ℍ-linéaire de la jacobienne en qs : projection de chaque bloc sur {L(w)} (base L(e_k) orthogonale, norme² 4)."""
    m = len(qs); J = np.zeros((m, m))
    for j in range(m):
        e = np.zeros(m); e[j] = h
        J[:, j] = (f(qs + e) - f(qs - e)) / (2 * h)
    P = np.zeros_like(J); n = m // 4
    Lb = [Lmat(np.eye(4)[k]) for k in range(4)]
    for i in range(n):
        for j in range(n):
            B = J[4 * i:4 * i + 4, 4 * j:4 * j + 4]
            P[4 * i:4 * i + 4, 4 * j:4 * j + 4] = sum(np.sum(B * Lb[k]) / 4 * Lb[k] for k in range(4))
    return float(np.linalg.norm(J - P) / np.linalg.norm(J))


def starts(n, qs, k, seed, rad):
    rng = np.random.default_rng(777 + seed)
    out = []
    for _ in range(k):
        v = rng.standard_normal(4 * n); v /= np.linalg.norm(v)
        out.append(qs + rad * np.sqrt(n) * v)
    return out


RAD = dict(hlin=1.5, bilat=1.5, poly=1.0, qnn=1.0, qres=1.0, gen=1.5)


def run(args):
    name, params, fam, n, seed, i, q0, bits, sigma, budget = args
    f, qs = problem(fam, n, seed)
    mode = "color" if name.startswith("sec-") else ("degraded" if name.endswith("-deg") else "exact")
    sens = [Sensor("s2", bits, "sph", sigma=(sigma,) * 3, seed=1000 * i + k) for k in range(n)]
    M = MeterN(f, n, mode, sens, budget)
    try:
        if name.startswith("sec-"):
            probe_hn(M, q0, name, n, **params)
        else:
            classic.broyden(M, q0, **params)
    except Hit:
        return True, M.cnt
    except (OutOfBudget, np.linalg.LinAlgError, FloatingPointError, ValueError, ZeroDivisionError, OverflowError):
        return False, budget
    return False, budget


PCFG = [dict(c=0.3, step=0.65, cap=1.5), dict(c=0.3, step=0.8, cap=3.0)]
METHODS = {"sec-full": PCFG, "sec-H": PCFG, "sec-Hinit": PCFG, "sec-H2full": PCFG,
           "broyden": [dict(hfd=1e-6), dict(hfd=1e-4)], "broyden-deg": [dict(hfd=1e-2), dict(hfd=1e-1)]}


def evaluate(name, params, fam, n, bits, sigma, seeds, nstart, budget, pool):
    tasks = []
    for seed in seeds:
        f, qs = problem(fam, n, seed)
        for i, q0 in enumerate(starts(n, qs, nstart, seed, RAD[fam])):
            tasks.append((name, params, fam, n, seed, i, q0, bits, sigma, budget))
    res = pool.map(run, tasks, chunksize=4)
    return np.array([r[0] for r in res]), np.array([r[1] for r in res], float)


if __name__ == "__main__":
    smoke = "--smoke" in sys.argv
    NS = [2] if smoke else [2, 3, 4]
    CASES = [(24, 0.0)] if smoke else [(24, 0.0), (24, 0.02)]
    seeds, nstart, budget = ([0], 4, 400) if smoke else ([0, 1, 2], 20, 600)
    out, defects = [], {}
    with Pool(2) as pool:
        for fam in FAMS:
            for n in NS:
                defects[f"{fam}/{n}"] = [hdefect(*problem(fam, n, s)) for s in seeds]
                print(fam, n, "défaut ℍ-linéaire au zéro :", np.round(defects[f"{fam}/{n}"], 3), flush=True)
                for bits, sg in CASES:
                    for name, cfgs in METHODS.items():
                        best = None
                        for p in cfgs:
                            ok, c = evaluate(name, p, fam, n, bits, sg, [100], 4 if smoke else 8, budget, pool)
                            key = (ok.mean(), -np.median(np.where(ok, c, budget)))
                            if best is None or key > best[0]: best = (key, p)
                        ok, c = evaluate(name, best[1], fam, n, bits, sg, seeds, nstart, budget, pool)
                        out.append(dict(fam=fam, n=n, bits=bits, sigma=sg, method=name, params=best[1], rate=float(ok.mean()),
                                        cost=float(np.median(c[ok])) if ok.any() else None, ok=ok.tolist(), costs=c.tolist()))
                        print(f"  {fam} n={n} σ={sg} {name:12s} {100*ok.mean():4.0f}%  {out[-1]['cost']}", flush=True)
                json.dump(dict(results=out, defects=defects), open(os.path.join(HERE, "bench_smoke.json" if smoke else "bench_hors_structure_results.json"), "w"))
