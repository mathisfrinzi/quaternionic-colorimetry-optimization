"""Partie 3.1 : f : ℍⁿ → ℍⁿ (4n dimensions réelles), n = 1, 2, 3, 4.

Une lecture de f = n couleurs (une par composante quaternionique), donc n·B bits par lecture. Pour chaque composante i le quotient
ρ_i = F_i(q)^{-1} F_i(q+h) vaut 1 + (W h)_i avec W = diag(F)^{-1} J ; l'équation de Newton normalisée est W h = -𝟙 (𝟙 = (1,0,0,0) par bloc).
Acceptation d'un pas : φ = moyenne des |ρ_i|² doit baisser (échelle absolue de chaque F_i inconnue : seul φ se lit dans les couleurs).

Fonctions : f(d) = A d (1 + |d_j|²/2 par bloc) + ε G d, d = q - q*, A ∈ ℍ^{n×n} (multiplication à GAUCHE : jacobienne ℍ-linéaire, 4n² paramètres au lieu de 16n²),
G réelle 4n×4n aléatoire. ε = 0 : presque ℍ-linéaire ; ε = 0.3, 1 : on dégrade vers le générique (même zéro).

Méthodes (S2 sphérique 24 bits par composante) :
  fd-full   : 4n sondes + base par pas, W réelle 4n×4n
  fd-H      : n sondes + base par pas, W = multiplication à gauche par une matrice quaternionique (structure ℍ-linéaire)
  sec-full  : Broyden réel (1 lecture par pas après initialisation par 4n sondes)
  sec-H     : Broyden ℍ-linéaire, w_ij += (d_i - (Wh)_i) conj(h_j)/|h|² (initialisation par n sondes)
  classiques (lectures exactes de F) : broyden, newton-FD, bfgs-FD, + broyden dégradé (même mesure dégradée par composante).
Usage : python bench_hn.py [--smoke]"""
import os, sys, json
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "..", "src"))
from multiprocessing import Pool
from quat import qmul, qinv, qabs, Lmat, rand_unit, E4
from sensors import Sensor, LN10
import classic
from probes import Hit, OutOfBudget

E0 = E4[0]


# ------------------------------------------------------------------------------------------------ problèmes
def problem(n, eps, seed):
    rng = np.random.default_rng(5000 + 97 * seed + 13 * n)
    A = np.zeros((n, n, 4))
    for i in range(n):
        for j in range(n):
            A[i, j] = 0.35 / np.sqrt(n) * rng.standard_normal(4)
        A[i, i] += rand_unit(rng) * 1.0
    G = rng.standard_normal((4 * n, 4 * n)) / np.sqrt(4 * n)
    qs = 0.5 * rng.standard_normal(4 * n)
    def f(q):
        d = (np.asarray(q, float) - qs).reshape(n, 4)
        s = 1.0 + 0.5 * np.sum(d * d, axis=1)
        out = np.zeros((n, 4))
        for i in range(n):
            for j in range(n):
                out[i] += qmul(A[i, j], d[j] * s[j])
        return out.reshape(-1) + eps * (G @ d.reshape(-1))
    return f, qs


def starts(n, qs, k, seed):
    rng = np.random.default_rng(777 + seed)
    out = []
    for _ in range(k):
        v = rng.standard_normal(4 * n); v /= np.linalg.norm(v)
        out.append(qs + 1.5 * np.sqrt(n) * v)
    return out


# ------------------------------------------------------------------------------------------------ compteur
class MeterN:
    """mode 'color' : n lectures Rd ; 'exact' : F ∈ R^{4n} exact ; 'degraded' : F̂ par composante (même résolution/bruit)."""
    def __init__(self, f, n, mode, sens, budget, thr=1e-3):
        self.f, self.n, self.mode, self.sens, self.budget, self.thr, self.cnt = f, n, mode, sens, budget, thr, 0

    def read(self, q):
        self.cnt += 1
        if self.cnt > self.budget:
            raise OutOfBudget
        F = self.f(q)
        if not np.all(np.isfinite(F)):
            raise OutOfBudget
        if np.linalg.norm(F) < self.thr:
            raise Hit
        if self.mode == "exact":
            return F
        Fb = F.reshape(self.n, 4)
        if self.mode == "color":
            return [self.sens[i].read(Fb[i]) for i in range(self.n)]
        return np.concatenate([self.sens[i].absolute(Fb[i]) for i in range(self.n)])


def quot(sens, r0, r1, pred):
    """ρ_i par composante ; amb = au moins une ambiguë (dent de scie)."""
    rho, amb = [], False
    for i, (a, b) in enumerate(zip(r0, r1)):
        r, am = sens[i].quotient(a, b, pred[i])
        rho.append(r); amb |= am
    return np.array(rho), amb


# ------------------------------------------------------------------------------------------------ sondes ℍⁿ
def Lblock(w):
    """w : (n,n,4) coefficients quaternioniques → matrice réelle 4n×4n de h ↦ (Σ_j w_ij h_j)_i."""
    n = w.shape[0]
    W = np.zeros((4 * n, 4 * n))
    for i in range(n):
        for j in range(n):
            W[4 * i:4 * i + 4, 4 * j:4 * j + 4] = Lmat(w[i, j])
    return W


def unit_vec(n):
    return np.tile(E0, n)


def probe_hn(M, q0, model, n, c=0.3, step=0.65, cap=1.5, eta0=0.3, acc2=0.92 ** 2, bt_max=5):
    """model : 'fd-full' | 'fd-H' | 'sec-full' | 'sec-H' | 'sec-Hinit' | 'sec-H2full'.
    Deux choix indépendants : comment la jacobienne est LUE (struct : n sondes, modèle ℍ-linéaire ; sinon 4n sondes) et comment
    elle est MISE À JOUR (upd_struct : Broyden quaternionique qui garde la structure ; sinon Broyden réel, qui peut apprendre la
    partie non ℍ-linéaire).
      sec-Hinit  : « projeter puis affiner » — lecture ℍ-linéaire (n + 1 lectures), mises à jour générales ensuite (partie 9).
      sec-H2full : mises à jour ℍ-linéaires tant que les pas sont acceptés ; au premier pas refusé, on passe définitivement aux
                   mises à jour générales en gardant la W courante (partie 9).
      sec-Hid    : a priori « jacobienne ≈ diagonale par blocs » (le bon a priori pour un point fixe, J = I − Dσ·W, comme Broyden
                   initialisé à l'identité) : W₀ = diag(Fᵢ⁻¹), la direction de Fᵢ lue dans la couleur, sa décade fixée par UNE
                   sonde simultanée η·(1, …, 1) ; 2 lectures au lieu de n + 1, puis Broyden quaternionique (partie 10).""" 
    sens = M.sens; q = np.array(q0, float); dim = 4 * n
    secant = model.startswith("sec")
    struct = model in ("fd-H", "sec-H", "sec-Hinit", "sec-H2full", "sec-Hid", "sec-Hid0")      # lecture
    upd_struct = model in ("sec-H", "sec-H2full", "sec-Hid", "sec-Hid0")                        # mise à jour
    diag_init = model in ("sec-Hid", "sec-Hid0")
    calib0 = model == "sec-Hid0"       # pas de lecture de calibrage : la décade est fixée par le premier pas lui-même
    need_calib = False
    pfloor = np.log(0.05) if calib0 else -np.inf
    r0 = M.read(q); eta = eta0; W = None; fails = 0
    one = unit_vec(n)
    nonlocal_flag = [False]                # sec-Hid0 : vrai tant que la décade de chaque |Fᵢ| n'est pas fixée
    def read_W(q, r0, eta):
        for _ in range(40):
            amb = False
            if calib0:
                # W₀ = diag(Fᵢ⁻¹) avec la décade lue telle quelle (|Fᵢ| supposé dans [1, 10[) : AUCUNE lecture ; le premier pas sert de sonde
                w = np.zeros((n, n, 4))
                for i in range(n):
                    P = np.concatenate([[np.cos(r0[i].th)], r0[i].u * np.sin(r0[i].th)])
                    w[i, i] = P * np.array([1, -1, -1, -1]) / np.exp(r0[i].ell)
                nonlocal_flag[0] = True
                return Lblock(w), eta
            if diag_init:
                # a priori J ≈ I, donc W ≈ diag(F)⁻¹ : la couleur donne la direction de chaque Fᵢ et ln|Fᵢ| à une décade près ;
                # une seule sonde simultanée η·(1, …, 1) donne l'ordre de grandeur |Fᵢ| ≈ η/|ρᵢ − 1| qui fixe la décade.
                e = np.zeros(dim); e[0::4] = eta
                rho, a = quot(sens, r0, M.read(q + e), [0.0] * n)
                if a: eta /= 2; continue
                w = np.zeros((n, n, 4))
                for i in range(n):
                    P = np.concatenate([[np.cos(r0[i].th)], r0[i].u * np.sin(r0[i].th)])          # direction de Fᵢ
                    m = eta / max(np.linalg.norm(rho[i] - E0), 1e-300)                              # ≈ |Fᵢ|
                    k = np.round((np.log(m) - r0[i].ell) / LN10)
                    mod = np.exp(r0[i].ell + k * LN10)
                    w[i, i] = P * np.array([1, -1, -1, -1]) / mod                                   # Fᵢ⁻¹
                return Lblock(w), eta
            if struct:                                   # n sondes : un quaternion unité (1) dans le bloc j
                w = np.zeros((n, n, 4))
                for j in range(n):
                    e = np.zeros(dim); e[4 * j] = eta
                    rho, a = quot(sens, r0, M.read(q + e), [0.0] * n)
                    if a: amb = True; break
                    w[:, j, :] = (rho - E0) / eta
                if amb: eta /= 2; continue
                return Lblock(w), eta
            cols = np.zeros((dim, dim))
            for j in range(dim):
                e = np.zeros(dim); e[j] = eta
                rho, a = quot(sens, r0, M.read(q + e), [0.0] * n)
                if a: amb = True; break
                cols[:, j] = (rho.reshape(-1) - one) / eta
            if amb: eta /= 2; continue
            return cols, eta
        raise OutOfBudget
    while True:
        if W is None:
            W, eta = read_W(q, r0, eta); fails = 0
        t, ok = 1.0, False
        for _ in range(bt_max + 1):
            h = -step * t * np.linalg.lstsq(W, one, rcond=1e-3)[0]
            nh = np.linalg.norm(h)
            if nh > cap: h *= cap / nh
            if nh < 1e-13: raise OutOfBudget
            pred = [max(float(np.log(max(np.linalg.norm((one + W @ h).reshape(n, 4)[i]), 1e-12))), pfloor) for i in range(n)]
            r1 = M.read(q + h)
            rho, amb = quot(sens, r0, r1, pred)
            if nonlocal_flag[0]:
                # le premier pas fait office de sonde : avec J ≈ I, |Fᵢ| ≈ |hᵢ|/|ρᵢ − 1| fixe la décade de chaque composante
                w = np.zeros((n, n, 4)); hb = h.reshape(n, 4)
                for i in range(n):
                    P = np.concatenate([[np.cos(r0[i].th)], r0[i].u * np.sin(r0[i].th)])
                    m = np.linalg.norm(hb[i]) / max(np.linalg.norm(rho[i] - E0), 1e-300)
                    k = np.round((np.log(max(m, 1e-300)) - r0[i].ell) / LN10)
                    w[i, i] = P * np.array([1, -1, -1, -1]) / np.exp(r0[i].ell + k * LN10)
                W = Lblock(w); nonlocal_flag[0] = False
                pred = [max(float(np.log(max(np.linalg.norm((one + W @ h).reshape(n, 4)[i]), 1e-12))), pfloor) for i in range(n)]
                rho, amb = quot(sens, r0, r1, pred)
            D = rho.reshape(-1) - one
            if secant and not amb:
                if upd_struct:
                    wmat = np.zeros((n, n, 4))              # coefficients courants de W
                    for i in range(n):
                        for j in range(n):
                            wmat[i, j] = W[4 * i:4 * i + 4, 4 * j:4 * j + 4][:, 0]    # colonne 0 de L(w) = w
                    Wh = (W @ h).reshape(n, 4); hb = h.reshape(n, 4); Db = D.reshape(n, 4); hn2 = h @ h
                    for i in range(n):
                        for j in range(n):
                            wmat[i, j] += qmul(Db[i] - Wh[i], hb[j] * np.array([1, -1, -1, -1])) / hn2
                    W = Lblock(wmat)
                else:
                    W = W + np.outer(D - W @ h, h) / (h @ h)
            if (not amb) and np.mean(np.sum(rho * rho, axis=1)) < acc2:
                for i in range(n):                          # rebase par composante : W_i ← ρ_i^{-1} W_i
                    Li = Lmat(qinv(rho[i]))
                    W[4 * i:4 * i + 4, :] = Li @ W[4 * i:4 * i + 4, :]
                q, r0, ok = q + h, r1, True
                break
            t *= 0.5
        if not ok:
            fails += 1
            if model == "sec-H2full":
                upd_struct = False                          # la structure ne suffit plus : on affine en général
            if fails >= 2 or not secant:
                W = None; eta = max(eta / 2, 1e-9)
                if eta <= 1e-9: raise OutOfBudget
            continue
        fails = 0
        if not secant:
            W = None
        elif not np.all(np.isfinite(W)):
            W = None


# ------------------------------------------------------------------------------------------------ banc
def run(args):
    name, params, n, eps, seed, i, q0, bits, sigma, budget = args
    f, qs = problem(n, eps, seed)
    mode = "color" if name.startswith(("fd", "sec")) else ("degraded" if name.endswith("-deg") else "exact")
    sens = [Sensor("s2", bits, "sph", sigma=(sigma,) * 3, seed=1000 * i + k) for k in range(n)]
    M = MeterN(f, n, mode, sens, budget)
    try:
        if name.startswith(("fd-", "sec-")):
            probe_hn(M, q0, name, n, **params)
        elif name.startswith("broyden"):
            classic.broyden(M, q0, **params)
        elif name == "newton-FD":
            classic.newton_fd(M, q0, **params)
        elif name == "bfgs-FD":
            classic.bfgs_fd(M, q0, **params)
    except Hit:
        return True, M.cnt
    except (OutOfBudget, np.linalg.LinAlgError, FloatingPointError, ValueError, ZeroDivisionError, OverflowError):
        return False, budget
    return False, budget


METHODS = {"fd-full": [dict(c=0.3, step=0.65, cap=1.5), dict(c=0.3, step=0.8, cap=3.0)],
           "fd-H": [dict(c=0.3, step=0.65, cap=1.5), dict(c=0.3, step=0.8, cap=3.0)],
           "sec-full": [dict(c=0.3, step=0.65, cap=1.5), dict(c=0.3, step=0.8, cap=3.0)],
           "sec-H": [dict(c=0.3, step=0.65, cap=1.5), dict(c=0.3, step=0.8, cap=3.0)],
           "broyden": [dict(hfd=1e-6), dict(hfd=1e-4)], "newton-FD": [dict(hfd=1e-6), dict(hfd=1e-4)], "bfgs-FD": [dict(hfd=1e-6)],
           "broyden-deg": [dict(hfd=1e-2), dict(hfd=1e-1)]}


def evaluate(name, params, n, eps, bits, sigma, nseeds, nstart, budget, pool, tune=False):
    tasks = []
    for seed in range(nseeds):
        f, qs = problem(n, eps, seed + (100 if tune else 0))
        for i, q0 in enumerate(starts(n, qs, nstart, seed + (100 if tune else 0))):
            tasks.append((name, params, n, eps, seed + (100 if tune else 0), i, q0, bits, sigma, budget))
    res = pool.map(run, tasks, chunksize=4)
    return np.array([r[0] for r in res]), np.array([r[1] for r in res], float)


if __name__ == "__main__":
    smoke = "--smoke" in sys.argv
    NS = [1, 2] if smoke else [1, 2, 3, 4]
    EPSS = [0.0, 1.0] if smoke else [0.0, 0.3, 1.0]
    CASES = [(24, 0.0)] if smoke else [(24, 0.0), (24, 0.02)]
    nseeds, nstart, budget = (1, 4, 300) if smoke else (3, 20, 600)
    out = []
    with Pool(2) as pool:
        for n in NS:
            for eps in EPSS:
                for bits, sg in CASES:
                    for name, cfgs in METHODS.items():
                        if name.endswith("-deg") and sg == 0.0 and False:
                            continue
                        best = None
                        for p in cfgs:                                 # réglage sur d'autres instances/départs
                            ok, c = evaluate(name, p, n, eps, bits, sg, 1, 8 if not smoke else 4, budget, pool, tune=True)
                            key = (ok.mean(), -np.median(np.where(ok, c, budget)))
                            if best is None or key > best[0]: best = (key, p)
                        ok, c = evaluate(name, best[1], n, eps, bits, sg, nseeds, nstart, budget, pool)
                        out.append(dict(n=n, eps=eps, bits=bits, sigma=sg, method=name, params=best[1], rate=float(ok.mean()),
                                        cost=float(np.median(c[ok])) if ok.any() else None, ok=ok.tolist(), costs=c.tolist()))
                        print(n, eps, bits, sg, f"{name:12s} {100*ok.mean():4.0f}%  {out[-1]['cost']}", flush=True)
    json.dump(out, open(os.path.join(HERE, "bench_hn_smoke.json" if smoke else "bench_hn_results.json"), "w"))
