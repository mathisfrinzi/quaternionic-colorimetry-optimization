"""Méthodes classiques de référence, sur le même compteur de lectures que les sondes (`Meter`).

Une « lecture » = une évaluation de F (un quaternion = 4 réels lus ensemble). Conventions de coût, comme dans le repo complexe :
  * différences finies avant : 4 lectures par jacobienne 4x4 (+ la lecture au point, déjà connue en général) ;
  * `newton_exactJ` : la jacobienne exacte est facturée 1 évaluation (borne basse, très favorable à Newton) ;
  * tout appel à F, y compris dans une recherche linéaire ou dans scipy, est compté.
Toutes lèvent `Hit` (succès) ou `OutOfBudget` via le Meter, ou retournent en cas d'échec propre.
"""
import numpy as np
from scipy.optimize import minimize, least_squares
from quat import E4, jac_num
from probes import Hit, OutOfBudget

nrm = np.linalg.norm


def fd_jac(M, q, F0, h):
    """Jacobienne par différences finies avant, en dimension len(q) (4 pour ℍ, 4n pour ℍⁿ : partie 8)."""
    n = len(q)
    J = np.empty((len(F0), n))
    for j in range(n):
        e = np.zeros(n); e[j] = 1.0
        J[:, j] = (M.read(q + h * e) - F0) / h
    return J


def _line(M, q, F, d, bt_max):
    """Recherche linéaire par divisions par 2 : accepte le premier pas qui baisse |F|."""
    t = 1.0
    for _ in range(bt_max + 1):
        Fn = M.read(q + t * d)
        if nrm(Fn) < nrm(F):
            return t, Fn
        t /= 2
    return None, None


def newton_fd(M, q0, hfd=1e-5, bt_max=6):
    """Newton réel 4D, jacobienne par différences finies avant (4 lectures par pas)."""
    q = np.array(q0, float); F = M.read(q)
    while True:
        J = fd_jac(M, q, F, hfd)
        d = -np.linalg.lstsq(J, F, rcond=1e-8)[0]
        t, Fn = _line(M, q, F, d, bt_max)
        if t is None:
            raise OutOfBudget
        q, F = q + t * d, Fn


def newton_exactJ(M, q0, bt_max=6):
    """Newton réel 4D avec jacobienne exacte facturée 1 évaluation (borne basse du coût de Newton)."""
    q = np.array(q0, float); F = M.read(q)
    while True:
        M.charge(1)
        J = jac_num(M.f, q)
        d = -np.linalg.lstsq(J, F, rcond=1e-8)[0]
        t, Fn = _line(M, q, F, d, bt_max)
        if t is None:
            raise OutOfBudget
        q, F = q + t * d, Fn


def broyden(M, q0, hfd=1e-5, bt_max=6, restart=2):
    """Broyden « bon » (mise à jour de l'inverse), jacobienne initiale par différences finies, 1 lecture par pas ensuite."""
    q = np.array(q0, float); F = M.read(q)
    H, fails = None, 0
    while True:
        if H is None:
            H = np.linalg.pinv(fd_jac(M, q, F, hfd), rcond=1e-8); fails = 0
        d = -H @ F
        t, Fn = _line(M, q, F, d, bt_max)
        if t is None:
            fails += 1
            if fails >= restart:
                H = None
            else:
                H = None
            continue
        s, y = t * d, Fn - F
        Hy = H @ y; den = s @ Hy
        if abs(den) > 1e-14:
            H = H + np.outer(s - Hy, s @ H) / den
        q, F = q + s, Fn
        fails = 0


def bfgs_fd(M, q0, hfd=1e-6):
    """BFGS (scipy) sur |F|^2, gradient 2 J^T F par différences finies : 4 lectures par gradient."""
    cache = {}
    def fun(x):
        F = M.read(x); cache[x.tobytes()] = F
        return float(F @ F)
    def jac(x):
        F0 = cache.get(x.tobytes())
        if F0 is None:
            F0 = M.read(x)
        return 2 * fd_jac(M, x, F0, hfd).T @ F0
    minimize(fun, np.array(q0, float), jac=jac, method="BFGS", options=dict(gtol=1e-14, maxiter=10 ** 6))
    raise OutOfBudget


def lm_fd(M, q0, diff_step=1e-6):
    """Levenberg-Marquardt (MINPACK) avec jacobienne par différences finies."""
    least_squares(lambda x: M.read(x), np.array(q0, float), method="lm", diff_step=diff_step, max_nfev=10 ** 6)
    raise OutOfBudget


def _loss(F, loss):
    n2 = float(F @ F)
    return float(np.log(np.sqrt(n2) + 1e-300)) if loss == "log" else n2


def nelder_mead(M, q0, simplex=0.3, loss="sq"):
    q0 = np.array(q0, float)
    sim = np.vstack([q0] + [q0 + simplex * E4[j] for j in range(4)])
    minimize(lambda x: _loss(M.read(x), loss), q0, method="Nelder-Mead", options=dict(initial_simplex=sim, xatol=1e-14, fatol=1e-14, maxiter=10 ** 6, maxfev=10 ** 6))
    raise OutOfBudget


def powell(M, q0, loss="sq"):
    minimize(lambda x: _loss(M.read(x), loss), np.array(q0, float), method="Powell", options=dict(maxiter=10 ** 6, maxfev=10 ** 6))
    raise OutOfBudget


def cma_es(M, q0, sigma=0.5, pop=8, loss="sq"):
    import cma
    es = cma.CMAEvolutionStrategy(list(q0), sigma, {"popsize": pop, "verbose": -9, "seed": 1})
    while True:
        xs = es.ask()
        es.tell(xs, [_loss(M.read(np.array(x)), loss) for x in xs])
