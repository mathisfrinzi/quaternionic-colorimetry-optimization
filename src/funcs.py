"""Fonctions de test.

Un `Problem` regroupe f : R^4 -> R^4 (quaternion -> quaternion), le rayon de tirage des départs, le TYPE DE CIBLE (point isolé,
sphère, courbe, surface, aucun zéro : on le signale dans chaque tableau, cf. fausse piste 5) et le seuil de succès.

Trois familles :
  * slice-régulières (partie 1) : g holomorphe à coefficients réels, étendue par f(x + v) = A(x,|v|) + (v/|v|) B(x,|v|) où
    g(x + i y) = A + i B. Sur la tranche C_i elle redonne exactement g ; les zéros non réels de g deviennent des SPHÈRES de zéros.
  * génériques (partie 2) : f lisse sans structure, zéro isolé / courbe / surface / aucun.
  * Möbius -> générique (partie 3) : f_ε = M + ε G(q - q*), M une transformation de Möbius, G une perturbation lisse aléatoire
    qui s'annule en q* (le zéro reste le même). ε règle la distance à la conformité, mesurée par `conf_defect` de la jacobienne.
"""
import numpy as np
from quat import qmul, qinv, qabs, mobius, jac_num, conf_defect, rand_unit, E4, Lmat, Rmat

np.seterr(all="ignore")


class Problem:
    def __init__(self, name, f, R, target, thr=1e-3, center=None, zero=None, meta=None, g_complex=None):
        self.name, self.f, self.R, self.target, self.thr = name, f, R, target, thr
        self.center = np.zeros(4) if center is None else np.asarray(center, float)
        self.zero, self.meta = zero, (meta or {})
        self.g_complex = g_complex

    def __call__(self, q):
        return self.f(q)

    def starts(self, n, seed):
        """n départs uniformes dans la boule de rayon R autour de `center`."""
        rng = np.random.default_rng(seed)
        out = []
        for _ in range(n):
            d = rand_unit(rng)
            out.append(self.center + self.R * rng.uniform() ** 0.25 * d)
        return np.array(out)

    def slice_starts(self, n, seed):
        """n départs dans la tranche C_i (x + i y), carré [-R, R]^2 comme dans le repo complexe."""
        rng = np.random.default_rng(seed)
        xy = rng.uniform(-self.R, self.R, (n, 2))
        return np.column_stack([xy, np.zeros((n, 2))])

    def jac(self, q):
        return jac_num(self.f, q)

    def defect(self, q):
        return conf_defect(self.jac(q))


# ------------------------------------------------------------------------------------------------ slice-régulières
def slice_extend(g):
    """Extension slice-régulière d'une fonction holomorphe g à coefficients réels."""
    def f(q):
        q = np.asarray(q, float)
        v = q[1:]
        y = float(np.sqrt(v @ v))
        w = complex(g(complex(q[0], y)))
        if y < 1e-300:
            return np.array([w.real, 0.0, 0.0, 0.0])
        return np.concatenate([[w.real], v / y * w.imag])
    return f


def _real_poly(roots_pairs, reals):
    """Polynôme réel à partir de paires de racines conjuguées et de racines réelles."""
    c = np.array([1.0])
    for r in reals:
        c = np.convolve(c, [1.0, -r])
    for r in roots_pairs:
        c = np.convolve(c, [1.0, -2 * r.real, abs(r) ** 2])
    return lambda z: np.polyval(c, z)


SLICE_SPECS = {
    # nom : (g holomorphe réelle, rayon, type de cible)
    "q2+1": (lambda z: z ** 2 + 1, 2.0, "sphère"),
    "q2-1": (lambda z: z ** 2 - 1, 2.0, "2 points réels"),
    "z3-1": (lambda z: z ** 3 - 1, 1.8, "point + sphère"),
    "quad": (_real_poly([-0.5 + 1.1j], []), 2.5, "sphère"),
    "poly5": (_real_poly([0.4 + 1.1j, -0.8 + 0.7j], [1.3]), 2.5, "point + 2 sphères"),
    "tanh": (np.tanh, 4.0, "point + sphères"),
    "sinz-z": (lambda z: np.sin(z) - z, 6.0, "point + sphères"),
    "exp-z-1": (lambda z: np.exp(z) - z - 1, 8.0, "point + sphères"),
    "rationnelle": (lambda z: (z ** 2 + 1) / (z ** 2 - 1), 2.5, "sphère"),
    "sigmoide-0.5": (lambda z: 1 / (1 + np.exp(-z)) - 0.5, 5.0, "point + sphères"),
}


def slice_problem(name, tol=1e-3):
    g, R, target = SLICE_SPECS[name]
    return Problem("slice:" + name, slice_extend(g), R, target, thr=tol, g_complex=g)


# ------------------------------------------------------------------------------------------------ génériques
def _rand_cond_matrix(rng, cond, n=4):
    Q1, _ = np.linalg.qr(rng.standard_normal((n, n)))
    Q2, _ = np.linalg.qr(rng.standard_normal((n, n)))
    return Q1 @ np.diag(np.geomspace(1.0, cond, n)) @ Q2


GENERIC_KINDS = ["lin", "quad", "cub", "ill", "curve", "surf", "nozero", "nozero_ill"]


def generic_problem(kind, seed, R=2.0):
    """Fonction lisse générique de R^4 dans R^4. Le zéro (ou le minimum) est en xs."""
    rng = np.random.default_rng(1000 + 7919 * seed)
    xs = 0.5 * rng.standard_normal(4)
    B = rand_unit(rng)
    if kind == "lin":
        A = _rand_cond_matrix(rng, 10.0)
        f = lambda q: A @ (np.asarray(q) - xs)
        target, thr = "point", 1e-3
    elif kind == "quad":
        A = _rand_cond_matrix(rng, 5.0)
        f = lambda q: A @ (np.asarray(q) - xs) + 0.5 * qmul(qmul(np.asarray(q) - xs, B), np.asarray(q) - xs)
        target, thr = "point (autres zéros possibles)", 1e-3
    elif kind == "cub":
        A = _rand_cond_matrix(rng, 5.0)
        f = lambda q: A @ (np.asarray(q) - xs) + 0.4 * (np.asarray(q) - xs) ** 3
        target, thr = "point", 1e-3
    elif kind == "ill":
        A = _rand_cond_matrix(rng, 100.0)
        f = lambda q: A @ (np.asarray(q) - xs) + 0.3 * np.tanh(np.asarray(q) - xs)
        target, thr = "point (mal conditionné)", 1e-3
    elif kind == "curve":                          # 3 équations en 4 inconnues : zéros = une courbe
        A = _rand_cond_matrix(rng, 4.0)
        def f(q):
            d = np.asarray(q) - xs
            y = A @ d + 0.4 * np.roll(d, 1) * d
            return np.array([y[0], y[1], y[2], 0.0])
        target, thr = "courbe (1D)", 1e-3
    elif kind == "surf":                           # 2 équations : zéros = une surface
        A = _rand_cond_matrix(rng, 4.0)
        def f(q):
            d = np.asarray(q) - xs
            y = A @ d + 0.4 * np.roll(d, 1) * d
            return np.array([y[0], y[1], 0.0, 0.0])
        target, thr = "surface (2D)", 1e-3
    elif kind in ("nozero", "nozero_ill"):         # |f| >= 1, minimum 1 en xs (nozero_ill : vallée étroite, conditionnement 100)
        A = _rand_cond_matrix(rng, 3.0 if kind == "nozero" else 100.0)
        def f(q):
            d = np.asarray(q) - xs
            y = A @ d
            return np.array([1.0 + d @ d, y[1], y[2], y[3]])
        target, thr = "aucun zéro (min |f| = 1)", 1.05
    else:
        raise ValueError(kind)
    return Problem(f"gen:{kind}", f, R, target, thr=thr, center=xs, zero=xs, meta=dict(seed=seed))


# ------------------------------------------------------------------------------------------------ Möbius -> générique
def mobius_problem(seed, eps, R=1.8, pole_dist=3.0, qstar_norm=1.0, tol=1e-3):
    """f_ε(q) = M(q) + ε G(q - q*), M(q) = (a q + b)(c q + d)^{-1}, zéro en q* = -a^{-1} b.
    La perturbation G(h) = A_g h + 0.15 · h B h est normalisée pour que ‖A_g‖_F = ‖J_M(q*)‖_F : ε = 1 met les deux sur
    le même pied au zéro. Le pôle de M est à distance ≈ pole_dist de l'origine (hors de la boule de départ)."""
    rng = np.random.default_rng(5000 + 104729 * seed)
    qs = qstar_norm * rand_unit(rng)
    a = rand_unit(rng)
    b = -qmul(a, qs)
    c = 0.35 * rand_unit(rng)
    d = pole_dist * 0.35 * rand_unit(rng)
    M = mobius(a, b, c, d)
    Jm = jac_num(M, qs)
    Ag = rng.standard_normal((4, 4))
    Ag *= np.linalg.norm(Jm) / np.linalg.norm(Ag)
    Bq = rand_unit(rng)
    def G(h):
        h = np.asarray(h, float)
        return Ag @ h + 0.15 * qmul(qmul(h, Bq), h)
    f = lambda q: M(q) + eps * G(np.asarray(q) - qs)
    J0 = jac_num(f, qs)
    P = Problem(f"mobius:eps={eps}", f, R, "point (isolé)", thr=tol, center=qs, zero=qs,
                meta=dict(seed=seed, eps=eps, defect_zero=conf_defect(J0)))
    # défaut moyen sur 16 points de la boule de départ
    pts = P.starts(16, 7)
    P.meta["defect_mean"] = float(np.mean([conf_defect(jac_num(f, p)) for p in pts]))
    return P


def build(spec):
    """Fabrique un Problem à partir d'un tuple picklable (pour les processus)."""
    kind = spec[0]
    if kind == "slice":
        return slice_problem(spec[1])
    if kind == "gen":
        return generic_problem(spec[1], spec[2])
    if kind == "mobius":
        return mobius_problem(spec[1], spec[2])
    raise ValueError(spec)
