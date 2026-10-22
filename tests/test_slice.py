"""Garde-fous de non-commutativité et de la tranche : python tests/test_slice.py"""
import os, sys
import numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "src"))
from quat import qmul, qinv
from funcs import slice_extend


def test_slice_matches_complex():
    g = lambda z: z * z + 1
    f = slice_extend(g)
    for z in (0.3 + 0.8j, -1.2 + 0.1j, 2.0 + 0.5j):
        F = f(np.array([z.real, z.imag, 0, 0]))
        assert np.allclose(F[0] + 1j * F[1], g(z)) and abs(F[2]) + abs(F[3]) < 1e-12


def test_slice_rotation_covariance():
    """f(q) = A + (v/|v|) B : l'image reste dans la même tranche que l'argument (le quotient reste dans ℂ_u)."""
    f = slice_extend(lambda z: z ** 3 - 1)
    q = np.array([0.4, 0.2, -0.5, 0.3]); h = np.array([0.0, 1e-6, 0, 0]) * 0
    F = f(q)
    assert np.allclose(F[1:] / np.linalg.norm(F[1:]), q[1:] / np.linalg.norm(q[1:]))


def test_hamilton_order_matters():
    rng = np.random.default_rng(0)
    a, b = rng.standard_normal(4), rng.standard_normal(4)
    assert not np.allclose(qmul(a, b), qmul(b, a))
    # quotient à gauche : F^{-1}(F ρ) = ρ, mais F (F^{-1} ρ) != ρ F^{-1} F en général
    assert np.allclose(qmul(qinv(a), qmul(a, b)), b)
    assert not np.allclose(qmul(qmul(a, b), qinv(a)), b)


if __name__ == "__main__":
    for n, fn in list(globals().items()):
        if n.startswith("test_"):
            fn(); print("ok", n)
