"""Tests de l'algèbre : python tests/test_quat.py  (ou pytest depuis quaternionic/)."""
import os, sys
import numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "src"))
from quat import *

rng = np.random.default_rng(0)


def test_hamilton():
    i, j, k = E4[1], E4[2], E4[3]
    assert np.allclose(qmul(i, j), k) and np.allclose(qmul(j, i), -k)       # non commutatif
    assert np.allclose(qmul(i, i), -E4[0])
    a, b, c = rng.standard_normal((3, 4))
    assert np.allclose(qmul(qmul(a, b), c), qmul(a, qmul(b, c)))            # associatif
    assert np.allclose(qmul(a, qinv(a)), E4[0])
    assert np.isclose(qabs(qmul(a, b)), qabs(a) * qabs(b))                  # le module est multiplicatif


def test_L_R():
    a, b = rng.standard_normal((2, 4))
    assert np.allclose(Lmat(a) @ b, qmul(a, b))
    assert np.allclose(Rmat(b) @ a, qmul(a, b))


def test_basis_orthogonal():
    B = BASIS.reshape(16, 16)
    assert np.allclose(B @ B.T, 4 * np.eye(16))


def test_rank_one_is_conformal():
    a, b = rng.standard_normal((2, 4))
    M = Lmat(a) @ Rmat(b)
    assert np.linalg.matrix_rank(to_Z(M), tol=1e-10) == 1
    assert conf_defect(M) < 1e-12
    assert np.allclose(conf_project(M), M)
    G = M.T @ M                                                             # conforme : M^T M = λ I
    assert np.allclose(G, G[0, 0] * np.eye(4))
    assert np.linalg.det(M) > 0


def test_generic_matrix_has_defect_and_projection_is_nearest():
    M = rng.standard_normal((4, 4))
    d = conf_defect(M)
    P = conf_project(M)
    assert 0.2 < d < 1
    assert np.isclose(np.linalg.norm(M - P) / np.linalg.norm(M), d)         # défaut = distance relative
    for _ in range(200):                                                    # aucune conforme tirée au hasard n'est plus proche
        a, b = rng.standard_normal((2, 4))
        assert np.linalg.norm(M - Lmat(a) @ Rmat(b)) >= np.linalg.norm(M - P) - 1e-9


def test_orientation_reversing_is_not_direct_conformal():
    M = np.diag([1.0, -1, -1, -1])                                          # conjugaison : conforme, mais inverse
    assert conf_defect(M) > 0.5


def test_mobius_jacobian_is_conformal_with_7_dof():
    a, b, c, d = rng.standard_normal((4, 4))
    f = mobius(a, b, c, d)
    q = rng.standard_normal(4)
    J = jac_num(f, q)
    assert conf_defect(J) < 1e-6
    al = a - qmul(f(q), c)                                                  # J h = (a - f c) h (c q + d)^{-1}
    be = qinv(qmul(c, q) + d)
    assert np.allclose(J, Lmat(al) @ Rmat(be), atol=1e-5)


def test_cross3():
    a, b, c = rng.standard_normal((3, 4))
    n = cross3(a, b, c)
    assert abs(n @ a) < 1e-10 and abs(n @ b) < 1e-10 and abs(n @ c) < 1e-10
    assert np.linalg.det(np.stack([a, b, c, n], 1)) > 0


if __name__ == "__main__":
    for name, fn in list(globals().items()):
        if name.startswith("test_"):
            fn(); print("ok", name)
