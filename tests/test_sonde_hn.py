"""La sonde générale (src/sonde_hn.py) reproduit exactement la sonde de la partie 3.4 (probe34) et vérifie ses modèles.
python tests/test_sonde_hn.py"""
import os, sys
import numpy as np
R = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
sys.path.insert(0, os.path.join(R, "src"))
for d in ("3.1_jacobienne_H_lineaire", "3.2_grande_dimension", "3.4_ambiguite"):
    sys.path.insert(0, os.path.join(R, "parties", "3_dimension_4n", d))
from sensors import Sensor
from probes import Hit
from quat import qmul, Lmat, Rmat, rand_unit
import sonde_hn as S
from bench_hn import MeterN
import bench_hn_grand as G
from sonde34 import probe34


def _run34(n, sigma, seed, i, opts):
    f, qs = G.problem(n, 0.0, seed); q0 = G.starts(n, qs, i + 1, seed)[i]
    sens = [Sensor("s2", 24, "sph", sigma=(sigma,) * 3, seed=1000 * i + k) for k in range(n)]
    M = MeterN(f, n, "color", sens, 600)
    try:
        probe34(M, q0, n, c=0.3, step=0.65, cap=1.5, **opts)
    except Hit:
        return True, M.cnt
    except S.ERRS:
        return False, 600
    return False, 600


def test_reproduit_probe34():
    for n, sigma, opts in [(4, 0.0, {}), (4, 0.02, {}), (8, 0.0, dict(lum="prior0")),
                           (8, 0.02, dict(verif_k=2, tau=0.25, double=True))]:
        for i in range(3):
            f, qs = G.problem(n, 0.0, 100); q0 = G.starts(n, qs, i + 1, 100)[i]
            a = _run34(n, sigma, 100, i, opts)
            ok, c, *_ = S.lancer(f, q0, n, "sonde", budget=600, sigma=sigma, i=i, c=0.3, step=0.65, cap=1.5, **opts)
            assert (ok, c) == a, (n, sigma, opts, i, (ok, c), a)


def test_quotient_droite():
    rng = np.random.default_rng(1)
    F0, F1 = rng.standard_normal((2, 4))
    s = Sensor("s2", 24, "sph", exact=True)
    r, _ = S.quotient(s.read(F0), s.read(F1), np.log(np.linalg.norm(F1) / np.linalg.norm(F0)), "R")
    assert np.allclose(r, qmul(F1, S.qinv(F0)), atol=1e-9)


def test_structures_conservees():
    """Normaliser à gauche garde les blocs L ; à droite, les blocs R ; à gauche, les blocs conformes."""
    rng = np.random.default_rng(2)
    a, c, F = rng.standard_normal((3, 4))
    LFi, RFi = Lmat(S.qinv(F)), Rmat(S.qinv(F))
    assert np.allclose(LFi @ Lmat(a), Lmat(qmul(S.qinv(F), a)))
    assert np.allclose(RFi @ Rmat(c), Rmat(qmul(c, S.qinv(F))))
    from quat import conf_defect
    assert conf_defect(LFi @ Lmat(a) @ Rmat(c)) < 1e-9


def test_deux_echelles_decade():
    rng = np.random.default_rng(3)
    s = S.CapteurDeuxEchelles(24)
    for _ in range(200):
        F = rng.standard_normal(4) * 10 ** rng.uniform(-5, 2)
        assert abs(s.read(F).lg - np.log10(np.linalg.norm(F))) < 0.05


def test_couche_par_sortie():
    """Partie 6.4 : (a) jacobienne de la couche exacte ; (b) le modèle Dbloc lit, avec un capteur exact, la vraie jacobienne
    normalisée à droite, en p sondes ; (c) avec m = 1, la sonde Dbloc fait exactement les mêmes lectures que la sonde D."""
    sys.path.insert(0, os.path.join(R, "parties", "6_reseaux"))
    from bench_couche import Couche
    from quat import qinv
    L = Couche(3, 2, "lin", 0); q = L.qs + 0.05                      # σ = identité : blocs exactement R(·)
    Jn = np.column_stack([(L.f(q + 1e-6 * e) - L.f(q - 1e-6 * e)) / 2e-6 for e in np.eye(24)])
    assert np.abs(Jn - L.J(q)).max() < 1e-7
    n = 6; sens = [Sensor("s2", 24, "sph", exact=True) for _ in range(n)]
    M = S.MeterLog(L.f, n, "color", sens, 10 ** 6)
    mod = S.MODELES["Dbloc"](n, m=3); r0 = M.read(q); eta = 1e-6
    dirs = mod.sondes(eta); assert len(dirs) == 2
    rhos = [np.array([S.quotient(a, b, 0.0, "R")[0] for a, b in zip(r0, M.read(q + e))]) for e in dirs]
    mod.lire(dirs, rhos, eta)
    F = L.f(q).reshape(n, 4); J = L.J(q)
    Wt = np.vstack([Rmat(qinv(F[r])) @ J[4 * r:4 * r + 4] for r in range(n)])
    assert np.linalg.norm(mod.mat() - Wt) / np.linalg.norm(Wt) < 1e-4
    L1 = Couche(1, 4, "radial-faible", 1)
    for i in range(3):
        q0 = L1.qs + 0.3 * np.random.default_rng(i).standard_normal(16)
        a = S.lancer(L1.f, q0, 4, "sonde", budget=600, i=i, modele="D")[:2]
        b = S.lancer(L1.f, q0, 4, "sonde", budget=600, i=i, modele="Dbloc", m=1)[:2]
        assert a == b, (a, b)


def test_broyden_q():
    """Partie 5.5 : le Broyden quaternionique converge sur un problème ℍ-linéaire avec n lectures de départ (au lieu de 4n),
    à gauche, et sur un neurone linéaire à droite ; le suivi (sans recherche linéaire) converge aussi."""
    sys.path.insert(0, os.path.join(R, "parties", "3_dimension_4n", "3.2_grande_dimension"))
    sys.path.insert(0, os.path.join(R, "parties", "3_dimension_4n", "3.1_jacobienne_H_lineaire"))
    sys.path.insert(0, os.path.join(R, "parties", "6_reseaux"))
    import bench_hn_grand as G
    from bench_entrainement import Neurone
    f, qs = G.problem(4, 0.0, 100); q0 = G.starts(4, qs, 1, 100)[0]
    ok, c, *_ = S.lancer(f, q0, 4, "broyden-q", budget=500, side="L"); assert ok and c < 40, (ok, c)
    ok2, c2, *_ = S.lancer(f, q0, 4, "broyden-q", budget=500, side="L", suivi=True); assert ok2
    N = Neurone(4, "lin", 0)
    ok3, c3, *_ = S.lancer(N.f, N.qs + 0.3, 4, "broyden-q", budget=500, side="R"); assert ok3, c3


if __name__ == "__main__":
    for name, fn in list(globals().items()):
        if name.startswith("test_"):
            fn(); print("ok", name)
