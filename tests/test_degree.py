"""Tests du degré topologique : python tests/test_degree.py"""
import os, sys
import numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "src"))
from quat import qmul, qconj, E4
from degree import *
from sensors import Sensor

exact = Sensor("s2", 24, exact=True)
ID = lambda q: np.asarray(q, float)
CONJ = lambda q: qconj(q)
SQ = lambda q: qmul(q, q)
CONJSQ = lambda q: qmul(qconj(q), qconj(q))
CUBE = lambda q: qmul(q, qmul(q, q))


def test_triangulation_counts_and_orientation():
    for d, (nv, nt) in {0: (8, 16), 1: (32, 128), 2: (192, 1024)}.items():
        V, T, o = triangulation(d)
        assert len(V) == nv and len(T) == nt, (d, len(V), len(T))
        assert np.allclose(np.linalg.norm(V, axis=1), 1)
    V, T, o = triangulation(2)
    assert np.all(o != 0)
    # le degré de l'identité vaut +1 si toutes les orientations sont cohérentes : somme des signes sur un point = 1
    G = V.copy()
    assert degree(G, 2)[0] == 1


def test_known_degrees_exact():
    for name, f, d in [("q", ID, 1), ("conj", CONJ, -1), ("q^2", SQ, 2), ("conj^2", CONJSQ, -2), ("q^3", CUBE, 3)]:
        got, agree, nr = certificate(f, exact, np.zeros(4), 1.0, depth=2)
        assert got == d and agree == 1.0, (name, got, d, agree)


def test_no_zero_gives_degree_zero():
    f = lambda q: np.asarray(q) + np.array([3.0, 0, 0, 0])
    assert certificate(f, exact, np.zeros(4), 1.0, depth=2)[0] == 0


def test_off_center_and_radius():
    f = lambda q: np.asarray(q) - np.array([0.3, 0.1, -0.2, 0.0])          # zéro en xs, |xs| < 1
    assert certificate(f, exact, np.zeros(4), 1.0, depth=2)[0] == 1
    assert certificate(f, exact, np.zeros(4), 0.2, depth=2)[0] == 0         # sphère trop petite : le zéro est dehors


def test_quantized_sensor_keeps_degree():
    s = Sensor("s2", 36, enc="sph")
    got, agree, nr = certificate(CONJSQ, s, np.zeros(4), 1.0, depth=2)
    assert got == -2 and agree > 0.9


if __name__ == "__main__":
    for name, fn in list(globals().items()):
        if name.startswith("test_"):
            fn(); print("ok", name)
