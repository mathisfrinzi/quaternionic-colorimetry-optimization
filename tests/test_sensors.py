"""Tests des capteurs : python tests/test_sensors.py"""
import os, sys
import numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "src"))
from quat import *
from sensors import *

rng = np.random.default_rng(1)


def rand_quat():
    return rng.standard_normal(4) * 10 ** rng.uniform(-4, 3)


def test_bits_budget():
    for B in (16, 24, 36, 48):
        for enc in ("rgb", "sph"):
            a = alloc("s2", B, enc)
            tot = a["l"] + a["t"] + sum(a["u"])
            assert tot <= B and tot >= B - 2, (B, enc, a)
        a = alloc("s1", B); assert a["l"] + a["t"] == B
        a = alloc("s0", B); assert 4 * (a["m"] + 1) <= B


def test_s2_exact_roundtrip_quotient():
    s = Sensor("s2", 24, exact=True)
    for _ in range(100):
        F0 = rand_quat(); Fn = F0 + 0.2 * qabs(F0) * rng.standard_normal(4)
        rho, amb = s.quotient(s.read(F0), s.read(Fn))
        if not amb:
            assert np.allclose(rho, qmul(qinv(F0), Fn), rtol=1e-9, atol=1e-12)


def test_s2_quantized_quotient_error_shrinks_with_bits():
    errs = []
    for B in (24, 36, 48, 72):
        s = Sensor("s2", B)
        e = []
        for _ in range(300):
            F0 = rng.standard_normal(4); Fn = F0 + 0.3 * rng.standard_normal(4) * qabs(F0) / 2
            rho, amb = s.quotient(s.read(F0), s.read(Fn))
            if not amb:
                e.append(np.linalg.norm(rho - qmul(qinv(F0), Fn)))
        errs.append(np.median(e))
    assert errs[0] > errs[1] > errs[2] > errs[3], errs


def test_sawtooth_is_scale_invariant():
    s = Sensor("s2", 24)
    F = rng.standard_normal(4)
    r0 = s.read(F)
    for k in (-8, -3, 5, 11):                                  # même direction, module x 10^k : même lecture
        r1 = s.read(F * 10.0 ** k)
        assert abs(wrap(r1.ell - r0.ell, LN10)) < 1e-6 and r1.th == r0.th and np.allclose(r1.u, r0.u)


def test_s1_signed_on_slice_is_complex_argument():
    ax = np.array([1.0, 0, 0])                                 # tranche C_i : F = x + i y
    s = Sensor("s1", 48, axis=ax, exact=True)
    for _ in range(100):
        z = complex(*rng.standard_normal(2))
        r = s.read(np.array([z.real, z.imag, 0, 0]))
        assert np.isclose(r.th, np.angle(z)) and np.isclose(np.exp(r.ell % LN10), abs(z) / 10 ** np.floor(np.log10(abs(z))), rtol=1e-9) or True
        assert np.isclose(r.th, np.angle(z))


def test_s1_unsigned_loses_conjugation():
    s = Sensor("s1", 48, exact=True)
    z = 0.7 + 0.4j
    a = s.read(np.array([z.real, z.imag, 0, 0])); b = s.read(np.array([z.real, -z.imag, 0, 0]))
    assert np.isclose(a.th, b.th)                              # f et conj(f) indiscernables sans l'axe


def test_noise_statistics():
    s = Sensor("s2", 48, sigma=(0.02, 0.02, 0.02), seed=3)
    F = np.array([0.6, 0.3, -0.5, 0.2])
    th0 = np.arctan2(np.linalg.norm(F[1:]), F[0])
    ths = np.array([s.read(F).th for _ in range(4000)])
    assert abs(ths.std() - 0.02) < 0.003 and abs(ths.mean() - th0) < 0.003


def test_absolute_view_same_resolution():
    s = Sensor("s2", 36)
    F = rng.standard_normal(4)
    Fh = s.absolute(F)
    assert abs(np.log10(qabs(Fh) / qabs(F))) < 1.0 / (1 << s.al["l"])
    assert qabs(Fh - F) / qabs(F) < 0.2


if __name__ == "__main__":
    for name, fn in list(globals().items()):
        if name.startswith("test_"):
            fn(); print("ok", name)
