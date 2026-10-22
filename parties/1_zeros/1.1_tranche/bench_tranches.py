"""Partie 1.1 (M1) : sanity check sur la tranche C_i.

Fonctions slice-régulières à coefficients réels, départs dans la tranche. Sur C_i, f est une fonction holomorphe ordinaire ;
le capteur S1 SIGNÉ (θ autour de l'axe i) se comporte comme le capteur complexe. On compare :
  - le code complexe du repo, tel quel (vendu dans vendor/), alimenté par la restriction à C_i de la fonction quaternionique ;
  - la même sonde sur le capteur quaternionique S1 (signé), à 16 et 24 bits ; S1 non signé (le signe de θ est perdu) ;
  - la sonde 4D générique sur S2 (elle ne sait pas qu'on est sur une tranche) ;
  - les méthodes classiques, qui lisent la valeur exacte de f (donc favorisées).
Critère de succès (décision documentée) : atteindre UN point de zéro (|f| < 1e-3), pas identifier la sphère entière.
Usage : python bench_tranches.py [--smoke]
"""
import os, sys, json
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", "..", ".."))
sys.path[:0] = [os.path.join(ROOT, "src"), os.path.join(ROOT, "vendor", "complex_colorimetry")]
import bench
from probes import probe_slice, probe_fd, Hit, OutOfBudget
import classic
from funcs import SLICE_SPECS

pick = lambda rng, xs: xs[int(rng.integers(len(xs)))]
AXIS_I = (1.0, 0.0, 0.0)          # axe de la tranche C_i dans l'espace des parties vectorielles (i = e1)


def complex_original(M, q0, bits=8, c=0.3, step_frac=0.7, cap=1.0, blind=True):
    """Le code complexe du repo, sans modification, sur la restriction de f à C_i."""
    from color_probe import descend_probe
    f_c = lambda z: complex(*M.f(np.array([z.real, z.imag, 0.0, 0.0]))[:2])
    z, it, ev = descend_probe(f_c, complex(q0[0], q0[1]), bits=bits, blind_rule=blind, c=c, step_frac=step_frac, cap=cap)
    cost = ev + 1
    if abs(f_c(z)) < M.thr and cost <= M.budget:
        M.n = cost
        raise Hit
    raise OutOfBudget


def newton_complex(M, q0):
    """Newton complexe : f(z) et f'(z) (facturée 1 évaluation) à chaque pas, comme dans le repo complexe."""
    z = complex(q0[0], q0[1])
    fc = lambda zz: complex(*M.f(np.array([zz.real, zz.imag, 0.0, 0.0]))[:2])
    while True:
        F = M.read(np.array([z.real, z.imag, 0.0, 0.0])); fz = complex(F[0], F[1])
        M.charge(1)
        h = 1e-6
        dz = (fc(z + h) - fc(z - h)) / (2 * h)
        if dz == 0 or not np.isfinite(dz):
            raise OutOfBudget
        z = z - fz / dz
        if not np.isfinite(z) or abs(z) > 1e6:
            raise OutOfBudget


def probe_sample(rng):
    return dict(c=pick(rng, [0.1, 0.2, 0.3, 0.45]), step_frac=pick(rng, [0.3, 0.5, 0.7, 0.9]), cap=pick(rng, [0.5, 1.0, 2.0]),
                blind_rule=bool(pick(rng, [True, False])))


def orig_sample(rng):
    p = probe_sample(rng); p["blind"] = p.pop("blind_rule"); return p


def fd_sample(rng):
    return dict(c=pick(rng, [0.15, 0.3, 0.45]), step=pick(rng, [0.5, 0.65, 0.8]), cap=pick(rng, [0.5, 1.0, 2.0]))


def register_all():
    R = bench.register
    for b in (8, 12):
        R(f"complexe-original-{b}b", complex_original, orig_sample, ncfg=8, starts="slice", fixed=dict(bits=b))
    for B in (16, 24):
        R(f"S1-signé-{B}b", probe_slice, probe_sample, ncfg=8, starts="slice", sensor=dict(kind="s1", bits=B, axis=AXIS_I))
    R("S1-non-signé-24b", probe_slice, probe_sample, ncfg=8, starts="slice", sensor=dict(kind="s1", bits=24))
    R("S2-4D-24b", probe_fd, fd_sample, ncfg=6, starts="slice", sensor=dict(kind="s2", bits=24), fixed=dict(model="full"))
    R("newton-complexe", newton_complex, lambda rng: {}, ncfg=1, starts="slice")
    R("newton-4D-FD", classic.newton_fd, lambda rng: dict(hfd=pick(rng, [1e-6, 1e-4])), ncfg=2, starts="slice")
    R("newton-4D-Jexacte", classic.newton_exactJ, lambda rng: {}, ncfg=1, starts="slice")
    R("broyden-4D", classic.broyden, lambda rng: dict(hfd=pick(rng, [1e-6, 1e-4])), ncfg=2, starts="slice")
    R("bfgs-4D", classic.bfgs_fd, lambda rng: {}, ncfg=1, starts="slice")
    R("nelder-mead-4D", classic.nelder_mead, lambda rng: dict(simplex=pick(rng, [0.3, 1.0]), loss=pick(rng, ["log", "sq"])), ncfg=3, starts="slice")
    R("cma-es-4D", classic.cma_es, lambda rng: dict(sigma=pick(rng, [0.3, 1.0]), pop=pick(rng, [4, 6]), loss=pick(rng, ["log", "sq"])), ncfg=3, starts="slice")


if __name__ == "__main__":
    smoke = "--smoke" in sys.argv
    register_all()
    methods = list(bench.REGISTRY)
    specs = [("slice", n) for n in SLICE_SPECS]
    res = bench.run_bench(methods, specs, specs, ntune=2 if smoke else 8, ntest=4 if smoke else 100, budget=1000,
                          ncfg=1 if smoke else None, out=None)
    res["specs"] = specs
    json.dump(res, open(os.path.join(HERE, "bench_tranches_smoke.json" if smoke else "bench_tranches_results.json"), "w"))
    for m in methods:
        s = bench.summarize(res, m)
        print(f"{m:22s} {100*s['rate']:5.1f}% ({100*s['lo']:.0f}-{100*s['hi']:.0f})  coût médian {s['cost']:.0f}")
