"""Partie 5 : armes égales. TOUTES les méthodes lisent la MÊME mesure dégradée (S2 sphérique, B bits, bruit σ).
Les méthodes couleur voient la lecture périodique (luminosité en dent de scie) ; les méthodes classiques voient la vue
`absolute` (même résolution, même bruit, mais luminosité non périodique : avantage donné au classique).
Familles : lin, quad, cub (génériques) + Möbius exact (ε=0). 3 instances × 30 départs, budget 600, Usage : python bench_armes_egales.py [--smoke]"""
import os, sys, json
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "src"))
import bench, classic
from probes import probe_fd, probe_secant

pick = lambda rng, xs: xs[int(rng.integers(len(xs)))]
BITS = [24, 36]; SIG = [0.0, 0.005, 0.02]
fd_s = lambda rng: dict(c=pick(rng, [0.15, 0.3, 0.45]), step=pick(rng, [0.5, 0.65, 0.8]), cap=pick(rng, [0.5, 1.0, 2.0]))
fd_h = lambda rng: dict(hfd=pick(rng, [1e-3, 1e-2, 1e-1, 3e-1]))


def register_all():
    R = bench.register
    for B in BITS:
        for s in SIG:
            sen = dict(kind="s2", bits=B, enc="sph", sigma=(s, s, s))
            t = f"@{B}/{s}"
            R("sonde-rank1" + t, probe_secant, fd_s, ncfg=6, sensor=sen, fixed=dict(model="rank1"))
            R("sonde-fd" + t, probe_fd, fd_s, ncfg=6, sensor=sen, fixed=dict(model="full"))
            R("broyden" + t, classic.broyden, fd_h, ncfg=4, sensor=sen, absolute=True)
            R("newton-FD" + t, classic.newton_fd, fd_h, ncfg=4, sensor=sen, absolute=True)
            R("bfgs-FD" + t, classic.bfgs_fd, lambda rng: dict(hfd=pick(rng, [1e-3, 1e-2, 1e-1])), ncfg=3, sensor=sen, absolute=True)
            R("nelder-mead" + t, classic.nelder_mead, lambda rng: dict(simplex=pick(rng, [0.3, 1.0]), loss=pick(rng, ["log", "sq"])), ncfg=3, sensor=sen, absolute=True)


if __name__ == "__main__":
    smoke = "--smoke" in sys.argv
    register_all(); methods = list(bench.REGISTRY)
    fam = [("gen", k, s) for k in ("lin", "quad", "cub") for s in (0, 1, 2)] + [("mobius", s, 0.0) for s in (0, 1, 2)]
    tune = [("gen", k, 10) for k in ("lin", "quad", "cub")]
    res = bench.run_bench(methods, tune, fam, ntune=2 if smoke else 5, ntest=3 if smoke else 30, budget=600, ncfg=1 if smoke else None)
    res.update(bits=BITS, sig=SIG)
    json.dump(res, open(os.path.join(HERE, "bench_smoke.json" if smoke else "bench_armes_egales_results.json"), "w"))
    for m in methods:
        s = bench.summarize(res, m); print(f"{m:26s} {100*s['rate']:5.1f}%  coût médian {s['cost']:.0f}")
