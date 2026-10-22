"""Partie 1.4 : trois ajouts à la sonde sécante (`sonde-rank1`), testés un par un puis ensemble, dans le cadre « armes égales »
de la partie 6 (même mesure dégradée pour toutes les méthodes, même départs, même graine de capteur → comparaisons appariées).

  dith : dithering soustractif connu (capteur) — l'erreur de quantification devient sans biais et indépendante d'une lecture
         à l'autre ; aussi offert à Broyden (`broyden+dith`) pour rester à armes égales.
  conf : démarrage conforme — la première jacobienne normalisée est lue avec 3 sondes et projetée sur L(α)R(β) (4 lectures au
         lieu de 5), puis mises à jour générales de rang 1.
  tr   : région de confiance pilotée par la luminosité au lieu du retour en arrière.

Familles : lin, quad, cub, ill (génériques) + Möbius exact (ε = 0) + Möbius dégradé (ε = 1), 3 instances × 30 départs, budget 600.
Capteur S2 sphérique, B ∈ {24, 36} bits, bruit σ ∈ {0, 0.02} sur chaque canal.
Usage : python bench_ameliorations.py [--smoke]"""
import os, sys, json
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "..", "src"))
import bench, classic
from probes import probe_secant

pick = lambda rng, xs: xs[int(rng.integers(len(xs)))]
BITS = [24, 36]; SIG = [0.0, 0.02]
VARIANTS = {"base": {}, "dith": {}, "conf": dict(init="conf"), "tr": dict(tr=True), "dith+conf+tr": dict(init="conf", tr=True)}
fd_s = lambda rng: dict(c=pick(rng, [0.15, 0.3, 0.45]), step=pick(rng, [0.5, 0.65, 0.8]), cap=pick(rng, [0.5, 1.0, 2.0]))
fd_h = lambda rng: dict(hfd=pick(rng, [1e-3, 1e-2, 1e-1, 3e-1]))


def register_all():
    R = bench.register
    for B in BITS:
        for s in SIG:
            t = f"@{B}/{s}"
            for v, kw in VARIANTS.items():
                sen = dict(kind="s2", bits=B, enc="sph", sigma=(s, s, s), dither="dith" in v)
                R(f"sonde-{v}{t}", probe_secant, fd_s, ncfg=6, sensor=sen, fixed=dict(model="rank1", **kw))
            for dith in (False, True):
                sen = dict(kind="s2", bits=B, enc="sph", sigma=(s, s, s), dither=dith)
                R(("broyden+dith" if dith else "broyden") + t, classic.broyden, fd_h, ncfg=4, sensor=sen, absolute=True)


if __name__ == "__main__":
    smoke = "--smoke" in sys.argv
    register_all(); methods = list(bench.REGISTRY)
    fam = [("gen", k, s) for k in ("lin", "quad", "cub", "ill") for s in (0, 1, 2)] + [("mobius", s, e) for e in (0.0, 1.0) for s in (0, 1, 2)]
    tune = [("gen", k, 10) for k in ("lin", "quad", "cub", "ill")] + [("mobius", 10, 0.0)]
    res = bench.run_bench(methods, tune, fam, ntune=2 if smoke else 5, ntest=3 if smoke else 30, budget=600, ncfg=1 if smoke else None)
    res.update(bits=BITS, sig=SIG, variants=list(VARIANTS))
    json.dump(res, open(os.path.join(HERE, "bench_smoke.json" if smoke else "bench_ameliorations_results.json"), "w"))
    for m in methods:
        s = bench.summarize(res, m); print(f"{m:28s} {100*s['rate']:5.1f}%  coût médian {s['cost']:.0f}", flush=True)
