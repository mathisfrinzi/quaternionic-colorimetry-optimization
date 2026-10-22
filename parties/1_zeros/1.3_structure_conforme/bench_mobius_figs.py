"""Partie 1.3 : figures et tableaux à partir de bench_mobius_results.json. Comparaisons APPARIÉES (mêmes départs, mêmes instances)."""
import os, sys, json
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "..", "src"))
import matplotlib.pyplot as plt
from report import PAL, style, md_table, INK2
from funcs import build
from bench import wilson

R = json.load(open(os.path.join(HERE, "bench_mobius_results.json")))
EPS = R["eps"]; D = R["data"]; BUD = R["budget"]
SEEDS = sorted({json.loads(k)[1] for k in D["broyden"]})
DEF = {e: float(np.mean([build(("mobius", s, e)).meta["defect_zero"] for s in SEEDS])) for e in EPS}
os.makedirs(os.path.join(HERE, "figures"), exist_ok=True)


def arr(m, e):
    ok = np.concatenate([D[m][json.dumps(["mobius", s, e])]["ok"] for s in SEEDS]).astype(bool)
    c = np.concatenate([D[m][json.dumps(["mobius", s, e])]["cost"] for s in SEEDS]).astype(float)
    return ok, c


def rate(m, e): return arr(m, e)[0].mean()
def med(m, e):
    ok, c = arr(m, e); return np.median(c[ok]) if ok.any() else np.nan


def paired(a, b, e, nb=2000, seed=0):
    """Différence de coût (b - a) sur les départs réussis par les deux, IC bootstrap ; + différence de succès."""
    oa, ca = arr(a, e); ob, cb = arr(b, e); both = oa & ob
    rng = np.random.default_rng(seed)
    if both.sum() < 5: return np.nan, np.nan, np.nan, both.sum(), (oa.mean() - ob.mean())
    d = np.log(cb[both] / ca[both])
    bs = [np.mean(rng.choice(d, len(d))) for _ in range(nb)]
    return float(np.mean(d)), *np.percentile(bs, [2.5, 97.5]), int(both.sum()), float(oa.mean() - ob.mean())


style()
X = [DEF[e] for e in EPS]
groups = [("fd-full", "fd-conf3", "fd-confproj"), ("sec-rank1", "sec-multi", "sec-conf")]
cols = [PAL["blue"], PAL["orange"], PAL["aqua"]]
fig, ax = plt.subplots(2, 2, figsize=(10, 7), sharex=True)
for j, B in enumerate((24, 36)):
    for gi, g in enumerate(groups):
        for k, m in enumerate(g):
            ax[0, j].plot(X, [100 * rate(f"{m}@{B}", e) for e in EPS], "o-", ms=4, color=cols[k], ls=["-", "--"][gi], label=m)
            ax[1, j].plot(X, [med(f"{m}@{B}", e) for e in EPS], "o-", ms=4, color=cols[k], ls=["-", "--"][gi])
    for m, c in (("broyden", PAL["violet"]), ("newton-Jexacte", PAL["red"])):
        ax[0, j].plot(X, [100 * rate(m, e) for e in EPS], ":", color=c, label=m + " (exact)")
        ax[1, j].plot(X, [med(m, e) for e in EPS], ":", color=c)
    ax[0, j].set_title(f"Succès (%), S2 sphérique {B} bits"); ax[1, j].set_title("Lectures médianes")
    ax[1, j].set_xlabel("défaut de conformité au zéro")
ax[0, 0].legend(fontsize=7, ncol=2)
fig.suptitle("Partie 1.3 — Möbius dégradé : conforme (traits) contre générique", x=0.01, ha="left")
fig.tight_layout(); fig.savefig(os.path.join(HERE, "figures", "01_degradation.png")); plt.close(fig)

pairs = [("fd-full", "fd-conf3"), ("fd-full", "fd-confproj"), ("sec-rank1", "sec-conf"), ("sec-rank1", "sec-multi")]
fig, ax = plt.subplots(1, 2, figsize=(10, 3.8), sharey=True)
tab = []
for j, B in enumerate((24, 36)):
    for (a, b), c in zip(pairs, [PAL["blue"], PAL["orange"], PAL["aqua"], PAL["violet"]]):
        v = [paired(f"{a}@{B}", f"{b}@{B}", e) for e in EPS]
        y = [100 * (1 - np.exp(x[0])) for x in v]      # gain en % de lectures de b par rapport à a (positif = b moins cher)
        lo = [100 * (1 - np.exp(x[2])) for x in v]; hi = [100 * (1 - np.exp(x[1])) for x in v]
        ax[j].errorbar(X, y, yerr=[np.array(y) - lo, np.array(hi) - y], fmt="o-", ms=4, color=c, capsize=2, label=f"{b} vs {a}")
        for e, x in zip(EPS, v): tab.append((B, f"{b} vs {a}", e, DEF[e], x))
    ax[j].axhline(0, color=INK2, lw=0.8); ax[j].set_title(f"{B} bits"); ax[j].set_xlabel("défaut de conformité au zéro")
ax[0].set_ylabel("gain en lectures (%), départs réussis par les deux"); ax[0].legend(fontsize=7)
fig.tight_layout(); fig.savefig(os.path.join(HERE, "figures", "02_gain_apparie.png")); plt.close(fig)

L = ["# Partie 1.3 — tableaux (générés par bench_mobius_figs.py)\n", f"Budget {BUD} lectures ; {len(SEEDS)} instances × {R['ntest']} départs par ε. Défaut = moyenne du défaut de conformité au zéro.\n"]
for B in (24, 36):
    L.append(f"\n## Succès % / lectures médianes, {B} bits\n")
    ms = [f"{m}@{B}" for g in groups for m in g] + ["broyden", "newton-FD", "newton-Jexacte", "bfgs-FD", "nelder-mead"]
    L.append(md_table(["méthode"] + [f"ε={e} (déf {DEF[e]:.2f})" for e in EPS],
                      [[m] + [f"{100*rate(m,e):.0f} / {med(m,e):.0f}" for e in EPS] for m in ms]))
L.append("\n## Gain apparié (en % de lectures, + = la 2ᵉ méthode est moins chère), IC95 % bootstrap, n = départs réussis par les deux\n")
L.append(md_table(["bits", "comparaison", "ε", "défaut", "gain %", "IC95", "n", "Δ succès (pts)"],
                  [[B, c, e, f"{d:.2f}", "—" if np.isnan(x[0]) else f"{100*(1-np.exp(x[0])):.0f}",
                    "—" if np.isnan(x[0]) else f"[{100*(1-np.exp(x[2])):.0f}, {100*(1-np.exp(x[1])):.0f}]", x[3], f"{-100*x[4]:+.0f}"]
                   for B, c, e, d, x in tab]))
open(os.path.join(HERE, "bench_mobius_tables.md"), "w").write("\n".join(L))
print("\n".join(L))
