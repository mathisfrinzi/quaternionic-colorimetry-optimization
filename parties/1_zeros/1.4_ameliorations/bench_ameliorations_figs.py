"""Partie 1.4 : tableaux et figure, comparaisons appariées (mêmes départs, mêmes instances) contre `sonde-base` et contre Broyden."""
import os, sys, json
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "..", "src"))
import matplotlib.pyplot as plt
from report import PAL, style, md_table, set_noms
set_noms(**{"broyden": "Broyden-mêmes-mesures", "broyden+dith": "Broyden-mêmes-mesures + dithering", "newton-FD": "Newton-mêmes-mesures", "bfgs-FD": "BFGS-mêmes-mesures", "nelder-mead": "Nelder-Mead-mêmes-mesures"})  # ici, les classiques lisent la mesure dégradée
R = json.load(open(os.path.join(HERE, "bench_ameliorations_results.json"))); D = R["data"]; style()
BITS, SIG, VAR = R["bits"], R["sig"], R["variants"]
os.makedirs(os.path.join(HERE, "figures"), exist_ok=True)
def arr(m):
    d = D[m]; k = list(d)
    return np.concatenate([d[x]["ok"] for x in k]).astype(bool), np.concatenate([d[x]["cost"] for x in k]).astype(float)
def gain(a, b, rng=np.random.default_rng(0)):
    """gain en % de lectures de b par rapport à a (+ : b moins cher), sur les départs réussis par les deux ; IC bootstrap."""
    oa, ca = arr(a); ob, cb = arr(b); both = oa & ob
    d = np.log(cb[both] / ca[both]); bs = [np.mean(rng.choice(d, len(d))) for _ in range(1000)]
    lo, hi = np.percentile(bs, [2.5, 97.5])
    return 100 * (1 - np.exp(d.mean())), 100 * (1 - np.exp(hi)), 100 * (1 - np.exp(lo)), int(both.sum())
names = [f"sonde-{v}" for v in VAR] + ["broyden", "broyden+dith"]
L = ["# Partie 1.4 — tableaux\n"]
rows = []
for B in BITS:
    for s in SIG:
        for n in names:
            ok, c = arr(f"{n}@{B}/{s}")
            rows.append([B, s, n, f"{100*ok.mean():.0f}", f"{np.median(c[ok]):.0f}" if ok.any() else "—"])
L += ["## Succès % et lectures médianes (zéros isolés + Möbius ε ∈ {0, 1}, 3 instances × 30 départs)\n", md_table(["bits", "σ", "méthode", "succès %", "lectures méd."], rows)]
rows = []
for B in BITS:
    for s in SIG:
        t = f"@{B}/{s}"
        for v in VAR[1:]:
            g = gain("sonde-base" + t, f"sonde-{v}" + t); rows.append([B, s, f"sonde-{v} contre sonde-base", f"{g[0]:+.0f} [{g[1]:+.0f}, {g[2]:+.0f}]", g[3]])
        for n in ("sonde-base", "sonde-dith"):
            g = gain("broyden" + t, n + t); rows.append([B, s, f"{n} contre broyden", f"{g[0]:+.0f} [{g[1]:+.0f}, {g[2]:+.0f}]", g[3]])
L += ["\n## Gains appariés en % de lectures (+ = la 1ʳᵉ méthode nommée est moins chère), IC 95 % bootstrap\n", md_table(["bits", "σ", "comparaison", "gain %", "n"], rows)]
open(os.path.join(HERE, "bench_ameliorations_tables.md"), "w").write("\n".join(L)); print("\n".join(L))
col = dict(zip(names, [PAL["blue"], PAL["aqua"], PAL["orange"], PAL["violet"], PAL["magenta"], PAL["red"], PAL["yellow"]]))
fig, ax = plt.subplots(1, 2, figsize=(11, 4), sharey=False)
x = np.arange(4); w = 0.11
cases = [(B, s) for B in BITS for s in SIG]
for i, n in enumerate(names):
    ax[0].bar(x + i * w, [100 * arr(f"{n}@{B}/{s}")[0].mean() for B, s in cases], w, color=col[n], label=n)
    ax[1].bar(x + i * w, [np.median(arr(f"{n}@{B}/{s}")[1][arr(f"{n}@{B}/{s}")[0]]) for B, s in cases], w, color=col[n])
for a, t in zip(ax, ("Succès (%)", "Lectures médianes (succès)")):
    a.set_xticks(x + 3 * w); a.set_xticklabels([f"{B} b, σ={s}" for B, s in cases]); a.set_title(t)
ax[0].set_ylim(40, 100); ax[0].legend(fontsize=6, ncol=2); fig.tight_layout(); fig.savefig(os.path.join(HERE, "figures", "01_ameliorations.png"))
