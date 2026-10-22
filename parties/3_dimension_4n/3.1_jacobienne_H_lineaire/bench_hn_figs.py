"""Partie 3.1 : figures et tableaux (ℍⁿ). Comparaisons appariées : mêmes instances, mêmes départs, même ordre."""
import os, sys, json
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "..", "src"))
import matplotlib.pyplot as plt
from report import PAL, style, md_table
R = json.load(open(os.path.join(HERE, "bench_hn_results.json"))); style()
os.makedirs(os.path.join(HERE, "figures"), exist_ok=True)
NS = sorted({r["n"] for r in R}); EPS = sorted({r["eps"] for r in R}); CASES = sorted({(r["bits"], r["sigma"]) for r in R})
MS = ["fd-full", "fd-H", "sec-full", "sec-H", "broyden", "broyden-deg", "newton-FD", "bfgs-FD"]
col = dict(zip(MS, [PAL["blue"], PAL["aqua"], PAL["violet"], PAL["green"], PAL["orange"], PAL["red"], PAL["yellow"], PAL["magenta"]]))
ls = {"fd-full": "--", "fd-H": "-", "sec-full": "--", "sec-H": "-", "broyden": ":", "broyden-deg": "-.", "newton-FD": ":", "bfgs-FD": ":"}
def get(n, e, c, m):
    for r in R:
        if r["n"] == n and r["eps"] == e and (r["bits"], r["sigma"]) == c and r["method"] == m: return r
def gain(a, b, rng=np.random.default_rng(0)):
    oa, ca = np.array(a["ok"]), np.array(a["costs"]); ob, cb = np.array(b["ok"]), np.array(b["costs"]); both = oa & ob
    if both.sum() < 5: return None
    d = np.log(cb[both] / ca[both]); bs = [np.mean(rng.choice(d, len(d))) for _ in range(1000)]
    lo, hi = np.percentile(bs, [2.5, 97.5]); return 100 * (1 - np.exp(d.mean())), 100 * (1 - np.exp(hi)), 100 * (1 - np.exp(lo))
for c in CASES:
    fig, ax = plt.subplots(2, len(EPS), figsize=(4 * len(EPS), 7), sharex=True)
    for k, e in enumerate(EPS):
        for m in MS:
            rs = [get(n, e, c, m) for n in NS]
            ax[0, k].plot(NS, [r["cost"] if r and r["cost"] else np.nan for r in rs], "o", ls=ls[m], color=col[m], ms=4, label=m)
            ax[1, k].plot(NS, [100 * r["rate"] if r else np.nan for r in rs], "o", ls=ls[m], color=col[m], ms=4)
        ax[0, k].set_yscale("log"); ax[0, k].set_title(f"ε = {e} — lectures médianes"); ax[1, k].set_title("succès (%)")
        ax[1, k].set_xlabel("n  (dimension réelle 4n)"); ax[1, k].set_xticks(NS)
    ax[0, 0].legend(fontsize=7, ncol=2)
    fig.suptitle(f"Partie 3.1 — ℍⁿ → ℍⁿ, S2 sphérique {c[0]} bits par composante, bruit {c[1]}", x=0.01, ha="left")
    fig.tight_layout(); fig.savefig(os.path.join(HERE, "figures", f"01_hn_{c[0]}b_sigma{c[1]}.png")); plt.close(fig)
L = ["# Partie 3.1 — tableaux (succès % / lectures médianes)\n"]
for c in CASES:
    L.append(f"\n## {c[0]} bits par composante, bruit {c[1]}\n")
    rows = []
    for e in EPS:
        for m in MS:
            rows.append([e, m] + [f"{100*get(n,e,c,m)['rate']:.0f} / {get(n,e,c,m)['cost']:.0f}" if get(n, e, c, m)["cost"] else f"{100*get(n,e,c,m)['rate']:.0f} / —" for n in NS])
    L.append(md_table(["ε", "méthode"] + [f"n={n}" for n in NS], rows))
    L.append(f"\n### Gains appariés ({c[0]} bits, σ = {c[1]}) en % de lectures (+ = la 1ʳᵉ méthode est moins chère), IC 95 %\n")
    rows = []
    for e in EPS:
        for a, b in (("fd-full", "fd-H"), ("sec-full", "sec-H"), ("broyden-deg", "sec-H"), ("broyden", "sec-H")):
            cells = []
            for n in NS:
                g = gain(get(n, e, c, a), get(n, e, c, b))
                cells.append("—" if g is None else f"{g[0]:+.0f} [{g[1]:+.0f}, {g[2]:+.0f}]")
            rows.append([e, f"{b} contre {a}"] + cells)
    L.append(md_table(["ε", "comparaison"] + [f"n={n}" for n in NS], rows))
open(os.path.join(HERE, "bench_hn_tables.md"), "w").write("\n".join(L)); print("\n".join(L))
