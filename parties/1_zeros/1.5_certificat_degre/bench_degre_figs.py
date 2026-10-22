"""Partie 1.5 : figures à partir de bench_degre_results.json."""
import os, sys, json
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "..", "src"))
import matplotlib.pyplot as plt
from report import PAL, style, md_table
R = json.load(open(os.path.join(HERE, "bench_degre_results.json")))
style(); os.makedirs(os.path.join(HERE, "figures"), exist_ok=True)
A = [r for r in R["a"] if r["sigma"] == 0]
bases = list(A[0]["per_base"])
fig, ax = plt.subplots(1, 3, figsize=(13, 3.8))
for d, c in ((1, PAL["orange"]), (2, PAL["blue"])):
    for enc, ls in (("sph", "-"), ("rgb", "--")):
        r = sorted([x for x in A if x["depth"] == d and x["enc"] == enc], key=lambda x: x["B"])
        ax[0].plot([x["B"] for x in r], [100 * x["rate"] for x in r], "o", ls=ls, color=c, ms=4, label=f"profondeur {d} ({r[0]['reads']} lectures), {enc}")
ax[0].set_xlabel("bits par lecture"); ax[0].set_ylabel("degré exact (%)"); ax[0].set_title("Degré retrouvé, sans bruit"); ax[0].legend(fontsize=7)
N = [r for r in R["a"] if r["B"] == 48]
for d, c in ((1, PAL["orange"]), (2, PAL["blue"])):
    for enc, ls in (("sph", "-"), ("rgb", "--")):
        r = sorted([x for x in N if x["depth"] == d and x["enc"] == enc], key=lambda x: x["sigma"])
        ax[1].plot([x["sigma"] for x in r], [100 * x["rate"] for x in r], "o", ls=ls, color=c, ms=4)
ax[1].set_xlabel("bruit σ par canal (rad / népers)"); ax[1].set_title("Robustesse au bruit (48 bits)")
rb = R["b"]
for i, (r, c) in enumerate(zip(rb, [PAL["blue"], PAL["aqua"], PAL["violet"], PAL["orange"], PAL["magenta"]])):
    ax[2].plot(r["d"], [100 * a for a in r["res"]["q"]["acc"]], "o-", ms=4, color=c, label=f"{r['enc']} {r['B']}b σ={r['sigma']} prof. {r['depth']}")
ax[2].axvline(1, color="#8a8984", lw=0.8); ax[2].set_xlabel("distance |x*| du zéro au centre (rayon 1)"); ax[2].set_title("Zéro dedans / dehors : bonne classification (%)"); ax[2].legend(fontsize=6)
fig.tight_layout(); fig.savefig(os.path.join(HERE, "figures", "01_degre.png")); plt.close(fig)
rows = [[r["enc"], r["B"], r["sigma"], r["depth"], r["reads"], f"{100*r['rate']:.0f}", *[f"{100*r['per_base'][b]:.0f}" for b in bases]] for r in R["a"]]
open(os.path.join(HERE, "bench_degre_tables.md"), "w").write("# Partie 1.5 — degré exact (%) par fonction de base\n\n" + md_table(["axe", "bits", "σ", "prof.", "lectures", "total"] + bases, rows))
