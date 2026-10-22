"""Partie 2 (suite) : figure et tableau de l'ablation des canaux."""
import os, sys, json
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "src"))
import matplotlib.pyplot as plt
from report import PAL, style, md_table
D = json.load(open(os.path.join(HERE, "bench_ablation_results.json"))); style()
os.makedirs(os.path.join(HERE, "figures"), exist_ok=True)
fams = ["nozero", "nozero_ill"]; col = {"luminosité": PAL["blue"], "teinte": PAL["orange"], "axe": PAL["aqua"], "équilibré": PAL["violet"]}
fig, ax = plt.subplots(2, 2, figsize=(10, 6.5))
for i, fam in enumerate(fams):
    for ch in ("luminosité", "teinte", "axe"):
        r = sorted([x for x in D if x["fam"] == fam and x["kind"] == "bruit" and x["canal"] == ch], key=lambda x: x["sigma"])
        ax[0, i].plot([x["sigma"] for x in r], [100 * x["rate"] for x in r], "o-", ms=4, color=col[ch], label=f"bruit sur {ch}")
    ax[0, i].set_title(f"{fam} : succès (%) vs bruit d'un seul canal"); ax[0, i].set_xlabel("σ")
    b = [x for x in D if x["fam"] == fam and x["kind"] == "bits"]
    ax[1, i].bar([x["canal"] for x in b], [100 * x["rate"] for x in b], color=[col[x["canal"]] for x in b])
    ax[1, i].set_title(f"{fam} : succès (%) quand un canal est réduit à 2 bits (36 bits au total)")
ax[0, 0].legend(fontsize=7); fig.tight_layout(); fig.savefig(os.path.join(HERE, "figures", "01_ablation_canaux.png"))
rows = [[x["fam"], x["kind"], x.get("canal", x.get("label", "")), x.get("sigma", ""), f"{100*x['rate']:.0f}", "—" if x["cost"] is None else f"{x['cost']:.0f}", x["eta"]] for x in D]
open(os.path.join(HERE, "bench_ablation_tables.md"), "w").write("# Partie 2 (suite) — ablation des canaux (succès %, lectures médianes, η retenu)\n\n" + md_table(["famille", "type", "canal / cas", "σ", "succès %", "lectures", "η"], rows))
