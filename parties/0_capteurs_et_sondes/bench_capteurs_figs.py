"""Partie 0 : figures des capteurs (erreur relative sur une colonne de W vs bits ; bruit par canal)."""
import os, sys, json
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "src"))
import matplotlib.pyplot as plt
from report import PAL, style
d = json.load(open(os.path.join(HERE, "bench_capteurs_results.json"))); style()
os.makedirs(os.path.join(HERE, "figures"), exist_ok=True)
fig, ax = plt.subplots(1, 2, figsize=(10, 3.8))
for (k, v), c in zip(d["curves"].items(), [PAL["violet"], PAL["blue"], PAL["aqua"], PAL["orange"]]):
    ax[0].semilogy(d["bits"][:len(v)], [max(x[0], 1e-12) for x in v], "o-", ms=4, color=c, label=k)
ax[0].set_xlabel("bits par lecture"); ax[0].set_ylabel("erreur relative sur une colonne de W"); ax[0].set_title("Précision de lecture vs bits"); ax[0].legend(fontsize=7)
for (k, v), c in zip(d["noise"].items(), [PAL["blue"], PAL["orange"], PAL["aqua"]]):
    ax[1].loglog(d["noise_sigma"], v, "o-", ms=4, color=c, label="bruit sur " + k)
ax[1].set_xlabel("σ du canal bruité"); ax[1].set_title("Sensibilité au bruit (48 bits)"); ax[1].legend(fontsize=7)
fig.tight_layout(); fig.savefig(os.path.join(HERE, "figures", "01_capteurs.png"))
