"""Partie 4.2 : d'où vient le défaut ℍ-linéaire d'un réservoir, et comment le contrôler (sans aucune sonde).
  1. défaut de la dérivée de l'activation en un quaternion z, selon |z| : split (moyenne et pire direction), radial (formule exacte)
  2. défaut de la jacobienne J = I − Dσ·L(W) au point fixe, selon la force de l'entrée β et le rayon spectral, pour les 3 activations
  3. défaut selon la dimension n (activation globale : décroît en 1/√(4n))
  4. défaut selon la taille t d'une pièce réelle (non quaternionique) ajoutée à W
Usage : python defaut_theorie.py"""
import os, sys, json
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "..", "src"))
import matplotlib.pyplot as plt
from report import PAL, style, INK2
from reservoir import Reservoir, dact_block, block_defect
style()
COL = dict(split=PAL["orange"], radial=PAL["blue"], glob=PAL["aqua"])
out = {}
# 1
rs = np.linspace(0.01, 4, 80); rng = np.random.default_rng(0)
U = rng.standard_normal((400, 4)); U /= np.linalg.norm(U, axis=1, keepdims=True)
sp_mean = [np.mean([block_defect(dact_block(r * u, "split")) for u in U]) for r in rs]
sp_max = [np.max([block_defect(dact_block(r * u, "split")) for u in U]) for r in rs]
rad = [block_defect(dact_block(r * U[0], "radial")) for r in rs]
out["activation"] = dict(r=rs.tolist(), split_mean=sp_mean, split_max=sp_max, radial=rad)
# 2
betas = [0.25, 0.5, 1.0, 2.0, 4.0]; rhos = [0.5, 0.8, 0.95]
grid = {k: [[np.mean([Reservoir(8, rho, b, k, s).defect() for s in range(3)]) for b in betas] for rho in rhos] for k in ("split", "radial", "global")}
out["grille"] = dict(betas=betas, rhos=rhos, defect=grid)
# 3
ns = [1, 2, 4, 8, 16, 32]
dn = {k: [np.mean([Reservoir(n, 0.95, 1.0, k, s).defect() for s in range(3)]) for n in ns] for k in ("split", "radial", "global")}
out["dimension"] = dict(n=ns, defect=dn)
# 4
ts = [0.0, 0.05, 0.1, 0.2, 0.3, 0.6, 1.0]
dt = {k: [np.mean([Reservoir(8, 0.95, 1.0, k, s, t).defect() for s in range(3)]) for t in ts] for k in ("split", "radial", "global")}
out["piece_reelle"] = dict(t=ts, defect=dt)
json.dump(out, open(os.path.join(HERE, "defaut_theorie.json"), "w"))

fig, ax = plt.subplots(1, 4, figsize=(15, 3.8))
ax[0].plot(rs, sp_mean, color=COL["split"], label="split (moyenne des directions)")
ax[0].plot(rs, sp_max, color=COL["split"], ls="--", label="split (pire direction)")
ax[0].plot(rs, rad, color=COL["radial"], label="radiale (exact)")
ax[0].axhline(0.5, color=INK2, lw=0.6, ls=":"); ax[0].axhline(np.sqrt(3) / 2, color=INK2, lw=0.6, ls=":")
ax[0].text(3.9, 0.51, "½", ha="right", color=INK2); ax[0].text(3.9, 0.875, "√3/2", ha="right", color=INK2, va="bottom")
ax[0].set_xlabel("|z| (pré-activation d'un quaternion)"); ax[0].set_title("défaut de la dérivée de l'activation"); ax[0].legend(fontsize=7)
for k, c in (("split", COL["split"]), ("radial", COL["radial"]), ("global", COL["glob"])):
    ax[1].plot(betas, grid[k][2], "o-", color=c, ms=4, label=k)
    ax[2].plot(ns, dn[k], "o-", color=c, ms=4, label=k)
    ax[3].plot(ts, dt[k], "o-", color=c, ms=4, label=k)
ax[2].plot(ns, [dn["global"][3] * np.sqrt(8 / n) for n in ns], ":", color=COL["glob"], lw=1, label="∝ 1/√n")
ax[1].set_xscale("log"); ax[1].set_xticks(betas); ax[1].set_xticklabels([str(b) for b in betas]); ax[1].minorticks_off(); ax[1].set_xlabel("force de l'entrée β"); ax[1].set_title("défaut de J au point fixe (n = 8, rayon 0,95)")
ax[2].set_xscale("log", base=2); ax[2].set_xticks(ns); ax[2].set_xticklabels([str(n) for n in ns]); ax[2].minorticks_off(); ax[2].set_xlabel("n (dimension 4n)"); ax[2].set_title("défaut selon la dimension (β = 1)")
ax[3].set_xlabel("taille t de la pièce réelle ajoutée à W"); ax[3].set_title("pièce apprise non quaternionique (n = 8)")
for a in ax[1:]: a.legend(fontsize=7)
for a in ax: a.set_ylim(0, 0.92)
fig.suptitle("Partie 4.2 — le défaut ℍ-linéaire d'un réservoir quaternionique vient de l'activation, de son point de fonctionnement et des pièces non quaternioniques",
             x=0.01, ha="left", fontsize=10.5)
fig.tight_layout(); fig.savefig(os.path.join(HERE, "figures", "01_defaut_theorie.png"))
print(json.dumps({k: v for k, v in out.items() if k != "activation"}, indent=0)[:1500])
