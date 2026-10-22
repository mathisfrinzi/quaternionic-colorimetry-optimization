"""Partie 4.2 : tableaux et figures des réservoirs (point fixe F(h)=h-σ(Wh+b))."""
import os, sys, json
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "..", "src"))
import matplotlib.pyplot as plt
from report import PAL, style, md_table, nom
R = json.load(open(os.path.join(HERE, "bench_reservoirs_results.json"))); style()
os.makedirs(os.path.join(HERE, "figures"), exist_ok=True)
A = [r for r in R if r["part"] == "A"]; B = [r for r in R if r["part"] == "B"]
MS = ["sec-H", "sec-Hid", "sec-full", "broyden", "broyden-id", "picard"]
COL = dict(zip(MS, [PAL["green"], PAL["violet"], PAL["aqua"], PAL["orange"], PAL["red"], PAL["blue"]]))
LAB = {m: nom(m) for m in MS}
KIND = {"split": "tanh par composante", "radial": "tanh radial", "global": "saturation globale"}
def cell(r): return "—" if r is None else f"{100*r['rate']:.0f} % / {r['cost']:.0f}" if r["cost"] == r["cost"] and r["cost"] is not None else f"{100*r['rate']:.0f} % / —"
def get(P, **kw):
    for r in P:
        if all(r[k] == v for k, v in kw.items()): return r
L = ["# Partie 4.2 — tableaux\n", "Cellules : succès % / lectures médianes (succès seulement). Défaut = défaut ℍ-linéaire de la jacobienne au point fixe.\n"]
for rho in (0.8, 0.95):
    L.append(f"\n## A. ρ = {rho}, β = 1\n")
    rows = []
    for n in (2, 8, 16, 32):
        for k in KIND:
            d = get(A, n=n, rho=rho, kind=k, method="broyden"); 
            if d is None: continue
            rows.append([n, KIND[k], f"{d['defect']:.2f}"] + [cell(get(A, n=n, rho=rho, kind=k, method=m)) for m in MS])
    L.append(md_table(["n", "activation", "défaut"] + [LAB[m] for m in MS], rows))
L.append("\n## B. défaut contrôlé (n = 8, ρ = 0,95)\n")
rows = []
for key in sorted({(r["kind"], r["beta"], r["t_real"]) for r in B}, key=lambda x: (x[2], x[0], x[1])):
    d = get(B, kind=key[0], beta=key[1], t_real=key[2], method="broyden")
    rows.append([KIND[key[0]], key[1], key[2], f"{d['defect']:.2f}"] + [cell(get(B, kind=key[0], beta=key[1], t_real=key[2], method=m)) for m in MS])
L.append(md_table(["activation", "β", "part réelle t", "défaut"] + [LAB[m] for m in MS], rows))
open(os.path.join(HERE, "bench_reservoirs_tables.md"), "w").write("\n".join(L))
# Figure 1 : lectures en fonction de n (ρ=0.95, β=1)
fig, axs = plt.subplots(1, 3, figsize=(11, 3.6), sharey=True)
for ax, k in zip(axs, KIND):
    for m in MS:
        xs = [n for n in (2, 8, 16, 32) if get(A, n=n, rho=0.95, kind=k, method=m)]
        ys = [get(A, n=n, rho=0.95, kind=k, method=m)["cost"] for n in xs]
        ax.plot(xs, ys, "o-", color=COL[m], label=LAB[m], lw=1.6, ms=4)
    ax.set_xscale("log", base=2); ax.set_yscale("log"); ax.set_xticks([2, 8, 16, 32]); ax.set_xticklabels([2, 8, 16, 32])
    ax.set_title(KIND[k]); ax.set_xlabel("n (dimension réelle 4n)")
axs[0].set_ylabel("lectures de couleur (médiane)"); axs[0].legend(fontsize=7.5, loc="upper left")
fig.suptitle("Réservoir à point fixe : coût en lectures (ρ = 0,95, β = 1)", x=0.01, ha="left", y=1.02)
fig.tight_layout(); fig.savefig(os.path.join(HERE, "figures", "02_lectures_vs_n.png"), bbox_inches="tight"); plt.close(fig)
# Figure 2 : lectures et succès vs défaut contrôlé (n=8)
fig, axs = plt.subplots(1, 2, figsize=(10, 3.8))
for m in ["sec-H", "sec-Hid", "broyden-id", "picard"]:
    pts = sorted((r["defect"], r["cost"] if r["cost"] == r["cost"] else np.nan, 100 * r["rate"]) for r in B if r["method"] == m)
    x, c, s = zip(*pts)
    axs[0].plot(x, c, "o", color=COL[m], label=LAB[m], ms=4, alpha=.8); axs[1].plot(x, s, "o", color=COL[m], ms=4, alpha=.8)
axs[0].set_yscale("log"); axs[0].set_xlabel("défaut ℍ-linéaire mesuré"); axs[0].set_ylabel("lectures (médiane)"); axs[0].legend(fontsize=8)
axs[1].set_xlabel("défaut ℍ-linéaire mesuré"); axs[1].set_ylabel("succès (%)")
fig.suptitle("Défaut contrôlé (β, activation, part réelle de W) : n = 8, ρ = 0,95", x=0.01, ha="left", y=1.02)
fig.tight_layout(); fig.savefig(os.path.join(HERE, "figures", "03_defaut_controle.png"), bbox_inches="tight"); plt.close(fig)
print("ok")
