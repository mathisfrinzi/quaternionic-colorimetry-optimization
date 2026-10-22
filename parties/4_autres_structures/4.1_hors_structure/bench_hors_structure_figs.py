"""Partie 4.1 : figures et tableaux. Gains appariés (mêmes instances, mêmes départs, même ordre) en fonction du défaut ℍ-linéaire."""
import os, sys, json
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "..", "src"))
import matplotlib.pyplot as plt
from report import PAL, style, md_table, INK2
D = json.load(open(os.path.join(HERE, "bench_hors_structure_results.json"))); R, DEF = D["results"], D["defects"]; style()
os.makedirs(os.path.join(HERE, "figures"), exist_ok=True)
FAMS = [f for f in ["hlin", "bilat", "poly", "qnn", "qres", "gen"] if any(r["fam"] == f for r in R)]
NS = sorted({r["n"] for r in R}); CASES = sorted({(r["bits"], r["sigma"]) for r in R})
MS = ["sec-full", "sec-H", "sec-Hinit", "sec-H2full", "broyden", "broyden-deg"]
LAB = dict(hlin="ℍ-linéaire (témoin)", bilat="bilatère a·d·b", poly="polynôme non commutatif", qnn="inversion réseau quaternionique",
           qres="point fixe réservoir quaternionique", gen="générique (témoin)")
col = dict(zip(MS, [PAL["violet"], PAL["green"], PAL["aqua"], PAL["blue"], PAL["orange"], PAL["red"]]))
def get(f, n, c, m):
    for r in R:
        if r["fam"] == f and r["n"] == n and (r["bits"], r["sigma"]) == c and r["method"] == m: return r
def gain(a, b, rng=np.random.default_rng(0)):
    if a is None or b is None: return None
    oa, ca = np.array(a["ok"]), np.array(a["costs"]); ob, cb = np.array(b["ok"]), np.array(b["costs"]); both = oa & ob
    if both.sum() < 5: return None
    d = np.log(cb[both] / ca[both]); bs = [np.mean(rng.choice(d, len(d))) for _ in range(1000)]
    lo, hi = np.percentile(bs, [2.5, 97.5]); return 100 * (1 - np.exp(d.mean())), 100 * (1 - np.exp(hi)), 100 * (1 - np.exp(lo))
L = ["# Partie 4.1 — tableaux\n", "Défaut ℍ-linéaire de la jacobienne au zéro (moyenne des instances) :\n",
     md_table(["famille"] + [f"n={n}" for n in NS], [[LAB[f]] + [f"{np.mean(DEF[f'{f}/{n}']):.2f}" if f"{f}/{n}" in DEF else "—" for n in NS] for f in FAMS])]
for c in CASES:
    L.append(f"\n## {c[0]} bits par composante, bruit {c[1]} — succès % / lectures médianes\n")
    rows = []
    for f in FAMS:
        for m in MS:
            cells = []
            for n in NS:
                r = get(f, n, c, m)
                cells.append("—" if r is None else (f"{100*r['rate']:.0f} / {r['cost']:.0f}" if r["cost"] else f"{100*r['rate']:.0f} / —"))
            rows.append([LAB[f], m] + cells)
    L.append(md_table(["famille", "méthode"] + [f"n={n}" for n in NS], rows))
    L.append(f"\n### Gains appariés contre Broyden à lectures EXACTES ({c[0]} bits, σ = {c[1]}) ; + = la sonde est moins chère ; Δ succès en points\n")
    rows = []
    for f in FAMS:
        for m in MS[:4]:
            cells = []
            for n in NS:
                a, b = get(f, n, c, "broyden"), get(f, n, c, m); g = gain(a, b)
                ds = "" if a is None or b is None else f" ({100*(b['rate']-a['rate']):+.0f} pts)"
                cells.append("—" if g is None else f"{g[0]:+.0f} [{g[1]:+.0f}, {g[2]:+.0f}]{ds}")
            rows.append([LAB[f], m] + cells)
    L.append(md_table(["famille", "sonde"] + [f"n={n}" for n in NS], rows))
open(os.path.join(HERE, "bench_hors_structure_tables.md"), "w").write("\n".join(L)); print("\n".join(L))
# figure 1 : gain contre Broyden exact en fonction du défaut (points, pas de lignes ; ▼ = hors cadre ou échec)
SHORT = dict(hlin="ℍ-lin", bilat="bilatère", poly="polynôme", qnn="réseau", qres="réservoir", gen="générique")
off = dict(zip(MS[:4], [-0.018, -0.006, 0.006, 0.018]))
for c in CASES:
    fig, ax = plt.subplots(1, len(NS), figsize=(4.4 * len(NS), 4.3), sharey=True)
    ax = np.atleast_1d(ax)
    for k, n in enumerate(NS):
        for f in FAMS:
            if f"{f}/{n}" not in DEF: continue
            x0 = np.mean(DEF[f"{f}/{n}"])
            ax[k].text(x0 + dict(hlin=-0.035, poly=0.035).get(f, 0), 72, SHORT[f], rotation=90, fontsize=7.5, color=INK2, ha="center", va="bottom")
            for m in MS[:4]:
                g = gain(get(f, n, c, "broyden"), get(f, n, c, m)); x = x0 + off[m]
                if g is None or g[0] < -200:
                    ax[k].plot(x, -195, "v", color=col[m], ms=6)
                else:
                    ax[k].errorbar(x, g[0], yerr=[[g[0] - g[1]], [g[2] - g[0]]], fmt="o", ms=5, capsize=2, color=col[m])
        ax[k].axhline(0, color=INK2, lw=0.8); ax[k].set_title(f"n = {n}  (dimension {4*n})", pad=34)
        ax[k].set_xlabel("défaut ℍ-linéaire de la jacobienne au zéro"); ax[k].set_ylim(-205, 70); ax[k].set_xlim(-0.06, 1.0)
    for m in MS[:4]: ax[0].plot([], [], "o", color=col[m], label=m)
    ax[0].plot([], [], "v", color=INK2, label="< −200 % ou échec"); ax[0].legend(fontsize=7, loc="lower left")
    ax[0].set_ylabel("gain en lectures contre Broyden exact (%)")
    fig.suptitle(f"Partie 4.1 — la sonde de couleur contre Broyden à lectures exactes, selon le défaut ℍ-linéaire ({c[0]} bits, σ = {c[1]})", x=0.01, ha="left")
    fig.tight_layout(); fig.savefig(os.path.join(HERE, "figures", f"01_gain_vs_defaut_sigma{c[1]}.png")); plt.close(fig)
# figure 2 : lectures médianes par famille au plus grand n
c = CASES[0]; n = NS[-1]
WRAP = dict(hlin="ℍ-linéaire\n(témoin)", bilat="bilatère\na·d·b", poly="polynôme\nnon commutatif", qnn="inversion\nréseau", qres="point fixe\nréservoir", gen="générique\n(témoin)")
fig, ax = plt.subplots(figsize=(11, 4.4)); w = 0.13; x = np.arange(len(FAMS))
for i, m in enumerate(MS):
    for xi, f in zip(x, FAMS):
        r = get(f, n, c, m)
        if r is None: continue
        if not r["cost"]:
            ax.text(xi + i * w, 11, "0 %", fontsize=6, ha="center", rotation=90, color=col[m]); continue
        ax.bar(xi + i * w, r["cost"], w, color=col[m], label=m if xi == 0 else None)
        if r["rate"] < 0.9: ax.text(xi + i * w, r["cost"] * 1.05, f"{100*r['rate']:.0f} %", fontsize=6, ha="center", rotation=90)
ax.set_xticks(x + 2.5 * w); ax.set_xticklabels([WRAP[f] for f in FAMS], fontsize=8.5); ax.set_yscale("log"); ax.set_ylim(10, 500)
ax.legend(fontsize=7, ncol=6, loc="upper left")
ax.set_title(f"Lectures médianes par famille, n = {n} (dimension {4*n}), {c[0]} bits, sans bruit (succès affiché s'il est < 90 %)")
fig.tight_layout(); fig.savefig(os.path.join(HERE, "figures", "02_par_famille.png")); plt.close(fig)
