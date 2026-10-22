"""Partie 1.2 : tableaux et figures (bench_generique_results.json, bench_bits_results.json)."""
import os, sys, json
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", "..", ".."))
sys.path[:0] = [os.path.join(ROOT, "src"), HERE]
import bench, report
from report import PAL, plt
from funcs import GENERIC_KINDS, generic_problem

report.style()
res = json.load(open(os.path.join(HERE, "bench_generique_results.json")))
methods = list(res["data"])
GENERIC_KINDS = [k for k in GENERIC_KINDS if any(json.loads(kk)[1] == k for kk in res["data"][methods[0]])]   # familles présentes dans le JSON
N = report.nom
keys_by = {k: [kk for kk in res["data"][methods[0]] if json.loads(kk)[1] == k] for k in GENERIC_KINDS}
TARGET = {k: generic_problem(k, 0).target for k in GENERIC_KINDS}
isprobe = lambda m: m.startswith("sonde")
S = {m: bench.summarize(res, m) for m in methods}
iso = [k for k in GENERIC_KINDS if k in ("lin", "quad", "cub", "ill")]
iso_keys = sum((keys_by[k] for k in iso), [])

# tableau
rows = []
for m in methods:
    s_iso = bench.summarize(res, m, iso_keys)
    fam = [100 * bench.summarize(res, m, keys_by[k])["rate"] for k in GENERIC_KINDS]
    rows.append([N(m), report.pct(s_iso), report.fmt_cost(s_iso["cost"])] + [f"{v:.0f}" for v in fam])
hdr = ["méthode", "zéro isolé : succès (IC 95 %)", "lectures méd."] + [f"{k}<br>({TARGET[k]})" for k in GENERIC_KINDS]
open(os.path.join(HERE, "bench_generique_tables.md"), "w").write(
    f"Budget {res['budget']} lectures ; 7 familles × 3 instances × {res['ntest']} départs. « zéro isolé » regroupe lin, quad, cub, ill "
    f"(480 départs). Colonnes de droite : succès (%) par famille.\n\n" + report.md_table(hdr, rows))

# figure 1 : carte par famille
M = np.array([[100 * bench.summarize(res, m, keys_by[k])["rate"] for k in GENERIC_KINDS] for m in methods])
fig, ax = plt.subplots(figsize=(10, 7))
report.heatmap(ax, M, [N(m) for m in methods], [f"{k}\n{TARGET[k]}" for k in GENERIC_KINDS])
ax.set_title("Succès (%) par famille de fonctions quelconques : sondes (capteur quantifié) puis méthodes classiques (valeur exacte)")
ax.axhline(sum(isprobe(m) for m in methods) - 0.5, color="k", lw=1.2)
fig.tight_layout(); fig.savefig(os.path.join(HERE, "figures", "01_familles.png")); plt.close(fig)

# figure 2 : succès / coût sur les zéros isolés
order = sorted(methods, key=lambda m: -bench.summarize(res, m, iso_keys)["rate"])
fig, (a, b) = plt.subplots(1, 2, figsize=(11, 6), sharey=True)
for i, m in enumerate(order):
    s = bench.summarize(res, m, iso_keys)
    col = PAL["blue"] if isprobe(m) else "#52514e"
    a.errorbar(100 * s["rate"], i, xerr=[[100 * (s["rate"] - s["lo"])], [100 * (s["hi"] - s["rate"])]], fmt="o", color=col, ms=5, capsize=2)
    if np.isfinite(s["cost"]):
        b.plot(s["cost"], i, "o", color=col, ms=5)
a.set_yticks(range(len(order))); a.set_yticklabels([N(m) for m in order]); a.invert_yaxis()
a.set_xlabel("succès (%), zéros isolés"); b.set_xlabel("lectures médianes (succès)"); b.set_xscale("log")
a.plot([], [], "o", color=PAL["blue"], label="sonde, capteur quantifié"); a.plot([], [], "o", color="#52514e", label="classique, lecture exacte")
fig.legend(*a.get_legend_handles_labels(), loc="lower center", ncol=2, fontsize=9)
fig.suptitle("Partie 1.2 : zéros isolés dans ℍ ≅ ℝ⁴, budget 600 lectures", x=0.01, ha="left")
fig.tight_layout(rect=(0, 0.05, 1, 0.96)); fig.savefig(os.path.join(HERE, "figures", "02_bilan.png")); plt.close(fig)

# figure 3 : balayage des bits
bp = os.path.join(HERE, "bench_bits_results.json")
if os.path.exists(bp):
    bd = json.load(open(bp)); B = bd["bits"]
    fig, (a, b) = plt.subplots(1, 2, figsize=(11, 4.4))
    sty = {"S0": (PAL["aqua"], "capteur-valeurs (4 composantes)"), "S1": (PAL["orange"], "couleur-sans-axe (ℓ, θ)"),
           "S2rgb": (PAL["blue"], "couleur-complète-RVB"), "S2sph": (PAL["violet"], "couleur-complète-angles")}
    for fam, (c, lab) in sty.items():
        a.plot(B, [100 * s["rate"] for s in bd[fam]], "-o", color=c, label=lab, ms=4)
        cs = [s["cost"] for s in bd[fam]]
        b.plot(B, cs, "-o", color=c, ms=4)
    e = bd["exacte"]
    a.axhline(100 * e["rate"], color="#52514e", ls="--", lw=1); a.text(B[0], 100 * e["rate"] + 1.5, "même algorithme, capteur exact", color="#52514e", fontsize=8)
    a.set_xlabel("bits par lecture (total)"); a.set_ylabel("succès (%)"); a.set_xscale("log", base=2); a.set_xticks(B); a.set_xticklabels(B)
    b.set_xlabel("bits par lecture (total)"); b.set_ylabel("lectures médianes (succès)"); b.set_xscale("log", base=2); b.set_xticks(B); b.set_xticklabels(B)
    a.legend(fontsize=8, loc="lower right")
    a.set_title("Succès"); b.set_title("Coût")
    fig.suptitle("Partie 1.2 : à bits égaux, le canal d'axe paie-t-il ses bits ? (zéros isolés, réglage fixe)", x=0.01, ha="left")
    fig.tight_layout(); fig.savefig(os.path.join(HERE, "figures", "03_bits.png")); plt.close(fig)
    rows = [[f"{B_}"] + [f"{100*bd[f][i]['rate']:.0f} % / {report.fmt_cost(bd[f][i]['cost'])}" for f in ("S0", "S1", "S2rgb", "S2sph")] for i, B_ in enumerate(B)]
    open(os.path.join(HERE, "bench_bits_tables.md"), "w").write(
        "Succès / lectures médianes par capteur, à B bits par lecture (zéros isolés, 480 départs, paramètres fixes).\n\n"
        + report.md_table(["bits", "capteur-valeurs", "couleur-sans-axe", "couleur-complète-RVB", "couleur-complète-angles"], rows) + f"\nMême algorithme, capteur exact : {report.pct(e)} en {report.fmt_cost(e['cost'])} lectures.\n")
print(open(os.path.join(HERE, "bench_generique_tables.md")).read())
