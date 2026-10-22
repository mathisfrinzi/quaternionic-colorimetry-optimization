"""Partie 1.1 : tableaux et figures à partir de bench_tranches_results.json."""
import os, sys, json
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", "..", ".."))
sys.path[:0] = [os.path.join(ROOT, "src"), HERE]
import bench, report
from report import PAL, INK2, plt
from funcs import SLICE_SPECS, slice_problem
from sensors import Sensor
from matplotlib.colors import hsv_to_rgb

report.style()
res = json.load(open(os.path.join(HERE, "bench_tranches_results.json")))
keys = list(res["data"]["S1-signé-24b"])
fnames = [json.loads(k)[1] for k in keys]
methods = list(res["data"])
INFO = {  # bits par lecture, lectures par pas, famille
    "complexe-original-8b": (24, "2", "complexe"), "complexe-original-12b": (36, "2", "complexe"),
    "S1-signé-16b": (16, "2", "S1"), "S1-signé-24b": (24, "2", "S1"), "S1-non-signé-24b": (24, "2", "S1"),
    "S2-4D-24b": (24, "5", "S2"), "newton-complexe": ("exacte", "2 (f + f')", "classique"),
    "newton-4D-FD": ("exacte", "5", "classique"), "newton-4D-Jexacte": ("exacte", "2 (F + J)", "classique"),
    "broyden-4D": ("exacte", "1 (+5 init.)", "classique"), "bfgs-4D": ("exacte", "5", "classique"),
    "nelder-mead-4D": ("exacte", "1", "classique"), "cma-es-4D": ("exacte", "1", "classique")}
COL = {"complexe": PAL["violet"], "S1": PAL["orange"], "S2": PAL["blue"], "classique": "#52514e"}
TARGET = {n: SLICE_SPECS[n][2] for n in SLICE_SPECS}

# ------------------------------------------------------------------ tableau principal
rows, S = [], {}
for m in methods:
    s = bench.summarize(res, m); S[m] = s
    b, rp, fam = INFO[m]
    rows.append([m, b, rp, report.pct(s), report.fmt_cost(s["cost"])])
tab = report.md_table(["méthode", "bits / lecture", "lectures / pas", "succès (IC 95 %)", "lectures médianes (succès)"], rows)
common = ["tanh", "z3-1", "sinz-z", "exp-z-1", "rationnelle", "sigmoide-0.5"]
ck = [k for k, n in zip(keys, fnames) if n in common]
# résultats publiés du dépôt complexe, s'il est cloné à côté de celui-ci (variable COMPLEX_REPO pour un autre emplacement)
_ref = os.path.join(os.environ.get("COMPLEX_REPO", os.path.join(ROOT, "..", "..", "complex-colorimetry-optimization")),
                    "parties", "2_zeros_en_2d", "bench2d_results.json")
ref = json.load(open(_ref)) if os.path.exists(_ref) else None
rows2 = []
for m in ["complexe-original-8b", "complexe-original-12b", "S1-signé-16b", "S1-signé-24b"]:
    s = bench.summarize(res, m, ck); rows2.append([m, report.pct(s), report.fmt_cost(s["cost"])])
if ref:
    for m, lab in [("sonde-aveugle-8b", "repo complexe : sonde-aveugle-8b (publié)"), ("sonde-aveugle-12b", "repo complexe : sonde-aveugle-12b (publié)")]:
        ok = np.concatenate([ref["data"][m][f]["ok"] for f in common]); c = np.concatenate([ref["data"][m][f]["cost"] for f in common])
        k, n = int(np.sum(ok)), len(ok); lo, hi = bench.wilson(k, n)
        rows2.append([lab, f"{100*k/n:.0f} % ({100*lo:.0f}–{100*hi:.0f})", f"{np.median(c[np.array(ok)]):.0f}"])
tab2 = report.md_table(["méthode (6 fonctions identiques au repo complexe)", "succès (IC 95 %)", "lectures médianes"], rows2)
M = np.array([[100 * np.mean(res["data"][m][k]["ok"]) for k in keys] for m in methods])
tab3 = report.md_table(["méthode"] + [f"{n}<br>({TARGET[n]})" for n in fnames], [[m] + [f"{v:.0f}" for v in M[i]] for i, m in enumerate(methods)])
open(os.path.join(HERE, "bench_tranches_tables.md"), "w").write(
    f"## Tableau principal (budget {res['budget']} lectures, 10 fonctions × {res['ntest']} départs)\n\n{tab}\n"
    f"## Comparaison à fonctions identiques avec le repo complexe\n\n{tab2}\n"
    f"## Succès (%) par fonction — type de cible entre parenthèses\n\n{tab3}")

# ------------------------------------------------------------------ figure 1 : bilan
order = sorted(methods, key=lambda m: -S[m]["rate"])
fig, (a, b) = plt.subplots(1, 2, figsize=(11, 5.2), sharey=True)
for i, m in enumerate(order):
    s, fam = S[m], INFO[m][2]
    a.errorbar(100 * s["rate"], i, xerr=[[100 * (s["rate"] - s["lo"])], [100 * (s["hi"] - s["rate"])]], fmt="o", color=COL[fam], ms=6, capsize=2)
    b.plot(s["cost"], i, "o", color=COL[fam], ms=6)
a.set_yticks(range(len(order))); a.set_yticklabels(order); a.invert_yaxis()
a.set_xlabel("succès (%), IC de Wilson 95 %"); b.set_xlabel("lectures médianes parmi les succès"); b.set_xscale("log")
a.set_title("Succès"); b.set_title("Coût")
for fam, lab in [("complexe", "code complexe du repo"), ("S1", "capteur S1 (quaternionique)"), ("S2", "capteur S2 (sonde 4D)"), ("classique", "classique, lecture exacte")]:
    a.plot([], [], "o", color=COL[fam], label=lab)
h_, l_ = a.get_legend_handles_labels()
fig.suptitle("Partie 1.1 — tranche C_i : la sonde complexe se transpose sur le capteur quaternionique S1", x=0.01, ha="left", fontsize=12)
fig.legend(h_, l_, loc="lower center", ncol=4, fontsize=8)
fig.tight_layout(rect=(0, 0.05, 1, 0.95)); fig.savefig(os.path.join(HERE, "figures", "01_bilan_tranches.png")); plt.close(fig)

# ------------------------------------------------------------------ figure 2 : par fonction
fig, ax = plt.subplots(figsize=(11, 5))
report.heatmap(ax, M, methods, [f"{n}\n{TARGET[n]}" for n in fnames])
ax.set_title("Succès (%) par fonction et par méthode")
fig.tight_layout(); fig.savefig(os.path.join(HERE, "figures", "02_par_fonction.png")); plt.close(fig)

# ------------------------------------------------------------------ figure 3 : la tranche cache la sphère (visualisation 3D)
from probes import probe_fd, Meter, Hit, OutOfBudget
P1 = slice_problem("q2+1"); P2 = slice_problem("z3-1")
chosen = res["chosen"]["S2-4D-24b"]
def finals(P, n=500, seed=3, bits=48):
    out = []
    for q0 in P.starts(n, seed):
        M_ = Meter(P.f, Sensor("s2", bits), budget=800, thr=1e-3)
        try:
            probe_fd(M_, q0, model="full", c=0.3, step=0.5, cap=1.0)
        except Hit:
            out.append(M_.last)
        except OutOfBudget:
            pass
    return np.array(out)
F1 = finals(P1); F2 = finals(P2)
fig = plt.figure(figsize=(12.5, 4.6))
# (a) domain coloring de la tranche
ax = fig.add_subplot(1, 3, 1)
xs = np.linspace(-2, 2, 400); X, Y = np.meshgrid(xs, xs)
Z = X + 1j * Y; Fv = Z ** 2 + 1
hue = (np.angle(Fv) / (2 * np.pi)) % 1; band = np.log10(np.abs(Fv) + 1e-300) % 1
img = hsv_to_rgb(np.stack([hue, 0.85 * np.ones_like(hue), 0.14 + 0.72 * band], -1))
ax.imshow(img, extent=[-2, 2, -2, 2], origin="lower"); ax.grid(False)
ax.plot([0, 0], [1, -1], "wo", ms=7, mec="k")
ax.set_title("(a) la tranche C_i de q²+1 :\ndeux points, la sphère est invisible"); ax.set_xlabel("Re q"); ax.set_ylabel("coef. de i")
# (b) sphère de zéros en 3D (q^2+1 : Re = 0, |v| = 1)
ax = fig.add_subplot(1, 3, 2, projection="3d")
u, v = np.meshgrid(np.linspace(0, 2 * np.pi, 30), np.linspace(0, np.pi, 15))
ax.plot_wireframe(np.cos(u) * np.sin(v), np.sin(u) * np.sin(v), np.cos(v), color="#c9c8c3", lw=0.5)
ax.scatter(F1[:, 1], F1[:, 2], F1[:, 3], s=8, color=PAL["blue"], alpha=0.8)
ax.scatter([1, -1], [0, 0], [0, 0], s=60, color="k", marker="*", label="les 2 points de C_i")
ax.set_title(f"(b) q²+1 : {len(F1)} arrivées de la sonde 4D\nsur la sphère |Im q| = 1 (Re q ≈ {np.abs(F1[:,0]).max():.0e})", fontsize=10)
ax.set_xlabel("i"); ax.set_ylabel("j"); ax.set_zlabel("k"); ax.legend(fontsize=7, loc="upper left")
# (c) z^3 - 1 : un point réel + une sphère (Re = -1/2, |v| = sqrt(3)/2)
ax = fig.add_subplot(1, 3, 3, projection="3d")
near_pt = np.linalg.norm(F2 - np.array([1, 0, 0, 0]), axis=1) < 0.05
r = np.sqrt(3) / 2
ax.plot_wireframe(r * np.cos(u) * np.sin(v), r * np.sin(u) * np.sin(v), r * np.cos(v), color="#c9c8c3", lw=0.5)
ax.scatter(F2[~near_pt, 1], F2[~near_pt, 2], F2[~near_pt, 3], s=8, color=PAL["orange"], alpha=0.8, label=f"sphère ({(~near_pt).sum()})")
ax.scatter([0], [0], [0], s=80, color=PAL["blue"], marker="o", label=f"point q = 1 ({near_pt.sum()} arrivées, hors du plan Re q = -1/2)")
ax.set_title("(c) q³−1 : projection sur Im q\n(le zéro réel q = 1 projette à l'origine)", fontsize=10)
ax.set_xlabel("i"); ax.set_ylabel("j"); ax.set_zlabel("k"); ax.legend(fontsize=7, loc="upper left")
fig.tight_layout(); fig.savefig(os.path.join(HERE, "figures", "03_spheres_de_zeros.png")); plt.close(fig)
json.dump(dict(q2p1_n=int(len(F1)), q2p1_re_max=float(np.abs(F1[:, 0]).max()), q2p1_norm_dev=float(np.abs(np.linalg.norm(F1[:, 1:], axis=1) - 1).max()),
               z3_point=int(near_pt.sum()), z3_sphere=int((~near_pt).sum())), open(os.path.join(HERE, "sphere_stats.json"), "w"))
print(tab); print(tab2); print(open(os.path.join(HERE, "sphere_stats.json")).read())
