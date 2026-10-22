"""Partie 5 : figure + tableaux + comparaison appariée sonde-rank1 vs Broyden (mêmes départs)."""
import os, sys, json
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "src"))
import matplotlib.pyplot as plt
from report import PAL, style, md_table, set_noms
set_noms(**{"broyden": "Broyden-mêmes-mesures", "broyden+dith": "Broyden-mêmes-mesures + dithering", "newton-FD": "Newton-mêmes-mesures", "bfgs-FD": "BFGS-mêmes-mesures", "nelder-mead": "Nelder-Mead-mêmes-mesures"})  # ici, les classiques lisent la mesure dégradée
R = json.load(open(os.path.join(HERE, "bench_armes_egales_results.json"))); D = R["data"]; style()
BITS, SIG = R["bits"], R["sig"]
M = ["sonde-rank1", "sonde-fd", "broyden", "newton-FD", "bfgs-FD", "nelder-mead"]
col = dict(zip(M, [PAL["blue"], PAL["aqua"], PAL["orange"], PAL["red"], PAL["violet"], PAL["yellow"]]))
def arr(m, B, s):
    d = D[f"{m}@{B}/{s}"]; k = list(d)
    return np.concatenate([d[x]["ok"] for x in k]).astype(bool), np.concatenate([d[x]["cost"] for x in k]).astype(float)
def stat(m, B, s):
    ok, c = arr(m, B, s); return 100 * ok.mean(), (np.median(c[ok]) if ok.any() else np.nan)
fig, ax = plt.subplots(2, 2, figsize=(10, 6.5), sharex=True)
for j, B in enumerate(BITS):
    for m in M:
        ax[0, j].plot(SIG, [stat(m, B, s)[0] for s in SIG], "o-", ms=4, color=col[m], label=m)
        ax[1, j].plot(SIG, [stat(m, B, s)[1] for s in SIG], "o-", ms=4, color=col[m])
    ax[0, j].set_title(f"Succès (%), {B} bits"); ax[1, j].set_title("Lectures médianes"); ax[1, j].set_xlabel("bruit σ par canal")
ax[0, 0].legend(fontsize=7); fig.suptitle("Partie 5 — armes égales : même mesure dégradée pour toutes les méthodes", x=0.01, ha="left")
fig.tight_layout(); fig.savefig(os.path.join(HERE, "figures", "01_armes_egales.png")) if os.makedirs(os.path.join(HERE, "figures"), exist_ok=True) is None else None
rows = [[f"{m}", B, s, f"{stat(m,B,s)[0]:.0f}", f"{stat(m,B,s)[1]:.0f}"] for B in BITS for s in SIG for m in M]
L = ["# Partie 5 — tableaux\n", md_table(["méthode", "bits", "σ", "succès %", "lectures médianes"], rows), "\n## Apparié : sonde-rank1 contre Broyden (gain en % de lectures, + = la sonde est moins chère ; IC95 % bootstrap)\n"]
rng = np.random.default_rng(0); pr = []
for B in BITS:
    for s in SIG:
        oa, ca = arr("broyden", B, s); ob, cb = arr("sonde-rank1", B, s); both = oa & ob
        d = np.log(cb[both] / ca[both]); bs = [np.mean(rng.choice(d, len(d))) for _ in range(2000)]
        lo, hi = np.percentile(bs, [2.5, 97.5]); pr.append([B, s, int(both.sum()), f"{100*(1-np.exp(d.mean())):.0f}", f"[{100*(1-np.exp(hi)):.0f}, {100*(1-np.exp(lo)):.0f}]"])
L.append(md_table(["bits", "σ", "n", "gain %", "IC95"], pr))
open(os.path.join(HERE, "bench_armes_egales_tables.md"), "w").write("\n".join(L)); print(md_table(["bits", "σ", "n", "gain %", "IC95"], pr))
