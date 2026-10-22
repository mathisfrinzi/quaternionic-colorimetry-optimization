"""Partie 3.4 : tableaux et figure à partir de bench34_results.json (et bench34_calib.json s'il existe).
Intervalles à 95 % : Wilson pour les taux de succès, bootstrap (2000 tirages) pour les médianes de lectures.
Différence appariée : même instance, même départ, les deux méthodes réussissent ; médiane des écarts de lectures, IC bootstrap."""
import os, json
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
HERE = os.path.dirname(os.path.abspath(__file__))
RNG = np.random.default_rng(0)
ORDRE = ["sonde-garde-ambiguite", "sonde-depliage-vers-zero", "R1-relecture-naive", "R1-relecture", "R2-par-composante",
         "R1+R2", "sonde-oracle-decade", "broyden-sans-couleur"]
NOMS = {"sonde-garde-ambiguite": "sonde-garde-ambiguïté", "sonde-depliage-vers-zero": "sonde-dépliage-vers-zéro",
        "R1-relecture-naive": "sonde-relecture-naïve", "R1-relecture": "sonde-relecture", "R2-par-composante": "sonde-par-composante",
        "R1+R2": "sonde-relecture-par-composante", "sonde-oracle-decade": "sonde-oracle-décade (triche)",
        "broyden-sans-couleur": "Broyden-sans-couleur"}


def wilson(k, m, z=1.96):
    if m == 0:
        return (np.nan, np.nan)
    p = k / m; d = 1 + z * z / m
    c = (p + z * z / (2 * m)) / d; w = z * np.sqrt(p * (1 - p) / m + z * z / (4 * m * m)) / d
    return (max(0.0, c - w), min(1.0, c + w))


def boot_median(x, B=2000):
    x = np.asarray(x, float)
    if len(x) == 0:
        return (np.nan, np.nan, np.nan)
    if len(x) == 1:
        return (x[0], np.nan, np.nan)
    b = np.median(RNG.choice(x, (B, len(x))), axis=1)
    return (float(np.median(x)), float(np.percentile(b, 2.5)), float(np.percentile(b, 97.5)))


def f1(x):
    return "–" if x is None or not np.isfinite(x) else (f"{x:.0f}" if abs(x) >= 10 else f"{x:.1f}")


def table(res, ref="sonde-oracle-decade"):
    lines = []
    cases = sorted({(r["n"], r["sigma"]) for r in res})
    for n, sg in cases:
        lines += [f"\n### n = {n}, σ = {sg:g}\n",
                  "| méthode | succès [IC 95 %] | lectures médianes [IC 95 %] | écart apparié à l'oracle [IC 95 %] (paires) | rejets ambiguïté / descente par lecture de pas | relectures par pas | pas acceptés avec une décade fausse muette |",
                  "|---|---|---|---|---|---|---|"]
        byk = {}
        for r in res:
            if r["n"] == n and r["sigma"] == sg:
                byk.setdefault(r["meth"], {})[(r["seed"], r["i"])] = r
        for m in [m for m in ORDRE if m in byk]:
            rr = list(byk[m].values()); ok = [r for r in rr if r["ok"]]
            lo, hi = wilson(len(ok), len(rr))
            med, a, b = boot_median([r["cost"] for r in ok])
            if m != ref and ref in byk:
                d = [byk[m][key]["cost"] - byk[ref][key]["cost"] for key in byk[m]
                     if key in byk[ref] and byk[m][key]["ok"] and byk[ref][key]["ok"]]
                dm, da, db = boot_median(d)
                ecart = f"{dm:+.0f} [{f1(da)} ; {f1(db)}] ({len(d)})" if d else "–"
            else:
                ecart = "réf."
            pas = sum(r.get("pas", 0) for r in rr)
            if pas:
                rej = f"{sum(r['rej_amb'] for r in rr)/pas:.2f} / {sum(r['rej_desc'] for r in rr)/pas:.2f}"
                ver = f"{sum(r['verif'] for r in rr)/pas:.2f}"
                acc = sum(r["acc"] for r in rr)
                fm = f"{100*sum(r['faux_muet_acc'] for r in rr)/max(acc,1):.1f} %"
            else:
                rej = ver = fm = "–"
            lines.append(f"| {NOMS.get(m, m)} | {100*len(ok)/len(rr):.0f} % [{100*lo:.0f} ; {100*hi:.0f}] ({len(rr)}) | "
                         f"{f1(med)} [{f1(a)} ; {f1(b)}] | {ecart} | {rej} | {ver} | {fm} |")
    return "\n".join(lines)


def figure(res, path):
    cases = sorted({(r["n"], r["sigma"]) for r in res})
    meths = [m for m in ORDRE if any(r["meth"] == m for r in res)]
    fig, axes = plt.subplots(2, len(cases), figsize=(4.2 * len(cases), 7), sharey="row", squeeze=False)
    cols = plt.cm.tab10(np.arange(len(meths)))
    for c, (n, sg) in enumerate(cases):
        for k, m in enumerate(meths):
            rr = [r for r in res if r["n"] == n and r["sigma"] == sg and r["meth"] == m]
            if not rr:
                continue
            ok = [r for r in rr if r["ok"]]
            p = len(ok) / len(rr); lo, hi = wilson(len(ok), len(rr))
            axes[0, c].bar(k, 100 * p, color=cols[k]); axes[0, c].errorbar(k, 100 * p, [[max(0.0, 100 * (p - lo))], [max(0.0, 100 * (hi - p))]], color="k", capsize=3)
            med, a, b = boot_median([r["cost"] for r in ok])
            if np.isfinite(med):
                axes[1, c].bar(k, med, color=cols[k])
                if np.isfinite(a):
                    axes[1, c].errorbar(k, med, [[max(0.0, med - a)], [max(0.0, b - med)]], color="k", capsize=3)
        axes[0, c].set_title(f"n = {n}, σ = {sg:g}")
        for ax in axes[:, c]:
            ax.set_xticks(range(len(meths))); ax.set_xticklabels([NOMS[m] for m in meths], rotation=60, ha="right", fontsize=8)
    axes[0, 0].set_ylabel("succès (%)"); axes[1, 0].set_ylabel("lectures médianes (succès)")
    fig.tight_layout(); fig.savefig(path, dpi=130); plt.close(fig)


def calib_table(cal):
    rows = {}
    for r in cal:
        rows.setdefault(r["meth"], {}).setdefault(r["sigma"], []).append(r)
    out = ["| configuration | σ = 0 : succès / lectures méd. | σ = 0,02 : succès / lectures méd. |", "|---|---|---|"]
    def cell(rr):
        ok = [r["cost"] for r in rr if r["ok"]]
        return f"{100*len(ok)/len(rr):.0f} % / {f1(np.median(ok)) if ok else '–'}"
    for m, d in sorted(rows.items(), key=lambda kv: -sum(np.mean([r["ok"] for r in v]) for v in kv[1].values())):
        out.append(f"| {m} | {cell(d.get(0.0, [])) if 0.0 in d else '–'} | {cell(d.get(0.02, [])) if 0.02 in d else '–'} |")
    return "\n".join(out)


if __name__ == "__main__":
    md = ["# Partie 3.4 : tableaux (générés par bench34_figs.py)"]
    pc = os.path.join(HERE, "bench34_calib.json")
    if os.path.exists(pc):
        md += ["\n## Calibration (n = 32, graines 300–304 × 6 départs = 30 départs, jamais réutilisées)\n", calib_table(json.load(open(pc)))]
    pr = os.path.join(HERE, "bench34_results.json")
    if os.path.exists(pr):
        res = json.load(open(pr))
        md += ["\n## Banc final (graines 100–109 × 10 départs = 100 départs par case)", table(res)]
        os.makedirs(os.path.join(HERE, "figures"), exist_ok=True)
        figure(res, os.path.join(HERE, "figures", "01_remedes_ambiguite.png"))
    open(os.path.join(HERE, "bench34_tables.md"), "w").write("\n".join(md) + "\n")
    print("\n".join(md))
