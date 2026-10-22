"""Tableaux et figures des bancs des parties 3.5, 4.3, 5.5 et 6 (succès avec IC de Wilson, lectures médianes, gains appariés).
Usage : python src/figs_hn.py [35] [43] [55] [61] [62] [63] [64]   (sans argument : tout ce dont le JSON existe)"""
import os, sys, json
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from report import PAL, style, INK, INK2, MUTED, nom, md_table
import banc_hn as B
style()
P = os.path.join(HERE, "..", "parties")
COL = [PAL["orange"], PAL["blue"], PAL["aqua"], PAL["violet"], PAL["yellow"], PAL["magenta"], PAL["green"], PAL["red"], MUTED]


def fr(x, d=2):
    return f"{x:.{d}f}".replace(".", ",")


def cell(s):
    if s is None: return "—"
    t = f"{100 * s['rate']:.0f} % [{max(0, 100 * s['lo']):.0f} ; {min(100, 100 * s['hi']):.0f}]"
    return t + (f" / {s['cost']:.0f}" if s["cost"] is not None else "")


def barres(ax, cases, meths, res_of, titre, budget_txt=""):
    """Deux barres par case : succès (haut) et lectures médianes (bas) — ici une figure à deux axes."""
    pass


def figure_cases(path, cases, meths, res_of, titre):
    """cases : liste de (étiquette, sélecteur) ; res_of(case, meth) -> résumé. Succès en haut, lectures médianes en bas."""
    fig, (a1, a2) = plt.subplots(2, 1, figsize=(max(7, 1.5 * len(cases) + 3), 6.2), sharex=True)
    w = 0.8 / len(meths); x = np.arange(len(cases))
    for k, m in enumerate(meths):
        S = [res_of(c, m) for c in cases]
        r = [np.nan if s is None else 100 * s["rate"] for s in S]
        lo = [0 if s is None else max(0.0, 100 * (s["rate"] - s["lo"])) for s in S]
        hi = [0 if s is None else max(0.0, 100 * (s["hi"] - s["rate"])) for s in S]
        cst = [np.nan if (s is None or s["cost"] is None) else s["cost"] for s in S]
        a1.bar(x + (k - (len(meths) - 1) / 2) * w, r, w * 0.92, color=COL[k % len(COL)], yerr=[lo, hi], error_kw=dict(lw=0.8, ecolor=INK2), label="modèle étendu (selon la famille)" if m == "étendu" else nom(m))
        a2.bar(x + (k - (len(meths) - 1) / 2) * w, cst, w * 0.92, color=COL[k % len(COL)])
    a1.set_ylabel("succès (%)"); a1.set_ylim(0, 105); a2.set_ylabel("lectures médianes (succès)")
    a2.set_xticks(x); a2.set_xticklabels([c[0] for c in cases], fontsize=8.5)
    a1.legend(fontsize=7.5, ncol=min(3, len(meths)), loc="lower left", bbox_to_anchor=(0, 1.02))
    a1.set_title(titre, pad=48 if len(meths) > 3 else 30)
    fig.tight_layout(); os.makedirs(os.path.dirname(path), exist_ok=True); fig.savefig(path); plt.close(fig)


# ------------------------------------------------------------------------------------------------ 3.5
def p35():
    d = os.path.join(P, "3_dimension_4n", "3.5_deux_echelles")
    pj = os.path.join(d, "bench35_results.json")
    if not os.path.exists(pj): return
    new = json.load(open(pj))
    old = json.load(open(os.path.join(P, "3_dimension_4n", "3.4_ambiguite", "bench34_results.json")))
    keep = {"sonde-garde-ambiguite", "sonde-depliage-vers-zero", "R1-relecture", "broyden-sans-couleur"}
    allr = [dict(r) for r in old if r["meth"] in keep] + new
    for r in allr:
        r.setdefault("seed", 0); r.setdefault("i", 0)
    meths = ["sonde-garde-ambiguite", "sonde-depliage-vers-zero", "R1-relecture", "sonde-relecture-econome", "sonde-vote",
             "sonde-deux-echelles", "sonde-deux-echelles-geo", "sonde-deux-echelles-pond", "sonde-oracle-decade",
             "sonde-oracle-decade-geo", "sonde-oracle-decade-pond", "broyden-sans-couleur"]
    cases = [(f"n = {n}, σ = {s}".replace(".", ","), (n, s)) for n in (16, 32) for s in (0.0, 0.02)]
    sel = lambda c, m: B.resume([r for r in allr if (r["n"], r["sigma"]) == c[1]], m)
    rows = [[nom(m)] + [cell(sel(c, m)) for c in cases] for m in meths]
    txt = "# Partie 3.5 : tableaux\n\nSuccès [IC 95 % de Wilson] / lectures médianes (succès seulement), 100 départs par case.\n\n"
    txt += md_table(["méthode"] + [c[0] for c in cases], rows)
    # diagnostic du test de descente
    txt += "\n## Test de descente : refus et refus à tort\n\nPart des lectures de pas refusées pour ambiguïté / par le test de descente ; entre parenthèses, part des refus du test de descente alors que la vraie moyenne des ‖ρᵢ‖² passait.\n\n"
    rows = []
    for m in meths:
        rr0 = [r for r in new if r["meth"] == m]
        if not rr0: continue
        line = [nom(m)]
        for c in cases:
            rr = [r for r in rr0 if (r["n"], r["sigma"]) == c[1]]
            if not rr: line.append("—"); continue
            pas = sum(r["pas"] for r in rr) or 1; rd = sum(r["rej_desc"] for r in rr)
            line.append(f"{fr(sum(r['rej_amb'] for r in rr)/pas)} / {fr(rd/pas)} ({fr(sum(r['rej_tort'] for r in rr)/max(rd,1))})")
        rows.append(line)
    txt += md_table(["méthode"] + [c[0] for c in cases], rows)
    # gains appariés
    txt += "\n## Gains appariés (+ : la première est moins chère ; IC 95 % bootstrap)\n\n"
    rows = []
    for a, b in (("sonde-deux-echelles", "R1-relecture"), ("sonde-deux-echelles", "sonde-depliage-vers-zero"),
                 ("sonde-deux-echelles-geo", "sonde-deux-echelles"), ("sonde-deux-echelles-pond", "sonde-deux-echelles"),
                 ("sonde-deux-echelles-geo", "broyden-sans-couleur"), ("sonde-relecture-econome", "R1-relecture"),
                 ("sonde-oracle-decade-geo", "sonde-oracle-decade")):
        line = [f"{nom(a)} contre {nom(b)}"]
        for c in cases:
            g = B.gain_apparie([r for r in allr if (r["n"], r["sigma"]) == c[1]], a, b)
            line.append("—" if g is None else f"{g['gain']:+.0f} [{g['lo']:+.0f} ; {g['hi']:+.0f}] ({g['paires']})")
        rows.append(line)
    txt += md_table(["comparaison"] + [c[0] for c in cases], rows)
    open(os.path.join(d, "bench35_tables.md"), "w").write(txt)
    fm = ["sonde-garde-ambiguite", "sonde-depliage-vers-zero", "R1-relecture", "sonde-deux-echelles", "sonde-deux-echelles-geo",
          "sonde-oracle-decade", "broyden-sans-couleur"]
    figure_cases(os.path.join(d, "figures", "01_deux_echelles.png"), cases, fm, sel,
                 "Lever l'ambiguïté à la source : capteur à deux échelles, à 24 bits par composante")


# ------------------------------------------------------------------------------------------------ 4.3
def p43():
    d = os.path.join(P, "4_autres_structures", "4.3_structures_etendues")
    pj = os.path.join(d, "bench_structures_results.json")
    if not os.path.exists(pj): return
    J = json.load(open(pj)); res, df = J["results"], J["defauts"]
    ext = {"bilat": "sonde-blocconf", "hrang": "sonde-Hrang", "qnn": "sonde-grise"}
    lab = {"bilat": "bilatère", "hrang": "ℍ-linéaire + rang 2", "qnn": "inversion réseau (tanh)"}
    txt = "# Partie 4.3 : tableaux\n\nSuccès [IC 95 %] / lectures médianes, 3 instances × 20 départs, budget 600.\n\n"
    for sg in (0.0, 0.02):
        txt += f"\n## σ = {sg}\n\n".replace(".", ",")
        rows = []
        for prob in ext:
            for n in (2, 3, 4):
                rr = [r for r in res if r["prob"] == prob and r["n"] == n and r["sigma"] == sg]
                if not rr: continue
                rows.append([lab[prob], n, fr(np.mean(df[f'{prob}/{n}']))] +
                            [cell(B.resume(rr, m)) for m in (ext[prob], "sonde-H", "sonde-full", "broyden", "broyden-deg")])
        txt += md_table(["famille", "n", "défaut ℍ-lin.", "modèle étendu", "sonde-H", "sonde-full", "broyden", "broyden-deg"], rows)
    txt += "\n## Gains appariés du modèle étendu (σ = 0)\n\n"
    rows = []
    for prob in ext:
        for b in ("sonde-H", "sonde-full", "broyden"):
            line = [lab[prob], f"contre {nom(b)}"]
            for n in (2, 3, 4):
                g = B.gain_apparie([r for r in res if r["prob"] == prob and r["n"] == n and r["sigma"] == 0.0], ext[prob], b)
                line.append("—" if g is None else f"{g['gain']:+.0f} [{g['lo']:+.0f} ; {g['hi']:+.0f}] ({g['paires']})")
            rows.append(line)
    txt += md_table(["famille", "comparaison", "n = 2", "n = 3", "n = 4"], rows)
    open(os.path.join(d, "bench_structures_tables.md"), "w").write(txt)
    cases = [(f"{lab[p]}\nn = {n}", (p, n)) for p in ext for n in (2, 3, 4)]
    sel = lambda c, m: B.resume([r for r in res if (r["prob"], r["n"]) == c[1] and r["sigma"] == 0.0 and
                                 (r["meth"] == m or (m == "étendu" and r["meth"] == ext[c[1][0]]))], ext[c[1][0]] if m == "étendu" else m)
    figure_cases(os.path.join(d, "figures", "01_structures_etendues.png"), cases, ["étendu", "sonde-H", "sonde-full", "broyden"], sel,
                 "Modèle étendu (conforme par blocs, + rang 2, boîte grise) contre ℍ-linéaire, sans structure et Broyden, σ = 0")


# ------------------------------------------------------------------------------------------------ 6
def p62():
    d = os.path.join(P, "6_reseaux"); pj = os.path.join(d, "bench_inversion_results.json")
    if not os.path.exists(pj): return
    J = json.load(open(pj)); res, df = J["results"], J["defauts"]
    lab = {"split": "tanh par composante", "radial": "tanh radial", "global": "saturation globale"}
    meths = ["sonde-H", "sonde-full", "broyden", "broyden-deg"]
    txt = "# Partie 6.2 : inversion d'un réseau, tableaux\n\nSuccès [IC 95 %] / lectures médianes, 3 instances × 20 départs.\n\n"
    for sg in (0.0, 0.02):
        txt += f"\n## σ = {sg}\n\n".replace(".", ",")
        rows = []
        for k in lab:
            for n in (2, 3, 4, 8):
                rr = [r for r in res if r["kind"] == k and r["n"] == n and r["sigma"] == sg]
                if rr: rows.append([lab[k], n, fr(np.mean(df[f'{k}/{n}']))] + [cell(B.resume(rr, m)) for m in meths])
        txt += md_table(["activation", "n", "défaut ℍ-lin.", *meths], rows)
    txt += "\n## Gain apparié de la sonde-couleur-mémoire-quaternionique contre Broyden-sans-couleur (σ = 0)\n\n"
    rows = []
    for k in lab:
        line = [lab[k]]
        for n in (2, 3, 4, 8):
            g = B.gain_apparie([r for r in res if r["kind"] == k and r["n"] == n and r["sigma"] == 0.0], "sonde-H", "broyden")
            line.append("—" if g is None else f"{g['gain']:+.0f} [{g['lo']:+.0f} ; {g['hi']:+.0f}] ({g['paires']})")
        rows.append(line)
    txt += md_table(["activation", "n = 2", "n = 3", "n = 4", "n = 8"], rows)
    open(os.path.join(d, "bench_inversion_tables.md"), "w").write(txt)
    cases = [(f"{lab[k]}\nn = {n}", (k, n)) for k in lab for n in (2, 4, 8)]
    sel = lambda c, m: B.resume([r for r in res if (r["kind"], r["n"]) == c[1] and r["sigma"] == 0.0], m)
    figure_cases(os.path.join(d, "figures", "02_inversion_activations.png"), cases, meths, sel,
                 "Inverser un réseau quaternionique : l'activation décide, σ = 0")


def p63():
    d = os.path.join(P, "6_reseaux"); pj = os.path.join(d, "bench_entrainement_results.json")
    if not os.path.exists(pj): return
    J = json.load(open(pj)); res, df = J["results"], J["defauts"]
    lab = {"lin": "linéaire", "radial-faible": "tanh radial, petit signal", "radial": "tanh radial", "split": "tanh par composante"}
    meths = ["sonde-D", "sonde-H", "sonde-full", "broyden", "broyden-deg"]
    txt = "# Partie 6.3 : entraîner un neurone, tableaux\n\nSuccès [IC 95 %] / lectures médianes, 3 instances × 20 départs.\n\n"
    for sg in (0.0, 0.02):
        txt += f"\n## σ = {sg}\n\n".replace(".", ",")
        rows = []
        for k in lab:
            for n in (2, 4, 8, 16):
                rr = [r for r in res if r["kind"] == k and r["n"] == n and r["sigma"] == sg]
                if rr:
                    dd = df[f"{k}/{n}"]
                    rows.append([lab[k], n, fr(np.mean(dd['R'])), fr(np.mean(dd['L']))] + [cell(B.resume(rr, m)) for m in meths])
        txt += md_table(["activation", "n", "défaut à droite", "défaut à gauche", *meths], rows)
    txt += "\n## Gains appariés de la sonde à droite (σ = 0)\n\n"
    rows = []
    for k in lab:
        for b in ("broyden", "sonde-full"):
            line = [lab[k], f"contre {nom(b)}"]
            for n in (2, 4, 8, 16):
                g = B.gain_apparie([r for r in res if r["kind"] == k and r["n"] == n and r["sigma"] == 0.0], "sonde-D", b)
                line.append("—" if g is None else f"{g['gain']:+.0f} [{g['lo']:+.0f} ; {g['hi']:+.0f}] ({g['paires']})")
            rows.append(line)
    txt += md_table(["activation", "comparaison", "n = 2", "n = 4", "n = 8", "n = 16"], rows)
    open(os.path.join(d, "bench_entrainement_tables.md"), "w").write(txt)
    cases = [(f"{lab[k]}\nn = {n}", (k, n)) for k in lab for n in ((4, 16) if k in ("lin", "radial-faible") else (4, 8))]
    sel = lambda c, m: B.resume([r for r in res if (r["kind"], r["n"]) == c[1] and r["sigma"] == 0.0], m)
    figure_cases(os.path.join(d, "figures", "03_entrainement_droite.png"), cases, meths, sel,
                 "Entraîner un neurone quaternionique : le quotient à droite garde la structure, σ = 0")


def p61():
    d = os.path.join(P, "6_reseaux"); pj = os.path.join(d, "bench_points_fixes_results.json")
    if not os.path.exists(pj): return
    J = json.load(open(pj))
    txt = "# Partie 6.1 : points fixes de RNN quaternioniques, tableaux\n"
    if "A" in J:
        A = J["A"]
        meths = ["sonde-H", "sonde-Hid", "sonde-H-relecture", "sonde-full", "broyden", "broyden-id"]
        for sg in (0.0, 0.02):
            txt += f"\n## A. Pièce apprise, σ = {sg}\n\nSuccès / lectures médianes (3 réseaux × 10 départs ; erreur de mémoire = écart relatif entre le point fixe atteint et la cible apprise).\n\n".replace("σ = 0.", "σ = 0,")
            rows = []
            for n in (4, 8):
                for kind in ("split", "radial", "global"):
                    for piece in ("aucune", "H", "penalise-1", "penalise-0.1", "libre"):
                        rr = [r for r in A if r["n"] == n and r["kind"] == kind and r["piece"] == piece and r["sigma"] == sg]
                        if not rr: continue
                        rows.append([n, kind, piece, fr(rr[0]['defaut_P']), fr(rr[0]['defaut_J']), fr(rr[0]['err_memoire'], 3)] +
                                    [B.cellule(B.resume(rr, m)) for m in meths])
            txt += md_table(["n", "activation", "pièce", "défaut de P", "défaut de J", "erreur mémoire", *[nom(m) for m in meths]], rows)
    if "B" in J:
        txt += "\n## B. Plusieurs points fixes (40 départs uniformes dans [−1,2 ; 1,2]^(4n), budget 600)\n\n"
        rows = []
        for b in J["B"]:
            c = b["certificats"]
            cs = f"{sum(x['deg_lu'] == x['deg_vrai'] for x in c)}/{len(c)}" if c else "—"
            rows.append([b["n"], b["seed"], b["n_vrais"], fr(b['defaut_moyen']), nom(b["meth"]), f"{100*b['succes']:.0f} %",
                         b["distincts"], b["lectures_total"], f"{b['lectures_med']:.0f}", cs])
        txt += md_table(["n", "réseau", "points fixes (Newton exact, 400 départs)", "défaut moyen", "méthode", "succès", "points distincts trouvés",
                         "lectures totales", "lectures médianes", "degré lu = signe de det J"], rows)
    open(os.path.join(d, "bench_points_fixes_tables.md"), "w").write(txt)
    if "A" in J:
        A = J["A"]
        cases = [(f"{k}\n{p}", (k, p)) for k in ("split", "radial", "global") for p in ("aucune", "H", "libre")]
        sel = lambda c, m: B.resume([r for r in A if (r["kind"], r["piece"]) == c[1] and r["n"] == 8 and r["sigma"] == 0.0], m)
        figure_cases(os.path.join(d, "figures", "01_points_fixes.png"), cases, ["sonde-H", "sonde-Hid", "broyden", "broyden-id"], sel,
                     "Point fixe d'un RNN quaternionique à pièce apprise, n = 8, σ = 0")


def p64():
    d = os.path.join(P, "6_reseaux"); pj = os.path.join(d, "bench_couche_results.json")
    if not os.path.exists(pj): return
    J = json.load(open(pj)); res, df = J["results"], J["defauts"]
    lab = {"lin": "linéaire", "radial-faible": "tanh radial, petit signal", "radial": "tanh radial"}
    meths = ["sonde-Dbloc", "sonde-D", "sonde-full", "broyden", "broyden-blocs", "broyden-deg"]
    cas = sorted({(r["m"], r["p"]) for r in res}, key=lambda x: (x[1], x[0]))
    txt = "# Partie 6.4 : entraîner une couche entière, tableaux\n\nSuccès [IC 95 %] / lectures médianes, 3 instances × 10 départs.\n\n"
    for sg in (0.0, 0.02):
        txt += f"\n## σ = {sg}\n\n".replace(".", ",")
        rows = []
        for k in lab:
            for m, p in cas:
                rr = [r for r in res if r["kind"] == k and (r["m"], r["p"]) == (m, p) and r["sigma"] == sg]
                if rr: rows.append([lab[k], m, p, 4 * m * p, fr(np.mean(df[f"{k}/{m}/{p}"]["R"]))] + [cell(B.resume(rr, x)) for x in meths])
        txt += md_table(["activation", "m (sorties)", "p (poids par sortie)", "dimension", "défaut à droite", *[nom(x) for x in meths]], rows)
    txt += "\n## Gains appariés de la sonde par sortie (σ = 0)\n\n"
    rows = []
    for k in lab:
        for b in ("broyden-blocs", "broyden", "sonde-D"):
            line = [lab[k], f"contre {nom(b)}"]
            for m, p in cas:
                g = B.gain_apparie([dict(r, seed=r["seed"]) for r in res if r["kind"] == k and (r["m"], r["p"]) == (m, p) and r["sigma"] == 0.0], "sonde-Dbloc", b)
                line.append("—" if g is None else f"{g['gain']:+.0f} [{g['lo']:+.0f} ; {g['hi']:+.0f}]")
            rows.append(line)
    txt += md_table(["activation", "comparaison", *[f"m = {m}, p = {p}" for m, p in cas]], rows)
    open(os.path.join(d, "bench_couche_tables.md"), "w").write(txt)
    # lectures médianes selon m, à p fixé
    fig, axs = plt.subplots(1, 2, figsize=(12, 4.4))
    for ax, p in zip(axs, (4, 8)):
        ms = sorted({r["m"] for r in res if r["p"] == p})
        for k2, x in enumerate(["sonde-Dbloc", "sonde-D", "broyden-blocs", "broyden"]):
            ys = []
            for m in ms:
                s_ = B.resume([r for r in res if r["kind"] == "radial-faible" and (r["m"], r["p"]) == (m, p) and r["sigma"] == 0.0], x)
                ys.append(np.nan if (s_ is None or s_["cost"] is None) else s_["cost"])
            ax.plot(ms, ys, "o-", color=COL[k2], lw=2, label=nom(x))
        ax.set_xscale("log", base=2); ax.set_xticks(ms); ax.set_xticklabels(ms)
        ax.set_xlabel("m (nombre de sorties)"); ax.set_ylabel("lectures médianes (succès)")
        ax.set_title(f"p = {p} poids par sortie, tanh radial petit signal, σ = 0")
    axs[0].legend(fontsize=8)
    fig.tight_layout(); fig.savefig(os.path.join(d, "figures", "04_couche_lectures_vs_m.png")); plt.close(fig)


def _best(rr, base):
    """Meilleur des pas de différence finie 0,5 et 1 pour la famille `base` (« bq-deg » ou « bq-deg2 ») : succès puis lectures."""
    best = None
    for h in ("0.5", "1"):
        s_ = B.resume(rr, f"{base}-{h}")
        if s_ and (best is None or (s_["rate"], -(s_["cost"] or 1e9)) > (best[1]["rate"], -(best[1]["cost"] or 1e9))): best = (h, s_)
    return [dict(r, meth=base) for r in rr if best and r["meth"] == f"{base}-{best[0]}"]


def _ligne(ax, xs, series, cols, ylab, logy=False):
    for (lab, ys), c in zip(series.items(), cols):
        ax.plot(xs, ys, "o-", color=c, lw=2, label=lab)
    ax.set_ylabel(ylab)
    if logy: ax.set_yscale("log")


def p56():
    """Partie 5.6 : boîte noire ℍ-linéaire instable (conditionnement) et explosive (plage dynamique), à bits égaux."""
    d = os.path.join(P, "5_armes_egales"); pj = os.path.join(d, "bench_boite_instable_results.json")
    if not os.path.exists(pj): return
    res = json.load(open(pj))
    meths = ["sonde", "sonde-oracle", "sonde-2e-geo", "bq-deg2", "bq-deg", "bq", "broyden"]
    cas = sorted({(r["k"], r["p"]) for r in res}, key=lambda c: (c[1] > 0, c[0], c[1]))
    def sel(n, k, p, sg):
        rr = [r for r in res if (r["n"], r["k"], r["p"], r["sigma"]) == (n, k, p, sg)]
        return rr + _best(rr, "bq-deg") + _best(rr, "bq-deg2")
    txt = ("# Partie 5.6 : boîte noire ℍ-linéaire instable, tableaux\n\nSuccès [IC 95 %] / lectures médianes ; « cond » : conditionnement 10^k ; "
           "« explosif » : exposant p (F croît comme ‖d‖^{2p+1}). Fenêtre de la mesure absolue à bits égaux : 10⁻⁴ à 10⁶. "
           "Broyden-quaternionique-mêmes-mesures et à bits égaux : meilleur des deux pas de différence finie (0,5 et 1) par case.\n")
    for sg in (0.0, 0.02):
        for n in sorted({r["n"] for r in res}):
            rows = []
            for k, p in cas:
                rr = sel(n, k, p, sg)
                if rr: rows.append([f"cond k = {k}" if p == 0 and k else ("référence (k = 0, p = 0)" if p == 0 else f"explosif p = {p}")] + [cell(B.resume(rr, m)) for m in meths])
            if rows:
                txt += f"\n## n = {n} (dimension {4 * n}), σ = {fr(sg)}\n\n" + md_table(["problème", *[nom(m) for m in meths]], rows)
    open(os.path.join(d, "bench_boite_instable_tables.md"), "w").write(txt)
    n0 = 8 if 8 in {r["n"] for r in res} else sorted({r["n"] for r in res})[0]
    fig, axs = plt.subplots(2, 2, figsize=(12.5, 7.6), sharex="col")
    fam = [("cond", [c for c in cas if c[1] == 0], lambda c: c[0], "conditionnement 10^k (k)"),
           ("explosif", [c for c in cas if c[0] == 0], lambda c: c[1], "exposant p de l'explosion")]
    noms = {"sonde": nom("sonde-Dbloc").replace("-à-droite-par-sortie", ""), "sonde-oracle": nom("sonde-oracle"), "sonde-2e-geo": nom("sonde-2e-geo"),
            "bq-deg2": nom("bq-deg2"), "bq-deg": nom("bq-deg"), "bq": nom("bq")}
    cols = {"sonde": PAL["orange"], "sonde-oracle": PAL["yellow"], "sonde-2e-geo": PAL["magenta"], "bq-deg2": PAL["blue"], "bq-deg": PAL["aqua"], "bq": MUTED}
    for c, (f_, cs, xf, xl) in enumerate(fam):
        xs = [xf(x) for x in cs]
        for m in noms:
            sr = [B.resume(sel(n0, k, p, 0.0), m) for k, p in cs]
            axs[0, c].plot(xs, [np.nan if x is None else 100 * x["rate"] for x in sr], "o-", color=cols[m], lw=2, label=noms[m])
            axs[1, c].plot(xs, [np.nan if (x is None or x["cost"] is None or x["rate"] < 0.3) else x["cost"] for x in sr], "o-", color=cols[m], lw=2)
        axs[1, c].set_xlabel(xl); axs[1, c].set_yscale("log")
        axs[0, c].set_title(("problème « cond » : A mal conditionné" if f_ == "cond" else "problème « explosif » : bien conditionné, F traverse des décades") + f", n = {n0}, σ = 0", fontsize=9.5)
    axs[0, 0].set_ylabel("succès (%)"); axs[1, 0].set_ylabel("lectures médianes (succès ≥ 30 %)")
    axs[0, 0].legend(fontsize=7, loc="lower left")
    fig.tight_layout(); fig.savefig(os.path.join(d, "figures", "03_boite_instable.png")); plt.close(fig)


def p57():
    """Partie 5.7 : capteur d'amplitude peu fiable (gain aléatoire par lecture)."""
    d = os.path.join(P, "5_armes_egales"); pj = os.path.join(d, "bench_gain_amplitude_results.json")
    if not os.path.exists(pj): return
    res = json.load(open(pj)); ns = sorted({r["n"] for r in res}); sls = sorted({r["sl"] for r in res})
    meths = ["sonde", "sonde-2e-geo", "bq-deg", "bq", "broyden"]
    def sel(n, sl):
        rr = [r for r in res if (r["n"], r["sl"]) == (n, sl)]
        return rr + _best(rr, "bq-deg")
    txt = ("# Partie 5.7 : capteur d'amplitude peu fiable, tableaux\n\nSuccès [IC 95 %] / lectures médianes ; σℓ : écart-type du bruit sur ln r, en népers "
           "(angle et axe sans bruit). Broyden-quaternionique-sans-couleur et Broyden-sans-couleur lisent la valeur exacte (ni bruit ni quantification) : bornes.\n")
    for n in ns:
        rows = [[fr(sl)] + [cell(B.resume(sel(n, sl), m)) for m in meths] for sl in sls if sel(n, sl)]
        txt += f"\n## n = {n} (dimension {4 * n})\n\n" + md_table(["σℓ (népers)", *[nom(m) for m in meths]], rows)
    open(os.path.join(d, "bench_gain_amplitude_tables.md"), "w").write(txt)
    fig, axs = plt.subplots(1, len(ns), figsize=(4.6 * len(ns), 4.2), sharey=True)
    axs = np.atleast_1d(axs)
    cols = {"sonde": PAL["orange"], "sonde-2e-geo": PAL["magenta"], "bq-deg": PAL["aqua"], "bq": MUTED}
    for ax, n in zip(axs, ns):
        for m in cols:
            sr = [B.resume(sel(n, sl), m) for sl in sls]
            ax.plot(range(len(sls)), [np.nan if x is None else 100 * x["rate"] for x in sr], "o-", color=cols[m], lw=2, label=nom(m))
        ax.set_xticks(range(len(sls))); ax.set_xticklabels([fr(x, 2) for x in sls], fontsize=8)
        ax.set_xlabel("bruit de gain σℓ (népers ; 1 ≈ facteur 2,7)"); ax.set_title(f"n = {n}", fontsize=10)
    axs[0].set_ylabel("succès (%)"); axs[0].legend(fontsize=7.5, loc="lower left")
    fig.suptitle("Amplitude peu fiable : le gain aléatoire de chaque lecture, angle et axe précis", x=0.02, ha="left", fontsize=11)
    fig.tight_layout(rect=(0, 0, 1, 0.93)); fig.savefig(os.path.join(d, "figures", "04_gain_amplitude.png")); plt.close(fig)


def p58():
    """Partie 5.8 : suivi d'un zéro mobile."""
    d = os.path.join(P, "5_armes_egales"); pj = os.path.join(d, "bench_suivi_results.json")
    if not os.path.exists(pj): return
    res = json.load(open(pj)); ns = sorted({r["n"] for r in res}); vs = sorted({r["v"] for r in res}); rs = sorted({r["r"] for r in res})
    meths = ["sonde", "sonde-2e", "sonde-oracle", "bq", "bq-deg", "bq-deg2", "bq-sans-suivi"]
    def med(n, v, r, sg, m):
        x = [q["retard"] for q in res if (q["n"], q["v"], q["r"], q["sigma"], q["meth"]) == (n, v, r, sg, m) and q["retard"] is not None]
        nb = sum(1 for q in res if (q["n"], q["v"], q["r"], q["sigma"], q["meth"]) == (n, v, r, sg, m))
        return (float(np.median(x)) if x else np.nan), len(x), nb
    txt = ("# Partie 5.8 : suivi d'un zéro mobile, tableaux\n\nRetard médian de suivi (en pas de temps) = distance médiane du meilleur point lu au zéro du moment, "
           "divisée par la vitesse v, sur les 40 derniers pas ; entre parenthèses : nombre de runs dont le suivi est défini / lancés (un run sans lecture exploitable "
           "n'a pas de retard). Un suivi parfait au rythme des lectures donnerait ≈ 1,7 / r.\n")
    for sg in sorted({r["sigma"] for r in res}):
        for n in ns:
            rows = []
            for v in vs:
                for r in rs:
                    cell_ = []
                    for m in meths:
                        x, k, nb = med(n, v, r, sg, m); cell_.append("—" if np.isnan(x) else f"{fr(x)} ({k}/{nb})")
                    rows.append([fr(v), r, *cell_])
            txt += f"\n## n = {n}, σ = {fr(sg)}\n\n" + md_table(["vitesse v", "lectures par pas r", *[nom(m) for m in meths]], rows)
    open(os.path.join(d, "bench_suivi_tables.md"), "w").write(txt)
    fig, axs = plt.subplots(1, len(vs), figsize=(4.3 * len(vs), 4.2), sharey=True); axs = np.atleast_1d(axs)
    cols = {"sonde": PAL["orange"], "sonde-2e": PAL["magenta"], "sonde-oracle": PAL["yellow"], "bq": MUTED, "bq-deg": PAL["aqua"], "bq-deg2": PAL["blue"], "bq-sans-suivi": PAL["red"]}
    for ax, v in zip(axs, vs):
        for m in cols:
            ax.plot(rs, [med(ns[-1], v, r, 0.0, m)[0] for r in rs], "o-", color=cols[m], lw=2, label=nom(m))
        ax.set_xscale("log", base=2); ax.set_yscale("log"); ax.set_xticks(rs); ax.set_xticklabels(rs)
        ax.set_xlabel("lectures par pas de temps r"); ax.set_title(f"vitesse v = {fr(v)}, n = {ns[-1]}, σ = 0", fontsize=9.5)
    axs[0].set_ylabel("retard de suivi (pas de temps)"); axs[0].legend(fontsize=7, loc="upper right")
    fig.tight_layout(); fig.savefig(os.path.join(d, "figures", "05_suivi.png")); plt.close(fig)


def p65():
    """Partie 6.5 : ablations de la sonde en régime saturé."""
    d = os.path.join(P, "6_reseaux"); pj = os.path.join(d, "bench_diag_sature_results.json")
    if not os.path.exists(pj): return
    res = json.load(open(pj)); noms = [m for m in ["base", "pond", "geo", "pas-1", "acc-1", "pond+pas-1", "pond+pas-1+acc-1", "rcond-1e-10", "eta-fixe", "pond+pas-1+rcond+eta-fixe", "aucun-test", "aucun-test-pas-065", "bq-sans-rl", "bq"] if any(r["meth"] == m for r in res)]
    lab = {"base": "sonde-oracle-décade (réglage de la partie)", "pond": "+ test pondéré", "geo": "+ test géométrique", "pas-1": "+ pas de Newton entier (1,0)",
           "acc-1": "+ seuil d'acceptation 1,0", "pond+pas-1": "+ test pondéré + pas 1,0", "pond+pas-1+acc-1": "+ test pondéré + pas 1,0 + seuil 1,0",
           "bq": "Broyden-quaternionique-sans-couleur",
           "aucun-test": "aucun test de descente (pas 1,0)", "aucun-test-pas-065": "aucun test de descente (pas 0,65)", "bq-sans-rl": "Broyden-quaternionique sans recherche linéaire (aucun test)",
           "rcond-1e-10": "+ coupure des valeurs singulières à 10⁻¹⁰ (écartée)", "eta-fixe": "+ taille des sondes fixe après échec (écartée)",
           "pond+pas-1+rcond+eta-fixe": "+ test pondéré + pas 1,0 + les deux précédentes"}
    txt = "# Partie 6.5 : pourquoi la sonde s'effondre en régime saturé, tableaux\n\nSuccès [IC 95 %] / lectures médianes ; décade donnée (oracle), σ = 0, 3 instances × 20 départs.\n"
    for kind in ("radial", "radial-faible"):
        ns = sorted({r["n"] for r in res if r["kind"] == kind})
        if not ns: continue
        rows = [[lab[m]] + [cell(B.resume([r for r in res if r["kind"] == kind and r["n"] == n], m)) for n in ns] for m in noms]
        txt += f"\n## neurone {'tanh radial (défaut 0,27)' if kind == 'radial' else 'tanh radial petit signal (défaut 0,05)'}\n\n" + md_table(["variante", *[f"n = {n}" for n in ns]], rows)
    open(os.path.join(d, "bench_diag_sature_tables.md"), "w").write(txt)
    ns = sorted({r["n"] for r in res if r["kind"] == "radial"})
    if not ns: return
    fig, ax = plt.subplots(figsize=(9, 4.2)); w = 0.8 / len(ns); x = np.arange(len(noms))
    for c, n in enumerate(ns):
        S_ = [B.resume([r for r in res if r["kind"] == "radial" and r["n"] == n], m) for m in noms]
        ax.bar(x + c * w, [0 if s_ is None else 100 * s_["rate"] for s_ in S_], w, color=[PAL["blue"], PAL["orange"]][c % 2], label=f"n = {n}")
    ax.set_xticks(x + w * (len(ns) - 1) / 2); ax.set_xticklabels([lab[m].replace("sonde-oracle-décade ", "base\n") for m in noms], rotation=28, ha="right", fontsize=7.5)
    ax.set_ylabel("succès (%)"); ax.set_title("Neurone tanh radial : on retire un défaut de la sonde à la fois (décade donnée)", fontsize=10); ax.legend()
    fig.tight_layout(); fig.savefig(os.path.join(d, "figures", "05_diag_sature.png")); plt.close(fig)


ENT55 = {"sonde": "sonde de la partie d'origine (quaternionique / à droite / par sortie)", "sonde-oracle": "même sonde, oracle-décade"}


def p55():
    """Partie 5.5 : même structure pour tous. Fusionne le banc 5.5 avec les résultats des parties 3.4, 3.5, 6.3 et 6.4."""
    d = os.path.join(P, "5_armes_egales"); pj = os.path.join(d, "bench_meme_structure_results.json")
    if not os.path.exists(pj): return
    res = [dict(r) for r in json.load(open(pj))]
    d3 = os.path.join(P, "3_dimension_4n")
    ren34 = {"sonde-garde-ambiguite": "sonde", "broyden-sans-couleur": "broyden"}
    for r in json.load(open(os.path.join(d3, "3.4_ambiguite", "bench34_results.json"))):
        if r["meth"] in ren34: res.append(dict(fam="H", kind="", n=r["n"], m=1, seed=r["seed"], i=r["i"], meth=ren34[r["meth"]], sigma=r["sigma"], ok=r["ok"], cost=r["cost"]))
    for r in json.load(open(os.path.join(d3, "3.5_deux_echelles", "bench35_results.json"))):
        if r["meth"] == "sonde-deux-echelles-geo": res.append(dict(fam="H", kind="", n=r["n"], m=1, seed=r["seed"], i=r["i"], meth="sonde-2e-geo", sigma=r["sigma"], ok=r["ok"], cost=r["cost"]))
    d6 = os.path.join(P, "6_reseaux")
    for r in json.load(open(os.path.join(d6, "bench_entrainement_results.json")))["results"]:
        x = {"sonde-D": "sonde", "broyden": "broyden"}.get(r["meth"])
        if x and r["kind"] in ("lin", "radial-faible", "radial"):
            res.append(dict(fam="neurone", kind=r["kind"], n=r["n"], m=1, seed=r["seed"], i=r["i"], meth=x, sigma=r["sigma"], ok=r["ok"], cost=r["cost"]))
    for r in json.load(open(os.path.join(d6, "bench_couche_results.json")))["results"]:
        x = {"sonde-Dbloc": "sonde", "broyden": "broyden", "broyden-blocs": "broyden-blocs"}.get(r["meth"])
        if x: res.append(dict(fam="couche", kind=r["kind"], n=r["n"], m=r["m"], seed=r["seed"], i=r["i"], meth=x, sigma=r["sigma"], ok=r["ok"], cost=r["cost"]))
    new = {(r["fam"], r["kind"], r["n"], r["m"]) for r in json.load(open(pj))}
    def sel(fam, kind, n, m, sg):
        rr = [r for r in res if (r["fam"], r["kind"], r["n"], r["m"], r["sigma"]) == (fam, kind, n, m, sg)]
        # Broyden-quaternionique-mêmes-mesures : le meilleur des deux pas (succès, puis lectures médianes) ; réglage favorable
        best = None
        for h in ("0.5", "1"):
            s_ = B.resume(rr, "bq-deg-" + h)
            if s_ and (best is None or (s_["rate"], -(s_["cost"] or 1e9)) > (best[1]["rate"], -(best[1]["cost"] or 1e9))): best = (h, s_)
        if best: rr = rr + [dict(r, meth="bq-deg") for r in rr if r["meth"] == "bq-deg-" + best[0]]
        return rr
    lab = {"": "ℍ-linéaire (3.2)", "lin": "identité", "radial-faible": "tanh radial, petit signal", "radial": "tanh radial"}
    fams = [("H", "Problème ℍ-linéaire des parties 3.1–3.5 (structure à gauche), 10 instances × 10 départs, départs à 1,5·√n"),
            ("neurone", "Entraîner un neurone (partie 6.3, structure à droite), 3 instances × 20 départs"),
            ("couche", "Entraîner une couche (partie 6.4, à droite + par sortie), 3 instances × 10 départs")]
    meths = ["sonde", "sonde-oracle", "sonde-2e-geo", "bq-deg", "bq", "broyden"]
    txt = ("# Partie 5.5 : même structure pour tous, tableaux\n\nSuccès [IC 95 %] / lectures médianes. Broyden-quaternionique-"
           "mêmes-mesures : meilleur des deux pas de différence finie (0,5 et 1) par case.\n")
    gains = []
    for fam, titre in fams:
        cas = sorted({(r["kind"], r["n"], r["m"]) for r in res if r["fam"] == fam and (fam, r["kind"], r["n"], r["m"]) in new},
                     key=lambda c: (list(lab).index(c[0]), c[1], c[2]))
        for sg in (0.0, 0.02):
            txt += f"\n## {titre}, σ = {fr(sg)}\n\n"
            rows = []
            for k, n, m in cas:
                rr = sel(fam, k, n, m, sg)
                taille = f"n = {n}" if fam != "couche" else f"{m} × {n // m}"
                rows.append([lab[k], taille, 4 * n] + [cell(B.resume(rr, x)) for x in meths])
                if sg == 0.0 or True:
                    for a, b in (("bq", "sonde"), ("bq-deg", "sonde-oracle"), ("bq", "broyden")):
                        g = B.gain_apparie(rr, a, b)
                        gains.append((fam, lab[k], taille, sg, a, b, g, B.resume(rr, a), B.resume(rr, b)))
            txt += md_table(["structure / activation", "taille", "dimension", *[ENT55.get(x, nom(x)) for x in meths]], rows)
    txt += "\n## Gains appariés (+ : la première méthode est moins chère ; IC 95 % bootstrap ; départs réussis par les deux)\n\n"
    rows = [[{"H": "ℍ-linéaire", "neurone": "neurone", "couche": "couche"}[f], k, t, fr(sg), f"{ENT55.get(a, nom(a))} contre {ENT55.get(b, nom(b))}",
             "—" if g is None else f"{g['gain']:+.0f} [{g['lo']:+.0f} ; {g['hi']:+.0f}] ({g['paires']} paires)"]
            for f, k, t, sg, a, b, g, sa, sb in gains]
    txt += md_table(["famille", "activation", "taille", "σ", "comparaison", "gain %"], rows)
    open(os.path.join(d, "bench_meme_structure_tables.md"), "w").write(txt)
    json.dump([dict(fam=f, kind=k, taille=t, sigma=sg, a=a, b=b, gain=g, ra=sa, rb=sb) for f, k, t, sg, a, b, g, sa, sb in gains],
              open(os.path.join(d, "bench_meme_structure_gains.json"), "w"))
    # figure : lectures médianes selon la dimension, σ = 0, une colonne par famille (cas presque linéaires)
    fig, axs = plt.subplots(1, 3, figsize=(15, 4.6))
    show = [("H", ""), ("neurone", "radial-faible"), ("couche", "radial-faible")]
    noms = {"sonde": "sonde-couleur-mémoire-quaternionique (la sonde de la partie)", "sonde-oracle": nom("sonde-oracle"),
            "bq-deg": nom("bq-deg"), "bq": nom("bq"), "broyden": nom("broyden")}
    for ax, (fam, k) in zip(axs, show):
        ns = sorted({r["n"] for r in res if r["fam"] == fam and r["kind"] == k and (fam, k, r["n"], r["m"]) in new})
        for c, x in enumerate(["sonde", "sonde-oracle", "bq-deg", "bq", "broyden"]):
            ys = []
            for n in ns:
                ms = sorted({r["m"] for r in res if r["fam"] == fam and r["kind"] == k and r["n"] == n and (fam, k, n, r["m"]) in new})
                s_ = B.resume(sel(fam, k, n, ms[-1], 0.0), x)
                ys.append(np.nan if (s_ is None or s_["cost"] is None or s_["rate"] < 0.5) else s_["cost"])
            ax.plot([4 * n for n in ns], ys, "o-" if x != "broyden" else "o:", color=[PAL["orange"], PAL["yellow"], PAL["blue"], PAL["aqua"], MUTED][c], lw=2, label=noms[x])
        ax.set_xscale("log", base=2); ax.set_xticks([4 * n for n in ns]); ax.set_xticklabels([4 * n for n in ns])
        ax.set_xlabel("dimension réelle"); ax.set_ylabel("lectures médianes (succès ≥ 50 %)")
        ax.set_title({"H": "ℍ-linéaire (parties 3.1–3.5)", "neurone": "un neurone, tanh petit signal (6.3)", "couche": "une couche, tanh petit signal (6.4)"}[fam], fontsize=10)
    axs[0].legend(fontsize=7.5, loc="upper left")
    fig.suptitle("À structure égale, σ = 0 : couleur (orange, jaune) contre valeur (bleu : mêmes mesures, vert : exacte)", fontsize=11)
    fig.tight_layout(); fig.savefig(os.path.join(d, "figures", "02_meme_structure.png")); plt.close(fig)


if __name__ == "__main__":
    todo = sys.argv[1:] or ["35", "43", "55", "56", "57", "58", "61", "62", "63", "64", "65"]
    for t in todo:
        globals()["p" + t]()
        print("ok", t)
