"""Schémas de synthèse : le fil de l'étude (06) et le bilan en un graphe (07). Usage : python schemas/make_syntheses.py
Les chiffres du bilan sont recopiés des tableaux des parties (voir parties/*/README.md) ; la partie 4.1 est relue dans son JSON."""
import os, sys, json
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "src"))
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
from report import PAL, style, INK, INK2, MUTED, nom
style()

# Le fil suit l'ordre chronologique de l'étude ; chaque boîte porte le numéro de la partie dans le rapport réorganisé.
STEPS = [
    ("0  capteurs", "que lit-on avec B bits ?", "#eef3fb"),
    ("1.1  tranche ℂ ⊂ ℍ", "la mécanique complexe\npasse-t-elle ?", "#e6f6ef"),
    ("1.2  ℍ → ℍ général", "4 sondes : Broyden gagne", "#fdeee6"),
    ("1.3  Möbius", "7 paramètres au lieu de 16", "#fdeee6"),
    ("2  minimums", "la luminosité s'aplatit", "#f1f1ee"),
    ("1.5  certificat", "ce que l'axe apporte en plus", "#e6f6ef"),
    ("5  armes égales", "mêmes mesures pour tous", "#fdeee6"),
    ("1.4  améliorations", "le bruit aidait la sonde ?", "#f1f1ee"),
    ("3.1–3.2  ℍⁿ", "structure ℍ-linéaire\nen dimension 4n", "#e6f6ef"),
    ("4  autres structures", "le gain survit-il ?\nle défaut se contrôle-t-il ?", "#eef3fb"),
    ("3.3  décrochage", "pourquoi le gain\nplafonne vers n = 8 ?", "#fdeee6"),
    ("3.4–3.5  ambiguïté", "relire, ou lire\nla décade ?", "#e6f6ef"),
    ("4.3  étendre", "plus de paramètres\nquand ça manque", "#fdeee6"),
    ("6  réseaux", "points fixes, inversion,\nneurone et couche (à droite)", "#e6f6ef"),
    ("5.5  même structure", "Broyden quaternionique :\nla valeur est moins chère", "#fdeee6"),
    ("5.6–5.8  niches", "instable, amplitude,\nzéro mobile : rien", "#fdeee6"),
    ("6.5  saturation", "le test de descente\nexplique l'essentiel", "#eef3fb"),
]
LIMITS = ["ça se lit : essayons", "ça marche, mais l'axe est inutile ici", "il faut une structure", "gain étroit : et les minima ?",
          "chercher autre chose que la vitesse", "et si tout le monde lisait pareil ?", "un indice : le bruit aide", "en dim 4, rien à économiser",
          "le gain est-il construit ?", "et en grande dimension ?", "lever l'ambiguïté", "et si la structure manque ?",
          "et sur de vrais réseaux ?", "la structure demande-t-elle la couleur ?", "une niche pour la couleur ?",
          "pourquoi la sonde s'effondre ?"]
PER_ROW = 9


def fil():
    L = 0.3 + 2.15 * PER_ROW
    fig, ax = plt.subplots(figsize=(L, 5.0)); ax.set_xlim(0, L); ax.set_ylim(0.8, 6.3); ax.set_axis_off()
    nb = len(STEPS); W_BOX = 1.85
    pos = [(0.15 + 2.15 * (k % PER_ROW), 4.4 if k < PER_ROW else 1.4) for k in range(nb)]
    pos = pos[:PER_ROW] + list(reversed(pos[PER_ROW:]))
    if nb - PER_ROW < PER_ROW:                                   # seconde ligne incomplète : alignée à droite, sous la fin
        pos = pos[:PER_ROW] + [(0.15 + 2.15 * (PER_ROW - 1 - j), 1.4) for j in range(nb - PER_ROW)]
    for (x, y), (t, s, c) in zip(pos, STEPS):
        ax.add_patch(FancyBboxPatch((x, y), W_BOX, 1.15, boxstyle="round,pad=0.02,rounding_size=0.06", fc=c, ec=INK2, lw=1))
        ax.text(x + W_BOX / 2, y + 0.8, t, ha="center", fontsize=9.2, color=INK, weight="bold")
        ax.text(x + W_BOX / 2, y + 0.4, s, ha="center", va="center", fontsize=7.2, color=INK2)
    for k in range(nb - 1):
        (x0, y0), (x1, y1) = pos[k], pos[k + 1]
        if y0 == y1:
            a, b = ((x0 + W_BOX, y0 + 0.57), (x1, y1 + 0.57)) if x1 > x0 else ((x0, y0 + 0.57), (x1 + W_BOX, y1 + 0.57))
            ty = y0 + 1.35 if y0 > 3 else y0 - 0.32
        else:
            a, b = (x0 + W_BOX / 2, y0), (x1 + W_BOX / 2, y1 + 1.15); ty = 3.0
        ax.add_patch(FancyArrowPatch(a, b, arrowstyle="-|>", mutation_scale=13, color=MUTED, lw=1.2))
        tx = (a[0] + b[0]) / 2 + (0.15 if y0 != y1 else 0)
        ax.text(tx, ty, LIMITS[k], ha="center" if y0 == y1 else "left", fontsize=6.9, color=PAL["orange"], style="italic")
    ax.text(0.2, 6.0, "Le fil de l'étude (ordre chronologique) : chaque étape naît d'une limite mesurée de la précédente", fontsize=12, color=INK)
    fig.savefig(os.path.join(HERE, "06_fil_de_l_etude.png")); plt.close(fig)


# (situation, meilleure sonde couleur (succès, lectures), meilleure référence classique (succès, lectures), ★ = mêmes mesures, partie)
ROWS = [
    ("tranche ℂ ⊂ ℍ", "sonde-couleur-sans-axe-signée", (100, 17), "Broyden-sans-couleur", (92, 18), False, "1.1"),
    ("ℍ → ℍ quelconque, zéro isolé", "sonde-couleur-mémoire 48 bits", (100, 15), "Broyden-sans-couleur", (100, 8), False, "1.2"),
    ("Möbius (défaut ≤ 0,1), 24 bits", "sonde-couleur-mémoire-conforme", (100, 17), "Broyden-sans-couleur", (100, 13), False, "1.3"),
    ("minimum non nul, facile*", "sonde-couleur-à-cible 36 bits", (100, 11), "Levenberg-Marquardt-sans-couleur", (100, 17), False, "2"),
    ("minimum non nul, vallée étroite*", "sonde-couleur-à-cible 36 bits", (100, 30), "Levenberg-Marquardt-sans-couleur", (100, 17), False, "2"),
    ("ℍ⁴ ℍ-linéaire (dim 16)", "sonde-couleur-mémoire-quaternionique", (97, 23), "Broyden-sans-couleur", (100, 29), False, "3.1"),
    ("ℍ⁴ dégradé ε = 1", "sonde-couleur-mémoire", (65, 86), "Broyden-sans-couleur", (100, 44), False, "3.1"),
    ("ℍ⁸ ℍ-linéaire (dim 32)", "sonde-couleur-mémoire-quaternionique", (90, 47), "Broyden-sans-couleur", (100, 49), False, "3.2"),
    ("ℍ³² ℍ-linéaire (dim 128)", "sonde-couleur-mémoire-quaternionique", (63, 167), "Broyden-sans-couleur", (100, 151), False, "3.2"),
    ("ℍ³² ℍ-linéaire, dépliage vers zéro", "sonde-dépliage-vers-zéro", (100, 125.5), "Broyden-sans-couleur", (100, 151), False, "3.3"),
    ("ℍ³² bruit 0,02, relecture", "sonde-relecture", (59, 219), "Broyden-sans-couleur", (100, 152), False, "3.4"),
    ("ℍ³² ℍ-linéaire, deux échelles + test géom. (24 bits)", "sonde-deux-échelles + test géométrique", (99, 122), "Broyden-sans-couleur", (100, 152), False, "3.5"),
    ("réservoir ℍ⁸, tanh radial (défaut 0,14)", "sonde-couleur-mémoire-quaternionique-départ-identité", (100, 18), "Broyden-sans-couleur-départ-identité", (100, 12), False, "4.2"),
    ("réservoir ℍ³², saturation globale (défaut 0,03)", "sonde-couleur-mémoire-quaternionique-départ-identité", (100, 27), "Broyden-sans-couleur-départ-identité", (100, 15), False, "4.2"),
    ("mêmes mesures, 24 bits, bruit 0,02", "sonde-couleur-mémoire", (100, 16), "Broyden-mêmes-mesures", (99, 26), True, "5"),
    ("mêmes mesures, 36 bits, sans bruit", "sonde-couleur-mémoire", (93, 13), "Broyden-mêmes-mesures", (100, 12), True, "5"),
    ("mêmes mesures, ℍ⁴ ℍ-linéaire", "sonde-couleur-mémoire-quaternionique", (97, 23), "Broyden-mêmes-mesures", (100, 61), True, "5"),
]
ROWS_APRES = [   # parties 5.5 et 6, dans l'ordre du rapport, après la partie 5
    ("ℍ¹⁶ ℍ-linéaire, même structure pour tous", "sonde-couleur-mémoire-quaternionique", (82, 100), "Broyden-quaternionique-sans-couleur", (100, 34), False, "5.5"),
    ("ℍ⁶⁴ couche 8 × 8 petit signal, même structure", "sonde-couleur-mémoire-quaternionique-à-droite-par-sortie", (100, 18), "Broyden-quaternionique-sans-couleur", (100, 14), False, "5.5"),
    ("ℍ⁸ point fixe RNN, pièce apprise", "sonde-couleur-mémoire-quaternionique", (100, 17), "Broyden-sans-couleur-départ-identité", (100, 9), False, "6.1"),
    ("ℍ⁴ inversion réseau, saturation globale", "sonde-couleur-mémoire", (48, 40), "Broyden-sans-couleur", (93, 29), False, "6.2"),
    ("ℍ⁸ entraîner un neurone (quotient à droite)", "sonde-couleur-mémoire-quaternionique-à-droite", (100, 17), "Broyden-sans-couleur", (100, 34), False, "6.3"),
    ("ℍ¹⁶ entraîner un neurone, tanh petit signal", "sonde-couleur-mémoire-quaternionique-à-droite", (93, 34), "Broyden-sans-couleur", (100, 68), False, "6.3"),
    ("ℍ⁶⁴ entraîner une couche 8 × 8, petit signal", "sonde-couleur-mémoire-quaternionique-à-droite-par-sortie", (100, 18), "Broyden-sans-couleur-par-sortie", (100, 37), False, "6.4"),
]
ROWS_43 = [
    ("ℍ⁴ bilatère, blocs conformes", "sonde-couleur-mémoire-conforme-par-blocs", (97, 104), "Broyden-sans-couleur", (100, 24), False, "4.3"),
    ("ℍ⁴ ℍ-linéaire + rang 2", "sonde-couleur-mémoire-quaternionique+rang-2", (97, 72), "Broyden-sans-couleur", (100, 31), False, "4.3"),
]


def rows9():
    p = os.path.join(HERE, "..", "parties", "4_autres_structures", "4.1_hors_structure", "bench_hors_structure_results.json")
    if not os.path.exists(p): return []
    R = json.load(open(p))["results"]; n = max(r["n"] for r in R); out = []
    lab = dict(poly="polynôme non commutatif", qnn="inversion réseau quaternionique", qres="point fixe réservoir quaternionique", bilat="bilatère a·d·b")
    for fam in ("poly", "qnn", "qres", "bilat"):
        rs = [r for r in R if r["fam"] == fam and r["n"] == n and r["sigma"] == 0.0]
        if not rs: continue
        probes = [r for r in rs if r["method"].startswith("sec-") and r["cost"]]
        best = max(probes, key=lambda r: (r["rate"] >= 0.9, -r["cost"] if r["rate"] >= 0.9 else r["rate"]))
        ref = [r for r in rs if r["method"] == "broyden"][0]
        out.append((f"ℍ{'⁴' if n == 4 else n} {lab[fam]}", nom(best["method"]), (round(100 * best["rate"]), best["cost"]), "Broyden-sans-couleur",
                    (round(100 * ref["rate"]), ref["cost"]), False, "4.1"))
    return out


def bilan():
    rows = ROWS[:-5] + rows9() + ROWS[-5:-3] + ROWS_43 + ROWS[-3:] + ROWS_APRES     # dans l'ordre du rapport : 1 à 6
    fig, ax = plt.subplots(1, 2, figsize=(17, 0.45 * len(rows) + 1.8), gridspec_kw=dict(width_ratios=[0.75, 1.5]))
    y = np.arange(len(rows))[::-1]
    for yi, r in zip(y, rows):
        ax[0].barh(yi + 0.18, r[2][0], 0.34, color=PAL["blue"]); ax[0].barh(yi - 0.18, r[4][0], 0.34, color=PAL["orange"])
        ratio = r[2][1] / r[4][1]
        ax[1].barh(yi, np.log2(ratio), 0.55, color=PAL["aqua"] if ratio < 1 else PAL["red"])
        ax[1].text(max(np.log2(ratio), 0) + 0.06, yi, f"×{ratio:.2f}  ({r[2][1]:.0f} contre {r[4][1]:.0f})",
                   va="center", ha="left", fontsize=7.5, color=INK2)
    ax[0].set_yticks(y); ax[0].set_yticklabels([("★ " if r[5] else "") + f"{r[0]}  [{r[6]}]" for r in rows], fontsize=8.5)
    ax[0].set_xlim(0, 105); ax[0].set_title("succès (%) : sonde (bleu), référence (orange)")
    ax[1].set_yticks(y); ax[1].set_yticklabels([f"{r[1]} / {r[3]}" for r in rows], fontsize=7.5)
    ax[1].axvline(0, color=INK2, lw=0.8); ax[1].set_xlim(-1.5, 3.6)
    ax[1].set_xticks([-1, 0, 1, 2]); ax[1].set_xticklabels(["÷2", "1", "×2", "×4"])
    ax[1].set_title("coût relatif en lectures, sonde / référence (vert : sonde moins chère)")
    fig.suptitle("Le bilan en un graphe : la meilleure sonde couleur face à la meilleure référence classique   "
                 "(★ mêmes mesures pour tous ; * |f| connu au point courant ; [n] = partie)", x=0.01, ha="left", fontsize=10)
    fig.tight_layout(); fig.savefig(os.path.join(HERE, "07_bilan.png")); plt.close(fig)


if __name__ == "__main__":
    fil(); bilan(); print("ok")
