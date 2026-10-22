"""Schéma 8 : de Broyden à la sonde, un ingrédient à la fois (problème ℍ-linéaire des parties 3.2 et 5.5, sans bruit).

Chaque barre ajoute ou retire une seule chose par rapport à la précédente :
  Broyden-sans-couleur                          valeur exacte, aucune structure (4n lectures de départ)
  + structure ℍ-linéaire                        Broyden-quaternionique-sans-couleur (n lectures de départ)
  + mesure à 24 bits, décade donnée             Broyden-quaternionique-mêmes-mesures
  + lecture par quotients de couleurs           sonde-oracle-décade (même information que la ligne précédente)
  + décade lue par le capteur à deux échelles   sonde-deux-échelles + test géométrique (bits égaux, plus d'oracle)
  décade inconnue (dent de scie)                sonde-couleur-mémoire-quaternionique de base
Les chiffres sont ceux de parties/5_armes_egales/bench_meme_structure_tables.md (σ = 0, mêmes 100 départs pour toutes les lignes).
Usage : python schemas/make_conclusion.py   → schemas/08_de_broyden_a_la_sonde.png"""
import os, sys
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "src"))
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from report import PAL, style, INK, INK2, MUTED, GRID
style()

# (étiquette, méthode du rapport, {n: (succès %, lectures médianes)}) ; tableaux de la partie 5.5, σ = 0, 100 départs
ETAPES = [
    ("valeur exacte, sans structure", "Broyden-sans-couleur", {16: (100, 84), 32: (100, 152)}),
    ("+ structure ℍ-linéaire", "Broyden-quaternionique-sans-couleur", {16: (100, 34), 32: (100, 53)}),
    ("+ mesure quantifiée à 24 bits\n   (décade donnée)", "Broyden-quaternionique-mêmes-mesures", {16: (100, 45), 32: (100, 70)}),
    ("+ quotients de couleurs\n   (même information)", "sonde-oracle-décade", {16: (100, 53), 32: (95, 89)}),
    ("+ décade lue par le capteur\n   à deux échelles (bits égaux)", "sonde-deux-échelles + test géométrique", {16: (97, 71), 32: (99, 122)}),
    ("décade inconnue\n   (dent de scie)", "sonde-couleur-mémoire-quaternionique", {16: (82, 100), 32: (86, 172)}),
]
COUL = [MUTED, PAL["blue"], PAL["blue"], PAL["orange"], PAL["orange"], PAL["orange"]]

fig, axes = plt.subplots(1, 2, figsize=(12.5, 5.2), sharey=True)
for ax, n in zip(axes, (16, 32)):
    y = list(range(len(ETAPES)))[::-1]
    vals = [e[2][n][1] for e in ETAPES]
    ax.barh(y, vals, color=COUL, height=0.62, edgecolor="white", linewidth=2)
    for yy, (lab, meth, d), v in zip(y, ETAPES, vals):
        s = d[n][0]
        ax.text(v + max(vals) * 0.015, yy, f"{v} lectures" + ("" if s == 100 else f"  ({s} % de succès)"), va="center", fontsize=9, color=INK)
        ax.text(max(vals) * 0.015, yy - 0.43, meth, va="center", fontsize=7.2, color=INK2)
    ax.axvline(vals[0], color=MUTED, lw=0.8, ls=":")
    ax.set_xlim(0, max(vals) * 1.45); ax.set_xlabel("lectures médianes jusqu'au zéro")
    ax.set_title(f"n = {n} quaternions (dimension {4 * n})", fontsize=10.5)
    ax.grid(axis="x", color=GRID); ax.grid(axis="y", visible=False)
axes[0].set_yticks(list(range(len(ETAPES)))[::-1]); axes[0].set_yticklabels([e[0] for e in ETAPES], fontsize=9)
fig.text(0.01, 0.015, "Gris : référence classique. Bleu : ce que gagne la structure (valeur lue). Orange : ce que coûte le passage à la couleur.\n"
         "Problème ℍ-linéaire des parties 3.2 et 5.5, sans bruit, 100 départs ; chaque ligne change un seul ingrédient par rapport à la précédente.",
         fontsize=8, color=INK2)
fig.suptitle("Conclusion : la structure fait le gain, la couleur le conserve en partie", x=0.01, ha="left", fontsize=12.5, color=INK)
fig.tight_layout(rect=(0, 0.07, 1, 0.95))
fig.savefig(os.path.join(HERE, "08_de_broyden_a_la_sonde.png")); plt.close(fig)
print("ok")
