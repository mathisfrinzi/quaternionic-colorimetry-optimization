"""Aides de rapport : style des figures (palette catégorielle validée, surface claire), tableaux Markdown."""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

PAL = dict(blue="#2a78d6", orange="#eb6834", aqua="#1baf7a", yellow="#eda100", magenta="#e87ba4",
           green="#008300", violet="#6250d6", red="#e34948")
INK, INK2, MUTED, GRID, SURF = "#0b0b0b", "#52514e", "#8a8984", "#e6e5e1", "#fcfcfb"


def style():
    plt.rcParams.update({
        "figure.facecolor": SURF, "axes.facecolor": SURF, "savefig.facecolor": SURF,
        "axes.edgecolor": GRID, "axes.labelcolor": INK2, "xtick.color": INK2, "ytick.color": INK2,
        "text.color": INK, "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.8,
        "axes.spines.top": False, "axes.spines.right": False, "font.size": 10, "axes.titlesize": 11,
        "axes.titlelocation": "left", "legend.frameon": False, "figure.dpi": 130, "savefig.dpi": 160,
    })
    installer_noms()


def pct(s):
    return f"{100 * s['rate']:.0f} % ({100 * s['lo']:.0f}–{100 * s['hi']:.0f})"


def md_table(headers, rows):
    out = ["| " + " | ".join(traduire(str(h)) for h in headers) + " |", "|" + "|".join("---" for _ in headers) + "|"]
    for r in rows:
        out.append("| " + " | ".join(traduire(str(c)) for c in r) + " |")
    return "\n".join(out) + "\n"


def fmt_cost(c):
    return "—" if (c is None or not np.isfinite(c)) else f"{c:.0f}"


def heatmap(ax, M, rows, cols, vmin=0, vmax=100, fmt="{:.0f}", cmap="Blues"):
    im = ax.imshow(M, aspect="auto", vmin=vmin, vmax=vmax, cmap=cmap)
    ax.set_xticks(range(len(cols))); ax.set_xticklabels(cols, rotation=35, ha="right")
    ax.set_yticks(range(len(rows))); ax.set_yticklabels(rows)
    ax.grid(False)
    for i in range(M.shape[0]):
        for j in range(M.shape[1]):
            v = M[i, j]
            if np.isfinite(v):
                ax.text(j, i, fmt.format(v), ha="center", va="center", fontsize=8,
                        color="white" if (v - vmin) / (vmax - vmin) > 0.6 else INK)
    return im


# ------------------------------------------------------------------------------------------------ noms du rapport
# Le code garde ses noms courts (sec-H, fd-full, broyden-deg…) ; les figures et les tableaux affichent des noms en mots.
# La correspondance complète est donnée dans parties/0_capteurs_et_sondes/README.md (glossaire).
# Règle : « sonde-… » = notre algorithme ; « …-sans-couleur » = méthode classique qui lit la valeur exacte de f ;
# « …-mêmes-mesures » = méthode classique qui lit la même mesure dégradée que la sonde (partie 5).
NOMS = {
    # capteurs
    "S0": "capteur-valeurs", "S1": "couleur-sans-axe", "S2": "couleur-complète",
    "S2rgb": "couleur-complète-RVB", "S2sph": "couleur-complète-angles", "exacte": "valeur-exacte",
    # dimension 4 : sondes
    "sonde-S0": "sonde-sans-couleur", "sonde-S1": "sonde-couleur-sans-axe",
    "sonde-S2rgb": "sonde-couleur-RVB", "sonde-S2sph": "sonde-couleur-angles", "sonde-exacte": "sonde-valeur-exacte",
    "fd-full": "sonde-couleur", "fd-conf3": "sonde-couleur-conforme-3-sondes", "fd-confproj": "sonde-couleur-conforme-projetée",
    "sec-rank1": "sonde-couleur-mémoire", "sec-multi": "sonde-couleur-mémoire-multiple", "sec-conf": "sonde-couleur-mémoire-conforme",
    "sonde-rank1": "sonde-couleur-mémoire", "sonde-fd": "sonde-couleur",
    "sonde-base": "sonde-couleur-mémoire", "sonde-dith": "sonde-couleur-mémoire + dithering",
    "sonde-conf": "sonde-couleur-mémoire + démarrage-conforme", "sonde-tr": "sonde-couleur-mémoire + région-de-confiance",
    "sonde-dith+conf+tr": "sonde-couleur-mémoire + les trois",
    "complexe-original": "sonde-complexe (dépôt d'origine)", "S1-signé": "sonde-couleur-sans-axe-signée",
    "S1-non-signé": "sonde-couleur-sans-axe", "S2-4D": "sonde-couleur (4D)",
    "cible": "sonde-couleur-à-cible",
    # dimension 4n
    "fd-H": "sonde-couleur-quaternionique", "sec-full": "sonde-couleur-mémoire",
    "sec-H": "sonde-couleur-mémoire-quaternionique", "sec-H2full": "sonde-couleur-mémoire-quaternionique-puis-libre",
    "sec-Hinit": "sonde-projeter-puis-affiner", "sec-Hid": "sonde-couleur-mémoire-quaternionique-départ-identité",
    # partie 3.3
    "base": "sonde-couleur-mémoire-quaternionique", "lum_exacte": "sonde-oracle-luminosité", "desamb": "sonde-oracle-décade",
    "W_vraie": "sonde-oracle-jacobienne", "prior0": "sonde-dépliage-vers-zéro", "noflag": "sonde-sans-drapeau",
    # méthodes classiques, valeur exacte
    "broyden": "Broyden-sans-couleur", "broyden-4D": "Broyden-sans-couleur", "broyden-id": "Broyden-sans-couleur-départ-identité",
    "newton-FD": "Newton-sans-couleur", "newton-4D-FD": "Newton-sans-couleur", "newton-Jexacte": "Newton-jacobienne-exacte",
    "newton-4D-Jexacte": "Newton-jacobienne-exacte", "newton-complexe": "Newton-complexe-sans-couleur",
    "bfgs-FD": "BFGS-sans-couleur", "bfgs-4D": "BFGS-sans-couleur", "lm-FD": "Levenberg-Marquardt-sans-couleur",
    "nelder-mead": "Nelder-Mead-sans-couleur", "nelder-mead-4D": "Nelder-Mead-sans-couleur", "powell": "Powell-sans-couleur",
    "cma-es": "CMA-ES-sans-couleur", "cma-es-4D": "CMA-ES-sans-couleur", "picard": "Picard-sans-couleur",
    # parties 3.5, 4.3 et 6 (src/sonde_hn.py)
    "sonde-H": "sonde-couleur-mémoire-quaternionique", "sonde-full": "sonde-couleur-mémoire",
    "sonde-Hid": "sonde-couleur-mémoire-quaternionique-départ-identité", "sonde-H-relecture": "sonde-relecture",
    "sonde-blocconf": "sonde-couleur-mémoire-conforme-par-blocs", "sonde-Hrang": "sonde-couleur-mémoire-quaternionique+rang-2",
    "sonde-grise": "sonde-couleur-mémoire-boîte-grise", "sonde-D": "sonde-couleur-mémoire-quaternionique-à-droite",
    "sonde-relecture-econome": "sonde-relecture-économe", "sonde-vote": "sonde-vote",
    "sonde-deux-echelles": "sonde-deux-échelles", "sonde-deux-echelles-geo": "sonde-deux-échelles + test géométrique",
    "sonde-deux-echelles-pond": "sonde-deux-échelles + test pondéré",
    "sonde-garde-ambiguite": "sonde-garde-ambiguïté", "sonde-depliage-vers-zero": "sonde-dépliage-vers-zéro",
    "R1-relecture": "sonde-relecture", "sonde-oracle-decade": "sonde-oracle-décade",
    "sonde-oracle-decade-geo": "sonde-oracle-décade + test géométrique", "sonde-oracle-decade-pond": "sonde-oracle-décade + test pondéré",
    "broyden-sans-couleur": "Broyden-sans-couleur",
    "bq-deg2": "Broyden-quaternionique-mesure-absolue-à-bits-égaux", "sonde-2e": "sonde-deux-échelles", "bq-sans-suivi": "Broyden-quaternionique-sans-couleur (sans variante suivi)", "bq": "Broyden-quaternionique-sans-couleur", "bq-deg": "Broyden-quaternionique-mêmes-mesures", "sonde-oracle": "sonde-oracle-décade", "sonde-2e-geo": "sonde-deux-échelles + test géométrique", "sonde-Dbloc": "sonde-couleur-mémoire-quaternionique-à-droite-par-sortie", "broyden-blocs": "Broyden-sans-couleur-par-sortie",
    # méthodes classiques, mêmes mesures
    "broyden-deg": "Broyden-mêmes-mesures",
}
_SURCHARGE = {}


def set_noms(**surcharge):
    """Surcharge locale (ex. partie 5 : `broyden` y lit la mesure dégradée, donc « Broyden-mêmes-mesures »)."""
    _SURCHARGE.update(surcharge)


import re as _re
_SUFFIXE = _re.compile(r"^(.*?)(-(\d+)b?)?(@.*)?$")


def nom(code):
    """Nom du rapport pour un code de méthode, avec ses suffixes éventuels (-24b, @24, @24/0.02, -Broyden)."""
    if not isinstance(code, str):
        return code
    table = {**NOMS, **_SURCHARGE}
    if code in table:
        return table[code]
    base, broy = (code[:-len("-Broyden")], True) if code.endswith("-Broyden") else (code, False)
    m = _SUFFIXE.match(base)
    tete, bits, arob = m.group(1), m.group(3), m.group(4)
    if tete not in table:
        return code
    out = table[tete]
    if broy:
        out = out.replace("sonde-couleur", "sonde-couleur-mémoire", 1)
    extra = []
    if bits: extra.append(f"{bits} bits")
    if arob:
        a = arob[1:]
        if "/" in a:
            b, s = a.split("/"); extra.append(f"{b} bits, σ = {s}")
        else:
            extra.append(f"{a} bits")
    return out + (f" ({', '.join(extra)})" if extra else "")


# mots ordinaires qui sont aussi des codes : jamais traduits dans un texte libre (seulement par nom())
_PAS_EN_TEXTE = {"base", "exacte", "cible", "S0", "S1", "S2"}
_PHRASES = [("S2 sphérique", "couleur-complète-angles"), ("S2 RGB", "couleur-complète-RVB"), ("S2 sph", "couleur-complète-angles")]
_JETON = _re.compile(r"(?<![\w-])([A-Za-z][\w+-]*(?:-\d+b)?(?:-Broyden)?(?:@[\d./]+)?)(?![\w-])")


def traduire(texte):
    """Remplace dans un texte libre chaque code de méthode connu par son nom de rapport."""
    if not isinstance(texte, str):
        return texte
    for a, b in _PHRASES:
        texte = texte.replace(a, b)
    def rep(m):
        t = m.group(1)
        if t in _PAS_EN_TEXTE:
            return t
        n = nom(t)
        return n if n != t else t
    return _JETON.sub(rep, texte)


def installer_noms():
    """Toute chaîne affichée par matplotlib (légendes, titres, graduations) passe par `traduire`."""
    import matplotlib.text as mt
    if getattr(mt.Text, "_noms_installes", False):
        return
    orig = mt.Text.set_text
    def set_text(self, s):
        return orig(self, traduire(s) if isinstance(s, str) else s)
    mt.Text.set_text = set_text
    mt.Text._noms_installes = True


installer_noms()          # actif dès l'import : tout script de figures affiche les noms du rapport
