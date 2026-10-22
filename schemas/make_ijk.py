"""Figures pédagogiques « avec / sans » dans l'espace (i, j, k).

Pour chaque élément ajouté à la sonde au fil du rapport, on lance la sonde SANS puis AVEC cet élément, sur les mêmes départs, et
on dessine un départ où l'élément aide (choisi parmi ceux où il aide : le gain médian de ces départs ; le nombre de départs où il
aide est écrit sous la figure, pour ne pas faire passer un cas choisi pour la règle).
  à gauche : la trajectoire dans l'espace (i, j, k), parties imaginaires de q − q* (en ℍⁿ : de la composante quaternionique la plus
             loin du zéro au départ). Chaque point est un pas accepté (une lecture qui fait baisser ‖f‖) ; la couleur du point
             est la partie réelle (axe 1), le zéro est l'étoile.
  à droite : ‖f‖ (meilleure valeur lue jusque-là) en fonction du nombre de lectures, avec la même référence classique grise.
Usage : python schemas/make_ijk.py [nom …]   (sans argument : toutes les figures)
Sortie : parties/<partie>/figures/ijk_<nom>.png, et schemas/ijk_resume.json (départs dessinés, nombre de départs où l'élément aide)."""
import os, sys, json
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
R = os.path.join(HERE, "..")
sys.path.insert(0, os.path.join(R, "src"))
for d in ("4_autres_structures/4.1_hors_structure", "4_autres_structures/4.2_reservoirs", "6_reseaux",
          "3_dimension_4n/3.2_grande_dimension", "3_dimension_4n/3.1_jacobienne_H_lineaire"):
    sys.path.insert(0, os.path.join(R, "parties", d))
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from report import PAL, style, INK, INK2, MUTED, GRID, SURF
from sensors import Sensor
from funcs import build
from probes import Meter, Hit, OutOfBudget, probe_fd, probe_secant, probe_s1
import classic
import sonde_hn as S
style()
ERRS = S.ERRS
PARTIES = os.path.join(R, "parties")


# ================================================================================================ lancer en dimension 4
class Meter4(Meter):
    def __init__(self, *a, **k):
        super().__init__(*a, **k); self.log = []

    def read(self, q):
        F = self.f(q)
        if np.all(np.isfinite(F)):
            self.log.append((self.n + 1, float(np.linalg.norm(F)), np.array(q, float)))
        return super().read(q)


def lancer4(P, q0, fn, sensor=None, absolute=False, budget=600, **kw):
    M = Meter4(P.f, sensor, budget=budget, thr=P.thr, absolute=absolute)
    try:
        fn(M, q0, **kw); ok = False
    except Hit:
        ok = True
    except ERRS:
        ok = False
    return ok, (M.n if ok else budget), M.log


# ================================================================================================ dessin
def chemin(log):
    """Pas acceptés ≈ lectures qui battent le record de ‖f‖ (une sonde n'accepte que des pas qui font baisser ‖f‖)."""
    best, pts = np.inf, []
    for c, nf, q in log:
        if nf < best:
            best = nf; pts.append((c, nf, q))
    return pts


def courbe(log):
    best, xs, ys = np.inf, [], []
    for c, nf, q in log:
        best = min(best, nf); xs.append(c); ys.append(best)
    return np.array(xs), np.array(ys)


def dessiner(path, titre, sans, avec, ref, qs, n, note, labels, seuil=1e-3):
    """sans / avec / ref : (ok, lectures, log)."""
    fig = plt.figure(figsize=(11.5, 4.6))
    ax = fig.add_subplot(1, 2, 1, projection="3d")
    q0 = sans[2][0][2]
    d0 = (q0 - qs).reshape(n, 4)
    j = int(np.argmax(np.linalg.norm(d0, axis=1)))
    allw = []
    for run, col, lab in ((sans, PAL["blue"], labels[0]), (avec, PAL["orange"], labels[1])):
        pts = chemin(run[2])
        D = np.array([(p[2] - qs).reshape(n, 4)[j] for p in pts])
        allw += list(D[:, 0])
        ax.plot(D[:, 1], D[:, 2], D[:, 3], color=col, lw=1.8, alpha=0.9,
                label=f"{lab} : {'succès' if run[0] else 'échec'}, {run[1]} lectures")
        ax.scatter(D[:, 1], D[:, 2], D[:, 3], c=D[:, 0], cmap="coolwarm", s=16, edgecolors=col, linewidths=0.6, depthshade=False)
    ax.scatter([d0[j, 1]], [d0[j, 2]], [d0[j, 3]], color=INK, s=40, marker="o", label="départ")
    ax.scatter([0], [0], [0], color=PAL["yellow"], s=170, marker="*", edgecolors=INK, linewidths=0.6, label="zéro q*" if seuil == 1e-3 else "minimum q*")
    ax.set_xlabel("i"); ax.set_ylabel("j"); ax.set_zlabel("k")
    for a in (ax.xaxis, ax.yaxis, ax.zaxis):
        a.pane.set_facecolor(SURF); a.pane.set_edgecolor(GRID)
    sub = "" if n == 1 else f" (composante {j + 1} sur {n})"
    ax.set_title("trajectoire dans (i, j, k) de q − q*" + sub + "\ncouleur des points : partie réelle (bleu < 0 < rouge)", fontsize=9.5)
    ax.legend(loc="upper left", fontsize=7.5, bbox_to_anchor=(-0.05, 1.02))
    ax.view_init(elev=24, azim=-58)
    from matplotlib.ticker import MaxNLocator
    for a in (ax.xaxis, ax.yaxis, ax.zaxis):
        a.set_major_locator(MaxNLocator(4)); a.set_tick_params(labelsize=7.5)
    ax2 = fig.add_subplot(1, 2, 2)
    for run, col, lab, ls in ((ref, MUTED, labels[2], ":"), (sans, PAL["blue"], labels[0], "-"), (avec, PAL["orange"], labels[1], "-")):
        if run is None: continue
        x, y = courbe(run[2])
        ax2.step(x, y, where="post", color=col, lw=2 if ls == "-" else 1.6, ls=ls, label=f"{lab} ({run[1]} lectures{'' if run[0] else ', échec'})")
    ax2.axhline(seuil, color=INK2, lw=0.8, ls="--")
    ax2.text(1, seuil * 1.3, "succès : ‖f‖ < 10⁻³" if seuil == 1e-3 else f"succès : ‖f‖ < {seuil:g} (minimum 1)".replace(".", ","), color=INK2, fontsize=8)
    ax2.set_yscale("log"); ax2.set_xlabel("lectures"); ax2.set_ylabel("meilleur ‖f‖ lu")
    ax2.set_title("convergence", fontsize=10); ax2.legend(fontsize=8, loc="best", frameon=True, facecolor=SURF, edgecolor=GRID, framealpha=0.92)
    fig.suptitle(titre, x=0.02, ha="left", fontsize=12, color=INK)
    fig.text(0.02, 0.01, note, fontsize=8, color=INK2)
    fig.tight_layout(rect=(0, 0.04, 1, 0.94))
    os.makedirs(os.path.dirname(path), exist_ok=True)
    fig.savefig(path); plt.close(fig)


def choisir(paires):
    """paires : liste de (sans, avec) = ((ok, c, log), (ok, c, log)). Retourne (indice, nb où l'élément aide, nb total).
    « aide » : avec réussit et (sans échoue ou coûte plus). On dessine le départ de gain médian parmi ceux-là."""
    aide = []
    for k, (s, a) in enumerate(paires):
        if a[0] and (not s[0] or a[1] < s[1]):
            gain = (s[1] - a[1]) / s[1] if s[0] else 1.0 + (600 - a[1]) / 600
            aide.append((gain, k))
    if not aide:
        return None, 0, len(paires)
    aide.sort()
    return aide[len(aide) // 2][1], len(aide), len(paires)


RESUME = {}


def element(nom, partie, titre, labels, gen, k=12):
    """gen(i) -> dict(sans=…, avec=…, ref=…, qs=…, n=…)."""
    runs = [gen(i) for i in range(k)]
    idx, na, nt = choisir([(r["sans"], r["avec"]) for r in runs])
    if idx is None:
        print(f"{nom}: l'élément n'aide sur aucun départ, figure non faite"); return
    r = runs[idx]
    gm = [((a["sans"][1] - a["avec"][1]) / a["sans"][1]) for a in runs if a["sans"][0] and a["avec"][0]]
    ss = sum(a["sans"][0] for a in runs); sa = sum(a["avec"][0] for a in runs)
    note = (f"Départ dessiné : n° {idx} sur {nt}. L'élément aide sur {na} départs sur {nt} ; succès sans / avec : {ss}/{nt} et {sa}/{nt} ; "
            f"gain médian en lectures quand les deux réussissent : {100*np.median(gm):+.0f} %." if gm else
            f"Départ dessiné : n° {idx} sur {nt}. L'élément aide sur {na} départs sur {nt} ; succès sans / avec : {ss}/{nt} et {sa}/{nt}.")
    path = os.path.join(PARTIES, partie, "figures", f"ijk_{nom}.png")
    dessiner(path, titre, r["sans"], r["avec"], r.get("ref"), r["qs"], r["n"], note, labels, r.get("seuil", 1e-3))
    RESUME[nom] = dict(partie=partie, depart=idx, aide=na, total=nt, succes_sans=ss, succes_avec=sa,
                       gain_median=float(np.median(gm)) if gm else None,
                       lectures=dict(sans=r["sans"][1], avec=r["avec"][1], ref=None if r.get("ref") is None else r["ref"][1]))
    print(f"{nom}: départ {idx}, aide {na}/{nt}, lectures sans {r['sans'][1]} avec {r['avec'][1]}", flush=True)


# ================================================================================================ éléments, dimension 4
def sen(seed, **kw):
    return Sensor(seed=seed, **{**dict(kind="s2", bits=24, enc="sph"), **kw})


def ref4(P, q0):
    return lancer4(P, q0, classic.broyden, hfd=1e-6)


def g_axe(i):
    """Partie 0/1.2 : l'axe. couleur-sans-axe (2 équations sur 4) → couleur-complète (quotient complet), problème mal conditionné."""
    P = build(("gen", "ill", i % 3)); q0 = P.starts(4, 100 + i)[i // 3]
    return dict(sans=lancer4(P, q0, probe_s1, sen(i, kind="s1", bits=48)), avec=lancer4(P, q0, probe_fd, sen(i, bits=48, enc="rgb")),
                ref=ref4(P, q0), qs=P.zero, n=1)


def g_memoire(i):
    P = build(("gen", ["lin", "quad", "cub", "ill"][i % 4], i // 4)); q0 = P.starts(3, 200 + i)[0]
    return dict(sans=lancer4(P, q0, probe_fd, sen(i, bits=48, enc="rgb")), avec=lancer4(P, q0, probe_secant, sen(i, bits=48, enc="rgb")),
                ref=ref4(P, q0), qs=P.zero, n=1)


def g_conforme(i):
    P = build(("mobius", i % 4, 0.0)); q0 = P.starts(3, 300 + i)[0]
    return dict(sans=lancer4(P, q0, probe_fd, sen(i), model="full"), avec=lancer4(P, q0, probe_fd, sen(i), model="conf3"),
                ref=ref4(P, q0), qs=P.zero, n=1)


def g_dither(i):
    """Sur des fonctions faciles, le dithering ne change presque rien (mêmes lectures à une près) : on prend la famille mal
    conditionnée, où l'erreur de quantification pèse."""
    P = build(("gen", "ill", i % 3)); q0 = P.starts(6, 400 + i // 3)[i // 3]
    return dict(sans=lancer4(P, q0, probe_secant, sen(i, bits=24)), avec=lancer4(P, q0, probe_secant, sen(i, bits=24, dither=True)),
                ref=ref4(P, q0), qs=P.zero, n=1)


def g_cible(i):
    """Partie 2 : minimum non nul (‖f‖ ≥ 1). Sonde à cible scalaire du dépôt complexe (luminosité de ‖f‖² − m₀) → sonde à cible
    quaternionique (couleur complète de F̃ = R(‖R‖ − m₀)/‖R‖). Mêmes hypothèses que la partie 2 (‖f‖ connu au point courant)."""
    sys.path.insert(0, os.path.join(PARTIES, "2_minimums")); sys.path.insert(0, os.path.join(R, "vendor", "complex_colorimetry"))
    from bench_minimum import cible
    import bench_ablation as BA
    P = build(("gen", "nozero_ill", i % 3)); q0 = P.starts(4, 1)[i // 3]
    sans = lancer4(P, q0, cible, bits=36, c=0.1, L=1.0)
    log = []
    class Pl:                                            # même problème, journal des lectures COMPTÉES par la sonde
        thr = P.thr
        @staticmethod
        def f(q):
            return P.f(q)
    sen0 = Sensor(kind="s2", enc="sph", bits=36, seed=1000 + i)
    class SenLog:
        def __init__(self, s): self.s = s
        def read(self, Ft):
            log.append((len(log) + 1, float(np.linalg.norm(P.f(SenLog.q))), np.array(SenLog.q))); return self.s.read(Ft)
        def quotient(self, *a, **k): return self.s.quotient(*a, **k)
    orig = BA.Sensor
    def run_log(q0):
        # BA.run lit via read(q, m0) → sen.read(F̃) ; on mémorise q avant chaque lecture en enveloppant Rtrue
        Pw = type("Pw", (), {})(); Pw.thr = P.thr
        def fw(q):
            SenLog.q = np.array(q, float); return P.f(q)
        Pw.f = fw
        return BA.run(Pw, q0, SenLog(sen0), 0.1)
    ok, c = run_log(q0)
    if ok:                                               # la lecture du succès s'arrête avant le capteur : on l'ajoute au journal
        log.append((c, float(np.linalg.norm(P.f(SenLog.q))), np.array(SenLog.q)))
    # BA.run évalue aussi f gratuitement (‖R‖ supposé connu) : seules les lectures passées au capteur sont dans `log`
    return dict(sans=sans, avec=(ok, c, log), ref=lancer4(P, q0, classic.lm_fd, diff_step=1e-6), qs=P.zero, n=1, seuil=1.05)


# ================================================================================================ éléments, dimension 4n
def lancerN(f, q0, n, genre="sonde", sigma=0.0, i=0, budget=2000, **kw):
    if genre == "sonde":
        kw = {**S.__dict__.get("P0", {}), **dict(c=0.3, step=0.65, cap=1.5 * max(1.0, np.sqrt(n / 4))), **kw}
    ok, c, log, *_ = S.lancer(f, q0, n, genre, budget=budget, sigma=sigma, i=i, log=True, **kw)
    return ok, c, log


def refN(f, q0, n):
    return lancerN(f, q0, n, "broyden", hfd=1e-6)


def _hn(n, i, eps=0.0):
    import bench_hn_grand as G
    seed = 100 + i // 10
    f, qs = G.problem(n, eps, seed); q0 = G.starts(n, qs, i % 10 + 1, seed)[i % 10]
    return f, qs, q0


def g_quaternionique(i):
    n = 4; f, qs, q0 = _hn(n, i)
    return dict(sans=lancerN(f, q0, n, modele="full", i=i, budget=600), avec=lancerN(f, q0, n, modele="H", i=i, budget=600),
                ref=refN(f, q0, n), qs=qs, n=n)


def g_depliage(i):
    n = 16; f, qs, q0 = _hn(n, i)
    return dict(sans=lancerN(f, q0, n, i=i), avec=lancerN(f, q0, n, lum="prior0", i=i), ref=refN(f, q0, n), qs=qs, n=n)


def g_relecture(i):
    n = 16; f, qs, q0 = _hn(n, i)
    return dict(sans=lancerN(f, q0, n, sigma=0.02, i=i), avec=lancerN(f, q0, n, sigma=0.02, i=i, verif_k=2, tau=0.25, double=True),
                ref=refN(f, q0, n), qs=qs, n=n)


def g_deux_echelles(i):
    n = 16; f, qs, q0 = _hn(n, i)
    return dict(sans=lancerN(f, q0, n, sigma=0.02, i=i, verif_k=2, tau=0.25, double=True),
                avec=lancerN(f, q0, n, sigma=0.02, i=i, deux_echelles=True), ref=refN(f, q0, n), qs=qs, n=n)


def g_test_geo(i):
    n = 16; f, qs, q0 = _hn(n, i)
    return dict(sans=lancerN(f, q0, n, sigma=0.02, i=i, deux_echelles=True),
                avec=lancerN(f, q0, n, sigma=0.02, i=i, deux_echelles=True, test="geo"), ref=refN(f, q0, n), qs=qs, n=n)


def g_identite(i):
    from reservoir import Reservoir
    n = 8; res = Reservoir(n, rho=0.95, beta=1.0, kind="global", seed=i % 3)
    rng = np.random.default_rng(8000 + i); v = rng.standard_normal(4 * n); v /= np.linalg.norm(v)
    q0 = res.hs + 0.5 * np.sqrt(n) * v
    return dict(sans=lancerN(res.F, q0, n, modele="H", i=i), avec=lancerN(res.F, q0, n, modele="Hid", i=i),
                ref=refN(res.F, q0, n), qs=res.hs, n=n)


def _p43(prob, n, i):
    import banc_hn as B
    sys.path.insert(0, os.path.join(PARTIES, "4_autres_structures", "4.3_structures_etendues"))
    import bench_structures  # enregistre les problèmes
    f, qs, extra = B.PROBLEMES[prob](n, i % 3)
    return f, qs, B.depart(qs, n, i, i % 3, bench_structures.RAD[prob]), extra


def g_blocconf(i):
    n = 3; f, qs, q0, _ = _p43("bilat", n, i)
    return dict(sans=lancerN(f, q0, n, modele="H", i=i, budget=600), avec=lancerN(f, q0, n, modele="blocconf", i=i, budget=600),
                ref=refN(f, q0, n), qs=qs, n=n)


def g_rang(i):
    n = 4; f, qs, q0, _ = _p43("hrang", n, i)
    return dict(sans=lancerN(f, q0, n, modele="H", i=i, budget=600), avec=lancerN(f, q0, n, modele="Hrang", k=2, i=i, budget=600),
                ref=refN(f, q0, n), qs=qs, n=n)


def g_grise(i):
    n = 3; f, qs, q0, extra = _p43("qnn", n, i)
    return dict(sans=lancerN(f, q0, n, modele="H", i=i, budget=600), avec=lancerN(f, q0, n, modele="grise", i=i, budget=600, **extra),
                ref=refN(f, q0, n), qs=qs, n=n)


def g_droite(i):
    from bench_entrainement import Neurone
    import banc_hn as B
    n = 4; m = Neurone(n, "radial-faible", i % 3); q0 = B.depart(m.qs, n, i, i % 3, 0.5)
    return dict(sans=lancerN(m.f, q0, n, modele="H", i=i, budget=600), avec=lancerN(m.f, q0, n, modele="D", i=i, budget=600),
                ref=refN(m.f, q0, n), qs=m.qs, n=n)


def g_par_sortie(i):
    """Partie 6.4 : couche de m = 4 sorties, p = 4 poids : structure à droite seule (m·p sondes) → + diagonale par sortie (p sondes)."""
    sys.path.insert(0, os.path.join(PARTIES, "6_reseaux"))
    from bench_couche import Couche
    import banc_hn as B
    m, p = 4, 4; n = m * p; L = Couche(m, p, "radial-faible", i % 3); q0 = B.depart(L.qs, n, i, i % 3, 0.5)
    return dict(sans=lancerN(L.f, q0, n, modele="D", i=i), avec=lancerN(L.f, q0, n, modele="Dbloc", m=m, i=i),
                ref=refN(L.f, q0, n), qs=L.qs, n=n)


def g_activation(i):
    """Partie 6.2 : ce n'est pas la sonde qui change mais le réseau : même sonde, tanh par composante → saturation globale."""
    from rnn import QNN
    import banc_hn as B
    n = 4; a, b = QNN(n, "split", i % 3), QNN(n, "global", i % 3)
    q0 = B.depart(a.qs, n, i, i % 3, 1.0)
    return dict(sans=lancerN(a.f, q0, n, modele="H", i=i, budget=600), avec=lancerN(b.f, q0, n, modele="H", i=i, budget=600),
                ref=refN(b.f, q0, n), qs=a.qs, n=n)


def g_sature(i):
    """Partie 6.5 : même sonde-oracle-décade sur un neurone tanh radial (n = 8), réglage de la partie → meilleur réglage du diagnostic
    (lu dans bench_diag_sature_results.json : le plus haut taux de succès à n = 8, à égalité le moins cher)."""
    sys.path.insert(0, os.path.join(PARTIES, "6_reseaux"))
    from bench_entrainement import Neurone
    import bench_diag_sature as D
    import banc_hn as B
    pj = os.path.join(PARTIES, "6_reseaux", "bench_diag_sature_results.json")
    R_ = [r for r in json.load(open(pj)) if r["kind"] == "radial" and r["n"] == 8 and r["meth"] in D.VARIANTES]
    sc = {m: (np.mean([r["ok"] for r in R_ if r["meth"] == m]), -np.median([r["cost"] for r in R_ if r["meth"] == m and r["ok"]] or [1e9])) for m in D.VARIANTES}
    meilleur = max(sc, key=lambda m: sc[m]); n = 8; seed = i % 3
    N = Neurone(n, "radial", seed); q0 = B.depart(N.qs, n, i, seed, 0.5)
    kw0 = {**D.P0, "cap": D.P0["cap"] * np.sqrt(n / 4)}
    kw1 = {**kw0, **D.VARIANTES[meilleur]}
    G_SATURE["meilleur"] = meilleur
    return dict(sans=lancerN(N.f, q0, n, modele="D", lum="desamb", f_vrai=N.f, i=i, **kw0),
                avec=lancerN(N.f, q0, n, modele="D", lum="desamb", f_vrai=N.f, i=i, **kw1),
                ref=lancerN(N.f, q0, n, "broyden-q", side="R", hfd=1e-6), qs=N.qs, n=n)


# ================================================================================================ partie 5 : la couleur face à la valeur
def g_broyden_q(i):
    """Partie 5.5 : la structure ℍ-linéaire sans la couleur. Broyden-sans-couleur (4n lectures de départ, mise à jour réelle) →
    Broyden-quaternionique-sans-couleur (n lectures de départ, mise à jour qui garde les blocs L(a)). Référence : la sonde."""
    n = 4; f, qs, q0 = _hn(n, i)
    return dict(sans=refN(f, q0, n), avec=lancerN(f, q0, n, "broyden-q", side="L", hfd=1e-6, budget=600),
                ref=lancerN(f, q0, n, modele="H", i=i, budget=600), qs=qs, n=n)


def g_instable(i):
    """Partie 5.6 : conditionnement 10³ (ℍ⁴). Même information de décade : la sonde-oracle-décade (couleur) → le Broyden quaternionique
    qui lit le capteur à deux échelles comme une valeur (mêmes 24 bits). Référence : Broyden-quaternionique-sans-couleur."""
    sys.path.insert(0, os.path.join(PARTIES, "5_armes_egales"))
    import bench_boite_instable as BI
    n, seed = 4, 100 + i % 3; f, qs, _ = BI.probleme(n, 3, 0, seed); q0 = BI.depart(qs, n, i, seed)
    return dict(sans=lancerN(f, q0, n, modele="H", lum="desamb", f_vrai=f, i=i, budget=1500),
                avec=lancerN(f, q0, n, "broyden-q-deg2", side="L", hfd=0.5, fenetre=BI.FEN, i=i, budget=1500),
                ref=lancerN(f, q0, n, "broyden-q", side="L", hfd=1e-6, budget=1500), qs=qs, n=n)


def g_amplitude(i):
    """Partie 5.7 : gain aléatoire indépendant à chaque lecture (σℓ = 0,2 népers, angle et axe exacts), ℍ⁴. La sonde-deux-échelles
    (quotients de deux lectures bruitées) → le Broyden quaternionique qui lit la valeur absolue bruitée. Référence : valeur exacte."""
    n = 4; f, qs, q0 = _hn(n, i); sg = (0.2, 0.0, 0.0)
    return dict(sans=lancerN(f, q0, n, sigma=sg, i=i, deux_echelles=True, test="geo", budget=1500),
                avec=lancerN(f, q0, n, "broyden-q-deg", sigma=sg, side="L", hfd=1.0, i=i, budget=1500),
                ref=lancerN(f, q0, n, "broyden-q", side="L", hfd=1e-6, budget=1500), qs=qs, n=n)


G_SATURE = {}


def figure_suivi(path):
    """Partie 5.8 : un zéro qui bouge. À gauche, dans (i, j, k), le zéro (courbe grise), le point lu par la sonde-deux-échelles et par la sonde
    en dent de scie ; à droite, la distance au zéro du moment à chaque lecture. Mêmes départs, mêmes lectures par pas (r = 4)."""
    sys.path.insert(0, os.path.join(PARTIES, "5_armes_egales"))
    import bench_suivi as BS
    n, v, r, seed, i = 4, 0.1, 4, 0, 1
    runs = {}
    for meth, lab, col in (("sonde", "sonde-couleur-mémoire-quaternionique (dent de scie)", PAL["blue"]),
                           ("sonde-2e", "sonde-deux-échelles (la décade se lit)", PAL["orange"]),
                           ("bq", "Broyden-quaternionique-sans-couleur", MUTED)):
        W = BS.Mobile(n, seed, v, r); budget = BS.T * r
        rng = np.random.default_rng(555 + 17 * i + seed); d = rng.standard_normal(4 * n)
        q0 = W.cible(0) + 1.5 * np.sqrt(n / 4) * 2 * d / np.linalg.norm(d)
        sens = None if meth == "bq" else S.capteurs(n, 0.0, i, deux_echelles=(meth == "sonde-2e"))
        M = S.MeterLog(W.f, n, "exact" if meth == "bq" else "color", sens, budget, thr=-1.0, log=True)
        suivi = dict(BS.P0, cap=BS.P0["cap"] * max(1.0, np.sqrt(n / 4)), step=1.0, acc=10.0)
        try:
            if meth == "bq": S.broyden_q(M, q0, n, side="L", hfd=0.5, suivi=True)
            else: S.sonde(M, q0, n, modele="H", **suivi)
        except ERRS:
            pass
        runs[meth] = (lab, col, W, M.log)
    W = runs["bq"][2]; j = int(np.argmax(np.linalg.norm((W.cible(0) - W.q0).reshape(n, 4), axis=1)))
    fig = plt.figure(figsize=(11.5, 4.6)); ax = fig.add_subplot(1, 2, 1, projection="3d")
    ts = np.linspace(0, BS.T, 200); C = np.array([W.cible(t).reshape(n, 4)[j] for t in ts])
    ax.plot(C[:, 1], C[:, 2], C[:, 3], color=INK2, lw=3, alpha=0.6, label="zéro q*(t) qui dérive")
    for meth, (lab, col, W_, log) in runs.items():
        if not log: continue
        pts = np.array([q.reshape(n, 4)[j] for c, nf, q in log])
        ax.plot(pts[:, 1], pts[:, 2], pts[:, 3], color=col, lw=1.2, alpha=0.9, label=f"{lab} : {len(log)} lectures")
    ax.set_xlabel("i"); ax.set_ylabel("j"); ax.set_zlabel("k")
    for a in (ax.xaxis, ax.yaxis, ax.zaxis): a.pane.set_facecolor(SURF); a.pane.set_edgecolor(GRID)
    ax.set_title(f"points lus dans (i, j, k) (composante {j + 1} sur {n}) et zéro mobile", fontsize=9.5)
    ax.legend(loc="upper left", fontsize=7, bbox_to_anchor=(-0.08, 1.02)); ax.view_init(elev=24, azim=-58)
    ax2 = fig.add_subplot(1, 2, 2)
    for meth, (lab, col, W_, log) in runs.items():
        if not log: continue
        xs = [c for c, nf, q in log]; ys = [np.linalg.norm(q - W_.cible((c - 1) / r)) for c, nf, q in log]
        ax2.plot(np.array(xs) / r, ys, color=col, lw=1.6, label=lab)
    ax2.set_yscale("log"); ax2.set_xlabel("temps (pas de temps, r = 4 lectures par pas)"); ax2.set_ylabel("distance au zéro du moment")
    ax2.set_title("suivi : la sonde en dent de scie décroche (ambiguïté de décade)", fontsize=9.5); ax2.legend(fontsize=7.5)
    fig.suptitle("Suivre un zéro mobile : même fonction ℍ-linéaire qui dérive, 4 lectures par pas de temps (ℍ⁴)", x=0.02, ha="left", fontsize=12, color=INK)
    fig.tight_layout(rect=(0, 0, 1, 0.94)); os.makedirs(os.path.dirname(path), exist_ok=True); fig.savefig(path); plt.close(fig)


ELEMENTS = {
    "cible": ("2_minimums", "Minimum non nul : sonde à cible scalaire → sonde à cible quaternionique (couleur complète, 36 bits, vallée étroite)",
              ("sonde-couleur-à-cible scalaire", "sonde-couleur-à-cible quaternionique", "Levenberg-Marquardt-sans-couleur"), g_cible),
    "axe": ("1_zeros", "Ajouter l'axe : couleur-sans-axe → couleur-complète, 48 bits (fonction mal conditionnée, ℍ)",
            ("sonde-couleur-sans-axe", "sonde-couleur (avec l'axe)", "Broyden-sans-couleur"), g_axe),
    "memoire": ("1_zeros", "Ajouter la mémoire : relire W à chaque pas → la mettre à jour (Broyden), ℍ",
                ("sonde-couleur", "sonde-couleur-mémoire", "Broyden-sans-couleur"), g_memoire),
    "conforme": ("1_zeros", "Ajouter la structure conforme : 4 sondes → 3 sondes + L(α)R(β) (Möbius, 24 bits)",
                 ("sonde-couleur", "sonde-couleur-conforme-3-sondes", "Broyden-sans-couleur"), g_conforme),
    "dithering": ("1_zeros", "Ajouter le dithering au capteur (24 bits, sonde-couleur-mémoire, fonction mal conditionnée)",
                  ("sans dithering", "avec dithering", "Broyden-sans-couleur"), g_dither),
    "quaternionique": ("3_dimension_4n", "Ajouter la structure ℍ-linéaire : 4n sondes → n sondes (ℍ⁴, n = 4)",
                       ("sonde-couleur-mémoire", "sonde-couleur-mémoire-quaternionique", "Broyden-sans-couleur"), g_quaternionique),
    "depliage": ("3_dimension_4n", "Déplier la luminosité vers zéro au lieu de la prédiction (ℍ¹⁶)",
                 ("sonde-garde-ambiguïté", "sonde-dépliage-vers-zéro", "Broyden-sans-couleur"), g_depliage),
    "relecture": ("3_dimension_4n", "Ajouter la relecture à mi-pas et le double dépliage (ℍ¹⁶, bruit 0,02)",
                  ("sonde-garde-ambiguïté", "sonde-relecture", "Broyden-sans-couleur"), g_relecture),
    "deux_echelles": ("3_dimension_4n", "Remplacer la relecture par un canal de luminosité à deux échelles (ℍ¹⁶, bruit 0,02, 24 bits)",
                      ("sonde-relecture", "sonde-deux-échelles", "Broyden-sans-couleur"), g_deux_echelles),
    "test_geo": ("3_dimension_4n", "Changer le test de descente : moyenne des ‖ρᵢ‖² → moyenne géométrique (ℍ¹⁶, bruit 0,02)",
                 ("sonde-deux-échelles", "sonde-deux-échelles + test géométrique", "Broyden-sans-couleur"), g_test_geo),
    "identite": ("4_autres_structures", "Changer le départ : jacobienne lue en n sondes → a priori « J ≈ identité » (réservoir, ℍ⁸)",
                 ("sonde-couleur-mémoire-quaternionique", "…-départ-identité", "Broyden-sans-couleur"), g_identite),
    "blocconf": ("4_autres_structures", "Étendre la structure : ℍ-linéaire → conforme par blocs (famille bilatère, ℍ³)",
                 ("sonde-couleur-mémoire-quaternionique", "sonde-…-conforme-par-blocs", "Broyden-sans-couleur"), g_blocconf),
    "rang": ("4_autres_structures", "Étendre la structure : ℍ-linéaire → ℍ-linéaire + rang 2 (ℍ⁴)",
             ("sonde-couleur-mémoire-quaternionique", "sonde-…-quaternionique+rang-2", "Broyden-sans-couleur"), g_rang),
    "grise": ("4_autres_structures", "Étendre la structure : ℍ-linéaire → boîte grise L(V)·Δ·L(U) (inversion d'un réseau, ℍ³)",
              ("sonde-couleur-mémoire-quaternionique", "sonde-…-boîte-grise", "Broyden-sans-couleur"), g_grise),
    "par_sortie": ("6_reseaux", "Entraîner une couche : structure à droite seule → + diagonale par sortie (m = 4 sorties, p = 4, ℍ¹⁶)",
                   ("sonde-…-à-droite", "sonde-…-à-droite-par-sortie", "Broyden-sans-couleur"), g_par_sortie),
    "activation": ("6_reseaux", "Changer l'activation du réseau à inverser : tanh par composante → saturation globale (ℍ⁴)",
                   ("tanh par composante", "saturation globale", "Broyden-sans-couleur (saturation globale)"), g_activation),
    "sature": ("6_reseaux", "Régime saturé : on retire un défaut de la sonde (neurone tanh radial, ℍ⁸, décade donnée)",
               ("sonde-oracle-décade (réglage de la partie)", "meilleur réglage du diagnostic", "Broyden-quaternionique-sans-couleur"), g_sature),
    "broyden_q": ("5_armes_egales", "La structure sans la couleur : Broyden-sans-couleur → Broyden-quaternionique-sans-couleur (ℍ⁴, n = 4)",
                  ("Broyden-sans-couleur", "Broyden-quaternionique-sans-couleur", "sonde-couleur-mémoire-quaternionique"), g_broyden_q),
    "instable": ("5_armes_egales", "Conditionnement 10³ : le même capteur à 24 bits lu comme une couleur, puis comme une valeur (ℍ⁴)",
                 ("sonde-oracle-décade (couleur)", "Broyden-quaternionique-mesure-absolue-à-bits-égaux", "Broyden-quaternionique-sans-couleur"), g_instable),
    "amplitude": ("5_armes_egales", "Gain aléatoire à chaque lecture (σℓ = 0,2) : quotients de couleurs → valeur absolue bruitée (ℍ⁴)",
                  ("sonde-deux-échelles + test géométrique", "Broyden-quaternionique-mêmes-mesures", "Broyden-quaternionique-sans-couleur"), g_amplitude),
    "droite": ("6_reseaux", "Entraîner un neurone : quotient à gauche → quotient à droite (ℍ⁴, tanh radial faible)",
               ("sonde-…-quaternionique (à gauche)", "sonde-…-quaternionique-à-droite", "Broyden-sans-couleur"), g_droite),
}

if __name__ == "__main__":
    noms = sys.argv[1:] or [*ELEMENTS, "suivi"]
    pj = os.path.join(HERE, "ijk_resume.json")
    if os.path.exists(pj):
        RESUME.update(json.load(open(pj)))
    for nom in noms:
        if nom == "suivi":
            figure_suivi(os.path.join(PARTIES, "5_armes_egales", "figures", "ijk_suivi.png")); print("suivi: ok"); continue
        partie, titre, labels, gen = ELEMENTS[nom]
        element(nom, partie, titre, labels, gen)
        json.dump(RESUME, open(pj, "w"), indent=1, ensure_ascii=False)
