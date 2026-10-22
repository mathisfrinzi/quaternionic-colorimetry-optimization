"""Schémas pédagogiques et représentations dans l'espace des quaternions. Usage : python schemas/make_schemas.py"""
import os, sys
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "src"))
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch, Circle, Arc
from report import PAL, style, INK, INK2, MUTED, GRID, SURF
from quat import qmul, qconj, rand_unit, to_Z, Lmat, Rmat
from sensors import Sensor
from funcs import build
from probes import Meter, Hit, OutOfBudget, probe_secant
import classic
style()
OUT = HERE


def axis_rgb(u): return (np.asarray(u) + 1) / 2


# ------------------------------------------------------------------ 1. couleur d'un quaternion
def fig_couleur():
    fig = plt.figure(figsize=(12, 3.9))
    ax = fig.add_subplot(1, 3, 1, projection="3d")
    n = 60; ph, la = np.meshgrid(np.linspace(0, np.pi, n), np.linspace(0, 2 * np.pi, 2 * n))
    X, Y, Z = np.sin(ph) * np.cos(la), np.sin(ph) * np.sin(la), np.cos(ph)
    cols = np.stack([axis_rgb(np.array([x, y, z])) for x, y, z in zip(X.ravel(), Y.ravel(), Z.ravel())]).reshape(X.shape + (3,))
    ax.plot_surface(X, Y, Z, facecolors=cols, rstride=1, cstride=1, shade=False, linewidth=0, antialiased=False)
    u = np.array([0.6, 0.3, 0.74]); u /= np.linalg.norm(u)
    ax.quiver(*(1.0 * u), *(0.9 * u), color=INK, lw=3, arrow_length_ratio=0.25)
    ax.text(*(2.05 * u), "u", fontsize=14)
    ax.view_init(elev=22, azim=25)
    ax.set_title("l'axe u ∈ S² → couleur RGB (u+1)/2\n(ou deux angles)"); ax.set_axis_off()
    ax = fig.add_subplot(1, 3, 2)
    th, r = 0.9, 1.0
    ax.add_patch(Arc((0, 0), 0.9, 0.9, theta1=0, theta2=np.degrees(th), color=PAL["orange"], lw=2))
    ax.annotate("", xy=(r * np.cos(th), r * np.sin(th)), xytext=(0, 0), arrowprops=dict(arrowstyle="-|>", color=INK, lw=2))
    ax.plot([0, 1.35], [0, 0], color=MUTED, lw=1); ax.plot([0, 0], [0, 1.3], color=MUTED, lw=1)
    ax.text(1.38, -0.05, "Re F", color=INK2); ax.text(0.03, 1.32, "‖Im F‖ (dans la direction u)", color=INK2)
    ax.text(0.55, 0.13, "θ", color=PAL["orange"], fontsize=13); ax.text(0.18, 0.62, "r = |F|", color=INK)
    ax.set_xlim(-0.2, 1.9); ax.set_ylim(-0.2, 1.5); ax.set_aspect("equal"); ax.set_axis_off()
    ax.set_title(r"$F = r\,e^{u\theta} = r\,(\cos\theta + u\,\sin\theta)$" + "\nluminosité ↔ r, teinte ↔ θ, axe ↔ u")
    ax = fig.add_subplot(1, 3, 3)
    r = np.logspace(-1, 2.2, 600); ax.semilogx(r, np.mod(np.log10(r), 1.0), color=PAL["blue"], lw=2)
    ax.set_xlabel("|F| (échelle log)"); ax.set_ylabel("luminosité lue = log₁₀|F| mod 1")
    ax.set_title("la luminosité est une dent de scie :\ntous les bits servent à la précision fine")
    fig.tight_layout(); fig.savefig(os.path.join(OUT, "01_couleur_quaternion.png")); plt.close(fig)


# ------------------------------------------------------------------ 2. pipeline d'une sonde
def box(ax, x, y, w, h, text, fc):
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.02,rounding_size=0.05", fc=fc, ec=INK2, lw=1))
    ax.text(x + w / 2, y + h / 2, text, ha="center", va="center", fontsize=9, color=INK)


def arrow(ax, x0, y0, x1, y1, text=None):
    ax.add_patch(FancyArrowPatch((x0, y0), (x1, y1), arrowstyle="-|>", mutation_scale=14, color=INK2, lw=1.4))
    if text: ax.text((x0 + x1) / 2, (y0 + y1) / 2 - 0.22, text, ha="center", fontsize=8, color=INK2)


def fig_pipeline():
    fig, ax = plt.subplots(figsize=(12, 4.3)); ax.set_xlim(0, 12); ax.set_ylim(0, 4.3); ax.set_axis_off()
    box(ax, 0.1, 2.9, 1.9, 0.9, "point q\n(4 réels)", "#e9f0fb")
    box(ax, 2.6, 2.9, 1.9, 0.9, "fonction f\n(boîte noire)", "#f1f1ee")
    box(ax, 5.1, 2.9, 2.2, 0.9, "capteur S2\n(quantifie, bruite)", "#fdeee6")
    box(ax, 7.9, 2.9, 2.0, 0.9, "couleur lue\n(ℓ, θ, u) en B bits", "#e6f6ef")
    arrow(ax, 2.0, 3.35, 2.6, 3.35); arrow(ax, 4.5, 3.35, 5.1, 3.35); arrow(ax, 7.3, 3.35, 7.9, 3.35)
    ax.text(10.1, 3.35, "l'algorithme\nne voit QUE ceci", fontsize=9, color=PAL["red"], va="center")
    box(ax, 0.1, 1.2, 2.8, 1.0, "quotient   ρ = F(q)⁻¹ F(q+h)\n(deux couleurs → un quaternion)", "#e6f6ef")
    box(ax, 3.4, 1.2, 3.0, 1.0, "jacobienne normalisée\nW = F⁻¹J ≈ (ρ − 1)/η", "#e6f6ef")
    box(ax, 6.9, 1.2, 2.6, 1.0, "pas de Newton\nW h = −1  →  h = −W⁺ e₀", "#e6f6ef")
    box(ax, 9.9, 1.2, 1.9, 1.0, "test lu :\n|ρ_nouveau| < 0,92", "#e6f6ef")
    arrow(ax, 2.9, 1.7, 3.4, 1.7); arrow(ax, 6.4, 1.7, 6.9, 1.7); arrow(ax, 9.5, 1.7, 9.9, 1.7)
    arrow(ax, 8.9, 2.9, 1.5, 2.2, "deux lectures : en q puis en q + η eⱼ")
    ax.text(0.1, 0.55, "Les trois familles de sondes ne diffèrent que par la façon d'obtenir W :", fontsize=9.5, color=INK)
    ax.text(0.1, 0.2, "FD : 4 sondes (η e₀…η e₃) à chaque pas   ·   SÉCANTE : 1 lecture par pas, W mise à jour   ·   CONFORME : 3 sondes, W = L(α)R(β)", fontsize=9, color=INK2)
    fig.suptitle("Ce que fait une sonde, étape par étape", x=0.01, ha="left", fontsize=12)
    fig.savefig(os.path.join(OUT, "02_pipeline_sonde.png")); plt.close(fig)


# ------------------------------------------------------------------ 3. S³ colorée (projection stéréographique)
def stereo(p):
    return p[:, 1:] / (1 - p[:, :1] + 1e-9)


def fig_s3():
    """Coloration de domaine sur S³ : on paramètre S³ par les tores de Hopf q = (cos η e^{iξ₁}) + (sin η e^{iξ₂}) j
    (η fixé = un tore plat, ξ₁, ξ₂ ∈ [0, 2π[) et on peint chaque point de la couleur de f(q)/|f(q)| (axe → RGB, θ → luminosité),
    comme la coloration de domaine du plan complexe."""
    n = 160
    x1, x2 = np.meshgrid(np.linspace(0, 2 * np.pi, n), np.linspace(0, 2 * np.pi, n))
    fns = [("f(q) = q  (degré +1)", lambda a: a), ("f(q) = q̄  (degré −1)", lambda a: qconj(a)),
           ("f(q) = q²  (degré +2)", lambda a: qmul(a, a)), ("f(q) = q̄²  (degré −2)", lambda a: qmul(qconj(a), qconj(a)))]
    etas = [np.pi / 8, np.pi / 4, 3 * np.pi / 8]
    fig, ax = plt.subplots(len(etas), 4, figsize=(12, 8.6))
    for r, eta in enumerate(etas):
        Q = np.stack([np.cos(eta) * np.cos(x1), np.cos(eta) * np.sin(x1), np.sin(eta) * np.cos(x2), np.sin(eta) * np.sin(x2)], -1)
        for k, (t, f) in enumerate(fns):
            F = np.apply_along_axis(f, -1, Q); F /= np.linalg.norm(F, axis=-1, keepdims=True)
            u = F[..., 1:] / (np.linalg.norm(F[..., 1:], axis=-1, keepdims=True) + 1e-12)
            th = np.arccos(np.clip(F[..., 0], -1, 1))
            img = axis_rgb(u) * (0.25 + 0.75 * (th / np.pi))[..., None]
            ax[r, k].imshow(np.clip(img, 0, 1), origin="lower", extent=(0, 2 * np.pi, 0, 2 * np.pi))
            ax[r, k].set_xticks([0, np.pi, 2 * np.pi]); ax[r, k].set_xticklabels(["0", "π", "2π"])
            ax[r, k].set_yticks([0, np.pi, 2 * np.pi]); ax[r, k].set_yticklabels(["0", "π", "2π"]); ax[r, k].grid(False)
            if r == 0: ax[r, k].set_title(t, fontsize=10)
            if k == 0: ax[r, k].set_ylabel(f"η = {['π/8', 'π/4', '3π/8'][r]}\nξ₂")
            if r == len(etas) - 1: ax[r, k].set_xlabel("ξ₁")
    fig.suptitle("Coloration de domaine sur la sphère S³ : trois tores de Hopf, chaque point peint avec la couleur de f/|f|\n"
                 "(teinte = axe u en RGB, sombre = θ petit). Un degré ±d se lit comme d enroulements des couleurs, dans un sens ou dans l'autre.",
                 x=0.01, ha="left", fontsize=10.5)
    fig.tight_layout(); fig.savefig(os.path.join(OUT, "03_S3_couleurs_degre.png")); plt.close(fig)


# ------------------------------------------------------------------ 4. conformité : cercles et ellipses
def fig_conformite():
    rng = np.random.default_rng(3)
    a, b = rand_unit(rng) * 1.3, rand_unit(rng) * 0.8
    Wc = Lmat(a) @ Rmat(b)                              # conforme (Möbius)
    Wg = rng.standard_normal((4, 4)); Wg = Wg / np.linalg.svd(Wg)[1].mean()    # générique
    fig, ax = plt.subplots(1, 3, figsize=(12, 3.8))
    t = np.linspace(0, 2 * np.pi, 200)
    for k, (W, nm, c) in enumerate(((Wc, "conforme  h ↦ α h β  (7 paramètres)", PAL["blue"]), (Wg, "générique  (16 paramètres)", PAL["orange"]))):
        for j in range(5):
            Q = np.linalg.qr(rng.standard_normal((4, 2)))[0]            # plan aléatoire de l'espace d'entrée
            Im = W @ Q; B = np.linalg.qr(Im)[0]; ell = B.T @ Im @ np.vstack([np.cos(t), np.sin(t)])
            ax[k].plot(ell[0], ell[1], color=c, lw=1.6, alpha=0.85)
        ax[k].set_aspect("equal"); ax[k].set_title(nm); ax[k].set_xlim(-2.6, 2.6); ax[k].set_ylim(-2.6, 2.6)
    ax[0].set_xlabel("5 plans de ℍ : chaque cercle unité reste un cercle (superposés)"); ax[1].set_xlabel("… une ellipse")
    for W, nm, c, o in ((Wc, "conforme", PAL["blue"], -0.15), (Wg, "générique", PAL["orange"], 0.15)):
        ax[2].bar(np.arange(4) + o, np.linalg.svd(W)[1], 0.3, color=c, label=nm)
    ax[2].set_title("valeurs singulières de W"); ax[2].legend(); ax[2].set_xticks(range(4)); ax[2].set_xticklabels(["σ₁", "σ₂", "σ₃", "σ₄"])
    fig.suptitle("Pourquoi la structure conforme aide : une jacobienne L(α)R(β) a 4 valeurs singulières égales (rotation + dilatation)", x=0.01, ha="left", fontsize=10.5)
    fig.tight_layout(); fig.savefig(os.path.join(OUT, "04_conformite.png")); plt.close(fig)


# ------------------------------------------------------------------ 5. trajectoires dans ℍ (4D) : sonde contre Broyden
def trace(method, P, q0, sensor=None, **kw):
    pts, vals = [], []
    M = Meter(P.f, sensor, budget=300, thr=P.thr, absolute=(sensor is not None and method == "broyden"))
    orig = M.read
    def rd(q):
        pts.append(np.array(q, float)); vals.append(float(np.linalg.norm(P.f(q)))); return orig(q)
    M.read = rd
    try:
        if method == "sonde": probe_secant(M, q0, model="rank1", **kw)
        else: classic.broyden(M, q0, **kw)
    except (Hit, OutOfBudget):
        pass
    return np.array(pts), np.array(vals)


def fig_trajectoires():
    P = build(("gen", "quad", 0)); q0 = P.starts(1, 4)[0]
    ps, vs = trace("sonde", P, q0, Sensor("s2", 36, "sph", seed=1), c=0.3, step=0.65, cap=1.0)
    pb, vb = trace("broyden", P, q0, Sensor("s2", 36, "sph", seed=1), hfd=1e-1)
    zero = P.zero
    fig = plt.figure(figsize=(12.5, 4))
    for k, (cols, nm) in enumerate((((1, 2, 3), "parties imaginaires (i, j, k)"), ((0, 1, 2), "(1, i, j)"))):
        ax = fig.add_subplot(1, 3, k + 1, projection="3d")
        for pts, c, lab in ((ps, PAL["blue"], "sonde sécante (couleurs)"), (pb, PAL["orange"], "Broyden (mêmes mesures)")):
            ax.plot(*[pts[:, i] for i in cols], "-", color=c, lw=1, alpha=0.8, label=lab); ax.scatter(*[pts[:, i] for i in cols], s=8, color=c)
        ax.scatter(*[[zero[i]] for i in cols], marker="*", s=160, color=PAL["red"], label="zéro de f")
        ax.scatter(*[[q0[i]] for i in cols], marker="s", s=40, color=INK, label="départ")
        ax.set_title(f"trajectoire dans ℍ ≅ ℝ⁴ — {nm}", fontsize=9)
        if k == 0: ax.legend(fontsize=7, loc="upper left")
    ax = fig.add_subplot(1, 3, 3)
    ax.semilogy(np.arange(1, len(vs) + 1), vs, "o-", ms=3, color=PAL["blue"], label=f"sonde ({len(vs)} lectures)")
    ax.semilogy(np.arange(1, len(vb) + 1), vb, "o-", ms=3, color=PAL["orange"], label=f"Broyden ({len(vb)} lectures)")
    ax.axhline(P.thr, color=PAL["red"], lw=1, ls="--"); ax.set_xlabel("lecture"); ax.set_ylabel("|f| (vraie valeur, jamais vue par la sonde)"); ax.legend(fontsize=8)
    fig.tight_layout(); fig.savefig(os.path.join(OUT, "05_trajectoires_4D.png")); plt.close(fig)


if __name__ == "__main__":
    for fn in (fig_couleur, fig_pipeline, fig_s3, fig_conformite, fig_trajectoires):
        fn(); print("ok", fn.__name__)
