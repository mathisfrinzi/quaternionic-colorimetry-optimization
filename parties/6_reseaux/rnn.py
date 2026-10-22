"""Partie 6 : réseaux quaternioniques (idées 1.1 à 1.3). Objets communs.

RNN quaternionique à entrée constante :  h ← σ(W·h + U·x + b),  h ∈ ℍⁿ, x ∈ ℍᵐ ; ses points fixes sont les zéros de
    F(h) = h − σ(W·h + U·x + b),        J = I − Dσ(z)·W.
W = W₀ + P : W₀ quaternionique (ℍ-linéaire) GELÉ, P la PIÈCE APPRISE (esprit du dépôt half-frozen-rnn). Trois façons d'apprendre P :
    "libre"      P réelle 4n×4n quelconque (moindres carrés sans contrainte)
    "H"          P ℍ-linéaire (blocs L(pᵢⱼ)) : la structure est garantie
    "penalise"   P réelle, avec une pénalité λ‖P − proj_H(P)‖² sur sa partie non ℍ-linéaire
Tâche apprise : MÉMOIRE ASSOCIATIVE. Pour K entrées xₖ, on veut que le point fixe soit une cible h*ₖ donnée, i.e.
    (W₀ + P)·h*ₖ = σ⁻¹(h*ₖ) − U·xₖ − b,  k = 1…K
qui est LINÉAIRE en P : on la résout exactement par moindres carrés (régularisés). C'est un apprentissage réel (P change ce que
le réseau calcule), sans rétropropagation.

Activations (comme en partie 4.2) : "split" (tanh par composante), "radial" (tanh(‖z‖) z/‖z‖ par quaternion), "global".
"""
import os, sys
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "src"))
sys.path.insert(0, os.path.join(HERE, "..", "4_autres_structures", "4.2_reservoirs"))
from quat import Lmat, Rmat, qmul, rand_unit, E4
from reservoir import act, dact_block, dact_global, hdefect
from sonde_hn import Lblock, proj_H

LB = [Lmat(E4[k]) for k in range(4)]


def act_inv(y, kind):
    """σ⁻¹ (cibles choisies dans l'image de σ)."""
    if kind == "split":
        return np.arctanh(np.clip(y, -0.999, 0.999))
    if kind == "global":
        n = len(y) // 4; r = np.linalg.norm(y) / np.sqrt(n)
        return y * (np.arctanh(min(r, 0.999)) / r if r > 1e-12 else 1.0)
    z = y.reshape(-1, 4); r = np.linalg.norm(z, axis=1, keepdims=True)
    s = np.where(r > 1e-12, np.arctanh(np.minimum(r, 0.999)) / np.maximum(r, 1e-300), 1.0)
    return (s * z).reshape(-1)


def qmat(rng, n, m=None, scale=1.0):
    m = n if m is None else m
    return scale * rng.standard_normal((n, m, 4)) / np.sqrt(4 * m)


def Lrect(A):
    n, m = A.shape[:2]; M = np.zeros((4 * n, 4 * m))
    for i in range(n):
        for j in range(m):
            M[4 * i:4 * i + 4, 4 * j:4 * j + 4] = Lmat(A[i, j])
    return M


class RNN:
    def __init__(self, n, kind="radial", rho0=0.6, K=None, piece="H", lam=10.0, seed=0, m=2, beta=0.8, cible=0.5):
        rng = np.random.default_rng(61000 + 1013 * seed + 7 * n)
        self.n, self.kind, self.m = n, kind, m
        W0 = Lblock(qmat(rng, n)); W0 *= rho0 / max(abs(np.linalg.eigvals(W0)))
        self.W0 = W0
        self.LU = Lrect(qmat(rng, n, m, 1.5))
        self.b = 0.3 * rng.standard_normal(4 * n)
        K = K if K is not None else max(2, n // 2)
        self.X = [beta * rng.standard_normal(4 * m) for _ in range(K)]
        # cibles : le point fixe du réseau gelé pour l'entrée xₖ, déplacé de `cible` par quaternion (ce que P doit apprendre)
        self.W = W0; self.P = np.zeros_like(W0)
        self.H = []
        for k in range(K):
            h0 = self.point_fixe(k)
            if h0 is None: h0 = np.zeros(4 * n)
            self.H.append(np.clip(h0 + cible * rng.standard_normal(4 * n) / 2, -0.9, 0.9))
        self.piece = piece
        self.P = self._apprendre(piece, lam) if piece != "aucune" else np.zeros_like(W0)
        self.W = W0 + self.P

    def _apprendre(self, piece, lam):
        n = self.n
        Hs = np.column_stack(self.H)                                                 # (4n, K)
        T = np.column_stack([act_inv(h, self.kind) - self.LU @ x - self.b for h, x in zip(self.H, self.X)]) - self.W0 @ Hs
        ridge = 1e-3
        if piece == "libre":
            return T @ Hs.T @ np.linalg.inv(Hs @ Hs.T + ridge * np.eye(4 * n))
        # base de l'espace des P : blocs L(e_k) (ℍ-linéaires, 4n²) puis le complément (12n²) ; P = Σ c_a B_a
        # résolu dans la base canonique de ℝ^{16n²} avec une pénalité sur la partie non ℍ-linéaire
        d = 4 * n
        # opérateur P ↦ P·Hs, vectorisé (vec en lignes) : (P Hs)_{r,k} = Σ_c P_{r,c} Hs_{c,k}
        A = np.kron(np.eye(d), Hs.T)                                                 # (d·K, d²)
        y = T.reshape(-1)
        Pi = self._proj_op(n)                                                        # projection orthogonale sur ℍ-lin (d², d²)
        if piece == "H":
            Bq = self._base_H(n)                                                     # (d², 4n²)
            c = np.linalg.solve((A @ Bq).T @ (A @ Bq) + ridge * np.eye(Bq.shape[1]), (A @ Bq).T @ y)
            return (Bq @ c).reshape(d, d)
        Q = np.eye(d * d) - Pi
        p = np.linalg.solve(A.T @ A + lam * Q + ridge * np.eye(d * d), A.T @ y)
        return p.reshape(d, d)

    @staticmethod
    def _base_H(n):
        d = 4 * n; cols = []
        for i in range(n):
            for j in range(n):
                for k in range(4):
                    M = np.zeros((d, d)); M[4 * i:4 * i + 4, 4 * j:4 * j + 4] = LB[k] / 2.0     # norme 1
                    cols.append(M.reshape(-1))
        return np.column_stack(cols)

    @classmethod
    def _proj_op(cls, n):
        B = cls._base_H(n)
        return B @ B.T

    # ------------------------------------------------------------------------------------------- dynamique
    def F(self, h, k=0):
        h = np.asarray(h, float)
        return h - act(self.W @ h + self.LU @ self.X[k] + self.b, self.kind)

    def J(self, h, k=0):
        z = self.W @ h + self.LU @ self.X[k] + self.b
        if self.kind == "global":
            return np.eye(4 * self.n) - dact_global(z) @ self.W
        zz = z.reshape(-1, 4); D = np.zeros((4 * self.n, 4 * self.n))
        for i in range(self.n):
            D[4 * i:4 * i + 4, 4 * i:4 * i + 4] = dact_block(zz[i], self.kind)
        return np.eye(4 * self.n) - D @ self.W

    def point_fixe(self, k=0, h0=None):
        h = np.zeros(4 * self.n) if h0 is None else np.array(h0, float)
        if h0 is None:
            for _ in range(200):
                h = act(self.W @ h + self.LU @ self.X[k] + self.b, self.kind)
        for _ in range(50):
            r = self.F(h, k)
            if np.linalg.norm(r) < 1e-13: break
            h = h - np.linalg.solve(self.J(h, k), r)
        return h if np.linalg.norm(self.F(h, k)) < 1e-9 else None

    def defauts(self):
        """Défaut ℍ-linéaire de J au point fixe de chaque entrée, et défaut de la pièce apprise ‖P − proj_H P‖/‖P‖."""
        dj = []
        for k in range(len(self.X)):
            hs = self.point_fixe(k, self.H[k])
            if hs is not None:
                dj.append(hdefect(self.J(hs, k)))
        nP = np.linalg.norm(self.P)
        dP = 0.0 if nP == 0 else float(np.linalg.norm(self.P - Lblock(proj_H(self.P))) / nP)
        return float(np.mean(dj)) if dj else np.nan, dP

    def erreur_memoire(self):
        """Erreur relative moyenne entre le point fixe atteint et la cible h*ₖ (qualité de l'apprentissage)."""
        e = []
        for k in range(len(self.X)):
            hs = self.point_fixe(k, self.H[k])
            e.append(np.nan if hs is None else np.linalg.norm(hs - self.H[k]) / np.linalg.norm(self.H[k]))
        return float(np.nanmean(e))


# ================================================================================================ inversion d'un réseau (6.2)
def qmat41(rng, n, diag=1.0, off=0.35):
    """Même tirage que bench_hors_structure.qmat (partie 4.1)."""
    A = off / np.sqrt(n) * rng.standard_normal((n, n, 4))
    for i in range(n):
        A[i, i] += diag * rand_unit(rng)
    return A


class QNN:
    """Réseau à une couche cachée y = V·σ(U·x + b) ; on l'inverse : f(x) = V·σ(U·x + b) − y*, zéro en x*.
    Même tirage que la famille « inversion réseau » de la partie 4.1 (graine 9000 + 131·seed + 17·n + 7·3), activation au choix."""
    def __init__(self, n, kind="split", seed=0):
        rng = np.random.default_rng(9000 + 131 * seed + 17 * n + 7 * 3)
        self.qs = 0.5 * rng.standard_normal(4 * n)
        self.U = qmat41(rng, n, 1.0, 0.6); self.V = qmat41(rng, n, 1.0, 0.6); self.b = 0.8 * rng.standard_normal(4 * n)
        self.LU, self.LV, self.n, self.kind = Lblock(self.U), Lblock(self.V), n, kind
        self.ys = self.net(self.qs)

    def net(self, x):
        return self.LV @ act(self.LU @ np.asarray(x, float) + self.b, self.kind)

    def f(self, x):
        return self.net(x) - self.ys

    def J(self, x):
        z = self.LU @ x + self.b
        if self.kind == "global":
            D = dact_global(z)
        else:
            zz = z.reshape(-1, 4); D = np.zeros((4 * self.n, 4 * self.n))
            for i in range(self.n):
                D[4 * i:4 * i + 4, 4 * i:4 * i + 4] = dact_block(zz[i], self.kind)
        return self.LV @ D @ self.LU
