"""Réservoirs quaternioniques et défaut ℍ-linéaire (partie 10).

Un réservoir quaternionique à entrée constante a pour points fixes les zéros de  F(h) = h − σ(W·h + b),  h ∈ ℍⁿ,
où W est une matrice n×n de quaternions (h ↦ W·h agit par multiplication à GAUCHE : ℍ-linéaire) et b = U·x le terme d'entrée.
Sa jacobienne vaut  J = I − Dσ(z)·L(W),  z = W·h + b,  Dσ = blocs 4×4 de la dérivée de l'activation.
J est ℍ-linéaire dès que chaque bloc de Dσ est une multiplication à gauche ; le défaut vient donc tout entier de l'activation
(ou d'une partie réelle non quaternionique ajoutée à W, ici `t_real`, qui imite une pièce apprise sans contrainte).

Activations :
  split  : tanh appliquée à chaque composante réelle (réseaux quaternioniques usuels)
  radial : σ(z) = tanh(|z|)·z/|z|  (ne touche que le module de chaque quaternion)
  global : σ(z) = tanh(ρ)·z/ρ, ρ = ‖z‖/√n la norme moyenne par quaternion (une seule saturation pour tout le vecteur, comme une
           normalisation de couche) : la déviation à l'identité est de rang 1 sur tout ℍⁿ, son poids relatif décroît en 1/√(4n)
"""
import os, sys
import numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "..", "src"))
from quat import Lmat, rand_unit

LB = [Lmat(np.eye(4)[k]) for k in range(4)]


def act(z, kind):
    if kind == "global":
        n = len(z) // 4; r = np.linalg.norm(z) / np.sqrt(n)
        return z * (np.tanh(r) / r if r > 1e-12 else 1.0)
    z = z.reshape(-1, 4)
    if kind == "split":
        return np.tanh(z).reshape(-1)
    r = np.linalg.norm(z, axis=1, keepdims=True)
    s = np.where(r > 1e-12, np.tanh(r) / np.maximum(r, 1e-300), 1.0)
    return (s * z).reshape(-1)


def dact_block(z, kind):
    """Dérivée 4×4 de l'activation en un quaternion z."""
    if kind == "split":
        return np.diag(1 - np.tanh(z) ** 2)
    r = np.linalg.norm(z)
    if r < 1e-8:
        return np.eye(4)
    s = np.tanh(r) / r
    ds = (1 - np.tanh(r) ** 2) / r - np.tanh(r) / r ** 2       # dérivée de s par rapport à r
    u = z / r
    return s * np.eye(4) + r * ds * np.outer(u, u)


def dact_global(z):
    n = len(z) // 4; r = np.linalg.norm(z) / np.sqrt(n)
    if r < 1e-8:
        return np.eye(len(z))
    s = np.tanh(r) / r; ds = (1 - np.tanh(r) ** 2) / r - np.tanh(r) / r ** 2
    return s * np.eye(len(z)) + ds * np.outer(z, z) / (n * r)        # dr/dz = z / (n r)


def block_defect(B):
    """Défaut ℍ-linéaire d'un bloc 4×4 : distance relative aux multiplications à gauche."""
    P = sum(np.sum(B * L) / 4 * L for L in LB)
    return float(np.linalg.norm(B - P) / np.linalg.norm(B))


def hdefect(J):
    n = J.shape[0] // 4
    P = np.zeros_like(J)
    for i in range(n):
        for j in range(n):
            B = J[4 * i:4 * i + 4, 4 * j:4 * j + 4]
            P[4 * i:4 * i + 4, 4 * j:4 * j + 4] = sum(np.sum(B * L) / 4 * L for L in LB)
    return float(np.linalg.norm(J - P) / np.linalg.norm(J))


class Reservoir:
    def __init__(self, n, rho=0.8, beta=1.0, kind="split", seed=0, t_real=0.0):
        rng = np.random.default_rng(31000 + 1009 * seed + 17 * n)
        self.n, self.kind = n, kind
        Wq = rng.standard_normal((n, n, 4)) / np.sqrt(4 * n)
        LW = np.zeros((4 * n, 4 * n))
        for i in range(n):
            for j in range(n):
                LW[4 * i:4 * i + 4, 4 * j:4 * j + 4] = Lmat(Wq[i, j])
        LW *= rho / max(abs(np.linalg.eigvals(LW)))
        if t_real > 0:                                          # pièce réelle (non quaternionique), ‖P‖ = t ‖W‖
            P = rng.standard_normal((4 * n, 4 * n)); P *= t_real * np.linalg.norm(LW) / np.linalg.norm(P)
            LW = LW + P
            LW *= rho / max(abs(np.linalg.eigvals(LW)))
        self.LW = LW
        bq = np.array([rand_unit(rng) for _ in range(n)]) * np.abs(rng.normal(1, 0.3, (n, 1)))
        self.b = beta * bq.reshape(-1)
        self.hs = self._fixed_point()

    def F(self, h):
        return np.asarray(h, float) - act(self.LW @ np.asarray(h, float) + self.b, self.kind)

    def J(self, h):
        if self.kind == "global":
            return np.eye(4 * self.n) - dact_global(self.LW @ h + self.b) @ self.LW
        z = (self.LW @ h + self.b).reshape(-1, 4)
        D = np.zeros((4 * self.n, 4 * self.n))
        for i in range(self.n):
            D[4 * i:4 * i + 4, 4 * i:4 * i + 4] = dact_block(z[i], self.kind)
        return np.eye(4 * self.n) - D @ self.LW

    def _fixed_point(self):
        h = np.zeros(4 * self.n)
        for _ in range(3000):
            h = act(self.LW @ h + self.b, self.kind)
        for _ in range(30):                                     # Newton exact pour finir
            r = self.F(h)
            if np.linalg.norm(r) < 1e-13: break
            h = h - np.linalg.solve(self.J(h), r)
        assert np.linalg.norm(self.F(h)) < 1e-9, "point fixe non trouvé"
        return h

    def defect(self):
        return hdefect(self.J(self.hs))
