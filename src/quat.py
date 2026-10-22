"""Algèbre des quaternions et structure conforme de R^4.

Convention : q = w + x i + y j + z k, stocké comme tableau (..., 4) = (w, x, y, z).

Deux faits servent partout dans le projet :

1. Multiplication à gauche / à droite. Pour des quaternions a, b : a·b = L(a) b = R(b) a, avec L(a), R(b) deux matrices
   4x4 réelles. Une application conforme de R^4 qui conserve l'orientation s'écrit h -> α h β (α, β quaternions),
   c'est-à-dire la matrice L(α) R(β) : 7 paramètres réels (8 moins l'échelle commune), contre 16 pour une matrice quelconque.

2. Les 16 matrices B[k,l] = L(e_k) R(e_l) sont orthogonales deux à deux (norme de Frobenius² = 4) : elles forment une base de
   M_4(R). Toute matrice M s'écrit M = Σ Z[k,l] B[k,l], et M est conforme (et directe) si et seulement si Z est de RANG 1
   (Z = α β^T). La matrice conforme la plus proche de M (en Frobenius) s'obtient donc en tronquant le SVD de Z à son premier
   terme : c'est la correspondance de Van Loan - Pitsianis (produit de Kronecker le plus proche) appliquée à ℍ ⊗ ℍ ≅ M_4(R).
   Le « défaut de conformité » sqrt(1 - σ1²/Σσ²) est exactement la distance relative de M aux matrices conformes.
"""
import numpy as np

E4 = np.eye(4)
_CONJ = np.array([1.0, -1.0, -1.0, -1.0])


def qmul(a, b):
    """Produit de Hamilton (l'ordre compte : i j = k, j i = -k)."""
    a = np.asarray(a, float); b = np.asarray(b, float)
    a0, a1, a2, a3 = a[..., 0], a[..., 1], a[..., 2], a[..., 3]
    b0, b1, b2, b3 = b[..., 0], b[..., 1], b[..., 2], b[..., 3]
    return np.stack([a0 * b0 - a1 * b1 - a2 * b2 - a3 * b3,
                     a0 * b1 + a1 * b0 + a2 * b3 - a3 * b2,
                     a0 * b2 - a1 * b3 + a2 * b0 + a3 * b1,
                     a0 * b3 + a1 * b2 - a2 * b1 + a3 * b0], -1)


def qconj(a):
    return np.asarray(a, float) * _CONJ


def qn2(a):
    a = np.asarray(a, float)
    return np.sum(a * a, -1)


def qabs(a):
    return np.sqrt(qn2(a))


def qinv(a):
    a = np.asarray(a, float)
    return qconj(a) / qn2(a)[..., None]


def qexp(v):
    """exp d'un quaternion v = (s, vec) : e^s (cos|vec| + vec/|vec| sin|vec|)."""
    v = np.asarray(v, float)
    s, w = v[..., 0], v[..., 1:]
    nw = np.sqrt(np.sum(w * w, -1))
    sinc = np.where(nw > 1e-12, np.sin(nw) / np.maximum(nw, 1e-300), 1.0 - nw ** 2 / 6)
    return np.exp(s)[..., None] * np.concatenate([np.cos(nw)[..., None], w * sinc[..., None]], -1)


def Lmat(a):
    """Matrice 4x4 de h -> a h."""
    a0, a1, a2, a3 = a
    return np.array([[a0, -a1, -a2, -a3], [a1, a0, -a3, a2], [a2, a3, a0, -a1], [a3, -a2, a1, a0]])


def Rmat(b):
    """Matrice 4x4 de h -> h b."""
    b0, b1, b2, b3 = b
    return np.array([[b0, -b1, -b2, -b3], [b1, b0, b3, -b2], [b2, -b3, b0, b1], [b3, b2, -b1, b0]])


# base orthogonale de M_4(R) : B[k, l] = L(e_k) R(e_l)
BASIS = np.array([[Lmat(E4[k]) @ Rmat(E4[l]) for l in range(4)] for k in range(4)])


def to_Z(M):
    """Coordonnées de M dans la base B : M = Σ Z[k,l] B[k,l]. M conforme directe <=> Z de rang 1."""
    return 0.25 * np.einsum("klij,ij->kl", BASIS, M)


def from_Z(Z):
    return np.einsum("kl,klij->ij", Z, BASIS)


def conf_project(M):
    """Matrice conforme directe la plus proche de M (Frobenius)."""
    U, S, Vt = np.linalg.svd(to_Z(M))
    return from_Z(S[0] * np.outer(U[:, 0], Vt[0]))


def conf_factors(M):
    """(α, β) tels que L(α) R(β) soit la matrice conforme la plus proche de M (|β| = 1)."""
    U, S, Vt = np.linalg.svd(to_Z(M))
    return np.sqrt(S[0]) * U[:, 0], np.sqrt(S[0]) * Vt[0]


def conf_defect(M):
    """Distance relative de M aux matrices conformes directes, dans [0, 1[ : 0 = exactement conforme."""
    S = np.linalg.svd(to_Z(M), compute_uv=False)
    return float(np.sqrt(max(0.0, 1.0 - S[0] ** 2 / np.sum(S ** 2))))


def cross3(a, b, c):
    """Vecteur de R^4 orthogonal à a, b, c (produit vectoriel généralisé), de norme = volume du parallélépipède.
    Le signe est choisi pour que det[a b c cross] > 0."""
    M = np.stack([a, b, c], 1)                      # 4x3
    out = np.empty(4)
    for i in range(4):
        out[i] = (-1) ** (i + 1) * np.linalg.det(np.delete(M, i, axis=0))
    return out


def mobius(a, b, c, d):
    """q -> (a q + b)(c q + d)^{-1}. L'ordre des produits compte."""
    a, b, c, d = map(lambda x: np.asarray(x, float), (a, b, c, d))
    def f(q):
        q = np.asarray(q, float)
        return qmul(qmul(a, q) + b, qinv(qmul(c, q) + d))
    return f


def jac_num(f, q, h=1e-6):
    """Jacobienne 4x4 par différences centrées (usage interne, non facturé)."""
    q = np.asarray(q, float)
    J = np.empty((4, 4))
    for j in range(4):
        e = E4[j] * h
        J[:, j] = (f(q + e) - f(q - e)) / (2 * h)
    return J


def rand_unit(rng, n=4):
    v = rng.standard_normal(n)
    return v / np.linalg.norm(v)
