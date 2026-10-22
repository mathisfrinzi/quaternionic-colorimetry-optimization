"""Sonde générale en ℍⁿ (dimension réelle 4n) : la sonde-couleur-mémoire-quaternionique des parties 3.1 et 3.4, avec des
modèles de jacobienne et des remèdes à l'ambiguïté interchangeables. Le code des parties 3.1 à 3.4 n'est pas modifié ; ce module
le reproduit (vérifié par `tests/test_sonde_hn.py` : mêmes lectures que `probe34` sur des départs fixés) et l'étend.

Rappel. Une lecture donne n couleurs. Pour chaque composante, le quotient de deux lectures vaut, selon le côté :
  à gauche   ρᵢ = Fᵢ(q)⁻¹ Fᵢ(q + h) = 1 + (W h)ᵢ,  W = diag(F)⁻¹ J                (parties 3 et 4)
  à droite   ρᵢ = Fᵢ(q + h) Fᵢ(q)⁻¹ = 1 + (W h)ᵢ,  (W h)ᵢ = (J h)ᵢ Fᵢ⁻¹             (partie 6.3)
et le pas de Newton résout W h = −𝟙. Le côté du quotient doit être celui qui CONSERVE la structure de J : une jacobienne dont les
blocs sont des multiplications à gauche L(a) reste de cette forme après normalisation à gauche (L(F⁻¹)L(a) = L(F⁻¹a)) ; une
jacobienne à blocs R(c) (multiplication à droite) reste à blocs R après normalisation à droite (h·c·F⁻¹ = R(cF⁻¹)h).

Modèles de W (`modele=`) :
  "H"        blocs L(wᵢⱼ), n sondes ; Broyden quaternionique (partie 3.1)              4n² paramètres
  "D"        blocs R(wᵢⱼ), n sondes, quotient à droite ; Broyden quaternionique à droite : wᵢⱼ += h̄ⱼ·rᵢ/‖h‖²
  "full"     réelle 4n×4n, 4n sondes ; Broyden réel                                      16n²
  "blocconf" blocs conformes L(αᵢⱼ)R(βᵢⱼ), 3n sondes (directions 1, i, j de chaque bloc) ; Broyden réel puis projection de
             chaque bloc sur les conformes (correspondance de Van Loan-Pitsianis de src/quat.py)          7n² (8n² moins l'échelle)
  "Hrang"    L(A) + C avec C réelle de rang ≤ k, n sondes ; A et C ajustés par moindres carrés alternés sur les K dernières
             sécantes (A au plus près de la précédente, C = troncature de rang k)                        4n² + 8nk
  "grise"    L(V)·Δ·L(U) avec U CONNUE (boîte grise : on connaît la première couche d'un réseau), Δ diagonale réelle ;
             2n sondes (directions L(U)⁻¹e₄ⱼ et L(U)⁻¹(e₄ⱼ₊₁+e₄ⱼ₊₂+e₄ⱼ₊₃)) ; mise à jour de norme minimale sur (V, Δ)  4n² + 3n

Luminosité (`lum=`) : "couleur" (garde l'ambiguïté), "prior0" (dépliage vers zéro, partie 3.3), "desamb" (oracle de décade,
partie 3.3, TRICHE). Avec un `CapteurDeuxEchelles` la décade est lue directement : plus d'ambiguïté.
Remèdes :
  verif_k, tau, double     relecture à pas réduit de la partie 3.4 (réglage retenu : k = 2, τ = 0,25, double dépliage)
  econome = True           ne relire que si le double dépliage est en désaccord (pas sur le seul écart à la prédiction)
  vote = True              les ambiguïtés sont corrélées : on estime le décalage commun c = médiane des écarts (dépliés,
                           en décades) des composantes sûres, et on redéplie les composantes douteuses autour de pred + c
Test de descente (`test=`) :
  "moy2"  moyenne des ‖ρᵢ‖² < 0,92² (parties 3 et 4)
  "geo"   moyenne des ln‖ρᵢ‖ < ln 0,92 (moyenne géométrique des rapports : une composante qui remonte pèse moins)
  "pond"  Σ mᵢ²‖ρᵢ‖² / Σ mᵢ² < 0,92², mᵢ = ‖Fᵢ‖ estimé en cumulant les variations de luminosité depuis le départ
          (à une décade commune près, inconnue mais sans effet sur le rapport) : on teste la vraie perte ‖F‖², pas la moyenne
          des rapports.

`trace` (liste) reçoit (lectures, q) à chaque pas accepté : c'est le chemin dessiné dans les figures (i, j, k).
`stats` (dict) reçoit des compteurs ; avec `f_vrai`, des diagnostics calculés avec la vraie fonction (sans effet sur
l'algorithme) : décades fausses muettes, et rejets « à tort » (pas refusé alors que la vraie moyenne des ‖ρᵢ‖² passait le test).
"""
import numpy as np
from quat import qmul, qinv, Lmat, Rmat, E4, conf_project, cross3
from sensors import Sensor, Rd, LN10, wrap, _uq
from probes import Hit, OutOfBudget

E0 = E4[0]
CONJ = np.array([1.0, -1.0, -1.0, -1.0])


# ================================================================================================ capteur à deux échelles
class Rd2(Rd):
    __slots__ = ("lg",)


class CapteurDeuxEchelles(Sensor):
    """Couleur-complète-angles avec un canal de luminosité grossier en plus du canal fin (partie 3.5, « deux échelles »).
    À bits égaux : 24 bits = luminosité fine 5 + angle 5 + axe 5 + 5 + luminosité grossière 4 (log10 r sur [lo, hi], 16 niveaux,
    pas de 0,67 décade pour 10 décades). La décade se lit donc directement : k = arrondi(grossier − fin)."""
    def __init__(self, bits=24, sigma=(0.0, 0.0, 0.0), seed=0, coarse_bits=4, lo=-7.0, hi=3.0):
        bf = (bits - coarse_bits) // 4
        al = dict(l=bits - coarse_bits - 3 * bf, t=bf, u=(bf, bf))
        super().__init__("s2", bits, "sph", sigma=sigma, seed=seed, al=al)
        self.cb, self.lo, self.hi = coarse_bits, lo, hi
        self.crng = np.random.default_rng(55_555 + seed)

    def read(self, F):
        r = super().read(F)
        lg = np.log10(max(float(np.linalg.norm(F)), 1e-300))
        if self.sigma[0]:
            lg += self.sigma[0] * self.crng.standard_normal() / LN10
        c = float(_uq(np.clip(lg, self.lo, self.hi), self.lo, self.hi, self.cb))
        frac = r.ell / LN10
        out = Rd2(r.ell, r.th, r.u)
        out.lg = np.round(c - frac) + frac
        return out


# ================================================================================================ compteur avec journal
class MeterLog:
    """Compteur ℍⁿ (comme bench_hn.MeterN) qui journalise chaque lecture : (numéro, ‖F‖, q). mode 'color' | 'exact' | 'degraded'."""
    def __init__(self, f, n, mode, sens, budget, thr=1e-3, log=False):
        self.f, self.n, self.mode, self.sens, self.budget, self.thr, self.cnt = f, n, mode, sens, budget, thr, 0
        self.log = [] if log else None

    def read(self, q):
        self.cnt += 1
        if self.cnt > self.budget:
            raise OutOfBudget
        F = self.f(q)
        if not np.all(np.isfinite(F)):
            raise OutOfBudget
        nF = float(np.linalg.norm(F))
        if self.log is not None:
            self.log.append((self.cnt, nF, np.array(q, float)))
        if nF < self.thr:
            raise Hit
        if self.mode == "exact":
            return F
        Fb = F.reshape(self.n, 4)
        if self.mode == "color":
            return [self.sens[i].read(Fb[i]) for i in range(self.n)]
        if self.mode == "degraded2":      # mesure ABSOLUE à bits égaux : le capteur à deux échelles lu comme une valeur (plage finie)
            out = []
            for i in range(self.n):
                r = self.sens[i].read(Fb[i])
                out.append(10.0 ** r.lg * np.concatenate([[np.cos(r.th)], r.u * np.sin(r.th)]))
            return np.concatenate(out)
        return np.concatenate([self.sens[i].absolute(Fb[i]) for i in range(self.n)])


def capteurs(n, sigma=0.0, i=0, deux_echelles=False, bits=24, fenetre=None):
    """sigma : un réel (même bruit sur les 3 canaux) ou un triplet (σ luminosité en népers, σ angle, σ axe).
    fenetre : (lo, hi) en décades pour la luminosité grossière du capteur à deux échelles (défaut (-7, 3))."""
    sg = tuple(sigma) if isinstance(sigma, (tuple, list)) else (sigma,) * 3
    if deux_echelles:
        kw = {} if fenetre is None else dict(lo=fenetre[0], hi=fenetre[1])
        return [CapteurDeuxEchelles(bits, sigma=sg, seed=1000 * i + k, **kw) for k in range(n)]
    return [Sensor("s2", bits, "sph", sigma=sg, seed=1000 * i + k) for k in range(n)]


# ================================================================================================ quotient
def _dir(r):
    return np.concatenate([[np.cos(r.th)], r.u * np.sin(r.th)])


def quotient(r0, r1, pred, side="L"):
    """ρ (à gauche F0⁻¹F1 ou à droite F1F0⁻¹) depuis deux lectures. Retourne (ρ, ambigu)."""
    if hasattr(r0, "lg"):
        dl, amb = (r1.lg - r0.lg) * LN10, False
    else:
        dl = pred + wrap(r1.ell - r0.ell - pred, LN10)
        amb = abs(dl - pred) > 0.4 * LN10
    P0, P1 = _dir(r0), _dir(r1)
    q = qmul(P0 * CONJ, P1) if side == "L" else qmul(P1, P0 * CONJ)
    return np.exp(dl) * q, amb


# ================================================================================================ modèles de jacobienne
def Lblock(w):
    n = w.shape[0]; W = np.zeros((4 * n, 4 * n))
    for i in range(n):
        for j in range(n):
            W[4 * i:4 * i + 4, 4 * j:4 * j + 4] = Lmat(w[i, j])
    return W


def Rblock(w):
    n = w.shape[0]; W = np.zeros((4 * n, 4 * n))
    for i in range(n):
        for j in range(n):
            W[4 * i:4 * i + 4, 4 * j:4 * j + 4] = Rmat(w[i, j])
    return W


def proj_H(W):
    """Projection orthogonale de chaque bloc sur les multiplications à gauche : coefficients w (n, n, 4)."""
    n = W.shape[0] // 4; w = np.zeros((n, n, 4))
    LB = [Lmat(E4[k]) for k in range(4)]
    for i in range(n):
        for j in range(n):
            B = W[4 * i:4 * i + 4, 4 * j:4 * j + 4]
            w[i, j] = [np.sum(B * LB[k]) / 4 for k in range(4)]
    return w


class Modele:
    """Interface : side ('L'/'R'), nb de sondes, lire(cols) → W, mat(), maj(h, D, rows), rebase(rho)."""
    side = "L"

    def __init__(self, n, **kw):
        self.n = n

    def sondes(self, eta):
        """Directions des sondes (liste de vecteurs de ℝ^{4n})."""
        n = self.n; out = []
        for j in range(n):
            e = np.zeros(4 * n); e[4 * j] = eta; out.append(e)
        return out

    def rebase(self, rho):
        for i in range(self.n):
            Li = Lmat(qinv(rho[i]))
            self.W[4 * i:4 * i + 4, :] = Li @ self.W[4 * i:4 * i + 4, :]

    def mat(self):
        return self.W

    def ok(self):
        return bool(np.all(np.isfinite(self.W)))


class ModeleH(Modele):
    def lire(self, dirs, rhos, eta):
        n = self.n; w = np.zeros((n, n, 4))
        for j in range(n):
            w[:, j, :] = (rhos[j] - E0) / eta
        self.W = Lblock(w)

    def maj(self, h, D, rows):
        n = self.n; hn2 = h @ h
        if hn2 <= 0: return
        A = (D - self.W @ h).reshape(n, 4); hb = h.reshape(n, 4) * CONJ
        Lb = np.hstack([Lmat(hb[j]) for j in range(n)]) / hn2
        for i in np.nonzero(rows)[0]:
            self.W[4 * i:4 * i + 4, :] += Lmat(A[i]) @ Lb


class ModeleD(Modele):
    """Blocs R(wᵢⱼ) : (W h)ᵢ = Σⱼ hⱼ wᵢⱼ. Quotient à droite."""
    side = "R"

    def lire(self, dirs, rhos, eta):
        n = self.n; self.w = np.zeros((n, n, 4))
        for j in range(n):
            self.w[:, j, :] = (rhos[j] - E0) / eta
        self.W = Rblock(self.w)

    def maj(self, h, D, rows):
        n = self.n; hn2 = h @ h
        if hn2 <= 0: return
        A = (D - self.W @ h).reshape(n, 4); hb = h.reshape(n, 4) * CONJ
        for i in np.nonzero(rows)[0]:
            for j in range(n):
                self.w[i, j] += qmul(hb[j], A[i]) / hn2
        self.W = Rblock(self.w)

    def rebase(self, rho):
        for i in range(self.n):
            ri = qinv(rho[i])
            for j in range(self.n):
                self.w[i, j] = qmul(self.w[i, j], ri)
        self.W = Rblock(self.w)


class ModeleDbloc(Modele):
    """Couche entière (partie 6.4) : m sorties, p poids quaternioniques par sortie, n = m·p. Le résidu (i, t) (sortie i, exemple t)
    ne dépend que des poids de la sortie i : la jacobienne est DIAGONALE PAR BLOCS (m blocs p × p), chaque bloc à structure
    « à droite » (blocs R(·), quotient à droite). Conséquence : une sonde sur le poids j de TOUTES les sorties à la fois lit la
    colonne j de tous les blocs ; p sondes au lieu de n = m·p. Mise à jour : Broyden à droite bloc par bloc (chaque bloc ne voit
    que la part h_i du pas qui le concerne, comme le Broyden creux de Schubert). Ordre des indices : résidu i·p + t, poids i·p + j."""
    side = "R"

    def __init__(self, n, m=1, **kw):
        self.n, self.m = n, m; self.p = n // m
        assert self.m * self.p == n

    def sondes(self, eta):
        out = []
        for j in range(self.p):
            e = np.zeros(4 * self.n)
            for i in range(self.m):
                e[4 * (i * self.p + j)] = eta
            out.append(e)
        return out

    def _mk(self):
        m, p = self.m, self.p; W = np.zeros((4 * self.n, 4 * self.n))
        for i in range(m):
            for t in range(p):
                for j in range(p):
                    r, c = i * p + t, i * p + j
                    W[4 * r:4 * r + 4, 4 * c:4 * c + 4] = Rmat(self.w[i, t, j])
        self.W = W

    def lire(self, dirs, rhos, eta):
        m, p = self.m, self.p; self.w = np.zeros((m, p, p, 4))
        for j in range(p):
            rho = rhos[j].reshape(m, p, 4)
            self.w[:, :, j, :] = (rho - E0) / eta
        self._mk()

    def maj(self, h, D, rows):
        m, p = self.m, self.p
        A = (D - self.W @ h).reshape(m, p, 4); hb = h.reshape(m, p, 4) * CONJ
        rws = np.asarray(rows).reshape(m, p)
        for i in range(m):
            hn2 = float(np.sum(hb[i] ** 2))
            if hn2 <= 0: continue
            for t in range(p):
                if not rws[i, t]: continue
                for j in range(p):
                    self.w[i, t, j] += qmul(hb[i, j], A[i, t]) / hn2
        self._mk()

    def rebase(self, rho):
        m, p = self.m, self.p; rr = np.asarray(rho).reshape(m, p, 4)
        for i in range(m):
            for t in range(p):
                ri = qinv(rr[i, t])
                for j in range(p):
                    self.w[i, t, j] = qmul(self.w[i, t, j], ri)
        self._mk()


class ModeleFull(Modele):
    def sondes(self, eta):
        return [eta * np.eye(4 * self.n)[j] for j in range(4 * self.n)]

    def lire(self, dirs, rhos, eta):
        self.W = np.column_stack([(r.reshape(-1) - np.tile(E0, self.n)) / eta for r in rhos])

    def maj(self, h, D, rows):
        if h @ h <= 0: return
        upd = np.outer(D - self.W @ h, h) / (h @ h)
        m = np.repeat(rows, 4)
        self.W[m] += upd[m]


class ModeleBlocConf(Modele):
    """Blocs conformes L(α)R(β) : 3 sondes par colonne de blocs (directions 1, i, j), 4ᵉ colonne par la structure."""
    def sondes(self, eta):
        n = self.n; out = []
        for j in range(n):
            for m in range(3):
                e = np.zeros(4 * n); e[4 * j + m] = eta; out.append(e)
        return out

    def lire(self, dirs, rhos, eta):
        n = self.n; self.W = np.zeros((4 * n, 4 * n))
        for j in range(n):
            cols = [(rhos[3 * j + m] - E0) / eta for m in range(3)]           # (n, 4) chacune
            for i in range(n):
                c = [cols[m][i] for m in range(3)]
                v = cross3(*c)
                lam = np.exp(np.mean([np.log(max(np.linalg.norm(x), 1e-300)) for x in c]))
                v = v / max(np.linalg.norm(v), 1e-300) * lam                      # orientation directe (comme conf3, m = 3)
                B = np.column_stack(c + [v])
                self.W[4 * i:4 * i + 4, 4 * j:4 * j + 4] = conf_project(B)
        self._proj()

    def _proj(self):
        n = self.n
        for i in range(n):
            for j in range(n):
                B = self.W[4 * i:4 * i + 4, 4 * j:4 * j + 4]
                if np.all(np.isfinite(B)):
                    self.W[4 * i:4 * i + 4, 4 * j:4 * j + 4] = conf_project(B)

    def maj(self, h, D, rows):
        if h @ h <= 0: return
        upd = np.outer(D - self.W @ h, h) / (h @ h)
        m = np.repeat(rows, 4)
        self.W[m] += upd[m]
        self._proj()


class ModeleHrang(Modele):
    """W = L(A) + C, C réelle de rang ≤ k, ajustés sur les K dernières sécantes."""
    def __init__(self, n, k=2, K=None, lam=1.0, **kw):
        self.n, self.k, self.K, self.lam = n, k, K or 3 * k, lam
        self.S, self.Y = [], []

    def lire(self, dirs, rhos, eta):
        n = self.n; self.A = np.zeros((n, n, 4))
        for j in range(n):
            self.A[:, j, :] = (rhos[j] - E0) / eta
        self.C = np.zeros((4 * n, 4 * n)); self.S, self.Y = [], []
        self.W = Lblock(self.A)

    def maj(self, h, D, rows):
        if h @ h <= 0: return
        n = self.n
        self.S.append(h.copy()); self.Y.append(D.copy())
        self.S, self.Y = self.S[-self.K:], self.Y[-self.K:]
        S = np.column_stack(self.S); Y = np.column_stack(self.Y)
        for _ in range(2):
            # A : au plus près de A courant, expliquant Y − C S (moindres carrés régularisés, ligne de blocs par ligne)
            Rr = Y - self.C @ S
            for i in range(n):
                # inconnues δa_ij (4n réels) : Σ_j δa_ij h_tj = Σ_j R(h_tj) δa_ij
                M = np.vstack([np.hstack([Rmat(S[4 * j:4 * j + 4, t]) for j in range(n)]) for t in range(S.shape[1])])
                tgt = np.concatenate([Rr[4 * i:4 * i + 4, t] - sum(qmul(self.A[i, j], S[4 * j:4 * j + 4, t]) for j in range(n))
                                      for t in range(S.shape[1])])
                sc = np.mean(np.sum(S * S, 0))
                da = np.linalg.solve(M.T @ M + self.lam * sc * np.eye(4 * n), M.T @ tgt)
                self.A[i] += da.reshape(n, 4)
            LA = Lblock(self.A)
            R = Y - LA @ S
            Cf = R @ np.linalg.pinv(S, rcond=1e-6)
            U, s, Vt = np.linalg.svd(Cf)
            self.C = (U[:, :self.k] * s[:self.k]) @ Vt[:self.k]
        self.W = Lblock(self.A) + self.C

    def rebase(self, rho):
        n = self.n
        for i in range(n):
            ri = qinv(rho[i]); Li = Lmat(ri)
            for j in range(n):
                self.A[i, j] = qmul(ri, self.A[i, j])
            self.C[4 * i:4 * i + 4, :] = Li @ self.C[4 * i:4 * i + 4, :]
            for t in range(len(self.Y)):
                self.Y[t][4 * i:4 * i + 4] = Li @ self.Y[t][4 * i:4 * i + 4]
        self.W = Lblock(self.A) + self.C


class ModeleGrise(Modele):
    """W = L(V)·Δ·L(U), U connue (n, n, 4), Δ diagonale réelle 4n. Normalisation : Δ₄ⱼ = 1 (l'échelle va dans V)."""
    def __init__(self, n, U=None, **kw):
        self.n = n; self.LU = Lblock(U); self.LUi = np.linalg.inv(self.LU)

    def sondes(self, eta):
        n = self.n; out = []
        for j in range(n):
            e = np.zeros(4 * n); e[4 * j] = 1.0
            g = np.zeros(4 * n); g[4 * j + 1:4 * j + 4] = 1.0
            for v in (e, g):
                d = self.LUi @ v
                out.append(eta * d / np.linalg.norm(d))
        return out

    def lire(self, dirs, rhos, eta):
        n = self.n; self.V = np.zeros((n, n, 4)); self.d = np.ones(4 * n)
        for j in range(n):
            for s_, (h, rho) in enumerate(zip(dirs[2 * j:2 * j + 2], rhos[2 * j:2 * j + 2])):
                x = self.LU @ h                                                      # = c·e₄ⱼ ou c·(e₄ⱼ₊₁+e₄ⱼ₊₂+e₄ⱼ₊₃)
                y = (rho - E0).reshape(n, 4)                                         # ≈ W h, blocs i
                if s_ == 0:
                    self.V[:, j, :] = y / x[4 * j]                                   # L(v)·e₀·c = c·v
                else:
                    # Σ_m d_m x_m v·e_m = y_i pour tout i : 3 inconnues d_m, 4n équations
                    Mm = np.column_stack([np.concatenate([qmul(self.V[i, j], E4[m]) * x[4 * j + m] for i in range(n)])
                                          for m in (1, 2, 3)])
                    self.d[4 * j + 1:4 * j + 4] = np.linalg.lstsq(Mm, y.reshape(-1), rcond=None)[0]
        self._mk()

    def _mk(self):
        self.W = Lblock(self.V) @ (self.d[:, None] * self.LU)

    def maj(self, h, D, rows):
        if h @ h <= 0: return
        n = self.n; x = self.LU @ h; r = D - self.W @ h
        # linéarisation : δW h = L(δV)·Δx + L(V)·diag(δd)·x ; inconnues (δV : 4n², δd : 4n), norme minimale
        dx = (self.d * x).reshape(n, 4)
        cols = []
        for i in range(n):
            for j in range(n):
                for m in range(4):
                    c = np.zeros(4 * n); c[4 * i:4 * i + 4] = qmul(E4[m], dx[j]); cols.append(c)
        LV = Lblock(self.V)
        for k in range(4 * n):
            cols.append(LV[:, k] * x[k])
        Amat = np.column_stack(cols)
        for j in range(n):                                                           # Δ₄ⱼ fixé à 1
            Amat[:, 4 * n * n + 4 * j] = 0.0
        m = np.repeat(rows, 4)
        sol = np.linalg.lstsq(Amat[m], r[m], rcond=1e-8)[0]
        self.V += sol[:4 * n * n].reshape(n, n, 4); self.d += sol[4 * n * n:]
        self._mk()

    def rebase(self, rho):
        for i in range(self.n):
            ri = qinv(rho[i])
            for j in range(self.n):
                self.V[i, j] = qmul(ri, self.V[i, j])
        self._mk()


class ModeleHid(ModeleH):
    """A priori « J ≈ identité » (point fixe) : W₀ = diag(Fᵢ⁻¹), direction de Fᵢ lue dans la couleur, décade fixée par UNE sonde
    simultanée η·(1, …, 1) ; puis Broyden quaternionique (sonde-…-départ-identité de la partie 4.2)."""
    def sondes(self, eta):
        e = np.zeros(4 * self.n); e[0::4] = eta
        return [e]

    def lire(self, dirs, rhos, eta):
        n = self.n; w = np.zeros((n, n, 4)); rho = rhos[0]
        for i in range(n):
            r = self.r0[i]
            m = eta / max(np.linalg.norm(rho[i] - E0), 1e-300)
            ell = getattr(r, "lg", None)
            if ell is not None:
                mod = 10 ** ell
            else:
                k = np.round((np.log(m) - r.ell) / LN10); mod = np.exp(r.ell + k * LN10)
            w[i, i] = _dir(r) * CONJ / mod
        self.W = Lblock(w)


MODELES = dict(Hid=ModeleHid, Dbloc=ModeleDbloc, H=ModeleH, D=ModeleD, full=ModeleFull, blocconf=ModeleBlocConf, Hrang=ModeleHrang, grise=ModeleGrise)


# ================================================================================================ sonde
def sonde(M, q0, n, modele="H", c=0.3, step=0.65, cap=1.5, eta0=0.3, acc=0.92, bt_max=5, rcond=1e-3, eta_fixe=False, lum="couleur",
          verif_k=0, tau=0.25, double=False, econome=False, vote=False, test="moy2", f_vrai=None, stats=None, trace=None,
          **kmod):
    sens = M.sens; q = np.array(q0, float); dim = 4 * n
    one = np.tile(E0, n)
    mod = MODELES[modele](n, **kmod); side = mod.side
    st = stats if stats is not None else {}
    for key in ("pas", "acc", "rej_amb", "rej_desc", "rej_tort", "verif", "faux_muet", "faux_muet_acc", "vote"):
        st.setdefault(key, 0)
    oracle = lum == "desamb"
    lacc2, lacc = acc ** 2, np.log(acc)

    def Fex(x):
        return f_vrai(x).reshape(n, 4)

    def quot(r0, r1, pred, F0=None, F1=None):
        rho = np.empty((n, 4)); amb = np.zeros(n, bool); dev = np.zeros(n)
        for i in range(n):
            p = 0.0 if lum == "prior0" else pred[i]
            r, am = quotient(r0[i], r1[i], p, side)
            if oracle:
                dl_vrai = np.log(np.linalg.norm(F1[i]) / np.linalg.norm(F0[i]))
                dl_lu = np.log(max(np.linalg.norm(r), 1e-300))
                kk = np.round((dl_vrai - dl_lu) / LN10)
                r = r * (np.exp(dl_lu + kk * LN10) / max(np.linalg.norm(r), 1e-300)); am = False
            rho[i] = r; amb[i] = am
            dev[i] = (np.log(max(np.linalg.norm(r), 1e-300)) - p) / LN10
        return rho, amb, dev

    def merite(rho, logm):
        if test == "geo":
            return float(np.mean(np.log(np.maximum(np.linalg.norm(rho, axis=1), 1e-300)))) < lacc
        if test == "pond":
            w = np.exp(2 * (logm - logm.max()))
            return float(np.sum(w * np.sum(rho * rho, axis=1)) / np.sum(w)) < lacc2
        return float(np.mean(np.sum(rho * rho, axis=1))) < lacc2

    r0 = M.read(q); eta = eta0; fails = 0; have = False
    F0 = Fex(q) if oracle else None
    logm = np.array([getattr(r, "lg", r.ell / LN10) * LN10 for r in r0])      # ln‖Fᵢ‖ à une décade commune près
    if trace is not None:
        trace.append((M.cnt, q.copy()))

    def lire_W(q, r0, eta, F0):
        for _ in range(40):
            dirs = mod.sondes(eta); rhos = []; bad = False; mod.r0 = r0
            for e in dirs:
                r1 = M.read(q + e)
                F1 = Fex(q + e) if oracle else None
                rho, a, _ = quot(r0, r1, [0.0] * n, F0, F1)
                if a.any():
                    bad = True; break
                rhos.append(rho)
            if bad:
                eta /= 2; continue
            mod.lire(dirs, rhos, eta)
            return eta
        raise OutOfBudget

    while True:
        if not have:
            eta = lire_W(q, r0, eta, F0); fails = 0; have = True
        t, ok = 1.0, False
        for _ in range(bt_max + 1):
            W = mod.mat()
            h = -step * t * np.linalg.lstsq(W, one, rcond=rcond)[0]
            nh = np.linalg.norm(h)
            if nh > cap: h *= cap / nh
            if nh < 1e-13: raise OutOfBudget
            pred = np.log(np.maximum(np.linalg.norm((one + W @ h).reshape(n, 4), axis=1), 1e-12))
            r1 = M.read(q + h)
            F1 = Fex(q + h) if oracle else None
            rho, amb, dev = quot(r0, r1, pred, F0, F1)
            st["pas"] += 1
            if vote and amb.any() and (~amb).sum() >= 2 and lum != "prior0":
                st["vote"] += 1
                cshift = float(np.median(dev[~amb]))
                for i in np.nonzero(amb)[0]:
                    r, am = quotient(r0[i], r1[i], pred[i] + cshift * LN10, side)
                    rho[i], amb[i] = r, am
            if verif_k and not oracle:
                doute = np.zeros(n, bool) if econome else ((np.abs(dev) > tau) | amb)
                if double or econome:
                    for i in range(n):
                        r0i, a0i = quotient(r0[i], r1[i], 0.0, side)
                        if a0i or abs(np.log(max(np.linalg.norm(r0i), 1e-300)) - np.log(max(np.linalg.norm(rho[i]), 1e-300))) > 0.5 * LN10:
                            doute[i] = True
                    if econome:
                        doute |= amb
                if doute.any():
                    st["verif"] += 1
                    hs = h / verif_k
                    rs = M.read(q + hs)
                    for i in np.nonzero(doute)[0]:
                        rsi, asi = quotient(r0[i], rs[i], 0.0, side)
                        pb = float(np.log(max(np.linalg.norm(E0 + verif_k * (rsi - E0)), 1e-12)))
                        r, am = quotient(r0[i], r1[i], pb, side)
                        rho[i] = r; amb[i] = bool(am or asi)
            if f_vrai is not None:                     # DIAGNOSTIC SEULEMENT
                Fq, Fh = Fex(q), Fex(q + h)
                dlv = np.log(np.linalg.norm(Fh, axis=1) / np.linalg.norm(Fq, axis=1))
                dll = np.log(np.maximum(np.linalg.norm(rho, axis=1), 1e-300))
                faux = np.abs(np.round((dlv - dll) / LN10)) >= 1
                st["faux_muet"] += int((faux & ~amb).sum())
                vrai_ok = float(np.mean(np.exp(2 * dlv))) < lacc2
                _fm = bool((faux & ~amb).any())
            any_amb = bool(amb.any())
            D = rho.reshape(-1) - one
            accept = (not any_amb) and merite(rho, logm)
            if accept:
                st["acc"] += 1
                if f_vrai is not None: st["faux_muet_acc"] += int(_fm)
            elif any_amb:
                st["rej_amb"] += 1
            else:
                st["rej_desc"] += 1
                if f_vrai is not None and vrai_ok: st["rej_tort"] += 1
            if not any_amb:
                mod.maj(h, D, np.ones(n, bool))
            if accept:
                logm = logm + np.log(np.maximum(np.linalg.norm(rho, axis=1), 1e-300))
                mod.rebase(rho)
                q, r0, ok = q + h, r1, True
                if oracle: F0 = F1
                if trace is not None: trace.append((M.cnt, q.copy()))
                break
            t *= 0.5
        if not ok:
            fails += 1
            if fails >= 2:
                have = False
                if not eta_fixe:                      # (diagnostic 6.5) ne pas réduire la taille des sondes après un échec de descente
                    eta = max(eta / 2, 1e-9)
                    if eta <= 1e-9: raise OutOfBudget
                else:
                    st["relectures"] = st.get("relectures", 0) + 1
                    if st["relectures"] > 400: raise OutOfBudget
            continue
        fails = 0
        if not mod.ok(): have = False


ERRS = (OutOfBudget, np.linalg.LinAlgError, FloatingPointError, ValueError, ZeroDivisionError, OverflowError)


def lancer(f, q0, n, methode, budget=600, sigma=0.0, i=0, log=False, **kw):
    """Lance une méthode sur (f, q0). methode : 'sonde' (options dans kw), 'broyden' (valeur exacte), 'broyden-deg' (mêmes
    mesures), 'broyden-id' (valeur exacte, départ identité),
    'broyden-q' / 'broyden-q-deg' / 'broyden-q-deg2' (Broyden quaternionique classique : même structure que la sonde ; valeur
    exacte / mêmes mesures à luminosité absolue / mesure absolue à bits égaux, capteur à deux échelles de plage finie). Retourne (succès, lectures, journal, trace, stats)."""
    import classic
    deux = kw.pop("deux_echelles", False); fen = kw.pop("fenetre", None)
    trace, stats = [], {}
    if methode == "sonde":
        M = MeterLog(f, n, "color", capteurs(n, sigma, i, deux, fenetre=fen), budget, log=log)
        call = lambda: sonde(M, q0, n, trace=trace, stats=stats, **kw)
    elif methode == "broyden-deg":
        M = MeterLog(f, n, "degraded", capteurs(n, sigma, i), budget, log=log)
        call = lambda: classic.broyden(M, q0, **kw)
    elif methode == "broyden-id":
        M = MeterLog(f, n, "exact", None, budget, log=log)
        call = lambda: broyden_id(M, q0, **kw)
    elif methode in ("broyden-q", "broyden-q-deg", "broyden-q-deg2"):
        exact = methode == "broyden-q"
        mode = {"broyden-q": "exact", "broyden-q-deg": "degraded", "broyden-q-deg2": "degraded2"}[methode]
        M = MeterLog(f, n, mode, None if exact else capteurs(n, sigma, i, mode == "degraded2", fenetre=fen), budget, log=log)
        kw.setdefault("hfd", 1e-6 if exact else 1e-2)
        call = lambda: broyden_q(M, q0, n, **kw)
    else:
        M = MeterLog(f, n, "exact", None, budget, log=log)
        call = lambda: classic.broyden(M, q0, **kw)
    try:
        call(); ok = False
    except Hit:
        ok = True
    except ERRS:
        ok = False
    return ok, (M.cnt if ok else budget), M.log, trace, stats


def broyden_id(M, q0, bt_max=6, **kw):
    """Broyden « bon », jacobienne initiale = identité (0 lecture) : le bon a priori pour un point fixe h − σ(Wh + b)."""
    q = np.array(q0, float); F = M.read(q); H = np.eye(len(q))
    nrm = np.linalg.norm
    while True:
        d = -H @ F; t = 1.0; Fn = None
        for _ in range(bt_max + 1):
            Fn = M.read(q + t * d)
            if nrm(Fn) < nrm(F): break
            t /= 2; Fn = None
        if Fn is None:
            H = np.eye(len(q)); continue
        s, y = t * d, Fn - F
        Hy = H @ y; den = s @ Hy
        if abs(den) > 1e-14:
            H = H + np.outer(s - Hy, s @ H) / den
        q, F = q + s, Fn


def broyden_q(M, q0, n, side="L", m=1, hfd=1e-6, bt_max=6, suivi=False, **kw):
    """Broyden QUATERNIONIQUE classique : la référence qui connaît la même structure que la sonde, sans la couleur.
    Jacobienne à blocs quaternioniques : L(a) (side='L', h ↦ a·h) ou R(b) (side='R', h ↦ h·b), et, si m > 1, diagonale par
    sortie (m groupes de p = n/m inconnues et de p résidus, comme ModeleDbloc). Départ : une différence finie le long de la
    direction réelle 1 par inconnue quaternionique (L(a)·1 = a, R(b)·1 = b), en perturbant l'inconnue j de tous les groupes à la
    fois : p lectures (n si m = 1), au lieu de 4n pour Broyden-sans-couleur. Ensuite, pas de Newton avec recherche linéaire
    (comme broyden_blocs) et mise à jour de Broyden qui garde la structure, groupe par groupe :
      à gauche  a_rj += (y − Js)_r · s̄_j / ‖s_g‖²      à droite  b_rj += s̄_j · (y − Js)_r / ‖s_g‖²
    (les deux vérifient l'équation sécante J⁺s = y). suivi=True : zéro mobile (le pas est toujours accepté, pas de recherche linéaire). Lit la valeur f (exacte, ou la mesure dégradée : 'broyden-q-deg')."""
    p = n // m; nrm = np.linalg.norm
    q = np.array(q0, float); F = M.read(q)
    blk = np.zeros((m, p, p, 4))                               # blk[g, r, j] : bloc (résidu g·p+r, inconnue g·p+j)
    mk = Lmat if side == "L" else Rmat

    def lire():
        for j in range(p):
            e = np.zeros((n, 4)); e[j::p, 0] = hfd
            dF = ((M.read(q + e.reshape(-1)) - F) / hfd).reshape(m, p, 4)
            blk[:, :, j, :] = dF

    def reel(g):
        Jg = np.zeros((4 * p, 4 * p))
        for r in range(p):
            for j in range(p):
                Jg[4 * r:4 * r + 4, 4 * j:4 * j + 4] = mk(blk[g, r, j])
        return Jg

    lire()
    while True:
        Fg = F.reshape(m, 4 * p)
        d = np.concatenate([-np.linalg.lstsq(reel(g), Fg[g], rcond=1e-10)[0] for g in range(m)])
        t, Fn = 1.0, None
        if suivi:                                       # suivi de zéro mobile : pas toujours accepté, aucune recherche linéaire
            Fn = M.read(q + d)
        else:
            for _ in range(bt_max + 1):
                Fn = M.read(q + t * d)
                if nrm(Fn) < nrm(F): break
                t /= 2; Fn = None
        if Fn is None:
            lire(); continue
        s = (t * d).reshape(m, p, 4); y = (Fn - F).reshape(m, p, 4)
        for g in range(m):
            Js = np.array([sum((qmul(blk[g, r, j], s[g, j]) if side == "L" else qmul(s[g, j], blk[g, r, j])) for j in range(p))
                           for r in range(p)])
            res, ss = y[g] - Js, float(np.sum(s[g] * s[g]))
            if ss <= 0: continue
            for r in range(p):
                for j in range(p):
                    sj = s[g, j] * CONJ
                    blk[g, r, j] += (qmul(res[r], sj) if side == "L" else qmul(sj, res[r])) / ss
        q, F = q + t * d, Fn
