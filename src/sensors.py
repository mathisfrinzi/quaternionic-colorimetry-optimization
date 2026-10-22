"""Capteurs quaternioniques : ce que l'algorithme voit de f(q) (jamais la valeur elle-même).

Un quaternion non nul se factorise F = r e^{uθ} = r (cos θ + u sin θ) : module r, angle θ, axe u ∈ S².

  S0  référence réelle : les 4 composantes de F, chacune en « virgule flottante logarithmique » (1 bit de signe + log|c| quantifié
      sur une plage fixe de décades). Échelle absolue, pas de dent de scie.
  S1  luminosité + teinte, SANS l'axe : luminosité = log10 r modulo 1 décade (dent de scie : tous les bits servent à la précision
      fine, l'ambiguïté de décade se lève par continuité), teinte = θ. Si `axis` est donné, θ est SIGNÉ autour de cet axe
      (θ = atan2(v·axis, F0) ∈ ]-π, π], exactement l'argument complexe quand F reste dans la tranche C_axis) ; sinon θ ∈ [0, π]
      et le signe est perdu.
  S2  luminosité + teinte + axe u. L'axe est encodé soit en RGB (u+1)/2 sur 3 canaux (Sangwine ; redondant puisque |u| = 1),
      soit en 2 angles (colatitude, longitude).

Budget : `bits` est le nombre TOTAL de bits par lecture ; `alloc` le répartit (S2-rgb 24 bits : 6 + 6 + 3×4).

Bruit (avant quantification), un paramètre par canal : σℓ sur ln r (népers), σθ sur θ (radians), σu rotation aléatoire de l'axe
(radians). Même modèle que le bruit multiplicatif exp(σ g1 + iσ g2) du repo complexe.
"""
import numpy as np
from quat import qmul, qconj, qinv, qabs

LN10 = float(np.log(10.0))
TWO_PI = 2 * np.pi


def wrap(x, p):
    """Ramène x modulo p dans [-p/2, p/2[."""
    return (x + p / 2) % p - p / 2


def alloc(kind, B, enc="rgb"):
    """Répartition des B bits entre les canaux."""
    if kind == "s2":
        if enc == "rgb":
            bu = max(1, int(round(B / 6)))
            bl = max(1, (B - 3 * bu) // 2)
            return dict(l=bl, t=B - 3 * bu - bl, u=(bu, bu, bu))
        b = max(1, B // 4)
        return dict(l=b, t=b, u=(b, b))
    if kind == "s1":
        bl = B // 2
        return dict(l=bl, t=B - bl)
    if kind == "s0":
        return dict(m=max(1, B // 4 - 1))
    raise ValueError(kind)


def _uq(x, lo, hi, bits):
    """Quantification uniforme de x ∈ [lo, hi] sur 2^bits niveaux (bornes comprises)."""
    L = (1 << bits) - 1
    i = np.clip(np.rint((x - lo) / (hi - lo) * L), 0, L)
    return lo + i / L * (hi - lo)


def _dq(x, lo, hi, bits, d):
    """Quantification uniforme avec dithering soustractif : on ajoute un décalage connu d ∈ [-1/2, 1/2[ (en pas de
    quantification), on quantifie, puis on retire le décalage. L'erreur devient uniforme, sans biais et indépendante d'une
    lecture à l'autre (au lieu d'être une fonction déterministe de x). La grille est prolongée d'un pas de chaque côté."""
    L = (1 << bits) - 1
    st = (hi - lo) / L
    i = np.clip(np.rint((x + d * st - lo) / st), -1, L + 1)
    return float(lo + i * st - d * st)


def _dq_periodic(x01, bits, d):
    n = 1 << bits
    return ((np.rint(x01 * n + d) - d) / n) % 1.0


def _uq_periodic(x01, bits):
    """Quantification d'une grandeur périodique de période 1 sur 2^bits niveaux."""
    n = 1 << bits
    return (np.rint(x01 * n) % n) / n


class Rd:
    """Une lecture quantifiée. ell : ln r modulo ln10 (S1/S2) ; th : angle ; u : axe (ou None) ; F : composantes (S0)."""
    __slots__ = ("ell", "th", "u", "F")

    def __init__(self, ell=0.0, th=0.0, u=None, F=None):
        self.ell, self.th, self.u, self.F = ell, th, u, F


class Sensor:
    def __init__(self, kind="s2", bits=24, enc="rgb", axis=None, sigma=(0.0, 0.0, 0.0), seed=0,
                 s0_range=(-12.0, 6.0), exact=False, al=None, dither=False):
        self.kind, self.bits, self.enc, self.exact = kind, bits, enc, exact
        self.dither = dither                                   # dithering soustractif connu (partie 7)
        self.drng = np.random.default_rng(10_000_019 + seed)   # suite pseudo-aléatoire connue de l'algorithme
        self.axis = None if axis is None else np.asarray(axis, float) / np.linalg.norm(axis)
        self.sigma = tuple(sigma)
        self.noisy = any(s > 0 for s in self.sigma)
        self.rng = np.random.default_rng(seed)
        self.s0_range = s0_range
        self.al = al if al is not None else alloc(kind, bits, enc)     # al : répartition des bits imposée
        if kind == "s1" and axis is None:
            self.signed = False
        else:
            self.signed = kind == "s1"
        # pas de quantification (pour fixer le seuil « aveugle » des sondes)
        if exact:
            self.step = 1e-12
        elif kind == "s0":
            self.step = (s0_range[1] - s0_range[0]) * LN10 / ((1 << self.al["m"]) - 1)        # pas en log (népers)
        else:
            st_l = LN10 / (1 << self.al["l"])                                                  # népers
            if kind == "s1":
                st_t = (TWO_PI / (1 << self.al["t"])) if self.signed else np.pi / ((1 << self.al["t"]) - 1)
                st_u = 0.0
            else:
                st_t = np.pi / ((1 << self.al["t"]) - 1)
                if enc == "rgb":
                    st_u = 2.0 / ((1 << min(self.al["u"])) - 1)                                # pas par composante ~ angle
                else:
                    st_u = max(np.pi / ((1 << self.al["u"][0]) - 1), TWO_PI / (1 << self.al["u"][1]))
            self.step = max(st_l, st_t, st_u)
        self.thr = 2.0 * self.step                                                              # signal minimal lisible

    # ------------------------------------------------------------------------------------------- lecture
    def read(self, F):
        F = np.asarray(F, float)
        if self.kind == "s0":
            return self._read_s0(F)
        r = float(np.sqrt(F @ F))
        v = F[1:]
        nv = float(np.sqrt(v @ v))
        lg = np.log10(max(r, 1e-300))
        if self.kind == "s1" and self.signed:
            th = float(np.arctan2(v @ self.axis, F[0]))                  # ]-π, π] (axe de référence)
            u = None
        else:
            th = float(np.arctan2(nv, F[0]))                             # [0, π]
            u = v / nv if nv > 1e-300 else np.array([0.0, 0.0, 1.0])
            if self.kind == "s1":
                u = None
        if self.noisy:
            sl, st, su = self.sigma
            if sl: lg += sl * self.rng.standard_normal() / LN10
            if st:
                th += st * self.rng.standard_normal()
                if not (self.kind == "s1" and self.signed):
                    th = abs(th); th = th if th <= np.pi else TWO_PI - th
                else:
                    th = wrap(th, TWO_PI)
            if su and u is not None:
                u = u + su * self.rng.standard_normal(3)
                u = u / np.linalg.norm(u)
        if self.exact:
            return Rd((lg % 1.0) * LN10, th, u)
        # quantification
        if self.dither:
            d = self.drng.uniform(-0.5, 0.5, 4)
            ell = _dq_periodic(lg % 1.0, self.al["l"], d[0]) * LN10
            if self.kind == "s1" and self.signed:
                th = _dq_periodic(((th + np.pi) / TWO_PI) % 1.0, self.al["t"], d[1]) * TWO_PI - np.pi
            else:
                th = _dq(th, 0.0, np.pi, self.al["t"], d[1])
            if self.kind == "s2":
                u = self._quant_axis(u, d[2:])
            return Rd(ell, th, u)
        band = _uq_periodic(lg % 1.0, self.al["l"])
        ell = band * LN10
        if self.kind == "s1" and self.signed:
            th = _uq_periodic(((th + np.pi) / TWO_PI) % 1.0, self.al["t"]) * TWO_PI - np.pi
        else:
            th = float(_uq(th, 0.0, np.pi, self.al["t"]))
        if self.kind == "s2":
            u = self._quant_axis(u)
        return Rd(ell, th, u)

    def _quant_axis(self, u, d=None):
        if self.enc == "rgb":
            c = (u + 1) / 2
            if d is not None:
                dd = [d[0], d[1], self.drng.uniform(-0.5, 0.5)]
                q = np.array([_dq(ci, 0.0, 1.0, b, di) for ci, b, di in zip(c, self.al["u"], dd)])
            else:
                q = np.array([_uq(ci, 0.0, 1.0, b) for ci, b in zip(c, self.al["u"])])
            uq = 2 * q - 1
            n = np.linalg.norm(uq)
            return uq / n if n > 1e-12 else np.array([0.0, 0.0, 1.0])
        phi = np.arccos(np.clip(u[2], -1, 1))
        lam = np.arctan2(u[1], u[0])
        if d is not None:
            phi = float(np.clip(_dq(phi, 0.0, np.pi, self.al["u"][0], d[0]), 0.0, np.pi))
            lam = _dq_periodic(((lam + np.pi) / TWO_PI) % 1.0, self.al["u"][1], d[1]) * TWO_PI - np.pi
        else:
            phi = _uq(phi, 0.0, np.pi, self.al["u"][0])
            lam = _uq_periodic(((lam + np.pi) / TWO_PI) % 1.0, self.al["u"][1]) * TWO_PI - np.pi
        return np.array([np.sin(phi) * np.cos(lam), np.sin(phi) * np.sin(lam), np.cos(phi)])

    def _read_s0(self, F):
        lo, hi = self.s0_range
        out = np.zeros(4)
        for i, c in enumerate(F):
            a = abs(c)
            lg = np.log10(max(a, 1e-300))
            if self.noisy:
                lg += self.sigma[0] * self.rng.standard_normal() / LN10
            if lg < lo:
                continue
            if self.exact:
                out[i] = np.sign(c) * 10 ** lg
                continue
            lgq = _uq(min(lg, hi), lo, hi, self.al["m"])
            out[i] = np.sign(c) * 10 ** lgq
        return Rd(F=out)

    # ------------------------------------------------------------------------------------------- quotient
    def quotient(self, r0, r1, pred=0.0):
        """ρ = F0^{-1} F1 (quaternion) calculé À PARTIR DES LECTURES. Retourne (ρ, ambigu).
        La luminosité est une dent de scie (ln r modulo ln10) : la variation Δℓ = ln|ρ| est dépliée au plus près de la valeur
        PRÉDITE `pred` (0 pour une petite sonde ; ln(1 - pas) pour un pas de Newton). ambigu = l'écart à la prédiction dépasse
        0.4 décade : on ne sait plus de combien |F| a varié. S0 (échelle absolue) n'a pas cette ambiguïté."""
        if self.kind == "s0":
            F0, F1 = r0.F, r1.F
            n0 = F0 @ F0
            if n0 < 1e-280:
                return np.array([1.0, 0, 0, 0]), True
            rho = qmul(qconj(F0) / n0, F1)
            return rho, False
        dl = pred + wrap(r1.ell - r0.ell - pred, LN10)
        amb = abs(dl - pred) > 0.4 * LN10
        P0 = np.concatenate([[np.cos(r0.th)], r0.u * np.sin(r0.th)])
        P1 = np.concatenate([[np.cos(r1.th)], r1.u * np.sin(r1.th)])
        return np.exp(dl) * qmul(qconj(P0), P1), amb

    def partial(self, r0, r1):
        """S1 : (Δℓ, Δθ) entre deux lectures, avec drapeau d'ambiguïté. Pas d'axe, donc pas de quotient quaternionique."""
        dl = wrap(r1.ell - r0.ell, LN10)
        dth = wrap(r1.th - r0.th, TWO_PI) if self.signed else (r1.th - r0.th)
        amb = abs(dl) > 0.4 * LN10 or abs(dth) > 0.8 * np.pi
        return dl, dth, amb

    # ------------------------------------------------------------------------------------------- vue « classique »
    def absolute(self, F):
        """Même résolution et même bruit, mais luminosité NON périodique (échelle absolue conservée) : ce que lit une
        méthode classique à armes égales (favorable à la méthode classique : plage de bits en plus). Retourne F̂ ∈ R^4."""
        F = np.asarray(F, float)
        if self.kind == "s0":
            return self._read_s0(F).F
        if self.kind == "s1":
            raise ValueError("S1 sans axe ne permet pas de reconstruire F")
        r = float(np.sqrt(F @ F))
        v = F[1:]; nv = float(np.sqrt(v @ v))
        lg = np.log10(max(r, 1e-300)); th = float(np.arctan2(nv, F[0]))
        u = v / nv if nv > 1e-300 else np.array([0.0, 0.0, 1.0])
        if self.noisy:
            sl, st, su = self.sigma
            if sl: lg += sl * self.rng.standard_normal() / LN10
            if st:
                th += st * self.rng.standard_normal(); th = abs(th); th = th if th <= np.pi else TWO_PI - th
            if su:
                u = u + su * self.rng.standard_normal(3); u = u / np.linalg.norm(u)
        if not self.exact and self.dither:
            d = self.drng.uniform(-0.5, 0.5, 4)
            n = 1 << self.al["l"]
            lg = (np.rint(lg * n + d[0]) - d[0]) / n
            th = float(np.clip(_dq(th, 0.0, np.pi, self.al["t"], d[1]), 0.0, np.pi))
            u = self._quant_axis(u, d[2:])
        elif not self.exact:
            n = 1 << self.al["l"]
            lg = np.rint(lg * n) / n
            th = float(_uq(th, 0.0, np.pi, self.al["t"]))
            u = self._quant_axis(u)
        return 10 ** lg * np.concatenate([[np.cos(th)], u * np.sin(th)])

    def describe(self):
        return f"{self.kind}{'-' + self.enc if self.kind == 's2' else ''}@{self.bits}b {self.al}"
