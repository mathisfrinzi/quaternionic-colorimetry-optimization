"""Sondes : algorithmes qui ne voient que des lectures de capteur (jamais la valeur de f).

Idée (héritée du repo complexe). Deux lectures F(q) et F(q + h) donnent le quotient ρ = F(q)^{-1} F(q + h) = 1 + W h + ..., où
W = F(q)^{-1} J est la jacobienne NORMALISÉE : les modules s'éliminent. La condition de zéro F + J h = 0 devient
    W h = -1     (1 = le quaternion réel unité = e0),
et la lecture d'un quotient suffit : ni f ni son échelle ne sont nécessaires.

En complexe, W est un nombre (2 paramètres) : une seule sonde donne tout. En quaternions :
  * générique : W est 4x4 réelle (16 paramètres) -> 4 sondes + 1 lecture de base = 5 lectures ;
  * conforme (Möbius) : W = L(α)R(β) (7 paramètres) -> 3 sondes suffisent (4 lectures), ou 4 sondes avec projection conforme
    (5 lectures, mais le bruit est moyenné sur la contrainte) ;
  * avec réutilisation des lectures (type Broyden) : 1 lecture par pas après l'initialisation ; le modèle conforme est
    ajusté sur une fenêtre de K pas par moindres carrés alternés (α, β).

Les méthodes de ce fichier :
  probe_fd       lecture de W par différences à chaque pas ; model = 'full' | 'conf3' | 'confproj'
  probe_secant   réutilise les lectures ; model = 'rank1' | 'multi' | 'conf'
  probe_s1       capteur S1 (luminosité + teinte, pas d'axe) : 2 équations sur 4 par pas (Newton partiel)
  probe_slice    sonde complexe à 2 lectures par pas, sur la tranche C_i (partie 1)

Toutes lisent via un `Meter` : il compte les lectures, applique le budget, et lève `Hit` dès qu'une lecture tombe sur un point
où |f| < seuil (critère de succès commun à toutes les méthodes, classiques comprises).
"""
import numpy as np
from quat import (qmul, qinv, qabs, qconj, E4, Lmat, Rmat, conf_project, conf_factors, cross3)
from sensors import wrap, LN10


class Hit(Exception):
    pass


class OutOfBudget(Exception):
    pass


class Meter:
    """Compteur de lectures. sensor=None : lecture exacte de F (méthodes classiques, parties 2-5).
    sensor=Sensor, absolute=False : lecture de couleur quantifiée (sondes).
    sensor=Sensor, absolute=True : mêmes mesures dégradées (même résolution, même bruit) mais à échelle absolue, pour les
    méthodes classiques de la partie 6 (« à armes égales »)."""
    def __init__(self, f, sensor=None, budget=600, thr=1e-3, absolute=False):
        self.f, self.sensor, self.budget, self.thr, self.n = f, sensor, budget, thr, 0
        self.absolute = absolute
        self.last = None                       # dernier point lu (pour localiser un succès)

    def read(self, q):
        self.n += 1
        if self.n > self.budget:
            raise OutOfBudget
        self.last = np.array(q, float)
        F = self.f(q)
        if not np.all(np.isfinite(F)):
            raise OutOfBudget
        nF = float(np.sqrt(F @ F))
        if nF < self.thr:
            raise Hit
        if self.sensor is None:
            return F
        return self.sensor.absolute(F) if self.absolute else self.sensor.read(F)

    def charge(self, k=1):
        """Facture k évaluations sans lire (ex. jacobienne exacte facturée 1 évaluation, convention du repo complexe)."""
        self.n += k
        if self.n > self.budget:
            raise OutOfBudget


E0 = E4[0]


def _solve_step(W, scale):
    """h = -scale * W^+ e0 (moindres carrés si W est (presque) singulière)."""
    sol = np.linalg.lstsq(W, E0, rcond=1e-3)[0]
    return -scale * sol


# ------------------------------------------------------------------------------------------------ modèle local par sondes
def _assemble(model, cols, it):
    """cols : dict {j: colonne de W}. Retourne W (4x4)."""
    if model == "full":
        return np.column_stack([cols[j] for j in range(4)])
    if model == "confproj":
        return conf_project(np.column_stack([cols[j] for j in range(4)]))
    if model == "conf3":
        m = [j for j in range(4) if j not in cols][0]
        others = [cols[j] for j in range(4) if j != m]
        n = cross3(*others)
        nn = np.linalg.norm(n)
        lam = np.exp(np.mean([np.log(max(np.linalg.norm(c), 1e-300)) for c in others]))
        w = (-1) ** (3 - m) * n / max(nn, 1e-300) * lam
        full = {**cols, m: w}
        return conf_project(np.column_stack([full[j] for j in range(4)]))
    raise ValueError(model)


def read_model(M, q, r0, eta, model, it, c, eta_max, thr_scale=1.0):
    """Lit W en q (r0 = lecture en q). Ajuste eta si une sonde est ambiguë (trop grand) ou aveugle (trop petit).
    Retourne (W, eta)."""
    sen = M.sensor
    idx = list(range(4)) if model in ("full", "confproj") else [j for j in range(4) if j != it % 4]
    for _ in range(60):
        cols, amb = {}, False
        for j in idx:
            rj = M.read(q + eta * E4[j])
            rho, a = sen.quotient(r0, rj)
            if a:
                amb = True
                break
            cols[j] = (rho - E0) / eta
        if amb:
            eta /= 2
            if eta < 1e-10:
                raise OutOfBudget
            continue
        sig = max(np.linalg.norm(v) for v in cols.values()) * eta
        if sig < sen.thr * thr_scale:                        # aveugle : le signal est sous le pas de quantification
            if eta >= eta_max:
                raise OutOfBudget
            eta = min(2 * eta, eta_max)
            continue
        return _assemble(model, cols, it), eta
    raise OutOfBudget


def probe_fd(M, q0, model="full", c=0.3, step=0.65, cap=1.0, eta0=0.3, eta_max=4.0, acc=0.92, bt_max=5):
    """Lecture de W à chaque pas (5 lectures par pas en 'full'/'confproj', 4 en 'conf3'), pas de Newton amorti, test de
    descente lu dans la couleur (|ρ| = |F_nouveau|/|F_ancien| doit baisser)."""
    q = np.array(q0, float)
    sen = M.sensor
    r0 = M.read(q)
    eta, it = eta0, 0
    while True:
        it += 1
        W, eta = read_model(M, q, r0, eta, model, it, c, eta_max)
        t, ok = 1.0, False
        for _ in range(bt_max + 1):
            h = _solve_step(W, step * t)
            nh = np.linalg.norm(h)
            if nh > cap:
                h *= cap / nh
            if nh < 1e-14:
                raise OutOfBudget
            pred = float(np.log(max(np.linalg.norm(E0 + W @ h), 1e-12)))
            r1 = M.read(q + h)
            rho, amb = sen.quotient(r0, r1, pred)
            if not amb and qabs(rho) < acc:
                q, r0, ok = q + h, r1, True
                break
            t *= 0.5
        if not ok:
            eta = max(eta / 2, 1e-10)
            if eta <= 1e-10:
                raise OutOfBudget
            continue
        eta = float(np.clip(c / max(np.linalg.norm(W, 2), 1e-12), 1e-10, eta_max))


# ------------------------------------------------------------------------------------------------ ajustement conforme
def conf_als(H, D, a, b, ridge=1e-3, iters=4):
    """Moindres carrés alternés pour D_k ≈ a h_k b (produits de Hamilton), partant de (a, b) avec un terme de rappel
    vers l'estimation précédente (la fenêtre peut avoir moins de données que de paramètres)."""
    y = np.concatenate(D)
    for _ in range(iters):
        A = np.vstack([Rmat(qmul(h, b)) for h in H])                 # D_k = R(h_k b) a
        lam = ridge * np.trace(A.T @ A) / 4
        a = np.linalg.solve(A.T @ A + lam * np.eye(4), A.T @ y + lam * a)
        A = np.vstack([Lmat(qmul(a, h)) for h in H])                 # D_k = L(a h_k) b
        lam = ridge * np.trace(A.T @ A) / 4
        b = np.linalg.solve(A.T @ A + lam * np.eye(4), A.T @ y + lam * b)
        nb = np.linalg.norm(b)
        b, a = b / nb, a * nb                                         # jauge |b| = 1
    return a, b


def probe_secant(M, q0, model="rank1", K=3, c=0.3, step=0.65, cap=1.0, eta0=0.3, eta_max=4.0, acc=0.92, bt_max=5,
                 init="full", ridge=1e-3, tr=False):
    """Réutilise les lectures : après l'initialisation (lecture de W), 1 lecture par pas.
    model : 'rank1' (Broyden, 1 sécante) | 'multi' (K sécantes, variation minimale de W) | 'conf' (K sécantes, W = L(α)R(β)).
    Après un pas accepté de quotient ρ, le changement de base multiplie W (ou α) à gauche par ρ^{-1} ; les sécantes
    mémorisées sont ramenées dans la nouvelle base de la même façon.
    init : 'full' (4 sondes + base) | 'conf' (démarrage conforme : 3 sondes, W projeté sur L(α)R(β), puis mises à jour
           générales ; partie 7).
    tr   : région de confiance pilotée par la luminosité (partie 7). Le pas reste le pas de Newton amorti (step), mais sa
           longueur maximale n'est plus le plafond fixe `cap` : c'est un rayon Δ adaptatif : on compare la
           baisse de log|F| lue (1 - |ρ|, |ρ| se lit dans la luminosité) à la baisse prédite par le modèle (1 - |1 + W h|) ; rapport > 0,75
           et pas au bord → Δ × 2 ; rapport < 0,25 ou pas refusé → Δ réduit. Le rayon est GARDÉ d'un pas à l'autre."""
    q = np.array(q0, float)
    sen = M.sensor
    r0 = M.read(q)
    eta = eta0
    W, pairs, a, b = None, [], None, None
    fails, it = 0, 0
    init_model = "conf3" if init == "conf" else "full"
    rad = cap
    while True:
        if W is None:
            it += 1
            W, eta = read_model(M, q, r0, eta, init_model, it, c, eta_max)
            pairs = []
            if model == "conf":
                a, b = conf_factors(W)
                W = Lmat(a) @ Rmat(b)
            fails = 0
        t = 1.0
        ok = False
        for _ in range(bt_max + 1):
            h = _solve_step(W, step * t)
            nh = np.linalg.norm(h)
            lim = rad if tr else cap
            if nh > lim:
                h *= lim / nh
            if nh < 1e-14:
                raise OutOfBudget
            pred = float(np.log(max(np.linalg.norm(E0 + W @ h), 1e-12)))
            r1 = M.read(q + h)
            rho, amb = sen.quotient(r0, r1, pred)
            if amb:
                t *= 0.5
                if tr:
                    rad = max(0.5 * min(rad, np.linalg.norm(h)), cap / 64)
                continue
            if tr:
                # baisses relatives de |F| (lue : 1 - |ρ| ; prédite : 1 - |1 + W h|), bornées : la prédiction d'un pas de Newton
                # complet vaut 1, le logarithme serait infini
                actual = 1.0 - qabs(rho)
                predicted = 1.0 - float(np.linalg.norm(E0 + W @ h))
                ratio = actual / predicted if predicted > 1e-9 else 0.0
                nh2 = np.linalg.norm(h)
                if ratio > 0.75 and nh2 > 0.9 * rad:
                    rad = min(2 * rad, 4 * cap)
                elif ratio < 0.25:
                    rad = max(0.35 * nh2, cap / 64)
            D = rho - E0
            # mise à jour du modèle avec la sécante (h, D), valable que le pas soit accepté ou non
            if model == "rank1":
                W = W + np.outer(D - W @ h, h) / (h @ h)
            else:
                pairs.append((h.copy(), D.copy()))
                pairs = pairs[-K:]
                Hs = [p[0] for p in pairs]; Ds = [p[1] for p in pairs]
                if model == "multi":
                    Hm = np.column_stack(Hs); Dm = np.column_stack(Ds)
                    G = Hm.T @ Hm
                    W = W + (Dm - W @ Hm) @ np.linalg.solve(G + 1e-9 * np.trace(G) / len(Hs) * np.eye(len(Hs)), Hm.T)
                else:
                    a, b = conf_als(Hs, Ds, a, b, ridge=ridge)
                    W = Lmat(a) @ Rmat(b)
            if qabs(rho) < acc:
                # rebase : W <- ρ^{-1} W ; sécantes et α idem
                ri = qinv(rho)
                W = Lmat(ri) @ W
                if model == "conf":
                    a = qmul(ri, a)
                if pairs:
                    pairs = [(hh, qmul(ri, dd)) for hh, dd in pairs]
                q, r0, ok = q + h, r1, True
                break
            t *= 0.5
            if tr:
                rad = max(0.5 * min(rad, np.linalg.norm(h)), cap / 64)
        if not ok:
            fails += 1
            if fails >= 2:
                W = None
                rad = cap                                             # nouvelle jacobienne : rayon remis à zéro
                eta = max(eta / 2, 1e-10)
                if eta <= 1e-10:
                    raise OutOfBudget
            continue
        fails = 0
        if not np.all(np.isfinite(W)):
            W = None


# ------------------------------------------------------------------------------------------------ capteur S1 (sans axe)
def probe_s1(M, q0, c=0.3, step=0.65, cap=1.0, eta0=0.3, eta_max=4.0, acc=0.92, bt_max=5, reg=1e-8):
    """S1 ne donne ni l'axe ni donc le quotient quaternionique, seulement (Δℓ, Δθ). Au premier ordre :
    Δℓ ≈ Re(W h) et Δθ ≈ (composante de W h sur l'axe courant). L'équation W h = -1 impose Re = -1 et composante radiale = 0 :
    deux équations sur quatre. h = solution de norme minimale de [∇ℓ ; ∇θ] h = (-step, 0). L'axe tourne d'un pas à
    l'autre, donc les équations utilisées changent (type Kaczmarz)."""
    q = np.array(q0, float)
    sen = M.sensor
    r0 = M.read(q)
    eta = eta0
    while True:
        for _ in range(60):
            rows_l, rows_t, amb = np.zeros(4), np.zeros(4), False
            for j in range(4):
                rj = M.read(q + eta * E4[j])
                dl, dth, a = sen.partial(r0, rj)
                if a:
                    amb = True
                    break
                rows_l[j], rows_t[j] = dl / eta, dth / eta
            if amb:
                eta /= 2
                if eta < 1e-10:
                    raise OutOfBudget
                continue
            if max(np.abs(rows_l).max(), np.abs(rows_t).max()) * eta < sen.thr:
                if eta >= eta_max:
                    raise OutOfBudget
                eta = min(2 * eta, eta_max)
                continue
            break
        else:
            raise OutOfBudget
        G = np.vstack([rows_l, rows_t])
        t, ok = 1.0, False
        for _ in range(bt_max + 1):
            rhs = np.array([-step * t, 0.0])
            h = G.T @ np.linalg.solve(G @ G.T + reg * np.trace(G @ G.T) * np.eye(2), rhs)
            nh = np.linalg.norm(h)
            if nh > cap:
                h *= cap / nh
            if nh < 1e-14:
                raise OutOfBudget
            pred = float(np.log(max(1 - step * t, 1e-6)))
            r1 = M.read(q + h)
            dl, dth, a = sen.partial(r0, r1)
            dl_unf = pred + wrap(dl - pred, LN10)
            if abs(dl_unf - pred) <= 0.4 * LN10 and np.exp(dl_unf) < acc:
                q, r0, ok = q + h, r1, True
                break
            t *= 0.5
        if not ok:
            eta = max(eta / 2, 1e-10)
            if eta <= 1e-10:
                raise OutOfBudget
            continue
        eta = float(np.clip(c / max(np.linalg.norm(G, 2), 1e-12), 1e-10, eta_max))


# ------------------------------------------------------------------------------------------------ partie 1 : tranche
def probe_slice(M, z0, c=0.3, step_frac=0.7, cap=1.0, eta0=0.3, eta_max=8.0, blind_rule=True, maxit=10 ** 6):
    """Sonde complexe du repo (descend_probe), portée sur le capteur S1 SIGNÉ et la tranche C_i : q = (x, y, 0, 0).
    Deux lectures par pas ; L = Δℓ + iΔθ ≈ η f'/f (dérivée logarithmique complexe) ; descente -conj(f'/f).
    Assertion de tranche : q ne quitte jamais C_i (fausse piste 2)."""
    sen = M.sensor
    z = complex(z0[0], z0[1]) if np.ndim(z0) else complex(z0)      # accepte un quaternion de départ (sur C_i) ou un complexe
    eta, u = eta0, 1.0 + 0j
    thr = sen.thr
    def rd(zz):
        return M.read(np.array([zz.real, zz.imag, 0.0, 0.0]))
    while True:
        r0 = rd(z)
        r1 = rd(z + eta * u)
        dl, dth, amb = sen.partial(r0, r1)
        if abs(dl) > 0.4 * LN10 or abs(dth) > 0.8 * np.pi:
            eta /= 2
            if eta < 1e-12:
                raise OutOfBudget
            continue
        L = dl + 1j * dth
        if abs(L) < thr and blind_rule:
            for k in (1, 2, 3):
                u2 = u * 1j ** k
                r2 = rd(z + eta * u2)
                d2l, d2t, _ = sen.partial(r0, r2)
                if abs(d2l) > 0.4 * LN10 or abs(d2t) > 0.8 * np.pi:
                    continue
                if abs(d2l + 1j * d2t) >= thr:
                    L, u = d2l + 1j * d2t, u2
                    break
        if abs(L) < thr:
            if eta >= eta_max:
                raise OutOfBudget
            eta = min(2 * eta, eta_max)
            continue
        w = L / (eta * u)
        z = z - min(step_frac / abs(w), cap) * np.conj(w) / abs(w)
        eta = float(np.clip(c / abs(w), 1e-9, eta_max))
