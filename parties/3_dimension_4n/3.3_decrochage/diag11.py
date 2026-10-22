"""Partie 3.3 : pourquoi et quand `sec-H` décroche quand n grandit.

Copie instrumentée de `probe_hn(model='sec-H')` (partie 8) : même algorithme, mêmes réglages, mais
  - chaque lecture de pas est journalisée (cause de rejet, nombre de composantes ambiguës, taille du pas,
    erreur de prédiction de la luminosité, erreur du modèle W par rapport à la vraie jacobienne normalisée) ;
  - des interrupteurs « oracle », qui lèvent UNE cause à la fois :
      lum="exact"   : la variation de luminosité ln|ρ_i| est lue exacte (plus de dent de scie, donc plus d'ambiguïté) ;
      W="vraie"     : W est remplacée à chaque pas par diag(F)^-1 J exacte (plus d'erreur de modèle sécant) ;
      partiel=True  : une composante ambiguë n'invalide plus le pas ; elle est exclue du test d'acceptation et de la
                      mise à jour sécante (la mise à jour de la ligne i de W ne dépend que de la composante i).
  amorti=a      : la mise à jour sécante divise par |h|² + a·thr² (thr = pas de quantification du capteur) au lieu de |h|² :
                  un pas plus petit que le bruit de lecture ne doit pas écraser W.
Le code de la partie 8 n'est pas modifié."""
import os, sys
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "..", "src")); sys.path.insert(0, os.path.join(HERE, "..", "3.1_jacobienne_H_lineaire"))
from quat import qmul, qinv, Lmat, E4
from sensors import Sensor, LN10
from probes import Hit, OutOfBudget
from bench_hn import MeterN, Lblock, unit_vec

E0 = E4[0]
CONJ = np.array([1, -1, -1, -1])


def jac_vraie_normalisee(f, q, n, h=1e-6):
    """W = diag(F)^-1 J par différences finies centrées sur f (lectures NON comptées : c'est un oracle)."""
    d = 4 * n
    J = np.empty((d, d))
    for k in range(d):
        e = np.zeros(d); e[k] = h
        J[:, k] = (f(q + e) - f(q - e)) / (2 * h)
    F = f(q).reshape(n, 4)
    W = np.empty((d, d))
    for i in range(n):
        Fi = F[i]
        W[4 * i:4 * i + 4] = Lmat(qinv(Fi)) @ J[4 * i:4 * i + 4]
    return W


def probe_diag(M, q0, n, c=0.3, step=0.65, cap=1.5, eta0=0.3, acc2=0.92 ** 2, bt_max=5,
               lum="couleur", W_oracle=False, partiel=False, log=None, f_vrai=None, amorti=0.0):
    """sec-H instrumenté. `log` : liste où l'on ajoute un dict par lecture de pas (et par lecture de sonde)."""
    sens = M.sens; q = np.array(q0, float); dim = 4 * n
    one = unit_vec(n)
    ev = log if log is not None else []

    def quot(r0, r1, pred, F0=None, F1=None):
        rho = np.empty((n, 4)); amb = np.zeros(n, bool); perr = np.zeros(n)
        for i in range(n):
            r, am = sens[i].quotient(r0[i], r1[i], 0.0 if lum == "prior0" else pred[i])
            if lum == "noflag":
                am = False        # on garde la décade « la plus proche de la prédiction » même si elle est douteuse
            if lum == "desamb":
                # on garde la lecture bruitée/quantifiée de la luminosité, mais la décade est connue : seule la dent de scie
                # est levée (le bruit et la quantification de ℓ restent)
                dl_vrai = np.log(np.linalg.norm(F1[i]) / np.linalg.norm(F0[i]))
                dl_lu = np.log(max(np.linalg.norm(r), 1e-300))
                k = np.round((dl_vrai - dl_lu) / LN10)
                r = r * (np.exp(dl_lu + k * LN10) / max(np.linalg.norm(r), 1e-300)); am = False
            if lum == "exact":
                # luminosité exacte : on garde la direction lue dans la couleur, mais |ρ_i| vient de |F_i| exacts
                dl = np.log(np.linalg.norm(F1[i]) / np.linalg.norm(F0[i]))
                nr = np.linalg.norm(r)
                r = r * (np.exp(dl) / max(nr, 1e-300)); am = False
            rho[i] = r; amb[i] = am
            perr[i] = abs(np.log(max(np.linalg.norm(r), 1e-300)) - pred[i])
        return rho, amb, perr

    def Fexact(x):
        return f_vrai(x).reshape(n, 4)

    r0 = M.read(q); eta = eta0; W = None; fails = 0
    F0 = Fexact(q) if lum in ("exact", "desamb") else None

    def read_W(q, r0, eta, F0):
        for _ in range(40):
            amb = False
            w = np.zeros((n, n, 4))
            for j in range(n):
                e = np.zeros(dim); e[4 * j] = eta
                r1 = M.read(q + e)
                F1 = Fexact(q + e) if lum in ("exact", "desamb") else None
                rho, a, _ = quot(r0, r1, [0.0] * n, F0, F1)
                if a.any():
                    amb = True; break
                w[:, j, :] = (rho - E0) / eta
            ev.append(dict(kind="sonde", n_amb=int(amb), ok=not amb))
            if amb:
                eta /= 2; continue
            return Lblock(w), eta
        raise OutOfBudget

    while True:
        if W is None:
            W, eta = read_W(q, r0, eta, F0); fails = 0
        t, ok = 1.0, False
        for _ in range(bt_max + 1):
            if W_oracle:
                W = jac_vraie_normalisee(f_vrai, q, n)
            h = -step * t * np.linalg.lstsq(W, one, rcond=1e-3)[0]
            nh = np.linalg.norm(h)
            if nh > cap: h *= cap / nh; nh = cap
            if nh < 1e-13: raise OutOfBudget
            pred = [float(np.log(max(np.linalg.norm((one + W @ h).reshape(n, 4)[i]), 1e-12))) for i in range(n)]
            r1 = M.read(q + h)
            F1 = Fexact(q + h) if lum in ("exact", "desamb") else None
            rho, amb, perr = quot(r0, r1, pred, F0, F1)
            D = rho.reshape(-1) - one
            ok_row = ~amb if partiel else np.ones(n, bool)
            any_amb = bool(amb.any())
            phi = float(np.mean(np.sum(rho[ok_row] ** 2, axis=1))) if ok_row.any() else np.inf
            if partiel:
                accept = (ok_row.sum() >= max(1, n // 2)) and phi < acc2
            else:
                accept = (not any_amb) and phi < acc2
            # erreur de modèle (si l'oracle est disponible) : ‖W − W_vrai‖ / ‖W_vrai‖, sur la base de q
            werr = None
            if f_vrai is not None and not W_oracle and len(ev) % 1 == 0:
                Wv = jac_vraie_normalisee(f_vrai, q, n)
                werr = float(np.linalg.norm(W - Wv) / np.linalg.norm(Wv))
            cause = "accepte" if accept else ("ambigu" if (any_amb and not partiel) else "descente")
            ev.append(dict(kind="pas", n_amb=int(amb.sum()), nh=float(nh), phi=phi, perr=float(perr.mean()),
                           werr=werr, cause=cause, t=t))
            if (not W_oracle) and (not any_amb or partiel):
                if True:                                    # mise à jour sécante ℍ-linéaire, ligne par ligne
                    wmat = np.zeros((n, n, 4))
                    for i in range(n):
                        for j in range(n):
                            wmat[i, j] = W[4 * i:4 * i + 4, 4 * j:4 * j + 4][:, 0]
                    Wh = (W @ h).reshape(n, 4); hb = h.reshape(n, 4); Db = D.reshape(n, 4)
                    # amortissement : |h|² + λ au dénominateur, λ = amorti × (pas de quantification du capteur)²
                    hn2 = h @ h + amorti * sens[0].thr ** 2
                    for i in range(n):
                        if not ok_row[i]:
                            continue
                        for j in range(n):
                            wmat[i, j] += qmul(Db[i] - Wh[i], hb[j] * CONJ) / hn2
                    W = Lblock(wmat)
            if accept:
                for i in range(n):
                    # rebase de TOUTES les composantes ; en mode partiel, une composante ambiguë est rebasée avec la valeur
                    # dépliée au plus près de la prédiction (juste tant que l'écart réel est < 0,5 décade)
                    Li = Lmat(qinv(rho[i]))
                    W[4 * i:4 * i + 4, :] = Li @ W[4 * i:4 * i + 4, :]
                q, r0, ok = q + h, r1, True
                if lum in ("exact", "desamb"): F0 = F1
                break
            t *= 0.5
        if not ok:
            fails += 1
            if fails >= 2:
                W = None; eta = max(eta / 2, 1e-9)
                if eta <= 1e-9: raise OutOfBudget
            continue
        fails = 0
        if not np.all(np.isfinite(W)): W = None
