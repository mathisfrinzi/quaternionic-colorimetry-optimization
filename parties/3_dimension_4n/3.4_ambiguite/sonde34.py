"""Partie 3.4 : lever l'ambiguïté de la luminosité sans oracle.

Même algorithme que la sonde-couleur-mémoire-quaternionique (`probe_hn(model='sec-H')`, partie 3.1) et que sa copie instrumentée
de la partie 3.3 (`diag11.probe_diag`), mêmes réglages, sans la mesure coûteuse de l'erreur de W (inutile ici). Les variantes :

  lum = "couleur"  : sonde qui garde l'ambiguïté (référence) ;
  lum = "prior0"   : sonde-dépliage-vers-zéro (partie 3.3) ;
  lum = "desamb"   : sonde-oracle-décade (partie 3.3, triche : la décade vient de la vraie fonction) ;

  verif_k = k > 0  : RELECTURE À PAS RÉDUIT. Après la lecture en q + h, on calcule pour chaque composante l'écart (en décades)
                     entre la variation de luminosité dépliée et la prédiction. Si une composante dépasse le seuil `tau`, on relit
                     UNE fois en q + h/k (une lecture couvre les n composantes). Sur ce petit pas on déplie vers 0, et le quaternion
                     ρˢᵢ = Fᵢ(q)⁻¹ Fᵢ(q + h/k) donne une dérivée directionnelle fraîche : on prédit le grand pas par
                     ln|1 + k(ρˢᵢ − 1)| (au lieu du modèle W, souvent faux) et on redéplie les composantes douteuses autour.
                     Une composante reste ambiguë si l'un des deux dépliages l'est.
  verif_maj        : la paire (h/k, ρˢ − 1) sert aussi de mise à jour de la mémoire (Broyden) avant celle du grand pas.
  tau              : seuil de doute, en décades. Le drapeau d'origine se lève à 0,4 ; en dessous, on vérifie aussi des
                     composantes non signalées (une erreur réelle de plus de 0,5 décade se replie sans drapeau).

  partiel = True   : MISE À JOUR PAR COMPOSANTE (mode `partiel` de diag11, jamais lancé en 3.3). Une composante ambiguë n'invalide
                     plus le pas : le test de descente porte sur les composantes sûres (au moins n/2), seules leurs lignes de W
                     reçoivent la mise à jour sécante ; au rebasage, une composante ambiguë prend la décade la plus proche de la
                     prédiction.
Le code des parties 3.1 et 3.3 n'est pas modifié."""
import os, sys
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "..", "src")); sys.path.insert(0, os.path.join(HERE, "..", "3.1_jacobienne_H_lineaire"))
from quat import qmul, qinv, Lmat, E4
from sensors import LN10
from probes import OutOfBudget
from bench_hn import Lblock, unit_vec

E0 = E4[0]
CONJ = np.array([1, -1, -1, -1])


def _secant(W, h, D, rows, n):
    """Mise à jour ℍ-linéaire w_ij += (D_i − (Wh)_i) conj(h_j) / |h|², sur les lignes `rows` seulement.
    Version matricielle : le bloc (i, j) de W vaut L(w_ij) et L(ab) = L(a) L(b), donc la ligne de blocs i reçoit
    L(D_i − (Wh)_i) · [L(conj h_1) … L(conj h_n)] / |h|² (même calcul que la double boucle de diag11, en O(n·d) produits 4×4)."""
    hn2 = h @ h
    if hn2 <= 0:
        return W
    W = W.copy()
    A = (D - W @ h).reshape(n, 4); hb = h.reshape(n, 4) * CONJ
    Lb = np.hstack([Lmat(hb[j]) for j in range(n)]) / hn2
    for i in np.nonzero(rows)[0]:
        W[4 * i:4 * i + 4, :] += Lmat(A[i]) @ Lb
    return W


def probe34(M, q0, n, c=0.3, step=0.65, cap=1.5, eta0=0.3, acc2=0.92 ** 2, bt_max=5,
            lum="couleur", verif_k=0, tau=0.4, verif_maj=False, partiel=False, double=False, repli=False, f_vrai=None, stats=None):
    """Sonde-couleur-mémoire-quaternionique avec remèdes à l'ambiguïté. `stats` (dict) reçoit des compteurs."""
    sens = M.sens; q = np.array(q0, float); dim = 4 * n
    one = unit_vec(n)
    st = stats if stats is not None else {}
    for key in ("pas", "acc", "rej_amb", "rej_desc", "verif", "verif_acc", "comp_amb", "comp_doute", "sonde_rej",
                "faux_muet", "faux_signale", "fausse_alerte", "faux_muet_verifie", "faux_muet_acc", "repli"):
        st.setdefault(key, 0)
    oracle = lum == "desamb"

    def quot(r0, r1, pred, F0=None, F1=None):
        rho = np.empty((n, 4)); amb = np.zeros(n, bool); dev = np.zeros(n)
        for i in range(n):
            p = 0.0 if lum == "prior0" else pred[i]
            r, am = sens[i].quotient(r0[i], r1[i], p)
            if oracle:
                dl_vrai = np.log(np.linalg.norm(F1[i]) / np.linalg.norm(F0[i]))
                dl_lu = np.log(max(np.linalg.norm(r), 1e-300))
                kk = np.round((dl_vrai - dl_lu) / LN10)
                r = r * (np.exp(dl_lu + kk * LN10) / max(np.linalg.norm(r), 1e-300)); am = False
            rho[i] = r; amb[i] = am
            dev[i] = abs(np.log(max(np.linalg.norm(r), 1e-300)) - p) / LN10
        return rho, amb, dev

    def Fexact(x):
        return f_vrai(x).reshape(n, 4)

    r0 = M.read(q); eta = eta0; W = None; fails = 0
    F0 = Fexact(q) if oracle else None

    def read_W(q, r0, eta, F0):
        for _ in range(40):
            amb = False
            w = np.zeros((n, n, 4))
            for j in range(n):
                e = np.zeros(dim); e[4 * j] = eta
                r1 = M.read(q + e)
                F1 = Fexact(q + e) if oracle else None
                rho, a, _ = quot(r0, r1, [0.0] * n, F0, F1)
                if a.any():
                    amb = True; break
                w[:, j, :] = (rho - E0) / eta
            if amb:
                st["sonde_rej"] += 1
                eta /= 2; continue
            return Lblock(w), eta
        raise OutOfBudget

    while True:
        if W is None:
            W, eta = read_W(q, r0, eta, F0); fails = 0
        t, ok = 1.0, False
        for _ in range(bt_max + 1):
            h = -step * t * np.linalg.lstsq(W, one, rcond=1e-3)[0]
            nh = np.linalg.norm(h)
            if nh > cap: h *= cap / nh; nh = cap
            if nh < 1e-13: raise OutOfBudget
            pred = np.log(np.maximum(np.linalg.norm((one + W @ h).reshape(n, 4), axis=1), 1e-12)).tolist()
            r1 = M.read(q + h)
            F1 = Fexact(q + h) if oracle else None
            rho, amb, dev = quot(r0, r1, pred, F0, F1)
            st["pas"] += 1
            verified = False
            if verif_k and not oracle:
                doute = (dev > tau) | amb
                if double:
                    # second prior : déplier aussi vers 0 ; si les deux dépliages ne tombent pas sur la même décade, on doute
                    for i in range(n):
                        r0i, a0i = sens[i].quotient(r0[i], r1[i], 0.0)
                        if a0i or abs(np.log(max(np.linalg.norm(r0i), 1e-300)) - np.log(max(np.linalg.norm(rho[i]), 1e-300))) > 0.5 * LN10:
                            doute[i] = True
                if doute.any():
                    st["verif"] += 1; st["comp_doute"] += int(doute.sum()); verified = True
                    hs = h / verif_k
                    rs = M.read(q + hs)
                    rho_s = np.empty((n, 4)); amb_s = np.zeros(n, bool)
                    for i in range(n):
                        rho_s[i], amb_s[i] = sens[i].quotient(r0[i], rs[i], 0.0)
                    for i in np.nonzero(doute)[0]:
                        pb = float(np.log(max(np.linalg.norm(E0 + verif_k * (rho_s[i] - E0)), 1e-12)))
                        r, am = sens[i].quotient(r0[i], r1[i], pb)
                        rho[i] = r; amb[i] = bool(am or amb_s[i])
                    if verif_maj and (~amb_s).any():
                        W = _secant(W, hs, rho_s.reshape(-1) - one, ~amb_s, n)
            st["comp_amb"] += int(amb.sum())
            if f_vrai is not None:
                # DIAGNOSTIC SEULEMENT (n'influence aucune décision) : la décade retenue est-elle la bonne ?
                Fq, Fh = Fexact(q), Fexact(q + h)
                dlv = np.log(np.linalg.norm(Fh, axis=1) / np.linalg.norm(Fq, axis=1))
                dll = np.log(np.maximum(np.linalg.norm(rho, axis=1), 1e-300))
                faux = np.abs(np.round((dlv - dll) / LN10)) >= 1
                st["faux_muet"] += int((faux & ~amb).sum()); st["faux_signale"] += int((faux & amb).sum())
                st["fausse_alerte"] += int((~faux & amb).sum())
                if verified: st["faux_muet_verifie"] += int((faux & ~amb & doute).sum())
                _faux_muet = bool((faux & ~amb).any())
            else:
                _faux_muet = False
            D = rho.reshape(-1) - one
            ok_row = ~amb if partiel else np.ones(n, bool)
            any_amb = bool(amb.any())
            phi = float(np.mean(np.sum(rho[ok_row] ** 2, axis=1))) if ok_row.any() else np.inf
            if partiel:
                accept = (ok_row.sum() >= max(1, n // 2)) and phi < acc2
            else:
                accept = (not any_amb) and phi < acc2
            if accept:
                st["acc"] += 1; st["verif_acc"] += int(verified); st["faux_muet_acc"] += int(_faux_muet)
            elif any_amb and not partiel:
                st["rej_amb"] += 1
            else:
                st["rej_desc"] += 1
            if not any_amb or partiel:
                W = _secant(W, h, D, ok_row, n)
            if (not accept) and repli and verified and not amb_s.any():
                # repli : le grand pas échoue, mais la relecture en q + h/k est une vraie lecture, sans ambiguïté ;
                # si elle passe le même test de descente, on l'accepte au lieu de relire plus près (rebroussement gratuit)
                phi_s = float(np.mean(np.sum(rho_s ** 2, axis=1)))
                if phi_s < acc2:
                    if not verif_maj:
                        W = _secant(W, hs, rho_s.reshape(-1) - one, np.ones(n, bool), n)
                    st["repli"] += 1; st["acc"] += 1; st["rej_amb" if (any_amb and not partiel) else "rej_desc"] -= 1
                    rho, h, r1, accept = rho_s, hs, rs, True
                    if oracle: F1 = Fexact(q + h)
            if accept:
                for i in range(n):
                    Li = Lmat(qinv(rho[i]))
                    W[4 * i:4 * i + 4, :] = Li @ W[4 * i:4 * i + 4, :]
                q, r0, ok = q + h, r1, True
                if oracle: F0 = F1
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
