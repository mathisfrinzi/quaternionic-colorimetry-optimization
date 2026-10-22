"""Partie 4c : comme 4b, mais |R| n'est plus lu aux points suivants : estimé par les quotients (échelle initiale seule donnée)."""
"""Partie 2 (suite) : ablation des canaux de couleur sur un minimum NON NUL, avec bruit.
Sonde « cible quaternionique » : l'instrument rend F̃ = R (|R| - m0)/|R| (m0 = (1-κ)|R(q)|, niveau cible choisi par l'algorithme comme
dans la sonde à cible du repo complexe : |R| au point COURANT est supposé connu pour fixer m0 et accepter un pas). La couleur de F̃
(luminosité ln|F̃|, teinte θ, axe u) plonge comme pour un vrai zéro (pas d'aplatissement près du minimum). Les quotients ρ_j = F̃0^{-1} F̃_j
donnent R0^{-1}R_j = (|R_j|/|R_0|) ρ_j/|ρ_j| avec |R_j| = m0 + s0|ρ_j| : jacobienne normalisée W en 4 sondes, pas de Levenberg-Marquardt
h = -(WᵀW + μI)^{-1}Wᵀe0 (indépendant de l'échelle de R). 5 lectures par itération (le point d'essai sert de base suivante).
Ablations : (i) bruit sur UN canal à la fois (σℓ, σθ, σu), (ii) bits famés par canal (2 bits au canal testé, le reste réparti).
Succès : |R| < 1.05 en vraie valeur (min |R| = 1). 3 instances × 30 départs par cas, budget 600. Famille nozero (cond. 3) et nozero_ill (cond. 100)."""
import os, sys, json
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "src"))
from funcs import build
from sensors import Sensor
from quat import E4, qmul, qconj

KAPPA = 0.5


def run_suivi(P, q0, sen, eta, mu0=1e-2, budget=600):
    n = [0]
    def Rtrue(q): return P.f(q)
    def read(q, m0):
        n[0] += 1
        R = Rtrue(q); nR = float(np.linalg.norm(R))
        if nR < P.thr:
            raise StopIteration(True)
        if n[0] > budget:
            raise StopIteration(False)
        Ft = R * ((nR - m0) / nR)
        return sen.read(Ft)
    q = np.array(q0, float); mu = mu0
    try:
        nR = float(np.linalg.norm(Rtrue(q)))                    # SEULE information donnée : l'échelle initiale de |R|
        m0 = (1 - KAPPA) * nR; s0 = KAPPA * nR
        r0 = read(q, m0)
        while True:
            W = np.zeros((4, 4)); bad = False
            for j in range(4):
                rj = read(q + eta * E4[j], m0)
                rho, amb = sen.quotient(r0, rj, 0.0)
                nrm = np.linalg.norm(rho)
                if amb or nrm < 1e-9:
                    bad = True; break
                Rj_over_R0 = ((m0 + s0 * nrm) / (m0 + s0)) * rho / nrm
                W[:, j] = (Rj_over_R0 - np.array([1.0, 0, 0, 0])) / eta
            if bad:
                eta = eta / 2 if eta > 1e-4 else eta
                r0 = read(q, m0); continue
            e0 = np.array([1.0, 0, 0, 0])
            while True:
                h = -np.linalg.solve(W.T @ W + mu * np.eye(4), W.T @ e0)
                nh = np.linalg.norm(h)
                if nh > 1.0: h *= 1.0 / nh
                qn = q + h
                rn = read(qn, m0)                               # lecture au point d'essai, même niveau cible m0
                rho_n, amb_n = sen.quotient(r0, rn, 0.0)
                nRn = m0 + s0 * float(np.linalg.norm(rho_n))    # |R| estimé (signe de |R|-m0 supposé positif), sans lire |R|
                if (not amb_n) and nRn < nR:                    # pas accepté (estimation)
                    q, nR, mu = qn, nRn, max(mu / 3, 1e-8)
                    m0 = (1 - KAPPA) * nR; s0 = KAPPA * nR
                    r0 = read(q, m0)
                    break
                mu *= 4
                if mu > 1e8:
                    return False, n[0]
    except StopIteration as e:
        return bool(e.args[0]), n[0]



def case(spec, sen_kw, eta, ntest=30):
    ok, cost = [], []
    for s in (0, 1, 2):
        P = build((spec[0], spec[1], s))
        for i, q0 in enumerate(P.starts(ntest, 1)):
            sen = Sensor(kind="s2", enc="sph", seed=1000 + i, **sen_kw)
            o, c = run_suivi(P, q0, sen, eta); ok.append(o); cost.append(c)
    ok = np.array(ok); cost = np.array(cost)
    return dict(rate=float(ok.mean()), cost=float(np.median(cost[ok])) if ok.any() else None)

if __name__ == "__main__":
    out = []
    for fam in ("nozero", "nozero_ill"):
        for B in (24, 36):
            for sg in (0.0, 0.02):
                r = max([case(("gen", fam), dict(bits=B, sigma=(sg, sg, sg)), e) | dict(eta=e) for e in (0.03, 0.1, 0.3)], key=lambda x: (x["rate"], -(x["cost"] or 1e9)))
                out.append(dict(fam=fam, bits=B, sigma=sg, **r)); print(out[-1], flush=True)
    json.dump(out, open(os.path.join(HERE, "bench_suivi_results.json"), "w"))
