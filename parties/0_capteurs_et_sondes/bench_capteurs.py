"""Partie 0 (M0) : que lit-on avec B bits ? Caractérisation des capteurs S0 / S1 / S2, sans aucun algorithme d'optimisation.

Mesure : une sonde idéale lit F0 puis F1 = F0 ρ avec ρ = 1 + s n (n unitaire, s = 0.3 : amplitude typique d'une sonde, |ρ - 1| = s).
L'erreur est l'écart relatif entre le quotient reconstruit à partir des lectures quantifiées (éventuellement bruitées) et ρ,
‖ρ̂ - ρ‖ / s : c'est l'erreur relative sur une colonne de la jacobienne normalisée W = F^{-1} J.
  1. erreur vs B pour S0 (4 composantes en log-flottant), S2-rgb, S2-sphérique (4 composantes du quotient) et S1 (2 composantes) ;
  2. répartition optimale des bits entre canaux pour S2-rgb à B = 24 et 36 ;
  3. sensibilité au bruit canal par canal (σℓ, σθ, σu seuls) à 48 bits.
Usage : python bench_capteurs.py
"""
import os, sys, json, itertools
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "src"))
from quat import qmul, qabs, rand_unit
from sensors import Sensor, alloc, LN10, wrap

S = 0.3
N = 3000


def sample_pairs(n, seed):
    rng = np.random.default_rng(seed)
    out = []
    for _ in range(n):
        F0 = rand_unit(rng) * 10 ** rng.uniform(-3, 3)
        rho = np.array([1.0, 0, 0, 0]) + S * rand_unit(rng)
        out.append((F0, rho, qmul(F0, rho)))
    return out


PAIRS = sample_pairs(N, 0)


def err_full(sen):
    """Erreur relative médiane sur le quotient quaternionique (S0, S2)."""
    e = []
    for F0, rho, F1 in PAIRS:
        rh, amb = sen.quotient(sen.read(F0), sen.read(F1))
        e.append(min(np.linalg.norm(rh - rho) / S, 10.0) if not amb else 10.0)
    return float(np.median(e)), float(np.mean(np.array(e) >= 10.0))


def err_s1(sen):
    """Erreur relative médiane sur (Δℓ, Δθ), les deux seules composantes que S1 reconstruit."""
    e = []
    for F0, rho, F1 in PAIRS:
        dl, dth, amb = sen.partial(sen.read(F0), sen.read(F1))
        a0, a1 = F0, F1
        th0 = np.arctan2(np.linalg.norm(a0[1:]), a0[0]); th1 = np.arctan2(np.linalg.norm(a1[1:]), a1[0])
        dl_t = np.log(np.linalg.norm(a1) / np.linalg.norm(a0)); dth_t = th1 - th0
        e.append(min(np.hypot(dl - dl_t, dth - dth_t) / S, 10.0) if not amb else 10.0)
    return float(np.median(e)), float(np.mean(np.array(e) >= 10.0))


if __name__ == "__main__":
    out = {}
    Bs = [16, 24, 32, 36, 48, 64, 96]
    curves = {"S0 (4 composantes, log-flottant)": [], "S2 rgb": [], "S2 sphérique": [], "S1 (Δℓ, Δθ) (Δℓ, Δθ)": []}
    for B in Bs:
        curves["S0 (4 composantes, log-flottant)"].append(err_full(Sensor("s0", B)))
        curves["S2 rgb"].append(err_full(Sensor("s2", B, "rgb")))
        curves["S2 sphérique"].append(err_full(Sensor("s2", B, "sph")))
        curves["S1 (Δℓ, Δθ) (Δℓ, Δθ)"].append(err_s1(Sensor("s1", B)))
    out["bits"] = Bs; out["curves"] = curves
    print("erreur relative médiane (fraction d'ambiguïtés) :")
    for k, v in curves.items():
        print(f"  {k:36s}", "  ".join(f"B={B}:{e:.3f}({a:.2f})" for B, (e, a) in zip(Bs, v)))

    # 2. répartition des bits S2-rgb
    alloc_res = {}
    for B in (24, 36):
        rows = []
        for bu in range(2, B // 3):
            for bl in range(2, B):
                bt = B - 3 * bu - bl
                if bt < 2:
                    continue
                e, a = err_full(Sensor("s2", B, "rgb", al=dict(l=bl, t=bt, u=(bu, bu, bu))))
                rows.append(dict(l=bl, t=bt, u=bu, err=e, amb=a))
        rows.sort(key=lambda r: r["err"])
        alloc_res[B] = rows
        d = alloc("s2", B, "rgb")
        dflt = next(r for r in rows if r["l"] == d["l"] and r["t"] == d["t"] and r["u"] == d["u"][0])
        print(f"B={B}: meilleure répartition (ℓ, θ, u×3) = ({rows[0]['l']}, {rows[0]['t']}, {rows[0]['u']}) erreur {rows[0]['err']:.3f} ; "
              f"défaut ({d['l']}, {d['t']}, {d['u'][0]}) erreur {dflt['err']:.3f} ; pire {rows[-1]['err']:.3f}")
    out["alloc"] = alloc_res

    # 3. bruit canal par canal, 48 bits (la quantification devient négligeable)
    sig = [0.002, 0.005, 0.01, 0.02, 0.05]
    noise = {}
    for chan, idx in (("ℓ (luminosité)", 0), ("θ (teinte)", 1), ("u (axe)", 2)):
        v = []
        for s_ in sig:
            sg = [0.0, 0.0, 0.0]; sg[idx] = s_
            v.append(err_full(Sensor("s2", 48, "rgb", sigma=tuple(sg), seed=5))[0])
        noise[chan] = v
        print(f"bruit seul sur {chan:16s}", "  ".join(f"σ={s_}:{e:.3f}" for s_, e in zip(sig, v)))
    out["noise_sigma"] = sig; out["noise"] = noise
    json.dump(out, open(os.path.join(HERE, "bench_capteurs_results.json"), "w"))
