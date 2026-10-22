"""Banc commun aux parties 3.5, 4.3 et 6 : des problèmes ℍⁿ → ℍⁿ nommés, des méthodes nommées, des départs appariés.

Un problème est une fonction `construire(seed) -> (f, q*, infos)` enregistrée dans PROBLEMES par les scripts de banc ; une méthode
est (genre, options) avec genre ∈ {sonde, broyden, broyden-deg, broyden-id} (voir sonde_hn.lancer). Toutes les méthodes voient
les mêmes départs (graine du départ = graine de l'instance + indice), d'où des comparaisons appariées.
"""
import numpy as np
import sonde_hn as S

PROBLEMES = {}
P0 = dict(c=0.3, step=0.65, cap=1.5)


def depart(qs, n, i, seed, rad):
    rng = np.random.default_rng(777 + 1009 * seed + i)
    v = rng.standard_normal(4 * n); v /= np.linalg.norm(v)
    return qs + rad * np.sqrt(n) * v


def run(a):
    """a = (probleme, n, seed, i, methode, (genre, options), sigma, budget, rad)."""
    prob, n, seed, i, meth, (genre, opts), sigma, budget, rad = a
    f, qs, _ = PROBLEMES[prob](n, seed)
    q0 = depart(qs, n, i, seed, rad)
    opts = dict(opts)
    if genre == "sonde":
        opts = {**P0, **opts}
        opts["cap"] = opts["cap"] * max(1.0, np.sqrt(n / 4))
    elif genre in ("broyden", "broyden-deg"):
        opts.setdefault("hfd", 1e-6 if genre == "broyden" else 1e-2)
    ok, c, *_ = S.lancer(f, q0, n, genre, budget=budget, sigma=sigma, i=i, **opts)
    return dict(prob=prob, n=n, seed=seed, i=i, meth=meth, sigma=sigma, ok=bool(ok), cost=int(c))


def wilson(k, m, z=1.96):
    if m == 0: return (np.nan, np.nan)
    p = k / m; d = 1 + z * z / m; c = (p + z * z / (2 * m)) / d; h = z * np.sqrt(p * (1 - p) / m + z * z / (4 * m * m)) / d
    return (c - h, c + h)


def gain_apparie(res, a, b, nboot=2000, seed=0):
    """Gain apparié de a sur b en % de lectures (+ : a moins cher), sur les départs réussis par les deux ; IC 95 % bootstrap."""
    da = {(r["seed"], r["i"]): r for r in res if r["meth"] == a}
    db = {(r["seed"], r["i"]): r for r in res if r["meth"] == b}
    pairs = [(da[k]["cost"], db[k]["cost"]) for k in da if k in db and da[k]["ok"] and db[k]["ok"]]
    if len(pairs) < 3:
        return None
    ca, cb = np.array(pairs, float).T
    g = lambda x, y: 100 * (1 - np.median(x) / np.median(y))
    rng = np.random.default_rng(seed)
    bs = [g(ca[ix], cb[ix]) for ix in (rng.integers(0, len(ca), len(ca)) for _ in range(nboot))]
    return dict(gain=float(g(ca, cb)), lo=float(np.percentile(bs, 2.5)), hi=float(np.percentile(bs, 97.5)), paires=len(pairs))


def resume(res, meth):
    rr = [r for r in res if r["meth"] == meth]
    if not rr: return None
    ok = np.array([r["ok"] for r in rr]); c = np.array([r["cost"] for r in rr], float)
    lo, hi = wilson(int(ok.sum()), len(ok))
    return dict(rate=float(ok.mean()), lo=lo, hi=hi, cost=float(np.median(c[ok])) if ok.any() else None, m=len(ok))


def cellule(s):
    if s is None: return "—"
    return f"{100 * s['rate']:.0f} %" + (f" / {s['cost']:.0f}" if s["cost"] is not None else "")
