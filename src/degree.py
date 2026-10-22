"""Degré topologique lu dans les couleurs (partie 5).

Pour un zéro isolé de f : R^4 -> R^4, l'application g = f/|f| envoie une petite sphère S³ autour du zéro dans S³ ; son degré
(entier, signé) est l'indice du zéro. Un degré non nul PROUVE qu'il y a un zéro dans la boule (théorème de Kronecker-Poincaré),
même si l'algorithme n'a jamais lu |f| en valeur absolue : seule la direction (θ, u) de f, que le capteur S2 donne, est utilisée.
(La luminosité n'est pas utilisée : une lecture de couleur ne dit pas si |f| est petit, elle ne le dit que par continuité.)

Méthode. On triangule S³ par le bord du polytope croisé (16 tétraèdres) raffiné par la subdivision de Bey (chaque tétraèdre en 8),
profondeur d : 16·8^d tétraèdres, et 8, 32, 192, ... sommets (= autant de lectures). On lit la direction de f en chaque sommet,
puis, pour une direction cible w tirée au hasard, on compte les préimages de w de l'application linéaire par morceaux (cône sur
l'image de chaque tétraèdre) : degré = Σ signes. On répète pour plusieurs w et on prend la valeur majoritaire (les lectures sont
déjà faites, c'est gratuit).
"""
import itertools
import numpy as np
from quat import E4


def _bey(tets):
    """Subdivision de Bey : chaque tétraèdre (4 points de R^4) donne 8 tétraèdres de même orientation relative."""
    out = []
    for x0, x1, x2, x3 in tets:
        m = lambda a, b: (a + b) / 2
        x01, x02, x03, x12, x13, x23 = m(x0, x1), m(x0, x2), m(x0, x3), m(x1, x2), m(x1, x3), m(x2, x3)
        out += [(x0, x01, x02, x03), (x01, x1, x12, x13), (x02, x12, x2, x23), (x03, x13, x23, x3),
                (x02, x13, x01, x03), (x02, x13, x03, x23), (x02, x13, x23, x12), (x02, x13, x12, x01)]
    return out


_CACHE = {}


def triangulation(depth):
    """(sommets (N,4) sur S³, tétraèdres (T,4) d'indices, orientation (T,) = signe de det[sommets du tétraèdre])."""
    if depth in _CACHE:
        return _CACHE[depth]
    tets = []
    for s in itertools.product((1.0, -1.0), repeat=4):
        tets.append(tuple(s[i] * E4[i] for i in range(4)))
    for _ in range(depth):
        tets = _bey(tets)
    index, verts, T = {}, [], []
    for tet in tets:
        ids = []
        for p in tet:
            key = tuple(np.round(p, 9))
            if key not in index:
                index[key] = len(verts)
                verts.append(p / np.linalg.norm(p))
            ids.append(index[key])
        T.append(ids)
    V = np.array(verts); T = np.array(T)
    orient = np.sign(np.array([np.linalg.det(V[t].T) for t in T]))
    _CACHE[depth] = (V, T, orient)
    return V, T, orient


def n_reads(depth):
    return len(triangulation(depth)[0])


def degree(G, depth, n_targets=24, seed=0):
    """G : (N, 4) directions (vecteurs unitaires approchés) lues aux sommets de la triangulation de profondeur `depth`.
    Retourne (degré majoritaire, accord = fraction des cibles w qui donnent ce degré)."""
    V, T, orient = triangulation(depth)
    GT = G[T]                                           # (T, 4 sommets, 4)
    A = np.transpose(GT, (0, 2, 1))                     # colonnes = images des sommets
    det = np.linalg.det(A)
    good = np.abs(det) > 1e-10                          # tétraèdres image non dégénérés
    Ag, dg, og = A[good], det[good], orient[good]
    rng = np.random.default_rng(seed)
    vals = []
    for _ in range(n_targets):
        w = rng.standard_normal(4); w /= np.linalg.norm(w)
        lam = np.linalg.solve(Ag, np.broadcast_to(w, (len(Ag), 4))[..., None])[..., 0]  # w = Σ λ_i g_i
        inside = np.all(lam >= 0, axis=1)
        vals.append(int(np.sum(np.sign(dg[inside]) * og[inside])))
    vals = np.array(vals)
    u, c = np.unique(vals, return_counts=True)
    k = int(np.argmax(c))
    return int(u[k]), float(c[k] / len(vals))


def read_direction(sensor, F):
    """Direction (vecteur unitaire de R^4) de F telle que le capteur la restitue : (cos θ, u sin θ), quantifiée, éventuellement bruitée.
    Sensor S2 (θ + axe) ou S0 (composantes) : la direction est reconstructible ; S1 (sans axe) : non."""
    if sensor.kind == "s0":
        c = sensor.read(F).F
        n = np.linalg.norm(c)
        return c / n if n > 0 else np.array([1.0, 0, 0, 0])
    if sensor.kind != "s2":
        raise ValueError("la direction de f exige l'axe : capteur S2 ou S0")
    r = sensor.read(F)
    return np.concatenate([[np.cos(r.th)], r.u * np.sin(r.th)])


def _rotation(seed):
    Q, R = np.linalg.qr(np.random.default_rng(10_000 + seed).standard_normal((4, 4)))
    Q = Q * np.sign(np.diag(R))
    return Q if np.linalg.det(Q) > 0 else Q * np.array([-1.0, 1, 1, 1])


def certificate(f, sensor, center, radius, depth=2, n_targets=24, seed=0):
    """Lit f aux n_reads(depth) sommets d'une sphère (center, radius) ; retourne (degré, accord, nombre de lectures).
    La sphère est lue dans un repère tourné au hasard (une rotation propre garde les orientations) : les sommets d'axe
    sont dégénérés pour des fonctions symétriques comme q² (les trois axes imaginaires ont la même image)."""
    V, T, orient = triangulation(depth)
    Rot = _rotation(seed)
    G = np.array([read_direction(sensor, f(np.asarray(center) + radius * (Rot @ v))) for v in V])
    d, agree = degree(G, depth, n_targets, seed)
    return d, agree, len(V)
