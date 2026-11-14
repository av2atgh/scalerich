"""Reference lattices in the unit square (same frame model): square, triangular, honeycomb.
Returns nodes, seg; boundary nodes = nodes within tol of the square's sides
(lattices are built so the outermost nodes lie on the sides)."""
import numpy as np


def _bonds(P, a, tol=1e-6):
    from scipy.spatial import cKDTree
    pairs = np.array(sorted(cKDTree(P).query_pairs(a * (1 + 1e-6))))
    return pairs


def _finish(P, seg):
    used = np.unique(seg)
    idx = -np.ones(len(P), int); idx[used] = np.arange(len(used))
    P = P[used]; seg = idx[seg]
    tol = 1e-9
    bnd = np.where((P[:, 0] < tol) | (P[:, 0] > 1 - tol) | (P[:, 1] < tol) | (P[:, 1] > 1 - tol))[0]
    return P, seg, bnd


def square(n=12, angle=0.0):
    x = np.linspace(0, 1, n + 1)
    P = np.array([(u, v) for u in x for v in x])
    return _finish(P, _bonds(P, 1.0 / n))


def triangular(n=12):
    a = 1.0 / n; h = a * np.sqrt(3) / 2
    rows = int(np.floor(1 / h))
    P = []
    for j in range(rows + 1):
        off = 0.5 * a * (j % 2)
        for i in range(-1, n + 2):
            P.append((off + i * a, j * h))
    P = np.array(P)
    seg = _bonds(P, a)
    keep = np.all((P[seg][:, :, 0] > -1e-9) & (P[seg][:, :, 0] < 1 + 1e-9), axis=1)
    seg = seg[keep]
    P[:, 1] /= rows * h          # stretch to unit height (tiny distortion)
    return _finish(P, seg)


def honeycomb(n=8):
    """armchair-free: zigzag rows; bonds of length a; clipped to the square."""
    a = 1.0 / (n * np.sqrt(3))
    a1 = np.array([np.sqrt(3) * a, 0]); a2 = np.array([np.sqrt(3) * a / 2, 1.5 * a])
    P = []
    for i in range(-3 * n, 3 * n):
        for j in range(-3 * n, 3 * n):
            o = i * a1 + j * a2
            P += [o, o + np.array([0, a])]
    P = np.array(P)
    P = P[(P[:, 0] > -1e-9) & (P[:, 0] < 1 + 1e-9) & (P[:, 1] > -1e-9) & (P[:, 1] < 1 + 1e-9)]
    seg = _bonds(P, a)
    # boundary: dangling (degree < 3) nodes near the sides, projected onto them
    deg = np.bincount(seg.ravel(), minlength=len(P))
    d = np.c_[P[:, 0], 1 - P[:, 0], P[:, 1], 1 - P[:, 1]]
    side = d.argmin(1); near = d.min(1) < 1.6 * a
    bnd = np.where((deg < 3) & near)[0]
    for k in bnd:
        s = side[k]
        if s == 0: P[k, 0] = 0
        elif s == 1: P[k, 0] = 1
        elif s == 2: P[k, 1] = 0
        else: P[k, 1] = 1
    used = np.unique(seg)
    idx = -np.ones(len(P), int); idx[used] = np.arange(len(used))
    bnd = idx[np.intersect1d(bnd, used)]
    return P[used], idx[seg], bnd


def voronoi(n=300, seed=0):
    """Poisson-Voronoi tessellation of the unit square (seeds mirrored in the four sides,
    so the tessellation is cut exactly by the sides, which become edges)."""
    from scipy.spatial import Voronoi
    rng = np.random.default_rng(seed)
    p = rng.random((n, 2))
    allp = np.vstack([p, p * [-1, 1], p * [-1, 1] + [2, 0], p * [1, -1], p * [1, -1] + [0, 2]])
    vor = Voronoi(allp)
    V = vor.vertices
    seg = []
    for (a, b) in vor.ridge_vertices:
        if a < 0 or b < 0:
            continue
        pa, pb = V[a], V[b]
        if np.all((pa > -1e-9) & (pa < 1 + 1e-9)) and np.all((pb > -1e-9) & (pb < 1 + 1e-9)):
            seg.append((a, b))
    seg = np.array(seg)
    P = V.copy()
    P[np.abs(P) < 1e-9] = 0; P[np.abs(P - 1) < 1e-9] = 1
    return _finish(P, seg)
