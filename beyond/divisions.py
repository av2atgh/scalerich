"""Recursive divisions of the unit square (zero-width topology).

A division rule has three parts, applied at every step t = 0..T-1:
  pick   which cell is cut:  'area'    cell chosen with probability ∝ area (SR)
                             'uniform' cell chosen uniformly among the leaves
                             'largest' the largest leaf (deterministic order)
  where  the cut passes:     'point'   a uniform point of the chosen cell
                             'centroid' the centroid of the chosen cell
  angle  of the cut:         'uniform' θ uniform in [0, π)
                             'rect'    θ ∈ {0, π/2}
                             'perp'    perpendicular to the cell's major axis
                                       (principal axis of its second moment)
SR (Both et al.) = (area, point, uniform). Cowan's iterated cell division
covers the whole family. Ligament t has thickness λ0 (t+1)^(-α) (painting by
birth order, as in SR). Output matches rg/srtree.run so rg/srtree.network works.
"""
import math
import heapq
import numpy as np

RULES = {
    'SR':        ('area', 'point', 'uniform'),
    'uniform':   ('uniform', 'point', 'uniform'),
    'rect':      ('area', 'point', 'rect'),
    'centroid':  ('area', 'centroid', 'uniform'),
    'balanced':  ('area', 'centroid', 'perp'),
    'regular':   ('largest', 'centroid', 'perp'),
}


def _clip(poly, p, d):
    v = poly
    e = np.roll(v, -1, axis=0) - v
    w = v - p
    den = d[0] * e[:, 1] - d[1] * e[:, 0]
    with np.errstate(divide='ignore', invalid='ignore'):
        s = (w[:, 0] * e[:, 1] - w[:, 1] * e[:, 0]) / den
        u = (w[:, 0] * d[1] - w[:, 1] * d[0]) / den
    ok = (u >= 0) & (u < 1) & np.isfinite(s)
    sn = np.where(ok & (s < 0), s, -np.inf)
    sp = np.where(ok & (s > 0), s, np.inf)
    i = int(np.argmax(sn)); j = int(np.argmin(sp))
    return i, sn[i], j, sp[j]


def _moments(poly):
    """area, centroid, second-moment tensor about the centroid (fan triangulation)."""
    x, y = poly[:, 0], poly[:, 1]
    x1, y1 = np.roll(x, -1), np.roll(y, -1)
    cr = x * y1 - x1 * y
    A = 0.5 * cr.sum()
    cx = ((x + x1) * cr).sum() / (6 * A); cy = ((y + y1) * cr).sum() / (6 * A)
    Ixx = ((y * y + y * y1 + y1 * y1) * cr).sum() / 12 - A * cy * cy
    Iyy = ((x * x + x * x1 + x1 * x1) * cr).sum() / 12 - A * cx * cx
    Ixy = ((x * y1 + 2 * x * y + 2 * x1 * y1 + x1 * y) * cr).sum() / 24 - A * cx * cy
    return A, np.array([cx, cy]), np.array([[Iyy, Ixy], [Ixy, Ixx]])


def _sample(poly, rng):
    p0 = poly[0]
    a = poly[1:-1] - p0; b = poly[2:] - p0
    w = 0.5 * np.abs(a[:, 0] * b[:, 1] - a[:, 1] * b[:, 0])
    k = rng.choice(len(w), p=w / w.sum())
    r1, r2 = rng.random(2)
    if r1 + r2 > 1:
        r1, r2 = 1 - r1, 1 - r2
    return p0 + r1 * a[k] + r2 * b[k]


def run(T, seed=0, rule='SR', box=(1.0, 1.0)):
    pick, where, angle = RULES[rule] if isinstance(rule, str) else rule
    rng = np.random.default_rng(seed)
    Wb, Hb = box
    polys = [np.array([[0., 0.], [Wb, 0.], [Wb, Hb], [0., Hb]])]
    labs = [[-1, -2, -3, -4]]
    areas = [1.0]
    leaves = [0]; pos = {0: 0}           # uniform pick: O(1) removal
    heap = [(-1.0, 0)]                   # largest pick
    cum = None
    A = np.zeros((T, 2)); B = np.zeros((T, 2))
    hostA = np.zeros(T, int); hostB = np.zeros(T, int)
    split_of = {}
    depth = {0: 0}
    pts = np.zeros((T, 2)); cells = []
    frac = np.zeros(T); cut_depth = np.zeros(T, int); parent_n = np.zeros(T, int)
    for t in range(T):
        # --- pick a cell
        if pick == 'area':
            p = rng.random(2) * np.array([Wb, Hb]); c = 0
            while c in split_of:
                q, nrm, k0, k1 = split_of[c]
                c = k0 if (p - q) @ nrm >= 0 else k1
        elif pick == 'uniform':
            c = leaves[rng.integers(len(leaves))]
            p = _sample(polys[c], rng)
        else:
            _, c = heapq.heappop(heap)
            p = _sample(polys[c], rng)
        poly = polys[c]; lab = labs[c]; m = len(poly)
        if where == 'centroid' or angle == 'perp' or isinstance(angle, tuple):
            _, cen, I2 = _moments(poly)
            if where == 'centroid':
                p = cen
        if angle == 'uniform':
            th = rng.random() * np.pi
        elif angle == 'rect':
            th = 0.5 * np.pi * rng.integers(2)
        elif isinstance(angle, tuple) and angle[0] == 'perpmix_field':
            w, V = np.linalg.eigh(I2)
            ax = V[:, 1]
            kap = float(angle[1](p[0], p[1]))
            th = math.atan2(ax[1], ax[0]) + 0.5 * np.pi + kap * np.pi * (rng.random() - 0.5) \
                + 0.02 * rng.standard_normal()
        elif isinstance(angle, tuple) and angle[0] == 'perpmix':
            # perpendicular to the major axis plus a uniform tilt in [-kappa pi/2, kappa pi/2]:
            # kappa = 0 -> 'balanced', kappa = 1 -> uniform angle ('centroid')
            w, V = np.linalg.eigh(I2)
            ax = V[:, 1]
            th = math.atan2(ax[1], ax[0]) + 0.5 * np.pi + angle[1] * np.pi * (rng.random() - 0.5) \
                + 0.02 * rng.standard_normal()
        else:
            w, V = np.linalg.eigh(I2)          # I2 ~ ∫(x-c)(x-c)^T, largest eig = major axis
            ax = V[:, 1]
            th = math.atan2(ax[1], ax[0]) + 0.5 * np.pi + 0.02 * rng.standard_normal()
        d = np.array([math.cos(th), math.sin(th)])
        i, sP, j, sQ = _clip(poly, p, d)
        P = p + sP * d; Q = p + sQ * d
        idx0 = [(i + 1 + k) % m for k in range((j - i) % m)]
        p0 = np.vstack([P] + [poly[k] for k in idx0] + [Q])
        l0 = [lab[i]] + [lab[k] for k in idx0] + [t]
        idx1 = [(j + 1 + k) % m for k in range((i - j) % m)]
        p1 = np.vstack([Q] + [poly[k] for k in idx1] + [P])
        l1 = [lab[j]] + [lab[k] for k in idx1] + [t]
        a0 = _moments(p0)[0]; a1 = _moments(p1)[0]
        n0 = len(polys); n1 = n0 + 1
        frac[t] = a0 / (a0 + a1); cut_depth[t] = depth[c]; parent_n[t] = m
        depth[n0] = depth[n1] = depth[c] + 1
        polys += [p0, p1]; labs += [l0, l1]; areas += [a0, a1]
        nrm = np.array([-d[1], d[0]])
        if (p0.mean(axis=0) - p) @ nrm < 0:
            nrm = -nrm
        split_of[c] = (p.copy(), nrm, n0, n1)
        # leaf bookkeeping
        k = pos.pop(c); last = leaves[-1]; leaves[k] = last; leaves.pop()
        if last != c:
            pos[last] = k
        for n_ in (n0, n1):
            pos[n_] = len(leaves); leaves.append(n_)
            heapq.heappush(heap, (-areas[n_], n_))
        pts[t] = p; cells.append(poly)
        A[t] = P; B[t] = Q; hostA[t] = lab[i]; hostB[t] = lab[j]
    return dict(T=T, seed=seed, A=A, B=B, hostA=hostA, hostB=hostB,
                leaf_area=np.array([areas[c] for c in leaves]),
                leaf_polys=[polys[c] for c in leaves],
                leaf_depth=np.array([depth[c] for c in leaves]),
                frac=frac, cut_depth=cut_depth, parent_n=parent_n,
                point=pts, cell=cells)
