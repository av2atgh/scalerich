"""SR model as a nested binary hierarchy of convex cells (zero-width topology).

The SR step "uniform point + uniform angle, grow both ways to the first hit" is
exactly: locate the convex cell containing the point (descend the splitting
tree: each internal node is a chord = half-plane test) and cut it with the
chord through the point.  Ligament t (t = 0..T-1, birth index t+1) has
thickness lam0 * (t+1)**(-alpha); thickness paints, it does not change the
topology (README: topology independent of (alpha, lam0) in the tunable phase).

Optional mode='rect' restricts angles to {0, pi/2} (Krapivsky-Ben-Naim
rectangular cracks, area-weighted) as a control for multiscaling.

run(T, seed) -> dict of arrays:
  A, B        (T,2) chord endpoints
  hostA/B     (T,) ligament id the endpoint lies on (-1..-4: square side
              bottom,right,top,left)
  parent      (T,) cell id split at step t;  kids (T,2) child cell ids
  depth       (T,) depth of the split cell in the tree (root = 0)
  area        (T,) area of the split cell; r (T,) area fraction of child 0
  ncorner     (T,) corners of the split cell
  cell_*      per cell: area, depth, ncorner, born (step created, -1 root),
              leaf flag (unsplit at T)
"""
import math
import numpy as np


def _clip(poly, p, d):
    """Intersections of line p + s d with convex CCW polygon poly (m,2).
    Returns (i, sP, j, sQ): edge i hit at s<0, edge j at s>0."""
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


def _area(poly):
    x, y = poly[:, 0], poly[:, 1]
    return 0.5 * float(np.dot(x, np.roll(y, -1)) - np.dot(np.roll(x, -1), y))


def run(T, seed=0, mode='iso'):
    rng = np.random.default_rng(seed)
    polys = [np.array([[0., 0.], [1., 0.], [1., 1.], [0., 1.]])]
    labs = [[-1, -2, -3, -4]]
    cdepth = [0]; cborn = [-1]; carea = [1.0]
    # tree: for internal cell c: chord point, normal, kids
    split_of = {}  # cell -> (point, normal, kid0, kid1)
    A = np.zeros((T, 2)); B = np.zeros((T, 2))
    hostA = np.zeros(T, int); hostB = np.zeros(T, int)
    parent = np.zeros(T, int); kids = np.zeros((T, 2), int)
    depth = np.zeros(T, int); area = np.zeros(T); r = np.zeros(T)
    ncorner = np.zeros(T, int)
    for t in range(T):
        p = rng.random(2)
        c = 0
        while c in split_of:
            q, nrm, k0, k1 = split_of[c]
            c = k0 if (p - q) @ nrm >= 0 else k1
        if mode == 'rect':
            th = 0.5 * np.pi * rng.integers(2)
        else:
            th = rng.random() * np.pi
        d = np.array([math.cos(th), math.sin(th)])
        poly = polys[c]; lab = labs[c]; m = len(poly)
        i, sP, j, sQ = _clip(poly, p, d)
        P = p + sP * d; Q = p + sQ * d
        # child 0: P, v_{i+1}..v_j, Q  (left of P->Q? it's the side reached from P going CCW)
        idx0 = [(i + 1 + k) % m for k in range((j - i) % m)]
        p0 = np.vstack([P] + [poly[k] for k in idx0] + [Q])
        # edges P->v(i+1) [lab i], v(k)->v(k+1) [lab k], ..., v(j)->Q [lab j], Q->P [t]
        l0 = [lab[i]] + [lab[k] for k in idx0] + [t]
        idx1 = [(j + 1 + k) % m for k in range((i - j) % m)]
        p1 = np.vstack([Q] + [poly[k] for k in idx1] + [P])
        l1 = [lab[j]] + [lab[k] for k in idx1] + [t]
        a0 = _area(p0); a1 = _area(p1)
        n0 = len(polys); n1 = n0 + 1
        polys += [p0, p1]; labs += [l0, l1]
        cdepth += [cdepth[c] + 1] * 2; cborn += [t, t]; carea += [a0, a1]
        # normal pointing into child 0
        nrm = np.array([-d[1], d[0]])
        cen0 = p0.mean(axis=0)
        if (cen0 - p) @ nrm < 0:
            nrm = -nrm
        split_of[c] = (p.copy(), nrm, n0, n1)
        A[t] = P; B[t] = Q; hostA[t] = lab[i]; hostB[t] = lab[j]
        parent[t] = c; kids[t] = (n0, n1); depth[t] = cdepth[c]
        area[t] = carea[c]; r[t] = a0 / carea[c]; ncorner[t] = m
    ncell = len(polys)
    leaf = np.ones(ncell, bool)
    for c in split_of:
        leaf[c] = False
    return dict(T=T, seed=seed, A=A, B=B, hostA=hostA, hostB=hostB,
                parent=parent, kids=kids, depth=depth, area=area, r=r,
                ncorner=ncorner, cell_area=np.array(carea),
                cell_depth=np.array(cdepth), cell_born=np.array(cborn),
                cell_ncorner=np.array([len(x) for x in polys]), leaf=leaf,
                polys=polys, labs=labs)


def thickness(T, lam0, alpha):
    return lam0 * np.arange(1, T + 1, dtype=float) ** (-alpha)


def network(g, frame=None):
    """Segment network.  Node 2t = endpoint A of ligament t, 2t+1 = endpoint B.
    A junction node is the endpoint of the stem ligament, lying inside the host.
    Returns nodes (N,2), seg (S,2) node pairs, seg_lig (S,) owning ligament,
    bnode: dict side(-1..-4) -> node ids on that side.
    frame: if not None, the square's sides are added as ligaments with ids
    T..T+3 (sides -1..-4) and 4 corner nodes 2T..2T+3."""
    T = g['T']
    nodes = np.empty((2 * T, 2))
    nodes[0::2] = g['A']; nodes[1::2] = g['B']
    on = {}  # host -> list of node ids lying on it
    for t in range(T):
        for k, h in ((2 * t, g['hostA'][t]), (2 * t + 1, g['hostB'][t])):
            on.setdefault(int(h), []).append(k)
    seg = []; seg_lig = []
    for t in range(T):
        ids = [2 * t, 2 * t + 1] + on.get(t, [])
        a = g['A'][t]; d = g['B'][t] - a
        s = (nodes[ids] - a) @ d
        o = np.argsort(s)
        ids = [ids[k] for k in o]
        for u, v in zip(ids[:-1], ids[1:]):
            seg.append((u, v)); seg_lig.append(t)
    bnode = {sd: on.get(sd, []) for sd in (-1, -2, -3, -4)}
    if frame is not None:
        corners = np.array([[0., 0.], [1., 0.], [1., 1.], [0., 1.]])
        nodes = np.vstack([nodes, corners]); c0 = 2 * T
        # side -1 bottom: corner0->corner1, -2 right: c1->c2, -3 top: c2->c3, -4 left: c3->c0
        for k, sd in enumerate((-1, -2, -3, -4)):
            ca, cb = c0 + k, c0 + (k + 1) % 4
            ids = [ca, cb] + on.get(sd, [])
            a = nodes[ca]; d = nodes[cb] - a
            s = (nodes[ids] - a) @ d
            ids = [ids[q] for q in np.argsort(s)]
            for u, v in zip(ids[:-1], ids[1:]):
                seg.append((u, v)); seg_lig.append(T + k)
            bnode[sd] = ids
    return nodes, np.array(seg), np.array(seg_lig), bnode


if __name__ == '__main__':
    import sys, time
    t0 = time.time()
    g = run(int(sys.argv[1]) if len(sys.argv) > 1 else 2000, 1)
    print('time', time.time() - t0)
    la = g['cell_area'][g['leaf']]
    print('leaf area sum', la.sum(), 'ncells', g['leaf'].sum())
    nodes, seg, sl, bn = network(g)
    print('nodes', len(nodes), 'segs', len(seg))
    L = np.linalg.norm(g['B'] - g['A'], axis=1)
    print('mean first ligament length check (single realisation)', L[0])
