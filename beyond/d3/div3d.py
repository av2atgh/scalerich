"""Recursive divisions of the unit cube by planar plates (zero-thickness topology).

Rule = (pick, where, normal):
  pick    'volume'  cell ∝ volume (SR: uniform point in the cube, then its cell)
          'uniform' cell uniform among leaves
          'largest' largest leaf first
  where   'point'   uniform point of the cell;  'centroid'
  normal  'iso'     isotropic unit normal;  'axes' random coordinate axis
          'perp'    major axis of the cell's inertia (cut ⟂ to the longest direction),
                    with 0.02 rad noise to break exact coincidences
SR in 3D (Both et al. Fig. 6a(v); d3/sr3d.c) = (volume, point, iso).

Plate t: convex polygon P[t] (m,3), ordered; host[t][i] = label of the cell face
that edge i (P[i] -> P[i+1]) lies on: a plate id >= 0, or a cube face -1..-6
(x=0, x=1, y=0, y=1, z=0, z=1).
"""
import numpy as np

RULES = {
    'SR':       ('volume', 'point', 'iso'),
    'uniform':  ('uniform', 'point', 'iso'),
    'rect':     ('volume', 'point', 'axes'),
    'centroid': ('volume', 'centroid', 'iso'),
    'balanced': ('volume', 'centroid', 'perp'),
    'regular':  ('largest', 'centroid', 'perp'),
}


def cube():
    V = np.array([[x, y, z] for x in (0., 1.) for y in (0., 1.) for z in (0., 1.)])
    idx = lambda x, y, z: 4 * x + 2 * y + z
    F = [[idx(0, 0, 0), idx(0, 0, 1), idx(0, 1, 1), idx(0, 1, 0)],   # x=0
         [idx(1, 0, 0), idx(1, 1, 0), idx(1, 1, 1), idx(1, 0, 1)],   # x=1
         [idx(0, 0, 0), idx(1, 0, 0), idx(1, 0, 1), idx(0, 0, 1)],   # y=0
         [idx(0, 1, 0), idx(0, 1, 1), idx(1, 1, 1), idx(1, 1, 0)],   # y=1
         [idx(0, 0, 0), idx(0, 1, 0), idx(1, 1, 0), idx(1, 0, 0)],   # z=0
         [idx(0, 0, 1), idx(1, 0, 1), idx(1, 1, 1), idx(0, 1, 1)]]   # z=1
    return V, F, [-1, -2, -3, -4, -5, -6]


def tets(V, F):
    """fan tetrahedra from the vertex mean: returns (k,4,3)"""
    c = V.mean(axis=0)
    T = []
    for f in F:
        for a, b in zip(f[1:-1], f[2:]):
            T.append((c, V[f[0]], V[a], V[b]))
    return np.array(T)


def moments(V, F):
    T = tets(V, F)
    e = T[:, 1:] - T[:, :1]
    vol = np.abs(np.einsum('ij,ij->i', e[:, 0], np.cross(e[:, 1], e[:, 2]))) / 6
    cen_t = T.mean(axis=1)
    Vt = vol.sum()
    cen = (vol[:, None] * cen_t).sum(0) / Vt
    r = T - cen
    S = r.sum(axis=1)
    cov = (vol[:, None, None] / 20) * (np.einsum('kai,kaj->kij', r, r) + np.einsum('ki,kj->kij', S, S))
    return Vt, cen, cov.sum(0), T, vol


def sample(T, vol, rng):
    k = rng.choice(len(vol), p=vol / vol.sum())
    w = rng.dirichlet(np.ones(4))
    return w @ T[k]


def clip(V, F, L, n, c, newlab):
    """split polyhedron by n.x = c. Returns (V1,F1,L1), (V2,F2,L2), section polygon (m,3), hosts."""
    d = V @ n - c
    nv = len(V)
    Vn = [v for v in V]
    cut = {}

    def xpt(a, b):
        k = (min(a, b), max(a, b))
        if k not in cut:
            s = d[a] / (d[a] - d[b])
            Vn.append(V[a] + s * (V[b] - V[a])); cut[k] = len(Vn) - 1
        return cut[k]

    Fp, Lp, Fm, Lm = [], [], [], []
    sec_edges = []                    # (i, j, host label)
    for f, lab in zip(F, L):
        pp, mm, xs = [], [], []
        m = len(f)
        for q in range(m):
            a, b = f[q], f[(q + 1) % m]
            if d[a] >= 0:
                pp.append(a)
            else:
                mm.append(a)
            if (d[a] >= 0) != (d[b] >= 0):
                x = xpt(a, b); pp.append(x); mm.append(x); xs.append(x)
        if len(pp) >= 3:
            Fp.append(pp); Lp.append(lab)
        if len(mm) >= 3:
            Fm.append(mm); Lm.append(lab)
        if len(xs) == 2:
            sec_edges.append((xs[0], xs[1], lab))
    Vn = np.array(Vn)
    # order the section polygon by chaining edges
    adj = {}
    for i, j, lab in sec_edges:
        adj.setdefault(i, []).append((j, lab)); adj.setdefault(j, []).append((i, lab))
    start = sec_edges[0][0]
    poly = [start]; hosts = []
    prev = None; cur = start
    while True:
        nxt = [(j, lab) for j, lab in adj[cur] if j != prev]
        if prev is None:
            nxt = adj[cur][:1]
        j, lab = nxt[0]
        hosts.append(lab)
        if j == start:
            break
        poly.append(j); prev, cur = cur, j
    poly = list(poly)
    P = Vn[poly]
    # orient so the polygon normal is +n
    nrm = np.cross(P[1] - P[0], P[2] - P[0])
    if nrm @ n < 0:
        poly = poly[::-1]
    emap = {}
    for i, j, lab in sec_edges:
        emap[(i, j)] = lab; emap[(j, i)] = lab
    hosts = [emap[(poly[q], poly[(q + 1) % len(poly)])] for q in range(len(poly))]
    Fp.append(poly[::-1]); Lp.append(newlab)     # outward normal of + cell is -n
    Fm.append(poly); Lm.append(newlab)

    def compact(Fs):
        used = sorted({v for f in Fs for v in f})
        mp = {v: k for k, v in enumerate(used)}
        return Vn[used], [[mp[v] for v in f] for f in Fs]
    V1, F1 = compact(Fp); V2, F2 = compact(Fm)
    return (V1, F1, Lp), (V2, F2, Lm), Vn[poly], hosts


def run(T, seed=0, rule='SR'):
    pick, where, normal = RULES[rule] if isinstance(rule, str) else rule
    rng = np.random.default_rng(seed)
    V, F, L = cube()
    cells = [(V, F, L)]
    vols = [1.0]
    P, H, N = [], [], []
    for t in range(T):
        va = np.array(vols)
        if pick == 'volume':
            c = int(rng.choice(len(cells), p=va / va.sum()))
        elif pick == 'uniform':
            c = int(rng.integers(len(cells)))
        else:
            c = int(np.argmax(va))
        V, F, L = cells[c]
        vol, cen, cov, Tt, tv = moments(V, F)
        p = cen if where == 'centroid' else sample(Tt, tv, rng)
        if normal == 'iso':
            n = rng.standard_normal(3); n /= np.linalg.norm(n)
        elif normal == 'axes':
            n = np.eye(3)[rng.integers(3)]
        else:
            w, Q = np.linalg.eigh(cov)
            n = Q[:, 2] + 0.02 * rng.standard_normal(3); n /= np.linalg.norm(n)
        A, B, poly, hosts = clip(V, F, L, n, float(n @ p), t)
        cells[c] = A; cells.append(B)
        vols[c] = moments(*A[:2])[0]; vols.append(moments(*B[:2])[0])
        P.append(poly); H.append(hosts); N.append(n)
    return dict(T=T, P=P, host=H, normal=np.array(N), vols=np.array(vols))


if __name__ == '__main__':
    import time
    for r in RULES:
        t0 = time.time(); g = run(300, 1, r)
        area = sum(0.5 * np.linalg.norm(np.cross(p[1:-1] - p[0], p[2:] - p[0]), axis=1).sum() for p in g['P'])
        print(r, 'vol sum %.6f' % g['vols'].sum(), 'area %.2f' % area,
              'mean m %.2f' % np.mean([len(p) for p in g['P']]), '%.1fs' % (time.time() - t0))
