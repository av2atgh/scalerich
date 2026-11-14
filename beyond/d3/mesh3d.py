"""Conforming triangle mesh of a plate division (div3d.run output).

Every plate edge is a segment lying on its host (another plate or a cube face).
On plate t the PSLG is: t's own boundary edges + the edges of later plates
hosted by t (T-junction lines, from either side; lines from opposite sides may
cross). Split points of each segment (crossings, T-junction ends) are found in
both PSLGs it belongs to and merged, then each piece is subdivided to length ≤ h,
so host and abutting plate share the same nodes along the junction. Triangle
('pYY': no Steiner points on any segment) meshes each plate.

mesh(g, h) -> X (N,3), tri (M,3), tplate (M,)
"""
import numpy as np
import triangle as tr
from scipy.spatial import cKDTree


def _frame(P):
    e1 = P[1] - P[0]; e1 /= np.linalg.norm(e1)
    n = np.cross(P[1] - P[0], P[2] - P[0]); n /= np.linalg.norm(n)
    return P[0], e1, np.cross(n, e1), n


def segment_splits(g):
    """plate-edge segments, the PSLG membership of each plate, and split parameters."""
    P, H, T = g['P'], g['host'], g['T']
    segs = []                                   # (plate, edge) -> endpoints
    owner = {}
    S = [[] for _ in range(T)]                  # segment ids in each plate's PSLG
    for s in range(T):
        m = len(P[s])
        for i in range(m):
            k = len(segs)
            segs.append((P[s][i], P[s][(i + 1) % m]))
            S[s].append(k)
            if H[s][i] >= 0:
                S[H[s][i]].append(k)
    split = [[0.0, 1.0] for _ in segs]
    frames = [_frame(p) for p in P]
    for t in range(T):
        ids = np.array(S[t])
        o, e1, e2, _ = frames[t]
        A = np.array([segs[k][0] for k in ids]) - o
        B = np.array([segs[k][1] for k in ids]) - o
        a = np.c_[A @ e1, A @ e2]; b = np.c_[B @ e1, B @ e2]
        r = b - a
        # pairwise: a_i + u r_i = a_j + v r_j
        den = r[:, None, 0] * r[None, :, 1] - r[:, None, 1] * r[None, :, 0]
        q = a[None, :, :] - a[:, None, :]
        with np.errstate(divide='ignore', invalid='ignore'):
            u = (q[..., 0] * r[None, :, 1] - q[..., 1] * r[None, :, 0]) / den
            v = (q[..., 0] * r[:, None, 1] - q[..., 1] * r[:, None, 0]) / den
        tol = 1e-7
        ok = (np.abs(den) > 1e-14) & (u > -tol) & (u < 1 + tol) & (v > -tol) & (v < 1 + tol)
        I, J = np.nonzero(ok)
        for i_, j_ in zip(I, J):
            split[ids[i_]].append(float(np.clip(u[i_, j_], 0, 1)))
    return segs, S, split, frames


def mesh(g, h=0.05, qmin=25):
    P, H, T = g['P'], g['host'], g['T']
    segs, S, split, frames = segment_splits(g)
    # subdivide
    pts = []; seg_pts = []
    for k, (p0, p1) in enumerate(segs):
        L = np.linalg.norm(p1 - p0)
        s = np.unique(np.round(np.array(split[k]), 12))
        keep = [s[0]]
        for x in s[1:]:
            if (x - keep[-1]) * L > 1e-9:
                keep.append(x)
        keep[-1] = 1.0
        fine = [keep[0]]
        for x0, x1 in zip(keep[:-1], keep[1:]):
            nd = max(1, int(np.ceil((x1 - x0) * L / h)))
            fine += list(x0 + (x1 - x0) * np.arange(1, nd + 1) / nd)
        fine = np.array(fine)
        seg_pts.append(np.arange(len(pts), len(pts) + len(fine)))
        pts += list(p0 + fine[:, None] * (p1 - p0))
    pts = np.array(pts)
    # merge coincident points
    tree = cKDTree(pts)
    pairs = tree.query_pairs(1e-8, output_type='ndarray')
    parent = np.arange(len(pts))

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]; x = parent[x]
        return x
    for i, j in pairs:
        ri, rj = find(i), find(j)
        if ri != rj:
            parent[max(ri, rj)] = min(ri, rj)
    root = np.array([find(i) for i in range(len(pts))])
    uniq, gid = np.unique(root, return_inverse=True)
    X = list(pts[uniq])
    tri_all = []; tplate = []
    for t in range(T):
        o, e1, e2, _ = frames[t]
        vids = []; sg = []
        for k in S[t]:
            g_ = gid[seg_pts[k]]
            vids += list(g_)
            sg += [(a_, b_) for a_, b_ in zip(g_[:-1], g_[1:]) if a_ != b_]
        vids = np.unique(vids)
        loc = {v: i for i, v in enumerate(vids)}
        Q = np.array([X[v] for v in vids]) - o
        xy = np.c_[Q @ e1, Q @ e2]
        sgl = np.unique(np.sort(np.array([(loc[a_], loc[b_]) for a_, b_ in sg]), axis=1), axis=0)
        area = 0.5 * abs(np.cross(P[t][1:-1] - P[t][0], P[t][2:] - P[t][0]).sum(0) @ frames[t][3])
        out = tr.triangulate({'vertices': xy, 'segments': sgl}, f'pYYQq{qmin}a{0.5 * h * h:.8f}')
        V = out['vertices']; nin = len(xy)
        mp = np.empty(len(V), int); mp[:nin] = vids
        for k in range(nin, len(V)):
            mp[k] = len(X); X.append(o + V[k, 0] * e1 + V[k, 1] * e2)
        tri_all.append(mp[out['triangles']]); tplate.append(np.full(len(out['triangles']), t))
    return np.array(X), np.vstack(tri_all), np.concatenate(tplate)


if __name__ == '__main__':
    import time, div3d
    for r in div3d.RULES:
        g = div3d.run(200, 1, r)
        t0 = time.time()
        X, tri, tp = mesh(g, 0.05)
        A = 0.5 * np.linalg.norm(np.cross(X[tri[:, 1]] - X[tri[:, 0]], X[tri[:, 2]] - X[tri[:, 0]]), axis=1)
        area = sum(0.5 * np.linalg.norm(np.cross(p[1:-1] - p[0], p[2:] - p[0]).sum(0)) for p in g['P'])
        # conformity: every mesh edge used by 1 triangle must lie on the cube boundary
        E = np.sort(np.vstack([tri[:, [0, 1]], tri[:, [1, 2]], tri[:, [2, 0]]]), axis=1)
        ue, cnt = np.unique(E, axis=0, return_counts=True)
        lone = ue[cnt == 1]
        onb = lambda x: np.any((x < 1e-9) | (x > 1 - 1e-9), axis=-1)
        bad = ~(onb(X[lone[:, 0]]) & onb(X[lone[:, 1]]))
        print(r, 'N', len(X), 'M', len(tri), 'area %.4f vs %.4f' % (A.sum(), area),
              'open interior edges', int(bad.sum()), 'max edge mult', cnt.max(), '%.1fs' % (time.time() - t0))
