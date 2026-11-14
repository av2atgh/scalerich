"""Open-cell 3D skeleton of a plate division: struts along every plate edge.

Each plate edge (a segment on its host) is split at the junction points found by
mesh3d.segment_splits; consecutive points are joined by a straight strut that
belongs to the plate whose edge it is, with square section lam_t x lam_t,
lam_t = lam0 t^-alpha, rho = sum lam^2 l. 3D Euler-Bernoulli frame (6 dof/node),
E0 = rho0 = 1, nu0 = 0.3, J = 0.1406 lam^4, rigid joints, L_eff = max(L, lam).
C (6x6) from kinematic tests on the cube faces, Christoffel n, clamped omega1.
"""
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spl
from scipy.spatial import cKDTree
import mesh3d, shell3d

NU = 0.3


def skeleton(g):
    segs, S, split, frames = mesh3d.segment_splits(g)
    owner = []
    for s, P in enumerate(g['P']):
        owner += [s] * len(P)
    pts, A, B, own = [], [], [], []
    for k, (p0, p1) in enumerate(segs):
        L = np.linalg.norm(p1 - p0)
        u = np.unique(np.round(np.array(split[k]), 12))
        keep = [u[0]]
        for x in u[1:]:
            if (x - keep[-1]) * L > 1e-9:
                keep.append(x)
        keep[-1] = 1.0
        base = len(pts)
        pts += list(p0 + np.array(keep)[:, None] * (p1 - p0))
        for i in range(len(keep) - 1):
            A.append(base + i); B.append(base + i + 1); own.append(owner[k])
    pts = np.array(pts)
    pairs = cKDTree(pts).query_pairs(1e-8, output_type='ndarray')
    parent = np.arange(len(pts))

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]; x = parent[x]
        return x
    for i, j in pairs:
        a, b = find(i), find(j)
        if a != b:
            parent[max(a, b)] = min(a, b)
    root = np.array([find(i) for i in range(len(pts))])
    uniq, gid = np.unique(root, return_inverse=True)
    X = pts[uniq]
    seg = np.c_[gid[A], gid[B]]
    own = np.array(own)
    ok = seg[:, 0] != seg[:, 1]
    seg, own = seg[ok], own[ok]
    key = np.sort(seg, axis=1)
    _, first = np.unique(key, axis=0, return_index=True)
    return X, seg[first], own[first]


def element_k(X, seg, lam):
    d = X[seg[:, 1]] - X[seg[:, 0]]
    L = np.linalg.norm(d, axis=1); ex = d / L[:, None]
    ref = np.where(np.abs(ex[:, 2:3]) < 0.9, np.array([[0, 0, 1.]]), np.array([[1., 0, 0]]))
    ey = np.cross(ref, ex); ey /= np.linalg.norm(ey, axis=1)[:, None]
    ez = np.cross(ex, ey)
    Le = np.maximum(L, lam)
    A = lam ** 2; I = lam ** 4 / 12; J = 0.1406 * lam ** 4; G = 1 / (2 * (1 + NU))
    m = len(L)
    k = np.zeros((m, 12, 12))
    a = A / Le; t = G * J / Le
    k[:, 0, 0] = k[:, 6, 6] = a; k[:, 0, 6] = k[:, 6, 0] = -a
    k[:, 3, 3] = k[:, 9, 9] = t; k[:, 3, 9] = k[:, 9, 3] = -t
    for (v, r, sgn) in ((1, 5, 1), (2, 4, -1)):     # bending in local xy (v, theta_z) and xz (w, theta_y)
        b1 = 12 * I / Le ** 3; b2 = 6 * I / Le ** 2; b3 = 4 * I / Le; b4 = 2 * I / Le
        V = [v, r, v + 6, r + 6]
        kb = np.array([[b1, sgn * b2, -b1, sgn * b2], [sgn * b2, b3, -sgn * b2, b4],
                       [-b1, -sgn * b2, b1, -sgn * b2], [sgn * b2, b4, -sgn * b2, b3]])
        for i in range(4):
            for j in range(4):
                k[:, V[i], V[j]] = kb[i, j]
    R3 = np.stack([ex, ey, ez], 1)                   # (m,3,3) rows = local axes
    T = np.zeros((m, 12, 12))
    for o in range(4):
        T[:, 3 * o:3 * o + 3, 3 * o:3 * o + 3] = R3
    return np.einsum('mji,mjk,mkl->mil', T, k, T), L


def analyse(X, seg, lam_s, resonance=True):
    N = len(X)
    Kg, L = element_k(X, seg, lam_s)
    dof = np.c_[6 * seg[:, 0][:, None] + np.arange(6), 6 * seg[:, 1][:, None] + np.arange(6)]
    K = sp.coo_matrix((Kg.ravel(), (np.repeat(dof, 12, axis=1).ravel(), np.tile(dof, (1, 12)).ravel())),
                      shape=(6 * N, 6 * N)).tocsr()
    bnd = np.any((X < 1e-9) | (X > 1 - 1e-9), axis=1)
    fixed = np.zeros(6 * N, bool)
    for c in range(3):
        fixed[6 * np.where(bnd)[0] + c] = True
    free = ~fixed
    Kff = K[free][:, free].tocsc()
    Kff = Kff + 1e-12 * abs(Kff.diagonal()).max() * sp.eye(Kff.shape[0], format='csc')
    lu = spl.splu(Kff)
    Ua3 = shell3d.strain_fields(X)                  # (6, 3N) translational
    Ua = np.zeros((6, 6 * N))
    for c in range(3):
        Ua[:, c::6] = Ua3[:, c::3]
    Kfx = K[free][:, fixed]
    U = Ua.copy()
    for I in range(6):
        U[I, free] = lu.solve(-(Kfx @ Ua[I, fixed]))
    C = U @ (K @ U.T); C = 0.5 * (C + C.T)
    rho = float((lam_s ** 2 * L).sum())
    n = shell3d.christoffel(C, rho)
    out = dict(rho=rho, n_x=float(n[0]), n_mean=float(n[1:].mean()), n_min=float(n.min()),
               n_max=float(n.max()), eta_n=float(n.min() / n.max()), C=C.tolist(), N=N)
    if resonance:
        m = lam_s ** 2 * L
        M = np.zeros(6 * N)
        for e in (0, 1):
            for c in range(3):
                np.add.at(M, 6 * seg[:, e] + c, 0.5 * m)
                np.add.at(M, 6 * seg[:, e] + 3 + c, m * L ** 2 / 24)
        Mf = sp.diags(M[free]).tocsc()
        w2 = spl.eigsh(Kff, k=1, M=Mf, sigma=0, which='LM', return_eigenvectors=False)
        out['omega1'] = float(np.sqrt(w2.min()))
    return out


if __name__ == '__main__':
    import time, div3d
    for rule in ('SR', 'balanced'):
        g = div3d.run(300, 100, rule)
        X, seg, own = skeleton(g)
        L = np.linalg.norm(X[seg[:, 1]] - X[seg[:, 0]], axis=1)
        deg = np.bincount(seg.ravel(), minlength=len(X))
        print(rule, 'nodes', len(X), 'struts', len(seg), 'mean z(interior) %.2f' % deg[~np.any((X < 1e-9) | (X > 1 - 1e-9), axis=1)].mean())
        for rho in (0.01, 0.03, 0.1):
            base = np.ones(300); lam_t = base * np.sqrt(rho / np.sum(base[own] ** 2 * L))
            t0 = time.time(); o = analyse(X, seg, lam_t[own])
            print('  rho %.2f n_x %.2f n_mean %.2f eta %.2f w1 %.3f %.1fs' % (o['rho'], o['n_x'], o['n_mean'], o['eta_n'], o['omega1'], time.time() - t0))
