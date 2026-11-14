"""Framed periodic tiling of 3D plate divisions (2x2x2 mirror supercell).

The six cube faces become plates (generation 0, thickness lam0) shared between
neighbouring tiles; the plate division is meshed with them (mesh3d), mirrored in
x = 1, y = 1, z = 1, coincident nodes merged, duplicate triangles on the internal
mirror planes removed, and the faces on x, y, z = 2 dropped (images of x, y, z = 0).
Opposite faces of the supercell are paired; fluctuations are periodic, one node pinned.
C per unit volume (supercell volume 8). -> compare with the kinematic shell3d result.
"""
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spl
from scipy.spatial import cKDTree
import mesh3d, shell3d

CUBE = [np.array(p, float) for p in (
    [[0, 0, 0], [0, 0, 1], [0, 1, 1], [0, 1, 0]], [[1, 0, 0], [1, 1, 0], [1, 1, 1], [1, 0, 1]],
    [[0, 0, 0], [1, 0, 0], [1, 0, 1], [0, 0, 1]], [[0, 1, 0], [0, 1, 1], [1, 1, 1], [1, 1, 0]],
    [[0, 0, 0], [0, 1, 0], [1, 1, 0], [1, 0, 0]], [[0, 0, 1], [1, 0, 1], [1, 1, 1], [0, 1, 1]])]


def with_frame(g):
    T = g['T']
    P = list(g['P']) + CUBE
    H = [[h if h >= 0 else T + (-h - 1) for h in hs] for hs in g['host']] + [[-1] * 4 for _ in range(6)]
    return dict(T=T + 6, P=P, host=H)


def _merge(X, tol=1e-9):
    pairs = cKDTree(X).query_pairs(tol, output_type='ndarray')
    parent = np.arange(len(X))

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]; x = parent[x]
        return x
    for i, j in pairs:
        a, b = find(i), find(j)
        if a != b:
            parent[max(a, b)] = min(a, b)
    root = np.array([find(i) for i in range(len(X))])
    uniq, gid = np.unique(root, return_inverse=True)
    return X[uniq], gid


def supercell(X, tri, lam_tri):
    Xs, Ts, Ls = [], [], []
    N = len(X)
    k = 0
    for sx in (0, 1):
        for sy in (0, 1):
            for sz in (0, 1):
                x = X.copy()
                for ax, sgn in enumerate((sx, sy, sz)):
                    if sgn:
                        x[:, ax] = 2 - x[:, ax]
                Xs.append(x); Ts.append(tri + k * N); Ls.append(lam_tri); k += 1
    Xa = np.vstack(Xs); Ta = np.vstack(Ts); La = np.concatenate(Ls)
    Xm, gid = _merge(Xa)
    Tm = gid[Ta]
    key = np.sort(Tm, axis=1)
    _, first = np.unique(key, axis=0, return_index=True)
    keep = np.zeros(len(Tm), bool); keep[first] = True
    for ax in range(3):
        keep &= ~np.all(np.abs(Xm[Tm][:, :, ax] - 2) < 1e-9, axis=1)
    return Xm, Tm[keep], La[keep]


def pairing(X, P=2.0, tol=1e-9):
    N = len(X)
    master = np.arange(N)
    for ax in range(3):
        lo = np.where(np.abs(X[:, ax]) < tol)[0]; hi = np.where(np.abs(X[:, ax] - P) < tol)[0]
        o = [c for c in range(3) if c != ax]
        d, idx = cKDTree(X[lo][:, o]).query(X[hi][:, o])
        assert np.all(d < 1e-7), 'unmatched periodic nodes %g' % d.max()
        master[hi] = lo[idx]
    for i in range(N):
        while master[master[i]] != master[i]:
            master[i] = master[master[i]]
    return master


def periodic_C(X, tri, lam_tri):
    K, A = shell3d.assemble(X, tri, lam_tri)
    master = pairing(X)
    N = len(X)
    ms = np.unique(master)
    red = -np.ones(N, int); red[ms] = np.arange(len(ms))
    P = sp.csr_matrix((np.ones(3 * N), (np.arange(3 * N), (3 * red[master][:, None] + np.arange(3)).ravel())),
                      shape=(3 * N, 3 * len(ms)))
    Kr = (P.T @ K @ P).tocsc()
    n = Kr.shape[0]
    keep = np.ones(n, bool); keep[:3] = False
    Kk = Kr[keep][:, keep]
    Kk = Kk + 1e-12 * abs(Kk.diagonal()).max() * sp.eye(Kk.shape[0], format='csc')
    lu = spl.splu(Kk.tocsc())
    Ua = shell3d.strain_fields(X - 0.0)
    U = []
    for I in range(6):
        rhs = -(P.T @ (K @ Ua[I]))
        w = np.zeros(n); w[keep] = lu.solve(rhs[keep])
        U.append(Ua[I] + P @ w)
    U = np.array(U)
    C = U @ (K @ U.T) / 8.0
    rho = float((lam_tri * A).sum() / 8.0)
    return 0.5 * (C + C.T), rho


def symmetric_normal(X, tri, lam_tri):
    """Framed cell with symmetry BCs (normal displacement prescribed on each cube face,
    tangential free): equivalent to the mirror-periodic tiling for normal strains.
    Returns C11, C22, C33 (uniaxial strain) and rho."""
    K, A = shell3d.assemble(X, tri, lam_tri)
    N = len(X)
    on = [(np.abs(X[:, ax]) < 1e-9) | (np.abs(X[:, ax] - 1) < 1e-9) for ax in range(3)]
    fixed = np.zeros(3 * N, bool)
    for ax in range(3):
        fixed[3 * np.where(on[ax])[0] + ax] = True
    free = ~fixed
    Kff = K[free][:, free].tocsc()
    Kff = Kff + 1e-12 * abs(Kff.diagonal()).max() * sp.eye(Kff.shape[0], format='csc')
    lu = spl.splu(Kff)
    Kfx = K[free][:, fixed]
    out = []
    for ax in range(3):
        u = np.zeros(3 * N)
        u[ax::3] = X[:, ax] - 0.5          # uniaxial strain along ax; other normal comps 0 on faces
        uf = u.copy(); uf[~fixed] = 0
        uf[free] = lu.solve(-(Kfx @ uf[fixed]))
        out.append(float(uf @ (K @ uf)))
    rho = float((lam_tri * A).sum())
    return np.array(out), rho
