"""Periodic unit cells by mirror tiling, periodic homogenisation and Bloch dispersion (2D).

The 2x2 mirror image of a square cell (reflections in x = 1 and y = 1) is exactly
periodic with period 2: nodes on x = 0 and x = 2 (y = 0 and y = 2) coincide in y (x).
Mirror symmetry makes the supercell orthotropic (C16 = C26 = 0).

periodic_C(...)  : u = eps.x + w, w periodic (translations and rotations), one node pinned.
                   C per unit area (area 4).  n(theta) via Christoffel.
bloch_x(...)     : Bloch waves along x, u_slave = u_master exp(i k (x_s - x_m)); for each k
                   the modes nearest a guessed frequency, the longitudinal one picked by its
                   plane-wave projection p_L = |sum m u_x e^{-ikx}|^2 / (sum m)(sum m|u|^2).
                   Gives n(omega) = k c0 / omega and the homogenisation edge.
Units: E0 = rho0 = 1 (c0 = 1), lengths in units of the original cell L.
"""
import os, sys
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spl
from scipy.spatial import cKDTree

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', 'rg'))
import srtree, frame  # noqa: E402
import waves  # noqa: E402


def mirror_supercell(nodes, seg, lam_s):
    Ns = len(nodes)
    X, S, Lm = [], [], []
    for k, (sx, sy) in enumerate(((0, 0), (1, 0), (0, 1), (1, 1))):
        x = nodes.copy()
        if sx: x[:, 0] = 2 - x[:, 0]
        if sy: x[:, 1] = 2 - x[:, 1]
        X.append(x); S.append(seg + k * Ns); Lm.append(lam_s)
    X = np.vstack(X); S = np.vstack(S); Lm = np.concatenate(Lm)
    pairs = cKDTree(X).query_pairs(1e-10, output_type='ndarray')
    parent = np.arange(len(X))
    for i, j in pairs:
        a, b = parent[i], parent[j]
        while parent[a] != a: a = parent[a]
        while parent[b] != b: b = parent[b]
        if a != b: parent[max(a, b)] = min(a, b)
    for i in range(len(X)):
        r = i
        while parent[r] != r: r = parent[r]
        parent[i] = r
    uniq, gid = np.unique(parent, return_inverse=True)
    return X[uniq], gid[S], Lm


def pairing(X, tol=1e-9, P=2.0):
    """slave -> (master, shift vector).  x = P paired with x = 0, y = P with y = 0."""
    N = len(X)
    master = np.arange(N); shift = np.zeros((N, 2))
    for ax in (0, 1):
        lo = np.where(np.abs(X[:, ax]) < tol)[0]; hi = np.where(np.abs(X[:, ax] - P) < tol)[0]
        o = 1 - ax
        tl = cKDTree(X[lo][:, [o]])
        d, idx = tl.query(X[hi][:, [o]])
        assert np.all(d < 1e-8), 'unmatched periodic nodes'
        for h, l in zip(hi, lo[idx]):
            master[h] = l; shift[h, ax] = P
    # resolve chains (a node on both x=P and y=P)
    for i in range(N):
        while master[master[i]] != master[i]:
            shift[i] += shift[master[i]]; master[i] = master[master[i]]
    return master, shift


def _P(master, shift, k=None, X=None):
    """dof map full(3N) <- reduced (3 * #masters); Bloch phases if k given (kx, ky)."""
    N = len(master)
    ms = np.unique(master)
    red = -np.ones(N, int); red[ms] = np.arange(len(ms))
    rows = np.arange(3 * N)
    cols = (3 * red[master][:, None] + np.arange(3)).ravel()
    if k is None:
        vals = np.ones(3 * N)
    else:
        ph = np.exp(1j * (shift @ np.asarray(k)))
        vals = np.repeat(ph, 3)
    return sp.csr_matrix((vals, (rows, cols)), shape=(3 * N, 3 * len(ms))), ms


def mass_vector(X, seg, lam_s):
    L = np.linalg.norm(X[seg[:, 1]] - X[seg[:, 0]], axis=1)
    m = lam_s * L
    M = np.zeros(3 * len(X))
    for e in (0, 1):
        np.add.at(M, 3 * seg[:, e], 0.5 * m); np.add.at(M, 3 * seg[:, e] + 1, 0.5 * m)
        np.add.at(M, 3 * seg[:, e] + 2, m * L ** 2 / 24)
    return M, L


def build(g, lam, frame_lam=None):
    """frame_lam: thickness of the square frame (generation 0), shared between neighbouring
    tiles; None = no frame (ligament ends meet their mirror images at bare kinks)."""
    if frame_lam is None:
        nodes, seg, sl, bn = srtree.network(g)
        X, S, Ls = mirror_supercell(nodes, seg, lam[sl])
    else:
        nodes, seg, sl, bn = srtree.network(g, frame=True)
        lam_ext = np.r_[lam, np.full(4, frame_lam)]
        X, S, Ls = mirror_supercell(nodes, seg, lam_ext[sl])
        # one copy of each frame segment: drop duplicates and the images on x = 2, y = 2
        key = np.sort(S, axis=1)
        _, first = np.unique(key, axis=0, return_index=True)
        keep = np.zeros(len(S), bool); keep[first] = True
        on2 = lambda ax: (np.abs(X[S[:, 0], ax] - 2) < 1e-9) & (np.abs(X[S[:, 1], ax] - 2) < 1e-9)
        keep &= ~on2(0) & ~on2(1)
        S, Ls = S[keep], Ls[keep]
    K = frame.assemble(X, S, Ls)[0]
    master, shift = pairing(X)
    M, L = mass_vector(X, S, Ls)
    rho = float((Ls * L).sum() / 4.0)
    return dict(X=X, S=S, lam=Ls, K=K, M=M, master=master, shift=shift, rho=rho)


def periodic_C(sc):
    X, K = sc['X'], sc['K']
    P, ms = _P(sc['master'], sc['shift'])
    Kr = (P.T @ K @ P).tocsc()
    n = Kr.shape[0]
    keep = np.ones(n, bool); keep[:2] = False         # pin translations of the first master
    Kk = Kr[keep][:, keep]
    Kk = Kk + 1e-13 * abs(Kk.diagonal()).max() * sp.eye(Kk.shape[0], format='csc')
    lu = spl.splu(Kk.tocsc())
    U = []
    x = X - 1.0
    for e in ((1, 0, 0), (0, 1, 0), (0, 0, 1)):
        exx, eyy, gm = e
        ua = np.zeros(3 * len(X))
        ua[0::3] = exx * x[:, 0] + 0.5 * gm * x[:, 1]
        ua[1::3] = eyy * x[:, 1] + 0.5 * gm * x[:, 0]
        rhs = -(P.T @ (K @ ua))
        w = np.zeros(n); w[keep] = lu.solve(rhs[keep])
        U.append(ua + P @ w)
    U = np.array(U)
    C = U @ (K @ U.T) / 4.0
    return 0.5 * (C + C.T)


def bloch_x(sc, n0, kmax=None, nk=40, nev=12, direction=0.0):
    """Longitudinal branch along direction (radians). Returns arrays k, omega, pL."""
    X, K, M = sc['X'], sc['K'], sc['M']
    d = np.array([np.cos(direction), np.sin(direction)])
    mtr = M[0::3]
    if kmax is None:
        kmax = 2 * np.pi            # wavelength = one original cell
    ks = np.linspace(kmax / nk, kmax, nk)
    out = []
    wguess = ks[0] / n0
    for kk in ks:
        kv = kk * d
        P, ms = _P(sc['master'], sc['shift'], kv)
        Kr = (P.conj().T @ K @ P).tocsc(); Mr = (P.conj().T @ sp.diags(M) @ P).tocsc()
        Kr = 0.5 * (Kr + Kr.conj().T); Mr = 0.5 * (Mr + Mr.conj().T)
        sig = wguess ** 2
        try:
            w2, V = spl.eigsh(Kr, k=nev, M=Mr, sigma=sig, which='LM')
        except Exception:
            out.append((kk, np.nan, np.nan)); continue
        U = (P @ V)                                   # full complex displacements
        ph = np.exp(-1j * (X @ kv))
        ux = U[0::3]; uy = U[1::3]
        ul = d[0] * ux + d[1] * uy
        num = np.abs((mtr * ph) @ ul) ** 2
        den = mtr.sum() * (mtr[:, None] * (np.abs(ux) ** 2 + np.abs(uy) ** 2)).sum(0)
        pL = num / den
        j = int(np.argmax(pL))
        w = np.sqrt(max(w2[j], 0))
        out.append((kk, w, pL[j]))
        if pL[j] > 0.3:
            wguess = w * (kk + ks[0]) / kk
        else:
            wguess = (kk + ks[0]) / n0
    return np.array(out)
