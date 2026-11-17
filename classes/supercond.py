"""Superconducting wire networks on the divisions: mean-field Tc(B) from the linearised
Ginzburg-Landau equation on the network (de Gennes-Alexander), solved with linear finite
elements along each wire (element length <= h) and Peierls phases, Landau gauge A = (-B y, 0),
flux quantum Phi0 = 1. Lowest generalised eigenvalue k^2 of K psi = k^2 M psi gives
(Tc0 - Tc(B))/Tc0 = xi0^2 k^2. Field as frustration f = B <A_cell> / Phi0, <A_cell> = 1/(T+1).
Open sample, frame included; X-crossings are single junctions. -> results/supercond.json"""
import sys, json, time, numpy as np
import scipy.sparse as sp, scipy.sparse.linalg as spl
from scipy.spatial import cKDTree
sys.path[:0] = ['../rg', '../beyond']
import qdiv, divisions, srtree, periodic

def wires(g):
    """junction network: unique nodes, segments (frame included), coincident nodes merged"""
    nodes, seg, sl, bn = srtree.network(g, frame=True)
    pairs = cKDTree(nodes).query_pairs(1e-10, output_type='ndarray')
    par = np.arange(len(nodes))
    for i, j in pairs:
        a, b = par[i], par[j]
        while par[a] != a: a = par[a]
        while par[b] != b: b = par[b]
        if a != b: par[max(a, b)] = min(a, b)
    for i in range(len(par)):
        r = i
        while par[r] != r: r = par[r]
        par[i] = r
    u, gid = np.unique(par, return_inverse=True)
    S = gid[seg]; S = S[S[:, 0] != S[:, 1]]
    S = np.unique(np.sort(S, axis=1), axis=0)
    return nodes[u], S

def discretise(X, S, h):
    """split every segment into elements of length <= h; returns node coords, element pairs"""
    P = [X]; E = []; n = len(X)
    for a, b in S:
        L = np.linalg.norm(X[b] - X[a]); k = max(1, int(np.ceil(L / h)))
        idx = [a]
        if k > 1:
            t = np.arange(1, k)[:, None] / k
            P.append(X[a] + t * (X[b] - X[a])); idx += list(range(n, n + k - 1)); n += k - 1
        idx.append(b)
        E += list(zip(idx[:-1], idx[1:]))
    return np.vstack(P), np.array(E)

def matrices(Y, E, B):
    a, b = E[:, 0], E[:, 1]; d = Y[b] - Y[a]; h = np.linalg.norm(d, axis=1)
    ymid = 0.5 * (Y[a, 1] + Y[b, 1])
    gam = 2 * np.pi * (-B * ymid * d[:, 0])          # phase from a to b
    e = np.exp(1j * gam); N = len(Y)
    rows = np.r_[a, b, a, b]; cols = np.r_[a, b, b, a]
    K = sp.csc_matrix((np.r_[1 / h, 1 / h, -np.conj(e) / h, -e / h], (rows, cols)), shape=(N, N))
    M = sp.csc_matrix((np.r_[h / 3, h / 3, h * np.conj(e) / 6, h * e / 6], (rows, cols)), shape=(N, N))
    return K, M

def k2min(Y, E, B):
    K, M = matrices(Y, E, B)
    w = spl.eigsh(K, k=1, M=M, sigma=-1e-3, which='LM', return_eigenvectors=False)
    return float(np.real(w[0]))


def torus(g, h):
    """mirror-periodic 2x2 supercell (frame shared, duplicates on x=2, y=2 removed), discretised;
    returns coords, elements, and for each element the boundary wrap needed for magnetic
    periodicity: node on y=2 -> master on y=0 with psi_top = exp(-2 pi i B Ly x) psi_bottom."""
    nodes, seg, sl, bn = srtree.network(g, frame=True)
    X, S, _ = periodic.mirror_supercell(nodes, seg, np.ones(len(sl)))
    S = S[S[:, 0] != S[:, 1]]; S = np.unique(np.sort(S, axis=1), axis=0)
    on2 = lambda ax: (np.abs(X[S[:, 0], ax] - 2) < 1e-9) & (np.abs(X[S[:, 1], ax] - 2) < 1e-9)
    S = S[~on2(0) & ~on2(1)]
    Y, E = discretise(X, S, h)
    # identify nodes on x=2 with x=0 and on y=2 with y=0
    key = lambda P: np.round(P, 9)
    idx = {tuple(k): i for i, k in enumerate(key(Y))}
    wrapy = np.zeros(len(Y))                       # x coordinate for phase if node was on y=2
    master = np.arange(len(Y)); flag = np.zeros(len(Y), bool)
    for i, (x, y) in enumerate(Y):
        xm, ym = x, y; wy = False
        if abs(x - 2) < 1e-9: xm = 0.0
        if abs(y - 2) < 1e-9: ym = 0.0; wy = True
        if (xm, ym) != (x, y):
            master[i] = idx[tuple(key(np.array([xm, ym])))]
            if wy: flag[i] = True; wrapy[i] = xm
    used = np.unique(master); ren = -np.ones(len(Y), int); ren[used] = np.arange(len(used))
    return Y, E, master, flag, wrapy, ren, len(used)

def tmatrices(Y, E, master, flag, wrapy, ren, n, B, Ly=2.0):
    a, b = E[:, 0], E[:, 1]; d = Y[b] - Y[a]; h = np.linalg.norm(d, axis=1)
    ymid = 0.5 * (Y[a, 1] + Y[b, 1])
    gam = 2 * np.pi * (-B * ymid * d[:, 0])
    chi = lambda i: np.where(flag[i], -2 * np.pi * B * Ly * wrapy[i], 0.0)
    gam = gam - chi(b) + chi(a)
    A, Bn = ren[master[a]], ren[master[b]]
    e = np.exp(1j * gam)
    rows = np.r_[A, Bn, A, Bn]; cols = np.r_[A, Bn, Bn, A]
    K = sp.csc_matrix((np.r_[1 / h, 1 / h, -np.conj(e) / h, -e / h], (rows, cols)), shape=(n, n))
    M = sp.csc_matrix((np.r_[h / 3, h / 3, h * np.conj(e) / 6, h * e / 6], (rows, cols)), shape=(n, n))
    return K, M

def framemask(TT):
    """nodes off the sample frame (x, y not in {0, 1}) of the mirror supercell"""
    Y, E, master, flag, wrapy, ren, n = TT
    pos = np.zeros((n, 2)); pos[ren[master]] = Y
    on = np.zeros(n, bool)
    for ax in (0, 1):
        for c in (0.0, 1.0): on |= np.abs(pos[:, ax] - c) < 1e-9
    return ~on

def k2torus(TT, B, mask=None):
    K, M = tmatrices(*TT, B)
    if mask is not None: K = K[mask][:, mask]; M = M[mask][:, mask]   # psi = 0 on the frame
    # k=6, ncv=40: commensurate fields give degenerate ground levels on which k=1 ARPACK stalls
    lu = spl.splu((K + 1e-3 * M).tocsc())
    op = spl.LinearOperator(K.shape, matvec=lu.solve, dtype=complex)
    for k, ncv in ((6, 40), (24, 100), (64, 220)):   # widen the block if the level is more degenerate
        try:
            w, v = spl.eigsh(K, k=k, M=M, sigma=-1e-3, OPinv=op, ncv=ncv, maxiter=3000, tol=1e-9)
            break
        except spl.ArpackNoConvergence as err:
            if len(err.eigenvalues): w, v = err.eigenvalues, err.eigenvectors
            else: continue
            break
    i = np.argmin(w.real); w = w[i:i + 1]
    p = np.abs(v[:, i]) ** 2; p /= p.sum()
    return float(np.real(w[0])), float(1 / (p ** 2).sum() / len(p))

TAU = qdiv.TAU
CASES = {'periodic': lambda T: qdiv.run(T, ('largest', 0.5, 'perp')),
         'golden-aligned': lambda T: qdiv.run(T, ('largest', TAU ** -2, 'perp')),
         'golden-anti': lambda T: qdiv.run(T, ('largest', ('avoid', TAU ** -2, 'low'), 'perp')),
         'r0.3-aligned': lambda T: qdiv.run(T, ('largest', 0.3, 'perp')),
         'uniform': lambda T: qdiv.run(T, ('largest', 'uniform', 'perp'), seed=0),
         'SR': lambda T: divisions.run(T, 0, 'SR')}

if __name__ == '__main__':
    T = int(sys.argv[1]) if len(sys.argv) > 1 else 255
    df = float(sys.argv[2]) if len(sys.argv) > 2 else 0.01
    names = sys.argv[3].split(',') if len(sys.argv) > 3 else list(CASES)
    hfac = float(sys.argv[4]) if len(sys.argv) > 4 else 0.25
    fout = sys.argv[5] if len(sys.argv) > 5 else 'results/supercond_torus.json'
    try: out = json.load(open(fout))
    except Exception: out = {}
    # torus of area 4: total flux 4B must be an integer -> f = n / (4 (T+1))
    nmax = int(round(3.0 * 4 * (T + 1))); step = max(1, int(round(df * 4 * (T + 1))))
    ns = np.arange(0, nmax + 1, step); fs = ns / (4 * (T + 1))
    for nm in names:
        t0 = time.time(); g = CASES[nm](T)
        TT = torus(g, hfac / np.sqrt(T + 1))
        bulk = len(sys.argv) > 6 and sys.argv[6] == 'bulk'
        mask = framemask(TT) if bulk else None
        res = [k2torus(TT, n / 4.0, mask) for n in ns]
        out['%s_T%d' % (nm, T)] = dict(rule=nm, T=T, f=fs.tolist(), k2=[r[0] for r in res],
                                       pr=[r[1] for r in res], Acell=1.0 / (T + 1), nodes=TT[-1], frame='removed' if bulk else 'kept')
        print(nm, T, 'nodes', TT[-1], 'k2(0)=%.3g' % res[0][0], 'time %.0fs' % (time.time() - t0), flush=True)
        json.dump(out, open(fout, 'w'))
