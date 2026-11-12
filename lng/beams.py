"""Linear-elastic Euler-Bernoulli frame on an SR realisation; per-ligament energy law.

Frame: every ligament is a straight beam (E = 1, unit depth, A = lam, I = lam^3/12)
cut into elements at its junction points; T-junctions are rigid (the ending
ligament shares the node, including the rotation, of the point on the centreline
of the ligament it hits).  Boundary nodes (ligament ends on the square) get the
affine displacement of a uniaxial strain u = (eps x, 0); rotations free.
Per-ligament energy U_i = sum of its element energies.

Fits (log-linear, orientation absorbed by 12 orientation-bin dummies):
  U ~ k^a lam^b f(theta)               (paper: a = 3/2, b = 1/2)
  U ~ l^a lam^b f(theta)
  U ~ k^a l^b lam^c (birth index)^d f(theta)
  U vs the affine stretch energy  U_aff = 1/2 lam l eps^2 cos^4(theta)
Output: results/beams_fit.txt, results/beams_<tag>.csv (per-ligament data)
"""
import sys, os, math
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import srcells, graphstats as gs
from exp_topology import edges_from_hits

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'results')


def build_frame(o):
    T = o['T']
    ends = o['ends']; hit = o['hit']; lam = o['lam'].copy()
    # direction and start point of each centreline (start = minus end)
    P0 = ends[:, 0:2]; P1 = ends[:, 2:4]
    L = np.linalg.norm(P1 - P0, axis=1)
    D = (P1 - P0) / L[:, None]
    nodes = []          # coordinates
    pts = [[] for _ in range(T)]   # (s, node) along each ligament
    end_node = np.zeros((T, 2), int)
    bnd = []
    for t in range(T):
        for col, P in ((1, P0[t]), (0, P1[t])):   # hit[:,1] is the minus end, hit[:,0] the plus end
            h = hit[t, col]
            if h >= 0:
                s = float(np.dot(P - P0[h], D[h]))
                s = min(max(s, 0.0), L[h])
                X = P0[h] + s * D[h]
                nid = len(nodes); nodes.append(X)
                pts[h].append((s, nid))
            else:
                nid = len(nodes); nodes.append(P.copy()); bnd.append(nid)
            end_node[t, 0 if col == 1 else 1] = nid
    elems = []   # (lig, n1, n2)
    for i in range(T):
        lst = [(0.0, end_node[i, 0])] + pts[i] + [(L[i], end_node[i, 1])]
        lst.sort()
        for (s1, a), (s2, b) in zip(lst[:-1], lst[1:]):
            if a != b:
                elems.append((i, a, b))
    return np.array(nodes), np.array(elems), np.array(bnd), L, D, lam


def element_K(x1, x2, A, I):
    dx, dy = x2 - x1
    l = math.hypot(dx, dy)
    if l < 1e-12:
        l = 1e-12
    c, s = dx / l, dy / l
    # a stub shorter than its thickness is not a slender beam: cap its stiffness
    # at that of a square block (l_eff = max(l, thickness)); avoids 1/l^3 blow-up
    l = max(l, math.sqrt(12.0 * I / max(A, 1e-300)))
    EA, EI = A / l, I
    k = np.array([
        [EA, 0, 0, -EA, 0, 0],
        [0, 12 * EI / l**3, 6 * EI / l**2, 0, -12 * EI / l**3, 6 * EI / l**2],
        [0, 6 * EI / l**2, 4 * EI / l, 0, -6 * EI / l**2, 2 * EI / l],
        [-EA, 0, 0, EA, 0, 0],
        [0, -12 * EI / l**3, -6 * EI / l**2, 0, 12 * EI / l**3, -6 * EI / l**2],
        [0, 6 * EI / l**2, 2 * EI / l, 0, -6 * EI / l**2, 4 * EI / l]])
    R = np.zeros((6, 6))
    r = np.array([[c, s, 0], [-s, c, 0], [0, 0, 1]])
    R[:3, :3] = r; R[3:, 3:] = r
    return R.T @ k @ R


def solve(o, eps=1e-3, direction=0.0):
    X, E, bnd, L, D, lam = build_frame(o)
    N = len(X)
    rows, cols, vals, Ke_list = [], [], [], []
    for (i, a, b) in E:
        A = lam[i]; I = lam[i] ** 3 / 12.0
        Ke = element_K(X[a], X[b], A, I)
        Ke_list.append(Ke)
        dofs = [3 * a, 3 * a + 1, 3 * a + 2, 3 * b, 3 * b + 1, 3 * b + 2]
        for p in range(6):
            for q in range(6):
                rows.append(dofs[p]); cols.append(dofs[q]); vals.append(Ke[p, q])
    K = sp.coo_matrix((vals, (rows, cols)), shape=(3 * N, 3 * N)).tocsr()
    u = np.zeros(3 * N)
    # uniaxial strain along `direction`: u = eps (n.x) n, origin at centre
    n = np.array([math.cos(direction), math.sin(direction)])
    fixed = []
    for v in bnd:
        disp = eps * np.dot(X[v] - 0.5, n) * n
        u[3 * v] = disp[0]; u[3 * v + 1] = disp[1]
        fixed += [3 * v, 3 * v + 1]
    fixed = np.array(sorted(set(fixed)))
    free = np.setdiff1d(np.arange(3 * N), fixed)
    Kff = K[free][:, free].tocsc()
    rhs = -K[free][:, fixed] @ u[fixed]
    # tiny regularisation for isolated rotational modes
    Kff = Kff + sp.identity(len(free), format='csc') * 1e-14 * abs(Kff.diagonal()).max()
    u[free] = spla.spsolve(Kff, rhs)
    T = o['T']
    U = np.zeros(T)
    for (i, a, b), Ke in zip(E, Ke_list):
        dofs = [3 * a, 3 * a + 1, 3 * a + 2, 3 * b, 3 * b + 1, 3 * b + 2]
        ue = u[dofs]
        U[i] += 0.5 * ue @ Ke @ ue
    return U, L, D


def fit(y, Xcols, names, theta, nb=12):
    """OLS of log y on log X columns + orientation-bin dummies; returns coefs, R2."""
    bins = np.minimum((theta % math.pi) / math.pi * nb, nb - 1).astype(int)
    Dm = np.zeros((len(y), nb)); Dm[np.arange(len(y)), bins] = 1
    M = np.column_stack(Xcols + [Dm])
    c, *_ = np.linalg.lstsq(M, y, rcond=None)
    res = y - M @ c
    R2 = 1 - res.var() / y.var()
    return dict(zip(names, c[:len(names)])), R2


def analyse(o, U, k, L, tag):
    """age law, k- vs l- laws, within-cohort partials."""
    Td = o['T']; th = o['theta']; lam = o['lam']; idx = np.arange(1, Td + 1)
    m = (o['hit'] >= 0).all(1) & (U > 0) & (k >= 2)
    y = np.log(U[m]); lk, ll, li = np.log(k[m]), np.log(L[m]), np.log(idx[m])
    lla = np.log(lam[m]); one = np.ones(m.sum())
    fmt = lambda r: ' '.join(f"{a}={b:.3f}" for a, b in r.items())
    out = [f"{tag}: rho={1-o['freeA'][-1]:.3f} n={m.sum()}  corr(lk,ll)={np.corrcoef(lk,ll)[0,1]:.3f} corr(lk,li)={np.corrcoef(lk,li)[0,1]:.3f}"]
    for names, cols in [(['i'], [li]), (['k'], [lk]), (['l'], [ll]), (['k', 'lam'], [lk, lla]),
                        (['l', 'lam'], [ll, lla]), (['k', 'i'], [lk, li]), (['l', 'i'], [ll, li]),
                        (['k', 'l', 'i'], [lk, ll, li])]:
        r, R2 = fit(y, [one] + cols, ['c'] + names, th[m])
        r.pop('c')
        out.append(f"   U ~ {'*'.join(names):8s}: {fmt(r)}  R2={R2:.3f}")
    rp, Rp = fit(y, [one, 1.5 * lk + 0.5 * lla], ['c', 'slope_k1.5_lam0.5'], th[m]); rp.pop('c')
    rq, Rq = fit(y, [one, 1.5 * ll + 0.5 * lla], ['c', 'slope_l1.5_lam0.5'], th[m]); rq.pop('c')
    out.append(f"   paper comb: {fmt(rp)} R2={Rp:.3f} | same with l for k: {fmt(rq)} R2={Rq:.3f}")
    # within-cohort: residualise on age (log i, 10 cohort dummies) then partial r of k and l
    nb = 10
    cb = np.minimum((li - li.min()) / (li.max() - li.min() + 1e-12) * nb, nb - 1).astype(int)
    Dc = np.zeros((len(y), nb)); Dc[np.arange(len(y)), cb] = 1
    def resid(v):
        M = np.column_stack([Dc, li])
        c, *_ = np.linalg.lstsq(M, v, rcond=None); return v - M @ c
    ry, rk, rl = resid(y), resid(lk), resid(ll)
    out.append(f"   within-cohort partial corr: U-k {np.corrcoef(ry,rk)[0,1]:.3f}  U-l {np.corrcoef(ry,rl)[0,1]:.3f}  k-l {np.corrcoef(rk,rl)[0,1]:.3f}")
    # U_aff
    Uaff = 0.5 * lam * L * 1e-6 * np.cos(th) ** 4
    out.append(f"   Utot={U.sum():.3e} Voigt(affine stretch, 0.5*E*eps^2*sum lam l cos^4)={Uaff.sum():.3e}")
    return '\n'.join(out)


if __name__ == '__main__':
    cases = [(1000, 0.0, 0.03), (3000, 0.25, 0.02), (3000, 0.5, 0.06), (3000, 0.75, 0.08),
             (3000, 1.0, 0.1), (5000, 0.5, 0.06), (5000, 1.0, 0.1)]
    fo = open(os.path.join(OUT, 'beams_fit.txt'), 'w')
    for (T, alpha, lam0) in cases:
        for seed in (1, 2):
            try:
                o = srcells.run(T, alpha=alpha, lam0=lam0, seed=seed)
                Td = o['T']
                U, L, D = solve(o, eps=1e-3, direction=0.0)
                A = gs.adjacency(edges_from_hits(o['hit']), Td)
                k = gs.degrees(A)
                tag = f"T{T}_a{alpha}_l{lam0}_s{seed}"
                txt = analyse(o, U, k, L, tag)
                print(txt, flush=True); fo.write(txt + '\n'); fo.flush()
                if seed == 1:
                    np.savetxt(os.path.join(OUT, f'beams_{tag}.csv'),
                               np.column_stack([np.arange(1, Td + 1), k, L, o['lam'], o['theta'], U]), delimiter=',',
                               header='i,k,l,lam,theta,U', fmt='%.6g', comments='')
            except Exception as ex:
                fo.write(f"# FAILED T{T} a{alpha} l{lam0} s{seed}: {ex!r}\n"); fo.flush()
                print('FAILED', ex, flush=True)
    fo.close()
