"""Maxwell fish-eye lens built from recursive-division cells (homogenised device simulation).

Cells: framed periodic tilings, alpha = 0 (results/psweep_T1000.json). Each rule's mean stiffness
is isotropised (2D bulk K and shear G) and tabulated against density; the isotropic index
n = sqrt(rho/(K+G)) is monotone in rho, so a target n fixes (rho, K, G) by log interpolation.
Band edge of each cell: Bloch omega_h(n) of its rule (results/bloch2_T1000.json, interpolated at
fixed n as in Fig. 2b); the device band is the lowest omega_h among the cells it uses.

Designs (lens radius R cells, profile n(r) = n0 * 2 / (1 + r^2/R^2) inside, n0 outside):
  mixed     balanced cells for n < n_switch, centroid cells above (n0 = 2.6, centre 5.2)
  centroid  / balanced / SR alone: n0 = lowest index of the rule, centre clipped at its highest.
Solver: 2D elastodynamics, bilinear quads (h = L/3), element-wise isotropic (K, G, rho),
harmonic point force (x) at the rim point (-R, 0); absorbing layer (mass-proportional damping,
quadratic ramp) around the domain. Image metric: mean |u|^2 within one cell of (R, 0) over the
mean |u|^2 on the circle r = R away from the source (|theta - pi| > pi/6).
Units: E0 = rho0 = 1, lengths in cells L, omega in c0/L.
"""
import json, os, sys, time
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spl

B = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..')


def cell_tables():
    p = json.load(open(os.path.join(B, 'results', 'psweep_T1000.json')))
    tab = {}
    for r in ('SR', 'centroid', 'balanced', 'regular'):
        rows = []
        for q in (0.1, 0.2, 0.3, 0.4):
            v = [x for x in p if x['rule'] == r and x['alpha'] == 0 and x['rho_t'] == q]
            C = np.mean([np.array(x['C']) for x in v], axis=0); rho = np.mean([x['rho'] for x in v])
            K = (C[0, 0] + C[1, 1] + 2 * C[0, 1]) / 4
            G = ((C[0, 0] + C[1, 1] - 2 * C[0, 1]) / 4 + C[2, 2]) / 2
            rows.append((np.sqrt(rho / (K + G)), rho, K, G))
        rows = np.array(sorted(rows))
        tab[r] = rows
    return tab


def band_tables():
    bl = json.load(open(os.path.join(B, 'results', 'bloch2_T1000.json')))
    curves = {}
    for x in bl:
        if x['resolved'] and x['omega_h'] > 0:
            curves.setdefault((x['rule'], x['seed'], x['dir']), []).append((x['n0'], x['omega_h']))
    grid = np.arange(2.0, 8.01, 0.1); out = {}
    for r in ('SR', 'centroid', 'balanced'):
        G = []; M = []
        for n in grid:
            v = []
            for (rr, _, _), cv in curves.items():
                if rr != r: continue
                cv = sorted(cv); xx = np.array([q[0] for q in cv]); yy = np.log([q[1] for q in cv])
                if xx.min() <= n <= xx.max(): v.append(np.interp(n, xx, yy))
            if len(v) >= 3: G.append(n); M.append(np.mean(v))
        out[r] = (np.array(G), np.array(M))
    return out


def props(tab, rule, n):
    t = tab[rule]
    ln = np.log(np.clip(n, t[0, 0], t[-1, 0]))
    f = lambda col: np.exp(np.interp(ln, np.log(t[:, 0]), np.log(t[:, col])))
    return f(1), f(2), f(3)


def omega_h(bt, rule, n):
    G, M = bt[rule]
    return float(np.exp(np.interp(n, G, M)))


def design(name, tab, R):
    lo = {r: tab[r][0, 0] for r in tab}; hi = {r: tab[r][-1, 0] for r in tab}
    if name == 'mixed':
        n0 = 2.6; n_sw = 3.5
        rule_of = lambda n: np.where(n < n_sw, 'balanced', 'centroid')
    else:
        n0 = lo[name]; rule_of = lambda n: np.full(np.shape(n), name)
    def nprof(r):
        n = np.where(r < R, n0 * 2 / (1 + (r / R) ** 2), n0)
        if name != 'mixed':
            n = np.clip(n, lo[name], hi[name])
        else:
            n = np.clip(n, lo['balanced'], hi['centroid'])
        return n
    return n0, nprof, rule_of


def solve(name, tab, bt, R, omega, h=1 / 3, margin=4, W=None):
    n0, nprof, rule_of = design(name, tab, R)
    if W is None:
        W = max(8.0, 2 * np.pi / (n0 * omega))            # absorbing layer about one background wavelength
    H = R + margin + W
    N = int(round(2 * H / h)); xs = np.linspace(-H, H, N + 1)
    X, Y = np.meshgrid(xs, xs, indexing='ij')
    # element centres
    xc = 0.5 * (xs[:-1] + xs[1:]); XC, YC = np.meshgrid(xc, xc, indexing='ij')
    rc = np.hypot(XC, YC)
    n_e = nprof(rc); rules = rule_of(n_e)
    rho = np.zeros_like(n_e); K = np.zeros_like(n_e); G = np.zeros_like(n_e)
    for r in np.unique(rules):
        m = rules == r
        rho[m], K[m], G[m] = props(tab, r, n_e[m])
    # Q4 reference matrices (2x2 Gauss), D = K*DK + G*DG with DK = [[1,1,0],[1,1,0],[0,0,0]], DG = [[1,-1,0],[-1,1,0],[0,0,1]]
    gp = np.array([-1, 1]) / np.sqrt(3)
    DK = np.array([[1, 1, 0], [1, 1, 0], [0, 0, 0.]]); DG = np.array([[1, -1, 0], [-1, 1, 0], [0, 0, 1.]])
    KK = np.zeros((8, 8)); KG = np.zeros((8, 8))
    for xi in gp:
        for eta in gp:
            dN = 0.25 * np.array([[-(1 - eta), (1 - eta), (1 + eta), -(1 + eta)],
                                  [-(1 - xi), -(1 + xi), (1 + xi), (1 - xi)]]) * (2 / h)
            Bm = np.zeros((3, 8))
            Bm[0, 0::2] = dN[0]; Bm[1, 1::2] = dN[1]; Bm[2, 0::2] = dN[1]; Bm[2, 1::2] = dN[0]
            w = (h / 2) ** 2
            KK += w * Bm.T @ DK @ Bm; KG += w * Bm.T @ DG @ Bm
    nid = lambda i, j: i * (N + 1) + j
    I, J = np.meshgrid(np.arange(N), np.arange(N), indexing='ij')
    nodes = np.stack([nid(I, J), nid(I + 1, J), nid(I + 1, J + 1), nid(I, J + 1)], -1).reshape(-1, 4)
    dof = np.repeat(2 * nodes, 2, axis=1) + np.tile([0, 1], 4)
    Ke = K.ravel()[:, None, None] * KK + G.ravel()[:, None, None] * KG
    nd = 2 * (N + 1) ** 2
    Kg = sp.coo_matrix((Ke.ravel(), (np.repeat(dof, 8, axis=1).ravel(), np.tile(dof, (1, 8)).ravel())), shape=(nd, nd)).tocsr()
    me = rho.ravel() * h * h / 4
    Mn = np.bincount(nodes.ravel(), weights=np.repeat(me, 4), minlength=(N + 1) ** 2)
    M = np.repeat(Mn, 2)
    d = np.clip((np.maximum(np.abs(X), np.abs(Y)) - (H - W)) / W, 0, 1).ravel()
    eta = np.repeat(3.0 * omega * d ** 2, 2)
    A = (Kg - omega ** 2 * sp.diags(M) + 1j * omega * sp.diags(eta * M)).tocsc()
    f = np.zeros(nd, complex)
    src = np.exp(-((X + R) ** 2 + Y ** 2) / (2 * 0.5 ** 2)).ravel()
    f[0::2] = src / src.sum()
    u = spl.spsolve(A, f)
    I2 = (np.abs(u[0::2]) ** 2 + np.abs(u[1::2]) ** 2).reshape(N + 1, N + 1)
    rn = np.hypot(X, Y); th = np.arctan2(Y, X)
    img = np.hypot(X - R, Y) < 1.0
    ring = (np.abs(rn - R) < 0.5) & (np.abs(np.abs(th) - np.pi) > np.pi / 6)
    gain = float(I2[img].mean() / I2[ring].mean())
    used = np.unique(rules[rc < R + 1e-9]) if np.any(rc < R) else np.unique(rules)
    wband = min(omega_h(bt, r, float(n_e[(rules == r) & (rc <= R)].max())) for r in used)
    # also the background cells
    wband = min(wband, omega_h(bt, str(rule_of(np.array(n0))), n0))
    return dict(gain=gain, omega_band=wband, n0=n0, xs=xs, I=I2, n_map=n_e, rule_map=rules)


if __name__ == '__main__':
    tab = cell_tables(); bt = band_tables()
    R = float(sys.argv[1]) if len(sys.argv) > 1 else 20.0
    out = {}
    for name in ('mixed', 'centroid', 'balanced', 'SR'):
        res = []
        for Rl in (1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 5.0):     # lens radius in background wavelengths
            n0 = design(name, tab, R)[0]
            om = 2 * np.pi * Rl / (n0 * R)
            t0 = time.time()
            s = solve(name, tab, bt, R, om)
            res.append(dict(R_over_lambda0=Rl, omega=om, gain=s['gain'], omega_band=s['omega_band'], n0=s['n0']))
            print(name, 'R/lambda0 %.1f omega %.3f gain %.2f band %.3f %s %.0fs' % (
                Rl, om, s['gain'], s['omega_band'], 'IN' if om < s['omega_band'] else 'out', time.time() - t0), flush=True)
            if Rl in (2.0, 3.0):
                np.savez_compressed(os.path.join(B, 'results', f'fisheye_{name}_{Rl}.npz'), xs=s['xs'], I=s['I'],
                                    n_map=s['n_map'], rule_map=(s['rule_map'] == 'centroid').astype(np.int8))
        out[name] = res
        json.dump(out, open(os.path.join(B, 'results', f'fisheye_R{R:g}.json'), 'w'), indent=1)
