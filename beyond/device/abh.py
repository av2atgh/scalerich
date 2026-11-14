"""Acoustic-black-hole (ABH) vibration trap in a 2.5D division lattice (homogenised device simulation).

A strip of width W cells carries in-plane waves from a uniform host section into a trap of length
L_t that ends at a free edge. In the trap the cells are graded along a path in (kappa, rho):
kappa = cut-angle tilt (0: perpendicular to the long axis, 'balanced'; 1: uniform angle,
'centroid'), rho = relative density; the path parameter follows t(x) = (x/L_t)^2 (ABH-like taper).
Cell properties (isotropised 2D bulk K and shear G, density) come from results/anglegrade.json
(framed periodic tilings, 4 realisations), interpolated linearly in kappa and in ln rho.
Every cell has a polymer-like loss factor eta (complex moduli K(1 + i eta), G(1 + i eta)).

Designs (all start from the host cell kappa = 0, rho = 0.2):
  uniform   no grading (lossy strip of host cells)                       -- reference
  step      abrupt switch to kappa = 1, rho = 0.1 over the whole trap
  density   kappa = 0, rho 0.2 -> 0.1
  angle     rho = 0.2, kappa 0 -> 1                (constant density)
  both      kappa 0 -> 1 and rho 0.2 -> 0.1
Excitation: harmonic line force (x) across the strip in the host section; absorbing layer at the
far left. Reflection: |R|^2 = <|u - u_inc|^2> / <|u_inc|^2> over the host section between source
and trap, u_inc from the same strip with the trap replaced by an absorbing continuation.
Absorbed fraction A = 1 - |R|^2. Validity: below the lowest Bloch band edge among the cells used
(omega_h estimated as 0.80 omega1, the median ratio at rho = 0.2 in anglegrade).
Units: E0 = rho0 = 1, lengths in cells L, omega in c0/L.
"""
import json, os, sys, time
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spl

B = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..')
KAPPAS = np.array([0.0, 0.2, 0.4, 0.6, 0.8, 1.0]); RHOS = np.array([0.1, 0.2, 0.4])


def tables():
    r = json.load(open(os.path.join(B, 'results', 'anglegrade.json')))
    T = {}
    for q in ('K', 'G', 'rho', 'n_iso', 'omega1'):
        T[q] = np.array([[np.exp(np.mean(np.log([o[q] for o in r if o['kappa'] == k and o['rho_t'] == p])))
                          for p in RHOS] for k in KAPPAS])
    v = [o for o in r if o['rho_t'] == 0.2 and o.get('omega_h_x')]
    ratio = np.median([0.5 * (o['omega_h_x'] + o['omega_h_y']) / o['omega1'] for o in v])
    return T, float(ratio)


def cell(T, kappa, rho):
    """bilinear in (kappa, ln rho) of ln K, ln G, ln rho_actual, ln omega1"""
    kappa = np.clip(kappa, 0, 1); lr = np.log(np.clip(rho, 0.1, 0.4))
    out = []
    for q in ('K', 'G', 'rho', 'omega1'):
        A = np.log(T[q])
        a = np.array([np.interp(lr, np.log(RHOS), A[i]) for i in range(len(KAPPAS))])
        out.append(np.exp(np.array([np.interp(k, KAPPAS, a[:, j]) for j, k in enumerate(np.atleast_1d(kappa))])
                          if np.ndim(kappa) else np.interp(kappa, KAPPAS, a)))
    return out


def path(design, t):
    t = np.asarray(t, float)
    if design == 'uniform': return 0 * t, 0.2 + 0 * t
    if design == 'step':    return np.ones_like(t), 0.1 + 0 * t
    if design == 'density': return 0 * t, 0.2 * (0.5 ** t)
    if design == 'angle':   return t, 0.2 + 0 * t
    if design == 'both':    return t, 0.2 * (0.5 ** t)
    raise ValueError(design)


def q4_refs(h):
    gp = np.array([-1, 1]) / np.sqrt(3)
    DK = np.array([[1, 1, 0], [1, 1, 0], [0, 0, 0.]]); DG = np.array([[1, -1, 0], [-1, 1, 0], [0, 0, 1.]])
    KK = np.zeros((8, 8)); KG = np.zeros((8, 8))
    for xi in gp:
        for eta in gp:
            dN = 0.25 * np.array([[-(1 - eta), (1 - eta), (1 + eta), -(1 + eta)],
                                  [-(1 - xi), -(1 + xi), (1 + xi), (1 - xi)]]) * (2 / h)
            Bm = np.zeros((3, 8))
            Bm[0, 0::2] = dN[0]; Bm[1, 1::2] = dN[1]; Bm[2, 0::2] = dN[1]; Bm[2, 1::2] = dN[0]
            KK += (h / 2) ** 2 * Bm.T @ DK @ Bm; KG += (h / 2) ** 2 * Bm.T @ DG @ Bm
    return KK, KG


def solve(T, design, omega, loss=0.03, W=4.0, L_abs=None, L_host=24.0, L_t=30.0, h=0.25, incident=False):
    if L_abs is None:
        L_abs = max(15.0, 2.5 * 2 * np.pi / (2.9 * omega))
    Lx = L_abs + L_host + L_t
    if incident:
        Lx = L_abs + L_host + L_t + L_abs                 # trap replaced by host + absorber
    Nx = int(round(Lx / h)); Ny = int(round(W / h))
    xc = (np.arange(Nx) + 0.5) * h
    x0 = L_abs + L_host                                   # trap start
    t = np.clip((xc - x0) / L_t, 0, 1) ** 2
    if incident:
        kap, rho_t = path('uniform', t)
    else:
        kap, rho_t = path(design, t)
        kap = np.where(xc < x0, 0.0, kap); rho_t = np.where(xc < x0, 0.2, rho_t)
    K, G, rho, w1 = cell(T, kap, rho_t)
    KK, KG = q4_refs(h)
    nid = lambda i, j: i * (Ny + 1) + j
    I, J = np.meshgrid(np.arange(Nx), np.arange(Ny), indexing='ij')
    nodes = np.stack([nid(I, J), nid(I + 1, J), nid(I + 1, J + 1), nid(I, J + 1)], -1).reshape(-1, 4)
    dof = np.repeat(2 * nodes, 2, axis=1) + np.tile([0, 1], 4)
    Ke_ = np.repeat(K, Ny)[:, None, None] * KK + np.repeat(G, Ny)[:, None, None] * KG
    nd = 2 * (Nx + 1) * (Ny + 1)
    Kg = sp.coo_matrix(((1 + 1j * loss) * Ke_.ravel(), (np.repeat(dof, 8, axis=1).ravel(), np.tile(dof, (1, 8)).ravel())),
                       shape=(nd, nd)).tocsr()
    me = np.repeat(rho, Ny) * h * h / 4
    M = np.repeat(np.bincount(nodes.ravel(), weights=np.repeat(me, 4), minlength=(Nx + 1) * (Ny + 1)), 2)
    xn = np.repeat(np.arange(Nx + 1) * h, Ny + 1)
    d = np.clip((L_abs - xn) / L_abs, 0, 1)
    if incident:
        d = np.maximum(d, np.clip((xn - (Lx - L_abs)) / L_abs, 0, 1))
    eta = np.repeat(3.0 * omega * d ** 2, 2)
    A = (Kg - omega ** 2 * sp.diags(M) + 1j * omega * sp.diags(eta * M)).tocsc()
    xs = L_abs + 2.0
    f = np.zeros(nd, complex)
    src = np.abs(xn - xs) < h / 2
    f[0::2][src] = 1.0 / src.sum()
    u = spl.spsolve(A, f)
    ux = u[0::2].reshape(Nx + 1, Ny + 1).mean(axis=1)    # width-averaged axial displacement
    xg = np.arange(Nx + 1) * h
    used_w1 = float(w1[xc >= x0].min()) if np.any(xc >= x0) else float(w1.min())
    mass_trap = float((rho[xc >= x0]).sum() * h * W)
    return xg, ux, used_w1, mass_trap


def reflection(T, design, omega, **kw):
    xg, u, w1, mass = solve(T, design, omega, **kw)
    xg0, u0, _, _ = solve(T, design, omega, incident=True, **kw)
    L_abs = kw.get('L_abs') or max(15.0, 2.5 * 2 * np.pi / (2.9 * omega))
    L_host = kw.get('L_host', 24.0)
    sel = (xg > L_abs + 6.0) & (xg < L_abs + L_host - 2.0)
    n = min(len(u), len(u0))
    R2 = float(np.mean(np.abs(u[sel] - u0[:len(u)][sel]) ** 2) / np.mean(np.abs(u0[:len(u)][sel]) ** 2))
    return R2, w1, mass


if __name__ == '__main__':
    T, ratio = tables()
    print('omega_h/omega1 at rho 0.2 (median) %.2f' % ratio)
    for k in (0.0, 0.5, 1.0):
        K, G, rho, w1 = cell(T, np.array([k]), np.array([0.2]))
        print('kappa %.1f rho 0.2: n %.2f' % (k, float(np.sqrt(rho / (K + G))[0])))
    loss = float(sys.argv[1]) if len(sys.argv) > 1 else 0.03
    out = {'ratio': ratio, 'loss': loss, 'runs': {}}
    omegas = np.round(np.geomspace(0.03, 0.8, 16), 4)
    for design in ('uniform', 'step', 'density', 'angle', 'both'):
        res = []
        for om in omegas:
            t0 = time.time()
            R2, w1, mass = reflection(T, design, float(om), loss=loss)
            res.append(dict(omega=float(om), R2=R2, A=1 - R2, band=ratio * w1, mass_trap=mass))
            print(design, 'omega %.3f  A %.3f  band %.3f %s  mass %.2f  %.1fs' % (
                om, 1 - R2, ratio * w1, 'IN' if om < ratio * w1 else 'out', mass, time.time() - t0), flush=True)
        out['runs'][design] = res
        json.dump(out, open(os.path.join(B, 'results', f'abh_loss{loss:g}.json'), 'w'), indent=1)
