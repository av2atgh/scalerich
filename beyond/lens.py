"""Lens step 2: Luneburg lens tiled from one division geometry per rule (2D frames).

Tiles: the calibrated cell (results/lens_calib.json), mirror-tiled so nodes match across
tile edges; every tile carries its frame (generation 0), shared edges take the mean of the
two tiles' thicknesses. Tile (i, j) has uniform ligament thickness set by the density that
gives n(r) = n_b sqrt(2 - (r/R)^2) inside the lens (r = distance of the tile centre) and n_b
outside, by log-interpolation of the calibration n(rho).
Drive: line force F_x ∝ mass on a strip at x = xs (plane longitudinal wave), harmonic omega.
Absorbing layers: mass-proportional damping eta(d) = eta_max (d/W)^2 in the W outer tiles
on all four sides. Field measured as the mass-weighted mean complex displacement in
half-tile bins (the effective wave, local ligament motion averaged out). Output: |u|^2 per node and the focal gain
G = I(lens)/I(no lens) at the focus (the back edge of the lens on its axis).
Units: E0 = rho0 = 1, lengths in tiles; omega in c0 / tile.
"""
import json, sys, time
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spl
from scipy.spatial import cKDTree
import divisions, periodic
import os
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'rg'))
import srtree, frame  # noqa: E402

NB = {'SR': 3.4, 'balanced': 2.9}


def tile_geometry(rule, TT, seed):
    g = divisions.run(TT, seed, rule)
    nodes, seg, sl, bn = srtree.network(g, frame=True)
    L = np.linalg.norm(g['B'] - g['A'], axis=1)
    return nodes, seg, sl, L, g['T']


def rho_for_n(cal, n_target):
    rho = np.array(cal['rho_int']); n = np.array(cal['n'])
    o = np.argsort(n)
    return float(np.exp(np.interp(n_target, n[o], np.log(rho[o]))))


def build_domain(rule, cal, TT, Nx, Ny, xc, yc, R, lens=True):
    nodes, seg, sl, Llig, T = tile_geometry(rule, TT, cal['seed'])
    nb = NB[rule]
    X, S, LAM, TILE = [], [], [], []
    N0 = len(nodes); k = 0
    rho_map = np.zeros((Nx, Ny))
    for i in range(Nx):
        for j in range(Ny):
            x = nodes.copy()
            if i % 2: x[:, 0] = 1 - x[:, 0]
            if j % 2: x[:, 1] = 1 - x[:, 1]
            x += [i, j]
            r = np.hypot(i + 0.5 - xc, j + 0.5 - yc)
            nt = nb * np.sqrt(2 - (r / R) ** 2) if (lens and r < R) else nb
            rho = rho_for_n(cal, nt); rho_map[i, j] = rho
            lam = rho / Llig.sum()
            X.append(x); S.append(seg + k * N0); LAM.append(np.full(len(seg), lam)); k += 1
    X = np.vstack(X); S = np.vstack(S); LAM = np.concatenate(LAM)
    pairs = cKDTree(X).query_pairs(1e-10, output_type='ndarray')
    parent = np.arange(len(X))

    def find(a):
        while parent[a] != a:
            parent[a] = parent[parent[a]]; a = parent[a]
        return a
    for a, b in pairs:
        ra, rb = find(a), find(b)
        if ra != rb: parent[max(ra, rb)] = min(ra, rb)
    root = np.array([find(a) for a in range(len(X))])
    uniq, gid = np.unique(root, return_inverse=True)
    X = X[uniq]; S = gid[S]
    key = np.sort(S, axis=1)
    uk, inv = np.unique(key, axis=0, return_inverse=True)
    lam_u = np.bincount(inv, weights=LAM) / np.bincount(inv)
    return X, uk, lam_u, rho_map


def solve(X, S, lam, omega, Nx, Ny, xs, W=2.5, eta_max=None):
    """line force F_x ∝ mass on the strip |x - xs| < 0.25 (plane wave both ways);
    absorbing layers (mass-proportional damping, quadratic ramp) on all four sides."""
    K = frame.assemble(X, S, lam)[0]
    M, L = periodic.mass_vector(X, S, lam)
    d = np.zeros(len(X))
    d = np.maximum(d, (W - X[:, 0]) / W)
    d = np.maximum(d, (X[:, 0] - (Nx - W)) / W)
    d = np.maximum(d, (W - X[:, 1]) / W)
    d = np.maximum(d, (X[:, 1] - (Ny - W)) / W)
    d = np.clip(d, 0, 1)
    if eta_max is None:
        eta_max = 2.0 * omega
    eta = np.repeat(eta_max * d ** 2, 3)
    A = (K - omega ** 2 * sp.diags(M) + 1j * omega * sp.diags(eta * M)).tocsc()
    f = np.zeros(3 * len(X), complex)
    src = np.abs(X[:, 0] - xs) < 0.25
    f[3 * np.where(src)[0]] = M[3 * np.where(src)[0]]
    return spl.spsolve(A, f)


def intensity(X, u):
    return np.abs(u[0::3]) ** 2 + np.abs(u[1::3]) ** 2


def coarse(X, u, m, Nx, Ny, b=0.5):
    """mass-weighted mean complex displacement in b x b bins -> effective-wave intensity map"""
    ix = np.minimum((X[:, 0] / b).astype(int), int(Nx / b) - 1)
    iy = np.minimum((X[:, 1] / b).astype(int), int(Ny / b) - 1)
    nb = (int(Nx / b), int(Ny / b))
    idx = ix * nb[1] + iy
    W = np.bincount(idx, weights=m, minlength=nb[0] * nb[1])
    ux = np.bincount(idx, weights=m * u[0::3].real, minlength=W.size) + 1j * np.bincount(idx, weights=m * u[0::3].imag, minlength=W.size)
    uy = np.bincount(idx, weights=m * u[1::3].real, minlength=W.size) + 1j * np.bincount(idx, weights=m * u[1::3].imag, minlength=W.size)
    W = np.where(W > 0, W, 1)
    I = (np.abs(ux / W) ** 2 + np.abs(uy / W) ** 2).reshape(nb)
    return I


if __name__ == '__main__':
    TT = 200
    cal = json.load(open('results/lens_calib.json'))
    R = float(sys.argv[1]) if len(sys.argv) > 1 else 4.0
    Nx, Ny = int(2 * R + 9), int(2 * R + 7)
    xs = 3.0                       # source line, just inside the left absorbing layer
    xc, yc = xs + 1.0 + R, Ny / 2
    out = {}
    for rule in ('SR', 'balanced'):
        nb = NB[rule]
        Xl, Sl, ll, rmap = build_domain(rule, cal[rule], TT, Nx, Ny, xc, yc, R, True)
        Xr, Sr, lr, _ = build_domain(rule, cal[rule], TT, Nx, Ny, xc, yc, R, False)
        print(rule, 'nodes', len(Xl), 'rho range %.3f-%.3f' % (rmap.min(), rmap.max()), flush=True)
        res = []
        for Dl in (0.75, 1.0, 1.25, 1.5, 1.75, 2.0, 2.25, 2.5, 2.75, 3.0, 3.5, 4.0):
            lamb = 2 * R / Dl     # background wavelength in tiles
            om = 2 * np.pi / (nb * lamb)
            t0 = time.time()
            ul = solve(Xl, Sl, ll, om, Nx, Ny, xs); ur = solve(Xr, Sr, lr, om, Nx, Ny, xs)
            ml = periodic.mass_vector(Xl, Sl, ll)[0][0::3]; mr = periodic.mass_vector(Xr, Sr, lr)[0][0::3]
            Il = coarse(Xl, ul, ml, Nx, Ny); Ir = coarse(Xr, ur, mr, Nx, Ny)
            bx = int((xc + R) / 0.5)                      # focal column (bins of 0.5 tile)
            by0, by1 = int((yc - 0.5) / 0.5), int((yc + 0.5) / 0.5)
            col = slice(int(3.0 / 0.5), int((Ny - 3.0) / 0.5))
            I_f = float(Il[bx, by0:by1].mean())
            G = I_f / float(Ir[bx, col].mean())
            contrast = I_f / float(Il[bx, col].mean())
            ys = (np.arange(Il.shape[1]) + 0.5) * 0.5
            prof = Il[bx].tolist()
            res.append(dict(lambda_b=lamb, omega=om, D_over_lambda=2 * R / lamb, gain=G, contrast=contrast,
                            ys=ys.tolist(), profile=prof, profile_ref=Ir[bx].tolist()))
            np.savez_compressed(f'results/lens_R{R:g}_{rule}_{Dl}.npz', Il=Il, Ir=Ir)
            print(rule, 'lambda_b %.2f omega %.3f D/lambda %.2f gain %.2f contrast %.2f %.0fs' % (lamb, om, 2 * R / lamb, G, contrast, time.time() - t0), flush=True)
        out[rule] = dict(nb=nb, R=R, Nx=Nx, Ny=Ny, xc=xc, yc=yc, rho_map=rmap.tolist(), runs=res)
        json.dump(out, open(f'results/lens_R{R:g}.json', 'w'))
