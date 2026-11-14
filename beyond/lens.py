"""Lens step 2: Luneburg lens tiled from one division geometry per rule (2D frames).

Tiles: the calibrated cell (results/lens_calib.json), mirror-tiled so nodes match across
tile edges; every tile carries its frame (generation 0), shared edges take the mean of the
two tiles' thicknesses. Tile (i, j) has uniform ligament thickness set by the density that
gives n(r) = n_b sqrt(2 - (r/R)^2) inside the lens (r = distance of the tile centre) and n_b
outside, by log-interpolation of the calibration n(rho).
Drive: u_x = 1 on the left boundary nodes (plane longitudinal wave), harmonic omega.
Absorbing layers: mass-proportional damping eta(d) = eta_max (d/W)^2 in the W outer tiles
on the right, top and bottom. Output: |u|^2 per node and the focal gain
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


def solve(X, S, lam, omega, Nx, Ny, W=2.5, eta_max=None):
    K = frame.assemble(X, S, lam)[0]
    M, L = periodic.mass_vector(X, S, lam)
    N = len(X)
    d = np.zeros(N)
    d = np.maximum(d, (X[:, 0] - (Nx - W)) / W)
    d = np.maximum(d, (W - X[:, 1]) / W)
    d = np.maximum(d, (X[:, 1] - (Ny - W)) / W)
    d = np.clip(d, 0, 1)
    if eta_max is None:
        eta_max = 2.0 * omega
    eta = np.repeat(eta_max * d ** 2, 3)
    A = (K - omega ** 2 * sp.diags(M) + 1j * omega * sp.diags(eta * M)).tocsc()
    left = np.where(X[:, 0] < 1e-9)[0]
    fixed = np.zeros(3 * N, bool); fixed[3 * left] = True
    u = np.zeros(3 * N, complex); u[3 * left] = 1.0
    free = ~fixed
    Aff = A[free][:, free].tocsc()
    rhs = -(A[free][:, fixed] @ u[fixed])
    u[free] = spl.spsolve(Aff, rhs)
    return u


def intensity(X, u):
    return np.abs(u[0::3]) ** 2 + np.abs(u[1::3]) ** 2


if __name__ == '__main__':
    TT = 200
    cal = json.load(open('results/lens_calib.json'))
    Nx, Ny, R = 12, 12, 3.0
    xc, yc = 1.5 + R, Ny / 2
    out = {}
    for rule in ('SR', 'balanced'):
        nb = NB[rule]
        Xl, Sl, ll, rmap = build_domain(rule, cal[rule], TT, Nx, Ny, xc, yc, R, True)
        Xr, Sr, lr, _ = build_domain(rule, cal[rule], TT, Nx, Ny, xc, yc, R, False)
        print(rule, 'nodes', len(Xl), 'rho range %.3f-%.3f' % (rmap.min(), rmap.max()), flush=True)
        res = []
        for lamb in (8.0, 6.0, 4.5, 3.5, 2.75, 2.25, 1.85, 1.5):     # background wavelength in tiles
            om = 2 * np.pi / (nb * lamb)
            t0 = time.time()
            ul = solve(Xl, Sl, ll, om, Nx, Ny); ur = solve(Xr, Sr, lr, om, Nx, Ny)
            Il = intensity(Xl, ul); Ir = intensity(Xr, ur)
            fx, fy = xc + R, yc
            sel_l = np.hypot(Xl[:, 0] - fx, Xl[:, 1] - fy) < 0.5
            sel_r = np.hypot(Xr[:, 0] - fx, Xr[:, 1] - fy) < 0.5
            G = float(np.mean(Il[sel_l]) / np.mean(Ir[sel_r]))
            # transverse profile on the focal line
            line = np.abs(Xl[:, 0] - fx) < 0.5
            ys = np.arange(2.75, Ny - 2.75, 0.5)
            prof = [float(np.mean(Il[line & (np.abs(Xl[:, 1] - y) < 0.25)])) for y in ys]
            res.append(dict(lambda_b=lamb, omega=om, D_over_lambda=2 * R / lamb, gain=G, ys=ys.tolist(), profile=prof))
            np.savez_compressed(f'results/lens_{rule}_{lamb}.npz', X=Xl, I=Il, Xr=Xr, Ir=Ir)
            print(rule, 'lambda_b %.2f omega %.3f D/lambda %.2f gain %.2f %.0fs' % (lamb, om, 2 * R / lamb, G, time.time() - t0), flush=True)
        out[rule] = dict(nb=nb, R=R, Nx=Nx, Ny=Ny, xc=xc, yc=yc, rho_map=rmap.tolist(), runs=res)
        json.dump(out, open('results/lens.json', 'w'))
