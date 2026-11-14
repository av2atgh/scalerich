"""Shock routing in a full lattice: one recursive division of the whole insert (no cell frames,
no homogenisation).

Insert W x H (units: the area holding n_per ligaments is 1), divided by the cut-angle rule
(area, centroid, ('perpmix_field', kappa(x, y))) with W*H*n_per ligaments; the outer boundary is a
skin of the same thickness. All ligaments have the same thickness, set by rho = 0.2.
Euler-Bernoulli frame, lumped mass, stiffness-proportional damping (loss eta at w_ref).
Bottom skin fixed (rigid base); half-sine force pulse (y) on the top skin over |x| < a; force on the
component = vertical reaction on bottom nodes with |x| < xc. Designs: kappa fields
  U uniform 0 | S uniform 1 | D min(1, (|x|/X)^2) (diverging) | C 1 - min(1, (|x|/X)^2) (converging).
Two realisations per design. Snapshots: mass-weighted |u| in 0.5 x 0.5 bins.
-> results/globalshock.json, results/globalshock_<design>_<seed>.npz
"""
import json, os, sys, time
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spl
HERE = os.path.dirname(os.path.abspath(__file__)); B = os.path.join(HERE, '..')
sys.path.insert(0, B); sys.path.insert(0, os.path.join(B, '..', 'rg'))
import divisions, srtree, frame


def build(kfield, W, H, n_per, seed, rho=0.2):
    T = int(W * H * n_per)
    g = divisions.run(T, seed, ('area', 'centroid', ('perpmix_field', kfield)), box=(W, H))
    nodes, seg, sl, bn = srtree.network(g, frame=True)
    L = np.linalg.norm(nodes[seg[:, 1]] - nodes[seg[:, 0]], axis=1)
    lam = rho * W * H / L.sum()
    X = nodes.copy(); X[:, 0] -= W / 2
    return X, seg, np.full(len(seg), lam), lam


def simulate(X, S, LAM, H, a=3.0, xc=4.0, tau=80.0, t_end=350.0, dt=1.0, eta=0.03, w_ref=0.1,
             snaps=(), W=40.0, b=0.5):
    K = frame.assemble(X, S, LAM)[0]
    L = np.linalg.norm(X[S[:, 1]] - X[S[:, 0]], axis=1); m = LAM * L
    N = len(X); M = np.zeros(3 * N)
    for e in (0, 1):
        np.add.at(M, 3 * S[:, e], 0.5 * m); np.add.at(M, 3 * S[:, e] + 1, 0.5 * m)
        np.add.at(M, 3 * S[:, e] + 2, m * L ** 2 / 24)
    M = np.maximum(M, 1e-12 * M.max())
    Cd = (eta / w_ref) * K
    bottom = np.where(X[:, 1] < 1e-9)[0]
    fixed = np.zeros(3 * N, bool)
    for c in range(3): fixed[3 * bottom + c] = True
    free = ~fixed
    top = np.where((np.abs(X[:, 1] - H) < 1e-9) & (np.abs(X[:, 0]) < a))[0]
    fdir = np.zeros(3 * N); fdir[3 * top + 1] = -1.0 / len(top)
    Kf = K[free][:, free]; Cf = Cd[free][:, free]; Mf = sp.diags(M[free])
    lu = spl.splu((Kf + (2 / dt) * Cf + (4 / dt ** 2) * Mf).tocsc())
    nf = free.sum(); u = np.zeros(nf); v = np.zeros(nf); acc = np.zeros(nf)
    fx = np.where(fixed)[0]; Kfx = K[fixed][:, free]; Cfx = Cd[fixed][:, free]
    comp = np.isin(fx, 3 * bottom[np.abs(X[bottom, 0]) < xc] + 1)
    ts = np.arange(0, t_end + dt / 2, dt); Fc = np.zeros(len(ts))
    F = lambda t: (np.sin(np.pi * t / tau) if t < tau else 0.0)
    mn = M[0::3]
    ix = np.clip(((X[:, 0] + W / 2) / b).astype(int), 0, int(W / b) - 1)
    iy = np.clip((X[:, 1] / b).astype(int), 0, int(H / b) - 1)
    idx = ix * int(H / b) + iy; nb = int(W / b) * int(H / b)
    wsum = np.bincount(idx, weights=mn, minlength=nb)
    images = {}
    for i, t in enumerate(ts[1:], start=1):
        rhs = fdir[free] * F(t) + Mf @ ((4 / dt ** 2) * u + (4 / dt) * v + acc) + Cf @ ((2 / dt) * u + v)
        un = lu.solve(rhs); vn = (2 / dt) * (un - u) - v; an = (4 / dt ** 2) * (un - u) - (4 / dt) * v - acc
        u, v, acc = un, vn, an
        Fc[i] = -(Kfx @ u + Cfx @ v)[comp].sum()
        if any(abs(t - s_) < dt / 2 for s_ in snaps):
            uf = np.zeros(3 * N); uf[free] = u
            amp = np.hypot(uf[0::3], uf[1::3])
            img = np.bincount(idx, weights=mn * amp, minlength=nb) / np.where(wsum > 0, wsum, 1)
            images[int(round(t))] = img.reshape(int(W / b), int(H / b))
    return ts, Fc, images


if __name__ == '__main__':
    W, H, n_per = 40.0, 20.0, 100
    tau = float(sys.argv[1]) if len(sys.argv) > 1 else 80.0
    Xg = 5.0
    fields = {'U': lambda x, y: 0.0, 'S': lambda x, y: 1.0,
              'D': lambda x, y: min(1.0, (abs(x - W / 2) / Xg) ** 2),
              'C': lambda x, y: 1.0 - min(1.0, (abs(x - W / 2) / Xg) ** 2)}
    OUT = os.path.join(B, 'results', 'globalshock.json')
    out = json.load(open(OUT)) if os.path.exists(OUT) else {}
    for seed in (1, 2):
        for name in ('U', 'D', 'C', 'S'):
            key = f'{name}_{seed}'
            if key in out:
                continue
            t0 = time.time()
            X, S, LAM, lam = build(fields[name], W, H, n_per, seed)
            tb = time.time() - t0
            ts, Fc, imgs = simulate(X, S, LAM, H, tau=tau, snaps=(60, 90, 120), W=W)
            out[key] = dict(design=name, seed=seed, t=ts.tolist(), F=Fc.tolist(), peak=float(np.abs(Fc).max()),
                            impulse=float(np.trapezoid(np.clip(-Fc, 0, None), ts)), lam=lam, nodes=len(X))
            np.savez_compressed(os.path.join(B, 'results', f'globalshock_{key}.npz'), **{f's{k}': v for k, v in imgs.items()})
            json.dump(out, open(OUT, 'w'))
            print(key, 'nodes %d  peak %.4f  impulse %.3f  build %.0fs total %.0fs' % (
                len(X), out[key]['peak'], out[key]['impulse'], tb, time.time() - t0), flush=True)
    for seed in (1, 2):
        pu = out[f'U_{seed}']['peak']; iu = out[f'U_{seed}']['impulse']
        print('seed', seed, ' '.join('%s peak %.2f impulse %.2f |' % (n, out[f'{n}_{seed}']['peak'] / pu,
                                                                       out[f'{n}_{seed}']['impulse'] / iu) for n in 'UDCS'))
