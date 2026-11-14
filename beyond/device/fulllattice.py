"""Full-lattice check of shock routing (Fig. 4) on a reduced device, against the homogenised model
on the identical setup.

Insert: Nx x Ny cells of unit size; cell (i, j) is an independent realisation of the cut-angle
division with T ligaments (rule (area, centroid, ('perpmix', kappa_i))), interior thickness set by
rho_t = 0.2; cell frames are the global grid lines (thickness = mean interior thickness), split at
every ligament end that lands on them, so neighbouring cells of different kappa join exactly.
Euler-Bernoulli frame (rg/frame.py), lumped mass, stiffness-proportional damping (loss eta at w_ref).
Bottom grid line fixed (u, v, rotation); half-sine force pulse (y) on the top grid line over |x| < a;
force on the component = vertical reaction on the bottom nodes with |x| < xc.
Homogenised twin: device/shock.run with the T-ligament cell tables (results/anglegrade_T{T}.json).
-> results/fulllattice_T{T}.json"""
import json, os, sys, time
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spl
from scipy.spatial import cKDTree
HERE = os.path.dirname(os.path.abspath(__file__)); B = os.path.join(HERE, '..')
sys.path.insert(0, B); sys.path.insert(0, os.path.join(B, '..', 'rg')); sys.path.insert(0, HERE)
import divisions, frame, abh, shock


def build(kappa_of_x, Nx, Ny, T, rho_t=0.2, seed0=0):
    P = []; S = []; LAM = []
    lam_all = []
    for i in range(Nx):
        xm = i + 0.5 - Nx / 2
        k = float(kappa_of_x(xm))
        for j in range(Ny):
            g = divisions.run(T, seed0 + 1000 * i + j, ('area', 'centroid', ('perpmix', k)))
            L = np.linalg.norm(g['B'] - g['A'], axis=1); lam = rho_t / L.sum()
            lam_all.append(lam)
            # ligament segments split at junctions (no frame): reuse srtree-style splitting
            import srtree
            nodes, seg, sl, bn = srtree.network(g)
            off = np.array([i - Nx / 2, j])
            base = sum(len(p) for p in P)
            P.append(nodes + off); S.append(seg + base); LAM.append(np.full(len(seg), lam))
    X = np.vstack(P); S = np.vstack(S); LAM = np.concatenate(LAM)
    lam_f = float(np.mean(lam_all))
    # frame: global grid lines, split at all nodes lying on them
    fr_nodes = []
    extra = []
    for gx in range(Nx + 1):
        x = gx - Nx / 2
        on = np.where(np.abs(X[:, 0] - x) < 1e-9)[0]
        ys = list(range(Ny + 1))
        extra += [(x, y) for y in ys]
        fr_nodes.append(('v', x, on))
    for gy in range(Ny + 1):
        on = np.where(np.abs(X[:, 1] - gy) < 1e-9)[0]
        extra += [(i - Nx / 2, gy) for i in range(Nx + 1)]
        fr_nodes.append(('h', gy, on))
    E = np.unique(np.array(extra), axis=0)
    base = len(X); X = np.vstack([X, E])
    # merge coincident nodes
    pairs = cKDTree(X).query_pairs(1e-9, output_type='ndarray')
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
    Xm = X[uniq]; S = gid[S]
    FS = []
    for gx in range(Nx + 1):
        x = gx - Nx / 2
        on = np.where(np.abs(Xm[:, 0] - x) < 1e-9)[0]
        on = on[np.argsort(Xm[on, 1])]
        FS += list(zip(on[:-1], on[1:]))
    for gy in range(Ny + 1):
        on = np.where(np.abs(Xm[:, 1] - gy) < 1e-9)[0]
        on = on[np.argsort(Xm[on, 0])]
        FS += list(zip(on[:-1], on[1:]))
    FS = np.array(FS)
    S = np.vstack([S, FS]); LAM = np.concatenate([LAM, np.full(len(FS), lam_f)])
    keep = S[:, 0] != S[:, 1]
    S, LAM = S[keep], LAM[keep]
    key = np.sort(S, axis=1); _, first = np.unique(key, axis=0, return_index=True)
    return Xm, S[first], LAM[first]


def simulate(Xm, S, LAM, Nx, Ny, a=3.0, xc=4.0, tau=80.0, t_end=400.0, dt=1.0, eta=0.03, w_ref=0.1):
    K = frame.assemble(Xm, S, LAM)[0]
    L = np.linalg.norm(Xm[S[:, 1]] - Xm[S[:, 0]], axis=1); m = LAM * L
    N = len(Xm); M = np.zeros(3 * N)
    for e in (0, 1):
        np.add.at(M, 3 * S[:, e], 0.5 * m); np.add.at(M, 3 * S[:, e] + 1, 0.5 * m)
        np.add.at(M, 3 * S[:, e] + 2, m * L ** 2 / 24)
    M = np.maximum(M, 1e-12 * M.max())
    Cd = (eta / w_ref) * K
    bottom = np.where(Xm[:, 1] < 1e-9)[0]
    fixed = np.zeros(3 * N, bool)
    for c in range(3): fixed[3 * bottom + c] = True
    free = ~fixed
    top = np.where((np.abs(Xm[:, 1] - Ny) < 1e-9) & (np.abs(Xm[:, 0]) < a))[0]
    fdir = np.zeros(3 * N); fdir[3 * top + 1] = -1.0 / len(top)
    Kf = K[free][:, free]; Cf = Cd[free][:, free]; Mf = sp.diags(M[free])
    t0 = time.time()
    lu = spl.splu((Kf + (2 / dt) * Cf + (4 / dt ** 2) * Mf).tocsc())
    print('  factorised %d dofs in %.0fs' % (free.sum(), time.time() - t0), flush=True)
    nf = free.sum(); u = np.zeros(nf); v = np.zeros(nf); acc = np.zeros(nf)
    fx = np.where(fixed)[0]; Kfx = K[fixed][:, free]; Cfx = Cd[fixed][:, free]
    comp = np.isin(fx, 3 * bottom[np.abs(Xm[bottom, 0]) < xc] + 1)
    ts = np.arange(0, t_end + dt / 2, dt); Fc = np.zeros(len(ts))
    F = lambda t: (np.sin(np.pi * t / tau) if t < tau else 0.0)
    for i, t in enumerate(ts[1:], start=1):
        rhs = fdir[free] * F(t) + Mf @ ((4 / dt ** 2) * u + (4 / dt) * v + acc) + Cf @ ((2 / dt) * u + v)
        un = lu.solve(rhs); vn = (2 / dt) * (un - u) - v; an = (4 / dt ** 2) * (un - u) - (4 / dt) * v - acc
        u, v, acc = un, vn, an
        Fc[i] = -(Kfx @ u + Cfx @ v)[comp].sum()
    return ts, Fc


if __name__ == '__main__':
    T = int(sys.argv[1]) if len(sys.argv) > 1 else 100
    Nx, Ny = int(sys.argv[2]) if len(sys.argv) > 2 else 40, int(sys.argv[3]) if len(sys.argv) > 3 else 20
    tau = float(sys.argv[4]) if len(sys.argv) > 4 else 80.0
    X5 = 5.0
    designs = {'U': lambda x: 0.0, 'D': lambda x: min(1.0, (abs(x) / X5) ** 2), 'C': lambda x: 1 - min(1.0, (abs(x) / X5) ** 2)}
    Tab, ratio = abh.tables(f'anglegrade_T{T}.json')
    out = {'T': T, 'Nx': Nx, 'Ny': Ny, 'tau': tau, 'full': {}, 'hom': {}}
    for name, kf in designs.items():
        t0 = time.time()
        Xm, S, LAM = build(kf, Nx, Ny, T)
        print(name, 'nodes', len(Xm), 'segs', len(S), 'built %.0fs' % (time.time() - t0), flush=True)
        ts, Fc = simulate(Xm, S, LAM, Nx, Ny, tau=tau)
        out['full'][name] = dict(t=ts.tolist(), F=Fc.tolist(), peak=float(np.abs(Fc).max()),
                                 impulse=float(np.trapezoid(np.clip(-Fc, 0, None), ts)))
        r = shock.run(Tab, ratio, name, Wd=float(Nx), Hd=float(Ny), h=0.25, X=X5, tau=tau, t_end=400.0, dt=1.0, w_ref=0.1)
        out['hom'][name] = dict(t=r['t'].tolist(), F=r['Fcomp'].tolist(), peak=float(np.abs(r['Fcomp']).max()),
                                impulse=float(np.trapezoid(np.clip(-r['Fcomp'], 0, None), r['t'])),
                                band=r['band'], frac_out=r['frac_out'])
        print(name, 'full peak %.4f  hom peak %.4f  band %.3f  spectrum above band %.3f  %.0fs' % (
            out['full'][name]['peak'], out['hom'][name]['peak'], r['band'], r['frac_out'], time.time() - t0), flush=True)
        json.dump(out, open(os.path.join(B, 'results', f'fulllattice_T{T}.json'), 'w'))
    for key in ('full', 'hom'):
        pu = out[key]['U']['peak']; iu = out[key]['U']['impulse']
        print(key, ' '.join('%s peak %.2f impulse %.2f |' % (n, out[key][n]['peak'] / pu, out[key][n]['impulse'] / iu) for n in designs))
