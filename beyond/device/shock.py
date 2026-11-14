"""Shock routing in a 2.5D division-lattice insert protecting a component (homogenised, time domain).

Slab of width Wd and depth Hd (cells), on a rigid base (bottom fixed); the component occupies the
bottom centre |x| < xc. A half-sine force pulse of duration tau acts on the top surface over
|x| < a, directly above the component. Cells: cut-angle family at fixed density rho = 0.2
(kappa = 0 'balanced' ... 1 'centroid'; n = 2.94 ... 4.41), isotropised K, G, rho from
results/anglegrade.json (abh.cell). Designs (kappa field):
  U  uniform kappa = 0                      S  uniform kappa = 1
  D  diverging: kappa = min(1, (|x|/X)^2)   (low index above the component, high at the sides)
  C  converging: kappa = 1 - min(1, (|x|/X)^2)   (control)
Time integration: Newmark average acceleration, lumped mass, stiffness-proportional damping
(loss eta at omega_ref). Output: force transmitted to the component (bottom reaction over |x| < xc),
its peak, arrival time (first crossing of 10% of the U peak), impulse; total base force; snapshots.
Validity: pulse spectrum vs the band edge of the cells used (0.8 omega1, abh.tables ratio).
Units: E0 = rho0 = 1, lengths in cells L, time in L/c0.
"""
import json, os, sys, time
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spl
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import abh

B = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..')


def kappa_field(design, x, X):
    s = np.minimum(1.0, (np.abs(x) / X) ** 2)
    return {'U': 0 * x, 'S': 1 + 0 * x, 'D': s, 'C': 1 - s}[design]


def run(T, ratio, design, Wd=60.0, Hd=30.0, h=0.5, X=15.0, a=3.0, xc=4.0, tau=40.0, t_end=500.0,
        dt=0.5, eta=0.03, w_ref=0.15, snaps=(60, 100, 140, 180), fix_rot=None):
    Nx = int(round(Wd / h)); Ny = int(round(Hd / h))
    xc_e = (np.arange(Nx) + 0.5) * h - Wd / 2
    kap = kappa_field(design, xc_e, X)
    K_, G_, rho_, w1 = abh.cell(T, kap, 0.2 + 0 * kap)
    KK, KG = abh.q4_refs(h)
    nid = lambda i, j: i * (Ny + 1) + j
    I, J = np.meshgrid(np.arange(Nx), np.arange(Ny), indexing='ij')
    nodes = np.stack([nid(I, J), nid(I + 1, J), nid(I + 1, J + 1), nid(I, J + 1)], -1).reshape(-1, 4)
    dof = np.repeat(2 * nodes, 2, axis=1) + np.tile([0, 1], 4)
    Ke = np.repeat(K_, Ny)[:, None, None] * KK + np.repeat(G_, Ny)[:, None, None] * KG
    nn = (Nx + 1) * (Ny + 1); nd = 2 * nn
    K = sp.coo_matrix((Ke.ravel(), (np.repeat(dof, 8, axis=1).ravel(), np.tile(dof, (1, 8)).ravel())), shape=(nd, nd)).tocsr()
    me = np.repeat(rho_, Ny) * h * h / 4
    M = np.repeat(np.bincount(nodes.ravel(), weights=np.repeat(me, 4), minlength=nn), 2)
    Cd = (eta / w_ref) * K
    X_ = np.repeat(np.arange(Nx + 1) * h - Wd / 2, Ny + 1); Y_ = np.tile(np.arange(Ny + 1) * h, Nx + 1)
    bottom = np.where(Y_ < 1e-9)[0]
    fixed = np.zeros(nd, bool); fixed[2 * bottom] = True; fixed[2 * bottom + 1] = True
    free = ~fixed
    comp = bottom[np.abs(X_[bottom]) < xc]
    top = np.where((np.abs(Y_ - Hd) < 1e-9) & (np.abs(X_) < a))[0]
    fdir = np.zeros(nd); fdir[2 * top + 1] = -1.0 / len(top)          # downward, unit total force
    Kf = K[free][:, free]; Cf = Cd[free][:, free]; Mf = sp.diags(M[free])
    A = (Kf + (2 / dt) * Cf + (4 / dt ** 2) * Mf).tocsc()
    lu = spl.splu(A)
    nf = free.sum()
    u = np.zeros(nf); v = np.zeros(nf); acc = np.zeros(nf)
    Kfx = K[fixed][:, free]; Cfx = Cd[fixed][:, free]
    fx_idx = np.where(fixed)[0]
    comp_y = np.isin(fx_idx, 2 * comp + 1)
    all_y = (fx_idx % 2 == 1)
    ts = np.arange(0, t_end + dt / 2, dt)
    Fcomp = np.zeros(len(ts)); Fbase = np.zeros(len(ts))
    F = lambda t: (np.sin(np.pi * t / tau) if t < tau else 0.0)
    snapshots = {}
    for i, t in enumerate(ts[1:], start=1):
        f = fdir[free] * F(t)
        rhs = f + Mf @ ((4 / dt ** 2) * u + (4 / dt) * v + acc) + Cf @ ((2 / dt) * u + v)
        un = lu.solve(rhs)
        vn = (2 / dt) * (un - u) - v
        an = (4 / dt ** 2) * (un - u) - (4 / dt) * v - acc
        u, v, acc = un, vn, an
        R = Kfx @ u + Cfx @ v                                         # reaction on fixed dofs
        Fcomp[i] = -R[comp_y].sum(); Fbase[i] = -R[all_y].sum()
        if any(abs(t - s_) < dt / 2 for s_ in snaps):
            uf = np.zeros(nd); uf[free] = u
            snapshots[int(round(t))] = np.hypot(uf[0::2], uf[1::2]).reshape(Nx + 1, Ny + 1)
    # spectrum check: fraction of pulse energy above the lowest band edge used
    band = ratio * float(w1.min())
    om = np.linspace(1e-3, 3.0, 3000)
    tt = np.linspace(0, tau, 400)
    S = np.abs(np.trapezoid(np.sin(np.pi * tt / tau)[None, :] * np.exp(-1j * om[:, None] * tt[None, :]), tt, axis=1)) ** 2
    frac_out = float(S[om > band].sum() / S.sum())
    return dict(t=ts, Fcomp=Fcomp, Fbase=Fbase, band=band, frac_out=frac_out, snaps=snapshots,
                n_min=float(np.sqrt(rho_ / (K_ + G_)).min()), n_max=float(np.sqrt(rho_ / (K_ + G_)).max()),
                mass=float(rho_.sum() * h * Hd))


if __name__ == '__main__':
    T, ratio = abh.tables()
    tau = float(sys.argv[1]) if len(sys.argv) > 1 else 40.0
    out = {}
    ref_peak = None
    for design in ('U', 'S', 'D', 'C'):
        t0 = time.time()
        r = run(T, ratio, design, tau=tau)
        Fc = r['Fcomp']; pk = float(np.abs(Fc).max())
        if ref_peak is None: ref_peak = pk
        arr = float(r['t'][np.argmax(np.abs(Fc) > 0.1 * ref_peak)])
        imp = float(np.trapezoid(np.clip(-Fc, 0, None), r['t']))   # compressive impulse
        out[design] = dict(peak=pk, arrival=arr, impulse_pos=imp, band=r['band'], frac_out=r['frac_out'],
                           n_range=[r['n_min'], r['n_max']], mass=r['mass'],
                           t=r['t'][::4].tolist(), Fcomp=Fc[::4].tolist(), Fbase=r['Fbase'][::4].tolist())
        np.savez_compressed(os.path.join(B, 'results', f'shock_{design}_tau{tau:g}.npz'),
                            **{f's{k}': v for k, v in r['snaps'].items()})
        print(design, 'peak %.4f (%.2f of U)  arrival %.1f  impulse+ %.3f  band %.3f  spectrum above band %.3f  n %.2f-%.2f  %.0fs' % (
            pk, pk / ref_peak, arr, imp, r['band'], r['frac_out'], r['n_min'], r['n_max'], time.time() - t0), flush=True)
    json.dump(out, open(os.path.join(B, 'results', f'shock_tau{tau:g}.json'), 'w'))
