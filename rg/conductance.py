"""Effective conductance of SR networks vs the affine (stiff-frame) value (Goal 3a).

Ligament segment conductance g = lam/len (unit sheet conductivity, unit depth).
Left side V=0, right side V=1 (and the y analogue); top/bottom free.
sigma_x = dissipated power; affine sigma_aff_x = sum_lig lam*l*cos^2(theta).
Phi = sigma/sigma_aff.  Per-generation phi(t): dissipation / affine in log bins
of birth index.

Generation-flow (Thevenin) prediction, see NOTES.md:
  per-chord phi = 1/(1 + kappa*g_c/sigma(t)),  g_c = lam_t/l_t,
  alpha<1/2: Phi = 1 - kappa (1-2 alpha)/c_l^2  (T-independent),
  alpha>1/2: phi(t) -> 1 at late t.
Writes results/conductance.json
"""
import json, sys, time
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spl
import srtree


def solve(nodes, seg, g, left, right):
    N = len(nodes)
    i, j = seg[:, 0], seg[:, 1]
    L = sp.coo_matrix((np.r_[-g, -g], (np.r_[i, j], np.r_[j, i])), shape=(N, N)).tocsr()
    L = L - sp.diags(np.asarray(L.sum(1)).ravel())
    V = np.full(N, np.nan); V[left] = 0.0; V[right] = 1.0
    fixed = ~np.isnan(V); free = ~fixed
    used = np.zeros(N, bool); used[i] = True; used[j] = True
    free &= used
    Lff = L[free][:, free]; Lfb = L[free][:, fixed]
    V[free] = spl.spsolve(Lff.tocsc(), -Lfb @ V[fixed])
    V[~used] = 0
    dV = V[i] - V[j]
    return V, g * dV ** 2


def run_one(T, seed, alphas):
    g0 = srtree.run(T, seed)
    nodes, seg, sl, bn = srtree.network(g0)
    Ls = np.linalg.norm(nodes[seg[:, 0]] - nodes[seg[:, 1]], axis=1)
    d = g0['B'] - g0['A']; l = np.linalg.norm(d, axis=1); c2 = (d[:, 0] / l) ** 2
    out = {}
    bins = np.unique(np.r_[0, np.round(np.logspace(0, np.log10(T), 10)).astype(int)])
    for a in alphas:
        lam = srtree.thickness(T, 1.0, a)
        g = lam[sl] / Ls
        r = {}
        for ax, (lf, rt), cc in (('x', (bn[-4], bn[-2]), c2), ('y', (bn[-1], bn[-3]), 1 - c2)):
            V, P = solve(nodes, seg, g, lf, rt)
            Pl = np.bincount(sl, P, minlength=T)
            aff = lam * l * cc
            phib = [float(Pl[b0:b1].sum() / aff[b0:b1].sum()) for b0, b1 in zip(bins[:-1], bins[1:])]
            r[ax] = dict(sigma=float(P.sum()), aff=float(aff.sum()), phi_bins=phib)
        out[a] = r
    return out, bins.tolist()


if __name__ == '__main__':
    Ts = [int(x) for x in sys.argv[1].split(',')] if len(sys.argv) > 1 else [1000, 3000, 10000, 30000]
    NS = int(sys.argv[2]) if len(sys.argv) > 2 else 4
    alphas = [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.8, 1.0]
    res = {}
    for T in Ts:
        acc = {a: [] for a in alphas}; phis = {a: [] for a in alphas}
        t0 = time.time()
        for s in range(NS):
            try:
                o, bins = run_one(T, 500 + s, alphas)
            except Exception as e:
                print('FAILED', T, s, e); continue
            for a in alphas:
                for ax in 'xy':
                    acc[a].append(o[a][ax]['sigma'] / o[a][ax]['aff'])
                    phis[a].append(o[a][ax]['phi_bins'])
        res[T] = {str(a): dict(Phi=float(np.mean(acc[a])), Phi_se=float(np.std(acc[a]) / np.sqrt(len(acc[a]))),
                               phi_bins=np.nanmean(np.array(phis[a]), 0).round(4).tolist(), bins=bins)
                  for a in alphas}
        print(T, 'time %.1f' % (time.time() - t0), {a: round(res[T][str(a)]['Phi'], 4) for a in alphas}, flush=True)
        json.dump(res, open('results/conductance.json', 'w'), indent=1)
