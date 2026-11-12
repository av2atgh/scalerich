"""Exact generation flow of the conductance (Goal 2a/3a).

Adding chord t (conductance g_t = lam_t/l_t between boundary points P, Q of
its cell) to the network of generations < t changes the dissipated power at
unit applied field by the rank-one (Thevenin / Sherman-Morrison) increment

    dP_t = g_t dV_t^2 / (1 + g_t R_t),

dV_t = potential difference between P and Q before the chord, R_t = two-point
resistance between P and Q with the electrodes grounded.  Hence exactly

    sigma(T) = sum_t dP_t = sum_t psi_t * aff_t,  aff_t = lam_t l_t cos^2 theta_t.

psi_t factorises into a field-transmission factor F_t = g dV^2/aff_t (how
much of the applied field the cell frame passes to the chord) and the
Thevenin factor 1/(1+g R).  The RG statement: the coupling g_t R_t is
stationary in ln t for alpha < 1/2 (marginal, set by alpha) and flows to 0
for alpha > 1/2; psi(t) converges to a fixed-line value psi(alpha).

Measures, at log-spaced checkpoints tau of one realisation:
  sigma(tau), sigma_aff(tau)    (x field)
  for chords tau..tau+K-1: dV, R, g, dP (checked against sigma increments)
Writes results/flow.json
"""
import json, sys, time
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spl
import srtree


def trunc(g, n):
    h = dict(T=n)
    for k in ('A', 'B', 'hostA', 'hostB'):
        h[k] = g[k][:n]
    return h


def lap(N, seg, gs):
    i, j = seg[:, 0], seg[:, 1]
    L = sp.coo_matrix((np.r_[-gs, -gs], (np.r_[i, j], np.r_[j, i])), shape=(N, N)).tocsr()
    return L - sp.diags(np.asarray(L.sum(1)).ravel())


def analyse(g, n, lam, K):
    """network of the first n+K ligaments; returns sigma(n'), for n'=n..n+K and
    Thevenin data for chords n..n+K-1."""
    h = trunc(g, n + K)
    nodes, seg, sl, bn = srtree.network(h)
    N = len(nodes)
    Ls = np.linalg.norm(nodes[seg[:, 0]] - nodes[seg[:, 1]], axis=1)
    fixed = np.zeros(N, bool); fixed[bn[-4]] = True; fixed[bn[-2]] = True
    Vb = np.zeros(N); Vb[bn[-2]] = 1.0
    rows = []
    sig = []
    for m in range(K + 1):
        # chords with index >= n+m are switched off
        on = sl < n + m
        gs = np.where(on, lam[sl] / Ls, 0.0)
        used = np.zeros(N, bool); used[seg[on, 0]] = True; used[seg[on, 1]] = True
        # chord endpoints on a ligament are already in the system; an endpoint on
        # a free side (-1, -3) is isolated: R = inf, dP = 0
        free = used & ~fixed
        L = lap(N, seg, gs)
        Lff = L[free][:, free].tocsc()
        lu = spl.splu(Lff)
        V = Vb.copy()
        V[free] = lu.solve(-(L[free][:, fixed] @ Vb[fixed]))
        dV = V[seg[:, 0]] - V[seg[:, 1]]
        sig.append(float(np.sum(gs * dV ** 2)))
        if m < K:
            c = n + m
            fidx = -np.ones(N, int); fidx[free] = np.arange(free.sum())
            b = np.zeros(free.sum())
            pa, pb = 2 * c, 2 * c + 1
            if fidx[pa] >= 0 and used[pa]: b[fidx[pa]] += 1
            if fidx[pb] >= 0 and used[pb]: b[fidx[pb]] -= 1
            x = lu.solve(b)
            xa = x[fidx[pa]] if fidx[pa] >= 0 and used[pa] else 0.0
            xb = x[fidx[pb]] if fidx[pb] >= 0 and used[pb] else 0.0
            R = xa - xb
            if not ((used[pa] or fixed[pa]) and (used[pb] or fixed[pb])):
                R = np.inf
            d = g['B'][c] - g['A'][c]; l = np.linalg.norm(d)
            gc = lam[c] / l
            dVc = V[pa] - V[pb] if np.isfinite(R) else 0.0
            rows.append(dict(t=c, gc=gc, R=float(R), dV=float(dVc), aff=float(lam[c] * d[0] ** 2 / l),
                             dP=float(gc * dVc ** 2 / (1 + gc * R)) if np.isfinite(R) else 0.0))
    return sig, rows


if __name__ == '__main__':
    T = int(sys.argv[1]) if len(sys.argv) > 1 else 20000
    NS = int(sys.argv[2]) if len(sys.argv) > 2 else 4
    K = 4
    alphas = [0.0, 0.2, 0.4, 0.5, 0.6, 0.8, 1.0]
    taus = np.unique(np.round(np.logspace(1, np.log10(T - K), 16)).astype(int))
    out = dict(T=T, NS=NS, K=K, taus=taus.tolist(), data={})
    for s in range(NS):
        g = srtree.run(T, 900 + s)
        for a in alphas:
            lam = srtree.thickness(T, 1.0, a)
            t0 = time.time()
            recs = []
            for tau in taus:
                try:
                    sig, rows = analyse(g, int(tau), lam, K)
                except Exception as e:
                    print('FAILED', s, a, tau, e); continue
                inc = np.diff(sig)
                for r_, di in zip(rows, inc):
                    r_['dP_direct'] = float(di)
                recs.append(dict(tau=int(tau), sigma=sig[0], rows=rows))
            out['data'].setdefault(str(a), []).append(recs)
            print('seed', s, 'alpha', a, 'time %.1f' % (time.time() - t0), flush=True)
        json.dump(out, open('results/flow.json', 'w'))
