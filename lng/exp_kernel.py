"""Attachment kernel and rate equation (book: Local network growth, 'one rate
equation for every model').

Measures the mean-field rate r_k(t) at which a node of degree k is hit,
kappa(k) = t * <hits on degree-k nodes at step t> / n_k(t), exactly from the growth
history (n_k(t) tracked step by step), in 2D (SR ligaments) and 3D (SR plates).
A rate equation with kernel (a + b k)/t and birth-degree law q_k then gives
  p_k (1 + a + b k) = (a + b (k-1)) p_{k-1} + q_k,
which is compared with the simulated degree distribution.

usage: python exp_kernel.py [2d|3d|all]
"""
import os
import sys
import numpy as np
from numba import njit

HERE = os.path.dirname(os.path.abspath(__file__))
RES = os.path.join(HERE, 'results')
KMAX = 4000


@njit(cache=True)
def kernel_scan(T, birth_deg, ev_t, ev_h, t0, kmax):
    """birth_deg[t]: degree of node t at birth (node t born at step t).
    ev_t, ev_h: hits (step, host), sorted by step.  Over steps t >= t0 returns
    H[k] = number of hits on nodes of degree k and D[k] = sum_t n_k(t)/t (exact),
    so that kappa(k) = H[k]/D[k] is the mean-field kernel t*r_k."""
    deg = np.zeros(T, np.int64)
    cnt = np.zeros(kmax + 2, np.float64)
    H = np.zeros(kmax + 2, np.float64)
    D = np.zeros(kmax + 2, np.float64)
    last = np.zeros(kmax + 2, np.int64)
    # cumulative harmonic weights: hw[t] = sum_{t0<=s<t} 1/(s+1)
    hw = np.zeros(T + 1, np.float64)
    for s in range(T):
        hw[s + 1] = hw[s] + (1.0 / (s + 1.0) if s >= t0 else 0.0)
    e = 0
    ne = len(ev_t)
    for t in range(T):
        e1 = e
        while e < ne and ev_t[e] == t:
            k = deg[ev_h[e]]
            if t >= t0 and k <= kmax:
                H[k] += 1.0
            e += 1
        for j in range(e1, e):
            h = ev_h[j]
            k = deg[h]
            for kk in (k, k + 1):
                if kk <= kmax:
                    D[kk] += cnt[kk] * (hw[t + 1] - hw[last[kk]])
                    last[kk] = t + 1
            if k <= kmax:
                cnt[k] -= 1
            if k + 1 <= kmax:
                cnt[k + 1] += 1
            deg[h] = k + 1
        b = birth_deg[t]
        deg[t] = b
        if b <= kmax:
            D[b] += cnt[b] * (hw[t + 1] - hw[last[b]])
            last[b] = t + 1
            cnt[b] += 1
    for kk in range(kmax + 1):
        D[kk] += cnt[kk] * (hw[T] - hw[last[kk]])
    return H, D, deg


def rate_equation_pk(a, b, q, kmax=KMAX):
    """stationary p_k of n_k(t+1)=n_k+[(a+b(k-1))n_{k-1}-(a+bk)n_k]/t + q_k."""
    p = np.zeros(kmax + 1)
    for k in range(kmax + 1):
        prev = (a + b * (k - 1)) * p[k - 1] if k >= 1 else 0.0
        qk = q[k] if k < len(q) else 0.0
        p[k] = (prev + qk) / (1.0 + a + b * k)
    return p


def ccdf(p):
    return p[::-1].cumsum()[::-1]


def fit_kernel(kappa, w, kmin, kmax):
    """weighted linear fit kappa(k) = a + b k over kmin<=k<=kmax."""
    k = np.arange(len(kappa))
    m = (k >= kmin) & (k <= kmax) & (w > 0)
    W = w[m]
    A = np.vstack([np.ones(m.sum()), k[m]]).T
    coef, *_ = np.linalg.lstsq(A * np.sqrt(W)[:, None], kappa[m] * np.sqrt(W), rcond=None)
    return coef


def run2d(seeds=(1, 2, 3), T=1_000_000):
    sys.path.insert(0, HERE)
    import srcells
    rows = []
    for s in seeds:
        o = srcells.run(T, lam0=0.0, seed=s)
        h = o['hit']
        Tn = len(h)
        birth = (h >= 0).sum(1)
        ev_t = np.repeat(np.arange(Tn), 2)
        ev_h = h.reshape(-1)
        m = ev_h >= 0
        ev_t, ev_h = ev_t[m].astype(np.int64), ev_h[m].astype(np.int64)
        H, D, deg = kernel_scan(Tn, birth.astype(np.int64), ev_t, ev_h, Tn // 10, KMAX)
        rows.append(dict(H=H, D=D, deg=deg, birth=birth, nhit=len(ev_t) / Tn))
    return rows


def run3d(names=('iso_s1', 'iso_s2', 'iso_s3')):
    sys.path.insert(0, os.path.join(HERE, '..', 'd3'))
    import common
    rows = []
    for nm in names:
        S_, E, C = common.load(nm)
        P = common.plate_arrays(S_, E)
        T = P['T'] + 1                     # plates labelled 1..T-1, index 0 unused
        birth = P['out'].astype(np.int64)
        birth[~P['ok']] = 0
        order = np.argsort(E['t'], kind='stable')
        ev_t = E['t'][order].astype(np.int64)
        ev_h = E['h'][order].astype(np.int64)
        m = ev_h >= 0
        ev_t, ev_h = ev_t[m], ev_h[m]
        H, D, deg = kernel_scan(T, birth, ev_t, ev_h, T // 10, KMAX)
        rows.append(dict(H=H, D=D, deg=deg[P['ok']], birth=birth[P['ok']],
                         nhit=len(ev_t) / P['ok'].sum()))
    return rows


def report(tag, rows, b_book, out):
    print(f'== {tag}', file=out)
    kap_all, fits = [], []
    for r in rows:
        kap = np.where(r['D'] > 0, r['H'] / np.maximum(r['D'], 1e-300), 0.0)
        w = np.where(kap > 0, r['H'] / np.maximum(kap, 1e-300) ** 2, 0.0)  # 1/var (Poisson H)
        kap_all.append(kap)
        nb = np.bincount(r['birth'], minlength=KMAX + 1)[:KMAX + 1] / len(r['birth'])
        obar = (np.arange(KMAX + 1) * nb).sum()
        ab = fit_kernel(kap, w, 2 if tag == '2D' else 4, 200)
        fits.append((ab, obar, nb, r))
        print(f'  realisation: hits/step {r["nhit"]:.3f}  <o> {obar:.3f}  '
              f'fit kappa = {ab[0]:+.3f} + {ab[1]:.4f} k  (k<=200)', file=out)
    kap = sum(r['H'] for r in rows) / np.maximum(sum(r['D'] for r in rows), 1e-300)
    Hs = sum(r['H'] for r in rows)
    ks = [2, 3, 4, 6, 8, 12, 16, 24, 32, 48, 64, 96, 128]
    print('  k      kappa(k)   kappa/k   +-      hits', file=out)
    for k in ks:
        if Hs[k] > 0:
            print(f'  {k:4d}  {kap[k]:8.3f}  {kap[k]/k:7.4f}  {kap[k]/k/np.sqrt(Hs[k]):.4f}  {int(Hs[k])}', file=out)
    # figure data: pooled simulated CCDF and rate equation with the deduced b
    degs = np.concatenate([f[3]['deg'] for f in fits])
    nb = np.bincount(np.concatenate([f[3]['birth'] for f in fits]), minlength=KMAX + 1)[:KMAX + 1]
    nb = nb / nb.sum(); obar = (np.arange(KMAX + 1) * nb).sum()
    pre = ccdf(rate_equation_pk(obar * (1 - 2 * b_book), b_book, nb))
    Kg = np.unique(np.logspace(np.log10(2), np.log10(min(degs.max(), KMAX)), 40).astype(int))
    np.savez(os.path.join(RES, f'kernel_{tag}.npz'), kappa=kap, H=Hs, K=Kg,
             Psim=np.array([(degs >= x).mean() for x in Kg]), Pre=pre[Kg], obar=obar, b=b_book)
    print(f'  deduced constants: b={b_book:.4f}, a=<o>(1-2b)={obar*(1-2*b_book):+.3f}', file=out)
    # rate equations: (i) deduced constants, (ii) fitted kernel; birth law measured
    for label, bsel in (('deduced b', b_book), ('fitted b', None)):
        rat = {K: [] for K in (8, 16, 32, 64, 128, 256)}
        for ab, obar, nb, r in fits:
            b = ab[1] if bsel is None else bsel
            a = obar * (1 - 2 * b)         # sum of rates = hits per step
            p = rate_equation_pk(a, b, nb)
            P = ccdf(p)
            deg = r['deg']
            for K in rat:
                sim = (deg >= K).mean()
                rat[K].append(sim / P[K] if P[K] > 0 else np.nan)
        bb = b_book if bsel is not None else np.mean([f[0][1] for f in fits])
        s = '  '.join(f'K={K}: {np.mean(v):.3f}' for K, v in rat.items())
        print(f'  rate equation [{label}={bb:.4f}, gamma=1+1/b={1+1/bb:.3f}]  '
              f'sim/theory CCDF: {s}', file=out)


def main(which):
    os.makedirs(RES, exist_ok=True)
    with open(os.path.join(RES, 'kernel.txt'), 'a') as out:
        if which in ('2d', 'all'):
            report('2D', run2d(), 0.5, out)
        if which in ('3d', 'all'):
            report('3D', run3d(), 2 / 3, out)
    print(open(os.path.join(RES, 'kernel.txt')).read())


if __name__ == '__main__':
    main(sys.argv[1] if len(sys.argv) > 1 else 'all')
