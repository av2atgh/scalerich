"""Which mechanism sets the 3D degree law: reinforcement (incidence PA) or a
fitness x age law (hit rate set by the plate area, fixed at birth)?

Three networks with the same newcomer out-degrees m_t (plate hits of the SR run):
  SR    the 3D plate simulation (data/<run>);
  FIT   fitness surrogate: each newcomer picks m_t distinct old plates with
        probability proportional to their birth area a_i (no reinforcement);
  PA    incidence-PA surrogate: picks m_t distinct old plates with probability
        proportional to s_i = 2 + hits received (pure reinforcement), b = m/(2+m).
Tests on each:
  (a) degree vs birth area x age (area-rule prediction);
  (b) dispersion of received hits around the best available predicted mean
      (Poisson = 1 for a fitness law; grows with the mean for PA);
  (c) reinforcement: partial correlation of hits in (2i,T] with hits in (i,2i]
      at fixed birth area, within cohorts;
  (d) fixed-plate growth exponent, MLE vs k_min and K^1.5 P(k>=K): shape
      comparison against the two controls instead of a bare tail slope.

usage: python3 exp_mechanism.py run -> results/mechanism_<run>.txt
"""
import sys, os
import numpy as np
from numba import njit
from common import load, plate_arrays, gamma_mle, RES


@njit(cache=True)
def grow_pa(mt, T, seed):
    np.random.seed(seed)
    tot = 2 * T + int(mt.sum()) + 10
    lst = np.empty(tot, np.int64); n = 0
    ht = np.empty(int(mt.sum()) + 10, np.int64); hh = np.empty_like(ht); ne = 0
    tg = np.empty(64, np.int64)
    for t in range(1, T + 1):
        m = mt[t]
        nd = t - 1
        if m > nd:
            m = nd
        c = 0
        tries = 0
        while c < m and tries < 10000:
            tries += 1
            x = lst[np.random.randint(0, n)]
            dup = False
            for j in range(c):
                if tg[j] == x:
                    dup = True
            if not dup:
                tg[c] = x; c += 1
        for j in range(c):
            ht[ne] = t; hh[ne] = tg[j]; ne += 1
            lst[n] = tg[j]; n += 1
        lst[n] = t; lst[n + 1] = t; n += 2
    return ht[:ne], hh[:ne]


@njit(cache=True)
def grow_fit(mt, a, T, seed):
    np.random.seed(seed)
    cum = np.cumsum(a)          # cum[i] = sum_{j<=i} a_j, a[0] = 0
    ht = np.empty(int(mt.sum()) + 10, np.int64); hh = np.empty_like(ht); ne = 0
    tg = np.empty(64, np.int64)
    for t in range(2, T + 1):
        m = mt[t]
        if m > t - 1:
            m = t - 1
        c = 0; tries = 0
        tot = cum[t - 1]
        while c < m and tries < 10000:
            tries += 1
            x = np.searchsorted(cum[:t], np.random.random() * tot, side='right')
            if x < 1:
                x = 1
            if x > t - 1:
                x = t - 1
            dup = False
            for j in range(c):
                if tg[j] == x:
                    dup = True
            if not dup:
                tg[c] = x; c += 1
        for j in range(c):
            ht[ne] = t; hh[ne] = tg[j]; ne += 1
    return ht[:ne], hh[:ne]


def analyse(label, Et, Eh, a, mt_plate, T, out, pred_area=True):
    recv = np.bincount(Eh, minlength=T + 1)
    k = mt_plate + recv
    i = np.arange(T + 1)
    # area-rule prediction (rate per unit area h(t) = hits_t / sum of areas so far)
    cum = np.cumsum(a)
    h = np.zeros(T + 1); h[2:] = mt_plate[2:] / cum[1:-1]
    H = np.cumsum(h[::-1])[::-1]; Hn = np.zeros(T + 1); Hn[:-1] = H[1:]
    pa = a * Hn
    sel = (i >= 1000) & (i < T // 100)
    cb = np.floor(np.log(np.maximum(i, 1)) / np.log(1.25)).astype(int)
    # (a)
    r = recv[sel].sum() / pa[sel].sum()
    idx = np.where(sel)[0]
    lr = np.log(recv[idx] + 0.5); la = np.log(a[idx])
    for arr in (lr, la):
        mu = np.bincount(cb[idx], weights=arr) / np.maximum(np.bincount(cb[idx]), 1)
        arr -= mu[cb[idx]]
    sl = np.polyfit(la, lr, 1)[0]
    out(f'[{label}] (a) i in [1e3,T/100): sum recv/sum pred_area = {r:.4f}; corr(recv, pred) = {np.corrcoef(recv[sel], pa[sel])[0, 1]:.4f}; '
        f'within-cohort slope of log recv on log a = {sl:.3f} (fitness x age: 1; PA: 0)')
    # (b) dispersion around the cohort mean and around the area-x-cohort prediction
    def cohort_mean(y, s):
        mu = np.bincount(cb[s], weights=y[s], minlength=cb.max() + 1) / np.maximum(np.bincount(cb[s], minlength=cb.max() + 1), 1)
        return mu[cb]
    pc = cohort_mean(recv.astype(float), sel)
    ratio = cohort_mean(recv.astype(float), sel) / np.maximum(cohort_mean(pa, sel), 1e-12)
    pac = pa * ratio
    rows = []
    for lo, hi in ((10, 30), (30, 100), (100, 300), (300, 1000)):
        for nm, pr in (('cohort', pc), ('area', pac)):
            b = sel & (pr >= lo) & (pr < hi)
            if b.sum() > 30:
                rows.append(f'{nm}[{lo},{hi}) {np.mean((recv[b] - pr[b])**2) / np.mean(pr[b]):.2f}')
    out(f'[{label}] (b) dispersion <(recv-pred)^2>/<pred> by predicted mean: ' + '; '.join(rows))
    # (c) reinforcement at fixed area
    s2 = (i >= 1000) & (i < 10000)
    idx = np.where(s2)[0]
    inwin = (Eh >= 1000) & (Eh < 10000)
    e_h = Eh[inwin]; e_t = Et[inwin]
    early = np.bincount(e_h[e_t <= 2 * e_h], minlength=T + 1)[idx].astype(float)
    late = np.bincount(e_h[e_t > 2 * e_h], minlength=T + 1)[idx].astype(float)
    v = [np.log(late + 0.5), np.log(early + 0.5), np.log(a[idx])]
    c2 = np.floor(np.log(idx) / np.log(1.1)).astype(int)
    for arr in v:
        mu = np.bincount(c2, weights=arr) / np.maximum(np.bincount(c2), 1)
        arr -= mu[c2]
    def resid(y, x):
        return y - np.polyfit(x, y, 1)[0] * x
    pl = np.corrcoef(resid(v[0], v[2]), resid(v[1], v[2]))[0, 1]
    co = np.linalg.lstsq(np.c_[v[1], v[2]], v[0], rcond=None)[0]
    out(f'[{label}] (c) i in [1e3,1e4): corr(late, early) within cohort = {np.corrcoef(v[0], v[1])[0, 1]:.4f}; '
        f'partial corr at fixed area = {pl:.4f}; log late = {co[0]:.3f} log early + {co[1]:.3f} log a')
    # (d) growth exponent of fixed plates and degree tail
    ts = np.unique(np.logspace(4, np.log10(T), 16).astype(int))
    mk = (Eh >= 1000) & (Eh < 10000); hh = Eh[mk]; tt = Et[mk]
    ls = [np.mean(np.log(2 + np.bincount(hh[tt <= tc], minlength=10000)[1000:10000])) for tc in ts]
    loc = np.diff(ls) / np.diff(np.log(ts))
    out(f'[{label}] (d) fixed plates i in [1e3,1e4): local d<ln s>/d ln t, last 4 intervals mean {loc[-4:].mean():.4f}')
    kk = k[1:]
    out(f'[{label}] (d) MLE k>=20,40,80,160,320,640: ' + ' '.join(f'{gamma_mle(kk, q)[0]:.3f}' for q in (20, 40, 80, 160, 320, 640))
        + ' | K^1.5 P(k>=K) K=20..640: ' + ' '.join(f'{K**1.5 * np.mean(kk >= K):.2f}' for K in (20, 40, 80, 160, 320, 640))
        + ' | K^(2+2/m-1) P: ' + ' '.join(f'{K**1.44 * np.mean(kk >= K):.2f}' for K in (20, 40, 80, 160, 320, 640)))
    return k


if __name__ == '__main__':
    name = sys.argv[1]
    lines = []

    def out(s):
        print(s, flush=True); lines.append(s)

    S, E, C = load(name)
    P = plate_arrays(S, E); T = P['T']
    a = P['area'].copy()
    mtp = P['out'].astype(np.int64)
    out(f'== {name}: T={T}; FIT and PA surrogates use the SR plate out-degrees m_t (mean {mtp[1:].mean():.4f}) and SR birth areas')
    pe = E['h'] >= 0
    kS = analyse('SR ', E['t'][pe].astype(np.int64), E['h'][pe].astype(np.int64), a, mtp, T, out)
    ht, hh = grow_fit(mtp, a, T, 11)
    kF = analyse('FIT', ht, hh, a, mtp, T, out)
    ht, hh = grow_pa(mtp, T, 12)
    kP = analyse('PA ', ht, hh, a, mtp, T, out)
    # full-shape comparison of the degree CCDF over k in [10, 2000]
    K = np.unique(np.logspace(1, np.log10(2000), 25).astype(int))
    cc = lambda k: np.array([np.mean(k[1:] >= q) for q in K])
    cS, cF, cP = cc(kS), cc(kF), cc(kP)
    out('[shape] CCDF ratio SR/FIT at K=' + ','.join(str(q) for q in K[::4]) + ': '
        + ' '.join(f'{x:.3f}' for x in (cS / cF)[::4]))
    out('[shape] CCDF ratio SR/PA  at same K: ' + ' '.join(f'{x:.3f}' for x in (cS / cP)[::4]))
    out(f'[shape] max |ln(SR/FIT)| = {np.max(np.abs(np.log(cS / cF))):.3f}; max |ln(SR/PA)| = {np.max(np.abs(np.log(cS / cP))):.3f}')
    open(os.path.join(RES, f'mechanism_{name}.txt'), 'w').write('\n'.join(lines) + '\n')
