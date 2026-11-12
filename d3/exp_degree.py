"""Goal 3: degree theory for 3D plates.

Two candidate rate equations for the hits received by plate i:
  (I) incidence PA: rate proportional to the number s_i of faces of i
      (s_i = 2 + hits received)  -> b = m/(2+m), gamma = 2 + 2/m;
  (A) area rule: rate proportional to the plate area a_i (both sides), i.e. a
      uniform hit rate per unit plate area h(t) = m_t / S(t),
      S = 2 sum_j a_j + 6 (cube)  -> k_i ~ (T/i)^(2/3), gamma = 5/2.
Tests: fixed-plate growth exponent, cohort exponent, MLE, rate per incidence and
per area by age (snapshot), regression on area and perimeter, host-age law,
tail amplitude.

usage: python3 exp_degree.py run1 [run2 ...] -> results/degree_<run>.txt
"""
import sys, os
import numpy as np
from common import load, plate_arrays, gamma_mle, cohort_beta, RES

for name in sys.argv[1:]:
    lines = []

    def out(s):
        print(s, flush=True); lines.append(s)

    S, E, C = load(name)
    P = plate_arrays(S, E)
    T = P['T']; k = P['k']; inter = P['interior']; a = P['area']; per = P['perim']
    recv = P['recv']; outd = P['out']
    g = S[:, 0] > 0
    mall = S[g, 3].mean(); mp = (S[g, 3] - S[g, 4]).mean()
    out(f'== {name}: T={T}; <m> all faces {mall:.4f}, plate faces {mp:.4f}; <k>={k[1:].mean():.4f} '
        f'(2<m_p> = {2 * mp:.4f})')
    out(f'   predictions: incidence PA b=m/(2+m)={mp / (2 + mp):.4f}, gamma=2+2/m={2 + 2 / mp:.4f}; '
        f'area rule b=2/3, gamma=5/2')
    # --- MLE. Interior-only selection biases the tail (large old plates touch the cube), so report both.
    for lab, sel in (('all plates', P['ok']), ('interior plates', inter)):
        kk = k[sel]
        s = ' '.join(f'k>={q}:{gamma_mle(kk, q)[0]:.3f}(n={gamma_mle(kk, q)[1]})' for q in (20, 40, 80, 160, 320, 640, 1280))
        out(f'[MLE] {lab}: {s}; kmax={kk.max()}')
    for lo in (10, 100, 1000, 10000):
        out(f'[interior fraction] i in [{lo},{10 * lo}): {inter[lo:10 * lo].mean():.3f}; <k> interior '
            f'{k[lo:10 * lo][inter[lo:10 * lo]].mean():.1f} vs touching cube {k[lo:10 * lo][~inter[lo:10 * lo] & P["ok"][lo:10 * lo]].mean():.1f}')
    # cohort exponent
    for lo, hi in ((1e2, T / 1e2), (1e3, T / 1e2), (10, T / 1e3)):
        out(f'[cohort] beta over i in [{lo:.0e},{hi:.0e}]: total k {cohort_beta(k, inter, lo, hi):.4f}, '
            f'received {cohort_beta(recv.astype(float) + 1e-9, inter, lo, hi):.4f}  (2/3=0.6667, inc-PA {mp / (2 + mp):.4f})')
    # --- fixed-plate growth: s_i(t) for old plates
    ts = np.unique(np.logspace(4, np.log10(T), 16).astype(int))
    Et = E['t']; Eh = E['h']
    for lo, hi in ((10, 100), (100, 1000), (1000, 10000)):
        sel = np.arange(lo, hi); sel = sel[inter[lo:hi]]
        mask_h = (Eh >= lo) & (Eh < hi)
        hh = Eh[mask_h]; tt = Et[mask_h]
        ls = []
        for tc in ts:
            cnt = np.bincount(hh[tt <= tc], minlength=hi)[sel]
            ls.append(np.mean(np.log(2 + cnt)))
        sl = np.polyfit(np.log(ts), ls, 1)[0]
        loc = np.diff(ls) / np.diff(np.log(ts))
        out(f'[fixed plates] i in [{lo},{hi}) (n={len(sel)}): d<ln s_i>/d ln t over [1e4,T] = {sl:.4f}; '
            f'local slopes over the last 4 intervals {np.round(loc[-4:], 4).tolist()} mean {loc[-4:].mean():.4f}  '
            f'(area rule 2/3, incidence PA {mp / (2 + mp):.4f})')
    # --- area rule prediction for each plate
    mt = np.zeros(T + 1); mt[S[g, 0].astype(int)] = S[g, 3]
    Scum = 6 + 2 * np.cumsum(a)          # S(t) after step t
    h = np.zeros(T + 1); h[1:] = mt[1:] / Scum[:-1]   # hit rate per unit side area at step t
    H = np.cumsum(h[::-1])[::-1]          # H[i] = sum_{t>=i} h(t)
    Hn = np.zeros(T + 1); Hn[:-1] = H[1:]   # sum over t > i
    pred = 2 * a * Hn
    for lo, hi in ((1e2, 1e3), (1e3, 1e4), (1e4, 1e5), (1e5, 1e6), (1e6, T)):
        sel = np.zeros(T + 1, bool); sel[int(lo):int(hi)] = True; sel &= inter
        if sel.sum() > 10:
            out(f'[area rule] i in [{lo:.0e},{hi:.0e}): sum recv / sum pred = {recv[sel].sum() / pred[sel].sum():.4f}, '
                f'corr(recv,pred) = {np.corrcoef(recv[sel], pred[sel])[0, 1]:.3f}, <recv>={recv[sel].mean():.2f}')
    sel = inter.copy(); sel[:100] = False; sel[T // 10:] = False
    out(f'[area rule] all interior i in [100, T/10): sum recv/sum pred = {recv[sel].sum() / pred[sel].sum():.4f}, '
        f'corr = {np.corrcoef(recv[sel], pred[sel])[0, 1]:.4f}, corr(log) = '
        f'{np.corrcoef(np.log(recv[sel] + 1), np.log(pred[sel] + 1))[0, 1]:.4f}')
    # cohort-normalised regression: log(recv/<recv>_c) on log(a/<a>_c), log(p/<p>_c)
    sel2 = sel.copy(); sel2[T // 100:] = False   # recv >~ 30: log not compressed by Poisson noise
    idx = np.where(sel2)[0]
    cb = np.floor(np.log(idx) / np.log(1.25)).astype(int)
    lr = np.log(recv[idx] + 0.5); la = np.log(a[idx]); lp = np.log(per[idx])
    for arr in (lr, la, lp):
        mu = np.bincount(cb, weights=arr) / np.maximum(np.bincount(cb), 1)
        arr -= mu[cb]
    X = np.c_[la, lp]
    coef = np.linalg.lstsq(X, lr, rcond=None)[0]
    out(f'[regression] i in [100,T/100), within-cohort log recv ~ {coef[0]:.3f} log a + {coef[1]:.3f} log p '
        f'(area rule: 1, 0; p ~ a^(1/2) collinear: corr(la,lp)={np.corrcoef(la, lp)[0, 1]:.3f}); '
        f'single-variable slopes: a {np.polyfit(la, lr, 1)[0]:.3f}, p {np.polyfit(lp, lr, 1)[0]:.3f}')
    # --- snapshot test: rate per incidence and per area by age h/t0
    for t0 in (T // 4, T // 2):
        t1 = int(t0 * 1.1)
        before = Et <= t0
        win = (Et > t0) & (Et <= t1) & (Eh >= 0)
        s0 = 2 + np.bincount(Eh[before & (Eh >= 0)], minlength=T + 1)[:t0 + 1]
        Hw = np.bincount(Eh[win], minlength=T + 1)[:t0 + 1]
        aa = a[:t0 + 1]
        ok = inter[:t0 + 1].copy(); ok[0] = False
        ages = np.arange(t0 + 1) / t0
        bins = [1e-5, 1e-4, 1e-3, 1e-2, 0.03, 0.1, 0.3, 0.6, 0.9, 1.0001]
        r_inc_all = Hw[ok].sum() / s0[ok].sum(); r_area_all = Hw[ok].sum() / aa[ok].sum()
        row_i, row_a = [], []
        for lo, hi in zip(bins[:-1], bins[1:]):
            b = ok & (ages >= lo) & (ages < hi)
            row_i.append(Hw[b].sum() / s0[b].sum() / r_inc_all)
            row_a.append(Hw[b].sum() / aa[b].sum() / r_area_all)
        out(f'[snapshot t0={t0}, window 10%] age bins h/t0 edges {bins}')
        out('    hit rate per incidence / mean: ' + ' '.join(f'{x:.3f}' for x in row_i))
        out('    hit rate per unit area / mean: ' + ' '.join(f'{x:.3f}' for x in row_a))
        # perimeter-of-face test: incidences per unit area (faces get smaller with age)
        dens = [s0[ok & (ages >= lo) & (ages < hi)].sum() / aa[ok & (ages >= lo) & (ages < hi)].sum() / (s0[ok].sum() / aa[ok].sum())
                for lo, hi in zip(bins[:-1], bins[1:])]
        out('    faces per unit area / mean:     ' + ' '.join(f'{x:.3f}' for x in dens))
    # --- host-age law
    sel = (Et > T // 10) & (Eh > 0)
    x = Eh[sel] / Et[sel]
    out('[host age] P(h/t <= y) for y = 1e-3, 1e-2, 0.1, 0.3, 0.5, 0.8: '
        + ' '.join(f'{np.mean(x <= y):.4f}/{y**(1/3):.4f}' for y in (1e-3, 1e-2, 0.1, 0.3, 0.5, 0.8))
        + '  (measured / U^3 law)')
    out('[host age] <(h/t)^al> for al = 1/6, 1/3, 1/2, 1: '
        + ' '.join(f'{np.mean(x**al):.4f}/{1 / (1 + 3 * al):.4f}' for al in (1/6, 1/3, 1/2, 1))
        + '  (measured / 1/(1+3 al))')
    out(f'[host age] incidence-PA alternative P(h/t<=y)=y^(1-b) with b={mp / (2 + mp):.4f}: at y=0.01 {0.01**(1 - mp / (2 + mp)):.4f}')
    # --- tail amplitude
    i = np.arange(T + 1)
    ca = np.mean((a * i ** (2 / 3))[T // 10:])
    xx = a[1:] * i[1:] ** (2 / 3) / ca
    allp = P['ok'][1:]
    x15 = np.mean(xx[allp] ** 1.5)
    amp = (mall / 2) ** 1.5 * x15
    ka = k[P['ok']]
    out(f'[tail] c_a={ca:.4f}, <x^1.5>={x15:.4f} (x = a_i i^(2/3)/c_a, all plates); predicted K^1.5 P(k>=K) -> (m/2)^1.5<x^1.5> = {amp:.3f}; '
        f'measured (all plates) at K=20,40,80,160,320,640,1280: ' + ' '.join(f'{K**1.5 * np.mean(ka >= K):.3f}' for K in (20, 40, 80, 160, 320, 640, 1280)))
    open(os.path.join(RES, f'degree_{name}.txt'), 'w').write('\n'.join(lines) + '\n')
