"""Statistics of the SR splitting tree (Goal 1).

- split-fraction distribution r (area fraction of a child), size-biased
  mu = E[-ln r'] and s2 = Var[-ln r'] where r' is the child keeping a tagged point
- depth of the cell split at step t vs ln t (renewal prediction: mean ln t/mu,
  variance (s2/mu^3) ln t)
- moments of leaf area and perimeter, M_q(t) = sum_leaves X^q ~ t^-gamma(q):
  ordinary scaling gamma_A(q) = q-1, gamma_P(q) = q/2-1; control: the
  Krapivsky-Ben-Naim rectangular model (mode='rect'), where the exact moment
  equation dM(a,b)/dt = [1/(a+1)+1/(b+1)-1] M(a+1,b+1) has the conserved family
  ab = 1 and gives gamma(a,0) = (a - sqrt(a^2+4))/2  (multiscaling of lengths)
- corner-number distribution, cell shape anisotropy stationarity
- thickness ratio chord/host: x = ln((t+1)/(h+1)) for host ligament h of chord t
  (stationary in t => lam_t/lam_host = e^{-alpha x} is a marginal, t-independent
  coupling)
- segment length of ligament i at final T vs birth index i
Writes results/tree_stats.json.
"""
import json, sys
import numpy as np
import srtree

T = int(sys.argv[1]) if len(sys.argv) > 1 else 30000
NS = int(sys.argv[2]) if len(sys.argv) > 2 else 12
checkpoints = np.unique(np.round(np.logspace(2, np.log10(T), 9)).astype(int))
qs = [0.5, 1.5, 2.0, 3.0]
aq = [0.5, 1.0, 2.0, 3.0, 4.0]
out = {'T': T, 'nseed': NS, 'checkpoints': checkpoints.tolist()}


def perim(poly):
    return float(np.linalg.norm(np.roll(poly, -1, 0) - poly, axis=1).sum())


def shape_ratio(poly):
    c = poly.mean(0); X = poly - c
    ev = np.linalg.eigvalsh(X.T @ X)
    return ev[0] / ev[1]


def ensemble(mode):
    res = {}
    R = []; D = []; ncorn = []; X = []; segL = []; shp = {}
    MA = np.zeros((NS, len(checkpoints), len(qs)))
    MP = np.zeros((NS, len(checkpoints), len(aq)))
    MX = np.zeros((NS, len(checkpoints), len(aq)))  # rect: sum x^a (width)
    for s in range(NS):
        g = srtree.run(T, 100 + s, mode)
        R.append(g['r']); D.append(g['depth'])
        ncorn.append(g['cell_ncorner'][g['leaf']])
        ncell = len(g['cell_area'])
        split_t = np.full(ncell, T + 1); split_t[g['parent']] = np.arange(T)
        born = g['cell_born']
        P = np.array([perim(p) for p in g['polys']])
        W = np.array([np.ptp(p[:, 0]) for p in g['polys']])
        for ci, tau in enumerate(checkpoints):
            alive = (born <= tau - 1) & (split_t >= tau)
            a = g['cell_area'][alive]
            MA[s, ci] = [np.sum(a ** q) for q in qs]
            MP[s, ci] = [np.sum(P[alive] ** q) for q in aq]
            MX[s, ci] = [np.sum(W[alive] ** q) for q in aq]
            if mode == 'iso' and s < 4:
                idx = np.where(alive)[0]
                shp.setdefault(int(tau), []).extend(
                    [shape_ratio(g['polys'][k]) for k in idx[:3000]])
        h = np.concatenate([g['hostA'], g['hostB']])
        tt = np.concatenate([np.arange(T)] * 2)
        m = h >= 0
        X.append(np.c_[tt[m], np.log((tt[m] + 1.0) / (h[m] + 1.0))])
        nodes, seg, sl, _ = srtree.network(g)
        Ls = np.linalg.norm(nodes[seg[:, 0]] - nodes[seg[:, 1]], axis=1)
        segL.append(np.c_[sl, Ls])
    R = np.concatenate(R); D = np.array(D); X = np.concatenate(X)
    segL = np.concatenate(segL)
    # split fractions: tagged point keeps child 0 with prob r
    rr = np.concatenate([R, 1 - R]); w = rr / rr.sum()
    mu = float(np.sum(w * -np.log(rr))); s2 = float(np.sum(w * np.log(rr) ** 2) - mu ** 2)
    hist, edges = np.histogram(np.minimum(R, 1 - R), bins=10, range=(0, 0.5), density=True)
    res['split'] = dict(mean_min_r=float(np.minimum(R, 1 - R).mean()), mu=mu, s2=s2,
                        hist_min_r=hist.round(3).tolist())
    # depth vs ln t
    tb = np.unique(np.round(np.logspace(1, np.log10(T - 1), 12)).astype(int))
    dm = [float(D[:, max(1, b // 2):b].mean()) for b in tb]
    dv = [float(D[:, max(1, b // 2):b].var()) for b in tb]
    lt = np.log([np.sqrt(b * max(1, b // 2)) for b in tb])
    sl_m = np.polyfit(lt[3:], dm[3:], 1); sl_v = np.polyfit(lt[3:], dv[3:], 1)
    res['depth'] = dict(t=tb.tolist(), mean=np.round(dm, 3).tolist(), var=np.round(dv, 3).tolist(),
                        slope_mean=float(sl_m[0]), pred_slope_mean=1 / mu,
                        slope_var=float(sl_v[0]), pred_slope_var=s2 / mu ** 3)
    # moments
    lc = np.log(checkpoints)
    mA = np.log(MA.mean(0)); mP = np.log(MP.mean(0)); mX = np.log(MX.mean(0))
    fit = lambda M, k: -float(np.polyfit(lc[3:], M[3:, k], 1)[0])
    res['gamma_area'] = {str(q): fit(mA, k) for k, q in enumerate(qs)}
    res['gamma_area_scaling'] = {str(q): q - 1 for q in qs}
    res['gamma_perim'] = {str(a): fit(mP, k) for k, a in enumerate(aq)}
    res['gamma_perim_scaling'] = {str(a): a / 2 - 1 for a in aq}
    res['gamma_width'] = {str(a): fit(mX, k) for k, a in enumerate(aq)}
    res['gamma_width_KBN'] = {str(a): (a - np.sqrt(a * a + 4)) / 2 for a in aq}
    cn = np.concatenate(ncorn)
    res['corners'] = {int(k): float(np.mean(cn == k)) for k in range(3, 10)}
    res['corners_mean'] = float(cn.mean())
    if shp:
        res['shape_ratio_mean_vs_t'] = {k: float(np.mean(v)) for k, v in shp.items()}
    # chord/host log-age gap, binned by t
    xb = {}
    for lo in [10, 30, 100, 300, 1000, 3000, 10000]:
        if lo * 3 > T: break
        m = (X[:, 0] >= lo) & (X[:, 0] < 3 * lo)
        x = X[m, 1]
        xb[lo] = dict(mean=float(x.mean()), sd=float(x.std()),
                      E_ratio_a0_25=float(np.mean(np.exp(-0.25 * x))),
                      E_ratio_a0_5=float(np.mean(np.exp(-0.5 * x))),
                      E_ratio_a1=float(np.mean(np.exp(-1.0 * x))))
    res['host_gap'] = xb
    # segment length vs birth index: mean over segments of ligament i, times sqrt(T)
    sb = {}
    for lo in [1, 10, 100, 1000, 10000]:
        if lo >= T: break
        m = (segL[:, 0] >= lo - 1) & (segL[:, 0] < min(T, 3 * lo) - 1 + (lo == 1))
        sb[lo] = float(np.mean(segL[m, 1]) * np.sqrt(T))
    res['segment_len_sqrtT_by_birth'] = sb
    return res


for mode in ('iso', 'rect'):
    out[mode] = ensemble(mode)
    print(mode, json.dumps(out[mode], indent=1)[:4000])
json.dump(out, open('results/tree_stats.json', 'w'), indent=1)
