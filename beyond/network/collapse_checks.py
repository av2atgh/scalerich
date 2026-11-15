"""Referee 1: leave-one-rule-out prediction, residuals against density, rms against p,
fixed per-rule pin ceilings.  -> results/collapse_checks.json"""
import json, os, sys
import numpy as np
from scipy.optimize import least_squares
HERE = os.path.dirname(os.path.abspath(__file__)); B = os.path.join(HERE, '..')
sys.path[:0] = [HERE]
import collapse as C

FR = ['SR', 'centroid', 'balanced', 'regular', 'rect'] + ['k%g' % k for k in (0, 0.2, 0.4, 0.6, 0.8, 1)]
lp = lambda ls, p: ((ls ** (p + 1)).sum() / ls.sum()) ** (1 / p)
D = [d for d in C.data if d[0] in FR]
lab = np.array([d[0] for d in D]); rho = np.array([d[1] for d in D]); y = np.log([d[2] for d in D])
X = np.log(np.array([d[4] / lp(d[3], 10) for d in D]))

def model(th, x, s=None):
    s_ = np.exp(th[0]) if s is None else s
    return -0.5 * np.log(s_ + np.exp(th[1] + th[2] * x))

def fit(m, s=None):
    th0 = [np.log(0.02), 0.0, 1.0] if s is None else [0.0, 0.0, 1.0]
    r = least_squares(lambda th: model(th, X[m], None if s is None else s[m]) - y[m], th0)
    return r.x

out = {}
# leave one rule out (the kappa family counts as one group)
grp = np.array(['kappa' if l.startswith('k') else l for l in lab])
loo = {}
for g in ('SR', 'centroid', 'balanced', 'regular', 'rect', 'kappa'):
    th = fit(grp != g); res = model(th, X[grp == g]) - y[grp == g]
    loo[g] = dict(rms=float(np.sqrt(np.mean(res ** 2))), bias=float(res.mean()), q=float(th[2]), s=float(np.exp(th[0])))
    print('leave out %-8s predicted rms %.3f bias %+.3f (fit q %.2f, n_max %.2f)' % (g, loo[g]['rms'], loo[g]['bias'], th[2], np.exp(-th[0] / 2)))
out['loo'] = loo
th = fit(np.ones(len(y), bool)); res = model(th, X) - y
out['all'] = dict(rms=float(np.sqrt(np.mean(res ** 2))))
# residuals against density, per rule
byrho = {}
for g in ('SR', 'centroid', 'balanced', 'regular', 'rect', 'kappa'):
    byrho[g] = {str(r): float(res[(grp == g) & (rho == r)].mean()) for r in (0.1, 0.2, 0.4)}
    print('residual vs rho %-8s' % g, ' '.join('%s:%+.3f' % kv for kv in byrho[g].items()))
out['res_rho'] = byrho
# without rect
m = grp != 'rect'; th2 = fit(m); r2 = model(th2, X[m]) - y[m]
out['no_rect'] = dict(rms=float(np.sqrt(np.mean(r2 ** 2))), q=float(th2[2]), nmax=float(np.exp(-th2[0] / 2)))
print('without rect: rms %.3f q %.2f n_max %.2f' % (out['no_rect']['rms'], th2[2], np.exp(-th2[0] / 2)))
# rms against p
rp = {}
for p in (1, 2, 4, 6, 8, 10, 15, 20, 40):
    Xp = np.log(np.array([d[4] / lp(d[3], p) for d in D]))
    rr = least_squares(lambda th: model(th, Xp) - y, [np.log(0.02), 0, 1]).fun
    rp[p] = float(np.sqrt(np.mean(rr ** 2)))
Xm = np.log(np.array([d[4] / d[3].max() for d in D]))
rp['max'] = float(np.sqrt(np.mean(least_squares(lambda th: model(th, Xm) - y, [np.log(0.02), 0, 1]).fun ** 2)))
Xc = np.log(np.array([d[4] for d in D]))
rp['none (lambda only)'] = float(np.sqrt(np.mean(least_squares(lambda th: model(th, Xc) - y, [np.log(0.02), 0, 1]).fun ** 2)))
print('rms vs p', rp)
out['rms_p'] = rp
# per-rule pin ceilings fixed (direction-averaged pinned index from twoterm.json, rule mean)
tt = json.load(open(os.path.join(B, 'results', 'twoterm.json')))
pin = {r: np.mean([1 / o['nmax_mean'] ** 2 for o in tt if o['rule'] == r]) for r in ('SR', 'centroid', 'balanced', 'regular')}
pin['kappa'] = None
m = np.isin(grp, ['SR', 'centroid', 'balanced', 'regular'])
s_fix = np.array([pin[g] if g in pin and pin[g] is not None else np.nan for g in grp])
r3 = least_squares(lambda th: model(th, X[m], s_fix[m]) - y[m], [0, 0, 1])
r4 = least_squares(lambda th: model(th, X[m]) - y[m], [np.log(0.02), 0, 1])
out['pin_fixed'] = dict(rms=float(np.sqrt(np.mean(r3.fun ** 2))), q=float(r3.x[2]),
                        rms_common=float(np.sqrt(np.mean(r4.fun ** 2))), q_common=float(r4.x[2]),
                        nmax_pin={k: float(v ** -0.5) for k, v in pin.items() if v})
print('per-rule pinned ceilings fixed (4 rules): rms %.3f q %.2f | common fitted s: rms %.3f q %.2f' % (
    out['pin_fixed']['rms'], r3.x[2], out['pin_fixed']['rms_common'], r4.x[2]), out['pin_fixed']['nmax_pin'])
json.dump(out, open(os.path.join(B, 'results', 'collapse_checks.json'), 'w'), indent=1)
