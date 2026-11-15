"""Does a span moment collapse the index of every rule onto one law n^-2 = s + c (lam/l_p)^q ?
l_p = (<l^(p+1)> / <l>)^(1/p): length-weighted p-th moment of the span lengths (frame graph).
Data: psweep (alpha = 0, rho_t 0.1-0.4, seeds 100-103) and anglegrade (kappa family)."""
import json, os, sys
import numpy as np
from scipy.optimize import least_squares
HERE = os.path.dirname(os.path.abspath(__file__)); B = os.path.join(HERE, '..')
sys.path[:0] = [B, os.path.join(B, '..', 'rg')]
import divisions, srtree

ps = json.load(open(os.path.join(B, 'results', 'psweep_T1000.json')))
ag = json.load(open(os.path.join(B, 'results', 'anglegrade.json')))
spans = {}
def get(rule, seed):
    key = (str(rule), seed)
    if key not in spans:
        g = divisions.run(1000, 100 + seed, rule)
        nodes, seg, sl, bn = srtree.network(g, frame=True)
        ls = np.linalg.norm(nodes[seg[:, 1]] - nodes[seg[:, 0]], axis=1); ls = ls[ls > 1e-12]
        L = np.linalg.norm(g['B'] - g['A'], axis=1)
        spans[key] = (ls, L.sum())
    return spans[key]

rows = []
for x in ps:
    if x['alpha'] == 0 and x['rule'] in ('SR', 'centroid', 'balanced', 'regular', 'rect', 'uniform'):
        rows.append((x['rule'], x['seed'], x['rho_t'], x['n_mean'], x['rule']))
for x in ag:
    rows.append((('area', 'centroid', ('perpmix', x['kappa'])), x['seed'], x['rho_t'], x['n_mean'], 'k%g' % x['kappa']))
data = []
for rule, seed, rho, n, lab in rows:
    ls, Ltot = get(rule, seed)
    data.append((lab, rho, n, ls, rho / Ltot))

def fit(p, labs_used, per_rule=False):
    D = [d for d in data if d[0] in labs_used]
    lp = np.array([((d[3] ** (p + 1)).sum() / d[3].sum()) ** (1 / p) if p > 0 else np.exp((d[3] * np.log(d[3])).sum() / d[3].sum()) for d in D])
    x = np.log(np.array([d[4] for d in D]) / lp); y = np.array([d[2] for d in D])
    def res(th):
        s, lc, q = th
        return np.log(1 / np.sqrt(np.exp(s) + np.exp(lc + q * x))) - np.log(y)
    r = least_squares(res, [np.log(0.02), 0.0, 1.0])
    return np.sqrt(np.mean(r.fun ** 2)), r.x, D, r.fun

if __name__ == '__main__':
    frame_rules = ['SR', 'centroid', 'balanced', 'regular', 'rect'] + ['k%g' % k for k in (0, 0.2, 0.4, 0.6, 0.8, 1)]
    for p in (0, 1, 2, 3, 4, 6):
        rms, th, D, fun = fit(p, frame_rules)
        lab = np.array([d[0] for d in D])
        per = {r: np.sqrt(np.mean(fun[lab == r] ** 2)) for r in ('SR', 'centroid', 'balanced', 'regular', 'rect')}
        bias = {r: np.mean(fun[lab == r]) for r in ('SR', 'centroid', 'balanced', 'regular', 'rect', 'k0', 'k1')}
        print('p=%d rms ln n %.3f  s %.4f c %.3g q %.2f  bias %s' % (p, rms, np.exp(th[0]), np.exp(th[1]), th[2],
              ' '.join('%s %+.2f' % kv for kv in bias.items())))
    # baseline: per-rule fits with the mean span
    tot = []
    for r in frame_rules:
        rms, th, D, fun = fit(1, [r]); tot += list(fun)
    print('per-rule fits (3 params each): rms %.3f' % np.sqrt(np.mean(np.square(tot))))
    rms, th, D, fun = fit(1, ['uniform']); print('uniform alone rms %.3f q %.2f' % (rms, th[2]))


def dump(p=10):
    # fitted on the rules with a common finite ceiling; rect (no common ceiling) and uniform (frame
    # model invalid) are shown but not fitted
    frame_rules = ['SR', 'centroid', 'balanced', 'regular'] + ['k%g' % k for k in (0, 0.2, 0.4, 0.6, 0.8, 1)]
    rms, th, D, fun = fit(p, frame_rules)
    lp = lambda ls: ((ls ** (p + 1)).sum() / ls.sum()) ** (1 / p)
    pts = [dict(rule=d[0], rho=d[1], n=d[2], x=float(d[4] / lp(d[3]))) for d in data]
    rms0, th0, D0, fun0 = fit(1, frame_rules)
    base = {}
    for d0, f0 in zip(D0, fun0): base.setdefault(d0[0], []).append(f0)
    per = {}
    for d1, f1 in zip(D, fun): per.setdefault(d1[0], []).append(f1)
    json.dump(dict(p=p, rms=rms, s=float(np.exp(th[0])), c=float(np.exp(th[1])), q=float(th[2]), pts=pts,
                   bias={k: float(np.mean(v)) for k, v in per.items()},
                   bias_mean_span={k: float(np.mean(v)) for k, v in base.items()}, rms_mean_span=rms0),
              open(os.path.join(B, 'results', 'collapse.json'), 'w'), indent=1)
    print('p', p, 'rms', rms, 'q', th[2], 's', np.exp(th[0]))
