"""Referee round 2, item A.3(ii): collapse with s fixed at each cell's own pinned value s_pin(T),
only c and q fitted; per-rule offsets, leave-one-out, and out-of-sample T = 300, 3000.
Also the thin-frame tiling (frame 1000x thinner): bending term without the frame lines.
-> results/collapse_pin.json"""
import json, os, sys
import numpy as np
from scipy.optimize import least_squares
HERE = os.path.dirname(os.path.abspath(__file__)); B = os.path.join(HERE, '..')
P = json.load(open(os.path.join(B, 'results', 'pinned_all.json')))
ps = json.load(open(os.path.join(B, 'results', 'psweep_T1000.json')))
ag = json.load(open(os.path.join(B, 'results', 'anglegrade.json')))
ts = json.load(open(os.path.join(B, 'results', 'tsize.json')))
pin = {(o['rule'], o['T'], o['seed']): o for o in P}
rows = []   # (group, T, rho, n, lam/l*, s_pin, n_thin)
for x in ps:
    if x['alpha'] == 0 and x['rule'] in ('SR', 'centroid', 'balanced', 'regular', 'rect') and x['rho_t'] in (0.1, 0.2, 0.4):
        o = pin[(x['rule'], 1000, x['seed'])]
        rows.append((x['rule'], 1000, x['rho_t'], x['n_mean'], x['rho_t'] / o['Ltot'] / o['l_star'], o['npin'] ** -2, o['n_thin'][str(x['rho_t'])]))
for x in ag:
    if x['rho_t'] in (0.1, 0.2, 0.4):
        o = pin[('kappa%g' % x['kappa'], 1000, x['seed'])]
        rows.append(('kappa', 1000, x['rho_t'], x['n_mean'], x['rho_t'] / o['Ltot'] / o['l_star'], o['npin'] ** -2, o['n_thin'][str(x['rho_t'])]))
for x in ts:
    o = pin[(x['rule'], x['T'], x['seed'])]
    rows.append((x['rule'], x['T'], x['rho_t'], x['n_mean'], x['rho_t'] / o['Ltot'] / o['l_star'], o['npin'] ** -2, o['n_thin'][str(x['rho_t'])]))
g = np.array([r[0] for r in rows]); T = np.array([r[1] for r in rows]); rho = np.array([r[2] for r in rows])
n = np.array([r[3] for r in rows]); X = np.log([r[4] for r in rows]); sp = np.array([r[5] for r in rows]); nt = np.array([r[6] for r in rows])
y = np.log(n)
def model(th, m, s=None, free_s=False):
    s_ = np.exp(th[2]) if free_s else s[m]
    return -0.5 * np.log(s_ + np.exp(th[0] + th[1] * X[m]))
def fit(m, s=sp, free_s=False, yy=None):
    yy = y if yy is None else yy
    r = least_squares(lambda th: model(th, m, s, free_s) - yy[m], [0, 1, np.log(0.02)] if free_s else [0, 1])
    return r.x, r.fun
out = {}
train = (T == 1000) & (g != 'rect')
for label, free in (('s = s_pin(cell)', False), ('common s', True)):
    th, f = fit(train, free_s=free)
    res = model(th, np.ones(len(y), bool), sp, free) - y
    d = dict(c=float(np.exp(th[0])), q=float(th[1]), rms=float(np.sqrt(np.mean(f ** 2))))
    if free: d['nmax'] = float(np.exp(-th[2] / 2))
    d['offset'] = {r: float(res[train & (g == r)].mean()) for r in ('SR', 'centroid', 'balanced', 'regular', 'kappa')}
    for TT in (300, 3000):
        d['T%d' % TT] = {r: [float(np.sqrt(np.mean(res[(T == TT) & (g == r)] ** 2))), float(res[(T == TT) & (g == r)].mean())]
                         for r in ('SR', 'centroid', 'balanced', 'regular')}
    out[label] = d
    print(label, 'c %.2f q %.2f rms %.3f' % (d['c'], d['q'], d['rms']), d.get('nmax', ''), 'offsets', {k: round(v, 3) for k, v in d['offset'].items()})
    for TT in (300, 3000): print('   out-of-sample T=%d rms/bias' % TT, {k: tuple(round(z, 3) for z in v) for k, v in d['T%d' % TT].items()})
# leave one out with s_pin
loo = {}
for r in ('SR', 'centroid', 'balanced', 'regular', 'kappa'):
    th, _ = fit(train & (g != r)); res = model(th, train & (g == r), sp) - y[train & (g == r)]
    loo[r] = float(np.sqrt(np.mean(res ** 2)))
out['loo_spin'] = loo; print('LOO (s_pin):', {k: round(v, 3) for k, v in loo.items()})
# thin frame: bending term alone, s from thin pinned (tiny)
yt = np.log(nt)
th, f = fit(train, s=np.zeros(len(y)), yy=yt)
res = model(th, np.ones(len(y), bool), np.zeros(len(y))) - yt
out['thin'] = dict(c=float(np.exp(th[0])), q=float(th[1]), rms=float(np.sqrt(np.mean(f ** 2))),
                   offset={r: float(res[train & (g == r)].mean()) for r in ('SR', 'centroid', 'balanced', 'regular', 'kappa')},
                   T3000={r: float(np.sqrt(np.mean(res[(T == 3000) & (g == r)] ** 2))) for r in ('SR', 'centroid', 'balanced', 'regular')},
                   T300={r: float(np.sqrt(np.mean(res[(T == 300) & (g == r)] ** 2))) for r in ('SR', 'centroid', 'balanced', 'regular')})
print('thin frame, pure power law n^-2 = c (lam/l*)^q: c %.2f q %.2f rms %.3f' % (out['thin']['c'], out['thin']['q'], out['thin']['rms']),
      {k: round(v, 3) for k, v in out['thin']['offset'].items()}, 'T300', {k: round(v, 3) for k, v in out['thin']['T300'].items()}, 'T3000', {k: round(v, 3) for k, v in out['thin']['T3000'].items()})
# m per rule and T
mm = {}
for r in ('SR', 'centroid', 'balanced', 'regular', 'rect'):
    mm[r] = {TT: [float(np.mean([o['m'] for o in P if o['rule'] == r and o['T'] == TT])), float(np.mean([o['m_thin'] for o in P if o['rule'] == r and o['T'] == TT])),
                  float(np.mean([o['npin'] for o in P if o['rule'] == r and o['T'] == TT]))] for TT in (300, 1000, 3000)}
    print(r, 'm (framed, thin) and <n_pin> by T:', {k: tuple(round(z, 2) for z in v) for k, v in mm[r].items()})
out['m'] = mm
json.dump(out, open(os.path.join(B, 'results', 'collapse_pin.json'), 'w'), indent=1)
# points for Fig. 4(b): bending term n^-2 - s_pin against lam/l*, T = 1e3
th, _ = fit(train)
out['fig'] = dict(c=float(np.exp(th[0])), q=float(th[1]),
                  pts=[dict(rule=str(g[i]), T=int(T[i]), x=float(np.exp(X[i])), b=float(n[i] ** -2 - sp[i])) for i in range(len(y))])
json.dump(out, open(os.path.join(B, 'results', 'collapse_pin.json'), 'w'), indent=1)
bad = [p for p in out['fig']['pts'] if p['b'] <= 0]
print('points with n above own pinned ceiling:', len(bad), set(p['rule'] for p in bad))
