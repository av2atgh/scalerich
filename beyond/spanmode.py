"""Parameter-free prediction of omega1: lowest clamped-clamped beam mode over all spans,
omega_s = 4.730^2 * sqrt(EI/(rho A)) / l^2 = 22.37 * lam / (sqrt(12) l^2).  2D, T = 1000."""
import json, sys, os, numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'rg'))
import divisions, srtree
rows = json.load(open('results/sweep_T1000.json'))
out = []
cache = {}
for r in rows:
    if r['rule'] not in divisions.RULES:
        continue
    key = (r['rule'], r['seed'])
    if key not in cache:
        g = divisions.run(1000, 100 + r['seed'], r['rule'])
        nodes, seg, sl, bn = srtree.network(g)
        L = np.linalg.norm(nodes[seg[:, 1]] - nodes[seg[:, 0]], axis=1)
        Lg = np.linalg.norm(g['B'] - g['A'], axis=1)
        cache[key] = (L, sl, Lg)
    L, sl, Lg = cache[key]
    base = np.arange(1, 1001.) ** -r['alpha']
    lam = base * r['rho_t'] / np.sum(base * Lg)
    ls = lam[sl]
    ok = L > ls                      # spans longer than their thickness
    w = 4.730 ** 2 * ls[ok] / (np.sqrt(12) * L[ok] ** 2)
    out.append((r['rule'], r['alpha'], r['rho_t'], r['omega1'], float(w.min())))
json.dump(out, open('results/spanmode.json', 'w'))
a = np.array([(o[3], o[4]) for o in out])
lr = np.log(a[:, 0] / a[:, 1])
print('all cells: omega1/omega_span  median %.2f  IQR %.2f-%.2f  corr(log) %.3f' % (
    np.exp(np.median(lr)), np.exp(np.percentile(lr, 25)), np.exp(np.percentile(lr, 75)),
    np.corrcoef(np.log(a[:, 0]), np.log(a[:, 1]))[0, 1]))
for rule in divisions.RULES:
    v = np.array([(o[3], o[4]) for o in out if o[0] == rule])
    q = np.log(v[:, 0] / v[:, 1])
    print(f'{rule:9s} ratio median {np.exp(np.median(q)):.2f}  corr {np.corrcoef(np.log(v[:,0]), np.log(v[:,1]))[0,1]:.3f}')
