"""Weakest-span prediction for the dense Bloch cells (framed tilings, alpha = 0):
omega_span = min over spans of 4.730^2 lam / (sqrt(12) l^2), spans longer than lam,
frame spans included. -> results/spanmode_bloch2.json"""
import json, os, sys, numpy as np, divisions
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'rg'))
import srtree
b = json.load(open('results/bloch2_T1000.json'))
cache = {}; out = []
for o in b:
    key = (o['rule'], o['seed'])
    if key not in cache:
        g = divisions.run(1000, 100 + o['seed'], o['rule'])
        nodes, seg, sl, bn = srtree.network(g, frame=True)
        Ls = np.linalg.norm(nodes[seg[:, 1]] - nodes[seg[:, 0]], axis=1)
        Lg = np.linalg.norm(g['B'] - g['A'], axis=1)
        cache[key] = (Ls, Lg.sum())
    Ls, Ltot = cache[key]
    lam = o['rho_t'] / Ltot
    ok = Ls > lam
    w = 4.730 ** 2 * lam / (np.sqrt(12) * Ls[ok] ** 2)
    out.append(dict(rule=o['rule'], seed=o['seed'], rho_t=o['rho_t'], dir=o['dir'], omega_span=float(w.min()),
                    omega_h=o['omega_h'], n0=o['n0'], resolved=o['resolved']))
json.dump(out, open('results/spanmode_bloch2.json', 'w'))
v = [x for x in out if x['resolved'] and x['omega_h'] > 0]
x = np.log([q['omega_span'] for q in v]); y = np.log([q['omega_h'] for q in v])
print(len(v), 'corr', round(np.corrcoef(x, y)[0, 1], 3), 'median ratio', round(float(np.exp(np.median(y - x))), 2))
