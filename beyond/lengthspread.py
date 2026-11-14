"""Index-bandwidth product vs segment-length geometry (2D, alpha = 0)."""
import json, sys, os, numpy as np, collections
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'rg'))
import divisions, srtree
rows = json.load(open('results/sweep_T1000.json'))
geo = collections.defaultdict(list)
for rule in divisions.RULES:
    for s in range(4):
        g = divisions.run(1000, 100 + s, rule)
        nodes, seg, sl, bn = srtree.network(g)
        L = np.linalg.norm(nodes[seg[:, 1]] - nodes[seg[:, 0]], axis=1)
        # mass-weighted mean span (uniform thickness) and the longest span
        lw = (L * L).sum() / L.sum()
        geo[rule].append((lw, L.max(), np.mean(L)))
print('rule      rho  n*w1   <l>_w/lmax^2   ratio')
out = {}
for rule in divisions.RULES:
    lw, lm, l1 = np.mean(geo[rule], axis=0)
    G = lw / lm ** 2
    for r in (0.1, 0.2, 0.4):
        v = [x for x in rows if x['rule'] == rule and x['alpha'] == 0 and x['rho_t'] == r]
        p = np.mean([x['n_mean'] * x['omega1'] for x in v])
        print(f'{rule:9s} {r:4.1f} {p:6.3f} {G:8.2f} {p/G:8.4f}')
    out[rule] = dict(lw=lw, lmax=lm, lmean=l1)
json.dump(out, open('results/lengthspread.json', 'w'), indent=1)
