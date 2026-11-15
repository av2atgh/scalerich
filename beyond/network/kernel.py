"""Referee 8: measured attachment kernel kappa(k) = t * (hit rate per node of degree k) for each rule,
T = 2e4, seeds 1-2, steps t > 2000.  Linear fit kappa = a + b k over 2 <= k <= 20.
-> results/kernel.json"""
import json, os, sys
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); B = os.path.join(HERE, '..')
sys.path[:0] = [B, HERE]
import divisions, netstats
KM = 60
out = {}
for name in ('SR', 'uniform', 'rect', 'centroid', 'balanced', 'regular'):
    hits = np.zeros(KM); expo = np.zeros(KM)
    for sd in (1, 2):
        g = divisions.run(20000, sd, netstats.rule_of(name)); T = g['T']
        deg = np.zeros(T, int); cnt = np.zeros(KM)     # cnt[k] = number of nodes with degree k
        for t in range(T):
            if t > 2000:
                expo += cnt / t                          # rate per node = kappa(k)/t
            for h in (g['hostA'][t], g['hostB'][t]):
                if h >= 0:
                    k = min(deg[h], KM - 1)
                    if t > 2000: hits[k] += 1
                    cnt[k] -= 1; deg[h] += 1; cnt[min(deg[h], KM - 1)] += 1
            d0 = int(g['hostA'][t] >= 0) + int(g['hostB'][t] >= 0)
            deg[t] = d0; cnt[d0] += 1
    kap = hits / np.maximum(expo, 1e-30)
    k = np.arange(KM); m = (k >= 2) & (k <= 20) & (hits > 50)
    b, a = np.polyfit(k[m], kap[m], 1)
    out[name] = dict(k=k[m].tolist(), kappa=kap[m].tolist(), a=float(a), b=float(b))
    print(name, 'kappa = %.3f + %.3f k' % (a, b), ' '.join('%d:%.2f' % (kk, v) for kk, v in zip(k[m][:10], kap[m][:10])), flush=True)
json.dump(out, open(os.path.join(B, 'results', 'kernel.json'), 'w'), indent=1)
