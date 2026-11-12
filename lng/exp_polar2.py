"""Mechanism of gamma < 3 under orientation fields: k vs length exponent for old
ligaments, and field wavelength dependence (wave field theta = pi*(x+y)/w).
Output: results/polar2.csv"""
import sys, os, math
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import srcells, graphstats as gs
from exp_topology import edges_from_hits
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'results')
circ = lambda x, y: math.atan2(y - 0.5, x - 0.5) + math.pi / 2
f = open(os.path.join(OUT, 'polar2.csv'), 'w')
f.write('rule,T,seed,gmle10,gmle20,beta,slope_logk_logl_old,mean_misalign_old_vs_new\n')
runs = [('uniform', dict()), ('circ', dict(angle='field', field=circ))]
for w in (1.0, 0.5, 0.25, 0.125):
    runs.append((f'wave_w{w}', dict(angle='field', field=(lambda w: (lambda x, y: math.pi * (x + y) / w))(w))))
for name, kw in runs:
    for T in (30000, 100000):
        s = 1
        o = srcells.run(T, seed=s, **kw)
        h = o['hit']; A = gs.adjacency(edges_from_hits(h), T); k = gs.degrees(A)
        ii = np.unique(np.logspace(np.log10(max(10, T / 3000)), np.log10(T / 10), 12).astype(int))
        mk = [k[a:int(a * 1.3) + 1].mean() for a in ii]
        beta = -np.polyfit(np.log(ii), np.log(mk), 1)[0]
        old = (np.arange(T) < T // 20) & (k > 4)
        sl = np.polyfit(np.log(o['length'][old]), np.log(k[old] - 2), 1)[0]
        # angle between each hit ligament and the newcomer that hit it (|sin|)
        th = o['theta']; nb = h[:, 0] >= 0
        mis = np.mean(np.abs(np.sin(th[nb] - th[h[nb, 0]])))
        f.write(f"{name},{T},{s},{gs.gamma_mle(k,10)[0]:.3f},{gs.gamma_mle(k,20)[0]:.3f},{beta:.3f},{sl:.3f},{mis:.3f}\n"); f.flush()
        print(name, T, flush=True)
f.close()
