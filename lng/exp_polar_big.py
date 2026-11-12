"""Asymptotics of gamma under orientation fields: T up to 1e6 (wave w=0.125, circ, uniform control).
Output: results/polar_big.csv"""
import sys, os, math, time
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import srcells, graphstats as gs
from exp_topology import edges_from_hits
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'results')
circ = lambda x, y: math.atan2(y - 0.5, x - 0.5) + math.pi / 2
wave = lambda x, y: math.pi * (x + y) / 0.125
f = open(os.path.join(OUT, 'polar_big.csv'), 'w')
f.write('rule,T,seed,gmle10,gmle20,gmle40,beta,secs\n')
for name, kw in [('uniform', {}), ('wave_w0.125', dict(angle='field', field=wave)), ('circ', dict(angle='field', field=circ))]:
    for T in (300000, 1000000):
        t0 = time.time(); o = srcells.run(T, seed=7, **kw)
        A = gs.adjacency(edges_from_hits(o['hit']), T); k = gs.degrees(A)
        ii = np.unique(np.logspace(np.log10(max(10, T / 3000)), np.log10(T / 10), 12).astype(int))
        mk = [k[a:int(a * 1.3) + 1].mean() for a in ii]
        beta = -np.polyfit(np.log(ii), np.log(mk), 1)[0]
        f.write(f"{name},{T},7,{gs.gamma_mle(k,10)[0]:.3f},{gs.gamma_mle(k,20)[0]:.3f},{gs.gamma_mle(k,40)[0]:.3f},{beta:.3f},{time.time()-t0:.0f}\n"); f.flush()
f.close()
