"""Pooled (across alpha) per-ligament energy fit from results/beams_T3000_*.csv:
U = C_alpha k^a lam^b f(theta)  vs  U = C_alpha l^a lam^b f(theta)  vs  C_alpha k^a i^c f.
alpha-specific intercepts absorb lam0 and density; b is identified only by the
variation of alpha.  Output: results/beams_pool.txt"""
import os, glob, math
import numpy as np
from beams import fit
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'results')
rows = []
for j, fn in enumerate(sorted(glob.glob(os.path.join(OUT, 'beams_T3000_a*_s1.csv')))):
    d = np.loadtxt(fn, delimiter=',', skiprows=1)
    i, k, l, lam, th, U = d.T
    m = (U > 0) & (k >= 2) & (i > 1)
    rows.append((j, i[m], k[m], l[m], lam[m], th[m], U[m], os.path.basename(fn)))
y = np.concatenate([np.log(r[6]) for r in rows])
n = len(y); G = len(rows)
Dg = np.zeros((n, G)); o = 0
for r in rows:
    Dg[o:o + len(r[1]), r[0]] = 1; o += len(r[1])
cat = lambda q: np.concatenate([np.log(r[q]) for r in rows])
th = np.concatenate([r[5] for r in rows])
lk, ll, lla, li = cat(2), cat(3), cat(4), cat(1)
out = ['files: ' + ', '.join(r[7] for r in rows)]
def f(names, cols):
    c, R2 = fit(y, [Dg[:, g] for g in range(G)] + cols, [f'g{g}' for g in range(G)] + names, th)
    return ' '.join(f"{a}={c[a]:.3f}" for a in names) + f"  R2={R2:.3f}"
out.append('U ~ k^a lam^b     : ' + f(['k', 'lam'], [lk, lla]))
out.append('U ~ l^a lam^b     : ' + f(['l', 'lam'], [ll, lla]))
out.append('U ~ k^a i^c       : ' + f(['k', 'i'], [lk, li]))
out.append('U ~ k^a l^b lam^c : ' + f(['k', 'l', 'lam'], [lk, ll, lla]))
out.append('U ~ k l lam i     : ' + f(['k', 'l', 'lam', 'i'], [lk, ll, lla, li]))
# paper law imposed: residual slope
comb = 1.5 * lk + 0.5 * lla
out.append('paper comb k^1.5 lam^0.5 (slope, should be 1): ' + f(['comb'], [comb]))
txt = '\n'.join(out); print(txt); open(os.path.join(OUT, 'beams_pool.txt'), 'w').write(txt + '\n')
