"""Orientation fields: can a polarized angle rule tune gamma?
Fields: circ (theta perpendicular to the radius from the centre: chord line at
distance |p-c| from the centre), wave (theta = pi*(x+y)), fraction pfield of
ligaments follow the field, rest uniform.  T-scaling for circ, pfield = 1.
Output: results/polar.csv; also degree vs distance of the line from the centre.
"""
import sys, os, time, math
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import srcells, graphstats as gs
from exp_topology import edges_from_hits
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'results')
circ = lambda x, y: math.atan2(y - 0.5, x - 0.5) + math.pi / 2
wave = lambda x, y: math.pi * (x + y)
radial = lambda x, y: math.atan2(y - 0.5, x - 0.5)
f = open(os.path.join(OUT, 'polar.csv'), 'w')
f.write('field,pfield,T,seed,gmle6,gmle10,gmle20,beta,c2,kmax,Cmean,P_tri,corr_k_l_old\n')
runs = [('circ', p, 30000) for p in (0.2, 0.5, 0.8, 1.0)] + [('wave', p, 30000) for p in (0.5, 1.0)] \
     + [('radial', 1.0, 30000)] + [('circ', 1.0, T) for T in (10000, 100000, 200000)]
fields = dict(circ=circ, wave=wave, radial=radial)
for name, pf, T in runs:
    for s in (1, 2):
        try:
            o = srcells.run(T, seed=s, angle='field', field=fields[name], pfield=pf)
            Td = o['T']; h = o['hit']
            A = gs.adjacency(edges_from_hits(h), Td); k = gs.degrees(A)
            C, _ = gs.local_clustering(A, k)
            nb = (h >= 0).all(1)
            pt = (np.asarray(A[h[nb, 0], h[nb, 1]]).ravel() > 0).mean()
            ii = np.unique(np.logspace(np.log10(max(10, Td / 3000)), np.log10(Td / 10), 12).astype(int))
            mk = [k[a:int(a * 1.3) + 1].mean() for a in ii]
            beta = -np.polyfit(np.log(ii), np.log(mk), 1)[0]
            old = np.arange(Td) < Td // 10
            cr = np.corrcoef(np.log(k[old & (k > 2)]), np.log(o['length'][old & (k > 2)]))[0, 1]
            f.write(f"{name},{pf},{Td},{s},{gs.gamma_mle(k,6)[0]:.3f},{gs.gamma_mle(k,10)[0]:.3f},{gs.gamma_mle(k,20)[0]:.3f},"
                    f"{beta:.3f},{(k*(k-1.)).mean()/k.mean():.2f},{k.max()},{C[k>=2].mean():.4f},{pt:.4f},{cr:.3f}\n")
            f.flush()
            if name == 'circ' and pf == 1.0 and T == 100000 and s == 1:
                # where are the hubs? degree vs line distance from centre and vs length
                m = o['mid']; e = o['ends']
                # distance of the line from the centre
                d = np.abs((e[:, 2] - e[:, 0]) * (e[:, 1] - 0.5) - (e[:, 3] - e[:, 1]) * (e[:, 0] - 0.5)) / np.maximum(o['length'], 1e-15)
                top = np.argsort(k)[::-1][:15]
                with open(os.path.join(OUT, 'polar_hubs.txt'), 'w') as g:
                    g.write('rank node k length dist_line_centre birth\n')
                    for r, i in enumerate(top):
                        g.write(f"{r} {i} {k[i]} {o['length'][i]:.4f} {d[i]:.4f} {i}\n")
                    # k / length for old nodes binned by line distance
                    for lo, hi in [(0, .02), (.02, .05), (.05, .1), (.1, .2), (.2, .3), (.3, .5), (.5, .8)]:
                        s_ = old & (d >= lo) & (d < hi) & (o['length'] > 0)
                        if s_.sum():
                            g.write(f"line-dist [{lo},{hi}) n={s_.sum()} mean (k-2)/length={np.mean((k[s_]-2)/o['length'][s_]):.1f}\n")
            print(name, pf, T, s, flush=True)
        except Exception as ex:
            f.write(f"# FAILED {name},{pf},{T},{s}: {ex!r}\n"); f.flush()
f.close()
