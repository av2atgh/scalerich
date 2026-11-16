"""Referee controls: (i) area pick with fixed r (golden, aligned): N_s, I_max, m, chi;
(ii) 2D peak search for plastic (3,2), (5,1) and anti-aligned golden at T = 16383 (brightest peak
on an axis?); (iii) diffuse baseline: I_max of a Poisson point set with the same N.
-> results/controls.json"""
import json, numpy as np, qdiv
from diffraction import vertices, bragg
from pisot import axis_top
from oned import theta
from mech import lines
from netq import xfrac
out = {}
TAU = qdiv.TAU
# (i)
rows = []
for T in (1000, 2000, 4000):
    for seed in (0, 1):
        g = qdiv.run(T, ('area', 1 / TAU ** 2, 'perp'), seed=seed)
        wh = np.sort(g['leaf_wh'], axis=1) * np.sqrt(T + 1)
        sim = len(set(np.round(wh[:, 1] / wh[:, 0], 6)))                    # shapes up to similarity
        I, k = axis_top(vertices(g), T); ln = lines(g)
        A = g['leaf_area'] * (T + 1)
        rows.append(dict(T=T, seed=seed, Ns_cong=len(set(map(tuple, np.round(wh, 6)))), Ns_sim=sim,
                         area_ratio=float(A.max() / A.min()), Imax=I, m=int((ln > 1 - 1e-6).sum()), chi=xfrac(g)))
        print('area-pick golden', rows[-1], flush=True)
out['area_pick'] = rows
# (ii)
rows = []
for nm, sp in [('(3,2)', theta(3, 2)[0] ** -3), ('(5,1)', theta(3, 2)[0] ** -5), ('golden-avoid', ('avoid', 1 / TAU ** 2, 'low'))]:
    T = 16383; P = vertices(qdiv.run(T, ('largest', sp, 'perp'))); I, q = bragg(P, T, ncand=60)
    rows.append(dict(rule=nm, T=T, I2d=I, q=q)); print('2D', nm, round(I, 4), np.round(q, 3), flush=True)
out['peak2d'] = rows
# (iii)
rng = np.random.default_rng(0); rows = []
for T in (4095, 65535, 1048575):
    N = int(1.2 * T); P = rng.random((N, 2)); I, k = axis_top(P, T)
    rows.append(dict(T=T, N=N, Imax=I)); print('Poisson', T, N, I, flush=True)
out['poisson'] = rows
# (iv) shapes up to similarity for the aligned commensurate splits
rows = []
for nm, r in [('(1,1)', 0.5), ('(2,1)', TAU ** -2), ('(3,1)', theta(3, 1)[0] ** -3), ('(3,2)', theta(3, 2)[0] ** -3),
              ('(5,1)', theta(3, 2)[0] ** -5), ('(4,1)', theta(4, 1)[0] ** -4), ('(4,3)', theta(4, 3)[0] ** -4), ('r=0.3', 0.3)]:
    for T in (1000, 100000):
        g = qdiv.run(T, ('largest', r, 'perp')); wh = np.sort(g['leaf_wh'], axis=1)
        rows.append(dict(rule=nm, T=T, Ns_sim=len(set(np.round(wh[:, 1] / wh[:, 0], 6))),
                         Ns_cong=len(set(map(tuple, np.round(wh * np.sqrt(T + 1), 6))))))
        print('shapes', rows[-1], flush=True)
out['shapes_sim'] = rows
json.dump(out, open('results/controls.json', 'w'), indent=1)
