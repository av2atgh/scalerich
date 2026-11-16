"""Mechanics of aligned vs anti-aligned splits: rigid and pinned joints, T = 300 ... 1e4,
frame thickness = ligament thickness ('frame1') or 1e-3 of it ('frame0'). -> results/mech2.json"""
import sys, json, time, numpy as np
sys.path[:0] = ['../beyond', '../rg']
import qdiv, periodic, waves
from mech import build, lines
TAU = qdiv.TAU
RULES = {'golden': ('largest', 1 / TAU ** 2, 'perp'),
         'golden-avoid': ('largest', ('avoid', 1 / TAU ** 2, 'low'), 'perp'),
         'golden-side': ('largest', ('side', 1 / TAU ** 2), 'perp'),
         'uniform': ('largest', 'uniform', 'perp')}
out = []; t0 = time.time()
for T in (300, 1000, 3000, 10000):
    for nm, rule in RULES.items():
        seeds = (0,) if nm in ('golden', 'golden-avoid') else ((0, 1, 2) if T < 10000 else (0,))
        for seed in seeds:
            g = qdiv.run(T, rule, seed=seed); L = np.linalg.norm(g['B'] - g['A'], axis=1)
            m = int((lines(g) > 1 - 1e-6).sum())
            for rho in ((0.1, 0.2, 0.4) if T < 10000 else (0.1,)):
                lam = rho / L.sum()
                for fr in ((1.0, 1e-3) if T <= 3000 else (1.0,)):
                    for bend in (1.0, 0.0):
                        sc = build(g, np.full(T, lam), lam * fr, bend=bend)
                        C = periodic.periodic_C(sc); th, n = waves.christoffel_n(C, sc['rho'])
                        o = dict(rule=nm, T=T, seed=seed, rho_t=rho, frame=fr, bend=bend, m=m,
                                 n_mean=float(n.mean()), n_x=float(n[0]), eta=float(n.min() / n.max()))
                        out.append(o)
                        print(nm, T, seed, rho, fr, bend, 'm=%d n=%.2f eta=%.2f' % (m, o['n_mean'], o['eta']), round(time.time() - t0), flush=True)
            json.dump(out, open('results/mech2.json', 'w'), indent=1)
