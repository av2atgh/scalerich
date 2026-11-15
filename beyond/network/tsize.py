"""Finite-T test (referee 5): index, isotropy and omega1 at T = 300 and 3000 (T = 1000 is psweep),
alpha = 0, rho_t = 0.1, 0.2, 0.4, seeds 100-103, with the span moment l_* of every cell, so that the
T dependence can be checked against the collapse n^-2 = s + c (lam/l_*)^q.  -> results/tsize.json"""
import json, os, sys, time
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); B = os.path.join(HERE, '..')
sys.path[:0] = [B, os.path.join(B, '..', 'rg')]
import divisions, srtree, periodic, waves
OUT = os.path.join(B, 'results', 'tsize.json')
out = json.load(open(OUT)) if os.path.exists(OUT) else []
done = {(o['rule'], o['T'], o['seed'], o['rho_t']) for o in out}
for T in (300, 3000):
    for rule in ('SR', 'centroid', 'balanced', 'regular', 'rect'):
        for s in range(4):
            g = divisions.run(T, 100 + s, rule)
            nodes, seg, sl, bn = srtree.network(g, frame=True)
            ls = np.linalg.norm(nodes[seg[:, 1]] - nodes[seg[:, 0]], axis=1); ls = ls[ls > 1e-12]
            lstar = float(((ls ** 11).sum() / ls.sum()) ** 0.1)
            L = np.linalg.norm(g['B'] - g['A'], axis=1)
            for r in (0.1, 0.2, 0.4):
                if (rule, T, s, r) in done: continue
                t0 = time.time()
                lam = np.full(T, r / L.sum())
                sc = periodic.build(g, lam, frame_lam=lam[0])
                C = periodic.periodic_C(sc); th, n = waves.christoffel_n(C, sc['rho'])
                try: w1 = float(waves.analyse(g, lam)['omega1'])
                except Exception: w1 = float('nan')
                o = dict(rule=rule, T=T, seed=s, rho_t=r, rho=sc['rho'], n_mean=float(n.mean()), n_x=float(n[0]),
                         eta_n=float(n.min() / n.max()), omega1=w1, lam=float(lam[0]), l_star=lstar,
                         l_max=float(ls.max()))
                out.append(o); json.dump(out, open(OUT, 'w'))
                print(rule, T, s, r, 'n %.2f eta %.2f w1 %.3g  %.0fs' % (o['n_mean'], o['eta_n'], w1, time.time() - t0), flush=True)
