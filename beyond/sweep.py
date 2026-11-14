"""Sweep: rules x alpha x rho x seeds, plus reference lattices. -> results/sweep.json"""
import json, sys, time
import numpy as np
import divisions, waves, lattices

T = int(sys.argv[1]) if len(sys.argv) > 1 else 1000
NS = 4
ALPHAS = (0.0, 0.25, 0.5, 0.75, 1.0)
RHOS = (0.1, 0.2, 0.3, 0.4)
rows = []
t0 = time.time()
for rule in divisions.RULES:
    for s in range(NS):
        g = divisions.run(T, 100 + s, rule)
        L = np.linalg.norm(g['B'] - g['A'], axis=1)
        for a in ALPHAS:
            base = np.arange(1, T + 1.) ** -a
            for r in RHOS:
                lam = base * r / np.sum(base * L)
                try:
                    o = waves.analyse(g, lam)
                except Exception as e:
                    print('FAILED', rule, s, a, r, e, flush=True); continue
                o.update(rule=rule, seed=s, alpha=a, rho_t=r, lam_max=float(lam[0]))
                rows.append(o)
    print(rule, round(time.time() - t0), flush=True)
for name, fn in (('square', lattices.square), ('triangular', lattices.triangular), ('honeycomb', lattices.honeycomb)):
    P, seg, bnd = fn()
    Ls = np.linalg.norm(P[seg[:, 1]] - P[seg[:, 0]], axis=1)
    for r in RHOS:
        lam_s = np.full(len(seg), r / Ls.sum())
        o = waves.analyse_frame(*waves.frame_system(P, seg, lam_s, bnd))
        o.update(rule=name, seed=0, alpha=0.0, rho_t=r)
        rows.append(o)
    print(name, flush=True)
json.dump(rows, open(f'results/sweep_T{T}.json', 'w'))
