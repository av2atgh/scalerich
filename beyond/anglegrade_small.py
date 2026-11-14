"""Cell tables for the cut-angle family with T ligaments per cell (for the full-lattice check).
Same format as results/anglegrade.json; rho 0.1, 0.2, 0.4; 4 seeds; no Bloch.
-> results/anglegrade_T{T}.json"""
import json, sys, numpy as np, divisions, periodic, waves
T = int(sys.argv[1]) if len(sys.argv) > 1 else 100
out = []
for kappa in (0.0, 0.2, 0.4, 0.6, 0.8, 1.0):
    for s in range(4):
        g = divisions.run(T, 700 + s, ('area', 'centroid', ('perpmix', kappa)))
        L = np.linalg.norm(g['B'] - g['A'], axis=1)
        for r in (0.1, 0.2, 0.4):
            lam = np.full(T, r / L.sum())
            sc = periodic.build(g, lam, frame_lam=lam[0])
            C = periodic.periodic_C(sc)
            K = (C[0, 0] + C[1, 1] + 2 * C[0, 1]) / 4; G = ((C[0, 0] + C[1, 1] - 2 * C[0, 1]) / 4 + C[2, 2]) / 2
            w1 = waves.analyse(g, lam)['omega1']
            out.append(dict(kappa=kappa, seed=s, rho_t=r, rho=sc['rho'], K=float(K), G=float(G),
                            n_iso=float(np.sqrt(sc['rho'] / (K + G))), omega1=w1,
                            omega_h_x=0.8 * w1, omega_h_y=0.8 * w1))
json.dump(out, open(f'results/anglegrade_T{T}.json', 'w'))
for k in (0.0, 0.5, 1.0):
    v = [o for o in out if abs(o['kappa'] - k) < 0.11 and o['rho_t'] == 0.2]
    if v: print('kappa ~%.1f rho .2: n %.2f +- %.2f  omega1 %.3f' % (k, np.mean([o['n_iso'] for o in v]), np.std([o['n_iso'] for o in v]), np.mean([o['omega1'] for o in v])))
