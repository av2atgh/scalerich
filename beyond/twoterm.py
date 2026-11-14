"""Two-term test.  C11 = rho [s + b (lam/l)^2]  ->  n^2 = 1/(s + b (lam/l)^2),  n_max = s^-1/2.
s from the pin-jointed periodic stiffness (bending switched off, tiny regulariser), per cell.
Then omega(n) = Omega sqrt(n^-2 - n_max^-2): one fitted Omega per rule, checked on omega1 (psweep)
and on the Bloch edge.  -> results/twoterm.json"""
import json, numpy as np, divisions, periodic, waves
p = json.load(open('results/psweep_T1000.json'))
out = []
for rule in ('SR', 'centroid', 'balanced', 'regular', 'uniform'):
    for s in range(4):
        g = divisions.run(1000, 100 + s, rule)
        L = np.linalg.norm(g['B'] - g['A'], axis=1)
        lam = np.full(1000, 0.2 / L.sum())
        sc = periodic.build(g, lam, frame_lam=lam[0], bend=1e-9)
        C = periodic.periodic_C(sc)
        th, n = waves.christoffel_n(C, sc['rho'])
        s_x = C[0, 0] / sc['rho']
        out.append(dict(rule=rule, seed=s, s_x=float(s_x), nmax_x=float(n[0]), nmax_mean=float(n.mean())))
        print(rule, s, 's_x %.4f  n_max_x %.2f' % (s_x, n[0]), flush=True)
json.dump(out, open('results/twoterm.json', 'w'), indent=1)
