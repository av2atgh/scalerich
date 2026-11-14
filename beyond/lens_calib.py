"""Lens step 1: tile calibration.  For each rule, T_tile ligaments, alpha = 0, several seeds:
framed periodic n_mean(rho_int) with the frame thickness equal to the tile's ligament thickness
(shared tile edges take the mean in the lens), isotropy, and omega1 of the clamped tile.  Picks the most isotropic seed per rule.
-> results/lens_calib.json"""
import json, sys, numpy as np, divisions, periodic, waves

TT = int(sys.argv[1]) if len(sys.argv) > 1 else 200
RHOS = np.array([0.02, 0.03, 0.05, 0.07, 0.1, 0.14, 0.2, 0.28, 0.4])
out = {}
for rule in ('SR', 'balanced'):
    best = None
    for s in range(8):
        g = divisions.run(TT, 500 + s, rule)
        L = np.linalg.norm(g['B'] - g['A'], axis=1)
        ns, etas, w1s = [], [], []
        for r in RHOS:
            lam = np.full(TT, r / L.sum())
            sc = periodic.build(g, lam, frame_lam=lam[0])   # frame scales with the tile
            C = periodic.periodic_C(sc)
            th, n = waves.christoffel_n(C, sc['rho'])
            ns.append(float(n.mean())); etas.append(float(n.min() / n.max()))
            w1s.append(waves.analyse(g, lam)['omega1'])
        rec = dict(seed=500 + s, rho_int=RHOS.tolist(), n=ns, eta=etas, omega1=w1s,
                   eta_min=float(min(etas)))
        print(rule, s, 'eta_min %.2f n %.2f-%.2f w1(0.1) %.3f' % (rec['eta_min'], min(ns), max(ns), w1s[4]), flush=True)
        if best is None or rec['eta_min'] > best['eta_min']:
            best = rec
    out[rule] = best
json.dump(out, open('results/lens_calib.json', 'w'), indent=1)
