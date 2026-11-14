"""Bloch dispersion of the longitudinal branch along x on framed periodic tilings (2D).
For each cell: n0 from the periodic C, then k in (0, 2pi] (wavelength down to one cell).
Homogenisation edge omega_h: first k at which the plane-wave projection p_L < 0.3 or
|n(omega)/n0 - 1| > 0.1; omega_h is the last omega before it.  Compared with omega1.
-> results/bloch_T1000.json"""
import json, time, numpy as np, divisions, periodic, waves

RULES = ('SR', 'balanced', 'regular', 'centroid', 'uniform')
out = []; t0 = time.time()
for rule in RULES:
    for s in (0, 1):
        g = divisions.run(1000, 100 + s, rule)
        L = np.linalg.norm(g['B'] - g['A'], axis=1)
        for a in (0.0, 0.5):
            base = np.arange(1, 1001.) ** -a
            for r in (0.05, 0.1, 0.2, 0.4):
                lam = base * r / np.sum(base * L)
                sc = periodic.build(g, lam, frame_lam=lam[0])
                C = periodic.periodic_C(sc)
                th, n = waves.christoffel_n(C, sc['rho'])
                n0 = float(n[0])
                b = periodic.bloch_x(sc, n0, nk=48, nev=12)
                k, w, p = b[:, 0], b[:, 1], b[:, 2]
                nk = k / w
                bad = (p < 0.3) | (np.abs(nk / n0 - 1) > 0.1) | ~np.isfinite(w)
                j = int(np.argmax(bad)) if bad.any() else len(k)
                wh = float(w[j - 1]) if j > 0 else 0.0
                kh = float(k[j - 1]) if j > 0 else 0.0
                out.append(dict(rule=rule, seed=s, alpha=a, rho_t=r, rho=sc['rho'], n0=n0,
                                n_mean=float(n.mean()), eta_n=float(n.min() / n.max()),
                                omega_h=wh, k_h=kh, k=k.tolist(), omega=w.tolist(), pL=p.tolist()))
                json.dump(out, open('results/bloch_T1000.json', 'w'))
        print(rule, s, round(time.time() - t0), flush=True)
