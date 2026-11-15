"""Dense Bloch sweep for the index-bandwidth curves (framed periodic tilings, alpha = 0).
Rules SR, centroid, balanced; 4 seeds; 8 densities; directions x and y; k grid up to
k_max = min(2 pi, 3 n0 omega1) with 24 points (resolves the band edge at high n).
Band edge omega_h as in blochsweep.py (p_L < 0.3 or |n/n0 - 1| > 0.1).
-> results/bloch2_T1000.json (resumable)"""
import json, os, time, numpy as np, divisions, periodic, waves
OUT = 'results/bloch2_T1000.json'
out = json.load(open(OUT)) if os.path.exists(OUT) else []
done = {(o['rule'], o['seed'], o['rho_t'], o['dir']) for o in out}
t0 = time.time()
import sys
RULES = sys.argv[1:] or ['SR', 'centroid', 'balanced']
for rule in RULES:
    for s in range(4):
        g = divisions.run(1000, 100 + s, rule)
        L = np.linalg.norm(g['B'] - g['A'], axis=1)
        for r in (0.03, 0.05, 0.07, 0.1, 0.14, 0.2, 0.28, 0.4):
            lam = np.full(1000, r / L.sum())
            sc = periodic.build(g, lam, frame_lam=lam[0])
            C = periodic.periodic_C(sc)
            th, n = waves.christoffel_n(C, sc['rho'])
            w1 = waves.analyse(g, lam)['omega1']
            for dname, ang, idx in (('x', 0.0, 0), ('y', np.pi / 2, 36)):
                if (rule, s, r, dname) in done:
                    continue
                n0 = float(n[idx])
                kmax = min(2 * np.pi, 3 * n0 * w1)
                b = periodic.bloch_x(sc, n0, kmax=kmax, nk=24, nev=10, direction=ang)
                k, w, p = b[:, 0], b[:, 1], b[:, 2]
                bad = (p < 0.3) | (np.abs(k / w / n0 - 1) > 0.1) | ~np.isfinite(w)
                j = int(np.argmax(bad)) if bad.any() else len(k)
                out.append(dict(rule=rule, seed=s, rho_t=r, dir=dname, rho=sc['rho'], n0=n0,
                                n_mean=float(n.mean()), eta_n=float(n.min() / n.max()), omega1=w1,
                                omega_h=float(w[j - 1]) if j > 0 else 0.0, k_h=float(k[j - 1]) if j > 0 else 0.0,
                                resolved=bool(j > 1), saturated=bool(not bad.any()),
                                k=k.tolist(), omega=w.tolist(), pL=p.tolist()))
                json.dump(out, open(OUT, 'w'))
        print(rule, s, round(time.time() - t0), flush=True)
