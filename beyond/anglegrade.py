"""Cut-angle grading at fixed density: rule (area, centroid, ('perpmix', kappa)).
kappa = 0: balanced; kappa = 1: centroid. Framed periodic n, eta, clamped omega1 (rho 0.1, 0.2, 0.4,
4 seeds); Bloch band edge along x and y at rho = 0.2. -> results/anglegrade.json (resumable)"""
import json, os, time, numpy as np, divisions, periodic, waves
OUT = 'results/anglegrade.json'
out = json.load(open(OUT)) if os.path.exists(OUT) else []
done = {(o['kappa'], o['seed'], o['rho_t']) for o in out}
t0 = time.time()
for kappa in (0.0, 0.2, 0.4, 0.6, 0.8, 1.0):
    rule = ('area', 'centroid', ('perpmix', kappa))
    for s in range(4):
        g = divisions.run(1000, 100 + s, rule)
        L = np.linalg.norm(g['B'] - g['A'], axis=1)
        for r in (0.1, 0.2, 0.4):
            if (kappa, s, r) in done:
                continue
            lam = np.full(1000, r / L.sum())
            sc = periodic.build(g, lam, frame_lam=lam[0])
            C = periodic.periodic_C(sc)
            th, n = waves.christoffel_n(C, sc['rho'])
            K = (C[0, 0] + C[1, 1] + 2 * C[0, 1]) / 4; G = ((C[0, 0] + C[1, 1] - 2 * C[0, 1]) / 4 + C[2, 2]) / 2
            w1 = waves.analyse(g, lam)['omega1']
            o = dict(kappa=kappa, seed=s, rho_t=r, rho=sc['rho'], n_x=float(n[0]), n_y=float(n[36]),
                     n_mean=float(n.mean()), n_iso=float(np.sqrt(sc['rho'] / (K + G))), K=float(K), G=float(G),
                     eta_n=float(n.min() / n.max()), omega1=w1)
            if r == 0.2:
                for dname, ang, idx in (('x', 0.0, 0), ('y', np.pi / 2, 36)):
                    n0 = float(n[idx])
                    b = periodic.bloch_x(sc, n0, kmax=min(2 * np.pi, 3 * n0 * w1), nk=24, nev=10, direction=ang)
                    k, w, p = b[:, 0], b[:, 1], b[:, 2]
                    bad = (p < 0.3) | (np.abs(k / w / n0 - 1) > 0.1) | ~np.isfinite(w)
                    j = int(np.argmax(bad)) if bad.any() else len(k)
                    o['omega_h_' + dname] = float(w[j - 1]) if j > 0 else 0.0
                    o['n0_' + dname] = n0
            out.append(o); json.dump(out, open(OUT, 'w'))
        print(kappa, s, round(time.time() - t0), flush=True)
