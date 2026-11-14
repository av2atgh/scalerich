"""3D plates: framed periodic (symmetry BCs) vs kinematic, n along the axes, T=300
-> ../results/pcheck3d.json"""
import json, time, numpy as np, div3d, mesh3d, shell3d, periodic3d
T = 300; out = []
for rule in ('SR', 'balanced', 'regular', 'centroid', 'rect'):
    for s in (0, 1):
        g = div3d.run(T, 100 + s, rule)
        X, tri, tp = mesh3d.mesh(g, 0.07)
        A = 0.5 * np.linalg.norm(np.cross(X[tri[:, 1]] - X[tri[:, 0]], X[tri[:, 2]] - X[tri[:, 0]]), axis=1)
        g2 = periodic3d.with_frame(g)
        X2, tri2, tp2 = mesh3d.mesh(g2, 0.07)
        for a in (0.0, 0.5):
            lam, _ = shell3d.thickness(g, tp, A, 0.1, a)
            t0 = time.time()
            ok = shell3d.analyse(X, tri, lam[tp], resonance=False)
            Ck = np.array(ok['C'])
            nk = np.sqrt(ok['rho'] / np.diag(Ck)[:3])
            lam2 = np.r_[lam, np.full(6, 0.5 * lam[0])]   # face plate shared by two tiles
            Cp, rho = periodic3d.symmetric_normal(X2, tri2, lam2[tp2])
            npr = np.sqrt(rho / Cp)
            o = dict(rule=rule, seed=s, alpha=a, rho_kin=ok['rho'], rho_per=rho,
                     n_axes_kin=nk.tolist(), n_axes_per=npr.tolist())
            out.append(o); json.dump(out, open('../results/pcheck3d.json', 'w'))
            print(rule, s, a, 'kin', nk.round(3), 'per', npr.round(3), 'rho %.3f/%.3f' % (ok['rho'], rho),
                  '%.0fs' % (time.time() - t0), flush=True)
