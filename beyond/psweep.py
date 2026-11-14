"""Periodic (mirror-supercell) effective stiffness for every 2D cell of sweep_T1000 + lattices.
-> results/psweep_T1000.json (framed periodic tiling) (same keys as sweep + n_x_kin, n_mean_kin)"""
import json, time, numpy as np, divisions, waves, periodic, lattices
import frame
rows = json.load(open('results/sweep_T1000.json'))
out = []; cache = {}; t0 = time.time()
for r in rows:
    rule = r['rule']
    if rule in divisions.RULES:
        key = (rule, r['seed'])
        if key not in cache:
            cache.clear(); cache[key] = divisions.run(1000, 100 + r['seed'], rule)
        g = cache[key]
        L = np.linalg.norm(g['B'] - g['A'], axis=1)
        base = np.arange(1, 1001.) ** -r['alpha']
        lam = base * r['rho_t'] / np.sum(base * L)
        sc = periodic.build(g, lam, frame_lam=lam[0])
    else:
        fn = getattr(lattices, rule)
        P, seg, bnd = fn()
        Ls = np.linalg.norm(P[seg[:, 1]] - P[seg[:, 0]], axis=1)
        lam_s = np.full(len(seg), r['rho_t'] / Ls.sum())
        try:
            X, S, Lm = periodic.mirror_supercell(P, seg, lam_s)
            K = frame.assemble(X, S, Lm)[0]
            master, shift = periodic.pairing(X)
            M, Lx = periodic.mass_vector(X, S, Lm)
            sc = dict(X=X, S=S, lam=Lm, K=K, M=M, master=master, shift=shift, rho=float((Lm * Lx).sum() / 4))
        except AssertionError as e:
            print('skip', rule, e); continue
    C = periodic.periodic_C(sc)
    th, n = waves.christoffel_n(C, sc['rho'])
    o = dict(rule=rule, seed=r['seed'], alpha=r['alpha'], rho_t=r['rho_t'], rho=sc['rho'],
             n_x=float(n[0]), n_mean=float(n.mean()), n_min=float(n.min()), n_max=float(n.max()),
             eta_n=float(n.min() / n.max()), C=C.tolist(), omega1=r['omega1'],
             n_x_kin=r['n_x'], n_mean_kin=r['n_mean'], eta_kin=r['eta_n'])
    out.append(o)
json.dump(out, open('results/psweep_T1000.json', 'w'))
print('done', len(out), round(time.time() - t0))
