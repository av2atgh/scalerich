"""Referee round 2, item A. For every cell used in the collapse (T = 300, 1e3, 3000; rules and the
cut-angle family; seeds 100-103): the direction-averaged pinned index of the framed tiling
(s_pin = <n_pin>^-2, independent of rho at alpha = 0), the pinned C11/lambda (number of straight
spanning lines), the same with a frame 1000x thinner, and the full index with that thin frame.
-> results/pinned_all.json"""
import json, os, sys, time
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); B = os.path.join(HERE, '..')
sys.path[:0] = [B, os.path.join(B, '..', 'rg'), HERE]
import divisions, periodic, waves, netstats, srtree
OUT = os.path.join(B, 'results', 'pinned_all.json')
out = json.load(open(OUT)) if os.path.exists(OUT) else []
done = {(o['rule'], o['T'], o['seed']) for o in out}
names = ['SR', 'centroid', 'balanced', 'regular', 'rect'] + ['kappa%g' % k for k in (0, 0.2, 0.4, 0.6, 0.8, 1)]
for T in (1000, 300, 3000):
    for name in names:
        if T != 1000 and name.startswith('kappa'): continue
        for s in range(4):
            if (name, T, s) in done: continue
            t0 = time.time()
            g = divisions.run(T, 100 + s, netstats.rule_of(name)); L = np.linalg.norm(g['B'] - g['A'], axis=1)
            nodes, seg, sl, bn = srtree.network(g, frame=True)
            ls = np.linalg.norm(nodes[seg[:, 1]] - nodes[seg[:, 0]], axis=1); ls = ls[ls > 1e-12]
            o = dict(rule=name, T=T, seed=s, Ltot=float(L.sum()), l_star=float(((ls ** 11).sum() / ls.sum()) ** 0.1))
            lam = np.full(T, 0.2 / L.sum())
            for tag, fl in (('', 1.0), ('_thin', 1e-3)):
                sc = periodic.build(g, lam, frame_lam=lam[0] * fl, bend=1e-9); C = periodic.periodic_C(sc)
                th, n = waves.christoffel_n(C, sc['rho'])
                o['npin' + tag] = float(n.mean()); o['m' + tag] = float(0.5 * (C[0, 0] + C[1, 1]) / lam[0])
            o['n_thin'] = {}
            for rho in (0.1, 0.2, 0.4):
                lam = np.full(T, rho / L.sum())
                sc = periodic.build(g, lam, frame_lam=lam[0] * 1e-3); C = periodic.periodic_C(sc)
                th, n = waves.christoffel_n(C, sc['rho']); o['n_thin'][str(rho)] = float(n.mean())
            out.append(o); json.dump(out, open(OUT, 'w'))
            print(name, T, s, 'npin %.2f m %.2f | thin: npin %.1f m %.3f | n_thin(.2) %.2f  %.0fs' % (
                o['npin'], o['m'], o['npin_thin'], o['m_thin'], o['n_thin']['0.2'], time.time() - t0), flush=True)
