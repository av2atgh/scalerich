"""Within-rule test: across 16 realisations of each rule (T = 1e3, alpha = 0, rho = 0.2), does the
longest span explain the realisation-to-realisation spread of the index and of omega1?
-> results/seedvar.json"""
import json, os, sys
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); B = os.path.join(HERE, '..')
sys.path[:0] = [B, os.path.join(B, '..', 'rg'), HERE]
import divisions, srtree, periodic, waves, netstats
OUT = os.path.join(B, 'results', 'seedvar.json')
out = json.load(open(OUT)) if os.path.exists(OUT) else []
done = {(o['rule'], o['seed']) for o in out}
for name in ('SR', 'centroid', 'balanced', 'regular', 'rect', 'kappa0.4'):
    for s in range(200, 216):
        if (name, s) in done: continue
        g = divisions.run(1000, s, netstats.rule_of(name))
        L = np.linalg.norm(g['B'] - g['A'], axis=1); lam = np.full(1000, 0.2 / L.sum())
        sc = periodic.build(g, lam, frame_lam=lam[0])
        C = periodic.periodic_C(sc); th, n = waves.christoffel_n(C, sc['rho'])
        try: w1 = waves.analyse(g, lam)['omega1']
        except Exception: w1 = float('nan')
        nodes, seg, sl, bn = srtree.network(g, frame=True)
        ls = np.linalg.norm(nodes[seg[:, 1]] - nodes[seg[:, 0]], axis=1); ls = ls[ls > 1e-12]
        cell_ext = []
        for p in g['leaf_polys']:
            d = p[:, None, :] - p[None, :, :]; cell_ext.append(np.sqrt((d ** 2).sum(-1)).max())
        o = dict(rule=name, seed=s, n_mean=float(n.mean()), n_x=float(n[0]), eta=float(n.min() / n.max()),
                 omega1=float(w1), lam=float(lam[0]), span_max=float(ls.max()),
                 span_p10=float(((ls ** 11).sum() / ls.sum()) ** 0.1), span_mean=float(ls.mean()),
                 cell_diam_max=float(max(cell_ext)), area_max=float(g['leaf_area'].max()))
        out.append(o); json.dump(out, open(OUT, 'w'))
        print(name, s, 'n %.2f w1 %.3f lmax %.3f lp10 %.3f dmax %.3f' % (o['n_mean'], w1, o['span_max'], o['span_p10'], o['cell_diam_max']), flush=True)
