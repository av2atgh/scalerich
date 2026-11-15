"""Network characterisation of every recursive-division rule (zero thickness).

Mappings of one division g (divisions.run):
  cell tree    binary tree of cells; split fraction a0/(a0+a1), depth of every leaf, leaf areas
  cells        corner number (vertices of the leaf polygon), elongation sqrt(I_max/I_min)
  dual graph   node = ligament, edge = T-junction (stem -- host); boundary sides excluded
  frame graph  node = junction, edge = span between consecutive junctions on a ligament
  geometry     ligament lengths and orientations; junction angle between stem and host
Usage: python netstats.py asym   (T = 2e4, 4 seeds, all rules + cut-angle family)  -> results/net_asym.json
       python netstats.py phys   (T = 1e3, seeds 100..103 = the cells of psweep)   -> results/net_phys.json
"""
import json, os, sys, time
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); B = os.path.join(HERE, '..')
sys.path[:0] = [B, os.path.join(B, '..', 'rg'), os.path.join(B, '..', 'lng')]
import divisions, srtree, graphstats

KAPPAS = (0.0, 0.2, 0.4, 0.6, 0.8, 1.0)


def rule_of(name):
    if name.startswith('kappa'):
        return ('area', 'centroid', ('perpmix', float(name[5:])))
    return name


def features(g, paths=True, nsrc=40):
    T = g['T']; f = {}
    # --- cell tree
    fr = g['frac']; fr = np.minimum(fr, 1 - fr)          # smaller daughter fraction in (0, 1/2]
    f['frac_min_mean'] = float(fr.mean()); f['frac_min_q05'] = float(np.quantile(fr, 0.05))
    f['frac_hist'] = np.histogram(fr, bins=20, range=(0, 0.5))[0].tolist()
    d = g['leaf_depth']; f['depth_mean'] = float(d.mean()); f['depth_var'] = float(d.var())
    f['depth_max'] = int(d.max())
    A = g['leaf_area'] * (T + 1)                         # mean 1
    f['area_cv'] = float(A.std()); f['area_max'] = float(A.max()); f['area_min'] = float(A.min())
    f['lnA_var'] = float(np.var(np.log(np.maximum(A, 1e-300))))
    f['sumA2'] = float((g['leaf_area'] ** 2).sum())
    # --- cells
    nc = np.array([len(p) for p in g['leaf_polys']])
    f['corner_frac'] = {int(k): float((nc == k).mean()) for k in range(3, 9)}
    f['corner_mean'] = float(nc.mean())
    el = []
    for p in g['leaf_polys']:
        a, _, I2 = divisions._moments(p)
        w = np.linalg.eigvalsh(I2)
        el.append(np.sqrt(w[1] / max(w[0], 1e-300)) if a > 0 else np.nan)
    el = np.array(el); el = el[np.isfinite(el)]
    f['elong_mean'] = float(el.mean()); f['elong_q95'] = float(np.quantile(el, 0.95))
    # --- ligament geometry
    v = g['B'] - g['A']; L = np.linalg.norm(v, axis=1); th = np.arctan2(v[:, 1], v[:, 0])
    t1 = np.arange(1, T + 1)
    f['lsqrt_t'] = float(np.mean((L * np.sqrt(t1))[T // 10:]))
    f['L_total_sqrtT'] = float(L.sum() / np.sqrt(T))
    f['fabric2'] = float(abs((L * np.exp(2j * th)).sum()) / L.sum())
    f['fabric4'] = float(abs((L * np.exp(4j * th)).sum()) / L.sum())
    f['fabric4_axis'] = float(((L * np.cos(4 * th)).sum()) / L.sum())   # >0: axis-aligned, <0: diagonal
    # junction angle: |sin| between stem t and host h
    s = []
    for hosts in (g['hostA'], g['hostB']):
        m = hosts >= 0
        vh = v[hosts[m]]
        cr = np.abs(v[m, 0] * vh[:, 1] - v[m, 1] * vh[:, 0]) / (L[m] * np.linalg.norm(vh, axis=1))
        s.append(cr)
    s = np.concatenate(s)
    f['junction_sin_mean'] = float(s.mean()); f['junction_oblique_frac'] = float((s < np.sin(np.pi / 6)).mean())
    # host age h/t for t > T/10
    h = []; tt = []
    for hosts in (g['hostA'], g['hostB']):
        m = (hosts >= 0) & (t1 - 1 > T // 10)
        h.append(hosts[m] + 1.0); tt.append(t1[m].astype(float))
    ht = np.concatenate(h) / np.concatenate(tt)
    f['ht_moments'] = [float(np.mean(ht ** a)) for a in (0.25, 0.5, 1.0)]   # U^2 law: 2/3, 1/2, 1/3
    # --- dual graph
    E = []
    for hosts in (g['hostA'], g['hostB']):
        m = hosts >= 0
        E.append(np.stack([np.nonzero(m)[0], hosts[m]], 1))
    E = np.concatenate(E)
    Adj = graphstats.adjacency(E, T)
    res, k, C, kn = graphstats.summary(Adj, nsrc=nsrc, paths=paths)
    for key in ('meank', 'kmax', 'c2', 'Cmean', 'transitivity', 'r', 'meand', 'lnmu', 'mu', 'Ck', 'knnk'):
        if key in res: f[key] = res[key]
    f['k2_over_k'] = float((k ** 2).mean() / k.mean())
    f['ccdf'] = {int(K): float((k >= K).mean()) for K in (3, 4, 6, 8, 12, 16, 24, 32, 48, 64, 96, 128)}
    f['pk'] = np.bincount(k, minlength=40)[:40].tolist()
    # --- frame graph: spans
    nodes, seg, sl, bn = srtree.network(g)
    ls = np.linalg.norm(nodes[seg[:, 1]] - nodes[seg[:, 0]], axis=1)
    ls = ls[ls > 1e-12] * np.sqrt(T)
    f['span_mean'] = float(ls.mean()); f['span_cv'] = float(ls.std() / ls.mean())
    f['span_max'] = float(ls.max()); f['span_q99'] = float(np.quantile(ls, 0.99))
    f['span_w'] = float((ls ** 2).sum() / ls.sum())      # length-weighted mean span
    f['span_hist'] = np.histogram(np.log10(ls), bins=30, range=(-3, 1.5))[0].tolist()
    f['nspan_per_lig'] = float(len(ls) / T)
    return f


def run_set(names, T, seeds, paths=True):
    out = []
    for name in names:
        for sd in seeds:
            t0 = time.time()
            g = divisions.run(T, sd, rule_of(name))
            f = features(g, paths=paths); f.update(rule=name, seed=sd, T=T)
            out.append(f)
            print(name, sd, 'k2/k %.2f C %.3f r %.3f  areaCV %.2f  depth %.1f  span_cv %.2f max %.2f  sin %.2f  %.0fs' % (
                f['k2_over_k'], f['Cmean'], f['r'], f['area_cv'], f['depth_mean'], f['span_cv'], f['span_max'],
                f['junction_sin_mean'], time.time() - t0), flush=True)
    return out


if __name__ == '__main__':
    import warnings; warnings.filterwarnings('ignore')
    mode = sys.argv[1]
    names = list(divisions.RULES) + ['kappa%g' % k for k in KAPPAS]
    if mode == 'asym':
        T = int(sys.argv[2]) if len(sys.argv) > 2 else 20000
        out = run_set(names, T, (1, 2, 3, 4))
        json.dump(out, open(os.path.join(B, 'results', f'net_asym_T{T}.json'), 'w'))
    else:
        out = run_set(names, 1000, (100, 101, 102, 103))
        json.dump(out, open(os.path.join(B, 'results', 'net_phys.json'), 'w'))
