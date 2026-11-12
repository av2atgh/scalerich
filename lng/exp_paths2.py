"""Referee comment 5: multiplicity with error bars and a controlled test.
A) size scaling, 6 seeds, T = 1e3..1e5, 100 sources: SR, SR* (degree-preserving rewiring), DM.
   Also the per-hop growth rate r = d<ln mu>/dd over d in [2, <d>].
B) control at T = 3e4, 6 seeds: excess multiplicity at FIXED degree sequence,
   Delta = <ln mu>(model) - <ln mu>(its rewiring), against the fraction of newcomers
   that close a square (1 - triangle closure), for SR rule variants and surrogates with varied p4.
Output: results/paths2_scaling.csv, results/paths2_control.csv"""
import sys, os, math, time
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import srcells, surrogate, graphstats as gs
from exp_topology import edges_from_hits
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'results')

def tri_frac(h):
    E = edges_from_hits(h); S = set(map(tuple, E.tolist()))
    m = (h[:, 0] >= 0) & (h[:, 1] >= 0)
    a, b = np.minimum(h[m, 0], h[m, 1]), np.maximum(h[m, 0], h[m, 1])
    return np.mean([(x, y) in S for x, y in zip(a.tolist(), b.tolist())])

def pstats(A, s, nsrc=100):
    k = gs.degrees(A); ps = gs.path_stats(A, nsrc=nsrc, seed=s)
    per = np.array([(d, n, lm) for d, n, m, lm in ps['per']])
    sel = (per[:, 0] >= 2) & (per[:, 0] <= ps['meand'])
    r = np.polyfit(per[sel, 0], per[sel, 2], 1, w=np.sqrt(per[sel, 1]))[0] if sel.sum() >= 2 else np.nan
    c2 = (k * (k - 1.0)).mean() / k.mean()
    return c2, ps['meand'], ps['mu'], ps['lnmu'], r

f = open(os.path.join(OUT, 'paths2_scaling.csv'), 'w')
f.write('model,T,seed,c2,meand,mu,lnmu,rate\n')
for T in (1000, 3000, 10000, 30000, 100000):
    for s in range(1, 7):
        try:
            for name, h in (('SR', srcells.run(T, seed=s)['hit']), ('DM', surrogate.run_dm(T, seed=s)['hit'])):
                A = gs.adjacency(edges_from_hits(h), T)
                todo = [(name, A)]
                if name == 'SR':
                    E2, _ = gs.configuration_rewire(A, gs.degrees(A), seed=s, nswap_per_edge=5)
                    todo.append(('SR*', gs.adjacency(E2, T)))
                for nm, AA in todo:
                    v = pstats(AA, s)
                    f.write(f"{nm},{T},{s}," + ",".join(f"{x:.4f}" for x in v) + "\n"); f.flush()
        except Exception as ex:
            f.write(f"# FAILED {T},{s}: {ex}\n"); f.flush()
    print('A', T, flush=True)
f.close()

g = open(os.path.join(OUT, 'paths2_control.csv'), 'w')
g.write('variant,seed,tri,lnmu,lnmu_rewired,delta,meand,meand_rewired,C\n')
T = 30000
variants = {
    'SR': lambda s: srcells.run(T, seed=s)['hit'],
    'SR_q0': lambda s: srcells.run(T, q=0.0, seed=s)['hit'],
    'SR_0_60_120': lambda s: srcells.run(T, angle='set', angles=[0, math.pi / 3, 2 * math.pi / 3], seed=s)['hit'],
    'SR_0_90': lambda s: srcells.run(T, angle='set', angles=[0, math.pi / 2], seed=s)['hit'],
}
for p4 in (0.1, 0.29, 0.6, 0.9):
    variants[f'surr_p4_{p4}'] = (lambda p4: (lambda s: surrogate.run(T, face='uniform', side='padj', padj=p4, seed=s)['hit']))(p4)
for name, gen in variants.items():
    for s in range(1, 7):
        try:
            h = gen(s); A = gs.adjacency(edges_from_hits(h), T); k = gs.degrees(A)
            C, _ = gs.local_clustering(A, k)
            v = pstats(A, s); E2, _ = gs.configuration_rewire(A, k, seed=s, nswap_per_edge=5)
            w = pstats(gs.adjacency(E2, T), s)
            g.write(f"{name},{s},{tri_frac(h):.4f},{v[3]:.4f},{w[3]:.4f},{v[3]-w[3]:.4f},{v[1]:.3f},{w[1]:.3f},{C[k>=2].mean():.4f}\n"); g.flush()
        except Exception as ex:
            g.write(f"# FAILED {name},{s}: {ex}\n"); g.flush()
    print('B', name, flush=True)
g.close()
