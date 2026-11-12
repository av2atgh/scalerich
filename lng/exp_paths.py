"""Size scaling of distances and shortest-path multiplicity.

Models: SR (geometric), SR* (degree-preserving rewiring of the SR dual),
U-padj (topological face surrogate, p4 = 0.29), DM (BB(n,1)).
Output: results/paths_scaling.csv  (model,T,seed,c2,meand,mu,lnmu,Cmean)
        results/paths_mu_of_d.csv  (model,T,d,npairs,mu,lnmu) for the largest T
Baseline (book, distance.tex): <mu> ~ c2/e for configuration networks.
"""
import sys, os, time
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import srcells, surrogate, graphstats as gs
from exp_topology import edges_from_hits

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'results')
Ts = [1000, 3000, 10000, 30000, 100000, 300000]
gens = {
    'SR': lambda T, s: srcells.run(T, seed=s)['hit'],
    'U-padj': lambda T, s: surrogate.run(T, face='uniform', side='padj', padj=0.29, seed=s)['hit'],
    'DM': lambda T, s: surrogate.run_dm(T, seed=s)['hit'],
}
f = open(os.path.join(OUT, 'paths_scaling.csv'), 'w')
f.write('model,T,seed,c2,meand,mu,lnmu,Cmean\n')
g = open(os.path.join(OUT, 'paths_mu_of_d.csv'), 'w')
g.write('model,T,d,npairs,mu,lnmu\n')
for T in Ts:
    seeds = [1, 2, 3] if T <= 30000 else [1, 2]
    nsrc = 60 if T <= 100000 else 30
    for name, gen in gens.items():
        for s in seeds:
            t0 = time.time()
            try:
                h = gen(T, s)
                E = edges_from_hits(h)
                A = gs.adjacency(E, T)
                todo = [(name, A)]
                if name == 'SR' and T <= 100000:
                    k = gs.degrees(A)
                    E2, acc = gs.configuration_rewire(A, k, seed=s, nswap_per_edge=5)
                    todo.append(('SR*', gs.adjacency(E2, T)))
                for nm, AA in todo:
                    k = gs.degrees(AA)
                    C, _ = gs.local_clustering(AA, k)
                    ps = gs.path_stats(AA, nsrc=nsrc, seed=s)
                    c2 = (k * (k - 1.0)).mean() / k.mean()
                    f.write(f"{nm},{T},{s},{c2:.4f},{ps['meand']:.4f},{ps['mu']:.4f},{ps['lnmu']:.4f},{C[k>=2].mean():.4f}\n")
                    f.flush()
                    if T == Ts[-1] or (T == 100000 and nm == 'SR*'):
                        if s == 1:
                            for d, n, m, lm in ps['per']:
                                g.write(f"{nm},{T},{d},{n},{m:.4f},{lm:.4f}\n")
                            g.flush()
                print(name, T, s, f"{time.time()-t0:.1f}s", flush=True)
            except Exception as ex:  # record and move on
                f.write(f"# FAILED {name},{T},{s}: {ex}\n")
                print('FAILED', name, T, s, ex, flush=True)
f.close(); g.close()
