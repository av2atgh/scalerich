"""SR dual vs topological surrogates: P(k), C(k), k_nn(k), paths, multiplicity.

Usage: python exp_topology.py  -> results/topology_*.txt
"""
import sys, os, json, time
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import srcells, surrogate, graphstats as gs

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'results')
os.makedirs(OUT, exist_ok=True)


def edges_from_hits(h):
    T = len(h)
    e = []
    for col in (0, 1):
        m = h[:, col] >= 0
        e.append(np.stack([np.arange(T)[m], h[m, col]], 1))
    E = np.concatenate(e)
    return np.unique(np.sort(E, 1), axis=0)


def measure(h, adjflag, T, nsrc=40, seed=0, paths=True):
    E = edges_from_hits(h)
    A = gs.adjacency(E, T)
    r, k, C, kn = gs.summary(A, nsrc=nsrc, seed=seed, paths=paths)
    nb = (h >= 0).all(1)
    r['P_celladj'] = float(adjflag[nb].mean())
    r['P_tri'] = float((np.asarray(A[h[nb, 0], h[nb, 1]]).ravel() > 0).mean())
    # CCDF amplitude K^2 P(k-2 >= K) at a few K
    for K in (8, 16, 32):
        r[f'K2ccdf{K}'] = float(K * K * ((k - 2) >= K).mean())
    return r, k


MODELS = {
    'SR':            lambda T, s: srcells.run(T, seed=s),
    'U-uniform':     lambda T, s: surrogate.run(T, face='uniform', side='uniform', seed=s),
    'U-padj':        lambda T, s: surrogate.run(T, face='uniform', side='padj', padj=0.29, seed=s),
    'N-padj':        lambda T, s: surrogate.run(T, face='nsides', side='padj', padj=0.29, seed=s),
    'A-padj':        lambda T, s: surrogate.run(T, face='area', side='padj', padj=0.29, seed=s),
    'DM':            lambda T, s: surrogate.run_dm(T, seed=s),
    'A-length':      lambda T, s: surrogate.run(T, face='area', side='length', seed=s),
}

if __name__ == '__main__':
    T = int(sys.argv[1]) if len(sys.argv) > 1 else 100000
    seeds = [1, 2, 3]
    rows = []
    with open(os.path.join(OUT, f'topology_T{T}.txt'), 'w') as f:
        for name, fn in MODELS.items():
            acc = []
            for s in seeds:
                t0 = time.time()
                o = fn(T, s)
                r, k = measure(o['hit'], o['adj'], o['T'], nsrc=40, seed=s)
                r['time'] = time.time() - t0
                acc.append(r)
            keys = ['meank', 'kmax', 'c2', 'Cmean', 'transitivity', 'r', 'gmle6', 'gmle10',
                    'K2ccdf8', 'K2ccdf16', 'K2ccdf32', 'P_celladj', 'P_tri', 'meand', 'mu', 'lnmu']
            line = name + '  ' + '  '.join(f"{kk}={np.mean([a[kk] for a in acc]):.3f}" for kk in keys)
            print(line, flush=True)
            f.write(line + '\n')
            # C(k), knn(k) of seed 1
            f.write('   C(k): ' + ' '.join(f"[{lo},{hi}):{v:.3f}" for lo, hi, n, v in acc[0]['Ck']) + '\n')
            f.write('   knn(k): ' + ' '.join(f"[{lo},{hi}):{v:.2f}" for lo, hi, n, v in acc[0]['knnk']) + '\n')
            f.write('   mu(d): ' + ' '.join(f"{d}:{m:.2f}" for d, n, m, lm in acc[0]['per'] if n > 100) + '\n')
            f.flush()
