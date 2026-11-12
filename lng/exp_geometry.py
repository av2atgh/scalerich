"""What geometry adds: thickness/jamming, angle rules, nucleation rules, cell weight q.

For each configuration (T = 30000 or until the free area is exhausted, 2 seeds):
  Tdone, final fill fraction rho = 1 - free area, gamma (MLE k>=6, k>=10),
  K^2 P(k-2 >= K) at K = 8, 16, cohort exponent beta (mean k of birth cohorts
  ~ i^-beta), C, P(corner cut), P(triangle), c2, <d>, <mu>.
Output: results/geometry.csv
"""
import sys, os, time, math
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import srcells, graphstats as gs
from exp_topology import edges_from_hits

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'results')
T = 30000
circ = lambda x, y: math.atan2(y - 0.5, x - 0.5) + math.pi / 2
CONFIGS = [
    ('SR0', dict()),
    # tunable-density phase
    ('a1.0_l0.1', dict(alpha=1.0, lam0=0.1)),
    ('a0.75_l0.03', dict(alpha=0.75, lam0=0.03)),
    ('a0.75_l0.1', dict(alpha=0.75, lam0=0.1)),
    ('a0.6_l0.03', dict(alpha=0.6, lam0=0.03)),
    ('a0.5_l0.01', dict(alpha=0.5, lam0=0.01)),
    ('a0.5_l0.03', dict(alpha=0.5, lam0=0.03)),
    # jammed phase
    ('a0.4_l0.01', dict(alpha=0.4, lam0=0.01)),
    ('a0.4_l0.02', dict(alpha=0.4, lam0=0.02)),
    ('a0.25_l0.005', dict(alpha=0.25, lam0=0.005)),
    ('a0.25_l0.01', dict(alpha=0.25, lam0=0.01)),
    ('a0_l0.002', dict(alpha=0.0, lam0=0.002)),
    ('a0_l0.005', dict(alpha=0.0, lam0=0.005)),
    ('a0_l0.02', dict(alpha=0.0, lam0=0.02)),
    ('a0_l0.05', dict(alpha=0.0, lam0=0.05)),
    # angle rules
    ('set_0_90', dict(angle='set', angles=[0.0, math.pi / 2])),
    ('set_0_60_120', dict(angle='set', angles=[0.0, math.pi / 3, 2 * math.pi / 3])),
    ('vm_k1', dict(angle='vonmises', kappa=1.0)),
    ('vm_k4', dict(angle='vonmises', kappa=4.0)),
    ('vm_k16', dict(angle='vonmises', kappa=16.0)),
    ('field_circ', dict(angle='field', field=circ, kappa=0.0)),
    ('field_circ_k4', dict(angle='field', field=circ, kappa=4.0)),
    # nucleation / cell weight
    ('q0', dict(q=0.0)),
    ('q0.5', dict(q=0.5)),
    ('q2', dict(q=2.0)),
    ('grad_x', dict(dens=lambda x, y: 0.05 + x)),
    ('grad_r', dict(dens=lambda x, y: 0.02 + math.hypot(x - 0.5, y - 0.5))),
    ('a0.5_l0.03_set0_90', dict(alpha=0.5, lam0=0.03, angle='set', angles=[0.0, math.pi / 2])),
]

if __name__ == '__main__':
    only = sys.argv[1:]  # optional subset of names
    path = os.path.join(OUT, 'geometry.csv')
    new = not os.path.exists(path) or not only
    f = open(path, 'w' if new else 'a')
    if new:
        f.write('name,seed,Tdone,rho,gmle6,gmle10,K2c8,K2c16,beta,Cmean,P_corner,P_tri,c2,meand,mu,lnmu,kmax,t_rho99\n')
    for name, kw in CONFIGS:
        if only and name not in only:
            continue
        for s in (1, 2):
            t0 = time.time()
            try:
                o = srcells.run(T, seed=s, **kw)
                Td = o['T']
                h = o['hit']
                E = edges_from_hits(h)
                A = gs.adjacency(E, Td)
                k = gs.degrees(A)
                C, _ = gs.local_clustering(A, k)
                nb = (h >= 0).all(1)
                pc = o['adj'][nb].mean()
                pt = (np.asarray(A[h[nb, 0], h[nb, 1]]).ravel() > 0).mean()
                g6 = gs.gamma_mle(k, 6)[0]; g10 = gs.gamma_mle(k, 10)[0]
                K8 = 64 * ((k - 2) >= 8).mean(); K16 = 256 * ((k - 2) >= 16).mean()
                # cohort exponent
                ii = np.unique(np.logspace(np.log10(max(10, Td / 3000)), np.log10(Td / 10), 12).astype(int))
                mk = [k[a:int(a * 1.3) + 1].mean() for a in ii]
                beta = -np.polyfit(np.log(ii), np.log(np.maximum(mk, 1e-9)), 1)[0]
                ps = gs.path_stats(A, nsrc=20, seed=s)
                rho = 1 - o['freeA'][-1] if Td else float('nan')
                jj = np.nonzero(o['freeA'] < 0.01)[0]
                t99 = int(jj[0]) if len(jj) and o['freeA'][-1] > -1 and kw.get('lam0', 0) > 0 else -1
                c2 = (k * (k - 1.0)).mean() / k.mean()
                f.write(f"{name},{s},{Td},{rho:.4f},{g6:.3f},{g10:.3f},{K8:.3f},{K16:.3f},{beta:.3f},"
                        f"{C[k>=2].mean():.4f},{pc:.4f},{pt:.4f},{c2:.3f},{ps['meand']:.3f},{ps['mu']:.3f},{ps['lnmu']:.3f},{k.max()},{t99}\n")
                f.flush()
                print(name, s, Td, f"{time.time()-t0:.1f}s", flush=True)
            except Exception as ex:
                f.write(f"# FAILED {name},{s}: {ex!r}\n"); f.flush()
                print('FAILED', name, s, repr(ex), flush=True)
    f.close()
