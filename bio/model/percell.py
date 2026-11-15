"""Per-compartment clocks: cell picked uniformly (each compartment divides at the
same rate, regardless of area) vs SR (picked prop. to area), plus balanced/centroid
(Errera-like). Dual graph: segment t joins its two host segments (hostA, hostB >= 0)."""
import sys, numpy as np
sys.path.insert(0,'../../beyond'); sys.path.insert(0,'../../lng')
import divisions, graphstats as gs
T=20000
for rule in ['SR','uniform','centroid','balanced']:
    out=[]
    for sd in range(2):
        r=divisions.run(T,seed=sd,rule=rule)
        E=[(t,h) for t in range(T) for h in (r['hostA'][t],r['hostB'][t]) if h>=0]
        A=gs.adjacency(E,T); k=gs.degrees(A); C,_=gs.local_clustering(A,k)
        L=np.linalg.norm(r['B']-r['A'],axis=1)
        out.append((gs.gamma_mle(k,10)[0],gs.gamma_mle(k,6)[0],k.max(),C[k>=2].mean(),[float((k>=q).mean()) for q in (8,16,32,64)],np.std(np.log(L))))
    o=np.array([x[:4] for x in out]).mean(0)
    print(f"{rule:9s} gamma10={o[0]:.2f} gamma6={o[1]:.2f} kmax={o[2]:.0f} C={o[3]:.3f} ccdf(8,16,32,64)={np.round(out[0][4],4)} sd(lnL)={out[0][5]:.2f}",flush=True)
