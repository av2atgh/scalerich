"""Sanity checks: Pi->0 reproduces SR statistics; geometry consistency."""
import sys, time, numpy as np
sys.path.insert(0, '.')
from pargrow import run
for Pi in [1e-3, 1.0, 1e5]:
    t=time.time(); r=run(3000, Pi, seed=1); dt=time.time()-t
    E=r['E']; n=len(r['P']); k=np.bincount(E.ravel(), minlength=n)
    free=(r['hitp']==-2).sum()+(r['hitm']==-2).sum()
    # geometric check: endpoint of tip lies on hit segment within its extent
    bad=0
    for i in range(n):
        for sd,h,s in ((1,r['hitp'][i],r['sp'][i]),(-1,r['hitm'][i],r['sm'][i])):
            if h<0: continue
            x=r['P'][i]+sd*s*r['D'][i]; w=x-r['P'][h]; u=w@r['D'][h]
            perp=abs(w[0]*r['D'][h][1]-w[1]*r['D'][h][0])
            ext=r['sp'][h] if u>=0 else r['sm'][h]
            if perp>1e-9 or abs(u)>ext+1e-9: bad+=1
    # crossings: any pair of final segments crossing in their interiors = X junction (should be 0)
    A=r['P']-r['sm'][:,None]*r['D']; B=r['P']+r['sp'][:,None]*r['D']
    X=0
    idx=np.random.default_rng(0).choice(n,300,replace=False)
    for i in idx:
        e=B-A; d=B[i]-A[i]; w=A-A[i]
        den=d[0]*e[:,1]-d[1]*e[:,0]
        with np.errstate(all='ignore'):
            s=(w[:,0]*e[:,1]-w[:,1]*e[:,0])/den; u=(w[:,0]*d[1]-w[:,1]*d[0])/den
        X+=((s>1e-7)&(s<1-1e-7)&(u>1e-7)&(u<1-1e-7)).sum()
    print(f"Pi={Pi:g} time={dt:.1f}s steps={r['steps']} <k>={k.mean():.3f} kmax={k.max()} free tips={free} bad hits={bad} X-crossings(300 sampled)={X}")
