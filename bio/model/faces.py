"""Reference face statistics for iterated division rules: geometric corners per cell
(interior angle < 180 - tol) vs topological neighbours; areas."""
import sys, numpy as np
sys.path.insert(0,'../../beyond')
import divisions
def corners(poly,tol=10):
    v=poly; a=np.roll(v,1,0)-v; b=np.roll(v,-1,0)-v
    na=np.linalg.norm(a,axis=1); nb=np.linalg.norm(b,axis=1); ok=(na>1e-12)&(nb>1e-12)
    cos=(a*b).sum(1)/np.where(ok,na*nb,1)
    ang=np.degrees(np.arccos(np.clip(cos,-1,1)))
    return int((ok&(ang<180-tol)).sum()), int(ok.sum())
for rule in ['SR','uniform','centroid','balanced','regular']:
    C=[];V=[];A=[]
    for sd in range(2):
        r=divisions.run(5000,seed=sd,rule=rule)
        for p,a in zip(r['leaf_polys'],r['leaf_area']):
            if a<=0: continue
            c,v=corners(np.asarray(p)); C.append(c); V.append(v); A.append(a)
    C=np.array(C);V=np.array(V);A=np.array(A); la=np.log(A)
    print(f"{rule:9s} corners {C.mean():.2f}  vertices(incl T-flat) {V.mean():.2f}  P(corners=4) {np.mean(C==4):.2f} sd(ln area) {la.std():.2f}")
