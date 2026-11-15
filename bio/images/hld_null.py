"""Width-ordered loop fusion (HLD) split fractions: SR (widths by birth, +- noise) vs Voronoi (random widths)."""
import sys, warnings, numpy as np; warnings.filterwarnings('ignore')
sys.path.insert(0,'.'); sys.path.insert(0,'../../beyond')
import divisions
from skimage.draw import polygon
from scipy.spatial import cKDTree
from planar import label_graph
def hld(walls, area):
    par={c:c for c in area}; A=dict(area); sp=[]
    def f(x):
        while par[x]!=x: par[x]=par[par[x]]; x=par[x]
        return x
    for a,b,w in sorted(walls,key=lambda x:x[2]):
        ra,rb=f(a),f(b)
        if ra==rb: continue
        sp.append(min(A[ra],A[rb])/(A[ra]+A[rb])); par[ra]=rb; A[rb]+=A[ra]
    return np.array(sp)
def report(name,sp):
    print(f"{name:28s} mean {sp.mean():.3f} quartiles {np.percentile(sp,[25,50,75]).round(3)} hist {np.histogram(sp,bins=10,range=(0,.5))[0].tolist()}",flush=True)
M=1536; rng=np.random.default_rng(0); T=800
r=divisions.run(T,seed=2,rule='SR')
Lab=np.zeros((M,M),np.int32)
for i,p in enumerate(r['leaf_polys']):
    p=np.asarray(p)*(M-1); rr,cc=polygon(p[:,1],p[:,0],(M,M)); Lab[rr,cc]=i+1
V,E,vl=label_graph(Lab)
A=r['A']*(M-1); B=r['B']*(M-1)
area={c:float(n) for c,n in zip(*np.unique(Lab[Lab>0],return_counts=True))}
births=[]
for u,v,pts,(a,b) in E:
    if a==0: continue
    m=pts[len(pts)//2][::-1]  # (x,y)
    d=B-A; t=np.clip(((m-A)*d).sum(1)/(d*d).sum(1),0,1); dist=np.hypot(*(A+t[:,None]*d-m).T)
    births.append((a,b,int(np.argmin(dist))))
for sig in [0,0.5,1,2]:
    w=[(a,b,-(np.log(t+1)) + sig*rng.standard_normal()) for a,b,t in births]   # older = wider
    report(f"SR widths by birth, noise {sig}",hld(w,area))
# Voronoi null
seeds=rng.random((T,2))*(M-1); yy,xx=np.mgrid[0:M,0:M]
_,lab=cKDTree(seeds).query(np.c_[xx.ravel(),yy.ravel()]); Lv=(lab+1).reshape(M,M).astype(np.int32)
V,E,vl=label_graph(Lv); area={c:float(n) for c,n in zip(*np.unique(Lv,return_counts=True))}
w=[(a,b,rng.random()) for _,_,_,(a,b) in E if a>0]
report("Voronoi, random widths",hld(w,area))
w=[(a,b,np.hypot(*(pts[0]-pts[-1]))) for _,_,pts,(a,b) in E if a>0]
report("Voronoi, width = wall length",hld(w,area))
