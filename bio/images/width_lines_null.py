"""Null for width-continuation lines: same construction on Voronoi (random / length widths)
and on SR (birth widths, +noise), from rasterised label images; boundary lines excluded."""
import sys, warnings, numpy as np; warnings.filterwarnings('ignore')
sys.path.insert(0,'.'); sys.path.insert(0,'../../beyond')
import divisions
from skimage.draw import polygon
from scipy.spatial import cKDTree
from planar import label_graph, lines, dual_stats
rng=np.random.default_rng(1); M=1536; T=800
def run(name,V,E,w):
    idx=[i for i,e in enumerate(E)]
    n1=np.array([e[0] for e in E]); n2=np.array([e[1] for e in E]); bnd=np.array([e[3][0]==0 for e in E])
    L=np.array([np.hypot(*np.diff(e[2],axis=0).T).sum() for e in E])
    inc={}
    for e in range(len(E)):
        if bnd[e]: continue
        inc.setdefault(n1[e],[]).append(e); inc.setdefault(n2[e],[]).append(e)
    recs=[]
    for j,es in inc.items():
        if len(es)<2: continue
        o=sorted(es,key=lambda e:-w[e]); recs.append(dict(j=j,m=len(es),through=[(o[0],o[1])],stems=o[2:]))
    li=lines(len(E),recs,width=w,length=L)
    lb=np.bincount(li['lid'],weights=bnd.astype(float),minlength=li['nL'])>0
    # lines touching a boundary vertex also excluded
    bv=set(n1[bnd])|set(n2[bnd]); tb=np.array([(n1[e] in bv) or (n2[e] in bv) for e in range(len(E))])
    lb|=np.bincount(li['lid'],weights=tb.astype(float),minlength=li['nL'])>0
    st,k=dual_stats(li['nL'],li['D'],li['L'],li['W']); m=~lb; kk=k[m]
    s=(kk>0); X=np.c_[np.ones(s.sum()),np.log(li['L'][m][s])]
    rk=np.log(kk[s])-X@np.linalg.lstsq(X,np.log(kk[s]),rcond=None)[0]; rw=np.log(li['W'][m][s])-X@np.linalg.lstsq(X,np.log(li['W'][m][s]),rcond=None)[0]
    print(f"{name:30s} lines={m.sum()} ccdf",{q:round(float((kk>=q).mean()),3) for q in (2,4,8,16,32)},"C",round(st['C'],3),"C(k)",st['Ck'][:5],f"pc(k,W|L)={np.corrcoef(rk,rw)[0,1]:.2f}",flush=True)
r=divisions.run(T,seed=2,rule='SR'); Lab=np.zeros((M,M),np.int32)
for i,p in enumerate(r['leaf_polys']):
    p=np.asarray(p)*(M-1); rr,cc=polygon(p[:,1],p[:,0],(M,M)); Lab[rr,cc]=i+1
pad=np.zeros((M+8,M+8),np.int32); pad[4:-4,4:-4]=Lab
V,E,vl=label_graph(pad)
A=r['A']*(M-1)+4; B=r['B']*(M-1)+4; birth=[]
for u,v,pts,(a,b) in E:
    m=pts[len(pts)//2][::-1]; d=B-A; t=np.clip(((m-A)*d).sum(1)/(d*d).sum(1),0,1)
    birth.append(np.argmin(np.hypot(*(A+t[:,None]*d-m).T)))
birth=np.array(birth)
for sig in (0,1):
    run(f"SR birth widths noise {sig}",V,E,np.exp(-0.5*np.log(birth+1)+sig*rng.standard_normal(len(E))))
run("SR random widths",V,E,rng.random(len(E)))
seeds=rng.random((T,2))*(M-1); yy,xx=np.mgrid[0:M,0:M]
_,lab=cKDTree(seeds).query(np.c_[xx.ravel(),yy.ravel()]); Lv=np.zeros((M+8,M+8),np.int32); Lv[4:-4,4:-4]=(lab+1).reshape(M,M)
V,E,vl=label_graph(Lv)
run("Voronoi random widths",V,E,rng.random(len(E)))
Lw=np.array([np.hypot(*np.diff(e[2],axis=0).T).sum() for e in E]); run("Voronoi width=length",V,E,Lw)
