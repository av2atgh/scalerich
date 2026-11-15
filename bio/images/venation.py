"""Leaf venation (LeafVeinCNN, 44 SLF leaves) as an iterated-division partition.
Tests: T share at 3-way nodes; at T, is the stem the thinnest arm; free ends;
lines (veins through collinear pairs) -> dual p_k, C(k), k vs width at fixed length;
HLD (width-ordered loop fusion) area split fractions vs SR (uniform, mean 0.25)."""
import sys, glob, json, numpy as np
sys.path.insert(0,'/Users/avazquez/av2atg/scalerich/bio/images')
from planar import classify_arms, lines, dual_stats
files=sorted(glob.glob('/Users/avazquez/av2atg/scalerich/bio/data/venation/*.npz'))
tol=float(sys.argv[1]) if len(sys.argv)>1 else 25.0
pool=dict(k=[],L=[],W=[],split=[],Tstem_thinnest=0,T=0,Y=0,X=0,deg=np.zeros(8,int),free_int=0,n1=0)
per=[]
for f in files:
    d=np.load(f)
    n1=d['e_EndNodes_1'].astype(int)-1; n2=d['e_EndNodes_2'].astype(int)-1
    W=d['e_Width']; L=d['e_Length']; bnd=(d['e_Ai']==0)|(d['e_Aj']==0)
    W=np.where(W>0,W,np.nan)
    def u(a): a=np.radians(a); return np.c_[-np.sin(a),np.cos(a)]   # (dy,dx) image coords
    ui=u(d['e_Or_ij']); uj=u(d['e_Or_ji'])
    arms={}
    for e in range(len(n1)):
        if n1[e]==n2[e]: continue
        arms.setdefault(n1[e],[]).append((e,ui[e])); arms.setdefault(n2[e],[]).append((e,uj[e]))
    deg=np.array([len(arms.get(i,[])) for i in range(len(d['n_node_ID']))])
    pool['deg']+=np.bincount(np.minimum(deg,7),minlength=8)
    recs=classify_arms(arms,tol)
    # degree-2 nodes: pass-through
    for j,a in arms.items():
        if len(a)==2: recs.append(dict(j=j,m=2,kind='2',through=[(a[0][0],a[1][0])],stems=[]))
    nodebnd=np.zeros(len(deg),bool); nodebnd[n1[bnd]]=True; nodebnd[n2[bnd]]=True
    for r in recs:
        if r['m']==3 and not nodebnd[r['j']]:
            if r['kind']=='T':
                pool['T']+=1; (p,q),s=r['through'][0],r['stems'][0]
                if W[s]<min(W[p],W[q]): pool['Tstem_thinnest']+=1
            else: pool['Y']+=1
        if r['kind']=='X': pool['X']+=1
    ends=np.where(deg==1)[0]; pool['n1']+=len(ends); pool['free_int']+=int((~nodebnd[ends]).sum())
    li=lines(len(n1),recs,width=np.nan_to_num(W,nan=np.nanmedian(W)),length=L)
    lb=np.bincount(li['lid'],weights=bnd.astype(float),minlength=li['nL'])>0
    st,k=dual_stats(li['nL'],li['D'],li['L'],li['W'])
    m=~lb
    pool['k']+=list(k[m]); pool['L']+=list(li['L'][m]); pool['W']+=list(li['W'][m])
    # HLD splits
    hid=d['h_node_ID'].astype(int); A=d['h_node_Area']; B=d['h_Boundary']; P=len(d['p_polygon_ID'])
    ch={}
    for c,p in zip(d['he_EndNodes_1'].astype(int),d['he_EndNodes_2'].astype(int)): ch.setdefault(p,[]).append(c)
    for p,cs in ch.items():
        if len(cs)==2 and B[cs[0]-1]==0 and B[cs[1]-1]==0:
            a1,a2=A[cs[0]-1],A[cs[1]-1]; pool['split'].append(min(a1,a2)/(a1+a2))
    per.append(dict(leaf=f.split('/')[-1][:-4],lines=int(m.sum()),C=st['C'],kmax=int(k[m].max())))
k=np.array(pool['k']); L=np.array(pool['L']); W=np.array(pool['W']); sp=np.array(pool['split'])
st,_=dual_stats(len(k),np.zeros((0,2),int))
out=dict(tol=tol,leaves=len(files),node_degree_counts=pool['deg'].tolist(),
  free_end_fraction_of_nodes=pool['n1']/pool['deg'].sum(),interior_free_ends=pool['free_int'],
  T=pool['T'],Y=pool['Y'],X=pool['X'],T_frac=pool['T']/(pool['T']+pool['Y']),
  T_stem_thinnest=pool['Tstem_thinnest']/pool['T'],
  interior_lines=len(k),kmean=float(k.mean()),ccdf={q:float((k>=q).mean()) for q in (1,2,4,8,16,32,64)},
  kmax=int(k.max()),C_mean_per_leaf=float(np.nanmean([p['C'] for p in per])),
  gamma4=float(1+(k>=4).sum()/np.sum(np.log(k[k>=4]/3.5))),
  corr_k_L=float(np.corrcoef(k,L)[0,1]))
s=(k>0)&(W>0)&(L>0)
X_=np.c_[np.ones(s.sum()),np.log(L[s])]
rk=np.log(k[s])-X_@np.linalg.lstsq(X_,np.log(k[s]),rcond=None)[0]
rw=np.log(W[s])-X_@np.linalg.lstsq(X_,np.log(W[s]),rcond=None)[0]
out['partial_logk_logW_given_logL']=float(np.corrcoef(rk,rw)[0,1])
out['slope_logk_logL']=float(np.polyfit(np.log(L[s]),np.log(k[s]),1)[0])
out['HLD_splits']=dict(n=len(sp),mean=float(sp.mean()),quartiles=np.percentile(sp,[25,50,75]).round(3).tolist(),
  hist_0_05=np.histogram(sp,bins=10,range=(0,0.5))[0].tolist())
print(json.dumps(out,indent=1))
np.savez('/Users/avazquez/av2atg/scalerich/bio/results/venation_pool.npz',k=k,L=L,W=W,split=sp)
json.dump(dict(agg=out,per=per),open(f'/Users/avazquez/av2atg/scalerich/bio/results/venation_tol{tol:g}.json','w'),indent=1)
