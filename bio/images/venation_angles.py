import glob, numpy as np, sys
sys.path.insert(0,'../../beyond')
files=sorted(glob.glob('/Users/avazquez/av2atg/scalerich/bio/data/venation/*.npz'))
H={'loop':[],'tree':[]}; thin={'loop':[0,0],'tree':[0,0]}
for f in files:
    d=np.load(f); n1=d['e_EndNodes_1'].astype(int)-1; n2=d['e_EndNodes_2'].astype(int)-1
    bnd=(d['e_Ai']==0)|(d['e_Aj']==0); typ=d['e_Type']; W=d['e_Width']
    def u(a): a=np.radians(a); return np.c_[-np.sin(a),np.cos(a)]
    ui=u(d['e_Or_ij']); uj=u(d['e_Or_ji']); arms={}
    for e in range(len(n1)):
        if n1[e]==n2[e]: continue
        arms.setdefault(n1[e],[]).append((e,ui[e])); arms.setdefault(n2[e],[]).append((e,uj[e]))
    for j,a in arms.items():
        if len(a)!=3: continue
        es=[x[0] for x in a]
        if bnd[es].any(): continue
        V=np.array([x[1] for x in a]); C=V@V.T; iu=np.triu_indices(3,1)
        p=np.argmin(C[iu]); amax=np.degrees(np.arccos(np.clip(C[iu][p],-1,1)))
        key='loop' if typ[es].all() else 'tree'
        H[key].append(amax)
        stem=3-iu[0][p]-iu[1][p]
        thin[key][1]+=1; thin[key][0]+= W[es[stem]]<=min(W[es[iu[0][p]]],W[es[iu[1][p]]])
for k,v in H.items():
    v=np.array(v); print(k,'n',len(v),'median max-angle',np.median(v).round(1),'frac>155',(v>155).mean().round(3),'frac>165',(v>165).mean().round(3),
      'hist 120..180 by 10:',np.histogram(v,bins=[0,120,130,140,150,160,170,180])[0].tolist(),
      'stem(opposite widest gap) thinnest', round(thin[k][0]/thin[k][1],3))
# model split fractions
import divisions
for rule in ['SR','centroid','balanced']:
    fr=[]
    for sd in range(3):
        r=divisions.run(3000,seed=sd,rule=rule); x=np.asarray(r['frac']); x=x[np.isfinite(x)]; fr+=list(np.minimum(x,1-x))
    fr=np.array(fr); print(rule,'split mean',fr.mean().round(3),'quartiles',np.percentile(fr,[25,50,75]).round(3),'hist',np.histogram(fr,bins=10,range=(0,.5))[0].tolist())
