"""Venation dual graph with lines defined by width continuation: at each node the two
widest arms continue (through), the rest are stems. Variants: all edges / loop edges only."""
import sys, glob, json, numpy as np
sys.path.insert(0,'.')
from planar import lines, dual_stats
files=sorted(glob.glob('/Users/avazquez/av2atg/scalerich/bio/data/venation/*.npz'))
for mode in ['all','loops']:
    K=[];LL=[];WW=[];Cs=[];CK={}
    for f in files:
        d=np.load(f); n1=d['e_EndNodes_1'].astype(int)-1; n2=d['e_EndNodes_2'].astype(int)-1
        keep=(n1!=n2)&(d['e_Width']>0)
        if mode=='loops': keep&=d['e_Type']==1
        bnd=(d['e_Ai']==0)|(d['e_Aj']==0)
        idx=np.where(keep)[0]; W=d['e_Width']; L=d['e_Length']
        inc={}
        for e in idx: inc.setdefault(n1[e],[]).append(e); inc.setdefault(n2[e],[]).append(e)
        recs=[]
        for j,es in inc.items():
            if len(es)<2: continue
            o=sorted(es,key=lambda e:-W[e])
            recs.append(dict(j=j,m=len(es),kind='w',through=[(o[0],o[1])],stems=o[2:]))
        li=lines(len(n1),recs,width=W,length=L)
        used=np.zeros(li['nL'],bool); used[li['lid'][idx]]=True
        lb=np.bincount(li['lid'],weights=bnd.astype(float),minlength=li['nL'])>0
        st,k=dual_stats(li['nL'],li['D'],li['L'],li['W'])
        m=used&~lb
        K+=list(k[m]); LL+=list(li['L'][m]); WW+=list(li['W'][m]); Cs.append(st['C'])
        for lo,n,c in st['Ck']: CK.setdefault(lo,[]).append(c)
    k=np.array(K); L=np.array(LL); W=np.array(WW); s=(k>0)&(L>0)&(W>0)
    X=np.c_[np.ones(s.sum()),np.log(L[s])]
    rk=np.log(k[s])-X@np.linalg.lstsq(X,np.log(k[s]),rcond=None)[0]; rw=np.log(W[s])-X@np.linalg.lstsq(X,np.log(W[s]),rcond=None)[0]
    x=k[k>=4].astype(float)
    print(mode, f"lines={len(k)} <k>={k.mean():.2f} kmax={k.max()} gamma4={1+len(x)/np.sum(np.log(x/3.5)):.2f} gamma8={1+(k>=8).sum()/np.sum(np.log(k[k>=8]/7.5)):.2f} C={np.nanmean(Cs):.3f}",
          "ccdf", {q:round(float((k>=q).mean()),5) for q in (1,2,4,8,16,32,64)},
          "C(k)", {lo:round(float(np.mean(v)),3) for lo,v in sorted(CK.items())},
          f"corr(logk,logW|logL)={np.corrcoef(rk,rw)[0,1]:.3f} slope logk~logL={np.polyfit(np.log(L[s]),np.log(k[s]),1)[0]:.2f}")
    np.savez(f'/Users/avazquez/av2atg/scalerich/bio/results/venation_width_{mode}.npz',k=k,L=L,W=W)
