"""Predictions of the plant model (Errera + no corner cut + tension maturation t0=0.7, tau=36 h)
vs no relaxation and SR; P(n), dual degree law, clustering. N cells, seeds."""
import sys, json, numpy as np
from plantdiv import Tissue
N=int(sys.argv[1]); S=int(sys.argv[2]); unit=36/np.log(2)
cfg={'plant (errera+relax, size pick)':dict(rule='errera',relax=True,t0=0.7,tau=36/unit,pick='size'),
     'plant, area pick':dict(rule='errera',relax=True,t0=0.7,tau=36/unit,pick='area'),
     'errera, no relax':dict(rule='errera',relax=False,pick='size'),
     'SR (no relax)':dict(rule='sr',relax=False,pick='area')}
out={}
for name,kw in cfg.items():
    Pn=[];K=[];C=[];cc=0;nd=0
    for sd in range(S):
        T=Tissue(seed=sd,**kw).grow(N)
        Pn+=list(T.sides()); k,Cl=T.dual(); K+=list(k); C+=Cl; cc+=T.corner_cuts; nd+=len(T.birth)
    Pn=np.array(Pn); K=np.array(K); x=K[K>=6].astype(float)
    Ck={}
    for lo,hi in [(2,3),(3,5),(5,8),(8,13),(13,21),(21,999)]:
        v=[c for kk,c in C if lo<=kk<hi]
        if len(v)>=5: Ck[lo]=round(float(np.mean(v)),3)
    o=dict(Pn={n:round(float((Pn==n).mean()),3) for n in range(3,11)},corner_cut=cc/nd,
           C=float(np.mean([c for _,c in C])),Ck=Ck,kmax=int(K.max()),gamma6=float(1+len(x)/np.sum(np.log(x/5.5))),
           ccdf={q:round(float((K>=q).mean()),4) for q in (4,8,16,32,64)})
    out[name]=o; print(name,o,flush=True)
json.dump(out,open(f'../results/plant_predictions_N{N}.json','w'),indent=1)
