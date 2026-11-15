"""Pipeline validation on rendered model networks with known dual graphs."""
import sys, json, numpy as np
sys.path.insert(0,'.'); sys.path.insert(0,'../model')
from pipeline import analyse, summary, render
from pargrow import run
out={}
def truth_k(r):
    n=len(r['P']); return np.bincount(r['E'].ravel(),minlength=n)
for name,Pi in [('SR',1e-3),('Gilbert',1e7)]:
    for N in (300,800):
        r=run(N,Pi,seed=1); A=r['P']-r['sm'][:,None]*r['D']; B=r['P']+r['sp'][:,None]*r['D']
        rank=np.arange(N); w=np.maximum(1.5,8*(rank+1)**-0.35)
        img=render(A,B,w,M=1024,psf=1.5,noise=0.05)
        s=summary(analyse(img)); kt=truth_k(r)
        s['truth']={f'P(k>={q})':float((kt>=q).mean()) for q in (2,4,8,16,32)}; s['truth']['kmax']=int(kt.max()); s['truth']['n']=N
        out[f'{name}_N{N}']=s
# Mikado: random fixed-length rods, X crossings, no T by construction (except chance)
rng=np.random.default_rng(2)
for N in (300,800):
    P=rng.random((N,2)); th=rng.random(N)*np.pi; d=np.c_[np.cos(th),np.sin(th)]*0.08
    A=(P-d).clip(0,1); B=(P+d).clip(0,1)
    img=render(A,B,np.full(N,2.0),M=1024); out[f'Mikado_N{N}']=summary(analyse(img))
for k,v in out.items():
    print(k, {a:(round(b,3) if isinstance(b,float) else b) for a,b in v.items()})
json.dump(out,open('../results/pipeline_validation.json','w'),indent=1)
