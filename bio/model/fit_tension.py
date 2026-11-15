"""Fit new-wall tension maturation (t0, tau) to the meristem through-angle vs age curve.
Model time unit = 1/(per-cell division rate) = Tcc/ln2, Tcc = 36 h (meristem, 4-h division fraction ~0.08)."""
import sys, json, numpy as np, itertools, time
from plantdiv import Tissue
Tcc=36.0; unit=Tcc/np.log(2)
data_age=np.array([0,4,8,12,16,20,24,28,32,36])+2.0   # observed-frame age + mean offset
data_ang=np.array([143.3,140.0,138.4,136.5,133.6,131.6,131.0,130.0,127.5,128.0])  # pooled medians (6 plants)
def curve(T):
    ja=np.array([(a*unit,g) for a,g,_ in T.junction_angles()])
    out=[]
    for a in data_age:
        s=np.abs(ja[:,0]-a)<2.0
        out.append(np.median(ja[s,1]) if s.sum()>=5 else np.nan)
    return np.array(out), ja
res=[]
N=int(sys.argv[1]) if len(sys.argv)>1 else 600
for t0,tau_h in itertools.product([0.3,0.45,0.6,0.75],[8,24,48,96]):
    cs=[]
    for sd in range(2):
        T=Tissue(seed=sd,rule='errera',relax=True,t0=t0,tau=tau_h/unit).grow(N); c,_=curve(T); cs.append(c)
    c=np.nanmean(cs,0); rms=float(np.sqrt(np.nanmean((c-data_ang)**2)))
    res.append(dict(t0=t0,tau_h=tau_h,curve=c.round(1).tolist(),rms=rms)); print(t0,tau_h,c.round(1),f"rms {rms:.2f}",flush=True)
json.dump(res,open('../results/fit_tension.json','w'),indent=1)
best=min(res,key=lambda r:r['rms']); print('best',best)
