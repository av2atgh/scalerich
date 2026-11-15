"""Turnover: segments removed at rate delta; steady state n_ss = R/delta.
Statistics on the segments present at the end (after ~6 lifetimes)."""
import sys, json, time, numpy as np
sys.path.insert(0, '.'); sys.path.insert(0, '../../lng')
from pargrow import run
import graphstats as gs
nss=int(sys.argv[1]); seeds=int(sys.argv[2]); Pis=[float(x) for x in sys.argv[3].split(',')]
stab=float(sys.argv[4]) if len(sys.argv)>4 else 0.0
res={}
for Pi in Pis:
    R=Pi; delta=R/nss*{0:1.0,1:1.45,2:1.70,3:1.71}.get(int(stab),1.0); Tend=8/(R/nss); N=int(R*Tend*1.05)+10
    rows=[]
    for sd in range(seeds):
        t=time.time(); r=run(N,Pi,seed=sd,delta=delta,Tmax=Tend,stab=stab)
        keep=r['present']; n=len(keep)
        A=gs.adjacency(r['E'],n); k=gs.degrees(A); C,_=gs.local_clustering(A,k)
        kk=k[keep]; free=((r['hitp']==-2)&keep).sum()+((r['hitm']==-2)&keep).sum()
        age=Tend-r['t0'][keep]
        o=dict(n=int(keep.sum()),kmean=float(kk.mean()),kmax=int(kk.max()),
               gamma6=float(gs.gamma_mle(kk,6)[0]),C=float(C[keep&(k>=2)].mean()),
               free_tips_per_seg=float(free/keep.sum()),
               corr_k_age=float(np.corrcoef(kk,age)[0,1]),
               ccdf=[[q,float((kk>=q).mean())] for q in [2,4,8,16,32,64]],sec=time.time()-t)
        rows.append(o)
    res[Pi]=rows
    print(f"stab={stab} nss={nss} Pi={Pi:g} n={np.mean([o['n'] for o in rows]):.0f} <k>={np.mean([o['kmean'] for o in rows]):.2f} kmax={np.mean([o['kmax'] for o in rows]):.0f} gamma6={np.mean([o['gamma6'] for o in rows]):.2f} C={np.mean([o['C'] for o in rows]):.3f} free/seg={np.mean([o['free_tips_per_seg'] for o in rows]):.3f} corr(k,age)={np.mean([o['corr_k_age'] for o in rows]):.2f} ccdf={rows[0]['ccdf']} t={np.mean([o['sec'] for o in rows]):.0f}s",flush=True)
json.dump(res,open(f'../results/turnover_nss{nss}_stab{stab:g}.json','w'),indent=1)
