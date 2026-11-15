"""Pi scan: dual-graph statistics of parallel iterated division (N segments)."""
import sys, json, time, numpy as np
sys.path.insert(0, '.'); sys.path.insert(0, '../../lng')
from pargrow import run
import graphstats as gs

def stats(r):
    n=len(r['P']); keep=r['present']; E=r['E']
    A=gs.adjacency(E, n); k=gs.degrees(A); C,_=gs.local_clustering(A,k)
    L=r['sp']+r['sm']; rank=np.argsort(np.argsort(r['t0']))
    out=dict(n=int(keep.sum()), kmean=float(k[keep].mean()), kmax=int(k.max()))
    out['gamma6']=float(gs.gamma_mle(k[keep],6)[0]); out['gamma10']=float(gs.gamma_mle(k[keep],10)[0])
    out['C']=float(C[k>=2].mean())
    bins=[2,3,5,8,13,21,34,55,89,1000]; ck=[]
    for lo,hi in zip(bins[:-1],bins[1:]):
        s=(k>=lo)&(k<hi)
        if s.sum()>=5: ck.append([lo,hi,int(s.sum()),float(C[s].mean())])
    out['Ck']=ck
    # degree slaved to length x age: k-2 vs L*(sqrt(n)-sqrt(rank))
    x=L*(np.sqrt(n)-np.sqrt(rank+1)); s=(rank>50)&(rank<0.8*n)
    out['corr_k_Lage']=float(np.corrcoef(k[s]-2,x[s])[0,1])
    out['corr_logk_rank']=float(np.corrcoef(np.log(k[k>0]),np.log(rank[k>0]+1))[0,1])
    # cohort mean degree vs birth rank
    out['cohort']=[[lo,hi,float(k[(rank>=lo)&(rank<hi)].mean())] for lo,hi in [(10,30),(100,300),(1000,3000),(10000,20000)] if hi<=n]
    # triangle closure fraction among segments with two hits
    adj=[set() for _ in range(n)]
    for a,b in E: adj[a].add(b); adj[b].add(a)
    tri=nb=0
    for i in range(n):
        a,b=r['hitp'][i],r['hitm'][i]
        if a>=0 and b>=0: nb+=1; tri+= b in adj[a]
    out['tri_closure']=tri/nb
    out['ccdf']=[[int(q),float((k>=q).mean())] for q in [2,4,8,16,32,64,128,256] ]
    return out, k

if __name__=='__main__':
    N=int(sys.argv[1]); seeds=int(sys.argv[2])
    Pis=[float(x) for x in sys.argv[3].split(',')]
    res={}
    for Pi in Pis:
        rows=[]
        for sd in range(seeds):
            t=time.time(); r=run(N,Pi,seed=sd); o,k=stats(r); o['sec']=time.time()-t; rows.append(o)
            np.save(f'../results/k_N{N}_Pi{Pi:g}_s{sd}.npy',k)
        res[Pi]=rows
        g=np.mean([o['gamma10'] for o in rows]); print(f"Pi={Pi:g} gamma10={g:.2f} gamma6={np.mean([o['gamma6'] for o in rows]):.2f} kmax={np.mean([o['kmax'] for o in rows]):.0f} C={np.mean([o['C'] for o in rows]):.3f} tri={np.mean([o['tri_closure'] for o in rows]):.3f} corrLage={np.mean([o['corr_k_Lage'] for o in rows]):.2f} cohorts={rows[0]['cohort']} Ck={[(a,round(c,3)) for a,b,m,c in rows[0]['Ck']]} t={np.mean([o['sec'] for o in rows]):.0f}s",flush=True)
    json.dump(res,open(f'../results/scan_N{N}.json','w'),indent=1)
