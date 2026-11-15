import numpy as np, json
from plantdiv import Tissue
unit=36/np.log(2); data={4:.04,5:.27,6:.41,7:.23,8:.045,9:.002}
dage=np.array([0,8,16,24,36])+2.0; dang=np.array([143.3,138.4,133.6,131.0,128.0])
for K in (20,100,500,2000):
    Pn=[];ja=[];tri=nb=0
    for sd in range(2):
        T=Tissue(seed=sd,rule='errera',relax=True,t0=0.7,tau=36/unit,K=K).grow(1000)
        Pn+=list(T.sides()); ja+=[(a*unit,g) for a,g,_ in T.junction_angles()]
        adj=[set() for _ in T.birth]
        for lid,hs in enumerate(T.hosts):
            a,b=hs
            if a>=0 and b>=0: nb+=1; tri+= b in adj[a]
            for h in hs:
                if h>=0: adj[lid].add(h); adj[h].add(lid)
    Pn=np.array(Pn); ja=np.array(ja)
    tv=0.5*sum(abs((Pn==n).mean()-data.get(n,0)) for n in range(3,12))
    cur=[np.median(ja[np.abs(ja[:,0]-a)<2,1]) for a in dage]
    print(f"K={K} TV(P(n))={tv:.3f} P(n)",{n:round(float((Pn==n).mean()),3) for n in range(4,10)},"angles",np.round(cur,1),"data",dang,f"rms {np.sqrt(np.mean((np.array(cur)-dang)**2)):.1f} host-linked frac {tri/nb:.3f}",flush=True)
