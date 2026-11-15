"""Model counterparts of the Refahi statistics (Tcc = 31 h): host-linked fraction restricted to
hosts born within a film window before the new wall, mean degree vs line age, angle vs age."""
import numpy as np, json
from plantdiv import Tissue
Tcc=31.0; unit=Tcc/np.log(2); film=150.0
out={}
for name,kw in [('relaxed K=20',dict(relax=True,K=20)),('relaxed K=100',dict(relax=True,K=100)),('unrelaxed',dict(relax=False))]:
    L=[];KA=[];ANG=[]
    for sd in range(2):
        T=Tissue(seed=sd,rule='errera',t0=0.7,tau=36/unit,**kw).grow(2000)
        b=np.array(T.birth)*unit; tend=T.t*unit; t0=tend-film
        adj={}
        for lid,hs in enumerate(T.hosts):
            for h in hs:
                if h>=0: adj.setdefault(lid,set()).add(h); adj.setdefault(h,set()).add(lid)
        for lid,(h1,h2) in enumerate(T.hosts):
            if h1<0 or h2<0 or b[lid]<t0: continue
            if b[h1]<t0 or b[h2]<t0: continue
            L.append(h2 in T.hosts[h1] or h1 in T.hosts[h2])
        k=np.array([len(adj.get(i,())) for i in range(len(b))])
        for i in range(len(b)):
            if b[i]>=t0: KA.append((tend-b[i],k[i]))
        ANG+=[(a*unit,g) for a,g,_ in T.junction_angles()]
    L=np.array(L); KA=np.array(KA); A=np.array(ANG)
    o=dict(linked=float(L.mean()),n=len(L))
    o['k_vs_age']={f'{lo}-{hi}':round(float(KA[(KA[:,0]>=lo)&(KA[:,0]<hi),1].mean()),2) for lo,hi in [(0,20),(20,40),(40,70),(70,110),(110,150)]}
    o['angle']={f'{lo}-{hi}':round(float(np.median(A[(A[:,0]>=lo)&(A[:,0]<hi),1])),1) for lo,hi in [(0,8),(8,16),(16,24),(24,36),(36,48),(48,72),(72,120)]}
    out[name]=o; print(name,o,flush=True)
json.dump(out,open('../results/refahi_model.json','w'),indent=1)
