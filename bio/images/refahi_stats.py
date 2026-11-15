import json, glob, numpy as np
rows=[]; ang=[]; linked=[]; K=[]; rates=[]
for f in sorted(glob.glob('/Users/avazquez/av2atg/scalerich/bio/results/refahi_FM*.json')):
    d=json.load(open(f)); birth={int(k):v for k,v in d['birth'].items()}; H={int(k):v for k,v in d['hostrel'].items()}
    tend=d['times'][-1]
    for ta,tb,fr,l2 in d['fracs']: rates.append((tb-ta,l2))
    ang+= [tuple(x) for x in d['relax']]
    for new,hs in H.items():
        if len(hs)!=2 or None in hs: continue
        h1,h2=hs; b=birth[new]
        if not all(isinstance(birth.get(h),int) and birth[h]<b for h in hs): continue
        if not all(h in H and len(H[h])==2 and None not in H[h] for h in hs): continue
        linked.append((h2 in H[h1]) or (h1 in H[h2]))
    # degree of film-born lines with both own hosts known
    deg={}
    for new,hs in H.items():
        for h in hs:
            if h is not None: deg[h]=deg.get(h,0)+1
    for new,hs in H.items():
        if len(hs)==2 and None not in hs and isinstance(birth[new],int):
            K.append((tend-birth[new], deg.get(new,0)+2, d['fm']))
    rows.append((d['fm'],len(H)))
r=np.array(rates); Tcc=r[:,0].sum()/r[:,1].sum()
print('per-cell doubling time from log2(daughters) per interval: Tcc =',round(Tcc,1),'h')
a=np.array(ang)
for lo,hi in [(0,8),(8,16),(16,24),(24,36),(36,48),(48,72),(72,120)]:
    s=(a[:,0]>=lo)&(a[:,0]<hi)
    if s.sum(): print(f'age [{lo},{hi}) h n {s.sum()} through-angle median {np.median(a[s,1]):.1f}')
L=np.array(linked); print('host-linked: n',len(L),'fraction',L.mean().round(3) if len(L) else None,'+-',np.sqrt(L.mean()*(1-L.mean())/len(L)).round(3) if len(L) else None)
K=np.array([(x[0],x[1]) for x in K],float)
for lo,hi in [(0,20),(20,40),(40,70),(70,110),(110,200)]:
    s=(K[:,0]>=lo)&(K[:,0]<hi)
    if s.sum(): print(f'line age [{lo},{hi}) h n {s.sum()} mean k {K[s,1].mean():.2f} ccdf',{q:round(float((K[s,1]>=q).mean()),3) for q in (3,4,6,8)},'max',int(K[s,1].max()))
json.dump(dict(Tcc=Tcc,linked=[int(x) for x in L],K=K.tolist(),ang=a.tolist()),open('/Users/avazquez/av2atg/scalerich/bio/results/refahi_stats.json','w'))
