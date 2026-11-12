import numpy as np, sys
rng=np.random.default_rng(int(sys.argv[2])); T=int(sys.argv[1])
A=np.zeros((T,2)); B=np.zeros((T,2)); born=np.zeros(T,int)
edges=[]; adj=[set() for _ in range(T)]; tri=0; nonb=0
for t in range(T):
    p=rng.random(2); th=rng.random()*np.pi; d=np.array([np.cos(th),np.sin(th)])
    # boundary
    with np.errstate(divide='ignore'):
        sx=np.array([(0-p[0])/d[0],(1-p[0])/d[0]]); sy=np.array([(0-p[1])/d[1],(1-p[1])/d[1]])
    ss=np.concatenate([sx,sy]); smax=ss[ss>0].min(); smin=ss[ss<0].max(); hp=hm=-1
    if t:
        a=A[:t]; e=B[:t]-a; w=a-p
        den=d[0]*e[:,1]-d[1]*e[:,0]
        s=(w[:,0]*e[:,1]-w[:,1]*e[:,0])/den; u=(w[:,0]*d[1]-w[:,1]*d[0])/den
        ok=(u>=0)&(u<=1)
        i=np.where(ok&(s>0)&(s<smax))[0]
        if len(i): j=i[np.argmin(s[i])]; smax=s[j]; hp=j
        i=np.where(ok&(s<0)&(s>smin))[0]
        if len(i): j=i[np.argmax(s[i])]; smin=s[j]; hm=j
    A[t]=p+smin*d; B[t]=p+smax*d
    for h in (hp,hm):
        if h>=0: adj[t].add(h); adj[h].add(t)
    if hp>=0 and hm>=0:
        nonb+=1; tri+= hm in adj[hp]
L=np.linalg.norm(B-A,axis=1); k=np.array([len(x) for x in adj])
np.savez(sys.argv[3],L=L,k=k)
print("tri",tri,"nonb",nonb);print("T",T,"<k>",k.mean(),"frac new nodes closing triangle",tri/nonb)
m=np.arange(T)>T//10
print("corr log k vs log L (k>0, old):",np.corrcoef(np.log(k[m&(k>0)]),np.log(L[m&(k>0)]))[0,1])
# slope of k on L for interior ligaments
sl=np.polyfit(np.log(L[m&(k>2)]),np.log(k[m&(k>2)]),1)[0]; print("slope log k vs log L",sl)
kk=np.sort(k)[::-1]; r=np.arange(1,T+1)
sel=(kk>=8)&(kk<=kk[20]); g=1-np.polyfit(np.log(kk[sel]),np.log(r[sel]),1)[0]; print("gamma from rank plot",g)
# mean degree vs birth time
for lo,hi in [(10,30),(100,300),(1000,3000)]: print("birth",lo,hi,"mean k",k[lo:hi].mean(),"expected ratio (t/i)^0.5")
# clustering of dual
C=[]
for i in range(T):
    n=list(adj[i]); ki=len(n)
    if ki<2: continue
    c=sum(1 for x in range(ki) for y in range(x+1,ki) if n[y] in adj[n[x]]); C.append(c/(ki*(ki-1)/2))
print("mean clustering",np.mean(C))
i=np.arange(T); y=k-2.0; x=L*(np.sqrt(T)-np.sqrt(i+1)); s=(i>50)&(i<T*0.8)
print("corr (k-2) vs L*(sqrt T - sqrt i):",np.corrcoef(y[s],x[s])[0,1], " ratio median",np.median(y[s&(x>0.5)]/x[s&(x>0.5)]))
Ci=np.zeros(T); knn=np.zeros(T)
for i in range(T):
    n=list(adj[i]); ki=len(n)
    if ki: knn[i]=np.mean([len(adj[j]) for j in n])
    if ki>=2: Ci[i]=sum(1 for x in range(ki) for y in range(x+1,ki) if n[y] in adj[n[x]])/(ki*(ki-1)/2)
for lo,hi in [(2,3),(3,5),(5,10),(10,20),(20,40),(40,80),(80,1000)]:
    s=(k>=lo)&(k<hi); print(f"k[{lo},{hi}) n={s.sum()} C={Ci[s].mean():.3f} knn={knn[s].mean():.2f}")
print("c2=<k(k-1)>/<k> =",(k*(k-1)).mean()/k.mean())
