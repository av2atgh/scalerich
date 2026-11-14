"""3D open-cell skeleton sweep: rules x alpha x rho x seeds (T=300) -> ../results/skel3d_T300.json"""
import json, time, numpy as np, div3d, skeleton3d
T=300; out=[]; t0=time.time()
for rule in div3d.RULES:
    for s in range(3):
        g=div3d.run(T,100+s,rule); X,seg,own=skeleton3d.skeleton(g)
        L=np.linalg.norm(X[seg[:,1]]-X[seg[:,0]],axis=1)
        for a in (0.0,0.5,1.0):
            base=np.arange(1,T+1.)**-a
            for r in (0.01,0.03,0.1):
                lam=base*np.sqrt(r/np.sum(base[own]**2*L)); ls=lam[own]
                try: o=skeleton3d.analyse(X,seg,ls)
                except Exception as e: print('FAILED',rule,s,a,r,e,flush=True); continue
                o.update(rule=rule,seed=s,alpha=a,rho_t=r,slender_min=float(np.min(L/ls)),
                         frac_stubby=float(np.mean(L<ls)))
                out.append(o); json.dump(out,open('../results/skel3d_T300.json','w'))
        print(rule,s,round(time.time()-t0),flush=True)
