"""3D sweep: rules x alpha x rho x seeds (T=300, h=0.07) -> ../results/sweep3d_T300.json (appends per run)"""
import json, os, sys, time, numpy as np, div3d, mesh3d, shell3d
T=int(sys.argv[1]) if len(sys.argv)>1 else 300
H=0.07; NS=3; ALPHAS=(0.0,0.5,1.0); RHOS=(0.05,0.1,0.2)
out=f'../results/sweep3d_T{T}.json'
rows=json.load(open(out)) if os.path.exists(out) else []
done={(r['rule'],r['seed'],r['alpha'],r['rho_t']) for r in rows}
t0=time.time()
for rule in div3d.RULES:
    for s in range(NS):
        if all((rule,s,a,r) in done for a in ALPHAS for r in RHOS): continue
        try:
            g=div3d.run(T,100+s,rule); X,tri,tp=mesh3d.mesh(g,H)
        except Exception as e:
            print('FAILED build',rule,s,e,flush=True); continue
        A=0.5*np.linalg.norm(np.cross(X[tri[:,1]]-X[tri[:,0]],X[tri[:,2]]-X[tri[:,0]]),axis=1)
        for a in ALPHAS:
            for r in RHOS:
                if (rule,s,a,r) in done: continue
                lam,ap=shell3d.thickness(g,tp,A,r,a)
                try: o=shell3d.analyse(X,tri,lam[tp])
                except Exception as e: print('FAILED',rule,s,a,r,e,flush=True); continue
                slim=lam/np.sqrt(np.maximum(ap,1e-12))
                o.update(rule=rule,seed=s,alpha=a,rho_t=r,thin_max=float(slim.max()),thin_med=float(np.median(slim)))
                rows.append(o); json.dump(rows,open(out,'w'))
        print(rule,s,round(time.time()-t0),flush=True)
