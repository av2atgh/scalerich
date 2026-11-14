"""3D plates: T=300 vs T=600, SR and balanced, alpha 0, rho 0.1, 2 seeds -> ../results/tcheck3d.json"""
import json, numpy as np, div3d, mesh3d, shell3d
out=[]
for T in (300,600):
    for rule in ('SR','balanced','regular'):
        for s in (0,1):
            g=div3d.run(T,300+s,rule); X,tri,tp=mesh3d.mesh(g,0.07)
            A=0.5*np.linalg.norm(np.cross(X[tri[:,1]]-X[tri[:,0]],X[tri[:,2]]-X[tri[:,0]]),axis=1)
            for a in (0.0,0.5):
                lam,_=shell3d.thickness(g,tp,A,0.1,a)
                o=shell3d.analyse(X,tri,lam[tp]); o.update(T=T,rule=rule,seed=s,alpha=a)
                out.append(o); json.dump(out,open('../results/tcheck3d.json','w'))
                print(T,rule,s,a,'n %.3f eta %.2f w1 %.3f'%(o['n_mean'],o['eta_n'],o['omega1']),flush=True)
