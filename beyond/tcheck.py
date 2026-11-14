import json, numpy as np, divisions, waves
rows=[]
for T in (1000, 3000):
    for rule in ('SR','balanced','uniform'):
        for s in range(3):
            g=divisions.run(T,200+s,rule); L=np.linalg.norm(g['B']-g['A'],axis=1)
            for a in (0.0,0.5):
                base=np.arange(1,T+1.)**-a
                for r in (0.1,0.2):
                    o=waves.analyse(g,base*r/np.sum(base*L)); o.update(rule=rule,T=T,alpha=a,rho_t=r,seed=s); rows.append(o)
json.dump(rows,open('results/tcheck.json','w'))
import collections; G=collections.defaultdict(list)
for o in rows: G[(o['rule'],o['alpha'],o['rho_t'],o['T'])].append(o)
for k in sorted(G):
    v=G[k]; print(k, 'n_mean %.2f eta %.2f w1/sqrtT %.4f'%(np.mean([x['n_mean'] for x in v]),np.mean([x['eta_n'] for x in v]),np.mean([x['omega1'] for x in v])/np.sqrt(k[3])))
