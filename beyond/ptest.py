import time, numpy as np, divisions, waves, periodic
for rule in ('SR','balanced'):
    g=divisions.run(1000,100,rule); L=np.linalg.norm(g['B']-g['A'],axis=1)
    lam=np.ones(1000); lam*=0.2/np.sum(lam*L)
    t0=time.time()
    o=waves.analyse(g,lam)
    sc=periodic.build(g,lam)
    C=periodic.periodic_C(sc); th,n=waves.christoffel_n(C,sc['rho'])
    print(rule,'rho %.3f/%.3f'%(o['rho'],sc['rho']),'kin n_x %.2f n_mean %.2f | per n_x %.2f n_mean %.2f eta %.2f'%(o['n_x'],o['n_mean'],n[0],n.mean(),n.min()/n.max()),'w1 %.3f'%o['omega1'],'%.1fs'%(time.time()-t0),flush=True)
    t0=time.time(); b=periodic.bloch_x(sc,n[0],nk=16)
    for k,w,p in b: print('   k %.2f w %.3f n %.2f pL %.2f'%(k,w,k/w if w>0 else np.nan,p))
    print('  bloch %.1fs'%(time.time()-t0),flush=True)
