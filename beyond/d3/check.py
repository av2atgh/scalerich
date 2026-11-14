"""checks: affine bound, rigid-body energy, mesh refinement, cubic-grid reference"""
import time, numpy as np, div3d, mesh3d, shell3d
def tri_area(X,tri): return 0.5*np.linalg.norm(np.cross(X[tri[:,1]]-X[tri[:,0]],X[tri[:,2]]-X[tri[:,0]]),axis=1)
for rule in ('SR','balanced'):
    g=div3d.run(200,1,rule)
    for h in (0.07,0.05,0.035):
        t0=time.time(); X,tri,tp=mesh3d.mesh(g,h); A=tri_area(X,tri)
        lam,_=shell3d.thickness(g,tp,A,0.1,0.0)
        K,_=shell3d.assemble(X,tri,lam[tp])
        # rigid rotation energy
        w=np.array([0.3,-0.2,0.5]); u=np.cross(w,X-0.5).ravel()
        o=shell3d.analyse(X,tri,lam[tp])
        print(rule,h,'N',len(X),'rigid E %.1e'%(u@(K@u)),'n_x %.3f n_mean %.3f eta %.3f n_aff %.3f w1 %.3f'%(o['n_x'],o['n_mean'],o['eta_n'],o['n_aff_mean'],o['omega1']),'%.1fs'%(time.time()-t0),flush=True)
