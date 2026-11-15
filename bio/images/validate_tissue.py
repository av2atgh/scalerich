import sys, warnings, numpy as np; warnings.filterwarnings('ignore')
sys.path.insert(0,'.'); sys.path.insert(0,'../../beyond')
import divisions
from skimage.draw import polygon
from planar import label_graph, classify, lines, dual_stats, kinds, corners_per_cell
M=int(sys.argv[2]) if len(sys.argv)>2 else 2048; T=int(sys.argv[1]) if len(sys.argv)>1 else 1000
for rule in ['SR','balanced']:
    r=divisions.run(T,seed=1,rule=rule)
    Lab=np.zeros((M,M),np.int32)
    for i,p in enumerate(r['leaf_polys']):
        p=np.asarray(p)*(M-1); rr,cc=polygon(p[:,1],p[:,0],(M,M)); Lab[rr,cc]=i+1
    pad=np.zeros((M+8,M+8),np.int32); pad[4:-4,4:-4]=Lab; Lab=pad   # background ring = 0
    V,E,vl=label_graph(Lab)
    recs=classify(V,[e[:3] for e in E],r1=3,r2=12)
    c,tf=kinds(recs)
    length=np.array([np.hypot(*np.diff(e[2],axis=0).T).sum() for e in E])
    li=lines(len(E),recs,length=length)
    st,k=dual_stats(li['nL'],li['D'],li['L'])
    cp=corners_per_cell(V,E,recs,vl); cn=np.array([x[0] for x in cp.values()]); nv=np.array([x[1] for x in cp.values()])
    kt=np.bincount(np.r_[r['hostA'],r['hostB']][np.r_[r['hostA'],r['hostB']]>=0],minlength=T)+ (np.c_[r['hostA'],r['hostB']]>=0).sum(1)
    print(rule, 'kinds',{a:c[a] for a in ('T','Y','X') if a in c},f"Tfrac={tf:.3f} corners={cn.mean():.2f} vertices={nv.mean():.2f} lines={li['nL']} (true {T})")
    print('   image dual', {a:b for a,b in st.items() if a!='Ck'}, st['Ck'])
    print('   truth      ', {f'P(k>={q})':round(float((kt>=q).mean()),4) for q in (2,4,8,16,32)}, 'kmax',kt.max())
