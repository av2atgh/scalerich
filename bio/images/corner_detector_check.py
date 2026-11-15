"""Calibration of the corner-cut detector: fraction of a cell's neighbour pairs that meet at a
vertex of the cell (expected ~2/(n-1)), Willis plant1, 20 h frame, first 200 cells."""
import sys, glob, numpy as np
sys.path.insert(0,'/Users/avazquez/av2atg/scalerich/bio/images')
import meristem as m
f=[x for x in glob.glob(m.ROOT+'/plant1/segmentation_tiffs/*') if '/20hrs_' in x][0]
lab,h=m.load_proj(f); V,E,vl,walls=m.frame(lab)
nb={}
for (a,b) in walls:
    if a: nb.setdefault(a,set()).add(b); nb.setdefault(b,set()).add(a)
vs=[set(l) for l in vl.values()]; adj=tot=0
for M,ns in list(nb.items())[:200]:
    ns=[x for x in ns if x]
    for i in range(len(ns)):
        for j in range(i+1,len(ns)):
            tot+=1; adj+= any({M,ns[i],ns[j]}<=s for s in vs)
print('neighbour pairs adjacent around M:',round(adj/tot,4),' mean neighbours',round(np.mean([len(v) for v in nb.values()]),3))
