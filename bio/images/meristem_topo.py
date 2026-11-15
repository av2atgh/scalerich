"""Meristem topology: P(n) of interior low-slope L1 cells, n_d1 + n_d2 - n_M at divisions,
division fraction per 4 h (cell cycle), host lines already linked (dual clustering proxy)."""
import sys, glob, re, json, numpy as np
sys.path.insert(0,'/Users/avazquez/av2atg/scalerich/bio/images')
import meristem as m
from scipy import ndimage as ndi
out={}
for plant in sys.argv[1:]:
    files=glob.glob(f'{m.ROOT}/{plant}/segmentation_tiffs/*hrs_*.tif.gz')
    times=sorted({int(re.match(r'(\d+)hrs',f.split('/')[-1]).group(1)) for f in files})
    nb={}; inner={}
    for t in times:
        f=[x for x in files if x.split('/')[-1].startswith(f'{t}hrs_')][0]
        lab,h=m.load_proj(f)
        # neighbours from 4-connectivity contacts
        a=np.r_[lab[:,:-1].ravel(),lab[:-1,:].ravel()]; b=np.r_[lab[:,1:].ravel(),lab[1:,:].ravel()]
        s=a!=b; pairs=set(zip(np.minimum(a[s],b[s]),np.maximum(a[s],b[s])))
        N={}
        for x,y in pairs: N.setdefault(x,set()).add(y); N.setdefault(y,set()).add(x)
        # keep cells away from rim: not touching 0 and no neighbour touching 0, low slope
        gy,gx=np.gradient(ndi.gaussian_filter(np.nan_to_num(h,nan=0.0),6)); sl=np.hypot(gx,gy)
        cm=ndi.center_of_mass(np.ones_like(lab),lab,list(N.keys()))
        ok={}
        for c,(cy,cx) in zip(N.keys(),cm):
            if c==0: continue
            ok[c]= (0 not in N[c]) and all(0 not in N.get(x,()) for x in N[c]) and sl[int(cy),int(cx)]<0.7
        nb[t]=N; inner[t]=ok
    Pn=[len(nb[t][c]) for t in times for c,ok in inner[t].items() if ok]
    dsum=[]; frac=[]
    for ta,tb in zip(times[:-1],times[1:]):
        if tb-ta!=4: continue
        try: mom,kids=m.tracking(plant,ta,tb)
        except FileNotFoundError: continue
        frac.append(np.mean([len(v)==2 for v in kids.values()]))
        for M,ds in kids.items():
            if len(ds)==2 and inner[ta].get(M) and all(inner[tb].get(d) for d in ds):
                # neighbour divisions landing on M's walls in the same interval would add sides; count them
                dsum.append(len(nb[tb][ds[0]])+len(nb[tb][ds[1]])-len(nb[ta][M]))
    Pn=np.array(Pn); dsum=np.array(dsum)
    out[plant]=dict(Pn={int(n):float((Pn==n).mean()) for n in range(3,11)},nPn=len(Pn),mean_n=float(Pn.mean()),
        dsum=np.bincount(dsum-dsum.min()).tolist(),dsum_min=int(dsum.min()),frac_div_4h=float(np.mean(frac)))
    print(plant,out[plant],flush=True)
json.dump(out,open('/Users/avazquez/av2atg/scalerich/bio/results/meristem_topo.json','w'),indent=1)
