"""Pipeline on BioSR F-actin linear SIM (COS-7, Lifeact; 31.3 nm/px). Run: ../.venv/bin/python -I run_biosr.py"""
import sys, glob, json, warnings, numpy as np
warnings.filterwarnings('ignore')
sys.path.insert(0,'/Users/avazquez/av2atg/scalerich/bio/images')
import mrcfile
from scipy import ndimage as ndi
from skimage.filters import gaussian, threshold_otsu
from pipeline import analyse, summary
files=sorted(glob.glob('/Users/avazquez/av2atg/scalerich/bio/data/raw/biosr/F-actin/Cell_*/SIM_gt.mrc'))
out={}; pool_k=[]; pool_len=[]; pool_W=[]; pool_I=[]
for f in files:
    with mrcfile.open(f, permissive=True) as m:
        img=np.asarray(m.data,float)
    if img.ndim==3: img=img.max(0)
    sm=gaussian(img,8); mask=sm>threshold_otsu(sm)*0.5
    mask=ndi.binary_fill_holes(ndi.binary_opening(mask,iterations=3))
    try:
        res=analyse(img,mask=mask,sigmas=(1,2,3),lmin=8,r1=6,r2=18,merge_len=8)
    except Exception as e:
        print(f,'fail',e); continue
    s=summary(res); name=f.split('/')[-2]; out[name]=s
    pool_k+=list(res['k']); pool_len+=list(res['flen']); pool_W+=list(res['W']); pool_I+=list(res['I'])
    print(name, s['kinds'].get('T',0), s['kinds'].get('Y',0), s['kinds'].get('X',0), f"Tfrac={s['T_frac']:.2f} X/J={s['X_per_junction']:.3f} nF={s['nF']} kmax={s['kmax']}",flush=True)
k=np.array(pool_k); L=np.array(pool_len); W=np.array(pool_W); I=np.array(pool_I)
np.savez('/Users/avazquez/av2atg/scalerich/bio/results/biosr_pool.npz',k=k,L=L,W=W,I=I)
tot={}
for v in out.values():
    for a,b in v['kinds'].items(): tot[a]=tot.get(a,0)+b
x=k[k>=4].astype(float)
agg=dict(cells=len(out),kinds=tot,T_frac=tot.get('T',0)/(tot.get('T',0)+tot.get('Y',0)),
  X_per_junction=tot.get('X',0)/sum(tot.values()),nF=int(len(k)),
  ccdf={q:float((k>=q).mean()) for q in (2,4,8,16,32)},gamma4=float(1+len(x)/np.sum(np.log(x/3.5))),
  corr_k_len=float(np.corrcoef(k,L)[0,1]),corr_k_W=float(np.corrcoef(k,W)[0,1]),corr_k_I=float(np.corrcoef(k,I)[0,1]))
# partial: k vs width at fixed length (residuals of log-log)
s=(k>0)&(L>0)
rk=np.log(k[s])-np.polyval(np.polyfit(np.log(L[s]),np.log(k[s]),1),np.log(L[s]))
rw=W[s]-np.polyval(np.polyfit(np.log(L[s]),W[s],1),np.log(L[s]))
agg['partial_k_W_given_L']=float(np.corrcoef(rk,rw)[0,1])
print(json.dumps(agg,indent=1))
json.dump(dict(cells=out,agg=agg),open('/Users/avazquez/av2atg/scalerich/bio/results/biosr_summary.json','w'),indent=1,default=str)
