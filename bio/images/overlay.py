import sys, warnings, numpy as np; warnings.filterwarnings('ignore')
sys.path.insert(0,'/Users/avazquez/av2atg/scalerich/bio/images')
import mrcfile, matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
from scipy import ndimage as ndi
from skimage.filters import gaussian, threshold_otsu
from pipeline import analyse
c=sys.argv[1]
img=np.asarray(mrcfile.open(f'/Users/avazquez/av2atg/scalerich/bio/data/raw/biosr/F-actin/{c}/SIM_gt.mrc',permissive=True).data,float)
if img.ndim==3: img=img.max(0)
sm=gaussian(img,8); mask=ndi.binary_fill_holes(ndi.binary_opening(sm>threshold_otsu(sm)*0.5,iterations=3))
r=analyse(img,mask=mask)
fig,ax=plt.subplots(1,2,figsize=(16,8))
ax[0].imshow(img,cmap='gray',vmax=np.percentile(img,99.5)); ax[0].set_title(c)
ax[1].imshow(img,cmap='gray',vmax=np.percentile(img,99.5)); ax[1].imshow(np.ma.masked_where(~r['S'],r['S']),cmap='autumn',alpha=0.9)
col={'T':'c','Y':'m','X':'y'}
for rec in r['recs']:
    p=r['G']['Jc'][rec['j']]; ax[1].plot(p[1],p[0],'o',ms=3,color=col.get(rec['kind'],'w'))
for a in ax: a.axis('off')
plt.tight_layout(); plt.savefig(f'/Users/avazquez/av2atg/scalerich/bio/results/overlay_{c}.png',dpi=70)
