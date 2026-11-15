"""Crossover rank n*: cohorts born before n* have mean degree below SR's.
n* = rank where cohort k / SR cohort k crosses 1/2 (and 0.8). Predicted n* ~ Pi^(2/3)."""
import numpy as np
N=20000; edges=np.unique(np.logspace(0,np.log10(N),30).astype(int))
def cohort(Pi):
    ks=np.mean([np.load(f'../results/k_N{N}_Pi{Pi}_s{s}.npy') for s in range(3)],axis=0)
    return np.array([ks[lo:hi].mean() for lo,hi in zip(edges[:-1],edges[1:])])
c=np.sqrt(edges[:-1]*edges[1:]); ref=cohort('0.001')
print("Pi     n*(1/2)  n*(0.8)  n*(1/2)/Pi^(2/3)")
for Pi in ['30','100','300','1000','10000','1e+06']:
    r=cohort(Pi)/ref
    def cross(th):
        j=np.where(r>=th)[0]
        return c[j[0]] if len(j) else np.nan
    print(f"{Pi:>6} {cross(0.5):8.0f} {cross(0.8):8.0f} {cross(0.5)/float(Pi)**(2/3):8.3f}   ratios:{np.round(r[:12],2)}")
