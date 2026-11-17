"""Gaps in dense spin-wave spectra (results/spinwave.json). Unfold by the Weyl count
x = L_tot k(omega)/pi with k^2 = [-(a+b) + sqrt((a-b)^2 + 4 omega^2)]/2; list the largest gaps
in k (a band gap keeps its width in k as T grows) with the IDOS per cell, N/(T+1), at the gap."""
import json, sys, numpy as np
a, b = 1.0, 4.0
def kof(om): return np.sqrt(np.maximum(0, (-(a + b) + np.sqrt((a - b) ** 2 + 4 * om ** 2)) / 2))
d = json.load(open('results/spinwave.json'))
for key, v in d.items():
    om = np.array(v['om']); k = kof(om); T = v['T']
    sel = (k > 0.3) & (k < 4.0); idx = np.where(sel)[0]
    dk = np.diff(k[idx]); x = v['Ltot'] * dk / np.pi
    top = np.argsort(dk)[::-1][:6]
    print('%-22s modes %5d  max gap (mean spacings) %.0f' % (key, sel.sum(), x.max()))
    print('    ' + '  '.join('k %.3f dk %.3f N/(T+1) %.3f' % (k[idx[j]], dk[j], (idx[j] + 1) / (T + 1)) for j in sorted(top)))
