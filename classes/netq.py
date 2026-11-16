"""Network features (beyond/network/netstats.features) for the split-law family, T = 2e4,
plus the fraction of junction points shared by two stems (X-crossings)."""
import sys, json, numpy as np
sys.path[:0] = ['../beyond', '../beyond/network', '../rg', '../lng']
import qdiv, divisions, netstats
from diffraction import RULES

def xfrac(g):
    P = np.round(np.vstack([g['A'], g['B']]), 10)
    _, c = np.unique(P, axis=0, return_counts=True)
    return float((c[c > 1]).sum() / len(P))

if __name__ == '__main__':
    T = 20000; out = []
    names = ['regular', 'golden', 'plastic', 'r0.3', 'golden-side', 'largest-uniform', 'SR', 'balanced']
    for name in names:
        for seed in (0, 1):
            g = qdiv.run(T, RULES[name], seed=seed) if name in RULES else divisions.run(T, seed=seed, rule=name)
            f = netstats.features(g, paths=False); f['xfrac'] = xfrac(g); f['rule'] = name; f['seed'] = seed
            out.append(f)
            print(name, seed, 'C=%.3f k2/k=%.1f kmax=%d K2P(16)=%.2f K2P(64)=%.2f depthvar=%.2f spanmax=%.2f spancv=%.2f x=%.3f tri? ht=%s' % (
                f['Cmean'], f['k2_over_k'], f['kmax'], 256 * f['ccdf'][16], 4096 * f['ccdf'][64], f['depth_var'],
                f['span_max'], f['span_cv'], f['xfrac'], np.round(f['ht_moments'], 3)), flush=True)
            if name in RULES and RULES[name][1] in (0.5,) or (name in ('golden', 'plastic') ): break
    json.dump(out, open('results/net.json', 'w'), indent=1)
