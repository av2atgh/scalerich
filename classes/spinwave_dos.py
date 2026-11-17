"""Coarse-grained spin-wave DOS by Sylvester inertia on a k grid (k in KGRID, omega(k) from the wire
dispersion). Bin counts divided by the Weyl count L_tot dk/pi; roughness = rms deviation from 1.
argv: names, Ts, seeds, out -> results/spinwave_dos.json"""
import sys, json, time, numpy as np
import spinwave as sw, spinwave_scale as ss, spinwave_count as sc_
a, b = 1.0, 4.0
KGRID = np.round(np.arange(1.0, 3.5001, 0.05), 4)
def omof(k): return np.sqrt((k ** 2 + a) * (k ** 2 + b))
if __name__ == '__main__':
    names = sys.argv[1].split(','); Ts = [int(x) for x in sys.argv[2].split(',')]
    seeds = [int(x) for x in sys.argv[3].split(',')]; fout = sys.argv[4]
    for nm in names:
        for T in Ts:
            for sd in (seeds if nm in ('uniform', 'SR') else [0]):
                try: out = json.load(open(fout))
                except Exception: out = {}
                key = '%s_T%d_s%d' % (nm, T, sd)
                if key in out: continue
                t0 = time.time()
                try:
                    X, S = sw.network(ss.build(nm, T, sd), T); Ku, Kw, G, info = sw.assemble(X, S)
                    N = sc_.counter(Ku, Kw, G); cnt = np.array([N(omof(k)) for k in KGRID])
                except Exception as err:
                    print(key, 'FAILED', repr(err), flush=True); continue
                ratio = np.diff(cnt) / (info['Ltot'] * np.diff(KGRID) / np.pi)
                out[key] = dict(rule=nm, T=T, seed=sd, k=KGRID.tolist(), N=cnt.tolist(), Ltot=info['Ltot'],
                                ratio=ratio.tolist(), rough=float(np.sqrt(np.mean((ratio - 1) ** 2))))
                json.dump(out, open(fout, 'w'))
                print(key, 'rough %.3f' % out[key]['rough'], '%.0fs' % (time.time() - t0), flush=True)
