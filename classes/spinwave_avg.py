"""Frequency-averaged spin-wave mode extent: for each omega0 in OMS, the 60 modes nearest omega0
(shift-invert); per window the median participation fraction; summary = geometric mean over
windows. argv: names, Ts, seeds, out -> results/spinwave_avg.json"""
import sys, json, time, numpy as np
import spinwave as sw, spinwave_scale as ss

OMS = (4.0, 4.5, 5.0, 5.5, 6.0, 6.5, 7.0, 8.0)

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
                    med = []
                    for w0 in OMS:
                        om, U, Wv = sw.spectrum_window(Ku, Kw, G, w0, 60)
                        med.append(float(np.median(sw.density(U, None, om, info, Wv=Wv))))
                except Exception as err:
                    print(key, 'FAILED', repr(err), flush=True); continue
                gm = float(np.exp(np.mean(np.log(med))))
                out[key] = dict(rule=nm, T=T, seed=sd, Ltot=info['Ltot'], oms=list(OMS), med=med, gm=gm)
                json.dump(out, open(fout, 'w'))
                print(key, 'gm %.5f' % gm, ' '.join('%.4f' % m for m in med), '%.0fs' % (time.time() - t0), flush=True)
