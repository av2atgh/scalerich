"""Spin-wave mode extent vs size: participation fraction of the nev modes nearest omega0
(shift-invert), for omega0 in OMS. argv: names, Ts, seeds, out. -> results/spinwave_scale.json"""
import sys, json, time, numpy as np
import spinwave as sw, supercond as sc

OMS = (5.0, 10.0, 20.0)

def build(nm, T, sd):
    if nm == 'uniform': return sc.qdiv.run(T, ('largest', 'uniform', 'perp'), seed=sd)
    if nm == 'SR': return sc.divisions.run(T, sd, 'SR')
    return sw.CASES[nm](T)

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
                    X, S = sw.network(build(nm, T, sd), T)
                    Ku, Kw, G, info = sw.assemble(X, S)
                    res = dict(rule=nm, T=T, seed=sd, nu=Ku.shape[0], Ltot=info['Ltot'])
                    for w0 in OMS:
                        om, U, Wv = sw.spectrum_window(Ku, Kw, G, w0, 150)
                        pr = sw.density(U, None, om, info, Wv=Wv)
                        res['om%g' % w0] = om.tolist(); res['pr%g' % w0] = pr.tolist()
                except Exception as err:
                    print(key, 'FAILED', repr(err), flush=True); continue
                out[key] = res; json.dump(out, open(fout, 'w'))
                print(key, 'nu', res['nu'], ' '.join('pr%g %.4f' % (w0, np.median(res['pr%g' % w0])) for w0 in OMS),
                      '%.0fs' % (time.time() - t0), flush=True)
