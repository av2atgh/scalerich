"""Summarise results/flow.json -> results/flow_summary.json (+ printed table).
psi = sum dP / sum aff over chords in a tau window (pooled over seeds),
F = sum g dV^2 / sum aff (field transmission), gR median and mean of 1/(1+gR)
weighted by g dV^2, Phi(tau) = sigma(tau)/sigma_aff(tau)."""
import json
import numpy as np
import srtree

d = json.load(open('results/flow.json'))
T = d['T']; NS = d['NS']; taus = d['taus']
gs = [srtree.run(T, 900 + s) for s in range(NS)]
out = {}
for a, runs in d['data'].items():
    a_ = float(a)
    lam = srtree.thickness(T, 1.0, a_)
    rows = {}
    for s, recs in enumerate(runs):
        dd = gs[s]['B'] - gs[s]['A']; l = np.linalg.norm(dd, axis=1)
        aff = np.cumsum(lam * dd[:, 0] ** 2 / l)
        for r in recs:
            tau = r['tau']
            e = rows.setdefault(tau, dict(dP=0., aff=0., F=0., th=0., gR=[], Phi=[]))
            e['Phi'].append(r['sigma'] / aff[tau - 1])
            for x in r['rows']:
                e['dP'] += x['dP']; e['aff'] += x['aff']
                if np.isfinite(x['R']):
                    gd = x['gc'] * x['dV'] ** 2
                    e['F'] += gd; e['gR'].append(x['gc'] * x['R'])
    # pool tau into 4 log windows
    ts = sorted(rows)
    win = np.array_split(np.array(ts), 5)
    res = []
    for w in win:
        dP = sum(rows[t]['dP'] for t in w); af = sum(rows[t]['aff'] for t in w)
        F = sum(rows[t]['F'] for t in w)
        gR = np.concatenate([rows[t]['gR'] for t in w])
        Phi = np.mean([rows[t]['Phi'] for t in w])
        res.append(dict(tau=[int(w[0]), int(w[-1])], psi=dP / af, F=F / af, thev=dP / F,
                        gR_median=float(np.median(gR)), Phi=float(Phi), n=int(len(gR))))
    out[a] = res
    print('alpha', a)
    for r in res:
        print('  tau %5d-%5d  psi %.3f  F %.3f  thevenin %.3f  median gR %.3f  Phi(tau) %.3f  (n=%d)' % (
            r['tau'][0], r['tau'][1], r['psi'], r['F'], r['thev'], r['gR_median'], r['Phi'], r['n']))
json.dump(out, open('results/flow_summary.json', 'w'), indent=1)

# window increments psi = d sigma / d sigma_aff between consecutive checkpoints (exact, all chords)
print('\nwindow increments psi = (sigma(t2)-sigma(t1))/(aff(t2)-aff(t1)), pooled over seeds')
inc = {}
for a, runs in d['data'].items():
    lam = srtree.thickness(T, 1.0, float(a))
    num = np.zeros(len(taus) - 1); den = np.zeros(len(taus) - 1)
    for s, recs in enumerate(runs):
        dd = gs[s]['B'] - gs[s]['A']; l = np.linalg.norm(dd, axis=1)
        aff = np.cumsum(lam * dd[:, 0] ** 2 / l)
        sig = {r['tau']: r['sigma'] for r in recs}
        for k in range(len(taus) - 1):
            t1, t2 = taus[k], taus[k + 1]
            if t1 in sig and t2 in sig:
                num[k] += sig[t2] - sig[t1]; den[k] += aff[t2 - 1] - aff[t1 - 1]
    inc[a] = (num / den).round(3).tolist()
    print(a, inc[a])
out['_window_psi'] = dict(taus=taus, psi=inc)
json.dump(out, open('results/flow_summary.json', 'w'), indent=1)
