"""Summaries of the elastic runs.
(1) E2 collapse test: window increment psi_el = dE/daff (aff = sum lam l cos^4,
    x/y averaged) against the newest-chord slenderness s = c_l tau^-1/2 / lam_tau
    at the window's geometric midpoint, all (alpha, rho_end) pooled.
(2) T-route energy growth x_T = dlnE/dlnT over the last 6 checkpoints, against
    the paper's 1 - alpha/2 (+p) and against the density route m_T.
Writes results/elastic_summary.json
"""
import json
import numpy as np

d = json.load(open('results/elastic_e2.json'))
taus = np.array(d['taus'])
out = {'collapse': [], 'xT': {}}
print('alpha rho_end | s_mid -> psi_el')
for key, recs in d.items():
    if key == 'taus':
        continue
    a, re = map(float, key.split('_'))
    seeds = sorted({r['seed'] for r in recs})
    E = np.array([[r['E'] for r in recs if r['seed'] == s] for s in seeds])
    A = np.array([[r['aff'] for r in recs if r['seed'] == s] for s in seeds])
    R = np.array([[r['rho'] for r in recs if r['seed'] == s] for s in seeds])
    S = np.array([[r['s_new'] for r in recs if r['seed'] == s] for s in seeds]).mean(0)
    psi = np.diff(E, axis=1).sum(0) / np.diff(A, axis=1).sum(0)
    smid = np.sqrt(S[1:] * S[:-1])
    for s_, p_ in zip(smid, psi):
        out['collapse'].append((a, re, float(s_), float(p_)))
    lE = np.log(E.mean(0)); lR = np.log(R.mean(0)); lT = np.log(taus)
    xT = float(np.polyfit(lT[-6:], lE[-6:], 1)[0])
    mT = float(np.polyfit(lR[-6:], lE[-6:], 1)[0])
    out['xT'][key] = dict(xT=xT, paper_1_minus_a2=1 - a / 2, mT=mT,
                          paper_m_p0=(1 - a / 2) / (0.5 - a) if a < 0.5 else None,
                          E_end=float(E.mean(0)[-1]), rho_end=float(R.mean(0)[-1]),
                          E_over_rho_end=float(E.mean(0)[-1] / R.mean(0)[-1]))
    print(key, ' '.join('%.2g:%.3f' % (s_, p_) for s_, p_ in zip(smid[4:], psi[4:])))
print('\nT-route: x_T = dlnE/dlnT and m_T = dlnE/dlnrho over tau in [%d, %d]' % (taus[-6], taus[-1]))
for k, v in out['xT'].items():
    print(k, {q: (round(w, 3) if isinstance(w, float) else w) for q, w in v.items()})
# collapse quality: bin pooled points by log s
C = np.array(out['collapse'])
C = C[C[:, 2] < 200]
bins = np.logspace(np.log10(C[:, 2].min()), np.log10(C[:, 2].max()), 8)
print('\npooled psi_el by slenderness bin (mean, sd across (alpha,rho) points)')
tab = []
for b0, b1 in zip(bins[:-1], bins[1:]):
    m = (C[:, 2] >= b0) & (C[:, 2] < b1)
    if m.sum():
        tab.append((float(np.sqrt(b0 * b1)), float(C[m, 3].mean()), float(C[m, 3].std()), int(m.sum()),
                    sorted({float(x) for x in C[m, 0]})))
        print('s ~ %.3g  psi %.3f +- %.3f  (n=%d, alphas %s)' % tab[-1])
out['collapse_table'] = tab
json.dump(out, open('results/elastic_summary.json', 'w'), indent=1)
