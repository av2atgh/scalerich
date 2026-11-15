"""Causal test: does the index of the balanced divisions come from quasi-crossings?
For every quasi-crossing (two stems on opposite sides of a host, offset < 0.1 of the local span,
angle < 15 deg), one partner that hosts no junction itself is moved along the host by a fraction
`shift` of the distance to its nearest neighbouring junction, away from its partner. Lengths change
by < 1%; topology is unchanged. Control: the same number of stems, chosen among those NOT in a
quasi-crossing, moved the same way. Index from the framed periodic tiling, alpha = 0, rho = 0.2.
-> results/breakx.json"""
import json, os, sys
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); B = os.path.join(HERE, '..')
sys.path[:0] = [B, os.path.join(B, '..', 'rg')]
import divisions, periodic, waves


def ends_on_hosts(g):
    T = g['T']; v = g['B'] - g['A']; L = np.linalg.norm(v, axis=1); u = v / L[:, None]
    ends = {}
    for e, (P, hosts) in enumerate(((g['A'], g['hostA']), (g['B'], g['hostB']))):
        for t in range(T):
            h = hosts[t]
            if h < 0: continue
            d = u[t] if e == 0 else -u[t]
            s = (P[t] - g['A'][h]) @ u[h]
            side = np.sign(d @ np.array([-u[h, 1], u[h, 0]]))
            ends.setdefault(int(h), []).append(dict(s=s, side=side, d=d, t=t, e=e))
    return ends, u, L


def classify(g, delta=0.1, phi=np.deg2rad(15)):
    ends, u, L = ends_on_hosts(g)
    hosted = np.zeros(g['T'], int)
    for h, lst in ends.items(): hosted[h] = len(lst)
    X = []; nonX = []
    for h, lst in ends.items():
        lst.sort(key=lambda q: q['s'])
        s = np.array([q['s'] for q in lst]); pts = np.concatenate([[0.0], s, [L[h]]])
        for i, q in enumerate(lst):
            q['gap'] = min(pts[i + 1] - pts[i], pts[i + 2] - pts[i + 1])
            q['host'] = h
            partner = None
            for j, r in enumerate(lst):
                if j == i or r['side'] == q['side']: continue
                if abs(r['s'] - q['s']) < delta * (pts[i + 2] - pts[i]) and \
                        np.arccos(np.clip(-(q['d'] @ r['d']), -1, 1)) < phi:
                    partner = r; break
            if hosted[q['t']] == 0:
                (X if partner is not None else nonX).append((q, partner))
    return X, nonX, u


def move(g, items, u, shift):
    g = dict(g); A = g['A'].copy(); Bp = g['B'].copy(); done = set()
    for q, partner in items:
        if q['t'] in done or (partner is not None and partner['t'] in done): continue
        sgn = 1.0 if partner is None or q['s'] >= partner['s'] else -1.0
        dp = sgn * shift * q['gap'] * u[q['host']]
        if q['e'] == 0: A[q['t']] += dp
        else: Bp[q['t']] += dp
        done.add(q['t'])
    g['A'] = A; g['B'] = Bp
    return g, len(done)


def index(g, rho=0.2):
    L = np.linalg.norm(g['B'] - g['A'], axis=1)
    lam = np.full(g['T'], rho / L.sum())
    sc = periodic.build(g, lam, frame_lam=lam[0])
    C = periodic.periodic_C(sc)
    th, n = waves.christoffel_n(C, sc['rho'])
    return float(n.mean()), float(n.min() / n.max())


if __name__ == '__main__':
    out = []
    rules = {'balanced': 'balanced', 'regular': 'regular', 'centroid': 'centroid', 'SR': 'SR'}
    for name, rule in rules.items():
        for s in range(4):
            g = divisions.run(1000, 100 + s, rule)
            X, nonX, u = classify(g)
            rng = np.random.default_rng(s)
            n0, e0 = index(g)
            row = dict(rule=name, seed=s, n0=n0, eta0=e0, nX=len(X))
            for shift in (0.25, 0.45):
                gX, mX = move(g, X, u, shift)
                ctrl = [nonX[i] for i in rng.choice(len(nonX), size=min(mX, len(nonX)), replace=False)]
                gC, mC = move(g, ctrl, u, shift)
                row[f'nX_{shift}'], row[f'etaX_{shift}'] = index(gX)
                row[f'nC_{shift}'], row[f'etaC_{shift}'] = index(gC)
                row[f'moved_{shift}'] = (mX, mC)
            out.append(row)
            print(name, s, 'n0 %.2f | break X (%d moved): %.2f / %.2f | control: %.2f / %.2f' % (
                n0, row['moved_0.25'][0], row['nX_0.25'], row['nX_0.45'], row['nC_0.25'], row['nC_0.45']), flush=True)
            json.dump(out, open(os.path.join(B, 'results', 'breakx.json'), 'w'), indent=1)
