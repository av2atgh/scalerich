"""Ligament (dual) graph of the divisions: ligaments = maximal straight chains of collinear
segments; two ligaments are linked if they meet at a junction. Frame ligaments excluded (they
span every sample). Reports <k>, <k^2>/<k> and max degree vs T.
-> results/ligaments.json"""
import sys, json, numpy as np
import supercond as sc

def ligament_graph(g):
    X, S = sc.wires(g)
    d = X[S[:, 1]] - X[S[:, 0]]; H = np.abs(d[:, 1]) < 1e-12
    par = np.arange(len(S))
    def find(a):
        while par[a] != a: par[a] = par[par[a]]; a = par[a]
        return a
    inc = [[] for _ in range(len(X))]
    for e, (a, b) in enumerate(S): inc[a].append(e); inc[b].append(e)
    u = d / np.linalg.norm(d, axis=1)[:, None]
    for v in range(len(X)):                            # join segments that continue straight through v
        es = inc[v]
        for i in range(len(es)):
            for j in range(i + 1, len(es)):
                e, f = es[i], es[j]
                se = u[e] if S[e, 0] == v else -u[e]; sf = u[f] if S[f, 0] == v else -u[f]
                if se @ sf < -1 + 1e-9: par[find(e)] = find(f)
    lig = np.array([find(e) for e in range(len(S))])
    L = np.linalg.norm(d, axis=1)
    onf = np.zeros(len(S), bool)                       # frame segments
    for ax in (0, 1):
        for c in (0.0, 1.0):
            onf |= (np.abs(X[S[:, 0], ax] - c) < 1e-9) & (np.abs(X[S[:, 1], ax] - c) < 1e-9)
    frame = set(lig[onf])
    edges = set()
    for v in range(len(X)):
        ls = sorted({lig[e] for e in inc[v]})
        for i in range(len(ls)):
            for j in range(i + 1, len(ls)):
                if ls[i] not in frame and ls[j] not in frame: edges.add((ls[i], ls[j]))
    ids = [l for l in np.unique(lig) if l not in frame]
    deg = {l: 0 for l in ids}
    for a, b in edges: deg[a] += 1; deg[b] += 1
    k = np.array(list(deg.values()), float)
    ln = {l: 0.0 for l in ids}
    for e in range(len(S)):
        if lig[e] in ln: ln[lig[e]] += L[e]
    return k, np.array(list(ln.values()))

if __name__ == '__main__':
    Ts = [250, 500, 1000, 2000, 4000, 8000, 16000]
    out = {}
    for nm in ['periodic', 'golden-aligned', 'golden-anti', 'r0.3-aligned', 'uniform', 'SR']:
        for T in Ts:
            for sd in ((0, 1) if nm in ('uniform', 'SR') else (0,)):
                if nm == 'uniform': g = sc.qdiv.run(T, ('largest', 'uniform', 'perp'), seed=sd)
                elif nm == 'SR': g = sc.divisions.run(T, sd, 'SR')
                else: g = sc.CASES[nm](T)
                k, ln = ligament_graph(g)
                r = dict(rule=nm, T=T, seed=sd, n=len(k), k1=k.mean(), k2k1=(k ** 2).mean() / k.mean(),
                         kmax=int(k.max()))
                out['%s_T%d_s%d' % (nm, T, sd)] = r
                print('%-15s T=%-6d s%d n=%-6d <k>=%.2f <k2>/<k>=%.2f kmax=%d ' %
                      (nm, T, sd, r['n'], r['k1'], r['k2k1'], r['kmax']), flush=True)
    json.dump(out, open('results/ligaments.json', 'w'))
    for nm in ['periodic', 'golden-aligned', 'golden-anti', 'r0.3-aligned', 'uniform', 'SR']:
        rs = [v for v in out.values() if v['rule'] == nm]
        T = np.array([v['T'] for v in rs]); y = np.array([v['k2k1'] for v in rs]); m = np.array([v['kmax'] for v in rs])
        sel = T >= 1000
        print('%-15s exponent of <k2>/<k> %.2f, of kmax %.2f (T>=1000)' % (nm, np.polyfit(np.log(T[sel]), np.log(y[sel]), 1)[0], np.polyfit(np.log(T[sel]), np.log(m[sel]), 1)[0]))
