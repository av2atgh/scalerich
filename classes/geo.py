"""Second referee round, concern 1(b) and 2: slenderness and near-crossings.
Lengths in units of 1/L_tot (L_tot = total ligament length), so that a length ell is ell/lambda =
ell L_tot / rho_bar thicknesses at areal density rho_bar (lambda = rho_bar / L_tot).
 seg  : beam segments between consecutive joints (zero-length X-crossing pieces dropped)
 delta: offset between a stem and the nearest stem on the other side of the same host (> 0)
 snap : finite-thickness test -- parallel ligaments whose centre lines are closer than lambda
        overlap; collinear coordinates within lambda are merged (single linkage) and the
        spanning lines m_snap(rho_bar) recounted (frame excluded).
 longest: longest straight line other than the first cut (minor point 2).
-> results/geo.json"""
import json, numpy as np, qdiv
from mech3 import RULES
from mech import lines

def segs(g):
    T = g['T']; A, B = g['A'], g['B']; out = []
    on = {}
    for t in range(T):
        for P, h in ((A[t], g['hostA'][t]), (B[t], g['hostB'][t])):
            on.setdefault(int(h), []).append(P)
    for t in range(T):
        d = B[t] - A[t]; L = np.linalg.norm(d); s = [0.0, L] + [float((P - A[t]) @ d / L) for P in on.get(t, [])]
        s = np.unique(np.round(np.sort(s), 12)); out += list(np.diff(s))
    return np.array(out), on

def offsets(g):
    """for every host t: stems from below/left vs above/right; nearest opposite distance."""
    T = g['T']; A, B = g['A'], g['B']; side = {}
    for t in range(T):
        for P, Q, h in ((A[t], B[t], g['hostA'][t]), (B[t], A[t], g['hostB'][t])):
            h = int(h)
            if h < 0: continue
            vert_host = abs(A[h, 0] - B[h, 0]) < 1e-12
            pos = P[1] if vert_host else P[0]; sd = (Q[0] > P[0]) if vert_host else (Q[1] > P[1])
            side.setdefault(h, ([], []))[int(sd)].append(pos)
    d = []
    for h, (a, b) in side.items():
        if not a or not b: continue
        a = np.sort(a); b = np.sort(b)
        for x, y in ((a, b), (b, a)):
            i = np.searchsorted(y, x); c = np.minimum(np.abs(x - y[np.clip(i, 0, len(y) - 1)]), np.abs(x - y[np.clip(i - 1, 0, len(y) - 1)]))
            d += list(c[c > 1e-9])
    return np.array(d)

def m_snap(g, tol):
    A, B = g['A'], g['B']; n = 0
    for ax in (0, 1):
        k = np.abs(A[:, ax] - B[:, ax]) < 1e-12
        c = A[k, ax].copy(); lo = np.minimum(A[k, 1 - ax], B[k, 1 - ax]); hi = np.maximum(A[k, 1 - ax], B[k, 1 - ax])
        u = np.unique(c); lab = np.r_[0, np.cumsum(np.diff(u) > tol)]; cl = lab[np.searchsorted(u, c)]
        for q in np.unique(cl):
            s = cl == q; iv = sorted(zip(lo[s], hi[s])); a, e = iv[0]
            for x, y in iv[1:]:
                if x <= e + tol: e = max(e, y)
                else: n += e - a > 1 - 2 * tol; a, e = x, y
            n += e - a > 1 - 2 * tol
    return n

if __name__ == '__main__':
    out = []
    for nm in ('golden-avoid', '(5,1)-avoid', '(3,2)-avoid', 'golden'):
        for T in (1000, 10000, 100000):
            g = qdiv.run(T, RULES[nm]); Lt = np.linalg.norm(g['B'] - g['A'], axis=1).sum()
            sg, _ = segs(g); sg = sg[sg > 1e-9] * Lt; dl = offsets(g) * Lt
            ms = {str(r): int(m_snap(g, r / Lt)) for r in (0.025, 0.05, 0.1, 0.2, 0.4)}
            ln = np.sort(lines(g))[::-1]
            o = dict(rule=nm, T=T, Ltot=float(Lt), seg_q=np.quantile(sg, [0, .01, .05, .5]).tolist(),
                     delta_q=(np.quantile(dl, [0, .01, .05, .5]).tolist() if len(dl) else None), n_delta=int(len(dl)),
                     frac_delta_lt={str(r): float((dl < r).mean()) for r in (0.1, 0.2, 0.4)} if len(dl) else None,
                     frac_seg_lt3={str(r): float((sg < 3 * r).mean()) for r in (0.05, 0.1, 0.2, 0.4)},
                     m=int((ln > 1 - 1e-6).sum()), m_snap=ms, longest2=float(ln[1]))
            out.append(o); print(json.dumps(o), flush=True)
    json.dump(out, open('results/geo.json', 'w'), indent=1)
