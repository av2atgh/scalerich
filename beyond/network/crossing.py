"""Quasi-crossings: two stems ending on opposite sides of the same host at nearly the same point and
nearly collinear, so that axial load passes through the host without bending it.
f_X = fraction of stem ends that have such a partner (offset < delta * local span, angle < phi_max)."""
import numpy as np


def quasi_crossings(g, delta=0.1, phi_max=np.deg2rad(15)):
    T = g['T']; A = g['A']; Bp = g['B']; v = Bp - A
    L = np.linalg.norm(v, axis=1); u = v / L[:, None]
    ends = {}
    for e, (P, hosts) in enumerate(((A, g['hostA']), (Bp, g['hostB']))):
        for t in range(T):
            h = hosts[t]
            if h < 0: continue
            # direction pointing from the junction into the stem
            d = u[t] if e == 0 else -u[t]
            s = (P[t] - A[h]) @ u[h]
            side = np.sign(d @ np.array([-u[h, 1], u[h, 0]]))
            ends.setdefault(h, []).append((s, side, d, t))
    nX = 0; ntot = 0
    for h, lst in ends.items():
        s = np.array([q[0] for q in lst]); order = np.argsort(s); s = s[order]
        lst = [lst[i] for i in order]
        # local span: distance to neighbouring junctions on h
        pts = np.concatenate([[0.0], s, [L[h]]])
        for i, (si, sd, di, ti) in enumerate(lst):
            loc = min(pts[i + 1] - pts[i], pts[i + 2] - pts[i + 1]) if len(lst) > 1 else L[h]
            ntot += 1
            for j in range(len(lst)):
                if j == i or lst[j][1] == sd: continue
                sj, _, dj, _ = lst[j]
                ang = np.arccos(np.clip(-(di @ dj), -1, 1))   # collinear through-path: dj = -di
                span_here = max(abs(sj - si), 0)
                ref = max(pts[i + 2] - pts[i], 1e-12)
                if span_here < delta * ref and ang < phi_max:
                    nX += 1; break
    return nX / max(ntot, 1)
