"""Goal 6: 3D beams (rods).  Zero-width lines in 3D never meet, so the model
needs a capture radius.  Natural version: rod t has radius R_t = lam0 t^-alpha / 2
(half the thickness); nucleation point uniform in the cube, rejected if it lies
within R_j + R_t of an existing rod axis; the rod grows both ways along an
isotropic direction until its tip comes within R_j + R_t of the axis of rod j
(surfaces touch) or reaches the cube.  Dual: rod = node, touch = edge (<= 2
links per newcomer, as in 2D).

Mean field: a rod is hit at a rate proportional to its capture cross-section
l_i (R_i + R_t) ~ l_i R_i (old rods are thicker).  The mean free path is
l_t ~ 1/W(t) with W = sum_j l_j R_j, so W ~ t^((1-alpha)/2), l_t ~ t^-((1-alpha)/2),
and k_i ~ (t/i)^((1+alpha)/2): beta = (1+alpha)/2, gamma = (3+alpha)/(1+alpha).
Rod volume sum_t l_t R_t^2 ~ sum t^-((1+3 alpha)/2) converges iff alpha > 1/3.

usage: python3 beams3d.py -> results/beams3d.txt  (brute force O(T^2), numba)
"""
import os, time
import numpy as np
from numba import njit
from common import gamma_mle, cohort_beta, RES


@njit(cache=True)
def seg_dist2(p, a, b):
    d0 = b[0] - a[0]; d1 = b[1] - a[1]; d2 = b[2] - a[2]
    L2 = d0 * d0 + d1 * d1 + d2 * d2
    s = ((p[0] - a[0]) * d0 + (p[1] - a[1]) * d1 + (p[2] - a[2]) * d2) / L2 if L2 > 0 else 0.0
    s = min(max(s, 0.0), 1.0)
    x = a[0] + s * d0 - p[0]; y = a[1] + s * d1 - p[1]; z = a[2] + s * d2 - p[2]
    return x * x + y * y + z * z


@njit(cache=True)
def ray_sphere(p, u, c, r2):
    """first s >= 0 with |p + s u - c|^2 <= r2 (u unit), inf if none."""
    w0 = p[0] - c[0]; w1 = p[1] - c[1]; w2 = p[2] - c[2]
    b = u[0] * w0 + u[1] * w1 + u[2] * w2
    cc = w0 * w0 + w1 * w1 + w2 * w2 - r2
    disc = b * b - cc
    if disc < 0:
        return np.inf
    s = -b - np.sqrt(disc)
    if s < 0:
        s2 = -b + np.sqrt(disc)
        return 0.0 if s2 >= 0 else np.inf
    return s


@njit(cache=True)
def ray_capsule(p, u, a, b, r):
    """first s >= 0 at which the point p + s u is within r of segment [a, b]."""
    r2 = r * r
    best = min(ray_sphere(p, u, a, r2), ray_sphere(p, u, b, r2))
    d = np.empty(3)
    d[0] = b[0] - a[0]; d[1] = b[1] - a[1]; d[2] = b[2] - a[2]
    L = np.sqrt(d[0] ** 2 + d[1] ** 2 + d[2] ** 2)
    if L > 0:
        d /= L
        w = np.empty(3); w[0] = p[0] - a[0]; w[1] = p[1] - a[1]; w[2] = p[2] - a[2]
        ud = u[0] * d[0] + u[1] * d[1] + u[2] * d[2]
        wd = w[0] * d[0] + w[1] * d[1] + w[2] * d[2]
        # perpendicular components
        up0 = u[0] - ud * d[0]; up1 = u[1] - ud * d[1]; up2 = u[2] - ud * d[2]
        wp0 = w[0] - wd * d[0]; wp1 = w[1] - wd * d[1]; wp2 = w[2] - wd * d[2]
        A = up0 * up0 + up1 * up1 + up2 * up2
        B = 2 * (up0 * wp0 + up1 * wp1 + up2 * wp2)
        C = wp0 * wp0 + wp1 * wp1 + wp2 * wp2 - r2
        if A > 1e-18:
            disc = B * B - 4 * A * C
            if disc >= 0:
                sq = np.sqrt(disc)
                for s in ((-B - sq) / (2 * A), (-B + sq) / (2 * A)):
                    if s >= 0:
                        h = wd + s * ud
                        if 0 <= h <= L and s < best:
                            best = s
                        break
    return best


@njit(cache=True)
def ray_box(p, u):
    best = np.inf
    for k in range(3):
        if u[k] > 0:
            s = (1 - p[k]) / u[k]
        elif u[k] < 0:
            s = -p[k] / u[k]
        else:
            continue
        if s < best:
            best = s
    return best


@njit(cache=True)
def grow(T, lam0, alpha, seed, maxtries):
    np.random.seed(seed)
    A = np.zeros((T + 1, 3)); B = np.zeros((T + 1, 3)); R = np.zeros(T + 1)
    lo = np.zeros((T + 1, 3)); hi = np.zeros((T + 1, 3))
    ea = np.full((T + 1, 2), -1, np.int64)
    vol = 0.0
    p = np.empty(3); u = np.empty(3); mu = np.empty(3)
    n = 0
    tries_tot = 0
    for t in range(1, T + 1):
        Rt = 0.5 * lam0 * t ** (-alpha)
        ok = False
        for tr in range(maxtries):
            tries_tot += 1
            p[0] = np.random.random(); p[1] = np.random.random(); p[2] = np.random.random()
            ok = True
            for j in range(1, n + 1):
                rr = R[j] + Rt
                if p[0] < lo[j, 0] - rr or p[0] > hi[j, 0] + rr or p[1] < lo[j, 1] - rr or p[1] > hi[j, 1] + rr \
                        or p[2] < lo[j, 2] - rr or p[2] > hi[j, 2] + rr:
                    continue
                if seg_dist2(p, A[j], B[j]) <= rr * rr:
                    ok = False
                    break
            if ok:
                break
        if not ok:
            return A[:n + 1], B[:n + 1], R[:n + 1], ea[:n + 1], n, tries_tot
        g = np.random.normal(0, 1, 3)
        g /= np.sqrt(g[0] ** 2 + g[1] ** 2 + g[2] ** 2)
        ends = np.empty((2, 3))
        for side in range(2):
            sg = 1.0 if side == 0 else -1.0
            u[0] = sg * g[0]; u[1] = sg * g[1]; u[2] = sg * g[2]
            best = ray_box(p, u); who = -1
            for j in range(1, n + 1):
                rr = R[j] + Rt
                # slab test of the ray segment [0, best] against the padded box
                s0 = 0.0; s1 = best; skip = False
                for k in range(3):
                    l = lo[j, k] - rr; h = hi[j, k] + rr
                    if abs(u[k]) < 1e-15:
                        if p[k] < l or p[k] > h:
                            skip = True
                            break
                    else:
                        q0 = (l - p[k]) / u[k]; q1 = (h - p[k]) / u[k]
                        if q0 > q1:
                            q0, q1 = q1, q0
                        s0 = max(s0, q0); s1 = min(s1, q1)
                        if s0 > s1:
                            skip = True
                            break
                if skip:
                    continue
                s = ray_capsule(p, u, A[j], B[j], rr)
                if s < best:
                    best = s; who = j
            for k in range(3):
                ends[side, k] = p[k] + best * u[k]
            ea[t, side] = who
        n = t
        for k in range(3):
            A[t, k] = ends[0, k]; B[t, k] = ends[1, k]
            lo[t, k] = min(ends[0, k], ends[1, k]); hi[t, k] = max(ends[0, k], ends[1, k])
        R[t] = Rt
    return A, B, R, ea, n, tries_tot


if __name__ == '__main__':
    lines = []

    def out(s):
        print(s, flush=True); lines.append(s)

    for T, lam0, alpha in ((100000, 0.02, 0.5), (100000, 0.05, 0.8), (100000, 0.004, 0.0)):
        t0 = time.time()
        A, B, R, ea, n, tries = grow(T, lam0, alpha, 1, 2000)
        L = np.linalg.norm(A - B, axis=1)
        ii = np.arange(n + 1)
        vfrac = np.sum(np.pi * R[1:] ** 2 * L[1:])
        src = np.repeat(ii, 2); dst = ea[:n + 1].ravel()
        ok = dst > 0
        k = np.bincount(src[ok], minlength=n + 1) + np.bincount(dst[ok], minlength=n + 1)
        recv = np.bincount(dst[ok], minlength=n + 1)
        sel = np.ones(n + 1, bool); sel[0] = False
        lslope = -np.polyfit(np.log(ii[100:]), np.log(L[100:] + 1e-12), 1)[0]
        bc = cohort_beta(k.astype(float), sel, 30, n / 30)
        # fixed-rod growth exponent: hits received by rods 100..1000 up to time t
        hit_t = np.repeat(ii, 2)[ok]; hit_h = dst[ok]
        ts = np.unique(np.logspace(np.log10(3000), np.log10(n), 8).astype(int))
        mk = (hit_h >= 100) & (hit_h < 1000)
        ls = [np.mean(np.log(1 + np.bincount(hit_h[mk][hit_t[mk] <= tc], minlength=1000)[100:1000])) for tc in ts]
        fs = np.polyfit(np.log(ts), ls, 1)[0]
        # hits vs capture cross-section l_i R_i
        s2 = (ii > 100) & (ii < n // 20)
        cbn = np.floor(np.log(np.maximum(ii, 1)) / np.log(1.25)).astype(int)
        y = np.log(recv[s2] + 0.5); xl = np.log(L[s2]); xr = np.log(R[s2])
        for arr in (y, xl):
            mu = np.bincount(cbn[s2], weights=arr) / np.maximum(np.bincount(cbn[s2]), 1)
            arr -= mu[cbn[s2]]
        el = np.polyfit(xl, y, 1)[0]
        out(f'== rods lam0={lam0} alpha={alpha}: rods placed {n} of {T} (nucleation tries {tries}), volume fraction {vfrac:.4f}, '
            f'[{time.time() - t0:.0f}s]')
        out(f'  length exponent l_t ~ t^-{lslope:.3f} (theory (1-alpha)/2 = {(1 - alpha) / 2:.3f}); <k>={k[1:].mean():.3f}; '
            f'fraction of ends on rods {ok.mean():.3f}')
        out(f'  cohort beta [30,n/30] = {bc:.3f}; fixed rods 100-1000 d<ln(1+recv)>/d ln t = {fs:.3f}  (theory (1+alpha)/2 = {(1 + alpha) / 2:.3f}); '
            f'within-cohort elasticity of hits on length = {el:.3f} (1)')
        out(f'  MLE k>=10,20,40,80: ' + ' '.join(f'{gamma_mle(k[1:], q)[0]:.3f}' for q in (10, 20, 40, 80))
            + f'  (theory (3+alpha)/(1+alpha) = {(3 + alpha) / (1 + alpha):.3f}); kmax={k.max()}')
    open(os.path.join(RES, 'beams3d.txt'), 'w').write('\n'.join(lines) + '\n')
