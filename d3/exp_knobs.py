"""Goal 5: knobs in 3D -- thickness/jamming and the three-axis orientation set.

Fill-fraction theory (zero-thickness area profile c(t) = <a_t> t^(2/3) from
the iso runs): a plate removes lam_t a_t of free volume and the free volume
F spread over t cells gives a_t = c(t) (F/t)^(2/3), so
   dF/dt = -lam0 c(t) F^(2/3) t^(-alpha-2/3),
   F_inf^(1/3) = 1 - (lam0/3) sum_t c(t) t^(-alpha-2/3).
The sum diverges for alpha <= 1/3 (jammed); for alpha > 1/3 it jams only if
lam0 > 3 / [c zeta(alpha+2/3)], and otherwise rho = 1 - F_inf ~ lam0 c zeta(alpha+2/3).

usage: python3 exp_knobs.py  -> results/knobs.txt
"""
import os
import numpy as np
from common import load, plate_arrays, gamma_mle, cohort_beta, RES, DATA
from exp_cluster import adjacency

lines = []


def out(s):
    print(s, flush=True); lines.append(s)


# c(t) from the isotropic zero-thickness runs
cs = []
for nm in ('iso_s1', 'iso_s2', 'iso_s3'):
    S, _, _ = load(nm)
    cs.append(S[:, 5] * S[:, 0] ** (2 / 3))
cbar = np.mean(cs, axis=0)
T0 = len(cbar)
edges = np.unique(np.logspace(0, np.log10(T0), 60).astype(int))
csm = np.zeros(T0)
for a, b in zip(edges[:-1], edges[1:]):
    csm[a - 1:b - 1] = cbar[a - 1:b - 1].mean()
csm[edges[-1] - 1:] = cbar[edges[-1] - 1:].mean()
c_inf = cbar[T0 // 2:].mean()
out(f'c(t)=<a_t>t^(2/3): t=1 {cbar[0]:.3f} (one draw per seed), t~1e2 {cbar[90:110].mean():.3f}, '
    f't~1e4 {cbar[9000:11000].mean():.3f}, asymptote {c_inf:.4f}')


def F_pred(lam0, alpha, T):
    F = 1.0; Fs = np.zeros(T)
    for t in range(1, T + 1):
        c = csm[t - 1] if t <= T0 else c_inf
        F = max(F - lam0 * c * F ** (2 / 3) * t ** (-alpha - 2 / 3), 0.0)
        Fs[t - 1] = F
    return Fs


def zeta_tail(s, n=10**6):
    t = np.arange(1, n + 1, dtype=float)
    return np.sum(t ** -s) + n ** (1 - s) / (s - 1)


for nm, (lam0, alpha) in (('lam_a0.8_l0.05', (0.05, 0.8)), ('lam_a0.5_l0.02', (0.02, 0.5)),
                          ('lam_a0.4_l0.02', (0.02, 0.4)), ('lam_a0.25_l0.02', (0.02, 0.25)),
                          ('lam_a0.0_l0.01', (0.01, 0.0))):
    S, E, C = load(nm)
    g = S[:, 0] > 0
    T = int(S[g, 0].max())
    Fp = F_pred(lam0, alpha, len(S))
    chk = [10**2, 10**3, 10**4, 10**5, len(S)]
    fs = ' '.join(f't={c}: {S[c - 1, 9]:.4f}/{Fp[c - 1]:.4f}' for c in chk)
    jam = np.where(S[:, 0] == -2)[0]
    jt = jam[0] + 1 if len(jam) else None
    lin = lam0 * c_inf * zeta_tail(alpha + 2 / 3) if alpha > 1 / 3 else float('inf')
    out(f'== {nm} (lam0={lam0}, alpha={alpha}): plates {g.sum()}, jammed at t={jt}; free volume sim/pred: {fs}; '
        f'linear rho estimate lam0 c zeta(alpha+2/3) = {lin:.4f}')
    P = plate_arrays(S, E)
    k = P['k']; ok = P['ok']
    mm = S[g, 3].mean()
    st = f'  <m>={mm:.3f}; MLE all plates k>=20,40,80: ' + ' '.join(f'{gamma_mle(k[ok], q)[0]:.3f}' for q in (20, 40, 80))
    if g.sum() > 20000:
        st += f'; cohort beta [1e2,T/1e2] {cohort_beta(k.astype(float), ok, 100, g.sum() / 100):.3f}'
    pe = E['h'] >= 0
    A = adjacency(E['t'][pe], E['h'][pe], len(S) + 1)
    kk = np.asarray(A.sum(1)).ravel()
    tri = np.asarray((A @ A).multiply(A).sum(1)).ravel() / 2
    m2 = kk >= 2
    Cl = 2 * tri[m2] / (kk[m2] * (kk[m2] - 1))
    hub = kk >= 50
    st += f'; <C>={Cl.mean():.4f}; <k C(k)> for k>=50: {np.mean(2 * tri[hub] / (kk[hub] - 1)) if hub.sum() else float("nan"):.3f}'
    out(st)

# three orthogonal orientations -> boxes
for nm in ('axes_s1', 'axes_s2', 'axes_s3'):
    S, E, C = load(nm)
    P = plate_arrays(S, E); T = P['T']; k = P['k']; ok = P['ok']
    m = S[:, 3]
    out(f'== {nm} (normals in {{x,y,z}}): T={T}; m values {np.unique(m).astype(int).tolist()}; cell faces {np.unique(C["nf"]).tolist()}; '
        f'<k>={k[1:].mean():.3f}; incidence PA b = 4/6 = area rule 2/3')
    out('  MLE all plates k>=20,40,80,160,320: ' + ' '.join(f'{gamma_mle(k[ok], q)[0]:.3f}' for q in (20, 40, 80, 160, 320))
        + f'; cohort beta [1e2,T/1e2] = {cohort_beta(k.astype(float), ok, 100, T / 100):.3f}'
        + ' ; K^1.5 P(k>=K) K=20..320: ' + ' '.join(f'{K**1.5 * np.mean(k[ok] >= K):.2f}' for K in (20, 40, 80, 160, 320)))
    if T <= 300000:
        pe = E['h'] >= 0
        A = adjacency(E['t'][pe], E['h'][pe], T + 1)
        kk = np.asarray(A.sum(1)).ravel()
        tri = np.asarray((A @ A).multiply(A).sum(1)).ravel() / 2
        m2 = kk >= 2
        hub = kk >= 50
        out(f'  <C>={np.mean(2 * tri[m2] / (kk[m2] * (kk[m2] - 1))):.4f}; <k C(k)> for k>=50: {np.mean(2 * tri[hub] / (kk[hub] - 1)):.3f}; '
            f'triangles per received link (fit, k>=20): {np.polyfit(P["recv"][kk >= 20], tri[kk >= 20], 1)[0]:.3f}')
open(os.path.join(RES, 'knobs.txt'), 'w').write('\n'.join(lines) + '\n')
