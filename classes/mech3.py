"""Second referee round, concern 1: rigid-joint index of the anti-aligned splits vs T at fixed
density, with the energy split. For the macroscopic strain e, the rigid-joint displacement U
gives C = U K U^T / A and its axial part C_a = U K_a U^T / A (K_a = bend 0 on the same supercell);
s = (C_a11 + C_a22) / (C11 + C22) is the share of the stiffness carried by stretching.
Also the pinned C (bend 0 solve) and a Timoshenko check (shear-deformable beams, nu = 0.3,
kappa = 5/6) at the densest point. Framed mirror-periodic tiling (mech.build).
-> results/mech3.json"""
import sys, json, time, numpy as np
import scipy.sparse as sp, scipy.sparse.linalg as spl
sys.path[:0] = ['../beyond', '../rg']
import qdiv, periodic, waves, frame
from mech import build, lines
from pisot import theta
TAU = qdiv.TAU; RHO = theta(3, 2)[0]
_ek = frame.element_k
def element_k_timo(L, c, s, EA, EI, nu=0.3, kap=5 / 6):
    """Timoshenko frame element (shear parameter phi = 12 EI / (kappa G A L^2))."""
    phi = 12 * EI / (kap * EA / (2 * (1 + nu)) * L ** 2)
    kg, R, kl = _ek(L, c, s, EA, EI)
    f = 1 / (1 + phi)
    kl = kl.copy()
    b = 12 * EI / L ** 3 * f; cc = 6 * EI / L ** 2 * f; d = (4 + phi) * EI / L * f; e = (2 - phi) * EI / L * f
    kl[:, 1, 1] = kl[:, 4, 4] = b; kl[:, 1, 4] = kl[:, 4, 1] = -b
    kl[:, 1, 2] = kl[:, 2, 1] = kl[:, 1, 5] = kl[:, 5, 1] = cc
    kl[:, 4, 2] = kl[:, 2, 4] = kl[:, 4, 5] = kl[:, 5, 4] = -cc
    kl[:, 2, 2] = kl[:, 5, 5] = d; kl[:, 2, 5] = kl[:, 5, 2] = e
    return np.einsum('nji,njk,nkl->nil', R, kl, R), R, kl

def C_split(sc, Ka):
    X, K = sc['X'], sc['K']
    P, ms = periodic._P(sc['master'], sc['shift'])
    Kr = (P.T @ K @ P).tocsc(); n = Kr.shape[0]
    keep = np.ones(n, bool); keep[:2] = False
    Kk = Kr[keep][:, keep]; Kk = Kk + 1e-13 * abs(Kk.diagonal()).max() * sp.eye(Kk.shape[0], format='csc')
    lu = spl.splu(Kk.tocsc()); U = []; x = X - 1.0
    for exx, eyy, gm in ((1, 0, 0), (0, 1, 0), (0, 0, 1)):
        ua = np.zeros(3 * len(X)); ua[0::3] = exx * x[:, 0] + 0.5 * gm * x[:, 1]; ua[1::3] = eyy * x[:, 1] + 0.5 * gm * x[:, 0]
        w = np.zeros(n); w[keep] = lu.solve(-(P.T @ (K @ ua))[keep]); U.append(ua + P @ w)
    U = np.array(U); C = U @ (K @ U.T) / 4.0; Ca = U @ (Ka @ U.T) / 4.0
    return 0.5 * (C + C.T), 0.5 * (Ca + Ca.T)

RULES = {'golden-avoid': ('largest', ('avoid', TAU ** -2, 'low'), 'perp'),
         '(5,1)-avoid': ('largest', ('avoid', RHO ** -5, 'low'), 'perp'),
         '(3,2)-avoid': ('largest', ('avoid', RHO ** -3, 'low'), 'perp'),
         'golden': ('largest', TAU ** -2, 'perp'),
         'uniform': ('largest', 'uniform', 'perp')}
if __name__ == '__main__':
    out = []; t0 = time.time()
    for T in (300, 1000, 3000, 10000):
        for nm, rule in RULES.items():
            for seed in ((0, 1, 2) if nm == 'uniform' and T < 10000 else (0,)):
                g = qdiv.run(T, rule, seed=seed); L = np.linalg.norm(g['B'] - g['A'], axis=1)
                m = int((lines(g) > 1 - 1e-6).sum())
                for rho in (0.05, 0.1, 0.2, 0.4):
                    lam = np.full(T, rho / L.sum())
                    for model in (('EB', 'timo') if rho >= 0.2 else ('EB',)):
                        frame.element_k = element_k_timo if model == 'timo' else _ek
                        sc = build(g, lam, lam[0], bend=1.0)
                        Ka = build(g, lam, lam[0], bend=0.0)['K']
                        C, Ca = C_split(sc, Ka)
                        Cp = periodic.periodic_C(dict(sc, K=Ka))
                        th, n = waves.christoffel_n(C, sc['rho']); thp, npn = waves.christoffel_n(Cp, sc['rho'])
                        o = dict(rule=nm, T=T, seed=seed, rho_t=rho, model=model, m=m, n_mean=float(n.mean()),
                                 n_x=float(n[0]), npin_x=float(npn[0]), npin_mean=float(npn.mean()),
                                 s=float((Ca[0, 0] + Ca[1, 1]) / (C[0, 0] + C[1, 1])),
                                 pin_ratio=float((Cp[0, 0] + Cp[1, 1]) / (C[0, 0] + C[1, 1])), C=C.tolist())
                        out.append(o)
                        print(nm, T, seed, rho, model, 'm=%d n=%.2f npin_x=%.2f s=%.3f pin/rig=%.3f' %
                              (m, o['n_mean'], o['npin_x'], o['s'], o['pin_ratio']), round(time.time() - t0), flush=True)
                frame.element_k = _ek
        json.dump(out, open('results/mech3.json', 'w'), indent=1)
