"""Long-wavelength elastic-wave properties of a recursive-division frame.

Frame: rg/frame.py (Euler-Bernoulli, rigid joints, E0 = rho0 = 1, unit depth,
EA = lam, EI = lam^3/12), no outer frame. Effective in-plane stiffness C
(Voigt 11,22,12 with engineering shear) from three kinematic (affine)
boundary tests: every boundary node gets u = eps.(x - x_c), rotations free;
C_ij = u_i^T K u_j (one factorisation). Kinematic BCs bound the periodic C
from above, so n is a lower estimate; the comparison between rules is what matters.

Waves (homogenised, ρ_eff = ρ̄): Christoffel Γ(θ); the quasi-longitudinal
velocity is v(θ)^2 = λ_max Γ(θ) / ρ̄, and n(θ) = c0 / v(θ) with c0 = 1. The paper's
n_eff is n(0). Isotropy of n: η_n = min_θ n / max_θ n.

Bandwidth: lowest local resonance ω1 of the frame with all boundary nodes
clamped (lumped mass λL per element, rotary inertia λL^3/24 per end), in
units c0 / L_cell. The homogenised n is valid only for ω well below ω1.
"""
import os, sys
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spl

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', 'rg'))
import srtree, frame  # noqa: E402


def system(g, lam):
    nodes, seg, sl, bn = srtree.network(g)
    lam_s = lam[sl]
    K, dof, R, kl, L = frame.assemble(nodes, seg, lam_s)
    bnd = np.unique(np.concatenate([bn[s] for s in (-1, -2, -3, -4)]))
    rho = float(np.sum(lam * np.linalg.norm(g['B'] - g['A'], axis=1)))
    return nodes, seg, lam_s, L, K, bnd, rho


def stiffness_tensor(nodes, K, bnd):
    N = len(nodes)
    x = nodes - 0.5
    fixed = np.zeros(3 * N, bool)
    fixed[3 * bnd] = fixed[3 * bnd + 1] = True
    free = ~fixed
    Kff = K[free][:, free].tocsc()
    Kff = Kff + 1e-13 * abs(Kff.diagonal()).max() * sp.eye(Kff.shape[0], format='csc')
    lu = spl.splu(Kff)
    Kfx = K[free][:, fixed]
    U = []
    for e in ((1, 0, 0), (0, 1, 0), (0, 0, 1)):        # eps_xx, eps_yy, gamma_xy
        exx, eyy, g = e
        u = np.zeros(3 * N)
        u[0::3] = exx * x[:, 0] + 0.5 * g * x[:, 1]
        u[1::3] = eyy * x[:, 1] + 0.5 * g * x[:, 0]
        uf = u.copy()
        uf[free] = lu.solve(-(Kfx @ u[fixed]))
        U.append(uf)
    U = np.array(U)
    C = U @ (K @ U.T)
    return 0.5 * (C + C.T), Kff, free


def christoffel_n(C, rho, nth=72):
    th = np.linspace(0, np.pi, nth, endpoint=False)
    c, s = np.cos(th), np.sin(th)
    C11, C22, C12, C16, C26, C66 = C[0, 0], C[1, 1], C[0, 1], C[0, 2], C[1, 2], C[2, 2]
    G11 = C11 * c * c + 2 * C16 * c * s + C66 * s * s
    G22 = C66 * c * c + 2 * C26 * c * s + C22 * s * s
    G12 = C16 * c * c + (C12 + C66) * c * s + C26 * s * s
    lmax = 0.5 * (G11 + G22) + np.sqrt(0.25 * (G11 - G22) ** 2 + G12 ** 2)
    return th, np.sqrt(rho / lmax)


def lowest_resonance(nodes, seg, lam_s, L, Kff, free, k=1):
    N = len(nodes)
    m = lam_s * L
    M = np.zeros(3 * N)
    for e in (0, 1):
        np.add.at(M, 3 * seg[:, e], 0.5 * m)
        np.add.at(M, 3 * seg[:, e] + 1, 0.5 * m)
        np.add.at(M, 3 * seg[:, e] + 2, m * L ** 2 / 24)
    Mf = sp.diags(M[free]).tocsc()
    w2 = spl.eigsh(Kff, k=k, M=Mf, sigma=0, which='LM', return_eigenvectors=False)
    return float(np.sqrt(np.min(w2)))


def analyse(g, lam, resonance=True):
    nodes, seg, lam_s, L, K, bnd, rho = system(g, lam)
    return analyse_frame(nodes, seg, lam_s, L, K, bnd, rho, resonance)


def frame_system(nodes, seg, lam_s, bnd):
    K, dof, R, kl, L = frame.assemble(nodes, seg, lam_s)
    return nodes, seg, lam_s, L, K, np.asarray(bnd), float(np.sum(lam_s * L))


def analyse_frame(nodes, seg, lam_s, L, K, bnd, rho, resonance=True):
    C, Kff, free = stiffness_tensor(nodes, K, bnd)
    th, n = christoffel_n(C, rho)
    out = dict(rho=rho, n_x=float(n[0]), n_mean=float(n.mean()), n_min=float(n.min()),
               n_max=float(n.max()), eta_n=float(n.min() / n.max()), C=C.tolist())
    if resonance:
        out['omega1'] = lowest_resonance(nodes, seg, lam_s, L, Kff, free)
    return out
