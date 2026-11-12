"""Linear-elastic beam frame on SR networks (Goal 3, elastic).

Euler-Bernoulli frame elements, rigid joints, E_s = 1, unit depth:
  axial EA = lam, bending EI = lam^3/12 (thickness lam of the owning ligament).
Uniaxial test along x: nodes on the left side ux = 0, on the right side
ux = -eps (eps = 1, linear); one left node also uy = 0; rotations free;
top/bottom sides free.  E_x = 2U/eps^2 (unit square).  y analogous.
Energy split per segment into axial and bending parts.

stiffness(nodes, seg, lam_seg, bn, axis) -> E, U_axial_seg, U_bend_seg
"""
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spl


def element_k(L, c, s, EA, EI):
    """global 6x6 stiffness for each element (vectorised): returns (n,6,6)"""
    n = len(L)
    k = np.zeros((n, 6, 6))
    a = EA / L; b = 12 * EI / L ** 3; cc = 6 * EI / L ** 2; d = 4 * EI / L; e = 2 * EI / L
    kl = np.zeros((n, 6, 6))
    kl[:, 0, 0] = kl[:, 3, 3] = a; kl[:, 0, 3] = kl[:, 3, 0] = -a
    kl[:, 1, 1] = kl[:, 4, 4] = b; kl[:, 1, 4] = kl[:, 4, 1] = -b
    kl[:, 1, 2] = kl[:, 2, 1] = kl[:, 1, 5] = kl[:, 5, 1] = cc
    kl[:, 4, 2] = kl[:, 2, 4] = kl[:, 4, 5] = kl[:, 5, 4] = -cc
    kl[:, 2, 2] = kl[:, 5, 5] = d; kl[:, 2, 5] = kl[:, 5, 2] = e
    R = np.zeros((n, 6, 6))
    for o in (0, 3):
        R[:, o, o] = c; R[:, o, o + 1] = s; R[:, o + 1, o] = -s; R[:, o + 1, o + 1] = c
        R[:, o + 2, o + 2] = 1
    return np.einsum('nji,njk,nkl->nil', R, kl, R), R, kl


def assemble(nodes, seg, lam_seg, bend=1.0, cap=True):
    """cap: stiffness length L_eff = max(L, lam) -- a segment shorter than its
    thickness lies inside a junction block; Euler-Bernoulli 12EI/L^3 and 4EI/L
    would otherwise diverge as L -> 0 and break the Voigt bound."""
    N = len(nodes)
    d = nodes[seg[:, 1]] - nodes[seg[:, 0]]
    L = np.linalg.norm(d, axis=1); c = d[:, 0] / L; s = d[:, 1] / L
    EA = lam_seg; EI = bend * lam_seg ** 3 / 12
    Le = np.maximum(L, lam_seg) if cap else L
    kg, R, kl = element_k(Le, c, s, EA, EI)
    dof = np.c_[3 * seg[:, 0], 3 * seg[:, 0] + 1, 3 * seg[:, 0] + 2,
                3 * seg[:, 1], 3 * seg[:, 1] + 1, 3 * seg[:, 1] + 2]
    I = np.repeat(dof, 6, axis=1).ravel(); J = np.tile(dof, (1, 6)).ravel()
    K = sp.coo_matrix((kg.ravel(), (I, J)), shape=(3 * N, 3 * N)).tocsr()
    return K, dof, R, kl, L


def stiffness(nodes, seg, lam_seg, bn, axis='x', bend=1.0, Kpack=None, cap=True):
    N = len(nodes)
    K, dof, R, kl, L = Kpack if Kpack is not None else assemble(nodes, seg, lam_seg, bend, cap)
    u = np.full(3 * N, np.nan)
    if axis == 'x':
        lo, hi, comp, other = bn[-4], bn[-2], 0, 1
    else:
        lo, hi, comp, other = bn[-1], bn[-3], 1, 0
    if len(lo) == 0 or len(hi) == 0:
        return np.nan, None, None
    u[3 * np.array(lo) + comp] = 0.0
    u[3 * np.array(hi) + comp] = -1.0
    u[3 * lo[0] + other] = 0.0
    fixed = ~np.isnan(u); free = ~fixed
    Kff = K[free][:, free].tocsc()
    reg = 1e-14 * abs(Kff.diagonal()).max()
    Kff = Kff + reg * sp.eye(Kff.shape[0], format='csc')
    rhs = -(K[free][:, fixed] @ u[fixed])
    u[free] = spl.spsolve(Kff, rhs)
    ue = u[dof]                                   # (S,6) global
    ul = np.einsum('nij,nj->ni', R, ue)           # local
    ua = np.zeros_like(ul); ua[:, [0, 3]] = ul[:, [0, 3]]
    ub = ul - ua
    Ua = 0.5 * np.einsum('ni,nij,nj->n', ua, kl, ua)
    Ub = 0.5 * np.einsum('ni,nij,nj->n', ub, kl, ub)
    E = 2 * (Ua.sum() + Ub.sum())
    return E, Ua, Ub
