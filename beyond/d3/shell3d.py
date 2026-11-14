"""Linear thin-shell model of a plate division; long-wavelength waves.

Plates of thickness lam_t (E0 = rho0 = 1, nu0 = 0.3), meshed by mesh3d.
  membrane: CST plane-stress triangles in each plate's plane;
  bending:  discrete hinges (linearised dihedral, stiffness 3|e|^2/(A1+A2) * D,
            D = E lam^3 / 12(1-nu^2), lam = min of the pair) for every pair of
            triangles sharing a mesh edge -- within a plate (flat) and across a
            T-junction line (rigid joint).
Translational dofs only. Effective C (6x6 Voigt, engineering shear) from six
kinematic tests (all nodes on the cube faces displaced affinely), one factorisation.
Waves: Christoffel over the sphere, n(dir) = c0/v_qL with c0 = sqrt(E0/rho0) = 1;
n_x along x (the paper's n_eff); eta_n = min n / max n.
omega1: lowest resonance with the cube faces clamped (lumped mass lam*A/3).
"""
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spl

NU = 0.3
VOIGT = {(0, 0): 0, (1, 1): 1, (2, 2): 2, (1, 2): 3, (2, 1): 3, (0, 2): 4, (2, 0): 4, (0, 1): 5, (1, 0): 5}


def assemble(X, tri, lam_tri):
    N = len(X); M = len(tri)
    x1, x2, x3 = X[tri[:, 0]], X[tri[:, 1]], X[tri[:, 2]]
    e1 = x2 - x1; e1 /= np.linalg.norm(e1, axis=1)[:, None]
    nn = np.cross(x2 - x1, x3 - x1); A2 = np.linalg.norm(nn, axis=1); nn /= A2[:, None]
    e2 = np.cross(nn, e1)
    A = 0.5 * A2
    loc = np.stack([np.c_[np.einsum('ij,ij->i', x - x1, e1), np.einsum('ij,ij->i', x - x1, e2)] for x in (x1, x2, x3)], 1)
    xs, ys = loc[..., 0], loc[..., 1]
    y23, y31, y12 = ys[:, 1] - ys[:, 2], ys[:, 2] - ys[:, 0], ys[:, 0] - ys[:, 1]
    x32, x13, x21 = xs[:, 2] - xs[:, 1], xs[:, 0] - xs[:, 2], xs[:, 1] - xs[:, 0]
    Bm = np.zeros((M, 3, 6))
    Bm[:, 0, [0, 2, 4]] = np.c_[y23, y31, y12]
    Bm[:, 1, [1, 3, 5]] = np.c_[x32, x13, x21]
    Bm[:, 2, [0, 2, 4]] = np.c_[x32, x13, x21]
    Bm[:, 2, [1, 3, 5]] = np.c_[y23, y31, y12]
    Bm /= (2 * A)[:, None, None]
    D = np.array([[1, NU, 0], [NU, 1, 0], [0, 0, (1 - NU) / 2]]) / (1 - NU ** 2)
    Kl = (lam_tri * A)[:, None, None] * np.einsum('mki,kl,mlj->mij', Bm, D, Bm)
    Tm = np.zeros((M, 6, 9))
    for a in range(3):
        Tm[:, 2 * a, 3 * a:3 * a + 3] = e1
        Tm[:, 2 * a + 1, 3 * a:3 * a + 3] = e2
    Kg = np.einsum('mki,mkl,mlj->mij', Tm, Kl, Tm)
    dof = (3 * tri[:, :, None] + np.arange(3)).reshape(M, 9)
    rows = [np.repeat(dof, 9, axis=1).ravel()]; cols = [np.tile(dof, (1, 9)).ravel()]; vals = [Kg.ravel()]
    # hinges
    E = np.vstack([tri[:, [0, 1, 2]], tri[:, [1, 2, 0]], tri[:, [2, 0, 1]]])   # (edge a, edge b, apex)
    tid = np.tile(np.arange(M), 3)
    key = np.sort(E[:, :2], axis=1)
    order = np.lexsort((key[:, 1], key[:, 0]))
    key, E, tid = key[order], E[order], tid[order]
    brk = np.r_[0, np.nonzero(np.any(np.diff(key, axis=0) != 0, axis=1))[0] + 1, len(key)]
    P1, P2, P3, P4, TA, TB = [], [], [], [], [], []
    for s, e in zip(brk[:-1], brk[1:]):
        if e - s < 2:
            continue
        for i in range(s, e):
            for j in range(i + 1, e):
                P1.append(key[i, 0]); P2.append(key[i, 1]); P3.append(E[i, 2]); P4.append(E[j, 2])
                TA.append(tid[i]); TB.append(tid[j])
    P1, P2, P3, P4, TA, TB = map(np.array, (P1, P2, P3, P4, TA, TB))
    xa, xb = X[P1], X[P2]
    ev = xb - xa; Le = np.linalg.norm(ev, axis=1); ah = ev / Le[:, None]
    G = []
    for xp in (X[P3], X[P4]):
        s = np.einsum('ij,ij->i', xp - xa, ah) / Le
        dv = xp - (xa + (s * Le)[:, None] * ah)
        hh = np.linalg.norm(dv, axis=1)
        G.append((s, dv / hh[:, None], hh))
    (s3, d3, h3), (s4, d4, h4) = G
    nA = np.cross(ah, d3); nB = -np.cross(ah, d4)
    g3 = nA / h3[:, None]; g4 = nB / h4[:, None]
    g1 = -(1 - s3)[:, None] * g3 - (1 - s4)[:, None] * g4
    g2 = -s3[:, None] * g3 - s4[:, None] * g4
    g = np.concatenate([g1, g2, g3, g4], axis=1)
    lam_h = np.minimum(lam_tri[TA], lam_tri[TB])
    kh = 3 * Le ** 2 / (A[TA] + A[TB]) * lam_h ** 3 / (12 * (1 - NU ** 2))
    hd = np.concatenate([3 * P1[:, None] + np.arange(3), 3 * P2[:, None] + np.arange(3),
                         3 * P3[:, None] + np.arange(3), 3 * P4[:, None] + np.arange(3)], axis=1)
    Kh = kh[:, None, None] * g[:, :, None] * g[:, None, :]
    rows.append(np.repeat(hd, 12, axis=1).ravel()); cols.append(np.tile(hd, (1, 12)).ravel()); vals.append(Kh.ravel())
    K = sp.coo_matrix((np.concatenate(vals), (np.concatenate(rows), np.concatenate(cols))), shape=(3 * N, 3 * N)).tocsr()
    return K, A


def strain_fields(X):
    x = X - 0.5
    U = []
    for I in range(6):
        eps = np.zeros((3, 3))
        for (i, j), k in VOIGT.items():
            if k == I:
                eps[i, j] = 1.0 if i == j else 0.5
        U.append((x @ eps.T).ravel())
    return np.array(U)


def christoffel(Cv, rho, ndir=600):
    k = np.arange(ndir) + 0.5
    z = 1 - k / ndir                        # upper hemisphere (n and -n equivalent)
    phi = np.pi * (1 + 5 ** 0.5) * k
    r = np.sqrt(1 - z * z)
    dirs = np.c_[r * np.cos(phi), r * np.sin(phi), z]
    dirs = np.vstack([[1, 0, 0], dirs])
    C4 = np.zeros((3, 3, 3, 3))
    for (i, j), I in VOIGT.items():
        for (k_, l), J in VOIGT.items():
            C4[i, j, k_, l] = Cv[I, J]
    G = np.einsum('ijkl,dj,dl->dik', C4, dirs, dirs)
    lmax = np.linalg.eigvalsh(G)[:, -1]
    return np.sqrt(rho / lmax)


def analyse(X, tri, lam_tri, resonance=True):
    K, A = assemble(X, tri, lam_tri)
    N = len(X)
    bnd = np.any((X < 1e-9) | (X > 1 - 1e-9), axis=1)
    fixed = np.repeat(bnd, 3); free = ~fixed
    Kff = K[free][:, free].tocsc()
    Kff = Kff + 1e-12 * abs(Kff.diagonal()).max() * sp.eye(Kff.shape[0], format='csc')
    lu = spl.splu(Kff)
    Ua = strain_fields(X)
    Kfx = K[free][:, fixed]
    U = Ua.copy()
    for I in range(6):
        U[I, free] = lu.solve(-(Kfx @ Ua[I, fixed]))
    C = U @ (K @ U.T); C = 0.5 * (C + C.T)
    Caff = Ua @ (K @ Ua.T); Caff = 0.5 * (Caff + Caff.T)
    rho = float((lam_tri * A).sum())
    n = christoffel(C, rho); na = christoffel(Caff, rho)
    out = dict(rho=rho, n_x=float(n[0]), n_mean=float(n[1:].mean()), n_min=float(n.min()), n_max=float(n.max()),
               eta_n=float(n.min() / n.max()), n_aff_mean=float(na[1:].mean()), C=C.tolist(), Caff=Caff.tolist(),
               N=N, M=len(tri))
    if resonance:
        m = np.zeros(N)
        for a in range(3):
            np.add.at(m, tri[:, a], lam_tri * A / 3)
        Mf = sp.diags(np.repeat(m, 3)[free]).tocsc()
        w2 = spl.eigsh(Kff, k=1, M=Mf, sigma=0, which='LM', return_eigenvectors=False)
        out['omega1'] = float(np.sqrt(w2.min()))
    return out


def thickness(g, tplate, A_tri, rho, alpha):
    T = g['T']
    base = np.arange(1, T + 1.) ** -alpha
    area_p = np.bincount(tplate, weights=A_tri, minlength=T)
    lam = base * rho / np.sum(base * area_p)
    return lam, area_p
