"""Referee 6: is the low index of the balanced divisions an Euler-Bernoulli artefact at high density?
Variants of the frame element: EB (default; spans shorter than lambda capped at lambda);
TIM: Timoshenko shear (nu0 = 0.3, k = 5/6), softer; TIM+RJ: Timoshenko with rigid joint zones,
approximated by a flexible length L - lambda (floor 0.5 lambda), stiffer.
Framed periodic tilings, T = 1e3, alpha = 0, rho = 0.2 and 0.4, seeds 100-103. -> results/beamcheck.json"""
import json, os, sys
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); B = os.path.join(HERE, '..')
sys.path[:0] = [B, os.path.join(B, '..', 'rg')]
import divisions, periodic, waves, frame
_orig = frame.element_k
MODE = {'m': 'EB'}

def element_k(L, c, s, EA, EI):
    if MODE['m'] == 'EB':
        return _orig(L, c, s, EA, EI)
    lam = np.sqrt(12 * EI / np.maximum(EA, 1e-300))
    Lf = np.maximum(L - lam, 0.5 * lam) if MODE['m'] == 'TIM+RJ' else L
    G = 1 / (2 * 1.3); phi = 12 * EI / ((5 / 6) * G * EA * Lf ** 2)
    n = len(L); kl = np.zeros((n, 6, 6))
    a = EA / L; f = EI / (1 + phi)
    b = 12 * f / Lf ** 3; cc = 6 * f / Lf ** 2; d = (4 + phi) * f / Lf; e = (2 - phi) * f / Lf
    kl[:, 0, 0] = kl[:, 3, 3] = a; kl[:, 0, 3] = kl[:, 3, 0] = -a
    kl[:, 1, 1] = kl[:, 4, 4] = b; kl[:, 1, 4] = kl[:, 4, 1] = -b
    kl[:, 1, 2] = kl[:, 2, 1] = kl[:, 1, 5] = kl[:, 5, 1] = cc
    kl[:, 4, 2] = kl[:, 2, 4] = kl[:, 4, 5] = kl[:, 5, 4] = -cc
    kl[:, 2, 2] = kl[:, 5, 5] = d; kl[:, 2, 5] = kl[:, 5, 2] = e
    R = np.zeros((n, 6, 6))
    for o in (0, 3):
        R[:, o, o] = c; R[:, o, o + 1] = s; R[:, o + 1, o] = -s; R[:, o + 1, o + 1] = c; R[:, o + 2, o + 2] = 1
    return np.einsum('nji,njk,nkl->nil', R, kl, R), R, kl
frame.element_k = element_k

out = []
for rule in ('SR', 'centroid', 'balanced', 'regular'):
    for s in range(4):
        g = divisions.run(1000, 100 + s, rule); L = np.linalg.norm(g['B'] - g['A'], axis=1)
        for rho in (0.2, 0.4):
            lam = np.full(1000, rho / L.sum()); row = dict(rule=rule, seed=s, rho_t=rho)
            for m in ('EB', 'TIM', 'TIM+RJ'):
                MODE['m'] = m
                sc = periodic.build(g, lam, frame_lam=lam[0]); C = periodic.periodic_C(sc)
                th, n = waves.christoffel_n(C, sc['rho']); row[m] = float(n.mean())
            out.append(row); print(rule, s, rho, ' '.join('%s %.2f' % (m, row[m]) for m in ('EB', 'TIM', 'TIM+RJ')), flush=True)
json.dump(out, open(os.path.join(B, 'results', 'beamcheck.json'), 'w'), indent=1)
for rule in ('SR', 'centroid', 'balanced', 'regular'):
    for rho in (0.2, 0.4):
        v = [o for o in out if o['rule'] == rule and o['rho_t'] == rho]
        print('%-9s rho %.1f  EB %.2f  TIM %.2f  TIM+RJ %.2f' % (rule, rho, *[np.mean([o[m] for o in v]) for m in ('EB', 'TIM', 'TIM+RJ')]))
