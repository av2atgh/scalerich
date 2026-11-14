# beyond/ — are SR the best recursive divisions for elastic-wave control?

Started 2026-11-14. Question: Both et al. claim "programmable elastic-wave
control" for SR (Cowan iterated division: cell ∝ area, uniform point, uniform
angle; thickness λ0 t^-α). Can other recursive divisions of the square do better?

## Reduction
The paper works in the long-wavelength regime (Bloch analysis of RVEs, GRIN lens at
f well below pixel resonances). There the wave response is fixed by (C_eff, ρ̄):
n(θ) = c0 / v_qL(θ), v_qL² = λ_max Γ(θ) / ρ̄ (Christoffel). Paper's n_eff = n(0).
Three figures of merit for a GRIN designer:
1. range of n at fixed density window,
2. isotropy of n (η_n = min n / max n) and the yield of quasi-isotropic RVEs,
3. bandwidth of the homogenised regime: lowest local resonance ω1 (boundary clamped).
Bounds: n ≥ 1 in every direction (Voigt); no upper bound (n → ∞ as C → 0), so
"widest range" alone is not a well-posed optimum: anything approaching a mechanism wins it.
Longitudinal waves see C11 ~ K + G, not E: honeycomb (stretch-dominated K) stays at n ≈ 1.9.

## Code
- `divisions.py` — rule = (pick, where, angle); SR, uniform-pick, rect, centroid,
  balanced (centroid, ⟂ major axis), regular (largest first, centroid, ⟂ major axis).
- `waves.py` — rg/frame.py beams, kinematic affine BCs on all boundary nodes → C;
  Christoffel n(θ); lowest clamped resonance ω1 (lumped mass), units c0/L_RVE.
- `lattices.py` — square, triangular, honeycomb references, same model.
- `sweep.py` (T = 1000, 4 seeds, α ∈ {0,…,1}, ρ̄ ∈ {0.1,…,0.4}) → results/sweep_T1000.json,
  `summarize.py` → results/summary_T1000.txt; `tcheck.py` (T = 1000 vs 3000) → results/tcheck.json.

## Results (2026-11-14)
- Calibration: SR n_x spans 1.66–7.60 over ρ̄ = 0.1–0.4, α = 0–1 (paper: 1.5–8.3).
  Square n_x = 1.41, triangular 1.71, honeycomb 1.9–2.1 (paper: 1.0–1.9 for square/hex).
- Density, not α, carries the n range. SR at ρ̄ = 0.2: n_x = 4.12, 3.06, 2.60, 2.56, 2.68
  for α = 0…1 (non-monotone). α mainly lowers ω1 (0.44 → 0.02) and raises anisotropy.
- Per-sample n_x ranges: SR 1.66–7.60, uniform 1.81–6.43, centroid 1.81–7.03,
  balanced 1.84–4.58, regular 1.77–3.53, rect 1.19–18.4 (n_max at 45°, η_n ≈ 0.27,
  ω1 ≤ 0.05: a near-mechanism in shear).
- Isotropic RVEs (η_n ≥ 0.8): SR 9/80 with n 3.06–4.66; balanced 16/80 with n 2.38–4.77;
  uniform-pick and centroid 0/80.
- Equal n, best ω1 among η_n ≥ 0.7 samples: n ∈ [2.5,3): balanced 1.39 vs SR none;
  [3,3.5): balanced 1.19 vs SR 0.95; [4,5): balanced 0.76 vs SR 0.59.
  Matched means: SR (ρ̄ 0.2, α 0) n 4.36, η 0.56, ω1 0.44; balanced (ρ̄ 0.1, α 0)
  n 4.28, η 0.77, ω1 0.74.
- More scales is worse: uniform-pick (broader cell sizes than SR) has lower n range,
  η_n ≈ 0.4–0.6 and ω1 3–5× below SR.
- T = 3000 check: n drifts up 5–10% for all rules, ranking unchanged; ω1/√T falls
  faster for SR than for balanced (longest slender members grow with the hierarchy).

## Reading
SR's advantage is n per unit density: at fixed ρ̄ it is softer in C11 than a balanced
division. At fixed n (what a GRIN lens needs) the balanced division of Cowan's
family is more isotropic, yields more isotropic RVEs, and has 1.3–3× the
homogenisation bandwidth; it reaches the same isotropic n window (up to ≈ 4.8) by
going to lower density. The bandwidth cost of SR is the breadth of its length
distribution: the first local resonance is set by the longest slender member,
ω1 ~ λ/ℓ_max², while n is set by the typical load-bearing slenderness.

## Caveats / next
- 2D Euler–Bernoulli beams; kinematic BCs (upper bound on C, so n is a lower
  estimate; boundary layer ~ largest cell, larger for SR). Periodic RVEs would be cleaner.
- Long-wavelength only: no band structure, no band gaps (the other "wave control" route).
- Not tried: rules designed for the GRIN merit directly (e.g. balanced backbone +
  painted thickness by depth), or density-at-fixed-n optimisation.

## 3D: plate divisions of the cube (2026-11-14, `d3/`)
- `d3/div3d.py` — same six rules on convex polyhedra (normals iso / axes / major axis);
  SR checks: m̄ = 4.57, total area 20.9 at T = 300 (d3/NOTES: 4.58, 3.3 T^{1/3}).
- `d3/mesh3d.py` — conforming triangle mesh across T-junction lines (Triangle 'pYY').
- `d3/shell3d.py` — CST membranes + discrete hinges (rigid joints), ν0 = 0.3, 6x6 C from
  kinematic tests, Christoffel on the sphere, clamped ω1. Rigid-body energy 1e-13;
  h = 0.07/0.05/0.035 → n within 1-2% (`d3/check.py`).
- `d3/sweep3d.py` T = 300, α ∈ {0, .5, 1}, ρ̄ ∈ {0.05, 0.1, 0.2}, 3 seeds → results/sweep3d_T300.json.
- 'uniform' invalid in 3D: plates tiny (total area ≈ 2), λ/√a up to 157 → excluded.
- Near-affine for every rule: n_aff = 1.30-1.32 (isotropic normals); C11/C11_aff = 0.5-0.67
  at α = 0, 0.9-1 at α = 1. SR ⟨n⟩ 1.36-1.73 over the scan; all valid rules 1.32-1.84;
  rect n_x up to 3.1 (η 0.44). SR at α = 0: n 1.71 → 1.62 for ρ̄ 0.05 → 0.2 (stretch-dominated).
- Isotropic (η ≥ 0.8): SR 0/27, centroid 2/27, balanced 9/27, regular 15/27.
- Equal n ≈ 1.5: SR (α .5, ρ̄ .1) η .64, ω1 0.20; balanced (α 0, ρ̄ .1) η .86, ω1 1.30.
- Hinge joints add energy under affine shear, so C can exceed the membrane Voigt value
  slightly (rect α = 1: n_x min 0.89 < 0.954).

## Length-spread check (2D, `lengthspread.py`)
n·ω1 at α = 0 ordered by G = ⟨ℓ⟩_w/ℓ_max² for isotropic-angle rules (uniform < SR <
centroid < balanced), saturating for regular; no collapse (ratio drifts with ρ̄).

## Manuscript (`manuscript/manuscript.tex`, 2026-11-14; formatted for PRL)
revtex4-2 prl, pdflatex+bibtex, figs.py → fig1-3.pdf. 4 pp, ~1600 words, abstract 518 chars.
Message: density (bending) sets the n range, not scale richness; at equal n a balanced
division wins isotropy and bandwidth; in 3D plates are near-affine and n collapses to 1.3-1.8.

## Span-mode prediction of ω1 (`spanmode.py`, 2026-11-14)
ω_span = min over spans of 4.730² λ_s/(√12 ℓ_s²) (clamped-clamped, no free parameter).
corr(ln ω1, ln ω_span) = 0.979 over 480 2D cells; ω1/ω_span median 0.48, IQR 0.44-0.80;
per rule 0.44 (uniform, rect, regular) to 0.70 (SR, centroid). 0.44 = (π/4.730)², the pinned limit.

## Periodic BCs, Bloch, 3D checks, lens (2026-11-14 evening)
- Bare mirror tiling (no frame) is an artefact: ligament ends meet their images at bare kinks
  (SR n_x up to 33, rect 189). Physical reference = framed periodic tiling (square frame =
  generation 0, thickness lambda0, shared between tiles): `periodic.build(g, lam, frame_lam)`.
- Framed periodic / kinematic n_mean: median 0.92 (10-90%: 0.82-1.00); `psweep.py` -> psweep_T1000.json.
  Periodic n_x ranges: SR 1.70-6.76, centroid 1.75-5.83, balanced 1.68-4.04, regular 1.64-3.28.
  eta>=0.8: SR 9/80, balanced 12/80. Voronoi reference n ≈ 2.1, eta ≈ 0.95 (results/voronoi.json).
- Bloch (`blochsweep.py`, 2 seeds, alpha 0): band edge omega_h vs clamped omega1: corr(ln) 0.82,
  omega_h/omega1 median 0.75. Index-bandwidth curves CROSS: balanced wins at low n (n≈2.7:
  1.04 vs SR 0.45), SR at high n (n≈4.5: SR 0.49 at rho 0.4 vs balanced 0.11 at rho 0.05).
  Index from architecture (SR, high density) keeps bandwidth; index from slenderness loses it.
  Dense sweep `blochsweep2.py` (SR/centroid/balanced, 4 seeds, 8 rho, x and y) running.
- 3D plates periodic (symmetry BCs on framed cell, `periodic3d.symmetric_normal`) vs kinematic:
  within 1-3% at alpha 0, 5-15% lower at alpha 0.5. T = 600 vs 300: n within 1%.
- 3D open-cell skeleton (`skeleton3d.py`, `skelsweep.py`): n 2.3-3.7 for every rule, weak rho
  dependence (near-affine); at equal n balanced/regular omega1 ≈ 3.4x SR. No high-n window in 3D.
- Lens (`lens.py`, `lens_calib.py`): at R = 4 tiles of T = 200, both lenses sit outside the
  homogenised band (D/lambda ≈ 0.9 at the band edge); no decisive result; paused.
- Framing decided by the user: report both regimes as they are (balanced at low n, SR at high n).
