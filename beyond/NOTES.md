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
- D = 3 not done (d3/beams3d.py is the starting point).
- Not tried: rules designed for the GRIN merit directly (e.g. balanced backbone +
  painted thickness by depth), or density-at-fixed-n optimisation.
