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

## Analytical picture (2026-11-14, `twoterm.py`, shell3d bend flag)
1. omega1 = weakest-span flexural mode (parameter-free; corr 0.979; Bloch edge = 0.75 omega1).
2. Index ceiling n <= n_max = (rho/C_pin)^1/2, C_pin with bending off (pinned 2D / membranes 3D).
   3D plates: n_max 1.731 (SR) vs n 1.709 (rho .05), 1.624 (rho .2); balanced 1.580 vs 1.536, 1.454.
   2D pinned framed tilings: n_max_x SR 4.75-8.09 (seed), centroid 5.4-7.9, balanced 5.60-5.70,
   regular 5.73-5.78, uniform 3.3-4.7 (results/twoterm.json).
3. n^-2 = s + c (lam/l)^q with s = n_max^-2: measured q ≈ 0.7-1.1 (not 2: stretching-bending
   coupling, sub-isostatic regime); balanced c ≈ 6x SR at rho 0.4 -> balanced needs low rho for
   high n -> near ceiling -> small band. Fits omega1(n) tight for balanced/regular (scatter
   0.01-0.09), loose for SR/centroid (0.3). omega1 ~ rho^0.79 (SR), rho^0.51 (balanced).
   Crossover: scaling argument only, no closed form.

## Dense Bloch sweep, final (2026-11-14, `blochsweep2.py`, 192 runs, 188 resolved)
Direction-resolved interpolation at fixed n (mean ln omega_h, SE 15-25%):
  n=3.0: balanced 0.65, centroid 0.65, SR 0.37 | n=3.5: centroid 0.50, balanced 0.35, SR 0.22
  n=4.0: centroid 0.33, balanced 0.25, SR 0.15 | n=5.0: SR 0.22, centroid 0.20, balanced 0.11
  n=5.5: centroid 0.23, SR 0.17, balanced 0.09. Coverage: balanced 2.06-6.3, SR/centroid 2.7-7.8.
=> balanced below n≈2.7 (and ties centroid at 3); centroid best from 3.5 up; SR > balanced only
above n≈4.5, ≈ centroid near 5. Per-density averaging over x/y (earlier fig) was misleading.
omega_h vs weakest span: corr 0.97/0.94/0.97 (SR/centroid/balanced); prefactor 0.63/0.79/0.31
(balanced below the pinned limit: uniform spans resonate collectively).
Manuscript updated (abstract 494 chars, ~2270 words, 4 pp): random angles give index,
centroid cut removes slender spans; mixing rules (fish-eye, waveguide) with impedance caveat.

## Fish-eye device test (2026-11-14, `device/fisheye.py`, R = 40 cells) — NOT for the paper
Homogenised continuum of rule cells (isotropised K, G, rho from psweep), point source on the rim,
image gain at the antipode. Gains non-monotone 0.2-5 for all designs (interference, lens only
3.5-4 background wavelengths). Mixed (balanced rim / centroid core, 2:1) max in-band gain 2.0;
balanced alone (clipped 1.66:1) 5.3 at R/lambda0 = 3. Impedance step at the rule switch
(balanced rho≈0.13 vs centroid ≈0.38 at n = 3.5) outweighs the contrast. Band edges 0.20-0.25 c0/L.
User: fish-eye is not a real application of division lattices (existing fish-eyes: graded plates,
Lefebvre 2015; microwave cells, Ma 2012). Searching for a real high-contrast lattice application.

## Real high-contrast lattice application? (search 2026-11-14) — none found
No elastic/acoustic lattice device mixes cell families to exceed ~1.6x index contrast. Built lattice
GRIN devices: one family, n <= 2 (Luneburg 1.41; Su et al. JASA 2017 pentamode ~2x, unverified).
High-n devices avoided: Fu et al., Phys. Rev. Materials 2, 105202 (2018), mirrored Luneburg
retroreflector chosen because the Eaton lens needs "extremely high" index (to verify before citing).
Mixing families practised in microwave GRIN (Qi, Tang & Cui, Sci. Rep. 5, 11732 (2015): 3 cell kinds,
to cross the substrate index) and in graded lattices for stiffness (Al-Ketan et al., JMBBM 102,
103520 (2020)). Unchecked: 2025 MSSP SC/BCC/ABH graded lattice (S0888327025014074).
Ma et al. NJP 14, 025001 (2012) RETRACTED (NJP 18, 059701 (2016)) — removed from refs.bib.

## Cut-angle grading (2026-11-14, `anglegrade.py`, rule (area, centroid, ('perpmix', kappa)))
kappa = 0 balanced ... 1 centroid; 4 seeds; framed periodic, isotropised n:
  rho 0.1: 3.98 -> 5.45 (1.37x); rho 0.2: 2.94 -> 4.41 (1.50x); rho 0.4: 2.37 -> 3.57 (1.51x).
  Monotone in kappa, seed spread 0.01-0.18. kappa x rho: 2.37-5.45 = 2.3x in ONE family
  (vs 2.45x by mixing rules) -> no rule switch, no impedance jump. Isotropy peaks at
  intermediate kappa (eta 0.94 at kappa 0.4, rho 0.4). Bloch edge at rho 0.2: 1.00, 0.87, 0.77,
  0.77, 0.61, 0.29 for kappa 0 ... 1. Replaces the mixing material in the paper.

## Shock routing (2026-11-14, `device/shock.py`) — Fig. 4 of the manuscript
Homogenised 60x30-cell insert on rigid base, half-sine pulse tau = 20 L/c0 above a component
(|x| < 4); cut-angle grading at rho = 0.2 (n 2.94-4.41), polymer loss 0.03; <1% of the pulse
spectrum above the band edge. Diverging grading (low-n column, high-n within X of the axis):
peak force 0.96 (X 25) ... 0.54 (X 5) ... 0.49 (X 3) of uniform; impulse 0.69-0.79. Converging
control: impulse 1.2-1.8x. Uniform high-n slab: arrival +50%, peak +16% (tau 40). X < 5 is beyond
safe homogenisation (cell-size gradients). Full-lattice check not done.
Vibration (ABH) trap (`device/abh.py`): grading +30-60% absorption, capped by tip band edge; not used.
Market/feasibility scan: helmets no (ms pulses, liner << wavelength); electronics drop/shock inserts
and blast pads plausible; prior art Radi et al. TAML 14, 100553 (2024). Manuscript now 5 pp.

## Full-lattice check of Fig. 4 FAILED the homogenised model (2026-11-14, `device/fulllattice.py`)
Reduced device (40x20 cells, T = 100 per cell, shared grid-line frames): full lattice U/D/C peak
1.00/0.71/0.73, impulse 1.00/0.66/0.86; homogenised (isotropised framed cells) 1.00/1.21/0.99,
1.00/1.20/0.88; absolute U peak 1.38 (full) vs 0.59 (hom). Cause: the shared square frame forms
straight vertical/horizontal load paths (strongly anisotropic), lost by isotropising K, G.
Fig. 4 (homogenised, framed T = 1000 tables) is therefore not reliable.
Replacement: `device/globalshock.py` — ONE division of the whole insert (rectangular domain,
position-dependent cut-angle tilt 'perpmix_field'), outer skin only, full lattice, no homogenisation.
Also fixed: Fig. 4(c) impulse was the tensile rebound; now compressive impulse.

## Global full-lattice shock result (2026-11-14, `device/globalshock.py`, 2 realisations)
40x20 insert, 80k ligaments, one division, skin only, tau = 80. Relative to uniform kappa = 0:
  D (diverging, X = 5): peak 0.79/0.94, impulse 0.82/0.82
  C (converging):       peak 1.21/1.09, impulse 1.10/0.95
  S (uniform kappa = 1): peak 0.84/0.92, impulse 0.84/0.83
=> D ≈ S: the reduction comes from the softer/slower material, not routing; directional effect
(D vs C) comparable to realisation scatter. Homogenised Fig. 4 (peak halved) not confirmed.
Recommendation: drop the application section and Fig. 4; keep cut-angle n-tuning as physics.

## Regular article: network mappings per rule (2026-11-15, `network/`)
User: drop the application; characterise the networks behind each Cowan division; for each structure
network mapping / physics / relation; network matters only where it explains a physical difference.
- `network/netstats.py` (T=2e4 asym, T=1e3 phys): cell tree (split fraction, depth, areas), corners,
  elongation, dual graph (degree, C, r, paths), frame spans, fabric, junction angle, host age.
- SAME for all rules: <corners>=4, <k>->4, K^2P(k>=K) flat 3.2-6.1 (gamma=3, b=1/2), h/t ~ U^2
  (moments within 0.06), <d> 7.1-8.3. => cannot explain physics.
- DIFFERENT: split law (uniform for SR/uniform/rect; [4/9,1/2] for centroid cuts = Grunbaum, min 0.4444),
  leaf-depth var 15-16 vs 0.9 vs 0.23, triangles 30/28/0/20/5/4 %, C .38/.45/0/.28/.06/.05,
  elongation e95 9.4/10.6/918/6.2/2.0/1.9, longest span*sqrt(T) 8.1/41/127/6.6/3.1/2.2 (T=2e4).
- KEY RELATION (`network/collapse.py`, results/collapse.json): with l_* = (sum l^11/sum l)^(1/10)
  (frame spans), n^-2 = s + c (lam/l_*)^q fits every frame rule + kappa family at every density with
  one (s=0.023, q=1.12): rms ln n 0.10 vs 0.087 per-rule fits; rule offsets <= 0.06 (mean span: 0.24,
  none: 0.28). uniform excluded (89% spans < 2 lam at rho .2: frame model invalid).
- Within-rule (`network/seedvar.py`, 16 seeds): n spread follows total length (mean span) for
  balanced/regular (corr .94/.97); SR/centroid spread too small.
- Quasi-crossings (`network/crossing.py`, `breakx.py`): correlate with low n but moving them changes
  n < 0.01 (moves ~ thickness) -> not used.
- Isotropy: kappa family eta peaks at 0.4 where fabric a4 -> 0 while cells compact; eta_aff .81/.97/.89.
  rect anisotropy = fabric (eta_aff .60).
- omega_h * l_max at fixed n NOT universal: centroid ~1.6x SR; balanced ~ SR (collective prefactor).
- Bloch for regular, rect running (`blochsweep2.py regular rect`, log results/bloch2b.log).
- Manuscript rewritten as regular article (revtex pre), 9 pp, figs 1-5; tbd: regular/rect Bloch.

## Potential applications beyond GRIN (discussion 2026-11-15; logged on user request, untested unless stated)
- Stiffness per weight: n^-2 = C11/(rho E0) is a normalised specific modulus. balanced/regular ~2.2x
  stiffer than SR at rho=0.2 (set by l_*) -> lightweight cores; 3D plates near-affine.
  SR/centroid/rect bending-dominated, soft at equal weight -> cushioning, packaging, impact absorption;
  rect: anisotropic compliance.
- Vibration: rect/uniform very low local resonances (rect w1 ~100x lower) -> low-frequency isolation;
  SR broad span spectrum -> broadband damping with lossy constituent (untested); balanced uniform spans
  -> narrow clean local-resonance gap (untested).
- Damage tolerance: path multiplicity <ln mu> 2.8-3.6 (rect, balanced) vs 1.9 (SR) -> redundant load
  paths; test with rg/failure.py (buckling cascades) across rules.
- Transport: conduction through ligaments (rg/conductance.py, exact); crack-templated transparent
  metal-mesh electrodes are T-junction divisions -> rule choice may matter for uniformity/transparency.
  Pore space: scale-rich pore sizes (SR, uniform) -> multiscale filtration, wicking, broadband porous
  absorbers; regular -> narrow pore window (scaffolds, catalyst supports).
- Biology: Errera's rule (shortest wall through the centre) ~ balanced; balanced maximises stiffness per
  weight and isotropy with compact cells. Wing venation and crack patterns elsewhere in the family.
- Proposed: "Beyond wave control" paragraph (stiffness reading + Errera); optional quick checks of
  conductance and buckling across the six rules (~1 h each), after the Bloch runs.

## Bloch for regular and rect (2026-11-15, results/bloch2_T1000.json, results/bloch2_rect.json)
- regular omega_h at n = 2.5/3/3.5/4: 0.89/0.42/0.27/0.16, below balanced (1.08/0.65/0.35/0.25):
  most uniform spans -> strongest collective resonance; most uniform division is not the best.
- rect: n from 1.7 to 12.8 across realisations and x/y at the same densities; omega_h < 0.1 everywhere,
  at the k-grid resolution limit (identical values across seeds) -> upper bounds only.
- Note: the first blochsweep2 process (regular rect) went on to rect after regular, writing into
  bloch2_T1000.json; figs.py replaces those rect rows by bloch2_rect.json.

## Referee report (PRE Polymers, 2026-11-15; ~/Downloads/pre_polymers_review.md) — response ~/Downloads/response_pre_polymers.md
User: disagree with point 1 (collapse holds); address the rest, hold back where correct.
New scripts/results (network/): collapse_checks.py (LOO 0.057-0.085, rms(p) saturates p>=6,
per-rule pinned s fixed rms .055 vs common .056), collapse refit WITHOUT rect: n_max 7.8 (= direction-
averaged pinned 6.7-8.3), q .97, c 1.2, rms .060 (per-rule .047, mean span .131, lambda only .152);
rect half within 15%, no common ceiling (pinned <n> 5.7-10.4). kinematic_eta.py: single-cell isotropy
keeps the kappa maximum (.77/.90/.90/.88/.75/.52); SR/centroid drop to .56/.55 (off-axis .3-.4).
kernel.py: slope b .47-.50 all but regular (non-stationary); balanced flat k<=4. tsize.py: n rises with
T (balanced 2.68/2.97/3.10 at T 300/1e3/3e3), ranking unchanged, law predicts T=3000 rms .05-.10.
beamcheck.py: Timoshenko +2-6%, +rigid joints -10-12%; balanced/SR ratio .64-.68.
Bloch: 8 realisations SR/centroid/balanced (results/bloch2_more_*.json, bloch_fixed_n.json);
grid step 15-18% of omega_h; zone folding not setting omega_h (20/247 near pi/2n).
Error found: rect band edges saturated in half the scans (k_h = k_max = 3 omega1): lower bounds only.
Pinned regulariser 1e-7..1e-11: no change; z = 3, 1/2 floppy mode per junction.
Manuscript 10 pp. Title kept. New refs (from memory, please double-check): heussinger2006,
broedersz2014, silva1995, zhu2000, graff1975, connelly1993.
