# rg/ — the SR material as a nested chygraph; renormalisation along the splitting hierarchy

Status 2026-11-12. All numbers below come from the scripts in this folder;
results in `results/` (json + log). Single-threaded, one process at a time.

## Summary

- **New and exact.** The SR splitting hierarchy, labelled by areas, *is a
  random binary search tree* (uniform keys). A chord through a uniform point
  at any angle cuts any convex cell into area fractions that are exactly
  Uniform(0,1), and cells are chosen ∝ area. So the T+1 cell areas are uniform
  spacings (Dirichlet(1,…,1)), split depth = BST insertion depth (2 ln t,
  variance 2 ln t), no multiscaling of areas — shape plays no role. Measured
  to 3 decimals.
- **New, asymptotically exact.** A chord's ends land on existing ligaments ∝
  length, and ℓ_h ∝ h^(−1/2), so the host's birth index is h = tU²
  (U uniform). The chord-to-host thickness ratio λ_t/λ_h = U^(2α) has a
  distribution that does not depend on t: an *exactly marginal* coupling, set
  by α. ⟨(h/t)^α⟩ = 1/(1+2α), measured 0.667/0.498/0.330 at α = 0.25/0.5/1.
- **Exact generation flow for transport.** σ(T) = Σ_t ψ_t aff_t with the
  rank-one (Thevenin) increment ψ_t aff_t = g_t ΔV_t²/(1+g_t R_t). For α < 1/2
  the increments are stationary in ln t, so σ/σ_aff tends to a fixed-line value
  Φ(α) = 0.59, 0.70, 0.82 at α = 0, 0.2, 0.4. α is relevant for density
  (eigenvalue 1/2−α) and exactly marginal for the non-affine prefactor.
- **Not new (translation).** The exact chygraph contraction along the tree is
  planar nested dissection: width w ≈ 5 n^0.53, root width ~ √T, work T^1.53.
  The diamond-like surrogate (a random sequential laminate) reproduces the
  marginality and the α-trend of Φ but sits 0.2 below the network.
- **Elasticity.** Pin-jointed SR networks are floppy, so stiffness comes from
  bending. The per-ligament energy law U_i ∝ i^(−x) is reproduced
  (x = 0.78/0.83/1.04 vs the paper's 0.75/0.90/1.05). The paper's m(α) of
  eq. (3) is **not** reproduced: m is between 1 and 3 when λ0 is varied (a
  theorem on that route), and m_T = 1.7–2.5 when T is varied, with no
  divergence as α → 1/2.
- **Failure proxy.** Euler-buckling cascades on SR frames are delocalised
  (strip fraction 0.16–0.20 against a random 0.16; square grid 0.47; mean
  avalanche 2 vs 11). A thickness-shuffled control is equally delocalised, so
  the thickness hierarchy is *not* the cause, consistent with the paper using
  α = 0 for this test.

## 1. The object (`srtree.py`)

Zero-width topology: thickness paints, it does not change the network (README).
The SR step "uniform point + uniform angle, grow to the first hit" is: descend
the splitting tree to the convex cell holding the point (each internal node is
a half-plane test), then cut that cell with the chord through the point.
Checked against the direct first-hit algorithm (40 seeds, T = 400):
⟨ℓ_1⟩ 0.917 vs 0.935, ⟨ℓ_t √t⟩ 1.1648 vs 1.1634.

Chygraph reading. Complex = cell. Members = its two child cells and the
chord's segments. Atoms (boundary vertices) = the junctions on the cell's
boundary that interior chords touch. The hierarchy is the splitting tree. Every
complex has two children and on average 4 corners (exactly: each split adds 4
corners). Leaves follow P(3,4,5,6) = 0.301, 0.452, 0.198, 0.042.

## 2. The splitting tree is a random BST (`tree_stats.py`, `exact_checks.py`)

**Lemma.** K is a convex body, x is uniform in K, and θ has any distribution
independent of x. The area fraction on one side of the line through x at
angle θ is Uniform(0,1).
*Proof.* For fixed θ, the offset s of x has density ℓ(s)/|K|, where ℓ(s) is
the chord length at offset s. The fraction on one side is
F(s) = ∫_{−∞}^s ℓ/|K|, the CDF of s, so F(s) is uniform. ∎

Cells are chosen ∝ area (uniform point in the square), so the area-labelled tree
is a Markov chain with the transitions of a BST built from uniform keys in
[0,1]: the 1D interval containing the key splits at the key. Consequences:

| quantity | exact (BST) | measured |
|---|---|---|
| split fraction min(r,1−r) | uniform on [0,½], density 2 | 2.017…1.988 (10 bins) |
| μ = E_sb[−ln r], σ² | ½, ¼ | 0.5001, 0.2497 |
| d⟨depth of split at t⟩/d ln t | 1/μ = 2 | 1.983 |
| d Var(depth)/d ln t | σ²/μ³ = 2 | 1.980 |
| Σ A² at T = 2·10⁴ | 2/(T+2) = 9.999e−5 | 9.998e−5 |
| P(TA > 0.5, 1, 2, 4) | e^−x = 0.607, 0.368, 0.135, 0.018 | 0.606, 0.368, 0.136, 0.018 |
| γ_A(q), M_q ~ t^−γ | q − 1 | −0.499, 0.499, 0.998, 1.997 (q = ½, 3/2, 2, 3) |
| tree height, T = 2·10⁴ | 4.311 ln T − 1.953 ln ln T + O(1) = 38.2 + O(1) | 34.0 |

T = 3·10⁴, 12 seeds, unless stated otherwise.

*Krapivsky–Ben-Naim.* The area process is identical for rectangular cracks
(`mode='rect'`), so areas never multiscale. Multiscaling lives in the lengths.
For rectangles, the moment equation dM(a,b)/dt = [1/(a+1)+1/(b+1)−1] M(a+1,b+1)
conserves M(a,b) on ab = 1, which gives M(a,0) ~ t^−γ with
γ = (a − √(a²+4))/2. Measured γ_width (rect): −0.778, −0.612, −0.402, −0.291,
−0.229, against −0.781, −0.618, −0.414, −0.303, −0.236 (a = ½, 1, 2, 3, 4).
With isotropic angles, the perimeter moments scale ordinarily, γ_P(a) = a/2 − 1
within 0.013. The mean shape ratio (vertex-covariance eigenvalues) drifts
slowly, 0.205 → 0.180 between t = 400 and 3·10⁴; this is not resolved.

*Host-age law.* Chord t ends on ligament h with h/t → U². Measured x = ln(t/h):
mean 1.28 → 1.96, sd 0.93 → 1.84 for t from 10 to 10⁴ (Exp(2): 2, 2).
Measured ⟨(h/t)^α⟩ = 0.667, 0.498, 0.330 at α = ¼, ½, 1, against 1/(1+2α) =
0.667, 0.500, 0.333. The rectangular model gives the same values. Hence
**λ_t/λ_host = U^(2α) is stationary**: in RG language the coupling between a
generation and the frame it attaches to is exactly marginal. Its moments
⟨(λ_h/λ_t)^p⟩ = 1/(1−2pα) diverge at α = 1/(2p) (p = 3, bending of the host:
α = 1/6).

*Segment length.* The segments of ligament i at the end of the run have length
≈ 0.57–0.66 T^(−1/2), the same for every birth decade except the youngest
(0.97). Every generation is cut down to the current mesh.

*Density bookkeeping.* ⟨ℓ_1⟩ = 0.9454 (the paper's ℓ0 = 0.946);
c_ℓ = lim ⟨ℓ_t⟩√t = 1.161. ρ(t) ≈ λ0 c_ℓ t^(1/2−α)/(1/2−α) + const. The
zero-width threshold for α > ½ is λ0* = 1/Σ⟨ℓ_t⟩t^(−α) = 0.084, 0.200, 0.369,
0.617 at α = 0.6, 0.75, 1, 1.5, close to 1/(c_ℓ ζ(α+½)) = 0.081, 0.187, 0.330,
0.524. The paper's (ℓ0 ζ(2α))^(−1) gives 0.189, 0.405, 0.643, 0.879. Both
diverge at α = ½. Away from it they differ by about 2×, which overlaps (finite
width, sibling `lng/srcells.py`) may or may not account for. Not settled here.

## 3. Exact contraction = nested dissection (`nested.py`)

Message of a cell = the Dirichlet-to-Neumann (Schur complement) matrix of its
interior conductance network onto its atoms. The parent sums its two
children's messages with the chord's segments and eliminates the junctions on
the chord. At the root, the electrodes are applied.

| T | rel. error vs sparse solve | w_max (root children) | Σ w³ |
|---|---|---|---|
| 300 | 1.4e−13 | 62 | 2.8e5 |
| 1000 | 2.2e−13 | 118 | 2.0e6 |
| 3000 | 2.3e−12 | 194 | 1.3e7 |
| 10000 | 1.4e−12 | 346 | 7.7e7 |

w ≈ 5 n^0.53 for a cell holding n chords (14.5, 52, 161, 287 at n ≈ 10, 100,
10³, 3·10³). The mean width is 7–8, but the top of the tree has w ~ √T. So the
complexes are narrow only at the bottom. The work grows as T^1.53, which is
Lipton–Tarjan nested dissection with the chords as separators. Unlike the
diamond lattice, the SR hierarchy does not seal loops into bounded separators:
chords end on ancestors' sides, so the junctions on a cell's boundary grow with
everything inside it.

## 4. Transport: exact generation flow and the fixed line Φ(α)

Ligament segment conductance λ/ℓ_seg. Field along x (and y); top and bottom
free. σ_aff = Σ λℓcos²θ is the affine (Voigt) upper bound, ½ρ when isotropic.

**Identity** (`flow.py`). Adding chord t (conductance g_t = λ_t/ℓ_t between P
and Q) changes the dissipated power by g_t ΔV_t²/(1 + g_t R_t), where ΔV_t is
the pre-existing potential difference and R_t the PQ resistance with the
electrodes grounded. Checked against direct solves to 4e−7. So
σ(T) = Σ_t ψ_t aff_t exactly, with ψ_t = (field transmission F_t) × (Thevenin
factor 1/(1+g_tR_t)).

Window increments ψ = Δσ/Δσ_aff (T = 2·10⁴, 4 seeds, from t ≈ 50 on):

| α | ψ (stationary windows) | Φ(T) = σ/σ_aff at T = 10³, 3·10³, 10⁴, 3·10⁴ |
|---|---|---|
| 0 | 0.58–0.61 | 0.584, 0.595, 0.592, 0.597 |
| 0.2 | 0.67–0.72 | 0.658, 0.677, 0.683, 0.691 |
| 0.4 | 0.75–0.86 | 0.710, 0.732, 0.746, 0.757 |
| 0.5 | 0.78–0.95, rising | 0.724, 0.746, 0.762, 0.774 |
| 0.6 | 0.82–1.07, rising | 0.731, 0.751, 0.767, 0.779 |
| 1.0 | 1.05–1.80, rising | 0.690, 0.697, 0.701, 0.704 |

(Φ(T) from `conductance.py`, x and y averaged, 4 seeds.)

Reading:
- For α < ½ the flow is stationary. Each new generation adds the same fraction
  ψ(α) of its affine weight, and Φ(T) → ψ(α) with a T^(−(½−α)) transient from
  the first ligaments, which explains the slow rise at α = 0.4. This is a
  *line of fixed points* labelled by α. The chord/frame coupling g_tR_t is
  stationary (median 0.92 at α = 0, about 0.6 at 0.2).
- For α > ½ the late increments exceed 1. The old thick wires are nearly
  equipotential, so the field is squeezed into the gaps the late chords bridge
  (field focusing, F > 1). The late weight is summable, so Φ(T) converges to a
  value set by the first few ligaments and the free edges (0.70 at α = 1).
- The median g_tR_t falls only slowly for α > ½ (0.37 → 0.17 at α = 1). R_t is
  dominated by the local host segments, so the coupling is set by the marginal
  host-thickness ratio rather than by 1/σ(t). My first guess, a continuum
  Thevenin argument (R ~ 1/σ(t), giving Φ = 1 − κ(1−2α)/c_ℓ² and Φ(½) = 1),
  is falsified: Φ(½) ≈ 0.8–0.9.
- The per-ligament dissipation/affine at the end of the run is a function of
  relative age t/T. Old ligaments are > 1 (they collect the current dumped by
  stems): 1.3–1.5 at t/T ~ 10⁻³, falling to 0.37 for the youngest decade
  (α = 0, T = 3·10⁴). That is the "dressing" of a ligament by its
  descendants. The ψ_t flow is the cleaner object.

**Diamond-like surrogate** (`laminate.py`). Smear each chord into its cell as a
rank-one medium (λℓ/A) nnᵀ and laminate the children across it. This is bond
moving in continuum form, and exact on a hierarchical laminate. Φ_lam = 0.366,
0.480, 0.619, 0.684, 0.738, 0.815 at α = 0, 0.2, 0.4, 0.5, 0.6, 1 (T = 3·10⁴).
It is stationary for α < ½ and has the right trend, but runs about 0.2 below the
network (empty leaves act as insulating layers), and it misses the decrease for
α > ½ (free edges, focusing). It is a lower-side bracket, not a quantitative
RG.

Density exponent for transport: m = 1 for every α (σ = Φ(α) ρ/2).

## 5. Elasticity (`frame.py`, `elastic.py`, `elastic_summary.py`)

Euler–Bernoulli frame: EA = λ, EI = λ³/12, rigid joints, uniaxial x and y with
free lateral sides. The stiffness length is capped at max(L, λ). Without the
cap, the 8–16% of segments shorter than their thickness (inside a junction
block) push E above the Voigt bound (0.82 vs 0.14 at ρ = 0.4) and make spurious
jumps.

- *Pin-jointed SR is floppy.* Interior junctions have degree 3, and the collinear
  host segments give free transverse modes. E(bend = 1e−8)/E(bend = 1) ≈ 2e−7.
  All stiffness is bending-stabilised, but at ρ ≳ 0.2 most of the energy is
  axial (fraction 0.4–0.8).
- *Per-ligament energy* (T = 3000, ρ = 0.4, 3 seeds). U_i ~ i^(−x) with
  x = 0.78, 0.83, 1.04 at α = 0, 0.3, 0.6; the paper gives 3/4 + α/2 = 0.75,
  0.90, 1.05. Each factor-3 birth window holds 10–20% of the energy (all scales
  contribute). The split of U ~ k^a λ^b into a and b is not identifiable inside
  one ensemble, because k_i ~ (T/i)^½ and λ_i ~ i^(−α) are both functions of i
  (fits give a = 1.64, 1.28, 1.43 at α = 0, 0.3, 0.6 and b = 1.25, 0.88 at α = 0.3, 0.6).
- *λ0 route* (fixed T, all thicknesses scaled). 1 ≤ m ≤ 3 is a theorem: stiffness
  is cK_axial + c³K_bend. Measured m = dlnE/dlnρ goes from 2.9 at ρ = 0.01 to
  1.70, 1.40, 1.16, 1.05, 1.01 at ρ = 0.6 (α = 0, 0.2, 0.4, 0.6, 0.8; T = 3000).
  Larger α reaches the stretching regime at lower density.
- *T route* (fixed λ0, 4 seeds, τ = 10³…10⁴). m_T = dlnE/dlnρ = 2.42, 2.41,
  2.44, 2.49, 2.47 (ρ_end = 0.1) and 2.16, 2.08, 2.00, 1.90, 1.73
  (ρ_end = 0.4) at α = 0, 0.1, 0.2, 0.3, 0.4. Paper eq. (3) with p = 0 gives
  2, 2.38, 3, 4.25, 8. x_T = dlnE/dlnT = 1.23, 1.00, 0.80, 0.62, 0.44
  (ρ_end = 0.1), against 1 − α/2 = 1, 0.95, 0.9, 0.85, 0.8. In the beam frame,
  E follows ρ^(≈2–2.5) along T for every α < ½. The 1/(½−α) divergence of
  eq. (3) only re-expresses ρ ∝ T^(½−α). At finite T, ρ(T) carries an O(1)
  early-ligament constant, and U_tot(T) grows much more slowly than T^(1−α/2).
  Unresolved: the paper's U_tot protocol (SI §3) and its fitted p(α).
- *RG reading.* Two couplings flow with ln t. (i) Slenderness of the newest
  generation, s_t = ℓ_t/λ_t ∝ t^(α−½). It has the opposite eigenvalue to density
  and s_t ρ(t) → 2c_ℓ²/(1−2α) = 2.70/(1−2α). At ρ = 0.4 that is 6.7 (α = 0),
  17 (0.3), 34 (0.4): larger α at fixed ρ means slenderer newcomers. (ii) The
  marginal host ratio U^(−2α). The bending increment from a junction on a host
  ∝ (λ_h/λ_t)³ has a mean of 1/(1−6α), which diverges at α = 1/6. Above that,
  bending increments are dominated by junctions on the oldest, thickest
  ligaments. In the pure-bending limit this predicts x_T = 3/2 − 3α (α < 1/6)
  and 1 (α > 1/6). The data (axial fraction 0.3–0.7) are not in that limit and
  show no kink. Window increments ψ_el do **not** collapse on s_t alone (at
  s ≈ 100, ψ_el = 0.006 at α = 0 vs 0.13–0.2 at α = 0.4). The host ratio has to
  enter as a second scaling variable. Collapsing on (s_t, α) is open.

## 6. Failure / compaction proxy (`failure.py`)

Euler buckling P_c = π²EI/L² on the same frames under x-compression. Extremal
dynamics with avalanches at fixed ε. A buckled segment keeps its bending and
loses its axial stiffness (×10⁻³). ρ = 0.2, T = 1000, first 150 events, 3 seeds:

| system | strip fraction (random 0.16) | nn-distance ratio | mean / max avalanche | E after / before |
|---|---|---|---|---|
| square grid, 5% noise | 0.467 | 0.87 | 11.3 / 57 | 0.006 |
| SR α = 0 | 0.169 | 0.75 | 2.1 / 13 | 0.16 |
| SR α = 0.3 | 0.198 | 0.83 | 1.9 / 23 | 0.21 |
| SR α = 0.6 | 0.198 | 0.85 | 1.7 / 21 | 0.35 |
| SR α = 0.6, λ shuffled | 0.162 | 0.86 | 2.0 / 17 | 0.18 |

So SR frames do not form bands, and the age–thickness hierarchy is not why:
the shuffled control is equally delocalised. The fibre-bundle intuition (thin
= young fail first and shed little load onto thicker hosts) holds for α > 0, but
it is not needed. Geometry (segment-length disorder, T-junction topology) is
enough. A hierarchical fibre-bundle RG on the tree was not built for that
reason.

## 7. What the chygraph/RG framing adds, and where it fails

Adds:
1. The exact identification of the hierarchy with a random BST. All depth,
   area, height and profile statistics are then imported from Devroye/Mahmoud,
   and multiscaling can be ruled out for areas regardless of shape.
2. The host-age law h = tU². It turns α into an exactly marginal coupling
   between generations, and gives the thresholds α = 1/(2p) where moments of
   the host-thickness ratio diverge.
3. The exact generation flow σ(T) = Σψ_t aff_t. Its stationarity for α < ½ is a
   fixed line; it says Φ is T-independent and gives its values. α enters as
   relevant for density and as marginal for the prefactor.

Fails / only translates:
- The exact contraction is nested dissection, with root width √T. No gain over
  a sparse solver. "Narrow complexes" holds only on average.
- Chord endpoints are non-local: they land on ancestors' sides, not on the
  cell's own two-terminal boundary. The cell is therefore not a w = 2 bond
  and no closed scalar recursion exists. The laminate surrogate that forces
  one is off by about 0.2 in Φ.
- Planarity and loops across cells: the frame of a cell is shared with its
  neighbours, so a cell's message depends on outside chords (the w ~ √n
  atoms).
- For elasticity, the RG gives the right variables (s_t, U^(2α)) but no
  collapse yet. For failure it was not needed.

## Next steps (ordered)

1. Get the paper's SI §3 (U_tot protocol, p(α)) and settle why eq. (3)
   diverges where the beam frame gives m_T ≈ 2. Rerun E2 under their protocol
   (fixed ρ across α, or continuum FE energy).
2. Elastic flow: collapse ψ_el on (s_t, ⟨(λ_h/λ_t)³⟩ truncated) and test the
   α = 1/6 change in x_T in the slender limit (ρ_end ≤ 0.03, T ≥ 3·10⁴).
3. Transport fixed line: an RDE over the BST for (F_t, g_tR_t), using the host
   law and the stationary cell-shape law, to compute ψ(α) without network
   solves (book Ch. 20 machinery). The laminate is the zeroth-order closure.
4. Field focusing for α > ½: does ψ_t keep growing (∝ t^(α−½)?) or saturate?
   Needs T ≥ 10⁵ with sparse checkpoints.
5. Finite-width topology (`lng/srcells.py`): does the λ0* boundary move from
   1/(c_ℓζ(α+½)) to the paper's ζ(2α) form?
6. Slow drift of the mean cell shape ratio (0.205 → 0.180): is there a slow
   elongation mode in the iso model?

## Files

`srtree.py` generator + segment network · `tree_stats.py` · `exact_checks.py` ·
`nested.py` · `conductance.py` · `flow.py`, `flow_analyse.py` · `laminate.py` ·
`frame.py`, `elastic.py`, `elastic_summary.py` · `failure.py` ·
`rg_section.tex` (draft section).
