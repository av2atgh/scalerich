# d3/ — the SR model with plates in three dimensions

Workstream notes, 2026-11-12. Every number below comes from a script in this
directory; outputs are in `results/`. Raw runs go to `data/` (not committed;
regenerate with `python3 run.py`, about 1 min in total). Everything is single-threaded.

## 0. Code

| script | what | output |
|---|---|---|
| `sr3d.c` | engine: iterated division of convex polyhedra (cell ∝ volume, plane through a uniform point, isotropic or {x,y,z} normals, optional thickness λ_t = λ0 t^-α). 2e6 plates in 7 s and 1.7 GB | `data/*.steps/.edges/.cells` |
| `run.py` | compiles the engine and runs the configurations (iso_s1–s3 at T=2e6, iso_s4 2e5, axes 1e6/2e5, five thickness runs at 3e5) | — |
| `common.py` | loaders, per-plate arrays, MLE, cohort exponent | — |
| `exp_cells.py` | goals 1–2: splitting, BST, Euler bookkeeping, m̄, areas | `cells.txt` |
| `exp_degree.py` | goal 3: rate equations, growth exponents, snapshot test, host age, tail | `degree_iso_s*.txt` |
| `exp_mechanism.py` | goal 3: reinforcement vs fitness×age, against FIT and PA surrogates | `mechanism_iso_s*.txt` |
| `exp_cluster.py` | goal 4: closures, C(k), k_nn, 3D random Apollonian network | `cluster_iso_s4.txt`, `cluster_ran3d.txt`, `knn_scaling.txt` |
| `exp_knobs.py` | goal 5: thickness/jamming, three-axis boxes | `knobs.txt` |
| `beams3d.py` | goal 6: rods with capture radius = thickness (brute force, numba) | `beams3d.txt` |
| `section3d.tex` | LaTeX section that could replace the "Three dimensions" paragraph | — |

The old `lng/sr3d.py` (Python, T ≤ 8e4) gives the same m̄ and C. Its
"γ_MLE(k≥40) on interior plates" is a biased estimator; see §3.

## 1. Volume splitting and the binary search tree

**Proof** (same as 2D, any dimension). For a fixed normal n, the offset s of
the plane through a uniform point x of a convex body K has density A(s)/V_K,
where A(s) is the section area. The volume on one side is ∫^s A = V_K F(s),
with F the distribution function of s. So the volume fraction is F(s) ~ U(0,1),
whatever the shape of K and the law of n. With volume-weighted selection this
gives a BST on uniform keys again: volumes are the spacings of T uniform points,
P(TV > x) → e^{-x}, Σ V² = 2/(T+2), and the split depth grows as 2 ln t with variance 2 ln t.

**Measured** (`cells.txt`, T = 2e6, 3 seeds):
- daughter fraction: mean 0.4997–0.5001, var 0.0832–0.0834 (1/12 = 0.0833), KS D = 0.0006–0.0010 (n = 2e6);
- depth slope 1.993–2.010, variance slope 1.92–2.06;
- P(NV > x) at x = 1, 3, 5: 0.3677–0.3679, 0.0497–0.0499, 0.0067–0.0068 (e^{-x}: 0.3679, 0.0498, 0.0067);
- Σ V² (T+2)/2 = 0.9992–1.0010.

## 2. Bookkeeping

**Euler identities.** A plane crossing m faces of a simple polyhedron (F, V, E) gives
daughters with F1+F2 = F+m+2, V1+V2 = V+2m, E1+E2 = E+3m. Euler holds:
(2m−3m+m+2) = 2 is added to the sum of V−E+F. Starting from the cube,
- ⟨F⟩ = 2 + m̄, ⟨V⟩ = 2m̄, ⟨E⟩ = 3m̄ (number averages over cells, up to O(1/T)),

where m̄ is the step average of m (with volume-weighted selection). Measured:
⟨F⟩ = 6.5851, ⟨V⟩ = 9.1702, ⟨E⟩ = 13.7554, equal to the identities to all
printed digits. This is the 3D counterpart of ⟨n⟩ = 4 (= 2+2) in 2D.

**m̄ is a dynamical constant, not a topological one.** In 2D a chord always
crosses m = 2 sides. In 3D m ranges over 3..F, and the Euler relations above
hold for any m̄. Measured: m̄ (all faces) = 4.5851–4.5855. For interior plates
(no cube face) it is 4.4562 (t ∈ [1e2,1e3)), 4.5175, 4.5633, 4.5772 ± 0.0011, and 4.5774 ± 0.0010 (t ≥ 1e6),
so it converges to 4.577. m grows by about 1/3 per face of the parent: ⟨m | F⟩ = 3.49, 3.83, 4.18,
4.50, 4.81 for F = 4..8. The volume-weighted parent has ⟨F⟩_V = 7.32, against a number average of 6.585.
The interior m distribution for m = 3..9 is 0.116, 0.388, 0.336, 0.129, 0.027, 0.004, 0.0003.
The cell-face number fractions for F = 4..10 are 0.044, 0.184, 0.282, 0.250, 0.148, 0.064, 0.021.
For comparison, the area-weighted cell of an isotropic Poisson line
tessellation (the section of a Poisson plane tessellation) has π²/2 = 4.93
vertices. No closed form for m̄ is known to us; it would need the stationary
shape law of the volume-weighted cell.

**Areas.**
- ⟨a_t⟩ t^{2/3} → c_a = 1.229–1.235; the fitted slope of log a_t is −0.6666 to −0.6671.
- Total plate area A(t)/t^{1/3} = 3.33 (1e3), 3.54, 3.63, 3.67, and 3.677 (2e6) → 3c_a = 3.69.
- Total perimeter is 7.45 t^{2/3}; ⟨p_t t^{1/3}⟩ → 4.96.
- Plate-area distribution P(a) ~ a^{-5/2}: the MLE is 2.490–2.510 above 10, 30 and 100 ⟨a_T⟩. These are the
  paper's "power-law area distributions", and the exponent follows from a_t ~ t^{-2/3} with t uniform.
- First plate in the unit cube (MC, 4e5 draws): ⟨a_1⟩ = 0.9413 ± 0.0005, ⟨p_1⟩ = 3.88,
  ⟨m_1⟩ = 4.41. This is the 3D analogue of ℓ0 = 0.946.
- Section shape: ⟨a/V^{2/3}⟩ = 0.813, ⟨p/V^{1/3}⟩ = 4.19.

## 3. Degree: γ = 5/2, by a fitness × age law (not incidence PA)

**Two candidate rate equations.** Out-degree at birth is the number of plate
faces cut, o_i = m_i − (cube faces). Each later hit adds one link and splits
one face of i, so incidences s_i = 2 + received.
- (I) Incidence PA: hits uniform per face. Plate incidences grow by 2+m̄ per
  step and links by m̄, so b = m̄/(2+m̄) = 0.694 and γ = 2+2/m̄ = 2.44.
- (A) Area rule: the hit rate per unit plate area (per side) is uniform at a
  given time, h(t) = m_t/S(t), with S = 2A(t) + 6. A plate's area is fixed at birth.
  Then received_i = 2a_i Σ_{t>i} h(t), and with a_i = c_a x_i i^{-2/3} and A = 3c_a t^{1/3}:

      k_i = o_i + (m̄/2) x_i [(T/i)^{2/3} − 1],   β = 2/3,   γ = 1 + 3/2 = 5/2,
      K^{3/2} P(k ≥ K) → (m̄/2)^{3/2} ⟨x^{3/2}⟩.

  This is the 3D form of the 2D law k−2 = 2x[(T/i)^{1/2}−1]. In general d,
  k − o = (m̄/(d−1)) x [(T/i)^{(d−1)/d} − 1] and γ = (2d−1)/(d−1).

**Why (I) fails, and why the two are not contradictory.** Bookkeeping alone
does not decide the matter: Σ_i s_i = (2+m̄)t holds under (A) for any m̄
(Σ_i 2a_i·(3m̄/(4c_a))(t^{2/3} − i^{2/3}) = m̄t). Geometry does decide it. Once a
plate is old, its faces are cells' faces of size t^{-1/3}, so s_i ∝ a_i t^{2/3},
and its growth exponent is 2/3, not 0.694. The rate per face is not uniform:
a newborn plate is one large face per side and is hit more often per face.
Snapshot at t0 = 1e6 (`degree_iso_s*.txt`):
- the hit rate per incidence is 0.94–0.95 × mean for old plates (h/t0 < 1e-3) and 1.02 for young plates (h/t0 > 0.3);
- the hit rate per unit area is 1.04 for old plates and 0.93–0.95 for young plates;
- faces per unit area are 1.10 (old) against 0.92 (young).

Neither rate is exactly flat, but for old plates both are constant in age, so the
exponent is set by face density ∝ t^{2/3}. The +4% per-area excess of old plates
(more fragmented faces, so more face perimeter per area) multiplies the
amplitude by about 1.035^{3/2} = 1.05.

In 2D, (I) and (A) coincide because m/(2+m) = 1/2 = (d−1)/d with m = 2. In 3D
they would coincide only for m̄ = 2(d−1) = 4, which is exactly the axis-aligned box
model (§5). The incidence-PA argument in the manuscript is therefore a 2D
coincidence, as the referee says.

**Tests** (T = 2e6, 3 seeds):
1. *Fixed-plate growth.* d⟨ln s_i⟩/d ln t over the last 4 log-intervals
   (t ≳ 4e5): 0.663–0.668 for i ∈ [10,100), [100,1e3) and [1e3,1e4) (2/3 = 0.667).
   The PA surrogate gives 0.688–0.690 (theory 0.694), and the FIT surrogate gives 0.658–0.661.
2. *Area-rule prediction per plate.* Σ recv/Σ pred = 1.015 and corr = 0.991–0.992
   (interior, i ∈ [100, T/10)). By cohort the ratio is 1.04 (i ~ 1e2–1e3) and 0.94–0.98 (young).
3. *Mechanism discrimination* (`exp_mechanism.py`; the FIT and PA surrogates use
   the SR out-degrees and birth areas):
   - (a) within-cohort elasticity of received hits on birth area: SR 0.942–0.948,
     FIT 0.988–0.997, PA 0.00. With area and perimeter jointly (i ∈ [100,T/100)), recv ~ a^{0.92} p^{0.06}:
     the rate is set by the area, not by the perimeter.
   - (b) dispersion ⟨(recv−pred)²⟩/⟨pred⟩ around the area×age prediction,
     for pred bins [10,30), [30,100), [100,300), [300,1000): SR 1.47–1.53, 1.86–1.93,
     2.39–2.64, 3.94–4.03; FIT 0.96–1.01; PA 100–380. Around the cohort mean alone,
     all three give 35–190. SR is a near-Poisson fitness law with a weak extra
     per-plate factor. Its coefficient of variation is √((D−1)/μ) ≈ 0.16 → 0.07
     as μ grows: local cell-size environment and face fragmentation.
   - (c) reinforcement, i ∈ [1e3,1e4): log(hits in (2i,T]) on log(hits in (i,2i])
     at fixed area gives slope SR 0.025–0.030, FIT 0.003–0.004, PA 0.55–0.56.
     The partial correlations are 0.10–0.12, 0.01–0.02 and 0.56.
   **Verdict: 3D plates follow a fitness × age law (hit rate ∝ birth area), with
   no reinforcement beyond a weak persistent environment factor.**
4. *Degree shape, with controls instead of a bare tail slope.*
   - The CCDF ratio SR/FIT is 0.97–1.07 for K = 10..827 (max |ln| 0.07–0.10 up to K = 2000).
     SR/PA drifts from 1.06 to 0.77–0.80.
   - MLE on all plates, kmin = 80..640: SR 2.47–2.59 (mostly 2.48–2.51), FIT 2.48–2.53
     (true asymptote 2.5), PA 2.40–2.46 (asymptote 2.44). At kmin = 20 all three
     read about 0.05 high (2.52–2.58), which is the finite-k curvature the referee warns about.
   - K^{3/2}P(k≥K) for K = 80..1280: SR 4.25–4.67, flat. The prediction (m̄/2)^{3/2}⟨x^{3/2}⟩ = 4.12 (⟨x^{3/2}⟩ = 1.188),
     and with the 1.05 fragmentation factor it is 4.33.
   - *Do not restrict to interior plates for the tail.* Old large plates touch the cube
     (interior fraction 0.18 for i < 100 and 0.53 for i < 1000), and the restriction
     biases the MLE upward to 2.7–3.6 at kmin ≥ 160. Boundary plates receive hits by
     area like the others.
   - Cohort β: 0.63–0.68 (noisy; it is pulled down by o_i and the −1 term at moderate T/i).
5. *Host age.* Under (A), P(h) = a_h/A(t), so h/t → U³: P(h/t ≤ y) = y^{1/3} and
   ⟨(h/t)^α⟩ = 1/(1+3α). Measured P at y = 1e-3, 1e-2, 0.1, 0.5: 0.094, 0.217, 0.474, 0.803
   (U³: 0.100, 0.215, 0.464, 0.794). The moments for α = 1/6, 1/3, 1/2, 1 are 0.665, 0.495,
   0.393, 0.243 (0.667, 0.500, 0.400, 0.250). Incidence PA would give y^{0.306}, which is close; this test does
   not discriminate.

**Thickness consequences** (from h/t = U³). The thickness ratio is λ_t/λ_h = U^{3α}, and
⟨(λ_h/λ_t)^p⟩ = 1/(1−3αp) diverges at α = 1/(3p). For the host's bending
stiffness per unit width (plates, ∝ λ³, p = 3) the divergence is at α = 1/9.

**Fill fraction and jamming** (`knobs.txt`). A plate removes λ_t a_t, and with free
volume F spread over t cells a_t = c(t)(F/t)^{2/3}, so
dF/dt = −λ0 c(t) F^{2/3} t^{−α−2/3}, and

    F_∞^{1/3} = 1 − (λ0/3) Σ_t c(t) t^{−α−2/3}.

- Jammed for α ≤ 1/3 (the 3D boundary; 1/2 in 2D).
- For α > 1/3, jammed if λ0 ≳ 3/[c ζ(α+2/3)].
- Otherwise ρ ≈ λ0 c ζ(α+2/3) for small λ0, the analogue of the paper's ρ ≈ ℓ0 λ0 ζ(α+β) with β = 2/3.

The simulated free volume matches the recursion (with the measured zero-thickness c(t)) at t = 1e2..3e5:
- (α, λ0) = (0.8, 0.05): 0.8656 vs 0.8647;
- (0.5, 0.02): 0.8716 vs 0.8725;
- (0.4, 0.02): 0.7995 vs 0.8011;
- (0.25, 0.02): 0.5377 vs 0.5417, still decreasing;
- (0, 0.01): 0.0610 vs 0.0842 at 1e5, and 0.0008 vs 0.0063 at 3e5, nearly jammed.

The linear estimate λ0 c ζ overestimates ρ for large λ0 (0.17 vs 0.13 at (0.8, 0.05)), because it neglects
the F^{2/3} depletion. The paper's main text states only that 3D gives "an analytically defined phase
diagram"; whether their boundary is α = 1/3 needs SI §1.6.

## 4. Clustering, k_nn, and the random Apollonian network

**Closures** (`cluster_iso_s4.txt`, T = 2e5). A newcomer's section polygon is a
cycle of faces, and consecutive faces meet along a cell edge, so their plates are in
contact. Measured: 4.388 consecutive plate pairs per newcomer, linked with
fraction 1.00000. There are 3.94 non-consecutive pairs, 64% of them already linked (contacts
outside the cell), adding 2.52. So 6.91 triangles are closed per newcomer
(7.17 for interior newcomers).

**C(k).** When i is hit, its two polygon neighbours are linked to i (+2), and on
average 1.26 more hit plates are linked to i (non-local). The fit over recv ≥ 20 is
e_i = 6.3 + 3.26 recv_i, so

    C(k) = 2[e0 + q(k − o)]/[k(k−1)] → 2q/k ≈ 6.5/k.

The measured k·C(k) rises to 6.2–6.3 for hubs. With e0 = 7.17 and q = 3.26 the formula matches
within 1% for k ≥ 16 (k ~ 16: 0.3607 vs 0.3601; k ~ 64: 0.1007 vs 0.1009; k ~ 256: 0.0257 vs 0.0258). At k ~ 8 it gives 0.606
against 0.644. ⟨C⟩ = 0.717, and the transitivity is 0.031.

**3D random Apollonian network** (pick a tetrahedron, join to its 4 vertices; N = 2e5):
γ_MLE = 2.49–2.50 (k ≥ 40..160), triangles exactly 3 per new link, so C(k) → 6/k;
⟨C⟩ = 0.824; k_nn of hubs 51.5. The SR plate dual is RAN-like in 3D. It has the same γ
(both are fitness×age / face-counting with b = 2/3; for RAN b = 2/3 exactly, since a
node gains 2 tetrahedra per hit while the total grows by 3). C(k) has nearly the same
prefactor (6.3–6.5 vs 6). The mean C is lower (0.72 vs 0.82) because a young plate's
neighbourhood is an m-cycle (C = 2/(m−1) plus non-local links), not a K4.
For contrast, the 2D SR has ⟨C⟩ = 0.38 against 0.77 for 2D RAN.

**k_nn and assortativity.** r = −0.036 (2e5) and −0.019 (2e6). k_nn(k) falls with k
and scales as T^{1/3}: k_nn(k~4)/T^{1/3} = 3.83 (2e5), 3.70 and 3.64 (2e6, two
seeds). At 2e6, k_nn at k ~ 4..256 is 467, 383, 280, 204, 147, 105 and 80.
Mean field: the out-neighbours of a node born at t are hosts with P(h) ∝ h^{-2/3},
whose degree is ∝ (T/h)^{2/3}. The sum converges at small h, so the k_nn of a young node is
≈ (m̄/6) ζ(4/3) T^{2/3} t^{-1/3}, i.e. ∝ T^{1/3} at t ~ T. For a node of degree k
(born at i ~ T k^{-3/2}) this gives k_nn ≈ c T^{1/3} k^{-1/2} + (in-neighbour term). The
measured local slope is −0.3 to −0.5, consistent with an additive constant. In 2D,
k_nn saturates at ln T instead; with γ < 3, k_nn diverges with T.

## 5. Knobs

- **Thickness, tunable and jammed** (`knobs.txt`, T = 3e5). MLE (all plates, k ≥ 40, 80):
  2.49–2.55. ⟨C⟩ = 0.713–0.717, k·C(k) for hubs 6.1–6.2, m̄ = 4.52–4.58. Thickness only
  paints, as in 2D. At α = 0 with near-complete jamming, m̄ falls to 4.00, ⟨C⟩ to 0.678 and kC to 5.3.
  The growth is truncated, not reshaped (MLE 2.46–2.62 over the available range).
- **Three orthogonal orientations → boxes.** Every cell is a box (F = 6) and every
  plate a rectangle (m = 4), exactly. Incidence PA would give b = 4/6 = 2/3, the
  same as the area rule: this is the 3D case where the two pictures coincide.
  MLE at k ≥ 40..320 is 2.47–2.60 (mostly 2.50–2.53), and cohort β is 0.64–0.67. K^{3/2}P(k≥K) = 3.3–3.5,
  lower than isotropic because m̄ = 4.
  - ⟨C⟩ = 0.548, kC(k) → 3.8, with 1.94 triangles per received link. These are only the two polygon
    neighbours: boxes have no non-local closures.
  - In 2D the axis-aligned version has C = 0, because a chord always hits opposite sides. In 3D the
    four sides of the rectangle are consecutively adjacent, so boxes cluster.
  - **Orientation sets change C, not γ.**
- Nothing found moves γ off 5/2 for plates, except the rod geometry (§6). Coherent orientation
  fields (the 2D lever) were not tried in 3D.

## 6. Beams (rods) in 3D

Zero-width lines in 3D meet with probability zero, so the beam model needs a
capture rule. The natural one: rod t has radius R_t = λ_t/2 and grows from a
free nucleation point until its tip comes within R_j + R_t of rod j's axis, or
reaches the cube.

Mean field: rod i is hit at a rate ∝ its capture cross-section ℓ_i(R_i + R_t) ≈ ℓ_i R_i,
and the mean free path is ℓ_t ~ 1/W with W = Σ ℓ_j R_j. This gives
W ~ t^{(1−α)/2}, ℓ_t ~ t^{−(1−α)/2}, and k_i ~ (t/i)^{(1+α)/2}, so

    β = (1+α)/2,  γ = (3+α)/(1+α)   (3 at α = 0, 2.33 at α = 1/2, → 2 at α → 1).

Rod volume Σ ℓ_t R_t² ~ Σ t^{−(1+3α)/2} converges iff α > 1/3. So for beams **thickness
tunes γ**, because the target size of a rod is its thickness. For plates and 2D ligaments
the target is the area or length, which thickness does not change.

Simulation (`beams3d.py`, brute force, T = 1e5, one seed each):
- α = 0, λ0 = 0.004 (volume fraction 0.088): ℓ_t ~ t^{−0.53} (theory 0.5); 93% of ends land on rods;
  the fixed-rod growth exponent is 0.556 (theory 0.5); MLE 3.0 at k ≥ 10, with kmax = 74 too small for a tail.
- α = 0.5 (λ0 = 0.02) and α = 0.8 (λ0 = 0.05): still boundary dominated. 32% and 62% of ends hit the
  cube; ℓ_t exponents are 0.21 and 0.07 against 0.25 and 0.10; fixed-rod growth is 0.61 and 0.62 against 0.75 and 0.90.
  The MLE is 2.9–3.1 and 2.2–2.35, and cohort β is 0.53 and 0.77.
- In all runs the hits scale with rod length at fixed age (elasticity 0.92–1.04), so the fitness picture holds.
  The (3+α)/(1+α) law is **not yet tested**: it needs T ≳ 1e6 with a spatial grid, or a larger λ0 to leave the
  boundary-dominated regime sooner.

**The paper's beam rule is unknown to us.** Fig. 6a(vi) shows beams, but the
capture or termination rule, and whether thickness enters it, is in SI §1.7.

## 7. What needs the paper's SI (§1.6–1.7)

- The 3D phase diagram: is the boundary α = 1/3, and is the prefactor
  a0 = ⟨a_1⟩ = 0.941 or c_a = 1.23 with ζ(α+2/3)?
- Which "area distribution" is power law: plate areas (we get a^{-5/2}) or cell volumes (exponential here)?
- Their degree exponent for plates and beams, and the beam termination rule.
- Whether nucleation is uniform in the domain with rejection in material (assumed here).

## Next steps (ordered)

1. Rods: add a uniform grid, run T = 1e6 at α = 0.5 and 0.8 with λ0 large enough to
   leave the boundary regime, and test β = (1+α)/2 and γ = (3+α)/(1+α).
2. Manuscript: replace the incidence-PA wording for 2D by the fitness×age
   (metric) law, with incidence PA as the 2D coincidence m/(2+m) = (d−1)/d. Use
   `section3d.tex` for Sec. VII.
3. Closed form or bound for m̄ = 4.577 (stationary shape law of the volume-weighted cell);
   the 3D analogue of the corner-cut integral for the per-face hit probability.
4. Coherent orientation fields in 3D (does the 2D lever, γ ≈ 2.45, have a 3D analogue below 2.5?).
5. Figure for the 3D section (CCDF SR/FIT/PA, fixed-plate growth, C(k) vs RAN).
