# lng/ — the SR model as local network growth

Workstream notes, 2026-11-12. Every number below comes from a script in this
directory, and the output files are in `results/`. All runs are single-threaded.

## 0. Code

| script | what | output |
|---|---|---|
| `srcells.py` | SR model run as iterated cell division (exact for zero thickness; thickness by half-plane offsets). Knobs: alpha, lam0, area^q selection, angle rules (set / von Mises / field with prob. pfield), nucleation density | — |
| `surrogate.py` | purely topological face-split rules (no geometry), plus DM = BB(n,1) | — |
| `graphstats.py` | degrees, C(k), k_nn, assortativity, BFS path counts, MS rewiring | — |
| `exp_topology.py` | SR vs surrogates at T=1e5 (3 seeds) | `topology_T100000.txt` |
| `exp_ck.py` | clustering mechanics, total length, k-l law | `ck_T100000.txt` |
| `exp_cornercut.py` | integral-geometry formula for the triangle closure | `cornercut.txt` |
| `exp_paths.py`, `paths_fit.py` | <d>, <mu>, <ln mu> vs T (1e3–3e5), mu(d) | `paths_scaling.csv`, `paths_mu_of_d.csv`, `paths_fit.txt` |
| `exp_geometry.py` | thickness/jamming, angle sets, spread, q, gradients | `geometry.csv` |
| `exp_polar.py`, `exp_polar2.py` | orientation fields (polarized rules) | `polar.csv`, `polar2.csv`, `polar_hubs.txt` |
| `beams.py`, `beams_pool.py` | Euler–Bernoulli frame and per-ligament energy fits | `beams_fit.txt`, `beams_pool.txt`, `beams_*.csv` |
| `sr3d.py` | 3D plates (iterated polyhedron division) | `sr3d.txt` |
| `fig_lng.py` | 4-panel figure | `fig_lng.png` |

Speed: 2e5 zero-thickness SR steps take 3 s, so no spatial hash is needed.

## 1. The SR step is a face-local rule

**Exact reformulation.** With zero thickness, a point that is uniform in the
unit square lies in cell K with probability A_K. The ligament through it, grown
both ways to the first hit, is exactly the chord of K through that point.
Ligaments never leave the cell they nucleate in, because the cell is convex and
bounded by ligaments. So the SR model is iterated cell division with
**area-weighted selection** and the **"chord through a uniform point,
uniform angle" division rule**.

This is exactly Cowan's model, Fig. 2(a) of R. Cowan, Adv. Appl. Prob. 42, 26
(2010) (verified from the paper). Cowan's dilation does not change the topology.
Known results from Cowan: daughter area fractions are uniform, and the limiting
cell areas are exponential (his Corollary 2). STIT uses perimeter-weighted
selection with uniformly random lines instead. The Nature paper cites Gilbert,
Kiang (1966) and Rosvall et al. (2005), but not Cowan. Whether the SI does is
not checked.

In the dual (ligament = node, junction = edge), the step is: pick a cell (a
face of the T-tessellation), and join a new node to the two ligaments that
carry the two boundary sides the chord hits. Exact bookkeeping:
- A cell's corner number n is fixed at birth, because a T-junction on a
  neighbour's side is a 180° vertex for that neighbour. A split gives
  n1 + n2 = n + 4, so the number-averaged mean is <n> = 4 exactly (measured
  4.0000).
- A triangle cell can only split into (3,4). So the number fraction of
  triangles is f3 = P_corner − w3, where w3 is the area-weighted triangle
  fraction. Measured: 0.30 = 0.449 − 0.147.
- Each node has 2 out-links (fewer if it ends on the boundary), so <k> → 4.
  Each later ligament ending on node i adds one link and cuts one of i's face
  sides into two. The number of face incidences s_i therefore equals k_i for
  interior nodes.

**Degree distribution (rate equation).** The hit rate on node i is the sum
over its face incidences. Per step there are 2 hits, and the total number of
incidences is 4t. If every incidence is equally likely to be hit, we get
linear preferential attachment with b = 2/<k> = 1/2: dk/dt = k/(2t). This is
the book's affine-rate equation (`distance.tex`) with a=0, b=1/2, c=1, d=0:

   p_k = 12 / [k(k+1)(k+2)],  P(k>=K) = 6/[K(K+1)],  gamma = 1 + 1/b = 3,

which is the BA(m=2) / BB(n,1) distribution.

Test (`topology_T100000.txt`, 3 seeds, T=1e5): K²P(k−2>=K) at K=8/16/32 is
3.53/4.45/5.19 for SR, 3.53/4.45/4.85 for the uniform-face surrogate, and
3.49/4.49/5.16 from the formula (fig_lng (a)). The geometric "spatial
preferential attachment" (rate ∝ length) and incidence-PA (rate ∝ k) are the
same statement. The portions of ligament i have the local cell size
ell_t ~ t^{-1/2}, so k_i ≈ l_i/ell_t.

**Geometric version (the probe's law, now with its constant).** The hit rate
per unit length per step is h(t) = 2/L(t), where L(t) = sum of lengths =
c_L √t. Measured c_L = 2.33 and <l_t>√t = 1.187. Hence

   k_i − 2 = (4 l_i / c_L)(√T − √i),

with mean ratio 1.006 and corr 0.954 (`ck_T100000.txt`).

In this version the noise sits in the birth length, x = l_i/<l_i>, with
<x²> = 1.39. In BA it sits in the Yule history (Gamma(2) gives <x²> = 1.5).
Asymptotic amplitudes are 4<x²> = 5.56 versus 6, so SR has about 7% fewer
far-tail hubs. Hub degree is predictable at birth from length; in BA it is
not.

**Triangle closure.** Integral geometry: chords through uniform points at
uniform angles have line density sigma dh dtheta/(pi A). A chord hits two
adjacent sides iff it cuts off one corner. For fixed theta, the integral of
sigma over the corner-cutting range is the area of the cut-off corner. So

   P_corner(K) = (1/(pi A)) ∫_0^pi [A_lo(theta) + A_hi(theta)] dtheta.

Values: square 2 ln2/pi = 0.4413 exactly, triangle 1, 1x4 rectangle 0.303,
regular pentagon 0.232, regular hexagon 0.136. The area-weighted average over
the cell population is 0.447–0.449, which matches the measured next-step
frequency 0.447–0.449 (`cornercut.txt`). By n=3..7: 1.000 0.461 0.254 0.152
0.10.

The probe's "0.495 triangle closure" is therefore **not 1/2**. It is
0.448 (corner cut, local) + 0.048 (the two opposite sides of the cell are
already joined elsewhere: lines meet outside the cell), totalling 0.496. The
other 0.50 of newcomers close a 4-cycle around the cell.

**Clustering.** Triangles at i form only at birth (prob. p0) or when i is hit
and the other hit side is dual-adjacent to i (prob. p per hit). Existing nodes
are never linked to each other. So e_i = p0 + p(k_i−2) and

   C(k) = 2[p0 + p(k−2)]/[k(k−1)] → 2p/k,

with p0 = 0.494 and fitted p = 0.53 (0.56 for k>7). The prediction is within
2% for k>=15 (`ck_T100000.txt`), and C(k)·k → 1.06–1.1, the probe's 1.1/k.
At k=2 the measured value is 0.536 against 0.494, because corner-cut ligaments
are short and so rarely hit again. This is the triadic-closure clustering law
of `friends.tex`, <c>_k ≈ 2mu_1/k, with 2mu_1 → 2p.

**Degree correlations.** r ≈ −0.02 in all face models.
k_nn(k) falls from 19.5 (k=2) to a plateau of 11.1–11.5 ≈ ln T for hubs
(T=1e5), the BA value (m/2) ln N with m = 2. Newcomers land on length-biased,
and hence degree-biased, sides; hubs collect neighbours of all ages
(fig_lng (c)).

**Place in the LNG parameter space.** SR dual =
- with prob. ~0.45, a corner cut: the newcomer joins both ends of a dual edge
  (the cell corner). This is BB(n,1) / Dorogovtsev–Mendes edge attachment, but
  with the edge drawn ∝ local cell size rather than uniformly.
- with prob. ~0.55, a face split across: the newcomer joins two nodes at
  distance 2 on a face cycle, closing a square.

It sits with BA(n,2), BB(n,1) and DS(n,1/3) on the gamma=3 line of the book's
master equation. Its clustering is half of DM's (C 0.38 vs 0.74), and the
missing half went into squares. It is also the 2D cousin of random Apollonian
networks (faces chosen, degree ∝ incident faces). In d dimensions RAN gives
gamma = (2d−1)/(d−1) (Kolossváry–Komjáthy–Vágó), the same formula SR gives
(§2.4).

**Topological surrogate (no geometry).** Pick a face (uniform, ∝ sides, or ∝ a
random-split "area"), then two sides.
- Every face-local variant reproduces P(k): gamma_MLE(k>=10) = 2.88–2.91,
  identical K² amplitudes. The variant with length-like side weights does not
  (gamma 3.6; it breaks incidence-PA).
- Clustering needs one parameter, the corner-cut probability for n>=4 cells:
  p4 = 0.29 gives C = 0.386 (SR 0.377), C(k) on the same curve
  (fig_lng (b)) and k_nn on the same curve.
- With uniform side pairs, P_corner = 0.74 and C = 0.59.
- **What the surrogate does not reproduce is distances and multiplicity.**
  At T=1e5, <d> = 14.1 vs 8.8 for SR, growing as 1.43 ln T vs 0.60 ln T, and
  <mu> explodes (1e5 vs 29). The surrogate creates arbitrarily deep "ladders"
  of squares. In geometry, chords favour long sides, which belong to old long
  ligaments, so geometry keeps the dual small-world. **Geometry matters only for
  metric/path properties, not for P(k), C(k) or k_nn(k).**

## 2. What geometry adds, and what can tune gamma

`geometry.csv`, T=3e4, 2 seeds; read gamma_MLE(k>=10) against SR0 2.90,
beta = cohort exponent.

1. **Thickness / density (tunable phase)** (alpha 0.5–1, lam0 0.01–0.1,
   rho = 0.12–0.45): gamma 2.85–2.92, beta 0.50–0.52, C 0.376–0.378,
   P_tri 0.49–0.50. The topology is identical; alpha only paints thickness.
2. **Jammed phase** (alpha 0–0.4): gamma 2.89–2.94, C 0.35–0.37,
   P_tri 0.48–0.49.
   - alpha=0, lam0=0.02 reaches rho=0.99 at t≈1425 and exhausts the free
     area at t≈1970. alpha=0, lam0=0.05 does so at t≈250/350.
   - Jamming truncates growth (it sets T); it does not reshape the topology.
   - The "corner-cut" fraction of free polygons rises to 0.69–0.74 (slivers),
     but dual triangle closure stays at 0.48.
3. **Angle sets.**
   - {0°, 90°}: all cells are rectangles, so P_tri = 0 and C = 0 exactly, with
     gamma unchanged at 2.87 (Krapivsky–Ben-Naim-type rectangular
     fragmentation; Cowan Fig. 7(b)). Squares only, so <ln mu> doubles
     (3.55 vs 1.52).
   - {0, 60, 120}: C = 0.17.
   - Cell selection ∝ area^q: q=0 gives C=0.45, q=2 gives C=0.35, gamma
     unchanged.
   - **The angle set and the selection weight tune clustering (0 to 0.45),
     not gamma.**
4. **Global alignment** (von Mises kappa=1/4/16 around one axis): gamma 2.86–2.91.
   Distances grow (<d> 7.7 → 15.3 at kappa=16: near-parallel ligaments).
5. **Density gradients** (linear in x, radial): nothing changes.
6. **Polarized orientation fields do tune gamma.** This is the paper's Fig. 6
   "polarized" rule (`polar.csv`, `polar2.csv`).
   - theta = circumferential field: gamma(k>=10) = 2.48–2.51, beta 0.62–0.67,
     stable over T = 1e4, 3e4, 1e5, 2e5. c2 grows to 87 at 2e5.
   - Wave field theta = pi(x+y)/w with w = 0.125–1: gamma 2.44–2.50.
   - Fraction pfield following the field: 0.2 → 2.87, 0.5 → 2.75–2.80,
     0.8 → 2.57–2.60, 1.0 → 2.49.
   - Mechanism: newcomers are locally parallel and cannot hit each other. They
     hit the ligaments that are misaligned with the local field. Misalignment
     grows with length for straight lines in a curved field, so the hit rate per
     unit length increases with l: k ~ l^{1.5–1.6} instead of l^{1.03}.
   - Junction angles become shallow (<|sin|> 0.17–0.4 vs 0.77).
   - Lines passing near the centre of the circumferential field have
     (k−2)/l = 570 vs 100–130 for distant lines (`polar_hubs.txt`).
   - Hit-rate homogeneity per incidence is the assumption that breaks.
   - Open: whether gamma → 3 as T→∞ for fields with a scale (w). No drift seen
     up to 1e5.
7. **Dimension** (`sr3d.txt`, plates in a cube, isotropic). A newcomer plate
   hits m plates (sides of its section polygon); incidences grow by 2+m per
   step and degree by 2m. So beta = m/(2+m) and gamma = 2 + 2/m.
   - Measured m = 4.47–4.57, predicting 2.44.
   - The dimensional argument k ~ (l_i/ell_T)^{d−1} gives gamma =
     (2d−1)/(d−1) = 2.5, the same as d-dim RAN.
   - Measured at T=8e4: beta 0.66, gamma_MLE(k>=40) 2.50 (2.55 at k>=20).
   - So **3D SR plates have gamma ≈ 2.5, not 3**. Check against their SI
     §1.6–1.7.

Rule-design summary (local-rule language): gamma = 1 + 1/b, with b = (links
received by old nodes per step)/(incidences created per step).
- 2D: any rule where each newcomer makes 2 junctions and old nodes are hit in
  proportion to their face incidences is pinned at b = 1/2. That covers
  thickness, jamming, angle sets, selection weights and gradients.
- Ways to move gamma:
  1. make the per-incidence hit rate node-dependent (coherent orientation
     fields: gamma down to ~2.45);
  2. go to plates in 3D (2.5);
  3. hypothetically, rules where a newcomer's own degree differs from the
     links it hands out (e.g. ligaments that cross r others: X-junctions keep
     b = 1/2; branching from a ligament keeps it too). Not tried.

## 3. Energy law (`beams_fit.txt`, `beams_pool.txt`)

**Identifiability (analytic).** At fixed alpha, lam_i = lam0 i^{-alpha}
exactly, and l_i, k_i ~ i^{-1/2} in the mean. The paper's marginals —
U ~ i^{-(3/4+alpha/2)}, k^{3/2+alpha}, ... — are all one age law rewritten.
U ~ k^a lam^b can match the age law for every alpha only if
a/2 + alpha b = 3/4 + alpha/2, i.e. a=3/2, b=1/2. With k ∝ l (k_i ≈ 2 l_i √T
at fixed T), **U ~ l^{3/2} lam^{1/2} f(theta) is the same law**. The data in
the paper cannot tell degree from length or age.

**Frame model** (Euler–Bernoulli, E=1, A=lam, I=lam³/12, rigid T-junctions,
affine uniaxial strain eps=1e-3 on boundary nodes; stubs with l<thickness are
capped). T = 1000–5000.
- Per-ligament energies span about 8 decades.
- R² with 12 orientation dummies is 0.3–0.45 per run (0.63 pooled with
  per-alpha intercepts).
- Age alone gives R² within 0.03–0.08 of the best fit.
- Within birth cohorts, the partial correlations U–k (0.22–0.30) and U–l
  (0.17–0.34) are comparable, with k slightly ahead at alpha>=0.5. Degree adds
  a little beyond length (more junctions means more loaded segments).
- Pooled over alpha = 0.25–1 (T=3000): U ∝ k^{1.64} lam^{1.61}
  (R² 0.634) ≈ l^{0.82} lam^{2.02} (R² 0.631).
- The paper's combination k^{1.5} lam^{0.5} has slope 1.63, not 1.
- The frame gives a lam exponent of 1.6–2.0, as expected between stretching
  (lam) and bending (lam³). It does not give 1/2.
- Age-law exponents: 1.15 (alpha 0, 0.25), 1.35–1.43 (0.5), 1.9–2.0 (0.75),
  2.4–2.5 (1.0), against 0.75–1.25 from 3/4 + alpha/2.
- Utot is 0.07–0.35 of the affine (Voigt) stretch bound, so the frame is
  bending-dominated for alpha >= 0.75.
- Caveat: a slender-beam frame is not their continuum FEM with thick struts,
  so the lam exponent is model-dependent. The robust statements are (i) the
  identifiability argument and (ii) k and l being statistically
  interchangeable.

## 4. Shortest-path multiplicity (`paths_fit.txt`, fig_lng (d))

Dual shortest paths are fewest-turn routes (load paths that change ligament
the fewest times).
- Degree-preserving rewired SR* follows the book's baseline: <mu>/(c2/e) =
  1.03, 0.98, 0.95, 0.94, 0.89 for T = 1e3 ... 1e5. c2 grows as ln T
  (7.9 → 15.6 for T = 1e3 → 3e5, as gamma=3 requires).
- SR: <mu> = 6.5, 9.7, 15.6, 22.3, 25.1, 39.1 for T = 1e3 ... 3e5, a boost of
  2.3 → 6.8 over c2/e, growing. The log-quadratic fit has a positive quadratic
  coefficient (0.57); <ln mu> = 0.18 ln T.
- DM / BB(n,1) stays near the baseline (1.1–1.2 × c2/e) in this range: its
  triangles add no redundancy.
- mu(d) in SR grows exponentially, about 0.45–0.6 per hop. That is the
  square-lattice signature (fewest-turn routes through a grid of ligaments).
  SR* saturates near c2.
- <d>_SR = 0.60 ln T + const (5.9 → 9.5), so the dual is small-world.
- So the SR boost comes from the ~55% square-closing splits (4-cycles), as in
  `distance.tex` "squares multiply paths", while geometry keeps the squares
  from stacking into ladders (cf. the surrogate, whose mu explodes).

## 5. Prior art (checked = fetched or seen in the paper's references)

- **Cowan 2010** (checked): this exact model (area-weighted + uniform-point
  chord), with exponential cell areas. Topology/dual degree not studied there,
  as far as read.
- **Gilbert 1967** (in the paper's refs): simultaneous crack growth from
  Poisson seeds. A different process, cited by the paper as a comparison
  architecture.
- **STIT** (Nagel & Weiss 2005; overview Redenbach & Jung, arXiv:2501.05166,
  checked): perimeter-weighted selection and uniformly random lines.
- **Road duals**: Rosvall et al. PRL 94, 028701 (2005) (in the paper);
  Kalapala et al. PRE 73, 026130 (2006) (checked abstract: dual degree
  exponent 2.2–2.4 for US/England/Denmark).
- **Random Apollonian networks**: Zhou, Yan, Wang, PRE 71, 046141 (2005),
  gamma=3, C≈0.74 (checked via arXiv listing). Exact form ∝ 1/[k(k+1)(k+2)]
  per the Wu–Xu–Wang comment. Kolossváry–Komjáthy–Vágó give d-dim
  gamma = (2d−1)/(d−1).
- Krapivsky & Ben-Naim (1994) rectangular fragmentation multiscaling: **not
  verified this session**.
- DM edge rule: the citation is to be pinned (the book calls it BB(n,1)).

**New here (as far as checked):**
1. The SR dual degree distribution is the BA(m=2) form 12/[k(k+1)(k+2)],
   explained as face-incidence PA, together with the exact length law and its
   constant.
2. Triangle closure 0.45 + 0.05, with an exact integral-geometry formula
   (square 2 ln2/pi) and the identity f3 = P_corner − w3.
3. C(k) = 2[p0+p(k−2)]/[k(k−1)].
4. gamma=3 is pinned by b = 1/2 against all of the paper's knobs except
   polarized fields (gamma ≈ 2.45–2.5, new) and 3D (gamma ≈ 2.5).
5. The energy law is unidentifiable between degree, length and age.
6. Multiplicity is boosted by squares; geometry is needed only for distances.

## Next steps (ordered)

1. Polarized fields: T → 1e6 for the wave field at w=0.125 to settle whether
   gamma → 3 asymptotically. Then derive gamma(pfield) from a two-population
   (aligned / misaligned) rate equation.
2. Read the paper's SI §1.4 (P(k) fit), §1.6–1.7 (3D exponent), §3 (energy-law
   fits; whether they fit within alpha or pooled), and check whether Cowan is
   cited.
3. Exact theory of p (per-hit triangle probability at hubs, 0.53–0.56) and of
   the non-local 0.048. Both need the joint law of cell shape next to long
   ligaments.
4. Multiplicity: decide between log-quadratic and power law for <mu>_SR with
   T to 1e6 (100 sources, <ln mu> as the stable observable). Compare with
   BB(n,1) and DS(n,1/3) at matched c2. Count squares per node.
5. Energy law with a continuum-like model (Timoshenko beams or 2D FE on thick
   struts) to see whether lam^{1/2} appears. The identifiability conclusion
   does not depend on this.
6. Confirm the Krapivsky–Ben-Naim and DM citations. Search for any prior
   statement of a T-tessellation dual degree distribution (STIT I-segment
   internal-vertex counts are the natural place).
7. 3D: larger T with interior-only statistics to pin 2.44 vs 2.5 (the
   m-formula vs the dimensional argument).

## Area selection vs degree (2026-11-12, `exp_area_degree.py`, `results/area_degree.txt`)
Question: does area-weighted cell selection route hits to high-degree ligaments?
Poisson regression of hits in [t,2t) on ln l, ln k(t), ln i, ligaments born in [t/40, t/4),
T = 1e5, three seeds, t = T/8, T/4, T/2:
- q = 1 (SR): a(ln l) = 1.00–1.10, b(ln k) = −0.09…+0.04 (mean −0.03, se 0.02–0.04).
  Degree adds nothing at fixed length and age; hits follow length.
- Catchment area per unit length vs k at fixed (age, length): r = −0.06…−0.07 (weakly negative).
- q = 0 (cells chosen uniformly, geometry kept): a = 0.01–0.02, b = 0.99–1.01.
  Hits follow degree exactly: every incidence equally likely = effective PA from local sampling.
Reading: area weighting links hits to length, not to degree. Without it, the same face-local
rule gives degree-proportional hits. Both are local; gamma = 3 either way (q=0 item above).
