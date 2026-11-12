# scalerich

Network-science readings of the "scale-rich metamaterials" model
(Both, Chen, ..., Portela & Barabási, Nature 659, 345 (2026),
doi:10.1038/s41586-026-11191-2; PDF in ~/Downloads, text in `paper/`).

Model: random nucleation point + random angle, ligament grows both ways to first
hit (T-junctions only), thickness λ_t = λ0 t^-α. Dual graph: ligaments = nodes,
junctions = edges. Claims P(k)~k^-3 via "spatial preferential attachment",
U_i ~ k^{3/2} λ^{1/2} f(θ).

## Probe results (2026-11-12, `common/sr_probe.py`, T = 20k, 3 seeds)
- γ ≈ 2.9; mean k of birth cohorts follows (T/i)^{1/2} (BA).
- k_i − 2 ≈ c ℓ_i (√T − √i), corr 0.94 → degree slaved to length × age.
- New ligament closes a dual triangle with prob 0.495 (adjacent cell sides),
  else a 4-cycle; mean clustering 0.38; C(k) ≈ 1.1/k; k_nn(k) decreasing.
- c2 = <k(k−1)>/<k> ≈ 12 at T = 20k (log-divergent for γ = 3).
- Topology independent of (α, λ0) in the tunable phase: α only paints thickness.

## Workstreams
- `lng/` — local network growth: the SR step as a face-local rule
  (pick cell ∝ area, join two boundary members); exact rate equations,
  mapping onto the LNG local-rule / triadic-closure families, what α and
  nucleation rules can and cannot tune, energy law refit on (ℓ, λ, age).
- `rg/` — chygraph / renormalisation: the cell-splitting hierarchy as nested
  complexes; real-space RG in the spirit of the diamond hierarchical lattice
  (Migdal–Kadanoff engine in chygraph/statmech/probe/rg_chygraph.py), for
  elastic/rigidity, percolation and transport on the SR structure.

Status and next steps: see `lng/NOTES.md`, `rg/NOTES.md`.

## Manuscript (2026-11-12)
`manuscript.tex` (PRE, revtex4-2, pdflatex + bibtex), figures from `figures.py`.
Anchored on local growth + Cowan's recursive cell division. The chygraph/RG
component (`rg/`) is dropped from the paper; only the Cowan-level facts from it
(uniform split → BST, host-age law h/t → U²) are used. `rg/` stays as a record.

## Pending revisions (apply at resubmission; manuscript frozen until then)
1. **"Non-local closures" is a misnomer (raised 2026-11-12).** Every SR step is
   face-local: the new node joins two sides of the chosen face. The 0.048 of
   steps that close a triangle without cutting a corner hit two non-adjacent
   sides whose ligaments already meet outside the cell. Both targets are on the
   face and already linked, so this is ordinary triadic closure (attach to both
   ends of a link), as local as a corner cut. Only the link is off the face.
   Changes:
   - Sec. IV paragraph title "Non-local closures" → "Closures through a link off
     the face"; add: "These steps are as local as the corner cuts: both targets
     lie on the chosen face and are already linked, so the new node closes a
     triangle on an existing link; only that link is not a side of the face."
   - Local-rule summary (end of Sec. IV): two moves, both confined to one face.
     ≈0.50 edge attachment (DMS), of which 0.45 on a corner link of the face and
     0.05 on a link joining two of its sides outside it; ≈0.50 attachment at
     distance ≥ 2 around the face, mostly closing a square. Drop "closes a
     triangle non-locally".
   - Conclusions, open points: "the non-local closure fraction 0.048" → "the
     fraction 0.048 of closures through a link off the face".
   - Fig. 1 caption already says "close a triangle with a ligament outside the
     cell"; fine as is. The 3D section does not use the term.
2. **"Fitness model" framing overreaches (raised 2026-11-12).** A plain fitness
   model (Bianconi–Barabási, targets drawn globally) gives C → 0; it cannot be
   the classification. Mechanisms A and B of Sec. IV are both face-local; they
   differ only in how the two sides of the chosen face are picked: A uniformly
   per side (each hit adds a side → reinforcement → BA), B in proportion to side
   length (length fixed at birth → no reinforcement). Clustering, C(k), k_nn come
   from face locality and do not discriminate (the A surrogate matches them with
   p₄ = 0.29); the tests of Sec. IV select B only for the degree law. Changes:
   - Name the mechanisms by the within-face choice: "A: uniform choice of sides
     (reinforcing)", "B: length-weighted choice of sides (no reinforcement)".
   - Mechanism B paragraph: replace "It is a fitness model [bianconi2001] in
     which the fitness is the birth length" by "Within the face-local rule, the
     birth length plays the role of a fitness [bianconi2001] for the degree
     only; targets are still confined to one face, which is what produces the
     clustering of Sec. V. A fitness model with global attachment would not."
   - Introduction: do not set local growth and fitness as alternatives. Replace
     "grows by a face-local rule, as local-growth models do, but that its degree
     law is of the fitness type" by "grows by a face-local rule in which the two
     targets within the face are chosen by length, so the degree is set by a
     birth attribute rather than by reinforcement".
   - Abstract: "fluctuation tests select the fitness mechanism" → "fluctuation
     tests select length weighting over reinforcement". Same in Sec. VIII (3D).
   - Conclusions: say explicitly that the local rule fixes the clustering and the
     length weighting fixes the degree fluctuations; both are needed.
