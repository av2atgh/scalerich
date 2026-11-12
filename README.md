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
