# Reporter-vector fixture provenance

## `synthetic_mcs_backbone.gb`

- Origin: hand-authored, repository-owned synthetic DNA fixture.
- Authored: 2026-08-27.
- Deterministic recreation: concatenate the 70 bp synthetic MCS-layout prefix
  recorded in the GenBank `misc_feature`, a 29 bp `ACGT` spacer, the 81 bp
  marker `ATG + GCC x 25 + TAA`, and `ACGT x 15`; retain the feature intervals
  and circular topology shown in the file.
- GENtle use: exact helper-vector validation tests in
  `src/engine/ops/reporter_ops.rs` exercise positive generic validation and
  derived restriction-site reporting without downloading a commercial vector.
- Limitation: this record is not Promega pGL4.10[luc2], contains no functional
  luciferase CDS, and must be rejected when the pGL4.10 catalog identity is
  requested.

The separate tutorial backbone at
`data/tutorial_inputs/gentle_mammalian_luciferase_backbone_v1.gb` remains the
offline VKORC1 teaching asset. Its own provenance likewise states that it is
not an exact deposited plasmid sequence.
