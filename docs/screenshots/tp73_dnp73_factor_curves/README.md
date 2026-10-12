# TP73 / DeltaNp73 local-factor tutorial evidence

This directory retains publication-safe manual/hybrid evidence for tutorial
08.17. The run used only the committed public `test_files/tp73.ncbi.gb` RefSeq
excerpt and an exact-revision Linux build at
`8652ecf29cf970e2d6c5df7e544c33f2bec116d4`.

- `01-tss-inventory.raw.png` shows all five annotated TP73 transcript starts.
  The first, six-transcript P1/TAp73 row is visibly refused as
  `missing_flanks`; GENtle neither pads nor hides it.
- `02-annotated-window.raw.png` shows the approved internal DeltaNp73 window
  before scoring: 701 bp, six transcript lanes, aligned local/genomic/TSS
  rulers, and no local motif curves.
- `03-local-factor-curves.raw.png` shows the same window after local scoring for
  TP73 (`MA0861.2`), E2F1 (`MA0024.3`) and PATZ1 (`MA1961.2`) with
  `llr_background_tail_log10`.

The score lanes are separate, strand-aware, non-cross-calibrated displays. They
are not occupancy, binding, reporter-activity or biological-regulation claims.
The GUI correctly refused an unbound Agent Assistant scoring request instead of
guessing a project sequence; this capture used the equivalent explicit GUI
controls after that safety check. The headless `promoters tss-view-svg` route
was separately replayed on the same sequence and settings.

The adjacent semantic snapshots retain native-window geometry.
[`evidence.json`](evidence.json) binds revision, binary, public input, inventory,
headless SVG and screenshots by SHA-256. This is Linux/X11 manual/hybrid
evidence, not Windows/macOS package acceptance or human scientific approval.
