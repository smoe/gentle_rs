# TP73 Genome-Anchored Evidence Viewer Runbook

This runbook is the public release proof for the TP73 evidence viewer. It uses
only committed public/local assets and tiny synthetic proof fixtures. It does
not download a full genome, full UCSC `rmsk`, raw CEL files, or SRA reads.

## Inputs

- TP73 locus: `test_files/tp73.ncbi.gb`
  - GRCh38.p14 Primary Assembly
  - chromosome 1, 3652516..3736201
- Proof fixtures: `test_files/fixtures/evidence_viewer/`
  - UCSC `rmsk`-style repeat rows and generated interval index
  - Clariom D-style microarray track manifest and per-contrast TSVs
  - CUT&RUN-style BED6 intervals

All proof fixture provenance and regeneration notes are in
`test_files/fixtures/evidence_viewer/README.md`.

## Regenerate The Repeat Sidecars

These commands are needed only after editing
`tp73_evidence_viewer.rmsk.hg38.txt`:

```bash
cargo run --quiet --bin gentle_cli -- resources sync-ucsc-rmsk \
  test_files/fixtures/evidence_viewer/tp73_evidence_viewer.rmsk.hg38.txt \
  test_files/fixtures/evidence_viewer/tp73_evidence_viewer.rmsk.hg38.json \
  --assembly hg38

cargo run --quiet --bin gentle_cli -- resources prepare-ucsc-rmsk-index \
  test_files/fixtures/evidence_viewer/tp73_evidence_viewer.rmsk.hg38.json \
  test_files/fixtures/evidence_viewer/tp73_evidence_viewer.rmsk.hg38.interval-index.json
```

## Headless Proof

Run from the repository root:

```bash
cargo run --quiet --bin gentle_cli -- --state /tmp/tp73_evidence_viewer.state.json \
  workflow @docs/examples/workflows/tp73_genome_evidence_viewer_release_proof.json
```

Expected artifacts are written under the git-ignored `artifacts/` directory:

- `artifacts/tp73_evidence_viewer/tp73_evidence_viewer.linear.svg`
- `artifacts/tp73_evidence_viewer/tp73_evidence_viewer.splicing.expert.svg`
- `artifacts/tp73_evidence_viewer/tp73_evidence_viewer.tfbs_score_tracks.svg`
- `artifacts/tp73_evidence_viewer/tp73_evidence_viewer.repeat_materialization.json`
- `artifacts/tp73_evidence_viewer/tp73_evidence_viewer.tfbs_score_tracks.json`

The workflow should preserve the TP73 genome anchor as `GRCh38.p14`, chromosome
`1`, 3652516..3736201, and should materialize:

- three repeat features from the tiny rmsk index,
- four projected array interval features from the two proof contrasts,
- two overlapping CUT&RUN-style BED track features,
- TFBS features plus a TFBS score-track SVG over the first 1200 bp.

## Manual GUI Smoke

Open the generated state:

```bash
cargo run --quiet --bin gentle -- --project /tmp/tp73_evidence_viewer.state.json
```

Smoke checklist:

- Confirm the active TP73 sequence shows a genome anchor/build status for
  `GRCh38.p14`.
- Toggle repeat features off and on.
- Toggle array tracks off and on.
- Keep TFBS display visible and inspect the first 1200 bp viewport.
- Select or hover one exon and confirm exon length modulo 3 / frame hints.
- Select one transcript or intron-adjacent transcript feature and confirm
  identity and genomic context are understandable.
- Select one repeat and confirm name, class/family, score/divergence, and
  genomic interval are shown.
- Select one array row and confirm contrast, `logFC`, `adj.P.Val`, probeset,
  transcript cluster, exon id, and assembly/projection status are shown.
- Select one CUT&RUN-style BED row and confirm track source/name/file, genomic
  interval, score, strand, and note are shown.
- Copy details for one repeat, one array row, and one BED track into a scratch
  note to confirm the details are usable outside the GUI.

## Deterministic Tests

Run the targeted proof test:

```bash
cargo test -q workflow_examples_tp73_evidence_viewer_release_proof_writes_artifacts_and_features
```

Useful neighboring checks:

```bash
cargo test -q genome_track_feature_details_include_provenance_and_coordinates
cargo test -q set_display_visibility_controls_array_features
cargo test -q set_display_visibility_controls_repeat_features
cargo test -q project_microarray_track_forward_anchor_materializes_array_features
cargo test -q materialize_repeat_features_creates_repeat_regions_on_reverse_anchor
```

Before handoff or release:

```bash
cargo test workflow_examples -- --test-threads=1
cargo check -q
git diff --check
```
