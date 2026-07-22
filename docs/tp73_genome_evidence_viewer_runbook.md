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
  - a directly imported CUT&RUN-style BED6 compatibility interval set
- Offline CUT&RUN V1-V3 fixture: `assets/cutrun.json` and
  `assets/cutrun/tp73_release_smoke/`
  - prepared/projected synthetic BED6 evidence
  - four deterministic synthetic paired FASTA fragments derived from the
    committed TP73 sequence
- Vendor-derived Clariom D TP73 subset:
  `test_files/fixtures/microarray_tracks/clariomd.tp73_vendor_subset.manifest.json`
  - uses real Thermo Fisher Clariom D Human na36 hg38 TP73 probeset IDs,
    transcript-cluster IDs, exon IDs, and hg38 coordinates from the derived
    gene-panel fixture
  - keeps expression/statistical values synthetic so no CEL-derived biological
    conclusion is implied

All proof fixture provenance and regeneration notes are in
`test_files/fixtures/evidence_viewer/README.md` and
`test_files/fixtures/microarray_tracks/README.md`, with the CUT&RUN-specific
record in `assets/cutrun/tp73_release_smoke/README.md`.

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
- `artifacts/tp73_evidence_viewer/tp73_evidence_viewer.cutrun_regulatory_support.json`
- `artifacts/tp73_evidence_viewer/tp73_evidence_viewer.cutrun_coverage.tsv`
- `artifacts/tp73_evidence_viewer/tp73_evidence_viewer.cutrun_cut_sites.tsv`
- `artifacts/tp73_evidence_viewer/tp73_evidence_viewer.cutrun_fragments.tsv`

The workflow should preserve the TP73 genome anchor as `GRCh38.p14`, chromosome
`1`, 3652516..3736201, and should materialize:

- three repeat features from the tiny rmsk index,
- four projected array interval features from the two proof contrasts,
- two overlapping prepared/projected CUT&RUN BED track features,
- four concordant CUT&RUN paired-read fragments in two support clusters,
- one V3 report combining the prepared dataset and saved V2 read report,
- TFBS features plus a TFBS score-track SVG over the first 1200 bp.

The zero-flank V2 step maps directly against the imported anchored TP73
sequence, so this proof intentionally does not require a prepared whole-genome
reference. All CUT&RUN values are synthetic software-test evidence and must not
be interpreted as experimental TP73 occupancy.

## Headless/Agent Smoke Queries

The inner Agent Assistant, MCP tools, or command-line shell can verify the same
evidence fields without clicking the GUI by using `features query` with
qualifiers included:

```bash
cargo run --quiet --bin gentle_cli -- --state /tmp/tp73_evidence_viewer.state.json \
  shell 'features query tp73_evidence_viewer --qual-contains gentle_generated=ucsc_rmsk --include-qualifiers --limit 10'

cargo run --quiet --bin gentle_cli -- --state /tmp/tp73_evidence_viewer.state.json \
  shell 'features query tp73_evidence_viewer --qual-contains gentle_track_source=Array --include-qualifiers --limit 10'

cargo run --quiet --bin gentle_cli -- --state /tmp/tp73_evidence_viewer.state.json \
  shell 'features query tp73_evidence_viewer --qual-contains gentle_track_source=BED --include-qualifiers --limit 10'
```

Expected machine-readable checks:

- repeat query: `matched_count=3`, with `rmsk_class`, `rmsk_family`, `score`,
  `rmsk_divergence_percent`, and genomic coordinate qualifiers.
- array query: `matched_count=4`, with `gentle_array_contrast`, `logFC`,
  `adj_P_Val`, `feature_id`, `transcript_cluster_id`, `exon_id`, and
  `gentle_array_projection_status`.
- BED query: `matched_count=2`, with `gentle_track_name`,
  `gentle_track_file`, `score`, `bed_strand`, and genomic interval qualifiers.

The proof array fixture contains three TSV rows per contrast. One row per
contrast is deliberately outside the anchored TP73 sequence, so the projected
viewer contains four array interval features: two rows from each contrast.

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
- In `Engine Ops -> CUT&RUN regulatory support`, use dataset
  `tp73_release_smoke_synthetic`, read report `tp73_release_smoke_reads`,
  catalog `assets/cutrun.json`, and cache
  `artifacts/tp73_evidence_viewer/cutrun_cache`; confirm two evidence sources
  and two support windows, then export the cached JSON.
- Copy details for one repeat, one array row, and one BED track into a scratch
  note to confirm the details are usable outside the GUI.

## Deterministic Tests

Run the targeted proof test:

```bash
cargo test -q workflow_examples_tp73_cutrun_release_proof_writes_artifacts_and_features
```

Useful neighboring checks:

```bash
cargo test -q genome_track_feature_details_include_provenance_and_coordinates
cargo test -q set_display_visibility_controls_array_features
cargo test -q set_display_visibility_controls_repeat_features
cargo test -q project_microarray_track_forward_anchor_materializes_array_features
cargo test -q project_microarray_track_uses_vendor_subset_on_tp73_genbank_anchor
cargo test -q materialize_repeat_features_creates_repeat_regions_on_reverse_anchor
cargo test -q cutrun
```

Before handoff or release:

```bash
cargo test workflow_examples -- --test-threads=1
cargo check -q
git diff --check
```
