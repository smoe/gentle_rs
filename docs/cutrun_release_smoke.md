# CUT&RUN Release Smoke Runbook

This runbook is the pre-release smoke/proof path for CUT&RUN support. It is
deliberately gene-agnostic in implementation: every step operates on a prepared
catalog dataset, a genome-anchored sequence id, and shared engine-owned report
records. TP73 is only a suggested proof target because it is experimentally
convenient for the near-term release.

## Scope

Use only existing contracts:

- `genomes status|prepare|extract-gene|extract-region`
- `reads acquire status|prepare|inspect` when a CUT&RUN catalog entry declares
  `reads_sra_accession`
- `cutrun list|status|prepare|project`
- `cutrun interpret|list-read-reports|show-read-report|export-coverage`
- `cutrun inspect-regulatory-support`
- `features tfbs-scan`
- `features tfbs-score-tracks-svg`
- optional DNA-window GUI inspection of the same regulatory-support report

Do not treat this runbook as a place to add de-novo motif discovery, new motif
scoring models, expression-stratified comparisons, or GUI-only biology logic.
Those are post-release work. The GUI path is limited to viewing/exporting the
engine-owned `gentle.cutrun_regulatory_support.v1` report.

If a catalog smoke dataset is SRA-backed, keep the smoke deterministic by
running `cutrun prepare` or `reads acquire prepare` explicitly;
`cutrun interpret --dataset ...` should consume the prepared manifest and not
perform surprise downloads during interpretation.

## Inputs

Choose these inputs before running the smoke:

- `STATE`: project/state file for the proof run.
- `GENOME_ID`: prepared reference genome catalog id.
- `SEQ_ID`: output id for the genome-anchored ROI sequence.
- `TARGET_GENE` or explicit `CHROM`, `START_1BASED`, `END_1BASED`: how to make
  the anchored ROI.
- `CUTRUN_DATASET_ID`: processed peaks/signal CUT&RUN catalog dataset id.
- `CUTRUN_READ_DATASET_ID` or `READ_R1`/`READ_R2`: optional raw reads for ROI
  interpretation.
- `TARGET_MOTIF`: motif/factor token used for TFBS scan and score tracks.
- `OUT`: output directory for JSON, TSV, and SVG artifacts.

TP73 example choices:

```bash
STATE=release_cutrun_tp73.gentle.json
GENOME_ID="Human GRCh38 Ensembl 116"
SEQ_ID=release_tp73_roi
TARGET_GENE=TP73
CUTRUN_CATALOG=assets/cutrun.d
CUTRUN_CACHE=data/cutrun
CUTRUN_DATASET_ID=replace_with_catalog_dataset_id
CUTRUN_READ_DATASET_ID=rostock_p73_sra_err15695857_p73_tap73alpha
TARGET_MOTIF=CTCF
OUT=exports/release_cutrun_tp73
mkdir -p "$OUT"
```

The same commands work for any other target by changing those variables.
The built-in CUT&RUN shard `assets/cutrun.d/rostock_p73_sra.json` records the
public `E-MTAB-15709` / `PRJEB100610` paired-end SRA runs. `cutrun status` is
safe for release smoke checks; `cutrun prepare` for these entries intentionally
acquires full raw reads through the shared SRA/read-acquisition path and should
only be run when disk and SRA Toolkit availability are explicit.

## Generic Smoke Path

1. Confirm the target genome is visible and prepared:

```bash
gentle_cli --state "$STATE" genomes status "$GENOME_ID" \
  --catalog assets/genomes.json \
  --cache-dir data/genomes

gentle_cli --state "$STATE" genomes prepare "$GENOME_ID" \
  --catalog assets/genomes.json \
  --cache-dir data/genomes
```

2. Create a genome-anchored ROI. For a gene-centered proof, use:

```bash
gentle_cli --state "$STATE" genomes extract-gene "$GENOME_ID" "$TARGET_GENE" \
  --occurrence 1 \
  --output-id "$SEQ_ID" \
  --annotation-scope core \
  --catalog assets/genomes.json \
  --cache-dir data/genomes
```

For a locus-only request, use `genomes extract-region` instead:

```bash
gentle_cli --state "$STATE" genomes extract-region "$GENOME_ID" "$CHROM" \
  "$START_1BASED" "$END_1BASED" \
  --output-id "$SEQ_ID" \
  --annotation-scope core \
  --catalog assets/genomes.json \
  --cache-dir data/genomes
```

3. Inspect and prepare the CUT&RUN processed dataset:

```bash
gentle_cli --state "$STATE" cutrun list \
  --catalog "$CUTRUN_CATALOG" \
  --filter "$TARGET_MOTIF"

gentle_cli --state "$STATE" cutrun status "$CUTRUN_DATASET_ID" \
  --catalog "$CUTRUN_CATALOG" \
  --cache-dir "$CUTRUN_CACHE"

gentle_cli --state "$STATE" cutrun prepare "$CUTRUN_DATASET_ID" \
  --catalog "$CUTRUN_CATALOG" \
  --cache-dir "$CUTRUN_CACHE"
```

4. Project prepared CUT&RUN peaks/signal onto the anchored ROI:

```bash
gentle_cli --state "$STATE" cutrun project "$SEQ_ID" "$CUTRUN_DATASET_ID" \
  --catalog "$CUTRUN_CATALOG" \
  --cache-dir "$CUTRUN_CACHE" \
  --clear-existing
```

5. If ROI raw reads are available, create a saved ROI read report. Dataset-backed
   reads are preferred when the catalog entry has `reads_r1`/`reads_r2` assets:

```bash
gentle_cli --state "$STATE" cutrun prepare "$CUTRUN_READ_DATASET_ID" \
  --catalog "$CUTRUN_CATALOG" \
  --cache-dir "$CUTRUN_CACHE"

gentle_cli --state "$STATE" cutrun interpret "$SEQ_ID" \
  --dataset "$CUTRUN_READ_DATASET_ID" \
  --catalog "$CUTRUN_CATALOG" \
  --cache-dir "$CUTRUN_CACHE" \
  --report-id release_cutrun_reads
```

For ad hoc paired FASTQ files:

```bash
gentle_cli --state "$STATE" cutrun interpret "$SEQ_ID" "$READ_R1" "$READ_R2" \
  --format fastq \
  --layout paired_end \
  --flank-bp 150 \
  --report-id release_cutrun_reads \
  --seed-kmer-len 9 \
  --min-seed-matches 2 \
  --max-mismatches 4 \
  --min-identity 0.9 \
  --max-fragment-span-bp 800
```

6. Export read-derived coverage summaries when a read report exists:

```bash
gentle_cli --state "$STATE" cutrun list-read-reports "$SEQ_ID"
gentle_cli --state "$STATE" cutrun show-read-report release_cutrun_reads
gentle_cli --state "$STATE" cutrun export-coverage release_cutrun_reads \
  "$OUT/cutrun.coverage.tsv" \
  --kind coverage
gentle_cli --state "$STATE" cutrun export-coverage release_cutrun_reads \
  "$OUT/cutrun.cut_sites.tsv" \
  --kind cut_sites
gentle_cli --state "$STATE" cutrun export-coverage release_cutrun_reads \
  "$OUT/cutrun.fragments.tsv" \
  --kind fragments
```

7. Inspect regulatory support from projected datasets and, when present, read
   reports:

```bash
gentle_cli --state "$STATE" cutrun inspect-regulatory-support "$SEQ_ID" \
  --dataset "$CUTRUN_DATASET_ID" \
  --read-report release_cutrun_reads \
  --neighbor-window-bp 150 \
  --species-filter human \
  --path "$OUT/cutrun.regulatory_support.json"
```

If no read report exists, omit `--read-report`; the report remains a processed
evidence smoke over peaks/signal.

8. Optional GUI inspection of the same regulatory-support report:

- Open the extracted ROI sequence in the DNA window.
- In `Engine Ops`, expand `CUT&RUN regulatory support`.
- Enter the same dataset id(s) and/or saved read-report id(s), or use
  `Use latest for sequence` when the read report was just generated.
- Keep promoter span blank for whole-sequence inspection, or enter the same
  0-based half-open interval used for CLI promoter-scoped inspection.
- Click `Inspect CUT&RUN regulatory support` and verify that the displayed
  source count, support-window count, TFBS confirmation rows, motif-absent
  windows, recurring motif context, and warnings match the CLI JSON.
- Use `Export cached CUT&RUN JSON...` when a GUI-side proof artifact is needed;
  the export reruns the shared engine operation with a path and does not use a
  GUI-only serializer.

9. Inspect motif/TFBS evidence with existing score surfaces:

```bash
gentle_cli --state "$STATE" features tfbs-scan "$SEQ_ID" \
  --motif "$TARGET_MOTIF" \
  --min-llr-quantile 0.95 \
  --max-hits 500 \
  --path "$OUT/tfbs_scan.json"

gentle_cli --state "$STATE" features tfbs-score-tracks-svg "$SEQ_ID" \
  "$OUT/tfbs_score_tracks.svg" \
  --motif "$TARGET_MOTIF" \
  --score-kind llr_background_tail_log10
```

## Signoff Checklist

- Genome status is `ready`, or prepare succeeds and records a manifest.
- The ROI sequence is genome-anchored and carries extraction provenance.
- CUT&RUN dataset status is `ready` before projection.
- Projection reports `gentle.cutrun_dataset_projection.v1` and nonzero projected
  peak and/or signal features for datasets expected to overlap the ROI.
- If reads are included, `show-read-report` returns
  `gentle.cutrun_read_report.v1`, a nonzero `total_units`, and clear mapped vs
  unmapped counts.
- Coverage, cut-site, and fragment TSVs are written for read-report proofs.
- Regulatory support JSON has schema `gentle.cutrun_regulatory_support.v1` and
  records evidence sources, support windows, confirmed/unconfirmed TFBS rows,
  motif-absent windows, and warnings.
- Optional GUI inspection displays the same regulatory-support report rather
  than producing a separate GUI-only interpretation.
- TFBS scan JSON and score-track SVG are produced from the same anchored
  sequence without requiring GUI-only state.
- Release notes record whether the proof was processed-evidence-only or also
  included raw ROI read interpretation.

## Deterministic Test Coverage

The committed deterministic tests use synthetic catalogs and toy anchored
sequences for read/projection behavior, plus a metadata-only assertion that the
built-in Rostock p73 SRA catalog shard is discoverable. The release smoke should
run:

```bash
cargo test -q cutrun
cargo test -q tfbs_score
cargo check -q
```

These tests exercise catalog discovery/status/prepare, stale prepare recovery,
anchored peaks/signal projection, ROI read interpretation, coverage export,
regulatory-support inspection, and the TFBS score surfaces used by this
runbook. A tiny TP73 raw-read subset should only be committed later if it is
derived directly from the public SRA runs with provenance.
