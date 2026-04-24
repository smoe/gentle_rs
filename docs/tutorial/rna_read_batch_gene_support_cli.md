# Batch-Compare cDNA FASTA.GZ Samples for One Target Gene (CLI Tutorial)

> Type: `CLI walkthrough + shared-shell parity`
> Status: `manual/hybrid`
> Drift note: this page is hand-written, but it is intentionally tied to the
> shared `rna-reads ...` command family and the engine-owned sample-sheet /
> gene-support export contracts.

This tutorial shows a practical batch workflow for the current RNA-read
mapping/export path when you have many cDNA read files and want one combined
table for downstream comparison.

The target scenario is:

- about `20` local `fasta.gz` or `fa.gz` cDNA read files,
- one target gene of interest, and
- one final TSV that records:
  - overall read abundance counters per sample,
  - target-gene assigned-read abundance,
  - exon-pair co-presence for the target gene,
  - overall mean cDNA length, and
  - mean cDNA length among reads assigned to the target gene.

You will:

1. choose one dedicated state file and output directory,
2. point GENtle at an already prepared target locus,
3. interpret each gzipped FASTA file into one saved RNA-read report,
4. run phase-2 alignment for each saved report,
5. export one combined sample sheet with target-gene cohort columns, and
6. optionally inspect one sample in more detail with the summary/audit tools.

The important current boundary is:

- the engine owns report interpretation, alignment, and cohort aggregation,
- `rna-reads export-sample-sheet` now owns the per-sample target-gene summary
  columns, but
- there is not yet a single one-shot manifest command that ingests all raw
  FASTA files directly, so the outer batch driver is still a shell loop.

## Inputs and Assumptions

This tutorial assumes all of the following are already true:

- you have one local directory containing your cDNA read files,
- each file is `*.fa.gz` or `*.fasta.gz`,
- you already have a target-locus sequence in one GENtle state file, and
- you know which `mRNA` feature should seed the RNA-read mapping run.

For a concrete TP53 starting point, first work through:

- [`docs/tutorial/generated/chapters/12_tp53_multi_gene_sparse_mapping_online.md`](./generated/chapters/12_tp53_multi_gene_sparse_mapping_online.md)

If you do not yet know the seed `FEATURE_ID`, query the loaded sequence first:

```bash
cargo run --bin gentle_cli -- --state "$STATE" features query "$SEQ_ID" --kind mRNA --label "$GENE_ID"
```

What to look for:

- a stable `feature_id` for the transcript/gene you want to seed from,
- the transcript label you expect, and
- the correct sequence id before you start the batch run.

## What You Will Verify

By the end, you should have confirmed all of these:

- GENtle accepts gzipped FASTA input directly for RNA-read interpretation,
- each input file becomes one persisted RNA-read report with a stable
  `report_id`,
- phase-2 alignment can be rerun for all retained rows with `selection=all`,
- `rna-reads export-sample-sheet` can aggregate many saved reports into one TSV,
- the sample sheet now includes:
  - `mean_read_length_bp`,
  - `gene_support_accepted_target_count`,
  - `gene_support_accepted_target_fraction_total`,
  - `gene_support_accepted_target_fraction_aligned`,
  - `gene_support_fragment_count`,
  - `gene_support_complete_count`,
  - `gene_support_mean_assigned_read_length_bp`,
  - `gene_support_exon_pair_support_json`, and
  - `gene_support_direct_transition_support_json`,
- and one saved report can still be inspected in detail via
  `summarize-gene-support` and `inspect-gene-support`.

## Step 1: Pick a Dedicated State and Output Directory

Use one dedicated state file so the batch export can safely use the exact
report ids you generate in this tutorial.

```bash
STATE=/tmp/gentle_rna_batch_tp53.gentle.json
OUT=/tmp/gentle_rna_batch_tp53
READ_DIR=/path/to/your/20_cdna_fastas
SEQ_ID=grch38_tp53
FEATURE_ID=0   # replace with the target mRNA feature_id
GENE_ID=TP53

mkdir -p "$OUT"
rm -f "$STATE" "$OUT"/report_ids.txt
```

If you are not working on TP53, substitute your own values for `SEQ_ID`,
`FEATURE_ID`, and `GENE_ID`.

Quick sanity check that you really have the expected file count:

```bash
find "$READ_DIR" -maxdepth 1 \( -name '*.fa.gz' -o -name '*.fasta.gz' \) | wc -l
```

Expected outcome for this scenario:

- the command prints `20`

## Step 2: Decide Whether You Need `single_gene` or `multi_gene_sparse`

Use `single_gene` when your target locus is biologically isolated enough that
you do not need nearby paralog rescue/contrast.

Use `multi_gene_sparse` when close relatives matter and you want the alignment
step to choose among a small requested gene family instead of only one target.

Recommended first-pass choice:

- `single_gene` if you are just quantifying one clean locus
- `multi_gene_sparse` if you are working in a TP53/TP63/TP73-style family

This tutorial shows the simpler `single_gene` path first, then the sparse
variant immediately after.

## Step 3: Interpret Each cDNA FASTA.GZ File into One Saved Report

The loop below:

- scans both `*.fa.gz` and `*.fasta.gz`,
- assigns one stable `report_id` per sample,
- stores all report ids in `report_ids.txt`, and
- leaves each report persisted in the shared state file for later alignment and
  export.

```bash
find "$READ_DIR" -maxdepth 1 \( -name '*.fa.gz' -o -name '*.fasta.gz' \) -print0 \
  | sort -z \
  | while IFS= read -r -d '' path; do
      base=$(basename "$path")
      sample=${base%.fasta.gz}
      sample=${sample%.fa.gz}
      report_id="${GENE_ID,,}_${sample}"
      echo "$report_id" >> "$OUT/report_ids.txt"

      cargo run --bin gentle_cli -- --state "$STATE" rna-reads interpret \
        "$SEQ_ID" "$FEATURE_ID" "$path" \
        --report-id "$report_id" \
        --profile nanopore_cdna_v1 \
        --format fasta \
        --scope all_overlapping_any_strand \
        --origin-mode single_gene
    done
```

What to verify:

- each run succeeds,
- `report_ids.txt` now has one line per input file,
- report ids are stable and human-readable, and
- you did not need to decompress the gzip files first.

If you need sparse family-aware indexing instead, change only the
interpretation call:

```bash
cargo run --bin gentle_cli -- --state "$STATE" rna-reads interpret \
  "$SEQ_ID" "$FEATURE_ID" "$path" \
  --report-id "$report_id" \
  --profile nanopore_cdna_v1 \
  --format fasta \
  --scope all_overlapping_any_strand \
  --origin-mode multi_gene_sparse \
  --target-gene TP53 \
  --target-gene TP63 \
  --target-gene TP73
```

That is the better choice when family members could steal or rescue assignments.

## Step 4: Run Phase-2 Alignment for Every Saved Report

Now rerun phase-2 on all retained rows for every report:

```bash
while IFS= read -r report_id; do
  cargo run --bin gentle_cli -- --state "$STATE" rna-reads align-report \
    "$report_id" \
    --selection all
done < "$OUT/report_ids.txt"
```

Why `--selection all` matters:

- it lets phase-2 score every retained saved row rather than only the stricter
  `seed_passed` subset
- that is usually the better first pass when the goal is cohort annotation and
  not only conservative filtering

## Step 5: Export One Combined Target-Gene Sample Sheet

This version uses the exact saved report ids rather than `--seq-id`, which is
safer if the state file contains anything else.

```bash
cmd=(
  cargo run --bin gentle_cli -- --state "$STATE"
  rna-reads export-sample-sheet
  "$OUT/${GENE_ID}_batch_sample_sheet.tsv"
  --gene "$GENE_ID"
  --complete-rule near
)

while IFS= read -r report_id; do
  cmd+=(--report-id "$report_id")
done < "$OUT/report_ids.txt"

"${cmd[@]}"
```

What to verify:

- the command succeeds with `report_count=20`,
- the written TSV has one header row plus one data row per sample,
- the output path is exactly:
  `"$OUT/${GENE_ID}_batch_sample_sheet.tsv"`

If you kept this state file dedicated to this batch run, the shorter variant is
also valid:

```bash
cargo run --bin gentle_cli -- --state "$STATE" rna-reads export-sample-sheet \
  "$OUT/${GENE_ID}_batch_sample_sheet.tsv" \
  --seq-id "$SEQ_ID" \
  --gene "$GENE_ID" \
  --complete-rule near
```

## Step 6: Read the Columns That Matter

The final sample sheet is designed so you can open it in R, Python, Excel,
LibreOffice, or a downstream stats pipeline without re-deriving the core
cohorts.

Columns to focus on first:

- `read_count_total`
  - total reads seen in the input file
- `read_count_aligned`
  - retained reads with a phase-2 best mapping
- `mean_read_length_bp`
  - overall mean cDNA length for the sample
- `gene_support_accepted_target_count`
  - number of aligned retained reads assigned to the requested target gene
- `gene_support_accepted_target_fraction_total`
  - accepted target reads divided by all reads in that sample
- `gene_support_accepted_target_fraction_aligned`
  - accepted target reads divided by aligned retained reads
- `gene_support_fragment_count`
  - assigned target-gene reads below the chosen complete rule
- `gene_support_complete_count`
  - assigned target-gene reads meeting the chosen complete rule
- `gene_support_mean_assigned_read_length_bp`
  - mean cDNA length among reads assigned to the target gene
- `gene_support_exon_pair_support_json`
  - ordered exon-pair co-presence rows for the accepted target-gene cohort
- `gene_support_direct_transition_support_json`
  - neighboring exon-transition rows only

Important interpretation detail:

- `gene_support_exon_pair_support_json` includes skipped ordered pairs such as
  `1->3`
- `gene_support_direct_transition_support_json` includes only neighboring exon
  steps such as `1->2`

That separation is what lets you distinguish broad exon co-presence from direct
junction evidence.

## Step 7: Audit One Interesting Sample in Detail

Once the sample sheet tells you which sample is interesting, drill into that
one report without redoing the whole batch.

Pick one saved report id, then write both the aggregate summary and the
row-level audit:

```bash
REPORT_ID=$(head -n 1 "$OUT/report_ids.txt")

cargo run --bin gentle_cli -- --state "$STATE" rna-reads summarize-gene-support \
  "$REPORT_ID" \
  --gene "$GENE_ID" \
  --complete-rule near \
  --output "$OUT/${REPORT_ID}_gene_summary.json"

cargo run --bin gentle_cli -- --state "$STATE" rna-reads inspect-gene-support \
  "$REPORT_ID" \
  --gene "$GENE_ID" \
  --complete-rule near \
  --cohort accepted \
  --output "$OUT/${REPORT_ID}_accepted_rows.json"
```

What these two files give you:

- `*_gene_summary.json`
  - compact cohort counts plus exon/exon-pair/direct-transition tables
- `*_accepted_rows.json`
  - exact read membership, full-length flags, mapped exon ordinals, exon pairs,
    direct transitions, and score/identity/coverage fields

## Optional: Compare `near` vs `strict` vs `exact`

If you want to separate looser fragments from more complete molecules, rerun
the sample-sheet export with a stricter complete rule:

```bash
cargo run --bin gentle_cli -- --state "$STATE" rna-reads export-sample-sheet \
  "$OUT/${GENE_ID}_batch_sample_sheet.strict.tsv" \
  --seq-id "$SEQ_ID" \
  --gene "$GENE_ID" \
  --complete-rule strict
```

Interpretation:

- `near`
  - default practical cohort split
- `strict`
  - `near` plus both transcript ends close and identity above threshold
- `exact`
  - `100%` template coverage only

## Engine / Shell Mapping

This tutorial exercises these shared routes and engine operations:

| Tutorial step | Shared route | Engine operation |
| --- | --- | --- |
| Interpret one gzipped FASTA | `rna-reads interpret ...` | `InterpretRnaReads` |
| Align retained report rows | `rna-reads align-report ...` | `AlignRnaReadReport` |
| Export combined sample sheet | `rna-reads export-sample-sheet ...` | `ExportRnaReadSampleSheet` |
| Summarize one sample | `rna-reads summarize-gene-support ...` | `SummarizeRnaReadGeneSupport` |
| Audit one sample | `rna-reads inspect-gene-support ...` | `InspectRnaReadGeneSupport` |

If you prefer the shared shell form, the same command family is also available
through:

```bash
cargo run --bin gentle_cli -- --state "$STATE" shell 'rna-reads ...'
```

## Why This Tutorial Matters

This is the current practical bridge between raw multi-sample cDNA input files
and downstream gene-centric comparison.

It gives you one deterministic route for:

- preserving one saved report per raw input file,
- aligning them through the shared engine path,
- exporting one consistent cohort table across all samples, and
- keeping the target-gene read set auditable per sample without reimplementing
  the cohort logic outside GENtle.

## Related Reading

- TP53 sparse-origin setup chapter:
  [`docs/tutorial/generated/chapters/12_tp53_multi_gene_sparse_mapping_online.md`](./generated/chapters/12_tp53_multi_gene_sparse_mapping_online.md)
- CLI reference for RNA-read commands:
  [`docs/cli.md`](../cli.md)
- tutorial landing page:
  [`docs/tutorial/README.md`](./README.md)
