# TP73 Pancreas Nanopore cDNA Benchmark Runbook

Purpose: run one clean, reproducible TP73 benchmark on the pancreatic cancer
Nanopore cDNA data from `GSE293661` / `SRR32957126`, using TP73 positive
transcript controls and TP53/TP63 negative controls to derive seed thresholds
before RNA-read mapping.

This runbook is written for `gentle.functional.domains` and assumes a local
GENtle checkout at `/home/clawbio/GENtle`.

For a fixed-parameter multi-accession cohort run after this single-run proof,
see `docs/tp73_pancreas_cohort_batch_runbook.md`.

The preferred execution path is the guarded helper script:

```bash
cd /home/clawbio/GENtle
scripts/tp73_pancreas_rna_mapping.sh strict-start
```

Run it inside `screen` or `tmux`. The script defines all paths, writes
`run.env`, refuses to start while another RNA-read mapping process is active,
and refuses to overwrite any state/checkpoint/log/report path that already
exists for the selected run. The manual blocks below document what the script
does internally and remain useful for debugging or adapting the workflow.

## 1. Environment

```bash
set -euo pipefail

cd /home/clawbio/GENtle

export RUN="${RUN:-SRR32957126}"
export GENE_ID="${GENE_ID:-TP73}"
export SEQ_ID="${SEQ_ID:-tp73_ncbi}"
export PAN_ROOT="${PAN_ROOT:-/home/clawbio/data/SRA/Pancreas_Carcinoma_Nanopore_GSE293661}"
export WORK="${WORK:-/home/clawbio/work/tp73_pancreas_benchmark_from_scratch}"
export STATE="${STATE:-$WORK/tp73_pancreas_from_scratch.gentle.json}"
export REPORT_ID="${REPORT_ID:-tp73_pancreas_${RUN}_from_scratch}"
export THREADS="${THREADS:-8}"

mkdir -p \
  "$PAN_ROOT" \
  "$WORK"/{fastq,logs,tmp,reports,manifests,checkpoints,post_interpret/{json,tsv,svg,fasta}}

cargo build --release --bin gentle_cli
GENTLE=(/home/clawbio/GENtle/target/release/gentle_cli --state "$STATE" --progress-stderr)
```

Recommended tool preflight:

```bash
for tool in prefetch vdb-validate fasterq-dump pigz seqkit jq /usr/bin/time; do
  command -v "$tool" >/dev/null || { echo "Missing required tool: $tool" >&2; exit 1; }
done
```

On Debian, `prefetch`, `vdb-validate`, and `fasterq-dump` are from
`sra-toolkit`; `pigz` and `time` may need explicit installation.

## 2. Retrieve And Validate SRA

```bash
prefetch --progress --resume yes --verify yes --max-size u "$RUN" \
  --output-directory "$PAN_ROOT" \
  2>&1 | tee "$WORK/logs/$RUN.prefetch.log"

vdb-validate "$PAN_ROOT/$RUN/$RUN.sra" \
  > "$WORK/logs/$RUN.vdb-validate.log" \
  2>&1
```

## 3. Convert To FASTA

`fasterq-dump --stdout` cannot be used with split-output modes. Write FASTQ to
disk first, then compress and convert.

```bash
if [ ! -s "$WORK/fastq/$RUN.split-spot.fastq.gz" ]; then
  /usr/bin/time -v fasterq-dump "$PAN_ROOT/$RUN/$RUN.sra" \
    --threads "$THREADS" \
    --temp "$WORK/tmp" \
    --outdir "$WORK/fastq" \
    --split-spot \
    > "$WORK/logs/$RUN.fasterq.stdout.log" \
    2> "$WORK/logs/$RUN.fasterq.stderr.log"

  if [ -s "$WORK/fastq/$RUN.fastq" ]; then
    mv "$WORK/fastq/$RUN.fastq" "$WORK/fastq/$RUN.split-spot.fastq"
  fi

  pigz -p "$THREADS" "$WORK/fastq/$RUN.split-spot.fastq"
fi

export READ_FASTQ="$WORK/fastq/$RUN.split-spot.fastq.gz"
export READ_FASTA="$WORK/fastq/$RUN.split-spot.fasta.gz"

gzip -t "$READ_FASTQ"
seqkit stats -a -T "$READ_FASTQ" \
  > "$WORK/logs/$RUN.split-spot.seqkit-stats.tsv"

if [ ! -s "$READ_FASTA" ]; then
  seqkit fq2fa -j "$THREADS" "$READ_FASTQ" | pigz -p "$THREADS" > "$READ_FASTA"
fi

gzip -t "$READ_FASTA"
seqkit stats -a -T "$READ_FASTA" \
  > "$WORK/logs/$RUN.split-spot.fasta.seqkit-stats.tsv"
```

## 4. Initialize Fresh GENtle State

```bash
cat > "$WORK/manifests/load_tp73_ncbi.workflow.json" <<'JSON'
{
  "run_id": "load_tp73_ncbi_for_pancreas_benchmark",
  "ops": [
    {
      "LoadFile": {
        "path": "test_files/tp73.ncbi.gb",
        "as_id": "tp73_ncbi"
      }
    }
  ]
}
JSON

"${GENTLE[@]}" workflow @"$WORK/manifests/load_tp73_ncbi.workflow.json" \
  > "$WORK/reports/load_tp73_ncbi.command.json" \
  2> "$WORK/logs/load_tp73_ncbi.stderr.log"

export SEED_FEATURE_ID="$(
  "${GENTLE[@]}" features query "$SEQ_ID" \
    --kind gene --label "$GENE_ID" --limit 1 \
    2> "$WORK/logs/query_${GENE_ID}_gene.stderr.log" \
    | jq -r '.rows[0].feature_id'
)"

test -n "$SEED_FEATURE_ID"
test "$SEED_FEATURE_ID" != "null"
echo "Using SEED_FEATURE_ID=$SEED_FEATURE_ID"
```

## 5. Derive TP73 Thresholds From Positive And Control Genes

Positive controls are hard must-pass transcript variants. They include Ensembl
TP73 transcripts plus the curated TP73 DeltaN and D_{Ex2,3}Np73 supplements so
non-standard TP73 isoforms participate in seed-filter optimization. TP53 and
TP63 are negative controls; optimized seed filters are rejected if either
control family passes above the configured ceiling.

```bash
"${GENTLE[@]}" rna-reads preflight-isoforms "$SEQ_ID" "$SEED_FEATURE_ID" \
  --scope all_overlapping_any_strand \
  --must-pass-transcript-fasta test_files/fixtures/mapping/ensembl_human_tp73_all.fasta \
  --must-pass-transcript-fasta test_files/fixtures/mapping/ensembl_chimp_tp73_all.fasta \
  --must-pass-transcript-fasta data/resources/tp73_dn_ena_transcripts.fasta \
  --must-pass-transcript-fasta data/resources/tp73_delta_ex2_3_refseq_derived_transcripts.fasta \
  --control-transcript-fasta test_files/fixtures/mapping/ensembl_human_tp53_all.fasta \
  --control-transcript-fasta test_files/fixtures/mapping/ensembl_human_tp63_all.fasta \
  --optimize-parameters \
  --max-control-match-probability 0.0 \
  > "$WORK/reports/$REPORT_ID.preflight.json" \
  2> "$WORK/logs/$REPORT_ID.preflight.stderr.log"

jq '{
  target: [.target_passed_transcript_count, .target_transcript_count],
  positives: [.positive_control_passed_transcript_count, .positive_control_transcript_count],
  controls: [.control_summaries[] | {
    control_id,
    transcript_count,
    passed_transcript_count,
    weighted_pass_probability,
    best_transcript_id,
    best_raw_hit_fraction,
    best_weighted_hit_fraction,
    best_unique_matched_kmers
  }],
  threshold_recommendation
}' "$WORK/reports/$REPORT_ID.preflight.json" \
  > "$WORK/reports/$REPORT_ID.preflight.summary.json"

SEED_FRAGMENT="$(
  jq -r '.threshold_recommendation.seed_filter_cli_fragment' \
    "$WORK/reports/$REPORT_ID.preflight.json"
)"

if [ -z "$SEED_FRAGMENT" ] || [ "$SEED_FRAGMENT" = "null" ]; then
  echo "No threshold recommendation was produced; inspect $WORK/reports/$REPORT_ID.preflight.summary.json" >&2
  exit 1
fi

read -r -a SEED_ARGS <<< "$SEED_FRAGMENT"
printf 'Seed args: %q ' "${SEED_ARGS[@]}" | tee "$WORK/logs/$REPORT_ID.seed_args.log"
printf '\n' | tee -a "$WORK/logs/$REPORT_ID.seed_args.log"
```

## 6. Interpret Reads

This phase stores only seed-passed rows but keeps full counters and checkpoints.

```bash
/usr/bin/time -v "${GENTLE[@]}" rna-reads interpret "$SEQ_ID" "$SEED_FEATURE_ID" "$READ_FASTA" \
  --report-id "$REPORT_ID" \
  --report-mode seed_passed_only \
  --checkpoint-path "$WORK/checkpoints/$REPORT_ID.checkpoint.json" \
  --checkpoint-every-reads 100000 \
  --profile nanopore_cdna_v1 \
  --format fasta \
  --scope all_overlapping_any_strand \
  --origin-mode single_gene \
  "${SEED_ARGS[@]}" \
  --align-band-bp 24 \
  --align-min-identity 0.55 \
  --max-secondary-mappings 3 \
  > "$WORK/reports/$REPORT_ID.interpret.command.json" \
  2> "$WORK/logs/$REPORT_ID.interpret.stderr.log"
```

Progress monitor from another shell:

```bash
while sleep 60; do
  date
  jq '{
    reads_processed,
    read_count_seed_passed,
    read_count_aligned,
    input_bytes_processed,
    input_bytes_total,
    cumulative_matched_kmers,
    cumulative_tested_kmers
  }' "$WORK/checkpoints/$REPORT_ID.checkpoint.json" 2>/dev/null || true
done
```

## 7. Align Retained Hits And Export Evidence

```bash
"${GENTLE[@]}" rna-reads align-report "$REPORT_ID" \
  --selection all \
  --align-band-bp 24 \
  --align-min-identity 0.55 \
  --max-secondary-mappings 3 \
  > "$WORK/post_interpret/json/$REPORT_ID.align_report.json" \
  2> "$WORK/logs/$REPORT_ID.align_report.stderr.log"

"${GENTLE[@]}" rna-reads show-report "$REPORT_ID" \
  > "$WORK/post_interpret/json/$REPORT_ID.show_report.json" \
  2> "$WORK/logs/$REPORT_ID.show_report.stderr.log"

"${GENTLE[@]}" rna-reads inspect-alignments "$REPORT_ID" \
  --selection aligned --limit 200 --sort score \
  > "$WORK/post_interpret/json/$REPORT_ID.aligned.top200.json" \
  2> "$WORK/logs/$REPORT_ID.aligned.top200.stderr.log"

"${GENTLE[@]}" rna-reads summarize-gene-support "$REPORT_ID" \
  --gene "$GENE_ID" --complete-rule near \
  --output "$WORK/post_interpret/json/$REPORT_ID.${GENE_ID}.gene_support.near.json" \
  > "$WORK/post_interpret/json/$REPORT_ID.${GENE_ID}.gene_support.near.command.json" \
  2> "$WORK/logs/$REPORT_ID.${GENE_ID}.gene_support.near.stderr.log"

"${GENTLE[@]}" rna-reads export-alignments-tsv "$REPORT_ID" \
  "$WORK/post_interpret/tsv/$REPORT_ID.alignments.aligned.tsv" \
  --selection aligned --subset-spec "from_scratch_tp73_pancreas_${RUN}" \
  > "$WORK/post_interpret/json/$REPORT_ID.export_alignments.command.json" \
  2> "$WORK/logs/$REPORT_ID.export_alignments.stderr.log"

"${GENTLE[@]}" rna-reads export-paths-tsv "$REPORT_ID" \
  "$WORK/post_interpret/tsv/$REPORT_ID.paths.aligned.tsv" \
  --selection aligned --subset-spec "from_scratch_tp73_pancreas_${RUN}" \
  > "$WORK/post_interpret/json/$REPORT_ID.export_paths.command.json" \
  2> "$WORK/logs/$REPORT_ID.export_paths.stderr.log"

"${GENTLE[@]}" rna-reads export-abundance-tsv "$REPORT_ID" \
  "$WORK/post_interpret/tsv/$REPORT_ID.abundance.aligned.tsv" \
  --selection aligned --subset-spec "from_scratch_tp73_pancreas_${RUN}" \
  > "$WORK/post_interpret/json/$REPORT_ID.export_abundance.command.json" \
  2> "$WORK/logs/$REPORT_ID.export_abundance.stderr.log"

"${GENTLE[@]}" rna-reads export-target-quality "$REPORT_ID" \
  "$WORK/post_interpret/svg/$REPORT_ID.${GENE_ID}.target_quality.svg" \
  --gene "$GENE_ID" --complete-rule near \
  > "$WORK/post_interpret/json/$REPORT_ID.${GENE_ID}.target_quality.command.json" \
  2> "$WORK/logs/$REPORT_ID.${GENE_ID}.target_quality.stderr.log"
```

## 8. Final Summary For Review

```bash
jq '{
  report_id: .report.report_id,
  read_count_total: .report.read_count_total,
  read_count_seed_passed: .report.read_count_seed_passed,
  read_count_aligned: .report.read_count_aligned,
  retained_count_msa_eligible: .report.retained_count_msa_eligible,
  read_lengths: .report.read_length_distributions
}' "$WORK/post_interpret/json/$REPORT_ID.show_report.json" \
  > "$WORK/reports/$REPORT_ID.final_summary.json"

{
  echo "# TP73 pancreas benchmark evidence bundle"
  echo
  echo "- run: $RUN"
  echo "- state: $STATE"
  echo "- report: $REPORT_ID"
  echo "- read FASTA: $READ_FASTA"
  echo "- preflight: $WORK/reports/$REPORT_ID.preflight.json"
  echo "- final summary: $WORK/reports/$REPORT_ID.final_summary.json"
  echo "- top alignments: $WORK/post_interpret/json/$REPORT_ID.aligned.top200.json"
  echo "- alignments TSV: $WORK/post_interpret/tsv/$REPORT_ID.alignments.aligned.tsv"
  echo "- target quality SVG: $WORK/post_interpret/svg/$REPORT_ID.${GENE_ID}.target_quality.svg"
} > "$WORK/reports/$REPORT_ID.evidence_bundle.md"
```

The two most important proof points are:

- the preflight summary shows all TP73 positives pass while TP53/TP63 controls
  stay below the configured pass ceiling
- the final report shows how many pancreatic cDNA reads seed-pass and align
  coherently to TP73 under those derived thresholds
