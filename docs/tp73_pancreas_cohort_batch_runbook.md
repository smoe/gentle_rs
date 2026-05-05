# TP73 Pancreas Nanopore cDNA Cohort Batch Runbook

Purpose: compare TP73 abundance and terminal isoform/path evidence across a
small pancreatic cancer Nanopore cDNA cohort while keeping the method fixed
across samples. TP73 is the benchmark target because it is experimentally
actionable for the current GENtle team; the batch pattern itself should remain
gene-agnostic.

This runbook builds on `docs/tp73_pancreas_benchmark_runbook.md`. Use the
single-run runbook first when validating a new machine. Use this cohort runbook
once the toolchain, SRA retrieval, FASTA conversion, RNA-read interpretation,
and phase-2 alignment are known to work for one accession.

## Scripted Path

Prefer the repository helper for the live pancreatic TP73 cohort run:

```bash
cd /home/clawbio/GENtle

scripts/tp73_pancreas_cohort.sh init

# Optional: fold the already harvested SRR32957126 strict run into the cohort
# table without rerunning the expensive mapping step.
scripts/tp73_pancreas_cohort.sh adopt-existing \
  /home/clawbio/work/tp73_pancreas_benchmark/strict_SRR32957126_20260503T003857Z \
  SRR32957126 SRR32957126

scripts/tp73_pancreas_cohort.sh status
scripts/tp73_pancreas_cohort.sh all
```

The serial `all` path walks sample-by-sample. For each manifest row it creates
or reuses:

- `"$WORK_ROOT/cohort/manifests/pancreas_runs.tsv"` with `SRR32957124` through
  `SRR32957129`
- SRA files under `"$PAN_ROOT"`
- shared FASTA files under `"$WORK_ROOT/fastq"`
- per-run states/logs under `"$WORK_ROOT/cohort/runs/$RUN"`
- a figure-ready cohort table at
  `"$WORK_ROOT/cohort/reports/tp73_pancreas_cohort.summary.tsv"`
- merged GENtle sample-sheet output at
  `"$WORK_ROOT/cohort/reports/tp73_pancreas_cohort.sample_sheet.near.tsv"`

The helper keeps the RNA mapping itself serial and deterministic. If disk space
becomes tight, set `DROP_FASTQ_AFTER_FASTA=1` only after confirming that the
FASTA conversion and `seqkit` stats are already written for that run.

After `summarize`, create a figure-source TSV and render the cohort overview:

```bash
COHORT_ROOT="${COHORT_ROOT:-/home/clawbio/work/tp73_pancreas_benchmark/cohort}"
mkdir -p "$COHORT_ROOT/figures"

jq -r '
  [
    "run_accession","sample_id","sample_name",
    "total_reads",
    "all_q0_bp","all_q25_bp","all_q50_bp","all_q75_bp","all_q100_bp","all_mean_bp",
    "strict_seed_passed_reads",
    "strict_seed_passed_max_read_bp",
    "strict_seed_passed_mean_read_bp",
    "accepted_tp73_reads",
    "accepted_tp73_per_million",
    "accepted_tp73_max_read_bp",
    "accepted_tp73_mean_read_bp"
  ],
  (.[] | [
    .run_accession,
    .sample_id,
    .sample_name,
    .read_count_total,
    .read_length_summaries.all_reads.quantiles_bp.q0,
    .read_length_summaries.all_reads.quantiles_bp.q25,
    .read_length_summaries.all_reads.quantiles_bp.q50,
    .read_length_summaries.all_reads.quantiles_bp.q75,
    .read_length_summaries.all_reads.quantiles_bp.q100,
    .read_length_summaries.all_reads.mean_bp,
    .strict_seed_passed_reads,
    .read_length_summaries.strict_seed_passed.quantiles_bp.q100,
    .read_length_summaries.strict_seed_passed.mean_bp,
    .gene_support.accepted_target_count,
    (if .read_count_total > 0 then (.gene_support.accepted_target_count * 1000000 / .read_count_total) else null end),
    .gene_support.accepted_target_read_lengths.max_length_bp,
    .gene_support.accepted_target_read_lengths.mean_length_bp
  ])
  | @tsv
' "$COHORT_ROOT/reports/tp73_pancreas_cohort.final_summaries.json" \
  > "$COHORT_ROOT/figures/tp73_pancreas_figure_source.tsv"

python3 scripts/plot_tp73_pancreas_cohort.py \
  "$COHORT_ROOT/figures/tp73_pancreas_figure_source.tsv" \
  --output "$COHORT_ROOT/figures/tp73_pancreas_cohort_overview.svg"
```

The SVG plot is dependency-free. If PNG/PDF conversion is needed on Debian,
install `librsvg2-bin` and run `rsvg-convert`, or use Inkscape.

## Generic Follow-On Gene Screens

Once the shared pancreatic FASTA files already exist, use the gene-agnostic
helper for TP53/TP63/E2F1/POU2F1-style follow-up screens instead of cloning the
TP73 wrapper again. It resolves the HUGO symbol through NCBI Gene, downloads a
compact GenBank locus, loads that locus into a base GENtle state, then copies
that state per sample so independent workers can run safely in parallel.

```bash
cd /home/clawbio/GENtle

scripts/pancreas_gene_rna_screen.sh run E2F1 --jobs 4
scripts/pancreas_gene_rna_screen.sh run POU2F1 --jobs 4
```

The helper prints a monitor command as soon as the run directory is initialized.
The same status table can also be queried manually:

```bash
RUN_ROOT=/home/clawbio/work/tp73_pancreas_benchmark/gene_screens/e2f1_pancreas_screen_YYYYMMDDTHHMMSSZ

scripts/pancreas_gene_rna_screen.sh status "$RUN_ROOT" | column -t -s $'\t'
scripts/pancreas_gene_rna_screen.sh monitor "$RUN_ROOT" 120
scripts/pancreas_gene_rna_screen.sh summarize "$RUN_ROOT"

python3 scripts/plot_pancreas_gene_screen.py \
  "$RUN_ROOT/figures/e2f1_pancreas_figure_source.tsv" \
  --gene E2F1 \
  --output "$RUN_ROOT/figures/e2f1_pancreas_overview.svg"
```

Useful defaults:

- `--jobs N` runs up to `N` sample workers at once; one worker is one
  single-core `gentle_cli` RNA-read mapping process.
- `--genbank-path PATH` bypasses NCBI if a curated locus file is already
  available.
- `--must-pass-transcript-fasta`, `--positive-transcript-fasta`, and
  `--control-transcript-fasta` are repeatable and feed the same isoform
  preflight contract used for TP73 controls.
- If no gene-specific transcript fixture exists, the helper still records the
  preflight report but falls back to the conservative default seed-filter
  fragment.
- `--concatemer-limit N` is off by default because first-pass abundance screens
  should not spend hours on transcriptome-scale concatemer audits unless that
  evidence is specifically needed.

Primary outputs live under the timestamped run root:

- `manifests/*_pancreas_inputs.tsv`
- `reports/*.preflight.json` and `reports/*.preflight.summary.json`
- per-sample `runs/$RUN/checkpoints/*.checkpoint.json`
- per-sample `runs/$RUN/post_interpret/{json,tsv,svg}/...`
- `reports/*_pancreas.summary.tsv`
- `figures/*_pancreas_figure_source.tsv`
- `figures/*_pancreas_overview.svg` when rendered with
  `scripts/plot_pancreas_gene_screen.py`

## 1. Batch Principle

For abundance variation, do not tune thresholds separately per sample after
seeing the sample results.

Use one fixed seed/alignment method for every run:

- derive TP73 seed filters once from transcript controls, not from one sample's
  observed abundance
- use the same `SEQ_ID`, `SEED_FEATURE_ID`, `SEED_ARGS`, profile, scope,
  origin mode, alignment band, identity threshold, and secondary-mapping limit
- keep one separate state file per run when running in parallel
- compare normalized counts and path evidence only after every sample was
  processed with the same settings

This makes low TP73 abundance interpretable as biology or library behavior
rather than a side effect of per-sample parameter choices.

## 2. Environment

```bash
set -euo pipefail

cd /home/clawbio/GENtle

export GENE_ID="${GENE_ID:-TP73}"
export SEQ_ID="${SEQ_ID:-tp73_ncbi}"
export PAN_ROOT="${PAN_ROOT:-/home/clawbio/data/SRA/Pancreas_Carcinoma_Nanopore_GSE293661}"
export WORK="${WORK:-/home/clawbio/work/tp73_pancreas_cohort_batch}"
export BASE_STATE="${BASE_STATE:-$WORK/base_tp73.gentle.json}"
export THREADS="${THREADS:-8}"

mkdir -p "$PAN_ROOT" "$WORK"/{logs,manifests,reports,runs}

cargo build --release --bin gentle_cli
GENTLE_BASE=(/home/clawbio/GENtle/target/release/gentle_cli --state "$BASE_STATE" --progress-stderr)
```

Tool preflight:

```bash
for tool in prefetch vdb-validate fasterq-dump pigz seqkit jq /usr/bin/time; do
  command -v "$tool" >/dev/null || { echo "Missing required tool: $tool" >&2; exit 1; }
done
```

## 3. Cohort Manifest

Keep the manifest intentionally small at first. Add cell-line/sample labels when
they are known; until then, use the SRA run accession as the sample id.

```bash
cat > "$WORK/manifests/pancreas_runs.tsv" <<'TSV'
run_accession	sample_id	note
SRR32957124	SRR32957124	pancreas_nanopore_cdna
SRR32957125	SRR32957125	pancreas_nanopore_cdna
SRR32957126	SRR32957126	pancreas_nanopore_cdna
SRR32957127	SRR32957127	pancreas_nanopore_cdna
SRR32957128	SRR32957128	pancreas_nanopore_cdna
SRR32957129	SRR32957129	pancreas_nanopore_cdna
TSV
```

If one accession is still partially downloaded or locked, remove it from the
first manifest and add it back later. A smaller completed cohort is better than
a blocked batch.

For new GENtle-native preparation runs, prefer the shared acquisition contract:
make a TSV with `sample_id` and `sra_accession`, then run
`reads acquire prepare MANIFEST.tsv --cache-dir "$PAN_ROOT" --work-dir "$WORK/read_acquisition"`.
`rna-reads batch-map --prepare-sra --read-cache-dir "$PAN_ROOT" --read-work-dir "$WORK/read_acquisition"`
uses the same layer before mapping. The older shell snippets below remain useful
when a cluster already has bespoke SRA/FASTQ/FASTA staging.

## 4. Prepare One Base State And One Seed Filter

The base state contains only the target sequence and shared preflight result.
Per-run states are copied from this file so cluster jobs never write to the
same mutable state concurrently.

```bash
cat > "$WORK/manifests/load_tp73_ncbi.workflow.json" <<'JSON'
{
  "run_id": "load_tp73_ncbi_for_pancreas_batch",
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

"${GENTLE_BASE[@]}" workflow @"$WORK/manifests/load_tp73_ncbi.workflow.json" \
  > "$WORK/reports/load_tp73_ncbi.command.json" \
  2> "$WORK/logs/load_tp73_ncbi.stderr.log"

export SEED_FEATURE_ID="$(
  "${GENTLE_BASE[@]}" features query "$SEQ_ID" \
    --kind gene --label "$GENE_ID" --limit 1 \
    2> "$WORK/logs/query_${GENE_ID}_gene.stderr.log" \
    | jq -r '.rows[0].feature_id'
)"

test -n "$SEED_FEATURE_ID"
test "$SEED_FEATURE_ID" != "null"
echo "Using SEED_FEATURE_ID=$SEED_FEATURE_ID"
```

Derive the fixed seed filter once:

```bash
"${GENTLE_BASE[@]}" rna-reads preflight-isoforms "$SEQ_ID" "$SEED_FEATURE_ID" \
  --scope all_overlapping_any_strand \
  --must-pass-transcript-fasta test_files/fixtures/mapping/ensembl_human_tp73_all.fasta \
  --must-pass-transcript-fasta test_files/fixtures/mapping/ensembl_chimp_tp73_all.fasta \
  --must-pass-transcript-fasta data/resources/tp73_dn_ena_transcripts.fasta \
  --must-pass-transcript-fasta data/resources/tp73_delta_ex2_3_refseq_derived_transcripts.fasta \
  --control-transcript-fasta test_files/fixtures/mapping/ensembl_human_tp53_all.fasta \
  --control-transcript-fasta test_files/fixtures/mapping/ensembl_human_tp63_all.fasta \
  --optimize-parameters \
  --max-control-match-probability 0.0 \
  > "$WORK/reports/tp73_cohort.preflight.json" \
  2> "$WORK/logs/tp73_cohort.preflight.stderr.log"

jq -r '.threshold_recommendation.seed_filter_cli_fragment' \
  "$WORK/reports/tp73_cohort.preflight.json" \
  > "$WORK/manifests/tp73_seed_args.txt"

if [ ! -s "$WORK/manifests/tp73_seed_args.txt" ] || \
   [ "$(cat "$WORK/manifests/tp73_seed_args.txt")" = "null" ]; then
  echo "No seed recommendation; inspect $WORK/reports/tp73_cohort.preflight.json" >&2
  exit 1
fi
```

## 5. Per-Run Worker Function

This function is written so it can be pasted into an interactive shell, run
sequentially, or wrapped by a cluster scheduler. It keeps all per-run artifacts
under `"$WORK/runs/$RUN"`.

```bash
run_one_pancreas_tp73() {
  local run="$1"
  local sample_id="${2:-$run}"
  local run_work="$WORK/runs/$run"
  local state="$run_work/$run.gentle.json"
  local safe_sample_id
  local report_id
  local read_fastq="$run_work/fastq/$run.split-spot.fastq.gz"
  local read_fasta="$run_work/fastq/$run.split-spot.fasta.gz"
  local seed_args

  safe_sample_id="$(printf '%s' "$sample_id" | tr -c 'A-Za-z0-9_.-' '_')"
  if [ "$safe_sample_id" = "$run" ]; then
    report_id="tp73_pancreas_${run}"
  else
    report_id="tp73_pancreas_${safe_sample_id}_${run}"
  fi

  mkdir -p "$run_work"/{fastq,logs,tmp,reports,checkpoints,post_interpret/{json,tsv,svg,fasta}}
  cp "$BASE_STATE" "$state"
  GENTLE=(/home/clawbio/GENtle/target/release/gentle_cli --state "$state" --progress-stderr)
  read -r -a seed_args < "$WORK/manifests/tp73_seed_args.txt"

  if [ ! -s "$PAN_ROOT/$run/$run.sra" ]; then
    prefetch --progress --resume yes --verify yes --max-size u "$run" \
      --output-directory "$PAN_ROOT" \
      2>&1 | tee "$run_work/logs/$run.prefetch.log"
  fi

  vdb-validate "$PAN_ROOT/$run/$run.sra" \
    > "$run_work/logs/$run.vdb-validate.log" \
    2>&1

  if [ ! -s "$read_fastq" ]; then
    /usr/bin/time -v fasterq-dump "$PAN_ROOT/$run/$run.sra" \
      --threads "$THREADS" \
      --temp "$run_work/tmp" \
      --outdir "$run_work/fastq" \
      --split-spot \
      > "$run_work/logs/$run.fasterq.stdout.log" \
      2> "$run_work/logs/$run.fasterq.stderr.log"

    if [ -s "$run_work/fastq/$run.fastq" ]; then
      mv "$run_work/fastq/$run.fastq" "$run_work/fastq/$run.split-spot.fastq"
    fi

    pigz -p "$THREADS" "$run_work/fastq/$run.split-spot.fastq"
  fi

  gzip -t "$read_fastq"
  seqkit stats -a -T "$read_fastq" \
    > "$run_work/logs/$run.split-spot.seqkit-stats.tsv"

  if [ ! -s "$read_fasta" ]; then
    seqkit fq2fa -j "$THREADS" "$read_fastq" | pigz -p "$THREADS" > "$read_fasta"
  fi

  gzip -t "$read_fasta"
  seqkit stats -a -T "$read_fasta" \
    > "$run_work/logs/$run.split-spot.fasta.seqkit-stats.tsv"

  /usr/bin/time -v "${GENTLE[@]}" rna-reads interpret "$SEQ_ID" "$SEED_FEATURE_ID" "$read_fasta" \
    --report-id "$report_id" \
    --report-mode seed_passed_only \
    --checkpoint-path "$run_work/checkpoints/$report_id.checkpoint.json" \
    --checkpoint-every-reads 100000 \
    --profile nanopore_cdna_v1 \
    --format fasta \
    --scope all_overlapping_any_strand \
    --origin-mode single_gene \
    "${seed_args[@]}" \
    --align-band-bp 24 \
    --align-min-identity 0.55 \
    --max-secondary-mappings 3 \
    > "$run_work/reports/$report_id.interpret.command.json" \
    2> "$run_work/logs/$report_id.interpret.stderr.log"

  "${GENTLE[@]}" rna-reads align-report "$report_id" \
    --selection all \
    --align-band-bp 24 \
    --align-min-identity 0.55 \
    --max-secondary-mappings 3 \
    > "$run_work/post_interpret/json/$report_id.align_report.json" \
    2> "$run_work/logs/$report_id.align_report.stderr.log"

  "${GENTLE[@]}" rna-reads show-report "$report_id" \
    > "$run_work/post_interpret/json/$report_id.show_report.json" \
    2> "$run_work/logs/$report_id.show_report.stderr.log"

  "${GENTLE[@]}" rna-reads summarize-gene-support "$report_id" \
    --gene "$GENE_ID" --complete-rule near \
    --output "$run_work/post_interpret/json/$report_id.${GENE_ID}.gene_support.near.json" \
    > "$run_work/post_interpret/json/$report_id.${GENE_ID}.gene_support.near.command.json" \
    2> "$run_work/logs/$report_id.${GENE_ID}.gene_support.near.stderr.log"

  "${GENTLE[@]}" rna-reads export-sample-sheet \
    "$run_work/post_interpret/tsv/$report_id.sample_sheet.near.tsv" \
    --report-id "$report_id" \
    --gene "$GENE_ID" \
    --complete-rule near \
    > "$run_work/post_interpret/json/$report_id.sample_sheet.near.command.json" \
    2> "$run_work/logs/$report_id.sample_sheet.near.stderr.log"

  "${GENTLE[@]}" rna-reads export-paths-tsv "$report_id" \
    "$run_work/post_interpret/tsv/$report_id.paths.aligned.tsv" \
    --selection aligned --subset-spec "pancreas_cohort_${run}" \
    > "$run_work/post_interpret/json/$report_id.export_paths.command.json" \
    2> "$run_work/logs/$report_id.export_paths.stderr.log"

  "${GENTLE[@]}" rna-reads export-abundance-tsv "$report_id" \
    "$run_work/post_interpret/tsv/$report_id.abundance.aligned.tsv" \
    --selection aligned --subset-spec "pancreas_cohort_${run}" \
    > "$run_work/post_interpret/json/$report_id.export_abundance.command.json" \
    2> "$run_work/logs/$report_id.export_abundance.stderr.log"

  "${GENTLE[@]}" rna-reads export-target-quality "$report_id" \
    "$run_work/post_interpret/svg/$report_id.${GENE_ID}.target_quality.svg" \
    --gene "$GENE_ID" --complete-rule near \
    > "$run_work/post_interpret/json/$report_id.${GENE_ID}.target_quality.command.json" \
    2> "$run_work/logs/$report_id.${GENE_ID}.target_quality.stderr.log"

  jq --arg run "$run" --arg sample_id "$sample_id" '{
    run: $run,
    sample_id: $sample_id,
    report_id: .report.report_id,
    read_count_total: .report.read_count_total,
    read_count_seed_passed: .report.read_count_seed_passed,
    read_count_aligned: .report.read_count_aligned,
    retained_count_msa_eligible: .report.retained_count_msa_eligible,
    read_lengths: .report.read_length_distributions
  }' "$run_work/post_interpret/json/$report_id.show_report.json" \
    > "$run_work/reports/$report_id.final_summary.json"
}
```

## 6. Run Sequentially First

Start sequentially until the full path is boring. This is kinder to disk space
and easier to debug than launching every accession at once.

```bash
tail -n +2 "$WORK/manifests/pancreas_runs.tsv" \
  | while IFS=$'\t' read -r run sample_id note; do
      echo "=== $run / $sample_id / $note ==="
      run_one_pancreas_tp73 "$run" "$sample_id"
    done
```

For a one-sample smoke test:

```bash
run_one_pancreas_tp73 SRR32957126 SRR32957126
```

## 7. Cluster/Parallel Use

Each run gets its own copied state file under `"$WORK/runs/$RUN"`, so parallel
jobs do not write to the same GENtle state. The shared resources are:

- SRA cache under `"$PAN_ROOT"`
- fixed base state `"$BASE_STATE"` copied at job start
- fixed seed args `"$WORK/manifests/tp73_seed_args.txt"`

On a cluster, submit one manifest row per job and call:

```bash
run_one_pancreas_tp73 "$RUN" "$SAMPLE_ID"
```

Avoid concurrent jobs for the same `RUN`, because `prefetch` and `fasterq-dump`
would compete for the same SRA/FASTQ artifacts. If a run is already downloaded
and converted, reruns reuse the existing `.sra`, `.fastq.gz`, and `.fasta.gz`
files.

## 8. Merge Cohort Outputs

Merge per-run sample sheets into one cohort table:

```bash
first_sheet="$(
  find "$WORK/runs" -path '*post_interpret/tsv/*.sample_sheet.near.tsv' \
    | sort \
    | head -n 1
)"

{
  head -n 1 "$first_sheet"
  find "$WORK/runs" -path '*post_interpret/tsv/*.sample_sheet.near.tsv' \
    | sort \
    | while read -r sheet; do
        tail -n +2 "$sheet"
      done
} > "$WORK/reports/tp73_pancreas_cohort.sample_sheet.near.tsv"
```

Merge compact JSON summaries:

```bash
jq -s '.' "$WORK"/runs/*/reports/*.final_summary.json \
  > "$WORK/reports/tp73_pancreas_cohort.final_summaries.json"
```

Collect path/abundance tables for downstream inspection:

```bash
find "$WORK/runs" -path '*post_interpret/tsv/*.paths.aligned.tsv' \
  | sort > "$WORK/manifests/path_tables.txt"

find "$WORK/runs" -path '*post_interpret/tsv/*.abundance.aligned.tsv' \
  | sort > "$WORK/manifests/abundance_tables.txt"
```

## 9. First Metrics To Compare

Use the merged sample sheet as the primary comparison table.

Start with these columns:

- `read_count_total`
- `mean_read_length_bp`
- `read_count_seed_passed`
- `seed_pass_fraction`
- `read_count_aligned`
- `aligned_fraction`
- `gene_support_accepted_target_count`
- `gene_support_accepted_target_fraction_total`
- `gene_support_accepted_target_fraction_aligned`
- `gene_support_fragment_count`
- `gene_support_complete_count`
- `gene_support_mean_assigned_read_length_bp`
- `gene_support_exon_pair_support_json`
- `gene_support_direct_transition_support_json`

For TP73 abundance, the most robust first normalization is accepted TP73 target
reads per total reads:

```text
TP73 accepted-per-million = gene_support_accepted_target_fraction_total * 1e6
```

For isoform questions, treat the current output as terminal-exon/path evidence,
not definitive isoform quantification. With median read lengths around 700 bp
and expected 3' bias, alpha/beta/gamma differences may be visible as terminal
path differences, but many reads will remain fragments.

## 10. Stop Rule

For the release proof, this batch is successful when it produces:

- one fixed preflight report and seed-argument file
- one completed state/report set per accession
- one merged sample sheet
- one compact JSON summary list
- path/abundance TSVs for every completed accession
- enough target-quality SVGs to manually inspect outliers

Do not let this batch become a new feature gate. De-novo motif discovery,
improved scoring, expression-stratified interpretation, and richer isoform
deconvolution belong after the release cut.
