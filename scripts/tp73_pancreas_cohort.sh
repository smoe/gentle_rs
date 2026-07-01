#!/usr/bin/env bash
# Orchestrate the TP73 pancreatic Nanopore cDNA cohort benchmark.
#
# The workflow is intentionally serial and state-isolated:
# - one shared base state holds TP73 and the fixed preflight thresholds,
# - each SRA run gets its own copied GENtle state,
# - SRA, FASTQ, FASTA, logs, reports, and final tables live in predictable paths,
# - reruns skip completed samples instead of silently duplicating work.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
DEFAULT_REPO_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"

GENE_ID="${GENE_ID:-TP73}"
SEQ_ID="${SEQ_ID:-tp73_ncbi}"
THREADS="${THREADS:-6}"
GENTLE_REPO="${GENTLE_REPO:-${GENTLE_RS:-$DEFAULT_REPO_ROOT}}"
GENTLE_BIN="${GENTLE_BIN:-$GENTLE_REPO/target/release/gentle_cli}"
WORK_ROOT="${WORK_ROOT:-${WORK:-/home/clawbio/work/tp73_pancreas_benchmark}}"
COHORT_ROOT="${COHORT_ROOT:-$WORK_ROOT/cohort}"
PAN_ROOT="${PAN_ROOT:-/home/clawbio/data/SRA/Pancreas_Carcinoma_Nanopore_GSE293661}"
FASTA_ROOT="${FASTA_ROOT:-$WORK_ROOT/fastq}"
MANIFEST="${MANIFEST:-$COHORT_ROOT/manifests/pancreas_runs.tsv}"
BASE_STATE="${BASE_STATE:-$COHORT_ROOT/base_tp73.gentle.json}"
SEED_ARGS_FILE="${SEED_ARGS_FILE:-$COHORT_ROOT/manifests/tp73_seed_args.txt}"
SRA_MAX_SIZE="${SRA_MAX_SIZE:-u}"
COMPLETE_RULE="${COMPLETE_RULE:-near}"
ALIGN_BAND_BP="${ALIGN_BAND_BP:-24}"
ALIGN_MIN_IDENTITY="${ALIGN_MIN_IDENTITY:-0.55}"
MAX_SECONDARY_MAPPINGS="${MAX_SECONDARY_MAPPINGS:-3}"
DROP_FASTQ_AFTER_FASTA="${DROP_FASTQ_AFTER_FASTA:-0}"

usage() {
  cat <<'EOF'
Usage:
  scripts/tp73_pancreas_cohort.sh init
  scripts/tp73_pancreas_cohort.sh download
  scripts/tp73_pancreas_cohort.sh prepare-fasta
  scripts/tp73_pancreas_cohort.sh run
  scripts/tp73_pancreas_cohort.sh one RUN [SAMPLE_ID]
  scripts/tp73_pancreas_cohort.sh adopt-existing RUN_WORK [RUN] [SAMPLE_ID]
  scripts/tp73_pancreas_cohort.sh summarize
  scripts/tp73_pancreas_cohort.sh status
  scripts/tp73_pancreas_cohort.sh all

Default paths:
  GENTLE_REPO=/home/clawbio/GENtle
  GENTLE_BIN=/home/clawbio/GENtle/target/release/gentle_cli
  WORK_ROOT=/home/clawbio/work/tp73_pancreas_benchmark
  COHORT_ROOT=$WORK_ROOT/cohort
  PAN_ROOT=/home/clawbio/data/SRA/Pancreas_Carcinoma_Nanopore_GSE293661
  FASTA_ROOT=$WORK_ROOT/fastq
  MANIFEST=$COHORT_ROOT/manifests/pancreas_runs.tsv

Notes:
  - The default manifest contains SRR32957124..SRR32957129.
  - The pipeline is serial by default. Use one terminal/screen and let it walk.
  - Set DROP_FASTQ_AFTER_FASTA=1 only if you explicitly want to delete
    intermediate split-spot FASTQ.gz files after FASTA.gz creation.
EOF
}

die() {
  echo "ERROR: $*" >&2
  exit 1
}

info() {
  echo "INFO: $*" >&2
}

require_executable() {
  local tool="$1"
  if [[ "$tool" == */* ]]; then
    [ -x "$tool" ] || die "Missing executable: $tool"
  else
    command -v "$tool" >/dev/null 2>&1 || die "Missing executable in PATH: $tool"
  fi
}

require_file() {
  local path="$1"
  [ -s "$path" ] || die "Missing required file: $path"
}

ensure_dirs() {
  mkdir -p \
    "$PAN_ROOT" \
    "$FASTA_ROOT" \
    "$COHORT_ROOT"/{logs,manifests,reports,runs,tmp}
}

require_repo_inputs() {
  cd "$GENTLE_REPO" || die "Cannot cd to GENTLE_REPO: $GENTLE_REPO"
  require_executable "$GENTLE_BIN"
  require_executable jq
  require_executable gzip
  require_executable seqkit
  require_executable /usr/bin/time
  require_file test_files/tp73.ncbi.gb
  require_file test_files/fixtures/mapping/ensembl_human_tp73_all.fasta
  require_file test_files/fixtures/mapping/ensembl_chimp_tp73_all.fasta
  require_file test_files/fixtures/mapping/ensembl_human_tp53_all.fasta
  require_file test_files/fixtures/mapping/ensembl_human_tp63_all.fasta
  require_file data/resources/tp73_dn_ena_transcripts.fasta
  require_file data/resources/tp73_delta_ex2_3_refseq_derived_transcripts.fasta
}

require_sra_tools() {
  require_executable prefetch
  require_executable vdb-validate
}

require_conversion_tools() {
  require_executable fasterq-dump
  require_executable pigz
}

create_default_manifest() {
  ensure_dirs
  if [ -s "$MANIFEST" ]; then
    info "Manifest already exists: $MANIFEST"
    return 0
  fi
  cat > "$MANIFEST" <<'TSV'
run_accession	sample_id	sample_name	note
SRR32957124	SRR32957124	SRR32957124	pancreas_nanopore_cdna_GSE293661
SRR32957125	SRR32957125	SRR32957125	pancreas_nanopore_cdna_GSE293661
SRR32957126	SRR32957126	SRR32957126	pancreas_nanopore_cdna_GSE293661
SRR32957127	SRR32957127	SRR32957127	pancreas_nanopore_cdna_GSE293661
SRR32957128	SRR32957128	SRR32957128	pancreas_nanopore_cdna_GSE293661
SRR32957129	SRR32957129	SRR32957129	pancreas_nanopore_cdna_GSE293661
TSV
  info "Wrote manifest: $MANIFEST"
}

for_manifest_rows() {
  local callback="$1"
  require_file "$MANIFEST"
  while IFS=$'\t' read -r run sample_id sample_name note rest; do
    if [ -z "${run:-}" ] || [[ "$run" == \#* ]] || [ "$run" = "run_accession" ]; then
      continue
    fi
    sample_id="${sample_id:-$run}"
    sample_name="${sample_name:-$sample_id}"
    note="${note:-}"
    "$callback" "$run" "$sample_id" "$sample_name" "$note"
  done < "$MANIFEST"
}

safe_token() {
  printf '%s' "$1" | tr -c 'A-Za-z0-9_.-' '_'
}

report_id_for() {
  local run="$1"
  local sample_id="$2"
  local safe_sample_id
  safe_sample_id="$(safe_token "$sample_id")"
  if [ "$safe_sample_id" = "$run" ]; then
    printf 'tp73_pancreas_%s\n' "$run"
  else
    printf 'tp73_pancreas_%s_%s\n' "$safe_sample_id" "$run"
  fi
}

run_work_for() {
  printf '%s/runs/%s\n' "$COHORT_ROOT" "$1"
}

state_for() {
  local run="$1"
  printf '%s/%s.gentle.json\n' "$(run_work_for "$run")" "$run"
}

read_fastq_for() {
  local run="$1"
  printf '%s/%s.split-spot.fastq.gz\n' "$FASTA_ROOT" "$run"
}

read_fasta_for() {
  local run="$1"
  printf '%s/%s.split-spot.fasta.gz\n' "$FASTA_ROOT" "$run"
}

sra_path_for() {
  local run="$1"
  printf '%s/%s/%s.sra\n' "$PAN_ROOT" "$run" "$run"
}

init_shared() {
  ensure_dirs
  create_default_manifest
  require_repo_inputs

  local -a gentle_base
  gentle_base=("$GENTLE_BIN" --state "$BASE_STATE" --progress-stderr)

  if [ ! -s "$BASE_STATE" ]; then
    cat > "$COHORT_ROOT/manifests/load_tp73_ncbi.workflow.json" <<'JSON'
{
  "run_id": "load_tp73_ncbi_for_pancreas_cohort",
  "ops": [
    { "LoadFile": { "path": "test_files/tp73.ncbi.gb", "as_id": "tp73_ncbi" } }
  ]
}
JSON
    "${gentle_base[@]}" workflow @"$COHORT_ROOT/manifests/load_tp73_ncbi.workflow.json" \
      > "$COHORT_ROOT/reports/load_tp73_ncbi.command.json" \
      2> "$COHORT_ROOT/logs/load_tp73_ncbi.stderr.log"
  else
    info "Base state already exists: $BASE_STATE"
  fi

  local seed_feature_id
  seed_feature_id="$(
    "${gentle_base[@]}" features query "$SEQ_ID" \
      --kind gene --label "$GENE_ID" --limit 1 \
      2> "$COHORT_ROOT/logs/query_${GENE_ID}_gene.stderr.log" \
      | jq -r '.rows[0].feature_id'
  )"
  [ -n "$seed_feature_id" ] && [ "$seed_feature_id" != "null" ] \
    || die "Could not resolve $GENE_ID seed feature in $BASE_STATE"
  printf '%s\n' "$seed_feature_id" > "$COHORT_ROOT/manifests/seed_feature_id.txt"
  info "Seed feature ID: $seed_feature_id"

  if [ ! -s "$SEED_ARGS_FILE" ]; then
    "${gentle_base[@]}" rna-reads preflight-isoforms "$SEQ_ID" "$seed_feature_id" \
      --scope all_overlapping_any_strand \
      --must-pass-transcript-fasta test_files/fixtures/mapping/ensembl_human_tp73_all.fasta \
      --must-pass-transcript-fasta test_files/fixtures/mapping/ensembl_chimp_tp73_all.fasta \
      --must-pass-transcript-fasta data/resources/tp73_dn_ena_transcripts.fasta \
      --must-pass-transcript-fasta data/resources/tp73_delta_ex2_3_refseq_derived_transcripts.fasta \
      --control-transcript-fasta test_files/fixtures/mapping/ensembl_human_tp53_all.fasta \
      --control-transcript-fasta test_files/fixtures/mapping/ensembl_human_tp63_all.fasta \
      --optimize-parameters \
      --max-control-match-probability 0.0 \
      > "$COHORT_ROOT/reports/tp73_cohort.preflight.json" \
      2> "$COHORT_ROOT/logs/tp73_cohort.preflight.stderr.log"

    jq -r '.threshold_recommendation.seed_filter_cli_fragment' \
      "$COHORT_ROOT/reports/tp73_cohort.preflight.json" \
      > "$SEED_ARGS_FILE"
  else
    info "Seed args already exist: $SEED_ARGS_FILE"
  fi

  require_file "$SEED_ARGS_FILE"
  if [ "$(cat "$SEED_ARGS_FILE")" = "null" ]; then
    die "Seed args are null; inspect $COHORT_ROOT/reports/tp73_cohort.preflight.json"
  fi
}

download_one() {
  local run="$1"
  local sample_id="$2"
  local sample_name="$3"
  local note="$4"
  local run_work sra
  require_sra_tools
  run_work="$(run_work_for "$run")"
  sra="$(sra_path_for "$run")"
  mkdir -p "$run_work/logs"

  if [ ! -s "$sra" ]; then
    info "Downloading $run ($sample_id; $sample_name; $note)"
    prefetch --progress --resume yes --verify yes --max-size "$SRA_MAX_SIZE" "$run" \
      --output-directory "$PAN_ROOT" \
      2>&1 | tee "$run_work/logs/$run.prefetch.log"
  else
    info "SRA already present for $run: $sra"
  fi

  vdb-validate "$sra" > "$run_work/logs/$run.vdb-validate.log" 2>&1
}

download_manifest() {
  ensure_dirs
  create_default_manifest
  require_sra_tools
  for_manifest_rows download_one
}

prepare_fasta_one() {
  local run="$1"
  local sample_id="$2"
  local sample_name="$3"
  local note="$4"
  local run_work sra read_fastq read_fasta raw_fastq tmp_fasta
  require_executable gzip
  require_executable seqkit
  require_conversion_tools
  run_work="$(run_work_for "$run")"
  sra="$(sra_path_for "$run")"
  read_fastq="$(read_fastq_for "$run")"
  read_fasta="$(read_fasta_for "$run")"
  raw_fastq="$FASTA_ROOT/$run.split-spot.fastq"
  mkdir -p "$run_work/logs" "$run_work/tmp" "$FASTA_ROOT"

  if [ -s "$read_fasta" ]; then
    gzip -t "$read_fasta"
    seqkit stats -a -T "$read_fasta" \
      > "$run_work/logs/$run.split-spot.fasta.seqkit-stats.tsv"
    info "FASTA already present for $run: $read_fasta"
    return 0
  fi

  require_file "$sra"
  if [ ! -s "$read_fastq" ]; then
    info "Converting SRA to split-spot FASTQ for $run ($sample_id; $sample_name; $note)"
    /usr/bin/time -v fasterq-dump "$sra" \
      --threads "$THREADS" \
      --temp "$run_work/tmp" \
      --outdir "$FASTA_ROOT" \
      --split-spot \
      > "$run_work/logs/$run.fasterq.stdout.log" \
      2> "$run_work/logs/$run.fasterq.stderr.log"

    if [ -s "$FASTA_ROOT/$run.fastq" ] && [ ! -e "$raw_fastq" ]; then
      mv "$FASTA_ROOT/$run.fastq" "$raw_fastq"
    fi
    [ -s "$raw_fastq" ] || die "fasterq-dump did not produce $raw_fastq"
    pigz -p "$THREADS" "$raw_fastq"
  fi

  gzip -t "$read_fastq"
  seqkit stats -a -T "$read_fastq" \
    > "$run_work/logs/$run.split-spot.seqkit-stats.tsv"

  info "Converting FASTQ.gz to FASTA.gz for $run"
  tmp_fasta="$read_fasta.tmp.$$"
  seqkit fq2fa -j "$THREADS" "$read_fastq" | pigz -p "$THREADS" > "$tmp_fasta"
  gzip -t "$tmp_fasta"
  mv "$tmp_fasta" "$read_fasta"
  seqkit stats -a -T "$read_fasta" \
    > "$run_work/logs/$run.split-spot.fasta.seqkit-stats.tsv"

  if [ "$DROP_FASTQ_AFTER_FASTA" = "1" ]; then
    rm -f "$read_fastq"
    info "Deleted intermediate FASTQ.gz for $run because DROP_FASTQ_AFTER_FASTA=1"
  fi
}

prepare_fasta_manifest() {
  ensure_dirs
  create_default_manifest
  require_sra_tools
  require_conversion_tools
  for_manifest_rows prepare_fasta_one
}

report_exists() {
  local state="$1"
  local report_id="$2"
  "$GENTLE_BIN" --state "$state" rna-reads show-report "$report_id" >/dev/null 2>&1
}

show_report_aligned_count() {
  local show_json="$1"
  jq -r '.report.read_count_aligned // 0' "$show_json"
}

write_final_summary() {
  local run="$1"
  local sample_id="$2"
  local sample_name="$3"
  local note="$4"
  local report_id="$5"
  local show_json="$6"
  local gene_json="$7"
  local out_json="$8"
  local empty_gene_json

  empty_gene_json="$COHORT_ROOT/tmp/empty_gene_support.json"
  mkdir -p "$(dirname "$out_json")" "$COHORT_ROOT/tmp"
  if [ ! -s "$empty_gene_json" ]; then
    printf '{}\n' > "$empty_gene_json"
  fi
  if [ ! -s "$gene_json" ]; then
    gene_json="$empty_gene_json"
  fi

  jq -n \
    --arg run "$run" \
    --arg sample_id "$sample_id" \
    --arg sample_name "$sample_name" \
    --arg note "$note" \
    --arg report_id "$report_id" \
    --slurpfile show "$show_json" \
    --slurpfile gene "$gene_json" '
    def len_at($counts; $idx): (($counts // [])[$idx] // 0);
    def count_sum($counts): reduce ($counts // [])[] as $count (0; . + ($count // 0));
    def weighted_sum($counts):
      reduce range(0; (($counts // []) | length)) as $idx
        (0; . + ($idx * len_at($counts; $idx)));
    def first_observed($counts):
      reduce range(0; (($counts // []) | length)) as $idx
        (null; if . == null and len_at($counts; $idx) > 0 then $idx else . end);
    def last_observed($counts):
      reduce range(0; (($counts // []) | length)) as $idx
        (null; if len_at($counts; $idx) > 0 then $idx else . end);
    def quantile_len($counts; $fraction):
      (count_sum($counts)) as $total
      | if $total == 0 then null
        else ((($total - 1) * $fraction) | floor) as $rank
        | reduce range(0; (($counts // []) | length)) as $idx
            ({cumulative: 0, value: null};
             if .value != null then .
             else (.cumulative + len_at($counts; $idx)) as $next
             | if $next > $rank
               then {cumulative: $next, value: $idx}
               else {cumulative: $next, value: null}
               end
             end)
        | .value
        end;
    def length_summary($counts):
      (count_sum($counts)) as $total
      | {
          sample_count: $total,
          mean_bp: (if $total == 0 then null else (weighted_sum($counts) / $total) end),
          quantiles_bp: {
            q0: first_observed($counts),
            q25: quantile_len($counts; 0.25),
            q50: quantile_len($counts; 0.50),
            q75: quantile_len($counts; 0.75),
            q100: last_observed($counts)
          }
        };
    ($show[0].report) as $report
    | ($gene[0] // {}) as $gene_support
    | ($report.hits // []) as $hits
    | {
        run_accession: $run,
        sample_id: $sample_id,
        sample_name: $sample_name,
        note: $note,
        report_id: ($report.report_id // $report_id),
        read_count_total: ($report.read_count_total // 0),
        strict_seed_passed_reads: ($report.read_count_seed_passed // 0),
        retained_report_rows: ($hits | length),
        retained_aligned_rows: ($report.read_count_aligned // 0),
        retained_msa_eligible_rows: ($report.retained_count_msa_eligible // 0),
        strict_seed_passed_aligned_rows: (
          $hits
          | map(select((.passed_seed_filter == true) and (.best_mapping != null)))
          | length
        ),
        strict_seed_passed_msa_eligible_rows: (
          $hits
          | map(select((.passed_seed_filter == true) and (.msa_eligible == true)))
          | length
        ),
        read_length_summaries: {
          all_reads: length_summary($report.read_length_counts_all),
          strict_seed_passed: length_summary($report.read_length_counts_seed_passed),
          retained_aligned: length_summary($report.read_length_counts_aligned),
          full_length_exact: length_summary($report.read_length_counts_full_length_exact),
          full_length_near: length_summary($report.read_length_counts_full_length_near),
          full_length_strict: length_summary($report.read_length_counts_full_length_strict)
        },
        gene_support: {
          complete_rule: ($gene_support.complete_rule // null),
          aligned_base_count: ($gene_support.aligned_base_count // null),
          accepted_target_count: ($gene_support.accepted_target_count // null),
          accepted_target_fraction_total: (
            if (($report.read_count_total // 0) > 0) and (($gene_support.accepted_target_count // null) != null)
            then ($gene_support.accepted_target_count / $report.read_count_total)
            else null end
          ),
          accepted_target_fraction_aligned: (
            if (($gene_support.aligned_base_count // 0) > 0) and (($gene_support.accepted_target_count // null) != null)
            then ($gene_support.accepted_target_count / $gene_support.aligned_base_count)
            else null end
          ),
          fragment_count: ($gene_support.fragment_count // null),
          complete_near_count: ($gene_support.complete_count // null),
          complete_strict_count: ($gene_support.complete_strict_count // null),
          complete_exact_count: ($gene_support.complete_exact_count // null),
          accepted_target_read_lengths: ($gene_support.accepted_target_read_lengths // null),
          accepted_target_fragment_lengths: ($gene_support.accepted_target_fragment_lengths // null),
          accepted_target_query_coverage: ($gene_support.accepted_target_query_coverage // null)
        }
      }
  ' > "$out_json"
}

run_analysis_one() {
  local run="$1"
  local sample_id="$2"
  local sample_name="$3"
  local note="$4"
  local run_work state report_id seed_feature_id read_fasta show_json gene_json summary_json
  run_work="$(run_work_for "$run")"
  state="$(state_for "$run")"
  report_id="$(report_id_for "$run" "$sample_id")"
  read_fasta="$(read_fasta_for "$run")"
  show_json="$run_work/post_interpret/json/$report_id.show_report.json"
  gene_json="$run_work/post_interpret/json/$report_id.$GENE_ID.gene_support.$COMPLETE_RULE.json"
  summary_json="$run_work/reports/$report_id.final_summary.json"

  if [ -s "$summary_json" ]; then
    info "Final summary already exists for $run: $summary_json"
    return 0
  fi

  init_shared
  mkdir -p "$run_work"/{logs,reports,checkpoints,post_interpret/{json,tsv,svg,fasta}}

  if [ ! -s "$read_fasta" ]; then
    download_one "$run" "$sample_id" "$sample_name" "$note"
    prepare_fasta_one "$run" "$sample_id" "$sample_name" "$note"
  fi
  require_file "$read_fasta"
  gzip -t "$read_fasta"

  if [ ! -s "$state" ]; then
    cp "$BASE_STATE" "$state"
  fi
  seed_feature_id="$(cat "$COHORT_ROOT/manifests/seed_feature_id.txt")"
  require_file "$SEED_ARGS_FILE"

  local -a gentle seed_args
  gentle=("$GENTLE_BIN" --state "$state" --progress-stderr)
  read -r -a seed_args < "$SEED_ARGS_FILE"

  if ! report_exists "$state" "$report_id"; then
    info "Interpreting $run as $report_id"
    /usr/bin/time -v "${gentle[@]}" rna-reads interpret "$SEQ_ID" "$seed_feature_id" "$read_fasta" \
      --report-id "$report_id" \
      --report-mode seed_passed_only \
      --checkpoint-path "$run_work/checkpoints/$report_id.checkpoint.json" \
      --checkpoint-every-reads 100000 \
      --profile nanopore_cdna_v1 \
      --format fasta \
      --scope all_overlapping_any_strand \
      --origin-mode single_gene \
      "${seed_args[@]}" \
      --align-band-bp "$ALIGN_BAND_BP" \
      --align-min-identity "$ALIGN_MIN_IDENTITY" \
      --max-secondary-mappings "$MAX_SECONDARY_MAPPINGS" \
      > "$run_work/reports/$report_id.interpret.command.json" \
      2> >(tee "$run_work/logs/$report_id.interpret.stderr.log" >&2)
  else
    info "Report already exists in state for $run: $report_id"
  fi

  "${gentle[@]}" rna-reads show-report "$report_id" \
    > "$show_json" \
    2> "$run_work/logs/$report_id.show_report.before_align.stderr.log"

  if [ "$(show_report_aligned_count "$show_json")" = "0" ]; then
    info "Harvesting alignments for $run"
    "${gentle[@]}" rna-reads align-report "$report_id" \
      --selection all \
      --align-band-bp "$ALIGN_BAND_BP" \
      --align-min-identity "$ALIGN_MIN_IDENTITY" \
      --max-secondary-mappings "$MAX_SECONDARY_MAPPINGS" \
      > "$run_work/post_interpret/json/$report_id.align_report.json" \
      2> >(tee "$run_work/logs/$report_id.align_report.stderr.log" >&2)
  else
    info "Report already has aligned retained rows for $run"
  fi

  "${gentle[@]}" rna-reads show-report "$report_id" \
    > "$show_json" \
    2> "$run_work/logs/$report_id.show_report.stderr.log"

  "${gentle[@]}" rna-reads summarize-gene-support "$report_id" \
    --gene "$GENE_ID" --complete-rule "$COMPLETE_RULE" \
    --output "$gene_json" \
    > "$run_work/post_interpret/json/$report_id.$GENE_ID.gene_support.$COMPLETE_RULE.command.json" \
    2> "$run_work/logs/$report_id.$GENE_ID.gene_support.$COMPLETE_RULE.stderr.log"

  "${gentle[@]}" rna-reads export-sample-sheet \
    "$run_work/post_interpret/tsv/$report_id.sample_sheet.$COMPLETE_RULE.tsv" \
    --report-id "$report_id" \
    --gene "$GENE_ID" \
    --complete-rule "$COMPLETE_RULE" \
    > "$run_work/post_interpret/json/$report_id.sample_sheet.$COMPLETE_RULE.command.json" \
    2> "$run_work/logs/$report_id.sample_sheet.$COMPLETE_RULE.stderr.log"

  "${gentle[@]}" rna-reads export-paths-tsv "$report_id" \
    "$run_work/post_interpret/tsv/$report_id.paths.aligned.tsv" \
    --selection aligned --subset-spec "pancreas_cohort_${run}" \
    > "$run_work/post_interpret/json/$report_id.export_paths.command.json" \
    2> "$run_work/logs/$report_id.export_paths.stderr.log"

  "${gentle[@]}" rna-reads export-abundance-tsv "$report_id" \
    "$run_work/post_interpret/tsv/$report_id.abundance.aligned.tsv" \
    --selection aligned --subset-spec "pancreas_cohort_${run}" \
    > "$run_work/post_interpret/json/$report_id.export_abundance.command.json" \
    2> "$run_work/logs/$report_id.export_abundance.stderr.log"

  "${gentle[@]}" rna-reads export-target-quality "$report_id" \
    "$run_work/post_interpret/svg/$report_id.$GENE_ID.target_quality.svg" \
    --gene "$GENE_ID" --complete-rule "$COMPLETE_RULE" \
    > "$run_work/post_interpret/json/$report_id.$GENE_ID.target_quality.command.json" \
    2> "$run_work/logs/$report_id.$GENE_ID.target_quality.stderr.log"

  write_final_summary "$run" "$sample_id" "$sample_name" "$note" \
    "$report_id" "$show_json" "$gene_json" "$summary_json"
}

run_manifest() {
  ensure_dirs
  create_default_manifest
  for_manifest_rows run_analysis_one
}

one_run() {
  local run="$1"
  local sample_id="${2:-$run}"
  local sample_name="${3:-$sample_id}"
  local note="${4:-manual_one_run}"
  download_one "$run" "$sample_id" "$sample_name" "$note"
  prepare_fasta_one "$run" "$sample_id" "$sample_name" "$note"
  run_analysis_one "$run" "$sample_id" "$sample_name" "$note"
  summarize_outputs
}

adopt_existing() {
  local original_run_work="${1:-}"
  local run="${2:-}"
  local sample_id="${3:-}"
  [ -n "$original_run_work" ] || die "adopt-existing requires RUN_WORK"
  require_file "$original_run_work/run.env"

  local original_report_id original_run original_gene original_show original_gene_json
  original_report_id="$(bash -c 'source "$1"; printf "%s" "$REPORT_ID"' bash "$original_run_work/run.env")"
  original_run="$(bash -c 'source "$1"; printf "%s" "${RUN:-}"' bash "$original_run_work/run.env")"
  original_gene="$(bash -c 'source "$1"; printf "%s" "${GENE_ID:-TP73}"' bash "$original_run_work/run.env")"
  run="${run:-$original_run}"
  sample_id="${sample_id:-$run}"
  original_show="$original_run_work/post_interpret/json/$original_report_id.show_report.json"
  original_gene_json="$original_run_work/post_interpret/json/$original_report_id.$original_gene.gene_support.$COMPLETE_RULE.json"
  require_file "$original_show"

  local run_work summary_json
  run_work="$(run_work_for "$run")"
  summary_json="$run_work/reports/$original_report_id.final_summary.json"
  mkdir -p "$run_work/reports"
  write_final_summary "$run" "$sample_id" "$sample_id" "adopted_existing_run" \
    "$original_report_id" "$original_show" "$original_gene_json" "$summary_json"
  info "Adopted existing run into cohort summary: $summary_json"
}

summarize_outputs() {
  ensure_dirs
  local summaries_json="$COHORT_ROOT/reports/tp73_pancreas_cohort.final_summaries.json"
  local summary_tsv="$COHORT_ROOT/reports/tp73_pancreas_cohort.summary.tsv"
  local sample_sheet="$COHORT_ROOT/reports/tp73_pancreas_cohort.sample_sheet.$COMPLETE_RULE.tsv"
  local summary_files_file="$COHORT_ROOT/manifests/final_summary_files.txt"
  local sheets_file="$COHORT_ROOT/manifests/sample_sheet_files.txt"

  find "$COHORT_ROOT/runs" -path "*/reports/*.final_summary.json" | sort > "$summary_files_file"
  if [ ! -s "$summary_files_file" ]; then
    die "No final summaries found under $COHORT_ROOT/runs"
  fi

  jq -s '.' $(cat "$summary_files_file") > "$summaries_json"

  {
    printf 'run_accession\tsample_id\tsample_name\treport_id\tread_count_total\tstrict_seed_passed_reads\tstrict_seed_passed_aligned_rows\tstrict_seed_passed_msa_eligible_rows\tretained_report_rows\tretained_aligned_rows\taccepted_target_count\taccepted_target_fraction_total\taccepted_target_fraction_aligned\tfragment_count\tcomplete_near_count\tcomplete_strict_count\tcomplete_exact_count\tall_mean_bp\tall_q25_bp\tall_q50_bp\tall_q75_bp\tall_q100_bp\taccepted_mean_bp\taccepted_fragment_mean_bp\taccepted_query_coverage_median_fraction\n'
    jq -r '
      .[]
      | [
          .run_accession,
          .sample_id,
          .sample_name,
          .report_id,
          .read_count_total,
          .strict_seed_passed_reads,
          .strict_seed_passed_aligned_rows,
          .strict_seed_passed_msa_eligible_rows,
          .retained_report_rows,
          .retained_aligned_rows,
          .gene_support.accepted_target_count,
          .gene_support.accepted_target_fraction_total,
          .gene_support.accepted_target_fraction_aligned,
          .gene_support.fragment_count,
          .gene_support.complete_near_count,
          .gene_support.complete_strict_count,
          .gene_support.complete_exact_count,
          .read_length_summaries.all_reads.mean_bp,
          .read_length_summaries.all_reads.quantiles_bp.q25,
          .read_length_summaries.all_reads.quantiles_bp.q50,
          .read_length_summaries.all_reads.quantiles_bp.q75,
          .read_length_summaries.all_reads.quantiles_bp.q100,
          .gene_support.accepted_target_read_lengths.mean_length_bp,
          .gene_support.accepted_target_fragment_lengths.mean_length_bp,
          .gene_support.accepted_target_query_coverage.median_fraction
        ]
      | @tsv
    ' "$summaries_json"
  } > "$summary_tsv"

  find "$COHORT_ROOT/runs" -path "*/post_interpret/tsv/*.sample_sheet.$COMPLETE_RULE.tsv" \
    | sort > "$sheets_file"
  if [ -s "$sheets_file" ]; then
    {
      head -n 1 "$(head -n 1 "$sheets_file")"
      while read -r sheet; do
        tail -n +2 "$sheet"
      done < "$sheets_file"
    } > "$sample_sheet"
  fi

  find "$COHORT_ROOT/runs" -path "*/post_interpret/tsv/*.paths.aligned.tsv" \
    | sort > "$COHORT_ROOT/manifests/path_tables.txt"
  find "$COHORT_ROOT/runs" -path "*/post_interpret/tsv/*.abundance.aligned.tsv" \
    | sort > "$COHORT_ROOT/manifests/abundance_tables.txt"

  info "Wrote cohort JSON: $summaries_json"
  info "Wrote cohort TSV:  $summary_tsv"
  if [ -s "$sample_sheet" ]; then
    info "Wrote sample sheet: $sample_sheet"
  fi
}

status_one() {
  local run="$1"
  local sample_id="$2"
  local sample_name="$3"
  local note="$4"
  local run_work sra fasta summary checkpoint reads seed aligned
  run_work="$(run_work_for "$run")"
  sra="$(sra_path_for "$run")"
  fasta="$(read_fasta_for "$run")"
  summary="$(find "$run_work/reports" -maxdepth 1 -name "*.final_summary.json" 2>/dev/null | sort | head -n 1 || true)"
  checkpoint="$(find "$run_work/checkpoints" -maxdepth 1 -name "*.checkpoint.json" 2>/dev/null | sort | head -n 1 || true)"
  reads=""
  seed=""
  aligned=""
  if [ -s "$summary" ]; then
    reads="$(jq -r '.read_count_total // ""' "$summary")"
    seed="$(jq -r '.strict_seed_passed_reads // ""' "$summary")"
    aligned="$(jq -r '.retained_aligned_rows // ""' "$summary")"
  elif [ -s "$checkpoint" ]; then
    reads="$(jq -r '.reads_processed // ""' "$checkpoint")"
    seed="$(jq -r '.read_count_seed_passed // ""' "$checkpoint")"
    aligned="$(jq -r '.read_count_aligned // ""' "$checkpoint")"
  fi
  printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
    "$run" \
    "$sample_id" \
    "$([ -s "$sra" ] && echo yes || echo no)" \
    "$([ -s "$fasta" ] && echo yes || echo no)" \
    "$([ -s "$summary" ] && echo complete || echo pending)" \
    "$reads" \
    "$seed" \
    "$aligned" \
    "$note"
}

status_manifest() {
  ensure_dirs
  create_default_manifest
  printf 'run_accession\tsample_id\tsra_present\tfasta_present\tanalysis_status\treads_or_processed\tseed_passed\taligned_or_retained_aligned\tnote\n'
  for_manifest_rows status_one
}

all_steps() {
  init_shared
  run_manifest
  summarize_outputs
}

main() {
  local command="${1:-}"
  shift || true
  case "$command" in
    init)
      init_shared
      ;;
    download)
      download_manifest
      ;;
    prepare-fasta|prepare)
      prepare_fasta_manifest
      ;;
    run)
      run_manifest
      ;;
    one)
      [ $# -ge 1 ] || die "one requires RUN [SAMPLE_ID]"
      one_run "$@"
      ;;
    adopt-existing)
      adopt_existing "$@"
      ;;
    summarize)
      summarize_outputs
      ;;
    status)
      status_manifest
      ;;
    all)
      all_steps
      ;;
    -h|--help|help|"")
      usage
      ;;
    *)
      usage >&2
      die "Unknown command: $command"
      ;;
  esac
}

main "$@"
