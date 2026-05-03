#!/usr/bin/env bash
# Run and harvest the strict TP73 pancreatic Nanopore cDNA RNA-mapping proof.
#
# This script is intentionally conservative:
# - it defines every path it needs,
# - it writes a run.env file for recovery/status commands,
# - it refuses to start while another RNA-read mapping process is active unless
#   ALLOW_PARALLEL_RNA_MAPPING=1 is set,
# - it runs the primary mapping in the foreground so a screen/tmux session owns
#   exactly one job.

set -uo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
DEFAULT_REPO_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"

RUN="${RUN:-SRR32957126}"
GENE_ID="${GENE_ID:-TP73}"
SEQ_ID="${SEQ_ID:-tp73_ncbi}"
THREADS="${THREADS:-6}"
GENTLE_REPO="${GENTLE_REPO:-${GENTLE_RS:-$DEFAULT_REPO_ROOT}}"
GENTLE_BIN="${GENTLE_BIN:-$GENTLE_REPO/target/release/gentle_cli}"
WORK_ROOT="${WORK_ROOT:-${WORK:-/home/clawbio/work/tp73_pancreas_benchmark}}"
PAN_ROOT="${PAN_ROOT:-/home/clawbio/data/SRA/Pancreas_Carcinoma_Nanopore_GSE293661}"
READ_FASTA="${READ_FASTA:-$WORK_ROOT/fastq/$RUN.split-spot.fasta.gz}"

usage() {
  cat <<'EOF'
Usage:
  scripts/tp73_pancreas_rna_mapping.sh strict-start
  scripts/tp73_pancreas_rna_mapping.sh status [RUN_WORK]
  scripts/tp73_pancreas_rna_mapping.sh harvest [RUN_WORK]
  scripts/tp73_pancreas_rna_mapping.sh stop [RUN_WORK]

Environment overrides:
  RUN=SRR32957126
  GENE_ID=TP73
  SEQ_ID=tp73_ncbi
  THREADS=6
  GENTLE_REPO=/home/clawbio/GENtle
  GENTLE_BIN=/home/clawbio/GENtle/target/release/gentle_cli
  WORK_ROOT=/home/clawbio/work/tp73_pancreas_benchmark
  PAN_ROOT=/home/clawbio/data/SRA/Pancreas_Carcinoma_Nanopore_GSE293661
  READ_FASTA=/home/clawbio/work/tp73_pancreas_benchmark/fastq/SRR32957126.split-spot.fasta.gz

Safety:
  strict-start refuses to run when any gentle_cli rna-reads interpret or
  align-report process is active. Set ALLOW_PARALLEL_RNA_MAPPING=1 to override.
EOF
}

die() {
  echo "ERROR: $*" >&2
  exit 1
}

info() {
  echo "INFO: $*" >&2
}

require_file() {
  local path="$1"
  [ -s "$path" ] || die "Missing required file: $path"
}

require_executable() {
  local tool="$1"
  if [[ "$tool" == */* ]]; then
    [ -x "$tool" ] || die "Missing executable: $tool"
  else
    command -v "$tool" >/dev/null 2>&1 || die "Missing executable in PATH: $tool"
  fi
}

latest_strict_run_work() {
  local latest
  latest="$(find "$WORK_ROOT" -maxdepth 1 -type d -name "strict_${RUN}_*" 2>/dev/null | sort | tail -n 1)"
  [ -n "$latest" ] || die "No strict_${RUN}_* run directory found under $WORK_ROOT"
  printf '%s\n' "$latest"
}

load_run_env() {
  local run_work="${1:-}"
  if [ -z "$run_work" ]; then
    run_work="$(latest_strict_run_work)"
  fi
  [ -d "$run_work" ] || die "Run directory does not exist: $run_work"
  [ -s "$run_work/run.env" ] || die "Missing run.env in: $run_work"
  # shellcheck disable=SC1090
  source "$run_work/run.env"
}

active_mapping_processes() {
  ps -axo pid,ppid,etime,pcpu,pmem,command \
    | grep -E 'gentle_cli .*rna-reads (interpret|align-report)' \
    | grep -v grep || true
}

refuse_if_mapping_active() {
  local active
  active="$(active_mapping_processes)"
  if [ -n "$active" ] && [ "${ALLOW_PARALLEL_RNA_MAPPING:-0}" != "1" ]; then
    echo "$active" >&2
    die "RNA mapping is already active. Stop it first or set ALLOW_PARALLEL_RNA_MAPPING=1."
  fi
}

refuse_existing_path() {
  local path="$1"
  if [ -e "$path" ]; then
    die "Refusing to overwrite existing path: $path"
  fi
}

refuse_existing_strict_outputs() {
  refuse_existing_path "$STATE"
  refuse_existing_path "$RUN_WORK/run.env"
  refuse_existing_path "$RUN_WORK/manifests/load_tp73_ncbi.workflow.json"
  refuse_existing_path "$RUN_WORK/reports/load_tp73_ncbi.command.json"
  refuse_existing_path "$RUN_WORK/logs/load_tp73_ncbi.stderr.log"
  refuse_existing_path "$RUN_WORK/logs/query_${GENE_ID}_gene.stderr.log"
  refuse_existing_path "$RUN_WORK/reports/$REPORT_ID.preflight.json"
  refuse_existing_path "$RUN_WORK/logs/$REPORT_ID.preflight.stderr.log"
  refuse_existing_path "$RUN_WORK/reports/$REPORT_ID.preflight.summary.json"
  refuse_existing_path "$RUN_WORK/checkpoints/$REPORT_ID.checkpoint.json"
  refuse_existing_path "$RUN_WORK/reports/$REPORT_ID.interpret.command.json"
  refuse_existing_path "$RUN_WORK/logs/$REPORT_ID.interpret.stderr.log"
}

refuse_existing_harvest_outputs() {
  refuse_existing_path "$RUN_WORK/reports/list_reports.before_harvest.json"
  refuse_existing_path "$RUN_WORK/logs/list_reports.before_harvest.stderr.log"
  refuse_existing_path "$RUN_WORK/post_interpret/json/$REPORT_ID.align_report.json"
  refuse_existing_path "$RUN_WORK/logs/$REPORT_ID.align_report.stderr.log"
  refuse_existing_path "$RUN_WORK/post_interpret/json/$REPORT_ID.show_report.json"
  refuse_existing_path "$RUN_WORK/logs/$REPORT_ID.show_report.stderr.log"
  refuse_existing_path "$RUN_WORK/post_interpret/json/$REPORT_ID.aligned.top200.json"
  refuse_existing_path "$RUN_WORK/logs/$REPORT_ID.aligned.top200.stderr.log"
  refuse_existing_path "$RUN_WORK/post_interpret/json/$REPORT_ID.$GENE_ID.gene_support.near.json"
  refuse_existing_path "$RUN_WORK/post_interpret/json/$REPORT_ID.$GENE_ID.gene_support.near.command.json"
  refuse_existing_path "$RUN_WORK/logs/$REPORT_ID.$GENE_ID.gene_support.near.stderr.log"
  refuse_existing_path "$RUN_WORK/post_interpret/tsv/$REPORT_ID.alignments.aligned.tsv"
  refuse_existing_path "$RUN_WORK/post_interpret/json/$REPORT_ID.export_alignments.command.json"
  refuse_existing_path "$RUN_WORK/logs/$REPORT_ID.export_alignments.stderr.log"
  refuse_existing_path "$RUN_WORK/post_interpret/tsv/$REPORT_ID.paths.aligned.tsv"
  refuse_existing_path "$RUN_WORK/post_interpret/json/$REPORT_ID.export_paths.command.json"
  refuse_existing_path "$RUN_WORK/logs/$REPORT_ID.export_paths.stderr.log"
  refuse_existing_path "$RUN_WORK/post_interpret/tsv/$REPORT_ID.abundance.aligned.tsv"
  refuse_existing_path "$RUN_WORK/post_interpret/json/$REPORT_ID.export_abundance.command.json"
  refuse_existing_path "$RUN_WORK/logs/$REPORT_ID.export_abundance.stderr.log"
  refuse_existing_path "$RUN_WORK/post_interpret/svg/$REPORT_ID.$GENE_ID.target_quality.svg"
  refuse_existing_path "$RUN_WORK/post_interpret/json/$REPORT_ID.$GENE_ID.target_quality.command.json"
  refuse_existing_path "$RUN_WORK/logs/$REPORT_ID.$GENE_ID.target_quality.stderr.log"
  refuse_existing_path "$RUN_WORK/reports/$REPORT_ID.final_summary.json"
  refuse_existing_path "$RUN_WORK/reports/$REPORT_ID.evidence_bundle.md"
}

write_run_env() {
  local run_work="$1"
  cat > "$run_work/run.env" <<EOF
export RUN="$RUN"
export GENE_ID="$GENE_ID"
export SEQ_ID="$SEQ_ID"
export THREADS="$THREADS"
export GENTLE_REPO="$GENTLE_REPO"
export GENTLE_BIN="$GENTLE_BIN"
export WORK_ROOT="$WORK_ROOT"
export PAN_ROOT="$PAN_ROOT"
export READ_FASTA="$READ_FASTA"
export RUN_WORK="$RUN_WORK"
export STATE="$STATE"
export REPORT_ID="$REPORT_ID"
export SEED_FEATURE_ID="${SEED_FEATURE_ID:-}"
export SEED_FRAGMENT="${SEED_FRAGMENT:-}"
EOF
}

check_common_inputs() {
  cd "$GENTLE_REPO" || die "Cannot cd to GENTLE_REPO: $GENTLE_REPO"

  require_executable "$GENTLE_BIN"
  require_executable jq
  require_executable gzip
  require_executable seqkit
  require_executable /usr/bin/time

  require_file "$READ_FASTA"
  require_file "$PAN_ROOT/$RUN/$RUN.sra"
  require_file test_files/tp73.ncbi.gb
  require_file test_files/fixtures/mapping/ensembl_human_tp73_all.fasta
  require_file test_files/fixtures/mapping/ensembl_chimp_tp73_all.fasta
  require_file test_files/fixtures/mapping/ensembl_human_tp53_all.fasta
  require_file test_files/fixtures/mapping/ensembl_human_tp63_all.fasta
  require_file data/resources/tp73_dn_ena_transcripts.fasta
  require_file data/resources/tp73_delta_ex2_3_refseq_derived_transcripts.fasta
}

strict_start() {
  refuse_if_mapping_active
  check_common_inputs

  local stamp
  stamp="$(date -u +%Y%m%dT%H%M%SZ)"
  RUN_WORK="${RUN_WORK:-$WORK_ROOT/strict_${RUN}_${stamp}}"
  STATE="${STATE:-$RUN_WORK/$RUN.strict.gentle.json}"
  REPORT_ID="${REPORT_ID:-tp73_strict_${RUN}_${stamp}}"

  refuse_existing_path "$RUN_WORK"
  mkdir -p "$RUN_WORK"/{logs,reports,manifests,checkpoints,post_interpret/{json,tsv,svg,fasta}}
  refuse_existing_strict_outputs
  write_run_env "$RUN_WORK"

  info "Run directory: $RUN_WORK"
  info "State: $STATE"
  info "Report ID: $REPORT_ID"

  gzip -t "$READ_FASTA" || die "gzip integrity check failed for $READ_FASTA"
  seqkit stats -a -T "$READ_FASTA" | tee "$RUN_WORK/reports/$RUN.fasta.seqkit-stats.tsv"
  git rev-parse HEAD | tee "$RUN_WORK/reports/gentle_git_commit.txt" >/dev/null || true

  local -a gentle
  gentle=("$GENTLE_BIN" --state "$STATE" --progress-stderr)

  cat > "$RUN_WORK/manifests/load_tp73_ncbi.workflow.json" <<'JSON'
{
  "run_id": "load_tp73_ncbi_for_strict_tp73_pancreas",
  "ops": [
    { "LoadFile": { "path": "test_files/tp73.ncbi.gb", "as_id": "tp73_ncbi" } }
  ]
}
JSON

  "${gentle[@]}" workflow @"$RUN_WORK/manifests/load_tp73_ncbi.workflow.json" \
    > "$RUN_WORK/reports/load_tp73_ncbi.command.json" \
    2> "$RUN_WORK/logs/load_tp73_ncbi.stderr.log" \
    || die "Could not load TP73 GenBank fixture; see $RUN_WORK/logs/load_tp73_ncbi.stderr.log"

  SEED_FEATURE_ID="$(
    "${gentle[@]}" features query "$SEQ_ID" \
      --kind gene --label "$GENE_ID" --limit 1 \
      2> "$RUN_WORK/logs/query_${GENE_ID}_gene.stderr.log" \
      | jq -r '.rows[0].feature_id'
  )"
  [ -n "$SEED_FEATURE_ID" ] && [ "$SEED_FEATURE_ID" != "null" ] \
    || die "Could not resolve seed feature ID; see $RUN_WORK/logs/query_${GENE_ID}_gene.stderr.log"
  write_run_env "$RUN_WORK"
  info "Seed feature ID: $SEED_FEATURE_ID"

  "${gentle[@]}" rna-reads preflight-isoforms "$SEQ_ID" "$SEED_FEATURE_ID" \
    --scope all_overlapping_any_strand \
    --must-pass-transcript-fasta test_files/fixtures/mapping/ensembl_human_tp73_all.fasta \
    --must-pass-transcript-fasta test_files/fixtures/mapping/ensembl_chimp_tp73_all.fasta \
    --must-pass-transcript-fasta data/resources/tp73_dn_ena_transcripts.fasta \
    --must-pass-transcript-fasta data/resources/tp73_delta_ex2_3_refseq_derived_transcripts.fasta \
    --control-transcript-fasta test_files/fixtures/mapping/ensembl_human_tp53_all.fasta \
    --control-transcript-fasta test_files/fixtures/mapping/ensembl_human_tp63_all.fasta \
    --optimize-parameters \
    --max-control-match-probability 0.0 \
    > "$RUN_WORK/reports/$REPORT_ID.preflight.json" \
    2> "$RUN_WORK/logs/$REPORT_ID.preflight.stderr.log" \
    || die "Preflight failed; see $RUN_WORK/logs/$REPORT_ID.preflight.stderr.log"

  SEED_FRAGMENT="$(
    jq -r '.threshold_recommendation.seed_filter_cli_fragment' \
      "$RUN_WORK/reports/$REPORT_ID.preflight.json"
  )"
  [ -n "$SEED_FRAGMENT" ] && [ "$SEED_FRAGMENT" != "null" ] \
    || die "Preflight did not produce seed_filter_cli_fragment"
  write_run_env "$RUN_WORK"

  jq '{
    target: [.target_passed_transcript_count, .target_transcript_count],
    positives: [.positive_control_passed_transcript_count, .positive_control_transcript_count],
    controls: [.control_summaries[] | {
      control_id,
      transcript_count,
      passed_transcript_count,
      weighted_pass_probability,
      best_transcript_id
    }],
    threshold_recommendation
  }' "$RUN_WORK/reports/$REPORT_ID.preflight.json" \
    | tee "$RUN_WORK/reports/$REPORT_ID.preflight.summary.json"

  local -a seed_args
  read -r -a seed_args <<< "$SEED_FRAGMENT"

  info "Starting strict mapping in the foreground. Leave this screen running."
  /usr/bin/time -v "${gentle[@]}" rna-reads interpret "$SEQ_ID" "$SEED_FEATURE_ID" "$READ_FASTA" \
    --report-id "$REPORT_ID" \
    --report-mode seed_passed_only \
    --checkpoint-path "$RUN_WORK/checkpoints/$REPORT_ID.checkpoint.json" \
    --checkpoint-every-reads 100000 \
    --profile nanopore_cdna_v1 \
    --format fasta \
    --scope all_overlapping_any_strand \
    --origin-mode single_gene \
    "${seed_args[@]}" \
    --align-band-bp 24 \
    --align-min-identity 0.55 \
    --max-secondary-mappings 3 \
    > "$RUN_WORK/reports/$REPORT_ID.interpret.command.json" \
    2> >(tee "$RUN_WORK/logs/$REPORT_ID.interpret.stderr.log" >&2)
}

status_run() {
  load_run_env "${1:-}"
  echo "RUN_WORK=$RUN_WORK"
  echo "STATE=$STATE"
  echo "REPORT_ID=$REPORT_ID"
  echo
  echo "Active mapping processes:"
  active_mapping_processes
  echo
  if [ -s "$RUN_WORK/checkpoints/$REPORT_ID.checkpoint.json" ]; then
    jq '{
      reads_processed,
      read_count_seed_passed,
      read_count_aligned,
      input_bytes_processed,
      input_bytes_total,
      percent_bytes: (100 * .input_bytes_processed / .input_bytes_total)
    }' "$RUN_WORK/checkpoints/$REPORT_ID.checkpoint.json"
  else
    echo "No checkpoint yet: $RUN_WORK/checkpoints/$REPORT_ID.checkpoint.json"
  fi
  echo
  if [ -s "$STATE" ]; then
    "$GENTLE_BIN" --state "$STATE" rna-reads list-reports "$SEQ_ID" || true
  fi
}

stop_run() {
  load_run_env "${1:-}"
  local pids
  pids="$(
    ps -axo pid,command \
      | grep -F "$REPORT_ID" \
      | grep -E 'gentle_cli|/usr/bin/time' \
      | grep -v grep \
      | awk '{print $1}'
  )"
  if [ -z "$pids" ]; then
    info "No active processes found for $REPORT_ID"
    return 0
  fi
  echo "$pids" | xargs kill -INT
  info "Sent SIGINT to: $pids"
}

harvest_run() {
  load_run_env "${1:-}"
  cd "$GENTLE_REPO" || die "Cannot cd to GENTLE_REPO: $GENTLE_REPO"

  local active
  active="$(active_mapping_processes | grep -F "$REPORT_ID" || true)"
  [ -z "$active" ] || die "Mapping still active for $REPORT_ID; wait or stop it before harvest."

  local -a gentle
  gentle=("$GENTLE_BIN" --state "$STATE" --progress-stderr)
  refuse_existing_harvest_outputs

  "${gentle[@]}" rna-reads list-reports "$SEQ_ID" \
    > "$RUN_WORK/reports/list_reports.before_harvest.json" \
    2> "$RUN_WORK/logs/list_reports.before_harvest.stderr.log" \
    || die "Could not list reports from $STATE"

  jq -e --arg report_id "$REPORT_ID" '.reports[] | select(.report_id == $report_id)' \
    "$RUN_WORK/reports/list_reports.before_harvest.json" >/dev/null \
    || die "Report $REPORT_ID is not stored in $STATE yet. Inspect checkpoint/log before harvest."

  "${gentle[@]}" rna-reads align-report "$REPORT_ID" \
    --selection all \
    --align-band-bp 24 \
    --align-min-identity 0.55 \
    --max-secondary-mappings 3 \
    > "$RUN_WORK/post_interpret/json/$REPORT_ID.align_report.json" \
    2> "$RUN_WORK/logs/$REPORT_ID.align_report.stderr.log"

  "${gentle[@]}" rna-reads show-report "$REPORT_ID" \
    > "$RUN_WORK/post_interpret/json/$REPORT_ID.show_report.json" \
    2> "$RUN_WORK/logs/$REPORT_ID.show_report.stderr.log"

  "${gentle[@]}" rna-reads inspect-alignments "$REPORT_ID" \
    --selection aligned --limit 200 --sort score \
    > "$RUN_WORK/post_interpret/json/$REPORT_ID.aligned.top200.json" \
    2> "$RUN_WORK/logs/$REPORT_ID.aligned.top200.stderr.log"

  "${gentle[@]}" rna-reads summarize-gene-support "$REPORT_ID" \
    --gene "$GENE_ID" --complete-rule near \
    --output "$RUN_WORK/post_interpret/json/$REPORT_ID.$GENE_ID.gene_support.near.json" \
    > "$RUN_WORK/post_interpret/json/$REPORT_ID.$GENE_ID.gene_support.near.command.json" \
    2> "$RUN_WORK/logs/$REPORT_ID.$GENE_ID.gene_support.near.stderr.log"

  "${gentle[@]}" rna-reads export-alignments-tsv "$REPORT_ID" \
    "$RUN_WORK/post_interpret/tsv/$REPORT_ID.alignments.aligned.tsv" \
    --selection aligned --subset-spec "strict_tp73_pancreas_${RUN}" \
    > "$RUN_WORK/post_interpret/json/$REPORT_ID.export_alignments.command.json" \
    2> "$RUN_WORK/logs/$REPORT_ID.export_alignments.stderr.log"

  "${gentle[@]}" rna-reads export-paths-tsv "$REPORT_ID" \
    "$RUN_WORK/post_interpret/tsv/$REPORT_ID.paths.aligned.tsv" \
    --selection aligned --subset-spec "strict_tp73_pancreas_${RUN}" \
    > "$RUN_WORK/post_interpret/json/$REPORT_ID.export_paths.command.json" \
    2> "$RUN_WORK/logs/$REPORT_ID.export_paths.stderr.log"

  "${gentle[@]}" rna-reads export-abundance-tsv "$REPORT_ID" \
    "$RUN_WORK/post_interpret/tsv/$REPORT_ID.abundance.aligned.tsv" \
    --selection aligned --subset-spec "strict_tp73_pancreas_${RUN}" \
    > "$RUN_WORK/post_interpret/json/$REPORT_ID.export_abundance.command.json" \
    2> "$RUN_WORK/logs/$REPORT_ID.export_abundance.stderr.log"

  "${gentle[@]}" rna-reads export-target-quality "$REPORT_ID" \
    "$RUN_WORK/post_interpret/svg/$REPORT_ID.$GENE_ID.target_quality.svg" \
    --gene "$GENE_ID" --complete-rule near \
    > "$RUN_WORK/post_interpret/json/$REPORT_ID.$GENE_ID.target_quality.command.json" \
    2> "$RUN_WORK/logs/$REPORT_ID.$GENE_ID.target_quality.stderr.log"

  jq '
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
    .report as $report
    | ($report.hits // []) as $hits
    | {
        report_id: $report.report_id,
        total_reads: $report.read_count_total,
        strict_seed_passed_reads: $report.read_count_seed_passed,
        retained_report_rows: ($hits | length),
        retained_aligned_rows: $report.read_count_aligned,
        retained_msa_eligible_rows: $report.retained_count_msa_eligible,
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
        note: "retained_aligned_rows is the alignment count over retained report rows after harvest --selection all; strict_seed_passed_* are the stringent TP73 seed-pass subset."
      }
  ' "$RUN_WORK/post_interpret/json/$REPORT_ID.show_report.json" \
    | tee "$RUN_WORK/reports/$REPORT_ID.final_summary.json"

  {
    echo "# TP73 pancreas strict RNA mapping evidence bundle"
    echo
    echo "- run: $RUN"
    echo "- state: $STATE"
    echo "- report: $REPORT_ID"
    echo "- read FASTA: $READ_FASTA"
    echo "- preflight: $RUN_WORK/reports/$REPORT_ID.preflight.json"
    echo "- final summary: $RUN_WORK/reports/$REPORT_ID.final_summary.json"
    echo "- top alignments: $RUN_WORK/post_interpret/json/$REPORT_ID.aligned.top200.json"
    echo "- alignments TSV: $RUN_WORK/post_interpret/tsv/$REPORT_ID.alignments.aligned.tsv"
    echo "- target quality SVG: $RUN_WORK/post_interpret/svg/$REPORT_ID.$GENE_ID.target_quality.svg"
  } > "$RUN_WORK/reports/$REPORT_ID.evidence_bundle.md"
}

main() {
  local command="${1:-}"
  shift || true
  case "$command" in
    strict-start)
      strict_start "$@"
      ;;
    status)
      status_run "${1:-}"
      ;;
    harvest)
      harvest_run "${1:-}"
      ;;
    stop)
      stop_run "${1:-}"
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
