#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
DEFAULT_OUT_DIR="${ROOT_DIR}/data/resources/jaspar_interest_catalog"
DEFAULT_RANDOM_LENGTH=10000
DEFAULT_RANDOM_SEED=123

usage() {
  cat <<'EOF'
benchmark_jaspar_interest_catalog.sh

Cluster-friendly helper for incrementally filling a stored JASPAR benchmark
catalog using GENtle's own internal motif machinery.

This script intentionally benchmarks one motif per task and stores one JSON
report per entry. A later merge step rebuilds a partial
`gentle.jaspar_registry_benchmark.v1` catalog from the stored entry reports.

Subcommands:

  run
    Benchmark one motif and write one stored entry report.

    Options:
      --motif MOTIF              Benchmark this motif id directly.
      --motifs-file FILE         Read motif ids from FILE (blank lines and '#'
                                 comments are ignored). When used with a
                                 cluster array, one task benchmarks one line.
      --task-index N             1-based index into --motifs-file. If omitted,
                                 the script tries one of:
                                   SLURM_ARRAY_TASK_ID
                                   SGE_TASK_ID
                                   PBS_ARRAY_INDEX
                                   LSB_JOBINDEX
      --out-dir DIR              Output directory
                                 (default: data/resources/jaspar_interest_catalog)
      --random-length N          Background DNA length (default: 10000)
      --seed N                   Deterministic random seed (default: 123)
      --gentle-cli-bin PATH      Use this compiled gentle_cli binary instead of
                                 `cargo run --quiet --bin gentle_cli --`.
      --force                    Recompute even if the stored JSON already exists.

  merge
    Merge stored per-entry reports into one partial registry benchmark snapshot.

    Options:
      --out-dir DIR              Output directory
      --output FILE              Merged benchmark JSON path
                                 (default: OUT_DIR/jaspar.registry_benchmark.partial.json)

Examples:

  # One local motif
  scripts/benchmark_jaspar_interest_catalog.sh run --motif SP1

  # One Slurm array task over a prepared motif list
  #SBATCH --array=1-25
  scripts/benchmark_jaspar_interest_catalog.sh run \
    --motifs-file tp_family_and_friends.txt \
    --out-dir /scratch/$USER/jaspar_interest

  # Rebuild the merged catalog after tasks complete
  scripts/benchmark_jaspar_interest_catalog.sh merge \
    --out-dir /scratch/$USER/jaspar_interest

Notes:

  - The default 10000 bp background is intentional. It is large enough for the
    low tail probabilities we care about when multiple-testing interpretation
    matters.
  - The primary storage layout is one JSON per motif under:
      OUT_DIR/entries/<motif>.presentation.json
    This avoids cluster write races on one shared catalog file.
  - If you already have a compiled binary, prefer:
      --gentle-cli-bin target/release/gentle_cli
    to avoid one Cargo invocation per task.
EOF
}

die() {
  echo "error: $*" >&2
  exit 1
}

sanitize_filename() {
  local raw="$1"
  raw="${raw//\//_}"
  raw="${raw//:/_}"
  raw="${raw// /_}"
  echo "${raw}"
}

run_gentle_cli() {
  local gentle_cli_bin="$1"
  shift
  if [[ -n "${gentle_cli_bin}" ]]; then
    "${gentle_cli_bin}" "$@"
  else
    (
      cd "${ROOT_DIR}"
      cargo run --quiet --bin gentle_cli -- "$@"
    )
  fi
}

load_motif_from_file() {
  local motifs_file="$1"
  local task_index="$2"
  [[ -f "${motifs_file}" ]] || die "motif file not found: ${motifs_file}"
  [[ -n "${task_index}" ]] || die "--motifs-file requires --task-index or a cluster array task id"
  [[ "${task_index}" =~ ^[0-9]+$ ]] || die "task index must be a positive integer"

  mapfile -t MOTIF_LINES < <(grep -v '^[[:space:]]*#' "${motifs_file}" | sed '/^[[:space:]]*$/d')
  local zero_based=$((task_index - 1))
  (( zero_based >= 0 )) || die "task index must be 1-based"
  (( zero_based < ${#MOTIF_LINES[@]} )) || die "task index ${task_index} exceeds motif list size ${#MOTIF_LINES[@]}"
  printf '%s\n' "${MOTIF_LINES[${zero_based}]}"
}

cmd_run() {
  local motif=""
  local motifs_file=""
  local out_dir="${DEFAULT_OUT_DIR}"
  local random_length="${DEFAULT_RANDOM_LENGTH}"
  local random_seed="${DEFAULT_RANDOM_SEED}"
  local task_index="${SLURM_ARRAY_TASK_ID:-${SGE_TASK_ID:-${PBS_ARRAY_INDEX:-${LSB_JOBINDEX:-}}}}"
  local gentle_cli_bin="${GENTLE_CLI_BIN:-}"
  local force=0

  while [[ $# -gt 0 ]]; do
    case "$1" in
      --motif)
        [[ $# -ge 2 ]] || die "--motif requires a value"
        motif="$2"
        shift 2
        ;;
      --motifs-file)
        [[ $# -ge 2 ]] || die "--motifs-file requires a value"
        motifs_file="$2"
        shift 2
        ;;
      --task-index)
        [[ $# -ge 2 ]] || die "--task-index requires a value"
        task_index="$2"
        shift 2
        ;;
      --out-dir)
        [[ $# -ge 2 ]] || die "--out-dir requires a value"
        out_dir="$2"
        shift 2
        ;;
      --random-length)
        [[ $# -ge 2 ]] || die "--random-length requires a value"
        random_length="$2"
        shift 2
        ;;
      --seed)
        [[ $# -ge 2 ]] || die "--seed requires a value"
        random_seed="$2"
        shift 2
        ;;
      --gentle-cli-bin)
        [[ $# -ge 2 ]] || die "--gentle-cli-bin requires a value"
        gentle_cli_bin="$2"
        shift 2
        ;;
      --force)
        force=1
        shift
        ;;
      -h|--help)
        usage
        exit 0
        ;;
      *)
        die "unknown run option: $1"
        ;;
    esac
  done

  if [[ -n "${motif}" && -n "${motifs_file}" ]]; then
    die "choose either --motif or --motifs-file, not both"
  fi
  if [[ -z "${motif}" && -z "${motifs_file}" ]]; then
    die "run requires either --motif or --motifs-file"
  fi
  if [[ -n "${motifs_file}" ]]; then
    motif="$(load_motif_from_file "${motifs_file}" "${task_index}")"
  fi

  mkdir -p "${out_dir}/entries"
  local safe_motif
  safe_motif="$(sanitize_filename "${motif}")"
  local output_path="${out_dir}/entries/${safe_motif}.presentation.json"

  if [[ -f "${output_path}" && "${force}" -ne 1 ]]; then
    echo "skip: ${motif} already stored at ${output_path}" >&2
    return 0
  fi

  echo "benchmark: ${motif} -> ${output_path}" >&2
  run_gentle_cli \
    "${gentle_cli_bin}" \
    resources summarize-jaspar \
    --motif "${motif}" \
    --random-length "${random_length}" \
    --seed "${random_seed}" \
    --output "${output_path}"
}

cmd_merge() {
  local out_dir="${DEFAULT_OUT_DIR}"
  local output_path=""

  while [[ $# -gt 0 ]]; do
    case "$1" in
      --out-dir)
        [[ $# -ge 2 ]] || die "--out-dir requires a value"
        out_dir="$2"
        shift 2
        ;;
      --output)
        [[ $# -ge 2 ]] || die "--output requires a value"
        output_path="$2"
        shift 2
        ;;
      -h|--help)
        usage
        exit 0
        ;;
      *)
        die "unknown merge option: $1"
        ;;
    esac
  done

  [[ -d "${out_dir}/entries" ]] || die "missing entries directory: ${out_dir}/entries"
  if [[ -z "${output_path}" ]]; then
    output_path="${out_dir}/jaspar.registry_benchmark.partial.json"
  fi
  mkdir -p "$(dirname "${output_path}")"

  python3 - "${out_dir}/entries" "${output_path}" <<'PY'
import glob
import json
import os
import sys
import time

entries_dir = sys.argv[1]
output_path = sys.argv[2]
paths = sorted(glob.glob(os.path.join(entries_dir, "*.json")))
if not paths:
    raise SystemExit(f"error: no stored JASPAR entry reports found under {entries_dir}")

rows = []
registry_counts = set()
random_lengths = set()
random_seeds = set()
background_models = set()

for path in paths:
    with open(path, "r", encoding="utf-8") as handle:
        report = json.load(handle)
    if report.get("schema") != "gentle.jaspar_entry_presentation.v1":
        raise SystemExit(f"error: {path} is not a gentle.jaspar_entry_presentation.v1 report")
    report_rows = report.get("rows") or []
    if len(report_rows) != 1:
        raise SystemExit(f"error: {path} expected exactly one JASPAR entry row, found {len(report_rows)}")
    rows.append(report_rows[0])
    registry_counts.add(report.get("registry_entry_count"))
    random_lengths.add(report.get("random_sequence_length_bp"))
    random_seeds.add(report.get("random_seed"))
    background_models.add(report.get("background_model"))

if len(registry_counts) != 1 or len(random_lengths) != 1 or len(random_seeds) != 1 or len(background_models) != 1:
    raise SystemExit(
        "error: stored entry reports disagree on registry size, random length, seed, or background model"
    )

rows.sort(key=lambda row: row.get("motif_id", ""))

def average(values):
    materialized = list(values)
    return sum(materialized) / len(materialized) if materialized else 0.0

def median(values):
    if not values:
        return 0.0
    ordered = sorted(values)
    mid = len(ordered) // 2
    if len(ordered) % 2:
        return ordered[mid]
    return (ordered[mid - 1] + ordered[mid]) / 2.0

def top_rows(rows, value_fn):
    ordered = sorted(
        rows,
        key=lambda row: (value_fn(row), row.get("motif_id", "")),
        reverse=True,
    )
    return [
        {
            "motif_id": row.get("motif_id", ""),
            "motif_name": row.get("motif_name"),
            "value": value_fn(row),
        }
        for row in ordered[:10]
    ]

def summary(score_kind, label, dist_key, max_key):
    distributions = [row[dist_key] for row in rows]
    return {
        "score_kind": score_kind,
        "label": label,
        "motif_count": len(rows),
        "global_min_score": min(dist["min_score"] for dist in distributions),
        "global_max_score": max(dist["max_score"] for dist in distributions),
        "mean_of_mean_scores": average(dist["mean_score"] for dist in distributions),
        "mean_of_stddev_scores": average(dist["stddev_score"] for dist in distributions),
        "median_of_p50_scores": median(dist["p50_score"] for dist in distributions),
        "mean_positive_fraction": average(dist["positive_fraction"] for dist in distributions),
        "top_max_score_rows": top_rows(rows, lambda row: row[max_key]),
        "top_positive_fraction_rows": top_rows(rows, lambda row: row[dist_key]["positive_fraction"]),
    }

merged = {
    "schema": "gentle.jaspar_registry_benchmark.v1",
    "generated_at_unix_ms": int(time.time() * 1000),
    "op_id": None,
    "run_id": None,
    "registry_entry_count": next(iter(registry_counts)),
    "benchmarked_entry_count": len(rows),
    "random_sequence_length_bp": next(iter(random_lengths)),
    "random_seed": next(iter(random_seeds)),
    "background_model": next(iter(background_models)),
    "score_family_summaries": [
        summary("llr_bits", "LLR bits", "llr_bits_distribution", "maximizing_llr_bits"),
        summary(
            "true_log_odds_bits",
            "True log-odds bits",
            "true_log_odds_bits_distribution",
            "maximizing_true_log_odds_bits",
        ),
    ],
    "rows": rows,
}

with open(output_path, "w", encoding="utf-8") as handle:
    json.dump(merged, handle, indent=2, sort_keys=False)
    handle.write("\n")

print(output_path)
PY
}

main() {
  local subcommand="${1:-}"
  if [[ -z "${subcommand}" ]]; then
    usage
    exit 1
  fi
  shift || true
  case "${subcommand}" in
    run)
      cmd_run "$@"
      ;;
    merge)
      cmd_merge "$@"
      ;;
    -h|--help|help)
      usage
      ;;
    *)
      die "unknown subcommand: ${subcommand}"
      ;;
  esac
}

main "$@"
