#!/usr/bin/env bash
# Screen an existing pancreatic Nanopore cDNA FASTA cohort for one HUGO gene.
#
# The helper is deliberately gene-agnostic:
# - it resolves a HUGO symbol to a compact NCBI GenBank genomic locus,
# - loads that locus into one base GENtle state,
# - copies the base state per sample,
# - runs RNA-read seed interpretation in parallel across samples, and
# - keeps interim checkpoints, logs, audit JSON, TSV exports, and figure-ready
#   summaries under one timestamped run directory.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
DEFAULT_REPO_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
SCRIPT_PATH="$SCRIPT_DIR/$(basename "$0")"

GENTLE_REPO="${GENTLE_REPO:-${GENTLE_RS:-$DEFAULT_REPO_ROOT}}"
GENTLE_BIN="${GENTLE_BIN:-$GENTLE_REPO/target/release/gentle_cli}"
WORK_ROOT="${WORK_ROOT:-${WORK:-/home/clawbio/work/tp73_pancreas_benchmark}}"
COHORT_ROOT="${COHORT_ROOT:-$WORK_ROOT/cohort}"
FASTA_ROOT="${FASTA_ROOT:-$WORK_ROOT/fastq}"
DEFAULT_MANIFEST="${DEFAULT_MANIFEST:-$COHORT_ROOT/manifests/pancreas_runs.tsv}"
OUT_ROOT="${OUT_ROOT:-$WORK_ROOT/gene_screens}"

JOBS="${JOBS:-1}"
NCBI_FLANK_BP="${NCBI_FLANK_BP:-0}"
SCOPE="${SCOPE:-all_overlapping_any_strand}"
PROFILE="${PROFILE:-nanopore_cdna_v1}"
INPUT_FORMAT="${INPUT_FORMAT:-fasta}"
ORIGIN_MODE="${ORIGIN_MODE:-single_gene}"
REPORT_MODE="${REPORT_MODE:-seed_passed_only}"
COMPLETE_RULE="${COMPLETE_RULE:-near}"
SCREEN_PHASE="${SCREEN_PHASE:-seed_only}"
ALIGN_SELECTION="${ALIGN_SELECTION:-seed_passed}"
ALIGN_BAND_BP="${ALIGN_BAND_BP:-24}"
ALIGN_MIN_IDENTITY="${ALIGN_MIN_IDENTITY:-0.60}"
MAX_SECONDARY_MAPPINGS="${MAX_SECONDARY_MAPPINGS:-5}"
CHECKPOINT_EVERY_READS="${CHECKPOINT_EVERY_READS:-100000}"
CONCATEMER_LIMIT="${CONCATEMER_LIMIT:-0}"
OPTIMIZE_PARAMETERS="${OPTIMIZE_PARAMETERS:-1}"
MAX_CONTROL_MATCH_PROBABILITY="${MAX_CONTROL_MATCH_PROBABILITY:-0.0}"
AUTO_FIXTURE="${AUTO_FIXTURE:-1}"
AUTO_FETCH_FIXTURES="${AUTO_FETCH_FIXTURES:-0}"
FIXTURE_DIR="${FIXTURE_DIR:-$GENTLE_REPO/test_files/fixtures/mapping}"
FIXTURE_SPECIES="${FIXTURE_SPECIES:-homo_sapiens}"
FIXTURE_SPECIES_LABEL="${FIXTURE_SPECIES_LABEL:-human}"
DEFAULT_SEED_FRAGMENT="${DEFAULT_SEED_FRAGMENT:---kmer-len 10 --seed-stride-bp 1 --min-seed-hit-fraction 0.300 --min-weighted-seed-hit-fraction 0.050 --min-unique-matched-kmers 12 --min-chain-consistency-fraction 0.40 --max-median-transcript-gap 4.00 --min-confirmed-transitions 1 --min-transition-support-fraction 0.05 --cdna-poly-t-flip --poly-t-prefix-min-bp 18}"

GENBANK_PATH=""
RUN_ROOT=""
MANIFEST="$DEFAULT_MANIFEST"
GENE_ID=""
SEQ_ID=""
GENE_LOWER=""
GENE_SAFE=""
BASE_STATE=""
SEED_FEATURE_ID=""
SEED_ARGS_FILE=""

declare -a MUST_PASS_FASTAS=()
declare -a POSITIVE_FASTAS=()
declare -a CONTROL_FASTAS=()
declare -a MUST_PASS_GENES=()
declare -a POSITIVE_GENES=()
declare -a CONTROL_GENES=()

usage() {
  cat <<'EOF'
Usage:
  scripts/pancreas_gene_rna_screen.sh run GENE --jobs N [options]
  scripts/pancreas_gene_rna_screen.sh group-plan GROUP [options]
  scripts/pancreas_gene_rna_screen.sh status RUN_ROOT
  scripts/pancreas_gene_rna_screen.sh summarize RUN_ROOT
  scripts/pancreas_gene_rna_screen.sh monitor RUN_ROOT [SECONDS]
  scripts/pancreas_gene_rna_screen.sh print-monitor RUN_ROOT

Examples:
  scripts/pancreas_gene_rna_screen.sh run E2F1 --jobs 4
  scripts/pancreas_gene_rna_screen.sh run POU2F1 --jobs 4
  scripts/pancreas_gene_rna_screen.sh run TP73 --jobs 4 --with-alignment --align-selection seed_passed
  scripts/pancreas_gene_rna_screen.sh group-plan p53_family --jobs 4
  scripts/pancreas_gene_rna_screen.sh status /home/clawbio/work/tp73_pancreas_benchmark/gene_screens/e2f1_pancreas_screen_...

Default inputs:
  GENTLE_REPO=/home/clawbio/GENtle
  GENTLE_BIN=$GENTLE_REPO/target/release/gentle_cli
  WORK_ROOT=/home/clawbio/work/tp73_pancreas_benchmark
  FASTA_ROOT=$WORK_ROOT/fastq
  MANIFEST=$WORK_ROOT/cohort/manifests/pancreas_runs.tsv
  OUT_ROOT=$WORK_ROOT/gene_screens

Options for run:
  --jobs N                         Number of sample workers to run in parallel.
  --manifest PATH                  Cohort TSV with run_accession/sample columns.
  --fasta-root PATH                Directory containing RUN.split-spot.fasta.gz.
  --out-root PATH                  Parent directory for timestamped run roots.
  --run-root PATH                  Explicit run root; must not already exist.
  --genbank-path PATH              Use an existing GenBank locus instead of NCBI.
  --seq-id ID                      GENtle sequence id, default lower(GENE)_ncbi.
  --flank-bp N                     Expand NCBI locus by N bp on both sides.
  --must-pass-transcript-fasta P   Repeatable positive isoform panel for preflight.
  --positive-transcript-fasta P    Repeatable softer positive panel for preflight.
  --control-transcript-fasta P     Repeatable control panel for preflight.
  --must-pass-gene GENE            Fetch/use Ensembl cDNA fixture as must-pass panel.
  --positive-gene GENE             Fetch/use Ensembl cDNA fixture as softer positive panel.
  --control-gene GENE              Fetch/use Ensembl cDNA fixture as control panel.
  --auto-fetch-fixtures            Fetch missing Ensembl cDNA fixtures through GENtle.
  --fixture-dir DIR                Fixture FASTA directory; default test_files/fixtures/mapping.
  --no-auto-fixture                Do not auto-use test_files/.../ensembl_human_gene_all.fasta.
  --no-optimize-parameters         Run preflight for reporting but use default seed args.
  --seed-args "..."                Explicit seed-filter fragment.
  --seed-only                      Interpret reads and summarize seed evidence only (default).
  --with-alignment                 Run phase-2 retained-read alignment and downstream reports.
  --align-selection MODE           all|seed_passed|aligned for --with-alignment; default seed_passed.
  --align-min-identity F           Default 0.60.
  --max-secondary-mappings N       Default 5.
  --concatemer-limit N             Optional per-sample concatemer audit; default 0.

Options for group-plan:
  --catalog PATH                   Optional gene-group catalog/overlay path.
  --out-dir PATH                   Directory for plan artifacts; default OUT_ROOT/group_plans.
  --jobs N                         Worker count embedded in generated per-gene commands.
  --manifest PATH                  Cohort TSV embedded in generated commands.
  --fasta-root PATH                FASTA root embedded in generated commands.
  --out-root PATH                  Parent for generated per-gene run roots.
  --seed-only                      Generate seed-only commands (default).
  --with-alignment                 Generate phase-2 alignment commands.
  --align-selection MODE           all|seed_passed|aligned for --with-alignment.

By default, gene screens stop after the seed phase. This keeps cohort-scale
triage fast and avoids accidentally aligning every retained row for abundant
or broadly seeded genes. Use --with-alignment only after a seed-only run has
identified a small enough evidence set worth phase-2 alignment.

The script uses existing FASTA.gz files. It intentionally does not download SRA
or recreate FASTA; use the TP73 cohort helper for those data-management steps.

Gene-group plans are review helpers, not automatic trusted curation. They
consume the shared GENtle gene-group catalog layer, write the resolved group
record beside the generated commands, and keep execution explicit so project
overlays, Gene Ontology mappings, and lab-facing local groups remain auditable.
EOF
}

timestamp() {
  date '+[%Y-%m-%d %H:%M:%S %Z]'
}

log() {
  echo "$(timestamp) $*" >&2
}

die() {
  echo "$(timestamp) ERROR: $*" >&2
  exit 1
}

trap_rm_rf_on_exit() {
  local path="$1"
  local quoted
  printf -v quoted '%q' "$path"
  trap "rm -rf -- $quoted" EXIT
}

quote_command_line() {
  local arg
  for arg in "$@"; do
    printf '%q ' "$arg"
  done
  printf '\n'
}

gentle_shell_quoted_token() {
  local token="$1"
  token="${token//\\/\\\\}"
  token="${token//\"/\\\"}"
  printf '"%s"' "$token"
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

safe_token() {
  printf '%s' "$1" | tr -c 'A-Za-z0-9_.-' '_'
}

lower_token() {
  printf '%s' "$1" | tr '[:upper:]' '[:lower:]'
}

fixture_token() {
  printf '%s' "$1" | tr '[:upper:]' '[:lower:]' | tr -c 'a-z0-9_' '_'
}

upper_token() {
  printf '%s' "$1" | tr '[:lower:]' '[:upper:]'
}

init_gene_tokens() {
  GENE_ID="$(upper_token "$GENE_ID")"
  GENE_LOWER="$(lower_token "$GENE_ID")"
  GENE_SAFE="$(safe_token "$GENE_LOWER")"
  SEQ_ID="${SEQ_ID:-${GENE_SAFE}_ncbi}"
}

parse_run_args() {
  [ $# -ge 1 ] || die "run requires a HUGO gene symbol"
  GENE_ID="$1"
  shift
  init_gene_tokens

  while [ $# -gt 0 ]; do
    case "$1" in
      --jobs)
        JOBS="${2:-}"
        shift 2
        ;;
      --manifest)
        MANIFEST="${2:-}"
        shift 2
        ;;
      --fasta-root)
        FASTA_ROOT="${2:-}"
        shift 2
        ;;
      --out-root)
        OUT_ROOT="${2:-}"
        shift 2
        ;;
      --run-root)
        RUN_ROOT="${2:-}"
        shift 2
        ;;
      --genbank-path)
        GENBANK_PATH="${2:-}"
        shift 2
        ;;
      --seq-id)
        SEQ_ID="${2:-}"
        shift 2
        ;;
      --flank-bp)
        NCBI_FLANK_BP="${2:-}"
        shift 2
        ;;
      --must-pass-transcript-fasta)
        MUST_PASS_FASTAS+=("${2:-}")
        shift 2
        ;;
      --positive-transcript-fasta)
        POSITIVE_FASTAS+=("${2:-}")
        shift 2
        ;;
      --control-transcript-fasta)
        CONTROL_FASTAS+=("${2:-}")
        shift 2
        ;;
      --must-pass-gene|--must-pass-transcript-gene)
        MUST_PASS_GENES+=("${2:-}")
        shift 2
        ;;
      --positive-gene|--positive-transcript-gene)
        POSITIVE_GENES+=("${2:-}")
        shift 2
        ;;
      --control-gene|--control-transcript-gene)
        CONTROL_GENES+=("${2:-}")
        shift 2
        ;;
      --auto-fetch-fixtures)
        AUTO_FETCH_FIXTURES=1
        shift
        ;;
      --fixture-dir)
        FIXTURE_DIR="${2:-}"
        shift 2
        ;;
      --no-auto-fixture)
        AUTO_FIXTURE=0
        shift
        ;;
      --no-optimize-parameters)
        OPTIMIZE_PARAMETERS=0
        shift
        ;;
      --seed-args)
        DEFAULT_SEED_FRAGMENT="${2:-}"
        OPTIMIZE_PARAMETERS=0
        shift 2
        ;;
      --seed-only)
        SCREEN_PHASE="seed_only"
        shift
        ;;
      --with-alignment)
        SCREEN_PHASE="with_alignment"
        shift
        ;;
      --align-selection)
        ALIGN_SELECTION="${2:-}"
        shift 2
        ;;
      --align-min-identity)
        ALIGN_MIN_IDENTITY="${2:-}"
        shift 2
        ;;
      --max-secondary-mappings)
        MAX_SECONDARY_MAPPINGS="${2:-}"
        shift 2
        ;;
      --concatemer-limit)
        CONCATEMER_LIMIT="${2:-}"
        shift 2
        ;;
      -h|--help|help)
        usage
        exit 0
        ;;
      *)
        die "Unknown run option: $1"
        ;;
    esac
  done

  case "$JOBS" in
    ''|*[!0-9]*)
      die "--jobs must be a positive integer"
      ;;
  esac
  [ "$JOBS" -ge 1 ] || die "--jobs must be at least 1"
  case "$NCBI_FLANK_BP" in
    ''|*[!0-9]*)
      die "--flank-bp must be a non-negative integer"
      ;;
  esac
  case "$CONCATEMER_LIMIT" in
    ''|*[!0-9]*)
      die "--concatemer-limit must be a non-negative integer"
      ;;
  esac
  case "$SCREEN_PHASE" in
    seed_only|with_alignment)
      ;;
    *)
      die "SCREEN_PHASE must be seed_only or with_alignment"
      ;;
  esac
  case "$ALIGN_SELECTION" in
    all|seed_passed|aligned)
      ;;
    *)
      die "--align-selection must be one of: all, seed_passed, aligned"
      ;;
  esac

  if [ -z "$RUN_ROOT" ]; then
    local stamp
    stamp="$(date -u +%Y%m%dT%H%M%SZ)"
    RUN_ROOT="$OUT_ROOT/${GENE_SAFE}_pancreas_screen_${stamp}"
  fi
  BASE_STATE="$RUN_ROOT/${GENE_SAFE}_pancreas.base.gentle.json"
  SEED_ARGS_FILE="$RUN_ROOT/manifests/${GENE_SAFE}_seed_args.txt"
}

write_gene_group_plan() {
  local group_token="$1"
  local catalog_path="$2"
  local out_dir="$3"
  local group_safe show_json genes_txt commands_sh plan_md group_arg shell_cmd group_label group_id source_path source_scope curation_status member_count gene group_run_out_root
  group_safe="$(safe_token "$(lower_token "$group_token")")"
  show_json="$out_dir/${group_safe}.gene_group.show.json"
  genes_txt="$out_dir/${group_safe}.genes.txt"
  commands_sh="$out_dir/${group_safe}.commands.sh"
  plan_md="$out_dir/${group_safe}.plan.md"
  group_run_out_root="$OUT_ROOT/gene_groups/$group_safe"

  mkdir -p "$out_dir"
  group_arg="$(gentle_shell_quoted_token "$group_token")"
  shell_cmd="gene-groups show $group_arg"
  if [ -n "$catalog_path" ]; then
    shell_cmd="$shell_cmd --catalog $(gentle_shell_quoted_token "$catalog_path")"
  fi

  log "Resolve gene group through shared GENtle catalog layer: $group_token"
  "$GENTLE_BIN" shell "$shell_cmd" > "$show_json"

  jq -r '
    .group.members[]?
    | select((.status // "included") == "included")
    | .symbol
    | select(. != null and . != "")
  ' "$show_json" | awk '!seen[$0]++' > "$genes_txt"

  if [ ! -s "$genes_txt" ]; then
    die "Gene group resolved but had no included member symbols: $group_token"
  fi

  group_label="$(jq -r '.group.label // .group.id // empty' "$show_json")"
  group_id="$(jq -r '.group.id // empty' "$show_json")"
  source_path="$(jq -r '.source_path // ""' "$show_json")"
  source_scope="$(jq -r '.source_scope // ""' "$show_json")"
  curation_status="$(jq -r '.group.curation_status // ""' "$show_json")"
  member_count="$(wc -l < "$genes_txt" | tr -d ' ')"

  {
    printf '#!/usr/bin/env bash\n'
    printf 'set -euo pipefail\n\n'
    printf 'cd %q\n\n' "$GENTLE_REPO"
    printf '# Gene-group pancreas screen plan generated by %s\n' "$SCRIPT_PATH"
    printf '# group_id=%s label=%s members=%s source=%s (%s)\n\n' \
      "$group_id" "$group_label" "$member_count" "$source_path" "$source_scope"
    while IFS= read -r gene; do
      [ -n "$gene" ] || continue
      local -a cmd
      cmd=("$SCRIPT_PATH" run "$gene" \
        --jobs "$JOBS" \
        --manifest "$MANIFEST" \
        --fasta-root "$FASTA_ROOT" \
        --out-root "$group_run_out_root")
      if [ "$SCREEN_PHASE" = "with_alignment" ]; then
        cmd+=(--with-alignment --align-selection "$ALIGN_SELECTION")
      else
        cmd+=(--seed-only)
      fi
      quote_command_line "${cmd[@]}"
    done < "$genes_txt"
  } > "$commands_sh"
  chmod +x "$commands_sh"

  {
    printf '# Pancreas gene-group screen plan: %s\n\n' "${group_label:-$group_token}"
    printf 'This is a reviewable command plan for the existing gene-agnostic pancreas RNA-read screen helper. It does not run the screen automatically.\n\n'
    printf -- '- group query: `%s`\n' "$group_token"
    printf -- '- resolved group id: `%s`\n' "${group_id:-n/a}"
    printf -- '- curation status: `%s`\n' "${curation_status:-n/a}"
    printf -- '- source: `%s` (`%s`)\n' "${source_path:-n/a}" "${source_scope:-n/a}"
    printf -- '- member count used for commands: `%s`\n' "$member_count"
    printf -- '- member list: `%s`\n' "$genes_txt"
    printf -- '- command script: `%s`\n' "$commands_sh"
    printf -- '- resolved group record: `%s`\n\n' "$show_json"
    printf '## Members\n\n'
    while IFS= read -r gene; do
      [ -n "$gene" ] || continue
      printf -- '- `%s`\n' "$gene"
    done < "$genes_txt"
    printf '\n## Execution\n\n'
    printf 'Inspect the group record and member list first. Then run individual commands from `%s`, or execute the script when the cohort FASTA inputs and resource choices are ready.\n\n' "$commands_sh"
    printf 'The generated commands use the reviewed GENtle group membership, but each per-gene screen still resolves its own locus/transcript resources flexibly through the helper options and local overlays.\n'
  } > "$plan_md"

  log "Wrote gene-group record: $show_json"
  log "Wrote member list:       $genes_txt"
  log "Wrote command script:    $commands_sh"
  log "Wrote review plan:       $plan_md"
}

parse_group_plan_args() {
  [ $# -ge 1 ] || die "group-plan requires a gene-group id or alias"
  local group_token="$1"
  local catalog_path=""
  local out_dir=""
  shift

  while [ $# -gt 0 ]; do
    case "$1" in
      --catalog)
        catalog_path="${2:-}"
        shift 2
        ;;
      --out-dir)
        out_dir="${2:-}"
        shift 2
        ;;
      --jobs)
        JOBS="${2:-}"
        shift 2
        ;;
      --manifest)
        MANIFEST="${2:-}"
        shift 2
        ;;
      --fasta-root)
        FASTA_ROOT="${2:-}"
        shift 2
        ;;
      --out-root)
        OUT_ROOT="${2:-}"
        shift 2
        ;;
      --seed-only)
        SCREEN_PHASE="seed_only"
        shift
        ;;
      --with-alignment)
        SCREEN_PHASE="with_alignment"
        shift
        ;;
      --align-selection)
        ALIGN_SELECTION="${2:-}"
        shift 2
        ;;
      -h|--help|help)
        usage
        exit 0
        ;;
      *)
        die "Unknown group-plan option: $1"
        ;;
    esac
  done

  case "$JOBS" in
    ''|*[!0-9]*)
      die "--jobs must be a positive integer"
      ;;
  esac
  [ "$JOBS" -ge 1 ] || die "--jobs must be at least 1"
  case "$SCREEN_PHASE" in
    seed_only|with_alignment)
      ;;
    *)
      die "SCREEN_PHASE must be seed_only or with_alignment"
      ;;
  esac
  case "$ALIGN_SELECTION" in
    all|seed_passed|aligned)
      ;;
    *)
      die "--align-selection must be one of: all, seed_passed, aligned"
      ;;
  esac

  out_dir="${out_dir:-$OUT_ROOT/group_plans}"
  require_executable "$GENTLE_BIN"
  require_executable jq
  write_gene_group_plan "$group_token" "$catalog_path" "$out_dir"
}

require_common_tools() {
  cd "$GENTLE_REPO" || die "Cannot cd to GENTLE_REPO: $GENTLE_REPO"
  require_executable "$GENTLE_BIN"
  require_executable jq
  require_executable curl
  require_executable gzip
  require_executable /usr/bin/time
}

curl_json_get() {
  local output="$1"
  shift
  curl -fsSL --get "$@" > "$output"
}

resolve_ncbi_gene_locus() {
  local resources_dir="$RUN_ROOT/resources"
  local reports_dir="$RUN_ROOT/reports"
  local search_json="$reports_dir/${GENE_SAFE}.ncbi_gene_search.json"
  local summary_json="$reports_dir/${GENE_SAFE}.ncbi_gene_summary.json"
  local locus_json="$reports_dir/${GENE_SAFE}.ncbi_locus.json"
  local gene_uid chraccver chrstart chrstop lo hi seq_start seq_stop
  mkdir -p "$resources_dir" "$reports_dir"

  if [ -n "$GENBANK_PATH" ]; then
    require_file "$GENBANK_PATH"
    cp "$GENBANK_PATH" "$resources_dir/${GENE_SAFE}.ncbi.gb"
    GENBANK_PATH="$resources_dir/${GENE_SAFE}.ncbi.gb"
    log "Using provided GenBank locus: $GENBANK_PATH"
    return 0
  fi

  GENBANK_PATH="$resources_dir/${GENE_SAFE}.ncbi.gb"
  if [ -s "$GENBANK_PATH" ]; then
    log "Reusing downloaded GenBank locus: $GENBANK_PATH"
    return 0
  fi

  log "Resolve $GENE_ID via NCBI Gene"
  curl_json_get "$search_json" \
    "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi" \
    --data-urlencode "db=gene" \
    --data-urlencode "retmode=json" \
    --data-urlencode "retmax=10" \
    --data-urlencode "term=Homo sapiens[Organism] AND ${GENE_ID}[Gene Name] AND alive[property]"

  gene_uid="$(jq -r '.esearchresult.idlist[0] // empty' "$search_json")"
  [ -n "$gene_uid" ] || die "NCBI Gene search returned no ids for $GENE_ID"

  curl_json_get "$summary_json" \
    "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi" \
    --data-urlencode "db=gene" \
    --data-urlencode "retmode=json" \
    --data-urlencode "id=$gene_uid"

  gene_uid="$(
    jq -r --arg gene "$GENE_ID" '
      first(
        .result.uids[] as $uid
        | .result[$uid] as $row
        | select((($row.nomenclaturesymbol // $row.name // "") | ascii_upcase) == $gene)
        | $uid
      ) // (.result.uids[0] // empty)
    ' "$summary_json"
  )"
  [ -n "$gene_uid" ] || die "Could not choose an NCBI Gene id for $GENE_ID"

  chraccver="$(
    jq -r --arg id "$gene_uid" '
      ([.result[$id].genomicinfo[]? | select((.chraccver // "") | startswith("NC_"))][0]
       // .result[$id].genomicinfo[0]
       // {}) .chraccver // empty
    ' "$summary_json"
  )"
  chrstart="$(
    jq -r --arg id "$gene_uid" '
      ([.result[$id].genomicinfo[]? | select((.chraccver // "") | startswith("NC_"))][0]
       // .result[$id].genomicinfo[0]
       // {}) .chrstart // empty
    ' "$summary_json"
  )"
  chrstop="$(
    jq -r --arg id "$gene_uid" '
      ([.result[$id].genomicinfo[]? | select((.chraccver // "") | startswith("NC_"))][0]
       // .result[$id].genomicinfo[0]
       // {}) .chrstop // empty
    ' "$summary_json"
  )"
  [ -n "$chraccver" ] || die "NCBI Gene summary has no primary chromosome accession for $GENE_ID"
  [ -n "$chrstart" ] || die "NCBI Gene summary has no chrstart for $GENE_ID"
  [ -n "$chrstop" ] || die "NCBI Gene summary has no chrstop for $GENE_ID"

  if [ "$chrstart" -le "$chrstop" ]; then
    lo="$chrstart"
    hi="$chrstop"
  else
    lo="$chrstop"
    hi="$chrstart"
  fi
  seq_start=$((lo + 1 - NCBI_FLANK_BP))
  seq_stop=$((hi + 1 + NCBI_FLANK_BP))
  [ "$seq_start" -ge 1 ] || seq_start=1

  jq -n \
    --arg gene "$GENE_ID" \
    --arg gene_uid "$gene_uid" \
    --arg chraccver "$chraccver" \
    --argjson chrstart "$chrstart" \
    --argjson chrstop "$chrstop" \
    --argjson seq_start "$seq_start" \
    --argjson seq_stop "$seq_stop" \
    --argjson flank_bp "$NCBI_FLANK_BP" \
    '{
      gene: $gene,
      ncbi_gene_id: $gene_uid,
      chraccver: $chraccver,
      ncbi_chrstart_raw: $chrstart,
      ncbi_chrstop_raw: $chrstop,
      efetch_seq_start_1based: $seq_start,
      efetch_seq_stop_1based: $seq_stop,
      flank_bp: $flank_bp,
      note: "NCBI Gene chrstart/chrstop are converted conservatively to 1-based efetch coordinates."
    }' > "$locus_json"

  log "Fetch $GENE_ID GenBank locus: $chraccver:$seq_start..$seq_stop"
  curl -fsSL --get "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi" \
    --data-urlencode "db=nuccore" \
    --data-urlencode "id=$chraccver" \
    --data-urlencode "seq_start=$seq_start" \
    --data-urlencode "seq_stop=$seq_stop" \
    --data-urlencode "rettype=gbwithparts" \
    --data-urlencode "retmode=text" \
    > "$GENBANK_PATH"

  require_file "$GENBANK_PATH"
  grep -q '^ORIGIN' "$GENBANK_PATH" || die "Downloaded GenBank file has no ORIGIN: $GENBANK_PATH"
  grep -q "/gene=\"$GENE_ID\"" "$GENBANK_PATH" \
    || log "WARNING: Downloaded GenBank file does not contain an exact /gene=\"$GENE_ID\" qualifier"
  log "Downloaded GenBank locus: $GENBANK_PATH ($(wc -c < "$GENBANK_PATH") bytes)"
}

write_load_workflow() {
  local workflow="$RUN_ROOT/manifests/load_${GENE_SAFE}_ncbi.workflow.json"
  jq -n \
    --arg path "$GENBANK_PATH" \
    --arg seq_id "$SEQ_ID" \
    --arg run_id "load_${GENE_SAFE}_ncbi_for_pancreas_gene_screen" \
    '{run_id: $run_id, ops: [{LoadFile: {path: $path, as_id: $seq_id}}]}' \
    > "$workflow"
  printf '%s\n' "$workflow"
}

load_base_state() {
  local workflow
  local -a gentle
  gentle=("$GENTLE_BIN" --state "$BASE_STATE" --progress-stderr)

  if [ -s "$BASE_STATE" ]; then
    log "Base state already exists: $BASE_STATE"
    return 0
  fi

  workflow="$(write_load_workflow)"
  log "Load $GENE_ID locus into base GENtle state"
  "${gentle[@]}" workflow @"$workflow" \
    > "$RUN_ROOT/reports/load_${GENE_SAFE}_ncbi.command.json" \
    2> "$RUN_ROOT/logs/load_${GENE_SAFE}_ncbi.stderr.log"
}

resolve_seed_feature() {
  local -a gentle
  gentle=("$GENTLE_BIN" --state "$BASE_STATE" --progress-stderr)

  log "Resolve $GENE_ID seed feature id"
  "${gentle[@]}" features query "$SEQ_ID" \
    --kind gene --label "$GENE_ID" --limit 20 \
    > "$RUN_ROOT/reports/${GENE_SAFE}.gene_features.json" \
    2> "$RUN_ROOT/logs/${GENE_SAFE}.gene_features.stderr.log"

  SEED_FEATURE_ID="$(
    jq -r --arg gene "$GENE_ID" '
      first(.rows[]? | select((.label // "" | ascii_upcase) == $gene) | .feature_id)
      // (.rows[0].feature_id // empty)
    ' "$RUN_ROOT/reports/${GENE_SAFE}.gene_features.json"
  )"
  [ -n "$SEED_FEATURE_ID" ] && [ "$SEED_FEATURE_ID" != "null" ] \
    || die "Could not resolve seed feature id for $GENE_ID; inspect $RUN_ROOT/reports/${GENE_SAFE}.gene_features.json"
  printf '%s\n' "$SEED_FEATURE_ID" > "$RUN_ROOT/manifests/seed_feature_id.txt"
  log "SEED_FEATURE_ID=$SEED_FEATURE_ID"
}

prepare_preflight_fastas() {
  ensure_transcript_fixture() {
    local gene="$1"
    local role="$2"
    local species_token gene_token fixture_path fixture_fetcher fixture_work stdout_log stderr_log
    species_token="$(fixture_token "$FIXTURE_SPECIES_LABEL")"
    gene_token="$(fixture_token "$gene")"
    fixture_path="$FIXTURE_DIR/ensembl_${species_token}_${gene_token}_all.fasta"
    if [ -s "$fixture_path" ]; then
      printf '%s\n' "$fixture_path"
      return 0
    fi
    [ "$AUTO_FETCH_FIXTURES" = "1" ] \
      || die "Missing $role transcript fixture for $gene: $fixture_path (use --auto-fetch-fixtures to retrieve it)"
    fixture_fetcher="$SCRIPT_DIR/fetch_ensembl_cdna_fixtures.sh"
    [ -x "$fixture_fetcher" ] || die "Missing executable fixture fetcher: $fixture_fetcher"
    fixture_work="$RUN_ROOT/resources/ensembl_cdna_fixtures"
    stdout_log="$RUN_ROOT/logs/fetch_fixture_$(fixture_token "$gene").stdout.log"
    stderr_log="$RUN_ROOT/logs/fetch_fixture_$(fixture_token "$gene").stderr.log"
    log "Auto-fetch Ensembl cDNA fixture for $gene ($role): $fixture_path"
    if ! "$fixture_fetcher" \
      --gene "$gene" \
      --species "$FIXTURE_SPECIES" \
      --species-label "$FIXTURE_SPECIES_LABEL" \
      --out-dir "$FIXTURE_DIR" \
      --work-dir "$fixture_work" \
      --gentle-bin "$GENTLE_BIN" \
      > "$stdout_log" \
      2> "$stderr_log"; then
      log "Fixture fetch for $gene failed; stdout=$stdout_log stderr=$stderr_log"
      tail -n 40 "$stderr_log" >&2 || true
      die "Could not fetch Ensembl cDNA fixture for $gene"
    fi
    if [ ! -s "$fixture_path" ]; then
      log "Fixture fetch for $gene completed but did not write $fixture_path"
      log "Fetcher stdout: $stdout_log"
      tail -n 40 "$stdout_log" >&2 || true
      log "Fetcher stderr: $stderr_log"
      tail -n 40 "$stderr_log" >&2 || true
      log "Matching fixture files currently under $FIXTURE_DIR:"
      find "$FIXTURE_DIR" -maxdepth 1 -type f -iname "*${gene_token}*" -print >&2 || true
      die "Missing fetched Ensembl cDNA fixture for $gene"
    fi
    printf '%s\n' "$fixture_path"
  }

  local auto_fixture gene fasta resolved_fixture
  auto_fixture="$FIXTURE_DIR/ensembl_$(fixture_token "$FIXTURE_SPECIES_LABEL")_$(fixture_token "$GENE_ID")_all.fasta"
  if [ "$AUTO_FIXTURE" = "1" ]; then
    if [ -s "$auto_fixture" ] || [ "$AUTO_FETCH_FIXTURES" = "1" ]; then
      resolved_fixture="$(ensure_transcript_fixture "$GENE_ID" "must-pass")"
      MUST_PASS_FASTAS+=("$resolved_fixture")
      log "Auto-added transcript fixture for $GENE_ID: $resolved_fixture"
    fi
  fi

  for gene in "${MUST_PASS_GENES[@]}"; do
    [ -n "$gene" ] || continue
    MUST_PASS_FASTAS+=("$(ensure_transcript_fixture "$gene" "must-pass")")
  done
  for gene in "${POSITIVE_GENES[@]}"; do
    [ -n "$gene" ] || continue
    POSITIVE_FASTAS+=("$(ensure_transcript_fixture "$gene" "positive")")
  done
  for gene in "${CONTROL_GENES[@]}"; do
    [ -n "$gene" ] || continue
    CONTROL_FASTAS+=("$(ensure_transcript_fixture "$gene" "control")")
  done

  for fasta in "${MUST_PASS_FASTAS[@]}" "${POSITIVE_FASTAS[@]}" "${CONTROL_FASTAS[@]}"; do
    [ -n "$fasta" ] || continue
    require_file "$fasta"
  done
}

run_preflight() {
  local -a gentle cmd
  local fasta seed_fragment optimize_preflight preflight_json preflight_stderr preflight_command
  gentle=("$GENTLE_BIN" --state "$BASE_STATE" --progress-stderr)
  prepare_preflight_fastas
  optimize_preflight="$OPTIMIZE_PARAMETERS"
  if [ "$OPTIMIZE_PARAMETERS" = "1" ] && [ "${#CONTROL_FASTAS[@]}" -eq 0 ]; then
    optimize_preflight=0
    log "No control transcript FASTA provided for $GENE_ID; running one-pass preflight without parameter optimization"
  fi

  cmd=("${gentle[@]}" rna-reads preflight-isoforms "$SEQ_ID" "$SEED_FEATURE_ID" --scope "$SCOPE")
  for fasta in "${MUST_PASS_FASTAS[@]}"; do
    cmd+=(--must-pass-transcript-fasta "$fasta")
  done
  for fasta in "${POSITIVE_FASTAS[@]}"; do
    cmd+=(--positive-transcript-fasta "$fasta")
  done
  for fasta in "${CONTROL_FASTAS[@]}"; do
    cmd+=(--control-transcript-fasta "$fasta")
  done
  if [ "$optimize_preflight" = "1" ]; then
    cmd+=(--optimize-parameters --max-control-match-probability "$MAX_CONTROL_MATCH_PROBABILITY")
  fi

  log "Run $GENE_ID isoform preflight"
  preflight_json="$RUN_ROOT/reports/${GENE_SAFE}.preflight.json"
  preflight_stderr="$RUN_ROOT/logs/${GENE_SAFE}.preflight.stderr.log"
  preflight_command="$RUN_ROOT/reports/${GENE_SAFE}.preflight.command.txt"
  quote_command_line "${cmd[@]}" > "$preflight_command"
  log "Preflight command: $preflight_command"
  if ! "${cmd[@]}" > "$preflight_json" 2> "$preflight_stderr"; then
    log "ERROR: $GENE_ID preflight failed; stderr tail follows ($preflight_stderr)"
    tail -n 80 "$preflight_stderr" >&2 || true
    die "$GENE_ID isoform preflight failed"
  fi

  seed_fragment="$(
    jq -r '.threshold_recommendation.seed_filter_cli_fragment // empty' \
      "$RUN_ROOT/reports/${GENE_SAFE}.preflight.json"
  )"
  if [ -z "$seed_fragment" ] || [ "$optimize_preflight" != "1" ]; then
    seed_fragment="$DEFAULT_SEED_FRAGMENT"
    log "Using default/explicit seed args for $GENE_ID"
  else
    log "Using optimized seed args for $GENE_ID"
  fi

  [ -n "$seed_fragment" ] || die "Seed argument fragment is empty"
  printf '%s\n' "$seed_fragment" > "$SEED_ARGS_FILE"

  jq '{
    target: [.target_passed_transcript_count, .target_transcript_count],
    positives: [.positive_control_passed_transcript_count, .positive_control_transcript_count],
    controls: [(.control_summaries // [])[] | {
      control_id,
      transcript_count,
      passed_transcript_count,
      weighted_pass_probability,
      best_transcript_id
    }],
    threshold_recommendation
  }' "$RUN_ROOT/reports/${GENE_SAFE}.preflight.json" \
    > "$RUN_ROOT/reports/${GENE_SAFE}.preflight.summary.json"
  log "Seed args: $(cat "$SEED_ARGS_FILE")"
}

report_id_for() {
  local run="$1"
  local sample_id="$2"
  local safe_sample
  safe_sample="$(safe_token "$sample_id")"
  if [ "$safe_sample" = "$run" ]; then
    printf '%s_pancreas_%s\n' "$GENE_SAFE" "$run"
  else
    printf '%s_pancreas_%s_%s\n' "$GENE_SAFE" "$safe_sample" "$run"
  fi
}

read_fasta_for() {
  local run="$1"
  printf '%s/%s.split-spot.fasta.gz\n' "$FASTA_ROOT" "$run"
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

write_sample_manifest_row() {
  local run="$1"
  local sample_id="$2"
  local sample_name="$3"
  local note="$4"
  local fasta report_id
  fasta="$(read_fasta_for "$run")"
  report_id="$(report_id_for "$run" "$sample_id")"
  if [ -s "$fasta" ]; then
    printf '%s\t%s\t%s\t%s\t%s\t%s\n' \
      "$run" "$sample_id" "$sample_name" "$note" "$fasta" "$report_id"
  else
    log "Skipping $run for $GENE_ID because FASTA is missing: $fasta"
  fi
}

build_sample_manifest() {
  local sample_manifest="$RUN_ROOT/manifests/${GENE_SAFE}_pancreas_inputs.tsv"
  {
    printf 'run_accession\tsample_id\tsample_name\tnote\tinput_path\treport_id\n'
    for_manifest_rows write_sample_manifest_row
  } > "$sample_manifest"
  if [ "$(tail -n +2 "$sample_manifest" | wc -l | tr -d ' ')" = "0" ]; then
    die "No usable FASTA rows found in $MANIFEST under $FASTA_ROOT"
  fi
  validate_sample_manifest "$sample_manifest"
  log "Wrote sample manifest: $sample_manifest"
  column -t -s $'\t' "$sample_manifest" >&2 || cat "$sample_manifest" >&2
}

validate_sample_manifest() {
  local sample_manifest="$1"
  local duplicate_runs
  duplicate_runs="$(
    awk -F '\t' '
      NR == 1 { next }
      $1 == "" { next }
      {
        count[$1] += 1
        rows[$1] = rows[$1] sprintf("  line %d: run=%s sample_id=%s sample_name=%s report_id=%s input=%s\n", NR, $1, $2, $3, $6, $5)
      }
      END {
        for (run in count) {
          if (count[run] > 1) {
            printf "%s\t%d\n%s", run, count[run], rows[run]
          }
        }
      }
    ' "$sample_manifest"
  )"
  if [ -n "$duplicate_runs" ]; then
    log "ERROR: duplicate run_accession values in generated sample manifest: $sample_manifest"
    printf '%s\n' "$duplicate_runs" >&2
    die "Duplicate run_accession values are not supported because each run uses one runs/RUN_ACCESSION work directory"
  fi
}

write_run_env() {
  cat > "$RUN_ROOT/run.env" <<EOF
export GENE_ID="$GENE_ID"
export GENE_LOWER="$GENE_LOWER"
export GENE_SAFE="$GENE_SAFE"
export SEQ_ID="$SEQ_ID"
export GENTLE_REPO="$GENTLE_REPO"
export GENTLE_BIN="$GENTLE_BIN"
export WORK_ROOT="$WORK_ROOT"
export COHORT_ROOT="$COHORT_ROOT"
export FASTA_ROOT="$FASTA_ROOT"
export MANIFEST="$MANIFEST"
export RUN_ROOT="$RUN_ROOT"
export BASE_STATE="$BASE_STATE"
export SEED_FEATURE_ID="$SEED_FEATURE_ID"
export SEED_ARGS_FILE="$SEED_ARGS_FILE"
export JOBS="$JOBS"
export PROFILE="$PROFILE"
export INPUT_FORMAT="$INPUT_FORMAT"
export SCOPE="$SCOPE"
export ORIGIN_MODE="$ORIGIN_MODE"
export REPORT_MODE="$REPORT_MODE"
export COMPLETE_RULE="$COMPLETE_RULE"
export SCREEN_PHASE="$SCREEN_PHASE"
export ALIGN_SELECTION="$ALIGN_SELECTION"
export ALIGN_BAND_BP="$ALIGN_BAND_BP"
export ALIGN_MIN_IDENTITY="$ALIGN_MIN_IDENTITY"
export MAX_SECONDARY_MAPPINGS="$MAX_SECONDARY_MAPPINGS"
export CHECKPOINT_EVERY_READS="$CHECKPOINT_EVERY_READS"
export CONCATEMER_LIMIT="$CONCATEMER_LIMIT"
EOF
}

sample_step() {
  local run_dir="$1"
  local step="$2"
  printf '%s\n' "$step" > "$run_dir/CURRENT_STEP.txt"
  log "$step"
}

run_optional() {
  local description="$1"
  local stderr_path="$2"
  shift 2
  if ! "$@" 2> "$stderr_path"; then
    log "WARNING: optional step failed: $description (see $stderr_path)"
    return 0
  fi
}

write_final_summary() {
  local run="$1"
  local sample_id="$2"
  local sample_name="$3"
  local note="$4"
  local report_id="$5"
  local show_json="$6"
  local gene_json="$7"
  local audit_json="$8"
  local out_json="$9"
  jq -n \
    --arg run "$run" \
    --arg sample_id "$sample_id" \
    --arg sample_name "$sample_name" \
    --arg note "$note" \
    --arg gene_id "$GENE_ID" \
    --arg report_id "$report_id" \
    --arg analysis_phase "$SCREEN_PHASE" \
    --arg align_selection "$ALIGN_SELECTION" \
    --slurpfile show "$show_json" \
    --slurpfile gene "$gene_json" \
    --slurpfile audit "$audit_json" '
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
            q90: quantile_len($counts; 0.90),
            q95: quantile_len($counts; 0.95),
            q99: quantile_len($counts; 0.99),
            q100: last_observed($counts)
          }
        };
    def accepted_status:
      .status == "accepted_fragment" or .status == "accepted_complete";
    ($show[0].report // {}) as $report
    | ($gene[0] // {}) as $gene_support
    | ($audit[0].rows // []) as $audit_rows
    | ($report.hits // []) as $hits
    | {
        run_accession: $run,
        sample_id: $sample_id,
        sample_name: $sample_name,
        note: $note,
        gene_id: $gene_id,
        analysis_phase: $analysis_phase,
        align_selection: $align_selection,
        report_id: ($report.report_id // $report_id),
        read_count_total: ($report.read_count_total // 0),
        strict_seed_passed_reads: ($report.read_count_seed_passed // 0),
        retained_report_rows: ($hits | length),
        retained_aligned_rows: ($report.read_count_aligned // 0),
        retained_msa_eligible_rows: ($report.retained_count_msa_eligible // 0),
        strict_seed_passed_aligned_rows: (
          $hits | map(select((.passed_seed_filter == true) and (.best_mapping != null))) | length
        ),
        strict_seed_passed_msa_eligible_rows: (
          $hits | map(select((.passed_seed_filter == true) and (.msa_eligible == true))) | length
        ),
        strict_seed_accepted_target_rows: (
          $audit_rows | map(select((.passed_seed_filter == true) and accepted_status)) | length
        ),
        strict_seed_complete_target_rows: (
          $audit_rows | map(select((.passed_seed_filter == true) and (.status == "accepted_complete"))) | length
        ),
        strict_seed_full_length_near_rows: (
          $audit_rows | map(select((.passed_seed_filter == true) and (.full_length_near == true))) | length
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

run_one_sample_impl() {
  local run="$1"
  local sample_id="$2"
  local sample_name="$3"
  local note="$4"
  local read_fasta="$5"
  local report_id="$6"
  local run_dir="$RUN_ROOT/runs/$run"
  local state="$run_dir/$run.${GENE_SAFE}.gentle.json"
  local show_json="$run_dir/post_interpret/json/$report_id.show_report.json"
  local gene_json="$run_dir/post_interpret/json/$report_id.$GENE_ID.gene_support.$COMPLETE_RULE.json"
  local audit_json="$run_dir/post_interpret/json/$report_id.$GENE_ID.gene_support.$COMPLETE_RULE.audit.json"
  local summary_json="$run_dir/reports/$report_id.final_summary.json"
  local seed_only_note_json="$run_dir/post_interpret/json/$report_id.seed_only.json"
  local lock_dir="$run_dir/.worker.lock"
  local -a gentle seed_args

  mkdir -p "$run_dir"/{logs,reports,checkpoints,post_interpret/{json,tsv,svg,fasta}}
  if ! mkdir "$lock_dir" 2>/dev/null; then
    die "Worker lock exists for $run: $lock_dir"
  fi
  trap_rm_rf_on_exit "$lock_dir"

  if [ -s "$summary_json" ]; then
    log "Summary already exists for $run: $summary_json"
    return 0
  fi
  [ ! -e "$state" ] || die "Refusing to overwrite existing state: $state"
  require_file "$read_fasta"
  gzip -t "$read_fasta"
  cp "$BASE_STATE" "$state"
  gentle=("$GENTLE_BIN" --state "$state" --progress-stderr)
  read -r -a seed_args < "$SEED_ARGS_FILE"

  sample_step "$run_dir" "[$run] interpret $GENE_ID as $report_id"
  /usr/bin/time -v "${gentle[@]}" rna-reads interpret "$SEQ_ID" "$SEED_FEATURE_ID" "$read_fasta" \
    --report-id "$report_id" \
    --report-mode "$REPORT_MODE" \
    --checkpoint-path "$run_dir/checkpoints/$report_id.checkpoint.json" \
    --checkpoint-every-reads "$CHECKPOINT_EVERY_READS" \
    --profile "$PROFILE" \
    --format "$INPUT_FORMAT" \
    --scope "$SCOPE" \
    --origin-mode "$ORIGIN_MODE" \
    "${seed_args[@]}" \
    --align-band-bp "$ALIGN_BAND_BP" \
    --align-min-identity "$ALIGN_MIN_IDENTITY" \
    --max-secondary-mappings "$MAX_SECONDARY_MAPPINGS" \
    > "$run_dir/reports/$report_id.interpret.command.json" \
    2> "$run_dir/logs/$report_id.interpret.stderr.log"

  if [ "$SCREEN_PHASE" = "seed_only" ]; then
    sample_step "$run_dir" "[$run] export seed-only report for $GENE_ID"
    "${gentle[@]}" rna-reads show-report "$report_id" \
      > "$show_json" \
      2> "$run_dir/logs/$report_id.show_report.stderr.log"
    printf '{}\n' > "$gene_json"
    printf '{"rows":[]}\n' > "$audit_json"
    jq -n \
      --arg mode "$SCREEN_PHASE" \
      --arg selection "$ALIGN_SELECTION" \
      --arg reason "cohort-scale screen stops after seed interpretation; run with --with-alignment for phase-2 mapping" \
      '{analysis_phase: $mode, align_selection_if_enabled: $selection, note: $reason}' \
      > "$seed_only_note_json"

    sample_step "$run_dir" "[$run] write seed-only final summary"
    write_final_summary "$run" "$sample_id" "$sample_name" "$note" \
      "$report_id" "$show_json" "$gene_json" "$audit_json" "$summary_json"
    printf '0\n' > "$run_dir/worker.exit"
    sample_step "$run_dir" "[$run] complete"
    return 0
  fi

  sample_step "$run_dir" "[$run] align retained rows for $GENE_ID"
  "${gentle[@]}" rna-reads align-report "$report_id" \
    --selection "$ALIGN_SELECTION" \
    --align-band-bp "$ALIGN_BAND_BP" \
    --align-min-identity "$ALIGN_MIN_IDENTITY" \
    --max-secondary-mappings "$MAX_SECONDARY_MAPPINGS" \
    > "$run_dir/post_interpret/json/$report_id.align_report.json" \
    2> "$run_dir/logs/$report_id.align_report.stderr.log"

  sample_step "$run_dir" "[$run] export aligned report and gene-support audit"
  "${gentle[@]}" rna-reads show-report "$report_id" \
    > "$show_json" \
    2> "$run_dir/logs/$report_id.show_report.stderr.log"

  "${gentle[@]}" rna-reads summarize-gene-support "$report_id" \
    --gene "$GENE_ID" --complete-rule "$COMPLETE_RULE" \
    --output "$gene_json" \
    > "$run_dir/post_interpret/json/$report_id.$GENE_ID.gene_support.$COMPLETE_RULE.command.json" \
    2> "$run_dir/logs/$report_id.$GENE_ID.gene_support.$COMPLETE_RULE.stderr.log"

  "${gentle[@]}" rna-reads inspect-gene-support "$report_id" \
    --gene "$GENE_ID" --complete-rule "$COMPLETE_RULE" --cohort all \
    --output "$audit_json" \
    > "$run_dir/post_interpret/json/$report_id.$GENE_ID.gene_support.$COMPLETE_RULE.audit.command.json" \
    2> "$run_dir/logs/$report_id.$GENE_ID.gene_support.$COMPLETE_RULE.audit.stderr.log"

  "${gentle[@]}" rna-reads inspect-alignments "$report_id" \
    --selection aligned --limit 200 --sort score \
    > "$run_dir/post_interpret/json/$report_id.aligned.top200.json" \
    2> "$run_dir/logs/$report_id.aligned.top200.stderr.log"

  "${gentle[@]}" rna-reads export-sample-sheet \
    "$run_dir/post_interpret/tsv/$report_id.sample_sheet.$COMPLETE_RULE.tsv" \
    --report-id "$report_id" \
    --gene "$GENE_ID" \
    --complete-rule "$COMPLETE_RULE" \
    > "$run_dir/post_interpret/json/$report_id.sample_sheet.$COMPLETE_RULE.command.json" \
    2> "$run_dir/logs/$report_id.sample_sheet.$COMPLETE_RULE.stderr.log"

  "${gentle[@]}" rna-reads export-alignments-tsv "$report_id" \
    "$run_dir/post_interpret/tsv/$report_id.alignments.aligned.tsv" \
    --selection aligned --subset-spec "${GENE_SAFE}_pancreas_${run}" \
    > "$run_dir/post_interpret/json/$report_id.export_alignments.command.json" \
    2> "$run_dir/logs/$report_id.export_alignments.stderr.log"

  "${gentle[@]}" rna-reads export-paths-tsv "$report_id" \
    "$run_dir/post_interpret/tsv/$report_id.paths.aligned.tsv" \
    --selection aligned --subset-spec "${GENE_SAFE}_pancreas_${run}" \
    > "$run_dir/post_interpret/json/$report_id.export_paths.command.json" \
    2> "$run_dir/logs/$report_id.export_paths.stderr.log"

  "${gentle[@]}" rna-reads export-abundance-tsv "$report_id" \
    "$run_dir/post_interpret/tsv/$report_id.abundance.aligned.tsv" \
    --selection aligned --subset-spec "${GENE_SAFE}_pancreas_${run}" \
    > "$run_dir/post_interpret/json/$report_id.export_abundance.command.json" \
    2> "$run_dir/logs/$report_id.export_abundance.stderr.log"

  run_optional "target-quality export for $run" \
    "$run_dir/logs/$report_id.$GENE_ID.target_quality.stderr.log" \
    "${gentle[@]}" rna-reads export-target-quality "$report_id" \
      "$run_dir/post_interpret/svg/$report_id.$GENE_ID.target_quality.svg" \
      --gene "$GENE_ID" --complete-rule "$COMPLETE_RULE" \
      > "$run_dir/post_interpret/json/$report_id.$GENE_ID.target_quality.command.json"

  if [ "$CONCATEMER_LIMIT" -gt 0 ]; then
    sample_step "$run_dir" "[$run] optional concatemer audit, limit $CONCATEMER_LIMIT"
    "${gentle[@]}" rna-reads inspect-concatemers "$report_id" \
      --selection aligned --limit "$CONCATEMER_LIMIT" \
      > "$run_dir/post_interpret/json/$report_id.concatemers.limit_${CONCATEMER_LIMIT}.json" \
      2> "$run_dir/logs/$report_id.concatemers.limit_${CONCATEMER_LIMIT}.stderr.log"
  fi

  sample_step "$run_dir" "[$run] write final summary"
  write_final_summary "$run" "$sample_id" "$sample_name" "$note" \
    "$report_id" "$show_json" "$gene_json" "$audit_json" "$summary_json"
  printf '0\n' > "$run_dir/worker.exit"
  sample_step "$run_dir" "[$run] complete"
}

run_one_sample() {
  local run="$1"
  local sample_id="$2"
  local sample_name="$3"
  local note="$4"
  local read_fasta="$5"
  local report_id="$6"
  local run_dir="$RUN_ROOT/runs/$run"
  local rc
  mkdir -p "$run_dir/logs"
  set +e
  (
    # Backgrounded shell functions run in a subshell. Reset inherited EXIT
    # traps so worker shutdown cannot touch the parent run lock, then restore
    # fail-fast behavior inside the actual sample pipeline.
    trap - EXIT
    set -euo pipefail
    run_one_sample_impl "$run" "$sample_id" "$sample_name" "$note" "$read_fasta" "$report_id"
  ) > "$run_dir/logs/$report_id.worker.stdout.log" \
    2> "$run_dir/logs/$report_id.worker.stderr.log"
  rc=$?
  set -e
  printf '%s\n' "$rc" > "$run_dir/worker.exit"
  if [ "$rc" -ne 0 ]; then
    printf '%s\n' "failed" > "$run_dir/CURRENT_STEP.txt"
  fi
  return "$rc"
}

print_failed_worker_summaries() {
  local root="$1"
  local run_dir run_name exit_code step stderr_log
  for run_dir in "$root"/runs/*; do
    [ -d "$run_dir" ] || continue
    run_name="$(basename "$run_dir")"
    exit_code="$(cat "$run_dir/worker.exit" 2>/dev/null || printf 'missing')"
    [ "$exit_code" != "0" ] || continue
    step="$(cat "$run_dir/CURRENT_STEP.txt" 2>/dev/null || printf 'unknown step')"
    log "FAILED worker: run=$run_name exit=$exit_code step=$step"
    for stderr_log in "$run_dir"/logs/*.worker.stderr.log "$run_dir"/logs/*.interpret.stderr.log "$run_dir"/logs/*.align_report.stderr.log; do
      [ -s "$stderr_log" ] || continue
      log "stderr tail: $stderr_log"
      tail -n 30 "$stderr_log" >&2 || true
    done
  done
}

active_jobs_count() {
  jobs -pr | wc -l | tr -d ' '
}

write_monitor_script() {
  local monitor="$RUN_ROOT/monitor.sh"
  cat > "$monitor" <<EOF
#!/usr/bin/env bash
set -euo pipefail
RUN_ROOT="$RUN_ROOT"
SCRIPT="$SCRIPT_PATH"
INTERVAL="\${1:-120}"

while true; do
  date
  "\$SCRIPT" status "\$RUN_ROOT" | column -t -s \$'\\t' || "\$SCRIPT" status "\$RUN_ROOT"
  echo
  echo "Active GENtle RNA-read processes:"
  ps -axo pid,etime,pcpu,pmem,command \\
    | grep -E 'gentle_cli .*rna-reads|pancreas_gene_rna_screen' \\
    | grep -v grep || true
  echo
  echo "Recent checkpoints/summaries:"
  find "\$RUN_ROOT/runs" -maxdepth 3 -type f \\
    \\( -name '*.checkpoint.json' -o -name '*.final_summary.json' -o -name 'worker.exit' \\) \\
    -printf '%TY-%Tm-%Td %TH:%TM %s %p\\n' 2>/dev/null \\
    | sort | tail -n 30 || true
  echo
  sleep "\$INTERVAL"
done
EOF
  chmod +x "$monitor"
  printf '%s\n' "$monitor"
}

print_monitor_hint() {
  local monitor="$RUN_ROOT/monitor.sh"
  cat <<EOF

Monitor this run from another shell/screen with:
  $monitor 120

One-shot status table:
  $SCRIPT_PATH status "$RUN_ROOT" | column -t -s \$'\\t'

Run root:
  $RUN_ROOT
EOF
}

status_from_manifest_row() {
  local run="$1"
  local sample_id="$2"
  local sample_name="$3"
  local note="$4"
  local read_fasta="$5"
  local report_id="$6"
  local run_dir="$RUN_ROOT/runs/$run"
  local summary="$run_dir/reports/$report_id.final_summary.json"
  local checkpoint="$run_dir/checkpoints/$report_id.checkpoint.json"
  local exit_file="$run_dir/worker.exit"
  local step_file="$run_dir/CURRENT_STEP.txt"
  local status reads seed aligned strict_accepted step exit_code
  status="pending"
  reads=""
  seed=""
  aligned=""
  strict_accepted=""
  step=""
  exit_code=""
  if [ -s "$summary" ]; then
    status="complete"
    reads="$(jq -r '.read_count_total // ""' "$summary")"
    seed="$(jq -r '.strict_seed_passed_reads // ""' "$summary")"
    aligned="$(jq -r '.retained_aligned_rows // ""' "$summary")"
    strict_accepted="$(jq -r '.strict_seed_accepted_target_rows // ""' "$summary")"
  elif [ -s "$checkpoint" ]; then
    status="running"
    reads="$(jq -r '.reads_processed // ""' "$checkpoint")"
    seed="$(jq -r '.read_count_seed_passed // ""' "$checkpoint")"
    aligned="$(jq -r '.read_count_aligned // ""' "$checkpoint")"
  fi
  if [ -s "$exit_file" ]; then
    exit_code="$(cat "$exit_file")"
    if [ "$exit_code" != "0" ] && [ "$status" != "complete" ]; then
      status="failed"
    fi
  fi
  if [ -s "$step_file" ]; then
    step="$(cat "$step_file")"
  fi
  printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
    "$run" "$sample_id" "$sample_name" \
    "$([ -s "$read_fasta" ] && echo yes || echo no)" \
    "$status" "$reads" "$seed" "$aligned" "$strict_accepted" "$exit_code" "$step"
}

load_existing_run_env() {
  RUN_ROOT="$1"
  require_file "$RUN_ROOT/run.env"
  # shellcheck disable=SC1090,SC1091
  source "$RUN_ROOT/run.env"
  # Run roots created before SCREEN_PHASE existed always ran phase-2 alignment.
  if ! grep -q '^export SCREEN_PHASE=' "$RUN_ROOT/run.env"; then
    SCREEN_PHASE="with_alignment"
  fi
  if ! grep -q '^export ALIGN_SELECTION=' "$RUN_ROOT/run.env"; then
    ALIGN_SELECTION="all"
  fi
}

status_run_root() {
  load_existing_run_env "$1"
  local sample_manifest="$RUN_ROOT/manifests/${GENE_SAFE}_pancreas_inputs.tsv"
  require_file "$sample_manifest"
  printf 'run_accession\tsample_id\tsample_name\tfasta_present\tanalysis_status\treads_or_processed\tstrict_seed_passed\tretained_aligned\tstrict_seed_accepted_target\texit_code\tcurrent_step\n'
  while IFS=$'\t' read -r run sample_id sample_name note read_fasta report_id rest; do
    if [ -z "${run:-}" ] || [ "$run" = "run_accession" ]; then
      continue
    fi
    status_from_manifest_row "$run" "$sample_id" "$sample_name" "$note" "$read_fasta" "$report_id"
  done < "$sample_manifest"
}

refresh_missing_read_length_tail_final_summaries() {
  local sample_manifest="$RUN_ROOT/manifests/${GENE_SAFE}_pancreas_inputs.tsv"
  local run sample_id sample_name note read_fasta report_id rest
  local run_dir show_json gene_json audit_json summary_json
  require_file "$sample_manifest"

  while IFS=$'\t' read -r run sample_id sample_name note read_fasta report_id rest; do
    if [ -z "${run:-}" ] || [ "$run" = "run_accession" ]; then
      continue
    fi
    run_dir="$RUN_ROOT/runs/$run"
    show_json="$run_dir/post_interpret/json/$report_id.show_report.json"
    gene_json="$run_dir/post_interpret/json/$report_id.$GENE_ID.gene_support.$COMPLETE_RULE.json"
    audit_json="$run_dir/post_interpret/json/$report_id.$GENE_ID.gene_support.$COMPLETE_RULE.audit.json"
    summary_json="$run_dir/reports/$report_id.final_summary.json"
    [ -s "$show_json" ] || continue
    [ -s "$gene_json" ] || continue
    [ -s "$audit_json" ] || continue
    if [ ! -s "$summary_json" ] \
      || ! jq -e '.read_length_summaries.all_reads.quantiles_bp | has("q90") and has("q95") and has("q99")' "$summary_json" >/dev/null 2>&1
    then
      log "Refresh final summary with q90/q95/q99 read-length quantiles: $summary_json"
      write_final_summary "$run" "$sample_id" "$sample_name" "$note" \
        "$report_id" "$show_json" "$gene_json" "$audit_json" "$summary_json"
    fi
  done < "$sample_manifest"
}

summarize_run_root() {
  load_existing_run_env "$1"
  local reports_dir="$RUN_ROOT/reports"
  local figures_dir="$RUN_ROOT/figures"
  local summaries_json="$reports_dir/${GENE_SAFE}_pancreas.final_summaries.json"
  local summary_tsv="$reports_dir/${GENE_SAFE}_pancreas.summary.tsv"
  local gene_screen_tsv="$reports_dir/${GENE_SAFE}_pancreas.gene_screen_summary.tsv"
  local figure_tsv="$figures_dir/${GENE_SAFE}_pancreas_figure_source.tsv"
  local summary_files="$RUN_ROOT/manifests/final_summary_files.txt"
  mkdir -p "$reports_dir" "$figures_dir"
  refresh_missing_read_length_tail_final_summaries

  find "$RUN_ROOT/runs" -path "*/reports/*.final_summary.json" | sort > "$summary_files"
  if [ ! -s "$summary_files" ]; then
    die "No final summaries found under $RUN_ROOT/runs"
  fi

  local -a summary_paths
  while IFS= read -r summary_path; do
    summary_paths+=("$summary_path")
  done < "$summary_files"
  jq -s '.' "${summary_paths[@]}" > "$summaries_json"

  {
    printf 'run_accession\tsample_id\tsample_name\treport_id\tread_count_total\tstrict_seed_passed_reads\tstrict_seed_accepted_target_rows\tstrict_seed_complete_target_rows\tstrict_seed_full_length_near_rows\tretained_report_rows\tretained_aligned_rows\taccepted_target_count\taccepted_target_fraction_total\taccepted_target_fraction_aligned\tfragment_count\tcomplete_near_count\tcomplete_strict_count\tcomplete_exact_count\tall_mean_bp\tall_q0_bp\tall_q25_bp\tall_q50_bp\tall_q75_bp\tall_q90_bp\tall_q95_bp\tall_q99_bp\tall_q100_bp\tstrict_seed_mean_bp\tstrict_seed_q90_bp\tstrict_seed_q95_bp\tstrict_seed_q99_bp\tstrict_seed_q100_bp\taccepted_mean_bp\taccepted_max_bp\n'
    jq -r '
      .[]
      | [
          .run_accession,
          .sample_id,
          .sample_name,
          .report_id,
          .read_count_total,
          .strict_seed_passed_reads,
          .strict_seed_accepted_target_rows,
          .strict_seed_complete_target_rows,
          .strict_seed_full_length_near_rows,
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
          .read_length_summaries.all_reads.quantiles_bp.q0,
          .read_length_summaries.all_reads.quantiles_bp.q25,
          .read_length_summaries.all_reads.quantiles_bp.q50,
          .read_length_summaries.all_reads.quantiles_bp.q75,
          .read_length_summaries.all_reads.quantiles_bp.q90,
          .read_length_summaries.all_reads.quantiles_bp.q95,
          .read_length_summaries.all_reads.quantiles_bp.q99,
          .read_length_summaries.all_reads.quantiles_bp.q100,
          .read_length_summaries.strict_seed_passed.mean_bp,
          .read_length_summaries.strict_seed_passed.quantiles_bp.q90,
          .read_length_summaries.strict_seed_passed.quantiles_bp.q95,
          .read_length_summaries.strict_seed_passed.quantiles_bp.q99,
          .read_length_summaries.strict_seed_passed.quantiles_bp.q100,
          .gene_support.accepted_target_read_lengths.mean_length_bp,
          .gene_support.accepted_target_read_lengths.max_length_bp
        ]
      | @tsv
    ' "$summaries_json"
  } > "$summary_tsv"

  {
    printf 'run_accession\tsample_id\tsample_name\ttotal_reads\tall_q0_bp\tall_q25_bp\tall_q50_bp\tall_q75_bp\tall_q90_bp\tall_q95_bp\tall_q99_bp\tall_q100_bp\tall_mean_bp\tstrict_seed_passed_reads\tstrict_seed_accepted_target_reads\tstrict_seed_passed_q90_read_bp\tstrict_seed_passed_q95_read_bp\tstrict_seed_passed_q99_read_bp\tstrict_seed_passed_max_read_bp\tstrict_seed_passed_mean_read_bp\taccepted_target_reads\taccepted_target_per_million\taccepted_target_max_read_bp\taccepted_target_mean_read_bp\n'
    jq -r '
      .[]
      | [
          .run_accession,
          .sample_id,
          .sample_name,
          .read_count_total,
          .read_length_summaries.all_reads.quantiles_bp.q0,
          .read_length_summaries.all_reads.quantiles_bp.q25,
          .read_length_summaries.all_reads.quantiles_bp.q50,
          .read_length_summaries.all_reads.quantiles_bp.q75,
          .read_length_summaries.all_reads.quantiles_bp.q90,
          .read_length_summaries.all_reads.quantiles_bp.q95,
          .read_length_summaries.all_reads.quantiles_bp.q99,
          .read_length_summaries.all_reads.quantiles_bp.q100,
          .read_length_summaries.all_reads.mean_bp,
          .strict_seed_passed_reads,
          .strict_seed_accepted_target_rows,
          .read_length_summaries.strict_seed_passed.quantiles_bp.q90,
          .read_length_summaries.strict_seed_passed.quantiles_bp.q95,
          .read_length_summaries.strict_seed_passed.quantiles_bp.q99,
          .read_length_summaries.strict_seed_passed.quantiles_bp.q100,
          .read_length_summaries.strict_seed_passed.mean_bp,
          .gene_support.accepted_target_count,
          (if .read_count_total > 0 and (.gene_support.accepted_target_count != null)
           then (.gene_support.accepted_target_count * 1000000 / .read_count_total)
           else null end),
          .gene_support.accepted_target_read_lengths.max_length_bp,
          .gene_support.accepted_target_read_lengths.mean_length_bp
        ]
      | @tsv
    ' "$summaries_json"
  } > "$figure_tsv"

  {
    printf 'schema\tgene\trun_accession\tsample_id\tsample_name\tsource_kind\tsource_path\tanalysis_phase\treport_id\tseq_id\tseed_feature_id\tinput_path\ttotal_reads\tall_q0_bp\tall_q25_bp\tall_q50_bp\tall_q75_bp\tall_q90_bp\tall_q95_bp\tall_q99_bp\tall_q100_bp\tall_mean_bp\tseed_passed_reads\tseed_passed_per_million\tseed_passed_q90_bp\tseed_passed_q95_bp\tseed_passed_q99_bp\tseed_passed_max_bp\tseed_passed_mean_bp\taccepted_target_count\taccepted_target_per_million\taccepted_target_max_bp\taccepted_target_mean_bp\talign_selection\n'
    jq -r \
      --arg schema 'gentle.rna_read_gene_screen_summary.v1' \
      --arg gene "$GENE_ID" \
      --arg source_path "$summaries_json" \
      --arg seq_id "$SEQ_ID" \
      --arg seed_feature_id "$SEED_FEATURE_ID" '
      .[]
      | [
          $schema,
          $gene,
          .run_accession,
          .sample_id,
          .sample_name,
          "pancreas_gene_rna_screen",
          $source_path,
          .analysis_phase,
          .report_id,
          $seq_id,
          $seed_feature_id,
          null,
          .read_count_total,
          .read_length_summaries.all_reads.quantiles_bp.q0,
          .read_length_summaries.all_reads.quantiles_bp.q25,
          .read_length_summaries.all_reads.quantiles_bp.q50,
          .read_length_summaries.all_reads.quantiles_bp.q75,
          .read_length_summaries.all_reads.quantiles_bp.q90,
          .read_length_summaries.all_reads.quantiles_bp.q95,
          .read_length_summaries.all_reads.quantiles_bp.q99,
          .read_length_summaries.all_reads.quantiles_bp.q100,
          .read_length_summaries.all_reads.mean_bp,
          .strict_seed_passed_reads,
          (if (.read_count_total // 0) > 0
           then ((.strict_seed_passed_reads // 0) * 1000000 / .read_count_total)
           else null end),
          .read_length_summaries.strict_seed_passed.quantiles_bp.q90,
          .read_length_summaries.strict_seed_passed.quantiles_bp.q95,
          .read_length_summaries.strict_seed_passed.quantiles_bp.q99,
          .read_length_summaries.strict_seed_passed.quantiles_bp.q100,
          .read_length_summaries.strict_seed_passed.mean_bp,
          .gene_support.accepted_target_count,
          (if (.read_count_total // 0) > 0 and (.gene_support.accepted_target_count != null)
           then (.gene_support.accepted_target_count * 1000000 / .read_count_total)
           else null end),
          .gene_support.accepted_target_read_lengths.max_length_bp,
          .gene_support.accepted_target_read_lengths.mean_length_bp,
          .align_selection
        ]
      | @tsv
    ' "$summaries_json"
  } > "$gene_screen_tsv"

  log "Wrote summary JSON: $summaries_json"
  log "Wrote summary TSV:  $summary_tsv"
  log "Wrote canonical gene-screen TSV: $gene_screen_tsv"
  log "Wrote figure TSV:   $figure_tsv"
}

run_screen() {
  require_common_tools
  [ ! -e "$RUN_ROOT" ] || die "Refusing to reuse existing run root: $RUN_ROOT"
  mkdir -p "$RUN_ROOT"/{logs,manifests,reports,resources,runs,tmp,figures}
  local run_lock="$RUN_ROOT/.run.lock"
  mkdir "$run_lock" || die "Could not acquire run lock: $run_lock"
  trap_rm_rf_on_exit "$run_lock"

  log "Run root: $RUN_ROOT"
  resolve_ncbi_gene_locus
  load_base_state
  resolve_seed_feature
  run_preflight
  build_sample_manifest
  write_run_env
  write_monitor_script >/dev/null
  print_monitor_hint
  if [ "$SCREEN_PHASE" = "seed_only" ] && [ "$CONCATEMER_LIMIT" -gt 0 ]; then
    log "NOTE: --concatemer-limit is ignored in seed-only mode; use --with-alignment for concatemer audits"
  fi

  local sample_manifest="$RUN_ROOT/manifests/${GENE_SAFE}_pancreas_inputs.tsv"
  local pids=""
  local launched=0
  local run sample_id sample_name note read_fasta report_id rest

  log "Start $GENE_ID screen with $JOBS parallel worker(s); phase=$SCREEN_PHASE; report_mode=$REPORT_MODE; align_selection=$ALIGN_SELECTION"
  while IFS=$'\t' read -r run sample_id sample_name note read_fasta report_id rest; do
    if [ -z "${run:-}" ] || [ "$run" = "run_accession" ]; then
      continue
    fi
    while [ "$(active_jobs_count)" -ge "$JOBS" ]; do
      sleep 5
    done
    log "Launch $run ($sample_name) -> $report_id"
    run_one_sample "$run" "$sample_id" "$sample_name" "$note" "$read_fasta" "$report_id" &
    pids="$pids $!"
    launched=$((launched + 1))
  done < "$sample_manifest"

  local failures=0
  local pid
  for pid in $pids; do
    if ! wait "$pid"; then
      failures=$((failures + 1))
    fi
  done

  log "Workers finished: launched=$launched failures=$failures"
  summarize_run_root "$RUN_ROOT" || true
  if [ "$failures" -ne 0 ]; then
    print_failed_worker_summaries "$RUN_ROOT"
    die "$failures sample worker(s) failed; inspect $RUN_ROOT/runs/*/logs"
  fi
  log "Completed $GENE_ID screen"
}

monitor_run_root() {
  local interval="${2:-120}"
  while true; do
    date
    status_run_root "$1" | column -t -s $'\t' || status_run_root "$1"
    echo
    echo "Active GENtle RNA-read processes:"
    # shellcheck disable=SC2009
    ps -axo pid,etime,pcpu,pmem,command \
      | grep -E 'gentle_cli .*rna-reads|pancreas_gene_rna_screen' \
      | grep -v grep || true
    echo
    sleep "$interval"
  done
}

main() {
  local command="${1:-}"
  shift || true
  case "$command" in
    run)
      parse_run_args "$@"
      run_screen
      ;;
    group-plan)
      parse_group_plan_args "$@"
      ;;
    status)
      [ $# -eq 1 ] || die "status requires RUN_ROOT"
      status_run_root "$1"
      ;;
    summarize)
      [ $# -eq 1 ] || die "summarize requires RUN_ROOT"
      summarize_run_root "$1"
      ;;
    monitor)
      [ $# -ge 1 ] || die "monitor requires RUN_ROOT [SECONDS]"
      monitor_run_root "$1" "${2:-120}"
      ;;
    print-monitor)
      [ $# -eq 1 ] || die "print-monitor requires RUN_ROOT"
      load_existing_run_env "$1"
      print_monitor_hint
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
