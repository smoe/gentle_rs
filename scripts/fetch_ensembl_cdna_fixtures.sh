#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"

usage() {
  cat <<'EOF'
Fetch Ensembl cDNA transcript FASTA fixtures through GENtle.

Usage:
  scripts/fetch_ensembl_cdna_fixtures.sh [OPTIONS] [GENE ...]

Default genes:
  E2F1 E2F2 E2F3 E2F4 E2F5 E2F6 E2F7 E2F8

Options:
  --gene GENE            Add one gene symbol. Repeatable.
  --genes CSV            Add comma-separated gene symbols.
  --species NAME         Ensembl REST species name (default: homo_sapiens).
  --species-label LABEL  Filename label for species (default: human).
  --out-dir DIR          Output directory (default: test_files/fixtures/mapping).
  --work-dir DIR         GENtle state/report workspace (default: target/ensembl_cdna_fixtures).
  --gentle-bin PATH      gentle_cli binary (default: GENTLE_BIN, release, then debug).
  --force                Overwrite existing FASTA/provenance files.
  --dry-run              Print planned GENtle commands without executing them.
  -h, --help             Show this help.

Each gene writes:
  <out-dir>/ensembl_<species-label>_<gene-token>_all.fasta
  <out-dir>/ensembl_<species-label>_<gene-token>_all.provenance.json

The script deliberately uses GENtle for Ensembl gene retrieval, locus import,
transcript derivation, and FASTA export. The sidecar JSON plus the aggregate
manifest record exact commands and report/state paths for fixture provenance.
EOF
}

default_gentle_bin() {
  if [ -n "${GENTLE_BIN:-}" ]; then
    printf '%s\n' "$GENTLE_BIN"
  elif [ -x "$REPO_ROOT/target/release/gentle_cli" ]; then
    printf '%s\n' "$REPO_ROOT/target/release/gentle_cli"
  elif [ -x "$REPO_ROOT/target/debug/gentle_cli" ]; then
    printf '%s\n' "$REPO_ROOT/target/debug/gentle_cli"
  else
    printf '%s\n' "$REPO_ROOT/target/release/gentle_cli"
  fi
}

csv_to_genes() {
  local csv="$1"
  local old_ifs="$IFS"
  IFS=','
  # shellcheck disable=SC2206
  local parts=($csv)
  IFS="$old_ifs"
  for part in "${parts[@]}"; do
    part="${part//[[:space:]]/}"
    if [ -n "$part" ]; then
      GENES+=("$part")
    fi
  done
}

require_tool() {
  local name="$1"
  command -v "$name" >/dev/null 2>&1 || {
    echo "ERROR: required tool '$name' was not found in PATH" >&2
    exit 1
  }
}

quote_cmd() {
  printf '%q ' "$@"
}

run_cmd() {
  local log="$1"
  shift
  {
    printf '# '
    quote_cmd "$@"
    printf '\n'
  } > "$log.command.txt"
  "$@" > "$log"
}

sanitize_token() {
  printf '%s' "$1" | tr '[:upper:]' '[:lower:]' | tr -c 'a-z0-9_' '_'
}

write_manifest_header_if_needed() {
  local manifest="$1"
  if [ ! -s "$manifest" ]; then
    printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
      generated_at_utc gene species species_label entry_id locus_seq_id \
      output_fasta transcript_count gentle_version state_path fetch_report \
      derive_report provenance_json > "$manifest"
  fi
}

GENES=()
SPECIES="${SPECIES:-homo_sapiens}"
SPECIES_LABEL="${SPECIES_LABEL:-human}"
OUT_DIR="$REPO_ROOT/test_files/fixtures/mapping"
WORK_DIR="$REPO_ROOT/target/ensembl_cdna_fixtures"
GENTLE_BIN_SELECTED="$(default_gentle_bin)"
FORCE=0
DRY_RUN=0

while [ "$#" -gt 0 ]; do
  case "$1" in
    -h|--help)
      usage
      exit 0
      ;;
    --gene)
      if [ "$#" -lt 2 ]; then
        echo "ERROR: Missing GENE after --gene" >&2
        exit 1
      fi
      GENES+=("$2")
      shift 2
      ;;
    --genes)
      if [ "$#" -lt 2 ]; then
        echo "ERROR: Missing CSV after --genes" >&2
        exit 1
      fi
      csv_to_genes "$2"
      shift 2
      ;;
    --species)
      if [ "$#" -lt 2 ]; then
        echo "ERROR: Missing NAME after --species" >&2
        exit 1
      fi
      SPECIES="$2"
      shift 2
      ;;
    --species-label)
      if [ "$#" -lt 2 ]; then
        echo "ERROR: Missing LABEL after --species-label" >&2
        exit 1
      fi
      SPECIES_LABEL="$2"
      shift 2
      ;;
    --out-dir)
      if [ "$#" -lt 2 ]; then
        echo "ERROR: Missing DIR after --out-dir" >&2
        exit 1
      fi
      OUT_DIR="$2"
      shift 2
      ;;
    --work-dir)
      if [ "$#" -lt 2 ]; then
        echo "ERROR: Missing DIR after --work-dir" >&2
        exit 1
      fi
      WORK_DIR="$2"
      shift 2
      ;;
    --gentle-bin)
      if [ "$#" -lt 2 ]; then
        echo "ERROR: Missing PATH after --gentle-bin" >&2
        exit 1
      fi
      GENTLE_BIN_SELECTED="$2"
      shift 2
      ;;
    --force)
      FORCE=1
      shift
      ;;
    --dry-run)
      DRY_RUN=1
      shift
      ;;
    --)
      shift
      while [ "$#" -gt 0 ]; do
        GENES+=("$1")
        shift
      done
      ;;
    -*)
      echo "ERROR: Unknown option: $1" >&2
      exit 1
      ;;
    *)
      GENES+=("$1")
      shift
      ;;
  esac
done

if [ "${#GENES[@]}" -eq 0 ]; then
  GENES=(E2F1 E2F2 E2F3 E2F4 E2F5 E2F6 E2F7 E2F8)
fi

if [ "$DRY_RUN" = "1" ]; then
  echo "Dry run: no files will be written."
  echo "GENtle binary: $GENTLE_BIN_SELECTED"
  echo "Species: $SPECIES ($SPECIES_LABEL)"
  echo "Output directory: $OUT_DIR"
  echo "Work directory: $WORK_DIR"
  for gene in "${GENES[@]}"; do
    gene_token="$(sanitize_token "$gene")"
    species_token="$(sanitize_token "$SPECIES_LABEL")"
    entry_id="ensembl_${species_token}_${gene_token}"
    locus_seq_id="${entry_id}_locus"
    out_fasta="$OUT_DIR/ensembl_${species_token}_${gene_token}_all.fasta"
    state_path="$WORK_DIR/$entry_id/${entry_id}.gentle.json"
    echo
    echo "Gene: $gene -> $out_fasta"
    quote_cmd "$GENTLE_BIN_SELECTED" --state "$state_path" ensembl-gene fetch "$gene" --species "$SPECIES" --entry-id "$entry_id"; echo
    quote_cmd "$GENTLE_BIN_SELECTED" --state "$state_path" ensembl-gene import-sequence "$entry_id" --output-id "$locus_seq_id"; echo
    quote_cmd "$GENTLE_BIN_SELECTED" --state "$state_path" transcripts derive "$locus_seq_id" --output-prefix "${entry_id}__cdna"; echo
    echo "Then export each derived transcript sequence via GENtle SaveFile and concatenate."
  done
  exit 0
fi

require_tool jq

if [ ! -x "$GENTLE_BIN_SELECTED" ]; then
  echo "ERROR: gentle_cli binary is not executable: $GENTLE_BIN_SELECTED" >&2
  echo "Build it with: cargo build --release --bin gentle_cli" >&2
  exit 1
fi

mkdir -p "$OUT_DIR" "$WORK_DIR"
MANIFEST="$OUT_DIR/ensembl_cdna_fixture_manifest.tsv"
write_manifest_header_if_needed "$MANIFEST"

GENTLE_VERSION="$("$GENTLE_BIN_SELECTED" --version 2>/dev/null || true)"
GENERATED_AT_UTC="$(date -u '+%Y-%m-%dT%H:%M:%SZ')"

for gene in "${GENES[@]}"; do
  gene_token="$(sanitize_token "$gene")"
  species_token="$(sanitize_token "$SPECIES_LABEL")"
  entry_id="ensembl_${species_token}_${gene_token}"
  locus_seq_id="${entry_id}_locus"
  out_fasta="$OUT_DIR/ensembl_${species_token}_${gene_token}_all.fasta"
  provenance_json="$OUT_DIR/ensembl_${species_token}_${gene_token}_all.provenance.json"
  gene_work="$WORK_DIR/$entry_id"
  state_path="$gene_work/$entry_id.gentle.json"
  fetch_json="$gene_work/fetch_ensembl_gene.json"
  import_json="$gene_work/import_ensembl_gene_sequence.json"
  derive_json="$gene_work/derive_transcripts.json"
  parts_dir="$gene_work/fasta_parts"
  tmp_fasta="$out_fasta.tmp.$$"

  if [ "$FORCE" != "1" ] && { [ -e "$out_fasta" ] || [ -e "$provenance_json" ]; }; then
    echo "ERROR: Refusing to overwrite existing output for $gene: $out_fasta" >&2
    echo "Use --force to regenerate." >&2
    exit 1
  fi

  mkdir -p "$gene_work" "$parts_dir"
  rm -f "$tmp_fasta"
  rm -f "$parts_dir"/*.fa

  echo "[$(date '+%Y-%m-%d %H:%M:%S')] Fetch $gene via GENtle"
  run_cmd "$fetch_json" \
    "$GENTLE_BIN_SELECTED" --state "$state_path" \
    ensembl-gene fetch "$gene" --species "$SPECIES" --entry-id "$entry_id"

  echo "[$(date '+%Y-%m-%d %H:%M:%S')] Import $gene locus into GENtle state"
  run_cmd "$import_json" \
    "$GENTLE_BIN_SELECTED" --state "$state_path" \
    ensembl-gene import-sequence "$entry_id" --output-id "$locus_seq_id"

  echo "[$(date '+%Y-%m-%d %H:%M:%S')] Derive transcript cDNAs for $gene"
  run_cmd "$derive_json" \
    "$GENTLE_BIN_SELECTED" --state "$state_path" \
    transcripts derive "$locus_seq_id" --output-prefix "${entry_id}__cdna"

  mapfile -t derived_ids < <(jq -r '.transcripts[]?.seq_id' "$derive_json")
  if [ "${#derived_ids[@]}" -eq 0 ]; then
    echo "ERROR: GENtle did not derive any transcript sequences for $gene" >&2
    exit 1
  fi

  part_paths=()
  idx=0
  for derived_id in "${derived_ids[@]}"; do
    idx=$((idx + 1))
    part_path="$parts_dir/$(printf '%04d' "$idx").fa"
    save_json="$parts_dir/$(printf '%04d' "$idx").save.json"
    op_json="$(jq -n --arg seq_id "$derived_id" --arg path "$part_path" \
      '{"SaveFile":{"seq_id":$seq_id,"path":$path,"format":"Fasta"}}')"
    run_cmd "$save_json" "$GENTLE_BIN_SELECTED" --state "$state_path" op "$op_json"
    part_paths+=("$part_path")
  done

  cat "${part_paths[@]}" > "$tmp_fasta"
  mv "$tmp_fasta" "$out_fasta"

  transcript_ids_json="$(jq -c '.transcripts[]?.seq_id' "$derive_json" | jq -sc '.')"
  jq -n \
    --arg schema "gentle.ensembl_cdna_fixture_provenance.v1" \
    --arg generated_at_utc "$GENERATED_AT_UTC" \
    --arg gene "$gene" \
    --arg species "$SPECIES" \
    --arg species_label "$SPECIES_LABEL" \
    --arg entry_id "$entry_id" \
    --arg locus_seq_id "$locus_seq_id" \
    --arg output_fasta "$out_fasta" \
    --arg state_path "$state_path" \
    --arg gentle_bin "$GENTLE_BIN_SELECTED" \
    --arg gentle_version "$GENTLE_VERSION" \
    --arg fetch_report "$fetch_json" \
    --arg import_report "$import_json" \
    --arg derive_report "$derive_json" \
    --argjson transcript_seq_ids "$transcript_ids_json" \
    '{
      schema: $schema,
      generated_at_utc: $generated_at_utc,
      gene: $gene,
      species: $species,
      species_label: $species_label,
      entry_id: $entry_id,
      locus_seq_id: $locus_seq_id,
      output_fasta: $output_fasta,
      transcript_count: ($transcript_seq_ids | length),
      transcript_seq_ids: $transcript_seq_ids,
      gentle: {
        binary: $gentle_bin,
        version: $gentle_version,
        state_path: $state_path,
        commands: {
          fetch_gene: ["ensembl-gene","fetch",$gene,"--species",$species,"--entry-id",$entry_id],
          import_sequence: ["ensembl-gene","import-sequence",$entry_id,"--output-id",$locus_seq_id],
          derive_transcripts: ["transcripts","derive",$locus_seq_id,"--output-prefix",($entry_id + "__cdna")]
        },
        reports: {
          fetch_gene: $fetch_report,
          import_sequence: $import_report,
          derive_transcripts: $derive_report
        }
      }
    }' > "$provenance_json"

  printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
    "$GENERATED_AT_UTC" "$gene" "$SPECIES" "$SPECIES_LABEL" "$entry_id" \
    "$locus_seq_id" "$out_fasta" "${#derived_ids[@]}" "$GENTLE_VERSION" \
    "$state_path" "$fetch_json" "$derive_json" "$provenance_json" >> "$MANIFEST"

  echo "[$(date '+%Y-%m-%d %H:%M:%S')] Wrote $out_fasta (${#derived_ids[@]} transcript records)"
done

echo "Manifest: $MANIFEST"
