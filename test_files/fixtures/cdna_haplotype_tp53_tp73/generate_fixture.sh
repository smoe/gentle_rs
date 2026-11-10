#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
OUT_ROOT=${1:-$SCRIPT_DIR}
DATA_DIR="$OUT_ROOT/data"

STORAGE_ROOT=${GENTLE_CDNA_STORAGE_ROOT:-/mnt/storage-box-1-gentle/Glen/IEGT/nanopore}
COHORT_ROOT="$STORAGE_ROOT/analysis/gene_counts_v1/results/cohort_v1_1"
REFERENCE_ROOT="$STORAGE_ROOT/reference/ensembl-116-grch38"
REFERENCE="$REFERENCE_ROOT/Homo_sapiens.GRCh38.dna.primary_assembly.fa"
ANNOTATION="$REFERENCE_ROOT/Homo_sapiens.GRCh38.116.gtf"
SAMPLE_SHEET="$STORAGE_ROOT/metadata/samples.tsv"

TP53_GTF=${GENTLE_TP53_GTF:-/home/clawbio/.openclaw/workspace/projects/long-read-sequencing-presentation/analysis/isoforms-enriched/references/TP53.ensembl116.gtf}
TP73_GTF=${GENTLE_TP73_GTF:-/home/clawbio/.openclaw/workspace/projects/long-read-sequencing-presentation/analysis/isoforms-enriched/references/TP73.ensembl116.gtf}

for tool in samtools gzip sha256sum awk sort stat; do
  command -v "$tool" >/dev/null || {
    echo "missing required tool: $tool" >&2
    exit 2
  }
done

for path in "$REFERENCE" "$REFERENCE.fai" "$ANNOTATION" "$SAMPLE_SHEET" "$TP53_GTF" "$TP73_GTF"; do
  test -r "$path" || {
    echo "missing or unreadable input: $path" >&2
    exit 2
  }
done

if test -d "$DATA_DIR" && find "$DATA_DIR" -mindepth 1 -print -quit | grep -q .; then
  echo "refusing to overwrite non-empty generated tree: $DATA_DIR" >&2
  exit 2
fi

mkdir -p "$DATA_DIR/annotation" "$DATA_DIR/reference" "$DATA_DIR/reads"

printf '%s\n' \
  $'gene\tgene_id\tchromosome\tstart_1based\tend_1based\tstrand\tassembly\tannotation_release' \
  $'TP73\tENSG00000078900.18\t1\t3652307\t3736201\t+\tGRCh38\tEnsembl 116' \
  $'TP53\tENSG00000141510.21\t17\t7661779\t7687546\t-\tGRCh38\tEnsembl 116' \
  > "$DATA_DIR/loci.tsv"

printf '%s\n' \
  $'barcode\tsample_id\tcell_line\tlibrary_content\tdescription\tcell_line_status' \
  $'barcode15\tsample1\tSK-MEL-28\tcDNA_regular\tcDNA regular; no TP53/TP73 target primer\tconfirmed' \
  $'barcode16\tsample2\tWM793\tcDNA_regular\tcDNA regular; no TP53/TP73 target primer\tconfirmed' \
  > "$DATA_DIR/samples.tsv"

{
  printf '# TP53 and TP73 records copied from the Ensembl 116 per-gene extracts.\n'
  cat "$TP53_GTF" "$TP73_GTF"
} > "$DATA_DIR/annotation/ensembl116.tp53_tp73.gtf"

samtools faidx "$REFERENCE" '1:3652307-3736201' \
  > "$DATA_DIR/reference/ensembl116.tp53_tp73.genomic_regions.fa"
samtools faidx "$REFERENCE" '17:7661779-7687546' \
  >> "$DATA_DIR/reference/ensembl116.tp53_tp73.genomic_regions.fa"

printf '%s\n' \
  $'barcode\tcell_line\tgene\tread_id\tflag\tchromosome\tstart_1based\tmapq\tcigar\tquery_length\tmean_base_quality' \
  > "$DATA_DIR/read_manifest.tsv"

for barcode in barcode15 barcode16; do
  case "$barcode" in
    barcode15) cell_line=SK-MEL-28 ;;
    barcode16) cell_line=WM793 ;;
  esac

  cram="$COHORT_ROOT/$barcode/$barcode.sorted.cram"
  crai="$cram.crai"
  receipt="$COHORT_ROOT/$barcode/receipt.json"
  test -r "$cram" && test -s "$crai" && test -r "$receipt" || {
    echo "missing CRAM, CRAI, or receipt for $barcode" >&2
    exit 2
  }
  samtools quickcheck -v "$cram"

  for gene_spec in 'TP73:1:3652307-3736201' 'TP53:17:7661779-7687546'; do
    IFS=: read -r gene chromosome span <<< "$gene_spec"
    bam="$DATA_DIR/reads/${barcode}.${gene}.primary.bam"
    fastq="$DATA_DIR/reads/${barcode}.${gene}.primary.fastq.gz"

    # Do not append a path-bearing @PG record: the source receipt already
    # binds the extraction tools, while --no-PG keeps regenerated BAM bytes
    # independent of the chosen output directory.
    samtools view --no-PG -b -F 2308 -T "$REFERENCE" "$cram" "$chromosome:$span" -o "$bam"
    samtools index "$bam"
    samtools quickcheck -v "$bam"
    samtools fastq -n -0 - "$bam" 2>/dev/null | gzip -n > "$fastq"

    samtools view "$bam" | awk -v barcode="$barcode" -v cell_line="$cell_line" -v gene="$gene" '
      BEGIN { OFS="\t" }
      {
        seq=$10; qual=$11; qlen=(seq=="*" ? 0 : length(seq));
        if (qual=="*" || qlen==0) meanq="NA";
        else {
          sum=0;
          for (i=1; i<=length(qual); i++) sum += index("!\"#$%&\047()*+,-./0123456789:;<=>?@ABCDEFGHIJKLMNOPQRSTUVWXYZ[\\]^_`abcdefghijklmnopqrstuvwxyz{|}~", substr(qual,i,1))-1;
          meanq=sprintf("%.3f", sum/length(qual));
        }
        print barcode,cell_line,gene,$1,$2,$3,$4,$5,$6,qlen,meanq;
      }
    ' >> "$DATA_DIR/read_manifest.tsv"

    gzip -dc "$fastq" | awk '
      NR%4==1 { records++ }
      NR%4==2 { seq=$0 }
      NR%4==0 {
        if (length(seq) != length($0)) {
          print "FASTQ sequence/quality length mismatch" > "/dev/stderr";
          exit 3;
        }
      }
      END {
        if (NR%4 != 0 || records == 0) exit 4;
      }
    '
  done
done

{
  printf 'role\tpath\tsize_bytes\tsha256_or_bound_digest\tnote\n'
  for barcode in barcode15 barcode16; do
    cram="$COHORT_ROOT/$barcode/$barcode.sorted.cram"
    receipt="$COHORT_ROOT/$barcode/receipt.json"
    case "$barcode" in
      barcode15) cram_sha256=b4f88f27d7cda404f07da36e95e4f106292966ed2d9eed9441b605ce5b902a9d ;;
      barcode16) cram_sha256=60fa43d8ba99a2b3132d4c581895a35292115ebe12d2c48795cb28827dab926e ;;
    esac
    printf 'source_cram\t%s\t%s\t%s\t%s\n' \
      "$cram" "$(stat -c %s "$cram")" "$cram_sha256" \
      "indexed archived full mapper output; full SHA-256 and samtools quickcheck verified 2026-11-10"
    printf 'cohort_receipt\t%s\t%s\t%s\t%s\n' \
      "$receipt" "$(stat -c %s "$receipt")" "$(sha256sum "$receipt" | awk '{print $1}')" \
      "source FASTQ hashes, tool versions, mapping command, and reference digests"
  done
  printf 'reference_fasta\t%s\t%s\t%s\t%s\n' \
    "$REFERENCE" "$(stat -c %s "$REFERENCE")" \
    '1e74081a49ceb9739cc14c812fbb8b3db978eb80ba8e5350beb80d8ad8dfef3b' \
    'digest from verified reference SHA256SUMS and cohort receipts'
  printf 'annotation_gtf\t%s\t%s\t%s\t%s\n' \
    "$ANNOTATION" "$(stat -c %s "$ANNOTATION")" \
    '474626b424caf2f2f512c0038d70d985832198551f6f2818b9528db9d1f4d9e7' \
    'digest from verified reference SHA256SUMS and cohort receipts'
  printf 'sample_sheet\t%s\t%s\t%s\t%s\n' \
    "$SAMPLE_SHEET" "$(stat -c %s "$SAMPLE_SHEET")" \
    "$(sha256sum "$SAMPLE_SHEET" | awk '{print $1}')" \
    'canonical corrected sample identities'
} > "$DATA_DIR/source_bindings.tsv"

(
  cd "$DATA_DIR"
  find . -type f ! -name SHA256SUMS -print0 \
    | LC_ALL=C sort -z \
    | xargs -0 sha256sum > SHA256SUMS
  sha256sum -c SHA256SUMS
)

echo "fixture generated at $DATA_DIR"
