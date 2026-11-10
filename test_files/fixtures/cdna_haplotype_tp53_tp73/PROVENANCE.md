# Provenance and scientific boundary

## Origin

- Delivery: ONT run `26R266`, staged under the IEGT nanopore delivery.
- Biological mapping (confirmed metadata, 2026-11-02):
  - `barcode15` = SK-MEL-28 / `sample1` / regular cDNA.
  - `barcode16` = WM793 / `sample2` / regular cDNA.
- Neither library used a TP53/TP73 target primer.
- Alignment: minimap2 `2.27-r1193`, preset `splice`, against Ensembl 116
  GRCh38 primary assembly; archived as CRAM by samtools `1.24`.
- Annotation: Ensembl 116 GTF.

The canonical local sample sheet is:

`/mnt/storage-box-1-gentle/Glen/IEGT/nanopore/metadata/samples.tsv`

The source CRAMs are:

- `/mnt/storage-box-1-gentle/Glen/IEGT/nanopore/analysis/gene_counts_v1/results/cohort_v1_1/barcode15/barcode15.sorted.cram`
- `/mnt/storage-box-1-gentle/Glen/IEGT/nanopore/analysis/gene_counts_v1/results/cohort_v1_1/barcode16/barcode16.sorted.cram`

The exact reference and annotation are:

- `/mnt/storage-box-1-gentle/Glen/IEGT/nanopore/reference/ensembl-116-grch38/Homo_sapiens.GRCh38.dna.primary_assembly.fa`
- `/mnt/storage-box-1-gentle/Glen/IEGT/nanopore/reference/ensembl-116-grch38/Homo_sapiens.GRCh38.116.gtf`

The full source assets remain outside Git because they total tens of gigabytes.
Their complete SHA-256 digests, frozen reference/annotation hashes and cohort
receipt bindings are recorded in `data/source_bindings.tsv`.

## Deterministic extraction

`generate_fixture.sh` performs the following for each sample and gene:

1. before creating output, authenticate `source_bindings.tsv` against the
   committed fixture manifest and compare the actual size/full streamed SHA-256
   of every bound CRAM, cohort receipt, reference, annotation and sample sheet;
2. authenticate the concatenated per-gene GTF bytes against the committed
   combined-GTF checksum, then require indexed CRAMs and run `samtools quickcheck`;
3. select primary mapped records overlapping the full Ensembl 116 gene span
   with `samtools view -F 2308` and no MAPQ floor;
4. retain the selected alignments as indexed BAM;
5. export the same records as gzip-normalized FASTQ (`gzip -n`) with original
   base calls and qualities;
6. emit a readwise manifest and validate sequence/quality length parity;
7. emit the computed, verified source bindings and checksum every generated file.

The archival bindings and all 18 retained data hashes are unchanged by the
provenance/checkout repair. Storage-prefix relocation is explicit; the selected
source bytes are not replaceable with new data under historical hashes. Keep
the sources immutable throughout extraction. Preflight does not lock files or
authenticate an auxiliary index independently of its archived source.

Coordinates are GRCh38 / Ensembl 116, one-based closed in `loci.tsv`:

- TP73: `1:3652307-3736201`, plus strand.
- TP53: `17:7661779-7687546`, minus strand.

## Usage in GENtle

Primary intended roles:

- development and integration tests for cDNA SNP evidence and local phasing;
- read-level association of phased allele evidence with observed junction
  chains;
- comparison of SK-MEL-28 and WM793 under one frozen mapping/reference
  contract;
- fail-closed tests for low depth, uncallable loci, disconnected phase blocks,
  sequencing errors, and extra-cluster review flags.

## Non-claims

- These are RNA/cDNA reads, not matched genomic DNA.
- Additional RNA haplotypes do not by themselves establish copy-number change;
  editing, mapping, RT, chimera, subclone, contamination, and sequencing-error
  explanations remain open.
- Apparent allele switches across splice segments do not by themselves prove
  trans-splicing or cross-allelic splicing.
- Locus overlap is not a complete-transcript call.
- `barcode16` has sparse TP73 coverage, so absence or imbalance claims require
  explicit depth and callability reporting.
