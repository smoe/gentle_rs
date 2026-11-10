# Real ONT cDNA fixture for TP53/TP73 allele and splice work

This owner-authorized fixture contains a bounded subset of real Oxford
Nanopore cDNA reads from two melanoma cell lines:

| Barcode | Cell line | Library |
| --- | --- | --- |
| `barcode15` | SK-MEL-28 | native regular cDNA, no TP53/TP73 target primer |
| `barcode16` | WM793 | native regular cDNA, no TP53/TP73 target primer |

The fixture is intended for developing and reviewing read-level SNP evidence,
local haplotypes, and the association between allele evidence and physical
splice chains. It is not a clinical dataset and does not establish genomic
copy number, trans-splicing, RNA editing, or complete full-length isoforms.

## Contents

- `data/annotation/ensembl116.tp53_tp73.gtf`: TP53 and TP73 Ensembl 116
  annotation used by the original cohort.
- `data/reference/ensembl116.tp53_tp73.genomic_regions.fa`: exact GRCh38
  genomic locus sequences for the two annotated gene spans. FASTA headers
  retain the source regions; consumers must not mistake them for complete
  chromosome contigs.
- `data/reads/*.primary.bam`: primary mapped alignments overlapping the named
  gene span, with no MAPQ floor. The cohort-compatible exclusion mask is
  `-F 2308` (unmapped, secondary, and supplementary alignments excluded).
- `data/reads/*.primary.fastq.gz`: the same selected records as FASTQ,
  preserving base calls and qualities carried by the CRAM/BAM record.
- `data/read_manifest.tsv`: read ID, sample, gene, flag, coordinate, MAPQ,
  CIGAR, query length, and mean base quality.
- `data/samples.tsv`, `data/loci.tsv`, and `data/source_bindings.tsv`: sample,
  coordinate, and source provenance.
- `data/SHA256SUMS`: checksums for every generated fixture file except the
  checksum file itself.

The existing transcript FASTAs in `test_files/fixtures/mapping/` provide the
compact TP53/TP73 transcript comparison set. They are not duplicated here.

## Regeneration on Glen's host

The immutable full inputs live on the mounted Hetzner Storage Box. Generate
into a fresh directory; the script refuses to overwrite a non-empty `data/`
tree:

```bash
bash test_files/fixtures/cdna_haplotype_tp53_tp73/generate_fixture.sh \
  /tmp/cdna_haplotype_tp53_tp73_regenerated
```

Compare the regenerated `data/SHA256SUMS` and files with the committed
fixture. See `PROVENANCE.md` for exact source bindings and interpretation
limits.

## Intended use

This fixture is an integration/development corpus, not a tiny unit-test oracle.
Tests should use deterministic projections from it and keep these distinctions:

1. reference support, alternate support, missing coverage, and uncallable
   evidence are different states;
2. phase labels are local to connected phase blocks unless reads link them;
3. a haplotype/splice association requires evidence on the same read;
4. apparent within-read allele switches must remain candidates for sequencing,
   alignment, RT, or chimera artifacts;
5. more than two RNA haplotype clusters are a review trigger, not a copy-number
   call without matched DNA evidence.
