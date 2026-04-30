# Example Asset Provenance

## `cdna_assay_demo.gb`

- Origin: synthetic, hand-crafted GENtle demo locus.
- Recreation: create a 32 bp linear GenBank record with one two-exon `TEST1`
  transcript (`join(1..12,21..32)`) and the sequence
  `ATGAAACCCGGGTTTTTTTTCCCAAATTTGGG`.
- Used by: cDNA PCR/qPCR assay-test workflow and ClawBio direct request
  examples.

## `cdna_assay_nonspecific_demo.gb`

- Origin: synthetic, hand-crafted GENtle demo locus.
- Recreation: create a 72 bp linear GenBank record with two `NONSPEC`
  transcripts sharing exon 1 (`1..12`) and using alternative second exons
  (`21..32` and `41..72`) so the same `AAACCC` / `CCCAAA` primer pair detects
  two transcript-derived cDNA products of different sizes.
- Used by: cDNA PCR/qPCR product-materialization and product-gel workflow and
  ClawBio examples for nonspecific assay product visualization.
