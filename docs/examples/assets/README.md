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

## `tp73_promoter_artifact_demo.gb`

- Origin: synthetic, hand-crafted GENtle promoter-signoff demo locus.
- Recreation: create a 249 bp linear GenBank record with one `TP73`-labeled
  gene, two transcript features sharing the same 5' boundary
  (`ENSTTP73DEMO1`, `ENSTTP73DEMO2`), one alternative-start transcript
  (`ENSTTP73DEMO3`), one exact shared promoter annotation (`61..116`), one
  synthetic `SP1` TFBS (`81..89`), one synthetic promoter variant (`95`), one
  synthetic LINE/L1 repeat interval (`105..116`), one TP73-like TFBS
  (`150..158`), and one CUT&RUN-labeled track interval (`130..152`).
- Exact sequence:
  `AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAACCCCCCCCCCCCCCCCCCCCGGGGCGGGGTTTTTTTTTTAAAAAAAAAACCCCCCCCCCCCCCCCCCCCGGGGGGGGGGGGGGGGGGGGATGTGTAACTTTTTTTTTTTTTTTTTTTTGGGGGGGGGGGGGGGGGGGGCCCCCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA`
- Used by: offline promoter-design artifact slice workflow for
  alternative-promoter grouping, promoter evidence matrix review, TFBS
  score-track SVG export, and TFBS similarity JSON export.
