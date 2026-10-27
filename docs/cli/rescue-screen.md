# Target Rescue Screen

<!-- manual-counterparts:start -->
No direct GUI counterpart: This chapter documents an automation interface, example or command-only workflow; no dedicated GUI chapter is paired with it.
<!-- manual-counterparts:end -->

### RNA-seq target rescue screen

`gentle_cli rescue-screen` is a standalone pilot calibration command for
screening RNA-seq reads against a small requested gene set after external
quantification. Typical inputs are reads that Salmon left unassigned or reads
that Salmon assigned to transcripts of interest.

This adapter is intentionally separate from the shared `rna-reads` engine
routes. It identifies redundant canonical k-mer support for requested genes;
it does **not** claim to call fusions, resolve splice isoforms, or establish a
definitive transcript of origin.

Minimal shape:

```sh
cargo run --bin gentle_cli -- rescue-screen \
  --transcript-fasta transcripts.fa \
  --gene TP73 \
  --gene SRSF1 \
  --reads salmon_unassigned.fastq.gz \
  --output-prefix artifacts/rescue/target_rescue
```

Accepted inputs include:

- one or more Ensembl-style transcript FASTA files through repeated
  `--transcript-fasta PATH`
- requested genes through `--genes SYM[,SYM,...]` and/or repeated
  `--gene SYM`
- streamed FASTA or FASTQ reads, plain or gzip-compressed, through repeated
  `--reads PATH`
- paired-end FASTA or FASTQ files through repeated `--read-pair R1,R2`; mates
  are read in lockstep and screened as one fragment without creating k-mers
  across the mate boundary
- an optional read-id allowlist
- optional Salmon bridge files through `--salmon-unmapped-names` and
  `--salmon-mappings-sam`
- optional local structural-evidence FASTA catalogs through `--exon-fasta`,
  `--junction-fasta`, `--exon-intron-boundary-fasta`, `--intron-fasta`, and
  `--genomic-region-fasta`; records use the same gene-symbol header tags as the
  transcript FASTA
- configurable `--kmer-len`, `--min-kmer-hits`, and repeated transcript-header
  `--gene-symbol-tag KEY` values

The output prefix produces deterministic reports:

- `<prefix>.gene_manifest.tsv`
- `<prefix>.gene_hits.tsv`
- `<prefix>.read_hits.tsv`
- `<prefix>.structural_hits.tsv`
- `<prefix>.run_metadata.json`
- `<prefix>.summary.json`
- `<prefix>.comparison.json` when Salmon bridge inputs define multiple read
  universes

Use this command as a transparent target-enrichment screen whose hit tables
can be audited, not as a substitute for an aligner or transcript quantifier.
Structural matches are labelled as known exon/junction evidence or as
retained-intron, intronic, or genomic-region candidates. They are not novel
exon, retained-intron, fusion, or isoform calls.

Annotated-junction and exon-intron-boundary FASTA records must include a
`boundary_after:N` header token, where `N` is the number of sequence bases to
the left of the junction or boundary. Only k-mers containing bases from both
sides are indexed, so a read matching only an adjacent exon or intron flank is
not counted as boundary evidence. Exon, intron, and explicitly supplied
genomic-region catalogs continue to use their full sequences.

For retained-intron and intronic candidates, `rna_anchored_hits` means that the
same read or paired fragment also passed the annotated-junction threshold for
the same gene. This makes an RNA-derived explanation more plausible, but does
not prove retention; unanchored hits remain visible because they may reflect
pre-mRNA, while also remaining compatible with genomic-DNA contamination. A
genomic-region catalog is never treated as RNA-specific and can simply be
omitted from RNA-focused runs.

The JSON and TSV reports distinguish physical input reads from evaluated
evidence units. A single-end record contributes one read and one evidence unit;
a paired fragment contributes two input reads and one evidence unit. Legacy v1
fields such as `selected_reads` retain their evidence-unit meaning, while the
additive `selected_input_reads` and `*_evidence_units` fields make the counting
unit explicit.

<!-- manual-page-nav:start -->

---

[CLI Manual](../cli.md) | [Previous: Observed Transcript Structures](observed-transcripts.md) | [Next: RNA Secondary Structure](rna-structure.md)
<!-- manual-page-nav:end -->
