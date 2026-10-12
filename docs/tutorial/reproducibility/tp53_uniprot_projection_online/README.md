# TP53 UniProt projection acceptance fixture

This directory retains the small locus and structured read-backs used to
review tutorial 06.04 on 2026-10-01. It does not replace the chapter's full
online `PrepareGenome` / `ExtractGenomeGene` workflow.

## Provenance

- `tp53_grch38_ensembl116.gb` combines the Ensembl 116 TP53 transcript/CDS
  geometry already retained in `docs/figures/tp53_ensembl116_panel_source.gb`
  with GRCh38 reference bases from UCSC `hg38`,
  `chr17:7661778-7687546` (0-based, end-exclusive).
- The UCSC request was
  `https://api.genome.ucsc.edu/getData/sequence?genome=hg38;chrom=chr17;start=7661778;end=7687546`
  on 2026-10-01. Its 25,768-base DNA payload has SHA-256
  `31d07cd16e14c12ae6a8049f8e9614be31713428d717594e8927cdf34f97dd6d`.
- The complete GenBank fixture has SHA-256
  `a6694752cd58bd16e6e5565725cc6d2fa3b127606f8ce3c6ee7890e327fecec3`.
- UniProt `P04637` was fetched by GENtle from
  `https://rest.uniprot.org/uniprotkb/P04637.txt` on 2026-10-01.
- Ensembl's gene and region sequence REST endpoints both returned HTTP 500
  during this review. The UCSC GRCh38 sequence therefore enables bounded
  local acceptance without claiming that the full Ensembl download path ran.

## Reproduce the retained read-backs

Use a fresh state path from the repository root:

```bash
STATE=/tmp/tp53-uniprot-review.gentle.json
gentle_cli --state "$STATE" shell \
  '/open file docs/tutorial/reproducibility/tp53_uniprot_projection_online/tp53_grch38_ensembl116.gb --id grch38_tp53'
gentle_cli --state "$STATE" shell \
  'uniprot fetch P04637 --entry-id P04637'
gentle_cli --state "$STATE" shell \
  'uniprot map P04637 grch38_tp53 --projection-id tp53_uniprot_p04637'
gentle_cli --state "$STATE" shell \
  'uniprot feature-coding-dna tp53_uniprot_p04637 DNA-binding --mode both --speed-profile human'
```

The review mapped six UniProt-linked transcripts and retained 16 explicit
unmapped-transcript warnings. The feature query produced one match on
`ENST00000269305.9`: amino acids 368-387, exon 11, 60 genomic coding bases and
a separately labelled 60-base translation-speed-oriented alternative.

Exact captured outputs:

- `map-result.json` SHA-256
  `bbdb7dad51f305061d24770f9e5bc5554ef34a1e330b8ee9556d7086a6028be5`
- `feature-coding-dna.json` SHA-256
  `d12fa76661ec541f77fc4c3bfda15078edbc159a39d6d278e756ad67adc45f24`

The fast checkout regression
`scripts.test_tutorial_checkouts.TutorialCheckoutTests.test_uniprot_review_evidence_hashes_survive_both_checkout_modes`
copies these three retained files into disposable Git repositories and verifies
the recorded hashes under LF and CRLF checkout settings. The scoped LF rules
are contractual; the test also proves that removing them changes the hashes.
