# Integrated locus context and selected-TSS profiles

`scripts/compose_locus_tss_profile_pdf.py` combines two already verified report
families without rescoring either one:

1. one complete locus-context SVG (transcripts, reporter architectures,
   regulatory regions, occupancy/chromatin lanes and promoter similarity);
2. the detailed TFBS pages for only the TSSs marked selected in an
   accession-pinned `gentle.tss_tfbs_profiles.v1` report.

The output is one PDF per gene. The context SVG is page 1. Selected TSS pages
follow in the deterministic order of the scored report. Other annotated TSSs
remain in the complete TSS supplement. Every page is 1400 pixels wide and every
detailed score axis uses the locus renderer's exact x=255..1050 plot frame, so
context and detail pages remain horizontally registered in a PDF viewer.

The composer also writes one selected-TSS FASTA per gene. It contains exactly
the transcript-oriented -500/+200 records shown on the detailed pages, in the
same order and with their original bound headers. It is a filtered projection of
the verified input bundle, not a new genome extraction.

The composer fails before publication unless:

- the locus SVG matches its locus-report receipt;
- the TSS report, index and every selected SVG match the TSS receipt;
- the index and report hashes agree;
- selected promoter IDs belong to the requested gene and occur exactly once;
- every selected TSS has gene, chromosome, strand, coordinate, transcript and
  sequence-digest metadata and bound selection evidence;
- every selected TSS coordinate falls inside a background band in the bound
  locus SVG;
- context and selected-detail SVGs share the 1400-pixel page and x=255..1050
  plot frame;
- the input manifest, `SHA256SUMS`, source FASTA and record-level sequence
  digests all agree with the TSS-profile receipt.

The evidence-window and displayed-window sequence hashes are deliberately not
joined: the earlier reporter selection used a -2000/+200 window, while the TFBS
detail uses -500/+200. Their stable join is the selected `promoter_id`, with
gene and physical TSS geometry checked independently.

Example for one gene:

```sh
python3 scripts/compose_locus_tss_profile_pdf.py \
  --gene TGFB1 \
  --locus-svg /path/to/TGFB1_with_TSS_similarity.svg \
  --locus-receipt /path/to/TGFB1_with_TSS_similarity.receipt.json \
  --tss-report /path/to/tss-profiles/report.json \
  --tss-index /path/to/tss-profiles/index.json \
  --tss-receipt /path/to/tss-profiles/receipt.json \
  --tss-bundle-manifest /path/to/target-tss-bundle/manifest.json \
  --gentle-cli /path/to/committed/gentle_cli \
  --producer-revision 0.1.0-internal.10+git.REVISION \
  --output-pdf /fresh/output/TGFB1_integrated_locus_selected_TSS_TFBS.pdf \
  --output-fasta /fresh/output/TGFB1_selected_TSS_minus500_plus200.fasta \
  --output-receipt /fresh/output/TGFB1_integrated_locus_selected_TSS_TFBS.receipt.json
```

The output parent may be created, but all three named outputs must be absent. The
composer stages the PDF, FASTA and receipt before promotion. The receipt binds
both input families, the ordered pages, exact selected-TSS joins, producer
script and executable, PDF/FASTA hashes, horizontal frame, score policy,
reference statement and scientific non-claims.

The composite is a presentation artifact. JASPAR tracks remain sequence-model
predictions, not measured binding, affinity, occupancy, promoter activity or
functional regulation. Reporter-selection evidence likewise does not prove TSS
usage, direct binding or construct activity.
