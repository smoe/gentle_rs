# Integrated locus context and selected-TSS profiles

`scripts/compose_locus_tss_profile_pdf.py` combines two already verified report
families without rescoring either one:

1. one complete locus-context SVG (transcripts, reporter architectures,
   regulatory regions, occupancy/chromatin lanes and promoter similarity);
2. the detailed TFBS pages for only the TSSs marked selected in an
   accession-pinned `gentle.tss_tfbs_profiles.v1` report.

The output is one PDF per gene. The context SVG is page 1. Selected TSS pages
follow in the deterministic order of the scored report. Other annotated TSSs
remain in the complete TSS supplement.

The composer fails before publication unless:

- the locus SVG matches its locus-report receipt;
- the TSS report, index and every selected SVG match the TSS receipt;
- the index and report hashes agree;
- selected promoter IDs belong to the requested gene and occur exactly once;
- every selected TSS has gene, chromosome, strand, coordinate, transcript and
  sequence-digest metadata and bound selection evidence;
- every selected TSS coordinate falls inside a background band in the bound
  locus SVG.

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
  --gentle-cli /path/to/committed/gentle_cli \
  --producer-revision 0.1.0-internal.10+git.REVISION \
  --output-pdf /fresh/output/TGFB1_integrated_locus_selected_TSS_TFBS.pdf \
  --output-receipt /fresh/output/TGFB1_integrated_locus_selected_TSS_TFBS.receipt.json
```

The output parent may be created, but both named outputs must be absent. The
composer stages the PDF and receipt before atomic promotion. The receipt binds
both input families, the ordered pages, exact selected-TSS joins, producer
script and executable, output hash, score policy, reference statement and
scientific non-claims.

The composite is a presentation artifact. JASPAR tracks remain sequence-model
predictions, not measured binding, affinity, occupancy, promoter activity or
functional regulation. Reporter-selection evidence likewise does not prove TSS
usage, direct binding or construct activity.
