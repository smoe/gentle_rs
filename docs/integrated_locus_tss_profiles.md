# Integrated locus context and selected-TSS profiles

`scripts/compose_locus_tss_profile_pdf.py` combines two already verified report
families without rescoring either one:

1. one complete locus-context SVG (transcripts, reporter architectures,
   regulatory regions, occupancy/chromatin lanes and promoter similarity);
2. the detailed TFBS pages for only the TSSs marked selected in an
   accession-pinned `gentle.tss_tfbs_profiles.v1` report.

The output is one PDF per gene. The context SVG is page 1. Selected TSS pages
follow in the deterministic order of the scored report, retaining all continuation
pages for a tall TSS. Each selected TSS still contributes only one FASTA record.
Other annotated TSSs
remain in the complete TSS supplement. Every page is 1400 pixels wide and every
detailed score axis uses the locus renderer's exact x=255..1050 plot frame, so
context and detail pages remain horizontally registered in a PDF viewer.

The composer also writes one selected-TSS FASTA per gene. It contains exactly
the transcript-oriented -500/+200 records shown on the detailed pages, in the
same order and with their original bound headers. It is a filtered projection of
the verified input bundle, not a new genome extraction.

The composer fails before publication unless:

- the locus SVG and source locus JSON match their locus-report receipt;
- the source locus JSON and SVG share a panel identity, and the locus and TSS
  report agree on reference assembly, chromosome and strand;
- every exact matrix accession in the detail panel occurs once in the locus
  overview's resolved JASPAR tracks (factor names cannot substitute for accessions);
- the TSS report, index and every selected SVG match the TSS receipt;
- the index and report hashes agree;
- selected promoter IDs belong to the requested gene and occur once in the
  selected-window list; continuation pages may repeat IDs but not SVG files;
- every selected TSS has gene, chromosome, strand, coordinate, transcript and
  sequence-digest metadata and bound selection evidence;
- every selected TSS coordinate falls inside a background band in the bound
  locus SVG;
- context and selected-detail SVGs share the 1400-pixel page and x=255..1050
  plot frame;
- the input manifest, `SHA256SUMS`, source FASTA and record-level sequence
  digests all agree with the TSS-profile receipt;
- the digest calculated from actual FASTA bases matches the header, manifest
  and report, without allowing header fields to overwrite the calculated value;
- the FASTA header, manifest and report have the same transcript membership
  and strand-aware -500/+200 genomic interval, not merely the same 701-bp length;
- the PDF renderer reports the actual font identities/hashes used for every
  page, from the same rasterization whose pixels enter the PDF.

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
  --locus-report /path/to/TGFB1_locus.report.json \
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

The output parents may be created, but all three named outputs and their `.partial`
paths must be absent. The composer stages the PDF, FASTA and receipt before
promotion. Same-directory hard links publish without replacing existing files;
the output filesystem must support them. The receipt is published last as the
completion marker. Caught errors and interruptions roll back this invocation's
outputs and partial files, preserving pre-existing/racing files. This is not a
single filesystem transaction across three paths: after process termination or
power loss, files without a receipt are incomplete and must not be distributed.
The receipt binds
both input families, the ordered pages, exact selected-TSS joins, producer
script and executable, PDF/FASTA hashes, horizontal frame, score policy,
reference statement and scientific non-claims.

## Refreshing The September 10 Bundle

The bundle produced at `8fd677c2` has internally matching PDF/FASTA hashes but
its CD44/TGFB1 overviews repeat `MA0815.1` under three TFAP2C labels. Do not fix
those PDFs by relabeling curves. Shared TF query expansion now retains resolved
matrix accessions rather than converting them back to ambiguous factor names.

Glen should regenerate the locus score reports/SVGs from the original requests
with explicit `MA0524.3`, `MA0814.3` and `MA0815.1` source IDs, then rerun the
bound similarity composition and this compositor with the matching locus JSON.
Keep the existing CUT&RUN/H3K4me3 inputs, track completeness, TSS bands and both
strand orientations. Retain the correct detail scores unless their input/policy
changes; their existing SVG pages can be re-used for composition and rasterized
by the updated audited `gentle_cli`.

Recheck all five genes, 13 selected FASTAs and their coordinates/digests, every
overview/detail matrix accession, continuation order, used-font hashes and
final output receipts. Rebind the web manifest to the new producer revision;
do not mix old PDFs with new receipts or present a hash check as independent
genome re-extraction. The original locus JSON/SVG and TSS export inputs remain
on Glen's host, so repository tests alone cannot certify the replacement bundle.

The composite is a presentation artifact. JASPAR tracks remain sequence-model
predictions, not measured binding, affinity, occupancy, promoter activity or
functional regulation. Reporter-selection evidence likewise does not prove TSS
usage, direct binding or construct activity.
