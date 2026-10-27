# Source-Coherent Transcripts

<!-- manual-counterparts:start -->
Related GUI reference: [Splicing Support And Transcript Actions](../gui/splicing-support-and-transcript-actions.md).
These topic links do not imply identical feature coverage across interfaces.
<!-- manual-counterparts:end -->

## Source-Coherent Transcripts

Existing `gene-locus prepare @request.json` and TSS `--context-manifest` inputs
accept optional `transcript_annotation_sources` for bound Ensembl/RefSeq
structure presentation. No network call or preferred-TSS selection is implied.
Use the enriched locus JSON for both overview and detailed pages. The existing
PDF compositor can also emit ordered SVG/HTML pages and a receipt-bound ZIP.
See [input examples and limitations](../transcript_source_presentation.md).

<!-- manual-page-nav:start -->

---

[CLI Manual](../cli.md) | [Previous: Gene Groups and Gene Sets](gene-sets.md) | [Next: Cryptic Splicing](splicing.md)
<!-- manual-page-nav:end -->
