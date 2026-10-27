### Accession-Pinned TSS Profiles

<!-- manual-counterparts:start -->
Related CLI reference: [Accession-Pinned TSS Profiles](../cli/tss-profiles.md).
These topic links do not imply identical feature coverage across interfaces.
<!-- manual-counterparts:end -->

The existing Shell accepts `features tss-tfbs-profiles` and
`features tss-tfbs-profiles-export`, using the shared engine and exporting
readable SVG/PNG/PDF documents. See the [offline example](../tss_tfbs_profiles.md).
This first development slice has no dedicated native setup wizard. It does
not replace the existing locus-evidence viewer or change its historical figures.
The same Shell commands accept `--context-manifest FILE` to align CUT&RUN,
exon/CDS structure, annotated translation boundaries and TATA evidence above the
detail-page scores. No new wizard is required; see the
[offline context example](../tss_tfbs_profiles.md#cutrun-gene-structure-and-tata-context).
SVG hover details remain available alongside printable coordinate summaries.
The same Shell routes accept repeated `--genomic-motif-evidence FILE` for
[imported-hit triangles](../tss_tfbs_profiles.md#imported-duckdb-motif-hits) in
detailed exported plots. This does not replace the native viewer's strip overlay.
`--formats svg,genbank,embl` also exports [annotated TSS windows](../tss_tfbs_profiles.md#annotated-tss-sequences)
for the normal sequence-file importer; this does not claim a new native TSS viewer.
Use `--scale-mode shared_across_tss` for the same accession's common Y range
across all supplied TSSs and both strands. Each page states its scale policy;
independent scaling remains the default. Older background-tail score reports
need [rescoring, not just re-export](../tss_tfbs_profiles.md#rescoring-existing-reports).

<!-- manual-page-nav:start -->

---

[GUI Manual](../gui.md) | [Previous: Promoter Cofactors](promoter-cofactors.md) | [Next: TATA-Box Evidence In The DNA Viewer](tata-box-evidence-in-the-dna-viewer.md)
<!-- manual-page-nav:end -->
