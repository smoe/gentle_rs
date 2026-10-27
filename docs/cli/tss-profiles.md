# Accession-Pinned TSS Profiles

<!-- manual-counterparts:start -->
Related GUI reference: [Accession-Pinned TSS Profiles](../gui/accession-pinned-tss-profiles.md) | [TATA-Box Evidence In The DNA Viewer](../gui/tata-box-evidence-in-the-dna-viewer.md).
These topic links do not imply identical feature coverage across interfaces.
<!-- manual-counterparts:end -->

## Accession-Pinned TSS Profiles

`features tss-tfbs-profiles` validates a pinned FASTA bundle and exact JASPAR
panel, computes both orientations and publishes audited profile documents.
`features tss-tfbs-profiles-export` consumes an existing report without rescoring.
See [complete offline commands and input requirements](../tss_tfbs_profiles.md).
These shared-parser routes also work in the GUI Shell and through typed MCP
`op`; output directories are fresh and writes require normal confirmation.
Both accept `--context-manifest FILE` for sequence-verified CUT&RUN/chromatin,
exon/CDS/translation-marker and TATA rows on detailed pages, without rescoring.
See the [context example and input contract](../tss_tfbs_profiles.md#cutrun-gene-structure-and-tata-context).
Use `--formats svg,vector_pdf` for an interactive SVG peer plus a static vector
PDF with selectable embedded text. The existing `pdf` token remains the
lossless raster-backed compatibility format; PDF viewers do not portably retain
SVG title-hover behavior.
Repeat `--genomic-motif-evidence FILE` to attach saved DuckDB-query reports as
separate [strand/span/score triangle lanes](../tss_tfbs_profiles.md#imported-duckdb-motif-hits)
and a complete imported-hit TSV, without rescoring or querying the package.
Both commands accept `--scale-mode independent|shared|shared_across_tss`.
`--formats svg,genbank,embl` additionally exports annotated, sequence-verified TSS
windows; [requirements and legacy replay](../tss_tfbs_profiles.md#annotated-tss-sequences).
The last mode shares one accession's range across all supplied TSSs and both
strands, without cross-matrix calibration. See [scores, scaling and required
rescoring of older tail reports](../tss_tfbs_profiles.md#scores-and-geometry).

The CLI prints an export summary, not the full scored report, for these routes
(also through `shell`, direct JSON `op`, and JSON `workflow`). It retains artifact
paths, hashes, counts, warnings and verification limits. Put `--full-report`
**before the command** to restore full result JSON on stdout. Exported
`report.json` is always complete; compute-only operations without an export,
GUI Shell and MCP retain their full result contract. See the
[stdout compatibility details](../tss_tfbs_profiles.md#cli-output).

<!-- manual-page-nav:start -->

---

[CLI Manual](../cli.md) | [Previous: Promoter Cofactors](promoter-cofactors.md) | [Next: TFBS and JASPAR Scans](tfbs.md)
<!-- manual-page-nav:end -->
