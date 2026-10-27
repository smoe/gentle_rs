### Splicing Structure And Evidence Views

<!-- manual-counterparts:start -->
Related CLI reference: [Cryptic Splicing](../cli/splicing.md).
These topic links do not imply identical feature coverage across interfaces.
<!-- manual-counterparts:end -->

- The Splicing Expert has `Structure`, `Evidence`, `Locus figure`, and
  `Cryptic screen` tabs.
  `Evidence` is a thin GUI over `FeatureExpertTarget::IsoformEvidence`:
  - isoform-panel choices come from panels imported for the active sequence.
    One compatible panel is seeded directly; several panels require an
    explicit choice; no compatible panel leaves composition blocked with an
    import hint. An explicit choice is stored with the project per sequence and
    restored after reopening it. `Use panel editor ID` remains a visibly
    separate advanced action and does not make an incompatible global editor
    default valid for the current locus
  - choose an imported isoform panel and optionally supply annotation-release
    text, persisted RNA-read/qPCR report ids, probe-evidence JSON, cDNA/EST
    resource JSON, an expression TSV, and projected occupancy track names
  - `Inspect evidence` caches the shared `gentle.gene_isoform_evidence.v2`
    report for the current splice group; opening a different group clears that
    cached result
  - transcript rows show biological 5'->3' exon order beside genomic-ascending
    order, while junction rows show stable assembly-local ids and
    transcript-oriented donor/acceptor coordinates
  - specificity, abundance, responsiveness, and assayability remain separate;
    the evidence ledger and provenance can be expanded for audit. Repeated
    contrast/source measurements remain individually inspectable; incompatible
    units are not collapsed. Recommendation tiers are rule-based guidance, not
    weighted evidence scores
  - selected BED/BigWig occupancy tracks are summarized as source-specific
    lanes and rendered beneath transcript models in the shared SVG with one
    common score scale; `*` selects every projected track in the gene span.
    Occupancy remains locus-level evidence and is not assigned to an isoform
  - `Copy report JSON` copies the machine-readable report and `Render SVG`
    delegates to the shared `RenderFeatureExpertSvg` operation
  - existing qPCR rows are candidates only. New assay design remains an
    explicit action in the Structure tab, and the GUI does not promote array
    overlap or missing evidence into biological validation
  - evidence status labels use one shared vocabulary: `observed evidence`,
    `candidate association`, `design constraint`, `not evaluated`, and
    `unresolved evidence`. These labels interpret the existing stable wire
    values; they do not promote candidates or overlaps into validation
- The standard linear DNA map and `Structure` tab also offer `Source annotation
  comparison (read-only)`. Load an enriched locus JSON there or reuse the report
  from `Locus figure`. These bounded, scrollable lanes use one engine projection,
  the same source filter, and the DNA-map viewport. Click highlights a comparison
  row; double-click fits its span without assigning PCR selection. Selecting a
  loaded transcript highlights full-chain/strand matches, not merely shared
  exons or matching labels. Thin/thick boxes show exon/CDS geometry; ticks show
  annotated starts. Hover exposes exact source provenance. Stale bindings and
  circular sequences do not receive live overlays. Comparison records never
  enter the actionable project transcript list or change primer-design scope.

<!-- manual-page-nav:start -->

---

[GUI Manual](../gui.md) | [Previous: Splicing Support And Transcript Actions](splicing-support-and-transcript-actions.md) | [Next: Locus Figure And Evidence Tracks](locus-figure-and-evidence-tracks.md)
<!-- manual-page-nav:end -->
