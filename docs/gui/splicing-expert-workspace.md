### Splicing Expert Workspace

<!-- manual-counterparts:start -->
Related CLI reference: [Cryptic Splicing](../cli/splicing.md).
These topic links do not imply identical feature coverage across interfaces.
<!-- manual-counterparts:end -->

- The splicing expert window is a free-standing top-level window and can be
  moved outside the DNA sequence window bounds.
- The splicing expert window opens slightly slimmer by default and its content
  is scrollable in both directions, so wide RNA-read tables can be inspected
  without stretching the window and very large transcript sets no longer hide
  the `Transcript x exon matrix` section.
- The first Splicing Expert frame now opens in a lightweight warm-up state
  before painting the full expert body on the next repaint, which keeps
  context-menu launches responsive on heavier loci.
- The same lightweight warm-up behavior now also applies to the dedicated
  `RNA-read Mapping` workspace so heavy loci do not try to paint the full
  mapping control surface in the same frame that opens the window.
- The same lightweight warm-up behavior now also applies to the dedicated
  `Promoter design` workspace so its full control surface is not rendered in
  the same frame that opens the window.
- Very large `Transcript x exon` or `Exon -> exon` matrices may open collapsed
  by default to reduce idle CPU usage; expand the section header to render the
  full table.
- Hover a transcript label in `Transcript x exon` to inspect its UniProt
  reference status and source xrefs. Evidence comes from entries already loaded
  into this project (`uniprot import-swissprot` or explicit `uniprot fetch` in
  the shared shell); refresh/reopen the expert view after loading evidence.
  An entry is relevant only through an exact normalized transcript xref or an
  exact stable Ensembl gene-ID xref. Names, labels, protein lengths and unrelated
  projections do not count. With no exact relation, status is `not_evaluated`.
  Gene-only relevance does not make a transcript `referenced`; `not_referenced`
  is scoped to the relevant loaded entries, not a global absence claim. The
  shared JSON and exported SVG hover titles retain the same source hashes and
  identifiers. Biological geometry and engine payload order are unchanged.
  The canvas and transcript/exon matrix share a display-only grouping: exact
  UniProt-referenced isoforms, other evaluated transcripts, then not evaluated.
  Target first within each group, then normalized transcript ID and feature ID.
  Headers are omitted when all lanes are unevaluated. One compact accession/
  review-status badge per matched lane carries `+N` for further entries and `?`
  for a loaded-evidence diagnostic, which can concern gene-level relevance
  rather than transcript ambiguity; the shared legend explains both symbols.
  Full evidence remains in hovers. Header clicks select nothing.
  Missing/duplicate lane-to-matrix joins are explicitly unevaluated; orphan
  matrix records remain listed, and missing exon cells show `n/a`, not absence.
  Every boundary tick remains visible, with rose exceptional markers. Full
  motif words appear once per oriented boundary/pair in the shared summary;
  expand **Per-transcript boundary details** for all original annotations.
  Summary labels, readable source-note hovers, per-transcript motif rows and
  class counts reuse the content-keyed presentation cache. SVG uses the same
  summary wording, not internal Rust debug fields. Non-canonical pairs retain
  exceptional ticks/dots even when their own dinucleotide is canonical.
  Unevaluated legacy views keep their 26/28-px lane spacing; evaluated groups
  reserve 36 px for headings/badges. Dense SVG matrices retain all cells within
  the canvas with nonoverlapping positive widths; a legend diagnostic discloses
  sub-2-pixel compression rather than implying missing exons or transitions.

<!-- manual-page-nav:start -->

---

[GUI Manual](../gui.md) | [Previous: Gene Set Inspector And Collection Operations](gene-set-inspector-and-collection-operations.md) | [Next: Splicing Support And Transcript Actions](splicing-support-and-transcript-actions.md)
<!-- manual-page-nav:end -->
