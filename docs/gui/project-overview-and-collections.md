### Project Overview And Collections

<!-- manual-counterparts:start -->
Related CLI reference: [Collections](../cli/collections.md).
These topic links do not imply identical feature coverage across interfaces.
<!-- manual-counterparts:end -->

The project main window (lineage page) supports two views:

- `Table`: tabular lineage view with per-sequence actions
- `Graph`: node/edge lineage visualization; this is the default for fresh
  projects and sessions without a saved lineage workspace preference
- clickable overview counts above the view toggle summarize sequences, pools,
  analyses, containers, and arrangements:
  - sequence/pool/analysis counts focus the lineage graph/table area
  - container and arrangement counts resize the existing project split to bring
    the relevant lower section into view
- analysis artifacts (dotplots, flexibility tracks, persisted primer/qPCR
  design reports, persisted transcript-native protein-derivation reports,
  persisted UniProt genome projections, and
  sequencing-confirmation reports) appear in lineage as dedicated analysis
  nodes linked to their source sequences
- SVG export operations (sequence, dotplot, feature-expert/splicing, gel) are
  also materialized as lineage analysis nodes, linked by operation edges for
  provenance/script parity
- `Containers`: container list with kind/member-count, contents mode, open
  actions, per-container gel export, and pool-artifact export for declared
  exhaustive containers
- container contents mode is editable directly from the lineage `Containers` panel:
  - `Declared only`: the listed members are intended to be the full known contents of a clean vial/tube
  - `Known subset`: the listed members are known/measured constituents of a more complex sample that may also contain unlisted material
- `Arrangements`: ordered layouts across containers, with arrangement-level
  gel preview/export for `Serial` layouts and linked physical rack drafts for
  serial or plate layouts. `Plate` rows keep gel actions disabled with the
  shared engine explanation that only serial arrangements can render gels.
- Rack placement is a linked physical layer on top of arrangements, not a replacement for them:
  - arrangements keep the semantic sample order
  - racks/plates keep the physical A1-slot placement

<!-- manual-page-nav:start -->

---

[GUI Manual](../gui.md) | [Previous: Main window layout](main-window-layout.md) | [Next: Physical Arrangements And Project Layout](physical-arrangements-and-project-layout.md)
<!-- manual-page-nav:end -->
