### Lineage Node Actions

<!-- manual-counterparts:start -->
No direct CLI counterpart: This chapter describes desktop presentation or interaction; no dedicated command reference is paired with it.
<!-- manual-counterparts:end -->

Node click behavior in lineage `Graph` view:

- Single-click: selects a node (highlight only).
- Double-click on a single-sequence node: opens that sequence window.
- Double-click on a pool node: opens a pool-context window (Engine Ops visible,
  pool member distribution available).
- Double-click on an analysis node:
  - `Dotplot` / `FlexibilityTrack` / `PrimerDesign` / `QpcrDesign` /
    `SequencingConfirmation`: opens the source sequence window
  - `RnaReadInterpretation`: opens the RNA-read Mapping workspace on the saved
    report
  - `ProteinDerivation`: opens the Derived Protein Expert on the saved report
  - `ReverseTranslation`: opens the created coding-sequence window referenced
    by the saved report
    - the selected-node detail panel also shows a compact codon-choice
      diagnostics summary (`preferred`, fallback count, GC, realized suffix Tm)
  - `UniprotProjection`: opens the UniProt Protein Expert on the saved
    projection
- Single-click on the dedicated `Gibson cloning` operation node in graph view,
  or on the `Op` cell for Gibson-created outputs in the table: reopens the
  Gibson specialist with the saved plan loaded again for review.
- Single-click on a PCR-related operation glyph in graph view, or on the `Op`
  cell for PCR-/primer-design-created outputs in the table: opens the
  dedicated `PCR Designer` on the stored template sequence context.
- Single-click on a `RestrictionCloningPcrHandoff` analysis node: reopens the
  same `PCR Designer` on the originating template and restores the saved
  cloning-handoff selection (vector, enzymes, leaders, and saved handoff id).
- Right-click context menu (graph and table node-id cells):
  - starts with the clicked item's display label, kind and primary ID; keeps that
    subject while the pointer moves into the popup, without retargeting to another
    hovered/selected item. If the project or row disappears, reopen the menu
  - `Inspect node details`: selects that item to show its existing details and
    provenance in the project view
  - `Open sequence` / `Open pool`: opens the clicked item in its existing viewer,
    without combining molecules or running analysis
  - `Open saved analysis`: opens supported saved reports in their specialist
    viewer without rerunning them; other node types retain detail inspection
  - `PCR Designer...`: opens primer-pair setup for the clicked single DNA
    sequence, even if another sequence or PCR workspace was previously selected.
    RNA/protein and missing DNA remain disabled with an explanation. A pool is
    not silently reduced to one template. Running design remains a separate action
  - `Copy node ID`: copies the stable lineage node identifier.
  - `Copy sequence/report ID`: copies the primary row identifier, such as a
    sequence ID, qPCR/primer report ID, macro instance ID, or arrangement ID.
  - `Copy display label`: copies the row label exactly as shown, falling back
    to the primary ID when no display label exists.
  - `Copy row summary`: copies a compact multi-line text summary with the node
    ID, primary ID, operation provenance, and row-type-specific context.
  - `Rename (leaf only)`: updates the node display name (sequence name).
  - `Remove (leaf only)`: opens a confirmation dialog, then removes that
    sequence/node from project state on confirm.
  - both are disabled for non-leaf nodes in this first pass.
- Graph blank space explains how to inspect a node or add data via
  `File > Open sequence`.
- Right-click context menu on the lineage graph canvas also includes
  `Save Graph as SVG...`.
  - the exported SVG follows the same visible grouped/hub-projected graph model
    used by the main window, so Gibson apply operations save as the same single
    `Gibson cloning` operation node users see on screen
  - exported Gibson hub connector edges stay unlabeled because the operation is
    named by the hub box itself, which keeps hero/tutorial figures readable
  - the suggested filename defaults to the project name stem plus
    `.lineage.svg`
- Retrieval-pattern badges are shown for retrieval-derived sequence nodes:
  - `GENE` for `ExtractGenomeGene`
  - `REGION` for `ExtractGenomeRegion`
  - `GB` for `FetchGenBankAccession`
  - `RS` for `FetchDbSnpRegion`
  - clicking the badge reopens the matching retrieval dialog with saved inputs
    prefilled
  - if a genome anchor exists, it is shown in the same tooltip; GenBank source
    provenance and genome-anchor provenance are complementary
- Macro-instance nodes are rendered as dedicated box nodes.
  - Double-click behavior remains no-op (informational node type).
  - Click or use table `Inspect` to open a persistent macro detail panel below
    the lineage graph/table:
    - status/template/routine summary
    - typed input/output bindings
    - emitted operation IDs with operation-label summaries
    - quick-open buttons for sequence outputs
- Analysis nodes are rendered as dedicated rounded boxes.
  - Table/tooltip details expose:
    - artifact id
    - source sequence (+ reference sequence for pairwise dotplots)
    - mode/model/backend depending on artifact family
    - point/bin count for dotplot/flexibility artifacts
    - pair/assay count for persisted primer/qPCR reports
  - reasoning-guided dotplots are linked from the construct-reasoning graph
    that recommended them; citation-free dotplots are labelled as manual/no
    reasoning citation rather than failed
  - selected dotplot details show the citation as visually distinct `PASS`,
    `FAIL`, or `UNKNOWN`, plus rationale, action id, driving evidence ids,
    verification reasons, and a shortcut back to the recommending graph

<!-- manual-page-nav:start -->

---

[GUI Manual](../gui.md) | [Previous: Status History And Background Jobs](status-history-and-background-jobs.md) | [Next: Lineage Node Groups](lineage-node-groups.md)
<!-- manual-page-nav:end -->
