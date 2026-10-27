### Conservation Workspace

<!-- manual-counterparts:start -->
No direct CLI counterpart: This chapter describes desktop presentation or interaction; no dedicated command reference is paired with it.
<!-- manual-counterparts:end -->

- `Conservation...` on a saved region opens a separate, resizable workspace:
  - the `Genomic Region Conservation` command-palette entry opens the saved-region
    manager for the active or first project sequence, where you choose the region
    and select `Conservation...`. It shares the Saved Genomic Regions opener:
    existing viewers are focused and pending viewers retain the manager request.
    An empty project opens sequence retrieval first
  - `Search request` edits the shared request: query resource, catalog/cache,
    explicit target genome IDs, required/optional status, expected ortholog
    loci and evidence sources, and every search-policy threshold. Empty targets
    still mean all validated local indexes. No species or orthology is guessed
  - `Import request...` and `Export request...` exchange the same JSON used by
    CLI/workflows. Imports must refer to the selected saved region; stale digest
    bindings are rejected and omitted bindings are pinned to its current digest.
    Editing/importing does not run BLAST or change an existing report. Request
    controls are disabled during a run; the worker uses the exact submitted copy
  - `Run local screen` snapshots the engine and performs preflight, local BLAST,
    projection, exact-support block calling, and composition away from the UI
    thread. The status line reports phase, target ordinal, elapsed time, and a
    final target/locus/block count
  - `Cancel` cooperatively stops the active BLAST child via the shared process
    supervisor. Closing the workspace or switching to another region also
    requests cancellation; cancelled or mismatched late results cannot replace
    the displayed report. Elapsed time continues updating during BLAST, which
    reports a heartbeat rather than inventing an intra-target percentage/ETA
  - `Open report...` validates and displays an existing content-bound homology
    report only when it still matches the exact saved query region. This allows
    a CLI/workflow result with explicit orthology declarations to be inspected
    without rerunning BLAST
  - no genome is downloaded or indexed from this window. Missing optional
    indexes are shown as unavailable and unrelated DNA-viewer work remains
    available
  - target readiness, conserved blocks, same-genome alternatives, and the
    query-referenced alignment are inspectable together. Target insertion bases
    stay in exported JSON and never create extra alignment columns
  - `Promoter similarity...` is a focused entry point into the same workspace.
    It requires the saved region's prepared-genome identity, selects that genome
    as the same-genome target, enables transcript-promoter annotation, and starts
    with 40 bp / 80% identity and a -2,000/+200 bp transcript window. Review or
    edit these declared parameters before `Run local screen`
  - the promoter recurrence matrix uses one row per distinct strand-aware
    genomic promoter window and retains every gene/transcript mapping for shared
    TSSs. Horizontal position is the selected query region, blue intensity is
    nucleotide identity, and the number inside each block is its 5'-to-3' order
    in the target promoter. A red outline begins a new structural segment when
    query order or alignment orientation changes; GENtle never joins across it
  - window, gene, and transcript frequencies are displayed separately. If the
    BLAST HSP or retained-locus budget is reached, the matrix labels frequencies
    as lower bounds. Similarity, rarity, and rearrangement remain reporter-design
    evidence, not promoter-activity or fragment-sufficiency claims
  - conserved-block and alignment lists instantiate only visible rows. The
    alignment has no former 50-locus display cap; selecting a block scrolls to
    its query tile without laying out the complete sequence on every repaint
  - other saved evidence regions fully contained in the query can be selected
    for the shared promoter-module assessment. Each PASS/NO rule is displayed;
    the result remains a reporter-design hypothesis rather than a functional
    conclusion
  - JSON and SVG export consume the shared report. Saving a selected block as a
    portable region is a separate explicit project mutation that retains the
    homology report digest and non-claims
  - semantic IDs cover the window, run/open/status, JSON/SVG export, module
    assessment, and block-save controls for Linux/Xvfb acceptance testing

<!-- manual-page-nav:start -->

---

[GUI Manual](../gui.md) | [Previous: Saved Genomic Regions](saved-genomic-regions.md) | [Next: Cryptic Splicing Screen](cryptic-splicing-screen.md)
<!-- manual-page-nav:end -->
