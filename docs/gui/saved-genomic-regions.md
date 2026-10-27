### Saved Genomic Regions

<!-- manual-counterparts:start -->
Related CLI reference: [Isoform Architecture and Genomic Regions](../cli/isoform-architecture.md).
These topic links do not imply identical feature coverage across interfaces.
<!-- manual-counterparts:end -->

- `Saved genomic regions` keeps reusable interval configuration out of the DNA
  map itself while remaining one click away from it:
  - open it with `Regions...` in the DNA toolbar, `Genome > Saved Genomic
    Regions...`, the command palette entry of the same name, or by selecting
    bases and choosing `Save/share selected genomic region...` from the map
    context menu
  - the same destination is reachable headlessly as the shared UI intent
    `ui open saved-genomic-regions` (aliases `genomic-regions`, `regions`),
    with `ui close saved-genomic-regions` reporting whether a manager was
    actually open. While a sequence is still loading the focus is queued and
    applies once the load completes
  - the resizable manager lists set, label, explicit assembly, 1-based
    inclusive coordinates, strand, selection method, and evidence availability.
    Its copy actions provide human coordinates, BED6 (0-based half-open), or
    the complete canonical ROI JSON
  - JSON is the lossless interchange. `BED + manifest...` writes BED6 plus a
    content-bound sidecar; importing bare BED requires an explicit reference in
    the shared request path rather than guessing from a chromosome label
  - CUT&RUN support-window and source-bound Ensembl Regulation tables expose
    capture actions that call the same `CaptureGenomicRegion` operation as the
    manager. Their evidence and conservative non-claims remain source-specific
  - `Locus figure` accepts saved region-set IDs and draws their current local
    projections on the shared coordinate axis. Stale projections remain typed
    in the report and are not silently treated as absent evidence. Selected,
    evidence-available occupancy regions additionally receive very thin
    vertical boundary guides across every coordinate-bearing track, so a saved
    dominant CUT&RUN support window remains visible beside transcripts,
    regulatory annotations, occupancy lanes, and motif scores. GENtle does not
    infer which window is dominant; save the reviewed support window and add
    its set to the composition
  - stable tutorial semantic IDs cover manager open/save/copy/import/export and
    CUT&RUN/Ensembl capture controls for Linux acceptance automation

<!-- manual-page-nav:start -->

---

[GUI Manual](../gui.md) | [Previous: Locus Figure And Evidence Tracks](locus-figure-and-evidence-tracks.md) | [Next: Conservation Workspace](conservation-workspace.md)
<!-- manual-page-nav:end -->
