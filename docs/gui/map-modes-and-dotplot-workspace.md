### Map Modes And Dotplot Workspace

<!-- manual-counterparts:start -->
Related CLI reference: [Dotplots, Flexibility and Alignment](../cli/dotplot.md).
These topic links do not imply identical feature coverage across interfaces.
<!-- manual-counterparts:end -->

Primary map modes (linear topology):

- `Standard map`
  - regular linear feature-map rendering
- `Splicing map`
  - transcript/exon lane rendering for selected
    `mRNA`/`ncRNA`/`misc_RNA`/`transcript`/`exon`/`gene`/`CDS` features
- `Dotplot map`
  - compact launcher/compute panel in sequence windows
  - the launcher fills the remaining map-pane height and scrolls within that
    pane at small sizes; its controls do not impose a 600-px minimum width
  - selecting `Dotplot map` opens a dedicated standalone `Dotplot` workspace
    window for full controls and rendering
  - opening/focusing a dotplot or flexibility window from an analysis,
    transcript/feature or construct-reasoning action preserves the current
    `Standard map`, `Splicing map` or other primary view; only explicit mode
    selection replaces it with this launcher
  - an already saved `Dotplot map` preference is not reset automatically;
    select `Standard map` to restore the sequence map while keeping Dotplot open
  - compact panel keeps quick operations in-window:
    - `Compute dotplot`
    - `Compute flexibility`
    - `Export Dotplot SVG...`
      - exports the currently loaded dotplot/flexibility rendering through engine
        operation `RenderDotplotSvg` (shared GUI/CLI/scripting parity path)
      - each export is recorded in operation history and appears in lineage
        table/graph as an export analysis node
      - default filename includes dotplot parameters (`mode`, query/reference
        spans, `word`, `step`, `mismatches`, optional `tile`) and display
        controls (`threshold`, `gain`; plus flexibility model/bin/smoothing when
        shown)
      - overlay exports switch to owner/reference-oriented naming and include
        the overlaid series count
      - pairwise, annotated self, and overlay dotplots draw a genome-context
        side rail beside the reference/genomic y-axis when the reference
        sequence already carries supported annotation features
      - the rail is not constrained to a selected gene: exon features keep
        sense/antisense strand cues in separate lanes, and materialized
        RepeatMasker/UCSC `rmsk` repeat features (`repeat_region` /
        `mobile_element`) render as a repeat lane
    - `dotplot_id` / `flex_track_id` selection
  - full parameter editing and plot inspection now live in the standalone window:
    - before a dotplot payload is loaded, the standalone window stays in a
      compact settings-first layout instead of reserving the full plot canvas
      height
    - self modes: `self_forward`, `self_reverse_complement`
    - pair modes: `pair_forward`, `pair_reverse_complement`
    - pair modes require `reference_seq_id`
      - reference controls now include:
        - `Use query as ref`
          - resets the reference to the current query sequence and expands the
            reference span to the full query length for pair-forward self-vs-self
            testing
        - `Use locus DNA`
          - restores the current splicing/RNA-mapping locus DNA sequence as the
            reference and resets `ref_start/ref_end` to the current locus ROI
        - `Annotated mRNA ref...`
          - derives one annotated transcript/mRNA sequence for the active locus
            and uses it as the pairwise reference without switching the active
            DNA/query window
      - changing the reference to a shorter sequence now automatically
        normalizes stale reference spans so old genomic `ref_end` values do not
        survive into RNA-read or transcript-reference computes
      - optional `ref_start` / `ref_end` set y-axis reference span
      - default behavior when `ref_start/ref_end` are empty:
        - run an initial full-reference pass
        - auto-fit to detected hit envelope (+padding)
        - recompute once on that fitted reference span
      - `Fit ref span to hits` remains available for manual re-fit when explicit
        reference spans are used
      - `Overlay transcript isoforms` becomes available when the active
        reference is the current splicing/RNA-mapping locus DNA:
        - multiple transcript lanes can be selected and overlaid against the
          same genomic/reference span in one canvas
        - `overlay x-axis` chooses how transcript queries share the x-axis:
          `% transcript` normalizes each isoform to 0..100%, `bp left-aligned`
          aligns transcript starts to emphasize different 3' splice endings,
          `bp right-aligned` aligns transcript ends to emphasize different
          5' splice starts, `shared exon anchored` shifts only the transcripts
          that contain the chosen shared exon so that exon start is shown at
          one common x-position, and `manual/domain anchor` uses explicit
          per-query anchors already stored in the loaded overlay payload
          (intended for curated cross-family/domain comparisons)
        - in `shared exon anchored` mode an explicit exon picker appears:
          - it lists only exons shared by at least two of the currently
            selected transcript lanes
          - the transcript checkbox state stays unchanged, but only the
            transcripts containing the picked exon are actually drawn/exported
          - the common anchor position is the maximal transcript-local exon
            start observed among the displayed transcripts for that exon
        - each isoform gets a deterministic color while the y-axis stays in
          shared genomic coordinates
        - genome-context annotation for the reference span is drawn beside the
          genomic axis when supported features are present on the reference
          sequence; exons keep sense/antisense lanes and materialized
          RepeatMasker/UCSC repeat features share the same side rail
        - manually anchored overlay payloads can also be re-opened here even if
          they were computed outside the current locus-transcript picker; the
          x-axis selector stays available so those persisted family/domain
          overlays can still be re-rendered/exported from the GUI
    - default parameters:
      - `half_window_bp`: defaults to the larger query/reference sequence length
        so the initial view spans the full comparison context
      - `word`: defaults to `7` (higher values are faster but less sensitive)
      - `step`: defaults to `1` (highest sampling density)
      - `mismatches`: defaults to `0` (exact-seed matching)
    - display controls:
      - `display threshold` (cell-density sensitivity)
      - `intensity gain` (contrast amplification for visible cells)
      - `Auto contrast` (estimates threshold/gain from payload density)
      - display-only controls update immediately; engine recompute is reserved
        for biological/compute input changes
      - if no cells pass visibility, the canvas shows an explicit message instead
        of a silent white panel
    - valid cheap requests now auto-compute after a short debounce; requests
      beyond the pair-evaluation budget stay dirty and wait for explicit
      `Compute dotplot`
    - optional paired flexibility-track panel (`AT richness` / `AT skew`)
    - boxplot summary panel:
      - rendered below the density map
      - each query bin shows reference-hit distribution
        (`min/q1/median/q3/max` + hit count)
      - useful for quickly spotting exon-band structure in cDNA-vs-genomic maps
    - linked crosshair:
      - hover for live `x/y` coordinates in the plotted span
      - click to lock crosshair
        - self modes: sync selection to x/y interval
        - pair modes: sync selection from x/query axis (y/reference stays informational)
      - right-click to clear the locked crosshair
      - overlay mode keeps hover-only crosshairs and reports one shared
        reference coordinate plus one per-isoform query coordinate using the
        selected overlay x-axis layout; query sync and locked crosshair stay
        disabled there because the overlay does not have one unambiguous shared
        query coordinate system
        - in `shared exon anchored` mode, transcripts that do support the
          chosen exon still report `outside` when the cursor lies beyond that
          shifted transcript extent
        - in `manual/domain anchor` mode, the same `outside` readout is used
          for series whose explicit anchored extent does not cover the hovered
          x position
    - status panel:
      - deterministic request diagnostics (mode, spans, seed parameters, estimated
        window counts / pair evaluations)
      - for `mismatches=0`, large pairwise requests use indexed exact-seed
        matching (not brute-force pair loops), so low-step runs can remain practical
      - loaded payload diagnostics (point count, estimated hit fraction, latest
        operation status/messages)
      - sparse pairwise payload guidance:
        - suggests `pair_reverse_complement` when `pair_forward` is sparse
          (common for cDNA vs genomic comparisons)
        - warns when strict settings (`word`/`step`/`mismatches`) are likely too
          aggressive for wide reference spans
        - warns when computed hits sit near reference-span edges and recommends
          ref-span fitting
      - interpretation notes:
        - `self_forward`: dominant diagonal is expected identity signal
        - reverse-complement self-pair zero-hits can be valid under strict/sparse
          sampling (not necessarily an error)
  - main lineage table/analysis details:
    - dotplot analysis rows include `Dotplot SVG` action (analogous to gel export)
      for direct export from the main project window

<!-- manual-page-nav:start -->

---

[GUI Manual](../gui.md) | [Previous: Strand Placement](strand-placement.md) | [Next: DNA Letters And Zoom](dna-letters-and-zoom.md)
<!-- manual-page-nav:end -->
