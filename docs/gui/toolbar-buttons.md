## Toolbar buttons

<!-- manual-counterparts:start -->
No direct CLI counterpart: This chapter describes desktop presentation or interaction; no dedicated command reference is paired with it.
<!-- manual-counterparts:end -->

The top toolbar in each DNA window provides these controls (left to right).
The exact set depends on `Mode`:

- `Genome`
  - read-focused; hides ROI extraction, sequence-derivation, and Engine Ops/Shell buttons
- `Chromosomal`
  - keeps ROI extraction (`Extract Sel`, `PCR ROI`), hides sequence-derivation and Engine Ops/Shell
- `Region`
  - full controls (default)
- `Gene`
  - full controls
- `cDNA`
  - full controls
  - defaults to intrinsic cDNA-focused display behavior:
    contextual transcript-projection features (`mRNA`/`exon`/`CDS`) are hidden
    unless explicitly enabled

Toolbar layout:

- The toolbar is arranged in stable left-anchored rows instead of one large
  wrap-everything flow:
  - map/navigation row
  - display/layer toggles row
    - includes a dedicated `Repeats (N)` toggle backed by shared
      `show_repeat_features` display state; the same setting is honored by
      linear/circular map rendering, feature-tree filtering, declutter restore,
      and SVG export
  - ROI/selection row
  - derivation/export/action row
- The `Selection formula` / `Apply Sel` controls now stay left-aligned on
  their own row instead of floating into the preceding button row when the
  window becomes narrower.
- Each toolbar row wraps to the available window width, including the `Go`
  coordinate fields and formula controls. The toolbar scrolls vertically when
  it reaches its height limit, leaving space for the sequence map.
- High-frequency map/display controls use existing icon assets where the
  meaning is clear, while every icon keeps the same hover/status description
  used for support and debugging. Text labels remain for controls without a
  clear existing icon.

Controls:

1. Circular/Linear toggle
   - Switches DNA map topology visualization between circular and linear.
   - For long circular sequences, switching into linear mode now opens a bounded
     initial linear viewport instead of forcing a full-sequence span immediately,
     to keep the window responsive; use `Fit Seq` if you explicitly want the
     full linear extent.
   - For debugging topology-switch hangs, launching GENtle with
     `GENTLE_TRACE_TOPOLOGY_TRANSITIONS=1` emits coarse stage tracing to the
     terminal and mirrors it into the window status line.
2. Show/Hide sequence panel
   - Shows or hides the sequence text panel.
3. Show/Hide map panel
   - Shows or hides the graphical DNA map.
4. Mode
   - Chooses biological presentation context (`Genome`, `Chromosomal`, `Region`, `Gene`, `cDNA`).
   - Toolbar groups are shown/hidden according to the selected context.
5. Splicing map (linear mode)
   - Toggles primary map rendering between standard feature-map mode and a
     read-only transcript/exon splicing-lane mode for the selected
     `mRNA`/`ncRNA`/`misc_RNA`/`transcript`/`exon`/`gene`/`CDS` feature.
   - Uses the same splicing payload/geometry contract as the dedicated
     Splicing Expert window.
   - Clicking a transcript lane focuses the corresponding transcript feature in
     the sequence view.
6. CDS
   - Shows or hides CDS feature rendering.
7. Gene
   - Shows or hides gene feature rendering.
8. mRNA
   - Shows or hides mRNA feature rendering.
9. Reasoning
   - Toggles the read-only construct-reasoning overlay for the active sequence.
   - Current overlay auto-refreshes the engine-derived graph and paints only
     annotation-grade sequence spans on the linear DNA map:
     - generated/inferred annotation-like spans such as promoter windows
     - cDNA-confirmed exon/splice/CDS/transcript spans
     - sequence intervals that directly support a stored fact or decision
   - Those visible spans are now backed by shared engine-owned annotation
     candidate records instead of adapter-local filtering alone.
   - Raw low-level evidence such as unrestricted restriction-site or TFBS spam
     stays in the graph but is no longer drawn by default just because it
     exists as `sequence_span` evidence.
   - Hovering an evidence span shows a compact "why is this here?" tooltip.
   - Clicking an evidence span routes the side-panel details view to that exact
     annotation candidate / reasoning record and now also shows:
     - whether it is generated, confirmed, or merely supporting
     - which stored facts and decisions it supports
     - transcript-context status and effect tags when that span participates in
       transcript-aware variant reasoning
   - The side-panel description area now also shows a dedicated
     `Construct reasoning` section for non-sequence facts and decision steps
     (for example host/helper/growth/selection context) so those records stay
     inspectable without being faked as coordinate spans.
   - The same side panel now also includes an `Annotation candidates` section
     listing the current graph-backed candidate annotations with their support
     summaries and transcript-aware warnings.
   - The same side panel now also includes an `Annotation summaries` section
     that collapses overlapping candidate spans into calmer locus-level rows
     such as overlapping exon or transcript-aware variant summaries.
     - clicked spans now also show which collapsed summary they belong to,
       including candidate count and review-state summary
     - this summary layer is intended to be the future integration point for
       similarity-derived annotation summaries as well, so large genomic views
       do not regress into raw overlay floods
   - Generated repeat/similarity annotations now also land here in the same
     calmer style:
     - low-complexity regions
     - homopolymer / tandem-repeat segments
     - direct-repeat and inverted-repeat clusters
     - cautious `Alu-like SINE candidate` mobile-element calls
     - operational-risk facts stay in the non-sequence inspector instead of
       duplicating every underlying seed/hit on the map
       - current predictor rows include:
         - `PCR/amplification review suggested`
         - `Nanopore/direct-sequencing review suggested`
         - `Repeat-driven mapping review suggested`
         - `Cloning stability review suggested`
       - task-aware repeat/similarity severity is shown as a compact structured
         grid under the same fact rows. Optional intrinsic evidence
         (`base_severity`/`base_score`) and objective-adjusted effective
         priority appear only when present; task applicability and its
         explicit/legacy/unspecified basis remain separate. This does not
         create extra map overlays
     - repeat/similarity summaries, candidates, and repeat-driven fact rows
       can now open the shared dotplot workspace directly:
       - `Dotplot` opens a self-forward dotplot centered on the implicated
         repeat/low-complexity span
       - `RevComp Dotplot` opens a self-reverse-complement dotplot for
         inverted-repeat / inversion-risk review
       - the linear viewport is focused onto the implicated span first so the
         dotplot uses the same region the reasoning row is describing
       - these buttons are backed by graph-level
         `inspection_actions[]` with deterministic action ids and
         `rationale`, `driving_evidence_ids[]`, mode, context tags, and focus
         ranges, not by GUI-local label matching
       - a collapsed `Evidence` pane beside each button presents the action's
         1-based focus range and length, readable dotplot mode, wrapped
         rationale, discrete evidence ids and context tags, and typed
         repeat-family provenance; the button and explanation therefore stay
         tied to the same protocol record
       - repeat-family rows enumerate every typed overlapping family with
         source, class/family agreement, and confidence; the plural provenance
         field is preferred while older singular provenance remains readable
       - the inspector reports graph snapshot freshness. Dotplot buttons are
         disabled, and the execution path independently rejects the action,
         when sequence/features, objective intent, or reasoning rules changed
         after graph generation; refresh is required before inspection
       - CLI/agent layers can list and run the same action objects with
         `construct-reasoning list-inspection-actions` and
         `construct-reasoning run-inspection-action`, so the GUI does not own a
         separate answer to "which dotplot should I inspect?"
  - Both the clicked-span detail view and the `Annotation candidates` section
    now expose shared-engine accept/reject/draft controls so reviewed
    automated-annotation candidates can be curated without inventing GUI-only
    state.
  - Accepted or locked generated candidates now also expose `Write Back`,
    which materializes the reviewed candidate as an ordinary sequence feature
    through the shared engine path and immediately refreshes the overlay and
    details panel.
  - Accepted candidates stay role-colored but gain a clear green review
    outline, while rejected candidates remain inspectable in a visibly muted
    form instead of disappearing from the graph.
   - That inspector now also surfaces host-route restriction/methylation review
     when route notes encode traits such as `hsdR- M+`, `dam+`, `dcm+`,
     `hsdR+`, or `MDRS+`.
   - Adapter/linker restriction-capture reasoning now also appears there as a
     non-sequence fact: blunt-end capture plans, minimal vs MCS-like adapter
     style, per-motif conflicts across the capture site plus extra retrieval
     sites, the stronger “all named adapter motifs already occur on the insert”
     case, and “planned insert methylation still needs enzyme-specific review”
     warnings stay visible without pretending they are coordinate spans.
   - Variant-aware reasoning now also appears there as non-sequence facts:
     promoter / TFBS / CDS / UTR / splice overlaps are summarized as likely
     effect candidates, and the same shared engine rules surface first assay
     suggestions such as paired promoter-luciferase reporters, allele-paired
     coding-expression comparison, minigene splicing, or UTR reporter assays.
   - The adjacent `Why` menu filters visible reasoning spans by role and
     evidence class without changing the underlying engine graph.
10. Ctx mRNA (shown in `cDNA` mode)
   - Explicit opt-in for contextual transcript-projection features derived from
     genome annotation/mapping context.
   - Default is off in `cDNA` mode (intrinsic cDNA evidence only).
11. Show/Hide TFBS
   - Toggles computed TFBS annotations.
   - Default is off.
12. Show/Hide restriction enzymes
   - Toggles restriction enzyme cut-site markers and labels.
   - Marker color now encodes cleavage geometry:
     - blunt = slate
     - `5'` overhang = blue
     - `3'` overhang = amber
   - The adjacent restriction mode control switches between:
     - `Preferred only`
     - `Preferred + unique`
     - `Unique only`
     - `All in view`
   - Non-`All` modes show `shown/total` counts so it is clear when cutters are
     hidden by the current filter instead of being absent.
   - Preferred enzymes are managed as a comma-separated REBASE-name list and
     include preset buttons that switch the filter to `Preferred only` for the
     canonical pUC MCS set
     (`EcoRI,SacI,KpnI,SmaI,BamHI,XbaI,SalI,PstI,SphI,HindIII`) and a Golden
     Gate Type IIS set (for example `BsaI`, `BsmBI`, `BbsI` when present in the
     active REBASE catalog).
13. Show/Hide GC content
   - Toggles GC-content visualization overlay.
   - Aggregation uses the configurable GC bin size (default `100 bp`).
14. Show/Hide ORFs
   - Toggles open reading frame overlays.
   - Predictions use ATG and the first in-frame TAA/TAG/TGA, on both strands.
     The existing display threshold is 100 codons including the stop (at least
     99 translated residues). These overlays are not annotation-backed CDSs.
   - Circular predictions may cross the origin, including a start or stop
     codon split by that origin, but never traverse more than one molecule
     length. Switching to the linear map splits such an overlay at the origin
     rather than displaying the complementary interval. SVG uses the same spans.
15. Show/Hide methylation sites
   - Toggles methylation-site markers.
16. Extract Sel
   - Extracts current map/text selection into a new sequence via
     `ExtractRegion`.
   - Disabled until a non-empty selection exists.
   - If the current viewport is the intended span, use `Select visible` first
     to promote the visible linear map window into the active selection.
     - When no explicit `output_id` is supplied, the default extracted sequence
       id now reflects the selected local `1-based` interval; dbSNP auto-fetched
       parents collapse back to the rsID stem instead of repeating the full
       fetched-locus coordinates.
17. Select visible
   - Copies the current visible linear map span into the active selection.
   - Intended as a quick bridge between viewport navigation and
     selection-first actions such as `Extract Sel`, `Queue PCR selection`, and
     `PCR ROI`.
18. Queue PCR selection
   - One-click shortcut beside `Extract Sel` that queues the current linear
     selection as one PCR region for later batch primer-pair design.
   - When a selection is present, the map area also shows an inline hint with
     the active span and the same PCR queue action.
19. RE scan
   - Toolbar menu for direct non-mutating restriction-site inspection on the
     current active DNA window.
   - Available entry points:
     - current selection
     - current visible linear span
     - whole active sequence
   - Uses the shared `FindRestrictionSites` engine operation rather than a
     GUI-only quick scan, so CLI/shell/agent calls see the same report schema.
   - Status feedback summarizes the top hits inline.
   - The same window now also caches the shared report in `Direct scan inspectors`,
     where you can inspect the portable hit table and export it again as JSON.
20. TFBS scan
   - Toolbar menu for direct non-mutating TFBS/JASPAR hit inspection on the
     current active DNA window.
   - Available entry points:
     - current selection
     - current visible linear span
     - whole active sequence
   - Reuses the current motif and threshold settings from the `TFBS annotation`
     panel instead of maintaining a second GUI-only parameter set.
   - Uses the shared `ScanTfbsHits` engine operation rather than a GUI-only
     quick scan, so CLI/shell/agent calls see the same report schema.
   - Status feedback summarizes the top hits inline; the quick GUI path caps
     the report at 200 hits to keep inspection responsive.
   - The same window now also caches the shared hit report in `Direct scan
     inspectors`, where you can inspect the portable hit table and export it
     again as JSON.
   - The same menu offers an optional `Precomputed hits` section for 1-64
     explicit JASPAR matrix IDs. It queries a configured finalized
     `jaspar-mapping` sparse package for the current selection, visible span,
     or whole anchored sequence through the shared
     `QueryGenomicMotifEvidence` operation.
   - This query runs in the existing background TFBS worker from an immutable
     engine snapshot. Its report and map strips are transient display state;
     they do not add features or dirty the project.
   - Package discovery is lazy. When
     `GENTLE_JASPAR_GENOME_SCAN_PACKAGE` or DuckDB is absent, the inspector
     reports that optional provider as unavailable while ordinary local TFBS
     hits, score tracks, and the rest of GENtle continue unchanged.
   - Package-native scores and local scores remain separate. The inspector
     shows per-motif retention-floor completeness and warns that a sparse
     package miss is not biological absence.
   - The configured package can also be the distinct **regulatory/TSS subset**.
     The inspector names that scope and exposes virtualized physical-TSS,
     transcript-owner and native regulatory-gene-link details. Empty subset
     coverage does not assert genome-wide motif absence. Use the shared Shell's
     `--inspect`, `--gene` or `--tss-id` options for package discovery/queries
     without loaded DNA; see the [offline walkthrough](../regulatory_motif_subset.md).
     Matching saved queries can be attached to TSS reports and then displayed
     through `Attach the quantitative report`, with a persistent subset label.
     No automatic annotation, rescoring or dedicated package wizard is added.
21. TFBS score tracks
   - Toolbar menu for direct non-mutating continuous TFBS/JASPAR score
     inspection on the current active DNA window.
   - Available entry points:
     - current selection
     - current visible linear span
     - whole active sequence
   - Reuses the same motif selection as `TFBS annotation` while keeping the
     chosen score-family/clipping mode in that panel, so no second GUI-only
     motif configuration is introduced.
   - Uses the shared `SummarizeTfbsScoreTracks` engine operation rather than a
     GUI-only painter, and the cached view can be exported from the same panel
     through `RenderTfbsScoreTracksSvg`.
22. TFBS similarity
   - Toolbar menu for direct non-mutating anchor-vs-candidate TFBS/JASPAR
     score-track ranking on the current active DNA window.
   - Available entry points:
     - current selection
     - current visible linear span
     - whole active sequence
   - Reuses the score-family/clipping mode from the `TFBS annotation` panel
     and keeps one explicit similarity subpanel there for:
     - anchor motif
     - candidate scope (`all JASPAR motifs` or an explicit/current subset)
     - ranking metric
     - optional cached remote-metadata enrichment
     - optional species filters
     - optional output row cap
   - Uses the shared `SummarizeTfbsTrackSimilarity` engine operation rather
     than a GUI-only ranking table, and the cached report can be exported from
     the same panel as JSON.
   - Current semantics are intentionally explicit: species filters restrict the
     candidate set before ranking; a future post-rank visibility filter is
     still tracked separately.
23. Selection formula
   - Toolbar input for formula-driven selection ranges in the DNA window.
   - Drag-selecting on the linear map now also updates this field to the
     current `=start .. end_exclusive` range, and existing map selections can
     be adjusted by dragging their left or right edge handles directly in the
     map.
   - Formula must start with `=` and define a range:
     - `=left .. right`
     - `=left to right`
   - Coordinate terms support:
     - numeric bp positions
     - feature-relative terms:
       - `KIND.start|end|middle|tss`
       - strand-aware upstream boundary (`KIND.upstream(N)`)
       - optional occurrence selector (`KIND[2]`, 1-based)
       - optional label selector (`KIND[label=TP73]`)
       - optional signed offsets (`+N`, `-N`)
   - Range formulas select the interval between the two resolved coordinates;
     this keeps `=gene.upstream(1000) .. gene.tss` usable on reverse-strand
     annotations where the upstream coordinate is numerically greater than the
     TSS.
   - `Apply Sel` resolves the formula and sets the active map/text selection,
     which can then be used directly by `Extract Sel`, `Queue PCR selection`,
     or `PCR ROI` actions.
   - To inspect that selection, use `Zoom to selection` in the navigation row
     or the map context menu. It fits the entire selected interval, even if it
     was off-screen. `Go to selection` centres the same interval without
     changing the zoom level. Neither action reapplies an edited formula:
     press `Apply Sel` first. Invalid formulas leave the previous selection
     unchanged and do not move the view.
   - The same `Selection formula` + `Apply Sel` control is also shown in the
     dedicated `PCR Designer` specialist window so pair-PCR setup can stay
     selection-first without switching back to the sequence toolbar.
20. PCR ROI
   - The DNA-window toolbar now also exposes a direct `Set PCR ROI` button for
     either the current map/text selection or the current painted ROI interval,
     alongside `Queue PCR selection`.
   - `PCR ROI` menu supports both single-ROI seeding and batch-queue capture:
     - seed Primer/qPCR ROI from current map/text selection or painted ROI
     - add current map/text selection to PCR region queue
     - seed Primer/qPCR ROI from selected feature bounds
     - add selected feature(s) to PCR region queue (one queued row per feature)
   - map context menu now also offers `Set PCR ROI from selection / painted ROI`
     for the active linear-map context.
   - map context menu now also offers `Simple PCR from selection` whenever a
     non-empty selection is active:
     - seeds the selected span as the core ROI
     - enables ROI flanking
     - applies simple left/right flank windows from the current beginner
       distance setting
     - opens or focuses the dedicated PCR Designer
   - when a feature is under the cursor or selected, the map context menu also
     exposes the same `Copy identifier`, `Copy description`, and
     `Copy tooltip text` actions as the feature tree.
- direct ROI seeding now opens or focuses the dedicated `PCR Designer`
  specialist so pair-PCR and qPCR setup can stay out of the crowded lower
  `Engine Ops` pane.
- the dedicated `PCR Designer` now has an explicit
  `Pair PCR | qPCR | Transcript panels | RT primer pool | Gene assay study` mode switch; ROI
  seeding is shared for the first two modes, transcript panels use a complete
  annotated Splicing Expert group, and RT-primer pools use explicitly ordered
  transcript targets.
   - after paint-dragging on the linear map, the post-drag chip now includes
     direct coordinate editing (`start..end`) for the painted interval, with
     explicit apply action (0-based, end-exclusive).
21. Export Seq
   - Exports the active sequence via engine `SaveFile`.
   - Output format is inferred from filename extension (`.gb/.gbk` => GenBank, `.embl/.emb` => EMBL, `.fa/.fasta` => FASTA).
   - GenBank and EMBL preserve the same annotated nucleotide record; FASTA is
     sequence-only. See [format parity and limits](../sequence_format_parity.md).
22. Export SVG
   - Exports the active sequence map via engine `RenderSequenceSvg`.
   - In linear mode, the SVG now honors the active bp viewport/crop instead of
     always exporting the full sequence length.
   - Single-base `variation` features are exported as direct DNA-baseline
     markers instead of generic feature blocks on separate lanes.
   - Linear exports now also add transcription-start hooked arrows for
     strand-bearing `gene`/`mRNA`/`CDS`/`promoter` features and suppress
     unlabeled fallback coordinate text that would otherwise clutter
     figure-oriented exports.
   - Accession-only transcript labels are demoted in favor of gene-style labels
     when possible, and nearby repeated non-gene labels are compacted so locus
     exports read more like figures than raw annotation dumps.
   - `mRNA`/`promoter` bars now render with pointed ends and use a short
     hooked transcription-start arrow so strand direction reads from shape as
     well as color/position.
   - In circular mode, the SVG now uses a transparent canvas, renders
     single-base `variation` features as explicit radial markers on the DNA
     ring, and marks transcription starts for strand-bearing
     `gene`/`mRNA`/`CDS`/`promoter` features with a short arrow shaft plus
     direction arrowhead. Circular figure exports also use a slightly larger
     ring and larger label fonts for readability.
   - Circular feature labels now anchor outward by circle side (`start` on the
     right, `end` on the left when possible) and nudge away from occupied
     restriction-label space rather than always sitting centered over the same
     radius.
   - When a `gene` and `CDS` describe the same named plasmid element, circular
     export collapses them into one bench-facing construct part; if the `CDS`
     is a contained coding subrange, export keeps the CDS-aligned geometry so
     promoter-plus-coding cassettes stay readable.
   - Circular cloning-map exports also use the same blue construct-part palette
     for `gene` and `CDS` tracks so plasmid parts read as one family unless a
     separate regulatory/variation color is needed.
   - Regulatory labels prefer interpretable named context (`standard_name` or
     `gene` plus regulatory class, for example `tac promoter` or
     `bla promoter`) when the annotation provides it.
   - RepeatMasker/UCSC `rmsk`-style `repeat_region` / `mobile_element`
     annotations with `repName` / `repClass` / `repFamily` (or
     `rmsk_*` / `repeat_*` / `rpt_*` aliases) use repeat-class labels,
     feature-tree grouping, and subtype colors instead of a single generic
     repeat style.
   - Generated UCSC-rmsk features materialized from the prepared interval
     index carry `gentle_generated=ucsc_rmsk` and remain ordinary
     `repeat_region` features, so the Repeat layer toggle controls them
     without adding a separate adapter-only track.
   - When imported sequence metadata provides a clearer title than the raw
     locus-like name, circular export prefers a bench-facing display title
     derived from `definition` plus accession/version when available.
   - `misc_feature` note text can also supply compact circular labels for
     common annotated vector landmarks such as `MCS` or `factor Xa`.
   - Circular cloning-map exports also add a compact narrative legend when
     annotation provides enough context, using fuller part names plus a short
     functional clause beside the figure instead of crowding that prose onto
     the ring itself.
   - Protein-facing `misc_feature` annotations such as a `factor Xa`
     recognition site now use a host-feature leader so the label reads as part
     of the expressed construct rather than as an unlabeled DNA-side tick.
21. Export View SVG
   - Exports the currently shown sequence-window view composition as SVG
     using the default `screen` profile (map panel + sequence panel extract).
   - When linear `Splicing map` mode is active, export uses the same splicing
     payload renderer as the on-screen splicing map (including transcript/exon
     lane context when an `mRNA`/`exon` is selected).
   - `View SVG ▾` exposes additional export profiles:
     - `wide-context`: larger canvas + expanded base-span context around the
       current linear viewport.
     - `print-a3`: print-oriented A3 landscape SVG with expanded context and
       physical-size metadata (`420mm x 297mm`) for direct print workflows.
   - SVG content and header layout do not depend on the Cargo build profile;
     internal `dev` packages do not insert debugging text or an extra header row.
22. Export RNA SVG (ssRNA only)
   - Exports RNA secondary-structure SVG via shared engine operation `RenderRnaStructureSvg`.
   - Shown only when active sequence is single-stranded RNA (`molecule_type` `RNA`/`ssRNA`).
23. Engine Ops
   - Shows/hides strict operation controls for explicit engine workflows.
   - Digest quick-fill actions:
     - `Use MCS enzymes`: fills digest enzyme list from MCS feature context
       (normalized to REBASE canonical names).
     - `Single-cutters (REBASE)`: fills digest list with enzymes cutting exactly
       once in the active sequence.
     - `Single-cutters in CDS`: same single-cutter filter constrained to
       cleavage positions inside CDS features.
24. Shell
   - Shows/hides the in-window GENtle Shell panel.
   - Uses the same shared command parser/executor as `gentle_cli shell`.

Hovering any button shows a tooltip in the UI.

<!-- manual-page-nav:start -->

---

[GUI Manual](../gui.md) | [Previous: Circular Map Labels](circular-map-labels.md) | [Next: Map interactions](map-interactions.md)
<!-- manual-page-nav:end -->
