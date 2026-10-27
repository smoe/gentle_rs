### Feature Tree And Selection

<!-- manual-counterparts:start -->
Related CLI reference: [Feature Annotation](../cli/features.md).
These topic links do not imply identical feature coverage across interfaces.
<!-- manual-counterparts:end -->

Feature tree grouping:

- The left-side feature tree includes a `Grouping` selector:
  - `Off`
  - `Auto (duplicates)`
  - `Always`
- Repeated labels (for example many RNA variants) can be grouped under each
  feature-kind section depending on that mode.
- `Auto (duplicates)` flattens singleton subgroups and only keeps grouped
  branches where duplicate labels exist.
- The feature-tree header includes focus presets:
  - `Cloning`: clears the tree filter and returns grouping to `Auto (duplicates)`
  - `Regulatory`: filters to regulatory annotations and uses `Always` grouping
  - `Repeats`: filters to repeat/rmsk annotations and uses `Always` grouping
  - `Tracks`: filters to imported track/array annotations and uses `Always`
    grouping
- A compact `Layers:` summary shows the current tree/filter mix (`core`,
  `regulatory`, `repeats`, `tracks`, `array`, `TFBS`, and `other` when present).
  In linear mode, those chips use the same `visible/total` convention as group
  headers.
- `gene` entries remain ungrouped (flat rows), even in grouped modes, because
  identifiers are expected to be unique.
- `mRNA` rows are grouped by their associated gene (when gene qualifiers are
  available), while transcript identifiers, spliced nucleotide lengths, and the
  best available transcript/isoform names are shown on individual mRNA entries.
- Regulatory entries (`enhancer`/`silencer`) are grouped in nested branches:
  - primary branch: regulatory class buckets (`enhancer`, `silencer`, `other`)
  - secondary branch:
    - `active region` when labels start with that prefix
    - enhancer marker token groups (for example histone-like markers such as
      `H3K4me1`/`H3K27ac-H3K4me1` and other protein-like tokens such as `CTCF`)
    - unclassified/non-enhancer+silencer regulatory rows remain in `other`
  - per-entry labels remove grouped prefixes to reduce redundancy.
- Group headings show counts as `visible/total` in linear mode (current stretch
  shown in the map), and as total count in circular mode.
- `Cmd` (macOS) / `Ctrl` (Windows/Linux) click on feature rows toggles
  multi-selection.
- In linear mode, multi-selected features are forced to use external labels
  (including short features), improving disambiguation.
- The feature-tree pane is horizontally resizable; when narrowed, the tree uses
  horizontal scrolling for long labels.
- In hosted windows, the feature tree and its split/scroll widgets use
  per-window ids so temporary relayouts such as toggling the sequence text
  panel do not collide with another hosted view of the same sequence.
- Lazy-opened sequence windows defer initial feature-tree/detail rendering to
  improve window-open responsiveness. Documents with at most 300 features load
  that browser automatically on the next frame; annotation-rich documents stay
  on the lightweight placeholder until `Load feature tree/details` is clicked,
  while the DNA map remains usable.
- Once the feature tree is visible, its derived grouping/filter/visibility
  model is cached across frames. Pan/zoom updates only exon-aware visibility
  and group counts, retaining labels and group membership. Sequence/annotation,
  topology, feature visibility settings or filter/grouping changes still rebuild
  the model. An intron-only viewport does not make a transcript visible merely
  because it overlaps that transcript's bounding span.
- Feature detail text remains below the feature tree in the left pane and uses
  the configurable feature-detail font size.
- The feature tree/details pane is top-aligned with the map and stretched to
  fill map-panel height, so no unused bottom gap is left above the sequence
  text panel.
- When multi-selection is active, the feature-tree header shows a
  `Multi-select active (N)` chip and a one-click `Clear multi-select` action.
- Right-clicking a feature row opens per-feature actions; the whole visible
  row, not only the text label, accepts the context menu gesture:
  - `Focus feature (current zoom)`
  - `Fit feature in view` (linear map mode)
  - `Copy identifier`, `Copy description`, and `Copy tooltip text`
  - `Copy spliced mRNA sequence`, `Copy spliced mRNA as FASTA`, and
    `Export spliced mRNA FASTA...` for `mRNA`/`transcript` rows; these use the
    same strand-aware exon assembly as the shared transcript derivation engine
    without adding a derived sequence to the project
  - `Use as promoter anchor (Engine Ops)` for `mRNA`/`transcript` rows
    (seeds anchored extraction defaults with strand-aware boundary/direction)
  - `Open Splicing Window`, `Open RNA-read Mapping`, and `Derive + Dotplot`
    for splicing-compatible feature kinds
    (`mRNA`/`transcript`/`ncRNA`/`misc_RNA`/`exon`/`gene`/`CDS`)
- Single-click on feature rows or feature glyphs now selects/focuses only; it
  no longer auto-opens the Splicing Expert window.
- Open splicing-linked expert windows deliberately via:
  - Agent Assistant or GUI Shell: `ui open splicing-expert SEQ_ID FEATURE_ID`,
    using a loaded sequence and a zero-based ID from `features query`.
    `ui focus` reuses that inspection without clearing evidence; `ui close`
    with the same operands closes only that expert, retaining project data.
    Opening waits for a lazy DNA viewer to finish loading. No active-window
    or canonical-transcript guess is made.
  - double-click a feature row
  - double-click a feature glyph on the DNA map
  - right-click context menu from feature tree or map
  - `Feature Details` pane buttons:
    `Open Splicing Window`, `Open RNA-read Mapping`, `Derive + Dotplot`
  - right-click the `Feature Details` title to reopen those same three
    destinations from the currently selected splicing-compatible feature, or
    to copy the selected feature identifier, description, or complete tooltip
    text
- `Filter` narrows all feature rows in the tree.
  - free text searches kind/label/range and selected qualifiers
  - scoped terms are supported:
    - `kind:mrna`
    - `label:TP73`
    - `range:6128..16430`
    - `source:BED` / `source:VCF`
    - `track:chip`
    - `path:peaks.bed` or `file:peaks.bed`
    - `note:enhancer`
    - `regulatory:enhancer`
    - `repeat_class:LINE` for RepeatMasker/UCSC `rmsk`-style repeat features
  - free-text `regulatory`, `repeat`, and `track` terms match their respective
    dense-layer families for quick focus presets
  - preset terms below the filter box are on/off toggle buttons

<!-- manual-page-nav:start -->

---

[GUI Manual](../gui.md) | [Previous: DNA Letters And Zoom](dna-letters-and-zoom.md) | [Next: Restriction-Site Inspection](restriction-site-inspection.md)
<!-- manual-page-nav:end -->
