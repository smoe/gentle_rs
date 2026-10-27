# Figure and Layout Exports

<!-- manual-counterparts:start -->
Related GUI reference: [Linear map conventions](../gui/linear-map-conventions.md) | [Circular map conventions](../gui/circular-map-conventions.md) | [Locus Figure And Evidence Tracks](../gui/locus-figure-and-evidence-tracks.md).
These topic links do not imply identical feature coverage across interfaces.
<!-- manual-counterparts:end -->

Feature-expert SVG parity status:

- `gentle_cli`: supported via direct commands and shared-shell command path
  (`inspect-feature-expert`, `render-feature-expert-svg`)
- `gentle_js`: supported via `render_feature_expert_svg(...)` wrapper over
  shared engine operation `RenderFeatureExpertSvg`
- `gentle_lua`: supported via `render_feature_expert_svg(...)` wrapper over
  shared engine operation `RenderFeatureExpertSvg`
- splicing expert SVG now includes:
  - UniProt-reference status/provenance in transcript-label hover titles, also
    exposed as `matrix_rows[].uniprot_reference` by `inspect-feature-expert`.
    Only an exact normalized transcript xref or an exact stable Ensembl gene-ID
    xref makes a loaded entry relevant. Names, labels, lengths and projections
    do not. No exact relation means `not_evaluated`, not absence; a gene-only
    match never establishes a transcript xref. Load evidence with
    `uniprot import-swissprot PATH [--entry-id ID]` or an explicitly requested
    `uniprot fetch` before inspecting; rendering itself does not access the network.
  - explicit junction transition support labels/table
  - display-only UniProt groups shared with the GUI (target first per group),
    with one accession/review badge per exact-xref lane and full source hovers.
    Missing/duplicate joins remain unevaluated; saved orphan matrix rows stay
    listed. Engine JSON order, fingerprints and command routes are unchanged.
  - lane boundary ticks and exceptional markers are retained; full motif words
    appear once per oriented boundary/pair in a shared summary. Repeated exon
    support labels are replaced by matrix counts/hovers. Narrow E# headers and
    arc numbers are deterministically suppressed on collision, never shrunk
    below 8 px; the target exon IDs appear only where they fit.
  - frequency-encoded transcript-vs-exon matrix cell coloring
  - predicted exon-to-exon transition matrix with frequency-encoded cell
    coloring
  - CDS flank phase coloring (`0/1/2`) on exon left/right edges when
    transcript `cds_ranges_1based` are available
  - exon-length modulo-3 (`len%3`) cues on exon headers (heuristic frame cue)
  - shared output semantics across CLI/JS/Lua because all routes call the same
    engine renderer


- Shared shell:
  - Commands:
    - `render-svg SEQ_ID linear|circular OUTPUT.svg`
    - `render-dotplot-svg SEQ_ID DOTPLOT_ID OUTPUT.svg [--flex-track ID] [--display-threshold N] [--intensity-gain N] [--overlay-x-axis percent_length|left_aligned_bp|right_aligned_bp|shared_exon_anchor|query_anchor_bp] [--overlay-anchor-exon START..END]`
    - `dotplot render-svg SEQ_ID DOTPLOT_ID OUTPUT.svg [--flex-track ID] [--display-threshold N] [--intensity-gain N] [--overlay-x-axis percent_length|left_aligned_bp|right_aligned_bp|shared_exon_anchor|query_anchor_bp] [--overlay-anchor-exon START..END]`
    - `dotplot overlay-compute OWNER_SEQ_ID [--reference-seq REF_SEQ_ID] --query-spec JSON_OR_@FILE [--query-spec JSON_OR_@FILE ...] [--ref-start N] [--ref-end N] [--word-size N] [--step N] [--max-mismatches N] [--tile-bp N] [--id DOTPLOT_ID]`
    - `render-rna-svg SEQ_ID OUTPUT.svg`
    - `rna-info SEQ_ID`
    - `render-lineage-svg OUTPUT.svg`
    - `protocol-cartoon list`
    - `protocol-cartoon render-svg PROTOCOL_ID OUTPUT.svg`
    - `protocol-cartoon render-template-svg TEMPLATE.json OUTPUT.svg`
    - `protocol-cartoon template-validate TEMPLATE.json`
    - `protocol-cartoon render-with-bindings TEMPLATE.json BINDINGS.json OUTPUT.svg`
    - `protocol-cartoon template-export PROTOCOL_ID OUTPUT.json`

Rendering export commands:

- `svg-png INPUT.svg OUTPUT.png [--scale N] [--drop-dotplot-metadata]`
  - Calls the shared deterministic `resvg` raster helper reused by the
    ClawBio/OpenClaw PNG-first graphics bundling path.
  - Intended for headless post-processing of existing SVG exports; it does not
    add new biology/render semantics by itself.
  - The JSON summary includes `font_face_count`.
  - SVGs containing `<text>` now fail early when no usable font faces are
    available, instead of silently producing label-free PNGs. Minimal
    containers should install a font package such as `fonts-dejavu-core` or
    `fonts-liberation`, or set `GENTLE_SVG_FONT_FILE` /
    `GENTLE_SVG_FONT_DIR` to readable TTF/OTF assets. Optional
    `GENTLE_SVG_MONOSPACE_FAMILY`, `GENTLE_SVG_SANS_SERIF_FAMILY`, and
    `GENTLE_SVG_SERIF_FAMILY` override generic family selection.
- `svg-pdf INPUT.svg OUTPUT.pdf [--scale N] [--drop-dotplot-metadata]`
  - Uses the same deterministic `resvg` path as `svg-png`, then embeds the
    rendered RGB image into a one-page PDF.
  - Intended for documentation and bench-facing handoff bundles that need a
    PDF copy of an existing SVG export without adding new drawing semantics.
- `svg-pdf-set OUTPUT.pdf INPUT.svg [INPUT.svg ...] [--scale N] [--drop-dotplot-metadata]`
  - Renders each ordered SVG as one page in a single lossless raster-backed PDF.
  - Pages retain their own dimensions; locus-context and detailed profile pages
    are not concatenated onto an unsafe extra-tall canvas.
  - The JSON summary records every input page, dimensions, available font-face
    count and actual used-font identities/hashes from that page's rasterization.
    It also declares lossless PDF image encoding. Failed font auditing prevents
    PDF publication; SVGs are never silently rendered through an unaudited fallback.
- `svg-vector-pdf-set OUTPUT.pdf INPUT.svg [INPUT.svg ...] [--scale N] [--drop-dotplot-metadata]`
  - Preserves each ordered SVG as a separate vector page with its own dimensions
    and embeds selectable, Unicode-mapped subset fonts.
  - The JSON summary records exact page geometry and used-font identities.
    This is a static PDF: SVG hover and URI interaction are not preserved.
    Unsupported SVG effects may be rasterized locally rather than converting
    the entire page into an image.
- `render-svg SEQ_ID linear|circular OUTPUT.svg`
  - Calls engine operation `RenderSequenceSvg`.
  - Linear exports honor the current stored linear viewport when one is set,
    so workflow-driven `SetLinearViewport` crops carry through into the SVG.
  - Single-base `variation` features render as DNA-baseline markers, which is
    useful for dbSNP-driven context figures.
  - Linear exports now also add transcription-start hooked arrows for
    strand-bearing `gene`/`mRNA`/`CDS`/`promoter` features and suppress
    unlabeled fallback coordinate labels that would otherwise clutter figure
    exports.
  - When transcript/gene features only carry accession-like ids, exports now
    prefer gene-style labels when available and suppress redundant nearby
    duplicates for cleaner communication figures.
  - Direction-bearing `mRNA`/`promoter` bars now render with arrowed ends, and
    the linear transcription-start cue uses a short hooked arrow rather than
    only a flat triangle marker.
  - Circular exports now use a transparent canvas, show single-base
    `variation` features as radial DNA-ring markers, and add
    transcription-start arrow shafts/arrowheads for strand-bearing
    `gene`/`mRNA`/`CDS`/`promoter` features. Circular figure exports also use
    a slightly larger ring and larger label fonts for readability.
- `render-dotplot-svg SEQ_ID DOTPLOT_ID OUTPUT.svg [--flex-track ID] [--display-threshold N] [--intensity-gain N] [--overlay-x-axis percent_length|left_aligned_bp|right_aligned_bp|shared_exon_anchor] [--overlay-anchor-exon START..END]`
  - Calls engine operation `RenderDotplotSvg`.
  - `DOTPLOT_ID` must exist in stored dotplot payloads (`dotplot compute ...` / GUI compute).
  - `SEQ_ID` is the owner sequence id of the stored payload; for reference-centered
    overlay payloads this is the shared reference/locus DNA sequence, not one
    particular query isoform.
  - `--flex-track` optionally overlays one stored flexibility track in the same SVG.
    Overlay dotplots suppress that panel because each query series uses its own
    normalized x-axis.
  - `--display-threshold` and `--intensity-gain` apply the same density/contrast controls as GUI dotplot display.
  - Annotated pairwise, annotated self, and overlay payload exports render a
    genome-context side rail next to the reference/genomic y-axis when the
    reference sequence already carries supported features; exon blocks keep
    sense/antisense strand cues in separate lanes, and materialized
    RepeatMasker/UCSC `rmsk` repeat features render in a repeat lane.
  - Dotplot rendering does not query a RepeatMasker index on its own; use
    `features materialize-repeats` first when the repeat context should travel
    with the project state and exported SVG.
  - `shared_exon_anchor` requires `--overlay-anchor-exon START..END` using the
    stored reference-exon coordinates from the overlay payload.
  - `query_anchor_bp` uses explicit `query_anchor_0based` values already stored
    in the overlay payload. This is the shared route for curated cross-family
    or domain-anchor figures; it does not require `--overlay-anchor-exon`.
  - Alias: `dotplot render-svg SEQ_ID DOTPLOT_ID OUTPUT.svg [--flex-track ID] [--display-threshold N] [--intensity-gain N] [--overlay-x-axis percent_length|left_aligned_bp|right_aligned_bp|shared_exon_anchor|query_anchor_bp] [--overlay-anchor-exon START..END]`.
- `dotplot overlay-compute OWNER_SEQ_ID [--reference-seq REF_SEQ_ID] --query-spec JSON_OR_@FILE [--query-spec JSON_OR_@FILE ...] [--ref-start N] [--ref-end N] [--word-size N] [--step N] [--max-mismatches N] [--tile-bp N] [--id DOTPLOT_ID]`
  - Calls engine operation `ComputeDotplotOverlay`.
  - `OWNER_SEQ_ID` identifies the stored payload owner and is also reused as
    the shared reference when `--reference-seq` is omitted.
  - Each `--query-spec` must deserialize to one `DotplotOverlayQuerySpec`
    record, for example:
    `{"seq_id":"tp53_iso_1","label":"Isoform 1","mode":"pair_forward","color_rgb":[29,78,216]}`
  - Query specs can optionally narrow one isoform/query span via
    `span_start_0based` / `span_end_0based` and switch per-series orientation
    with `mode=pair_forward|pair_reverse_complement`.
  - Cross-family/domain-anchor overlays can also carry
    `query_anchor_0based` plus optional `query_anchor_label`, for example:
    `{"seq_id":"tp63_family_query","label":"TP63","query_anchor_0based":1044,"query_anchor_label":"shared core motif CATGTGTAACAG","mode":"pair_forward","color_rgb":[220,38,38]}`
  - `--ref-start` / `--ref-end` constrain the shared reference span on the
    y-axis for every stored series.
- `render-rna-svg SEQ_ID OUTPUT.svg`
  - Calls engine operation `RenderRnaStructureSvg`.
  - Accepts only single-stranded RNA (`molecule_type` of `RNA`/`ssRNA`).
  - Uses external `RNAfold` for folding and `rnapkin` for drawing (set `GENTLE_RNAFOLD_BIN` / `GENTLE_RNAPKIN_BIN` to override executable paths).
- `render-lineage-svg OUTPUT.svg`
  - Calls engine operation `RenderLineageSvg`.
- `protocol-cartoon list`
  - Lists deterministic protocol-cartoon IDs currently supported by the engine.
- `protocol-cartoon render-svg PROTOCOL_ID OUTPUT.svg`
  - Calls engine operation `RenderProtocolCartoonSvg`.
  - Built-in protocol IDs currently include `gibson.two_fragment`,
    `gibson.single_insert_dual_junction`, `pcr.assay.pair`,
    `pcr.assay.pair.no_product`, `pcr.assay.pair.with_tail`,
    `pcr.oe.substitution`, and `pcr.assay.qpcr`.
- `protocol-cartoon render-template-svg TEMPLATE.json OUTPUT.svg`
  - Calls engine operation `RenderProtocolCartoonTemplateSvg`.
  - Loads template JSON schema `gentle.protocol_cartoon_template.v1` and resolves
    sparse defaults before rendering.
- `protocol-cartoon template-validate TEMPLATE.json`
  - Calls engine operation `ValidateProtocolCartoonTemplate`.
  - Validates template JSON schema `gentle.protocol_cartoon_template.v1` and resolved cartoon semantics without writing SVG output.
- `protocol-cartoon render-with-bindings TEMPLATE.json BINDINGS.json OUTPUT.svg`
  - Calls engine operation `RenderProtocolCartoonTemplateWithBindingsSvg`.
  - Applies binding JSON schema `gentle.protocol_cartoon_template_bindings.v1` to the template before resolving and rendering.
- `protocol-cartoon template-export PROTOCOL_ID OUTPUT.json`
  - Calls engine operation `ExportProtocolCartoonTemplateJson`.
  - Exports canonical built-in template JSON (`gentle.protocol_cartoon_template.v1`)
    for the requested protocol ID.
  - Protocol-cartoon commands intentionally use one canonical namespace with no
    extra alias routes, which keeps scripted and AI-facing command selection
    unambiguous.
- `gibson preview PLAN_JSON_OR_@FILE [--output OUTPUT.json]`
  - Calls the shared non-mutating Gibson preview derivation path.
  - Accepts the single-insert destination-first subset of
    `gentle.gibson_assembly_plan.v1`.
  - Optional plan input:
    `validation_policy.desired_unique_restriction_site_enzyme_name`
    requests one new unique REBASE site on a terminal overlap when the current
    defined-site single-insert plan can support it.
  - The plan payload is resolved against the current project state, so the
    referenced `seq_id` values must already exist in memory or in the loaded
    state file.
  - Returns one `gentle.gibson_assembly_preview.v1` payload with:
    - resolved opening/junction overlaps
    - Gibson primer suggestions (`5' overlap + 3' priming`)
    - validation/advisory findings
    - structured suggested design adjustments for the current priming-window
      blocker when the overlap side already resolves
    - optional structured requested-unique-site outcome
    - protocol-cartoon bindings for the factual Gibson strip
- `gibson apply PLAN_JSON_OR_@FILE`
  - Calls engine operation `ApplyGibsonAssemblyPlan`.
  - Consumes the same plan JSON accepted by `gibson preview`.
  - Creates new project sequences for the Gibson outputs:
    - left insert primer
    - right insert primer
    - assembled product
  - Also records one serial arrangement for gel-oriented review:
    - original vector lane
    - ordered insert lane(s)
    - assembled product lane
    - recommended DNA ladders flanking the samples on export
  - Transfers destination and insert features onto the assembled product
    deterministically.
  - Partially consumed destination annotations are trimmed or projected when a
    truthful product-side rewrite is available.
  - Surviving MCS annotations are revalidated against actual restriction sites
    on the assembled product, including new sites introduced by the insert.
  - Returns the standard mutating `OpResult`, including created sequence IDs,
    warnings, messages, and the operation ID used for provenance/reopen flows.
- `cache inspect [--references|--helpers|--both] [--cache-dir PATH ...]`
  - Calls the shared non-mutating prepared-cache inspection path.
  - Returns one `gentle.prepared_cache_inspection.v1` payload with:
    - selected cache roots
    - manifest-backed prepared installs
    - orphaned remnants under those roots
    - per-entry artifact groups, byte totals, and file counts
  - If no `--cache-dir` is given, default roots are:
    - `--references`: `data/genomes`
    - `--helpers`: `data/helper_genomes`
    - `--both`: both roots above
- `cache clear blast-db-only|derived-indexes-only|selected-prepared|all-prepared-in-cache [--references|--helpers|--both] [--cache-dir PATH ...] [--prepared-id ID ...] [--prepared-path PATH ...] [--include-orphans]`
  - Calls the shared prepared-cache cleanup path.
  - Returns one `gentle.prepared_cache_cleanup.v1` payload with exact affected
    paths, removed byte totals, and per-entry cleanup results.
  - `blast-db-only` and `derived-indexes-only` apply only to manifest-backed
    prepared installs and require one or more `--prepared-id` or
    `--prepared-path` selectors.
  - `selected-prepared` removes only the selected prepared installs and may
    also remove orphaned remnants when `--include-orphans` is set.
  - `--prepared-path` is the precise selector when duplicate prepared ids exist
    across multiple selected cache roots.
  - `all-prepared-in-cache` clears all prepared installs under the selected
    roots and ignores `--prepared-id` / `--prepared-path`.
  - Cleanup is conservative by design:
    - only selected known cache roots are touched
    - catalog JSON and project state files remain unchanged
    - `.gentle_state.json`, MCP/runtime files, backdrop/runtime caches, and
      `target/` are out of scope
- `render-pool-gel-svg IDS|'-' OUTPUT.svg [--ladders NAME[,NAME]] [--containers ID[,ID]] [--arrangement ARR_ID] [--agarose-pct FLOAT] [--buffer tae|tbe] [--topology-aware true|false] [--lane-label-layout auto|horizontal|wrapped|staggered|angled] [--band-label-layout auto|inline|panel] [--isoform-markers auto|off]`
  - Calls engine operation `RenderPoolGelSvg`.
  - Use `IDS` as a comma-separated sequence-id list, or pass `-`/`_` when using `--containers` or `--arrangement`.
  - `--containers` renders one lane per container ID.
  - `--arrangement` renders lanes from a stored serial arrangement in its
    declared lane and role-label order. Plate arrangements are rejected with a
    serial-only readiness error.
  - Collection lifting classifies an explicit container as `combine` and a
    serial arrangement as `arrange`. A project sequence selection is not
    automatically lifted because the direct `IDS` form deliberately
    co-migrates every supplied sequence in one lane; invoking that form is the
    caller's explicit pooling choice.
  - Optional `--ladders` overrides both auto mode and any saved arrangement
    ladder pair.
  - If `--ladders` is omitted for `--arrangement`, engine uses the saved
    arrangement ladder pair when present, otherwise falls back to auto ladder
    selection.
  - If `--ladders` is omitted without `--arrangement`, engine auto-selects one
    or two ladders based on pool bp range.
  - `--agarose-pct` / `--buffer` / `--topology-aware` override the shared
    deterministic gel-run profile for that render.
  - Defaults are `1.0%`, `TAE`, topology-aware `true`.
  - `--lane-label-layout auto` is the default. It leaves short labels
    horizontal, wraps isolated long labels over multiple lines, and angles
    labels when adjacent long names cannot be wrapped without collisions.
  - `horizontal`, `wrapped`, `staggered`, and `angled` force one deterministic
    placement strategy. The SVG grows downward when wrapped or angled labels
    need more room; gel and band coordinates do not change.
  - `--band-label-layout auto` draws an in-gel size annotation only when the
    complete text fits before the next lane or gel edge. `panel` keeps all band
    text in the fragment table, while `inline` restores the legacy always-inline
    behavior for callers that explicitly want it.
  - `--isoform-markers auto` (default) recognizes stable Ensembl transcript
    accessions (`ENST...` and species-specific `ENS...T...`) and RefSeq
    transcript accessions (`NM_`, `NR_`, `XM_`, `XR_`) in product identifiers.
    Version suffixes are ignored when matching the same isoform.
  - Every recognized isoform receives a deterministic color and relative square
    position that repeat on every band containing it. Merged bands may therefore
    carry several markers. The side legend and fragment rows also show a sorted,
    zero-based binary code using `O=0` and `I=1`, so identification does not rely
    on color perception. Use `--isoform-markers off` for a plain gel.
  - when topology-aware mode is on, explicit sequence hints like
    `supercoiled`, `relaxed circular`, or `nicked/open circular` refine
    apparent migration beyond the old generic circular-vs-linear split.
- `render-gel-svg ...`
  - Alias for `render-pool-gel-svg ...` with identical semantics.
- `arrange-serial CONTAINER_IDS [--id ARR_ID] [--name TEXT] [--ladders NAME[,NAME]]`
  - Calls engine operation `CreateArrangementSerial`.
  - Persists a serial lane setup that can be reused by `render-pool-gel-svg --arrangement`.
  - For a resolved gene set, use `op` with `CreateGeneSetArrangement`, explicit
    `resolution_id`, `member_bindings`, fresh `arrangement_id` and
    `confirm_physical_placement=true`; see [tutorial 03.10](../tutorial/03-10_gene_set_storage_projection.md).
    This binds 1..24 exclusive singleton tubes without pooling or cloning and
    creates the existing default rack draft. No new dedicated storage subcommand
    or curated MCP tool is added; physical stock is not inferred.
- `arrange-set-ladders ARR_ID [--ladders NAME[,NAME]]`
  - Calls engine operation `SetArrangementLadders`.
  - Persists one or two ladder names on an existing serial arrangement.
  - Omit `--ladders` to clear back to shared engine auto ladder selection.
- `containers set-exclusive CONTAINER_ID true|false`
  - Calls engine operation `SetContainerDeclaredContentsExclusive`.
  - `true` restores a clean-vial / declared-only interpretation.
  - `false` marks the listed members as a known/measured subset of a more
    complex sample.
- `racks create-from-arrangement ARR_ID [--rack-id ID] [--name TEXT] [--profile small_tube_4x6|plate_6|plate_96|plate_384]`
  - Calls engine operation `CreateRackFromArrangement`.
  - Creates one deterministic physical rack draft from one stored arrangement.
  - If `--profile` is omitted, engine chooses the smallest built-in profile
    that fits the arrangement payload plus ladder-reference slots.
- `racks place-arrangement ARR_ID --rack RACK_ID`
  - Calls engine operation `PlaceArrangementOnRack`.
  - Appends the arrangement as one contiguous block onto the next free region
    in the target rack.
- `racks move RACK_ID --from A1 --to B1 [--block]`
  - Calls engine operation `MoveRackPlacement`.
  - Without `--block`, moves one sample within its arrangement block using
    shift-neighbor semantics.
  - With `--block`, moves the whole arrangement block and shifts later
    occupied positions in fill order.
- `racks move-samples RACK_ID --from A1 [--from A2 ...] --to B1`
  - Calls engine operation `MoveRackSamples`.
  - Moves multiple selected samples together as one contiguous group within a
    single arrangement block.
  - Shared engine preserves the selected samples in current rack order even if
    the CLI request lists them in another order.
- `racks move-blocks RACK_ID --arrangement ARR_ID [--arrangement ARR_ID ...] --to B1`
  - Calls engine operation `MoveRackArrangementBlocks`.
  - Moves the selected arrangement blocks together as one contiguous group.
  - Shared engine preserves the selected blocks in current rack order even if
    the CLI request lists them in another order.
- `racks show RACK_ID`
  - Returns `gentle.rack_state.v1` JSON with the saved rack record and resolved
    occupied positions.
- `racks labels-svg RACK_ID OUTPUT.svg [--arrangement ARR_ID] [--preset compact_cards|print_a4|wide_cards]`
  - Calls engine operation `ExportRackLabelsSvg`.
  - Writes one deterministic label sheet for the whole rack or one selected
    arrangement block on that rack.
  - Built-in presets:
    - `compact_cards` (current compact two-column cards)
    - `print_a4` (denser print-oriented A4 sheet)
    - `wide_cards` (single-column wider cards)
- `racks fabrication-svg RACK_ID OUTPUT.svg [--template storage_pcr_tube_rack|pipetting_pcr_tube_rack|cell_culture_plate|safe_lock_tube_rack_collar_seated|safe_lock_tube_rack_floor_seated|pcr_tube_96_position_carrier|safe_lock_storage_box_derived]`
  - Calls engine operation `ExportRackFabricationSvg`.
  - Writes one top-view fabrication/planning SVG for the full saved rack.
  - Current built-in physical templates:
    - `storage_pcr_tube_rack`
    - `pipetting_pcr_tube_rack`
    - `cell_culture_plate`
    - `safe_lock_tube_rack_collar_seated`
    - `safe_lock_tube_rack_floor_seated`
    - `pcr_tube_96_position_carrier`
- `racks isometric-svg RACK_ID OUTPUT.svg [--template storage_pcr_tube_rack|pipetting_pcr_tube_rack|cell_culture_plate|safe_lock_tube_rack_collar_seated|safe_lock_tube_rack_floor_seated|pcr_tube_96_position_carrier|safe_lock_storage_box_derived]`
  - Calls engine operation `ExportRackIsometricSvg`.
  - Writes one pseudo-3D/isometric SVG for the full saved rack.
  - Intended for presentation-ready rack review and README/documentation
    figures while still using the same saved rack/template state as the other
    physical exports.
- `racks hero-svg RACK_ID OUTPUT.svg [--template storage_pcr_tube_rack|pipetting_pcr_tube_rack|cell_culture_plate|safe_lock_tube_rack_collar_seated|safe_lock_tube_rack_floor_seated|pcr_tube_96_position_carrier|safe_lock_storage_box_derived]`
  - Calls engine operation `ExportRackHeroSvg`.
  - Writes one README-facing top-down rack/plate SVG with a clipped upper-left
    orientation corner, row/column labels, empty circular slots, and
    saved-arrangement labels/rings when slots are occupied.
  - Uses the same saved rack/template state as the technical isometric,
    fabrication, simulation, and OpenSCAD exports.
- `racks openscad RACK_ID OUTPUT.scad [--template storage_pcr_tube_rack|pipetting_pcr_tube_rack|cell_culture_plate|safe_lock_tube_rack_collar_seated|safe_lock_tube_rack_floor_seated|pcr_tube_96_position_carrier|safe_lock_storage_box_derived]`
  - Calls engine operation `ExportRackOpenScad`.
  - Writes one parameterized OpenSCAD source file for the full saved rack.
- `racks carrier-labels-svg RACK_ID OUTPUT.svg [--arrangement ARR_ID] [--template storage_pcr_tube_rack|pipetting_pcr_tube_rack|cell_culture_plate|safe_lock_tube_rack_collar_seated|safe_lock_tube_rack_floor_seated|pcr_tube_96_position_carrier|safe_lock_storage_box_derived] [--preset front_strip_and_cards|front_strip_only|module_cards_only]`
  - Calls engine operation `ExportRackCarrierLabelsSvg`.
  - Writes one carrier-matched front-strip plus module-label SVG sheet for the
    full rack or one selected arrangement block on that rack.
  - Built-in carrier-label presets:
    - `front_strip_and_cards`
    - `front_strip_only`
    - `module_cards_only`
- `racks simulation-json RACK_ID OUTPUT.json [--template storage_pcr_tube_rack|pipetting_pcr_tube_rack|cell_culture_plate|safe_lock_tube_rack_collar_seated|safe_lock_tube_rack_floor_seated|pcr_tube_96_position_carrier|safe_lock_storage_box_derived]`
  - Calls engine operation `ExportRackSimulationJson`.
  - Writes one machine-readable rack geometry/placement JSON export for
    downstream simulation adapters.
- `racks set-profile RACK_ID small_tube_4x6|plate_6|plate_96|plate_384`
  - Calls engine operation `SetRackProfile`.
  - Reflows occupied positions onto another built-in rack/plate profile while
    preserving arrangement order.
- `racks apply-template RACK_ID bench_rows|plate_columns|plate_edge_avoidance`
  - Calls engine operation `ApplyRackTemplate`.
  - Applies one shared authoring shortcut on top of the existing rack snapshot.
  - Built-in templates:
    - `bench_rows`
      - row-major fill
      - clears blocked coordinates
    - `plate_columns`
      - column-major fill
      - clears blocked coordinates
    - `plate_edge_avoidance`
      - column-major fill
      - blocks the outer perimeter
- `racks set-fill-direction RACK_ID row_major|column_major`
  - Calls engine operation `SetRackFillDirection`.
  - Reflows occupied positions onto the same rack geometry under a different
    deterministic fill order.
- `racks set-carrier-pitch RACK_ID MM|--clear`
  - Calls `SetRackCarrierPitch { rack_id, pitch_um }`; `null` clears the override.
  - MM is a positive decimal with at most three fractional digits, converted
    exactly to whole micrometres. Exponents, signs and extra precision refuse.
  - Sets explicit compact Safe-Lock spacing, never an inferred cap clearance.
    Collar-seated minimum is 11.75 mm; floor-seated exports need at least
    12.05 mm; maximum is 20 mm. The setter admits collar-compatible spacing,
    and each export revalidates its own support mode. Old defaults stay 20 mm.
  - Other physical templates refuse a saved override; clear it before those
    exports. Clearing restores the original template output.
  - Grid format is separate: `racks set-custom-profile RACK_ID ROWS COLUMNS`
    retains the pitch. Variable grids/labels do not prove fit to a storage box.
- `racks set-custom-profile RACK_ID ROWS COLUMNS`
  - Calls engine operation `SetRackProfileCustom`.
  - Reflows occupied positions onto one custom A1-style rack geometry.
  - A1-style row labels continue beyond `Z` as `AA`, `AB`, ...
  - Dimensions must be positive and their checked product at most 1,536 slots.
    This is a processing budget, not a physical-capacity or fit claim. Larger
    historical rack snapshots still load, but grid processing and rack exports
    refuse them without truncating the saved data or overwriting an output.
- `racks set-blocked RACK_ID COORD[,COORD...]...`
  - Calls engine operation `SetRackBlockedCoordinates`.
  - Reserves one normalized blocked-coordinate set on the rack and reflows
    occupied positions across the remaining available slots.
- `racks set-blocked RACK_ID --clear`
  - Clears all blocked rack coordinates.

Tube reference models use the existing generic `op` route, also available to
the GUI Shell and inner agent; there is no second labware command parser:

```text
op '{"InspectTubeProfiles":{}}'
op '{"ExportTubeOpenScad":{"profile":"eppendorf_safe_lock_1_5ml","path":"tube.visualization.scad"}}'
op '{"ExportTubeOpenScad":{"profile":"brand_pcr_strip_0_2ml_781320","tube_count":5,"path":"pcr-strip-5.visualization.scad"}}'
```

Inspection returns `tube_profiles` with drawing source URLs/digests and per-field
`drawing`, `approximate_drawing` or `not_documented` basis. Available profile IDs
are `eppendorf_safe_lock_0_5ml`, `eppendorf_safe_lock_1_5ml` and
`eppendorf_safe_lock_2_0ml`, `eppendorf_pcr_tube_0_2ml` and
`brand_pcr_strip_0_2ml_781320`. Export writes the explicitly selected file and
returns the selected profile, but neither operation changes project state.
Normal adapter confirmation policies still apply to artifact writes. The SCAD
is an outer-body/collar visualization with labelled assumptions, not a printable
vessel; cap, hinge and cavity are omitted. No fit with existing rack/plate
templates is inferred. See [drawing provenance](../../data/labware/README.md).

The BRAND profile represents a standard **0.2 ml per tube** strip, article 781320.
`tube_count` selects **1 through 8 contiguous tubes** after cutting; omission
selects eight. Only the retained tubes and their internal bridges are drawn.
Individual-tube profiles default to one and reject counts other than one;
zero or oversized counts fail before overwriting an output. The strip uses the
manufacturer's 9-mm pitch and cited dimensions, but transition/bridge details
are explicitly illustrative. Caps, cavities, numbering and cut stubs/end tab
are omitted. The original Eppendorf individual PCR model still refuses its
undocumented seating collar. No rack binding or inventory change is made.

The existing rack export commands also accept
`safe_lock_tube_rack_collar_seated` and `safe_lock_tube_rack_floor_seated`.
Both are unverified 1.5/2.0 ml carriers with conservative 20-mm spacing,
not PCR racks. All new-template exports disclose support mode, drawing basis,
engineering assumptions and fit/material/temperature nonclaims. The three
older templates retain their geometry and omit the new metadata keys.
`pcr_tube_96_position_carrier` additionally requires an 8x12 grid and uses
nominal 127.76 x 85.48 mm / 9-mm pitch with explicit unverified cap/holder/thermal
nonclaims. It is not a plate or thermocycler block. The PCR reference profile
has no documented sufficient collar; `ExportTubeOpenScad` refuses it cleanly
without changing an existing file rather than fabricating that dimension.

`safe_lock_storage_box_derived` is a separate variable-grid template, also
available through every rack export and generic operation. For example:
`racks openscad RACK_ID derived-box.scad --template safe_lock_storage_box_derived`.
It derives square cells and a separate lid for floor-supported 1.5/2.0 ml
tubes, not a commercial-box profile or qualified fit. Clear a saved compact
pitch with `racks set-carrier-pitch RACK_ID --clear` first; otherwise exports
refuse before writing. Simulation JSON includes optional `storage_box` geometry;
SVGs visibly disclose their round position markers. Caps/hinges may collide;
fit, temperature, material and lid retention remain unverified. See the
[engineering assumptions](../../data/labware/README.md#derived-storage-box).

<!-- manual-page-nav:start -->

---

[CLI Manual](../cli.md) | [Previous: Dotplots, Flexibility and Alignment](dotplot.md) | [Next: Construct Reasoning](construct-reasoning.md)
<!-- manual-page-nav:end -->
