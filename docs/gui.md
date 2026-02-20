# GENtle GUI Manual

This page documents the current graphical interface of GENtle.

## Start the GUI

```bash
cargo run --bin gentle
cargo run --bin gentle -- path/to/project.gentle.json
cargo run --features screenshot-capture --bin gentle -- --allow-screenshots
```

The GUI opens an empty project unless a project path is passed on startup.

Screenshot capture policy:

- Default builds keep screenshot capture unavailable.
- Screenshot capture requires:
  - build feature `screenshot-capture`
  - runtime flag `--allow-screenshots`

## Configuration Window

GENtle provides a project-level configuration window (not per-sequence).
It opens in an independent window/viewport, separate from the main project window.

Access:

- Main window menu: `File -> Configuration...` or `Settings -> Configuration...`
- macOS app menu: `GENtle -> GENtle Settings...`
- Shortcut: `Cmd+,`

Tabs:

- `External Applications`
  - Configure external executable overrides for:
    - `rnapkin` (`GENTLE_RNAPKIN_BIN`)
    - `makeblastdb` (`GENTLE_MAKEBLASTDB_BIN`)
    - `blastn` (`GENTLE_BLASTN_BIN`)
    - `bigWigToBedGraph` (`GENTLE_BIGWIG_TO_BEDGRAPH_BIN`)
  - Validate executable availability/version from within the UI.
- `Graphics`
  - Configure project-level display visibility defaults (panels, feature layers, overlays).
  - Applies to sequence windows through the shared engine display state.
  - `Apply + Refresh Open Windows` forces immediate refresh of all currently open sequence windows.

Persistence:

- Configuration is persisted in an app settings file at `~/.gentle_gui_settings.json`.
- Saved settings are restored on app startup.

## Main window layout

A DNA window is split into two visual areas:

- DNA map panel: graphical map (circular or linear)
- Sequence panel: text-oriented sequence rows

Both panels can be shown/hidden from the toolbar.

The project main window (lineage page) supports two views:

- `Table`: tabular lineage view with per-sequence actions
- `Graph`: node/edge lineage visualization
- `Containers`: container list with kind/member-count and open actions

Node click behavior in lineage `Graph` view:

- Single-click: selects a node (highlight only).
- Double-click on a single-sequence node: opens that sequence window.
- Double-click on a pool node: opens a pool-context window (Engine Ops visible,
  pool member distribution available).

## Toolbar buttons

The top toolbar in each DNA window provides these controls (left to right):

1. Circular/Linear toggle
   - Switches DNA map topology visualization between circular and linear.
2. Show/Hide sequence panel
   - Shows or hides the sequence text panel.
3. Show/Hide map panel
   - Shows or hides the graphical DNA map.
4. CDS
   - Shows or hides CDS feature rendering.
5. Gene
   - Shows or hides gene feature rendering.
6. mRNA
   - Shows or hides mRNA feature rendering.
7. Show/Hide TFBS
   - Toggles computed TFBS annotations.
   - Default is off.
8. Show/Hide restriction enzymes
   - Toggles restriction enzyme cut-site markers and labels.
9. Show/Hide GC content
   - Toggles GC-content visualization overlay.
10. Show/Hide ORFs
   - Toggles open reading frame overlays.
11. Show/Hide methylation sites
   - Toggles methylation-site markers.
12. Export Seq
   - Exports the active sequence via engine `SaveFile`.
   - Output format is inferred from filename extension (`.gb/.gbk` => GenBank, `.fa/.fasta` => FASTA).
13. Export SVG
   - Exports the active sequence map via engine `RenderSequenceSvg`.
14. Export RNA SVG (ssRNA only)
   - Exports RNA secondary-structure SVG via shared engine operation `RenderRnaStructureSvg`.
   - Shown only when active sequence is single-stranded RNA (`molecule_type` `RNA`/`ssRNA`).
15. Engine Ops
   - Shows/hides strict operation controls for explicit engine workflows.
16. Shell
   - Shows/hides the in-window GENtle Shell panel.
   - Uses the same shared command parser/executor as `gentle_cli shell`.

Hovering any button shows a tooltip in the UI.

## RNA Structure (ssRNA)

When the active sequence is single-stranded RNA (`molecule_type` `RNA`/`ssRNA`),
the DNA window shows an `RNA Structure (rnapkin)` panel in the top area.

Features:

- `Refresh RNA Structure`
  - Runs `rnapkin` through shared engine APIs to fetch text output and refresh SVG preview.
- `Export RNA SVG`
  - Saves a chosen SVG path through engine operation `RenderRnaStructureSvg`.
- Textual report
  - Shows command metadata and `stdout`/`stderr` from `rnapkin`.
- SVG preview
  - Displays rendered structure image in-panel.

Runtime dependency:

- `rnapkin` must be installed and reachable in `PATH`, or
- set `GENTLE_RNAPKIN_BIN` to the `rnapkin` executable path.

## GENtle Shell (GUI)

The DNA-window toolbar includes a `Shell` button that opens a command panel.

Behavior:

- Uses shared command parsing/execution (`src/engine_shell.rs`).
- Same command set as CLI `gentle_cli shell`.
- Shows a command preview before execution.
- Maintains command output history in the panel.

Supported commands:

- `help`
- `capabilities`
- `state-summary`
- `load-project PATH`
- `save-project PATH`
- `screenshot-window OUTPUT.png` (only in screenshot-enabled builds)
- `render-svg SEQ_ID linear|circular OUTPUT.svg`
- `render-rna-svg SEQ_ID OUTPUT.svg`
- `rna-info SEQ_ID`
- `render-lineage-svg OUTPUT.svg`
- `render-pool-gel-svg IDS OUTPUT.svg [--ladders NAME[,NAME]]`
- `ladders list [--molecule dna|rna] [--filter TEXT]`
- `ladders export OUTPUT.json [--molecule dna|rna] [--filter TEXT]`
- `export-pool IDS OUTPUT.pool.gentle.json [HUMAN_ID]`
- `import-pool INPUT.pool.gentle.json [PREFIX]`
- `resources sync-rebase INPUT.withrefm_or_URL [OUTPUT.rebase.json] [--commercial-only]`
- `resources sync-jaspar INPUT.jaspar_or_URL [OUTPUT.motifs.json]`
- `genomes list [--catalog PATH]`
- `genomes validate-catalog [--catalog PATH]`
- `genomes status GENOME_ID [--catalog PATH] [--cache-dir PATH]`
- `genomes genes GENOME_ID [--catalog PATH] [--cache-dir PATH] [--filter REGEX] [--biotype NAME] [--limit N] [--offset N]`
- `genomes prepare GENOME_ID [--catalog PATH] [--cache-dir PATH] [--timeout-secs N]`
- `genomes blast GENOME_ID QUERY_SEQUENCE [--max-hits N] [--task blastn-short|blastn] [--catalog PATH] [--cache-dir PATH]`
- `genomes extract-region GENOME_ID CHR START END [--output-id ID] [--catalog PATH] [--cache-dir PATH]`
- `genomes extract-gene GENOME_ID QUERY [--occurrence N] [--output-id ID] [--catalog PATH] [--cache-dir PATH]`
- `helpers list [--catalog PATH]`
- `helpers validate-catalog [--catalog PATH]`
- `helpers status HELPER_ID [--catalog PATH] [--cache-dir PATH]`
- `helpers genes HELPER_ID [--catalog PATH] [--cache-dir PATH] [--filter REGEX] [--biotype NAME] [--limit N] [--offset N]`
- `helpers prepare HELPER_ID [--catalog PATH] [--cache-dir PATH] [--timeout-secs N]`
- `helpers blast HELPER_ID QUERY_SEQUENCE [--max-hits N] [--task blastn-short|blastn] [--catalog PATH] [--cache-dir PATH]`
- `helpers extract-region HELPER_ID CHR START END [--output-id ID] [--catalog PATH] [--cache-dir PATH]`
- `helpers extract-gene HELPER_ID QUERY [--occurrence N] [--output-id ID] [--catalog PATH] [--cache-dir PATH]`
- `tracks import-bed SEQ_ID PATH [--name NAME] [--min-score N] [--max-score N] [--clear-existing]`
- `tracks import-bigwig SEQ_ID PATH [--name NAME] [--min-score N] [--max-score N] [--clear-existing]`
- `tracks import-vcf SEQ_ID PATH [--name NAME] [--min-score N] [--max-score N] [--clear-existing]`
- `candidates list`
- `candidates delete SET_NAME`
- `candidates generate SET_NAME SEQ_ID --length N [--step N] [--feature-kind KIND] [--feature-label-regex REGEX] [--max-distance N] [--limit N]`
- `candidates show SET_NAME [--limit N] [--offset N]`
- `candidates metrics SET_NAME`
- `candidates score SET_NAME METRIC_NAME EXPRESSION`
- `candidates score-distance SET_NAME METRIC_NAME [--feature-kind KIND] [--feature-label-regex REGEX]`
- `candidates filter INPUT_SET OUTPUT_SET --metric METRIC_NAME [--min N] [--max N] [--min-quantile Q] [--max-quantile Q]`
- `candidates set-op union|intersect|subtract LEFT_SET RIGHT_SET OUTPUT_SET`
- `candidates macro [--transactional] [--file PATH | SCRIPT_OR_@FILE]`
- `set-param NAME JSON_VALUE`
- `op <operation-json-or-@file>`
- `workflow <workflow-json-or-@file>`

Screenshot command (shared shell, implemented):

- `screenshot-window OUTPUT.png`
  - compile-time gated: available only in builds with feature
    `screenshot-capture`
  - runtime gated: still requires startup opt-in flag `--allow-screenshots`
  - captures active/topmost GENtle window only
  - window lookup uses native AppKit in-process path (no AppleScript bridge)
  - saves to caller-provided output filename/path
- current backend support: macOS (`screencapture`); other platforms currently
  return unsupported

## Candidate-Set Workflow (Engine Ops + GUI Shell)

Candidate-set generation/scoring/filtering is backed by shared engine
operations and available in two GUI paths:

- `Engine Ops -> Candidate sets (scoring/filtering)` dedicated form panel
- GUI shell (`candidates ...` command family, including `candidates macro`)

Recommended flow in one sequence window:

1. Generate a seed set:
   - Engine Ops panel: fill `set`, `seq`, `length`, `step`, optional feature filters, then `Generate`
   - Shell equivalent: `candidates generate sgrnas my_seq --length 20 --step 1 --limit 10000`
2. Add derived metrics:
   - Engine Ops panel: `Score expr` and `Score distance`
   - Shell equivalents:
     - `candidates score sgrnas gc_bias "100 * (gc_fraction - at_fraction)"`
     - `candidates score-distance sgrnas dist_gene --feature-kind gene`
3. Filter into explicit subsets:
   - Engine Ops panel: `Filter` with value and/or quantile bounds
   - Shell equivalents:
     - `candidates filter sgrnas sgrnas_gc_ok --metric gc_bias --min -20 --max 20`
     - `candidates filter sgrnas_gc_ok sgrnas_top --metric dist_gene --max-quantile 0.25`
4. Combine subsets:
   - Engine Ops panel: `Apply set-op` (`union`/`intersect`/`subtract`)
   - Shell equivalent: `candidates set-op intersect sgrnas_top other_set final_set`
5. Inspect/paginate/export:
   - Engine Ops panel: page controls (`limit`/`offset`, `Prev`/`Next`), local sort key, and `Export selected set as JSON`
   - Shell equivalent: `candidates show final_set --limit 50 --offset 0`
6. Optional macro execution:
   - Engine Ops panel: `Run candidates macro` with multiline script
   - Shell equivalent: `candidates macro SCRIPT_OR_@FILE`

Persistence:

- In-memory candidate sets are tracked at
  `metadata["candidate_sets"]` (`gentle.candidate_sets.v1`).
- On save, candidate sets are externalized into a sidecar index + JSONL store;
  project metadata stores a reference schema (`gentle.candidate_sets.ref.v1`).
- On project load, sidecar-backed candidate sets are rehydrated automatically.

## Documentation automation status

- Auto-updated documentation with embedded graphics is explicitly postponed.
- Current policy is manual documentation updates with optional manual screenshot
  artifacts.

## About GENtle

Two About entries exist on macOS:

- `Help -> About GENtle`: custom GENtle About window (icon + version/build text)
- app menu `GENtle -> About GENtle`: standard macOS About panel

The custom window and CLI `--version` share the same text payload.

## Help manuals

The `Help` menu now includes:

- `GUI Manual`: opens `docs/gui.md` in an in-app markdown viewer
- `CLI Manual`: opens `docs/cli.md` in an in-app markdown viewer
- on macOS, app menu `GENtle -> GENtle Help...` opens the same help window
- help now opens in its own native window (separate viewport), not as an overlay in the project window

Help content loading behavior:

- if `docs/gui.md` / `docs/cli.md` exists at runtime, GUI loads those files
- otherwise GUI falls back to embedded copies compiled into the app binary

Markdown image support:

- image rendering is enabled in the help viewer
- use standard markdown image syntax (`![alt](path-or-url)`)
- relative image paths in `docs/*.md` are resolved relative to the markdown file location
- `Reload` in the help window reloads markdown + images from disk
- help viewer supports in-window text search (`Cmd/Ctrl+F` focuses search box)
- search UI includes match count and `Prev`/`Next` navigation

## Map interactions

The DNA map supports mouse interactions:

- Hover: highlights/inspects map elements
- Click: selects a feature
- Double-click: creates a sequence selection from the clicked feature or restriction site

### Linear DNA map: zoom and pan (mouse/touchpad)

When the sequence is in linear mode and your pointer is over the map:

- Zoom with vertical scroll:
  - mouse wheel up/down
  - touchpad two-finger vertical swipe
- Pan with horizontal scroll:
  - touchpad two-finger horizontal swipe
  - on many mice/systems: Shift + mouse wheel
- Zoom is centered around the current cursor position on the map.

Toolbar alternatives (linear mode):

- `-`: zoom out
- `+`: zoom in
- `Fit`: reset view to full sequence
- `Pan` slider: move the current viewport left/right

Notes:

- The map reacts to scroll only while hovered.
- Scroll direction follows your OS/input-device settings (for example natural
  scrolling).

### Lineage graph: zoom and pan (mouse/touchpad)

In `Main window -> Graph` view:

- Zoom:
  - hold `Cmd` (macOS) or `Ctrl` (Windows/Linux), then scroll vertically
    over the graph
  - or use `-`, `+`, `Reset`, `Fit`, and the `Zoom` slider above the graph
  - `Reset Layout` restores default node placement after manual moves
- Pan:
  - use scrollbars, mouse wheel, or touchpad scrolling in the graph viewport
  - hold `Space` and drag on empty graph background
- Node layout:
  - drag a node with the mouse to reposition it
  - double-click a node to open it (`pool` nodes open pool view)
  - single-click selects a node for visual focus
  - hover shows node details at pointer; pool nodes include pool range and
    auto-selected ladder hints
  - moved positions are persisted in project metadata and restored when the
    project is reopened
  - automatic graph layout uses DAG layering (parents left, children right)
    with crossing-minimizing ordering between layers
- Dense graph readability:
  - `Compact labels` reduces label density when the graph is crowded
  - operation labels are shortened and node detail labels are reduced
- Workspace persistence:
  - graph node positions, zoom, compact-label toggle, graph scroll/pan offset,
    and preferred graph/container panel heights are stored in project metadata
  - panel heights can be adjusted with `Graph h` and `Container h`

Why zoom controls differ between views:

- Linear DNA map uses plain vertical scroll for zoom because it has dedicated
  horizontal panning (scroll + pan slider) and a cursor-centered local viewport.
- Lineage graph uses `Cmd/Ctrl + scroll` for zoom so plain scrolling remains
  available for navigating the larger graph canvas.

## Linear map conventions

Current linear map conventions are:

- Forward-strand features are shown above the DNA backbone
- Reverse-strand features are shown below the DNA backbone
- Regulatory features (for example TFBS / `regulatory_region`) are grouped into
  dedicated upper lanes above forward-strand coding features to keep dense loci
  readable
- `REG@TOP` / `REG@DNA` toggle in the map toolbar switches regulatory-feature
  placement between dedicated top lanes and a single near-baseline/GC-strip
  lane
- Directional features use arrow-shaped ends
- Feature labels are lane-packed to reduce overlap
- Coordinate fallback labels are suppressed for unlabeled regulatory features
  to avoid clutter in dense tracks
- Restriction enzyme labels are lane-packed to reduce overlap

## Circular map conventions

Current circular map conventions include:

- Features arranged around the circular backbone
- Restriction enzyme labels around the perimeter
- Optional overlays for GC content, ORFs, and methylation sites

## Pool Distribution (Engine Ops)

When a strict engine operation creates multiple products (for example digest or
ligation), the Engine Ops panel shows a molecular-weight proxy distribution
based on sequence length in base pairs (`bp`).

Bin breaks are computed as follows:

- Compute `span = max_bp - min_bp` across created products.
- Choose `bucket_size` adaptively:
  - `50 bp` if `span <= 500`
  - `100 bp` if `span <= 2000`
  - `250 bp` otherwise
- Anchor buckets at `min_bp` (not around a center value).
- For each product length `bp`, compute:
  - `idx = (bp - min_bp) / bucket_size`
- Display bucket range:
  - `lo = min_bp + idx * bucket_size`
- `hi = lo + bucket_size - 1`

## Pool Gel (Ladder View)

When pool members are available, Engine Ops also shows a virtual agarose-gel
preview:

- one or two DNA ladders are auto-selected to span the pool bp range
- ladders are drawn in lanes next to a pool lane (bands from pool members)
- pool-band labels show bp, with multiplicity when multiple members have the
  same length

Ladder source:

- built-in ladder catalog: `assets/dna_ladders.json` (derived from historical
  GENtle ladder data; upstream legacy files used "marker" naming)
- built-in RNA ladder catalog: `assets/rna_ladders.json`
- historical references:
  - `https://github.com/GENtle-persons/gentle-m/blob/main/src/marker.txt`
  - `http://en.wikibooks.org/wiki/GENtle/DNA_markers`

Controls:

- `ladder preset`: quick selection for common ladder pairs (or Auto)
- `gel ladders`: optional comma-separated ladder names
  - if empty, ladder selection is automatic
- `Export Pool Gel SVG`: writes the current ladder + pool band view to SVG via
  shared engine operation `RenderPoolGelSvg`

## Engine Settings (Engine Ops)

Within `Region extraction and engine settings`, GUI provides:

- `Feature details font` slider (`8.0..24.0 px`)
  - controls the font size used in the feature tree/detail text
  - persists in project display settings (`feature_details_font_size`)
  - `Reset Font` restores default (`11.0 px`)

## Design Constraints Filter (Engine Ops)

The core operations panel includes `Filter Design` for practical design-constraint
screening.

Inputs:

- `Design inputs`: comma-separated sequence IDs
- `GC min` / `GC max`
- `max homopoly`
- `Reject ambiguous`
- `Avoid U6 TTTT`
- `Forbidden motifs` (comma-separated IUPAC motifs)
- `Unique`
- `prefix`

Defaults in the GUI form:

- `GC min=0.30`, `GC max=0.70`
- `max homopoly=4`
- `Reject ambiguous=true`
- `Avoid U6 TTTT=true`

Execution calls engine operation `FilterByDesignConstraints` and creates filtered
in-silico selection outputs.

## Anchored Region Extraction (Engine Ops)

The Engine Ops panel includes an `Extract Anchored` form for promoter-like
design constraints:

- fixed anchor:
  - absolute position (0-based), or
  - feature boundary (`kind`, optional `label`, `Start/End`, `occurrence`)
- direction:
  - `Upstream` or `Downstream`
- flexible target length:
  - `target len` with `tolerance`
- hard constraints:
  - required restriction enzyme sites
  - required TF motifs (ID/name or IUPAC)
  - optional forward/reverse primer constraints
- candidate controls:
  - `Unique`
  - `max candidates`
  - output prefix

Execution calls engine operation `ExtractAnchoredRegion` and creates one or
more candidate sequences as operation results.

## TFBS Annotation (Engine Ops)

The Engine Ops panel includes `TFBS annotation (log-likelihood ratio)`:

- motif selection:
  - `All known JASPAR motifs` mode (explicit all)
  - or selected motifs as comma-separated IDs/names/IUPAC
  - selection helper:
    - filter JASPAR entries
    - add entries to selected motif list with `+`
- global thresholds:
  - `min llr_bits` (absolute log-likelihood ratio score, base 2)
  - `min llr_quantile` (empirical quantile in the scanned sequence region; both strands)
- per-TF overrides:
  - `per-TF min llr_bits` as `TF=VALUE,...`
  - `per-TF min llr_quantile` as `TF=VALUE,...`
- `Clear previous TFBS`:
  - removes prior generated TFBS annotations before writing fresh results

Execution calls engine operation `AnnotateTfbs` and writes scored TFBS features
onto the active sequence. Generated TFBS qualifiers now include four scores:
`llr_bits`, `llr_quantile`, `true_log_odds_bits`, and
`true_log_odds_quantile`.

Safety behavior:

- GUI uses engine default TFBS cap (`500` accepted hits per operation) to keep
  the UI responsive on dense scans.
- CLI/JSON workflows can override this via `AnnotateTfbs.max_hits` (`0` means
  unlimited).

While TFBS annotation is running, GUI shows live progress indicators and keeps
repainting until completion:

- per-motif percentage (based on actual scan steps)
- total percentage across all selected motifs (fraction of TFs addressed)

TFBS display reduction (no recomputation needed):

- `TFBS display filter` offers four checkbox-enabled criteria:
  - `llr_bits`
  - `llr_quantile`
  - `true_log_odds_bits`
  - `true_log_odds_quantile`
- each enabled criterion applies its threshold live to visible TFBS features
- disabling all criteria shows all TFBS features
- sequence SVG export now uses the same TFBS filter settings as the GUI display

VCF display filtering (no recomputation needed):

- `VCF display filter` in Engine Ops offers live filtering criteria for imported
  variant features:
  - class toggles: `SNP`, `INS`, `DEL`, `SV`, `OTHER`
  - `PASS only`
  - `QUAL >= min` and `QUAL <= max` with enable checkboxes
  - required INFO keys (comma-separated)
- changes are written to shared display state and applied immediately.
- sequence SVG export uses the same VCF display filter state as the GUI display.

## Reference Genomes (Main Menu)

Reference-genome workflow is split into separate dialogs (masks) launched from
the main project window menu. `Prepare Reference Genome...` opens as its own
standalone window/viewport (not embedded in the project canvas):

- `File -> Prepare Reference Genome...` (or `Genome -> Prepare Reference Genome...`)
- `File -> Prepared References...` (or `Genome -> Prepared References...`)
- `File -> Retrieve Genome Sequence...` (or `Genome -> Retrieve Genome Sequence...`)
- `File -> BLAST Genome Sequence...` (or `Genome -> BLAST Genome Sequence...`)
- `File -> Prepare Helper Genome...` (or `Genome -> Prepare Helper Genome...`)
- `File -> Retrieve Helper Sequence...` (or `Genome -> Retrieve Helper Sequence...`)
- `File -> BLAST Helper Sequence...` (or `Genome -> BLAST Helper Sequence...`)

Recommended flow:

1. Prepare/cache a genome once:
   - open `Prepare Reference Genome...`
   - set `catalog` + `cache_dir`
   - select `genome` from dropdown values loaded from catalog JSON
   - only genomes that are not yet prepared in the selected cache are shown
   - click `Prepare Genome`
   - this runs in background, shows live progress, and builds local FASTA, gene, and BLAST indexes
   - a running prepare task can be cancelled via `Cancel Prepare`
   - optional `timeout_sec` can be set in GUI (or `--timeout-secs` in CLI/shell)
   - if a previous attempt already downloaded `sequence.fa` but annotation
     resolution fails (for example wrong GTF path), the next retry reuses the
     local FASTA instead of downloading it again
2. Extract a region when needed:
   - open `Retrieve Genome Sequence...`
   - select the same `genome` from dropdown
   - retrieval dropdown lists only genomes already prepared in the selected
     `catalog` + `cache_dir`
   - use `Gene filter` as case-insensitive regular expression
     (for example: `^TP53$`, `^HLA-.*`)
   - use `Biotype filter` checkboxes (auto-populated from parsed annotation data)
   - tune `Top matches` to cap candidate scanning on very large genomes
   - browse candidates page-by-page (`Prev`/`Next`) to avoid rendering huge lists
   - click a filtered gene to auto-fill `chr`, `start_1based`, `end_1based`
   - optionally edit coordinates and set `output_id`
   - click `Extract Selected Gene` (engine op `ExtractGenomeGene`) or
     `Extract Region` (explicit coordinate extraction)
   - coordinates are 1-based and inclusive
3. Run BLAST searches against prepared references:
   - open `BLAST Genome Sequence...`
   - select a prepared `genome` from the selected `catalog` + `cache_dir`
   - choose query source:
     - `Manual sequence` (paste sequence text)
     - `Project sequence` (blast one loaded sequence by id)
     - `Project pool` (blast all members of a selected pool/container)
   - set `max_hits` and `task` (`blastn-short` or `blastn`)
   - click `Run BLAST`
   - BLAST runs in background and keeps the UI responsive; pool mode returns one result set per member
4. Inspect prepared installations when needed:
   - open `Prepared References...`
   - review per-genome install size, readiness flags, source types, and short
     SHA-1 fingerprints
   - use `Retrieve` directly from an inspected row
5. Overlay signal tracks (BED, BigWig, or VCF) onto an extracted sequence:
   - open `Genome -> Import Genome Track...`
   - this opens in its own floating window
   - optionally select a genome-anchored sequence (for one-sequence import)
   - choose a track file (`.bed`, `.bed.gz`, `.bw`, `.bigWig`, `.vcf`, or `.vcf.gz`)
   - optionally set track name and score filters
   - click one of:
     - `Import To Selected`: import onto only the currently selected anchored sequence
     - `Import To All Anchored (One-Time)`: import onto all currently anchored sequences without saving a subscription
     - `Import To All Anchored + Track`: import onto all currently anchored sequences and save a tracked subscription for auto-sync to future anchored extracts
   - import runs in a background task with live progress and can be cancelled
     from the same dialog (`Cancel Import`).
   - tracked files are listed in the same window and can be managed:
     - `Apply now`: re-apply one tracked file to all currently anchored sequences
     - `Remove`: delete one tracked subscription (already imported features remain)
     - `Clear Tracked Files`: delete all tracked subscriptions (already imported features remain)

Equivalent workflow JSON (still supported via workflow runner):

```json
[
  {"PrepareGenome":{"genome_id":"Human GRCh38 Ensembl 116","catalog_path":"assets/genomes.json","cache_dir":"data/genomes"}},
  {"ExtractGenomeGene":{"genome_id":"Human GRCh38 Ensembl 116","gene_query":"TP53","occurrence":1,"output_id":"grch38_tp53","catalog_path":"assets/genomes.json","cache_dir":"data/genomes"}},
  {"ExtractGenomeRegion":{"genome_id":"Human GRCh38 Ensembl 116","chromosome":"1","start_1based":1000000,"end_1based":1001500,"output_id":"grch38_chr1_1000000_1001500","catalog_path":"assets/genomes.json","cache_dir":"data/genomes"}},
  {"ImportGenomeBedTrack":{"seq_id":"grch38_tp53","path":"data/chipseq/peaks.bed.gz","track_name":"H3K27ac","min_score":10.0,"max_score":null,"clear_existing":false}}
]
```

Notes:

- `Catalog...` opens a file picker for genome catalog JSON.
- `cache_dir` can be selected via folder picker; persistent project-local
  storage is recommended over `/tmp`.
- If `catalog` is empty, engine uses default `assets/genomes.json`.
- Bundled `assets/genomes.json` currently includes Human GRCh38 (Ensembl 113 and 116),
  Mouse GRCm39 Ensembl 116, Rat GRCr8 Ensembl 116, Saccharomyces cerevisiae
  S288c (Ensembl 113 and 116), and `LocalProject` (backed by
  `test_files/AB011549.2.fa` + `test_files/AB011549.2.gb`).
- A curated starter catalog for local helper systems is available at
  `assets/helper_genomes.json` (plasmid/lentivirus/adenovirus/AAV plus yeast/E. coli host references).
- Catalog entries may specify `ncbi_assembly_accession` + `ncbi_assembly_name`
  instead of explicit remote URLs. In that case GENtle derives NCBI GenBank/RefSeq
  FTP URLs for sequence (`*_genomic.fna.gz`) and annotation (`*_genomic.gff.gz`)
  during `Prepare Genome`.
- Catalog entries may also specify `genbank_accession` (for helper vectors and
  non-assembly records). If explicit URLs are absent, GENtle derives NCBI EFetch
  sources for FASTA sequence plus GenBank annotation (`gbwithparts`) and then
  indexes extracted feature records for search/retrieval.
- Malformed GTF/GFF lines are now reported as warnings with file/line context
  while valid gene records continue to be indexed.
- Genome track import accepts BED (`.bed` / `.bed.gz`) and BigWig (`.bw` / `.bigWig`) inputs.
- Genome track import also accepts VCF (`.vcf` / `.vcf.gz`) inputs.
- BigWig import uses `bigWigToBedGraph` and can be overridden via `GENTLE_BIGWIG_TO_BEDGRAPH_BIN`.
- VCF score filtering uses VCF `QUAL` for `min_score` / `max_score`.
- Genome track import requires a genome-anchored sequence (`ExtractGenomeRegion` or
  `ExtractGenomeGene`) to map genomic coordinates into local sequence coordinates.
- Imported GenBank files that include an NCBI `ACCESSION ... REGION: start..end`
  header are auto-anchored on load (for example `NC_000001 REGION: 3652516..3736201`),
  so BED tracks can be aligned immediately.
- `ACCESSION ... REGION: complement(start..end)` is interpreted as a reverse
  genomic anchor (`strand -`), and BED intervals are remapped into local
  sequence coordinates accordingly.
- `Import To All Anchored + Track` stores a tracked subscription in project
  metadata and auto-applies it to newly added anchored sequences.
- Imported BED/BigWig/VCF signals are materialized as generated `track` features.
  They render as dense signal tracks (lane-packed in linear view) and appear in
  the feature tree under a dedicated `Tracks` category, grouped by experiment
  (`gentle_track_name`).
- Prepare/Retrieve dialogs show resolved source types for the selected entry
  (`local`, `ncbi_assembly`, `genbank_accession`, `remote_http`).
- Retrieval fields are enabled only after the selected genome is prepared.
- During preparation, a persistent `genes.json` index is built in the genome
  cache to keep retrieval responsive.
- During preparation, GENtle also prepares a BLAST nucleotide index (if
  `makeblastdb` is available).
- BLAST command-line searches use `blastn` and can be configured in
  `Configuration -> External Applications`.
- HTTP downloads support resume/retry behavior and continue from partial files
  when possible.
- Manifest integrity fields (`sequence_sha1`, `annotation_sha1`) are captured
  after successful preparation and shown in `Prepared References...`.
- `start_1based` and `end_1based` fields are constrained to numeric input and
  up to 10 digits.
- Extracted regions are added to project state and opened in sequence windows.

## Engine Ops State Persistence

Engine Ops panel input state is persisted in project metadata per active
sequence id. This includes panel visibility and operation-form inputs.

The same persistence payload now also stores:

- shell panel visibility
- last shell command text

Metadata key format:

- `gui.engine_ops.<seq_id>`

## Loading sequence files

Use the top application menu:

- `File -> Open Sequence...`
- `File -> Open Project...`
- `File -> Open Recent Project...`
- `File -> Prepare Reference Genome...`
- `File -> Prepared References...`
- `File -> Retrieve Genome Sequence...`
- `File -> BLAST Genome Sequence...`
- `File -> Prepare Helper Genome...`
- `File -> BLAST Helper Sequence...`
- `File -> Retrieve Helper Sequence...`
- `File -> Import REBASE Data...`
- `File -> Import JASPAR Data...`
- `File -> Save Project...`
- `File -> Export DALG SVG...`

Shortcut:

- `Cmd+Shift+G` opens `Retrieve Genome Sequence...`.
- `Cmd+Shift+P` opens `Prepare Reference Genome...`.
- `Cmd+Shift+L` opens `BLAST Genome Sequence...`.

Resource import behavior:

- `Import REBASE Data...`
  - Parses a REBASE/Bairoch input file and updates
    `data/resources/rebase.enzymes.json`.
  - Loaded project sequences are refreshed so restriction-site tracks use the
    newly imported enzyme set immediately.
- `Import JASPAR Data...`
  - Parses a JASPAR PFM input file and updates
    `data/resources/jaspar.motifs.json`.
  - TF motif lookups for anchored extraction use the refreshed motif registry.

Supported in the current flow:

- GenBank files (`.gb`, `.gbk` and similar)
- FASTA files (`.fa`, `.fasta`)
  - default interpretation: synthetic blunt `dsDNA`
  - optional FASTA-header metadata tokens for synthetic oligos:
    - `molecule=ssdna` for single-stranded DNA
    - `molecule=rna` for RNA (input `T` is normalized to `U`)
    - `molecule=dsdna` plus optional overhangs:
      - `forward_5=...` (alias `f5=...`)
      - `forward_3=...` (alias `f3=...`)
      - `reverse_5=...` (alias `r5=...`)
      - `reverse_3=...` (alias `r3=...`)

Example FASTA headers:

- `>oligo_ss molecule=ssdna`
- `>oligo_rna molecule=rna`
- `>oligo_ds molecule=dsdna f5=GATC r5=CTAG`

## Notes and limitations

- The GUI is under active development.
- Some advanced feature-location edge cases may still require refinement.
- Visual behavior may continue to evolve as map parity between circular and linear modes improves.
