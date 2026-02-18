# GENtle GUI Manual

This page documents the current graphical interface of GENtle.

## Start the GUI

```bash
cargo run --bin gentle
cargo run --bin gentle -- path/to/project.gentle.json
```

The GUI opens an empty project unless a project path is passed on startup.

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

- Click on a single-sequence node: opens that sequence window.
- Click on a pool node: opens a pool-context window (Engine Ops visible, pool
  member distribution available).
- Double-click on a node: runs `SelectCandidate` for that node.

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

Hovering any button shows a tooltip in the UI.

## About GENtle

Two About entries exist on macOS:

- `Help -> About GENtle`: custom GENtle About window (icon + version/build text)
- app menu `GENtle -> About GENtle`: standard macOS About panel

The custom window and CLI `--version` share the same text payload.

## Help manuals

The `Help` menu now includes:

- `GUI Manual`: opens `docs/gui.md` in an in-app markdown viewer
- `CLI Manual`: opens `docs/cli.md` in an in-app markdown viewer
- help now opens in its own native window (separate viewport), not as an overlay in the project window

Help content loading behavior:

- if `docs/gui.md` / `docs/cli.md` exists at runtime, GUI loads those files
- otherwise GUI falls back to embedded copies compiled into the app binary

Markdown image support:

- image rendering is enabled in the help viewer
- use standard markdown image syntax (`![alt](path-or-url)`)

## Map interactions

The DNA map supports mouse interactions:

- Hover: highlights/inspects map elements
- Click: selects a feature
- Double-click: creates a sequence selection from the clicked feature or restriction site

## Linear map conventions

Current linear map conventions are:

- Forward-strand features are shown above the DNA backbone
- Reverse-strand features are shown below the DNA backbone
- Directional features use arrow-shaped ends
- Feature labels are lane-packed to reduce overlap
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

## Reference Genomes (Main Menu)

Reference-genome workflow is split into two separate dialogs (masks) in the
main project window menu:

- `File -> Prepare Reference Genome...` (or `Genome -> Prepare Reference Genome...`)
- `File -> Retrieve Genome Sequence...` (or `Genome -> Retrieve Genome Sequence...`)

Recommended flow:

1. Prepare/cache a genome once:
   - open `Prepare Reference Genome...`
   - set `catalog` + `cache_dir`
   - select `genome` from dropdown values loaded from catalog JSON
   - click `Prepare Genome`
   - this runs in background, shows live progress, and builds local FASTA index
2. Extract a region when needed:
   - open `Retrieve Genome Sequence...`
   - select the same `genome` from dropdown
   - use `Gene filter` (name/id text filter, analogous to TF interest filtering)
   - click a filtered gene to auto-fill `chr`, `start_1based`, `end_1based`
   - optionally edit coordinates and set `output_id`
   - click `Extract Region`
   - coordinates are 1-based and inclusive

Equivalent workflow JSON (still supported via workflow runner):

```json
[
  {"PrepareGenome":{"genome_id":"Human GRCh38 Ensembl 113","catalog_path":"assets/genomes.json","cache_dir":"data/genomes"}},
  {"ExtractGenomeRegion":{"genome_id":"Human GRCh38 Ensembl 113","chromosome":"1","start_1based":1000000,"end_1based":1001500,"output_id":"grch38_chr1_1000000_1001500","catalog_path":"assets/genomes.json","cache_dir":"data/genomes"}}
]
```

Notes:

- `Catalog...` opens a file picker for genome catalog JSON.
- `cache_dir` can be selected via folder picker; persistent project-local
  storage is recommended over `/tmp`.
- If `catalog` is empty, engine uses default `assets/genomes.json`.
- Retrieval fields are enabled only after the selected genome is prepared.
- Extracted regions are added to project state and opened in sequence windows.

## Engine Ops State Persistence

Engine Ops panel input state is persisted in project metadata per active
sequence id. This includes panel visibility and operation-form inputs.

Metadata key format:

- `gui.engine_ops.<seq_id>`

## Loading sequence files

Use the top application menu:

- `File -> Open Sequence...`
- `File -> Open Project...`
- `File -> Prepare Reference Genome...`
- `File -> Retrieve Genome Sequence...`
- `File -> Import REBASE Data...`
- `File -> Import JASPAR Data...`
- `File -> Save Project...`
- `File -> Export DALG SVG...`

Shortcut:

- `Cmd+Shift+G` opens `Retrieve Genome Sequence...`.
- `Cmd+Shift+P` opens `Prepare Reference Genome...`.

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

## Notes and limitations

- The GUI is under active development.
- Some advanced feature-location edge cases may still require refinement.
- Visual behavior may continue to evolve as map parity between circular and linear modes improves.
