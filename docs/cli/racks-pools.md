# Racks, Pools and Ladders

<!-- manual-counterparts:start -->
Related GUI reference: [Physical Arrangements And Project Layout](../gui/physical-arrangements-and-project-layout.md) | [Pool Distribution (Engine Ops)](../gui/pool-distribution-engine-ops.md) | [Serial Gel (Ladder View)](../gui/serial-gel-ladder-view.md) | [Rack View](../gui/rack-view.md).
These topic links do not imply identical feature coverage across interfaces.
<!-- manual-counterparts:end -->

- Shared shell:
  - Commands:
    - `render-pool-gel-svg IDS|'-' OUTPUT.svg [--ladders NAME[,NAME]] [--containers ID[,ID]] [--arrangement ARR_ID] [--lane-label-layout auto|horizontal|wrapped|staggered|angled] [--band-label-layout auto|inline|panel] [--isoform-markers auto|off]`
    - `render-gel-svg IDS|'-' OUTPUT.svg [--ladders NAME[,NAME]] [--containers ID[,ID]] [--arrangement ARR_ID] [--lane-label-layout auto|horizontal|wrapped|staggered|angled] [--band-label-layout auto|inline|panel] [--isoform-markers auto|off]`
    - `arrange-serial CONTAINER_IDS [--id ARR_ID] [--name TEXT] [--ladders NAME[,NAME]]`
    - `arrange-set-ladders ARR_ID [--ladders NAME[,NAME]]`
    - `racks create-from-arrangement ARR_ID [--rack-id ID] [--name TEXT] [--profile small_tube_4x6|plate_6|plate_96|plate_384]`
    - `racks place-arrangement ARR_ID --rack RACK_ID`
    - `racks move RACK_ID --from A1 --to B1 [--block]`
    - `racks move-blocks RACK_ID --arrangement ARR_ID [--arrangement ARR_ID ...] --to B1`
    - `racks show RACK_ID`
    - `racks labels-svg RACK_ID OUTPUT.svg [--arrangement ARR_ID] [--preset compact_cards|print_a4|wide_cards]`
    - `racks fabrication-svg RACK_ID OUTPUT.svg [--template storage_pcr_tube_rack|pipetting_pcr_tube_rack|cell_culture_plate|safe_lock_tube_rack_collar_seated|safe_lock_tube_rack_floor_seated|pcr_tube_96_position_carrier|safe_lock_storage_box_derived]`
    - `racks isometric-svg RACK_ID OUTPUT.svg [--template storage_pcr_tube_rack|pipetting_pcr_tube_rack|cell_culture_plate|safe_lock_tube_rack_collar_seated|safe_lock_tube_rack_floor_seated|pcr_tube_96_position_carrier|safe_lock_storage_box_derived]`
    - `racks hero-svg RACK_ID OUTPUT.svg [--template storage_pcr_tube_rack|pipetting_pcr_tube_rack|cell_culture_plate|safe_lock_tube_rack_collar_seated|safe_lock_tube_rack_floor_seated|pcr_tube_96_position_carrier|safe_lock_storage_box_derived]`
    - `racks openscad RACK_ID OUTPUT.scad [--template storage_pcr_tube_rack|pipetting_pcr_tube_rack|cell_culture_plate|safe_lock_tube_rack_collar_seated|safe_lock_tube_rack_floor_seated|pcr_tube_96_position_carrier|safe_lock_storage_box_derived]`
    - `racks set-profile RACK_ID small_tube_4x6|plate_6|plate_96|plate_384`
    - `racks set-custom-profile RACK_ID ROWS COLUMNS`
    - `ladders list [--molecule dna|rna] [--filter TEXT]`
    - `ladders export OUTPUT.json [--molecule dna|rna] [--filter TEXT]`
    - `export-pool IDS OUTPUT.pool.gentle.json [HUMAN_ID]`
    - `export-run-bundle OUTPUT.run_bundle.json [--run-id RUN_ID]`
    - `export-lab-instructions OUTPUT.{md|odt|docx} [--format markdown|odt|docx] [--run-id RUN_ID] [--title TEXT] [--audience TEXT]`
    - `import-pool INPUT.pool.gentle.json [PREFIX]`

Pool exchange commands:

- `export-pool IDS OUTPUT.pool.gentle.json [HUMAN_ID]`
- `collections run export-pool CONTAINER_ID --path OUTPUT.pool.gentle.json [--pool-id POOL_ID] [--human-id HUMAN_ID]`
  - Exports explicit sequence IDs (`IDS` is comma-separated) with topology and
    overhang fields into a versioned JSON pool artifact.
  - topology now preserves richer gel-oriented circular-form hints when GENtle
    can infer them (`supercoiled`, `relaxed_circular`, `nicked_circular`)
  - Adds `human_id` at pool level and per-member `human_id`.
- `import-pool INPUT.pool.gentle.json [PREFIX]`
  - Imports pool members into current state; generated IDs are prefixed.
  - richer topology strings such as `supercoiled` or `nicked_circular` are
    preserved as gel-topology hints on imported circular molecules.
- `export-run-bundle OUTPUT.run_bundle.json [--run-id RUN_ID]`
  - Calls engine operation `ExportProcessRunBundle`.
  - Writes a deterministic process/audit bundle JSON
    (`gentle.process_run_bundle.v1`) with:
    - extracted operation inputs
    - chronological parameter overrides
    - selected operation log records
    - output summaries (created/changed ids, exported artifact paths)
  - Omitting `--run-id` exports all operation-log rows.
- `export-lab-instructions OUTPUT.{md|odt|docx} [--format markdown|odt|docx] [--run-id RUN_ID] [--title TEXT] [--audience TEXT]`
  - Calls engine operation `ExportLabAssistantInstructions`.
  - Writes a bench handoff (`gentle.lab_assistant_instructions.v2`) for
    non-IT lab assistants from recorded design operations and current
    sequence/container state. The output format is inferred from the extension
    unless `--format` is supplied.
  - `odt` and `docx` reports embed a lineage overview graphic when SVG
    rasterization succeeds; `md` stays a lightweight text fallback.
  - Includes material IDs, designed outputs, container/arrangement/rack/gel
    references where available, design-derived bench steps, checkpoints, safety
    scope, and record-keeping notes.
  - Omitting `--run-id` summarizes all operation-log rows; with an empty or
    unrelated state, GENtle writes a scaffold and warns that no recorded design
    operations were available.
  - GENtle does not invent reagent volumes, incubation temperatures, or
    kit-specific timings; local SOPs and supervisor-approved conditions remain
    authoritative.
  - Offline demo:
    `workflow @docs/examples/workflows/gibson_arrangements_baseline.json`
    now designs the pGEX/insert Gibson example and writes
    `artifacts/gibson_lab_assistant_handoff.md` for ClawBio/Telegram previews.


DNA ladder catalog commands:

- `ladders list [--molecule dna|rna] [--filter TEXT]`
  - Returns engine-inspected ladder catalog (`schema`, `ladder_count`, `ladders`).
  - Default molecule is `dna`; use `--molecule rna` for RNA ladders.
  - Optional `--filter` keeps only ladders whose names contain `TEXT` (case-insensitive).
- `ladders export OUTPUT.json [--molecule dna|rna] [--filter TEXT]`
  - Calls engine operation `ExportDnaLadders` or `ExportRnaLadders`.
  - Writes structured ladder catalog JSON to disk.

<!-- manual-page-nav:start -->

---

[CLI Manual](../cli.md) | [Previous: Measured Gel Images](gel-images.md) | [Next: Panel Workflows](panels.md)
<!-- manual-page-nav:end -->
