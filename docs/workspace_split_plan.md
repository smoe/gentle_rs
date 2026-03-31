# Workspace / Crate Split Plan

Last updated: 2026-03-31

Purpose: describe the intended multi-crate end state for GENtle and the order
in which the current monocrate should be split without breaking the
"one deterministic engine, thin adapters" rule.

This is an implementation plan, not yet a statement that the split is already
in place. Durable architectural invariants remain in
[`docs/architecture.md`](architecture.md); execution status remains in
[`docs/roadmap.md`](roadmap.md).

Current status:

- the Cargo workspace scaffold and empty member crates now exist;
- production code still mostly lives in the root crate;
- Phase 1 has now started by extracting the first stable id aliases, shared
  enums, and `EngineError`/`ErrorCode` into `crates/gentle-protocol`;
- the extracted contract slice now also includes the dotplot view payloads,
  `SequenceAlignmentReport`, sequence-feature query payloads, and the small
  `TfbsProgress` / `GenomeTrackImportProgress` / `Capabilities` records;
- the extracted slice now also covers the portable RNA-read mapping contract
  layer:
  `SplicingScopePreset`, seed/alignment configs, inspection/export payloads,
  progress previews, and the persisted RNA-read report + summary records;
- the extracted slice now also covers the portable feature-expert vocabulary:
  expert target enums, TFBS/restriction/splicing/isoform view payloads, and
  the shared instruction strings those views expose;
- the root engine surface currently re-exports those extracted types so
  downstream callers do not need to change all at once.

## 1. Goals

The crate split should improve:

- compile-time isolation,
- dependency hygiene,
- headless reuse,
- future browser/WebAssembly optionality,
- and developer navigation.

The split must **not** create:

- adapter-specific biology logic,
- GUI-owned execution semantics,
- or a maze of tiny per-feature crates with unclear ownership.

## 2. Target workspace shape

Suggested end-state workspace members:

### `crates/gentle-protocol`

Shared machine-readable contracts only.

Contents:

- persistent state records,
- operation/workflow payload models,
- report/result/progress payloads,
- stable schema ids and helper enums,
- no execution logic,
- no egui/eframe dependency.

Primary consumers:

- `gentle-engine`
- `gentle-shell`
- `gentle-render`
- `gentle-gui`
- CLI/MCP/JS/Lua/Python adapters

This crate becomes the stable vocabulary layer for cross-adapter parity tests.

### `crates/gentle-engine`

Deterministic biology/workflow execution.

Contents:

- operation application and workflow execution,
- provenance and journal handling,
- candidate, guide, primer, Gibson, RNA-read, sequencing-confirmation, and
  planning logic,
- genome preparation / extraction / BLAST orchestration,
- shared state mutation helpers,
- no egui/eframe dependency.

Likely sources moved here first:

- `src/engine.rs`
- `src/engine/`
- `src/genomes.rs`
- `src/gibson_planning.rs`
- `src/resource_sync.rs`
- parser/import helpers that are execution-facing rather than UI-facing

### `crates/gentle-render`

Headless rendering and export helpers.

Contents:

- lineage SVG export,
- protocol-cartoon rendering,
- sequence SVG export,
- feature-expert SVG export,
- pool-gel rendering,
- other deterministic figure/export paths that should work from CLI/MCP
  without GUI code.

Rule:

- keep this crate headless and reusable;
- do not let egui types leak into it.

### `crates/gentle-shell`

Shared textual command grammar and executor.

Contents:

- shell parsing,
- shell help rendering,
- command execution over `gentle-engine`,
- shell-specific adapter glue that remains headless.

Likely sources:

- `src/engine_shell.rs`
- `src/engine_shell/`
- `src/shell_docs.rs`

This crate should remain the common base for:

- GUI shell,
- `gentle_cli shell`,
- MCP shell-like routes,
- future agent/voice textual command entry points.

### `crates/gentle-gui`

All egui/eframe application code.

Contents:

- `GENtleApp`,
- `MainAreaDna`,
- window management,
- help/settings dialogs,
- background task polling,
- GUI-specific presentation models,
- native-menu integration,
- icons/backdrops and other GUI resources/helpers.

Likely sources:

- `src/app.rs`
- `src/app/`
- `src/main_area_dna.rs`
- `src/main_area_dna/`
- `src/window*.rs`
- `src/dna_display.rs`
- GUI-only render/view helpers that genuinely depend on egui

### Thin adapter/binary layer

Entry points should stay small and explicit:

- `gentle` -> depends on `gentle-gui`
- `gentle_cli` -> depends on `gentle-shell`
- `gentle_mcp` -> depends on `gentle-shell` plus MCP transport glue
- JS/Lua wrappers -> depend on `gentle-engine` and/or `gentle-shell`
- Python wrapper stays out-of-tree as a thin subprocess/client layer

## 3. Intended dependency direction

The important thing is not only *which* crates exist, but which direction the
dependencies point.

```text
gentle-protocol
    ^
    |
gentle-engine ------> gentle-render
    ^                     ^
    |                     |
gentle-shell ------------+
    ^
    |
gentle-gui

thin bins/adapters -> gentle-gui or gentle-shell
```

Interpretation:

- `gentle-protocol` is the bottom shared vocabulary layer.
- `gentle-engine` owns execution and depends on protocol.
- `gentle-render` can depend on protocol and engine-owned read-only helpers, but
  not on GUI.
- `gentle-shell` depends on protocol/engine/render, but not on GUI.
- `gentle-gui` may depend on protocol/engine/render/shell, but nothing below
  should depend on it.

## 4. What should deliberately stay together at first

The first split should **not** produce many micro-crates.

These should stay inside `gentle-engine` initially:

- `genomes`
- `gibson_planning`
- RNA-read analysis
- sequencing confirmation
- planning/meta-layer
- candidate/guide design
- resource sync helpers tightly coupled to engine behavior

Reason:

- they share state/provenance/report semantics heavily,
- the current risk is adapter drift, not insufficient biological sub-crating,
- and too-fine crates would add dependency churn before the core execution
  seams are stable.

Similarly, GUI code should stay in one crate initially instead of being split
immediately by window type.

## 5. Migration order

Recommended order:

### Phase 0: finish monocrate seam discovery

Continue current module decomposition until ownership is obvious and
documentation/file maps make the large files navigable.

Current evidence this is already working:

- `src/engine/` extraction
- `src/engine_shell/` extraction
- `src/app/` extraction
- `src/main_area_dna/` extraction

### Phase 1: extract `gentle-protocol`

Move the most stable shared records first:

- operation/workflow ids and basic envelopes
- report summary rows
- progress records
- stable enums used across GUI/CLI/MCP/JS/Lua

Success criterion:

- adapters can compile against shared contracts without pulling in execution or
  GUI code.

### Phase 2: extract `gentle-engine`

Move deterministic execution next, keeping tests close.

Success criterion:

- direct engine tests still run without egui/eframe,
- GUI and CLI both depend on the extracted engine instead of monolithic
  `src/lib.rs` internals.

### Phase 3: extract `gentle-render`

Peel headless rendering/export paths away from the GUI.

Success criterion:

- CLI/MCP figure exports remain available with no GUI dependency chain.

### Phase 4: extract `gentle-shell`

Move shell parsing/execution after protocol and engine are stable.

Success criterion:

- GUI shell and CLI shell still share the exact same executor crate.

### Phase 5: extract `gentle-gui`

Move egui/eframe application code last.

Success criterion:

- GUI becomes a clean top-layer crate,
- lower crates become easier to test headlessly and eventually reuse in other
  frontends.

## 6. What not to do

Avoid these traps:

- Do not split by adapter first.
  Example: no separate "GUI engine" and "CLI engine" crates.

- Do not split each biology family into its own crate immediately.
  Example: no first-wave `gentle-gibson`, `gentle-genomes`, `gentle-rna-reads`
  crates.

- Do not move rendering into GUI-only crates if the same export is needed from
  CLI/MCP.

- Do not let JS/Lua feature-gating become a reason to fork business logic.

- Do not treat crate extraction as finished unless parity tests still prove the
  same behavior across adapters.

## 7. Near-term concrete moves

The first practical "crate-ready" extraction candidates are:

1. `engine` public contracts now concentrated in `src/engine.rs` plus
   `src/engine/protocol.rs`
2. `engine_shell` public contracts in `src/engine_shell.rs` plus
   `src/engine_shell/`
3. headless export helpers such as:
   - `src/lineage_export.rs`
   - `src/protocol_cartoon.rs`
   - `src/render_export.rs`
   - `src/pool_gel.rs`

The first things that should *not* be split out yet are the GUI monoliths:

- `src/app.rs`
- `src/main_area_dna.rs`

They are now better documented, but they should be extracted into a single
`gentle-gui` crate later rather than fragmented prematurely.

## 8. Definition of done for the split

The split is successful when:

- lower crates compile and test headlessly,
- egui/eframe dependencies sit only in `gentle-gui`,
- GUI/CLI/MCP/JS/Lua still execute the same deterministic engine contracts,
- shell behavior remains shared rather than reimplemented,
- and workspace boundaries make code ownership clearer instead of more opaque.

## 9. First extraction checklist against current `src/lib.rs`

The current public library surface in [`src/lib.rs`](../src/lib.rs) is the best
inventory for planning the first real moves. The mapping below is intentionally
pragmatic rather than perfectly pure.

### Likely future `gentle-protocol` sources

These are not yet clean whole-module moves; they are the contract slices that
should be extracted first:

- `src/engine/protocol.rs`
- public report/progress/state-summary records currently still living in
  `src/engine.rs`
- public genome/blast report records in `src/genomes.rs`
- shared expert-view payloads in `src/feature_expert.rs`

Immediate task:

- identify which public types can move without dragging execution helpers or
  egui types with them.

### Likely future `gentle-engine` modules

Modules that are primarily deterministic biology/state logic:

- `src/engine.rs`
- `src/engine/`
- `src/genomes.rs`
- `src/dna_sequence.rs`
- `src/feature_location.rs`
- `src/gibson_planning.rs`
- `src/resource_sync.rs`
- `src/ncbi_genbank_xml.rs`
- `src/uniprot.rs`
- `src/open_reading_frame.rs`
- `src/methylation_sites.rs`
- `src/tf_motifs.rs`
- `src/restriction_enzyme.rs`
- `src/iupac_code.rs`
- `src/pssm.rs`
- `src/protease.rs`
- `src/amino_acids.rs`
- `src/enzymes.rs`
- `src/dna_ladder.rs`
- `src/feature_expert.rs`

Immediate task:

- keep these together during the first extraction unless a module is clearly
  pure protocol or pure rendering.

### Likely future `gentle-render` modules

These are the current headless export/render candidates:

- `src/lineage_export.rs`
- `src/protocol_cartoon.rs`
- `src/pool_gel.rs`
- `src/render_export.rs`

Maybe later, if they remain headless:

- other export-only figure helpers that do not import egui types

Explicitly **not** in `gentle-render` at first:

- `src/render_dna.rs`
- `src/render_dna_linear.rs`
- `src/render_dna_circular.rs`
- `src/render_sequence.rs`
- `src/sequence_rows.rs`
- `src/sequence_rows_*.rs`

Reason:

- those are currently egui-bound rendering widgets and belong with the GUI
  crate until/unless a true headless abstraction exists.

### Likely future `gentle-shell` modules

- `src/engine_shell.rs`
- `src/engine_shell/`
- `src/shell_docs.rs`

Potential later neighbors, depending on cleanup:

- `src/mcp_server.rs` may stay a thin adapter on top of shell instead of moving
  into the shell crate itself
- `src/agent_bridge.rs` may remain a separate top-layer helper until its
  transport responsibilities settle

### Likely future `gentle-gui` modules

- `src/app.rs`
- `src/app/`
- `src/main_area_dna.rs`
- `src/main_area_dna/`
- `src/window.rs`
- `src/window_dna.rs`
- `src/dna_display.rs`
- `src/render_dna.rs`
- `src/render_dna_linear.rs`
- `src/render_dna_circular.rs`
- `src/render_sequence.rs`
- `src/sequence_rows.rs`
- `src/sequence_rows_blank.rs`
- `src/sequence_rows_dna.rs`
- `src/sequence_rows_restriction_enzymes.rs`
- `src/window_backdrop.rs`
- `src/icons.rs`
- `src/egui_compat.rs`
- `src/scroll_input_policy.rs`

Immediate task:

- keep GUI-specific presentation code together even if some pieces look
  "render-like"; the boundary we want first is headless vs egui-bound, not
  "drawing" vs "not drawing".

### Likely thin adapter/top-layer modules

These probably stay above the first-wave crate split or move only later:

- `src/mcp_server.rs`
- `src/js_interface.rs`
- `src/lua_interface.rs`
- `src/agent_bridge.rs`
- `src/workflow_examples.rs`
- `src/about.rs`
- `src/tool_overrides.rs`
- `src/test_support.rs`

They are important, but they should not distort the first extraction waves.

## 10. First practical extraction milestone

The first milestone I would actually implement is:

1. create a workspace with empty `gentle-protocol` and `gentle-engine` crates
2. move only obviously stable protocol records first
3. make the current root crate temporarily depend on those extracted crates
4. keep GUI/shell code in place until protocol + engine boundaries prove
   themselves in tests

That gives the project a real workspace without forcing a risky all-at-once GUI
move.
