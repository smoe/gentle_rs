# GENtle Roadmap and Status

Last updated: 2026-02-20

Purpose: shared implementation status, known gaps, and prioritized execution
order. Durable architecture constraints and decisions remain in
`docs/architecture.md`. Machine contracts remain in `docs/protocol.md`.

## 1. Current implementation snapshot

### Engine and adapter baseline in place

- Shared engine core in `src/engine.rs` with operation/workflow execution.
- Shared shell layer in `src/engine_shell.rs` reused by GUI Shell and
  `gentle_cli shell`.
- CLI adapter in `src/bin/gentle_cli.rs` with state/capability utilities and
  first-class command trees (`genomes`, `helpers`, `resources`, `tracks`,
  `ladders`, `candidates`, `import-pool`).
- JS and Lua adapters expose shared operation/workflow bridges and convenience
  wrappers over engine contracts.

### Biology/analysis capabilities already implemented

- Container-first model in engine state (`ProjectState.container_state`) with
  container-aware digest/merge/ligation/filter operations.
- Reference-genome preparation and extraction (`PrepareGenome`,
  `ExtractGenomeRegion`, `ExtractGenomeGene`) including catalog-backed helper
  flows and BLAST integration.
- Genome track import operations (BED, BigWig via conversion, VCF, BLAST hits)
  with anchor-aware coordinate remapping.
- Resource ingestion/update path for REBASE and JASPAR snapshots across GUI/CLI
  and scripting adapters.
- TFBS annotation guardrails (default cap, explicit unlimited mode), progress
  reporting, and persistent display-time filtering criteria.
- Candidate-set workflow (generate/score/filter/set operations + macro scripts)
  with persistent storage model.
- Ladder-aware virtual gel rendering and SVG export routes.

### GUI baseline in place

- Main lineage page with table/graph/container views.
- Graph contract:
  - single-click selects
  - double-click opens sequence/pool context
  - hover exposes node context (pool range and ladder hints for pool nodes)
- Dense rendering controls including regulatory placement policy and visibility
  toggles persisted through display state.
- Scrollable/resizable Engine Ops area and shared-shell panel.
- Context-sensitive hover descriptions on actionable controls.

### High-level parity snapshot

| Area | Status |
|---|---|
| Core cloning/editing operations across GUI/CLI/JS/Lua | Done |
| Export/render operations (sequence/lineage/pool gel) | Done |
| Reference genome + track import surfaces | Done |
| Shared shell parity across GUI/CLI | Done |
| Shared operation protocol usage | Partial |

Notes:

- Detailed per-operation parity is code-backed; use
  `cargo run --bin gentle_cli -- capabilities` and
  `cargo run --bin gentle_cli -- state-summary` for current runtime inventory.
- `import-pool` and some resource utilities remain adapter-level contracts.

## 2. Active known gaps

- View-model contract is not yet formalized as a frontend-neutral schema.
- Some utilities are still adapter-level rather than engine operations
  (notably `import-pool` and resource-sync utilities).
- No dedicated engine operation yet for exporting a full run/process as a
  technical-assistant protocol text artifact.
- guideRNA workflow is still incomplete (guide-candidate model, oligo
  generation/export, macro template flow; draft in `docs/rna_guides_spec.md`).
- Screenshot bridge execution is intentionally disabled by current security
  policy despite historical implementation work.
- Auto-updated documentation with embedded graphics remains postponed.
- Zoom/pan policy is not yet unified across canvases and should converge to a
  modifier-key-centric contract.

### Missing test coverage (current priority list)

- Add shell execution tests for `resources sync-rebase` and
  `resources sync-jaspar` using local fixture input files (no network dependency).
- Add a focused test that verifies motif reload side effects after
  `resources sync-jaspar`.
- Add JS adapter tests for `import_pool(...)` and resource-sync wrapper paths.
- Add Lua adapter tests for `import_pool(...)` and resource-sync wrapper paths.
- Add CLI integration tests that confirm `import-pool` / `resources` top-level
  dispatch uses the shared shell parser/executor path.

## 3. Recommended execution order

### Phase A: parity hardening

- Keep adapter-level helpers thin and aligned with engine operations.
- Continue parity checks as new operations are introduced.

### Phase B: GUI routing discipline

- Keep UI logic thin; route business logic through engine operations only.

### Phase C: shared view-model contract

- Define frontend-neutral schema for tracks, features, labels, overlays,
  selections, and interaction metadata.

### Phase D: protocol hardening

- Versioned schemas and compatibility policy.
- Richer error taxonomy and validation.
- Operation provenance metadata (engine version/timestamps/input references).
- Process-protocol export contract.
- Keep screenshot re-enable work as the last item in this phase and only after
  explicit endpoint-security exception/approval.

### Phase E: interpretation (later)

- Image/sketch to `StatePatchProposal` translation with confidence scoring and
  explicit confirmation before apply.

## 4. Practical resume checklist

1. Read `docs/architecture.md`.
2. Read `docs/protocol.md`.
3. Read this file (`docs/roadmap.md`).
4. Run quick sanity:
   - `cargo check -q`
   - `cargo run --bin gentle_cli -- capabilities`
   - `cargo run --bin gentle_cli -- state-summary`
5. Continue with highest-priority item from Section 2.
