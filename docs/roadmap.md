# GENtle Roadmap and Status

Last updated: 2026-02-23

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
  `ExtractGenomeRegion`, `ExtractGenomeGene`, `ExtendGenomeAnchor`) including
  catalog-backed helper flows and BLAST integration.
- Genome track import operations (BED, BigWig via conversion, VCF, BLAST hits)
  with anchor-aware coordinate remapping.
- Resource ingestion/update path for REBASE and JASPAR snapshots across GUI/CLI
  and scripting adapters.
- TFBS annotation guardrails (default cap, explicit unlimited mode), progress
  reporting, and persistent display-time filtering criteria.
- Feature expert-view pipeline for selected TFBS and restriction sites:
  - shared expert payload generation in engine
  - GUI-side expert panel rendering from engine payload
  - SVG export via `RenderFeatureExpertSvg`
  - shell/CLI/JS/Lua access (`inspect-feature-expert`,
    `render-feature-expert-svg`)
- VCF display filtering parity between GUI and SVG export (`SetParameter`/shared
  display criteria).
- Candidate-set workflow (generate/score/filter/set operations + macro scripts)
  with persistent storage model, including between-anchor window generation.
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
- Help-window shell command reference generated from `docs/glossary.json` with
  interface filter controls (`All`, GUI shell, CLI shell, CLI direct, JS, Lua).
- Experimental window backdrop styling path:
  - optional per-window-type accent tint (`main`, `sequence`, `pool`,
    `configuration`, `help`)
  - optional image watermark path per window type
  - persisted in app settings and live-applied from Configuration -> Graphics

### High-level parity snapshot

| Area | Status |
|---|---|
| Core cloning/editing operations across GUI/CLI/JS/Lua | Done |
| Export/render operations (sequence/lineage/pool gel) | Done |
| Reference genome + track import surfaces | Done |
| Shared shell parity across GUI/CLI | Done |
| Feature expert views (TFBS/restriction) via shared engine model | Done |
| Candidate strand-relation controls across adapters | Done |
| Cloning-mode macro presets (SnapGene-style workflows) | Planned |
| Gel simulation realism and arrangement modeling | Partial |
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
- Backdrop-image ingest and UX hardening are still incomplete:
  - no dedicated file picker yet for backdrop image paths
  - monochrome conversion currently relies on tinting/asset choice and needs a
    stricter renderer-side grayscale pass
  - per-window readability guardrails (contrast checks, auto-dimming) are not
    yet enforced
- Standardized cloning protocol macros are not yet packaged as first-class
  reusable templates (restriction-only, Gibson, Golden Gate, Gateway, TOPO,
  TA/GC, In-Fusion, NEBuilder HiFi).
- UI-intent command family exists for GUI routing, but mutating-intent guard
  policy is not yet fully hardened for agent/voice-driven invocation paths.
- Chromosomal-scale track overview is still missing: BED-derived features should
  also be visualized at chromosome level, including an optional density view for
  large regions.
- Feature expert-view scope is currently targeted to TFBS and restriction
  sites; future extension to additional feature classes should preserve the same
  `FeatureExpertTarget -> FeatureExpertView -> SVG` contract.
- XML sequence/annotation import is not yet integrated into the shared runtime
  import paths; current primary format remains GenBank (+ FASTA for
  sequence-only).

### XML import integration track (GenBank-first)

Goal: add XML import support without creating a second semantic model.

Execution order:

1. Parser scaffolding and format detection:
   - Add explicit format detection at import boundaries (`LoadFile`,
     annotation parser dispatch) with deterministic precedence:
     GenBank -> FASTA -> XML.
   - Introduce one normalized intermediate import record
     (sequence, topology, feature list, source metadata) used by all formats.
2. Sequence + annotation XML adapter (priority scope):
   - Implement `GBSet/GBSeq` parsing first (smallest useful scope for NCBI
     GenBank XML).
   - Map XML fields to existing `DNAsequence` + feature qualifier structures so
     downstream operations remain unchanged.
   - Reuse the same annotation-to-`GenomeGeneRecord` extraction rules currently
     used for GenBank feature normalization.
3. Annotation-index and genome/helper pipeline integration:
   - Extend annotation parser dispatch in `src/genomes.rs` so prepared
     genome/helper workflows can ingest XML annotation sources.
   - Keep GenBank as preferred catalog source and XML as optional, explicit
     fallback.
4. Adapter/UI exposure and diagnostics:
   - Add XML file filters in GUI open dialogs and document support in CLI/GUI
     manuals.
   - Emit clear unsupported-dialect errors (for example native Bioseq XML when
     only GBSeq adapter is enabled).
5. Parity tests and fixture governance:
   - Add cross-format parity tests based on tiny paired fixtures:
     `test_files/fixtures/import_parity/toy.small.fa`,
     `test_files/fixtures/import_parity/toy.small.gb`,
     `test_files/fixtures/import_parity/toy.small.gbseq.xml`.
   - Assert equal sequence length/content and equivalent mapped gene intervals
     across GenBank and XML imports.
   - Keep large exploratory XML samples out of committed default fixtures.

### Text/voice control track (new)

Goal: make GUI behavior addressable through deterministic text commands first,
then layer voice control on top of the same command plane.

Current baseline:

- `Prepared References...` exists and includes the chromosome line inspector.
- Shell/agent UI-intent routing is implemented:
  - `ui intents`
  - `ui open|focus TARGET ...`
  - `ui prepared-genomes ...`
  - `ui latest-prepared SPECIES ...`
- Prepared references supports one-shot disambiguation/open flow:
  - `ui open prepared-references --species human --latest`
  - explicit `--genome-id` still overrides query-based selection

Planned work:

1. Add guarded execution policy for mutating actions invoked via text/voice.
2. Extend intent grammar/aliasing so natural language maps more reliably to
   deterministic UI-intent commands without guessing.
3. Add optional voice adapter path (STT/TTS) that emits/consumes the same
   deterministic shell/UI intent contracts through the existing agent interface
   execution path.

Scope estimate:

- Moderate (multi-phase): core shell/parser infrastructure already exists;
  missing piece is UI-intent routing + capability discovery + voice adapter.

### SnapGene-parity benchmark track (new)

Goal: align GENtle’s operation-level workflows and visual outputs with the
major protocol families and usability expectations commonly seen in SnapGene,
while keeping GENtle’s shared-engine and open-protocol architecture.

#### A) Macro templates for major cloning modes (priority)

- Restriction cloning:
  - single-fragment insertion, multi-fragment insertion, linear ligation.
- Gibson Assembly:
  - one-insert, multi-insert, circularize-fragment workflows.
- NEBuilder HiFi:
  - one-insert and multi-insert assembly workflows.
- In-Fusion:
  - one-insert and multi-insert assembly workflows.
- Golden Gate:
  - Type IIS overhang-aware assembly templates.
- Gateway:
  - BP, LR, and BP+LR multi-insert workflows.
- TOPO:
  - TA TOPO, blunt TOPO, directional TOPO workflows.
- TA/GC cloning:
  - vector + PCR-fragment templates with orientation handling.

Implementation note:

- Each protocol ships as a macro template (`macros template-put`) backed only by
  shared engine/shell operations, then exposed identically in GUI/CLI/JS/Lua.
- Templates must emit auditable operation logs suitable for DALG-derived
  protocol export.

#### B) Agarose gel simulation improvements (high value)

Current state is strong for ladder-aware pool preview and export, but still
simplified relative to wet-lab interpretation.

Planned upgrades:

- Generalize virtual gels from pool-only usage to any container/tube:
  - support one-lane "single-product proof" use cases
  - support mixed cardinality lanes (single-sequence and pool containers)
- Introduce arrangement nodes in lineage/DALG:
  - arrangement is an explicit node type that groups multiple input tubes under
    one experimental setup
  - first arrangement mode: `serial` (gel lanes)
  - later arrangement mode: `plate` (plate-reader/assay modeling)
- Enforce same-setup semantics for one serial arrangement:
  - all lanes in one virtual gel share one run configuration and one ladder
    context
  - no implicit "spliced gel from different runs" output in a single
    arrangement
- Multi-lane digest conditions with per-lane enzyme sets and batch
  "apply-to-all" behavior.
- Optional uncut/topology-aware migration model for circular DNA
  (supercoiled/nicked/linearized lane behavior).
- Band intensity based on estimated DNA mass per band, not only multiplicity.
- Co-migration grouping thresholds and explicit merged-band annotation.
- Lane-side fragment table with bp, estimated mass, source fragments, and
  cut-site context.
- Optional gel-conditions parameters (agarose %, buffer/model preset) with a
  deterministic default profile.

#### C) Sequence/map clarity benchmarking (high value)

Use a small fixed corpus of examples to compare clarity of GENtle maps and SVG
exports versus established readability patterns.

Planned upgrades:

- Build visual benchmark fixtures (sparse, dense annotations, dense RE sites,
  promoter/regulatory-heavy loci).
- Add explicit label-placement modes where useful (for example,
  "inside-preferred" vs "outside-preferred" behavior by feature class).
- Keep linear/circular/SVG parity for:
  - label anchoring
  - overlap resolution
  - lane separation
  - feature-to-label traceability
- Add quantitative readability checks in snapshot tests:
  - unresolved label overlaps
  - clipped labels
  - hidden high-priority labels (genes/CDS/ORFs/selected RE sites)
- Continue node/table/sequence-window identity consistency
  (`node_id <-> seq_id <-> open-window title`) for lineage traceability.

### UX declutter and readability improvements (planned)

- Add one-click map view presets (`Anchored`, `Cloning`, `Annotation`,
  `Signal`) that apply curated layer visibility bundles.
- Add zoom-aware rendering policy that suppresses tiny low-value glyphs at wide
  zoom and progressively reveals detail when zoomed in.
- Add a one-click declutter action that temporarily disables low-value overlays
  when feature overlap/noise is high.
- Add per-layer counts in visibility controls (for example `ORF (N)`) so users
  can predict visual noise before enabling a layer.
- Separate visual lanes more strictly (annotation rectangles vs predictive
  overlays vs signal tracks) to avoid overlap collisions.
- Expand window visual identity from experimental to production-ready:
  - add backdrop asset validation/file picker
  - add preview panel for each window type
  - enforce readability floor (contrast + max opacity caps)

### Post-delivery hardening backlog (shared)

- Finish-button discoverability hardening:
  - Run a completeness audit so every actionable button has both a concise
    tooltip and a stable hover-status name in the status bar.
  - Add regression checks for newly introduced dialogs/panels so these labels
    do not silently disappear.
- Dense-feature readability hardening:
  - Continue smarter lane packing with stronger non-overlap constraints and
    reduced wasted vertical gaps.
  - Keep strict `REG@DNA` mode enforcing DNA-level placement only (never
    intruding into gene lanes) and add focused visual regression fixtures.
- Undo/redo reliability hardening:
  - Expand operation-level history tests around mixed GUI workflows and
    background-job gating.
  - Keep visible history and transitions consistent after imports/retries.
- Guided anchored-import hardening:
  - Keep preflight summary deterministic for GenBank/BED/BigWig/VCF
    (anchor detection, match status, projected tracks, apply-to-all behavior).
  - Add broader fixture coverage for multi-anchor projects and re-anchoring.
- Command-palette hardening:
  - Keep action naming/keywords stable and improve ranking as command count
    grows.
  - Add keyboard-navigation regressions (`Cmd/Ctrl+K`, arrows, enter, escape).
- Background-jobs robustness:
  - Implemented now:
    - typed job-event records with `kind/phase/job_id/timestamp/summary`
    - monotonic per-job IDs for prepare/track-import/BLAST workers
    - stale worker-message rejection by `job_id`
    - idempotent cancellation handlers reused by dialogs and jobs panel
    - explicit retry events from jobs panel actions
  - Next:
    - add BLAST cancellation support (worker-side interruption)
    - optionally persist recent job-event history across restarts
    - track retry argument snapshots for reproducibility/debugging

### Current branch blockers (must clear first)

- None currently blocking on this branch. Latest local run: `cargo test -q`
  passed (`263 passed, 1 ignored` in main suite; additional suites green).

### Stability TODO (queued, items 5-7)

- Add malformed-annotation reporting that summarizes non-fatal GTF/GFF parse
  issues and exposes file/line context in engine/CLI/GUI.
- Add external-binary preflight diagnostics (BLAST and related tools) with
  explicit "found/missing/version/path" reporting before long jobs start.
- Extend cancellation/timebox controls beyond current genome-track import flow
  so long-running genome-prepare/BLAST jobs can be stopped cleanly.

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
- Start with gel work:
  - add arrangement-node model
  - generalize gel rendering/operations from pools to any tube/container
  - keep one-run/one-setup semantics per serial arrangement
- Add protocol macro template packs for the cloning modes listed in Section 2.
- Add visual benchmark fixtures and readability regression gates for map export.

### Phase B: GUI routing discipline

- Keep UI logic thin; route business logic through engine operations only.
- Add `ui ...` intent routing path so shell/agent commands can open/focus GUI
  dialogs without bypassing existing `src/app.rs` dialog openers.
- Keep prepared-reference discoverability strong:
  - menu route (`File/Genome -> Prepared References...`)
  - command-palette route (`Prepared References`)
  - shell/agent route (`ui open prepared-references`)
  - shell/agent one-shot disambiguation route
    (`ui open prepared-references --species human --latest`).

### Phase C: shared view-model contract

- Define frontend-neutral schema for tracks, features, labels, overlays,
  selections, and interaction metadata.

### Phase D: protocol hardening

- Versioned schemas and compatibility policy.
- Richer error taxonomy and validation.
- Operation provenance metadata (engine version/timestamps/input references).
- Process-protocol export contract.
- UI-intent protocol/versioning and compatibility guarantees for agent/voice
  callers.
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
   - `cargo check -q` (expected green)
   - `cargo run --bin gentle_cli -- capabilities` (run after `cargo check`
     is green)
   - `cargo run --bin gentle_cli -- state-summary` (run after `cargo check`
     is green)
5. Continue with highest-priority item from Section 2.
