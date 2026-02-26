# GENtle Roadmap and Status

Last updated: 2026-02-26

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
- Feature expert-view pipeline for selected TFBS, restriction sites, and
  splicing groups:
  - shared expert payload generation in engine
  - GUI-side expert panel rendering from engine payload
  - SVG export via `RenderFeatureExpertSvg`
  - shell/CLI/JS/Lua access (`inspect-feature-expert`,
    `render-feature-expert-svg`)
- VCF display filtering parity between GUI and SVG export (`SetParameter`/shared
  display criteria).
- Candidate-set workflow (generate/score/filter/set operations + macro scripts)
  with persistent storage model, including between-anchor window generation.
- Workflow and candidate macro template catalogs now preserve optional external
  reference URLs (`details_url`) alongside template names/descriptions across
  engine and shared shell adapter surfaces.
- Ladder-aware virtual gel rendering and SVG export routes, including
  container-based and arrangement-based serial gel export surfaces.

### GUI baseline in place

- Main lineage page with table/graph/container views.
- Graph contract:
  - single-click selects
  - double-click opens sequence/pool context
  - hover exposes node context (pool range and ladder hints for pool nodes)
  - serial arrangements are rendered as dedicated graph nodes linked from lane-source sequences
- Dense rendering controls including regulatory placement policy and visibility
  toggles persisted through display state.
- Linear DNA-letter rendering controls now include:
  - configurable standard letter span threshold
  - optional helical-compressed letter mode up to higher spans (default 2000 bp)
  - configurable double-strand display with optional 180° reverse-letter rotation
  - optional sequence-panel auto-hide when map letters are visible
- Linear-map drag selection can now be extracted directly to a new sequence via
  `ExtractRegion`, preserving overlapping features in the derived fragment.
- Scrollable/resizable Engine Ops area and shared-shell panel.
- Context-sensitive hover descriptions on actionable controls.
- Help-window shell command reference generated from `docs/glossary.json` with
  interface filter controls (`All`, GUI shell, CLI shell, CLI direct, JS, Lua).
- Native macOS menu mirrors of open windows are available under:
  - `Window -> GENtle Open Windows…`
  - `GENtle -> GENtle Windows…`
  - entries now sync by stable viewport keys (not transient index position)
  - selecting an entry requests focus for the specific window key
  - active window state is mirrored with native menu checkmarks
- Circular sequence-map feature labels now use collision-aware placement:
  labels can slide within feature spans and avoid overlap with already rendered
  labels (for example restriction-site annotations).
- GC-content overlays now use a configurable bin size (`gc_content_bin_size_bp`)
  shared across linear view, circular view, and SVG export.
- Help viewer image handling:
  - markdown images render at constrained width
  - image captions are authored inline in markdown (`*Figure: ...*`)
- Experimental window backdrop styling path:
  - optional per-window-type accent tint (`main`, `sequence`, `splicing`,
    `pool`, `configuration`, `help`)
  - optional image watermark path per window type
  - configuration UI now includes tint color pickers, image-file pickers,
    and live path validation
  - persisted in app settings and live-applied from Configuration -> Graphics

### High-level parity snapshot

| Area | Status |
|---|---|
| Core cloning/editing operations across GUI/CLI/JS/Lua | Done |
| Export/render operations (sequence/lineage/pool gel) | Done |
| Reference genome + track import surfaces | Done |
| Shared shell parity across GUI/CLI | Done |
| Feature expert views (TFBS/restriction/splicing) via shared engine model | Done |
| Candidate strand-relation controls across adapters | Done |
| Alternative-splicing interpretation (lanes/boundaries/events/evidence) | Done (expert-view baseline) |
| Cloning-mode macro presets (SnapGene-style workflows) | Partial (starter templates only) |
| AI communication routes (agent bridge + MCP server) | Partial (agent bridge + guarded MCP op/workflow baseline implemented) |
| Gel simulation realism and arrangement modeling | Partial |
| Shared operation protocol usage | Partial |

Notes:

- Detailed per-operation parity is code-backed; use
  `cargo run --bin gentle_cli -- capabilities` and
  `cargo run --bin gentle_cli -- state-summary` for current runtime inventory.
- `import-pool` and some resource utilities remain adapter-level contracts.

## 2. Active known gaps (priority-ordered)

1. Cloning routine standardization is incomplete:
   - no first-class cloning-routine catalog exists yet to map protocol
     vocabulary to selectable macro entries with explicit family/status metadata
   - workflow macro templates do not yet declare typed input/output contracts
   - lineage graph does not yet render explicit macro-instance box nodes with
     input/output edges
   - protocol-family template packs are still incomplete (restriction-only,
     Gibson, Golden Gate, Gateway, TOPO, TA/GC, In-Fusion, NEBuilder HiFi)
2. MCP route now has guarded mutating execution (`op`/`workflow`), but parity
   breadth is still incomplete (for example UI-intent routing and broader
   adapter-equivalence coverage).
3. Mutating-intent safety policy is not yet fully hardened across agent, voice,
   and MCP invocation paths.
4. Core architecture parity gaps remain:
   - some utilities are still adapter-level rather than engine operations
     (notably `import-pool` and resource-sync utilities)
   - no dedicated engine operation yet for exporting a full run/process as a
     technical-assistant protocol text artifact
   - view-model contract is not yet formalized as a frontend-neutral schema
5. guideRNA workflow remains incomplete (guide-candidate model, oligo
   generation/export, macro template flow; draft in `docs/rna_guides_spec.md`).
6. XML import follow-up remains for `INSDSet/INSDSeq` dialect support.
7. Visualization and workflow UX gaps remain:
   - chromosomal-scale BED overview/density view is missing
   - dedicated primary map-mode splicing view is still pending
   - zoom/pan policy is not yet unified across canvases
   - UI-level snapshot tests for feature-tree grouping/collapse are pending
   - backdrop-image readability guardrails and stricter grayscale handling are
     incomplete
8. Cross-application clipboard interoperability for sequence + feature transfer
   is not yet implemented (current baseline is deterministic in-app extraction).
9. Screenshot bridge execution remains disabled by security policy.
10. Auto-updated documentation with embedded graphics remains postponed.

### MCP server communication track

Goal: add MCP as a first-class AI communication route while keeping one
deterministic engine contract across all adapters.

Current baseline:

- `gentle_mcp` stdio server is implemented.
- baseline tools implemented:
  - `capabilities`
  - `state_summary`
  - `op` (shared engine operation execution; requires explicit `confirm=true`)
  - `workflow` (shared engine workflow execution; requires explicit `confirm=true`)
  - `help`
- successful mutating calls persist state to the resolved `state_path`.
- shared shell + operation contracts already exist for deterministic execution.

Planned work:

1. Add safe UI-intent routing through existing `ui ...` contracts.
2. Add adapter-equivalence tests (CLI shell vs MCP tool invocations) for key
   cloning flows.
3. Keep structured schema compatibility clear across JSON-RPC envelopes and
   MCP tool result payloads.
4. Keep zero MCP-only biology logic branches.

UI-intent tool routine (implementation outline):

1. Capability discovery:
   - add MCP UI-intent discovery surface mirroring `ui intents` output
   - include action/target schema and argument contracts in structured payload
2. Deterministic resolution:
   - add helper routes equivalent to `ui prepared-genomes` /
     `ui latest-prepared` for target disambiguation before open/focus
   - return explicit disambiguation payloads when multiple targets match
3. Guarded execution:
   - route `open`/`focus` through shared `ui ...` parser/executor path
   - require explicit confirmation for mutating UI intents when introduced
4. Result contract:
   - return machine-readable execution report (`executed`, `resolved_target`,
     warnings/errors, follow-up choices)
5. Parity tests:
   - add deterministic equivalence tests asserting MCP UI-intent calls and CLI
     shell `ui ...` commands produce matching resolved targets and outcomes

### Alternative-splicing interpretation track (baseline implemented; follow-ups)

Goal: make exon/intron boundary interpretation and alternative-splicing
inspection explicit, readable, and adapter-equivalent.

Status:

1. Implemented baseline:
   - shared engine-owned splicing payload (`SplicingExpertView`) and selector
     (`FeatureExpertTarget::SplicingFeature`) are in place.
   - boundary markers, event summaries, transcript-vs-exon matrix, and
     junction-support arcs are rendered from the shared payload.
   - GUI expert panel and SVG export use the same payload; shell/CLI/JS/Lua
     can inspect/export via existing feature-expert commands.
2. Remaining follow-ups:
   - add dense-fixture snapshot tests for boundary visibility/non-overlap.
   - add explicit geometry invariants to ensure label text can never alter
     exon/intron footprints.
   - add a dedicated primary map-mode splicing view (beyond expert/detail
     panel embedding) for full-sequence workflows.

### XML import integration track (GenBank-first)

Goal: add XML import support without creating a second semantic model.

Status:

1. Implemented baseline:
   - deterministic runtime import detection order now includes XML fallback:
     `GenBank -> EMBL -> FASTA -> XML`.
   - `GBSet/GBSeq` sequence import is mapped to existing `DNAsequence` +
     feature qualifier structures (no format-specific biology logic branch).
   - genome annotation parser dispatch now ingests `.xml` annotation sources
     and normalizes through existing GenBank-like `GenomeGeneRecord` mapping.
   - GUI open-sequence dialog now exposes XML file filters and docs list XML
     support explicitly.
   - unsupported XML dialects (for example `INSDSet/INSDSeq`) return explicit,
     deterministic diagnostics.
2. Remaining follow-ups:
   - implement additive `INSDSet/INSDSeq` adapter without diverging semantic
     normalization.
   - extend cross-format parity fixtures/tests to include additional XML edge
     cases (multi-record, qualifier-only, interval-only locations).
   - keep large exploratory XML samples out of committed default fixtures.

### Cloning routine catalog + macro-box graph track (new)

Goal: map current cloning vocabulary to executable GENtle routines and represent
macro instances as explicit graph boxes with typed inputs/outputs.

Current baseline:

- shared workflow macro persistence and execution are implemented
  (`UpsertWorkflowMacroTemplate`, `macros template-*`, `macros run`).
- starter macro import pack exists (`assets/cloning_patterns.json`).
- routine-family coverage and graph visualization are still incomplete.

Planned work:

1. Add a versioned cloning-routine catalog manifest with family/tag/status and
   template bindings.
2. Extend macro templates with optional typed input/output port contracts.
3. Persist macro-run instance records (resolved bindings + emitted op ids).
4. Render macro instances as box nodes in lineage/workflow graph with explicit
   input/output edges; support multiple instances of the same routine per
   project.
5. Fill protocol-family packs incrementally (restriction, Gibson, Golden Gate,
   Gateway, TOPO, TA/GC, In-Fusion, NEBuilder HiFi).

Detailed plan and support crosswalk:

- `docs/cloning_routine_catalog_plan.md`

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
- Routine vocabulary crosswalk and graph-node rollout plan are maintained in
  `docs/cloning_routine_catalog_plan.md`.

#### B) Agarose gel simulation improvements (high value)

Current state is strong for ladder-aware pool preview and export, but still
simplified relative to wet-lab interpretation.

Planned upgrades:

- Keep universal gel support for any container/tube (now implemented) and
  harden workflows:
  - preserve one-lane "single-product proof" use cases
  - preserve mixed cardinality lanes (single-sequence and pool containers)
  - improve arrangement edit/authoring UX for lane setup reuse
- Extend arrangement nodes in lineage/DALG:
  - arrangement is already modeled and rendered as an explicit node type that
    groups multiple input tubes under one experimental setup
  - current arrangement mode: `serial` (gel lanes)
  - next arrangement mode: `plate` (plate-reader/assay modeling)
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
  passed (`335 passed, 1 ignored` in main suite; additional suites green).

### Stability TODO (queued, items 5-7)

- Add malformed-annotation reporting that summarizes non-fatal GTF/GFF parse
  issues and exposes file/line context in engine/CLI/GUI.
- Add external-binary preflight diagnostics (BLAST and related tools) with
  explicit "found/missing/version/path" reporting before long jobs start.
- Extend cancellation/timebox controls beyond current genome-track import flow
  so long-running genome-prepare/BLAST jobs can be stopped cleanly.

### Missing test coverage (current priority list)

- Done (2026-02-24): shared-shell execution tests now cover
  `resources sync-rebase` and `resources sync-jaspar` using local fixture
  inputs (no network dependency), including a focused assertion that motif
  reload side effects are applied after `resources sync-jaspar`.
- Add JS adapter tests for `import_pool(...)` and resource-sync wrapper paths.
- Add Lua adapter tests for `import_pool(...)` and resource-sync wrapper paths.
- Add CLI integration tests that confirm `import-pool` / `resources` top-level
  dispatch uses the shared shell parser/executor path.

## 3. Recommended execution order

### Phase A: AI communication + safety plane

- Expand MCP server from guarded op/workflow baseline to broader deterministic
  tool coverage over shared shell/engine routes.
- Keep handlers thin and adapter-equivalent (no MCP-only biology branches).
- Harden mutating-intent safety policy uniformly across agent, voice, and MCP
  invocation paths.
- Implement the UI-intent tool routine in this order:
  discovery -> deterministic resolution -> guarded execution -> structured
  result contract -> CLI-shell/MCP parity tests.
- Keep `ui ...` intent routing deterministic and continue discoverability paths
  (menu, command palette, shell/agent/MCP intent surfaces).

### Phase B: cloning routine standardization

- Execute cloning-routine catalog phases 1-4:
  - catalog schema + routine indexing
  - typed macro input/output contracts
  - macro-run instance recording
  - macro box nodes in lineage/workflow graph
- Add protocol macro template packs for the cloning modes listed in Section 2
  (restriction, Gibson, Golden Gate, Gateway, TOPO, TA/GC, In-Fusion,
  NEBuilder HiFi).

### Phase C: engine/protocol parity hardening

- Keep adapter-level helpers thin and aligned with engine operations.
- Promote remaining adapter-level utilities into first-class engine operations.
- Add process-protocol export contract and richer versioned schema/error
  policy.
- Define shared frontend-neutral view-model schema for tracks/features/overlays.
- Complete XML follow-up (`INSDSet/INSDSeq`) without semantic divergence.

### Phase D: visualization and workflow UX

- Continue alternative-splicing follow-ups:
  - dense-fixture regression tests for boundary visibility and label safety
  - coordinate-true geometry invariants
  - dedicated primary map-mode splicing view
- Continue gel work:
  - arrangement authoring/editing UX (create/update/reorder lanes)
  - `plate` arrangement mode as first-class engine + adapter entity
  - one-run/one-setup semantics per serial arrangement
  - realism upgrades (topology, intensity, co-migration, lane tables)
- Add visual benchmark fixtures and readability regression gates for map export.
- Unify zoom/pan policy and close remaining feature-tree UI snapshot gaps.

### Phase E: integration polish and deferred policy items

- Add cross-application clipboard interoperability through versioned contracts.
- Keep screenshot re-enable work as the final item and only after explicit
  endpoint-security exception/approval.
- Keep auto-updated documentation with embedded graphics postponed until above
  safety/contract priorities are complete.

### Phase F: interpretation (later)

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
