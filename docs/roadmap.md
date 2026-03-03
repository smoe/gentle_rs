# GENtle Roadmap and Status

Last updated: 2026-03-03

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
- Workflow macro pattern import now accepts hierarchical file catalogs:
  - directory-recursive import via `macros template-import PATH`
  - per-template schema support (`gentle.cloning_pattern_template.v1`)
  - starter hierarchy at `assets/cloning_patterns_catalog`
  - GUI `Patterns` menu mirrors folder hierarchy for in-app selection/import
  - `template-import PATH` accepts single template file, legacy pack file, or
    folder tree
  - template execution can be repeated per project with different bindings and
    persists resulting sequences/containers via regular operation persistence
- Typed cloning-routine catalog baseline is now available:
  - versioned manifest at `assets/cloning_routines.json`
  - shared-shell/CLI discovery command:
    `routines list [--catalog PATH] [--family NAME] [--status NAME] [--tag TAG] [--query TEXT]`
  - GUI `Patterns` menu now exposes routine browsing grouped by family/status
    and routine-level template import actions
  - workflow template preflight baseline is now wired into
    `macros template-run`:
    - typed routine input-port checks for supported kinds
      (`sequence`, `container`, `candidate_set`, `guide_set`, `string`,
      `number`, `bool`, `path`, `sequence_anchor`)
    - non-mutating inspection mode:
      `macros template-run ... --validate-only`
    - preflight payload is included in template-run output
  - workflow macro templates now support optional typed
    `input_ports`/`output_ports` metadata:
    - persisted in engine template state
    - import-compatible via `gentle.cloning_pattern_template.v1` files
    - shell authoring path via `macros template-put --input-port ... --output-port ...`
  - preflight now validates both input and output ports and reports
    `contract_source` (`template_ports` or `routine_catalog`) in
    `gentle.macro_template_preflight.v1`
  - preflight now also enforces baseline cross-port semantics:
    - output alias/collision diagnostics
    - sequence/container compatibility diagnostics
    - sequence-anchor semantic checks when one sequence context is bound
  - mutating `macros run` and `macros template-run` now append deterministic
    lineage macro-instance records for success and failure/cancel pathways
    (typed bound inputs/outputs + emitted op ids + status/status_message)
  - macro-instance introspection commands are available:
    `macros instance-list`, `macros instance-show`
- Ladder-aware virtual gel rendering and SVG export routes, including
  container-based and arrangement-based serial gel export surfaces.
- Primer-design report baseline:
  - engine operation `DesignPrimerPairs` now persists deterministic report
    payloads (`gentle.primer_design_report.v1`) in project metadata
  - shared-shell/CLI inspection/export commands are available:
    `primers design`, `primers list-reports`, `primers show-report`,
    `primers export-report`
  - backend selection is now available through engine parameters and shell
    command options:
    - `primer_design_backend=auto|internal|primer3`
    - `primer3_executable` path override
    - `primers design ... --backend ... --primer3-exec ...`
  - report metadata now records backend provenance
    (`requested`, `used`, optional fallback reason + Primer3 executable/version)
  - `auto` mode now falls back deterministically to internal scoring when
    Primer3 is unavailable
- Executable tutorial baseline is now integrated with canonical workflow
  examples:
  - tutorial manifest source:
    `docs/tutorial/manifest.json` (`gentle.tutorial_manifest.v1`)
  - tutorial chapters now include explicit narrative learning objectives and
    concept tags, with generated recurrence mapping ("where this concept
    reappears") in chapter pages and tutorial index
  - committed generated tutorial output:
    `docs/tutorial/generated/` (chapters + retained artifacts + report)
  - dedicated contributor onboarding tutorial chapter added:
    `contribute_to_gentle_development` (backed by executable
    `contribute_gentle_development_baseline` workflow example)
  - `gentle_examples_docs` now supports:
    `tutorial-generate` and `tutorial-check`
  - CI now gates tutorial drift and workflow runtime coverage through
    `tutorial-check` and `cargo test workflow_examples`
  - CI now includes a deterministic internal-preview smoke route that exercises:
    - workflow examples (`load_branch_reverse_complement_pgex_fasta`,
      `digest_ligation_extract_region_minimal`)
    - fast primer-pair report generation + persisted report inspection
  - CI now also runs the full Rust test suite (`cargo test -q`) to catch
    cross-module regressions that may not be covered by targeted smoke steps.

### GUI baseline in place

- Main lineage page with table/graph/container views.
- Main lineage page content is vertically scrollable so graph height does not
  hide container/arrangement sections.
- Main lineage node-group baseline:
  - disjoint node-group model (strict one-group-per-node membership)
  - table indentation under group representative rows
  - graph enclosure outlines for grouped nodes
  - collapsed groups project external edges to representative nodes
  - table/graph node context-menu actions for mark/draft/create workflows
- Graph contract:
  - single-click selects
  - double-click opens sequence/pool context
  - hover exposes node context (pool range and ladder hints for pool nodes)
  - operation transitions render as square intermediate glyphs between node
    circles, with compact operation symbols and optional edge-label text
  - projected parallel transitions (common with collapsed groups) are
    coalesced into one visible edge with grouped operation semantics
  - collapsed group representatives display hidden-internals badges (total
    hidden operation count + top hidden operation-family chips)
  - serial arrangements are rendered as dedicated graph nodes linked from lane-source sequences
  - macro instances are rendered as dedicated box nodes with explicit
    input/output edges (sequence/container-resolved where applicable)
- Dense rendering controls including regulatory placement policy and visibility
  toggles persisted through display state.
- Linear DNA-letter rendering controls now include:
  - adaptive viewport-density routing (default-on) with explicit mode selector:
    `auto-adaptive`, `force-standard`, `force-helical`, `force-condensed-10`
  - compressed-letter toggle semantics scoped to auto mode
    (`linear_sequence_helical_letters_enabled`)
  - deterministic tier routing in auto mode (`1.5x / 2x / 10x` capacity tiers)
    with explicit `OFF` when condensed capacity is exceeded
  - modulo-10 seam offset control for row alignment (`(bp + offset) % 10`)
  - condensed layout backbone replacement and outward feature-lane clearance
  - configurable double-strand display with optional 180° reverse-letter rotation
  - configurable helical strand geometry (`parallel` vs `mirrored`) for
    reverse-strand visual semantics
  - continuous-helical rendering now uses oscillating phase geometry with
    endpoint-anchored width projection (full viewport span retained while
    avoiding monotonic "ramp" artifacts in 1.5x..2.0x density range)
  - helical projection smoothing now interpolates inside compressed columns to
    avoid local x-plateaus/crowding at curve extrema and improve symmetry
  - continuous-helical rendering now applies cycle-level arc-length
    compensation for more uniform bp-to-bp spacing along the helical curve
  - configurable reverse-strand letter opacity now de-emphasizes reverse rows
    consistently across linear map and sequence panel views
  - when backbone hiding is active and letters are visible, baseline tick marks
    are suppressed alongside the backbone line for clearer dense views
  - optional sequence-panel auto-hide when map letters are visible
  - one shared routing helper used by renderer diagnostics and Sequence-window
    status/auto-hide decisions (UI-status parity)
- Linear-map drag selection can now be extracted directly to a new sequence via
  `ExtractRegion`, preserving overlapping features in the derived fragment.
- Sequence-window `Export View SVG` now includes profile-aware export routes:
  - default `screen` profile for current-window composition
  - `wide-context` profile for larger canvas and expanded linear bp context
  - `print-a3` profile with A3 landscape physical-size metadata and higher
    sequence text capacity
  - debug builds include routing-tier diagnostics in the SVG header block
- Scrollable/resizable Engine Ops area and shared-shell panel.
- DNA sequence windows now expose direct genome-anchor extension controls
  (`Extend 5'` / `Extend 3'` with bp + optional output ID) next to anchor
  status, using shared `ExtendGenomeAnchor` engine semantics (same as shell/CLI).
- Context-sensitive hover descriptions on actionable controls.
- Unified interaction policy baseline is implemented across primary canvases and panes:
  - default wheel/trackpad scroll pans/scrolls
  - zoom on canvases is `Shift + wheel`
  - hand-pan drag mode is `Option` (Alt) + drag
  - legacy graph aliases (`Cmd/Ctrl + wheel`, `Space + drag`) remain
    transitionally enabled
  - keyboard pane scrolling (`Arrow`, `PageUp/PageDown`, `Home/End`) is enabled
    for major scroll panes (help, lineage tables/lists, jobs/history,
    configuration, and sequence-window feature/detail panes)
  - linear DNA map now supports two-axis pan:
    - horizontal pan for bp viewport
    - vertical pan for rendered feature-lane stack
  - linear DNA toolbar fit actions are split:
    - `Fit Seq` for full-sequence horizontal fit
    - `Fit Features` for vertical recenter of current subsequence
- Help-window shell command reference generated from `docs/glossary.json` with
  interface filter controls (`All`, GUI shell, CLI shell, CLI direct, JS, Lua).
- Help menu now includes a dedicated `Reviewer Quickstart` guide with internal
  preview warnings, known limitations, and deterministic colleague walkthrough
  steps.
- BLAST UX/provenance hardening baseline:
  - BLAST dialog now emits live query heartbeat status (elapsed runtime while
    `blastn` is running) in addition to query-count progress.
  - status/result views show explicit invocation template and resolved command line.
  - BLAST-hit import operations now persist invocation metadata in
    `ImportBlastHitsTrack.blast_provenance` for operation history/lineage context.
  - agent-suggested shell commands can execute shared BLAST routes
    (`genomes/helpers blast`, `genomes/helpers blast-track`); recursion guardrail
    blocks only nested `agents ask`
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
- Configuration apply workflow:
  - staged settings now use one consistent apply model across external tools,
    graphics, and window styling
  - explicit `Unapplied changes` indicator is shown in Configuration
  - bottom `Cancel`/`Apply` actions are persistent and remain visible while tab
    content scrolls

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
   - typed cloning-routine catalog + template-port preflight baseline is now
     integrated into macro validation and lineage visualization, but semantic
     depth is still incomplete
   - richer preflight semantics are now baseline, but protocol-specific
     constraints remain incomplete (for example routine-family-specific
     compatibility constraints beyond generic alias/container/anchor checks)
   - macro-instance lineage recording now covers success and failure/cancel and
     supports shell introspection, but replay-oriented helpers are still pending
     (for example per-instance re-run from recorded bindings)
   - macro-node drill-down now has a persistent detail panel in GUI lineage
     view, but dense-view controls are still pending (edge-density/aggregation)
   - protocol-family template packs are still incomplete (restriction-only,
     Gibson, Golden Gate, Gateway, TOPO, TA/GC, In-Fusion, NEBuilder HiFi)
   - cross-tool benchmarking (Serial Cloner + MacVector + SnapGene synthesis)
     confirms additional repeated gaps are not yet first-class:
     primer design/validation workflows, auto-annotation library scans,
     sequencing-confirmation workflows, and interactive cloning clipboard/model
   - primer design backend parity is still incomplete:
     - Primer3 backend baseline is now available behind `DesignPrimerPairs`
       (with deterministic auto-fallback to internal backend), but deeper
       constraint-mapping parity still needs hardening
     - backend preflight/version diagnostics are now captured per report, but
       dedicated GUI preflight/status UX is still pending
     - no adapter-equivalence test matrix comparing internal vs Primer3-backed
       report normalization/provenance behavior
2. MCP route now has guarded mutating execution (`op`/`workflow`) and
   UI-intent parity baseline (`ui_intents`, `ui_intent`,
   `ui_prepared_genomes`, `ui_latest_prepared`), but broader parity breadth is
   still incomplete (additional shared-shell route coverage, richer result
   contracts for future mutating UI intents, and wider adapter-equivalence
   coverage).
3. Mutating-intent safety policy is not yet fully hardened across agent, voice,
   and MCP invocation paths.
4. Async long-running command orchestration is still incomplete:
   - BLAST async job-handle/progress/cancel baseline is now available through
     shared shell (`genomes/helpers blast-start|status|cancel|list`) and MCP
     (`blast_async_start|status|cancel|list`)
   - agent auto-execution still needs higher-level orchestration for polling and
     multi-step async flows (it currently executes one suggested command at a
     time)
   - upcoming primer-pair selection workflows are expected to fan out into
     multiple BLAST searches; this still needs workflow-level async orchestration
     on top of the baseline job primitives
5. Core architecture parity gaps remain:
   - some utilities are still adapter-level rather than engine operations
     (notably `import-pool` and resource-sync utilities)
   - no dedicated engine operation yet for exporting a full run/process as a
     technical-assistant protocol text artifact
   - view-model contract is not yet formalized as a frontend-neutral schema
6. guideRNA workflow remains incomplete (guide-candidate model, oligo
   generation/export, macro template flow; draft in `docs/rna_guides_spec.md`).
7. XML import follow-up remains for `INSDSet/INSDSeq` dialect support.
8. Visualization and workflow UX gaps remain:
   - chromosomal-scale BED overview/density view is missing
   - dedicated primary map-mode splicing view is still pending
   - adaptive linear DNA-letter routing baseline is implemented (shared
     renderer/UI decision path, seam-offset behavior, condensed backbone
     replacement, deterministic annotation-clearance tests); dense
     snapshot/readability benchmarking remains pending
   - any screenshot-based readability baseline artifacts require manual human
     contribution while agent screenshot execution remains policy-disabled
   - unified zoom/pan behavior is now implemented for map/graph/help/list panes;
    focused-region fallback behavior and wider regression coverage are pending
   - UI-level snapshot tests for feature-tree grouping/collapse are pending
   - backdrop-image readability guardrails and stricter grayscale handling are
     incomplete
   - optional all-open-window refresh behavior for applied display changes
     should be finalized consistently across relevant graphics-setting changes
   - sequence-panel amino-acid row rendering is intentionally deferred until
     engine-owned transcript/CDS-aware translation contracts are implemented
     (explicit codon-table resolution plus frame/phase context); current dormant
     AA-row path is maintained as safe no-op
   - publication/release visual-polish mode is still pending:
     - add a dedicated `Publication mode` preset for GUI + SVG export paths
     - include deterministic readability-focused defaults (backdrop strength,
       debug overlay visibility, reverse-strand emphasis, typography/spacing)
     - treat exact figure-style tuning as intentionally iterative while core
       functionality and workflows continue to evolve
9. Cross-application clipboard interoperability for sequence + feature transfer
   is not yet implemented (current baseline is deterministic in-app extraction).
10. Screenshot bridge execution remains disabled by security policy.
11. Auto-updated documentation with embedded graphics remains postponed.

### MCP server communication track (UI-intent parity baseline implemented)

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
  - `ui_intents`
  - `ui_intent`
  - `ui_prepared_genomes`
  - `ui_latest_prepared`
  - `blast_async_start`
  - `blast_async_status`
  - `blast_async_cancel`
  - `blast_async_list`
- successful mutating calls persist state to the resolved `state_path`.
- UI-intent MCP tools now route through shared shell parser/executor paths and
  are enforced as non-mutating.
- deterministic adapter-equivalence tests now assert MCP-vs-shared-shell parity
  for all current UI-intent tools plus async BLAST status routing.

Planned work:

1. Extend adapter-equivalence tests (CLI shell vs MCP tool invocations) for key
   cloning flows beyond UI-intent helpers.
2. Keep structured schema compatibility clear across JSON-RPC envelopes and
   MCP tool result payloads.
3. Keep zero MCP-only biology logic branches.
4. For future mutating UI-intent tools, require explicit confirmation and
   preserve adapter-equivalent execution reports (`executed`,
   `resolved_target`, warnings/errors, follow-up choices).

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

### Adaptive linear DNA letter routing track

Goal: route linear DNA-letter rendering by viewport capacity (not fixed bp
thresholds), while preserving strict base order, deterministic compressed
layouts, and renderer/UI parity.

Current baseline:

- shared routing helper (`linear_base_routing`) is now authoritative for both
  renderer mode selection and Sequence-window status/auto-hide decisions
- deterministic auto tiers use density capacity limits:
  - `<= 1.5x`: standard
  - `<= 2x`: helical (when compressed mode enabled)
  - `<= 10x`: condensed-10 (when compressed mode enabled)
  - `> 10x`: `OFF`
- explicit mode override is available:
  - `AutoAdaptive` (default)
  - `StandardLinear`
  - `ContinuousHelical`
  - `Condensed10Row`
- compressed-letter toggle now applies to auto mode only:
  when disabled, auto mode can use standard or `OFF` only
- legacy fixed-threshold `SetParameter` knobs are compatibility-accepted as
  deterministic deprecated no-op messages (no runtime routing effect)
- phase-offset mapping remains available (`(bp + offset_bp) % 10`, clamp `0..9`)
- condensed layout keeps backbone replacement/suppression and deterministic
  outward annotation-lane clearance
- upper-right renderer diagnostics now report adaptive state and metrics
  (active mode, route policy, density, cols-fit, glyph width, reason)
- persisted configuration migration baseline:
  - schema-versioned config
  - legacy helical toggle (`false`) migrates to `true`
  - legacy layout mode migrates to `AutoAdaptive`
  - migrated settings are rewritten once after successful load

Remaining follow-ups:

1. Add snapshot-style dense readability regression assets/benchmarks
   (manual screenshot contribution path remains policy-constrained).
2. Continue condensed/helical stress testing on very large annotation density.

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
- macro import supports both:
  - legacy pack file (`assets/cloning_patterns.json`)
  - per-template hierarchy (`assets/cloning_patterns_catalog/**/*.json`)
- GUI `Patterns` menu mirrors hierarchy and imports templates through shared
  shell command contracts.
- templates can be executed repeatedly with distinct bindings/output IDs in one
  project; resulting outputs persist through regular operation state.
- typed cloning-routine manifest baseline is implemented at
  `assets/cloning_routines.json` (`gentle.cloning_routines.v1`), including
  routine-family/status/tag metadata and typed input/output port declarations.
- shared-shell/CLI discovery now supports routine list/filter/search via
  `routines list ...`.
- GUI routine discovery baseline is now exposed under `Patterns` with grouped
  family/status browsing.
- workflow template preflight baseline is implemented in shared shell:
  - `macros template-run ... --validate-only` reports typed preflight without
    mutating state
  - normal `macros template-run` now carries preflight output and blocks
    execution on preflight errors
  - preflight now validates both input and output ports and reports
    `contract_source` (`template_ports` or `routine_catalog`)
- workflow macro templates now accept optional typed `input_ports` and
  `output_ports` contracts in template metadata.
- preflight now enforces baseline cross-port semantics (alias/collision,
  container/sequence compatibility, anchor checks with bound sequence context).
- mutating `macros run` and `macros template-run` now persist lineage
  macro-instance records for success/failure/cancel
  (bound inputs/outputs + emitted op ids + status/status_message).
- shared-shell introspection for recorded macro instances is now available via
  `macros instance-list` and `macros instance-show`.
- lineage graph and lineage SVG export now include explicit macro box nodes and
  input/output edge rendering.
- GUI lineage now includes a persistent selected-macro detail pane with
  inputs/outputs and emitted operation drill-down.
- routine-family coverage and deeper semantic validation are still incomplete.

Planned work:

1. Extend semantic preflight from baseline generic checks to
   routine-family-specific rules (protocol-aware compatibility models).
2. Expand macro-node dense-view controls (edge-density filtering/aggregation,
   compact operation summaries for large projects).
3. Fill protocol-family packs incrementally (restriction, Gibson, Golden Gate,
   Gateway, TOPO, TA/GC, In-Fusion, NEBuilder HiFi).
4. Postponed: add replay helpers for recorded macro instances only after
   routine-family preflight models and protocol-family packs are stabilized.

Postponed item detail (deferred): macro-instance replay helpers

1. Deferral reason
Replay semantics are currently too volatile while routine-family constraints and
template packs are still moving. Implementing replay now would lock in
incomplete behavior and create churn in replay contracts.

2. Target command surface (when resumed)
`macros instance-replay MACRO_INSTANCE_ID [--transactional] [--validate-only] [--allow-template-drift] [--output-prefix PREFIX]`
`macros instance-diff MACRO_INSTANCE_ID [--against LAST|INSTANCE_ID]`

3. Required replay metadata additions
Store `template_schema`, `template_name`, and resolved binding payload hash.
Store rendered macro script hash and optional rendered script snapshot.
Store operation-journal span boundaries and run fingerprint for deterministic
comparison.

4. Replay execution contract
Load macro instance by ID and resolve replay source.
If template-based and template changed, fail closed unless
`--allow-template-drift` is set.
Re-run preflight with recorded bindings.
Apply optional output-id prefix remapping to avoid collisions in current state.
Execute using the same shared shell parser/executor path as normal macros.
Return a deterministic replay report containing old/new op-id mapping.

5. Diff contract
Compare original and replayed runs by operation family sequence and
created/changed sequence IDs.
Report `equivalent` when semantic outputs match after optional ID remapping.
Report structured mismatch rows when operation families or outputs diverge.

6. Failure/cancel behavior
Replay failure/cancel must append a new macro-instance lineage row with
`status=failed|cancelled` and a status message.
Replay should never mutate state when `--validate-only` is used.

7. Acceptance tests (resume criteria)
Deterministic replay of a template-backed successful macro with
`equivalent=true`.
Deterministic replay with output-prefix remap avoids ID collisions.
Template drift is rejected by default and accepted only with
`--allow-template-drift`.
Replay failure path records failed macro-instance lineage with status message.

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

### Cross-tool parity synthesis (Serial Cloner + MacVector + SnapGene)

Goal: prioritize missing capabilities that recur across multiple external
cloning tools, instead of chasing one-off parity points.

Reference matrix:

- `docs/cloning_tool_gap_matrix.md`

Repeated multi-tool gaps to prioritize:

1. Primer design and validation workflow:
   - baseline now implemented:
     - first-class `DesignPrimerPairs` operation
     - persisted report contract + shell/CLI inspect/export routes
     - optional Primer3 backend selection (`auto|internal|primer3`) with
       deterministic fallback and backend provenance in reports
   - next:
     - expand Primer3 constraint-mapping parity and fixture-backed equivalence
       coverage versus internal backend normalization
     - add richer Primer3 preflight diagnostics/UI surfacing
       (binary/version/config-path checks + environment guidance)
     - pair interaction checks and richer thermodynamic scoring
     - saved/reusable primer sets with explicit versioning
     - async-capable batch off-target/specificity checks so primer-pair
       selection can run multiple BLAST searches through agent/MCP/CLI routes
       with progress/cancel parity
2. Auto-annotation library scan:
   - detect missing/plausible features from curated vector/feature libraries
   - keep machine-readable confidence/overlap diagnostics for automation
3. Sequencing confirmation workflow:
   - import read evidence (Sanger/NGS-aligned scope), map evidence to expected
     constructs, and emit deterministic pass/fail evidence summaries
4. Interactive cloning workspace model:
   - add explicit fragment/clipboard-style assembly workspace with end
     compatibility feedback, while execution still routes through shared ops
5. CRISPR workflow closure:
   - extend existing guide-design baseline to include practical
     screening/confirmation workflow outputs

Notes:

- If visual comparisons include screenshot/raster baselines, those artifacts
  remain manual contributions while screenshot execution is policy-disabled.

### Primer3 wrapper integration track (new)

Goal: add Primer3 tooling support without fragmenting GENtle's shared engine
contracts or adapter parity guarantees.

Phase 1 (wrapper + normalization baseline): implemented baseline

- Engine-owned Primer3 adapter is now wired behind `DesignPrimerPairs`.
- Shared report schema remains unchanged (`gentle.primer_design_report.v1`) and
  deterministic tie-break ordering remains backend-independent.
- Reports now include backend provenance fields
  (`requested`, `used`, optional fallback reason, executable/version).
- Shell/CLI backend override controls are now available:
  - `primers design ... --backend auto|internal|primer3`
  - `primers design ... --primer3-exec PATH`
  - `set-param primer_design_backend ...`
  - `set-param primer3_executable ...`

Phase 2 (tooling diagnostics + compatibility hardening): in progress

- Shell/CLI report payload now captures executable/version diagnostics when
  Primer3 is used (or fallback metadata in auto mode).
- Remaining:
  - dedicated GUI preflight/status views for Primer3 availability
  - additional config-path diagnostics for complex installations
- Keep failure modes deterministic and machine-readable
  (`Unsupported`/`Io`/`InvalidInput` with stable message contracts).
- Add fixture-backed adapter-equivalence tests that assert matching normalized
  report semantics between internal and Primer3 backends for representative
  inputs.

Phase 3 (async specificity tier + agent/MCP parity): baseline started

- Shared shell/CLI + MCP now expose async BLAST primitives:
  - `genomes/helpers blast-start|status|cancel|list`
  - MCP `blast_async_start|status|cancel|list`
- Deterministic parity test baseline now covers MCP-vs-shared-shell async
  BLAST status routing.
- Remaining:
  - primer-pair-specific multi-BLAST orchestration on top of these primitives
  - richer progress granularity and GUI binding for agent-triggered async jobs
  - broader cross-adapter integration tests for cancellation/progress semantics

### Unified BLAST abstraction + primer UI track (new)

Goal: unify BLAST execution semantics for standalone BLAST, single-primer
design, and primer-pair design through one engine-owned service contract.

Detailed execution plan: `docs/primer_design_specialist_plan.md`.

- Primer specialist window + persisted alternative-filter views are tracked in
  `docs/primer_design_specialist_plan.md`.

Phase 1 (engine abstraction + options layering):

- Add shared BLAST service/request/response/progress abstractions in engine,
  reusable by BLAST and primer workflows.
- Implement deterministic option layering:
  - built-in defaults
  - optional defaults file
  - project-level override operation
  - per-request JSON overrides
- Enforce strict option validation (unknown key/type/range failures).
- Persist both raw override JSON and resolved effective options in provenance.

Status (2026-03-01):

- Implemented baseline in engine:
  - layered option resolution (`built-in -> defaults file -> project override -> quick flags -> request JSON`)
  - strict request/default/project JSON object validation
  - threshold filtering (`max_evalue`, `min_identity_percent`,
    `min_query_coverage_percent`, `min_alignment_length_bp`, `min_bit_score`,
    `unique_best_hit`)
  - provenance/report fields for raw request override + resolved effective options
  - project-level parameter operations:
    `blast_options_override`, `blast_options_defaults_path`
- Shared-shell/CLI BLAST routes now accept:
  - `--options-json JSON_OR_@FILE`
  - `--options-file PATH`
- JS/Lua BLAST wrappers now accept optional `options_json` argument and route
  through the same request-options engine path.
- Remaining in Phase 1:
  - GUI advanced options editor/preset UX (Phase 2 UI work item)
  - parity tests that exercise identical option layering across GUI/CLI/JS/Lua.

Phase 2 (BLAST UI alignment):

- Refactor BLAST specialist window into stable sections:
  input, target, options, execution, results.
- Add advanced options JSON editor + preset selector while keeping quick
  controls for common fields.
- Keep heartbeat/query-count progress and invocation visibility wired through
  the shared BLAST progress/result contracts.

Status (2026-03-01):

- Implemented GUI BLAST sectioning baseline (`Target/Input/Options/Execution/Results`).
- Implemented BLAST options controls in GUI:
  - quick controls (`task`, `max_hits`)
  - preset selector
  - structured threshold controls (toggle + typed fields)
  - advanced JSON editor + JSON file loader
  - effective-options preflight preview via shared engine resolver.
- Result panel now displays both request override JSON and resolved effective
  options payloads.
- Implemented BLAST cancellation from both:
  - BLAST dialog execution controls (`Cancel BLAST`)
  - Background Jobs panel (`BLAST -> Cancel`)
  - worker-side interruption now propagates cancellation through shared engine/genomes BLAST paths.
- Remaining in Phase 2:
  - preset catalog persistence and named reusable GUI presets
  - polish structured threshold UX (units/tooltips/preset interop helpers).

Phase 3 (single-primer design UI):

- Add dedicated sequence-context primer-design specialist window:
  - target region selector
  - primer constraints
  - specificity policy (backed by shared BLAST abstraction)
  - ranked primer table with diagnostics
- Add actions for annotation/export and forwarding selected primers into PCR
  workflows.

Phase 4 (primer-pair UI and specificity tier):

- Expand GUI around `DesignPrimerPairs` with explicit amplicon intent, pair
  constraints, and specificity tiers.
- Route pair-specificity BLAST fan-out through the same async BLAST service
  contract used by standalone BLAST.
- Keep report/provenance schema adapter-equivalent.

Phase 5 (adapter parity + tests):

- Add parity tests for request normalization, option resolution, progress, and
  provenance across GUI/CLI/JS/Lua paths.
- Add regression fixtures for primer/specificity edge cases and deterministic
  tie-break behavior in ranked outputs.

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
- If benchmark packs include screenshot/raster artifacts, capture and curation
  are manual contributions (agent screenshot route remains disabled by policy).
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
    - optionally persist recent job-event history across restarts
    - track retry argument snapshots for reproducibility/debugging

### Current branch blockers (must clear first)

- None currently blocking on this branch. Latest local run: `cargo test -q`
  passed (`528 passed, 1 ignored` in main suite; additional suites green).

### Stability TODO (queued, items 5-7)

- Add malformed-annotation reporting that summarizes non-fatal GTF/GFF parse
  issues and exposes file/line context in engine/CLI/GUI.
- Add external-binary preflight diagnostics (BLAST and related tools) with
  explicit "found/missing/version/path" reporting before long jobs start.
- Extend cancellation/timebox controls beyond current genome-track import flow
  (BLAST now supported; continue hardening prepare + future long-running jobs).

### Missing test coverage (current priority list)

- Done (2026-02-24): shared-shell execution tests now cover
  `resources sync-rebase` and `resources sync-jaspar` using local fixture
  inputs (no network dependency), including a focused assertion that motif
  reload side effects are applied after `resources sync-jaspar`.
- Done (2026-02-27): MCP UI-intent parity tests now compare MCP tool outputs
  against direct shared-shell `ui ...` outputs for discovery, prepared-query,
  latest-selection, and open-intent resolution.
- Done (2026-03-03): JS adapter parity tests now compare `import_pool(...)`,
  `sync_rebase_resource(...)`, and `sync_jaspar_resource(...)` wrapper outcomes
  against direct shared-shell execution.
- Done (2026-03-03): Lua adapter parity tests now compare `import_pool(...)`,
  `sync_rebase(...)`, and `sync_jaspar(...)` wrapper outcomes against direct
  shared-shell execution.
- Done (2026-03-03): CLI forwarded-dispatch parity tests now verify
  `import-pool` / `resources sync-rebase` / `resources sync-jaspar` top-level
  routes produce adapter-equivalent outputs and state deltas to direct shared
  shell execution.

## 3. Recommended execution order

### Phase A: AI communication + safety plane

- Expand MCP server from guarded op/workflow baseline to broader deterministic
  tool coverage over shared shell/engine routes.
- Keep handlers thin and adapter-equivalent (no MCP-only biology branches).
- Harden mutating-intent safety policy uniformly across agent, voice, and MCP
  invocation paths.
- Keep the implemented UI-intent tool routine stable and extend parity to
  additional routed command families.
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
- Start repeated cross-tool cloning UX gaps after routine packs land:
  - primer design/validation workflow contracts
  - Primer3 wrapper integration + backend-equivalence test matrix
  - auto-annotation library scan contracts
  - interactive cloning workspace/clipboard model

### Phase C: engine/protocol parity hardening

- Keep adapter-level helpers thin and aligned with engine operations.
- Promote remaining adapter-level utilities into first-class engine operations.
- Add process-protocol export contract and richer versioned schema/error
  policy.
- Define shared frontend-neutral view-model schema for tracks/features/overlays.
- Complete XML follow-up (`INSDSet/INSDSeq`) without semantic divergence.
- Add sequencing-confirmation evidence contracts (read-aligned construct
  validation summaries) as a deterministic shared-engine path.

### Phase D: visualization and workflow UX

- Continue dense-case hardening for adaptive linear DNA letter routing:
  - visual benchmark fixtures and regression gates for crowded labels/features
  - snapshot-style stress coverage for condensed readability constraints
  - manual screenshot contribution for curated visual baselines where required
    (agent screenshot capture remains policy-disabled)
- Continue alternative-splicing follow-ups:
  - dense-fixture regression tests for boundary visibility and label safety
  - coordinate-true geometry invariants
  - dedicated primary map-mode splicing view
- Continue gel work:
  - arrangement authoring/editing UX (create/update/reorder lanes)
  - `plate` arrangement mode as first-class engine + adapter entity
  - one-run/one-setup semantics per serial arrangement
  - realism upgrades (topology, intensity, co-migration, lane tables)
- Add visual benchmark fixtures and readability regression gates for map
  export; treat screenshot/raster baseline assets as manual contributions.
- Add focused-region fallback + regression/snapshot coverage for unified
  scroll/zoom policy and close remaining feature-tree UI snapshot gaps.

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
