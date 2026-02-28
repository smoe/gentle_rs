# GENtle Roadmap and Status

Last updated: 2026-02-28

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
  - configurable standard letter span threshold
  - optional helical-compressed letter mode with explicit layout selector
    (`continuous-helical` / `condensed-10-row`)
  - dedicated max-span thresholds for continuous-helical and condensed layouts
    (default condensed target: 1500 bp)
  - modulo-10 seam offset control for row alignment (`(bp + offset) % 10`)
  - condensed layout backbone replacement and outward feature-lane clearance
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
2. MCP route now has guarded mutating execution (`op`/`workflow`) and
   UI-intent parity baseline (`ui_intents`, `ui_intent`,
   `ui_prepared_genomes`, `ui_latest_prepared`), but broader parity breadth is
   still incomplete (additional shared-shell route coverage, richer result
   contracts for future mutating UI intents, and wider adapter-equivalence
   coverage).
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
   - condensed 10-row helix-mimic DNA-letter layout baseline is implemented,
     including dedicated span threshold controls, seam-offset behavior,
     backbone replacement, and deterministic annotation-clearance tests;
     dense snapshot/readability benchmarking remains pending
   - any screenshot-based readability baseline artifacts require manual human
     contribution while agent screenshot execution remains policy-disabled
   - zoom/pan policy is not yet unified across canvases
   - UI-level snapshot tests for feature-tree grouping/collapse are pending
   - backdrop-image readability guardrails and stricter grayscale handling are
     incomplete
   - optional all-open-window refresh behavior for applied display changes
     should be finalized consistently across relevant graphics-setting changes
   - sequence-panel amino-acid row rendering is intentionally deferred until
     engine-owned transcript/CDS-aware translation contracts are implemented
     (explicit codon-table resolution plus frame/phase context); current dormant
     AA-row path is maintained as safe no-op
8. Cross-application clipboard interoperability for sequence + feature transfer
   is not yet implemented (current baseline is deterministic in-app extraction).
9. Screenshot bridge execution remains disabled by security policy.
10. Auto-updated documentation with embedded graphics remains postponed.

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
- successful mutating calls persist state to the resolved `state_path`.
- UI-intent MCP tools now route through shared shell parser/executor paths and
  are enforced as non-mutating.
- deterministic adapter-equivalence tests now assert MCP-vs-shared-shell parity
  for all current UI-intent tools.

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

### Linear DNA condensed 10-row layout track (new)

Goal: add a second linear DNA-letter layout optimized for sequences up to a
user-configured span (default target about 1500 bp), with strict base order,
strong X compression, and readability from stable row separation.

Desired layout contract:

- preserve strict nucleotide order (no reordering/skipping)
- render letters in 10 discrete rows (bottom row = modulo class `0`)
- keep X progression monotonic and tightly compressed per bp
- in condensed mode, DNA letters become the primary backbone representation
  (replace/suppress the current black DNA baseline line)
- shift feature annotations/labels outward from the sequence-text band so the
  condensed DNA letters remain readable
- expose manual seam shift via `offset_bp` so conserved columns can be aligned
- row mapping formula: `row = (bp + offset_bp) mod 10`, with `offset_bp` in
  `0..9` (`0` keeps current anchor convention)

Current baseline:

- helical-compressed linear-letter rendering already exists and is controlled by
  display settings/UI knobs (including user-adjustable max span)
- project display setting + UI controls for
  `linear_sequence_helical_phase_offset_bp` already exist with clamp `0..9`
- explicit layout mode selection (`continuous-helical` / `condensed-10-row`)
  is now available and persisted through project display settings
- condensed mode now applies modulo-offset row mapping
  (`(bp + offset_bp) % 10`) and suppresses the black backbone line while DNA
  letters are visible
- feature lanes are pushed outward in condensed mode with deterministic
  readability clearance from the sequence-letter band
- condensed row spacing and baseline band reservation were increased to prevent
  vertical letter overlap and preserve a dedicated readable letter band
- condensed row stacks now anchor outward from the baseline with a deterministic
  minimum per-row step so dense tracks retain full 10-row readability
- dedicated condensed layout span threshold is now persisted and configurable
  (default target 1500 bp) independently of continuous-helical max span

Implementation steps (phased):

1. Renderer mode + mapping:
   - add explicit layout mode (`continuous-helical` vs `condensed-10-row`)
   - implement `(bp + offset_bp) % 10` condensed row mapping
   - suppress/replace black backbone line with DNA letters in condensed mode
   - status: implemented baseline
2. Annotation reflow:
   - push feature annotation lanes/labels outward with deterministic minimum
     spacing from condensed DNA text rows
   - status: implemented baseline
3. Controls + persistence:
   - expose mode + span threshold controls (default target around 1500 bp) in
     Sequence window and Configuration -> Graphics
   - keep live apply/sync + project persistence behavior adapter-equivalent
   - status: implemented baseline
4. Tests + docs:
   - add deterministic tests for offsets `0/3/9`, seam behavior, backbone
     replacement, and annotation clearance
   - update help text/tooltips/docs for the new mode semantics
   - status: implemented baseline

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
   - first-class primer design/test flow (pair interaction checks, practical
     constraints, saved primer sets)
   - keep deterministic engine contracts so GUI/CLI/JS/Lua/MCP stay aligned
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
- Done (2026-02-27): MCP UI-intent parity tests now compare MCP tool outputs
  against direct shared-shell `ui ...` outputs for discovery, prepared-query,
  latest-selection, and open-intent resolution.
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

- Continue dense-case hardening for the linear DNA condensed 10-row layout:
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
