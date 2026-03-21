# GENtle Architecture (Working Draft)

Last updated: 2026-03-11

This document describes how GENtle is intended to work and the durable
architecture constraints behind implementation choices.

It has two goals:

- Human-readable architecture guide
- Stable record of architectural decisions and invariants

Related shared documents:

- `docs/protocol.md`: operation/state/result contracts
- `docs/roadmap.md`: implementation status, known gaps, and execution order

## 1. Product intent

GENtle is a DNA/cloning workbench with multiple access paths:

- GUI for interactive use
- JavaScript shell for scripted experiments
- Lua shell for scripted experiments
- Python module wrapper for scripted/notebook automation
- CLI for automation and AI tools
- MCP server route for tool-based AI integration, including capability
  discovery/negotiation (`tools/list`, `capabilities`)

The long-term requirement is strict behavioral parity:

- The same biological operations must run through the same core routines,
  regardless of entry point.
- This document does not describe any current protocol that could be used
  to change the human genome or that of animals or plants.

Wet-lab semantic rule (target model):

- A main DNA window represents a wet-lab container (tube/vial), not a single
  guaranteed molecule.
- A container may hold multiple candidate molecules/fragments.
- Filter-like steps (PCR, gel extraction, in silico selection) produce a new
  container with a narrower candidate set.

Strategic aims:

1. Keep one deterministic engine contract across GUI, CLI, JS, Lua, Python,
   and MCP.
2. Preserve provenance so every derived result can be traced and replayed.
3. Make every process exportable as a human-readable protocol text that a
   technical assistant can follow step by step (inputs, operations, expected
   outputs, and checkpoints).
4. Support optional, explicit screenshot artifact generation for documentation
   and progress communication without weakening default safety boundaries
   (compile-time feature gate + explicit runtime opt-in).
5. Provide context-sensitive hover descriptions for actionable GUI buttons so
   control intent is explicit and shared button references are unambiguous.
6. Require explicit provenance documentation for all test data fixtures so
   source/reconstruction and usage context are always clear.
7. Support optional per-window-type visual identity cues (subtle monochrome
   image backdrops + color accents) without harming scientific readability.
8. Make sequence/feature transfer pathways explicit and deterministic,
   including in-app selection extraction and planned cross-application
   clipboard interchange through the same engine operation contracts.

Test-data provenance rule:

- Every committed test fixture (for example under `test_files/` or `tests/`)
  must have clear provenance documentation.
- Documentation must include at least:
  - exact upstream/source reference when known, or an explicit statement that
    it is synthetic/hand-crafted
  - deterministic self-creation/retrieval steps
  - where in GENtle the fixture is used (runtime path, parser path, and/or
    test names)
- New fixture files are not considered done until this provenance note exists
  in a nearby README or equivalent fixture manifest.

Import-format policy (GenBank-first, XML additive):

- GenBank remains the canonical annotation import format in current production
  paths.
- XML support is additive and must normalize to the same internal structures
  used by GenBank import before touching business logic.
- The first XML scope is sequence + feature/annotation import only (no
  format-specific behavior divergence).
- XML parser adapters should target explicit NCBI dialects (`GBSet/GBSeq`
  first, `INSDSet/INSDSeq` second) and reject unknown dialects with clear
  diagnostics.
- Cross-format fixtures must remain small and paired (same biological content
  across `.fa`, `.gb`, `.xml`) so parity tests can detect semantic drift.

Information-preserving transfer default (decision):

- Default behavior for import/retrieval/transcode paths is to preserve and
  transfer relevant structured biological information (sequence + annotation)
  whenever source data contains it.
- For genome-region retrieval specifically, annotation projection defaults to
  `annotation_scope=core` when unset.
- `annotation_scope=full` is available for richer transfer and can be paired
  with `max_annotation_features` to keep payload size bounded.
- Opt-out is allowed when associated data volume is large and not immediately
  relevant to the active cloning task, but this must be explicit at the call
  site.
- Gene-level annotation is considered core cloning context and should not be
  omitted by default.

Tooltip coverage rule:

- Every actionable GUI button must include a concise, context-sensitive
  hover description (`on_hover_text`).
- This includes menu actions, toolbar buttons, dialog primary/secondary
  actions, shell controls, and strict engine operation triggers.
- Definition of done for new controls: button label + hover description
  are both present at implementation time.
- Maintenance check: scan button call sites and verify coverage before
  release (for example via `rg` over `.button`/`.small_button` call sites).

Visual consistency rule:

- Legends must map 1:1 to rendered glyph semantics (shape + color + label).
- If node/track glyph semantics change, legend text/colors must be updated in
  the same change.
- Interaction contracts should be explicit where ambiguity is costly
  (for example, click vs double-click behavior).

cDNA presentation rule:

- `cDNA` mode defaults to intrinsic cDNA evidence and must not automatically
  flood the map with contextual transcript projections.
- Contextual transcript-projection features should be tagged explicitly in
  qualifiers (for deterministic filtering across adapters/reloads).
- Showing all contextual transcripts in `cDNA` mode is an explicit opt-in
  action, not default behavior.

Interaction consistency rule:

- GUI scroll/zoom interaction must follow one cross-surface policy:
  - default wheel/trackpad motion scrolls or pans content
  - zoom is reserved for `Shift + wheel` on zoom-capable canvases
  - `Option` (Alt) + drag enters hand-pan mode on zoom-capable canvases
  - in linear DNA view, wheel/trackpad pan must support both axes:
    horizontal movement pans bp viewport; vertical movement pans rendered lane stack
- Cursor affordances must mirror mode intent:
  - `Option` hover uses hand (`Grab`)
  - `Option` drag uses grabbing hand (`Grabbing`)
  - `Shift` hover on zoom-capable canvases uses zoom cursor (`ZoomIn`/`ZoomOut`)
- Arrow keys pan active scroll surfaces; `Shift + arrows` map to zoom intents
  on zoom-capable canvases.
- Linear DNA toolbar `Fit` actions must be explicit and non-ambiguous:
  - `Fit Seq`: full sequence horizontal fit
  - `Fit Features`: geometry-based vertical fit for currently shown subsequence
- Non-zoom scroll panes should expose keyboard navigation (`Arrow`,
  `PageUp/PageDown`, `Home/End`) when the pane is active and text fields do not
  hold keyboard focus.
- Temporary compatibility aliases are allowed (for example legacy
  `Cmd/Ctrl + wheel` or `Space + drag`) but must be documented and treated as
  transitional.

Window visual identity rule:

- Backdrop styling is decorative and optional; it must never reduce readability
  of sequence/map labels, controls, or quantitative overlays.
- Window-type coloring must remain stable (`main`, `sequence`, `splicing`,
  `pool`, `configuration`, `help`) and should not encode scientific semantics.
- Backdrop images are expected to be monochrome or monochrome-tinted in
  rendering, with low opacity and a fast global disable path.

Discoverability rule:

- High-frequency actions must be reachable from both menus and searchable
  command surfaces (Command Palette).
- GUI status reporting should be centralized (status bar + background jobs
  panel) so long-running work is visible and cancellable.
- Hovered actionable controls should expose a stable human-readable name for
  shared debugging/support language.
- Command discoverability and documentation should be glossary-driven:
  - shell command help is generated from `docs/glossary.json` (single source of truth)
  - the help viewer should support interface/language filtering (GUI shell,
    CLI shell, CLI direct, JS, Lua, MCP, or all) to avoid duplicated manuals
    while keeping context-specific views.
- Menu information architecture should separate sequence ingress from reference
  lifecycle:
  - `File` may include sequence ingress/start actions (for example local open,
    accession fetch, and sequence retrieval) because users treat these as
    project-seeding entry points.
  - `Genome` remains the canonical home for prepared-reference lifecycle and
    heavy analysis workflows (prepare, BLAST, track import, inspection).
- Prepared-reference inspection must remain directly reachable from
  `Genome -> Prepared References...`, and searchable as
  `Prepared References` in Command Palette.
- Chromosome inspection currently lives inside `Prepared References...` as an
  embedded `Chromosome inspector` section (one proportional line per contig).
- Configuration dialogs with staged edits must keep primary commit actions
  (`Cancel`/`Apply`) persistently visible (not scroll-hidden) and expose an
  explicit unapplied-changes indicator.

## 2. Core architecture rule

GENtle should follow a single-engine architecture:

1. `Core engine` (single source of truth)
2. `Adapters/frontends` (GUI, JS, Lua, Python, CLI, MCP server)
3. `View model/renderers` (presentation layer only)

### Non-negotiable invariant

Frontends must not implement their own cloning/business logic.
They only translate user input into engine operations and display results.
This includes MCP tool handlers.

## 3. Implementation status and roadmap

Shared implementation status, parity matrices, known gaps, and execution order
are now maintained in `docs/roadmap.md`.

## 4. Engine model

### State

`ProjectState` stores:

- `sequences: HashMap<SeqId, DNAsequence>`
- current sequence windows are DNA/RNA-oriented. First-class protein sequence
  windows are intentionally deferred. UniProt integration is currently modeled
  as metadata/projection state, not as native sequence-window materialization.
- `metadata: HashMap<String, serde_json::Value>`
  - includes in-memory candidate scoring/filter sets at
    `metadata["candidate_sets"]` (`gentle.candidate_sets.v1`)
  - on disk, candidate metadata may be a sidecar reference
    (`gentle.candidate_sets.ref.v1`) pointing to JSONL-indexed records
  - candidate sidecar load/save now uses stricter safety semantics:
    - atomic project-file writes
    - staged/rollback-aware sidecar replacement
    - strict-load opt-in via env (`GENTLE_CANDIDATE_STORE_STRICT_LOAD`)
      with non-strict warning fallback metadata
- `lineage: LineageGraph`
  - `nodes` (sequence lineage nodes with origin + creation op)
  - `seq_to_node` (current sequence id -> latest lineage node)
  - `edges` (parent -> child derivation edges, including multi-parent ligation)
  - `macro_instances` (shared-shell macro execution records with typed bound
    ports, emitted op ids, and status)
- analysis artifacts (dotplots/flexibility tracks) are persisted in
  `metadata["dotplot_analysis"]` and projected into GUI lineage as dedicated
  analysis nodes linked from source sequence lineage nodes
- `container_state: ContainerState`
  - `containers` (explicit singleton/pool/selection containers)
  - `seq_to_latest_container` (latest container membership index)
  - `next_container_counter` (stable container id generation)

Container semantics now exist as first-class state.
`SelectCandidate` remains explicit in-silico disambiguation when multiple
possible products exist.

Display state (`ProjectState.display`) now also carries TFBS and VCF
display-filter criteria so interactive GUI filtering and SVG export use one
shared contract.

### Operation-based execution

Clients submit operations (or workflows of operations).
Engine returns `OpResult` with:

- deterministic operation id (`op_id`)
- created/changed sequence ids
- warnings/messages

This enables:

- reproducibility
- auditability
- machine-to-machine use
- operation-level undo/redo and replay support

Progress/cancellation contract:

- Engine progress callbacks are cooperative and return `bool`.
- `true` continues the running operation.
- `false` requests cancellation.
- Long-running operations that support this contract (genome-track imports and
  genome-prepare flow) can stop early while returning warnings/progress
  summaries.
- Genome-prepare additionally supports explicit timeboxing (`timeout_seconds`
  in operation payload, `--timeout-secs` in CLI/shell, `timeout_sec` in GUI).
- Annotation parsing during genome preparation is fault-tolerant: malformed
  tabular annotation lines are skipped, summarized, and reported with capped
  file/line examples in warnings.

Interactive orchestration contract:

- Long-running jobs are surfaced in one GUI background-jobs panel with progress
  snapshots and cancel/retry controls where the engine callback contract
  supports cancellation.
- Operation history is surfaced in a dedicated GUI panel backed by engine
  operation journal + checkpoint stacks.

### Provenance and DAG (recommended)

Yes, every state should carry provenance metadata, and yes, this naturally
forms a DAG.

Proposed model:

- `ProjectState.metadata["provenance"]` stores:
  - `state_id`: deterministic or UUID id of this state snapshot
  - `parents`: zero or more parent `state_id`s
  - `created_at`: timestamp
  - `source_context`: optional descriptor (GUI session, CLI run, API caller)
  - `last_run_id`: workflow/run identifier when applicable
- Every applied operation appends an operation record edge:
  - `from_state_id` -> `to_state_id`
  - `op_id`, `run_id`, `operation payload hash`, `result summary`

Why this works with your CLI concern:

- GUI: can provide rich `source_context` (window/session/project metadata).
- CLI/agents: can omit human context and still provide machine context
  (`run_id`, tool name, caller id) when available.
- If no external context is provided, provenance still remains valid through
  operation edges and parent state ids.
- Current implementation now appends extraction-level provenance entries at
  `ProjectState.metadata["provenance"]["genome_extractions"]` for
  `ExtractGenomeRegion` and `ExtractGenomeGene`, including source descriptors,
  checksum fields, and anchor-strand context (`anchor_strand`) when available.
- `ExtractGenomeGene` now also auto-derives an exon-concatenated synthetic
  companion sequence (`<seq_id>__exons`) when transcript exon annotation is
  available; this derived sequence is lineage-linked to the extracted genomic
  gene sequence and intentionally not recorded as one contiguous genome anchor.

Practical rule:

- Context fields are optional, but lineage fields (`state_id`, `parents`,
  operation edge) are mandatory for DAG integrity.

## 5. Adapter responsibilities

### GUI

- Convert UI events to engine operations
- Render state and selections
- Never duplicate biology logic
- Provide searchable action routing (Command Palette) that dispatches only to
  existing engine/adaptor actions (no duplicate logic path).
- Keep anchored-data imports preflighted in UI (detected anchor, matching
  status, projected targets), while execution remains engine-owned.
- Current GUI-only routing note:
  - reference/helper status dialogs (including `Prepared References...`) are
    opened via menu/command-palette actions in `src/app.rs`
  - shared shell commands can query/prepare/extract, but cannot yet directly
    open/focus GUI dialogs

### JavaScript/Lua shells

- Provide ergonomic function wrappers
- Internally call engine operations
- Return structured results where possible
- Keep helper parity for state/display mutations (`set_parameter`,
  `set_vcf_display_filter`) so automation does not depend on GUI-only controls

### CLI (`gentle_cli`)

- Stable JSON interface for scripted and AI-driven workflows
- State import/export and summary
- Capabilities reporting
- Shared shell command path (`shell`) aligned with GUI Shell panel
- Structured command glossary contract:
  - source file: `docs/glossary.json`
  - consumed by runtime `help` rendering (`text|json|markdown`)
  - intended as single machine-readable command semantics index for
    CLI/GUI-shell/JS/Lua documentation generation
- Protocol-first workflow example contract:
  - canonical source files: `docs/examples/workflows/*.json`
  - schema: `gentle.workflow_example.v1`
  - snippet generator binary: `gentle_examples_docs`
  - generated adapter snippets: `docs/examples/generated/*.md`
  - canonical human-facing tutorial entry page: `docs/tutorial/README.md`
    - aggregates generated executable chapters, GUI walkthroughs, and
      agent-facing/reference tutorials
    - keeps tutorial type/status labels explicit so executable and
      hand-written material are distinguishable at discovery time
  - tutorial manifest source: `docs/tutorial/manifest.json`
  - tutorial schema: `gentle.tutorial_manifest.v1`
  - committed tutorial runtime outputs: `docs/tutorial/generated/`
  - tutorial generator/check commands:
    - `gentle_examples_docs tutorial-generate`
    - `gentle_examples_docs tutorial-check`
  - test gating by example metadata:
    - `always`: execute in default test runs
    - `online`: execute only with `GENTLE_TEST_ONLINE=1`
    - `GENTLE_SKIP_REMOTE_TESTS=1` force-disables remote-resource tests
      regardless of `GENTLE_TEST_ONLINE`
    - `skip`: parse/validate only
  - tutorial information-architecture rule:
    - users should have one canonical landing page for tutorials
    - executable chapters may remain generated, but manual walkthroughs and
      reference guides must be linked from the same entry point
    - distinction between tutorial, recipe, and reference material should be
      explicit in page metadata/catalog text
- Screenshot bridge status (temporarily disabled by security policy):
  - historical implementation existed as a compile-time + runtime gated adapter
    bridge (`screenshot-capture` feature + `--allow-screenshots` startup flag)
  - observed field behavior:
    - endpoint security (Cortex XDR) flagged related local test/process activity
      as malware
    - repeated warnings were seen for targeted test invocations of screenshot
      paths
  - current enforced behavior:
    - `screenshot-window` is rejected with a deterministic
      "disabled by security policy" message
    - `--allow-screenshots` is rejected by CLI/GUI startup parsing
    - screenshot-specific tests/invocation entry points were removed
    - rendering snapshot tests are now behind opt-in feature flag
      `snapshot-tests` and are disabled by default
  - immediate policy:
    - keep screenshot/image-capture execution disabled by default in this repo
    - keep auto-updated documentation with graphics postponed
  - how to re-enable later (intentional opt-in process):
    - obtain endpoint-security approval/exception first (project path, Rust test
      binaries, and host screenshot utility execution path)
    - revert screenshot security-disable gates in:
      - `src/engine_shell.rs` (`screenshot-window` parse/execute branches)
      - `src/bin/gentle.rs` and `src/bin/gentle_cli.rs`
        (`--allow-screenshots` argument handling and help text)
    - optionally restore screenshot-specific tests only after security exception
      is confirmed stable
    - rebuild with feature `screenshot-capture`
      (example: `cargo run --features screenshot-capture --bin gentle -- --allow-screenshots`)
    - run with startup flag `--allow-screenshots` for any process that should
      accept screenshot commands
    - validate from GUI shell first (active GENtle window only)
    - if snapshot tests are intentionally desired again, run with explicit
      feature opt-in (example:
      `cargo test --features snapshot-tests -q render_export::tests::snapshot_`)

#### Cross-platform screenshot relaunch (deferred)

This track is intentionally deferred until GUI stability improves and should be
kept as a design target, not an active implementation item.

Fixed decisions for relaunch:

- active GENtle window only for raster screenshot capture
- all GENtle windows are eligible targets (main, sequence, and auxiliary
  specialist windows)
- capture content area only (no native OS title bar/frame pixels)
- raster screenshot capture is GUI-context only
- CLI/MCP headless workflows use explicit SVG export commands instead of raster
  screenshot capture

Implementation direction:

- use a shared `egui/eframe` viewport screenshot event pipeline for all desktop
  platforms
- avoid OS shell screenshot utilities in the primary path
- keep screenshot capture behind explicit runtime/feature opt-in policy gates
  even after relaunch

### MCP server (implemented parity-expanded baseline)

Current baseline:

- `gentle_mcp` stdio server route is available.
- MCP serves both:
  - tool execution (`tools/call`)
  - capability discovery/negotiation (`tools/list`, `capabilities`, `help`)
- Exposed tools:
  - `capabilities`
  - `state_summary`
  - `op` (shared engine operation execution; explicit `confirm=true` required)
  - `workflow` (shared engine workflow execution; explicit `confirm=true` required)
  - `help`
  - `ui_intents`
  - `ui_intent`
  - `ui_prepared_genomes`
  - `ui_latest_prepared`
- Handlers are thin and map to existing shared contracts
  (`GentleEngine::capabilities`, state summary, `Engine::apply`,
  `Engine::apply_workflow`, glossary-backed help, and shared shell
  parser/executor for `ui ...` routes).
- Mutating MCP tools persist project state to disk after successful execution.
- Structured JSON-RPC diagnostics are returned for invalid requests/params.
- UI-intent MCP tools are explicitly non-mutating and reject unexpected
  `state_changed=true` results from routed shared-shell commands.
- MCP state access is file-backed (`state_path`) rather than live GUI-memory
  inspection; unsaved in-memory GUI edits are not exposed through MCP until
  persisted.
- Pointing MCP at the same project state file as active GUI/CLI workflows is
  an intentional shared-state model and should be treated as a trusted-client
  boundary.
- Adapter-equivalence tests cover MCP-vs-shared-shell parity for:
  `ui_intents`, `ui_prepared_genomes`, `ui_latest_prepared`, and `ui_intent`.

Remaining expansion scope:

- extend mutating tool breadth and additional shell-route coverage over shared
  paths.
- extend adapter-equivalence tests for additional MCP vs CLI shell flows.
- keep zero MCP-only biology logic branches.

Feature expert-view command baseline (implemented):

- `inspect-feature-expert` (shared shell) returns a structured expert view for
  a selected feature target.
- `render-feature-expert-svg` renders the same expert view to SVG via engine
  operation `RenderFeatureExpertSvg`.
- Target kinds currently supported:
  - TFBS feature by feature id
  - restriction site by cut position (with optional enzyme/recognition span
    hints)
- Direct CLI (`gentle_cli`) forwards both commands, and JS/Lua wrappers expose
  equivalent calls.

### GUI intent command plane (implemented baseline)

Current baseline:

- Shared shell `ui ...` commands now exist:
  - `ui intents`
  - `ui open TARGET ...`
  - `ui focus TARGET ...`
  - `ui prepared-genomes ...`
  - `ui latest-prepared SPECIES ...`
- GUI-side intent handlers in `src/app.rs` now map `ui open|focus` intents to
  existing dialog openers (Prepared References, prepare/retrieve/blast, track
  import, agent assistant, helper-genome dialogs).
- UI-intent capability/introspection output is available via `ui intents`.
- Query helpers are implemented and can be composed with open/focus for
  prepared-reference selection, for example:
  - `ui open prepared-references --species human --latest`
  - explicit `--genome-id` overrides query-based selection

Remaining additions for full text/voice GUI control:

1. Keep destructive operation confirmation semantics explicit when driven by
   text/voice commands.
2. Continue intent-grammar hardening for broader natural-language-to-intent
   mapping while preserving deterministic execution.
3. Add optional voice transport (STT/TTS) that emits/consumes the same
   deterministic intent contracts through the existing agent interface.
4. Extend the intent plane with explicit clipboard/export intents so
   cross-application copy/paste can be routed through auditable engine
   operations rather than frontend-only ad hoc handlers.

UI-intent tool routine (target contract for MCP/agent/voice adapters):

1. Discover:
   - caller requests current intent capability catalog (`ui intents` equivalent)
   - adapter returns stable action/target list and required/optional arguments
2. Resolve:
   - caller resolves target selection deterministically (for example
     `ui prepared-genomes` / `ui latest-prepared`) before open/focus execution
   - ambiguous query results must return structured "needs disambiguation"
     payloads instead of guessing
3. Guard:
   - mutating/destructive intents require explicit confirmation field in the
     adapter request before execution
   - non-mutating intents (`open`, `focus`, list/query helpers) execute without
     mutating confirmation
4. Execute:
   - adapter maps request to existing shared `ui ...` command contracts
     (single parser/executor path; no frontend-local behavior forks)
   - GUI host performs only thin routing to existing window/dialog openers
5. Report:
   - adapter returns structured result with `executed`, `resolved_target`,
     warnings/errors, and optional follow-up suggestions (for example candidate
     target list for user choice)
   - result payloads stay machine-readable and adapter-equivalent

Voice path note:

- Biology-by-voice should route through the same agent interface used for text
  suggestions/execution, with voice as transport only
  (STT -> agent prompt -> deterministic shell/UI intents -> TTS).
- No biology logic should move into voice-specific code.

## 6. Protocol-first direction

GENtle should expose versioned machine contracts:

- operations
- state
- results
- errors
- capabilities

Why:

- makes AI/automation robust
- simplifies future integrations (including OpenClaw)
- reduces coupling to Rust internals

## 7. AI-tooling requirement

GENtle should be usable as a skill/tool by agents.

Minimum requirements:

1. Deterministic operation endpoints
2. Structured errors and warnings
3. State persistence and export
4. Capability negotiation
5. Stable versioning policy

ClawBio/OpenClaw compatibility direction:

- Present GENtle to ClawBio as a thin skill wrapper around deterministic
  `gentle_cli` command surfaces (not by duplicating biology logic in the skill).
- Wrapper outputs should remain reproducibility-first:
  `report.md`, `result.json`, and reproducibility artifacts
  (`commands.sh`, `environment.yml`, checksums).
- Keep command contracts explicit in request payloads so orchestration layers
  can route high-level asks while execution remains replayable.

Current work satisfies (1) through (4) for most operations; remaining gaps are
mainly view-model formalization and promoting remaining adapter-level utility
contracts into stable engine operations.

MCP direction:

- MCP is the preferred standardized route for external AI tool orchestration and
  capability negotiation.
- `agents ask` remains supported as a complementary route for chat-like
  assistance and suggested-command execution.

### Agent assistant bridge (implemented)

GENtle now provides a shared agent-assistance bridge across GUI and CLI shell:

- Shared shell commands:
  - `agents list [--catalog PATH]`
  - `agents ask SYSTEM_ID --prompt TEXT [--catalog PATH] [--base-url URL] [--model MODEL] [--timeout-secs N] [--connect-timeout-secs N] [--read-timeout-secs N] [--max-retries N] [--max-response-bytes N] [--allow-auto-exec] [--execute-all] [--execute-index N ...] [--no-state-summary]`
- Catalog source defaults to `assets/agent_systems.json`.
- Catalog entries describe transport and invocation details:
  - `builtin_echo` (offline/demo transport)
  - `external_json_stdio` (external adapter command over stdin/stdout JSON)
  - `native_openai` / `native_openai_compat` (built-in HTTP adapters with optional per-request base URL/model override)
  - `native_openai_compat` is the preferred path for OpenAI-compatible hosted
    providers (including [clawbio.ai](https://clawbio.ai/)) through explicit
    base-URL + model selection
  - `native_openai_compat` requires explicit model resolution (catalog model or
    per-request override; `unspecified` is rejected) and keeps endpoint host/port
    deterministic (no hidden host fallback probing)
  - runtime guardrails are overrideable per request (timeout/connect/read,
    retry budget, max response bytes) and surfaced in invocation telemetry
- Agent request payload includes:
  - system id and prompt
  - optional project state summary context
- Agent response payload supports mixed outcomes per reply:
  - plain assistant message (`chat`)
  - follow-up questions (`ask`)
  - suggested shell commands with explicit execution intent (`ask` or `auto`)

Execution safety model:

- There is no global always-execute mode.
- Execution is evaluated per returned suggestion:
  - explicit user-run by index (`--execute-index`, GUI per-row Run)
  - bulk explicit run (`--execute-all`)
  - auto-only when caller enables `--allow-auto-exec` and suggestion intent is
    `auto`
- Nested/recursive `agents ask` execution is blocked in suggested-command runs.
- Suggested commands are executed through the same shared shell parser/executor
  used by GUI shell and `gentle_cli shell`.
- This includes BLAST shell routes (`genomes blast`, `helpers blast`,
  `genomes blast-track`, `helpers blast-track`) when agents suggest them.
- Async BLAST job contract baseline is now available through shared shell:
  - `genomes/helpers blast-start`
  - `genomes/helpers blast-status`
  - `genomes/helpers blast-cancel`
  - `genomes/helpers blast-list`
- Async BLAST execution is now scheduler-backed (bounded FIFO queue with
  configurable max-concurrency), exposing explicit `queued`/`running` states
  before terminal outcomes.
- The same async BLAST surface is exposed through MCP tools
  (`blast_async_start|status|cancel|list`) with shared parser/executor parity.
- Remaining limitation:
  - agent auto-execution currently runs one suggested command at a time and does
    not yet orchestrate multi-step poll/wait loops automatically for async jobs.
  - planned primer-pair multi-BLAST specificity fan-out still needs dedicated
    workflow-level async orchestration on top of the baseline job primitives.
- Protocol-level JSON shapes and exact command contract are specified in
  `docs/protocol.md` (agent catalog/request/response schemas and execution
  intent semantics).

Minimal-success rollout profile (recommended):

- Start with catalog-driven stdio adapters (for example OpenAI bridge command)
  rather than ad-hoc RAG/training.
- Keep default interaction mode at `chat` + `ask`; enable `auto` only for
  tightly allowlisted low-risk commands.
- Include compact machine context (`state_summary`) instead of large free-form
  project dumps to keep prompts deterministic and understandable.
- For local/small models, provide a domain bootstrap package:
  - `docs/ai_cloning_primer.md`
  - `docs/ai_task_playbooks.md`
  - `docs/ai_prompt_contract.md`
  - `docs/examples/ai_cloning_examples.md`
  - optional term extension: `docs/ai_glossary_extensions.json`

### Candidate query/optimization command contract (current)

Current implementation is engine-level and exposed through:

- first-class `gentle_cli candidates ...`
- shared shell (`gentle_cli shell` + GUI shell panel)
- dedicated Engine Ops GUI candidate panel

Command surface:

- `candidates list`
- `candidates show SET_NAME [--limit N] [--offset N]`
- `candidates metrics SET_NAME`
- `candidates delete SET_NAME`
- `candidates generate SET_NAME SEQ_ID --length N ...`
- `candidates generate-between-anchors SET_NAME SEQ_ID --length N ...`
- `candidates score SET_NAME METRIC_NAME EXPRESSION`
- `candidates score-distance SET_NAME METRIC_NAME ...`
- `candidates score-weighted SET_NAME METRIC_NAME --term METRIC:WEIGHT[:max|min] ...`
- `candidates top-k INPUT_SET OUTPUT_SET --metric METRIC --k N ...`
- `candidates pareto INPUT_SET OUTPUT_SET --objective METRIC[:max|min] ...`
- `candidates filter INPUT_SET OUTPUT_SET --metric METRIC ...`
- `candidates set-op union|intersect|subtract LEFT RIGHT OUTPUT`
- `macros run [--transactional] [--file PATH | SCRIPT_OR_@FILE]`
- `macros template-list|template-show|template-put|template-delete|template-import|template-run ...`
- `routines list [--catalog PATH] [--family NAME] [--status NAME] [--tag TAG] [--query TEXT]`
- `routines explain ROUTINE_ID [--catalog PATH]`
- `routines compare ROUTINE_A ROUTINE_B [--catalog PATH]`
- `candidates macro [--transactional] [--file PATH | SCRIPT_OR_@FILE]`
- `candidates template-list|template-show|template-put|template-delete|template-run ...`
- `set-param NAME JSON_VALUE`
- Starter cloning-pattern pack: `assets/cloning_patterns.json` (`gentle.cloning_patterns.v1`)
- Hierarchical per-template cloning-pattern catalog:
  `assets/cloning_patterns_catalog/**/*.json`
  (`gentle.cloning_pattern_template.v1`)
- `macros template-import PATH` supports:
  - single template file
  - legacy pack file
  - recursive directory import
- Typed cloning-routine manifest baseline:
  `assets/cloning_routines.json` (`gentle.cloning_routines.v1`)
- First Gibson + restriction family baselines are now shipped:
  - routine: `gibson.two_fragment_overlap_preview`
  - template:
    `assets/cloning_patterns_catalog/gibson/overlap_assembly/gibson_two_fragment_overlap_preview.json`
  - routine: `restriction.digest_ligate_extract_sticky`
  - template:
    `assets/cloning_patterns_catalog/restriction/digest_ligation/digest_ligate_extract_sticky.json`
- GUI `Patterns` menu mirrors the catalog folder hierarchy and routes imports
  through the same shared-shell command path.
- GUI `Patterns` menu now also exposes routine catalog discovery grouped by
  family/status and routine-linked template import actions.
- Multiple template runs in one project are supported today via normal
  operation execution and persistence.
- Macro-instance lineage nodes are now first-class in graph/SVG rendering and
  record success/failure/cancel status with optional status message.
- Shared-shell macro introspection contracts are available via:
  - `macros instance-list`
  - `macros instance-show MACRO_INSTANCE_ID`
- Shared-shell preflight now supports routine-family semantic hooks in addition
  to generic typed-port checks:
  - Gibson routines validate adjacent sequence overlap compatibility against
    declared overlap length (`overlap_bp`) before mutation.
  - Restriction routines validate enzyme resolution and digest-sanity semantics
    (distinct enzyme inputs, recognition-site presence, fragment-index and
    extract-range checks) before mutation.

Concrete patch plan: routine-application assistant and alternative-awareness

- Extend routine metadata (`gentle.cloning_routines.v1`) with additive
  explanation fields used by both humans and AI:
  - `purpose`, `mechanism`, `requires[]`, `contraindications[]`
  - `confusing_alternatives[]` (routine ids)
  - `difference_matrix[]` (fixed comparison axes)
  - `disambiguation_questions[]`
  - `failure_modes[]`
- Add non-mutating shared-shell explainability routes:
  - `routines explain ROUTINE_ID`
  - `routines compare ROUTINE_A ROUTINE_B`
- Keep these routes adapter-equivalent (GUI/CLI/JS/Lua/MCP consume the same
  JSON payload; no adapter-only explanation logic branches).
- Add a GUI Routine Assistant that applies routines in deterministic stages:
  - goal capture -> ranked candidates -> compare (`why this / why not`) ->
    typed parameter form -> preflight preview -> transactional run -> export.
- Keep preflight semantics authoritative for execution gating; explanation data
  is guidance and decision context, not an execution shortcut.
- Extend protocol/run-bundle export with optional routine decision trace:
  selected routine, considered alternatives, disambiguation answers, and
  preflight summary at decision time.
- Add deterministic parity tests for explain/compare payloads and assistant
  stage mapping (`routine -> form -> preflight -> execution`).

Meta-planning layer (time/cost/local-fit) and lab-management agent constraints:

- Scope boundary (v1):
  - planning affects recommendation/ranking only,
  - execution semantics of cloning operations remain unchanged.
- Engine-owned planning contracts:
  - `gentle.planning_profile.v1`
  - `gentle.planning_objective.v1`
  - `gentle.planning_estimate.v1`
  - `gentle.planning_suggestion.v1`
  - `gentle.planning_sync_status.v1`
- Effective profile merge order is deterministic:
  - `global_profile -> confirmed_agent_overlay -> project_override`
- Routine-ranking payloads include planning estimates:
  - `estimated_time_hours`, `estimated_cost`, `local_fit_score`,
    `composite_meta_score`, and machine-readable explanation details.
- Purchasing simplification (explicit v1 design decision):
  - missing required local material class adds procurement delay,
  - default delay is `10` business days (`procurement_business_days_default`),
  - item-level override supported (`inventory.<class>.procurement_business_days`),
  - delay is applied once per missing required class (deduplicated),
  - business-day model currently excludes weekends only (no holiday calendar),
  - business-day delays map to `estimated_time_hours` using deterministic
    weekend-aware conversion (`24h * 7/5` per business day).
- Lab-management agent integration policy:
  - sync mode supports both pull and push suggestion capture,
  - all incoming/outgoing agent updates are persisted as `pending` suggestions
    with source/confidence/snapshot metadata,
  - updates are advisory; activation requires explicit user accept/reject
    action (no auto-apply).

Decision-trace capture/export contract (detailed plan):

- Why this exists:
  - Routine execution provenance currently captures what ran, but not why one
    routine was selected over alternatives.
  - `decision_trace` closes that gap by serializing user/agent intent,
    alternatives considered, preflight gating, and final execution outcome in
    one deterministic record.
- Non-goals for phase 1:
  - No autonomous recommendation engine.
  - No hidden ranking logic in adapters.
  - No mutation of execution behavior based on trace content.

- Scope (phase 1):
  - Capture routine-selection intent and gating outcomes from GUI
    `Routine Assistant` runs.
  - Persist traces in project metadata while a session is active.
  - Include the finalized trace payload in `ExportProcessRunBundle` output.
  - Keep all fields additive and optional for backward-compatible loading.

- Lifecycle capture points (ordered):
  - `assistant_started`: assistant session ID allocated, source adapter recorded.
  - `intent_captured`: goal text/query text plus routine-list filters recorded.
  - `candidates_shown`: deterministic ordered candidate IDs snapshot recorded.
  - `routine_selected`: selected routine ID + presented alternatives snapshot.
  - `comparison_opened` (0..N): each compare pair and comparison payload hash.
  - `disambiguation_answered` (0..N): ordered answer commits, no inferred values.
  - `bindings_committed`: typed input-port bindings frozen for preflight.
  - `preflight_evaluated` (1..N): latest preflight snapshot retained as
    canonical gate state, prior snapshots retained in history.
  - `execution_submitted` (0..1): transactional flag and submission timestamp.
  - `execution_finished` (0..1): success/failure plus macro-instance and op IDs.
  - `trace_exported` (0..N): run-bundle path(s) and export timestamps.

- Trace object schema (planned):
  - schema: `gentle.routine_decision_trace.v1`
  - identity:
    - `trace_id` (deterministic per assistant run)
    - `source` (`gui_routine_assistant`, then `cli/js/lua`)
    - `status` (`draft`, `preflight_failed`, `ready`, `executed`,
      `execution_failed`, `aborted`, `exported`)
    - `created_at_unix_ms`, `updated_at_unix_ms`
  - intent and search context:
    - `goal_text`
    - `query_text`
    - `candidate_snapshot`:
      - applied list/search filters
      - ordered candidate routine IDs displayed to the user
      - optional `catalog_version`/`catalog_hash`
  - selection and alternatives context:
    - `selected_routine` (id/title/family/status/template)
    - `alternatives_presented` (ordered routine IDs)
    - `comparisons` (ordered):
      - `left_id`, `right_id`
      - shared/left-only/right-only tags
      - aligned `difference_matrix` rows consumed by UI
      - optional payload hash for determinism/debug
  - disambiguation context:
    - `disambiguation_questions_presented` (ordered IDs/text)
    - `disambiguation_answers` (ordered by question ID)
    - free-text allowed now; typed answer schema is additive in later version
  - parameter and gating context:
    - `bindings_snapshot` (typed port -> bound value map)
    - `preflight_history` (ordered snapshots)
    - `preflight_snapshot` (last snapshot, canonical):
      - `can_execute`
      - `warnings[]`
      - `errors[]`
      - `contract_source`
      - optional checked-port summary hash/count
  - execution context:
    - `execution_attempted`
    - `execution_success`
    - `transactional`
    - `macro_instance_id` (when available)
    - `emitted_operation_ids` (ordered)
    - `execution_error` (structured summary if failed)
  - export context:
    - `export_events[]`:
      - `run_bundle_path`
      - `exported_at_unix_ms`
      - optional `bundle_hash`

- Run-bundle placement (planned):
  - `gentle.process_run_bundle.v1` includes top-level
    `decision_traces[]` (ordered by `created_at_unix_ms` then `trace_id`).
  - Bundle export can include:
    - all traces (default for reproducibility)
    - selected trace by ID (optional future export filter)

- Determinism and normalization rules:
  - Stable ordering is mandatory for candidates, comparisons, questions,
    answers, preflight snapshots, and emitted operation IDs.
  - Missing optional values are represented explicitly as `null` or empty arrays
    (no adapter-specific field dropping).
  - Text normalization at capture time:
    - trim outer whitespace
    - preserve user-entered interior whitespace
    - preserve case (no adapter-side rewriting)
  - Trace ID generation must be deterministic relative to assistant session
    identity and monotonic local counter, not wall-clock randomness.

- Failure/cancel behavior:
  - If execution is never submitted, export partial trace with
    `execution_attempted=false` and current status (`draft`,
    `preflight_failed`, or `aborted`).
  - If preflight fails, trace remains exportable and includes complete
    `preflight_snapshot` + history.
  - If execution fails, trace includes structured failure summary and emitted
    operation IDs up to failure boundary.

- Retention and redaction policy (phase 1):
  - Trace fields must avoid secret capture by default; only user-visible
    assistant inputs and deterministic engine outputs are recorded.
  - File paths may be stored when they are explicit binding inputs; no
    automatic expansion of unrelated environment paths.
  - Additive redaction flags may be introduced later, but default export is
    full trace for reproducibility.

- Adapter-parity target:
  - GUI is phase-1 producer.
  - CLI/JS/Lua should submit equivalent payloads through the same schema once
    routine-assistant flows are exposed there.
  - Adapters may differ in UX, but not in serialized field meaning or
    ordering semantics.

- Validation and test expectations:
  - deterministic snapshot tests for success/failure/aborted traces
  - round-trip load/save test preserving ordering and null/empty field shape
  - run-bundle export tests verifying trace inclusion and stable serialization
  - parity tests where the same staged decisions from different adapters yield
    schema-equivalent trace payloads.

This enables reusable query composition:

1. generate explicit candidate sets
2. attach base/derived scores
3. run explicit optimizer primitives (weighted objective, top-k, Pareto frontier)
4. filter by absolute threshold and/or quantile
5. intersect/union/subtract with other candidate sets

### Primer/qPCR design report command contract (baseline)

Primer-pair and qPCR design are first-class engine operations plus shared-shell
inspection/export paths:

- Engine operation:
  - `DesignPrimerPairs { template, roi_start_0based, roi_end_0based, forward, reverse, pair_constraints?, min_amplicon_bp, max_amplicon_bp, max_tm_delta_c?, max_pairs?, report_id? }`
  - `DesignQpcrAssays { template, roi_start_0based, roi_end_0based, forward, reverse, probe, pair_constraints?, min_amplicon_bp, max_amplicon_bp, max_tm_delta_c?, max_probe_tm_delta_c?, max_assays?, report_id? }`
  - `forward`/`reverse` side constraints now include optional sequence-level filters:
    `fixed_5prime`, `fixed_3prime`, `required_motifs[]`, `forbidden_motifs[]`,
    and `locked_positions[]` (offset/base locks, IUPAC-aware).
- Persisted metadata:
  - `primer_design_reports` (`gentle.primer_design_reports.v1`)
  - `qpcr_design_reports` (`gentle.qpcr_design_reports.v1`)
  - report payload schema: `gentle.primer_design_report.v1`
  - report payload schema: `gentle.qpcr_design_report.v1`
- Shared-shell commands:
  - `primers design REQUEST_JSON_OR_@FILE [--backend auto|internal|primer3] [--primer3-exec PATH]`
  - `primers design-qpcr REQUEST_JSON_OR_@FILE [--backend auto|internal|primer3] [--primer3-exec PATH]`
  - `primers seed-from-feature SEQ_ID FEATURE_ID`
  - `primers seed-from-splicing SEQ_ID FEATURE_ID`
  - `primers list-reports`
  - `primers show-report REPORT_ID`
  - `primers export-report REPORT_ID OUTPUT.json`
  - `primers list-qpcr-reports`
  - `primers show-qpcr-report REPORT_ID`
  - `primers export-qpcr-report REPORT_ID OUTPUT.json`

Adapter contract:

- `gentle_cli` exposes this surface directly (`gentle_cli primers ...`) by
  forwarding into the same shared-shell parser/executor path.
- JS/Lua adapters consume the same operation contract via `apply_operation`
  and can use the same shared-shell routes for report inspection/export.
- No adapter-specific primer-design business logic is permitted.
- ROI seed helpers are intentionally non-mutating and return seed payload
  schema `gentle.primer_seed_request.v1` with ready-to-run operation JSON for
  both pair-PCR and qPCR design.

Primer3 wrapper integration status (baseline implemented):

- `DesignPrimerPairs` remains the canonical external contract for
  GUI/CLI/JS/Lua/agent/MCP entry points; Primer3 integration is an engine
  backend implementation detail.
- Engine backend selection is now explicit and deterministic via parameters:
  - `primer_design_backend = auto|internal|primer3`
  - `primer3_executable` (default `primer3_core`)
- `auto` mode attempts Primer3 and deterministically falls back to the internal
  backend with explicit warning text when Primer3 is unavailable.
- Primer3 outputs are normalized into the existing
  `gentle.primer_design_report.v1` pair schema with shared deterministic
  ranking/tie-break behavior.
- Reports now carry backend provenance metadata:
  - requested backend
  - used backend
  - optional fallback reason
  - optional Primer3 executable/version details
- Remaining work:
  - deeper Primer3 constraint mapping parity for edge-case constraints
  - broader fixture-backed equivalence matrix between internal and Primer3
  - dedicated preflight/status views in GUI configuration panels
  - multi-BLAST specificity tiers integrated into primer-pair post-filtering

### Primer/PCR/BLAST UI and internal BLAST abstraction plan (new)

Goal: present coherent user interfaces for single-primer design, primer-pair
design, and BLAST, while enforcing one engine-owned BLAST abstraction reused by
all three pathways.

Detailed implementation plan: `docs/primer_design_specialist_plan.md`.

Current status (2026-03-01):

- Engine now implements layered BLAST option resolution and validation, with
  persisted request/effective option provenance in BLAST reports/import history.
- Shared-shell/CLI and JS/Lua routes can now submit per-request BLAST JSON
  option overrides.
- GUI now includes a BLAST options section with quick controls, presets,
  advanced JSON editing, and engine-backed effective-options preflight preview.

Internal abstraction (engine-owned, adapter-neutral):

- Add a shared BLAST service boundary used by GUI/CLI/JS/Lua/MCP-facing flows:
  - request shape includes:
    - query payload(s)
    - target genome/helper selection
    - `task`
    - JSON options object (`options_json`) for advanced thresholds/controls
  - response shape includes:
    - parsed hits
    - warnings/errors
    - invocation/provenance
    - resolved effective options (post-default merge)
  - progress shape includes:
    - phase/status text
    - query-level done/total counts
    - elapsed time
- Keep deterministic option resolution in engine:
  - built-in defaults
  - optional defaults file
  - project-level override operation
  - per-request JSON override
  - strict validation for unknown keys and type/range mismatches
- Primer and primer-pair scoring/specificity steps must call this service, not
  invoke BLAST binaries directly from adapter/UI code.

BLAST specialist UI direction:

- Keep the dedicated BLAST viewport and evolve it to a five-block structure:
  1. Input (`manual`, `project sequence`, `project pool`)
  2. Target (`prepared genome/helper`)
  3. Options (quick controls + advanced JSON editor + presets)
  4. Execution (run/cancel/retry with shared progress contract)
  5. Results (hits table, import-to-track, export JSON, invocation detail)
- Show both invocation template (pre-run) and resolved invocation (post-run)
  consistently in status/results.

Single-primer design UI direction:

- Add a dedicated primer-design specialist window, sequence-context aware.
- Required sections:
  - target region (selection/feature/manual coordinates)
  - primer constraints (length/GC/Tm/3' rules/homopolymer/ambiguity)
  - specificity policy (none/quick/strict) backed by shared BLAST service
  - ranked result list with per-primer diagnostics
- Provide actions to annotate selected primer(s), export sets, and forward
  selected primers into PCR workflows.

Primer-pair (PCR) design UI direction:

- Keep `DesignPrimerPairs` as canonical engine operation and expand GUI around
  explicit pair-design workflows (amplicon intent + pair constraints +
  specificity tier + ranked reports).
- Pair-specificity stages that fan out into multiple BLAST checks must use the
  same async progress/cancel semantics as standalone BLAST.
- Selected pair handoff into `Pcr`/`PcrAdvanced`/`PcrMutagenesis` should be
  explicit and deterministic (no adapter-local hidden transformations).

Parity/invariant requirements:

- No adapter-specific biology logic for primer/blast scoring.
- Same request/response/progress/provenance semantics across GUI/CLI/JS/Lua.
- Operation history and lineage summaries must preserve BLAST invocation and
  effective-option provenance for primer and primer-pair specificity checks.

Local anchor model (for candidate/extraction workflows):

- Local sequence anchors are an in-sequence abstraction (`SequenceAnchor`) and
  are intentionally separate from genome-provenance anchoring metadata used by
  `ExtractGenomeRegion`/`ExtractGenomeGene`/`ExtendGenomeAnchor`.
- Current local anchor forms: absolute `Position` and `FeatureBoundary`
  (`Start|End|Middle`).

Feature-distance geometry controls (engine + shell/CLI):

- `--feature-geometry feature_span|feature_parts|feature_boundaries`
- `--feature-boundary any|five_prime|three_prime|start|end`
- `--strand-relation any|same|opposite`
- defaults preserve prior behavior (`feature_span` + `any`)
- `--feature-boundary` is only meaningful with `feature_boundaries`

## 8. Rendering and interpretation direction

Two additional layers are needed for full human/AI collaboration:

1. `state -> view model -> rendering`
2. `human drawing/image -> interpreted state proposal` (later)

Important separation:

- Engine remains deterministic and biology-centric
- Rendering remains presentation-centric
- Image interpretation remains probabilistic and should output proposals,
  not direct mutations without confirmation

### Feature expert view contract (implemented baseline)

Expert plots for selected features are engine-owned view models, not GUI-local
recomputations. This keeps GUI/CLI/JS/Lua behavior aligned.

Shared view model:

- `FeatureExpertTarget`: stable target selector
  (TFBS feature, restriction site, splicing-seeded feature, or isoform panel).
- `FeatureExpertView`: typed payload with one of:
  - `TfbsExpertView`
  - `RestrictionSiteExpertView`
  - `SplicingExpertView`
  - `IsoformArchitectureExpertView`

TFBS expert semantics:

- One motif column per PSSM position.
- Column bar total height is information content:
  - `IC_j = 2 + sum_b (p_j,b * log2(p_j,b))`
  - convention: `0 * log2(0) = 0`
- Nucleotide segments stack within each bar by their relative frequencies
  (`p_j,A/C/G/T`).
- The matched sequence instance is overlaid as a polyline crossing the chosen
  base segment in each column.
- Instruction text is part of the shared model so GUI and exported SVG can show
  consistent explanatory wording.

Restriction-site expert semantics:

- Render a compact top/bottom strand view for the recognized site context.
- Mark cleavage position explicitly.
- Include enzyme/site metadata and shared instruction text.

Splicing expert semantics:

- Group transcripts by shared locus/gene context and strand.
- Render one transcript lane per isoform on a shared coordinate axis.
- Keep exon/intron geometry coordinate-true (labels never resize feature
  geometry).
- Mark donor/acceptor boundaries with canonical/non-canonical motif context.
- Summarize alternative-splicing events and transcript-vs-exon membership in a
  matrix.
- Render junction-support arcs from the same engine payload used by GUI and
  SVG export.

Isoform-architecture expert semantics:

- Import is panel-driven and explicit:
  - curated JSON resource (`gentle.isoform_panel_resource.v1`)
  - imported via `ImportIsoformPanel { seq_id, panel_path, panel_id?, strict }`
- Rendering model is two-section and deterministic:
  - transcript/exon architecture lanes on genomic coordinates
  - protein/domain architecture lanes on amino-acid coordinates
- One stable row order is shared between transcript and protein sections.
- Domain/transactivation annotations come from the curated panel resource, not
  from GUI-local inference.
- Transcript mapping warnings are part of the shared payload, so GUI, shell,
  CLI, and SVG export report identical mapping diagnostics.

Export semantics:

- `RenderFeatureExpertSvg` is the canonical engine export path.
- `RenderIsoformArchitectureSvg` is the canonical panel-specific export path.
- GUI detail panel and shell/CLI scripting consume the same expert view
  contract before/while exporting.
- SVG output is therefore adapter-equivalent by construction.

### Alternative-splicing interpretation contract (implemented baseline)

Alternative-splicing interpretation now follows the same deterministic
engine-owned view-model pattern used by other feature expert views.

Shared splicing view model:

- `FeatureExpertTarget::SplicingFeature`: target selector.
- `SplicingExpertView`: typed payload carrying:
  - transcript lanes over one shared genomic axis
  - exon/intron geometry per transcript
  - splice-boundary markers (donor/acceptor)
  - event summaries (for example exon skip, alt 5', alt 3', mutually exclusive,
    intron retention)
  - optional junction-evidence arcs and support values
  - optional transcript-vs-exon presence matrix.

Rendering invariants:

- Biological geometry is coordinate-driven:
  - exon rectangles are sized strictly by genomic coordinates.
  - labels must never expand exon/intron footprints or hide neighboring
    exon-intron structures.
- One lane per transcript in splicing mode with stable lane ordering.
- Introns are rendered as connectors only; exon presence/absence carries the
  primary interpretation.

Interaction and export semantics:

- GUI splicing mode, SVG export, and scripting interfaces must consume the same
  `SplicingExpertView` contract.
- Boundary/event summaries must therefore be adapter-equivalent by design.
- Implemented operations mirror the expert-view pattern:
  - inspect splicing view payload (`inspect-feature-expert ... splicing ...`)
  - render splicing view SVG from the same payload
    (`render-feature-expert-svg ... splicing ...`).

### RNA-seq evidence direction for cloning-candidate regions (planned)

Primary near-term scope is small genomic region interpretation for cloning
candidate workflows, not whole-transcriptome alignment inside GENtle.

Detailed implementation plan:
`docs/rna_seq_nanopore_cloning_regions_plan.md`.

Detailed follow-up plan for multi-gene origin classification and sparse
annotated seed indexing:
`docs/rna_read_origin_sparse_index_plan.md`.

Architecture constraints for this track:

- RNA evidence handling must remain engine-owned and deterministic; adapters
  only invoke/display shared contracts.
- Initial RNA-seq path is Nanopore cDNA evidence ingestion with deterministic
  filtering and partial ROI mapping overlays.
- Read-orientation handling must stay explicit and deterministic:
  cDNA-oriented poly-T prefix reverse-complement normalization is configurable,
  and direct-RNA mode must remain available without implicit strand mixing.
- For `all-overlap / both-strands` scope, RNA-read interpretation should remain
  a single combined run, but assignment decisions must include
  strand-partitioned diagnostics and deterministic tie-breaks so cross-strand
  seed reuse is always auditable.
- Deferred mapping pass for this track should remain seed-hash based (exon-only
  then exon-exon junction then intronic fallback) instead of mandatory
  dynamic-programming alignment, so behavior stays deterministic and scalable on
  large long-read batches.
- Batch-oriented outputs must remain engine-owned: sample-sheet export contracts
  (TSV + machine-readable frequency columns) are shared across GUI/CLI/agents
  and derived from persisted RNA-read reports.
- Per-read exon-path and exon/transition abundance TSV exports should remain
  first-class engine contracts (not GUI-only derivations) so downstream cohort
  analysis is adapter-equivalent.
- Fixed seed-hit thresholds are bootstrap defaults for phase-1 only; the
  architecture must support a transcriptome-scale signal-to-noise model that can
  derive dynamic acceptance thresholds from empirical background distributions.
- Suggestion-first curation applies: RNA-derived edit candidates are
  non-mutating until explicitly user-confirmed through shared engine edit
  operations.
- Cross-adapter parity (GUI/CLI/shared shell) is required for summary
  inspection/export surfaces.
- Long-running RNA-read interpretation should avoid exclusive engine-write
  lock residency during compute-heavy read parsing/filtering so adapters stay
  responsive while progress streams are active.
- SNR/background diagnostics must be persisted in engine-owned report payloads
  so GUI/CLI/shared shell/agent adapters can explain decisions consistently.
- RNA-read seed evidence should remain reusable outside interpretation-only
  views; primer-design workflows may consume persisted seed-support maps as
  deterministic guidance input (candidate-region ranking), but primer scoring
  contracts must remain engine-owned and adapter-equivalent.

Ownership/defer boundary:

- Transposon-specific interpretation and modeling decisions are intentionally
  deferred pending Anze Karlek's direction.
- Seed-capture workflow abstraction (biotech analogy for seed-hash filtering) is
  intentionally deferred for a follow-up track; when implemented, it should be
  a reusable engine operation usable in generic workflows, not only a splicing
  expert window action (tracked in
  `docs/rna_seq_nanopore_cloning_regions_plan.md`).

### Dotplot + promoter-flexibility contract (implemented baseline; follow-ups)

Goal:

- Add low-latency sequence-self and sequence-vs-sequence analysis in sequence
  windows without coupling biological logic to GUI rendering.
- Support promoter-oriented interpretation (repeat structure + local
  flexibility proxies) through engine-owned, exportable view models.

Design constraints:

- Dotplot must not block first paint of sequence windows.
- Dotplot/flexibility computation must be cancellable and progress-reported.
- GUI/CLI/JS/Lua/SVG must consume the same engine payloads.
- Feature geometry remains authoritative; analysis overlays are additive.

Status (2026-03-19):

- Implemented baseline:
  - engine operations `ComputeDotplot` and `ComputeFlexibilityTrack`
  - persisted analysis payloads in project metadata (`dotplot_analysis`)
  - shared-shell/CLI read/write command surfaces:
    - `dotplot compute|list|show|render-svg`
    - `render-dotplot-svg`
    - `flex compute|list|show`
  - engine SVG export parity:
    - `RenderDotplotSvg { seq_id, dotplot_id, path, flex_track_id?, display_density_threshold?, display_intensity_gain? }`
    - shared-shell/CLI command route: `render-dotplot-svg ...`
    - GUI `Export Dotplot SVG...` now uses the same engine operation (no
      GUI-only renderer path for persisted exports)
  - lineage projection parity:
    - SVG export operations are materialized as analysis nodes in the main
      lineage table/graph and linked to source sequences via operation edges
  - adapter convenience wrapper parity:
    - JS: `render_dotplot_svg(...)`
    - Lua: `render_dotplot_svg(...)`
    - Python: `render_dotplot_svg(...)`

Engine objects and operations:

- Dotplot payload schema:
  - `gentle.dotplot_view.v2`
  - typed fields:
    - `seq_id`
    - `reference_seq_id?` (self mode keeps this `null`)
    - `span_start_0based`, `span_end_0based`
    - `reference_span_start_0based`, `reference_span_end_0based`
    - `mode` (`self_forward`, `self_reverse_complement`, `pair_forward`,
      `pair_reverse_complement`)
    - `word_size`, `step_bp`, `max_mismatches`
    - sparse/packed match points in deterministic order
    - per-query-bin boxplot summary over reference-hit coordinates:
      - `boxplot_bin_count`
      - `boxplot_bins[]` with `min/q1/median/q3/max` and `hit_count`
    - optional aggregated density tiles (multiresolution)
- Flexibility payload schema:
  - `gentle.flexibility_track.v1`
  - typed fields:
    - `seq_id`, span
    - bin size
    - one or more score series (for example bendability, A/T-run burden,
      duplex-destabilization proxy)
    - normalization metadata and min/max ranges
- Operations:
  - implemented:
    - `ComputeDotplot { seq_id, reference_seq_id?, span_start_0based?, span_end_0based?, reference_span_start_0based?, reference_span_end_0based?, mode, word_size, step_bp, max_mismatches, tile_bp?, store_as? }`
    - `ComputeFlexibilityTrack { seq_id, span_start_0based?, span_end_0based?, model, bin_bp, smoothing_bp?, store_as? }`
    - `RenderDotplotSvg { seq_id, dotplot_id, path, flex_track_id?, display_density_threshold?, display_intensity_gain? }`
  - optional later:
    - explicit inspect operations (`InspectDotplot`, `InspectFlexibilityTrack`)
      if adapter-neutral retrieval routes need to be operation-based instead of
      shell read commands.

Low-latency strategy (engine-first):

- Use indexed seeds (k-mer/minimizer style) for candidate match discovery, not
  O(n^2) brute-force raster loops.
- Implemented baseline optimization:
  - when `ComputeDotplot.max_mismatches == 0`, engine uses indexed exact-seed
    matching over sampled windows instead of brute-force pairwise comparisons.
- Build multiresolution tiles so initial view is coarse+fast; refine only for
  visible region/zoom level.
- Cache results in project metadata keyed by `(seq_id, span, params_hash)`.
- Reuse cache across GUI redraws and SVG export; no adapter-side recomputation.
- Route long runs through background jobs with cooperative cancellation and
  deterministic progress phases (`index`, `seed-match`, `tile-aggregate`,
  `finalize`).

GUI contract (implemented baseline + follow-ups):

- Implemented:
  - sequence-window `Dotplot map` compact launcher (separate from map background)
  - dedicated standalone `Dotplot` workspace window (full controls + rendering)
  - bounded-span compute controls (compact launcher + workspace)
  - pair-mode helper to fit `ref_start/ref_end` to loaded hit envelope (+padding)
  - linked hover/locked crosshair behavior
  - pair-mode rendering with separate query/reference axis spans
  - query-axis selection sync in pair mode
  - sparse pairwise diagnostics (orientation hinting, strict-parameter warning,
    reference-edge warning)
- Follow-ups:
  - optional overlay mode with low alpha while keeping dedicated panel mode as
    primary for readability
  - compact parameter presets:
  - `fast-preview` (large word, coarse step)
  - `balanced`
  - `sensitive` (small word, fine step)
  - promoter-flexibility toggles as separate tracks shown beneath/alongside
    dotplot (not fused into dotplot colors by default).

Determinism and parity rules:

- Same parameter set must yield byte-stable JSON payloads (ordering and
  rounding rules fixed in engine).
- SVG export consumes only stored/returned engine payloads.
- CLI/JS/Lua interfaces call the same operations and can render/export the same
  artifacts without GUI-only logic.

Testing contract:

- Deterministic fixture tests with known direct repeats, inverted repeats, and
  low-complexity promoter motifs.
- Parameter-sensitivity tests (`word_size`, `max_mismatches`, `step_bp`) with
  stable expected deltas.
- Snapshot tests for `RenderDotplotSvg`.
- Regression tests for cache-key correctness and cancellation behavior.

### Amino-acid translation row contract (deferred)

The sequence-panel amino-acid row path is intentionally deferred and must remain
non-operative until translation semantics are engine-defined.

Deferred-scope rules:

- Frontend renderers must not infer amino-acid rows from raw DNA/exon context.
- Translation-table resolution must be explicit and deterministic (for example
  codon table identifiers, including contexts where table choice differs by
  species/organelle).
- Deterministic amino-acid output must be transcript/CDS-context aware (frame
  and phase semantics are required); exon-only translation requests are
  therefore not implicitly rendered.
- Until the contract is implemented in engine/view-model form, dormant AA-row
  renderer paths must no-op safely (never panic).

## 9. Decision log (concise)

- Adopt single shared engine for all frontends: accepted
- Use JSON operation/workflow protocol for automation: accepted
- Add capabilities endpoint for integration negotiation: accepted
- Keep rendering separate from core operations: accepted
- Promote pool export to engine operation (`ExportPool`): accepted
- Promote sequence/map SVG and lineage SVG export to engine operations
  (`RenderSequenceSvg`, `RenderLineageSvg`): accepted and implemented
- Promote DNA ladder catalog export to engine operation (`ExportDnaLadders`)
  and expose ladder inspection/export across CLI/JS/Lua/shared shell:
  accepted and implemented
- Add opt-in screenshot artifact bridge (`--allow-screenshots`) for
  documentation/progress image generation with active-window-only capture:
  accepted and implemented historically (adapter-level command path; macOS
  backend), then temporarily disabled in this branch due endpoint-security
  warnings until explicit exception/approval is in place
- Postpone auto-updated documentation with graphics:
  accepted and postponed (manual documentation updates remain in effect)
- Add protocol-grade process export as a strategic requirement:
  accepted and planned (engine-level contract pending)
- Promote container semantics to first-class engine state: accepted and
  implemented (`ProjectState.container_state`)
- Add container-first operations: accepted and implemented
  (`DigestContainer`, `MergeContainersById`, `LigationContainer`,
  `FilterContainerByMolecularWeight`)
- Add normalized engine state summary and expose via CLI/JS/Lua:
  accepted and implemented
- Add thin Python adapter wrapper over deterministic CLI contracts
  (`integrations/python/gentle_py`):
  accepted and implemented
- Add TFBS safety cap semantics (`default=500`, `0=unlimited`):
  accepted and implemented
- Add shared TFBS display filtering criteria and enforce parity in SVG export:
  accepted and implemented
- Add shared VCF display filtering criteria (class/FILTER/QUAL/INFO) and enforce
  parity in SVG export: accepted and implemented
- Add dedicated GUI BLAST viewport with background execution/progress and
  project-sequence/pool query modes: accepted and implemented
- Add BigWig genome-signal import operation (`ImportGenomeBigWigTrack`) with
  `bigWigToBedGraph` conversion and shared shell/adapter exposure:
  accepted and implemented
- Defer sequence-panel amino-acid row rendering until transcript/CDS-aware
  translation contracts (explicit table + phase/frame context) are engine-owned:
  accepted and deferred
- Add cancellable background genome-track import in GUI using shared engine
  progress callbacks: accepted and implemented
- Add shared shell `set-param` and JS/Lua helpers (`set_parameter`,
  `set_vcf_display_filter`) for display-parameter parity: accepted and
  implemented
- Add genome-anchor extension operation (`ExtendGenomeAnchor`) and expose it via
  shared shell/CLI/JS/Lua adapters: accepted and implemented
- Add candidate feature-distance geometry controls (`feature_geometry_mode`,
  `feature_boundary_mode`) across engine and shell/CLI: accepted and
  implemented
- Add candidate feature strand-relation control (`feature_strand_relation`) in
  engine operation schema and shell/CLI/GUI adapters: accepted and implemented
- Add tracked genome-signal subscriptions with auto-sync for newly anchored
  sequences via engine-managed subscription metadata:
  accepted and implemented
- Add agent-assistant bridge with catalog-driven transports, per-reply
  execution intents, shared-shell execution, and standalone GUI viewport:
  accepted and implemented
- Keep explicit compatibility target for OpenAI-compatible hosted providers
  (including [clawbio.ai](https://clawbio.ai/)) via `native_openai_compat`:
  accepted and implemented
- Add ClawBio/OpenClaw skill scaffold that wraps deterministic `gentle_cli`
  routes and emits reproducibility bundles:
  accepted and implemented (`integrations/clawbio/skills/gentle-cloning`)
- Add curated isoform-panel import + TP53-style transcript/protein expert view
  + deterministic SVG export route (`ImportIsoformPanel`,
  `RenderIsoformArchitectureSvg`): accepted and implemented
- Add async long-running command execution contract for agent-suggested BLAST
  and primer-pair multi-BLAST workflows: accepted and planned
- Add shared GUI intent command plane (`ui intents`, `ui open|focus`, prepared
  query helpers) and wire it into GUI-host dialog openers with deterministic
  prepared-reference selection: accepted and implemented (baseline)
- Add MCP server adapter that exposes shared deterministic command/operation
  surfaces to external AI tools without duplicating biology logic:
  accepted and implemented as a parity-expanded baseline
  (`capabilities`, `state_summary`, `op`, `workflow`, `help`,
  `ui_intents`, `ui_intent`, `ui_prepared_genomes`, `ui_latest_prepared`)
  with explicit confirmation semantics for mutating tools and routed shared
  parser/executor paths for UI-intent tools
- Replace runtime process-environment mutation for tool-path overrides with a
  process-local override registry (`tool_overrides`) to keep Rust 2024-safe
  behavior without unsafe env writes: accepted and implemented
- Prevent idle redraw churn by deduplicating window-title viewport commands in
  GUI update loop: accepted and implemented
- Defer for now: `import-pool` engine operation until first-class container
  semantics settle for stable import behavior across adapters.
- Add guideRNA-design base layer:
  accepted and partially implemented as first-class engine operations
  (`UpsertGuideSet`, `DeleteGuideSet`, `FilterGuidesPractical`,
  `GenerateGuideOligos`, `ExportGuideOligos`, `ExportGuideProtocolText`) with
  shared-shell/CLI exposure via `guides ...`; remaining phase covers off-target
  ranking and macro-template packaging (`docs/rna_guides_spec.md`).
- Add design-constraint filtering operation (`FilterByDesignConstraints`;
  GC bounds, homopolymer cap, U6 `TTTT` avoidance, forbidden motifs): accepted
  and implemented.
- Add feature expert-view model + SVG export contract (TFBS + restriction
  sites) and expose through GUI/shell/CLI/JS/Lua:
  accepted and implemented.
- Add first-class alternative-splicing interpretation view model (transcript
  lanes, splice boundaries, event summaries, junction overlays) with strict
  GUI/SVG parity and no label-driven geometry distortion:
  accepted and implemented baseline (follow-up hardening remains: dense-fixture
  regressions + dedicated primary map-mode integration).
- Add candidate-set scoring/filter/set-algebra workflow (derived expressions,
  value+quantile filters, and explicit set union/intersection/subtraction):
  accepted and implemented as first-class engine operations, exposed through
  first-class CLI (`gentle_cli candidates ...`), shared shell (`candidates`),
  and dedicated GUI Engine Ops candidate forms.
- Add first-class primer-pair design/report baseline (`DesignPrimerPairs`) with
  shared-shell report inspection/export (`primers design|list-reports|show-report|export-report`)
  and metadata persistence (`gentle.primer_design_report.v1`):
  accepted and implemented baseline (follow-up planned for richer thermodynamic
  scoring and off-target evaluation tiers).
- Add Primer3 wrapper backend integration for `DesignPrimerPairs` with
  deterministic request/response normalization, provenance capture, and
  adapter-equivalent behavior across GUI/CLI/JS/Lua/agent/MCP:
  accepted and planned.
- Add optional per-window-type visual identity backdrops with persisted app
  settings and subtle tint/watermark rendering:
  accepted and in progress (experimental GUI implementation).
- Unify configuration editing contract with explicit staged-change signaling and
  always-visible commit footer controls (`Cancel`/`Apply`) for configuration
  viewports: accepted and implemented.
- Adopt contract-first source documentation policy to support
  human/Codex maintainability without high-noise churn:
  accepted and in progress.
  - Add concise rustdoc comments (`//!`, `///`) on core modules/types/functions
    that define contracts, invariants, and non-obvious semantics.
  - Prefer "why/invariant/edge-case" documentation over line-by-line narration.
  - Keep documentation updates in focused low-noise passes so delivery and
    rebases stay stable.
