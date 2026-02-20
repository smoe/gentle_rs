# GENtle Architecture (Working Draft)

Last updated: 2026-02-20

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
- CLI for automation and AI tools

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

1. Keep one deterministic engine contract across GUI, CLI, JS, and Lua.
2. Preserve provenance so every derived result can be traced and replayed.
3. Make every process exportable as a human-readable protocol text that a
   technical assistant can follow step by step (inputs, operations, expected
   outputs, and checkpoints).
4. Support optional, explicit screenshot artifact generation for documentation
   and progress communication without weakening default safety boundaries
   (compile-time feature gate + explicit runtime opt-in).
5. Provide context-sensitive hover descriptions for actionable GUI buttons so
   control intent is explicit and shared button references are unambiguous.

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

## 2. Core architecture rule

GENtle should follow a single-engine architecture:

1. `Core engine` (single source of truth)
2. `Adapters/frontends` (GUI, JS, Lua, CLI)
3. `View model/renderers` (presentation layer only)

### Non-negotiable invariant

Frontends must not implement their own cloning/business logic.
They only translate user input into engine operations and display results.

## 3. Implementation status and roadmap

Shared implementation status, parity matrices, known gaps, and execution order
are now maintained in `docs/roadmap.md`.

## 4. Engine model

### State

`ProjectState` stores:

- `sequences: HashMap<SeqId, DNAsequence>`
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
- eventual undo/replay support

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

Practical rule:

- Context fields are optional, but lineage fields (`state_id`, `parents`,
  operation edge) are mandatory for DAG integrity.

## 5. Adapter responsibilities

### GUI

- Convert UI events to engine operations
- Render state and selections
- Never duplicate biology logic

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

Current work satisfies (1) through (4) for most operations; remaining gaps are
mainly view-model formalization and promoting remaining adapter-level utility
contracts into stable engine operations.

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
- `candidates score SET_NAME METRIC_NAME EXPRESSION`
- `candidates score-distance SET_NAME METRIC_NAME ...`
- `candidates filter INPUT_SET OUTPUT_SET --metric METRIC ...`
- `candidates set-op union|intersect|subtract LEFT RIGHT OUTPUT`
- `candidates macro [--transactional] [--file PATH | SCRIPT_OR_@FILE]`
- `set-param NAME JSON_VALUE`

This enables reusable query composition:

1. generate explicit candidate sets
2. attach base/derived scores
3. filter by absolute threshold and/or quantile
4. intersect/union/subtract with other candidate sets

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
- Prevent idle redraw churn by deduplicating window-title viewport commands in
  GUI update loop: accepted and implemented
- Defer for now: `import-pool` engine operation until first-class container
  semantics settle for stable import behavior across adapters.
- Add guideRNA-design base layer (practical design-constraint filters + oligo/export
  + macro-template workflow model): accepted and planned (spec draft in
  `docs/rna_guides_spec.md`).
- Add design-constraint filtering operation (`FilterByDesignConstraints`;
  GC bounds, homopolymer cap, U6 `TTTT` avoidance, forbidden motifs): accepted
  and implemented.
- Add candidate-set scoring/filter/set-algebra workflow (derived expressions,
  value+quantile filters, and explicit set union/intersection/subtraction):
  accepted and implemented as first-class engine operations, exposed through
  first-class CLI (`gentle_cli candidates ...`), shared shell (`candidates`),
  and dedicated GUI Engine Ops candidate forms.
