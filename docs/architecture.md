# GENtle Architecture (Working Draft)

Last updated: 2026-02-18

This document describes how GENtle is intended to work, what is already
implemented, and what should be built next.

It has two goals:

- Human-readable architecture guide
- Reliable handoff/resume note if work is interrupted

## 1. Product intent

GENtle is a DNA/cloning workbench with multiple access paths:

- GUI for interactive use
- JavaScript shell for scripted experiments
- Lua shell for scripted experiments
- CLI for automation and AI tools

The long-term requirement is strict behavioral parity:

- The same biological operations must run through the same core routines,
  regardless of entry point.

Wet-lab semantic rule (target model):

- A main DNA window represents a wet-lab container (tube/vial), not a single
  guaranteed molecule.
- A container may hold multiple candidate molecules/fragments.
- Filter-like steps (PCR, gel extraction, in silico selection) produce a new
  container with a narrower candidate set.

## 2. Core architecture rule

GENtle should follow a single-engine architecture:

1. `Core engine` (single source of truth)
2. `Adapters/frontends` (GUI, JS, Lua, CLI)
3. `View model/renderers` (presentation layer only)

### Non-negotiable invariant

Frontends must not implement their own cloning/business logic.
They only translate user input into engine operations and display results.

## 3. Current implementation status

### Already in place

- New shared engine module: `src/engine.rs`
- New automation CLI: `src/bin/gentle_cli.rs`
- Protocol draft doc: `docs/protocol.md`
- CLI docs updated: `docs/cli.md`
- Project lineage/provenance graph in engine state (`ProjectState.lineage`)
- Main window lineage view (branch/merge aware sequence history list)
- First-class container state in engine (`ProjectState.container_state`)
- Container-aware operations in engine (`DigestContainer`,
  `MergeContainersById`, `LigationContainer`,
  `FilterContainerByMolecularWeight`)
- Shared normalized state summary from engine (`EngineStateSummary`)
- Resource sync/auto-feed baseline:
  - CLI import of REBASE (`resources sync-rebase`) and JASPAR (`resources sync-jaspar`)
  - JS shell resource sync (`sync_rebase`, `sync_jaspar`)
  - Lua shell resource sync (`sync_rebase`, `sync_jaspar`)
  - GUI File menu import actions (`Import REBASE Data...`, `Import JASPAR Data...`)
  - Runtime loading of local REBASE enzyme snapshot and local JASPAR motif registry
- Built-in JASPAR motif snapshot (2026 CORE non-redundant derived) is shipped in `assets/jaspar.motifs.json`
    with runtime override at `data/resources/jaspar.motifs.json`
- TFBS runtime guardrails:
  - `AnnotateTfbs.max_hits` supports safe caps
  - default cap is `500` accepted hits per operation
  - `max_hits: 0` means unlimited (explicit opt-out for CLI/automation)
- TFBS progress + responsiveness hardening:
  - engine emits coarse-grained progress for long scans
  - GUI worker forwarding and polling are throttled to avoid UI starvation
- TFBS display filtering is first-class and persistent:
  - four optional criteria (`llr_bits`, `llr_quantile`,
    `true_log_odds_bits`, `true_log_odds_quantile`)
  - filter state is synchronized into `ProjectState.display`
- SVG export parity for TFBS filtering:
  - linear/circular SVG export uses the same TFBS display criteria as GUI
  - no recomputation required to reduce visual TFBS density

### Engine capabilities currently implemented

- State container (`ProjectState`)
- Operations:
  - `LoadFile`
  - `SaveFile`
  - `RenderSequenceSvg`
  - `RenderLineageSvg`
  - `ExportPool`
  - `DigestContainer`
  - `MergeContainersById`
  - `LigationContainer`
  - `FilterContainerByMolecularWeight`
  - `Digest`
  - `MergeContainers`
  - `Ligation` (protocol-driven: sticky/blunt)
  - `Pcr`
  - `PcrAdvanced`
  - `PcrMutagenesis`
  - `ExtractRegion`
  - `ExtractAnchoredRegion`
  - `SelectCandidate`
  - `FilterByMolecularWeight`
  - `Reverse`
  - `Complement`
  - `ReverseComplement`
  - `Branch`
  - `SetDisplayVisibility`
  - `SetTopology`
  - `RecomputeFeatures`
  - `SetParameter`
  - `AnnotateTfbs` (TFBS annotation with LLR + true-log-odds scoring, progress reporting,
    and safety cap support)
- Workflow execution (`Workflow`)
- Structured errors (`ErrorCode` + message)
- Capabilities negotiation (`GentleEngine::capabilities()`)
- Shared state summary (`GentleEngine::summarize_state()`)
- Persistent state through CLI JSON file

### GUI status relevant to architecture

- Linear renderer rewritten
- Added overlap management and strand-aware placement in linear mode
- Added button tooltips
- Added linear overlays for ORF/GC/methylation/ticks
- Main lineage view supports both table/graph and a container list with open
  actions (sequence/pool)
- Engine Ops panel supports asynchronous TFBS annotation with live per-motif and
  total progress bars
- TFBS display filtering is interactive (checkbox criteria + thresholds) and
  persistent per sequence
- Engine Ops panel area is scrollable and vertically resizable

### Status matrix

Legend:

- `Done`: implemented and available
- `Partial`: available but not yet routed through shared engine everywhere
- `Planned`: not implemented yet

| Area | GUI | JS | Lua | CLI (`gentle_cli`) | Core engine |
|---|---|---|---|---|---|
| Load sequence file | Done | Done | Done | Done | Done |
| Save sequence file | Done | Done | Done | Done | Done |
| Restriction digest | Done | Done | Done | Done | Done |
| Container state model | Done | Done | Done | Done | Done |
| Container-first operations | Partial | Done | Done | Done | Done |
| Topology change (linear/circular) | Done | Done | Done | Done | Done |
| Feature recomputation | Partial | Done | Done | Done | Done |
| Ligation | Done | Done | Done | Done | Done |
| PCR | Done | Done | Done | Done | Done |
| Region extraction/editing | Partial | Done | Done | Done | Done |
| View toggles (show/hide tracks and panels) | Done | Done | Done | Done | Done |
| TFBS annotation with progress | Done | Done | Done | Done | Done |
| TFBS display filtering (4 criteria) | Done | Done | Done | Done | Done |
| Sequence SVG export | Partial | Done | Done | Done | Done |
| Lineage SVG export | Partial | Done | Done | Done | Done |
| Pool export (overhang-aware) | Missing | Done | Done | Done | Done |
| State summary (seq + container) | Done | Done | Done | Done | Done |
| Shared operation protocol | Partial | Done | Done | Done | Done |

### Operation parity matrix (code-backed, current)

This matrix is stricter than the high-level status matrix above.
It describes whether a frontend currently executes the corresponding shared
engine operation.

Legend:

- `Wired`: frontend invokes `Operation` through `GentleEngine`
- `Exposed`: available through frontend API/command, but not necessarily bound in GUI
- `Missing`: not currently available in that adapter

| Operation | GUI | CLI (`gentle_cli`) | JS shell | Lua shell | Engine |
|---|---|---|---|---|---|
| `LoadFile` | Wired | Wired | Exposed | Exposed | Implemented |
| `SaveFile` | Wired | Wired | Exposed | Exposed | Implemented |
| `RenderSequenceSvg` | Missing | Wired | Exposed | Exposed | Implemented |
| `RenderLineageSvg` | Missing | Wired | Exposed | Exposed | Implemented |
| `ExportPool` | Missing | Wired | Exposed | Exposed | Implemented |
| `DigestContainer` | Missing | Wired | Exposed | Exposed | Implemented |
| `MergeContainersById` | Missing | Wired | Exposed | Exposed | Implemented |
| `LigationContainer` | Missing | Wired | Exposed | Exposed | Implemented |
| `FilterContainerByMolecularWeight` | Missing | Wired | Exposed | Exposed | Implemented |
| `Digest` | Wired | Wired | Exposed | Exposed | Implemented |
| `MergeContainers` | Wired | Wired | Exposed | Exposed | Implemented |
| `Ligation` | Wired | Wired | Exposed | Exposed | Implemented |
| `Pcr` | Wired | Wired | Exposed | Exposed | Implemented |
| `PcrAdvanced` | Wired | Wired | Exposed | Exposed | Implemented |
| `PcrMutagenesis` | Wired | Wired | Exposed | Exposed | Implemented |
| `ExtractRegion` | Missing | Wired | Exposed | Exposed | Implemented |
| `ExtractAnchoredRegion` | Wired | Wired | Exposed | Exposed | Implemented |
| `SelectCandidate` | Wired | Wired | Exposed | Exposed | Implemented |
| `FilterByMolecularWeight` | Wired | Wired | Exposed | Exposed | Implemented |
| `Reverse` | Wired | Wired | Exposed | Exposed | Implemented |
| `Complement` | Wired | Wired | Exposed | Exposed | Implemented |
| `ReverseComplement` | Wired | Wired | Exposed | Exposed | Implemented |
| `Branch` | Wired | Wired | Exposed | Exposed | Implemented |
| `SetDisplayVisibility` | Wired | Wired | Exposed | Exposed | Implemented |
| `SetTopology` | Wired | Wired | Exposed | Exposed | Implemented |
| `RecomputeFeatures` | Missing | Wired | Exposed | Exposed | Implemented |
| `SetParameter` | Missing | Wired | Exposed | Exposed | Implemented |
| `AnnotateTfbs` | Wired | Wired | Exposed | Exposed | Implemented |
| `Workflow` (multi-op run) | Missing | Wired | Exposed | Exposed | Implemented |

Notes from current code:

- GUI now routes most day-to-day cloning actions through engine operations
  (digest/merge/ligation/PCR variants/reverse-family/toggles/save).
- GUI TFBS feature density can be reduced at display time via score-based
  criteria without mutating sequence annotations.
- SVG export now honors TFBS display criteria from shared display state.
- GUI has a container-centric read/view layer, but container-first operations
  are not yet directly exposed as GUI actions.
- Remaining GUI operation gaps are mainly `ExtractRegion`, `RecomputeFeatures`,
  `SetParameter`, container-first operation actions, `RenderSequenceSvg`,
  `RenderLineageSvg`, and generic workflow execution.
- CLI exposes all implemented operations (`op`/`workflow`) and adds some
  adapter-level utilities (render and import helpers).
- REBASE/JASPAR snapshot update paths are now exposed across CLI, JS, Lua, and
  GUI (GUI through File-menu import actions).
- JS/Lua expose generic operation/workflow bridges through `GentleEngine`
  (`apply_operation`, `apply_workflow`) and therefore can reach the full
  operation set.
- CLI/JS/Lua now share engine-provided normalized state summaries
  (`summarize_state` / `state_summary`).

Immediate parity goal:

1. Finish remaining GUI operation gaps (`ExtractRegion`, `RecomputeFeatures`,
   `SetParameter`, container-first actions, optional generic workflow runner).
2. Move remaining CLI-only utilities that represent stable contracts into engine-level
   operations.
3. Keep JS/Lua convenience wrappers thin over shared operation contracts.

## 4. Engine model

### State

`ProjectState` stores:

- `sequences: HashMap<SeqId, DNAsequence>`
- `metadata: HashMap<String, serde_json::Value>`
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

Display state (`ProjectState.display`) now also carries TFBS display-filter
criteria so interactive GUI filtering and SVG export use one shared contract.

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

### CLI (`gentle_cli`)

- Stable JSON interface for scripted and AI-driven workflows
- State import/export and summary
- Capabilities reporting

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
mainly GUI wiring gaps described in section 11.

## 8. Rendering and interpretation direction

Two additional layers are needed for full human/AI collaboration:

1. `state -> view model -> rendering`
2. `human drawing/image -> interpreted state proposal` (later)

Important separation:

- Engine remains deterministic and biology-centric
- Rendering remains presentation-centric
- Image interpretation remains probabilistic and should output proposals,
  not direct mutations without confirmation

## 9. Roadmap (recommended order)

### Phase A: expand operation parity

Phase A core scope is largely complete.

Remaining Phase A items:

- close GUI operation gaps (`ExtractRegion`, `RecomputeFeatures`, `SetParameter`)
- expose container-first operations directly in GUI action surfaces
- wire GUI export actions to shared render operations (`RenderSequenceSvg`,
  `RenderLineageSvg`)
- keep adapter-level helpers as wrappers over engine operations

### Phase B: GUI migration

Route GUI actions through engine operations.
Keep UI state thin and avoid direct mutations where possible.

### Phase C: view model contract

Introduce a frontend-neutral view model schema:

- tracks
- features
- labels
- overlays (ORF/GC/methylation)
- selections

Use same view model for GUI and machine consumers.

### Phase D: protocol hardening

- versioned schema files
- compatibility policy
- richer error taxonomy and validation
- operation provenance metadata (engine version, timestamps, input references)

### Phase E: interpretation (later)

- parse sketches/images into `StatePatchProposal`
- confidence scores
- explicit user/agent confirmation before apply

## 10. Practical resume checklist

If work is interrupted, resume in this order:

1. Read this file (`docs/architecture.md`)
2. Read machine contract (`docs/protocol.md`)
3. Inspect engine (`src/engine.rs`)
4. Inspect CLI adapter (`src/bin/gentle_cli.rs`)
5. Run quick sanity:
   - `cargo check -q`
   - `cargo run --bin gentle_cli -- capabilities`
   - `cargo run --bin gentle_cli -- state-summary`
   - `cargo test -q render_export::tests::`
6. Continue with next Phase A operation implementation

## 11. Current known gaps

- GUI is close to parity but still missing a few engine operations
- GUI does not yet provide direct triggers for container-first operations
  (`DigestContainer`, `MergeContainersById`, `LigationContainer`,
  `FilterContainerByMolecularWeight`)
- GUI exports are not yet fully routed through engine render operations
  (`RenderSequenceSvg`, `RenderLineageSvg`) despite those operations being
  implemented and available in CLI/JS/Lua.
- View model contract is not yet formalized
- Some rendering/import/export utilities are still adapter-level contracts
  instead of engine operations

## 12. Decision log (concise)

- Adopt single shared engine for all frontends: accepted
- Use JSON operation/workflow protocol for automation: accepted
- Add capabilities endpoint for integration negotiation: accepted
- Keep rendering separate from core operations: accepted
- Promote pool export to engine operation (`ExportPool`): accepted
- Promote sequence/map SVG and lineage SVG export to engine operations
  (`RenderSequenceSvg`, `RenderLineageSvg`): accepted and implemented
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
- Defer for now: `import-pool` engine operation until first-class container
  semantics settle for stable import behavior across adapters.
