# GENtle Architecture (Working Draft)

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

### Engine capabilities currently implemented

- State container (`ProjectState`)
- Operations:
  - `LoadFile`
  - `SaveFile`
  - `Digest`
  - `SetTopology`
  - `RecomputeFeatures`
- Workflow execution (`Workflow`)
- Structured errors (`ErrorCode` + message)
- Capabilities negotiation (`GentleEngine::capabilities()`)
- Persistent state through CLI JSON file

### GUI status relevant to architecture

- Linear renderer rewritten
- Added overlap management and strand-aware placement in linear mode
- Added button tooltips
- Added linear overlays for ORF/GC/methylation/ticks

### Status matrix

Legend:

- `Done`: implemented and available
- `Partial`: available but not yet routed through shared engine everywhere
- `Planned`: not implemented yet

| Area | GUI | JS | Lua | CLI (`gentle_cli`) | Core engine |
|---|---|---|---|---|---|
| Load sequence file | Done | Done | Done | Done | Done |
| Save sequence file | Done | Done | Done | Done | Done |
| Restriction digest | Done | Done | Planned | Done | Done |
| Topology change (linear/circular) | Done | Partial | Partial | Done | Done |
| Feature recomputation | Done | Partial | Partial | Done | Done |
| Ligation | Planned | Planned | Planned | Done | Done |
| PCR | Planned | Planned | Planned | Done | Done |
| Region extraction/editing | Partial | Partial | Partial | Done | Done |
| View toggles (show/hide tracks and panels) | Done | Planned | Planned | Done | Done |
| Shared operation protocol | Partial | Partial | Partial | Done | Done |

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
| `SaveFile` | Missing | Wired | Exposed | Exposed | Implemented |
| `Digest` | Missing | Wired | Exposed | Exposed | Implemented |
| `MergeContainers` | Missing | Wired | Exposed | Exposed | Implemented |
| `Ligation` | Missing | Wired | Exposed | Exposed | Implemented |
| `Pcr` | Missing | Wired | Exposed | Exposed | Implemented |
| `PcrAdvanced` | Missing | Wired | Exposed | Exposed | Implemented |
| `PcrMutagenesis` | Missing | Wired | Exposed | Exposed | Implemented |
| `ExtractRegion` | Missing | Wired | Exposed | Exposed | Implemented |
| `SelectCandidate` | Missing | Wired | Exposed | Exposed | Implemented |
| `FilterByMolecularWeight` | Missing | Wired | Exposed | Exposed | Implemented |
| `Reverse` | Wired | Wired | Exposed | Exposed | Implemented |
| `Complement` | Wired | Wired | Exposed | Exposed | Implemented |
| `ReverseComplement` | Wired | Wired | Exposed | Exposed | Implemented |
| `Branch` | Wired | Wired | Exposed | Exposed | Implemented |
| `SetDisplayVisibility` | Wired | Wired | Exposed | Exposed | Implemented |
| `SetTopology` | Wired | Wired | Exposed | Exposed | Implemented |
| `RecomputeFeatures` | Missing | Wired | Exposed | Exposed | Implemented |
| `SetParameter` | Missing | Wired | Exposed | Exposed | Implemented |
| `Workflow` (multi-op run) | Missing | Wired | Exposed | Exposed | Implemented |

Notes from current code:

- GUI menu currently exposes file-open and DNA display toggles/topology through
  engine operations.
- GUI does not yet expose save/digest/ligation/PCR/extract/workflow actions as
  engine operations.
- CLI exposes all implemented operations and workflow execution.
- JS/Lua now expose generic operation/workflow bridges through `GentleEngine`
  (`apply_operation`, `apply_workflow`), with convenience helpers still present.

Immediate parity goal for Step 1 completion:

1. Enumerate each GUI action and bind it to exactly one `Operation` or
   `Workflow`.
2. Add missing GUI actions for core operations already available in engine/CLI.
3. Add JS/Lua wrappers that call engine operations instead of bespoke routines.

## 4. Engine model

### State

`ProjectState` stores:

- `sequences: HashMap<SeqId, DNAsequence>`
- `metadata: HashMap<String, serde_json::Value>`
- `lineage: LineageGraph`
  - `nodes` (sequence lineage nodes with origin + creation op)
  - `seq_to_node` (current sequence id -> latest lineage node)
  - `edges` (parent -> child derivation edges, including multi-parent ligation)

Planned extension for wet-lab semantics:

- Introduce explicit container/pool entities so operations can target
  "all candidates in container X" instead of one selected sequence id.
- Keep `SelectCandidate` as explicit in-silico disambiguation when multiple
  possible products exist.

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

Current work satisfies (1) through (4) at draft level for a subset of operations.

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

Add operations to engine and CLI:

- `MergeContainers` (or `Mix`)
- protocol-driven `Ligation` over merged candidate pools
- `PCR`
- region extraction/edit operations
- annotation operations

Then migrate JS/Lua wrappers to call engine instead of bespoke logic.

Ligation protocol rule (planned):

- `Ligation` must accept a `protocol` argument.
- End compatibility behavior is derived from protocol, not free-form booleans.
- Built-in minimal protocols:
  - `sticky`
  - `blunt`
- Additional established protocols can be added as named presets later.

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
6. Continue with next Phase A operation implementation

## 11. Current known gaps

- GUI is not yet fully migrated to the shared operation API surface
- Container/pool semantics are not yet first-class in engine state
- Ligation is not yet protocol-driven and sticky-end compatibility aware
- View model contract is not yet formalized
- Operation set is still a subset of required cloning workflows

## 12. Decision log (concise)

- Adopt single shared engine for all frontends: accepted
- Use JSON operation/workflow protocol for automation: accepted
- Add capabilities endpoint for integration negotiation: accepted
- Keep rendering separate from core operations: accepted
