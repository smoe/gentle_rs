# GENtle Architecture (Working Draft)

Last updated: 2026-02-19

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
   and progress communication without weakening default safety boundaries.

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
- New shared shell command layer: `src/engine_shell.rs`
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
- Reference genome prepare/extract baseline:
  - engine operation `PrepareGenome` (one-time local install/cache of sequence + annotation)
  - engine operation `ExtractGenomeRegion` (indexed region retrieval into project sequence state)
  - engine operation `ExtractGenomeGene` (gene-name/id-based retrieval from indexed annotations)
  - optional `cache_dir` override passed through both operations
  - progress events for background genome preparation in GUI
  - annotation-derived gene listing API for GUI gene-name filtering and gene extraction
  - persistent `genes.json` index built during prepare for fast retrieval
  - local manifest + FASTA index persisted per prepared genome
  - catalog support for NCBI assembly-derived downloads via
    `ncbi_assembly_accession` + `ncbi_assembly_name` (direct GenBank/RefSeq FTP)
  - catalog support for GenBank accession-derived helper downloads via
    `genbank_accession` (NCBI EFetch FASTA + GenBank annotation)
  - GenBank FEATURE normalization into indexed retrieval records (for vector
    elements such as promoter/CDS/origin/LTR/ITR in helper workflows)
  - additional starter helper-system catalog shipped as `assets/helper_genomes.json`
    for local vector inventories (plasmid/lenti/adeno/AAV + yeast/E. coli hosts)
  - catalog validation hardening for inconsistent source declarations
    (assembly vs GenBank fields, unpublished placeholders, missing source fields)
  - HTTP download hardening for genome preparation (timeouts, retry/backoff on
    transient status/network errors, clearer fetch error reporting)
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
- Linear map viewport state is first-class and persistent:
  - zoom/pan state (`linear_view_start_bp`, `linear_view_span_bp`)
  - display-density LOD hides cluttered tracks at coarse scales (notably restriction sites)
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
  - `RenderPoolGelSvg`
  - `ExportDnaLadders`
  - `ExportRnaLadders`
  - `ExportPool`
  - `PrepareGenome`
  - `ExtractGenomeRegion`
  - `ExtractGenomeGene`
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
  - `SetLinearViewport`
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
- Circular renderer now supports multipart feature locations (`Join`, `Order`,
  nested `Complement`/`External`): exon-like segments are rendered separately
  with intron arches between them
- Circular intron arches have feature-aware styling (notably emphasized for
  mRNA) with hover/selection accenting for readability
- Added button tooltips
- Added linear overlays for ORF/GC/methylation/ticks
- Main lineage view supports both table/graph and a container list with open
  actions (sequence/pool)
- Pool-context Engine Ops includes ladder-aware virtual gel preview and shared
  SVG export route
- GUI DNA window now includes `Shell` panel backed by shared shell command
  parsing/execution (same command semantics as `gentle_cli shell`)
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
| Container-first operations | Done | Done | Done | Done | Done |
| Topology change (linear/circular) | Done | Done | Done | Done | Done |
| Feature recomputation | Done | Done | Done | Done | Done |
| Ligation | Done | Done | Done | Done | Done |
| PCR | Done | Done | Done | Done | Done |
| Region extraction/editing | Done | Done | Done | Done | Done |
| View toggles (show/hide tracks and panels) | Done | Done | Done | Done | Done |
| TFBS annotation with progress | Done | Done | Done | Done | Done |
| TFBS display filtering (4 criteria) | Done | Done | Done | Done | Done |
| Sequence SVG export | Done | Done | Done | Done | Done |
| Lineage SVG export | Done | Done | Done | Done | Done |
| Pool gel SVG export (auto ladder selection) | Done | Done | Done | Done | Done |
| Pool export (overhang-aware) | Done | Done | Done | Done | Done |
| Reference genome prepare + gene/region extraction | Done | Done | Done | Done | Done |
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
| `RenderSequenceSvg` | Wired | Wired | Exposed | Exposed | Implemented |
| `RenderLineageSvg` | Wired | Wired | Exposed | Exposed | Implemented |
| `RenderPoolGelSvg` | Wired | Wired | Exposed | Exposed | Implemented |
| `ExportDnaLadders` | Exposed | Wired | Exposed | Exposed | Implemented |
| `ExportRnaLadders` | Exposed | Wired | Exposed | Exposed | Implemented |
| `ExportPool` | Wired | Wired | Exposed | Exposed | Implemented |
| `PrepareGenome` | Wired | Wired | Exposed | Exposed | Implemented |
| `ExtractGenomeRegion` | Wired | Wired | Exposed | Exposed | Implemented |
| `ExtractGenomeGene` | Wired | Wired | Exposed | Exposed | Implemented |
| `DigestContainer` | Wired | Wired | Exposed | Exposed | Implemented |
| `MergeContainersById` | Wired | Wired | Exposed | Exposed | Implemented |
| `LigationContainer` | Wired | Wired | Exposed | Exposed | Implemented |
| `FilterContainerByMolecularWeight` | Wired | Wired | Exposed | Exposed | Implemented |
| `Digest` | Wired | Wired | Exposed | Exposed | Implemented |
| `MergeContainers` | Wired | Wired | Exposed | Exposed | Implemented |
| `Ligation` | Wired | Wired | Exposed | Exposed | Implemented |
| `Pcr` | Wired | Wired | Exposed | Exposed | Implemented |
| `PcrAdvanced` | Wired | Wired | Exposed | Exposed | Implemented |
| `PcrMutagenesis` | Wired | Wired | Exposed | Exposed | Implemented |
| `ExtractRegion` | Wired | Wired | Exposed | Exposed | Implemented |
| `ExtractAnchoredRegion` | Wired | Wired | Exposed | Exposed | Implemented |
| `SelectCandidate` | Wired | Wired | Exposed | Exposed | Implemented |
| `FilterByMolecularWeight` | Wired | Wired | Exposed | Exposed | Implemented |
| `Reverse` | Wired | Wired | Exposed | Exposed | Implemented |
| `Complement` | Wired | Wired | Exposed | Exposed | Implemented |
| `ReverseComplement` | Wired | Wired | Exposed | Exposed | Implemented |
| `Branch` | Wired | Wired | Exposed | Exposed | Implemented |
| `SetDisplayVisibility` | Wired | Wired | Exposed | Exposed | Implemented |
| `SetLinearViewport` | Wired | Wired | Exposed | Exposed | Implemented |
| `SetTopology` | Wired | Wired | Exposed | Exposed | Implemented |
| `RecomputeFeatures` | Wired | Wired | Exposed | Exposed | Implemented |
| `SetParameter` | Wired | Wired | Exposed | Exposed | Implemented |
| `AnnotateTfbs` | Wired | Wired | Exposed | Exposed | Implemented |
| `Workflow` (multi-op run) | Wired | Wired | Exposed | Exposed | Implemented |

Notes from current code:

- GUI now routes most day-to-day cloning actions through engine operations
  (digest/merge/ligation/PCR variants/reverse-family/toggles/save).
- GUI TFBS feature density can be reduced at display time via score-based
  criteria without mutating sequence annotations.
- SVG export now honors TFBS display criteria from shared display state.
- GUI now exposes direct actions for `ExtractRegion`, `RecomputeFeatures`,
  `SetParameter`, container-first operations, sequence SVG rendering, lineage
  SVG rendering, and workflow execution through shared engine operations.
- GUI now exposes dedicated controls for `PrepareGenome` and
  `ExtractGenomeRegion` from the main-window menu as separate dialogs:
  `Prepare Reference Genome...` and `Retrieve Genome Sequence...`.
- GUI now exposes a third reference-genome dialog, `Prepared References...`,
  to inspect prepared installations (paths, readiness flags, source types,
  and checksum fingerprints).
- GUI also exposes helper-catalog shortcuts (`Prepare Helper Genome...`,
  `Retrieve Helper Sequence...`) that preselect helper catalog/cache defaults.
- GUI genome selection is catalog-backed dropdown (no free-text genome id
  entry); prepare dialogs only list unprepared entries and retrieval dialogs
  only list prepared entries. Retrieval provides paged/top-N regex filtering,
  biotype checkbox filtering from parsed annotations, and direct engine-backed
  extraction.
- GUI operation parity is near-complete for sequence/container workflows; some
  utility-style contracts (notably `ExportDnaLadders` / `ExportRnaLadders`) are currently exposed
  through shared shell/adapter surfaces rather than dedicated first-class GUI
  controls.
- CLI exposes all implemented operations (`op`/`workflow`) and adds some
  adapter-level utilities (render and import helpers).
- CLI now exposes a shared shell command path (`gentle_cli shell ...`) that
  reuses `src/engine_shell.rs` (also used by GUI `Shell` panel).
- Shared shell command coverage now includes `genomes`, `helpers`,
  `resources`, `ladders`, and `import-pool`, and `gentle_cli` top-level dispatch routes
  these trees through the same shared parser/executor used by GUI Shell.
- Shared shell/CLI `genomes status` and `helpers status` now include resolved
  source type reporting (`sequence_source_type`, `annotation_source_type`) in
  addition to prepared/not-prepared state.
- Shared shell/CLI now include `genomes validate-catalog` and
  `helpers validate-catalog` for preflight catalog validation.
- Genome preparation now persists SHA-1 integrity fields in per-genome manifest
  files and validates/backfills them on cache reuse; HTTP source downloads use
  resumable Range requests with retry/backoff.
- CLI includes dedicated `helpers` convenience subcommands (list/status/genes/
  prepare/extract-region/extract-gene) that default to
  `assets/helper_genomes.json`.
- CLI/JS/Lua now expose dedicated reference-genome helper surfaces in addition
  to raw operation bridges:
  - catalog listing
  - prepared-status checks
  - gene-index listing
  - prepare/extract convenience wrappers
- REBASE/JASPAR snapshot update paths are now exposed across CLI, JS, Lua, and
  GUI (GUI through File-menu import actions).
- JS/Lua expose generic operation/workflow bridges through `GentleEngine`
  (`apply_operation`, `apply_workflow`) and therefore can reach the full
  operation set.
- CLI/JS/Lua now share engine-provided normalized state summaries
  (`summarize_state` / `state_summary`).

Immediate parity goal:

1. Maintain operation parity as new operations are added.
2. Move remaining CLI-only utilities that represent stable contracts into
   engine-level operations.
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
- Current implementation now appends extraction-level provenance entries at
  `ProjectState.metadata["provenance"]["genome_extractions"]` for
  `ExtractGenomeRegion` and `ExtractGenomeGene`, including source descriptors
  and checksum fields when available.

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
- Shared shell command path (`shell`) aligned with GUI Shell panel
- Optional screenshot bridge guarded by `--allow-screenshots`:
  - default behavior: disabled
  - enabled behavior: allow invoking the host system screenshot utility from an
    explicit screenshot command path
  - capture scope: active GENtle window only (not full desktop)
  - output: caller-provided custom filename/path for saved image artifact
  - purpose: automate documentation refresh and produce image-backed progress
    updates for third-party communication

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

Phase A core scope is complete.

Remaining Phase A items:

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
- process-protocol export contract (plain-text, technical-assistant-friendly
  step list derived from workflow/operation provenance)
- define adapter contract for screenshot artifacts behind explicit opt-in
  (`--allow-screenshots`) with active-window-only capture semantics

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

- View model contract is not yet formalized
- Some rendering/import/export utilities are still adapter-level contracts
  instead of engine operations
- `import-pool` remains a shared shell/CLI utility contract (no dedicated
  engine operation yet)
- No dedicated engine operation yet for exporting a full run/process as a
  technical-assistant protocol text artifact

## 12. Decision log (concise)

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
  accepted and planned (adapter contract + command path pending)
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
- Defer for now: `import-pool` engine operation until first-class container
  semantics settle for stable import behavior across adapters.
