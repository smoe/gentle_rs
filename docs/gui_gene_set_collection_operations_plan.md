# GUI Gene Set And Collection Operations Plan

Last updated: 2026-06-18

## Current Gap

GENtle has several collection-shaped concepts, but they are not yet presented
as one coherent GUI subject:

- gene sets resolve catalog groups, explicit members, external mappings,
  genomic neighborhoods, or deterministic samples into auditable member rows
- containers and pools represent physical samples that may contain multiple
  molecules
- arrangements preserve semantic order for lanes, racks, plates, labels, and
  gel review
- analysis reports can derive promoter windows, fragments, amplicons, or
  neighboring loci from multiple source members
- alignments compare multiple members together and preserve correspondence

The missing layer is not a gene-set-only editor. The missing layer is a shared
collection subject model that lets GUI, CLI, MCP, JS, Lua, Python, and agent
routes agree on how single-sequence operations lift over a set of sequences.

## Scope

This plan covers the implementation path for:

- a GUI-visible collection subject mental model
- engine-owned operation lifting metadata
- first-class gene-set collection affordances in the GUI
- collection operation launch, readiness, result inspection, and export
- parity tests proving GUI affordances call the same operation/shell/capability
  paths as CLI/MCP/agent routes

## Non-Goals

- Do not add GUI-only biology logic or GUI-only per-member loops.
- Do not treat gene sets as physical samples unless the user explicitly creates
  a pool/container or arrangement from them.
- Do not make all collection operation families prominent in one release slice.
- Do not require live GO, Ensembl ortholog/paralog, or network lookup support
  for V1 collection flows.
- Do not replace existing containers, arrangements, racks, or pool-gel
  contracts with a new parallel model.

## Design Rules

- A gene set is a logical collection source, not a physical object by itself.
- A collection subject must preserve member identity, ordering semantics when
  present, source provenance, warnings, and typed failures.
- Every single-sequence operation that becomes visible on a collection must
  declare one lifting mode:
  - `map`: run once per member and return per-member reports plus an aggregate
    summary
  - `combine`: intentionally treat members as one pool/sample
  - `compare`: consume members together, for example alignment or cohort
    comparison
  - `arrange`: preserve or create order/placement
  - `derive`: create descendant members from each source member
  - `reject`: return a typed reason because collection input is invalid
- GUI controls expose readiness and result summaries for the selected lifting
  mode, then call the same named engine operation, shell command, or typed
  capability route that CLI/MCP/agent surfaces can call.

## Data Model Direction

Prefer extending portable protocol records before adding GUI state.

Planned protocol concepts:

- `CollectionSubjectRef`:
  - project sequence ids
  - container id
  - arrangement id
  - saved gene-set resolution id/path
  - inline gene-set request
  - future alignment/report ids
- `CollectionMemberRef`:
  - stable member id within the collection
  - optional `seq_id`, `container_id`, gene symbol/gene id, or derived region
  - source provenance row(s)
  - optional ordering/placement metadata
- `CollectionLiftPolicy`:
  - operation family/name
  - allowed lifting modes
  - default lifting mode
  - readiness requirements
  - result payload kind
  - typed rejection reasons
- `CollectionOperationReport`:
  - collection subject summary
  - selected lifting mode
  - per-member status rows
  - aggregate counts and warnings
  - links to created reports, sequences, containers, arrangements, or exports

These names are working names. The implementation should reuse existing
records where possible and introduce new schema only when a report needs to be
portable across adapters.

## Operation Lifting Inventory

Start with documented behavior before adding controls:

| Operation family | Initial lifting mode | First useful GUI surface | Engine/shell path |
|---|---|---|---|
| Export FASTA / sequence inspection | `map` or explicit multi-record export | collection action menu | existing sequence/export shell routes |
| Restriction digest | `map` | collection operation launcher | `Digest` / shared shell digest route |
| PCR / primer design | `map` | collection operation launcher, later primer specialist handoff | `DesignPrimerPairs`, qPCR assay routes |
| BLAST / feature scan / TFBS scan | `map` | collection operation launcher | existing scan operations and shell routes |
| Promoter or neighboring-region derivation | `derive` | gene-set resolve/promoter cohort panel | `ResolveGeneSet`, `BuildGeneSetPromoterCohort` |
| Pool / gel render | `combine` only after explicit pool action; otherwise `arrange` | containers/arrangements panels | `ExportPool`, `CreateArrangementSerial`, `RenderPoolGelSvg` |
| Multiple sequence alignment | `compare` | future alignment workspace | future MSA report route |
| Rack/freezer/inventory placement | `arrange` | arrangements/rack panels | existing rack/arrangement operations plus future inventory routes |
| CUT&RUN gene-set support | `derive` then `map`/aggregate over promoter windows | gene-set evidence panel | `InspectCutRunGeneSetRegulatorySupport` |

## Implementation Phases

### Phase 0: design contract

Status: documented; behavior implementation is still pending.

- Add the durable architecture rule for sequence collection subjects.
- Add protocol wording for collection semantics and operation lifting modes.
- Add GUI manual wording for the collection mental model.
- Add this implementation plan and cross-link it from architecture/protocol/GUI
  docs.

Acceptance:

- No code behavior changes.
- Docs clearly distinguish logical sets, pools, arrangements, alignments,
  derived collections, and storage projections.

### Phase 1: capability and command inventory

- Add a generated or hand-maintained inventory of single-sequence operations and
  their planned collection lifting mode.
- Prefer deriving the inventory from glossary, engine operation metadata, and
  existing capability descriptors when possible.
- Mark each row as one of:
  - supported now
  - ready to wire through an existing shared operation
  - needs engine report contract
  - deliberately rejected for collection input
- Keep this inventory adapter-neutral; GUI consumes it, but does not own it.

File targets:

- `docs/protocol.md`
- `docs/gui_cli_mcp_parity.md`
- `crates/gentle-protocol/src/lib.rs`
- `tests/capability_registry_parity.rs`

Tests:

- capability test that every prominent collection GUI affordance has a named
  operation/shell/capability route
- parity freshness test updates if the inventory feeds generated docs

### Phase 2: protocol report and typed rejection baseline

- Add the smallest portable collection report shape needed to return:
  - selected collection subject
  - lifting mode
  - per-member success/error status
  - aggregate warnings
  - links to created artifacts
- Add typed rejection reasons for collection-incompatible operations before any
  GUI collection launcher uses them.
- Avoid introducing a broad new schema until one smoke operation needs it.

File targets:

- `crates/gentle-protocol/src/lib.rs`
- `src/engine/protocol.rs`
- `docs/protocol.md`
- `src/engine/tests.rs`

Tests:

- serialization round-trip for collection report/rejection rows
- one engine unit test proving a rejected collection operation returns a typed
  reason, not a string-only GUI error

### Phase 3: gene-set collection inspector

- Add a GUI inspector for existing gene-set operations:
  - source kind: catalog group, explicit members, external mapping, neighbors,
    deterministic random
  - genome/catalog inputs
  - draft/deprecated gates
  - resolve action
  - resolved/unresolved members
  - warnings and provenance
  - export JSON
- Use `ResolveGeneSet` through the shared operation path.
- Do not duplicate gene-group catalog logic in the GUI.

File targets:

- `src/app.rs`
- `src/app/*_ui.rs` or a new focused `src/app/gene_set_ui.rs`
- `src/engine_shell.rs`
- `src/engine_shell/command_parsers.rs`
- `docs/gui.md`
- `docs/cli.md`

Tests:

- GUI/app helper test that collected form parameters create the same
  `ResolveGeneSet` operation as the shared shell parser
- shell parser test for each source kind used by the GUI
- capability/parity test proving the GUI affordance maps to `gene-sets resolve`
  or `ResolveGeneSet`

### Phase 4: promoter-cohort derivation and arrangement handoff

- From a resolved gene-set collection, offer promoter-window derivation through
  `BuildGeneSetPromoterCohort`.
- Show:
  - returned window count
  - unresolved members
  - promoter span defaults and overrides
  - relationship expectation
  - warnings
- Add explicit handoff choices:
  - materialize/export promoter cohort JSON
  - create/open an arrangement-like ordered view when windows are materialized
    as sequences
  - render one lane per derived member only after an arrangement/container
    exists

File targets:

- `src/app/gene_set_ui.rs` or nearby app UI module
- `src/main_area_dna.rs` only if sequence-window context launch is needed
- `src/engine/protocol.rs`
- `docs/gui.md`
- `docs/protocol.md`

Tests:

- engine test for relationship flag preservation in promoter cohort reports
- GUI helper test that promoter derivation uses `BuildGeneSetPromoterCohort`
- parity test for `gene-sets promoter-cohort` GUI affordance

### Phase 5: collection operation launcher

- Add a generic collection operation launcher that starts with a small curated
  operation set:
  - export FASTA / multi-record artifact
  - digest per member
  - BLAST per member
  - promoter derivation from gene set
  - pool/gel via explicit pool or arrangement action
- The launcher displays the selected lifting mode before execution.
- The launcher must not silently convert logical sets into physical pools.
- Each operation row must show whether it is ready, requires materialization,
  or is unsupported with a typed reason.

File targets:

- new focused GUI module, for example `src/app/collection_operations_ui.rs`
- shared operation metadata in protocol/engine layer
- `src/app/tests.rs`
- `docs/gui.md`

Tests:

- GUI helper tests for readiness states and selected operation payloads
- parity test that launcher operation ids match shared capability descriptors
- regression test preventing GUI-only per-member execution loops for supported
  shared routes

### Phase 6: containers, arrangements, and storage projections

- Allow explicit conversion of logical collections into:
  - a physical pool/container
  - a serial arrangement
  - a rack/storage projection
- Preserve original logical set identity and provenance.
- Require explicit user intent for physical pooling or freezer/rack placement.
- Reuse existing `ExportPool`, `CreateArrangementSerial`,
  `CreateRackFromArrangement`, and rack placement operations.

File targets:

- existing container/arrangement/rack UI modules
- `src/engine.rs`
- `src/engine/ops/operation_handlers.rs`
- `docs/gui.md`
- `docs/protocol.md`

Tests:

- operation tests for provenance links from logical set to container/arrangement
- GUI tests for explicit pool versus one-lane-per-member arrangement choices
- parity tests for all new visible actions

### Phase 7: comparison/alignment workspace

- Add MSA or multi-member comparison only after a shared engine report exists.
- Treat alignment as `compare`, not as a pool or arrangement.
- Preserve member order and aligned-column correspondence.

File targets:

- future alignment engine/report module
- future GUI alignment workspace
- docs/protocol.md
- docs/gui.md

Tests:

- alignment report serialization and deterministic fixture tests
- GUI affordance parity against the shared alignment operation

## Release-Oriented First Slice

The smallest useful implementation PR after this plan should be:

1. Add collection lifting inventory metadata for the three gene-set operations
   and one existing arrangement/pool operation.
2. Add GUI helper code that converts a gene-set resolve form into
   `Operation::ResolveGeneSet`.
3. Add tests proving that the GUI helper and shared shell parser produce the
   same operation payload for one catalog-group source and one explicit-members
   source.
4. Keep display minimal: status, warnings, resolved/unresolved member table,
   and JSON export.

This slice avoids broad collection launchers while proving the core parity
shape.

## Acceptance Criteria For The Full Plan

- Gene sets appear in the GUI as logical collection sources, not as physical
  samples.
- A resolved gene set can feed promoter-cohort derivation, CUT&RUN support,
  export, and later arrangement/materialization flows through shared engine
  operations.
- Prominent GUI collection actions are backed by named operations, shell
  commands, or typed capabilities.
- Tests prove GUI affordances and CLI/MCP/agent-reachable routes do not drift.
- Existing container, pool, arrangement, rack, and gel semantics remain intact.
- Unsupported collection inputs fail with typed reasons instead of hidden GUI
  conditionals.

## Manual Review Still Needed

- Whether the first GUI entry point belongs in the main project window, command
  palette, genome/evidence workspace, or a dedicated collection specialist
  window.
- Whether wet-lab users expect default one-lane-per-member arrangements or an
  explicit prompt before arranging resolved gene sets.
- How much biological context is appropriate in the first gene-set result table
  before the UI becomes a gene-catalog editor.
- Scientific review of relationship-expectation wording for co-regulated and
  anti-co-regulated cohorts.
