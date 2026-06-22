# Gene Set Retrieval Producers Plan

Last updated: 2026-06-22

Planning status: design document only. This pass does not implement protocol,
engine, shell, lineage, GUI, or adapter code. Claude review: not run; a
read-only Claude second opinion can be requested before implementation.

## Purpose

This plan designs the next implementation slice for retrieval-backed gene-set
producers: local-cache or explicitly imported retrieval functions that return
genes and normalize them into provenance-rich, graph-addressable
`GeneSetResolutionReport` artifacts.

The target artifact is a resolved logical gene set. It is not a transient
symbol list, not a DNA/protein sequence, and not a wet-lab container. This
extends the parent collection-subject direction in
[`docs/gui_gene_set_collection_operations_plan.md`](gui_gene_set_collection_operations_plan.md):
gene sets are first-class logical collection subjects whose operations must
flow through shared engine/shell/protocol contracts across GUI, CLI, MCP,
JS, Lua, Python, and agent paths.

## Scope Boundaries

Do not implement code in this planning pass.

Do not hide live network retrieval behind existing `gene-sets resolve`,
collection, promoter, or CUT&RUN routes. V1 producer execution is offline and
local-cache friendly. Live GO download, live Ensembl ortholog/paralog sync, or
other provider calls are later explicit opt-in importers.

Keep these four concepts distinct:

- ontology-assigned genes: genes a provider/cache assigns to an ontology term
- local catalog groups citing a term: GENtle gene-group records whose
  `external_mappings` include the term
- enrichment results: statistical comparisons that reference terms but do not
  by themselves define the provider's term membership set
- evidence-derived co-regulated cohorts: candidate sets derived from local
  evidence by explicit scoring/threshold rules, not proof of regulation

Do not force gene sets into lineage `Sequence` nodes or physical `Pool` nodes.
They need a typed logical artifact node when lineage support lands.

## Inspected Context

This plan was reconciled against:

- `crates/gentle-protocol/src/gene_sets.rs`
- `crates/gentle-protocol/src/gene_groups.rs`
- `src/engine/analysis/gene_sets.rs`
- `src/lineage_export.rs`
- `docs/gui_gene_set_collection_operations_plan.md`
- `docs/protocol.md`
- `docs/cli.md`
- `docs/roadmap.md`

The committed baseline at `HEAD` already has gene-set resolution, promoter
cohort, CUT&RUN support schemas, and the cohort-relationship surface
(`GeneSetCohortRelationship`, `relationship_flags`) described in the task. This
dirty worktree additionally contains broader uncommitted producer, importer,
lineage, and documentation edits. Those dirty changes should not be treated as
one merged implementation unit. If this plan is turned into PRs, split or
replay only the phase under review from a clean `main`, or explicitly carve
this worktree back into the phase-sized sequence below.

For clarity, this document distinguishes:

- `HEAD baseline`: what a future smallest safe first PR should assume
- `dirty worktree`: useful implementation evidence, but not accepted scope for
  this planning deliverable

## Current Baseline Inventory

### Gene-Set Resolution Shape

At the `HEAD baseline`, `GeneSetResolutionReport` already carries:

- `schema`, `generated_at_unix_ms`, optional `op_id`, optional `run_id`
- `request: GeneSetRequest`
- optional `genome_id`
- optional `gene_group_catalog_label`
- optional `genome_catalog_label`
- `contributing_group_ids`
- optional `random: GeneSetRandomProvenance`
- requested/resolved/unresolved counts
- `resolved_members`
- `unresolved_members`
- `warnings`
- report-level `provenance: Vec<GeneSetProvenanceRow>`

`GeneSetResolvedMember` already carries:

- `dedup_key`
- `symbol`
- optional `gene_id`
- aliases
- optional chromosome/start/end/strand
- optional `biotype`
- optional `member_status`
- optional `confidence`
- `contributing_group_ids`
- member-level provenance rows

`GeneSetUnresolvedMember` already carries query, reason, source kind, and
optional source id. This is the right place to preserve failed provider/cache
members instead of dropping them.

`GeneSetProvenanceRow` already carries source kind, source id, optional label,
optional path, and optional free-text note. It is useful for report and member
source fragments, but it is not structured enough to encode provider/cache
version, query filters, or scoring rules.

`GeneSetRandomProvenance` already carries genome id/build, gene-index source,
seed, universe size, and foreground exclusion count. Keep this random-specific
rather than generalizing it into producer provenance.

### Missing In The HEAD Resolution Report

The `HEAD baseline` report does not yet carry:

- review status
- provider id/name/version
- cache id/path/version/digest
- organism
- taxon id
- symbol namespace
- structured query metadata
- structured filter metadata
- co-regulated dataset, contrast, scoring, threshold, direction, and
  relationship metadata

`GeneGroupRecord` already carries `organism`, `taxon_id`, and
`symbol_namespace`. Retrieval producers should lift equivalent set-level
context onto `GeneSetResolutionReport`, because a provider/cache-produced
logical set can exist without being a local gene-group catalog record.

The dirty worktree already prototypes these missing report-level fields and
producer structs. Treat that as Phase 1 implementation evidence, not as a
reason to combine later engine, shell, importer, and lineage changes into the
same PR.

### Existing Source Taxonomy

`GeneSetRequest` already has:

- `CatalogGroup`
- `ExplicitMembers`
- `ExternalMapping { namespace, id }`
- `GenomicNeighbors`
- `Random`

`src/engine/analysis/gene_sets.rs::resolve_gene_set` resolves this taxonomy. It
loads optional genome indexes, applies gene-group curation/member gating,
resolves local `ExternalMapping` through catalog `external_mappings`, handles
neighbors/random, deduplicates by resolved identity, and returns one
`GeneSetResolutionReport`.

The producer families overlap existing sources:

- direct gene-list retrieval is retrieval-backed `ExplicitMembers`
- ontology assignment retrieval is not existing `ExternalMapping`; it asks a
  provider/cache for genes assigned to a term
- co-regulated cohort retrieval derives `ExplicitMembers` from evidence rows
  and connects to cohort-relationship metadata

### Relationship And Dirty Worktree Note

`GeneSetCohortRelationship` and `relationship_flags` are already present in the
baseline gene-set promoter/CUT&RUN report schemas, and the dirty worktree
threads them through producer and lineage experiments. That direction fits
producer family 3, but it should stay distinct:

- a co-regulated producer records an expectation and the evidence/filtering
  recipe used to build a candidate set
- promoter/CUT&RUN reports may later emit non-blocking flags that compare
  observed evidence against the declared expectation
- neither should claim that co-membership proves regulation

## Producer Layer Decision

Recommendation: add a producer layer above resolution, not one new
`GeneSetRequest` variant per producer family.

Rationale:

- `GeneSetRequest` is the normalized analysis-facing resolver taxonomy.
- Producers are acquisition and normalization steps. They read local caches or
  importer outputs, select candidate members, and then call existing resolution
  semantics.
- Keeping producers above resolution preserves `gene-sets resolve` behavior and
  avoids duplicate meanings for `ExplicitMembers` and `ExternalMapping`.
- The durable output remains a `GeneSetResolutionReport` with additive
  producer metadata.

Producer mapping:

| Producer family | Resolver input | Why it is not a new source variant |
|---|---|---|
| Direct retrieved gene list | `GeneSetRequest::ExplicitMembers` | The source is still an explicit candidate member list; metadata records provider/cache origin. |
| Ontology assignment cache | `GeneSetRequest::ExplicitMembers` after cache lookup | It asks "which genes are assigned to this term by this cache?", not "which local groups cite this term?" |
| Local catalog groups citing a term | Existing `GeneSetRequest::ExternalMapping` | This remains catalog curation over `external_mappings`; no producer required. |
| Evidence-derived co-regulated cohort | `GeneSetRequest::ExplicitMembers` after thresholding | The producer recipe and relationship expectation live in metadata. |
| Ensembl ortholog/paralog sync | Deferred opt-in importer that writes a cache or reviewable fragment | Not a hidden live resolver source in V1. |

If future replay requires referring to a large immutable producer artifact
without embedding all candidate members, add one generic request form such as
`ProducerArtifact { producer_id, artifact_id }`. Do not add one
`GeneSetRequest` variant per provider.

## Additive Protocol Deltas

All schema deltas should be additive and back-compatible:

- new containers use `#[serde(default)]`
- new optional fields use `Option<T>` plus `skip_serializing_if`
- new enums have conservative defaults
- old `gentle.gene_set_resolution.v1` payloads deserialize with default values
- `gentle.gene_set_resolution.v1` remains the schema id unless a breaking
  payload change is made

### New Types

Candidate types in `crates/gentle-protocol/src/gene_sets.rs`:

```rust
#[serde(rename_all = "snake_case")]
pub enum GeneSetProducerKind {
    DirectGeneList,
    OntologyAssignment,
    CoRegulatedCohort,
}

#[serde(rename_all = "snake_case")]
pub enum GeneSetResolutionReviewStatus {
    Unreviewed,
    Reviewed,
    Included,
    Draft,
    Deprecated,
}

pub struct GeneSetProducerProvenance {
    pub producer_kind: GeneSetProducerKind,
    pub provider_id: String,
    pub provider_label: Option<String>,
    pub provider_version: Option<String>,
    pub cache_id: Option<String>,
    pub cache_path: Option<String>,
    pub cache_version: Option<String>,
    pub cache_digest: Option<String>,
    pub import_op_id: Option<String>,
    pub imported_at_unix_ms: Option<u128>,
}

pub struct GeneSetProducerFilter {
    pub field: String,
    pub operator: String,
    pub value: String,
}

pub struct GeneSetProducerQueryMetadata {
    pub query_kind: String,
    pub query_id: Option<String>,
    pub query_label: Option<String>,
    pub organism: Option<String>,
    pub taxon_id: Option<String>,
    pub symbol_namespace: Option<String>,
    pub filters: Vec<GeneSetProducerFilter>,
}

pub struct GeneSetCoRegulatedProducerMetadata {
    pub dataset_ids: Vec<String>,
    pub contrast_labels: Vec<String>,
    pub condition_labels: Vec<String>,
    pub normalization_method: String,
    pub scoring_method: String,
    pub threshold_rule: String,
    pub sign_direction_rule: String,
    pub relationship: GeneSetCohortRelationship,
    pub interpretation_note: String,
}
```

Recommended defaults:

- `GeneSetProducerKind::DirectGeneList`
- `GeneSetResolutionReviewStatus::Unreviewed`
- `GeneSetCohortRelationship::Unspecified`
- `GeneSetCoRegulatedProducerMetadata::interpretation_note`:
  `"This evidence-derived cohort is a retrieval result and does not prove regulation."`

### Fields Added To GeneSetResolutionReport

Add these fields additively:

```rust
#[serde(default)]
pub review_status: GeneSetResolutionReviewStatus,

#[serde(default, skip_serializing_if = "Option::is_none")]
pub organism: Option<String>,

#[serde(default, skip_serializing_if = "Option::is_none")]
pub taxon_id: Option<String>,

#[serde(default, skip_serializing_if = "Option::is_none")]
pub symbol_namespace: Option<String>,

#[serde(default, skip_serializing_if = "Option::is_none")]
pub producer: Option<GeneSetProducerProvenance>,

#[serde(default, skip_serializing_if = "Option::is_none")]
pub query_metadata: Option<GeneSetProducerQueryMetadata>,

#[serde(default, skip_serializing_if = "Option::is_none")]
pub co_regulated_metadata: Option<GeneSetCoRegulatedProducerMetadata>,
```

Why these live on `GeneSetResolutionReport` rather than only in provenance
rows:

- organism/taxon/namespace describe set-wide interpretation, including
  unresolved rows
- provider/cache version is a replay property for the whole artifact
- query/filter metadata must be machine-readable for parity and UI review
- per-member `GeneSetProvenanceRow` remains appropriate for row-level source
  fragments and duplicate-collapse auditing

Do not add these producer fields to `GeneGroupRecord`. Gene groups remain the
reviewable local catalog/curation surface. A future importer may write
reviewable catalog fragments, but that is a separate promote/import step.

## Producer Families

### Direct Gene-List Retrieval

Meaning: read a local provider/cache artifact that already contains gene
symbols or gene ids for one organism/namespace.

Examples:

- local curated TSV/JSON from a collaborator
- cached provider export
- project resource bundle generated by an explicit importer

Semantics:

- candidate members normalize through existing explicit-member resolution
- unresolved provider rows remain in `unresolved_members`
- provider/cache/query/filter metadata is stored at report level
- member provenance records the cache row or list id when available

Recommended V1 default:

- accept local JSON/TSV only
- no live URL
- require explicit organism/taxon/namespace or a cache manifest supplying them
- reject unsupported cache major versions

### Ontology Assignment Retrieval

Meaning: read a local ontology-assignment cache and return genes assigned to a
term under explicit filters.

This is distinct from existing `ExternalMapping`:

- `GeneSetRequest::ExternalMapping { namespace: "GO", id }` asks which GENtle
  local catalog groups cite the GO term
- ontology assignment retrieval asks which genes the selected provider/cache
  assigns to the term

This is distinct from enrichment:

- ontology assignment retrieval produces a membership set
- enrichment compares foreground/background sets and returns statistical
  results
- enrichment results can reference terms but should not become a membership
  producer unless the user explicitly materializes a provider/cache term set

Recommended V1 default:

- implement only offline local GO/ontology assignment cache lookup
- require provider/cache version plus evidence/filter metadata
- keep live GO download/indexing as a later explicit route such as
  `resources import-ontology-assignment-cache` or
  `resources sync-go-assignments`

### Evidence-Derived Co-Regulated Cohort Retrieval

Meaning: derive candidate gene members from local evidence rows using explicit
dataset, contrast, scoring, threshold, and sign rules.

Required metadata:

- dataset ids
- contrast labels
- condition labels
- normalization method
- scoring method
- threshold rule
- sign/direction rule
- declared relationship expectation
- interpretation note stating that the cohort does not prove regulation

Recommended V1 default:

- `review_status=unreviewed`
- `relationship=co_regulated` only when the route or user flags explicitly
  declare that expectation; otherwise `manual` or `unspecified`
- no language that implies confirmed regulation from co-membership alone
- downstream promoter/CUT&RUN reports may later compare evidence against the
  declared expectation using `relationship_flags`

## Graph And Lineage Direction

At the `HEAD baseline`, `src/lineage_export.rs` has:

- `LineageSvgNodeKind::{Sequence, Pool, Arrangement, Macro, Analysis, OperationHub}`
- `EngineLineageRenderRowKind::{Sequence, Arrangement, Macro, Analysis, OperationHub}`

There is no logical gene-set node kind. Add one rather than overloading
sequence or pool nodes.

The dirty worktree already prototypes `LineageSvgNodeKind::GeneSet`,
`EngineLineageRenderRowKind::GeneSet`, and related subtitle metadata. Treat
that as Phase 5 implementation evidence. It should not be bundled into the
smallest protocol-metadata PR.

Recommended lineage additions:

- `LineageSvgNodeKind::GeneSet`
- `EngineLineageRenderRowKind::GeneSet`
- row metadata:
  - stable gene-set resolution artifact id
  - producer kind
  - resolved member count
  - unresolved member count
  - organism
  - taxon id
  - symbol namespace

Rendering default:

- a distinct logical-set shape, for example a hexagon or compact rounded
  rectangle with a color not used by sequences or physical pools
- subtitle:
  `members=N unresolved=M | producer=ontology_assignment | taxon=9606`

Lineage should connect:

- producer operation or importer operation -> gene-set artifact
- gene-set artifact -> promoter cohort, CUT&RUN support, and future collection
  operations
- gene-set artifact -> materialized sequences only when materialization is an
  explicit operation

Persisted or indexed gene-set artifacts need stable ids. A minimal id policy is
to derive ids deterministically from schema, producer/source context,
generated-at/op id when present, and a digest of resolved/unresolved member
identity. If a report already has a future `artifact_id`, use that instead.

## Candidate Shell And API Routes

Offline producer routes:

```bash
gene-sets produce direct-list --cache CACHE.json_or_tsv [--query LIST_ID] [--genome GENOME_ID] [--provider-id ID] [--provider-version VERSION] [--cache-version VERSION] [--organism NAME|--taxon-id N|--namespace NAMESPACE] [--filter FIELD=VALUE] [--output OUTPUT.json]

gene-sets produce ontology-assignment --cache CACHE.json_or_tsv --term GO:NNNNNNN [--ontology-namespace GO] [--evidence-code CODE] [--genome GENOME_ID] [--provider-id ID] [--provider-version VERSION] [--cache-version VERSION] [--organism NAME|--taxon-id N|--namespace NAMESPACE] [--filter FIELD=VALUE] [--output OUTPUT.json]

gene-sets produce co-regulated --cache CACHE.json_or_tsv --dataset DATASET_ID --contrast LABEL --score METHOD --threshold RULE --direction both|positive|negative [--relationship manual|co_regulated|anti_co_regulated] [--genome GENOME_ID] [--provider-id ID] [--provider-version VERSION] [--cache-version VERSION] [--organism NAME|--taxon-id N|--namespace NAMESPACE] [--filter FIELD=VALUE] [--output OUTPUT.json]
```

Explicit local importer routes:

```bash
resources import-gene-list-cache --input PATH --provider PROVIDER --version VERSION --output CACHE.json
resources import-ontology-assignment-cache --input PATH --namespace GO --provider PROVIDER --version VERSION --output CACHE.json
resources import-co-regulated-cache --input PATH --dataset DATASET_ID --normalization METHOD --output CACHE.json
```

Deferred opt-in live routes:

```bash
resources sync-go-assignments --provider PROVIDER --taxon-id N --output CACHE.json
resources sync-ensembl-homology --species SPECIES --output CACHE.json
```

If any shell surface grows, update:

- `docs/glossary.json`
- `docs/cli.md`
- `docs/protocol.md`
- `docs/agent_interface.md` when MCP dedicated tools are not added
- `docs/gui_cli_mcp_parity.md` through the parity generator
- shell parser tests
- capability/parity tests
- JS/Lua/Python/agent documentation or explicit non-surfacing reasons

## Phased Roadmap

### Phase 0: Planning Document

Scope:

- this file only
- no engine, shell, protocol, lineage, or docs parity changes beyond the plan

Acceptance:

- existing contracts and the dirty worktree are reconciled explicitly
- smallest safe first PR is described but not executed

### Phase 1: Additive Protocol Metadata

Smallest safe first PR.

Scope:

- add producer metadata structs/enums
- add optional/defaulted fields to `GeneSetResolutionReport`
- add serde compatibility tests for old payloads
- update `docs/protocol.md`

Do not include:

- shell routes
- engine producer behavior
- live retrieval
- lineage nodes
- GUI affordances

Green-build acceptance:

- old `gentle.gene_set_resolution.v1` JSON deserializes with default producer
  metadata
- new producer metadata serializes deterministically
- `cargo test -q -p gentle-protocol gene_set`
- `cargo check -q`

### Phase 2: Offline Direct-List Producer

Scope:

- parse local JSON/TSV direct-list cache rows
- convert selected rows to `GeneSetRequest::ExplicitMembers`
- call existing resolver semantics
- attach producer/query/filter metadata
- preserve unresolved provider rows

Acceptance:

- no network calls
- deterministic ordering and dedup match existing resolver behavior
- tests cover provider/cache version, organism/taxon/namespace, filters, and
  unresolved rows

### Phase 3: Offline Ontology Assignment Producer

Scope:

- parse local ontology assignment cache rows
- lookup by namespace/id, initially GO-compatible
- record evidence/filter metadata structurally
- output uses explicit-member resolution after cache lookup

Acceptance:

- tests prove ontology assignment retrieval is distinct from local
  `ExternalMapping`
- zero assignment rows emit warning/unresolved term context
- no live GO download/indexing

### Phase 4: Offline Co-Regulated Cohort Producer

Scope:

- parse local dataset/contrast cache rows
- apply explicit normalization/scoring/threshold/sign rules
- preserve relationship expectation with `GeneSetCohortRelationship`
- record "does not prove regulation" interpretation note

Acceptance:

- tests cover direction/sign filters
- tests cover `relationship=co_regulated` and
  `relationship=anti_co_regulated`
- downstream promoter/CUT&RUN reports can carry the relationship without
  recomputing producer membership

### Phase 5: Graph-Visible Gene-Set Artifacts

Scope:

- add `GeneSet` node kind to lineage SVG and internal render rows
- persist or index produced gene-set resolution artifacts with stable ids
- render logical-set nodes distinctly from sequences and pools
- link producer operation -> gene set -> downstream promoter/CUT&RUN analysis

Acceptance:

- lineage tests prove a gene-set node is not a sequence node
- graph links producer operation to the gene set and the gene set to downstream
  analysis
- GUI/CLI lineage export remains deterministic

### Phase 6: Opt-In Live Importers

Scope:

- explicit `resources sync-*` or `resources import-*` routes only
- importers write reviewable cache artifacts or catalog fragments
- existing `gene-sets resolve` remains offline/local

Acceptance:

- network tests are opt-in and skipped by default in CI
- offline CI uses tiny local fixtures only
- import artifacts carry provider/cache version and reproducibility metadata

## Smallest Safe First PR

Title direction:

`protocol(gene-sets): add retrieval producer metadata`

Files:

- `crates/gentle-protocol/src/gene_sets.rs`
- `crates/gentle-protocol/src/lib.rs` only if re-exports are needed
- `docs/protocol.md`

Scope:

- protocol structs/enums only
- optional/defaulted fields on `GeneSetResolutionReport`
- serde back-compat tests
- protocol note explaining producer layer vs resolver source taxonomy

Not included:

- no engine producer behavior
- no shell routes
- no glossary/parity changes
- no lineage nodes
- no live network code

Acceptance:

- old `gentle.gene_set_resolution.v1` payloads deserialize
- new producer metadata serializes deterministically
- `cargo test -q -p gentle-protocol gene_set`
- `cargo check -q`

## Required Tests

Resolution and provenance:

- producer metadata round-trips through JSON
- provider/cache version and digest are preserved
- organism/taxon/namespace default from cache manifest or explicit flags
- structured filters round-trip without free-text parsing
- dedup still uses existing resolved identity rules
- duplicate-collapse preserves first-seen member plus merged provenance

Offline/cache behavior:

- missing cache returns a typed error
- unsupported cache major version returns a typed error
- producer routes do not perform live network access by default
- live importer tests require explicit opt-in environment variables

Unresolved reporting:

- unknown direct-list symbols become `unresolved_members`
- ontology term with no assignments emits warning plus unresolved term context
- ambiguous genome mappings remain unresolved unless a future explicit policy
  is added

Producer/source reconciliation:

- direct-list output has `request=ExplicitMembers` plus producer metadata
- ontology assignment output does not use `ExternalMapping`
- existing local GO catalog lookup still uses `ExternalMapping`
- co-regulated output records relationship expectation but not regulatory proof

Lineage:

- future gene-set artifacts render as `GeneSet`, not `Sequence`
- downstream promoter cohort links back to the gene-set artifact
- CUT&RUN support links to the same logical set or to the promoter cohort
  derived from it

Parity:

- new shell paths are in `docs/glossary.json`
- MCP exclusions or dedicated tools are documented
- parity matrix is regenerated
- capability/parity tests remain green
- CLI/MCP/JS/Lua/Python/agent surfacing is implemented or explicitly marked
  with a reason

## File-By-File Targets

Phase 1:

- `crates/gentle-protocol/src/gene_sets.rs`
  - add metadata structs/enums
  - add defaulted fields to `GeneSetResolutionReport`
  - add serde back-compat tests
- `crates/gentle-protocol/src/lib.rs`
  - re-export new protocol types if needed
- `docs/protocol.md`
  - document producer metadata and the offline-first producer layer

Phases 2-4:

- `src/engine/analysis/gene_sets.rs`
  - add producer helpers above existing resolver
  - keep `resolve_gene_set` semantics intact
- `src/engine.rs`
  - add operation variants only when a producer route is implemented
- `src/engine/protocol.rs`
  - expose portable request/report records if root engine protocols need them
- `src/engine/ops/operation_handlers.rs`
  - dispatch producer operations through shared helpers
- `src/engine_shell.rs`
  - parse producer routes after protocol metadata exists
- `docs/cli.md`
  - add route docs once shell commands exist
- `docs/glossary.json`
  - add glossary rows for new shell paths
- `docs/agent_interface.md`
  - add MCP exclusions or tool mappings
- `docs/gui_cli_mcp_parity.md`
  - regenerate after glossary/capability changes
- `tests/capability_registry_parity.rs`
  - keep parity gate green
- `tests/mcp_capability_surface.rs`
  - keep MCP tools/exclusions synchronized

Phase 5:

- `src/lineage_export.rs`
  - add `GeneSet` row/node kind and render/subtitle behavior
- `src/engine.rs`
  - store or index stable gene-set artifacts if they become project-visible
- `src/engine/tests.rs`
  - assert producer -> gene-set -> analysis lineage links
- `src/app/*`
  - add future GUI entry points only after shared artifact semantics exist

Phase 6:

- `src/engine_shell.rs`
  - add explicit importer/sync commands only after offline producers are stable
- `src/engine/analysis/*` or a focused resource-import module
  - keep network/download logic outside the resolver
- `docs/cli.md`, `docs/protocol.md`, `docs/agent_interface.md`
  - document opt-in/network semantics and offline CI defaults

## Risks And Recommended Defaults

| Risk or question | Recommended default |
|---|---|
| New `GeneSetRequest` variants vs producer layer | Use a producer layer above resolution; keep `GeneSetRequest` normalized and small. |
| Where organism/taxon/namespace live | Add optional set-level fields to `GeneSetResolutionReport`; copy from cache, explicit flags, or gene group where applicable. |
| Review status vocabulary | Add `GeneSetResolutionReviewStatus` with default `Unreviewed`; retrieved sets are not trusted until reviewed. |
| Cache version compatibility | Require provider/cache version metadata; reject unsupported major versions. |
| Provider row provenance granularity | Store provider/cache metadata at report level and per-member row provenance in existing `GeneSetProvenanceRow`; add row fields later only if tests show ambiguity. |
| GO assignment vs local GO-citing groups | Keep `ExternalMapping` for local catalog groups; use ontology assignment producer metadata over explicit members for provider/cache term membership. |
| Enrichment results as producers | Do not treat enrichment results as producers in this slice; add a separate enrichment report later. |
| Co-regulated wording | Always include a "does not prove regulation" interpretation note. Default review status to `Unreviewed`. |
| Relationship expectation default | Use `Unspecified` unless the route or flags explicitly declare `manual`, `co_regulated`, or `anti_co_regulated`. |
| Graph support | Add `GeneSet` artifact nodes; never overload `Sequence` or physical `Pool`. |
| Live network pressure | Keep live sync as explicit opt-in resource importers; no hidden network calls from resolver or collection operations. |
| Dirty worktree already contains broader implementation | Split into phase-sized PRs or replay from clean main; do not merge producer routes, lineage storage, and relationship flags as one diff. |
