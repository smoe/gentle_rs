# Gene Set Retrieval Producers: Landed Status And Design History

Last verified: 2026-08-04

Implementation status: the offline direct-list, ontology-assignment, and
co-regulated producer slices are landed across the shared protocol, engine,
shell/direct CLI, persistence, lineage export, tests, and user documentation.
The producer commands remain GUI shell-only rather than prominent GUI
authoring workflows. Statistical enrichment and live provider synchronization
remain separate future work.

## Purpose

This document records the accepted design and current implementation status of
retrieval-backed gene-set producers. It is retained as design history, not as
an instruction to reimplement the shipped producer layer.

The durable output is a provenance-rich `GeneSetResolutionReport`: a logical
gene set, not a transient symbol list, sequence, or physical pool. Producers
acquire or derive candidate members from explicit local inputs and then reuse
the existing gene-set resolver.

## Status Summary

| Area | Current status |
|---|---|
| Additive producer protocol metadata | Landed |
| Offline direct-list producer | Landed |
| Offline ontology-assignment producer | Landed |
| Offline co-regulated producer | Landed |
| Explicit local cache import routes | Landed |
| Typed `GeneSet` lineage artifacts | Landed |
| Prominent GUI producer authoring | Not implemented; GUI shell access only |
| Statistical gene-set enrichment | Not implemented; roadmap item |
| Live GO download/indexing | Deferred explicit importer |
| Live Ensembl ortholog/paralog sync | Deferred; offline ortholog workflows exist |

The authoritative interface status is generated in
[`gui_cli_mcp_parity.md`](gui_cli_mcp_parity.md). The producer commands,
promoter-cohort command, and ortholog commands are GUI shell-only. By contrast,
`gene-sets resolve` and `gene-sets create-pool` have prominent GUI affordances.

## Accepted Design

### Producers Stay Above Resolution

`GeneSetRequest` remains the normalized resolver taxonomy:

- `CatalogGroup`
- `ExplicitMembers`
- `ExternalMapping`
- `GenomicNeighbors`
- `Random`

Retrieval producers are acquisition and normalization steps above that
taxonomy. They select candidate members from a local cache or evidence table,
then resolve those candidates through `ExplicitMembers`. This avoids one new
request variant per provider and preserves existing resolver behavior for
deduplication, unresolved rows, context binding, and provenance.

The three landed producer kinds are:

```rust
pub enum GeneSetProducerKind {
    DirectGeneList,
    OntologyAssignment,
    CoRegulatedCohort,
}
```

### Distinct GO-Shaped Paths

Two paths may mention the same ontology term but answer different questions:

- `GeneSetRequest::ExternalMapping` resolves reviewed local GENtle catalog
  groups whose `external_mappings` cite the term.
- `GeneSetProducerKind::OntologyAssignment` selects genes that a named offline
  provider/cache assigns to the term, then resolves those genes as explicit
  members.

Neither path is statistical enrichment. Ontology assignment supplies a term's
membership set. Enrichment compares a foreground set with an explicit
background universe and reports term-level statistics.

### Offline-First And Explicit Provenance

Producer execution is deterministic and local-cache based. It does not hide a
network fetch inside `gene-sets resolve`, collection operations, promoter
analysis, or CUT&RUN analysis.

Producer-backed reports preserve:

- producer kind and review status
- provider id, label, and version
- cache id, path, version, digest, and import identity
- organism, taxon id, and symbol namespace
- structured query and filter metadata
- unresolved candidate rows and warnings
- per-member source provenance
- for co-regulated cohorts, dataset/contrast/condition, normalization,
  scoring, threshold, sign/direction, and declared relationship metadata

Retrieved membership is not a trust shortcut. The conservative review default
is `unreviewed`, and a co-regulated cohort carries an interpretation note that
co-membership does not prove regulation.

### Logical Artifacts, Not Physical Samples

Produced reports are persisted as logical gene-set artifacts. Lineage export
uses `LineageSvgNodeKind::GeneSet`, keeping these artifacts distinct from
sequences and physical pools. Downstream promoter and CUT&RUN analyses can link
to the same gene-set artifact. Physical material is created only through an
explicit operation such as `gene-sets create-pool`.

## Landed Protocol Surface

The shared protocol defines:

- `GeneSetProducerKind`
- `GeneSetResolutionReviewStatus`
- `GeneSetProducerProvenance`
- `GeneSetProducerFilter`
- `GeneSetProducerQueryMetadata`
- `GeneSetCoRegulatedProducerMetadata`
- additive producer/query/context fields on `GeneSetResolutionReport`

The local cache schemas are:

- `gentle.gene_set_direct_list_cache.v1`
- `gentle.gene_set_ontology_assignment_cache.v1`
- `gentle.gene_set_co_regulated_cache.v1`

These additions retain the `gentle.gene_set_resolution.v1` schema and use
defaulted/optional fields so older reports remain readable.

## Landed Producer Families

### Direct Gene List

`gene-sets produce direct-list` reads a local JSON, TSV, or CSV cache containing
gene symbols or ids. It selects an optional list/query, resolves the candidate
members through existing explicit-member semantics, and preserves provider,
cache, context, query, filter, and unresolved-row provenance.

### Ontology Assignment

`gene-sets produce ontology-assignment` reads a local assignment cache and
selects genes assigned to a specified term under optional evidence and filter
constraints. A term with no assignment rows produces an unresolved term row
and warning rather than a silent successful empty set.

The implementation and tests keep this path distinct from local
`ExternalMapping` catalog lookup.

### Co-Regulated Cohort

`gene-sets produce co-regulated` reads local evidence-derived cohort rows,
filters by dataset/contrast/condition, applies an explicit score threshold and
sign-direction rule, and records the declared relationship expectation.

The report preserves the normalization and scoring recipe. Downstream evidence
may compare observations with the declared expectation through non-blocking
relationship flags, but the producer does not claim confirmed regulation.

## Landed Commands

Producer routes:

```text
gene-sets produce direct-list
gene-sets produce ontology-assignment
gene-sets produce co-regulated
```

Explicit local cache import routes:

```text
resources import-gene-list-cache
resources import-ontology-assignment-cache
resources import-co-regulated-cache
```

Exact arguments and output contracts are documented in [`cli.md`](cli.md) and
[`protocol.md`](protocol.md). The shell routes dispatch shared engine
operations rather than duplicating producer biology in an adapter.

## Validation Coverage

Deterministic tests cover:

- additive protocol serialization and backward-compatible defaults
- all three local cache schemas and unsupported-major-version rejection
- provider/cache/context/query/filter provenance
- unresolved and duplicate member handling
- ontology assignment versus local `ExternalMapping`
- co-regulated threshold, direction, and anti-relationship behavior
- producer shell parsing and engine dispatch
- typed gene-set lineage rendering
- glossary/capability interface parity

Tiny synthetic local fixtures are used; producer tests do not require live
network access.

## Historical Phase Outcomes

The original forward plan was implemented as the following slices:

| Original phase | Outcome |
|---|---|
| Phase 0: design | Superseded by this landed-status record |
| Phase 1: additive protocol metadata | Landed |
| Phase 2: direct-list producer | Landed |
| Phase 3: ontology-assignment producer | Landed |
| Phase 4: co-regulated producer | Landed |
| Phase 5: graph-visible gene-set artifacts | Landed |
| Phase 6: live importers | Local file importers landed; network sync remains deferred |

## Still Open

### Statistical Enrichment

Enrichment should be modeled as a separate analysis report, provisionally
`GeneSetEnrichmentReport`, not as `GeneSetProducerKind::Enrichment`.

Its inputs should include:

- a resolved foreground logical gene set
- an explicit background universe
- identifier namespace and mapping audit
- one ontology/assignment cache with provider and release provenance
- biological-context compatibility policy
- statistical and multiple-testing-correction parameters

Its outputs should include term-level counts, effect/statistical measures,
adjusted values, unresolved-member accounting, warnings, and exact input/cache
bindings. Enrichment must not redefine provider term membership or imply causal
regulation. A later explicit user operation may materialize a selected term's
provider-assigned members, but that remains an ontology-assignment producer
action rather than a side effect of enrichment.

### Live Provider Synchronization

Live GO download/indexing and live Ensembl ortholog/paralog synchronization
remain explicit, opt-in future importer work. They must write reviewable,
versioned local artifacts and must not become hidden resolver network calls.

### GUI Product Decision

Prominent producer authoring in the GUI is a separate product decision. The
current shell-only GUI status is accurate and should not be changed by
documentation cleanup alone.

## Durable Defaults

- Keep `GeneSetRequest` small; producers stay above resolution.
- Default retrieved sets to `unreviewed`.
- Require provider/cache and biological-context provenance.
- Preserve unresolved rows instead of silently dropping them.
- Keep ontology assignment, local catalog external mapping, and enrichment
  semantically distinct.
- Treat co-regulated membership as a declared evidence-derived cohort, not
  proof of regulation.
- Keep gene sets logical until an explicit materialization/pooling operation.
- Keep producer execution offline by default and live synchronization explicit.
