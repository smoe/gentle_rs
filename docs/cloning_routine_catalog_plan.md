# Cloning Routine Catalog: Current Contract And Residual Scope

Source review: 2026-10-07. The original catalog/typed-port/macro-graph foundation
is implemented; this file is not an instruction to rebuild it.

## Current Contract

`assets/cloning_routines.json` is the maintained `gentle.cloning_routines.v1`
catalog. Its 22 definitions cover restriction, PCR, sequence transforms, guide
workflows, Gibson, Golden Gate, Gateway, TOPO, TA/GC, In-Fusion and NEBuilder
HiFi. The current entries are marked `implemented`; that status describes the
declared template, not comprehensive protocol support or wet-lab validation.

Each definition supplies identity, family/status, vocabulary, explainability,
a backing template and typed input/output ports. Discover the actual inventory
with `routines list`, `routines explain ROUTINE_ID` and
`routines compare ROUTINE_A ROUTINE_B`; see [CLI contracts](cli.md) rather than
assuming every variant of a named family exists.

`assets/cloning_patterns_catalog/` stores the templates, including a separate
promoter-reporter template that is not a peer cloning-routine family.
[Pattern catalog](cloning_pattern_catalog.md) documents imports and inventory.

Macro execution already records typed `LineageMacroInstance` objects,
resolved inputs/outputs, expanded operations and status. Lineage graph macro
boxes and their edges are implemented. Repeated runs retain distinct instance
IDs; preflight and execution failures remain visible. Shared shell/engine
contracts, not GUI-local biology, own this behavior.

## Remaining Scope

The [roadmap](roadmap.md) Phase B owns future variants and protocol-pack
clarity. Do not mark a whole family missing because a directional, multi-insert
or screening variant lacks a template, or mark a whole family complete because
one template is present.

For each proposed variant, first identify an actual gap in `routines list`
and `macros template-list`. State its junction/product semantics, typed ports,
failure modes and shared-engine preflight. Require deterministic tests for
binding, repeated-instance recording and exact output geometry before adding
it. General vendor-protocol and wet-lab process modeling remains deferred;
catalog presence is not a laboratory success claim.

## Design Provenance

The original vocabulary crosswalk used review identifier
`S1016847825000585` / `PMID:40449799`. This is historical design provenance,
not a new literature assessment or evidence that a GENtle template reproduces
every method in that review.
