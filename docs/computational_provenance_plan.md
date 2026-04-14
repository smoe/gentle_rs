## Computational Provenance Parity

Last updated: 2026-04-14

This note turns GENtle's general provenance principle into a concrete rollout
plan for persisted computational contributions.

The core question is not only "was a report stored?", but also:

- can a user see that this computational result belongs to a sample?
- can they inspect what inputs and settings produced it?
- can GUI/CLI/agent adapters reopen it consistently from the same project
  state?
- does it appear in lineage/provenance views alongside wet-lab-style sequence
  derivations?

## Current State

GENtle is already strong in three places:

- sequence lineage and operation ids for mutating sequence/container workflows
- immutable operation journal rows and run-bundle export
- a small but real set of first-class analysis artifacts in lineage:
  - dotplots
  - flexibility tracks
  - primer-design reports
  - qPCR-design reports
  - UniProt genome projections
  - sequencing-confirmation reports

The main parity gap is that several other computational stores are preserved
only as metadata/report inventories rather than graph-visible provenance
artifacts. The result is that provenance exists, but users do not encounter it
with the same clarity everywhere.

## Architectural Goal

For persisted computational outputs that may affect sample understanding or
bench action, GENtle should converge on one common standard:

1. The result has a stable artifact id / report id.
2. The result records `op_id` and `run_id`.
3. The result records upstream sequence ids and relevant non-sequence inputs.
4. The result records enough request/effective parameter provenance to explain
   "why this conclusion looked the way it did."
5. The result is reopenable/exportable across GUI/CLI/agent adapters.
6. The result becomes lineage-visible as an analysis artifact node when it is
   persisted and sample-relevant.

This does not mean every helper visualization must be persisted. Transient
helper overlays are still acceptable when all of the following are true:

- they are clearly presented as helper state rather than a final conclusion
- they are deterministically regenerable from persisted inputs/reports
- no user would reasonably treat them as the sole audit trail

## Priority Rollout

### 1. Primer-design and qPCR-design reports

Why first:

- they are strong bench-facing computational conclusions
- users already reason about them like generated assay artifacts, not just
  hidden metadata

Shipped baseline:

- persisted primer/qPCR reports now carry `op_id` / `run_id`
- GUI lineage projects them as analysis artifacts linked from the template
  sequence
- lineage can reopen the PCR Designer on the selected report
- lineage details expose backend provenance plus pair/assay counts

Expected next step:

- extend the same artifact contract to richer primer-adjacent outputs such as
  backend-equivalence audits, nested-PCR follow-ons, or future batch summary
  artifacts when those become persisted review objects

### 2. Projection/annotation analysis artifacts

First candidates:

- UniProt genome projections
- isoform-panel renderable expert views
- grouped TFBS region summaries when persisted

Why:

- these outputs often shape how a sample is interpreted, annotated, or
  communicated
- they are more than raw exports, but they are not yet treated uniformly as
  project-visible computational artifacts

Shipped baseline:

- persisted UniProt genome projections now carry `op_id` / `run_id`
- GUI lineage projects them as analysis artifacts linked from the source
  sequence
- lineage can reopen the UniProt protein expert on the selected projection
- lineage details expose the upstream UniProt `entry_id`, transcript filter,
  and projected transcript count

Expected next step:

- extend the same artifact contract to curated isoform-panel expert views and
  persisted grouped TFBS summaries when those review objects become regular
  sample-facing handoff artifacts

### 3. Protein-side and self-translation outputs

Protein-side work is now starting to enter the same provenance model instead of
becoming a special case later.

Shipped baseline:

- persisted transcript-native protein-derivation reports now carry stable
  `report_id` plus stored `op_id` / `run_id`
- GUI lineage projects those reports as analysis artifacts linked from the
  source nucleotide sequence
- lineage can reopen the transcript-first `Open Derived Protein Expert` path
  directly from the persisted derivation report
- lineage details expose the derivation-mode summary and derived-protein count

Expected next step:

- extend the same artifact contract to reverse/self-translation planning outputs
  so protein-to-DNA handoff reasoning is visible with the same clarity as
  transcript-to-protein derivation

This includes future persisted outputs such as:

- protein-derived reasoning artifacts
- self/reverse-translation planning outputs
- protein-to-DNA handoff artifacts that influence cloning design

The important rule is that protein-side computational conclusions should not
silently bypass the same lineage/artifact visibility expected for nucleotide
workflows.

### 4. Construct reasoning and planning outputs

This track is already emerging through construct-reasoning graphs and
run-bundle exports, but it should converge with the same artifact model.

Shipped baseline:

- persisted construct-reasoning graphs now carry stable `graph_id` plus stored
  `op_id` / `run_id`
- GUI lineage projects them as analysis artifacts linked from the source
  sequence
- lineage can reopen the sequence window on the selected construct-reasoning
  graph
- lineage details expose the objective id/goal plus evidence/decision/candidate
  counts

Expected next step:

- reasoning graphs and planning/decision traces should be inspectable as
  sample-linked computational artifacts rather than only side-channel metadata
- users should be able to tell not only what ran, but why a route or
  interpretation was preferred

### 5. RNA-read interpretation reports

Why fifth:

- they already carry rich request/progress/hit provenance
- they can drive biologically meaningful conclusions
- they deserve first-class lineage presence once the lighter-weight report
  families establish the shared artifact vocabulary cleanly

Shipped baseline:

- persisted RNA-read reports now carry `op_id` / `run_id`
- GUI lineage projects them as analysis artifacts linked from the source
  sequence
- lineage can reopen the RNA-read Mapping workspace directly on the selected
  persisted report
- lineage details expose the stored interpretation profile,
  report-mode/origin-mode summary, total read count, seed-passed count, aligned
  read count, and target-gene count

Follow-on:

- expose linked audit artifacts such as gene-support summaries/audits as either
  subordinate exports or explicit companion artifacts when they become
  important enough for routine review

## Shared Contract Work

To avoid one-off provenance fields per report family, the project should
gradually introduce a reusable computational-artifact provenance shape for
persisted reports/stores.

Minimum desired fields:

- `artifact_id` / `report_id`
- `op_id`
- `run_id`
- `generated_at_unix_ms`
- `primary_seq_id`
- `related_seq_ids[]`
- `external_inputs[]`
  - file path / accession / imported record id
  - checksum when practical
  - optional human label
- `request_summary`
- `effective_settings_summary`
- `reopen_hint`
- `export_kinds[]`

This does not have to become one Rust struct immediately, but new persisted
artifact families should trend toward the same vocabulary instead of inventing
fresh field names each time.

## Acceptance Criteria For One Family

When a computational family is considered provenance-parity complete, it should
meet all of these:

- persisted result can be listed and shown from the shared engine
- result carries stable artifact/report identity plus operation/run linkage
- result records upstream sequence ids
- result records non-sequence inputs strongly enough for audit/replay
- GUI lineage shows it as an analysis artifact node
- lineage can reopen the relevant view/workspace directly
- run-bundle export includes enough information to explain how it was produced

## Non-Goals

- force every ephemeral helper preview into persistent project state
- duplicate the full operation payload into every report family
- treat exported SVG/TSV/JSON files alone as the canonical provenance layer

The canonical provenance should remain project state + operation journal +
run-bundle, with exports as downstream human-readable/materialized views of the
same underlying record.
