# Sequencing Confirmation Implementation Plan

Status: called-read + raw-trace + trace-aware engine/shell baseline plus GUI chromatogram review and sequencing-primer overlays landed (2026-04-02)

Current implementation snapshot:

- shipped:
  - shared-engine `ConfirmConstructReads` / `ListSequencingConfirmationReports`
    / `ShowSequencingConfirmationReport` / export operations
  - persisted sequencing-confirmation report store
  - shared shell/CLI `seq-confirm` run/list/show/export routes
  - deterministic regression coverage for forward, reverse-complement, and
    truncated-read evidence cases
  - raw ABI/AB1/SCF trace import plus shared `seq-trace import|list|show`
  - trace-aware confirmation refinement through the same
    `ConfirmConstructReads` contract using imported `trace_ids`
  - optional baseline/reference context through the same
    `ConfirmConstructReads` contract for intended-edit vs reversion
    classification
  - GUI sequencing-confirmation specialist now accepts imported trace ids and
    includes imported-trace metadata/called-base review plus
    variant-focused chromatogram inspection in the same window
  - sequencing-primer overlay suggestions through the same shared engine/shell
    path plus the GUI specialist
  - unresolved-locus guidance rows now rank the best existing primer hit for
    each non-confirmed target or problem variant
  - saved-report driven fresh-primer proposals can now be emitted through the
    same `SuggestSequencingPrimers` contract when no good existing hit exists
  - persisted confirmation reports now project into lineage/history views and
    reopen the sequencing-confirmation specialist from GUI lineage artifacts

Current follow-up after this plan:

- deeper whole-trace browsing and chromatogram-navigation polish

Purpose: define the first shared-engine sequencing-confirmation workflow for
GENtle so construct validation from read evidence becomes a first-class,
deterministic capability across GUI/CLI/JS/Lua rather than an informal
alignment-only exercise.

This plan refines the roadmap item in `docs/roadmap.md` and the legacy intake
priorities in `docs/legacy_gentle_m_intake.md`.

## Summary

Build sequencing confirmation in four stages:

1. called-read construct confirmation reports
2. GUI specialist and adapter parity
3. raw ABI/AB1/SCF trace import and inspection
4. trace-aware confirmation refinement

The first shipped milestone should answer one bench question well:

> Given an expected construct and one or more sequencing reads, which expected
> junctions/features/edits are confirmed, contradicted, or still uncovered?

## Problem Statement

GENtle can already align sequences and produce structured alignment reports, but
it does not yet expose a first-class workflow for using read evidence to
validate a construct.

Current gap:

- users can compute alignments, but GENtle does not produce a construct-level
  confirmation verdict
- there is no engine-owned report saying which expected regions are confirmed
- there is no GUI specialist oriented around confirmation evidence
- raw Sanger trace import exists only in legacy `gentle-m`, not in the Rust
  rewrite

## Goals

1. Add a deterministic shared-engine confirmation report contract.
2. Preserve one source of truth across GUI/CLI/JS/Lua.
3. Distinguish `contradicted` from `insufficient_evidence`.
4. Make expected construct checkpoints explicit:
   - full-span coverage
   - junction confirmation
   - feature presence
   - expected edit confirmation
   - optional restriction-site checks
5. Keep confirmation inspectable in lineage/artifact views rather than only in
   transient status text.

## Non-goals (initial track)

- chromatogram editing
- de novo contig assembly
- BAM/SAM/pileup workflows
- large-cohort NGS variant calling
- adapter-local confirmation logic

## Architecture Rules For This Track

- Reuse current pairwise-alignment machinery where possible:
  - `AlignSequences`
  - `SequenceAlignmentReport`
- Store confirmation as engine-owned reports, following existing report-store
  patterns used for primer and RNA-read workflows.
- Project confirmation outputs into lineage/artifact views per the
  graph-inspectability rule in `docs/architecture.md`.
- Keep raw traces and confirmation reports separate:
  - traces are evidence inputs
  - confirmation reports are interpreted validation outputs

## Phase 1: Called-Read Construct Confirmation

### Deliverable

A deterministic report workflow that accepts expected construct sequence(s) plus
one or more called read sequences and emits a construct-level confirmation
report.

### Engine Contracts

Add new engine schema constants and report/store types in `src/engine.rs`:

- `gentle.sequencing_confirmation_report.v1`
- `gentle.sequencing_confirmation_reports.v1`
- `gentle.sequencing_confirmation_support_tsv_export.v1`

Add operation enum variants in `src/engine.rs`:

- `ConfirmConstructReads`
- `ListSequencingConfirmationReports`
- `ShowSequencingConfirmationReport`
- `ExportSequencingConfirmationReport`
- `ExportSequencingConfirmationSupportTsv`

Suggested operation shape:

```rust
ConfirmConstructReads {
    expected_seq_id: SeqId,
    read_seq_ids: Vec<SeqId>,
    #[serde(default)]
    report_id: Option<String>,
    #[serde(default)]
    targets: Vec<ConfirmationTargetSpec>,
    #[serde(default)]
    allow_reverse_complement: bool,
    #[serde(default)]
    mode: PairwiseAlignmentMode,
    #[serde(default = "default_pairwise_match_score")]
    match_score: i32,
    #[serde(default = "default_pairwise_mismatch_score")]
    mismatch_score: i32,
    #[serde(default = "default_pairwise_gap_open")]
    gap_open: i32,
    #[serde(default = "default_pairwise_gap_extend")]
    gap_extend: i32,
    #[serde(default)]
    min_identity_fraction: Option<f64>,
    #[serde(default)]
    min_target_coverage_fraction: Option<f64>,
}
```

### Report Model

Add a new report model in `src/engine.rs`:

- `SequencingConfirmationReport`
- `SequencingConfirmationReportSummary`
- `SequencingConfirmationTargetRow`
- `SequencingConfirmationReadRow`
- `SequencingConfirmationDiscrepancyRow`
- `SequencingConfirmationUncoveredInterval`

Key report fields:

- `schema`, `report_id`, `expected_seq_id`, `generated_at_unix_ms`
- `overall_status`
  - `confirmed`
  - `contradicted`
  - `insufficient_evidence`
- `targets[]`
- `reads[]`
- `uncovered_intervals[]`
- `warnings[]`

Target kinds for v1:

- `full_span`
- `junction`
- `feature_presence`
- `expected_edit`
- `restriction_site`

Each target row should include:

- stable `target_id`
- `label`
- `kind`
- expected coordinate or sequence semantics
- `status`
- `supporting_read_ids`
- `contradicting_read_ids`
- short machine-readable reason

Each read row should include:

- `read_seq_id`
- selected orientation
- best embedded `SequenceAlignmentReport`
- covered target intervals
- discrepancy rows:
  - mismatch
  - insertion
  - deletion
- `usable_for_confirmation`

### Adapter Surfaces

CLI/shared shell:

- `seq-confirm run EXPECTED_SEQ_ID --reads ID[,ID...] [--report-id ID]`
- `seq-confirm list-reports [EXPECTED_SEQ_ID]`
- `seq-confirm show-report REPORT_ID`
- `seq-confirm export-report REPORT_ID OUTPUT.json`
- `seq-confirm export-support-tsv REPORT_ID OUTPUT.tsv`

JS/Lua/Python:

- expose via `apply_operation`
- avoid wrapper-only biology logic

### Likely File Touch Points

- `src/engine.rs`
- `src/engine/ops/operation_handlers.rs`
- `src/engine_shell/command_parsers.rs`
- `src/engine_shell.rs`
- `docs/cli.md`
- `src/engine/tests.rs`

Optional extraction target after scaffolding:

- `src/engine/analysis/sequencing_confirmation.rs`

### Lineage / Artifact Visibility

Per `docs/architecture.md`, confirmation must not live only in status text.

Definition of done for this phase:

- a report/artifact node is projected for the confirmation operation
- the artifact node links back to the expected construct and operation id
- CLI and GUI can reopen/show/export the same stored report

### Acceptance Criteria

1. One perfect forward read confirms a target junction.
2. One reverse-complement read confirms the same target when allowed.
3. Truncated coverage yields `insufficient_evidence`, not `contradicted`.
4. A covered unexpected indel or mismatch yields `contradicted`.
5. Two reads covering different checkpoints can produce overall `confirmed`.
6. JSON export is deterministic across repeated runs.

### Default Policy Decisions For Phase 1

To keep phase 1 buildable without reopening design on every edge case, use
these defaults unless later tests force a change:

- default alignment mode: `local`
- default orientation policy: try forward first, then reverse complement when
  `allow_reverse_complement=true`, and keep the better-scoring usable result
- default target generation in GUI:
  - full construct span
  - all explicit user-selected checkpoints
  - no automatic "all features" expansion until the report model is stable
- contradiction rule:
  - only covered contradictory evidence should mark a target as
    `contradicted`
  - missing coverage must remain `insufficient_evidence`
- overall report status:
  - `contradicted` if any required target is contradicted
  - `confirmed` if all required targets are confirmed
  - otherwise `insufficient_evidence`

### Example Phase-1 Operation Payload

```json
{
  "expected_seq_id": "p53_reporter_expected",
  "read_seq_ids": ["read_m13_fwd", "read_sv40_rev"],
  "report_id": "p53_reporter_seqconf_001",
  "targets": [
    {
      "target_id": "left_junction",
      "kind": "junction",
      "label": "vector-insert left junction",
      "left_end_0based_exclusive": 411,
      "right_start_0based": 412
    },
    {
      "target_id": "right_junction",
      "kind": "junction",
      "label": "insert-vector right junction",
      "left_end_0based_exclusive": 1488,
      "right_start_0based": 1489
    }
  ],
  "allow_reverse_complement": true,
  "mode": "local"
}
```

### Example Phase-1 Report Shape

```json
{
  "schema": "gentle.sequencing_confirmation_report.v1",
  "report_id": "p53_reporter_seqconf_001",
  "expected_seq_id": "p53_reporter_expected",
  "overall_status": "confirmed",
  "targets": [
    {
      "target_id": "left_junction",
      "kind": "junction",
      "status": "confirmed",
      "supporting_read_ids": ["read_m13_fwd"],
      "contradicting_read_ids": []
    },
    {
      "target_id": "right_junction",
      "kind": "junction",
      "status": "confirmed",
      "supporting_read_ids": ["read_sv40_rev"],
      "contradicting_read_ids": []
    }
  ],
  "reads": [
    {
      "read_seq_id": "read_m13_fwd",
      "selected_orientation": "forward",
      "usable_for_confirmation": true
    },
    {
      "read_seq_id": "read_sv40_rev",
      "selected_orientation": "reverse_complement",
      "usable_for_confirmation": true
    }
  ],
  "uncovered_intervals": [],
  "warnings": []
}
```

### Phase-1 File-By-File Checklist

1. `src/engine.rs`
   - add schema constants
   - add report/store structs
   - add operation enum variants
   - add normalization helpers for `report_id`
   - add list/get/export helpers mirroring primer/qPCR patterns
2. `src/engine/ops/operation_handlers.rs`
   - resolve reads and expected sequence
   - compute per-read best orientation/alignment
   - classify target rows and overall report status
   - persist report
   - emit artifact-facing messages/result payload
3. `src/engine_shell/command_parsers.rs`
   - parse `seq-confirm ...` commands
4. `src/engine_shell.rs`
   - add shared-shell execution and show/export output formatting
5. `src/engine/tests.rs`
   - add deterministic positive/negative/coverage-gap cases
   - add export determinism test
   - add report-store list/show coverage
6. lineage/artifact path
   - add graph-visible confirmation artifact node projection
7. `docs/cli.md`
   - add CLI contract and examples

### Suggested Phase-1 Test Names

- `test_confirm_construct_reads_confirms_single_junction_forward_read`
- `test_confirm_construct_reads_accepts_reverse_complement_support`
- `test_confirm_construct_reads_marks_truncated_target_as_insufficient_evidence`
- `test_confirm_construct_reads_marks_covered_unexpected_indel_as_contradicted`
- `test_confirm_construct_reads_combines_two_reads_into_overall_confirmed`
- `test_export_sequencing_confirmation_report_is_deterministic`
- `test_list_sequencing_confirmation_reports_orders_by_report_id`

### Suggested Phase-1 Synthetic Fixtures

- `test_files/fixtures/sequencing_confirmation/README.md`
- `test_files/fixtures/sequencing_confirmation/expected_construct_simple.gb`
- `test_files/fixtures/sequencing_confirmation/read_forward_perfect.fa`
- `test_files/fixtures/sequencing_confirmation/read_reverse_perfect.fa`
- `test_files/fixtures/sequencing_confirmation/read_truncated.fa`
- `test_files/fixtures/sequencing_confirmation/read_junction_indel.fa`

## Phase 2: GUI Confirmation Specialist

### Deliverable

A specialist window oriented around construct validation rather than generic
sequence alignment.

### GUI Scope

User flow:

1. choose expected construct
2. select one or more read sequences already present in project state
3. optionally add or edit targets
4. run confirmation
5. inspect:
   - overall verdict
   - target checklist
   - read table
   - discrepancy table
   - export actions

### Initial Placement

Preferred launch points:

- sequence-window context for active construct
- `Patterns` or `Analysis` menu entry
- command palette entry

### Behavior Rules

- GUI must not compute its own verdicts
- GUI may propose default targets from current features/junctions
- saved reports reopen with the same engine-owned content

### Likely File Touch Points

- `src/app.rs`
- `src/main_area_dna.rs`
- `docs/gui.md`
- `docs/tutorial/README.md` once stable enough for user-facing docs

## Phase 3: Raw Trace Import And Inspection

### Deliverable

Engine-owned import of ABI/AB1/SCF evidence with enough metadata to support
inspection now and confirmation integration later.

### Legacy Code Seeds

- `gentle-m` `src/ABItype.cpp`
- `gentle-m` `src/SCFtype.cpp`
- orientation sanity-check idea from `src/TSequencingAssistantDialog.cpp`

### Engine Contracts

Add trace import/report types, likely in `src/engine.rs`:

- `SequencingTraceRecord`
- `SequencingTraceSummary`
- optional `SequencingTraceImportReport`

Add operations:

- `ImportSequencingTrace`
- `ListSequencingTraces`
- `ShowSequencingTrace`

Initial import scope:

- parse called bases
- parse base-position / peak-position support when available
- keep original source path and format provenance
- import trace evidence without mutating construct sequences

### Important Constraint

Phase 3 should not block shipping phase 1.
Called-read confirmation must stand on its own first.

## Phase 4: Trace-Aware Confirmation

### Deliverable

Use trace support to refine evidence interpretation without changing the core
confirmation contract.

Status:

- shipped baseline:
  - `ConfirmConstructReads` accepts imported `trace_ids` alongside
    `read_seq_ids`
  - imported traces contribute their stored called-base strings as evidence
  - confirmation reports now expose evidence kind plus optional `trace_id` for
    each evidence row
  - optional `baseline_seq_id` now lets the engine classify
    intended edits vs reference reversions without changing expected-construct
    primacy
  - confidence-aware contradiction softening is now applied to low-confidence
    trace loci
  - GUI chromatogram review is now attached to trace-backed evidence rows
- still pending within this phase:
  - deeper ambiguous-base policies beyond the current low-confidence softening
  - lineage/artifact projection for confirmation reports

Possible extensions:

- flag low-confidence contradictory calls
- flag ambiguous-base positions
- distinguish covered-but-weak evidence from clean support
- connect one confirmation read row to one trace record by evidence id

### Rule

Do not fork the report model into “called-read confirmation” versus
“trace confirmation”.
Trace support should enrich the same confirmation report.

## Immediate Follow-up: Sequencing-Primer Overlays

This is not part of the minimum confirmation milestone, but it should follow
immediately after phase 1 or phase 2.

Legacy code seed:

- `gentle-m` `src/MiscDialogs.cpp`
  - `findBestMatch`
  - `matchToVector`
  - `addSequencingPrimer`

Planned behavior:

- scan primer collections against expected construct
- emit directional sequencing-primer candidate annotations
- use exact 3' anneal threshold semantics first
- let confirmation UI reuse these overlays to explain expected coverage gaps

## Fixtures And Provenance Plan

Committed fixture requirements from `docs/architecture.md` apply.

Add a dedicated fixture manifest near the new test data, for example:

- `test_files/fixtures/sequencing_confirmation/README.md`

Current shortlist/provenance note:

- `test_files/fixtures/sequencing_confirmation/README.md`
  - deterministic phase-1 baseline: public `U13852` / `pGEX-3X` derived
    synthetic reads
  - real-read exploratory benchmark: `PRJNA1066256` / `SRR27605537`
  - raw-trace parser seed fixtures: Biopython `Tests/Abi` pinned to a specific
    upstream commit

Fixture classes to prepare:

1. Synthetic plasmid + perfect forward read
2. Synthetic plasmid + reverse-complement read
3. Synthetic plasmid + truncated read
4. Synthetic plasmid + junction mismatch / indel read
5. Tiny ABI/AB1/SCF fixtures with exact origin and recreation notes

Preference order:

- synthetic tiny fixtures first for deterministic behavior
- only add real trace fixtures when provenance and redistribution status are
  clear

## Test Matrix

### Phase 1

- report normalization and deterministic serialization
- target-state classification:
  - confirmed
  - contradicted
  - insufficient_evidence
- CLI command parity with direct engine operation
- artifact/lineage projection coverage

### Phase 2

- GUI specialist smoke tests for loading and exporting stored reports
- no adapter-local verdict regressions because verdict logic stays in engine

### Phase 3

- ABI/AB1 parser tests on tiny committed fixtures
- SCF parser tests on tiny committed fixtures
- trace import provenance checks

### Phase 4

- trace-linked contradiction/ambiguity tests
- repeated-run determinism with identical trace input

## Recommended Build Order

1. `src/engine.rs` report schemas, enums, store types, and helpers
2. `src/engine/ops/operation_handlers.rs` confirmation execution path
3. `src/engine_shell/command_parsers.rs` and `src/engine_shell.rs` CLI/shared shell routes
4. `src/engine/tests.rs` deterministic report-classification tests
5. lineage/artifact projection for confirmation reports
6. GUI specialist
7. trace import
8. trace-aware confirmation enrichment

## Phase-1 Exit Criteria Before Starting Trace Import

Do not begin ABI/AB1/SCF intake until all of the following are true:

1. phase-1 synthetic fixtures are committed with provenance notes
2. report classification behavior is locked by deterministic tests
3. CLI show/export routes are stable enough to support fixture-driven review
4. GUI uses the stored engine report instead of recomputing confirmation logic
5. confirmation artifacts are visible in lineage/history

## Risks

- confusing lack-of-coverage with negative evidence
- overfitting the model to Sanger now in a way that blocks later NGS-aligned
  evidence
- letting GUI defaults for targets drift away from CLI/script defaults
- importing raw traces before the report contract is stable and ending up with
  two parallel workflows

## Definition Of Done For The First Useful Milestone

GENtle should let a user:

1. select an expected construct already in project state
2. select one or more sequencing reads already in project state
3. run one shared-engine confirmation operation
4. see a stable target-by-target verdict
5. export the same report from GUI and CLI
6. reopen that report later as a graph-visible artifact
