# Sequencing Confirmation: Current Contract And Remaining Work

Source/status reconciliation: 2026-10-07. Called-read, raw-trace and trace-aware
confirmation, the GUI specialist, primer guidance and lineage reopening are
implemented. The original staged implementation instructions are preserved in
Git history; this file is not a request to implement trace import again.

## Current Contract

The shared engine's `ConfirmConstructReads` compares an explicit expected
construct with loaded called reads and/or imported `trace_ids`.
`ListSequencingConfirmationReports`, show and export operations reuse
persisted `gentle.sequencing_confirmation_report.v1` evidence.
The `seq-confirm` command family and GUI specialist inspect the same report,
not independently computed verdicts.

The report distinguishes `confirmed`, `contradicted` and
`insufficient_evidence`. Coverage gaps remain inspectable and do not become
negative evidence. Expected target spans/junctions/features/edits must remain
explicit; do not substitute a plausible reference for the actual construct.
Optional baseline/reference context separates intended edits from reversions
without changing the expected construct's primacy.

`ImportSequencingTrace` and `seq-trace import|list|show` parse ABI/AB1/SCF
into separate engine-owned trace records. Trace-aware confirmation uses their
called bases and confidence/channel context. It enriches the same confirmation
report rather than creating a second trace-only verdict model. Current
low-confidence contradiction softening is not a complete ambiguity policy.

The GUI Sequencing Confirmation specialist includes per-evidence alignment and
trace/chromatogram inspection, construct-scale coverage, discrepancy/gap
selection, unresolved navigation and summary export. Stored report artifacts
already reopen from lineage/history. `SuggestSequencingPrimers` provides
existing-hit and saved-report-driven fresh-primer guidance with retained
annealing/read-coverage context; a suggestion is not evidence that a read will
actually cover or confirm a construct.

## Maintained Interfaces And Implementation

Canonical schemas, flags and examples belong in [protocol](protocol.md),
[CLI](cli.md) and [GUI](gui.md). Source ownership:

- `src/engine/io/sequencing_traces.rs`: trace input records and parsing.
- `src/engine/analysis/sequencing_confirmation.rs`: confirmation and primer
  guidance.
- `src/engine/ops/operation_handlers.rs`: shared operation execution.
- `src/main_area_dna/sequencing_confirmation_ui.rs`: thin report inspection.
- `src/engine_shell/command_parsers.rs` and `src/engine_shell.rs`: shared routes.

The [legacy intake reference](legacy_gentle_m_intake.md) retains upstream
algorithm provenance. Do not port wxWidgets dialogs or duplicate confirmation
biology in frontend adapters.

## Remaining Scope

- Define deeper ambiguous-base policies beyond confidence-based softening.
  State how uncertain calls affect support, contradiction and not-evaluable
  outcomes, including both orientations. No universal quality threshold or
  automatic base editing is implied.
- Measure and refine chromatogram navigation on actual retained traces.
  Keep raw channels/calls immutable and confidence/coverage visible.
- Recheck current-candidate native GUI, save/reopen, lineage and adapter parity
  against exact report/input hashes. Existing lineage projection is not pending
  implementation, but its current acceptance still needs evidence.

The [roadmap](roadmap.md) owns priority. De novo assembly, chromatogram editing,
large-cohort NGS/pileup workflows and automatic restriction-site introduction
are not implied by this baseline.

## Deterministic Verification And Provenance

Existing fixtures and recreation/use notes remain in
`test_files/fixtures/sequencing_confirmation/README.md`.
Relevant regression paths cover forward/reverse reads, truncation, junction
differences, exact report ordering/export and trace intake/confirmation.
For an ambiguity-policy change, add independent low-confidence, ambiguous-call
and uncovered controls; keep `contradicted` distinct from
`insufficient_evidence`. Test report persistence, trace-ID binding and
primer-guidance provenance through shared adapters.

This reconciliation did not run Rust tests, native GUI checks or laboratory
validation. Parser availability, a historical test result or a primer proposal
does not confirm a clone without the specific read evidence.
