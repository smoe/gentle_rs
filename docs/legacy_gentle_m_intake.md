# Legacy gentle-m Intake: Status And Residual Ideas

Source review: 2026-10-07. The original intake's sequencing baseline is present
in Rust. This is a provenance/status reference, not an instruction to port it
again or to recreate wxWidgets dialogs.

## Upstream Design Provenance

The original inspection used `GENtle-persons/gentle-m` at
`ad80adfd837b2eebc1bc0e678d49fb88e39f7ba6`, especially:

- `docs/manual.adoc`
- `src/TSequencingAssistantDialog.cpp`, `src/ABItype.cpp`, `src/SCFtype.cpp`
- `src/MiscDialogs.cpp`, `src/TSilmutDialog.cpp`, `src/AutoAnnotate.cpp`
- `src/TLigationDialog.cpp`

That pinned source is historical design evidence; it was not newly fetched
during this review. Retain algorithmic ideas only through deterministic shared
engine reports/operations, never legacy database plumbing or adapter-local
biology.

## Migrated Baseline

- Called-read and trace-aware construct confirmation:
  `ConfirmConstructReads`, persisted report inspection/export and the GUI
  Sequencing Confirmation specialist.
- ABI/AB1/SCF intake: `ImportSequencingTrace` and
  `seq-trace import|list|show`.
- Sequencing-primer overlays and saved-report-guided proposals:
  `SuggestSequencingPrimers` with explicit annealing and read-coverage context.
- Ligation primitives and typed cloning-routine/macro workflows; the legacy
  generic ligation dialog is not a missing Rust requirement.

Current code is in `src/engine/io/sequencing_traces.rs`,
`src/engine/analysis/sequencing_confirmation.rs` and
`src/main_area_dna/sequencing_confirmation_ui.rs`. Canonical interfaces:
[CLI](cli.md), [protocol](protocol.md), [GUI](gui.md).
The [sequencing plan](sequencing_confirmation_plan.md) retains trace review
and policy details; implementation does not imply a clone is confirmed without
its actual read evidence.

## Remaining Intake Ideas

- Synonymous restriction-site **introduction** with amino-acid preservation,
  before/after site counts and a verification-digest prediction is not the
  existing [sequence-design](dna_sequence_optimization_plan.md) motif-avoidance/
  GC contract. Do not mark it implemented just because synonymous editing exists.
- Report-first curated feature-library auto-annotation remains distinct from
  genome annotation import, feature editing and motif scanning. It needs source/
  version binding, overlap/ambiguity rules and a separate confirmed apply step.

These are parking-lot ideas in the [roadmap](roadmap.md), not release gates or
an approved implementation tranche. Keep ambiguous-base policy and real trace/
GUI acceptance separate from parser/template availability.
