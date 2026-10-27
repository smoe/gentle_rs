## Sequencing Confirmation

<!-- manual-counterparts:start -->
Related CLI reference: [Sequencing Confirmation](../cli/sequencing-confirmation.md).
These topic links do not imply identical feature coverage across interfaces.
<!-- manual-counterparts:end -->

The dedicated `Patterns -> Sequencing Confirmation...` window is the GUI front
end for construct verification from called reads and imported sequencing traces.

Tutorial:

- [`docs/tutorial/10-02_sequencing_confirmation_gui.md`](../tutorial/10-02_sequencing_confirmation_gui.md)

Current scope:

- expected construct is the active/opened sequence window
- optional baseline/reference sequence ID can be supplied to classify intended
  edits versus reference reversions
- read evidence is supplied as already-loaded project sequence IDs
- imported ABI/AB1/SCF evidence can be supplied either by:
  - importing local trace files directly in this window, or
  - reusing already-imported trace IDs from the project evidence store
- target coverage can be:
  - the default full construct span
  - one or more explicit junction checkpoints
  - automatically inferred expected-edit checkpoints when a baseline sequence is
    supplied
- orientation handling follows the shared engine option:
  forward only or forward + reverse-complement trial

Key controls:

- `Baseline/reference sequence ID`: optional sequence id used only to explain
  whether an observed supported difference is the intended edit or a reversion
- `Read sequence IDs`: comma-separated IDs of read sequences already present in
  the project
- `Imported trace IDs`: comma-separated IDs from the sequencing-trace evidence
  store
- `Raw Trace Import`:
  - `trace file`
  - optional `trace id`
  - `associate with expected construct`
  - `add imported trace to current run`
- `Include full construct span target`
- `Junction breakpoints (0-based)` plus `flank`
- alignment scoring knobs:
  `mode`, `match`, `mismatch`, `gap open`, `gap extend`
- confirmation thresholds:
  `min identity`, `min target coverage`
- `report id`

Actions:

- `Run confirmation`
- `Import trace`
- `Refresh reports`
- `Show selected`
- `Export JSON...`
- `Export TSV...`
- `Add selected trace to run` / `Use selected only` from the imported-trace
  review pane

Report review surface:

- overall verdict (`confirmed`, `contradicted`, `insufficient_evidence`)
- construct overview strip for the whole expected sequence:
  - target spans
  - evidence spans
  - uncovered coverage gaps
  - variant loci
  - selected-evidence discrepancy spans
  - clickable sync into the detailed review panes
- optional unresolved-first review mode for target and variant ordering
- selected coverage-gap summary for unsupported regions that are not yet tied
  to an explicit target
- per-target table with kind, span, support counts, and reason
- per-evidence table with sequence-vs-trace kind, chosen orientation,
  usability, identity, and coverage
- imported trace review pane with called-base preview plus confidence/peak
  summaries
- variant/checkpoint list for inferred or explicit expected-edit loci
- chromatogram pane for trace-backed loci:
  - overlaid `A/C/G/T` curves
  - called bases and peak positions
  - expected allele, optional baseline allele, observed allele
  - classification badge:
    `expected_match`, `intended_edit_confirmed`, `reference_reversion`,
    `unexpected_difference`, `low_confidence`, or `insufficient_evidence`
  - trace-base browser mode with called-base slider plus `First/Prev/Next/Last`
    and clip-window jumps
- selected-evidence alignment snapshot for quick inspection without dropping to
  shell
- saved-report evidence rows can now drive that snapshot directly, including
  `Prev` / `Next` navigation across multiple reads or traces
- selected-evidence discrepancy table with one row per mismatch/insertion/deletion

Current limitations:

- older stored traces without raw curve arrays stay usable for confirmation,
  but the GUI asks for re-import before chromatogram curve inspection
- no chromatogram editing or base re-calling yet
- current chromatogram browsing is inspectable by variant or called base, but
  it is still not a full whole-trace editor
- anneal `Tm/GC/hits` ignore non-annealing 5' tails; dimer/structure diagnostics still use full oligo sequence
- `Show report_id` includes top-candidate clamp/dimer diagnostics in the status line

<!-- manual-page-nav:start -->

---

[GUI Manual](../gui.md) | [Previous: Primer, qPCR, and Sequence-Specific RT Design](primer-qpcr-and-sequence-specific-rt-design.md) | [Next: Promoter Cofactors](promoter-cofactors.md)
<!-- manual-page-nav:end -->
