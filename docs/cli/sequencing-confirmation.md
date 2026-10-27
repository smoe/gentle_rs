# Sequencing Confirmation

<!-- manual-counterparts:start -->
Related GUI reference: [Sequencing Confirmation Workspace](../gui/sequencing-confirmation-workspace.md) | [Sequencing Confirmation](../gui/sequencing-confirmation.md).
These topic links do not imply identical feature coverage across interfaces.
<!-- manual-counterparts:end -->

Sequencing-confirmation capability status (called-read plus trace-aware engine baseline):

- `gentle_cli`: supported via shared-shell commands plus the current GUI
  specialist:
  - `seq-trace import PATH [--trace-id ID] [--seq-id ID]`
  - `seq-trace list [SEQ_ID]`
  - `seq-trace show TRACE_ID`
  - `seq-confirm run EXPECTED_SEQ_ID [--baseline BASELINE_SEQ_ID] --reads ID[,ID...] [--trace-id ID]... [--trace-ids ID[,ID...]]`
  - `seq-confirm list-reports [EXPECTED_SEQ_ID]`
  - `seq-confirm show-report REPORT_ID`
  - `seq-confirm export-report REPORT_ID OUTPUT.json`
  - `seq-confirm export-support-tsv REPORT_ID OUTPUT.tsv`
  - `seq-primer suggest EXPECTED_SEQ_ID [--primers ID[,ID...]] [--confirmation-report REPORT_ID] [--min-3prime-anneal-bp N] [--predicted-read-length-bp N]`
  raw-trace intake is backed by `ImportSequencingTrace`,
  `ListSequencingTraces`, and `ShowSequencingTrace`.
  backed by `ConfirmConstructReads`, `ListSequencingConfirmationReports`,
  `ShowSequencingConfirmationReport`, `ExportSequencingConfirmationReport`, and
  `ExportSequencingConfirmationSupportTsv`. Sequencing-primer overlays are
  backed by `SuggestSequencingPrimers`.
  Current scope:
  - expected construct plus already-loaded read sequences
  - imported ABI/AB1/SCF traces can now participate directly in
    `seq-confirm run` without first materializing a project sequence
  - default full-span confirmation when no explicit targets are supplied
  - explicit junction targets via `--junction LEFT_END_0BASED`
    and `--junction-flank N`
  - forward and reverse-complement read evaluation through the same shared
    engine report contract
  - raw ABI/AB1/SCF intake stores trace evidence separately from
    confirmation reports:
    - preserves file-supplied called bases
    - preserves called-base confidence arrays when available
    - preserves peak locations when available
    - preserves raw per-channel chromatogram curves and clip windows when
      available from the source trace
    - keeps one imported trace store that can be listed/shown without running
      confirmation
  - trace-aware confirmation keeps one shared report model:
    - sequence-backed and trace-backed evidence rows land in the same
      confirmation report
    - trace-backed rows expose their `trace_id` and evidence kind
    - target support/contradiction ids can now reflect imported trace ids
    - optional `--baseline` context lets confirmation classify one locus as:
      intended edit, reference reversion, unexpected difference, or
      low-confidence/insufficient evidence
    - `seq-confirm show-report` now exposes the persisted `baseline_seq_id`
      plus per-variant rows for expected edits inferred from the
      baseline-vs-expected diff
  - sequencing-primer suggestion now has two modes through the same
    `SuggestSequencingPrimers` operation:
    - overlay existing primer ids onto the expected construct
    - or, with only `--confirmation-report`, emit fresh primer proposals for
      unresolved loci that still lack a good existing hit
    - when existing primer hits do exist but sit outside the preferred
      sequencing window, the report now keeps both the best existing hit and
      the top fresh proposal
  - GUI parity is now available through `Patterns -> Sequencing Confirmation...`
    and command palette `Sequencing Confirmation`, backed by the same
    `ConfirmConstructReads` / report-store path:
    - optional baseline sequence input is available in the run form
    - the review pane now includes a variant list and chromatogram curve view
      for trace-backed loci
    - sequencing-primer overlays can now be suggested in the same specialist
      from already-loaded primer sequence ids
    - the overlay report can optionally annotate which saved-report targets and
      variant loci fall inside each predicted primer-derived read span
    - unresolved targets and variant loci now also get one shared
      recommendation row naming the best existing primer hit to clarify them
    - older traces without stored curve arrays stay usable for confirmation,
      but the GUI asks you to re-import them before chromatogram review
  Not yet included in this phase:
  - GUI-side raw trace import
  - chromatogram editing or base re-calling
  - full whole-trace browsing beyond the current variant-focused window
  - lineage/artifact projection for confirmation reports
- `gentle_js`: baseline support via `apply_operation` for the same operation
  family.
- `gentle_lua`: baseline support via `apply_operation` for the same operation
  family.
- `gentle_py`: baseline support via `op(...)` for the same operation family.

<!-- manual-page-nav:start -->

---

[CLI Manual](../cli.md) | [Previous: Guide Design](guides.md) | [Next: Promoter Design](promoters.md)
<!-- manual-page-nav:end -->
