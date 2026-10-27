### Sequencing Confirmation Workspace

<!-- manual-counterparts:start -->
Related CLI reference: [Sequencing Confirmation](../cli/sequencing-confirmation.md).
These topic links do not imply identical feature coverage across interfaces.
<!-- manual-counterparts:end -->

- `Patterns -> Sequencing Confirmation...`
  - opens a dedicated called-read construct-confirmation specialist for the
    active/opened expected construct sequence.
  - the same entry is also available from the command palette as
    `Sequencing Confirmation`.
  - current scope mirrors shared-shell `seq-confirm run`:
    - choose one expected construct sequence
    - enter one or more already-loaded read sequence IDs and/or imported
      sequencing-trace IDs
    - keep the default full-span target and/or add explicit junction
      breakpoints (`0-based` left-end positions plus flank bp)
    - run confirmation through the same shared `ConfirmConstructReads`
      operation used by CLI/shell
  - imported trace review is built in:
    - browse imported ABI/AB1/SCF trace records already stored in the project
    - inspect called-base previews, confidence/peak summaries, source path,
      and sample/run metadata
    - append the selected trace directly to the current confirmation inputs
  - baseline-aware variant review is built in:
    - optional baseline/reference sequence id can be supplied in the run form
    - baseline-vs-expected edits become first-class review checkpoints
    - selected trace-backed checkpoints render a chromatogram curve pane with
      expected, baseline, observed, and classification badges
  - whole-trace browsing is built in:
    - the same chromatogram pane can switch between `Variant locus` and
      `Trace base browser`
    - trace-base browsing supports `First`, `Prev`, `Next`, `Last`,
      `Clip start`, and `Clip end` jumps plus a called-base slider
    - the focused base summary shows called base, confidence, peak index, and
      clip-window status
  - saved-report review is built in:
    - choose a persisted report for the current expected construct
    - inspect overall status, target-level outcomes, evidence-level outcomes,
      per-variant rows, and a selected-evidence alignment snapshot
    - a construct overview strip now projects targets, evidence spans, coverage
      gaps, variant loci, and selected-evidence discrepancies onto the expected
      construct span so review can start from the construct rather than from
      one table
    - clicking the overview synchronizes the target, variant, evidence, and
      trace-backed inspection panes; clicking a coverage gap jumps to the most
      relevant insufficient/contradicted target when possible
    - a clicked coverage gap now also stays highlighted and opens a dedicated
      unsupported-region summary, even when no explicit target overlaps it yet
    - the unsupported-region summary now also shows the nearest retained
      evidence edges plus the most relevant existing primer overlays or fresh
      primer proposals that could close that gap next
    - `Prev unresolved` / `Next unresolved` buttons now step the saved-report
      review focus across unresolved targets, problem variants, and coverage
      gaps without requiring manual clicking in each table
    - `Copy summary` and `Export summary...` now hand off a compact Markdown
      unresolved-review snapshot for lab notes, issues, or AI-assisted
      follow-up, including the currently matching sequencing-primer overlay
      guidance when one is loaded for that report
    - an optional unresolved-first review mode keeps contradicted or
      insufficient rows at the top of the target/variant inspection queue
    - the evidence table itself is now the focus picker for alignment and
      discrepancy review, with `Prev` / `Next` navigation for multi-read runs
    - the selected-evidence snapshot now includes a full discrepancy table,
      not only a short text preview
    - export the selected report as JSON or support TSV
    - persisted reports now also appear in the project lineage graph/table as
      analysis artifacts and reopen this specialist on the selected report
  - sequencing-primer overlays are built in:
    - enter already-loaded primer sequence ids in the same specialist
    - suggest predicted read spans from exact 3' anneal hits on the expected
      construct
    - optionally annotate which saved-report targets or variant loci each
      primer-derived span would cover
    - leave primer ids empty to get report-driven fresh primer proposals for
      unresolved loci when you only have the saved confirmation report
    - unresolved targets and variant loci now get a `Recommended Next Primers`
      table that picks the best existing primer hit by 3' distance and coverage
      of additional unresolved loci
    - a `Proposed New Primers` table now suggests fresh primers when no
      existing hit is good enough for the preferred sequencing window
  - selection convenience:
    - the active sequence selection can append its `start`, `midpoint`, or
      `end` as explicit junction breakpoints
  - current limitation:
    - older traces without stored curve arrays must be re-imported before
      chromatogram curves can be reviewed
    - chromatogram review is now browseable by called base, but it is still
      not a full whole-trace editor or base-calling workbench
  - shared UI-intent parity:
    - `ui open sequencing-confirmation`
    - `ui focus sequencing-confirmation`

<!-- manual-page-nav:start -->

---

[GUI Manual](../gui.md) | [Previous: Gibson Workspace](gibson-workspace.md) | [Next: Experiment Planning](experiment-planning.md)
<!-- manual-page-nav:end -->
