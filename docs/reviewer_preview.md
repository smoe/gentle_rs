# Reviewer Quickstart (Internal Preview)

This page is embedded in `Help -> Reviewer Quickstart` and is for colleagues
evaluating GENtle `0.1.0-internal.13`. Reconciled 2026-10-15; it is an inspection
guide, not an acceptance result. Use `Help -> About GENtle` to retain the exact
source revision and build identity of the executable you actually test.

## Internal-Only Warning

Internal preview builds are engineering milestones, not validated assays.

- Not a public/stable release.
- Not intended for production, clinical, or regulated workflows.
- APIs, UI behavior, and file formats may still change.

## Short Offline Command Check

Use an already-built `gentle_cli` from the candidate being reviewed. The commands
below are for a POSIX shell, run from the repository root with that executable
on `PATH`; on Windows use the same workflow paths with separate fresh state
files and `gentle_cli.exe`. Compilation is not part of this walkthrough.

Each example loads its own input. Keep separate state files: combining the
examples would reuse the same sequence ID, not extend the first workflow.

```bash
WORK=$(mktemp -d "${TMPDIR:-/tmp}/gentle-review.XXXXXX")
gentle_cli --version
gentle_cli --state "$WORK/branch.json" workflow @docs/examples/workflows/load_branch_reverse_complement_pgex_fasta.json
gentle_cli --state "$WORK/digest.json" workflow @docs/examples/workflows/digest_ligation_extract_region_minimal.json
```

Expected outcomes:

- Branch workflow creates `pgex_fasta_branch` and `pgex_fasta_branch_rc`.
- Digest/ligation workflow creates `lig_extract`, a 120-base extract.
- Both commands return operation records and save independent state files in
  the fresh directory. Retain them with your feedback; no existing project is
  removed or overwritten.

## GUI Checks

1. Use `File -> Open Tutorial Project...` to load
   `Load FASTA, branch, and reverse-complement`. Inspect the graph, then
   double-click the original or derived sequence to open the DNA viewer.
2. Select a short interval and exercise `Go to selection` / `Zoom to selection`;
   the viewport changes, not the sequence. Resize the window and check that
   wrapped controls remain reachable. See [DNA viewer controls](gui.md).
3. For a separate primer-design check, follow
   [Simple PCR From a Selected Core Region](tutorial/04-01_simple_pcr_selection_gui.md)
   from its starter project, not a completed report. It opens the dedicated
   `PCR Designer`; inspect the selected pair, backend and saved report.
   Allow the tutorial's stated time rather than promising a ten-minute run.
4. Save to a new project path, close, reopen and compare the sequence and
   report identifiers. Report failures with the exact candidate revision.

The same persisted reports can be inspected through `primers list-reports`,
`primers show-report REPORT_ID` and `primers export-report REPORT_ID OUTPUT.json`.
CLI/headless UI intents acknowledge a request; they do not prove a GUI opened.

## Known Limits And Pending Acceptance

- Generated primer pairs are not order-ready without the required specificity,
  oligo-QC, coverage and approval evidence.
- Managed BLAST has asynchronous admission and a bounded queue; migration of
  other long-running commands is incomplete. Resource preparation may need
  explicit external tools and references; this walkthrough downloads none.
- Current-candidate native/independent-package, startup/latency and real-locus
  acceptance remain pending. Published `.12` CI package smokes, historical Linux
  screenshots and nearby revisions do not certify the executable under review.
- Some advanced operations remain shared-shell-first. Proposed persisted primer
  filtered views are not existing commands; see the
  [remaining specialist scope](primer_design_specialist_plan.md).

## Reporting Feedback

When reporting issues, include:

- exact executable version/source revision and OS,
- command/workflow or tutorial step used,
- sequence ID and view span,
- screenshots only through an explicitly approved capture path (if visual),
- exact warning/error text,
- whether behavior reproduced after restart.
