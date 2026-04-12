# Inspect an Imported Sequencing Trace and Confirm a Construct (GUI Tutorial)

> Type: `GUI walkthrough + shared-engine parity`
> Status: `manual/hybrid`
> Drift note: this page is hand-written, but it is intentionally tied to the
> `Patterns -> Sequencing Confirmation...` specialist, the same shared
> `ImportSequencingTrace` / `ConfirmConstructReads` engine routes used by the
> CLI, and deterministic local fixtures.

This tutorial shows the current GUI-first sequencing-confirmation workflow with
raw trace import, baseline-aware variant classification, and chromatogram
inspection.

You will:

1. load one tiny expected construct and one tiny baseline/reference sequence,
2. import one bundled ABI/AB1 trace directly inside the sequencing-confirmation
   specialist,
3. run confirmation using trace evidence only,
4. inspect the inferred intended-edit checkpoint in the variant list,
5. review the matching chromatogram curves at that locus, and
6. reopen the stored confirmation report from lineage.

Current boundary to keep in mind:

- the expected construct remains the primary truth,
- the baseline is used only to explain whether a supported difference is the
  intended edit or a reversion,
- chromatogram review is now available in the GUI,
- and the same pane now also supports trace-base browsing,
- but it is still not a full chromatogram editor.

## Inputs

This walkthrough uses only committed local files:

- expected construct:
  [`docs/tutorial/inputs/sequencing_confirmation_trace_demo_construct.fa`](./inputs/sequencing_confirmation_trace_demo_construct.fa)
- baseline/reference sequence:
  [`docs/tutorial/inputs/sequencing_confirmation_trace_demo_baseline.fa`](./inputs/sequencing_confirmation_trace_demo_baseline.fa)
- ABI/AB1 trace:
  [`test_files/fixtures/sequencing_confirmation/3100.ab1`](../../test_files/fixtures/sequencing_confirmation/3100.ab1)
- optional malformed negative-control trace:
  [`test_files/fixtures/sequencing_confirmation/fake.ab1`](../../test_files/fixtures/sequencing_confirmation/fake.ab1)

Why this pair works:

- `trace_demo_construct` is a short `48 bp` expected construct derived from the
  first `48` called bases of the bundled `3100.ab1` trace.
- `trace_demo_baseline` is the same sequence except for one intentional tutorial
  SNP at position `24` (1-based):
  - expected allele = `A`
  - baseline allele = `G`
- that lets one imported trace confirm both:
  - a full-span construct match, and
  - one inferred intended edit without needing a second binary trace fixture.

## What You Will Verify

By the end, you should have confirmed all of these:

- the sequencing-confirmation specialist can import a local `ABI/AB1` trace
  directly
- imported traces remain evidence records instead of mutating project sequences
- a baseline sequence automatically creates one expected-edit checkpoint
- the imported trace can confirm the expected construct without any preloaded
  read sequence IDs
- the variant row classifies the baseline-vs-expected SNP as
  `intended_edit_confirmed`
- the chromatogram pane shows live `A/C/G/T` curves and the expected/baseline
  alleles at the selected locus
- the saved confirmation report appears as a lineage artifact and reopens the
  specialist when selected

## Step 1: Load the Expected Construct and Baseline

GUI:

1. start GENtle
2. `File -> Open Sequence...`
3. load
   [`docs/tutorial/inputs/sequencing_confirmation_trace_demo_construct.fa`](./inputs/sequencing_confirmation_trace_demo_construct.fa)
4. `File -> Open Sequence...`
5. load
   [`docs/tutorial/inputs/sequencing_confirmation_trace_demo_baseline.fa`](./inputs/sequencing_confirmation_trace_demo_baseline.fa)

What to verify:

- both sequences open as ordinary DNA sequence windows
- the expected construct ID is `trace_demo_construct`
- the baseline/reference sequence ID is `trace_demo_baseline`

## Step 2: Save a Temporary Project

Do this before importing the trace so the stored report and lineage artifact are
easy to revisit.

GUI:

1. `File -> Save Project As...`
2. save to a temporary location, for example:
   `~/Desktop/trace_demo_gui.project.gentle.json`

## Step 3: Open the Sequencing Confirmation Specialist

Make sure the expected construct window is the active DNA window first.

GUI:

1. click the `trace_demo_construct` sequence window so it is focused
2. `Patterns -> Sequencing Confirmation...`

What to verify:

- the sequencing-confirmation specialist opens as a dedicated window
- the expected construct is already set to `trace_demo_construct`
- the window contains sections for:
  - evidence inputs
  - raw trace import
  - targets
  - report review
  - imported trace review
  - chromatogram

## Step 4: Fill the Confirmation Setup

GUI:

1. set `Baseline/reference sequence ID` to `trace_demo_baseline`
2. leave `Read sequence IDs` empty
3. keep `Include full construct span target` enabled
4. add one explicit junction breakpoint:
   - `breakpoint = 24`
   - `flank = 12`
5. set `report id` to `trace_demo_gui_confirm`

What this means:

- the report will include one full-span target,
- one explicit junction target centered on the middle of the construct, and
- one inferred expected-edit checkpoint from the baseline-vs-expected SNP.

## Step 5: Import the Bundled ABI/AB1 Trace

Use the built-in `Raw Trace Import` box inside the same specialist.

GUI:

1. in `Raw Trace Import`, choose the file:
   [`test_files/fixtures/sequencing_confirmation/3100.ab1`](../../test_files/fixtures/sequencing_confirmation/3100.ab1)
2. set optional `trace id` to `abi_demo_trace_gui`
3. keep `associate with expected construct` enabled
4. keep `add imported trace to current run` enabled
5. click `Import trace`

What to verify:

- the import succeeds without creating a new project sequence
- the imported-trace review pane now lists `abi_demo_trace_gui`
- the trace is already present in the current run inputs
- the review pane shows:
  - called-base preview
  - confidence summary
  - peak summary
  - channel summaries

## Step 6: Run Confirmation

GUI:

1. click `Run confirmation`

What to verify:

- overall verdict is `confirmed`
- the evidence table contains one usable trace-backed row
- the evidence row is selectable and becomes the alignment snapshot focus in
  the saved-report pane
- the per-target table shows:
  - full construct span confirmed
  - explicit junction target confirmed
  - one inferred expected-edit checkpoint confirmed

## Step 7: Inspect the Variant and Chromatogram

Use the variant/checkpoint list to jump directly to the tutorial SNP.

GUI:

1. in the variant/checkpoint list, select the inferred expected-edit row
2. inspect the chromatogram pane that focuses on that locus

What to verify:

- the selected row classifies as `intended_edit_confirmed`
- the expected allele is `A`
- the baseline allele is `G`
- the observed trace-supported allele is `A`
- the chromatogram pane shows:
  - overlaid `A/C/G/T` curves
  - called bases and peak positions
  - the selected locus centered in the local window

Why this matters:

- the difference from the baseline is not treated as suspicious here
- because it is the intended construct change, it is counted as positive
  support instead of contradiction

## Step 8: Reopen the Stored Report from Lineage

This step verifies the closing loop: sequencing confirmation is now a
lineage-visible project artifact, not just a transient specialist session.

GUI:

1. return to the project main window
2. open the `Lineage Graph` tab if it is not already visible
3. locate the sequencing-confirmation analysis artifact for
   `trace_demo_gui_confirm`
4. select it

What to verify:

- the artifact appears as a distinct sequencing-confirmation analysis node
- opening it returns you to the sequencing-confirmation specialist
- the specialist reloads the stored report rather than starting from scratch
- the saved report now also shows a construct overview strip:
  targets, evidence spans, uncovered gaps, and variant loci should all line up
  on the expected construct ruler
- clicking the inferred expected-edit marker or the evidence span in that
  overview should keep the chromatogram and alignment panes in sync
- clicking a coverage gap should jump the review focus to the closest
  insufficient-evidence target when one overlaps that unsupported region
- if no target overlaps the clicked gap, the saved-report pane should still
  open an explicit unsupported-region summary for that interval

## Optional Negative Control

The bundled `fake.ab1` file is intentionally malformed.

GUI:

1. keep the specialist open
2. in `Raw Trace Import`, choose:
   [`test_files/fixtures/sequencing_confirmation/fake.ab1`](../../test_files/fixtures/sequencing_confirmation/fake.ab1)
3. click `Import trace`

Expected outcome:

- the import fails with a deterministic format/input error
- the already imported `abi_demo_trace_gui` record remains intact
- the current run inputs are unchanged

## GUI / Shared-Engine Mapping

This tutorial exercises these shared operations through the GUI:

| Tutorial step | GUI action | Shared engine route |
| --- | --- | --- |
| Load expected/baseline | `File -> Open Sequence...` | `LoadFile` |
| Import trace | `Import trace` | `ImportSequencingTrace` |
| Run confirmation | `Run confirmation` | `ConfirmConstructReads` |
| Review stored report | `Show selected` / lineage reopen | `ShowSequencingConfirmationReport` |
| Export report | `Export JSON...` / `Export TSV...` | export report operations |

## Why This Tutorial Matters

This is the first GUI-first path that closes the whole sequencing-confirmation
loop inside GENtle:

- import raw trace evidence,
- classify intended edits versus reversions,
- inspect the chromatogram at the flagged locus, and
- revisit the saved confirmation artifact from lineage.

That makes sequencing confirmation feel like a project-native workflow instead
of a shell-only side path.

## Related Reading

- CLI/shared-shell parity walkthrough:
  [`docs/tutorial/sequencing_confirmation_trace_cli.md`](./sequencing_confirmation_trace_cli.md)
- GUI reference section:
  [`docs/gui.md`](../gui.md)
- sequencing-confirmation implementation plan:
  [`docs/sequencing_confirmation_plan.md`](../sequencing_confirmation_plan.md)
- fixture provenance and public benchmark shortlist:
  [`test_files/fixtures/sequencing_confirmation/README.md`](../../test_files/fixtures/sequencing_confirmation/README.md)
