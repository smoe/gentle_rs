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

This is an offline, source-reviewed walkthrough, not a recorded live GUI
acceptance run. No network or screenshot capture is required. The expected
sequence is derived from the same trace: this is a reproducibility exercise,
not independent biological validation of a construct or a real intended edit.
`confirmed` is a target-level software verdict, not proof of whole-plasmid
correctness, sample purity, clinical suitability, or experimental success.

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

- FASTA entry `trace_demo_construct` is a short `48 bp` expected construct derived from the
  first `48` called bases of the bundled `3100.ab1` trace.
- FASTA entry `trace_demo_baseline` is the same sequence except for one intentional tutorial
  SNP at position `24` (1-based):
  - expected allele = `A`
  - baseline allele = `G`
- that lets one imported trace demonstrate a full-span called-base match and
  one inferred intended-edit checkpoint without a second binary trace fixture.
- the supplied ABI `PCON` confidence at position 24 is `9`, below the engine's
  `20` threshold. Although the called allele matches `A`, this checkpoint must
  remain `low_confidence_or_ambiguous` / `insufficient_evidence`.

The importer recognizes ABI/AB1 (`ABIF` magic bytes) and SCF (`.scf` magic
bytes), not just extensions. Use `Raw Trace Import`, not `Open Sequence...`,
for trace evidence. The repository has no committed public SCF fixture; SCF
coverage is synthetic-only in the parser and shell tests. Your own local SCF
file can use the same import controls, but its content determines its verdict.
Do not rename `3100.ab1` to claim SCF coverage.

Confidence values are preserved from the file when available, not invented or
recalibrated. Full-span/junction confirmation uses called-base alignment and
coverage, not quality-weighted consensus. Variant review uses supplied
confidence/ambiguity checks, but missing confidence is not evidence of high
quality. There is no automatic quality trimming, peak re-calling, mixed-peak
allele-fraction inference, or chromatogram editing. See the
[CLI format and scientific limits](./10-01_sequencing_confirmation_trace_cli.md#format-and-scientific-limits)
for the ABI versus SCF field handling and current confidence threshold.

## What You Will Verify

By the end, you should have confirmed all of these:

- the sequencing-confirmation specialist can import a local `ABI/AB1` trace
  directly
- imported traces remain evidence records instead of mutating project sequences
- a baseline sequence automatically creates one expected-edit checkpoint
- the imported trace can support full-span and junction alignment targets
  without any preloaded read sequence IDs
- the variant row classifies the matching but low-confidence baseline-vs-expected
  SNP as `low_confidence_or_ambiguous`, keeping the overall report
  `insufficient_evidence` rather than silently claiming an intended edit
- the chromatogram pane shows live `A/C/G/T` curves and the expected/baseline
  alleles at the selected locus
- the saved confirmation report appears as a lineage artifact and reopens the
  specialist when selected

## Step 1: Load the Expected Construct and Baseline

GUI:

1. start an already-built GENtle GUI in a fresh project; save unrelated work
   before using `File -> New Project`
2. `File -> Open Sequence...`
3. load
   [`docs/tutorial/inputs/sequencing_confirmation_trace_demo_construct.fa`](./inputs/sequencing_confirmation_trace_demo_construct.fa)
4. `File -> Open Sequence...`
5. load
   [`docs/tutorial/inputs/sequencing_confirmation_trace_demo_baseline.fa`](./inputs/sequencing_confirmation_trace_demo_baseline.fa)

What to verify:

- both sequences open as ordinary DNA sequence windows
- the expected construct project ID is `sequencing_confirmation_trace_demo_construct`
- the baseline/reference project ID is `sequencing_confirmation_trace_demo_baseline`

`Open Sequence...` derives project IDs from filenames, not FASTA headers. The
short `trace_demo_construct` / `trace_demo_baseline` headers may still appear
as sequence names. Unlike the CLI guide's explicit `as_id`, this GUI import
does not assign those short IDs. In a nonempty project, collisions add suffixes;
use the actual imported IDs rather than guessing them.

## Step 2: Save a Temporary Project

Do this before importing the trace so the stored report and lineage artifact are
easy to revisit.

GUI:

1. `File -> Save Project...`
2. create a new empty tutorial directory with the save dialog's new-folder
   control and save `trace_demo_gui.project.gentle.json` there
3. do not overwrite an existing project or reuse another run's trace/report IDs

Keep this directory for the later report exports. Saving now establishes the
project location; save again after import and confirmation to persist them.

## Step 3: Open the Sequencing Confirmation Specialist

Make sure the expected construct window is the active DNA window first.

GUI:

1. focus the expected construct window imported from
   `sequencing_confirmation_trace_demo_construct.fa`
2. `Patterns -> Sequencing Confirmation...`

What to verify:

- the sequencing-confirmation specialist opens as a dedicated window
- the expected construct is already set to `sequencing_confirmation_trace_demo_construct`
- the window contains sections for:
  - evidence inputs
  - raw trace import
  - targets
  - report review
  - imported trace review
  - chromatogram

## Step 4: Fill the Confirmation Setup

GUI:

1. set `Baseline/reference sequence ID` to `sequencing_confirmation_trace_demo_baseline`
2. leave `Read sequence IDs` empty
3. keep `Include full construct span target` enabled
4. in `Junction breakpoints (0-based)`, enter `24`:
   - `flank = 12`
5. set `report id` to `trace_demo_gui_confirm`
6. leave alignment `mode` at `local`, `match` at `2`, `mismatch` at `-3`,
   `gap open` at `-5`, `gap extend` at `-1`, `min identity` at `0.80`,
   `min target coverage` at `1.0`, and reverse-complement trials enabled

What this means:

- the report will include one full-span target,
- one explicit junction target centered on the middle of the construct, and
- one inferred expected-edit checkpoint from the baseline-vs-expected SNP.

Boundary `24` is between 1-based bases 24 and 25, while the SNP itself is at
1-based base 24 (0-based index 23). The junction window is `[12, 36)` in
0-based, end-exclusive coordinates. These are deliberately different coordinate
conventions, not an off-by-one error in the expected edit.

## Step 5: Import the Bundled ABI/AB1 Trace

Use the built-in `Raw Trace Import` box inside the same specialist.

GUI:

1. in `Raw Trace Import`, click `Browse...` beside `trace file` and choose:
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

For this AB1 fixture the confidence, peak and curve arrays are populated; that
is not guaranteed for every accepted file. Inspect import warnings. Older
stored traces without curve arrays can still support confirmation, but must
be re-imported from the source file to display chromatogram curves.

## Step 6: Run Confirmation

GUI:

1. click `Run confirmation`

What to verify:

- overall verdict is `insufficient_evidence`
- the evidence table contains one usable trace-backed row
- the evidence row is selectable and becomes the alignment snapshot focus in
  the saved-report pane
- the per-target table shows:
  - full construct span confirmed
  - explicit junction target confirmed
  - one inferred expected-edit checkpoint remains `insufficient_evidence`

## Step 7: Inspect the Variant and Chromatogram

Use the variant/checkpoint list to jump directly to the tutorial SNP.

GUI:

1. in the variant/checkpoint list, select the inferred expected-edit row
2. inspect the chromatogram pane that focuses on that locus

What to verify:

- the selected row classifies as `low_confidence_or_ambiguous`
- the expected allele is `A`
- the baseline allele is `G`
- the observed trace-supported allele is `A`
- the confidence minimum is `9`; a matching called letter alone does not pass
  the expected-edit checkpoint's confidence check
- the chromatogram pane shows:
  - overlaid `A/C/G/T` curves
  - called bases and peak positions
  - the selected locus centered in the local window

Why this matters:

- the difference from the baseline is an intended change, not automatically a
  contradiction
- this trace's low-confidence call still cannot confirm that intended change;
  retain the unresolved checkpoint rather than changing quality values or
  lowering expectations to force a pass

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
  targets, evidence spans, and variant loci should all line up
  on the expected construct ruler
- clicking the inferred expected-edit marker or the evidence span in that
  overview should keep the chromatogram and alignment panes in sync

This 48 bp called-base match is not expected to have uncovered gaps, but its
low-confidence expected-edit target and variant remain unresolved. Use
`Next unresolved` / `Prev unresolved` to visit that checkpoint. Gap review is
an existing control for other evidence, not a success criterion exercised here.
No primer overlay has been loaded in this tutorial, so primer guidance is not
part of the expected output.

## Step 9: Export and Save the Reviewed Report

1. In `Saved report`, choose `trace_demo_gui_confirm`; use `Refresh reports`
   and `Show selected` if needed.
2. Click `Export JSON...` and choose a new `trace_demo_gui_confirm.json` path
   inside the tutorial directory.
3. Click `Export TSV...` and choose a new `trace_demo_gui_confirm.tsv` path.
4. Optionally use `Copy summary` or `Export summary...` for a Markdown review
   snapshot. With this fixture, expect the low-confidence checkpoint, not an
   empty unresolved list.
5. Save the project again with `File -> Save Project...`, then reopen it and the
   lineage artifact to check that the trace and report survived reload.

The JSON should contain `report_id = trace_demo_gui_confirm`,
`baseline_seq_id = sequencing_confirmation_trace_demo_baseline`,
`trace_ids = [abi_demo_trace_gui]`,
and the `low_confidence_or_ambiguous` variant. The TSV should include confirmed
full-span/junction rows and an insufficient-evidence expected-edit row. These
exports are local review artifacts, not a sequencing certificate.

## Optional Negative Control

The bundled `fake.ab1` file is intentionally malformed.

GUI:

1. keep the specialist open
2. in `Raw Trace Import`, choose:
   [`test_files/fixtures/sequencing_confirmation/fake.ab1`](../../test_files/fixtures/sequencing_confirmation/fake.ab1)
3. use a fresh `trace id`, such as `fake_trace_gui`, then click `Import trace`

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
  [`docs/tutorial/10-01_sequencing_confirmation_trace_cli.md`](./10-01_sequencing_confirmation_trace_cli.md)
- GUI reference section:
  [`docs/gui.md`](../gui.md)
- sequencing-confirmation implementation plan:
  [`docs/sequencing_confirmation_plan.md`](../sequencing_confirmation_plan.md)
- fixture provenance and public benchmark shortlist:
  [`test_files/fixtures/sequencing_confirmation/README.md`](../../test_files/fixtures/sequencing_confirmation/README.md)

## Feedback

If this tutorial is confusing, execution-stale, biologically suspect, or missing a useful figure, please open the matching tutorial issue template and include the context copied from GENtle Help -> Tutorial -> Copy Feedback Context.

- Tutorial title:
- Tutorial/chapter id:
- Step reached:
- Expected vs. actual:
- Interface used: GUI / CLI / Agent Assistant / ClawBio

Paste the Tutorial feedback context here:

```text

```
