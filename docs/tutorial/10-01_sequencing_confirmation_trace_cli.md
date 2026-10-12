# Confirm a Construct from an Imported Sequencing Trace (CLI Tutorial)

> Type: `CLI walkthrough + shared-shell parity`
> Status: `manual/hybrid`
> Drift note: this page is hand-written, but it is intentionally tied to the
> shared `seq-trace ...` and `seq-confirm ...` shell routes plus deterministic
> local fixtures.

This tutorial shows the shared shell/CLI route for the current raw-trace
sequencing-confirmation workflow.

You will:

1. load one tiny expected construct sequence,
2. import one bundled ABI/AB1 trace,
3. inspect the persisted trace record,
4. run trace-aware construct confirmation against one explicit junction target,
5. inspect the stored confirmation report, and
6. export JSON + TSV artifacts for handoff or regression checks.

The important current boundary is:

- imported traces are stored as evidence inputs,
- confirmation still uses the shared `ConfirmConstructReads` report model, and
- the GUI specialist now uses that same persisted report model; this page stays
  on the shell/CLI route for parity and regression checks.

## Inputs

This walkthrough uses only committed local files:

- expected construct:
  [`docs/tutorial/inputs/sequencing_confirmation_trace_demo_construct.fa`](./inputs/sequencing_confirmation_trace_demo_construct.fa)
- ABI/AB1 trace:
  [`test_files/fixtures/sequencing_confirmation/3100.ab1`](../../test_files/fixtures/sequencing_confirmation/3100.ab1)
- optional malformed negative-control trace:
  [`test_files/fixtures/sequencing_confirmation/fake.ab1`](../../test_files/fixtures/sequencing_confirmation/fake.ab1)

Why this construct works:

- `trace_demo_construct` is a short `48 bp` sequence derived from the first
  `48` called bases of the bundled `3100.ab1` trace.
- That keeps the tutorial deterministic, local, and easy to confirm from one
  real imported trace without adding a second binary tutorial fixture.

## What You Will Verify

By the end, you should have confirmed all of these:

- `seq-trace import` stores one trace record without mutating sequences
- `seq-trace show` exposes called bases, confidence counts, peak counts, and
  sample/run metadata
- `seq-confirm run` can use `--trace-id` directly
- the resulting report keeps trace-backed evidence in the same confirmation
  schema as sequence-backed evidence
- the target support ids now refer to the imported trace id
- JSON and TSV exports come from the same stored report

## Format and Scientific Limits

`seq-trace import` detects content, not the filename extension: ABI/AB1 uses
`ABIF` magic bytes and SCF uses `.scf` magic bytes. Renaming a file is not format
conversion. Use this evidence import route, not ordinary sequence `LoadFile`,
for chromatograms.

ABI import preserves file-supplied `PBAS` calls, optional `PCON` confidence
bytes, `PLOC` peaks, and available channel arrays. Missing confidence or peak
arrays produce warnings; import does not invent quality values. SCF import
also reads called bases, peaks, channels and clip metadata, selecting the
called base's confidence byte from the four supplied values (the maximum for
an ambiguous call). These are source values, not newly calibrated qualities.

There is no committed public SCF fixture in this repository. Existing SCF tests
construct synthetic bytes in
[`src/engine/io/sequencing_traces.rs`](../../src/engine/io/sequencing_traces.rs)
and [`src/engine_shell/tests.rs`](../../src/engine_shell/tests.rs). With your own
reviewed local SCF file, repeat Steps 3-4 with its path and a fresh trace ID;
expect `format = scf`, but do not expect the AB1 demo's confirmation verdict.

Confirmation aligns existing called bases; it does not re-call peaks, edit
chromatograms, automatically quality-trim reads, or interpret mixed peaks as
validated allele fractions. Full-span/junction support uses alignment identity
and target coverage, not a quality-weighted consensus. Variant/checkpoint
classification can mark ambiguous calls or supplied confidence minima below
`20` as `low_confidence_or_ambiguous`; absent confidence is not proof of high
quality and does not itself trigger that threshold.

This expected sequence was copied from the same trace: agreement is a software
plumbing check, not independent biological validation. `confirmed` means the
requested targets passed the selected rules, not whole-plasmid correctness,
sample purity, clinical suitability, or experimental success. Inspect coverage,
discrepancies and missing evidence before drawing scientific conclusions.

## Step 1: Pick a Dedicated State File

Run Bash commands from the repository root using an already-built CLI from
this checkout. Set `GENTLE_CLI` to an absolute binary path if needed. Stop if
the executable check fails; this walkthrough does not build it.

```bash
GENTLE_CLI="${GENTLE_CLI:-$PWD/target/debug/gentle_cli}"
test -x "$GENTLE_CLI"
RUN_DIR=$(mktemp -d "${TMPDIR:-/tmp}/gentle-seq-trace.XXXXXX") || exit 1
STATE="$RUN_DIR/tutorial.gentle.json"
printf 'Tutorial files: %s\n' "$RUN_DIR"
```

Keep this shell open. Everything below uses the same `--state "$STATE"` path.
All inputs are local and all outputs remain under this fresh directory. Do not
delete a fixed state file or reuse an export path to restart the tutorial.

## Step 2: Load the Expected Construct

Import the tiny expected construct under one stable sequence id:

```bash
"$GENTLE_CLI" --state "$STATE" op '{"LoadFile":{"path":"docs/tutorial/inputs/sequencing_confirmation_trace_demo_construct.fa","as_id":"trace_demo_construct"}}' --confirm
```

What to verify:

- the operation succeeds,
- `trace_demo_construct` is now present in the state, and
- no trace evidence has been imported yet.

Optional quick summary:

```bash
"$GENTLE_CLI" --state "$STATE" state-summary
```

## Step 3: Import the Bundled ABI/AB1 Trace

Import the public bundled AB1 fixture into the sequencing-trace evidence store:

```bash
"$GENTLE_CLI" --state "$STATE" shell 'seq-trace import test_files/fixtures/sequencing_confirmation/3100.ab1 --trace-id abi_demo_trace'
```

Read the `import_report` object in the command's JSON output (the full record
is also returned under `trace`). What to verify:

- `trace_id` is `abi_demo_trace`
- `format` is `abi_ab1`
- `called_base_count` is non-zero
- the trace import reports sample/run metadata from the file
- no project sequence was created or mutated just by importing the trace

## Step 4: Inspect the Persisted Trace Record

List the imported traces:

```bash
"$GENTLE_CLI" --state "$STATE" shell 'seq-trace list'
```

Then inspect the specific record:

```bash
"$GENTLE_CLI" --state "$STATE" shell 'seq-trace show abi_demo_trace'
```

What to look for under the output's `trace` object:

- `called_bases` is present and long
- `called_base_confidence_values` and `peak_locations` are populated
- `sample_name` is the file-derived ABI sample label
- `channel_summaries` reports four processed channels

The first 48 calls include low confidence values. At 1-based base 24 the file
supplies `9`, so the GUI companion's baseline-derived expected-edit checkpoint
is insufficient evidence even though this CLI's alignment-only junction can
be confirmed. Do not interpret populated confidence arrays as uniformly good
quality or replace their values to force a positive verdict.

This is the important separation point:

- the trace now exists as evidence input,
- but no confirmation report exists yet.

## Step 5: Run Trace-Aware Confirmation

Now run `seq-confirm` using the imported trace directly as evidence.

This tutorial uses one explicit junction target at boundary `24` (0-based
left-end position, between 1-based bases 24 and 25). Flank `12` covers
`[12, 36)` in 0-based, end-exclusive coordinates, or bases 13-36 in 1-based
coordinates. This is a synthetic checkpoint, not a demonstrated cloning join.

```bash
"$GENTLE_CLI" --state "$STATE" shell 'seq-confirm run trace_demo_construct --trace-id abi_demo_trace --junction 24 --junction-flank 12 --report-id trace_demo_confirm'
```

What to verify:

- the command succeeds without any `--reads` input
- the report id is `trace_demo_confirm`
- `report.overall_status` is `confirmed`
- the target status is `confirmed`
- `support_read_ids` contains `abi_demo_trace`

This is the new behavior the tutorial is meant to exercise:

- imported trace ids can now be used directly in `seq-confirm run`
- the confirmation report keeps trace-backed evidence rows in the same schema
  used for called-read confirmation

## Step 6: Inspect the Stored Confirmation Report

Because an explicit junction was supplied, this CLI run does not automatically
add a full-span target. Omit the junction flags in a separate run with a fresh
report ID to request the default full span. The GUI companion also adds a
baseline-derived intended-edit checkpoint; the reports need not have the same
target count.

Show the persisted report:

```bash
"$GENTLE_CLI" --state "$STATE" shell 'seq-confirm show-report trace_demo_confirm'
```

What to look for under the `report` object:

- `trace_ids` contains `abi_demo_trace`
- `read_seq_ids` is empty in this trace-only example
- `reads[0].evidence_kind` is `trace`
- `reads[0].trace_id` is `abi_demo_trace`
- `reads[0].confirmed_target_ids` contains `junction_1`

This is how GUI trace review reconnects the confirmation row back to its trace
record without inventing a second confirmation report family.

## Step 7: Export the Report Artifacts

Export the stored report as JSON:

```bash
"$GENTLE_CLI" --state "$STATE" shell "seq-confirm export-report trace_demo_confirm \"$RUN_DIR/trace_demo_confirm.json\""
```

Export the target-support TSV:

```bash
"$GENTLE_CLI" --state "$STATE" shell "seq-confirm export-support-tsv trace_demo_confirm \"$RUN_DIR/trace_demo_confirm.tsv\""
```

What to verify:

- `$RUN_DIR/trace_demo_confirm.json` contains the report itself, without the
  shell response's outer `report` wrapper, and the same report id and `trace_ids`
- `$RUN_DIR/trace_demo_confirm.tsv` includes the `junction_1` row
- the TSV support column contains `abi_demo_trace`

## Optional Negative Control

The bundled `fake.ab1` file is intentionally malformed. Importing it should
fail deterministically:

```bash
"$GENTLE_CLI" --state "$STATE" shell 'seq-trace import test_files/fixtures/sequencing_confirmation/fake.ab1 --trace-id fake_trace'
```

Expected outcome:

- the command fails with a deterministic input/format error
- the existing `abi_demo_trace` record remains intact

Run `seq-trace list` again: only `abi_demo_trace` should be present. An explicit
trace ID can replace a previous record on successful import, so use fresh IDs
for different evidence files rather than treating IDs as append-only.

## Engine / Shell Mapping

This tutorial exercises these shared operations and shell routes:

| Tutorial step | Shared route | Engine operation |
| --- | --- | --- |
| Load expected construct | `op '{"LoadFile":...}'` | `LoadFile` |
| Import trace | `seq-trace import ...` | `ImportSequencingTrace` |
| List traces | `seq-trace list` | `ListSequencingTraces` |
| Show trace | `seq-trace show ...` | `ShowSequencingTrace` |
| Confirm from imported trace | `seq-confirm run ... --trace-id ...` | `ConfirmConstructReads` |
| Show stored report | `seq-confirm show-report ...` | `ShowSequencingConfirmationReport` |
| Export report | `seq-confirm export-report ...` | `ExportSequencingConfirmationReport` |
| Export target support TSV | `seq-confirm export-support-tsv ...` | `ExportSequencingConfirmationSupportTsv` |

## Why This Tutorial Matters

Before raw-trace support, sequencing confirmation in GENtle was limited to
already-materialized read sequences.

This tutorial now gives one deterministic local route for checking the next
trust-building step:

- real imported trace evidence can participate in construct confirmation,
- the evidence store stays separate from project sequences, and
- the resulting report is still one shared engine-owned artifact usable across
  CLI, shell, GUI review, and future agent workflows.

## Related Reading

- sequencing-confirmation implementation plan:
  [`docs/sequencing_confirmation_plan.md`](../sequencing_confirmation_plan.md)
- fixture provenance and public benchmark shortlist:
  [`test_files/fixtures/sequencing_confirmation/README.md`](../../test_files/fixtures/sequencing_confirmation/README.md)
- tutorial landing page:
  [`docs/tutorial/README.md`](./README.md)

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
