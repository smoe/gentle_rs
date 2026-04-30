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

## Step 1: Pick a Dedicated State File

Use a temporary state path so the tutorial does not interfere with your normal
CLI project state.

```bash
STATE=/tmp/gentle_seq_trace_tutorial.gentle.json
rm -f "$STATE"
```

Everything below uses that same `--state "$STATE"` path.

## Step 2: Load the Expected Construct

Import the tiny expected construct under one stable sequence id:

```bash
cargo run --bin gentle_cli -- --state "$STATE" op '{"LoadFile":{"path":"docs/tutorial/inputs/sequencing_confirmation_trace_demo_construct.fa","as_id":"trace_demo_construct"}}' --confirm
```

What to verify:

- the operation succeeds,
- `trace_demo_construct` is now present in the state, and
- no trace evidence has been imported yet.

Optional quick summary:

```bash
cargo run --bin gentle_cli -- --state "$STATE" state-summary
```

## Step 3: Import the Bundled ABI/AB1 Trace

Import the public bundled AB1 fixture into the sequencing-trace evidence store:

```bash
cargo run --bin gentle_cli -- --state "$STATE" shell 'seq-trace import test_files/fixtures/sequencing_confirmation/3100.ab1 --trace-id abi_demo_trace'
```

What to verify:

- `trace_id` is `abi_demo_trace`
- `format` is `abi_ab1`
- `called_base_count` is non-zero
- the trace import reports sample/run metadata from the file
- no project sequence was created or mutated just by importing the trace

## Step 4: Inspect the Persisted Trace Record

List the imported traces:

```bash
cargo run --bin gentle_cli -- --state "$STATE" shell 'seq-trace list'
```

Then inspect the specific record:

```bash
cargo run --bin gentle_cli -- --state "$STATE" shell 'seq-trace show abi_demo_trace'
```

What to look for:

- `called_bases` is present and long
- `called_base_confidence_values` and `peak_locations` are populated
- `sample_name` is the file-derived ABI sample label
- `channel_summaries` reports four processed channels

This is the important separation point:

- the trace now exists as evidence input,
- but no confirmation report exists yet.

## Step 5: Run Trace-Aware Confirmation

Now run `seq-confirm` using the imported trace directly as evidence.

This tutorial uses one explicit junction target centered at base `24` of the
`48 bp` construct:

```bash
cargo run --bin gentle_cli -- --state "$STATE" shell 'seq-confirm run trace_demo_construct --trace-id abi_demo_trace --junction 24 --junction-flank 12 --report-id trace_demo_confirm'
```

What to verify:

- the command succeeds without any `--reads` input
- the report id is `trace_demo_confirm`
- overall status is `confirmed`
- the target status is `confirmed`
- `support_read_ids` contains `abi_demo_trace`

This is the new behavior the tutorial is meant to exercise:

- imported trace ids can now be used directly in `seq-confirm run`
- the confirmation report keeps trace-backed evidence rows in the same schema
  used for called-read confirmation

## Step 6: Inspect the Stored Confirmation Report

Show the persisted report:

```bash
cargo run --bin gentle_cli -- --state "$STATE" shell 'seq-confirm show-report trace_demo_confirm'
```

What to look for in the report payload:

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
cargo run --bin gentle_cli -- --state "$STATE" shell 'seq-confirm export-report trace_demo_confirm /tmp/trace_demo_confirm.json'
```

Export the target-support TSV:

```bash
cargo run --bin gentle_cli -- --state "$STATE" shell 'seq-confirm export-support-tsv trace_demo_confirm /tmp/trace_demo_confirm.tsv'
```

What to verify:

- `/tmp/trace_demo_confirm.json` contains the same report id and `trace_ids`
- `/tmp/trace_demo_confirm.tsv` includes the `junction_1` row
- the TSV support column contains `abi_demo_trace`

## Optional Negative Control

The bundled `fake.ab1` file is intentionally malformed. Importing it should
fail deterministically:

```bash
cargo run --bin gentle_cli -- --state "$STATE" shell 'seq-trace import test_files/fixtures/sequencing_confirmation/fake.ab1 --trace-id fake_trace'
```

Expected outcome:

- the command fails with a deterministic input/format error
- the existing `abi_demo_trace` record remains intact

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
