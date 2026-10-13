# Glen: Privilege-Free Post-BLAST Phase Attribution

## Target And Boundary

Steffen must name the full committed SHA containing this document and
`src/primer_import_trace.rs`; use that exact source in an isolated clean
checkout. If it is not yet available, request the pushed ref or a Git bundle.
Do not substitute a nearby main, the exercise tag or historical `f35a58e6`.
Codex has not compiled or executed this diagnostic patch. It is not a speedup
repair, a new release candidate, a batch scheduler or biological acceptance.

Keep the [75-file profiling evidence](audits/glen_followups_20261008/post_blast_profile_f35a58e6/README.md)
and the original 213-file audit immutable. Retain all admission requirements
and reference/input identities in the
[historical procedure](glen_post_blast_import_profiling_f35a58e6.md#admission-checks),
but attribute the new binary/checks to Steffen's new SHA, not `f35a58e6`.
Glen may reuse his already-restored, verified evidence tree read-only; never
overwrite it, rewrite bound paths/hashes, reprepare references or rerun BLAST.
Use fresh state files for every attempt. Stop if source/input/reference
admission cannot be satisfied. Preserve the dirty primary worktree.

## Focused Checks Before Real Inputs

On the approved external host, record source/tree/lockfile and toolchain/host
identity. Build only the ordinary dev/default-feature CLI with one job and
`CARGO_INCREMENTAL=0`, as in the preceding profile. Do not introduce opt-level
overrides, frame pointers, `gui-profiler` or `package-opt1`. Record effective
settings, binary size/hash and `--version`.

```bash
cargo test -q --locked --lib -j1 primer_import_trace -- --test-threads=1
cargo test -q --locked --lib -j1 engine::primer_specificity_multi::tests:: -- --test-threads=1
cargo test -q --locked --lib -j1 specificity_matrix_ -- --test-threads=1
cargo check -q --locked
```

Confirm the trace tests execute rather than matching zero tests. The synthetic
successful import replay must retain the exact summary and unchanged state with
tracing enabled. Missing-input refusal must preserve the error and state.
Recorder tests cover exact opt-in, inert/disabled spans, errors, counts, bounded
output, sink failure and nested/panic cleanup. The synthetic pending import
also exercises the persistence span without turning incomplete evidence into
a pass. The `specificity_matrix_` filter includes the shared-shell
parser test. Matrix execution also rechecks the merged MCP test's
UTF-8 state-path repair. Native Windows execution remains separately unrun if
the host is Linux; Unix mock-BLAST execution does not certify Windows.

Stop on a compile or focused-test failure and return the complete primary log;
do not repair source or skip an assertion. Rustfmt/source checks by Codex are
not substitutes for these runs.

## Bounded Real Imports

Use the unchanged ACTB/TP53 handoff/manifest pairs and their current prepared
Ensembl-116 cDNA/GRCh38 references. Run sequentially, one untraced and one traced
import per case, each with a 60-second deadline and five-second TERM grace.
Maximum four launches, total 300 seconds including cleanup. Stop on resource
pressure, insufficient evidence space, changed inputs or surviving run-owned
processes. No perf installation, sudo, security-policy changes or retries.

Use the fresh-directory/state setup in the historical procedure. `CMD` contains
the absolute CLI plus `--state`, the initially absent state path and
`primers specificity-multi-import HANDOFF MANIFEST`. These are separate templates;
allocate a new run directory for each. Preserve nonzero exits immediately:

```bash
/usr/bin/time -v -o "$R/time.txt" \
  timeout --signal=TERM --kill-after=5s 60s \
  env -u GENTLE_PRIMER_IMPORT_TRACE "${CMD[@]}" \
  >"$R/stdout.json" 2>"$R/stderr.txt"
rc=$?
printf '%s\n' "$rc" >"$R/exit-code.txt"
```

```bash
/usr/bin/time -v -o "$R/time.txt" \
  timeout --signal=TERM --kill-after=5s 60s \
  env GENTLE_PRIMER_IMPORT_TRACE=1 "${CMD[@]}" \
  >"$R/stdout.json" 2>"$R/stderr.txt"
rc=$?
printf '%s\n' "$rc" >"$R/exit-code.txt"
```

Unset or any value other than exactly `1` leaves tracing disabled. Trace records
are JSON lines on stderr with schema `gentle.primer_import_trace.v1`, numeric
import/span identities and fixed phases. The cap is 512 events per import,
including a final `truncated` marker when exceeded. No paths, sequences or
biological identifiers are recorded. Capture stderr to a regular local file,
not a potentially stalled pipe. Keep other stderr diagnostics intact.

Pair `begin` and `completed|failed|abandoned` by `(import, span)`;
`duration_us` is inclusive wall time, not CPU time. Nested durations overlap:
never sum them into total time or CPU percentages. A timed-out span may retain
only its begin event; identify it as unfinished, not zero/complete. `abandoned`
means early return or unwind without explicit completion, not cooperative
cancellation. `truncated` means later observations are unavailable.
Completed import/reference spans mean the code returned normally, not a
specificity pass; incomplete/optional evidence keeps its existing semantics.

Report annotation hash, read and JSON decode separately from database inspection
and repeated validation; then report lookup, forward/reverse hit evaluation and
amplicon joining if reached. Decode counts are annotation records, hit counts
are pre-compaction HSP interpretations, join counts are candidate comparisons,
output-validation counts are accepted command outputs and import counts are
reference rows. An enclosing interpretation count is retained post-compaction
hits. No per-hit logging or sampling bypass is added.

## Evidence To Return

Return per-check/case pass/fail/timeout/blocked/not-run, exact SHA/tree/lock and
binary hashes, profile/host/toolchain, full commands/environment, wall/user/system
time, RSS, signals/exit, process cleanup and state status/hash. Retain raw
stdout/stderr/timing, paired-span durations/counts, unmatched begins and any
truncation. Compare only same-source/same-profile baselines; tracing overhead is
not a speedup or native responsiveness claim.

If imports finish, retain full summaries and child reports, compare their
content/reference/input bindings and verdicts between trace modes, and reopen
saved states. Disclose any excluded volatile creation timestamps individually;
do not strip scientific fields or assume stdout/state equality after timeout.
Recheck original input hashes and both reference inspection identities after
all attempts. Package only new evidence with per-file hashes and archive
size/SHA/location, separately from the two historical bundles.

Identify the expensive or unfinished wall-clock phase, including repeated-call
counts and observational limits. No dominant CPU stack attribution follows.
No automatic caching/removal of validation, reference alias, terminal-policy
change, hit truncation, fan-out/readiness work, publication or gate closure.
