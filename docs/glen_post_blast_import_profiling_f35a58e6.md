# Glen: Bounded Post-BLAST Import Profiling

## Frozen Target And Scope

Owner-requested diagnostic, 2026-10-13. Read this procedure from the coordination
commit forwarded by Steffen, but execute GENtle at exactly
`f35a58e6d2bdd7e0a5551e4531c336d7949bfaf0`, tree
`f678fa505010119b0b546db116a82438e07f6405`. Do not substitute current main or the
exercise tag. Preserve the dirty primary checkout, `paper`, `output/`, prepared
references and [original evidence](audits/glen_followups_20261008/candidate_f35a58e6/README.md).
No source repair, push, CI dispatch, package/tag/release or scientific approval.

Use an existing exact-source `gentle_cli`, or rebuild only that entrypoint in a
clean isolated checkout on the approved external host. Match the original
ordinary `dev`/default-feature recipe, not the separate debug-zero GUI harness
or `package-opt1`. Record effective profile/environment and disclose any
toolchain or binary difference; do not introduce optimization, profiling
features, frame-pointer flags or another source revision to improve the result.
Codex runs no local Rust build/check/test. No repeat of the full test suite or
new BLAST search is requested; the existing focused passes remain source-bound.

## Admission Checks

1. Record clean HEAD/tree, `Cargo.lock` SHA-256, toolchain/platform, CLI
   `--version` and SHA-256, relevant Cargo/Rust overrides and profiler versions.
   Lock digest must be
   `ead4421545a125267e33af84aeeea39f0e7ea2f60717a59ef0271d0e71855f1b`.
2. Verify archive size 103,172,261 and SHA-256
   `35e000d086f17348ea76dbaf88712f44b808abea045a72509cabc6d6034a1dcd`.
   Inspect safe membership before extraction; verify all 213 manifest hashes.
   Never overwrite an existing run directory or execute its archived scripts.
3. Restore only the disposable archived `evidence/` tree at its original
   `/tmp/gentle-final-candidate-1Ece9l/evidence/` paths, and only if the entire
   old run root is absent. Catalog, query and TSV paths are content-bound.
   Do not relocate by rewriting JSON, fabricate hashes/execution records,
   reprepare databases, or use symlink aliases. If exact restoration or the
   original references are unavailable, report blocked; a newly generated
   handoff/search needs separate authorization and attribution.
4. Use `artifacts/primer-specificity/bundle` (ACTB: 5,999 HSP rows, 686,525
   bytes) and `bundle-unique` (TP53: 12,202 rows, 1,404,773 bytes). Exclude the
   redundant stricter-policy replay. Recheck queries, child/parent handoffs,
   four completed zero-exit command rows and every TSV size/SHA against the
   retained manifests. Both parent snapshots are null; use a distinct,
   initially absent disposable state path for each import.
5. Confirm current reference identities with exact-source
   `genomes inspect-blast-resource "Ensembl 116 cDNA" --catalog ORIGINAL_CATALOG`
   and `genomes status "Human GRCh38 Ensembl 116" --catalog ORIGINAL_CATALOG`
   (`components.blast_database`), using existing approved BLAST tools. Do not
   execute any preparation hint. Expected cDNA fingerprint:
   `sha256:3159519891bd8c0d4690e9c2f3de8eeb9621152b58c749e3d7d4920f4913263f`,
   465,769 sequences; annotation SHA-1
   `5f571d5775368aeb4e8a88c3ce4b1c8c138e3acf`.
   GRCh38 fingerprint:
   `sha256:90a64ae8ff24fdbab67286bab563c95a6b2a04be37f50ad811474106c1ba3e7e`,
   706 contigs. Keep catalog IDs, assemblies, masking and index kinds distinct.
   Missing/stale fingerprints are a refusal, not permission to bypass admission.

## Limits And Checks

Supervise these sequentially; no parallel imports or unrelated workload.
Maximum five import launches: one 60-second uninstrumented baseline and one
60-second CPU profile per case, plus at most one 10-second ACTB file-open trace.
Each has at most five seconds of termination grace; total import wall budget
is five minutes. Stop early on sufficient evidence, resource pressure, less
than 2 GiB free evidence space, or 100 MiB of evidence from any single attempt.
Do not extend deadlines or keep retrying permission failures.

The following are separate GNU/Linux command templates, not one script. Set
absolute `CLI`, `HANDOFF`, `MANIFEST` paths after admission; set `PROFILE_ROOT`
to approved scratch space outside the archive/reference trees. Choose `CASE`
(`ACTB` or `TP53`) and `MODE` (`baseline`, `cpu` or optional `opens`). Rerun this
setup before each template so it gets fresh state/output paths:

```bash
R="$(mktemp -d "$PROFILE_ROOT/$CASE-$MODE.XXXXXX")" || exit 1
STATE="$R/state.gentle.json"
test ! -e "$STATE" || exit 1
CMD=("$CLI" --state "$STATE" primers specificity-multi-import "$HANDOFF" "$MANIFEST")
```

Capture exit codes even for deadlines. In the supervising shell, temporarily
disable `set -e` for each attempted command, then retain its status immediately;
do not append `|| true` or combine the import with another GENtle command.

Baseline, once per case:

```bash
/usr/bin/time -v -o "$R/time.txt" \
  timeout --signal=TERM --kill-after=5s 60s "${CMD[@]}" \
  >"$R/stdout.json" 2>"$R/stderr.txt"
rc=$?
printf '%s\n' "$rc" >"$R/exit-code.txt"
```

CPU sample, once per case, only if unprivileged perf is already permitted:

```bash
/usr/bin/time -v -o "$R/time.txt" \
  timeout --signal=INT --kill-after=5s 60s \
  perf record -F 99 --call-graph dwarf,8192 -o "$R/perf.data" -- "${CMD[@]}" \
  >"$R/stdout.json" 2>"$R/stderr.txt"
rc=$?
printf '%s\n' "$rc" >"$R/exit-code.txt"
timeout 10s perf report -i "$R/perf.data" --stdio --percent-limit 1 \
  >"$R/perf-report.txt" 2>"$R/perf-report.stderr.txt"
report_rc=$?
printf '%s\n' "$report_rc" >"$R/perf-report-exit-code.txt"
```

Optional ACTB file-open trace, once; no file contents traced:

```bash
/usr/bin/time -v -o "$R/time.txt" \
  timeout --signal=TERM --kill-after=5s 10s \
  strace -f -e trace=openat -s 256 -o "$R/opens.log" -- "${CMD[@]}" \
  >"$R/stdout.json" 2>"$R/stderr.txt"
rc=$?
printf '%s\n' "$rc" >"$R/exit-code.txt"
```

Ordinary GNU timeout status 124 means deadline, not a
failed biological verdict; forced termination/137 must be explicit. Verify all
run-owned processes have exited before the next attempt. If a process survives
the bounded cleanup, stop and report it; do not kill unrelated processes.
Record state absence/presence and raw hashes before/after. A killed process is
not proof of cooperative cancellation or application-level rollback.

Use only already-authorized unprivileged profiler access. No `sudo`, password
prompt, `perf_event_paranoid`/ptrace-policy change, isolation bypass, credentials
or security-setting adjustment. If perf/strace is unavailable or denied, retain
the error and the baseline, mark that diagnostic blocked and stop that tool.
Instrumented times include profiler/timeout overhead; they are not comparable
to uninstrumented native latency. Unresolved stacks remain unresolved.

## Evidence To Return

Return pass/fail/timeout/blocked/not-run for each case/tool, with exact source,
tree, lock and binary/input hashes, host/tool versions, effective build recipe,
complete argument vectors, environment overrides, wall/user/system time, peak
RSS, deadline/signal/exit and process cleanup. Retain stdout/stderr, timing,
raw CPU profile/text summary, file-open counts by relevant path and state
before/after status/hash. If import completes, retain the full summary, child
reports, reference bindings, verdicts and saved/reopened disposable state.
If not, do not infer a specificity verdict from successful BLAST or an empty
output. Rehash the retained handoffs/queries/TSVs and reference identities after
the attempts; report any change rather than repairing it.

Identify the dominant sampled phase, with percentages/path counts and limits:
database fingerprint checks, subject-annotation loading, per-hit prepared
FASTA-index/sequence lookup, terminal mismatch evaluation or amplicon joins.
The repeated index-loading path is a hypothesis, not the required conclusion;
the archived default policy disables full alignment. Do not discard hits,
change policies or weaken provenance to get a faster pass.

Package new evidence separately with a per-file SHA-256 manifest, archive
size/hash and retained location. Preserve the old bundle unchanged. Return
findings before any repair; no publication, native/package acceptance or human
scientific approval is implied.
