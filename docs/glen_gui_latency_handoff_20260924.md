# Glen Handoff: GC Counts And Remaining GUI Latency

Prepared: 2026-09-24. This is an exact-revision diagnostic assignment under
[DEC-039](decisions.md#dec-039-external-auditor-owns-the-performance-verdict),
not a new implementation plan or release approval. The
[latency plan](dna_feature_rendering_latency_plan.md) owns scope and ordering.

## Request To Glen

Please compare these two clean revisions without changing the development
checkout, publishing, tagging or merging:

| Role | Exact revision |
| --- | --- |
| Before | `ad0338a7e5bd14442ca50722a16eb1677cb82790` |
| After | `ee9eb4836752b093764add824617416ea307435a` |

The change replaces whole-sequence GC calculation used only for visible-bin
counts with interval arithmetic. Actual GC values and annotation visibility
are unchanged. The lockfile, Cargo profiles, density fixtures, benchmark and
runner are byte-identical between these revisions. Preserve that controlled
pair even if `main` advances; later revisions need a separately named audit.

The native sections below follow the packaged-GUI measurement contract added
on `main` in [`0e3fd09a`](https://github.com/smoe/gentle_rs/blob/0e3fd09aea82546275efdc274622c97c0d6a1e93/docs/dna_feature_rendering_latency_plan.md#measurement-profiles).
Keep three questions separate: the fixed-revision GC code comparison, native
acceptance of an identified candidate package, and a same-SHA native profile
comparison. The last two do not replace either revision in the GC comparison.

This task has not pushed either revision. If either object is unavailable in
your checkout, request publication/transfer of the exact commit and stop rather
than substituting another tip. Source identity in every binary/receipt must
match its designated revision.

First establish whether the reduction in work has a measurable effect on pan
and zoom. Then attribute remaining startup/first-window and interaction delays
so Codex can select the next small fix. A negligible timing change is a useful
result, not a failed experiment to conceal.

## Scope And Non-Claims

- The nine fixtures at the fixed comparison revisions cover 20 kbp / 250 kbp / 2 Mbp crossed with
  100 / 1,000 / 10,000 features. Do not modify them for this comparison.
- The owner confirmed the interactive target of linear loci up to 250 kbp /
  5,000 features on 2026-09-24. The 2 Mbp and 10,000-feature cases remain stress
  cases. Circular maps require correctness and native interaction regressions,
  but `.12` sets no circular timing budget. The fixed-revision runs below may
  diagnose costs, but cannot certify the agreed target. The exact 5,000-feature
  boundary is absent from this historical ladder and must not be inferred by
  interpolation. Current `main` adds it in `d035912e`, with derived-layer
  inventories, 130 cases and 90 counter observations. Use that workload for a
  separately identified current baseline; do not backport it into this controlled
  pair or claim the 117-case comparison establishes current boundary coverage.
- Keep `.11` packaging/release acceptance separate. This is neither a `.12`
  release verdict nor permission to change scientific output or skip work.
- Developer verification on macOS passed nine GC/layer/cache Rust tests,
  twelve Python runner tests, locked offline `cargo check`, format and whitespace
  checks. It did not establish native GUI or Windows acceptance.

## Controlled Linux CPU Comparison

Use one stable host/toolchain with the same flags, power settings, thread limits
and local storage. Avoid concurrent builds or unrelated heavy jobs. Retain CPU,
RAM, OS, display/session details and cache conditions. Do not compare these
thin-LTO results with historical fat-LTO or Cargo-default release results.

The following Bash uses GNU `/usr/bin/time -v`, local temporary storage and
separate worktrees/targets. Set `repo` to your existing verification clone.
Ensure enough disk space before building; compile cost and peak RSS are evidence,
not startup/runtime measurements. Missing offline dependencies are setup gaps,
not product failures. Do not relax `--locked` or silently fetch during a run.

```bash
set -euo pipefail
repo=/absolute/path/to/gentle_rs
before=ad0338a7e5bd14442ca50722a16eb1677cb82790
after=ee9eb4836752b093764add824617416ea307435a
git -C "$repo" cat-file -e "$before^{commit}"
git -C "$repo" cat-file -e "$after^{commit}"
audit_root=$(mktemp -d "${TMPDIR:-/tmp}/gentle12-audit.XXXXXX")
mkdir -p "$audit_root/worktrees" "$audit_root/targets"
printf 'Retain all evidence under %s\n' "$audit_root"
df -h "$audit_root"
uname -a > "$audit_root/host.txt"

prepare_one() (
    label=$1
    revision=$2
    tree="$audit_root/worktrees/$label"
    git -C "$repo" worktree add --detach "$tree" "$revision"
    cd "$tree"
    test "$(git rev-parse HEAD)" = "$revision"
    test -z "$(git status --porcelain)"
    export CARGO_BUILD_JOBS=1
    export CARGO_TARGET_DIR="$audit_root/targets/$label"
    /usr/bin/time -v -o "$audit_root/$label-build.time.txt" \
        python3 scripts/dna_feature_latency.py prepare \
        --profile bench-audit --timeout 7200 \
        --output "$audit_root/$label-build"
)

replay_one() (
    label=$1
    name=$2
    mode=$3
    /usr/bin/time -v -o "$audit_root/$label-$name.time.txt" \
        python3 "$audit_root/worktrees/$label/scripts/dna_feature_latency.py" run \
        --receipt "$audit_root/$label-build/build.json" \
        --output "$audit_root/$label-$name" \
        --mode "$mode" --repeats 1 --timeout 3600
)

prepare_one before "$before"
prepare_one after "$after"
replay_one before smoke smoke
replay_one after smoke smoke
replay_one before audit-a audit
replay_one after audit-a audit
replay_one after audit-b audit
replay_one before audit-b audit
```

Every output directory must be new. Stop on errors; retain incomplete runs and
logs without reclassifying them as passes. A timeout does not prove an algorithmic
defect. After a timeout, check for surviving audit-owned compiler/benchmark
processes and stop those before another run; do not terminate unrelated work.
Report build wall time and maximum RSS exactly as the measurement tool
reports them, separately from sampling and native timing. Keep the target
directories intact until both binaries have finished their replays.

Each successful audit process must provide all 117 Criterion estimates and 81
counter observations. Results are under each output's `repeat-1/`: `run.json`,
`run.log`, `work.tsv` and `criterion/`. Match fixture hashes and interaction IDs
across revisions; retain both source revisions, binary hashes, effective profile settings,
toolchains and flags. The ABBA order reduces a simple order bias but is not a
statistical guarantee. The runner uses Criterion `--quick`; assess repeatability
and uncertainty before interpreting differences. Do not invent a speedup or a
release threshold from two successful quick runs.

Expected structural result: on a pan that rebuilds layer counts, the before
revision traverses GC bases again; the after revision contributes zero
`layer_gc_bases` in that path. Tree rebuilds may remain, and the separate GC
track still computes GC values when needed. Counters are not milliseconds.
Report effects per length, density and interaction rather than only a pooled
mean; retain regressions and inconclusive results too.

## Remaining Startup And Native Interaction Attribution

For native acceptance, use the packaged `--release` GUI at an explicitly named
candidate SHA, preferably extracted from its package. Retain artifact and binary
hashes, source SHA, build command, effective profile settings, features/flags,
toolchain and lockfile separately from the Criterion receipts. A local rebuild
must match the packaging recipe and be labelled rebuilt, not extracted-package
evidence. Do not substitute a `benchmark-support`, `gui-test-support`, semantic
snapshot-writer or profiler-specific build; label such diagnostics separately.
Do not silently use the older `ee9eb483` release recipe as current packaged-GUI
evidence. If no candidate package is available, report that gap explicitly.
Use the public PATZ1/TP73 fixtures from the
[existing runbook](../benches/README.md#gui-critical-operations), copied into an
isolated audit profile with bound project/report hashes. Do not use private
state, silently reopen a personal project or save over the supplied fixtures.

Set `GENTLE_GUI_STARTUP_TRACE` to a new path in an existing local evidence
directory, then exercise both empty startup and explicit `--project PATH`
startup. Open the intended DNA window and exit normally to retain the trace.
See [checkpoint semantics](../benches/README.md#startup-phase-checkpoints).
Record cold/warm conditions and external launch/visibly-confirmed-content times.
CPU markers alone do not prove that the right locus became visible.

Separate process launch, app initialization, project read/decode/install,
DNA dispatch, worker scheduling/lock/clone, hydration and first visible content.
Also exercise pan/zoom, resize, layer toggles and selection with recorded window
sizes. Keep harness waits, snapshot overhead and compositor/event delays distinct
from GENtle computation. Do not subtract optimized embedded-frame times from
debug Xvfb timings or sum nested/overlapping startup spans. Missing/lost spans
are unavailable measurements, not zero durations.

Prioritize attribution on the identified packaged candidate; a native GC
improvement/regression claim additionally needs before/after measurements with
the same effective profile and other build inputs. A local macOS
cross-check remains a separate follow-up before a general scheduling change;
Linux/Xvfb evidence alone does not establish a cross-platform root cause.

## Separate Same-SHA Native Profile Comparison

At the selected candidate SHA, compare the packaged-release GUI with a native
GUI built using `bench-audit`, not the Criterion harness. Prebuild both outside
timing, keeping features, toolchain and other build inputs identical. Use the
same host, public PATZ1 project/report, window sizes and tracing settings, with
equivalent isolated application settings and cache state. Repeat open and resize
operations, alternate run order, distinguish cold/warm conditions, and retain
individual samples, repetition counts, medians and spread with both binary hashes.

Record effective settings and any overrides, not only profile names. The release
recipe changed after the fixed GC pair; `bench-audit` retained its settings.
Neither a speedup nor a slowdown follows from that change without measurements.
This experiment measures the whole profile recipe, not LTO alone. Do not pool
it with the different-SHA GC experiment, derive a universal correction for older
timings, or substitute its outcome for packaged-native acceptance. It requires
no new installer-CI stage or change to production profiles.

## Requested Reply

Please return:

1. Exact revisions, binary/input/fixture hashes, host/toolchain/effective profiles/flags,
   retained evidence location and clean-worktree status.
2. Separate build wall time/RSS, smoke completeness and per-interaction CPU
   comparisons, including uncertainty, regressions and unchanged cases.
3. Confirmation or rejection of the GC work reduction; distinguish it from
   any measured responsiveness improvement. Do not promote historical
   2 Mbp/10,000-feature stress results into the agreed interactive envelope.
4. Packaged-native startup/interaction attribution and separately labelled
   same-SHA profile results, or explicitly why either is not yet available.
   Classify product defects, environment/build limits and harness
   gaps separately; keep missing evidence visible.
5. One recommended next implementation slice backed by the measured dominant
   cost, or an explicit request for the missing evidence needed to choose it.

No push, merge, release tag, changed acceptance threshold or automatic GUI
performance verdict is requested. Please retain the original evidence rather
than replacing previous audits.
