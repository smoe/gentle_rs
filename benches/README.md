# GENtle Performance Audit Targets

Developers verify that benchmark fixtures and scientific assertions still run
without collecting a statistical baseline:

```bash
cargo test -p gentle-benchmarks --bench gui_operations
cargo test -p gentle-benchmarks --bench specificity_finalization
```

Timed runs and release-facing comparisons belong to the external auditor.

## Audit profiles

GENtle provides two distinct optimized harness modes:

- Routine audit uses `--profile bench-audit`: `opt-level=3`, explicit stripping,
  thin LTO, 16 codegen units and `panic=unwind`. These effective settings are
  unchanged by the native release-profile changes in `e4c94bf3` and `ffe5c637`.
  Cargo forces unwind for benchmark targets; matching that setting throughout the
  audit dependency graph avoids a separate panic-strategy build.
- Release-profile harness runs use `cargo bench -p gentle-benchmarks` without a
  profile override and inherit the current release settings, including
  `lto="off"`, 16 codegen units and no stripping. Historical fat-LTO,
  one-codegen-unit results retain their original settings; do not relabel them
  or reuse them as baselines for the new recipe.

Cargo keeps compiled artifacts under `target/bench-audit/` and
`target/release/`, respectively. Criterion 0.8.2 does not infer the Cargo
profile and otherwise stores both modes under the shared `target/criterion/`
tree. Set `CRITERION_HOME` as shown below so measurements remain under the
matching profile directory. Do not pool these modes or use one as a code-change
baseline for the other: they differ in LTO and stripping. A controlled same-SHA
profile comparison is a separate diagnostic, not a reusable timing correction.
Retain the revision and effective settings, not just the profile name, which
spans different historical configurations.
The benchmarks live in the dedicated, non-published `gentle-benchmarks`
workspace package. It depends on GENtle as a library with default features
disabled and explicitly enables `desktop-gui` plus `benchmark-support`, because
the GUI target imports the feature-gated DNA-window modules it measures. Cargo
still does not prepare the application's GUI, CLI, MCP, documentation, or
publication-report binaries before sampling. The first library build and link
can remain substantial; subsequent runs reuse the profile cache. GENtle's
build fingerprint follows the current loose branch ref
and consults repository-wide `packed-refs` only when that loose ref is absent,
so Git maintenance in another worktree does not invalidate an otherwise
unchanged audit build.

Neither harness is the packaged GUI. Native release-facing acceptance uses the
packaged `--release` desktop binary, without the semantic snapshot writer;
retain extracted-package evidence separately from matching local builds.
Glen's repeated, same-SHA PATZ1 native open/resize comparison of release and
`bench-audit` GUI binaries follows the
[measurement-profile contract](../docs/dna_feature_rendering_latency_plan.md#measurement-profiles).
The headless ladder stays on `bench-audit` for historical continuity. No runtime
speed difference is assumed, and this audit adds no installer-CI benchmark stage.

## GUI-critical operations

The primary Criterion target exercises the real DNA-window constructor and
embedded egui render path with a feature-free 120 kbp control and the public
TP73 locus fixture:

```bash
CRITERION_HOME="$PWD/target/bench-audit/criterion" \
  cargo bench --profile bench-audit -p gentle-benchmarks \
  --bench gui_operations -- --quick --noplot
```

For the release-profile harness, omit `--profile` and change only the retained
result directory:

```bash
CRITERION_HOME="$PWD/target/release/criterion" \
  cargo bench -p gentle-benchmarks --bench gui_operations -- \
  --quick --noplot
```

It reports eager constructor, deferred UI-thread hydration, first-frame, and
steady-frame CPU costs. The hydration benchmark directly exercises the
`replace_loaded_sequence` boundary used after lazy background loading, while
keeping the background record clone outside the timed loop. Benchmark IDs
include sequence length, feature count, and the first 12 characters of the
input SHA-256. This headless target cannot measure native window creation, GPU
submission, compositor latency, input delivery, or whether interaction feels
responsive. It explicitly clears each returned egui texture delta because no
renderer exists to upload it. Those native concerns remain part of the external
GUI/Puffin acceptance described in `docs/testing.md`.

For the authentic public PATZ1 workload, first prepare a fresh offline project
with the exact CLI under review and then opt the benchmark into both the project
and its bound locus report:

```bash
audit_root="$(mktemp -d)"
python3 scripts/prepare_real_patz1_tutorial.py \
  --gentle-cli target/debug/gentle_cli \
  --output-dir "$audit_root/patz1"
GENTLE_GUI_BENCH_PATZ1_STATE="$audit_root/patz1/patz1.gentle.json" \
GENTLE_GUI_BENCH_PATZ1_REPORT="$audit_root/patz1/locus.report.json" \
CRITERION_HOME="$PWD/target/bench-audit/criterion" \
  cargo bench --profile bench-audit -p gentle-benchmarks \
  --bench gui_operations -- --quick --noplot
```

The PATZ1 setup fails closed unless it sees the expected 20,802-bp locus,
at least 75 loaded gene/transcript/exon features, 13 Ensembl transcript records
and 17 source-comparison records (13 Ensembl plus four RefSeq). It loads the
portable report through the same sequence-binding and presentation-cache path
as the GUI. First and steady frames are sampled at `820x520`, `1200x800`,
`1600x1000`, and `1920x1080`; separate cases measure the first frame after a
resize from `1200x800` to the other three sizes. The benchmark ID uses the
pinned fixture-manifest hash rather than the generated project bytes, whose
audit timestamp changes between equivalent preparations. Retain the exact
project/report hashes beside Criterion's results.

The two scroll groups, `dna_window_wheel_event_frame` and
`dna_window_first_frame_after_wheel`, dispatch an actual egui horizontal wheel
event over the painted DNA map, without modifiers. They use a 5-kbp viewport
starting at base 1,000, a six-point delta and fixed 60-Hz input timestamps.
Every iteration starts with a fresh warmed window outside the timed section,
so clamping, accumulated scrolling and window destruction cannot masquerade as
steady input performance. The second group excludes event dispatch and measures
only the following frame. These are CPU-side frames, not native delivery or
input-to-paint latency.

Before timing, both scroll directions must change the actual viewport without
changing its span or sequence length, emit paint shapes and leave sequence,
annotation, enzyme and methylation inputs byte-identical. The following frame
must retain the moved viewport. `GUI_SCROLL_OBSERVATION` JSON lines record these
checks against the full Git SHA and version. The pointer target comes from the
clipped egui interaction rectangle, not the potentially off-screen paint extent.
checks for every fixture and size, including the exact source revision and
content digest. A failed assertion invalidates the run. An omitted PATZ1 input
does not count as PATZ1 acceptance. Use `-- --test` instead of `-- --quick
--noplot` for deterministic smoke verification; retain its log and build receipt
separately from an auditor's timed run. No historical results are overwritten.

## DNA feature-density latency

The `.12` measurement foundation adds `dna_feature_latency` over the real
`MainAreaDna` presentation path. It separates length from feature density: all
nine combinations of 20,000 / 250,000 / 2,000,000 bp and 100 / 1,000 / 10,000
features, plus the exact proposed 250,000-bp / 5,000-feature interactive boundary.
The large cases are stress probes, not promised interactive limits. The boundary
is available for measurement; its inclusion is not owner approval or acceptance.

**Fixture provenance and recreation:**
`src/main_area_dna/latency_benchmark.rs::feature_density_fixture` generates exact
counts in memory from a fixed artificial DNA repeat, plus/minus joined mRNA and
CDS locations, exons, regulatory and repeat features. Half the features cluster
in the initial 5-kbp viewport; half are spread over the locus. There are no
biological claims, private inputs, downloaded annotations or RNG state. Running
the target recreates the fixtures and hashes the exact sequence/feature JSON.
The helper and view-only controls compile only with `benchmark-support`.

`feature_density_boundary_fixture` preserves that annotation geometry, and
inserts the artificial string `GAATTCGATCCCAGGATG` + 120 `GCC` codons + `TAA`
at position 100 and every 9,973 bp thereafter (zero-based). It configures one
synthetic `GAATTC` cutter (cut 1, overlap 4), both Dam/Dcm motif scans and visible
restriction/GC/ORF/methylation layers. This creates a nonempty derived workload
without changing the original nine fixtures or claiming real biology. The new
fixture digest additionally binds the enzyme definition and methylation mode.
One synthetic enzyme is not a full-catalog performance test: retain authentic
TP73/PATZ1 cases and their configured layer/catalog counts as separate evidence.

The 130 cases comprise constructor, hydration, first frame (deferred/loaded
tree), steady, one-base pan, zoom, mRNA toggle, selection, hover, and resize from
1200x800 to 820x520, 1600x1000 and 1920x1080. Each interaction also emits untimed
before/after cache counters, for 90 observations. Sequence/feature hashes must
remain unchanged. The synthetic run has no engine, locus-report hydration,
native window, actual X11 click, GPU, or compositor. Use `gui_operations` above
for engine-backed TP73/PATZ1 window cases and native acceptance for the rest.
Each observation also includes current before/after layer inventories: total
and viewport-eligible restriction groups, GC bins, ORFs and methylation sites,
their enabled flags, bin size and viewport. These are shared-toolbar counts,
not pixel visibility or readiness claims. Reading them does not populate a
cache or alter the measured counters. Missing/stale inventories fail the run;
empty or disabled derived layers invalidate the new boundary case.

Build **once**, separately from runtime, with no network access:

```bash
python3 scripts/dna_feature_latency.py prepare --profile dev \
  --output /tmp/gentle-dna-latency-dev-build
python3 scripts/dna_feature_latency.py run \
  --receipt /tmp/gentle-dna-latency-dev-build/build.json \
  --output /tmp/gentle-dna-latency-smoke --mode smoke
```

This runs two fresh processes by default (`--repeats 1` is available). Every
output directory must be new. Preparation uses `cargo test --locked --offline
--profile PROFILE -p gentle-benchmarks --bench dna_feature_latency --no-run`;
missing cached dependencies cause an explicit build failure, never a download.
The receipt retains build command/logs/time, source revision, dirty-diff and
untracked-build-input hashes, lockfile, toolchain, profile and binary SHA-256.
Only clean, frozen revisions are eligible for the auditor's baseline.

For Glen's timed run, prepare with `--profile bench-audit` into a different
directory and run with `--mode audit`. The latter executes the verified binary
directly with `--bench --quick --noplot`; it never invokes Cargo and rejects a
development-profile or dirty-source receipt. Baselines remain isolated in each run's
`criterion/` directory. Keep the full directory, not only the extracted means.
Two process repeats are evidence for the auditor to assess, not an automatic
stability or performance pass.

Each replay records `run.json` (including the runner's own hash), exact `run.log`
bytes/hash, `work.tsv` deltas and,
for an audit, raw Criterion artifacts and estimate hashes. The runner isolates
HOME/XDG/temp/cache paths, removes inherited `GENTLE_*` options, disables the
optional diagnostics pane/profiler, and records host/locale/thread settings.
It verifies all ten fixtures and all nine interactions per fixture; missing,
duplicate, wrong-revision or inconsistent-fixture records fail. Audit estimates
must match all 130 content-bound Criterion case IDs, not just a file count.
`work.tsv` includes viewport-eligible layer counts; full inventories remain in
`run.json`. Timeouts and
child failures retain their receipt and are never successes. Build or source
identity changes require a new build receipt. Do not edit the tree during
preparation. Native input-to-content latency is explicitly **not measured**.

The binary, build receipt and observations declare `density_boundary_v1`.
Older binaries/receipts without that field retain their `density_ladder_v1`
identity (117 cases, 81 observations). The runner can replay them, but records
`interactive_boundary_exercised=false`; they cannot certify the boundary by
interpolation. A receipt/binary workload mismatch fails before execution. Even
`interactive_boundary_exercised=true` is only a completed workload, not a timing
pass; the run mode distinguishes smoke from an audit.

Glen's [2026-09-27 current-main audit](../docs/glen_gui_latency_audit_20260927.md)
establishes the first clean `density_boundary_v1` CPU baseline at `0467867c`:
two complete `bench-audit` replays, 130/130 Criterion estimates and 90/90
counter observations. The 250-kbp/5,000-feature medians were 62.696 ms pan,
52.328 ms zoom and 51.209 ms mRNA toggle, with zero GC bases rescanned. This is
headless CPU evidence, not native-GUI or release acceptance.

The diagnostics pane and Puffin scopes are described in `docs/gui.md`;
`docs/dna_feature_rendering_latency_plan.md` retains the conditional optimization
sequence. Developers smoke-test; Glen provides timed and native acceptance.

For isolated tree-path attribution (not the Criterion ladder or native
acceptance), run the opt-in synthetic continuous-pan diagnostic:

```bash
cargo test --locked --offline --lib --no-default-features \
  --features desktop-gui,benchmark-support feature_tree_continuous_pan_diagnostic \
  -- --ignored --nocapture --test-threads=1
```

It reuses `feature_density_fixture` at 250 kbp with 1,000/5,000/10,000 features,
warms the tree, then records 50 one-base pan refreshes (without drawing). Raw
microsecond samples, debug-assertion mode, revision and work counters go to
stdout. Retain the exact source/diff, binary hash, command/toolchain and host
alongside the log; its revision alone does not bind an uncommitted experiment.
Construction, compilation and native display are outside those samples. There
is no timing assertion. `tree_viewport_updates` in current diagnostics separates
recounts from full `tree_builds`; older logs without it are unavailable on that
axis. These observations cannot certify a release budget or be pooled with
optimized ladder measurements.

### Startup phase checkpoints

The desktop binary can retain a small CPU timeline without `gui-test-support`,
Puffin, screenshots or per-frame disk writes. On a locally built, identified
binary, use an explicit **new file** in an existing local output directory:

```bash
GENTLE_GUI_STARTUP_TRACE=/tmp/gentle-startup-empty-01.json \
  /absolute/path/to/gentle
GENTLE_GUI_STARTUP_TRACE=/tmp/gentle-startup-project-01.json \
  /absolute/path/to/gentle --project /absolute/path/to/public-project.gentle.json
```

Open the intended DNA window, then exit normally. The destination is written
only after the native loop returns; an existing file is never overwritten.
Missing parent directories or write failures are reported on stderr without
changing project behavior. Forced termination can leave no file: absence is
not evidence of fast startup or a successful run. There is no live flush button
or hang-capture promise. Unset/empty `GENTLE_GUI_STARTUP_TRACE` is inert.

`gentle.gui_startup_trace.v1` binds the compiled source revision and reports
monotonic microseconds since Rust `main` entry. It excludes OS process loading
before that point. Each event has a fixed phase/kind and process-local ordinal
subject/span, never a sequence name, bases, path, credential or error text.
Matched span begin/end rows distinguish completion, failure and interruption.
The record is capped at 512 events; contention/cap losses increment
`dropped_events`. Any loss or unmatched span makes the affected duration
unavailable, not zero. Do not sum nested or overlapping worker spans.

Use these boundaries for attribution:

- `native_run` begin to `app_initialize` begin separates pre-constructor native
  setup from app defaults, help preparation, configuration and credential load.
  `native_run` ends on exit, not when startup finishes.
- `help_preparation` contains `help_manuals`, `help_shell_reference`,
  `help_tutorial_discovery` and `help_tutorial_selected_load`. Reloads can repeat
  these subspans; do not sum them with their enclosing span.
- `help_open` and `help_tutorial_open` measure the first actual help handlers,
  including payload preparation and search refresh, not rendering/presentation.
  A focus-only return on an already open manual does not consume `help_open`.
  `help_tutorial_menu_discovery` measures the first guided-walkthrough discovery
  called by File > Open tutorial project, not the rest of that menu's work.
  First-use spans are recorded once per process; later calls still execute.
  `help_tutorial_switch` measures the first changed tutorial selection while
  help is already open (including the Topic picker), excluding subsequent
  search refresh and painting. Startup's initial selection does not consume it.
- `project_read_decode` and `project_install` distinguish file/parsing work from
  installation and presentation-state reset. Failed loads have no install span.
- `root_first_frame` may be a splash; `root_workspace_frame` excludes it. Both
  precede the optional semantic snapshot writer, not GPU/compositor completion.
- `dna_open_dispatch` measures the common application dispatcher, including a
  focus/no-match return; completion does not imply a window opened. Constructor
  and worker subjects identify subsequent work without persistent sequence IDs.
- `dna_worker_scheduled`, `dna_engine_read_lock`, `dna_sequence_clone`,
  `dna_worker_result` and `dna_hydrate` separate scheduling, lock wait, copying
  and foreground hydration. Direct/eager routes legitimately lack worker spans.
- `dna_native_content_frame` / `dna_embedded_content_frame` occur once after
  loaded DNA rendering returns, never for loading/error placeholders. Deferred
  feature trees or evidence jobs may still be pending. These are CPU markers,
  **not subject-verified, fully ready or compositor-presented content**.

Glen should retain exact binary/lockfile/project/report hashes, toolchain,
profile, flags, isolated profile/cache setup and cold/warm conditions alongside
each trace. Keep packaged-release native macOS and Linux runs as separate
evidence classes; record process-launch and visibly confirmed content with the
external native audit. This file alone cannot certify an interaction budget,
and enabling tracing adds some overhead. The existing Criterion replay remains
unchanged and clears inherited `GENTLE_*` options.

The additive `help_image_work` object retains fixed, path-free **whole-session**
totals: `svg_references`, `cache_hits`, `preparation_failures` (cache-directory
creation), `rasterization_attempts`, `rasterization_completed` (including errors),
`rasterization_failures` and `rasterization_us`. References are visits, not unique
images. Rasterization time surrounds the existing file-to-PNG call, including
file I/O, font discovery, parsing, rendering and writing; it does not isolate
font cost. These totals include startup, reloads, help and menu work and overlap
their spans. They consume no events from the 512-event budget. Contention is
nonblocking and increments `dropped_help_image_observations`; overflow sets
`saturated`. An unfinished attempt, loss or saturation prevents a complete-total
claim. The reader retains partial counters but marks the complete duration
unavailable. Counter loss does not invalidate intact event spans, nor event loss
intact counters. Older v1 traces without this extension remain readable and show
image work as unavailable, not zero. No paths, image names or error text are
recorded; loading, fallback and cache behavior remain unchanged.

#### Read a saved startup trace

Summarize one existing trace without launching GENtle, Cargo, a GUI or a timed
benchmark. From the repository root:

```bash
python3 scripts/gui_startup_summary.py /tmp/gentle-startup-project-01.json \
  --output-dir /tmp/gentle-startup-project-01-summary
python3 -m unittest scripts.test_gui_startup_summary -v
```

The output directory must be new. It contains the untouched `trace.json`, a
machine-readable `summary.json` (`gentle.gui_startup_summary.v1`) and
`summary.md`. The report hashes the exact input bytes and reader script; LF and
CRLF inputs remain distinct byte evidence. This is not a redaction tool: the
raw input is retained verbatim, including any extra fields a producer supplied.
Review it before sharing. Invalid/oversized input is rejected before output is
created; prior output is never replaced. The input limit is 1 MiB / 512 events.

Read the CPU markers first, then the inclusive spans for each process-local
subject. Project decoding/installation, lock acquisition, cloning and hydration
remain separate. Same-subject boundaries show worker scheduling to lock attempt,
worker result to hydration, and construction/hydration to native or embedded
content markers. There is no inferred link from an open dispatcher to a DNA
subject, and no matching of subject ordinals across runs. Eager windows may
legitimately lack worker/hydration measurements.

Missing, ambiguous, unsuccessfully completed or inverted boundary endpoints
produce `null` / `unavailable`, never zero. Any recorded loss suppresses all
derived durations because losses cannot be localized; the observed timestamps
remain available for inspection. Failed/interrupted spans retain their outcome
and, when both endpoints exist without loss, their elapsed duration. Do not sum
inclusive/nested/concurrent spans. `native_run` includes the entire session until
exit, not just startup. `no_recorded_loss` is trace integrity, not successful
startup; the performance verdict always remains `not_assessed`.

The tool does not combine runs, pick a bottleneck automatically or verify binary,
profile or subject identity. Glen still binds the external evidence described
above and owns attribution and native acceptance. CI exercises synthetic parser
and real CLI/file-boundary regressions only, not GUI timings.

## Specificity finalization

Run the deterministic primer-specificity benchmark with:

```bash
CRITERION_HOME="$PWD/target/bench-audit/criterion" \
  cargo bench --profile bench-audit -p gentle-benchmarks \
  --bench specificity_finalization -- --quick --noplot
```

`specificity_finalization.rs` generates non-biological complete-query BLAST
rows in memory at 100, 1,000, and 6,600 total HSPs. Every synthetic subject has
one inward-facing primer pair; one product is declared intended and all others
are exact off-target products. Fixture construction and scientific accounting
assertions run outside the timed loops. The first 12 characters of the SHA-256
over both generated TSV streams are included in each Criterion benchmark ID.

The benchmark executes no BLAST process, downloads no reference, and uses no
private data. Its tiny temporary genome catalog exists only to make catalog
resolution independent of locally prepared references. The generator in the
benchmark source is the deterministic recreation procedure and the source
revision printed at startup binds a run to its GENtle checkout.

Compare performance only when runner, toolchain, build profile, fixture-hash
ID, and GENtle revision all match. Baseline labels and retained audit metadata
must name the profile. Shared CI runners are suitable for compilation or smoke
runs, not strict regression thresholds.

## External-auditor ownership

GENtle developers provide and smoke-test these targets but do not assign their
own performance verdict. Glen is the named external auditor for the current
release assessment. He should run an exact clean revision on a stable host and
retain:

```bash
git rev-parse HEAD
rustc -Vv
cargo -V
uname -a
CRITERION_HOME="$PWD/target/bench-audit/criterion" \
  /usr/bin/time -v cargo bench --profile bench-audit \
  -p gentle-benchmarks --bench gui_operations -- \
  --quick --noplot --save-baseline bench-audit_HOST_TOOLCHAIN_REVISION
```

Use `/usr/bin/time -l` instead of `-v` on macOS. Keep the console log and the
corresponding `target/bench-audit/criterion/` subtree with the audit record.
For a release-profile harness baseline, use `cargo bench -p gentle-benchmarks`
without `--profile`, set
`CRITERION_HOME="$PWD/target/release/criterion"`, and prefix the baseline label
with `release-profile` (do not rename historical evidence). Compare baseline
repeats only against a baseline produced by the same
profile, host, toolchain, fixture hash, and GENtle revision. The auditor decides
whether a difference is meaningful after also performing the real GUI
interaction check; GENtle does not turn a noisy shared runner result into a
release verdict.
