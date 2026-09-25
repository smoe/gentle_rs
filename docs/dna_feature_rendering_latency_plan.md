# DNA Feature Rendering Latency Plan (`.12`)

Planned: 2026-09-20. This is the first `.12` priority. Implementation provides
instrumentation, deterministic fixtures and bounded changes; the named external
auditor retains baselines and owns the verdict
([DEC-039](decisions.md#dec-039-external-auditor-owns-the-performance-verdict)).

Aim: make DNA-feature presentation in the sequence viewer respond immediately on
real annotated loci - opening a locus window, panning/zooming, resizing,
toggling feature layers, selecting and hovering. This plan changes *when*
features appear, never *what* is shown: identical features, coordinates, labels,
strands and scientific outputs.

Implementation update: S0's density tools landed in `5893aa35` and are
smoke-tested. Opt-in startup phase checkpoints (`d24c3fa7`) cover app
initialization, project loading and first root/DNA CPU frames; they do not
confirm native presentation. The exact proposed 250 kbp / 5,000-feature
boundary now includes explicit derived layers and content-bound workload checks.
Owner scope choices, the macOS cross-check and Glen's timed/native audit remain
pending. A bounded S2 simplification (`ee9eb483`) counts GC bins arithmetically
instead of computing GC values just to discard them. It preserves half-open
counts and does not claim a measured GUI speedup or completion of S2; broader
runtime slices remain pending.
B0 separates build feedback from runtime performance. The consolidated
hypothesis ledger and S0 section below retain both updates. This plan does not
change the `.11` candidate; counts of work are not timing evidence or performance
acceptance.

## Review Reconciliation (2026-09-21)

The owner supplied Claude's read-only review; these revisions were checked
against `4e858f3a66cec54f5bdc3eec76fb2f18683773e9`, with the same latency
paths at `70a3b038`. They are Codex's revised proposal, not a second Claude
review or owner approval of the open scope choices.

| Original plan | Claude feedback | Revised minimal plan |
| --- | --- | --- |
| Measure startup, but only S1-S5 own fixes | The largest observed delay has no implementation owner | Add S6 for startup/first usable content; run it first if attribution confirms the dominant product delay |
| Prebuilt runner; ambiguous workspace non-goal | Build feedback needs its own scope | Reuse the landed build-once runner; B0 permits evidence-led build-boundary work, independently of runtime hotspots, not an unconditional root-crate split |
| Workload ceiling and circular scope remain open | Decide these before the audit | Record explicit owner choices below; add exact boundary fixtures before a pass/fail audit |
| First-content budget lacks a dependency | Whole-sequence recomputation may exhaust it | Make S3/freshness-correct reuse a conditional prerequisite; do not extrapolate the whole PATZ1 constructor as recomputation time |
| Native evidence is Linux/Xvfb only | Cross-check on macOS before S5 | Require a local packaged-profile macOS trace and retain platform-specific conclusions |
| Three overlapping priority lists | Keep one implementation sequence | Roadmap owns release scope, this plan owns steps, and the historical acceptance report links here |

## Build Profile Context (Two Steps)

The native release recipe changed twice while this plan was being written, and
[DEC-039](decisions.md#dec-039-external-auditor-owns-the-performance-verdict)
keeps both steps as distinct evidence classes:

- `e4c94bf3` dropped the explicit release block,
  so native releases took Cargo's defaults instead of fat LTO, one codegen unit,
  `panic=abort` and stripping. Its default `lto=false` still permitted local
  thin LTO. `release-fast` retained its explicit panic/stripping settings.
- `ffe5c637` then set `[profile.release] lto = "off"`, disabling even
  within-crate LTO. That is the current tree state.

`bench-audit` inherits `opt-level=3` from `release` and explicitly sets thin LTO,
16 codegen units, `panic=unwind` and stripping. Making stripping explicit
preserved its previously inherited value, so the
audit profile's effective settings survived both steps unchanged. Three
consequences bind this plan:

- The PATZ1 means quoted below keep their original profile identity. Per DEC-039
  and the [audit-profile notes](../benches/README.md#audit-profiles), record the
  effective settings, not only the `release` profile name, which now spans three
  historical configurations.
- Today's native `release` binary and the thin-LTO audit binary use different
  optimization recipes. Neither the direction nor the perceptibility of a
  runtime difference is established. Record effective settings with every
  measurement and do not pool traces across either boundary.
- B0's build-cost premise changed: the native release path no longer performs a
  fat-LTO link. Re-baseline build cost on current `main` before proposing any
  dependency-boundary change.

These are build-configuration changes with their own acceptance evidence, still
under release validation. Their runtime effects remain unmeasured; neither
establishes a speedup or regression of the runtime hypotheses below.

## Audit Admission Decisions

These two owner decisions were requested on 2026-09-21 and remain pending.
Recommendations are not accepted performance promises:

- **Interactive envelope:** recommend linear loci up to 250 kbp and 5,000
  loaded features; retain 2 Mbp and 10,000 features as stress cases. The exact
  250 kbp / 5,000-feature boundary is now available with the same clustered
  annotation and deferred/loaded-tree modes, plus explicit synthetic restriction,
  GC, ORF and methylation layers and their counts. This closes the fixture gap,
  not the scope decision or timing gate. The original nine cases remain
  comparable; their results alone cannot certify this boundary by interpolation.
- **Circular scope:** recommend linear optimization first, with circular-map
  correctness and native interaction regression checks still required. If the
  owner instead requires circular optimization, agree its fixture envelope and
  budgets before running a dedicated pass/fail audit; linear results do not
  establish circular performance.

Confirm these choices before Glen's release-target audit. Preparatory traces
may still diagnose a problem, but cannot declare an undecided workload accepted.
Glen proposes host-bound numerical budgets from the baseline; the owner agrees
the release target before optimization. Do not relax it after seeing results
without an explicit recorded scope decision. Stress cases must remain correct,
inspectable and interruptible where work is cancellable, but are not silently
held to the interactive timing table.

### Measurement Profiles

For native acceptance, **measure the actual packaged GUI and declare its profile**,
not the `bench-audit` harness. Since 2026-09-25, internal installers use `dev`;
final releases use `release`. The historical "release-like" label still means
the optimized release profile, not an unoptimized interim package. Keep those
evidence classes and their agreed budgets separate; neither certifies the
other and slower development binaries do not silently relax latency gates.
Prefer the extracted candidate package and retain its
artifact and binary hashes, source SHA, lockfile, toolchain, features and
effective profile settings. A locally rebuilt counterpart must match the
packaging recipe and be identified as such; it is not extracted-package evidence.
Do not substitute `benchmark-support`, `gui-test-support`, a semantic snapshot
writer or a profiler-specific build for the shipped binary. Keep instrumented
diagnostics separately labelled and retain external input-to-content evidence.

Keep `bench-audit` for continuity of the headless Criterion ladder. The optional
release-profile Criterion harness also remains a harness: its target, features
and event/rendering environment differ from the packaged desktop application.
Neither harness can satisfy native acceptance, even when its CPU timings pass.

Glen should additionally measure a controlled **same-SHA native profile pair**:
the packaged GUI (explicitly `dev` or `release`) and a native GUI built with `bench-audit`, keeping
features and all other build inputs identical. Prebuild both outside timing;
use the same host, toolchain, public PATZ1 project/report, window sizes, tracing
settings and equivalent isolated application settings/cache state. Repeat open
and resize operations for each profile, alternate run order, separate cold and
warm conditions, and retain individual samples, repetition counts, medians and
spread alongside both binary hashes. No semantic snapshot writer is used.

This comparison measures the whole profile recipe, not LTO alone. Report any
observed difference for that workload and host; do not derive a universal offset
or correct historical timings with it. Ordinary code-change comparisons still
hold the profile constant. Native release acceptance must meet its own agreed
budgets regardless of the audit ladder's outcome.

This is an external audit requirement, **not an additional installer-CI build
or benchmark stage**. It does not explain past compiler/runner terminations;
build reliability and runtime performance require separate evidence.

## Execution Order And Ownership

The slice numbers identify work, not a mandatory S1-through-S6 order. Codex
provides instrumentation, smoke checks and bounded fixes; Glen owns timed and
native acceptance; the owner decides release scope and any explicit deferral.

1. Resolve audit admission, then finish S0 attribution. Reuse the prebuilt
   runner; pursue B0 only if build cost still impedes those audits.
2. Compare startup/first-content and ongoing-interaction costs. Prioritize S6
   if the large startup observation is a product delay; otherwise document its
   harness/environment contribution and select the largest confirmed H1-H4 path.
3. Apply S1-S3 according to measured contribution, not their numbering. S5 needs
   platform-specific native evidence; S4 needs a remaining drawing-budget miss
   within the agreed envelope. Neither is an automatic follow-up.
4. Re-run correctness and compatible before/after measurements after each
   change, then obtain one exact-candidate verdict. Preserve `.11` evidence
   separately; this plan does not reopen or silently waive its release gates.

## Why this leads `.12`

[GUI usability acceptance 2026-09-20](gui_usability_acceptance_20260920.md)
accepted the fixed TSS workflow but explicitly did **not** accept general GUI
usability or performance. The evidence already separates two things:

- Both isolated acceptance runs spent about 32 s in application startup plus
  opening the first locus, while the remaining operations usually took about
  1-3 s each.
- Native PATZ1 resize-to-confirmed-content took 409.6 / 63.1 / 229.4 / 163.7 /
  163.3 ms (run A) and 425.7 / 62.3 / 187.3 / 203.7 / 183.4 ms (run B), with the
  `gui-test-support` snapshot writer included.
- The optimized `bench-audit` Criterion means for the same authentic PATZ1
  project, taken before the release-profile steps above and valid only for their
  own revision and effective settings, were 27.852 ms eager DNA-window
  construction, 32.742 ms deferred
  UI-thread hydration, 28.1-29.5 ms first embedded frame at all four viewport
  sizes, 0.786-0.793 ms steady embedded frame, and 3.83-3.94 ms first frame
  after a resize.

Optimized steady painting of this small 20,802-bp PATZ1 locus with 75+ loaded
features is inexpensive in the embedded benchmark. It neither explains the
debug harness's startup observation nor proves all painting fast.
Construction/hydration, viewport-dependent work, native event-to-content delay
and process/project startup remain distinct
hypotheses. The roughly 32 seconds includes fixed waits and capture overhead;
do not call it 32 seconds of application startup or subtract the optimized
Criterion timings from it. These historical numbers are not density-ladder
acceptance; feature-count scaling still awaits the auditor.

## Scope

In scope: process/project startup through first usable root and DNA windows,
DNA-map feature presentation, the feature tree and layer-visibility panel,
per-window hydration, and repaint/resize scheduling. The proposed linear-first
optimization scope and circular regression floor await the owner decisions
above. B0 is a separate, bounded build-feedback workstream, not GUI redesign.

Out of scope for `.12`: visual redesign, new feature classes, engine-side
annotation work, sequence-text panel work beyond what the feature path forces,
and any change to displayed biology. Retain
[DEC-026](decisions.md#dec-026-long-gui-jobs-use-optimistic-engine-snapshots)
and
[DEC-048](decisions.md#dec-048-command-submission-and-observation-never-wait-for-execution);
this plan adds no new execution authority.

## Proposed interaction budgets

These proposals assume the recommended envelope, pending owner confirmation.
Glen establishes reproducible numerical targets before optimization; each is
bound to a host, toolchain, profile, fixture hash and GENtle revision. First
usable content means subject-correct selected feature content, not an empty
window or progress placeholder. Record pending derived work separately.

| Interaction | Subject | Proposed budget |
| --- | --- | ---: |
| First confirmed feature content after an open request in a usable application | <= 250 kbp, <= 5,000 features | <= 500 ms, no single frame > 100 ms; conditional S3 dependency below |
| Pan / zoom step | same | p95 <= 16.7 ms CPU prepare+paint |
| Live resize | 1920x1080 | every intermediate frame <= 33 ms; confirmed content <= 150 ms |
| Feature layer toggle (CDS/repeat/array/TFBS) | same | <= 100 ms to confirmed content |
| Selection / hover | same | <= 16.7 ms, no model rebuild |
| Tab switch / multiwindow focus | same | <= 100 ms to confirmed content |
| Process start to usable root window, and separately first subject-correct locus | empty profile and reference project, cold/warm recorded separately | S6 owns attribution and any fix; agree budgets after baseline and before optimization/acceptance |

The first-content target depends on S3 if H3/H4 leaves insufficient frame or
wall-time budget. Freshness-correct reuse or deferred work is a possible means,
not a presumed fix. The 27.852 ms figure includes the entire PATZ1 constructor;
its recomputation share is unknown, so multiplying it by a length ratio is not
a measured prediction. Test the selected boundary instead. The 2 Mbp stress
case has no 500 ms promise. Likewise, startup cannot receive a green verdict
while its budgets remain unspecified; an explicit owner deferral is not a pass.

## Hypotheses And Evidence Ledger

The evidence below separates verified work from unmeasured timing contributions.
Before optimization, record measurements across the fixture ladder for H1-H4
and native traces for H5-H6. A rejected hypothesis stays documented with its
evidence. Counters measure work, not time: a hit does not establish that a cache
is fast, and a miss does not establish that it dominates a frame.

**H1 - Full relayout on every changed rectangle.**
In `src/render_dna_linear.rs`, `RenderDnaLinear::render` recomputes
`layout_features` whenever the drawing rectangle changed or the layout dirty
flag is set. During a live drag-resize that is one full relayout per frame. The
relayout allocates a `Seed` per visible feature, including two `String`s (label,
upper-cased kind), then sorts and lane-allocates them. The rectangle/dirty guard
already exists; another signature containing the rectangle would still miss on
every resize. Layout/index counters and `draw_features` scopes separate these
costs. The renderer regression confirms resize increments layouts without
rebuilding the interval index; the density-dependent runtime share is unmeasured.

**H2 - Viewport-keyed model rebuilds.**
`FeatureTreeCacheKey` (`src/main_area_dna/feature_tree_ui.rs`) and
`LayerVisibilityCacheKey` (`src/main_area_dna.rs`) both contain the viewport,
and `active_linear_viewport_range` returns `Some` for every linear sequence. A
one-base pan invalidates both caches. At `ad0338a7`, layer-count misses also
computed GC values across the full sequence just to count overlapping bins:
O(all features + bases), not merely O(visible features). The bounded S2 repair
uses length/bin/viewport arithmetic for that count; GC values used for drawing
and scientific output are unchanged. The pan regression still requires tree
builds 1 -> 2, but now requires zero `layer_gc_bases` and exact counts before
and after crossing a bin boundary, without a project mutation. Historical
1,200 -> 2,400 GC-base counters belong to the pre-fix implementation.
Whole-view timing and the remaining tree/feature-count work still need Glen's
audit. Static grouping and viewport counts must be separated carefully:
the tree uses exact exon-piece overlap, whereas the existing feature interval
index uses transcript bounding spans, including introns. They are not
interchangeable visibility oracles.

**H3 - One-shot hydration inside a single UI frame.**
`WindowDna::poll_deferred_load` (`src/window_dna.rs`) calls
`replace_loaded_sequence` -> `replace_active_dna` (`src/main_area_dna.rs`),
which invalidates every derived cache, reconciles the viewport, updates the map
renderer and refreshes the construct-reasoning overlay in one frame (32.742 ms
on PATZ1). Background lock/clone, foreground hydration, viewport reconciliation,
map update and overlay refresh now have separate scopes; the background thread
still only clones the record. Their dense-locus and native costs await audit.

**H4 - Whole-sequence recompute at window construction.**
`MainAreaDna::new` calls `ensure_local_restriction_site_catalog_current`, which
runs `DNAsequence::update_computed_features` (`src/dna_sequence.rs`): restriction
sites, restriction groups, ORFs, methylation sites and GC content over the whole
sequence, regardless of what the first frame shows. Separate restriction, ORF,
methylation and GC scopes now attribute this work. Its contribution to the
27.852 ms eager construction on 20 kbp must be measured at 250 kbp and 2 Mbp.
Existing computed features alone are not proof that those values are current;
any reuse needs sequence, enzyme-catalog and parameter identity.

**H5 - The native gap.**
Debug native resize and optimized embedded-frame timings differ markedly, but
are not comparable measurements. Attribution across event delivery, viewport
synchronization, repaint scheduling, GPU upload, compositor work and the
`gui-test-support` snapshot writer does not exist yet. Cross-check on native
macOS before generalizing Xvfb/Openbox observations into a product fix. Direct
display setters in the headless harness measure presentation work, not X11
input, user navigation, window focus or compositor behavior. Scheduling changes
remain conditional on a native trace without the snapshot writer.

**H6 - Startup and first usable content.**
The harness combines process startup, project/window opening, fixed waits and
capture work. Use the bounded CPU phase markers around `src/bin/gentle.rs`,
app initialization, project load, first root/workspace and DNA frames alongside
independent native presentation timestamps and the existing profiler scopes.
Use an empty clean profile and the same public PATZ1/TP73 projects. A small
constructor benchmark cannot exonerate or explain process startup.

## Slices

### S0 - Measurement Foundation (Tools Implemented; Audit Pending)

**Implemented and smoke-tested; reuse these tools:**

- `dna_feature_latency` benchmarks the real `MainAreaDna` constructor, hydration,
  first frame with tree deferred/loaded, steady frame, one-base pan, zoom, mRNA
  layer toggle, feature selection, hover and three resize transitions. It emits
  130 cases and 90 counter observations across all nine combinations of
  20 kbp / 250 kbp / 2 Mbp and 100 / 1,000 / 10,000 features, plus the exact
  250 kbp / 5,000-feature boundary. Length and count vary independently; the
  largest case is a stress probe, not an interactive performance promise.
- Each fixture binds exact sequence and feature bytes. The generator uses
  structured, overlapping plus/minus transcripts, CDS, exons, regulatory and
  repeat features, with half clustered in the first 5 kbp. It uses no private
  annotation, random state, network, prepared genome or binary fixture blob.
- The additional boundary fixture has explicit synthetic cut-site/ORF/methylation
  inputs and enabled derived layers. Untimed inventories reuse current toolbar
  counts; absent/stale caches are unavailable, not zero. The old nine cases and
  their identities are unchanged. One synthetic enzyme does not represent the
  full restriction catalog; authentic workloads remain required.
- `scripts/dna_feature_latency.py prepare` builds offline once and binds the
  executable, source/diff, toolchain, profile and lockfile. `run` rechecks that
  binary and executes it directly with isolated profile/cache/temp directories.
  Build time is separate; failure, timeout, raw logs and work tables are retained.
  Timed audits reject dirty-source or development-profile receipts.
  Workload identity distinguishes legacy nine-fixture receipts from boundary
  runs; exact Criterion case IDs must agree with the observation content hashes.
- `GENTLE_DNA_CACHE_DIAGNOSTICS=1` exposes an opt-in DNA-viewer diagnostics pane.
  It is observation-only, reads renderer counters without waiting on its lock,
  contains no sequence/feature names, and writes nothing to project state.
- Additional Puffin scopes distinguish feature painting, interval indexing,
  construction computations and hydration substeps. Existing tree build/render,
  layer-count and display-sync scopes are retained.
- `GENTLE_GUI_STARTUP_TRACE` retains bounded, process-local CPU checkpoints
  without biological identifiers or per-frame I/O, writing a new JSON file
  only on native-loop exit. The [startup runbook](../benches/README.md#startup-phase-checkpoints)
  distinguishes splash, project decoding/installation, worker lock/clone,
  hydration and content-frame return; losses/failures remain explicit. Neither
  these markers nor a successful return prove visible or fully ready content.

See [the benchmark runbook](../benches/README.md#dna-feature-density-latency)
for exact generation, build and prebuilt replay commands. Keep the existing
`gui_operations` TP73/PATZ1 workload: the synthetic subtree harness does not
replace engine-backed window/report hydration or real annotations.

**Pending audit and exit criteria (Glen):**

1. Resolve the owner admission decisions above; if a different envelope is
   selected, add its exact boundary before acceptance. Freeze a clean SHA and run two
   prebuilt `bench-audit` repeats on a stable host. Retain the per-fixture,
   per-interaction tables, raw distributions and
   counter deltas; assess repeatability rather than treating two successful
   processes as a timing pass. Bind every result to source revision, profile,
   toolchain, fixture and binary hashes.
2. Capture native input-to-content traces at 820x520, 1200x800, 1600x1000 and
   1920x1080. Separate process startup, open, clone/hydration, first paint,
   event delivery, repaint and compositor delays. Use the actual packaged GUI
   under [Measurement Profiles](#measurement-profiles), without the semantic
   snapshot writer; capture locally, not on per-frame
   network storage. Retain project/report hashes, tool/environment identity and
   the binary's effective profile settings, which the profile name alone no
   longer fixes.
   Bind the CPU startup trace to external native observations; S6 owns any
   resulting product fixes, not a faster harness checkpoint.
3. Before a general S5 scheduling change, obtain a short local macOS native
   resize/open trace on the same source and public fixture, with its own binary,
   packaged profile, window sizes and environment recorded and no
   semantic snapshot writer. This is a diagnostic cross-check, not a macOS performance verdict.
   Do not pool Linux and macOS timings or attribute differences to an OS from
   unlike profiles. Windows native acceptance remains separate for Windows claims.
4. Record the confirmed/rejected/unresolved timing hypotheses and agree
   numerical budgets with the owner before optimization. Only then authorize
   the corresponding S1-S6 changes against their evidence requirements below.
   Record workload and circular scope; keep stress outcomes and release-target
   outcomes separate, without silently changing the admitted density ceiling.

The headless runner does not establish native acceptance, even when every smoke
case passes. Retain Glen's repeated native profile-pair comparison separately
from the ladder and native acceptance. S0's implementation is not to be
repeated; its measurement gate remains open.

### B0 - Build feedback (independent of runtime hypotheses)

The dedicated benchmark crate still depends on the root library with
`desktop-gui`; a cold or changed-source build can be expensive. However,
`scripts/dna_feature_latency.py prepare` already builds once and `run` replays
the hash-bound binary without Cargo. Use that path first. It solves repeat-run
rebuilds, not initial compilation or iteration after a source edit.

If build cost still blocks work, first re-baseline it on current `main`:
`e4c94bf3` removed fat LTO and `ffe5c637` disabled remaining local LTO in the
native release path, so earlier build timings no longer describe that recipe.
Then retain cold and incremental Cargo timings, compile versus LTO/link phases, wall time and
peak RSS. Evaluate the smallest
evidence-backed build/dependency-boundary change; a narrowly reviewed crate
extraction is eligible without proving a runtime hotspot. Follow DEC-006/007,
preserve production profiles and benchmark semantics, and remeasure build cost
at named revisions. Root source size alone does not prove which extraction
helps, and no broad engine/GUI split is authorized by this plan.

This is distinct from a single-root **GUI workspace** (replacing native child
windows with one application workspace), which remains a conditional UX/runtime
change. Neither faster linking nor a different audit profile establishes a
runtime speedup; build comparisons have their own acceptance evidence.

### S1 - Reduce Measured Layout Allocations And Dependencies (H1)

`RenderDnaLinear::render` already checks the rectangle and a layout dirty flag.
A new signature with the same rectangle cannot avoid work during real resizing.
First remove a measured redundant dependency or allocation, preserving dirty
invalidation and hit areas. Keep selection and hover from unnecessarily dirtying
layout. Reuse seed/lane buffers only if measured allocation cost warrants it;
avoid per-feature label/kind allocation where comparisons can borrow data.

Evidence required first: H1 confirmed on the ladder, with the relayout share of
a resize frame reported.

### S2 - Decouple tree and layer counts from the viewport (H2)

Bounded count-only repair: `GcContents::region_count_for_viewport` replaces the
full-sequence scan previously used solely for the GC layer's bin count. Oracle
tests compare it with materialized bins, including partial/zero-sized bins,
half-open boundaries and density-ladder lengths. This removes unnecessary work
without a new cache, worker, density policy or scientific change. It is not a
timed acceptance result; the broader S2 work below remains conditional.

Split both models into a viewport-independent part (grouping, filtering,
ordering, labels - keyed by feature generation and display revision) and a cheap
viewport pass that only marks visibility and recounts. Reuse an index only with
exact exon-piece and half-open overlap semantics, not transcript bounding spans
that include introns. Collapsed/offscreen models may defer work but must expose
correct counts when queried. Retain the arithmetic GC-bin count separately
from GC-value computation. Do not simply remove the viewport from existing keys.

Evidence required first: H2 confirmed, including the rebuild cost at 10^3 and
10^4 features during a continuous pan.

### S3 - Take one-shot work off the first frames (H3/H4)

- Move derivable presentation work out of the hydration frame: precompute what
  the background load can already produce, and time-slice the remainder across
  frames with a visible, bounded progress state rather than one long frame.
- Reuse restriction sites, ORFs, methylation and GC content only after checking
  sequence, enzyme-catalog and parameter identity; presence alone is not
  freshness. Where recomputation is needed, make it explicit, cancellable and
  off the first paint using existing immutable snapshots and owner checks.
  Never install stale/cancelled results or change scientific output.
- Keep the existing deferred feature-tree behaviour
  (`FEATURE_TREE_DEFERRED_AUTO_LOAD_MAX_FEATURES = 300`) and re-examine the
  threshold only with ladder evidence.

Evidence required first: H3/H4 confirmed, with construction and hydration split
by sub-step at 250 kbp and 2 Mbp.

### S4 - Density-aware feature drawing (only if S1-S3 leave a gap)

Here "gap" means an attributed drawing-budget miss inside the owner-approved
interactive envelope, not merely a slow 10,000-feature stress run.

The renderer already has explicit detail thresholds
(`FEATURE_LABEL_MAX_BP_PER_PX`, `RE_SITE_MAX_BP_PER_PX`, ORF and methylation
limits). If drawing still dominates at high feature density, extend that
mechanism - merged low-width features, label lanes bounded by measured cost -
rather than introducing implicit culling. Any density rule must be visible to
the user and documented; a feature must never disappear silently.

### S5 - Resize and wake-up path (H5)

Only after S0 attributes the native gap across event delivery, viewport
synchronization, repaint/wake-up scheduling and compositor work: coalesce
intermediate resize events so one relayout serves a burst, remove redundant
repaint wake-ups on the DNA-window path, and establish whether the snapshot
writer or the renderer explains the acceptance-harness timings. Keep CPU paint,
semantic-harness waits and release-binary timings distinct in every report.
Require the S0 macOS cross-check before a general scheduling change. If the
delay belongs to the Xvfb harness or snapshot writer, repair/label that path
rather than changing production scheduling; a Linux-specific product defect
can still warrant a narrowly evidenced fix. Worker migration, virtualization
or a single-root GUI workspace remain conditional on runtime profiling.

### S6 - Startup and first usable content (H6; may run first)

Codex owns the smallest product fix justified by S0's startup attribution;
Glen owns its repeatable acceptance. If this is the dominant confirmed delay,
prioritize S6 ahead of S1-S5 rather than leaving it until the end.

- Separate process launch to responsive root window, requested project load,
  DNA-window admission/hydration and first subject-correct feature content.
  Retain time spent in harness waits, external I/O, GPU/window initialization
  and project parsing separately; do not silently remove costs from totals.
- Compare empty-profile startup and explicit public-project opening, with
  cold/warm conditions documented. Do not change project restoration, skip
  validation or discard annotations merely to make startup look faster.
- Fix only the confirmed blocking phase. Where that phase is the existing
  constructor/hydration path, implement S3 once and report its contribution to
  both open-window and process-to-content totals. Preserve cancellation,
  stale-result checks and sequence identity for any background work.
- Agree numerical startup targets with the owner and Glen after baseline
  attribution, before optimization. Exit requires target-bound before/after
  native evidence on the selected projects, not a splash screen or a quicker
  harness checkpoint. If no product defect remains, retain that finding and
  test the actual startup against the agreed targets; do not infer acceptance.

## Correctness guardrails

- Extend the existing brute-force oracles in `src/render_dna_linear.rs` tests
  (feature interval index vs. exhaustive overlap, half-open edges, generation
  invalidation) to cover every new cache or reuse path introduced here.
- Caches are keyed by explicit generations and revisions, never by wall clock or
  frame counter; every new key gets an invalidation regression.
- Preserve subject bindings, cancellation semantics and scientific outputs. A
  performance change that alters an export, a report or a coordinate is a
  defect, not a trade-off.
- Keep visual-regression fixtures for the DNA map; compare before/after at
  identical viewport and window size.

## Acceptance

1. Implementation posts the S0 table plus the confirmed/rejected hypothesis
   ledger at one named SHA, with `cargo test -p gentle-benchmarks --bench
   gui_operations` green and the fixture regeneration command recorded.
2. The auditor repeats the headless ladder on a stable host with `bench-audit`
   and retains raw Criterion artifacts and environment metadata. Baseline
   repeats share a revision; before/after comparisons name both revisions and
   hold effective profile settings, host, toolchain and fixture identity
   constant, spanning no profile redefinition. Do not pool
   different profiles or treat a build-speed improvement as a runtime result.
3. Native acceptance uses the actual packaged GUI and covers open, pan, zoom,
   live resize, layer toggles and selection on PATZ1 and TP73 plus the exact agreed envelope boundary. Record
   startup/S6 targets, the macOS cross-check and circular regression checks;
   retain the largest ladder fixture as a separate stress outcome. Record the
   repeated same-SHA native profile comparison under Measurement Profiles;
   neither it nor the headless ladder replaces packaged-profile acceptance.
4. The performance verdict is the auditor's. An unresolved owner scope choice
   or unbudgeted startup cannot be a green `.12` result. Implementation-side
   numbers are evidence, not acceptance.

## Non-goals without evidence

Virtualization, worker migration, GPU renderer replacement, caching or culling
strategies, and single-root GUI workspace restructuring are all conditional on
a confirmed, attributed runtime hotspot. This plan explicitly forbids starting
with them. Build-boundary work follows B0 instead; that carve-out is neither a
runtime-performance claim nor approval for a broad root-crate rewrite.

## Open questions

- Owner: confirm the interactive envelope and circular scope under Audit
  Admission Decisions before the release-target audit.
- Owner/Glen: freeze startup and interaction budgets after baseline attribution,
  before optimization; proposals above are not accepted measurements.
- Is the native gap compositor-bound on all three platforms, or specific to the
  Xvfb/Openbox acceptance environment?
