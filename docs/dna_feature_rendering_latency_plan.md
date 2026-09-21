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

Implementation update: the existing S0 density tools are implemented and
smoke-tested. The added startup phase markers, exact owner-selected boundary
cases and macOS cross-check below remain pending, as do Glen's timed/native
audit and runtime slices S1-S6. B0 separates build feedback from runtime
performance. This plan does not change the `.11` candidate; counts of work are
not timing evidence or performance acceptance.

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
| Native evidence is Linux/Xvfb only | Cross-check on macOS before S5 | Require a local release-like macOS trace and retain platform-specific conclusions |
| Three overlapping priority lists | Keep one implementation sequence | Roadmap owns release scope, this plan owns steps, and the historical acceptance report links here |

## Audit Admission Decisions

These two owner decisions were requested on 2026-09-21 and remain pending.
Recommendations are not accepted performance promises:

- **Interactive envelope:** recommend linear loci up to 250 kbp and 5,000
  loaded features; retain 2 Mbp and 10,000 features as stress cases. If selected,
  add at least the exact 250 kbp / 5,000-feature boundary with the same clustered
  annotation and deferred/loaded-tree modes. The current 100/1,000/10,000 ladder
  cannot certify that boundary by interpolation. Include representative visible
  restriction/derived layers and record their counts, not only input features.
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

## Improvements To The Original Plan

1. **An explicit layout signature is not yet a fix.** `RenderDnaLinear::render`
   already checks the rectangle and layout dirty flag. A signature containing
   the same rectangle still misses on every resize. First separate actual
   layout builds, interval-index builds and paint cost; then remove a measured
   redundant dependency or allocation, with invalidation tests.
2. **Vary length and count independently.** Use all nine combinations of
   20 kbp / 250 kbp / 2 Mbp and 100 / 1,000 / 10,000 features. Otherwise
   whole-sequence recomputation can masquerade as feature-count scaling.
   Measure both deferred and explicitly loaded feature trees. The 10,000-feature
   case is a stress workload, not an interactive performance promise.
3. **Keep visibility honest.** Static grouping and viewport visibility are
   different products. A bounding-span interval index is not an exon-piece
   overlap oracle: an intron-only viewport must not make a transcript exon
   visible. Collapsed or offscreen groups still need correct visible/total
   counts when queried. Do not simply remove the viewport from existing keys.
4. **Measure another H2 dependency.** Layer-count cache misses currently call
   `GcContents::new_from_sequence_with_bin_size` across the full sequence.
   Record bases traversed separately from feature visits; the layer-count cost
   is not exclusively a tree problem.
5. **Existing computed features are not proof of freshness.** Before skipping
   construction-time recomputation, establish sequence, enzyme-catalog and
   parameter identity. Worker migration must retain DEC-026 snapshot checks,
   cancellation, exact results and a responsive observer, as required by DEC-048.
6. **CPU is not native latency.** Direct display setters in a headless egui
   benchmark measure presentation work, not X11 events, user navigation, GPU
   upload, wake-up scheduling, window focus or compositor behavior.

## Hypotheses And Evidence Ledger

| Hypothesis | Instrumented boundary | Current interpretation |
| --- | --- | --- |
| H1: rectangle changes trigger costly layout | `RenderDnaLinear::layout_features`, `draw_features`, layout/index counters | Existing dirty check verified; density/timing contribution awaits audit |
| H2: pan rebuilds whole-locus tree/count models | `FeatureTree::build_model`, tree/count hits and builds, feature visits and GC bases | Viewport-keyed rebuilds verified in source and a deterministic pan regression; runtime share awaits audit |
| H3: hydration performs presentation work in one UI frame | `WindowDna::poll_deferred_load.hydrate`, replacement, viewport reconciliation, map update, overlay refresh | Background lock/clone and foreground hydration remain separate; native cost awaits audit |
| H4: constructor recomputes whole-sequence derived features | Separate restriction, ORF, methylation and GC scopes in `DNAsequence::update_computed_features` | Called by the constructor; length-scaling contribution awaits audit |
| H5: native input-to-content gap | Existing public native acceptance plus a new auditor trace without snapshot writer | Not explained by the CPU harness; scheduling changes remain blocked |
| H6: startup/first-locus dominates perceived delay | Proposed process entry, app construction, project load, first usable root window, open request and first subject-correct paint markers | About 32 s observed only in the debug acceptance harness; product time versus fixed waits/I/O/snapshot overhead is unresolved |

Counters measure work, not time. A counter hit does not establish that the
cache is fast, and a miss does not establish that it dominates a frame.

## S0: Measurement Foundation

Implemented developer tools:

- `dna_feature_latency` benchmarks the real `MainAreaDna` constructor, hydration,
  first frame with tree deferred/loaded, steady frame, one-base pan, zoom, mRNA
  layer toggle, feature selection, hover and three resize transitions. It emits
  117 cases and 81 counter observations across nine synthetic fixtures.
- Each fixture binds exact sequence and feature bytes. The generator uses
  structured, overlapping plus/minus transcripts, CDS, exons, regulatory and
  repeat features, with half clustered in the first 5 kbp. It uses no private
  annotation, random state, network, prepared genome or binary fixture blob.
- `scripts/dna_feature_latency.py prepare` builds offline once and binds the
  executable, source/diff, toolchain, profile and lockfile. `run` rechecks that
  binary and executes it directly with isolated profile/cache/temp directories.
  Build time is separate; failure, timeout, raw logs and work tables are retained.
- `GENTLE_DNA_CACHE_DIAGNOSTICS=1` exposes an opt-in DNA-viewer diagnostics pane.
  It is observation-only, reads renderer counters without waiting on its lock,
  contains no sequence/feature names, and writes nothing to project state.
- Additional Puffin scopes distinguish feature painting, interval indexing,
  construction computations and hydration substeps. Existing tree build/render,
  layer-count and display-sync scopes are retained.

See [the benchmark runbook](../benches/README.md#dna-feature-density-latency)
for exact generation, build and prebuilt replay commands. Keep the existing
`gui_operations` TP73/PATZ1 workload: the new synthetic subtree harness does not
replace engine-backed window/report hydration or real annotations.

**S0 exit still requires Glen:** freeze a clean SHA, run two `bench-audit`
repeats on a stable host, compare distributions and counter deltas, and capture
native input-to-content traces at 820x520, 1200x800, 1600x1000 and 1920x1080.
Separate process startup, open, clone/hydration, first paint, event delivery,
repaint and compositor delays. Use a release-like binary without the semantic
snapshot writer; capture locally, not on per-frame network storage. Retain
project/report hashes and tool/environment identity. The headless runner does
not claim this acceptance, even when every smoke case passes.

Before a general S5 scheduling change, also obtain a short local macOS native
resize/open trace on the same source and public fixture, with its own binary,
profile, window sizes and environment recorded and no semantic snapshot writer.
One trace is a diagnostic cross-check, not a macOS performance verdict. Do not
pool Linux and macOS timings or attribute differences to an OS from unlike
profiles. Windows native acceptance remains separate when making Windows claims.

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
  project were 27.852 ms eager DNA-window construction, 32.742 ms deferred
  UI-thread hydration, 28.1-29.5 ms first embedded frame at all four viewport
  sizes, 0.786-0.793 ms steady embedded frame, and 3.83-3.94 ms first frame
  after a resize.

Optimized steady painting of this small PATZ1 locus is inexpensive in the
embedded benchmark. It neither explains the debug harness's startup observation
nor proves all painting fast. Construction/hydration, viewport-dependent work,
native event-to-content delay and process/project startup remain distinct
hypotheses. The roughly 32 seconds includes fixed waits and capture overhead;
do not call it 32 seconds of application startup or subtract the optimized
Criterion timings from it. Feature-count scaling still awaits the auditor.

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

## Hypotheses, to be confirmed or rejected before any optimization

Each hypothesis gets an instrumented scope or counter, a measurement across the
fixture ladder, and a recorded result. A rejected hypothesis is documented with
its evidence rather than silently dropped.

**H1 - Full relayout on every changed rectangle.**
In `src/render_dna_linear.rs`, `RenderDnaLinear::render` recomputes
`layout_features` whenever the drawing rectangle changed or the layout dirty
flag is set. During a live drag-resize that is one full relayout per frame. The
relayout allocates a `Seed` per visible feature, including two `String`s (label,
upper-cased kind), then sorts and lane-allocates them. Expected to scale with
visible feature count, which PATZ1 barely exercises.

**H2 - Viewport-keyed model rebuilds.**
`FeatureTreeCacheKey` (`src/main_area_dna/feature_tree_ui.rs`) and
`LayerVisibilityCacheKey` (`src/main_area_dna.rs`) both contain the viewport,
and `active_linear_viewport_range` returns `Some` for every linear sequence. A
one-base pan therefore invalidates both caches and rebuilds models that are
O(all features), not O(visible features), even while the tree is collapsed or
scrolled out of sight.

**H3 - One-shot hydration inside a single UI frame.**
`WindowDna::poll_deferred_load` (`src/window_dna.rs`) calls
`replace_loaded_sequence` -> `replace_active_dna` (`src/main_area_dna.rs`),
which invalidates every derived cache, reconciles the viewport, updates the map
renderer and refreshes the construct-reasoning overlay in one frame (32.742 ms
on PATZ1). The background thread only clones the record; all presentation work
lands on the UI thread.

**H4 - Whole-sequence recompute at window construction.**
`MainAreaDna::new` calls `ensure_local_restriction_site_catalog_current`, which
runs `DNAsequence::update_computed_features` (`src/dna_sequence.rs`): restriction
sites, restriction groups, ORFs, methylation sites and GC content over the whole
sequence. That is linear in sequence length and unrelated to what the first
frame shows. Its cost is part of the 27.852 ms eager construction on 20 kbp and
must be measured at 250 kbp and 2 Mbp.

**H5 - The native gap.**
Debug native resize and optimized embedded-frame timings differ markedly, but
are not comparable measurements. Attribution across event delivery, viewport
synchronization, repaint scheduling, GPU upload, compositor work and the
`gui-test-support` snapshot writer does not exist yet. Cross-check on native
macOS before generalizing Xvfb/Openbox observations into a product fix.

**H6 - Startup and first usable content.**
The harness combines process startup, project/window opening, fixed waits and
capture work. Add phase markers around `src/bin/gentle.rs`,
`GENtleApp::new_with_project`, project load and the first usable root/window
frames; connect the open request to existing hydration/paint scopes. Use an
empty clean profile and the same public PATZ1/TP73 projects. A small constructor
benchmark cannot exonerate or explain process startup.

## Slices

### S0 - Measurement harness (prerequisite for runtime changes S1-S6)

- Deterministic, offline, hash-bound **feature-density fixture ladder**:
  synthetic annotated sequences at roughly 10^2 / 10^3 / 10^4 features over
  20 kbp / 250 kbp / 2 Mbp, with exon-structured transcripts, regulatory tracks
  and overlapping lanes, plus the existing public PATZ1 and TP73 fixtures.
  Regeneration must be a documented command, not a stored binary blob.
- Extend `gui-profiler` Puffin scopes to the currently unscoped feature path:
  `draw_features`, feature-tree model build vs. render, layer-visibility count,
  hydration sub-steps and the construct-reasoning overlay refresh.
- Surface the existing hit/miss counters (restriction-site presentation, overlay
  presentation, GC, engine display sync, feature tree, layer visibility) through
  one diagnostic line or debug pane so cache thrash is observable without a
  profiler build.
- Prebuilt, quick-running benchmark/acceptance runner with isolated profile
  baselines; measure build and link cost separately so it is never reported as
  runtime latency.
- Native input-to-content timing on a release-like binary **without** the
  semantic snapshot writer, at the four PATZ1 viewport sizes (`820x520`,
  `1200x800`, `1600x1000`, `1920x1080`), capturing locally - never to per-frame
  network storage, since the discarded CIFS harness blocked in kernel I/O.
- Separate process startup, first DNA-window construction, deferred hydration
  and first paint, so the acceptance harness's ~32 s startup-plus-first-locus
  observation is attributed rather than restated. Retain exact source,
  toolchain, profile, project and report hashes plus native traces with every
  recorded measurement.
- Add the owner-selected envelope boundary where the current ladder lacks it,
  and retain a local macOS native cross-check before general S5 changes. Startup
  instrumentation and baseline collection belong here; S6 owns the resulting
  product fixes, not another open-ended profiling task.

Exit criteria: one command produces a per-fixture, per-interaction table; two
repeats on a stable host agree; every number carries source revision, profile,
toolchain and fixture hash. The workload, circular scope and budgets are
recorded; stress outcomes and release-target outcomes remain separate.

### B0 - Build feedback (independent of runtime hypotheses)

The dedicated benchmark crate still depends on the root library with
`desktop-gui`; a cold or changed-source build can be expensive. However,
`scripts/dna_feature_latency.py prepare` already builds once and `run` replays
the hash-bound binary without Cargo. Use that path first. It solves repeat-run
rebuilds, not initial compilation or iteration after a source edit.

If build cost still blocks work, retain cold and incremental Cargo timings,
compile versus LTO/link phases, wall time and peak RSS. Evaluate the smallest
evidence-backed build/dependency-boundary change; a narrowly reviewed crate
extraction is eligible without proving a runtime hotspot. Follow DEC-006/007,
preserve production profiles and benchmark semantics, and remeasure build cost
at named revisions. Root source size alone does not prove which extraction
helps, and no broad engine/GUI split is authorized by this plan.

This is distinct from a single-root **GUI workspace** (replacing native child
windows with one application workspace), which remains a conditional UX/runtime
change. Neither faster linking nor a different audit profile establishes a
runtime speedup; build comparisons have their own acceptance evidence.

### S1 - Stop relayouting what did not change (H1)

`RenderDnaLinear::render` already checks the rectangle and a layout dirty flag.
A new signature with the same rectangle cannot avoid work during real resizing.
First remove a measured redundant dependency or allocation, preserving dirty
invalidation and hit areas. Keep selection and hover from unnecessarily dirtying
layout. Reuse seed/lane buffers only if measured allocation cost warrants it;
avoid per-feature label/kind allocation where comparisons can borrow data.

Evidence required first: H1 confirmed on the ladder, with the relayout share of
a resize frame reported.

### S2 - Decouple tree and layer counts from the viewport (H2)

Split both models into a viewport-independent part (grouping, filtering,
ordering, labels - keyed by feature generation and display revision) and a cheap
viewport pass that only marks visibility and recounts. Reuse an index only with
exact exon-piece and half-open overlap semantics, not transcript bounding spans
that include introns. Collapsed/offscreen models may defer work but must expose
correct counts when queried. Investigate arithmetic GC-bin counts separately
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
2. The auditor repeats the ladder on a stable host with a release-like binary
   and retains raw Criterion artifacts and environment metadata. Baseline
   repeats share a revision; before/after comparisons name both revisions and
   hold profile, host, toolchain and fixture identity constant. Do not pool
   different profiles or treat a build-speed improvement as a runtime result.
3. Native acceptance covers open, pan, zoom, live resize, layer toggles and
   selection on PATZ1 and TP73 plus the exact agreed envelope boundary. Record
   startup/S6 targets, the macOS cross-check and circular regression checks;
   retain the largest ladder fixture as a separate stress outcome.
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
