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

Status reviewed: 2026-09-21 at `70a3b038`. S0's developer tools landed in
`5893aa35` and are smoke-tested; its timed/native audit and S1-S5 remain pending.
This consolidated plan has one hypothesis ledger and one S0 section below.
The `.11` release gate is unchanged; counts of work are not timing evidence or
performance acceptance.

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

Steady painting of this small 20,802-bp locus with 75+ loaded features is not the
measured hotspot. The investigation therefore separates (a) one-shot
construction and hydration work on the UI thread, (b) repeated viewport/resize
work, and (c) the unattributed distance between an OS event and confirmed
content, which the headless benchmark cannot observe. These historical numbers
are not density-ladder acceptance; its timing still awaits the external auditor.

## Scope

In scope: linear and circular DNA-map feature presentation, the feature tree and
layer-visibility panel that mirror it, the per-window hydration path that first
makes features visible, and the repaint/resize scheduling around them.

Out of scope for `.12`: visual redesign, new feature classes, engine-side
annotation work, sequence-text panel work beyond what the feature path forces,
and any change to displayed biology. Retain
[DEC-026](decisions.md#dec-026-long-gui-jobs-use-optimistic-engine-snapshots)
and
[DEC-048](decisions.md#dec-048-command-submission-and-observation-never-wait-for-execution);
this plan adds no new execution authority.

## Proposed interaction budgets

These are proposals to be confirmed, tightened or replaced by the auditor, and
are only meaningful bound to one host, toolchain, profile, fixture hash and
GENtle revision.

| Interaction | Subject | Proposed budget |
| --- | --- | ---: |
| First confirmed feature content after open | <= 250 kbp, <= 5,000 features | <= 500 ms, no single frame > 100 ms |
| Pan / zoom step | same | p95 <= 16.7 ms CPU prepare+paint |
| Live resize | 1920x1080 | every intermediate frame <= 33 ms; confirmed content <= 150 ms |
| Feature layer toggle (CDS/repeat/array/TFBS) | same | <= 100 ms to confirmed content |
| Selection / hover | same | <= 16.7 ms, no model rebuild |
| Tab switch / multiwindow focus | same | <= 100 ms to confirmed content |
| Process start to first usable window | reference project | attributed, not yet budgeted |

## Hypotheses And Evidence Ledger

The evidence below separates verified work from unmeasured timing contributions.
Before optimization, record measurements across the fixture ladder for H1-H4
and native traces for H5. A rejected hypothesis stays documented with its
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
one-base pan invalidates both caches. Layer-count misses also call
`GcContents::new_from_sequence_with_bin_size` across the full sequence, so this
is O(all features + bases), not merely O(visible features). Tree/count
hits/builds, feature visits and `layer_gc_bases` record the work separately.
The deterministic pan regression observes tree builds 1 -> 2 and GC bases
1,200 -> 2,400 without a project mutation; dominance at 2 Mbp remains a timing
hypothesis. Static grouping and viewport counts must be separated carefully:
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
Native resize-to-confirmed-content substantially exceeded the CPU-side first
frame after resize. Attribution across event delivery, viewport
synchronization, repaint scheduling, GPU upload, compositor work and the
`gui-test-support` snapshot writer does not exist yet. Until it does, no
conclusion about the renderer is warranted from the native numbers. Direct
display setters in the headless harness measure presentation work, not X11
input, user navigation, window focus or compositor behavior. Scheduling changes
remain conditional on a native trace without the snapshot writer.

## Slices

### S0 - Measurement Foundation (Tools Implemented; Audit Pending)

**Implemented and smoke-tested; reuse these tools:**

- `dna_feature_latency` benchmarks the real `MainAreaDna` constructor, hydration,
  first frame with tree deferred/loaded, steady frame, one-base pan, zoom, mRNA
  layer toggle, feature selection, hover and three resize transitions. It emits
  117 cases and 81 counter observations across all nine combinations of
  20 kbp / 250 kbp / 2 Mbp and 100 / 1,000 / 10,000 features. Length and count
  vary independently; the largest case is a stress probe, not an interactive
  performance promise.
- Each fixture binds exact sequence and feature bytes. The generator uses
  structured, overlapping plus/minus transcripts, CDS, exons, regulatory and
  repeat features, with half clustered in the first 5 kbp. It uses no private
  annotation, random state, network, prepared genome or binary fixture blob.
- `scripts/dna_feature_latency.py prepare` builds offline once and binds the
  executable, source/diff, toolchain, profile and lockfile. `run` rechecks that
  binary and executes it directly with isolated profile/cache/temp directories.
  Build time is separate; failure, timeout, raw logs and work tables are retained.
  Timed audits reject dirty-source or development-profile receipts.
- `GENTLE_DNA_CACHE_DIAGNOSTICS=1` exposes an opt-in DNA-viewer diagnostics pane.
  It is observation-only, reads renderer counters without waiting on its lock,
  contains no sequence/feature names, and writes nothing to project state.
- Additional Puffin scopes distinguish feature painting, interval indexing,
  construction computations and hydration substeps. Existing tree build/render,
  layer-count and display-sync scopes are retained.

See [the benchmark runbook](../benches/README.md#dna-feature-density-latency)
for exact generation, build and prebuilt replay commands. Keep the existing
`gui_operations` TP73/PATZ1 workload: the synthetic subtree harness does not
replace engine-backed window/report hydration or real annotations.

**Pending audit and exit criteria (Glen):**

1. Freeze a clean SHA and run two prebuilt `bench-audit` repeats on a stable
   host. Retain the per-fixture, per-interaction tables, raw distributions and
   counter deltas; assess repeatability rather than treating two successful
   processes as a timing pass. Bind every result to source revision, profile,
   toolchain, fixture and binary hashes.
2. Capture native input-to-content traces at 820x520, 1200x800, 1600x1000 and
   1920x1080. Separate process startup, open, clone/hydration, first paint,
   event delivery, repaint and compositor delays. Use a release-like binary
   without the semantic snapshot writer; capture locally, not on per-frame
   network storage. Retain project/report hashes and tool/environment identity.
3. Record the confirmed/rejected/unresolved timing hypotheses above and choose
   an interactive density ceiling from evidence. Only then authorize the
   corresponding S1-S5 changes against their evidence requirements below.

The headless runner does not establish native acceptance, even when every smoke
case passes. S0's implementation is not to be repeated; its measurement gate
remains open.

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
Worker migration, virtualization or a single-root workspace remain out of scope
without profiling evidence.

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
2. The auditor repeats the ladder on a stable host with a release-like binary,
   retains raw Criterion artifacts and environment metadata, and compares only
   matching profile/host/toolchain/fixture/revision baselines.
3. Native acceptance covers open, pan, zoom, live resize, layer toggles and
   selection on PATZ1 and TP73, plus the largest ladder fixture, with
   before/after bound to two SHAs.
4. The performance verdict is the auditor's. Implementation-side numbers are
   evidence, not acceptance.

## Non-goals without evidence

Virtualization, worker migration, GPU renderer replacement, caching or culling
strategies, and workspace restructuring are all conditional on a confirmed,
attributed hotspot. This plan explicitly forbids starting with them.

## Open questions

- Which feature-count ceiling does `.12` commit to supporting interactively?
- Does the circular map need the same treatment in `.12`, or only after the
  linear path is attributed?
- Is the native gap compositor-bound on all three platforms, or specific to the
  Xvfb/Openbox acceptance environment?
