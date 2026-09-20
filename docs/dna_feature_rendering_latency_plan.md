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

Steady painting of 20,802 bp with 75+ loaded features is therefore not the
problem. The remaining latency sits in (a) one-shot construction and hydration
work on the UI thread, (b) work that is redone whenever the viewport or the
window rectangle changes, and (c) the unattributed distance between an OS event
and confirmed content, which the headless benchmark cannot observe. PATZ1 is
also a *small* locus; feature-count scaling has never been measured.

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
Native resize-to-confirmed-content exceeded the CPU-side first frame after
resize by roughly 20x-100x. Attribution across event delivery, viewport
synchronization, repaint scheduling, GPU upload, compositor work and the
`gui-test-support` snapshot writer does not exist yet. Until it does, no
conclusion about the renderer is warranted from the native numbers.

## Slices

### S0 - Measurement harness (prerequisite, blocks S1-S5)

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

Exit criteria: one command produces a per-fixture, per-interaction table; two
repeats on a stable host agree; every number carries source revision, profile,
toolchain and fixture hash.

### S1 - Stop relayouting what did not change (H1)

Key `layout_features` results by an explicit signature (viewport, drawing
rectangle, vertical offset, display revision, feature generation) and reuse the
previous layout when the signature is unchanged; keep selection and hover out of
the signature, since they only affect drawing. Reuse the seed and lane buffers
across relayouts instead of reallocating, and avoid per-feature `String`
allocation for labels and kinds that only feed comparisons.

Evidence required first: H1 confirmed on the ladder, with the relayout share of
a resize frame reported.

### S2 - Decouple tree and layer counts from the viewport (H2)

Split both models into a viewport-independent part (grouping, filtering,
ordering, labels - keyed by feature generation and display revision) and a cheap
viewport pass that only marks visibility and recounts, reusing the existing
feature interval index rather than scanning all features. Do not rebuild either
model while the panel is collapsed or scrolled out of view.

Evidence required first: H2 confirmed, including the rebuild cost at 10^3 and
10^4 features during a continuous pan.

### S3 - Take one-shot work off the first frames (H3/H4)

- Move derivable presentation work out of the hydration frame: precompute what
  the background load can already produce, and time-slice the remainder across
  frames with a visible, bounded progress state rather than one long frame.
- Do not recompute restriction sites, ORFs, methylation and GC content on the UI
  thread at window construction when the engine record already carries them;
  where a recompute is genuinely needed, make it explicit, cancellable and off
  the first paint.
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
