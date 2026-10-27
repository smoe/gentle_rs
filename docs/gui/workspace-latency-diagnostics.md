### Workspace Latency Diagnostics

<!-- manual-counterparts:start -->
Related CLI reference: [Latency Diagnostics](../cli/diagnostics.md).
These topic links do not imply identical feature coverage across interfaces.
<!-- manual-counterparts:end -->

- For debugging Splicing Expert launch stalls, launching GENtle with
  `GENTLE_TRACE_SPLICING_EXPERT=1` emits coarse open/render stage tracing to
  the terminal and mirrors it into the sequence-window status line.
- For debugging RNA-read Mapping launch stalls, launching GENtle with
  `GENTLE_TRACE_RNA_READ_MAPPING=1` emits coarse open/render stage tracing to
  the terminal and mirrors it into the sequence-window status line.
- For debugging Promoter design launch stalls, launching GENtle with
  `GENTLE_TRACE_PROMOTER_DESIGN=1` emits coarse open/render stage tracing to
  the terminal and mirrors it into the sequence-window status line.
- For GUI latency investigations, build with `--features gui-profiler` and
  launch with `GENTLE_GUI_PROFILE=1 cargo run --features gui-profiler --bin
  gentle`. GENtle starts a localhost Puffin stream on `127.0.0.1:8585` by
  default; override it with `GENTLE_GUI_PROFILE_ADDR=127.0.0.1:PORT` and inspect
  the frame spans with `puffin_viewer`. The profiling spans are coarse by
  design and focus on app-frame, DNA-window shell, DNA-map, sequence-text, and
  feature-tree costs. Gene-locus investigations additionally expose separate
  scopes for GUI composition/preview rendering and the engine's occupancy,
  motif, probe-effect, and complete report-building phases.
- DNA latency scopes additionally separate feature drawing/indexing, deferred
  record lock/clone, viewport reconciliation, map/overlay refresh, and the
  restriction, ORF, methylation and GC constructor computations. Launch with
  `GENTLE_DNA_CACHE_DIAGNOSTICS=1` to expose an optional **DNA cache diagnostics**
  pane without a profiler build. Its cumulative counters show tree/layer model
  builds, feature visits, GC bases traversed, layout/index builds, and existing
  restriction/overlay/GC/display-sync cache activity. A busy or circular renderer
  reports no linear counters, not zero work. Reading counters does not wait on
  the renderer, rebuild caches or alter project state; counters contain no
  sequence or feature names. They count work, not milliseconds. See the
  [latency plan](../dna_feature_rendering_latency_plan.md) and
  [prebuilt audit runbook](../../benches/README.md#dna-feature-density-latency).
  `tree_viewport_updates` counts in-place visibility recounts separately from
  `tree_builds`; these reuse hits are also included in `tree_hits`. Older logs
  without this counter do not establish a zero recount cost.
  GC layer counts use bin geometry without rescanning sequence bases on pan;
  `layer_gc_bases` therefore stays zero on that path. This does not suppress
  the separate GC-value computation needed for the actual GC track.
  For startup attribution, set `GENTLE_GUI_STARTUP_TRACE` to a new JSON output
  path before launching the desktop binary. The bounded, name-free CPU trace is
  written only on exit, not on the paint path; no screenshots, database queries
  or project writes are introduced. Failed loads and loading placeholders do
  not become content checkpoints. See the [startup runbook](../../benches/README.md#startup-phase-checkpoints)
  for missing-event handling and the distinction from native presentation and
  performance acceptance. The offline `scripts/gui_startup_summary.py` reader
  produces per-subject CPU tables and retains the exact trace without relaunching
  GENtle; missing durations stay unavailable, not zero. Help preparation has
  separate manual, shell-reference, tutorial-discovery and selected-load spans,
  plus first-use help/tutorial/menu-discovery spans and path-free, session-wide
  SVG/cache/rasterization totals. These are attribution only: they do not change
  help loading or prove that rasterization causes the observed startup delay.
  This diagnostic is off by default.

<!-- manual-page-nav:start -->

---

[GUI Manual](../gui.md) | [Previous: Configuration Window](configuration-window.md) | [Next: GENtle Shell (GUI)](gentle-shell-gui.md)
<!-- manual-page-nav:end -->
