# Glen Audit: Current DNA-Feature Boundary Baseline

Prepared: 2026-09-27. This is an external-auditor CPU baseline under
[DEC-039](decisions.md#dec-039-external-auditor-owns-the-performance-verdict),
not native-GUI acceptance or a `.12` release verdict.

## Identity And Scope

- Revision: `0467867cce153ab4c89c4093c32c003b9f70e513`, the exact fetched
  `origin/main` at preparation and replay time.
- Source state: clean detached worktree; no untracked build inputs.
- Workload: `density_boundary_v1`, ten fixtures, 130 Criterion cases and 90
  untimed work-counter observations.
- Profile: `bench-audit`; `rustc 1.95.0 (59807616e 2026-04-14)`; locked,
  offline preparation followed by direct binary replays.
- Benchmark binary SHA-256:
  `e25fc12f60b9badde60aa53124e4bcd32569aa8468aa337ca0c50ea730a1f2e2`.
- Retained raw evidence on Glen's host:
  `/home/clawbio/.openclaw/workspace/artifacts/gentle-gui-latency-audit-20260927`.
  The directory contains build logs/receipt and both complete Criterion trees,
  `run.json`, `run.log` and `work.tsv` files. The disposable Cargo target was
  removed after retention to recover disk space.

The two replay receipts have SHA-256
`74171d7017b4f1cff3b281bc2cbf66eb5bd5bb651b3acf2a124cb00c08b72f4c`
and `c6a31e7102819da850a75947fc0fe2869c079d35bb0bbcf6c61b859a048ef2b2`.
Both report `failure: null`, 130 statistics, 90 observations and
`interactive_boundary_exercised: true`.

## Build Observation

The cold audit preparation completed in 28:54.65 with exit status 0, no swap,
and 21,023,304 KiB maximum resident set size for the timed process tree. This
is the thin-LTO auditor recipe, not the internal `dev` package recipe. It
therefore does not contradict the lower package-build RSS recorded for the B0
extraction, and it remains unsuitable as evidence that a small hosted runner
can build the audit binary.

## Exact 250 kbp / 5,000-Feature Boundary

The two direct replays completed in 46.686 and 51.893 seconds. Times below are
Criterion mean point estimates in milliseconds; the final column is their
two-repeat median. Quick-mode repeats characterize the baseline but do not set
an acceptance threshold.

| Case | Replay 1 | Replay 2 | Median |
| --- | ---: | ---: | ---: |
| Constructor | 4.510 | 4.548 | 4.529 |
| Hydration | 0.002 | 0.002 | 0.002 |
| First frame, tree deferred | 31.338 | 32.283 | 31.810 |
| First frame, tree loaded | 69.615 | 71.806 | 70.711 |
| Steady frame | 2.468 | 2.357 | 2.413 |
| Pan by one base | 61.491 | 63.901 | 62.696 |
| Zoom | 51.700 | 52.956 | 52.328 |
| Toggle mRNA | 50.768 | 51.649 | 51.209 |
| Select | 3.595 | 3.516 | 3.555 |
| Hover | 2.318 | 2.363 | 2.340 |
| Resize 820x520 | 16.315 | 17.343 | 16.829 |
| Resize 1600x1000 | 17.046 | 17.558 | 17.302 |
| Resize 1920x1080 | 17.083 | 17.987 | 17.535 |

The boundary fixture hash is
`e61ecd4986ed34c3d072c76f0fa502c7a2825a9175a00601f983b19357e9d67f`.
Its enabled derived-layer inventory is nonempty and stable: 26 restriction
groups (1 viewport eligible), 2,500 GC bins (50), 30 ORFs (1), and 52
methylation sites (2), over viewport `[0, 5000]` with 100-bp GC bins.

Pan, zoom and mRNA toggle each rebuild the feature tree, layer counts and one
linear layout once. They scan zero GC bases. Steady, select and hover rebuild
none of those units. Each resize performs one layout and no tree/layer rebuild.
Both replays have the same counter and inventory results.

## Continuity With The Retained GC Baseline

The prior post-GC baseline at `ee9eb483` used the same host, toolchain and
`bench-audit` profile, but the legacy nine-fixture `density_ladder_v1`
workload. Only its 117 content-identical case IDs are compared here; the new
boundary is not interpolated from them.

Using the median of two quick repeats at each revision, the common cases show a
geometric aggregate shift of +1.72% and a median shift of +1.43%. Of 117 cases,
104 are within +/-5%, one is at least 5% faster, and twelve are at least 5%
slower. Most larger percentage changes are approximately two-microsecond
hydration cases. The pan/zoom/mRNA-toggle cases at 250 kbp and 2 Mbp across
100, 1,000 and 10,000 features range from +0.16% to +4.78%. These two quick
repeats do not establish a runtime regression; they support continuity of the
post-GC baseline while retaining the observed +6.53% 20-kbp/10,000-feature pan
sample as an inconclusive outlier.

## Verdict And Next Evidence

The current headless boundary workload is complete and internally consistent.
It confirms the zero-GC-rescan work reduction at the agreed 250-kbp/5,000-
feature boundary. It does not measure native event delivery, compositor delay,
window creation, project loading, circular-map behavior or scientific/visual
correctness.

The follow-up [native Linux audit](glen_native_gui_audit_20260927.md) exercises
a local exact-identity internal-`dev` rebuild. PATZ1 and the synthetic TSS path
are visually correct in that run, while fresh startup is dominated by a
repeatable 37.5--37.8 s help-preparation span. Exact downloaded-package,
nonempty DuckDB-hit and macOS evidence remain open. Circular maps keep
correctness coverage without a new `.12` timing promise.
