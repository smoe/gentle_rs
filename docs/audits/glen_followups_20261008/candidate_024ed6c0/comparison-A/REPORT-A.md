# Help/Startup Comparison A — Lazy Discovery

## Identity and method

- Before: `e723f46eb3fe58b15ba4bad5ff975024cd761437`, binary SHA-256
  `ec9e288692afb8dc34e2a34fb02bba3107811b5a029068fa6ba1924d32bb74a3`.
- After: `64baba4f1cd21c9df731b232d9f82ce55c9d6495`, binary SHA-256
  `346e7b034dfb214d67a806d75895b6d9b4a0fbc4c764e9605b3349611e021f1d`.
- Both: Rust/Cargo 1.96.1, default features, ordinary `dev`, one build job,
  Linux x86_64, Xvfb 1600x1000x24 plus Openbox, 266 installed font rows.
- Runtime resources and working directory were fixed to the After checkout.
  `Cargo.toml`, `Cargo.lock`, `docs/tutorial`, and `assets` were checked for
  the intended recipe/resource relationship before the builds.
- Public input only:
  `docs/tutorial/inputs/gibson-project.gentle.json`, SHA-256
  `9ae0a7f82c8cb6298693c6736896df2526815b8af12e78d7499522db05a4c766`.
- For each revision: three new-TMPDIR Cold runs and three Warm runs reusing
  one cache populated by an excluded fill. HOME and all XDG state were fresh
  per run. Every included run ended normally with 52 retained events, no
  dropped events, no failed spans, and no dropped image observations.
- Fixed interaction: launch the public project; wait for an externally painted
  usable main-window control; open the subject-correct
  `XXU13852 [n-1] (gibson_destination_pgex)` DNA window; open GUI Help; open
  Tutorial; switch to `01.02 Contribute to GENtle development`; open the
  File > Open Tutorial Project submenu; close normally.
- Values below are median `[min, max]` in milliseconds. Inclusive trace spans
  are reported independently and are never summed.

## Results

| Condition | Revision | usable main | subject window | help preparation | tutorial discovery |
| --- | --- | ---: | ---: | ---: | ---: |
| Cold | Before | 34457.376 [34305.855, 34486.570] | 35535.908 [35377.002, 35635.603] | 33714.684 [33549.240, 33730.462] | 33666.608 [33502.320, 33682.097] |
| Cold | After | 2671.788 [2602.472, 2674.802] | 3732.066 [3707.692, 3777.489] | 50.631 [48.758, 51.297] | 2.457 [2.436, 2.474] |
| Warm | Before | 2650.298 [2648.927, 2656.588] | 3762.287 [3696.341, 3763.885] | 82.548 [79.322, 82.651] | 34.000 [32.977, 34.102] |
| Warm | After | 2644.200 [2599.443, 2680.515] | 3751.480 [3621.532, 3772.256] | 49.192 [49.118, 51.018] | 2.422 [2.386, 2.517] |

Cold usable-main median improved by 31,785.588 ms and the subject-correct
window by 31,803.842 ms. Warm medians are effectively unchanged.

First-use work did not inherit the removed startup delay:

- first GUI Help window: Before Cold 181.172 ms; After Cold 121.782 ms;
- tutorial window: Before Cold 4887.558 ms; After Cold 4770.819 ms;
- Topic switch: Before Cold 257.310 ms; After Cold 257.870 ms;
- tutorial-project discovery span: Before Cold 3.049 ms; After Cold 1.567 ms.

Before Cold recorded 29 SVG references, 29 completed rasterizations, zero
failures, and median 33,629,831 microseconds of rasterization work. Before Warm
recorded the same 29 references as cache hits. After Cold and Warm recorded no
SVG references or rasterization work in this fixed sequence. The producer did
not emit a `help_tutorial_open` span for this UI route; that metric is
**unavailable**, not zero. The externally observed tutorial-window timing is
retained instead.

## Verdict

**Pass for the bounded lazy-discovery comparison.** The Cold improvement is
large, repeatable, and localized to tutorial discovery; it does not transfer
the removed delay to first Help, Tutorial, Topic, menu use, or the first public
DNA subject. Warm behavior remains unchanged.

The measurements do not justify additional shared-font loading, failed-SVG
caching, or menu caching: After performs no image work in the measured path,
there are no SVG failures, and menu discovery is already about 1.6 ms. This
Linux virtual-X11 result is not physical-platform or release acceptance.

Raw per-run values and spreads are in `analysis-A/runs.csv`, `runs.json`, and
`statistics.json`. Harness-development runs that lacked a requested boundary
remain explicitly retained under `rejected-*` and are excluded from statistics.
