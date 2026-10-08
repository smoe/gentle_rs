# Follow-Up Evidence, Not Historical Replacement

`scroll_cpu_smoke_d1bf3fd3.json` is copied byte-for-byte from the
`gui-scroll-cpu-smoke-d1bf3fd321c651483125ee6f588f8129490e0cd5` artifact of
[Actions run 37766941695](https://github.com/smoe/gentle_rs/actions/runs/37766941695).
Receipt SHA-256:
`60add8b059510d660c36a25c1b7115138b3fbc8693dcecf2b42fd2316ddceab5`.

Reproduction: on the exact source revision, build the locked `gui_operations`
target in the development profile, prepare PATZ1 with
`scripts/prepare_real_patz1_tutorial.py`, and run the retained executable with
`--test` and both `GENTLE_GUI_BENCH_PATZ1_STATE` and
`GENTLE_GUI_BENCH_PATZ1_REPORT` set to that prepared project's paths. The
workflow definition at the tested SHA records the exact build commands.

The Ubuntu 24.04/x86-64 run used Rust 1.99.0. All three fixtures, four sizes and
two directions passed (24 observations). The downloaded receipt, executable,
prepared inputs and all 25 listed artifact hashes were independently checked.
Executable SHA-256:
`d9a10afc01266334935ae83286acef417dc7d60af0be4dee19a062ebb811d6d8`.

This accepts deterministic wheel-event and next-frame **CPU benchmark cases**,
not performance limits, native input-to-paint latency, package-opt1 artifacts,
the full CI run or human scientific approval. The run's Linux/Windows test
compiles failed separately on a sibling-method visibility error, subsequently
repaired in `1aee220c`. Glen's original audit and evidence are not overwritten.
