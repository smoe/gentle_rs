# GENtle Performance Audit Targets

Developers verify that benchmark fixtures and scientific assertions still run
without collecting a statistical baseline:

```bash
cargo test --bench gui_operations --features benchmark-support
cargo test --bench specificity_finalization --features benchmark-support
```

Timed runs and release-facing comparisons belong to the external auditor.

## GUI-critical operations

The primary Criterion target exercises the real DNA-window constructor and
embedded egui render path with a feature-free 120 kbp control and the public
TP73 locus fixture:

```bash
cargo bench --bench gui_operations --features benchmark-support
```

It reports eager constructor, deferred UI-thread hydration, first-frame, and
steady-frame CPU costs. The hydration benchmark directly exercises the
`replace_loaded_sequence` boundary used after lazy background loading, while
keeping the background record clone outside the timed loop. Benchmark IDs
include sequence length, feature count, and the first 12 characters of the
input SHA-256. This headless target cannot measure native window creation, GPU
submission, compositor latency, input delivery, or whether interaction feels
responsive. Those remain part of the external GUI/Puffin acceptance described
in `docs/testing.md`.

## Specificity finalization

Run the deterministic primer-specificity benchmark with:

```bash
cargo bench --bench specificity_finalization --features benchmark-support
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

Compare performance only on the same runner, toolchain, build profile, and
fixture-hash ID. Shared CI runners are suitable for compilation or smoke runs,
not strict regression thresholds.

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
cargo bench --bench gui_operations --features benchmark-support -- --save-baseline AUDIT_LABEL
```

Keep the console log and the corresponding `target/criterion/` subtree with
the audit record. Compare only against a baseline produced by the same host and
toolchain. The auditor decides whether a difference is meaningful after also
performing the real GUI interaction check; GENtle does not turn a noisy shared
runner result into a release verdict.
