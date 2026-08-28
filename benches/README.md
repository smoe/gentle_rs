# GENtle Performance Audit Targets

Developers verify that benchmark fixtures and scientific assertions still run
without collecting a statistical baseline:

```bash
cargo test --bench gui_operations --features benchmark-support
cargo test --bench specificity_finalization --features benchmark-support
```

Timed runs and release-facing comparisons belong to the external auditor.

## Audit profiles

GENtle provides two deliberately non-comparable optimized modes:

- Routine audit uses `--profile bench-audit`. It keeps release `opt-level=3`
  and stripping, but uses thin LTO, 16 codegen units, and `panic=unwind` to
  reduce cold-build pressure and make cached audits practical. Cargo already
  forces unwind for benchmark targets; matching that setting for the binary
  targets built alongside an integration benchmark avoids compiling the root
  library once for each panic strategy.
- Exact release-like audit uses plain `cargo bench` and therefore preserves the
  existing fat-LTO, one-codegen-unit benchmark build. Use it when compilation
  cost and exact release-like code generation are part of the question.

Cargo keeps compiled artifacts under `target/bench-audit/` and
`target/release/`, respectively. Criterion 0.8.2 does not infer the Cargo
profile and otherwise stores both modes under the shared `target/criterion/`
tree. Set `CRITERION_HOME` as shown below so measurements remain under the
matching profile directory. Never compare results across these modes: the
routine mode deliberately differs in LTO, codegen units, and panic strategy.
The first routine build can still be substantial. Cargo automatically prepares
enabled binary targets for an integration benchmark so it can expose
`CARGO_BIN_EXE_*`; those targets therefore retain their own Thin-LTO link cost.
Subsequent runs reuse the profile cache. GENtle's build fingerprint follows the
current loose branch ref and consults repository-wide `packed-refs` only when
that loose ref is absent, so Git maintenance in another worktree does not
invalidate an otherwise unchanged audit build.

## GUI-critical operations

The primary Criterion target exercises the real DNA-window constructor and
embedded egui render path with a feature-free 120 kbp control and the public
TP73 locus fixture:

```bash
CRITERION_HOME=target/bench-audit/criterion \
  cargo bench --profile bench-audit --bench gui_operations \
  --features benchmark-support -- --quick --noplot
```

For the exact release-like mode, omit `--profile` and change only the retained
result directory:

```bash
CRITERION_HOME=target/release/criterion \
  cargo bench --bench gui_operations --features benchmark-support \
  -- --quick --noplot
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
CRITERION_HOME=target/bench-audit/criterion \
  cargo bench --profile bench-audit --bench specificity_finalization \
  --features benchmark-support -- --quick --noplot
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
CRITERION_HOME=target/bench-audit/criterion \
  /usr/bin/time -v cargo bench --profile bench-audit \
  --bench gui_operations --features benchmark-support -- \
  --quick --noplot --save-baseline bench-audit_HOST_TOOLCHAIN_REVISION
```

Use `/usr/bin/time -l` instead of `-v` on macOS. Keep the console log and the
corresponding `target/bench-audit/criterion/` subtree with the audit record.
For an exact release-like baseline, use plain `cargo bench`, set
`CRITERION_HOME=target/release/criterion`, and prefix the baseline label with
`release-like`. Compare only against a baseline produced by the same profile,
host, toolchain, fixture hash, and GENtle revision. The auditor decides whether
a difference is meaningful after also performing the real GUI interaction
check; GENtle does not turn a noisy shared runner result into a release verdict.
