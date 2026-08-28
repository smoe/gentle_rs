# GENtle Testing Strategy (Draft)

This document proposes automated testing for both operation correctness and
rendering correctness.

## 1. Test pyramid

1. Engine unit tests (fast, deterministic)
2. Protocol/CLI integration tests (JSON in/out, state persistence)
3. Rendering snapshot tests (linear + circular graphics)
4. Optional GUI smoke tests (window opens, toggles wired)

## 2. Engine tests (required)

Scope:

- operation semantics (`LoadFile`, `Digest`, `Ligation`, `Pcr`, `ExtractRegion`,
  `FilterByDesignConstraints`)
- structured errors/warnings
- deterministic operation log shape
- display toggle operations (`SetDisplayVisibility`)

Execution:

```bash
cargo test engine::tests::
```

## 3. CLI/protocol tests (required)

Scope:

- JSON contract compatibility
- `@file.json` input handling
- state persistence (`--state`)
- command parity (`op`, `workflow`, `state-summary`, import/export)

Recommended method:

- add integration tests that run `gentle_cli` as a subprocess
- compare normalized JSON output against golden files

### 3.1 Canonical workflow example tests

Canonical protocol examples in `docs/examples/workflows/*.json` are now part of
the test surface.

Execution:

```bash
cargo test workflow_examples -- --test-threads=1
```

Online-only examples are opt-in:

```bash
GENTLE_TEST_ONLINE=1 cargo test workflow_examples -- --test-threads=1
```

`test_mode` policy from each example file:

- `always`: parsed, validated, and executed in default runs
- `online`: executed only with `GENTLE_TEST_ONLINE=1`
- `skip`: parsed/validated only

Primer-design CI policy:

- Always-run workflow/tutorial examples that exercise primer design should pin
  `primer_design_backend=internal` explicitly so routine CI stays fast and does
  not depend on a locally installed `primer3_core`.
- Real external-binary primer-design checks are opt-in via:

```bash
GENTLE_TEST_EXTERNAL_BINARIES=1 cargo test real_primer3 -- --test-threads=1
```

- Those opt-in tests are intended to validate the actual `primer3_core`
  integration path separately from the default deterministic internal-backend
  CI surface.

### 3.2 Tutorial drift and runtime checks

Tutorial source + generated outputs are part of the test surface:

- source units: `docs/tutorial/sources/`
- generated runtime manifest: `docs/tutorial/manifest.json`
- generated output: `docs/tutorial/generated/`

Validation commands:

```bash
cargo run --bin gentle_examples_docs -- tutorial-check
cargo run --bin gentle_examples_docs -- tutorial-manifest-check
cargo run --bin gentle_examples_docs -- tutorial-catalog-check
cargo test workflow_examples -- --test-threads=1
```

`tutorial-check` regenerates tutorial outputs in a temp directory and compares
byte-for-byte with committed files under `docs/tutorial/generated`. On failure,
it prints a paste-ready tutorial feedback context with the chapter id, source
JSON, workflow JSON, generated chapter path, artifact directory, failing check,
and suggested GitHub issue-template category.

CI additionally runs a CLI smoke path for core tutorial chapters via:

- `cargo run --bin gentle_cli -- workflow @docs/examples/workflows/<core>.json`

### 3.3 Planned opt-in Mistral inner-agent conformance test

Add one live-provider integration test for the native Mistral Agent Assistant
transport. This is a provider conformance test, not a deterministic-output
snapshot and not part of default offline CI.

Activation and isolation:

- If `MISTRAL_API_KEY` is missing or blank, print one explicit skip line and
  return before model discovery or generation. Do not use `#[ignore]`, because
  that would also suppress the test when the key is deliberately supplied.
- When the key is present, use the existing `mistral_large_native` transport
  with `--model mistral-small-latest` by default. Allow
  `GENTLE_TEST_MISTRAL_MODEL` to override the model without changing the
  checked-in agent catalog.
- Use a temporary project state and report directory. Never print the API key,
  persist it in project state, or include it in failure snapshots.
- Disable retries and use a bounded request timeout so a provider outage cannot
  stall the suite indefinitely.

Test routine:

1. Run native-Mistral preflight and require successful authentication, model
   discovery, and selected-model resolution.
2. Submit the fixed, offline restriction-map prompt from the Agent Interfaces
   tutorial. Request one `execution=ask` suggestion and no network-backed
   biological operation.
3. Validate `gentle.agent_response.v1`, require exactly one non-recursive
   command, and parse it through `parse_shell_line`.
4. Resolve the suggested command's canonical capability descriptor and require
   `annotation_status=fact_annotated`; reject a command not represented by
   introspection.
5. Bind its arguments and require `introspect readiness` to report `ready`
   against the temporary state before execution.
6. Execute through the shared shell path with auto-execution disabled, then
   verify declared hard effects through `introspect verify-effects`.
7. For the restriction scan, compare normalized scientific rows with a direct
   shared-shell run: enzyme, recognition interval, recognition sequence, cut
   geometry, and matched-site count. Ignore prose, provider token usage,
   elapsed time, generated timestamps, and option ordering.
8. Run a second state-aware prompt that refers to the now-loaded `agent_toy`
   sequence. Require the suggestion to use that real sequence id rather than
   inventing another project object, then repeat parse/readiness checks.

The pass condition is deterministic containment: a variable model proposal is
accepted only when it lands inside GENtle's fixed schema, parser, capability,
readiness, execution, and effect-verification contracts. Repeating the provider
call and demanding byte-identical prose or JSON would test provider sampling,
not GENtle determinism.

## 4. Rendering tests (required)

Graphics are part of functionality. Visibility toggles and biology overlays must
be tested at the rendered output level.

### 4.1 Proposed renderer test contract

Introduce deterministic rendering export functions from a shared view model:

- `render_linear_svg(view_model) -> String`
- `render_circular_svg(view_model) -> String`

Use fixed viewport size, fonts, and deterministic layout order.

### 4.2 Snapshot workflow

For each test case:

1. Build state through engine operations
2. Build view model with selected display settings
3. Export linear and circular SVG
4. Compare against approved snapshot files

You can export from persisted engine state with:

```bash
cargo run --bin gentle_cli -- --state state.json render-svg SEQ_ID linear out.linear.svg
cargo run --bin gentle_cli -- --state state.json render-svg SEQ_ID circular out.circular.svg
```

Suggested snapshot folders:

- `tests/snapshots/linear/*.svg`
- `tests/snapshots/circular/*.svg`

### 4.3 Cases to include first

- baseline map with all overlays on
- each visibility toggle off/on individually
- dense feature labels
- linear strand direction correctness
- restriction site label crowding
- circular zero-point crossing features
- design-constraint filter pass/fail cases (GC bounds, homopolymer cap, U6
  `TTTT` rejection, forbidden motifs)

## 5. Regression gates

Recommended CI gate (minimum):

1. `cargo check -q`
2. engine unit tests
3. CLI integration tests
4. rendering snapshot tests

If snapshot output changes, require reviewer approval and explicit snapshot update.

### 5.1 External performance audit

Criterion benchmarks characterize CPU-bound, in-process GENtle algorithms;
they do not replace functional assertions or end-to-end acceptance. The
primary target exercises the actual eager DNA-window constructor, deferred
UI-thread hydration boundary, and first and steady embedded egui frames over a
feature-free control and the public TP73 locus:

```bash
cargo bench --bench gui_operations --features benchmark-support
```

Implementation-side verification stops at the non-statistical smoke path:

```bash
cargo test --bench gui_operations --features benchmark-support
cargo test --bench specificity_finalization --features benchmark-support
```

This is a CPU-side proxy for GUI-critical preparation and painting. It cannot
measure native viewport creation, GPU/compositor behavior, event delivery, or
perceived interaction stalls. Glen is the named external auditor for the
current release assessment and therefore combines those results with an
exact-revision GUI run using the existing `gui-profiler` feature/Puffin scopes
and the DNA-viewer resize/repaint harness. The auditor, not the implementation
author, retains baselines and assigns the performance verdict.

The secondary target covers complete-query BLAST HSP parsing and
primer-specificity interpretation at 100, 1,000, and 6,600 total rows:

```bash
cargo bench --bench specificity_finalization --features benchmark-support
```

The benchmark constructs its synthetic input outside timed loops, verifies the
expected hit, pairing, intended-product, and off-target counts once, reports
throughput in HSP rows, and includes a fixture SHA-256 token in every benchmark
ID. It launches no BLAST process and reads no real reference dataset. Criterion
results should only be compared for the same runner, toolchain, fixture hash,
and optimized profile. Pull-request CI may compile or smoke the targets, but
strict timing thresholds belong to repeated external-auditor runs on a stable
host after several baselines exist. The exact commands and retained evidence
are documented in [`../benches/README.md`](../benches/README.md).

Full Ensembl FASTA retrieval, peak RSS, real BLAST/Primer3 execution, GUI frame
behavior, and Xvfb/macOS acceptance remain system-level measurements rather
than unit-test claims.

## 6. Semantic GUI acceptance identifiers

GENtle exposes a bounded, read-only semantic widget snapshot only when built
with `--features gui-test-support`. The feature is disabled by default. Set
`GENTLE_GUI_TEST_SNAPSHOT` to an explicit output path when launching `gentle`;
the app atomically replaces that JSON file after each root frame. There is no
listener, input-injection API, screenshot permission bypass, or command route.

Semantic IDs are lowercase dotted purpose names such as
`dna.splitter.info_width` and do not derive from translated labels or screen
position. Repeated domain rows share one semantic ID and carry an opaque
`subject_scope` hash derived from stable domain fields rather than list index.
Snapshots never include visible labels, sequence strings, paths, credentials,
or arbitrary user text. Rectangles are screen-relative egui logical points,
including the native viewport origin; multiply them by the item's
`pixels_per_point` for physical pixels.

Example launch:

```bash
cargo build --locked --features gui-test-support --bin gentle
GENTLE_GUI_TEST_SNAPSHOT=/tmp/gentle-gui.json \
  target/debug/gentle --project /path/to/tp73-proof.project.gentle.json
jq '.generation, .settled, .items[] | select(.semantic_id == "dna.splitter.info_width")' \
  /tmp/gentle-gui.json
```

The snapshot schema is `gentle.gui_semantic_snapshot.v1`. `generation`
advances monotonically and `settled` is false while known deferred GUI work,
including initial feature-tree hydration, is pending. Harnesses should wait on
those fields rather than a fixed sleep. Semantic IDs locate GUI affordances;
scientific correctness still comes from the typed project state and reports.
The Linux example [`../scripts/gui_semantic_xvfb_smoke.sh`](../scripts/gui_semantic_xvfb_smoke.sh)
self-skips when X11 tools or a display are unavailable and accepts explicit
opaque repeat/array row scopes plus a typed-state verification command.
than Criterion microbenchmarks.

## 6. Practical implementation order

1. Keep extending engine tests alongside new operations
2. Add CLI integration tests for current protocol
3. Implement deterministic SVG exporter from view model
4. Add snapshot suite and CI wiring

## 7. Notes

- PNG snapshot testing is possible but less stable than SVG/text snapshots.
- Prefer view-model-to-SVG tests over full GUI pixel tests for determinism.
- Keep one small GUI smoke test suite for wiring confidence only.
