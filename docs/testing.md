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

## 6. Practical implementation order

1. Keep extending engine tests alongside new operations
2. Add CLI integration tests for current protocol
3. Implement deterministic SVG exporter from view model
4. Add snapshot suite and CI wiring

## 7. Notes

- PNG snapshot testing is possible but less stable than SVG/text snapshots.
- Prefer view-model-to-SVG tests over full GUI pixel tests for determinism.
- Keep one small GUI smoke test suite for wiring confidence only.
