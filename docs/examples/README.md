# GENtle Protocol-First Examples

This folder defines canonical workflow examples independent of any one adapter
syntax.

## Source of truth

- Canonical examples live in `docs/examples/workflows/*.json`.
- Schema: `gentle.workflow_example.v1`.
- Each file carries:
  - metadata (`id`, `title`, `summary`, `test_mode`, `required_files`)
  - canonical `workflow` JSON payload

`test_mode` values:

- `always`: validated and executed by default test runs
- `online`: executed only when `GENTLE_TEST_ONLINE=1`
- `skip`: parsed and documented, but not executed in automated tests

## Generate adapter snippets

Generate Markdown snippets for CLI/shell/JS/Lua:

```bash
cargo run --bin gentle_examples_docs -- generate
```

Generated output is written to `docs/examples/generated`.

Validation only:

```bash
cargo run --bin gentle_examples_docs -- --check
```

## Test examples

Run default (offline-safe) example tests:

```bash
cargo test workflow_examples -- --test-threads=1
```

Run online examples (explicit opt-in):

```bash
GENTLE_TEST_ONLINE=1 cargo test workflow_examples -- --test-threads=1
```
