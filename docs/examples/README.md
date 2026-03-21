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

## Draft design-resource examples

- Draft, non-executable planning/design artifacts may live outside
  `docs/examples/workflows/`.
- Current draft resource example:
  - `docs/examples/plans/gibson_destination_first_single_insert.json`
  - schema: `gentle.gibson_assembly_plan.v1`
- These files are documentation resources for protocol/data-model design and
  are not yet executed by workflow-example tests.

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

## Generate tutorial output

Generated tutorial pages and retained runtime artifacts are committed under:

- `docs/tutorial/generated`

Source manifest:

- `docs/tutorial/manifest.json` (schema `gentle.tutorial_manifest.v1`)

Generate tutorial output:

```bash
cargo run --bin gentle_examples_docs -- tutorial-generate
```

Validate committed tutorial output against fresh generation:

```bash
cargo run --bin gentle_examples_docs -- tutorial-check
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
