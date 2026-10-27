# Protocol Examples and Tutorials

<!-- manual-counterparts:start -->
Related GUI reference: [Detailed Reference](../gui/detailed-reference.md).
These topic links do not imply identical feature coverage across interfaces.
<!-- manual-counterparts:end -->

Structured workflow examples:

- canonical source files: `docs/examples/workflows/*.json`
- schema: `gentle.workflow_example.v1`
- includes GUI-first parity skeleton:
  - `docs/examples/workflows/vkorc1_rs9923231_promoter_luciferase_assay_planning.json`
  - merged canonical tutorial:
    `docs/tutorial/08-04_vkorc1_warfarin_promoter_luciferase_gui.md`
- on-demand snippet generation (CLI/shell/JS/Lua):
  - `cargo run --bin gentle_examples_docs -- generate`
- validation only:
  - `cargo run --bin gentle_examples_docs -- --check`
- tutorial generation:
  - `cargo run --bin gentle_examples_docs -- tutorial-generate`
- tutorial drift check:
  - `cargo run --bin gentle_examples_docs -- tutorial-check`
- test gating:
  - `always` examples run in default tests
  - `online` examples run only with `GENTLE_TEST_ONLINE=1`
  - set `GENTLE_SKIP_REMOTE_TESTS=1` to force-skip all remote-resource tests
    (takes precedence over `GENTLE_TEST_ONLINE`)
- `skip` examples are syntax-checked only


## Protocol-first example source

Canonical workflow examples are adapter-neutral JSON files:

- source: `docs/examples/workflows/*.json`
- generated adapter snippets: `docs/examples/generated/*.md`
- tutorial runtime manifest: `docs/tutorial/manifest.json`
  (generated from `docs/tutorial/sources/`)
- committed generated tutorial output: `docs/tutorial/generated/`

Regenerate snippets on demand:

```bash
cargo run --bin gentle_examples_docs -- generate
```

Validate example files and schema without writing output:

```bash
cargo run --bin gentle_examples_docs -- --check
```

Generate tutorial markdown + retained artifacts:

```bash
cargo run --bin gentle_examples_docs -- tutorial-generate
```

Validate committed tutorial-generated output:

```bash
cargo run --bin gentle_examples_docs -- tutorial-check
```

Before replay, this command validates every source unit declaring
`agent_parity`: the documented command must occur in its guide, parse through
the shared shell, and respect the existing automatic-execution mutation and
undo/redo guards. A guide may call an inspect/export command non-mutating even
when the parser conservatively requires review; that is allowed with `ask`.

Learner-created `@file` requests can bind a `parser_payload` JSON template;
see [the source contract](../testing.md#tutorial-agent-parity-contracts). Guide,
template and `@file` paths resolve under `--repo-root`, not the process working
directory. They must be portable repository-relative paths; absolute paths,
raw parent traversal and symlink escapes are rejected. Approval placeholders
are parsed unchanged, never replaced with invented approval values.

Parity failures exit 1 on stderr with the prefix
`tutorial agent parity failed (N findings):`, followed by one JSON finding per
line (`code`, `tutorial_id`, `case_id`, `detail`). Success adds an
`agent_parity` summary with `tutorials`, `cases`, `declared_mutating` and
`parser_state_mutating` counts to the existing JSON summary. The parity step
never executes commands; **the rest of `tutorial-check` still replays workflows
in temporary directories**. This is not a model, GUI or scientific acceptance
verdict. Custom legacy manifests without source units keep their existing
behavior and report zero declared parity contracts.

<!-- manual-page-nav:start -->

---

[CLI Manual](../cli.md) | [Previous: Introspection, Facts and Readiness](discovery.md) | [Next: Independent Command Examples](command-examples.md)
<!-- manual-page-nav:end -->
