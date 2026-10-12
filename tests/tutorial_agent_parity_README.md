# Tutorial agent-parity fixtures

Origin: entirely hand-authored synthetic JSON and Markdown created in temporary
directories by `tutorial_agent_parity_cli.rs` and
`src/workflow_examples/agent_parity/tests.rs`. No external data or model calls.

Recreate with:

```sh
cargo test --lib tutorial_agent -- --test-threads=1
cargo test --lib workflow_examples::agent_parity -- --test-threads=1
cargo test --test tutorial_agent_parity_cli
```

The CLI test uses the actual `gentle_examples_docs tutorial-check` executable,
an unrelated working directory and deliberately absent workflow/manifest/output
inputs. Invalid contracts must be reported before replay loading or execution.
Both LF and CRLF fixture text are tested. Unit cases exercise semantic JSON
binding, retained unapproved placeholders, path quoting, raw traversal/Windows
paths, and Unix symlink escapes. Unix symlink tests do not claim native Windows
symlink acceptance. All files are disposable; no project or biological result
is changed. These tests validate contracts, not inner-agent language quality.
