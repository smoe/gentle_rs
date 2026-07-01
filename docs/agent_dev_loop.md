# Agent Development Loop

Pattern: an agent drives the working-tree GENtle binary directly, keeps state in
an explicit session file, and validates each code edit before continuing.

## Short Wrapper

Use the in-tree wrapper instead of spelling out `cargo run` each time:

```sh
scripts/dev-gentle-cli capabilities
```

The wrapper runs `cargo run -q --bin gentle_cli -- ...`, passes arguments to
`gentle_cli` unchanged, and preserves stdout/stderr/exit status. If
`CARGO_TARGET_DIR` is already set, it respects it. Otherwise it uses a shared
workspace target at `../../.gentle_target_shared` so parallel worktrees do not
all rebuild from cold. Set `GENTLE_SHARED_TARGET_DIR` to choose a different
fallback without overriding Cargo for the whole shell.

For benchmark-style calls, put `--release` before the GENtle command:

```sh
scripts/dev-gentle-cli --release capabilities
```

Every invocation tees stdout to
`$GENTLE_DEV_LOG_DIR/gentle-cli-*.json`, defaulting to `/tmp/gentle-dev`, while
still printing the same stdout to the terminal:

```sh
rg '"schema"|"ok"' /tmp/gentle-dev
```

## Pre-flight

Before a longer agent session, ask the CLI to check the local loop:

```sh
scripts/dev-gentle-cli --state /tmp/gentle-agent-session.json doctor --agent
```

The JSON report covers cargo availability, whether the current binary is newer
than source mtimes, whether the `--state` path is writable, and whether
`gentle_examples_docs --check` passes. The command exits `0` only when every
reported check is healthy.

## Discover Operations

Use capabilities as the first discovery step:

```sh
scripts/dev-gentle-cli capabilities
```

The payload includes `capability_registry`, projected from
`gentle_protocol`, so CLI/MCP/JavaScript/Lua discovery share the same operation
names, descriptions, schemas, mutability flags, and adapter exposure metadata.

For MCP parity without installing anything, drive the in-tree binary over
stdio:

```sh
printf 'Content-Length: 58\r\n\r\n{"jsonrpc":"2.0","id":1,"method":"tools/list","params":{}}' \
  | cargo run -q --bin gentle_mcp --
```

## Thread State

Use one explicit state file for a session and pass it to every mutating call:

```sh
STATE=/tmp/gentle-agent-session.json
scripts/dev-gentle-cli --state "$STATE" state-summary
scripts/dev-gentle-cli --state "$STATE" op @operation.json
scripts/dev-gentle-cli --state "$STATE" workflow @workflow.json
scripts/dev-gentle-cli --state "$STATE" export-state /tmp/session-copy.gentle.json
```

This keeps the agent loop stateless from the shell's point of view while still
letting GENtle persist project context when the biology or workflow requires it.

## Inspect Artifacts

Most render and export routes write explicit files. Keep paths under `/tmp` or
the current branch's artifact area, then inspect by file type:

```sh
scripts/dev-gentle-cli --state "$STATE" panels render-isoform-svg SEQ PANEL /tmp/panel.svg
scripts/dev-gentle-cli --state "$STATE" protocol-cartoon render-svg PROTOCOL_ID /tmp/cartoon.svg
scripts/dev-gentle-cli --state "$STATE" state-summary > /tmp/state-summary.json
```

Use `rg`/`jq` for JSON, a browser or image viewer for SVG/PNG, and the logged
stdout files under `/tmp/gentle-dev` when comparing prior calls.

## Extend From Inside The Loop

After each source edit, run:

```sh
cargo check -q
```

For the current minimal validation slice from `CONTRIBUTING.md`, run:

```sh
cargo check -q
cargo test -q workflow_examples -- --test-threads=1
cargo run --bin gentle_examples_docs -- --check
cargo run --bin gentle_examples_docs -- tutorial-check
```

Use narrower targeted tests first when working on one adapter or operation, then
run this slice before handoff.
