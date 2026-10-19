# Glen: Parity And Expression Inspector Follow-Up

After Steffen's push, fetch `smoe/gentle_rs main` once and record its full SHA
and tree. Use a clean detached checkout and an explicitly set, absolute,
run-owned `CARGO_TARGET_DIR`; preserve primary checkouts, shared caches and
historical evidence. No interim release freeze or recurring bundle is needed.

## Only The Pending Checks

Verify the `8f37453b` tube parity repair and the newer expression inspector
(`5893ebb0`) plus MCP/agent/locator follow-up (`b3b3e8a9`) are included.
Do not repeat the completed `394dd5d` README, sequence-design/actual-CLI,
expression/storage starter, workflow, rack/labware or nine-file carrier audits.
Their [recorded passes](glen_audit_followups_20261008.md#focused-verification-at-394dd5d-2026-10-18)
remain source-bound; normal CI may naturally exercise them at newer source.

On Linux, run from the checkout root, with no regeneration or source edits.
Set the fresh run-owned target before any Cargo invocation:

```sh
: "${CARGO_TARGET_DIR:?Set an absolute run-owned target directory}"
case "$CARGO_TARGET_DIR" in /*) ;; *) printf '%s\n' 'CARGO_TARGET_DIR must be absolute' >&2; exit 1 ;; esac
cargo build --locked --no-default-features -j1 --bin gentle_examples_docs
DOCS="$CARGO_TARGET_DIR/debug/gentle_examples_docs"
"$DOCS" parity-matrix-check
"$DOCS" tutorial-catalog-check
"$DOCS" tutorial-manifest-check
cargo test --locked -j1 --test parity_matrix_freshness
cargo test --locked --lib -j1 protein_expression_inspector_ -- --test-threads=1
cargo test --locked -p gentle-shell -j1 ui_intent
cargo test --locked --lib -j1 mcp_server::protein_expression::tests:: -- --test-threads=1
cargo test --locked --lib -j1 protein_expression_inner_agent_suggestion_ -- --test-threads=1
cargo test --locked --lib -j1 protein_expression_agent_card_ -- --test-threads=1
cargo test --locked --lib --features gui-test-support -j1 protein_expression_capture_targets_ -- --test-threads=1
python3 -m unittest scripts.test_tutorial_checkouts scripts.test_specialist_tutorials.SpecialistSourceTests -v
```

Require seven inspector tests, three named MCP tests, one inner-agent suggestion
test, one control-card test and one feature-gated locator test; the shared
UI-intent filter must also select tests. A zero count/skip is not acceptance.
Retain the actual commands, counts and any environment/resource limitations.

Check both tube rows/counts survive the new MCP registry projection. Compare
inspector/shared-shell/MCP/agent-command reports on the same saved synthetic
state and inputs, excluding only `generated_at_unix_ms`. Missing inputs must
stay unanswered; explicit provider projection/withholding, unknown/path/null
refusals, copy-only actions, stale/closed-worker rejection and unchanged
project/history/files must hold. MCP cannot see unsaved GUI changes.

The locator test is an in-memory semantic check, not a desktop capture. Four
[retained macOS teaching images](screenshots/tutorial_gui_acceptance/protein_expression_handoff_cli/README.md)
remain `captured_not_accepted` at their recorded source, without retained typed
GUI reports/snapshots. Do not relabel them as this revision's runtime proof.

## Return Through Steffen

Send one concise result for forwarding to this main Codex chat: per-check
pass/fail/not-run/blocked, exact SHA/tree, OS/toolchain/features/profile,
lockfile/binary hashes, commands/counts/skips and the normal same-source CI
URL/outcome. Distinguish primary failures from cascade skips. Supply a small
manifested folder with raw streams, command/exit receipts and any concrete
defect; retain bulky binaries and captures on Glen's host. A later main tip
does not turn this scoped result into a new synchronization task.

No push, dispatch/rerun, provider/live-model call, new desktop capture,
package/container build, tag, publication or scientific approval is requested.
Normal CI can be observed without restarting it. The separate six-gate/native,
scroll, package and scientific boundaries stay open; new gel acceptance is
not an additional manual task in this request. Codex ran no local Rust.
