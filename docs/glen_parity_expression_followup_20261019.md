# Glen: Parity And Expression Inspector Follow-Up

Updated after Glen's `33e67234` return: the original parity/expression campaign
is finished. This file now requests only the two bounded repairs below after
the owner's next ordinary push, not another full campaign or source freeze.

After Steffen's push, fetch `smoe/gentle_rs main` once and record its full SHA
and tree. Use a clean detached checkout and an explicitly set, absolute,
run-owned `CARGO_TARGET_DIR`; preserve primary checkouts, shared caches and
historical evidence. No interim release freeze or recurring bundle is needed.

## Only The Pending Checks

The [verified return and Windows failure](glen_audit_followups_20261008.md#parity-and-expression-verification-at-33e67234-2026-10-19)
close parity freshness, catalog/manifest, UI intents, MCP report/refusal,
deterministic agent-command/control-card, locator and Python groups at
`33e67234`. Do not repeat those groups or the completed `394dd5d` audits.
Normal CI may naturally exercise them at newer source; preserve their original
source-bound passes without upgrading the whole failed Windows job to green.

Recheck the corrected hosted-layer test (`7bc3158f`) and the separate expression
MCP fact-annotation repair found in full Windows CI:

On Linux, run from the checkout root, with no regeneration or source edits.
Set the fresh run-owned target before any Cargo invocation:

```sh
: "${CARGO_TARGET_DIR:?Set an absolute run-owned target directory}"
case "$CARGO_TARGET_DIR" in /*) ;; *) printf '%s\n' 'CARGO_TARGET_DIR must be absolute' >&2; exit 1 ;; esac
cargo test --locked --lib -j1 app::protein_expression_ui::tests::protein_expression_inspector_small_hosted_window_paints_without_mutation -- --exact --test-threads=1
cargo test --locked --lib -j1 protein_expression_inspector_ -- --test-threads=1
cargo test --locked --lib -j1 execute_introspect_protein_expression_mcp_contract_ -- --test-threads=1
cargo test --locked --lib -j1 execute_introspect_capabilities_projects_full_registry_with_fact_annotations -- --test-threads=1
cargo test --locked --lib -j1 execute_introspect_all_capability_section_has_full_fact_annotations -- --test-threads=1
```

Require actual counts 1/7/1/1/1, with no zero-test or skipped pass. The window
must be foreground on opening, then middle after consuming the focus request,
with the report and engine state preserved in both frames. The direct and
aggregate introspection routes must retain the shared read-only/no-project
contract, seven optional MCP fields and no state effects; inspection readiness
is not product or expression readiness. No business-logic or focus-policy
change is intended. Retain failures and environment limitations honestly.

Four
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
