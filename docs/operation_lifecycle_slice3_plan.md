# Lifecycle Slice 3: Construction, Workflows And Implicit Persistence

Status, 2026-10-27: owner-authorized revised implementation authored; Rust,
actual-process, scripting-runtime and native verification remain with Glen/CI.
This replaces the original proposal, not the historical executed evidence for
slices 1--2. The owner supplied Claude's review and explicitly excluded the
new payload-resolution layer. No further Claude call was made. See the
[historical review prompt](operation_lifecycle_slice3_claude_review_prompt.md),
[parent proposal](operation_lifecycle_policy_proposal.md#revised-slice-order)
and [verification handoff](glen_next_verification_handoff.md#lifecycle-slice-3-2026-10-27).
No release gate, S2 authorization or scientific acceptance follows from this slice.

## Source And Scope

Implementation began on `gentle_rs_3` in
`/Users/u005069/.codex/worktrees/faae/gentle_rs` at `ce4d33b5` (then also local
main). The original inspection was at `320d5803`, later `7d7519ad`, in the
47dd review worktree. These are provenance, not an instruction to freeze main.
Unrelated `outputs/` and the previously paused goal stay untouched; no push.

Only `PlanDnaSequenceDesign` (with or without `path`) and
`EvaluateGenomeEditIntent` (S1) are migrated. `ApplyDnaSequenceDesign`, explicit
`SaveFile`, ordinary loading and unknown operations remain compatibility controls.
No biological, schema, parser, agent-approval, scheduler or writer change.

## One Engine-Owned Fact

`engine::LifecycleRequest::{Operation,Workflow}::preserves_project_content()` is
a pure selector over already decoded requests. It performs no file reads,
parsing, model calls or template expansion. False means the consumer retains its
existing fallback, not a universal instruction to save. In particular MCP's
existing RNA file-only no-save list remains intact.

| Request | Construction/apply repair | Checkpoints | Implicit project save |
| --- | --- | --- | --- |
| Plan, absent or present path; S1 | Skip persistent reconciliation | Preserve undo/redo | Suppress on migrated raw routes |
| Nonempty workflow of only those operations | Same | Existing per-step semantics | Suppress |
| Empty workflow, mixed or unmigrated request | Legacy construction; covered children still skip their own repair | Legacy per-step semantics | Existing host fallback |

`GentleEngine::from_state_for` retains fresh instance identity, history defaults
and operation-counter reseeding without repairing supplied state. Ordinary
`from_state` is unchanged. This is not a sticky engine mode: a subsequent
ordinary mutation still reconciles normally. `apply_internal` consults the same
fact before even a refused inspection can repair state.

Successful inspections still journal and advance execution bookkeeping. They do
not acquire checkpoints or advance mutation/structural revisions. Whole
`history_summary` equality is therefore the wrong oracle; compare checkpoint
contents and mutation/structural revisions separately from operation-log count.

## Adapter Boundaries

- Raw CLI op/workflow and raw MCP op/workflow use the decoded request for
  construction and implicit-save policy. MCP `confirm=true` still precedes
  payload parsing and project loading, including no-path previews and S1.
- Generic JavaScript/Lua operation/workflow bindings use the same factory.
  Compare decoded state-in/state-out semantically, not serialized byte order.
  Python remains a thin CLI delegate; it gains no separate policy.
- Shared-shell/named/forwarded CLI and MCP-shell helpers are **not refactored**.
  Their existing change flags already suppress inspection saves; regressions
  cover their file-backed spellings. They may construct an ordinarily repaired
  engine before resolving a payload; no general state-out promise is added.
- No resolve-once envelope, public payload helper, second grammar, early
  `@file` read or new captured-payload identity is introduced. CLI and shared
  shell retain their existing different file-intent/error rules. Conservative
  wrapper/candidate confirmation remains additive.
- Raw MCP's returned `state_summary` now describes the preserved post-loader
  state, not constructor-repaired containers/arrangements. Counts can differ
  from older versions without sequence/report differences. The standalone
  `state_summary` tool remains an ordinary construction route.

Missing project files remain missing. Empty, noncanonical, LF/CRLF project files
and candidate sidecar bytes remain untouched when no explicit writer targets
them. In-memory preservation starts **after** the existing loader, including
sidecar hydration and warning metadata. Permissive loading can still remove a
broken candidate reference in memory; unmigrated autosaving routes can persist
that loss. This residual is not repaired globally here. Strict/permissive load
behavior and explicit saves are unchanged.

## Workflows, Managed Jobs And Failures

Both workflow executors retain order, `run_id`, child journals/checkpoints and
cooperative cancellation. No nested `apply`, rollback or error-path save is
added. In-process failures retain an executed prefix; raw CLI/MCP save only on
overall success; scripting may raise without state-out; managed failures discard
unpublished deltas. Already written reports survive those failures.

Managed generic Op/Workflow remains a **journal-publishing project worker**,
even when all operations are covered inspections. Keep the loaded-target
snapshot, owner identity, structural/history stale checks and guarded commit.
Do not move it to RNA's independent observation lane or discard its journal.
Persisted content/checkpoints/dirty revisions are unchanged on inspection commit.
[DEC-026](decisions.md#dec-026-long-gui-jobs-use-optimistic-engine-snapshots)
distinguishes this from observation-only workers that never commit.

A concurrent inspection journal append can still make another managed job stale.
The test pins that boundary rather than weakening it. Cancellation before commit
can discard a computed journal while leaving a requested report on disk; report
publication is not proof of job commit or scientific completion.

Opaque agent-plan and macro/template wrappers that construct before discovering
the selected request remain legacy. JS/Lua regressions explicitly pin that
repair rather than claiming universal wrapper parity. Nested-agent guards,
approval requirements and managed admission remain unchanged.

## Authored Regressions And External Checks

New synthetic fixtures live inline in `src/engine/lifecycle_slice3_tests.rs` and
the Python companion, reusing tutorial 06.07's hand-crafted literal CDS and a
six-base invented S1 request. No private input or new byte-bound fixture.

- Engine: genuinely unreconciled construction, unique identity/counter reseeding,
  later mutation, state/revision/checkpoint-pointer and redo preservation on
  success/refusal/bad output/cancellation/infeasible/exhausted outcomes; both
  workflow executors, order/run IDs, failure/cancel prefixes and legacy controls.
- Service: loaded-target op/workflow, guarded journal publication with no dirty
  state/history loss; concurrent-journal stale rejection and cancellation after
  explicit report writing. Existing owner-replacement/stale tests remain controls.
- MCP: preview must leave state absent, denied Apply leaves it absent, approved
  Apply still creates exact DNA; shared-shell `@file` helper bytes stay unchanged.
- Real JS/Lua VM tests: generic operation/workflow state-out, inline/loaded
  targets, refusal/mutation/empty controls, plus opaque-wrapper legacy controls.
- Opt-in actual CLI/MCP/Python companion: absent/empty/unreconciled projects,
  LF/CRLF payloads, apostrophe paths, all-covered and legacy workflows, explicit
  report/failure prefixes, pre-I/O MCP confirmation, state-summary counts and
  candidate sidecar preservation with isolated strict/permissive subprocess env.

No test may claim execution from source inspection or a skipped process/VM suite.
Run nonzero selected tests at one normally integrated source; reuse existing
same-source CI results rather than repeating audit campaigns:

```sh
cargo check -q --locked --tests
cargo check -q --locked --no-default-features --lib --bin gentle_cli --bin gentle_mcp
cargo test --locked --no-default-features -j1 --lib lifecycle_slice3_ -- --test-threads=1
cargo test --locked --no-default-features -j1 --lib sequence_design -- --test-threads=1
cargo test --locked --no-default-features -j1 --lib genome_edit_intent -- --test-threads=1
cargo test --locked --no-default-features -j1 --lib command_execution -- --test-threads=1
cargo test --locked --no-default-features --features js-interface -j1 --lib lifecycle_slice3_ -- --test-threads=1
cargo test --locked --no-default-features --features lua-interface -j1 --lib lifecycle_slice3_ -- --test-threads=1
cargo test --locked -j1 --test capability_registry_parity --test parity_matrix_freshness
cargo fmt --check
GENTLE_LIFECYCLE_BIN_DIR=/absolute/same-source/bin python3 -m unittest scripts.test_operation_lifecycle_slice3 -v
GENTLE_TUTORIAL_BIN_DIR=/absolute/same-source/bin python3 -m unittest scripts.test_sequence_design_tutorial -v
```

Missing V8/Lua prerequisites are blocked/unrun evidence, not parity. Glen also
owns native GUI dirty/redo checks and Windows file/handle/path acceptance. Current
local checks are Python source/metadata/checkout checks and formatting only:
focused lifecycle/tutorial suite 8 passed, 9 actual-binary skips; full LF/CRLF
checkout suite 75 passed; documentation-status/maturity group 23 passed (includes
11 documentation tests also in the checkout suite). Scoped rustfmt, documentation
status and whitespace checks passed. Session-close: 4 OK, 2 warnings (intentional
dirty tree/manual plan-fidelity reminder), 0 failures. None is Rust/runtime or
native acceptance.

## Explicit Remaining Boundaries

- No captured-payload layer or broad opaque-wrapper migration.
- No global sidecar-loader/implicit-save repair or journal concurrency redesign.
- Slice 3's autosave suppression alone does not secure explicit report paths.
  The separately authorized [Slice 4](operation_lifecycle_slice4_plan.md) now
  authors Plan-only fresh output and bound project/sidecar exclusions, with
  external verification pending. It is not generic exporter migration, overwrite
  permission, hostile parent-race protection or crash durability.
- S2, new editor biology and release/native/scientific acceptance stay separate.
