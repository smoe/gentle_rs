# Claude Review Prompt: Lifecycle Slice 3

Historical prompt prepared 2026-10-26 for Steffen to forward. The owner supplied
Claude's review and authorized the revised Slice 3 on 2026-10-26, explicitly
**without the new payload-resolution layer**. The reconciled
[implementation plan/status](operation_lifecycle_slice3_plan.md) supersedes the
proposal below, including its resolve-once and baseline-before-Apply suggestions.
Rust/runtime/native verification remains pending with Glen/CI. No direct Claude
call was made here; the historical review instructions below are not a new task.

---

Please give a read-only software-contract review of GENtle's next lifecycle
slice: request-aware construction, workflows and implicit project autosave.
The approval defect and static metadata from your previous review are already
implemented. Do not re-plan those fixes, S1 arithmetic, editor chemistry or S2.
Find gaps in the proposed next boundary and return a revised minimal plan.

## Repository And Restrictions

Repo/worktree: `/Users/u005069/.codex/worktrees/47dd/gentle_rs`.
Branch: `gentle_rs_2_main`.
HEAD/local main when inspection began: `320d580376dbcec41f9e1f25d463bd7d6afb60f6`.
During preparation, local main advanced to
`7d7519ad33e2ca932d98f82799e051656f21e513` (CHANGELOG/Glen handoff only).
Runtime inspection predates the owner's subsequent commit/rebase request;
the recorded SHA is a source baseline, not a frozen-HEAD constraint. No fetch
was performed. Read the current files, including working-copy edits if present:

- `docs/operation_lifecycle_slice3_plan.md` (the full original Codex plan).
- `docs/operation_lifecycle_policy_proposal.md`, especially revised slice order
  and the source/feature-scoped executed verification.
- `AGENTS.md`, relevant architecture/decisions/roadmap, protocol/CLI contracts
  and testing portability rules. The roadmap owns release scope, not this prompt.

Use Read/Grep/Glob only. State verifiable source drift; if allowed tools cannot
verify HEAD, say so instead of inventing a Git verdict. Do not silently switch
to the primary checkout. No shell/build/test/model/network calls, edits, fetch,
checkout/rebase, commits/push, CI dispatch, credentials, private sequences,
`data/`, `paper/` or unrelated `output/`. Give feedback in your response, not a
new file. Neither this request nor a positive verdict authorizes implementation.

## Already Landed, Not Current Missing Work

`src/engine/operation_lifecycle.rs::operation_approval_policy` derives Plan's
confirmation from the typed optional path. Every supplied path requires a write-
attempt approval; absence cannot waive explicit candidate/tool gates. Static
Plan raw-op/named-command metadata is External. SaveFile/S1/Apply classifications
and independent raw MCP confirmation remain unchanged.

Codex executed desktop registry 21/21, matrix freshness 1/1, canonical matrix/
catalog/manifest and actual tutorial CLI companion 9/9 at `ab80743e`. At
`98c9331a`, default Cargo check/formatting plus headless sequence-design 32/32
and S1 14/14 passed, including the five new approval regressions, four MCP
preview/apply checks and the unreconciled-state/MCP controls. Those passes are
not current-source, JS/Lua-runtime, native GUI, Windows or full-suite acceptance.
Full source/toolchain/profile/binary identities are in the parent proposal.

## Original Codex Plan To Critique

1. Extend the pure engine-owned typed lifecycle facts only for Plan and S1:
   preserve project content, skip persistent construction/apply repair, preserve
   checkpoints and forbid implicit project save. Keep session journals, op-ID
   reseeding, fresh instance identity, execution bookkeeping and the existing
   approval facet distinct. `None` retains each consumer's legacy behavior,
   including existing MCP RNA file-only exceptions, not "always save".
2. Add explicitly selected request-aware factories for decoded Op/Workflow;
   leave ordinary `from_state`, `new`, project opening, deserialization and
   shared file loading unchanged. Generic JS/Lua direct state-out must preserve
   the supplied decoded state. File-backed baselines are post-existing-loader:
   candidate-sidecar hydration and warning handling must not be silently removed.
3. A non-empty workflow of only covered inspections gets preserve/no-autosave;
   any mutation/unmigrated child or an empty workflow keeps legacy aggregate
   construction/save fallback. Preserve child order, `run_id`, checkpoint and
   journal semantics in both workflow implementations. Cover failures and
   cancellation without adding rollback, error-path autosave or a new transaction.
4. Reuse shared parsing to resolve/capture the effective Op/Workflow before
   construction. Do not inspect `@file` for classification and re-read a different
   payload at dispatch. Resolution I/O stays in execution/worker phases, never
   the pure policy or GUI admission. Preserve approved/raw request and receipt
   identities; use no second command grammar, registry or scientific executor.
   Original wrapper/tool/candidate gates remain additive: resolving a conservative
   raw `op @file` into a no-path Plan must not silently waive its confirmation.
5. Wire raw CLI/MCP, named/shared-shell typed routes and generic JS/Lua direct
   op/workflow coherently. Python remains a thin CLI delegate. Direct CLI intent
   and all MCP/tool/candidate/nested-agent gates stay unchanged. Keep raw MCP
   confirmation ahead of payload parsing/resolution and state loading;
   do not introduce new pre-confirmation reads to classify a file-backed request.
   Keep covered inspection workflow `state_changed=false` typed, not a new whole-project
   serialization in GUI frames; preserve mixed/unmigrated fallbacks.
6. Keep GUI/managed snapshots, project-instance/history/stale/cancellation guards
   and existing journal publication. Do not reuse RNA's no-project observation
   mode for a loaded-sequence preview. Proposed cut line: keep the existing
   guarded session-journal commit without changing persisted state/checkpoints;
   do not advertise "no commit of anything". Challenge this against architecture.
7. Opaque agent-plan and macro/template wrappers remain explicitly legacy unless
   already resolved through the shared boundary. This is not universal adapter
   preservation or permission to unblock S2. Challenge whether this residual
   is acceptable or belongs in the same minimal slice.
8. Add fail-first real state/file/runtime/worker tests, explicitly initialize a
   project before Apply in fixtures that rely on preview-created files, and
   preserve scientific report/DNA/digest/approval/undo assertions. Only after
   owner-approved implementation reconcile contracts/tutorial wording. Alias/
   destination and writer hardening remain separately gated Slice 4; forbidding
   autosave does not prevent an explicit report path overwriting a project.

## Source Anchors Already Inspected

Please verify actual symbols/callers; line numbers are deliberately omitted.

- `src/engine.rs::GentleEngine::{from_state,new,reseed_op_counter_from_state}`:
  normal reconstruction reconciles before dispatch. A new inspection factory
  still needs fresh process-local identity and noncolliding operation IDs.
- `src/engine/ops/operation_handlers.rs::apply_internal`: only S1 skips repairs
  today; Plan writes through the existing handler/writer, including cancelled
  reports. Check post-handler effects as well as the initial preamble.
- `src/engine/state/lineage_containers.rs::{add_lineage_edges,add_container_from_result}`:
  tail helpers return for empty created-sequence results; retain tests for this.
- `src/engine.rs::{operation_checkpoint_kind,apply_with_progress,apply_workflow_with_progress}`
  and `impl Engine::{apply,apply_workflow}`: journals/execution revisions differ
  from project mutation and undo/redo. Do not solve preservation by dropping logs.
- `src/engine.rs::ProjectState::load_from_path` and
  `src/cli_support.rs::load_state_or_default`: sidecar hydration/load diagnostics
  can alter in-memory metadata before any constructor. Keep byte and semantic
  preservation claims separate; loading errors must not become empty-state success.
- `src/bin/gentle_cli.rs::{execute_direct_operation,parse_forwarded_shell_command,run}`:
  raw op's S1 save exemption, raw workflow's unconditional post-success save,
  and parsed outer shell commands whose payload is resolved later.
- `src/engine_shell.rs::{resolve_operation_payload,execute_op_command,execute_workflow_command,execute_parsed_workflow_command}`:
  shared decoding, explicit false flags for Plan/S1 op and serialized workflow
  comparison. `resolve_operation_payload` is crate-visible today; propose a
  minimal public Rust transport boundary, not frontend-owned policy copies.
- `ShellCommand::is_state_mutating` handles typed inline JSON but conservatively
  falls back for file-backed Op payloads. New resolution must not remove an
  existing envelope confirmation gate or enable check/re-read substitution.
- `src/mcp_server.rs::{op_tool_result,workflow_tool_result,run_shell_tool_with_optional_persist}`:
  retain raw confirmation and the existing RNA save exemptions outside this slice.
- `src/js_interface.rs::{apply_operation_impl,apply_workflow,execute_agent_plan_impl}`;
  `src/lua_interface.rs::LuaInterface` registration of `apply_operation`,
  `apply_workflow`, `execute_agent_plan`: some wrappers reconstruct before
  resolving candidate payloads. Distinguish direct generic routes from opaque ones.
- `integrations/python/gentle_py/client.py::GentleClient::{op,workflow,shell}`:
  delegates to CLI. Mocks do not establish actual-process lifecycle parity.
- `src/command_execution.rs::{fork_admitted_command,execute_command,submit_work_with_cancel}`;
  `src/command_execution/rna_artifacts.rs::{admission,freeze,publication}`;
  `src/engine.rs::{clone_without_history,commit_detached_execution}`: generic
  Plan/S1/Workflow currently use project-bound snapshot/publication. A no-project
  observation shortcut can lose loaded targets or stale-owner checks. The current
  commit preserves redo when there was no checkpoint/project mutation.
- `src/mcp_server.rs::tests::assert_sequence_design_mcp_parity`: reads a
  preview-created project before Apply. Preserve false-confirmation and exact
  science checks while separating preview from explicit project initialization.
- `src/engine/analysis/genome_editing/tests.rs::genome_edit_intent_preserves_unreconciled_state_on_success_and_refusal`,
  sequence-design tests, CLI forwarded/shared-shell helpers, optional JS/Lua
  runtime tests, `src/command_execution.rs` stale/cancel/failure tests and
  `scripts/test_sequence_design_tutorial.py`: reuse controls/provenance and
  identify missing ingress coverage instead of broad fixture regeneration.

## Requested Output

Return concise findings ordered by severity, citing file plus symbol and a
concrete counterexample. Distinguish incorrect current-source claims, risks in
the proposed contract and intentionally deferred work. Please answer explicitly:

1. Are the factory/policy and captured-payload boundaries sufficient, minimal
   and correctly ordered? Any other repair before dispatch or after the handler?
2. Is the empty/mixed/unknown workflow fallback conservative without changing
   failure-prefix, progress, cancellation, journaling or save semantics per host?
3. Can current guarded managed journal publication remain, or does the read-only
   worker invariant require a small explicit distinction now? Do not drop
   project snapshots, execution receipts or owner/stale checks to make it pass.
4. Is the opaque wrapper cut line honest enough, or must selected agent/macro
   resolution join this slice? Explain the smallest compatible alternative.
5. Do loader hydration, explicit output aliases or old preview/Apply fixtures
   invalidate any preservation claim or proposed test oracle?

Then give a revised minimal implementation order, required fail-first tests and
remaining owner decisions. No global constructor change, all-readonly migration,
writer redesign, scheduling rewrite or S2 by implication. No execution or
scientific/native/platform acceptance is requested or granted by this review.
