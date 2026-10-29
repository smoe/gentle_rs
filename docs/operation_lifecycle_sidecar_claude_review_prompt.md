# Claude Review Prompt: Degraded Candidate-Sidecar Load/Save Safety

Prepared 2026-10-29 for Steffen to forward manually. **Read-only review, not
implementation authorization.** No Claude invocation or verdict is claimed.
The previous goal stays paused. Rust, actual-process/VM and native verification
remain with Glen/CI; no push or new release gate.

---

Please critique the proposed next lifecycle slice below. Verify the current
source, distinguish facts from proposed compatibility changes, and recommend
the smallest useful repair. Do not implement anything.

## Repository And Restrictions

Repository: `/Users/u005069/.codex/worktrees/faae/gentle_rs`.
Branch: `gentle_rs_3`. Inspected HEAD:
`186937bd926449930b69a5d08d7bb33fcf5174d6`.
Recheck drift through source inspection; this SHA is provenance, not an audit
freeze or authority to reset/rebase another checkout.

Use Read/Grep/Glob only. No shell, builds, tests, native/VM execution, network or
provider calls, file writes, commit/rebase/push or CI dispatch. Do not inspect
credentials, private sequencing/project data, `output/`, `outputs/`, `data/` or
`paper/`. Read tracked source and explicit synthetic test definitions only.
Return feedback in your answer, not a new file. A positive review does not
authorize implementation or resume the previous goal.

Read AGENTS.md and relevant sections of `docs/architecture.md`,
`docs/decisions.md`, `docs/roadmap.md`, `docs/protocol.md`, `docs/gui.md`,
`docs/cli.md`, `docs/CHANGELOG.md` and the cross-platform rules in
`docs/testing.md`.

## Existing Lifecycle Context

Slices 1--5 implement bounded approval/discovery, decoded inspection
construction/no-autosave, Plan-only fresh report writing, and snapshot-bound
managed no-path previews. Those contracts are not being reopened here.

`docs/operation_lifecycle_slice3_plan.md`, under Adapter Boundaries and Explicit
Remaining Boundaries, explicitly defers permissive candidate-sidecar loss on
unmigrated saves. Its preservation promise starts after the loader. The parent
`docs/operation_lifecycle_policy_proposal.md` and roadmap retain loader semantics
as separately gated work. The proposed slice is not an already approved Slice 6.

Commit `e39469bd` records a scoped rehearsal at binary-reported source
`284e5b70`: 25 source/process checks passed, plus limited macOS GUI observations.
Native cancellation/concurrency remained inconclusive. See
`docs/glen_next_verification_handoff.md#synthetic-lifecycle-rehearsal-at-284e5b70-2026-10-28`.
Neither that evidence nor the owner's separate local CI report proves the
degraded-load/save boundary safe. Do not require repeating unrelated acceptance.

## Source Evidence Already Checked

These are static observations, not an executed loss reproduction:

1. `src/engine.rs::ProjectState::load_from_path` parses the project and calls
   `hydrate_candidate_store_from_external_ref`. On failure, strict mode returns
   the error. Otherwise it removes `metadata["candidate_sets"]` and writes
   `candidate_sets_load_warning` containing the error message and strict-env
   name. Successful hydration removes that warning.
2. `strict_candidate_store_load_enabled` implements the existing opt-in
   `GENTLE_CANDIDATE_STORE_STRICT_LOAD`. Do not assume a default-policy change
   is a neutral implementation detail.
3. `ProjectState::save_to_path` clones state, prepares a candidate-store
   transaction, commits that sidecar transaction, then atomically writes the
   project, with existing rollback/finalization handling.
4. `prepare_candidate_store_sidecar_transaction` chooses Remove when candidate
   metadata is absent but the derived save-owned sidecar directory exists.
   Thus saving a degraded state can persist the lost reference and remove that
   directory. This is a source-derived risk, not a claim every missing-sidecar
   load deletes files.
5. `read_candidate_store_from_metadata` deserializes `CandidateStore` or falls
   back to its default. Merely retaining an external reference may therefore
   still present an empty store. Inspect candidate readers and writers, not only
   the load function, before proposing that one-line fix.
6. `src/cli_support.rs::load_state_or_default` preserves its missing/empty-file
   defaults and delegates real loads to `ProjectState::load_from_path`.
   GUI `save_project_to_file`, CLI, shared-shell LoadProject/SaveProject, MCP,
   and JS/Lua file load/save bindings use the shared persistence methods.
   Verify the full caller inventory and identify any bypass; Python delegates
   to CLI, not a separate persistence implementation.
7. Existing tests in `src/engine/tests.rs` explicitly require permissive load
   to drop the reference and record a warning, and strict load to fail:
   `test_project_state_load_degrades_when_candidate_sidecar_is_missing` and
   `test_project_state_load_strict_mode_errors_when_candidate_sidecar_is_missing`.
   Nearby tests cover healthy round-trips and sidecar replacement/removal.

Additional source anchors: `CandidateStoreReference`,
`load_candidate_store_from_ref`, `load_candidate_store_for_externalization`,
`read_candidate_store`, `write_candidate_store`,
`commit_candidate_store_transaction`, `rollback_candidate_store_transaction`
and `finalize_candidate_store_transaction` in `src/engine.rs`;
`src/app.rs`, `src/bin/gentle_cli.rs`, `src/mcp_server.rs`,
`src/engine_shell.rs`, `src/js_interface.rs`, `src/lua_interface.rs`;
`scripts/test_operation_lifecycle_slice3.py` and
`docs/tutorial/06-07_synthetic_sequence_design.md`.

## Original Codex Plan To Critique

Proposed scope: **protect saving after a degraded candidate-store load**.

1. Retain the unresolved original candidate-store reference and a structured
   diagnostic. Unavailable candidates must not be represented as known-empty.
2. Keep unaffected project inspection available. Candidate-dependent operations
   must not infer biological absence or silently replace unresolved content.
3. Refuse unsafe explicit saves and autosaves through shared persistence, before
   any project/sidecar/staging/parent-directory write or removal, while the load
   remains unresolved. No adapter-local save exception.
4. Add deterministic synthetic load -> inspect -> save regressions for missing
   and malformed sidecars, healthy controls, legacy cases, LF/CRLF and quoted
   paths. Preserve healthy atomic replacement/rollback and lifecycle controls.
5. Document refusal and the safe recovery boundary in tutorial 06.07 and the
   relevant interfaces/decisions. Preserve the distinction between inspection,
   in-memory mutation, persisted state and already published external artifacts.

This is an initial proposal, not a chosen DTO or blanket ban on all project
work. Challenge the scope if retaining references plus blocking saves requires
an unnecessarily broad redesign. No new generic payload-resolution layer.

## Questions Requiring Review

- What is the smallest safe representation? Compare retaining the exact raw
  reference plus diagnostic, guarding the existing degraded state, and a stricter
  write-capable load boundary. State compatibility costs and recommend one;
  do not introduce a generic loader framework by default.
- Where should unresolved status live so clone, detached execution, history,
  serde state-in/state-out and reload cannot accidentally clear it? An
  unpersisted flag can be lost; warning metadata alone is not authenticated
  evidence. Define the supported threat model without promising protection
  against arbitrary caller edits to project JSON.
- How should legacy projects already containing a warning but no recoverable
  reference behave? Do not invent the original reference or clear the warning
  merely because no hydration was attempted. Distinguish ordinary known-empty
  stores and intentional candidate-set deletion from unresolved prior content.
- Which candidate accessors/mutators need a structured unavailable result or
  refusal? Does fixing unknown-versus-empty require a small additive inspection
  result, or does it force a separate slice? Avoid a new registry/catalog system.
- Cover Save, autosave, Save As and explicit project-copy/export routes. A
  relative source reference must not be reinterpreted relative to the new
  destination. Do not retry hydration against a different project's sidecar or
  silently copy/discard unresolved data. Separate the source binding from the
  destination before claiming preservation.
- If a command mutates memory and its implicit save is refused, what must its
  result say? Do not claim rollback of in-memory work or external files merely
  because persisted project bytes stayed unchanged. Check GUI and agent feedback
  as well as CLI exit status and MCP errors.
- Recommend the smallest usable recovery policy. Restoring the source sidecar
  and reloading may suffice; automatic repair, force-save, discard confirmation
  and a recovery GUI are not authorized. If safe editing needs more, identify
  that owner decision explicitly rather than silently adding commands.
- Identify existing schema/digest/constructor assumptions affected by any new
  field. Healthy legacy payloads must not gain serialized defaults or altered
  scientific hashes without an explicit reason. Preserve strict-load semantics,
  raw-tool/candidate approval and the existing nested-agent ban.

## Proposed Acceptance Boundaries For Glen/CI

- Actual missing/malformed index and record failures: unaffected sequence
  inspection works; unavailable candidate data never becomes an asserted empty
  result. Required strict-load refusal remains a control.
- Attempt same-path Save, autosave, Save As and applicable project-copy routes
  after degraded load. Assert original project bytes and all sentinel sidecar
  bytes survive; fresh destinations and parents remain absent on refusal.
- Healthy inline/external and genuinely empty stores retain existing save,
  intentional deletion, replacement and rollback behavior. Include a negative
  control proving the preservation oracle detects a destructive save.
- Cover clones, checkpoints/undo/redo and scripting serialization where relevant;
  old warning-only state must not become silently clean. Restore valid source
  files and reload to exercise the explicitly chosen recovery path.
- Exercise shared engine/GUI helper, CLI/shared shell, MCP and JS/Lua file
  bindings, plus Python delegation and deterministic agent error feedback as
  applicable. Do not call a live inner model or author parallel persistence.
- Use isolated strict/permissive environment fixtures, LF/CRLF text cases,
  shared shell quoting, spaces/apostrophes/non-ASCII names and native Windows
  handle/path checks. Never weaken hashes or normalize contractual raw bytes.
- Retain Plan/S1 no-autosave, Apply approval, fresh report writer, snapshot/owner,
  scientific-output and unrelated legacy controls. New saved fixtures need
  provenance; prefer inline hand-crafted synthetic cases.

Recommend exact nonzero-count targeted Rust/VM/actual-process filters and locked
default/headless test-target checks for Glen/CI; run none. Distinguish authored
tests from executed results and platform acceptance. No full-suite or native
verdict is inherited from source inspection.

## Non-Goals And Requested Answer

No implementation, automatic recovery/discard, strict-default switch, generic
payload resolver, opaque-wrapper migration, queue, journal merge, exporter
rewrite, S2 editor biology, production-data replay or release-gate expansion.
Preserve unrelated work and the paused goal; no push.

Return: (1) verified facts/drift; (2) severity-ordered concerns with source
symbols; (3) revised minimal plan and explicit compatibility/recovery decisions;
(4) tests and Glen/CI commands; (5) implement/revise/defer recommendation.
Call out any proposal larger than the evidence warrants. Implementation requires
separate owner authorization after this review.
