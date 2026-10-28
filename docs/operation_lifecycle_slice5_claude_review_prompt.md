# Claude Review Prompt: Scoped Lifecycle Slice 5

Prepared 2026-10-28 for Steffen to forward manually. **Read-only review, not
implementation authorization.** No Claude invocation or verdict is claimed.
Slice 5 is a proposed follow-up, not an already approved fifth slice. The
previous goal stays paused; Rust/native verification remains with Glen/CI.

The owner has now supplied Claude's source review at `550e6e0a`. Its feedback
and Codex's qualifications are recorded in the [revised minimal plan](operation_lifecycle_slice5_plan.md).
This original prompt remains historical; neither the review nor that plan
authorizes changing DEC-026 or removing managed-preview History rows. Claude
has not reviewed Codex's subsequent qualifications.
The owner subsequently accepted the revised DEC-026/History exception and
authorized implementation on 2026-10-28. The linked plan records authored work
and pending Glen/CI checks; this prompt is historical, not another review claim.

---

Please critique the following narrow lifecycle proposal. Return concise,
source-grounded findings and a revised minimal plan, or recommend deferral if
the benefit does not justify the contract change. Do not implement anything.

## Repository And Restrictions

Repo: `/Users/u005069/.codex/worktrees/faae/gentle_rs`, branch `gentle_rs_3`.
Inspected HEAD and local main:
`550725695437e602c001bd9c881a2f66f861e33c`.
Lifecycle Slice 4 is `1adad235a0aaf43d21f7666b5aa5114b7937d77f`, already an
ancestor. Verify current drift rather than requesting a new audit freeze.

Use Read/Grep/Glob only. No shell, builds, tests, native/VM execution, model or
network calls, file writes, checkout/rebase/commit/push or CI dispatch. Do not
inspect credentials, private data, `outputs/`, `output/`, `data/` or `paper/`.
Return feedback in your answer, not another file. A positive review does not
authorize implementation or resume the previous goal.

Read AGENTS.md, relevant architecture/decisions/roadmap sections, protocol/GUI
contracts and `docs/testing.md`'s cross-platform regression rules. Glen/CI owns
eventual Rust, actual-process, scripting-runtime and native verification. Do not
claim that existing source checks or skipped tests establish runtime acceptance.

## Context Already Delivered

- Slices 1/2 implement request-aware report-write approval and conservative
  discovery. Every present Plan path requires write-attempt approval; absent
  paths do not remove independent candidate/tool requirements.
- Slice 3 preserves project content, dirty revisions and undo/redo for covered
  decoded inspections. Successful operations still journal. Managed generic
  Op/Workflow jobs retain snapshots and guarded project commits. Opaque wrappers,
  ordinary constructors and mixed/empty workflows retain legacy behavior.
- Slice 4 adds Plan-only fresh-file output, project/sidecar exclusions and a
  serde-skipped host project binding. Do not reopen its writer, dependencies,
  Unicode support limits or approval contract in this proposal.
- Slice 3/4's remaining integrated-source Rust/process/VM/native acceptance is
  tracked in `docs/glen_next_verification_handoff.md`, not established here.
- The parent proposal contains historical future-work wording. Use its revised
  slice order and the dedicated Slice 3/4 contracts to distinguish implemented
  work from deliberately deferred boundaries. No S2 editor work is authorized.

## Verified Motivation And Source Anchors

1. `src/engine/operation_lifecycle.rs` separates project preservation from
   approval. A preservation fact is not a no-journal/no-commit permission.
2. `src/engine.rs::fork_detached_execution` retains owner, structural revision
   and journal baseline; `commit_detached_execution` rejects changed journal
   length as well as owner/structural changes. Do not simply weaken that guard.
3. `src/command_execution/lifecycle_tests.rs` explicitly pins current behavior
   in `lifecycle_slice3_managed_inspection_keeps_concurrent_journal_stale_guard`:
   a completed inspection can make an otherwise held managed job stale.
4. `src/command_execution.rs::{submit,submit_parsed,execute_command,
   submit_work_with_cancel}` separates admission, snapshot, execution and
   publication. Ordinary work uses `Publication::Project`; successful guarded
   publication appends the private journal. Cancellation before commit discards
   that delta. Independent RNA publications are named exceptions, not a generic
   read-only precedent. Reuse infrastructure only where its contract fits.
5. `src/command_execution/rna_artifacts.rs::{admission,freeze,publication}`
   already distinguishes inline admission from worker-only payload resolution.
   Do not introduce another payload-resolution layer, early request-file read,
   grammar or generic capability-based exemption.
6. `src/app/routine_and_agent_assistant_ui.rs::{poll_agent_commands,
   cancel_agent_commands_for_project_change}` currently conditions some owner
   handling on `receipt.requires_project_commit`. That is not automatically
   correct for a new snapshot-bound, noncommitting result.
7. `src/engine/analysis/sequence_design.rs` binds the report to its request,
   source DNA/features and exact apply validation. `src/agent_bridge.rs` and
   `src/agent_feedback.rs` distinguish receipts, scientific outcomes and explicit
   user sharing. Receipt hashes do not disclose DNA or authorize application.

Other references: `docs/operation_lifecycle_policy_proposal.md`,
`docs/operation_lifecycle_slice3_plan.md`,
`docs/operation_lifecycle_slice4_plan.md`, DEC-026 and
`docs/tutorial/06-07_synthetic_sequence_design.md`.

## Original Codex Plan To Critique

Proposed name: **standalone no-file preview publication**.

1. Cover only an already decoded, standalone
   `Operation::PlanDnaSequenceDesign { path: None, .. }` through the managed
   command service. Include inline and loaded-sequence targets, and named/raw
   commands when they lower to the same operation. Select from the actual typed
   request, never model flags, names, `state_changed=false` or registry ReadOnly.
2. Keep the existing immutable project snapshot, project admission/owner/source
   checks and execution capacity. This is not project-independent RNA work and
   does not gain an observation-lane slot. Reuse worker resolution after existing
   approval; execute the exact resolved request once. No new payload envelope.
3. Publish the bounded report and source-bound audit through the command result/
   receipt without swapping the live engine or merging the private journal.
   Keep raw engine/CLI/MCP/scripting operation semantics and journals unchanged.
   Managed-preview history becomes a command audit, not a project operation-log
   entry. This is an explicit behavior change requiring owner authorization.
4. Preserve scientific DNA, report/digest fields, search statuses/budgets and
   exact subsequent Apply approval. Bind the managed result to the executed
   request, source snapshot/owner and available source hashes. Never silently
   attribute a stale/replaced-project result to the current project or model
   turn. Do not turn retained historical evidence into automatic Apply or share
   it with the inner agent without the existing user action.
5. Preserve cooperative cancellation and terminal receipt ordering. No-path
   previews have no published file to justify copying RNA's late-cancel rule.
   Specify a small, race-safe cancellation/publication boundary; completion of
   computation is not necessarily successful service publication. A published
   result must not later change scientific status or disappear through relabeling.
6. Leave all other requests on their current paths: `path=Some(...)` including
   empty/invalid values, S1, Apply, workflows even if inspection-only, macros,
   templates, opaque agent plans and arbitrary `submit_work` closures. No queue,
   generic journal merge, loader repair, global exporter migration or S2.
7. Add deterministic synthetic concurrency, stale-owner, cancellation and
   noncoverage tests. Update decisions/architecture and the existing lifecycle
   status docs, with a small tutorial 06.07 explanation of managed previews,
   audit/history and explicit Apply. Record authored tests as unrun here.

## Questions That Need A Decision, Not An Assumption

- Is retaining a per-command audit instead of the live journal the smallest
  useful change? Which existing history/agent/result consumer would regress?
  Is deferring this preferable to inventing a general journal-merge scheme?
- Where can the existing worker freeze select this exact standalone operation
  without duplicate resolution, premature `@file` I/O or parser/approval drift?
  Check both `submit` and `submit_parsed`, not just one GUI call site.
- Separate **needs a project snapshot**, **requires a live project commit** and
  **result belongs to a project/turn**. Can existing fields express this safely,
  or is one small additive distinction necessary? Do not mark needs_project=false
  just to avoid a commit. Do not drop result-owner checks as a side effect.
- Recommend exact source-freshness behavior: reject before result publication
  when source/project changed, or retain explicitly historical evidence. Explain
  how changes after publication but before consumption/sharing are detected.
  Include loaded-source annotations, project replacement, display-only changes
  and output-ID collisions; avoid whole-project serialization on the UI thread.
- Specify cancellation before snapshot, during search, after a computed result
  but before publication, and after publication. Distinguish the report's search
  status from the service receipt. Do not infer success from zero committed ops.
- What minimal receipt/result binding is necessary without changing scientific
  schemas or treating a content digest as authenticated search provenance?
  Retain bounded output/refusal, consumption and eviction behavior.

## Proposed Fail-First Acceptance Matrix

- A no-path preview finishes while a held genuine mutation remains current;
  two no-path previews finish without invalidating one another. Per-command
  audits remain inspectable while live journal/state/revisions/undo/redo stay
  unchanged by preview publication.
- Actual structural/source edits and project replacement are not made safe by
  removing the journal append. Test the chosen stale-result rule at execution,
  publication and delayed GUI/agent consumption. Ordinary viewport/display edits
  must not become biological edits merely to simplify the test.
- Cover named/raw forms, existing worker-only `@file` inputs, inline/loaded
  targets, malformed requests and result-size refusal. Approval/candidate/raw MCP
  requirements and nested-agent guards remain intact.
- Use barriers/channels, not sleeps, for cancellation windows and delayed
  consumption. Retain exact search status/DNA/digests and explicit Apply checks.
- With-path previews, S1, covered/mixed workflows, opaque wrappers and trusted
  arbitrary work remain conservative controls. Preserve their journal and stale
  guards; do not merely replace the existing stale-guard regression with a pass.
- Exercise host/UI and inner-agent receipt/result handling as authored
  deterministic tests, not a live-model call. Scientific parity with synchronous
  raw adapters does not imply identical host-local audit placement.

Use existing synthetic tutorial fixtures and service/agent test helpers. Glen/CI
will run focused command-execution/lifecycle/sequence-design/agent tests, locked
default/headless test-target checks, applicable parity/actual-process checks and
native consumption acceptance at the integrated source. Give precise recommended
commands if source review changes that inventory; run none yourself.

## Requested Answer

Return (1) verified current facts and any drift; (2) severity-ordered concerns
with source symbols; (3) a revised minimal plan, including the exact audit,
staleness and cancellation decisions; (4) missing tests and Glen/CI checks;
(5) a clear implement/revise/defer recommendation. Distinguish facts from
inferences and proposed contract changes. Do not expand the scope to close all
lifecycle boundaries, and do not claim implementation or runtime acceptance.
