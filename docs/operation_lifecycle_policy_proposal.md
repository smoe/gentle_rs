# Request-Aware Operation Lifecycle Policy

Status: proposal for owner-forwarded Claude review, 2026-10-25. Not implemented,
not an adopted decision, and not authorization for editing S2. Review prompt:
[read-only Claude request](operation_lifecycle_claude_review_prompt_20261025.md).

## Source And Scope

Repository: `/Users/u005069/.codex/worktrees/47dd/gentle_rs`.
Branch: `gentle_rs_2_main`; local main and preparation HEAD:
`81b675048ab0c886185d55e101551790d8c3bea5`. The clean tracked branch was rebased
onto local main before this inspection; no fetch was performed. These documents
were prepared as working-copy drafts; committing them does not adopt the policy.
Untracked `output/` is unrelated and untouched.

Start with `PlanDnaSequenceDesign`, distinguish its absent/present output path,
and use `EvaluateGenomeEditIntent`, `ApplyDnaSequenceDesign` and `SaveFile` as
controls. Do not migrate every operation or infer policy from names, prefixes,
the registry's coarse `CapabilityMutation`, or model-supplied declarations.
The [S1 review reconciliation](base_prime_editing_preparation.md#landed-s1-review-reconciliation-2026-10-25)
and [existing safety decision](decisions.md#read-only-genome-edit-intent) explain
why a single read-only boolean is insufficient.

## Observed Behavior, Not Desired Guarantees

This is source inspection, not an executed regression or security exploit.

| Boundary | Current source fact |
| --- | --- |
| Shared solver | `engine/analysis/sequence_design.rs::plan_dna_sequence_design` takes `&self` and returns a transient report. Apply is separate. |
| Shared parser | `engine_shell/command_parsers.rs::parse_sequence_design_command` lowers both named commands to `ShellCommand::Op`; `--path` becomes `PlanDnaSequenceDesign.path`. There is not a second named-command executor. |
| Registry | `gentle-protocol::infer_engine_operation_mutation` falls through to `Mutating` for `PlanDnaSequenceDesign`; `SaveFile` is explicitly `ReadOnly`. Neither label describes all lifecycle effects. |
| Detailed introspection | `engine_shell.rs::sequence_design_capability_descriptor` says planning is nonmutating and needs no confirmation, while also advertising an optional external artifact effect. |
| Parsed-command confirmation | `ShellCommand::is_state_mutating` exempts every `PlanDnaSequenceDesign`, irrespective of `path`. Stored-plan execution relies on this policy in addition to candidate flags. |
| Apply preamble | `apply_internal` reconciles lineage/containers before planning, including before a refused request. Only `EvaluateGenomeEditIntent` currently skips these repairs. |
| History | `operation_checkpoint_kind` excludes planning, but successful `apply`/`apply_with_progress` still append an operation receipt and advance execution bookkeeping. No undo checkpoint does not mean no journal activity. |
| Shared-shell result | `execute_parsed_operation_command` reports `state_changed=false` for planning without serializing the project. On unreconciled state this can hide the preamble's repairs. |
| Explicit report writer | The plan handler calls `write_pretty_json_file` when `path` is present. That helper creates parents and uses `fs::write`; it is not atomic/no-clobber publication. Cancellation is a report status, so this handler can also write a cancelled report. |
| CLI | Forwarded/shared-shell commands save on `state_changed`; direct `gentle_cli op` instead saves after success, except for the S1 exemption. Thus the planning spellings have different implicit project-file behavior. |
| MCP | Raw `op` requires `confirm=true` independently, then saves after success except for S1. `sequence_design_mcp_confirmation_preview_and_exact_apply_parity` currently relies on the plan-created state file. |
| Construction/workflows | `GentleEngine::from_state` reconciles before dispatch; skipping only the apply preamble is not sufficient to preserve the loaded in-memory state. Direct CLI/MCP workflow paths also save after success. |
| Scripting | `js_interface.rs::apply_operation_impl`/`apply_workflow` and Lua's `apply_operation`/`apply_workflow` also construct with `from_state` and return the resulting state. Python's `GentleClient.op`/`workflow` delegate to CLI. They need the same covered-request construction policy, not a separate script implementation. |

The existing shell parity test explicitly expects a plan with `--path` to have
`is_state_mutating()==false`. Preserve its scientific assertions, but do not
bless that boolean as proof that an agent may write without confirmation.
The new S1 unreconciled-state regression in `ced03a22` remains relevant and
requires actual Glen/CI execution, not another source-only acceptance claim.

## Proposed Policy Dimensions

Use one request-aware, engine-owned policy selector, tentatively
`operation_lifecycle(&Operation) -> Option<OperationLifecyclePolicy>`, in a small
headless helper under `src/engine/`. Names below are design vocabulary, not new
serialized API fields. It must examine the typed request, not serialize the
whole operation just to recover its variant name. No file/network I/O occurs
during policy selection.
`Operation` currently lives in the root engine. Keep the selector headless;
do not introduce a dependency from `gentle-protocol` back to the root crate or
from the engine to GUI/CLI/MCP modules. The protocol's static family projection
must be explicitly aligned and cross-checked, not call a root payload selector.

| Dimension | Meaning |
| --- | --- |
| Project effect | Preserve persisted project content, or allow the existing mutation handler. This is distinct from external files. |
| Explicit artifact effect | No outputs, or the caller-selected report destination. Preserve actual paths; a possible effect is not successful publication. |
| Reconciliation | Skip persistent lineage/container repair, or use the existing legacy/mutation setup. |
| History | Preserve undo/redo, or use the existing full checkpoint. Keep operation journals and execution counters as separately disclosed session bookkeeping. |
| Implicit project persistence | Forbid adapter autosave, or retain the host's existing mutation-save behavior. This is not an instruction to write a requested report. |
| Conversation approval | No operation-derived approval, or approval required for the exact artifact/mutation intent. Adapter/tool and candidate gates remain additive. |

Return `None` for an unmigrated operation: each caller retains its existing
conservative behavior. `None` must never mean no effects, no confirmation or no
autosave. Do not change generic `from_state`, all checkpoint rules or all export
permissions merely to remove match lists. Existing mutation/stale-owner guards,
parsed-command confirmation and the nested-agent ban remain intact.

### Recommended Target For The First Migration

| Typed request | Project content | Explicit artifact | Reconcile | Undo/redo | Implicit project save | Operation-derived conversation approval |
| --- | --- | --- | --- | --- | --- | --- |
| `EvaluateGenomeEditIntent` | Preserve | None | Skip | Preserve | Forbid | None additional; preserve raw MCP tool gate |
| Plan, `path=None` | Preserve | None | Skip | Preserve | Forbid | None additional |
| Plan, `path=Some(...)` | Preserve | Requested preview only | Skip | Preserve | Forbid | Required for exact report-write intent |
| `ApplyDnaSequenceDesign` | Existing guarded mutation | No added artifact | Existing | Existing full checkpoint | Existing host behavior | Existing confirmation plus exact approval digest |
| `SaveFile` and other operations | Unmigrated | Unmigrated | Unmigrated | Unmigrated | Unmigrated | Existing policy; no new exemption |

Empty, malformed or unsafe supplied paths are not equivalent to `None`: retain
the write classification and refuse them before execution. Loading `@request`
remains an explicit input read, not a filesystem-write permission. The policy
does not authorize an implicit reference fetch, sequence creation or S2 design.

## Concrete Consumer Rules

1. **Engine construction and apply:** add an explicitly selected inspection
   construction path that retains project content without persistent repairs;
   leave ordinary `from_state` unchanged. Preserve the necessary session ID,
   history defaults and op-counter reseeding. Both engine apply paths and
   workflows must consult the same reconciliation/history facts. The pure plan
   reads sequence/annotation data, not repaired container/lineage indexes.
2. **Shell and agents:** use the request policy both for confirmation and for
   `state_changed`. Keep `is_state_mutating` as a compatibility wrapper if needed,
   but document that the existing boolean also enforces external-write approval.
   A forged candidate with both flags false cannot run a plan with an output
   path without confirmation. Named command and generic `op` must agree.
3. **CLI/MCP persistence:** suppress implicit project-file creation/rewrite for
   both plan forms, not the explicitly requested report. Preserve direct CLI
   invocation as explicit user intent without inventing a universal `--confirm`
   flag; delegated execution still uses conversational approval. Raw MCP `op`
   and `workflow` retain `confirm=true`, including a no-path preview and S1.
4. **Workflows and hosts:** aggregate facts, not operation names. A workflow
   containing only covered inspections must preserve project content and its
   file. One artifact requires artifact approval; one mutation/unmigrated op
   retains the existing conservative host behavior. Keep workflow tool gates,
   execution order, partial-side-effect reporting, progress and cancellation.
   Do not bypass the managed-command stale-project/admission safeguards.
   Cover the generic JS/Lua state-in/state-out constructors and Python's CLI
   routes in that inventory; an unchanged disk file does not prove that a
   scripting response preserved unreconciled project content.
5. **Metadata:** recommend projecting Plan's durable-project effect as
   `ReadOnly`, with explicit conditional artifact details; the coarse label is
   never an authorization source. A static family descriptor should conservatively
   disclose that report output can require approval; a concrete parsed no-path
   request can have no operation-derived approval. Explicit candidate/tool
   requirements are never removed. Reconcile canonical registry, detailed
   introspection, effects and generated projections together; ask Claude to
   challenge this compatibility projection before adopting it.

### Output And Cancellation Boundaries Requiring Review

Forbid *implicit* project autosave does not prevent a caller from supplying
the project filename as the explicit report path. Proposed rule: reject a
preview destination aliasing the selected project file before writing, including
existing aliases supported by the repository's cross-platform path helpers.
Inspect those helpers before choosing the implementation. Do not promise
symlink/hard-link or race safety without an executed design/test path. If this
cannot be enforced narrowly, defer the with-path migration rather than claim
project-file preservation that it cannot provide.

Do not globally replace `write_pretty_json_file` or silently change report
overwrite semantics. A stronger fresh/atomic/no-clobber contract is a separate
owner decision unless required for the approved collision rule. Path validation
must precede normalization; use the existing shared quoting and portability
rules, including spaced/Unicode/apostrophe paths and native Windows cases.

Preserve report bytes and outcome semantics: a search result of
`status=cancelled` with a supplied path is different from cancellation before
dispatch and may retain the existing explicitly requested report. No code may
call either case a successfully published artifact without an actual write.
Failure before dispatch writes nothing; do not introduce workflow rollback of
already written files or infer that undo removes them.

## Minimal Implementation Order After Review And Approval

1. **Characterize and add failing boundary tests.** Pin the scientific/report
   bytes and current journal/undo contract; expose the target permission,
   unreconciled-state and autosave mismatches without weakening expectations.
   Inventory the actual managed/scripting call sites before selecting them.
2. **Migrate the covered requests coherently.** Add the scoped policy selector
   and wire construction, reconciliation, checkpoints, shell confirmation,
   state-change reporting and CLI/MCP persistence together. Reuse it in covered
   workflow/host paths; leave all other operations on explicit legacy fallback.
   Do not ship an intermediate change that suppresses approval but still writes.
3. **Align metadata and user contracts.** Update the canonical classification,
   detailed discovery and affected generated artifacts with LF/CRLF coverage.
   Correct tests that assume MCP planning creates a project file: explicitly
   initialize a project for subsequent apply instead of loosening apply checks.
   Explain the new artifact approval and autosave behavior in existing CLI/
   protocol/tutorial 06.07 text. No new active editing tutorial or release gate.

Keep S1 regressions as a control, not a reason to repeat its scientific review.
Do not expand this to every read-only operation before evidence warrants it.
S2 still needs separate implementation authorization after this prerequisite.

## Required Deterministic Verification

| Boundary | Required expectation |
| --- | --- |
| Policy admission | Exact typed covered variants and both path forms; malformed paths stay writes until refused; an unknown/legacy operation cannot gain an exemption from its name or registry label. |
| State and refusal | Inline and loaded input, both normal and deliberately unreconciled projects, success/refusal/cancellation: project bytes, mutation/structural revisions and undo/redo unchanged for plan. Independent repair controls must actually change their fixtures. |
| Session bookkeeping | Successful operations retain their existing journal/op-ID/execution receipt contract; refused operations do not acquire successful receipts. No undo checkpoint or redo clearing for preview. |
| Files | Missing state stays missing; noncanonical existing state JSON stays byte-identical through direct CLI, named/shared shell and raw MCP. Requested report content equals the same direct-engine report. Protect project/report alias collisions and preserve sentinel files on refusal. |
| Approval | Unconfirmed and falsely read-only agent candidates cannot write reports or mutate via apply, whether named, `op`, workflow or macro-wrapped. No-path requests retain existing applicable candidate/tool gates; MCP raw tools still require confirmation. |
| Workflow/host | All-inspection, mixed artifact and mixed mutation cases; no nested apply/checkpoint; progress/cancel/stale-owner behavior preserved, and earlier file effects are disclosed on later failure. |
| Metadata/fallback | Canonical/detailed metadata agrees on Plan's project effect and conditional artifacts; unrelated operations including `SaveFile` keep their previous behavior. |
| Scientific identity | Same solver strategy/budgets, hashes, statuses, search claims, exact apply approval/receipt, lineage/container derivation and undo. No preview report migration or new solver run at apply. |
| Portability | LF/CRLF request and generated-projection checks; shared quoted paths; native Windows file handles/alias cases and no source-file rewrite on refusal. |

Reuse `genome_edit_intent`, `sequence_design`, `sequence_design_*_parity`,
`capability_registry_parity`, `parity_matrix_freshness`, the relevant agent/
command-execution tests and canonical generator checks; select exact commands
against the eventual implementation. Include normal locked Cargo check and
`scripts.test_sequence_design_tutorial`, `scripts.test_tutorial_checkouts`.
Current no-local-Rust restrictions remain in force; Glen/CI runs belong in the
existing consolidated handoff only after an implementation exists. No runtime,
native, scientific or performance acceptance is claimed for this proposal.

## Questions For Claude And The Owner

- Is the scope small enough while covering the constructor, workflow and raw
  adapter save paths that could otherwise undermine a strict preview?
- Is `ReadOnly` plus conditional artifact evidence the right compatibility
  projection, or should the coarse family be `External` until metadata can
  express those dimensions? Neither answer may weaken request-time approval.
- Can the project/report collision rule reuse a sound existing helper without
  expanding this into a generic writer rewrite? Identify limits explicitly.
- Are journal/cancel semantics and static-versus-concrete confirmation stated
  precisely enough to preserve existing host and scripting behavior?

These are review questions, not accepted invariants. After Claude's source-
grounded reply, reconcile it with the original proposal, revise the minimal
plan and obtain implementation authorization before changing runtime behavior.
