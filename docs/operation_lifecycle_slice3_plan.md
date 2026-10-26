# Lifecycle Slice 3: Construction, Workflows And Implicit Persistence

Status, 2026-10-26: proposed software contract, **not implemented and not yet
Claude-reviewed**. The owner requested this plan and a
[forwardable review prompt](operation_lifecycle_slice3_claude_review_prompt.md),
not implementation. Prior approval of slices 1--2 does not authorize this slice.
The [parent proposal](operation_lifecycle_policy_proposal.md#revised-slice-order)
and [roadmap](roadmap.md) retain scope/order authority; S2 remains separate.

## Baseline And Scope

- Worktree: `/Users/u005069/.codex/worktrees/47dd/gentle_rs`.
- Branch: `gentle_rs_2_main`; HEAD and local main when inspection began:
  `320d580376dbcec41f9e1f25d463bd7d6afb60f6`. No fetch was performed.
- During preparation, local main advanced to
  `7d7519ad33e2ca932d98f82799e051656f21e513` in CHANGELOG/Glen-handoff text only.
  Runtime inspection predates the owner's subsequent commit/rebase request.
  These SHAs record inspection provenance, not a frozen-HEAD requirement;
  check actual source drift when reviewing the committed/rebased documents.
- Tracked tree was clean before these planning edits; unrelated `output/` is
  untouched. Private data, the primary checkout's untracked prompts and `paper/`
  are not inputs to this plan.
- Slices 1--2 already provide Plan's request-aware write-attempt approval and
  External family metadata. Preserve their
  [source-scoped executed evidence](operation_lifecycle_policy_proposal.md#executed-local-verification-2026-10-26),
  including headless 32/32 sequence-design and 14/14 S1 at `98c9331a` and earlier
  desktop parity checks at `ab80743e`. Those are not new-source execution here.
- Cover `PlanDnaSequenceDesign` with either path form and
  `EvaluateGenomeEditIntent` (S1). Do not migrate all nominally read-only ops.
  `ApplyDnaSequenceDesign`, `SaveFile` and unknown/unmigrated operations are
  compatibility controls, not newly exempt operations.

Goal: a covered inspection must not repair persisted project content merely
because a host reconstructs an engine, runs a workflow, or saves after success.
Keep explicit requested artifacts distinct from implicit project persistence.

## Verified Source Map

Use file plus symbol rather than historical line numbers. These are inspected
facts at the baseline above, not claims that the proposed fixes exist.

| Source and anchor | Current boundary |
| --- | --- |
| `src/engine/operation_lifecycle.rs::operation_approval_policy` | Typed Plan approval only; does not control construction, history or saving. |
| `src/engine.rs::GentleEngine::{new,from_state,reseed_op_counter_from_state}` | `from_state` sets session/history defaults, repairs lineage/containers and reseeds counters. Ordinary `new` does not perform those repairs. |
| `src/engine/ops/operation_handlers.rs::apply_internal` | Only S1 skips pre-dispatch repair today. Plan can repair even on refusal. Its handler can write a cancelled report through the existing writer. |
| `src/engine/state/lineage_containers.rs::{add_lineage_edges,add_container_from_result}` | Post-handler helpers return without creating lineage/containers when no sequence was created, as for Plan/S1. Keep these controls rather than assume all post-handler code is inert. |
| `src/engine.rs::{operation_checkpoint_kind,apply_with_progress,apply_workflow_with_progress}` and `impl Engine` | Plan/S1 have no undo checkpoint; successful operations still append journals and advance execution bookkeeping. Both workflow implementations retain per-step semantics and `run_id`. |
| `src/engine.rs::ProjectState::load_from_path`; `src/cli_support.rs::load_state_or_default` | Loading can hydrate candidate sidecars and change load-warning metadata before engine construction. Missing/empty files yield default state; existing JSON is not rewritten by loading. |
| `src/bin/gentle_cli.rs::{execute_direct_operation,parse_forwarded_shell_command,run}` | Raw op suppresses saving only for S1. Raw workflow saves after success. Named/shared-shell paths already parse an outer `ShellCommand`, construct normally, then use `state_changed`; Op/Workflow payloads may still require resolution. |
| `src/engine_shell.rs::{resolve_operation_payload,execute_op_command,execute_workflow_command,execute_parsed_workflow_command}` | Shared payload decoding/execution already exists. Plan/S1 raw-op results use explicit false project-change flags; workflow uses before/after state serialization. |
| `src/mcp_server.rs::{op_tool_result,workflow_tool_result,run_shell_tool_with_optional_persist}` | Raw tools retain mandatory confirmation. Op's no-save list includes S1 and named RNA file-only operations; workflow saves after success. Shell construction happens before effective payload resolution. |
| `src/js_interface.rs::{apply_operation_impl,apply_workflow}`; `src/lua_interface.rs::LuaInterface` registration of `apply_operation`/`apply_workflow` | Generic state-in/state-out routes reconstruct normally and can return repaired state despite a covered operation. |
| `integrations/python/gentle_py/client.py::GentleClient::{op,workflow,shell}` | Thin CLI delegation, not a second engine or persistence implementation. |
| `src/command_execution.rs::{fork_admitted_command,execute_command,submit_work_with_cancel}`; `src/command_execution/rna_artifacts.rs::{admission,freeze,publication}` | Generic Plan/S1/Workflow currently use the project-bound detached path. A false `state_changed` alone does not prevent its guarded commit/journal publication. Named RNA observations instead may execute without a project snapshot. |
| `src/js_interface.rs::execute_agent_plan_impl`; Lua `execute_agent_plan`; CLI/MCP generic agent/macro wrappers | Some reconstruct an engine before resolving the selected candidate or expansion. Do not advertise universal wrapper preservation from fixing direct Op/Workflow. |
| `src/mcp_server.rs::tests::assert_sequence_design_mcp_parity` | Explicitly reads a project file currently created by preview before testing Apply. This fixture must change deliberately when preview stops autosaving. |

## Proposed Contract

Keep a pure, typed, engine-owned lifecycle selector beside the existing approval
facet. Suggested internal facts are project-content preservation, pre-dispatch
reconciliation, checkpoint policy and implicit-save eligibility. Names are
provisional Rust design vocabulary, not new serialized schemas or endpoints.
Existing approval behavior must remain equivalent, including additive host gates.

| Request | Construction/apply repair | Undo/redo | Implicit project save | Explicit artifact and approval |
| --- | --- | --- | --- | --- |
| S1 | Skip persistent repair | Preserve checkpoints | Forbid | None; preserve independent MCP/tool confirmation. |
| Plan, absent path | Skip persistent repair | Preserve checkpoints | Forbid | None; do not remove explicit candidate/tool requirements. |
| Plan, any present path | Skip persistent repair | Preserve checkpoints | Forbid | Existing report write attempt and approval, including bad/empty paths and cancelled reports. |
| Apply, SaveFile, any unmigrated operation | Existing behavior | Existing behavior | Existing host-specific behavior | Existing confirmation, digest, effects and writer. |

`None` from the selector means **use that consumer's existing fallback**, not
"always save" or "no effects". In particular, do not erase MCP's existing RNA
file-only save exemptions. Do not infer lifecycle from registry `ReadOnly`,
`state_changed`, a command prefix or model-provided flags.

### Define What Is Preserved

- Direct state-in/state-out: compare the supplied decoded `ProjectState` with
  returned state semantically, including unreconciled lineage, containers and
  metadata. This is not byte identity of arbitrarily formatted input JSON.
- File-backed routes: absent files remain absent; empty/noncanonical/LF/CRLF
  existing project bytes remain unchanged. No implicit parent/sidecar writes.
  The engine baseline is the state returned by the existing shared loader;
  sidecar hydration and its diagnostics are not disabled or relabelled repair.
- Live engines: preserve dispatch-time project content, mutation/structural
  revisions and existing undo/redo checkpoints. Ordinary project opening may
  already have reconciled data; this slice does not undo that loading behavior.
- Preserve session semantics: fresh engine-instance identity, history defaults,
  counter reseeding, successful operation journals, workflow `run_id`, execution
  revisions and existing failed-attempt bookkeeping. Whole `history_summary`
  equality is not appropriate when its operation-log fields legitimately change.

### Constructor And Apply Seams

Prefer request-aware engine factories (tentative `from_state_for_operation` /
`from_state_for_workflow`) over changing ordinary `from_state` or introducing a
sticky "read-only engine" mode. A covered factory initializes normal session
identity/defaults and reseeds counters from existing state without modifying it.
It may later execute an ordinary mutation, which must retain its legacy repair.
Expose only the small Rust boundary needed by binaries; keep policy decisions
in the engine, with no dependency from protocol/engine back to frontend modules.

Consult the same facts at the start of `apply_internal`, before validation or
refusal can cause repair. Retain successful journaling and no-checkpoint behavior
through both apply methods and both workflow methods; do not rewrite workflows
as nested `apply` calls that change `run_id`, double-journal or add checkpoints.
Do not remove journals to make a persistence test pass.

### Workflow Aggregation And Failure Semantics

| Workflow | Recommended aggregate |
| --- | --- |
| Non-empty, every op covered | Preserve construction/project content, forbid implicit save; explicit report intents combine without becoming project mutations. |
| Any mutation or unmigrated op, in either order | Legacy aggregate construction/save fallback; each covered child still skips its own repair. This is not a whole-workflow no-change guarantee. |
| Empty workflow | Preserve legacy behavior for compatibility; do not gain an exemption through vacuous `all()`. Ask Claude to challenge this default. |
| Invalid/unresolved payload | No guessed inspection exemption; preserve existing structured refusal and explicit input-read errors. |

Keep operation order, child checkpoint boundaries, workflow journals, progress
and cancellation intact. Preserve each host's error behavior: in-process failed
workflows can retain an executed prefix; CLI/MCP raw file routes currently save
only after overall success; JS/Lua can raise without returning state; managed
failure discards unpublished engine deltas. Do not make these transactional or
add error-path autosaves. Already written explicit artifacts are not rolled back.
A cancelled Plan report returned successfully is not the same as workflow
cancellation/error, and a written file is not a completed scientific search.

### Resolve Once, Then Select Construction

Thin Op/Workflow adapters can select policy from their decoded requests before
constructing the engine. Named/shared-shell routes need an existing-parser-based
resolution boundary: retain typed `Operation`/`Workflow` plus captured payload
identity for execution, rather than read `@file` for classification and read it
again during dispatch. A small Rust-only resolved envelope/helper is acceptable;
a second command grammar, operation registry or scientific executor is not.
Use the captured request for lifecycle, applicable approval checks and execution.
Where a host already checks an original wrapper's conservative confirmation,
keep that gate additive: raw `op @file` must not become unconfirmed merely
because capture reveals a no-path Plan. This is a prevention requirement, not
a claim that the current file-backed route already bypasses approval.
Retain denials that currently precede payload I/O; resolution must not become
a reason to read files before a conservative wrapper/tool gate can reject.

Resolution I/O belongs in the existing execution/worker phase, not the pure
selector, GUI frame or nonblocking admission. Preserve original approved request
bytes/receipt hashes and bind the policy to the actual executed payload. Reuse
shared parser/error rules, not normalization that changes digest-bound inputs.

## Host Scope And Explicit Limits

1. Wire raw CLI op/workflow, forwarded typed commands and `shell op|workflow`
   through request-aware construction and the shared implicit-save fact. No new
   universal CLI confirmation flag; direct invocation remains user intent.
2. Wire raw MCP op/workflow and applicable shared-shell tools coherently. Keep
   raw `confirm=true` checks before payload parsing/resolution, project loading
   and construction, not merely before dispatch; retain all legacy exemptions
   for operations outside this slice. Do not infer successful artifact writes.
3. Cover generic JS/Lua direct operation/workflow state-out using the same
   factories, with real runtime tests under their optional features. Python
   inherits the CLI contract and needs actual-process verification, not a new
   Python policy implementation. Explicit `save_project` is not autosave.
4. GUI live engines use apply policy without reconstruction. Keep dirty-state
   checks revision-based; covered shell workflows can return false from typed
   policy without two full-project serializations. Mixed/unmigrated routes keep
   their existing change-detection fallback.
5. Managed Op/Workflow must keep required project snapshots, owner/history/
   structural checks, cancellation and short guarded publication. Recommended
   default: retain current journal/execution publication while asserting unchanged
   `ProjectState`, mutation revision and checkpoints. Do not route Plan/S1 into
   RNA's no-project `Observe` path, drop loaded targets, or claim no commit of
   any session bookkeeping. Ask Claude whether this bounded approach is
   compatible with the architecture's read-only worker rule or needs a smaller
   explicit journal-publication distinction within this slice.

Opaque `agents execute-plan`, macro/template expansion and other wrappers that
reconstruct before discovering their effective operation remain legacy unless
already resolved through the shared boundary above. This is a deliberate proposed
cut line, **not** a promise that every spelling/state-out wrapper is fixed.
Claude must assess whether this residual makes the slice misleading and must
be included now. Never expand templates or call models during policy selection;
never change nested-agent, mutation-approval or managed admission guards.

## Required Fail-First Tests

Reuse existing explicitly synthetic sequence-design/S1 fixtures and document
any new fixture consumer. Generate temporary projects locally; no private DNA,
reference downloads, model calls or native screenshots are needed.

| Boundary | Required positive and negative controls |
| --- | --- |
| Construction | Sequence-only state missing lineage/containers is genuinely repairable under ordinary `from_state`, but unchanged under covered factories. Verify fresh instance IDs, reseeded unique op IDs and normal defaults; later mutation still repairs. |
| Apply/history | Plan without/with a report and S1 preserve state, mutation/structural revisions, checkpoint contents and redo on success/refusal. Check invalid binding, invalid path, exhausted/infeasible/cancelled reports; retain existing journal/counter semantics. |
| Workflows | S1/Plan in both orders; multiple requested reports; a failing later step; cancellation before first, after a written prefix and before final result publication. Empty, Apply, SetParameter and SaveFile/unmigrated mixed controls retain legacy behavior and `run_id`. |
| CLI/MCP actual ingress | Raw op, named plan, shared-shell op/workflow and raw workflow; missing/empty/existing noncanonical project bytes; explicit report only. Keep raw MCP false/missing confirmation refused before payload/state I/O. Never use a helper-only pass as file-boundary proof. |
| Resolution | Inline and `@file` payloads execute the captured request, including a deterministic file-change-between-resolution-and-dispatch control. Keep original command/approval/receipt identity, conservative file-wrapper confirmation and existing invalid-JSON errors. |
| JS/Lua | Execute real generic op/workflow runtime bindings against unreconciled state; compare decoded state-in/state-out. Cover loaded and inline targets, refusals and mutation controls. Source grep or mocked clients are not VM execution. |
| Python | Actual supplied CLI for op/workflow/shell, missing and byte-preserved existing projects; mocked wrapper tests alone do not close parity. |
| Managed GUI/shared service | Valid loaded-target snapshot, preserved live state/history/dirty revisions, retained journal semantics, replaced-project and stale-history refusals, cancellation before/after compute and explicit artifact survival. No weakened lock/admission guard. |
| Load/fallback | Candidate sidecar hydration and strict/permissive load diagnostics stay unchanged. Existing RNA file-only exemptions, ordinary loads, Apply digests, SaveFile behavior and explicit save requests remain controls. |

`assert_sequence_design_mcp_parity` currently reads a preview-created project.
First assert preview does not create that file; explicitly initialize a baseline
project before the separate Apply phase, without manufacturing an operation to
seed counters. Keep confirmation/refusal, scientific report equality, exact DNA,
hashes and undo assertions intact. Audit tutorial 06.07 and its actual CLI
companion for similar reliance; do not weaken tests or automatically bless new
scientific output/approval hashes to make the lifecycle change pass.

## Proposed Implementation Order After Review And Owner Approval

1. Engine facts/factories, apply reconciliation and workflow aggregation, plus
   fail-first state/history controls. Keep the slices 1--2 approval tests green.
2. Shared resolution and thin CLI/MCP/JS/Lua/Python consumers together with
   direct file/state-out tests and deliberate preview/Apply fixture separation.
3. Managed/GUI controls and the reviewed wrapper cut line. Stop for an owner
   decision if preserving journals requires a broader scheduler/wire change.
4. Reconcile protocol/CLI/GUI/tutorial wording and the existing decisions only
   after implementation is tested. No generated/maturity churn unless a
   canonical check proves it is required; use the consolidated Glen handoff,
   not a separate branch freeze or repeated audit campaign.

Suggested new test prefix: `lifecycle_slice3_` (not present yet). Future runs
must prove the requested tests executed, not report success for zero matches:

```sh
cargo test --locked --offline --no-default-features -j 1 --lib lifecycle_slice3_ -- --test-threads=1
cargo test --locked --offline --no-default-features -j 1 --lib sequence_design -- --test-threads=1
cargo test --locked --offline --no-default-features -j 1 --lib genome_edit_intent -- --test-threads=1
cargo test --locked --offline --no-default-features --features js-interface -j 1 --lib lifecycle_slice3_ -- --test-threads=1
cargo test --locked --offline --no-default-features --features lua-interface -j 1 --lib lifecycle_slice3_ -- --test-threads=1
cargo test --locked --offline -j 1 --test capability_registry_parity --test parity_matrix_freshness
cargo check -q --locked --offline -j 1
cargo fmt --check
GENTLE_TUTORIAL_BIN_DIR=/absolute/path/to/same-source/binaries python3 -m unittest scripts.test_sequence_design_tutorial -q
```

Also require a new actual CLI/MCP process boundary test and actual Python client
coverage, with their final test commands documented at implementation time.
Missing optional JS/V8/Lua dependencies are blocked/unrun evidence, not a skip
count converted into parity. No dependency downloads/builds are authorized by
this planning document. Preserve existing LF/CRLF tests; add scoped checkout
coverage/provenance if any new byte-bound fixture or freshness gate is introduced.
Use shared shell quoting and native Windows handle/path tests, without hash
weakening, production timeout changes or fixture identity normalization.

## Review Questions And Exit Boundary

- Is the explicit-factory plus typed policy boundary sufficient and minimal?
  Should shell resolution return an engine/run pair or a small resolved request?
- Is legacy fallback for empty/mixed workflows the right compatibility choice?
  Verify failure-prefix behavior per host rather than imposing one global rule.
- Must opaque agent/template wrappers join this slice, or can the stated cut
  line remain without claiming universal readonly behavior or unblocking S2?
- Can managed journal publication stay on the existing guarded path without
  weakening the read-only-worker invariant or silently discarding receipts?
- Are the decoded-state/post-loader baselines and sidecar limitations honest?
  Is any additional persisted mutation performed after the identified preamble?

Explicit report/project-path aliases, symlink/hard-link/race safety, writer
overwrite/atomicity/no-clobber changes remain Slice 4. **Forbid autosave is not
a guarantee that a user-selected report path cannot overwrite the project.**
Use disjoint report/project paths in preservation tests; do not claim alias
refusal exists. No S2, new editor/protocol biology, broad registry reclassification,
global constructor/loader change, asynchronous-service redesign or new release
gate is authorized. Review feedback and owner implementation approval are still
required; this plan itself grants neither implementation nor acceptance.
