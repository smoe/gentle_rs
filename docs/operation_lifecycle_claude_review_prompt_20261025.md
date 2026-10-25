# Claude Review Prompt: PlanDnaSequenceDesign Lifecycle

Prepared 2026-10-25 for Steffen to forward. No direct Claude invocation or
review has occurred. This is a design-only request before implementation.

---

Please give a read-only software-contract review of GENtle's proposed operation
lifecycle policy, starting with `PlanDnaSequenceDesign`. Your previous S1 review
accepted the rebuttal of deriving all exemptions from `CapabilityMutation` and
recommended per-dimension policy before editing S2. This proposal addresses
that prerequisite; it does not implement S2 or reopen completed S1 arithmetic.

## Repository And Safety

Worktree: `/Users/u005069/.codex/worktrees/47dd/gentle_rs`.
Branch: `gentle_rs_2_main`.
Preparation HEAD and local main:
`81b675048ab0c886185d55e101551790d8c3bea5`.
No fetch was performed. Read the **working-copy** proposal, not only HEAD:
`docs/operation_lifecycle_policy_proposal.md` and this prompt. They may be
uncommitted when forwarded. State any source drift or identity limit you can
verify with read-only tools; do not silently substitute the primary checkout.

Read `AGENTS.md`, relevant architecture/decision/roadmap sections, and existing
sequence-design/S1 contracts before judging the proposal. Do not edit files,
run shell commands/builds/tests, fetch/rebase/switch branches, invoke another
agent, inspect credentials, access private sequences or unrelated `paper/`
and `output/`, upload data, commit/push, dispatch CI or change release status.
Use Read/Grep/Glob only. No literature/editor-chemistry review is requested.

## Original Codex Proposal

1. Add one typed, request-aware engine policy selector for the covered preview/
   apply operations, with separate project effects, explicit artifacts,
   reconciliation, undo/redo, implicit project autosave and conversation approval.
   Keep successful operation journals/execution bookkeeping distinct from
   project mutation. Never derive permission from one registry label.
2. Plan without a path: preserve project content/history, skip persistent
   repair, emit no file, and prevent implicit CLI/MCP project saves. Plan with
   a path: the same project policy, but explicitly write the requested report
   and require artifact approval in delegated execution. Empty/bad paths never
   become a no-path exemption. Preserve explicit direct CLI intent and raw
   MCP's independent `confirm=true` gate.
3. Preserve S1's current no-reconcile/no-save semantics and separately approved
   apply's mutation/checkpoint/digest behavior. Unknown/unmigrated operations
   retain legacy conservative behavior; `SaveFile` gains no new exemption.
4. Wire constructor, apply/checkpoint, shell confirmation/state-change and raw
   CLI/MCP save seams coherently. Handle covered workflows and verify managed/
   scripting consumers without replacing their stale-owner or tool gates.
5. Recommend projecting Plan as project-read-only with conditional external
   artifact details. Static family metadata must conservatively disclose the
   possible approval requirement; a concrete no-path request can have no
   operation-derived requirement, but explicit candidate/tool gates stay binding.
   Challenge this metadata projection rather than assuming it is settled.
6. Reject explicit preview output aliasing the selected project file before
   writing, using existing portable helpers only if sound. Do not conflate
   autosave suppression with destination safety or silently change every writer's
   overwrite/atomicity contract. Flag a narrower alternative if this is too broad.
7. Characterize scientific/report/journal bytes, add fail-first boundary tests,
   migrate the covered lifecycle coherently, then align metadata/generated
   projections/docs. Preserve current cancelled-report semantics, explicit
   approval and unrelated operation behavior. No runtime change before review.

Full proposal and test matrix: `docs/operation_lifecycle_policy_proposal.md`.
All proposed behavior changes must be distinguished from current facts.

## Evidence Already Gathered

- `crates/gentle-protocol/src/lib.rs::infer_engine_operation_mutation`:
  `SaveFile` is `ReadOnly`; Plan falls through to `Mutating`. The detailed
  shell descriptor contradicts that coarse classification.
- `src/engine_shell/command_parsers.rs::parse_sequence_design_command` lowers
  named plan/apply commands to the same typed generic `Op`; do not propose a
  second runner for named commands.
- `src/engine_shell.rs::sequence_design_capability_descriptor` calls planning
  nonmutating/unconfirmed yet lists an optional external artifact effect.
  `ShellCommand::is_state_mutating` exempts Plan regardless of output path.
- `src/agent_execution.rs::execute_agent_plan_candidate` (verify the current
  symbol/call chain) checks parsed command policy as well as candidate flags.
  A request-sensitive correction must reach stored/forged-flag proposals, not
  just change tutorial or discovery wording.
- `src/engine.rs::{from_state,operation_checkpoint_kind,apply,apply_with_progress}`
  and workflow counterparts; `src/engine/ops/operation_handlers.rs::apply_internal`:
  construction and preamble can reconcile, while planning has no undo checkpoint.
  Success still appends operation records. S1 alone skips the preamble repair.
- Plan's operation-handler arm writes via `write_pretty_json_file`, which
  creates parent directories and uses `fs::write`; it is not no-clobber/atomic.
  The pure `src/engine/analysis/sequence_design.rs::plan_dna_sequence_design`
  takes `&self`. Cancellation can return a report, not just an error.
- `src/bin/gentle_cli.rs::{main/run dispatcher,execute_direct_operation}` and
  `src/mcp_server.rs::{op_tool_result,workflow_tool_result}`: raw operations
  currently save project files after success except for S1; shared-shell paths
  instead use `state_changed`. Verify exact symbols rather than trusting names.
- `src/js_interface.rs::apply_operation_impl`/`apply_workflow` and the generic
  `apply_operation`/`apply_workflow` functions in `src/lua_interface.rs` construct
  via `from_state` and return state. Python's `integrations/python/gentle_py/client.py`
  delegates op/workflow to CLI. Inspect these before claiming adapter parity;
  keep the policy headless and avoid a protocol-to-root dependency cycle.
- `src/engine/analysis/sequence_design/tests.rs::assert_sequence_design_shell_parity`
  expects the with-path command to be nonmutating; preserve scientific assertions
  while reviewing what this boolean authorizes.
- `src/mcp_server.rs::assert_sequence_design_mcp_parity` reads a state file
  created by planning before applying. A no-autosave migration must change that
  fixture deliberately, not permit apply without its normal approval.
- `src/engine/analysis/genome_editing/tests.rs`, CLI and MCP S1 file tests:
  `ced03a22` added explicit unreconciled-state controls. That follow-up has not
  had an executed Rust verdict from Codex/Claude; do not claim otherwise.

Relevant documents: `docs/decisions.md` (S1, synthetic design, DEC-013, DEC-026,
DEC-048), `docs/architecture.md`, `docs/dna_sequence_optimization_plan.md`,
`docs/base_prime_editing_preparation.md`, `docs/cli.md`, `docs/protocol.md`,
`docs/gui.md`, `docs/testing.md` portability rules, `docs/roadmap.md`, and
`docs/glen_next_verification_handoff.md`. No active tutorial, metadata, generated
artifact or durable decision has been changed to pretend the policy exists.

## Requested Review Output

Please return a concise, source-grounded critique and revised minimal plan:

1. Verify/correct the current-state evidence and identify missed constructor,
   workflow, macro, host, direct-method or state-save seams. Distinguish project
   files, requested artifacts and transient journals; show concrete paths.
2. Judge the proposed target matrix. Specifically test the reasoning for both
   Plan path forms, S1, Apply and unchanged `SaveFile`; do not suggest one
   boolean or prefix heuristic that weakens existing gates.
3. Resolve the metadata question: project-`ReadOnly` plus explicit artifacts
   versus conservative `External`. Explain static versus request-specific
   approval, candidate-required confirmation and raw MCP tool gates.
4. Review output/project path collisions, cancellation and any partial file
   effects. State what existing helpers can actually guarantee and which
   claims require native platform tests or separate authorization.
5. List only concrete defects/risks, ranked by severity with file/symbol evidence;
   say explicitly where no new defect exists. Propose the smallest coherent
   implementation order and independent negative controls, including LF/CRLF
   and native Windows behavior. Flag unnecessary scope or unsafe fallbacks.

Source inspection is not executed Rust, filesystem, native GUI, scientific or
release acceptance. End with the exact scope/identity and what you did not run.
