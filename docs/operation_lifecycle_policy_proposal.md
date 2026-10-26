# Request-Aware Operation Lifecycle Policy

Status, 2026-10-27: slices 1--2 retain the source-scoped executed macOS approval/
parity evidence below. After a further owner-supplied Claude review, the owner
authorized [revised Slice 3](operation_lifecycle_slice3_plan.md), explicitly
without the new payload-resolution layer. Construction/workflow/no-autosave
changes and regressions are authored; Rust/actual-process/VM/native verification
remains with Glen/CI. The original review/source observations below are historical;
the revised slice order records current scope. Alias/writer rules and S2 remain
separately gated. Original review prompt:
[read-only Claude request](operation_lifecycle_claude_review_prompt_20261025.md).

## Source And Scope

Repository: `/Users/u005069/.codex/worktrees/47dd/gentle_rs`.
Proposal preparation baseline on `gentle_rs_2_main` (then local main/HEAD):
`81b675048ab0c886185d55e101551790d8c3bea5`. The clean tracked branch was rebased
onto local main before this inspection; no fetch was performed. These documents
were prepared as working-copy drafts; committing them does not adopt the policy.
Untracked `output/` is unrelated and untouched.

The review was forwarded from Claude's inspection of `adacb2bc`; Codex rechecked
the load-bearing paths at local main `b544a21b81b0d98896599f4ae21cdf60eea59b62`
before editing. The primary checkout's untracked
`docs/operation_lifecycle_codex_prompt_20261025.md` was read, not copied or
modified. No fetch, Rust execution, model call or private-data use was performed
during that source review. Later authorized execution is recorded below.

## Review Reconciliation And Approved Scope

Original Codex plan: migrate optional-path planning through a full lifecycle
selector and construction/persistence coverage before S2; it incorrectly framed
the missing approval as an invariant to preserve and preferred static ReadOnly.
Claude verified the source claims and identified these corrections:

| Finding | Revised Codex action |
| --- | --- |
| F1: ordinary planner report writes need no confirmation | Accept as a present permission defect, not just forged metadata. Implement the engine-owned request approval facet and reuse the shared parsed-command gate. Add ordinary planner and forged/wrapped regressions. Source repair is not an executed acceptance claim. |
| F2: static external precedent already exists | Use `External` in the canonical registry and `external`/confirmation in detailed discovery, like TSS SVG export. Only a concrete absent-path request adds no operation-derived approval; explicit candidate/tool requirements stay additive. |
| F3: S1 decision overstates scripting preservation | Narrow the decision now: apply-time no-reconcile and named CLI/MCP file-save suppression do not preserve unreconciled JS/Lua constructor state-out. Leave constructor migration separate. |
| F4: approval and destination restrictions are independent | Defer alias/collision and writer changes explicitly; they must not delay confirmation. Do not claim project/report alias safety or implicit-save suppression for Plan. |
| F5: one existing assertion must flip | Change `assert_sequence_design_shell_parity`'s with-path confirmation assertion only. Preserve saved/direct/report equality checks and exact apply/undo assertions. |
| F6: approval covers cancelled-report writes too | Approve the write attempt, not search success. Retain a saved report's `cancelled` status. A completed file write is publication evidence, not a completed search; do not conflate the two. |

Two qualifications to the review: planner candidates derive their mutation flag
from the parsed command; static descriptor changes alone cannot protect forged
flags. Also, actually writing a cancelled report can succeed: the scientific
outcome must remain cancelled, not be relabelled completed because publication
succeeded. The existing writer/report bytes are unchanged.

Static alignment also exposed a missing source fact: neither existing raw
sequence-design operation ID was in `PUBLIC_ENGINE_OPERATION_NAMES`, despite
the named glossary commands. Register both existing variants, with Plan External
and Apply's unchanged Mutating classification. This is required discovery/parity
alignment, not a new operation or adjacent biological feature. Add their two
shell-only matrix rows and counts using the existing projection rules; canonical
Rust byte-freshness is now checked at the recorded source below. Existing
LF/CRLF matrix coverage is retained, with no new freshness gate or maturity
regeneration.

The implemented helper is `engine/operation_lifecycle.rs::operation_approval_policy`.
It selects only `PlanDnaSequenceDesign`, inspecting `path.is_some()` through
typed variants without I/O or registry inference. Every present path requires
confirmation, even empty/malformed strings. `None` for any other operation means
legacy behavior, not an exemption. This is only the approval facet, not the full
policy type or the construction/persistence guarantees proposed below.

Start with `PlanDnaSequenceDesign`, distinguish its absent/present output path,
and use `EvaluateGenomeEditIntent`, `ApplyDnaSequenceDesign` and `SaveFile` as
controls. Do not migrate every operation or infer policy from names, prefixes,
the registry's coarse `CapabilityMutation`, or model-supplied declarations.
The [S1 review reconciliation](base_prime_editing_preparation.md#landed-s1-review-reconciliation-2026-10-25)
and [existing safety decision](decisions.md#read-only-genome-edit-intent) explain
why a single read-only boolean is insufficient.

## Pre-Repair Behavior, Not Desired Guarantees

These observations were rechecked before the bounded repair at `b544a21b`.
They are source inspection, not an executed regression or security exploit.

| Boundary | Pre-repair source fact |
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

The existing shell parity test expected a plan with `--path` to have
`is_state_mutating()==false`. Slice 1 flips that assertion; its scientific
report equalities remain untouched. The compatibility boolean includes external
write approval, not only persistent project mutation.
The S1 unreconciled-state regression in `ced03a22` remains a control; its later
headless macOS execution is recorded below, not inherited from source review.

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

### Full Lifecycle Target Beyond Slices 1--2

Only Plan's final approval column is implemented by the current bounded repair.
The construction, reconciliation, history and autosave columns remain future
targets, not guarantees of this partial migration.

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

## Full Migration Consumer Rules

The approval/static-metadata parts of rules 2 and 5 are implemented in source.
Other lifecycle changes remain deferred pending separate owner authorization.

1. **Engine construction and apply:** add an explicitly selected inspection
   construction path that retains project content without persistent repairs;
   leave ordinary `from_state` unchanged. Preserve the necessary session ID,
   history defaults and op-counter reseeding. Both engine apply paths and
   workflows must consult the same reconciliation/history facts. The pure plan
   reads sequence/annotation data, not repaired container/lineage indexes.
2. **Shell and agents:** use the request policy both for confirmation and for
   `state_changed`. Keep `is_state_mutating` as a compatibility wrapper if needed,
   but document that the existing boolean also enforces external-write approval.
   Repair the existing normal-planner defect as well as forged flags: neither
   may run a plan with an output path without confirmation. Named command and
   generic `op` must agree. Approval is for the write attempt, including a
   cancelled report, not an assurance of successful search or publication.
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
5. **Metadata:** project Plan's static family as `External` with confirmation,
   retaining explicit optional-artifact details; the coarse label is never an
   authorization source. A concrete parsed no-path
   request can have no operation-derived approval. Explicit candidate/tool
   requirements are never removed. Reconcile canonical registry, detailed
   introspection, effects and generated projections together. This replaces the
   original ReadOnly recommendation following Claude's source-grounded review.

### Output And Cancellation Boundaries Requiring Review

Forbid *implicit* project autosave does not prevent a caller from supplying
the project filename as the explicit report path. Proposed rule: reject a
preview destination aliasing the selected project file before writing, including
existing aliases supported by the repository's cross-platform path helpers.
Inspect those helpers before choosing the implementation. Do not promise
symlink/hard-link or race safety without an executed design/test path. This is
explicitly deferred from slices 1--2 and must not delay approval repair. The
current approval facet makes no destination or strict project-file-preservation
guarantee for Plan.

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

## Revised Slice Order

1. **Confirmation repair plus decision-record correction.** Owner-authorized;
   source implemented and covered by the focused macOS execution below. The
   engine-owned approval facet covers only Plan's optional path and the shared
   parsed-command consumers.
   Ordinary planner, forged Shell/Op/Workflow and macro-wrapped cases have
   executed headless regressions, including approved cancelled-report publication
   and additive no-path candidate gates. Narrow S1's scripting claim now rather
   than imply its constructor is already covered. Flip only the historical
   with-path boolean assertion; preserve scientific/report/apply equalities.
2. **Conservative static metadata.** Owner-authorized; source implemented with
   canonical Rust freshness/parity checked at the recorded source. Plan's
   operation and named-command families are External; detailed discovery requires
   confirmation and explains concrete absent-path handling. Do not change `SaveFile`, S1 or
   Apply classification, glossary routes or maturity attestations. Existing
   matrix surfacing columns do not encode mutation, but the missing raw IDs
   require two shell-only rows/counts. No new byte-exact artifact gate or
   unrelated regeneration is introduced. Retain existing LF/CRLF checks.
3. **Construction, scripting and workflow lifecycle coverage.** Owner authorized
   the Claude-reviewed revision on 2026-10-26 **without a new payload-resolution
   layer**. The [reconciled Slice 3](operation_lifecycle_slice3_plan.md) is authored;
   Rust/actual-process/VM/native verification remains with Glen/CI. One explicit
   decoded-request factory/selector covers Plan/S1 and nonempty all-covered
   workflows; ordinary construction, mixed/empty/opaque wrappers and parsers
   stay legacy. Preserve journals and guarded managed commits, not just DNA.
   The MCP fixture now asserts preview/denied Apply leave state absent; approved
   Apply creates it. No artificial baseline or new request/digest is needed.
4. **Alias/collision and writer contract.** Separate owner decision; not a
   prerequisite for slices 1--3. Inspect existing portable helpers before any
   promises about symlinks, hard links, races, atomicity or no-clobber behavior.

Keep S1 regressions as a control, not a reason to repeat its scientific review.
Do not expand this to every read-only operation before evidence warrants it.
S2 still needs separate implementation authorization after this prerequisite.

## Full Target Verification

Approval/metadata and unchanged scientific-report comparisons apply to the
current repair. Constructor, autosave, strict project-content preservation and
alias/writer rows below remain requirements for the separately authorized later
slices, not claims that slices 1--2 implement them.

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
The owner-requested local checks below are a bounded exception to the earlier
no-local-Rust restriction, not authorization to implement the remaining policy.
Glen/CI verification stays in the existing consolidated handoff. The full target,
native, scientific and performance acceptance are not established by these checks.

## Remaining Owner Decisions

Claude's source review is received; External versus ReadOnly and write-attempt
approval are no longer unanswered design questions for slices 1--2. The executed
evidence below supports the bounded confirmation repair on the named macOS
source/features. The `ced03a22` S1 control also ran headlessly; neither result
transfers to unrun constructors, platforms or later integrated source.

Separate authorization is still required for construction/scripting/autosave
migration, the destination-alias/writer policy, and S2 itself. Strict project
preservation cannot be claimed from a confirmation guard or unchanged disk
bytes alone. There is no new release gate, writer redesign or S2 editor work.

## Authoring Checks Before Rust Execution

macOS, 2026-10-26, based on `b544a21b` with the scoped working-copy changes:

```sh
python3 -m unittest scripts.test_sequence_design_tutorial scripts.test_primer_trace_contract scripts.test_tutorial_checkouts scripts.test_capability_maturity -q
python3 scripts/check_documentation_status.py
python3 scripts/generate_capability_maturity.py --check
rustfmt --edition 2024 --check --config skip_children=true src/engine/operation_lifecycle.rs src/engine.rs src/engine_shell.rs src/engine_shell/tests.rs src/agent_planner.rs src/agent_execution.rs src/agent_bridge.rs src/engine/analysis/sequence_design/tests.rs crates/gentle-protocol/src/lib.rs tests/capability_registry_parity.rs
git diff --check
```

The combined Python run discovers 103 tests: 99 pass, four explicitly configured
binary-dependent CLI replays skip. A later focused source/maturity run discovers
28 tests: 24 pass, the same four skip. These are source/teaching/checkout checks,
not execution of the Rust permission gate. LF/CRLF matrix/fixture coverage and
matrix count consistency pass. No fresh Rust canonical matrix generation or
byte-freshness result is claimed. Maturity projections and `Cargo.lock` remain
unchanged. Documentation status, formatting and whitespace checks pass.

At this authoring stage, six new Rust regressions, Cargo check and the S1 control
were still unrun; the subsequent execution record follows. No new scientific
fixture was committed; the existing teaching input's README records its
approval-test consumers.
The source diff leaves solver, writer, operation handler, CLI/MCP persistence
and JS/Lua constructors untouched. The added raw registry IDs and two-row matrix
projection are the necessary static-alignment extension, not an adjacent feature.
Session-close reports four OK, two expected dirty-tree/manual-plan warnings and
zero failures; unrelated untracked `output/` is untouched. Committing these
changes is not execution evidence. No push, workflow dispatch, model call,
native or scientific acceptance is performed here.

## Executed Local Verification (2026-10-26)

The owner requested parity execution and rebasing, superseding the earlier local
Rust restriction for this bounded verification. No runtime repair or generated
artifact change was needed. Sources and features must remain distinguished:

| Tested source | Features | Executed results |
| --- | --- | --- |
| `ab80743e300a42083a24de2c63dcabf680b5335b`, on local main `07389910` | Defaults: `desktop-gui`, `screenshot-capture` | Registry/adapter integration 21/21; matrix freshness 1/1; canonical matrix check OK, catalog 69 entries OK, manifest 30 chapters OK; configured tutorial companion 9/9, including its four actual-CLI cases. |
| `98c9331aea3bb4062447d91db0052ceb832580a5`, rebased onto local main `f1a60ae4` | Defaults | `cargo check` and whole-workspace formatting pass. |
| `98c9331aea3bb4062447d91db0052ceb832580a5` | `--no-default-features` | Sequence-design library filter 32/32, including all five new shared-engine approval regressions and four existing MCP preview/apply parity tests; S1 filter 14/14, including unreconciled-state and MCP state-file controls. No ignored tests in these filters. |

The sixth new regression is the registry-family test in the first row. Rebase
added RNA payload repairs, not changes to the approval/registry/parity files;
this does not turn the earlier desktop binary into a newly built-source pass.
The default-feature library-test build was interrupted for the requested rebase
before running tests (exit 130), not counted as a failure or pass.

Exact commands, from the checkout root, with no profile/RUSTFLAGS overrides:

```sh
cargo test --locked --offline -j 1 --test capability_registry_parity --test parity_matrix_freshness
target/debug/gentle_examples_docs parity-matrix-check
target/debug/gentle_examples_docs tutorial-catalog-check
target/debug/gentle_examples_docs tutorial-manifest-check
GENTLE_TUTORIAL_BIN_DIR="$PWD/target/debug" python3 -m unittest scripts.test_sequence_design_tutorial -q
cargo check -q --locked --offline -j 1
cargo fmt --check
cargo test --locked --offline --no-default-features -j 1 --lib sequence_design -- --test-threads=1
cargo test --locked --offline --no-default-features -j 1 --lib genome_edit_intent -- --test-threads=1
```

Host: Darwin 27.0.0 arm64 / `aarch64-apple-darwin`; rustc
`1.100.0-beta.5 (8ea38001b 2026-10-13)`; Cargo
`1.100.0-beta.5 (7552a885e 2026-09-29)`. Cargo test/dev uses opt-level 0,
debug 1, incremental true, 256 codegen units, dependency debug 0, and default
debug/overflow checks. The checkout config supplies a 16-MiB Rust thread stack.

SHA-256 identities:

```text
Cargo.lock (unchanged): 97671e924a77d61de6dad2e6c6fe3afbee440b8606d79ca0ad2b494fd32d66ad
ab80743e gentle_cli: 6b90b57bba45c9ffdef03200586ca84f5d8622f8c401f51f73409bb28c42fe32
ab80743e gentle_examples_docs: a3e6e3ea0680994b88f7679a3cdf6bb80676cd0fa33a1972b21b41301bdd199f
98c9331a headless library-test binary: c13309a2386c616cda28d0d3e896a88f42b82fe5fc5a765c72c4e3c22267204d
```

The fast Python run at `ab80743e` passes 99/103, with four optional CLI cases
skipped; the configured nine-test run above executes those four without skips.
The post-rebase focused source/maturity run passes 24/28 with the same four
optional skips. LF/CRLF coverage, documentation/maturity checks and whitespace
pass; no lockfile or maturity output regeneration. The linker reports a large
`__eh_frame` warning; headless compilation also reports three unused existing
GUI diagnostics helpers. These warnings were not suppressed or repaired.

Only the two owner-approved inactive incremental caches were removed after an
open-file check; the active build, sources and unrelated `output/` were preserved.
Session-close: four OK, two warnings (intentional documentation edits/existing
`output/`, and the manual plan-fidelity reminder), no failures. The final diff
only records verification and reconciles status; no additional policy is adopted.
Full-suite, native GUI, Windows/Linux execution, live-model, package, performance
and scientific acceptance remain unrun here. No push, fetch or workflow dispatch.
Construction/reconciliation/autosave, destination safety and S2 remain separate.
