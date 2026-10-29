# Lifecycle Slice 5: Snapshot-Bound Preview Publication

Status, 2026-10-29: **implemented in `b0521123`, with follow-up tests in
`98f69995`; the owner-forwarded Claude report records 28 passing Rust tests at
`284e5b70`**. See the [scoped reported result](#claude-reported-rust-verification-at-284e5b70),
not a new-source or native acceptance claim. The later
owner-authorized [284e5b70 rehearsal](glen_next_verification_handoff.md#synthetic-lifecycle-rehearsal-at-284e5b70-2026-10-28)
passes 25 source/process checks and scoped macOS GUI observations. Native
cancellation/concurrency is inconclusive; remaining platform/VM/full-suite
acceptance stays separate. No Rust build, model call or push in that rehearsal;
the prior goal stays paused. Authoring-time unrun lists below remain historical,
not a request to repeat the explicitly executed subset.

The [original review prompt](operation_lifecycle_slice5_claude_review_prompt.md)
records Codex's proposal. The owner supplied Claude's review against local main
`550e6e0ad67fbb301e9a10bb4cc996924db3388f`; Codex checked the implicated source
at branch HEAD `550725695437e602c001bd9c881a2f66f861e33c`. Those runtime paths
have no HEAD-to-main diff. Main's intervening observed-transcript S2 planning
is unrelated. These are source observations, not a new acceptance freeze.

## Original Plan And Review

Original Codex scope: only standalone no-path `PlanDnaSequenceDesign` through
the managed service, using a project snapshot but publishing a per-command
result instead of merging its private operation log. Other operations, report
writers, workflows and raw adapter semantics remain unchanged.

Claude recommends revise, then implement. Accept its three essential safeguards:
a preview still needs the live project's snapshot; publication must check owner
and structural freshness; the GUI must not infer project ownership solely from
whether a result commits. Also accept a separate publication branch, explicit
cancel ordering, existing bounded receipts instead of a new audit store, and
unchanged conservative S1/workflow controls. No generic journal merge.

Codex qualifications at authoring, before the later shipped-code review:

- `dna_design_source` already rejects an occupied output ID during preview as
  well as Apply. Keep that validation; a later collision is caught by structural
  freshness or fresh Apply validation, not by removing the preview check.
- Historical-result labelling is required, not optional. Carry the source
  binding into retained GUI output and the explicit sharing path; a one-time
  summary string at polling cannot detect a later edit. Existing
  `AgentCommandOutput` has no such binding. Keep any addition host-local and
  limited to this preview, not a scientific report field or payload framework.
- `CommandReceipt` is bounded, local and nonpersistent. Consumed terminal
  receipts become eligible for capacity eviction; consumption itself does not
  immediately delete them. Do not promise a durable replacement for History.
- A held engine write lock alone does not establish the computed-before-
  publication cancellation window. Tests must also synchronize on actual
  completion of work, not timing. Never sleep and assume publication was reached.
- Cargo accepts one positional filter before `--`; multiple filters after it
  are test-binary arguments. Separate Cargo-filter invocations remain valid.
  Require nonzero counts and record which filters actually executed.

The owner-forwarded shipped-code review accepts request/report binding and the
stronger delayed display/sharing safeguards. Do not retain a blanket claim that
those additions remain Claude-unreviewed, or infer an audit of every separate
qualification from that source review.

## Owner Decisions Accepted, 2026-10-28

1. Add a narrow DEC-026 exception: managed standalone no-path previews retain
   command results/receipts without publishing their private journal to the live
   project. Do not relax `commit_detached_execution` or its generic stale guard.
2. Accept the user-visible consequence: these managed previews stop adding
   project History/operation-log rows. Results and bounded local receipts remain;
   they are not persisted project history. Raw engine/CLI/MCP/JS/Lua behavior is
   not changed, and exact Apply still requires separate approval.

The owner explicitly accepted both decisions and requested implementation plus
tutorial 06.07. Claude's recommendation alone did not authorize these changes.
This slice adds no editor biology or release gate.

## Revised Minimal Implementation

### Selection And Snapshot

Select only typed `PlanDnaSequenceDesign { path: None, .. }` in the existing
worker classification/freezing path. A small `Policy::Preview`-like variant
uses the Execution lane and explicitly requires a project snapshot, just like
Project. Do not inherit the existing default `needs_project=false` of other
non-project variants. Keep admission's owner/structural/journal baseline checks.

`submit` and production `submit_parsed` must converge. Inline admission can
recognize the variant; `op @file` remains conservatively Project until existing
worker-only freezing, with identical baseline/lane requirements. Named commands
already lower to Op. No new resolver, early file read, repeated resolution,
operation schema, grammar, generic ReadOnly inference or capacity exemption.

### Publication, Freshness And Cancellation

Use a dedicated preview publication branch, not RNA's independent branch:

1. Require a typed sequence-design result for the selected request. Retain the
   bounded full output and its hash; no compact RNA fallback or silent truncation.
2. Take the live engine read lock, then the job lock, matching existing lock
   order. Before publishing, reject cancellation, replaced owner, changed
   structural revision or an unexpected `run.state_changed`.
3. Do not compare journal length at this publication boundary: this branch does
   not merge project history. Ordinary display/auxiliary metadata edits and Save
   As need not invalidate the preview. Accept conservative refusal on unrelated
   structural edits even for inline DNA. Do not change the admission guard.
4. Under the same locks mark publication as `preview`, retain result owner and
   snapshot structural revision plus result hash, set zero committed operations
   and `requires_project_commit=false`. A successful private apply still records
   its computed operation count. Drop the detached engine outside the locks.

Before-snapshot cancellation stays unchanged. Cancellation during search or
after computation but before publication refuses service publication, giving a
Cancelled receipt; no-file previews have no retained artifact to protect.
Cancelled search DNA must never be relabelled feasible/completed. Conversely,
`SearchExhausted`/`ProvenInfeasible` can be successfully published reports under
a Completed service receipt, without claiming feasibility or an applicable
approval. After publication, `cancel()` returns false; do not relabel the result.

Keep the existing result hashes and scientific report schema unchanged. A new
`publication="preview"` receipt value is a host-local additive distinction, not
authenticated search provenance or application approval.

### GUI And Inner Agent

At consumption, apply owner checks to results carrying a result owner as well
as results requiring commit. A result from a replaced project is not attributed
to its replacement or another turn; retain its receipt. Keep existing session/
turn handling. The project-change cancellation rule still covers this job while
it is pending publication; do not accidentally treat it as independent RNA work.

A same-owner result published before a later structural edit may be inspected
as **historical snapshot evidence**. Preserve owner/revision context in the
retained local output and reevaluate freshness on display and explicit sharing,
not just when polling. Do not modify the scientific JSON to insert a warning.
Use a small local binding or existing receipt association; no general context
registry, project serialization or source rehashing in the UI frame. A historical
share must explicitly identify its snapshot status. A different-owner share is
refused or cleared. Preserve bounded, user-initiated sharing and fresh Apply
source/feature/code/output-ID validation; no automatic send or Apply.

### Noncoverage Controls

Leave `path=Some(...)`, including invalid/empty paths, S1, all workflows,
macros/templates, opaque agent plans, arbitrary trusted `submit_work`, RNA lanes
and other operations unchanged. No generic state/journal merge, loader repair,
queue, broad observation migration or output-writer change. Slices 3/4 and raw
adapter contracts retain their own pending verification.

## Authored Test Scope

- Managed named/raw inline and loaded-target no-path previews: publication,
  owner/revision, one computed/zero committed operations, exact direct report/
  approval equality, unchanged live state/journal/revisions/undo/redo.
- Held genuine mutation survives a completed preview; two synchronized previews
  both publish. Keep S1's existing concurrent-journal stale rejection unchanged.
- Structural edit or replacement before publication refuses; display/metadata
  edits and Save As allow publication without restoring snapshot state/bindings.
- Cancel before snapshot, during search, after computed work but before
  publication, and after publication. Use barriers/channels and actual receipt
  observations to establish the intended stage. Check scientific and service
  statuses independently, including exhaustion/infeasibility.
- Existing worker `@file` path, admission/freeze agreement and Execution/project
  snapshot requirements. Retain with-path, workflow, S1 and arbitrary-work
  controls; split only the old Slice 3 test's standalone Op expectation.
- Publication followed by edit/replacement before polling, and edit/replacement
  after polling but before display/sharing. Historical status must reach the
  explicit agent draft, while wrong-owner results are not reattributed/shared.
- Preserve existing bounded-output errors. If realistic admitted previews cannot
  reach 16 MiB, record the end-to-end coverage limit; test the existing boundary
  with synthetic data where feasible, without a new production budget override.

Update DEC-026, architecture, lifecycle status, protocol/GUI docs where changed,
tutorial 06.07 and changelog. Adjust roadmap only for real status changes, with
no Done entry. Teach ephemeral command audit versus project History and exact
Apply approval; do not claim live-agent/native acceptance from authored tests.

## Implementation Ledger

`command_execution/rna_artifacts.rs` adds the closed Preview selection with
project snapshot/Execution admission, existing worker freezing and exact typed
request/report binding. `command_execution.rs` adds a separate preview branch
with owner/structural/cancel checks and no live journal merge. Generic commit,
admission and raw operation contracts are unchanged. Output limits have no
compact preview fallback; an infrastructure-only synthetic oversize test does
not claim a reachable 16 MiB solver case.

`app.rs` retains a tiny host-local owner/revision stamp in `AgentCommandOutput`.
The Agent Assistant consumer checks owner independently of commit, retains the
stamp after receipt consumption, labels historical output at display/sharing,
and refuses wrong-owner/busy-project sharing. Other output constructors retain
no stamp. Scientific JSON is unchanged and user sharing remains explicit.

Synthetic tests live in `command_execution/preview_tests.rs`, the existing
policy module and Agent Assistant command tests. Split only the Slice 3 Op
expectation; keep its S1/workflow controls and Slice 4 tests. Tutorial 06.07,
canonical GUI chapters, protocol, lifecycle docs, architecture/DEC-026,
roadmap/changelog and consolidated handoff describe authored behavior and
pending external verification, not new-source runtime acceptance.

Local non-Rust verification: 25 focused lifecycle/tutorial/maturity checks pass
with 12 explicit actual-process skips (37 total), plus all 83 LF/CRLF-aware
checkout tests. Manual projection freshness, documentation status, maturity
freshness and whitespace pass. Scoped rustfmt passes for the changed service/
GUI/test/agent modules. `app.rs` retains three existing Help formatting groups
from HEAD; stdin formatter comparisons confirm Slice 5 adds no new formatting
differences. Preserve that unrelated code rather than silently reformat it.
These checks execute no GENtle binary, Rust tests, VM, native GUI or model.
Session-close: 4 OK, 2 warnings (intentional working changes/unrelated `outputs/`
and the manual plan-fidelity reminder), 0 failures. Roadmap remains 400 lines.

## Claude-Reported Rust Verification At 284e5b70

Received through the owner on 2026-10-29: Claude reports running
`cargo test --locked --lib -- lifecycle_slice5 lifecycle_slice3` on macOS with
default features at `284e5b7049016e42ea19dd4bd50d7e93516d7ca7`: **28 passed,
zero failed**, comprising 17 Slice 5 tests (three GUI consumer unit tests) and
11 Slice 3 controls. `Cargo.lock` was reported unchanged. This is attributed
external evidence, not a Codex rerun; raw logs, toolchain and binary hashes were
not supplied in the forwarded comment. The three GUI tests are not native GUI
acceptance. See the [consolidated return](glen_next_verification_handoff.md#claude-reported-lifecycle-rust-results-at-284e5b70-2026-10-29)
for source review, formatting limits and outstanding checks.

The authoring-time unrun statements above/below are historical for this executed
subset. Current-source headless, other suites, VM/parity and native/platform
acceptance remain separate. Preserve the inconclusive long-job/cancellation
rehearsal, rather than treating these unit results as its completion.

## Remaining Glen/CI Verification

Owner-requested continuation, 2026-10-28: a further source pass at integrated
`a52c7b74` identified two missing test boundaries, not a confirmed runtime bug.
New synthetic regressions explicitly retain a read-only S1 journal append
while a preview finishes, and feed a published report unchanged into explicit
Apply after an unrelated edit, changed source DNA, changed annotations or an
occupied output ID. An unrelated edit can make a preview historical without
invalidating its actual source; the other three controls must refuse. Tutorial
06.07 now explains that distinction. Runtime, policy and scientific contracts
are unchanged. They were unrun at authoring; the reported 17-test Slice 5 return
above includes the two follow-up tests at its historical source.
Continuation checks pass 13 Python source/teaching cases with 12 process skips
(25 total), the targeted real LF/CRLF checkout test, scoped rustfmt and
documentation-status/whitespace checks. Session-close: 4 OK, 2 intentional
working-tree/manual-review warnings, 0 failures. No Rust/native result is claimed.

Use the ordinary integrated source and record full SHA/features/platform and
nonzero counts. Reuse unaffected same-source evidence; no new audit campaign.

```sh
cargo check -q --locked --tests
cargo check -q --locked --no-default-features --tests
cargo test --locked -j1 --lib command_execution:: -- --test-threads=1
cargo test --locked -j1 --lib routine_and_agent_assistant_ui::command_tests -- --test-threads=1
cargo test --locked -j1 --lib lifecycle_slice3_ -- --test-threads=1
cargo test --locked -j1 --lib lifecycle_slice4_ -- --test-threads=1
cargo test --locked -j1 --lib sequence_design -- --test-threads=1
cargo test --locked -j1 --lib operation_lifecycle -- --test-threads=1
cargo test --locked -j1 --lib agent_execution::tests::sequence_design -- --test-threads=1
cargo test --locked -j1 --test capability_registry_parity --test parity_matrix_freshness
cargo fmt --check
```

Claude reports scoped Slice 5 rustfmt passing but repository-wide formatting
failing with about 78 unrelated diffs at `284e5b70`. Do not treat that estimate
as current-source evidence or mass-format unrelated files. Retain the repository
gate, separate baseline drift from changed-file regressions and report scoped
formatting against the source actually checked; a global failure is not a pass.

Keep existing raw CLI/MCP/scripting scientific parity and source-matched
lifecycle/tutorial companions passing. Native checks: preview alongside a long
job, structural versus display edit, project replacement, delayed consumption/
sharing, visible snapshot warning and absence from project History. Preserve
normal report-with-path and Apply behavior. No live-model execution is required
for the deterministic host contract, and no native/biological verdict is claimed.
