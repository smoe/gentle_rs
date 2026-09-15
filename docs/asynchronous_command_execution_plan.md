# Responsive Commands And Workflows

Status: reviewed plan; rollback/status prerequisites and the first managed
Agent Assistant command tranche implemented on 2026-09-15. Broader asynchronous
migration and live acceptance remain pending. The user supplied Claude's
review on 2026-09-14. The original proposal and review remain below;
the [reconciled Codex plan](#reconciled-codex-plan-2026-09-14) supersedes their
implementation ordering and identifies findings not supported by this checkout.

Baseline: `a5b7028157493d44c1812fc8e1f39e1ccbea6f6a` on `gentle_rs_3`, plus
uncommitted architecture/DEC-048/roadmap/changelog documentation. Unrelated
`outputs/` and any running genome preparation must remain untouched.

## At A Glance

1. **One way to submit.** A prepared Run button submits the same exact command
   as typed command input, visibly records it and preserves an unfinished chat
   draft. It does not ask the model to reinterpret the command.
2. **Quick work stays quiet.** Show its result directly. Longer work gets one
   activity with status, progress when known, and cancellation. Neither case
   blocks the GUI or prevents conversation.
3. **A workflow is also work.** Many quick steps must remain observable and
   cancellable as one activity, with completed steps distinct from saved results.
4. **Keep the safety rules.** Approval, result validation, project identity,
   undo and atomicity remain intact. Cancellation is a request until confirmed.
5. **Deliver in stages.** First fix confirmed rollback/status side effects and
   verify existing snapshot protections. Reuse the existing worker/status
   mechanisms, then prove genome preparation from Run and a long workflow of
   small steps. Migrate other command families afterward. No new biological logic.

## 1. User Contract

GENtle remains interactive while commands execute. A quick command displays
its result directly, without queued/completed message noise. A longer command
or workflow becomes one inspectable, controllable activity. The user can inspect
other windows, continue talking to the agent and cancel work. Neither human nor
agent must choose a special background command to keep GENtle responsive.

A prepared Run button and an explicitly typed command use the same submission
path. Both execute the exact GENtle command, not a new language-model
interpretation. Natural-language conversation remains a separate input mode.

The first acceptance case is genome preparation initiated by an Agent Assistant
button. The second is a long workflow consisting of individually quick steps.
Passing only the first does not establish architectural correctness.

## 2. Verified Starting Points

- `src/app/routine_and_agent_assistant_ui.rs`: Run calls
  `execute_agent_suggestion` during rendering. It converges with direct prompt
  commands in `execute_agent_shell_command_from_ui`, which takes the engine
  write lock, sets `progress_callback: None` and waits for execution.
- `agent_prompt_direct_shell_command` recognises slash commands and a few
  exact control words, not every valid shell command. Simply inserting an
  arbitrary command into the current prose prompt can route it to the model.
- The same UI disables prompt submission while `agent_task` exists. Model
  invocation and local command execution need distinct identities and controls.
- `src/engine_shell.rs`: expanded-stack workers are immediately joined. They
  solve stack depth, not asynchronous submission. `execute_parsed_workflow_command`
  serializes state before/after `apply_workflow` and waits for the entire result.
- `src/app.rs`: genome preparation already runs off-thread, with progress,
  cancellation and genome/catalog/cache/scope ownership in Background Jobs.
- `src/genomes.rs`: preparation has activity records, resource locks, duplicate
  run handling and component-level readiness. Reuse these contracts.
- BLAST has start/status/cancel/list routes, but `blast-start` probes tools and
  clones the full engine before returning. Its status path refreshes the
  scheduler and persists a store in engine metadata. It is not a fully
  compliant generic supervisor to copy unchanged.
- `src/engine.rs`: `clone_without_history`, `fork_detached_execution` and
  `commit_detached_execution` implement DEC-026. History-free snapshots still
  clone state/journal and require measurement; off-thread does not imply cheap.
- `src/agent_feedback.rs`: typed receipts distinguish dispatched/running from
  completed for existing BLAST envelopes. Model-facing feedback is bounded and
  redacted; job handling must not introduce raw paths, sequences or errors there.

## 3. Original Codex Plan

### A. One Submission Surface

Use one explicit submission function for prepared buttons, command-mode text,
GUI shell and adapter calls. Do not simulate keystrokes or execute shell strings
through an OS shell. Preserve existing parser and approval/recursion guards.

The Agent Assistant offers Message/Command modes using the same text control.
Run submits in Command mode and records the exact command visibly. Preserve an
unfinished message as a separate draft rather than silently overwriting it.
Existing slash/control shortcuts remain compatible. Editing a prepared command
invalidates any approval that no longer binds the exact request.

A short presentation grace period, initially about 250 ms and test-configurable,
suppresses job-indicator noise. It is not an execution deadline or a UI-thread
wait. If execution finishes quickly, show only its result. Otherwise show one
activity with elapsed time, phase, progress and inspect/cancel actions.

Only provably bounded in-memory admission/lookup work may run inline. No slow
operation is run on the UI thread with the hope of moving it after a threshold.
Submission IDs distinguish retries/double-clicks from intentional fresh runs.

### B. Shared Execution Controller, Separate Observation And Control

Introduce a headless shared execution controller below the GUI adapters. Keep
deterministic operation implementations reusable and synchronous inside their
worker; do not convert every biological primitive into its own task/thread.
The controller owns submitted activities, bounded scheduling, worker messages,
latest snapshots, cancellation tokens and terminal receipts.

Provisional public concept: a versioned command-execution envelope, separate
from the scientific Operation/request. Final schema/route names require review.
It binds execution/parent IDs, origin, project-instance or resource scope,
request identity, effective input/approval bindings, owner lifetime, state,
phase, snapshot revision/time, cancellation support and result/error references.
Large file resolution, snapshot construction, hashing and preflight happen in
explicit preparation phases, not synchronously before acknowledgment.

Separate running state from result availability and scientific readiness. States
must distinguish accepted/queued, preparing/validating, waiting or blocked with
reasons, running, cancelling, finalizing and terminal completed/failed/cancelled/
interrupted outcomes. Commit conflicts preserve a result with explicit stale
applicability; they are not successful application. These names are provisional.

Status reads use bounded cached snapshots and never refresh files, run probes,
or wait on the execution lock. Result retrieval is separate, reports not-ready,
and uses bounded projections or artifact references for large payloads. Direct
cancel targets the controller, not the busy command queue. A cancel request is
acknowledged promptly; actual termination remains cooperative/capability-bound.

GUI, shell, MCP, JS/Lua/Python consume the same controller contract. Existing
blocking CLI/script APIs may wait as explicit documented compatibility wrappers,
without blocking the shared controller. Actual async submission requires a live
owner (GUI, persistent MCP/shell/scripting session); a one-shot process must not
claim its worker survives exit. Detached services and cross-process control are
not prerequisites for this first fix and must not be implied by its job ID.

### C. Separate Resource Work From Project Mutation

Genome preparation uses the existing resource executor without taking a live
project write lock. Bind genome ID, source/catalog identity, cache and mode once;
changing the UI selection must not change ownership or progress attribution.
Reuse current preparation locks/readiness checks. Equivalent preparation can
attach to an existing compatible job, but identical scientific mutations are
not automatically deduplicated or reused without their existing authorization.

Read-only work uses the relevant bounded immutable snapshot. Project mutations
use DEC-026 detached execution and guarded commit, including a project-instance
identity so a result cannot apply to a different reopened project. Snapshot
creation, commit validation and dropping old large state must be inspected for
lock duration; a background worker holding a long write lock is still a defect.
Keep observation/control available even while preparation waits for a resource
or project lock. Serialize conflicting writers; do not enable blanket parallel
mutations merely because execution is asynchronous.

Keep lifecycle/progress receipts outside scientific project metadata. Preserve
existing operation journal/undo semantics for actual commits. A prepared genome
can remain valid when project attachment fails or is stale. Success receipts,
resource publication and partial-file cleanup retain existing validation rules.

### D. Workflows Stay One Managed Activity

The same controller handles `workflow`, `op`, macros and existing approved study
wrappers, rather than merely backgrounding their individual commands. Execute
the original ordered payload through existing engine paths; do not regenerate
plans, change digests, weaken transactional flags or infer new dependencies.

Expose the parent activity and child-step ordinal/operation/phase, counts of
computed versus durably committed results, pending steps and blockers. Show
one expandable workflow in the GUI, not hundreds of queued messages. A sequence
of cheap steps gets cooperative cancellation checks and scheduling opportunities
between steps; total runtime is not hidden by resetting a per-step grace period.
Preserve declared order and atomicity; this is not a new general DAG executor.

Cancellation stops launching further steps after the controller observes the
request and requests cancellation of the active step where supported. Explain
what committed state remains, what is rolled back and which external artifacts
remain. Never label uninterruptible work cancelled before it stops. Keep exact
failed/partial/rolled-back receipts without turning intermediate results into
durable scientific state. Existing approved checkpoint/reuse policies remain.

Progress is phase-local measured work, not an invented workflow percentage.
Step counts do not imply equal cost. Show last-progress age; ETA is optional,
qualified and derived only where the implementation has a meaningful basis.

### E. Conversation Is Not Execution

Local work must not occupy `agent_task` or disable new conversation. Feed current
bounded typed activity snapshots into the agent's context, with original
session/project ownership and omission metadata. A new prompt can ask what is
happening, request cancellation of a specific activity, or discuss unrelated
work. Exact status/cancel controls also remain available without any model call.

Keep model-request cancellation separate from command cancellation. Do not
silently enable overlapping model requests with unordered conversational
history; one active model response may retain its own explicit busy/cancel state.
Cancelling command work must remain possible even during that model response.
Ambiguous natural-language cancellation identifies/clarifies its target before
acting. Clearing conversation does not silently cancel or reassign executions.
No automatic recursive `agents ask` is introduced on completion events.

### F. Inventory And Migration

Seed an audit from existing parser/capability metadata and `docs/glossary.json`
(573 entries / 526 distinct catalog paths at the baseline, including bindings
and aliases; not proof of complete parser coverage). Map aliases to handlers
and operations, then compare against parser/dispatcher routes. Do not maintain
a second handwritten scientific command registry.

For each route record admission cost/I/O, execution path, locks, state/resource
effects, progress/cancellation, owner lifetime, result binding, migration verdict
and regression. Use explicit unknown/not-assessed rather than inferred safety.
Generic wrappers derive requirements from their actual contained operations.
Never classify solely by command name, duration observed on a tiny fixture, or
the presence of a worker thread. New routes must have a declared execution class
and test/justification; exact automated coverage mechanism needs review.

First priorities: genome/helper preparation; status/preflight and existing BLAST
async routes; downloads and resource indexing; primer/BLAST/RNA/CUT&RUN/array
work; rendering and large imports/exports; compound workflows. Quick inspections
also need proof that they do not wait behind execution. Unsupported background
paths remain explicit migration gaps, never silently run inline in a GUI handler.

## 4. Delivery Slices And Acceptance

1. Freeze the lifecycle/ownership and migration inventory after review. Introduce
   the shared controller plus fake workers, observer/control path and input-mode
   contract. No production genome or model needed.
2. Route prepared buttons and direct command submission through it for genome
   preparation, and adapt the existing preparation dialog to the same service.
   Preserve explicit payload/options and terminal preparation result semantics;
   do not require an LLM round trip or a special background command.
3. Add managed ordinary/transactional workflows through the existing executor,
   including many-small-step responsiveness, cancellation and guarded commits.
   This slice is required before claiming the reported architectural gap closed.
4. Adapt other command families in audited batches, starting with status/probes
   and BLAST. Expose consistent host-owned execution controls to other adapters;
   keep blocking compatibility wrappers documented. Do not call the universal
   contract implemented until inventory coverage and parity are demonstrated.

Deterministic tests use barriers/channels or an injected clock, not timing-only
assertions or real downloads. Include:

- Actual Agent Assistant Run-button click and typed Command submission produce
  the same bound execution (except origin), preserve drafts/approval and make
  no model call. Natural prose remains prose.
- A held worker cannot prevent acknowledgment, view changes, snapshot inspection,
  another conversation request or direct cancellation. Slow preflight/snapshot
  construction is exercised too, not only the main computation.
- Fast completion causes no transient job-message spam; a long chain of fast
  steps becomes one activity and remains cancellable before it finishes.
- Duplicate delivery versus intentional rerun, resource conflicts, explicit
  dependencies, failing children and existing atomic/non-atomic modes.
- Cancel during queued/running/finalizing phases; stale project-instance/structural
  baselines; project switch/close; owner death; orphan/interrupted status.
- Exact terminal results/journal/undo behavior match the synchronous reference
  execution under the same inputs. No completed claim precedes validated result
  publication/commit; unavailable progress/readiness is never a pass.
- Same lifecycle/receipt facts across GUI and persistent shell/MCP/scripting
  adapters; legacy CLI waits only on its own result. Separate one-shot process
  lifetime tests prevent false promises of detached execution.
- Offline GUI semantic checkpoints prove inspection and Cancel remain usable.
  A fake provider suffices for conversation responsiveness; no live inner-agent
  or external model tests are needed.

## 5. Explicit Non-Goals And Review Questions

No new biology, relaxed approvals, global automatic result cache, rewritten
workflow ordering/atomicity, server deployment, implicit restart/resume, or
interference with the running T2T preparation. A safe retry is not permission
to recompute or import prior scientific results silently.

Claude should check:

1. Is the proposed controller boundary compatible with root-engine extraction
   and current session ownership, or is a smaller existing runtime API better?
2. Which exact mutation/journal semantics must the genome-dialog/shared-shell
   convergence preserve, given their different present result paths?
3. Can DEC-026 snapshot/commit costs preserve responsiveness as proposed, and
   which existing UI readers could still block during these phases?
4. Is explicit Message/Command mode the smallest safe unification, including
   drafts, slash compatibility and prepared-command approval binding?
5. What is the minimal workflow cancellation hook that preserves existing
   atomicity and external-effect policy? Which wrappers need separate adapters?
6. Can existing BLAST/runtime/feedback contracts be extended without copying
   their known synchronous probe/clone/status-persistence limitations?
7. Which inventory coverage check and tests are necessary for the first two
   acceptance cases, and which requirements should be deferred explicitly?

## Review Record

The original Codex plan above is preserved for comparison. A read-only Claude
consultation was attempted using `-p --permission-mode plan --allowed-tools
"Read Grep Glob" --output-format text`. The first invocation returned no output
and was stopped; the network-enabled retry returned an HTTP 401 authentication
failure reporting an expired OAuth access token. No critique was received, no
credentials were inspected, and no reviewed/agreed schema is claimed.

Subsequently the user supplied Claude's review, appended below. The failed local
invocation is retained as history, not the current review status. Codex's
source-checked reconciliation follows the supplied review; neither review is
implementation or release acceptance.

### Additional Codex Check Before Implementation

Audit the entire command lifecycle, including admission and completion helpers,
not only `execute_shell_command_with_options`. In particular,
`finish_agent_shell_run` calls `import_agent_ensembl_gene_fetch_result`, which
performs another engine operation under a live write lock. Backgrounding the
initial fetch alone would leave this path blocking. State snapshots, approval
checks, imports and result formatting must receive the same bounded-work/lock
review; GUI-only presentation stays on the UI thread. Preserve the existing
fetch/import result distinction and authorization rather than hiding an extra
mutation inside a completion callback.

This adds one explicit audit/test obligation to slices 1 and 4; it does not
expand the first genome-preparation delivery into an Ensembl rewrite.

## Claude Review, 2026-09-14

Read-only review of the plan above against the tree. Reviewed at `main`
`48141673`, which has the plan's baseline `a5b70281` as an ancestor (four
subsequent `fix(genomes)` commits, none touching the paths below). Every claim
here cites a file and line; nothing was executed or modified.

The plan's diagnosis is right and its verified starting points check out. The
main disagreement is about size: three of the four things Section B proposes to
build already exist, and four correctness defects that async execution would
amplify are not in the plan at all. Those defects should be fixed before slice
1, not during slice 4.

### Q1 - Controller boundary, or a smaller existing API?

A smaller existing API. The plan would duplicate three subsystems:

- `src/runtime_status.rs` already defines the envelope Section B calls
  "provisional". `RuntimeStatusFrame` (`runtime_status.rs:158`) carries
  `frame_id`, `parent_frame_id`, `kind`, `label`, `phase`, `detail`, `state`
  in {Running, Waiting, Completed, Failed, Cancelled}, started/updated
  timestamps, `progress_percent`, `bytes_done`/`bytes_total`, `thread`.
  `RuntimeStatusActivity` (`runtime_status.rs:197`) adds `activity_id`,
  `source`, `scope` in {ProcessLocal, PersistedActivity, ProjectAsyncRegistry},
  `lifecycle_status`, `observation` in {Live, CrossProcess, Completed, Failed,
  Cancelled, Stale, Unknown}, `stale_reason`, `origin_process_id`. Snapshots are
  generation-versioned (`runtime_status.rs:398`) and taken without the *engine*
  lock - `snapshot_with_generation` still takes the registry mutex and clones
  active frames, so it is engine-lock-independent, not lock-free. `parent_frame_id` is already the parent/child structure Section D wants.
- `src/background_engine.rs:11` `execute_on_engine_snapshot` already implements
  Section C: fork under a read lock, run unlocked, commit under a short guarded
  write lock, drop the prior engine afterwards.
- `src/app.rs` already repeats one task handle eight times - `GenomePrepareTask`
  (`app.rs:1681`), `TutorialProjectTask`, `GenomeTrackImportTask`,
  `SequenceIngressTask`, `JasparBackgroundTask`, `DbSnpFetchTask`,
  `GenomeBlastTask`, `AgentAskTask` - each `{job_id, started, cancel_requested:
  Arc<AtomicBool>, runtime_frame: RuntimeStatusGuard, receiver:
  mpsc::Receiver<..>}`, over `BackgroundJobKind` (`app.rs:1364`) and
  `BackgroundJobEventPhase` {Started, CancelRequested, Completed, Failed,
  Retried, IgnoredStale} (`app.rs:1509`), with retry snapshots and stale-message
  filtering already handled.

Exactly two pieces are missing, and both are small. The registry holds no cancel
*request*: `RuntimeStatusGuard::cancel` (`runtime_status.rs:516`) only marks a
frame cancelled after the fact. And there is no result slot or terminal receipt
binding on a frame. Adding a cancel flag and a typed result handle to
`runtime_status` closes the gap without a second lifecycle vocabulary.

One hard constraint the plan does not state. `execute_reference_and_track_command_with_expanded_stack`
(`engine_shell.rs:64669`) and its three siblings pass `&mut GentleEngine` across
a thread boundary as a `usize` and rely on the immediate join for soundness.
The async boundary must therefore sit **above** `execute_shell_command_with_options`,
with the worker owning a detached engine outright; it can never move inside the
dispatcher. That is compatible with root-engine extraction - it is a call-site
rule, not a type rule - but it rules out a controller that owns worker messages
for individual operations beneath the GUI adapters.

On session ownership, the precedent not to copy is `BLAST_ASYNC_JOBS`: a
process-global `Mutex<HashMap>` hydrated from *project* metadata
(`engine_shell.rs:14886`), so identity is process-scoped while the store is
project-scoped. Two projects in one process share a namespace.

### Q2 - Which mutation/journal semantics must the convergence preserve?

The two paths are not two implementations of one contract; they are different
contracts, and the plan does not say which wins.

- Shell `ShellCommand::ReferencePrepare` (`engine_shell.rs:55031`) calls
  `engine.apply(Operation::PrepareGenome{..})`, which journals an
  `OperationRecord`, captures an undo checkpoint and bumps revisions
  (`engine.rs:30408`), and reports `state_changed: true`.
- The GUI dialog (`app.rs:12089`) calls the associated function
  `GentleEngine::prepare_reference_genome_once` (`engine.rs:9584`) - **no engine
  instance, no lock, no journal, no checkpoint** - synthesizes an `OpResult`
  with `op_id: "background-prepare-genome"` and empty `created_seq_ids`
  (`app.rs:12245`), and on completion calls only `invalidate_genome_genes()`
  (`app.rs:12527`).

So the dialog is a resource action with zero project effect; the shell command
is a journaled, undoable project operation. Converging them naively either makes
the Run button start dirtying the project and creating undo checkpoints where
the dialog never did, or drops the journal entry that `src/lua_interface.rs:1572`,
`src/bin/gentle_cli.rs:3726`, `src/workflow_examples.rs:3525` and roughly twenty
engine tests depend on.

Recommendation: do not converge the result paths. Background the *resource*
phase - which already needs no engine at all - and finish with the same
`engine.apply(Operation::PrepareGenome{..})` on the host thread as a short
guarded commit once the cache is populated. Whether that re-run is cheap against
a prepared cache must be measured, not assumed. Note also that the dialog's
`OpResult` is a hand-written literal of roughly 100 `None` fields
(`app.rs:12245`); any shared result type must not require a third copy of it.

### Q3 - DEC-026 costs, and which UI readers still block

The mechanism is sound, and `commit_detached_execution` is already careful to
leave evicted checkpoints in the prior engine so they drop outside the lock
(`engine.rs:7085-7113`). Two costs and one defect are unbudgeted; a third
claimed defect is **retracted** - see item 2.

1. **Fork cost is O(project + journal), and the journal is unbounded.**
   `clone_without_history` (`engine.rs:7024`) clones `state`, `journal` and
   `state.metadata` under the read lock. `history_limit` bounds only
   `undo_stack` (`engine.rs:10077`); nothing ever trims `journal`, and each
   `OperationRecord` stores a full `OpResult`. Every fork gets more expensive as
   a session runs. Measure this against a T2T-scale project before slice 2.

2. **RETRACTED - concurrent display and metadata edits are already
   preserved.** The original review claimed `commit_detached_execution` could
   silently discard them, because it validates staleness on `structural_revision`
   and `journal.len()` only. That reading stopped at the staleness check and
   missed the merge. Before the swap, `commit_detached_execution` calls
   `rebase_detached_commit_revisions` (`engine.rs:7186`), which carries the live
   `state.display` forward unconditionally (`:7235`), three-way merges
   `state.metadata` per key against the fork baseline - adopting a live-only
   change, rejecting a genuine both-sides conflict with a stale error
   (`:7194-7234`) - and patches live display and metadata into every undo
   checkpoint the detached execution created (`:7237-7258`). Revisions are then
   rebased onto live and bumped only per actual change class (`:7259-7269`).
   `src/background_engine.rs` covers each behavior directly:
   `detached_commit_preserves_concurrent_display_changes`,
   `detached_commit_preserves_concurrent_auxiliary_metadata_through_undo`,
   `detached_commit_rejects_conflicting_auxiliary_metadata` and
   `detached_workflow_preserves_live_display_in_each_new_undo_checkpoint`.

   The correction also explains the design the review misread as an oversight:
   mutation-only changes are *mergeable* and structural ones are not, which is
   exactly why the baseline check discriminates between them. Adding
   `mutation_revision` to that check would reject commits the rebase already
   handles correctly, turning an ordinary viewport change into a stale
   biological result. **Do not make that change.** Retain these protections and
   their tests as a regression boundary for every later slice.

3. **A BLAST poll that observes a job transition invalidates every in-flight
   detached execution.** *(RESOLVED in `c1eb0ef7`; both the insert and the
   remove now go through `auxiliary_metadata_mut()`, the byte-identical early
   return is intact, and a steady-state poll stays inert.)* `persist_blast_async_jobs_to_engine`
   (`engine_shell.rs:14854`) writes the job store through `state_mut()`, which
   bumps `structural_revision` (`engine.rs:7139`) - the one baseline the rebase
   in item 2 does **not** merge, and the exact value
   `commit_detached_execution` rejects on. Scope correction: the original review
   said "every poll". It does not. The helper returns early when the serialized
   store is byte-identical, and a still-running job's `try_recv` yields `Empty`
   without touching any status field, so a steady-state poll is inert. The bump
   occurs when a poll observes a transition - a job completing, failing, being
   cancelled or orphaned - or when `prune_blast_async_jobs_locked` drops a
   terminal job. That is a narrower window than claimed, but it is also the
   window in which long concurrent work is most likely to be mid-flight.

   Routing the write through `auxiliary_metadata_mut()` is the minimal fix, and
   item 2 sharpens why it is the right one: the job store is an ordinary
   `state.metadata` key, so at mutation level it would participate in the
   three-way merge rather than hit structural rejection. It is an interim fix,
   not DEC-048 compliance - a poll would still persist state and drive the
   scheduler. Separating observation from persistence and scheduling, per the
   plan's own Section C, is the real repair.

UI readers that still block: `agent_execution_revision()`
(`routine_and_agent_assistant_ui.rs:2830`, a read lock taken twice per command);
`blast_external_binary_preflight_report()` at `app.rs:12108` - a read lock *plus*
an external binary probe, run before the genome job is queued, which is the
"slow preflight before acknowledgment" DEC-048 forbids; and
`clone_without_history()` on the UI thread at
`routine_and_agent_assistant_ui.rs:2351` and
`main_area_dna/rna_read_mapping_ui.rs:5043` and `:5104`.

### Q4 - Is explicit Message/Command mode the smallest safe unification?

No. It is larger than what exists, and it would create the draft problem it then
proposes to solve. Mode is already *derived* from the text:
`agent_prompt_direct_shell_command` (`routine_and_agent_assistant_ui.rs:1834`)
inspects the single `agent_prompt` buffer, and the button relabels itself
"Run command" versus "Ask agent" (`:5957`). One buffer, one submit path, no
draft to preserve, no mode to desynchronize. An explicit toggle adds a second
buffer and a new failure mode - submitting in the wrong mode - for no safety
gain.

Prepared-command approval binding is likewise already a non-issue: suggestion
text is rendered as a read-only `egui::Label` (`:6188`), so there is nothing to
edit and no approval to invalidate. If editing is ever added, the pattern to
copy is `agent_screenshot_consent_binding_error` (`:274`), which already
invalidates on system id, project generation and originating-turn presence.

The two real defects are in the same code and are much smaller:

- `can_submit_prompt = !running && ..` (`:5945`) gates the **local command**
  branch on the model request being idle. That single expression is the Section E
  violation verbatim. The command branch should not consult `running` at all.
- The detector accepts only a leading `/` or the literals
  `capabilities|help|state-summary`. Widening it is the actual work, and it must
  stay conservative: prose that accidentally parses is worse than prose that
  reaches the model.

Smallest safe change: keep derived mode, drop `running` from the command branch,
widen the detector under test, leave the suggestion Run path alone.

### Q5 - Minimal workflow cancellation hook, and which wrappers need adapters

The hook already exists and is one call away. `apply_workflow_with_progress`
(`engine.rs:10177`) documents that "returning `false` requests cancellation",
and `ShellExecutionOptions::progress_callback` is already
`Arc<Mutex<Box<dyn FnMut(OperationProgress) -> bool + Send>>>`
(`engine_shell.rs:274`) forwarded by `forward_shell_progress` (`:57243`). Two
things stop it working:

- `execute_parsed_workflow_command` (`engine_shell.rs:64875`) calls the
  **non**-progress `engine.apply_workflow`, so no callback ever reaches a
  workflow.
- Neither `apply_workflow` nor `apply_workflow_with_progress` checks anything
  **between** ops. A workflow of five hundred fast ops that emit no progress is
  uncancellable by construction - which is exactly the plan's second acceptance
  case.

Minimal hook: emit one synthetic step-boundary `OperationProgress` before each
op in `apply_workflow_with_progress`, and on `false` stop launching further ops
and return the results computed so far as an explicit cancelled-partial outcome
rather than `Err`. Declared order and per-op atomicity are untouched; only
"launch the next op" is gated. Then switch `execute_parsed_workflow_command` to
the progress variant.

**Blocker before slice 3.** *(RESOLVED in `c1eb0ef7` - see the correction
below.)* The existing transactional rollback is unsound. `run_candidates_macro`
and `run_workflow_macro` roll back with
`*engine = GentleEngine::from_state(state)` (`engine_shell.rs:48991`, `:49086`),
and `from_state` (`engine.rs:6931`) builds with
`..Self::default()`. A rollback therefore erases the **entire** operation
journal and both history stacks, and resets execution, mutation and structural
revisions to zero - not just the macro's own effects. That is already a history
bug synchronously; with detached commits in flight it destroys the revision
monotonicity `commit_detached_execution` relies on for staleness detection. Fix
rollback to restore `state` while preserving journal and history and bumping
revisions forward.

*Correction and resolution.* The review also cited `engine_shell.rs:52477` and
`:67507` as rollback sites. They are not: `:67507` is `LoadProject` and `:52477`
is a sequence-pool import, where rebuilding from state is intended. Only the two
macro runners were the defect, and `c1eb0ef7` replaced both with one engine-owned
`with_rollback_on_error` guard (`engine.rs:7071`). Its `Drop`
(`engine.rs:6867`) restores state, journal and both history stacks, and sets each
revision to `max(current) + 1` with `op_counter` likewise `max`-ed - strictly
monotonic, so an in-flight detached fork is correctly rejected as stale rather
than silently passing. Being `Drop`-based, it also covers an early `?` and panic
unwinding. That is more complete than the review asked for.

`LoadProject` does remain relevant to a *different* point: `from_state` resets
revisions to zero, so a fork taken on a fresh project (structural revision 0,
empty journal) would pass the staleness check against a newly loaded project.
Revisions alone cannot carry project identity across a reopen. This is the
concrete case for the owner/project-instance binding in reconciled slice 1, not
a rollback issue.

Separately: `execute_parsed_workflow_command` computes `state_changed` by
serializing the whole project with `serde_json::to_value(engine.snapshot())`
twice, before and after, purely to diff (`engine_shell.rs:64879-64893`). That is
an O(project) cost inside the write-locked path, and worth removing.

Correction to the original review: it proposed comparing `mutation_revision`
instead and called the swap equivalent. It is not. The serialized diff reports
**net** state; the revision counts **edits**, including ones a later operation
reverses. A workflow that creates a sequence and then deletes it reports
`state_changed: false` today and would report `true` under a revision
comparison. That is a contract change for every `ShellRunResult` consumer, so it
needs an explicit decision and its own tests - and it is an optimization, not a
prerequisite for the first responsiveness fix. Sequence it behind the confirmed
defects.

Wrappers that need their own adapters: the macro runners, which re-enter
`execute_shell_command_with_options` per statement and already forbid nesting
(`engine_shell.rs:48984`); and the digest-verified study wrappers
(`verify_approved_gene_isoform_assay_study_workflow`, `engine_shell.rs:64907`),
whose `approved_workflow_sha256` / `operation_batch_sha256` binding must be
checked at admission **and** re-checked at commit, since approved bytes must
still bind the same request across an async gap.

### Q6 - Extending BLAST/runtime/feedback without copying their limits

Runtime and feedback: extend them. `runtime_status.rs` needs only a cancel-request
flag and a result slot; `agent_feedback.rs` receipts already separate dispatched
and running from completed.

BLAST: reuse the shape, not the code. Four limitations must not carry forward.

1. `blast-start` (`engine_shell.rs:55258`) runs
   `blast_external_binary_preflight_report()` - an external tool probe - and then
   `engine.clone()`, a **full** clone including both history stacks rather than
   `clone_without_history`, both before acknowledging. Both belong after
   admission, and the clone should be `fork_detached_execution`.
2. The job store lives in project metadata through `state_mut()`, so status
   polls dirty the project and bump structural revision (Q3.3).
3. `with_blast_async_registry` (`engine_shell.rs:14886`) runs refresh, dispatch,
   prune and persist on **every** status read. "Status" is a scheduler tick and a
   write, not the bounded snapshot read DEC-048 requires.
4. `BLAST_ASYNC_JOBS` is process-global but hydrated per-engine (Q1).

### Q7 - Necessary coverage and tests, and what to defer

Necessary for the first two acceptance cases:

- Parity: the derived-command detector routes `/genomes prepare ..` identically
  from the Run button and from typed input, makes no model call under a fake
  provider, and is not gated by `running`.
- Held worker, using a barrier rather than timing: acknowledgment,
  `runtime_status_snapshot` and cancel all work while the worker is blocked -
  run once with the barrier in preflight, not only in the main computation.
- Genome-prepare equivalence: async submission yields the same journal entry,
  undo availability and `OpResult` as `engine.apply(Operation::PrepareGenome{..})`
  under identical inputs.
- Regressions for the two confirmed defects: a BLAST poll that observes a job
  transition does not invalidate an in-flight detached execution; transactional
  macro rollback preserves the journal, both history stacks and monotonic
  revisions. No new display/metadata test is needed - the four existing
  `src/background_engine.rs` tests named in Q3 item 2 already cover it, and they
  are the regression boundary that must keep passing through every slice.

`docs/glossary.json` is a reasonable seed - 573 commands over 526 distinct
catalog paths, confirmed - but it is a documentation inventory, not a parser
inventory, and cannot prove coverage. Use it to enumerate, then reconcile
mechanically against the `ShellCommand` variants.

Defer explicitly: cross-process and detached-service lifetime, restart/resume,
a general DAG executor, migration of the roughly seventy
`execute_shell_command*` call sites across fifteen files (`lua_interface.rs` 21,
`mcp_server.rs` 13, `main_area_dna.rs` 10, `bin/gentle_cli.rs` 7, and the rest),
and any automated enforcement that new routes declare an execution class.

### Revised Minimal Plan

Five slices. Slice 0 is new and is the main change to the delivery order: these
are existing defects that asynchrony amplifies from latent to routine, and each
is small enough to land before any new machinery.

0. **Prerequisites, no async.** Make transactional macro rollback preserve the
   journal, both history stacks and monotonic revisions. Stop BLAST job-store
   persistence from advancing the structural revision. Measure
   `clone_without_history` against a synthetic project and a grown journal, and
   record size, build and methodology - not against the running T2T study.
   Run the existing detached display/metadata/undo tests as the entry gate.
   **Not** in this slice: the retracted `mutation_revision` baseline change
   (Q3 item 2), and the workflow `state_changed` optimization, which needs its
   own contract decision (Q5).
1. **Generalize, do not invent.** Fold the eight repeated task handles into one
   `BackgroundTask<T>`; add a cancel-request flag and a result slot to
   `runtime_status`. No new envelope schema, no second lifecycle vocabulary,
   fake workers only.
2. **Genome preparation from the Run button.** Background the resource phase,
   which already needs no engine; keep the journaled `apply` as a short guarded
   commit on the host. Move the makeblastdb preflight after admission. Drop
   `running` from the local-command gate. Keep the dialog and shell result paths
   distinct and document the difference.
3. **Workflows as one activity.** Step-boundary progress and cancellation in
   `apply_workflow_with_progress`; cancelled-partial as a distinct outcome from
   `Err`; parent and child frames through the existing `parent_frame_id`.
4. **Audited migration batches.** Status and probe routes first, then BLAST.
   Blocking wrappers stay, and stay documented.

Sections A through F of the original plan remain the correct requirements. The
disagreement is scope and order, not direction: most of Section B exists, Section
C's engine primitive exists, Section D's hook exists but is unwired, and Section
E is one boolean. Three correctness defects are missing from the plan entirely
and should lead; a fourth claimed defect is retracted in Q3 item 2.

## Reconciled Codex Plan, 2026-09-14

This section supersedes the original delivery sequence, not the user contract.
Reuse existing mechanisms and address confirmed correctness prerequisites first.
No runtime changes are authorized or claimed by this documentation update.

### Findings Accepted And Corrected

- **Reuse is the right direction.** `background_engine::execute_on_engine_snapshot`
  already owns detached execution, and `runtime_status` already provides active
  frames, generation snapshots and parent IDs. Extend them rather than introduce
  a competing lifecycle. Factor task-handle boilerplate only where the first
  migrated paths need it; converting every existing GUI task is not a prerequisite.
- **The claimed display/metadata loss is not present as described.** Before
  swapping, `commit_detached_execution` calls `rebase_detached_commit_revisions`
  (`src/engine.rs:7186`). That helper preserves live display state, three-way
  merges metadata, rejects conflicting keys and rebases new undo checkpoints.
  Existing `src/background_engine.rs` tests exercise each behavior. The engine
  source blob is identical at this checkout's `a5b70281` and Claude's reviewed
  `48141673`; this is not explained by a later fix. Do not add a blanket
  `mutation_revision` rejection: that would turn ordinary viewport changes into
  stale biological results, contrary to DEC-026. Retain these protections/tests.
- **BLAST polling can invalidate detached results.** When its serialized job
  store changes, `persist_blast_async_jobs_to_engine` uses `state_mut`, advancing
  the structural revision. An unchanged poll does not. Moving writes to the
  auxiliary accessor is an interim fix for structural invalidation, not complete
  DEC-048 compliance: polling would still persist state, schedule work and risk
  metadata conflicts. Separate observation from persistence/scheduling, preserving
  existing persisted-job recovery until a documented replacement exists.
- **Rollback and cancellation findings are confirmed.** Transactional macro
  failure reconstructs the engine with `from_state`, discarding the prior journal,
  history and revisions. Workflow shell execution omits the existing progress
  callback, and the engine's progress variant lacks between-operation checks.
  Repair these before expanding asynchronous mutation.
- **Runtime snapshots are engine-lock-independent, not lock-free.**
  `snapshot_with_generation` takes the registry mutex and clones active frames.
  Frames disappear when their guards drop; parent inheritance is thread-local.
  A cancel flag/result slot alone does not supply execution ownership, approval
  binding, terminal retention, admission limits or explicit cross-worker parents.
  Add only those missing pieces in a host-owned execution record linked to the
  existing frame. Keep bulky results outside the observation lock.
- **Keep the synchronous stack helpers synchronous.** Their pointer transfer is
  safe only with the immediate join. The asynchronous boundary sits above
  `execute_shell_command_with_options`, in a worker owning its execution state;
  never remove the inner join to obtain background behavior.
- **Keep the existing derived input mode.** An explicit Message/Command toggle
  is unnecessary for the initial fix. Slash commands are unambiguous; Run can
  supply a typed exact-command submission directly to the same handler without
  inserting into or overwriting the prompt buffer. Keep the visible command
  record and suggestion guards. Do not broaden prose detection merely because
  text parses as a command. Relax only the model-busy condition for local
  commands; this alone does not eliminate their blocking execution or the locks
  needed to start a conversation while work runs.
- **Preserve distinct genome completion semantics.** The dialog prepares a
  resource without a project journal entry; shell `PrepareGenome` journals an
  operation. Share resource execution but retain these caller contracts. Do not
  run `apply(PrepareGenome)` on the UI thread a second time assuming a warm cache
  makes it bounded. Validate resources off-thread; use an internal shared
  finalization path for exactly-once journal/history handling, or retain fully
  detached execution until that split is proven safe. Normal `apply` must use
  the same result construction and validation rules, not an adapter-made result.
- **Do not turn cancellation into successful partial execution.**
  `Ok(Vec<OpResult>)` currently means workflow success. A managed outcome must
  distinguish cancelled/failed/completed and computed/committed/rolled-back
  children. Preserve existing legacy error behavior; expose partial receipts
  separately. Approval wrappers must not mistake cancellation for success.
- **Revision comparison is an optimization requiring a contract check.** The
  current workflow `state_changed` compares net serialized state. Revisions
  record edits, including changes later reversed. Replacing the former with the
  latter is not automatically equivalent. Prove or explicitly document semantics
  before changing it; it is not a prerequisite for the first responsiveness fix.

### Revised Delivery Sequence

0. **Correctness and measurement.** Run the existing detached display/metadata/
   undo tests. Add regressions before repairing transactional rollback: retain
   the pre-transaction journal and undo/redo baseline, discard only transaction
   effects, keep revisions monotonic, and reject stale in-flight commits. Cover
   parsing/nesting errors and cancellation as well as execution failures; an
   early `?` must not escape a transaction after earlier statements mutated it.
   Keep external-file/resource side effects explicit, not falsely rolled back.
   Repair BLAST observation-induced structural mutation with durability tests.
   Measure fork/commit/lock cost against a realistically sized, non-private
   synthetic project and growing journal; record size, build and methodology.
   Do not inspect or reuse the running T2T study for this benchmark.
1. **Connect existing infrastructure.** Build one headless admission/control
   service around the existing runtime frames, task channels/cancel tokens and
   detached executor. Bind request, owner/project instance, approval and result
   identity; retain terminal receipts after frame removal. Explicitly propagate
   parents across workers. Bound queue/result retention and expose cached
   status/cancellation independently of execution locks. Prove this with fake
   workers before wiring biology. No second vocabulary or sweeping task refactor.
2. **First visible fix: genome preparation.** Prepared Run and typed slash input
   share exact submission, with existing approval/recursion guards and no LLM
   round trip. Preserve prompt drafts. Move tool probes/resource checks after
   acknowledgment, preserve dialog versus shell journal semantics, and keep
   completion processing off the UI thread except bounded presentation/commit.
   Quick results remain quiet; long work exposes an activity. Prove navigation,
   status, another conversation and direct Cancel while preflight/work is held.
3. **Second required case: long workflows of small steps.** Forward progress
   through the existing workflow executor and add between-step cancellation
   checks before launching the next operation. Use explicit managed terminal
   outcomes and parent/child frames without altering approved ordered payloads.
   Test atomic and non-atomic failure/cancellation, including the commit race:
   cancellation observed before publication prevents commit; a completed commit
   is reported as completed, not retroactively labelled cancelled. Freeze and
   validate approved input bytes off-thread; perform bounded binding/revision
   checks at publication without re-reading large files under a commit lock.
4. **Audited migration.** Status/probe paths, BLAST and the remaining families
   follow in batches. Include admission, snapshots and completion helpers in
   each audit, not just the central executor. Preserve blocking CLI/script
   compatibility as caller-side waiting. Record unassessed routes explicitly;
   broad adapter migration, cross-process services and restart/resume remain
   separate from the first two acceptance cases.

### Required Acceptance Evidence

Retain section 4's controlled-worker tests, with these clarifications:

- Run versus typed slash command binds the same request without touching an
  existing prose draft or invoking a provider; local control also works during
  a held model response. Use a fake provider, not a live agent test.
- Existing detached display/metadata preservation and conflict tests remain
  green; no new blanket mutation-revision rejection is introduced.
- Transaction failure preserves pre-existing undo/redo and journal, advances
  invalidation identity monotonically, and leaves no aborted child as committed.
- Changed BLAST job observations do not invalidate a detached scientific result;
  cached status reads do not mutate project state or perform scheduler/I/O work.
  Recovery of existing persisted BLAST jobs is separately verified.
- Resource work executes once; shell journal/result semantics match synchronous
  execution while dialog-only preparation does not acquire project side effects.
- Cancellation during a sequence of non-progress-emitting operations has an
  explicit terminal outcome; legacy callers never receive ordinary partial success.
- Runtime frames retain correct parents across worker threads; terminal receipt
  availability, owner/project scope, bounded retention and stale-result rejection
  are tested independently of whether a frame is still visible.

Implementation status is recorded below; the historical reconciliation was planning only.

### First Managed-Command Tranche, 2026-09-15

The initial service uses runtime frame states and history-free detached execution,
not a second execution engine. Admission performs bounded text hashing and a
nonblocking owner/revision check; snapshot copying, probes and execution run on
workers. The admitted structural revision and journal length are checked again
before copying, so a queued worker never silently adopts later sequence edits.
Project-instance identities close the reopen/revision-zero hole without rejecting
mergeable display/metadata edits. Status/cancel use an independent control lock;
result publication and cancellation share a final guarded boundary. Terminal
receipts retain request/output hashes and computed versus committed counts.

The first GUI migration covers `genomes prepare` (including its helper alias),
`op`, `workflow`, workflow/candidate macros/templates, and the verified single
gene-study workflow route. Run and typed submission use the same path; the draft
is untouched. Direct commands no longer depend on the model being idle. Progress
and Cancel appear after a 200 ms quiet interval. Completion retains its original
turn/session, not the conversation turn current when the worker finishes.
Its state-change flag describes that command's result, not unrelated edits made
while it ran. Running status text stays bounded even for a large JSON workflow.

Workflow and macro boundaries now check cancellation and expose counts. Managed
errors/cancellation discard the detached engine delta; standalone synchronous
nontransactional execution keeps its prior prefix semantics. External file,
resource and independently running child-job effects are not transactional.
The shell still journals preparation once; dialog-only resource preparation is
unchanged. No warm-cache re-execution is performed on the UI thread.

Limits: four active workers (no implicit waiting queue), 32 retained records,
1 MiB command text and 16 MiB interactive result. Busy admission is explicit;
consumed terminal receipts may be evicted at capacity. This is process-local,
not restart recovery or a new CLI job server. Concurrent structural changes can
still reject a result, including changes made by another admitted worker.

Remaining: pure cached BLAST observation, BLAST-start full-clone/probe ordering,
audited batch/per-operation failure receipts, other GUI command families and
adapters, cross-process control/recovery, and live T2T acceptance. The synthetic
1-Mbp/growing-journal cost probe is not a bound for real annotated projects.
No live model, private genome preparation or production report is used here.

#### Focused Verification

The first-tranche working tree passed 45 focused tests (one timing probe ignored
in that suite); the probe then passed separately. Reproduction commands:

```bash
cargo test -q -j 1 --lib -- command_execution::tests:: routine_and_agent_assistant_ui::command_tests:: background_engine::tests:: workflow_progress_ transactional_rollback_ shared_history_rebase_ blast_async_store_ agent_prompt_direct_shell_command_ runtime_status::tests:: --test-threads=2
cargo test -q -j 1 --lib detached_execution_snapshot_cost_probe -- --ignored --nocapture
```

The GUI tests exercise shared Run/prompt handlers, held-worker inspection,
cancellation and project replacement programmatically, not live desktop clicks.
The local preparation test uses generated FASTA/GTF files, missing-tool overrides
and two cold/warm runs; it does not download or prepare a production genome.

The debug-build probe on this development machine used a synthetic 1-Mbp state
and three samples per journal size. Concurrent builds were present; these are
observations, not a CI timing threshold or representative biological benchmark.

| Journal records | Snapshot/read-lock range (ms) | No-op commit-lock range (ms) |
| --- | --- | --- |
| 0 | 0.051-0.126 | 0.009-0.094 |
| 100 | 2.022-2.595 | 0.041-0.057 |
| 1,000 | 19.419-24.922 | 0.048-0.094 |
| 10,000 | 189.222-271.554 | 0.050-0.053 |

Journal-copy cost is demonstrably nonconstant and remains an optimization target.
This probe excludes expensive real annotation payloads and nontrivial commit
rebasing; it does not close the live responsiveness acceptance gate. The linker
emitted the known macOS large `__eh_frame` warning but exited successfully.

### Reconciliation Verification

Checked on `a5b7028157493d44c1812fc8e1f39e1ccbea6f6a`, with documentation edits only:

- `cargo test -q --lib background_engine::tests::`: 16 passed, including live
  display preservation, metadata merging through undo and conflict rejection.
  The macOS linker emitted its large `__eh_frame` warning; linking and tests
  still completed successfully.
- `cargo check -q`: passed.
- `git diff --check`: passed.
- `python3 scripts/maintenance_chore.py session-close --plan
  docs/asynchronous_command_execution_plan.md`: 4 OK, 2 warnings, no failures.
  Warnings are the intentional uncommitted documentation/unrelated `outputs/`
  and the manual plan-fidelity reminder. Scope remains planning, not runtime work.

Rollback, BLAST polling and workflow cancellation findings were inspected in
source, not fixed or exercised by new regressions in this reconciliation.
No private data, live model, genome preparation or production report was run.

## Claude Cross-Read Of The Reconciled Plan, 2026-09-14

Read-only cross-read of the corrections above, performed before implementation.
The Claude review section has been corrected in place: Q3 item 2 is retracted,
Q3 item 3 is narrowed, Q1's "lock-free" wording is fixed, Q5's revision-comparison
claim is withdrawn as an equivalence, and the revised Slice 0 no longer contains
the retracted change. Both earlier proposals remain intact.

**Corrections to Claude, verified and accepted.**

- *Display/metadata preservation.* Accepted without reservation. The original
  finding stopped at the staleness check in `commit_detached_execution` and never
  read `rebase_detached_commit_revisions` (`engine.rs:7186`), which is where the
  protection lives. The four named `src/background_engine.rs` tests exercise it
  directly. The design is also deliberate rather than incidental: mutation-only
  changes are mergeable and structural ones are not, which is precisely why the
  baseline check discriminates. A blanket `mutation_revision` rejection would
  regress working behavior and is withdrawn.
- *"The engine source blob is identical at both commits."* Verified:
  `a5b70281:src/engine.rs` and `48141673:src/engine.rs` are both blob
  `97124eb727423ae3aa00799b2cc876bb877e6c70`. The mistake was a reading error,
  not a version skew, and the review was corrected rather than excused.
- *"An unchanged BLAST poll does not advance the structural revision."*
  Verified. `persist_blast_async_jobs_to_engine` returns early on a
  byte-identical store, and a still-running job's `try_recv` yields `Empty`
  without touching any status field. "Every poll" was wrong; the bump needs an
  observed transition or a prune. Q3 item 3 now says so.
- *"Runtime snapshots are engine-lock-independent, not lock-free."* Verified;
  `snapshot_with_generation` takes the registry mutex and clones active frames.
  The distinction matters for the admission service's observation path.
- *"Revision comparison is an optimization requiring a contract check."*
  Accepted, and the strongest of the corrections. A serialized diff reports net
  state; a revision counts edits including reversed ones, so create-then-delete
  changes answer. Treating the swap as equivalent would have altered
  `ShellRunResult` semantics for every consumer. Correctly resequenced out of
  the prerequisite slice.

**Cross-read findings on the reconciled plan itself.**

- The delivery sequence is sound and nothing from the review that survived
  verification was dropped. Slice 0's "growing journal" benchmark picks up the
  unbounded-journal fork cost; the expanded-stack pointer constraint, derived
  input mode, distinct genome journal semantics, unwired workflow progress and
  cancellation-is-not-success all carry forward intact.
- Two review details are not explicit in the reconciled slices and should not be
  lost in slice 4's BLAST batch: `blast-start` clones with `engine.clone()` -
  a **full** clone including both history stacks - where `fork_detached_execution`
  is the correct primitive (`engine_shell.rs:55260`); and it runs
  `blast_external_binary_preflight_report()`, an external tool probe, before
  acknowledging (`:55258`). Both are the same admission-ordering defect slice 2
  fixes for genome preparation.
- The two named safeguards are the right ones to hold. "A warm cache does not
  justify blocking execution" is the correct reading of the genome-completion
  risk, and is stricter than the review's "measure it" - that is the safer
  default. "Cancelled partial workflows must never appear successful" is load
  bearing given that `Ok(Vec<OpResult>)` currently means success, and it is what
  makes the between-step cancellation hook safe to add.

No further disagreement. The prerequisite slice is well-formed and the two
confirmed defects - transactional rollback and BLAST observation-induced
structural mutation - are the right place to start.

## Implementation Progress (2026-09-15)

The first two correctness prerequisites are implemented, not the complete
asynchronous lifecycle:

- Both transactional macro runners use one engine-owned rollback guard. It
  restores pre-run state, journal, undo/redo and history limits on error or
  unwind, advances revisions and retains the larger baseline/final operation
  counter. Parsing/nesting failures after an executed prefix no longer
  escape rollback. Public workflow failure lineage remains outside rollback;
  external resource/filesystem/job effects are explicitly not rolled back.
- Private history checkpoints use immutable shared ownership so retaining the
  transaction baseline does not deep-copy all historical projects. Undo/redo
  can still consume checkpoints, including within a transaction, and detached
  rebasing uses copy-on-write. This is the only additional internal ownership
  change required here; persisted formats and public operation payloads do
  not change. Current state and the growing journal still need cloning.
- Changed/pruned BLAST stores use auxiliary metadata mutation. Identical
  polls remain no-ops; real transitions still dirty the persisted state and
  invalidate redo, but do not bump structural revision. Existing three-way
  metadata merging preserves live status through detached commit and undo;
  genuine conflicts remain rejected. Restart recovery is unchanged.

The new regressions first reproduced loss of an older journal on execution
error and retention of an executed prefix on parse error. Coverage now also
includes forbidden nesting, macro undo/redo, checkpoint eviction, successful
transactions, nontransactional failure, panic unwinding, failure receipts,
BLAST transition/prune/no-op and restart interruption/cancellation.

Still pending: the growing-journal snapshot-cost benchmark, pure cached BLAST
observation separated from scheduling/persistence, reusable activity control
and result delivery, genome Run/slash submission and workflow cancellation.
BLAST-start still needs a history-free fork and tool probing moved after
admission. Shared checkpoint ownership does not make its full-engine clone an
appropriate worker boundary. No UI responsiveness or live T2T acceptance is
claimed by these prerequisite fixes.

Verification for this slice:

```bash
cargo test -q --lib -- transactional_macro_ transaction_rollback_ \
  shared_history_rebase_ background_engine::tests:: blast_async_store_ \
  macro_transactional execute_async_blast_ execute_macros_run_records \
  execute_macros_run_failed execute_macros_template_run_records \
  engine::tests::test_engine_history engine::tests::test_display_history \
  engine::tests::test_data_mutation_history \
  engine_shell::tests::execute_history_commands --test-threads=1
cargo check -q
cargo fmt --all -- --check
git diff --check
scripts/maintenance_chore.py session-close --plan docs/asynchronous_command_execution_plan.md
```

- Focused suite: **41 passed**, including all 16 existing detached-engine
  protections and the shared-checkpoint isolation regression. Synthetic inputs
  only; no private study, genome preparation, model call or GUI acceptance run.
- Check, formatting and whitespace gates pass. The macOS linker still warns
  about the large `__eh_frame` section; linking and the unwind test succeed.
- Maintenance: 4 OK, 2 warnings, 0 failures. Warnings are the intentionally
  uncommitted worktree (including untouched unrelated `outputs/`) and manual
  plan-fidelity review. Only the prerequisites above are implemented; existing
  planning edits and Claude's corrected review were preserved.
