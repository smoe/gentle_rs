# Lifecycle Slice 4: Preview Destination And Writer Safety

Status, 2026-10-27: **proposed, not implemented or Claude-reviewed**. The owner
requested a scoped plan for Claude, not a runtime change. Review/implementation
approval and Glen/CI verification are separate. No direct Claude call, Rust or
native execution, push or fetch. Source observations below predate the owner's
subsequent commit/rebase request; their SHAs are inspection provenance, not a
frozen-HEAD requirement. The previous goal remains paused.
Use the [forwardable review prompt](operation_lifecycle_slice4_claude_review_prompt.md).
The [parent lifecycle proposal](operation_lifecycle_policy_proposal.md#revised-slice-order)
and roadmap retain scope authority; this is not S2 or a new release gate.

## Baseline And Verified Evidence

Worktree `/Users/u005069/.codex/worktrees/faae/gentle_rs`, branch `gentle_rs_3`,
HEAD `f0b8371f4a99a455c60fc71a7769f3ea0380e7f3`. Local main has moved to
`a2637079f5a2ee1b49310f988243497562a65ccc` via three commits. Inspection found
no diff between these refs in the engine/handler/CLI/MCP/lifecycle policy files
listed below. Main's `a2637079` repairs test history calls to the existing
`undo_last_operation`/`redo_last_operation` API; its commit records 11 headless
lifecycle passes on macOS. That is source-scoped evidence recorded by another
thread, not execution here or Slice 4 acceptance. Before implementation, use
current main including that integration repair; do not reintroduce old calls.
Unrelated `outputs/` is untouched. No private data is an input to this plan.

| Source symbol | Observed boundary |
| --- | --- |
| `src/engine/ops/operation_handlers.rs::apply_internal`, Plan arm | Compute preview, then call the general JSON writer if `path` is present; only afterwards return/journal the successful operation. Cancelled reports can also be written. |
| `src/engine.rs::write_pretty_json_file` | Serialize pretty JSON, call `ensure_output_parent_dir`, then `std::fs::write`; existing files can be truncated and destination symlinks followed. Many unrelated operations use this helper. |
| `src/engine/ops/candidate_guides.rs::ensure_output_parent_dir` | Calls `create_dir_all`; even a later failed report write can leave directories. |
| `src/engine.rs::ProjectState::{load_from_path,save_to_path,candidate_store_sidecar_dir}` | Project filename is a host argument, not a persisted `ProjectState` field; loading hydrates candidate references and can replace them with warnings. Saving uses its own replacement/sidecar transaction. Do not change it in this slice. |
| `src/engine.rs::write_text_file_atomically` | Stage, flush/sync and `persist` with replacement: useful staging precedent, **not** a no-clobber writer. |
| `src/engine/io/{gel_images,gel_workbench}.rs` | Existing single-file exports stage with `NamedTempFile`, sync and `persist_noclobber`. Prefer that pattern over inventing a generic export framework. |
| `src/tss_profile_export.rs::{checked_directory,destination,is_parent_component}` | Checks raw traversal before joining/normalizing; includes Windows verbatim `Normal("..")`. Its fresh-directory policy is not a drop-in single-file rule. |
| `src/app.rs::current_project_path`, `set_current_project_path`; CLI/MCP state-path inputs | Hosts know file identity that plain engine/script state does not. Unsaved GUI projects and pure state-in/state-out scripts have no active project file. |
| `src/engine_shell.rs::ShellExecutionOptions`; `src/command_execution.rs` | No current protected-path context. Managed execution writes reports inside the worker before the guarded project/journal commit. |
| `src/engine/analysis/sequence_design/tests.rs::sequence_design_tutorial...` | Both strategies currently reuse one preview path. New no-clobber policy requires distinct fresh paths, not dropping either strategy/assertions. |

`tempfile = 3.27.0` is already a direct dependency. Its locally cached source
(`src/file/mod.rs::NamedTempFile::persist_noclobber`, Unix/Windows `imp::persist`)
explicitly promises no replacement but disclaims universal atomicity/security.
Unix can fall back to hard-link then unlink, leaving an extra temporary link;
Windows uses a non-replacing move. Inspect the resolved dependency again during
implementation. Do not turn this source inspection into a cross-platform pass.

## Recommended Scope And Behavior

Only migrate the optional output of `PlanDnaSequenceDesign`. No-path Plan and S1
remain file-free. Keep scientific request/report schemas, JSON bytes, approval
digest, search outcomes, progress, journals, checkpoints and Slice 3 construction/
autosave behavior unchanged. Writer error codes may become more specific.

**Default: a preview report creates one fresh ordinary file in an existing
directory, never replaces anything.** No overwrite/retry/delete-existing flag,
automatic filename, directory creation, backup system or filename extension rule.
If a report already exists, choose another explicit filename and obtain any
required command approval. Identical bytes are still an existing destination,
not permission to replace. Approval of a write attempt is not overwrite consent.

| Destination | Proposed outcome |
| --- | --- |
| Fresh ordinary leaf, existing eligible parent, outside protected paths | Stage exact JSON and publish with a no-clobber primitive. |
| Existing ordinary file, directory, symlink (including dangling), hard link or special file | Refuse without opening it for writing; retain its bytes and link identity. No separate same-file library is needed merely to reject all existing leaves. |
| Active project filename, even absent/empty, or its known sidecar root/descendants | Refuse using the host-bound protection context. Fresh-file-only alone does **not** protect an absent project path. |
| Missing parent, raw traversal/control characters, unsupported/ambiguous path semantics | Refuse before staging; no parent creation or guessed normalization. |
| Competitor creates the destination after preflight | Final no-clobber publication must lose safely; preserve competitor bytes. No delete/retry fallback. |

## Minimal Shared Boundary To Review

### Protected Project Binding, Not Payload Resolution

Propose one small engine-owned **nonserialized I/O context** for this writer.
It describes an explicitly file-backed host invocation (active project path,
protected candidate sidecar roots/locators, captured working-directory base), or
an explicitly stateless/unsaved invocation. Names and concrete API remain open
for review; do not build a general effect registry or permission framework.

- CLI/MCP use their actual state-path argument, including defaults and paths
  that do not exist. GUI uses its current project binding; unsaved means no file,
  not a fabricated default filename. Python delegates to CLI.
- Generic JS/Lua state-only calls are genuinely stateless: no project pathname
  can be inferred from DNA, metadata labels or a serialized report. They still
  receive universal no-clobber rules. File-loading wrappers that know a backing
  project must carry that binding, or state their exclusion and fail closed
  rather than claim project-alias coverage.
- Capture sidecar locator facts from the existing load, before hydration or
  permissive-warning fallback discards references; also protect the normal
  derived sidecar root used by a later save. Reuse loader/path rules, not an
  adapter JSON parser, filesystem crawl, second project read or private scan.
  Ask Claude whether a minimal internal load result/context is necessary or a
  smaller existing seam suffices. Do not persist context in project metadata.
- Bind/propagate context through direct apply, workflows, shared shell,
  forwarded commands and managed snapshots without changing their payload
  grammar, approval bytes, raw command hashes or resolution timing. No new
  operation-payload resolution layer. Do not eagerly read `@file` to classify it.
- Keep raw MCP confirmation ahead of context I/O/loading and all worker I/O
  outside GUI/admission locks. Test default/state/project aliases and host paths.
- Session context must not be imported from an untrusted request/report, mutate
  scientific digests, leak through project serialization, or be accidentally
  committed over a newer host binding by a detached result. Define its lifetime
  across clone/fork, project replacement, Save As, errors and unwinding. A sticky
  stale context or forgotten host adapter is a correctness defect, not a fallback.

The minimum guarantee covers the binding of the invocation, not unrelated files
or paths first selected by another actor afterwards. Ask Claude to challenge
Save As / concurrent binding changes before choosing a context API. A general
live binding registry or second artifact-commit phase is **not** pre-approved;
stop and narrow scope if needed rather than silently adding either.

### Path Admission And Alias Identity

Inspect raw components before absolute joining/canonicalization so `..` is not
erased; use the existing Windows regression as precedent. Reject empty, control,
drive-relative and unsupported device/stream forms with stable reasons. Do not
globally change shell escaping, trim into another path, or lowercase every path.

For a missing leaf, resolve/check the existing parent and compare the leaf in
that directory against protected names. Parent directory identity and native
case/normalization semantics matter: string equality, `canonicalize` on a missing
leaf, extension tests, or chromosome-style aliases are insufficient. Protect
known sidecar subtrees as subtrees, not just existing `index.json` files. Unknown
identity/unsupported filesystem semantics must refuse, not pass as non-alias.
Avoid promising a universal Unicode/case-fold implementation from source tests.

Choose the smallest demonstrably safe parent policy. Strictly rejecting all
symlink ancestors would reject common macOS `/tmp` and `/var` aliases; resolving
them can be acceptable only with an explicit race/identity boundary. Windows
junctions/reparse points are not covered merely by `is_symlink()`. Do not copy
the TSS helper wholesale or bless a canonical path without examining these cases.

### Stage And Publish, Not Truncate Then Repair

Use a Plan-specific shared writer; do **not** change `write_pretty_json_file`
for its unrelated callers. Validate destination/context before expensive search,
then recheck immediately before publication; both checks occur after approval
in the worker/execution phase. Missing-path Plan does not touch filesystem state.

Serialize with the same pretty-JSON bytes (including current newline behavior),
write into a uniquely created same-directory temporary file, flush/sync the
write-capable handle, then publish through a no-clobber primitive. Keep the
original error on failure; clean up only files owned by this invocation. Never
fall back to truncation, unlinking the destination, overwrite-rename, copy-over,
or treating an I/O error as successful report publication. Do not extend report
success into whole-workflow/managed-commit success.

Explicitly separate guarantees:

1. No replacement of an existing destination, including a final-name race.
2. Complete staged bytes before the final name is exposed on supported targets.
3. Temporary cleanup, which can fail or leave an extra link on some filesystems.
4. Crash/power-loss durability and parent-directory namespace races, which are
   not established by `flush`, file `sync_all`, canonicalization or `tempfile`
   alone. No universal atomicity or hostile-directory sandbox claim.

Rechecking paths is not a proof against parent swaps. If stronger adversarial
parent-directory protection is required, review handle-relative/native APIs
and dependencies separately or reject unsupported parents/filesystems. Do not
add bespoke unsafe cross-platform filesystem code under this small slice.

Cancelled search reports remain writable through this same approved safe path.
Writer refusal/failure returns an operation error and no new success journal;
preserve source state, revisions, redo and existing reports. Already published
files may survive a later workflow failure, service cancellation or stale commit,
as in Slice 3. No external rollback or search-status relabelling is added.

Keep deterministic machine-readable failure reasons in the existing
`EngineError.code`/`cause_chain`, not OS text or adapter-only strings. Proposed
categories distinguish invalid/unsupported path, protected path, existing
destination, unavailable identity/parent, staging write/sync and final publication
failure. Native error details remain diagnostic. Do not say every failed publish
leaves no final name: uncertain OS outcomes require inspecting only the owned
publication safely, never deleting a destination that may belong to another
writer. Ask review whether the existing error surface suffices; a new public
artifact receipt/schema is not assumed by this plan.

## Minimal Implementation Sequence After Review And Approval

1. Shared protected-path/context contract and fail-first Plan-only admission/
   writer tests. Confirm every required host has a bounded binding seam before
   changing behavior; implementation starts on then-current main, not this old
   branch. No new dependency unless review proves the existing stack inadequate.
2. Thin host propagation and staged no-clobber Plan writer together; test direct,
   named/forwarded/shared-shell, workflow, MCP, scripting and managed paths. Do
   not ship a GUI-only guard or an unrestricted raw-op bypass.
3. Reconcile protocol/CLI/GUI/decisions and tutorial 06.07 in the same behavioral
   change. Explain fresh outputs, existing-parent requirement, approval and scope.
   Adjust repeated synthetic previews to use distinct fresh paths, retaining both
   strategies and exact scientific/DNA/digest assertions. Changelog records work;
   roadmap tracks remaining verification, never a Done entry or new release gate.

## Fail-First Verification Matrix

Use temporary synthetic projects and sidecars with exact before/after bytes;
reuse tutorial 06.07's documented hand-crafted request. No private files or new
retained fixtures are necessary. Add fixture provenance if that changes.

| Boundary | Required controls |
| --- | --- |
| Existing destinations | Sentinel file, empty file, prior identical/different preview, directory, FIFO/special file where supported; bytes/metadata/link identities unchanged on refusal. |
| Alias protection | Project existing/empty/absent; absolute/relative/dot spelling, alternate parent alias, hard link, destination and dangling symlink; sidecar index/record plus unused new leaf under protected root. |
| Native path semantics | Windows drive-relative/UNC/verbatim traversal, ADS/device/reparse/case aliases; macOS case/Unicode variants and `/tmp`/`/var`; supported Linux aliases. Test actual filesystem behavior, not assumed platform names. |
| Sidecar loading | Normal and nonstandard reference roots, hydration, missing/corrupt references, strict/permissive loading; preserve project/sidecar bytes and original diagnostics. No global env races in parallel tests. |
| Bound context | Loaded and absent CLI/MCP state, unsaved GUI, stateless script; actual defaults, Save As/project replacement, snapshot context propagation and no stale binding overwrite. Never serialize context. |
| Staging failures | Inject serialize/write/flush/sync/publish failures; no truncated final file, no false success/journal, owned-temp cleanup and explicit residual on cleanup failure. Prefer deterministic test-only hooks, not sleeps or permissions-only failures. |
| Final-name races | Barrier between preflight and publication: another file/symlink/hard link wins; destination is preserved. Parent-swap test must either prove refusal or pin the expressly weaker supported-parent boundary, never pretend a second stat eliminates TOCTOU. |
| State/history | Plan with/without path; feasible/exhausted/infeasible/cancelled report; source/report/digest equality, unchanged project/undo/redo/dirty revisions; later mutation and explicit SaveFile/Apply retain legacy controls. |
| Workflows/service | Two fresh reports; duplicate target fails without replacing first; later failure/cancel/stale commit preserves published prefix; no new project/autosave/queue/observation semantics. |
| Agent/adapters | Existing report-write ask, wrapped/raw confirmation before I/O, denied MCP, inner-agent mutation/nesting guards; actual CLI/MCP/Python and JS/Lua bindings, LF/CRLF hash-bound project/input bytes and quoted paths. No model call required. |

Glen/CI owns all Rust/native execution. Proposed prefix `lifecycle_slice4_`;
require nonzero selected counts and record source/toolchain/features/platform.
After implementation, extend the existing actual-process lifecycle companion,
not a second tool runner. Suggested external groups:

```sh
cargo check -q --locked --tests
cargo check -q --locked --no-default-features --lib --bin gentle_cli --bin gentle_mcp
cargo test --locked --no-default-features -j1 --lib lifecycle_slice4_ -- --test-threads=1
cargo test --locked --no-default-features -j1 --lib sequence_design -- --test-threads=1
cargo test --locked --no-default-features -j1 --lib lifecycle_slice3_ -- --test-threads=1
cargo test --locked --no-default-features -j1 --lib command_execution -- --test-threads=1
cargo test --locked --no-default-features --features js-interface -j1 --lib lifecycle_slice4_ -- --test-threads=1
cargo test --locked --no-default-features --features lua-interface -j1 --lib lifecycle_slice4_ -- --test-threads=1
cargo test --locked -j1 --test capability_registry_parity --test parity_matrix_freshness
cargo fmt --check
GENTLE_LIFECYCLE_BIN_DIR=/absolute/same-source/bin python3 -m unittest scripts.test_operation_lifecycle_slice3 -v
GENTLE_TUTORIAL_BIN_DIR=/absolute/same-source/bin python3 -m unittest scripts.test_sequence_design_tutorial -v
```

Reuse same-source CI rather than repeat completed groups. Native Windows and
macOS filesystem checks cannot be inferred from LF/CRLF source tests or Linux
passes. Any missing privileges/filesystem support is an explicit unrun boundary.

## Review Decisions Still Needed

- Is fresh-file-only plus existing-parent policy the smallest useful default?
  No replacement option is proposed; challenge compatibility rather than silently
  restoring overwrite or automatic directory creation.
- What is the smallest nonserialized binding seam that covers absent projects,
  loader-lost sidecar references and all hosts without a payload-resolution layer?
  Are stateless calls and Save As/concurrent binding changes correctly bounded?
- Which parent/case/reparse policies are supportable with current dependencies?
  Is the stated cooperative-filesystem threat boundary sufficient, or must a
  smaller stronger implementation be chosen before promising alias safety?
- Does `persist_noclobber` plus explicit limited claims meet the needed writer
  contract? Do not equate no replacement, complete publication and crash durability.

Until those decisions and owner approval, this is a draft, not an adopted
filesystem invariant. S2, generic exporter migration, project-save/sidecar
transaction redesign, overwrite approval, opaque lifecycle reconciliation,
bounded queues and adversarial filesystem sandboxing remain outside scope.
