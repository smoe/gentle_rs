# Lifecycle Slice 4: Preview Destination And Writer Safety

Status, 2026-10-27: **owner authorized the reconciled implementation; authored,
with Rust/actual-process/VM/native verification pending with Glen/CI**. The
Plan writer, host binding, tutorial and tests now follow this contract. No new
scientific schema or payload resolver. No direct Claude invocation, local Rust/
native execution, push, fetch or goal resumption. The earlier goal remains paused.
The [original review prompt](operation_lifecycle_slice4_claude_review_prompt.md)
is historical; this document owns the revised minimal contract. The
[parent lifecycle proposal](operation_lifecycle_policy_proposal.md#revised-slice-order)
and roadmap retain scope authority. This is not S2 or a new release gate.

## Source And Evidence

Review reconciliation on branch `gentle_rs_3`, worktree
`/Users/u005069/.codex/worktrees/faae/gentle_rs`, HEAD
`b8e0118feb44726392702398be6b6d52218eec66`; local main observed at
`a6ba0da6053cd6e514021b86d38f2169429fd400`. The original proposal at `b8e0118f`
records inspection of `f0b8371f` and newer main `a2637079`; those are historical
source facts, not a frozen-HEAD requirement. Main's Slice 3 test-API repair is
already an ancestor here. Its recorded 11 headless lifecycle passes are not
Slice 4 execution or scripting/native acceptance. Unrelated `outputs/` is untouched.
During implementation the owner requested another rebase: the branch advanced
to main `1804605f1ba42369f943213414d79fd9a40013cb` and restored the scoped
Slice 4 working changes without conflicts. The new RNA S1 changes are retained;
their recorded passes are not Slice 4 verification. A recovery stash
`617a89ee3d08cf57ecde2e5866b279507b9bc44d` retains the pre-rebase work only;
unrelated outputs and older stashes were not applied or changed.

The table records the pre-implementation inspection, not current writer behavior.

| Source symbol | Pre-implementation relevance |
| --- | --- |
| `src/engine/ops/operation_handlers.rs::apply_internal`, Plan arm | Compute preview, then call the shared JSON writer before successful journaling. Cancelled reports can be written too. |
| `src/engine.rs::write_pretty_json_file`; `src/engine/ops/candidate_guides.rs::ensure_output_parent_dir` | Current writer serializes pretty JSON, creates parents and uses truncating `fs::write`. Many unrelated operations use it; do not change their behavior here. |
| `src/engine.rs::ProjectState`, `src/cli_support.rs::load_state_or_default` | Missing/empty state paths load as default. Project JSON requires `sequences` and `metadata`; a preview saved there is not a valid project. |
| `ProjectState::{candidate_store_sidecar_dir,prepare_candidate_store_sidecar_transaction,commit_candidate_store_transaction,finalize_candidate_store_transaction}` | Saving owns the derived `<project filename>.candidates` root and can replace/remove that entire tree, including otherwise fresh reports placed inside it. This root depends only on the project path. |
| `src/engine.rs::GentleEngine`, `clone_without_history`, `fork_detached_execution`, `commit_detached_execution`, `rebase_detached_commit_revisions` | Session fields can be serde-skipped, but detached snapshots use an explicit struct literal, not automatic derived Clone. Commit swaps engines and currently preserves live display/metadata, not an as-yet nonexistent host path. |
| `src/app.rs::current_project_path`; CLI/MCP state arguments | Hosts know the backing/output project filename. Unsaved GUI and pure script state have no such file. |
| `src/engine/io/{gel_images,gel_workbench}.rs`; `src/engine/analysis/rna_panel/tools.rs` | Existing small staged-file/sync/no-clobber patterns to reuse, without a general export framework. |
| `Cargo.lock`; cached `tempfile 3.27.0` documentation/implementation | `same-file 1.0.6` is locked transitively; adding a direct dependency would still be an explicit reviewed manifest change. Tempfile's no-clobber API disclaims universal atomicity/security and can leave a temporary link. |

## Review Reconciliation, 2026-10-27

Original Codex plan: fresh-only Plan output plus a nonserialized protected-path
context, preservation of loader-lost sidecar facts, broad native alias admission,
preflight/recheck and staged publication. Claude's supplied review recommends a
smaller single host-path field, derived sidecar exclusion and no-clobber writer.

| Feedback | Revised Codex decision |
| --- | --- |
| Fresh-only already protects all existing files/links | Accept. Refuse any existing destination, including dangling symlinks, before opening it for writing. Final no-clobber publication handles a new competitor. |
| No loader-fact capture is needed | Accept after tracing the save transaction. Protect the derived save-owned sidecar tree; other existing files are protected by fresh-only. Do not change loading, retain discarded references or reread/scan sidecars. Do not claim a subtree exclusion for arbitrary external directories. |
| One serde-skipped `host_project_path` is enough | Accept the field over a context type/registry, but explicitly propagate it in `clone_without_history`. A derived Clone alone does not reach managed jobs. Preserve the current live binding across detached commit. |
| Save As during computation need not trigger live rebinding | Accept the publication boundary: a later explicit project save may replace a report, or no-clobber may refuse it. Do not promise report survival against later explicit saves. Separately prevent a detached commit restoring its old host path over the live binding. |
| Use parent identity and conservative ignore-case leaves | Accept the smaller strategy, but not a universal Unicode claim. Case folding is not Unicode normalization. Native alias gaps must be conservatively refused or explicitly resolved before support is claimed. Only NotFound permits skipping an absent project-parent check; permission/identity errors must not become non-alias. |
| Drop blanket traversal rejection and a separate pre-publish path recheck | Accept under the stated non-hostile-parent boundary. Resolve the existing directory using native path semantics; final no-clobber, not another stat, handles destination races. Unsupported/malformed forms still fail closed. No stronger directory-swap claim. |
| Reuse writer/error patterns | Accept: Plan-specific helper, existing InvalidInput/Io plus stable reason prefixes, exact JSON bytes and original publication errors. No new public error enum variants or receipt schema. |
| Separate no replacement, completeness, cleanup and durability | Retain. No universal atomicity, crash durability or adversarial filesystem sandbox claim. |

## Revised Minimal Contract

Only `PlanDnaSequenceDesign.path` changes. No-path Plan and S1 remain file-free.
Approval of every present path remains a write-attempt gate, including cancelled
reports; raw MCP confirmation still precedes payload/state reads. No new payload
resolver, parser, command grammar, operation field or scientific digest change.

### 1. Fresh-Only Plan Writer

- Parent directory must already exist. Never create it or an implicit filename.
- `symlink_metadata(destination)` success means refuse every file type. Errors
  other than NotFound also refuse; do not convert access failures to absence.
- Stage identical `to_string_pretty` bytes in a uniquely created same-directory
  temp file, using the existing writer pattern. Flush/sync a write-capable handle
  and publish with `persist_noclobber`, retaining its original error.
- No overwrite flag, mkdir, delete-existing, retry-with-replacement or global
  change to `write_pretty_json_file`. Identical prior reports are still existing
  destinations: use a fresh explicit filename and obtain required command review.
- Perform destination/protection admission after approval, in execution/worker
  phases and before the search for fail-fast behavior. No-clobber publication
  itself is the final-name race guard; no separate path recheck is treated as proof.
- Cancelled reports use the same writer. Refusal/write failure creates no success
  journal, leaves project/undo/redo/dirty revisions unchanged and never relabels
  search completion. Published files can survive later workflow/service failure,
  cancellation or stale commit, as in Slice 3.

### 2. One Nonserialized Host Binding

Session engine field: `#[serde(skip)] host_project_path: Option<PathBuf>`.
Bind an absolute spelling of the host's actual project path, even when missing;
do not canonicalize the absent project file or infer it from project metadata.

- CLI/MCP bind the actual state/project argument, including defaults. Python
  needs no new implementation because it delegates to CLI.
- GUI binds its current project on load/Save As and clears it on a new unsaved
  project. File I/O/identity resolution stays outside GUI/admission locks.
- Generic JS/Lua state-in/state-out calls remain unbound. Loading state into a
  script does not serialize a filename into it. Universal no-clobber still
  applies, but no unknown absent backing file is inferred or protected.
- Copy the field in `clone_without_history` and thus `fork_detached_execution`.
  During guarded commit retain the **live** host binding, even if a worker used
  an older path for its preflight. Do not add a new stale predicate/registry or
  artifact-commit phase. Existing owner/history/structural checks stay intact.
- Audit resets/reconstruction, including shell `LoadProject`, and rollback:
  do not lose the host's output binding or accidentally reinterpret an import/
  export as a Save As. Keep each host's existing current-project semantics.
- The field never enters ProjectState, operation payloads, report/approval/command
  hashes, serialized engine data, undo checkpoints or dirty-state accounting.

No loader result/context API, sidecar reference preservation or second file read
is needed. Reuse the existing pure sidecar-name derivation rather than duplicate
its suffix/fallback behavior in adapters.

### 3. Small Protected-Path Check

Resolve the output's existing parent using native filesystem semantics. Compare
its identity to the bound project's existing parent, then conservatively compare
the leaf to the project leaf and derived sidecar name ignoring case. For an
existing derived sidecar root, compare identities while walking the resolved
output parent's ancestors; reject the whole tree, including new unused leaves.
Refuse a file named as the derived sidecar root even when that root is absent.

If the project parent truly does not exist, no report with an existing parent
can be created there, so that comparison can be skipped. Do not skip errors
such as permission denied or unavailable identity. Avoid hand-written textual
normalization, extension rules or treating a failed identity comparison as false.

`same-file = "1.0.6"` is now a direct root dependency, reusing its already locked
package/checksum. Its portable native directory handles avoid bespoke unsafe
Windows identity code. Only GENtle's dependency edge changes in Cargo.lock;
no package version/feature refresh or zero-build-effect claim. Glen/CI must
check headless/desktop targets and the locked manifest after this explicit change.

Ordinary resolvable `..`, `/tmp`/`/var` aliases and junctions need not be banned
merely for being aliases; use the actual resolved directory. Do not claim that
this protects against someone swapping a parent during the operation. Unknown
device/stream/verbatim forms must fail, not silently bypass the existing-parent
or protected-name checks.

Unicode remains a precise support decision, not a reason for a general path
framework. ASCII ignore-case is insufficient for composed/decomposed aliases.
The small safe fallback is to refuse ambiguous same-parent protected-name
comparisons involving non-ASCII names until coverage exists; describe the
over-refusal rather than silently claiming all native normalization is handled.
Claude did not independently review this Codex qualification; the owner has now
authorized the revised plan including conservative refusal. Stronger comparison
and actual native-filesystem acceptance are not implied.

## Implementation And Pending Verification

`src/engine/io/sequence_design_preview.rs` owns admission/publication; its tests
are inline-synthetic fixtures in `sequence_design_preview/tests.rs`. The handler
preflights before search and leaves other JSON writers alone. Raw/named/shared
CLI/MCP bind host paths after existing gates. GUI project binding follows load/
Save As/reset; shell LoadProject and pool import retain the host output binding.
Explicit clones, serialized checkpoint import and detached commit keep the correct
session binding. Generic JS/Lua need no new runtime API; tests exercise their
existing bindings. The existing CLI/MCP/Python companion is extended, not replaced.

Tutorial 06.07 and its source/agent/discovery guidance explain fresh outputs and
support limits. Both strategies retain separate previews and all biological,
digest and undo assertions. Engine/service/GUI/VM/MCP tests are authored, not run
here. Source checks, scoped formatting and documentation checks are not Rust or
native acceptance. Use the [consolidated handoff](glen_next_verification_handoff.md#lifecycle-slice-4-2026-10-27)
for exact external commands and changed files.

Post-rebase local checks pass 22 focused Python source/teaching/maturity checks
with 12 explicit actual-process skips (34 total), plus all 75 LF/CRLF-aware
checkout tests. Canonical maturity freshness required only main's new glossary-
source hash in its four projections, not changed capability definitions or
verdicts. Documentation status, scoped rustfmt and whitespace pass. None of
these checks execute GENtle, Rust tests, VMs or a native GUI.

Use temporary synthetic files and tutorial 06.07's hand-crafted CDS. Preserve
exact before/after bytes, source/report/approval hashes and existing scientific
assertions. Proposed Rust prefix `lifecycle_slice4_`; require nonzero counts.

| Boundary | Required fail-first controls |
| --- | --- |
| Writer | Existing/empty/identical report, directory, dangling symlink, hard link and supported special file; missing parent; no directories created and no destination opened for truncation. |
| Project/sidecar | Absent/default/empty/existing state, relative spelling and parent alias, conservative case variant, derived-root leaf and new file inside that root. Loader warnings stay unchanged; no captured custom-root metadata. |
| Binding lifetime | Actual CLI/MCP defaults/aliases, GUI load/Save As/new project, explicit snapshot copy, reconstruction/rollback, and a managed commit that keeps a newer live binding. JS/Lua stateless controls and serde exclusion. |
| Faults/races | Deterministic staged write/sync/persist failure and a barrier-controlled final-name competitor; preserve competitor bytes, original error and owned-temp cleanup. Do not assert hostile-parent or power-loss guarantees. |
| Workflow/service | Two distinct fresh reports; duplicate path refuses without replacing the first; failure/cancel/stale prefix survives; loaded targets, journals, undo/redo and owner checks unchanged. |
| Adapters/agents | Direct/named/forwarded/shared-shell CLI, raw/shared MCP, workflow, managed execution, actual Python and generic JS/Lua; approval before I/O, ask/nested-agent guards, quoted paths, LF/CRLF input/project bytes. No live model call. |
| Native support | Windows handles/case/verbatim/UNC/reparse behavior, macOS aliases/case/Unicode, Linux identities and symlink/hard-link cases. Unsupported privilege/filesystem cases are explicit unrun evidence, not passes. |

Update protocol/CLI/GUI/decisions and tutorial 06.07 in the same behavior change.
Teach fresh paths and existing parents. Give each strategy/repeated preview its
own output path; do not unlink previously retained evidence or remove either
strategy/scientific assertion to make tests pass. Keep the original JSON bytes.

Glen/CI owns all Rust, actual-process/VM and native execution. Use same-source
normal CI rather than a new audit campaign. The existing handoff records
default `cargo check --locked
--tests`, headless CLI/MCP check, `lifecycle_slice4_`, sequence-design, Slice 3,
command-execution and JS/Lua filters, registry/parity freshness and formatting.
Extend/run the existing actual-process lifecycle and tutorial companions with
explicit same-source `GENTLE_LIFECYCLE_BIN_DIR`/`GENTLE_TUTORIAL_BIN_DIR`; never
allow Python to fall back to a local Cargo build.

Claim only no replacement of existing destinations, complete staged bytes before
final-name publication on supported filesystems, and best-effort cleanup. A
later explicit Save As can overwrite a report; parent swaps, crash durability
and universal atomicity remain unclaimed. No loader/save-transaction redesign,
new context layer, overwrite approval, generic exporter migration, S2, queue or
new release gate is included in this implementation.
