# Claude Review Prompt: Lifecycle Slice 4

Prepared 2026-10-27 for Steffen to forward manually. **Read-only review, not
implementation. No direct Claude invocation or verdict is claimed.**

---

Please review the scoped Plan preview destination/writer proposal in
`docs/operation_lifecycle_slice4_plan.md`. Give concise critique and a revised
minimal plan. Challenge unsafe guarantees and unnecessary machinery; do not
redesign GENtle's entire exporter, loader, lifecycle or job system.

## Repository, Task And Restrictions

Repo `/Users/u005069/.codex/worktrees/faae/gentle_rs`, branch `gentle_rs_3`,
inspected HEAD `f0b8371f4a99a455c60fc71a7769f3ea0380e7f3`. Local main observed at
`a2637079f5a2ee1b49310f988243497562a65ccc`. No fetch/rebase was performed during
inspection; these source observations predate the owner's commit/rebase request.
Check current drift rather than treating those SHAs as a frozen-HEAD requirement.
Main's three newer commits include the observed-transcript kernel and repair
of Slice 3 test calls to `undo_last_operation`/`redo_last_operation`; do not
mistake the older branch's test API gap for a new writer finding. Core writer/
handler/lifecycle/CLI/MCP files inspected have no HEAD-to-main diff. Verify
current drift with your available read-only tools; report anything unverifiable.

Use Read/Grep/Glob only. No shell/build/test/model/network calls, file writes,
checkout/rebase/commit/push, CI dispatch, credential/private data inspection,
`outputs/`, `data/` or `paper/`. Read AGENTS.md, relevant architecture/decisions/
roadmap, protocol/CLI/GUI and testing cross-platform rules. Keep the earlier
goal paused. Return feedback in your answer, not another file. A positive
review does not authorize implementation. The owner reserves Rust/native
verification to Glen/CI and forbids a new payload-resolution layer.

## Context Already Delivered

Lifecycle Slice 3 (`f0b8371f`) preserves Plan/S1 construction, project content,
undo/redo and raw CLI/MCP implicit-save behavior, with existing journal-
publishing managed commits. Empty/mixed/opaque wrappers remain legacy.
Report-write approval is already required for every present Plan path,
including cancelled reports and invalid destinations; raw MCP confirmation
is additive and precedes loading. Do not re-plan those as missing features.

The unresolved problem is explicit writing: Plan still calls the generic JSON
writer, which creates parents then uses `fs::write`. Suppressing autosave does
not prevent truncating the project through a chosen report path. No-path Plan
and S1 are not writers. Approval of a report write is not overwrite consent.

## Original Codex Plan To Critique

1. Migrate **only** `PlanDnaSequenceDesign.path` to fresh-file-only output in
   an existing parent, with no overwrite flag, auto-name, mkdir, delete/retry or
   broad generic-writer change. Refuse all existing leaves, including identical
   prior reports, hard links and dangling symlinks. Keep exact JSON/digests and
   scientific/search semantics; requested cancelled reports still use this path.
2. Add the smallest engine-owned, nonserialized host I/O binding needed to
   exclude an active project filename **even when absent**, plus known/default
   candidate sidecar roots. CLI/MCP bind their real/default state path; GUI its
   current file or explicitly unsaved state; pure JS/Lua state calls are stateless;
   Python delegates. Preserve loader-lost reference facts during existing loading
   rather than rereading/scanning or changing ProjectState. No operation-payload
   envelope, parser change, early `@file` read or scientific schema/metadata field.
   Exact API is not fixed: ask whether this is too much and identify a smaller seam.
3. Validate raw paths before normalization; compare existing parent identity and
   native leaf/subtree aliases for protected absent paths. Do not canonicalize
   away traversal, globally lowercase, infer from extensions, or assume
   `is_symlink` detects Windows reparse points. Distinguish normal macOS `/tmp`/
   `/var` aliases from adversarial parent replacement; reject unsupported/unknown
   semantics rather than claim protection. Review Save As and snapshot context
   lifetime without silently adding a registry or second artifact-commit phase.
4. Preflight after approval in the worker, compute, then recheck and stage exact
   pretty JSON in a unique same-directory temp file; flush/sync a write-capable
   handle and use no-clobber publication. Reject a final-name competitor safely.
   Preserve original errors; cleanup only owned temps. Never fall back to a
   truncating writer. Separate no replacement, staged completeness, cleanup,
   crash durability and hostile parent-directory races in claims/tests.
5. Keep managed owner/history/stale/cancel behavior and workflow prefixes:
   published reports can survive a later failure/cancel/stale commit. No report
   status relabelling or external rollback. Bind host paths without changing
   command/approval hashes, no-project observation lanes or agent guards.
6. Add fail-first synthetic state/file/race/fault and actual adapter tests;
   tutorial 06.07 teaches fresh destinations. Existing two-strategy tests must
   choose separate filenames rather than overwrite one report. Retain exact
   reports, biological assertions, approvals and both strategy runs. Document
   behavior/limits in existing contracts/decisions/changelog, with verification
   pending on the roadmap. Leave Rust/native execution to Glen/CI.

## Evidence To Inspect

- `src/engine/ops/operation_handlers.rs::apply_internal`, Plan arm.
- `src/engine.rs::write_pretty_json_file`, `write_text_file_atomically`,
  `ProjectState::{load_from_path,save_to_path,candidate_store_sidecar_dir}` and
  candidate reference hydration/transaction helpers.
- `src/engine/ops/candidate_guides.rs::ensure_output_parent_dir`.
- `src/engine/io/gel_images.rs`, `src/engine/io/gel_workbench.rs`: existing
  `NamedTempFile` + `sync_all` + `persist_noclobber` export precedents.
- `src/tss_profile_export.rs::{checked_directory,destination,is_parent_component}`:
  raw Windows traversal precedent, not a universal single-file/TOCTOU guarantee.
- `src/app.rs::current_project_path` and Save As/load/reset handling;
  `src/app/routine_and_agent_assistant_ui.rs` managed command submission.
- `src/engine_shell.rs::{ShellExecutionOptions,resolve_operation_payload}`,
  `src/command_execution.rs`; raw/named/forwarded CLI and MCP execution helpers;
  generic JS/Lua and CLI-backed Python routes. Identify a binding bypass.
- `src/engine/analysis/sequence_design/tests.rs` tutorial test reuses one report
  path across strategies; Slice 3 engine/service tests and
  `scripts/test_operation_lifecycle_slice3.py` also contain repeated targets.
- Tutorial `docs/tutorial/06-07_synthetic_sequence_design.md`, input provenance
  README and `scripts/test_sequence_design_tutorial.py`.
- `Cargo.toml`/`Cargo.lock` resolve `tempfile 3.27.0`. If accessible, inspect its
  cached `src/file/mod.rs::persist_noclobber` and `file/imp/{unix,windows}.rs`.
  Its docs explicitly disclaim universal atomicity/security; Unix can use
  link/unlink fallback. Do not assume library naming proves stronger guarantees.

## Questions And Required Reply

1. Is the fresh-file-only default acceptable and minimal? Name all compatibility
   changes rather than reintroducing overwrite under a hidden exception.
2. Can one smaller existing seam supply absent-project and sidecar protection?
   Does proposed context preservation through loading/fork/Save As create hidden
   project mutations, lost host binding or expensive admission work?
3. Are the alias/native-filesystem and race claims defensible? Distinguish a
   real safety requirement from a stronger deferred sandbox guarantee. Propose
   fail-closed support boundaries or a smaller safe cut, not hand-wavy canonicalization.
4. Is cached tempfile behavior sufficient with the stated limits? Which tests
   must fail under today's writer and exercise the actual adapter/file boundary?
5. What belongs in the first implementation and what must remain separately
   approved? Avoid new dependencies/unsafe OS layers without evidence.

Return source baseline and drift, concise severity-ordered findings with file/
symbol evidence, original-plan points accepted/rejected, and a revised minimal
implementation order plus Glen/CI acceptance matrix. Mark unresolved owner
decisions and unverified platform assumptions. Claim no execution or acceptance.
