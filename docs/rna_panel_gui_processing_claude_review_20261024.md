# Claude Review: GUI-Only FASTQ Processing And Monitoring

Status, 2026-10-24: review delivered and reconciled with the owner. The authorized
implementation ledger is [RNA panel processing and monitoring](rna_panel_gui_processing_plan.md).
The original review brief below is retained as context, not a second active plan.

## Task And Boundaries

Review this plan read-only and return a concise critique plus a revised minimal
implementation plan. Do not implement it or edit this document.

Repository: `/Users/u005069/.codex/worktrees/faae/gentle_rs`.
Verified starting revision: `7c87bf4378677283fa65e1c81167bdd988e24dc5`
(`feat(rna): prepare gene-set panels with shared graph provenance`). Recheck
HEAD and relevant files first; distinguish already-landed work from remaining
gaps. At preparation, local `main` and `gentle_rs_3` point here; this is not a
claim that remote `main` contains it. Unrelated untracked `outputs/` is present.

Keep the owner's previous goal paused. Do not build/test Rust, run GENtle/native
GUI, invoke models/providers, access private sequencing data, inspect `outputs/`,
capture screenshots, regenerate artifacts, commit, rebase or push. Rust/native
execution and acceptance belong to Glen/CI. Return findings in your reply only;
identify static evidence versus proposed tests and unverified hypotheses.

Read `AGENTS.md`, relevant sections of `docs/architecture.md`,
`docs/decisions.md`, `docs/roadmap.md`, `docs/gui.md`, `docs/cli.md`,
`docs/protocol.md`, and `docs/testing.md` (cross-platform rules).
Consult `docs/asynchronous_command_execution_plan.md`, especially the existing
managed-command tranche, rather than treating its historical proposal as open
implementation work. This is a bounded follow-up, not a new release gate.

## User Outcome

A wet-lab user should be able to prepare a reference panel, select completed
FASTQ.gz chunks, review sample/checksum/policy bindings, start a screen, continue
using other windows and the Agent Assistant, monitor/cancel it, inspect retained
evidence, and export summaries without typing shell commands or hand-editing
JSON. CLI, MCP, scripting and the inner agent must retain equivalent scientific
operations, observation and control. "GUI-only" describes user convenience,
not GUI-owned biology or a requirement that agents drive pixels.

Reuse GENtle's existing gene sets. The current requested panel has ten names:
TP53, TP63, TP73, PATZ1, POU2F1, POU2F2, FUS, TDP43, MDM2 and MDM4, with growth
to about twenty enriched targets. Resolve aliases against the explicit source;
do not silently infer the stable identity of TDP43 from its label. The initial
libraries are owner-bound `26R266-1 / barcode13` and `26R266-2 / barcode14`, both
custom-primer/adapter-only; filename numbering is not an identity oracle.
Review synthetic examples only, not those actual downloads.

## Original Codex Plan To Review

1. **Native screen setup.** Extend the existing RNA-panel workspace with saved
   panel selection, explicit completed-file selection, supplier checksum-manifest
   input, per-file sample/barcode/priming bindings, policy controls and fresh
   output/control paths. Provide request review/copy/save and approved Run.
   Parse gzip directly; no upfront full checksum pass, full decompression/copy,
   growing-directory watcher or implicit resource preparation.
2. **Nonmodal run monitor.** Submit through existing command admission/runtime
   machinery, not another scheduler. Show job identity, phase, file ordinal,
   reads and compressed-byte progress, with cooperative cancellation. Other
   windows and conversation remain usable; closing the monitor is not cancel.
   Admission/capacity failures are explicit, not a silently invented queue.
3. **Saved-result inspector.** Project shared panel/screen/summary records into
   sample coverage, verified/missing inputs, transcript support and bounded
   nomination inspection. Preserve ambiguity and all outcome distinctions.
   Sorting/filtering rows is presentation, never reassignment of reads.
4. **Summary export and graph.** Invoke the existing summary exporter, reopen
   outputs and explicitly link panel -> screen -> summary in project/SVG graphs.
   Keep portable immutable evidence separate from transient runtime activity;
   graph registration must have a shared engine contract and clear undo/file
   semantics, not GUI-only metadata writes.
5. **End-to-end tutorial acceptance.** Extend tutorial 10.04 through these
   native controls, including cancellation, partial/control-error inspection,
   explicit verification and summary export. Pair each GUI action with its
   CLI/MCP/inner-agent route. Glen captures source-bound synthetic screenshots
   and verifies native/hosted interaction; authored tests alone do not prove it.

This proposes reusing shared contracts before adding new ones. A bounded queue
may be a later shared-service extension if justified, but is not already shipped
and is not permission to introduce a general workflow/scheduler abstraction.

## Evidence Already Gathered

- `src/app/rna_panel_ui.rs`: the new native preparation editor handles refresh,
  gene-set binding, opening and building panels; it does not yet have a screen
  form or screen/summary inspector. It has snapshot-bound background work,
  shared operations, UI intent and graph reopening. Review which pieces should
  be reused rather than duplicating them.
- `src/engine/protocol/rna_panel.rs`: portable build, screen, summary, progress,
  file-receipt and outcome records. `src/engine/analysis/rna_panel.rs` and its
  `input.rs`, `control.rs`, `summary.rs`, `verification.rs`, `tools.rs` and
  `project.rs` children own biology, streaming intake, retained units, exports
  and optional panel registration. No new classifier/aligner is proposed.
- Existing shared routes are `rna-reads panel-demo`, `panel-build`,
  `panel-screen`, `panel-show`, `panel-cancel`, `panel-export-summary`, and
  `panel-build --list-sources [GENE_SET_ID]`. See
  `src/engine_shell/command_parsers.rs`, `src/engine_shell.rs`,
  `src/engine/ops/operation_handlers.rs` and `docs/protocol.md`.
- `src/command_execution.rs`: existing process-local `CommandExecutionService`
  with submit/status/cancel/take_result, detached guarded publication, four
  active workers, 32 retained receipts and a 16 MiB interactive-output limit.
  It handles typed RNA-panel progress. Result consumption is one-time; receipts
  may be evicted. These handles are not a durable queue or cross-process server.
- `src/runtime_status.rs`, `src/background_engine.rs`, and
  `src/app/routine_and_agent_assistant_ui.rs`: runtime observations, owner/
  revision checks, managed command admission and existing UI controls. Runtime
  observation is not the scientific completion receipt. Registry snapshots use
  synchronization; do not describe them as lock-free.
- `src/agent_feedback.rs` distinguishes screen cancellation/control errors from
  successful inspection of those reports. `src/agent_bridge.rs`,
  `src/agent_execution.rs`, `src/mcp_server.rs` and
  `crates/gentle-shell/src/ui_intent.rs` cover agent approval and adapter access.
- Tutorial source is `docs/tutorial/10-04_targeted_nanopore_fastq_panel.md`, with
  source/catalog metadata under `docs/tutorial/sources/` and `catalog.json`.
  Its synthetic demo already exercises two small libraries, shared-exon
  ambiguity, a distinguishing junction that does not eliminate another
  candidate, same-pass checksums, inspection and summary export. New preparation
  acceptance is listed in `docs/glen_next_verification_handoff.md`.

## Questions Requiring A Decision

1. **Actual parity and ownership.** What can use existing routes unchanged,
   and what minimal shared typed contracts are missing for preflight,
   start/observation/control, saved-run lookup or graph registration? Supply a
   small GUI / CLI / MCP / inner-agent matrix, including confirmation/read-only
   behavior. A Rust host API alone is not proof of a usable external-agent route.
   Explain same-host versus separate-process access: process-local numeric job
   IDs cannot silently become portable IDs or evidence of restart recovery.
2. **Cancellation and external artifacts.** Trace managed-command cancellation
   through the screen's cooperative callback/marker, shell result, agent
   feedback and result consumption. The command service can return an error at
   its final cancellation/staleness gate after external files were written.
   Can a valid cancelled/control-error report become inaccessible in the GUI?
   Recommend the smallest fix, if needed, that exposes verified retained units
   and bound paths without publishing stale engine mutations or success claims.
   Do not infer a confirmed defect from one wrapper without tracing callers.
3. **Three separate statuses.** Define command execution, scientific run
   outcome and project attachment separately. A completed command or successful
   `panel-show` does not turn `control_error` into successful analysis. A marker
   failure is not authorized cancellation, even when every file verified.
   Inspection after a project switch may remain possible as historical evidence
   without being attached to the new project. How should controls show this?
4. **Intake and resource boundaries.** What shared preflight can run without
   reading terabytes? Preserve explicit file-completion attestation, ambiguous
   manifest-name refusals, original/effective request binding, no-overwrite
   output rules, path/handle safety and same-pass raw compressed-byte hashing.
   Bytes/read-ahead and read counters must not imply verified EOF or safe ETA.
   Do not load large sidecars into RAM or exceed interactive-result limits;
   default to metadata and explicitly bounded rows/verification.
5. **Lifetime and recovery.** What happens when a view closes, a project is
   replaced/edited, a worker fails, capacity is reached, a result is consumed or
   evicted, the process restarts, or the external drive disappears? Reopening a
   published report is not worker reattachment. Lost activity stays unknown/
   unassessed until supported by evidence; no automatic retry, rescan, in-file
   resume, changed-input resume or invented successful completion. Explain the
   boundary between an explicit remaining-file run and unsupported continuation.
6. **Graph and evidence durability.** Reuse panel registration patterns but
   justify the exact stored receipt/reference boundary for screens/summaries.
   Avoid embedding huge nominal-read tables in the project or refreshing hashes
   from external disks during paint. Retain sample/panel/policy/hash provenance,
   partial diagnostics, deliberate relocation checks and undo semantics. New
   panel versions must not silently inherit old screening counts.
7. **Minimal UI and delivery order.** Can setup/monitor/results/export be tabs
   or views of the existing workspace instead of independent overlapping tools?
   Recommend the smallest complete first tranche and ordered follow-ups, not
   just another setup form that still requires shell/JSON for the rest. Identify
   any correctness prerequisite before adding controls. Keep queueing bounded
   and separate unless the end-to-end requirement genuinely needs it now.

## Scientific And Approval Invariants

Preserve class-bound partial co-nomination sets and overlapping per-member
support; no forced winner, fractional allocation, whole-read compatibility,
exclusive isoform call, expression or molecule-copy claim. Unavailable evidence
is not numerical zero. One candidate's junction support does not erase another
candidate's shared-exon nomination. Keep legacy class association unresolved.

Summary export must retain its same-panel/full-policy and duplicate-input
guards. Different file hashes do not prove independent reads: a facility's
merged FASTQ and constituent chunks can overlap. Display filters must not alter
the competitor universe; source/transcript changes require a new panel/run.

Existing artifact writes, project replacements and registration require explicit
approval. New paths cannot bypass parser-classified mutation guards, transfer
approval after request edits, or permit recursive agent invocation. Monitoring
must remain read-only and must not journal polling that invalidates active work;
control requests retain their existing distinct approval and lifecycle rules.

## Required Review Output

- Verified HEAD, inspected paths and what is already implemented.
- Severity-ordered findings with concrete file/line evidence; distinguish
  demonstrated defects from static risks and owner decisions.
- Revised minimal plan with shared types/routes, adapter ownership, mutation
  rules, non-goals and focused commit boundaries. Name anything unnecessary in
  the original plan, not just extra features to add.
- A verification matrix for Glen/CI: synthetic success, ambiguity, no matching
  reads versus missing evidence, cancellation before/during/after verified
  units, control errors, hard checksum/budget failures, stale owners, output
  conflicts, result consumption/eviction, save/restart, LF/CRLF and quoted
  cross-platform paths. Compare canonical scientific payloads across adapters,
  excluding documented timestamps/operation IDs rather than weakening hashes.
- Tutorial 10.04 changes needed for a complete no-shell rehearsal, with semantic
  GUI targets, typed checkpoints and later source-bound screenshots. Include
  deterministic inner-agent route/approval tests; live model usability and
  native behavior remain separate Glen/CI checks, not claims of this review.

Do not expand into POD5, alignment-based isoform resolution, biological
calibration, throughput claims, generic mapping primitives or a second planning
system. Conclude with recommended next implementation slice and unresolved
questions only. No implementation or execution is authorized by this prompt.
