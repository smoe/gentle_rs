# Lifecycle Slice 6: Degraded Candidate-Store Save Safety

Status, 2026-10-29: owner authorized the Claude-revised implementation, including
the recommended unknown-schema and legacy-warning decisions. Runtime changes
and synthetic regressions are authored; Rust, actual-process/VM and native
verification remain with Glen/CI. No push or resumption of the previous goal.

## Review And Reconciliation

The [original read-only prompt](operation_lifecycle_sidecar_claude_review_prompt.md)
proposed retaining unresolved references, guarding saves and candidate use, and
keeping unaffected inspection available. The owner supplied Claude's review of
`186937bd`, then requested implementation. No further Claude invocation was made.

Claude confirmed the destructive derived-sidecar removal path and required
guards at candidate read/write boundaries as well as the saver. It also found
silent GUI save errors, legacy warnings cleared on reload, and Save As resolving
retained relative references against another destination. Accept these repairs.

Reuse the overlapping main-checkout work through a scoped source snapshot,
without editing that checkout. Remove its extra post-load strict check, count/
duplicate validation and diagnostic source/cause duplication. The existing
reference/index hydration grammar and count defaults remain unchanged. This
worktree starts from `5411ad10`; later-source integration needs its own checks.
During commit/rebase the overlapping foundation landed as `1fdeb1d9`. Retain
its historical authoring evidence separately; the reconciled runtime/scripts
are unchanged from this worktree's reviewed pre-rebase implementation.

## Approved Minimal Contract

1. Failed permissive hydration keeps the exact raw `candidate_sets` value and
   adds `candidate_sets_load_warning`, schema
   `gentle.candidate_store_load_diagnostic.v1`, status `unresolved`, stable code
   `candidate_store_unavailable`, message and existing strict-env hint. Only
   successful existing hydration of an original reference clears a prior warning.
2. One in-memory `ensure_candidate_store_available` check rejects a warning,
   unhydrated reference, malformed inline data or unsupported candidate schema.
   It never retries I/O or resolves a reference relative to a save destination.
3. Shared readers/writers/listing use that check. Listing returns a Result,
   shell/routine/GUI consumers report unavailable, and project fact graphs add
   `unavailable_facts` only when nonempty. Candidate absence is unknown, not false;
   deterministic inner-agent context omits the zero count and explains the gap.
4. `save_to_path` checks before cloning, staging, mkdir, replacement or cleanup.
   This protects Save, autosave, Save As and existing project-copy routes through
   shared persistence. Healthy absent/inline/empty-store cleanup is unchanged.
5. GUI save failures are visible without clearing dirty state, rebinding a path
   or continuing a pending project action. CLI/MCP one-shot failures explain
   discarded transient state and non-rollback of separately published files.

Decisions: unknown schemas keep existing load compatibility (including strict
mode), but candidate use and saves fail closed. Legacy warning-only states also
remain loadable but cannot be saved; only a pre-loss backup can restore a lost
reference. General pre-mutation refusal remains deferred. Strict mode still
returns the original error on actual hydration failure; its default is unchanged.

## Recovery And Non-Goals

Restore the source sidecar or an intact project/sidecar backup, then reopen the
original project. The existing degraded session is not silently unlocked. Keep
unsaved GUI work separate; no new recovery merge, discard, force-save or repair
command is provided. Diagnostics are not authenticated; intentionally removing
both the reference and warning in raw JSON is outside this persistence contract.
Successful hydration establishes readability under the existing parser, not
authenticated provenance or proof that replacement files are the lost originals.

No automatic repair, new payload resolver, journal/queue migration, wider
exporter rewrite, count/duplicate validation, Windows-separator migration,
staging cleanup refactor, crash-durability or hostile-parent guarantee. Parseable
truncation is not newly detected. File publication is not reversed by a later
project-save refusal. Prior lifecycle/S2/scientific/release gates stay separate.

## Authored Verification And Remaining Checks

Inline synthetic fixtures recreate saver-produced indexes/JSONL in temporary
folders. They cover actual missing/malformed failures, unknown schemas, warning-
only state, same-path/Save As refusal with intact siblings and a foreign decoy,
LF/CRLF, restoration/reopen, intentional deletion as a negative control, history/
snapshot/serde transfer and unavailable facts. VM tests load warning-only files
through public loaders and refuse both same-path saves and Save As. The existing
process companion adds CLI/shared-shell/MCP/Python refusal and report-survival
cases without auto-building binaries. GUI helper tests do not establish native
dialog or platform acceptance. All new Rust/process/VM checks are authored,
not executed here; source/Python checks alone do not prove runtime behavior.

Use the [bounded Glen handoff](glen_next_verification_handoff.md#lifecycle-slice-6-2026-10-29)
for exact filters, source/binary/platform receipts and remaining native checks.
Preserve the separately reported Slice 3/5 results and rehearsal; do not create
another old-source audit campaign or claim a new scientific verdict.

Local authoring checks: 19 source/teaching checks pass; 15 optional actual-binary
cases skip explicitly. Two real LF/CRLF checkout tests, documentation status,
manual freshness and whitespace pass. Formatter output is applied only to new/
changed Rust lines, retaining unrelated baseline formatting; this is not Cargo
type checking. Session-close: 4 OK, 2 expected working-tree/manual-plan warnings,
0 failures. The roadmap remains 400 lines. Runtime/VM/native tests and Cargo
checks stay unrun; all working changes are scoped to this slice and unrelated
`outputs/` is preserved.
