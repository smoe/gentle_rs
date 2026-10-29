# Lifecycle Slice 6: Degraded Candidate-Store Save Safety

Status, 2026-10-29: owner-authorized implementation and synthetic regressions
authored on main, following source inspection at `186937bd`. Rust/actual-process,
VM and native Windows acceptance remain with Glen/CI. No local Rust build or
push; the previous goal stays paused. This is a separate extension, not an
unfinished promise of Slice 5. The plan was not Claude-reviewed.

## Shared Boundary

Permissive `ProjectState::load_from_path` preserves the exact unresolved
`candidate_sets` reference and adds `candidate_sets_load_warning`, schema
`gentle.candidate_store_load_diagnostic.v1`, status `unresolved`, stable reason
`candidate_store_unavailable`, original source-project spelling and cause.
Index/set and record counts must agree; incomplete records cannot become an
apparently healthy empty store. Strict loading still refuses failed hydration.

An old warning-only state remains unresolved: no path/reference is guessed,
and neither JSON/VM round-tripping nor opening an inline/empty replacement
silently clears its warning. Only reopening an original external reference
with successful complete hydration clears a prior warning. An unhydrated raw
reference also cannot be saved relative to a new destination.
Hydration establishes readability/schema/count consistency, not authenticated
provenance or proof that replacement files are the lost originals.

`ProjectState::save_to_path` checks availability before cloning, staging,
directory creation, replacement or sidecar cleanup. This covers GUI Save/Save As,
CLI/shared-shell/MCP project saves and implicit saves, plus JS/Lua save bindings
and CLI-backed Python. Candidate access also refuses; listing returns a Result
instead of presenting unavailable data as zero rows. The GUI displays that error.
Facts mark `candidate_set.exists` unavailable, so closed-world evaluation returns
unknown; agent context omits its zero count and discloses the missing evidence.
Unrelated sequence/state inspection remains possible. Healthy absent/empty
stores and intentional empty-store cleanup retain their existing behavior.

This protects the save operation, not all preceding work. An autosaving command
may already have changed its transient engine or written a separately approved
external artifact before save refusal. There is no new rollback or checkpoint
policy, parser, approval bypass, payload resolver, queue or exporter migration.

## Recovery And Exclusions

Restore the original sidecar or an intact project/sidecar backup, then reopen
the original project. Restoring a file does not silently authorize an already
degraded in-memory session to save. Keep that session's unsaved work separate;
this slice provides no merge or repair command. Warning-only files whose
reference was already lost require an intact backup, not an inferred filename.

No automatic repair, discard command, recovery overwrite, external retrieval,
crash-durability guarantee or hostile parent-directory race protection is added.
Opaque wrappers, generic journals/concurrency and base/prime-editing S2 remain
separate. Arbitrary trusted code that edits metadata or writes raw JSON itself
is outside the shared persistence contract.

## Verification

Inline, explicitly synthetic Rust fixtures recreate real saver-produced indexes
and JSONL records in temporary folders. Tests cover missing/corrupt/truncated
files, LF/CRLF, quoted paths, exact references/diagnostics, unchanged project and
sidecar trees, new/existing Save As destinations, unknown facts, read/mutation
refusals, legacy warning-only JSON, strict loading and healthy recovery/cleanup.
The companion extends existing opt-in CLI/MCP tests; no binary is auto-built.

Authoring checks: 37 Python/source tests passed, including real disposable
LF/CRLF checkouts. The two new actual-process checks were explicitly skipped
without same-source binaries. Scoped rustfmt, documentation/manual freshness
and whitespace checks passed; none constitutes Rust/runtime or native acceptance.

Follow the [bounded Glen handoff](glen_next_verification_handoff.md#lifecycle-slice-6-2026-10-29).
Require nonzero Rust/process counts and retain native Windows path/handle and
GUI Save/Save As refusal verdicts separately. Python/source checks alone are
not runtime, cross-platform, package or scientific acceptance.
