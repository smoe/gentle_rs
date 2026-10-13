# GENtle `v0.1.0-internal.13` Release Notes

Status: development draft after the owner-published
[`.12` prerelease](release_notes_v0.1.0-internal.12.md#publication-evidence)
at `5e6201a69a599cba87491c86cfd7311cb6fee173` on 2026-10-13.
**Release aim: undecided.** No `.13` candidate SHA, tag, packages, publication,
performance improvement or scientific acceptance is selected or claimed.

## Post-Release Development

- The version rollover itself changes only application/workspace metadata and release
  projections. Dependencies, `package-opt1`, `bench-audit`, `release-fast` and
  the independent `gentle-sequence-design` version remain unchanged.
- Main already contains post-tag read-only primer-panel specificity-matrix
  inspection (`dee6033d`, merged through `a76b8418`) and opt-in import-phase
  diagnostics (`52908b2f`), with subsequent fixture/identity repairs. These are
  `.13` development content, not features of the released `.12` packages.
- Diagnostics are not a demonstrated performance repair. Current-source CI,
  real-reference measurements and native acceptance retain distinct scopes.

## Next Bounded Optimization Step

The owner's 2026-10-13 request selects `package-opt2` for future native/container
builds, raising only opt-level 1 to 2. `package-opt1` remains unchanged for
comparisons. Keep LTO off, 256 codegen units, assertions/overflow checks on,
unwinding, no stripping, incremental/debug off, one build job and unoptimized
helpers. Build paths, archive suffixes and receipts name `package-opt2`.
No dependency or application behavior changes. Exact-source package builds,
extracted smokes and measured resource/runtime comparison remain not-run;
no speedup, native acceptance, push or dispatch follows from implementation.

## Carried-Forward Gates

1. Validate S6 cold/warm startup and first use after lazy help discovery using
   the [latency plan](../dna_feature_rendering_latency_plan.md). CPU checkpoints
   and headless frame tests do not measure native presentation latency.
2. Have Glen compare `dev`, `package-opt1` and `package-opt2` at one SHA/toolchain/host: build
   resources, repeated startup, PATZ1 open/resize and unchanged UI/export output.
   Successful `.12` hosted packaging does not establish a speedup.
3. Resolve the [real-reference primer-import performance block](../glen_post_blast_import_phase_trace.md)
   from bounded phase/CPU evidence, preserving terminal rules, HSP coverage,
   reference bindings, fingerprints and atomicity. Do not infer biological
   failure from a timeout or authorize panel fan-out before this prerequisite.
4. Complete the [remaining native/human tutorial checks](../glen_final_candidate_execution_handoff.md#remaining-native-and-human-checks),
   independent extracted-package review, 08.04 host-policy prerequisites and
   live-agent work only at an explicitly authorized source and host.
5. Keep FLNA native orientation/readability, historical-release and experimental
   checks pending. The human teaching-exon choice and transcript-version label
   addendum are already recorded; do not reopen them or activate the draft.
6. Preserve the separate `.10`/`.11` scientific ledgers. New checks at a new
   SHA do not backfill acceptance of those released artifacts.

## Verification Boundary

Normal [CI run 38432492850](https://github.com/smoe/gentle_rs/actions/runs/38432492850)
finished successfully at `8b6f8186c81636491dadb446e94ca68429ba0860`, final update
2026-10-13 19:05:33 UTC. Windows, macOS, headless CLI/MCP, release safety and
the selected-platform gate passed; full Linux and optional audit jobs skipped.
This resolves that repaired-fixture run, not later main `024ed6c0` or this
rollover. See the [dated audit update](../glen_audit_followups_20261008.md#repaired-source-audit-authorization-at-8b6f8186).

Local Rust execution is prohibited for this session. `cargo metadata`, the
Rust version guard, Cargo check/build/test and a freshly built tutorial replay
are **not run**. The [rollover handoff](../internal_13_rollover_glen_handoff.md)
names the required external commands and receipt boundaries. Historical `.11`
tutorial bytes/hashes remain intact; do not regenerate merely to print `.13`.

All existing tags/assets stay unchanged. No push, workflow dispatch, benchmark
run or publication follows from this draft.
