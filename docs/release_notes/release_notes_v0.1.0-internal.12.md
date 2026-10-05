# GENtle `v0.1.0-internal.12` Release Notes

Status: development draft following the owner-confirmed `.11` publication at
`a51bbc064ea29d4e7a03c559c5a6e0a6f8e1bd29` on 2026-09-30.
No `.12` tag, package publication, performance improvement or scientific
acceptance is claimed. The [`.11` notes](release_notes_v0.1.0-internal.11.md)
retain the successful artifact cycle and earlier failed-run history separately.

## First Optimization Step

Native installers and the headless container select **`package-opt1`**, a
dedicated Cargo profile inheriting `dev`. It sets `opt-level=1` and `lto="off"`;
`lto=false` would allow local ThinLTO and is not equivalent. This is a controlled
first step, not a return to the former fully optimized release recipe.

- Keep 256 codegen units, no incremental compilation, no debug information,
  debug assertions and overflow checks enabled, and panic unwinding. Local
  `dev` builds remain unoptimized. Normal GUI/export content remains identical
  across profiles by contract; package smokes still need to verify the candidate.
- Keep one Cargo build job, the existing job allowance and resource diagnostics.
  Native packages contain the same five binaries; the container contains the
  same three headless binaries. JS/Lua remain optional source interfaces.
- Keep `cargo-bundle` and `rnapkin` on their pinned, unoptimized helper recipes.
  Preserve RNAPKIN's corrected lockfile and SVG/PNG smokes.
- Name archives with `-package-opt1`, stage from `target/package-opt1`, and bind
  the explicit optimization and safety settings into build receipts. Reject
  missing/mismatched settings rather than combining different recipes.
- Keep `bench-audit` and `release-fast` unchanged for historical comparisons;
  neither selects the new package build. Higher optimization, thin LTO, reduced
  codegen-unit counts, stripping and panic-abort are not enabled by this change.

## Learning And Agent Workflows

- Tutorial 04.08 now teaches **Find Primer-Pair Markers for PATZ1 Isoforms**:
  compare annotation sources, review the request, and explain why seven pairs
  cover 13 cDNA classes yet leave nine class pairs unresolved in the retained
  example. Planning/readiness/dossier steps are optional continuations. The
  endpoint-only feasibility command is no longer prescribed for this SYBR case.
- Approval-bound outer routing for TSS collections remains on the shared
  command/proposal path. Retained branch-local evidence does not imply live
  inner-agent or scientific acceptance of this new candidate.
- Tutorial review projections now use Git history and repository-relative
  paths, never checkout timestamps. Unavailable history is an explicit check
  error; stored catalog age follows the revision date, while live age warnings
  remain current-date-based. No human review or historical panel is relabelled.

## Acceptance Before Publication

**Validation pending:** historical synthetic transcript-panel reports retain
their true `.11` provenance. That is no longer a rollover error: the revised
contract checks historical integrity separately from a fresh `.12` result
replay, allowing only the two scoped generator fields and their verified ledger
hashes to differ. No stored report bytes or hashes were changed. The
[updated Glen handoff](../internal_12_tutorial_regeneration_glen_prompt.md)
requires the post-integration contract change and records the comparison, path
and checkout checks to complete. The replay evidence below covers only the
reported checks, not package or scientific acceptance.

**Reported replay evidence, 2026-10-02:** Glen reproduced the checkout failure
at `f1c904a4` and repaired the generated 08.03 chapter, tutorial hub and ledger
in `e462de1f`, now integrated. With Rust 1.96.1 he reports successful exact-commit
LF and CRLF replays, each covering parity, 51 examples and 29 chapters, plus
17/17 Python checkout regressions. The historical `.11` panel bytes and hashes
remain unchanged. The two stale human-review warnings remain expected. This
is externally reported evidence, not a local rerun or native macOS/Windows
Actions, `package-opt1` build, or scientific acceptance. The remaining handoff
checks and native-platform results still need to be recorded.

**CI follow-up, 2026-10-05:** the October 4 run `37235516663` at `f523a757`
finished with Linux and Windows catalog-test failures; their preceding drift
checks passed, Windows's full suite was skipped and macOS was not selected.
`19889dae` already corrects that catalog projection. The separate Git-history
guard, checkout mock repair and twentieth TSS case counts are locally validated;
they do not establish current native Actions or package acceptance. The
[candidate-preparation handoff](../internal_12_candidate_preparation.md)
records exact historical evidence, permitted local checks, build-only commands
and open compatibility/runner decisions. No workflow was dispatched.

1. Build the `.12` generator and verify replay against the unchanged historical
   baselines. Regenerate only for a reviewed result change, never merely to
   relabel `.11` provenance. Run version, tutorial, LF/CRLF and packaging-policy
   gates on the committed candidate.
2. Run build-only macOS, Windows, Linux and container jobs at one frozen SHA.
   Retain recipe-bound receipts and extracted-package smokes before considering
   publication. No workflow dispatch, push or tag mutation is implied here.
3. Have Glen compare `dev` and `package-opt1` at the same SHA/toolchain/host,
   recording build time/peak memory/package size plus repeated startup and
   PATZ1 open/resize timings. Keep `bench-audit` evidence separate. Decide on
   `opt-level=2` or thin LTO only after this measured result.
4. Carry forward the [roadmap](../roadmap.md)'s exact-candidate GUI, biological
   and external-tool checks. The existing linear-locus envelope (250 kbp /
   5,000 features) and circular correctness checks are not newly passed gates.
