# Glen: Validate The Internal.13 Development Rollover

Prepared 2026-10-13. This is a handoff, not an executed audit or publication
authorization. The owner confirmed `.12` at `5e6201a6`; its
[published ledger](release_notes/release_notes_v0.1.0-internal.12.md#publication-evidence)
remains separate from this development-version validation.

**2026-10-15 scope update:** the owner selected early-feature/3D review and
minimal expression/storage tutorial paths for `.13`. Use the current
[preparation handoff](internal_13_preparation_handoff.md) and release-note ledger
for ordered candidate gates; the version commands below remain required, not
the complete release scope. Performance/live-agent/historical-science work is
non-blocking carry-forward, not a newly selected release aim.

**2026-10-13 follow-up:** the owner subsequently requests a bounded opt-level
increase. Future packaging now selects `package-opt2`; `package-opt1` remains
unchanged. The rollover-only checks below still validate version/tutorial
metadata; do not describe the combined source as profile-unchanged. Recipe
acceptance and same-SHA opt1/opt2 comparison follow the
[release process](release.md) and
[measurement contract](dna_feature_rendering_latency_plan.md#measurement-profiles),
only with separately authorized execution. No Rust/package build ran on the Mac.

## Bind The Source

Wait for an **owner-provided pushed full SHA** containing the `.13` rollover
and its carry-forward documentation. Local Codex commits have not been pushed;
do not substitute a moving main, `8b6f8186`, `f35a58e6`, `.12` binaries or a
nearby source. Use an isolated clean checkout on your approved external host,
not Steffen's Mac or a private/dirty worktree.

Record repository/ref, full SHA/tree, clean status, raw Cargo.lock SHA-256,
Rust/Cargo versions, OS/architecture, effective profile/environment and resulting
executable identities/hashes. Read AGENTS.md, the
[post-release checklist](release.md#post-release-development-version), the
[`.13` draft](release_notes/release_notes_v0.1.0-internal.13.md) and the historical
baseline contract in [the prior handoff](internal_12_tutorial_regeneration_glen_prompt.md).

## Required Checks: Not Run Locally

The rollover changes root workspace version plus seven local inherited
lockfile package versions to `0.1.0-internal.13`. Verify the independently
versioned `gentle-sequence-design` stays `0.1.0`, with no dependency, profile
or scientific-fixture changes. Use a new explicit local target directory:

```bash
export CARGO_TARGET_DIR=/absolute/local/path/gentle-internal13-validation
export CARGO_BUILD_JOBS=1
export CARGO_INCREMENTAL=0
export CARGO_PROFILE_DEV_DEBUG=0
cargo metadata --locked --offline --no-deps --format-version 1
cargo test --locked -j1 --test release_version_consistency
cargo check -q --locked -j1
cargo build --locked -j1 --bin gentle_examples_docs
BIN="$CARGO_TARGET_DIR/debug/gentle_examples_docs"
python3 -m unittest scripts.test_release_candidate -v
python3 -m unittest scripts.test_tutorial_checkouts -v
"$BIN" --check
"$BIN" tutorial-catalog-check
"$BIN" tutorial-manifest-check
"$BIN" tutorial-check
python3 scripts/check_tutorial_checkouts.py --binary "$BIN" --timeout-seconds 1800
git diff --exit-code -- Cargo.lock docs/tutorial/generated docs/tutorial/catalog.json docs/tutorial/manifest.json
git diff --check
```

The historical panel baseline contract allows only its two scoped generator
fields and validated associated checksum projections to differ during replay.
Original reports retain their true `.11` provenance and raw hashes. Do not
rewrite baseline bytes, broaden exclusions or reuse an old generator just to
pass `.13`. For real result drift, retain the temporary replay and scientific
diff for review rather than regenerating committed evidence automatically.
Git-history-dependent tutorial metadata must be checked on the final committed
source, not an uncommitted worktree or a binary from an earlier source.

## Delivery

Return commands, exit codes, exact source/binary/input/lock hashes, platform,
profile, warnings, failures/skips and whether baseline bytes stayed unchanged.
Rebuild and rerun affected gates after any source repair. Keep an unexecuted
check marked not-run, including unavailable dependencies or unsupported hosts.
The existing two stale human-review warnings are not new scientific approvals.

This request is version/tutorial validation only, not an installer build,
native GUI/benchmark acceptance, new model call or a repeat of Glen's completed
`f35a58e6` learner checks. No tag mutation, dependency update, dispatch or push
is authorized here. Native Windows/macOS, actual packaged-profile performance
and the remaining scientific gates stay distinct from a Linux development build.

## Current CI Evidence

[Normal CI 38432492850](https://github.com/smoe/gentle_rs/actions/runs/38432492850)
completed successfully at `8b6f8186c81636491dadb446e94ca68429ba0860` (final update
2026-10-13 19:05:33 UTC). Windows, macOS, headless CLI/MCP, release safety and
the aggregate selected-platform gate passed. Full Linux and optional native
replay/scroll/audit-CLI jobs skipped. No failure repair is needed for that run.
It does not validate the subsequent local baseline
`024ed6c0520d38c4221517219d44c792084de4e5` or the rollover commits. No new run was
dispatched. Retain the [dated update and prior failure history](glen_audit_followups_20261008.md#repaired-source-audit-authorization-at-8b6f8186).

## Historical Candidate Suggestions: Not The Selected Aim

The October 13 suggestions below are retained for context. The October 15 owner
selection above supersedes them; they are not current blocking alternatives.

1. **Measured startup and interactive responsiveness.** The
   [native baseline](glen_native_gui_audit_20260927.md) exposed the
   37.5--37.8-second help-preparation span; lazy discovery needs S6 cold/warm
   and first-use validation. Use the [latency plan](dna_feature_rendering_latency_plan.md)
   and repeated same-SHA `dev` / `package-opt1` comparison, preserving the
   250 kbp / 5,000-feature envelope, separate stress cases and circular
   correctness. This is a measured acceptance aim, not a promised speedup.
2. **Complete a bounded real-reference primer evidence loop.** Post-tag
   matrix inspection and opt-in phase diagnostics exist, but real-reference
   imports still need [bounded CPU/phase measurements](glen_post_blast_import_phase_trace.md).
   Identify the bottleneck before a repair; preserve biological rules and
   reference binding. Conditional panel fan-out is a separate approval, not
   the proposed starting task; follow the [specialist plan](primer_design_specialist_plan.md).
3. **Consolidate native and scientific acceptance.** The
   [final-candidate handoff](glen_final_candidate_execution_handoff.md#remaining-native-and-human-checks)
   separates completed source-bound learner checks from outstanding native
   readability/orientation, TSS/assay, real microarray/gel and live-agent
   evidence. Preserve FLNA's completed human exon choice and version-label
   addendum, keep its draft inactive and retain `.10`/`.11` ledgers. This aim
   records evidence at a chosen source, not automatic scientific approval.

The previous `.12` responsiveness aim is historical, not automatically `.13`'s
aim. Preserve the detailed [headless baseline](glen_gui_latency_audit_20260927.md),
[fixed comparison](glen_gui_latency_handoff_20260924.md) and latency plan rather
than treating their older SHA-bound passes as current acceptance. Selected-assay
cDNA overview/detail still needs merged-source checking; Glen's `18ae846e`
visual review is separate, and genomic projection needs a bound cDNA/exon map.
No broad GUI extraction, higher optimization or new biological contract is
authorized by any of these choices.

## Owner Decisions

- The owner selected A/B/C; complete the current preparation/readiness checks.
- Integrate and push the separate local `pin-ubuntu-2404` commit before
  the [October 19 rolling migration](release.md#ubuntu-runner-decision). It
  pins only the ten currently rolling hosts and updates that inventory; it
  does not qualify Ubuntu 26 or dispatch CI.
- Authorize any push of the local release and/or pin commits. Nothing is
  pushed and no tag, release asset or publication is changed in this session.
- Supply Glen a pushed exact SHA and authorize the external rollover checks
  above. A successful older CI run is not a substitute for those checks.

## Local Non-Rust Verification

At the documentation handoff, these checks pass:

- `env PYTHONDONTWRITEBYTECODE=1 python3 -m unittest scripts.test_release_candidate -v`:
  41 tests, including unchanged historical tutorial provenance checks.
- `git diff --check`: no whitespace errors.
- `env PYTHONDONTWRITEBYTECODE=1 python3 scripts/maintenance_chore.py release-gate`:
  four checks pass, no warnings or failures.
- `env PYTHONDONTWRITEBYTECODE=1 python3 scripts/maintenance_chore.py session-close --plan docs/roadmap.md`:
  no failures; manual plan-fidelity review remains required, and the pre-commit
  dirty-tree warning names only the six intended carry-forward documents.
- A read-only Python `tomllib` comparison against `024ed6c0` confirms exactly
  seven inherited version changes, unchanged independent version/dependencies/
  profiles and a 400-line roadmap. No generated tutorial/scientific artifact
  changes are included.

These are static/policy checks, not Rust execution or fresh tutorial replay.
The separate runner-pin proposal requires its own workflow/policy checks.
