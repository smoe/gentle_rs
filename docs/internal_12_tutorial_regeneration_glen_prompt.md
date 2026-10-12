# Glen Prompt: Verify Internal.12 Against Historical Tutorial Baselines

Updated 2026-10-02. This replaces the previous mandatory-regeneration request.
The filename remains stable for existing handoff links. This is a prompt, not
an executed job or an acceptance record.

## Request And Source

Please verify the revised tutorial-baseline contract in `smoe/gentle_rs` on
your Linux machine, in a separate bot branch/worktree. Do not build on
Steffen's Mac. No installer build, benchmarks or unrelated tutorial repairs.

The source must include both:

- Integration `a2a88995616e64795a7ec024951d97617f5e6c1b`, which merged your
  `ae21fc0542a030132658a828ad43ffcac0006839` series onto the `.12` rollover.
- The subsequent historical-baseline contract change described here, including
  `tutorial_panel_comparison_value`, `tutorial_generation_ledger_equal` and the
  renamed historical-provenance Python/Rust guards. The integration SHA alone
  is no longer sufficient.

Codex has not pushed this work. Obtain the exact owner-provided/published ref
including the follow-up, verify the integration is its ancestor, and record
the full actual source SHA. If unavailable, report that prerequisite; do not
substitute stale upstream main, your old review branch or an old executable.

Read `AGENTS.md`, the historical tutorial baseline rule in
`docs/architecture.md`, `docs/testing.md`, the post-release checklist in
`docs/release.md`, the roadmap and the `.12` release-note draft.

## What Changed

These reports are historical synthetic references, not real PATZ1 acceptance:

```text
docs/tutorial/generated/artifacts/patz1_transcript_assay_panels_cli/artifacts/patz1_routine_common_region_screen.report.json
docs/tutorial/generated/artifacts/patz1_transcript_assay_panels_cli/artifacts/patz1_sybr_juc_panel.report.json
docs/tutorial/generated/artifacts/patz1_transcript_assay_panels_cli/artifacts/patz1_endpoint_end_matrix.report.json
```

Their true `.11` generator identities and original raw hashes remain intact.
The previous current-version assertion conflated historical provenance with
current compatibility. The new contract separates them:

- The fast gate checks historical schemas, consistent nonempty generator
  identities and the original raw hashes; it does not demand today's version.
- Fresh replay output must identify the currently executing package and source
  revision. Only `selection_audit_generator_revision` and
  `provenance.gentle_version` inside assay `primer_pair_summary` records in
  `selected_assays` and `short_sybr_junction_assays` may differ for these three
  allowlisted files. Everything else must match, including scientific data,
  input hashes, schemas, settings, warnings and selection decisions.
- Before projecting the three associated checksum entries for ledger
  comparison, validate each original checksum against its own file bytes.
  Other ledger fields and checksums are not excluded.
- Fresh exports preserve their real source revision. Comparison is in-memory;
  it must not rewrite either baseline or fresh output. No general recursive
  ignore list, provenance relabeling or hash normalization is authorized.

## Build And Check

Use a clean worktree and an explicit local target directory, not network
storage. Record Cargo.lock SHA-256, Rust/Cargo versions, source/binary hashes,
profile and environment overrides. Do not run `cargo update` or accept lock
drift. A cold dev generator is sufficient for this scope and is not packaged
`package-opt1` acceptance:

```bash
export CARGO_TARGET_DIR=/absolute/local/path/gentle-internal12-baseline-target
export CARGO_BUILD_JOBS=1
export CARGO_INCREMENTAL=0
export CARGO_PROFILE_DEV_DEBUG=0
cargo build --locked -j1 --bin gentle_examples_docs
BIN="$CARGO_TARGET_DIR/debug/gentle_examples_docs"
python3 -m unittest scripts.test_release_candidate -v
python3 -m unittest scripts.test_tutorial_checkouts -v
cargo test --locked -j1 --test release_version_consistency
cargo test --locked -j1 --lib historical_panel_
cargo test --locked -j1 --lib retained_panel_baselines_keep_valid_historical_provenance_and_hashes
cargo test --locked -j1 --lib retained_tutorial_artifact_normalization_preserves_generator_identity
cargo test --locked -j1 --lib workflow_examples_patz1_endpoint_and_sybr_panels_are_explicit_and_primer_only
cargo test --locked -j1 --lib export_promoter_artifact_manifest_
cargo test --locked -j1 --lib promoter_manifest_paths_match_in_direct_and_tutorial_replay
cargo test --locked -j1 --lib rewrite_example_paths_handles_promoter_and_handoff_outputs
cargo check -q --locked -j1
"$BIN" --check
"$BIN" tutorial-catalog-check
"$BIN" tutorial-manifest-check
"$BIN" tutorial-check
git diff --check
```

The expected outcome is a successful current replay against unchanged `.11`
baselines. Do not regenerate them merely to make version numbers agree.
If replay exposes actual result drift, generate into a temporary directory,
retain the before/after diff and explain the scientific changes before proposing
replacement baselines. Never broaden the comparison exclusion to hide them.

The integration also fixed LF protection for four TP53/TP73 evidence files and
manifest-relative promoter paths. Run the real canonical promoter workflow in
a disposable directory through its normal CLI route with documented fixture
prerequisites. Require six present components, zero missing required components
and references resolving beside the manifest. A cwd file must not mask a missing
manifest-relative component. The synthetic path unit test is not biological
artifact acceptance.

On the final committed candidate, run without an attributes overlay:

```bash
python3 scripts/check_tutorial_checkouts.py --binary "$BIN" --timeout-seconds 1800
python3 scripts/maintenance_chore.py session-close --plan docs/internal_12_tutorial_regeneration_glen_prompt.md
```

## Delivery And Boundaries

- Report the exact source/executable/final-commit SHAs, lock hash, profile and
  environment, commands/results, failures/skips and whether baseline bytes
  stayed unchanged. A stale executable is not current-replay acceptance.
- If code changes while resolving a failure, rebuild and rerun affected gates.
  A later documentation-only commit is not a fresh executable build; identify
  both revisions rather than claiming exact-binary equivalence silently.
- Only after successful execution, update the pending replay verdict in the
  roadmap/release notes and add evidence to the changelog on a reviewable bot
  branch. Return the commit and comparison link; never push Steffen's main.
- Leave `.11` tags/assets and `.12` `package-opt1` unchanged. No dependency or
  optimization changes, release publication, tag moves or workflow dispatch.
- Native Windows/macOS, package smoke, GUI/scientific acceptance and unrelated
  tutorial/GUI issues remain separate. Linux LF/CRLF simulation is not native
  Windows acceptance.

Codex has only made static source, syntax, formatting, attribute and hash
checks. No local tests, builds, Cargo check, GUI runs or benchmarks were run.
