# Glen Prompt: Close The Internal.12 Tutorial Regeneration Blocker

Prepared 2026-10-02, after integrating the tutorial review into local `main`.
This is a handoff prompt, not an executed job or an acceptance record.

## Request

Please close the version-bound tutorial regeneration gap in `smoe/gentle_rs`
on your Linux machine, in a separate bot branch/worktree. Do not run builds on
Steffen's Mac. Keep this task limited to regeneration, the integration
regressions below and their evidence; no installer build or benchmarking.

Use the integrated source, not your old review branch or a stale upstream main:

- Integration: `a2a88995616e64795a7ec024951d97617f5e6c1b`.
- Rollover parent: `e89a6932c655b9af93e88c363447d77ba3514755`.
- Your merged review head: `ae21fc0542a030132658a828ad43ffcac0006839`.
- Workspace version: `0.1.0-internal.12`.

Codex has not pushed this integration. First establish that the integration
commit is available through an owner-provided/published ref. If it is not,
report that prerequisite rather than substituting the old `4e7adf6e` main or
your `ae21fc05` branch. Record the actual full source SHA and verify that the
integration is its ancestor. A subsequent documentation-only handoff commit
may be included; do not silently incorporate other new work.

Read `AGENTS.md`, `docs/release.md` (Post-Release Development Version),
`docs/testing.md`, `docs/roadmap.md` and the `.12` release-note draft first.

## Known Failure

Cargo and the release metadata already declare `.12`, but these retained
reports still contain `.11` in `selection_audit_generator_revision` and
`gentle_version` (including nested occurrences):

```text
docs/tutorial/generated/artifacts/patz1_transcript_assay_panels_cli/artifacts/patz1_routine_common_region_screen.report.json
docs/tutorial/generated/artifacts/patz1_transcript_assay_panels_cli/artifacts/patz1_sybr_juc_panel.report.json
docs/tutorial/generated/artifacts/patz1_transcript_assay_panels_cli/artifacts/patz1_endpoint_end_matrix.report.json
```

The fast guard
`scripts.test_release_candidate.WorkflowWiringTests.test_replayed_tutorial_reports_use_current_package_version`
rejects this intermediate tree. This is a known preparation blocker, not
evidence that the new build profile or the scientific calculation failed.
Confirm and retain that initial failure before changing generated files.

## Regenerate, Do Not Relabel

1. Record source SHA, clean worktree status, Cargo.lock SHA-256, Rust/Cargo
   versions, build profile and explicit environment overrides. Use a dedicated
   local target directory, not network storage. Do not run `cargo update` or
   accept lockfile drift.
2. Build `gentle_examples_docs` from the integrated `.12` source, locked and
   with one job. A cold dev build is sufficient for this narrowly scoped
   regeneration; it is not `package-opt1` acceptance. Do not reuse an old `.11`
   executable. For example, with an explicitly chosen absolute target path:

```bash
export CARGO_TARGET_DIR=/absolute/local/path/gentle-internal12-regeneration-target
export CARGO_BUILD_JOBS=1
export CARGO_INCREMENTAL=0
export CARGO_PROFILE_DEV_DEBUG=0
cargo build --locked -j1 --bin gentle_examples_docs
BIN="$CARGO_TARGET_DIR/debug/gentle_examples_docs"
GENERATED_DIR=$(mktemp -d)
"$BIN" tutorial-generate --tutorial-output "$GENERATED_DIR"
```

3. Compare the complete generated tree with `docs/tutorial/generated` before
   copying anything back. Copy the legitimately regenerated affected reports
   and generator-produced `report.json` checksum ledger. No search-and-replace
   of version strings, manual provenance edits, changed expected hashes to
   conceal drift, or normalization of raw hash-bound inputs.
4. Compare primer sequences/coordinates, scores, class coverage and selection
   decisions before and after regeneration. Explain any scientific difference;
   do not dismiss it as version churn. Investigate any additional generated-file
   changes and keep only those required by the integrated source. Preserve the
   reviewed tutorial prose, screenshots and their historical evidence.

## Verify The Integration Repairs Too

The merge corrected two review findings which have not been execution-tested
locally. Keep these fixes rather than reintroducing the former behavior:

- Four TP53/TP73 GenBank/JSON/SVG evidence files now have scoped `text eol=lf`
  rules. Their original hashes are unchanged. The fast checkout regression
  verifies LF/CRLF and includes a missing-rule negative control.
- Promoter component paths are relative only to the output manifest directory
  (or explicitly absolute). The canonical workflow now uses sibling filenames;
  the tutorial path adapter preserves these references. There is no fallback
  to a same-named file under cwd. Missing artifacts must remain missing.

Run and retain the following focused checks, keeping the same build environment:

```bash
python3 -m unittest scripts.test_tutorial_checkouts -v
python3 -m unittest scripts.test_release_candidate -v
cargo test --locked -j1 --test release_version_consistency
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

Also execute the real canonical promoter workflow through its normal CLI route
in a disposable working directory, with the documented fixture prerequisites.
Confirm six present components, zero missing required components, and that each
returned relative path resolves beside the saved manifest. Keep this separate
from the new synthetic path regression, which uses placeholder content and
does not validate the biological artifacts.

After committing the regenerated files, run the exact-committed checkout gate,
without an attributes overlay:

```bash
python3 scripts/check_tutorial_checkouts.py --binary "$BIN" --timeout-seconds 1800
python3 scripts/maintenance_chore.py session-close --plan docs/internal_12_tutorial_regeneration_glen_prompt.md
```

Record both the executable's source SHA and the final artifact-commit SHA. If
executable sources change while fixing a failure, rebuild before claiming
acceptance. A generated-files/docs-only commit is not a new binary build; say
so explicitly. LF/CRLF simulation on Linux is not native Windows acceptance.

## Scope And Delivery

- Preserve the `.12` `package-opt1` recipe: opt-level 1, LTO off, 256 codegen
  units, assertions and overflow checks on, unwinding, no stripping,
  incremental=false, debug=0. No further optimization or dependency updates.
- Leave the published `.11` tag, assets and receipts unchanged. Do not push to
  Steffen's main, move tags, publish a release or dispatch packaging workflows.
- Do not fold in 08.02, glossary/decision changes, GUI stack-overflow fixes,
  sibling worktree edits or unrelated tutorial improvements.
- Only after the relevant gates pass, replace the pending blocker in the
  roadmap and `.12` release notes with the precise regeneration verdict, and
  add a changelog entry. Leave package/native GUI/scientific acceptance pending.
- Return a reviewable bot-branch commit and comparison link, the complete
  regeneration diff summary, exact commands/results, source/binary/lock hashes,
  profile/environment, remaining failures/skips and native-platform gaps.

Codex's integration checks were static only: formatting, Python/JSON syntax,
Git attributes, the four unchanged evidence hashes and whitespace. No local
Rust/Python test suite, Cargo check/build, GUI run or benchmark was executed.
