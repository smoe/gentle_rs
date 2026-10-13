# Glen: Validate The Internal.13 Development Rollover

Prepared 2026-10-13. This is a handoff, not an executed audit or publication
authorization. The owner confirmed `.12` at `5e6201a6`; its
[published ledger](release_notes/release_notes_v0.1.0-internal.12.md#publication-evidence)
remains separate from this development-version validation.

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
