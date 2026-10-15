# Internal.13 Preparation Handoff

Prepared 2026-10-15. This records local implementation and future owner/Glen
actions, **not** dispatch, acceptance or publication authority. Work is on
`gentle_rs_3`, based on local main `22b06a9b`, not stale origin/main. Main is
checked out by another session; unrelated work there and `outputs/` here are
not changed. The previous unrelated goal stays paused.

## Local Commits And Evidence

| Commit | Scope | Verification boundary |
| --- | --- | --- |
| `9edba9e7` | Owner-selected .13 aim and gates | Documentation; no acceptance |
| `bb18cb44` | OpenSCAD floor/cap correction and three example sources | Six actual renders; Python semantic checks; Rust regression unrun |
| `71b63601` | Early-feature audit and native review instructions | Static inventory; native checks pending |
| `30bc1b9d` | Expression starter, review fixture, tutorial and MCP parity test | Offline source checks; Rust/CLI/MCP execution pending |
| `a5cfe580` | Explicit gene-set arrangement/storage projection and tutorial | Offline source checks; four Rust tests and runtime pending |

Full changed-file inventories are in those commits. The
[early-feature audit](audits/early_feature_review_20261015.md) records 21 feature
families, review-time test locations and all template/geometry limitations.
No inspected family is orphaned; the well-floor/cap defect is corrected, not
claimed independently accepted. Temporary reference PNGs remain outside Git.

### Runner Pin And Separate Fallback

- `pin-ubuntu-2404`: `5e12565acad7486c50b12a2d8c37313d74a8aa2f`, rebased onto
  local main `420d913e`, integrated into main at the owner's request on
  2026-10-15. Pins the ten rolling jobs; no Ubuntu 26 qualification. No push
  performed here; owner-approved push is required before candidate selection.
- `package-opt1-fallback`: `5f9f9e4a5b9aaafc3f51e74fd20c9a8ca991938d`, based on
  this task's `a5cfe580`. One optional commit switches native/container selection,
  receipts, paths, archive suffixes and policy tests together; both Cargo recipes
  and historical release bytes stay intact. All 103 offline packaging-policy
  tests pass. **Do not merge merely because it exists**: main remains opt2.

Use the fallback commit, not a blanket merge of an outdated branch, after
reconciling with the final integrated source. Retain the actual hosted opt2
failure; a new profile means a new candidate and fresh gates.

## External Completion Before Candidate Selection

No Cargo, rustc, GENtle binary, model, provider or real Docker build ran here.
Formatting is not compilation. Glen/CI must compile the new paths and execute
the authored regressions plus existing release/version and full CI gates.
Use a disposable clean checkout and target directory as described in the
[rollover handoff](internal_13_rollover_glen_handoff.md), preserving private data.

Focused commands for Glen/CI (not run locally):

```bash
cargo test --locked -j1 --lib rack_openscad_
cargo test --locked -j1 --lib gene_set_arrangement_
cargo test --locked -j1 --lib protein_expression_
cargo test --locked -j1 --test protein_expression_handoff_cli
cargo test --locked -j1 --test capability_registry_parity
cargo test --locked -j1 --test parity_matrix_freshness
cargo test --locked -j1 --test workflow_examples
cargo test --locked -j1 --test release_version_consistency
cargo check -q --locked -j1
cargo build --locked -j1 --bin gentle_cli --bin gentle_examples_docs
```

Before finalizing a candidate, use that freshly built docs binary to register
tutorial 03.10 in source/catalog metadata, keeping it manual/hybrid like 03.08;
review 09.03's source metadata as well. Generate/check workflow snippets for
both new executable starters, catalog, manifest and parity matrix. Neither
starter is a completed oracle for the user's task. Reconcile generated drift
and LF/CRLF checks in a focused follow-up commit, then select the candidate
containing it. Do **not** claim the unchanged old generated catalog covers 03.10.
Do not overwrite historical tutorial/scientific evidence to pass freshness.

Replay tutorials [09.03](tutorial/09-03_protein_expression_handoff_cli.md) and
[03.10](tutorial/03-10_gene_set_storage_projection.md) literally in fresh states.
Verify neutral/in-house/explicit-provider review and no state writes/provider
contact; for storage verify confirmation refusals, order/provenance, save/reopen,
same-session undo/redo and all three exports. Render the **newly exported** SCAD;
the reference renders are not a substitute. Exercise existing generic MCP/agent
confirmation without a live model. Native rack/label/3D checks follow the
[final handoff](glen_final_candidate_execution_handoff.md#internal13-early-features).

B remains minimal: no inspector/ranking, commercial contact or ordering. C
remains a digital placement draft: no freezer-box inventory, barcode workflow,
new GUI mapper or implicit biological source binding. Broader refinement is
separate from making these two literal paths pass.

## Owner Procedure

Commands below are instructions for the owner/main-maintaining thread; none
were executed here. Integrate this task into main without losing unrelated
work, complete the external/generated follow-up above, then:

1. Push approved main. Preserve the push run's Windows and sampled Unix results.
2. Confirm the integrated `5e12565a` pin is present in pushed main before
   October 19 and before candidate selection. An older published source still
   uses unpinned workflow definitions; the local merge is not hosted evidence.
3. Record the clean full candidate SHA, matching pushed main, lock hash and
   profile. Keep main at that tip while dispatching; if it moves, reselect and
   rerun rather than combining different workflow/source results.
4. Dispatch full Linux; retain Windows and **both** Unix verdicts at this SHA.
   If the push sampled Linux rather than macOS, also dispatch macOS explicitly.
5. Run build-only native installers/container; check each logged candidate and
   workflow revision, effective profile, package hashes, receipts and extracted
   smokes. Do not use a nearby SHA or stale published artifact.

```bash
# In the clean, integrated main checkout; owner authorization required.
git merge-base --is-ancestor 5e12565acad7486c50b12a2d8c37313d74a8aa2f HEAD
git push origin main
CANDIDATE_SHA=$(git rev-parse HEAD)
test "$CANDIDATE_SHA" = "$(git rev-parse origin/main)"
test -z "$(git status --porcelain)"
gh workflow run ci.yml -R smoe/gentle_rs --ref main -f platform=linux
# Run when macOS is not already covered at this exact SHA:
gh workflow run ci.yml -R smoe/gentle_rs --ref main -f platform=macos
gh workflow run release.yml -R smoe/gentle_rs --ref main \
  -f tag=v0.1.0-internal.13 -f candidate_sha="$CANDIDATE_SHA" -F publish=false
gh workflow run container.yml -R smoe/gentle_rs --ref main \
  -f tag=v0.1.0-internal.13 -f candidate_sha="$CANDIDATE_SHA" -F publish=false
```

6. If hosted opt2 packaging fails, retain run/logs and distinguish infrastructure
   from a product failure. Only with owner selection, integrate the **single**
   fallback commit, preserve the Ubuntu pin, push, and select a new SHA. Repeat
   steps 3-5 and all affected acceptance on opt1; an old opt1 pass is insufficient.
7. Give Glen that SHA for rollover/version/tutorial, early-feature/A/B/C checks
   and a scoped readiness verdict. The [release-note ledger](release_notes/release_notes_v0.1.0-internal.13.md#internal13-acceptance-ledger)
   starts entirely Pending. Latency, import performance, profile comparisons,
   live-agent and historical scientific ledgers carry forward without blocking
   .13, not as silently waived passes.
8. Only after all blockers and explicit owner approval: tag/publish using the
   [release process](release.md#standard-tagged-release). No tag command here
   authorizes publication. Preserve all existing tags/assets.
9. After owner-confirmed publication, perform the [.14 development rollover](release.md#post-release-development-version)
   immediately; synchronize versions/lock/readme/roadmap/changelog/new notes and
   run the external version guard. Do not advance .13 before confirmation.

## Local Verification

Executed offline, without invoking a real GENtle/Cargo executable:

```bash
PYTHONDONTWRITEBYTECODE=1 python3 -m unittest \
  scripts.test_storage_projection_tutorial scripts.test_early_feature_review \
  scripts.test_specialist_tutorials scripts.test_synthetic_biology_tutorial \
  scripts.test_tutorial_checkouts -q
PYTHONDONTWRITEBYTECODE=1 python3 -m unittest scripts.test_release_candidate \
  scripts.test_container scripts.test_package_desktop -q
git diff --check
python3 scripts/maintenance_chore.py session-close --plan docs/roadmap.md
```

First suite: **80 tests, 77 passed, three expected optional real-binary skips**.
Second: **103/103 passed** on the active opt2 branch. The same three packaging
modules pass **103/103** on the separate opt1 fallback; the pin branch's
`scripts.test_release_candidate` passes **46/46**. Package-smoke log lines in
these tests use synthetic executables, not actual extracted GENtle packages.
Local documentation link targets/Bash example syntax and whitespace pass;
the roadmap remains 400 lines. Session-close reports **4 OK, 2 warnings,
0 failures**: intentional working-tree/untracked `outputs/` entries and the
manual plan-fidelity reminder. Scope was reviewed manually; unrelated files,
historical evidence and lockfiles remain untouched.
