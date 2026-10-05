# Internal .12 Candidate Preparation

This is a handoff, not a selected candidate or release verdict. Inspection began
at `2829eec7` on 2026-10-05. No push, dispatch, tag, upload, publication, runner
change or optimization increase is authorized. Preserve private outputs and
the published `.11` artifacts.

## Changed Files

- CI logic/regressions: `src/workflow_examples.rs`,
  `scripts/test_tutorial_checkouts.py`, `src/workflow_examples/agent_parity.rs`
  and its `tests.rs`.
- Contract/testing guidance: `docs/architecture.md`, `docs/protocol.md` and
  `docs/testing.md`.
- Status/handoff: `docs/CHANGELOG.md`, `docs/roadmap.md`,
  `docs/release_notes/release_notes_v0.1.0-internal.12.md`,
  `docs/clawbio_gentle_integration_onepager.md` and this document.

## Corrected CI Evidence

- A later [run 37292637621](https://github.com/smoe/gentle_rs/actions/runs/37292637621)
  at `462b5ffe848e53bd287545ded11154bb9bee0288` is also red. Its fast checkout
  gate reports the unbound mock, and both selected desktop jobs report the
  outdated 19-versus-20 TSS case count. `42867017` already fixes those failures.
  Windows additionally fails
  `rejects_raw_absolute_traversal_and_windows_paths_before_normalizing`:
  `/tmp/request.json` produces `file_unavailable` instead of `unsafe_path`.
  That rooted path lacks a drive prefix, so Windows `Path::is_absolute()` is
  false. Reject roots with `Path::has_root()` before joining/canonicalizing,
  covering guide, template and payload bindings as well as direct `@file`.
  Preserve raw traversal, Windows spelling and symlink guards; do not change
  command grammar or regenerate the catalog. A new native Windows verdict is
  still required; this red run is not evidence about the unpushed repair.
  The additional repair is local commit
  `a086bbda5d372f1caf00c6034d88a3862f029910`, not a selected candidate.
- [Run 37235516663](https://github.com/smoe/gentle_rs/actions/runs/37235516663)
  at `f523a757cd2a2861196b1a5c0b36df67983ce0e6` completed with Linux and Windows
  failing the two catalog tests; their preceding tutorial-drift checks passed.
  macOS was skipped. Windows completed at 2026-10-04 22:38:51 UTC; its full suite
  was skipped, not passed.
- The review reversed the catalog diff: `git show f523a757:docs/tutorial/catalog.json`
  contains the **graphic** reason. The generator/log warning selects **source
  JSON**. Both dependencies last changed in `ae21fc05` at the same UTC commit
  date; source priority wins. `19889dae` already repairs that projection.
  There is no evidence that CI used a newer PNG mtime in that failure.
- The separate latent defect was real: failed/empty Git lookup could fall back
  to checkout mtime. Review projections now reject unavailable history,
  preserve missing-graphic/newest-date/priority ordering and use relative paths.
  The committed catalog's age reference is HEAD's commit date, while live age
  warnings remain current-date-based. No human review is relabelled fresh.
- The checkout early-failure regression introduced in `19889dae` referenced an
  unbound mock. Bind the mock and retain its exact stop-after-failure assertion.
  The TSS collection export added the twentieth agent-parity case; update the
  count assertions, not the mutation guard.
- `42867017184c6fd65ebbd479c69bc0e73382ba20` is the local CI repair, not a
  selected or pushed candidate. The old age-only test fixture also needed a
  committed synthetic workflow; repair its evidence instead of restoring the
  silent-history fallback.
- [Windows job 109794599955](https://github.com/smoe/gentle_rs/actions/runs/36686854112/job/109794599955)
  at `a51bbc064ea29d4e7a03c559c5a6e0a6f8e1bd29` completed the unfiltered
  `cargo test -q --locked --workspace` step successfully under normal parallelism.
  This covers the non-ignored `gentle_cli` forwarded-services tests and closes
  the historical JASPAR-lock recheck. It is **not** current `.12` acceptance.
  The separate CertUtil SHA-1 sharing-violation note remains open.

## Local Verification

Host: macOS, Rust `1.100.0-beta.1 (e3feeb59c 2026-09-27)`. These are functional
working-tree checks, not package, GUI, performance or frozen-candidate receipts.
No catalog, historical report or retained screenshot was regenerated.

### Integrated-Head Recheck

The follow-up started at `7387beb4` and fast-forwarded this clean branch to
local `main` `fb9e6e00e4c6b6d7994b7979326b5e999f047ab6`, which already includes
the two scoped commits `42867017` and `d1883dde`. A fetch still placed public
`origin/main` at `462b5ffe`; no new push or dispatch was performed.

At the unchanged integrated HEAD, both original catalog tests pass separately
and from a fresh full clone at
`/private/tmp/gentle-internal12-current-fb9e6e00` (not shallow, not a worktree
copy). Reuse the just-built test executable by absolute path with the clone
as CWD and `RUST_MIN_STACK=16777216`; its SHA-256 before the rooted-path repair
is `02fc0c69ed00787d934e8eaa874b080d93b79c944ed21f4bd5c4891c6c1443a0`.
The lockfile SHA-256 remains
`3dfe44c7a08bf779e32289cf1312103ea1fb5883e44bba5f6e1912ce3ed77745`.
The combined checkout/walkthrough/publication Python checks run 28 tests:
25 pass and three have expected platform skips. Formatting and locked Cargo
check also pass. These are pre-repair checks, not evidence that Windows has
accepted the new path guard.
Post-repair results are recorded separately below; do not merge their identity
with the earlier 92-test verification.

### Earlier Scoped Verification

```bash
python3 -m unittest scripts.test_tutorial_checkouts -q
cargo test -q --locked --lib \
  workflow_examples::tests::tutorial_catalog_check_passes_on_committed_tree \
  -- --exact --test-threads=1
cargo test -q --locked --lib workflow_examples -- --test-threads=1
cargo test -q --locked --lib \
  workflow_examples::tests::tutorial_review_manifest_stale_date_warns \
  -- --exact --test-threads=1
cargo check -q --locked
rustfmt --edition 2024 --check --config skip_children=true src/workflow_examples.rs
git diff --check
```

- Python checkout regressions: **18 passed**, including LF/CRLF date ties and
  stop-after-first-failure assertions. Catalog checks and independent history/
  age regressions pass.
- The first 92-test serial workflow run had **91 passes / 1 failure**: the
  age-warning fixture had no Git-bound workflow. Its repaired fixture passes
  independently; the corrected full replay passes **92/92**, with no ignored
  tests. This is the requested workflow-example subset, not the full workspace.
- A direct harness launch without Cargo's environment aborted on the existing
  default-stack boundary. Reuse `.cargo/config.toml`'s **unchanged**
  `RUST_MIN_STACK=16777216`; no new stack policy was introduced. The full replay
  uses the rebuilt executable with SHA-256
  `5926b2b7b9b4c84e765641acaa1628e7e810d50e36a76c5a03831a114ecee46c`:

```bash
RUST_MIN_STACK=16777216 \
  target/debug/build/GENtle/39c6fd870ce74af2/out/gentle-39c6fd870ce74af2 \
  workflow_examples --test-threads=1
```

Fresh-clone verification used a full clone, **not** a worktree copy:

```bash
git clone --no-hardlinks --quiet /Users/u005069/.codex/worktrees/faae/gentle_rs \
  /private/tmp/gentle-internal12-review-42867017
git -C /private/tmp/gentle-internal12-review-42867017 \
  rev-parse --is-shallow-repository HEAD
```

The clone reports `false` and the full local repair SHA above. From that clone's
working directory, invoke the same executable by absolute path, with
`RUST_MIN_STACK=16777216`, once for each exact test filter below and
`--exact --test-threads=1`. Both pass with fresh checkout timestamps:

- `workflow_examples::tests::tutorial_catalog_check_passes_on_committed_tree`
- `workflow_examples::tests::tutorial_source_units_generate_committed_catalog`

The clone checks use the same compiled generator logic, not an older `.11`
binary and not a second monolithic build. Authoritative generation requires
available Git history; CI retains `fetch-depth: 0`. No new shallow-history
qualification is claimed. The macOS linker emitted the known non-blocking
`__eh_frame` compact-unwind warning; it did not fail the focused test build.

Session-close: `python3 scripts/maintenance_chore.py session-close --plan
docs/roadmap.md` reports **4 OK / 2 warnings / 0 failures**. The warnings are
pre-commit documentation edits/untracked artifacts and the manual plan-fidelity
reminder. The post-commit check has the same totals, with only unrelated
`outputs/` remaining untracked. Preserve that directory; all task edits are the CI repair or
this documentation handoff. The final locked check and 92-test serial replay
are green. Manual plan-fidelity review confirms these scoped changes; the
historical CI diff correction, extra mock/count repairs and age-fixture repair
are recorded above, not silently treated as the original timestamp hypothesis.

### Rooted-Path Repair Verification

On Darwin 27 arm64 with Rust/Cargo 1.100.0-beta.1, the patched tree based on
`fb9e6e00` passes **93/93 workflow-example library tests**, including all
**11 agent-parity tests**, with no ignored cases. The latter tests also pass
as a separate targeted run. This uses CI's serial library-test order; the
full workspace/all-target invocation was not run locally.

```bash
cargo test -q --locked --offline -j 1 --lib workflow_examples::agent_parity \
  -- --test-threads=1
cargo test -q --locked --offline -j 1 --lib workflow_examples \
  -- --test-threads=1
python3 -m unittest scripts.test_tutorial_checkouts \
  scripts.test_tutorial_walkthroughs \
  scripts.test_publish_tutorial_gui_screenshots -q
cargo check -q --locked --offline -j 1
cargo fmt --check
git diff --check
```

The Python run has **25 passes / 3 expected platform skips**. Locked Cargo
check, formatting and whitespace checks pass. The linker emits the existing
non-blocking `__eh_frame` compact-unwind warning. The rebuilt test executable
SHA-256 is `510c1e1157ba303d411d7ceb62096b2114bef49d4fb64e5ae5460886adb6dc05`.
Both original catalog tests pass separately in the worktree and again from
the same fresh clone, using this rebuilt executable by absolute path and
the existing 16 MiB test-thread stack. Catalog and manifest bytes, historical
reports, owner prompt files, dependencies and workflows remain unchanged.

These are pre-commit functional checks, not exact packaged-candidate receipts:
the clone's data remains bound to `fb9e6e00`, and the executable was built from
the patched working tree before `a086bbda`. Native Windows/Linux reruns, the
full workspace suite, package/container builds, benchmarks and live GUI or
scientific acceptance remain unrun. Rerun frozen-candidate gates after the
owner selects and pushes the final SHA; do not relabel this evidence.

## First package-opt1 Build-Only Runs

Only after green CI on one clean, pushed, owner-approved SHA, use the commands
from [Build-Only Candidate Verification](release.md#build-only-candidate-verification).
`git rev-parse HEAD` resolves the full 40-character source SHA, not a branch or
abbreviation. Freeze that source and the workflow ref; retain their identities
separately. These commands are instructions for the owner and have **not** run:

```bash
CANDIDATE_SHA=$(git rev-parse HEAD)
WORKFLOW_REF=main # already pushed at the frozen candidate
gh workflow run release.yml -R smoe/gentle_rs --ref "$WORKFLOW_REF" \
  -f tag=v0.1.0-internal.12 -f candidate_sha="$CANDIDATE_SHA" -F publish=false
gh workflow run container.yml -R smoe/gentle_rs --ref "$WORKFLOW_REF" \
  -f tag=v0.1.0-internal.12 -f candidate_sha="$CANDIDATE_SHA" -F publish=false
```

Retain both run URLs/repository identities, logged candidate/workflow SHAs,
toolchain/lock hashes and recipe settings. Do not merge nearby revisions' passes.

| Job evidence | What it can close |
| --- | --- |
| Three native archives ending `-package-opt1.dmg`, `-package-opt1.zip`, `-package-opt1.tar.gz`; per-platform `.build.json` and release-attributes inventory; exact archive hashes; successful extracted-package smokes | Roadmap's fresh `.12` native package recipe/entrypoint proof and release-notes acceptance item 2's native portion |
| `gentle.release_candidate.v1` and `gentle.container_build.v1`, local image ID/digest, three headless entrypoints, runtime-library/assets checks, unprivileged/no-network checks and actual nonempty RNAPKIN SVG/PNG smokes | Roadmap's changed `runtime-cli` recipe gate and acceptance item 2's container portion |
| Compiler peak RSS where actually measured, pressure/swap/cgroup logs, build time and package sizes | Resources for diagnosing this recipe; never a performance or GUI pass |

If a metric is absent, record **unavailable** and retain the logs; successful
packaging does not measure compiler peak RSS by itself. Glen's same-SHA repeated
`dev`/`package-opt1` startup/PATZ1 comparison (release-notes item 3), scientific
checks and live GUI acceptance (item 4) remain open after these jobs. Do not
dispatch a benchmark, change timeouts/parallelism or start B0 extraction in
response to an unmeasured memory hypothesis. A later hosted kill gets a separate
log-based diagnosis and owner-approved plan.

The current container smoke creates the RNAPKIN files inside a removed container;
it does not upload those SVG/PNG bytes or compiler peak-RSS measurements. Retain
the successful smoke log, image identity and pinned helper-lock receipt, and
record those retention gaps explicitly. Do not invent file hashes or treat
absent resource measurements as zero.

## ClawBio Compatibility Decision

No mode is removed. Earliest-removal `.11` was an eligibility marker, not an
automatic deadline; the owner must decide retention/removal for `.12`.
The 25 deprecated modes in `DEPRECATED_REQUEST_MODE_EQUIVALENTS` are:

| Family | Modes |
| --- | --- |
| Primer | `primer-preflight`, `primer-seed-from-feature`, `primer-seed-from-splicing`, `primer-design`, `primer-report-list`, `primer-report-show`, `primer-report-export` |
| qPCR | `qpcr-seed-from-feature`, `qpcr-seed-from-splicing`, `qpcr-design`, `qpcr-report-list`, `qpcr-report-show`, `qpcr-report-export` |
| cDNA/transcript panel | `cdna-pcr-test`, `cdna-qpcr-test`, `transcript-qpcr-panel` |
| Restriction handoff | `restriction-cloning-pcr-handoff`, `restriction-cloning-pcr-handoff-seed`, `restriction-cloning-vector-suggestions`, `restriction-cloning-handoff-list`, `restriction-cloning-handoff-show`, `restriction-cloning-handoff-export` |
| Other shell normalizers | `pcr-protocol-cartoon`, `exon-skip-plan`, `exon-skip-materialize` |

Call sites and regression scope:

- `integrations/clawbio/skills/gentle-cloning/gentle_cloning.py`: the deprecation
  dictionary, `_request_mode_warnings`, validation and shared-shell request
  routing (including the exon-skip, cDNA and restriction handoff branches).
- `integrations/clawbio/skills/gentle-cloning/INTENTS.json` and its `examples/`
  still expose compatibility requests; review every affected route before any
  removal rather than changing the warning alone.
- Documentation: the skill `SKILL.md`/`README.md` and
  [mode/deprecation ledger](clawbio_gentle_integration_onepager.md#deprecation-table).
- `integrations/clawbio/skills/gentle-cloning/tests/test_gentle_cloning.py`:
  supported-mode inventory plus exon-skip, qPCR seed, cDNA and PCR-cartoon
  command/confirmation regressions. Descriptor/example parity remains required.

## Ubuntu Runner Decision

The [Linux annotation](https://github.com/smoe/gentle_rs/actions/runs/37235516663/job/111533721645)
states that `ubuntu-latest` starts migrating to Ubuntu 26 on 2026-10-19 and
links [runner-images issue 14748](https://github.com/actions/runner-images/issues/14748).
Current consumers are:

| Workflow | `ubuntu-latest` jobs |
| --- | --- |
| `ci.yml` | `release-policy`, `select-platform`, `headless-linux`, `linux`, `ci-summary` |
| `container.yml` | `container-build-check`, `container-publish` |
| `release.yml` | `collect-packages`, `publish-release` |
| `release-candidate.yml` | `resolve` |

The Linux installer matrix already pins **`ubuntu-24.04`**; do not misdescribe it
as `ubuntu-latest`. A rolling host can change compiler/system dependencies,
graphics/font behavior, Docker/BuildKit or runtime libraries, even when source
and lockfile agree. Proposed owner decision: pin the listed hosts to
`ubuntu-24.04` for the `.12` candidate, then separately qualify Ubuntu 26.
No workflow file has been changed here.

## Remaining Decisions And Evidence

Owner decisions remain: push the local fixes; select/freeze a green candidate;
authorize the two build-only dispatches; choose ClawBio compatibility policy;
decide the runner pin. No tag, upload, package/performance acceptance or Glen
scientific/GUI verdict follows from this handoff. The full workspace suite,
native Windows/Linux reruns, macOS Actions/GUI checks, package/container builds
and benchmarks were not run; local macOS Rust checks are the scoped evidence
above, not a platform-release verdict.
