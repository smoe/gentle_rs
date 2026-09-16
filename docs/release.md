# Release Process

This project publishes installable desktop packages for:

- macOS: `.dmg`
- Windows: `.zip` (contains `gentle.exe`)
- Linux: `.tar.gz` for x86-64, built on Ubuntu 24.04; GUI plus CLI/MCP/script
  entrypoints and tracked resources. See [Linux tarball quick start](linux_tarball.md).

Explicitly approved releases also publish container images through GitHub
Container Registry (GHCR); a tag push alone only runs build checks:

- headless CLI image:
  `ghcr.io/<owner>/<repo>:cli` and `ghcr.io/<owner>/<repo>:<tag>-cli`
- browser-served GUI image:
  `ghcr.io/<owner>/<repo>:gui` and `ghcr.io/<owner>/<repo>:<tag>`
- `latest` is updated only by an explicitly approved publish run and remains a GUI tag
- current image platform: `linux/amd64`

The current release workflow adds an actual Linux tarball; Debian, RPM and AppImage
packaging remain deferred. `linux_distribution=tarball` records the artifact
actually built, not a future intention. Until that workflow has passed on the
tag candidate, Linux download packaging remains an unverified release gate.

## Actions Runtime

CI, candidate resolution, installers and containers use action versions that
declare Node 24, rather than forcing Node 20 actions onto a different runtime.
The selected action versions require Actions Runner 2.327.1 or newer; the
GitHub-hosted runner observed during `.10` validation was 2.337.0.
Self-hosted runners must also satisfy the Node 24 OS/architecture requirements.
These are CI runner requirements, not new requirements for the GENtle binaries.
See [GitHub's Node 20 retirement notice](https://github.blog/changelog/2025-09-19-deprecation-of-node-20-on-github-actions-runners/).

This runtime-only migration uses `checkout@v5`, `setup-python@v6`, `cache@v5`,
`upload-artifact@v6`, `download-artifact@v7`, Docker setup/login v4, build/push v7,
metadata v6 and `softprops/action-gh-release@v3`. Artifact upload v5 and download
v6 only had preliminary Node 24 support and still declared Node 20. The existing
attestation v3 already delegates to Node 24 actions; the Rust-toolchain action
uses shell steps. Both remain unchanged. Artifact names, extraction layout,
cache keys, build profiles, candidate checks and publication approval remain
unchanged. The offline workflow-wiring test guards these reviewed references;
actual execution and package validation remain GitHub's responsibility.

## Candidate Approval

`v0.1.0-internal.10` remains unreleased pending Glen's exact-candidate readiness
verdict and the release owner's approval. A Git tag, draft release, successful
build or individual passing test is not release sign-off. Do not advance the
candidate to `.11` merely because a `.10` tag exists.

The existing `.10` tag points to `052cf125`, not the current development
candidate. Leave that tag unchanged unless the release owner explicitly
authorizes reconciliation after reviewing the accepted SHA. Packaging builds
the tag's revision; an older tag must not stand in for the candidate evaluated
by Glen.

## Build-Only Candidate Verification

After the tutorial fix and candidate changes are merged, record one clean
candidate commit and use it for both workflows. Run these commands only when
that commit is available on GitHub and includes the candidate-verification
helper. The `--ref main` selects the workflow definition, while `candidate_sha`
selects the exact source to build; both revisions are retained in the receipts.
Keep `main` frozen at the candidate while dispatching the acceptance runs.

```bash
CANDIDATE_SHA=$(git rev-parse HEAD)
gh workflow run release.yml --ref main \
  -f tag=v0.1.0-internal.10 -f candidate_sha="$CANDIDATE_SHA" -F publish=false
gh workflow run container.yml --ref main \
  -f tag=v0.1.0-internal.10 -f candidate_sha="$CANDIDATE_SHA" -F publish=false
```

Manual runs default to **build-only**. They require the full 40-character commit
SHA and a version label matching `Cargo.toml`; branch names and abbreviated
SHAs are rejected. The version label need not be an existing tag. In particular,
the older `.10` tag does not prevent evaluating the new `.10` candidate, and no
tag is created or moved. Build/check jobs have read-only repository permission;
the write-capable publication jobs are skipped.

Download the Actions artifacts from those specific run IDs, not from a generic
"latest" run. Installer runs retain three actual packages, per-platform build
receipts and an aggregate metadata receipt. Container checks build/load both
runtime targets, run their CLI/MCP entrypoints with container networking disabled,
and retain local image IDs/digests and a container receipt. Docker build records
remain with the Actions run; container validation does not upload images to GHCR.
These are packaging/entrypoint checks, not graphical or scientific acceptance.

The shared `gentle.release_candidate.v1` receipt records candidate SHA, lockfile
hash, workflow revision, version label and `validate_only`/`publish` mode.
Installer receipts also bind the actual archive digest, toolchain, release
profile and script features. Collection compares every receipt with the selected
candidate, not merely with the other receipts; missing, stale, mixed-mode or
modified packages fail closed. Docker retains its existing `release-fast`
profile and Debian `forky` build arguments, distinct from the installers'
`release` profile. Neither profile nor production optimization is changed here.
`.gitattributes` keeps `Cargo.lock` byte-identical under LF and CRLF checkouts;
candidate receipts continue to hash actual bytes, not normalized text. The
offline release-policy tests exercise both Git checkout policies and reject
unprotected conversion or subsequent lockfile edits.

Only after Glen's readiness verdict and release-owner approval may an owner
explicitly request `publish=true`, with the same SHA and an already-existing tag
that points to it. Publishing a GitHub Release is also an explicit publication
event. The tag is checked before building and again in the publication job;
an older/moved tag fails rather than silently publishing another revision.
Publication is a new run, not promotion of previously retained validation
artifacts: its own checks/receipts must pass, and external dependency changes
can produce different bytes even for the same source SHA. Do not treat earlier
artifact hashes as the hashes of the published packages.

No workflow dispatch, tag change, GHCR push or GitHub Release publication is
authorized merely by these instructions. Both real CI workflows still need to
run successfully on the selected candidate.

## Source Archive Exclusions

Tutorial runtime-generated files are committed for GitHub browsing under:

- `docs/tutorial/generated/`

They are excluded from `git archive` source bundles via:

- `.gitattributes`:
  - `docs/tutorial/generated export-ignore`
  - `docs/tutorial/generated/** export-ignore`

Local guard check:

```bash
archive_path=/tmp/gentle-src.tar
git archive --format=tar HEAD > "$archive_path"
tar -tf "$archive_path" | grep '^docs/tutorial/generated/' && echo "unexpected"
```

## Workflows

- CI validation workflow: `.github/workflows/ci.yml`
  - Runs checks/tests on pushes to `main` and pull requests.
  - Does not publish release assets.
- Container workflow: `.github/workflows/container.yml`
  - Runs no-push build checks for both Debian-first runtime targets
    (`runtime-cli`, `runtime-gui`) on tag pushes, published-release events and
    manual dispatch, using the resolved immutable candidate commit.
  - Publishes `linux/amd64` GHCR images only on a published-release event or
    an explicit manual `publish=true` run after matching tag/SHA checks:
    `:cli` / `:<tag>-cli` for headless use and `:gui` / `:<tag>` for GUI use.
  - A tag push alone never logs into GHCR or moves `latest`.
- Release workflow: `.github/workflows/release.yml`
  - Triggered automatically when a GitHub Release is published.
  - Can also be run manually via `workflow_dispatch` with inputs:
    - `tag`
    - `candidate_sha` (full commit SHA, required for manual runs)
    - `publish` (boolean, defaults to `false`)
    - `linux_distribution` (currently only `tarball`)
  - Builds macOS and Windows installers and a Linux tarball, runs smoke checks,
    and retains Actions artifacts in both modes. Publication is separate and
    explicit, never a side effect of a build-only run.
  - Resolves the candidate once and checks out that SHA in every subsequent
    job. Release events bind to the event SHA and verify the tag still matches.
  - Forces `CARGO_TARGET_DIR=target` for the release job so installer artifact
    discovery stays inside the checked-out workspace even if a developer or CI
    wrapper has overridden the default local target directory.
  - Caches the Cargo registry, but intentionally builds installer artifacts
    from a fresh per-run local `target/` tree rather than restoring a cached
    `target/` directory, so stale tag-build outputs cannot masquerade as a
    successful package build.
  - Logs the immediate `${CARGO_TARGET_DIR}/release` output layout after each
    platform build so missing bundle/binary regressions fail close to the build
    step rather than only during later packaging collection.
  - Also publishes a release-attributes JSON file:
    - `gentle-<tag>-release-attributes.json`
    - schema marker: `gentle.release_attributes.v1`
    - includes actual `linux_distribution`, common revision and lockfile hash,
      plus artifact names, sizes and SHA-256 digests.
  - Checks tag/package identity and builds the locked script-enabled binaries
    before packaging. Per-platform build receipts bind tag, full revision,
    lockfile hash, toolchain and profile; publication rejects mismatched or
    missing platform receipts and missing archive formats.
- Shared identity workflow: `.github/workflows/release-candidate.yml` calls
  `scripts/release_candidate.py`. Its fail-closed policy is tested offline by
  `python3 -m unittest scripts.test_release_candidate -v` on every CI push/PR.

## Artifact Naming

Release assets are normalized to:

- `gentle-<tag>-macos-<arch>.dmg`
- `gentle-<tag>-windows-<arch>.zip`
- `gentle-<tag>-linux-x64.tar.gz`
- `gentle-<tag>-<platform>-<arch>.build.json`
- `gentle-<tag>-release-attributes.json`

Example:

- `gentle-v0.1.0-macos-arm64.dmg`
- `gentle-v0.1.0-windows-x64.zip`
- `gentle-v0.1.0-linux-x64.tar.gz`
- `gentle-v0.1.0-release-attributes.json`

## Local Pre-Tag Smoke Checklist

Before pushing an internal or public release tag, run a release-shaped local
smoke pass that matches the packaging feature set rather than the lean default
developer build.

Required local matrix:

```bash
cargo check -q --locked
cargo test --locked --workspace
cargo test --locked -q --test release_version_consistency
cargo build --locked --release --features script-interfaces --bins
target/release/gentle --version
target/release/gentle_cli capabilities
target/release/gentle_js --version
target/release/gentle_lua --version
target/release/gentle_examples_docs --check
target/release/gentle_examples_docs tutorial-check
target/release/gentle_examples_docs tutorial-manifest-check
target/release/gentle_examples_docs tutorial-catalog-check
target/release/gentle_mcp --help
target/release/gentle_publication_report --help
```

Record the full candidate SHA, clean-tree status, lockfile hash, toolchain and
features/profile before running gates. Use the same checkout and built
binaries for the whole evidence set. Results from nearby commits cannot be
combined into an exact-candidate pass; the `.10` release notes contain the
pending scientific, GUI, tutorial and packaging gate ledger. Glen owns Linux
acceptance; Windows/macOS and both Docker runtime targets require CI evidence.

In particular, all three Linux tutorial smoke chapters must finish. The compact
Simple-PCR input bounds the work; it does not waive the ten-minute compute
budget, the non-empty primer report or the ROI-flanking checks.

Release-note expectations for that smoke pass:

- record pass/fail per command in the versioned internal release-notes document
  under `docs/release_notes/`
- call out any intentionally skipped entrypoint or known failure explicitly
- when exercising the pre-release CUT&RUN proof path, follow
  [`docs/cutrun_release_smoke.md`](cutrun_release_smoke.md) and record whether
  the proof used processed evidence only or also included ROI read
  interpretation, plus whether the optional DNA-window GUI report inspector was
  used to review the same regulatory-support payload

Release-workflow assumptions to re-check before tagging:

- macOS installer output remains `.dmg`
- Windows installer output remains `.zip`
- Linux tarball is actually built, extracted and checked before publication

## Standard Tagged Release

1. Ensure `main` is green in CI.
2. Create and push a version tag:
   - `git tag vX.Y.Z`
   - `git push smoe vX.Y.Z`
3. Publish a GitHub Release for that tag.
4. Wait for `Release Installers` workflow completion.
5. Verify GitHub Release contains all three desktop artifacts, their build
   receipts and the release-attributes inventory. Check the separate container
   workflow for both runtime images from the same tag.

## Manual Release Re-run

Use Actions → `Release Installers` → `Run workflow` and provide:

- `tag`: the version label matching the candidate's `Cargo.toml`
- `candidate_sha`: the full 40-character commit to rebuild
- `publish`: leave `false` for build-only downloadable Actions artifacts
- `linux_distribution`: `tarball` (the only implemented Linux download format)

Only an explicitly approved `publish=true` run updates release assets; its tag
must already point to the candidate SHA. A build-only run does not require a
tag to exist and never creates, moves or publishes one.

A tag push triggers container build checks, not `Release Installers`. To rebuild
Windows/macOS/Linux packages after retagging, dispatch this workflow explicitly
with the intended SHA and mode. Container publication remains separate and
requires its own approved publication event.

## Smoke Checks in Release Workflow

- macOS:
  - mounts the `.dmg`
  - verifies an `.app` bundle exists with `Contents/Info.plist`
- Windows:
  - extracts the ZIP package
  - verifies `gentle.exe` is present and non-empty
- Linux:
  - extracts the tarball into a fresh directory and verifies its file checksums
  - verifies the retained full revision matches the tag checkout
  - runs every packaged entrypoint's version/help/capabilities smoke from the
    extracted directory and checks required resource files
  - does not claim live graphical, scientific or cross-distribution acceptance

## Rollback / Recovery

If a release artifact is broken:

1. Fix on `main` and rerun the new candidate's gates.
2. Publish a new version/tag containing the fix. Do not move an existing tag
   or attach binaries from a different commit to it.
3. A same-tag workflow rerun is appropriate only for transient infrastructure
   failures: it rebuilds that tag, not the newer `main` fix. If needed, remove
   broken assets before republishing verified replacements for that same tag.

## Internal Release Notes

For internal tags, keep versioned release-note documents under
`docs/release_notes/`.
