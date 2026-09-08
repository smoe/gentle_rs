# Release Process

This project publishes installable desktop packages for:

- macOS: `.dmg`
- Windows: `.zip` (contains `gentle.exe`)
- Linux: `.tar.gz` for x86-64, built on Ubuntu 24.04; GUI plus CLI/MCP/script
  entrypoints and tracked resources. See [Linux tarball quick start](linux_tarball.md).

Release tags also publish GitHub-downloadable container images through GitHub
Container Registry (GHCR):

- headless CLI image:
  `ghcr.io/<owner>/<repo>:cli` and `ghcr.io/<owner>/<repo>:<tag>-cli`
- browser-served GUI image:
  `ghcr.io/<owner>/<repo>:gui` and `ghcr.io/<owner>/<repo>:<tag>`
- `latest` is updated only from release-tag publishes and remains a GUI tag
- current image platform: `linux/amd64`

The `.11` workflow adds an actual Linux tarball; Debian, RPM and AppImage
packaging remain deferred. `linux_distribution=tarball` records the artifact
actually built, not a future intention. Until that workflow has passed on the
tag candidate, Linux download packaging remains an unverified release gate.

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
    (`runtime-cli`, `runtime-gui`) on tag pushes and manual dispatch.
  - Publishes `linux/amd64` GHCR images only when the workflow ref is a tag
    matching `v*`:
    `:cli` / `:<tag>-cli` for headless use and `:gui` / `:<tag>` for GUI use.
  - Moves `latest` only on release-tag publishes, as a GUI compatibility tag.
- Release workflow: `.github/workflows/release.yml`
  - Triggered automatically when a GitHub Release is published.
  - Can also be run manually via `workflow_dispatch` with inputs:
    - `tag`
    - `linux_distribution` (currently only `tarball`)
  - Builds macOS and Windows installers and a Linux tarball, runs smoke checks, and publishes assets
    to the corresponding GitHub Release.
  - Uses the release tag itself as the checkout/build target on `release`
    events, so publishing a release for an already-pushed tag still produces
    the desktop installers without needing a second tag push.
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
combined into an exact-candidate pass; the `.11` release notes contain the
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

- `tag`: an existing tag (for example `v0.1.0`)
- `linux_distribution`: `tarball` (the only implemented Linux download format)

This rebuilds installers and updates assets on that tag’s release.

Container publishing remains tag-driven through `.github/workflows/container.yml`
rather than the desktop-installer release workflow, so the tag push should
still happen before or alongside release publication.

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
