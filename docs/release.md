# Release Process

This project publishes installable desktop packages for:

- macOS: `.dmg`
- Windows: `.zip` (contains `gentle.exe`)

Release tags also publish GitHub-downloadable container images through GitHub
Container Registry (GHCR):

- headless CLI image:
  `ghcr.io/<owner>/<repo>:cli` and `ghcr.io/<owner>/<repo>:<tag>-cli`
- browser-served GUI image:
  `ghcr.io/<owner>/<repo>:gui` and `ghcr.io/<owner>/<repo>:<tag>`
- `latest` is updated only from release-tag publishes and remains a GUI tag
- current image platform: `linux/amd64`

Linux installable packaging is intentionally deferred. Until Debian packaging is
ready, releases should default Linux distribution metadata to `tarball`.

Linux distribution intent is now tracked as a release attribute (`deferred`,
`deb`, `appimage`, `rpm`, or `tarball`) in a release-metadata JSON asset so
each release records the planned Linux channel explicitly. Current default:
`tarball`.

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
  - Builds both Debian-first runtime targets (`runtime-cli`, `runtime-gui`) on
    PRs / `main`.
  - Publishes `linux/amd64` GHCR images from tag pushes matching `v*`:
    `:cli` / `:<tag>-cli` for headless use and `:gui` / `:<tag>` for GUI use.
  - Moves `latest` only on release-tag publishes, as a GUI compatibility tag.
- Release workflow: `.github/workflows/release.yml`
  - Triggered by tag pushes matching `v*`.
  - Can also be run manually via `workflow_dispatch` with inputs:
    - `tag`
    - `linux_distribution` (`deferred|deb|appimage|rpm|tarball`)
  - Builds macOS and Windows installers, runs smoke checks, and publishes assets
    to the corresponding GitHub Release.
  - Forces `CARGO_TARGET_DIR=target` for the release job so installer artifact
    discovery stays inside the checked-out workspace even though normal local
    development builds default to the shared worktree target directory from
    `.cargo/config.toml`.
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
    - includes selected `linux_distribution`.

## Artifact Naming

Release assets are normalized to:

- `gentle-<tag>-macos-<arch>.dmg`
- `gentle-<tag>-windows-<arch>.zip`
- `gentle-<tag>-release-attributes.json`

Example:

- `gentle-v0.1.0-macos-arm64.dmg`
- `gentle-v0.1.0-windows-x64.zip`
- `gentle-v0.1.0-release-attributes.json`

## Local Pre-Tag Smoke Checklist

Before pushing an internal or public release tag, run a release-shaped local
smoke pass that matches the packaging feature set rather than the lean default
developer build.

Required local matrix:

```bash
cargo build --release --features script-interfaces
cargo run --release --bin gentle -- --version
cargo run --release --bin gentle_cli -- capabilities
cargo run --release --features js-interface --bin gentle_js -- --version
cargo run --release --features lua-interface --bin gentle_lua -- --version
cargo run --release --bin gentle_examples_docs -- --check
cargo run --release --bin gentle_examples_docs -- tutorial-check
cargo run --release --bin gentle_mcp -- --help
```

Release-note expectations for that smoke pass:

- record pass/fail per command in the versioned root release-notes document
- call out any intentionally skipped entrypoint or known failure explicitly

Release-workflow assumptions to re-check before tagging:

- macOS installer output remains `.dmg`
- Windows installer output remains `.zip`
- Linux release metadata defaults to `tarball`

## Standard Tagged Release

1. Ensure `main` is green in CI.
2. Create and push a version tag:
   - `git tag vX.Y.Z`
   - `git push smoe vX.Y.Z`
3. Wait for `Release Installers` workflow completion.
4. Verify GitHub Release contains both installer assets.

## Manual Release Re-run

Use Actions → `Release Installers` → `Run workflow` and provide:

- `tag`: an existing tag (for example `v0.1.0`)
- `linux_distribution`: planned Linux channel for this release metadata

This rebuilds installers and updates assets on that tag’s release.

Container publishing remains tag-driven through `.github/workflows/container.yml`
rather than the desktop-installer release workflow.

## Smoke Checks in Release Workflow

- macOS:
  - mounts the `.dmg`
  - verifies an `.app` bundle exists with `Contents/Info.plist`
- Windows:
  - extracts the ZIP package
  - verifies `gentle.exe` is present and non-empty

## Rollback / Recovery

If a release artifact is broken:

1. Fix on `main`.
2. Re-run release workflow for the same tag via manual dispatch, or
   create a new patch tag (`vX.Y.Z+1`) and publish normally.
3. If needed, remove incorrect assets from the GitHub Release UI before
   republishing.

## Internal Release Notes

For internal tags, keep a versioned release-notes document at repository root,
e.g.:

- `release_notes_v0.1.0-internal.2.md`
