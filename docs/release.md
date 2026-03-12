# Release Process

This project publishes installable desktop packages for:

- macOS: `.dmg`
- Windows: `.zip` (contains `gentle.exe`)

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
- Release workflow: `.github/workflows/release.yml`
  - Triggered by tag pushes matching `v*`.
  - Can also be run manually via `workflow_dispatch` with inputs:
    - `tag`
    - `linux_distribution` (`deferred|deb|appimage|rpm|tarball`)
  - Builds macOS and Windows installers, runs smoke checks, and publishes assets
    to the corresponding GitHub Release.
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

## Standard Tagged Release

1. Ensure `main` is green in CI.
2. Create and push a version tag:
   - `git tag vX.Y.Z`
   - `git push origin vX.Y.Z`
3. Wait for `Release Installers` workflow completion.
4. Verify GitHub Release contains both installer assets.

## Manual Release Re-run

Use Actions → `Release Installers` → `Run workflow` and provide:

- `tag`: an existing tag (for example `v0.1.0`)
- `linux_distribution`: planned Linux channel for this release metadata

This rebuilds installers and updates assets on that tag’s release.

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
