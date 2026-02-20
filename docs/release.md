# Release Process

This project publishes installable desktop packages for:

- macOS: `.dmg`
- Windows: `.msi`

Linux installable packaging is intentionally deferred. The target is Debian
packaging once dependency packaging is finalized.

## Workflows

- CI validation workflow: `.github/workflows/ci.yml`
  - Runs checks/tests on pushes to `main` and pull requests.
  - Does not publish release assets.
- Release workflow: `.github/workflows/release.yml`
  - Triggered by tag pushes matching `v*`.
  - Can also be run manually via `workflow_dispatch` with a `tag` input.
  - Builds macOS and Windows installers, runs smoke checks, and publishes assets
    to the corresponding GitHub Release.

## Artifact Naming

Release assets are normalized to:

- `gentle-<tag>-macos-<arch>.dmg`
- `gentle-<tag>-windows-<arch>.msi`

Example:

- `gentle-v0.1.0-macos-arm64.dmg`
- `gentle-v0.1.0-windows-x64.msi`

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

This rebuilds installers and updates assets on that tag’s release.

## Smoke Checks in Release Workflow

- macOS:
  - mounts the `.dmg`
  - verifies an `.app` bundle exists with `Contents/Info.plist`
- Windows:
  - performs MSI administrative extract (`msiexec /a`)
  - verifies `gentle.exe` is present and non-empty

## Rollback / Recovery

If a release artifact is broken:

1. Fix on `main`.
2. Re-run release workflow for the same tag via manual dispatch, or
   create a new patch tag (`vX.Y.Z+1`) and publish normally.
3. If needed, remove incorrect assets from the GitHub Release UI before
   republishing.

