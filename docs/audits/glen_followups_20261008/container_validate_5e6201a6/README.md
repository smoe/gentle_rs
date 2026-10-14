# Internal.12 Non-Publishing Container Evidence

This supplements gate 6; it does not replace any historical report or accept
newer main, `.13`/opt2, native installers, physical GUI or scientific results.
Independently retrieved and checked on 2026-10-14 without building or pulling
an image, executing a model, pushing or dispatching a workflow.

## Origin And Exact Bytes

`container-build.json` is the original 2,425-byte CI receipt from
[tag-push run 38322028567](https://github.com/smoe/gentle_rs/actions/runs/38322028567),
attempt 1, exact source/workflow
`5e6201a69a599cba87491c86cfd7311cb6fee173`. The
[container job 115052172732](https://github.com/smoe/gentle_rs/actions/runs/38322028567/job/115052172732)
completed successfully; the GHCR publication job was skipped. This is distinct
from the later release-event publication run `38330177090`.

Artifact `11781064566`, named
`container-evidence-5e6201a69a599cba87491c86cfd7311cb6fee173`, contains only this
receipt. The 1,154-byte ZIP independently matches GitHub's artifact digest
`0128405dd13a7e5fa85a0c8262f0ac7fd22364e466dcce9e5a5ede22c394be9d`.
The extracted receipt SHA-256 is
`682d7d10df1b907dd80121c7055c62b3335952147fc7a3bf41988e7f13ffc42d`.
No fields, spacing or line endings were changed. Its root lockfile digest
matches `git show 5e6201a6:Cargo.lock`; the recipe matches that source's
`[profile.package-opt1]`, not the current default package profile.

Recreate into a fresh evidence directory, verify the ZIP digest/membership
before extracting the single named file, and compare its raw digest:

```sh
gh api --allow-escape-sequences repos/smoe/gentle_rs/actions/artifacts/11781064566/zip > container-evidence.zip
shasum -a 256 container-evidence.zip
unzip -l container-evidence.zip
unzip -p container-evidence.zip container-build.json > container-build.json
shasum -a 256 container-build.json
```

## Verified Scope And Limits

The actual `linux/amd64` `runtime-cli` build uses `package-opt1`, no desktop
features, `mode=validate_only` and `pushed=false`. Its workflow executes the
default entrypoint and all three help routes with `--network none`, checks
runtime libraries and desktop-binary exclusion, imports the Python adapter,
and exercises the synthetic RNAPKIN SVG/PNG helper smoke. CI executes each
packaged binary's `--version` and computes its SHA-256 inside the loaded image.
All three retained identities bind to `.12` and the exact source above.

The receipt is independently hash-checked here. Image and executable digests
remain CI-observed, not independently rehashed from a locally pulled image.
This record does not establish native archive membership/extraction, downloaded
installer integrity, native readability, a timed performance budget, genuine
live-agent behavior or final-source six-gate completion. Independent native
package acceptance and new-source/profile checks remain open.

GENtle consumes this as audit/provenance evidence, not a biological runtime
input. `scripts.test_package_desktop` checks its frozen bytes and declared
scope. The fast `scripts.test_tutorial_checkouts` gate tests LF/CRLF preservation
and refusal of a missing scoped `.gitattributes` rule. No existing frozen
evidence, `paper` or `output/` was changed.
