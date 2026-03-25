# Release Notes: `v0.1.0-internal.2`

Release baseline: changes since tag `v0.1.0-internal.1`.

This internal milestone cut focuses on getting the current GUI/CLI/tutorial
surface into a release-shaped, demonstrable state rather than broad new
feature expansion.

## Highlights

- Gibson specialist baseline is now internally releasable:
  - destination-first planning and preview for one or more ordered inserts
  - Gibson-specific overlap/Tm controls, primer suggestions, review blocks,
    and cartoon export from one shared preview contract
  - apply creates a single lineage-visible `Gibson cloning` operation with
    primer/product outputs that can reopen the specialist for review
- GUI lineage export now matches what users actually see:
  - `File -> Export DALG SVG...` and graph-canvas `Save Graph as SVG...`
    export the same grouped / Gibson-hub graph model shown in the main window
  - save/export dialogs now suggest filenames derived from the current project
    name rather than generic placeholders
- Tutorial/help copy was tightened for the interim release:
  - Gibson landing/tutorial pages now call out the current multi-insert apply
    guardrail explicitly
  - release docs now include a release-shaped pre-tag smoke checklist using
    `--features script-interfaces`

## Known Limitations

- Multi-insert Gibson execution currently requires a defined destination
  opening; `existing_termini` remains the single-fragment handoff path.
- The manual Gibson tutorial remains intentionally centered on the shipped
  single-insert apply path, with ordered multi-insert preview coverage called
  out as a guardrailed capability rather than a broader tutorial flow.
- Linux installer packaging remains deferred; release metadata still defaults
  the Linux distribution channel to `tarball`.

## Smoke Check Results

Local pre-tag matrix run on 2026-03-25:

- `cargo build --release --features script-interfaces`: pass
- `cargo run --release --bin gentle -- --version`: pass
  - reported `GENtle 0.1.0-internal.2`
- `cargo run --release --bin gentle_cli -- capabilities`: pass
- `cargo run --release --features js-interface --bin gentle_js -- --version`: pass
- `cargo run --release --features lua-interface --bin gentle_lua -- --version`: pass
- `cargo run --release --bin gentle_examples_docs -- --check`: pass
  - reported `example_count = 17`
- `cargo run --release --bin gentle_examples_docs -- tutorial-check`: pass
  - reported `chapter_count = 15`
- `cargo run --release --bin gentle_mcp -- --help`: pass

Execution note:

- after the successful release-shaped build, the entrypoint probes for
  `gentle_cli`, `gentle_js`, `gentle_lua`, `gentle_examples_docs`, and
  `gentle_mcp` were run against the freshly built `target/release/*` binaries
  to avoid redundant relinks while validating the same release artifacts

## Install / Package Notes

- verified in `.github/workflows/release.yml`:
  - macOS release artifact: `.dmg`
  - Windows release artifact: `.zip` containing `gentle.exe`
  - Linux release metadata default: `tarball`

## Pre-Tag Checklist

- Run the local release-shaped smoke matrix from [`docs/release.md`](docs/release.md).
- Confirm the Gibson specialist help/status text still warns that multi-insert
  apply requires a defined destination opening.
- Verify README/tutorial wording still matches current menu entries:
  - `Patterns -> Gibson...`
  - `File -> Export DALG SVG...`
  - graph-canvas `Save Graph as SVG...`
