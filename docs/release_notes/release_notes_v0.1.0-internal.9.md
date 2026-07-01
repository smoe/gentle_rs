# Release Notes / Changelog: `v0.1.0-internal.9` (draft)

This draft internal release is a focused follow-up to
`v0.1.0-internal.8`. It keeps the biology engine contract stable and
concentrates on four practical areas:

- smoother macOS GUI window behavior by defaulting back to native child
  viewports where they work better than hosted/nested windows;
- opt-in GUI profiling and small DNA-rendering responsiveness fixes;
- tighter ClawBio/OpenClaw action-envelope semantics, including stable
  `skill_id` reporting from the `gentle-cloning` skill;
- clearer rack, plate, SVG/PDF, and OpenSCAD figure exports for presentation
  and bench-layout communication.

The release is still an interim cut. It should be read as "more usable and
better instrumented" rather than as a final GUI architecture.

## Highlights

- **macOS native child viewports are the default again** on
  `egui`/`eframe` `0.34.3`. Manual testing showed native OS windows are much
  smoother than the hosted/nested fallback for moving, resizing, and interacting
  with DNA sequence windows.
- **Hosted/nested windows remain available as an explicit fallback** via
  `GENTLE_MACOS_HOSTED_CHILD_VIEWPORTS=1 cargo run --bin gentle`.
- **Root-window fullscreen restoration is disabled on macOS** because opening
  child windows from a fullscreen/maximized root can trigger unstable
  Split View-style flicker. GENtle now delays new macOS native sequence windows
  until the root has returned to a regular window.
- **GUI latency investigation is now instrumentable** through the optional
  `gui-profiler` feature and `GENTLE_GUI_PROFILE=1`, with coarse Puffin spans
  around app, DNA-window, DNA-map, sequence-text, and feature-tree rendering.
- **ClawBio/OpenClaw response contracts are clearer**: `gentle-cloning` now
  reports `skill_id` and keeps suggested-action routing explicit so downstream
  OpenClaw-style retention layers can invoke GENtle actions without inventing
  a parallel command surface.
- **Rack/plate figure export is richer**: the six-well cell-culture plate
  profile/template, culture-well SVG rendering, arrangement-labelled top-down
  SVGs, `svg-pdf` conversion, README figure assets, and OpenSCAD `.scad`
  sources are now part of the presentation path.
- **The startup splash screen uses the transparent mascot asset** at
  `assets/mascots/Mascot_transparent.png`, with a larger splash presentation.
- **CI parity drift was corrected** by refreshing
  `docs/gui_cli_mcp_parity.md` after protocol/catalog changes.

## Notable Changes by Area

### 1) macOS GUI Window Behavior

- macOS defaults to native OS child viewports on `egui`/`eframe` `0.34.3`.
- `GENTLE_MACOS_HOSTED_CHILD_VIEWPORTS=1` remains the explicit hosted-window
  fallback for regression testing or machines where native viewports misbehave.
- The root GENtle window no longer restores persisted fullscreen/maximized
  state on macOS.
- New macOS native sequence windows are delayed until the root window leaves
  fullscreen/maximized state, avoiding the flickering Split View/tab-like path
  observed when opening a child window from fullscreen.
- Hosted-window ownership and resize handling from `.8` remains in place for
  the fallback path, but hosted mode is no longer the default macOS experience.

### 2) GUI Profiling and Responsiveness

- Added optional Puffin-based GUI frame profiling behind
  `--features gui-profiler`.
- Runtime opt-in uses `GENTLE_GUI_PROFILE=1`; normal builds and normal GUI
  startup remain unchanged.
- Coarse spans cover app update/render, DNA windows, DNA map rendering,
  sequence text rendering, and feature-tree rendering so future latency work
  can be driven from measurements instead of hunches.
- The linear DNA renderer received a small responsiveness pass, scoped to
  reducing unnecessary work without changing sequence-rendering semantics.
- Beta compiler compatibility fixes landed alongside the profiling pass.

### 3) ClawBio/OpenClaw Integration

- The `integrations/clawbio/skills/gentle-cloning` wrapper now exposes a
  stable `skill_id` in its response/catalog metadata.
- ClawBio/OpenClaw suggested-action envelopes were tightened so action ids,
  labels, routing hints, preferred artifacts, and scope remain explicit.
- The `gentle-cloning` skill is now wired toward existing GENtle reporter
  catalog, reporter recommender, corpus export, and reporter construct handoff
  routes, so synthetic-biology follow-up can quote GENtle reports instead of
  inventing prose.
- The boundary remains intentionally deterministic: OpenClaw can retain and
  present suggested actions, but GENtle still owns the actual shell/engine
  route semantics.

### 4) Rack, Plate, and Figure Export

- Added a six-well cell-culture plate rack profile/template and example
  workflow at `docs/examples/workflows/cell_culture_plate.json`.
- Added culture-well SVG rendering and arrangement-labelled top-down
  `racks hero-svg` output.
- Added `svg-pdf` conversion support for figure workflows.
- Added README-facing SVG/PDF assets for the cell-culture plate layout.
- Added OpenSCAD sources for the cell-culture plate and Gibson rack/storage
  figures under `docs/figures/`, keeping printable/presentation geometry
  inspectable as source rather than only as rendered assets.

### 5) Documentation, CI, and Dependency Hygiene

- Refreshed `Cargo.lock` after the local compiler/dependency update pass.
- Updated GUI and CLI docs for native macOS window defaults, profiling
  switches, and rack/plate export routes.
- Regenerated `docs/gui_cli_mcp_parity.md` after ClawBio/protocol catalog
  changes so the parity freshness test passes again.
- Preserved the release packaging story from `.8`: macOS `.dmg`, Windows
  `.zip`, and Linux/container paths unchanged.

## Release-Facing Known Limitations

- macOS native child viewports are currently the preferred default, but opening
  child windows from a fullscreen root remains unstable enough that GENtle
  deliberately avoids restoring fullscreen/maximized root-window state on
  macOS.
- Hosted/nested windows remain useful as a fallback and for comparison, but
  dense overlapping hosted-window setups are still not the target experience.
- The GUI profiler is opt-in developer instrumentation, not a user-facing
  performance panel.
- The `.9` notes do not imply a package-version bump by themselves; update
  `Cargo.toml` separately when the tag is actually cut.
- `proc-macro-error2 v2.0.1` can emit a future-incompatibility warning on
  newer beta toolchains through the `bio -> getset -> proc-macro-error2`
  dependency chain. This is external to GENtle; avoid a local patch unless the
  warning becomes a hard build failure.

## Suggested Pre-Tag Smoke Checklist

At minimum before tagging `v0.1.0-internal.9`, run:

```bash
cargo check -q
cargo test -q --test parity_matrix_freshness
cargo test -q --lib egui_compat --no-fail-fast
cargo run --bin gentle_examples_docs -- --check
cargo run --bin gentle_examples_docs -- tutorial-check
cargo run --release --bin gentle -- --version
cargo run --release --bin gentle_cli -- capabilities
```

Manual GUI smoke:

- On macOS, launch `cargo run --bin gentle` from a regular, non-fullscreen
  root window.
- Open project, DNA sequence, Splicing Expert, and Promoter design windows.
- Confirm native windows move and resize smoothly enough for ordinary use.
- Open a DNA sequence window and confirm DNA selection latency is acceptable.
- Confirm opening a child window from a formerly fullscreen/maximized root does
  not enter the earlier Split View-style flicker path.
- Optionally compare hosted fallback with
  `GENTLE_MACOS_HOSTED_CHILD_VIEWPORTS=1 cargo run --bin gentle`, but treat
  native child viewports as the expected macOS default.

Release-shaped smoke check results: _pending_.
