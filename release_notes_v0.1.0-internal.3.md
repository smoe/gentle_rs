# Release Notes: `v0.1.0-internal.3` (draft)

This draft internal release focuses on making the DNA-window workflows feel
more deliberate, clarifying RNA/splicing analysis entry points, tightening
cloning-oriented restriction-site ergonomics, and carrying the release story
forward with a GitHub-downloadable container image alongside the desktop
artifacts.

These notes intentionally follow the same “show both what is already stable and
what is clearly being prepared” style as the previous internal cut, so readers
can see where GENtle is becoming more usable for biologists even before every
path has the same end-to-end validation depth.

## Highlights

- DNA-window UX moved toward more explicit intent:
  - Splicing no longer auto-opens on single-click selection.
  - The Splicing Expert now opens deliberately by double-click, context-menu
    action, or the description-pane button.
  - Painted PCR intervals can now be edited numerically (`start..end`) after a
    drag instead of forcing a repaint.
- RNA/splicing workflows became easier to understand and reach:
  - the dedicated cDNA / RNA-read mapping workspace is easier to open from the
    DNA window
  - former scope wording was clarified to `Region / gene scope`, with more
    explicit descriptions of broad vs target-focused transcript selection
- Restriction/cloning ergonomics improved:
  - restriction-site display is now decoupled from the underlying site
    computation
  - preferred enzyme lists are preserved as first-class user choices
  - a Golden Gate Type IIS preset is now available from the DNA window and
    graphics configuration
- Prepared-cache cleanup is now path-based end to end:
  - duplicate prepared ids across multiple cache roots can be inspected,
    previewed, and cleaned independently
  - shared shell/CLI parity now includes `--prepared-path`
- Release-tag container distribution is now part of the story:
  - GitHub publishes a Debian-first multi-arch GHCR image for release tags
  - the same image carries GUI, CLI, MCP, JS, Lua, and Python-entrypoint
    workflows
- Deferred-window placement behavior was hardened:
  - one-shot viewport positioning is now sent explicitly for deferred/embedded
    windows, which should reduce the recent “maximisation” / misplaced-window
    drift

## Notable Changes by Area

### 1) DNA Window, Splicing, and PCR UX

- Feature selection and expert opening behavior are now separated:
  - select/focus stays lightweight on single-click
  - deeper splicing inspection requires explicit user intent
- The Splicing Expert is now reachable from:
  - feature-tree double-click
  - DNA-map double-click
  - feature/map context menus (`Open Splicing Window`)
  - description-pane action button
- Painted PCR intervals now keep typed start/end editor fields in the post-drag
  chip, with an explicit apply action.

### 2) RNA-Read Mapping Clarity

- RNA-read mapping language was tightened around the actual decision being made:
  - `Region / gene scope` now describes how broad the transcript-template pool
    becomes before seeding/alignment
  - target-group vs all-overlapping wording is more explicit about strand and
    competition tradeoffs
- The DNA window now offers more direct access into the dedicated cDNA
  read-mapping workspace.

### 3) Restriction / Cloning Ergonomics

- Restriction-site display modes now behave more like a bench-facing filter
  choice than a recomputation toggle.
- Preferred-enzyme handling has a clearer baseline:
  - pUC MCS defaults
  - Golden Gate Type IIS defaults
  - manual comma-separated preferred list
- The Golden Gate preset is drawn from the active REBASE-derived catalog rather
  than being hard-coded only as a UI label.

### 4) Prepared Cache Cleanup and Rebuild Handoff

- Cache cleanup now keys prepared installs by path rather than only by
  `entry_id`, avoiding accidental co-selection when the same prepared genome id
  exists under more than one selected cache root.
- Shared shell/CLI parity now exposes the same exact-targeting concept via
  `--prepared-path`.
- The GUI cleanup flow keeps its rebuild handoff behavior after partial cleanup
  while clearing selection state more deterministically.

### 5) Release / Distribution Story

- In addition to macOS `.dmg` and Windows `.zip` release artifacts, release tags
  now publish a GitHub-downloadable GHCR container image.
- The image is Debian-first and multi-arch (`linux/amd64`, `linux/arm64`).
- It is intended as the shared Linux/macOS container route for:
  - GUI via browser/noVNC
  - CLI
  - MCP server
  - embedded JavaScript and Lua shells
  - thin Python wrapper access

## Interim Release Readiness

- The current work continues the shift toward explicit, shared-engine behavior:
  less automatic UI magic, clearer handoff points, and more stable cross-adapter
  semantics.
- The new container distribution path strengthens the “download and run a real
  GENtle environment” story even while native Linux installer packaging remains
  deferred.

## Release-Facing Known Limitations

- Multi-insert Gibson execution still currently requires a defined destination
  opening; `existing_termini` remains the single-fragment handoff path.
- Linux native desktop installer packaging remains deferred; Linux release
  metadata still defaults to `tarball` even though the GHCR container image is
  now a real release-tag distribution route.
- Several areas noted above remain stronger in the GUI and shared shell than in
  polished tutorial coverage; this draft should be refined again once the next
  actual cut is smoke-tested.

## Install / Package Notes

- macOS release artifact remains `.dmg`
- Windows release artifact remains `.zip` containing `gentle.exe`
- release tags also publish a GitHub-downloadable GHCR image:
  - `ghcr.io/<owner>/<repo>:<tag>`
  - `latest` moves on release-tag publishes
  - current published platforms: `linux/amd64`, `linux/arm64`
- Linux release metadata still defaults to `tarball`

## Release-Shaped Smoke Check Results

Pending for this draft. Before tagging, run the local release-shaped matrix from
[`docs/release.md`](docs/release.md) and record pass/fail here.
