# Release Notes: `v0.1.0-internal.3` (draft)

This draft internal release focuses on making the DNA-window workflows feel
more deliberate, clarifying RNA/splicing analysis entry points, tightening
cloning-oriented restriction-site ergonomics, adding promoter-aware genome-gene
retrieval, and carrying the release story forward with a GitHub-downloadable
container image alongside the desktop artifacts.

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
  - the dedicated `RNA-read Mapping` window now exposes its main mapping
    parameters directly instead of hiding them behind a collapsed
    Nanopore-specific section
  - former scope wording was clarified to `Region / gene scope`, with more
    explicit descriptions of broad vs target-focused transcript selection
  - saved-report alignment detail now reconstructs the exact `rust-bio`
    pairwise alignment, including explicit reverse-complement query
    orientation when that is how the read fits the transcript template
- Promoter-oriented genome retrieval is now first-class:
  - `Extract Selected Gene` can retrieve either the full gene span or
    `CDS + promoter`
  - users can request an explicit upstream promoter length without manually
    calculating coordinates around the coding region
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
- Tutorial/release metadata hardening moved closer to release shape:
  - committed tutorial catalog/manifest/generated outputs are checked more
    strictly against their source units in CI
  - this reduces the chance of shipping a tag whose downloadable tutorial/docs
    assets are out of sync with the committed generators
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
- That dedicated workspace now behaves more like a true specialist surface:
  - primary run controls are visible immediately on first open
  - the window is no longer visually framed as a mostly hidden
    “Nanopore cDNA interpretation” foldout
  - the action language is less tied to one sequencing vendor profile even
    though `nanopore_cdna_v1` remains the currently built-in profile
- Saved-report aligned-read inspection became more trustworthy and easier to
  interpret:
  - selected rows now reconstruct the exact `rust-bio` pairwise alignment
    strings for the retained phase-2 hit
  - the alignment midline explicitly means:
    `|` exact match, `.` mismatch, blank indel/gap
  - reverse-complemented read orientation is surfaced explicitly instead of
    leaving complementary-strand reads looking like a confusing wrong-way
    alignment
- Score-density/read-retention behavior was also hardened:
  - post-run histogram inspection can focus exact retained saved-report rows
    from a chosen score bin instead of only a capped live preview
  - `composite gate` and retained-only `retained replay` populations make it
    clearer whether the histogram shows all scored rows, recorded gate-passing
    rows, or a fast replay under current controls
  - retained saved reports now keep obvious high-score outliers through rescue
    rules beyond the baseline top-rank cap, so later inspection/alignment is
    less likely to lose them

### 3) Genome Retrieval and Promoter-Focused Extraction

- Gene retrieval now better matches cloning-oriented questions:
  - GUI `Extract Selected Gene` supports either the full gene span or
    `CDS + promoter`
  - an explicit upstream promoter length can be entered directly in base pairs
- This uses one shared engine contract (`extract_mode=gene|coding_with_promoter`
  plus `promoter_upstream_bp`) rather than a GUI-only coordinate helper.
- The intended use case is promoter/reporter and promoter-plus-coding-region
  work where users should not have to hand-calculate strand-aware upstream
  spans around CDS coordinates.

### 4) Restriction / Cloning Ergonomics

- Restriction-site display modes now behave more like a bench-facing filter
  choice than a recomputation toggle.
- Preferred-enzyme handling has a clearer baseline:
  - pUC MCS defaults
  - Golden Gate Type IIS defaults
  - manual comma-separated preferred list
- The Golden Gate preset is drawn from the active REBASE-derived catalog rather
  than being hard-coded only as a UI label.

### 5) Prepared Cache Cleanup and Rebuild Handoff

- Cache cleanup now keys prepared installs by path rather than only by
  `entry_id`, avoiding accidental co-selection when the same prepared genome id
  exists under more than one selected cache root.
- Shared shell/CLI parity now exposes the same exact-targeting concept via
  `--prepared-path`.
- The GUI cleanup flow keeps its rebuild handoff behavior after partial cleanup
  while clearing selection state more deterministically.

### 6) Release / Distribution Story

- In addition to macOS `.dmg` and Windows `.zip` release artifacts, release tags
  now publish a GitHub-downloadable GHCR container image.
- The image is Debian-first and multi-arch (`linux/amd64`, `linux/arm64`).
- It is intended as the shared Linux/macOS container route for:
  - GUI via browser/noVNC
  - CLI
  - MCP server
  - embedded JavaScript and Lua shells
  - thin Python wrapper access
- Tutorial/documentation artifacts are now more explicitly part of the release
  hygiene story:
  - tutorial chapter generation/checking is treated as a release-facing CI gate
  - tutorial source ordering and committed generated outputs are checked for
    deterministic consistency before a tag can be considered clean

## Interim Release Readiness

- The current work continues the shift toward explicit, shared-engine behavior:
  less automatic UI magic, clearer handoff points, and more stable cross-adapter
  semantics.
- This draft now also tells a more complete release-facing story for promoter
  retrieval and RNA-read interpretation:
  - promoter-plus-coding extraction is now a real shared feature, not a manual
    coordinate workaround
  - RNA-read aligned-hit review now explains what the alignment glyphs mean and
    when a complementary-strand read was reverse-complemented for scoring
- The new container distribution path strengthens the “download and run a real
  GENtle environment” story even while native Linux installer packaging remains
  deferred.

## Release-Facing Known Limitations

- Multi-insert Gibson execution still currently requires a defined destination
  opening; `existing_termini` remains the single-fragment handoff path.
- RNA-read `retained replay` is intentionally a retained-report replay under
  current controls; it does not revisit reads that were never retained by the
  original run.
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
