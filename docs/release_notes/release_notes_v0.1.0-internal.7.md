# Release Notes / Changelog: `v0.1.0-internal.7` (draft)

This draft internal release is the **pre-refactor stability cut**. Compared
with `v0.1.0-internal.6`, this cut does not introduce a new biology theme on
its own — instead it consolidates the regulatory, RNA-mapping, qPCR/cDNA, and
ClawBio surfaces that landed during the v6 cycle, hardens them with broader
test coverage and CI fixes, and tidies the repository so the upcoming
`gentle-engine` workspace extraction can run against a clean baseline.

Tagging now captures the green-CI snapshot that has the full Phase 1
(`gentle-protocol`) extraction in place, plus several follow-on slices in
`gentle-render`, `gentle-shell`, and `gentle-gui`, while keeping the root
crate's biology execution behavior untouched. Subsequent tags are expected to
begin moving execution code out of the root crate.

## Highlights

- RNA-read mapping became a routinely usable cohort/preflight surface:
  - paralog negative controls (TP53/TP63) and positive controls (TP73) are now
    integrated into the deterministic threshold/preflight path
  - gene-group catalog records (new in this cut) link directly into RNA-read
    mapping so cohort comparisons can use named families
  - cohort/batch automation, status outputs, family-plot graphics, and
    `--seed-only` defaults consolidated around the gene-agnostic pancreas
    screen helper script
- ClawBio-facing surface has its largest single-cut expansion since the skill
  scaffold was introduced:
  - qPCR/TaqMan, transcript-derived cDNA PCR/qPCR assay testing, gel-ready
    cDNA product visualisation, and lab-assistant handoff exports all gained
    typed ClawBio modes
  - exon-skip planning/materialisation is now a two-phase ClawBio workflow with
    return-channel parity (`exon-skip-plan` then `exon-skip-materialize`)
  - command-palette UI-intent entries are now rebuilt from shared metadata and
    handed off through the ClawBio `capabilities` surface so chat planners and
    GUI surface the same names
  - PCR protocol cartoons are now reachable from ClawBio
- Helper-construct vocabulary is now an engine-owned, catalog-extensible
  semantics layer rather than ad-hoc strings:
  - normalized helper-construct semantics, catalog-extensible vocabulary
    overlays, doctor command, and guide fallback/artifact continuation
- TFBS / promoter / regulatory reasoning continued to deepen:
  - first promoter evidence matrix (draft), reasoned promoter-region
    definition, multi-gene promoter window comparisons, and collapsed
    identical-promoter views by transcript count
  - sequence-logo presentation fixes and JASPAR identifier CI hardening
- Repeat features became a real shared track:
  - `gentle-protocol` parsing and display contracts for UCSC RepeatMasker
    (`rmsk`), resource-prep and integration into the joint API, RMSK record
    tightening, and lightweight GUI repeat-details
- Coordinate-system clarity:
  - explicit mapping between genome, transcriptome, and proteome coordinates
  - a clean separation between 0- and 1-based positions across the engine
    surface
  - transcript-position graphics for cDNA PCR/qPCR assay testing
- Dotplot inspection gained genome-context awareness:
  - genome-context rail, exon side-block annotation, genomic carry-over
- Primer3 hardening:
  - first pre-release Primer3 hardening pass, Primer3 provenance on oligo QC
    reports, Tm recalibration, self-hybridisation tests in qPCR, and a local
    Primer-BLAST-style specificity check
- Repository tidies that prepare for the workspace split:
  - obsolete uncompiled source tree removed
  - stray root `lib.rs` and duplicate root generated SVGs removed
  - release notes consolidated under `docs/release_notes/`
  - app/main-area module decomposition continued, shell-first modularisation
    advanced
- First external-service-provider implementation:
  - Metabion as the initial provider through the provider-neutral handoff
    contract; GUI layer for external services; service-provider handoff
    reorganisation; tutorial for the Metabion handoff route

## Notable Changes by Area

### 1) RNA-Read Mapping, Gene Groups, and Cohort Screens

- TP53/TP63 are now first-class **paralog negative controls** for RNA-read
  mapping; TP73 remains the positive control. Threshold/preflight tests
  consume both groups deterministically.
- New deterministic **gene-group layer** with catalog-extensible records and a
  link directly into RNA-read mapping so cohorts and panels can be described
  by group rather than ad-hoc symbol lists.
- **External knowledge wiring**: GO is now an external source feeding the
  internal gene-group catalog; agent-driven gene catalogs landed; a small
  helper script `scripts/fetch_ensembl_cdna_fixtures.sh` retrieves cDNA
  fixtures for paralog controls.
- `--seed-only` is now the **default for cohort runs**; phase-2 alignment is
  explicit opt-in via `--with-alignment --align-selection seed_passed`.
- RNA-mapping presentation:
  - 90th percentile is also shown for cohort plots
  - q90/q100 derivation rewritten for the new plot abstraction
  - family-picture plots rebuilt around a shared plotting abstraction
  - status outputs for batch/cohort runs are now stable across list/show
- Auto-download path for missing gene sequences is corrected for the small-set
  edge case.
- Concatemer review, cDNA vs direct-RNA reporting, and read-report
  list/show/export are now explicit and stable across CLI surfaces.

### 2) ClawBio Skill Surface, UI Intents, and Lab-Assistant Handoff

- **Two-phase exon-skip isoform workflow**: `exon-skip-plan` produces a stored
  plan; `exon-skip-materialize` consumes a `plan_id` with `confirm=true`. The
  return-channel parity contract makes the materialised product inspectable
  on the same ClawBio response shape used by the plan step.
- New ClawBio modes / parity slices for:
  - qPCR / TaqMan design and direct assay testing
  - transcript-derived cDNA PCR/qPCR assay testing, including gel-ready
    materialisation (idempotent) of non-specific products into a vial/container
    and product-pool gel rendering
  - lab-assistant handoff export (assistant-ready wet-lab cloning instructions
    from current design/run history)
  - PCR protocol cartoon
  - direct Ensembl ROI fetching
- **UI-intent discoverability**: command-palette entries are now rebuilt from
  shared metadata; that same metadata is merged into the ClawBio
  `capabilities` payload so chat planners and the GUI palette agree on action
  names.
- ClawBio polish + skill-description revisions; a "GENtle classic vs local
  rewrite" disambiguation pass clarifies which runtime the wrapper is talking
  to (`mode: "version"` / `request_runtime_version.json`).
- Gibson tutorial reworked into a complete ClawBio demo path.
- Several `fix: ClawBio …` entries removed redundant continuation steps and
  fixed the suggested-action follow-up in demos.

### 3) Helper-Construct Vocabulary and Catalog Semantics

- The **helper-construct semantics vocabulary** is now an engine-owned,
  catalog-extensible surface rather than ad-hoc strings. New layers:
  - normalized helper-construct semantics
  - catalog-extensible overlay support
  - `helper-construct-vocab doctor` command for inspecting resolved records,
    overlay provenance, duplicate canonical terms, alias collisions, and
    routine-hint/routine-family diagnostics
  - guide fallback and artifact continuation for GENtle-side answers when a
    catalog entry is missing

### 4) TFBS, JASPAR, and Promoter Reasoning

- Draft **promoter evidence matrix** record.
- First reasoned **promoter-region definition** linking back to the design
  window for downstream planning.
- Multi-gene promoter window comparisons; identical promoter windows are now
  collapsed by transcript count to reduce visual noise.
- Sequence-logo presentation fixes; CI JASPAR identifier hardening (resilient
  to upstream identifier surface drift).

### 5) Repeats (RMSK) as a Shared Track

- Engine/protocol-level parsing for UCSC RepeatMasker (`rmsk`) — moves repeat
  features off ad-hoc strings onto typed records.
- Resource-prep contracts for `rmsk` tables; rmsk repeat-display preparation;
  joint-API integration so other tracks can co-reference repeat context.
- Lightweight GUI repeat-details panel; tightened RMSK record formatting and
  display hardening for malformed lines.

### 6) Coordinate Mapping and cDNA Transcript Context

- Explicit **genome ↔ transcriptome ↔ proteome coordinate mapping** records
  and helpers.
- Clear distinction between **0- and 1-based** sequence positions across the
  engine surface; affected reports and request schemas now state the basis.
- **Transcript-position graphics** for cDNA PCR/qPCR assay testing; exon
  architecture and primer-legend refinements in cDNA transcript-map SVGs.
- New `feature: emit more TaqMan-like qPCR requests` rounds out the qPCR/TaqMan
  emission path for typed automation.
- Genomic-aligned cDNA assay maps so common primer loci, exon identity, and
  spliced-source gaps line up across isoforms when ClawBio needs comparison.

### 7) Dotplot Genome Context

- **Genome-context rail** on dotplots.
- **Exon side-block annotation** for dotplots.
- Genomic carry-over awareness so context is preserved across nearby loci
  rather than dropping at panel edges.

### 8) Primer3 Hardening and Specificity

- First pre-release Primer3 hardening pass.
- **Primer3 provenance** is recorded on oligo QC reports.
- Tm recalibration.
- qPCR self-hybridisation tests added to the design QC path.
- **Local Primer-BLAST-style specificity** check landed for primer/qPCR design.
- New `tests: cover qPCR specificity fallback`.

### 9) Restriction Sites, Sequence Inspection, and Selection

- **SnapGene-style restriction-site inspector** in the GUI.
- First **restriction-site tooltip/inspector** with map-level context menu.
- Same restriction-enzyme detail surface is now also reachable through CLI and
  MCP routes (parity between GUI hover/popover and headless automation).
- **Formula-based selection** in the GUI (Enter-key interaction covered by
  tests); shared-shell / CLI formula parity.
- Tooltip copy → clipboard, including names/descriptions in the DNA sequence
  viewer.

### 10) Visual SVG Guardrails

- Machine-readable **SVG role markers** for visual benchmark fixtures
  (`data-gentle-role`, `data-gentle-feature-kind`).
- Visual benchmark fixtures for dense SVG exports.
- **Visual SVG lint guardrail** and quantitative SVG readability checks for
  marked labels.

### 11) GenBank / XML Import Parity

- **INSDSet / INSDSeq import** support.
- XML parity work alongside a compact multi-record GBSeq fixture for
  cross-format sequence parity tests.

### 12) External Service Providers and Order-Form Concept

- First **provider implementation** (Metabion) on the provider-neutral
  service-request/status/artifact contract.
- First **GUI layer for external services** for inspecting handoff state.
- Service-provider handoff path reorganised; tutorial for the Metabion
  external-service handoff (GUI + CLI).
- **Order-form concept** sketched (`planning: describing link to order forms`,
  `planning: GeneArt roadmap.md`) — provider-neutral order-batch artifacts
  with line-item provenance back to primer/qPCR reports.

### 13) Microarray, Publication Data, and Genome Projection

- Publication-resource downloads are now reproducible across runs.
- Generic **manifest-driven batch helper** for bulk acquisitions.
- **Genome-assembly projection** support (cross-build coordinate projection).
- Microarray + publication data wired into the broader regulatory grouping
  coverage.

### 14) SRA Toolkit Integration

- **SRA toolkit wrapper** as an executable resource.
- SRA toolkit user-handling improvements (credential / cache handling, not
  the credential body itself).

### 15) Agent Bridge and Safety

- Agent Assistant live-setup path.
- **LLM safety hardening** pass (general guardrails on agent-invoked
  commands).
- JS and Lua agent safety review and tightening.

### 16) Architecture, Workspace Split, and Internal Engineering

- **Continued shell-first modularisation**: more shell commands moved off the
  large root dispatcher onto crate-owned executors, executor remains
  stack-safe for nested workflow/macro commands.
- **GUI module-decomposition** for CUT&RUN; CUT&RUN regulatory-support GUI
  inspection slice now lives in its own module.
- `src/app.rs` reorganisation continued.
- `gentle_cli` source refactoring: filtering, candidates, pool helpers,
  rendering/export, and protocol-cartoon CLI families split into smaller
  modules.
- CLI adapter split: small legacy services extracted; shared `cli_support`
  module landed.
- **Two-way CLI documentation guard**: glossary <-> CLI command surface stays
  in sync.
- Genome-prepare lifecycle hardening (stale-lease / heartbeat refinements).
- Multiple `fix: stack overflow` recurrences in the shell executor are
  addressed; the executor continues to run on an explicitly sized worker
  stack with recursion guard.

### 17) Repository Tidies (pre-tag cleanups)

- **Obsolete uncompiled source tree** under `src/obsolete/` removed.
- **Stray root `lib.rs`** removed (was a stale orphan duplicate of the
  `gentle-protocol` crate root and was not referenced by any Cargo target).
- **Duplicate root generated SVGs** removed (TP73/E2F1 UniProt projections
  now exist only in their proper locations).
- Internal release notes moved under `docs/release_notes/`.

### 18) GUI Reliability

- **GUI theme** layer implemented.
- DNA sequence viewer / RNA mapping **window titles** stabilised, including
  splicing windows and context-menu launched windows.
- **Window stacking order** fix and "drag moves inner frame" fix.
- Saved-changes indicator visibility fixes.
- Hardened SVG → PNG rasterization (font availability path; relative path
  resolution).
- Consistent **window-open behavior** across launch points.

### 19) CI and Test Robustness

- Multiple `fix: CI stack overflow` rounds for the shell-test surface, all
  resolved.
- Machine-independence fixes for tutorial generation under CI.
- **Reduced-precision float comparisons** for CI compatibility across hosts.
- CI test URL specification race resolved.
- **Daily Python bug-scan** for `--help` coverage of helper scripts.
- New scrolling/egui regression tests; visual benchmark tests; agent-safety
  tests for JS/Lua wrappers.

## Release-Facing Known Limitations

- The `gentle-engine` workspace member exists as a 5-line placeholder.
  Execution code still lives in the root crate. The next planned post-tag
  work moves the first execution-side modules into `gentle-engine`, starting
  with the file-loader seam currently still attached to `GENtleApp`.
- CI currently runs on `macos-latest` only. A `ubuntu-latest` parity job is
  the recommended post-tag follow-up given the Debian-first container
  policy in `docs/architecture.md`.
- Twelve top-level dependencies in `Cargo.toml` are still on wildcard
  versions (`version = "*"`). A pre-extraction PR should consolidate these
  into `[workspace.dependencies]` with explicit pins captured from the
  current resolved `Cargo.lock`.
- The ClawBio wrapper's `SUPPORTED_REQUEST_MODES` list has grown to ~40
  entries; a planned post-tag audit will classify each as thin pass-through
  vs deprecation candidate.
- No GUI parity for several newer CLI/MCP routes (helper-construct doctor,
  some external-service provider flows). These remain CLI/MCP-first by
  design for this release.

## Install / Package Notes

- macOS release artifact remains `.dmg`.
- Windows release artifact remains `.zip` containing `gentle.exe`.
- Linux native installer packaging remains deferred; release metadata
  continues to default Linux distribution intent to `tarball`.
- Container/runtime story unchanged:
  - `ghcr.io/smoe/gentle_rs:cli` for headless/agent/ClawBio use
  - `ghcr.io/smoe/gentle_rs:gui` or `:latest` for interactive GUI container use
- ViennaRNA dependency in the Docker image now correctly resolves from the
  Debian non-free channel.

## Suggested Pre-Tag Smoke Checklist

At minimum before tagging `v0.1.0-internal.7`, run:

```bash
cargo check -q
cargo check -q --features js-interface
cargo check -q --features lua-interface
cargo check --all-features -q
cargo test -q workflow_examples -- --test-threads=1
cargo test -q rna_reads
cargo test -q "engine::tests::"
cargo run --bin gentle_examples_docs -- --check
cargo run --bin gentle_examples_docs -- tutorial-check
cargo run --release --bin gentle -- --version
cargo run --release --bin gentle_cli -- capabilities
```

Because this release touched RNA mapping, qPCR/cDNA assay testing,
ClawBio-facing modes, helper-construct semantics, restriction-site UI, and
the gene-group layer, one short manual smoke pass should also cover:

- launching the GUI, opening a known TP73 sequence/project
- opening the RNA-read Mapping workspace and inspecting a saved report
- exporting at least one evidence artifact (TSV, SVG, or evidence bundle)
- running one `gentle_cli` ClawBio request example end-to-end
  (`integrations/clawbio/skills/gentle-cloning/examples/request_workflow_simple_pcr_primer_design_offline.json`
  is the smallest safe choice)
- running one transcript-aware qPCR design from the GUI
- opening the new SnapGene-style restriction-site inspector and confirming
  the same detail is reachable from `gentle_cli`

## Post-Tag Refactor Direction

The tag is intentionally cut before larger reorganisation work. The next
planned PRs against `v0.1.0-internal.7` are:

1. Repo-root tidies were completed pre-tag; the equivalent post-tag task is
   to keep `lib.rs`, root SVGs, and `src/obsolete/` from reappearing
   (CI guard).
2. Extract `GENtleApp::load_from_file` + the five `load_dna_from_*_file`
   helpers off the GUI top-level type so the engine stops depending on
   `crate::app::GENtleApp`. This unblocks `gentle-engine` extraction.
3. Seed `gentle-engine` with its first real module (recommended:
   `dna_sequence` or `feature_location`) once (2) lands.
4. Consolidate root dependencies into `[workspace.dependencies]` with
   explicit pins.
5. Add a `ubuntu-latest` parity job to `.github/workflows/ci.yml`.
6. ClawBio thin-wrapper consolidation: runtime `intents` capability,
   parity tests between `INTENTS.json` / `examples/*.json` / `SKILL.md`
   keyword list, and removal of the literal `"clawbio"` default in
   `services handoff --scope`.

Each of these is intended to be one independent PR so the tag-to-tag diff
remains reviewable.

## Release-Shaped Smoke Check Results

*Pending — to be populated when the pre-tag smoke runs are executed.*

Local release-prep checks expected on the tag date:

- `git diff --check`: _pending_
- `cargo check -q`: _pending_
- `cargo test -q workflow_examples -- --test-threads=1`: _pending_
- `cargo test -q rna_reads`: _pending_
- `cargo build --release --features script-interfaces`: _pending_
  - run with `CARGO_TARGET_DIR=/tmp/gentle_release_smoke`
  - run with `CARGO_INCREMENTAL=0`
- release-artifact checks from the same built target:
  - `gentle --version`: _pending_
  - `gentle_cli capabilities`: _pending_
  - `gentle_js --version`: _pending_
  - `gentle_lua --version`: _pending_
  - `gentle_examples_docs --check`: _pending_
  - `gentle_examples_docs tutorial-check`: _pending_
  - `gentle_mcp --help`: _pending_

Manual GUI smoke pass: _pending_.
