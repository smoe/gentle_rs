# Release Notes: `v0.1.0-internal.2`

Release baseline: changes since tag `v0.1.0-internal.1`.

Commit-range reference (`v0.1.0-internal.1..HEAD`):

- `e895bf2` Multi-index for gene similarity
- `a156e88` Preparation of sequence alignments - bio::alignment::pairwise::banded
- `5737c5a` Code reorganisations, modularisation of indices, tutorials
- `c7970a9` Code restructuring
- `7c4bdd1` Progress indicator triggered every two seconds
- `0f3b0eb` Fixed CPU usage in Expert Splicing Window
- `213aeba` Improvements to cDNA mapping, Preparing for priorisation of alternative cloning paths
- `6b42736` cDNA mapping defaults to 10mers (from 9mers)
- `6f09fad` DNA feature grouping activated by default
- `0d45084` Adding test data to the archive
- `676b8a0` Splicing cDNA filter, layout fixes
- `70a557b` RNA-seq mapping initiated, removed translations
- `5df59f6` Removed some sqlite cruft from original GENtle
- `ebcd0f5` Fixes
- `bb06fca` ClawBio, hardening
- `6d46246` More towards cloning patterns
- `e0a391e` More consistent export
- `2fafee1` Hardening blast searches
- `d2a91bd` Improvements to Splicing Expert window - coloring
- `6870f2f` Scrolling in splicing expert window
- `f92e60b` Fixed genomes::tests::
- `fb4a909` More hardening, more tutorials, also on Gibson
- `897b858` Uniprot import. More towards primer design
- `52bf3bf` Implemented mapping of isoforms as show case

This internal release focuses on deterministic agent integration, cloning
workflow standardization, expert-view improvements, and release-readiness
hardening.

## Final Cut Additions (latest bits before `v0.1.0-internal.2`)

- RNA-read interpretation and sparse-origin work moved forward substantially:
  - default seed length moved to 10-mers
  - multi-index/similarity and alignment-preparation groundwork landed
  - cDNA mapping responsiveness improved (progress cadence and CPU behavior)
- GUI defaults and clarity improved:
  - feature grouping now defaults to active
  - splicing/cDNA layout and expert window behavior hardened
- Tutorial screenshots were added for upcoming visual docs integration:
  - TP73/luciferase and genomic retrieval-related captures:
    - `docs/screenshots/screenshot_GUI_main_with_TP73_retrieved.png`
    - `docs/screenshots/screenshot_GUI_DNA_sequence_with_TP73_retrieved.png`
    - `docs/screenshots/screenshot_GUI_prepare_reference_genome.png`
    - `docs/screenshots/screenshot_GUI_retrieve_genome_sequence.png`
    - `docs/screenshots/screenshot_GUI_retrieve_genomic_sequence_2.png`
  - Dotplot/cDNA-vs-genomic tutorial captures:
    - `docs/screenshots/tutorial_cdna_genomic_01_main_start.png`
    - `docs/screenshots/tutorial_cdna_genomic_02_fetch_cdna_dialog.png`
    - `docs/screenshots/tutorial_cdna_genomic_06_cdna_sequence_window.png`

## Highlights

- Python interface added (`gentle-py`) as a thin deterministic adapter over
  `gentle_cli`:
  - package path: `integrations/python/gentle_py/`
  - methods: `capabilities`, `state_summary`, `op`, `workflow`, `shell`
  - structured `GentleCliError` for automation-grade failure handling
  - packaging metadata added via
    `integrations/python/pyproject.toml`
- ClawBio/OpenClaw integration scaffold added:
  - `integrations/clawbio/skills/gentle-cloning/`
  - includes `SKILL.md`, wrapper runner, request examples, tests, and
    ready-to-paste `catalog_entry.json`
- Splicing expert usability and interpretability improved:
  - standalone window with vertical scrolling
  - transition-frequency-aware matrix coloring
  - predicted exon->exon transition matrix
  - CDS phase/modulo cues and expanded overlay semantics
  - deterministic SVG export parity for expert views
- BLAST robustness and async orchestration expanded:
  - stronger async job handling (`start/status/cancel/list`)
  - hardened invocation/provenance tracking
  - clearer integration across shell/CLI/MCP/agent paths
- Genome extraction annotation projection expanded and hardened:
  - `ExtractGenomeRegion` now supports explicit
    `annotation_scope=none|core|full` policy
  - optional `max_annotation_features` deterministic safety cap for heavy
    intervals (`full -> core -> none` fallback sequence)
  - structured operation telemetry emitted via
    `OpResult.genome_annotation_projection` (requested/effective scope, attached
    vs dropped feature counts, fallback reason)
  - shell/CLI flags added:
    `--annotation-scope` / `--max-annotation-features`
  - legacy flags retained for compatibility:
    `--include-genomic-annotation` / `--no-include-genomic-annotation`
- Agent adapter/runtime resilience expanded:
  - OpenAI native + OpenAI-compatible endpoints (including local runtimes)
  - adapter availability checks and clearer operator-facing error guidance
  - stricter schema/error handling and configurable timeout behavior
- Cloning pattern/routine coverage significantly expanded:
  - richer routine catalog entries and pattern templates across
    Gibson/Golden Gate/Gateway/TOPO/TA-GC/In-Fusion/NEBuilder families
  - stronger preflight/validation and run-bundle consistency
- Primer backend checks are now unified before running primer design:
  - new command:
    `primers preflight [--backend ...] [--primer3-exec ...]`
  - this runs a dry check first (tool available, version readable, setup looks
    valid) before execution
  - GUI now uses the same preflight check/report path as shell/CLI
  - deterministic tests now cover CLI forwarding for `primers preflight`
- Tutorial and workflow assets expanded:
  - new canonical online/tutorial examples (including TP53/TP63 tracks)
  - TP73 promoter luciferase planning skeleton workflow and GUI tutorial
    materials added
  - TP73 showcase/tutorial content consolidated into one canonical page:
    `docs/tutorial/tp73_promoter_luciferase_gui.md`
- UniProt mapping functionality expanded:
  - import/fetch/projection routes and shell/CLI exposure improved
- Architecture/protocol/docs hardening:
  - clearer schema/adapter guidance for agent interfaces, MCP, and external
    orchestration
  - stricter machine-readable interface posture in protocol docs

## Notable Changes by Area

### 1) Python and External Agent Interfaces

- Added thin Python adapter module + tests + packaging:
  - `integrations/python/gentle_py/client.py`
  - `integrations/python/tests/test_client.py`
  - `integrations/python/pyproject.toml`
- Added ClawBio/OpenClaw skill bridge:
  - `integrations/clawbio/skills/gentle-cloning/*`
- Updated docs for agent/tutorial/interface clarity:
  - `docs/agent_interfaces_tutorial.md`
  - `docs/architecture.md`
  - `docs/protocol.md`
  - `docs/cli.md`

### 2) Splicing and Feature Expert Views

- Extended splicing expert model + rendering:
  - matrix and transition-frequency enhancements
  - additional evidence visualization and phase/modulo cues
  - improved viewport behavior for dense transcript sets

### 3) BLAST, Genomes, and Track Handling

- Hardened BLAST orchestration and provenance metadata in shell/engine paths.
- Improved genome-path reliability and test stabilization around genome
  operations.
- Added scoped annotation transfer controls for region extraction:
  - core/default transfer keeps overlapping gene/transcript context
  - full transfer mode includes richer annotation payloads where available
  - deterministic fallback behavior for capped feature budgets
- Added extract-region telemetry fields for adapter/UI reporting:
  - enables GUI status and machine clients to inspect projection outcomes
    without parsing free-text messages
- Extended shell/CLI parser coverage with deterministic tests for new
  extract-region options and defaults.

### 3a) Primer Backend Preflight and Diagnostics

- Preflight means a dry-run check before execution: verify backend/tool
  availability and report what would be used, without changing project state.
- Added one shared preflight report path for Primer3 availability/version
  checks (engine + shell + CLI + GUI).
- GUI preflight status now shows the same shared report payload used by
  shell/CLI.
- Added deterministic tests for:
  - shell parsing/execution of `primers preflight`
  - CLI forwarded `primers preflight` parser and dispatch parity

### 4) Cloning Routine Catalog and Macro Workflows

- Expanded template catalog and routine metadata coverage.
- Improved preflight and execution consistency for multi-step cloning plans.

### 5) Tutorial / Demonstration Assets

- Expanded tutorial manifest and generated chapters.
- Added TP73 promoter luciferase planning skeleton and GUI-focused walkthrough.
- Added two dedicated GUI tutorials:
  - `docs/tutorial/tp73_promoter_luciferase_gui.md`
  - `docs/tutorial/two_sequence_dotplot_gui.md`
- Added screenshot assets for these tutorial tracks (for release/docs embedding):
  - TP73/luciferase and genomic retrieval screenshots under `docs/screenshots/`
  - cDNA vs genomic dotplot screenshots under `docs/screenshots/tutorial_cdna_genomic_*.png`

## Developer Notes

- CI/build guardrails now include broader validation of deterministic behavior
  in key adapter routes.
- Runtime and docs now explicitly track compatibility with OpenAI-compatible
  provider endpoints (including ClawBio-relevant integration pathways).

## Breaking Changes

No intentional protocol-major breaking changes were introduced in this internal
release. Existing project states and primary CLI operation/workflow routes
remain on the same schema family (`v1` contracts).

## Migration and Compatibility Notes

- `ExtractGenomeRegion` now defaults to `annotation_scope=core` when
  annotation scope is not explicitly set.
- Legacy flags are preserved for compatibility and map to scope semantics:
  - `--include-genomic-annotation` maps to `annotation_scope=core`.
  - `--no-include-genomic-annotation` maps to `annotation_scope=none`.
- For deterministic automation behavior, prefer explicit
  `--annotation-scope` and optional `--max-annotation-features` over legacy
  include/omit flags.

## Validated Commands (this release candidate)

- `cargo check -q`
- `cargo test -q`
- `cargo test -q test_extract_genome_region_include_annotation_attaches_features_and_sets_name`
- `cargo test -q test_extract_genome_region_full_scope_feature_cap_falls_back_to_core`
- `cargo test -q parse_genomes_extract_region_with_annotation_flag`
- `cargo test -q parse_genomes_extract_region_with_scope_and_cap`
- `cargo test -q execute_genomes_extract_region_default_scope_core_with_telemetry`
- `cargo test -q execute_genomes_extend_anchor_creates_sequence`

## Known Gaps (still open after this release)

Priority-ordered tracking and execution detail live in `docs/roadmap.md`
(see section "2. Active known gaps (priority-ordered)").

- Primer3 deep parity and advanced backend normalization still require
  additional hardening.
- Some adapter-level helper routes remain to be promoted into first-class
  engine operations.
- Additional routine-family semantic depth and cross-tool benchmark parity are
  still in progress.
- GUI still exposes the simple include/omit annotation toggle for region
  extraction; direct in-GUI `core|full|none` selector can be added later as an
  enhancement (engine/shell/CLI support is already in place).

## Suggested Internal Release Validation Checklist

1. Run `cargo check -q`.
2. Run selected high-value tests (engine shell/parity/genome/feature expert).
3. Smoke-test:
   - GUI splicing expert rendering + SVG export
   - BLAST async flow from shell/CLI
   - extract-region scope/cap behavior from shell/CLI:
     - default (`core`) extraction
     - `--annotation-scope full --max-annotation-features N` fallback path
   - Python adapter (`integrations/python`)
   - ClawBio scaffold demo (`integrations/clawbio`)
4. Verify tutorial assets and TP73 luciferase planning docs are accessible.
