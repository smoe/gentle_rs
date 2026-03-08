# Release Notes: `v0.1.0-internal.2`

Release baseline: changes since tag `v0.1.0-internal.1`.

Commit-range reference (`v0.1.0-internal.1..HEAD`):

- `52bf3bf` Implemented mapping of isoforms as show case
- `897b858` Uniprot import. More towards primer design.
- `fb4a909` More hardening, more tutorials, also on Gibson
- `f92e60b` Fixed genomes::tests::
- `6870f2f` Scrolling in splicing expert window
- `d2a91bd` Improvements to Splicing Expert window - coloring
- `2fafee1` Hardening blast searches
- `e0a391e` More consistent export.
- `6d46246` More towards cloning patterns

This internal release focuses on deterministic agent integration, cloning
workflow standardization, expert-view improvements, and release-readiness
hardening.

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
- Tutorial and workflow assets expanded:
  - new canonical online/tutorial examples (including TP53/TP63 tracks)
  - TP73 promoter luciferase planning skeleton workflow and GUI tutorial
    materials added
  - web-ready narrative showcase page added:
    `docs/tutorial/tp73_promoter_luciferase_showcase.md`
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

### 4) Cloning Routine Catalog and Macro Workflows

- Expanded template catalog and routine metadata coverage.
- Improved preflight and execution consistency for multi-step cloning plans.

### 5) Tutorial / Demonstration Assets

- Expanded tutorial manifest and generated chapters.
- Added TP73 promoter luciferase planning skeleton and GUI-focused walkthrough.

## Developer Notes

- CI/build guardrails now include broader validation of deterministic behavior
  in key adapter routes.
- Runtime and docs now explicitly track compatibility with OpenAI-compatible
  provider endpoints (including ClawBio-relevant integration pathways).

## Breaking Changes

No intentional protocol-major breaking changes were introduced in this internal
release. Existing project states and primary CLI operation/workflow routes
remain on the same schema family (`v1` contracts).

## Known Gaps (still open after this release)

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
