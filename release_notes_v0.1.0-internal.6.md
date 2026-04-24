# Release Notes / Changelog: `v0.1.0-internal.6` (draft)

This draft internal release is centered on regulatory evidence, motif-aware
reasoning, and more complete automation handoff. Compared with
`v0.1.0-internal.5`, this cut moves GENtle from protein/construct planning
toward evidence-backed regulatory interpretation: CUT&RUN can now be organized
as a local evidence family, TFBS score landscapes can be inspected and compared
more directly, RNA-read mapping is more usable from CLI and GUI paths, and
variant/promoter follow-up has a stronger planning story.

The release is also a practical integration cut. More features now have shared
engine/shell/CLI records rather than GUI-only behavior, and ClawBio-facing
request examples now cover a wider span of GENtle's regulatory, protein,
sequence-inspection, and service-readiness surfaces.

## Highlights

- CUT&RUN support has advanced from planning to a real shared evidence track:
  - catalog-backed dataset discovery, status, preparation, and projection
  - anchored BED/BigWig projection for processed peaks and signal
  - ROI-first FASTA/FASTQ read interpretation with saved read reports
  - per-base coverage, cut-site, and fragment export
  - regulatory-support reasoning over peaks, signal islands, and read clusters
  - motif-context summaries for strong supported windows, including cases where
    the assayed factor has no direct local motif
- TFBS/JASPAR analysis became much more useful for regulatory review:
  - stateless sequence scans and score-track routes are now shared
  - GUI cached inspectors exist for TFBS hit scans, score tracks, and
    similarity ranking
  - TFBS score-track similarity and correlation reports gained richer
    strand-aware output
  - JASPAR metadata, expert views, and motif-group/resource summaries are
    available to both local and ClawBio-facing workflows
- Promoter and variant follow-up gained a more explicit interpretation path:
  - multi-gene promoter TFBS characterization landed
  - VKORC1 / rs9923231 promoter-reporter planning artifacts were expanded
  - nf-core-inspired decision context and experimental follow-up request
    catalogs were added for variant/regulatory reasoning
- RNA-read mapping became more reviewable and more scriptable:
  - CLI parity now covers retained-read/report operations that previously felt
    GUI-biased
  - read-length statistics and result-export paths were added
  - transcript-catalog indexing supports concatemer and artifact review
  - several GUI stability and window-placement fixes landed for the RNA mapping
    workspace
- Protein and peptide workflows continued to broaden:
  - first-class protease digestion is now a shared operation path for protein
    and peptide sequences
  - protein gel and 2D gel workflows were added to ClawBio/demo surfaces
- Interoperability and automation improved:
  - SnapGene `.dna` parsing is now available through the bundled reusable
    parser package
  - service/resource status reporting and handoff/doctor flows were expanded
    for ClawBio and agent-driven execution
  - Ensembl gene download/import paths and direct sequence-inspection examples
    broaden the "ask for context, inspect locally" workflow

## Notable Changes by Area

### 1) CUT&RUN Evidence and Regulatory Reasoning

- Added the shared CUT&RUN catalog/cache/report family:
  - `ListCutRunDatasets`
  - `ShowCutRunDatasetStatus`
  - `PrepareCutRunDataset`
  - `ProjectCutRunDataset`
- Processed CUT&RUN datasets can now be projected onto genome-anchored
  sequences using existing BED/BigWig track semantics instead of a separate
  one-off projection path.
- Dataset preparation now uses lease/heartbeat status records so duplicate
  prepare calls can reuse an active install and stale attempts are visible.
- Added ROI-first read interpretation:
  - `InterpretCutRunReads`
  - `ListCutRunReadReports`
  - `ShowCutRunReadReport`
  - `ExportCutRunReadCoverage`
- V2 read reports retain explicit read-unit status, fragment spans, coverage,
  cut-site counts, and compact support clusters.
- Added V3 regulatory support inspection:
  - `InspectCutRunRegulatorySupport`
  - shared shell/CLI route:
    `cutrun inspect-regulatory-support ...`
- V3 aggregates prepared peaks, prepared signal islands, and saved ROI read
  reports into support windows.
- Theoretical TFBS rows are split into `confirmed` and `unconfirmed`.
- Strong supported windows without the assayed target motif are reported as
  motif-context windows and classified as:
  - `context_supported_by_other_motifs`
  - `motif_poor_supported`
- The latest reasoning pass now handles non-DNA-binding/cofactor-style evidence
  more honestly: when a selected CUT&RUN target factor cannot be resolved to a
  local motif, strong supported windows are still retained and reported through
  context-only motif reasoning.
- Motif-context scans use a high-confidence default
  `motif_context_min_llr_quantile = 0.95` and omit resolved target motifs from
  the "other motifs" context list.

### 2) TFBS, JASPAR, ATtRACT, and Promoter Analysis

- Stateless sequence inspection now covers direct TFBS/JASPAR hit scanning and
  score-track summarization without forcing project-state materialization.
- TFBS score tracks now use the same state-optional `SequenceScanTarget`
  contract as other quick sequence scans.
- GUI support expanded:
  - direct TFBS scan actions from the DNA window
  - cached TFBS hit inspector
  - cached TFBS score-track inspector/export
  - cached TFBS similarity ranking
  - click-through navigation from hit tables back to sequence spans
- TFBS score-track comparison broadened:
  - strand-specific similarity/correlation axes
  - Spearman support for promoter-correlation views
  - positive-only score-track visualization
  - imported BED overlays can be shown under TFBS score tracks
- JASPAR support is stronger:
  - motif metadata includes species/class/family context
  - JASPAR expert presentation is available
  - one-entry JASPAR expert display was tightened
  - benchmark/helper tooling exists for motif-interest catalogs
- ATtRACT support moved from exploratory to a conservative PWM-backed layer:
  - deterministic snapshot fingerprinting
  - provenance-first PWM records
  - strict-vs-windowed comparison behavior documented
  - splice-site motif annotation and intron-signal heuristics improved
- Multi-gene promoter TFBS characterization landed, including SVG/rendering
  support and ClawBio request examples for TERT/TP73-style promoter comparison.

### 3) RNA Mapping, Concatemer Review, and Read Artifacts

- RNA-read mapping gained CLI parity for retained-report/result workflows and
  related command parsing.
- Read-length statistics are now reported for RNA mapping, improving quick
  diagnosis of suspicious input distributions.
- RNA-mapping result export paths were added.
- A reusable transcript-catalog index supports concatemer/audit workflows.
- GUI cleanup around naming, report IDs, and concatemer review reduced
  confusion when moving between active runs and saved reports.
- Several RNA-mapping window fixes landed:
  - stale detached hosts removed
  - RNA mapping no longer closes erroneously when a parent closes
  - window arrangement and foreground behavior were corrected
  - redundant scrolling was removed
  - stability fixes addressed recent GUI instability
- Nanopore adapter/barcode resources and artifact-oriented feature handling
  were added as groundwork for better long-read diagnostics.

### 4) Variant, Promoter, and Experimental Follow-Up Planning

- Variant Follow-up now has a fuller GUI story, including tutorial coverage
  and a more explicit expert handoff.
- VKORC1 / rs9923231 promoter-reporter planning artifacts were expanded with
  reproducibility commands, candidate reports, context JSON, and SVG output.
- Annotation-candidate curation and write-back continued to mature, allowing
  reasoning-derived suggestions to become project-visible annotations.
- The experimental follow-up layer broadened:
  - machine-readable request catalog
  - follow-up catalog graph generator
  - perturbation-planning contract groundwork
  - nf-core-inspired decision-context additions for genetic variants
- qPCR follow-up became more concrete:
  - qPCR reports
  - qPCR protocol cartoons
  - qPCR seed helpers
  - a dedicated place in the PCR Designer

### 5) Protein, Protease, Gels, and SnapGene Interoperability

- Added a shared protease-digest path for first-class protein/peptide
  sequences, with engine, shell, CLI, GUI, and tests moving onto one contract.
- Protein-gel and 2D protein-gel render/export paths are now available and used
  in tutorial/ClawBio examples.
- Ensembl/UniProt audit workflows continued to expand, including direct
  Ensembl gene download/import routes and TP73 UniProt projection audit
  material.
- SnapGene `.dna` import support landed through a reusable parser package under
  `packages/snapgene-reader`.
- Reverse joined SnapGene feature locations were fixed after the initial parser
  landing.

### 6) ClawBio, Services, and Automation Readiness

- ClawBio request coverage expanded substantially:
  - stateless sequence inspection
  - Ensembl preflight and gene extraction/import
  - TFBS hit scans, summaries, score tracks, and similarity
  - JASPAR resource summaries and TF-group resolution
  - VKORC1 context/planning workflows
  - protein gels and TP73 score-track workflows
- Service/resource readiness reporting is now more explicit:
  - service status
  - resource status
  - handoff/doctor checks
  - clearer local-checkout wrapper failure reporting
- GENtle version/status information is available through ClawBio-facing paths.
- A deterministic graphical GENtle demo path and PNG presentation support were
  added for ClawBio demonstrations.

### 7) GUI and Window Reliability

- Online tutorial project opening now runs in the background instead of
  blocking the main event loop.
- DNA sequence-window widget IDs were tightened to avoid collisions.
- Transition behavior from circular to linear representation was fixed and
  documented.
- Hosted-window title bars and RNA-mapping foreground behavior were corrected.
- Configuration/help/window work is still an active stabilization area, but
  this cut contains several concrete fixes that reduce window drift and stale
  hosted-window state.

### 8) Architecture, Release, and Internal Engineering

- The workspace split continued to absorb shared contracts into
  `gentle-protocol` and shared renderers into `gentle-render`.
- `src/engine/cutrun.rs`, motif-statistics helpers, promoter-design analysis,
  variant-promoter analysis, service-readiness, and resource-status paths now
  carry more of the newer shared behavior.
- Shell dispatch was split further to avoid stack overflows in dense command
  families, including arrangement/rack/export and CUT&RUN paths.
- Release workflow hardening continued:
  - release-installer workflow now targets published GitHub releases
  - release builds avoid stale target-cache reuse
  - release attributes record Linux distribution intent
- The local workspace is currently versioned as `0.1.0-internal.6` in
  `Cargo.toml` and `Cargo.lock`.

## Release-Facing Known Limitations

- CUT&RUN V3 is still a baseline reasoning surface:
  - support-window strength thresholds are not yet user-tunable
  - GUI inspection for V3 reports is still deferred
  - full FASTQ-scale indexing/alignment remains intentionally lightweight
    compared with dedicated genomics pipelines
- TFBS/JASPAR reasoning is much richer, but motif presence remains a
  computational hypothesis; CUT&RUN and other evidence should still be treated
  as contextual support, not definitive binding proof by itself.
- RNA-read mapping is more stable and scriptable, but long-read artifact and
  isoform interpretation remain under active tuning.
- Variant/follow-up planning now has more machinery, but practical experiment
  ranking by reagent availability, cost, and lab-specific feasibility is still
  early.
- SnapGene support is import-focused; richer SnapGene-specific presentation
  metadata and export remain future work.
- GUI window behavior has improved, but recent configuration and hosted-window
  edge cases remain worth manual smoke testing before tagging.

## Install / Package Notes

- macOS release artifact remains `.dmg`.
- Windows release artifact remains `.zip` containing `gentle.exe`.
- Linux native installer packaging remains deferred; release metadata should
  continue to default Linux distribution intent to `tarball`.
- Container/runtime story should continue to distinguish:
  - `ghcr.io/smoe/gentle_rs:cli` for headless/agent/ClawBio use
  - `ghcr.io/smoe/gentle_rs:gui` or `:latest` for interactive GUI container use
- For automation-driven work, prefer the headless CLI/MCP routes over the GUI
  container unless the task is explicitly interactive.

## Suggested Pre-Tag Smoke Checklist

At minimum before tagging `v0.1.0-internal.6`, run:

```bash
cargo check -q
cargo test -q inspect_cutrun_regulatory_support
cargo test -q rna_reads
cargo run --release --bin gentle -- --version
cargo run --release --bin gentle_cli -- capabilities
cargo run --release --bin gentle_examples_docs -- --check
```

Because this release touched regulatory reasoning, RNA mapping, GUI-hosted
windows, and ClawBio-facing automation, one short manual smoke pass should also
cover:

- opening the GUI and the Configuration window
- opening the RNA-read mapping workspace from a DNA window
- running one TFBS hit scan and opening its cached inspector
- running one CUT&RUN regulatory-support CLI/shell command on a small fixture
- invoking one ClawBio local-checkout request wrapper against `gentle_cli`

## Release-Shaped Smoke Check Results

Partial local release-prep checks run on 2026-04-24:

- `git diff --check`: pass
- `cargo check -q`: pass
- `cargo test -q inspect_cutrun_regulatory_support`: pass

The broader release-shaped matrix from `docs/release.md` remains pending before
final tagging.
