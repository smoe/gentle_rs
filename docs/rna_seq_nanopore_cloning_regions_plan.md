# RNA-seq Evidence for Cloning-Candidate Regions (Nanopore cDNA First)

Status: in progress (phase-1 baseline implementation underway)

## Goal and scope

Near-term objective: support cloning-oriented interpretation for small genomic
regions of interest (ROI) using deterministic RNA-seq evidence overlays.

In scope:

- small-region genomic workflows tied to cloning-candidate decisions
- exon-aware evidence display in ROI context
- suggestion-first curation support with explicit user confirmation before any
  feature edit

Out of scope (initial phase):

- full whole-transcriptome alignment/quantification pipeline inside GENtle
- transposon interpretation policy/modeling decisions

## Implementation phases

### Phase 1: initial Nanopore cDNA evidence baseline

- Ingest externally generated Nanopore cDNA-derived evidence for a selected
  genomic ROI.
- Apply deterministic filtering suitable for long-read error profiles.
- Support partial mapping overlays against genomic region/exon context.
- Produce deterministic region-level summary outputs for inspection.
- Persist per-report exon-support and exon-exon junction-support frequency
  summaries so cohort-level downstream interpretation can reuse the same
  engine-authored metrics.
- Support batch sample-sheet export (TSV) for many runs/files, including
  machine-readable frequency columns for each sample/report.
- Keep the initial seed-hit gate deterministic (`min_seed_hit_fraction=0.30`)
  as a bootstrap default, with explicit plan to replace/tune it once
  transcriptome-scale background/noise estimates are available.

### Phase 2: engine-owned summarized evidence model

- Promote imported evidence into engine-owned summarized views for ROI-level
  interpretation.
- Include coverage/junction/exon-support style summary metrics scoped to the
  active ROI.
- Add deterministic sample-sheet merge/append semantics so large cohorts can be
  assembled without adapter-specific post-processing logic.
- Expose inspect/export parity for those summaries across GUI/shared shell/CLI
  routes.

### Phase 3: suggestion-to-curation handoff

- Add RNA-evidence-based curation suggestions integrated with the feature-edit
  track.
- Keep mutation paths explicit: apply-suggestion actions remain user-confirmed,
  deterministic engine operations only (no GUI-only mutation logic).

### Phase 4: transcriptome-scale signal-to-noise calibration

- Add a formal signal-to-noise model for full-transcriptome inputs where
  cloning-ROI signal is expected to be a small minority.
- Compute an empirical background seed-hit distribution per run/sample (or
  cached cohort baseline) and derive dynamic acceptance bands from that
  distribution.
- Promote the fixed `30%` rule to a configurable bootstrap prior, then support
  SNR-normalized decision thresholds (for example z-score/quantile gates) for
  higher specificity on transcriptome-scale data.
- Persist SNR diagnostics in report metadata so GUI/shell/CLI/agents can inspect
  why reads were accepted/rejected under dynamic thresholds.

## Deferred transposon direction

Transposon-specific modeling and expression interpretation are explicitly
tracked as deferred and pending Anze-led direction.

## Acceptance criteria

- Deterministic ROI evidence summaries exist for the Nanopore cDNA input path.
- Batch sample-sheet export exists and preserves deterministic per-report
  exon/junction frequency summaries.
- Summary inspect/export semantics are consistent across GUI/shell/CLI.
- Engine remains the biology source of truth (no adapter-local biology logic).

## Risks and constraints

- Long-read error profile handling is filter-sensitive and must be explicit.
- Initial mapping is intentionally partial and ROI-oriented.
- Full-scale mapping/quantification remains a later track.
- Transcriptome-scale SNR estimation requires representative background data;
  early calibration may need iterative retuning by dataset type/species.

## Test and validation expectations (implementation phase)

- Deterministic tests for ROI evidence import/filtering and summary stability.
- Cross-adapter parity tests for summary visibility and export contracts.
- Documentation consistency checks:
  - `docs/architecture.md` references this plan.
  - `docs/roadmap.md` references this plan.
  - no conflicting demo-first framing in this track.
