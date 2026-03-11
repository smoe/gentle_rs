# RNA-seq Evidence for Cloning-Candidate Regions (Nanopore cDNA First)

Status: planned

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

### Phase 2: engine-owned summarized evidence model

- Promote imported evidence into engine-owned summarized views for ROI-level
  interpretation.
- Include coverage/junction/exon-support style summary metrics scoped to the
  active ROI.
- Expose inspect/export parity for those summaries across GUI/shared shell/CLI
  routes.

### Phase 3: suggestion-to-curation handoff

- Add RNA-evidence-based curation suggestions integrated with the feature-edit
  track.
- Keep mutation paths explicit: apply-suggestion actions remain user-confirmed,
  deterministic engine operations only (no GUI-only mutation logic).

## Deferred transposon direction

Transposon-specific modeling and expression interpretation are explicitly
tracked as deferred and pending Anze-led direction.

## Acceptance criteria

- Deterministic ROI evidence summaries exist for the Nanopore cDNA input path.
- Summary inspect/export semantics are consistent across GUI/shell/CLI.
- Engine remains the biology source of truth (no adapter-local biology logic).

## Risks and constraints

- Long-read error profile handling is filter-sensitive and must be explicit.
- Initial mapping is intentionally partial and ROI-oriented.
- Full-scale mapping/quantification remains a later track.

## Test and validation expectations (implementation phase)

- Deterministic tests for ROI evidence import/filtering and summary stability.
- Cross-adapter parity tests for summary visibility and export contracts.
- Documentation consistency checks:
  - `docs/architecture.md` references this plan.
  - `docs/roadmap.md` references this plan.
  - no conflicting demo-first framing in this track.
