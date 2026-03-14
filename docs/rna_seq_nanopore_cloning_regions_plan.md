# RNA-seq Evidence for Cloning-Candidate Regions (Nanopore cDNA First)

Status: in progress (phase-1 baseline implementation underway)

Related follow-up plan:

- `docs/rna_read_origin_sparse_index_plan.md` (multi-gene sparse seed index,
  origin classes, antisense-ready ROI seed-capture model)

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
- Keep phase-1 execution seed-filter-only (no inline dynamic-programming
  alignment) so streaming progress remains responsive on large FASTA inputs.
- Keep read-orientation mode explicit in UI/CLI:
  - cDNA mode: optional poly-T prefix reverse-complement normalization
  - direct-RNA mode: no automatic reverse-complement normalization
- Produce deterministic region-level summary outputs for inspection.
- Persist per-report exon-support and exon-exon junction-support frequency
  summaries so cohort-level downstream interpretation can reuse the same
  engine-authored metrics.
- Support batch sample-sheet export (TSV) for many runs/files, including
  machine-readable frequency columns for each sample/report.
- Export additional deterministic TSV artifacts:
  - per-read exon-path table
  - exon/transition abundance table
- Keep the initial seed-hit gate deterministic (`min_seed_hit_fraction=0.30`)
  as a bootstrap default, with explicit plan to replace/tune it once
  transcriptome-scale background/noise estimates are available.

### Phase 2: engine-owned summarized evidence model

- Promote imported evidence into engine-owned summarized views for ROI-level
  interpretation.
- Include coverage/junction/exon-support style summary metrics scoped to the
  active ROI.
- Implement seed-hash mapping cascade for shortlisted reads (no DP requirement):
  - hash full reads in coding orientation and first map against exon-only seed
    space for the scoped ROI
  - unresolved seeds then map against exon-exon junction seed space
  - remaining unresolved seeds optionally map against intronic seed space
  - remaining unmatched seeds are retained as noise/error indicators and
    potential trans-splicing signal for downstream clustering/inspection
- Add deterministic sample-sheet merge/append semantics so large cohorts can be
  assembled without adapter-specific post-processing logic.
- Keep seed-capture workflow abstraction documented:
  model seed-hash filtering as a reusable biotech-style enrichment operation
  that can be composed in standard workflows (deferred implementation).
- Expose inspect/export parity for those summaries across GUI/shared shell/CLI
  routes.

### Phase 2b: dual-strand joint-run hardening (all-overlap scope)

Goal: keep one combined run over `+/-` transcript lanes while making
cross-strand evidence assignment explicit, auditable, and deterministic.

- Keep `all-overlap / both-strands` as one run contract:
  no independent per-strand reruns required.
- Add strand-partitioned seed diagnostics per read:
  - best `+` transcript-chain score tuple
  - best `-` transcript-chain score tuple
  - final selected strand and deterministic reason string
- Add strand-partitioned aggregate diagnostics per isoform row:
  - assigned reads from same-strand evidence
  - competing opposite-strand evidence counts
  - ambiguity counters (near-tie rows) for review
- Add deterministic tie-break policy for read assignment:
  1. higher confirmed transition count
  2. higher exon support count
  3. lower median transcript gap
  4. strand-consistent chain preference (if configured strict)
  5. lower transcript feature id
- Add optional strict mode:
  require chain assignment and transition-confirmation support from the same
  strand partition before a read can pass seed gate.
- Keep default mode permissive but transparent:
  mixed-strand seed reuse remains allowed, but diagnostics must disclose when
  final assignment depends on cross-strand shared seeds.
- Extend report schema and exports:
  include per-read strand-assignment diagnostics and aggregate ambiguity
  counters so GUI/shell/CLI/agents remain contract-equivalent.

### Phase 3: suggestion-to-curation handoff

- Add RNA-evidence-based curation suggestions integrated with the feature-edit
  track.
- Keep mutation paths explicit: apply-suggestion actions remain user-confirmed,
  deterministic engine operations only (no GUI-only mutation logic).
- Add evidence-to-primer handoff:
  seed-support and origin diagnostics should be reusable as deterministic
  guidance input for primer-target ranking, while final primer scoring remains
  owned by shared primer-design engine contracts.

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
- Dual-strand joint-run assignment is deterministic and auditable:
  - one combined run can include both strands,
  - per-read strand decision diagnostics are exported,
  - strict same-strand mode is available when higher specificity is needed.

## Risks and constraints

- Long-read error profile handling is filter-sensitive and must be explicit.
- Initial mapping is intentionally partial and ROI-oriented.
- Seed-space design must stay strand/direction aware so exon-only matching does
  not over-count reverse/non-coding artifacts.
- Full-scale mapping/quantification remains a later track.
- Transcriptome-scale SNR estimation requires representative background data;
  early calibration may need iterative retuning by dataset type/species.

## Test and validation expectations (implementation phase)

- Deterministic tests for ROI evidence import/filtering and summary stability.
- Cross-adapter parity tests for summary visibility and export contracts.
- Dual-strand joint-run tests:
  - deterministic winner selection in `+/-` near-tie scenarios
  - strict mode rejects mixed-strand-supported reads
  - report/schema parity for strand diagnostics in GUI/shell/CLI exports
- Documentation consistency checks:
  - `docs/architecture.md` references this plan.
  - `docs/roadmap.md` references this plan.
  - no conflicting demo-first framing in this track.
