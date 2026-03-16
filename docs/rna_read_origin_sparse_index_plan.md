# RNA Read Origin Classification + Sparse Multi-Gene Seed Index Plan

Status: in progress (phase-1 schema/operation scaffolding landed; list/sample-sheet sparse-origin provenance surfaced; multi-gene sparse index + ROI capture still pending)

## Summary

This plan extends the current RNA-read seed-filter workflow with a deterministic
multi-gene origin model:

- keep current strict TP73-oriented behavior as default-compatible baseline
- add optional multi-gene targeting in one run
- replace per-gene index fan-out with one sparse annotated seed index
- classify reads by coherent origin class instead of only pass/fail
- prepare antisense and annotation-incomplete ROI inspection workflows

## Problem statement

Current filtering can produce reads with strong local TP73 similarity but weak
global coherence across full read span. These reads are biologically relevant
for interpretation but are hidden by a single strict gate.

At the same time, conserved regions (for example domain-rich sequence blocks)
can contribute seed support across multiple genes/transcripts. We need
explicit origin diagnostics, not only a single scalar score.

## Goals

1. Preserve current strict TP73 gate as a stable baseline (no regression).
2. Add deterministic origin classes for interpretation and review.
3. Support optional multi-gene target sets in one combined run.
4. Keep strand evidence explicit and auditable in both-strands mode.
5. Use one sparse annotated seed index (scalable, no duplicated per-gene
   indexes).

## Non-goals (this track)

- full whole-transcriptome alignment/quantification inside GENtle
- automatic online retrieval during scoring loops
- adapter-local biology logic

## Core design decision

Use one annotated sparse inverted index:

- `seed_bits -> postings[]`

Each posting carries:

- `gene_id`
- `transcript_id`
- `strand`
- coordinate (`genomic` and optional `exonic_compact`)
- context (`exon`, `junction`, optional `intron`)
- occurrence metadata for weighting

Optional aggregate sparse matrix views (CSR/CSC style) may be derived for fast
gene/transcript score accumulation but are secondary representations, not
independent sources of truth.

## Input and catalog model

1. Splicing Expert accepts a target-gene set (one or many genes).
2. Engine builds transcript templates from locally available annotation first.
3. Additional transcript resources can be imported explicitly and versioned in
   metadata (source, timestamp, retrieval path).
4. Scoring always runs against persisted deterministic templates, never against
   live remote state.

## Read interpretation pipeline

### Stage A: seed accumulation (fast pass)

- full-read seed hashing (existing policy)
- accumulate posting hits across all targeted genes/transcripts
- compute:
  - raw seed fraction
  - weighted seed fraction (occurrence-aware)
  - unique matched seed count
  - strand-partitioned support totals

### Stage B: coherence refinement (top candidates only)

- select top-K transcript/gene candidates from Stage A
- compute chain/coherence metrics and transition support per candidate:
  - chain support fraction
  - median transcript gap and gap count
  - confirmed transitions / total transitions
- optional stricter confirmation with longer seeds on top candidates
  (`k=11/13`) to suppress conserved-domain bleed-through

### Stage C: origin class assignment

Assign deterministic class labels with reason metadata:

1. `TP73_coherent`
2. `TP73_partial_local_block`
3. `ROI_same_strand_local_block`
4. `ROI_reverse_strand_local_block`
5. `TP_family_ambiguous`
6. `Background_likely`

Also emit:

- `strand_confidence` (`high` or `low`)
- `origin_confidence` (numeric + reason string)
- top competing candidate diagnostics

## Strand policy

- Keep one combined run for both strands when requested.
- Do not hide strand ambiguity:
  - emit per-read `selected_strand`, `competing_opposite_strand`,
    `ambiguous_near_tie`
  - expose strand-partitioned support counters in aggregate views
- If read orientation is uncertain, classify as strand-ambiguous rather than
  forcing same/reverse labels without confidence.

## Antisense and annotation-incomplete support

Add an annotation-independent ROI seed-capture layer in parallel to transcript
templates:

- genomic ROI seeds, strand-tagged, context-aware
- supports antisense discovery even when full antisense transcript models are
  incomplete
- preserves "DNA-coated magnetic bead" mental model (seed as capture handle)

## Engine/API additions (planned)

- extend RNA-read report schema with:
  - `origin_class`
  - `origin_confidence`
  - `strand_confidence`
  - candidate contribution table (top gene/transcript contributors)
- add optional operation params:
  - `target_gene_ids[]`
  - `origin_mode=single_gene|multi_gene_sparse`
  - `roi_seed_capture_enabled`

## GUI/CLI/agent surfaces (planned)

GUI:

- target-gene multiselect
- origin-class summary panel and filters
- per-read contributor breakdown

CLI/shared shell/MCP:

- parity fields in list/show/export routes
- deterministic class-aware filtering in export commands

## Performance strategy

1. Batch-based read processing with deterministic merge order.
2. Parallel Stage A/Stage B compute on worker pool.
3. Keep bounded top-hit heaps and bounded candidate lists.
4. Persist per-batch timing counters for `seed`, `infer`, and candidate fan-out
   diagnostics.

## Phased rollout

### Phase 1

- schema scaffolding for class labels + candidate contribution diagnostics
- non-breaking defaults (`single_gene` behavior unchanged)

### Phase 2

- sparse multi-gene index builder
- multi-gene target-set UI + CLI/shell params

### Phase 3

- deterministic origin-class assignment + GUI class summaries
- strand-confidence and ambiguity reporting

### Phase 4

- ROI seed-capture layer for antisense/annotation-incomplete regions
- performance tuning and deterministic benchmark fixtures

## Acceptance criteria

1. Existing TP73 strict mode remains behavior-compatible when multi-gene mode
   is off.
2. Multi-gene mode emits deterministic origin classes and class reasons.
3. Same input yields adapter-equivalent class/report/export output across
   GUI/CLI/shell/MCP.
4. Strand ambiguity is explicit in report payloads (never implicit).

## Test matrix (planned)

1. Regression:
   - current TP73 strict tests remain green under default mode.
2. Positive controls:
   - compact committed pack under `test_files/fixtures/mapping/` classifies as
     coherent TP73 support for TP73-family sets.
3. Negative controls:
   - TP53 set in `test_files/fixtures/mapping/` should not classify as
     `TP73_coherent`.
   - optional extended decoy runs may use
     `test_files/mapping/False_TP73/*` locally.
4. Cross-family controls:
   - TP53/TP63 sets produce ambiguous/family/background classes as expected.
5. Determinism:
   - identical results under repeated runs and different worker counts.

## Risks and constraints

- conserved-domain seed reuse can inflate ambiguous classes without proper
  weighting and coherence gating.
- multi-gene scope expansion can increase compute cost if candidate pruning is
  not bounded.
- external annotation imports must remain provenance-tracked and versioned to
  keep analyses reproducible.
