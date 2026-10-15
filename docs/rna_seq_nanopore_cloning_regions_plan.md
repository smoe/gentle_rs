# RNA-Read Evidence For Cloning Regions: Status And Remaining Work

Source review: 2026-10-07. The Nanopore cDNA mapping/report baseline is
implemented; it is not merely phase-1 scaffolding. This residual plan does not
authorize a second mapping engine.

## Existing Baseline

The shared RNA-read engine supports FASTA intake, transcript seed filtering,
explicit read-orientation handling, exon-path/junction support, multi-gene
targeting and strand/origin diagnostics. Interpretation's phase-1 streaming pass
is seed-filter-only; the separate mapping/alignment path supplies retained-read
alignment. Do not promise that interpretation itself performs inline DP.

Persisted reports drive the RNA Mapping GUI and shared inspect/export routes.
Batch sample sheets, per-read exon paths, exon/transition abundance, score-density
and alignment exports already exist. Samples retain their own report/resource
context; aggregate counts are not automatically normalized expression.

The [origin plan](rna_read_origin_sparse_index_plan.md) distinguishes actual
multi-gene targeting from the still-inactive annotation-independent ROI capture
option. Canonical interfaces are in [CLI](cli.md), [protocol](protocol.md) and
[GUI](gui.md); implementation is in `src/engine/analysis/rna_reads.rs` and
`src/main_area_dna/rna_read_mapping_ui.rs`.

## Remaining Work

### Calibration And Background

Positive/control preflight and explicit threshold recommendations already
exist. They do not establish a calibrated transcriptome-wide background model.
Before replacing bootstrap thresholds with SNR-normalized bands, retain
representative background/control data, species/library/error-profile context,
a denominator and reproducible threshold rules. Validate false-positive and
false-negative behavior; do not infer trans-splicing from unmatched seeds or
expression changes from unnormalized abundance.

### Suggestion-To-Curation And Primer Guidance

A dedicated RNA-evidence curation-suggestion contract remains a proposal.
Define a content-bound report of proposed feature changes, evidence and reasons
before a separate, confirmed shared-engine mutation. Existing general feature
editing is not evidence that this automatic suggestion pipeline exists.
RNA support may guide candidate choice, but cannot replace complete primer
specificity, coverage and oligo-QC gates.

### Cohort Presentation And Scale Acceptance

The [roadmap](roadmap.md) retains virtualized multi-sample/gene-cohort alignment
display. Use existing engine-authored support and sample-sheet fields rather
than recomputing biology in the GUI. Large-cohort throughput, calibrated origin
confidence and live GUI acceptance require their own measured evidence.

Transposon-specific modeling remains deferred pending a defined biological
question. Full whole-transcriptome alignment/quantification, unsolicited remote
retrieval and automatic annotation changes are not this track.

## Verification For Future Slices

Preserve default single-gene behavior, explicit missing-target reporting,
any-strand near-tie diagnostics and exact selected-subset export semantics.
For each new proposal, require deterministic engine tests, adapter-equivalent
inspection and provenance-bound fixtures. Keep real-data calibration and native
GUI acceptance separate from synthetic regression results.
