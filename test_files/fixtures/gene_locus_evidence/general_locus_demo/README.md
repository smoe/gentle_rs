# General gene-locus evidence demo fixture

## Origin

Every file in this directory is hand-crafted synthetic data created for GENtle.
It contains no private sequence, experimental measurement, model output,
credential, or network-derived payload. The names `LOCUSDEMO`,
`ENSGDEMO0001`, and the `Synthetic1` chromosome are deliberately fictional.

## Contents

- `ensembl_gene_entry.json` is an offline Ensembl-shaped entry for a 1200 bp
  plus-strand locus. It has two versioned protein-coding transcripts with
  distinct 5' starts and distinct annotated translation starts. Its sequence
  begins with 170 `A` bases, the 26 bp motif island
  `GACATGTCTGGACATGTGGGGCGGGG`, and 104 `A` bases. The island contains a
  near-maximal TP73 site and a maximal SP1 site overlapping by one base. The
  remainder is `ATG`, 48 `A` bases, `ATG`, 543 `A` bases, `TAA`, and 300 `A`
  bases.
- `tp73_occupancy.bed` is synthetic locus-level TP73 occupancy evidence.
- `h3k4me3_context.bed` is a separately declared synthetic chromatin-context
  lane. Neither BED file represents a biological experiment.
- `occupancy_layout.json` assigns stable source ids, assembly, assay,
  mark/factor, roles, condition, cell line, batch, and independent scales. It
  also requests one deliberately absent SERPINE1 lane so `not_prepared` is
  visible rather than silently omitted. No role is inferred from a filename.
- `external_model_scores.json` is a provider-neutral, coordinate-bound set of
  synthetic forward/reverse model scores and two site calls with the typed
  `provider_declared_uncalibrated` state. No model was run.
- `preparation_request.json` composes those sources with four local JASPAR
  requests for TP73, PATZ1, E2F1, and SP1, a fixed 1000 bp scale bar, the
  canonical three-way reporter architecture comparison for both transcripts
  (six design rows total), and SVG, PNG, and PDF outputs. The combined figure
  includes an explicit material key for genomic source intervals, orange
  spliced exon/cDNA segments, and schematic red `LUC` blocks.

## Deterministic recreation

1. Recreate the sequence with the segments listed above and verify a length of
   1200 bp.
2. Verify its digest is
   `sha256:4813bc39b42785cf6731d92c07fc65bc60e646248d323d22f4525c9a599b2e64`.
3. Canonicalize the external score payload as the JSON object containing only
   `track_start_0based`, `forward_scores`, `reverse_scores`, and `sites`, using
   GENtle's serializer. Its declared score-payload digest is
   `sha256:984272df50f25143f47e8eddcbf53dbc8b003c7e68de43bef50e495d6b18aed7`.
   The request-schema and request bindings are SHA-256 of the literal synthetic
   strings `gentle.synthetic.regulatory.request_schema.v1` and
   `general_locus_demo external score request v1`, respectively.
4. BED coordinates are zero-based half-open and map to the fixture's
   `GRCh38/Synthetic1:100001..101200` one-based inclusive anchor.

## GENtle use

`prepare_gene_locus_evidence_runs_fully_offline_with_normalized_score_sources`
executes the typed `PrepareGeneLocusEvidence` operation directly against the
offline entry without network access or pre-existing project state. The test
verifies the canonical genome anchor, automatic isoform panel, two
transcript/TSS and CDS-start classes, explicit occupancy/chromatin roles, four
normalized JASPAR and external-model score tracks, one explicit unavailable
occupancy row, six canonical reporter architecture rows, and non-empty combined
SVG/PNG/PDF exports. The architecture report embeds the same
normalized locus-evidence report rendered in the figure; no reporter geometry
is re-derived by the display adapter. Genomically sourced architecture segments
remain to scale; each red `LUC` block is an explicitly schematic synthetic
continuation and does not claim a genomic span. The runbook uses the same
request as the public offline demonstration.

Local source paths are omitted from the portable report and receipt by default;
their SHA-256 bindings remain available for audit. Set
`include_local_source_paths=true` only when a deliberately machine-local report
is required.

This fixture proves data-flow, coordinate, provenance, and rendering behavior
only. PWM and external-model values are predictions, CUT&RUN-style rows are
synthetic occupancy evidence, H3K4me3 is chromatin context, and no layer proves
binding, affinity, isoform-specific regulation, or causality.
