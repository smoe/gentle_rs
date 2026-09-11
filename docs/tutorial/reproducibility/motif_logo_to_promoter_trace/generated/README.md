# TSS TFBS Profile Export

This export projects the supplied report without rescoring or reading source sequences.

## Interpretation Limits

JASPAR tracks are sequence-model predictions, not measured binding, affinity, occupancy, cofactor interaction, promoter activity or functional regulation. Cross-matrix magnitudes are not comparable without a documented calibration. Annotated transcript starts are TSS candidates, not experimentally established initiation sites. Correlation compares model outputs, not biological correctness or co-regulation. A factor also appearing as a target gene does not establish autoregulation.

Bundle consistency and hashes establish internal consistency, not independent reference authenticity.
The report's verification status is retained verbatim in report.json. Nothing here adds biological validation.

## Source Provenance

Input manifest SHA-256: 4224a9f99a6a5e77125f30a8f3f43fc692e6bcfde35eabe74ae68e4ac92d38e8

Input source schema: gentle.tss_fasta_bundle.v1

Input source revision: not_assessed

Input dataset ID: not_assessed

Bundle-source producer SHA-256 (source-declared): not_assessed

Separately declared annotation release: synthetic-annotation-v1

Profile producer executable SHA-256: 4c7fe3968561168083bb54a7d09f4b4657857cb8a76f555b3d95da850f6703e4

The named receipt input_manifest_sha256 is copied from the report's explicit bundle_manifest
input binding (legacy manifest role also accepted), which binds the original manifest bytes.
When report.source is present its manifest_sha256 must match that binding exactly; its
source_revision is retained separately from the scoring producer and exporter build revisions.
A legacy report without source metadata still binds its input manifest, but source_revision
remains null. Missing optional provenance is not_assessed (JSON/TSV null), not a guessed value.
In particular, annotation release is never extracted from words in a genome label.
The optional producer_executable receipt input binds the supplied profile producer binary
digest; bundle_source_producer binds the separate source-declared producer digest when present.
Neither substitutes for receipt.executable_sha256, which hashes the running exporter. These
identities are copied from the report, not independently reverified against original input files.

## Files And Coordinates

report.json is the complete computational report. Gene JSON files are protocol-compatible subsets,
including their comparisons, original raw scores, exact panel, normalization metadata, source
provenance, optional producer binary digest, selection evidence and input bindings.
Each resolved_matrix digest is lowercase SHA-256 of serde_json::to_vec of the tuple
(accession, Some(exact factor name), matrix_counts), serialized as [accession,name,counts].
Counts retain A/C/G/T column order and f64 JSON serialization; this is not a counts-only digest.
index.json lists genes, TSS memberships, pages and data-file SHA-256 values. Gene order and TSS order
follow the supplied report; filenames combine safe symbol/ID fragments and a full identity digest.
The index repeats named source provenance and selected-panel descriptions for cross-gene inspection.
comparisons.tsv is the cross-gene table; each gene also has its own comparison TSV.

Score TSV is long-form: one row per sequence-base window start, matrix and local motif strand.
local_window_start_0based is a start, never a motif center or an inferred motif 5-prime endpoint.
tss_relative_window_start_bp and genomic_window_start_1based use protocol TssGeometry methods.
The common start coordinate is identical for both local motif strands. Genomic positions decrease
along a minus-transcript window. Motif intervals are inclusive, 1-based and stored ascending;
motif_5prime_genomic_1based is the opposite interval endpoint for the local reverse motif.
Both local and genomic motif strands are explicit. An unavailable raw score is literal null, not zero.
Non-fitting starts (including the trailing motif-length-minus-one bases) have null scores and motif
interval/endpoint fields, with availability=outside_window. Internal nulls remain unavailable.
display_score alone applies panel.clip_negative (max(raw_score, 0)); raw_score is unchanged.
Score units and all computation policies are in the panel and score_policy, not inferred from heights.
Comparison values are copied, not recomputed: their method retains strand, smoothing and raw-input
policy; null correlations carry undefined_reason. Display clipping does not alter these values.
Every unordered same-factor matrix pair must have one comparison on each local strand. Paired
counts must match starts where both raw vectors are available; excluded counts cover the longer
vector's start domain minus paired starts. Validation checks these masks, not correlation formulas.
Every TSV begins with # comment lines carrying canonical non-claims and coordinate/method policies.
Provenance and selected-panel legends are also included in those comments. Selection columns
preserve the descriptive label, legend, criterion and factor; unselected rows use null for all four.
Skip these lines before reading the column header. TSV text escapes backslash, tab, CR and LF
as \\, \t, \r and \n respectively; a leading # in a text cell is escaped as \#.

## Rendering And Replay

SVG pages use gentle_render::tss_profiles::render_tss_profile_pages. Optional PNG uses in-process
resvg at scale 1 with metadata stripping disabled. PDF is a single-page raster-backed RGB image,
not a vector PDF. Its FlateDecode/zlib compression preserves every RGB pixel and the resolution.
Receipt metadata records that encoding, actual dimensions, available font-face counts and the
root's locked resvg version. Audited PNG/PDF helpers inspect positioned glyphs in the same parsed
rendering tree, including fallback and nested SVG trees. Each font_identities entry records its
family names, PostScript name, face_index and lowercase SHA-256 of the complete font source/container
bytes. The face index is recorded separately, not mixed into that digest. used_font_face_count counts
these distinct used source/face identities; font_face_count counts available faces, not used faces.
An unbound glyph-used face fails export instead of leaving a partial audit. This is selected
layout-font evidence, not proof that every glyph produced visible pixels or that glyph coverage is
complete. Font files and host paths are not copied into the export. For font-matched replay, supply
the exact recorded font bytes and compare the newly selected identities, backend and options.
Font pinning is not enforced; cross-host byte-identical replay is NOT verified by recording these
identities alone. SVG-only pages retain font-family declarations and viewer-dependent font selection;
their font_identities remain null. PNG/PDF font bindings do not pin how another viewer renders SVG.

For report-only replay, load report.json as TssProfileReport and export-request.json as
ExportTssProfilesRequest. Replace output_dir (stored portably as '.') with a fresh output directory,
then invoke export_tss_profiles through the shared export operation. No scores are recomputed.
To recompute instead, separately recover the exact hash-bound manifest, FASTA, panel, optional
selection and registry sources. The original inputs are not copied into this export. Use the
producer revision and producer lockfile digest in report.json, not the exporter revision as a substitute.

Producer revision: 0.1.0-internal.10+git.8b7b6b472c0d2df3177b31901ceab69932ab90b1
Exporter revision: 0.1.0-internal.10+git.8b7b6b472c0d2df3177b31901ceab69932ab90b1
Receipt executable_sha256 hashes the running exporter executable. Receipt lockfile_sha256 hashes
the exporter build's embedded Cargo.lock; index.json separately retains the producer lockfile digest.
No timestamps or staging/host paths enter the result identity. report_sha256 hashes the exact
compact report.json bytes. The receipt binds all published files except receipt.json itself.
The index hashes neither itself nor the receipt, avoiding a circular hash. Verification checks
the complete inventory, file hashes and report bindings; it is not a digital signature or proof
of authenticity against an attacker who can rewrite both files and the receipt.

## Selected TSS Descriptions

Selection labels, legends, criteria and optional factors below are descriptive report metadata,
not an independent evidence assessment. Selection does not establish TSS usage, direct binding
or promoter activity. A selected TSS without recorded evidence gets only a generic label and
legend; its criterion and factor remain not_assessed. No CUT&RUN support is inferred.

### SYNMINUS / synthetic-gene-minus / synthetic-minus

Selected-panel label: Selected TSS

Legend: Selected according to the supplied report. No descriptive evidence criterion was recorded; this does not establish TSS usage, direct binding or promoter activity.

Recorded criterion: not_assessed

Recorded factor: not_assessed

