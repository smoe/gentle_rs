# Dotplots, Flexibility and Alignment

<!-- manual-counterparts:start -->
Related GUI reference: [Map Modes And Dotplot Workspace](../gui/map-modes-and-dotplot-workspace.md).
These topic links do not imply identical feature coverage across interfaces.
<!-- manual-counterparts:end -->

Dotplot/flexibility capability status:

- `gentle_cli`: supported via shared-shell/direct commands:
  - `dotplot compute|overlay-compute|list|show|render-svg`
  - `render-dotplot-svg`
  - `transcripts derive`
  - `flex compute|list|show`
  - `splicing-refs derive`
  - `align compute`
  backed by `ComputeDotplot`, `ComputeDotplotOverlay`, `DeriveTranscriptSequences`,
  `ComputeFlexibilityTrack`, `DeriveSplicingReferences`,
  `AlignSequences`, and `RenderDotplotSvg`.
  - `dotplot compute` supports self and pairwise modes via
    `--mode self_forward|self_reverse_complement|pair_forward|pair_reverse_complement`
    with optional `--reference-seq`, `--ref-start`, and `--ref-end`.
  - `dotplot overlay-compute` accepts repeated `--query-spec JSON_OR_@FILE`
    payloads for shared-reference multi-query overlays; if `--reference-seq` is
    omitted, the owner sequence id is reused as the shared reference.
  - `dotplot render-svg` / `render-dotplot-svg` accept
    `--overlay-x-axis percent_length|left_aligned_bp|right_aligned_bp|shared_exon_anchor`
    so overlay exports can switch between normalized transcript length,
    left-aligned bp, right-aligned bp, and shared-exon-anchor layouts
    without recomputing the dotplot.
  - `transcripts derive` supports:
    - full-sequence derivation (all `mRNA`/`transcript` features)
    - selected feature derivation (`--feature-id`)
    - splicing-scope constrained derivation from one seed feature
      (`--scope ...` with exactly one `--feature-id`)
    - additive synthetic local `CDS` + translated-protein qualifiers when the
      admitted transcript has resolvable CDS context and translation-table
      metadata
  - exon-skipped isoform creation is a two-phase shared-shell route:
    - `transcripts exon-skip-plan SEQ_ID --feature-id N [--skip exon_2|START..END ...] [--overlap START..END ...] [--length-mod3 0|1|2] [--coding-mod3 0|1|2] [--coding-context utr-only|cds-only|mixed-utr-cds] [--phase-entry codon-boundary|split-codon|split-codon-1|split-codon-2|unavailable] [--feature-query-json JSON] [--plan-id ID]`
      where exon-skip `START..END` intervals are 1-based inclusive genomic
      coordinates
    - `transcripts exon-skip-materialize PLAN_ID [--candidate-id ID ...] [--output-prefix PREFIX] [--return genbank|cdna_fasta|amino_acid_sequence|amino_acid_fasta ...]`
    - phase 1 stores `gentle.exon_skip_selection_plan.v1`; phase 2 consumes the
      stored plan and creates both a genomic annotation product and retained
      exon cDNA/mRNA sequence
    - introspection exposes stored plans as `exon_skip_plan.exists(PLAN_ID)`;
      materialization readiness requires that fact, while planning can verify it
      after a deterministic `--plan-id` is supplied
    - phase-1 candidate rows include exon position (`ordinal`,
      `transcript_exon_count`, `transcript_position`), support frequency
      (`support_transcript_count`, `support_transcript_total`,
      `support_fraction`), full-exon frame hints (`length_bp`, `length_mod3`,
      `frame_neutral_length`), coding-only frame hints (`coding_skip_bp`,
      `coding_skip_mod3`, `frame_neutral_coding_skip`, `coding_context`),
      flanking intron lengths, optional `left_cds_phase` / `right_cds_phase`,
      stable `cds_phase_entry_kind`, and any `cds_phase_warning` /
      `coding_frame_note`; `--length-mod3`, `--coding-mod3`,
      `--coding-context`, and `--phase-entry` select by those engine-owned
      attributes
    - `--return` lets automation callers request a compact handoff payload in
      the materialization report, such as the adjusted GenBank entry or just
      the amino-acid sequence
    - coding hints are specific to the selected transcript, not all CDS records
      in the locus. Ambiguous associations are rejected; absent CDS context is
      `unknown`, not evidence for a UTR-only exon or a frame-neutral coding skip.
      Materialization projects the original CDS and warns when a coding end is
      lost; it does not replace an annotated fragment with a newly found ORF.
      Rebuild older stored plans before relying on their coding hints.
    - for the committed TP73 `NM_001126241.3` (DeltaNp73beta) annotation, skipping
      the terminal exon `80232..83686` removes 16 coding bases from the 1,353-bp
      CDS. The result is 1,337 coding bp within 1,571-bp spliced cDNA (including
      234 bp of retained 5' UTR), not a complete beta isoform: the native stop is
      lost and two trailing coding bases form an incomplete codon. These are
      local coordinates in `test_files/tp73.ncbi.gb`, not chromosome positions.
  - `transcripts residue-genomic-coordinates SEQ_ID RESIDUE_START [RESIDUE_END]`
    maps transcript-native protein residue(s) back to 1-based genomic codon
    nucleotide positions; `--transcript ID` narrows the query to a transcript
    id, label, or `n-N` feature id, and split codons across exon junctions are
    reported base-by-base in coding order.
- `gentle_js`: baseline support via `apply_operation` (`ComputeDotplot`,
  `DeriveTranscriptSequences`, `ComputeFlexibilityTrack`,
  `DeriveSplicingReferences`, `AlignSequences`) plus
  `render_dotplot_svg(...)` convenience wrapper over `RenderDotplotSvg`.
- `gentle_lua`: baseline support via `apply_operation` (`ComputeDotplot`,
  `DeriveTranscriptSequences`, `ComputeFlexibilityTrack`,
  `DeriveSplicingReferences`, `AlignSequences`) plus
  `render_dotplot_svg(...)` convenience wrapper over `RenderDotplotSvg`.
- `gentle_py`: baseline support via `op(...)` plus
  `render_dotplot_svg(...)` convenience wrapper over `RenderDotplotSvg`.


- Shared shell:
  - Commands:
    - `dotplot compute SEQ_ID [--reference-seq REF_SEQ_ID] [--start N] [--end N] [--ref-start N] [--ref-end N] [--mode self_forward|self_reverse_complement|pair_forward|pair_reverse_complement] [--word-size N] [--step N] [--max-mismatches N] [--tile-bp N] [--id DOTPLOT_ID]`
    - `dotplot overlay-compute OWNER_SEQ_ID [--reference-seq REF_SEQ_ID] --query-spec JSON_OR_@FILE [--query-spec JSON_OR_@FILE ...] [--ref-start N] [--ref-end N] [--word-size N] [--step N] [--max-mismatches N] [--tile-bp N] [--id DOTPLOT_ID]`
    - `dotplot list [SEQ_ID]`
    - `dotplot show DOTPLOT_ID`
    - `dotplot render-svg SEQ_ID DOTPLOT_ID OUTPUT.svg [--flex-track ID] [--display-threshold N] [--intensity-gain N] [--overlay-x-axis percent_length|left_aligned_bp|right_aligned_bp|shared_exon_anchor] [--overlay-anchor-exon START..END]`
    - `transcripts derive SEQ_ID [--feature-id N ...] [--scope all_overlapping_any_strand|target_group_any_strand|all_overlapping_target_strand|target_group_target_strand] [--output-prefix PREFIX]`
    - `transcripts residue-genomic-coordinates SEQ_ID RESIDUE_START [RESIDUE_END] [--transcript ID]`
    - `flex compute SEQ_ID [--start N] [--end N] [--model at_richness|at_skew] [--bin-bp N] [--smoothing-bp N] [--id TRACK_ID]`
    - `flex list [SEQ_ID]`
    - `flex show TRACK_ID`
    - `splicing-refs derive SEQ_ID START_0BASED END_0BASED [--seed-feature-id N] [--scope all_overlapping_any_strand|target_group_any_strand|all_overlapping_target_strand|target_group_target_strand] [--output-prefix PREFIX]`
    - `align compute QUERY_SEQ_ID TARGET_SEQ_ID or --query-sequence-text DNA --target-sequence-text DNA [--query-id-hint TEXT] [--target-id-hint TEXT] [--query-range START..END|--query-start N --query-end N] [--target-range START..END|--target-start N --target-end N] [--mode global|local] [--match N] [--mismatch N] [--gap-open N] [--gap-extend N]`
    - `reads acquire status MANIFEST.tsv --cache-dir DIR --work-dir DIR`
    - `reads acquire prepare MANIFEST.tsv --cache-dir DIR --work-dir DIR [--analysis-format fasta|fastq] [--read-layout single_end|paired_end|split_spot] [--threads N] [--max-size SIZE] [--min-free-gb N] [--drop-intermediate-fastq] [--continue-on-error]`
    - `reads acquire inspect RUN_ACCESSION --cache-dir DIR --work-dir DIR`
    - `reads acquire cancel RUN_ACCESSION --cache-dir DIR --work-dir DIR`

<!-- manual-page-nav:start -->

---

[CLI Manual](../cli.md) | [Previous: Feature Annotation](features.md) | [Next: Figure and Layout Exports](rendering.md)
<!-- manual-page-nav:end -->
