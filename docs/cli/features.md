# Feature Annotation

<!-- manual-counterparts:start -->
Related GUI reference: [Strand Placement](../gui/strand-placement.md) | [Feature Tree And Selection](../gui/feature-tree-and-selection.md).
These topic links do not imply identical feature coverage across interfaces.
<!-- manual-counterparts:end -->

- Shared shell:
  - Commands:
    - `features formula SEQ_ID EXPR`
    - `features edit-location SEQ_ID FEATURE_INDEX [--segment-index INDEX] --start-1based N --end-1based-inclusive M [--dry-run] [--expected-feature-fingerprint-sha256 SHA] [--path OUT.json]`
    - `features create SEQ_ID --kind KIND --start-1based N --end-1based-inclusive M [--strand forward|reverse] [--qualifier KEY[=VALUE] ...] [--dry-run] [--expected-annotation-state-fingerprint-sha256 SHA] [--path OUT.json]`
    - `features delete SEQ_ID FEATURE_INDEX [--dry-run] [--expected-feature-fingerprint-sha256 SHA] [--expected-annotation-state-fingerprint-sha256 SHA] [--path OUT.json]`
    - `features split SEQ_ID FEATURE_INDEX --split-before-1based N [--dry-run] [--expected-feature-fingerprint-sha256 SHA] [--expected-annotation-state-fingerprint-sha256 SHA] [--path OUT.json]`
    - `features merge SEQ_ID FIRST_FEATURE_INDEX SECOND_FEATURE_INDEX [--dry-run] [--expected-first-feature-fingerprint-sha256 SHA] [--expected-second-feature-fingerprint-sha256 SHA] [--expected-annotation-state-fingerprint-sha256 SHA] [--path OUT.json]`
    - `features query SEQ_ID [--kind KIND] [--kind-not KIND] [--range START..END|--start N --end N] [--overlap|--within|--contains] [--strand any|forward|reverse] [--label TEXT] [--label-regex REGEX] [--qual KEY] [--qual-contains KEY=VALUE] [--qual-regex KEY=REGEX] [--min-len N] [--max-len N] [--nearest-to POSITION] [--limit N] [--offset N] [--sort feature_id|start|end|kind|length] [--desc] [--include-source] [--include-qualifiers]`
    - `features export-bed SEQ_ID OUTPUT.bed [--coordinate-mode auto|local|genomic] [--include-restriction-sites] [--restriction-enzyme NAME] [--kind KIND] [--kind-not KIND] [--range START..END|--start N --end N] [--overlap|--within|--contains] [--strand any|forward|reverse] [--label TEXT] [--label-regex REGEX] [--qual KEY] [--qual-contains KEY=VALUE] [--qual-regex KEY=REGEX] [--min-len N] [--max-len N] [--nearest-to POSITION] [--limit N] [--offset N] [--sort feature_id|start|end|kind|length] [--desc] [--include-source] [--include-qualifiers]`
    - `features tfbs-summary SEQ_ID --focus START..END [--context START..END] [--min-focus-count N] [--min-context-count N] [--limit N]`
    - `features genomic-motif-evidence [SEQ_ID|--region-set ID] --motif JASPAR_ID [--motif JASPAR_ID ...|--motifs CSV] [--range START..END] [--region [ID=]CHR:START..END ...] [--package DIR] [--database FILE] [--duckdb FILE] [--genome-id ID] [--min-score VALUE] [--min-pwm-relative-score VALUE] [--max-rows N] [--max-payload-files N] [--timeout-seconds N] [--path FILE.json]`
    - Regulatory/TSS subset packages additionally support `features genomic-motif-evidence --inspect [--search TEXT] [--catalog-offset N] [--catalog-limit N] --package DIR [--path FILE.json]`, or `--gene EXACT_ID_OR_NAME` / `--tss-id ID` in place of the sequence/region target, with 1-64 explicit motifs. See the [offline walkthrough](../regulatory_motif_subset.md).
    - `features repeat-query GENOME_ID --rmsk PATH [--rep-class CLASS] [--rep-family FAMILY] [--rep-name NAME] [--alias ALIAS] [--chromosome CHR] [--range START..END] [--limit N] [--path FILE.json]`
    - `features repeat-overlaps SEQ_ID --index RMSK_INTERVAL_INDEX.json [--range START..END] [--limit N] [--path FILE.json]`
    - `features materialize-repeats SEQ_ID --index RMSK_INTERVAL_INDEX.json [--max-features N] [--append] [--path FILE.json]`
    - `features encode-ccre-overlaps SEQ_ID --index INDEX.json [--bed BED] [--range START..END] [--class pELS|dELS] [--limit N] [--path FILE.json]`
    - `features materialize-encode-ccres SEQ_ID --index INDEX.json [--bed BED] [--class pELS|dELS] [--max-features N] [--append] [--path FILE.json]`
    - `features ensembl-regulation-overlaps SEQ_ID --index INDEX.json [--intervals PATH] [--range START..END] [--feature-type TYPE] [--limit N] [--path FILE.json]`
    - `features materialize-ensembl-regulation SEQ_ID --index INDEX.json [--intervals PATH] [--feature-type TYPE] [--max-features N] [--append] [--path FILE.json]`
    - `features repeat-cohort GENOME_ID --rmsk PATH [--rep-class CLASS] [--rep-family FAMILY] [--rep-name NAME] [--alias ALIAS] [--geometry repeat_midpoint|transcript_5utr_start|pol2_promoter_upstream|cds_stop_context] [--upstream-bp N] [--downstream-bp N] [--limit N] [--catalog PATH] [--cache DIR] [--path FILE.json]`
    - `features window-cohort-tfbs COHORT_JSON --motif TOKEN [--motif TOKEN ...] [--motifs CSV] [--score-kind llr_bits|llr_quantile|llr_background_quantile|llr_background_tail_log10|true_log_odds_bits|true_log_odds_quantile|true_log_odds_background_quantile|true_log_odds_background_tail_log10] [--allow-negative] [--catalog PATH] [--cache DIR] [--path FILE.json]`
    - `features promoter-evidence-matrix SEQ_ID [--gene-label LABEL] [--transcript-id ID] [--promoter-upstream-bp N] [--promoter-downstream-bp N] [--no-feature-overlaps] [--path FILE.json]`
    - `features promoter-isoform-comparison SEQ_ID [--gene-label LABEL] [--transcript-id ID] [--promoter-upstream-bp N] [--promoter-downstream-bp N] [--no-feature-overlaps] [--path FILE.json]`
    - `features promoter-expression-evidence SEQ_ID [--gene-label LABEL] [--transcript-id ID] [--promoter-upstream-bp N] [--promoter-downstream-bp N] [--expression-json JSON] [--source-label LABEL] [--path FILE.json]`
    - `features promoter-artifact-manifest SEQ_ID --artifact-json JSON --path FILE.json [--gene-label LABEL]`

- Shared shell:
  - Commands:
    - `features tfbs-score-tracks-svg SEQ_ID OUTPUT.svg --motif TOKEN [--motif TOKEN ...] [--motifs CSV] [--range START..END|--start N --end N] [--score-kind llr_bits|llr_quantile|llr_background_quantile|llr_background_tail_log10|true_log_odds_bits|true_log_odds_quantile|true_log_odds_background_quantile|true_log_odds_background_tail_log10] [--allow-negative]`
    - `features tfbs-score-tracks-svg --sequence-text DNA --output OUTPUT.svg [--topology linear|circular] [--id-hint TEXT] --motif TOKEN [--motif TOKEN ...] [--motifs CSV] [--range START..END|--start N --end N] [--score-kind llr_bits|llr_quantile|llr_background_quantile|llr_background_tail_log10|true_log_odds_bits|true_log_odds_quantile|true_log_odds_background_quantile|true_log_odds_background_tail_log10] [--allow-negative]`
    - `features tfbs-track-similarity SEQ_ID --anchor-motif TOKEN [--candidate-motif TOKEN ...|--candidate-motifs CSV|--candidate-motif ALL] [--range START..END|--start N --end N] [--ranking-metric raw_pearson|smoothed_pearson|raw_spearman|smoothed_spearman] [--score-kind llr_bits|llr_quantile|llr_background_quantile|llr_background_tail_log10|true_log_odds_bits|true_log_odds_quantile|true_log_odds_background_quantile|true_log_odds_background_tail_log10] [--allow-negative] [--species TEXT] [--include-remote-metadata] [--limit N] [--path FILE.json]`
    - `features tfbs-track-similarity --sequence-text DNA [--topology linear|circular] [--id-hint TEXT] --anchor-motif TOKEN [--candidate-motif TOKEN ...|--candidate-motifs CSV|--candidate-motif ALL] [--range START..END|--start N --end N] [--ranking-metric raw_pearson|smoothed_pearson|raw_spearman|smoothed_spearman] [--score-kind llr_bits|llr_quantile|llr_background_quantile|llr_background_tail_log10|true_log_odds_bits|true_log_odds_quantile|true_log_odds_background_quantile|true_log_odds_background_tail_log10] [--allow-negative] [--species TEXT] [--include-remote-metadata] [--limit N] [--path FILE.json]`
      - the background-normalized score kinds suppress everything below the
        `0.95` modeled random-background quantile and can show either the
        surviving mid-rank percentile itself or `-log10(background tail)`
      - those background-normalized views are no longer limited by the old
        finite-sample ceiling alone: GENtle now calibrates them against a
        quantized IID random-DNA window model while still carrying a larger
        deterministic random sample for summary statistics
      - the shared JSON payload now also includes `top_peaks` per motif so
        agents can reason over a short ranked list of unusual windows instead
      - the same JSON payload now also includes pairwise raw/smoothed
        correlation summaries and primary-peak offsets so agent/CLI workflows
        can judge whether strong motif neighborhoods are actually synchronized
      - motif tokens accept exact ids/names, aliases like `OCT4`, built-in
        groups like `Yamanaka factors` / `stemness`, and family-like queries
        like `KLF family`
      - omitting `--range` uses the full selected/stored sequence span
      - stacked SVG labels now render as `TF name (JASPAR id)` when possible,
        and TSS context uses one shared dashed line with a single top-level
        kinked arrow if feature-derived or promoter-provenance TSS metadata is
        available
    - `features tfbs-score-track-correlation-svg SEQ_ID OUTPUT.svg --motif TOKEN [--motif TOKEN ...] [--motifs CSV] [--range START..END|--start N --end N] [--score-kind llr_bits|llr_quantile|llr_background_quantile|llr_background_tail_log10|true_log_odds_bits|true_log_odds_quantile|true_log_odds_background_quantile|true_log_odds_background_tail_log10] [--correlation-metric pearson|spearman] [--signal-source max_strands|forward_only|reverse_only] [--allow-negative]`
      - `max_strands` now renders one strand-expanded all-vs-all matrix:
        every motif appears twice on each axis, ordered `F` then `R`, so each
        TF pair naturally reads as one adjacent `F-F / F-R / R-F / R-R` block
      - `forward_only` and `reverse_only` keep the single-signal projections
      - renders the same shared synchrony sidecar as a dedicated dual heatmap:
        smoothed selected-metric correlation on the left, raw selected-metric
        correlation on the right
      - `--correlation-metric spearman` is often the better promoter-reading
        choice when the displayed score tracks are clearly non-normal or
        heavily tied by clipping/background thresholding
      - the footer ranks the strongest synchronized TF pairs and keeps signed
        primary-peak offsets visible for quick promoter interpretation
    - `features tfbs-scan SEQ_ID --motif TOKEN [--motif TOKEN ...] [--motifs CSV] [--range START..END|--start N --end N] [--min-llr-bits VALUE] [--min-llr-quantile VALUE] [--per-tf-min-llr-bits TF=VALUE] [--per-tf-min-llr-quantile TF=VALUE] [--max-hits N] [--path FILE.json]`
    - `features tfbs-scan --sequence-text DNA [--topology linear|circular] [--id-hint TEXT] --motif TOKEN [--motif TOKEN ...] [--motifs CSV] [--range START..END|--start N --end N] [--min-llr-bits VALUE] [--min-llr-quantile VALUE] [--per-tf-min-llr-bits TF=VALUE] [--per-tf-min-llr-quantile TF=VALUE] [--max-hits N] [--path FILE.json]`
      - the same motif-token expansion rules apply here too, so promoter scans
        can use aliases, family-style queries, or functional groups
    - `features restriction-scan SEQ_ID [--range START..END|--start N --end N] [--enzyme NAME] [--max-sites-per-enzyme N] [--no-cut-geometry] [--path FILE.json]`
    - `features restriction-scan --sequence-text DNA [--topology linear|circular] [--id-hint TEXT] [--range START..END|--start N --end N] [--enzyme NAME] [--max-sites-per-enzyme N] [--no-cut-geometry] [--path FILE.json]`

- Shared shell:
  - Commands:
    - Feature formula helper notes (`features formula`):
      - non-mutating structured result schema:
        `gentle.feature_coordinate_formula_resolution.v1`
      - resolves the same coordinate/range formulas used by GUI selection and
        primer/qPCR ROI fields, including `KIND.tss` and strand-aware
        `KIND.upstream(N)`
      - range expressions return normalized `start_0based` /
        `end_0based_exclusive` intervals, so reverse-strand upstream/TSS
        formulas are directly usable by agents
    - Feature-location edit notes (`features edit-location`):
      - use `--dry-run` first; applying requires the returned
        `before_feature_fingerprint_sha256`
      - without `--segment-index`, accepts an exact simple GenBank range and
        preserves its strand
      - `--segment-index` selects one zero-based child in an exact flat
        `Join`/`Order`, optionally under one outer `Complement`; the operation
        preserves operator kind, wrapper, child count/order, qualifiers, and
        every unedited child
      - nested, fuzzy, between-base, external, gap, bond, one-of,
        non-monotonic, and circular cross-origin compounds remain read-only
      - compound reports distinguish stored child index from biological
        5-prime-to-3-prime segment number and emit review warnings for overlap,
        changed stored direction, or a CDS length delta not divisible by three
      - the report lists boundary-sharing annotations for review but modifies
        only the selected feature
      - a successful apply is a normal undoable project mutation
    - Feature-record curation notes (`features create`, `features delete`,
      `features split`, `features merge`):
      - always preview first; apply uses the returned annotation-state lock,
        while Delete/Split additionally lock their selected feature and Merge
        locks both selected features
      - Create appends one exact forward or reverse range; repeat
        `--qualifier` in the required stored order and omit `=VALUE` for a
        valueless qualifier
      - Delete accepts any existing location shape and reports the complete
        deleted record; later feature indices shift down by one
      - Split accepts only an exact simple `Range` or `Complement(Range)`;
        `--split-before-1based N` creates genomic-left and genomic-right
        records at the original index and the next index, preserving kind,
        strand, and ordered qualifiers on both explicit outputs
      - Merge accepts only two exactly touching simple records with identical
        kind, strand, and ordered qualifiers; it keeps the lower table index
        and removes the higher one. Gaps, overlaps, compounds, fuzzy
        boundaries, and metadata conflicts are errors rather than guessed
        repairs
      - overlap and matching `locus_tag`, `gene`, `protein_id`, or
        `transcript_id` values are informational review candidates, not
        inferred dependencies or automatic propagation
      - every successful apply uses the normal undo/redo history
    - Feature query helper notes (`features query`):
      - non-mutating structured result schema:
        `gentle.sequence_feature_query_result.v1`
      - deterministic filters over feature kind, range relation, strand, labels,
        qualifiers, and length
      - deterministic ordering (`feature_id|start|end|kind|length`) with
        `offset`/`limit` paging suitable for agent iteration
      - `--nearest-to POSITION` accepts a 0-based local sequence coordinate,
        orders matching intervals by nearest-base distance, and adds
        `distance_to_query_bp` to each returned row; stable ties use feature
        start, end, and id
      - nearest projected microarray interval example:
        `features query grch38_tp73 --qual-contains gentle_track_source=Array --nearest-to 42000 --limit 1 --include-qualifiers`
    - Feature BED export helper notes (`features export-bed`):
      - non-mutating structured result schema:
        `gentle.sequence_feature_bed_export.v1`
      - writes BED6 plus deterministic extra columns:
        `kind`, `row_id`, `coordinate_source`, `qualifiers_json`
      - exports all matching rows by default when `--limit` is omitted; use
        `--offset`/`--limit` only when you want a bounded subset
      - `--coordinate-mode auto` prefers genomic BED coordinates when features
        carry `chromosome` + `genomic_start_1based` +
        `genomic_end_1based`; otherwise it falls back to local `SEQ_ID`
        coordinates
      - `--include-restriction-sites` appends deterministic REBASE-derived
        `restriction_site` rows, and `--restriction-enzyme NAME` can be
        repeated to keep only selected enzymes

<!-- manual-page-nav:start -->

---

[CLI Manual](../cli.md) | [Previous: Macros and Routines](macros-routines.md) | [Next: Dotplots, Flexibility and Alignment](dotplot.md)
<!-- manual-page-nav:end -->
