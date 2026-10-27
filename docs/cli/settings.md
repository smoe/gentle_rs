# Engine Settings

<!-- manual-counterparts:start -->
Related GUI reference: [Configuration Window](../gui/configuration-window.md) | [Engine Settings (Engine Ops)](../gui/engine-settings-engine-ops.md).
These topic links do not imply identical feature coverage across interfaces.
<!-- manual-counterparts:end -->

- Shared shell:
  - Commands:
    - `set-param NAME JSON_VALUE`
      - anchor-related examples:
        - `set-param require_verified_genome_anchor_for_extension true`
        - `set-param genome_anchor_prepared_fallback_policy "single_compatible|always_explicit|off"`
    - promoter-design TF score-track examples:
      - structured summary through the raw operation bridge:
        - `gentle_cli op '{"type":"SummarizeTfbsScoreTracks","target":{"kind":"seq_id","seq_id":"tp73_context","span_start_0based":0,"span_end_0based_exclusive":4000},"motifs":["TP73","SP1","BACH2","PATZ1"],"score_kind":"llr_background_tail_log10","clip_negative":true}'`
      - shared SVG export through the shell/CLI route:
        - `gentle_cli shell 'features tfbs-score-tracks-svg tp73_context docs/figures/tp73_upstream_tfbs_score_tracks.svg --motif TP53 --motif TP63 --motif TP73 --motif PATZ1 --motif SP1 --motif BACH2 --motif REST --range 15564..16764 --score-kind llr_background_tail_log10 --allow-negative'`
      - state-optional inline export through the same shell/CLI route:
        - `gentle_cli shell 'features tfbs-score-tracks-svg --sequence-text ACGTACGTACGT --id-hint inline_promoter --output /tmp/inline_tfbs_tracks.svg --motif SP1 --range 0..12 --score-kind llr_background_tail_log10'`
      - anchor-vs-candidate similarity ranking through the same shell/CLI
        family:
        - `gentle_cli shell 'features tfbs-track-similarity tp73_context --anchor-motif TP73 --candidate-motif ALL --range 15564..16764 --ranking-metric smoothed_spearman --score-kind llr_background_tail_log10 --species "Homo sapiens" --include-remote-metadata --limit 25 --path /tmp/tp73_tfbs_similarity.json'`
      - promoter evidence ledger for later GUI/ClawBio interpretation:
        - `gentle_cli shell 'features promoter-evidence-matrix tp73_context --gene-label TP73 --promoter-upstream-bp 1000 --promoter-downstream-bp 200 --path /tmp/tp73_promoter_evidence.json'`
      - isoform-aware promoter comparison and expression evidence assignment:
        - `gentle_cli shell 'features promoter-isoform-comparison tp73_context --gene-label TP73 --promoter-upstream-bp 1000 --promoter-downstream-bp 200 --path /tmp/tp73_isoform_promoters.json'`
        - `gentle_cli shell 'features promoter-expression-evidence tp73_context --gene-label TP73 --expression-json {"transcript_id":"ENST...","sample_id":"case_1","condition":"case","value":12.4,"unit":"TPM"} --source-label rna_seq_table_1 --path /tmp/tp73_promoter_expression.json'`
      - dynamic promoter workflow example:
        - `gentle_cli shell 'resources resolve-tf-query "Yamanaka factors" OCT4 "KLF family" --output /tmp/tf_queries.json'`
        - `gentle_cli genomes extract-promoter "Human GRCh38 Ensembl 116" TERT --output-id grch38_tert_promoter --upstream-bp 1000 --downstream-bp 200`
        - `gentle_cli shell 'features tfbs-scan grch38_tert_promoter --motif "Yamanaka factors" --motif SP1 --min-llr-quantile 0.95 --max-hits 100 --path /tmp/tert_tfbs_scan.json'`
        - `gentle_cli shell 'features tfbs-score-tracks-svg grch38_tert_promoter /tmp/tert_yamanaka_sp1_tracks.svg --motif "Yamanaka factors" --motif SP1 --score-kind llr_background_tail_log10'`
        - `gentle_cli shell 'features tfbs-track-similarity grch38_tert_promoter --anchor-motif SP1 --candidate-motif "Yamanaka factors" --ranking-metric smoothed_spearman --score-kind llr_background_tail_log10 --species "Homo sapiens" --include-remote-metadata --limit 25 --path /tmp/tert_sp1_yamanaka_similarity.json'`
      - multi-gene promoter comparison through the same reference shell family:
        - `gentle_cli shell 'genomes promoter-tfbs-summary "Human GRCh38 Ensembl 116" --gene TERT --gene TP73 --motif "Yamanaka factors" --motif SP1 --upstream-bp 1000 --downstream-bp 200 --score-kind llr_background_tail_log10 --path /tmp/tert_tp73_promoter_tfbs.summary.json'`
        - `gentle_cli shell 'genomes promoter-cohort-comparison "Human GRCh38 Ensembl 116" --cohort-label tert_tp73_review --cohort-kind manual --gene TERT --gene TP73 --motif "Yamanaka factors" --motif SP1 --expression-json {"gene_label":"TP73","sample_id":"case_1","condition":"case","value":12.4,"unit":"TPM"} --expression-source-label reviewed_table --path /tmp/tert_tp73_promoter_cohort.json'`
        - `gentle_cli shell 'genomes promoter-tfbs-svg "Human GRCh38 Ensembl 116" --gene TERT --gene TP73 --motif "Yamanaka factors" --motif SP1 --upstream-bp 1000 --downstream-bp 200 --score-kind llr_background_tail_log10 /tmp/tert_tp73_promoter_tfbs.svg'`
      - `--score-kind` now lets you inspect raw bit scores, in-window
        quantiles, random-background percentiles, or `-log10(background tail)`
        views; the background-tail modes are the best first pass when random
        sequence already produces visually noisy low-level activity
      - clip-negative mode remains the intended default for presentation-oriented
        promoter plots on the bit-based score kinds because it suppresses
        negative-only windows and leaves only positive motif support;
        `--allow-negative` is available when a figure should show the full
        continuous score landscape even when only a subset of factors crosses
        into positive support

<!-- manual-page-nav:start -->

---

[CLI Manual](../cli.md) | [Previous: Services and Resource Catalogs](catalogs.md) | [Next: Latency Diagnostics](diagnostics.md)
<!-- manual-page-nav:end -->
