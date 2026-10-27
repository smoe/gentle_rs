# Tracks and Arrays

<!-- manual-counterparts:start -->
Related GUI reference: [Array Tracks In DNA Viewer](../gui/array-tracks-in-dna-viewer.md).
These topic links do not imply identical feature coverage across interfaces.
<!-- manual-counterparts:end -->

- Shared shell:
  - Commands:
    - `helpers list [--catalog PATH] [--filter TEXT]`
    - `helpers ensembl-available [--collection all|vertebrates|metazoa] [--filter TEXT]`
    - `helpers install-ensembl SPECIES_DIR [--collection vertebrates|metazoa] [--catalog PATH] [--output-catalog PATH] [--genome-id ID] [--cache-dir PATH] [--timeout-secs N]`
    - `helpers validate-catalog [--catalog PATH]`
    - `hosts list [--catalog PATH] [--filter TEXT]`
    - `helpers status HELPER_ID [--catalog PATH] [--cache-dir PATH]`
    - `helpers genes HELPER_ID [--catalog PATH] [--cache-dir PATH] [--filter REGEX] [--biotype NAME] [--limit N] [--offset N]`
    - `helpers prepare HELPER_ID [--catalog PATH] [--cache-dir PATH] [--timeout-secs N]`
    - `helpers blast HELPER_ID QUERY_SEQUENCE [--max-hits N] [--task blastn-short|blastn] [--options-json JSON_OR_@FILE | --options-file PATH] [--catalog PATH] [--cache-dir PATH]`
    - `helpers blast-start HELPER_ID QUERY_SEQUENCE [--max-hits N] [--task blastn-short|blastn] [--options-json JSON_OR_@FILE | --options-file PATH] [--catalog PATH] [--cache-dir PATH]`
    - `helpers blast-status JOB_ID [--with-report]`
    - `helpers blast-cancel JOB_ID`
    - `helpers blast-list`
    - `helpers blast-track HELPER_ID QUERY_SEQUENCE TARGET_SEQ_ID [--max-hits N] [--task blastn-short|blastn] [--options-json JSON_OR_@FILE | --options-file PATH] [--track-name NAME] [--clear-existing] [--catalog PATH] [--cache-dir PATH]`
    - `helpers extract-region HELPER_ID CHR START END [--output-id ID] [--annotation-scope none|core|full] [--max-annotation-features N] [--include-genomic-annotation|--no-include-genomic-annotation] [--rmsk-index PATH] [--max-repeat-features N] [--append-repeat-features] [--catalog PATH] [--cache-dir PATH]`
    - `helpers extract-gene HELPER_ID QUERY [--occurrence N] [--output-id ID] [--extract-mode gene|coding_with_promoter] [--promoter-upstream-bp N] [--annotation-scope none|core|full] [--max-annotation-features N] [--include-genomic-annotation|--no-include-genomic-annotation] [--rmsk-index PATH] [--max-repeat-features N] [--append-repeat-features] [--catalog PATH] [--cache-dir PATH]`
    - `helpers extract-promoter HELPER_ID QUERY [--occurrence N] [--transcript-id ID] [--output-id ID] [--upstream-bp N] [--downstream-bp N] [--annotation-scope none|core|full] [--max-annotation-features N] [--include-genomic-annotation|--no-include-genomic-annotation] [--rmsk-index PATH] [--max-repeat-features N] [--append-repeat-features] [--catalog PATH] [--cache-dir PATH]`
    - `proteases list [--filter TEXT] [--output PATH]`
    - `proteases show QUERY [--output PATH]`
    - `proteases digest SEQ_ID PROTEASE[,PROTEASE...] [--output-prefix PREFIX] [--min-length-aa N] [--predict-only]`
    - `proteases digest-gel-svg SEQ_ID PROTEASE[,PROTEASE...] OUTPUT.svg [--min-length-aa N] [--ladder NAME] [--ladders CSV]`
    - `tracks import-bed SEQ_ID PATH [--name NAME] [--min-score N] [--max-score N] [--clear-existing]`
    - `tracks import-bigwig SEQ_ID PATH [--name NAME] [--min-score N] [--max-score N] [--clear-existing]`
    - `tracks import-vcf SEQ_ID PATH [--name NAME] [--min-score N] [--max-score N] [--clear-existing]`
    - `tracks tracked list`
    - `tracks tracked add PATH [--source auto|bed|bigwig|vcf] [--name NAME] [--min-score N] [--max-score N] [--clear-existing]`
    - `tracks tracked remove INDEX`
    - `tracks tracked clear`
    - `tracks tracked apply [--index N] [--only-new-anchors]`
    - `arrays inspect-microarray-track MANIFEST`
    - `arrays project-microarray-track SEQ_ID MANIFEST [--contrasts CSV] [--level probeset] [--min-abs-logfc N] [--max-adj-p N] [--max-features N] [--clear-existing]`
    - `arrays inspect-probe-region-output OUTPUT_DIR`
    - `arrays import-apt-probe-region-output SUMMARY.tsv ANNOTATION.csv OUTPUT_DIR [--metadata PATH] [--condition-column NAME] [--sample-column NAME] [--probe-intensity PATH] [--probe-id-column NAME] [--platform NAME] [--normalization NAME] [--coordinate-system ID] [--genome-build ID]`
    - `arrays run-probe-region-backend PLAN.json [--backend NAME] --allow-external-execution`
    - `arrays render-probe-region-output-svg OUTPUT_DIR OUTPUT.svg`
    - `arrays render-probe-region-evidence-svg REPORT.json OUTPUT.svg`
    - `arrays project-probe-region-output SEQ_ID OUTPUT_DIR [--contrasts CSV] [--level probe_region|pm_probe] [--min-abs-logfc N] [--max-features N] [--clear-existing]`
    - `arrays interpret-probe-region-evidence SEQ_ID [--gene LABEL] [--level all|probe_region|pm_probe] [--min-abs-logfc N] [--threshold-source TEXT] [--policy-sha256 SHA256] [--path FILE]`
    - `arrays probe-regions (--cel PATH ... | --dataset ID) (--gene SYMBOL|--genes CSV|--locus LOCUS|--loci CSV|--transcript-cluster-id ID|--probeset-id ID ...) [--metadata PATH] [--platform NAME] [--annotation-library PATH] [--r-library-path PATH ...] [--condition-column NAME] [--sample-column NAME] [--block-column NAME] [--paired-by-replicate-suffix] [--normalization rma|quantile-feature|none] [--plot] [--output DIR] [--cache-dir DIR] [--dry-run]`

## Probe-Region Preparation and Evidence

- `Rscript scripts/analyze_p73_clariomd_probe_level.R [DATASET_DIR] [OUTPUT_DIR]`
  - External-analysis helper for the Rostock `E-MTAB-14704` Clariom D CEL
    resource.
  - Requires R/Bioconductor packages `oligo`, `limma`, and
    `pd.clariom.d.human`.
  - Reads the prepared SDRF/CEL files, exports sample metadata,
    transcript-cluster RMA matrices, probeset-level RMA matrices when
    supported by the platform design package, limma contrasts, and first-pass
    probeset heterogeneity tables for splice-variant triage.
- `Rscript scripts/probe_regions_oligo.R --r-library-path .r-lib --cel sample1.CEL --cel sample2.CEL --metadata samples.tsv --gene PATZ1 --platform-package pd.clariom.d.human --coordinate-system hg38 --genome-build GRCh38 --output analysis/probe_regions`
  - Generic external R/oligo backend helper for the `arrays probe-regions`
    preflight contract.
  - Currently supports `--normalization rma`; other normalization modes remain
    plan-only until a matching backend is added.
  - Requires explicit CEL paths and at least one selector unless
    `--allow-all-features` is used intentionally.
  - Requires R/Bioconductor packages `oligo`, `limma`, `Biobase`, `DBI`,
    `RSQLite`, and the platform design package such as `pd.clariom.d.human`.
    The platform package is installed externally and normally contains its own
    `extdata/<package-name>.sqlite`; this database is not a separate GENtle
    input. Probeset output stops with a diagnostic if the selected package has
    no unambiguous SQLite database. Transcript-cluster-only output can continue
    from `netaffxTranscript.rda` with an explicit warning.
  - `--r-library-path PATH` is repeatable. It prepends explicit agent-local,
    user, or system package roots to R's `.libPaths()` for both preflight and
    execution. With no explicit flag, an existing workspace `.r-lib` remains
    the compatibility default. Missing-package errors print every effective
    library path and ask the user to check this flag when the result differs
    from the R environment in which packages were installed.
  - Writes `region_intensity_chrom_order.csv`, expression/feature TSVs,
    limma contrast TSVs when metadata defines conditions, a normalized matrix
    manifest, provenance JSON, and `sessionInfo.txt`. Supplying
    `--coordinate-system` and `--genome-build` records the coordinate basis
    needed before any later genome-anchored projection can be accepted.
    Provenance also records the R version, exact versions of `oligo`, `limma`,
    `Biobase`, `DBI`, `RSQLite`, and the platform package, plus the helper
    method version and input fingerprints. The script only checks dependencies;
    it does not install them. A direct R invocation records path, size, and
    modification time; `arrays run-probe-region-backend` additionally binds
    those inputs with Rust-computed SHA-256 fingerprints in the finalized
    provenance.
- `arrays inspect-probe-region-output analysis/probe_regions`
  - Read-only GENtle inspection of a completed `probe_regions_oligo.R` output
    directory.
  - Emits `gentle.probe_region_output_inspection.v1` JSON summarizing the
    region table, optional `sample_table.tsv`, normalized matrix manifest,
    provenance JSON, sample/condition/logFC columns, chromosomes, gene symbols,
    bounded preview rows, row counts, feature counts, coordinate/build
    declarations, declared coordinate projection maps, projection-readiness
    blockers, and required-column problems.
  - The same inspection includes `r_version`, `package_versions`,
    `analysis_method_version`, and `input_fingerprints` when newer helper
    provenance supplies them. Older output directories remain readable with
    those optional fields absent.
- `arrays import-apt-probe-region-output apt.summary.tsv annotation.csv analysis/probe_regions --metadata samples.csv --condition-column condition --sample-column file --probe-intensity probe_intensity.tsv --probe-id-column probe_id --platform Clariom_D_Human --normalization rma-sketch --coordinate-system hg38 --genome-build GRCh38`
  - Converts explicit APT summary output plus an explicit annotation/NetAffx
    coordinate table into GENtle's helper-output directory contract. The
    annotation table must provide probeset/region id, chromosome, start, and
    stop columns, with optional strand, transcript-cluster, probe-count, and
    gene-symbol columns. When optional metadata is supplied, matching sample
    rows add `mean_log2_*`, `sd_log2_*`, and default `log2FC_*` columns to the
    imported region table. If annotation rows also include `probe_id` and
    PM probe `probe_start`/`probe_stop` coordinates, GENtle writes
    `probe_intensity_chrom_order.csv`. With `--probe-intensity`, that table
    uses rows from the supplied probe-level matrix and marks them
    `probe_level_input`; without it, GENtle preserves the compatibility
    fallback by using parent-probeset summary intensities marked as
    `parent_probeset_summary`. GENtle writes
    `region_intensity_chrom_order.csv`,
    `normalized_feature_matrix_manifest.json`, and `provenance.json`; it does
    not run APT.
- `arrays render-probe-region-output-svg analysis/probe_regions analysis/probe_regions/probe_region_plot.svg`
  - Read-only deterministic SVG export for inspected helper output. The SVG
    plots `mean_log2_*` condition tracks in the upper panel and `log2FC_*`
    tracks in the lower panel, using the existing chromosome-ordered helper
    table and reporting projection blockers without running R/APT.
- `arrays run-probe-region-backend analysis/probe_regions/plan.json --allow-external-execution`
  - Explicit local execution of a persisted `gentle.probe_region_plan.v1`.
    GENtle refuses to run without `--allow-external-execution`, refuses plans
    whose preflight/backend readiness failed, captures stdout/stderr/exit
    status, validates the resulting four-file helper-output contract, and
    writes hardened `gentle.probe_region_backend_provenance.v1` provenance. It
    does not download or install CEL files, vendor resources, R packages, or
    APT.
- `arrays render-probe-region-evidence-svg analysis/probe_regions/tp73_interpretation.json analysis/probe_regions/tp73_probe_geometry.svg`
  - Read-only deterministic SVG export for
    `gentle.probe_region_evidence_interpretation.v2`. The SVG draws
    transcript lanes, only the exon ranges and junction spans present in the
    report, parent probeset spans, and PM probe intervals. It is a review-only
    transcript/exon-geometry constraint visualization; it does not infer
    isoform support, probe sequence specificity, multi-hit status, or full gene
    models.
- `arrays project-probe-region-output grch38_tp73 analysis/probe_regions --contrasts TAp73-AdGFP --level pm_probe --min-abs-logfc 0.5`
  - Projects direct-coordinate-compatible or explicitly mapped `log2FC_*`
    helper-output rows into genome-anchored array features on the target
    sequence. Default `--level probe_region` projects
    `region_intensity_chrom_order.csv`; `--level pm_probe` projects only true
    PM probe rows from `probe_intensity_chrom_order.csv` marked as
    `probe_level_input`. The route refuses projection unless the helper output
    declares `coordinate_system`/`genome_build` compatible with the target
    sequence anchor, or a `coordinate_projections[]`/`projection_maps[]`
    interval map whose source and target builds match the helper output and
    anchor.
- `arrays interpret-probe-region-evidence grch38_tp73 --gene TP73 --level pm_probe --min-abs-logfc 0.5 --threshold-source gene_isoform_assay_study_policy.min_abs_regional_effect --policy-sha256 sha256:... --path analysis/probe_regions/tp73_interpretation.json`
  - Compares projected probe/probeset-region array features with the
    sequence's transcript/exon annotations and writes
    `gentle.probe_region_evidence_interpretation.v2`. The report preserves
    shared-transcript overlaps, parent probeset context, structured
    exon/junction/transcript mappings, conservative geometry scores,
    score-basis guardrails, transcript-level review labels, probe sequence
    multi-hit-not-assessed status, and coordinate-projection ambiguity rather
    than converting array intervals directly into isoform support calls.
  - To locate the nearest projected interval first, compose this operation with
    `features query SEQ_ID --qual-contains gentle_track_source=Array
    --nearest-to POSITION --limit 1 --include-qualifiers`. Array projections
    retain genomic coordinate qualifiers; transcript/exon/junction association
    remains the responsibility of this interpretation report rather than being
    inferred by the distance query.

<!-- manual-page-nav:start -->

---

[CLI Manual](../cli.md) | [Previous: Repeats and Regulatory Annotations](repeats-regulation.md) | [Next: Multi-Gene Publication Bundles](publication-bundle.md)
<!-- manual-page-nav:end -->
