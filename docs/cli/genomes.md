# Reference Genomes

<!-- manual-counterparts:start -->
Related GUI reference: [Reference Genomes (Main Menu)](../gui/reference-genomes-main-menu.md).
These topic links do not imply identical feature coverage across interfaces.
<!-- manual-counterparts:end -->

Reference genome capability status:

- `gentle_cli`: supported via shared engine operations (`PrepareGenome`, `ExtractGenomeRegion`, `ExtractGenomeGene`) and shell-level `genomes/helpers blast`
- `gentle_js`: supported via dedicated helpers (`list_reference_genomes`, `list_reference_catalog_entries`, `list_helper_catalog_entries`, `list_helper_semantics_vocabulary`, `list_host_profile_catalog_entries`, `list_ensembl_installable_genomes`, `is_reference_genome_prepared`, `list_reference_genome_genes`, `blast_reference_genome`, `blast_helper_genome`, `prepare_genome`, `extract_genome_region`, `extract_genome_gene`) and `apply_operation`
- `gentle_lua`: supported via dedicated helpers (`list_reference_genomes`, `list_reference_catalog_entries`, `list_helper_catalog_entries`, `list_helper_semantics_vocabulary`, `list_host_profile_catalog_entries`, `list_ensembl_installable_genomes`, `is_reference_genome_prepared`, `list_reference_genome_genes`, `blast_reference_genome`, `blast_helper_genome`, `prepare_genome`, `extract_genome_region`, `extract_genome_gene`) and `apply_operation`


- Shared shell:
  - Commands:
    - `genomes list [--catalog PATH] [--filter TEXT]`
    - `genomes ensembl-available [--collection all|vertebrates|metazoa] [--filter TEXT]`
    - `genomes install-ensembl SPECIES_DIR [--collection vertebrates|metazoa] [--catalog PATH] [--output-catalog PATH] [--genome-id ID] [--cache-dir PATH] [--timeout-secs N]`
    - `genomes validate-catalog [--catalog PATH]`
    - `genomes status GENOME_ID [--catalog PATH] [--cache-dir PATH]`
    - `genomes genes GENOME_ID [--catalog PATH] [--cache-dir PATH] [--filter REGEX] [--biotype NAME] [--limit N] [--offset N]`
    - `genomes prepare GENOME_ID [--catalog PATH] [--cache-dir PATH] [--timeout-secs N]`
    - `genomes prepare-blast-resource RESOURCE_ID [--catalog PATH] [--cache-dir PATH] [--force-rebuild]`
    - `genomes import-blast-resource RESOURCE_ID [--catalog PATH] [--cache-dir PATH]` (alias)
    - `genomes adopt-blast-resource RESOURCE_ID --sequence-path PREPARED.fa --blast-db-prefix PREFIX --expected-content-fingerprint sha256:... [--expected-sequence-sha1 SHA1] [--verify-source-identities] [--catalog PATH] [--cache-dir PATH]`
    - `genomes inspect-blast-resource RESOURCE_ID [--catalog PATH] [--cache-dir PATH]`
    - `genomes blast GENOME_ID QUERY_SEQUENCE [--max-hits N] [--task blastn-short|blastn] [--options-json JSON_OR_@FILE | --options-file PATH] [--catalog PATH] [--cache-dir PATH]`
    - `genomes blast-start GENOME_ID QUERY_SEQUENCE [--max-hits N] [--task blastn-short|blastn] [--options-json JSON_OR_@FILE | --options-file PATH] [--catalog PATH] [--cache-dir PATH]`
    - `genomes blast-status JOB_ID [--with-report]`
    - `genomes blast-cancel JOB_ID`
    - `genomes blast-list`
    - `genomes blast-track GENOME_ID QUERY_SEQUENCE TARGET_SEQ_ID [--max-hits N] [--task blastn-short|blastn] [--options-json JSON_OR_@FILE | --options-file PATH] [--track-name NAME] [--clear-existing] [--catalog PATH] [--cache-dir PATH]`
    - `genomes extract-region GENOME_ID CHR START END [--output-id ID] [--annotation-scope none|core|full] [--max-annotation-features N] [--include-genomic-annotation|--no-include-genomic-annotation] [--rmsk-index PATH] [--max-repeat-features N] [--append-repeat-features] [--catalog PATH] [--cache-dir PATH]`
    - `genomes extract-gene GENOME_ID QUERY [--occurrence N] [--output-id ID] [--extract-mode gene|coding_with_promoter] [--promoter-upstream-bp N] [--annotation-scope none|core|full] [--max-annotation-features N] [--include-genomic-annotation|--no-include-genomic-annotation] [--rmsk-index PATH] [--max-repeat-features N] [--append-repeat-features] [--catalog PATH] [--cache-dir PATH]`
    - `genomes extract-promoter GENOME_ID QUERY [--occurrence N] [--transcript-id ID] [--output-id ID] [--upstream-bp N] [--downstream-bp N] [--annotation-scope none|core|full] [--max-annotation-features N] [--include-genomic-annotation|--no-include-genomic-annotation] [--rmsk-index PATH] [--max-repeat-features N] [--append-repeat-features] [--catalog PATH] [--cache-dir PATH]`

- Shared shell:
  - Commands:
    - `genomes promoter-cohort-comparison GENOME_ID --cohort-label LABEL --cohort-kind manual|co_regulated|anti_co_regulated --gene QUERY [--gene ...|--gene-json JSON] --motif TOKEN [--motif TOKEN ...|--motifs CSV] [--source-seq-id SEQ_ID ...] [--upstream-bp N] [--downstream-bp N] [--score-kind KIND] [--allow-negative] [--expression-json JSON] [--expression-source-label LABEL] [--cutrun-dataset ID ...] [--cutrun-read-report ID ...] [--catalog PATH] [--cache-dir PATH] [--path FILE.json]`

Genome convenience commands:

- `genomes list [--catalog PATH]`
  - Lists available genomes in the catalog.
- `genomes ensembl-available [--collection all|vertebrates|metazoa] [--filter TEXT]`
  - Lists Ensembl species directories that currently look installable because both FASTA and GTF listings are present.
  - Returns a read-only discovery report with current listing URLs and latest release numbers seen per collection.
  - Vertebrate discovery uses `https://ftp.ensembl.org/pub/current/fasta/` and
    `https://ftp.ensembl.org/pub/current/gtf/`; installation and spec updates pin
    files under `pub/release-N/fasta/` and `pub/release-N/gtf/`. Metazoa retains
    its separate collection root. Discovery does not rewrite existing catalogs.
- `genomes install-ensembl SPECIES_DIR [--collection vertebrates|metazoa] [--catalog PATH] [--output-catalog PATH] [--genome-id ID] [--cache-dir PATH] [--timeout-secs N]`
  - Resolves concrete current Ensembl FASTA/GTF URLs for one species directory, writes a real catalog entry, and immediately runs the normal prepare workflow.
  - When the active catalog is a writable single JSON file, GENtle updates that file in place.
  - When the active catalog comes from discovery/overlay mode, GENtle writes only the new entry into an overlay fragment instead of copying the whole merged catalog.
  - `--output-catalog PATH` can point to a catalog JSON file or to a directory where a fragment file should be created.
  - `--genome-id ID` overrides the default deterministic id derived from species + assembly + release.
- `genomes validate-catalog [--catalog PATH]`
  - Verifies catalog JSON schema/entry rules and that each entry resolves usable
    sequence/annotation source definitions.
- `genomes update-ensembl-specs [--catalog PATH] [--output-catalog PATH]`
  - Refreshes explicit pinned Ensembl URLs/specs for catalog rows that carry
    `ensembl_template` metadata.
  - Older pinned release rows stay in the catalog.
  - If a newer current species release exists, GENtle adds or refreshes the
    newest pinned row rather than overwriting the old one.
  - If the active catalog file is not writable, `--output-catalog PATH` is
    required so the updated catalog can be written as a copy.
- `genomes status GENOME_ID [--catalog PATH] [--cache-dir PATH]`
  - Shows whether the genome cache is prepared/indexed and validates its BLAST
    component with `blastdbcmd`.
  - Also reports stable lifecycle fields:
    - `resource_key`
    - `display_name`
    - `prepared`
    - `lifecycle_status = missing|running|ready|failed|cancelled|stale`
    - `current_activity` when a recent or active prepare marker exists
    - `component_ready` and `components`, including
      `components.blast_database` (`gentle.blast_database_inspection.v1`) with
      index kind, assembly/release, masking, prefix, BLAST database/tool
      versions, sequence/base counts, content fingerprint, validation status,
      warnings, and compatible operations
  - Also reports `effective_cache_dir`, `requested_catalog_key`,
    `requested_family`, and `compatible_prepared_options` so "not prepared"
    states explain exactly which cache root was inspected and whether a
    same-family prepared install already exists.
  - A valid component set is authoritative for reuse even when an older
    whole-prepare activity record is `stale`; in that case no redundant
    `prepare_command` is offered.
  - When lifecycle is `missing`, `failed`, `cancelled`, or `stale` and the
    required components are not valid, output includes a ready-to-run
    `prepare_command` hint. When lifecycle is `running`, `ready`, or components
    independently validate, that hint is suppressed.
  - On Unix, a `running` activity whose `owner_pid` no longer exists is treated
    as `stale` immediately, even if the heartbeat is younger than the regular
    stale timeout.
  - Also reports resolved source details: `sequence_source_type`,
    `annotation_source_type`, `sequence_source`, `annotation_source`.
  - Also reports optional mass metadata:
    `nucleotide_length_bp`, `molecular_mass_da`, `molecular_mass_source`.
- `genomes genes GENOME_ID [--catalog PATH] [--cache-dir PATH] [--filter REGEX] [--biotype NAME] [--limit N] [--offset N]`
  - Lists indexed genes from prepared cache (paged by `--limit`/`--offset`).
  - `--filter` is a case-insensitive regular expression.
  - `--biotype` can be repeated to constrain matches to selected biotypes.
- `genomes prepare GENOME_ID [--catalog PATH] [--cache-dir PATH] [--timeout-secs N]`
  - Runs engine `PrepareGenome`.
  - Output now includes `binary_preflight` (`gentle.blast_external_binary_preflight.v1`)
    with `makeblastdb`/`blastn` found/missing/version/path diagnostics captured
    before prepare work starts.
  - If an active prepare is reused, the human message reports the effective
    prepared-install/activity-status path when available, not only the raw
    `--cache-dir` argument. This matters for relative cache directories because
    they resolve relative to the catalog file or fragment that supplied the
    entry.
- `genomes prepare-blast-resource RESOURCE_ID [--catalog PATH] [--cache-dir PATH] [--force-rebuild]`
  - Prepares a catalog-declared, sequence-only nucleotide BLAST resource.
    `genomes import-blast-resource` is an exact alias for local collections.
  - The initial supported `blast_index_kind` is `transcriptome_cdna`. The route
    resolves `sequence_local` or `sequence_remote`, materializes the FASTA
    (converting sequence-bearing GenBank records when supplied), runs
    `makeblastdb`, validates the database with `blastdbcmd`, fingerprints the
    content, and writes a dedicated resource manifest. It also parses
    available Ensembl FASTA-defline or GenBank record identity once into a
    checksum-bound `gentle.blast_subject_annotation_index.v1` sidecar. It does
    not invent missing symbols or gene ids, nor create FASTA-index, gene-index,
    or transcript-index files.
  - A minimal catalog row is:

    ```json
    {
      "Human GRCh38 Ensembl 116 cDNA": {
        "description": "Ensembl 116 Homo sapiens whole-cDNA collection",
        "sequence_local": "/path/to/Homo_sapiens.GRCh38.cdna.all.fa.gz",
        "blast_index_kind": "transcriptome_cdna",
        "reference_name": "Homo_sapiens.GRCh38.cdna.all",
        "reference_release": "Ensembl 116",
        "blast_masking": "unmasked",
        "cache_dir": "data/genomes"
      }
    }
    ```

    Use `sequence_remote` instead of `sequence_local` for a catalogued remote
    FASTA. GENtle records the actual source and digest; users do not hand-write
    a prepared-resource manifest.
  - An existing content-bound database is protected when the catalog source or
    index content at the same prefix changes. Ordinary preparation refuses to
    accept or delete it; use adoption for independently reviewed replacement
    content, or `--force-rebuild` when GENtle should deliberately replace it.
- `genomes adopt-blast-resource RESOURCE_ID --sequence-path PREPARED.fa --blast-db-prefix PREFIX --expected-content-fingerprint sha256:... [--expected-sequence-sha1 SHA1] [--verify-source-identities] [--catalog PATH] [--cache-dir PATH]`
  - Registers corrected catalog metadata for an already-valid sequence-only
    database without invoking `makeblastdb` or altering index files.
  - The expected database fingerprint is mandatory. The report preserves the
    fingerprint algorithm (`sha256:index-file-full-or-edge-sampled-v1`) rather
    than presenting an edge-sampled large-file identity as a full content hash.
  - The prepared FASTA SHA-1 is always recorded; `prepared_sequence_verified`
    is true only when `--expected-sequence-sha1` was supplied and matched. Its
    relationship to the database is `attested` unless
    `--verify-source-identities` compares the complete accession/length set and
    succeeds; this check does not claim a residue-by-residue comparison.
  - Adoption regenerates the descriptive subject-annotation sidecar, writes a
    dedicated `blast-resource-manifest.json`, and leaves any legacy
    `manifest.json` untouched.
- `genomes inspect-blast-resource RESOURCE_ID [--catalog PATH] [--cache-dir PATH]`
  - Returns whether the resource is prepared plus
    `gentle.blast_database_inspection.v1`: reference name/release, index kind,
    masking, sequence count, `blastdbcmd` validation, content fingerprint, and
    compatible operations. Subject-annotation status, total/annotated record
    counts, source counts, sidecar path, and fingerprint let users verify that
    gene labels are available before running primer specificity.
  - A prepared `transcriptome_cdna` resource can be used directly as
    `--target-genome` for `primers specificity`,
    `primers specificity-plan`, and whole-panel specificity planning.
- `genomes remove-prepared GENOME_ID [--catalog PATH] [--cache-dir PATH]`
  - Deletes only the prepared install directory for one genome from the selected cache.
  - Catalog JSON is left unchanged.
- `genomes remove-catalog-entry GENOME_ID [--catalog PATH] [--output-catalog PATH]`
  - Deletes only one row from the selected genome catalog JSON.
  - Prepared cache files are left unchanged until `remove-prepared` is run explicitly.
  - If the active catalog file is not writable, `--output-catalog PATH` is
    required so the edited catalog can be written as a copy.
- `genomes blast GENOME_ID QUERY_SEQUENCE [--max-hits N] [--task blastn-short|blastn] [--options-json JSON_OR_@FILE | --options-file PATH] [--catalog PATH] [--cache-dir PATH]`
  - Runs `blastn` against prepared genome cache/index.
  - `--task` defaults to `blastn-short`; accepted values: `blastn-short`, `blastn`.
  - `--options-json` / `--options-file` accept a JSON object that can override quick options and include thresholds (`max_evalue`, `min_identity_percent`, `min_query_coverage_percent`, `min_alignment_length_bp`, `min_bit_score`, `unique_best_hit`).
  - Output includes `binary_preflight` with explicit BLAST tool diagnostics.
- `genomes blast-start GENOME_ID QUERY_SEQUENCE [--max-hits N] [--task blastn-short|blastn] [--options-json JSON_OR_@FILE | --options-file PATH] [--catalog PATH] [--cache-dir PATH]`
  - Starts one async BLAST job and returns a stable `job_id`.
  - Jobs run through a bounded FIFO scheduler and may return initial state
    `queued` or `running` depending on current slot availability.
  - Start schema is `gentle.blast_async_start.v2`: `binary_preflight` is null
    and `binary_preflight_status` is `pending`, never an implied readiness pass.
    Tool probes run in the worker; status returns their diagnostics.
  - Up to 64 jobs can wait beyond the configured worker slots
    (`GENTLE_BLAST_ASYNC_MAX_CONCURRENT`). Excess submissions are rejected.
    Queued jobs preserve submitted options; cancellation before dispatch does
    not launch BLAST. Queue progression does not depend on polling.
- `genomes blast-status JOB_ID [--with-report]`
  - Polls async BLAST job status; `--with-report` includes final report payload when available.
  - Status includes scheduler metadata (`max_concurrent_jobs`, `running_jobs`,
    `queued_jobs`, and `queue_position` while queued).
  - This is a cached, non-mutating observation, not a refresh or wait. Phase,
    update time, request/result hashes and preflight diagnostics are available.
    Retrieve/export the full result while the issuing process is alive: new
    saved project receipts retain status/hashes, not full BLAST reports.
- `genomes blast-cancel JOB_ID`
  - Requests cooperative cancellation for one async BLAST job.
- `genomes blast-list`
  - Lists known async genome-BLAST jobs in the current process.
- `genomes extract-region GENOME_ID CHR START END [--output-id ID] [--annotation-scope none|core|full] [--max-annotation-features N] [--include-genomic-annotation|--no-include-genomic-annotation] [--rmsk-index PATH] [--max-repeat-features N] [--append-repeat-features] [--catalog PATH] [--cache-dir PATH]`
  - Runs engine `ExtractGenomeRegion`.
  - `--annotation-scope` selects projection policy (`none`, `core`, `full`).
  - Default scope is `core` when neither scope nor legacy flags are set.
  - `--max-annotation-features N` applies a deterministic safety cap
    (`0` disables the cap for explicit unrestricted transfer).
  - legacy `--include-genomic-annotation` maps to `--annotation-scope core`.
  - legacy `--no-include-genomic-annotation` maps to `--annotation-scope none`.
  - `--rmsk-index PATH` materializes overlapping UCSC rmsk repeats on the
    extracted anchored sequence; `--max-repeat-features N` caps the generated
    repeat annotations, and `--append-repeat-features` keeps any pre-existing
    generated rmsk annotations.
  - Result payload includes `genome_annotation_projection` telemetry.
- `genomes extract-gene GENOME_ID QUERY [--occurrence N] [--output-id ID] [--extract-mode gene|coding_with_promoter] [--promoter-upstream-bp N] [--annotation-scope none|core|full] [--max-annotation-features N] [--include-genomic-annotation|--no-include-genomic-annotation] [--rmsk-index PATH] [--max-repeat-features N] [--append-repeat-features] [--catalog PATH] [--cache-dir PATH]`
  - Runs engine `ExtractGenomeGene`.
- `genomes extract-promoter GENOME_ID QUERY [--occurrence N] [--transcript-id ID] [--output-id ID] [--upstream-bp N] [--downstream-bp N] [--annotation-scope none|core|full] [--max-annotation-features N] [--include-genomic-annotation|--no-include-genomic-annotation] [--rmsk-index PATH] [--max-repeat-features N] [--append-repeat-features] [--catalog PATH] [--cache-dir PATH]`
  - Runs engine `ExtractGenomePromoterSlice`.
  - Derives one unclipped promoter slice directly from transcript TSS geometry.
  - When `--transcript-id` is omitted, GENtle deterministically uses the
    outermost 5' transcript among the matched gene’s transcript records and
    warns when multiple transcript candidates existed.
  - Defaults are `--upstream-bp 1000` and `--downstream-bp 200`.
- `genomes promoter-tfbs-summary GENOME_ID [--gene QUERY[::OCCURRENCE][@TRANSCRIPT_ID][#DISPLAY_LABEL] ...|--gene-json JSON|--gene-set GROUP_OR_GO|--gene-set-resolution RESOLUTION.json] --motif TOKEN [--motif TOKEN ...|--motifs CSV] [--upstream-bp N] [--downstream-bp N] [--score-kind llr_bits|llr_quantile|llr_background_quantile|llr_background_tail_log10|true_log_odds_bits|true_log_odds_quantile|true_log_odds_background_quantile|true_log_odds_background_tail_log10] [--allow-negative] [--catalog PATH] [--gene-group-catalog PATH] [--cache-dir PATH] [--path FILE.json]`
  - Runs engine `SummarizeMultiGenePromoterTfbs`.
  - Resolves repeated `--gene` tokens into promoter-aligned windows over one
    prepared reference, then scores the requested motif set over every gene.
  - Minus-strand promoters are reverse-complemented before scoring so all gene
    panels share the same transcription-oriented x-axis.
  - Returns portable schema `gentle.multi_gene_promoter_tfbs.v1`.
- `genomes promoter-cohort-comparison GENOME_ID --cohort-label LABEL --cohort-kind manual|co_regulated|anti_co_regulated --gene QUERY[::OCCURRENCE][@TRANSCRIPT_ID][#DISPLAY_LABEL] [--gene ...|--gene-json JSON] --motif TOKEN [--motif TOKEN ...|--motifs CSV] [--source-seq-id SEQ_ID ...] [--upstream-bp N] [--downstream-bp N] [--score-kind llr_bits|llr_quantile|llr_background_quantile|llr_background_tail_log10|true_log_odds_bits|true_log_odds_quantile|true_log_odds_background_quantile|true_log_odds_background_tail_log10] [--allow-negative] [--expression-json JSON] [--expression-source-label LABEL] [--cutrun-dataset ID ...] [--cutrun-read-report ID ...] [--catalog PATH] [--cache-dir PATH] [--path FILE.json]`
  - Runs engine `SummarizePromoterCohortComparison`.
  - First slice supports only manual, co-regulated, and anti-co-regulated
    same-genome cohorts. Use `orthologs resolve-promoter-cohort` and
    `orthologs promoter-comparison` for offline cross-species promoter
    reasoning.
  - Resolves each gene/transcript into a strand-aware promoter window, reuses
    multi-gene TFBS score summaries and TFBS similarity ranking, and emits
    pairwise similarity, shared/common peaks, cohort-specific/outlier peaks,
    optional expression associations, expectation-based `relationship_flags[]`
    for co-/anti-co-regulated review cues, and unresolved-gene warnings.
    Relationship flags are evidence-triage hints, not biological verdicts.
  - Returns portable schema `gentle.promoter_cohort_comparison.v1`.
- `genomes promoter-tfbs-svg GENOME_ID --gene QUERY[::OCCURRENCE][@TRANSCRIPT_ID][#DISPLAY_LABEL] [--gene ...|--gene-json JSON] --motif TOKEN [--motif TOKEN ...|--motifs CSV] [--upstream-bp N] [--downstream-bp N] [--score-kind llr_bits|llr_quantile|llr_background_quantile|llr_background_tail_log10|true_log_odds_bits|true_log_odds_quantile|true_log_odds_background_quantile|true_log_odds_background_tail_log10] [--allow-negative] [--catalog PATH] [--cache-dir PATH] OUTPUT.svg`
  - Runs engine `RenderMultiGenePromoterTfbsSvg`.
  - Exports one small-multiples SVG with one promoter-aligned score-track panel
    per requested gene.
- `genomes extend-anchor SEQ_ID 5p|3p LENGTH_BP [--output-id ID] [--catalog PATH] [--cache-dir PATH] [--prepared-genome GENOME_ID]`
  - Runs engine `ExtendGenomeAnchor`.
  - Extends an already genome-anchored sequence in-silico on contextual `5'` or `3'`.
  - Reprojects core gene/transcript annotation from the prepared reference onto
    the new sequence, including reverse-anchor location/local-CDS transforms.
    Keep the annotation telemetry/warnings; this does not copy the parent's
    imported tracks or authorize reusing an earlier TSS approval.
  - If exact anchor genome id is not prepared but one compatible assembly-family
    cache exists, extension auto-uses that cache and emits a warning.
  - If multiple compatible prepared caches exist, command fails and lists
    explicit options.
- `genomes verify-anchor SEQ_ID [--catalog PATH] [--cache-dir PATH] [--prepared-genome GENOME_ID]`
  - Runs engine `VerifyGenomeAnchor`.
  - Re-checks one anchored sequence against prepared reference sequence at
    recorded coordinates and records `anchor_verified` in provenance.

<!-- manual-page-nav:start -->

---

[CLI Manual](../cli.md) | [Previous: Read-Only SNV Intent](snv-intent.md) | [Next: Helper Sequences](helpers.md)
<!-- manual-page-nav:end -->
