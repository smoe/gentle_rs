# RNA Read Interpretation

<!-- manual-counterparts:start -->
Related GUI reference: [RNA-Read Mapping Workspace](../gui/rna-read-mapping-workspace.md).
These topic links do not imply identical feature coverage across interfaces.
<!-- manual-counterparts:end -->

RNA-read interpretation capability status (Nanopore cDNA phase-1):

- `gentle_cli`: supported via shared-shell/direct commands:
  - `rna-reads interpret`
  - `rna-reads batch-map`
  - `rna-reads align-report`
  - `rna-reads preflight-isoforms`
  - `rna-reads list-reports`
  - `rna-reads show-report`
  - `rna-reads show-alignment`
  - `rna-reads show-alignments`
  - `rna-reads summarize-gene-support`
  - `rna-reads inspect-gene-support`
  - `rna-reads inspect-alignments`
  - `rna-reads inspect-concatemers`
  - `rna-reads build-transcript-index`
  - `rna-reads allele-hash-screen`
  - `rna-reads materialize-hits`
  - `rna-reads export-report`
  - `rna-reads export-hits-fasta`
  - `rna-reads export-sample-sheet`
  - `rna-reads export-target-quality`
  - `rna-reads export-paths-tsv`
  - `rna-reads export-abundance-tsv`
  - `rna-reads export-dexseq-annotation-gff`
  - `rna-reads export-dexseq-counts-tsv`
  - `rna-reads verify-dexseq`
  - `rna-reads export-score-density-svg`
  - `rna-reads export-alignments-tsv`
  - `rna-reads export-isoform-triage-tsv`
  - `rna-reads export-alignment-dotplot-svg`
  backed by `InterpretRnaReads`, `AlignRnaReadReport`,
  `RunRnaReadBatchMap`, `PreflightRnaReadIsoforms`,
  `ListRnaReadReports`, `ShowRnaReadReport`,
  `SummarizeRnaReadGeneSupport`, `InspectRnaReadGeneSupport`,
  `MaterializeRnaReadHitSequences`, `ExportRnaReadReport`,
  `ExportRnaReadHitsFasta`, `ExportRnaReadSampleSheet`,
  `ExportRnaReadTargetQuality`,
  `ExportRnaReadExonPathsTsv`, `ExportRnaReadExonAbundanceTsv`,
  `ExportRnaReadDexseqAnnotationGff`, `ExportRnaReadDexseqCountsTsv`,
  `VerifyRnaReadDexseqExports`,
  `ExportRnaReadScoreDensitySvg`,
  `ExportRnaReadAlignmentsTsv`, `ExportRnaReadIsoformTriageTsv`, and
  `ExportRnaReadAlignmentDotplotSvg`, plus
  the shared engine alignment-display helper used by the GUI `Show alignment`
  pane.
  Input supports FASTA plus gzipped FASTA (`.fa/.fasta` and `.fa.gz/.fasta.gz`).
  Concatenated gzip members are accepted for gzipped FASTA input as well.
  `rna-reads allele-hash-screen` and the standalone
  `gentle_cli allele-hash-screen` command run a deterministic, offline
  allele-aware k-mer screen from an explicit transcript FASTA, transcript-
  coordinate variant table or locally projected reviewed VCF. Reads may come
  from explicit FASTA/FASTQ files or directly from the target-gene accepted
  cohort of a persisted RNA-read report via `--from-rna-report`. Optional
  `--salmon-unmapped-names` and `--salmon-mappings-sam` inputs select Salmon
  cohorts from explicit read files; they do not replace those sequence files.
  Paired inputs are evaluated as fragments without forming k-mers across the
  mate boundary. The output is sequence evidence
  (`gentle.rna_allele_hash_screen.v3`) for haplotype-biased read support. Each
  variant, transcript, phase block, and the gene roll-up carries an auditable
  `hap1_preferred`, `hap2_preferred`, `balanced`,
  `inconclusive_low_depth`, or `unphased_allele_level_only` verdict. Verdicts
  use `--min-informative-reads` (default 10) and the inclusive
  `--balanced-band-lo`/`--balanced-band-hi` (defaults 0.40/0.60); the reported
  two-sided binomial value is advisory and never gates the threshold call.
  RNA-report retained gene/transcript support is labeled expression context,
  with a low-support caveat below the configured informative-depth floor, while
  the operative weight remains allele-informative depth. The report does not
  infer functional allelic imbalance, significance, or causation. Per-read
  `source_origins[]` and aggregate `source_provenance[]` record the input basis. Provenance tags may
  overlap, for example an explicit-file read selected as Salmon-unassigned.
  Existing v1/v2 reports remain deserializable without retrospective verdicts.
  Unphased genotypes remain allele-level evidence, and separate phase sets are
  reported as separate blocks rather than being combined into a fabricated
  haplotype. A transcript/gene roll-up across multiple phased blocks carries a
  `block_local_phase_labels_aggregated` caveat because those labels do not by
  themselves establish one linked gene-wide haplotype.
  Long runs can opt into live progress with the global `--progress-stderr`
  option before the `rna-reads` command; this keeps the final JSON payload on
  stdout while periodic `progress rna-reads ...` lines are written to stderr.
  A complete headless TP73 pancreatic cancer cDNA benchmark recipe, including
  SRA retrieval and TP53/TP63 preflight controls, is maintained in
  `docs/tp73_pancreas_benchmark_runbook.md`.
  The fixed-parameter multi-accession cohort variant for comparing TP73
  abundance across pancreatic cancer runs is maintained in
  `docs/tp73_pancreas_cohort_batch_runbook.md`.
  Shared SRA-backed setup is now available through `reads acquire ...` and uses
  the external SRA Toolkit (`prefetch`, `vdb-validate`, `fasterq-dump`) behind
  one GENtle-owned manifest/status/log contract. RNA and CUT&RUN consume the
  prepared FASTA/FASTQ paths from that shared layer instead of each workflow
  inventing its own downloader.
  `reads acquire prepare` is explicit setup work, emits live
  `progress read-acquisition ...` lines when the global `--progress-stderr`
  option is enabled, records the cooperative cancel-marker path in its activity
  JSON, and re-checks `--min-free-gb` while external SRA Toolkit phases are
  running. Use `reads acquire cancel RUN_ACCESSION ...` to request cancellation
  of an active run from another shell/session.
  Phase split:
  - `interpret`: seed-filter pass (Nanopore phase-1)
  - `batch-map`: manifest-driven wrapper over the same per-sample
    `interpret` + `align-report` engine paths.
    - manifest v1 is TSV with `sample_id` plus either `input_path` or
      `sra_accession`.
    - FASTA rows produce persisted RNA-read reports plus
      `batch_report.json`, `batch_summary.tsv`, `sample_sheet.tsv`,
      `gene_screen_summary.tsv`, `isoform_support.tsv`,
      `concatemer_partner_summary.tsv`, and per-sample gene-support/concatemer
      JSON artifacts under `--out-dir`.
    - `gene_screen_summary.tsv` uses schema
      `gentle.rna_read_gene_screen_summary.v1` and is the preferred
      figure/statistics substrate for cross-gene family plots because it
      records the conservative seed-passed count, all-read q90/q95/q99 length
      quantiles, seed-passed read-length summaries, accepted-target fallback
      fields, and source provenance in one consistent table.
    - by default, SRA-only rows are planned, not fetched: they are marked
      `needs_preparation` and written to `sra_preparation_plan.tsv` plus
      `sra_preparation_commands.sh`.
    - with `--prepare-sra --read-cache-dir DIR --read-work-dir DIR`,
      `batch-map` runs the shared read-acquisition prepare first and maps from
      the prepared FASTA path; `--drop-intermediate-fastq` may remove only the
      conversion FASTQ when the final analysis input is FASTA.
    - default batch settings are `report_mode=full`, `align-selection=all`,
      `complete-rule=near`, `max-secondary-mappings=5`, and
      `continue-on-error=true`.
  - `align-report`: retained-hit alignment pass (phase-2) that updates mapping
    fields, MSA-eligibility counters, exon-transition rows, and exon/junction
    abundance frequencies in the persisted report; retained hits are re-ranked
    by alignment-aware retention rank after alignment.
  - `preflight-isoforms`: non-mutating seed-filter preflight for target
    transcript representation against repeatable positive and negative
    control transcript FASTAs;
    with `--optimize-parameters`, candidate filters must preserve target
    support, pass every `--positive-transcript-fasta` /
    `--must-pass-transcript-fasta` record, and keep each negative control
    group at or below
    `--max-control-match-probability`. The JSON payload includes
    `threshold_recommendation.seed_filter_cli_fragment` and
    `threshold_recommendation.interpret_command_fragment` so control-derived
    thresholds can be pasted into `rna-reads interpret` runs.
    - supports explicit row filtering via
      `--record-indices i,j,k` (0-based stored `record_index` values);
      when provided, this overrides `--selection`.
    - default adapter behavior now uses `selection=all` so rescued retained
      rows also receive round-2 similarity/coverage scores.
  - `summarize-gene-support`: non-mutating target-gene cohort summary over one
    saved aligned RNA-read report.
    - `--gene GENE` is required and repeatable; matches are case-insensitive
      against the same transcript group-label logic used by the splicing view.
    - optional `--record-indices i,j,k` restricts the base aligned cohort to an
      exact saved-report subset before target-gene filtering.
    - `--complete-rule near|strict|exact` controls the fragment-vs-complete
      cohort split; strict and exact counts are still reported explicitly in the
      payload.
    - output is JSON on stdout by default; `--output PATH` writes the exact same
      payload to disk.
    - payload tables include:
      - per-cohort exon support
      - ordered exon-pair support (including skipped pairs like `1->3`)
      - neighboring direct-transition support only (for example `1->2`)
  - `inspect-gene-support`: non-mutating row-level audit of the same
    target-gene cohort logic used by `summarize-gene-support`.
    - `--gene GENE` is required and repeatable, with the same
      case-insensitive group-label matching.
    - optional `--record-indices i,j,k` restricts the evaluation universe
      before cohort classification.
    - `--complete-rule near|strict|exact` controls whether accepted rows land
      in `accepted_fragment` vs `accepted_complete`.
    - `--cohort all|accepted|fragment|complete|rejected` filters the returned
      row list without changing the grouped top-level record-index arrays.
    - output is JSON on stdout by default; `--output PATH` writes the exact
      same payload to disk.
    - rows include resolved gene/transcript identity when available, status +
      machine-readable `status_reason`, full-length flags/class, mapped exon
      ordinals, ordered exon pairs, neighboring direct transitions, phase-2
      score/identity/coverage, and `passed_seed_filter`.
  - `inspect-alignments`: non-mutating ranked alignment inspection
  - `show-alignment`: non-mutating exact phase-2 pairwise alignment detail for
    one saved `record_index`, returning the same machine-readable
    `RnaReadAlignmentDisplay` contract used by the GUI detail pane.
  - `show-alignments`: non-mutating batch wrapper over the same display builder
    for explicit record indices or all rows in one gene-support cohort; records
    without `best_mapping` are reported as skipped instead of aborting the
    whole batch.
  - `inspect-concatemers`: non-mutating fragment/concatemer suspicion audit
    over persisted report hits, with optional structured subset controls:
    `--effect-filter`, `--sort`, `--search`, and `--record-indices`.
  - `build-transcript-index`: prepares a reusable external transcript catalog
    JSON for later concatemer audits, with precomputed per-template k-mer
    positions and preserved source-path provenance.
  - `materialize-hits`: materializes saved-report read sequences back into
    project state for downstream shell/CLI dotplots or other sequence-centric
    analysis without reopening the original FASTA input.
  - `export-target-quality`: exports the target-fragment quality summary for
    one saved report and one or more requested genes/groups.
    - `.json` output creates or extends a comparison bundle so repeated runs
      can accumulate side-by-side across target genes, GENtle versions, and
      parameter changes.
    - `.svg` output writes a comparison figure plus a `.bundle.json` sidecar
      that stores the comparable entries.
    - if the requested SVG path already exists without a reusable GENtle
      sidecar bundle, the export preserves that file and writes a sibling
      `*_compare.svg` instead of overwriting it.
  - `export-alignments-tsv`: non-mutating ranked alignment-row TSV export for
    downstream table-based triage.
  - `export-isoform-triage-tsv`: non-mutating conservative read-level
    classification into known-isoform, ambiguous, gene-supported, and
    off-target/bad-seed bins; it does not call novel isoforms.
  - `export-alignment-dotplot-svg`: non-mutating dotplot-like scatter export
    (coverage vs identity, score-colored points) for aligned report hits.
- `gentle_js`: baseline support via `apply_operation` for the same operation
  family.
- `gentle_lua`: baseline support via `apply_operation` for the same operation
  family.
- GUI workflow parity:
  - Splicing Expert `Nanopore cDNA interpretation` now exposes
    `Prepare Workflow Op` and `Copy Workflow JSON` so the exact
    `InterpretRnaReads` payload can be executed via regular workflow routes.
  - Example:
    - `gentle_cli workflow @docs/examples/workflows/load_branch_reverse_complement_pgex_fasta.json`
    - or paste copied JSON directly:
      `gentle_cli workflow '{"run_id":"workflow_rna_reads_tp73_ncbi_cdna_srr32957124","ops":[{"InterpretRnaReads":{"seq_id":"tp73.ncbi","seed_feature_id":0,"profile":"nanopore_cdna_v1","input_path":"reads.fa.gz","input_format":"fasta","scope":"all_overlapping_any_strand","seed_filter":{"kmer_len":10,"seed_stride_bp":1,"min_seed_hit_fraction":0.30,"min_weighted_seed_hit_fraction":0.05,"min_unique_matched_kmers":12,"max_median_transcript_gap":4.0,"min_chain_consistency_fraction":0.40,"min_confirmed_exon_transitions":1,"min_transition_support_fraction":0.05,"cdna_poly_t_flip_enabled":true,"poly_t_prefix_min_bp":18},"align_config":{"band_width_bp":24,"min_identity_fraction":0.60,"max_secondary_mappings":0},"report_id":"cdna_reads"}}]}'`


- Shared shell:
  - Commands:

## Shared RNA-Read Commands

- Shared shell:
  - Commands:
    - `rna-reads interpret SEQ_ID FEATURE_ID INPUT.fa[.gz] [--report-id ID] [--report-mode full|seed_passed_only] [--checkpoint-path PATH] [--checkpoint-every-reads N] [--resume-from-checkpoint|--no-resume-from-checkpoint] [--profile nanopore_cdna_v1] [--format fasta] [--scope all_overlapping_any_strand|target_group_any_strand|all_overlapping_target_strand|target_group_target_strand] [--origin-mode single_gene|multi_gene_sparse] [--target-gene GENE_ID]... [--roi-seed-capture|--no-roi-seed-capture] [--kmer-len N] [--seed-stride-bp N] [--min-seed-hit-fraction F] [--min-weighted-seed-hit-fraction F] [--min-unique-matched-kmers N] [--min-chain-consistency-fraction F] [--max-median-transcript-gap F] [--min-confirmed-transitions N] [--min-transition-support-fraction F] [--cdna-poly-t-flip|--no-cdna-poly-t-flip] [--poly-t-prefix-min-bp N] [--align-band-bp N] [--align-min-identity F] [--max-secondary-mappings N]`
    - `rna-reads batch-map MANIFEST.tsv --seq-id SEQ_ID --seed-feature-id FEATURE_ID --gene GENE_ID [--gene GENE_ID ...] --out-dir OUT [--target-gene GENE_ID]... [--origin-mode single_gene|multi_gene_sparse] [--report-mode full|seed_passed_only] [--align-selection all|seed_passed|aligned] [--complete-rule near|strict|exact] [--max-secondary-mappings N] [--continue-on-error|--fail-fast] [--prepare-sra] [--read-cache-dir DIR] [--read-work-dir DIR] [--drop-intermediate-fastq] [--transcript-fasta PATH]... [--transcript-index PATH]...`
    - `rna-reads align-report REPORT_ID [--selection all|seed_passed|aligned] [--record-indices i,j,k] [--align-band-bp N] [--align-min-identity F] [--max-secondary-mappings N]`
    - `rna-reads preflight-isoforms SEQ_ID FEATURE_ID [--scope all_overlapping_any_strand|target_group_any_strand|all_overlapping_target_strand|target_group_target_strand] [--positive-transcript-fasta PATH ...] [--must-pass-transcript-fasta PATH ...] [--control-transcript-fasta PATH ...] [--optimize-parameters] [--max-control-match-probability F] [seed filter options]`
    - `rna-reads list-reports [SEQ_ID]`
    - `rna-reads show-report REPORT_ID`
    - `rna-reads show-alignment REPORT_ID RECORD_INDEX`
    - `rna-reads show-alignments REPORT_ID [--gene GENE_ID ... --cohort all|accepted|fragment|complete|rejected --complete-rule near|strict|exact | --record-indices i,j,k] [--limit N] [--output PATH]`
    - `rna-reads summarize-gene-support REPORT_ID --gene GENE_ID [--gene GENE_ID ...] [--record-indices i,j,k] [--complete-rule near|strict|exact] [--output PATH]`
    - `rna-reads inspect-gene-support REPORT_ID --gene GENE_ID [--gene GENE_ID ...] [--record-indices i,j,k] [--complete-rule near|strict|exact] [--cohort all|accepted|fragment|complete|rejected] [--output PATH]`
    - `rna-reads inspect-alignments REPORT_ID [--selection all|seed_passed|aligned] [--limit N] [--effect-filter all_aligned|confirmed_only|disagreement_only|reassigned_only|no_phase1_only|selected_only] [--sort rank|identity|coverage|score] [--search TEXT] [--record-indices i,j,k] [--score-bin-variant all_scored|composite_seed_gate] [--score-bin-index N] [--score-bin-count M]`
    - `rna-reads inspect-concatemers REPORT_ID [--selection all|seed_passed|aligned] [--limit N] [--record-indices i,j,k] [--internal-homopolymer-min-bp N] [--end-margin-bp N] [--max-primary-query-cov F] [--min-secondary-identity F] [--max-secondary-query-overlap F] [--adapter-fasta PATH] [--adapter-min-match-bp N] [--fragment-min-bp N] [--fragment-max-parts N] [--fragment-min-identity F] [--fragment-min-query-cov F] [--transcript-fasta PATH]... [--transcript-index PATH]...`
    - `rna-reads build-transcript-index OUTPUT.json [--kmer-len N] --transcript-fasta PATH [--transcript-fasta PATH ...]`
    - `rna-reads allele-hash-screen --gene GENE --transcript-fasta PATH (--variant-table PATH | --vcf PATH --transcript-map PATH [--vcf-sample SAMPLE]) [--from-rna-report REPORT_ID] [--read-file PATH ...] [--read-pair R1,R2 ...] [--salmon-unmapped-names PATH] [--salmon-mappings-sam PATH] [--read-id-allowlist PATH] [--kmer-len N] [--min-unique-kmer-hits N] [--min-informative-reads N] [--balanced-band-lo F] [--balanced-band-hi F] [--max-inline-read-calls N] --out OUT_DIR`
    - `rna-reads materialize-hits REPORT_ID [--selection all|seed_passed|aligned] [--record-indices i,j,k] [--output-prefix PREFIX]`
    - `rna-reads export-report REPORT_ID OUTPUT.json`
    - `rna-reads export-hits-fasta REPORT_ID OUTPUT.fa [--selection all|seed_passed|aligned] [--record-indices i,j,k] [--subset-spec TEXT]`
    - `rna-reads export-sample-sheet OUTPUT.tsv [--seq-id ID] [--report-id ID]... [--gene GENE_ID]... [--complete-rule near|strict|exact] [--append]`
    - `rna-reads export-target-quality REPORT_ID OUTPUT.{svg|json} --gene GENE_ID [--gene GENE_ID ...] [--complete-rule near|strict|exact]`
    - `rna-reads export-paths-tsv REPORT_ID OUTPUT.tsv [--selection all|seed_passed|aligned] [--record-indices i,j,k] [--subset-spec TEXT]`
    - `rna-reads export-abundance-tsv REPORT_ID OUTPUT.tsv [--selection all|seed_passed|aligned] [--record-indices i,j,k] [--subset-spec TEXT]`
    - `rna-reads export-dexseq-annotation-gff REPORT_ID OUTPUT.gff`
    - `rna-reads export-dexseq-counts-tsv REPORT_ID OUTPUT.tsv [--selection all|seed_passed|aligned] [--record-indices i,j,k] [--subset-spec TEXT]`
    - `rna-reads verify-dexseq REPORT_ID OUTPUT.gff OUTPUT.tsv [--selection all|seed_passed|aligned] [--record-indices i,j,k] [--subset-spec TEXT] [--r-library-path PATH ...]`
    - `rna-reads export-score-density-svg REPORT_ID OUTPUT.svg [--scale linear|log] [--variant all_scored|composite_seed_gate]`
    - `rna-reads export-alignments-tsv REPORT_ID OUTPUT.tsv [--selection all|seed_passed|aligned] [--limit N] [--record-indices i,j,k] [--subset-spec TEXT]`
    - `rna-reads export-isoform-triage-tsv REPORT_ID OUTPUT.tsv [--selection all|seed_passed|aligned] [--limit N] [--record-indices i,j,k] [--subset-spec TEXT] [--min-identity F] [--min-query-coverage F] [--min-confirmed-transition-fraction F] [--max-secondary-mappings N]`
    - `rna-reads export-alignment-dotplot-svg REPORT_ID OUTPUT.svg [--selection all|seed_passed|aligned] [--max-points N]`
    - `rna-reads show-alignments` returns
      `gentle.rna_read_alignment_display_batch.v1` with final
      `selected_record_indices`, per-read `entries[]` containing the same
      display records as `show-alignment`, and `skipped_records[]` for reads
      without `best_mapping`.
    - `rna-reads align-report` re-ranks retained hits by alignment-aware
      retention rank after mapping refresh.
    - `rna-reads inspect-alignments` returns ranked aligned rows suitable for
      read-level inspection without mutating report payloads, and now echoes
      the formal structured subset spec (`effect_filter`, `sort_key`,
      `search`, `selected_record_indices`, `score_density_variant`,
      `score_bin_index`, `score_bin_count`) plus `subset_match_count` in the
      JSON payload so agent-driven inspection stays reproducible. Rows now also
      include full-length flags:
      - `full_length_exact`
      - `full_length_near`
      - `full_length_strict`
    - `rna-reads inspect-concatemers` returns a conservative
      `gentle.rna_read_concatemer_inspection.v1` audit over one saved report.
      It does not claim to prove library chimeras by itself; instead it ranks
      retained reads by fragment-fusion-like signals so suspicious reads can be
      reviewed before stronger conclusions are drawn.
      - current signals are:
        - low primary query coverage
        - internal poly(A) bridge away from read ends
        - internal poly(T) bridge away from read ends
        - disjoint secondary mappings with limited query overlap
        - phase-1 local-block / partial-origin classification
        - optional internal adapter-like matches from an external oligo FASTA
        - optional iterative fragment decomposition over the admitted
          transcript-template set, or one or more appended external transcript
          FASTA catalogs via repeated `--transcript-fasta`, reporting how many distinct
          genes/groups a suspicious read can be explained by
      - committed default adapter signatures live at
        `data/resources/nanopore_direct_cdna_kit14_adapters.fasta`
      - a bundled multiplex/barcoding signature set is also available at
        `data/resources/nanopore_native_barcodes_nbd114_24_v14.fasta`
      - if your study used a different ONT kit, point `--adapter-fasta` at
        the correct oligo FASTA instead of editing engine code
      - when one or more `--transcript-fasta` arguments are provided, the
        top-level payload also
        summarizes repeated non-primary partner genes/transcripts across the
        suspicious cohort so recurring associations can be reviewed directly
      - optional `--record-indices i,j,k` narrows the audit to an exact saved
        aligned-read subset and the returned payload now echoes both
        `selected_record_indices` and `subset_match_count` so transcriptome
        reruns stay reproducible
      - fragment decomposition is BLAST-like but engine-local: the best
        transcript-template fragment is found, masked, and the remaining
        unmasked sequence is queried again until no sufficiently good fragment
        remains or `--fragment-max-parts` is reached
      - `strong` rows currently require a disjoint secondary mapping plus at
        least one additional signal
      - if the aligned report was built with `max_secondary_mappings=0`, the
        payload warns that disjoint-secondary evidence is unavailable and the
        audit becomes weaker
      - practical recommendation:
        rerun `rna-reads align-report REPORT_ID --max-secondary-mappings 3`
        (or higher) before using this audit to argue for concatemer/fusion
        artifacts in real data
      - practical transcript-catalog recipes:
        - coding transcripts only:
          `rna-reads inspect-concatemers REPORT_ID --transcript-fasta /Users/u005069/GitHub/gentle_rs/data/transcriptomes/ensembl/release-116/Homo_sapiens.GRCh38.cdna.all.fa.gz`
        - coding + noncoding transcripts in one run:
          `rna-reads inspect-concatemers REPORT_ID --transcript-fasta /Users/u005069/GitHub/gentle_rs/data/transcriptomes/ensembl/release-116/Homo_sapiens.GRCh38.cdna.all.fa.gz --transcript-fasta /Users/u005069/GitHub/gentle_rs/data/transcriptomes/ensembl/release-116/Homo_sapiens.GRCh38.ncrna.fa.gz`
        - TP73 DeltaN ENA alpha/beta/gamma supplement for focused isoform
          review:
          `rna-reads inspect-concatemers REPORT_ID --selection aligned --record-indices i,j,k --transcript-fasta data/resources/tp73_dn_ena_transcripts.fasta`
          The companion curation/evaluation resource is
          `assets/panels/tp73_dn_isoforms_v1.json`.
        - TP73 D_{Ex2,3}Np73 alpha/beta RefSeq-derived supplement for
          focused exon-skip review:
          `rna-reads inspect-concatemers REPORT_ID --selection aligned --record-indices i,j,k --transcript-fasta data/resources/tp73_delta_ex2_3_refseq_derived_transcripts.fasta`
          The companion curation/evaluation resource is
          `assets/panels/tp73_delta_ex2_3_isoforms_v1.json`.
      - practical TP73-style two-pass cohort review:
          1. run a cheap local suspicion pass first:
             `rna-reads inspect-concatemers REPORT_ID --selection aligned --limit 200 --fragment-max-parts 0`
          2. collect the suspicious `record_index` values from that first
             payload
          3. rerun only those exact rows with the external transcript
             catalogs:
             `rna-reads inspect-concatemers REPORT_ID --selection aligned --record-indices i,j,k --transcript-fasta /Users/u005069/GitHub/gentle_rs/data/transcriptomes/ensembl/release-116/Homo_sapiens.GRCh38.cdna.all.fa.gz --transcript-fasta /Users/u005069/GitHub/gentle_rs/data/transcriptomes/ensembl/release-116/Homo_sapiens.GRCh38.ncrna.fa.gz`
        - reusable prepared transcript catalog for repeated runs:
          `rna-reads build-transcript-index /Users/u005069/GitHub/gentle_rs/data/transcriptomes/ensembl/release-116/Homo_sapiens.GRCh38.cdna_plus_ncrna.rna_read_index.json --transcript-fasta /Users/u005069/GitHub/gentle_rs/data/transcriptomes/ensembl/release-116/Homo_sapiens.GRCh38.cdna.all.fa.gz --transcript-fasta /Users/u005069/GitHub/gentle_rs/data/transcriptomes/ensembl/release-116/Homo_sapiens.GRCh38.ncrna.fa.gz`
        - reuse the prepared index in later audits:
          `rna-reads inspect-concatemers REPORT_ID --transcript-index /Users/u005069/GitHub/gentle_rs/data/transcriptomes/ensembl/release-116/Homo_sapiens.GRCh38.cdna_plus_ncrna.rna_read_index.json`
      - current practical note:
        for first real human transcriptome-scale runs, prefer passing the
        Ensembl FASTAs directly with repeated `--transcript-fasta`. The new
        reusable JSON transcript-catalog index is correct on focused test
        cases, but the current whole-transcriptome JSON format is still heavy
        enough that building a full human shared artifact may be slower and
        larger than is pleasant in day-to-day use.
      - practical memory note:
        when transcriptome-wide reruns are too heavy for a whole saved cohort,
        first use the local audit to collect suspicious `record_index` values,
        then rerun only that narrowed subset with `--record-indices`; using
        `--fragment-max-parts 0` on the first pass disables fragment
        decomposition completely and keeps that initial triage much cheaper
      - implementation note:
        `build-transcript-index` now prepares a reusable transcript catalog
        JSON with per-template k-mer positions so repeated concatemer audits do
        not need to reparse every FASTA on every run. This is the first shared
        prepared-index step for the external transcriptome path. It is not yet
        the heavier fully joined global inverted seed index one might build for
        even larger transcriptome-scale scans.
    - `rna-reads summarize-gene-support` returns a machine-readable
      `gentle.rna_read_gene_support_summary.v1` payload for one or more target
      genes from a saved aligned report:
      - base cohort = aligned retained rows, optionally narrowed by
        `--record-indices`
      - accepted target cohort = base rows whose best mapping resolves to one
        of the requested genes
      - cohort tables are emitted for `all_target`, `fragments`, and
        `complete`
      - `exon_pair_support` includes skipped ordered pairs like `1->3`
      - `direct_transition_support` includes neighboring exon steps only
      - `--output PATH` writes the same JSON returned on stdout
    - `rna-reads export-target-quality` materializes a comparison-friendly
      target-quality artifact from that same shared-engine summary:
      - `.json` output writes a cumulative comparison bundle that appends or
        replaces one matching entry rather than blindly overwriting prior runs
      - `.svg` output writes a comparison figure plus a `.bundle.json`
        sidecar so later exports for another gene, GENtle version, or
        parameter set can be compared directly
      - if a requested SVG path already exists without a reusable GENtle
        sidecar bundle, the engine preserves it and writes `*_compare.svg`
        instead
    - `rna-reads inspect-gene-support` returns a row-level
      `gentle.rna_read_gene_support_audit.v1` payload for the same saved-report
      target-gene evaluation:
      - `status` is one of `unaligned`, `aligned_other_gene`,
        `accepted_fragment`, or `accepted_complete`
      - grouped top-level arrays echo the exact accepted/fragment/complete
        record indices regardless of the chosen row `--cohort` filter
      - rows carry machine-readable `status_reason`, full-length flags/class,
        mapped exon ordinals, exon pairs, direct transitions, and phase-2
        alignment score/identity/coverage fields
    - `rna-reads export-alignments-tsv` writes the same ranked alignment rows
      in TSV form for downstream filtering/sorting; `--record-indices`
      exports an exact saved-report subset and overrides coarse `--selection`;
      `--subset-spec` records the formal subset definition that produced it.
      TSV rows now include:
      - `target_coverage_fraction`
      - `full_length_exact`
      - `full_length_near`
      - `full_length_strict`
      - `full_length_class` (`exact` / `strict_end` / `near` / `partial`)
    - `rna-reads export-paths-tsv` and `rna-reads export-abundance-tsv` now
      accept the same `--record-indices` exact-subset override plus optional
      `--subset-spec` provenance. Abundance exports count per-bin support from
      `exon_path` and exclude `_` intra-exon adjacencies from junction rows.
    - `rna-reads export-dexseq-annotation-gff` writes the report's persisted,
      selection-independent exonic-part partition in the flattened GFF shape
      consumed by DEXSeq. `rna-reads export-dexseq-counts-tsv` writes the
      matching plain two-column HTSeq-style count file; it deliberately emits
      no metadata comments or extra columns. The join key is
      `<aggregate_gene_id>:E<three-digit-part-number>`.
      - reuse one annotation GFF across samples made from the same sequence
        annotation and exonic-part partition
      - each selected retained report hit contributes at most one count per
        touched global ordinal, even if that ordinal repeats in `exon_path`
      - `_empty` counts selected rows with an empty `exon_path`, `_lowaqual`
        counts selected rows that failed the seed gate, and `_notaligned`
        counts selected rows without a best alignment; these independently
        observed diagnostics can overlap
      - `_ambiguous` and `_ambiguous_readpair` are currently always `0`
        because GENtle does not persist those HTSeq assignment categories
      - command selection and `subset_spec` provenance are returned as JSON,
        not inserted into the strict two-column count file
      - reports that predate persisted `exonic_part_bins` are rejected; re-run
        alignment rather than reconstructing a possibly different partition at
        export time
      - pass the per-sample count files and shared GFF to R as
        `DEXSeqDataSetFromHTSeq(countfiles, sampleData, design, flattenedfile)`
      - `rna-reads verify-dexseq` writes both artifacts and constructs a real
        `DEXSeqDataSetFromHTSeq()` with two pseudo-samples to verify that the
        GFF and count keys join. It never installs R packages; its JSON result
        reports `Rscript`/`DEXSeq` readiness, effective R library paths, the
        generated command, and a `DEXSEQ_OK` summary when verification runs.
    - `rna-reads export-isoform-triage-tsv` writes a conservative read-level
      TSV with bins `known_isoform_confirmed`, `known_isoform_ambiguous`,
      `gene_supported_no_isoform_call`, and `off_target_or_bad_seed`; recurrent
      unresolved patterns are left for a later aggregate novelty workflow.
    - `rna-reads export-alignment-dotplot-svg` emits a dotplot-style alignment
      scatter (coverage vs identity) with score-based point coloring.
    - `rna-reads export-hits-fasta` headers include seed metrics and exon-path
      annotations:
      - `--record-indices i,j,k` exports the exact saved-report subset and
        overrides coarse `--selection`
      - `--subset-spec TEXT` records the formal subset definition in the FASTA
        headers
      - `exon_path_tx=<transcript_id|none>`
      - `transcript_exon_path=<ordinal_path|none>` uses transcript-local exon
        numbering and is the primary human-facing path; `:` marks
        hash-confirmed adjacent exon transitions and `-` marks unconfirmed
        adjacency
      - `exon_path=<ordinal_path|none>` is the DEXSeq-style genomic projection
        into disjoint exonic-part bins, numbered by ascending genomic
        coordinate; `_` marks adjacent bins inside the same transcript-local
        exon (no splice), while `:` / `-` retain confirmed / unconfirmed splice
        junction semantics
      - schema note: these global projection ordinals are renumbered relative
        to pre-exonic-parts releases; `transcript_exon_path` is unchanged and
        remains the primary human-facing field
      - `exon_transitions=<confirmed>/<total>`
      - `rc_applied=<true|false>` indicating automatic cDNA poly-T
        reverse-complement normalization was applied
      - `origin_class=<...>` with `origin_conf` / `strand_conf` summary values
        from deterministic origin classification

## Normalization and Seed Span

- Shared shell:
  - Commands:
    - cDNA/direct-RNA seed-normalization semantics:
      - default (`--cdna-poly-t-flip`): reads with a T-rich 5' head are
        reverse-complement normalized before hashing/scoring
        - head detection tolerates small interruptions/mismatches in the
          poly-T stretch near the 5' end
      - direct RNA (`--no-cdna-poly-t-flip`): reads are scored as provided
      - `--poly-t-prefix-min-bp` controls the minimum T support threshold used
        by the 5' head detector for automatic cDNA flip (default `18`)
    - phase-1 seed-span policy:
      - all reads are hashed across their full span (no 3-window sampling)
      - `--seed-stride-bp` is the real phase-1 density knob: `1` means one
        ordered k-mer start per base, while larger values trade sensitivity for
        speed
    - sparse-origin options:
      - `--origin-mode` accepts `single_gene` (default) and
        `multi_gene_sparse`
      - `--target-gene` (repeatable) and `--roi-seed-capture` are persisted in
        the report payload for deterministic follow-up runs
      - `multi_gene_sparse` expands template indexing with local-annotation
        matches from `--target-gene`
      - `--roi-seed-capture` remains planned and is surfaced as a deterministic
        warning until implemented
    - `rna-reads list-reports` rows expose sparse-origin request provenance:
      `origin_mode`, `input_orientation_mode`, `input_orientation_label`,
      `target_gene_count`, and `roi_seed_capture_enabled`
      - shell/CLI JSON output also includes `summary_rows[]` for quick
        human-readable triage (`mode`, cDNA/direct-RNA orientation, `origin`,
        target count, ROI-capture flag, and read counters)
    - `rna-reads show-report` output includes a `summary` field with the same
      compact provenance framing and top-level `input_orientation_mode` /
      `input_orientation_label` JSON markers, plus aligned/full-length
      percentages and compact read-length histogram summaries for:
      - all reads
      - seed-passed reads
      - aligned reads
      - full-length exact / near / strict subsets
    - tutorial reference (TP53 basis + multi-gene sparse mapping):
      `docs/tutorial/generated/chapters/07-02_tp53_multi_gene_sparse_mapping_online.md`
    - batch target-gene cohort tutorial:
      `docs/tutorial/07-01_rna_read_batch_gene_support_cli.md`
    - `rna-reads export-sample-sheet` includes sparse-origin provenance columns
      (`report_mode`, `input_orientation_mode`, `input_orientation_label`,
      `origin_mode`, `target_gene_count`, `target_gene_ids_json`,
      `roi_seed_capture_enabled`), `mean_read_length_bp`, and
      `origin_class_counts_json` alongside exon/junction frequency JSON fields
    - when `rna-reads export-sample-sheet` is given one or more `--gene`
      arguments, each report row also carries target-gene cohort metrics:
      - accepted-target counts/fractions
      - fragment vs complete counts (`near|strict|exact`)
      - `gene_support_mean_assigned_read_length_bp`
      - JSON-serialized exon support, ordered exon-pair co-presence, and
        neighboring direct-transition support

## Report Compaction and Checkpoints

- Shared shell:
  - Commands:
    - report compaction and checkpoint options:
      - `--report-mode full` (default): persist retained top hits as ranked
      - `--report-mode seed_passed_only`: persist only retained hits that
        passed the seed gate (stream counters are unchanged)
      - `--checkpoint-path PATH --checkpoint-every-reads N`: persist
        deterministic resume snapshots during streaming
      - `--resume-from-checkpoint`: load checkpoint snapshot and continue from
        the saved processed-read index
    - phase-1 seed-pass gate policy:
      - `pass = raw_hit_fraction >= min_seed_hit_fraction AND weighted_hit_fraction >= min_weighted_seed_hit_fraction AND unique_matched_kmers >= min(min_unique_matched_kmers, tested_kmers) AND chain_consistency_fraction >= min_chain_consistency_fraction AND median_transcript_gap <= max_median_transcript_gap AND confirmed_transitions >= min_confirmed_transitions AND confirmed_transition_fraction >= min_transition_support_fraction`
      - `chain_consistency_fraction` is the fraction of matched seed observations
        explained by one coherent transcript-offset chain
      - `weighted_hit_fraction` is occurrence-normalized inside the scoped seed
        index: `sum(1/occurrence_count(seed_bits))/tested_kmers`
      - defaults:
        - `--kmer-len 10`
        - `--min-seed-hit-fraction 0.30`
        - `--min-weighted-seed-hit-fraction 0.05`
        - `--min-unique-matched-kmers 12`
        - `--min-chain-consistency-fraction 0.40`
        - `--max-median-transcript-gap 4.0`
        - `--min-confirmed-transitions 1`
        - `--min-transition-support-fraction 0.05`

## Evidence Scope

- Shared shell:
  - Commands:
    - RNA-read scope semantics:
      - `all_overlapping_any_strand`: index all overlapping transcripts on any
        strand, including antisense/opposite-strand genes relative to the
        selected target gene/group
      - `target_group_any_strand`: index only target-group transcripts, any
        annotated strand allowed
      - `all_overlapping_target_strand`: index all overlapping transcripts on
        the selected target gene/group's annotated strand only
      - `target_group_target_strand`: index only target-group transcripts on the
        selected target gene/group's annotated strand
      - strand-scoring note:
        any-strand modes score against the union of admitted target-gene-strand
        and antisense/opposite-strand templates; target-gene-strand modes
        exclude antisense/opposite-strand templates from scoring.
      - seed-index note:
        indexed seeds include exon-body and exon-exon transition k-mers for
        admitted transcripts.

<!-- manual-page-nav:start -->

---

[CLI Manual](../cli.md) | [Previous: RNA Read Panel](rna-read-panel.md) | [Next: Primer-Anchored Read Trace](primer-trace.md)
<!-- manual-page-nav:end -->
