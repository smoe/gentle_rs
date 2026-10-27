# Primer Design and Handoff

<!-- manual-counterparts:start -->
Related GUI reference: [Primer, qPCR, and Sequence-Specific RT Design](../gui/primer-qpcr-and-sequence-specific-rt-design.md).
These topic links do not imply identical feature coverage across interfaces.
<!-- manual-counterparts:end -->

Primer-design report capability status:

- `gentle_cli`: supported via shared-shell `primers ...` commands and direct forwarding (`gentle_cli primers ...`), backed by `DesignPrimerPairs`, `DesignQpcrAssays`, `BuildTranscriptQpcrPanel`, `AssessPrimerPairSpecificity`, and the post-design cloning handoff operation `PrepareRestrictionCloningPcrHandoff`, plus non-mutating ROI seed helpers (`primers seed-from-feature`, `primers seed-from-splicing`), a local BLAST specificity confirmation route (`primers specificity`), a non-mutating restriction-cloning handoff request seeder (`primers seed-restriction-cloning-handoff`), vector suggestion helpers (`primers restriction-cloning-vector-suggestions`), and persisted report inspect/export helpers for primer, qPCR, and restriction-cloning handoff reports. `--progress` now also streams `progress primers ...` lines for shared-shell primer/qPCR design commands, not only raw `op` / `workflow` JSON execution.
- Saved primer/qPCR reports and materialized primer/probe sequences are
  in-silico design instantiations: they satisfy a PCR/qPCR design constraint,
  but they do not mean the physical oligo has been ordered, received, or found
  in local inventory.
- `gentle-cloning` ClawBio skill: exposes typed request modes over the same
  shared PCR/qPCR/TaqMan surface, including Primer3/backend preflight, PCR and
  qPCR seed helpers, `DesignPrimerPairs`/`DesignQpcrAssays` payload execution,
  cDNA PCR/qPCR assay tests, primer/qPCR report list/show/export, restriction-
  cloning PCR handoffs, vector suggestions, and PCR-family protocol cartoons.
- `gentle_js`: supported via `apply_operation` (`DesignPrimerPairs`, `DesignQpcrAssays`, `BuildTranscriptQpcrPanel`, `DesignTranscriptAssayPanel`, `PrepareRestrictionCloningPcrHandoff`); the returned project state and operation result carry the persisted transcript-panel report
- `gentle_lua`: supported via `apply_operation` (`DesignPrimerPairs`, `DesignQpcrAssays`, `BuildTranscriptQpcrPanel`, `DesignTranscriptAssayPanel`, `PrepareRestrictionCloningPcrHandoff`); the returned project state and operation result carry the persisted transcript-panel report


- Shared shell:
  - Commands:
    - `primers design REQUEST_JSON_OR_@FILE [--backend auto|internal|primer3] [--primer3-exec PATH]`
    - `primers design-transcript-capture-pool REQUEST_JSON_OR_@FILE`
      - offline exact-site discovery for explicit transcript sets, sense or
        antisense roles, UTR/range windows, sharing permissions and fixed-oligo
        stage-aware pool diagnostics; never an approved assay or order
      - use `primers list-reports`, `show-report`, and `export-report` for the
        saved result; see [request examples and limits](../transcript_capture_pools.md)
    - `primers design-terminal-exon-rt-pool REQUEST_JSON_OR_@FILE`
      - the optional request object `genomic_specificity` accepts
        `target_genome_id`, optional `catalog_path`/`cache_dir`, and
        `hit_warning_threshold` (default `20000`); GENtle then runs exhaustive
        local `blastn-short` checks for every selected variable segment and
        persists exact genomic hits plus database-content provenance in the
        same report
      - accepts `gentle.terminal_exon_rt_primer_pool_request.v1`
      - resolves each ordered target through the annotated Splicing Expert
        transcript model and searches exact-length sequence-specific segments
        near the transcript-oriented start of its terminal exon
      - the complete oligo is the fixed 5-prime adapter followed by the reverse
        complement of the selected mature-transcript segment
      - caller order identifies the requested pool members. GENtle retains a
        bounded candidate list for each member and selects the final set
        jointly, protecting complete-oligo/variable self-complementarity and
        homopolymer quality before minimizing mutual 3-prime and general
        complementarity and distance; Tm is reported but does not affect
        ranking
      - example request file:
        ```json
        {
          "schema": "gentle.terminal_exon_rt_primer_pool_request.v1",
          "fixed_adapter_5prime": "ACTTGCCTGTCGCTCTATCTTC",
          "variable_length_bp": 22,
          "terminal_exon_search_window_bp": 250,
          "max_candidates_per_target": 10,
          "targets": [
            {
              "seq_id": "tp73_locus",
              "source_feature_id": 0,
              "transcript_id": "NM_005427.4",
              "label": "TP73"
            },
            {
              "seq_id": "patz1_locus",
              "source_feature_id": 0,
              "transcript_id": "NM_032052.4",
              "label": "PATZ1"
            }
          ],
          "genomic_specificity": {
            "target_genome_id": "human_grch38_ensembl_116",
            "hit_warning_threshold": 20000
          },
          "report_id": "tp73_patz1_terminal_exon_rt_pool"
        }
        ```
        The sequence, feature, and transcript identifiers are illustrative;
        use the exact identifiers shown for the loaded annotations in GENtle's
        Splicing Expert or introspection output.
      - run with
        `gentle_cli primers design-terminal-exon-rt-pool @terminal_exon_rt_pool.json`
      - without `genomic_specificity`, this report evaluates terminal-exon
        geometry and within-request pool compatibility only. With it, the
        report additionally contains exhaustive, content-fingerprint-bound
        genomic-DNA exact-hit evidence; it still does not claim
        whole-transcriptome specificity
      - self/cross-complementarity is reported as complementary-base runs,
        not thermodynamic hairpin/dimer predictions; see
        [self- and cross-interaction review](../transcript_capture_pools.md#self--and-cross-interaction-review)
        for the existing checks and the separate Primer3 follow-up
    - `primers design-qpcr REQUEST_JSON_OR_@FILE [--backend auto|internal|primer3] [--primer3-exec PATH]`
    - `primers specificity REPORT_ID --pair-rank N --target-genome GENOME_ID [--max-target-amplicon-bp N | --readiness-max-amplicon-bp N --exploratory-max-amplicon-bp N] [--report-detail compact|full] [--min-primer-coverage-fraction F] [--max-3prime-mismatches N] [--three-prime-window-bp N] [--min-total-mismatches-to-unintended-target N] [--max-hits-per-primer N] [--path OUTPUT.json]`
    - `primers specificity --forward SEQ --reverse SEQ --target-genome GENOME_ID [--max-target-amplicon-bp N | --readiness-max-amplicon-bp N --exploratory-max-amplicon-bp N] [--report-detail compact|full] [--min-primer-coverage-fraction F] [--max-3prime-mismatches N] [--three-prime-window-bp N] [--min-total-mismatches-to-unintended-target N] [--max-hits-per-primer N] [--path OUTPUT.json]`

- Shared shell:
  - Commands:
    - `primers specificity-plan REPORT_ID --pair-rank N --target-genome GENOME_ID --output-dir DIR [same policy/catalog/cache options as specificity]`
    - `primers specificity-plan --forward SEQ --reverse SEQ --target-genome GENOME_ID --output-dir DIR [same policy/catalog/cache options as specificity]`
    - `primers specificity-import HANDOFF.json [--path OUTPUT.json]`
    - `primers transcript-assay-specificity-plan PANEL_REPORT_ID --target-genome GENOME_ID --output-dir DIR [same policy/catalog/cache options as specificity]`
    - `primers transcript-assay-specificity-finalize HANDOFF.json EXECUTION_MANIFEST_JSON_OR_@FILE [--path ACCEPTANCE.json]`
    - `primers transcript-assay-specificity-redesign REQUEST_JSON_OR_@FILE [--path OUTPUT.json]`

## PrimerBank and Transcript Checks

- Shared shell:
  - Commands:
    - PrimerBank lookup and transcript checks:
      - `primers primerbank search QUERY [--by gene-symbol|gene-id|genbank|protein|primerbank-id|keyword] [--species human|mouse|all] [--html SAVED.html] [--path OUTPUT.json]`
      - `primers primerbank show PRIMERBANK_ID [--species human|mouse|all] [--html SAVED.html] [--path OUTPUT.json]`
      - `primers primerbank test-cdna SEQ_ID FEATURE_ID PRIMERBANK_ID --species human|mouse [--html SAVED.html] [--transcript-id ID] [--min-amplicon-bp N] [--max-amplicon-bp N] [--max-mismatches N] [--require-3prime-exact-bases N] [--transcript-order transcript_id|genomic_first_exon|genomic_last_exon|antisense_first_exon] [--map-coordinate-mode cdna|genomic_aligned] [--path OUTPUT.json] [--svg OUTPUT.svg]`
      - live lookup submits one query to PrimerBank's public HTML search form;
        `--html` parses a saved result instead for offline/reproducible work.
        GENtle preserves PrimerBank coding-sequence coordinates and source links
        but does not reinterpret catalog presence as current-transcript
        compatibility, genomic specificity, or experimental validation.
        `test-cdna` is the explicit compatibility continuation through the
        existing transcript-aware PCR route; whole-genome specificity remains
        a separate `primers specificity-plan`/`specificity-import` step.
        Search/show JSON includes `species_check` and per-gene match status.
        `test-cdna` requires a concrete expected species and refuses mismatch
        or unresolved PrimerBank records before running the transcript assay.
        It also compares an available `organism` annotation on the selected
        project sequence/transcript with that expected species. A known target
        mismatch is refused; absent target-species metadata is reported as
        `unresolved` with a warning rather than silently presented as verified.
        `--path` preserves the complete catalog-plus-cDNA report, while
        `--svg` writes its transcript map.
    - `primers import-external-pairs INPUT.json|tsv SEQ_ID FEATURE_ID [--format auto|json|tsv] [--report-id ID] [--transcript-id ID] [--transcript-order transcript_id|genomic_first_exon|genomic_last_exon|antisense_first_exon] [--map-coordinate-mode cdna|genomic_aligned] [--min-amplicon-bp N] [--max-amplicon-bp N] [--max-mismatches N] [--require-3prime-exact-bases N] [--specificity-target-genome GENOME_ID] [--specificity-catalog PATH] [--specificity-cache-dir DIR] [--artifact-output-dir DIR] [--materialize-products] [--product-gel-ladder NAME]... [--path OUTPUT.json]`
      - imports exact 5-prime-to-3-prime forward/reverse sequences with provider,
        catalogue/source/accession provenance and aliases. Commercial rows require
        both provider and catalogue id
      - whitespace and copied position digits are removed, RNA `U` is normalized
        to DNA `T`, and the remaining sequence must be valid IUPAC
      - duplicate oriented sequence pairs share stable sequence-derived oligo/pair
        ids while every distinct source row and claim remains attached
      - GENtle recomputes cDNA products, transcript maps, genomic carryover, oligo
        QC, Tm, GC, and pair delta-Tm. Source targeting/validation claims never
        establish transcript coverage or specificity
      - specificity is `not_run` unless a prepared target is requested. For a
        requested run, GENtle derives intended transcript products from its
        computed cDNA assay and projects genomic geometry only from a
        provenance-bearing project genome anchor; provider claims never supply
        missing geometry.
        `--materialize-products` makes the import state-changing and, with
        `--artifact-output-dir`, writes the corresponding product gel
    - `primers screen-variants REQUEST_JSON_OR_@FILE [--path OUTPUT.json] [--evidence-dir DIR]`
      - accepts `gentle.primer_variant_screen_request.v1` with one or more
        primer pairs, explicit current-assembly binding segments, and either a
        local VCF/VCF.gz or local variant-resource manifest
      - scans the VCF once for all physical pairs. Junction-spanning oligos use
        separate exonic segments; reverse-strand positions and the configurable
        critical 3-prime window are reported in oligo coordinates
      - verifies assembly, contig, and overlapping REF alleles. Missing AF is
        unknown rather than zero, and indels/MNVs remain conservative evidence
        rather than being silently normalized into a haplotype claim
      - optional `source.annotation_info_fields` (or the corresponding manifest
        field) retains selected local VCF INFO values for inspection, including
        dbSNP-derived or splice-annotation context, without interpreting those
        values as causal evidence
      - setting `degenerate_rescue_minimum_frequency` emits a deterministic
        mixed-IUPAC proposal only for passing, reference-verified simple primer
        SNVs at or above that known population frequency. Critical 3-prime
        changes additionally require
        `allow_critical_3prime_degenerate_rescue=true`
      - every proposal reports its synthesis-mixture complexity and receives a
        new physical pair id. It is not an order-ready genotype claim: rerun
        specificity, thermodynamic, cDNA, and experimental validation for the
        adjusted pair. cDNA PCR/qPCR already interprets IUPAC symbols as allele
        sets, equivalent to regular-expression character classes
      - `--path` writes the aggregate `gentle.primer_variant_screen.v1` report;
        `--evidence-dir` writes one `gentle.primer_variant_evidence.v1` JSON per
        pair for direct use with `primers experimental-handoff
        ... --variant-evidence PATH`
      - candidate source rows remain provenance only. PrimerBank or commercial
        catalog presence is not GENtle validation, and this route never
        downloads dbSNP or another variant resource
    - `primers test-cdna-pcr SEQ_ID FEATURE_ID --forward SEQ --reverse SEQ [--transcript-id ID] [--transcript-order transcript_id|genomic_first_exon|genomic_last_exon|antisense_first_exon] [--map-coordinate-mode cdna|genomic_aligned] [--min-amplicon-bp N] [--max-amplicon-bp N] [--max-mismatches N] [--require-3prime-exact-bases N] [--path OUTPUT.json] [--svg OUTPUT.svg] [--materialize-products] [--product-output-prefix PREFIX] [--product-gel-svg OUTPUT.svg] [--product-gel-ladder NAME]...`
    - `primers test-cdna-qpcr SEQ_ID FEATURE_ID --forward SEQ --reverse SEQ --probe SEQ [--transcript-id ID] [--transcript-order transcript_id|genomic_first_exon|genomic_last_exon|antisense_first_exon] [--map-coordinate-mode cdna|genomic_aligned] [--min-amplicon-bp N] [--max-amplicon-bp N] [--max-mismatches N] [--require-3prime-exact-bases N] [--path OUTPUT.json] [--svg OUTPUT.svg] [--materialize-products] [--product-output-prefix PREFIX] [--product-gel-svg OUTPUT.svg] [--product-gel-ladder NAME]...`
    - `primers transcript-qpcr-panel SEQ_ID FEATURE_ID SHARED_QPCR_REPORT_ID [--path OUTPUT.json]`
    - `primers design-transcript-assay-panel OPERATION_JSON_OR_@FILE [--backend auto|internal|primer3] [--primer3-exec PATH]`
    - `primers design-group-target REQUEST_JSON_OR_@FILE [--path OUTPUT.json] [--backend auto|internal|primer3] [--primer3-exec PATH]`
    - `primers design-transcript-assay-panel SEQ_ID FEATURE_ID [--assay-kind endpoint-rt-pcr|sybr-qpcr|taqman-qpcr] [--cdna-synthesis oligo-dt|random-hexamers|gene-specific|mixed] [--objective pan-transcript|one-per-class|minimal-discrimination-panel|isoform-end-matrix] [--coverage-policy require-all|best-effort] [--coverage-universe JSON_OR_@FILE] [--assay-tier routine-common-region-screen|isoform-discrimination|long-range-structure-discovery] [--preferred-min-amplicon-bp N --preferred-max-amplicon-bp N] [--junctions JSON_OR_@FILE] [--junction-evidence PATH ...] [--junction-evidence-priority required|preferred] [--min-3prime-junction-overlap-bp N] [--min-5prime-junction-overlap-bp N] [--annotation-release TEXT] [--min-amplicon-bp N] [--max-amplicon-bp N] [--max-assays-per-class N] [--max-mismatches N] [--require-3prime-exact-bases N] [--oligo-dt-5prime-risk-threshold-bp N] [--report-id ID] [--path OUTPUT.json] [--backend auto|internal|primer3] [--primer3-exec PATH]`
      - The group-target route accepts `gentle.primer_group_target_design.v1`
        request fields as JSON. It aligns two or more loaded related sequences,
        derives exact shared primer intervals, refuses excessive search before
        Primer3, and reports every retained pair against every member.
        `require_all` is strict; incomplete coverage is returned only when
        `coverage_policy` is explicitly `best_effort`.
      - The same `DesignPrimerGroupTarget` operation is reachable through
        generic `op`/workflow, MCP, JavaScript, and Lua adapters. The dedicated
        command is the discoverable shell projection, not a second algorithm.
      - for oligo-dT cDNA, matrix cells always report each predicted product's
        required reach from the annotated transcript 3-prime end. Supplying
        `--oligo-dt-5prime-risk-threshold-bp` classifies that geometry against
        the experiment-specific threshold; omitting it leaves the distance
        visible without inventing a universal safe cutoff. A missing primer
        target remains structurally distinct from an indeterminate no-product
        result.
    - `primers compose-gene-assay-routine REQUEST_JSON_OR_@FILE [--path OUTPUT.json]`
      - joins one exported isoform-evidence ledger with persisted transcript
        assay panels into `gentle.gene_transcript_assay_routine.v1`
      - accepts either the bare `gentle.gene_isoform_evidence.v1|v2` report or
        the tagged `isoform_evidence` feature-expert JSON returned by the
        shared inspection route
      - accepts an optional `expected_isoform_evidence_sha256`; mismatch is a
        hard error, while missing or stale panel specificity remains visibly
        unaccepted
      - this command only composes existing reports and digests. It does not
        rerun design or specificity and does not change project state
      - example request:
        ```json
        {
          "label": "PATZ1 common, junction, and endpoint review",
          "isoform_evidence_path": "analysis/patz1_isoform_evidence.json",
          "expected_isoform_evidence_sha256": "sha256:<digest>",
          "transcript_assay_panel_report_ids": [
            "patz1_common_control",
            "patz1_junction_sybr",
            "patz1_endpoint_matrix"
          ]
        }
        ```
        Save this as `analysis/patz1_assay_routine_request.json`, then run
        `primers compose-gene-assay-routine @analysis/patz1_assay_routine_request.json --path analysis/patz1_assay_routine.json`.
    - `primers experimental-handoff PANEL_REPORT_ID [--policy JSON_OR_@FILE] [--variant-evidence PATH ...] [--order-form-id ID] [--path OUTPUT.json] [--order-table OUTPUT.tsv] [--gel-svg OUTPUT.svg] [--gel-ladder NAME ...]`
    - `primers test-cdna-qpcr-fasta CDNA_FASTA[.gz] [CDNA_FASTA[.gz] ...] --forward SEQ --reverse SEQ --probe SEQ [--transcript-id ID] [--min-amplicon-bp N] [--max-amplicon-bp N] [--max-mismatches N] [--require-3prime-exact-bases N] [--path OUTPUT.json] [--svg OUTPUT.svg]`
    - `primers preflight [--backend auto|internal|primer3] [--primer3-exec PATH]`
      - returns `gentle.primer3_preflight.v1` with configured/effective
        executable, resolved path, working directory, and
        reachability/version diagnostics plus tri-state `progress_supported`
        capability detection from `--help`
    - `primers seed-from-feature SEQ_ID FEATURE_ID`
    - `primers seed-from-splicing SEQ_ID FEATURE_ID`
    - `primers seed-qpcr-from-feature SEQ_ID FEATURE_ID`
    - `primers seed-qpcr-from-splicing SEQ_ID FEATURE_ID [--mode shared_gene|distinguish_transcript] [--transcript-id ID] [--specificity-evidence junction_only|unique_exon_or_chain|either_prefer_junction]`
    - `primers list-reports`
      - lists persisted primer-design reports, terminal-exon RT-primer-pool
        reports, and primer-specificity artifacts in separate arrays
    - `primers show-report REPORT_ID`
      - accepts all three kinds; design reports include `simple_pcr_pairs` with
        per-pair left/right distance from the core ROI, overlap flags, and
        flanking labels for quick CLI inspection. The full design report also
        exposes additive score terms, bounded evaluated rejection rows,
        capture completeness/accounting, and its construct-reasoning graph id.
    - `primers export-report REPORT_ID OUTPUT.json`
      - exports a primer-design report, terminal-exon RT-primer-pool report, or
        persisted primer-specificity artifact without rerunning analysis
    - `primers list-qpcr-reports`
    - `primers show-qpcr-report REPORT_ID`
      - persisted qPCR report output now includes `best_assay_summary` plus
        machine-readable `best_assay_probe_placement`
    - `primers export-qpcr-report REPORT_ID OUTPUT.json`
    - `primers list-transcript-assay-panels`
    - `primers show-transcript-assay-panel REPORT_ID`
    - `primers export-transcript-assay-panel REPORT_ID OUTPUT.json`
    - `primers list-transcript-assay-fallbacks`
    - `primers show-transcript-assay-fallback EXECUTION_ID`
    - `primers export-transcript-assay-fallback EXECUTION_ID OUTPUT.json`
    - `primers oligo-order create REQUEST_JSON_OR_@FILE`
      - creates and persists a first-class `gentle.oligo_order_form.v1` from
        supplied line items; duplicate line items are kept and reported
    - `primers oligo-order from-primer-report REPORT_ID --pair-rank N[,N...] [--form-id ID] [--scale TEXT] [--purification TEXT] [--modification TEXT ...]`
      - expands each selected primer pair into forward/reverse order lines with
        report-rank provenance
    - `primers oligo-order from-qpcr-report REPORT_ID --assay-rank N[,N...] [--include-probe true|false] [--form-id ID] [--scale TEXT] [--purification TEXT] [--modification TEXT ...]`
      - expands each selected qPCR assay into forward/reverse/probe order lines
        by default
    - `primers oligo-order list`
    - `primers oligo-order show FORM_ID`
    - `primers oligo-order export FORM_ID OUTPUT.json`
    - `primers oligo-order route FORM_ID`
      - classifies the stored form through the existing
        `services delivery-route` logic using `source_target.kind =
        oligo_order_form`
    - `primers oligo-order quote FORM_ID [--provider metabion|geneart] [--service-kind KIND] [--output-dir DIR]`
      - prepares the existing external-service quote handoff bundle from a
        duplicate-reviewed order form; omitted provider/service kind values use
        the route recommendation, and `--output-dir` writes the normalized
        JSON/CSV/Markdown bundle files
    - `primers oligo-order review-dedup FORM_ID [--reviewer NAME] [--duplicate-action keep-separate] [--note TEXT]`
      - records review of exact procurement duplicates without merging or
        removing lines

- Shared shell:
  - Commands:

## Primer Request Payloads

- Shared shell:
  - Commands:
    - Primer request payload notes (`primers design` / `primers design-qpcr`):
      - `forward`/`reverse`/`probe` side constraints support optional
        sequence filters:
        `non_annealing_5prime_tail`, `fixed_5prime`, `fixed_3prime`,
        `required_motifs[]`, `forbidden_motifs[]`, `locked_positions[]`.
      - `non_annealing_5prime_tail` is added to the oligo sequence while
        anneal `tm_c`/`gc_fraction`/`anneal_hits` remain computed on the
        template-binding segment only.
      - to design a partner around one complete supplied canonical primer,
        set that side's `min_length == max_length == fixed_5prime.len()` and
        put the full 5'-to-3' sequence in `fixed_5prime`. GENtle passes this to
        Primer3 as `SEQUENCE_PRIMER` or `SEQUENCE_PRIMER_REVCOMP`; a shorter
        `fixed_5prime` remains only a prefix filter.
      - `pair_constraints` is optional and supports:
        `require_roi_flanking`, amplicon motif filters, and fixed amplicon
        start/end coordinates.
      - `pair_constraints.rejected_near_miss_limit` controls deterministic
        evaluated rejection retention (`null`/omitted = `20`,
        `0` = disabled, maximum `100`). This does not change candidate
        selection, rank, or score. Primer-pair and insertion-primer reports
        retain evaluated pair-level rows; qPCR reports retain evaluated
        pair/probe assay-level rows.
      - retained pairs from both the internal and Primer3 proposal paths use
        the same additive `gentle_primer_pair_rank_v1` score decomposition.
        Primer3 near-miss capture is explicitly `incomplete`: GENtle can retain
        only Primer3-returned pairs that its post-filters reject, not Primer3's
        hidden internal rejection space.
      - retained qPCR assays use additive model
        `gentle_qpcr_assay_rank_v1` (`higher_is_better`): the primer-pair terms
        plus probe Tm-offset and amplicon-midpoint contributions. Probe
        self/primer complementarity diagnostics are first-class zero-weight
        observational terms in v1, so provenance records them without changing
        the established rank score.
      - transcript-aware qPCR capture is explicitly `incomplete`: aggregate
        transcript-local rejections are retained in the census, but rejected
        coordinates are not projected back to the source sequence.
      - every primer-pair or qPCR design report links a report-content-fingerprinted
        construct-reasoning graph. Its bounded rejected intervals reach the
        existing linear-map overlay as non-verdict `ContextEvidence`.
      - design does not currently consult repeat, variant, or
        paralogue-similarity tracks for region-level exclusion. Their absence
        from a graph means `not_run`, not checked and clear. New reports carry
        that omission explicitly as `excluded_region_analysis_status` and
        `excluded_region_analysis_reason`; legacy reports omit both.

## ROI and qPCR Seeds

- Shared shell:
  - Commands:
    - Primer ROI seed helper notes (`primers seed-from-feature` / `primers seed-from-splicing`):
      - returns non-mutating schema `gentle.primer_seed_request.v1`
      - includes `template`, source metadata, `roi_start_0based`,
        `roi_end_0based_exclusive`
      - includes ready-to-run operations:
        `operations.design_primer_pairs` (`DesignPrimerPairs`) and
        `operations.design_qpcr_assays` (`DesignQpcrAssays`)
    - qPCR-only seed helper notes
      (`primers seed-qpcr-from-feature` /
      `primers seed-qpcr-from-splicing`):
      - returns non-mutating schema `gentle.qpcr_seed_request.v1`
      - includes `template`, source metadata, `roi_start_0based`,
        `roi_end_0based_exclusive`
      - `primers seed-qpcr-from-splicing` additionally accepts:
        `--mode shared_gene|distinguish_transcript`
      - `--transcript-id ID` is required when
        `--mode distinguish_transcript`
      - `--specificity-evidence junction_only|unique_exon_or_chain|either_prefer_junction`
        is available when `--mode distinguish_transcript`
      - splicing-seeded qPCR requests now carry engine-owned
        `transcript_targeting` metadata inside the emitted
        `DesignQpcrAssays` operation, so the same transcript-aware intent is
        preserved across shell, CLI, GUI, ClawBio, JS, and Lua
      - includes `rationale.summary`, `rationale.why_this_roi`, and
        `rationale.recommended_defaults`
      - seeded qPCR defaults are probe-based/TaqMan-like: 80-200 bp amplicons,
        18-24 bp primers at 63-73 C with pair delta <= 3 C, and a 20-30 bp
        probe at 71-80 C so retained probes normally sit above the primer mean
      - includes one ready-to-run `operation` (`DesignQpcrAssays`)
      - `primers design-qpcr` now accepts either that full
        `gentle.qpcr_seed_request.v1` payload directly or the nested bare
        `{"DesignQpcrAssays": ...}` operation
      - includes built-in protocol-cartoon metadata for `pcr.assay.qpcr`, so
        shell/CLI/ClawBio flows can promote the same qPCR strip without hard-
        coding that protocol id elsewhere

- Shared shell:
  - Commands:

## Experimental Handoff

- Shared shell:
  - Commands:
    - Experimental handoff notes (`primers experimental-handoff`):
      - consumes one persisted transcript-panel report and emits one
        deterministic card per selected primer pair; it automatically runs the
        shared cDNA assay test for each pair and links it by canonical `pair_id`
        plus instance-specific `assay_test_id`
      - the default policy requires critical oligo QC, annotation provenance,
        and separate passing `genomic_carryover` and
        `transcriptome_specificity` gates. RT-PCR/qPCR therefore cannot become
        order-ready from genomic evidence alone. Missing variant evidence is
        shown as `not_evaluated` but is optional unless the supplied policy
        requires it. Evaluated failures remain visible blockers
      - `--policy` accepts inline JSON or `@FILE`; `--variant-evidence` may be
        repeated; `--order-form-id` projects existing modifications, scale,
        purification, and duplicate-review state without merging order lines
      - `--path` writes the complete JSON package and `--order-table` writes a
        compact TSV. Neither route submits an order or records a wet-lab result
      - `--gel-svg` writes a GENtle-attributed virtual gel with one sample lane
        per primer pair under shared conditions; `--gel-ladder` may be repeated.
        Lanes with no predicted cDNA product remain visible. Unless the policy
        supplies named gel conditions, the figure is illustrative and does not
        satisfy the separate gel-resolution readiness gate
      - exact cDNA product digests define sequence classes only. Gel separation
        is assessed separately and only when named gel conditions are present
        in the policy
      - the typed `BuildExperimentalAssayHandoff` operation is also available
        through generic MCP `op`, workflows, JavaScript, and Lua; the GUI Shell
        surface accepts the same command

<!-- manual-page-nav:start -->

---

[CLI Manual](../cli.md) | [Previous: Protein and Structure Evidence](protein-evidence.md) | [Next: Reference-Bound Primer Specificity](primer-specificity.md)
<!-- manual-page-nav:end -->
