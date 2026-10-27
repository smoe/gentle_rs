# Introspection, Facts and Readiness

<!-- manual-counterparts:start -->
No direct GUI counterpart: This chapter documents an automation interface, example or command-only workflow; no dedicated GUI chapter is paired with it.
<!-- manual-counterparts:end -->

- Shared shell:
  - Commands:
    - `help`
    - `capabilities`
    - `state-summary`
    - `facts graph [--evidence SCAN.json ...]`
    - `facts eval FACT_EXPR_JSON_OR_@FILE [--evidence SCAN.json ...]`
      - evaluates `gentle.fact_expression.v1`-style preconditions against the
        deterministic `gentle.project_fact_graph.v1`; persisted
        reverse-translation, protein-derivation, primer/qPCR design,
        restriction-cloning handoff, sequencing-confirmation, CUT&RUN read,
        and RNA-read interpretation reports are projected as `report.exists`;
        imported sequencing traces are projected as `sequencing_trace.exists`;
        and restriction-site absence requires explicit scan evidence rather
        than missing hits in state

## Facts and Capability Descriptors

- Shared shell:
  - Commands:
    - `introspect facts [--domain project|view|host|config] [--seq-id SEQ_ID] [--evidence SCAN.json ...] [--ui-host true|false]`
      - returns `gentle.introspection.v1` fact read-back grouped by domain
      - `--seq-id` filters facts to one concrete loaded sequence subject
      - headless CLI contexts project `ui.host_available=false` unless
        `--ui-host true` is supplied by a GUI-attached caller/test harness
    - `introspect capabilities [--kind operation|view_intent|host_config]`
      - returns the shared capability registry plus fact-annotated descriptors
        for the validated introspection slice; registry-only rows remain
        discoverable but do not invent preconditions
      - the top-level `introspect facts`, `introspect capabilities`,
        `introspect readiness`, `introspect verify-effects`,
        `introspect runtime`, and `introspect all` shell routes are themselves
        fact-annotated as
        read-only self-description/status routes with no project-state
        preconditions

## Live Runtime Inspection

- Shared shell:
  - Commands:
    - `introspect runtime`
      - returns `gentle.runtime_status.v1`, with process-local live `frames[]`
        for the current GENtle process plus observed `activities[]` from
        existing persisted/project ledgers
      - active shell commands appear as `shell_command` frames; genome
        preparation additionally mirrors current prepare phase/progress as a
        `resource_prepare` frame
      - observed activities include default-catalog genome prepare markers,
        CUT&RUN shared-asset markers, and project-backed BLAST async jobs,
        tagged as live, cross-process, completed, failed, cancelled, stale, or
        unknown
      - no runtime-status file is written; completed frames drop out of the
        stack rather than becoming a history log
      - on Unix, sending SIGUSR1 to `gentle`, `gentle_cli`, or `gentle_mcp`
        prints the process-local snapshot to STDERR

## Readiness Evaluation

- Shared shell:
  - Commands:
    - `introspect readiness [CAPABILITY_ID] [--arg NAME=VALUE ...] [--seq-id SEQ_ID] [--readiness ready|blocked|unknown] [--evidence SCAN.json ...] [--ui-host true|false]`
      - evaluates bound or unbound capability readiness through the shared fact
        evaluator; unbound templated args report `unknown`
      - fact-annotated descriptors are evaluated through their full
        `precondition_expr`, including `any` branches, not only flat `all`
        atom lists
      - partial alternatives retain three-valued truth: false OR an unbound
        alternative is `unknown`; true OR an unbound alternative is `ready`.
        Binding diagnostics remain visible but are not an independent veto
      - TSS SVG export deliberately selects the supplied exclusive `SEQ_ID`
        or `COLLECTION_ID`: a missing selected target stays `blocked`, while
        no target selection stays `unknown`
      - `--seq-id` binds `SEQ_ID`, and `--readiness` filters evaluated rows

## Readiness: Reports, Sequencing and Agents

- Shared shell:
  - Commands:
    - Readiness annotation details:
      - persisted-report inspection commands such as `primers show-report`,
        `primers show-qpcr-report`, and
        `primers show-restriction-cloning-handoff`, plus
        `seq-confirm show-report`, `cutrun show-read-report`, and
        `rna-reads show-report`, become ready when the matching `report.exists`
        fact is projected
      - raw persisted-report show operation rows
        `ShowSequencingConfirmationReport`, `ShowCutRunReadReport`, and
        `ShowRnaReadReport` mirror their shell `show-*` routes: they require
        the matching `report.exists` fact and declare no side effects
      - sequencing-trace commands project `sequencing_trace.exists(TRACE_ID)`;
        `seq-trace list` / `ListSequencingTraces` are catalog-ready,
        `seq-trace show` / `ShowSequencingTrace` require the trace fact, and
        `seq-trace import` / `ImportSequencingTrace` can verify the same fact
        when the caller supplied an explicit `TRACE_ID`
      - sequencing-confirmation run rows `seq-confirm run` and
        `ConfirmConstructReads` require the expected construct sequence
        (`sequence.exists(EXPECTED_SEQ_ID)`) and can verify
        `report.exists(REPORT_ID) == sequencing_confirmation` when the caller
        supplied a deterministic report id. Optional repeated read sequences
        and trace ids are declared as readable evidence inputs; callers should
        repeat readiness checks for each concrete id until list-valued
        introspection binding is promoted.
      - sequencing-primer overlay rows `seq-primer suggest` and
        `SuggestSequencingPrimers` require the expected construct sequence.
        Optional primer sequence ids and sequencing-confirmation reports are
        declared as readable inputs so bound clients can check them explicitly.
      - agent dispatch rows `agents ask`, `ask_agent_system`, `agents plan`,
        and `agent_plan` require
        `host.tool_available(SYSTEM_ID) == true`, matching the same host fact
        used by `agents preflight`. `agents execute-plan` and
        `agent_execute_plan` are payload-ready without project facts because
        their real precondition is a valid supplied agent-plan payload and
        candidate id; command-dependent project effects are advertised only as
        `may_on_success`.
      - read-acquisition routes `reads acquire status`, `reads acquire prepare`,
        `reads acquire inspect`, `reads acquire cancel`, and raw
        `ReadAcquire*` rows are payload/path-ready without project facts.
        Prepare and cancel rows model the external work/cache updates as
        `artifact.written(WORK_DIR)` external handoffs; readiness does not
        assert that the manifest, SRA accession, SRA Toolkit, or paths are
        valid until execution.
      - persisted-report export commands such as
        `reverse-translate export-report`, `primers export-report`, and
        `primers export-qpcr-report`, and
        `primers export-restriction-cloning-handoff`, plus
        `seq-confirm export-report`, `seq-confirm export-support-tsv`,
        `cutrun export-coverage`, `rna-reads export-report`, and specialized
        RNA-read artifact exports such as `rna-reads export-hits-fasta`,
        `rna-reads export-target-quality`, `rna-reads export-paths-tsv`,
        `rna-reads export-abundance-tsv`,
        `rna-reads export-dexseq-annotation-gff`,
        `rna-reads export-dexseq-counts-tsv`,
        `rna-reads verify-dexseq`,
        `rna-reads export-score-density-svg`,
        `rna-reads export-alignments-tsv`,
        `rna-reads export-isoform-triage-tsv`, and
        `rna-reads export-alignment-dotplot-svg`, use the same `report.exists`
        readiness and model output files as `artifact.written` external
        handoffs. The raw engine operation rows
        `ExportPrimerDesignReport`,
        `ExportSequencingConfirmationReport`,
        `ExportSequencingConfirmationSupportTsv`,
        `ExportCutRunReadCoverage`, `ExportRnaReadReport`,
        `ExportRnaReadHitsFasta`, `ExportRnaReadTargetQuality`,
        `ExportRnaReadExonPathsTsv`, `ExportRnaReadExonAbundanceTsv`,
        `ExportRnaReadDexseqAnnotationGff`, `ExportRnaReadDexseqCountsTsv`,
        `VerifyRnaReadDexseqExports`,
        `ExportRnaReadScoreDensitySvg`, `ExportRnaReadAlignmentsTsv`,
        `ExportRnaReadIsoformTriageTsv`, and
        `ExportRnaReadAlignmentDotplotSvg` expose the same model for
        registry-driven adapters. Raw RNA-read gene-support analysis rows
        `SummarizeRnaReadGeneSupport` and `InspectRnaReadGeneSupport`
        also require `report.exists(REPORT_ID) == rna_read` and model their
        optional JSON `path` as an `artifact.written` external handoff.

## Readiness: RNA and Primer Workflows

- Shared shell:
  - Commands:
    - Readiness annotation details:
      - RNA-read processing commands `rna-reads align-report` /
        `AlignRnaReadReport` require `report.exists(REPORT_ID) == rna_read`
        and verify that the report remains present after the mutating alignment
        update. `rna-reads materialize-hits` /
        `MaterializeRnaReadHitSequences` use the same report readiness but do
        not yet declare a hard sequence-creation effect because created ids are
        selection-dependent. `rna-reads preflight-isoforms` /
        `PreflightRnaReadIsoforms` require `sequence.exists(SEQ_ID)` and have
        no project-state effect.
      - primer helper readbacks are fact-annotated where the readiness can be
        expressed with existing facts: `primers preflight` is ready without
        project state, `primers seed-from-feature`,
        `primers seed-from-splicing`, and
        `primers restriction-cloning-vector-suggestions` require
        `sequence.exists(SEQ_ID)`, and
        `primers seed-restriction-cloning-handoff` requires both
        `report.exists(PRIMER_REPORT_ID) == primer_design` and
        `sequence.exists(DESTINATION_VECTOR_SEQ_ID)`. cDNA assay test routes
        remain registry-only until their conditional product-materialization
        effects are modeled.

## Readiness: Catalogs, Services and Planning

- Shared shell:
  - Commands:
    - Readiness annotation details:
      - no-project local catalog/report operations such as
        `SummarizeJasparEntries`, `BenchmarkJasparRegistry`,
        `ListJasparCatalog`, `ResolveTfQueries`, `ListReporterCatalog`, and
        `RecommendReporters` are ready in catalog mode and model optional JSON
        `path` outputs as `artifact.written` external handoffs
      - reporter shell routes `reporters list` and `reporters recommend`
        mirror the raw reporter catalog operations as catalog-ready optional
        JSON artifact handoffs; `reporters export-corpus` and raw
        `ExportReporterCorpus` are also catalog-ready and model their required
        JSON/JSONL output path as an `artifact.written` external handoff
      - service status/provider catalog routes `services status`,
        `services providers list`, `services providers doctor`,
        `services delivery-route`, `services project-preflight`,
        `services project-quote`, `services handoff`, and `services guide` are
        fact-annotated as ready without project state. Doctor, quote, and
        handoff output paths are modeled as `artifact.written` external
        handoffs. These routes prepare local classification, validation, guide,
        or handoff records only; they do not submit vendor orders.
        `services route-project-source` remains separate because its
        preconditions depend on the selected project object kind.
      - planning read-only routes `planning consult cloning`,
        `planning protein-expression-handoff`, `planning profile show`,
        `planning objective show`, `planning suggestions list`, and
        `planning sync status` are fact-annotated as ready without project
        state and declare no side effects. Planning profile/objective
        set/clear, suggestion accept/reject, and sync pull/push are also ready
        without project facts, but declare only a non-verifiable
        `may_on_success` planning-state effect because profile/objective/
        suggestion facts are not projected yet.
      - shell-level resource/catalog inspection routes such as
        `resources summarize-jaspar`, `resources status`,
        `resources sync-rebase`, `resources sync-jaspar`,
        `resources sync-ucsc-rmsk`,
        `resources import-gene-list-cache`,
        `resources import-ontology-assignment-cache`,
        `resources import-co-regulated-cache`,
        `resources install-ucsc-rmsk`,
        `resources prepare-ucsc-rmsk-index`,
        `resources sync-jaspar-remote-metadata`,
        `resources benchmark-jaspar`,
        `resources suggest-ucsc-rmsk-index`, `resources list-jaspar`,
        `resources inspect-jaspar`, `resources resolve-tf-query`,
        `resources list-publication-datasets`,
        `resources status-publication-dataset`, `genomes validate-catalog`,
        `helpers validate-catalog`, `helpers vocabulary list`, and
        `helpers vocabulary doctor` are also fact-annotated as ready without
        project state. Optional JSON output paths are modeled as
        `artifact.written` external handoffs; catalog, motif, dataset, and file
        path validation remains part of execution.
      - local cache/resource mutation routes `cache clear` and
        `resources prepare-publication-dataset` are ready without project
        facts, but declare only a non-verifiable `may_on_success` effect
        because their concrete filesystem/cache changes depend on supplied
        paths, download choices, and execution-time validation.
      - host/helper/protease/microRNA catalog helper routes such as
        `hosts list`, `list_host_profile_catalog_entries`,
        `host_profile_catalog_entries`, `list_helper_catalog_entries`,
        `helper_catalog_entries`, `helper_semantics_vocabulary`,
        `helper_interpretation`, `proteases list`, `proteases show`,
        `mirna explain-seed`, and `mirna catalog-show` are fact-annotated as
        ready without project state. Protease list/show JSON output paths are
        modeled as optional `artifact.written` handoffs; catalog lookup
        validation remains part of execution.

## Readiness: Genomes and Helper Resources

- Shared shell:
  - Commands:
    - Readiness annotation details:
      - reference/helper genome cache inspection routes (`genomes status`,
        `helpers status`, `genomes genes`, `helpers genes`) and JS/Lua helper
        adapters (`list_reference_genomes`, `list_reference_catalog_entries`,
        `is_reference_genome_prepared`, `list_reference_genome_genes`) are
        fact-annotated as ready without project state. This means the command
        can be evaluated from the current project-fact graph; concrete genome
        ids, prepared caches, catalog paths, regex filters, and biotype filters
        are still validated by the command when it runs.
      - reference/helper genome extraction routes (`genomes|helpers
        extract-region`, `extract-gene`, and `extract-promoter`) plus raw
        `ExtractGenome*` operation rows and JS/Lua helper aliases are
        fact-annotated as project-sequence creation routes. The prepared cache
        is still validated during execution; introspection can verify the
        resulting `sequence.exists(OUTPUT_ID)` fact when the caller supplied a
        deterministic output id. Anchor extension/verification routes require
        the source `sequence.exists(SEQ_ID)` fact.
      - Ensembl discovery/catalog-maintenance routes
        (`genomes ensembl-available`, `helpers ensembl-available`,
        `ensembl_installable_genomes`, `list_ensembl_installable_genomes`,
        `genomes preview-ensembl-specs`, `helpers preview-ensembl-specs`,
        `genomes update-ensembl-specs`, and
        `helpers update-ensembl-specs`) are fact-annotated as ready without
        project state. Update rows model the external catalog output as
        `artifact.written(OUTPUT_CATALOG_PATH)`; concrete catalog paths and
        Ensembl/network availability remain execution-time concerns.
      - async BLAST status/list routes (`genomes blast-status`,
        `helpers blast-status`, `genomes blast-list`, `helpers blast-list`,
        `blast_async_status`, `blast_async_list`) are fact-annotated with no
        project-state preconditions. They observe cached live or saved receipts
        without dispatching work or rewriting project metadata; job ids and
        optional terminal reports remain execution-time
        validation concerns. Async start/cancel and synchronous BLAST routes
        remain registry-only in this slice.
      - built-in ladder catalog routes such as `ladders list`,
        `inspect_dna_ladders`, `inspect_rna_ladders`, `list_dna_ladders`, and
        `list_rna_ladders` are ready without project state; ladder export
        routes `ladders export`, `export_dna_ladders`, `export_rna_ladders`,
        `ExportDnaLadders`, and `ExportRnaLadders` additionally model the JSON
        output path as an `artifact.written` external handoff
      - agent-system catalog routes `agents list`, `agent_systems`, and
        `list_agent_systems` are ready without project state because they only
        enumerate configured systems
      - `agents preflight`, `agents discover-models`, `agent_preflight`, and
        `agent_models` are fact-annotated host-config routes that require
        `host.tool_available(SYSTEM_ID) == true`. Preflight may emit an
        external-handoff report; model discovery declares no project effects.
        Asking, planning, and plan execution remain separate adapter/transport
        operations.
      - adapter parity aliases `state_summary`, `reference_catalog_entries`,
        `ui_intents`, `ui_prepared_genomes`, and `ui_latest_prepared`, plus
        shell routes `ui prepared-genomes` and `ui latest-prepared`, are
        fact-annotated with the same no-project readiness as their shared
        contracts (`state-summary`, reference catalog listing, and deterministic
        UI catalog/prepared-genome query routes).

## Readiness: Views, Tracks and Analysis Catalogs

- Shared shell:
  - Commands:
    - Readiness annotation details:
      - generic GUI intent requests (`ui open`, `ui focus`, `ui close`, and
        `ui_intent`) are fact-annotated as view intents that require
        `ui.host_available == true`. In headless command-line use they can
        still resolve to a structured intent payload, but readiness reports that
        the intent cannot be applied until a GUI host is attached.
      - protocol-cartoon catalog/render/template routes are ready without
        project state; render/export rows model SVG or JSON outputs as
        `artifact.written` external handoffs, while template input path
        validation remains part of execution rather than a project fact
      - prepared-cache inspection, CUT&RUN dataset catalog/status inspection,
        and array helper inspection routes are ready without project state;
        external file/directory/catalog validation remains part of execution,
        and `arrays render-probe-region-output-svg` models its SVG output as an
        `artifact.written` external handoff
      - genome-track import and array projection routes (`tracks import-bed`,
        `tracks import-bigwig`, `tracks import-vcf`, their raw/MCP operation
        rows, `genomes|helpers blast-track`,
        `arrays project-microarray-track`, `ProjectMicroarrayTrack`, and
        `arrays project-probe-region-output`) require the loaded target
        sequence fact. They currently remain readiness-only because feature
        freshness/track-update facts are not projected yet
      - tracked genome signal-file subscription routes
        (`tracks tracked list|add|remove|clear|apply`) have no project-fact
        preconditions. List is read-only; add/remove/clear/apply declare only a
        non-verifiable `may_on_success` project-state effect because imported
        features depend on current anchors and external file validation
      - sequence-scan reporting/rendering rows (`FindRestrictionSites`,
        `features tfbs-score-tracks-svg`, `RenderTfbsScoreTracksSvg`,
        `SummarizeTfbsScoreTracks`, `features tfbs-track-similarity`, and
        `SummarizeTfbsTrackSimilarity`) require
        `sequence.exists(SEQ_ID)` when the scan target is a loaded project
        sequence. Inline sequence text and motif/enzyme validation remain
        execution-time checks; JSON/SVG output paths are modeled as
        `artifact.written` external handoffs where present
      - catalog/list routes for candidate sets, guide sets, workflow macros,
        candidate macro templates, and routine catalogs are ready without
        project state; routes that inspect a named persisted set/template still
        remain registry-only until those object-existence facts are projected
      - construct-reasoning graph list routes
        `construct-reasoning list-graphs` and `construct_reasoning_graphs` are
        ready without project state; their optional `SEQ_ID` is treated as a
        filter. Named graph inspection/action routes require
        `construct_reasoning_graph.exists(GRAPH_ID)`; status/writeback actions
        verify that the graph remains present, and graph export models the JSON
        output as an `artifact.written` external handoff
      - persisted analysis-payload list routes `dotplot list` and `flex list`
        are ready without project state; optional `SEQ_ID` is treated as a
        filter, while show/render-by-id routes remain registry-only until
        dotplot/flex-track existence facts are projected
      - local metadata/catalog routes `genomes list`, `helpers list`,
        `ensembl-gene list`, `ensembl-protein list`, and
        `gene-groups list|show|resolve|doctor` are ready without project
        state; catalog/group/token validation remains part of execution, and
        gene-group JSON outputs are modeled as optional `artifact.written`
        external handoffs
      - report-list commands such as `primers list-reports`,
        `primers list-qpcr-reports`,
        `primers list-restriction-cloning-handoffs`, and
        `reverse-translate list-reports`, plus `cutrun list-read-reports` and
        `rna-reads list-reports`, have no project preconditions and are ready
        in catalog mode. The raw engine operation rows
        `ListSequencingConfirmationReports`, `ListCutRunReadReports`, and
        `ListRnaReadReports` follow the same catalog-ready convention.

## Readiness: Read Reports, Facts and Settings

- Shared shell:
  - Commands:
    - Readiness annotation details:
      - RNA-read report inspection routes (`rna-reads show-alignment`,
        `rna-reads show-alignments`, `rna-reads summarize-gene-support`,
        `rna-reads inspect-gene-support`, `rna-reads inspect-alignments`, and
        `rna-reads inspect-concatemers`) require a persisted `rna_read`
        `report.exists(REPORT_ID)`. Routes with an optional `--output` path
        model it as an `artifact.written` external handoff.
      - self-description/status commands such as `help`, `capabilities`,
        `state-summary`, and `history status` also have no project
        preconditions and are ready in catalog mode
      - fact-layer commands `facts graph` and `facts eval` are also
        fact-annotated, no-mutation routes; `facts eval` validates its supplied
        expression at execution time rather than as a project-state
        precondition
      - engine-owned configuration parameters are projected as closed-world
        `config.param` facts under `introspect facts --domain config`; the
        subject id is the parameter name and the value is the current JSON
        value
      - `set-param` and the lower-case adapter alias `set_parameter` are
        fact-annotated `host_config` capabilities with no project preconditions
        and a `must_on_success` effect that verifies
        `config.param(PARAM_NAME) == PARAM_VALUE`; the raw `SetParameter`
        operation row exposes the same effect for registry-driven adapters
      - view/display commands such as `display show|hide|visibility` and the
        raw `SetDisplayVisibility` operation row have no project preconditions
        and declare a `view.visible_tracks` `view_session` effect;
        `introspect facts --domain view` projects the current display-layer
        booleans as closed-world read-back state

## Readiness: Sequence Design and Physical Layouts

- Shared shell:
  - Commands:
    - Readiness annotation details:
      - direct sequence-derivation engine operations such as `Reverse`,
        `Complement`, `ReverseComplement`, `Branch`, and `ExtractRegion` are
        fact-annotated through their raw operation rows: they require
        `sequence.exists(INPUT_SEQ_ID)` and verify
        `sequence.exists(OUTPUT_ID)` when a deterministic output id is supplied
      - `MaterializeVariantAllele` mirrors `variant materialize-allele`:
        it requires `sequence.exists(INPUT_SEQ_ID)` and verifies
        `sequence.exists(OUTPUT_ID)` for deterministic allele materialization
      - `AlignSequences` mirrors `align compute`: it requires both
        `sequence.exists(QUERY_SEQ_ID)` and `sequence.exists(TARGET_SEQ_ID)`;
        it declares no project-state effects because the alignment result is
        returned directly
      - raw primer/qPCR design operation rows `DesignPrimerPairs`,
        `DesignInsertionPrimerPairs`, and `DesignQpcrAssays` mirror their shell
        design routes: they require `sequence.exists(TEMPLATE_SEQ_ID)` and
        verify the requested `report.exists(REPORT_ID)` kind when a
        deterministic report id is supplied
      - dotplot and flexibility routes project stored payload facts:
        `dotplot.exists(DOTPLOT_ID)` and
        `flexibility_track.exists(TRACK_ID)`. `dotplot compute`,
        `ComputeDotplot`, `dotplot overlay-compute`,
        `ComputeDotplotOverlay`, `flex compute`, and
        `ComputeFlexibilityTrack` require loaded sequence inputs and can verify
        those payload facts when deterministic ids are supplied. `dotplot show`,
        `flex show`, and `RenderDotplotSvg` become ready only once the stored
        payload fact exists; `RenderDotplotSvg` models the SVG path as an
        `artifact.written` external handoff
      - candidate-window optimization routes project
        `candidate_set.exists(SET_NAME)`. Candidate generation requires a
        loaded source sequence and verifies the named set when deterministic
        `SET_NAME` is supplied. Candidate show/metrics/score/delete routes
        require the existing set; filter/top-k/pareto/set-op routes require
        input set facts and verify the named output set. Deletion currently has
        no positive hard effect because the fact language does not yet model
        absence effects
      - guide-design routes project `guide_set.exists(GUIDE_SET_ID)`,
        `guide_filter_report.exists(GUIDE_SET_ID)`, and
        `guide_oligo_set.exists(OLIGO_SET_ID)`. Guide-set upsert and oligo
        generation can verify deterministic output ids; practical filtering
        verifies the filter report and, when supplied, the filtered output guide
        set. Guide exports require the guide set and model CSV/FASTA/protocol
        files as `artifact.written` external handoffs
      - container/rack authoring routes project
        `container.exists(CONTAINER_ID)`, `arrangement.exists(ARRANGEMENT_ID)`,
        and `rack.exists(RACK_ID)`. Container exclusivity updates require the
        container fact. Serial arrangement and rack creation can verify supplied
        deterministic ids. Rack move/profile/template/block commands require
        and preserve `rack.exists`; rack SVG/OpenSCAD/simulation exports require
        the rack and model files as `artifact.written` external handoffs.
        Single-container transform operations (`DigestContainer`,
        `LigationContainer`, and `FilterContainerByMolecularWeight`) require
        `container.exists(CONTAINER_ID)` and declare `may_on_success` product
        effects because concrete product ids are derived during execution
      - macro-template routes project
        `workflow_macro_template.exists(TEMPLATE_NAME)`,
        `candidate_macro_template.exists(TEMPLATE_NAME)`, and
        `macro_instance.exists(MACRO_INSTANCE_ID)`. Template upsert can verify
        deterministic template names. Template show/delete/run require the
        matching template fact; delete routes do not yet declare absence
        effects. `macros instance-show` requires a recorded macro-instance fact

## Readiness: Project and External Sequence Operations

- Shared shell:
  - Commands:
    - Readiness annotation details:
      - project/session utility routes are fact-annotated without existing
        project facts: `history undo`, `history redo`, `load-project`, and
        `load_project` declare `may_on_success` effects because the concrete
        state changes depend on the undo/redo stack or loaded project file.
        `save-project` and `save_project` model their external project-state
        output as `artifact.written(OUTPUT_PATH)`. `load_dna` mirrors
        `LoadFile` for JS/Lua adapters and can verify
        `sequence.exists(OUTPUT_ID)` only when the adapter supplies a stable
        output id.
      - external sequence creation routes are fact-annotated for both shell and
        raw operation callers: `LoadFile`, `load_dna`, `genbank fetch`,
        `FetchGenBankAccession`, `ensembl-region fetch`,
        `FetchEnsemblRegion`, `dbsnp fetch`, `FetchDbSnpRegion`,
        and `FetchUniprotLinkedGenBank`. These fetch/import routes have no
        existing-project preconditions; `introspect verify-effects` can verify
        `sequence.exists(OUTPUT_ID)` only when the caller supplied the explicit
        `--as-id`, `as_id`, `--output-id`, or `output_id` that execution created
      - protein/gene metadata routes project
        `uniprot_entry.exists(ENTRY_ID)`,
        `ensembl_gene_entry.exists(ENTRY_ID)`, and
        `ensembl_protein_entry.exists(ENTRY_ID)`. UniProt and Ensembl
        fetch/import routes can verify explicit entry ids; show routes and
        metadata-backed sequence imports require the matching stored metadata
        fact. `ImportUniprotEntrySequence`, `ensembl-gene import-sequence`,
        `ImportEnsemblGeneSequence`, `ensembl-protein import-sequence`, and
        `ImportEnsemblProteinSequence` can also verify
        `sequence.exists(OUTPUT_ID)` when a deterministic output id was supplied
      - UniProt projection routes project closed-world
        `uniprot_projection.exists(PROJECTION_ID)` facts. `uniprot map` and
        `ProjectUniprotToGenome` require the stored UniProt-entry fact plus the
        target `sequence.exists(SEQ_ID)` fact and can verify a supplied
        deterministic `PROJECTION_ID`. Projection-show, feature-coding,
        Ensembl-link resolution, transcript accounting, Ensembl comparison, and
        audit-generation routes require the projection fact. UniProt projection
        audit and audit-parity reports still project as
        `report.exists(REPORT_ID)` facts, so their show/export routes can be
        readiness-checked.

## Readiness: Core and Derived Sequence Operations

- Shared shell:
  - Commands:
    - Readiness annotation details:
      - raw core sequence operation rows are fact-annotated where the state
        contract is deterministic: `SaveFile` requires `sequence.exists(SEQ_ID)`
        and declares an `artifact.written(OUTPUT_PATH)` external handoff;
        `Digest` and the glossary `digest` alias require
        `sequence.exists(INPUT_SEQ_ID)` but are readiness-only because fragment
        ids are prefix/index-derived; `ExtractAnchoredRegion` requires the same
        input sequence fact and is readiness-only because candidate ids are
        output-prefix/rank-derived; `SelectCandidate` also requires
        `sequence.exists(INPUT_SEQ_ID)` and can verify
        `sequence.exists(OUTPUT_ID)` when a deterministic selected-candidate id
        was supplied; `Pcr`, `PcrAdvanced`, and `PcrMutagenesis` require
        `sequence.exists(TEMPLATE_SEQ_ID)` and can verify
        `sequence.exists(OUTPUT_ID)` for deterministic product ids; and
        `PcrOverlapExtensionMutagenesis` requires the template sequence but is
        readiness-only until generated candidate ids are projected as facts
      - `render-dotplot-svg` is fact-annotated with the same readiness contract
        as `RenderDotplotSvg`: it requires `sequence.exists(SEQ_ID)` and
        `dotplot.exists(DOTPLOT_ID)`, then models the SVG path as an
        `artifact.written(OUTPUT_PATH)` external handoff
      - raw cDNA assay-test rows `TestCdnaPcr` and `TestCdnaQpcr` are
        fact-annotated over `sequence.exists(SEQ_ID)` and advertise only
        conservative optional-output effects. The shell routes
        `primers test-cdna-pcr` and `primers test-cdna-qpcr` expose the same
        readiness. `TestCdnaQpcrFasta`, `primers test-cdna-qpcr-fasta`, and
        `ProjectGenomeInterval` have no project-state preconditions because
        their FASTA/projection-map paths are external files validated during
        execution
      - `ReverseTranslateProteinSequence` mirrors `reverse-translate run`:
        the input sequence must exist and have `sequence.kind == protein`, and
        a deterministic `OUTPUT_ID` can be verified as the generated coding-DNA
        sequence
      - raw transcript/protein/splicing derivation rows
        `DeriveTranscriptSequences`, `DeriveProteinSequences`, and
        `DeriveSplicingReferences` require `sequence.exists(SEQ_ID)` and
        declare conservative `may_on_success` sequence-creation effects because
        output ids depend on feature ids, prefixes, and uniqueness rules.
        `DeriveProteinSequences` can additionally verify
        `report.exists(REPORT_ID) == protein_derivation` when a deterministic
        report id was supplied.
      - protease catalog and digest routes are fact-annotated:
        `proteases list` and `proteases show` are no-precondition catalog
        reads; `proteases digest`, `ProteaseDigestProteinSequence`, and
        `proteases digest-gel-svg` require an existing protein-kind sequence.
        Digest materialization advertises only a `may_on_success` sequence
        effect because peptide ids are prefix/index-derived. Protein gel SVG
        routes require `report.exists(REPORT_ID) == protein_derivation`, while
        `RenderProteaseDigestGelSvg` can also be ready from an explicit
        protein sequence. SVG outputs are modeled as `artifact.written`
        external handoffs.

## Readiness: Topology, Views and Export

- Shared shell:
  - Commands:
    - Readiness annotation details:
      - `SetTopology` is fact-annotated through its raw operation row: it
        requires `sequence.exists(SEQ_ID)` and verifies
        `sequence.circular(SEQ_ID) == CIRCULAR`
      - `RecomputeFeatures` is fact-annotated through its raw operation row as
        readiness-only: it requires `sequence.exists(SEQ_ID)`, but it does not
        yet declare a hard fact effect because computed-feature freshness is
        not projected as a fact
      - `SetLinearViewport` is fact-annotated through its raw operation row:
        it has no project preconditions and verifies
        `view.viewport(linear_sequence) == {"start_bp": START_BP, "span_bp": SPAN_BP}`
      - `ui intents` is a no-precondition catalog route for discovering GUI
        intent targets and commands
      - external-render commands such as `render-svg`, `render-rna-svg`, and
        `render-lineage-svg` model output files as `artifact.written` external
        handoffs rather than saved project state; sequence renderers become
        ready when their input `sequence.exists` fact is projected. The raw
        engine operation rows `RenderSequenceSvg`, `RenderRnaStructureSvg`,
        `RenderFeatureExpertSvg`, `RenderTfbsScoreTrackCorrelationSvg`, and
        `RenderLineageSvg` expose the same readiness and external-handoff model
        for registry-driven adapters.
      - read-only sequence inspectors such as `rna-info`, `features query`, and
        `features tfbs-summary` become ready when the inspected
        `sequence.exists` fact is projected. The raw `SummarizeTfbsRegion`
        and `QueryProteinResidueGenomicCoordinates` operation rows expose the
        same sequence-readiness model. Variant promoter/reporter inspectors
        (`variant promoter-context`, `variant reporter-fragments`, and their
        raw operation rows) use the same sequence-readiness model and treat
        optional JSON paths as `artifact.written` external handoffs, not
        persistent project reports. `variant annotate-promoters`,
        `AnnotatePromoterWindows`, and `AnnotateTfbs` are sequence-gated
        mutating annotations without a hard fact effect until feature
        freshness/count facts exist;
        `inspect-feature-expert` uses the same sequence-level readiness until
        target-specific facts exist
      - external feature exports such as `features export-bed` use the same
        sequence readiness and model the BED path as an `artifact.written`
        external handoff. The raw `ExportFeaturesBed` operation row exposes
        the same sequence-readiness and BED-handoff model; `render-feature-expert-svg`
        does the same for its SVG path
      - raw sequence-context operations `InspectSequenceContextView` and
        `ExportSequenceContextBundle` use the same sequence readiness; bundle
        export models its output directory as an `artifact.written` external
        handoff

## Effect Verification and Complete Introspection

- Shared shell:
  - Commands:
    - `introspect verify-effects CAPABILITY_ID [--arg NAME=VALUE ...] [--seq-id SEQ_ID] [--evidence SCAN.json ...] [--ui-host true|false]`
      - verifies `effect_kind: "must_on_success"` rows for a fact-annotated
        capability against current project facts and supplied evidence
    - `introspect all [--evidence SCAN.json ...] [--ui-host true|false]`
      - returns facts, registry-backed capabilities, and catalog-mode readiness
        in one payload

<!-- manual-page-nav:start -->

---

[CLI Manual](../cli.md) | [Previous: Shared Shell](shell.md) | [Next: Protocol Examples and Tutorials](examples.md)
<!-- manual-page-nav:end -->
