# Transcript Assay Panels

<!-- manual-counterparts:start -->
Related GUI reference: [Primer, qPCR, and Sequence-Specific RT Design](../gui/primer-qpcr-and-sequence-specific-rt-design.md).
These topic links do not imply identical feature coverage across interfaces.
<!-- manual-counterparts:end -->

- Shared shell:
  - Commands:

## cDNA Assay Tests

- Shared shell:
  - Commands:
    - cDNA PCR/qPCR assay test notes
      (`primers test-cdna-pcr` / `primers test-cdna-qpcr` /
      `primers test-cdna-qpcr-fasta`):
      - returns non-mutating schema `gentle.cdna_assay_test_report.v1`
      - derives transcript cDNA templates from the selected splicing group
        using the zero-based `FEATURE_ID`
      - `primers test-cdna-qpcr-fasta` instead streams one or more transcript
        FASTA/FASTA.gz files directly, so callers can screen complete Ensembl
        cDNA and ncRNA catalogs without importing them into project state
      - `--transcript-id ID` narrows the test to one transcript product
      - transcript-derived tests default to Ensembl transcript-id row order;
        `--transcript-order genomic_first_exon` orders rows by the first
        transcript exon in strand-aware genomic direction,
        `genomic_last_exon` orders by the last exon genomic coordinate, and
        `antisense_first_exon` inverts the first-exon order for antisense
        interpretation
      - transcript-derived tests default to `--map-coordinate-mode cdna`, where
        each row uses its transcript-local cDNA length as the axis; use
        `--map-coordinate-mode genomic_aligned` to draw every row against the
        shared source/genomic coordinate range so common primer loci and
        skipped introns line up across isoforms
      - `--max-mismatches` defaults to exact/IUPAC-compatible matching
      - `--require-3prime-exact-bases` gates primer hits only
      - qPCR probes may bind either cDNA orientation but must fall inside the
        primer-bounded amplicon interior
      - returned reports classify genomic-DNA carryover risk using junction
        primer/probe evidence, mapped genomic-equivalent amplicon length, and
        the configured max amplicon window
      - FASTA screens count every scanned record but omit non-detected
        transcript rows from the returned report unless `--transcript-id` is
        supplied
      - every report carries a `transcript_map` SVG
        (`gentle.cdna_assay_transcript_map.v1`) that places functional
        products, symbolic forward/reverse primer hits, probe hits,
        group/source exon labels (`E1`, `E2`, ...) ordered by transcript
        strand, source-exon identity colors/patterns, and exon-junction ticks
        on each transcript cDNA coordinate axis; requested forward/reverse
        primer and probe sequences are shown once in the legend, forward primer
        glyphs are one-sided above the cDNA axis, reverse primer glyphs are
        one-sided below it, probe glyphs are one-sided on the detected
        probe-binding strand, crowded maps auto-wrap into multiple columns, and
        SVG tooltips retain the transcript-local exon ordinal when it differs
        from the group label
      - in `genomic_aligned` mode the same SVG uses source/genomic coordinates
        for exon blocks, primer/probe glyphs, and amplicon source spans; dashed
        product bridges mark spliced cDNA products whose source ranges skip
        genomic sequence while length labels remain cDNA amplicon lengths
      - `--path OUTPUT.json` writes the same structured report returned on
        stdout
      - `--svg OUTPUT.svg` writes the embedded transcript-map SVG as a
        ClawBio/OpenClaw-friendly graphical artifact; shell output also
        advertises it in `preferred_artifacts[]`
      - transcript-derived `test-cdna-pcr` / `test-cdna-qpcr` remain
        non-mutating by default; add `--materialize-products` to turn detected
        products into deterministic linear sequence entries and one
        singleton/pool container
      - materialization is idempotent: repeating the same assay with the same
        prefix reuses matching product sequences/container and reports created
        vs reused ids in the structured result instead of duplicating product
        vials
      - `--product-gel-svg OUTPUT.svg` implies product materialization and
        renders the product container through the existing pool-gel renderer,
        so non-specific cDNA PCR/qPCR products appear as multiple bands in one
        lane; repeat `--product-gel-ladder NAME` to add ladder lanes. The same
        result includes `gel_band_rows[]` and `gel_summary_lines[]` for
        terminal/Telegram displays that should not rely on SVG text rendering
        alone.
      - report, transcript-map SVG, and product-gel SVG paths create parent
        directories automatically when needed
      - when no products are detected, the assay still succeeds with
        `not_detected`, writes no product vial/gel, and reports that outcome in
        the materialization summary

## Transcript qPCR Panels

- Shared shell:
  - Commands:
    - Transcript qPCR panel notes
      (`primers transcript-qpcr-panel`):
      - this is the compatibility route for fitting a characteristic forward
        primer to one stored shared reverse/probe assay; use
        `primers design-transcript-assay-panel` for matrix-centered all-transcript
        design
      - consumes a stored shared-gene qPCR report and the zero-based source
        `FEATURE_ID`
      - returns `gentle.transcript_qpcr_panel.v1` with shared forward/reverse/
        probe records plus one row per transcript
      - per-transcript rows try to keep the shared reverse primer and probe
        while selecting a characteristic forward primer; exon-junction primers
        are preferred when uniquely informative, but single-exon/exon-chain
        primers are allowed and reported explicitly
      - source coordinates are always local 0-based/exclusive plus display
        1-based/inclusive; genomic 1-based/inclusive coordinates are included
        when the source sequence has a genome anchor
      - byte-identical mature cDNAs receive
        `not_distinguishable_between_members`; other unresolved transcripts
        receive deterministic `not_found` rows instead of silent omission

## Transcript Assay Panel V2

- Shared shell:
  - Commands:
    - Transcript assay panel v2 notes
      (`primers design-transcript-assay-panel`):
      - preflight an exact endpoint operation with
        `primers inspect-transcript-assay-feasibility OPERATION_JSON_OR_@FILE
        [--path OUTPUT.json]`. The read-only command emits
        `gentle.transcript_assay_panel_feasibility.v1`, keeps the operation
        digest, stable end-class/reaction ids, effective primer/product
        constraints, and deterministic geometry blockers, and never starts
        Primer3. MCP exposes the same non-mutating route as
        `transcript_assay_feasibility`, accepting the exact operation object or
        its JSON string
      - the same preflight now derives
        `gentle.transcript_assay_primer_search_plan.v1`: narrow annotation-
        anchored left/right intervals, rejected sequence regions with reasons,
        deterministic candidate/pair estimates, and the exact
        `SEQUENCE_PRIMER_PAIR_OK_REGION_LIST`, junction/target, exclusion, and
        product-range fields that a subsequent Primer3 run will receive. It
        returns `search_space_too_broad` or `no_admissible_regions` before
        launching Primer3 when a required target cannot fit the declared
        budget. The complete operation-JSON form can set the optional
        `search_policy`; concise CLI flags retain its documented defaults. An
        advanced request may set
        `search_policy.cdna_similarity_map = {"path":"...json",
        "expected_sha256":"sha256:..."}` to bind a
        `gentle.transcript_assay_cdna_similarity_map.v1` report. GENtle
        validates the map against current mature-cDNA digests and uses
        `deprioritize` overlaps only to run lower-risk bounded records first;
        emitted Primer3 constraints and mandatory post-design complete-cDNA
        and genomic specificity gates are unchanged
      - `search_policy.build_cdna_similarity_map` instead asks GENtle to build
        and content-bind that advisory map from a prepared transcriptome/cDNA
        resource before Primer3. It is mutually exclusive with a map file;
        read-only feasibility never launches this external BLAST step
      - `search_policy.interval_evidence` accepts explicit hard exclusions,
        soft deprioritization, and informative intervals. Optional project
        projection can add known-variant and RepeatMasker evidence. The same
        plan records every exact `SEQUENCE_EXCLUDED_REGION` emitted to Primer3
      - `search_policy.primer3_chemistry` selects executable defaults, an
        explicitly emitted standard-PCR computational baseline, or restricted
        custom Primer3 solution/Tm overrides. It records calculation inputs and
        does not claim the wet-lab master mix uses them
      - `search_policy.max_junction_single_side_match_bp` can cap either
        contiguous side of a junction-spanning primer in addition to the
        existing minimum 3-prime and 5-prime overlaps
      - `search_policy.min_intervening_intron_bp` can require either a
        junction-spanning primer or an annotated intervening intron of at least
        the requested size. This reduces genomic-carryover risk by design but
        does not replace the genomic specificity gate
      - the command accepts either the concise `SEQ_ID FEATURE_ID` flag form or
        a complete externally tagged `DesignTranscriptAssayPanel` or
        `DesignTranscriptAssayPanelWithFallback` operation as inline JSON /
        `@FILE`. The operation form exposes every shared primer, probe, pair,
        junction, assay, fallback-policy, and provenance field without adding
        a CLI-only request schema
      - that same operation object is accepted unchanged by `gentle_cli op`,
        workflow JSON, MCP tool `op` (`confirm=true`), JavaScript
        `apply_operation`, and Lua `apply_operation`; all routes execute the
        shared engine operation
      - `--assay-kind` selects `endpoint-rt-pcr`, primer-only `sybr-qpcr`, or
        `taqman-qpcr`; omission preserves the original TaqMan-compatible
        forward/reverse/probe behavior for existing callers
      - endpoint mode requires `--objective isoform-end-matrix`, constructs
        only annotated first-end x terminal-end combinations, defaults to a
        10,000 bp ceiling, and reports the predicted band sizes for every
        selected assay against every transcript. When
        `--max-assays-per-class` is omitted, endpoint searches ask Primer3 for
        at most four candidate pairs per reaction; callers can still approve a
        different explicit limit
      - SYBR mode never fabricates an internal probe; its selected assays keep
        primer Tm/structure checks and the same full transcript product matrix
      - explicit junctions can be supplied with `--junctions` in transcript-
        local, genomic-intron-span, or exon-ordinal coordinates; repeat
        `--junction-evidence` to turn every Clariom-style JUC mapping in an
        interpretation report into an audited required/preferred target
      - Primer3 receives `SEQUENCE_OVERLAP_JUNCTION_LIST` plus the configured
        3-prime/5-prime overlap minima; GENtle independently records whether
        the selected forward primer, reverse primer, or neither spans each
        requested junction
      - defaults to `--objective pan-transcript --coverage-policy require-all`;
        an assay that omits any exact cDNA class is rejected rather than saved
      - `--coverage-universe JSON_OR_@FILE` independently limits the mandatory
        annotation set without removing any objective. Omission keeps all
        annotated cDNA classes. An explicit universe names transcript ids; a
        UniProt-supported universe names content-bound
        `gentle.uniprot_projection_audit.v1` sources from which GENtle derives
        the complete `(entry id, named protein-isoform id)` inventory. All
        mapped transcripts are retained beneath each protein target, but only
        one audit-ranked representative per target seeds automatic searches.
        A named isoform with no current annotation mapping stops before
        Primer3, including under `best-effort`; an entry with no named
        alternative products contributes one canonical protein target
      - for authoritative linked-transcript accounting, include both source
        kinds in the UniProt universe: the projection audit defines the
        protein-isoform targets and a
        `gentle.uniprot_linked_transcript_inventory.v1` preserves every linked
        transcript across primary, patch, and alternate loci. GENtle validates
        both the file hash and the inventory's internal content hash, then
        rejects assembly or annotation-release mismatches before Primer3.
        Exact mature-cDNA digest identity can establish cDNA coverage for an
        off-locus record, but its genomic carryover specificity remains
        unassessed until that locus is searched. Audit-only requests remain
        compatible but cannot claim an authoritative linked-record denominator
      - example UniProt universe file:
        `{"kind":"uniprot_supported_isoforms","sources":[{"source_kind":"uniprot_projection_audit","path":"reports/G1.uniprot_projection_audit.json","expected_sha256":"sha256:..."},{"source_kind":"uniprot_linked_transcript_inventory","path":"reports/G1.uniprot_linked_transcript_inventory.json","expected_sha256":"sha256:..."}]}`.
        Study-request normalization fills and binds the observed hash and the
        exact `required_uniprot_isoforms` list; direct design execution requires
        both. Regenerate legacy projection audits that lack
        `protein_isoform_inventory`. The report shows target coverage and
        distinction separately, and states explicitly that UniProt
        prioritization is not protein-expression evidence. GenBank/Ensembl CDS
        disagreements remain audit evidence rather than being silently merged
      - under this universe the existing objectives use protein targets as the
        mandatory units. Every mapped transcript/exact-cDNA class remains in
        the detection matrix, so a complete protein-target panel may still list
        informational uncovered transcript classes. Use an all-annotation or
        explicit-transcript plan when those models must each be mandatory
      - preserve alternatives as separate named requests: all-annotation,
        UniProt-prioritized, and manually explicit universes can each use
        pan-transcript, one-per-class, discrimination, or endpoint objectives
      - `--coverage-policy best-effort` is an explicit opt-in that saves a
        visibly partial report with uncovered class ids and warnings
      - `--objective maximally-informative-panel` uses the versioned default
        informative policy in the concise form. For pre-approved automatic
        recovery, use the complete wrapper operation instead: keep the nested
        strict operation at `coverage_policy=require_all`, set
        `fallback_submission.mode=preapproved_informative_partial`, and approve
        its whole-panel `max_assays`, minimum resolvable product-size difference,
        nonzero-increment rule, long-range permission, and distinct fallback
        report id. Only typed strict coverage infeasibility can submit that
        second operation. The command returns the panel when one completed and
        the linked `fallback_execution`; list/show/export inspect the persisted
        audit without rerunning design
      - anomalous wall-clock behavior does not silently select biological
        `best-effort`: GENtle may stop one bounded Primer3 record and continue
        other already planned records, but the requested coverage policy is
        unchanged and `require-all` still rejects any resulting uncovered
        class or endpoint reaction
      - strict endpoint matrices now reject annotation/geometry impossibilities
        before Primer3. Explicit best-effort matrices skip only those proven
        impossible reactions and retain each reaction, equivalence class, and
        blocker in the canonical report; feasible reactions can still end as
        `primer_search_exhausted`, which is distinct from structural
        impossibility and timeout/interruption
      - `--assay-tier` records experimental purpose independently of the
        objective. `routine-common-region-screen` requires `pan-transcript`
        and accepts commonality only from transcript annotation; Clariom PSR
        intensity can support the selected region but cannot establish that it
        is shared
      - the paired preferred-range flags add a routine-length preference inside
        the allowed `--min-amplicon-bp/--max-amplicon-bp` range. Biological
        coverage remains primary; short allowed products are explicit
        `allowed_nonpreferred` choices, while only products above the preferred
        maximum are described as `long_range_fallback`
      - pair summaries retain separate PSR/JUC evidence, one concise selection
        explanation, and a bounded deterministic list of considered
        alternatives. Endpoint gel intensity is rough/semi-quantitative, not
        a quantitative expression measurement
      - exact cDNA equivalence means byte identity only; near matches remain
        separate classes even when a single assay is difficult to design
      - each selected assay is shown against every transcript in a typed
        product matrix, and exact-hit absence is used only as a certain-negative
        shortcut when mismatch tolerance is zero
      - reports carry `gentle.transcript_assay_panel.v2` and persist through the
        existing project save/load path; list/show/export use the commands above
      - every newly selected assay also carries a redundant
        `gentle.primer_pair_summary.v2` communication view. It copies the
        forward/reverse sequences in 5-prime-to-3-prime orientation, oligo
        lengths, GENtle/Primer3 melting temperatures (`tm_c`), GC fraction and
        unrounded percent, pair `tm_delta_c`, design-transcript binding
        coordinates, the canonical designed amplicon coordinates/length,
        transcript product rows, concise oligo-QC reasons derived from the
        stored pair rule flags and metrics, junction matches,
        whole-genome specificity status, and GENtle/backend/Primer3 provenance.
        `predicted_amplicon_lengths_bp` is instead the deduplicated union from
        cross-transcript detection and may be empty even when the designed
        amplicon is present. Summary generation does not rerun sequence or
        thermodynamic analysis.
        Existing panel reports are enriched when shown or exported
      - summary v2 adds immutable sequence-derived primer ids, transcript/exon-
        aware display labels, optional aliases, typed de-novo/literature/lab
        origin, optional anchor/companion role, the satisfied design objective,
        structured selection-reason rows (`code`, `message`, `related_ids`),
        and the annotation release/reference transcript used for exon numbering
      - primer-junction overlap, amplicon-junction overlap, and selection due
        to junction evidence are three separate fields. Clariom/JUC inputs are
        preserved as repeatable `probe_region_influenced` rows with platform,
        feature/region/junction, contrast/statistic, priority, and source hash;
        they are never described as exact probe-sequence reuse
      - an evidence-free panel says
        `de_novo_no_external_selection_evidence`. A report created before the
        v2 context was stored asks for a rerun rather than inventing exon labels
        or external provenance
      - the summary does not invent a PCR annealing temperature. In the absence
        of an explicit chemistry/polymerase protocol model, experimental Ta is
        undetermined. Per-transcript genomic carryover geometry remains in the
        cDNA assay-test report; the panel summary reports `not_evaluated` rather
        than inferring a pass from junction placement
      - reports also expose order-ready primer rows, endpoint reaction/band
        matrices, a separate short-SYBR junction table, unresolved targets with
        reasons, annotation/JUC provenance, and prepared-genome specificity
        follow-up commands
      - with `--cdna-synthesis oligo-dt`, long endpoint reports warn that 5-prime
        reverse-transcription completeness can be limiting even when the PCR
        polymerase accepts 10 kb products; a missing long band is not called
        isoform absence without short end-specific assays or 5-prime RACE
      - deterministic worked example:
        `docs/examples/workflows/patz1_endpoint_sybr_transcript_assay_panel_offline.json`
        (synthetic sequence; not orderable human PATZ1 primers)

## Gene Isoform-Assay Studies

- Shared shell:
  - Commands:
    - Gene isoform-assay study planner:
      - PCR Designer's `Gene assay study` mode uses these same commands; it
        does not add a second planner or publication narrative. The bounded
        [GUI/Shell tutorial](../tutorial/04-08_gene_assay_study_gui.md) includes an
        offline preparation helper, explicitly pending dossier and two reviews.
      - `primers plan-gene-isoform-study REQUEST_JSON_OR_@FILE [--normalize-only] [--normalized-request OUTPUT.json] [--path PLAN.json] [--workflow WORKFLOW.json]`
      - normalization content-binds the complete request and all effective
        defaults before the first approval. With `--normalize-only`, stdout is
        the normalized request itself and can be persisted verbatim; planning emits a transparent
        recommendation plus exact ordered assay operations for a separate
        batch approval, but does not execute design
      - optional `fallback_submission` is part of that normalized first-stage
        approval. `mode=never` is the legacy/default behavior. With
        `preapproved_informative_partial`, eligible strict discrimination steps
        are emitted as exact `DesignTranscriptAssayPanelWithFallback` objects;
        their policy and derived fallback report id are also covered by the
        operation-batch and workflow digests used for second-stage approval
      - an optional request-level `coverage_universe` is resolved before the
        first approval and copied unchanged into every emitted assay operation.
        Keep separate named plans when comparing all annotated cDNA classes,
        UniProt-supported isoforms, and a manual transcript set; GENtle does
        not discard the other planning alternatives
      - the plan's `approved_workflow_sha256` binds the exact emitted workflow
        bytes. After the second approval, run
        `primers execute-gene-isoform-study-workflow PLAN.json
        OPERATIONS.workflow.json`; GENtle also rechecks
        `operation_batch_sha256` before executing any operation and returns
        `gentle.gene_isoform_assay_study_workflow_execution.v1`
      - for one second-stage approval across multiple studies, create a request
        with ordered `{plan_path, workflow_path}` entries and run
        `primers compose-gene-isoform-study-workflow-batch REQUEST_JSON_OR_@FILE`.
        Review and approve the emitted
        `gentle.gene_isoform_assay_study_workflow_batch.v1`, then run
        `primers execute-gene-isoform-study-workflow-batch BATCH_JSON_OR_@FILE`.
        GENtle validates the complete ordered set before the first operation,
        executes against a detached project snapshot, and commits the project
        only after the entire batch succeeds; the single-study route remains
        available
      - add `--checkpoint-dir DIR` to retain a content-bound checkpoint after
        each successful approved operation even if a later operation fails.
        Checkpoints do not alter the live project and are never reused
        automatically
      - `--on-gene-failure abort|continue` (default `abort`) selects what a
        failing gene does to the rest of the batch. `continue` moves atomicity
        down one level rather than removing it: the failing gene is rolled back
        to the boundary at which it started, so a gene is never committed
        half-applied, and GENtle then continues with the other precomputed
        genes. The report names every skipped gene in `gene_failures[]`, sets
        `batch_complete: false`, and reports `execution_atomicity: "per_gene"`.
        A batch in which every gene failed still commits nothing
      - `continue` is available only when no approved operation declares an
        external output path. GENtle checks this before the first operation and
        refuses the batch otherwise, because an in-memory gene rollback cannot
        revoke a file already written. Remove operation-level output paths or
        use the default `abort` policy for an all-or-nothing detached run
      - after a failed gene is rolled back, `last_checkpoint_manifest` points
        to the last checkpoint from a fully retained gene; it never names an
        operation checkpoint whose state was subsequently discarded
      - with `continue`, checkpoints are written only after a complete gene
        succeeds. This prevents the run from adding a discoverable checkpoint
        for a partial state from a later rolled-back gene. The default `abort`
        policy retains per-operation checkpoints for failure recovery
      - `continue` governs execution failures only. A batch-basis, plan,
        workflow, operation-digest, or endpoint-feasibility rejection still
        refuses the complete batch under either policy, because those are
        approval failures rather than results. The requested
        `--coverage-policy` is also unchanged: a gene that fails strict
        coverage is skipped, never downgraded
      - after a `continue` run with failures the live baseline has changed, so
        that run's frozen checkpoint no longer matches and its reuse proposal
        is correctly refused. Compose a new batch containing only the failed
        genes against the committed baseline; do not re-run the original batch
        with `--reuse-proposal`
      - inspect an unchanged approved prefix with
        `primers inspect-gene-isoform-study-reuse BATCH_JSON_OR_@FILE
        --checkpoint-dir DIR [--path PROPOSAL.json]`. A matching
        `gentle.gene_isoform_assay_study_reuse_proposal.v1` binds the baseline
        project, exact ordered prefix, checkpoint bytes, and running GENtle
        executable. Approve that exact `proposal_sha256`, then execute with
        `--reuse-proposal @PROPOSAL.json --approve-reuse-sha256 SHA256`.
        Rebuilding GENtle, changing an approved operation, or changing the
        baseline invalidates the proposal instead of silently recomputing or
        importing work

        ```json
        {
          "schema": "gentle.gene_isoform_assay_study_workflow_batch_request.v1",
          "batch_id": "tp73_association_second_stage",
          "entries": [
            {"plan_path": "analysis/PATZ1.plan.json", "workflow_path": "analysis/PATZ1.workflow.json"},
            {"plan_path": "analysis/TFAP2C.plan.json", "workflow_path": "analysis/TFAP2C.workflow.json"}
          ]
        }
        ```
      - policy decisions keep exact-cDNA complexity, annotation coverage,
        Clariom support/abundance/contrast responsiveness, assayability, prior
        priority, missing evidence, and explicit overrides as separate fields
      - `primers publish-gene-isoform-study REQUEST.json OUTPUT_DIR [--profile ID] [--blocks ID,ID] [--pdf]` renders only declared blocks from
        content-addressed plan/handoff/order reports. It emits a meta-page, one
        HTML page per gene, `print.html`, order TSVs, and a projection receipt
      - a dossier can be published while part of a study is still running or
        has failed. Each request gene may declare `status` (`resolved`,
        `pending`, `unresolved`) and a `status_reason`. An omitted status is
        derived: a gene with at least one assay handoff is `resolved`, a gene
        without one is `pending`, so an existing request keeps working and
        still reports honestly. `unresolved` requires a `status_reason` stating
        why the established automatism could not address the gene, and a gene
        declared `resolved` without any handoff is rejected rather than
        published as complete
      - the index lists every gene with its status and reason, the canonical
        report carries `resolved_gene_count`, `pending_gene_count`,
        `unresolved_gene_count`, and `complete`, and each incomplete gene page
        and the print document open with a visible notice. Biologists can start
        on the finished genes immediately; regenerate the dossier from the same
        request once the outstanding results are in
      - all standard profiles include a detailed quality-assurance block with
        the exact handoff policy and panel hash, every per-assay gate and
        required flag, GENtle's recorded status/summary, evidence-report ids,
        blockers, and warnings. Annotated-transcript product testing and
        whole-transcriptome/cDNA specificity are separate gates; missing,
        incomplete, and `not_evaluated` evidence is not presented as a pass
      - all standard profiles also include an explanatory provenance table for
        the bound study plan, experimental handoff, and oligo order form. It
        reports each source's role, schema/report id, path, SHA-256, and the
        relevant operation-batch, panel, policy, handoff, and duplicate-review
        bindings. Provenance describes processing and identity; it does not
        turn user-supplied/external evidence into a GENtle result
      - `primers oligo-order from-experimental-handoff HANDOFF.json --expected-sha256 SHA256 [...]` creates a form only from rows that passed the
        handoff's exact readiness policy and retains the handoff/policy/card/
        assay/pair/oligo provenance on every line

<!-- manual-page-nav:start -->

---

[CLI Manual](../cli.md) | [Previous: Reference-Bound Primer Specificity](primer-specificity.md) | [Next: Restriction Cloning](restriction-cloning.md)
<!-- manual-page-nav:end -->
