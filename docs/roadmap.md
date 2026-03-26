# GENtle Roadmap and Status

Last updated: 2026-03-26

Purpose: shared implementation status, known gaps, and prioritized execution
order. Durable architecture constraints and decisions remain in
`docs/architecture.md`. Machine contracts remain in `docs/protocol.md`.

## 1. Current implementation snapshot

### Engine and adapter baseline in place

- Shared engine core in `src/engine.rs` with operation/workflow execution.
- Engine module decomposition is now underway:
  - `src/engine/tests.rs` (engine test suite extracted from monolithic engine file)
  - `src/engine/ops/candidate_guides.rs` (candidate-set, guide-design, and macro-template operations)
  - `src/engine/analysis/candidate_metrics.rs` (candidate metrics, feature-distance helpers, and expression scoring)
  - `src/engine/analysis/feature_expert_ops.rs` (feature-expert, splicing, and isoform helper routines)
  - `src/engine/io/genome_tracks.rs` (genome track parsing/import helpers for BED/BigWig/VCF/BLAST overlays)
  - `src/engine/io/import_anchors.rs` (GenBank-derived anchor parsing and import-origin helpers)
  - `src/engine/state/lineage_containers.rs` (lineage-node and container/arrangement helper routines)
  - `src/engine/ops/operation_handlers.rs` (core operation-dispatch handler extracted from `apply_internal`)
  - `src/engine/analysis/rna_reads.rs` (RNA-read report parsing, scoring, and export helper routines)
  - `src/engine/state/sequence_ops.rs` (sequence digestion, enzyme resolution, FASTA/pool export, and primer/overhang utilities)
- Shared shell layer in `src/engine_shell.rs` reused by GUI Shell and
  `gentle_cli shell`.
- Engine-shell decomposition is now underway:
  - `src/engine_shell/command_parsers.rs` (candidate/guide/macro/routine/planning command parser helpers)
  - `src/engine_shell/tests.rs` (engine-shell test suite extracted from monolithic shell file)
- CLI adapter in `src/bin/gentle_cli.rs` with state/capability utilities and
  first-class command trees (`genomes`, `helpers`, `resources`, `tracks`,
  `ladders`, `candidates`, `import-pool`).
- Planning meta-layer shell/CLI baseline is now available:
  - direct/shared-shell command family: `planning ...`
  - profile/objective/suggestion/sync contracts:
    `gentle.planning_profile.v1`,
    `gentle.planning_objective.v1`,
    `gentle.planning_estimate.v1`,
    `gentle.planning_suggestion.v1`,
    `gentle.planning_sync_status.v1`
  - deterministic profile merge precedence:
    `global -> confirmed_agent_overlay -> project_override`
  - advisory agent sync lifecycle with explicit user
    `accept`/`reject` activation policy
  - strict schema compatibility checks for planning profile/objective payloads
    (mismatched schema ids fail fast)
  - procurement/queue business-day delays are converted to elapsed
    `estimated_time_hours` with deterministic weekend-aware mapping
    (`24h * 7/5` per business day)
  - GUI planning baseline is now available:
    - standalone `Planning` window in `Patterns -> Planning...`
    - command palette action: `Planning`
    - in-window profile/objective editing + pull/push suggestion registration
      + explicit accept/reject resolution
- Gibson specialist baseline is now available:
  - standalone `Gibson` window in `Patterns -> Gibson...`
  - command palette action: `Gibson`
  - shared non-mutating preview path:
    `gibson preview PLAN_JSON_OR_@FILE [--output OUTPUT.json]`
  - shared mutating apply path:
    `gibson apply PLAN_JSON_OR_@FILE`
  - current scope:
    - destination-first Gibson preview/apply for one or more ordered inserts
    - resolved junction chain (`n + 1` junctions for `n` inserts)
    - Gibson-specific primer suggestions (`5' overlap + 3' priming`) with two
      primers per insert fragment
    - validation/advisory findings
    - destination opening suggestions start with unique cutters named by the
      MCS annotation, with an explicit option to reveal other unique cutters
      on the destination
    - specific cutter suggestions now show REBASE-derived recognition/cut
      information and fill the actual cleavage window rather than the whole
      recognition span
    - shared Tₘ-model assumptions surfaced in preview notes and Gibson-window
      help so GUI/CLI users see the same explanation
    - dedicated `Tₘ Model` box in the Gibson window
    - apply-time sequence creation for all insert primers plus the assembled
      product
    - Gibson apply outputs now materialize as separate singleton containers
      (one vial per primer plus the assembled product) instead of one
      synthetic pool container
    - protocol-cartoon preview now shows both destination-insert junctions
      explicitly for the current single-insert flow instead of collapsing to
      one representative overlap
    - multi-insert preview now carries an ordered dynamic cartoon spec through
      the same export/review path, including strand-specific 5' chew-back,
      exposed overlap/tail geometry, and downstream fill-in/seal states
    - deterministic feature transfer onto the assembled Gibson product
    - lineage reopen path: clicking a Gibson operation reopens the specialist
      with the saved plan loaded again
    - partially consumed destination annotations are now trimmed/projected when
      the product-side rewrite is honest enough to preserve them
    - surviving MCS annotations are cross-checked against actual assembled
      product restriction sites and uniqueness, including insert-introduced
      sites
    - protocol-cartoon preview/export
    - Routine Assistant handoff when the destination already uses existing
      linear termini
  - current limitation:
    - multi-insert execution currently requires a defined destination opening;
      `existing_termini` remains the single-fragment handoff path
- CLI state loading now treats empty/whitespace `--state` files as
  uninitialized and starts from `ProjectState::default()` instead of failing
  parse.
- Primer-design CI split is now explicit:
  - always-run workflow/tutorial examples pin `primer_design_backend=internal`
    so routine CI remains deterministic and does not depend on a local
    `primer3_core`
  - real external-binary Primer3 smoke tests are opt-in via
    `GENTLE_TEST_EXTERNAL_BINARIES=1`
- Remote-resource test execution now supports explicit skip policy via
  `GENTLE_SKIP_REMOTE_TESTS=1` (overrides `GENTLE_TEST_ONLINE=1`) for
  packaging/offline environments.
- JS and Lua adapters expose shared operation/workflow bridges and convenience
  wrappers over engine contracts.
- Embedded JS/Lua adapters are now compile-time optional build targets:
  - default Cargo builds avoid `deno_core`/`mlua` unless scripting features are
    explicitly enabled
  - Cargo features:
    - `js-interface`
    - `lua-interface`
    - `script-interfaces`
  - `gentle_js` / `gentle_lua` binaries require their matching feature
  - release packaging builds enable `script-interfaces` so tagged release
    artifacts are compiled against the broader embedded scripting surface
- Python adapter baseline is now available as a thin `gentle_cli` wrapper:
  - path: `integrations/python/gentle_py/`
  - deterministic methods for `capabilities`, `state-summary`, `op`,
  `workflow`, and `shell`
  - no biology logic duplication (subprocess bridge only)
- Debian-first container baseline is now available:
  - authoritative image definition in repository `Dockerfile`
  - GUI served through `Xvfb` + `openbox` + `x11vnc` + `noVNC`
  - CLI/MCP/JS/Lua/Python wrapper available from the same image
  - helper-tool coverage in container:
    - Debian packages for `primer3`, `ncbi-blast+`, `python3-pybigwig`,
      `novnc`, `x11vnc`, `xvfb`
    - explicit container compatibility wrapper for `bigWigToBedGraph`
    - explicit non-Debian exception for `rnapkin`
  - user documentation lives in `docs/container.md`
  - GHCR workflow now lives in `.github/workflows/container.yml`

### Biology/analysis capabilities already implemented

- Container-first model in engine state (`ProjectState.container_state`) with
  container-aware digest/merge/ligation/filter operations.
- Reference-genome preparation and extraction (`PrepareGenome`,
  `ExtractGenomeRegion`, `ExtractGenomeGene`, `ExtendGenomeAnchor`) including
  catalog-backed helper flows and BLAST integration.
- Genome-menu specialist dialogs now preserve independent Reference vs Helper
  setup state (`catalog`/`cache_dir`) and use scope-specific window titles +
  BLAST target messaging so helper actions no longer inherit reference
  interface/context labels.
- `ExtractGenomeRegion` now supports explicit annotation projection scopes
  (`annotation_scope=none|core|full`), defaulting to `core` when unset so
  extracted coordinate slices include overlapping gene/transcript features by
  default.
- `ExtractGenomeRegion` now emits structured annotation projection telemetry in
  `OpResult` (`genome_annotation_projection`) with requested/effective scope,
  attached/dropped counts, and fallback reason fields.
- `ExtractGenomeRegion` now supports `max_annotation_features` caps with
  deterministic fallback (`full -> core -> none`) and warning guidance for
  follow-up unrestricted transfer.
- `ExtractGenomeGene` now supports interval selection via
  `extract_mode=gene|coding_with_promoter` plus `promoter_upstream_bp`, so
  GUI/CLI/scripting callers can retrieve CDS-spanning genomic context with an
  explicit strand-aware 5' promoter flank instead of manually calculating
  coordinates.
- Shell/CLI `genomes/helpers extract-region` now expose
  `--annotation-scope`/`--max-annotation-features` while retaining legacy
  `--include-genomic-annotation` / `--no-include-genomic-annotation` flags for
  compatibility.
- GUI `Retrieve Genomic Sequence` now exposes matching extraction controls for
  annotation scope/cap (`none|core|full`, optional `max_annotation_features`)
  and surfaces structured projection outcomes in status output (requested/effective
  scope, attached/dropped counts, optional fallback reason).
- Deterministic tests now cover:
  - default core projection telemetry,
  - full-scope cap fallback behavior,
  - shell parse/execute paths for new scope/cap flags.
- Imported GenBank region anchors now record verification status against catalog
  sequence (`verified`/`unverified`/`verification n/a`) and extension now
  hardens against stale/invalid anchor `catalog_path` by falling back to the
  default genome catalog when no explicit catalog override is provided.
- Helper genome extraction now auto-attaches a canonical `misc_feature` MCS
  annotation for pUC18/pUC19 when source annotation does not already provide an
  MCS feature (motif-detected, deterministic qualifiers, uniqueness required).
- Engine Ops digest panel now includes MCS/REBASE quick-fill actions:
  `Use MCS enzymes`, `Single-cutters (REBASE)`, and
  `Single-cutters in CDS`.
- Next UI polish target for the DNA sequence window:
  make MCS annotations more prominent and directly legible in the sequence
  view, so cloning-oriented users can recognize the intended insertion/cut
  region without needing to expose `misc_feature` detail manually first.
- New strict policy switch is available:
  `SetParameter(require_verified_genome_anchor_for_extension=true)` enforces
  verified anchors for `ExtendGenomeAnchor` (unverified/unknown anchors fail).
- Genome-anchor extension/read paths now resolve assembly-family compatible
  prepared caches (for example `GRCh38.p14` -> prepared `Human GRCh38 ...`)
  when the match is unique; ambiguous compatible matches return explicit
  options instead of guessing.
- Genome-anchor resolution policy now supports deterministic modes via
  `SetParameter(genome_anchor_prepared_fallback_policy=off|single_compatible|always_explicit)`.
- Sequence-window anchored controls now include:
  - explicit prepared-genome chooser popup when multiple compatible prepared
    references exist,
  - `Re-verify anchor` action (engine `VerifyGenomeAnchor`),
  - visible `anchor check: verified|unverified|n/a` badge.
- Sequence-window open path is now lazy:
  - sequence payload is loaded asynchronously on window open,
  - feature-tree/details rendering can be deferred per window and enabled
    interactively (`Load feature tree now`) to keep first paint responsive on
    very feature-dense entries.
- `PrepareGenome` now detects source-path/URL drift for existing manifests and
  rebuilds from configured sources instead of silently reusing stale caches.
- Prepared-genome reindexing now has an explicit cached-files-first path:
  reindex keeps cached local sequence/annotation files, progress appears
  immediately when launched, and deleting cached downloads for a full refresh
  requires an explicit confirmation choice.
- Genome-catalog maintenance now has a shared baseline:
  - bundled `assets/genomes.json` now includes additional Ensembl templates for
    zebrafish, chimp, dog, and Drosophila, while fixing the rat pinned release
  - catalog rows can carry `ensembl_template` metadata so GENtle can refresh
    explicit pinned Ensembl URLs on demand without preparing genomes
  - `Prepared References...` now separates `Remove Prepared...` from
    `Remove Catalog Entry...` so cache deletion and catalog editing are explicit
- Prepare/reindex/refresh dialogs now expose a structured step checklist
  up front (per-step bars, completed checkmarks, compact jobs-panel summary)
  instead of a single raw phase bar, so long-running indexing work stays
  visible and auditable.
- Genome track import operations (BED, BigWig via conversion, VCF, BLAST hits)
  with anchor-aware coordinate remapping.
- Resource ingestion/update path for REBASE and JASPAR snapshots across GUI/CLI
  and scripting adapters.
- TFBS annotation guardrails (default cap, explicit unlimited mode), progress
  reporting, and persistent display-time filtering criteria.
- Feature expert-view pipeline for selected TFBS, restriction sites, and
  splicing groups, plus curated isoform architecture panels:
  - shared expert payload generation in engine
  - GUI-side expert panel rendering from engine payload
  - SVG export via `RenderFeatureExpertSvg` and
    `RenderIsoformArchitectureSvg`
  - shell/CLI/JS/Lua access:
    - `inspect-feature-expert` / `render-feature-expert-svg`
    - `panels import-isoform` / `panels inspect-isoform` /
      `panels render-isoform-svg`
- Transcript derivation baseline is now additive and parity-wired:
  - engine operation `DeriveTranscriptSequences` derives cDNA/transcript
    sequences from `mRNA`/`transcript` features
  - shared-shell/CLI command:
    `transcripts derive SEQ_ID [--feature-id N ...] [--scope ...] [--output-prefix PREFIX]`
  - Splicing Expert quick actions in GUI:
    - `Derive group transcripts`
    - `Derive all mRNA`
    - `Derive + Dotplot` (selected transcript; strand-aware pair mode)
- New sequence-analysis operation baseline:
  - `DeriveSplicingReferences` derives a sequence-window DNA slice, transcript-
    oriented mRNA isoforms, and an exon-consecutive artificial reference in one
    deterministic engine operation.
  - `AlignSequences` adds deterministic local/global pairwise alignment
    reporting (score/identity/coverage/CIGAR-like ops) as structured
    `OpResult.sequence_alignment` payloads.
  - shared-shell/direct command routing is available for both operations:
    - `splicing-refs derive ...`
    - `align compute ...`
- VCF display filtering parity between GUI and SVG export (`SetParameter`/shared
  display criteria).
- Candidate-set workflow (generate/score/filter/set operations + macro scripts)
  with persistent storage model, including between-anchor window generation.
- Workflow and candidate macro template catalogs now preserve optional external
  reference URLs (`details_url`) alongside template names/descriptions across
  engine and shared shell adapter surfaces.
- Workflow macro pattern import now accepts hierarchical file catalogs:
  - directory-recursive import via `macros template-import PATH`
  - per-template schema support (`gentle.cloning_pattern_template.v1`)
  - starter hierarchy at `assets/cloning_patterns_catalog`
  - GUI `Patterns` menu mirrors folder hierarchy for in-app selection/import
  - `template-import PATH` accepts single template file, legacy pack file, or
    folder tree
  - template execution can be repeated per project with different bindings and
    persists resulting sequences/containers via regular operation persistence
- Typed cloning-routine catalog baseline is now available:
  - versioned manifest at `assets/cloning_routines.json`
  - additive explainability metadata baseline is now supported in catalog rows:
    - `purpose`, `mechanism`, `requires`, `contraindications`
    - `confusing_alternatives`, `difference_matrix`
    - `disambiguation_questions`, `failure_modes`
  - shared-shell/CLI discovery command:
    `routines list [--catalog PATH] [--family NAME] [--status NAME] [--tag TAG] [--query TEXT]`
    - routine rows now include planning estimate fields when planning context
      exists (`estimated_time_hours`, `estimated_cost`, `local_fit_score`,
      `composite_meta_score`, full `planning_estimate` explanation payload)
  - shared-shell/CLI explainability commands:
    - `routines explain ROUTINE_ID [--catalog PATH]`
    - `routines compare ROUTINE_A ROUTINE_B [--catalog PATH]`
    - compare matrix now includes planning-estimate axes
      (`estimated_time_hours`, `estimated_cost`, `local_fit_score`,
      `composite_meta_score`)
  - GUI staged Routine Assistant baseline is now available
    (`Patterns -> Routine Assistant...`) with:
    - goal/candidate selection
    - compare stage (`routines compare`)
    - typed binding form
    - preflight stage (`macros template-run --validate-only`)
    - transactional run + run-bundle export
    - Gibson topology guard:
      - circular sequence bindings block preflight/execute
      - one-click `Linearize Vector...` action branches + linearizes + rebinds
  - GUI `Patterns` menu now exposes routine browsing grouped by family/status
    and routine-level template import actions
  - workflow template preflight baseline is now wired into
    `macros template-run`:
    - typed routine input-port checks for supported kinds
      (`sequence`, `container`, `candidate_set`, `guide_set`, `string`,
      `number`, `bool`, `path`, `sequence_anchor`)
    - non-mutating inspection mode:
      `macros template-run ... --validate-only`
    - preflight payload is included in template-run output
  - workflow macro templates now support optional typed
    `input_ports`/`output_ports` metadata:
    - persisted in engine template state
    - import-compatible via `gentle.cloning_pattern_template.v1` files
    - shell authoring path via `macros template-put --input-port ... --output-port ...`
  - preflight now validates both input and output ports and reports
    `contract_source` (`template_ports` or `routine_catalog`) in
    `gentle.macro_template_preflight.v1`
  - preflight now also enforces baseline cross-port semantics:
    - output alias/collision diagnostics
    - sequence/container compatibility diagnostics
    - sequence-anchor semantic checks when one sequence context is bound
  - routine-family preflight baseline now includes:
    - Gibson overlap checks:
      - validates adjacent fragment suffix/prefix homology against configured
        overlap length
      - emits explicit mismatch diagnostics before mutation
    - Restriction digest/ligation checks:
      - validates enzyme-name resolution and duplicate-enzyme misuse
      - validates recognition-site presence across bound sequence inputs
      - validates digest parameter sanity (`left_fragment`/`right_fragment`,
        `extract_from`/`extract_to`)
  - first shipped Gibson + restriction routine/template baselines:
    - routine id: `gibson.two_fragment_overlap_preview`
    - template:
      `assets/cloning_patterns_catalog/gibson/overlap_assembly/gibson_two_fragment_overlap_preview.json`
    - routine id: `restriction.digest_ligate_extract_sticky`
    - template:
      `assets/cloning_patterns_catalog/restriction/digest_ligation/digest_ligate_extract_sticky.json`
  - mutating `macros run` and `macros template-run` now append deterministic
    lineage macro-instance records for success and failure/cancel pathways
    (typed bound inputs/outputs + emitted op ids + status/status_message)
  - macro-instance introspection commands are available:
    `macros instance-list`, `macros instance-show`
- Ladder-aware virtual gel rendering and SVG export routes, including
  container-based and arrangement-based serial gel export surfaces.
- Protocol-cartoon SVG generation baseline is now available:
  - engine operation `RenderProtocolCartoonSvg { protocol, path }`
  - engine operation `RenderProtocolCartoonTemplateSvg { template_path, path }`
  - engine operation `ValidateProtocolCartoonTemplate { template_path }`
  - top-level README now embeds a deterministic built-in Gibson strip rendered
    from the shared protocol-cartoon engine route
  - engine operation
    `RenderProtocolCartoonTemplateWithBindingsSvg { template_path, bindings_path, path }`
  - engine operation `ExportProtocolCartoonTemplateJson { protocol, path }`
  - shared-shell/CLI routes:
    - `protocol-cartoon list`
    - `protocol-cartoon render-svg PROTOCOL_ID OUTPUT.svg`
    - `protocol-cartoon render-template-svg TEMPLATE.json OUTPUT.svg`
    - `protocol-cartoon template-validate TEMPLATE.json`
    - `protocol-cartoon render-with-bindings TEMPLATE.json BINDINGS.json OUTPUT.svg`
    - `protocol-cartoon template-export PROTOCOL_ID OUTPUT.json`
    - command surface intentionally stays canonical under
      `protocol-cartoon ...` without extra aliases
  - shipped built-in protocol ids now include:
    - `gibson.two_fragment`
    - `gibson.single_insert_dual_junction`
    - `pcr.assay.pair`
    - `pcr.assay.pair.no_product`
    - `pcr.assay.pair.with_tail`
    - `pcr.oe.substitution`
    - `pcr.assay.qpcr`
  - renderer contract is now abstraction-first:
    - protocol strip = ordered event sequence
    - event payload = one or more DNA molecules
    - molecule payload supports:
      - topology: `linear|circular`
      - feature fragments with explicit lengths/colors
      - optional split strand coloring and per-strand lengths per feature
        (top/bottom independently)
      - linear-end styles:
        `NotShown`, `Continuation`, `Blunt`, or
        `Sticky { polarity: FivePrime|ThreePrime, nt }`
    - template schema baseline:
      - `gentle.protocol_cartoon_template.v1`
      - sparse template inputs are expanded with deterministic defaults
        (action/caption/topology/end styles/feature lengths/palette)
    - bindings schema baseline:
      - `gentle.protocol_cartoon_template_bindings.v1`
      - deterministic overrides target defaults and event/molecule/feature ids
    - validation rejects malformed cartoons (empty ids/events, zero-length
      features, circular molecules with linear end styles, zero-nt sticky ends)
  - Gibson specialist experience now confirms the intended growth model:
    extend protocol cartoons by template families + bindings (including
    repeated events/molecules for multi-fragment cases), with shared figure
    building blocks under the built-ins, not protocol-specific renderer
    branches
  - internal shared protocol-cartoon building-block baseline is now in place
    for:
    - duplex feature spans
    - strand-specific single-strand tails
    - linear molecule rows with deterministic end styles
    - event rows reused across Gibson variants and intended for future PCR
      variants
  - PCR-assay family baseline is now shipped through `pcr.assay.*`:
    - `pcr.assay.pair`
    - `pcr.assay.pair.no_product`
    - `pcr.assay.pair.with_tail`
    - `pcr.assay.qpcr`
    - pair/qPCR renders now show explicit primer glyphs with oriented 5'/3'
      end labels
    - tailed pair-PCR renders now show non-annealing 5' primer tails and
      their deterministic carry-through into the product lane
  - overlap-extension substitution mutagenesis baseline is now shipped through
    `pcr.oe.substitution`:
    - six-panel OE-PCR strip with primer set `a..f`
    - strand-specific anneal-gap lane plus explicit polymerase fill lane
    - deterministic geometry bindings (`flank_bp`, `overlap_bp`, `insert_bp`)
      can be applied through template bindings
  - next protocol-cartoon family expansion for PCR:
    - batch, nested, inverse, and richer artifact/readout variants through the
      same template/binding abstraction
- Top-level README showcase expansion is now underway:
  - shipped showcase figures now cover distinct GENtle usage modes:
    - built-in Gibson protocol cartoon
    - single-insert dual-junction Gibson mechanism strip
    - Gibson lineage graph export showing one operation plus concrete outputs
    - TP53/p53 isoform architecture from Ensembl 116 geometry plus curated
      isoform panel overlays
    - TP73 cDNA-vs-genomic dotplot derived locally from `test_files/tp73.ncbi.gb`
      and rasterized to README PNG via `gentle_examples_docs svg-png`
      (`resvg`) while preserving bp axis labels
  - current README dotplot path is intentionally offline:
    - derive `NM_001126241.3` from local TP73 GenBank feature geometry
    - compute dotplot through shared `DeriveTranscriptSequences` +
      `ComputeDotplot` + `RenderDotplotSvg` engine routes
    - convert SVG -> PNG in docs tooling to keep README rendering lightweight
  - next showcase candidates queued on the same track:
    - end-to-end cloning assembly
      (`docs/examples/workflows/digest_ligation_extract_region_minimal.json`)
    - TP53-family sparse RNA mapping
      (`docs/examples/workflows/tp53_multi_gene_sparse_mapping_online.json`)
    - batch primer design
      (`docs/examples/workflows/pcr_selection_batch_primer_pairs_offline.json`)
    - overlap-extension substitution mutagenesis baseline
      (`docs/examples/workflows/pcr_overlap_extension_substitution_offline.json`)
    - guide export and protocol handoff
      (`docs/examples/workflows/guides_export_csv_and_protocol.json`)
    - TP73 promoter luciferase assay planning
      (`docs/examples/workflows/tp73_promoter_luciferase_assay_planning.json`)
- Deterministic process run-bundle export baseline is now implemented:
  - engine operation `ExportProcessRunBundle { path, run_id? }`
  - shared-shell/CLI command
    `export-run-bundle OUTPUT.run_bundle.json [--run-id RUN_ID]`
  - exported schema `gentle.process_run_bundle.v1` includes:
    - extracted operation inputs
    - chronological `SetParameter` overrides
    - selected immutable operation log records
    - output summaries (created/changed entities + exported artifacts)
    - full parameter snapshot
- Primer-design report baseline:
  - engine operation `DesignPrimerPairs` now persists deterministic report
    payloads (`gentle.primer_design_report.v1`) in project metadata
  - engine operation `DesignInsertionPrimerPairs` now provides an
    insertion-first wrapper over pair design (requested anchors + extension
    sequences + anchor-adjacent windows) while persisting the same report
    schema with optional `insertion_context` compensation rows
  - engine operation `PcrOverlapExtensionMutagenesis` now supports
    overlap-extension insertion/deletion/replacement workflows with
    deterministic inner-overlap tail synthesis and graph-visible staged
    artifacts (outer/inner primers, stage-1 left/right products, stage-2 final
    mutant product, and stage-specific containers)
    - insertion/replacement runs now also emit
      `OpResult.protocol_cartoon_preview` metadata for
      `pcr.oe.substitution` (geometry + deterministic template bindings)
    - GUI primer panel now surfaces the last captured
      `OpResult.protocol_cartoon_preview` payload, including one-click
      render/export of the bound protocol-cartoon SVG for operation-level
      inspection/handoff
  - engine operation `DesignQpcrAssays` now persists deterministic qPCR report
    payloads (`gentle.qpcr_design_report.v1`) with forward/reverse/probe assay
    records in project metadata
  - side-sequence constraints are now supported for forward/reverse/probe:
    - `non_annealing_5prime_tail` (5' add-on excluded from anneal Tm/GC/hit)
    - `fixed_5prime`, `fixed_3prime`
    - `required_motifs[]`, `forbidden_motifs[]`
    - `locked_positions[]`
  - pair-level amplicon constraints are now supported through
    `pair_constraints`:
    - `require_roi_flanking`
    - required/forbidden amplicon motifs
    - fixed amplicon start/end coordinates
  - shared-shell/CLI inspection/export commands are available:
    `primers design`, `primers list-reports`, `primers show-report`,
    `primers export-report`, `primers design-qpcr`,
    `primers list-qpcr-reports`, `primers show-qpcr-report`,
    `primers export-qpcr-report`
  - GUI Engine Ops now exposes dedicated primer/qPCR forms for those operations,
    including explicit side-sequence constraints and pair constraints (no raw
    JSON required for common interactive use)
  - GUI Engine Ops now also exposes report-management helpers in the same
    primer/qPCR panel:
    - list persisted report IDs
    - show report summary for current `report_id`
    - export current `report_id` to JSON via save dialog
  - DNA-window `PCR ROI` menu now supports selection-first queue capture:
    - add current selection to PCR queue
    - seed/queue from current visible linear map span (selection-independent fallback)
    - add selected feature(s) to PCR queue (one region per feature)
  - DNA-window toolbar now also exposes a visible `Queue PCR selection`
    action beside `Extract Sel` so the current linear selection can be queued
    without opening the menu first; the map area also shows an inline
    selection-ready hint when a non-empty span exists
  - paint-first PCR coordinate UX baseline is now in place (pair-PCR v1):
    - fixed semantic paint roles on linear DNA map:
      ROI (green), upstream primer window (red), downstream primer window (blue)
    - drag paints one interval for the selected role and shows live
      `start..end (len bp)` overlay text
    - `Shift+drag` on ROI immediately appends queue row metadata; `Option/Alt+drag`
      preserves pan behavior with higher priority than paint
    - post-drag action chip now exposes:
      `Set PCR ROI`, `Add ROI to Queue`, `Open PCR Designer`
    - pair-PCR live geometry preview in the primer panel now resolves through
      `pcr_assay_pair_geometry_bindings(...)` with deterministic clamping to
      baseline template maxima
  - GUI primer panel now includes queued PCR batch execution:
    - queue table (`source`, `template`, `start/end/len`) with row remove/clear
    - explicit `Queue current ROI spec` action in-panel, clarifying that queue
      rows are ROI specifications (`template + start + end`) rather than
      immediate Primer3 jobs
    - explicit `Queue visible map span` in-panel shortcut for map-first workflows
    - batch run action `Design Primer Pairs for queued regions`
    - deterministic batch report suffixing (`{base}_rNN`)
    - optional per-region `ExtractRegion` copy artifacts
    - batch-results table with per-region status and quick actions:
      `Show` / `Export` / `Open` (copy-first fallback to template)
  - GUI primer/qPCR forms now include field-level hover help for ROI,
    side/pair constraints, and motif-format expectations (IUPAC + comma-separated
    examples), including `require ROI flanking`
  - dedicated PCR Designer specialist window is now available:
    - `Patterns -> PCR Designer...` and command palette `PCR Designer`
    - sequence-context aware dedicated viewport with paint controls + map +
      queue summary on left and pair-PCR constraints/run/report panel on right
    - qPCR remains in the existing Engine Ops panel for this v1 scope
    - shared-shell UI intents now include:
      `ui open pcr-design` and `ui focus pcr-design`
  - single-run GUI `Design Primer Pairs` execution now runs asynchronously in a
    background worker to keep the sequence window responsive while Primer3 is
    running
  - queued GUI batch `Design Primer Pairs for queued regions` now runs through
    the same asynchronous background worker path (no UI-thread blocking while
    Primer3 executes across queued regions)
  - internal primer-pair candidate cross-product evaluation now uses a
    deterministic bounded budget (scaled from requested `max_pairs` with a hard
    cap) so pathological ROI/constraint combinations cannot appear
    non-terminating; skipped combinations are surfaced in rejection summaries
    and operation warnings (PCR + qPCR paths)
  - backend selection is now available through engine parameters and shell
    command options:
    - `primer_design_backend=auto|internal|primer3`
    - `primer3_executable` path override
    - `primers design ... --backend ... --primer3-exec ...`
  - report metadata now records backend provenance
    (`requested`, `used`, optional fallback reason + Primer3 executable/version)
  - Primer3-backed reports now persist explain diagnostics and the exact
    Boulder-IO request payload used for that run, enabling deterministic local
    reruns of failed/empty designs from the exported request file
  - GUI primer report summaries now explicitly flag zero-pair outcomes as
    `NO ACCEPTED PRIMER PAIRS` with rejection counters + Primer3 explain text,
    and the primer panel includes `Export Primer3 input...` for report-level
    request export
  - `auto` mode now falls back deterministically to internal scoring when
    Primer3 is unavailable
  - `DesignPrimerPairs` is now lineage/graph-inspectable as a mutating
    operation:
    - each accepted pair materializes three derived sequences:
      - `..._fwd`
      - `..._rev`
      - `..._amplicon` (predicted product, including configured 5' tails)
    - one dedicated container is created per primer pair
      (forward + reverse + predicted amplicon)
    - lineage edges are emitted from the template to each created primer
      or amplicon sequence with the operation id
    - generic aggregate container auto-creation is skipped for this operation
      so `seq_to_latest_container` stays pair-scoped
  - pair-ranking/scoring now includes explicit primer-quality heuristics:
    - preferred length window (`20..30 bp`)
    - 3' GC-clamp preference
    - secondary-structure penalties (homopolymer/self-complementary runs)
    - primer-pair dimer-risk penalties (global and 3'-anchored)
  - report payloads now include per-primer and pair diagnostics for those
    heuristics (clamp status, run lengths, dimer run metrics)
  - deterministic engine tests now cover side/pair/probe constraint pass/fail
    behavior and shared-shell primer/qPCR request parsing with the extended
    operation payload
- UniProt/SWISS-PROT import + genome-projection baseline:
  - engine operations:
    - `ImportUniprotSwissProt` (offline SWISS text ingestion)
    - `FetchUniprotSwissProt` (online accession/id fetch)
    - `ProjectUniprotToGenome` (feature-to-genome mapping via transcript/CDS)
    - `FetchGenBankAccession` (online GenBank accession fetch + direct import)
    - `ImportUniprotEntrySequence` remains present for protocol compatibility
      but is intentionally disabled (`Unsupported`) until first-class
      protein-sequence windows are implemented.
  - persisted metadata stores:
    - `gentle.uniprot_entries.v1`
    - `gentle.uniprot_genome_projections.v1`
  - shared-shell commands:
    - `genbank fetch ACCESSION [--as-id ID]`
    - `uniprot fetch`, `uniprot import-swissprot`, `uniprot list`,
      `uniprot show`, `uniprot map`, `uniprot projection-list`,
      `uniprot projection-show`
  - next UI follow-up:
    - add a dedicated projection viewer for protein -> transcript -> genome
      mapping results (for example TP53/p53), reusing the persisted
      `gentle.uniprot_genome_projections.v1` payload and keeping the GUI as a
      thin presentation layer over shared projection state rather than adding
      first-class protein sequence windows
- Executable tutorial baseline is now integrated with canonical workflow
  examples:
  - canonical tutorial landing page now exists at `docs/tutorial/README.md`
    and links:
    - generated executable chapter hub
    - hand-written GUI walkthroughs
    - agent/reference tutorial material
  - canonical machine-readable discovery catalog now exists at:
    `docs/tutorial/catalog.json` (`gentle.tutorial_catalog.v1`)
    - current role:
      - unified discovery metadata for generated and hand-written tutorial pages
      - consumed by GUI help/tutorial discovery for curated ordering and
        inclusion of agent/reference pages
      - used as the root for resolving executable tutorial manifest location in
        `Open Tutorial Project...`
    - intended future role: broader shared source for tutorial discovery
      surfaces beyond the current GUI/help path
  - tutorial discovery authoring now has an explicit source layer:
    - `docs/tutorial/sources/catalog_meta.json`
    - `docs/tutorial/sources/*.json`
    - generation/check commands:
      - `gentle_examples_docs tutorial-catalog-generate`
      - `gentle_examples_docs tutorial-catalog-check`
      - `gentle_examples_docs tutorial-manifest-generate`
      - `gentle_examples_docs tutorial-manifest-check`
  - tutorial pages now distinguish their role more explicitly:
    `generated+checked`, `manual`, `manual/hybrid`, `manual/reference`
  - generated tutorial runtime manifest:
    `docs/tutorial/manifest.json` (`gentle.tutorial_manifest.v1`)
  - tutorial chapters now include explicit narrative learning objectives and
    concept tags, with generated recurrence mapping ("where this concept
    reappears") in chapter pages and tutorial index
  - committed generated tutorial output:
    `docs/tutorial/generated/` (chapters + retained artifacts + report)
  - dedicated contributor onboarding tutorial chapter added:
    `contribute_to_gentle_development` (backed by executable
    `contribute_gentle_development_baseline` workflow example)
  - dedicated manual Gibson specialist testing tutorial added:
    `docs/tutorial/gibson_specialist_testing_gui.md`
    - local committed inputs only
    - explicit GUI checkpoints for preview, exports, and cartoon visibility
    - explicit shared CLI parity check through `gibson preview`
  - executable `Open Tutorial Project...` baseline added for Gibson specialist
    testing:
    - workflow example:
      `docs/examples/workflows/gibson_specialist_testing_baseline.json`
    - generated chapter:
      `docs/tutorial/generated/chapters/15_gibson_specialist_testing_baseline.md`
    - role:
      - preloads stable destination/insert IDs for the manual Gibson testing
        walkthrough
      - now opens the matching Help/Tutorial guide automatically when the
        starter project is opened from the GUI
  - dedicated online TP53 isoform architecture chapter added:
    `tp53_isoform_architecture_online` (backed by canonical
    `tp53_isoform_architecture_online` workflow example and curated panel
    resource `assets/panels/tp53_isoforms_v1.json`)
  - dedicated online TP63 anchor-extension chapter added:
    `tp63_anchor_extension_online` (backed by canonical
    `tp63_extend_anchor_online` workflow example for coordinate retrieval +
    +/-2 kb extension)
  - dedicated offline PCR selection-first chapter added:
    `pcr_selection_batch_primer_pairs_offline` (backed by canonical
    `pcr_selection_batch_primer_pairs_offline` workflow example with
    deterministic multi-region primer report IDs)
    - chapter walkthrough now includes explicit mixed-source ROI queue capture
      checklist (current selection + selected features), plus row-level
      batch-results action guidance (`Show`/`Export`/`Open`)
    - TP73 promoter GUI tutorial Step 5 now references the dedicated PCR chapter
      and mirrors the same queue-first flow
  - dedicated Gibson planning chapter added:
    `gibson_two_fragment_overlap_preview` (backed by canonical
    `gibson_two_fragment_overlap_preview` workflow example and routine-catalog
    template import/run path)
  - GUI now exposes `File -> Open Tutorial Project...`:
    chapter entries are loaded from tutorial manifest/examples, grouped by tier
    (`Core`/`Advanced`/`Online`), materialized through shared workflow engine
    execution, and opened as inspectable project states
  - `gentle_examples_docs` now supports:
    `tutorial-generate` and `tutorial-check`
  - CI now gates tutorial drift and workflow runtime coverage through
    `tutorial-check` and `cargo test workflow_examples`
  - CI now includes a deterministic internal-preview smoke route that exercises:
    - workflow examples (`load_branch_reverse_complement_pgex_fasta`,
      `digest_ligation_extract_region_minimal`)
    - fast primer-pair report generation + persisted report inspection
  - CI now also runs the full Rust test suite (`cargo test -q`) to catch
    cross-module regressions that may not be covered by targeted smoke steps.
- ClawBio/OpenClaw external skill scaffold is now available for deterministic
  GENtle execution from a skill runtime:
  - path: `integrations/clawbio/skills/gentle-cloning/`
  - wrapper request/result schemas:
    `gentle.clawbio_skill_request.v1` /
    `gentle.clawbio_skill_result.v1`
  - outputs include reproducibility bundle artifacts (`commands.sh`,
    `environment.yml`, checksums)

### GUI baseline in place

- Main lineage page with table/graph/container views.
- Main lineage page content is vertically scrollable so graph height does not
  hide container/arrangement sections.
- Main lineage pane now supports a draggable split between the
  `Table`/`Graph` area and the `Containers`/`Arrangements` area.
- Main lineage node-group baseline:
  - disjoint node-group model (strict one-group-per-node membership)
  - table indentation under group representative rows
  - graph enclosure outlines for grouped nodes
  - collapsed groups project external edges to representative nodes
  - table/graph node context-menu actions for mark/draft/create workflows
- Graph contract:
  - single-click selects
  - double-click opens sequence/pool context
  - hover exposes node context (pool range and ladder hints for pool nodes)
  - operation transitions render as square intermediate glyphs between node
    circles, with compact operation symbols and optional edge-label text
  - projected parallel transitions (common with collapsed groups) are
    coalesced into one visible edge with grouped operation semantics
  - collapsed group representatives display hidden-internals badges (total
    hidden operation count + top hidden operation-family chips)
  - serial arrangements are rendered as dedicated graph nodes linked from lane-source sequences
  - macro instances are rendered as dedicated box nodes with explicit
    input/output edges (sequence/container-resolved where applicable)
  - dotplot/flexibility analysis artifacts are rendered as dedicated analysis
    nodes linked from source sequence nodes
  - lineage SVG export from the main GUI now reuses the same grouped +
    Gibson-hub projection model shown in the graph view, so saved DALG output
    matches the visible `Gibson cloning` operation node instead of the older
    raw multi-edge serialization
  - engine/CLI `RenderLineageSvg` now also projects Gibson apply operations
    through that same dedicated `Gibson cloning` hub so command-line SVGs no
    longer fall back to the impossible raw multi-edge form
  - projected Gibson hub SVG export now centers the hub label inside the box
    and suppresses redundant connector-edge labels for cleaner tutorial/hero
    figures
  - graph-canvas context menu now includes `Save Graph as SVG...`
  - `Save Project...` / lineage-SVG save dialogs now default filenames from the
    current project name/path rather than fixed placeholders
  - interim release-readiness pass now also covers:
    - explicit Gibson specialist/help/tutorial guardrail text for the current
      multi-insert `defined opening` requirement
    - versioned root release notes for `v0.1.0-internal.2`
    - `docs/release.md` local pre-tag smoke checklist matching the
      `script-interfaces` packaging build
    - task-oriented readiness matrices in `README.md`, `docs/gui.md`, and
      `docs/tutorial/README.md` so biologist-facing users can distinguish
      recommended, caveated, and exploratory paths without reading the roadmap
- Dense rendering controls including regulatory placement policy and visibility
  toggles persisted through display state.
- Linear DNA-letter rendering controls now include:
  - adaptive viewport-density routing (default-on) with explicit mode selector:
    `auto-adaptive`, `force-standard`, `force-helical`, `force-condensed-10`
  - compressed-letter toggle semantics scoped to auto mode
    (`linear_sequence_helical_letters_enabled`)
  - deterministic tier routing in auto mode (`1.5x / 2x / 10x` capacity tiers)
    with explicit `OFF` when condensed capacity is exceeded
  - modulo-10 seam offset control for row alignment (`(bp + offset) % 10`)
  - condensed layout backbone replacement and outward feature-lane clearance
  - configurable double-strand display with optional 180° reverse-letter rotation
  - configurable helical strand geometry (`parallel` vs `mirrored`) for
    reverse-strand visual semantics
  - continuous-helical rendering now uses oscillating phase geometry with
    endpoint-anchored width projection (full viewport span retained while
    avoiding monotonic "ramp" artifacts in 1.5x..2.0x density range)
  - helical projection smoothing now interpolates inside compressed columns to
    avoid local x-plateaus/crowding at curve extrema and improve symmetry
  - continuous-helical rendering now applies cycle-level arc-length
    compensation for more uniform bp-to-bp spacing along the helical curve
  - configurable reverse-strand letter opacity now de-emphasizes reverse rows
    consistently across linear map and sequence panel views
  - when backbone hiding is active and letters are visible, baseline tick marks
    are suppressed alongside the backbone line for clearer dense views
  - optional sequence-panel auto-hide when map letters are visible
  - one shared routing helper used by renderer diagnostics and Sequence-window
    status/auto-hide decisions (UI-status parity)
- Linear-map drag selection can now be extracted directly to a new sequence via
  `ExtractRegion`, preserving overlapping features in the derived fragment.
- Sequence-window `Export View SVG` now includes profile-aware export routes:
  - default `screen` profile for current-window composition
  - `wide-context` profile for larger canvas and expanded linear bp context
  - `print-a3` profile with A3 landscape physical-size metadata and higher
    sequence text capacity
  - debug builds include routing-tier diagnostics in the SVG header block
- Scrollable/resizable Engine Ops area and shared-shell panel.
- DNA sequence windows now expose direct genome-anchor extension controls
  (`Extend 5'` / `Extend 3'` with bp + optional output ID) next to anchor
  status, using shared `ExtendGenomeAnchor` engine semantics (same as shell/CLI).
- File menu now includes a dedicated `UniProt Mapping...` specialist dialog
  (fetch/import/project) backed directly by shared UniProt engine operations
  (`FetchUniprotSwissProt`, `ImportUniprotSwissProt`,
  `ProjectUniprotToGenome`) and a recent-entry table for quick `entry_id`
  reuse.
- File menu now also includes `Fetch GenBank Accession...` specialist dialog
  backed by shared engine operation `FetchGenBankAccession` (accession + optional
  `as_id`, with imported sequence window auto-open).
- `File -> Open Sequence...` now supports multi-file selection and imports the
  chosen files sequentially through the same per-file `LoadFile` path, opening
  one sequence window per successful import.
- Linear sequence windows now include a dedicated primary `Splicing map` mode
  (read-only) for selected `mRNA`/`exon` features.
  - Primary map splicing lanes reuse the same `SplicingExpertView` payload and
    lane-geometry renderer used by the splicing expert window.
  - Primary splicing lanes are clickable and focus transcript features in the
    sequence view.
  - `Export View SVG` now preserves this mode via the same splicing payload
    renderer (with deterministic placeholder output when no splicing feature is selected).
  - No duplicate biology logic was introduced; the mode is presentation-only.
- Context-sensitive hover descriptions on actionable controls.
- Unified interaction policy baseline is implemented across primary canvases and panes:
  - default wheel/trackpad scroll pans/scrolls
  - zoom on canvases is `Shift + wheel`
  - hand-pan drag mode is `Option` (Alt) + drag
  - legacy graph aliases (`Cmd/Ctrl + wheel`, `Space + drag`) remain
    transitionally enabled
  - keyboard pane scrolling (`Arrow`, `PageUp/PageDown`, `Home/End`) is enabled
    for major scroll panes (help, lineage tables/lists, jobs/history,
    configuration, and sequence-window feature/detail panes)
  - linear DNA map now supports two-axis pan:
    - horizontal pan for bp viewport
    - vertical pan for rendered feature-lane stack
    - explicit map-side vertical pan slider + quick `0` fit action
  - linear DNA toolbar fit actions are split:
    - `Fit Seq` for full-sequence horizontal fit
    - `Fit Features` for vertical recenter of current subsequence
- Help-window shell command reference generated from `docs/glossary.json` with
  interface filter controls (`All`, GUI shell, CLI shell, CLI direct, JS, Lua).
- Help open-path latency hardening:
  - opening help now reuses already-loaded markdown/glossary payloads
    (no unconditional disk reload/parse on each open)
  - when Help is already open on the requested tab, re-open now takes a
    focus-only fast path (no redundant search refresh)
  - explicit `Reload` in help remains the deterministic path for on-disk doc refresh
  - help search/render paths no longer clone full markdown buffers on each
    refresh/frame
- Window-backdrop load path now uses a stable texture-size hint so opening
  differently sized windows does not trigger per-window-size image reload paths.
- Slow-open diagnostics now emit status-bar timing hints when:
  - native macOS `GENtle Help...` menu dispatch to app-loop consumption exceeds
    a threshold
  - Configuration runtime sync exceeds a threshold
  - Help payload load exceeds a threshold
  - viewport focus acquisition exceeds a threshold
  - first-frame Help/Configuration render exceeds a threshold
  - total Help/Configuration window activation exceeds a threshold
  - native macOS open-window menu synchronization exceeds a threshold
  - native macOS window-menu synchronization is now deferred while an
    in-progress Help/Configuration open probe is active (sync resumes on next frame)
- Help menu now includes a dedicated `Reviewer Quickstart` guide with internal
  preview warnings, known limitations, and deterministic colleague walkthrough
  steps.
- BLAST UX/provenance hardening baseline:
  - BLAST dialog now emits live query heartbeat status (elapsed runtime while
    `blastn` is running) in addition to query-count progress.
  - status/result views show explicit invocation template and resolved command line.
  - BLAST-hit import operations now persist invocation metadata in
    `ImportBlastHitsTrack.blast_provenance` for operation history/lineage context.
  - Async BLAST scheduler now uses bounded FIFO dispatch with explicit
    `queued`/`running` job-state transitions and scheduler metadata in
    `blast-start|status|list` responses.
  - agent-suggested shell commands can execute shared BLAST routes
    (`genomes/helpers blast`, `genomes/helpers blast-track`); recursion guardrail
    blocks only nested `agents ask`
- Native macOS menu mirrors of open windows are available under:
  - `Window -> GENtle Open Windows…`
  - `GENtle -> GENtle Windows…`
  - entries now sync by stable viewport keys (not transient index position)
  - selecting an entry requests focus for the specific window key
  - active window state is mirrored with native menu checkmarks
- Circular sequence-map feature labels now use collision-aware placement:
  labels can slide within feature spans and avoid overlap with already rendered
  labels (for example restriction-site annotations).
- GC-content overlays now use a configurable bin size (`gc_content_bin_size_bp`)
  shared across linear view, circular view, and SVG export.
- Help viewer image handling:
  - markdown images render at constrained width
  - image captions are authored inline in markdown (`*Figure: ...*`)
  - CommonMark rendering path is re-enabled on the current `egui/eframe` stack
    (`egui_commonmark 0.22` aligned with `egui 0.33`), replacing the temporary
    plain-text fallback
- Experimental window backdrop styling path:
  - optional per-window-type accent tint (`main`, `sequence`, `splicing`,
    `pool`, `configuration`, `help`)
  - optional image watermark path per window type
  - configuration UI now includes tint color pickers, image-file pickers,
    and live path validation
  - persisted in app settings and live-applied from Configuration -> Graphics
- Configuration apply workflow:
  - staged settings now use one consistent apply model across external tools,
    graphics, and window styling
  - explicit `Unapplied changes` indicator is shown in Configuration
  - bottom `Cancel`/`Apply` actions are persistent and remain visible while tab
    content scrolls

### High-level parity snapshot

| Area | Status |
|---|---|
| Core cloning/editing operations across GUI/CLI/JS/Lua | Done |
| Export/render operations (sequence/lineage/pool gel) | Done |
| Reference genome + track import surfaces | Done |
| Shared shell parity across GUI/CLI | Done |
| Feature expert views (TFBS/restriction/splicing) via shared engine model | Done |
| Candidate strand-relation controls across adapters | Done |
| Alternative-splicing interpretation (lanes/boundaries/events/evidence) | Done (expert + primary-map read-only baseline) |
| Cloning-mode macro presets (SnapGene-style workflows) | Partial (starter templates only) |
| AI communication routes (agent bridge + MCP server) | Partial (agent bridge + guarded MCP op/workflow baseline implemented; ClawBio skill scaffold available) |
| Gel simulation realism and arrangement modeling | Partial |
| Shared operation protocol usage | Partial |

Notes:

- Detailed per-operation parity is code-backed; use
  `cargo run --bin gentle_cli -- capabilities` and
  `cargo run --bin gentle_cli -- state-summary` for current runtime inventory.
- `import-pool` and some resource utilities remain adapter-level contracts.

## 2. Active known gaps (priority-ordered)

1. Cloning routine standardization is incomplete:
   - typed cloning-routine catalog + template-port preflight baseline is now
     integrated into macro validation and lineage visualization, but semantic
     depth is still incomplete
   - richer preflight semantics are now baseline (generic + Gibson +
     restriction + Golden Gate + Gateway + TOPO + TA/GC + In-Fusion +
     NEBuilder HiFi), but deeper protocol-specific constraints/edge cases
     remain incomplete (including deeper restriction variants)
   - macro-instance lineage recording now covers success and failure/cancel and
     supports shell introspection, but replay-oriented helpers are still pending
     (for example per-instance re-run from recorded bindings)
   - macro-node drill-down now has a persistent detail panel in GUI lineage
     view, but dense-view controls are still pending (edge-density/aggregation)
   - protocol-family template packs now include baseline coverage for
     restriction, Gibson, Golden Gate, Gateway, TOPO, TA/GC, In-Fusion, and
     NEBuilder HiFi; advanced variants and richer family-specific templates are
     still incomplete
   - cross-tool benchmarking (Serial Cloner + MacVector + SnapGene synthesis)
     confirms additional repeated gaps are not yet first-class:
     primer design/validation workflows, auto-annotation library scans,
     sequencing-confirmation workflows, and interactive cloning clipboard/model
   - PCRtools parity intake (paper-first review, 2026-03-20; DOI:
     `10.1016/j.omtn.2025.102716`) highlights additional assay-family gaps not
     yet first-class in shared engine contracts:
     - inverse PCR on circular templates
     - bisulfite-aware PCR/LAMP assay design modes
     - KASP/PACE/ASQ allele-specific genotyping assay workflows
     - LAMP primer-set design (FIP/BIP/F3/B3 with optional loop primers)
     - custom multiplex tiling panel design + primer pooling strategy outputs
     - virtual PCR/off-target search with mismatch-tolerance reporting
   - primer design backend parity is still incomplete:
     - Primer3 backend baseline is now available behind `DesignPrimerPairs`
       (with deterministic auto-fallback to internal backend), but deeper
       constraint-mapping parity still needs hardening
     - backend preflight/version diagnostics are now shared across engine/shell
       (`primers preflight`) and GUI Engine Ops preflight controls, but broader
       cross-tool diagnostics for additional external binaries remain pending
     - adapter-equivalence coverage is still incomplete:
       - shared-shell parity baseline exists for internal-vs-Primer3 report
         normalization/provenance
       - CLI forwarded-route parity now covers `primers preflight`
       - wider matrix coverage is still pending across additional primer routes
         and adapters
2. MCP route now has guarded mutating execution (`op`/`workflow`) and
   UI-intent parity baseline (`ui_intents`, `ui_intent`,
   `ui_prepared_genomes`, `ui_latest_prepared`), but broader parity breadth is
   still incomplete (additional shared-shell route coverage, richer result
   contracts for future mutating UI intents, and wider adapter-equivalence
   coverage).
3. Mutating-intent safety policy is not yet fully hardened across agent, voice,
   and MCP invocation paths.
4. Async long-running command orchestration is still incomplete:
   - BLAST async job-handle/progress/cancel baseline is now available through
     shared shell (`genomes/helpers blast-start|status|cancel|list`) and MCP
     (`blast_async_start|status|cancel|list`)
   - async BLAST queueing/concurrency baseline is now implemented (bounded FIFO
     scheduler with configurable max concurrency)
   - agent auto-execution still needs higher-level orchestration for polling and
     multi-step async flows (it currently executes one suggested command at a
     time)
   - upcoming primer-pair selection workflows are expected to fan out into
     multiple BLAST searches; this still needs workflow-level async orchestration
     on top of the baseline job primitives
5. Core architecture parity gaps remain:
   - some utilities are still adapter-level rather than engine operations
     (notably `import-pool` and resource-sync utilities)
   - no dedicated engine operation yet for exporting a full run/process as a
     technical-assistant protocol text artifact
   - view-model contract is not yet formalized as a frontend-neutral schema
6. guideRNA workflow remains incomplete (guide-candidate model, oligo
   generation/export, macro template flow; draft in `docs/rna_guides_spec.md`).
7. XML import follow-up remains for `INSDSet/INSDSeq` dialect support.
8. Visualization and workflow UX gaps remain:
   - chromosomal-scale BED overview/density view is missing
   - genome-extract failure diagnostics now include alias-aware guidance
     (`17` vs `chr17`) plus prepared-contig previews/suggestions from the
     active reference; retrieval status now exposes a dedicated
     one-click `Apply suggested contig` action
   - GUI-driven feature editing is not yet first-class:
     - no explicit edit workflow yet for exon/intron/transcript boundary curation
       informed by RNA-seq interpretation
     - feature edits should be performed in sequence/splicing-specialist
       contexts (not in the main lineage/project window)
     - engine-backed edit operations + provenance model are still pending
   - adaptive linear DNA-letter routing baseline is implemented (shared
     renderer/UI decision path, seam-offset behavior, condensed backbone
     replacement, deterministic annotation-clearance tests); dense
     snapshot/readability benchmarking remains pending
   - any screenshot-based readability baseline artifacts require manual human
     contribution while agent screenshot execution remains policy-disabled
   - unified zoom/pan behavior is now implemented for map/graph/help/list panes;
    focused-region fallback behavior and wider regression coverage are pending
   - UI-level snapshot tests for feature-tree grouping/collapse are pending
   - backdrop-image readability guardrails and stricter grayscale handling are
     incomplete
   - optional all-open-window refresh behavior for applied display changes
     should be finalized consistently across relevant graphics-setting changes
   - sequence-panel amino-acid row rendering is intentionally deferred until
     engine-owned transcript/CDS-aware translation contracts are implemented
     (explicit codon-table resolution plus frame/phase context); current dormant
     AA-row path is maintained as safe no-op
   - translation-speed modeling is not yet implemented:
     - it must become genome-aware/organism-aware rather than relying only on
       the translation table used to decode codons
     - translation of a sequence should be able to consult reference codon-usage
       and translational-speed context for the relevant genome/host
     - target use cases:
       - bacterial protein-production work where codon choice can affect
         expression/translation pace
       - human interpretation workflows where presumed coding effects should be
         compared against a human reference translation-speed baseline
   - tutorial/executable-guidance UX still needs an explicit per-step checklist
     with state-verifiable completion markers inside GUI walkthrough flows
     (including promoter cloning walkthroughs such as TP73)
     - PCR selection-first chapter and TP73 Step 5 now include explicit
       queue-capture and batch-results checklists; broader chapter coverage
       remains pending, although the Gibson specialist testing tutorial now
       provides that style of checklist for the new destination-first Gibson
       flow
    - tutorial authoring is still structurally split:
      executable chapter metadata now comes from `docs/tutorial/sources/`, while
      GUI and agent/reference tutorials still keep their narrative content in
      hand-written markdown
      - current mitigation is the canonical landing page
        `docs/tutorial/README.md`
        plus generated discovery metadata in `docs/tutorial/catalog.json`
        sourced from `docs/tutorial/sources/`
      - before first release, run one tutorial-consistency sweep across manual
        and generated tutorial pages to verify current GUI menu/control names,
        GUI/CLI parity wording, and cloning/design tutorial cross-links remain
        aligned
      - future improvement should reduce drift between source units and the
        hand-written narrative bodies, rather than continuing to depend on a
        single large manifest
   - TP73 cDNA-vs-genomic dotplot tutorial now has an explicit screenshot
     coverage checklist in `docs/tutorial/two_sequence_dotplot_gui.md`; pending
     additions are focused on dotplot-stage captures (`10..11`) to complete the
     visual walkthrough
   - canonical reusable workflow baseline is now added for this route:
     `docs/examples/workflows/tp73_cdna_genomic_dotplot_online.json`
     and linked as chapter `13` in `docs/tutorial/manifest.json`
   - tutorial quality gate now includes deterministic image-link existence
     coverage for TP73 cDNA-vs-genomic tutorial markdown
     (`tutorial_tp73_cdna_genomic_markdown_image_links_exist` in
     `src/workflow_examples.rs`)
   - publication/release visual-polish mode is still pending:
     - add a dedicated `Publication mode` preset for GUI + SVG export paths
     - include deterministic readability-focused defaults (backdrop strength,
       debug overlay visibility, reverse-strand emphasis, typography/spacing)
     - treat exact figure-style tuning as intentionally iterative while core
       functionality and workflows continue to evolve
9. Cross-application clipboard interoperability for sequence + feature transfer
   is not yet implemented (current baseline is deterministic in-app extraction).
10. Screenshot bridge execution remains disabled by security policy.
11. Auto-updated documentation with embedded graphics remains postponed.
12. RNA mapping benchmark-fixture curation is still incomplete:
   - compact committed pack now exists at `test_files/fixtures/mapping/` and
     is used by deterministic TP73 seed-filter tests
   - legacy `test_files/mapping/True_TP73/` and
     `test_files/mapping/False_TP73/` corpora still
     need full pinned-source provenance capture before broader CI adoption
   - follow-up: add a TP53-forward benchmark profile (higher-abundance,
     TP73-family-adjacent) while keeping fixture footprint small
13. cDNA presentation semantics are partially enforced (follow-up remains):
   - implemented:
     - cDNA mode now defaults to intrinsic evidence by hiding contextual
       transcript-projection features (`mRNA`/`exon`/`CDS`) unless explicitly
       enabled
     - explicit cDNA-only opt-in toggle is available in toolbar
       (`Ctx mRNA`), default off
     - genomic annotation projection now tags generated context features with
       deterministic qualifiers
       (`gentle_generated=genome_annotation_projection`,
       `gentle_context_layer=contextual_transcript|contextual_gene`)
   - remaining:
     - continue refining intrinsic-vs-context split for broader mapped layers
       beyond transcript projection
     - wire dotplot/flex outputs into the same intrinsic-vs-context visibility
       contract
     - keep "all contextual transcripts" contract stable across reopen/reload
       and expand fixture coverage for legacy projects without context tags

### MCP server communication track (UI-intent parity baseline implemented)

Goal: add MCP as a first-class AI communication and capability-negotiation
route while keeping one deterministic engine contract across all adapters.

Current baseline:

- `gentle_mcp` stdio server is implemented.
- baseline tools implemented:
  - `capabilities`
  - `state_summary`
  - `op` (shared engine operation execution; requires explicit `confirm=true`)
  - `workflow` (shared engine workflow execution; requires explicit `confirm=true`)
  - `help`
  - `ui_intents`
  - `ui_intent`
  - `ui_prepared_genomes`
  - `ui_latest_prepared`
  - `blast_async_start`
  - `blast_async_status`
  - `blast_async_cancel`
  - `blast_async_list`
- successful mutating calls persist state to the resolved `state_path`.
- UI-intent MCP tools now route through shared shell parser/executor paths and
  are enforced as non-mutating.
- deterministic adapter-equivalence tests now assert MCP-vs-shared-shell parity
  for all current UI-intent tools plus async BLAST status routing.

Planned work:

1. Extend adapter-equivalence tests (CLI shell vs MCP tool invocations) for key
   cloning flows beyond UI-intent helpers.
2. Keep structured schema compatibility clear across JSON-RPC envelopes and
   MCP tool result payloads.
3. Keep zero MCP-only biology logic branches.
4. For future mutating UI-intent tools, require explicit confirmation and
   preserve adapter-equivalent execution reports (`executed`,
   `resolved_target`, warnings/errors, follow-up choices).

### Alternative-splicing interpretation track (baseline implemented; follow-ups)

Goal: make exon/intron boundary interpretation and alternative-splicing
inspection explicit, readable, and adapter-equivalent.

Status:

1. Implemented baseline:
   - shared engine-owned splicing payload (`SplicingExpertView`) and selector
     (`FeatureExpertTarget::SplicingFeature`) are in place.
   - boundary markers, event summaries, transcript-vs-exon matrix, and
     junction-support arcs are rendered from the shared payload.
   - predicted exon->exon transition matrix is rendered in GUI and SVG, with
     frequency-coded support cells and exon `len%3` header cues (heuristic
     frame signal).
   - CDS-aware exon flank phase coloring (`0/1/2`) is rendered in GUI and SVG
     from shared splicing payload data when transcript `cds_ranges_1based`
     qualifiers are present.
   - sequence-window primary `Splicing map` mode (read-only) is now available
     and uses the same shared payload/geometry path as the expert view.
   - GUI expert panel and SVG export use the same payload; shell/CLI/JS/Lua
     can inspect/export via existing feature-expert commands.
2. Remaining follow-ups:
   - add dense-fixture snapshot tests for boundary visibility/non-overlap.
   - Add deterministic tests for crowded transcript sets to ensure geometry/readability stays stable (no exon/intron distortion, no clipping).
   - add explicit geometry invariants to ensure label text can never alter
     exon/intron footprints.

### GUI-driven feature editing track (new)

Goal: make annotation curation (especially exon/intron/transcript structure)
editable through deterministic engine operations with full provenance and
cross-adapter parity.

Status:

1. Current baseline:
   - strong read-only interpretation exists (feature map + splicing expert).
   - no first-class engine operation set for editing existing feature geometry
     or qualifiers.
2. Required implementation:
   - define engine operations for feature editing (add/update/delete) with
     explicit geometry and qualifier contracts.
   - include reversible provenance payloads in operation history:
     - before/after feature snapshots
     - reason/evidence metadata (for example RNA-seq-supported curation note)
   - expose identical edit operations in GUI, shared shell/CLI, JS, and Lua.
3. GUI scope rule:
   - editing must live in sequence/splicing-specialist windows where local
     coordinate context is visible.
   - the main lineage/project window remains orchestration/navigation focused
     and should not become a feature geometry editor.
4. RNA-seq alignment path (later in this track):
   - import/overlay RNA-seq-derived splice evidence as suggestions.
   - offer user-confirmed "apply suggestion as feature edit" actions via the
     same engine edit operations (no GUI-only mutation path).

### RNA-seq evidence for cloning-candidate regions (Nanopore cDNA first)

Goal: prioritize cloning-oriented interpretation for small genomic regions of
interest, with deterministic RNA-seq evidence overlays and explicit
suggestion-first editing handoff.

Detailed execution plan: `docs/rna_seq_nanopore_cloning_regions_plan.md`.

Detailed origin-classification + sparse-index follow-up plan:
`docs/rna_read_origin_sparse_index_plan.md`.

Status:

- In progress (phase-1 baseline implementation started).
- This replaces earlier demo-oriented framing for this direction; immediate
  focus is production-aligned ROI workflows for cloning-candidate regions.
- Implemented in current baseline:
  - Splicing Expert `Nanopore cDNA interpretation` run path is now asynchronous
    (non-blocking UI) with live progress updates.
  - Splicing Expert now exposes a regular workflow access route for the same
    cDNA mapping payload:
    `Prepare Workflow Op` stages `run_id + ops` in Engine Ops workflow runner,
    and `Copy Workflow JSON` exports a complete workflow object for CLI/shell
    execution.
  - Running seed-confirmation histogram (genomic-position bins) is shown during
    execution and updates every 1000 reads (`+` strand up, `-` strand down).
  - Histogram bars now use sqrt-scaled heights so low-frequency bins remain
    visible instead of collapsing under one dominant bin.
  - Seed-hash catalog preview is now available directly in Splicing Expert and
    overlaid as red position dots with hover sequence/hash detail.
  - Live top-read preview is now available during streaming; selected rows
    highlight their supported seed positions in green with recompute-time
    telemetry for in-window hash-speed sanity checks.
  - Streaming status now shows ETA estimated from gzip/plain input bytes
    consumed and elapsed runtime; ETA refresh is synchronized with read-count
    update cadence.
  - Seed-hit ranking now includes inverse-occurrence weighting to dampen
    low-complexity repetitive-seed dominance in top-read previews.
  - Phase-1 seed-pass logic now uses a composite gate:
    `raw >= min_hit AND weighted >= min_weighted AND unique >= min(min_unique, tested_kmers)`,
    with GUI/shell controls for all thresholds.
  - Added transcript-chain spacing guard to reduce dispersed false positives:
    reads must also satisfy `median transcript seed-gap <= max_median_transcript_gap`
    (default `4.0`), with exported per-read gap diagnostics.
  - Added transition-aware specificity guard:
    reads must also satisfy `confirmed_transitions >= min_confirmed_transitions`
    and `confirmed_transition_fraction >= min_transition_support_fraction`.
  - Added isoform-support ranking panel during runs:
    one row per known transcript (label/id/strand) with assigned reads,
    seed-pass counts, transition coverage, and auto-picked best isoform.
  - Added strand-audit diagnostics to the joint-run isoform/read display:
    isoform rows now expose chain-same/opposite-competition/ambiguity counts,
    and top-read rows expose selected strand plus opposite-strand competition
    and ambiguity flags.
  - Added origin-classification schema scaffolding (non-breaking):
    per-read origin class/reason/confidence fields, candidate-contribution
    hints, and running/final origin-class count summaries in reports/progress.
  - Activated sparse multi-gene template expansion in engine runtime:
    `InterpretRnaReads` with `origin_mode=multi_gene_sparse` now expands
    transcript-template indexing using local-annotation matches from
    `target_gene_ids[]` (deterministic, shared across GUI/CLI/shell/agent
    routes).
  - `roi_seed_capture_enabled` remains a documented deterministic no-op with
    explicit warning until the ROI capture layer is implemented.
  - Added non-breaking report-compaction controls:
    `InterpretRnaReads` and `rna-reads interpret` now accept `report_mode`
    (`full|seed_passed_only`) so persisted retained-hit payloads can be
    compacted to seed-pass rows only while preserving full run counters.
  - Added deterministic checkpoint/resume scaffolding for long runs:
    `checkpoint_path`, `checkpoint_every_reads`, and
    `resume_from_checkpoint` are now supported in engine and shell routes;
    checkpoint snapshots are serialized to
    `gentle.rna_read_interpret_checkpoint.v1` and resume restores counters,
    support tables, retained hits, and score-density bins.
  - Splicing Expert advanced Nanopore controls expose active sparse-origin
    settings (`origin mode`, `target genes`, `ROI seed capture`) and pass them
    through unchanged to engine/shell contracts.
  - Deterministic TP73 seed-filter tests now use compact committed mapping
    fixtures (`test_files/fixtures/mapping/`) and cover:
    - expected TP73 positive behavior (TP73-derived reads with 30% deletions
      still pass)
    - close-family negative control (human TP53 set rejected for TP73 model)
    - generic negative control (same-length deterministic random reads fail)
  - Histogram guide overlays are now user-toggleable (`Exons`, `Introns`) for
    clearer exon-context interpretation during filtering runs.
  - Histogram coordinate mode now supports genomic axis and exonic-only compact
    axis (merged exons adjacent, introns removed from coordinate span).
  - Empty RNA report IDs now default deterministically from input file names
    (`cdna_<filename_stem>`).
  - Inline per-read alignment is disabled for phase-1 Nanopore runs so
    streaming progress stays responsive.
  - Phase-2 retained-read alignment is now implemented via shared operation
    `AlignRnaReadReport` (shell: `rna-reads align-report`; GUI:
    `Run alignment phase (retained report)`), refreshing per-hit mappings,
    `msa_eligible` diagnostics, transition/isoform support, and exon/junction
    abundance frequencies.
  - Phase-2 alignment now deterministically re-ranks retained hits by
    alignment-aware retention rank (alignment identity/coverage/score first,
    seed metrics as tie-breakers) so downstream candidate ranking reflects
    mapping quality.
  - Added non-mutating alignment inspection route:
    `rna-reads inspect-alignments REPORT_ID [--selection ...] [--limit N]
    [--effect-filter ...] [--sort ...] [--search TEXT]
    [--record-indices i,j,k]` returns ranked aligned-hit rows for report triage
    without mutating state, with a structured subset spec echoed in JSON for
    reproducible agent-driven inspection.
  - `Mapped cDNA` inspection in Splicing Expert is now read-first by default:
    saved-report aligned rows expose phase-1 vs phase-2 comparison,
    deterministic effect labels (`confirmed`, `reassigned`, `no phase-1 tx`),
    and per-read mapped exon/junction contribution details, while the old
    top-hit list remains a clearly labeled capped live preview.
  - Added ranked alignment TSV export for downstream tabular analysis:
    `ExportRnaReadAlignmentsTsv` /
    `rna-reads export-alignments-tsv REPORT_ID OUTPUT.tsv
    [--selection ...] [--limit N]`.
  - Added dotplot-style alignment inspection export:
    `ExportRnaReadAlignmentDotplotSvg` /
    `rna-reads export-alignment-dotplot-svg REPORT_ID OUTPUT.svg
    [--selection ...] [--max-points N]` renders coverage-vs-identity scatter
    with score-based coloring.
  - Phase-1 seed filtering now hashes full read span for every read (replacing
    prior sampled-window behavior) to improve filtering sensitivity.
  - Splicing Expert now exposes the active phase-1 hash-density control
    (`seed_stride_bp`, default `1`) so seed-start spacing can be tuned without
    diverging from engine/CLI/shell contracts.
  - Runtime panel now reports detailed compute breakdown
    (`seed`, `align`, `io`, `parse`, `norm`, `infer`, `emit`, `other`) and
    throughput metrics (`reads/s`, `bp/s`, mean/median/p95 read length) to
    make overhead diagnosis explicit.
  - Runtime/summary status now reports seed-pass percentages (not only raw
    counts) for faster triage on large batches.
  - Seed-hit score-density panel now supports `Linear`/`Log` rendering
    (`Log` default) so sparse higher-score bins remain visible in heavily
    skewed distributions while still allowing raw-count inspection.
  - Splicing Expert now exposes explicit cDNA/direct-RNA interpretation mode
    control (`Input is cDNA` checkbox) with configurable poly-T prefix minimum
    for automatic reverse-complement normalization.
  - Splicing Expert / `InterpretRnaReads` now admit `ncRNA` seed features in
    addition to `mRNA`/`transcript`/`exon`, so reverse-strand antisense loci
    such as `TP73-AS3` can use the shared cDNA filtering path without
    feature-kind workarounds.
  - cDNA normalization now uses tolerant 5' T-rich head detection (minor
    interruptions in poly-T tails are accepted) to better normalize real
    Nanopore cDNA tails before seed scoring.
  - Added coherent-chain specificity guard:
    reads must satisfy `chain_consistency_fraction >= min_chain_consistency_fraction`
    (default `0.40`) so dispersed local-seed matches are less likely to pass.
  - RNA-read seed-filter defaults now use `kmer_len=10` (was `9`) to reduce
    broad false-positive matches in large human transcriptome cDNA screens.
  - Exon-transition support counters now update only for reads that pass the
    full seed gate; rejected reads no longer inflate transition support tables.
  - RNA-read hit exports now include reverse-complement provenance markers
    (`rc_applied`) in FASTA headers and TSV path exports.
  - New RNA-read TSV exports are available for downstream analysis:
    - per-read exon-path table (`rna-reads export-paths-tsv` / GUI export)
    - exon/transition abundance table (`rna-reads export-abundance-tsv` / GUI export)
  - RNA-read score-density chart export is now available as SVG from GUI/CLI
    (`Export Score Density (SVG)...` /
    `rna-reads export-score-density-svg ... --scale linear|log`).
  - Splicing Expert now shows seed-confirmed exon-exon transition support
    tables with per-transition read counts/percentages and indexed
    junction-crossing seed-bit diagnostics.
  - RNA-read reports persist exon-support and exon-exon junction-support
    frequency schema fields; phase-1 seed-only runs populate placeholders until
    the explicit phase-2 alignment step is run.
  - Splicing Expert RNA-read support tables now separate `Seed` diagnostics
    from `Mapped` phase-2 support so exon/junction/isoform interpretation can
    follow retained-read alignments instead of reusing seed/path heuristics.
  - Mapped exon/junction support counting now follows aligned
    transcript-template offsets, reducing false support inflation from
    alternative exons or junctions that lie inside the same broad genomic span
    but were not actually traversed by the retained-read mapping.
  - `rna-reads export-hits-fasta` now includes exon-path annotations in FASTA
    headers (`:` confirmed adjacent transition by seeds, `-` unconfirmed).
  - `rna-reads export-sample-sheet` / GUI sample-sheet export produce TSV
    cohort summaries with per-report frequency JSON columns for downstream
    annotation workflows.
  - RNA-read report listing/sample-sheet exports now carry sparse-origin
    request provenance (`origin_mode`, target-gene counts/IDs, ROI-capture
    flag) so multi-gene scaffold runs remain auditable in cohort-level tables.
  - `rna-reads list-reports` / `show-report` shell outputs now include compact
    human-readable provenance summaries (`summary_rows[]` / `summary`) for
    faster triage without post-processing report JSON.

Track boundaries:

1. Initial data path is Nanopore cDNA evidence with deterministic
   filtering + partial ROI mapping (not full transcriptome alignment in
   GENtle at this stage), with FASTA-first ingest in phase 1.
2. Evidence modeling remains engine-owned and adapter-equivalent across
   GUI/shell/CLI.
3. Curation flow remains non-mutating by default and suggestion-first, with
   explicit user-confirmed apply actions only.
4. Transposon-specific interpretation/modeling is explicitly deferred and
   pending Anze-led direction.
5. Add a transcriptome-scale SNR calibration subtrack:
   - treat the current `30%` seed-hit threshold as bootstrap default only
   - empirically estimate background/noise seed-hit distributions for
     full-transcriptome input
   - add online common-seed suppression (`IDF`-like penalty/backprop on
     frequently observed seed bits across reads) as an optional scoring term
     to reduce broad-domain/repetitive false positives
   - evolve toward SNR-normalized acceptance thresholds with report-level
     diagnostics and deterministic cross-adapter parity.
6. Add SRA/FASTA ingestion + storage strategy subtrack:
   - baseline phase-1 contract remains external SRA conversion and FASTA ingest
     (`.sra` is not parsed directly by GENtle in phase 1)
   - record observed storage behavior in design assumptions:
     raw FASTA can be larger than `.sra` (for example, a ~9.4 GB `.sra` to
     ~11 GB FASTA conversion), so disk planning is mandatory
   - evaluate a streaming ingestion path via `fasterq-dump --fasta --stdout`
     as an optional execution mode for large runs
   - because the Nanopore flow is two-pass (seed filter then alignment),
     define deterministic replay/spool rules up front:
     - stream mode: keep deterministic shortlisted-hit spool entries
       (read id, ordinal index, seed metrics, sequence; byte offset when
       available) and run pass 2 from this spool
     - artifact mode: keep full FASTA artifact for exact replay/audit and pass 2
   - prefer compressed persisted exports (`.fasta.gz` plus JSON report sidecar)
     for retained hit sets and summaries.
7. Extend sample-sheet cohorts with user metadata:
   - add optional sample annotation fields (condition/timepoint/replicate,
     extraction notes) in shared engine-owned sample-sheet contracts
  - keep CLI/GUI/agent parity for writing/merging these annotations with RNA
    evidence metrics.
8. Defer seed-capture workflow abstraction (planned, not implemented):
   - model the seed-hash filtering stage as a reusable biotech-style
     "seed-capture/enrichment" operation suitable for standard workflows
     (instead of only a splicing-window action)
   - keep this track documented first, then implement through shared engine
     operation contracts so GUI/CLI/agents remain parity-safe.
   - detailed note is maintained in
     `docs/rna_seq_nanopore_cloning_regions_plan.md` (phase-2 follow-up).
9. Define phase-2 mapping strategy as seed-cascade (no DP requirement):
   - shortlist hits from phase-1 seed filtering, then run deterministic
     seed-space mapping in ordered tiers: exon-only -> exon-exon junctions ->
     intronic fallback
   - preserve unmatched-seed residuals in report payloads as candidate
     sequencing-error/trans-splicing indicators for downstream clustering.
10. Add dual-strand joint-run isoform assignment hardening (planned):
   - keep one run over both strands for `all-overlap / both-strands` scope,
     but add explicit strand-partitioned score diagnostics in the isoform table
     to make cross-strand contributions auditable.
   - add deterministic tie-break policy when `+/-` isoforms compete for the
     same read-chain evidence (strand-consistent chain preferred, then score).
   - expose optional strict mode requiring seed-chain transcript assignment and
     exon-transition confirmation to come from the same strand partition.
   - extend report/export schema with per-read strand-decision diagnostics and
     mixed-strand ambiguity counters so GUI/CLI/shell/agent surfaces stay
     parity-equivalent.
   - detailed execution notes: `docs/rna_seq_nanopore_cloning_regions_plan.md`
     (Phase 2b).
11. Add RNA-seed-guided primer-design handoff (planned):
   - reuse persisted seed-hash support maps to prioritize primerable regions
     with coherent transcript evidence.
   - down-rank primer targets dominated by high-frequency/common seed hashes
     to reduce cross-gene conserved-domain artifacts.
   - keep this as guidance input only; final primer scoring/reporting remains
     under shared engine primer-design contracts.
12. Defer focal-region similarity cohort analysis (planned, not implemented):
   - treat a well-understood focal ROI as a hypothesis anchor for genome-wide
     comparison.
   - discover the most similar loci sequence-first using existing
     BLAST/genome-extraction/alignment primitives.
   - materialize a ranked cohort of comparable loci as extracted anchored
     windows.
   - keep similarity discovery sequence-first in v1 and use a fixed matched
     window for comparable loci.
   - prepare downstream evidence tracks, including CUT&RUN, to run against
     either one ROI or one stored cohort.
   - inspection should reuse structured reports, TSV export, alignments,
     dotplots, and projected tracks before any new GUI trigger is added.
   - keep this engine-owned and adapter-equivalent; do not move similarity
     logic into routine macros or GUI-only flows.
   - preserve existing cDNA filtering/mapping behavior unchanged; this deferred
     track should not trigger premature refactors in the current RNA/cDNA path.

### Isoform-architecture panel track (baseline implemented; follow-ups)

Goal: support deterministic transcript/protein architecture rendering for curated
panels (for example TP53 Figure-1-style layouts) through one shared engine path.

Status:

1. Implemented baseline:
   - new operations:
     - `ImportIsoformPanel`
     - `RenderIsoformArchitectureSvg`
   - curated panel resource schema + starter asset:
     - `gentle.isoform_panel_resource.v1`
     - `assets/panels/tp53_isoforms_v1.json`
   - shared shell/CLI command family:
     - `panels import-isoform`
     - `panels inspect-isoform`
     - `panels render-isoform-svg`
     - `panels validate-isoform`
   - GUI Engine Ops section for panel import/inspect/export
   - tutorial/example integration:
     - `docs/examples/workflows/tp53_isoform_architecture_online.json`
     - chapter `tp53_isoform_architecture_online`
   - short hardening pass completed:
     - curator-facing panel validation/report API
       (`gentle.isoform_panel_validation_report.v1`) via
       `panels validate-isoform`
     - deterministic parity test for isoform SVG output equivalence across
       render routes (`RenderIsoformArchitectureSvg`,
       `panels render-isoform-svg`, and `render-feature-expert-svg ... isoform`)
2. Remaining follow-ups:
   - add optional publication-style presets (font/spacing/legend layout) while
     keeping the engine payload and geometry deterministic.

### Adaptive linear DNA letter routing track

Goal: route linear DNA-letter rendering by viewport capacity (not fixed bp
thresholds), while preserving strict base order, deterministic compressed
layouts, and renderer/UI parity.

Current baseline:

- shared routing helper (`linear_base_routing`) is now authoritative for both
  renderer mode selection and Sequence-window status/auto-hide decisions
- deterministic auto tiers use density capacity limits:
  - `<= 1.5x`: standard
  - `<= 2x`: helical (when compressed mode enabled)
  - `<= 10x`: condensed-10 (when compressed mode enabled)
  - `> 10x`: `OFF`
- explicit mode override is available:
  - `AutoAdaptive` (default)
  - `StandardLinear`
  - `ContinuousHelical`
  - `Condensed10Row`
- compressed-letter toggle now applies to auto mode only:
  when disabled, auto mode can use standard or `OFF` only
- legacy fixed-threshold `SetParameter` knobs are compatibility-accepted as
  deterministic deprecated no-op messages (no runtime routing effect)
- phase-offset mapping remains available (`(bp + offset_bp) % 10`, clamp `0..9`)
- condensed layout keeps backbone replacement/suppression and deterministic
  outward annotation-lane clearance
- upper-right renderer diagnostics now report adaptive state and metrics
  (active mode, route policy, density, cols-fit, glyph width, reason)
- persisted configuration migration baseline:
  - schema-versioned config
  - legacy helical toggle (`false`) migrates to `true`
  - legacy layout mode migrates to `AutoAdaptive`
  - migrated settings are rewritten once after successful load

Remaining follow-ups:

1. Add snapshot-style dense readability regression assets/benchmarks
   (manual screenshot contribution path remains policy-constrained).
2. Continue condensed/helical stress testing on very large annotation density.

### XML import integration track (GenBank-first)

Goal: add XML import support without creating a second semantic model.

Status:

1. Implemented baseline:
   - deterministic runtime import detection order now includes XML fallback:
     `GenBank -> EMBL -> FASTA -> XML`.
   - `GBSet/GBSeq` sequence import is mapped to existing `DNAsequence` +
     feature qualifier structures (no format-specific biology logic branch).
   - genome annotation parser dispatch now ingests `.xml` annotation sources
     and normalizes through existing GenBank-like `GenomeGeneRecord` mapping.
   - GUI open-sequence dialog now exposes XML file filters and docs list XML
     support explicitly.
   - unsupported XML dialects (for example `INSDSet/INSDSeq`) return explicit,
     deterministic diagnostics.
2. Remaining follow-ups:
   - implement additive `INSDSet/INSDSeq` adapter without diverging semantic
     normalization.
   - extend cross-format parity fixtures/tests to include additional XML edge
     cases (multi-record, qualifier-only, interval-only locations).
   - keep large exploratory XML samples out of committed default fixtures.

### Cloning routine catalog + macro-box graph track (new)

Goal: map current cloning vocabulary to executable GENtle routines and represent
macro instances as explicit graph boxes with typed inputs/outputs.

Current baseline:

- shared workflow macro persistence and execution are implemented
  (`UpsertWorkflowMacroTemplate`, `macros template-*`, `macros run`).
- macro import supports both:
  - legacy pack file (`assets/cloning_patterns.json`)
  - per-template hierarchy (`assets/cloning_patterns_catalog/**/*.json`)
- GUI `Patterns` menu mirrors hierarchy and imports templates through shared
  shell command contracts.
- templates can be executed repeatedly with distinct bindings/output IDs in one
  project; resulting outputs persist through regular operation state.
- typed cloning-routine manifest baseline is implemented at
  `assets/cloning_routines.json` (`gentle.cloning_routines.v1`), including
  routine-family/status/tag metadata and typed input/output port declarations.
- shared-shell/CLI discovery now supports routine list/filter/search via
  `routines list ...`.
- shared-shell/CLI explainability baseline is now implemented:
  - `routines explain ROUTINE_ID [--catalog PATH]`
  - `routines compare ROUTINE_A ROUTINE_B [--catalog PATH]`
  - response schemas:
    `gentle.cloning_routine_explain.v1`,
    `gentle.cloning_routine_compare.v1`
- routine catalog rows now support additive explainability metadata
  (`purpose`, `mechanism`, `requires`, `contraindications`,
  `confusing_alternatives`, `difference_matrix`,
  `disambiguation_questions`, `failure_modes`), with initial curation for core
  cloning families.
- GUI routine discovery baseline is now exposed under `Patterns` with grouped
  family/status browsing.
- GUI staged Routine Assistant baseline is implemented at
  `Patterns -> Routine Assistant...`, backed by shared-shell explain/compare +
  template-run command paths.
- workflow template preflight baseline is implemented in shared shell:
  - `macros template-run ... --validate-only` reports typed preflight without
    mutating state
  - normal `macros template-run` now carries preflight output and blocks
    execution on preflight errors
  - preflight now validates both input and output ports and reports
    `contract_source` (`template_ports` or `routine_catalog`)
- workflow macro templates now accept optional typed `input_ports` and
  `output_ports` contracts in template metadata.
- preflight now enforces baseline cross-port semantics (alias/collision,
  container/sequence compatibility, anchor checks with bound sequence context).
- mutating `macros run` and `macros template-run` now persist lineage
  macro-instance records for success/failure/cancel
  (bound inputs/outputs + emitted op ids + status/status_message).
- shared-shell introspection for recorded macro instances is now available via
  `macros instance-list` and `macros instance-show`.
- lineage graph and lineage SVG export now include explicit macro box nodes and
  input/output edge rendering.
- GUI lineage now includes a persistent selected-macro detail pane with
  inputs/outputs and emitted operation drill-down.
- GUI lineage node context menus now include `Rename (leaf only)` and
  `Remove (leaf only)` for sequence nodes, with strict leaf guards to avoid
  downstream-edge/container breakage in this first iteration.
- routine-family coverage and deeper semantic validation are still incomplete
  beyond the current family baselines (Gibson, restriction, Golden Gate,
  Gateway, TOPO, TA/GC, In-Fusion, NEBuilder HiFi).

Planned work:

1. Extend semantic preflight from the current baseline family checks to deeper
   routine-family-specific rules (protocol-aware compatibility models and edge
   cases).
2. Expand macro-node dense-view controls (edge-density filtering/aggregation,
   compact operation summaries for large projects).
3. Expand protocol-family packs beyond current baseline coverage (advanced
   family variants, richer constraints, and higher-depth templates).
4. Harden GUI Routine Assistant:
   - decision-trace capture/export integration
   - richer comparison guidance for missing/weak metadata rows
   - adapter-parity tests for staged flow output contracts.
5. Postponed: add replay helpers for recorded macro instances only after
   routine-family preflight models and protocol-family packs are stabilized.

Concrete patch plan (architecture-aligned): routine-application assistant +
self-describing alternatives

Status (2026-03-20): metadata schema extension + shared `routines explain` /
`routines compare` command surfaces are implemented; GUI staged Routine
Assistant baseline is implemented; decision-trace export wiring now includes
deterministic run-bundle inclusion plus `disambiguation_questions_presented`,
`disambiguation_answers`, and ordered `preflight_history` capture from GUI
Routine Assistant flows (including explicit compare-stage answer entry), with
deterministic lifecycle tests now covering `preflight_failed`,
`execution_failed`, and `aborted` trace statuses; JS/Lua adapter fixtures now
assert `ExportProcessRunBundle` decision-trace payload equivalence with direct
engine execution for equivalent staged decisions, and shared-shell
`export-run-bundle` fixtures now assert the same equivalence for CLI shell
routing.

1. Metadata contract extension (`assets/cloning_routines.json`)
- Add additive explanation fields for routine-selection clarity:
  - `purpose`, `mechanism`, `requires[]`, `contraindications[]`
  - `confusing_alternatives[]`
  - `difference_matrix[]`
  - `disambiguation_questions[]`
  - `failure_modes[]`
- Keep schema additive and backward-compatible.

2. Shared explanation commands (non-mutating, adapter-neutral)
- Add:
  - `routines explain ROUTINE_ID`
  - `routines compare ROUTINE_A ROUTINE_B`
- Output one deterministic JSON payload family used unchanged by
  GUI/CLI/JS/Lua/AI/MCP renderers.

3. GUI Routine Assistant (apply flow)
- Add a dedicated assistant panel with fixed stages:
  1) goal capture,
  2) candidate routines,
  3) alternative comparison (`why this / why not`),
  4) typed parameter form (from routine ports),
  5) preflight preview (`--validate-only` equivalent),
  6) transactional execute,
  7) protocol/run-bundle export.
- Keep execution routed through existing macro/template operations only.

4. Decision-trace export for reproducibility
- Deliver a first-class deterministic trace payload for routine choice:
  - schema:
    - `gentle.routine_decision_trace.v1` (additive metadata object)
    - embedded in `gentle.process_run_bundle.v1` as `decision_traces[]`
  - capture lifecycle (ordered):
    1) assistant start + source adapter,
    2) intent/query capture,
    3) candidate snapshot,
    4) selected routine + alternatives presented,
    5) compare events,
    6) disambiguation answers,
    7) binding snapshot,
    8) preflight snapshot/history,
    9) execution attempt/outcome,
    10) export event(s).
  - status model:
    - `draft`, `preflight_failed`, `ready`, `executed`,
      `execution_failed`, `aborted`, `exported`
  - deterministic serialization guarantees:
    - stable ordering for candidate lists/comparisons/questions/answers/op IDs
    - explicit `null`/`[]` field encoding (no adapter-specific omissions)
    - normalized text capture (trim outer whitespace, preserve case/content)
  - retention/redaction baseline:
    - include explicit user-facing inputs and engine outputs only
    - avoid implicit secret/env capture
  - integration rollout:
    - phase 1 producer: GUI Routine Assistant
    - phase 2 parity: CLI/JS/Lua assistant-equivalent flows
  - acceptance tests:
    - success trace snapshot (full lifecycle)
    - preflight-failed partial trace snapshot
    - execution-failed trace snapshot
    - aborted trace snapshot (no execution attempted)
    - run-bundle inclusion + deterministic ordering snapshot
    - adapter-parity fixture asserting schema-equivalent payloads for
      equivalent staged decisions.

5. Milestone order and acceptance gates
- M1: Golden Gate pack + family preflight + explain/compare rows. (Done)
- M2: Gateway/TOPO/TA-GC packs + family preflight + confusion mapping. (Done)
- M3: In-Fusion/NEBuilder packs + overlap-family preflight + confusion mapping. (Done)
- M4: GUI Routine Assistant + decision-trace export wiring. (In progress)
- Gate per milestone:
  - one validate-only success + one deterministic failure per family/mode
  - one transactional run creating outputs + macro-instance lineage rows
  - `routines list`/`explain`/`compare` payload consistency checks.

Postponed item detail (deferred): macro-instance replay helpers

1. Deferral reason
Replay semantics are currently too volatile while routine-family constraints and
template packs are still moving. Implementing replay now would lock in
incomplete behavior and create churn in replay contracts.

2. Target command surface (when resumed)
`macros instance-replay MACRO_INSTANCE_ID [--transactional] [--validate-only] [--allow-template-drift] [--output-prefix PREFIX]`
`macros instance-diff MACRO_INSTANCE_ID [--against LAST|INSTANCE_ID]`

3. Required replay metadata additions
Store `template_schema`, `template_name`, and resolved binding payload hash.
Store rendered macro script hash and optional rendered script snapshot.
Store operation-journal span boundaries and run fingerprint for deterministic
comparison.

4. Replay execution contract
Load macro instance by ID and resolve replay source.
If template-based and template changed, fail closed unless
`--allow-template-drift` is set.
Re-run preflight with recorded bindings.
Apply optional output-id prefix remapping to avoid collisions in current state.
Execute using the same shared shell parser/executor path as normal macros.
Return a deterministic replay report containing old/new op-id mapping.

5. Diff contract
Compare original and replayed runs by operation family sequence and
created/changed sequence IDs.
Report `equivalent` when semantic outputs match after optional ID remapping.
Report structured mismatch rows when operation families or outputs diverge.

6. Failure/cancel behavior
Replay failure/cancel must append a new macro-instance lineage row with
`status=failed|cancelled` and a status message.
Replay should never mutate state when `--validate-only` is used.

7. Acceptance tests (resume criteria)
Deterministic replay of a template-backed successful macro with
`equivalent=true`.
Deterministic replay with output-prefix remap avoids ID collisions.
Template drift is rejected by default and accepted only with
`--allow-template-drift`.
Replay failure path records failed macro-instance lineage with status message.

Detailed plan and support crosswalk:

- `docs/cloning_routine_catalog_plan.md`

### Text/voice control track (new)

Goal: make GUI behavior addressable through deterministic text commands first,
then layer voice control on top of the same command plane.

Current baseline:

- `Prepared References...` exists, includes the chromosome line inspector, and
  offers confirmed per-row reinstall for stale/partial prepared genomes.
- `Prepare Reference Genome...` now keeps a structured per-step checklist visible
  through completion, adds active-step ETA for determinate byte-based work, and
  surfaces explicit reinstall-from-sources recovery when reindex detects an
  inconsistent cached install.
- Shell/agent UI-intent routing is implemented:
  - `ui intents`
  - `ui open|focus TARGET ...`
  - `ui prepared-genomes ...`
  - `ui latest-prepared SPECIES ...`
- Prepared references supports one-shot disambiguation/open flow:
  - `ui open prepared-references --species human --latest`
  - explicit `--genome-id` still overrides query-based selection

Planned work:

1. Add guarded execution policy for mutating actions invoked via text/voice.
2. Extend intent grammar/aliasing so natural language maps more reliably to
   deterministic UI-intent commands without guessing.
3. Add optional voice adapter path (STT/TTS) that emits/consumes the same
   deterministic shell/UI intent contracts through the existing agent interface
   execution path.

Scope estimate:

- Moderate (multi-phase): core shell/parser infrastructure already exists;
  missing piece is UI-intent routing + capability discovery + voice adapter.

### SnapGene-parity benchmark track (new)

Goal: align GENtle’s operation-level workflows and visual outputs with the
major protocol families and usability expectations commonly seen in SnapGene,
while keeping GENtle’s shared-engine and open-protocol architecture.

#### A) Macro templates for major cloning modes (priority)

- Restriction cloning:
  - single-fragment insertion, multi-fragment insertion, linear ligation.
- Gibson Assembly:
  - baseline shipped: two-fragment overlap planning preview + Gibson-specific
    preflight overlap diagnostics.
  - destination-first planning model is now preview-consumed for the restricted
    single-insert v1 specialist flow:
    - schema: `gentle.gibson_assembly_plan.v1`
    - preview response: `gentle.gibson_assembly_preview.v1`
    - shared shell/direct surface:
      `gibson preview PLAN_JSON_OR_@FILE [--output OUTPUT.json]`
      `gibson apply PLAN_JSON_OR_@FILE`
    - GUI specialist:
      `Patterns -> Gibson...`
  - current single-insert v1 now covers:
    - destination opening resolution (`existing_termini` or `defined_site`)
    - destination opening suggestions from unique restriction cutters, with
      MCS-linked single-cutters prioritized first and other single-cutters
      available on demand
    - cutter-derived opening coordinates resolve to cleavage windows/cutpoints
      from shared restriction-enzyme geometry
    - terminal overlap derivation
    - Gibson-specific primer suggestions
      (`5' overlap + 3' gene-specific priming segment`)
    - blocking/advisory validation
    - deterministic creation of output sequence nodes for:
      - left insert primer
      - right insert primer
      - assembled product
    - lineage operation reopen behavior for Gibson apply
    - protocol-cartoon preview/export from the same resolved plan
  - next:
    - multi-fragment Gibson planning, preview, and cartoon generation
    - richer user influence over Gibson PCR/primer design while still keeping
      the specialist high-level and Gibson-specific
    - richer semantic rescue for interrupted annotations beyond the current
      trim/join/edited-locus projection rules, especially for CDS/gene models
    - specialist option to introduce one new unique restriction-endonuclease
      cleavage site during Gibson primer/design derivation
    - circularize-fragment workflows and richer Routine Assistant handoff
  - extend the protocol-cartoon baseline from `gibson.two_fragment` to
    multi-fragment Gibson assembly, including circular assembled products and
    canonical template/binding support for the generated explanatory strip.
- NEBuilder HiFi:
  - one-insert and multi-insert assembly workflows.
- In-Fusion:
  - one-insert and multi-insert assembly workflows.
- Golden Gate:
  - Type IIS overhang-aware assembly templates.
- Gateway:
  - BP, LR, and BP+LR multi-insert workflows.
- TOPO:
  - TA TOPO, blunt TOPO, directional TOPO workflows.
- TA/GC cloning:
  - vector + PCR-fragment templates with orientation handling.

Implementation note:

- Each protocol ships as a macro template (`macros template-put`) backed only by
  shared engine/shell operations, then exposed identically in GUI/CLI/JS/Lua.
- Templates must emit auditable operation logs suitable for DALG-derived
  protocol export.
- Routine vocabulary crosswalk and graph-node rollout plan are maintained in
  `docs/cloning_routine_catalog_plan.md`.

### Cross-tool parity synthesis (Serial Cloner + MacVector + SnapGene)

Goal: prioritize missing capabilities that recur across multiple external
cloning tools, instead of chasing one-off parity points.

Reference matrix:

- `docs/cloning_tool_gap_matrix.md`

Repeated multi-tool gaps to prioritize:

1. Primer design and validation workflow:
   - baseline now implemented:
     - first-class `DesignPrimerPairs` operation
     - persisted report contract + shell/CLI inspect/export routes
     - optional Primer3 backend selection (`auto|internal|primer3`) with
       deterministic fallback and backend provenance in reports
   - next:
     - extend the shipped shared PCR-assay protocol-cartoon family on top of
       the existing template/binding renderer and shared figure building
       blocks:
       - batch ROI queue and nested-PCR extensions through bindings, repeated
         events/molecules, and shared blocks rather than new renderer
         semantics
       - additional explicit artifact lanes once the renderer grows beyond
         DNA-only rows
     - add nested-PCR primer design workflow contracts:
       - outer + inner primer-pair design in one deterministic request/report
       - explicit nesting constraints (inner amplicon must lie within outer
         amplicon with configurable margin)
       - adapter-equivalent export/inspection routes for outer/inner pair sets
     - add explicit long-range PCR workflow mode:
       - support primer-pair design where target amplicons can exceed one
         viewport span in sequence GUI (far-apart anchors/primers)
       - preserve deterministic reporting and GUI guidance for very long
         amplicons.
     - add chromosomal-translocation PCR workflow mode:
       - support primer-pair intent where forward and reverse primers target
         different chromosomes in the reference genome model
       - produce deterministic fusion-aware report payloads for breakpoint
         tests.
     - add deferred genome-translocation modeling track (depends on genome-view):
       - allow users to declare/simulate translocation events in-genome view.
       - define how synthetic fusion loci are represented for primer design.
       - evaluate BLAST indexing implications for fusion references and
         breakpoint-spanning primer specificity checks.
     - expand Primer3 constraint-mapping parity and fixture-backed equivalence
       coverage versus internal backend normalization
     - add richer Primer3 preflight diagnostics/UI surfacing
       (binary/version/config-path checks + environment guidance)
     - add hash-guided primer targeting mode:
       - use persisted RNA seed-support maps to bias candidate primer anchors
         toward transcript-coherent regions
       - penalize anchors dominated by common/reused hashes that show broad
         cross-gene ambiguity
       - keep this as deterministic guidance input; final pair scoring still
         runs through shared primer report contracts
     - pair interaction checks and richer ion-/structure-aware thermodynamic
       scoring beyond the current shared nearest-neighbor baseline
     - saved/reusable primer sets with explicit versioning
     - async-capable batch off-target/specificity checks so primer-pair
       selection can run multiple BLAST searches through agent/MCP/CLI routes
       with progress/cancel parity
2. Auto-annotation library scan:
   - detect missing/plausible features from curated vector/feature libraries
   - keep machine-readable confidence/overlap diagnostics for automation
   - legacy `gentle-m` confirms this is still valuable, but intake should be
     report-first/apply-second rather than direct annotation side effects
   - use `gentle-m` `src/AutoAnnotate.cpp` only as heuristic seed material; do
     not port dialog/database behavior directly
3. Sequencing confirmation workflow:
   - import read evidence (Sanger/NGS-aligned scope), map evidence to expected
     constructs, and emit deterministic pass/fail evidence summaries
   - detailed implementation plan is tracked in
     `docs/sequencing_confirmation_plan.md`
   - legacy `gentle-m` intake priority: do this before resurrecting broader
     legacy cloning conveniences because it most directly improves user trust
   - stage the work as:
     - called-sequence Sanger/amplicon confirmation first
     - ABI/AB1/SCF raw-trace import and inspection second
     - trace-aware confirmation integration third
   - likely legacy code seeds:
     - `gentle-m` `src/ABItype.cpp`
     - `gentle-m` `src/SCFtype.cpp`
     - `gentle-m` `src/TSequencingAssistantDialog.cpp`
   - related follow-up from the same legacy intake:
     sequencing-primer suggestion overlays should come immediately after the
     core confirmation report path, not as an isolated UI feature
4. Interactive cloning workspace model:
   - add explicit fragment/clipboard-style assembly workspace with end
     compatibility feedback, while execution still routes through shared ops
   - the legacy ligation dialog remains informative as algorithm/reference
     material (`gentle-m` `src/TLigationDialog.cpp`), but it is lower priority
     than typed confirmation, primer, and annotation workflows
5. CRISPR workflow closure:
   - extend existing guide-design baseline to include practical
     screening/confirmation workflow outputs
6. PCR assay modality expansion (PCRtools-derived parity set):
   - add inverse-PCR mode for circular templates with deterministic amplicon
     boundary semantics and report parity across adapters
   - add bisulfite-aware primer-design modes (PCR and optional LAMP branch)
     with explicit conversion policy fields in operation contracts
   - add allele-specific genotyping workflow contracts (KASP/PACE/ASQ) for
     SNP/InDel assay design and machine-readable assay-output bundles
   - add LAMP primer-set design workflow (`F3/B3/FIP/BIP` + optional `LF/LB`)
     with deterministic scoring/reporting and GUI assist
   - add multiplex tiling panel-design contracts (target ranges -> tiled primer
     pools) with deterministic pool-assignment/export payloads
   - add virtual-PCR/off-target route with mismatch-tolerance controls and
     structured multi-hit reporting for multiplex preflight

Notes:

- If visual comparisons include screenshot/raster baselines, those artifacts
  remain manual contributions while screenshot execution is policy-disabled.
- External-code intake decision (2026-03-20, paper-first intake for
  `10.1016/j.omtn.2025.102716`):
  - do not vendor/import upstream JavaScript directly into GENtle engine paths
  - rationale:
    - primary feature intake source is the peer-reviewed article discussion and
      methods/data-availability statements, with repository content treated as
      secondary context only
    - architecture mismatch (PCRtools places core assay logic in browser-page
      JavaScript, while GENtle requires shared engine-owned contracts)
    - maintainability/auditability risk from obfuscated/minified upstream
      scripts in core assay modules
    - GPLv3 upstream licensing is technically consumable via GENtle's `GPL-2+`
      policy, but would still require explicit maintainership/legal sign-off
      before any direct code derivation
  - preferred path:
    - re-spec required behavior in `docs/protocol.md`
    - reimplement in Rust engine with deterministic fixtures/tests
    - keep optional output-level comparison harnesses against PCRtools as
      non-shipping validation aids

### Primer3 wrapper integration track (new)

Goal: add Primer3 tooling support without fragmenting GENtle's shared engine
contracts or adapter parity guarantees.

Phase 1 (wrapper + normalization baseline): implemented baseline

- Engine-owned Primer3 adapter is now wired behind `DesignPrimerPairs`.
- Shared report schema remains unchanged (`gentle.primer_design_report.v1`) and
  deterministic tie-break ordering remains backend-independent.
- Reports now include backend provenance fields
  (`requested`, `used`, optional fallback reason, executable/version).
- Shell/CLI backend override controls are now available:
  - `primers design ... --backend auto|internal|primer3`
  - `primers design ... --primer3-exec PATH`
  - `set-param primer_design_backend ...`
  - `set-param primer3_executable ...`

Phase 2 (tooling diagnostics + compatibility hardening): in progress

- Shell/CLI report payload now captures executable/version diagnostics when
  Primer3 is used (or fallback metadata in auto mode).
- GUI Engine Ops now includes dedicated Primer3 preflight/status controls:
  - backend selector (`auto|internal|primer3`)
  - executable path field (`primer3_executable`)
  - explicit probe action with structured status text
- shared-shell non-mutating preflight command is now available:
  - `primers preflight [--backend auto|internal|primer3] [--primer3-exec PATH]`
- Fixture-backed adapter-equivalence coverage now includes deterministic
  internal-vs-Primer3 normalization parity tests (fake local Primer3 script +
  committed key/value fixture).
- Remaining:
  - additional config-path diagnostics for complex installations
- Keep failure modes deterministic and machine-readable
  (`Unsupported`/`Io`/`InvalidInput` with stable message contracts).

Phase 3 (async specificity tier + agent/MCP parity): baseline started

- Shared shell/CLI + MCP now expose async BLAST primitives:
  - `genomes/helpers blast-start|status|cancel|list`
  - MCP `blast_async_start|status|cancel|list`
- Async BLAST now runs through a bounded FIFO scheduler with explicit
  queue/running metadata in job status payloads.
- Deterministic parity test baseline now covers MCP-vs-shared-shell async
  BLAST status routing.
- Remaining:
  - primer-pair-specific multi-BLAST orchestration on top of these primitives
  - richer progress granularity and GUI binding for agent-triggered async jobs
  - broader cross-adapter integration tests for cancellation/progress semantics
  - optional shared-shell scheduler control command
    (`genomes/helpers blast-scheduler`) to inspect effective async BLAST
    concurrency and set per-process overrides without relying on environment
    variables

### Unified BLAST abstraction + primer UI track (new)

Goal: unify BLAST execution semantics for standalone BLAST, single-primer
design, and primer-pair design through one engine-owned service contract.

Detailed execution plan: `docs/primer_design_specialist_plan.md`.

- Primer specialist window + persisted alternative-filter views are tracked in
  `docs/primer_design_specialist_plan.md`.

Phase 1 (engine abstraction + options layering):

- Add shared BLAST service/request/response/progress abstractions in engine,
  reusable by BLAST and primer workflows.
- Implement deterministic option layering:
  - built-in defaults
  - optional defaults file
  - project-level override operation
  - per-request JSON overrides
- Enforce strict option validation (unknown key/type/range failures).
- Persist both raw override JSON and resolved effective options in provenance.

Status (2026-03-01):

- Implemented baseline in engine:
  - layered option resolution (`built-in -> defaults file -> project override -> quick flags -> request JSON`)
  - strict request/default/project JSON object validation
  - threshold filtering (`max_evalue`, `min_identity_percent`,
    `min_query_coverage_percent`, `min_alignment_length_bp`, `min_bit_score`,
    `unique_best_hit`)
  - provenance/report fields for raw request override + resolved effective options
  - project-level parameter operations:
    `blast_options_override`, `blast_options_defaults_path`
- Shared-shell/CLI BLAST routes now accept:
  - `--options-json JSON_OR_@FILE`
  - `--options-file PATH`
- JS/Lua BLAST wrappers now accept optional `options_json` argument and route
  through the same request-options engine path.
- Remaining in Phase 1:
  - GUI advanced options editor/preset UX (Phase 2 UI work item)
  - parity tests that exercise identical option layering across GUI/CLI/JS/Lua.

Phase 2 (BLAST UI alignment):

- Refactor BLAST specialist window into stable sections:
  input, target, options, execution, results.
- Add advanced options JSON editor + preset selector while keeping quick
  controls for common fields.
- Keep heartbeat/query-count progress and invocation visibility wired through
  the shared BLAST progress/result contracts.

Status (2026-03-01):

- Implemented GUI BLAST sectioning baseline (`Target/Input/Options/Execution/Results`).
- Implemented BLAST options controls in GUI:
  - quick controls (`task`, `max_hits`)
  - preset selector
  - structured threshold controls (toggle + typed fields)
  - advanced JSON editor + JSON file loader
  - effective-options preflight preview via shared engine resolver.
- Result panel now displays both request override JSON and resolved effective
  options payloads.
- Implemented BLAST cancellation from both:
  - BLAST dialog execution controls (`Cancel BLAST`)
  - Background Jobs panel (`BLAST -> Cancel`)
  - worker-side interruption now propagates cancellation through shared engine/genomes BLAST paths.
- Remaining in Phase 2:
  - preset catalog persistence and named reusable GUI presets
  - polish structured threshold UX (units/tooltips/preset interop helpers).

Phase 3 (single-primer design UI):

- Add dedicated sequence-context primer-design specialist window:
  - target region selector
  - primer constraints
  - specificity policy (backed by shared BLAST abstraction)
  - ranked primer table with diagnostics
- Add actions for annotation/export and forwarding selected primers into PCR
  workflows.

Phase 4 (primer-pair UI and specificity tier):

- Expand GUI around `DesignPrimerPairs` with explicit amplicon intent, pair
  constraints, and specificity tiers.
- Add nested-PCR mode in the same primer-pair specialist flow:
  - design/report outer and inner primer pairs together
  - expose nesting margin/amplicon-span constraints with deterministic
    validation
- Route pair-specificity BLAST fan-out through the same async BLAST service
  contract used by standalone BLAST.
- Keep report/provenance schema adapter-equivalent.

Phase 5 (adapter parity + tests):

- Add parity tests for request normalization, option resolution, progress, and
  provenance across GUI/CLI/JS/Lua paths.
- Add regression fixtures for primer/specificity edge cases and deterministic
  tie-break behavior in ranked outputs.

#### B) Agarose gel simulation improvements (high value)

Current state is strong for ladder-aware pool preview and export, but still
simplified relative to wet-lab interpretation.

Planned upgrades:

- Keep universal gel support for any container/tube (now implemented) and
  harden workflows:
  - preserve one-lane "single-product proof" use cases
  - preserve mixed cardinality lanes (single-sequence and pool containers)
  - improve arrangement edit/authoring UX for lane setup reuse
- Extend arrangement nodes in lineage/DALG:
  - arrangement is already modeled and rendered as an explicit node type that
    groups multiple input tubes under one experimental setup
  - current arrangement mode: `serial` (gel lanes)
  - next arrangement mode: `plate` (plate-reader/assay modeling)
- Enforce same-setup semantics for one serial arrangement:
  - all lanes in one virtual gel share one run configuration and one ladder
    context
  - no implicit "spliced gel from different runs" output in a single
    arrangement
- Multi-lane digest conditions with per-lane enzyme sets and batch
  "apply-to-all" behavior.
- Optional uncut/topology-aware migration model for circular DNA
  (supercoiled/nicked/linearized lane behavior).
- Band intensity based on estimated DNA mass per band, not only multiplicity.
- Co-migration grouping thresholds and explicit merged-band annotation.
- Lane-side fragment table with bp, estimated mass, source fragments, and
  cut-site context.
- Optional gel-conditions parameters (agarose %, buffer/model preset) with a
  deterministic default profile.

#### C) Sequence/map clarity benchmarking (high value)

Use a small fixed corpus of examples to compare clarity of GENtle maps and SVG
exports versus established readability patterns.

Planned upgrades:

- Build visual benchmark fixtures (sparse, dense annotations, dense RE sites,
  promoter/regulatory-heavy loci).
- If benchmark packs include screenshot/raster artifacts, capture and curation
  are manual contributions (agent screenshot route remains disabled by policy).
- Add explicit label-placement modes where useful (for example,
  "inside-preferred" vs "outside-preferred" behavior by feature class).
- Keep linear/circular/SVG parity for:
  - label anchoring
  - overlap resolution
  - lane separation
  - feature-to-label traceability
- Add quantitative readability checks in snapshot tests:
  - unresolved label overlaps
  - clipped labels
  - hidden high-priority labels (genes/CDS/ORFs/selected RE sites)
- Continue node/table/sequence-window identity consistency
  (`node_id <-> seq_id <-> open-window title`) for lineage traceability.

#### D) Dotplot + promoter-flexibility view track (new, high value)

Status (2026-03-19):

- Implemented baseline:
  - engine operations:
    - `ComputeDotplot`
    - `ComputeFlexibilityTrack`
  - persisted metadata schemas:
    - `gentle.dotplot_view.v2`
    - `gentle.flexibility_track.v1`
  - shared-shell/CLI commands:
    - `dotplot compute|list|show|render-svg`
    - `render-dotplot-svg`
    - `flex compute|list|show`
  - deterministic tests:
    - engine operation storage/retrieval
    - shell parse + shell execute paths
  - pairwise extension (2026-03-15):
    - `ComputeDotplot` supports `pair_forward` and
      `pair_reverse_complement` with explicit reference sequence/span fields
    - GUI `Dotplot map` supports query-vs-reference rendering with separate
      x/y spans and query-axis selection sync
    - shell/CLI `dotplot compute` supports `--reference-seq`, `--ref-start`,
      `--ref-end`
  - lineage integration (2026-03-19):
    - main lineage table/graph now materializes persisted dotplot/flex artifacts
      as analysis nodes
    - analysis nodes are linked to query/reference source sequences by operation
      edges and carry artifact id + mode/model + count metadata in tooltips/details
  - mapping-summary refinement (2026-03-19):
    - dotplot payload now stores per-query-bin boxplot summaries of reference-hit
      distributions (`boxplot_bin_count`, `boxplot_bins`)
    - standalone Dotplot workspace renders these boxplots below the density map
      for rapid exon-band distribution inspection
  - GUI export baseline (2026-03-19):
    - compact Dotplot map panel and standalone Dotplot workspace now expose
      `Export Dotplot SVG...`
    - Dotplot SVG export now routes through engine operation
      `RenderDotplotSvg` (GUI/CLI/shared-shell parity path)
    - main lineage table/details now expose `Dotplot SVG` action on dotplot
      analysis rows (analogous to gel export)
    - SVG export operations are projected into lineage table/graph as analysis
      nodes linked to source sequences
    - default export filename now encodes active dotplot parameters (`mode`,
      spans, `word`, `step`, `mismatches`, optional `tile`) and display controls
      (`threshold`, `gain`; plus flexibility parameters when enabled)
  - gene-extraction refinement (2026-03-19):
    - `ExtractGenomeGene` now auto-creates an exon-concatenated synthetic
      companion sequence (`<seq_id>__exons`) with deterministic `N` spacers
      between merged exon blocks for cleaner cDNA-vs-exon-only dotplot workflows
- Remaining:
  - additional overlay controls beyond crosshair baseline
  - additional adapter convenience wrappers as new dotplot operations are added

Latest GUI baseline (2026-03-09):

- Sequence windows now include a primary `Dotplot map` mode.
- Sequence-window `Dotplot map` now acts as a compact launcher/compute panel
  and opens a dedicated standalone `Dotplot` workspace window for full
  parameter editing and plot inspection.
- Dotplot compute defaults to a full-span view by setting `half_window_bp` to
  the larger query/reference sequence length; default `max_mismatches` is `0`
  (exact-seed baseline).
- Dotplot mode now supports:
  - parameterized compute (`mode`, `word_size`, `step`, `max_mismatches`,
    optional `tile_bp`)
  - persisted payload selection (`dotplot_id`)
  - optional paired flexibility-track panel (`flex_track_id`, model/bin/smooth)
  - linked crosshair baseline (hover + click-to-lock + sequence-selection sync)
  - density-render safeguards (point sampling cap for responsive redraws)
  - pairwise sparse-result diagnostics (orientation hinting, strict-parameter
    warnings, reference-edge warnings)
  - pairwise default auto-fit for reference span:
    - when `ref_start/ref_end` are empty, dotplot compute performs one
      full-reference pass, auto-fits to hit envelope (+padding), and recomputes
  - `Fit ref span to hits` manual helper remains available for explicit-span
    re-fit workflows
  - expanded engine pair-evaluation guardrail (now `100,000,000`) for larger
    pairwise spans before requiring coarser sampling
  - exact-seed acceleration for `mismatches=0` requests (indexed k-mer matching
    path avoids brute-force pair loops at low step sizes)

Goal:

- Add a low-latency dotplot workflow in DNA sequence windows and pair it with
  promoter-oriented flexibility tracks for interpretation.
- Keep this engine-owned so GUI/CLI/JS/Lua/SVG remain deterministic and
  equivalent.

Why now:

- Promoter-heavy analysis needs quick repeat/inverted-repeat inspection.
- Dotplot alone is not enough for flexibility claims; it should be paired with
  explicit score tracks.
- Existing background-job and SVG/export contracts are mature enough to host
  this without frontend-only logic.

Phase 1 (engine payload + compute baseline):

- Add typed payload schemas:
  - `gentle.dotplot_view.v2`
  - `gentle.flexibility_track.v1`
- Add compute operations:
  - `ComputeDotplot` (self-forward, self-revcomp, parameterized seed/mismatch)
  - `ComputeFlexibilityTrack` (model + binning + smoothing options)
- Store outputs in project metadata with deterministic cache keys:
  `(seq_id, span, params_hash, model_version)`.
- Add cooperative progress phases + cancellation support.
- Status: implemented baseline; progress-phase callbacks and advanced cache-key
  versioning are follow-up hardening.

Phase 2 (GUI sequence-window workflow):

- Add a dedicated `Dotplot` mode in sequence windows (primary view).
- Status: implemented baseline (bounded compute controls + persisted selection
  + in-view rendering + optional flexibility track panel + linked crosshair).
- Add optional low-alpha overlay mode (secondary, opt-in).
- Add linked coordinate crosshair between dotplot and sequence map.
- Add parameter presets:
  - `Fast preview`
  - `Balanced`
  - `Sensitive`
- Add flexibility-track toggles in the same control group, but keep visual
  channels separate from dotplot symbols to avoid ambiguity.

Phase 3 (CLI/JS/Lua parity + export):

- Add shared-shell commands:
  - `dotplot compute ...`
  - `dotplot show ...`
  - `dotplot render-svg ...`
  - `render-dotplot-svg ...`
  - `flex compute ...`
  - `flex show ...`
- Add adapter wrappers in JS/Lua/Python over the same operation payloads.
- Add engine export operation:
  - `RenderDotplotSvg` (optional flexibility-track panels on same coordinate axis).
- Status: partial.
  - implemented: `dotplot compute|list|show`, `flex compute|list|show`, GUI
    `Export Dotplot SVG...`, engine op `RenderDotplotSvg`, shared-shell/CLI
    `render-dotplot-svg`, `dotplot render-svg`, JS/Lua/Python
    `render_dotplot_svg(...)` convenience wrappers

Phase 4 (latency hardening):

- Use indexed seed matching (k-mer/minimizer style), not brute-force full
  matrix fill.
- Add multiresolution tile generation so initial draw is fast and refinement is
  viewport-driven.
- Reuse cached tiles/payloads between GUI redraw and SVG export.
- Add explicit performance budgets for common promoter spans.

Acceptance gates:

- Deterministic fixtures:
  - direct-repeat promoter fragment
  - inverted-repeat promoter fragment
  - low-complexity fragment.
- Snapshot tests for dotplot/flexibility SVG output.
- Parity tests proving GUI/CLI/JS/Lua receive equivalent payloads for the same
  parameter set.
- Cancellation + resume tests for long spans.
- Cache-key regression tests (parameter changes must invalidate/recompute;
  unchanged inputs must reuse cache).

Post-baseline follow-ups:

- Evaluate adding additional promoter mechanics models after baseline
  reproducibility is stable.

### UX declutter and readability improvements (planned)

- Add one-click map view presets (`Anchored`, `Cloning`, `Annotation`,
  `Signal`) that apply curated layer visibility bundles.
- Add zoom-aware rendering policy that suppresses tiny low-value glyphs at wide
  zoom and progressively reveals detail when zoomed in.
- Add a one-click declutter action that temporarily disables low-value overlays
  when feature overlap/noise is high.
- Add per-layer counts in visibility controls (for example `ORF (N)`) so users
  can predict visual noise before enabling a layer.
- Separate visual lanes more strictly (annotation rectangles vs predictive
  overlays vs signal tracks) to avoid overlap collisions.
- Expand window visual identity from experimental to production-ready:
  - add backdrop asset validation/file picker
  - add preview panel for each window type
  - enforce readability floor (contrast + max opacity caps)

### Post-delivery hardening backlog (shared)

- Finish-button discoverability hardening:
  - Run a completeness audit so every actionable button has both a concise
    tooltip and a stable hover-status name in the status bar.
  - Add regression checks for newly introduced dialogs/panels so these labels
    do not silently disappear.
- Dense-feature readability hardening:
  - Continue smarter lane packing with stronger non-overlap constraints and
    reduced wasted vertical gaps.
  - Keep strict `REG@DNA` mode enforcing DNA-level placement only (never
    intruding into gene lanes) and add focused visual regression fixtures.
- Undo/redo reliability hardening:
  - Expand operation-level history tests around mixed GUI workflows and
    background-job gating.
  - Keep visible history and transitions consistent after imports/retries.
- Guided anchored-import hardening:
  - Keep preflight summary deterministic for GenBank/BED/BigWig/VCF
    (anchor detection, match status, projected tracks, apply-to-all behavior).
  - Add broader fixture coverage for multi-anchor projects and re-anchoring.
- Command-palette hardening:
  - Keep action naming/keywords stable and improve ranking as command count
    grows.
  - Add keyboard-navigation regressions (`Cmd/Ctrl+K`, arrows, enter, escape).
- Background-jobs robustness:
  - Implemented now:
    - typed job-event records with `kind/phase/job_id/timestamp/summary`
    - monotonic per-job IDs for prepare/track-import/BLAST workers
    - stale worker-message rejection by `job_id`
    - idempotent cancellation handlers reused by dialogs and jobs panel
    - explicit retry events from jobs panel actions
    - recent job-event history now persists in project metadata
      (`gui.background_job_history`, schema
      `gentle.gui_background_job_history.v1`) and is restored on project load
      / reset paths.
    - retry actions from the Jobs panel now capture deterministic argument
      snapshots (per job family) and persist them in the same metadata payload
      for reproducibility/debugging across restarts.
    - background-jobs retry snapshot list now supports kind/text filtering and
      filtered JSON export for larger-history triage handoff.
    - background-jobs retry snapshot retention controls now support retain-count
      cap, explicit oldest-first prune, and clear-all cleanup (persisted via the
      same metadata payload).
    - background-jobs retry snapshot workflows now include filtered bulk-delete
      and archive-and-delete (JSON artifact + removal) for targeted per-filter
      cleanup.
    - destructive filtered cleanup now uses staged confirm UX with preview
      summaries (match counts/kinds/origins/id-range) before delete/archive.
    - staged destructive cleanup now includes a dry-run diff panel that previews
      "would remove" vs "would remain" snapshot rows before confirm.
    - destructive cleanup confirm now requires an action-specific type-to-confirm
      phrase before confirm is enabled.
    - successful destructive cleanup actions now append persisted cleanup-audit
      entries (action/filter/counts/archive path) rendered in the Jobs panel.
    - cleanup-audit history now has a dedicated JSON report export action; this
      export is intentionally read-only and does not self-append audit entries.
    - cleanup-audit history now supports action/text filtering and independent
      retention controls (`retain newest N`, prune oldest, clear-all) so
      long-running sessions can keep bounded, searchable audit trails.
    - cleanup-audit report export now respects current audit filters (action +
      text), and audit `clear all` now uses staged type-to-confirm before
      removal.
  - Next:
    - optionally add dual-mode cleanup-audit export (`filtered` vs `full`) and
      quick-filter chips for common audit actions

### Current branch blockers (must clear first)

- None currently blocking on this branch. Latest local run: `cargo test -q`
  passed (`528 passed, 1 ignored` in main suite; additional suites green).

### Stability TODO (queued: none; items 4, 5, 6, 7 done)

- Done (2026-03-04): 5) malformed-annotation reporting now summarizes non-fatal
  GTF/GFF parse issues with file/line context and surfaces the warnings through
  prepare/report payloads consumed by engine/CLI/GUI paths.

- Done (2026-03-19): 4) async long-running BLAST job durability hardened:
  - BLAST async job snapshots are now persisted in project metadata
    (`gentle.blast_async_job_store.v1`) with deterministic ordering and
    counter carry-forward,
  - restart/reload recovery now normalizes orphaned non-terminal jobs
    deterministically (`queued|running` -> `failed` with explicit interruption
    reason, or `cancelled` when cancellation was requested),
  - cancellation semantics after restart/reload are deterministic (cancel can
    be applied to recovered non-terminal snapshots before normalization),
  - conformance tests now cover restart recovery and stale-snapshot
    race-transition normalization behavior.

- Done (2026-03-19): 7) cancellation/timebox hardening extended to prepare jobs:
  - prepare worker progress now routes timeout enforcement through one source
    of truth (engine-level `timeout_seconds`) instead of duplicate app-side
    timebox checks,
  - prepare completion status now uses explicit cancellation-request +
    configured-timebox context to classify terminal outcomes
    (`cancelled` vs `timed out` vs `failed`) deterministically,
  - prepare progress/completion UI/status summaries now include
    cancellation-request context, including success-after-cancel-request cases,
  - conformance tests now cover cancel-request classification, timebox
    precedence over cancellation wording, and completion reporting after a
    cancel request.
- Done (2026-03-19): 6) external-binary preflight diagnostics added for BLAST
  and related tooling before long jobs start:
  - shared engine/genomes preflight schema
    `gentle.blast_external_binary_preflight.v1` now reports per-tool
    `found/missing`, `version`, configured executable token, and resolved path
    for `blastn` + `makeblastdb`,
  - shared-shell routes now include deterministic `binary_preflight` payloads
    for `prepare`, `blast`, `blast-track`, and async `blast-start`,
  - GUI prepare/BLAST background-job startup status now surfaces the same
    preflight diagnostics before long-running execution proceeds,
  - conformance tests now cover missing-tool and configured-tool preflight
    reporting plus shell-output payload coverage.

### Missing test coverage (current priority list)

- Done (2026-02-24): shared-shell execution tests now cover
  `resources sync-rebase` and `resources sync-jaspar` using local fixture
  inputs (no network dependency), including a focused assertion that motif
  reload side effects are applied after `resources sync-jaspar`.
- Done (2026-02-27): MCP UI-intent parity tests now compare MCP tool outputs
  against direct shared-shell `ui ...` outputs for discovery, prepared-query,
  latest-selection, and open-intent resolution.
- Done (2026-03-03): JS adapter parity tests now compare `import_pool(...)`,
  `sync_rebase_resource(...)`, and `sync_jaspar_resource(...)` wrapper outcomes
  against direct shared-shell execution.
- Done (2026-03-03): Lua adapter parity tests now compare `import_pool(...)`,
  `sync_rebase(...)`, and `sync_jaspar(...)` wrapper outcomes against direct
  shared-shell execution.
- Done (2026-03-03): CLI forwarded-dispatch parity tests now verify
  `import-pool` / `resources sync-rebase` / `resources sync-jaspar` top-level
  routes produce adapter-equivalent outputs and state deltas to direct shared
  shell execution.

## 3. Recommended execution order

### Phase A: AI communication + safety plane

- Expand MCP server from guarded op/workflow baseline to broader deterministic
  tool coverage over shared shell/engine routes.
- Keep handlers thin and adapter-equivalent (no MCP-only biology branches).
- Harden mutating-intent safety policy uniformly across agent, voice, and MCP
  invocation paths.
- Keep the implemented UI-intent tool routine stable and extend parity to
  additional routed command families.
- Keep `ui ...` intent routing deterministic and continue discoverability paths
  (menu, command palette, shell/agent/MCP intent surfaces).

### Phase B: cloning routine standardization

- Execute cloning-routine catalog phases 1-4:
  - catalog schema + routine indexing
  - typed macro input/output contracts
  - macro-run instance recording
  - macro box nodes in lineage/workflow graph
- Add protocol macro template packs for the cloning modes listed in Section 2
  (restriction + Gibson baseline already shipped; next focus on Golden Gate,
  Gateway, TOPO, TA/GC, In-Fusion, NEBuilder HiFi).
- Implement routine-application assistant surfaces on top of the same routine
  catalog/preflight contracts:
  - shared `routines explain` / `routines compare` payloads
  - GUI staged apply flow with explicit alternative disambiguation
  - optional decision-trace export in process protocol/run-bundle artifacts.
- Guardrail for future genomic-context routines:
  - do not assume "one focal locus only" if that would block later
    cohort-based comparative evidence workflows.
  - future routine work may consume focal-ROI / matched-locus-cohort /
    assay-evidence-comparison outputs, but should not invent a parallel
    similarity/cohort model in macros, catalog metadata, or adapter code.
  - any new routine touching genomic-context ranking should preserve the path
    toward focal ROI -> matched-locus cohort -> assay evidence comparison.
- Start repeated cross-tool cloning UX gaps after routine packs land:
  - primer design/validation workflow contracts
  - nested-PCR primer workflow contract + UI/report parity
  - Primer3 wrapper integration + backend-equivalence test matrix
  - PCRtools-derived assay-modality expansion:
    - inverse PCR on circular templates
    - bisulfite-aware PCR/LAMP design modes
    - KASP/PACE/ASQ allele-specific genotyping workflows
    - LAMP primer-set workflow contracts
    - multiplex tiling panel design + pool assignment exports
    - virtual-PCR/off-target reporting routes
  - auto-annotation library scan contracts
  - interactive cloning workspace/clipboard model

### Phase C: engine/protocol parity hardening

- Keep adapter-level helpers thin and aligned with engine operations.
- Promote remaining adapter-level utilities into first-class engine operations.
- Extend process-protocol export beyond the shipped run-bundle baseline:
  - workflow-scoped artifact packaging (optional input/output file copies and
    checksums)
  - stricter schema evolution and compatibility reporting for protocol bundles.
- Define shared frontend-neutral view-model schema for tracks/features/overlays.
- Complete XML follow-up (`INSDSet/INSDSeq`) without semantic divergence.
- Add sequencing-confirmation evidence contracts (read-aligned construct
  validation summaries) as a deterministic shared-engine path.
- Detailed sequencing-confirmation design note:
  - `docs/sequencing_confirmation_plan.md`
- Legacy `gentle-m` intake note:
  - priority/rationale/code-mining plan is tracked in
    `docs/legacy_gentle_m_intake.md`
  - best near-term legacy carry-forward targets are:
    sequencing confirmation, sequencing-primer overlays, silent-mutation
    verification, and auto-annotation

### Phase D: visualization and workflow UX

- Continue dense-case hardening for adaptive linear DNA letter routing:
  - visual benchmark fixtures and regression gates for crowded labels/features
  - snapshot-style stress coverage for condensed readability constraints
  - manual screenshot contribution for curated visual baselines where required
    (agent screenshot capture remains policy-disabled)
- Continue alternative-splicing follow-ups:
  - dense-fixture regression tests for boundary visibility and label safety
  - coordinate-true geometry invariants
  - primary splicing map-mode polishing (dense labels + optional interactions)
- Start GUI-driven feature editing track:
  - engine operation contract for add/update/delete feature edits
  - provenance-rich edit history entries and reversible snapshots
  - sequence/splicing-window editing UX (main window stays orchestration-only)
- Continue gel work:
  - arrangement authoring/editing UX (create/update/reorder lanes)
  - `plate` arrangement mode as first-class engine + adapter entity
  - one-run/one-setup semantics per serial arrangement
  - realism upgrades (topology, intensity, co-migration, lane tables)
- Add visual benchmark fixtures and readability regression gates for map
  export; treat screenshot/raster baseline assets as manual contributions.
- Add focused-region fallback + regression/snapshot coverage for unified
  scroll/zoom policy and close remaining feature-tree UI snapshot gaps.

### Phase E: integration polish and deferred policy items

- Add cross-application clipboard interoperability through versioned contracts.
- Continue OCI distribution follow-up for the shipped Debian-first container:
  - validate Linux/Apptainer pull path against the published OCI tag in CI or
    release smoke documentation
  - pin the non-Debian `rnapkin` install source/version in the Dockerfile once
    release policy is fixed
  - decide whether release attachments should also include an exported
    `docker-archive` or `.sif`, or whether OCI-only remains sufficient
  - keep a dedicated `.def` file deferred unless HPC-specific divergence
    appears
- Keep screenshot re-enable work as the final item and only after explicit
  endpoint-security exception/approval.
- Keep first-class protein-sequence support deferred:
  - UniProt remains metadata/projection input (fetch/import/list/show/map).
  - do not reopen protein sequence-window import until a dedicated protein
    window model exists (layout, feature toggles, exports, provenance UX).
  - when resumed, ship as a separate export/analysis channel rather than a
    DNA-window fallback.
- Keep auto-updated documentation with embedded graphics postponed until above
  safety/contract priorities are complete.
- Add deferred documentation track: **GUI screenshot atlas automation (blocked
  by stability gate)**.
  - Stability gate:
    - do not start screenshot-atlas implementation until GUI stability baseline
      is reached.
  - Staged deliverables after gate:
    1. re-enable `screenshot-window` in GUI shell context via shared viewport
       screenshot capture path.
    2. add `docs/screenshots/manifest.json` and generated
       `docs/gui_screenshots.md`.
    3. add `screenshot-atlas list` / `screenshot-atlas run` command routes.
    4. keep CLI/MCP raster screenshot capture disabled outside GUI context and
       use explicit SVG export commands for headless workflows.

### Phase F: interpretation (later)

- Image/sketch to `StatePatchProposal` translation with confidence scoring and
  explicit confirmation before apply.

## 4. Practical resume checklist

1. Read `docs/architecture.md`.
2. Read `docs/protocol.md`.
3. Read this file (`docs/roadmap.md`).
4. Run quick sanity:
   - `cargo check -q` (expected green)
   - `cargo run --bin gentle_cli -- capabilities` (run after `cargo check`
     is green)
   - `cargo run --bin gentle_cli -- state-summary` (run after `cargo check`
     is green)
5. Continue with highest-priority item from Section 2.

## 5. Speculative parking lot (not prioritized)

Purpose:

- capture plausible future ideas that are worth remembering,
- avoid losing them during focused delivery work,
- keep them clearly separate from committed execution order and current status.

Current parking-lot ideas:

- Cross-platform findings/artifact inspection for agent-driven work:
  - treat user-facing "findings" (interpretation notes, rationale summaries,
    evidence callouts, suggested next steps, and explanation artifacts) as
    engine-owned portable records rather than GUI-only session text.
  - likely placement if pursued:
    - compact structured summaries in `ProjectState.metadata`
    - optional larger sidecar/run-bundle artifacts for rich attachments
      (for example exports, screenshots, or reproducibility payloads)
    - lineage links back to the operation/report/sequence context that produced
      the finding
  - intended benefit:
    - the same finding can be inspected from GUI, CLI, MCP, JS/Lua, Python, or
      a future browser/WebAssembly frontend without adapter-specific re-entry or
      loss of provenance/explanation context.
  - note:
    - this is intentionally a parking-lot idea, not a scheduled roadmap item
      yet; if adopted, it should preserve the single-engine / machine-readable
      contract described in `docs/architecture.md`.
