# GENtle Changelog

This file is the auditable record of completed work that used to accumulate in
[`roadmap.md`](roadmap.md). Keep entries brief and outcome-focused. Durable
implementation rules belong in [`decisions.md`](decisions.md), not here.

Maintenance rule:

- Do not add completed-work bullets to [`roadmap.md`](roadmap.md).
- At the end of each session, move completed roadmap items here.
- Prefer a date and one short outcome sentence. Include command names,
  document names, schemas, or feature names only when they help a reader
  understand what changed.

## 2026-06-27

- Added fact-aware introspection for sequence-gated variant promoter/reporter
  routes and raw operation rows: `variant annotate-promoters`,
  `AnnotatePromoterWindows`, `variant promoter-context`,
  `SummarizeVariantPromoterContext`, `variant reporter-fragments`,
  `SuggestPromoterReporterFragments`, and `AnnotateTfbs`. Optional JSON output
  paths are modeled as `artifact.written` external handoffs rather than
  invented persistent project reports.
- Added closed-world `dotplot.exists` and `flexibility_track.exists` project
  facts plus fact-aware introspection for dotplot/flex compute, show, overlay,
  and SVG render routes. Deterministic `DOTPLOT_ID`/`TRACK_ID` bindings can now
  be verified after compute operations, while dotplot SVG export remains an
  external handoff.
- Added closed-world `candidate_set.exists` introspection for persisted
  candidate-window sets and fact-aware readiness/effect descriptors for
  candidate generation, show/metrics, score, filter, top-k, Pareto, set algebra,
  and delete routes. Delete rows are readiness-gated but do not yet claim a
  positive hard absence effect.
- Added closed-world `guide_set.exists`, `guide_filter_report.exists`, and
  `guide_oligo_set.exists` introspection for guide-design workflows, including
  fact-aware descriptors for guide-set upsert/show/filter/delete, oligo
  generation/list/show, and guide oligo/protocol exports.
- Added closed-world `container.exists`, `arrangement.exists`, and
  `rack.exists` introspection for persisted wet-lab container/rack authoring,
  including readiness/effect descriptors for arrangement creation, rack
  placement/mutation, and rack SVG/OpenSCAD/simulation exports.
- Added closed-world `workflow_macro_template.exists`,
  `candidate_macro_template.exists`, and `macro_instance.exists`
  introspection for persisted macro templates and recorded macro lineage rows,
  including readiness/effect descriptors for template show/upsert/delete/run
  and macro-instance inspection.
- Added fact-aware external sequence creation introspection for `LoadFile`,
  `genbank fetch`, `FetchGenBankAccession`, `ensembl-region fetch`,
  `FetchEnsemblRegion`, `dbsnp fetch`, `FetchDbSnpRegion`,
  `FetchUniprotLinkedGenBank`, `ImportUniprotEntrySequence`, and Ensembl
  gene/protein entry sequence import routes. These routes have no project-state
  preconditions, and post-run verification checks `sequence.exists(OUTPUT_ID)`
  when a deterministic id was supplied.
- Added fact-aware introspection for raw core sequence operation rows:
  `SaveFile`, `Digest`, `Pcr`, `PcrAdvanced`, `PcrMutagenesis`, and
  `PcrOverlapExtensionMutagenesis`. PCR rows require a loaded template and can
  verify deterministic product ids; digest and overlap-extension rows are
  readiness-only until prefix/rank-derived products are projected as facts.
- Added closed-world `uniprot_entry.exists`,
  `ensembl_gene_entry.exists`, and `ensembl_protein_entry.exists`
  introspection for stored protein/gene metadata entries. UniProt and Ensembl
  fetch/import routes can verify explicit entry ids, show/import-sequence
  routes require stored metadata facts, and metadata-backed sequence imports can
  still verify deterministic `sequence.exists(OUTPUT_ID)` products.

## 2026-06-26

- Added the first shared `introspect facts|capabilities|readiness|all`
  implementation, including additive fact domains, explicit headless
  `ui.host_available=false`, deterministic agent-system
  `host.tool_available` facts, registry-backed capability discovery,
  fact-annotated readiness descriptors, and `--seq-id` / `--readiness`
  introspection scoping.
- Added `introspect verify-effects` to verify fact-annotated
  `must_on_success` effects against the current fact graph plus supplied
  evidence, starting with restriction-scan report verification.
- Added fact-aware introspection for `sequence create`, including readiness
  with no fact preconditions and post-run verification that the requested
  output sequence id exists.
- Added fact-aware introspection for `variant materialize-allele`, including
  input sequence readiness and post-run verification of the materialized output
  sequence id.
- Added fact-aware readiness introspection for `align compute`, checking that
  both loaded query and target sequence ids exist before alignment.
- Added matching fact-aware introspection for the raw `AlignSequences`
  operation row so registry-driven adapters can make the same query/target
  readiness check as the shell command.
- Added fact-aware readiness introspection for `render-svg`, `render-rna-svg`,
  and `render-lineage-svg`, checking sequence inputs where applicable and
  modeling external SVG files as `artifact.written` `external_handoff` effects.
- Added matching fact-aware introspection for the raw render operation rows
  `RenderSequenceSvg`, `RenderRnaStructureSvg`, `RenderFeatureExpertSvg`,
  `RenderTfbsScoreTrackCorrelationSvg`, and `RenderLineageSvg`, keeping
  registry-driven adapters aligned with the shell render commands.
- Added fact-aware introspection for unambiguous raw persisted-report operation
  rows: `ListSequencingConfirmationReports`, `ListCutRunReadReports`,
  `ListRnaReadReports`, `ExportSequencingConfirmationReport`,
  `ExportSequencingConfirmationSupportTsv`, `ExportCutRunReadCoverage`, and
  `ExportRnaReadReport`.
- Added fact-aware introspection for raw persisted-report show operations:
  `ShowSequencingConfirmationReport`, `ShowCutRunReadReport`, and
  `ShowRnaReadReport`.
- Updated `introspect readiness` to evaluate full descriptor
  `precondition_expr` trees, including `any` branches, and added fact-aware
  introspection for the shared raw `ExportPrimerDesignReport` operation across
  primer and qPCR design reports.
- Added fact-aware introspection for raw primer/qPCR design operations:
  `DesignPrimerPairs`, `DesignInsertionPrimerPairs`, and `DesignQpcrAssays`.
- Added fact-aware introspection for raw `ReverseTranslateProteinSequence`,
  including protein-kind input readiness and deterministic output sequence
  effect verification.
- Added fact-aware introspection for raw `MaterializeVariantAllele`, including
  input-sequence readiness and deterministic output sequence effect
  verification.
- Added fact-aware readiness introspection for the read-only `rna-info`
  sequence inspector.
- Added fact-aware introspection for `reverse-translate run`, including
  protein-kind input readiness and post-run verification of the requested
  coding-DNA output sequence id.
- Added `report.exists` fact projection for persisted reverse-translation
  reports and fact-aware readiness introspection for
  `reverse-translate show-report`.
- Added `report.exists` fact projection for persisted primer/qPCR design
  reports and fact-aware introspection for `primers design` and
  `primers design-qpcr`, including template-sequence readiness and post-run
  verification of the requested report id.
- Added `report.exists` fact projection and fact-aware introspection for
  `primers prepare-restriction-cloning`, checking template, primer-report, and
  destination-vector readiness and verifying the generated handoff report.
- Added fact-aware readiness introspection for persisted primer, qPCR, and
  restriction-cloning report inspection commands:
  `primers show-report`, `primers show-qpcr-report`, and
  `primers show-restriction-cloning-handoff`.
- Added fact-aware, no-precondition catalog readiness for report-list commands:
  `primers list-reports`, `primers list-qpcr-reports`,
  `primers list-restriction-cloning-handoffs`, and
  `reverse-translate list-reports`.
- Added fact-aware readiness introspection for `features tfbs-summary`, treating
  it as a read-only sequence inspector whose readiness depends on the inspected
  `sequence.exists` fact.
- Added matching fact-aware introspection for the raw `SummarizeTfbsRegion`
  operation row with the same sequence-readiness model.
- Added fact-aware introspection for the raw
  `QueryProteinResidueGenomicCoordinates` operation row as a read-only
  sequence inspector.
- Added fact-aware catalog readiness for self-description commands
  `capabilities` and `state-summary`, both with empty project preconditions.
- Added fact-aware readiness introspection for `features query` and
  `features export-bed`; feature queries are read-only sequence inspectors,
  while BED exports model the output path as an `artifact.written` external
  handoff.
- Added matching fact-aware introspection for the raw `ExportFeaturesBed`
  operation row, keeping registry-driven adapters aligned with the shell BED
  export command.
- Added fact-aware introspection for raw sequence-context operations:
  `InspectSequenceContextView` checks sequence readiness, and
  `ExportSequenceContextBundle` additionally models its output directory as an
  external artifact handoff.
- Added fact-aware readiness introspection for `inspect-feature-expert` and
  `render-feature-expert-svg`; both require the inspected sequence to exist,
  and the SVG renderer models its output path as an `artifact.written`
  external handoff.
- Registered `view.selection` and `view.visible_tracks` fact vocabulary entries
  and added fact-aware catalog readiness for the `display` visibility command
  with a `view.visible_tracks` view-session effect.
- Added fact-aware catalog readiness for `ui intents`, exposing GUI intent
  discovery as a no-precondition view-intent route.
- Added fact-aware catalog readiness for `help`, exposing shell help/topic
  lookup as a no-precondition command-catalog route.
- Added fact-aware catalog readiness for `history status`, exposing
  undo/redo availability as a no-precondition status route.
- Registered `facts graph` and `facts eval` in the shared glossary and added
  fact-aware catalog readiness for both fact-layer routes.
- Added fact-aware readiness introspection for `reverse-translate export-report`,
  checking the source report and modeling the output JSON path as an
  `artifact.written` external handoff.
- Added fact-aware readiness introspection for `primers export-report` and
  `primers export-qpcr-report`, using persisted report facts and external JSON
  handoff effects.
- Added fact-aware readiness introspection for
  `primers export-restriction-cloning-handoff`, using the persisted handoff
  report fact and an external JSON handoff effect.
- Added `report.exists` fact projection for persisted sequencing-confirmation
  reports and fact-aware readiness introspection for `seq-confirm list-reports`,
  `seq-confirm show-report`, `seq-confirm export-report`, and
  `seq-confirm export-support-tsv`.
- Added fact-aware, no-precondition catalog readiness for
  `cutrun list-read-reports` and `rna-reads list-reports`.
- Added `report.exists` fact projection for persisted CUT&RUN read and RNA-read
  interpretation reports, plus fact-aware readiness introspection for
  `cutrun show-read-report`, `cutrun export-coverage`,
  `rna-reads show-report`, and `rna-reads export-report`.
- Added fact-aware introspection for raw RNA-read gene-support operations:
  `SummarizeRnaReadGeneSupport` and `InspectRnaReadGeneSupport` require an
  existing `rna_read` report and treat optional JSON output paths as external
  artifact handoffs.
- Added fact-aware introspection for no-project local catalog/report operations
  `SummarizeJasparEntries`, `BenchmarkJasparRegistry`, `ListJasparCatalog`,
  `ResolveTfQueries`, `ListReporterCatalog`, and `RecommendReporters`, treating
  optional JSON output paths as external artifact handoffs.
- Added `config.param` fact projection for engine-owned configuration
  parameters and fact-aware `set-param` introspection, including JSON-value
  effect verification through `introspect verify-effects`.
- Added fact-aware introspection for direct sequence-derivation operation rows
  `Reverse`, `Complement`, `ReverseComplement`, `Branch`, and `ExtractRegion`,
  checking the input sequence and verifying deterministic output sequence ids.
- Added fact-aware introspection for the raw `SetTopology` operation row,
  checking the target sequence and verifying the projected `sequence.circular`
  boolean fact.
- Added readiness-only fact-aware introspection for the raw `RecomputeFeatures`
  operation row, checking that the target sequence exists without inventing a
  computed-feature freshness fact.
- Added the closed-world `view.viewport` fact projection and fact-aware
  introspection for the raw `SetLinearViewport` operation row, including nested
  argument binding for exact viewport verification.
- Added closed-world `view.visible_tracks` fact projection from persisted
  display-layer booleans so words-only clients can read back visible track
  state, not only see display intent descriptors.
- Added fact-aware introspection for raw `SetDisplayVisibility` and
  `SetParameter` operation rows, mirroring the existing `display` and
  `set-param` shell descriptors without introducing duplicate fact names.
- Added fact-aware introspection for specialized RNA-read artifact exports
  such as hit FASTA, target-quality, exon-path, exon-abundance, score-density,
  alignment TSV, isoform-triage, and alignment-dotplot exports, including raw
  engine operation rows.
- Added fact-aware introspection for RNA-read alignment, isoform preflight, and
  hit-sequence materialization routes, with conservative effect declarations
  for selection-dependent materialized sequence ids.
- Added fact-aware introspection for built-in DNA/RNA ladder catalog list and
  export routes, including JS/Lua helper names and raw ladder export operation
  rows.
- Added fact-aware introspection for agent-system catalog list routes
  `agents list`, `agent_systems`, and `list_agent_systems`, keeping configured
  system enumeration distinct from live adapter/model-discovery readiness.
- Added fact-aware introspection for the raw `agent_preflight` MCP/tool
  operation row, mirroring `agents preflight` through the shared
  `host.tool_available(SYSTEM_ID)` readiness fact.
- Added fact-aware introspection for protocol-cartoon catalog, render, template
  validation, and template export routes, including the raw engine operation
  rows and external SVG/JSON artifact handoff effects.
- Added fact-aware introspection for no-project external inspection routes:
  prepared-cache inspection, CUT&RUN dataset catalog/status inspection, and
  array helper inspection/rendering, with SVG artifact handoff modeling where
  applicable.
- Added fact-aware readiness descriptors for genome-track imports, BLAST-track
  imports, and array projection routes. These rows require the loaded target
  sequence and intentionally remain readiness-only until feature
  freshness/track-update facts are projected.
- Added fact-aware introspection for no-project catalog/list routes covering
  candidate sets, candidate macro templates, guide sets, workflow macro
  instances/templates, and routine catalog list/explain/compare operations.
- Added fact-aware introspection for construct-reasoning graph list routes,
  treating optional sequence ids as filters rather than readiness
  preconditions.
- Added closed-world `construct_reasoning_graph.exists` introspection for
  persisted construct-reasoning graphs, including readiness/effect descriptors
  for named graph inspection, inspection-action listing/running, annotation
  status/writeback routes, and graph JSON export.
- Added fact-aware introspection for persisted dotplot and flexibility-track
  list routes, treating optional sequence ids as filters rather than readiness
  preconditions.
- Added fact-aware introspection for local metadata/catalog routes covering
  reference/helper catalog listing, stored Ensembl gene/protein metadata lists,
  and gene-group list/show/resolve/doctor routes with optional JSON artifact
  handoffs.
- Added fact-aware introspection for shell-level resource/catalog inspection
  routes covering JASPAR summary/list/inspect/TF-query resolution, resource
  status, UCSC rmsk indexing suggestions, publication dataset list/status,
  reference/helper catalog validation, and helper vocabulary list/doctor
  commands.
- Added fact-aware introspection for reporter catalog shell routes
  `reporters list`, `reporters recommend`, `reporters export-corpus`, and raw
  `ExportReporterCorpus`, keeping shell and operation-level catalog readiness
  aligned.
- Added fact-aware introspection for service readiness/provider catalog routes
  `services status`, `services providers list`, and
  `services providers doctor`, with optional doctor JSON output modeled as an
  external artifact handoff.
- Added fact-aware introspection for planning read-back routes
  `planning profile show`, `planning objective show`, and
  `planning suggestions list`, keeping mutating planning updates separate.
- Added fact-aware introspection for host/helper/protease/microRNA catalog
  helper routes, including JS/Lua/MCP helper catalog row names and optional
  protease JSON artifact handoffs.
- Added fact-aware introspection for reference/helper genome cache inspection
  and adapter readbacks (`genomes status|genes`, `helpers status|genes`,
  `list_reference_genomes`, `list_reference_catalog_entries`,
  `is_reference_genome_prepared`, `list_reference_genome_genes`), keeping
  concrete catalog/cache validation as execution-time behavior.
- Added fact-aware descriptors for the top-level read-only introspection shell
  routes: `introspect facts`, `introspect capabilities`,
  `introspect readiness`, `introspect verify-effects`, and `introspect all`.
- Added fact-aware introspection for async BLAST status/list readbacks:
  `genomes blast-status`, `helpers blast-status`, `genomes blast-list`,
  `helpers blast-list`, `blast_async_status`, and `blast_async_list`, while
  leaving start/cancel/execution rows registry-only until their external
  runtime effects are modeled.
- Added fact-aware introspection for agent model-discovery routes
  `agents discover-models` and `agent_models`, using the same
  `host.tool_available(SYSTEM_ID)` readiness model as agent preflight while
  declaring no project effects for model enumeration.
- Added fact-aware introspection for adapter parity aliases:
  `state_summary`, `reference_catalog_entries`, `ui_intents`,
  `ui_prepared_genomes`, and `ui_latest_prepared`.
- Added fact-aware introspection for generic GUI intent requests:
  `ui open`, `ui focus`, `ui close`, and `ui_intent`, with readiness gated by
  `ui.host_available`.
- Added fact-aware introspection for sequence-scan report/render routes:
  `FindRestrictionSites`, `features tfbs-score-tracks-svg`,
  `RenderTfbsScoreTracksSvg`, `SummarizeTfbsScoreTracks`,
  `features tfbs-track-similarity`, and `SummarizeTfbsTrackSimilarity`.
- Added fact-aware introspection for local external-service handoff routes:
  `services delivery-route`, `services project-preflight`,
  `services project-quote`, `services handoff`, and `services guide`, while
  keeping `services route-project-source` separate until conditional
  project-object preconditions are modeled.
- Added closed-world `sequencing_trace.exists` introspection for imported
  sequencing-trace evidence records, including readiness/effect descriptors for
  trace import/list/show shell routes and raw operation rows.
- Added fact-aware introspection for protease catalog/digest and protein-gel
  rendering routes. Protease digest readiness now requires an existing
  protein-kind sequence, persisted protein-derivation reports project as
  `report.exists == protein_derivation`, and SVG render routes model their
  output paths as `artifact.written` handoffs.
- Added fact-aware introspection for non-mutating primer helper readbacks:
  Primer3/backend preflight, feature/splicing ROI seed helpers, restriction
  cloning vector suggestions, and restriction-cloning handoff request seeding.

## 2026-06-24

- Added shared Agent Assistant/shell controls for DNA sequence-window selection
  (`ui selection sequence-window ...`) and project display-layer visibility
  (`display show|hide|visibility`), keeping window/selection actions separate
  from sequence deletion.
- Added Promoter design GUI parity for offline ortholog promoter cohorts and
  comparisons, plus ortholog relationship expectations that emit non-blocking
  TFBS/CUT&RUN divergence or concordance flags.

## 2026-06-22

- Added retrieval-producer metadata to `gentle.gene_set_resolution.v1` and the
  offline `gene-sets produce direct-list` route for local JSON/TSV gene-list
  caches, resolving candidates through the existing explicit-member resolver.
- Added the offline `gene-sets produce ontology-assignment` route for local
  ontology assignment JSON/TSV caches, keeping provider term membership distinct
  from local gene-group `external_mapping` resolution.
- Added the offline `gene-sets produce co-regulated` route for local
  evidence-derived cohort caches with explicit dataset/contrast, score
  threshold, sign-direction, and relationship-expectation metadata.
- Made produced gene-set resolutions, promoter cohorts, and gene-set CUT&RUN
  support reports persist as logical lineage artifacts rendered as `GeneSet`
  nodes, linked from producer operation to downstream analysis.

## 2026-06-21

- Added shared operand metavariable conventions for glossary `usage` rows and
  updated Agent Assistant prompt guidance so inner helpers consult GENtle docs
  before proposing commands and ask instead of guessing ambiguous operands.
- Added JS and Lua adapter bindings for construct-reasoning inspection actions,
  so scripts can list and run the same portable dotplot recommendations exposed
  through CLI, ClawBio, and MCP shell routes.
- Made construct-reasoning repeat/similarity task severity quantitative and
  objective-specific with protocol `score` fields, explicit score-to-bucket
  thresholds, and visible objective boost/down-weight rationales.

## 2026-06-19

- Deepened construct-reasoning repeat-family interpretation with a shared
  Alu/SINE/LINE/LTR/satellite taxonomy, stricter internal-vs-curated
  corroboration, and confidence tiers for repeat-family provenance.
- Added rule-based `task_severities[]` to construct-reasoning facts so
  repeat/similarity warnings can report PCR, nanopore, read-mapping, cloning
  stability, and construct-maintenance severity without creating extra map
  overlays.
- Completed the Clariom/probe-region real-data-to-figure bridge slice with
  coordinate-consistent `gentle.probe_region_evidence_interpretation.v2`,
  gated `arrays run-probe-region-backend`, GUI shared-capability surfacing, and
  runbook/docs coverage for the explicit E-MTAB-14704 TP73 loop.

## 2026-06-18

- Integrated materialized RepeatMasker/UCSC `rmsk`-style repeat annotations
  into construct reasoning so overlapping curated Alu/SINE repeat-family rows
  back soft internal repeat/mobile-element calls without duplicating fact rows.
- Made native HTTP agent transports tolerant of local models that omit or
  mis-shape the `schema` field when the returned JSON otherwise matches
  `gentle.agent_response.v1`; external stdio adapters remain strict.
- Added an explicit Msty MLX OpenAI-compatible agent template for
  `http://localhost:11973/v1`, with GUI/CLI docs that distinguish it from the
  `11964` Msty gateway when that gateway reports no model ids.
- Added construct-reasoning inspection-action rationale to the portable action
  payload and surfaced action mode, focus, evidence ids, and rationale in the
  existing GUI inspector rows.
- Added offline-first ortholog promoter reasoning with local
  `gentle.ortholog_resource.v1` mapping resources,
  `ResolveOrthologPromoterCohort`, `SummarizeOrthologPromoterComparison`, and
  `orthologs resolve-promoter-cohort` / `orthologs promoter-comparison` shell
  routes, including species/genome-matched CUT&RUN motif-support states.
- Added additive `relationship_flags[]` to
  `gentle.promoter_cohort_comparison.v1`, surfacing unexpected TFBS-track
  divergence for declared co-regulated promoter cohorts and unexpected
  concordance for declared anti-co-regulated cohorts as non-blocking review
  cues.

## 2026-06-17

- Documented sequence collection subjects as the shared model for gene sets,
  pools, arrangements, alignments, derived collections, and storage projections,
  with a GUI implementation plan for adapter-parity operation lifting.
- Completed gene-set cohort relationship support by deriving evaluated-only
  CUT&RUN occupancy consistency flags for declared co-regulated and
  anti-co-regulated promoter cohorts.
- Extended `planning protein-expression-handoff --seq-id` with read-only
  sequence/readiness/CDS/tag context so high-yield protein-expression handoffs
  report product uncertainty before review-gated provider or construct work.
- Added optional gene-set promoter-cohort relationship expectations
  (`manual`, `co_regulated`, `anti_co_regulated`) to promoter cohort and
  CUT&RUN aggregate reports, with non-blocking expectation flags and
  `--relationship` shell parsing.
- Routed protein-expression `suggested_next_actions[]` from
  `product_definition.readiness.status`, distinguishing CDS candidates,
  protein-only targets, and sequences that need CDS/ORF boundary review.
- Added Promoter design GUI parity for promoter expression evidence and
  cached CUT&RUN regulatory-support reports, additive four-state TFBS support
  status, and promoter cohort comparison; `genomes promoter-cohort-comparison`
  emits `gentle.promoter_cohort_comparison.v1`.
- Added portable construct-reasoning inspection actions for repeat/similarity
  dotplot handoffs, carrying deterministic action ids and driving evidence ids
  in the graph payload instead of deriving GUI buttons from labels.
- Added `construct-reasoning list-inspection-actions` and
  `run-inspection-action` so CLI/ClawBio users can list the same portable
  repeat/similarity dotplot recommendations and compute/render the selected
  action through shared dotplot operations; listing now supports fact,
  annotation/candidate, evidence, sequence, action-kind, and summary filters,
  and ClawBio exposes typed request modes
  `construct-reasoning-list-inspections` and
  `construct-reasoning-run-inspection` over the same shell contract.
- Added optional true PM probe-intensity matrix input to
  `arrays import-apt-probe-region-output` and the Clariom D GUI panel, with
  probe-level sample, condition-summary, and logFC columns preserved in
  `probe_intensity_chrom_order.csv`.
- Added `--level pm_probe` projection for probe-region helper outputs, so true
  PM probe rows marked `probe_level_input` can become genome-anchored array
  features without promoting parent-summary fallback rows.
- Added `arrays interpret-probe-region-evidence` /
  `InterpretProbeRegionEvidence` to compare projected probe/probeset-region
  array features with transcript/exon geometry while preserving shared
  transcript, multi-hit, and coordinate-projection ambiguity.
- Added first-class sequence-window controls to run and preview
  `InterpretProbeRegionEvidence` reports from the Clariom D / probe-region
  evidence panel, including bounded evidence and transcript geometry tables.
- Extended `InterpretProbeRegionEvidence` with explicit per-evidence
  transcript mappings that record exon ordinals, exon ranges, junction spans,
  and overlap base counts without turning array evidence into isoform calls.
- Added conservative geometry scores and score-basis guardrails to
  `InterpretProbeRegionEvidence`, with transcript-level unique/shared/
  constraining score sums for review-only probeset evidence triage.
- Added transcript-level `review_status` labels to probe-region interpretation
  reports so GUI/CLI users can triage unique, shared, constraining, and absent
  geometry without treating the result as an isoform call.

## 2026-06-16

- Added `services route-project-source` plus External Services GUI helpers so
  selected sequences/spans, persisted oligo order forms, and primer report pair
  ranks can produce delivery-route candidates while preserving duplicate-review
  gates before quote handoff.
- Added first-class sequence-window GUI controls for importing explicit APT
  probe-region summaries plus annotation tables into inspectable GENtle helper
  output, reusing the shared shell import/inspect/export/project contracts.
- Added Phase B oligo order handoff routes: `primers oligo-order route` and
  `quote` now reuse external-service delivery/quote contracts, block unreviewed
  duplicate forms before quote handoff, and preserve oligo line-item
  provenance/modification fields in normalized JSON/CSV bundles.

## 2026-06-15

- Added Phase A first-class oligo order forms under `primers oligo-order`,
  with deterministic line ids, primer/qPCR report provenance, duplicate/reuse
  grouping, review marking, list/show/export, and persistence inside the
  existing primer-design report store.
- Added the command-palette `Evidence Preparation` assistant for the TP73
  evidence-viewer proof path, reusing shared operations for local repeat,
  array, BED, TFBS, and proof-export preparation while keeping CEL/R/vendor
  steps as explicit copy-command handoffs.
- Added `services delivery-route` and the
  `gentle.external_service_delivery_route.v1` contract so generic "deliver this
  sequence" wording is classified by sequence kind before selecting Metabion or
  GeneArt quote-handoff routes.
- Added ClawBio external-service intents and request examples for provider
  catalog/doctor checks plus review-only Metabion oligo/m-block and GeneArt
  cloned-gene/protein-expression preflight and quote handoff routes.
- Routed high-yield protein-expression ClawBio intents such as "maximal amount
  of protein" to the read-only `planning protein-expression-handoff` request
  example and added the same scenario to the experimental follow-up catalog,
  graph, and human guide, with a review-gated GeneArt protein-expression quote
  packet as the downstream provider handoff.
- Added an explicit `services project-quote @...` suggested next action to the
  protein-expression handoff report so GeneArt preflight and quote-packet
  preparation remain separate review stages.
- Added a native GUI inspector for completed probe-region helper output, plus
  coordinate/build provenance, bounded preview rows, and projection-readiness
  blockers in `gentle.probe_region_output_inspection.v1`.
- Added `arrays render-probe-region-output-svg` and a matching GUI export
  control for deterministic native SVG plots of inspected `mean_log2_*` and
  `log2FC_*` probe-region helper tracks.
- Added `arrays project-probe-region-output` / `ProjectProbeRegionOutput` for
  direct-coordinate-compatible projection of inspected helper `log2FC_*` rows
  into genome-anchored array features.
- Added first-class sequence-window controls for projecting inspected
  probe-region helper output into genome-anchored array features without
  leaving the Clariom D evidence panel.
- Added explicit `coordinate_projections[]` / `projection_maps[]` support for
  probe-region helper-output projection, reusing the existing interval-map
  contract to preserve native and displayed array coordinates.
- Added explicit Affymetrix Power Tools command planning to
  `arrays probe-regions` when user-supplied PGF/CLF and optional MPS library
  files are present.
- Added `arrays import-apt-probe-region-output` to convert explicit APT summary
  output plus an explicit annotation/NetAffx coordinate table into GENtle's
  probe-region helper-output directory contract.
- Extended explicit APT probe-region imports with optional sample metadata
  matching that writes `mean_log2_*`, `sd_log2_*`, and default `log2FC_*`
  tracks for native inspection, plotting, and projection.
- Added optional `probe_intensity_chrom_order.csv` output for explicit APT
  imports when annotation rows provide PM probe coordinates, marking values as
  parent probeset-summary intensities when no explicit PM probe matrix is
  supplied.

## 2026-06-14

- Added the read-only `planning protein-expression-handoff` route emitting
  `gentle.protein_expression_handoff.v1` with product context, chassis/route
  candidates, high-yield missing questions, and a GeneArt protein-expression
  preflight scaffold.
- Added tutorial `09.03` for using the protein-expression handoff route from
  CLI/agent workflows, plus a ClawBio example request.

## 2026-06-13

- Added `arrays inspect-probe-region-output` to validate completed
  `probe_regions_oligo.R` output directories and summarize region-table,
  sample, condition, logFC, chromosome/gene, manifest, and provenance status as
  `gentle.probe_region_output_inspection.v1`.

## 2026-06-12

- Added `arrays probe-regions` as a shared-shell, GUI-shell-visible
  Affymetrix CEL probe/probeset-region preflight, emitting
  `gentle.probe_region_plan.v1` with CEL, metadata, annotation/library,
  platform, dependency, backend-candidate, metadata/contrast, output, and
  cache-readiness checks.
- Added `scripts/probe_regions_oligo.R` as a generic explicit R/oligo helper
  for the `arrays probe-regions` plan, producing chromosome-ordered intensity
  CSVs, expression/feature TSVs, limma contrast tables, provenance, and an
  RMA-normalized matrix manifest from user-supplied CEL files.

## 2026-06-06

- Added `protein_expression_max_yield` as a planning consultation
  `biological_intent`, with explicit high-yield protein-expression questions
  before any construct/vector route is accepted.

## 2026-06-05

- Switched push/PR CI to one commit-sampled platform per run, with manual
  platform override, a PowerShell-native Windows validation pass matching the
  Unix test shape, clearer `macOS` naming, a single Linux runner when Linux is
  selected, and no multi-GB `target` cache restores.
- Limited the container-image workflow to release tags and manual dispatch so
  ordinary pushes and PRs no longer build Docker CLI/GUI images.
- Added explicit `biological_intent` fields to reporter recommendation and
  reporter construct handoff reports, and staged generated synthetic-biology
  chef background assets for a future synthetic-biology window.

## 2026-06-04

- Wired the ClawBio `gentle-cloning` skill toward the existing reporter
  catalog, recommender, corpus export, and reporter construct handoff routes
  so synthetic-biology follow-up quotes GENtle reports instead of improvising.

## 2026-05-29

- Added opt-in GUI frame profiling behind `--features gui-profiler` and
  `GENTLE_GUI_PROFILE=1`, with coarse Puffin spans for app, DNA-window,
  DNA-map, sequence-text, and feature-tree latency investigations.

## 2026-05-28

- Switched the startup splash screen to the app-ready transparent mascot asset
  at `assets/mascots/Mascot_transparent.png` and enlarged its splash
  presentation.
- Defaulted macOS back to native OS child viewports on `egui/eframe` `0.34.3`,
  with `GENTLE_MACOS_HOSTED_CHILD_VIEWPORTS=1` left as the explicit hosted
  fallback and root-window fullscreen restoration disabled on macOS.
- Delayed new macOS native sequence windows until the root window leaves
  fullscreen/maximized state, avoiding the unstable Split View-style flicker
  path when opening the first child window.
- Tightened hosted egui window frame-drag ownership so resize-edge drags keep
  priority over lower hosted-window body hits while ordinary hosted-window body
  interactions such as DNA selection are no longer treated as frame drags.
- Added a six-well cell-culture plate rack profile/template, culture-well SVG
  rendering, arrangement-labelled top-down `racks hero-svg`, `svg-pdf`
  conversion, and README figure assets for plate-layout presentation.

## 2026-05-27

- Added the `codex_local_stdio` Agent Assistant provider and
  `scripts/codex-agent-bridge`, allowing GENtle to route inner-agent requests
  through a logged-in local Codex CLI/App account without requiring
  `OPENAI_API_KEY`.
- Hardened hosted egui window drag ownership during embedded-window resize/move
  interactions so lower hosted windows do not react to a drag that began on the
  active window frame.
- Switched the egui/eframe family from the temporary upstream git override to
  the published `0.34.3` crates while preserving GENtle's hosted-window
  title-bar movement and lower-window drag lockout in the compatibility wrapper.
- Added the shared `mirna` target-site scan service with built-in
  `hsa-miR-96-5p` seed catalog support, JSON schema
  `gentle.mirna_target_scan.v1`, and CLI/shared-shell commands for seed
  explanation, catalog inspection, and annotated target scanning.
- Added `Patterns -> microRNA Target Scan...` as a graphical GUI wrapper over
  the shared scan command, including seed-pairing drawings, region-specific
  splicing interpretation, and side-by-side ortholog/candidate snippet scans.
- Hardened the microRNA scanner with typed evidence tags, a reverse-strand
  coordinate regression, and warnings when a known catalog name is paired with a
  non-canonical mature-sequence override.

## 2026-05-26

- Refreshed `Cargo.lock` after a pre-release `cargo update` and verified the
  release-facing examples/tutorial checks against the updated patch versions.
- Fixed linear DNA-map drag selection so drags that originate on splitters or
  hosted-window decorations cannot become sequence selections after entering
  the map canvas.
- Added `integrations/clawbio/local_agent_handoff.md` as a shared routing note
  for Codex, Claude, OpenClaw, and other local agents: use GENtle through the
  known ClawBio `gentle-cloning` runner, inspect `result.json`/`report.md` and
  artifacts, and avoid inventing a second GENtle command surface.
- Tightened DNA-map drag ownership checks so hosted auxiliary window drags
  cannot pass through into DNA selection.

## 2026-05-24

- Added `File -> New Sequence...` and `File -> New Sequence from Clipboard...`
  so typed or clipboard IUPAC DNA can become an ordinary project sequence
  through `CreateSequenceFromText`.
- Added a phase-1 GENtle-side ClawBio bridge: compact GUI context export,
  `CancelToken`-based subprocess transport, `result.json` parsing, and an
  updated federation/boundary contract in `docs/clawbio_gentle_integration_onepager.md`.
- Added the phase-2 ClawBio GUI panel under `Services -> ClawBio...`, with
  worker-thread dispatch, cancellation, stdout/stderr output-bundle logs,
  artifact/report links, and verbatim suggested-action request dispatch.
- Fixed the External Services provider/service picker so changing from
  metabion to GeneArt refreshes the editable request JSON and clears stale
  preflight/quote previews.

## 2026-05-23

- Started Phase 4 of the tutorial presentation overhaul by surfacing
  `review_manifest.json` status in `catalog.json`, generated tutorial
  front matter, generated and human tutorial hubs, and GUI tutorial hover text.
- Completed the second Phase 4 tutorial review slice: catalog/generated
  tutorial outputs now carry dependency-aware stale reasons, feedback issue
  template links, and the two May 18 human reviews are recorded as stale after
  the subsequent tutorial source changes.
- Began Phase 5 of the tutorial presentation overhaul by carrying declared
  `graphics[]` metadata through catalog/manifest/report outputs, embedding the
  promoter-design generated TFBS SVG beside its tutorial step, and tracking one
  existing screenshot-backed dotplot tutorial figure as review-staleness input.

## 2026-05-22

- Added a Phase 0 tutorial coverage audit for the decimal-numbering overhaul,
  identifying uncataloged tutorial-like Markdown pages and the TP73/VKORC1
  promoter-luciferase decision gate before schema or rename work begins.
- Began Phase 1 of the tutorial presentation overhaul with v2 catalog/manifest
  schemas, v4 source-unit grouping/graphics fields, a 10-group tutorial
  taxonomy, and two newly cataloged guided walkthroughs for stateless sequence
  inspection and TFBS similarity ranking.
- Completed the Phase 1 tutorial group-assignment review gate: all 40 source
  units now use `gentle.tutorial_source.v4`, 38 tutorial units carry derived
  decimal ids, and the generated hub plus landscape overview remain grouped but
  intentionally unnumbered.
- Began Phase 2 tutorial presentation regrouping: the human tutorial hub,
  generated reference hub, Help tutorial picker, and `File -> Open Tutorial
  Project...` menu now present tutorials by content group and decimal id rather
  than flat lists or tier-only buckets.
- Began Phase 3 tutorial filename convergence: numbered tutorial sources,
  generated chapters, and hand-written walkthroughs now use decimal
  `<GG>-<PP>_<unit_id>` filenames, while unnumbered reference units keep
  explicit `ref_*.json` source names.
- Archived superseded TP73 promoter/luciferase tutorial pointers and moved the
  planned VKORC1 / Factor X companion tutorial note into the roadmap parking
  lot instead of active tutorial numbering.
- Added deterministic helper vector-card and catalog-doctor inspection routes
  with GUI Planning-window parity for metadata-only vector candidates.

## 2026-05-21

- Added metadata-only Luebeck S1 vector/resource candidates to the helper
  catalog, with explicit sequence availability, redistribution, biosafety,
  procurement, host-compatibility, marker/origin, and empty-backbone fields.

## 2026-05-19

- Added GUI export for lab-assistant cloning handoff reports and upgraded
  `ExportLabAssistantInstructions` to `gentle.lab_assistant_instructions.v2`
  with ODT/DOCX editable reports plus embedded lineage graphics where
  rasterization is available.
- Made `gentle_examples_docs tutorial-check` failures include paste-ready
  tutorial feedback context with chapter/source/workflow/artifact paths and a
  suggested issue-template category.
- Added `gentle.tutorial_review_manifest.v1` review metadata with nonfatal
  tutorial-check warnings for missing, unknown, or stale tutorial review
  entries.
- Regenerated executable tutorial chapters with review/automation provenance
  front matter and a standardized tutorial feedback section.
- Added a GUI Tutorial help button that copies issue-ready feedback context
  with tutorial id, source JSON, workflow, artifact, version, platform, and
  current search context.
- Grouped the in-app Tutorial topic picker by catalog audience buckets and
  added standardized feedback sections to hand-written tutorial pages.
- Added the read-only `planning consult cloning` route and
  `gentle.planning_cloning_consultation.v1` report so Agent Assistant and
  ClawBio can quote deterministic cloning-strategy/vector advice instead of
  improvising prose.
- Pinned `planning consult cloning` v1 to the 11 catalogued routine families,
  kept `--seq-id` traceability-only, and left marker/promoter/MCS gaps as
  explicit questions instead of narrative ranking heuristics.
- Clarified the Agent Assistant tutorial's same-live-state relationship with
  the GUI and GUI Shell, and added `CreateSequenceFromText` /
  `sequence create --sequence-text ...` so agent-proposed sequence text can be
  materialized as an ordinary project sequence.
- Added direct OpenAI usage and billing hyperlinks to GUI quota-failure
  messages from the Agent Assistant, including `Test Setup` live-probe results
  whose static configuration is otherwise available.
- Tightened the native Agent Assistant system prompt so providers are told to
  suggest only GENtle shared-shell commands, avoid invented `fs.*`/gateway-style
  verbs, ask for local file paths when discovery is outside the GENtle contract,
  and use ASCII punctuation for egui-safe rendering.
- Documented the current non-OpenClaw command boundary for the in-app Agent
  Assistant, including Finder/Spotlight as the macOS path-discovery fallback
  before handing a selected file path to GENtle.
- Added a protocol-owned GENtle-local slash alias registry covering `/help`,
  `/list`, GUI file opening, exact file import, pasted sequence creation, and
  explicit external fetch aliases while keeping vague filesystem/OpenClaw-style
  commands rejected with typed alternatives.
- Added `Ctrl+Return` submission for the Agent Assistant prompt editor.
- Added a `Copy Response JSON` button for the latest Agent Assistant response.
- The Agent Assistant suggestion table now validates suggested commands before
  enabling `Run`, marking hallucinated or unsupported commands as invalid.

## 2026-05-18

- Split the tutorial hub into hand-written guided walkthroughs versus
  machine-generated executable reference chapters, with reciprocal links for
  generated chapters that have guided counterparts.
- Added native Mistral Agent Assistant integration with a quick-start GUI route,
  `MISTRAL_API_KEY` session/env handling, model discovery, and live setup
  preflight parity.
- Added tutorial stewardship rules so recurring maintenance checks tutorial
  implementation drift, wet-lab readability, biological teaching depth, and
  Codex-versus-human review sign-off boundaries.
- Kept the Agent Assistant and Agent Interfaces tutorial pinned in the
  in-app `Help -> Tutorials` picker even if catalog discovery falls back to
  markdown scanning.
- Enriched the first generated tutorial batch with v3 per-step CLI snippets,
  expected-result callouts, and prerequisite links.
- Enriched the Gibson generated tutorial batch with v3 per-step CLI snippets,
  expected-result callouts, and prerequisite links.
- Enriched the online generated tutorial batch with v3 per-step CLI snippets,
  expected-result callouts, and prerequisite links to reference preparation.
- Completed generated tutorial-source enrichment so every executable chapter
  now has v3 per-step CLI/expected-result teaching fields.
- Polished generated tutorial rendering with quieter provenance in Help,
  applied-concept wording, earlier parameter explanation, CLI state guidance,
  SVG label fallbacks, and capped tutorial-project progress before final open.
- Surfaced documentation-only guided walkthroughs, including the Agent
  Assistant and Agent Interfaces tutorial, from `File -> Open Tutorial
  Project...` so readers can find them from the tutorial-opening menu as well
  as `Help -> Tutorials`.
- Fixed `File -> Open Tutorial Project... -> Guided walkthroughs` to read all
  `manual/reference` tutorial entries from `docs/tutorial/catalog.json`
  instead of showing only the hard-coded agent tutorial entry.

## 2026-05-17

- Added v3 tutorial-source teaching fields for generated chapters, including
  per-step CLI/expected-outcome rendering, prerequisite links, and online
  local-execution callouts while keeping v2 sources loadable during migration.
- Added a reporter construct handoff tutorial that shows how to generate and
  inspect `gentle.reporter_construct_handoff.v1` from a saved
  promoter-reporter candidate set before running any macro commands.

## 2026-05-16

- Restored shared-shell reachability for the reporter catalog/recommender
  family so `reporters list`, `reporters recommend`, and
  `reporters export-corpus` now work through GUI Shell and
  `gentle_cli shell ...` as documented.

## 2026-05-15

- Regenerated the GUI/CLI/MCP parity matrix from `CapabilityDescriptor`
  surfacing metadata and added a freshness guard plus regeneration script.
- Replaced binary capability adapter exposure with per-adapter surfacing states
  and added `EngineError.cause_chain` preservation across adapter boundaries.
- Made engine-backed missing adapter routes render as parity gaps unless a
  curated not-applicable override supplies a real justification.
- Added `gui-menu` glossary surfacing and GUI menu/palette oracle checks so the
  parity matrix can distinguish first-class GUI affordances from GUI shell reachability.
- Added conservative RNA-read isoform triage TSV export for known-isoform,
  ambiguous, gene-supported/no-call, and off-target/bad-seed read bins without
  calling novel isoforms.
- Added Splicing Expert isoform read-support inspection using saved RNA-read
  mapped-isoform triage counts, with audit/export links back to contributing
  aligned reads.
- Fixed Agent/Routine Assistant hosted-window focus lookup so embedded macOS
  windows are raised by their stable hosted ids instead of stale title layers.
- Added an offline reporter-recommender V1 with a provenance-gated local
  reporter catalog, deterministic constraint ranking, rejected-candidate
  reasons, and JSON/JSONL corpus export for local AI retrieval or training prep.
- Added a read-only reporter construct handoff plan that joins saved
  promoter-reporter candidate reports, the offline luciferase recommender, and
  the existing `allele_paired_promoter_luciferase_reporter` macro without
  creating constructs automatically.

## 2026-05-14

- Made `AlignSequences` accept inline ASCII `SequenceScanTarget` operands for
  stateless pairwise alignment while preserving legacy `query_seq_id` /
  `target_seq_id` workflows.
- Added a protocol-owned capability registry for CLI/MCP/JavaScript/Lua
  discovery and normalized CLI/MCP adapter failures into structured
  `EngineError` payloads.
- Added a short in-tree agent development loop with `scripts/dev-gentle-cli`,
  `gentle_cli doctor --agent`, and a dedicated `docs/agent_dev_loop.md`.
- Hardened MCP capability discovery so `tools/list` carries glossary-backed
  descriptions, input/output schema metadata, mutating/external safety flags,
  and a documented exclusion ledger for shell commands without dedicated tools.
- Added generated tutorial chapters that still need explicit human functional
  confirmation to the tutorial catalog and Help -> Tutorials picker with the
  `generated+checked/human-pending` status.

## 2026-05-13

- Selected the TP73 genome-anchored evidence viewer as the next release aim and
  added the offline proof workflow, public runbook, fixture provenance, and
  generic genome-track detail polish needed for the first smoke path.
- Added a GUI/CLI/MCP parity audit and brought the glossary-backed help catalog
  in sync for implemented agent preflight/model/planning routes plus
  construct-reasoning graph routes.
- Fixed hosted Promoter design foregrounding so DNA-owned auxiliary windows now
  open above their parent sequence window instead of behind it.
- Added a dedicated macOS `egui` viewport investigation pack plus a separate
  source folder for the minimal repro harness.
- Added an optional `--expression-tsv` heatmap overlay to isoform architecture
  SVG exports, keyed by `isoform_id`/`sample_label` values.

## 2026-05-12

- Cleared the completed interim-release block from `docs/roadmap.md` and reset
  the roadmap to choosing the next release aim.

## 2026-05-11

- Revised the GUI egui stack from `0.34.1` to `0.34.2` across `eframe`,
  `egui`, `egui_extras`, and the `gentle-gui` helper crate.
- Added the ClawBio `mode=intents` runtime descriptor surface, trigger-keyword
  generation/drift checks, and a provider-neutral default scope for
  `services handoff`.
- Added `rna-reads show-alignments` for headless batch export of the same
  per-read alignment display payload used by `rna-reads show-alignment` and the
  GUI splicing-expert detail pane.

## 2026-05-10

- Moved sequence file-loading dispatch into `dna_sequence::load_from_file` so
  engine, JS, Lua, and tests no longer depend on `GENtleApp` for parser
  selection.
- Added a TP73 pancreas proof-bundle `audit` handoff to
  `scripts/tp73_pancreas_rna_mapping.sh`, including gene-group snapshot capture
  and expanded evidence-bundle paths for release review.

## 2026-05-09

- Aligned the `docs/release.md` Local Pre-Tag Smoke Checklist with the roadmap
  release gate by adding `cargo check -q` before the release build matrix.
- Split roadmap maintenance into `roadmap.md` for next work, this changelog for
  completed work, and `decisions.md` for durable implementation decisions.
- Added [`maintenance_chore_plan.md`](maintenance_chore_plan.md) to define
  five recurring repository hygiene chores: session close, daily bug scan,
  weekly drift/test/provenance scan, release-gate readiness, and monthly
  decisions/parity review.
- Added `scripts/maintenance_chore.py session-close` as the first runnable
  maintenance chore implementation, covering roadmap size, completed-history
  bullets, invariant links, conflict markers, worktree/artifact warnings, and
  plan-fidelity handoff evidence.
- Polished the session-close runner after external review: conflict-marker
  detection now uses the exact seven-equals Git separator, split-doc existence
  is checked explicitly, and plan fidelity is reported as a manual reminder.
- Added `scripts/maintenance_chore.py drift-scan --base-ref REF` as the first
  weekly drift/test/provenance automation slice, warning when shared
  engine/protocol/shell-contract files change without obvious test evidence in
  the same diff.
- Polished `drift-scan` after external review by narrowing the warning to
  shared-contract code without separate test-file evidence, recognizing plural
  Rust `_tests.rs` files, normalizing rename status display, and documenting
  inline-test and pure protocol-doc limits.
- Added `scripts/maintenance_chore.py release-gate` as the first
  release-readiness automation slice, checking compact roadmap release-gate
  structure, local release-gate links, release-facing document presence, and
  `cargo check` alignment with `docs/release.md`.
- Hardened the release-gate scan after external review so a missing roadmap
  returns a structured failure and titled Markdown links do not produce false
  broken-link findings.
- ClawBio/OpenClaw `gentle-cloning` capabilities replies gained shared
  UI-intent catalog handoff support through `ui_intent_catalog`,
  `ui_intent_catalog_error`, and `kind = ui_intent` suggested actions.
- Stash-pop conflicts around UI-intent roadmap/docs wording were resolved by
  preserving both shared GUI catalog consumption and ClawBio/OpenClaw handoff
  documentation.

## 2026-05-08

- Release planning narrowed pre-release work to TP73 pancreatic Nanopore cDNA
  proof artifacts, release-signoff tests, GUI smoke, and documentation
  alignment.
- Gene-agnostic pancreas RNA screening was routed through
  `scripts/pancreas_gene_rna_screen.sh` for seed-only cohort triage and
  optional retained-read alignment.

## 2026-05-02

- Restriction-site inspection gained pinned Description-panel details plus
  `Copy Summary` and `Copy Detail JSON` actions for protocol notes and agent
  handoff.
- Restriction-site context menus gained pin/copy actions for hovered or
  selected labels without changing the rendering layer.

## 2026-05-01

- Restriction-site hover tooltips began reusing the shared
  `RestrictionSiteExpertView`, including enzyme grouping, cut geometry,
  tooltip lines, and REBASE links.
- `restriction_site_detail` exposed the same restriction-site detail record
  through MCP for non-rendered agent inspection.
- GUI Agent Assistant suggested-command execution and shared plan execution
  reject nested `agents ...` invocations across the full ask/plan/execute
  family.

## 2026-04-28

- RepeatMasker/UCSC `rmsk`-style repeat annotations gained deterministic
  labels and subtype colors across GUI maps, feature-tree grouping/filtering,
  and SVG export.

## 2026-04-21

- Repeat/mobile-element reasoning gained direct Dotplot and RevComp Dotplot
  actions from the DNA-window inspector.

## 2026-04-18

- Construct-reasoning annotation candidates gained accepted/rejected/locked
  review state, shared shell/CLI mutation routes, GUI curation, and JS/Lua/MCP
  parity.
- Accepted or locked construct-reasoning annotation candidates can be written
  back as ordinary sequence features through a shared engine report.
- Construct-reasoning graphs gained portable annotation-candidate summaries so
  long genomic views can show collapsed context instead of raw overlap lists.
- Similarity/repeat/mobile-element reasoning began emitting generated
  repeat-region annotations and operational-risk facts for PCR, mapping,
  nanopore, and cloning review.

## 2026-04-17

- Construct-reasoning overlays were narrowed to annotation-grade and
  decision-linked spans so raw restriction/TFBS evidence no longer floods long
  genomic windows by default.
- Construct-reasoning graphs gained a portable annotation-candidates layer
  shared by the sequence overlay and inspector.

## 2026-04-16

- Adapter/linker restriction-capture reasoning was added to construct
  objectives, GUI summaries, run bundles, shared shell/CLI, JS, Lua, and MCP.

## 2026-04-13

- Variant markers became first-class construct-reasoning evidence with
  promoter/TFBS/CDS/UTR/splice effect hypotheses and first assay-family
  suggestions.
- Routine planning began consuming transcript-aware variant summaries and
  construct-reasoning context for `routines list`, `routines explain`, and
  `routines compare`.

## 2026-04-12

- Construct-reasoning contracts reserved host/helper/growth fields and began
  emitting non-sequence context facts and decisions.
- Host-route restriction/methylation, selection/complementation, helper-backed
  selection, and growth-condition interpretation were added to shared
  construct-reasoning outputs.
- Process run-bundles gained a portable `construct_reasoning` section for
  ClawBio/OpenClaw offline inspection.
- DNA-window inspection gained non-sequence construct-reasoning summaries while
  keeping overlays sequence-only.
- A starter host-profile catalog and GUI browser were added.

## 2026-04-10

- Construct-reasoning overlays gained hover/click inspection, side-panel
  evidence details, and GUI filters without mutating the engine-owned graph.

## 2026-04-09

- Construct-reasoning graph records and metadata storage were introduced for
  deterministic read-only reasoning over existing sequence evidence.
- The DNA window began auto-refreshing a read-only reasoning graph and painting
  selected reasoning spans through a `Reasoning` layer.

## 2026-03-19

- Async BLAST job durability was hardened with persisted job snapshots,
  deterministic restart/reload recovery, cancellation semantics, and
  conformance coverage.
- Prepare-job cancellation and timeout classification became engine-owned and
  deterministic.
- BLAST and related genome routes gained external-binary preflight diagnostics
  for `blastn` and `makeblastdb`.

## 2026-03-04

- Malformed GTF/GFF annotation reporting began summarizing non-fatal parse
  issues with file/line context in prepare/report payloads.

## 2026-03-03

- JS adapter parity tests began comparing shared-shell-backed import-pool,
  REBASE sync, and JASPAR sync outcomes.
- Lua adapter parity tests began comparing import-pool, REBASE sync, and
  JASPAR sync outcomes.
- CLI forwarded-dispatch parity tests began verifying import-pool and resource
  sync routes against shared-shell execution.

## 2026-02-27

- MCP UI-intent parity tests began comparing discovery, prepared-query,
  latest-selection, and open-intent outputs against shared-shell `ui ...`
  execution.

## 2026-02-24

- Shared-shell execution tests began covering local-fixture REBASE and JASPAR
  resource sync without network dependency.

## Undated Historical Baseline Imported From Roadmap

- Shared engine operation/workflow execution exists in `src/engine.rs` with
  adapters for GUI, CLI, shared shell, JS, Lua, Python wrapper, and MCP.
- Build tooling defaults to per-worktree Cargo target directories, while
  release-installer CI avoids stale `target/` caches for tag/manual builds.
- Multi-crate extraction has started around `gentle-protocol`,
  `gentle-engine`, `gentle-render`, `gentle-shell`, and `gentle-gui`, with
  portable contracts and selected render/shell/gui helpers already moved.
- Engine, engine-shell, and GUI module decomposition has begun, including
  extracted tests, helper modules, and navigation-oriented module docs.
- Shared shell and CLI cover capabilities, state summary, genomes, helpers,
  resources, tracks, ladders, candidates, services, protocol cartoons,
  RNA-read routes, CUT&RUN routes, sequencing confirmation, UI intents, and
  many direct rendering/export paths.
- Python integration exists as `integrations/python/gentle_py`, wrapping
  deterministic CLI contracts.
- ClawBio/OpenClaw integration exists as a copy-ready `gentle-cloning` skill
  scaffold with local checkout and Apptainer launchers, example requests,
  reproducibility bundles, PNG-first graphics bundling, lifecycle-aware
  suggested actions, and first-run bootstrap documentation.
- MCP server baseline exposes deterministic tool routes including
  `capabilities`, `state_summary`, guarded `op`/`workflow`, `help`,
  UI-intent tools, prepared-genome helpers, BLAST async tools, and
  restriction-site detail.
- GUI Agent Assistant and typed planner/execution boundaries exist; machine
  callers are expected to prefer typed `agents plan` / `agents execute-plan`
  style contracts over chat-only flows.
- GUI command palette, menus, shell, agent, MCP, and ClawBio discovery now
  consume the shared UI-intent target catalog instead of maintaining separate
  target lists.
- Feature expert views for TFBS, restriction, splicing, isoform panels, and
  related SVG exports are engine-owned and exposed through multiple adapters.
- Candidate-set scoring/filter/set algebra, guide-set operations, practical
  guide filtering, oligo/protocol export, and design-constraint filtering are
  first-class engine operations.
- Reference/helper genome catalog, preparation, extraction, BLAST, anchored
  track import, subscription, and cache cleanup surfaces are implemented across
  engine/shell/CLI/GUI paths.
- Genome-anchor resolution supports deterministic fallback policy controls.
- BigWig/VCF/BED/BLAST genome-track imports and display filters are shared
  engine/display contracts.
- SnapGene `.dna` import is supported through a reusable headless parser crate
  with synthetic fixtures alongside GenBank/EMBL/XML parity fixtures.
- XML import is additive and GenBank-like in semantics, with XML dialect
  rejection/parity tests.
- Process run-bundle export, protocol-cartoon rendering, lab-assistant
  instruction export, routine decision traces, and protocol template export
  have shared record/export baselines.
- Cloning routine catalog, macro templates, macro validation, routine explain
  and compare, Routine Assistant GUI baseline, macro lineage visualization, and
  baseline restriction/Gibson/Golden Gate/Gateway/TOPO/TA/GC/In-Fusion/
  NEBuilder HiFi packs are present.
- Gibson specialist, Gibson preview/apply, arrangement/rack/ladder state,
  rack-profile authoring, label export, fabrication/export, OpenSCAD export,
  and rack simulation JSON paths are present.
- Virtual pool gels, protein gels, protease digest gels, 2D protein gels,
  arrangement-aware gel exports, and text band rows are present.
- Primer design, qPCR assay design, primer specificity assessment, Primer3
  backend baseline, backend preflight, and progress parity exist, with richer
  constraint parity still tracked as future work.
- PCR, advanced PCR, mutagenesis, overlap-extension mutagenesis,
  insertion-primer design, restriction-cloning PCR handoff, transcript-derived
  cDNA assays, and qPCR setup/promotion helpers are present.
- UniProt, Ensembl gene/protein, GenBank, dbSNP, and related projection/import
  paths exist with persisted reports and provenance for several artifact
  families.
- Alternative-splicing interpretation baseline exists with transcript lanes,
  boundaries, event/evidence summaries, primary map read-only mode, and expert
  SVG parity.
- RNA-read interpretation includes seed/preflight workflows, retained-read
  alignment, saved report stores, GUI mapping workspace, CLI/ClawBio examples,
  target-gene cohort summaries/audits, score-density export, and batch mapping
  helpers.
- CUT&RUN support includes dataset status/prepare/project/read interpretation,
  ROI read reports, coverage/export, regulatory-support inspection, and GUI
  inspector baselines.
- Isoform architecture panels, curated TP53/TP73 panel resources, online/offline
  workflows, expert SVG exports, protein/domain lanes, and validation reports
  are present.
- TFBS/JASPAR resources include catalog listing, remote metadata sync,
  deterministic registry benchmarks, single-entry inspection, promoter-window
  summaries, multi-gene promoter views, score tracks, score-track similarity,
  and query resolution.
- Promoter and variant context work includes variant promoter summaries,
  reporter-fragment suggestions, allele materialization, luciferase planning
  examples, and regulatory evidence ledgers.
- Dense rendering controls, adaptive linear DNA letter routing, scroll/zoom
  policy, display visibility, topology, viewport controls, and GUI declutter
  baselines are present.
- Agent recursion guardrails block nested `agents ...` execution directly and
  through macro expansion paths.
- Test-data provenance policy, remote-resource skip policy, source-documentation
  pilot docs, and cargo-doc module documentation goals are established.
- Screenshot artifact support exists historically but remains disabled by
  security policy unless explicit opt-in/approval is restored.
