# GENtle Engine Protocol (Draft v1)

This document defines the draft machine-facing protocol for operating GENtle
through a shared core engine.

Goal:

- GUI, CLI, JavaScript, Lua, and Python wrappers call the same core routines.
- AI tools can run deterministic cloning workflows with reproducible logs.

## Design principles

- Protocol-first: versioned JSON request/response shapes
- Capability negotiation: clients discover supported operations and formats
- Deterministic operation log: each operation emits a stable op id and result
- Structured errors: machine-parseable error code + message

## Capabilities

`gentle_cli capabilities` returns:

- `protocol_version`
- `supported_operations`
- `supported_export_formats`
- `deterministic_operation_log`

## Protocol-first workflow examples

Canonical, adapter-independent examples are defined in:

- `docs/examples/workflows/*.json`
- schema: `gentle.workflow_example.v1`

Each example includes:

- metadata (`id`, `title`, `summary`)
- test policy (`test_mode`: `always|online|skip`)
- required local files (`required_files`)
- canonical `workflow` payload

Adapter snippets (CLI/shared shell/JavaScript/Lua) are generated on demand from
those canonical files:

```bash
cargo run --bin gentle_examples_docs -- generate
```

Validation only:

```bash
cargo run --bin gentle_examples_docs -- --check
```

Tutorial manifest + generated outputs:

- discovery catalog: `docs/tutorial/catalog.json`
- discovery schema: `gentle.tutorial_catalog.v1`
- shared tutorial source units:
  - `docs/tutorial/sources/catalog_meta.json`
  - `docs/tutorial/sources/*.json`
- source-unit schemas:
  - `gentle.tutorial_catalog_meta.v1`
  - `gentle.tutorial_source.v2`
- generated runtime manifest: `docs/tutorial/manifest.json`
- runtime manifest schema: `gentle.tutorial_manifest.v1`
- committed generated outputs: `docs/tutorial/generated/`

Catalog/manifest split:

- `docs/tutorial/catalog.json` is the canonical discovery layer for all
  tutorials, including hand-written walkthroughs and agent/reference guides.
- `docs/tutorial/sources/` is the authoring layer for both the discovery
  catalog and the executable tutorial runtime manifest.
- `docs/tutorial/manifest.json` is a generated runtime contract used for
  chapter output and tutorial runtime checks.
- GUI help/tutorial discovery may consume the catalog directly for curated
  ordering and metadata, while executable tutorial project materialization still
  resolves through the manifest/workflow example path.

Generate/check tutorial outputs:

```bash
cargo run --bin gentle_examples_docs -- tutorial-generate
cargo run --bin gentle_examples_docs -- tutorial-check
cargo run --bin gentle_examples_docs -- tutorial-catalog-generate
cargo run --bin gentle_examples_docs -- tutorial-catalog-check
cargo run --bin gentle_examples_docs -- tutorial-manifest-generate
cargo run --bin gentle_examples_docs -- tutorial-manifest-check
```

## Draft design resources

### `gentle.gibson_assembly_plan.v1`

Purpose:

- describe one Gibson cloning project in a destination-first way,
- separate user-specified plan inputs from derived design consequences,
- provide one canonical JSON artifact that future routines, primer design, and
  protocol-cartoon rendering can all read from.

Status:

- draft schema and documentation artifact only
- not yet accepted by a direct engine operation in this release

Canonical examples:

- `docs/examples/plans/gibson_destination_first_single_insert.json`
- `docs/examples/plans/gibson_destination_first_multi_insert.json`

Top-level structure:

- `schema`, `id`, `title`, `summary`
- `destination`
  - destination molecule (`seq_id`, prior topology)
  - explicit opening definition (`mode`, label, resulting left/right ends)
- `product`
  - intended output topology and output-id hint
- `fragments[]`
  - participating inserts or non-destination fragments
  - orientation plus per-end adaptation strategy
  - optional source-coordinate hints for deterministic ordering
- `assembly_order[]`
  - explicit left-to-right order of destination ends and inserts
  - supports future multi-fragment Gibson plans without changing the model
- `junctions[]`
  - one record per adjacent join
  - required overlap length
  - explicit overlap partition across the left/right adjacent members
  - whether overlap is derived from destination context or user-specified
  - explicit `distinct_from` constraints for terminal junctions
- `validation_policy`
  - hard requirements:
    - unique/unambiguous destination opening
    - distinct terminal junctions
    - adjacency-consistent overlaps
  - advisory checks:
    - overlap-length design range
    - fragment-count-aware overlap targets
    - overlap Tm
    - destination/fragment/reference uniqueness heuristics
- `derived_design`
  - derived overlap sequences
  - primer design suggestions
  - advisory notes and validation outcomes

Current draft value vocabulary:

- `destination.topology_before_opening`
  - `linear`
  - `circular`
- `destination.opening.mode`
  - `defined_site`
    - a user-selected opening site/window on an existing destination molecule
  - reserved future values:
    - `restriction_digest`
    - `pcr_linearization`
    - `inverse_pcr`
- `destination.opening.uniqueness_requirement`
  - `must_be_unambiguous`
    - opening ambiguity is a hard validation error
  - `advisory_only`
    - ambiguity is surfaced, but not automatically fatal
- `product.topology`
  - `linear`
  - `circular`
- `fragments[].role`
  - `insert`
  - reserved future values:
    - `backbone_fragment`
    - `bridge_fragment`
- `fragments[].orientation`
  - `forward`
  - `reverse`
- `fragments[].source_span_1based`
  - optional source-coordinate hint for plans that preserve source order
  - shape:
    - `source_seq_id`
    - `start`
    - `end`
- `fragments[].left_end_strategy.mode` / `right_end_strategy.mode`
  - `native_overlap`
    - fragment terminus is already expected to satisfy the required overlap
  - `primer_added_overlap`
    - overlap is expected to be introduced via a primer tail
    - this is the current draft bucket for overlap-extension / primer-stitching
      style adaptation
  - reserved future values:
    - `synthetic_terminal_sequence`
    - `library_defined_overlap`
- `assembly_order[].kind`
  - `destination_end`
  - `fragment`
- `junctions[].overlap_source`
  - `derive_from_destination_left_flank`
  - `derive_from_destination_right_flank`
  - `derive_from_adjacent_fragment_ends`
  - `designed_bridge_sequence`
    - internal junction overlap chosen as a synthetic bridge/adaptor sequence
  - reserved future value:
    - `user_specified_sequence`
- `junctions[].overlap_partition`
  - explicit contribution of the overlap region from the adjacent members
  - shape:
    - `left_member_bp`
    - `right_member_bp`
  - invariant:
    - `left_member_bp + right_member_bp == required_overlap_bp`
  - examples:
    - left-member only overlap: `30 + 0`
    - right-member only overlap: `0 + 30`
    - split overlap: `20 + 20`
- `validation_policy.adjacency_overlap_mismatch`
  - `error`
  - `warn`
- `validation_policy.uniqueness_checks.*`
  - `off`
  - `warn`
  - `error`
- `validation_policy.reference_contexts[].severity`
  - `warn`
  - `error`

Input-vs-derived boundary in the draft model:

- Intended user/planner inputs:
  - `destination`
  - `product`
  - `fragments`
  - `assembly_order`
  - `junctions[].required_overlap_bp`
  - `junctions[].overlap_partition`
  - `junctions[].distinct_from`
  - `validation_policy`
- Intended normalized/derived outputs:
  - `derived_design.junction_overlaps`
  - `derived_design.primer_design_suggestions`
  - `derived_design.notes`
- Transition fields that may begin as user hints and later become resolved
  values:
  - `junctions[].overlap_source`
  - fragment end strategies (`native_overlap` vs `primer_added_overlap`)

Interpretation:

- Gibson plans are modeled as explicit assembly junctions around an opened
  destination, not merely as an unordered bag of fragments.
- The destination opening defines two terminal junctions, and therefore two
  required overlap regions.
- Inserts may already satisfy those terminal overlaps or may require primer-tail
  adaptation.
- The overlap at one junction should be treated as a selection around the
  in-silico junction, not merely as one scalar length:
  - it may come entirely from the left member
  - entirely from the right member
  - or be split across both members
- For plans that preserve an existing source order, fragment ordering should
  follow ascending bp coordinates (low bp positions first) unless an explicit
  alternative order is requested.
- For multi-fragment plans, internal fragment-fragment junctions may be created
  through primer-added bridge overlaps rather than relying on pre-existing
  native overlap.
- Uniqueness is best treated in layers:
  - destination opening uniqueness: hard validation
  - left/right terminal overlap distinctness: hard validation
  - destination/fragment/genome uniqueness heuristics: advisory checks

Design intent:

- make the same JSON artifact useful for:
  - preflight Gibson validation
  - primer design derivation
  - workflow/macro instantiation
  - factual protocol-cartoon generation
  - reproducible AI-facing project context

Practical overlap heuristics (draft defaults):

- single-insert / two-fragment style assemblies often fit comfortably in the
  20-40 bp range
- multi-fragment assemblies should usually move toward longer overlaps
- a practical starting rule for the draft model is:
  - 1-2 assembled fragments: 20-40 bp overlaps
  - 3-5 assembled fragments: 40 bp overlaps
  - 6+ assembled fragments: 50-100 bp overlaps
- internal multi-fragment junctions introduced by primer-added bridge overlaps
  should follow the same fragment-count-aware guidance rather than being treated
  as exempt from overlap design heuristics

Primer design conventions (draft):

- Gibson primer suggestions should be modeled as two-part primers:
  - `overlap_5prime`
    - non-priming 5' overlap segment used for assembly of adjacent fragments
  - `priming_3prime`
    - gene-specific 3' priming segment used for PCR amplification from the
      source template
- Primer design should start from an in silico assembled product/junction view,
  then work backward to fragment-specific PCR primers.
- Overlap choice is best treated as Tm-aware rather than length-only:
  - simple PCR-fragment-into-vector assemblies may use shorter overlaps when
    overlap Tm is already adequate
  - more complex or multi-fragment assemblies often justify longer overlaps
- The overlap region may lie entirely within one adjacent member or be split
  across the two members around a junction.
- Two primers that implement the same overlap sequence can still belong to
  different PCR reactions, because each primes a different template fragment.

Assembly setup heuristics (draft advisory layer):

- linearized destination can be prepared by PCR amplification or by restriction
  digestion
- PCR cleanup is not always required, but carryover should stay modest
  relative to the final assembly reaction volume
- column purification is especially worth recommending for:
  - assemblies of three or more PCR fragments
  - assemblies involving fragments longer than ~5 kb
- direct vector + insert assemblies often benefit from insert concentration
  above vector concentration
- multi-fragment vector assemblies should generally move toward equimolar
  fragment usage
- some constructs may validate in silico but still perform poorly because of
  biological burden or instability in the propagation host
  (for example repeats or toxic products)

Planning implication:

- these factors should usually surface as `derived_design` advisories or future
  execution/setup guidance, not as hard failures in the core Gibson junction
  model

Normalization/derivation phases:

1. Resolve the destination opening into explicit `dest_left` / `dest_right`
   terminal context.
2. Normalize `assembly_order[]` into one adjacency chain.
   - when fragments carry compatible `source_span_1based` hints for one source
     context, default normalization should preserve ascending bp order
3. Materialize one `junction` per adjacent pair in that chain.
4. Derive required overlap sequences from destination flanks and/or adjacent
   fragment termini.
   - respect the junction-specific overlap partition when choosing the final
     overlap sequence around that adjacency
   - internal multi-fragment junctions may instead use designed bridge
     sequences introduced by primer-added overlaps
5. Detect whether each fragment end already satisfies its required overlap or
   requires adaptation (for example primer-added tails).
6. Run hard validation and advisory design checks.
7. Expose derived overlaps, primer design suggestions, and cartoon-ready event
   semantics through `derived_design`.
8. Attach reaction/setup advisories (cleanup, stoichiometry, host-risk notes)
   without conflating them with the hard overlap/junction logic.

Current invariants for the draft model:

- `assembly_order[]` defines the intended adjacency order explicitly.
- `junctions[]` should cover every adjacent pair in `assembly_order[]`.
- terminal junctions are the ones adjacent to the opened destination ends.
- `junctions[].overlap_partition.left_member_bp +
   junctions[].overlap_partition.right_member_bp` should equal
  `junctions[].required_overlap_bp`.
- terminal junction distinctness is a hard validation rule for opened
  destination-vector Gibson plans.
- when source-order hints are present and no contrary manual order is given,
  low bp positions should precede high bp positions in normalization.
- destination-opening uniqueness is a hard validation rule.
- broader destination/fragment/genome uniqueness checks are advisory unless a
  stricter policy is requested.
- `derived_design` may contain unresolved/null sequences at pure planning time;
  this allows the same schema to exist before sequence extraction or primer
  design has been run.
- `junctions[].distinct_from` is currently intended primarily for terminal
  destination-defined junctions, not as a requirement that every internal
  fragment-fragment junction be globally unique.
- `native_overlap` is an expectation about the fragment terminus; it still
  requires sequence confirmation in validation/derivation.
- `designed_bridge_sequence` should be treated as a designed internal overlap,
  suitable for primer-stitching style workflows; GENtle should validate it for
  distinctness/design heuristics rather than treating it as biologically
  privileged just because it was user-supplied.

## Core entities

### ProjectState

```json
{
  "sequences": {"seq_id": "DNAsequence object"},
  "metadata": {"any": "json"},
  "display": {"ui_visibility_tfbs_and_linear_viewport_state": "..."},
  "lineage": {"nodes": {}, "edges": []},
  "parameters": {"max_fragments_per_container": 80000},
  "container_state": {"containers": {}, "seq_to_latest_container": {}}
}
```

Semantic interpretation:

- In GUI terms, a project window represents a wet-lab container context.
- A container may map to multiple candidate sequences/fragments.
- Explicit container objects are first-class state (`container_state`) and are
  indexed from sequence ids via `seq_to_latest_container`.

### Operation

Current draft operations:

- `LoadFile { path, as_id? }`
- `SaveFile { seq_id, path, format }`
- `RenderSequenceSvg { seq_id, mode, path }`
- `RenderDotplotSvg { seq_id, dotplot_id, path, flex_track_id?, display_density_threshold?, display_intensity_gain? }`
- `RenderFeatureExpertSvg { seq_id, target, path }`
  - shared renderer contract across GUI/CLI/JS/Lua for TFBS/restriction/splicing/isoform expert exports
  - splicing SVG includes explicit junction-support counts, frequency-encoded transcript-vs-exon matrix coloring, predicted exon->exon transition matrix support coloring, exon `len%3` (genomic-length modulo 3) cues, and CDS flank phase edge coloring (`0/1/2`) when transcript `cds_ranges_1based` are available
- `RenderIsoformArchitectureSvg { seq_id, panel_id, path }`
- `RenderRnaStructureSvg { seq_id, path }`
- `RenderLineageSvg { path }`
- `RenderPoolGelSvg { inputs, path, ladders? }`
- `RenderProtocolCartoonSvg { protocol, path }`
- `RenderProtocolCartoonTemplateSvg { template_path, path }`
- `ValidateProtocolCartoonTemplate { template_path }`
- `RenderProtocolCartoonTemplateWithBindingsSvg { template_path, bindings_path, path }`
- `ExportProtocolCartoonTemplateJson { protocol, path }`
- `ExportDnaLadders { path, name_filter? }`
- `ExportRnaLadders { path, name_filter? }`
- `ExportPool { inputs, path, pool_id?, human_id? }`
- `ExportProcessRunBundle { path, run_id? }`
- `Digest { input, enzymes, output_prefix? }`
- `Ligation { inputs, circularize_if_possible, protocol, output_id?, output_prefix?, unique? }`
- `MergeContainers { inputs, output_prefix? }`
- `Pcr { template, forward_primer, reverse_primer, output_id?, unique? }`
- `PcrAdvanced { template, forward_primer, reverse_primer, output_id?, unique? }`
- `PcrMutagenesis { template, forward_primer, reverse_primer, mutations, output_id?, unique?, require_all_mutations? }`
- `DesignPrimerPairs { ... }` (implemented baseline)
- `DesignQpcrAssays { ... }` (implemented baseline; forward/reverse/probe)
- `ComputeDotplot { seq_id, reference_seq_id?, span_start_0based?, span_end_0based?, reference_span_start_0based?, reference_span_end_0based?, mode, word_size, step_bp, max_mismatches?, tile_bp?, store_as? }` (implemented baseline, self + pairwise)
- `ComputeFlexibilityTrack { seq_id, span_start_0based?, span_end_0based?, model, bin_bp, smoothing_bp?, store_as? }` (implemented baseline)
- `DeriveSplicingReferences { seq_id, span_start_0based, span_end_0based, seed_feature_id?, scope?, output_prefix? }` (implemented baseline; emits derived DNA window + mRNA isoforms + exon-reference sequence)
- `AlignSequences { query_seq_id, target_seq_id, query_span_start_0based?, query_span_end_0based?, target_span_start_0based?, target_span_end_0based?, mode?, match_score?, mismatch_score?, gap_open?, gap_extend? }` (implemented baseline; returns structured pairwise local/global report in `OpResult.sequence_alignment`)
- `InterpretRnaReads { seq_id, seed_feature_id, profile, input_path, input_format, scope, origin_mode?, target_gene_ids?, roi_seed_capture_enabled?, seed_filter, align_config, report_id?, report_mode?, checkpoint_path?, checkpoint_every_reads?, resume_from_checkpoint? }` (Nanopore cDNA phase-1 seed-filter pass; `multi_gene_sparse` expands local transcript-template indexing, while ROI capture remains planned)
- `AlignRnaReadReport { report_id, selection, align_config_override?, selected_record_indices? }` (Nanopore cDNA phase-2 retained-hit alignment pass; updates mapping/MSA/abundance report fields and re-ranks retained hits by alignment-aware retention rank)
- `ListRnaReadReports { seq_id? }`
- `ShowRnaReadReport { report_id }`
- `ExportRnaReadReport { report_id, path }`
- `ExportRnaReadHitsFasta { report_id, path, selection }`
- `ExportRnaReadSampleSheet { path, seq_id?, report_ids?, append? }`
- `ExportRnaReadExonPathsTsv { report_id, path, selection }`
- `ExportRnaReadExonAbundanceTsv { report_id, path, selection }`
- `ExportRnaReadScoreDensitySvg { report_id, path, scale }`
- `ExportRnaReadAlignmentsTsv { report_id, path, selection, limit? }`
- `ExportRnaReadAlignmentDotplotSvg { report_id, path, selection, max_points }`
- `ExtractRegion { input, from, to, output_id? }`
- `PrepareGenome { genome_id, catalog_path?, cache_dir?, timeout_seconds? }`
- `ExtractGenomeRegion { genome_id, chromosome, start_1based, end_1based, output_id?, annotation_scope?, max_annotation_features?, include_genomic_annotation?, catalog_path?, cache_dir? }`
  - `annotation_scope` accepts `none|core|full` and defaults to `core` when omitted.
  - `max_annotation_features` is an optional safety cap (0 or omitted = unlimited for explicit requests).
  - legacy `include_genomic_annotation` is still accepted (`true` -> `core`, `false` -> `none`) for compatibility.
  - operation results include `genome_annotation_projection` telemetry (requested/effective scope, feature counts, fallback metadata).
  - for helper genome IDs containing `pUC18`/`pUC19`, the engine applies a deterministic fallback MCS `misc_feature` annotation when source annotation does not already include an MCS feature and exactly one canonical MCS motif is found.
  - source-derived and fallback MCS features expose `mcs_expected_sites` with REBASE-normalized enzyme names when recognizable.
- `ExtractGenomeGene { genome_id, gene_query, occurrence?, output_id?, annotation_scope?, max_annotation_features?, include_genomic_annotation?, catalog_path?, cache_dir? }`
  - `annotation_scope` accepts `none|core|full` and defaults to `core` when omitted.
  - `max_annotation_features` is an optional safety cap (0 or omitted = unlimited for explicit requests).
  - legacy `include_genomic_annotation` is still accepted (`true` -> `core`, `false` -> `none`) for compatibility.
  - operation results include `genome_annotation_projection` telemetry (requested/effective scope, feature counts, fallback metadata).
  - for helper genome IDs containing `pUC18`/`pUC19`, the same deterministic MCS fallback annotation behavior applies when an MCS feature is missing; non-unique motif matches are warned and skipped.
- `ExtendGenomeAnchor { seq_id, side, length_bp, output_id?, catalog_path?, cache_dir?, prepared_genome_id? }`
- `VerifyGenomeAnchor { seq_id, catalog_path?, cache_dir?, prepared_genome_id? }`
- `ImportBlastHitsTrack { seq_id, hits[], track_name?, clear_existing?, blast_provenance? }`
  - optional `blast_provenance` payload preserves invocation context
    (`genome_id`, `query_label`, `query_length`, `max_hits`, `task`,
    `blastn_executable`, `blast_db_prefix`, raw `command[]`, `command_line`,
    `catalog_path?`, `cache_dir?`, `options_override_json?`,
    `effective_options_json?`) for sequence-history/audit views.
- `SelectCandidate { input, criterion, output_id? }`
- `ImportIsoformPanel { seq_id, panel_path, panel_id?, strict }`
- `ImportUniprotSwissProt { path, entry_id? }`
- `FetchUniprotSwissProt { query, entry_id? }`
- `ImportUniprotEntrySequence { entry_id, output_id? }`
  - currently returns `Unsupported`: protein sequence windows are deferred; use
    UniProt entries as metadata/projection sources in this release.
- `FetchGenBankAccession { accession, as_id? }`
- `ProjectUniprotToGenome { seq_id, entry_id, projection_id?, transcript_id? }`
- `GenerateCandidateSet { set_name, seq_id, length_bp, step_bp, feature_kinds[], feature_label_regex?, max_distance_bp?, feature_geometry_mode?, feature_boundary_mode?, feature_strand_relation?, limit? }`
- `GenerateCandidateSetBetweenAnchors { set_name, seq_id, anchor_a, anchor_b, length_bp, step_bp, limit? }`
- `DeleteCandidateSet { set_name }`
- `UpsertGuideSet { guide_set_id, guides[] }`
- `DeleteGuideSet { guide_set_id }`
- `FilterGuidesPractical { guide_set_id, config?, output_guide_set_id? }`
- `GenerateGuideOligos { guide_set_id, template_id, apply_5prime_g_extension?, output_oligo_set_id?, passed_only? }`
- `ExportGuideOligos { guide_set_id, oligo_set_id?, format: csv_table|plate_csv|fasta, path, plate_format? }`
- `ExportGuideProtocolText { guide_set_id, oligo_set_id?, path, include_qc_checklist? }`
- `ScoreCandidateSetExpression { set_name, metric, expression }`
- `ScoreCandidateSetDistance { set_name, metric, feature_kinds[], feature_label_regex?, feature_geometry_mode?, feature_boundary_mode?, feature_strand_relation? }`
- `FilterCandidateSet { input_set, output_set, metric, min?, max?, min_quantile?, max_quantile? }`
- `CandidateSetOp { op: union|intersect|subtract, left_set, right_set, output_set }`
- `ScoreCandidateSetWeightedObjective { set_name, metric, objectives[], normalize_metrics? }`
- `TopKCandidateSet { input_set, output_set, metric, k, direction?, tie_break? }`
- `ParetoFrontierCandidateSet { input_set, output_set, objectives[], max_candidates?, tie_break? }`
- `UpsertWorkflowMacroTemplate { name, description?, details_url?, parameters[], script }`
- `DeleteWorkflowMacroTemplate { name }`
- `UpsertCandidateMacroTemplate { name, description?, details_url?, parameters[], script }`
- `DeleteCandidateMacroTemplate { name }`
- `FilterByMolecularWeight { inputs, min_bp, max_bp, error, unique, output_prefix? }`
- `FilterByDesignConstraints { inputs, gc_min?, gc_max?, max_homopolymer_run?, reject_ambiguous_bases?, avoid_u6_terminator_tttt?, forbidden_motifs?, unique, output_prefix? }`
- `Reverse { input, output_id? }`
- `Complement { input, output_id? }`
- `ReverseComplement { input, output_id? }`
- `Branch { input, output_id? }`
- `SetDisplayVisibility { target, visible }`
- `SetLinearViewport { start_bp, span_bp }`
- `SetTopology { seq_id, circular }`
- `RecomputeFeatures { seq_id }`
- `SetParameter { name, value }` (purely in-silico project parameter change)

Isoform-panel operation semantics (current):

- `ImportIsoformPanel` loads curated panel resources with schema
  `gentle.isoform_panel_resource.v1` and binds them to one sequence context.
- `strict=true` enforces hard failure when panel transcript mapping fails;
  `strict=false` records warnings and keeps partial mappings.
- `RenderIsoformArchitectureSvg` emits a deterministic two-section architecture
  SVG (transcript/exon lanes + protein/domain lanes) derived from the same
  expert payload used by GUI/shell inspection.
  - when CDS ranges are available for mapped transcripts, SVG rendering uses
    dual coding in the top panel (faint full exons + solid CDS blocks) and
    adds a genome boundary rail with semi-transparent flank ribbons mapping
    boundary intervals to amino-acid spans on the shared protein reference axis
    (`1 aa ... max aa`), so mapping is readable across all protein lanes.
    Identical ribbons are merged and rendered once with support-weighted opacity.
- `gentle.isoform_panel_resource.v1` supports optional protein reference-span
  hints per isoform:
  - `reference_start_aa` (1-based inclusive)
  - `reference_end_aa` (1-based inclusive)
  - when present, protein lanes render and clip domains within this span while
    keeping one shared amino-acid axis across isoforms (useful for TP53
    N-terminus/C-terminus class overlays).
- `gentle.isoform_panel_resource.v1` also supports panel-level transcript
  geometry mode:
  - `transcript_geometry_mode: exon|cds` (default `exon`)
  - `cds` renders top-panel lanes from transcript CDS segments when available,
    falling back to exon geometry per transcript if CDS metadata is missing.

`LoadFile` import detection semantics (current):

- deterministic probe order: `GenBank -> EMBL -> FASTA -> XML`
- XML scope: `GBSet/GBSeq` is supported
- unsupported XML dialects (for example `INSDSet/INSDSeq`) return explicit
  schema/dialect diagnostics

`ExtendGenomeAnchor` side semantics:

- `side` accepts `five_prime` or `three_prime`.
- Direction is contextual to anchor strand.
- On anchor strand `-`, `five_prime` increases physical genomic position.
- If the anchor genome id is not prepared exactly, the engine can auto-resolve
  to one compatible prepared assembly-family entry (for example `GRCh38.p14`
  -> `Human GRCh38 Ensembl 116`).
- If multiple compatible prepared entries exist, extension fails with a
  deterministic options list so caller/GUI can choose explicitly.
- `prepared_genome_id` can be passed explicitly to force a specific prepared
  cache and bypass compatibility auto-selection.

`VerifyGenomeAnchor` semantics:

- Re-checks one anchored sequence against the selected prepared genome cache at
  recorded coordinates/strand.
- Writes one new provenance entry with `operation = VerifyGenomeAnchor` and
  `anchor_verified = true|false`.
- Returns an in-place state change (`changed_seq_ids`) for the same sequence id
  so GUI/CLI can refresh verification badges/status lines deterministically.

Local `SequenceAnchor` semantics (distinct from genome provenance anchoring):

- `SequenceAnchor` currently supports:
  - `Position { zero_based }`
  - `FeatureBoundary { feature_kind?, feature_label?, boundary, occurrence? }`
- `boundary` accepts `Start`, `End`, or `Middle`.
- This anchor model resolves in-sequence positions and is used for
  in-silico extraction/scoring workflows (`ExtractAnchoredRegion`,
  `GenerateCandidateSetBetweenAnchors`).

Adapter utility contracts (current, non-engine operations):

For narrative/operator guidance on when to use CLI, MCP, Agent Assistant, or an
external coding agent runtime, see:

- `docs/agent_interfaces_tutorial.md`

- `help [COMMAND ...] [--format text|json|markdown] [--interface ...]`
  - backed by structured glossary source `docs/glossary.json`
  - `--format text` renders human-readable help
  - `--format json` renders machine-readable help catalog/topic payload
  - `--format markdown` renders documentation-ready markdown
  - `--interface` accepts: `all|cli-direct|cli-shell|gui-shell|js|lua|mcp`
    (`mcp` currently aliases to shared shell command docs)
- shared-shell isoform panel routes:
  - `panels import-isoform SEQ_ID PANEL_PATH [--panel-id ID] [--strict]`
  - `panels inspect-isoform SEQ_ID PANEL_ID`
  - `panels render-isoform-svg SEQ_ID PANEL_ID OUTPUT.svg`
  - `panels validate-isoform PANEL_PATH [--panel-id ID]`
- shared-shell UniProt routes:
  - `uniprot fetch QUERY [--entry-id ID]`
  - `uniprot import-swissprot PATH [--entry-id ID]`
  - `uniprot list`
  - `uniprot show ENTRY_ID`
  - `uniprot map ENTRY_ID SEQ_ID [--projection-id ID] [--transcript ID]`
  - `uniprot projection-list [--seq SEQ_ID]`
  - `uniprot projection-show PROJECTION_ID`
- shared-shell GenBank route:
  - `genbank fetch ACCESSION [--as-id ID]`
- shared-shell protocol-cartoon routes:
  - `protocol-cartoon list`
  - `protocol-cartoon render-svg PROTOCOL_ID OUTPUT.svg`
  - `protocol-cartoon render-template-svg TEMPLATE.json OUTPUT.svg`
  - `protocol-cartoon template-validate TEMPLATE.json`
  - `protocol-cartoon render-with-bindings TEMPLATE.json BINDINGS.json OUTPUT.svg`
  - `protocol-cartoon template-export PROTOCOL_ID OUTPUT.json`
  - command surface is intentionally canonical: protocol-cartoon routes do not
    expose extra alias names

- Python adapter wrapper (`integrations/python/gentle_py`):
  - thin subprocess-based wrapper over `gentle_cli`
  - deterministic methods:
    - `capabilities()`
    - `state_summary()`
    - `op(operation)`
    - `workflow(workflow|workflow_path)`
    - `shell(line, expect_json=False)`
    - `render_dotplot_svg(seq_id, dotplot_id, output_svg, ...)`
  - raises structured `GentleCliError` with:
    - `code` (best-effort extracted stable code token)
    - `command`, `exit_code`, `stdout`, `stderr`
  - executable resolution order:
    1. constructor `cli_cmd`
    2. `GENTLE_CLI_CMD`
    3. `gentle_cli` on `PATH`
    4. repository fallback `cargo run --quiet --bin gentle_cli --`

- `gentle_mcp` (stdio MCP adapter, expanded UI-intent parity baseline)
  - MCP role:
    - request/response transport for tool execution (`tools/call`)
    - standardized capability discovery/negotiation (`tools/list`,
      `capabilities`, `help`)
  - current tools:
    - `capabilities`
    - `state_summary`
    - `op` (apply one `Operation`; requires explicit `confirm=true`)
    - `workflow` (apply one `Workflow`; requires explicit `confirm=true`)
    - `help`
    - `ui_intents` (shared `ui intents` catalog)
    - `ui_intent` (shared `ui open|focus ...` resolution path)
    - `ui_prepared_genomes` (shared `ui prepared-genomes ...` query path)
    - `ui_latest_prepared` (shared `ui latest-prepared ...` query path)
  - successful mutating calls (`op`, `workflow`) persist state to the resolved
    `state_path`
  - UI-intent tools route through the shared shell parser/executor
    (`parse_shell_tokens` + `execute_shell_command_with_options`) and are
    required to remain non-mutating (`state_changed = false`)
  - tool handlers are adapter wrappers over existing deterministic engine/shell
    contracts (no MCP-only biology logic branch)
  - stdio framing/validation hardening:
    - `Content-Length` is required, duplicate headers are rejected
    - maximum accepted frame size is `8 MiB`
    - parsed JSON nesting depth is capped at `96`
    - `tools/call` params are strict (`name`, optional `arguments` only)
    - `tools/call.arguments` must be a JSON object

MCP UI-intent tool contracts (current):

- `ui_intents`
  - arguments:
    - `state_path?` (optional; accepted for interface symmetry)
  - behavior:
    - executes shared shell command: `ui intents`
  - result:
    - structured payload schema: `gentle.ui_intents.v1`
    - includes stable `targets`, `commands`, and deterministic notes

- `ui_intent`
  - arguments:
    - required: `action` (`open|focus`), `target`
    - optional: `state_path`, `genome_id`, `helpers`, `catalog_path`,
      `cache_dir`, `filter`, `species`, `latest`
  - behavior:
    - executes shared shell command:
      - `ui open TARGET ...` or `ui focus TARGET ...`
    - for `target = prepared-references`, optional query flags can resolve
      `selected_genome_id` deterministically through the same helper path used
      by shared shell/CLI
    - parser guardrails are preserved:
      - query flags (`--helpers`, `--catalog`, `--cache-dir`, `--filter`,
        `--species`, `--latest`) are rejected for non-`prepared-references`
        targets
  - result:
    - structured payload schema: `gentle.ui_intent.v1`
    - fields include `ui_intent`, `selected_genome_id`, optional
      `prepared_query`, `applied=false`, and deterministic `message`

- `ui_prepared_genomes`
  - arguments:
    - optional: `state_path`, `helpers`, `catalog_path`, `cache_dir`, `filter`,
      `species`, `latest`
  - behavior:
    - executes shared shell command: `ui prepared-genomes ...`
  - result:
    - structured payload schema: `gentle.ui_prepared_genomes.v1`
    - includes `prepared_count`, sorted `genomes[]`, and `selected_genome_id`

- `ui_latest_prepared`
  - arguments:
    - required: `species`
    - optional: `state_path`, `helpers`, `catalog_path`, `cache_dir`
  - behavior:
    - executes shared shell command: `ui latest-prepared SPECIES ...`
  - result:
    - structured payload schema: `gentle.ui_latest_prepared.v1`
    - includes `selected_genome_id` and nested `prepared_query` payload

MCP UI-intent JSON-RPC example (abbreviated):

```json
{
  "jsonrpc": "2.0",
  "id": 7,
  "method": "tools/call",
  "params": {
    "name": "ui_intent",
    "arguments": {
      "action": "open",
      "target": "prepared-references",
      "catalog_path": "assets/genomes.json",
      "species": "human",
      "latest": true
    }
  }
}
```

Result envelope shape:

```json
{
  "jsonrpc": "2.0",
  "id": 7,
  "result": {
    "isError": false,
    "structuredContent": {
      "schema": "gentle.ui_intent.v1",
      "selected_genome_id": "Human GRCh38 Ensembl 116",
      "applied": false
    }
  }
}
```

Adapter-equivalence guarantee for UI-intent tools:

- deterministic parity tests compare MCP UI-intent tool outputs with direct
  shared shell `ui ...` command outputs for:
  - intent catalog (`ui_intents`)
  - prepared query (`ui_prepared_genomes`)
  - latest helper (`ui_latest_prepared`)
  - open/focus intent resolution (`ui_intent`)

- `macros run/instance-list/instance-show/template-list/template-show/template-put/template-delete/template-import/template-run`
  - shared-shell macro adapter family for full operation/workflow scripting
  - template persistence is backed by engine operations
    `UpsertWorkflowMacroTemplate`/`DeleteWorkflowMacroTemplate`
  - `template-put` supports optional typed port contracts:
    - `--input-port PORT_ID:KIND[:one|many][:required|optional][:description]`
    - `--output-port PORT_ID:KIND[:one|many][:required|optional][:description]`
  - `template-import PATH` accepts:
    - one pack JSON file (`gentle.cloning_patterns.v1`)
    - one single-template JSON file (`gentle.cloning_pattern_template.v1`)
    - one directory tree (recursive `*.json` import; files must use one of the
      schemas above)
  - imports are transactional; if one template fails validation, no imported
    template changes are kept
  - expanded scripts can execute `op ...` and `workflow ...` statements and
    optionally roll back via `--transactional`
  - `template-run` supports non-mutating preflight mode via `--validate-only`
  - template-run responses now include a preflight payload
    (`gentle.macro_template_preflight.v1`) with warnings/errors and typed
    input/output port validation rows (`contract_source` indicates whether
    checks came from template metadata or routine catalog)
  - preflight includes cross-port semantic checks (alias/collision checks,
    input sequence/container consistency, and sequence-anchor semantics when
    sequence context is unambiguous)
  - routine-family semantic checks are now supported:
    - Gibson routines validate adjacent fragment overlap compatibility against
      configured overlap length before execution
    - Restriction routines validate enzyme-name resolution, duplicate-enzyme
      misuse, enzyme-site presence across bound input sequences, and common
      digest parameter sanity (`left_fragment`/`right_fragment`,
      `extract_from`/`extract_to`)
  - mutating `macros run` / `macros template-run` executions always persist one
    lineage macro-instance record (`ok`/`failed`/`cancelled`)
  - successful runs return `macro_instance_id`; failed runs include
    `macro_instance_id=...` in error messages
  - `macros instance-list` and `macros instance-show` expose persisted lineage
    macro-instance records as first-class introspection contracts

- `routines list [--catalog PATH] [--family NAME] [--status NAME] [--tag TAG] [--query TEXT]`
  - shared-shell/CLI routine catalog discovery surface
  - default catalog path: `assets/cloning_routines.json`
  - typed catalog schema: `gentle.cloning_routines.v1`
  - response schema: `gentle.cloning_routines_list.v1`
  - filters are case-insensitive; query performs substring match across
    routine id/title/family/status/template/tags/summary plus explainability
    metadata fields
- `routines explain ROUTINE_ID [--catalog PATH]`
  - shared-shell/CLI routine explainability surface
  - response schema: `gentle.cloning_routine_explain.v1`
  - returns one routine definition plus normalized explanation payload
    (purpose/mechanism/requires/contraindications/disambiguation/failure modes)
    and resolved confusing alternatives
- `routines compare ROUTINE_A ROUTINE_B [--catalog PATH]`
  - shared-shell/CLI deterministic routine comparison surface
  - response schema: `gentle.cloning_routine_compare.v1`
  - returns both routine definitions plus comparison payload:
    shared/unique tags, cross-reference status, aligned difference-matrix rows,
    and merged disambiguation questions
  - includes planning-aware estimate rows in comparison payload:
    - `estimated_time_hours`
    - `estimated_cost`
    - `local_fit_score`
    - `composite_meta_score`

- Planning meta-layer contracts (shared shell/CLI, engine-owned):
  - profile schema: `gentle.planning_profile.v1`
  - objective schema: `gentle.planning_objective.v1`
  - estimate schema: `gentle.planning_estimate.v1`
  - suggestion schema: `gentle.planning_suggestion.v1`
  - sync-status schema: `gentle.planning_sync_status.v1`
  - merge precedence for effective profile:
    - `global_profile -> confirmed_agent_overlay -> project_override`
  - purchasing latency heuristic in v1:
    - each missing required material class adds default
      `procurement_business_days_default` (default `10`) to estimate
      (Monday-Friday business-day model; no holiday calendar yet)
    - business-day delays are converted to `estimated_time_hours` with a
      deterministic weekend-aware factor (`24h * 7/5` per business day)
  - schema compatibility rule:
    - profile/objective payloads with mismatched schema ids are rejected
      (`InvalidInput`) instead of silently coerced
- `planning profile show [--scope global|project_override|confirmed_agent_overlay|effective]`
  - inspect one planning profile scope or merged effective profile
- `planning profile set JSON_OR_@FILE [--scope global|project_override|confirmed_agent_overlay]`
  - set/replace selected planning profile scope
- `planning profile clear [--scope global|project_override|confirmed_agent_overlay]`
  - clear selected planning profile scope
- `planning objective show`
  - inspect current planning objective
- `planning objective set JSON_OR_@FILE`
  - set/replace planning objective
- `planning objective clear`
  - clear planning objective (engine defaults apply)
- `planning suggestions list [--status pending|accepted|rejected]`
  - list pending/resolved planning sync suggestions
- `planning suggestions accept SUGGESTION_ID`
  - accept suggestion and apply patch into confirmed overlay/objective
- `planning suggestions reject SUGGESTION_ID [--reason TEXT]`
  - reject suggestion with optional reason
- `planning sync status`
  - inspect planning sync lifecycle metadata
- `planning sync pull JSON_OR_@FILE [--source ID] [--confidence N] [--snapshot-id ID]`
  - register inbound advisory suggestion as pending
- `planning sync push JSON_OR_@FILE [--source ID] [--confidence N] [--snapshot-id ID]`
  - register outbound advisory suggestion as pending
  - payload for `planning sync pull|push`:
    - optional `profile_patch` (`gentle.planning_profile.v1`)
    - optional `objective_patch` (`gentle.planning_objective.v1`)
    - optional `message`
  - activation policy remains explicit user action (`accept`/`reject`);
    no auto-apply in v1

- `screenshot-window OUTPUT.png`
  - currently disabled by security policy
  - returns deterministic disabled message from shared shell/CLI/GUI command
    paths
  - kept as reserved adapter contract for future re-enable after explicit
    endpoint-security approval

- `agents list [--catalog PATH]`
  - Lists configured agent systems from catalog JSON.
  - Default catalog: `assets/agent_systems.json`.

- `agents ask SYSTEM_ID --prompt TEXT [--catalog PATH] [--base-url URL] [--model MODEL] [--timeout-secs N] [--connect-timeout-secs N] [--read-timeout-secs N] [--max-retries N] [--max-response-bytes N] [--allow-auto-exec] [--execute-all] [--execute-index N ...] [--no-state-summary]`
  - Invokes one configured agent system via catalog transport.
  - `--base-url` applies a per-request runtime base URL override for native
    transports (`native_openai`, `native_openai_compat`).
  - `--model` applies a per-request runtime model override for native
    transports (`native_openai`, `native_openai_compat`).
  - `--timeout-secs` applies a per-request timeout override for stdio/native
    transports (maps to `GENTLE_AGENT_TIMEOUT_SECS`).
  - `--connect-timeout-secs` applies a per-request HTTP connect timeout override
    for native transports (maps to `GENTLE_AGENT_CONNECT_TIMEOUT_SECS`).
  - `--read-timeout-secs` applies a per-request read timeout override for
    stdio/native transports (maps to `GENTLE_AGENT_READ_TIMEOUT_SECS`).
  - `--max-retries` applies a per-request transient retry budget override
    (maps to `GENTLE_AGENT_MAX_RETRIES`; `0` disables retries).
  - `--max-response-bytes` applies a per-request response body/output cap
    override (maps to `GENTLE_AGENT_MAX_RESPONSE_BYTES`).
  - `--no-state-summary` suppresses project context injection.
  - Suggested-command execution is per-suggestion only (no global always-execute).

Agent bridge catalog schema (`gentle.agent_systems.v1`):

```json
{
  "schema": "gentle.agent_systems.v1",
  "systems": [
    {
      "id": "openai_gpt5_stdio",
      "label": "OpenAI GPT-5 (stdio bridge)",
      "description": "Optional human-readable description",
      "transport": "external_json_stdio",
      "command": ["openai-agent-bridge", "--model", "gpt-5"],
      "env": {},
      "working_dir": null
    },
    {
      "id": "openai_gpt5_native",
      "label": "OpenAI GPT-5 (native HTTP)",
      "transport": "native_openai",
      "model": "gpt-5",
      "base_url": "https://api.openai.com/v1",
      "env": {}
    }
  ]
}
```

Transport notes:

- `builtin_echo`: offline/demo transport.
- `external_json_stdio`: requires local bridge executable from `command[0]`.
- `native_openai`: built-in OpenAI HTTP adapter; requires `OPENAI_API_KEY`
  (environment or system-level `env` override in catalog entry).
- `native_openai_compat`: built-in OpenAI-compatible local HTTP adapter
  (`/chat/completions`), intended for local services such as Jan/Msty/Ollama
  when they expose an OpenAI-compatible endpoint. API key is optional.
- `GENTLE_AGENT_BASE_URL` (or CLI `--base-url`) overrides catalog `base_url`
  per request for `native_openai` and `native_openai_compat`.
- `GENTLE_AGENT_MODEL` (or CLI `--model`) overrides catalog `model` per request
  for `native_openai` and `native_openai_compat`.
- `GENTLE_AGENT_TIMEOUT_SECS` (or CLI `--timeout-secs`) overrides request
  timeout per attempt for agent transports.
- `GENTLE_AGENT_CONNECT_TIMEOUT_SECS` (or CLI `--connect-timeout-secs`)
  overrides HTTP connect timeout for native transports.
- `GENTLE_AGENT_READ_TIMEOUT_SECS` (or CLI `--read-timeout-secs`) overrides
  read timeout for stdio/native transports.
- `GENTLE_AGENT_MAX_RETRIES` (or CLI `--max-retries`) overrides transient retry
  count (`0` disables retries).
- `GENTLE_AGENT_MAX_RESPONSE_BYTES` (or CLI `--max-response-bytes`) overrides
  response-size cap per attempt (stdout/stderr or HTTP body).
- `native_openai_compat` requires a concrete model name; value `unspecified`
  is treated as missing and the request is rejected until a model is provided.
- `native_openai_compat` does not silently switch host/port; it uses catalog
  `base_url` or explicit `GENTLE_AGENT_BASE_URL`.

Agent request payload schema (`gentle.agent_request.v1`):

```json
{
  "schema": "gentle.agent_request.v1",
  "system_id": "openai_gpt5_stdio",
  "prompt": "User request text",
  "sent_at_unix_ms": 1768860000000,
  "state_summary": {}
}
```

Agent response payload schema (`gentle.agent_response.v1`):

```json
{
  "schema": "gentle.agent_response.v1",
  "assistant_message": "Text response",
  "questions": ["Optional follow-up question"],
  "suggested_commands": [
    {
      "title": "Optional short label",
      "rationale": "Optional reason",
      "command": "state-summary",
      "execution": "ask"
    }
  ]
}
```

Agent execution intent semantics:

- `chat`: explain/ask only, never executed as shell command.
- `ask`: executable suggestion requiring explicit user confirmation.
- `auto`: executable suggestion eligible for automatic execution only when
  caller enables `--allow-auto-exec`.

Agent schema/compatibility policy:

- `schema` is mandatory for catalog/request/response JSON objects.
- Supported major versions (current): `gentle.agent_systems.v1`,
  `gentle.agent_request.v1`, `gentle.agent_response.v1`.
- Future incompatible major versions (for example `.v2`) are rejected with a
  deterministic schema-unsupported error.
- Response validation is strict for canonical fields:
  - top-level allowed: `schema`, `assistant_message`, `questions`,
    `suggested_commands` plus extension keys prefixed with `x_` or `x-`
  - `suggested_commands[]` allowed: `title`, `rationale`, `command`,
    `execution` plus extension keys prefixed with `x_` or `x-`
  - unsupported canonical fields (for example `commands`, `mode`) are rejected

Execution safety rules:

- There is no global always-execute mode.
- Execution is per suggestion:
  - explicit run (`--execute-index`, `--execute-all`, GUI row `Run`)
  - optional auto-run only for `execution = auto` + `--allow-auto-exec`
- Recursive `agents ask` execution from suggested commands is blocked.

Failure-handling policy for external adapters:

- Adapter invocations use bounded retry with exponential backoff for transient
  failures.
- OpenAI `429` with `insufficient_quota` is treated as non-transient (no retry)
  and returned with the original API error body plus billing/usage guidance.
- Missing/unreachable adapter binaries fail gracefully with deterministic
  adapter-unavailable errors.
- CLI/shell errors are stable and prefixed for scripting, e.g.:
  - `AGENT_INVALID_INPUT`
  - `AGENT_SCHEMA_VALIDATION`
  - `AGENT_SCHEMA_UNSUPPORTED`
  - `AGENT_ADAPTER_UNAVAILABLE`
  - `AGENT_ADAPTER_TRANSIENT`
  - `AGENT_ADAPTER_FAILED`
  - `AGENT_RESPONSE_PARSE`
  - `AGENT_RESPONSE_VALIDATION`

ClawBio/OpenClaw integration scaffold schemas:

- integration path:
  `integrations/clawbio/skills/gentle-cloning/`
- wrapper request schema: `gentle.clawbio_skill_request.v1`
  - `mode`: `capabilities|state-summary|shell|op|workflow|raw`
  - optional: `state_path`, `timeout_secs`
  - mode-specific:
    - `shell`: `shell_line`
    - `op`: `operation` (JSON object/string)
    - `workflow`: `workflow` or `workflow_path`
    - `raw`: `raw_args[]`
- wrapper result schema: `gentle.clawbio_skill_result.v1`
  - `status`: `ok|command_failed|timeout|failed|degraded_demo`
  - includes resolver details, executed command, exit code, stdout/stderr, and
    generated artifact paths
- reproducibility outputs:
  - `report.md`
  - `result.json`
  - `reproducibility/commands.sh`
  - `reproducibility/environment.yml`
  - `reproducibility/checksums.sha256`

Planned operation refinements:

- `MergeContainers { inputs, output_prefix? }`
  - Explicitly models wet-lab mixing of multiple tubes/pools.
- Protocol-based ligation:
  - `Ligation { input_container, protocol, output_container?, ... }`
  - `protocol` determines allowed end joins.
  - Initial protocol values:
    - `sticky`
    - `blunt`
  - Future protocol values may include established ligation workflows
    represented as named presets.

Current parameter support:

- `max_fragments_per_container` (default `80000`)
  - limits digest fragment output per operation
  - also serves as ligation product-count limit guard
- `require_verified_genome_anchor_for_extension` (default `false`)
  - when `true`, `ExtendGenomeAnchor` requires anchor provenance with
    `anchor_verified=true`
  - anchors with `anchor_verified=false` or missing verification status are
    rejected in strict mode
  - alias parameters accepted: `strict_genome_anchor_verification`,
    `strict_anchor_verification`
- `genome_anchor_prepared_fallback_policy` (default `single_compatible`)
  - controls how `ExtendGenomeAnchor` / `VerifyGenomeAnchor` resolve anchor
    genome ids when exact prepared cache id is not present.
  - accepted values:
    - `off` (no compatibility fallback; must match exact prepared id)
    - `single_compatible` (auto-fallback only when one compatible prepared
      cache exists)
    - `always_explicit` (never auto-fallback; require explicit selection even
      when only one compatible prepared cache exists)
  - alias parameters accepted: `genome_anchor_fallback_mode`,
    `genome_anchor_prepared_mode`
- primer-design backend controls:
  - `primer_design_backend` (default `auto`)
    - accepted values: `auto`, `internal`, `primer3`
    - `auto` tries Primer3 and falls back deterministically to internal scoring
      with explicit warning + fallback reason in report metadata
  - `primer3_executable` (default `"primer3_core"`)
    - executable path/name used when backend is `primer3` or `auto`
    - alias parameters accepted: `primer3_backend_executable`, `primer3_path`
- `feature_details_font_size` (default `9.0`, range `8.0..24.0`)
  - controls GUI font size for the feature tree entries and feature range details
- `regulatory_feature_max_view_span_bp` (default `50000`, range `>= 0`)
  - hides regulatory feature overlays in linear view when current view span
    exceeds this threshold (`0` disables regulatory overlays)
- `gc_content_bin_size_bp` (default `100`, range `>= 1`)
  - controls GC-content aggregation bin size for linear/circular rendering and
    SVG export
- Linear DNA-letter routing parameters:
  - `linear_sequence_letter_layout_mode` (default `AutoAdaptive`)
    - supported canonical modes:
      - `auto|adaptive|auto_adaptive`
      - `standard|standard_linear`
      - `helical|continuous_helical`
      - `condensed_10_row|condensed`
    - auto mode uses deterministic viewport-density tiers:
      - `<= 1.5x`: standard
      - `<= 2x`: helical (if compressed letters enabled)
      - `<= 10x`: condensed-10 (if compressed letters enabled)
      - `> 10x`: `OFF`
  - `linear_sequence_helical_letters_enabled` (default `true`)
    - applies to auto mode only (allows/disallows compressed auto tiers)
  - `linear_sequence_helical_phase_offset_bp` (range `0..9`)
    - seam offset used by helical/condensed row mapping
  - reverse/helical strand geometry controls:
    - `linear_show_double_strand_bases` / `linear_show_reverse_strand_bases`
      (bool alias pair; controls reverse-strand letter visibility)
    - `linear_helical_parallel_strands` (default `true`)
      - `true`: forward/reverse helical slant stays parallel
      - `false`: forward/reverse helical slant is mirrored (cross-over look)
    - `reverse_strand_visual_opacity` (range `0.2..1.0`, default `0.55`)
      - shared reverse-strand emphasis in linear map and sequence panel
- Legacy linear-letter threshold knobs are compatibility-only and return
  deterministic deprecated no-op messages (no routing effect):
  - `linear_sequence_base_text_max_view_span_bp`
  - `linear_sequence_helical_max_view_span_bp`
  - `linear_sequence_condensed_max_view_span_bp`
- VCF display filter parameters (shared GUI/SVG state):
  - `vcf_display_show_snp`
  - `vcf_display_show_ins`
  - `vcf_display_show_del`
  - `vcf_display_show_sv`
  - `vcf_display_show_other`
  - `vcf_display_pass_only`
  - `vcf_display_use_min_qual`
  - `vcf_display_min_qual`
  - `vcf_display_use_max_qual`
  - `vcf_display_max_qual`
  - `vcf_display_required_info_keys` (CSV string or string array)
- BLAST options-layer parameters:
  - `blast_options_override` (JSON object or `null`)
    - project-level BLAST option layer merged before per-command request JSON
    - supports the same keys as request JSON (`task`, `max_hits`, `thresholds`)
  - `blast_options_defaults_path` (string path or `null`)
    - optional defaults-file path used ahead of project/request layers
    - if unset, engine falls back to `assets/blast_defaults.json`

Current ligation protocol behavior:

- `protocol` is mandatory.
- If `protocol = Blunt`, ligation enumerates ordered input pairs with blunt-end
  compatibility checks.
- If `protocol = Sticky`, ligation enumerates ordered input pairs with sticky-end
  overhang compatibility checks.
- `unique = true` requires exactly one product.

`FilterByMolecularWeight` semantics:

- Applies a bp-range filter across provided input sequence ids.
- Effective accepted range is expanded by `error`:
  - `effective_min = floor(min_bp * (1 - error))`
  - `effective_max = ceil(max_bp * (1 + error))`
- `unique = true` requires exactly one match, otherwise the operation fails.

`FilterByDesignConstraints` semantics:

- Applies practical design-constraint filters across provided input sequence ids.
- Optional GC bounds:
  - `gc_min` and/or `gc_max` (fractional range `0.0..1.0`)
  - when both are provided, `gc_min <= gc_max` is required
- Optional homopolymer cap:
  - `max_homopolymer_run >= 1`
  - rejects candidates with a longer A/C/G/T run
- `reject_ambiguous_bases` (default `true`):
  - rejects sequences containing non-ACGT letters
- `avoid_u6_terminator_tttt` (default `true`):
  - rejects sequences containing `TTTT`
- Optional `forbidden_motifs`:
  - IUPAC motifs; reject when motif appears on either strand
- `unique = true` requires exactly one match, otherwise the operation fails.

Guide-design semantics:

- Guide sets persist in `ProjectState.metadata["guide_design"]`
  (`schema = gentle.guide_design.v1`) and include:
  - guide sets
  - practical-filter reports
  - oligo sets
  - audit log entries for guide operations/exports
- `UpsertGuideSet`:
  - normalizes guide fields and validates required properties
  - sorts by rank (then guide id) and rejects duplicate `guide_id` within one set
- `FilterGuidesPractical`:
  - applies deterministic practical filters over one guide set
  - supports GC bounds, global/per-base homopolymer limits, ambiguous-base
    rejection, U6 `TTTT` avoidance, dinucleotide repeat cap, forbidden motifs,
    and required 5' base checks
  - can emit a passed-only output guide set (`output_guide_set_id`)
  - always persists a structured per-guide report with reasons/warnings/metrics
- `GenerateGuideOligos`:
  - generates forward/reverse oligos using a named template
  - supports optional 5' G extension and passed-only mode
  - persists generated oligo records in named oligo sets
- `ExportGuideOligos`:
  - exports an oligo set as `csv_table`, `plate_csv` (96/384), or `fasta`
  - records export actions in the guide-design audit log
- `ExportGuideProtocolText`:
  - exports a deterministic human-readable protocol text artifact
  - optional QC checklist can be included/excluded

Candidate-set semantics:

- `GenerateCandidateSet` creates a persisted candidate window set over one source
  sequence and computes baseline metrics for each candidate.
- `GenerateCandidateSetBetweenAnchors` creates a persisted candidate window set
  constrained to the in-sequence interval between two local anchors.
- `ScoreCandidateSetExpression` computes a derived metric from an arithmetic
  expression over existing metrics.
- `ScoreCandidateSetDistance` computes feature-distance metrics against filtered
  feature targets.
- `FilterCandidateSet` keeps/drops candidates by absolute bounds and/or quantile
  bounds for a named metric.
- `CandidateSetOp` supports set algebra (`union`, `intersect`, `subtract`) over
  candidate identity (`seq_id`, `start_0based`, `end_0based`).
- `ScoreCandidateSetWeightedObjective` computes one metric from weighted
  objective terms (`maximize`/`minimize` per term, optional normalization).
- `TopKCandidateSet` selects an explicit top-k subset for one metric with a
  deterministic tie-break policy.
- `ParetoFrontierCandidateSet` keeps non-dominated candidates for multiple
  objectives (`maximize`/`minimize` per objective), with optional tie-break
  truncation.
- Workflow macro templates are persisted in project metadata:
  - `UpsertWorkflowMacroTemplate` stores/replaces named templates
  - `DeleteWorkflowMacroTemplate` removes templates
  - each template now carries `template_schema`
    (`gentle.cloning_macro_template.v1`) so cloning-operation macro intent is
    explicit at engine level
  - optional `details_url` can link to external protocol/reference material
  - optional typed `input_ports`/`output_ports` can be persisted directly in
    template metadata (same port shape as routine catalog ports)
  - template expansion/binding is exposed through adapter command surfaces
    (`macros template-*`, including `macros template-import PATH`)
  - expanded scripts run through shared shell execution (`macros run`) and can
    orchestrate full cloning operations via `op ...` or `workflow ...` payloads
  - shipped starter assets:
    - legacy pack: `assets/cloning_patterns.json` (`gentle.cloning_patterns.v1`)
    - hierarchical catalog: `assets/cloning_patterns_catalog/**/*.json`
      (`gentle.cloning_pattern_template.v1`, one template per file)
    - Gibson baseline template:
      `assets/cloning_patterns_catalog/gibson/overlap_assembly/gibson_two_fragment_overlap_preview.json`
    - Restriction baseline template:
      `assets/cloning_patterns_catalog/restriction/digest_ligation/digest_ligate_extract_sticky.json`
- Typed cloning-routine catalog baseline:
  - manifest: `assets/cloning_routines.json`
  - schema: `gentle.cloning_routines.v1`
  - typed routine metadata fields include routine family/status/tags, linked
    template name/path, and typed input/output port declarations
  - includes Gibson + restriction family baselines:
    - `gibson.two_fragment_overlap_preview`
    - `restriction.digest_ligate_extract_sticky`
  - adapter discovery surface:
    `routines list [--catalog PATH] [--family NAME] [--status NAME] [--tag TAG] [--query TEXT]`
  - explainability and comparison surfaces:
    - `routines explain ROUTINE_ID [--catalog PATH]`
    - `routines compare ROUTINE_A ROUTINE_B [--catalog PATH]`
- Macro-instance lineage baseline:
  - mutating `macros run` / `macros template-run` append one
    `LineageMacroInstance` record in project lineage state for success and
    failure pathways
  - records include deterministic `macro_instance_id`, optional
    `routine_id/template_name`, typed bound inputs/outputs, emitted `op_id`s,
    status, and optional `status_message`
  - lineage graph + lineage SVG consume these records as macro box nodes with
    explicit input/output edges where sequence/container references resolve
- Candidate macro templates are persisted in project metadata:
  - `UpsertCandidateMacroTemplate` stores/replaces named templates
  - `DeleteCandidateMacroTemplate` removes templates
  - optional `details_url` can link to external protocol/reference material
  - template expansion/binding is exposed through adapter command surfaces
    (`candidates template-*`)
- Between-anchor generation augments baseline metrics with anchor-aware fields
  (`distance_to_anchor_a_bp`, `distance_to_anchor_b_bp`,
  `distance_to_nearest_anchor_bp`, interval span metadata).

Feature-distance geometry controls (candidate generation and distance scoring):

- `feature_geometry_mode` (optional, default `feature_span`):
  - `feature_span`: one interval per feature using whole-feature bounds
  - `feature_parts`: one interval per explicit location part (ignores intronic
    gaps for multipart features)
  - `feature_boundaries`: boundary points of explicit location parts
- `feature_boundary_mode` (optional, default `any`):
  - `any`, `five_prime`, `three_prime`, `start`, `end`
  - only meaningful when `feature_geometry_mode = feature_boundaries`
- `feature_strand_relation` (optional, default `any`):
  - `any`, `same`, `opposite`
  - current engine interpretation is sequence-forward relative
    (`same = '+'`, `opposite = '-'` feature strand)
- Directed-boundary interpretation:
  - on `+` strand: `five_prime = start`, `three_prime = end`
  - on `-` strand: `five_prime = end`, `three_prime = start`
  - for unknown strand, `five_prime`/`three_prime` conservatively include both
    boundaries.

`RenderPoolGelSvg` semantics:

- Accepts explicit `inputs` (sequence ids) and an output `path`.
- Computes pool molecular-weight proxy from sequence bp lengths.
- Chooses one or two ladders to span pool range:
  - from explicit `ladders` list when provided
  - otherwise from built-in ladder catalog (auto mode)
- Renders ladder lanes plus pooled band lane as SVG artifact.

`RenderDotplotSvg` semantics:

- Inputs:
  - `seq_id` (owner/query sequence id for the stored dotplot payload)
  - `dotplot_id` (stored payload id from `ComputeDotplot`)
  - `path` (output SVG)
  - optional `flex_track_id` (adds flexibility panel in same SVG)
  - optional `display_density_threshold` and `display_intensity_gain` (display tuning)
- Ownership checks:
  - dotplot payload must belong to `seq_id`
  - optional flexibility track must also belong to `seq_id`
- Output:
  - deterministic SVG dotplot artifact; operation is non-mutating.

`RenderProtocolCartoonSvg` semantics:

- Inputs:
  - `protocol` (currently supported baseline: `gibson.two_fragment`)
  - `path` (output SVG)
- Behavior:
  - renders a deterministic protocol-cartoon strip through one engine route,
    independent of GUI/CLI entry point.
  - emits canonical conceptual step order for the requested protocol as an
    ordered event-sequence model.
  - template representation baseline is now available in engine internals:
    - schema id: `gentle.protocol_cartoon_template.v1`
    - sparse template rows (event/molecule/feature) are resolved with
      deterministic defaults into render-ready specs.
  - internal model used by renderer:
    - event -> molecules -> feature fragments
    - molecule topology supports `linear|circular`
    - linear molecules may carry end styles
      (`NotShown`, `Continuation`, `Blunt`, or
      `Sticky { polarity: FivePrime|ThreePrime, nt }`)
    - feature fragments can optionally render different top-strand and
      bottom-strand colors and lengths, plus strand-specific nicks after a
      segment boundary; this is useful for annealed overlaps, exonuclease
      chew-back cartoons with single-stranded tails, and polymerase-filled
      intermediates that still require ligase
  - malformed protocol cartoon specs fail validation and render deterministic
    invalid-spec SVG diagnostics instead of panicking.
- Output:
  - deterministic SVG artifact; operation is non-mutating.

`RenderProtocolCartoonTemplateSvg` semantics:

- Inputs:
  - `template_path` (JSON file path, schema `gentle.protocol_cartoon_template.v1`)
  - `path` (output SVG)
- Behavior:
  - reads template JSON from disk and parses it deterministically.
  - resolves sparse event/molecule/feature rows using deterministic defaults
    (action/caption/topology/end styles/feature length/palette).
  - validates resolved cartoon semantics before rendering.
- Output:
  - deterministic SVG artifact; operation is non-mutating.

`ValidateProtocolCartoonTemplate` semantics:

- Inputs:
  - `template_path` (JSON file path, schema `gentle.protocol_cartoon_template.v1`)
- Behavior:
  - reads and parses template JSON deterministically.
  - resolves sparse defaults and validates resolved cartoon semantics.
  - emits validation diagnostics through operation result messages; no SVG is
    written.
- Output:
  - non-mutating validation result suitable for pre-render checks in CLI/GUI
    flows.

`RenderProtocolCartoonTemplateWithBindingsSvg` semantics:

- Inputs:
  - `template_path` (JSON file path, schema `gentle.protocol_cartoon_template.v1`)
  - `bindings_path` (JSON file path, schema
    `gentle.protocol_cartoon_template_bindings.v1`)
  - `path` (output SVG)
- Behavior:
  - loads template and binding payloads.
  - applies deterministic ID-targeted overrides (defaults, event, molecule,
    feature) and then resolves the bound template.
  - validates resolved semantics before SVG rendering.
- Output:
  - deterministic SVG artifact; operation is non-mutating.

`ExportProtocolCartoonTemplateJson` semantics:

- Inputs:
  - `protocol` (built-in protocol cartoon id, for example `gibson.two_fragment`)
  - `path` (output JSON file)
- Behavior:
  - materializes the canonical built-in template
    (`gentle.protocol_cartoon_template.v1`) for the requested protocol.
  - writes deterministic pretty JSON suitable for user editing/tweaking.
- Output:
  - deterministic JSON artifact; operation is non-mutating.

`ExportProcessRunBundle` semantics:

- Exports a deterministic JSON run bundle artifact (`gentle.process_run_bundle.v1`)
  for reproducibility/audit.
- Inputs:
  - `path` (required): output JSON file
  - `run_id` (optional): when set, only operation-log rows for that `run_id`
    are exported; when omitted, all operation-log rows are exported.
- Payload sections:
  - `inputs`:
    - per-operation extracted input references
      (`sequence_ids`, `container_ids`, `arrangement_ids`, candidate/guide sets,
      genome ids, file inputs)
    - aggregated referenced ids and inferred `root_sequence_ids`
  - `parameter_overrides`:
    - chronological `SetParameter` overrides with `op_id`, `record_index`,
      parameter `name`, and exact JSON `value`
  - `decision_traces`:
    - optional routine-assistant trace rows
      (`gentle.routine_decision_trace.v1`) captured in project metadata and
      exported for routine-selection reproducibility (`trace_id`, selected
      routine/alternatives, disambiguation questions/answers, binding snapshot,
      ordered `preflight_history`, canonical `preflight_snapshot`, execution
      outcome, export events)
  - `operation_log`:
    - selected immutable operation records (`run_id`, operation payload, result)
  - `outputs`:
    - created/changed sequence ids
    - final sequence summaries for affected ids
    - container/arrangement ids created by selected operations
    - file artifact paths produced by selected operations
  - `parameter_snapshot`:
    - full current engine parameter snapshot at export time.
- Failure modes:
  - empty `path` => `InvalidInput`
  - unknown filtered `run_id` (no selected rows) => `NotFound`

RNA secondary-structure semantics:

- Inspection API:
  - `GentleEngine::inspect_rna_structure(seq_id)`
  - Runs `rnapkin -v -p <sequence>` and returns structured text report (`stdout`/`stderr` + command metadata).
- Export operation:
  - `RenderRnaStructureSvg { seq_id, path }`
  - Runs `rnapkin <sequence> <path>` and expects SVG output at `path`.
- Input constraints:
  - accepted only for single-stranded RNA (`molecule_type` `RNA` or `ssRNA`)
  - empty sequence is rejected
- Runtime dependency:
  - external `rnapkin` executable is required
  - executable path resolution order:
    1. env var `GENTLE_RNAPKIN_BIN`
    2. fallback executable name `rnapkin` in `PATH`

DNA ladder catalog semantics:

- Inspection API:
  - `GentleEngine::inspect_dna_ladders(name_filter?)`
  - Returns structured ladder metadata:
    - `schema` (`gentle.dna_ladders.v1`)
    - `ladder_count`
    - `ladders[]` (`name`, `loading_hint`, `min_bp`, `max_bp`, `band_count`, `bands`)
- Export operation:
  - `ExportDnaLadders { path, name_filter? }`
  - Writes the same structured payload to JSON at `path`.
  - Optional `name_filter` applies case-insensitive name matching before export.

RNA ladder catalog semantics:

- Inspection API:
  - `GentleEngine::inspect_rna_ladders(name_filter?)`
  - Returns structured ladder metadata:
    - `schema` (`gentle.rna_ladders.v1`)
    - `ladder_count`
    - `ladders[]` (`name`, `loading_hint`, `min_nt`, `max_nt`, `band_count`, `bands`)
- Export operation:
  - `ExportRnaLadders { path, name_filter? }`
  - Writes the same structured payload to JSON at `path`.
  - Optional `name_filter` applies case-insensitive name matching before export.

Historical screenshot artifact contract (currently disabled):

- Guardrail:
  - command is currently rejected by security policy even when
    `--allow-screenshots` is provided.
- Command surface:
  - direct CLI: `gentle_cli screenshot-window OUTPUT.png`
  - shared shell (CLI and GUI shell panel): `screenshot-window OUTPUT.png`
- Scope and safety:
  - captures only the active/topmost GENtle window
  - window lookup is native AppKit in-process (no AppleScript automation path)
  - command is primarily intended for GUI shell contexts with an active window
  - rejects full-desktop capture and non-GENtle targets
  - rejects request if no eligible active GENtle window is available
  - current backend support is macOS (`screencapture`); non-macOS returns
    unsupported
- Output:
  - writes an image file at caller-provided `OUTPUT` path (custom filename
    supported)
  - recommended default image format is inferred from extension (e.g. `.png`)
- Result payload shape:

```json
{
  "schema": "gentle.screenshot.v1",
  "path": "docs/images/gui-main.png",
  "window_title": "GENtle - pGEX-3X",
  "captured_at_unix_ms": 1768860000000,
  "pixel_width": 1680,
  "pixel_height": 1020,
  "backend": "macos.screencapture"
}
```

Operation progress/cancellation semantics:

- `apply_with_progress` and workflow progress callbacks receive
  `OperationProgress` updates.
- Callback return value:
  - `true`: continue
  - `false`: request cancellation
- Current event families:
  - `Tfbs`
  - `GenomePrepare`
  - `GenomeTrackImport`
  - `RnaReadInterpret`
- Current cancellation support:
  - genome preparation supports cooperative cancellation plus optional
    `timeout_seconds` timeboxing and reports deterministic cancellation/timeout
    outcomes.
  - genome-track imports support cooperative cancellation and return partial
    import warnings.
  - RNA-read interpretation uses cooperative callback checks while emitting
    periodic progress snapshots (including seed-confirmation histogram bins).

`Pcr` semantics (current):

- Exact primer matching on linear templates.
- Enumerates all valid amplicons formed by forward-primer matches and downstream
  reverse-primer binding matches.
- `unique = true` requires exactly one amplicon; otherwise fails.
- `output_id` may only be used when exactly one amplicon is produced.

`PcrAdvanced` semantics:

- Primer spec fields:
  - `sequence` (full primer, 5'->3')
  - `anneal_len` (3' suffix length used for template binding)
  - `max_mismatches` (allowed mismatches within anneal part)
  - `require_3prime_exact_bases` (hard exact-match requirement at primer 3' end)
  - `library_mode` (`Enumerate` or `Sample`) for degenerate/IUPAC primers
  - `max_variants` cap for primer-library expansion
  - `sample_seed` deterministic seed when `library_mode = Sample`
- Supports 5' tails and mismatch-mediated mutagenesis.
- Supports degenerate/randomized synthetic primers via IUPAC codes.
- Product is constructed from:
  - full forward primer sequence
  - template interior between forward and reverse anneal windows
  - reverse-complement of full reverse primer sequence

`PcrMutagenesis` semantics:

- Builds on `PcrAdvanced` primer behavior.
- Accepts explicit SNP intents:
  - `zero_based_position`
  - `reference`
  - `alternate`
- Validates reference bases against the template.
- Filters amplicons to those that introduce requested SNPs.
- `require_all_mutations` (default `true`) controls whether all or at least one
  mutation must be introduced.

`DesignPrimerPairs` contract (implemented baseline):

- Purpose:
  - propose ranked forward/reverse primer pairs for one linear template under
    explicit constraints
  - provide deterministic, machine-readable reports that can be consumed by
    GUI/CLI/scripting/agents
- Operation payload:

```json
{
  "DesignPrimerPairs": {
    "template": "seq_id",
    "roi_start_0based": 1000,
    "roi_end_0based": 1600,
    "forward": {
      "min_length": 20,
      "max_length": 30,
      "location_0based": null,
      "start_0based": null,
      "end_0based": null,
      "min_tm_c": 55.0,
      "max_tm_c": 68.0,
      "min_gc_fraction": 0.35,
      "max_gc_fraction": 0.70,
      "max_anneal_hits": 1,
      "non_annealing_5prime_tail": null,
      "fixed_5prime": null,
      "fixed_3prime": null,
      "required_motifs": [],
      "forbidden_motifs": [],
      "locked_positions": []
    },
    "reverse": {
      "min_length": 20,
      "max_length": 30,
      "location_0based": null,
      "start_0based": null,
      "end_0based": null,
      "min_tm_c": 55.0,
      "max_tm_c": 68.0,
      "min_gc_fraction": 0.35,
      "max_gc_fraction": 0.70,
      "max_anneal_hits": 1,
      "non_annealing_5prime_tail": null,
      "fixed_5prime": null,
      "fixed_3prime": null,
      "required_motifs": [],
      "forbidden_motifs": [],
      "locked_positions": []
    },
    "pair_constraints": {
      "require_roi_flanking": false,
      "required_amplicon_motifs": [],
      "forbidden_amplicon_motifs": [],
      "fixed_amplicon_start_0based": null,
      "fixed_amplicon_end_0based_exclusive": null
    },
    "min_amplicon_bp": 120,
    "max_amplicon_bp": 1200,
    "max_tm_delta_c": 2.0,
    "max_pairs": 200,
    "report_id": "tp73_roi_primers_v1"
  }
}
```

- `max_tm_delta_c`, `max_pairs`, `report_id`, and `pair_constraints` are optional in current
  implementation:
  - `max_tm_delta_c` default: `2.0`
  - `max_pairs` default: `200`
  - `report_id` default: auto-generated deterministic-safe id stem
  - `pair_constraints` default:
    `{"require_roi_flanking":false,"required_amplicon_motifs":[],"forbidden_amplicon_motifs":[],"fixed_amplicon_start_0based":null,"fixed_amplicon_end_0based_exclusive":null}`
- Side constraints (`forward`, `reverse`, and qPCR `probe`) accept optional
  sequence-level filters:
  - `non_annealing_5prime_tail` (added to the final oligo but excluded from
    anneal Tm/GC/hit calculations)
  - `fixed_5prime`, `fixed_3prime`
  - `required_motifs[]`, `forbidden_motifs[]`
  - `locked_positions[]` entries (`offset_0based`, single IUPAC `base`)
- Built-in primer-ranking heuristics (internal and Primer3 pair-ranking stage):
  - preferred primer length window: `20..30 bp` (outside window is penalized)
  - 3' GC clamp preference (`G/C` at terminal 3' base)
  - secondary-structure risk penalty (homopolymer and self-complementary runs)
  - primer-dimer risk penalty (global and 3'-anchored complementary runs)

- Report schema:
  - `gentle.primer_design_report.v1`
  - deterministic ordering by score then tie-break fields
  - backend metadata block:
    - `backend.requested` (`auto|internal|primer3`)
    - `backend.used` (`internal|primer3`)
    - optional `backend.fallback_reason`
    - optional `backend.primer3_executable`
    - optional `backend.primer3_version`
  - each pair includes:
    - forward/reverse primer sequence and genomic binding window
    - per-primer diagnostics:
      - `length_bp`
      - `anneal_length_bp`
      - `non_annealing_5prime_tail_bp`
      - `three_prime_base`
      - `three_prime_gc_clamp`
      - `longest_homopolymer_run_bp`
      - `self_complementary_run_bp`
    - estimated `tm_c` and `gc_fraction` for annealing segment only
    - anneal-hit counts per side
    - amplicon start/end/length
    - pair dimer diagnostics:
      - `primer_pair_complementary_run_bp`
      - `primer_pair_3prime_complementary_run_bp`
    - rule-pass flags and aggregate score
  - optional rejection summary buckets (for explainability):
    - out-of-window
    - GC/Tm out of bounds
    - non-unique anneal
    - primer sequence-constraint failure
    - pair constraint failure
    - amplicon-size or ROI-coverage failure

`DesignQpcrAssays` contract (implemented baseline):

- Purpose:
  - propose ranked qPCR assays with three oligos (forward primer, reverse
    primer, internal probe) for one linear template.
- Operation payload shape:
  - same core fields as `DesignPrimerPairs` plus:
    - `probe` (`PrimerDesignSideConstraint`)
    - `max_probe_tm_delta_c` (probe Tm distance to mean primer Tm)
    - `max_assays` (result cap)
  - `pair_constraints` is supported identically to `DesignPrimerPairs` and
    applies to forward/reverse pair proposal before probe selection.
- Current baseline behavior:
  - forward/reverse pair generation follows the same backend routing as
    `DesignPrimerPairs` (`auto|internal|primer3` for pair proposal).
  - probe selection is deterministic, constrained to amplicon interior, and
    reuses the same side sequence-constraint fields (`fixed_5prime`,
    `fixed_3prime`, motifs, locked positions).
  - probe Tm gating is enforced via `max_probe_tm_delta_c`.
- Report schema:
  - `gentle.qpcr_design_report.v1`
  - includes ranked `assays[]` with forward/reverse/probe oligos, amplicon
    window, and rule flags.
  - includes qPCR rejection summary with pair-level and probe-level counters.

Primer-design shell command family (implemented):

- Shared-shell family:
  - `primers design REQUEST_JSON_OR_@FILE [--backend auto|internal|primer3] [--primer3-exec PATH]`
  - `primers design-qpcr REQUEST_JSON_OR_@FILE [--backend auto|internal|primer3] [--primer3-exec PATH]`
  - `primers seed-from-feature SEQ_ID FEATURE_ID`
  - `primers seed-from-splicing SEQ_ID FEATURE_ID`
  - `primers list-reports`
  - `primers show-report REPORT_ID`
  - `primers export-report REPORT_ID OUTPUT.json`
  - `primers list-qpcr-reports`
  - `primers show-qpcr-report REPORT_ID`
  - `primers export-qpcr-report REPORT_ID OUTPUT.json`
- `primers design` expects an operation payload whose root variant is
  `{"DesignPrimerPairs": {...}}`.
- `primers design-qpcr` expects an operation payload whose root variant is
  `{"DesignQpcrAssays": {...}}`.
- `primers seed-from-feature` and `primers seed-from-splicing` are
  non-mutating helper commands that resolve an ROI and emit seeded operation
  payloads for both pair-PCR and qPCR design.
- Response schemas:
  - `gentle.primer_seed_request.v1`
  - `gentle.primer_design_report.v1`
  - `gentle.primer_design_report_list.v1`
  - `gentle.qpcr_design_report.v1`
  - `gentle.qpcr_design_report_list.v1`
- `gentle.primer_seed_request.v1` payload fields:
  - `template`
  - `source` (`kind=feature|splicing`, `feature_id`, and splicing metadata when available)
  - `roi_start_0based`
  - `roi_end_0based_exclusive`
  - `operations.design_primer_pairs` (`{"DesignPrimerPairs": ...}`)
  - `operations.design_qpcr_assays` (`{"DesignQpcrAssays": ...}`)

Dotplot + flexibility operation contract (implemented baseline):

- Dotplot operation:
  - `ComputeDotplot { seq_id, reference_seq_id?, span_start_0based?, span_end_0based?, reference_span_start_0based?, reference_span_end_0based?, mode, word_size, step_bp, max_mismatches?, tile_bp?, store_as? }`
  - `mode`: `self_forward | self_reverse_complement | pair_forward | pair_reverse_complement`
  - pair modes require `reference_seq_id` and use the optional
    `reference_span_start_0based` / `reference_span_end_0based` for the
    y/reference axis.
  - stores payload schema `gentle.dotplot_view.v2`
  - payload includes:
    - sparse match points (`points[]`)
    - per-query-bin reference-distribution boxplot summary
      (`boxplot_bin_count`, `boxplot_bins[]` with
      `min/q1/median/q3/max + hit_count`)
  - guardrails:
    - `word_size >= 1`
    - `step_bp >= 1`
    - query/reference spans must satisfy `0 <= start < end <= sequence_len`
    - pair-evaluation safety limit is enforced for latency protection
    - point count is capped with deterministic truncation warning
- Flexibility operation:
  - `ComputeFlexibilityTrack { seq_id, span_start_0based?, span_end_0based?, model, bin_bp, smoothing_bp?, store_as? }`
  - `model`: `at_richness | at_skew`
  - stores payload schema `gentle.flexibility_track.v1`
  - guardrails:
    - `bin_bp >= 1`
    - same span validation contract as dotplot
    - optional smoothing uses deterministic moving-average bins
- Metadata persistence:
  - metadata key: `dotplot_analysis`
  - store schema: `gentle.dotplot_analysis_store.v1`
  - both dotplots and flexibility tracks are persisted under this key
- Shared-shell command family:
  - `dotplot compute SEQ_ID [--reference-seq REF_SEQ_ID] [--start N] [--end N] [--ref-start N] [--ref-end N] [--mode self_forward|self_reverse_complement|pair_forward|pair_reverse_complement] [--word-size N] [--step N] [--max-mismatches N] [--tile-bp N] [--id DOTPLOT_ID]`
  - `dotplot list [SEQ_ID]`
  - `dotplot show DOTPLOT_ID`
  - `dotplot render-svg SEQ_ID DOTPLOT_ID OUTPUT.svg [--flex-track ID] [--display-threshold N] [--intensity-gain N]`
  - `render-dotplot-svg SEQ_ID DOTPLOT_ID OUTPUT.svg [--flex-track ID] [--display-threshold N] [--intensity-gain N]` (alias)
  - `flex compute SEQ_ID [--start N] [--end N] [--model at_richness|at_skew] [--bin-bp N] [--smoothing-bp N] [--id TRACK_ID]`
  - `flex list [SEQ_ID]`
  - `flex show TRACK_ID`

Splicing-reference derivation + pairwise alignment operation contract (implemented baseline):

- Splicing-reference derivation operation:
  - `DeriveSplicingReferences { seq_id, span_start_0based, span_end_0based, seed_feature_id?, scope?, output_prefix? }`
  - emits multiple derived sequence outputs from one genomic window:
    - DNA window (`..._dna`)
    - one mRNA sequence per transcript lane (`..._mrna_*`, transcript orientation, `T->U`)
    - exon-consecutive artificial reference sequence (`..._exon_reference`)
  - if `seed_feature_id` is omitted, engine selects one overlapping mRNA feature deterministically from the requested span
  - default `scope`: `target_group_target_strand`
- Pairwise alignment operation:
  - `AlignSequences { query_seq_id, target_seq_id, query_span_start_0based?, query_span_end_0based?, target_span_start_0based?, target_span_end_0based?, mode?, match_score?, mismatch_score?, gap_open?, gap_extend? }`
  - `mode`: `global | local` (default `global`)
  - scoring defaults: `match=2`, `mismatch=-3`, `gap_open=-5`, `gap_extend=-1`
  - returns structured payload `sequence_alignment` with spans, score, coverage, identity, and CIGAR-like compact operations string
  - non-mutating operation (no sequence/container state mutation)
- Shared-shell command family:
  - `splicing-refs derive SEQ_ID START_0BASED END_0BASED [--seed-feature-id N] [--scope all_overlapping_both_strands|target_group_any_strand|all_overlapping_target_strand|target_group_target_strand] [--output-prefix PREFIX]`
  - `align compute QUERY_SEQ_ID TARGET_SEQ_ID [--query-start N] [--query-end N] [--target-start N] [--target-end N] [--mode global|local] [--match N] [--mismatch N] [--gap-open N] [--gap-extend N]`

RNA-read interpretation contract (Nanopore cDNA phase-1 baseline):

- Operations:
  - `InterpretRnaReads { seq_id, seed_feature_id, profile, input_path, input_format, scope, origin_mode?, target_gene_ids?, roi_seed_capture_enabled?, seed_filter, align_config, report_id?, report_mode?, checkpoint_path?, checkpoint_every_reads?, resume_from_checkpoint? }`
  - `AlignRnaReadReport { report_id, selection, align_config_override?, selected_record_indices? }`
  - implemented profile: `nanopore_cdna_v1`
  - implemented input format: `fasta` (`.fa/.fasta`, optional `.fa.gz/.fasta.gz`; `.sra` must be converted externally in phase-1)
  - default seed/filter constants:
    - `kmer_len=10`
    - `min_seed_hit_fraction=0.30` (bootstrap default; future SNR calibration track can override policy)
    - `min_weighted_seed_hit_fraction=0.05`
    - `min_unique_matched_kmers=12`
    - `min_chain_consistency_fraction=0.40`
    - `max_median_transcript_gap=4.0`
    - `min_confirmed_exon_transitions=1`
    - `min_transition_support_fraction=0.05`
    - weighted-hit definition:
      - `weighted_hit_fraction = sum(1 / occurrence_count(seed_bits)) / tested_kmers`
      - `occurrence_count` is measured inside the active scoped seed index
    - seed pass gate:
      - `raw_hit_fraction >= min_seed_hit_fraction`
      - `AND weighted_hit_fraction >= min_weighted_seed_hit_fraction`
      - `AND unique_matched_kmers >= min(min_unique_matched_kmers, tested_kmers)`
      - `AND chain_consistency_fraction >= min_chain_consistency_fraction`
      - `AND median_transcript_gap <= max_median_transcript_gap`
      - `AND confirmed_transitions >= min_confirmed_exon_transitions`
      - `AND confirmed_transition_fraction >= min_transition_support_fraction`
  - phase-1 seed-span behavior:
    - full-read hashing is always used for every read
    - `short_full_hash_max_bp`, `long_window_bp`, and `long_window_count`
      remain compatibility fields and currently have no runtime effect
  - sparse-origin behavior:
    - `origin_mode` accepts `single_gene|multi_gene_sparse` (default
      `single_gene`)
    - `target_gene_ids[]` and `roi_seed_capture_enabled` are persisted in the
      report payload for deterministic follow-up runs
    - `multi_gene_sparse` expands local transcript-template indexing with
      transcripts matched from `target_gene_ids[]`
    - `roi_seed_capture_enabled=true` is currently a deterministic no-op with
      explicit warning in report `warnings[]` until the ROI capture layer is
      implemented
  - report compaction and resume behavior:
    - `report_mode=full` keeps retained top hits exactly as ranked
    - `report_mode=seed_passed_only` keeps only retained hits that passed the
      composite seed gate (counters remain based on the full stream)
    - `checkpoint_path` + `checkpoint_every_reads` writes deterministic JSON
      snapshots (`gentle.rna_read_interpret_checkpoint.v1`) during streaming
    - `resume_from_checkpoint=true` resumes from the checkpoint snapshot and
      fast-forwards already-processed records deterministically
  - phase-2 alignment behavior:
    - `AlignRnaReadReport` loads a persisted report and reprocesses a selected
      retained subset (`all|seed_passed|aligned`)
    - optional `selected_record_indices[]` (0-based stored `record_index`)
      overrides the selection preset and aligns only the explicit subset
    - aligner configuration uses `align_config_override` when supplied,
      otherwise the report-stored `align_config`
    - mapping backend uses `bio::alignment::pairwise::banded` with
      `align_band_bp` as band width (`w`) and transcript-seed `kmer_len` as
      seed length (`k`), plus deterministic dense fallback when the banded
      solver yields no mapping
    - updated report fields include:
      - per-hit mapping fields (`best_mapping`, `secondary_mappings`)
      - per-hit `msa_eligible` and `msa_eligibility_reason`
      - aggregate `read_count_aligned` and `retained_count_msa_eligible`
      - refreshed transition/isoform support rows and exon/junction abundance
        frequencies
      - deterministic retained-hit re-ranking by alignment-aware retention rank
- Report persistence:
  - report schema: `gentle.rna_read_report.v1`
  - metadata store schema: `gentle.rna_read_reports.v1`
  - metadata key: `rna_read_reports`
  - `rna-reads list-reports` summary rows include sparse-origin request
    provenance:
    - `origin_mode`
    - `target_gene_count`
    - `roi_seed_capture_enabled`
  - report payload now includes per-report:
    - `exon_support_frequencies[]`
    - `junction_support_frequencies[]`
    - storage/streaming controls:
      - `report_mode` (`full` or `seed_passed_only`)
      - `checkpoint_path` / `checkpoint_every_reads`
      - `resumed_from_checkpoint`
    - request provenance fields:
      - `origin_mode`
      - `target_gene_ids[]`
      - `roi_seed_capture_enabled`
    - `origin_class_counts` (running/final deterministic class tallies)
  - per-hit payload now includes:
    - `origin_class`
    - `origin_reason`
    - `origin_confidence`
    - `strand_confidence`
    - `origin_candidates[]` (selected/plus/minus/seed-chain candidate hints)
    - `best_mapping.alignment_mode` (`semiglobal` preferred, with deterministic
      local fallback when quality is better)
  - alignment inspection payload schema:
    - `gentle.rna_read_alignment_inspection.v1`
    - produced by non-mutating shared-shell inspection command
      `rna-reads inspect-alignments`
- Sample-sheet export:
  - operation: `ExportRnaReadSampleSheet { path, seq_id?, report_ids?, append? }`
  - export schema: `gentle.rna_read_sample_sheet_export.v1`
  - output: TSV with run/read metrics, sparse-origin request provenance
    (`report_mode`, `origin_mode`, `target_gene_count`,
    `target_gene_ids_json`, `roi_seed_capture_enabled`), JSON-serialized
    exon/junction frequency columns, and `origin_class_counts_json` for
    cohort-level downstream analysis.
- Shared-shell command family:
  - `rna-reads interpret SEQ_ID FEATURE_ID INPUT.fa[.gz] [--report-id ID] [--report-mode full|seed_passed_only] [--checkpoint-path PATH] [--checkpoint-every-reads N] [--resume-from-checkpoint|--no-resume-from-checkpoint] [--profile nanopore_cdna_v1] [--format fasta] [--scope all_overlapping_both_strands|target_group_any_strand|all_overlapping_target_strand|target_group_target_strand] [--origin-mode single_gene|multi_gene_sparse] [--target-gene GENE_ID]... [--roi-seed-capture|--no-roi-seed-capture] [--kmer-len N] [--short-max-bp N] [--long-window-bp N] [--long-window-count N] [--min-seed-hit-fraction F] [--min-weighted-seed-hit-fraction F] [--min-unique-matched-kmers N] [--min-chain-consistency-fraction F] [--max-median-transcript-gap F] [--min-confirmed-transitions N] [--min-transition-support-fraction F] [--cdna-poly-t-flip|--no-cdna-poly-t-flip] [--poly-t-prefix-min-bp N] [--align-band-bp N] [--align-min-identity F] [--max-secondary-mappings N]`
  - `rna-reads align-report REPORT_ID [--selection all|seed_passed|aligned] [--record-indices i,j,k] [--align-band-bp N] [--align-min-identity F] [--max-secondary-mappings N]`
  - `rna-reads list-reports [SEQ_ID]`
  - `rna-reads show-report REPORT_ID`
  - `rna-reads inspect-alignments REPORT_ID [--selection all|seed_passed|aligned] [--limit N]`
  - `rna-reads export-report REPORT_ID OUTPUT.json`
  - `rna-reads export-hits-fasta REPORT_ID OUTPUT.fa [--selection all|seed_passed|aligned]`
  - `rna-reads export-sample-sheet OUTPUT.tsv [--seq-id ID] [--report-id ID]... [--append]`
  - `rna-reads export-paths-tsv REPORT_ID OUTPUT.tsv [--selection all|seed_passed|aligned]`
  - `rna-reads export-abundance-tsv REPORT_ID OUTPUT.tsv [--selection all|seed_passed|aligned]`
  - `rna-reads export-score-density-svg REPORT_ID OUTPUT.svg [--scale linear|log]`
  - `rna-reads export-alignments-tsv REPORT_ID OUTPUT.tsv [--selection all|seed_passed|aligned] [--limit N]`
  - `rna-reads export-alignment-dotplot-svg REPORT_ID OUTPUT.svg [--selection all|seed_passed|aligned] [--max-points N]`
  - shell output convenience fields:
    - `rna-reads list-reports` includes `summary_rows[]` with concise
      human-readable provenance lines (`mode`, `origin`, target count,
      ROI-capture flag, read counters)
    - `rna-reads show-report` includes `summary` with the same provenance
      framing for one report
    - `rna-reads inspect-alignments` returns top aligned rows ranked by
      alignment-aware retention score (mapping + seed metrics)
- Alignment-TSV export:
  - operation:
    `ExportRnaReadAlignmentsTsv { report_id, path, selection, limit? }`
  - export schema: `gentle.rna_read_alignment_tsv_export.v1`
  - output: ranked alignment rows as TSV (`rank`, mapping metrics, and seed
    metrics) with optional top-`N` truncation via `limit`.
- Alignment-dotplot export:
  - operation:
    `ExportRnaReadAlignmentDotplotSvg { report_id, path, selection, max_points }`
  - export schema: `gentle.rna_read_alignment_dotplot_svg_export.v1`
  - output: SVG scatter of query coverage vs identity for aligned hits with
    score-colored points and report-config threshold guide.
- `rna-reads export-hits-fasta` header extensions:
  - `exon_path_tx=<transcript_id|none>`
  - `exon_path=<ordinal_path|none>` using `:` for hash-confirmed adjacent
    exon transitions and `-` for unconfirmed adjacency
  - `exon_transitions=<confirmed>/<total>`
  - `rc_applied=<true|false>` (automatic cDNA poly-T reverse-complement
    normalization marker)
  - `origin_class=<...>` plus `origin_conf=<...>` and `strand_conf=<...>`
- cDNA/direct-RNA normalization controls in `seed_filter`:
  - `cdna_poly_t_flip_enabled` (default `true`)
  - `poly_t_prefix_min_bp` (default `18`): minimum T support used by the
    tolerant 5' poly-T-head detector (minor interruptions in the head are
    accepted)
- Scope/strand semantics for `InterpretRnaReads`:
  - `all_overlapping_both_strands`: all overlapping transcripts on both strands
  - `target_group_any_strand`: target-group transcripts only, both strands
  - `all_overlapping_target_strand`: all overlapping transcripts on target
    strand only
  - `target_group_target_strand`: target-group transcripts on target strand only
  - scoring note:
    both-strand modes score against the union of admitted strand-specific
    templates; target-strand modes exclude opposite-strand templates.
  - seed-index note:
    indexed seeds include annotated exon-body and exon-exon transition k-mers
    for admitted transcripts.

Async BLAST shell contract (agent/MCP-ready baseline):

- Shared-shell families (both `genomes` and `helpers` scopes):
  - `blast-start GENOME_ID QUERY_SEQUENCE ...`
  - `blast-status JOB_ID [--with-report]`
  - `blast-cancel JOB_ID`
  - `blast-list`
- Deterministic job payload schemas:
  - `gentle.blast_async_start.v1`
  - `gentle.blast_async_status.v1`
  - `gentle.blast_async_cancel.v1`
  - `gentle.blast_async_list.v1`
- External-binary preflight payload:
  - `blast-start` responses now include `binary_preflight` with schema
    `gentle.blast_external_binary_preflight.v1`.
  - payload includes deterministic `blastn` and `makeblastdb` probe rows with:
    `found`, `version`, `executable`, and resolved `path` diagnostics.
  - equivalent preflight payload is also emitted by synchronous shared-shell
    routes `prepare`, `blast`, and `blast-track`.
- Job status contract:
  - `job_id` stable per process
  - non-terminal states: `queued | running`
  - terminal states: `completed | failed | cancelled`
  - scheduler metadata:
    - `max_concurrent_jobs`
    - `running_jobs`
    - `queued_jobs`
    - `queue_position` (present while state is `queued`)
  - optional final `report` on `blast-status --with-report`
- Durability/restart semantics:
  - BLAST async status snapshots are persisted in project metadata as
    `blast_async_jobs` (`gentle.blast_async_job_store.v1`).
  - On restart/reload, recovered jobs that were previously non-terminal but no
    longer have an active worker context are normalized deterministically:
    - `cancel_requested=true` -> `cancelled`
    - otherwise -> `failed` with explicit restart/reload interruption reason.
  - `blast-start`, `blast-status`, `blast-cancel`, and `blast-list` may mark
    shell state as changed when they persist updated async job snapshots.
- Scheduler policy:
  - async BLAST jobs are executed by a bounded FIFO scheduler (queue + worker slots)
  - default concurrency uses host CPU parallelism
  - optional override via environment variable
    `GENTLE_BLAST_ASYNC_MAX_CONCURRENT` (clamped to `1..256`)
- `gentle_mcp` exposes equivalent tool routes:
  - `blast_async_start`
  - `blast_async_status`
  - `blast_async_cancel`
  - `blast_async_list`

### Workflow

```json
{
  "run_id": "string",
  "ops": ["Operation", "Operation", "..."]
}
```

Notes:

- Splicing Expert `Nanopore cDNA interpretation` uses this same workflow shape
  when you click `Copy Workflow JSON`.
- `Prepare Workflow Op` in the same panel writes `run_id`/`ops` into the GUI
  workflow runner so the exact `InterpretRnaReads` payload can be rerun through
  the generic workflow path.

### OpResult

```json
{
  "op_id": "op-1",
  "created_seq_ids": ["..."],
  "changed_seq_ids": ["..."],
  "warnings": ["..."],
  "messages": ["..."],
  "genome_annotation_projection": null,
  "sequence_alignment": null
}
```

### Error

```json
{
  "code": "InvalidInput|NotFound|Unsupported|Io|Internal",
  "message": "human-readable explanation"
}
```

## State model

`gentle_cli` persists engine state in JSON (`.gentle_state.json` by default).

This supports:

- resumable multi-step workflows
- external inspection
- reproducibility and audit trails

## Recommended AI-agent flow

1. Query `capabilities`
2. Import or initialize state
3. Apply one operation at a time, checking warnings/errors
4. Save/export artifacts
5. Optionally export final state for handoff

## Planned next additions

- richer sequence-editing and annotation operation set
- ligation protocol presets with sticky/blunt compatibility derivation
- render/view model endpoint for frontend-independent graphical representation
- schema publication for strict client-side validation
- CRISPR guide-design next phase:
  - off-target search/ranking contracts
  - on-target efficacy model integration hooks
  - guide-design macro/template expansion into deterministic `Workflow` JSON
  - see draft: `docs/rna_guides_spec.md`
