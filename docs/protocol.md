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

## TFBS region summary contract

GENtle now exposes a portable grouped TFBS summary contract for comparing one
focus window against a wider context window on the same sequence.

Current shared-shell route:

```bash
gentle_cli shell 'features tfbs-summary SEQ_ID --focus START..END [--context START..END] [--min-focus-count N] [--min-context-count N] [--limit N]'
```

First-class operation route:

```json
{"SummarizeTfbsRegion":{"seq_id":"SEQ_ID","focus_start_0based":2900,"focus_end_0based_exclusive":3100,"context_start_0based":0,"context_end_0based_exclusive":6001,"min_focus_occurrences":1,"min_context_occurrences":0,"limit":25}}
```

Portable schema:

- `gentle.tfbs_region_summary.v1`

Request fields:

- `seq_id`
- `focus_start_0based`
- `focus_end_0based_exclusive`
- optional `context_start_0based`
- optional `context_end_0based_exclusive`
- `min_focus_occurrences`
- `min_context_occurrences`
- optional `limit`

Result fields:

- sequence/focus/context bounds and widths
- total TFBS hit counts in the focus and context spans
- grouped rows keyed by TF name with:
  - `motif_ids`
  - `focus_occurrences`
  - `context_occurrences`
  - `outside_focus_occurrences`
  - focus/context/outside densities per kb
  - focus-vs-context and focus-vs-outside density ratios

Grouping policy:

- prefer `bound_moiety`
- otherwise `standard_name`
- otherwise `gene`
- otherwise `name`
- otherwise `tf_id`

## Draft design resources

### `gentle.gibson_assembly_plan.v1`

Purpose:

- describe one Gibson cloning project in a destination-first way,
- separate user-specified plan inputs from derived design consequences,
- provide one canonical JSON artifact that future routines, primer design, and
  protocol-cartoon rendering can all read from.

Status:

- destination-first single-insert and ordered multi-insert plans are now
  consumed by the shared Gibson preview/apply path (`gibson preview ...`,
  `gibson apply`, and the `Patterns -> Gibson...` specialist window)
- current limit:
  - multi-insert execution currently assumes a defined destination opening
  - `existing_termini` remains the single-fragment handoff path for now

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
  - optional design request:
    - `desired_unique_restriction_site_enzyme_name`
      asks the shared preview/apply path to try introducing one new unique
      REBASE cutter site on one terminal overlap if the assembled product can
      still remain uniquely cut there
- `derived_design`
  - derived overlap sequences
  - primer design suggestions
  - advisory notes and validation outcomes

Current draft value vocabulary:

- `destination.topology_before_opening`
  - `linear`
  - `circular`
- `destination.opening.mode`
  - `existing_termini`
    - use the current termini of an already-linear destination sequence
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

### `gentle.gibson_assembly_preview.v1`

Purpose:

- provide one deterministic, non-mutating preview response for the current
  Gibson specialist flow,
- keep GUI, shared shell, and direct CLI on the same overlap/primer/cartoon
  derivation path.

Current shared entry points:

- `gibson preview PLAN_JSON_OR_@FILE [--output OUTPUT.json]`
- `gibson apply PLAN_JSON_OR_@FILE`
- GUI specialist window: `Patterns -> Gibson...`

Top-level structure:

- `schema`, `plan_id`, `title`, `summary`
- `can_execute`
- `destination`
  - resolved opening mode/span or cutpoint and actual topology
- `fragments[]`
  - resolved ordered insert rows (fragment id, template seq id, orientation,
    length)
- `insert`
  - compatibility mirror of the first insert row for older single-insert
    consumers
- `resolved_junctions[]`
  - overlap bp, left/right member contributions, overlap Tm, resolved overlap
    sequence, source note
- `primer_suggestions[]`
  - full primer sequence plus explicit `overlap_5prime` and
    `priming_3prime` segments
- `warnings[]`, `errors[]`, `notes[]`
  - includes the shared Tₘ-model note used by GUI/CLI so the assumptions stay
    visible to the user
  - notes also carry explicit design-review guidance that separates:
    - overlap-side success/failure
    - PCR 3' priming-side success/failure
    so adapters can explain when the current blocker is priming rather than
    Gibson overlap derivation
- `suggested_design_adjustments[]`
  - optional structured next-step relaxations when overlap derivation already
    succeeds and the remaining blocker is only the 3' priming window
  - current v1 targets:
    - increasing `priming_segment_max_length_bp`
    - lowering `priming_segment_tm_min_celsius`
  - intended for adapters to offer deterministic “apply and rerun preview”
    actions without parsing prose notes
- `unique_restriction_site`
  - optional structured outcome for a requested
    `validation_policy.desired_unique_restriction_site_enzyme_name`
  - reports whether the requested site was:
    - already unique in the assembled product
    - newly engineered on one terminal overlap
  - carries the enzyme name, terminal side/junction, engineered overlap
    sequence, motif offset, mutation count, and user-facing message so adapters
    do not have to infer this from notes/error prose
- `cartoon`
  - built-in protocol id plus template bindings for single-insert flows
  - multi-insert previews may instead carry one fully resolved
    `ProtocolCartoonSpec` directly
  - intended to stay mechanism-first:
    - show resolved fragment flow and achieved homology/overlap relationships
    - preserve strand-specific 5' chew-back / exposed-tail geometry rather
      than flattening the mechanism to duplex-only blocks
    - avoid drawing full primer objects or low-level PCR parameterization inside
      the cartoon itself
    - keep primer sequences, priming segments, Tm assumptions, and related PCR
      details in adjacent textual/review payloads instead
- `routine_handoff`
  - best-effort Routine Assistant handoff metadata for existing execution paths

Current v1 scope and limits:

- one or more insert fragments in an explicit ordered chain
- destination-first order:
  `destination_left -> insert_1 -> ... -> insert_n -> destination_right`
- the shared preview derives `n + 1` explicit Gibson junctions for `n` inserts
- terminal overlaps are derived from destination context; internal junctions
  are normalized from the adjacent fragment ends / partition rules
- user influence over PCR design stays high-level and Gibson-specific:
  overlap bp range, minimum overlap Tm, priming-segment Tm window, and
  priming-segment length window
- current execution limitation:
  - multi-insert apply currently requires `destination.opening.mode=defined_site`
  - `existing_termini` remains the single-fragment path used by the current
    Routine Assistant handoff
- current unique-site engineering limitation:
  - only the single-insert `defined_site` path is supported
  - only palindromic cutter recognition sequences are currently handled
  - overlap windows must be non-wrapping in the displayed destination sequence
- current Tₘ fields use the shared GENtle nearest-neighbor estimate with fixed
  assumptions:
  - exact complement
  - 50 mM monovalent salt
  - 250 nM total oligo concentration
  - no mismatch/dangling-end/Mg correction
  - fallback to the simple 2/4 estimate for ambiguous or very short sequences
- generic PCR/qPCR request editing is intentionally out of scope for this
  specialist flow
- mutating execution now exists as engine operation
  `ApplyGibsonAssemblyPlan`:
  - consumes the same plan JSON
  - creates deterministic sequence outputs for:
    - left insert primer
    - right insert primer
    - assembled product
  - creates one shared serial arrangement for downstream gel review:
    - original destination vector
    - ordered insert lane(s)
    - assembled product
    - recommended DNA ladders carried with the arrangement for flanking export
  - transfers destination and insert features onto the assembled product
    deterministically through the shared engine path
  - destination features intersecting the consumed opening are now projected
    when a truthful rewrite is available:
    - one-sided overlaps are trimmed to the surviving product span
    - simple spanning features can survive as multipart remnants
    - MCS-like annotations are projected to the edited locus and revalidated
      against actual restriction-enzyme sites on the assembled product
  - the MCS cross-check is product-aware:
    - `mcs_expected_sites` is rewritten to the currently unique cutter set for
      that annotated region on the assembled product
    - `mcs_expected_sites_original`, `mcs_region_sites`,
      `mcs_nonunique_sites`, `mcs_gained_unique_sites`, and
      `mcs_lost_or_nonunique_sites` preserve the cross-check result
    - insert-derived sequence may introduce new sites, and those new sites are
      considered during the same validation pass
  - records one operation-log row so GUI lineage/CLI state replay can reopen
    the specialist from the saved plan without silently re-running it

Normalization/derivation phases:

1. Resolve the destination opening into explicit `dest_left` / `dest_right`
   terminal context.
   - for cutter-derived openings, the resolved coordinates represent the actual
     cleavage window between the recessed termini rather than the whole
     recognition span
   - equal start/end is therefore valid and means a blunt cutpoint
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

### `gentle.prepared_cache_inspection.v1`

Purpose:

- provide one deterministic, non-mutating inspection payload for prepared
  reference/helper cache roots,
- let GUI, shared shell, and direct CLI report exactly what GENtle created
  locally before any deletion happens.

Current shared entry points:

- `cache inspect [--references|--helpers|--both] [--cache-dir PATH ...]`
- GUI specialist window: `Genome -> Clear Caches...`
- `Prepared References... -> Clear Caches...`

Top-level structure:

- `schema`, `cache_roots[]`
- `entries[]`
- `entry_count`, `total_size_bytes`, `total_file_count`

Entry structure:

- `entry_id`
- `classification`
  - `prepared_install`
  - `orphaned_remnant`
- `cache_root`, `path`
- `artifact_stats[]`
  - `group`
    - `cached_sources`
    - `derived_indexes`
    - `blast_db`
  - `total_size_bytes`
  - `file_count`
- `total_size_bytes`, `file_count`

Inspection rules:

- inspection stays rooted in the selected cache roots only
- default roots are adapter-facing conventions:
  - references: `data/genomes`
  - helpers: `data/helper_genomes`
- orphaned remnants are inspectable even when they are not backed by a
  manifest
- catalog JSON, project state files, MCP/runtime files, backdrop/runtime
  caches, and developer build artifacts are out of scope

### `gentle.prepared_cache_cleanup.v1`

Purpose:

- provide one deterministic cleanup result payload for conservative prepared
  cache deletion workflows,
- keep partial rebuild/reindex cleanup and full prepared-install deletion on
  the same shared contract across GUI/CLI/shell.

Current shared entry points:

- `cache clear blast-db-only|derived-indexes-only|selected-prepared|all-prepared-in-cache ...`
- GUI specialist window: `Genome -> Clear Caches...`

Top-level structure:

- `schema`, `mode`, `cache_roots[]`
- `selected_prepared_ids[]`
- `selected_prepared_paths[]`
- `include_orphaned_remnants`
- `results[]`
- `entry_count`, `removed_item_count`, `removed_bytes`, `removed_file_count`

Per-entry result structure:

- `entry_id`
- `classification`
  - `prepared_install`
  - `orphaned_remnant`
- `cache_root`, `path`
- `removed`
- `removed_artifact_groups[]`
- `removed_bytes`, `removed_file_count`
- `skipped_reason?`

Cleanup modes:

- `blast_db_only`
  - remove only BLAST DB sidecars for selected manifest-backed installs
- `derived_indexes_only`
  - remove BLAST DB sidecars plus `sequence.fa.fai` and `genes.json`
  - cached sources and manifests remain so reindex-from-cached-files still
    works
- `selected_prepared_installs`
  - remove only explicitly selected prepared installs
  - optional `include_orphaned_remnants` also removes orphaned remnants under
    the same selected roots
- `all_prepared_in_cache`
  - remove all prepared installs under the selected roots
  - optional `include_orphaned_remnants` extends that deletion to orphaned
    remnants

Cleanup rules:

- `blast_db_only` and `derived_indexes_only` apply only to manifest-backed
  prepared installs
- selective cleanup modes accept either `selected_prepared_ids[]` or
  `selected_prepared_paths[]`
- `selected_prepared_paths[]` are the precise selector when duplicate prepared
  ids exist across multiple selected cache roots
- orphaned remnants can only be deleted through the full-delete modes
- cleanup never scans the whole workspace; it only touches the selected roots
- cleanup does not treat catalog JSON, `.gentle_state.json`, MCP/runtime files,
  backdrop/runtime caches, or `target/` as cache

## Core entities

### ProjectState

```json
{
  "sequences": {"seq_id": "DNAsequence object"},
  "metadata": {"any": "json"},
  "display": {"ui_visibility_tfbs_and_linear_viewport_state": "..."},
  "lineage": {"nodes": {}, "edges": []},
  "parameters": {"max_fragments_per_container": 80000},
  "container_state": {
    "containers": {},
    "arrangements": {},
    "racks": {},
    "seq_to_latest_container": {}
  }
}
```

Semantic interpretation:

- In GUI terms, a project window represents a wet-lab container context.
- A container may map to multiple candidate sequences/fragments.
- Explicit container objects are first-class state (`container_state`) and are
  indexed from sequence ids via `seq_to_latest_container`.
- Arrangements stay the semantic experiment-order layer.
- Racks are the linked physical placement layer and may host one or more
  arrangements without changing arrangement identity.

### Rack placement entities

- `RackProfileKind`
  - built-in physical carriers:
    - `small_tube_4x6`
    - `plate_96`
    - `plate_384`
  - persisted custom snapshots use:
    - `custom`
- `RackProfileSnapshot`
  - persisted row/column/fill-direction/blocked-slot snapshot used by one saved rack
  - `fill_direction`
    - `row_major`
    - `column_major`
  - `blocked_coordinates[]`
    - normalized A1-style coordinate list
- `Rack`
  - one saved physical rack/plate draft
- `RackPlacementEntry`
  - one occupied A1-style coordinate on that rack
  - points back to:
    - `arrangement_id`
    - arrangement-local `order_index`
    - one `occupant`
- `RackOccupant`
  - `container`
  - `ladder_reference`

Rack-placement invariants:

- rack placement consumes arrangement order instead of duplicating experiment
  meaning in a second free-form list
- default placement is deterministic:
  - choose the smallest fitting built-in profile
  - fill row-major
  - use A1-style coordinates
- saved rack snapshots may then refine physical layout with:
  - `fill_direction = row_major|column_major`
  - `blocked_coordinates[]`
- A1-style row labels continue beyond `Z` as `AA`, `AB`, ...
- moving one sample or arrangement block is shift-neighbor by default; it
  preserves occupied order instead of creating arbitrary holes

### `gentle.rack_state.v1`

Purpose:

- provide one deterministic inspection payload for saved rack state
- keep GUI rack view and CLI/shell inspection on one shared state contract

Current shared entry point:

- `racks show RACK_ID`

Top-level structure:

- `schema`
- `rack`
- `placements[]`

Placement payload:

- `coordinate`
- `arrangement_id`
- `order_index`
- `role_label`
- `occupant`
  - `kind = container`
    - `container_id`
    - `container_name?`
    - `seq_id?`
  - `kind = ladder_reference`
    - `ladder_name`
  - `kind = empty`

### Operation

Current draft operations:

- `LoadFile { path, as_id? }`
- `SaveFile { seq_id, path, format }`
- `RenderSequenceSvg { seq_id, mode, path }`
  - linear exports honor the current stored linear viewport in `display`
    (`linear_view_start_bp` / `linear_view_span_bp`) when that viewport is a
    proper subsequence crop
  - single-base `variation` features render as baseline markers in linear SVG
    output rather than as generic detached feature blocks
  - linear exports now also mark transcription starts/directions for
    strand-bearing `gene`/`mRNA`/`CDS`/`promoter` features and suppress
    unlabeled fallback coordinate text that would otherwise clutter
    figure-oriented exports
  - linear exports also prefer gene-style labels over accession-only transcript
    ids when possible and compact nearby repeated non-gene labels
  - direction-bearing `mRNA`/`promoter` bars render with arrowed ends, and the
    linear TSS cue uses a short hooked arrow so direction survives
    figure-oriented contexts
  - circular exports now use a transparent canvas and render single-base
    `variation` features as explicit radial markers on the DNA ring
  - circular exports also mark transcription starts for strand-bearing
    `gene`/`mRNA`/`CDS`/`promoter` features with a short arrow shaft plus
    direction arrowhead
  - circular exports also use a slightly larger ring and larger label fonts so
    figure-oriented construct maps stay readable when embedded in docs
- `RenderDotplotSvg { seq_id, dotplot_id, path, flex_track_id?, display_density_threshold?, display_intensity_gain? }`
- `RenderFeatureExpertSvg { seq_id, target, path }`
  - shared renderer contract across GUI/CLI/JS/Lua for TFBS/restriction/splicing/isoform expert exports
  - splicing SVG includes explicit junction-support counts, frequency-encoded transcript-vs-exon matrix coloring, predicted exon->exon transition matrix support coloring, exon `len%3` (genomic-length modulo 3) cues, and CDS flank phase edge coloring (`0/1/2`) when transcript `cds_ranges_1based` are available
- `RenderIsoformArchitectureSvg { seq_id, panel_id, path }`
- `RenderRnaStructureSvg { seq_id, path }`
- `RenderLineageSvg { path }`
- `RenderPoolGelSvg { inputs, path, ladders?, container_ids?, arrangement_id? }`
- `CreateArrangementSerial { container_ids, arrangement_id?, name?, ladders? }`
- `SetArrangementLadders { arrangement_id, ladders? }`
- `CreateRackFromArrangement { arrangement_id, rack_id?, name?, profile? }`
- `PlaceArrangementOnRack { arrangement_id, rack_id }`
- `MoveRackPlacement { rack_id, from_coordinate, to_coordinate, move_block? }`
- `MoveRackArrangementBlocks { rack_id, arrangement_ids[], to_coordinate }`
- `SetRackProfile { rack_id, profile }`
- `ApplyRackTemplate { rack_id, template }`
- `SetRackFillDirection { rack_id, fill_direction }`
- `SetRackProfileCustom { rack_id, rows, columns }`
- `SetRackBlockedCoordinates { rack_id, blocked_coordinates }`
- `ExportRackLabelsSvg { rack_id, path, arrangement_id?, preset }`
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
- `PcrOverlapExtensionMutagenesis { ... }` (implemented baseline; insertion/deletion/replacement overlap-extension flow)
- `DesignQpcrAssays { ... }` (implemented baseline; forward/reverse/probe)
- `ComputeDotplot { seq_id, reference_seq_id?, span_start_0based?, span_end_0based?, reference_span_start_0based?, reference_span_end_0based?, mode, word_size, step_bp, max_mismatches?, tile_bp?, store_as? }` (implemented baseline, self + pairwise)
- `ComputeFlexibilityTrack { seq_id, span_start_0based?, span_end_0based?, model, bin_bp, smoothing_bp?, store_as? }` (implemented baseline)
- `DeriveSplicingReferences { seq_id, span_start_0based, span_end_0based, seed_feature_id?, scope?, output_prefix? }` (implemented baseline; emits derived DNA window + mRNA isoforms + exon-reference sequence)
- `AlignSequences { query_seq_id, target_seq_id, query_span_start_0based?, query_span_end_0based?, target_span_start_0based?, target_span_end_0based?, mode?, match_score?, mismatch_score?, gap_open?, gap_extend? }` (implemented baseline; returns structured pairwise local/global report in `OpResult.sequence_alignment`)
- `ImportSequencingTrace { path, trace_id?, seq_id? }` (implemented baseline; imports one ABI/AB1 or SCF evidence file into the shared sequencing-trace store without mutating construct sequences)
- `ListSequencingTraces { seq_id? }`
- `ShowSequencingTrace { trace_id }`
- `ConfirmConstructReads { expected_seq_id, baseline_seq_id?, read_seq_ids?, trace_ids?, targets?, alignment_mode?, match_score?, mismatch_score?, gap_open?, gap_extend?, min_identity_fraction?, min_target_coverage_fraction?, allow_reverse_complement?, report_id? }` (implemented baseline; accepts already-loaded read sequences and/or imported sequencing traces as evidence inputs into one shared confirmation report, with optional baseline context for intended-edit vs reversion classification)
- `InterpretRnaReads { seq_id, seed_feature_id, profile, input_path, input_format, scope, origin_mode?, target_gene_ids?, roi_seed_capture_enabled?, seed_filter, align_config, report_id?, report_mode?, checkpoint_path?, checkpoint_every_reads?, resume_from_checkpoint? }` (Nanopore cDNA phase-1 seed-filter pass; `multi_gene_sparse` expands local transcript-template indexing, while ROI capture remains planned)
- `AlignRnaReadReport { report_id, selection, align_config_override?, selected_record_indices? }` (Nanopore cDNA phase-2 retained-hit alignment pass; updates mapping/MSA/abundance report fields and re-ranks retained hits by alignment-aware retention rank)
- `ListRnaReadReports { seq_id? }`
- `ShowRnaReadReport { report_id }`
- `ExportRnaReadReport { report_id, path }`
- `ExportRnaReadHitsFasta { report_id, path, selection, selected_record_indices?, subset_spec? }`
- `ExportRnaReadSampleSheet { path, seq_id?, report_ids?, append? }`
- `ExportRnaReadExonPathsTsv { report_id, path, selection, selected_record_indices?, subset_spec? }`
- `ExportRnaReadExonAbundanceTsv { report_id, path, selection, selected_record_indices?, subset_spec? }`
- `ExportRnaReadScoreDensitySvg { report_id, path, scale, variant }`
- `ExportRnaReadAlignmentsTsv { report_id, path, selection, limit?, selected_record_indices?, subset_spec? }`
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
- `ExtractGenomeGene { genome_id, gene_query, occurrence?, output_id?, extract_mode?, promoter_upstream_bp?, annotation_scope?, max_annotation_features?, include_genomic_annotation?, catalog_path?, cache_dir? }`
  - `annotation_scope` accepts `none|core|full` and defaults to `core` when omitted.
  - `max_annotation_features` is an optional safety cap (0 or omitted = unlimited for explicit requests).
  - legacy `include_genomic_annotation` is still accepted (`true` -> `core`, `false` -> `none`) for compatibility.
  - operation results include `genome_annotation_projection` telemetry (requested/effective scope, feature counts, fallback metadata).
  - for helper genome IDs containing `pUC18`/`pUC19`, the same deterministic MCS fallback annotation behavior applies when an MCS feature is missing; non-unique motif matches are warned and skipped.
- `ExtendGenomeAnchor { seq_id, side, length_bp, output_id?, catalog_path?, cache_dir?, prepared_genome_id? }`
- `VerifyGenomeAnchor { seq_id, catalog_path?, cache_dir?, prepared_genome_id? }`

Sequencing-trace evidence notes:

- raw traces are stored separately from `SequencingConfirmationReport`
  payloads; importing a trace does not run confirmation and does not mutate any
  sequence entry
- `ImportSequencingTrace` currently auto-detects:
  - ABIF/AB1 via `ABIF` magic bytes
  - SCF via `.scf` magic bytes
- stored `SequencingTraceRecord` payloads preserve:
  - file-supplied called bases
  - called-base confidence arrays when available
  - peak locations when available
  - raw per-channel intensity arrays when available
  - compact per-channel trace-length summaries
  - optional clip window metadata when present in the source file
  - optional sample/run/machine metadata when present in the source file
- trace-aware confirmation now reuses the same `SequencingConfirmationReport`
  model:
  - `ConfirmConstructReads` accepts `trace_ids` in addition to `read_seq_ids`
  - `ConfirmConstructReads` accepts optional `baseline_seq_id` so the expected
    construct remains primary truth while baseline context can distinguish
    intended edits from reference reversions
  - per-evidence rows expose evidence kind plus optional `trace_id`
  - target support/contradiction ids may now refer to imported trace ids when
    traces provide the relevant evidence
  - report payloads now include:
    - `baseline_seq_id?`
    - per-target `expected_bases?` / `baseline_bases?` for expected-edit loci
    - `variants[]` rows with observed allele, evidence id, confidence summary,
      peak center, and classification:
      `expected_match|intended_edit_confirmed|reference_reversion|unexpected_difference|low_confidence_or_ambiguous|insufficient_evidence`
  - persisted confirmation reports now project as lineage analysis artifacts in
    both the GUI lineage workspace and shared `RenderLineageSvg` export:
    nodes are keyed by `report_id`, attach to the expected construct plus
    optional baseline/reference sequence, and reopen the sequencing-confirmation
    specialist on that stored report in GUI adapters
- `SuggestSequencingPrimers { expected_seq_id, primer_seq_ids[], confirmation_report_id?, min_3prime_anneal_bp, predicted_read_length_bp }`
  - non-mutating helper for sequencing-confirmation review and primer coverage
    planning
  - `primer_seq_ids[]` may be empty when `confirmation_report_id` is present:
    that mode proposes fresh sequencing primers for unresolved loci using the
    expected construct plus the saved report context
  - returns `SequencingPrimerOverlayReport` with per-hit orientation, anneal
    span, predicted read span, optional coverage annotations against a
    persisted sequencing-confirmation report, per-problem guidance rows naming
    the best existing primer hit for unresolved targets or variant loci, and
    `proposals[]` rows for fresh primer candidates when no good existing hit is
    available
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
- `FetchDbSnpRegion { rs_id, genome_id, flank_bp?, output_id?, annotation_scope?, max_annotation_features?, catalog_path?, cache_dir? }`
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
- shared-shell dbSNP route:
  - `dbsnp fetch RS_ID GENOME_ID [--flank-bp N] [--output-id ID] [--annotation-scope none|core|full] [--max-annotation-features N] [--catalog PATH] [--cache-dir PATH]`
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
  - current stable targets:
    - `prepared-references`
    - `prepare-reference-genome`
    - `retrieve-genome-sequence`
    - `blast-genome-sequence`
    - `import-genome-track`
    - `pcr-design`
    - `sequencing-confirmation`
    - `agent-assistant`
    - `prepare-helper-genome`
    - `retrieve-helper-sequence`
    - `blast-helper-sequence`
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
- Optional `container_ids` renders one lane per referenced stored container.
- Optional `arrangement_id` renders one lane per stored serial-arrangement lane.
- Computes pool molecular-weight proxy from sequence bp lengths.
- Chooses one or two ladders to span pool range:
  - from explicit `ladders` list when provided
  - otherwise from saved arrangement ladders when `arrangement_id` is present
    and the arrangement stores a ladder choice
  - otherwise from built-in ladder catalog (auto mode)
- Renders ladder lanes plus pooled band lane as SVG artifact.

`CreateArrangementSerial` semantics:

- Persists an ordered serial-lane setup over stored containers.
- Optional `ladders` can store one symmetric ladder or one left/right ladder
  pair for later gel preview/export reuse.
- Also materializes one default physical rack draft:
  - choose the smallest built-in rack/plate profile that fits the arrangement
    payload plus ladder-reference positions
  - place the arrangement block row-major in that rack
  - link the arrangement back to the resulting `default_rack_id`

`SetArrangementLadders` semantics:

- Mutates an existing serial arrangement in place.
- `ladders = null` clears back to shared engine auto ladder selection.
- one ladder name means the same ladder is used on both sides during
  arrangement-based gel preview/export.
- two ladder names mean explicit left/right ladder selection.

`CreateRackFromArrangement` semantics:

- Creates one new physical rack/plate draft from one stored arrangement.
- Optional `profile` overrides the default smallest-fitting profile choice.
- If `profile` is omitted, the engine chooses in this order:
  - `small_tube_4x6`
  - `plate_96`
  - `plate_384`
- Placement is row-major and preserves arrangement order.
- Ladder-bearing arrangements reserve left/right ladder-reference positions in
  the same contiguous block.

`PlaceArrangementOnRack` semantics:

- Places one arrangement onto an existing rack as one contiguous block at the
  next free region in fill order.
- Existing rack occupants stay in order; the appended arrangement does not
  reorder earlier blocks.
- Shared racks are therefore possible without losing arrangement identity.

`MoveRackPlacement` semantics:

- Moves one occupied rack coordinate within one saved rack.
- `move_block=false` means move one sample within its arrangement block and
  shift neighboring occupied positions to preserve order.
- `move_block=true` means move the whole arrangement block and shift later
  occupied blocks in fill order.
- The operation is order-preserving by design; it does not treat arbitrary
  holes as the primary editing model.

`MoveRackArrangementBlocks` semantics:

- Moves two or more selected arrangement blocks together within one saved rack.
- `arrangement_ids` are normalized against the rack's current occupied order;
  the shared engine preserves the rack-ordering of selected blocks even if the
  request lists them in another order.
- The selected blocks move as one contiguous combined group.
- Later occupied blocks shift in fill order to keep the rack contiguous.
- This is the shared engine contract behind rack-editor multi-select moves.

`SetRackProfile` semantics:

- Reprojects one saved rack onto another built-in profile.
- Existing arrangement order is preserved while coordinates are reflowed under
  the target profile geometry.
- Existing fill direction is preserved.
- Existing blocked coordinates are preserved when still in-bounds for the new
  geometry; out-of-bounds blocked coordinates are dropped deterministically.

`ApplyRackTemplate` semantics:

- Applies one engine-owned quick-authoring template on top of an existing rack
  snapshot.
- Built-in templates:
  - `bench_rows`
    - `fill_direction = row_major`
    - `blocked_coordinates = []`
  - `plate_columns`
    - `fill_direction = column_major`
    - `blocked_coordinates = []`
  - `plate_edge_avoidance`
    - `fill_direction = column_major`
    - `blocked_coordinates = outer perimeter of the current profile`
- Existing arrangement order is preserved while occupied coordinates are
  reflowed onto the resulting available slots.
- `plate_edge_avoidance` requires at least a `3 x 3` profile so an interior
  region remains after blocking the perimeter.

`SetRackFillDirection` semantics:

- Reprojects one saved rack onto the same geometry with a different fill order.
- Supported values:
  - `row_major`
  - `column_major`
- Existing arrangement order is preserved while occupied coordinates are
  reassigned under the new fill order.

`SetRackProfileCustom` semantics:

- Reprojects one saved rack onto one custom A1-style geometry.
- `rows` and `columns` are persisted directly in the rack profile snapshot.
- Existing fill direction is preserved.
- Existing blocked coordinates are preserved when still in-bounds for the new
  geometry; out-of-bounds blocked coordinates are dropped deterministically.
- Existing arrangement order is preserved while coordinates are reflowed under
  the custom geometry.
- A1-style row labels continue beyond `Z` as `AA`, `AB`, ...

`SetRackBlockedCoordinates` semantics:

- Persists one normalized blocked/reserved coordinate set on the rack profile.
- Blocked coordinates are excluded from placement capacity and fill-order
  reflow.
- Existing arrangement order is preserved while occupied coordinates are
  reassigned onto the remaining available positions.
- Duplicate blocked coordinates are removed deterministically.

`ExportRackLabelsSvg` semantics:

- Writes one deterministic SVG label sheet for a saved rack.
- Optional `arrangement_id` restricts output to labels belonging to one
  arrangement block on that rack.
- `preset` is engine-owned and defaults to `compact_cards`.
- Built-in presets:
  - `compact_cards`
  - `print_a4`
  - `wide_cards`
- Label rows currently include:
  - rack id
  - position
  - role
  - container/ladder display name
  - sequence id when sequence-backed
  - bp length/topology when sequence-backed

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
  - `protocol` (built-in ids now include `gibson.two_fragment`,
    `gibson.single_insert_dual_junction`, `pcr.assay.pair`,
    `pcr.assay.pair.no_product`, `pcr.assay.pair.with_tail`, and
    `pcr.assay.qpcr`)
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
    - built-in protocol families should now be composed from shared internal
      figure building blocks (feature spans, strand-specific tails, linear
      molecule rows, event rows) rather than ad-hoc per-protocol struct
      literals; this keeps future PCR/Gibson growth on one composition model
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
  - `protocol` (built-in protocol cartoon id, for example `gibson.two_fragment`
    or `pcr.assay.pair` / `pcr.assay.pair.with_tail` / `pcr.assay.qpcr`)
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

Protocol-cartoon family growth direction (planned):

- Gibson specialist work now validates the abstraction-first protocol-cartoon
  strategy:
  - one protocol family should be expressed as a canonical
    `gentle.protocol_cartoon_template.v1` template plus deterministic bindings
  - the renderer should grow a shared collection of reusable figure building
    blocks that protocol families compose, instead of embedding protocol-
    specific drawing code in each built-in cartoon
  - protocol growth in count/shape (for example multi-fragment Gibson) should
    prefer repeated events, repeated molecules, and binding-level overrides
    rather than renderer-specific special cases
- Implemented baseline:
  - built-in Gibson cartoons now compose from shared internal building blocks
    for duplex spans, strand-specific tails, linear molecule rows, and event
    rows
  - this is intentionally still mechanism-first: Gibson cartoons describe
    fragment flow and achieved homology relationships, not full primer objects
    or low-level PCR parameter details
- PCR-assay cartoons should follow the same rule. The shared renderer should
  remain chemistry-agnostic and continue to render only ordered events,
  molecules, and features; PCR-specific meaning belongs in template structure
  and bindings, not in new PCR-only drawing primitives.
- PCR cartoon purpose:
  - explain assay intent and artifact flow, not every thermocycler sub-step
  - stay aligned with engine-owned operations, reports, and lineage-visible
    artifacts
  - keep lower-level primer sequences/thermocycler details in adjacent textual
    reports or inspectors rather than the cartoon itself
- Canonical PCR assay scene vocabulary should stay stable across modalities:
  - source template/context event
  - target/ROI event (selected span, feature-derived span, or queued region)
  - assay setup event (forward/reverse pair and optional probe or staged
    inner/outer sets may be named, but do not need literal primer glyphs)
  - amplification event
  - product/artifact event (amplicon, extracted copy, report, export, or
    explicit no-accepted-pairs outcome)
- Implemented PCR baseline in the `pcr.assay.*` family:
  - `pcr.assay.pair`: base strip with one selected ROI, one assay-setup lane,
    one amplification step, one amplicon/report outcome, and explicit forward/
    reverse primer glyphs with 5'/3' orientation
  - `pcr.assay.pair.no_product`: same family with an explicit report-only
    terminal state when no accepted primer pair yields a product
  - `pcr.assay.pair.with_tail`: insertion-first strip with requested extension
    sequences + insertion anchors, anchor-adjacent primer windows, and carried-
    through inserted terminal tails in the final amplicon
  - `pcr.oe.substitution`: six-step overlap-extension substitution strip with
    primer set `a`..`f`, first-step product haplotypes (`AB`/`CD`/`EF`),
    strand-specific anneal-gap geometry, and polymerase fill
- Implemented qPCR baseline in the same `pcr.assay.*` family:
  - `pcr.assay.qpcr`: same base strip enriched with one internal probe window
    plus explicit forward/reverse primer glyphs, a retained probe marker, and
    one quantitative readout terminal state
- Planned PCR modality adaptation should continue through the same
  `pcr.assay.*` protocol-cartoon family:
  - nested PCR: same family with two amplification stages (outer -> inner)
    instead of one, reusing the same event vocabulary
  - inverse PCR: same family with circular-template bindings and outward-facing
    primer semantics
  - batch/multiplex/tiling: repeated assay groups or repeated output lanes in
    bindings, not new renderer semantics per assay count
  - empty/failure outcomes: report/artifact nodes can render without product
    nodes when no accepted assay is produced
- Recommended rollout order:
  - extend the shipped PCR/qPCR baseline to queued batch PCR without changing
    renderer semantics
  - add nested, inverse, long-range, and multiplex variants as further
    template/binding expansions
- Naming/design rule:
  - do not introduce one built-in protocol id per assay count or minor UI view
  - prefer one stable protocol family with bindings that carry assay modality,
    stage count, molecule presence, and repeated-lane structure
  - keep generated explanatory strips exportable through the existing
    `protocol-cartoon ...` routes

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
  - `DbSnpFetch`
  - `RnaReadInterpret`
- Current cancellation support:
  - genome preparation supports cooperative cancellation plus optional
    `timeout_seconds` timeboxing and reports deterministic cancellation/timeout
    outcomes.
  - genome-track imports support cooperative cancellation and return partial
    import warnings.
  - dbSNP fetch currently emits staged progress events (`validate_input`,
    `inspect_prepared_genome`, `contact_server`, `wait_response`,
    `parse_response`, `resolve_placement`, `extract_region`,
    `attach_variant_marker`) but does not yet expose cooperative cancellation.
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
  - mutating artifact materialization per accepted pair:
    - one forward-primer sequence (`..._fwd`)
    - one reverse-primer sequence (`..._rev`)
    - one predicted amplicon sequence (`..._amplicon`) built from:
      full forward primer + template interior + reverse-primer reverse-complement
      (including non-annealing 5' tails)
    - one per-pair pool container containing all three artifacts
  - optional insertion-anchored context:
    - `insertion_context` (present when `DesignInsertionPrimerPairs` is used)
    - records requested forward/reverse anchor positions, extension sequences,
      window constraints, shift budget, and per-pair compensation rows
      (`forward_anchor_shift_bp`, `reverse_anchor_shift_bp`,
      compensation segments, and compensated primer/tail strings)

`DesignInsertionPrimerPairs` contract (implemented MVP):

- Purpose:
  - insertion-first wrapper around deterministic pair-primer design when the
    user already knows insert extensions and requested insertion anchors.
- Operation payload shape:

```json
{
  "DesignInsertionPrimerPairs": {
    "template": "seq_id",
    "insertion": {
      "requested_forward_3prime_end_0based_exclusive": 620,
      "requested_reverse_3prime_start_0based": 700,
      "forward_extension_5prime": "GAATTC",
      "reverse_extension_5prime": "CTCGAG",
      "forward_window_start_0based": 560,
      "forward_window_end_0based_exclusive": 650,
      "reverse_window_start_0based": 660,
      "reverse_window_end_0based_exclusive": 760,
      "max_anchor_shift_bp": 12
    },
    "forward": {
      "min_length": 20,
      "max_length": 30
    },
    "reverse": {
      "min_length": 20,
      "max_length": 30
    },
    "pair_constraints": {
      "require_roi_flanking": false
    },
    "min_amplicon_bp": 120,
    "max_amplicon_bp": 1200,
    "max_tm_delta_c": 2.0,
    "max_pairs": 200,
    "report_id": "tp73_insert_v1"
  }
}
```

- MVP behavior:
  - the insertion block is normalized first (IUPAC extension validation +
    anchor/window bounds checks)
  - forward/reverse primer windows are enforced from insertion windows
  - forward/reverse non-annealing tails are set from insertion extensions
  - primer design backend selection remains identical to `DesignPrimerPairs`
    (`auto|internal|primer3`)
  - resulting report is the same primer-report schema with populated
    `insertion_context` rows for shift/compensation inspection
  - no dedicated GUI form yet; operation is available through `op`/workflow
    payloads.

`PcrOverlapExtensionMutagenesis` contract (implemented baseline):

- Purpose:
  - deterministic overlap-extension insertion/deletion/replacement mutagenesis
    planning + staged product materialization in the main operation graph.
- Operation payload shape:

```json
{
  "PcrOverlapExtensionMutagenesis": {
    "template": "seq_id",
    "edit_start_0based": 620,
    "edit_end_0based_exclusive": 640,
    "insert_sequence": "GGTACC",
    "constraints": {
      "overlap_bp": 24,
      "outer_forward": {
        "min_length": 20,
        "max_length": 30
      },
      "outer_reverse": {
        "min_length": 20,
        "max_length": 30
      },
      "inner_forward": {
        "min_length": 18,
        "max_length": 28
      },
      "inner_reverse": {
        "min_length": 18,
        "max_length": 28
      }
    },
    "output_prefix": "tp73_oe_mut"
  }
}
```

- Baseline behavior:
  - `edit_start_0based..edit_end_0based_exclusive` defines the replaced region
    on the original template.
    - insertion: `edit_start == edit_end` and `insert_sequence` non-empty
    - deletion: `insert_sequence` empty and `edit_end > edit_start`
    - replacement: both deletion and insertion are non-empty
  - inner primers are chosen upstream/downstream of the edit and receive dynamic
    5' overlap tails derived from the mutant sequence so stage-1 products share
    one explicit overlap segment (minimum `overlap_bp`).
  - outer primers amplify both stage-1 fragments and the stage-2 final mutant
    amplicon.
  - operation materializes graph-visible artifacts:
    - primers: `..._outer_fwd`, `..._outer_rev`, `..._inner_fwd`, `..._inner_rev`
    - stage-1 products: `..._stage1_left`, `..._stage1_right`
    - final stage-2 mutant: `..._mutant`
    - three per-stage pool containers (left, right, final)
  - operation warnings include deterministic candidate-search limit notices when
    the combinatorial search budget is exhausted.
  - insertion/replacement runs now also emit
    `OpResult.protocol_cartoon_preview` for built-in protocol
    `pcr.oe.substitution`, including deterministic
    `flank_bp`/`overlap_bp`/`insert_bp` geometry and bound template overrides
    (`gentle.protocol_cartoon_template_bindings.v1`) for adapter rendering.

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

Feature-query shell contract (implemented):

- Shared-shell command:
  - `features query SEQ_ID [--kind KIND] [--kind-not KIND] [--range START..END|--start N --end N] [--overlap|--within|--contains] [--strand any|forward|reverse] [--label TEXT] [--label-regex REGEX] [--qual KEY] [--qual-contains KEY=VALUE] [--qual-regex KEY=REGEX] [--min-len N] [--max-len N] [--limit N] [--offset N] [--sort feature_id|start|end|kind|length] [--desc] [--include-source] [--include-qualifiers]`
- Execution semantics:
  - non-mutating engine inspection over one sequence’s feature table
  - deterministic filter pipeline:
    kind include/exclude, optional range relation (`overlap|within|contains`),
    strand filter, label contains/regex, qualifier filters, and length bounds
  - deterministic ordering by requested sort key with stable tie-breaks +
    pagination (`offset`/`limit`)
- Response schema:
  - `gentle.sequence_feature_query_result.v1`
  - fields include:
    - `seq_id`, `sequence_length_bp`, `total_feature_count`
    - `matched_count`, `returned_count`, `offset`, `limit`
    - normalized `query`
    - `rows[]` with `feature_id`, `kind`, `start_0based`,
      `end_0based_exclusive`, `length_bp`, `strand`, `label`, `labels[]`, and
      optional qualifier maps when requested (`--include-qualifiers`)

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
  - `SummarizeRnaReadGeneSupport { report_id, gene_ids, selected_record_indices?, complete_rule?, path? }`
  - `InspectRnaReadGeneSupport { report_id, gene_ids, selected_record_indices?, complete_rule?, cohort_filter?, path? }`
  - `seed_feature_id` may reference an `mRNA`, `transcript`, `ncRNA`,
    `misc_RNA`, `exon`, `gene`, or `CDS` feature; transcript-template
    admission then follows the selected splicing-scope rules around that seed.
  - implemented profile: `nanopore_cdna_v1`
  - implemented input format: `fasta` (`.fa/.fasta`, optional `.fa.gz/.fasta.gz`; `.sra` must be converted externally in phase-1)
  - default seed/filter constants:
    - `kmer_len=10`
    - `seed_stride_bp=1`
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
    - seed-start density is controlled by `seed_stride_bp`
    - default density is one start per base (`seed_stride_bp=1`)
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
    - `report_mode=seed_passed_only` keeps a smaller retained subset for later
      inspection/alignment:
      - retained hits that passed the composite seed gate
      - retained hits at or above raw `min hit`
      - counters still remain based on the full stream
    - `checkpoint_path` + `checkpoint_every_reads` writes deterministic JSON
      snapshots (`gentle.rna_read_interpret_checkpoint.v1`) during streaming
    - `resume_from_checkpoint=true` resumes from the checkpoint snapshot and
      fast-forwards already-processed records deterministically
  - phase-2 alignment behavior:
    - `AlignRnaReadReport` loads a persisted report and reprocesses a selected
      retained subset (`all|seed_passed|aligned`)
    - phase-2 progress events now emit once per selected retained row
      (`update_every_reads=1`) so adapters can show visible row-by-row advance
    - GUI/shared-shell default selection is `seed_passed`
    - optional `selected_record_indices[]` (0-based stored `record_index`)
      overrides the selection preset and aligns only the explicit subset
    - `selection=all` remains available when you deliberately want the broader
      rescued-retained working set to receive round-2 similarity/coverage
      scores
    - `selection=aligned` means rerun phase 2 only on retained rows that
      already have a stored mapping from an earlier phase-2 pass
    - if `selection=seed_passed` matches no retained hits and no explicit
      record indices were supplied, the engine falls back to retained rows at
      or above raw `min hit`, and if that is still empty, to the highest
      phase-1 score retained row
    - aligner configuration uses `align_config_override` when supplied,
      otherwise the report-stored `align_config`
    - mapping backend uses `bio::alignment::pairwise::banded` with
      `align_band_bp` as band width (`w`) and transcript-seed `kmer_len` as
      seed length (`k`), plus deterministic dense fallback when the banded
      solver yields no mapping
    - phase-2 pairwise alignment evaluates both query orientations for every
      retained row (stored query plus reverse complement) and keeps the
      best-scoring deterministic candidate, preferring semiglobal over local
      and non-reversed over reversed only as later tie-breakers
    - selected retained rows are pairwise aligned regardless of whether their
      recomputed composite seed-pass flag remains true; the seed-pass result is
      still recomputed and stored independently for later inspection
    - updated report fields include:
      - per-hit mapping fields (`best_mapping`, `secondary_mappings`)
      - per-hit `msa_eligible` and `msa_eligibility_reason`
      - aggregate `read_count_aligned` and `retained_count_msa_eligible`
      - refreshed seed/path diagnostics
        (`transition_support_rows`, `isoform_support_rows`)
      - refreshed mapped support rows
        (`exon_support_frequencies`, `junction_support_frequencies`,
        `mapped_isoform_support_rows`)
      - mapped exon/junction support is derived from aligned transcript-template
        offsets first and falls back to coarse genomic-span overlap only for
        legacy mappings that do not carry template offsets
      - deterministic retained-hit re-ranking by alignment-aware retention rank
  - alignment inspection behavior:
    - `rna-reads inspect-alignments` accepts coarse `selection` plus a
      structured subset contract:
      - `effect_filter = all_aligned|confirmed_only|disagreement_only|reassigned_only|no_phase1_only|selected_only`
      - `sort_key = rank|identity|coverage|score`
      - `search = free-text match over read ids, transcript ids/labels,
        effect labels, and `#record_index` labels`
      - `selected_record_indices[]` provides the explicit subset for
        `selected_only`
      - `score_density_variant = all_scored|composite_seed_gate|retained_replay_current_controls`
      - optional `score_density_seed_filter_override` carries the current
        seed-gate controls when an adapter requests retained-only replay under
        current controls
      - `score_bin_index` + `score_bin_count` provide a formal
        score-density-bin subset for reproducible histogram-driven inspection
        within that chosen histogram population
    - inspection payload now includes:
      - `aligned_count`: aligned rows admitted by coarse `selection`
      - `subset_match_count`: aligned rows matching the structured subset
        before `limit`
      - `row_count`: returned rows after `limit`
      - `subset_spec`: normalized structured subset object echoed back in the
        response for deterministic replay
    - row `rank` remains the original alignment-aware retention rank even when
      subset sorting reorders the returned rows
  - on-demand pairwise-alignment detail behavior:
    - the engine can reconstruct the exact phase-2 read-vs-transcript-template
      alignment for one retained row from the saved report plus admitted
      transcript-template set
    - detail payload schema:
      - `gentle.rna_read_alignment_detail.v1`
    - payload includes:
      - selected retained row id (`record_index`, `header_id`)
      - transcript/template target identity
      - phase-2 `alignment_mode`
      - alignment backend (`banded` or `dense_fallback`)
      - aligned query/template spans, full template length, score, identity,
        query coverage, transcript-template coverage, and CIGAR
      - aligned `query / relation / target` text rows for manual inspection of
        low-complexity or partial confirmations
  - exact-subset export behavior:
    - `ExportRnaReadHitsFasta`, `ExportRnaReadExonPathsTsv`,
      `ExportRnaReadExonAbundanceTsv`, and `ExportRnaReadAlignmentsTsv` accept
      optional `selected_record_indices[]`
    - when present, the explicit 0-based stored `record_index` subset
      overrides the coarse `selection` preset
    - these exports also accept optional `subset_spec`, a human-readable formal
      description such as `filter=... | sort=... | search=...`; when provided,
      the exported artifact records both the explicit `record_index` subset and
      the subset definition that produced it
    - intended for exporting the exact contributor reads surfaced by mapped
      `Audit` actions in the GUI
  - target-gene cohort summary behavior:
    - `SummarizeRnaReadGeneSupport` is non-mutating and consumes one persisted
      aligned RNA-read report
    - required `gene_ids[]` are normalized/deduplicated and matched
      case-insensitively against the same splicing group-label logic already
      used for transcript grouping
    - output schema:
      - `gentle.rna_read_gene_support_summary.v1`
    - base cohort:
      - retained rows with `best_mapping` present
      - optionally intersected with explicit `selected_record_indices[]`
    - accepted target cohort:
      - base-cohort rows whose `best_mapping.transcript_feature_id` resolves to
        one of the requested matched genes/groups
    - complete/fragment split:
      - `complete_rule = near|strict|exact` controls which accepted rows land
        in the `complete` cohort
      - `fragment` is the remaining accepted-target cohort
      - summary still reports nested `complete_strict_count` and
        `complete_exact_count` regardless of the chosen `complete_rule`
    - support attribution is derived from phase-2 mapped support, not phase-1
      `exon_path`
    - per-cohort output blocks:
      - `all_target`
      - `fragments`
      - `complete`
    - each block includes:
      - `read_count`
      - `exon_support[]`
      - `exon_pair_support[]`
      - `direct_transition_support[]`
    - row semantics:
      - `exon_support[]`: each exon counted at most once per read
      - `exon_pair_support[]`: every ordered exon_i -> exon_j pair observed in
        the mapped exon order once per read, including skipped pairs like
        `1->3`
      - `direct_transition_support[]`: neighboring exon steps only, so
        `1->2` is counted but skipped pairs like `1->3` are not
      - all support fractions are normalized by the enclosing cohort size
      - exon and pair rows carry deterministic gene-level exon ordinals plus
        genomic coordinates for auditability
    - when `path` / shell `--output` is provided, the exact same JSON payload
      returned to the caller is also written to disk
  - target-gene cohort audit behavior:
    - `InspectRnaReadGeneSupport` is non-mutating and shares the same
      requested-gene matching, selected-record restriction, accepted-target
      logic, and `complete_rule` classification used by
      `SummarizeRnaReadGeneSupport`
    - output schema:
      - `gentle.rna_read_gene_support_audit.v1`
    - evaluation universe:
      - all selected saved-report rows after `selected_record_indices[]`
        filtering, including unaligned retained rows
    - grouped top-level subset handles:
      - `accepted_target_record_indices[]`
      - `fragment_record_indices[]`
      - `complete_record_indices[]`
      - `complete_strict_record_indices[]`
      - `complete_exact_record_indices[]`
    - row status values:
      - `unaligned`
      - `aligned_other_gene`
      - `accepted_fragment`
      - `accepted_complete`
    - row payload includes:
      - `record_index`, `header_id`
      - resolved `gene_id` when available
      - aligned transcript identity (`transcript_feature_id`,
        `transcript_id`, `transcript_label`)
      - machine-readable `status_reason`
      - `full_length_exact`, `full_length_near`, `full_length_strict`, and
        derived `full_length_class`
      - `mapped_exon_ordinals[]`
      - ordered `exon_pairs[]`
      - ordered `direct_transition_pairs[]`
      - phase-2 `score`, `identity_fraction`, `query_coverage_fraction`
      - `passed_seed_filter` as provenance only
    - `cohort_filter = all|accepted|fragment|complete|rejected` limits the
      returned `rows[]` set without changing the grouped top-level subset
      arrays
    - when `path` / shell `--output` is provided, the exact same JSON payload
      returned to the caller is also written to disk
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
    - `score_density_bins[]` (`all_scored` phase-1 histogram)
    - `seed_pass_score_density_bins[]` (`composite_seed_gate` histogram)
    - exact read-length histograms (`length_bp -> count`) for
      deterministic subset auditing:
      - `read_length_counts_all`
      - `read_length_counts_seed_passed`
      - `read_length_counts_aligned`
      - `read_length_counts_full_length_exact`
      - `read_length_counts_full_length_near`
      - `read_length_counts_full_length_strict`
      - checkpoint snapshots mirror these vectors so resume/restart keeps
        histogram accumulation deterministic
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
    - `best_mapping.query_reverse_complemented` (whether phase-2 had to
      reverse-complement the stored read to fit the chosen transcript-template
      mapping)
  - alignment inspection payload schema:
    - `gentle.rna_read_alignment_inspection.v1`
    - produced by non-mutating shared-shell inspection command
      `rna-reads inspect-alignments`
    - each row now carries:
      - phase-1 transcript-assignment fields
        (`phase1_primary_transcript_id`, `seed_chain_transcript_id`,
        `exon_path_transcript_id`, `exon_path`,
        `exon_transitions_confirmed/total`, `selected_strand`,
        `reverse_complement_applied`)
      - phase-2 best-mapping fields
        (`transcript_id`, `transcript_label`, `strand`, `alignment_mode`,
        `target_start_1based`, `target_end_1based`, `target_length_bp`,
        `identity_fraction`, `query_coverage_fraction`,
        `target_coverage_fraction`, `score`, `secondary_mapping_count`)
      - full-length classification flags (derived deterministically from
        transcript-template coverage and current alignment threshold):
        - `full_length_exact` (`100%` template coverage)
        - `full_length_near` (`>=95%` template coverage)
        - `full_length_strict`
          (`near` + both template ends within `15 bp` + identity above
          active alignment threshold)
      - deterministic comparison field `alignment_effect`
        (`confirmed_assignment`, `reassigned_transcript`,
        `aligned_without_phase1_assignment`)
      - mapped-support attribution arrays for the best mapping
        (`mapped_exon_support[]`, `mapped_junction_support[]`)
    - top-level payload now also carries:
      - `aligned_count`
      - `subset_match_count`
      - `row_count`
      - `limit`
      - normalized `subset_spec`
        (`effect_filter`, `sort_key`, `search`,
        `selected_record_indices[]`, `score_density_variant`,
        `score_bin_index`, `score_bin_count`)
- Sample-sheet export:
  - operation: `ExportRnaReadSampleSheet { path, seq_id?, report_ids?, append? }`
  - export schema: `gentle.rna_read_sample_sheet_export.v1`
  - output: TSV with run/read metrics, sparse-origin request provenance
    (`report_mode`, `origin_mode`, `target_gene_count`,
    `target_gene_ids_json`, `roi_seed_capture_enabled`), JSON-serialized
    exon/junction frequency columns, and `origin_class_counts_json` for
    cohort-level downstream analysis.
- Shared-shell command family:
  - `rna-reads interpret SEQ_ID FEATURE_ID INPUT.fa[.gz] [--report-id ID] [--report-mode full|seed_passed_only] [--checkpoint-path PATH] [--checkpoint-every-reads N] [--resume-from-checkpoint|--no-resume-from-checkpoint] [--profile nanopore_cdna_v1] [--format fasta] [--scope all_overlapping_both_strands|target_group_any_strand|all_overlapping_target_strand|target_group_target_strand] [--origin-mode single_gene|multi_gene_sparse] [--target-gene GENE_ID]... [--roi-seed-capture|--no-roi-seed-capture] [--kmer-len N] [--seed-stride-bp N] [--min-seed-hit-fraction F] [--min-weighted-seed-hit-fraction F] [--min-unique-matched-kmers N] [--min-chain-consistency-fraction F] [--max-median-transcript-gap F] [--min-confirmed-transitions N] [--min-transition-support-fraction F] [--cdna-poly-t-flip|--no-cdna-poly-t-flip] [--poly-t-prefix-min-bp N] [--align-band-bp N] [--align-min-identity F] [--max-secondary-mappings N]`
  - `rna-reads align-report REPORT_ID [--selection all|seed_passed|aligned] [--record-indices i,j,k] [--align-band-bp N] [--align-min-identity F] [--max-secondary-mappings N]`
  - `rna-reads list-reports [SEQ_ID]`
  - `rna-reads show-report REPORT_ID`
  - `rna-reads summarize-gene-support REPORT_ID --gene GENE_ID [--gene GENE_ID ...] [--record-indices i,j,k] [--complete-rule near|strict|exact] [--output PATH]`
  - `rna-reads inspect-gene-support REPORT_ID --gene GENE_ID [--gene GENE_ID ...] [--record-indices i,j,k] [--complete-rule near|strict|exact] [--cohort all|accepted|fragment|complete|rejected] [--output PATH]`
  - `rna-reads inspect-alignments REPORT_ID [--selection all|seed_passed|aligned] [--limit N] [--effect-filter all_aligned|confirmed_only|disagreement_only|reassigned_only|no_phase1_only|selected_only] [--sort rank|identity|coverage|score] [--search TEXT] [--record-indices i,j,k] [--score-bin-variant all_scored|composite_seed_gate] [--score-bin-index N] [--score-bin-count M]`
  - `rna-reads export-report REPORT_ID OUTPUT.json`
  - `rna-reads export-hits-fasta REPORT_ID OUTPUT.fa [--selection all|seed_passed|aligned] [--record-indices i,j,k] [--subset-spec TEXT]`
  - `rna-reads export-sample-sheet OUTPUT.tsv [--seq-id ID] [--report-id ID]... [--append]`
  - `rna-reads export-paths-tsv REPORT_ID OUTPUT.tsv [--selection all|seed_passed|aligned] [--record-indices i,j,k] [--subset-spec TEXT]`
  - `rna-reads export-abundance-tsv REPORT_ID OUTPUT.tsv [--selection all|seed_passed|aligned] [--record-indices i,j,k] [--subset-spec TEXT]`
  - `rna-reads export-score-density-svg REPORT_ID OUTPUT.svg [--scale linear|log] [--variant all_scored|composite_seed_gate]`
  - `rna-reads export-alignments-tsv REPORT_ID OUTPUT.tsv [--selection all|seed_passed|aligned] [--limit N] [--record-indices i,j,k] [--subset-spec TEXT]`
  - `rna-reads export-alignment-dotplot-svg REPORT_ID OUTPUT.svg [--selection all|seed_passed|aligned] [--max-points N]`
  - shell output convenience fields:
    - `rna-reads list-reports` includes `summary_rows[]` with concise
      human-readable provenance lines (`mode`, `origin`, target count,
      ROI-capture flag, read counters)
    - `rna-reads show-report` includes `summary` with the same provenance
      framing for one report
    - `rna-reads summarize-gene-support` returns the full
      `gentle.rna_read_gene_support_summary.v1` payload directly, including
      `requested_gene_ids`, `matched_gene_ids`, `missing_gene_ids`,
      selected-record echo fields, and per-cohort support tables
    - `rna-reads inspect-gene-support` returns the full
      `gentle.rna_read_gene_support_audit.v1` payload directly, including
      grouped cohort record-index arrays plus row-level `status`,
      `status_reason`, full-length fields, and mapped exon/junction audit data
    - `rna-reads inspect-alignments` returns aligned rows ranked by
      alignment-aware retention score (mapping + seed metrics), plus a
      structured `subset_spec` payload (`effect_filter`, `sort_key`, `search`,
      `selected_record_indices`, `score_density_variant`, `score_bin_index`,
      `score_bin_count`) and
      `subset_match_count`
- Alignment-TSV export:
  - operation:
    `ExportRnaReadAlignmentsTsv { report_id, path, selection, limit?, selected_record_indices?, subset_spec? }`
  - export schema: `gentle.rna_read_alignment_tsv_export.v1`
  - output: ranked alignment rows as TSV with:
    - leading `#` metadata lines for report provenance (`selection`, `limit`,
      `selected_record_indices`, `subset_spec`, `profile`, `scope`, `origin_mode`)
    - seed-screen sampling/gating context (`k`, `seed_stride_bp`,
      overlap/order-density wording, seed thresholds)
    - alignment config summary (`min_identity_fraction`,
      `max_secondary_mappings`)
    - phase-1 transcript/path diagnostics
    - phase-2 mapping metrics
    - `alignment_effect`
    - compact mapped exon/junction attribution columns
    - optional top-`N` truncation via `limit`
- Score-density SVG export:
  - `rna-reads export-score-density-svg` writes the same report summary used by
    the GUI plus seed-screen provenance in the SVG header:
    - `variant = all_scored|composite_seed_gate|retained_replay_current_controls`
    - `profile`, `report_mode`, `scope`, `origin_mode`
    - seed-filter summary with `k`, `seed_stride_bp`, thresholds, and
      overlap/order-density wording
    - optional `replay_seed_filter` summary when the export uses retained-only
      replay under current controls
    - whether bins were stored in the report or derived from retained hits
- Alignment-dotplot export:
  - operation:
    `ExportRnaReadAlignmentDotplotSvg { report_id, path, selection, max_points }`
  - export schema: `gentle.rna_read_alignment_dotplot_svg_export.v1`
  - output: SVG scatter of query coverage vs identity for aligned hits with
    score-colored points and report-config threshold guide.
- Read-sequence materialization:
  - operation:
    `MaterializeRnaReadHitSequences { report_id, selection, selected_record_indices?, output_prefix? }`
  - output:
    - creates one ordinary project sequence per selected retained RNA-read hit
    - exact `selected_record_indices` takes precedence over coarse `selection`
    - intended for downstream dotplots/manual inspection of saved-report
      outliers without re-reading the FASTA input
- `rna-reads export-hits-fasta` header extensions:
  - optional `selected_record_indices[]` overrides the coarse selection preset
    for exact saved-report subset export
  - optional `subset_spec` records the formal subset definition that produced
    that explicit `record_index` subset
  - `exon_path_tx=<transcript_id|none>`
  - `exon_path=<ordinal_path|none>` using `:` for hash-confirmed adjacent
    exon transitions and `-` for unconfirmed adjacency
  - `exon_transitions=<confirmed>/<total>`
  - `rc_applied=<true|false>` (automatic cDNA poly-T reverse-complement
    normalization marker)
  - `origin_class=<...>` plus `origin_conf=<...>` and `strand_conf=<...>`
- `rna-reads export-exon-paths-tsv` and `rna-reads export-exon-abundance-tsv`
  now begin with the same `#` report/seed-screen provenance block used by the
  alignment TSV export, minus alignment-only fields; optional `subset_spec`
  records the formal subset definition alongside `selected_record_indices`
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
