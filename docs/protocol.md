# GENtle Engine Protocol (Draft v1)

This document defines the draft machine-facing protocol for operating GENtle
through a shared core engine.

Goal:

- GUI, CLI, JavaScript, and Lua call the same core routines.
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
- `RenderRnaStructureSvg { seq_id, path }`
- `RenderLineageSvg { path }`
- `RenderPoolGelSvg { inputs, path, ladders? }`
- `ExportDnaLadders { path, name_filter? }`
- `ExportRnaLadders { path, name_filter? }`
- `Digest { input, enzymes, output_prefix? }`
- `Ligation { inputs, circularize_if_possible, protocol, output_id?, output_prefix?, unique? }`
- `MergeContainers { inputs, output_prefix? }`
- `Pcr { template, forward_primer, reverse_primer, output_id?, unique? }`
- `PcrAdvanced { template, forward_primer, reverse_primer, output_id?, unique? }`
- `PcrMutagenesis { template, forward_primer, reverse_primer, mutations, output_id?, unique?, require_all_mutations? }`
- `ExtractRegion { input, from, to, output_id? }`
- `ExtendGenomeAnchor { seq_id, side, length_bp, output_id?, catalog_path?, cache_dir? }`
- `SelectCandidate { input, criterion, output_id? }`
- `GenerateCandidateSet { set_name, seq_id, length_bp, step_bp, feature_kinds[], feature_label_regex?, max_distance_bp?, feature_geometry_mode?, feature_boundary_mode?, feature_strand_relation?, limit? }`
- `GenerateCandidateSetBetweenAnchors { set_name, seq_id, anchor_a, anchor_b, length_bp, step_bp, limit? }`
- `DeleteCandidateSet { set_name }`
- `ScoreCandidateSetExpression { set_name, metric, expression }`
- `ScoreCandidateSetDistance { set_name, metric, feature_kinds[], feature_label_regex?, feature_geometry_mode?, feature_boundary_mode?, feature_strand_relation? }`
- `FilterCandidateSet { input_set, output_set, metric, min?, max?, min_quantile?, max_quantile? }`
- `CandidateSetOp { op: union|intersect|subtract, left_set, right_set, output_set }`
- `ScoreCandidateSetWeightedObjective { set_name, metric, objectives[], normalize_metrics? }`
- `TopKCandidateSet { input_set, output_set, metric, k, direction?, tie_break? }`
- `ParetoFrontierCandidateSet { input_set, output_set, objectives[], max_candidates?, tie_break? }`
- `UpsertWorkflowMacroTemplate { name, description?, parameters[], script }`
- `DeleteWorkflowMacroTemplate { name }`
- `UpsertCandidateMacroTemplate { name, description?, parameters[], script }`
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

`ExtendGenomeAnchor` side semantics:

- `side` accepts `five_prime` or `three_prime`.
- Direction is contextual to anchor strand.
- On anchor strand `-`, `five_prime` increases physical genomic position.

Local `SequenceAnchor` semantics (distinct from genome provenance anchoring):

- `SequenceAnchor` currently supports:
  - `Position { zero_based }`
  - `FeatureBoundary { feature_kind?, feature_label?, boundary, occurrence? }`
- `boundary` accepts `Start`, `End`, or `Middle`.
- This anchor model resolves in-sequence positions and is used for
  in-silico extraction/scoring workflows (`ExtractAnchoredRegion`,
  `GenerateCandidateSetBetweenAnchors`).

Adapter utility contracts (current, non-engine operations):

- `help [COMMAND ...] [--format text|json|markdown] [--interface ...]`
  - backed by structured glossary source `docs/glossary.json`
  - `--format text` renders human-readable help
  - `--format json` renders machine-readable help catalog/topic payload
  - `--format markdown` renders documentation-ready markdown

- `macros run/template-list/template-show/template-put/template-delete/template-run`
  - shared-shell macro adapter family for full operation/workflow scripting
  - template persistence is backed by engine operations
    `UpsertWorkflowMacroTemplate`/`DeleteWorkflowMacroTemplate`
  - expanded scripts can execute `op ...` and `workflow ...` statements and
    optionally roll back via `--transactional`

- `screenshot-window OUTPUT.png`
  - currently disabled by security policy
  - returns deterministic disabled message from shared shell/CLI/GUI command
    paths
  - kept as reserved adapter contract for future re-enable after explicit
    endpoint-security approval

- `agents list [--catalog PATH]`
  - Lists configured agent systems from catalog JSON.
  - Default catalog: `assets/agent_systems.json`.

- `agents ask SYSTEM_ID --prompt TEXT [--catalog PATH] [--allow-auto-exec] [--execute-all] [--execute-index N ...] [--no-state-summary]`
  - Invokes one configured agent system via catalog transport.
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
    }
  ]
}
```

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
- `feature_details_font_size` (default `11.0`, range `8.0..24.0`)
  - controls GUI font size for the feature tree entries and feature range details
- `regulatory_feature_max_view_span_bp` (default `50000`, range `>= 0`)
  - hides regulatory feature overlays in linear view when current view span
    exceeds this threshold (`0` disables regulatory overlays)
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
  - template expansion/binding is exposed through adapter command surfaces
    (`macros template-*`)
  - expanded scripts run through shared shell execution (`macros run`) and can
    orchestrate full cloning operations via `op ...` or `workflow ...` payloads
- Candidate macro templates are persisted in project metadata:
  - `UpsertCandidateMacroTemplate` stores/replaces named templates
  - `DeleteCandidateMacroTemplate` removes templates
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
- Current cancellation support:
  - genome-track imports support cooperative cancellation and return partial
    import warnings.

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

### Workflow

```json
{
  "run_id": "string",
  "ops": ["Operation", "Operation", "..."]
}
```

### OpResult

```json
{
  "op_id": "op-1",
  "created_seq_ids": ["..."],
  "changed_seq_ids": ["..."],
  "warnings": ["..."],
  "messages": ["..."]
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
- CRISPR guide-design base layer:
  - guide-candidate domain model and ranking layer
  - oligo generation and export contracts (table/plate/protocol text)
  - macro/template expansion into deterministic `Workflow` JSON
  - design-constraint filtering is already available as `FilterByDesignConstraints`
  - see draft: `docs/rna_guides_spec.md`
