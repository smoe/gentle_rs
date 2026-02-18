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

## Core entities

### ProjectState

```json
{
  "sequences": {"seq_id": "DNAsequence object"},
  "metadata": {"any": "json"},
  "display": {"ui_visibility_and_tfbs_filter_state": "..."},
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
- `Digest { input, enzymes, output_prefix? }`
- `Ligation { inputs, circularize_if_possible, protocol, output_id?, output_prefix?, unique? }`
- `MergeContainers { inputs, output_prefix? }`
- `Pcr { template, forward_primer, reverse_primer, output_id?, unique? }`
- `PcrAdvanced { template, forward_primer, reverse_primer, output_id?, unique? }`
- `PcrMutagenesis { template, forward_primer, reverse_primer, mutations, output_id?, unique?, require_all_mutations? }`
- `ExtractRegion { input, from, to, output_id? }`
- `SelectCandidate { input, criterion, output_id? }`
- `FilterByMolecularWeight { inputs, min_bp, max_bp, error, unique, output_prefix? }`
- `Reverse { input, output_id? }`
- `Complement { input, output_id? }`
- `ReverseComplement { input, output_id? }`
- `Branch { input, output_id? }`
- `SetDisplayVisibility { target, visible }`
- `SetTopology { seq_id, circular }`
- `RecomputeFeatures { seq_id }`
- `SetParameter { name, value }` (purely in-silico project parameter change)

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
