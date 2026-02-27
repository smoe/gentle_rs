# Cloning Routine Catalog and Macro-Graph Plan

Last updated: 2026-02-26

Purpose: map current cloning vocabulary coverage to GENtle capabilities and
define an implementation plan for a first-class cloning-routine catalog where
macro instances render as workflow graph boxes with explicit inputs/outputs.

Reference vocabulary source: cloning-method review article
`S1016847825000585` (Elsevier/ScienceDirect) and linked PubMed/PMC record
`PMID: 40449799`.

## 1. Coverage crosswalk (current GENtle vs core cloning vocabulary)

| Routine family / vocabulary | Current GENtle support | Gap status |
|---|---|---|
| Restriction digest + ligation (`restriction enzyme`, `sticky end`, `blunt end`, `ligation`) | Implemented primitives: `Digest`, `Ligation` (`Sticky`/`Blunt`), container-aware variants, starter macro `digest_ligate_extract_sticky` | Partial: no curated full template family (single-insert, multi-insert, directional variants) |
| PCR cloning (`primer`, `amplicon`, site introduction) | Implemented primitives: `Pcr`, `PcrAdvanced`, `PcrMutagenesis`, starter macro `pcr_site_insertion_then_digest` | Partial: no protocol-labeled template catalog coverage and no workflow-graph macro nodes |
| TA/GC and TOPO cloning (`A-overhang`, `T-overhang`, `TOPO`) | No protocol-specific operation/template set | Missing |
| Gibson Assembly (`homology overlap`) | No dedicated overlap-assembly operation; only approximable through generic steps/macros | Missing |
| In-Fusion / NEBuilder HiFi (`recombination-style overlap assembly`) | No dedicated operation/template semantics for overlap validation and assembly intent | Missing |
| Golden Gate (`Type IIS`, modular assembly) | No Type IIS-aware protocol template or overhang compatibility planner | Missing |
| Gateway recombination (`BP`, `LR`, `att` sites) | No recombination-specific operation/template family | Missing |
| Clone screening/validation vocabulary (`colony PCR`, `Sanger verification`) | Generic sequence operations and render/export exist, but no protocol templates for screening flows | Partial/missing |

Current reusable infrastructure already in place:

- Shared operations in `src/engine.rs` and adapter parity via `src/engine_shell.rs`.
- Workflow macro persistence (`UpsertWorkflowMacroTemplate`,
  `DeleteWorkflowMacroTemplate`) with schema
  `gentle.cloning_macro_template.v1`.
- Importable starter pack `assets/cloning_patterns.json` with schema
  `gentle.cloning_patterns.v1`.
- Hierarchical per-template catalog at `assets/cloning_patterns_catalog`:
  - one macro template per file (`gentle.cloning_pattern_template.v1`)
  - directory hierarchy doubles as GUI selection hierarchy
  - shared shell import supports directory-recursive ingest via
    `macros template-import assets/cloning_patterns_catalog`
- Graph rendering currently supports sequence nodes (including pooled sequence
  render variants) and arrangement nodes, but not explicit macro-instance
  nodes.

Current starter hierarchy (implemented):

- `assets/cloning_patterns_catalog/restriction/digest_ligation/digest_ligate_extract_sticky.json`
- `assets/cloning_patterns_catalog/sequence/transform/branch_reverse_complement.json`
- `assets/cloning_patterns_catalog/pcr/site_insertion/pcr_site_insertion_then_digest.json`
- `assets/cloning_patterns_catalog/crispr/guides/candidate_scans/grna_candidate_priority_scan.json`
- `assets/cloning_patterns_catalog/crispr/guides/candidate_scans/grna_anchor_window_scan.json`
- `assets/cloning_patterns_catalog/crispr/guides/oligos/grna_practical_filter_and_oligos.json`

## 2. Target outcome

1. Users can choose cloning routines from a curated catalog (GUI/CLI/shell/JS/Lua parity).
2. Every routine has a typed input/output contract.
3. Each macro execution is represented as a box node in lineage/workflow graph
   with deterministic edges from inputs to outputs.
4. Multiple instances of the same routine can coexist in one project and remain
   distinguishable/auditable.

## 3. Catalog model (proposed)

Add a catalog manifest (`assets/cloning_routines.json`, versioned schema) with
entries like:

- `routine_id` (stable machine key; e.g. `restriction_single_insert`)
- `title` (user-facing label)
- `family` (`restriction`, `gibson`, `golden_gate`, `gateway`, `topo`, etc.)
- `vocabulary_tags` (search terms/synonyms)
- `summary`
- `details_url` (optional external reference)
- `template_name` (backing workflow macro template)
- `input_ports[]` (typed ports)
- `output_ports[]` (typed ports)
- `status` (`implemented`, `partial`, `planned`)

Typed port fields:

- `port_id`
- `kind` (`sequence`, `container`, `candidate_set`, `guide_set`, `string`, `number`, `bool`, `path`)
- `required` (`true/false`)
- `cardinality` (`one`, `many`)
- `description`

## 4. Macro-instance graph box model (proposed)

Add a macro-run record to project metadata/state:

- `macro_instance_id` (deterministic, e.g. `macro:<run_id>:<ordinal>`)
- `routine_id`
- `template_name`
- `run_id`
- `created_at_unix_ms`
- `bound_inputs` (resolved port -> concrete ids/values)
- `bound_outputs` (resolved port -> concrete ids/values)
- `expanded_op_ids` (ordered op ids emitted by the run)
- `status` (`ok`, `failed`, `cancelled`)

Graph behavior:

- render macro instances as rectangular box nodes.
- add edges `input sequence/container -> macro box -> output sequence/container`.
- preserve existing sequence and arrangement nodes; macro nodes are additive.
- allow repeated instances of same template by unique `macro_instance_id`.

## 5. Delivery plan

### Phase 1: Catalog foundation (read-only)

- Add schema + loader for cloning routine catalog.
- Add shared-shell command to list/search routines by family/tag/status.
- Seed catalog from existing starter macros and mark missing families as
  `planned`.

### Phase 2: Typed I/O contract for workflow macros

- Extend workflow macro template payload with optional explicit
  `input_ports`/`output_ports`.
- Validate `template-run` bindings against required inputs and type/cardinality
  constraints.
- Keep backward compatibility for existing script-only templates.

### Phase 3: Macro run recording

- Persist macro-run instance records after `macros run` and
  `macros template-run`.
- Capture concrete resolved inputs/outputs and emitted operation ids.
- Add deterministic tests for replay-safe recording.

### Phase 4: Graph box rendering

- Extend lineage graph model to include macro-instance node kind.
- Render macro nodes as boxes with label `routine/title + instance id`.
- Add hover/details panel with ports and emitted op list.
- Ensure multiple instances of same macro render independently.

### Phase 5: Catalog population and protocol packs

- Ship first complete routine families:
  - restriction cloning (single/multi/directional)
  - PCR-cloning starter pack
  - Golden Gate baseline
  - Gibson baseline
- Add remaining families (`Gateway`, `TOPO`, `TA/GC`, `In-Fusion`,
  `NEBuilder HiFi`) as additive packs.

## 6. Test and acceptance criteria

- Deterministic unit tests for:
  - catalog load/validate/sort/filter
  - macro template I/O contract validation
  - macro-run instance recording with stable IDs
  - graph-node/edge derivation from macro-run records
- At least one integration test that runs the same routine twice in one project
  and confirms two separate macro box instances with correct edges.
- `cargo check -q` and targeted tests for modified engine/shell/graph paths.
