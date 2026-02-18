# GENtle CLI Manual

This page documents command-line entry points for GENtle.

## Overview

GENtle currently provides three binaries:

- `gentle`: graphical desktop app
- `gentle_js`: interactive JavaScript shell
- `gentle_lua`: interactive Lua shell
- `gentle_cli`: JSON operation/workflow CLI for automation and AI tools

## Build and run

From the repository root:

```bash
cargo run --bin gentle
cargo run --bin gentle_js
cargo run --bin gentle_lua
cargo run --bin gentle_cli -- capabilities
cargo run --bin gentle -- --version
cargo run --bin gentle_cli -- --version
```

For optimized builds:

```bash
cargo run --release --bin gentle
cargo run --release --bin gentle_js
cargo run --release --bin gentle_lua
cargo run --release --bin gentle_cli -- capabilities
cargo run --release --bin gentle -- --version
cargo run --release --bin gentle_cli -- --version
```

Important: with `cargo run`, arguments for GENtle binaries must come after `--`.
Example: `cargo run --bin gentle -- --version`.

## `gentle` (GUI launcher)

`gentle` starts the graphical application.

```bash
cargo run --bin gentle
cargo run --bin gentle -- path/to/project.gentle.json
cargo run --bin gentle -- --version
```

If a project path is provided, GENtle opens that project at startup.

Current behavior:

- opens GUI windows
- opens an empty project by default
- optionally opens a project path passed on startup
- supports project/sequence actions from the File menu

## `gentle_js` (JavaScript shell)

`gentle_js` starts an interactive JavaScript REPL backed by GENtle data structures.

```bash
cargo run --bin gentle_js
```

Exit methods:

- type `exit`
- `Ctrl+C`
- `Ctrl+D`

### JavaScript shell functions

1. `load_dna(path)`
   - Loads a DNA sequence from file.
2. `write_gb(seq, path)`
   - Writes a sequence object to a GenBank file.
3. `load_project(path)` / `save_project(state, path)`
   - Loads/saves GENtle project JSON.
4. `capabilities()`
   - Returns shared-engine capabilities.
5. `apply_operation(state, op)`
   - Applies one engine operation to a project state.
   - `op` may be a JS object or JSON string.
   - Returns `{ state, result }`.
6. `apply_workflow(state, workflow)`
   - Applies a workflow to a project state.
   - `workflow` may be a JS object or JSON string.
   - Returns `{ state, results }`.

### JavaScript example

```javascript
state = { sequences: {}, metadata: {}, display: {}, lineage: {} };
op = { LoadFile: { path: "test_files/pGEX-3X.gb", as_id: "pgex" } };
r1 = apply_operation(state, op);
r2 = apply_operation(r1.state, {
  Digest: { input: "pgex", enzymes: ["BamHI", "EcoRI"], output_prefix: "frag" },
});
console.log(r2.result.created_seq_ids);
save_project(r2.state, "session.gentle.json");
```

## `gentle_lua` (Lua shell)

`gentle_lua` starts an interactive Lua REPL backed by GENtle data structures.

```bash
cargo run --bin gentle_lua
```

Exit methods:

- type `exit`
- `Ctrl+C`
- `Ctrl+D`

### Lua shell functions

1. `load_dna(filename)`
   - Loads a DNA sequence from file.
2. `write_gb(filename, seq)`
   - Writes a sequence object to a GenBank file.
3. `load_project(filename)` / `save_project(filename, project)`
   - Loads/saves GENtle project JSON.
4. `capabilities()`
   - Returns shared-engine capabilities.
5. `apply_operation(project, op)`
   - Applies one engine operation; `op` can be Lua table or JSON string.
   - Returns table with `state` and `result`.
6. `apply_workflow(project, workflow)`
   - Applies workflow; `workflow` can be Lua table or JSON string.
   - Returns table with `state` and `results`.

### Lua example

```lua
state = { sequences = {}, metadata = {}, display = {}, lineage = {} }
r1 = apply_operation(state, {
  LoadFile = { path = "test_files/pGEX-3X.gb", as_id = "pgex" }
})
r2 = apply_operation(r1.state, {
  Branch = { input = "pgex", output_id = "pgex_copy" }
})
print(r2.result.created_seq_ids[1])
save_project("session.gentle.json", r2.state)
```

## File format expectations

Current CLI workflows rely on sequence files supported by internal loaders:

- GenBank
- FASTA

## `gentle_cli` (automation/agent CLI)

`gentle_cli` exposes operation-based control with JSON input and output.
It is intended for scripted workflows and AI agents.

State is persisted to `.gentle_state.json` by default (override with `--state PATH`).

### Commands

```bash
cargo run --bin gentle_cli -- --version
cargo run --bin gentle_cli -- capabilities
cargo run --bin gentle_cli -- state-summary
cargo run --bin gentle_cli -- op '<operation-json>'
cargo run --bin gentle_cli -- workflow '<workflow-json>'
cargo run --bin gentle_cli -- export-state state.json
cargo run --bin gentle_cli -- import-state state.json
cargo run --bin gentle_cli -- save-project project.gentle.json
cargo run --bin gentle_cli -- load-project project.gentle.json
cargo run --bin gentle_cli -- render-svg pgex linear pgex.linear.svg
cargo run --bin gentle_cli -- render-svg pgex circular pgex.circular.svg
```

You can pass JSON from a file with `@file.json`.

Project aliases:

- `save-project PATH` aliases `export-state PATH`
- `load-project PATH` aliases `import-state PATH`

### Example operations

Load a file:

```json
{"LoadFile":{"path":"test_files/pGEX-3X.gb","as_id":"pgex"}}
```

Digest:

```json
{"Digest":{"input":"pgex","enzymes":["BamHI","EcoRI"],"output_prefix":"pgex_frag"}}
```

Ligation:

```json
{"Ligation":{"inputs":["pgex_frag_1","pgex_frag_2"],"circularize_if_possible":true,"output_id":"re_ligated","protocol":"Sticky","output_prefix":"lig","unique":true}}
```

Merge multiple containers/pools into one candidate pool namespace:

```json
{"MergeContainers":{"inputs":["tubeA_1","tubeB_3","tubeC_2"],"output_prefix":"merged_pool"}}
```

Note:

- Ligation is protocol-driven and requires explicit `protocol`
  (`"Sticky"` or `"Blunt"`).

PCR:

```json
{"Pcr":{"template":"pgex","forward_primer":"ATGGCT","reverse_primer":"CGTACC","output_id":"amplicon1","unique":true}}
```

Advanced PCR (tail insertion / mismatch-tolerant annealing):

```json
{"PcrAdvanced":{"template":"pgex","forward_primer":{"sequence":"GGATCCATGGCT","anneal_len":6,"max_mismatches":1,"require_3prime_exact_bases":4},"reverse_primer":{"sequence":"CGTACC","anneal_len":6,"max_mismatches":0,"require_3prime_exact_bases":4},"output_id":"amplicon_site_added","unique":true}}
```

Advanced PCR with degenerate/randomized primer library:

```json
{"PcrAdvanced":{"template":"pgex","forward_primer":{"sequence":"ATNAAA","anneal_len":6,"max_mismatches":1,"require_3prime_exact_bases":3,"library_mode":"Sample","max_variants":10,"sample_seed":42},"reverse_primer":{"sequence":"AAACCC","anneal_len":6,"max_mismatches":0,"require_3prime_exact_bases":4},"unique":false}}
```

PCR mutagenesis (explicit SNP intent):

```json
{"PcrMutagenesis":{"template":"pgex","forward_primer":{"sequence":"ATCAAA","anneal_len":6,"max_mismatches":1,"require_3prime_exact_bases":3},"reverse_primer":{"sequence":"AAACCC","anneal_len":6,"max_mismatches":0,"require_3prime_exact_bases":4},"mutations":[{"zero_based_position":2,"reference":"G","alternate":"C"}],"output_id":"snp_product","unique":true,"require_all_mutations":true}}
```

Extract region (`from` inclusive, `to` exclusive):

```json
{"ExtractRegion":{"input":"pgex","from":100,"to":900,"output_id":"insert_candidate"}}
```

Select one candidate in-silico (explicit provenance step):

```json
{"SelectCandidate":{"input":"pgex_frag_1","criterion":"band_size_range:450-550bp","output_id":"picked_band"}}
```

Filter candidates by molecular-weight proxy (bp length) with tolerance and uniqueness:

```json
{"FilterByMolecularWeight":{"inputs":["pgex_frag_1","pgex_frag_2","pgex_frag_3"],"min_bp":450,"max_bp":550,"error":0.10,"unique":true,"output_prefix":"mw_pick"}}
```

Create transformed or branched candidates:

```json
{"Reverse":{"input":"pgex","output_id":"pgex_rev"}}
{"Complement":{"input":"pgex","output_id":"pgex_comp"}}
{"ReverseComplement":{"input":"pgex","output_id":"pgex_revcomp"}}
{"Branch":{"input":"pgex","output_id":"pgex_branch"}}
```

Set visibility of a GUI-equivalent display toggle (example: features):

```json
{"SetDisplayVisibility":{"target":"Features","visible":false}}
```

Set an in-silico engine parameter (example: cap fragment/product combinatorics):

```json
{"SetParameter":{"name":"max_fragments_per_container","value":80000}}
```

Available `target` values:

- `SequencePanel`
- `MapPanel`
- `Features`
- `RestrictionEnzymes`
- `GcContents`
- `OpenReadingFrames`
- `MethylationSites`

Save as GenBank:

```json
{"SaveFile":{"seq_id":"pgex_frag_1","path":"frag1.gb","format":"GenBank"}}
```

### Current limitations in the new operation layer

- Ligation currently concatenates input sequences in order without sticky-end compatibility checks.
- PCR currently supports linear templates and exact primer matching only (no mismatch model yet).

## Error behavior

- REPLs print runtime errors without exiting automatically.
- Invalid file paths or unsupported content produce load/write errors.
- Some functions assume computed annotations are available after `load_dna`.
