# GENtle CLI Manual

This page documents command-line entry points for GENtle.

## Overview

GENtle currently provides three binaries:

- `gentle`: graphical desktop app
- `gentle_js`: interactive JavaScript shell
- `gentle_lua`: interactive Lua shell
- `gentle_cli`: JSON operation/workflow CLI for automation and AI tools

Resource update capability status:

- `gentle_cli`: supported (`resources sync-rebase`, `resources sync-jaspar`)
- `gentle_js`: supported (`sync_rebase`, `sync_jaspar`)
- `gentle_lua`: supported (`sync_rebase`, `sync_jaspar`)

Reference genome capability status:

- `gentle_cli`: supported via shared engine operations (`PrepareGenome`, `ExtractGenomeRegion`, `ExtractGenomeGene`) and shell-level `genomes/helpers blast`
- `gentle_js`: supported via dedicated helpers (`list_reference_genomes`, `is_reference_genome_prepared`, `list_reference_genome_genes`, `blast_reference_genome`, `blast_helper_genome`, `prepare_genome`, `extract_genome_region`, `extract_genome_gene`) and `apply_operation`
- `gentle_lua`: supported via dedicated helpers (`list_reference_genomes`, `is_reference_genome_prepared`, `list_reference_genome_genes`, `blast_reference_genome`, `blast_helper_genome`, `prepare_genome`, `extract_genome_region`, `extract_genome_gene`) and `apply_operation`

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
5. `state_summary(state)`
   - Returns normalized sequence/container/display summary.
6. `inspect_dna_ladders(name_filter)`
   - Returns built-in ladder catalog as structured JSON.
   - `name_filter` is optional (`null`/`""` means all ladders).
7. `export_dna_ladders(path, name_filter)`
   - Writes ladder catalog JSON to `path`.
   - Optional `name_filter` limits exported ladders by case-insensitive name match.
8. `apply_operation(state, op)`
   - Applies one engine operation to a project state.
   - `op` may be a JS object or JSON string.
   - Returns `{ state, result }`.
9. `apply_workflow(state, workflow)`
   - Applies a workflow to a project state.
   - `workflow` may be a JS object or JSON string.
   - Returns `{ state, results }`.
10. `sync_rebase(input, output, commercial_only)`
   - Parses REBASE/Bairoch input and writes a REBASE resource JSON snapshot.
   - `output` is optional (`null`/`""` uses default runtime resource path).
11. `sync_jaspar(input, output)`
   - Parses JASPAR PFM text and writes motif resource JSON snapshot.
   - `output` is optional (`null`/`""` uses default runtime resource path).
12. `list_reference_genomes(catalog_path)`
    - Lists genome IDs from the genome catalog.
    - `catalog_path` is optional (`null`/`""` uses default catalog).
13. `is_reference_genome_prepared(genome_id, catalog_path, cache_dir)`
    - Checks whether a genome cache is prepared and indexed.
14. `list_reference_genome_genes(genome_id, catalog_path, cache_dir)`
    - Lists indexed genes from prepared genome cache.
15. `prepare_genome(state, genome_id, catalog_path, cache_dir)`
    - Convenience wrapper around engine `PrepareGenome`.
16. `extract_genome_region(state, genome_id, chromosome, start_1based, end_1based, output_id, catalog_path, cache_dir)`
    - Convenience wrapper around engine `ExtractGenomeRegion`.
17. `extract_genome_gene(state, genome_id, gene_query, occurrence, output_id, catalog_path, cache_dir)`
    - Convenience wrapper around engine `ExtractGenomeGene`.
18. `blast_reference_genome(genome_id, query_sequence, max_hits, task, catalog_path, cache_dir)`
    - Runs BLAST (`blastn`/`blastn-short`) against a prepared reference genome.
19. `blast_helper_genome(helper_id, query_sequence, max_hits, task, catalog_path, cache_dir)`
    - Same as `blast_reference_genome`, but defaults to helper catalog context.

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
5. `state_summary(project)`
   - Returns normalized sequence/container/display summary.
6. `inspect_dna_ladders([name_filter])`
   - Returns built-in ladder catalog as structured table.
7. `export_dna_ladders(output_json, [name_filter])`
   - Writes ladder catalog JSON to disk.
8. `apply_operation(project, op)`
   - Applies one engine operation; `op` can be Lua table or JSON string.
   - Returns table with `state` and `result`.
9. `apply_workflow(project, workflow)`
   - Applies workflow; `workflow` can be Lua table or JSON string.
   - Returns table with `state` and `results`.
10. `sync_rebase(input, output, commercial_only)`
   - Parses REBASE/Bairoch input and writes a REBASE resource JSON snapshot.
   - `output` and `commercial_only` are optional.
11. `sync_jaspar(input, output)`
   - Parses JASPAR PFM text and writes motif resource JSON snapshot.
   - `output` is optional.
12. `list_reference_genomes(catalog_path)`
    - Lists genome IDs from the genome catalog.
13. `is_reference_genome_prepared(genome_id, catalog_path, cache_dir)`
    - Checks whether a genome cache is prepared and indexed.
14. `list_reference_genome_genes(genome_id, catalog_path, cache_dir)`
    - Lists indexed genes from prepared genome cache.
15. `prepare_genome(project, genome_id, catalog_path, cache_dir)`
    - Convenience wrapper around engine `PrepareGenome`.
16. `extract_genome_region(project, genome_id, chromosome, start_1based, end_1based, output_id, catalog_path, cache_dir)`
    - Convenience wrapper around engine `ExtractGenomeRegion`.
17. `extract_genome_gene(project, genome_id, gene_query, occurrence, output_id, catalog_path, cache_dir)`
    - Convenience wrapper around engine `ExtractGenomeGene`.
18. `blast_reference_genome(genome_id, query_sequence, [max_hits], [task], [catalog_path], [cache_dir])`
    - Runs BLAST (`blastn`/`blastn-short`) against a prepared reference genome.
19. `blast_helper_genome(helper_id, query_sequence, [max_hits], [task], [catalog_path], [cache_dir])`
    - Same as `blast_reference_genome`, but defaults to helper catalog context.

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
  - default interpretation: synthetic blunt `dsDNA`
  - optional FASTA-header metadata tokens for synthetic oligos:
    - `molecule=ssdna` for single-stranded DNA
    - `molecule=rna` for RNA (input `T` is normalized to `U`)
    - `molecule=dsdna` plus optional overhangs:
      - `forward_5=...` (alias `f5=...`)
      - `forward_3=...` (alias `f3=...`)
      - `reverse_5=...` (alias `r5=...`)
      - `reverse_3=...` (alias `r3=...`)

Example FASTA headers:

- `>oligo_ss molecule=ssdna`
- `>oligo_rna molecule=rna`
- `>oligo_ds molecule=dsdna f5=GATC r5=CTAG`

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
cargo run --bin gentle_cli -- --progress op '<operation-json>'
cargo run --bin gentle_cli -- --progress-stdout workflow '<workflow-json>'
cargo run --bin gentle_cli -- export-state state.json
cargo run --bin gentle_cli -- import-state state.json
cargo run --bin gentle_cli -- save-project project.gentle.json
cargo run --bin gentle_cli -- load-project project.gentle.json
cargo run --bin gentle_cli -- render-svg pgex linear pgex.linear.svg
cargo run --bin gentle_cli -- render-svg pgex circular pgex.circular.svg
cargo run --bin gentle_cli -- render-rna-svg rna_seq rna.secondary.svg
cargo run --bin gentle_cli -- rna-info rna_seq
cargo run --bin gentle_cli -- render-lineage-svg lineage.svg
cargo run --bin gentle_cli -- shell 'help'
cargo run --bin gentle_cli -- shell 'state-summary'
cargo run --bin gentle_cli -- shell 'op @op.json'
cargo run --bin gentle_cli -- render-pool-gel-svg frag_1,frag_2 digest.gel.svg
cargo run --bin gentle_cli -- render-pool-gel-svg frag_1,frag_2 digest.gel.svg --ladders "NEB 100bp DNA Ladder,NEB 1kb DNA Ladder"
cargo run --bin gentle_cli -- ladders list
cargo run --bin gentle_cli -- ladders list --filter NEB
cargo run --bin gentle_cli -- ladders list --molecule rna
cargo run --bin gentle_cli -- ladders export dna_ladders.snapshot.json
cargo run --bin gentle_cli -- ladders export dna_ladders.neb.json --filter NEB
cargo run --bin gentle_cli -- ladders export rna_ladders.snapshot.json --molecule rna
cargo run --bin gentle_cli -- export-pool frag_1,frag_2 digest.pool.gentle.json "BamHI+EcoRI digest pool"
cargo run --bin gentle_cli -- import-pool digest.pool.gentle.json imported
cargo run --bin gentle_cli -- resources sync-rebase rebase.withrefm data/resources/rebase.enzymes.json --commercial-only
cargo run --bin gentle_cli -- resources sync-jaspar JASPAR2026_CORE_non-redundant_pfms_jaspar.txt data/resources/jaspar.motifs.json
cargo run --bin gentle_cli -- op '{"PrepareGenome":{"genome_id":"ToyGenome","catalog_path":"catalog.json"}}'
cargo run --bin gentle_cli -- op '{"ExtractGenomeRegion":{"genome_id":"ToyGenome","chromosome":"chr1","start_1based":1001,"end_1based":1600,"output_id":"toy_chr1_1001_1600","catalog_path":"catalog.json"}}'
cargo run --bin gentle_cli -- op '{"ExtractGenomeGene":{"genome_id":"ToyGenome","gene_query":"MYGENE","occurrence":1,"output_id":"toy_mygene","catalog_path":"catalog.json"}}'
cargo run --bin gentle_cli -- genomes list --catalog assets/genomes.json
cargo run --bin gentle_cli -- genomes list --catalog assets/helper_genomes.json
cargo run --bin gentle_cli -- genomes validate-catalog --catalog assets/genomes.json
cargo run --bin gentle_cli -- genomes status "Human GRCh38 Ensembl 116" --catalog assets/genomes.json --cache-dir data/genomes
cargo run --bin gentle_cli -- genomes genes "Human GRCh38 Ensembl 116" --catalog assets/genomes.json --cache-dir data/genomes --filter "^TP53$" --biotype protein_coding --limit 20
cargo run --bin gentle_cli -- genomes prepare "Human GRCh38 Ensembl 116" --catalog assets/genomes.json --cache-dir data/genomes
cargo run --bin gentle_cli -- genomes blast "Human GRCh38 Ensembl 116" ACGTACGTACGT --task blastn-short --max-hits 10 --catalog assets/genomes.json --cache-dir data/genomes
cargo run --bin gentle_cli -- genomes extract-region "Human GRCh38 Ensembl 116" 1 1000000 1001500 --output-id grch38_chr1_slice --catalog assets/genomes.json --cache-dir data/genomes
cargo run --bin gentle_cli -- genomes extract-gene "Human GRCh38 Ensembl 116" TP53 --occurrence 1 --output-id grch38_tp53 --catalog assets/genomes.json --cache-dir data/genomes
cargo run --bin gentle_cli -- tracks import-bed grch38_tp53 data/chipseq/peaks.bed.gz --name H3K27ac --min-score 10 --clear-existing
cargo run --bin gentle_cli -- tracks import-bigwig grch38_tp53 data/chipseq/signal.bw --name ATAC --min-score 0.2 --clear-existing
cargo run --bin gentle_cli -- helpers list
cargo run --bin gentle_cli -- helpers validate-catalog
cargo run --bin gentle_cli -- helpers status "Plasmid pUC19 (local)"
cargo run --bin gentle_cli -- helpers prepare "Plasmid pUC19 (local)" --cache-dir data/helper_genomes
cargo run --bin gentle_cli -- helpers genes "Plasmid pUC19 (local)" --filter bla --limit 20
cargo run --bin gentle_cli -- helpers blast "Plasmid pUC19 (local)" ACGTACGTACGT --task blastn-short --max-hits 10 --cache-dir data/helper_genomes
```

You can pass JSON from a file with `@file.json`.

Global CLI options:

- `--state PATH`: use a non-default project state file
- `--progress` or `--progress-stderr`: print live progress events to `stderr`
- `--progress-stdout`: print live progress events to `stdout`
- `--allow-screenshots`: opt-in guard for window screenshot commands

Current progress events include TFBS annotation updates and genome-prepare
updates (download/index phases).
When `--progress-stdout` is used, progress lines are emitted before the final JSON output.

`state-summary` output includes:

- sequences
- containers
- display visibility flags

Project aliases:

- `save-project PATH` aliases `export-state PATH`
- `load-project PATH` aliases `import-state PATH`

Shared shell command:

- `shell '<command line>'`
  - Parses and executes commands through the same shared parser/executor used by
    the GUI `Shell` panel (`src/engine_shell.rs`).
  - Shell-supported commands:
    - `help`
    - `capabilities`
    - `state-summary`
    - `load-project PATH`
    - `save-project PATH`
    - `render-svg SEQ_ID linear|circular OUTPUT.svg`
    - `render-rna-svg SEQ_ID OUTPUT.svg`
    - `rna-info SEQ_ID`
    - `render-lineage-svg OUTPUT.svg`
    - `render-pool-gel-svg IDS OUTPUT.svg [--ladders NAME[,NAME]]`
    - `ladders list [--molecule dna|rna] [--filter TEXT]`
    - `ladders export OUTPUT.json [--molecule dna|rna] [--filter TEXT]`
    - `export-pool IDS OUTPUT.pool.gentle.json [HUMAN_ID]`
    - `import-pool INPUT.pool.gentle.json [PREFIX]`
    - `resources sync-rebase INPUT.withrefm_or_URL [OUTPUT.rebase.json] [--commercial-only]`
    - `resources sync-jaspar INPUT.jaspar_or_URL [OUTPUT.motifs.json]`
    - `genomes list [--catalog PATH]`
    - `genomes validate-catalog [--catalog PATH]`
    - `genomes status GENOME_ID [--catalog PATH] [--cache-dir PATH]`
    - `genomes genes GENOME_ID [--catalog PATH] [--cache-dir PATH] [--filter REGEX] [--biotype NAME] [--limit N] [--offset N]`
    - `genomes prepare GENOME_ID [--catalog PATH] [--cache-dir PATH]`
    - `genomes blast GENOME_ID QUERY_SEQUENCE [--max-hits N] [--task blastn-short|blastn] [--catalog PATH] [--cache-dir PATH]`
    - `genomes extract-region GENOME_ID CHR START END [--output-id ID] [--catalog PATH] [--cache-dir PATH]`
    - `genomes extract-gene GENOME_ID QUERY [--occurrence N] [--output-id ID] [--catalog PATH] [--cache-dir PATH]`
    - `helpers list [--catalog PATH]`
    - `helpers validate-catalog [--catalog PATH]`
    - `helpers status HELPER_ID [--catalog PATH] [--cache-dir PATH]`
    - `helpers genes HELPER_ID [--catalog PATH] [--cache-dir PATH] [--filter REGEX] [--biotype NAME] [--limit N] [--offset N]`
    - `helpers prepare HELPER_ID [--catalog PATH] [--cache-dir PATH]`
    - `helpers blast HELPER_ID QUERY_SEQUENCE [--max-hits N] [--task blastn-short|blastn] [--catalog PATH] [--cache-dir PATH]`
    - `helpers extract-region HELPER_ID CHR START END [--output-id ID] [--catalog PATH] [--cache-dir PATH]`
    - `helpers extract-gene HELPER_ID QUERY [--occurrence N] [--output-id ID] [--catalog PATH] [--cache-dir PATH]`
    - `tracks import-bed SEQ_ID PATH [--name NAME] [--min-score N] [--max-score N] [--clear-existing]`
    - `tracks import-bigwig SEQ_ID PATH [--name NAME] [--min-score N] [--max-score N] [--clear-existing]`
    - `op <operation-json-or-@file>`
    - `workflow <workflow-json-or-@file>`
    - `screenshot-window OUTPUT.png` (requires process startup with
      `--allow-screenshots`)
  - Use single quotes around JSON payloads to preserve whitespace:
    - `gentle_cli shell 'workflow {"run_id":"r1","ops":[]}'`

Screenshot bridge:

- Startup guard:
  - screenshot commands are accepted only when `gentle_cli` was started with
    `--allow-screenshots`.
- Command surface:
  - direct CLI: `gentle_cli --allow-screenshots screenshot-window OUTPUT.png`
  - shared shell: `gentle_cli --allow-screenshots shell 'screenshot-window OUTPUT.png'`
- Scope/safety:
  - captures only the active/topmost GENtle window
  - full-screen capture and non-GENtle window capture are out of contract
  - window lookup is native AppKit in-process (no AppleScript automation path)
  - in practice this is most useful from the GUI shell; a headless `gentle_cli`
    process usually has no eligible window
  - current backend support is macOS (`screencapture`); non-macOS returns an
    unsupported error for now
- Output:
  - caller chooses the full output path/filename (for docs auto-update workflows)
  - expected result schema: `gentle.screenshot.v1`

Pool exchange commands:

- `export-pool IDS OUTPUT.pool.gentle.json [HUMAN_ID]`
  - Exports explicit sequence IDs (`IDS` is comma-separated) with topology and
    overhang fields into a versioned JSON pool artifact.
  - Adds `human_id` at pool level and per-member `human_id`.
- `import-pool INPUT.pool.gentle.json [PREFIX]`
  - Imports pool members into current state; generated IDs are prefixed.

Rendering export commands:

- `render-svg SEQ_ID linear|circular OUTPUT.svg`
  - Calls engine operation `RenderSequenceSvg`.
- `render-rna-svg SEQ_ID OUTPUT.svg`
  - Calls engine operation `RenderRnaStructureSvg`.
  - Accepts only single-stranded RNA (`molecule_type` of `RNA`/`ssRNA`).
  - Uses external `rnapkin` executable (set `GENTLE_RNAPKIN_BIN` to override executable path).
- `render-lineage-svg OUTPUT.svg`
  - Calls engine operation `RenderLineageSvg`.
- `render-pool-gel-svg IDS OUTPUT.svg [--ladders NAME[,NAME]]`
  - Calls engine operation `RenderPoolGelSvg`.
  - `IDS` is a comma-separated sequence-id list.
  - Optional `--ladders` overrides auto ladder selection.
  - If `--ladders` is omitted, engine auto-selects one or two ladders based on
    pool bp range.

RNA secondary-structure text command:

- `rna-info SEQ_ID`
  - Returns a JSON report from `rnapkin -v -p <sequence>`.
  - Includes `stdout`/`stderr` and command metadata.
  - Accepts only single-stranded RNA (`molecule_type` of `RNA`/`ssRNA`).
  - Uses external `rnapkin` executable (set `GENTLE_RNAPKIN_BIN` to override executable path).

DNA ladder catalog commands:

- `ladders list [--molecule dna|rna] [--filter TEXT]`
  - Returns engine-inspected ladder catalog (`schema`, `ladder_count`, `ladders`).
  - Default molecule is `dna`; use `--molecule rna` for RNA ladders.
  - Optional `--filter` keeps only ladders whose names contain `TEXT` (case-insensitive).
- `ladders export OUTPUT.json [--molecule dna|rna] [--filter TEXT]`
  - Calls engine operation `ExportDnaLadders` or `ExportRnaLadders`.
  - Writes structured ladder catalog JSON to disk.

Resource sync commands:

- `resources sync-rebase INPUT.withrefm [OUTPUT.rebase.json] [--commercial-only]`
  - Parses REBASE/Bairoch-style records (`withrefm`) into GENtle restriction-enzyme JSON.
  - `INPUT` may be a local file path or an `https://...` URL.
  - Default output: `data/resources/rebase.enzymes.json`.
  - At runtime, GENtle auto-loads this file (if present) and overrides embedded restriction enzymes.
- `resources sync-jaspar INPUT.jaspar.txt [OUTPUT.motifs.json]`
  - Parses JASPAR PFM text into motif snapshot JSON with IUPAC consensus.
  - `INPUT` may be a local file path or an `https://...` URL.
  - Default output: `data/resources/jaspar.motifs.json`.
  - `ExtractAnchoredRegion` can resolve TF motifs by ID/name from this local registry.
  - GENtle ships with built-in motifs in `assets/jaspar.motifs.json` (currently generated from JASPAR 2026 CORE non-redundant); this command provides local updates/extensions.

Genome convenience commands:

- `genomes list [--catalog PATH]`
  - Lists available genomes in the catalog.
- `genomes validate-catalog [--catalog PATH]`
  - Verifies catalog JSON schema/entry rules and that each entry resolves usable
    sequence/annotation source definitions.
- `genomes status GENOME_ID [--catalog PATH] [--cache-dir PATH]`
  - Shows whether the genome cache is prepared/indexed.
  - Also reports resolved source details: `sequence_source_type`,
    `annotation_source_type`, `sequence_source`, `annotation_source`.
- `genomes genes GENOME_ID [--catalog PATH] [--cache-dir PATH] [--filter REGEX] [--biotype NAME] [--limit N] [--offset N]`
  - Lists indexed genes from prepared cache (paged by `--limit`/`--offset`).
  - `--filter` is a case-insensitive regular expression.
  - `--biotype` can be repeated to constrain matches to selected biotypes.
- `genomes prepare GENOME_ID [--catalog PATH] [--cache-dir PATH]`
  - Runs engine `PrepareGenome`.
- `genomes blast GENOME_ID QUERY_SEQUENCE [--max-hits N] [--task blastn-short|blastn] [--catalog PATH] [--cache-dir PATH]`
  - Runs `blastn` against prepared genome cache/index.
  - `--task` defaults to `blastn-short`; accepted values: `blastn-short`, `blastn`.
- `genomes extract-region GENOME_ID CHR START END [--output-id ID] [--catalog PATH] [--cache-dir PATH]`
  - Runs engine `ExtractGenomeRegion`.
- `genomes extract-gene GENOME_ID QUERY [--occurrence N] [--output-id ID] [--catalog PATH] [--cache-dir PATH]`
  - Runs engine `ExtractGenomeGene`.

Helper convenience commands:

- `helpers list [--catalog PATH]`
  - Same behavior as `genomes list`, but defaults to `assets/helper_genomes.json`.
- `helpers validate-catalog [--catalog PATH]`
  - Same behavior as `genomes validate-catalog`, with helper-catalog default.
- `helpers status HELPER_ID [--catalog PATH] [--cache-dir PATH]`
  - Same behavior as `genomes status`, with helper-catalog default.
- `helpers genes HELPER_ID [--catalog PATH] [--cache-dir PATH] [--filter REGEX] [--biotype NAME] [--limit N] [--offset N]`
  - Same behavior as `genomes genes`, with helper-catalog default.
- `helpers prepare HELPER_ID [--catalog PATH] [--cache-dir PATH]`
  - Same behavior as `genomes prepare`, with helper-catalog default.
- `helpers blast HELPER_ID QUERY_SEQUENCE [--max-hits N] [--task blastn-short|blastn] [--catalog PATH] [--cache-dir PATH]`
  - Same behavior as `genomes blast`, with helper-catalog default.
- `helpers extract-region HELPER_ID CHR START END [--output-id ID] [--catalog PATH] [--cache-dir PATH]`
  - Same behavior as `genomes extract-region`, with helper-catalog default.
- `helpers extract-gene HELPER_ID QUERY [--occurrence N] [--output-id ID] [--catalog PATH] [--cache-dir PATH]`
  - Same behavior as `genomes extract-gene`, with helper-catalog default.

Recommended suffixes:

- pool artifacts: `*.pool.gentle.json`
- project/state artifacts: `*.project.gentle.json` (or existing `*.gentle.json`)

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

Extract anchored region with flexible 5' boundary and fixed anchor side:

```json
{"ExtractAnchoredRegion":{"input":"pgex","anchor":{"FeatureBoundary":{"feature_kind":"CDS","feature_label":null,"boundary":"Start","occurrence":0}},"direction":"Upstream","target_length_bp":500,"length_tolerance_bp":100,"required_re_sites":["EcoRI"],"required_tf_motifs":["TATAAA"],"forward_primer":"GAATTC","reverse_primer":"CGTACC","output_prefix":"promoter","unique":false,"max_candidates":20}}
```

`required_tf_motifs` accepts either IUPAC motif strings or motif IDs/names from
the local JASPAR snapshot.

Annotate TFBS features using log-likelihood ratio thresholds:

```json
{"AnnotateTfbs":{"seq_id":"pgex","motifs":["MA0139.1","SP1","TATAAA"],"min_llr_bits":0.0,"min_llr_quantile":0.95,"per_tf_thresholds":[{"tf":"SP1","min_llr_bits":-1.0,"min_llr_quantile":0.80}],"clear_existing":true,"max_hits":500}}
```

`AnnotateTfbs` writes generated TFBS features with score qualifiers:

- `llr_bits`: absolute log-likelihood ratio score (base 2)
- `llr_quantile`: empirical score quantile in the scanned sequence region (both strands)
- `true_log_odds_bits`: smoothed true log-odds score (base 2)
- `true_log_odds_quantile`: empirical quantile of `true_log_odds_bits` in the scanned region
- compatibility aliases are also written: `log_odds_ratio_bits`, `log_odds_ratio_quantile`

Motif-selection shortcuts:

- `motifs: ["ALL"]` or `motifs: ["*"]` scans all motifs from the local JASPAR registry.
- `max_hits` controls safety capping of generated TFBS features:
  - omitted: default cap of `500` hits
  - `0`: unlimited (no cap)
  - `N > 0`: stop after `N` accepted hits

Select one candidate in-silico (explicit provenance step):

```json
{"SelectCandidate":{"input":"pgex_frag_1","criterion":"band_size_range:450-550bp","output_id":"picked_band"}}
```

Filter candidates by molecular-weight proxy (bp length) with tolerance and uniqueness:

```json
{"FilterByMolecularWeight":{"inputs":["pgex_frag_1","pgex_frag_2","pgex_frag_3"],"min_bp":450,"max_bp":550,"error":0.10,"unique":true,"output_prefix":"mw_pick"}}
```

Apply practical sequence-quality filters (GC bounds, homopolymer cap, U6 `TTTT`
avoidance, forbidden motifs):

```json
{"FilterBySequenceQuality":{"inputs":["g1","g2","g3"],"gc_min":0.30,"gc_max":0.70,"max_homopolymer_run":4,"reject_ambiguous_bases":true,"avoid_u6_terminator_tttt":true,"forbidden_motifs":["GAATTC"],"unique":false,"output_prefix":"sq_pick"}}
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

Set linear-map viewport (zoom/pan state):

```json
{"SetLinearViewport":{"start_bp":1000,"span_bp":5000}}
```

Set an in-silico engine parameter (example: cap fragment/product combinatorics):

```json
{"SetParameter":{"name":"max_fragments_per_container","value":80000}}
```

Set feature-details font size used in the feature tree/details panel (valid range `8.0..24.0`):

```json
{"SetParameter":{"name":"feature_details_font_size","value":10.5}}
```

Export selected pool members (engine operation):

```json
{"ExportPool":{"inputs":["frag_1","frag_2"],"path":"digest.pool.gentle.json","pool_id":"digest_1","human_id":"BamHI+EcoRI digest"}}
```

Export built-in DNA ladder catalog (engine operation):

```json
{"ExportDnaLadders":{"path":"dna_ladders.snapshot.json","name_filter":null}}
{"ExportDnaLadders":{"path":"dna_ladders.neb.json","name_filter":"NEB"}}
{"ExportRnaLadders":{"path":"rna_ladders.snapshot.json","name_filter":null}}
{"ExportRnaLadders":{"path":"rna_ladders.neb.json","name_filter":"NEB"}}
```

Render sequence SVG (engine operation):

```json
{"RenderSequenceSvg":{"seq_id":"pgex","mode":"Linear","path":"pgex.linear.svg"}}
{"RenderSequenceSvg":{"seq_id":"pgex","mode":"Circular","path":"pgex.circular.svg"}}
```

Render RNA secondary-structure SVG (engine operation):

```json
{"RenderRnaStructureSvg":{"seq_id":"rna_seq","path":"rna.secondary.svg"}}
```

Render lineage SVG (engine operation):

```json
{"RenderLineageSvg":{"path":"lineage.svg"}}
```

Render pool gel SVG with ladder selection (engine operation):

```json
{"RenderPoolGelSvg":{"inputs":["frag_1","frag_2","frag_3"],"path":"digest.gel.svg","ladders":["NEB 100bp DNA Ladder","NEB 1kb DNA Ladder"]}}
```

Render pool gel SVG with automatic ladder selection:

```json
{"RenderPoolGelSvg":{"inputs":["frag_1","frag_2","frag_3"],"path":"digest.auto.gel.svg","ladders":null}}
```

Prepare a whole reference genome once (download/copy sequence + annotation and
build local FASTA and BLAST indexes):

```json
{"PrepareGenome":{"genome_id":"Human GRCh38 Ensembl 116","catalog_path":"assets/genomes.json","cache_dir":"data/genomes"}}
```

Extract a genomic interval from the prepared cache (1-based inclusive
coordinates):

```json
{"ExtractGenomeRegion":{"genome_id":"Human GRCh38 Ensembl 116","chromosome":"1","start_1based":1000000,"end_1based":1001500,"output_id":"grch38_chr1_1000000_1001500","catalog_path":"assets/genomes.json","cache_dir":"data/genomes"}}
```

Extract by gene query (name/id) from the prepared gene index:

```json
{"ExtractGenomeGene":{"genome_id":"Human GRCh38 Ensembl 116","gene_query":"TP53","occurrence":1,"output_id":"grch38_tp53","catalog_path":"assets/genomes.json","cache_dir":"data/genomes"}}
```

Import BED intervals (plain `.bed` or gzipped `.bed.gz`) onto a genome-anchored
sequence:

```json
{"ImportGenomeBedTrack":{"seq_id":"grch38_tp53","path":"data/chipseq/peaks.bed.gz","track_name":"H3K27ac","min_score":10.0,"max_score":null,"clear_existing":false}}
```

Import BigWig signal tracks (`.bw` / `.bigWig`) onto a genome-anchored
sequence:

```json
{"ImportGenomeBigWigTrack":{"seq_id":"grch38_tp53","path":"data/chipseq/signal.bw","track_name":"ATAC","min_score":0.2,"max_score":null,"clear_existing":false}}
```

Notes:

- `PrepareGenome` is intended as a one-time setup step per genome and cache
  location.
- During prepare, GENtle also attempts to create a BLAST nucleotide index
  (`makeblastdb`) under the genome install directory.
- HTTP-based source downloads are resumable (Range requests with retry/backoff),
  and completed installs persist SHA-1 checksums for sequence/annotation files.
- A catalog entry can either define explicit URLs (`sequence_remote` /
  `annotations_remote`) or define `ncbi_assembly_accession` +
  `ncbi_assembly_name` to derive direct NCBI GenBank/RefSeq FTP downloads.
- Helper/vector entries can also define `genbank_accession`; when explicit
  URLs are absent, GENtle derives direct NCBI EFetch sources
  (`rettype=fasta` + `rettype=gbwithparts`) for one-time prepare/index.
- `ExtractGenomeRegion` expects the genome to have been prepared already.
- `ExtractGenomeGene` also expects prepared cache and gene index.
- `genomes/helpers blast` expects prepared cache and a BLAST index.
  If index files are missing, GENtle tries to build them on demand.
- BLAST executable overrides:
  - `GENTLE_MAKEBLASTDB_BIN` (default: `makeblastdb`)
  - `GENTLE_BLASTN_BIN` (default: `blastn`)
- `ImportGenomeBedTrack` expects `seq_id` to be a sequence created by
  `ExtractGenomeRegion` or `ExtractGenomeGene` (genome-anchored provenance).
- `ImportGenomeBigWigTrack` expects the same genome-anchored `seq_id`.
- BED import accepts local `.bed` and `.bed.gz` files.
- BigWig import accepts local `.bw` and `.bigWig` files and uses
  `bigWigToBedGraph` (override with `GENTLE_BIGWIG_TO_BEDGRAPH_BIN`).
- `ExtractGenomeRegion` and `ExtractGenomeGene` append extraction provenance
  records into `ProjectState.metadata["provenance"]["genome_extractions"]`
  (genome id, coordinates/query, source descriptors, and checksums when present).
- If `catalog_path` is omitted, engine default catalog is `assets/genomes.json`.
- Bundled `assets/genomes.json` currently includes Human GRCh38 (Ensembl 113 and 116),
  Mouse GRCm39 Ensembl 116, Rat GRCr8 Ensembl 116, Saccharomyces cerevisiae
  S288c (Ensembl 113 and 116), and `LocalProject` (backed by
  `test_files/AB011549.2.fa` + `test_files/AB011549.2.gb`).
- `cache_dir` is optional. If omitted, catalog/default cache settings are used.
- `chromosome` accepts exact contig names and also tolerates `chr` prefix
  differences (`1` vs `chr1`).
- For `ExtractGenomeGene`, `occurrence` is 1-based among matching records.

Available `target` values:

- `SequencePanel`
- `MapPanel`
- `Features`
- `CdsFeatures`
- `GeneFeatures`
- `MrnaFeatures`
- `Tfbs`
- `RestrictionEnzymes`
- `GcContents`
- `OpenReadingFrames`
- `MethylationSites`

Save as GenBank:

```json
{"SaveFile":{"seq_id":"pgex_frag_1","path":"frag1.gb","format":"GenBank"}}
```

### Current limitations in the new operation layer

- `import-pool` is currently a shared shell/CLI utility command and not yet an
  engine operation.

## Error behavior

- REPLs print runtime errors without exiting automatically.
- Invalid file paths or unsupported content produce load/write errors.
- Some functions assume computed annotations are available after `load_dna`.
