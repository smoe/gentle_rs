# GENtle CLI Manual

This page documents command-line entry points for GENtle.

## Overview

GENtle currently provides six binaries:

- `gentle`: graphical desktop app (GUI)
- `gentle_cli`: JSON operation/workflow CLI for automation and AI tools
- `gentle_js`: interactive JavaScript shell
- `gentle_lua`: interactive Lua shell
- `gentle_examples_docs`: generates adapter snippets from canonical protocol examples
- `gentle_mcp`: MCP stdio server (guarded mutating + UI-intent parity baseline)

In addition, the GUI includes an embedded `Shell` panel that uses the same
shared shell parser/executor as `gentle_cli shell`.

Structured command glossary:

- `docs/glossary.json` is the machine-readable command glossary used for
  per-command help rendering (`help ...`) and catalog export
  (`help --format json|markdown`).

Structured workflow examples:

- canonical source files: `docs/examples/workflows/*.json`
- schema: `gentle.workflow_example.v1`
- on-demand snippet generation (CLI/shell/JS/Lua):
  - `cargo run --bin gentle_examples_docs -- generate`
- validation only:
  - `cargo run --bin gentle_examples_docs -- --check`
- test gating:
  - `always` examples run in default tests
  - `online` examples run only with `GENTLE_TEST_ONLINE=1`
  - `skip` examples are syntax-checked only

Architecture invariant: all adapters/frontends above route cloning/business
behavior through the same shared engine.

Resource update capability status:

- `gentle_cli`: supported (`resources sync-rebase`, `resources sync-jaspar`)
- `gentle_js`: supported (`sync_rebase`, `sync_jaspar`)
- `gentle_lua`: supported (`sync_rebase`, `sync_jaspar`)

Reference genome capability status:

- `gentle_cli`: supported via shared engine operations (`PrepareGenome`, `ExtractGenomeRegion`, `ExtractGenomeGene`) and shell-level `genomes/helpers blast`
- `gentle_js`: supported via dedicated helpers (`list_reference_genomes`, `is_reference_genome_prepared`, `list_reference_genome_genes`, `blast_reference_genome`, `blast_helper_genome`, `prepare_genome`, `extract_genome_region`, `extract_genome_gene`) and `apply_operation`
- `gentle_lua`: supported via dedicated helpers (`list_reference_genomes`, `is_reference_genome_prepared`, `list_reference_genome_genes`, `blast_reference_genome`, `blast_helper_genome`, `prepare_genome`, `extract_genome_region`, `extract_genome_gene`) and `apply_operation`

Agent-assistant capability status:

- `gentle_cli`: supported via shared-shell command family (`agents list`, `agents ask`) and direct forwarding (`gentle_cli agents ...`)
- `gentle_js`: supported via helper wrappers (`list_agent_systems`, `ask_agent_system`) over shared shell execution
- `gentle_lua`: supported via helper wrappers (`list_agent_systems`, `ask_agent_system`) over shared shell execution
- GUI: supported via standalone `Agent Assistant` viewport using the same bridge and command executor

Candidate-set capability status:

- `gentle_cli`: supported as first-class `candidates ...` commands and shared-shell `candidates ...` commands, backed by shared engine operations (`GenerateCandidateSet`, `GenerateCandidateSetBetweenAnchors`, `DeleteCandidateSet`, `ScoreCandidateSetExpression`, `ScoreCandidateSetDistance`, `ScoreCandidateSetWeightedObjective`, `TopKCandidateSet`, `ParetoFrontierCandidateSet`, `FilterCandidateSet`, `CandidateSetOp`, `UpsertCandidateMacroTemplate`, `DeleteCandidateMacroTemplate`)
- `gentle_js`: supported via `apply_operation` with the same candidate-set operations
- `gentle_lua`: supported via `apply_operation` with the same candidate-set operations

Guide-design capability status:

- `gentle_cli`: supported as first-class `guides ...` commands and shared-shell `guides ...` commands, backed by shared engine operations (`UpsertGuideSet`, `DeleteGuideSet`, `FilterGuidesPractical`, `GenerateGuideOligos`, `ExportGuideOligos`, `ExportGuideProtocolText`)
- `gentle_js`: supported via `apply_operation` with the same guide-design operations
- `gentle_lua`: supported via `apply_operation` with the same guide-design operations

## Build and run

From the repository root:

```bash
cargo run --bin gentle
cargo run --bin gentle_js
cargo run --bin gentle_lua
cargo run --bin gentle_cli -- capabilities
cargo run --bin gentle_examples_docs -- --check
cargo run --bin gentle_mcp -- --help
cargo run --bin gentle -- --version
cargo run --bin gentle_cli -- --version
```

For optimized builds:

```bash
cargo run --release --bin gentle
cargo run --release --bin gentle_js
cargo run --release --bin gentle_lua
cargo run --release --bin gentle_cli -- capabilities
cargo run --release --bin gentle_examples_docs -- --check
cargo run --release --bin gentle_mcp -- --help
cargo run --release --bin gentle -- --version
cargo run --release --bin gentle_cli -- --version
```

Important: with `cargo run`, arguments for GENtle binaries must come after `--`.
Example: `cargo run --bin gentle -- --version`.

## Protocol-first example source

Canonical workflow examples are adapter-neutral JSON files:

- source: `docs/examples/workflows/*.json`
- generated adapter snippets: `docs/examples/generated/*.md`

Regenerate snippets on demand:

```bash
cargo run --bin gentle_examples_docs -- generate
```

Validate example files and schema without writing output:

```bash
cargo run --bin gentle_examples_docs -- --check
```

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
10. `import_pool(state, input_pool_json, prefix)`
   - Imports a `.pool.gentle.json` artifact into `state` via shared adapter logic.
   - `prefix` is optional (`null`/`""` defaults to `pool`).
   - Returns `{ state, state_changed, output }`.
11. `sync_rebase(input, output, commercial_only)`
   - Parses REBASE/Bairoch input and writes a REBASE resource JSON snapshot.
   - `output` is optional (`null`/`""` uses default runtime resource path).
12. `sync_jaspar(input, output)`
   - Parses JASPAR PFM text and writes motif resource JSON snapshot.
   - `output` is optional (`null`/`""` uses default runtime resource path).
13. `list_reference_genomes(catalog_path)`
    - Lists genome IDs from the genome catalog.
    - `catalog_path` is optional (`null`/`""` uses default catalog).
14. `is_reference_genome_prepared(genome_id, catalog_path, cache_dir)`
    - Checks whether a genome cache is prepared and indexed.
15. `list_reference_genome_genes(genome_id, catalog_path, cache_dir)`
    - Lists indexed genes from prepared genome cache.
16. `prepare_genome(state, genome_id, catalog_path, cache_dir)`
    - Convenience wrapper around engine `PrepareGenome`.
17. `extract_genome_region(state, genome_id, chromosome, start_1based, end_1based, output_id, catalog_path, cache_dir)`
    - Convenience wrapper around engine `ExtractGenomeRegion`.
18. `extract_genome_gene(state, genome_id, gene_query, occurrence, output_id, catalog_path, cache_dir)`
    - Convenience wrapper around engine `ExtractGenomeGene`.
19. `blast_reference_genome(genome_id, query_sequence, max_hits, task, catalog_path, cache_dir)`
    - Runs BLAST (`blastn`/`blastn-short`) against a prepared reference genome.
20. `blast_helper_genome(helper_id, query_sequence, max_hits, task, catalog_path, cache_dir)`
    - Same as `blast_reference_genome`, but defaults to helper catalog context.
21. `set_parameter(state, name, value)`
    - Convenience wrapper for engine `SetParameter`.
22. `set_vcf_display_filter(state, options)`
    - Convenience wrapper that updates one or more VCF display-filter parameters
      via `SetParameter` in one call.
23. `list_agent_systems(catalog_path)`
    - Lists configured agent systems from catalog JSON.
    - `catalog_path` is optional (`null`/`""` uses `assets/agent_systems.json`).
24. `ask_agent_system(state, system_id, prompt, options)`
    - Invokes one configured agent system through shared shell execution.
    - Returns `{ state, state_changed, output }`.
    - `state` may be `null` to use an empty/default project state.
    - `options` supports:
      - `catalog_path`
      - `allow_auto_exec` / `execute_all`
      - `execute_indices` (1-based)
      - `include_state_summary` (default `true`)

### JavaScript example

Use generated adapter snippets to stay synchronized with canonical workflow JSON:

- `docs/examples/generated/load_and_digest_pgex.md`
- `docs/examples/generated/load_branch_reverse_complement_pgex_fasta.md`
- `docs/examples/generated/guides_filter_and_generate_oligos.md`
- `docs/examples/generated/guides_export_csv_and_protocol.md`

## `gentle_mcp` (MCP stdio server)

`gentle_mcp` starts a Model Context Protocol server over stdio.

Current tools:

- `capabilities`
- `state_summary`
- `op` (apply one operation; requires `confirm=true`)
- `workflow` (apply one workflow; requires `confirm=true`)
- `help`
- `ui_intents` (discover deterministic UI-intent contracts)
- `ui_intent` (run deterministic `ui open|focus` intent resolution path)
- `ui_prepared_genomes` (run deterministic prepared-genome query path)
- `ui_latest_prepared` (resolve latest prepared genome for one species)

Tool-parity rule:

- MCP UI-intent tools execute through the same shared parser/executor route as
  CLI shell commands (`ui ...`), and return the same structured payload
  contracts.
- No MCP-only biology/UI logic branches are allowed.

Run:

```bash
cargo run --bin gentle_mcp
cargo run --bin gentle_mcp -- --state path/to/project.gentle.json
```

Minimum MCP JSON-RPC flow:

1. Initialize:

```json
{"jsonrpc":"2.0","id":1,"method":"initialize","params":{"protocolVersion":"2025-06-18"}}
```

2. List tools:

```json
{"jsonrpc":"2.0","id":2,"method":"tools/list","params":{}}
```

3. Call one tool:

```json
{
  "jsonrpc":"2.0",
  "id":3,
  "method":"tools/call",
  "params":{
    "name":"ui_intent",
    "arguments":{
      "action":"open",
      "target":"prepared-references",
      "catalog_path":"assets/genomes.json",
      "species":"human",
      "latest":true
    }
  }
}
```

`ui_intent` arguments:

- required:
  - `action`: `open|focus`
  - `target`: one of:
    - `prepared-references`
    - `prepare-reference-genome`
    - `retrieve-genome-sequence`
    - `blast-genome-sequence`
    - `import-genome-track`
    - `agent-assistant`
    - `prepare-helper-genome`
    - `retrieve-helper-sequence`
    - `blast-helper-sequence`
- optional:
  - `state_path`
  - `genome_id`
  - `helpers`
  - `catalog_path`
  - `cache_dir`
  - `filter`
  - `species`
  - `latest`

`ui_prepared_genomes` arguments:

- optional: `state_path`, `helpers`, `catalog_path`, `cache_dir`, `filter`,
  `species`, `latest`

`ui_latest_prepared` arguments:

- required: `species`
- optional: `state_path`, `helpers`, `catalog_path`, `cache_dir`

`ui_intents` arguments:

- optional: `state_path` (accepted for API symmetry)

MCP tool result envelope behavior:

- `result.isError = false`: success; inspect `result.structuredContent`.
- `result.isError = true`: deterministic tool-level error; inspect
  `result.content[0].text`.
- UI-intent structured payloads match shared shell schemas:
  - `gentle.ui_intents.v1`
  - `gentle.ui_intent.v1`
  - `gentle.ui_prepared_genomes.v1`
  - `gentle.ui_latest_prepared.v1`

Mutating tool safety:

- `op` and `workflow` reject execution unless `confirm=true` is explicitly set.
- UI-intent tools are currently non-mutating; if a routed command ever reports
  state mutation unexpectedly, MCP returns an explicit tool error.

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
2. `write_gb(seq, filename)`
   - Writes a sequence object to a GenBank file.
3. `load_project(filename)` / `save_project(state, filename)`
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
10. `import_pool(project, input_pool_json, [prefix])`
   - Imports a `.pool.gentle.json` artifact into `project` via shared adapter logic.
   - Returns table with `state`, `state_changed`, and `output`.
11. `sync_rebase(input, output, commercial_only)`
   - Parses REBASE/Bairoch input and writes a REBASE resource JSON snapshot.
   - `output` and `commercial_only` are optional.
12. `sync_jaspar(input, output)`
   - Parses JASPAR PFM text and writes motif resource JSON snapshot.
   - `output` is optional.
13. `list_reference_genomes(catalog_path)`
    - Lists genome IDs from the genome catalog.
14. `is_reference_genome_prepared(genome_id, catalog_path, cache_dir)`
    - Checks whether a genome cache is prepared and indexed.
15. `list_reference_genome_genes(genome_id, catalog_path, cache_dir)`
    - Lists indexed genes from prepared genome cache.
16. `prepare_genome(project, genome_id, catalog_path, cache_dir)`
    - Convenience wrapper around engine `PrepareGenome`.
17. `extract_genome_region(project, genome_id, chromosome, start_1based, end_1based, output_id, catalog_path, cache_dir)`
    - Convenience wrapper around engine `ExtractGenomeRegion`.
18. `extract_genome_gene(project, genome_id, gene_query, occurrence, output_id, catalog_path, cache_dir)`
    - Convenience wrapper around engine `ExtractGenomeGene`.
19. `blast_reference_genome(genome_id, query_sequence, [max_hits], [task], [catalog_path], [cache_dir])`
    - Runs BLAST (`blastn`/`blastn-short`) against a prepared reference genome.
20. `blast_helper_genome(helper_id, query_sequence, [max_hits], [task], [catalog_path], [cache_dir])`
    - Same as `blast_reference_genome`, but defaults to helper catalog context.
21. `set_parameter(project, name, value)`
    - Convenience wrapper for engine `SetParameter`.
22. `set_vcf_display_filter(project, opts)`
    - Convenience wrapper that updates one or more VCF display-filter parameters
      via `SetParameter` in one call.
23. `list_agent_systems([catalog_path])`
    - Lists configured agent systems from catalog JSON.
    - `catalog_path` is optional (defaults to `assets/agent_systems.json`).
24. `ask_agent_system(project, system_id, prompt, [catalog_path], [allow_auto_exec], [execute_all], [execute_indices], [include_state_summary])`
    - Invokes one configured agent system through shared shell execution.
    - Returns table with `state`, `state_changed`, and `output`.
    - `project` may be `nil` to use an empty/default project state.

### Lua example

Use generated adapter snippets to stay synchronized with canonical workflow JSON:

- `docs/examples/generated/load_and_digest_pgex.md`
- `docs/examples/generated/load_branch_reverse_complement_pgex_fasta.md`
- `docs/examples/generated/guides_filter_and_generate_oligos.md`
- `docs/examples/generated/guides_export_csv_and_protocol.md`

## File format expectations

Current CLI workflows rely on sequence files supported by internal loaders:

- GenBank
- EMBL
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
- NCBI GenBank XML (`GBSet/GBSeq`)
  - currently rejected with explicit diagnostics for unsupported XML dialects
    (for example `INSDSet/INSDSeq`)

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
cargo run --bin gentle_cli -- help
cargo run --bin gentle_cli -- help candidates generate
cargo run --bin gentle_cli -- help --format json
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
cargo run --bin gentle_cli -- render-pool-gel-svg - digest.gel.svg --containers container-3,container-8
cargo run --bin gentle_cli -- render-pool-gel-svg - digest.gel.svg --arrangement arrangement-2
cargo run --bin gentle_cli -- render-gel-svg - digest.gel.svg --arrangement arrangement-2
cargo run --bin gentle_cli -- arrange-serial container-3,container-8 --id arrangement-2 --name "Digest run A" --ladders "NEB 100bp DNA Ladder"
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
cargo run --bin gentle_cli -- agents list
cargo run --bin gentle_cli -- agents list --catalog assets/agent_systems.json
cargo run --bin gentle_cli -- agents ask builtin_echo --prompt "summarize current project state"
cargo run --bin gentle_cli -- agents ask builtin_echo --prompt "ask: Which sequence should I use?" --execute-index 1
cargo run --bin gentle_cli -- agents ask local_llama_compat --prompt "summarize project context" --base-url http://localhost:11964 --model deepseek-r1:8b
cargo run --bin gentle_cli -- op '{"PrepareGenome":{"genome_id":"ToyGenome","catalog_path":"catalog.json"}}'
cargo run --bin gentle_cli -- op '{"ExtractGenomeRegion":{"genome_id":"ToyGenome","chromosome":"chr1","start_1based":1001,"end_1based":1600,"output_id":"toy_chr1_1001_1600","catalog_path":"catalog.json"}}'
cargo run --bin gentle_cli -- op '{"ExtractGenomeGene":{"genome_id":"ToyGenome","gene_query":"MYGENE","occurrence":1,"output_id":"toy_mygene","catalog_path":"catalog.json"}}'
cargo run --bin gentle_cli -- genomes list --catalog assets/genomes.json
cargo run --bin gentle_cli -- genomes list --catalog assets/helper_genomes.json
cargo run --bin gentle_cli -- genomes validate-catalog --catalog assets/genomes.json
cargo run --bin gentle_cli -- genomes status "Human GRCh38 Ensembl 116" --catalog assets/genomes.json --cache-dir data/genomes
cargo run --bin gentle_cli -- genomes genes "Human GRCh38 Ensembl 116" --catalog assets/genomes.json --cache-dir data/genomes --filter "^TP53$" --biotype protein_coding --limit 20
cargo run --bin gentle_cli -- genomes prepare "Human GRCh38 Ensembl 116" --catalog assets/genomes.json --cache-dir data/genomes --timeout-secs 3600
cargo run --bin gentle_cli -- genomes blast "Human GRCh38 Ensembl 116" ACGTACGTACGT --task blastn-short --max-hits 10 --catalog assets/genomes.json --cache-dir data/genomes
cargo run --bin gentle_cli -- genomes extract-region "Human GRCh38 Ensembl 116" 1 1000000 1001500 --output-id grch38_chr1_slice --catalog assets/genomes.json --cache-dir data/genomes
cargo run --bin gentle_cli -- genomes extract-gene "Human GRCh38 Ensembl 116" TP53 --occurrence 1 --output-id grch38_tp53 --catalog assets/genomes.json --cache-dir data/genomes
cargo run --bin gentle_cli -- tracks import-bed grch38_tp53 data/chipseq/peaks.bed.gz --name H3K27ac --min-score 10 --clear-existing
cargo run --bin gentle_cli -- tracks import-bigwig grch38_tp53 data/chipseq/signal.bw --name ATAC --min-score 0.2 --clear-existing
cargo run --bin gentle_cli -- tracks import-vcf grch38_tp53 data/variants/sample.vcf.gz --name Variants --min-score 20 --clear-existing
cargo run --bin gentle_cli -- helpers list
cargo run --bin gentle_cli -- helpers validate-catalog
cargo run --bin gentle_cli -- helpers status "Plasmid pUC19 (local)"
cargo run --bin gentle_cli -- helpers prepare "Plasmid pUC19 (local)" --cache-dir data/helper_genomes --timeout-secs 600
cargo run --bin gentle_cli -- helpers genes "Plasmid pUC19 (local)" --filter bla --limit 20
cargo run --bin gentle_cli -- helpers blast "Plasmid pUC19 (local)" ACGTACGTACGT --task blastn-short --max-hits 10 --cache-dir data/helper_genomes
cargo run --bin gentle_cli -- candidates generate sgrnas chr1_window --length 20 --step 1 --feature-kind gene --max-distance 500 --limit 5000
cargo run --bin gentle_cli -- candidates score sgrnas gc_balance "100 * (gc_fraction - at_fraction)"
cargo run --bin gentle_cli -- candidates score-weighted sgrnas priority --term gc_fraction:0.7:max --term distance_to_seq_start_bp:0.3:min --normalize
cargo run --bin gentle_cli -- candidates top-k sgrnas sgrnas_top --metric priority --k 100 --direction max
cargo run --bin gentle_cli -- candidates pareto sgrnas sgrnas_front --objective gc_fraction:max --objective distance_to_seq_start_bp:min --max-candidates 200
cargo run --bin gentle_cli -- candidates filter sgrnas sgrnas_q95 --metric gc_balance --min-quantile 0.95
cargo run --bin gentle_cli -- candidates macro @candidate_flow.gsh
cargo run --bin gentle_cli -- candidates macro --transactional --file candidate_flow.gsh
cargo run --bin gentle_cli -- candidates template-put scan_tp53 --script 'generate ${set_name} ${seq_id} --length ${len} --step 1' --param set_name --param seq_id=grch38_tp53 --param len=20
cargo run --bin gentle_cli -- candidates template-run scan_tp53 --bind set_name=tp53_candidates --transactional
cargo run --bin gentle_cli -- guides list
cargo run --bin gentle_cli -- guides put tp73_guides --json '[{"guide_id":"g1","seq_id":"tp73","start_0based":100,"end_0based_exclusive":120,"strand":"+","protospacer":"GACCTGTTGACGATGTTCCA","pam":"AGG","nuclease":"SpCas9","cut_offset_from_protospacer_start":17}]'
cargo run --bin gentle_cli -- guides filter tp73_guides --config '{"gc_min":0.3,"gc_max":0.7,"avoid_u6_terminator_tttt":true}' --output-set tp73_guides_pass
cargo run --bin gentle_cli -- guides oligos-generate tp73_guides lenti_bsmbi_u6_default --apply-5prime-g-extension --output-oligo-set tp73_lenti --passed-only
cargo run --bin gentle_cli -- guides oligos-export tp73_guides exports/tp73_guides.csv --format csv_table --oligo-set tp73_lenti
cargo run --bin gentle_cli -- guides protocol-export tp73_guides exports/tp73_guides.protocol.txt --oligo-set tp73_lenti
cargo run --bin gentle_cli -- shell 'macros template-list'
cargo run --bin gentle_cli -- shell 'macros template-import assets/cloning_patterns.json'
cargo run --bin gentle_cli -- shell 'macros template-import assets/cloning_patterns_catalog'
cargo run --bin gentle_cli -- shell 'macros run --transactional --file cloning_flow.gsh'
cargo run --bin gentle_cli -- shell 'set-param vcf_display_pass_only true'
cargo run --bin gentle_cli -- shell 'set-param vcf_display_required_info_keys ["AF","DP"]'
cargo run --bin gentle_cli -- shell 'set-param tfbs_display_min_llr_quantile 0.95'
```

You can pass JSON from a file with `@file.json`.

Global CLI options:

- `--state PATH`: use a non-default project state file
- `--progress` or `--progress-stderr`: print live progress events to `stderr`
- `--progress-stdout`: print live progress events to `stdout`
- `--allow-screenshots`: currently rejected (screenshot bridge disabled by security policy)

Current progress events include TFBS annotation updates, genome-prepare
updates (download/index phases), and genome-track import updates.
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
    - `render-pool-gel-svg IDS|'-' OUTPUT.svg [--ladders NAME[,NAME]] [--containers ID[,ID]] [--arrangement ARR_ID]`
    - `render-gel-svg IDS|'-' OUTPUT.svg [--ladders NAME[,NAME]] [--containers ID[,ID]] [--arrangement ARR_ID]`
    - `arrange-serial CONTAINER_IDS [--id ARR_ID] [--name TEXT] [--ladders NAME[,NAME]]`
    - `ladders list [--molecule dna|rna] [--filter TEXT]`
    - `ladders export OUTPUT.json [--molecule dna|rna] [--filter TEXT]`
    - `export-pool IDS OUTPUT.pool.gentle.json [HUMAN_ID]`
    - `import-pool INPUT.pool.gentle.json [PREFIX]`
    - `resources sync-rebase INPUT.withrefm_or_URL [OUTPUT.rebase.json] [--commercial-only]`
    - `resources sync-jaspar INPUT.jaspar_or_URL [OUTPUT.motifs.json]`
    - `agents list [--catalog PATH]`
    - `agents ask SYSTEM_ID --prompt TEXT [--catalog PATH] [--base-url URL] [--model MODEL] [--timeout-secs N] [--connect-timeout-secs N] [--read-timeout-secs N] [--max-retries N] [--max-response-bytes N] [--allow-auto-exec] [--execute-all] [--execute-index N ...] [--no-state-summary]`
    - `genomes list [--catalog PATH]`
    - `genomes validate-catalog [--catalog PATH]`
    - `genomes status GENOME_ID [--catalog PATH] [--cache-dir PATH]`
    - `genomes genes GENOME_ID [--catalog PATH] [--cache-dir PATH] [--filter REGEX] [--biotype NAME] [--limit N] [--offset N]`
    - `genomes prepare GENOME_ID [--catalog PATH] [--cache-dir PATH] [--timeout-secs N]`
    - `genomes blast GENOME_ID QUERY_SEQUENCE [--max-hits N] [--task blastn-short|blastn] [--catalog PATH] [--cache-dir PATH]`
    - `genomes extract-region GENOME_ID CHR START END [--output-id ID] [--catalog PATH] [--cache-dir PATH]`
    - `genomes extract-gene GENOME_ID QUERY [--occurrence N] [--output-id ID] [--catalog PATH] [--cache-dir PATH]`
    - `helpers list [--catalog PATH]`
    - `helpers validate-catalog [--catalog PATH]`
    - `helpers status HELPER_ID [--catalog PATH] [--cache-dir PATH]`
    - `helpers genes HELPER_ID [--catalog PATH] [--cache-dir PATH] [--filter REGEX] [--biotype NAME] [--limit N] [--offset N]`
    - `helpers prepare HELPER_ID [--catalog PATH] [--cache-dir PATH] [--timeout-secs N]`
    - `helpers blast HELPER_ID QUERY_SEQUENCE [--max-hits N] [--task blastn-short|blastn] [--catalog PATH] [--cache-dir PATH]`
    - `helpers extract-region HELPER_ID CHR START END [--output-id ID] [--catalog PATH] [--cache-dir PATH]`
    - `helpers extract-gene HELPER_ID QUERY [--occurrence N] [--output-id ID] [--catalog PATH] [--cache-dir PATH]`
    - `tracks import-bed SEQ_ID PATH [--name NAME] [--min-score N] [--max-score N] [--clear-existing]`
    - `tracks import-bigwig SEQ_ID PATH [--name NAME] [--min-score N] [--max-score N] [--clear-existing]`
    - `tracks import-vcf SEQ_ID PATH [--name NAME] [--min-score N] [--max-score N] [--clear-existing]`
    - `macros run [--transactional] [--file PATH | SCRIPT_OR_@FILE]`
    - `macros template-list`
    - `macros template-show TEMPLATE_NAME`
    - `macros template-put TEMPLATE_NAME (--script SCRIPT_OR_@FILE|--file PATH) [--description TEXT] [--details-url URL] [--param NAME|NAME=DEFAULT ...]`
    - `macros template-delete TEMPLATE_NAME`
    - `macros template-import PATH`
    - `macros template-run TEMPLATE_NAME [--bind KEY=VALUE ...] [--transactional]`
    - `routines list [--catalog PATH] [--family NAME] [--status NAME] [--tag TAG] [--query TEXT]`
    - `candidates list`
    - `candidates delete SET_NAME`
    - `candidates generate SET_NAME SEQ_ID --length N [--step N] [--feature-kind KIND] [--feature-label-regex REGEX] [--max-distance N] [--feature-geometry feature_span|feature_parts|feature_boundaries] [--feature-boundary any|five_prime|three_prime|start|end] [--strand-relation any|same|opposite] [--limit N]`
    - `candidates generate-between-anchors SET_NAME SEQ_ID --length N (--anchor-a-pos N|--anchor-a-json JSON) (--anchor-b-pos N|--anchor-b-json JSON) [--step N] [--limit N]`
    - `candidates show SET_NAME [--limit N] [--offset N]`
    - `candidates metrics SET_NAME`
    - `candidates score SET_NAME METRIC_NAME EXPRESSION`
    - `candidates score-distance SET_NAME METRIC_NAME [--feature-kind KIND] [--feature-label-regex REGEX] [--feature-geometry feature_span|feature_parts|feature_boundaries] [--feature-boundary any|five_prime|three_prime|start|end] [--strand-relation any|same|opposite]`
    - `candidates score-weighted SET_NAME METRIC_NAME --term METRIC:WEIGHT[:max|min] [--term ...] [--normalize|--no-normalize]`
    - `candidates top-k INPUT_SET OUTPUT_SET --metric METRIC_NAME --k N [--direction max|min] [--tie-break seq_start_end|seq_end_start|length_ascending|length_descending|sequence_lexicographic]`
    - `candidates pareto INPUT_SET OUTPUT_SET --objective METRIC[:max|min] [--objective ...] [--max-candidates N] [--tie-break seq_start_end|seq_end_start|length_ascending|length_descending|sequence_lexicographic]`
    - `candidates filter INPUT_SET OUTPUT_SET --metric METRIC_NAME [--min N] [--max N] [--min-quantile Q] [--max-quantile Q]`
    - `candidates set-op union|intersect|subtract LEFT_SET RIGHT_SET OUTPUT_SET`
    - `candidates macro [--transactional] [--file PATH | SCRIPT_OR_@FILE]`
    - `candidates template-list`
    - `candidates template-show TEMPLATE_NAME`
    - `candidates template-put TEMPLATE_NAME (--script SCRIPT_OR_@FILE|--file PATH) [--description TEXT] [--details-url URL] [--param NAME|NAME=DEFAULT ...]`
    - `candidates template-delete TEMPLATE_NAME`
    - `candidates template-run TEMPLATE_NAME [--bind KEY=VALUE ...] [--transactional]`
    - `guides list`
    - `guides show GUIDE_SET_ID [--limit N] [--offset N]`
    - `guides put GUIDE_SET_ID (--json JSON|@FILE|--file PATH)`
    - `guides delete GUIDE_SET_ID`
    - `guides filter GUIDE_SET_ID [--config JSON|@FILE] [--config-file PATH] [--output-set GUIDE_SET_ID]`
    - `guides filter-show GUIDE_SET_ID`
    - `guides oligos-generate GUIDE_SET_ID TEMPLATE_ID [--apply-5prime-g-extension] [--output-oligo-set ID] [--passed-only]`
    - `guides oligos-list [--guide-set GUIDE_SET_ID]`
    - `guides oligos-show OLIGO_SET_ID`
    - `guides oligos-export GUIDE_SET_ID OUTPUT_PATH [--format csv_table|plate_csv|fasta] [--plate 96|384] [--oligo-set ID]`
    - `guides protocol-export GUIDE_SET_ID OUTPUT_PATH [--oligo-set ID] [--no-qc]`
    - `set-param NAME JSON_VALUE`
    - `op <operation-json-or-@file>`
    - `workflow <workflow-json-or-@file>`
    - `screenshot-window OUTPUT.png` (currently disabled by security policy)
    - `help [COMMAND ...] [--format text|json|markdown] [--interface all|cli-direct|cli-shell|gui-shell|js|lua|mcp]`
  - Use single quotes around JSON payloads to preserve whitespace:
    - `gentle_cli shell 'workflow {"run_id":"r1","ops":[]}'`
  - Structured help export for automation:
    - `gentle_cli help --format json`
    - `gentle_cli help --format markdown`

Screenshot bridge:

- Current status: disabled by security policy.
- `screenshot-window` currently returns a deterministic disabled-policy message.
- `--allow-screenshots` is rejected by argument parsing.
- Manual documentation updates remain the current approach for image artifacts.

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
- `render-pool-gel-svg IDS|'-' OUTPUT.svg [--ladders NAME[,NAME]] [--containers ID[,ID]] [--arrangement ARR_ID]`
  - Calls engine operation `RenderPoolGelSvg`.
  - Use `IDS` as a comma-separated sequence-id list, or pass `-`/`_` when using `--containers` or `--arrangement`.
  - `--containers` renders one lane per container ID.
  - `--arrangement` renders lanes from a stored serial arrangement.
  - Optional `--ladders` overrides auto ladder selection.
  - If `--ladders` is omitted, engine auto-selects one or two ladders based on
    pool bp range.
- `render-gel-svg ...`
  - Alias for `render-pool-gel-svg ...` with identical semantics.
- `arrange-serial CONTAINER_IDS [--id ARR_ID] [--name TEXT] [--ladders NAME[,NAME]]`
  - Calls engine operation `CreateArrangementSerial`.
  - Persists a serial lane setup that can be reused by `render-pool-gel-svg --arrangement`.

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

Agent bridge commands:

- `agents list [--catalog PATH]`
  - Lists configured agent systems from catalog JSON.
  - Default catalog path: `assets/agent_systems.json`.
  - Includes availability status (`available`, `availability_reason`) so callers
    can skip systems that are not currently runnable.
- `agents ask SYSTEM_ID --prompt TEXT [--catalog PATH] [--base-url URL] [--model MODEL] [--timeout-secs N] [--connect-timeout-secs N] [--read-timeout-secs N] [--max-retries N] [--max-response-bytes N] [--allow-auto-exec] [--execute-all] [--execute-index N ...] [--no-state-summary]`
  - Invokes one configured agent system and returns message/questions/suggested shell commands.
  - `--base-url` sets a per-request runtime endpoint override (maps to
    `GENTLE_AGENT_BASE_URL`) for native transports.
  - `--model` sets a per-request model override (maps to `GENTLE_AGENT_MODEL`)
    for native transports.
  - `--timeout-secs` sets per-request timeout override (maps to
    `GENTLE_AGENT_TIMEOUT_SECS`) for stdio/native transports.
  - `--connect-timeout-secs` sets per-request HTTP connect timeout (maps to
    `GENTLE_AGENT_CONNECT_TIMEOUT_SECS`) for native transports.
  - `--read-timeout-secs` sets per-request read timeout (maps to
    `GENTLE_AGENT_READ_TIMEOUT_SECS`) for stdio/native transports.
  - `--max-retries` sets per-request transient retry budget (maps to
    `GENTLE_AGENT_MAX_RETRIES`; `0` means no retries).
  - `--max-response-bytes` caps adapter response size per attempt (maps to
    `GENTLE_AGENT_MAX_RESPONSE_BYTES`).
  - for `native_openai_compat`, model must resolve to a concrete value (catalog
    model or `--model`); `unspecified` is rejected.
  - `--no-state-summary` disables project-context injection in the request.
  - External adapter responses must be strict `gentle.agent_response.v1` JSON.
  - Unknown canonical response fields are rejected (extensions must use `x_` or `x-` prefix).
  - Adapter calls retry transient failures with exponential backoff before returning an error.
  - Suggested commands are executed only when explicitly selected
    (`--execute-all`, `--execute-index`) or when `--allow-auto-exec` is enabled
    and the suggestion intent is `auto`.
  - Recursive `agents ask` execution from suggested commands is blocked by design.
  - Failures use deterministic error prefixes for scripting, e.g.
    `AGENT_ADAPTER_UNAVAILABLE`, `AGENT_ADAPTER_TRANSIENT`,
    `AGENT_RESPONSE_PARSE`, `AGENT_RESPONSE_VALIDATION`.
  - Supported transports:
    - `builtin_echo` (offline/demo)
    - `external_json_stdio` (local bridge executable)
    - `native_openai` (built-in OpenAI HTTP adapter, requires `OPENAI_API_KEY`)
    - `native_openai_compat` (built-in local OpenAI-compatible adapter for
      Jan/Msty/Ollama-style `/chat/completions` endpoints; key optional)
      - endpoint host/port come from catalog `base_url` (or `--base-url` if set);
        GENtle does not silently switch to a different host/port

Prompting users should send to agents:

- Prefer compact structured prompts so local models stay reliable.
- Use this template:

```text
Objective:
<one clear goal>

Context:
<sequence/genome/helper ids>

Inputs:
- anchors/coordinates:
- feature labels/kinds:

Constraints:
- length:
- GC range:
- motifs/sites:
- strand assumptions:

Output wanted:
- exact commands
- validation checklist

Execution policy:
chat-only | ask-before-run | allow-auto-exec
```

- For local models without cloning background, preload:
  - `docs/ai_cloning_primer.md`
  - `docs/ai_task_playbooks.md`
  - `docs/ai_prompt_contract.md`
  - `docs/examples/ai_cloning_examples.md`
  - optional: `docs/ai_glossary_extensions.json`

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
  - Also reports optional mass metadata:
    `nucleotide_length_bp`, `molecular_mass_da`, `molecular_mass_source`.
- `genomes genes GENOME_ID [--catalog PATH] [--cache-dir PATH] [--filter REGEX] [--biotype NAME] [--limit N] [--offset N]`
  - Lists indexed genes from prepared cache (paged by `--limit`/`--offset`).
  - `--filter` is a case-insensitive regular expression.
  - `--biotype` can be repeated to constrain matches to selected biotypes.
- `genomes prepare GENOME_ID [--catalog PATH] [--cache-dir PATH] [--timeout-secs N]`
  - Runs engine `PrepareGenome`.
- `genomes blast GENOME_ID QUERY_SEQUENCE [--max-hits N] [--task blastn-short|blastn] [--catalog PATH] [--cache-dir PATH]`
  - Runs `blastn` against prepared genome cache/index.
  - `--task` defaults to `blastn-short`; accepted values: `blastn-short`, `blastn`.
- `genomes extract-region GENOME_ID CHR START END [--output-id ID] [--catalog PATH] [--cache-dir PATH]`
  - Runs engine `ExtractGenomeRegion`.
- `genomes extract-gene GENOME_ID QUERY [--occurrence N] [--output-id ID] [--catalog PATH] [--cache-dir PATH]`
  - Runs engine `ExtractGenomeGene`.
- `genomes extend-anchor SEQ_ID 5p|3p LENGTH_BP [--output-id ID] [--catalog PATH] [--cache-dir PATH]`
  - Runs engine `ExtendGenomeAnchor`.
  - Extends an already genome-anchored sequence in-silico on contextual `5'` or `3'`.

Helper convenience commands:

- `helpers list [--catalog PATH]`
  - Same behavior as `genomes list`, but defaults to `assets/helper_genomes.json`.
- `helpers validate-catalog [--catalog PATH]`
  - Same behavior as `genomes validate-catalog`, with helper-catalog default.
- `helpers status HELPER_ID [--catalog PATH] [--cache-dir PATH]`
  - Same behavior as `genomes status`, with helper-catalog default
    (including length/mass metadata fields).
- `helpers genes HELPER_ID [--catalog PATH] [--cache-dir PATH] [--filter REGEX] [--biotype NAME] [--limit N] [--offset N]`
  - Same behavior as `genomes genes`, with helper-catalog default.
- `helpers prepare HELPER_ID [--catalog PATH] [--cache-dir PATH] [--timeout-secs N]`
  - Same behavior as `genomes prepare`, with helper-catalog default.
  - `--timeout-secs N`: optional prepare-job timebox.
- `helpers blast HELPER_ID QUERY_SEQUENCE [--max-hits N] [--task blastn-short|blastn] [--catalog PATH] [--cache-dir PATH]`
  - Same behavior as `genomes blast`, with helper-catalog default.
- `helpers extract-region HELPER_ID CHR START END [--output-id ID] [--catalog PATH] [--cache-dir PATH]`
  - Same behavior as `genomes extract-region`, with helper-catalog default.
- `helpers extract-gene HELPER_ID QUERY [--occurrence N] [--output-id ID] [--catalog PATH] [--cache-dir PATH]`
  - Same behavior as `genomes extract-gene`, with helper-catalog default.
- `helpers extend-anchor SEQ_ID 5p|3p LENGTH_BP [--output-id ID] [--catalog PATH] [--cache-dir PATH]`
  - Same behavior as `genomes extend-anchor`, with helper-catalog default.

Workflow macro commands (`gentle_cli shell 'macros ...'`):

- `macros run [--transactional] [--file PATH | SCRIPT_OR_@FILE]`
  - Executes semicolon/newline-separated shell statements.
  - Supports transactional rollback (`--transactional`) when any statement fails.
  - Designed for full cloning workflows through `op ...` and `workflow ...`
    statements (Digest/Ligation/PCR/ExtractRegion/container ops, etc.).
- `macros template-list`
  - Lists persisted workflow macro templates.
- `macros template-show TEMPLATE_NAME`
  - Shows one persisted workflow template definition.
- `macros template-put TEMPLATE_NAME (--script SCRIPT_OR_@FILE|--file PATH) [--description TEXT] [--details-url URL] [--param NAME|NAME=DEFAULT ...]`
  - Creates/updates a named workflow macro template in project metadata.
  - Placeholders in script use `${param_name}` and must be declared via `--param`.
  - Optional `--details-url URL` records external protocol/reference details for
    catalog display.
- `macros template-delete TEMPLATE_NAME`
  - Deletes one persisted workflow template.
- `macros template-import PATH`
  - Imports workflow macro templates from:
    - one pack JSON file (`gentle.cloning_patterns.v1`)
    - one single-template JSON file (`gentle.cloning_pattern_template.v1`)
    - one directory tree (recursive `*.json` import)
  - If one template fails validation, no imported template changes are kept.
- `macros template-run TEMPLATE_NAME [--bind KEY=VALUE ...] [--transactional]`
  - Expands a named template with provided bindings/defaults, then executes it as
    a workflow macro script.

Typed routine catalog command (`gentle_cli routines ...` or `gentle_cli shell 'routines ...'`):

- `routines list [--catalog PATH] [--family NAME] [--status NAME] [--tag TAG] [--query TEXT]`
  - Lists typed cloning routines from catalog JSON (`gentle.cloning_routines.v1`).
  - `--family`, `--status`, `--tag`: exact case-insensitive filters.
  - `--query`: case-insensitive substring match across id/title/family/status/template/tags/summary.
  - Default catalog path: `assets/cloning_routines.json`.

Shipped starter assets:

- Legacy pack:
  - `assets/cloning_patterns.json` (`gentle.cloning_patterns.v1`)
- Hierarchical catalog (one template file per macro):
  - `assets/cloning_patterns_catalog/**/*.json`
  - each file schema: `gentle.cloning_pattern_template.v1`
  - folder hierarchy is used by GUI `Patterns` menu hierarchy
- Typed routine manifest:
  - `assets/cloning_routines.json`
  - schema: `gentle.cloning_routines.v1`
- Import commands:
  - `gentle_cli shell 'macros template-import assets/cloning_patterns.json'`
  - `gentle_cli shell 'macros template-import assets/cloning_patterns_catalog'`

Candidate-set commands (`gentle_cli candidates ...` and `gentle_cli shell 'candidates ...'`):

- `candidates list`
  - Lists available candidate sets from project metadata.
- `candidates delete SET_NAME`
  - Removes one candidate set.
- `candidates generate SET_NAME SEQ_ID --length N [--step N] [--feature-kind KIND] [--feature-label-regex REGEX] [--max-distance N] [--feature-geometry feature_span|feature_parts|feature_boundaries] [--feature-boundary any|five_prime|three_prime|start|end] [--strand-relation any|same|opposite] [--limit N]`
  - Creates a candidate set from fixed-length windows on `SEQ_ID`.
  - `--feature-kind` can be repeated to constrain nearest-feature context.
  - `--feature-geometry`, `--feature-boundary`, and `--strand-relation`
    control how directed feature distance targets are selected.
- `candidates generate-between-anchors SET_NAME SEQ_ID --length N (--anchor-a-pos N|--anchor-a-json JSON) (--anchor-b-pos N|--anchor-b-json JSON) [--step N] [--limit N]`
  - Creates a candidate set from fixed-length windows constrained to the
    interval between two local sequence anchors on `SEQ_ID`.
  - Anchor JSON accepts the same local `SequenceAnchor` schema used by
    `ExtractAnchoredRegion` (`Position` or `FeatureBoundary` with
    `Start|End|Middle`).
- `candidates show SET_NAME [--limit N] [--offset N]`
  - Pages records (`sequence`, `coordinates`, `metrics`) from one set.
- `candidates metrics SET_NAME`
  - Lists available metric names in one set.
- `candidates score SET_NAME METRIC_NAME EXPRESSION`
  - Computes a derived metric expression for all records in a set.
- `candidates score-distance SET_NAME METRIC_NAME [--feature-kind KIND] [--feature-label-regex REGEX] [--feature-geometry feature_span|feature_parts|feature_boundaries] [--feature-boundary any|five_prime|three_prime|start|end] [--strand-relation any|same|opposite]`
  - Computes nearest-feature distance metric with optional feature filters.
- `candidates score-weighted SET_NAME METRIC_NAME --term METRIC:WEIGHT[:max|min] [--term ...] [--normalize|--no-normalize]`
  - Computes one weighted objective metric from existing metrics.
  - `--term` can be repeated; default direction is `max`.
  - `--normalize` is enabled by default (min-max scaling per term).
- `candidates top-k INPUT_SET OUTPUT_SET --metric METRIC_NAME --k N [--direction max|min] [--tie-break seq_start_end|seq_end_start|length_ascending|length_descending|sequence_lexicographic]`
  - Materializes explicit top-k selection for a metric.
  - deterministic tie-break policy avoids unstable ordering.
- `candidates pareto INPUT_SET OUTPUT_SET --objective METRIC[:max|min] [--objective ...] [--max-candidates N] [--tie-break seq_start_end|seq_end_start|length_ascending|length_descending|sequence_lexicographic]`
  - Materializes a Pareto frontier for multi-objective optimization.
  - optionally truncates with deterministic tie-break (`--max-candidates`).
- `candidates filter INPUT_SET OUTPUT_SET --metric METRIC_NAME [--min N] [--max N] [--min-quantile Q] [--max-quantile Q]`
  - Creates `OUTPUT_SET` by value and/or quantile thresholds.
- `candidates set-op union|intersect|subtract LEFT_SET RIGHT_SET OUTPUT_SET`
  - Creates set algebra output from two sets.
- `candidates macro SCRIPT_OR_@FILE`
  - Runs multiple candidate statements in order (semicolon/newline separated).
  - Nested macro calls are rejected.
- `candidates template-list`
  - Lists persisted candidate macro templates.
- `candidates template-show TEMPLATE_NAME`
  - Shows one persisted template definition.
- `candidates template-put TEMPLATE_NAME (--script SCRIPT_OR_@FILE|--file PATH) [--description TEXT] [--details-url URL] [--param NAME|NAME=DEFAULT ...]`
  - Creates/updates a named template in project metadata.
  - placeholders in script use `${param_name}` and must be declared via `--param`.
  - Optional `--details-url URL` records external protocol/reference details for
    catalog display.
- `candidates template-delete TEMPLATE_NAME`
  - Deletes one persisted template.
- `candidates template-run TEMPLATE_NAME [--bind KEY=VALUE ...] [--transactional]`
  - Expands a named template with provided bindings/defaults, then executes it as
    a candidate macro script.

Guide-design commands (`gentle_cli guides ...` and `gentle_cli shell 'guides ...'`):

- `guides list`
  - Lists persisted guide sets from guide-design metadata.
- `guides show GUIDE_SET_ID [--limit N] [--offset N]`
  - Pages guide rows from one guide set.
- `guides put GUIDE_SET_ID (--json JSON|@FILE|--file PATH)`
  - Creates/replaces one guide set from JSON array payload.
  - Input payload is `Vec<GuideCandidate>` objects.
- `guides delete GUIDE_SET_ID`
  - Deletes one guide set.
- `guides filter GUIDE_SET_ID [--config JSON|@FILE] [--config-file PATH] [--output-set GUIDE_SET_ID]`
  - Applies practical guide constraints and persists a filter report.
  - Optional `--output-set` writes passed guides into a new set.
- `guides filter-show GUIDE_SET_ID`
  - Returns persisted practical-filter report (reasons, warnings, metrics).
- `guides oligos-generate GUIDE_SET_ID TEMPLATE_ID [--apply-5prime-g-extension] [--output-oligo-set ID] [--passed-only]`
  - Generates guide oligos from selected template.
  - Built-ins include `lenti_bsmbi_u6_default` and `plain_forward_reverse`.
- `guides oligos-list [--guide-set GUIDE_SET_ID]`
  - Lists persisted oligo sets (optionally filtered by source guide set).
- `guides oligos-show OLIGO_SET_ID`
  - Shows one oligo set with generated forward/reverse oligos.
- `guides oligos-export GUIDE_SET_ID OUTPUT_PATH [--format csv_table|plate_csv|fasta] [--plate 96|384] [--oligo-set ID]`
  - Exports oligo rows for ordering or plate layouts.
- `guides protocol-export GUIDE_SET_ID OUTPUT_PATH [--oligo-set ID] [--no-qc]`
  - Exports a human-readable wet-lab protocol text for generated oligos.

Notes:

- In-memory candidate sets persist in
  `ProjectState.metadata["candidate_sets"]` (`gentle.candidate_sets.v1`).
- On save, candidate sets are externalized into a sidecar index + JSONL files;
  project metadata stores a reference schema (`gentle.candidate_sets.ref.v1`).
- On load, sidecar-backed candidate metadata is hydrated automatically.
- `list/show/metrics` are read-only commands.
- `macros run/template-put/template-delete/template-run` and
  `delete/generate/generate-between-anchors/score/score-distance/score-weighted/top-k/pareto/filter/set-op/macro/template-put/template-delete/template-run`
  mutate state and are available through CLI shell and GUI shell.

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

Note:

- `ExtractAnchoredRegion.anchor` is a local sequence anchor (in-sequence
  coordinate resolver), not genome-assembly provenance anchoring.

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

Apply practical design-constraint filters (GC bounds, homopolymer cap, U6 `TTTT`
avoidance, forbidden motifs):

```json
{"FilterByDesignConstraints":{"inputs":["g1","g2","g3"],"gc_min":0.30,"gc_max":0.70,"max_homopolymer_run":4,"reject_ambiguous_bases":true,"avoid_u6_terminator_tttt":true,"forbidden_motifs":["GAATTC"],"unique":false,"output_prefix":"design_pick"}}
```

Notes:

- `reject_ambiguous_bases` defaults to `true` if omitted.
- `avoid_u6_terminator_tttt` defaults to `true` if omitted.
- `forbidden_motifs` supports IUPAC and is checked on both strands.

Generate a candidate set from windows on one sequence:

```json
{"GenerateCandidateSet":{"set_name":"sgrna_windows","seq_id":"grch38_tp53","length_bp":20,"step_bp":1,"feature_kinds":["gene"],"feature_label_regex":"^TP53$","max_distance_bp":5000,"limit":10000}}
```

Generate a candidate set only between two local anchors:

```json
{"GenerateCandidateSetBetweenAnchors":{"set_name":"sgrna_between","seq_id":"grch38_tp53","anchor_a":{"FeatureBoundary":{"feature_kind":"gene","feature_label":"TP53","boundary":"Start","occurrence":0}},"anchor_b":{"FeatureBoundary":{"feature_kind":"gene","feature_label":"TP53","boundary":"End","occurrence":0}},"length_bp":20,"step_bp":1,"limit":10000}}
```

Add a derived metric expression:

```json
{"ScoreCandidateSetExpression":{"set_name":"sgrna_windows","metric":"gc_balance","expression":"100 * (gc_fraction - at_fraction)"}}
```

Add distance-to-feature score:

```json
{"ScoreCandidateSetDistance":{"set_name":"sgrna_windows","metric":"distance_to_cds_bp","feature_kinds":["CDS"],"feature_label_regex":null}}
```

Add a weighted multi-objective score:

```json
{"ScoreCandidateSetWeightedObjective":{"set_name":"sgrna_windows","metric":"priority_score","objectives":[{"metric":"gc_fraction","weight":0.7,"direction":"maximize"},{"metric":"distance_to_cds_bp","weight":0.3,"direction":"minimize"}],"normalize_metrics":true}}
```

Select top-k by one metric:

```json
{"TopKCandidateSet":{"input_set":"sgrna_windows","output_set":"sgrna_top20","metric":"priority_score","k":20,"direction":"maximize","tie_break":"seq_start_end"}}
```

Compute Pareto frontier for multiple objectives:

```json
{"ParetoFrontierCandidateSet":{"input_set":"sgrna_windows","output_set":"sgrna_frontier","objectives":[{"metric":"gc_fraction","direction":"maximize"},{"metric":"distance_to_cds_bp","direction":"minimize"}],"max_candidates":100,"tie_break":"seq_start_end"}}
```

Filter by absolute value and quantiles into a new explicit set:

```json
{"FilterCandidateSet":{"input_set":"sgrna_windows","output_set":"sgrna_windows_top","metric":"gc_balance","min":-20.0,"max":20.0,"min_quantile":0.10,"max_quantile":0.90}}
```

Intersect/union/subtract explicit candidate sets:

```json
{"CandidateSetOp":{"op":"intersect","left_set":"sgrna_windows_top","right_set":"near_tss","output_set":"sgrna_final"}}
{"CandidateSetOp":{"op":"union","left_set":"set_a","right_set":"set_b","output_set":"set_union"}}
{"CandidateSetOp":{"op":"subtract","left_set":"set_a","right_set":"set_b","output_set":"set_a_minus_b"}}
```

Delete a candidate set:

```json
{"DeleteCandidateSet":{"set_name":"sgrna_windows"}}
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

Set regulatory-overlay max linear view span threshold (`50000` recommended for anchored genome maps):

```json
{"SetParameter":{"name":"regulatory_feature_max_view_span_bp","value":50000}}
```

Set TFBS display filtering parameters shared by GUI and SVG export:

```json
{"SetParameter":{"name":"show_tfbs","value":true}}
{"SetParameter":{"name":"tfbs_display_use_llr_bits","value":true}}
{"SetParameter":{"name":"tfbs_display_min_llr_bits","value":0.0}}
{"SetParameter":{"name":"tfbs_display_use_llr_quantile","value":true}}
{"SetParameter":{"name":"tfbs_display_min_llr_quantile","value":0.95}}
{"SetParameter":{"name":"tfbs_display_use_true_log_odds_bits","value":false}}
{"SetParameter":{"name":"tfbs_display_min_true_log_odds_bits","value":0.0}}
{"SetParameter":{"name":"tfbs_display_use_true_log_odds_quantile","value":false}}
{"SetParameter":{"name":"tfbs_display_min_true_log_odds_quantile","value":0.95}}
```

Supported TFBS display parameter names:

- `show_tfbs`
- `tfbs_display_use_llr_bits`
- `tfbs_display_min_llr_bits`
- `tfbs_display_use_llr_quantile`
- `tfbs_display_min_llr_quantile` (quantile in range `0.0..1.0`)
- `tfbs_display_use_true_log_odds_bits`
- `tfbs_display_min_true_log_odds_bits`
- `tfbs_display_use_true_log_odds_quantile`
- `tfbs_display_min_true_log_odds_quantile` (quantile in range `0.0..1.0`)

Set VCF display filtering parameters shared by GUI and SVG export:

```json
{"SetParameter":{"name":"vcf_display_pass_only","value":true}}
{"SetParameter":{"name":"vcf_display_use_min_qual","value":true}}
{"SetParameter":{"name":"vcf_display_min_qual","value":30.0}}
{"SetParameter":{"name":"vcf_display_required_info_keys","value":["AF","DP"]}}
```

Supported VCF display parameter names:

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
- `vcf_display_required_info_keys` (string CSV or string array)

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

Extend an anchored sequence on the contextual 5' side by 250 bp:

```json
{"ExtendGenomeAnchor":{"seq_id":"grch38_tp53","side":"five_prime","length_bp":250,"output_id":"grch38_tp53_ext5","catalog_path":"assets/genomes.json","cache_dir":"data/genomes"}}
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

Import VCF variants (`.vcf` / `.vcf.gz`) onto a genome-anchored sequence:

```json
{"ImportGenomeVcfTrack":{"seq_id":"grch38_tp53","path":"data/variants/sample.vcf.gz","track_name":"Variants","min_score":20.0,"max_score":null,"clear_existing":false}}
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
- Catalog entries may also include optional physical parameters:
  `nucleotide_length_bp` and `molecular_mass_da`.
- If `molecular_mass_da` is omitted but nucleotide length is available
  (directly or from prepared FASTA index), GENtle estimates dsDNA molecular
  mass and labels it as `estimated_from_nucleotide_length`.
- `ExtractGenomeRegion` expects the genome to have been prepared already.
- `ExtractGenomeGene` also expects prepared cache and gene index.
- `genomes/helpers blast` expects prepared cache and a BLAST index.
  If index files are missing, GENtle tries to build them on demand.
- BLAST executable overrides:
  - `GENTLE_MAKEBLASTDB_BIN` (default: `makeblastdb`)
  - `GENTLE_BLASTN_BIN` (default: `blastn`)
- `ImportGenomeBedTrack` expects `seq_id` to be a sequence created by
  `ExtractGenomeRegion`, `ExtractGenomeGene`, or `ExtendGenomeAnchor`
  (genome-anchored provenance).
- `ImportGenomeBigWigTrack` expects the same genome-anchored `seq_id`.
- `ImportGenomeVcfTrack` expects the same genome-anchored `seq_id`.
- BED import accepts local `.bed` and `.bed.gz` files.
- BigWig import accepts local `.bw` and `.bigWig` files and uses
  `bigWigToBedGraph` (override with `GENTLE_BIGWIG_TO_BEDGRAPH_BIN`).
- VCF import accepts local `.vcf` and `.vcf.gz` files.
- For VCF import, `min_score` / `max_score` filter on VCF `QUAL`.
- `ExtractGenomeRegion`, `ExtractGenomeGene`, and `ExtendGenomeAnchor` append extraction provenance
  records into `ProjectState.metadata["provenance"]["genome_extractions"]`
  (genome id, coordinates/query, source descriptors, and checksums when present).
- If `catalog_path` is omitted, engine default catalog is `assets/genomes.json`.
- Bundled `assets/genomes.json` currently includes Human GRCh38 (Ensembl 113 and 116),
  Mouse GRCm39 Ensembl 116, Rat GRCr8 Ensembl 116, Saccharomyces cerevisiae
  S288c (Ensembl 113 and 116), and `LocalProject` (backed by
  `test_files/fixtures/genomes/AB011549.2.fa` +
  `test_files/fixtures/genomes/AB011549.2.gb`).
- `cache_dir` is optional. If omitted, catalog/default cache settings are used.
- `chromosome` accepts exact contig names and also tolerates `chr` prefix
  differences (`1` vs `chr1`).
- For `ExtractGenomeGene`, `occurrence` is 1-based among matching records.
- For `ExtendGenomeAnchor`, `side` is contextual to anchor strand.
  On anchor strand `-`, `5'` increases physical genomic position.

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

- `import-pool` is currently an adapter-level utility contract (CLI/GUI shared
  shell + JS/Lua wrappers) and not yet an engine operation.

## Error behavior

- REPLs print runtime errors without exiting automatically.
- Invalid file paths or unsupported content produce load/write errors.
- Some functions assume computed annotations are available after `load_dna`.
