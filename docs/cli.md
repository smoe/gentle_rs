# GENtle CLI Manual

This page documents command-line entry points for GENtle.

## Overview

GENtle currently provides six binaries:

- `gentle`: graphical desktop app (GUI)
- `gentle_cli`: JSON operation/workflow CLI for automation and AI tools
- `gentle_js`: interactive JavaScript shell (optional build target; requires
  feature `js-interface`)
- `gentle_lua`: interactive Lua shell (optional build target; requires feature
  `lua-interface`)
- `gentle_examples_docs`: generates adapter snippets and tutorial artifacts from canonical protocol examples
- `gentle_mcp`: MCP stdio server (guarded mutating + UI-intent parity baseline;
  includes standardized capability discovery via `tools/list`,
  `capabilities`, and `help`)

In addition, the GUI includes an embedded `Shell` panel that uses the same
shared shell parser/executor as `gentle_cli shell`.

Structured command glossary:

- `docs/glossary.json` is the machine-readable command glossary used for
  per-command help rendering (`help ...`) and catalog export
  (`help --format json|markdown`).
- Shell-only dry-run catalog preview routes are available through the shared
  shell, for example:
  - `gentle_cli shell 'genomes preview-ensembl-specs --catalog assets/genomes.json'`
  - `gentle_cli shell 'helpers preview-ensembl-specs --catalog assets/helper_genomes.json'`

Structured workflow examples:

- canonical source files: `docs/examples/workflows/*.json`
- schema: `gentle.workflow_example.v1`
- includes GUI-first parity skeleton:
  - `docs/examples/workflows/tp73_promoter_luciferase_assay_planning.json`
  - merged canonical tutorial: `docs/tutorial/tp73_promoter_luciferase_gui.md`
- on-demand snippet generation (CLI/shell/JS/Lua):
  - `cargo run --bin gentle_examples_docs -- generate`
- validation only:
  - `cargo run --bin gentle_examples_docs -- --check`
- tutorial generation:
  - `cargo run --bin gentle_examples_docs -- tutorial-generate`
- tutorial drift check:
  - `cargo run --bin gentle_examples_docs -- tutorial-check`
- test gating:
  - `always` examples run in default tests
  - `online` examples run only with `GENTLE_TEST_ONLINE=1`
  - set `GENTLE_SKIP_REMOTE_TESTS=1` to force-skip all remote-resource tests
    (takes precedence over `GENTLE_TEST_ONLINE`)
  - `skip` examples are syntax-checked only

Architecture invariant: all adapters/frontends above route cloning/business
behavior through the same shared engine.

## ClawBio/OpenClaw skill scaffold

GENtle ships a copy-ready ClawBio skill scaffold at:

- `integrations/clawbio/skills/gentle-cloning/`
- ready-to-paste catalog object:
  `integrations/clawbio/skills/gentle-cloning/catalog_entry.json`

Purpose:

- expose deterministic GENtle CLI routes as a ClawBio skill,
- keep reproducibility artifacts (`report.md`, `result.json`,
  `reproducibility/*`) for each run.

Quick usage (inside the scaffold directory):

```bash
python gentle_cloning.py --demo --output /tmp/gentle_clawbio_demo
python gentle_cloning.py --input examples/request_capabilities.json --output /tmp/gentle_clawbio_caps
```

When copied into a ClawBio checkout (`skills/gentle-cloning/`):

```bash
python clawbio.py run gentle-cloning --demo
python clawbio.py run gentle-cloning --input skills/gentle-cloning/examples/request_workflow_file.json --output /tmp/gentle_clawbio_run
```

Command resolution order used by the wrapper:

1. `--gentle-cli` explicit command string
2. `GENTLE_CLI_CMD` environment variable
3. `gentle_cli` on `PATH`
4. repository-local fallback: `cargo run --quiet --bin gentle_cli --`

## Python module wrapper (`gentle_py`)

GENtle also ships a thin Python wrapper over `gentle_cli`:

- module path: `integrations/python/gentle_py/`
- quick docs: `integrations/python/README.md`

Quick usage:

```bash
python3 - <<'PY'
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path("integrations/python").resolve()))
from gentle_py import GentleClient

client = GentleClient(state_path=".gentle_state.json")
print(json.dumps(client.capabilities(), indent=2))
PY
```

Deterministic API surface:

- `capabilities()`
- `state_summary()`
- `op(operation)`
- `workflow(workflow|workflow_path=...)`
- `shell(line, expect_json=False)`
- `render_dotplot_svg(seq_id, dotplot_id, output_svg, ...)`

CLI resolution order:

1. constructor `cli_cmd`
2. `GENTLE_CLI_CMD`
3. `gentle_cli` on `PATH`
4. repository fallback `cargo run --quiet --bin gentle_cli --`

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

Primer-design report capability status:

- `gentle_cli`: supported via shared-shell `primers ...` commands and direct forwarding (`gentle_cli primers ...`), backed by `DesignPrimerPairs` and `DesignQpcrAssays` plus non-mutating ROI seed helpers (`primers seed-from-feature`, `primers seed-from-splicing`) and persisted report inspect/export helpers
- `gentle_js`: supported via `apply_operation` (`DesignPrimerPairs`, `DesignQpcrAssays`) plus shared-shell execution for report listing/show/export
- `gentle_lua`: supported via `apply_operation` (`DesignPrimerPairs`, `DesignQpcrAssays`) plus shared-shell execution for report listing/show/export

Dotplot/flexibility capability status:

- `gentle_cli`: supported via shared-shell/direct commands:
  - `dotplot compute|list|show|render-svg`
  - `render-dotplot-svg`
  - `transcripts derive`
  - `flex compute|list|show`
  - `splicing-refs derive`
  - `align compute`
  backed by `ComputeDotplot`, `DeriveTranscriptSequences`,
  `ComputeFlexibilityTrack`, `DeriveSplicingReferences`,
  `AlignSequences`, and `RenderDotplotSvg`.
  - `dotplot compute` supports self and pairwise modes via
    `--mode self_forward|self_reverse_complement|pair_forward|pair_reverse_complement`
    with optional `--reference-seq`, `--ref-start`, and `--ref-end`.
  - `transcripts derive` supports:
    - full-sequence derivation (all `mRNA`/`transcript` features)
    - selected feature derivation (`--feature-id`)
    - splicing-scope constrained derivation from one seed feature
      (`--scope ...` with exactly one `--feature-id`)
- `gentle_js`: baseline support via `apply_operation` (`ComputeDotplot`,
  `DeriveTranscriptSequences`, `ComputeFlexibilityTrack`,
  `DeriveSplicingReferences`, `AlignSequences`) plus
  `render_dotplot_svg(...)` convenience wrapper over `RenderDotplotSvg`.
- `gentle_lua`: baseline support via `apply_operation` (`ComputeDotplot`,
  `DeriveTranscriptSequences`, `ComputeFlexibilityTrack`,
  `DeriveSplicingReferences`, `AlignSequences`) plus
  `render_dotplot_svg(...)` convenience wrapper over `RenderDotplotSvg`.
- `gentle_py`: baseline support via `op(...)` plus
  `render_dotplot_svg(...)` convenience wrapper over `RenderDotplotSvg`.

RNA-read interpretation capability status (Nanopore cDNA phase-1):

- `gentle_cli`: supported via shared-shell/direct commands:
  - `rna-reads interpret`
  - `rna-reads align-report`
  - `rna-reads list-reports`
  - `rna-reads show-report`
  - `rna-reads summarize-gene-support`
  - `rna-reads inspect-alignments`
  - `rna-reads export-report`
  - `rna-reads export-hits-fasta`
  - `rna-reads export-sample-sheet`
  - `rna-reads export-paths-tsv`
  - `rna-reads export-abundance-tsv`
  - `rna-reads export-score-density-svg`
  - `rna-reads export-alignments-tsv`
  - `rna-reads export-alignment-dotplot-svg`
  backed by `InterpretRnaReads`, `AlignRnaReadReport`,
  `ListRnaReadReports`, `ShowRnaReadReport`,
  `SummarizeRnaReadGeneSupport`, `ExportRnaReadReport`,
  `ExportRnaReadHitsFasta`, `ExportRnaReadSampleSheet`,
  `ExportRnaReadExonPathsTsv`, `ExportRnaReadExonAbundanceTsv`, `ExportRnaReadScoreDensitySvg`,
  `ExportRnaReadAlignmentsTsv`, and `ExportRnaReadAlignmentDotplotSvg`.
  Input supports FASTA plus gzipped FASTA (`.fa/.fasta` and `.fa.gz/.fasta.gz`).
  Concatenated gzip members are accepted for gzipped FASTA input as well.
  Progress output includes periodic `progress rna-reads ...` lines during
  `apply_with_progress` runs.
  Phase split:
  - `interpret`: seed-filter pass (Nanopore phase-1)
    - `align-report`: retained-hit alignment pass (phase-2) that updates
      mapping fields, MSA-eligibility counters, exon-transition rows, and
      exon/junction abundance frequencies in the persisted report; retained hits
      are re-ranked by alignment-aware retention rank after alignment.
    - supports explicit row filtering via
      `--record-indices i,j,k` (0-based stored `record_index` values);
      when provided, this overrides `--selection`.
    - default adapter behavior now uses `selection=all` so rescued retained
      rows also receive round-2 similarity/coverage scores.
  - `summarize-gene-support`: non-mutating target-gene cohort summary over one
    saved aligned RNA-read report.
    - `--gene GENE` is required and repeatable; matches are case-insensitive
      against the same transcript group-label logic used by the splicing view.
    - optional `--record-indices i,j,k` restricts the base aligned cohort to an
      exact saved-report subset before target-gene filtering.
    - `--complete-rule near|strict|exact` controls the fragment-vs-complete
      cohort split; strict and exact counts are still reported explicitly in the
      payload.
    - output is JSON on stdout by default; `--output PATH` writes the exact same
      payload to disk.
    - payload tables include:
      - per-cohort exon support
      - ordered exon-pair support (including skipped pairs like `1->3`)
      - neighboring direct-transition support only (for example `1->2`)
  - `inspect-alignments`: non-mutating ranked alignment inspection
    over persisted report hits, with optional structured subset controls:
    `--effect-filter`, `--sort`, `--search`, and `--record-indices`.
  - `export-alignments-tsv`: non-mutating ranked alignment-row TSV export for
    downstream table-based triage.
  - `export-alignment-dotplot-svg`: non-mutating dotplot-like scatter export
    (coverage vs identity, score-colored points) for aligned report hits.
- `gentle_js`: baseline support via `apply_operation` for the same operation
  family.
- `gentle_lua`: baseline support via `apply_operation` for the same operation
  family.
- GUI workflow parity:
  - Splicing Expert `Nanopore cDNA interpretation` now exposes
    `Prepare Workflow Op` and `Copy Workflow JSON` so the exact
    `InterpretRnaReads` payload can be executed via regular workflow routes.
  - Example:
    - `gentle_cli workflow @docs/examples/workflows/load_branch_reverse_complement_pgex_fasta.json`
    - or paste copied JSON directly:
      `gentle_cli workflow '{"run_id":"workflow_rna_reads_tp73_ncbi_cdna_srr32957124","ops":[{"InterpretRnaReads":{"seq_id":"tp73.ncbi","seed_feature_id":0,"profile":"nanopore_cdna_v1","input_path":"reads.fa.gz","input_format":"fasta","scope":"all_overlapping_both_strands","seed_filter":{"kmer_len":10,"seed_stride_bp":1,"min_seed_hit_fraction":0.30,"min_weighted_seed_hit_fraction":0.05,"min_unique_matched_kmers":12,"max_median_transcript_gap":4.0,"min_chain_consistency_fraction":0.40,"min_confirmed_exon_transitions":1,"min_transition_support_fraction":0.05,"cdna_poly_t_flip_enabled":true,"poly_t_prefix_min_bp":18},"align_config":{"band_width_bp":24,"min_identity_fraction":0.60,"max_secondary_mappings":0},"report_id":"cdna_reads"}}]}'`

Sequencing-confirmation capability status (called-read plus trace-aware engine baseline):

- `gentle_cli`: supported via shared-shell commands plus the current GUI
  specialist:
  - `seq-trace import PATH [--trace-id ID] [--seq-id ID]`
  - `seq-trace list [SEQ_ID]`
  - `seq-trace show TRACE_ID`
  - `seq-confirm run EXPECTED_SEQ_ID [--baseline BASELINE_SEQ_ID] --reads ID[,ID...] [--trace-id ID]... [--trace-ids ID[,ID...]]`
  - `seq-confirm list-reports [EXPECTED_SEQ_ID]`
  - `seq-confirm show-report REPORT_ID`
  - `seq-confirm export-report REPORT_ID OUTPUT.json`
  - `seq-confirm export-support-tsv REPORT_ID OUTPUT.tsv`
  - `seq-primer suggest EXPECTED_SEQ_ID [--primers ID[,ID...]] [--confirmation-report REPORT_ID] [--min-3prime-anneal-bp N] [--predicted-read-length-bp N]`
  raw-trace intake is backed by `ImportSequencingTrace`,
  `ListSequencingTraces`, and `ShowSequencingTrace`.
  backed by `ConfirmConstructReads`, `ListSequencingConfirmationReports`,
  `ShowSequencingConfirmationReport`, `ExportSequencingConfirmationReport`, and
  `ExportSequencingConfirmationSupportTsv`. Sequencing-primer overlays are
  backed by `SuggestSequencingPrimers`.
  Current scope:
  - expected construct plus already-loaded read sequences
  - imported ABI/AB1/SCF traces can now participate directly in
    `seq-confirm run` without first materializing a project sequence
  - default full-span confirmation when no explicit targets are supplied
  - explicit junction targets via `--junction LEFT_END_0BASED`
    and `--junction-flank N`
  - forward and reverse-complement read evaluation through the same shared
    engine report contract
  - raw ABI/AB1/SCF intake stores trace evidence separately from
    confirmation reports:
    - preserves file-supplied called bases
    - preserves called-base confidence arrays when available
    - preserves peak locations when available
    - preserves raw per-channel chromatogram curves and clip windows when
      available from the source trace
    - keeps one imported trace store that can be listed/shown without running
      confirmation
  - trace-aware confirmation keeps one shared report model:
    - sequence-backed and trace-backed evidence rows land in the same
      confirmation report
    - trace-backed rows expose their `trace_id` and evidence kind
    - target support/contradiction ids can now reflect imported trace ids
    - optional `--baseline` context lets confirmation classify one locus as:
      intended edit, reference reversion, unexpected difference, or
      low-confidence/insufficient evidence
    - `seq-confirm show-report` now exposes the persisted `baseline_seq_id`
      plus per-variant rows for expected edits inferred from the
      baseline-vs-expected diff
  - sequencing-primer suggestion now has two modes through the same
    `SuggestSequencingPrimers` operation:
    - overlay existing primer ids onto the expected construct
    - or, with only `--confirmation-report`, emit fresh primer proposals for
      unresolved loci that still lack a good existing hit
    - when existing primer hits do exist but sit outside the preferred
      sequencing window, the report now keeps both the best existing hit and
      the top fresh proposal
  - GUI parity is now available through `Patterns -> Sequencing Confirmation...`
    and command palette `Sequencing Confirmation`, backed by the same
    `ConfirmConstructReads` / report-store path:
    - optional baseline sequence input is available in the run form
    - the review pane now includes a variant list and chromatogram curve view
      for trace-backed loci
    - sequencing-primer overlays can now be suggested in the same specialist
      from already-loaded primer sequence ids
    - the overlay report can optionally annotate which saved-report targets and
      variant loci fall inside each predicted primer-derived read span
    - unresolved targets and variant loci now also get one shared
      recommendation row naming the best existing primer hit to clarify them
    - older traces without stored curve arrays stay usable for confirmation,
      but the GUI asks you to re-import them before chromatogram review
  Not yet included in this phase:
  - GUI-side raw trace import
  - chromatogram editing or base re-calling
  - full whole-trace browsing beyond the current variant-focused window
  - lineage/artifact projection for confirmation reports
- `gentle_js`: baseline support via `apply_operation` for the same operation
  family.
- `gentle_lua`: baseline support via `apply_operation` for the same operation
  family.
- `gentle_py`: baseline support via `op(...)` for the same operation family.

Feature-expert SVG parity status:

- `gentle_cli`: supported via direct commands and shared-shell command path
  (`inspect-feature-expert`, `render-feature-expert-svg`)
- `gentle_js`: supported via `render_feature_expert_svg(...)` wrapper over
  shared engine operation `RenderFeatureExpertSvg`
- `gentle_lua`: supported via `render_feature_expert_svg(...)` wrapper over
  shared engine operation `RenderFeatureExpertSvg`
- splicing expert SVG now includes:
  - explicit junction transition support labels/table
  - frequency-encoded transcript-vs-exon matrix cell coloring
  - predicted exon-to-exon transition matrix with frequency-encoded cell
    coloring
  - CDS flank phase coloring (`0/1/2`) on exon left/right edges when
    transcript `cds_ranges_1based` are available
  - exon-length modulo-3 (`len%3`) cues on exon headers (heuristic frame cue)
  - shared output semantics across CLI/JS/Lua because all routes call the same
    engine renderer

UniProt mapping capability status:

- shared shell (`gentle_cli shell`, GUI shell): supported via `uniprot ...` commands
  - `uniprot fetch QUERY [--entry-id ID]`
  - `uniprot import-swissprot PATH [--entry-id ID]`
  - `uniprot list` / `uniprot show ENTRY_ID`
  - `uniprot map ENTRY_ID SEQ_ID [--projection-id ID] [--transcript ID]`
  - `uniprot projection-list [--seq SEQ_ID]`
  - `uniprot projection-show PROJECTION_ID`
- shared shell (`gentle_cli shell`, GUI shell): GenBank accession import
  - `genbank fetch ACCESSION [--as-id ID]`
- shared shell (`gentle_cli shell`, GUI shell): dbSNP-guided annotated region extraction
  - `dbsnp fetch RS_ID GENOME_ID [--flank-bp N] [--output-id ID] [--annotation-scope none|core|full] [--max-annotation-features N] [--catalog PATH] [--cache-dir PATH]`
- engine operations behind those commands:
  - `FetchUniprotSwissProt`, `ImportUniprotSwissProt`, `ProjectUniprotToGenome`,
    `FetchGenBankAccession`, `FetchDbSnpRegion`
- `ImportUniprotEntrySequence` is currently disabled (`Unsupported`) because
  first-class protein sequence windows are deferred; use UniProt metadata +
  projection routes instead.

## Build and run

From the repository root:

```bash
cargo run --bin gentle
cargo run --features js-interface --bin gentle_js
cargo run --features lua-interface --bin gentle_lua
cargo run --bin gentle_cli -- capabilities
cargo run --bin gentle_examples_docs -- --check
cargo run --bin gentle_examples_docs -- tutorial-check
cargo run --bin gentle_mcp -- --help
cargo run --bin gentle -- --version
cargo run --bin gentle_cli -- --version
```

For optimized builds:

```bash
cargo run --release --bin gentle
cargo run --release --features js-interface --bin gentle_js
cargo run --release --features lua-interface --bin gentle_lua
cargo run --release --bin gentle_cli -- capabilities
cargo run --release --bin gentle_examples_docs -- --check
cargo run --release --bin gentle_examples_docs -- tutorial-check
cargo run --release --bin gentle_mcp -- --help
cargo run --release --bin gentle -- --version
cargo run --release --bin gentle_cli -- --version
```

Important: with `cargo run`, arguments for GENtle binaries must come after `--`.
Example: `cargo run --bin gentle -- --version`.

JavaScript and Lua shells are compile-time optional so default `cargo check`
and `cargo build` do not pull in their heavier runtime dependencies. The Python
wrapper in `integrations/python/gentle_py` remains separate from Cargo
features.

Release-build policy:

- default local builds stay lean unless JS/Lua support is requested explicitly
- release packaging builds enable `--features script-interfaces`
- for a local optimized build that matches release behavior, use
  `cargo build --release --features script-interfaces`

## Protocol-first example source

Canonical workflow examples are adapter-neutral JSON files:

- source: `docs/examples/workflows/*.json`
- generated adapter snippets: `docs/examples/generated/*.md`
- tutorial runtime manifest: `docs/tutorial/manifest.json`
  (generated from `docs/tutorial/sources/`)
- committed generated tutorial output: `docs/tutorial/generated/`

Regenerate snippets on demand:

```bash
cargo run --bin gentle_examples_docs -- generate
```

Validate example files and schema without writing output:

```bash
cargo run --bin gentle_examples_docs -- --check
```

Generate tutorial markdown + retained artifacts:

```bash
cargo run --bin gentle_examples_docs -- tutorial-generate
```

Validate committed tutorial-generated output:

```bash
cargo run --bin gentle_examples_docs -- tutorial-check
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
This binary is available when Cargo feature `js-interface` is enabled.

```bash
cargo run --features js-interface --bin gentle_js
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
18. `extract_genome_gene(state, genome_id, gene_query, occurrence, output_id, catalog_path, cache_dir, annotation_scope, max_annotation_features, extract_mode, promoter_upstream_bp)`
    - Convenience wrapper around engine `ExtractGenomeGene`.
19. `blast_reference_genome(genome_id, query_sequence, max_hits, task, catalog_path, cache_dir, options_json)`
    - Runs BLAST (`blastn`/`blastn-short`) against a prepared reference genome.
    - `options_json` is optional and may override any quick option (`max_hits`, `task`) plus threshold fields.
20. `blast_helper_genome(helper_id, query_sequence, max_hits, task, catalog_path, cache_dir, options_json)`
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
25. `render_dotplot_svg(state, seq_id, dotplot_id, output_svg, options)`
    - Convenience wrapper around engine `RenderDotplotSvg`.
    - `options` supports:
      - `flex_track_id` (or `flexTrackId`)
      - `display_density_threshold`
      - `display_intensity_gain`

### JavaScript example

Use generated adapter snippets to stay synchronized with canonical workflow JSON:

- `docs/examples/generated/load_and_digest_pgex.md`
- `docs/examples/generated/load_branch_reverse_complement_pgex_fasta.md`
- `docs/examples/generated/guides_filter_and_generate_oligos.md`
- `docs/examples/generated/digest_ligation_extract_region_minimal.md`
- `docs/examples/generated/contribute_gentle_development_baseline.md`
- `docs/examples/generated/guides_export_csv_and_protocol.md`

## `gentle_mcp` (MCP stdio server)

`gentle_mcp` starts a Model Context Protocol server over stdio.
It provides both tool execution (`tools/call`) and capability
discovery/negotiation (`tools/list`, `capabilities`, `help`).

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
- `blast_async_start` (start async BLAST job through shared shell route)
- `blast_async_status` (poll async BLAST job status and optional report)
- `blast_async_cancel` (request cancellation for async BLAST job)
- `blast_async_list` (list known async BLAST jobs in scope)

Tool-parity rule:

- MCP UI-intent tools execute through the same shared parser/executor route as
  CLI shell commands (`ui ...`), and return the same structured payload
  contracts.
- No MCP-only biology/UI logic branches are allowed.

Data-boundary and trust model:

- `gentle_mcp` reads/writes project state via the resolved `state_path`
  (default `.gentle_state.json`).
- It does not attach to another running GUI process's in-memory engine state.
- Unsaved GUI-only edits are therefore not visible to MCP until they are
  persisted to disk and read from that same state path.
- If MCP is pointed at the same project state file used by active GUI/CLI
  workflows, that shared visibility is intentional.
- Treat connected MCP clients as trusted local operators; `confirm=true`
  expresses explicit mutating intent but is not an authentication boundary.

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
    - `pcr-design`
    - `sequencing-confirmation`
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

`blast_async_start` arguments:

- required: `genome_id`, `query_sequence`
- optional: `state_path`, `helpers`, `max_hits`, `task`, `options_json`,
  `catalog_path`, `cache_dir`

`blast_async_status` arguments:

- required: `job_id`
- optional: `helpers`, `with_report`

`blast_async_cancel` arguments:

- required: `job_id`
- optional: `helpers`

`blast_async_list` arguments:

- optional: `helpers`

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
This binary is available when Cargo feature `lua-interface` is enabled.

```bash
cargo run --features lua-interface --bin gentle_lua
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
18. `extract_genome_gene(project, genome_id, gene_query, occurrence, output_id, catalog_path, cache_dir, annotation_scope, max_annotation_features, extract_mode, promoter_upstream_bp)`
    - Convenience wrapper around engine `ExtractGenomeGene`.
19. `blast_reference_genome(genome_id, query_sequence, [max_hits], [task], [catalog_path], [cache_dir], [options_json])`
    - Runs BLAST (`blastn`/`blastn-short`) against a prepared reference genome.
    - `options_json` is optional and may override any quick option (`max_hits`, `task`) plus threshold fields.
20. `blast_helper_genome(helper_id, query_sequence, [max_hits], [task], [catalog_path], [cache_dir], [options_json])`
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
25. `render_dotplot_svg(project, seq_id, dotplot_id, output_svg, [flex_track_id], [display_density_threshold], [display_intensity_gain])`
    - Convenience wrapper around engine `RenderDotplotSvg`.

### Lua example

Use generated adapter snippets to stay synchronized with canonical workflow JSON:

- `docs/examples/generated/load_and_digest_pgex.md`
- `docs/examples/generated/load_branch_reverse_complement_pgex_fasta.md`
- `docs/examples/generated/guides_filter_and_generate_oligos.md`
- `docs/examples/generated/digest_ligation_extract_region_minimal.md`
- `docs/examples/generated/contribute_gentle_development_baseline.md`
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
cargo run --bin gentle_cli -- op op.json
cargo run --bin gentle_cli -- workflow '<workflow-json>'
cargo run --bin gentle_cli -- workflow docs/examples/workflows/load_branch_reverse_complement_pgex_fasta.json
cargo run --bin gentle_cli -- workflow @docs/examples/workflows/load_branch_reverse_complement_pgex_fasta.json
cargo run --bin gentle_cli -- workflow @docs/examples/workflows/rna_reads_interpret_cdna_tp73_template.json
cargo run --bin gentle_cli -- --progress op '<operation-json>'
cargo run --bin gentle_cli -- --progress-stdout workflow '<workflow-json>'
cargo run --bin gentle_cli -- export-state state.json
cargo run --bin gentle_cli -- import-state state.json
cargo run --bin gentle_cli -- save-project project.gentle.json
cargo run --bin gentle_cli -- load-project project.gentle.json
cargo run --bin gentle_cli -- render-svg pgex linear pgex.linear.svg
cargo run --bin gentle_cli -- render-svg pgex circular pgex.circular.svg
cargo run --bin gentle_cli -- render-dotplot-svg tp73_cdna dotplot_primary tp73.dotplot.svg
cargo run --bin gentle_cli -- render-dotplot-svg tp73_cdna dotplot_primary tp73.dotplot.svg --flex-track flex_primary --display-threshold 0.2 --intensity-gain 1.8
cargo run --bin gentle_cli -- render-rna-svg rna_seq rna.secondary.svg
cargo run --bin gentle_cli -- rna-info rna_seq
cargo run --bin gentle_cli -- render-lineage-svg lineage.svg
cargo run --bin gentle_cli -- protocol-cartoon list
cargo run --bin gentle_cli -- protocol-cartoon render-svg gibson.two_fragment gibson.protocol.svg
cargo run --bin gentle_cli -- protocol-cartoon render-svg pcr.assay.pair pcr_pair.protocol.svg
cargo run --bin gentle_cli -- protocol-cartoon render-svg pcr.assay.pair.with_tail pcr_pair_with_tail.protocol.svg
cargo run --bin gentle_cli -- protocol-cartoon render-svg pcr.oe.substitution pcr_oe_substitution.protocol.svg
cargo run --bin gentle_cli -- protocol-cartoon render-svg pcr.assay.qpcr qpcr.protocol.svg
cargo run --bin gentle_cli -- protocol-cartoon render-template-svg docs/examples/protocol_cartoon/demo_template.json demo.protocol.svg
cargo run --bin gentle_cli -- protocol-cartoon template-validate docs/examples/protocol_cartoon/demo_template.json
cargo run --bin gentle_cli -- protocol-cartoon render-with-bindings docs/examples/protocol_cartoon/demo_template.json docs/examples/protocol_cartoon/demo_bindings.json demo.bound.protocol.svg
cargo run --bin gentle_cli -- protocol-cartoon template-export gibson.two_fragment gibson.template.json
cargo run --bin gentle_cli -- gibson preview @docs/examples/plans/gibson_destination_first_single_insert.json --output gibson.preview.json
cargo run --bin gentle_cli -- gibson apply @docs/examples/plans/gibson_destination_first_single_insert.json
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
cargo run --bin gentle_cli -- op '{"ExtractGenomeRegion":{"genome_id":"ToyGenome","chromosome":"chr1","start_1based":1001,"end_1based":1600,"output_id":"toy_chr1_1001_1600","annotation_scope":"core","catalog_path":"catalog.json"}}'
cargo run --bin gentle_cli -- op '{"ExtractGenomeGene":{"genome_id":"ToyGenome","gene_query":"MYGENE","occurrence":1,"output_id":"toy_mygene","catalog_path":"catalog.json"}}'
cargo run --bin gentle_cli -- op '{"DesignInsertionPrimerPairs":{"template":"toy_chr1_1001_1600","insertion":{"requested_forward_3prime_end_0based_exclusive":220,"requested_reverse_3prime_start_0based":300,"forward_extension_5prime":"GAATTC","reverse_extension_5prime":"CTCGAG","forward_window_start_0based":170,"forward_window_end_0based_exclusive":245,"reverse_window_start_0based":275,"reverse_window_end_0based_exclusive":360,"max_anchor_shift_bp":12},"forward":{"min_length":20,"max_length":30},"reverse":{"min_length":20,"max_length":30},"pair_constraints":{"require_roi_flanking":false},"min_amplicon_bp":120,"max_amplicon_bp":1200,"max_tm_delta_c":2.0,"max_pairs":50,"report_id":"insertion_demo_v1"}}'
cargo run --bin gentle_cli -- genomes list --catalog assets/genomes.json
cargo run --bin gentle_cli -- genomes list --catalog assets/helper_genomes.json
cargo run --bin gentle_cli -- genomes validate-catalog --catalog assets/genomes.json
cargo run --bin gentle_cli -- genomes update-ensembl-specs --catalog assets/genomes.json --output-catalog exports/genomes.updated.json
cargo run --bin gentle_cli -- genomes status "Human GRCh38 Ensembl 116" --catalog assets/genomes.json --cache-dir data/genomes
cargo run --bin gentle_cli -- genomes genes "Human GRCh38 Ensembl 116" --catalog assets/genomes.json --cache-dir data/genomes --filter "^TP53$" --biotype protein_coding --limit 20
cargo run --bin gentle_cli -- genomes prepare "Human GRCh38 Ensembl 116" --catalog assets/genomes.json --cache-dir data/genomes --timeout-secs 3600
cargo run --bin gentle_cli -- genomes remove-prepared "Human GRCh38 Ensembl 116" --catalog assets/genomes.json --cache-dir data/genomes
cargo run --bin gentle_cli -- genomes remove-catalog-entry "Human GRCh38 Ensembl 113" --catalog exports/genomes.custom.json
cargo run --bin gentle_cli -- genomes blast "Human GRCh38 Ensembl 116" ACGTACGTACGT --task blastn-short --max-hits 10 --options-json '{"thresholds":{"min_identity_percent":97.0,"min_query_coverage_percent":80.0}}' --catalog assets/genomes.json --cache-dir data/genomes
cargo run --bin gentle_cli -- genomes extract-region "Human GRCh38 Ensembl 116" 1 1000000 1001500 --output-id grch38_chr1_slice --annotation-scope core --catalog assets/genomes.json --cache-dir data/genomes
cargo run --bin gentle_cli -- genomes extract-gene "Human GRCh38 Ensembl 116" TP53 --occurrence 1 --output-id grch38_tp53 --catalog assets/genomes.json --cache-dir data/genomes
cargo run --bin gentle_cli -- tracks import-bed grch38_tp53 data/chipseq/peaks.bed.gz --name H3K27ac --min-score 10 --clear-existing
cargo run --bin gentle_cli -- tracks import-bigwig grch38_tp53 data/chipseq/signal.bw --name ATAC --min-score 0.2 --clear-existing
cargo run --bin gentle_cli -- tracks import-vcf grch38_tp53 data/variants/sample.vcf.gz --name Variants --min-score 20 --clear-existing
cargo run --bin gentle_cli -- helpers list
cargo run --bin gentle_cli -- helpers validate-catalog
cargo run --bin gentle_cli -- helpers update-ensembl-specs --catalog assets/helper_genomes.json --output-catalog exports/helper_genomes.updated.json
cargo run --bin gentle_cli -- helpers status "Plasmid pUC19 (online)"
cargo run --bin gentle_cli -- helpers prepare "Plasmid pUC19 (online)" --cache-dir data/helper_genomes --timeout-secs 600
cargo run --bin gentle_cli -- helpers remove-prepared "Plasmid pUC19 (online)" --cache-dir data/helper_genomes
cargo run --bin gentle_cli -- helpers genes "Plasmid pUC19 (online)" --filter bla --limit 20
cargo run --bin gentle_cli -- helpers blast "Plasmid pUC19 (online)" ACGTACGTACGT --task blastn-short --max-hits 10 --options-json '{"thresholds":{"min_identity_percent":95.0}}' --cache-dir data/helper_genomes
cargo run --bin gentle_cli -- cache inspect --references --cache-dir data/genomes
cargo run --bin gentle_cli -- cache clear derived-indexes-only --references --cache-dir data/genomes --prepared-id "Human GRCh38 Ensembl 116"
cargo run --bin gentle_cli -- cache clear selected-prepared --both --cache-dir data/genomes --cache-dir data/helper_genomes --prepared-path data/helper_genomes/localproject
cargo run --bin gentle_cli -- cache clear all-prepared-in-cache --both --include-orphans
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
cargo run --bin gentle_cli -- shell 'macros template-run gibson_two_fragment_overlap_preview --bind left_seq_id=gibson_left --bind right_seq_id=gibson_right --bind overlap_bp=20 --bind assembly_prefix=gibson_demo --bind output_id=gibson_demo_forward --validate-only'
cargo run --bin gentle_cli -- shell 'macros template-run gibson_two_fragment_overlap_preview --bind left_seq_id=gibson_left --bind right_seq_id=gibson_right --bind overlap_bp=20 --bind assembly_prefix=gibson_demo --bind output_id=gibson_demo_forward --transactional'
cargo run --bin gentle_cli -- shell 'macros run --transactional --file cloning_flow.gsh'
cargo run --bin gentle_cli -- shell 'set-param vcf_display_pass_only true'
cargo run --bin gentle_cli -- shell 'set-param vcf_display_required_info_keys ["AF","DP"]'
cargo run --bin gentle_cli -- shell 'set-param tfbs_display_min_llr_quantile 0.95'
cargo run --bin gentle_cli -- shell 'panels import-isoform grch38_tp53 assets/panels/tp53_isoforms_v1.json --panel-id tp53_isoforms_v1'
cargo run --bin gentle_cli -- shell 'panels inspect-isoform grch38_tp53 tp53_isoforms_v1'
cargo run --bin gentle_cli -- shell 'panels render-isoform-svg grch38_tp53 tp53_isoforms_v1 exports/tp53_isoform_architecture.svg'
cargo run --bin gentle_cli -- shell 'panels validate-isoform assets/panels/tp53_isoforms_v1.json --panel-id tp53_isoforms_v1'
cargo run --bin gentle_cli -- inspect-feature-expert grch38_tp53 isoform tp53_isoforms_v1
cargo run --bin gentle_cli -- render-feature-expert-svg grch38_tp53 isoform tp53_isoforms_v1 exports/tp53_isoform_architecture.svg
```

You can pass JSON from a file with `@file.json` or a bare existing file path.
When loading from file path, an initial shebang line (`#!...`) is ignored so
executable script files can embed JSON payloads directly.
`workflow` accepts both raw workflow payloads (`{"run_id":"...","ops":[...]}`)
and wrapped protocol example payloads (`{"workflow":{...}}`).

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
    - `render-dotplot-svg SEQ_ID DOTPLOT_ID OUTPUT.svg [--flex-track ID] [--display-threshold N] [--intensity-gain N]`
    - `dotplot render-svg SEQ_ID DOTPLOT_ID OUTPUT.svg [--flex-track ID] [--display-threshold N] [--intensity-gain N]`
    - `render-rna-svg SEQ_ID OUTPUT.svg`
    - `rna-info SEQ_ID`
    - `render-lineage-svg OUTPUT.svg`
    - `protocol-cartoon list`
    - `protocol-cartoon render-svg PROTOCOL_ID OUTPUT.svg`
    - `protocol-cartoon render-template-svg TEMPLATE.json OUTPUT.svg`
    - `protocol-cartoon template-validate TEMPLATE.json`
    - `protocol-cartoon render-with-bindings TEMPLATE.json BINDINGS.json OUTPUT.svg`
    - `protocol-cartoon template-export PROTOCOL_ID OUTPUT.json`
    - `render-pool-gel-svg IDS|'-' OUTPUT.svg [--ladders NAME[,NAME]] [--containers ID[,ID]] [--arrangement ARR_ID]`
    - `render-gel-svg IDS|'-' OUTPUT.svg [--ladders NAME[,NAME]] [--containers ID[,ID]] [--arrangement ARR_ID]`
    - `arrange-serial CONTAINER_IDS [--id ARR_ID] [--name TEXT] [--ladders NAME[,NAME]]`
    - `ladders list [--molecule dna|rna] [--filter TEXT]`
    - `ladders export OUTPUT.json [--molecule dna|rna] [--filter TEXT]`
    - `export-pool IDS OUTPUT.pool.gentle.json [HUMAN_ID]`
    - `export-run-bundle OUTPUT.run_bundle.json [--run-id RUN_ID]`
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
    - `genomes blast GENOME_ID QUERY_SEQUENCE [--max-hits N] [--task blastn-short|blastn] [--options-json JSON_OR_@FILE | --options-file PATH] [--catalog PATH] [--cache-dir PATH]`
    - `genomes blast-start GENOME_ID QUERY_SEQUENCE [--max-hits N] [--task blastn-short|blastn] [--options-json JSON_OR_@FILE | --options-file PATH] [--catalog PATH] [--cache-dir PATH]`
    - `genomes blast-status JOB_ID [--with-report]`
    - `genomes blast-cancel JOB_ID`
    - `genomes blast-list`
    - `genomes extract-region GENOME_ID CHR START END [--output-id ID] [--annotation-scope none|core|full] [--max-annotation-features N] [--include-genomic-annotation|--no-include-genomic-annotation] [--catalog PATH] [--cache-dir PATH]`
    - `genomes extract-gene GENOME_ID QUERY [--occurrence N] [--output-id ID] [--extract-mode gene|coding_with_promoter] [--promoter-upstream-bp N] [--annotation-scope none|core|full] [--max-annotation-features N] [--include-genomic-annotation|--no-include-genomic-annotation] [--catalog PATH] [--cache-dir PATH]`
    - `helpers list [--catalog PATH]`
    - `helpers validate-catalog [--catalog PATH]`
    - `helpers status HELPER_ID [--catalog PATH] [--cache-dir PATH]`
    - `helpers genes HELPER_ID [--catalog PATH] [--cache-dir PATH] [--filter REGEX] [--biotype NAME] [--limit N] [--offset N]`
    - `helpers prepare HELPER_ID [--catalog PATH] [--cache-dir PATH] [--timeout-secs N]`
    - `helpers blast HELPER_ID QUERY_SEQUENCE [--max-hits N] [--task blastn-short|blastn] [--options-json JSON_OR_@FILE | --options-file PATH] [--catalog PATH] [--cache-dir PATH]`
    - `helpers blast-start HELPER_ID QUERY_SEQUENCE [--max-hits N] [--task blastn-short|blastn] [--options-json JSON_OR_@FILE | --options-file PATH] [--catalog PATH] [--cache-dir PATH]`
    - `helpers blast-status JOB_ID [--with-report]`
    - `helpers blast-cancel JOB_ID`
    - `helpers blast-list`
    - `helpers extract-region HELPER_ID CHR START END [--output-id ID] [--annotation-scope none|core|full] [--max-annotation-features N] [--include-genomic-annotation|--no-include-genomic-annotation] [--catalog PATH] [--cache-dir PATH]`
    - `helpers extract-gene HELPER_ID QUERY [--occurrence N] [--output-id ID] [--extract-mode gene|coding_with_promoter] [--promoter-upstream-bp N] [--annotation-scope none|core|full] [--max-annotation-features N] [--include-genomic-annotation|--no-include-genomic-annotation] [--catalog PATH] [--cache-dir PATH]`
    - `tracks import-bed SEQ_ID PATH [--name NAME] [--min-score N] [--max-score N] [--clear-existing]`
    - `tracks import-bigwig SEQ_ID PATH [--name NAME] [--min-score N] [--max-score N] [--clear-existing]`
    - `tracks import-vcf SEQ_ID PATH [--name NAME] [--min-score N] [--max-score N] [--clear-existing]`
    - `macros run [--transactional] [--file PATH | SCRIPT_OR_@FILE]`
    - `macros instance-list`
    - `macros instance-show MACRO_INSTANCE_ID`
    - `macros template-list`
    - `macros template-show TEMPLATE_NAME`
    - `macros template-put TEMPLATE_NAME (--script SCRIPT_OR_@FILE|--file PATH) [--description TEXT] [--details-url URL] [--param NAME|NAME=DEFAULT ...] [--input-port PORT_ID:KIND[:one|many][:required|optional][:description]] [--output-port PORT_ID:KIND[:one|many][:required|optional][:description]]`
    - `macros template-delete TEMPLATE_NAME`
    - `macros template-import PATH`
    - `macros template-run TEMPLATE_NAME [--bind KEY=VALUE ...] [--transactional] [--validate-only]`
    - `routines list [--catalog PATH] [--family NAME] [--status NAME] [--tag TAG] [--query TEXT]`
    - `routines explain ROUTINE_ID [--catalog PATH]`
    - `routines compare ROUTINE_A ROUTINE_B [--catalog PATH]`
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
    - `features query SEQ_ID [--kind KIND] [--kind-not KIND] [--range START..END|--start N --end N] [--overlap|--within|--contains] [--strand any|forward|reverse] [--label TEXT] [--label-regex REGEX] [--qual KEY] [--qual-contains KEY=VALUE] [--qual-regex KEY=REGEX] [--min-len N] [--max-len N] [--limit N] [--offset N] [--sort feature_id|start|end|kind|length] [--desc] [--include-source] [--include-qualifiers]`
    - `primers design REQUEST_JSON_OR_@FILE [--backend auto|internal|primer3] [--primer3-exec PATH]`
    - `primers design-qpcr REQUEST_JSON_OR_@FILE [--backend auto|internal|primer3] [--primer3-exec PATH]`
    - `primers preflight [--backend auto|internal|primer3] [--primer3-exec PATH]`
    - `primers seed-from-feature SEQ_ID FEATURE_ID`
    - `primers seed-from-splicing SEQ_ID FEATURE_ID`
    - `primers list-reports`
    - `primers show-report REPORT_ID`
    - `primers export-report REPORT_ID OUTPUT.json`
    - `primers list-qpcr-reports`
    - `primers show-qpcr-report REPORT_ID`
    - `primers export-qpcr-report REPORT_ID OUTPUT.json`
    - `dotplot compute SEQ_ID [--reference-seq REF_SEQ_ID] [--start N] [--end N] [--ref-start N] [--ref-end N] [--mode self_forward|self_reverse_complement|pair_forward|pair_reverse_complement] [--word-size N] [--step N] [--max-mismatches N] [--tile-bp N] [--id DOTPLOT_ID]`
    - `dotplot list [SEQ_ID]`
    - `dotplot show DOTPLOT_ID`
    - `dotplot render-svg SEQ_ID DOTPLOT_ID OUTPUT.svg [--flex-track ID] [--display-threshold N] [--intensity-gain N]`
    - `transcripts derive SEQ_ID [--feature-id N ...] [--scope all_overlapping_both_strands|target_group_any_strand|all_overlapping_target_strand|target_group_target_strand] [--output-prefix PREFIX]`
    - `flex compute SEQ_ID [--start N] [--end N] [--model at_richness|at_skew] [--bin-bp N] [--smoothing-bp N] [--id TRACK_ID]`
    - `flex list [SEQ_ID]`
    - `flex show TRACK_ID`
    - `splicing-refs derive SEQ_ID START_0BASED END_0BASED [--seed-feature-id N] [--scope all_overlapping_both_strands|target_group_any_strand|all_overlapping_target_strand|target_group_target_strand] [--output-prefix PREFIX]`
    - `align compute QUERY_SEQ_ID TARGET_SEQ_ID [--query-start N] [--query-end N] [--target-start N] [--target-end N] [--mode global|local] [--match N] [--mismatch N] [--gap-open N] [--gap-extend N]`
    - `rna-reads interpret SEQ_ID FEATURE_ID INPUT.fa[.gz] [--report-id ID] [--report-mode full|seed_passed_only] [--checkpoint-path PATH] [--checkpoint-every-reads N] [--resume-from-checkpoint|--no-resume-from-checkpoint] [--profile nanopore_cdna_v1] [--format fasta] [--scope all_overlapping_both_strands|target_group_any_strand|all_overlapping_target_strand|target_group_target_strand] [--origin-mode single_gene|multi_gene_sparse] [--target-gene GENE_ID]... [--roi-seed-capture|--no-roi-seed-capture] [--kmer-len N] [--seed-stride-bp N] [--min-seed-hit-fraction F] [--min-weighted-seed-hit-fraction F] [--min-unique-matched-kmers N] [--min-chain-consistency-fraction F] [--max-median-transcript-gap F] [--min-confirmed-transitions N] [--min-transition-support-fraction F] [--cdna-poly-t-flip|--no-cdna-poly-t-flip] [--poly-t-prefix-min-bp N] [--align-band-bp N] [--align-min-identity F] [--max-secondary-mappings N]`
    - `rna-reads align-report REPORT_ID [--selection all|seed_passed|aligned] [--align-band-bp N] [--align-min-identity F] [--max-secondary-mappings N]`
    - `rna-reads list-reports [SEQ_ID]`
    - `rna-reads show-report REPORT_ID`
    - `rna-reads summarize-gene-support REPORT_ID --gene GENE_ID [--gene GENE_ID ...] [--record-indices i,j,k] [--complete-rule near|strict|exact] [--output PATH]`
    - `rna-reads inspect-alignments REPORT_ID [--selection all|seed_passed|aligned] [--limit N] [--effect-filter all_aligned|confirmed_only|disagreement_only|reassigned_only|no_phase1_only|selected_only] [--sort rank|identity|coverage|score] [--search TEXT] [--record-indices i,j,k] [--score-bin-variant all_scored|composite_seed_gate] [--score-bin-index N] [--score-bin-count M]`
    - `rna-reads export-report REPORT_ID OUTPUT.json`
    - `rna-reads export-hits-fasta REPORT_ID OUTPUT.fa [--selection all|seed_passed|aligned] [--record-indices i,j,k] [--subset-spec TEXT]`
    - `rna-reads export-sample-sheet OUTPUT.tsv [--seq-id ID] [--report-id ID]... [--append]`
    - `rna-reads export-paths-tsv REPORT_ID OUTPUT.tsv [--selection all|seed_passed|aligned] [--record-indices i,j,k] [--subset-spec TEXT]`
    - `rna-reads export-abundance-tsv REPORT_ID OUTPUT.tsv [--selection all|seed_passed|aligned] [--record-indices i,j,k] [--subset-spec TEXT]`
    - `rna-reads export-score-density-svg REPORT_ID OUTPUT.svg [--scale linear|log] [--variant all_scored|composite_seed_gate]`
    - `rna-reads export-alignments-tsv REPORT_ID OUTPUT.tsv [--selection all|seed_passed|aligned] [--limit N] [--record-indices i,j,k] [--subset-spec TEXT]`
    - `rna-reads export-alignment-dotplot-svg REPORT_ID OUTPUT.svg [--selection all|seed_passed|aligned] [--max-points N]`
    - `rna-reads align-report` re-ranks retained hits by alignment-aware
      retention rank after mapping refresh.
    - `rna-reads inspect-alignments` returns ranked aligned rows suitable for
      read-level inspection without mutating report payloads, and now echoes
      the formal structured subset spec (`effect_filter`, `sort_key`,
      `search`, `selected_record_indices`, `score_density_variant`,
      `score_bin_index`, `score_bin_count`) plus `subset_match_count` in the
      JSON payload so agent-driven inspection stays reproducible. Rows now also
      include full-length flags:
      - `full_length_exact`
      - `full_length_near`
      - `full_length_strict`
    - `rna-reads summarize-gene-support` returns a machine-readable
      `gentle.rna_read_gene_support_summary.v1` payload for one or more target
      genes from a saved aligned report:
      - base cohort = aligned retained rows, optionally narrowed by
        `--record-indices`
      - accepted target cohort = base rows whose best mapping resolves to one
        of the requested genes
      - cohort tables are emitted for `all_target`, `fragments`, and
        `complete`
      - `exon_pair_support` includes skipped ordered pairs like `1->3`
      - `direct_transition_support` includes neighboring exon steps only
      - `--output PATH` writes the same JSON returned on stdout
    - `rna-reads export-alignments-tsv` writes the same ranked alignment rows
      in TSV form for downstream filtering/sorting; `--record-indices`
      exports an exact saved-report subset and overrides coarse `--selection`;
      `--subset-spec` records the formal subset definition that produced it.
      TSV rows now include:
      - `target_coverage_fraction`
      - `full_length_exact`
      - `full_length_near`
      - `full_length_strict`
      - `full_length_class` (`exact` / `strict_end` / `near` / `partial`)
    - `rna-reads export-paths-tsv` and `rna-reads export-abundance-tsv` now
      accept the same `--record-indices` exact-subset override plus optional
      `--subset-spec` provenance.
    - `rna-reads export-alignment-dotplot-svg` emits a dotplot-style alignment
      scatter (coverage vs identity) with score-based point coloring.
    - `rna-reads export-hits-fasta` headers include seed metrics and exon-path
      annotations:
      - `--record-indices i,j,k` exports the exact saved-report subset and
        overrides coarse `--selection`
      - `--subset-spec TEXT` records the formal subset definition in the FASTA
        headers
      - `exon_path_tx=<transcript_id|none>`
      - `exon_path=<ordinal_path|none>` where `:` marks hash-confirmed adjacent
        exon transitions and `-` marks unconfirmed adjacency
      - `exon_transitions=<confirmed>/<total>`
      - `rc_applied=<true|false>` indicating automatic cDNA poly-T
        reverse-complement normalization was applied
      - `origin_class=<...>` with `origin_conf` / `strand_conf` summary values
        from deterministic origin classification
    - cDNA/direct-RNA seed-normalization semantics:
      - default (`--cdna-poly-t-flip`): reads with a T-rich 5' head are
        reverse-complement normalized before hashing/scoring
        - head detection tolerates small interruptions/mismatches in the
          poly-T stretch near the 5' end
      - direct RNA (`--no-cdna-poly-t-flip`): reads are scored as provided
      - `--poly-t-prefix-min-bp` controls the minimum T support threshold used
        by the 5' head detector for automatic cDNA flip (default `18`)
    - phase-1 seed-span policy:
      - all reads are hashed across their full span (no 3-window sampling)
      - `--seed-stride-bp` is the real phase-1 density knob: `1` means one
        ordered k-mer start per base, while larger values trade sensitivity for
        speed
    - sparse-origin options:
      - `--origin-mode` accepts `single_gene` (default) and
        `multi_gene_sparse`
      - `--target-gene` (repeatable) and `--roi-seed-capture` are persisted in
        the report payload for deterministic follow-up runs
      - `multi_gene_sparse` expands template indexing with local-annotation
        matches from `--target-gene`
      - `--roi-seed-capture` remains planned and is surfaced as a deterministic
        warning until implemented
    - `rna-reads list-reports` rows expose sparse-origin request provenance:
      `origin_mode`, `target_gene_count`, and `roi_seed_capture_enabled`
      - shell/CLI JSON output also includes `summary_rows[]` for quick
        human-readable triage (`mode`, `origin`, target count, ROI-capture
        flag, and read counters)
    - `rna-reads show-report` output includes a `summary` field with the same
      compact provenance framing, plus aligned/full-length percentages and
      compact read-length histogram summaries for:
      - all reads
      - seed-passed reads
      - aligned reads
      - full-length exact / near / strict subsets
    - tutorial reference (TP53 basis + multi-gene sparse mapping):
      `docs/tutorial/generated/chapters/12_tp53_multi_gene_sparse_mapping_online.md`
    - `rna-reads export-sample-sheet` includes sparse-origin provenance columns
      (`report_mode`, `origin_mode`, `target_gene_count`,
      `target_gene_ids_json`, `roi_seed_capture_enabled`) and
      `origin_class_counts_json` alongside exon/junction frequency JSON fields
    - report compaction and checkpoint options:
      - `--report-mode full` (default): persist retained top hits as ranked
      - `--report-mode seed_passed_only`: persist only retained hits that
        passed the seed gate (stream counters are unchanged)
      - `--checkpoint-path PATH --checkpoint-every-reads N`: persist
        deterministic resume snapshots during streaming
      - `--resume-from-checkpoint`: load checkpoint snapshot and continue from
        the saved processed-read index
    - phase-1 seed-pass gate policy:
      - `pass = raw_hit_fraction >= min_seed_hit_fraction AND weighted_hit_fraction >= min_weighted_seed_hit_fraction AND unique_matched_kmers >= min(min_unique_matched_kmers, tested_kmers) AND chain_consistency_fraction >= min_chain_consistency_fraction AND median_transcript_gap <= max_median_transcript_gap AND confirmed_transitions >= min_confirmed_transitions AND confirmed_transition_fraction >= min_transition_support_fraction`
      - `chain_consistency_fraction` is the fraction of matched seed observations
        explained by one coherent transcript-offset chain
      - `weighted_hit_fraction` is occurrence-normalized inside the scoped seed
        index: `sum(1/occurrence_count(seed_bits))/tested_kmers`
      - defaults:
        - `--kmer-len 10`
        - `--min-seed-hit-fraction 0.30`
        - `--min-weighted-seed-hit-fraction 0.05`
        - `--min-unique-matched-kmers 12`
        - `--min-chain-consistency-fraction 0.40`
        - `--max-median-transcript-gap 4.0`
        - `--min-confirmed-transitions 1`
        - `--min-transition-support-fraction 0.05`
    - RNA-read scope semantics:
      - `all_overlapping_both_strands`: index all overlapping transcripts on
        both strands
      - `target_group_any_strand`: index only target-group transcripts, both
        strands allowed
      - `all_overlapping_target_strand`: index all overlapping transcripts on
        target strand only
      - `target_group_target_strand`: index only target-group transcripts on
        target strand
      - strand-scoring note:
        both-strand modes score against the union of admitted `+/-` templates;
        target-strand modes exclude opposite-strand templates from scoring.
      - seed-index note:
        indexed seeds include exon-body and exon-exon transition k-mers for
        admitted transcripts.
    - Primer request payload notes (`primers design` / `primers design-qpcr`):
      - `forward`/`reverse`/`probe` side constraints support optional
        sequence filters:
        `non_annealing_5prime_tail`, `fixed_5prime`, `fixed_3prime`,
        `required_motifs[]`, `forbidden_motifs[]`, `locked_positions[]`.
      - `non_annealing_5prime_tail` is added to the oligo sequence while
        anneal `tm_c`/`gc_fraction`/`anneal_hits` remain computed on the
        template-binding segment only.
      - `pair_constraints` is optional and supports:
        `require_roi_flanking`, amplicon motif filters, and fixed amplicon
        start/end coordinates.
    - Primer ROI seed helper notes (`primers seed-from-feature` / `primers seed-from-splicing`):
      - returns non-mutating schema `gentle.primer_seed_request.v1`
      - includes `template`, source metadata, `roi_start_0based`,
        `roi_end_0based_exclusive`
      - includes ready-to-run operations:
        `operations.design_primer_pairs` (`DesignPrimerPairs`) and
        `operations.design_qpcr_assays` (`DesignQpcrAssays`)
    - Feature query helper notes (`features query`):
      - non-mutating structured result schema:
        `gentle.sequence_feature_query_result.v1`
      - deterministic filters over feature kind, range relation, strand, labels,
        qualifiers, and length
      - deterministic ordering (`feature_id|start|end|kind|length`) with
        `offset`/`limit` paging suitable for agent iteration
    - `panels import-isoform SEQ_ID PANEL_PATH [--panel-id ID] [--strict]`
    - `panels inspect-isoform SEQ_ID PANEL_ID`
    - `panels render-isoform-svg SEQ_ID PANEL_ID OUTPUT.svg`
    - `panels validate-isoform PANEL_PATH [--panel-id ID]`
    - `set-param NAME JSON_VALUE`
      - anchor-related examples:
        - `set-param require_verified_genome_anchor_for_extension true`
        - `set-param genome_anchor_prepared_fallback_policy "single_compatible|always_explicit|off"`
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

Isoform architecture panel workflow:

- canonical panel resource:
  - `assets/panels/tp53_isoforms_v1.json`
- operation route (JSON):
  - `ImportIsoformPanel` + `RenderIsoformArchitectureSvg`
- shared shell route:
  - `panels import-isoform ...`
  - `panels inspect-isoform ...`
  - `panels render-isoform-svg ...`
  - `panels validate-isoform ...`
- direct expert route:
  - `inspect-feature-expert SEQ_ID isoform PANEL_ID`
  - `render-feature-expert-svg SEQ_ID isoform PANEL_ID OUTPUT.svg`
  - same command family for splicing:
    - `inspect-feature-expert SEQ_ID splicing FEATURE_ID`
    - `render-feature-expert-svg SEQ_ID splicing FEATURE_ID OUTPUT.svg`

Pool exchange commands:

- `export-pool IDS OUTPUT.pool.gentle.json [HUMAN_ID]`
  - Exports explicit sequence IDs (`IDS` is comma-separated) with topology and
    overhang fields into a versioned JSON pool artifact.
  - Adds `human_id` at pool level and per-member `human_id`.
- `import-pool INPUT.pool.gentle.json [PREFIX]`
  - Imports pool members into current state; generated IDs are prefixed.
- `export-run-bundle OUTPUT.run_bundle.json [--run-id RUN_ID]`
  - Calls engine operation `ExportProcessRunBundle`.
  - Writes a deterministic process/audit bundle JSON
    (`gentle.process_run_bundle.v1`) with:
    - extracted operation inputs
    - chronological parameter overrides
    - selected operation log records
    - output summaries (created/changed ids, exported artifact paths)
  - Omitting `--run-id` exports all operation-log rows.

Rendering export commands:

- `render-svg SEQ_ID linear|circular OUTPUT.svg`
  - Calls engine operation `RenderSequenceSvg`.
- `render-dotplot-svg SEQ_ID DOTPLOT_ID OUTPUT.svg [--flex-track ID] [--display-threshold N] [--intensity-gain N]`
  - Calls engine operation `RenderDotplotSvg`.
  - `DOTPLOT_ID` must exist in stored dotplot payloads (`dotplot compute ...` / GUI compute).
  - `--flex-track` optionally overlays one stored flexibility track in the same SVG.
  - `--display-threshold` and `--intensity-gain` apply the same density/contrast controls as GUI dotplot display.
  - Alias: `dotplot render-svg SEQ_ID DOTPLOT_ID OUTPUT.svg [--flex-track ID] [--display-threshold N] [--intensity-gain N]`.
- `render-rna-svg SEQ_ID OUTPUT.svg`
  - Calls engine operation `RenderRnaStructureSvg`.
  - Accepts only single-stranded RNA (`molecule_type` of `RNA`/`ssRNA`).
  - Uses external `rnapkin` executable (set `GENTLE_RNAPKIN_BIN` to override executable path).
- `render-lineage-svg OUTPUT.svg`
  - Calls engine operation `RenderLineageSvg`.
- `protocol-cartoon list`
  - Lists deterministic protocol-cartoon IDs currently supported by the engine.
- `protocol-cartoon render-svg PROTOCOL_ID OUTPUT.svg`
  - Calls engine operation `RenderProtocolCartoonSvg`.
  - Built-in protocol IDs currently include `gibson.two_fragment`,
    `gibson.single_insert_dual_junction`, `pcr.assay.pair`,
    `pcr.assay.pair.no_product`, `pcr.assay.pair.with_tail`,
    `pcr.oe.substitution`, and `pcr.assay.qpcr`.
- `protocol-cartoon render-template-svg TEMPLATE.json OUTPUT.svg`
  - Calls engine operation `RenderProtocolCartoonTemplateSvg`.
  - Loads template JSON schema `gentle.protocol_cartoon_template.v1` and resolves
    sparse defaults before rendering.
- `protocol-cartoon template-validate TEMPLATE.json`
  - Calls engine operation `ValidateProtocolCartoonTemplate`.
  - Validates template JSON schema `gentle.protocol_cartoon_template.v1` and resolved cartoon semantics without writing SVG output.
- `protocol-cartoon render-with-bindings TEMPLATE.json BINDINGS.json OUTPUT.svg`
  - Calls engine operation `RenderProtocolCartoonTemplateWithBindingsSvg`.
  - Applies binding JSON schema `gentle.protocol_cartoon_template_bindings.v1` to the template before resolving and rendering.
- `protocol-cartoon template-export PROTOCOL_ID OUTPUT.json`
  - Calls engine operation `ExportProtocolCartoonTemplateJson`.
  - Exports canonical built-in template JSON (`gentle.protocol_cartoon_template.v1`)
    for the requested protocol ID.
  - Protocol-cartoon commands intentionally use one canonical namespace with no
    extra alias routes, which keeps scripted and AI-facing command selection
    unambiguous.
- `gibson preview PLAN_JSON_OR_@FILE [--output OUTPUT.json]`
  - Calls the shared non-mutating Gibson preview derivation path.
  - Accepts the single-insert destination-first subset of
    `gentle.gibson_assembly_plan.v1`.
  - Optional plan input:
    `validation_policy.desired_unique_restriction_site_enzyme_name`
    requests one new unique REBASE site on a terminal overlap when the current
    defined-site single-insert plan can support it.
  - The plan payload is resolved against the current project state, so the
    referenced `seq_id` values must already exist in memory or in the loaded
    state file.
  - Returns one `gentle.gibson_assembly_preview.v1` payload with:
    - resolved opening/junction overlaps
    - Gibson primer suggestions (`5' overlap + 3' priming`)
    - validation/advisory findings
    - structured suggested design adjustments for the current priming-window
      blocker when the overlap side already resolves
    - optional structured requested-unique-site outcome
    - protocol-cartoon bindings for the factual Gibson strip
- `gibson apply PLAN_JSON_OR_@FILE`
  - Calls engine operation `ApplyGibsonAssemblyPlan`.
  - Consumes the same plan JSON accepted by `gibson preview`.
  - Creates new project sequences for the Gibson outputs:
    - left insert primer
    - right insert primer
    - assembled product
  - Also records one serial arrangement for gel-oriented review:
    - original vector lane
    - ordered insert lane(s)
    - assembled product lane
    - recommended DNA ladders flanking the samples on export
  - Transfers destination and insert features onto the assembled product
    deterministically.
  - Partially consumed destination annotations are trimmed or projected when a
    truthful product-side rewrite is available.
  - Surviving MCS annotations are revalidated against actual restriction sites
    on the assembled product, including new sites introduced by the insert.
  - Returns the standard mutating `OpResult`, including created sequence IDs,
    warnings, messages, and the operation ID used for provenance/reopen flows.
- `cache inspect [--references|--helpers|--both] [--cache-dir PATH ...]`
  - Calls the shared non-mutating prepared-cache inspection path.
  - Returns one `gentle.prepared_cache_inspection.v1` payload with:
    - selected cache roots
    - manifest-backed prepared installs
    - orphaned remnants under those roots
    - per-entry artifact groups, byte totals, and file counts
  - If no `--cache-dir` is given, default roots are:
    - `--references`: `data/genomes`
    - `--helpers`: `data/helper_genomes`
    - `--both`: both roots above
- `cache clear blast-db-only|derived-indexes-only|selected-prepared|all-prepared-in-cache [--references|--helpers|--both] [--cache-dir PATH ...] [--prepared-id ID ...] [--prepared-path PATH ...] [--include-orphans]`
  - Calls the shared prepared-cache cleanup path.
  - Returns one `gentle.prepared_cache_cleanup.v1` payload with exact affected
    paths, removed byte totals, and per-entry cleanup results.
  - `blast-db-only` and `derived-indexes-only` apply only to manifest-backed
    prepared installs and require one or more `--prepared-id` or
    `--prepared-path` selectors.
  - `selected-prepared` removes only the selected prepared installs and may
    also remove orphaned remnants when `--include-orphans` is set.
  - `--prepared-path` is the precise selector when duplicate prepared ids exist
    across multiple selected cache roots.
  - `all-prepared-in-cache` clears all prepared installs under the selected
    roots and ignores `--prepared-id` / `--prepared-path`.
  - Cleanup is conservative by design:
    - only selected known cache roots are touched
    - catalog JSON and project state files remain unchanged
    - `.gentle_state.json`, MCP/runtime files, backdrop/runtime caches, and
      `target/` are out of scope
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

Conceptual/tutorial companion:

- `docs/agent_interfaces_tutorial.md` (role map + interface comparison:
  CLI/shared shell vs MCP vs Agent Assistant vs external coding agents).

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
- `genomes update-ensembl-specs [--catalog PATH] [--output-catalog PATH]`
  - Refreshes explicit pinned Ensembl URLs/specs for catalog rows that carry
    `ensembl_template` metadata.
  - Older pinned release rows stay in the catalog.
  - If a newer current species release exists, GENtle adds or refreshes the
    newest pinned row rather than overwriting the old one.
  - If the active catalog file is not writable, `--output-catalog PATH` is
    required so the updated catalog can be written as a copy.
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
  - Output now includes `binary_preflight` (`gentle.blast_external_binary_preflight.v1`)
    with `makeblastdb`/`blastn` found/missing/version/path diagnostics captured
    before prepare work starts.
- `genomes remove-prepared GENOME_ID [--catalog PATH] [--cache-dir PATH]`
  - Deletes only the prepared install directory for one genome from the selected cache.
  - Catalog JSON is left unchanged.
- `genomes remove-catalog-entry GENOME_ID [--catalog PATH] [--output-catalog PATH]`
  - Deletes only one row from the selected genome catalog JSON.
  - Prepared cache files are left unchanged until `remove-prepared` is run explicitly.
  - If the active catalog file is not writable, `--output-catalog PATH` is
    required so the edited catalog can be written as a copy.
- `genomes blast GENOME_ID QUERY_SEQUENCE [--max-hits N] [--task blastn-short|blastn] [--options-json JSON_OR_@FILE | --options-file PATH] [--catalog PATH] [--cache-dir PATH]`
  - Runs `blastn` against prepared genome cache/index.
  - `--task` defaults to `blastn-short`; accepted values: `blastn-short`, `blastn`.
  - `--options-json` / `--options-file` accept a JSON object that can override quick options and include thresholds (`max_evalue`, `min_identity_percent`, `min_query_coverage_percent`, `min_alignment_length_bp`, `min_bit_score`, `unique_best_hit`).
  - Output includes `binary_preflight` with explicit BLAST tool diagnostics.
- `genomes blast-start GENOME_ID QUERY_SEQUENCE [--max-hits N] [--task blastn-short|blastn] [--options-json JSON_OR_@FILE | --options-file PATH] [--catalog PATH] [--cache-dir PATH]`
  - Starts one async BLAST job and returns a stable `job_id`.
  - Jobs run through a bounded FIFO scheduler and may return initial state
    `queued` or `running` depending on current slot availability.
  - Start payload includes `binary_preflight` with explicit BLAST tool diagnostics.
- `genomes blast-status JOB_ID [--with-report]`
  - Polls async BLAST job status; `--with-report` includes final report payload when available.
  - Status includes scheduler metadata (`max_concurrent_jobs`, `running_jobs`,
    `queued_jobs`, and `queue_position` while queued).
- `genomes blast-cancel JOB_ID`
  - Requests cooperative cancellation for one async BLAST job.
- `genomes blast-list`
  - Lists known async genome-BLAST jobs in the current process.
- `genomes extract-region GENOME_ID CHR START END [--output-id ID] [--annotation-scope none|core|full] [--max-annotation-features N] [--include-genomic-annotation|--no-include-genomic-annotation] [--catalog PATH] [--cache-dir PATH]`
  - Runs engine `ExtractGenomeRegion`.
  - `--annotation-scope` selects projection policy (`none`, `core`, `full`).
  - Default scope is `core` when neither scope nor legacy flags are set.
  - `--max-annotation-features N` applies a deterministic safety cap
    (`0` disables the cap for explicit unrestricted transfer).
  - legacy `--include-genomic-annotation` maps to `--annotation-scope core`.
  - legacy `--no-include-genomic-annotation` maps to `--annotation-scope none`.
  - Result payload includes `genome_annotation_projection` telemetry.
- `genomes extract-gene GENOME_ID QUERY [--occurrence N] [--output-id ID] [--extract-mode gene|coding_with_promoter] [--promoter-upstream-bp N] [--annotation-scope none|core|full] [--max-annotation-features N] [--include-genomic-annotation|--no-include-genomic-annotation] [--catalog PATH] [--cache-dir PATH]`
  - Runs engine `ExtractGenomeGene`.
- `genomes extend-anchor SEQ_ID 5p|3p LENGTH_BP [--output-id ID] [--catalog PATH] [--cache-dir PATH] [--prepared-genome GENOME_ID]`
  - Runs engine `ExtendGenomeAnchor`.
  - Extends an already genome-anchored sequence in-silico on contextual `5'` or `3'`.
  - If exact anchor genome id is not prepared but one compatible assembly-family
    cache exists, extension auto-uses that cache and emits a warning.
  - If multiple compatible prepared caches exist, command fails and lists
    explicit options.
- `genomes verify-anchor SEQ_ID [--catalog PATH] [--cache-dir PATH] [--prepared-genome GENOME_ID]`
  - Runs engine `VerifyGenomeAnchor`.
  - Re-checks one anchored sequence against prepared reference sequence at
    recorded coordinates and records `anchor_verified` in provenance.

Helper convenience commands:

- `helpers list [--catalog PATH]`
  - Same behavior as `genomes list`, but defaults to `assets/helper_genomes.json`.
- `helpers validate-catalog [--catalog PATH]`
  - Same behavior as `genomes validate-catalog`, with helper-catalog default.
- `helpers update-ensembl-specs [--catalog PATH] [--output-catalog PATH]`
  - Same behavior as `genomes update-ensembl-specs`, with helper-catalog default.
- `helpers status HELPER_ID [--catalog PATH] [--cache-dir PATH]`
  - Same behavior as `genomes status`, with helper-catalog default
    (including length/mass metadata fields).
- `helpers genes HELPER_ID [--catalog PATH] [--cache-dir PATH] [--filter REGEX] [--biotype NAME] [--limit N] [--offset N]`
  - Same behavior as `genomes genes`, with helper-catalog default.
- `helpers prepare HELPER_ID [--catalog PATH] [--cache-dir PATH] [--timeout-secs N]`
  - Same behavior as `genomes prepare`, with helper-catalog default.
  - `--timeout-secs N`: optional prepare-job timebox.
- `helpers remove-prepared HELPER_ID [--catalog PATH] [--cache-dir PATH]`
  - Same behavior as `genomes remove-prepared`, with helper-catalog default.
- `cache inspect [--references|--helpers|--both] [--cache-dir PATH ...]`
  - Lists prepared installs and orphaned remnants under the selected cache
    roots with artifact-group byte totals.
  - Default roots are `data/genomes` for references and
    `data/helper_genomes` for helpers.
- `cache clear blast-db-only|derived-indexes-only|selected-prepared|all-prepared-in-cache [--references|--helpers|--both] [--cache-dir PATH ...] [--prepared-id ID ...] [--prepared-path PATH ...] [--include-orphans]`
  - Conservatively deletes derived cache artifacts or prepared installs inside
    the selected roots only.
  - Partial modes keep cached FASTA/annotation sources and manifests intact so
    reindex-from-cached-files remains possible.
  - Catalog JSON, `.gentle_state.json`, MCP/runtime files, backdrop/runtime
    caches, and `target/` are not treated as cache.
- `helpers remove-catalog-entry HELPER_ID [--catalog PATH] [--output-catalog PATH]`
  - Same behavior as `genomes remove-catalog-entry`, with helper-catalog default.
- `helpers blast HELPER_ID QUERY_SEQUENCE [--max-hits N] [--task blastn-short|blastn] [--options-json JSON_OR_@FILE | --options-file PATH] [--catalog PATH] [--cache-dir PATH]`
  - Same behavior as `genomes blast`, with helper-catalog default.
- `helpers blast-start HELPER_ID QUERY_SEQUENCE [--max-hits N] [--task blastn-short|blastn] [--options-json JSON_OR_@FILE | --options-file PATH] [--catalog PATH] [--cache-dir PATH]`
  - Same behavior as `genomes blast-start`, with helper-catalog default.
- `helpers blast-status JOB_ID [--with-report]`
  - Same behavior as `genomes blast-status`, scoped to helper jobs.
- `helpers blast-cancel JOB_ID`
  - Same behavior as `genomes blast-cancel`, scoped to helper jobs.
- `helpers blast-list`
  - Same behavior as `genomes blast-list`, scoped to helper jobs.
- `helpers extract-region HELPER_ID CHR START END [--output-id ID] [--annotation-scope none|core|full] [--max-annotation-features N] [--include-genomic-annotation|--no-include-genomic-annotation] [--catalog PATH] [--cache-dir PATH]`
  - Same behavior as `genomes extract-region`, with helper-catalog default.
  - For helper IDs containing `pUC18`/`pUC19`, GENtle auto-attaches a
    canonical MCS `misc_feature` when no MCS annotation is present in source
    annotation and exactly one canonical MCS motif is found.
  - MCS features expose `mcs_expected_sites` with REBASE-normalized enzyme names
    when recognizable from source/fallback annotation text.
- `helpers extract-gene HELPER_ID QUERY [--occurrence N] [--output-id ID] [--extract-mode gene|coding_with_promoter] [--promoter-upstream-bp N] [--annotation-scope none|core|full] [--max-annotation-features N] [--include-genomic-annotation|--no-include-genomic-annotation] [--catalog PATH] [--cache-dir PATH]`
  - Same behavior as `genomes extract-gene`, with helper-catalog default.
  - pUC18/pUC19 helper extractions apply the same automatic MCS fallback
    annotation behavior when applicable (non-unique motif matches are warned and
    skipped).
- `helpers extend-anchor SEQ_ID 5p|3p LENGTH_BP [--output-id ID] [--catalog PATH] [--cache-dir PATH] [--prepared-genome GENOME_ID]`
  - Same behavior as `genomes extend-anchor`, with helper-catalog default.
- `helpers verify-anchor SEQ_ID [--catalog PATH] [--cache-dir PATH] [--prepared-genome GENOME_ID]`
  - Same behavior as `genomes verify-anchor`, with helper-catalog default.

Workflow macro commands (`gentle_cli shell 'macros ...'`):

- `macros run [--transactional] [--file PATH | SCRIPT_OR_@FILE]`
  - Executes semicolon/newline-separated shell statements.
  - Existing file paths are auto-loaded even without `@` (shebang-friendly).
  - Supports transactional rollback (`--transactional`) when any statement fails.
  - Designed for full cloning workflows through `op ...` and `workflow ...`
    statements (Digest/Ligation/PCR/ExtractRegion/container ops, etc.).
  - All runs persist a lineage macro-instance record:
    - success: status `ok`
    - failure: status `failed` (or `cancelled` when cancellation-like error text is detected)
  - Successful runs return `macro_instance_id`; failed runs include
    `macro_instance_id=...` in the error message.
- `macros instance-list`
  - Lists recorded macro-instance lineage rows.
  - Response schema: `gentle.lineage_macro_instances.v1`.
- `macros instance-show MACRO_INSTANCE_ID`
  - Shows one recorded macro-instance lineage row.
  - Response schema: `gentle.lineage_macro_instance.v1`.
- `macros template-list`
  - Lists persisted workflow macro templates.
- `macros template-show TEMPLATE_NAME`
  - Shows one persisted workflow template definition.
- `macros template-put TEMPLATE_NAME (--script SCRIPT_OR_@FILE|--file PATH) [--description TEXT] [--details-url URL] [--param NAME|NAME=DEFAULT ...] [--input-port PORT_ID:KIND[:one|many][:required|optional][:description]] [--output-port PORT_ID:KIND[:one|many][:required|optional][:description]]`
  - Creates/updates a named workflow macro template in project metadata.
  - Placeholders in script use `${param_name}` and must be declared via `--param`.
  - Optional `--details-url URL` records external protocol/reference details for
    catalog display.
  - Optional typed port contracts (`--input-port` / `--output-port`) are stored
    in template metadata and used as first-class preflight contract source.
- `macros template-delete TEMPLATE_NAME`
  - Deletes one persisted workflow template.
- `macros template-import PATH`
  - Imports workflow macro templates from:
    - one pack JSON file (`gentle.cloning_patterns.v1`)
    - one single-template JSON file (`gentle.cloning_pattern_template.v1`)
    - one directory tree (recursive `*.json` import)
  - If one template fails validation, no imported template changes are kept.
- `macros template-run TEMPLATE_NAME [--bind KEY=VALUE ...] [--transactional] [--validate-only]`
  - Expands a named template with provided bindings/defaults, then executes it as
    a workflow macro script.
  - Executes typed preflight checks before mutation using template port
    contracts when present, otherwise routine catalog mapping.
  - Preflight now includes richer semantics:
    - cross-port alias/collision checks,
    - input sequence vs input container consistency checks,
    - sequence-anchor semantic checks against the bound input sequence when unambiguous,
    - Gibson-family overlap checks (adjacent suffix/prefix validation against
      configured overlap length on Gibson routines),
    - restriction-family digest checks (enzyme-name resolution, duplicate-enzyme
      misuse, recognition-site presence across bound input sequences, and
      digest parameter sanity for fragment indices/extract range).
  - `--validate-only` runs preflight only and never mutates state.
  - Mutating runs now always record a macro-instance lineage row:
    - success: `ok`
    - preflight/execute failure: `failed`/`cancelled`
  - Successful runs return `macro_instance_id`; failed runs include
    `macro_instance_id=...` in the error message.

Typed routine catalog command (`gentle_cli routines ...` or `gentle_cli shell 'routines ...'`):

- `routines list [--catalog PATH] [--family NAME] [--status NAME] [--tag TAG] [--query TEXT]`
  - Lists typed cloning routines from catalog JSON (`gentle.cloning_routines.v1`).
  - `--family`, `--status`, `--tag`: exact case-insensitive filters.
  - `--query`: case-insensitive substring match across id/title/family/status/template/tags/summary plus explainability metadata fields.
  - Default catalog path: `assets/cloning_routines.json`.
- `routines explain ROUTINE_ID [--catalog PATH]`
  - Returns structured explainability payload for one routine.
  - Response schema: `gentle.cloning_routine_explain.v1`.
  - Includes resolved alternatives plus purpose/mechanism/requirements, contraindications, disambiguation questions, and failure modes.
- `routines compare ROUTINE_A ROUTINE_B [--catalog PATH]`
  - Returns deterministic side-by-side comparison payload for two routines.
  - Response schema: `gentle.cloning_routine_compare.v1`.
  - Includes shared/unique vocabulary tags, difference-matrix rows, and merged disambiguation questions.
  - Also includes planning-aware estimate rows:
    `estimated_time_hours`, `estimated_cost`, `local_fit_score`,
    `composite_meta_score`.

Planning meta-layer commands (`gentle_cli planning ...` or `gentle_cli shell 'planning ...'`):

- Planning schemas:
  - `gentle.planning_profile.v1`
  - `gentle.planning_objective.v1`
  - `gentle.planning_estimate.v1`
  - `gentle.planning_suggestion.v1`
  - `gentle.planning_sync_status.v1`
- Effective profile merge precedence:
  - `global_profile -> confirmed_agent_overlay -> project_override`
- Schema compatibility:
  - payloads declaring a non-matching planning schema id are rejected
    (not auto-upgraded silently).
- `planning profile show [--scope global|project_override|confirmed_agent_overlay|effective]`
  - Shows one profile scope or merged effective profile.
- `planning profile set JSON_OR_@FILE [--scope global|project_override|confirmed_agent_overlay]`
  - Sets/replaces one profile scope.
- `planning profile clear [--scope global|project_override|confirmed_agent_overlay]`
  - Clears one profile scope.
- `planning objective show`
  - Shows current objective weights/guardrails.
- `planning objective set JSON_OR_@FILE`
  - Sets/replaces current planning objective.
- `planning objective clear`
  - Clears objective and falls back to engine defaults.
- `planning suggestions list [--status pending|accepted|rejected]`
  - Lists planning sync suggestions.
- `planning suggestions accept SUGGESTION_ID`
  - Accepts a pending suggestion and applies its patches.
- `planning suggestions reject SUGGESTION_ID [--reason TEXT]`
  - Rejects a pending suggestion.
- `planning sync status`
  - Shows sync metadata and pending count.
- `planning sync pull JSON_OR_@FILE [--source ID] [--confidence N] [--snapshot-id ID]`
  - Registers a pending inbound advisory suggestion.
- `planning sync push JSON_OR_@FILE [--source ID] [--confidence N] [--snapshot-id ID]`
  - Registers a pending outbound advisory suggestion.
- Sync payload shape (`pull`/`push`):
  - optional `profile_patch` (`gentle.planning_profile.v1`)
  - optional `objective_patch` (`gentle.planning_objective.v1`)
  - optional `message`

Planning estimate rule (v1 purchasing simplification):

- Missing required local material adds procurement delay using
  `procurement_business_days_default` (default `10`) unless overridden on the
  item (`inventory.<class>.procurement_business_days`).
- Delay is applied once per missing required material class (deduplicated).
- Business-day model in v1: Monday-Friday only (no holiday calendar).
- Business-day delays are converted to `estimated_time_hours` with a
  deterministic weekend-aware factor (`24h * 7/5` per business day).

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
  - includes `gibson.two_fragment_overlap_preview`
    -> `assets/cloning_patterns_catalog/gibson/overlap_assembly/gibson_two_fragment_overlap_preview.json`
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

Derive splicing references from a sequence window (DNA window + mRNA isoforms +
exon-reference):

```json
{"DeriveSplicingReferences":{"seq_id":"tp73_window","span_start_0based":1200,"span_end_0based":2600,"seed_feature_id":null,"scope":"target_group_target_strand","output_prefix":"tp73_refs"}}
```

Pairwise sequence alignment (global/local) with structured result payload:

```json
{"AlignSequences":{"query_seq_id":"tp73_refs_mrna_NM_001204186","target_seq_id":"tp73_refs_exon_reference","mode":"local","match_score":2,"mismatch_score":-3,"gap_open":-5,"gap_extend":-1}}
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

Enable strict anchor-verification policy for genome-anchor extension:

```json
{"SetParameter":{"name":"require_verified_genome_anchor_for_extension","value":true}}
```

Set prepared-genome fallback policy for anchor extension/verification:

```json
{"SetParameter":{"name":"genome_anchor_prepared_fallback_policy","value":"single_compatible"}}
{"SetParameter":{"name":"genome_anchor_prepared_fallback_policy","value":"always_explicit"}}
{"SetParameter":{"name":"genome_anchor_prepared_fallback_policy","value":"off"}}
```

Set feature-details font size used in the feature tree/details panel (valid range `8.0..24.0`):

```json
{"SetParameter":{"name":"feature_details_font_size","value":10.5}}
```

Set sequence text-panel maximum length (`200000` default, `0` means unlimited):

```json
{"SetParameter":{"name":"sequence_panel_max_text_length_bp","value":200000}}
```

Set regulatory-overlay max linear view span threshold (`50000` recommended for anchored genome maps):

```json
{"SetParameter":{"name":"regulatory_feature_max_view_span_bp","value":50000}}
```

Set primer-design backend controls (internal baseline + optional Primer3 backend):

```json
{"SetParameter":{"name":"primer_design_backend","value":"auto"}}
{"SetParameter":{"name":"primer3_executable","value":"primer3_core"}}
```

Set adaptive linear DNA-letter routing parameters (shared GUI/runtime semantics):

```json
{"SetParameter":{"name":"linear_sequence_letter_layout_mode","value":"auto"}}
{"SetParameter":{"name":"linear_sequence_helical_letters_enabled","value":true}}
{"SetParameter":{"name":"linear_sequence_helical_phase_offset_bp","value":3}}
{"SetParameter":{"name":"linear_show_reverse_strand_bases","value":true}}
{"SetParameter":{"name":"linear_helical_parallel_strands","value":true}}
{"SetParameter":{"name":"reverse_strand_visual_opacity","value":0.55}}
```

Supported `linear_sequence_letter_layout_mode` values:

- `auto`, `adaptive`, `auto_adaptive`
- `standard`, `standard_linear`
- `helical`, `continuous_helical`
- `condensed_10_row`, `condensed-10-row`, `condensed`

Compatibility notes:

- Legacy fixed-threshold parameters are still accepted for compatibility but are
  deprecated no-ops under adaptive routing:
  - `linear_sequence_base_text_max_view_span_bp`
  - `linear_sequence_helical_max_view_span_bp`
  - `linear_sequence_condensed_max_view_span_bp`
- `linear_sequence_helical_letters_enabled` applies to auto mode only
  (`linear_sequence_letter_layout_mode = auto*`).
- reverse-strand visibility aliases map to the same control:
  - `linear_show_double_strand_bases`
  - `linear_show_reverse_strand_bases`
- `linear_helical_parallel_strands=true` keeps forward/reverse helical slant
  in parallel; `false` uses mirrored slant.
- `reverse_strand_visual_opacity` (aliases:
  `linear_reverse_strand_visual_opacity`,
  `linear_reverse_strand_letter_opacity`) controls reverse-strand letter
  emphasis across linear map and sequence panel (`0.2..1.0`).

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

- Gibson apply operations export as the same dedicated `Gibson cloning` hub
  shown in the GUI, instead of a raw parent-to-each-output multi-edge graph.

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
{"ExtractGenomeRegion":{"genome_id":"Human GRCh38 Ensembl 116","chromosome":"1","start_1based":1000000,"end_1based":1001500,"output_id":"grch38_chr1_1000000_1001500","annotation_scope":"core","catalog_path":"assets/genomes.json","cache_dir":"data/genomes"}}
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
- BLAST progress:
  - BLAST+ does not provide a native percent-progress CLI flag for `blastn`.
  - GENtle surfaces deterministic progress at orchestration level (query counts,
    running elapsed time, and invocation context in adapters that support live status).
  - Async BLAST jobs additionally expose queue/scheduler state:
    `queued|running|completed|failed|cancelled` and scheduler counters.
  - Scheduler concurrency defaults to available CPU parallelism and can be
    overridden with `GENTLE_BLAST_ASYNC_MAX_CONCURRENT` (`1..256`).
- BLAST options layering (engine contract):
  - built-in defaults (`task=blastn-short`, `max_hits=25`)
  - optional defaults file (`assets/blast_defaults.json` or path set by parameter)
  - optional project override metadata (`set-param blast_options_override ...`)
  - quick flags (`--max-hits`, `--task`)
  - per-command JSON override (`--options-json` / `--options-file`)
  - project defaults file path can be set via
    `set-param blast_options_defaults_path '"path/to/blast_defaults.json"'`
- BLAST executable overrides:
  - `GENTLE_MAKEBLASTDB_BIN` (default: `makeblastdb`)
  - `GENTLE_BLASTN_BIN` (default: `blastn`)
- When BLAST hits are imported as features (`ImportBlastHitsTrack` /
  shell `genomes|helpers blast-track`), operation payload now includes
  `blast_provenance` with invocation details (`blastn_executable`, `blast_db_prefix`,
  raw `command` args, and compact `command_line`) so sequence history can
  trace how hits were generated.
- `ImportGenomeBedTrack` expects `seq_id` to be a sequence created by
  `ExtractGenomeRegion`, `ExtractGenomeGene`, or `ExtendGenomeAnchor`
  (genome-anchored provenance).
- `ImportGenomeBigWigTrack` expects the same genome-anchored `seq_id`.
- `ImportGenomeVcfTrack` expects the same genome-anchored `seq_id`.
- BED import accepts local `.bed` and `.bed.gz` files.
- Concatenated gzip members are accepted for `.bed.gz` track input.
- BigWig import accepts local `.bw` and `.bigWig` files and uses
  `bigWigToBedGraph` (override with `GENTLE_BIGWIG_TO_BEDGRAPH_BIN`).
- VCF import accepts local `.vcf` and `.vcf.gz` files.
- Concatenated gzip members are accepted for `.vcf.gz` track input.
- For VCF import, `min_score` / `max_score` filter on VCF `QUAL`.
- `ExtractGenomeRegion`, `ExtractGenomeGene`, and `ExtendGenomeAnchor` append extraction provenance
  records into `ProjectState.metadata["provenance"]["genome_extractions"]`
  (genome id, coordinates/query, source descriptors, and checksums when present).
- If `catalog_path` is omitted, engine default catalog is `assets/genomes.json`.
- Bundled `assets/genomes.json` currently includes Human GRCh38 (Ensembl 113 and 116),
  Mouse GRCm39 Ensembl 116, Rat GRCr8 Ensembl 115, Danio rerio GRCz11
  Ensembl 115, Pan troglodytes Pan_tro_3.0 Ensembl 115, Canis lupus familiaris
  ROS_Cfam_1.0 Ensembl 115, Drosophila melanogaster BDGP6.54 Ensembl Metazoa 62,
  Caenorhabditis elegans WBcel235 Ensembl 115, Saccharomyces cerevisiae S288c
  (Ensembl 113 and 115), and `LocalProject` (backed by
  `test_files/fixtures/genomes/AB011549.2.fa` +
  `test_files/fixtures/genomes/AB011549.2.gb`).
- `cache_dir` is optional. If omitted, catalog/default cache settings are used.
- `PrepareGenome` now validates that gene-bearing contigs parsed from the
  prepared annotation are present in the prepared FASTA index; truncated or
  mismatched installs fail during preparation instead of only later during
  extraction.
- Catalog rows with `ensembl_template` metadata can be refreshed on demand via
  `update-ensembl-specs`, which rewrites explicit pinned URLs without preparing
  any genomes.
- `chromosome` accepts exact contig names and also tolerates `chr` prefix
  differences (`1` vs `chr1`).
- Missing-contig extraction errors now also report the prepared `sequence.fa`
  and `sequence.fa.fai` paths for cache debugging.
- For `ExtractGenomeGene`, `occurrence` is 1-based among matching records.
- `ExtractGenomeGene.extract_mode=coding_with_promoter` resolves the CDS span
  from transcript `CDS` annotation and applies `promoter_upstream_bp`
  strand-aware on the gene's 5' side.
- For `ExtendGenomeAnchor`, `side` is contextual to anchor strand.
  On anchor strand `-`, `5'` increases physical genomic position.
- For anchor-extension reads, genome ids can resolve through assembly-family
  compatibility (for example `GRCh38.p14` using a prepared `Human GRCh38 ...`
  cache) when unique; ambiguous matches are rejected with options.

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
