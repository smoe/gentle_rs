# JavaScript Shell

<!-- manual-counterparts:start -->
No direct GUI counterpart: This chapter documents an automation interface, example or command-only workflow; no dedicated GUI chapter is paired with it.
<!-- manual-counterparts:end -->

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
14. `list_reference_catalog_entries(catalog_path, filter)`
    - Lists structured reference-catalog entries, including typed metadata.
    - `filter` is optional and matches the same search surface as `genomes list --filter`.
15. `list_helper_catalog_entries(catalog_path, filter)`
    - Lists structured helper-catalog entries, including optional normalized
      `interpretation` records.
    - `filter` is optional and matches the same search surface as `helpers list --filter`.
16. `list_helper_semantics_vocabulary(vocabulary_path, filter)`
    - Lists resolved helper semantics vocabulary terms, aliases, descriptions, sources, and routine hints.
    - `vocabulary_path` is optional (`null`/`""` uses default built-in/system/user/project discovery).
    - `filter` is optional and matches axis, value, aliases, descriptions, sources, and routine hints.
17. `list_host_profile_catalog_entries(catalog_path, filter)`
    - Lists structured host-profile catalog entries for construct-reasoning host/strain lookup.
    - `filter` is optional and matches ids, aliases, species, strain, genotype tags, phenotype tags, and notes.
18. `list_ensembl_installable_genomes(collection, filter)`
    - Lists Ensembl species directories that currently appear installable because both FASTA and GTF listings are present.
    - `collection` is optional (`all`, `vertebrates`, `metazoa`).
19. `list_construct_reasoning_graphs(state, seq_id)`
    - Lists stored construct-reasoning graphs plus compact shared-shell summary rows.
    - `seq_id` is optional and limits rows to one active sequence.
19. `show_construct_reasoning_graph(state, graph_id)`
    - Returns one stored construct-reasoning graph plus the same compact summary block exposed by `construct-reasoning show-graph`.
20. `list_construct_reasoning_inspection_actions(state, graph_id, options)`
    - Lists portable recommended inspection actions for a construct-reasoning graph through `construct-reasoning list-inspection-actions`.
    - `options.fact_id`, `options.annotation_id`, and `options.summary_id` are optional filters.
20. `run_construct_reasoning_inspection_action(state, graph_id, action_id, options)`
    - Runs a portable recommended inspection action through `construct-reasoning run-inspection-action` and returns the updated state plus shared-shell mutation output.
    - Dotplot options include `word_size`, `step_bp`, `max_mismatches`, `tile_bp`, `dotplot_id`, and `render_svg_path`.
20. `set_construct_reasoning_annotation_status(state, graph_id, annotation_id, status)`
    - Updates one stored annotation-candidate review status and returns the updated state plus shared-shell mutation output.
21. `is_reference_genome_prepared(genome_id, catalog_path, cache_dir)`
    - Checks whether a genome cache is prepared and indexed.
21. `list_reference_genome_genes(genome_id, catalog_path, cache_dir)`
    - Lists indexed genes from prepared genome cache.
22. `prepare_genome(state, genome_id, catalog_path, cache_dir)`
    - Convenience wrapper around engine `PrepareGenome`.
23. `extract_genome_region(state, genome_id, chromosome, start_1based, end_1based, output_id, catalog_path, cache_dir)`
    - Convenience wrapper around engine `ExtractGenomeRegion`.
24. `extract_genome_gene(state, genome_id, gene_query, occurrence, output_id, catalog_path, cache_dir, annotation_scope, max_annotation_features, extract_mode, promoter_upstream_bp)`
    - Convenience wrapper around engine `ExtractGenomeGene`.
25. `blast_reference_genome(genome_id, query_sequence, max_hits, task, catalog_path, cache_dir, options_json)`
    - Runs BLAST (`blastn`/`blastn-short`) against a prepared reference genome.
    - `options_json` is optional and may override any quick option (`max_hits`, `task`) plus threshold fields.
26. `blast_helper_genome(helper_id, query_sequence, max_hits, task, catalog_path, cache_dir, options_json)`
    - Same as `blast_reference_genome`, but defaults to helper catalog context.
27. `set_parameter(state, name, value)`
    - Convenience wrapper for engine `SetParameter`.
28. `set_vcf_display_filter(state, options)`
    - Convenience wrapper that updates one or more VCF display-filter parameters
      via `SetParameter` in one call.
29. `list_agent_systems(catalog_path)`
    - Lists configured agent systems from catalog JSON.
    - `catalog_path` is optional (`null`/`""` uses `assets/agent_systems.json`).
30. `ask_agent_system(state, system_id, prompt, options)`
    - Invokes one configured agent system through shared shell execution.
    - Returns `{ state, state_changed, output }`.
    - `state` may be `null` to use an empty/default project state.
    - `options` supports:
      - `catalog_path`
      - `allow_auto_exec` / `execute_all`
      - `execute_indices` (1-based)
      - `include_state_summary` (default `true`)
31. `plan_agent_system(state, system_id, prompt, options)`
    - Compiles prose into a typed `gentle.agent_plan_result.v1` plan through
      the shared `agents plan` route.
    - `state` may be `null`; options support `catalog_path`,
      `include_state_summary`, `max_candidates`, and
      `allow_mutating_candidates`.
32. `execute_agent_plan(state, plan, candidate_id, options)`
    - Executes one stored plan candidate through shared `agents execute-plan`.
    - `plan` may be a JSON string or JS object.
    - Returns `{ state, state_changed, output }`; options support `confirm`.
33. `render_dotplot_svg(state, seq_id, dotplot_id, output_svg, options)`
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
- `docs/examples/generated/tp73_isoform_protein_gel_offline.md`
- `docs/examples/generated/tp73_isoform_protein_2d_gel_offline.md`

<!-- manual-page-nav:start -->

---

[CLI Manual](../cli.md) | [Previous: ClawBio and OpenClaw](clawbio.md) | [Next: Lua Shell](lua.md)
<!-- manual-page-nav:end -->
