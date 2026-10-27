# Lua Shell

<!-- manual-counterparts:start -->
No direct GUI counterpart: This chapter documents an automation interface, example or command-only workflow; no dedicated GUI chapter is paired with it.
<!-- manual-counterparts:end -->

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
14. `list_reference_catalog_entries([catalog_path], [filter])`
    - Lists structured reference-catalog entries, including typed metadata.
15. `list_helper_catalog_entries([catalog_path], [filter])`
    - Lists structured helper-catalog entries, including optional normalized
      `interpretation` records.
16. `list_helper_semantics_vocabulary([vocabulary_path], [filter])`
    - Lists resolved helper semantics vocabulary terms.
17. `list_host_profile_catalog_entries([catalog_path], [filter])`
    - Lists structured host-profile catalog entries for construct-reasoning host/strain lookup.
18. `list_ensembl_installable_genomes([collection], [filter])`
    - Lists current Ensembl candidates where both FASTA and GTF species listings exist.
18. `list_construct_reasoning_graphs(project, [seq_id])`
    - Lists stored construct-reasoning graphs plus compact shared-shell summary rows.
19. `show_construct_reasoning_graph(project, graph_id)`
    - Returns one stored construct-reasoning graph plus the same compact summary block exposed by `construct-reasoning show-graph`.
20. `list_construct_reasoning_inspection_actions(project, graph_id, [fact_id], [annotation_id], [summary_id])`
    - Lists portable recommended inspection actions for a construct-reasoning graph through `construct-reasoning list-inspection-actions`.
20. `run_construct_reasoning_inspection_action(project, graph_id, action_id, [word_size], [step_bp], [max_mismatches], [tile_bp], [dotplot_id], [render_svg_path])`
    - Runs a portable recommended inspection action through `construct-reasoning run-inspection-action` and returns the updated project state plus shared-shell mutation output.
20. `set_construct_reasoning_annotation_status(project, graph_id, annotation_id, status)`
    - Updates one stored annotation-candidate review status and returns the updated project state plus shared-shell mutation output.
21. `is_reference_genome_prepared(genome_id, catalog_path, cache_dir)`
    - Checks whether a genome cache is prepared and indexed.
21. `list_reference_genome_genes(genome_id, catalog_path, cache_dir)`
    - Lists indexed genes from prepared genome cache.
22. `prepare_genome(project, genome_id, catalog_path, cache_dir)`
    - Convenience wrapper around engine `PrepareGenome`.
23. `extract_genome_region(project, genome_id, chromosome, start_1based, end_1based, output_id, catalog_path, cache_dir)`
    - Convenience wrapper around engine `ExtractGenomeRegion`.
24. `extract_genome_gene(project, genome_id, gene_query, occurrence, output_id, catalog_path, cache_dir, annotation_scope, max_annotation_features, extract_mode, promoter_upstream_bp)`
    - Convenience wrapper around engine `ExtractGenomeGene`.
25. `blast_reference_genome(genome_id, query_sequence, [max_hits], [task], [catalog_path], [cache_dir], [options_json])`
    - Runs BLAST (`blastn`/`blastn-short`) against a prepared reference genome.
    - `options_json` is optional and may override any quick option (`max_hits`, `task`) plus threshold fields.
26. `blast_helper_genome(helper_id, query_sequence, [max_hits], [task], [catalog_path], [cache_dir], [options_json])`
    - Same as `blast_reference_genome`, but defaults to helper catalog context.
27. `set_parameter(project, name, value)`
    - Convenience wrapper for engine `SetParameter`.
28. `set_vcf_display_filter(project, opts)`
    - Convenience wrapper that updates one or more VCF display-filter parameters
      via `SetParameter` in one call.
29. `list_agent_systems([catalog_path])`
    - Lists configured agent systems from catalog JSON.
    - `catalog_path` is optional (defaults to `assets/agent_systems.json`).
30. `ask_agent_system(project, system_id, prompt, [catalog_path], [allow_auto_exec], [execute_all], [execute_indices], [include_state_summary])`
    - Invokes one configured agent system through shared shell execution.
    - Returns table with `state`, `state_changed`, and `output`.
    - `project` may be `nil` to use an empty/default project state.
31. `plan_agent_system(project, system_id, prompt, [catalog_path], [include_state_summary], [max_candidates], [allow_mutating_candidates])`
    - Compiles prose into a typed `gentle.agent_plan_result.v1` plan through
      the shared `agents plan` route.
32. `execute_agent_plan(project, plan_json_or_table, candidate_id, [confirm])`
    - Executes one stored plan candidate through shared `agents execute-plan`
      and returns table with `state`, `state_changed`, and `output`.
33. `render_dotplot_svg(project, seq_id, dotplot_id, output_svg, [flex_track_id], [display_density_threshold], [display_intensity_gain])`
    - Convenience wrapper around engine `RenderDotplotSvg`.

### Lua example

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

[CLI Manual](../cli.md) | [Previous: JavaScript Shell](javascript.md) | [Next: Python Wrapper](python.md)
<!-- manual-page-nav:end -->
