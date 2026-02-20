use crate::app::GENtleApp;
use crate::dna_sequence::DNAsequence;
use crate::engine::{Engine, GentleEngine, Operation, ProjectState, Workflow};
use crate::enzymes::active_restriction_enzymes;
use crate::methylation_sites::MethylationMode;
use crate::{resource_sync, tf_motifs};
use mlua::prelude::*;
use mlua::{Error, Value};
use mlua::{IntoLuaMulti, Lua, MultiValue, Result as LuaResult};
use serde::Serialize;
use serde_json::json;

#[derive(Clone, Debug, Default)]
pub struct LuaInterface {
    lua: Lua,
}

impl LuaInterface {
    pub fn new() -> Self {
        Self { lua: Lua::new() }
    }

    pub fn lua(&self) -> &Lua {
        &self.lua
    }

    pub fn load_dna(path: &str) -> LuaResult<DNAsequence> {
        let mut dna =
            GENtleApp::load_from_file(path).map_err(|e| LuaError::RuntimeError(e.to_string()))?;

        // Add default enzymes and stuff
        *dna.restriction_enzymes_mut() = active_restriction_enzymes();
        dna.set_max_restriction_enzyme_sites(Some(2));
        dna.set_methylation_mode(MethylationMode::both());
        dna.update_computed_features();
        Ok(dna)
    }

    fn err(s: &str) -> Error {
        Error::RuntimeError(s.to_string())
    }

    pub fn help_main() {
        println!("Interactive Lua Shell (type 'exit' to quit)");
        println!("Available Rust functions:");
        println!("  - load_dna(filename): Loads a DNA sequence from a file");
        println!("  - write_gb(filename,seq): Writes a DNA sequence to a GenBank file");
        println!("  - load_project(filename): Loads a GENtle project JSON");
        println!("  - save_project(filename,project): Saves a GENtle project JSON");
        println!("  - capabilities(): Returns engine capabilities");
        println!("  - state_summary(project): Returns sequence/container summary");
        println!("  - inspect_dna_ladders([name_filter]): Returns built-in DNA ladder catalog");
        println!("  - export_dna_ladders(output_json, [name_filter]): Writes ladder catalog JSON");
        println!("  - inspect_rna_ladders([name_filter]): Returns built-in RNA ladder catalog");
        println!(
            "  - export_rna_ladders(output_json, [name_filter]): Writes RNA ladder catalog JSON"
        );
        println!(
            "  - apply_operation(project, op): Applies an operation (Lua table or JSON string)"
        );
        println!("  - apply_workflow(project, wf): Applies a workflow (Lua table or JSON string)");
        println!("  - list_reference_genomes([catalog_path]): Lists catalog genome names");
        println!(
            "  - is_reference_genome_prepared(genome_id, [catalog_path], [cache_dir]): Prepared check"
        );
        println!(
            "  - list_reference_genome_genes(genome_id, [catalog_path], [cache_dir]): Lists indexed genes"
        );
        println!(
            "  - blast_reference_genome(genome_id, query_sequence, [max_hits], [task], [catalog_path], [cache_dir]): BLAST query against prepared genome"
        );
        println!(
            "  - blast_helper_genome(helper_id, query_sequence, [max_hits], [task], [catalog_path], [cache_dir]): BLAST query against helper catalog genome"
        );
        println!(
            "  - prepare_genome(project, genome_id, [catalog_path], [cache_dir]): Engine op helper"
        );
        println!("  - extract_genome_region(project, genome_id, chr, start, end, [output_id], [catalog_path], [cache_dir]): Engine op helper");
        println!("  - extract_genome_gene(project, genome_id, gene_query, [occurrence], [output_id], [catalog_path], [cache_dir]): Engine op helper");
        println!("  - import_genome_bed_track(project, seq_id, path, [track_name], [min_score], [max_score], [clear_existing]): BED/BED.GZ overlay helper");
        println!("  - import_genome_bigwig_track(project, seq_id, path, [track_name], [min_score], [max_score], [clear_existing]): BigWig overlay helper");
        println!("  - import_genome_vcf_track(project, seq_id, path, [track_name], [min_score], [max_score], [clear_existing]): VCF/VCF.GZ variant overlay helper");
        println!("  - render_pool_gel_svg(project, ids_csv, output_svg, [ladders_csv]): Engine op helper");
        println!("  - sync_rebase(input, [output], [commercial_only]): Sync REBASE data");
        println!("  - sync_jaspar(input, [output]): Sync JASPAR motif data");
        println!("A sequence has the following properties:\n- seq.restriction_enzymes\n- seq.restriction_enzyme_sites\n- seq.open_reading_frames\n- seq.methylation_sites");
    }

    fn write_gb(seq: DNAsequence, filename: String) -> LuaResult<bool> {
        // lua.to_value(seq.restriction_enzyme_sites())
        seq.write_genbank_file(&filename)
            .map_err(|e| Self::err(&format!("{}", e)))?;
        Ok(true)
    }

    fn sync_rebase(
        input: String,
        output: Option<String>,
        commercial_only: Option<bool>,
    ) -> LuaResult<resource_sync::SyncReport> {
        let output = output.as_deref().map(str::trim).filter(|s| !s.is_empty());
        resource_sync::sync_rebase(&input, output, commercial_only.unwrap_or(false))
            .map_err(|e| Self::err(&e))
    }

    fn sync_jaspar(input: String, output: Option<String>) -> LuaResult<resource_sync::SyncReport> {
        let output = output.as_deref().map(str::trim).filter(|s| !s.is_empty());
        let report = resource_sync::sync_jaspar(&input, output).map_err(|e| Self::err(&e))?;
        tf_motifs::reload();
        Ok(report)
    }

    fn list_reference_genomes(catalog_path: Option<String>) -> LuaResult<Vec<String>> {
        let catalog_path = catalog_path
            .as_deref()
            .map(str::trim)
            .filter(|v| !v.is_empty());
        GentleEngine::list_reference_genomes(catalog_path).map_err(|e| Self::err(&e.to_string()))
    }

    fn inspect_dna_ladders(
        name_filter: Option<String>,
    ) -> LuaResult<crate::engine::DnaLadderCatalog> {
        let name_filter = name_filter
            .as_deref()
            .map(str::trim)
            .filter(|v| !v.is_empty());
        Ok(GentleEngine::inspect_dna_ladders(name_filter))
    }

    fn export_dna_ladders(
        output_json: String,
        name_filter: Option<String>,
    ) -> LuaResult<crate::engine::DnaLadderExportReport> {
        let name_filter = name_filter
            .as_deref()
            .map(str::trim)
            .filter(|v| !v.is_empty());
        GentleEngine::export_dna_ladders(&output_json, name_filter)
            .map_err(|e| Self::err(&e.to_string()))
    }

    fn inspect_rna_ladders(
        name_filter: Option<String>,
    ) -> LuaResult<crate::engine::RnaLadderCatalog> {
        let name_filter = name_filter
            .as_deref()
            .map(str::trim)
            .filter(|v| !v.is_empty());
        Ok(GentleEngine::inspect_rna_ladders(name_filter))
    }

    fn export_rna_ladders(
        output_json: String,
        name_filter: Option<String>,
    ) -> LuaResult<crate::engine::RnaLadderExportReport> {
        let name_filter = name_filter
            .as_deref()
            .map(str::trim)
            .filter(|v| !v.is_empty());
        GentleEngine::export_rna_ladders(&output_json, name_filter)
            .map_err(|e| Self::err(&e.to_string()))
    }

    fn is_reference_genome_prepared(
        genome_id: String,
        catalog_path: Option<String>,
        cache_dir: Option<String>,
    ) -> LuaResult<bool> {
        GentleEngine::is_reference_genome_prepared(
            catalog_path
                .as_deref()
                .map(str::trim)
                .filter(|v| !v.is_empty()),
            &genome_id,
            cache_dir
                .as_deref()
                .map(str::trim)
                .filter(|v| !v.is_empty()),
        )
        .map_err(|e| Self::err(&e.to_string()))
    }

    fn list_reference_genome_genes(
        genome_id: String,
        catalog_path: Option<String>,
        cache_dir: Option<String>,
    ) -> LuaResult<Vec<crate::genomes::GenomeGeneRecord>> {
        GentleEngine::list_reference_genome_genes(
            catalog_path
                .as_deref()
                .map(str::trim)
                .filter(|v| !v.is_empty()),
            &genome_id,
            cache_dir
                .as_deref()
                .map(str::trim)
                .filter(|v| !v.is_empty()),
        )
        .map_err(|e| Self::err(&e.to_string()))
    }

    fn blast_reference_genome(
        genome_id: String,
        query_sequence: String,
        max_hits: Option<usize>,
        task: Option<String>,
        catalog_path: Option<String>,
        cache_dir: Option<String>,
    ) -> LuaResult<crate::genomes::GenomeBlastReport> {
        GentleEngine::blast_reference_genome(
            catalog_path
                .as_deref()
                .map(str::trim)
                .filter(|v| !v.is_empty()),
            &genome_id,
            &query_sequence,
            max_hits.unwrap_or(25).max(1),
            task.as_deref().map(str::trim).filter(|v| !v.is_empty()),
            cache_dir
                .as_deref()
                .map(str::trim)
                .filter(|v| !v.is_empty()),
        )
        .map_err(|e| Self::err(&e.to_string()))
    }

    fn blast_helper_genome(
        genome_id: String,
        query_sequence: String,
        max_hits: Option<usize>,
        task: Option<String>,
        catalog_path: Option<String>,
        cache_dir: Option<String>,
    ) -> LuaResult<crate::genomes::GenomeBlastReport> {
        GentleEngine::blast_helper_genome(
            &genome_id,
            &query_sequence,
            max_hits.unwrap_or(25).max(1),
            task.as_deref().map(str::trim).filter(|v| !v.is_empty()),
            catalog_path
                .as_deref()
                .map(str::trim)
                .filter(|v| !v.is_empty()),
            cache_dir
                .as_deref()
                .map(str::trim)
                .filter(|v| !v.is_empty()),
        )
        .map_err(|e| Self::err(&e.to_string()))
    }

    // fn restriction_enzyme_digest(seq: DNAsequence, enzymes: String) -> LuaResult<Vec<DNAsequence>> {
    //     let enzymes = enzymes.split(',').map(|s| s.trim()).collect::<Vec<_>>();
    //     Ok(vec![])
    // }

    fn parse_or_decode<T>(lua: &Lua, input: Value) -> LuaResult<T>
    where
        T: serde::de::DeserializeOwned,
    {
        match input {
            Value::String(s) => {
                let s = s.to_str().map_err(|e| Self::err(&e.to_string()))?;
                serde_json::from_str(s.as_ref()).map_err(|e| Self::err(&format!("{e}")))
            }
            other => {
                let json_value = lua
                    .from_value::<serde_json::Value>(other)
                    .map_err(|e| Self::err(&format!("{e}")))?;
                serde_json::from_value(json_value).map_err(|e| Self::err(&format!("{e}")))
            }
        }
    }

    pub fn register_rust_functions(&self) -> LuaResult<()> {
        self.lua.globals().set(
            "load_dna",
            self.lua
                .create_function(|_lua, filename: String| Self::load_dna(&filename))?,
        )?;

        self.lua.globals().set(
            "write_gb",
            self.lua
                .create_function(|_lua, (filename, seq): (String, DNAsequence)| {
                    Self::write_gb(seq, filename)
                })?,
        )?;

        self.lua.globals().set(
            "load_project",
            self.lua.create_function(|lua, filename: String| {
                let state = ProjectState::load_from_path(&filename)
                    .map_err(|e| Self::err(&e.to_string()))?;
                lua.to_value(&state)
            })?,
        )?;

        self.lua.globals().set(
            "save_project",
            self.lua
                .create_function(|lua, (filename, state): (String, Value)| {
                    let state: ProjectState = lua
                        .from_value(state)
                        .map_err(|e| Self::err(&format!("Invalid project value: {e}")))?;
                    state
                        .save_to_path(&filename)
                        .map_err(|e| Self::err(&e.to_string()))?;
                    Ok(true)
                })?,
        )?;

        self.lua.globals().set(
            "capabilities",
            self.lua
                .create_function(|lua, _: ()| lua.to_value(&GentleEngine::capabilities()))?,
        )?;

        self.lua.globals().set(
            "state_summary",
            self.lua.create_function(|lua, state: Value| {
                let state: ProjectState = lua
                    .from_value(state)
                    .map_err(|e| Self::err(&format!("Invalid project value: {e}")))?;
                let engine = GentleEngine::from_state(state);
                lua.to_value(&engine.summarize_state())
            })?,
        )?;

        self.lua.globals().set(
            "inspect_dna_ladders",
            self.lua
                .create_function(|lua, name_filter: Option<String>| {
                    let ladders = Self::inspect_dna_ladders(name_filter)?;
                    lua.to_value(&ladders)
                })?,
        )?;

        self.lua.globals().set(
            "list_dna_ladders",
            self.lua
                .create_function(|lua, name_filter: Option<String>| {
                    let ladders = Self::inspect_dna_ladders(name_filter)?;
                    lua.to_value(&ladders)
                })?,
        )?;

        self.lua.globals().set(
            "export_dna_ladders",
            self.lua.create_function(
                |lua, (output_json, name_filter): (String, Option<String>)| {
                    let report = Self::export_dna_ladders(output_json, name_filter)?;
                    lua.to_value(&report)
                },
            )?,
        )?;

        self.lua.globals().set(
            "inspect_rna_ladders",
            self.lua
                .create_function(|lua, name_filter: Option<String>| {
                    let ladders = Self::inspect_rna_ladders(name_filter)?;
                    lua.to_value(&ladders)
                })?,
        )?;

        self.lua.globals().set(
            "list_rna_ladders",
            self.lua
                .create_function(|lua, name_filter: Option<String>| {
                    let ladders = Self::inspect_rna_ladders(name_filter)?;
                    lua.to_value(&ladders)
                })?,
        )?;

        self.lua.globals().set(
            "export_rna_ladders",
            self.lua.create_function(
                |lua, (output_json, name_filter): (String, Option<String>)| {
                    let report = Self::export_rna_ladders(output_json, name_filter)?;
                    lua.to_value(&report)
                },
            )?,
        )?;

        self.lua.globals().set(
            "list_reference_genomes",
            self.lua
                .create_function(|lua, catalog_path: Option<String>| {
                    let genomes = Self::list_reference_genomes(catalog_path)?;
                    lua.to_value(&genomes)
                })?,
        )?;

        self.lua.globals().set(
            "is_reference_genome_prepared",
            self.lua.create_function(
                |lua, (genome_id, catalog_path, cache_dir): (String, Option<String>, Option<String>)| {
                    let prepared =
                        Self::is_reference_genome_prepared(genome_id, catalog_path, cache_dir)?;
                    lua.to_value(&prepared)
                },
            )?,
        )?;

        self.lua.globals().set(
            "list_reference_genome_genes",
            self.lua.create_function(
                |lua, (genome_id, catalog_path, cache_dir): (String, Option<String>, Option<String>)| {
                    let genes = Self::list_reference_genome_genes(genome_id, catalog_path, cache_dir)?;
                    lua.to_value(&genes)
                },
            )?,
        )?;

        self.lua.globals().set(
            "blast_reference_genome",
            self.lua.create_function(
                |lua,
                 (genome_id, query_sequence, max_hits, task, catalog_path, cache_dir): (
                    String,
                    String,
                    Option<usize>,
                    Option<String>,
                    Option<String>,
                    Option<String>,
                )| {
                    let report = Self::blast_reference_genome(
                        genome_id,
                        query_sequence,
                        max_hits,
                        task,
                        catalog_path,
                        cache_dir,
                    )?;
                    lua.to_value(&report)
                },
            )?,
        )?;

        self.lua.globals().set(
            "blast_helper_genome",
            self.lua.create_function(
                |lua,
                 (genome_id, query_sequence, max_hits, task, catalog_path, cache_dir): (
                    String,
                    String,
                    Option<usize>,
                    Option<String>,
                    Option<String>,
                    Option<String>,
                )| {
                    let report = Self::blast_helper_genome(
                        genome_id,
                        query_sequence,
                        max_hits,
                        task,
                        catalog_path,
                        cache_dir,
                    )?;
                    lua.to_value(&report)
                },
            )?,
        )?;

        self.lua.globals().set(
            "apply_operation",
            self.lua
                .create_function(|lua, (state, op): (Value, Value)| {
                    let state: ProjectState = lua
                        .from_value(state)
                        .map_err(|e| Self::err(&format!("Invalid project value: {e}")))?;
                    let op: Operation = Self::parse_or_decode(lua, op)?;
                    let mut engine = GentleEngine::from_state(state);
                    let result = engine.apply(op).map_err(|e| Self::err(&e.to_string()))?;
                    #[derive(Serialize)]
                    struct Response {
                        state: ProjectState,
                        result: crate::engine::OpResult,
                    }
                    let response = Response {
                        state: engine.state().clone(),
                        result,
                    };
                    lua.to_value(&response)
                })?,
        )?;

        self.lua.globals().set(
            "prepare_genome",
            self.lua.create_function(
                |lua,
                 (state, genome_id, catalog_path, cache_dir): (
                    Value,
                    String,
                    Option<String>,
                    Option<String>,
                )| {
                    let state: ProjectState = lua
                        .from_value(state)
                        .map_err(|e| Self::err(&format!("Invalid project value: {e}")))?;
                    let mut engine = GentleEngine::from_state(state);
                    let result = engine
                        .apply(Operation::PrepareGenome {
                            genome_id,
                            catalog_path,
                            cache_dir,
                        })
                        .map_err(|e| Self::err(&e.to_string()))?;
                    #[derive(Serialize)]
                    struct Response {
                        state: ProjectState,
                        result: crate::engine::OpResult,
                    }
                    lua.to_value(&Response {
                        state: engine.state().clone(),
                        result,
                    })
                },
            )?,
        )?;

        self.lua.globals().set(
            "extract_genome_region",
            self.lua.create_function(
                |lua,
                 (
                    state,
                    genome_id,
                    chromosome,
                    start_1based,
                    end_1based,
                    output_id,
                    catalog_path,
                    cache_dir,
                ): (
                    Value,
                    String,
                    String,
                    usize,
                    usize,
                    Option<String>,
                    Option<String>,
                    Option<String>,
                )| {
                    let state: ProjectState = lua
                        .from_value(state)
                        .map_err(|e| Self::err(&format!("Invalid project value: {e}")))?;
                    let mut engine = GentleEngine::from_state(state);
                    let result = engine
                        .apply(Operation::ExtractGenomeRegion {
                            genome_id,
                            chromosome,
                            start_1based,
                            end_1based,
                            output_id,
                            catalog_path,
                            cache_dir,
                        })
                        .map_err(|e| Self::err(&e.to_string()))?;
                    #[derive(Serialize)]
                    struct Response {
                        state: ProjectState,
                        result: crate::engine::OpResult,
                    }
                    lua.to_value(&Response {
                        state: engine.state().clone(),
                        result,
                    })
                },
            )?,
        )?;

        self.lua.globals().set(
            "extract_genome_gene",
            self.lua.create_function(
                |lua, (state, genome_id, gene_query, occurrence, output_id, catalog_path, cache_dir): (Value, String, String, Option<usize>, Option<String>, Option<String>, Option<String>)| {
                    let state: ProjectState = lua
                        .from_value(state)
                        .map_err(|e| Self::err(&format!("Invalid project value: {e}")))?;
                    let mut engine = GentleEngine::from_state(state);
                    let result = engine
                        .apply(Operation::ExtractGenomeGene {
                            genome_id,
                            gene_query,
                            occurrence,
                            output_id,
                            catalog_path,
                            cache_dir,
                        })
                        .map_err(|e| Self::err(&e.to_string()))?;
                    #[derive(Serialize)]
                    struct Response {
                        state: ProjectState,
                        result: crate::engine::OpResult,
                    }
                    lua.to_value(&Response {
                        state: engine.state().clone(),
                        result,
                    })
                },
            )?,
        )?;

        self.lua.globals().set(
            "import_genome_bed_track",
            self.lua.create_function(
                |lua,
                 (state, seq_id, path, track_name, min_score, max_score, clear_existing): (
                    Value,
                    String,
                    String,
                    Option<String>,
                    Option<f64>,
                    Option<f64>,
                    Option<bool>,
                )| {
                    let state: ProjectState = lua
                        .from_value(state)
                        .map_err(|e| Self::err(&format!("Invalid project value: {e}")))?;
                    let mut engine = GentleEngine::from_state(state);
                    let result = engine
                        .apply(Operation::ImportGenomeBedTrack {
                            seq_id,
                            path,
                            track_name,
                            min_score,
                            max_score,
                            clear_existing,
                        })
                        .map_err(|e| Self::err(&e.to_string()))?;
                    #[derive(Serialize)]
                    struct Response {
                        state: ProjectState,
                        result: crate::engine::OpResult,
                    }
                    lua.to_value(&Response {
                        state: engine.state().clone(),
                        result,
                    })
                },
            )?,
        )?;

        self.lua.globals().set(
            "import_genome_bigwig_track",
            self.lua.create_function(
                |lua,
                 (state, seq_id, path, track_name, min_score, max_score, clear_existing): (
                    Value,
                    String,
                    String,
                    Option<String>,
                    Option<f64>,
                    Option<f64>,
                    Option<bool>,
                )| {
                    let state: ProjectState = lua
                        .from_value(state)
                        .map_err(|e| Self::err(&format!("Invalid project value: {e}")))?;
                    let mut engine = GentleEngine::from_state(state);
                    let result = engine
                        .apply(Operation::ImportGenomeBigWigTrack {
                            seq_id,
                            path,
                            track_name,
                            min_score,
                            max_score,
                            clear_existing,
                        })
                        .map_err(|e| Self::err(&e.to_string()))?;
                    #[derive(Serialize)]
                    struct Response {
                        state: ProjectState,
                        result: crate::engine::OpResult,
                    }
                    lua.to_value(&Response {
                        state: engine.state().clone(),
                        result,
                    })
                },
            )?,
        )?;

        self.lua.globals().set(
            "import_genome_vcf_track",
            self.lua.create_function(
                |lua,
                 (state, seq_id, path, track_name, min_score, max_score, clear_existing): (
                    Value,
                    String,
                    String,
                    Option<String>,
                    Option<f64>,
                    Option<f64>,
                    Option<bool>,
                )| {
                    let state: ProjectState = lua
                        .from_value(state)
                        .map_err(|e| Self::err(&format!("Invalid project value: {e}")))?;
                    let mut engine = GentleEngine::from_state(state);
                    let result = engine
                        .apply(Operation::ImportGenomeVcfTrack {
                            seq_id,
                            path,
                            track_name,
                            min_score,
                            max_score,
                            clear_existing,
                        })
                        .map_err(|e| Self::err(&e.to_string()))?;
                    #[derive(Serialize)]
                    struct Response {
                        state: ProjectState,
                        result: crate::engine::OpResult,
                    }
                    lua.to_value(&Response {
                        state: engine.state().clone(),
                        result,
                    })
                },
            )?,
        )?;

        self.lua.globals().set(
            "render_pool_gel_svg",
            self.lua.create_function(
                |lua,
                 (state, ids_csv, output_svg, ladders_csv): (
                    Value,
                    String,
                    String,
                    Option<String>,
                )| {
                    let state: ProjectState = lua
                        .from_value(state)
                        .map_err(|e| Self::err(&format!("Invalid project value: {e}")))?;
                    let inputs = ids_csv
                        .split(',')
                        .map(|s| s.trim())
                        .filter(|s| !s.is_empty())
                        .map(|s| s.to_string())
                        .collect::<Vec<_>>();
                    if inputs.is_empty() {
                        return Err(Self::err(
                            "render_pool_gel_svg requires at least one sequence id",
                        ));
                    }
                    let ladders = ladders_csv
                        .as_deref()
                        .map(|s| {
                            s.split(',')
                                .map(|v| v.trim())
                                .filter(|v| !v.is_empty())
                                .map(|v| v.to_string())
                                .collect::<Vec<_>>()
                        })
                        .filter(|v| !v.is_empty());
                    let mut engine = GentleEngine::from_state(state);
                    let result = engine
                        .apply(Operation::RenderPoolGelSvg {
                            inputs,
                            path: output_svg,
                            ladders,
                        })
                        .map_err(|e| Self::err(&e.to_string()))?;
                    #[derive(Serialize)]
                    struct Response {
                        state: ProjectState,
                        result: crate::engine::OpResult,
                    }
                    lua.to_value(&Response {
                        state: engine.state().clone(),
                        result,
                    })
                },
            )?,
        )?;

        self.lua.globals().set(
            "apply_workflow",
            self.lua
                .create_function(|lua, (state, workflow): (Value, Value)| {
                    let state: ProjectState = lua
                        .from_value(state)
                        .map_err(|e| Self::err(&format!("Invalid project value: {e}")))?;
                    let workflow: Workflow = Self::parse_or_decode(lua, workflow)?;
                    let mut engine = GentleEngine::from_state(state);
                    let results = engine
                        .apply_workflow(workflow)
                        .map_err(|e| Self::err(&e.to_string()))?;
                    #[derive(Serialize)]
                    struct Response {
                        state: ProjectState,
                        results: Vec<crate::engine::OpResult>,
                    }
                    let response = Response {
                        state: engine.state().clone(),
                        results,
                    };
                    lua.to_value(&response)
                })?,
        )?;

        self.lua.globals().set(
            "sync_rebase",
            self.lua.create_function(
                |lua, (input, output, commercial_only): (String, Option<String>, Option<bool>)| {
                    let report = Self::sync_rebase(input, output, commercial_only)?;
                    lua.to_value(&report)
                },
            )?,
        )?;

        self.lua.globals().set(
            "sync_jaspar",
            self.lua
                .create_function(|lua, (input, output): (String, Option<String>)| {
                    let report = Self::sync_jaspar(input, output)?;
                    lua.to_value(&report)
                })?,
        )?;

        Ok(())
    }

    pub fn format_lua_value(value: &Value) -> String {
        format!("{}", json!(value))
    }
}

impl IntoLuaMulti for DNAsequence {
    fn into_lua_multi(self, lua: &Lua) -> LuaResult<MultiValue> {
        // Convert the struct to a Lua table using serde
        let table = lua.to_value(&self)?;
        let mut ret = vec![table];
        json!(self).as_object().unwrap().iter().for_each(|(_k, v)| {
            let value = lua.to_value(v).unwrap();
            ret.push(value);
        });

        // Return both the table and individual values
        Ok(MultiValue::from_vec(ret))
    }
}

impl FromLuaMulti for DNAsequence {
    fn from_lua_multi(values: MultiValue, _lua: &Lua) -> LuaResult<Self> {
        let table = values.front().unwrap();
        let table = json!(table);
        let ret: DNAsequence = serde_json::from_value(table).unwrap();
        Ok(ret)
    }
}
