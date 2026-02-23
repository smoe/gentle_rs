use crate::app::GENtleApp;
use crate::dna_sequence::DNAsequence;
use crate::engine::{
    Engine, FeatureExpertTarget, GenomeAnchorSide, GentleEngine, Operation, ProjectState, Workflow,
};
use crate::engine_shell::{ShellCommand, execute_shell_command};
use crate::enzymes::active_restriction_enzymes;
use crate::methylation_sites::MethylationMode;
use crate::resource_sync;
use mlua::prelude::*;
use mlua::{Error, Value};
use mlua::{IntoLuaMulti, Lua, MultiValue, Result as LuaResult};
use serde::Serialize;
use serde_json::json;

#[derive(Clone, Debug, Default)]
pub struct LuaInterface {
    lua: Lua,
}

#[derive(Serialize)]
struct ShellUtilityResponse {
    state: ProjectState,
    state_changed: bool,
    output: serde_json::Value,
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

    fn sync_report_from_shell_output(
        output: serde_json::Value,
        context: &str,
    ) -> LuaResult<resource_sync::SyncReport> {
        let Some(report_value) = output.get("report").cloned() else {
            return Err(Self::err(&format!(
                "{context}: missing 'report' field in shell output"
            )));
        };
        serde_json::from_value(report_value)
            .map_err(|e| Self::err(&format!("{context}: invalid sync report output: {e}")))
    }

    pub fn help_main() {
        println!("Interactive Lua Shell (type 'exit' to quit)");
        println!("Available Rust functions:");
        println!("  - load_dna(filename): Loads a DNA sequence from a file");
        println!("  - write_gb(seq, filename): Writes a DNA sequence to a GenBank file");
        println!("  - load_project(filename): Loads a GENtle project JSON");
        println!("  - save_project(state, filename): Saves a GENtle project JSON");
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
        println!("  - set_parameter(project, name, value): Helper for Operation::SetParameter");
        println!(
            "  - set_vcf_display_filter(project, opts): Convenience helper for VCF display criteria updates"
        );
        println!("  - apply_workflow(project, wf): Applies a workflow (Lua table or JSON string)");
        println!(
            "  - inspect_feature_expert(project, seq_id, target): Returns TFBS/restriction expert JSON (target table or JSON string)"
        );
        println!(
            "  - render_feature_expert_svg(project, seq_id, target, output_svg): Writes expert-view SVG via shared engine operation"
        );
        println!("  - list_reference_genomes([catalog_path]): Lists catalog genome names");
        println!(
            "  - list_agent_systems([catalog_path]): Lists external/internal AI systems from agent catalog"
        );
        println!(
            "  - ask_agent_system(project_or_nil, system_id, prompt, [catalog_path], [allow_auto_exec], [execute_all], [execute_indices], [include_state_summary]): Asks one configured AI system"
        );
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
        println!(
            "  - extract_genome_region(project, genome_id, chr, start, end, [output_id], [catalog_path], [cache_dir]): Engine op helper"
        );
        println!(
            "  - extract_genome_gene(project, genome_id, gene_query, [occurrence], [output_id], [catalog_path], [cache_dir]): Engine op helper"
        );
        println!(
            "  - extend_genome_anchor(project, seq_id, side_5p_or_3p, length_bp, [output_id], [catalog_path], [cache_dir]): Engine op helper"
        );
        println!(
            "  - import_genome_bed_track(project, seq_id, path, [track_name], [min_score], [max_score], [clear_existing]): BED/BED.GZ overlay helper"
        );
        println!(
            "  - import_genome_bigwig_track(project, seq_id, path, [track_name], [min_score], [max_score], [clear_existing]): BigWig overlay helper"
        );
        println!(
            "  - import_genome_vcf_track(project, seq_id, path, [track_name], [min_score], [max_score], [clear_existing]): VCF/VCF.GZ variant overlay helper"
        );
        println!(
            "  - render_pool_gel_svg(project, ids_csv, output_svg, [ladders_csv]): Engine op helper"
        );
        println!(
            "  - import_pool(project, input_pool_json, [prefix]): Shared adapter helper (returns updated state + report)"
        );
        println!("  - sync_rebase(input, [output], [commercial_only]): Sync REBASE data");
        println!("  - sync_jaspar(input, [output]): Sync JASPAR motif data");
        println!(
            "A sequence has the following properties:\n- seq.restriction_enzymes\n- seq.restriction_enzyme_sites\n- seq.open_reading_frames\n- seq.methylation_sites"
        );
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
        let mut engine = GentleEngine::from_state(ProjectState::default());
        let output = output
            .as_deref()
            .map(str::trim)
            .filter(|s| !s.is_empty())
            .map(str::to_string);
        let command = ShellCommand::ResourcesSyncRebase {
            input,
            output,
            commercial_only: commercial_only.unwrap_or(false),
        };
        let run = execute_shell_command(&mut engine, &command)
            .map_err(|e| Self::err(&format!("resources sync-rebase failed: {e}")))?;
        Self::sync_report_from_shell_output(run.output, "resources sync-rebase")
    }

    fn sync_jaspar(input: String, output: Option<String>) -> LuaResult<resource_sync::SyncReport> {
        let mut engine = GentleEngine::from_state(ProjectState::default());
        let output = output
            .as_deref()
            .map(str::trim)
            .filter(|s| !s.is_empty())
            .map(str::to_string);
        let command = ShellCommand::ResourcesSyncJaspar { input, output };
        let run = execute_shell_command(&mut engine, &command)
            .map_err(|e| Self::err(&format!("resources sync-jaspar failed: {e}")))?;
        Self::sync_report_from_shell_output(run.output, "resources sync-jaspar")
    }

    fn import_pool(
        state: ProjectState,
        input: String,
        prefix: Option<String>,
    ) -> LuaResult<ShellUtilityResponse> {
        let mut engine = GentleEngine::from_state(state);
        let command = ShellCommand::ImportPool {
            input,
            prefix: prefix
                .as_deref()
                .map(str::trim)
                .filter(|v| !v.is_empty())
                .unwrap_or("pool")
                .to_string(),
        };
        let run = execute_shell_command(&mut engine, &command)
            .map_err(|e| Self::err(&format!("import-pool failed: {e}")))?;
        Ok(ShellUtilityResponse {
            state: engine.state().clone(),
            state_changed: run.state_changed,
            output: run.output,
        })
    }

    fn list_reference_genomes(catalog_path: Option<String>) -> LuaResult<Vec<String>> {
        let catalog_path = catalog_path
            .as_deref()
            .map(str::trim)
            .filter(|v| !v.is_empty());
        GentleEngine::list_reference_genomes(catalog_path).map_err(|e| Self::err(&e.to_string()))
    }

    fn list_agent_systems(catalog_path: Option<String>) -> LuaResult<serde_json::Value> {
        let mut engine = GentleEngine::from_state(ProjectState::default());
        let command = ShellCommand::AgentsList {
            catalog_path: catalog_path
                .as_deref()
                .map(str::trim)
                .filter(|v| !v.is_empty())
                .map(str::to_string),
        };
        let run = execute_shell_command(&mut engine, &command)
            .map_err(|e| Self::err(&format!("agents list failed: {e}")))?;
        Ok(run.output)
    }

    #[allow(clippy::too_many_arguments)]
    fn ask_agent_system(
        state: Option<ProjectState>,
        system_id: String,
        prompt: String,
        catalog_path: Option<String>,
        allow_auto_exec: Option<bool>,
        execute_all: Option<bool>,
        execute_indices: Option<Vec<usize>>,
        include_state_summary: Option<bool>,
    ) -> LuaResult<ShellUtilityResponse> {
        let mut engine = GentleEngine::from_state(state.unwrap_or_default());
        let command = ShellCommand::AgentsAsk {
            system_id,
            prompt,
            catalog_path: catalog_path
                .as_deref()
                .map(str::trim)
                .filter(|v| !v.is_empty())
                .map(str::to_string),
            base_url_override: None,
            model_override: None,
            timeout_seconds: None,
            connect_timeout_seconds: None,
            read_timeout_seconds: None,
            max_retries: None,
            max_response_bytes: None,
            include_state_summary: include_state_summary.unwrap_or(true),
            allow_auto_exec: allow_auto_exec.unwrap_or(false),
            execute_all: execute_all.unwrap_or(false),
            execute_indices: execute_indices.unwrap_or_default(),
        };
        let run = execute_shell_command(&mut engine, &command)
            .map_err(|e| Self::err(&format!("agents ask failed: {e}")))?;
        Ok(ShellUtilityResponse {
            state: engine.state().clone(),
            state_changed: run.state_changed,
            output: run.output,
        })
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
                .create_function(|lua, (seq_value, filename): (Value, String)| {
                    let seq: DNAsequence = lua
                        .from_value(seq_value)
                        .map_err(|e| Self::err(&format!("Invalid sequence value: {e}")))?;
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
                .create_function(|lua, (state, filename): (Value, String)| {
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
            "list_agent_systems",
            self.lua
                .create_function(|lua, catalog_path: Option<String>| {
                    let systems = Self::list_agent_systems(catalog_path)?;
                    lua.to_value(&systems)
                })?,
        )?;

        self.lua.globals().set(
            "ask_agent_system",
            self.lua.create_function(
                |lua,
                 (
                    state,
                    system_id,
                    prompt,
                    catalog_path,
                    allow_auto_exec,
                    execute_all,
                    execute_indices,
                    include_state_summary,
                ): (
                    Value,
                    String,
                    String,
                    Option<String>,
                    Option<bool>,
                    Option<bool>,
                    Option<Vec<usize>>,
                    Option<bool>,
                )| {
                    let state = if matches!(state, Value::Nil) {
                        None
                    } else {
                        Some(
                            lua.from_value(state)
                                .map_err(|e| Self::err(&format!("Invalid project value: {e}")))?,
                        )
                    };
                    let response = Self::ask_agent_system(
                        state,
                        system_id,
                        prompt,
                        catalog_path,
                        allow_auto_exec,
                        execute_all,
                        execute_indices,
                        include_state_summary,
                    )?;
                    lua.to_value(&response)
                },
            )?,
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
            "inspect_feature_expert",
            self.lua
                .create_function(|lua, (state, seq_id, target): (Value, String, Value)| {
                    let state: ProjectState = lua
                        .from_value(state)
                        .map_err(|e| Self::err(&format!("Invalid project value: {e}")))?;
                    let target: FeatureExpertTarget = Self::parse_or_decode(lua, target)?;
                    let engine = GentleEngine::from_state(state);
                    let view = engine
                        .inspect_feature_expert(&seq_id, &target)
                        .map_err(|e| Self::err(&e.to_string()))?;
                    lua.to_value(&view)
                })?,
        )?;

        self.lua.globals().set(
            "render_feature_expert_svg",
            self.lua.create_function(
                |lua, (state, seq_id, target, path): (Value, String, Value, String)| {
                    let state: ProjectState = lua
                        .from_value(state)
                        .map_err(|e| Self::err(&format!("Invalid project value: {e}")))?;
                    let target: FeatureExpertTarget = Self::parse_or_decode(lua, target)?;
                    let mut engine = GentleEngine::from_state(state);
                    let result = engine
                        .apply(Operation::RenderFeatureExpertSvg {
                            seq_id,
                            target,
                            path,
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
            "set_parameter",
            self.lua
                .create_function(|lua, (state, name, value): (Value, String, Value)| {
                    let state: ProjectState = lua
                        .from_value(state)
                        .map_err(|e| Self::err(&format!("Invalid project value: {e}")))?;
                    let value: serde_json::Value = lua
                        .from_value(value)
                        .map_err(|e| Self::err(&format!("Invalid parameter value: {e}")))?;
                    let mut engine = GentleEngine::from_state(state);
                    let result = engine
                        .apply(Operation::SetParameter { name, value })
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
                })?,
        )?;

        self.lua.globals().set(
            "set_vcf_display_filter",
            self.lua
                .create_function(|lua, (state, options): (Value, Value)| {
                    let state: ProjectState = lua
                        .from_value(state)
                        .map_err(|e| Self::err(&format!("Invalid project value: {e}")))?;
                    let options: serde_json::Value = lua
                        .from_value(options)
                        .map_err(|e| Self::err(&format!("Invalid options table: {e}")))?;
                    let obj = options.as_object().ok_or_else(|| {
                        Self::err(
                            "set_vcf_display_filter expects a Lua table/object with optional fields",
                        )
                    })?;
                    let mut engine = GentleEngine::from_state(state);
                    let mut results: Vec<crate::engine::OpResult> = vec![];
                    let mut apply_one = |name: &str, value: serde_json::Value| -> LuaResult<()> {
                        let result = engine
                            .apply(Operation::SetParameter {
                                name: name.to_string(),
                                value,
                            })
                            .map_err(|e| Self::err(&e.to_string()))?;
                        results.push(result);
                        Ok(())
                    };
                    if let Some(v) = obj.get("show_snp") {
                        apply_one("vcf_display_show_snp", v.clone())?;
                    }
                    if let Some(v) = obj.get("show_ins") {
                        apply_one("vcf_display_show_ins", v.clone())?;
                    }
                    if let Some(v) = obj.get("show_del") {
                        apply_one("vcf_display_show_del", v.clone())?;
                    }
                    if let Some(v) = obj.get("show_sv") {
                        apply_one("vcf_display_show_sv", v.clone())?;
                    }
                    if let Some(v) = obj.get("show_other") {
                        apply_one("vcf_display_show_other", v.clone())?;
                    }
                    if let Some(v) = obj.get("pass_only") {
                        apply_one("vcf_display_pass_only", v.clone())?;
                    }
                    if let Some(v) = obj.get("use_min_qual") {
                        apply_one("vcf_display_use_min_qual", v.clone())?;
                    }
                    if let Some(v) = obj.get("min_qual") {
                        apply_one("vcf_display_min_qual", v.clone())?;
                    }
                    if let Some(v) = obj.get("use_max_qual") {
                        apply_one("vcf_display_use_max_qual", v.clone())?;
                    }
                    if let Some(v) = obj.get("max_qual") {
                        apply_one("vcf_display_max_qual", v.clone())?;
                    }
                    if let Some(v) = obj.get("required_info_keys") {
                        apply_one("vcf_display_required_info_keys", v.clone())?;
                    }
                    #[derive(Serialize)]
                    struct Response {
                        state: ProjectState,
                        applied: usize,
                        results: Vec<crate::engine::OpResult>,
                    }
                    lua.to_value(&Response {
                        state: engine.state().clone(),
                        applied: results.len(),
                        results,
                    })
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
                            timeout_seconds: None,
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
            "extend_genome_anchor",
            self.lua.create_function(
                |lua,
                 (state, seq_id, side_raw, length_bp, output_id, catalog_path, cache_dir): (
                    Value,
                    String,
                    String,
                    usize,
                    Option<String>,
                    Option<String>,
                    Option<String>,
                )| {
                    let side = match side_raw.trim().to_ascii_lowercase().as_str() {
                        "5" | "5p" | "5prime" | "5'" | "five_prime" | "five-prime" => {
                            GenomeAnchorSide::FivePrime
                        }
                        "3" | "3p" | "3prime" | "3'" | "three_prime" | "three-prime" => {
                            GenomeAnchorSide::ThreePrime
                        }
                        _ => return Err(Self::err("extend_genome_anchor side must be 5p or 3p")),
                    };
                    let state: ProjectState = lua
                        .from_value(state)
                        .map_err(|e| Self::err(&format!("Invalid project value: {e}")))?;
                    let mut engine = GentleEngine::from_state(state);
                    let result = engine
                        .apply(Operation::ExtendGenomeAnchor {
                            seq_id,
                            side,
                            length_bp,
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
                            container_ids: None,
                            arrangement_id: None,
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

        self.lua.globals().set(
            "import_pool",
            self.lua.create_function(
                |lua, (state, input, prefix): (Value, String, Option<String>)| {
                    let state: ProjectState = lua
                        .from_value(state)
                        .map_err(|e| Self::err(&format!("Invalid project value: {e}")))?;
                    let response = Self::import_pool(state, input, prefix)?;
                    lua.to_value(&response)
                },
            )?,
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

#[cfg(test)]
mod tests {
    use super::*;
    use serde_json::json;
    use std::fs;
    use tempfile::tempdir;

    #[test]
    fn lua_sync_rebase_wrapper_writes_snapshot() {
        let td = tempdir().expect("tempdir");
        let input_path = td.path().join("rebase.withrefm");
        let output_path = td.path().join("rebase.json");
        fs::write(
            &input_path,
            "<1>EcoRI\n<2>EcoRI\n<3>GAATTC (1/5)\n<7>N\n//\n",
        )
        .expect("write rebase input");
        let report = LuaInterface::sync_rebase(
            input_path.to_string_lossy().to_string(),
            Some(output_path.to_string_lossy().to_string()),
            Some(true),
        )
        .expect("sync rebase");
        assert_eq!(report.resource, "rebase-commercial");
        assert_eq!(report.item_count, 1);
        assert!(output_path.exists());
    }

    #[test]
    fn lua_sync_jaspar_wrapper_writes_snapshot() {
        let td = tempdir().expect("tempdir");
        let input_path = td.path().join("motifs.pfm");
        let output_path = td.path().join("motifs.json");
        fs::write(
            &input_path,
            ">MA0001.1 TEST\nA [ 10 0 0 0 ]\nC [ 0 10 0 0 ]\nG [ 0 0 10 0 ]\nT [ 0 0 0 10 ]\n",
        )
        .expect("write jaspar input");
        let report = LuaInterface::sync_jaspar(
            input_path.to_string_lossy().to_string(),
            Some(output_path.to_string_lossy().to_string()),
        )
        .expect("sync jaspar");
        assert_eq!(report.resource, "jaspar");
        assert_eq!(report.item_count, 1);
        assert!(output_path.exists());
    }

    #[test]
    fn lua_import_pool_wrapper_loads_member_into_state() {
        let td = tempdir().expect("tempdir");
        let pool_path = td.path().join("demo.pool.gentle.json");
        let pool_json = json!({
            "schema": "gentle.pool.v1",
            "pool_id": "demo_pool",
            "human_id": "demo",
            "member_count": 1,
            "members": [
                {
                    "seq_id": "member_1",
                    "human_id": "member_1",
                    "name": "Member One",
                    "sequence": "ATGCATGC",
                    "length_bp": 8,
                    "topology": "linear",
                    "ends": {
                        "end_type": "blunt",
                        "forward_5": "",
                        "forward_3": "",
                        "reverse_5": "",
                        "reverse_3": ""
                    }
                }
            ]
        });
        fs::write(
            &pool_path,
            serde_json::to_string_pretty(&pool_json).expect("serialize pool json"),
        )
        .expect("write pool json");
        let out = LuaInterface::import_pool(
            ProjectState::default(),
            pool_path.to_string_lossy().to_string(),
            Some("lua_pool".to_string()),
        )
        .expect("import pool");
        assert!(out.state_changed);
        assert_eq!(out.output["pool_id"].as_str(), Some("demo_pool"));
        assert_eq!(out.output["member_count"].as_u64(), Some(1));
        let imported_ids = out
            .output
            .get("imported_ids")
            .and_then(|v| v.as_array())
            .expect("imported_ids array");
        assert_eq!(imported_ids.len(), 1);
        let imported_id = imported_ids[0].as_str().expect("imported id string");
        assert!(out.state.sequences.contains_key(imported_id));
    }

    #[test]
    fn lua_list_agent_systems_wrapper_reads_catalog() {
        let td = tempdir().expect("tempdir");
        let catalog_path = td.path().join("agents.json");
        fs::write(
            &catalog_path,
            r#"{
  "schema": "gentle.agent_systems.v1",
  "systems": [
    { "id": "builtin_echo", "label": "Built-in Echo", "transport": "builtin_echo" }
  ]
}"#,
        )
        .expect("write agent catalog");
        let out =
            LuaInterface::list_agent_systems(Some(catalog_path.to_string_lossy().to_string()))
                .expect("list agent systems");
        assert_eq!(out["schema"].as_str(), Some("gentle.agent_systems_list.v1"));
        assert_eq!(out["system_count"].as_u64(), Some(1));
        assert_eq!(out["systems"][0]["id"].as_str(), Some("builtin_echo"));
    }

    #[test]
    fn lua_ask_agent_system_wrapper_uses_builtin_echo() {
        let td = tempdir().expect("tempdir");
        let catalog_path = td.path().join("agents.json");
        fs::write(
            &catalog_path,
            r#"{
  "schema": "gentle.agent_systems.v1",
  "systems": [
    { "id": "builtin_echo", "label": "Built-in Echo", "transport": "builtin_echo" }
  ]
}"#,
        )
        .expect("write agent catalog");
        let out = LuaInterface::ask_agent_system(
            Some(ProjectState::default()),
            "builtin_echo".to_string(),
            "hello".to_string(),
            Some(catalog_path.to_string_lossy().to_string()),
            Some(false),
            Some(false),
            Some(vec![]),
            Some(true),
        )
        .expect("ask agent");
        assert!(!out.state_changed);
        assert_eq!(
            out.output["schema"].as_str(),
            Some("gentle.agent_ask_result.v1")
        );
        assert_eq!(
            out.output["invocation"]["system_id"].as_str(),
            Some("builtin_echo")
        );
    }

    #[test]
    fn lua_ask_agent_system_wrapper_accepts_nil_state() {
        let td = tempdir().expect("tempdir");
        let catalog_path = td.path().join("agents.json");
        fs::write(
            &catalog_path,
            r#"{
  "schema": "gentle.agent_systems.v1",
  "systems": [
    { "id": "builtin_echo", "label": "Built-in Echo", "transport": "builtin_echo" }
  ]
}"#,
        )
        .expect("write agent catalog");
        let out = LuaInterface::ask_agent_system(
            None,
            "builtin_echo".to_string(),
            "hello".to_string(),
            Some(catalog_path.to_string_lossy().to_string()),
            Some(false),
            Some(false),
            Some(vec![]),
            Some(true),
        )
        .expect("ask agent with nil state");
        assert!(!out.state_changed);
        assert_eq!(
            out.output["schema"].as_str(),
            Some("gentle.agent_ask_result.v1")
        );
        assert_eq!(
            out.output["invocation"]["system_id"].as_str(),
            Some("builtin_echo")
        );
    }
}
