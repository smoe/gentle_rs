use crate::{
    app::GENtleApp,
    dna_sequence::DNAsequence,
    engine::{Engine, EngineStateSummary, GentleEngine, Operation, ProjectState, Workflow},
    engine_shell::{execute_shell_command, ShellCommand},
    enzymes::active_restriction_enzymes,
    methylation_sites::MethylationMode,
    resource_sync,
};
use deno_core::*;
use serde::Serialize;

fn empty_to_none(value: &str) -> Option<&str> {
    let trimmed = value.trim();
    if trimmed.is_empty() {
        None
    } else {
        Some(trimmed)
    }
}

#[derive(Serialize)]
struct OperationApplyResponse {
    state: ProjectState,
    result: crate::engine::OpResult,
}

#[derive(Serialize)]
struct WorkflowApplyResponse {
    state: ProjectState,
    results: Vec<crate::engine::OpResult>,
}

#[derive(Serialize)]
struct GenomePreparedResponse {
    prepared: bool,
}

#[derive(Serialize)]
struct ShellUtilityApplyResponse {
    state: ProjectState,
    state_changed: bool,
    output: serde_json::Value,
}

fn sync_report_from_shell_output(
    output: serde_json::Value,
    context: &str,
) -> Result<resource_sync::SyncReport, deno_core::anyhow::Error> {
    let Some(report_value) = output.get("report").cloned() else {
        return Err(deno_core::anyhow::anyhow!(
            "{context}: missing 'report' field in shell output"
        ));
    };
    serde_json::from_value(report_value).map_err(|e| {
        deno_core::anyhow::anyhow!("{context}: could not parse sync report from shell output: {e}")
    })
}

#[op2]
#[serde]
fn load_dna(#[string] path: &str) -> Result<DNAsequence, deno_core::anyhow::Error> {
    let mut dna = GENtleApp::load_from_file(path)?;

    // Add default enzymes and stuff
    *dna.restriction_enzymes_mut() = active_restriction_enzymes();
    dna.set_max_restriction_enzyme_sites(Some(2));
    dna.set_methylation_mode(MethylationMode::both());
    dna.update_computed_features();
    Ok(dna)
}

#[op2]
#[serde]
fn sync_rebase_resource(
    #[string] input: &str,
    #[string] output: &str,
    commercial_only: bool,
) -> Result<resource_sync::SyncReport, deno_core::anyhow::Error> {
    let mut engine = GentleEngine::from_state(ProjectState::default());
    let command = ShellCommand::ResourcesSyncRebase {
        input: input.to_string(),
        output: empty_to_none(output).map(str::to_string),
        commercial_only,
    };
    let run = execute_shell_command(&mut engine, &command).map_err(|e| {
        deno_core::anyhow::anyhow!("resources sync-rebase failed: {e}")
    })?;
    sync_report_from_shell_output(run.output, "resources sync-rebase")
}

#[op2]
#[serde]
fn sync_jaspar_resource(
    #[string] input: &str,
    #[string] output: &str,
) -> Result<resource_sync::SyncReport, deno_core::anyhow::Error> {
    let mut engine = GentleEngine::from_state(ProjectState::default());
    let command = ShellCommand::ResourcesSyncJaspar {
        input: input.to_string(),
        output: empty_to_none(output).map(str::to_string),
    };
    let run = execute_shell_command(&mut engine, &command).map_err(|e| {
        deno_core::anyhow::anyhow!("resources sync-jaspar failed: {e}")
    })?;
    sync_report_from_shell_output(run.output, "resources sync-jaspar")
}

#[op2]
#[serde]
fn import_pool(
    #[serde] state: ProjectState,
    #[string] input: &str,
    #[string] prefix: &str,
) -> Result<ShellUtilityApplyResponse, deno_core::anyhow::Error> {
    let mut engine = GentleEngine::from_state(state);
    let prefix = if prefix.trim().is_empty() {
        "pool".to_string()
    } else {
        prefix.trim().to_string()
    };
    let command = ShellCommand::ImportPool {
        input: input.to_string(),
        prefix,
    };
    let run = execute_shell_command(&mut engine, &command)
        .map_err(|e| deno_core::anyhow::anyhow!("import-pool failed: {e}"))?;
    Ok(ShellUtilityApplyResponse {
        state: engine.state().clone(),
        state_changed: run.state_changed,
        output: run.output,
    })
}

#[op2]
fn write_gb(
    #[serde] seq: DNAsequence,
    #[string] path: &str,
) -> Result<(), deno_core::anyhow::Error> {
    seq.write_genbank_file(path)?;
    Ok(())
}

#[op2]
#[serde]
fn load_project(#[string] path: &str) -> Result<ProjectState, deno_core::anyhow::Error> {
    let state = ProjectState::load_from_path(path)?;
    Ok(state)
}

#[op2]
fn save_project(
    #[serde] state: ProjectState,
    #[string] path: &str,
) -> Result<(), deno_core::anyhow::Error> {
    state.save_to_path(path)?;
    Ok(())
}

#[op2]
#[serde]
fn capabilities() -> Result<crate::engine::Capabilities, deno_core::anyhow::Error> {
    Ok(GentleEngine::capabilities())
}

#[op2]
#[serde]
fn state_summary(
    #[serde] state: ProjectState,
) -> Result<EngineStateSummary, deno_core::anyhow::Error> {
    let engine = GentleEngine::from_state(state);
    Ok(engine.summarize_state())
}

#[op2]
#[serde]
fn inspect_dna_ladders(
    #[string] name_filter: &str,
) -> Result<crate::engine::DnaLadderCatalog, deno_core::anyhow::Error> {
    Ok(GentleEngine::inspect_dna_ladders(empty_to_none(
        name_filter,
    )))
}

#[op2]
#[serde]
fn inspect_rna_ladders(
    #[string] name_filter: &str,
) -> Result<crate::engine::RnaLadderCatalog, deno_core::anyhow::Error> {
    Ok(GentleEngine::inspect_rna_ladders(empty_to_none(
        name_filter,
    )))
}

#[op2]
#[serde]
fn export_dna_ladders(
    #[string] path: &str,
    #[string] name_filter: &str,
) -> Result<crate::engine::DnaLadderExportReport, deno_core::anyhow::Error> {
    GentleEngine::export_dna_ladders(path, empty_to_none(name_filter))
        .map_err(|e| deno_core::anyhow::anyhow!(e.to_string()))
}

#[op2]
#[serde]
fn export_rna_ladders(
    #[string] path: &str,
    #[string] name_filter: &str,
) -> Result<crate::engine::RnaLadderExportReport, deno_core::anyhow::Error> {
    GentleEngine::export_rna_ladders(path, empty_to_none(name_filter))
        .map_err(|e| deno_core::anyhow::anyhow!(e.to_string()))
}

#[op2]
#[serde]
fn list_reference_genomes(
    #[string] catalog_path: &str,
) -> Result<Vec<String>, deno_core::anyhow::Error> {
    GentleEngine::list_reference_genomes(empty_to_none(catalog_path))
        .map_err(|e| deno_core::anyhow::anyhow!(e.to_string()))
}

#[op2]
#[serde]
fn is_reference_genome_prepared(
    #[string] genome_id: &str,
    #[string] catalog_path: &str,
    #[string] cache_dir: &str,
) -> Result<GenomePreparedResponse, deno_core::anyhow::Error> {
    let prepared = GentleEngine::is_reference_genome_prepared(
        empty_to_none(catalog_path),
        genome_id,
        empty_to_none(cache_dir),
    )
    .map_err(|e| deno_core::anyhow::anyhow!(e.to_string()))?;
    Ok(GenomePreparedResponse { prepared })
}

#[op2]
#[serde]
fn list_reference_genome_genes(
    #[string] genome_id: &str,
    #[string] catalog_path: &str,
    #[string] cache_dir: &str,
) -> Result<Vec<crate::genomes::GenomeGeneRecord>, deno_core::anyhow::Error> {
    GentleEngine::list_reference_genome_genes(
        empty_to_none(catalog_path),
        genome_id,
        empty_to_none(cache_dir),
    )
    .map_err(|e| deno_core::anyhow::anyhow!(e.to_string()))
}

#[op2]
#[serde]
fn blast_reference_genome(
    #[string] genome_id: &str,
    #[string] query_sequence: &str,
    max_hits: u32,
    #[string] task: &str,
    #[string] catalog_path: &str,
    #[string] cache_dir: &str,
) -> Result<crate::genomes::GenomeBlastReport, deno_core::anyhow::Error> {
    GentleEngine::blast_reference_genome(
        empty_to_none(catalog_path),
        genome_id,
        query_sequence,
        max_hits.max(1) as usize,
        empty_to_none(task),
        empty_to_none(cache_dir),
    )
    .map_err(|e| deno_core::anyhow::anyhow!(e.to_string()))
}

#[op2]
#[serde]
fn blast_helper_genome(
    #[string] genome_id: &str,
    #[string] query_sequence: &str,
    max_hits: u32,
    #[string] task: &str,
    #[string] catalog_path: &str,
    #[string] cache_dir: &str,
) -> Result<crate::genomes::GenomeBlastReport, deno_core::anyhow::Error> {
    GentleEngine::blast_helper_genome(
        genome_id,
        query_sequence,
        max_hits.max(1) as usize,
        empty_to_none(task),
        empty_to_none(catalog_path),
        empty_to_none(cache_dir),
    )
    .map_err(|e| deno_core::anyhow::anyhow!(e.to_string()))
}

#[op2]
#[serde]
fn apply_operation(
    #[serde] state: ProjectState,
    #[string] op_json: &str,
) -> Result<OperationApplyResponse, deno_core::anyhow::Error> {
    let op: Operation = serde_json::from_str(op_json)?;
    let mut engine = GentleEngine::from_state(state);
    let result = engine.apply(op)?;
    Ok(OperationApplyResponse {
        state: engine.state().clone(),
        result,
    })
}

#[op2]
#[serde]
fn apply_workflow(
    #[serde] state: ProjectState,
    #[string] workflow_json: &str,
) -> Result<WorkflowApplyResponse, deno_core::anyhow::Error> {
    let workflow: Workflow = serde_json::from_str(workflow_json)?;
    let mut engine = GentleEngine::from_state(state);
    let results = engine.apply_workflow(workflow)?;
    Ok(WorkflowApplyResponse {
        state: engine.state().clone(),
        results,
    })
}

pub struct JavaScriptInterface {
    runtime: JsRuntime,
}

impl JavaScriptInterface {
    pub fn new() -> Self {
        // Build a deno_core::Extension providing custom ops
        const LOAD_DNA: OpDecl = load_dna();
        const WRITE_GB: OpDecl = write_gb();
        const LOAD_PROJECT: OpDecl = load_project();
        const SAVE_PROJECT: OpDecl = save_project();
        const CAPABILITIES: OpDecl = capabilities();
        const STATE_SUMMARY: OpDecl = state_summary();
        const INSPECT_DNA_LADDERS: OpDecl = inspect_dna_ladders();
        const INSPECT_RNA_LADDERS: OpDecl = inspect_rna_ladders();
        const EXPORT_DNA_LADDERS: OpDecl = export_dna_ladders();
        const EXPORT_RNA_LADDERS: OpDecl = export_rna_ladders();
        const LIST_REFERENCE_GENOMES: OpDecl = list_reference_genomes();
        const IS_REFERENCE_GENOME_PREPARED: OpDecl = is_reference_genome_prepared();
        const LIST_REFERENCE_GENOME_GENES: OpDecl = list_reference_genome_genes();
        const BLAST_REFERENCE_GENOME: OpDecl = blast_reference_genome();
        const BLAST_HELPER_GENOME: OpDecl = blast_helper_genome();
        const APPLY_OPERATION: OpDecl = apply_operation();
        const APPLY_WORKFLOW: OpDecl = apply_workflow();
        const SYNC_REBASE_RESOURCE: OpDecl = sync_rebase_resource();
        const SYNC_JASPAR_RESOURCE: OpDecl = sync_jaspar_resource();
        const IMPORT_POOL: OpDecl = import_pool();
        let ext = Extension {
            name: "my_ext",
            ops: std::borrow::Cow::Borrowed(&[
                LOAD_DNA,
                WRITE_GB,
                LOAD_PROJECT,
                SAVE_PROJECT,
                CAPABILITIES,
                STATE_SUMMARY,
                INSPECT_DNA_LADDERS,
                INSPECT_RNA_LADDERS,
                EXPORT_DNA_LADDERS,
                EXPORT_RNA_LADDERS,
                LIST_REFERENCE_GENOMES,
                IS_REFERENCE_GENOME_PREPARED,
                LIST_REFERENCE_GENOME_GENES,
                BLAST_REFERENCE_GENOME,
                BLAST_HELPER_GENOME,
                APPLY_OPERATION,
                APPLY_WORKFLOW,
                SYNC_REBASE_RESOURCE,
                SYNC_JASPAR_RESOURCE,
                IMPORT_POOL,
            ]),
            ..Default::default()
        };

        let mut ret = Self {
            runtime: JsRuntime::new(RuntimeOptions {
                extensions: vec![ext],
                ..Default::default()
            }),
        };
        let init_code = r#"
        	function load_dna(path) {return Deno.core.ops.load_dna(path)}
         	function write_gb(seq,path) {return Deno.core.ops.write_gb(seq,path)}
         	function load_project(path) {return Deno.core.ops.load_project(path)}
          	function save_project(state,path) {return Deno.core.ops.save_project(state,path)}
	          	function capabilities() {return Deno.core.ops.capabilities()}
	          	function state_summary(state) {return Deno.core.ops.state_summary(state)}
	          	function inspect_dna_ladders(name_filter) {
	          		return Deno.core.ops.inspect_dna_ladders(name_filter ?? "");
	          	}
		          	function list_dna_ladders(name_filter) {
		          		return inspect_dna_ladders(name_filter);
		          	}
		          	function inspect_rna_ladders(name_filter) {
		          		return Deno.core.ops.inspect_rna_ladders(name_filter ?? "");
		          	}
		          	function list_rna_ladders(name_filter) {
		          		return inspect_rna_ladders(name_filter);
		          	}
		          	function export_dna_ladders(path, name_filter) {
		          		return Deno.core.ops.export_dna_ladders(path, name_filter ?? "");
		          	}
		          	function export_rna_ladders(path, name_filter) {
		          		return Deno.core.ops.export_rna_ladders(path, name_filter ?? "");
		          	}
		          	function list_reference_genomes(catalog_path) {
		          		return Deno.core.ops.list_reference_genomes(catalog_path ?? "");
		          	}
          	function is_reference_genome_prepared(genome_id, catalog_path, cache_dir) {
          		const status = Deno.core.ops.is_reference_genome_prepared(genome_id, catalog_path ?? "", cache_dir ?? "");
          		return !!status.prepared;
          	}
          	function list_reference_genome_genes(genome_id, catalog_path, cache_dir) {
          		return Deno.core.ops.list_reference_genome_genes(genome_id, catalog_path ?? "", cache_dir ?? "");
          	}
          	function blast_reference_genome(genome_id, query_sequence, max_hits, task, catalog_path, cache_dir) {
          		return Deno.core.ops.blast_reference_genome(
          			genome_id,
          			query_sequence,
          			(max_hits === undefined ? 25 : max_hits),
          			task ?? "",
          			catalog_path ?? "",
          			cache_dir ?? ""
          		);
          	}
          	function blast_helper_genome(genome_id, query_sequence, max_hits, task, catalog_path, cache_dir) {
          		return Deno.core.ops.blast_helper_genome(
          			genome_id,
          			query_sequence,
          			(max_hits === undefined ? 25 : max_hits),
          			task ?? "",
          			catalog_path ?? "",
          			cache_dir ?? ""
          		);
          	}
	          	function apply_operation(state, op) {
	          		const payload = (typeof op === "string") ? op : JSON.stringify(op);
	          		return Deno.core.ops.apply_operation(state, payload);
	          	}
	          	function set_parameter(state, name, value) {
	          		return apply_operation(state, {
	          			SetParameter: {
	          				name: name,
	          				value: value
	          			}
	          		});
	          	}
	          	function set_vcf_display_filter(state, options) {
	          		const opts = options ?? {};
	          		let currentState = state;
	          		let lastResult = null;
	          		const applyOne = (name, value) => {
	          			const applied = set_parameter(currentState, name, value);
	          			currentState = applied.state;
	          			lastResult = applied.result;
	          		};
	          		if (opts.show_snp !== undefined) applyOne("vcf_display_show_snp", !!opts.show_snp);
	          		if (opts.show_ins !== undefined) applyOne("vcf_display_show_ins", !!opts.show_ins);
	          		if (opts.show_del !== undefined) applyOne("vcf_display_show_del", !!opts.show_del);
	          		if (opts.show_sv !== undefined) applyOne("vcf_display_show_sv", !!opts.show_sv);
	          		if (opts.show_other !== undefined) applyOne("vcf_display_show_other", !!opts.show_other);
	          		if (opts.pass_only !== undefined) applyOne("vcf_display_pass_only", !!opts.pass_only);
	          		if (opts.use_min_qual !== undefined) applyOne("vcf_display_use_min_qual", !!opts.use_min_qual);
	          		if (opts.min_qual !== undefined) applyOne("vcf_display_min_qual", Number(opts.min_qual));
	          		if (opts.use_max_qual !== undefined) applyOne("vcf_display_use_max_qual", !!opts.use_max_qual);
	          		if (opts.max_qual !== undefined) applyOne("vcf_display_max_qual", Number(opts.max_qual));
	          		if (opts.required_info_keys !== undefined) {
	          			applyOne("vcf_display_required_info_keys", opts.required_info_keys);
	          		}
	          		return { state: currentState, result: lastResult };
	          	}
	          	function apply_workflow(state, workflow) {
	          		const payload = (typeof workflow === "string") ? workflow : JSON.stringify(workflow);
	          		return Deno.core.ops.apply_workflow(state, payload);
	          	}
          	function sync_rebase(input, output, commercial_only) {
          		return Deno.core.ops.sync_rebase_resource(input, output ?? "", commercial_only ?? false);
          	}
          	function sync_jaspar(input, output) {
          		return Deno.core.ops.sync_jaspar_resource(input, output ?? "");
          	}
          	function import_pool(state, input, prefix) {
          		return Deno.core.ops.import_pool(state, input, prefix ?? "pool");
          	}
          	function digest(state, input, enzymes, output_id) {
          		const op = {
          			Digest: {
          				input: input,
          				enzymes: enzymes.split(",").map(s => s.trim()).filter(Boolean),
          				output_prefix: output_id ?? null
          			}
          		};
          		return apply_operation(state, op);
          	}
          	function prepare_genome(state, genome_id, catalog_path, cache_dir) {
          		return apply_operation(state, {
          			PrepareGenome: {
          				genome_id: genome_id,
          				catalog_path: catalog_path ?? null,
          				cache_dir: cache_dir ?? null
          			}
          		});
          	}
          	function extract_genome_region(state, genome_id, chromosome, start_1based, end_1based, output_id, catalog_path, cache_dir) {
          		return apply_operation(state, {
          			ExtractGenomeRegion: {
          				genome_id: genome_id,
          				chromosome: chromosome,
          				start_1based: start_1based,
          				end_1based: end_1based,
          				output_id: output_id ?? null,
          				catalog_path: catalog_path ?? null,
          				cache_dir: cache_dir ?? null
          			}
          		});
          	}
          	function extract_genome_gene(state, genome_id, gene_query, occurrence, output_id, catalog_path, cache_dir) {
          		return apply_operation(state, {
          			ExtractGenomeGene: {
          				genome_id: genome_id,
          				gene_query: gene_query,
          				occurrence: occurrence ?? null,
          				output_id: output_id ?? null,
          				catalog_path: catalog_path ?? null,
          				cache_dir: cache_dir ?? null
          			}
          		});
          	}
          	function extend_genome_anchor(state, seq_id, side, length_bp, output_id, catalog_path, cache_dir) {
          		const rawSide = String(side ?? "5p").trim().toLowerCase();
          		const sideValue = (rawSide === "3" || rawSide === "3p" || rawSide === "3prime" || rawSide === "3'" || rawSide === "three_prime" || rawSide === "three-prime")
          			? "three_prime"
          			: ((rawSide === "5" || rawSide === "5p" || rawSide === "5prime" || rawSide === "5'" || rawSide === "five_prime" || rawSide === "five-prime")
          				? "five_prime"
          				: null);
          		if (!sideValue) {
          			throw new Error("extend_genome_anchor side must be 5p or 3p");
          		}
          		return apply_operation(state, {
          			ExtendGenomeAnchor: {
          				seq_id: seq_id,
          				side: sideValue,
          				length_bp: Number(length_bp),
          				output_id: output_id ?? null,
          				catalog_path: catalog_path ?? null,
          				cache_dir: cache_dir ?? null
          			}
          		});
          	}
          	function import_genome_bed_track(state, seq_id, path, track_name, min_score, max_score, clear_existing) {
          		return apply_operation(state, {
          			ImportGenomeBedTrack: {
          				seq_id: seq_id,
          				path: path,
          				track_name: track_name ?? null,
          				min_score: (min_score === undefined ? null : min_score),
          				max_score: (max_score === undefined ? null : max_score),
          				clear_existing: (clear_existing === undefined ? null : !!clear_existing)
          			}
          		});
          	}
          	function import_genome_bigwig_track(state, seq_id, path, track_name, min_score, max_score, clear_existing) {
          		return apply_operation(state, {
          			ImportGenomeBigWigTrack: {
          				seq_id: seq_id,
          				path: path,
          				track_name: track_name ?? null,
          				min_score: (min_score === undefined ? null : min_score),
          				max_score: (max_score === undefined ? null : max_score),
          				clear_existing: (clear_existing === undefined ? null : !!clear_existing)
          			}
          		});
          	}
          	function import_genome_vcf_track(state, seq_id, path, track_name, min_score, max_score, clear_existing) {
          		return apply_operation(state, {
          			ImportGenomeVcfTrack: {
          				seq_id: seq_id,
          				path: path,
          				track_name: track_name ?? null,
          				min_score: (min_score === undefined ? null : min_score),
          				max_score: (max_score === undefined ? null : max_score),
          				clear_existing: (clear_existing === undefined ? null : !!clear_existing)
          			}
          		});
          	}
          	function render_pool_gel_svg(state, inputs, path, ladders) {
          		const seqIds = (Array.isArray(inputs) ? inputs : String(inputs).split(","))
          			.map(s => String(s).trim())
          			.filter(Boolean);
          		const ladderIds = (ladders == null)
          			? null
          			: (Array.isArray(ladders) ? ladders : String(ladders).split(","))
          				.map(s => String(s).trim())
          				.filter(Boolean);
          		return apply_operation(state, {
          			RenderPoolGelSvg: {
          				inputs: seqIds,
          				path: path,
          				ladders: (ladderIds && ladderIds.length > 0) ? ladderIds : null
          			}
          		});
          	}
        "#
        .to_string();
        ret.run(init_code);
        println!("Interactive JavaScript Shell (type 'exit' to quit)");
        ret
    }

    pub fn run_checked(&mut self, code: String) -> Result<(), String> {
        self.runtime
            .execute_script("<usage>", code)
            .map(|_| ())
            .map_err(|e| e.to_string())
    }

    pub fn run(&mut self, code: String) {
        if let Err(e) = self.run_checked(code) {
            eprintln!("{e}");
        }
    }
}

impl Default for JavaScriptInterface {
    fn default() -> Self {
        Self::new()
    }
}
