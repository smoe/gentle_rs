use crate::{
    app::GENtleApp,
    dna_sequence::DNAsequence,
    engine::{Engine, EngineStateSummary, GentleEngine, Operation, ProjectState, Workflow},
    enzymes::active_restriction_enzymes,
    methylation_sites::MethylationMode,
    resource_sync, tf_motifs,
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
    resource_sync::sync_rebase(
        input,
        if output.trim().is_empty() {
            None
        } else {
            Some(output)
        },
        commercial_only,
    )
    .map_err(|e| deno_core::anyhow::anyhow!(e))
}

#[op2]
#[serde]
fn sync_jaspar_resource(
    #[string] input: &str,
    #[string] output: &str,
) -> Result<resource_sync::SyncReport, deno_core::anyhow::Error> {
    let report = resource_sync::sync_jaspar(
        input,
        if output.trim().is_empty() {
            None
        } else {
            Some(output)
        },
    )
    .map_err(|e| deno_core::anyhow::anyhow!(e))?;
    tf_motifs::reload();
    Ok(report)
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
        const APPLY_OPERATION: OpDecl = apply_operation();
        const APPLY_WORKFLOW: OpDecl = apply_workflow();
        const SYNC_REBASE_RESOURCE: OpDecl = sync_rebase_resource();
        const SYNC_JASPAR_RESOURCE: OpDecl = sync_jaspar_resource();
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
                APPLY_OPERATION,
                APPLY_WORKFLOW,
                SYNC_REBASE_RESOURCE,
                SYNC_JASPAR_RESOURCE,
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
          	function apply_operation(state, op) {
          		const payload = (typeof op === "string") ? op : JSON.stringify(op);
          		return Deno.core.ops.apply_operation(state, payload);
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
