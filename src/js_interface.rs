//! JavaScript adapter wrappers over shared engine contracts.

use crate::{
    app::GENtleApp,
    dna_sequence::DNAsequence,
    engine::{
        Engine, EngineStateSummary, FeatureExpertTarget, GentleEngine, Operation, ProjectState,
        Workflow,
    },
    engine_shell::{ShellCommand, execute_shell_command},
    enzymes::active_restriction_enzymes,
    methylation_sites::MethylationMode,
    resource_sync,
};
use deno_core::*;
use serde::Serialize;
use std::borrow::Cow;

#[derive(Debug)]
struct JsAnyhow(deno_core::anyhow::Error);

impl From<deno_core::anyhow::Error> for JsAnyhow {
    fn from(value: deno_core::anyhow::Error) -> Self {
        Self(value)
    }
}

impl std::fmt::Display for JsAnyhow {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        std::fmt::Display::fmt(&self.0, f)
    }
}

impl std::error::Error for JsAnyhow {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        self.0.source()
    }
}

impl deno_error::JsErrorClass for JsAnyhow {
    fn get_class(&self) -> Cow<'static, str> {
        Cow::Borrowed("Error")
    }

    fn get_message(&self) -> Cow<'static, str> {
        Cow::Owned(self.0.to_string())
    }

    fn get_additional_properties(&self) -> deno_error::AdditionalProperties {
        Box::new(std::iter::empty())
    }

    fn get_ref(&self) -> &(dyn std::error::Error + Send + Sync + 'static) {
        self
    }
}

fn empty_to_none(value: &str) -> Option<&str> {
    let trimmed = value.trim();
    if trimmed.is_empty() {
        None
    } else {
        Some(trimmed)
    }
}

fn compose_blast_request_override(
    max_hits: u32,
    task: &str,
    options_json: &str,
) -> Result<serde_json::Value, JsAnyhow> {
    let mut request = serde_json::Map::new();
    request.insert(
        "max_hits".to_string(),
        serde_json::Value::from(max_hits.max(1) as usize),
    );
    if let Some(task) = empty_to_none(task) {
        request.insert("task".to_string(), serde_json::Value::from(task));
    }
    if let Some(raw) = empty_to_none(options_json) {
        let parsed: serde_json::Value = serde_json::from_str(raw)
            .map_err(|e| deno_core::anyhow::anyhow!("Invalid BLAST options_json payload: {e}"))?;
        let obj = parsed.as_object().ok_or_else(|| {
            deno_core::anyhow::anyhow!("BLAST options_json payload must decode to a JSON object")
        })?;
        for (k, v) in obj {
            request.insert(k.clone(), v.clone());
        }
    }
    Ok(serde_json::Value::Object(request))
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
) -> Result<resource_sync::SyncReport, JsAnyhow> {
    let Some(report_value) = output.get("report").cloned() else {
        return Err(deno_core::anyhow::anyhow!(
            "{context}: missing 'report' field in shell output"
        )
        .into());
    };
    serde_json::from_value(report_value)
        .map_err(|e| {
            deno_core::anyhow::anyhow!(
                "{context}: could not parse sync report from shell output: {e}"
            )
        })
        .map_err(Into::into)
}

fn sync_rebase_resource_impl(
    input: &str,
    output: &str,
    commercial_only: bool,
) -> Result<resource_sync::SyncReport, JsAnyhow> {
    let mut engine = GentleEngine::from_state(ProjectState::default());
    let command = ShellCommand::ResourcesSyncRebase {
        input: input.to_string(),
        output: empty_to_none(output).map(str::to_string),
        commercial_only,
    };
    let run = execute_shell_command(&mut engine, &command)
        .map_err(|e| deno_core::anyhow::anyhow!("resources sync-rebase failed: {e}"))?;
    sync_report_from_shell_output(run.output, "resources sync-rebase")
}

fn sync_jaspar_resource_impl(
    input: &str,
    output: &str,
) -> Result<resource_sync::SyncReport, JsAnyhow> {
    let mut engine = GentleEngine::from_state(ProjectState::default());
    let command = ShellCommand::ResourcesSyncJaspar {
        input: input.to_string(),
        output: empty_to_none(output).map(str::to_string),
    };
    let run = execute_shell_command(&mut engine, &command)
        .map_err(|e| deno_core::anyhow::anyhow!("resources sync-jaspar failed: {e}"))?;
    sync_report_from_shell_output(run.output, "resources sync-jaspar")
}

fn import_pool_impl(
    state: ProjectState,
    input: &str,
    prefix: &str,
) -> Result<ShellUtilityApplyResponse, JsAnyhow> {
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

fn list_agent_systems_impl(catalog_path: &str) -> Result<serde_json::Value, JsAnyhow> {
    let mut engine = GentleEngine::from_state(ProjectState::default());
    let command = ShellCommand::AgentsList {
        catalog_path: empty_to_none(catalog_path).map(str::to_string),
    };
    let run = execute_shell_command(&mut engine, &command)
        .map_err(|e| deno_core::anyhow::anyhow!("agents list failed: {e}"))?;
    Ok(run.output)
}

fn list_construct_reasoning_graphs_impl(
    state: ProjectState,
    seq_id: &str,
) -> Result<serde_json::Value, JsAnyhow> {
    let mut engine = GentleEngine::from_state(state);
    let command = ShellCommand::ConstructReasoningListGraphs {
        seq_id: empty_to_none(seq_id).map(str::to_string),
    };
    let run = execute_shell_command(&mut engine, &command)
        .map_err(|e| deno_core::anyhow::anyhow!("construct-reasoning list-graphs failed: {e}"))?;
    Ok(run.output)
}

fn show_construct_reasoning_graph_impl(
    state: ProjectState,
    graph_id: &str,
) -> Result<serde_json::Value, JsAnyhow> {
    let mut engine = GentleEngine::from_state(state);
    let command = ShellCommand::ConstructReasoningShowGraph {
        graph_id: graph_id.to_string(),
    };
    let run = execute_shell_command(&mut engine, &command)
        .map_err(|e| deno_core::anyhow::anyhow!("construct-reasoning show-graph failed: {e}"))?;
    Ok(run.output)
}

fn set_construct_reasoning_annotation_status_impl(
    state: ProjectState,
    graph_id: &str,
    annotation_id: &str,
    editable_status: &str,
) -> Result<ShellUtilityApplyResponse, JsAnyhow> {
    let graph_id = graph_id.trim().to_string();
    let annotation_id = annotation_id.trim().to_string();
    let editable_status = match editable_status.trim().to_ascii_lowercase().as_str() {
        "draft" => crate::engine::EditableStatus::Draft,
        "accepted" => crate::engine::EditableStatus::Accepted,
        "rejected" => crate::engine::EditableStatus::Rejected,
        "locked" => crate::engine::EditableStatus::Locked,
        other => {
            return Err(deno_core::anyhow::anyhow!(
                "Unsupported construct-reasoning annotation status '{other}' (expected draft|accepted|rejected|locked)"
            )
            .into())
        }
    };
    let mut engine = GentleEngine::from_state(state);
    let command = ShellCommand::ConstructReasoningSetAnnotationStatus {
        graph_id,
        annotation_id,
        editable_status,
    };
    let run = execute_shell_command(&mut engine, &command).map_err(|e| {
        deno_core::anyhow::anyhow!("construct-reasoning set-annotation-status failed: {e}")
    })?;
    Ok(ShellUtilityApplyResponse {
        state: engine.state().clone(),
        state_changed: run.state_changed,
        output: run.output,
    })
}

fn write_back_construct_reasoning_annotation_impl(
    state: ProjectState,
    graph_id: &str,
    annotation_id: &str,
) -> Result<ShellUtilityApplyResponse, JsAnyhow> {
    let mut engine = GentleEngine::from_state(state);
    let command = ShellCommand::ConstructReasoningWriteAnnotation {
        graph_id: graph_id.trim().to_string(),
        annotation_id: annotation_id.trim().to_string(),
    };
    let run = execute_shell_command(&mut engine, &command).map_err(|e| {
        deno_core::anyhow::anyhow!("construct-reasoning write-annotation failed: {e}")
    })?;
    Ok(ShellUtilityApplyResponse {
        state: engine.state().clone(),
        state_changed: run.state_changed,
        output: run.output,
    })
}

#[allow(clippy::too_many_arguments)]
fn ask_agent_system_impl(
    state: ProjectState,
    system_id: &str,
    prompt: &str,
    catalog_path: &str,
    allow_auto_exec: bool,
    execute_all: bool,
    execute_indices: &[u32],
    include_state_summary: bool,
) -> Result<ShellUtilityApplyResponse, JsAnyhow> {
    let mut engine = GentleEngine::from_state(state);
    let command = ShellCommand::AgentsAsk {
        system_id: system_id.to_string(),
        prompt: prompt.to_string(),
        catalog_path: empty_to_none(catalog_path).map(str::to_string),
        base_url_override: None,
        model_override: None,
        timeout_seconds: None,
        connect_timeout_seconds: None,
        read_timeout_seconds: None,
        max_retries: None,
        max_response_bytes: None,
        include_state_summary,
        allow_auto_exec,
        execute_all,
        execute_indices: execute_indices.iter().map(|idx| *idx as usize).collect(),
    };
    let run = execute_shell_command(&mut engine, &command)
        .map_err(|e| deno_core::anyhow::anyhow!("agents ask failed: {e}"))?;
    Ok(ShellUtilityApplyResponse {
        state: engine.state().clone(),
        state_changed: run.state_changed,
        output: run.output,
    })
}

#[op2]
#[serde]
fn load_dna(#[string] path: &str) -> Result<DNAsequence, JsAnyhow> {
    let mut dna = GENtleApp::load_from_file(path).map_err(deno_core::anyhow::Error::from)?;

    // Add default enzymes and stuff
    *dna.restriction_enzymes_mut() = active_restriction_enzymes();
    dna.set_max_restriction_enzyme_sites(None);
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
) -> Result<resource_sync::SyncReport, JsAnyhow> {
    sync_rebase_resource_impl(input, output, commercial_only)
}

#[op2]
#[serde]
fn sync_jaspar_resource(
    #[string] input: &str,
    #[string] output: &str,
) -> Result<resource_sync::SyncReport, JsAnyhow> {
    sync_jaspar_resource_impl(input, output)
}

#[op2]
#[serde]
fn import_pool(
    #[serde] state: ProjectState,
    #[string] input: &str,
    #[string] prefix: &str,
) -> Result<ShellUtilityApplyResponse, JsAnyhow> {
    import_pool_impl(state, input, prefix)
}

#[op2]
fn write_gb(#[serde] seq: DNAsequence, #[string] path: &str) -> Result<(), JsAnyhow> {
    seq.write_genbank_file(path)
        .map_err(deno_core::anyhow::Error::from)?;
    Ok(())
}

#[op2]
#[serde]
fn load_project(#[string] path: &str) -> Result<ProjectState, JsAnyhow> {
    let state = ProjectState::load_from_path(path).map_err(deno_core::anyhow::Error::from)?;
    Ok(state)
}

#[op2]
fn save_project(#[serde] state: ProjectState, #[string] path: &str) -> Result<(), JsAnyhow> {
    state
        .save_to_path(path)
        .map_err(deno_core::anyhow::Error::from)?;
    Ok(())
}

#[op2]
#[serde]
fn capabilities() -> Result<crate::engine::Capabilities, JsAnyhow> {
    Ok(GentleEngine::capabilities())
}

#[op2]
#[serde]
fn state_summary(#[serde] state: ProjectState) -> Result<EngineStateSummary, JsAnyhow> {
    let engine = GentleEngine::from_state(state);
    Ok(engine.summarize_state())
}

#[op2]
#[serde]
fn inspect_dna_ladders(
    #[string] name_filter: &str,
) -> Result<crate::engine::DnaLadderCatalog, JsAnyhow> {
    Ok(GentleEngine::inspect_dna_ladders(empty_to_none(
        name_filter,
    )))
}

#[op2]
#[serde]
fn inspect_rna_ladders(
    #[string] name_filter: &str,
) -> Result<crate::engine::RnaLadderCatalog, JsAnyhow> {
    Ok(GentleEngine::inspect_rna_ladders(empty_to_none(
        name_filter,
    )))
}

#[op2]
#[serde]
fn export_dna_ladders(
    #[string] path: &str,
    #[string] name_filter: &str,
) -> Result<crate::engine::DnaLadderExportReport, JsAnyhow> {
    GentleEngine::export_dna_ladders(path, empty_to_none(name_filter))
        .map_err(|e| deno_core::anyhow::anyhow!(e.to_string()))
        .map_err(Into::into)
}

#[op2]
#[serde]
fn export_rna_ladders(
    #[string] path: &str,
    #[string] name_filter: &str,
) -> Result<crate::engine::RnaLadderExportReport, JsAnyhow> {
    GentleEngine::export_rna_ladders(path, empty_to_none(name_filter))
        .map_err(|e| deno_core::anyhow::anyhow!(e.to_string()))
        .map_err(Into::into)
}

#[op2]
#[serde]
fn list_reference_genomes(#[string] catalog_path: &str) -> Result<Vec<String>, JsAnyhow> {
    GentleEngine::list_reference_genomes(empty_to_none(catalog_path))
        .map_err(|e| deno_core::anyhow::anyhow!(e.to_string()))
        .map_err(Into::into)
}

#[op2]
#[serde]
fn list_reference_catalog_entries(
    #[string] catalog_path: &str,
    #[string] filter: &str,
) -> Result<Vec<crate::genomes::GenomeCatalogListEntry>, JsAnyhow> {
    GentleEngine::list_reference_catalog_entries(empty_to_none(catalog_path), empty_to_none(filter))
        .map_err(|e| deno_core::anyhow::anyhow!(e.to_string()))
        .map_err(Into::into)
}

#[op2]
#[serde]
fn list_helper_catalog_entries(
    #[string] catalog_path: &str,
    #[string] filter: &str,
) -> Result<Vec<crate::genomes::GenomeCatalogListEntry>, JsAnyhow> {
    GentleEngine::list_helper_catalog_entries(empty_to_none(catalog_path), empty_to_none(filter))
        .map_err(|e| deno_core::anyhow::anyhow!(e.to_string()))
        .map_err(Into::into)
}

#[op2]
#[serde]
fn list_host_profile_catalog_entries(
    #[string] catalog_path: &str,
    #[string] filter: &str,
) -> Result<Vec<crate::engine::HostProfileRecord>, JsAnyhow> {
    GentleEngine::list_host_profile_catalog_entries(
        empty_to_none(catalog_path),
        empty_to_none(filter),
    )
    .map_err(|e| deno_core::anyhow::anyhow!(e.to_string()))
    .map_err(Into::into)
}

#[op2]
#[serde]
fn list_ensembl_installable_genomes(
    #[string] collection: &str,
    #[string] filter: &str,
) -> Result<crate::genomes::EnsemblInstallableGenomeCatalog, JsAnyhow> {
    GentleEngine::discover_ensembl_installable_genomes(
        empty_to_none(collection),
        empty_to_none(filter),
    )
    .map_err(|e| deno_core::anyhow::anyhow!(e.to_string()))
    .map_err(Into::into)
}

#[op2]
#[serde]
fn list_construct_reasoning_graphs(
    #[serde] state: ProjectState,
    #[string] seq_id: &str,
) -> Result<serde_json::Value, JsAnyhow> {
    list_construct_reasoning_graphs_impl(state, seq_id)
}

#[op2]
#[serde]
fn show_construct_reasoning_graph(
    #[serde] state: ProjectState,
    #[string] graph_id: &str,
) -> Result<serde_json::Value, JsAnyhow> {
    show_construct_reasoning_graph_impl(state, graph_id)
}

#[op2]
#[serde]
fn set_construct_reasoning_annotation_status(
    #[serde] state: ProjectState,
    #[string] graph_id: &str,
    #[string] annotation_id: &str,
    #[string] editable_status: &str,
) -> Result<ShellUtilityApplyResponse, JsAnyhow> {
    set_construct_reasoning_annotation_status_impl(state, graph_id, annotation_id, editable_status)
}

#[op2]
#[serde]
fn write_back_construct_reasoning_annotation(
    #[serde] state: ProjectState,
    #[string] graph_id: &str,
    #[string] annotation_id: &str,
) -> Result<ShellUtilityApplyResponse, JsAnyhow> {
    write_back_construct_reasoning_annotation_impl(state, graph_id, annotation_id)
}

#[op2]
#[serde]
fn list_agent_systems(#[string] catalog_path: &str) -> Result<serde_json::Value, JsAnyhow> {
    list_agent_systems_impl(catalog_path)
}

#[op2]
#[serde]
fn ask_agent_system(
    #[serde] state: Option<ProjectState>,
    #[string] system_id: &str,
    #[string] prompt: &str,
    #[string] catalog_path: &str,
    allow_auto_exec: bool,
    execute_all: bool,
    #[serde] execute_indices: Vec<u32>,
    include_state_summary: bool,
) -> Result<ShellUtilityApplyResponse, JsAnyhow> {
    ask_agent_system_impl(
        state.unwrap_or_default(),
        system_id,
        prompt,
        catalog_path,
        allow_auto_exec,
        execute_all,
        &execute_indices,
        include_state_summary,
    )
}

#[op2]
#[serde]
fn is_reference_genome_prepared(
    #[string] genome_id: &str,
    #[string] catalog_path: &str,
    #[string] cache_dir: &str,
) -> Result<GenomePreparedResponse, JsAnyhow> {
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
) -> Result<Vec<crate::genomes::GenomeGeneRecord>, JsAnyhow> {
    GentleEngine::list_reference_genome_genes(
        empty_to_none(catalog_path),
        genome_id,
        empty_to_none(cache_dir),
    )
    .map_err(|e| deno_core::anyhow::anyhow!(e.to_string()))
    .map_err(Into::into)
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
    #[string] options_json: &str,
) -> Result<crate::genomes::GenomeBlastReport, JsAnyhow> {
    let request = compose_blast_request_override(max_hits, task, options_json)?;
    GentleEngine::blast_reference_genome_with_request_options(
        empty_to_none(catalog_path),
        genome_id,
        query_sequence,
        Some(&request),
        empty_to_none(cache_dir),
    )
    .map_err(|e| deno_core::anyhow::anyhow!(e.to_string()))
    .map_err(Into::into)
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
    #[string] options_json: &str,
) -> Result<crate::genomes::GenomeBlastReport, JsAnyhow> {
    let request = compose_blast_request_override(max_hits, task, options_json)?;
    GentleEngine::blast_helper_genome_with_request_options(
        genome_id,
        query_sequence,
        Some(&request),
        empty_to_none(catalog_path),
        empty_to_none(cache_dir),
    )
    .map_err(|e| deno_core::anyhow::anyhow!(e.to_string()))
    .map_err(Into::into)
}

#[op2]
#[serde]
fn apply_operation(
    #[serde] state: ProjectState,
    #[string] op_json: &str,
) -> Result<OperationApplyResponse, JsAnyhow> {
    apply_operation_impl(state, op_json)
}

fn apply_operation_impl(
    state: ProjectState,
    op_json: &str,
) -> Result<OperationApplyResponse, JsAnyhow> {
    let op: Operation = serde_json::from_str(op_json).map_err(deno_core::anyhow::Error::from)?;
    let mut engine = GentleEngine::from_state(state);
    let result = engine.apply(op).map_err(deno_core::anyhow::Error::from)?;
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
) -> Result<WorkflowApplyResponse, JsAnyhow> {
    let workflow: Workflow =
        serde_json::from_str(workflow_json).map_err(deno_core::anyhow::Error::from)?;
    let mut engine = GentleEngine::from_state(state);
    let results = engine
        .apply_workflow(workflow)
        .map_err(deno_core::anyhow::Error::from)?;
    Ok(WorkflowApplyResponse {
        state: engine.state().clone(),
        results,
    })
}

#[op2]
#[serde]
fn inspect_feature_expert(
    #[serde] state: ProjectState,
    #[string] seq_id: &str,
    #[string] target_json: &str,
) -> Result<serde_json::Value, JsAnyhow> {
    let target: FeatureExpertTarget =
        serde_json::from_str(target_json).map_err(deno_core::anyhow::Error::from)?;
    let engine = GentleEngine::from_state(state);
    let view = engine
        .inspect_feature_expert(seq_id, &target)
        .map_err(deno_core::anyhow::Error::from)?;
    serde_json::to_value(view)
        .map_err(|e| deno_core::anyhow::anyhow!(e.to_string()))
        .map_err(Into::into)
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
        const LIST_REFERENCE_CATALOG_ENTRIES: OpDecl = list_reference_catalog_entries();
        const LIST_HELPER_CATALOG_ENTRIES: OpDecl = list_helper_catalog_entries();
        const LIST_HOST_PROFILE_CATALOG_ENTRIES: OpDecl = list_host_profile_catalog_entries();
        const LIST_ENSEMBL_INSTALLABLE_GENOMES: OpDecl = list_ensembl_installable_genomes();
        const LIST_CONSTRUCT_REASONING_GRAPHS: OpDecl = list_construct_reasoning_graphs();
        const SHOW_CONSTRUCT_REASONING_GRAPH: OpDecl = show_construct_reasoning_graph();
        const SET_CONSTRUCT_REASONING_ANNOTATION_STATUS: OpDecl =
            set_construct_reasoning_annotation_status();
        const WRITE_BACK_CONSTRUCT_REASONING_ANNOTATION: OpDecl =
            write_back_construct_reasoning_annotation();
        const LIST_AGENT_SYSTEMS: OpDecl = list_agent_systems();
        const ASK_AGENT_SYSTEM: OpDecl = ask_agent_system();
        const IS_REFERENCE_GENOME_PREPARED: OpDecl = is_reference_genome_prepared();
        const LIST_REFERENCE_GENOME_GENES: OpDecl = list_reference_genome_genes();
        const BLAST_REFERENCE_GENOME: OpDecl = blast_reference_genome();
        const BLAST_HELPER_GENOME: OpDecl = blast_helper_genome();
        const APPLY_OPERATION: OpDecl = apply_operation();
        const APPLY_WORKFLOW: OpDecl = apply_workflow();
        const INSPECT_FEATURE_EXPERT: OpDecl = inspect_feature_expert();
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
                LIST_REFERENCE_CATALOG_ENTRIES,
                LIST_HELPER_CATALOG_ENTRIES,
                LIST_HOST_PROFILE_CATALOG_ENTRIES,
                LIST_ENSEMBL_INSTALLABLE_GENOMES,
                LIST_CONSTRUCT_REASONING_GRAPHS,
                SHOW_CONSTRUCT_REASONING_GRAPH,
                SET_CONSTRUCT_REASONING_ANNOTATION_STATUS,
                WRITE_BACK_CONSTRUCT_REASONING_ANNOTATION,
                LIST_AGENT_SYSTEMS,
                ASK_AGENT_SYSTEM,
                IS_REFERENCE_GENOME_PREPARED,
                LIST_REFERENCE_GENOME_GENES,
                BLAST_REFERENCE_GENOME,
                BLAST_HELPER_GENOME,
                APPLY_OPERATION,
                APPLY_WORKFLOW,
                INSPECT_FEATURE_EXPERT,
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
                      function list_reference_catalog_entries(catalog_path, filter) {
                        return Deno.core.ops.list_reference_catalog_entries(catalog_path ?? "", filter ?? "");
                      }
                      function list_helper_catalog_entries(catalog_path, filter) {
                        return Deno.core.ops.list_helper_catalog_entries(catalog_path ?? "", filter ?? "");
                      }
                      function list_host_profile_catalog_entries(catalog_path, filter) {
                        return Deno.core.ops.list_host_profile_catalog_entries(catalog_path ?? "", filter ?? "");
                      }
                      function list_ensembl_installable_genomes(collection, filter) {
                        return Deno.core.ops.list_ensembl_installable_genomes(collection ?? "", filter ?? "");
                      }
                      function list_construct_reasoning_graphs(state, seq_id) {
                        return Deno.core.ops.list_construct_reasoning_graphs(state, seq_id ?? "");
                      }
                      function show_construct_reasoning_graph(state, graph_id) {
                        return Deno.core.ops.show_construct_reasoning_graph(state, graph_id);
                      }
                      function set_construct_reasoning_annotation_status(state, graph_id, annotation_id, editable_status) {
                        return Deno.core.ops.set_construct_reasoning_annotation_status(
                          state,
                          graph_id,
                          annotation_id,
                          editable_status
                        );
                      }
                      function write_back_construct_reasoning_annotation(state, graph_id, annotation_id) {
                        return Deno.core.ops.write_back_construct_reasoning_annotation(
                          state,
                          graph_id,
                          annotation_id
                        );
                      }
                      function list_agent_systems(catalog_path) {
                        return Deno.core.ops.list_agent_systems(catalog_path ?? "");
                      }
			          	function ask_agent_system(state, system_id, prompt, options) {
			          		const opts = options ?? {};
			          		return Deno.core.ops.ask_agent_system(
			          			(state === undefined ? null : state),
			          			system_id,
			          			prompt,
			          			opts.catalog_path ?? "",
			          			!!opts.allow_auto_exec,
			          			!!opts.execute_all,
			          			Array.isArray(opts.execute_indices) ? opts.execute_indices : [],
			          			(opts.include_state_summary === undefined) ? true : !!opts.include_state_summary
			          		);
			          	}
	          	function is_reference_genome_prepared(genome_id, catalog_path, cache_dir) {
	          		const status = Deno.core.ops.is_reference_genome_prepared(genome_id, catalog_path ?? "", cache_dir ?? "");
	          		return !!status.prepared;
	          	}
          	function list_reference_genome_genes(genome_id, catalog_path, cache_dir) {
          		return Deno.core.ops.list_reference_genome_genes(genome_id, catalog_path ?? "", cache_dir ?? "");
          	}
          	function blast_reference_genome(genome_id, query_sequence, max_hits, task, catalog_path, cache_dir, options_json) {
          		return Deno.core.ops.blast_reference_genome(
          			genome_id,
          			query_sequence,
          			(max_hits === undefined ? 25 : max_hits),
          			task ?? "",
          			catalog_path ?? "",
          			cache_dir ?? "",
          			options_json ?? ""
          		);
          	}
          	function blast_helper_genome(genome_id, query_sequence, max_hits, task, catalog_path, cache_dir, options_json) {
          		return Deno.core.ops.blast_helper_genome(
          			genome_id,
          			query_sequence,
          			(max_hits === undefined ? 25 : max_hits),
          			task ?? "",
          			catalog_path ?? "",
          			cache_dir ?? "",
          			options_json ?? ""
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
          	function inspect_feature_expert(state, seq_id, target) {
          		const payload = (typeof target === "string") ? target : JSON.stringify(target);
          		return Deno.core.ops.inspect_feature_expert(state, seq_id, payload);
          	}
          	function render_feature_expert_svg(state, seq_id, target, path) {
          		return apply_operation(state, {
          			RenderFeatureExpertSvg: {
          				seq_id: seq_id,
          				target: (typeof target === "string") ? JSON.parse(target) : target,
          				path: path
          			}
          		});
          	}
          	function render_dotplot_svg(state, seq_id, dotplot_id, path, options) {
          		const opts = options ?? {};
          		const flexTrack = opts.flex_track_id ?? opts.flexTrackId ?? null;
          		let density = null;
          		if (opts.display_density_threshold !== undefined && opts.display_density_threshold !== null) {
          			const parsed = Number(opts.display_density_threshold);
          			if (!Number.isFinite(parsed)) {
          				throw new Error("render_dotplot_svg options.display_density_threshold must be numeric");
          			}
          			density = parsed;
          		}
          		let gain = null;
          		if (opts.display_intensity_gain !== undefined && opts.display_intensity_gain !== null) {
          			const parsed = Number(opts.display_intensity_gain);
          			if (!Number.isFinite(parsed)) {
          				throw new Error("render_dotplot_svg options.display_intensity_gain must be numeric");
          			}
          			gain = parsed;
          		}
          		return apply_operation(state, {
          			RenderDotplotSvg: {
          				seq_id: seq_id,
          				dotplot_id: dotplot_id,
          				path: path,
          				flex_track_id: flexTrack,
          				display_density_threshold: density,
          				display_intensity_gain: gain
          			}
          		});
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
          	function extract_genome_region(state, genome_id, chromosome, start_1based, end_1based, output_id, catalog_path, cache_dir, annotation_scope, max_annotation_features) {
          		return apply_operation(state, {
          			ExtractGenomeRegion: {
          				genome_id: genome_id,
          				chromosome: chromosome,
          				start_1based: start_1based,
          				end_1based: end_1based,
          				output_id: output_id ?? null,
          				annotation_scope: annotation_scope ?? null,
          				max_annotation_features: (max_annotation_features === undefined || max_annotation_features === null) ? null : Number(max_annotation_features),
          				catalog_path: catalog_path ?? null,
          				cache_dir: cache_dir ?? null
          			}
          		});
          	}
          	function extract_genome_gene(state, genome_id, gene_query, occurrence, output_id, catalog_path, cache_dir, annotation_scope, max_annotation_features, extract_mode, promoter_upstream_bp) {
          		return apply_operation(state, {
          			ExtractGenomeGene: {
          				genome_id: genome_id,
          				gene_query: gene_query,
          				occurrence: occurrence ?? null,
          				output_id: output_id ?? null,
          				extract_mode: extract_mode ?? null,
          				promoter_upstream_bp: (promoter_upstream_bp === undefined || promoter_upstream_bp === null) ? null : Number(promoter_upstream_bp),
          				annotation_scope: annotation_scope ?? null,
          				max_annotation_features: (max_annotation_features === undefined || max_annotation_features === null) ? null : Number(max_annotation_features),
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

#[cfg(test)]
mod tests {
    use super::*;
    use crate::dna_sequence::DNAsequence;
    use crate::engine::{
        AdapterCaptureProtectionMode, AdapterCaptureStyle, AdapterRestrictionCapturePlan,
        ConstructObjective, GentleEngine, ProjectState,
    };
    use crate::engine_shell::execute_shell_command;
    use crate::test_support::{
        decision_trace_fixture_state, write_demo_jaspar_pfm, write_demo_pool_json,
        write_demo_rebase_withrefm,
    };
    use std::fs;
    use tempfile::tempdir;

    fn normalize_construct_reasoning_status_output(
        mut value: serde_json::Value,
    ) -> serde_json::Value {
        if let Some(graph) = value
            .get_mut("graph")
            .and_then(serde_json::Value::as_object_mut)
        {
            graph.remove("generated_at_unix_ms");
        }
        value
    }

    #[test]
    fn js_sync_rebase_resource_wrapper_writes_snapshot() {
        let td = tempdir().expect("tempdir");
        let input_path = write_demo_rebase_withrefm(td.path());
        let output_path = td.path().join("rebase.json");
        let report = sync_rebase_resource_impl(
            input_path.to_string_lossy().as_ref(),
            output_path.to_string_lossy().as_ref(),
            true,
        )
        .expect("sync rebase");
        assert_eq!(report.resource, "rebase-commercial");
        assert_eq!(report.item_count, 1);
        assert!(output_path.exists());
    }

    #[test]
    fn js_sync_jaspar_resource_wrapper_writes_snapshot() {
        let td = tempdir().expect("tempdir");
        let input_path = write_demo_jaspar_pfm(td.path());
        let output_path = td.path().join("motifs.json");
        let report = sync_jaspar_resource_impl(
            input_path.to_string_lossy().as_ref(),
            output_path.to_string_lossy().as_ref(),
        )
        .expect("sync jaspar");
        assert_eq!(report.resource, "jaspar");
        assert_eq!(report.item_count, 1);
        assert!(output_path.exists());
    }

    #[test]
    fn js_import_pool_wrapper_loads_member_into_state() {
        let td = tempdir().expect("tempdir");
        let pool_path = write_demo_pool_json(td.path());
        let state = ProjectState::default();
        let out = import_pool_impl(state, pool_path.to_string_lossy().as_ref(), "js_pool")
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
    fn js_sync_rebase_wrapper_matches_shared_shell_report() {
        let td = tempdir().expect("tempdir");
        let input_path = write_demo_rebase_withrefm(td.path());
        let output_path = td.path().join("rebase.json");

        let wrapper_report = sync_rebase_resource_impl(
            input_path.to_string_lossy().as_ref(),
            output_path.to_string_lossy().as_ref(),
            true,
        )
        .expect("wrapper sync rebase");

        let mut engine = GentleEngine::from_state(ProjectState::default());
        let shell_run = execute_shell_command(
            &mut engine,
            &ShellCommand::ResourcesSyncRebase {
                input: input_path.to_string_lossy().to_string(),
                output: Some(output_path.to_string_lossy().to_string()),
                commercial_only: true,
            },
        )
        .expect("shell sync rebase");
        let shell_report =
            sync_report_from_shell_output(shell_run.output, "resources sync-rebase shell")
                .expect("shell report");

        assert_eq!(
            serde_json::to_value(wrapper_report).expect("serialize wrapper report"),
            serde_json::to_value(shell_report).expect("serialize shell report")
        );
    }

    #[test]
    fn js_sync_jaspar_wrapper_matches_shared_shell_report() {
        let td = tempdir().expect("tempdir");
        let input_path = write_demo_jaspar_pfm(td.path());
        let output_path = td.path().join("motifs.json");

        let wrapper_report = sync_jaspar_resource_impl(
            input_path.to_string_lossy().as_ref(),
            output_path.to_string_lossy().as_ref(),
        )
        .expect("wrapper sync jaspar");

        let mut engine = GentleEngine::from_state(ProjectState::default());
        let shell_run = execute_shell_command(
            &mut engine,
            &ShellCommand::ResourcesSyncJaspar {
                input: input_path.to_string_lossy().to_string(),
                output: Some(output_path.to_string_lossy().to_string()),
            },
        )
        .expect("shell sync jaspar");
        let shell_report =
            sync_report_from_shell_output(shell_run.output, "resources sync-jaspar shell")
                .expect("shell report");

        assert_eq!(
            serde_json::to_value(wrapper_report).expect("serialize wrapper report"),
            serde_json::to_value(shell_report).expect("serialize shell report")
        );
    }

    #[test]
    fn js_import_pool_wrapper_matches_shared_shell_output_and_state() {
        let td = tempdir().expect("tempdir");
        let pool_path = write_demo_pool_json(td.path());

        let wrapper = import_pool_impl(
            ProjectState::default(),
            pool_path.to_string_lossy().as_ref(),
            "js_pool",
        )
        .expect("wrapper import pool");

        let mut shell_engine = GentleEngine::from_state(ProjectState::default());
        let shell_run = execute_shell_command(
            &mut shell_engine,
            &ShellCommand::ImportPool {
                input: pool_path.to_string_lossy().to_string(),
                prefix: "js_pool".to_string(),
            },
        )
        .expect("shell import pool");

        assert_eq!(wrapper.state_changed, shell_run.state_changed);
        assert_eq!(wrapper.output, shell_run.output);
        assert_eq!(
            wrapper
                .state
                .sequences
                .keys()
                .cloned()
                .collect::<std::collections::BTreeSet<_>>(),
            shell_engine
                .state()
                .sequences
                .keys()
                .cloned()
                .collect::<std::collections::BTreeSet<_>>(),
            "wrapper and shared-shell import_pool should yield same sequence ids"
        );
    }

    #[test]
    fn js_list_agent_systems_wrapper_reads_catalog() {
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
        let out = list_agent_systems_impl(catalog_path.to_string_lossy().as_ref())
            .expect("list agent systems");
        assert_eq!(out["schema"].as_str(), Some("gentle.agent_systems_list.v1"));
        assert_eq!(out["system_count"].as_u64(), Some(1));
        assert_eq!(out["systems"][0]["id"].as_str(), Some("builtin_echo"));
    }

    #[test]
    fn js_ask_agent_system_wrapper_uses_builtin_echo() {
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
        let out = ask_agent_system_impl(
            ProjectState::default(),
            "builtin_echo",
            "hello",
            catalog_path.to_string_lossy().as_ref(),
            false,
            false,
            &[],
            true,
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
    fn js_ask_agent_system_accepts_null_state() {
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
        let catalog_js = serde_json::to_string(catalog_path.to_string_lossy().as_ref())
            .expect("serialize catalog path");
        let mut js = JavaScriptInterface::default();
        js.run_checked(format!(
            r#"
                const out = ask_agent_system(null, "builtin_echo", "hello", {{
                    catalog_path: {catalog_js}
                }});
                if (out.output?.schema !== "gentle.agent_ask_result.v1") {{
                    throw new Error("unexpected schema");
                }}
            "#
        ))
        .expect("ask agent with null state");
    }

    #[test]
    fn js_dotplot_svg_wrapper_is_registered() {
        let mut js = JavaScriptInterface::default();
        js.run_checked(
            r#"
                if (typeof render_dotplot_svg !== "function") {
                    throw new Error("render_dotplot_svg wrapper is missing");
                }
            "#
            .to_string(),
        )
        .expect("render_dotplot_svg wrapper should be registered");
    }

    #[test]
    fn js_reference_and_helper_catalog_entry_wrappers_are_registered() {
        let mut js = JavaScriptInterface::default();
        js.run_checked(
            r#"
                if (typeof list_reference_catalog_entries !== "function") {
                    throw new Error("list_reference_catalog_entries wrapper is missing");
                }
                if (typeof list_helper_catalog_entries !== "function") {
                    throw new Error("list_helper_catalog_entries wrapper is missing");
                }
                if (typeof list_host_profile_catalog_entries !== "function") {
                    throw new Error("list_host_profile_catalog_entries wrapper is missing");
                }
                if (typeof list_ensembl_installable_genomes !== "function") {
                    throw new Error("list_ensembl_installable_genomes wrapper is missing");
                }
                if (typeof list_construct_reasoning_graphs !== "function") {
                    throw new Error("list_construct_reasoning_graphs wrapper is missing");
                }
                if (typeof show_construct_reasoning_graph !== "function") {
                    throw new Error("show_construct_reasoning_graph wrapper is missing");
                }
                if (typeof set_construct_reasoning_annotation_status !== "function") {
                    throw new Error("set_construct_reasoning_annotation_status wrapper is missing");
                }
                if (typeof write_back_construct_reasoning_annotation !== "function") {
                    throw new Error("write_back_construct_reasoning_annotation wrapper is missing");
                }
            "#
            .to_string(),
        )
        .expect("catalog entry wrappers should be registered");
    }

    #[test]
    fn js_helper_catalog_entry_wrapper_exposes_interpretation() {
        let td = tempdir().expect("tempdir");
        let catalog_path = td.path().join("helpers.json");
        fs::write(
            &catalog_path,
            r#"{
  "pGEX-demo": {
    "sequence_remote": "https://example.invalid/pgex.fa.gz",
    "annotations_remote": "https://example.invalid/pgex.gb.gz",
    "summary": "GST fusion expression vector",
    "helper_kind": "plasmid_vector",
    "host_system": "Escherichia coli",
    "search_terms": ["factor xa"],
    "semantics": {
      "schema": "gentle.helper_semantics.v1",
      "affordances": ["affinity_purification", "protease_tag_removal"],
      "components": [
        {"id": "gst", "kind": "fusion_tag", "label": "GST"},
        {"id": "mcs", "kind": "cloning_site", "label": "MCS"}
      ]
    }
  }
}"#,
        )
        .expect("write helpers catalog");

        let mut js = JavaScriptInterface::default();
        let catalog_js = serde_json::to_string(catalog_path.to_string_lossy().as_ref())
            .expect("serialize catalog path");
        js.run_checked(format!(
            r#"
                const rows = list_helper_catalog_entries({catalog_js}, "factor xa");
                if (rows.length !== 1) {{
                    throw new Error(`expected one helper row, got ${{rows.length}}`);
                }}
                const interpretation = rows[0].interpretation;
                if (!interpretation) {{
                    throw new Error("missing interpretation");
                }}
                if (interpretation.helper_kinds[0] !== "plasmid_vector") {{
                    throw new Error("missing helper kind");
                }}
                if (!interpretation.offered_functions.includes("fusion_tagging")) {{
                    throw new Error("missing derived function");
                }}
            "#
        ))
        .expect("helper catalog entries via js");
    }

    #[test]
    fn js_host_profile_catalog_entry_wrapper_lists_rows() {
        let td = tempdir().expect("tempdir");
        let catalog_path = td.path().join("host_profiles.json");
        fs::write(
            &catalog_path,
            r#"{
  "schema": "gentle.host_profile_catalog.v1",
  "profiles": [
    {
      "profile_id": "ecoli_dh5alpha",
      "species": "Escherichia coli",
      "strain": "DH5alpha",
      "aliases": ["DH5α"],
      "genotype_tags": ["deoR", "endA1"],
      "phenotype_tags": ["large_insert_friendly"],
      "notes": ["Common cloning host"],
      "source_notes": ["Synthetic regression fixture"]
    }
  ]
}"#,
        )
        .expect("write host profile catalog");

        let mut js = JavaScriptInterface::default();
        let catalog_js = serde_json::to_string(catalog_path.to_string_lossy().as_ref())
            .expect("serialize catalog path");
        js.run_checked(format!(
            r#"
                const rows = list_host_profile_catalog_entries({catalog_js}, "deoR");
                if (rows.length !== 1) {{
                    throw new Error(`expected one host row, got ${{rows.length}}`);
                }}
                if (rows[0].profile_id !== "ecoli_dh5alpha") {{
                    throw new Error("unexpected profile id");
                }}
            "#
        ))
        .expect("host profile catalog entries via js");
    }

    #[test]
    fn js_construct_reasoning_graph_wrappers_match_shared_shell_output() {
        let mut state = ProjectState::default();
        state.sequences.insert(
            "adapter_capture_js".to_string(),
            DNAsequence::from_sequence("GAATTCTCTAGAGCGGCCGCTTT").expect("sequence"),
        );
        let mut engine = GentleEngine::from_state(state.clone());
        let objective = engine
            .upsert_construct_objective(ConstructObjective {
                title: "JS adapter capture".to_string(),
                goal: "Expose adapter capture shell summary via JS".to_string(),
                adapter_restriction_capture_plans: vec![AdapterRestrictionCapturePlan {
                    capture_id: "mcs_capture".to_string(),
                    restriction_enzyme_name: "EcoRI".to_string(),
                    blunt_insert_required: true,
                    adapter_style: AdapterCaptureStyle::McsLike,
                    protection_mode: AdapterCaptureProtectionMode::InsertMethylation,
                    extra_retrieval_enzyme_names: vec!["XbaI".to_string(), "NotI".to_string()],
                    notes: vec![],
                }],
                ..ConstructObjective::default()
            })
            .expect("objective");
        let graph = engine
            .build_construct_reasoning_graph(
                "adapter_capture_js",
                Some(&objective.objective_id),
                Some("adapter_capture_js_graph"),
            )
            .expect("build graph");
        state = engine.state().clone();

        let wrapper_list =
            list_construct_reasoning_graphs_impl(state.clone(), "adapter_capture_js")
                .expect("js list wrapper");
        let wrapper_show = show_construct_reasoning_graph_impl(state.clone(), &graph.graph_id)
            .expect("js show wrapper");

        let mut shell_engine = GentleEngine::from_state(state);
        let shell_list = execute_shell_command(
            &mut shell_engine,
            &ShellCommand::ConstructReasoningListGraphs {
                seq_id: Some("adapter_capture_js".to_string()),
            },
        )
        .expect("shell list");
        let shell_show = execute_shell_command(
            &mut shell_engine,
            &ShellCommand::ConstructReasoningShowGraph {
                graph_id: graph.graph_id.clone(),
            },
        )
        .expect("shell show");

        assert_eq!(wrapper_list, shell_list.output);
        assert_eq!(wrapper_show, shell_show.output);
        assert!(
            wrapper_show["summary"]["fact_summaries"]
                .as_array()
                .map(|rows| rows.iter().any(|row| {
                    row["fact_type"].as_str() == Some("adapter_restriction_capture_context")
                }))
                .unwrap_or(false)
        );
    }

    #[test]
    fn js_construct_reasoning_annotation_status_wrapper_matches_shared_shell_output() {
        let mut dna = DNAsequence::from_sequence("ATGCGTATGCGT").expect("sequence");
        dna.features_mut().push(gb_io::seq::Feature {
            kind: "exon".into(),
            location: gb_io::seq::Location::simple_range(2, 10),
            qualifiers: vec![
                ("label".into(), Some("Confirmed exon".to_string())),
                ("evidence".into(), Some("supported by cDNA".to_string())),
            ],
        });
        dna.update_computed_features();

        let mut state = ProjectState::default();
        state
            .sequences
            .insert("construct_reasoning_js_status".to_string(), dna);
        let mut engine = GentleEngine::from_state(state.clone());
        let graph = engine
            .build_construct_reasoning_graph(
                "construct_reasoning_js_status",
                None,
                Some("construct_reasoning_js_status_graph"),
            )
            .expect("build graph");
        let annotation_id = graph
            .annotation_candidates
            .iter()
            .find(|candidate| candidate.role == crate::engine::ConstructRole::Exon)
            .map(|candidate| candidate.annotation_id.clone())
            .expect("annotation candidate id");
        state = engine.state().clone();

        let wrapper = set_construct_reasoning_annotation_status_impl(
            state.clone(),
            &graph.graph_id,
            &annotation_id,
            "accepted",
        )
        .expect("js status wrapper");

        let mut shell_engine = GentleEngine::from_state(state);
        let shell_run = execute_shell_command(
            &mut shell_engine,
            &ShellCommand::ConstructReasoningSetAnnotationStatus {
                graph_id: graph.graph_id.clone(),
                annotation_id,
                editable_status: crate::engine::EditableStatus::Accepted,
            },
        )
        .expect("shell status command");

        assert_eq!(wrapper.state_changed, shell_run.state_changed);
        assert_eq!(
            normalize_construct_reasoning_status_output(wrapper.output),
            normalize_construct_reasoning_status_output(shell_run.output)
        );
    }

    #[test]
    fn js_construct_reasoning_annotation_writeback_wrapper_matches_shared_shell_output() {
        let mut dna = DNAsequence::from_sequence(&"A".repeat(6001)).expect("sequence");
        dna.features_mut().push(gb_io::seq::Feature {
            kind: "mRNA".into(),
            location: gb_io::seq::Location::Complement(Box::new(
                gb_io::seq::Location::simple_range(2000, 2613),
            )),
            qualifiers: vec![
                ("gene".into(), Some("VKORC1".to_string())),
                ("transcript_id".into(), Some("ENSTVKORC1".to_string())),
                ("label".into(), Some("VKORC1-201".to_string())),
            ],
        });
        dna.update_computed_features();

        let mut state = ProjectState::default();
        state
            .sequences
            .insert("construct_reasoning_js_writeback".to_string(), dna);
        let mut engine = GentleEngine::from_state(state.clone());
        let graph = engine
            .build_construct_reasoning_graph(
                "construct_reasoning_js_writeback",
                None,
                Some("construct_reasoning_js_writeback_graph"),
            )
            .expect("build graph");
        let annotation_id = graph
            .annotation_candidates
            .iter()
            .find(|candidate| candidate.role == crate::engine::ConstructRole::Promoter)
            .map(|candidate| candidate.annotation_id.clone())
            .expect("annotation candidate id");
        engine
            .set_construct_reasoning_annotation_candidate_status(
                &graph.graph_id,
                &annotation_id,
                crate::engine::EditableStatus::Accepted,
            )
            .expect("accept annotation candidate");
        state = engine.state().clone();

        let wrapper = write_back_construct_reasoning_annotation_impl(
            state.clone(),
            &graph.graph_id,
            &annotation_id,
        )
        .expect("js writeback wrapper");

        let mut shell_engine = GentleEngine::from_state(state);
        let shell_run = execute_shell_command(
            &mut shell_engine,
            &ShellCommand::ConstructReasoningWriteAnnotation {
                graph_id: graph.graph_id.clone(),
                annotation_id,
            },
        )
        .expect("shell writeback command");

        assert_eq!(wrapper.state_changed, shell_run.state_changed);
        assert_eq!(
            normalize_construct_reasoning_status_output(wrapper.output),
            normalize_construct_reasoning_status_output(shell_run.output)
        );
    }

    #[test]
    fn js_apply_operation_export_run_bundle_matches_engine_decision_traces() {
        let td = tempdir().expect("tempdir");
        let wrapper_path = td.path().join("wrapper.run_bundle.json");
        let engine_path = td.path().join("engine.run_bundle.json");
        let operation = Operation::ExportProcessRunBundle {
            path: wrapper_path.to_string_lossy().to_string(),
            run_id: None,
        };
        let operation_json = serde_json::to_string(&operation).expect("serialize operation");
        let _wrapper = apply_operation_impl(decision_trace_fixture_state(), &operation_json)
            .expect("js apply operation export run bundle");

        let wrapper_bundle_text =
            fs::read_to_string(&wrapper_path).expect("read wrapper bundle output");
        let wrapper_bundle: crate::engine::ProcessRunBundleExport =
            serde_json::from_str(&wrapper_bundle_text).expect("parse wrapper bundle");

        let mut engine = GentleEngine::from_state(decision_trace_fixture_state());
        engine
            .apply(Operation::ExportProcessRunBundle {
                path: engine_path.to_string_lossy().to_string(),
                run_id: None,
            })
            .expect("engine export run bundle");
        let engine_bundle_text = fs::read_to_string(&engine_path).expect("read engine bundle");
        let engine_bundle: crate::engine::ProcessRunBundleExport =
            serde_json::from_str(&engine_bundle_text).expect("parse engine bundle");

        assert_eq!(
            serde_json::to_value(&wrapper_bundle.decision_traces)
                .expect("serialize wrapper traces"),
            serde_json::to_value(&engine_bundle.decision_traces).expect("serialize engine traces")
        );
        assert_eq!(wrapper_bundle.decision_traces.len(), 1);
        assert_eq!(
            wrapper_bundle.decision_traces[0].status,
            "preflight_failed".to_string()
        );
        assert_eq!(
            wrapper_bundle.decision_traces[0]
                .preflight_snapshot
                .as_ref()
                .expect("snapshot derived from history")
                .errors,
            vec!["missing sequence".to_string()]
        );
    }
}
