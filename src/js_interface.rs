use crate::{
    app::GENtleApp,
    dna_sequence::DNAsequence,
    engine::{Engine, GentleEngine, Operation, ProjectState, Workflow},
    methylation_sites::MethylationMode,
    ENZYMES,
};
use deno_core::*;
use serde::Serialize;

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

#[op2]
#[serde]
fn load_dna(#[string] path: &str) -> Result<DNAsequence, deno_core::anyhow::Error> {
    let mut dna = GENtleApp::load_from_file(path)?;

    // Add default enzymes and stuff
    ENZYMES
        .restriction_enzymes()
        .clone_into(dna.restriction_enzymes_mut());
    dna.set_max_restriction_enzyme_sites(Some(2));
    dna.set_methylation_mode(MethylationMode::both());
    dna.update_computed_features();
    Ok(dna)
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
        const APPLY_OPERATION: OpDecl = apply_operation();
        const APPLY_WORKFLOW: OpDecl = apply_workflow();
        let ext = Extension {
            name: "my_ext",
            ops: std::borrow::Cow::Borrowed(&[
                LOAD_DNA,
                WRITE_GB,
                LOAD_PROJECT,
                SAVE_PROJECT,
                CAPABILITIES,
                APPLY_OPERATION,
                APPLY_WORKFLOW,
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
          	function apply_operation(state, op) {
          		const payload = (typeof op === "string") ? op : JSON.stringify(op);
          		return Deno.core.ops.apply_operation(state, payload);
          	}
          	function apply_workflow(state, workflow) {
          		const payload = (typeof workflow === "string") ? workflow : JSON.stringify(workflow);
          		return Deno.core.ops.apply_workflow(state, payload);
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
        "#
        .to_string();
        ret.run(init_code);
        println!("Interactive JavaScript Shell (type 'exit' to quit)");
        ret
    }

    pub fn run(&mut self, code: String) {
        match self.runtime.execute_script("<usage>", code) {
            Ok(_) => {}
            Err(e) => eprintln!("{}", e),
        }
    }
}

impl Default for JavaScriptInterface {
    fn default() -> Self {
        Self::new()
    }
}
