use crate::app::GENtleApp;
use crate::dna_sequence::DNAsequence;
use crate::engine::{Engine, GentleEngine, Operation, ProjectState, Workflow};
use crate::methylation_sites::MethylationMode;
use crate::ENZYMES;
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
        ENZYMES
            .restriction_enzymes()
            .clone_into(dna.restriction_enzymes_mut());
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
        println!(
            "  - apply_operation(project, op): Applies an operation (Lua table or JSON string)"
        );
        println!("  - apply_workflow(project, wf): Applies a workflow (Lua table or JSON string)");
        println!("A sequence has the following properties:\n- seq.restriction_enzymes\n- seq.restriction_enzyme_sites\n- seq.open_reading_frames\n- seq.methylation_sites");
    }

    fn write_gb(seq: DNAsequence, filename: String) -> LuaResult<bool> {
        // lua.to_value(seq.restriction_enzyme_sites())
        seq.write_genbank_file(&filename)
            .map_err(|e| Self::err(&format!("{}", e)))?;
        Ok(true)
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
