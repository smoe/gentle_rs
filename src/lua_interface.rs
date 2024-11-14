use crate::app::GENtleApp;
use crate::dna_sequence::DNAsequence;
use crate::methylation_sites::MethylationMode;
use crate::ENZYMES;
use mlua::prelude::*;
use mlua::{Error, Value};
use mlua::{IntoLuaMulti, Lua, MultiValue, Result as LuaResult};
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

        // self.lua.globals().set(
        //     "digest",
        //     self.lua
        //         .create_function(|_lua, (enzymes, seq): (String, DNAsequence)| {
        //             Self::restriction_enzyme_digest(seq, enzymes)
        //         })?,
        // )?;

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
