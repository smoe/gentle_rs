use crate::dna_sequence::DNAsequence;
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
        let mut dna = if let Ok(dna) = DNAsequence::from_genbank_file(path) {
            Self::first_dna_sequence(dna)
        } else if let Ok(dna) = DNAsequence::from_fasta_file(path) {
            Self::first_dna_sequence(dna)
        } else {
            return Err(Self::err(&format!(
                "Could not load DNA from file: {}",
                path
            )));
        }?;
        // Add default enzymes
        ENZYMES
            .restriction_enzymes()
            .clone_into(dna.restriction_enzymes_mut());
        dna.update_computed_features();
        Ok(dna)
    }

    fn first_dna_sequence(dna: Vec<DNAsequence>) -> Result<DNAsequence, Error> {
        Ok(dna
            .first()
            .ok_or_else(|| Self::err("No sequence in file"))?
            .to_owned())
    }

    fn err(s: &str) -> Error {
        Error::RuntimeError(s.to_string())
    }

    pub fn help_main() {
        println!("Interactive Lua Shell (type 'exit' to quit)");
        println!("Available Rust functions:");
        println!("  - load_dna(filename): Loads a DNA sequence from a file");
        println!("  - restriction_sites(seq): Returns restriction sites for a DNA sequence");
    }

    fn restriction_sites(mut seq: DNAsequence, lua: &Lua) -> LuaResult<Value> {
        let v = lua.to_value(seq.restriction_enzyme_sites())?;
        Ok(v)
    }

    pub fn register_rust_functions(&self) -> LuaResult<()> {
        self.lua.globals().set(
            "load_dna",
            self.lua
                .create_function(|_, filename: String| Self::load_dna(&filename))?,
        )?;

        self.lua.globals().set(
            "restriction_sites",
            self.lua
                .create_function(|lua, seq: DNAsequence| Self::restriction_sites(seq, lua))?,
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
