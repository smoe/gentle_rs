use crate::dna_sequence::DNAsequence;
use mlua::prelude::*;
use mlua::{Error, Value};
use mlua::{IntoLuaMulti, Lua, MultiValue, Result as LuaResult};
use serde::{Deserialize, Serialize};
use serde_json::json;

#[derive(Clone, Debug, Default, Serialize, Deserialize)]
pub struct LuaInterface {}

impl LuaInterface {
    pub fn load_dna(path: &str) -> Result<DNAsequence, Error> {
        if let Ok(dna) = DNAsequence::from_genbank_file(path) {
            Self::first_dna_sequence(dna)
        } else if let Ok(dna) = DNAsequence::from_fasta_file(path) {
            Self::first_dna_sequence(dna)
        } else {
            Err(Self::err(&format!(
                "Could not load DNA from file: {}",
                path
            )))
        }
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
    }

    fn restriction_sites(
        mut seq: DNAsequence,
    ) -> Vec<crate::restriction_enzyme::RestrictionEnzymeSite> {
        seq.update_computed_features();
        seq.restriction_enzyme_sites().to_owned()
    }

    pub fn register_rust_functions(lua: &Lua) -> LuaResult<()> {
        lua.globals().set(
            "load_dna",
            lua.create_function(|_, filename: String| Self::load_dna(&filename))?,
        )?;

        // lua.globals().set(
        //     "restriction_sites",
        //     lua.create_function(|_, filename: D| Self::load_dna(&filename))?,
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
