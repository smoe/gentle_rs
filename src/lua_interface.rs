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

    pub fn register_rust_functions(lua: &Lua) -> LuaResult<()> {
        // Register the multiply function
        // lua.globals().set(
        //     "rust_multiply",
        //     lua.create_function(|_, (a, b): (i64, i64)| Ok(rust_multiply(a, b)))?,
        // )?;

        // Register the string reverse function
        // lua.globals().set(
        //     "rust_reverse",
        //     lua.create_function(|_, s: String| Ok(rust_reverse_string(s)))?,
        // )?;

        lua.globals().set(
            "load_dna",
            lua.create_function(|_, filename: String| Self::load_dna(&filename))?,
        )?;

        Ok(())
    }

    pub fn format_lua_value(value: &Value) -> String {
        match value {
            Value::Nil => "nil".to_string(),
            Value::Boolean(b) => b.to_string(),
            Value::Integer(i) => i.to_string(),
            Value::Number(n) => n.to_string(),
            Value::String(s) => {
                if let Ok(str) = s.to_str() {
                    str.to_string()
                } else {
                    "<binary string>".to_string()
                }
            }
            Value::Table(_) => "<table>".to_string(),
            Value::Function(_) => "<function>".to_string(),
            Value::Thread(_) => "<thread>".to_string(),
            Value::UserData(_) => "<userdata>".to_string(),
            Value::Error(e) => format!("<error: {}>", e),
            _ => "<unknown>".to_string(),
        }
    }
}

// impl UserData for LuaInterface {
//     fn add_fields<F: UserDataFields<Self>>(fields: &mut F) {
//         fields.add_field_method_get("val", |_, this| Ok(this.value));
//     }

//     fn add_methods<M: UserDataMethods<Self>>(methods: &mut M) {
//         methods.add_function("new", |_, value: i32| Ok(Self { value }));
//         methods.add_function("hello", |_, value: String| Ok(format!("Hello, {}!", value)));
//     }
// }

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
