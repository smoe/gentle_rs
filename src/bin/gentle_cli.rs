use mlua::{Lua, Result as LuaResult, Value};
use rustyline::error::ReadlineError;
use rustyline::DefaultEditor;
use std::error::Error;

// Example custom Rust function that will be callable from Lua
fn rust_multiply(a: i64, b: i64) -> i64 {
    a * b
}

// Example custom Rust function that works with strings
fn rust_reverse_string(s: String) -> String {
    s.chars().rev().collect()
}

fn register_rust_functions(lua: &Lua) -> LuaResult<()> {
    // Register the multiply function
    lua.globals().set(
        "rust_multiply",
        lua.create_function(|_, (a, b): (i64, i64)| Ok(rust_multiply(a, b)))?,
    )?;

    // Register the string reverse function
    lua.globals().set(
        "rust_reverse",
        lua.create_function(|_, s: String| Ok(rust_reverse_string(s)))?,
    )?;

    Ok(())
}

fn format_lua_value(value: &Value) -> String {
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

fn main() -> Result<(), Box<dyn Error>> {
    // Create a new Lua state
    let lua = Lua::new();

    // Register our custom Rust functions
    register_rust_functions(&lua)?;

    // Create line editor
    let mut rl = DefaultEditor::new()?;

    println!("Interactive Lua Shell (type 'exit' to quit)");
    println!("Available Rust functions:");
    println!("  - rust_multiply(a, b): Multiplies two numbers");
    println!("  - rust_reverse(s): Reverses a string");

    loop {
        // Read line from user
        let readline = rl.readline("lua> ");

        match readline {
            Ok(line) => {
                if line.trim() == "exit" {
                    break;
                }

                rl.add_history_entry(line.as_str())?;

                // First try to evaluate as an expression
                let result = lua.load(&line).eval::<Value>();

                // If that fails, try to execute as a statement
                let result = match result {
                    Ok(value) => Ok(value),
                    Err(_) => lua.load(&line).exec().map(|_| Value::Nil),
                };

                match result {
                    Ok(value) => {
                        // Only print non-nil values
                        if !matches!(value, Value::Nil) {
                            println!("{}", format_lua_value(&value));
                        }
                    }
                    Err(e) => println!("Error: {}", e),
                }
            }
            Err(ReadlineError::Interrupted) => {
                println!("CTRL-C");
                break;
            }
            Err(ReadlineError::Eof) => {
                println!("CTRL-D");
                break;
            }
            Err(err) => {
                println!("Error: {}", err);
                break;
            }
        }
    }

    Ok(())
}

// use gentle::dna_sequence::DNAsequence;
// use mlua::prelude::*;
// use mlua::Error;
// use mlua::{UserData, UserDataFields, UserDataMethods};
// use serde::{Deserialize, Serialize};
// use std::io::Write;

// #[derive(Clone, Debug, Default, Serialize, Deserialize)]
// struct LuaInterface {
//     value: i32,
// }

// impl LuaInterface {
//     fn load_dna(path: &str) -> Result<DNAsequence, Error> {
//         if let Ok(dna) = DNAsequence::from_genbank_file(path) {
//             Self::first_dna_sequence(dna)
//         } else if let Ok(dna) = DNAsequence::from_fasta_file(path) {
//             Self::first_dna_sequence(dna)
//         } else {
//             Err(Self::err(&format!(
//                 "Could not load DNA from file: {}",
//                 path
//             )))
//         }
//     }

//     fn first_dna_sequence(dna: Vec<DNAsequence>) -> Result<DNAsequence, Error> {
//         Ok(dna
//             .first()
//             .ok_or_else(|| Self::err("No sequence in file"))?
//             .to_owned())
//     }

//     fn err(s: &str) -> Error {
//         Error::RuntimeError(s.to_string())
//     }
// }

// impl UserData for LuaInterface {
//     fn add_fields<F: UserDataFields<Self>>(fields: &mut F) {
//         fields.add_field_method_get("val", |_, this| Ok(this.value));
//     }

//     fn add_methods<M: UserDataMethods<Self>>(methods: &mut M) {
//         methods.add_function("new", |_, value: i32| Ok(Self { value }));
//         methods.add_function("hello", |_, value: String| Ok(format!("Hello, {}!", value)));
//     }
// }

// fn main() -> LuaResult<()> {
//     let lua = Lua::new();

//     lua.globals()
//         .set("gentle", lua.create_proxy::<LuaInterface>()?)?;

//     println!("GENtle interactive Lua interface. Type `exit` or press Ctrl-C to exit.");
//     loop {
//         print!(">> ");
//         std::io::stdout().lock().flush().unwrap();
//         let mut input = String::new();
//         std::io::stdin().read_line(&mut input).unwrap();
//         let input = input.trim();
//         if input == "exit" {
//             break;
//         }

//         // let input = format!("function _() return {} end", input);
//         // let chunk: mlua::Function = lua.load(input).eval()?;
//         let chunk = lua.load(input);

//         match chunk.eval::<String>() {
//             Ok(result) => println!("=> {result:?}"),
//             Err(e) => eprintln!("Error: {}", e),
//         }
//     }

//     Ok(())
// }
