use gentle::lua_interface::LuaInterface;
use mlua::Value as LuaValue;
use rustyline::error::ReadlineError;
use rustyline::DefaultEditor;
use std::error::Error;

fn main() -> Result<(), Box<dyn Error>> {
    // Create a new Lua state
    let li = LuaInterface::new();

    // Register our custom Rust functions
    li.register_rust_functions()?;

    // Create line editor
    let mut rl = DefaultEditor::new()?;

    LuaInterface::help_main();

    loop {
        // Read line from user
        let readline = rl.readline("GENtle> ");

        match readline {
            Ok(line) => {
                if line.trim() == "exit" {
                    break;
                }

                rl.add_history_entry(line.as_str())?;

                // First try to evaluate as an expression
                let result = li.lua().load(&line).eval::<LuaValue>();

                // If that fails, try to execute as a statement
                let result = match result {
                    Ok(value) => Ok(value),
                    Err(_) => li.lua().load(&line).exec().map(|_| LuaValue::Nil),
                };

                match result {
                    Ok(value) => {
                        // Only print non-nil values
                        if !matches!(value, LuaValue::Nil) {
                            println!("{}", LuaInterface::format_lua_value(&value));
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
                // break;
            }
        }
    }

    Ok(())
}
