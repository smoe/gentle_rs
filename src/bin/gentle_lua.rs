//! Lua shell binary entry point over shared engine operations.

use gentle::lua_interface::LuaInterface;
use gentle::{
    about,
    cli_support::{SingleProjectCliOptions, parse_single_project_cli_args},
};
use mlua::Value as LuaValue;
use rustyline::DefaultEditor;
use rustyline::error::ReadlineError;
use std::{env, error::Error};

fn print_help() {
    println!(
        "Usage:\n  \
gentle_lua [--help|-h] [--version|-V]\n  \
gentle_lua [--project PATH] [PATH]\n\n  \
If PATH is provided, the project is loaded into global variable 'project'."
    );
}

fn main() -> Result<(), Box<dyn Error>> {
    let args: Vec<String> = env::args().skip(1).collect();
    let cli = match parse_single_project_cli_args(&args, SingleProjectCliOptions::standard()) {
        Ok(parsed) => parsed,
        Err(e) => {
            eprintln!("{e}");
            print_help();
            return Ok(());
        }
    };
    if cli.show_help {
        print_help();
        return Ok(());
    }
    if cli.show_version {
        println!("{}", about::version_cli_text());
        return Ok(());
    }

    // Create a new Lua state
    let li = LuaInterface::new();

    // Register our custom Rust functions
    li.register_rust_functions()?;

    if let Some(path) = cli.project_path {
        let path_json = serde_json::to_string(&path)?;
        li.lua()
            .load(format!("project = load_project({path_json})"))
            .exec()?;
        println!("Loaded project '{path}' into global variable 'project'");
    }

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
