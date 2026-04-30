//! Lua shell binary entry point over shared engine operations.

use gentle::about;
use gentle::lua_interface::LuaInterface;
use mlua::Value as LuaValue;
use rustyline::DefaultEditor;
use rustyline::error::ReadlineError;
use std::{env, error::Error};

#[derive(Debug, Default)]
struct CliArgs {
    show_help: bool,
    show_version: bool,
    project_path: Option<String>,
}

fn print_help() {
    println!(
        "Usage:\n  \
gentle_lua [--help|-h] [--version|-V]\n  \
gentle_lua [--project PATH] [PATH]\n\n  \
If PATH is provided, the project is loaded into global variable 'project'."
    );
}

fn parse_cli_args(args: &[String]) -> Result<CliArgs, String> {
    let mut parsed = CliArgs::default();
    let mut idx = 0usize;
    while idx < args.len() {
        match args[idx].as_str() {
            "--help" | "-h" => {
                parsed.show_help = true;
                idx += 1;
            }
            "--version" | "-V" => {
                parsed.show_version = true;
                idx += 1;
            }
            "--project" => {
                if idx + 1 >= args.len() {
                    return Err("Missing PATH after --project".to_string());
                }
                parsed.project_path = Some(args[idx + 1].clone());
                idx += 2;
            }
            arg if arg.starts_with('-') => {
                return Err(format!("Unknown option '{arg}'"));
            }
            path => {
                if parsed.project_path.is_some() {
                    return Err(format!(
                        "Multiple project paths provided ('{}' and '{}')",
                        parsed.project_path.as_deref().unwrap_or_default(),
                        path
                    ));
                }
                parsed.project_path = Some(path.to_string());
                idx += 1;
            }
        }
    }
    Ok(parsed)
}

fn main() -> Result<(), Box<dyn Error>> {
    let args: Vec<String> = env::args().skip(1).collect();
    let cli = match parse_cli_args(&args) {
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
