//! JavaScript shell binary entry point over shared engine operations.

use gentle::js_interface::JavaScriptInterface;
use gentle::{
    about,
    cli_support::{SingleProjectCliOptions, parse_single_project_cli_args},
};
use rustyline::DefaultEditor;
use rustyline::error::ReadlineError;
use std::{env, error::Error};

fn print_help() {
    println!(
        "Usage:\n  \
gentle_js [--help|-h] [--version|-V]\n  \
gentle_js [--project PATH] [PATH]\n\n  \
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

    // Create the JavaScript runtime
    let mut js = JavaScriptInterface::default();
    if let Some(path) = cli.project_path {
        let path_json = serde_json::to_string(&path)?;
        js.run_checked(format!("globalThis.project = load_project({path_json});"))
            .map_err(|e| format!("Could not load project '{path}': {e}"))?;
        println!("Loaded project '{path}' into global variable 'project'");
    }

    // Create line editor
    let mut rl = DefaultEditor::new()?;

    loop {
        // Read line from user
        let readline = rl.readline("GENtle> ");

        match readline {
            Ok(line) => {
                if line.trim() == "exit" {
                    break;
                }

                rl.add_history_entry(line.as_str())?;

                js.run(line);

                // let result = runtime.execute_script("<usage>", line);
                // match result {
                //     Ok(_state) => {
                //         // println!("!{:?}", &value);
                //     }
                //     Err(e) => println!("Error: {}", e),
                // }
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
