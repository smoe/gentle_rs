//! MCP stdio server binary for GENtle.

use gentle::{
    about,
    mcp_server::{DEFAULT_MCP_STATE_PATH, run_stdio_server},
};
use std::env;

fn usage() {
    println!(
        "Usage:\n  \
gentle_mcp [--state PATH|--project PATH] [--help|-h] [--version|-V]\n\n  \
Starts a MCP stdio server with tools:\n  \
  - capabilities\n  \
  - state_summary\n  \
  - op (requires confirm=true)\n  \
  - workflow (requires confirm=true)\n  \
  - help\n\n  \
Default state path: {default_state}\n",
        default_state = DEFAULT_MCP_STATE_PATH
    );
}

fn parse_state_path(args: &[String]) -> Result<Option<String>, String> {
    let mut state_path: Option<String> = None;
    let mut idx = 1usize;
    while idx < args.len() {
        match args[idx].as_str() {
            "--state" | "--project" => {
                if idx + 1 >= args.len() {
                    return Err(format!("Missing PATH after {}", args[idx]));
                }
                state_path = Some(args[idx + 1].clone());
                idx += 2;
            }
            "--help" | "-h" => return Ok(None),
            "--version" | "-V" => return Ok(None),
            other => {
                return Err(format!("Unknown argument '{other}'. Use --help for usage."));
            }
        }
    }
    Ok(state_path)
}

fn run() -> Result<(), String> {
    let args = env::args().collect::<Vec<_>>();
    if args.iter().any(|a| a == "--help" || a == "-h") {
        usage();
        return Ok(());
    }
    if args.iter().any(|a| a == "--version" || a == "-V") {
        println!("{}", about::version_cli_text());
        return Ok(());
    }
    let state_path = parse_state_path(&args)?.unwrap_or_else(|| DEFAULT_MCP_STATE_PATH.to_string());
    run_stdio_server(&state_path)
}

fn main() {
    if let Err(err) = run() {
        eprintln!("{err}");
        std::process::exit(1);
    }
}
