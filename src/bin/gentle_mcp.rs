//! MCP stdio server binary for GENtle.

use gentle::{
    about,
    cli_support::parse_state_path_cli_args,
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
  - restriction_site_detail\n  \
  - op (requires confirm=true)\n  \
  - workflow (requires confirm=true)\n  \
  - help\n\n  \
MCP role in GENtle:\n  \
  - tool execution via tools/call\n  \
  - capability discovery/negotiation via tools/list, capabilities, and help\n\n  \
Default state path: {default_state}\n",
        default_state = DEFAULT_MCP_STATE_PATH
    );
}

fn run() -> Result<(), String> {
    let args = env::args().skip(1).collect::<Vec<_>>();
    let parsed = parse_state_path_cli_args(&args)?;
    if parsed.show_help {
        usage();
        return Ok(());
    }
    if parsed.show_version {
        println!("{}", about::version_cli_text());
        return Ok(());
    }
    let state_path = parsed
        .state_path
        .unwrap_or_else(|| DEFAULT_MCP_STATE_PATH.to_string());
    run_stdio_server(&state_path)
}

fn main() {
    if let Err(err) = run() {
        eprintln!("{err}");
        std::process::exit(1);
    }
}
