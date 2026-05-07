//! Direct `services` CLI command handling kept out of `run()`.
//!
//! Public `services ...` routes are normally forwarded through the shared shell
//! parser before direct dispatch. This module preserves the small legacy direct
//! handler so the CLI adapter can keep shrinking one command family at a time.

use super::*;

pub(super) const USAGE: &str = "\
  gentle_cli services status
  gentle_cli services providers list
  gentle_cli services providers doctor [--catalog PATH] [--output PATH]
  gentle_cli services project-preflight REQUEST_JSON_OR_@FILE
  gentle_cli services project-quote REQUEST_JSON_OR_@FILE

  gentle_cli services handoff [--scope NAME] [--output PATH]
  gentle_cli services guide --channel telegram [--section overview|readiness|gene-context|tfbs|inline-dna|cloning|isoforms|follow-up] [--gene SYMBOL]

";

pub(super) fn handle_services_family(args: &[String], cmd_idx: usize) -> Result<(), String> {
    if args.len() <= cmd_idx + 1 {
        usage();
        return Err("services requires a subcommand: status".to_string());
    }
    match args[cmd_idx + 1].as_str() {
        "status" => print_json(&service_readiness_status()?),
        other => Err(format!(
            "Unknown services subcommand '{other}' (expected status)"
        )),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn argv(values: &[&str]) -> Vec<String> {
        values.iter().map(|value| value.to_string()).collect()
    }

    #[test]
    fn usage_block_keeps_forwarded_service_routes() {
        assert!(USAGE.contains("services providers list"));
        assert!(USAGE.contains("services providers doctor"));
        assert!(USAGE.contains("services project-preflight"));
        assert!(USAGE.contains("services handoff"));
        assert!(USAGE.contains("services guide --channel telegram"));
    }

    #[test]
    fn services_missing_subcommand_reports_legacy_message() {
        let args = argv(&["gentle_cli", "services"]);
        let err = handle_services_family(&args, 1).expect_err("missing subcommand should fail");
        assert_eq!(err, "services requires a subcommand: status");
    }

    #[test]
    fn services_unknown_subcommand_reports_legacy_message() {
        let args = argv(&["gentle_cli", "services", "providers"]);
        let err = handle_services_family(&args, 1).expect_err("unknown subcommand should fail");
        assert_eq!(
            err,
            "Unknown services subcommand 'providers' (expected status)"
        );
    }
}
