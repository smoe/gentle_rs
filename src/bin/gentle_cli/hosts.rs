//! Direct `hosts` CLI command handling kept out of `run()`.
//!
//! Public `hosts ...` invocations are forwarded through the shared shell parser
//! before direct dispatch. This module preserves the legacy direct handler so
//! the CLI adapter split can proceed one small command family at a time.

use super::*;

pub(super) fn handle_hosts_family(args: &[String], cmd_idx: usize) -> Result<(), String> {
    if args.len() <= cmd_idx + 1 {
        usage();
        return Err("hosts requires a subcommand".to_string());
    }
    match args[cmd_idx + 1].as_str() {
        "list" => {
            let mut catalog_path: Option<String> = None;
            let mut filter: Option<String> = None;
            let mut idx = cmd_idx + 2;
            while idx < args.len() {
                match args[idx].as_str() {
                    "--catalog" => {
                        if idx + 1 >= args.len() {
                            return Err("Missing PATH after --catalog for hosts list".to_string());
                        }
                        catalog_path = Some(args[idx + 1].clone());
                        idx += 2;
                    }
                    "--filter" => {
                        if idx + 1 >= args.len() {
                            return Err("Missing TEXT after --filter for hosts list".to_string());
                        }
                        filter = Some(args[idx + 1].clone());
                        idx += 2;
                    }
                    other => {
                        return Err(format!("Unknown option '{}' for hosts list", other));
                    }
                }
            }
            let entries = GentleEngine::list_host_profile_catalog_entries(
                explicit_catalog_arg(&catalog_path),
                filter.as_deref(),
            )
            .map_err(|e| e.to_string())?;
            let profile_ids = entries
                .iter()
                .map(|entry| entry.profile_id.clone())
                .collect::<Vec<_>>();
            print_json(&json!({
                "catalog_path": explicit_catalog_arg(&catalog_path)
                    .unwrap_or(DEFAULT_HOST_PROFILE_CATALOG_PATH),
                "filter": filter,
                "profile_count": profile_ids.len(),
                "profile_ids": profile_ids,
                "entries": entries,
            }))
        }
        other => Err(format!("Unknown hosts subcommand '{}'", other)),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn argv(values: &[&str]) -> Vec<String> {
        values.iter().map(|value| value.to_string()).collect()
    }

    #[test]
    fn hosts_missing_subcommand_reports_legacy_message() {
        let args = argv(&["gentle_cli", "hosts"]);
        let err = handle_hosts_family(&args, 1).expect_err("missing subcommand should fail");
        assert_eq!(err, "hosts requires a subcommand");
    }

    #[test]
    fn hosts_list_missing_catalog_path_reports_legacy_message() {
        let args = argv(&["gentle_cli", "hosts", "list", "--catalog"]);
        let err = handle_hosts_family(&args, 1).expect_err("missing catalog path should fail");
        assert_eq!(err, "Missing PATH after --catalog for hosts list");
    }

    #[test]
    fn hosts_unknown_subcommand_reports_legacy_message() {
        let args = argv(&["gentle_cli", "hosts", "show"]);
        let err = handle_hosts_family(&args, 1).expect_err("unknown subcommand should fail");
        assert_eq!(err, "Unknown hosts subcommand 'show'");
    }
}
