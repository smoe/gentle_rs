//! Direct `resources` CLI command handling kept out of `run()`.
//!
//! Most `resources ...` invocations are intentionally forwarded through the
//! shared shell parser before direct dispatch. This module preserves the legacy
//! direct handler as a small command-family island so future CLI cleanup can
//! move one family at a time without changing shell-forwarding behavior.

use super::*;

pub(super) const USAGE: &str = "\
  gentle_cli resources status
  gentle_cli resources sync-rebase INPUT.withrefm [OUTPUT.rebase.json] [--commercial-only]
  gentle_cli resources sync-jaspar INPUT.jaspar.txt [OUTPUT.motifs.json]

  gentle_cli resources sync-ucsc-rmsk INPUT.rmsk.txt_or_txt.gz [OUTPUT.rmsk.json] [--assembly DB] [--limit N]
  gentle_cli resources install-ucsc-rmsk [--assembly DB] [--input PATH_OR_URL] [--resource-output PATH] [--index-output PATH]
  gentle_cli resources prepare-ucsc-rmsk-index RESOURCE.rmsk.json [OUTPUT.interval-index.json]
  gentle_cli resources suggest-ucsc-rmsk-index [--assembly DB] [--output OUTPUT.json]

  gentle_cli resources sync-jaspar-remote-metadata [--motif TOKEN ...] [--motifs CSV] [--all] [--filter TOKEN] [--limit N] [--output OUTPUT.json]

  gentle_cli resources summarize-jaspar [--motif TOKEN ...] [--motifs CSV] [--all] [--random-length N] [--seed N] [--output OUTPUT.json]

  gentle_cli shell 'resources resolve-tf-query QUERY [QUERY ...] [--output OUTPUT.json]'

  gentle_cli resources benchmark-jaspar [--random-length N] [--seed N] [--output OUTPUT.json]

  gentle_cli resources list-jaspar [--filter TOKEN] [--limit N] [--fetch-remote] [--output OUTPUT.json]

  gentle_cli resources list-publication-datasets [--filter TEXT] [--catalog PATH] [--output OUTPUT.json]
  gentle_cli resources status-publication-dataset DATASET_ID [--catalog PATH] [--cache-dir PATH]
  gentle_cli resources prepare-publication-dataset DATASET_ID [--catalog PATH] [--cache-dir PATH] [--download-files] [--max-files N] [--category NAME|--categories CSV]

  gentle_cli resources inspect-jaspar MOTIF [--random-length N] [--seed N] [--fetch-remote] [--output OUTPUT.json]

";

pub(super) fn handle_resources_family(args: &[String], cmd_idx: usize) -> Result<(), String> {
    if args.len() <= cmd_idx + 1 {
        usage();
        return Err(
            "resources requires a subcommand: status, sync-rebase, sync-jaspar, sync-ucsc-rmsk, install-ucsc-rmsk, prepare-ucsc-rmsk-index, suggest-ucsc-rmsk-index, sync-jaspar-remote-metadata, summarize-jaspar, benchmark-jaspar, list-jaspar or inspect-jaspar"
                .to_string(),
        );
    }
    match args[cmd_idx + 1].as_str() {
        "status" => print_json(&resource_catalog_status()),
        "sync-rebase" => {
            if args.len() <= cmd_idx + 2 {
                usage();
                return Err("resources sync-rebase requires INPUT.withrefm path".to_string());
            }
            let input = &args[cmd_idx + 2];
            let mut output = DEFAULT_REBASE_RESOURCE_PATH.to_string();
            let mut commercial_only = false;
            for arg in args.iter().skip(cmd_idx + 3) {
                if arg == "--commercial-only" {
                    commercial_only = true;
                } else if !arg.starts_with("--") {
                    output = arg.clone();
                }
            }
            let text = read_text_input(input)?;
            let enzymes = parse_rebase_withrefm(&text, commercial_only);
            if enzymes.is_empty() {
                return Err(format!(
                    "No REBASE enzymes were parsed from '{input}'{}",
                    if commercial_only {
                        " (commercial-only filter active)"
                    } else {
                        ""
                    }
                ));
            }
            ensure_parent_dir(&output)?;
            let json = serde_json::to_string_pretty(&enzymes)
                .map_err(|e| format!("Could not serialize REBASE resource snapshot: {e}"))?;
            fs::write(&output, json)
                .map_err(|e| format!("Could not write REBASE output '{output}': {e}"))?;
            println!(
                "Synced {} REBASE enzymes to '{}'{}",
                enzymes.len(),
                output,
                if commercial_only {
                    " (commercial-only)"
                } else {
                    ""
                }
            );
            Ok(())
        }
        "sync-jaspar" => {
            if args.len() <= cmd_idx + 2 {
                usage();
                return Err("resources sync-jaspar requires INPUT.jaspar.txt path".to_string());
            }
            let input = &args[cmd_idx + 2];
            let output = args
                .get(cmd_idx + 3)
                .filter(|s| !s.starts_with("--"))
                .cloned()
                .unwrap_or_else(|| DEFAULT_JASPAR_RESOURCE_PATH.to_string());
            let text = read_text_input(input)?;
            let motifs = parse_jaspar_motifs(&text)?;
            if motifs.is_empty() {
                return Err(format!("No JASPAR motifs were parsed from '{input}'"));
            }
            let snapshot = JasparMotifSnapshot {
                schema: "gentle.tf_motifs.v2".to_string(),
                source: input.clone(),
                fetched_at_unix_ms: now_unix_ms(),
                motif_count: motifs.len(),
                motifs,
            };
            ensure_parent_dir(&output)?;
            let json = serde_json::to_string_pretty(&snapshot)
                .map_err(|e| format!("Could not serialize JASPAR resource snapshot: {e}"))?;
            fs::write(&output, json)
                .map_err(|e| format!("Could not write JASPAR output '{output}': {e}"))?;
            println!(
                "Synced {} JASPAR motifs to '{}'",
                snapshot.motif_count, output
            );
            Ok(())
        }
        "sync-jaspar-remote-metadata" => {
            let mut motifs: Vec<String> = vec![];
            let mut use_all = false;
            let mut filter: Option<String> = None;
            let mut limit: Option<usize> = None;
            let mut output: Option<String> = None;
            let mut idx = cmd_idx + 2;
            while idx < args.len() {
                match args[idx].as_str() {
                    "--motif" => {
                        if idx + 1 >= args.len() {
                            return Err(
                                "Missing TOKEN after --motif for resources sync-jaspar-remote-metadata"
                                    .to_string(),
                            );
                        }
                        motifs.push(args[idx + 1].clone());
                        idx += 2;
                    }
                    "--motifs" => {
                        if idx + 1 >= args.len() {
                            return Err(
                                "Missing CSV after --motifs for resources sync-jaspar-remote-metadata"
                                    .to_string(),
                            );
                        }
                        motifs.extend(
                            args[idx + 1]
                                .split(',')
                                .map(str::trim)
                                .filter(|value| !value.is_empty())
                                .map(str::to_string),
                        );
                        idx += 2;
                    }
                    "--all" => {
                        use_all = true;
                        idx += 1;
                    }
                    "--filter" => {
                        if idx + 1 >= args.len() {
                            return Err("Missing TOKEN after --filter".to_string());
                        }
                        filter = Some(args[idx + 1].clone());
                        idx += 2;
                    }
                    "--limit" => {
                        if idx + 1 >= args.len() {
                            return Err("Missing N after --limit".to_string());
                        }
                        limit = Some(args[idx + 1].parse().map_err(|e| {
                            format!(
                                "Invalid --limit '{}' for resources sync-jaspar-remote-metadata: {e}",
                                args[idx + 1]
                            )
                        })?);
                        idx += 2;
                    }
                    "--output" => {
                        if idx + 1 >= args.len() {
                            return Err("Missing PATH after --output".to_string());
                        }
                        output = Some(args[idx + 1].clone());
                        idx += 2;
                    }
                    other => {
                        return Err(format!(
                            "Unknown option '{}' for resources sync-jaspar-remote-metadata",
                            other
                        ));
                    }
                }
            }
            if use_all && !motifs.is_empty() {
                return Err(
                    "resources sync-jaspar-remote-metadata cannot combine --all with --motif/--motifs"
                        .to_string(),
                );
            }
            if use_all {
                motifs.clear();
            }
            let mut engine = GentleEngine::new();
            let op_result = engine
                .apply(Operation::SyncJasparRemoteMetadata {
                    motifs,
                    filter,
                    limit,
                    path: output.or_else(|| Some(DEFAULT_JASPAR_REMOTE_METADATA_PATH.to_string())),
                })
                .map_err(|e| e.to_string())?;
            print_json(&json!({ "result": op_result }))
        }
        "summarize-jaspar" => {
            let mut motifs: Vec<String> = vec![];
            let mut use_all = false;
            let mut random_sequence_length_bp =
                DEFAULT_JASPAR_PRESENTATION_RANDOM_SEQUENCE_LENGTH_BP;
            let mut random_seed = DEFAULT_JASPAR_PRESENTATION_RANDOM_SEED;
            let mut output: Option<String> = None;
            let mut idx = cmd_idx + 2;
            while idx < args.len() {
                match args[idx].as_str() {
                    "--motif" => {
                        if idx + 1 >= args.len() {
                            return Err(
                                "Missing TOKEN after --motif for resources summarize-jaspar"
                                    .to_string(),
                            );
                        }
                        motifs.push(args[idx + 1].clone());
                        idx += 2;
                    }
                    "--motifs" => {
                        if idx + 1 >= args.len() {
                            return Err(
                                "Missing CSV after --motifs for resources summarize-jaspar"
                                    .to_string(),
                            );
                        }
                        motifs.extend(
                            args[idx + 1]
                                .split(',')
                                .map(str::trim)
                                .filter(|value| !value.is_empty())
                                .map(str::to_string),
                        );
                        idx += 2;
                    }
                    "--all" => {
                        use_all = true;
                        idx += 1;
                    }
                    "--random-length" => {
                        if idx + 1 >= args.len() {
                            return Err("Missing N after --random-length".to_string());
                        }
                        random_sequence_length_bp = args[idx + 1].parse().map_err(|e| {
                            format!(
                                "Invalid --random-length '{}' for resources summarize-jaspar: {e}",
                                args[idx + 1]
                            )
                        })?;
                        idx += 2;
                    }
                    "--seed" => {
                        if idx + 1 >= args.len() {
                            return Err("Missing N after --seed".to_string());
                        }
                        random_seed = args[idx + 1].parse().map_err(|e| {
                            format!(
                                "Invalid --seed '{}' for resources summarize-jaspar: {e}",
                                args[idx + 1]
                            )
                        })?;
                        idx += 2;
                    }
                    "--output" => {
                        if idx + 1 >= args.len() {
                            return Err("Missing PATH after --output".to_string());
                        }
                        output = Some(args[idx + 1].clone());
                        idx += 2;
                    }
                    other => {
                        return Err(format!(
                            "Unknown option '{}' for resources summarize-jaspar",
                            other
                        ));
                    }
                }
            }
            if use_all && !motifs.is_empty() {
                return Err(
                    "resources summarize-jaspar cannot combine --all with --motif/--motifs"
                        .to_string(),
                );
            }
            if use_all {
                motifs.clear();
            }
            let mut engine = GentleEngine::new();
            let op_result = engine
                .apply(Operation::SummarizeJasparEntries {
                    motifs,
                    random_sequence_length_bp,
                    random_seed,
                    path: output,
                })
                .map_err(|e| e.to_string())?;
            print_json(&json!({ "result": op_result }))
        }
        "benchmark-jaspar" => {
            let mut random_sequence_length_bp =
                DEFAULT_JASPAR_PRESENTATION_RANDOM_SEQUENCE_LENGTH_BP;
            let mut random_seed = DEFAULT_JASPAR_PRESENTATION_RANDOM_SEED;
            let mut output: Option<String> = None;
            let mut idx = cmd_idx + 2;
            while idx < args.len() {
                match args[idx].as_str() {
                    "--random-length" => {
                        if idx + 1 >= args.len() {
                            return Err("Missing N after --random-length".to_string());
                        }
                        random_sequence_length_bp = args[idx + 1].parse().map_err(|e| {
                            format!(
                                "Invalid --random-length '{}' for resources benchmark-jaspar: {e}",
                                args[idx + 1]
                            )
                        })?;
                        idx += 2;
                    }
                    "--seed" => {
                        if idx + 1 >= args.len() {
                            return Err("Missing N after --seed".to_string());
                        }
                        random_seed = args[idx + 1].parse().map_err(|e| {
                            format!(
                                "Invalid --seed '{}' for resources benchmark-jaspar: {e}",
                                args[idx + 1]
                            )
                        })?;
                        idx += 2;
                    }
                    "--output" => {
                        if idx + 1 >= args.len() {
                            return Err("Missing PATH after --output".to_string());
                        }
                        output = Some(args[idx + 1].clone());
                        idx += 2;
                    }
                    other => {
                        return Err(format!(
                            "Unknown option '{}' for resources benchmark-jaspar",
                            other
                        ));
                    }
                }
            }
            let mut engine = GentleEngine::new();
            let op_result = engine
                .apply(Operation::BenchmarkJasparRegistry {
                    random_sequence_length_bp,
                    random_seed,
                    path: output.or_else(|| Some(DEFAULT_JASPAR_BENCHMARK_PATH.to_string())),
                })
                .map_err(|e| e.to_string())?;
            print_json(&json!({ "result": op_result }))
        }
        "list-jaspar" => {
            let mut filter: Option<String> = None;
            let mut limit: Option<usize> = None;
            let mut fetch_remote = false;
            let mut output: Option<String> = None;
            let mut idx = cmd_idx + 2;
            while idx < args.len() {
                match args[idx].as_str() {
                    "--filter" => {
                        if idx + 1 >= args.len() {
                            return Err("Missing TOKEN after --filter".to_string());
                        }
                        filter = Some(args[idx + 1].clone());
                        idx += 2;
                    }
                    "--limit" => {
                        if idx + 1 >= args.len() {
                            return Err("Missing N after --limit".to_string());
                        }
                        limit = Some(args[idx + 1].parse().map_err(|e| {
                            format!(
                                "Invalid --limit '{}' for resources list-jaspar: {e}",
                                args[idx + 1]
                            )
                        })?);
                        idx += 2;
                    }
                    "--fetch-remote" => {
                        fetch_remote = true;
                        idx += 1;
                    }
                    "--output" => {
                        if idx + 1 >= args.len() {
                            return Err("Missing PATH after --output".to_string());
                        }
                        output = Some(args[idx + 1].clone());
                        idx += 2;
                    }
                    other => {
                        return Err(format!(
                            "Unknown option '{}' for resources list-jaspar",
                            other
                        ));
                    }
                }
            }
            let mut engine = GentleEngine::new();
            let op_result = engine
                .apply(Operation::ListJasparCatalog {
                    filter,
                    limit,
                    include_remote_metadata: fetch_remote,
                    refresh_remote_metadata: fetch_remote,
                    path: output,
                })
                .map_err(|e| e.to_string())?;
            print_json(&json!({ "result": op_result }))
        }
        "inspect-jaspar" => {
            if args.len() <= cmd_idx + 2 {
                return Err(
                    "resources inspect-jaspar requires MOTIF [--random-length N] [--seed N] [--fetch-remote] [--output OUTPUT.json]"
                        .to_string(),
                );
            }
            let motif = args[cmd_idx + 2].clone();
            let mut random_sequence_length_bp =
                DEFAULT_JASPAR_PRESENTATION_RANDOM_SEQUENCE_LENGTH_BP;
            let mut random_seed = DEFAULT_JASPAR_PRESENTATION_RANDOM_SEED;
            let mut fetch_remote = false;
            let mut output: Option<String> = None;
            let mut idx = cmd_idx + 3;
            while idx < args.len() {
                match args[idx].as_str() {
                    "--random-length" => {
                        if idx + 1 >= args.len() {
                            return Err("Missing N after --random-length".to_string());
                        }
                        random_sequence_length_bp = args[idx + 1].parse().map_err(|e| {
                            format!(
                                "Invalid --random-length '{}' for resources inspect-jaspar: {e}",
                                args[idx + 1]
                            )
                        })?;
                        idx += 2;
                    }
                    "--seed" => {
                        if idx + 1 >= args.len() {
                            return Err("Missing N after --seed".to_string());
                        }
                        random_seed = args[idx + 1].parse().map_err(|e| {
                            format!(
                                "Invalid --seed '{}' for resources inspect-jaspar: {e}",
                                args[idx + 1]
                            )
                        })?;
                        idx += 2;
                    }
                    "--fetch-remote" => {
                        fetch_remote = true;
                        idx += 1;
                    }
                    "--output" => {
                        if idx + 1 >= args.len() {
                            return Err("Missing PATH after --output".to_string());
                        }
                        output = Some(args[idx + 1].clone());
                        idx += 2;
                    }
                    other => {
                        return Err(format!(
                            "Unknown option '{}' for resources inspect-jaspar",
                            other
                        ));
                    }
                }
            }
            let mut engine = GentleEngine::new();
            let op_result = engine
                .apply(Operation::InspectJasparEntry {
                    motif,
                    random_sequence_length_bp,
                    random_seed,
                    include_remote_metadata: fetch_remote,
                    refresh_remote_metadata: fetch_remote,
                    path: output,
                })
                .map_err(|e| e.to_string())?;
            print_json(&json!({ "result": op_result }))
        }
        other => Err(format!(
            "Unknown resources subcommand '{other}' (expected status, sync-rebase, sync-jaspar, sync-ucsc-rmsk, install-ucsc-rmsk, prepare-ucsc-rmsk-index, suggest-ucsc-rmsk-index, sync-jaspar-remote-metadata, summarize-jaspar, benchmark-jaspar, list-jaspar or inspect-jaspar)"
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
    fn usage_block_keeps_publication_dataset_routes() {
        assert!(USAGE.contains("resources list-publication-datasets"));
        assert!(USAGE.contains("resources prepare-publication-dataset"));
    }

    #[test]
    fn resources_missing_subcommand_reports_legacy_message() {
        let args = argv(&["gentle_cli", "resources"]);
        let err = handle_resources_family(&args, 1).expect_err("missing subcommand should fail");
        assert!(err.starts_with("resources requires a subcommand: status, sync-rebase"));
    }

    #[test]
    fn inspect_jaspar_missing_motif_reports_legacy_message() {
        let args = argv(&["gentle_cli", "resources", "inspect-jaspar"]);
        let err = handle_resources_family(&args, 1).expect_err("missing motif should fail");
        assert_eq!(
            err,
            "resources inspect-jaspar requires MOTIF [--random-length N] [--seed N] [--fetch-remote] [--output OUTPUT.json]"
        );
    }

    #[test]
    fn summarize_jaspar_rejects_all_with_explicit_motif_before_engine_dispatch() {
        let args = argv(&[
            "gentle_cli",
            "resources",
            "summarize-jaspar",
            "--all",
            "--motif",
            "MA0001.1",
        ]);
        let err = handle_resources_family(&args, 1).expect_err("conflicting motif flags fail");
        assert_eq!(
            err,
            "resources summarize-jaspar cannot combine --all with --motif/--motifs"
        );
    }
}
