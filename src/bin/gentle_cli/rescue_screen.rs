//! Direct `gentle_cli rescue-screen` adapter.

use super::*;
use gentle::target_rescue::{TargetRescueRequest, run_target_rescue_screen};

pub(super) fn handle_rescue_screen(args: &[String], cmd_idx: usize) -> Result<(), String> {
    let mut idx = cmd_idx + 1;
    let mut request = TargetRescueRequest::default();
    while idx < args.len() {
        match args[idx].as_str() {
            "--transcript-fasta" => {
                request
                    .transcript_fastas
                    .push(gentle_cli_args::required_value(
                        args,
                        idx,
                        "--transcript-fasta",
                        "PATH",
                        "rescue-screen",
                    )?);
                idx += 2;
            }
            "--genes" => {
                let raw = gentle_cli_args::required_value(
                    args,
                    idx,
                    "--genes",
                    "SYM[,SYM,...]",
                    "rescue-screen",
                )?;
                request.genes.extend(
                    raw.split(',')
                        .map(str::trim)
                        .filter(|gene| !gene.is_empty())
                        .map(str::to_string),
                );
                idx += 2;
            }
            "--gene" => {
                request.genes.push(gentle_cli_args::required_value(
                    args,
                    idx,
                    "--gene",
                    "SYM",
                    "rescue-screen",
                )?);
                idx += 2;
            }
            "--reads" => {
                request.reads.push(gentle_cli_args::required_value(
                    args,
                    idx,
                    "--reads",
                    "PATH",
                    "rescue-screen",
                )?);
                idx += 2;
            }
            "--read-id-allowlist" => {
                request.read_id_allowlist = Some(gentle_cli_args::required_value(
                    args,
                    idx,
                    "--read-id-allowlist",
                    "PATH",
                    "rescue-screen",
                )?);
                idx += 2;
            }
            "--salmon-unmapped-names" => {
                request.salmon_unmapped_names = Some(gentle_cli_args::required_value(
                    args,
                    idx,
                    "--salmon-unmapped-names",
                    "PATH",
                    "rescue-screen",
                )?);
                idx += 2;
            }
            "--salmon-mappings-sam" => {
                request.salmon_mappings_sam = Some(gentle_cli_args::required_value(
                    args,
                    idx,
                    "--salmon-mappings-sam",
                    "PATH",
                    "rescue-screen",
                )?);
                idx += 2;
            }
            "--kmer-len" => {
                request.kmer_len = gentle_cli_args::parse_required_value(
                    args,
                    idx,
                    "--kmer-len",
                    "N",
                    "rescue-screen",
                    Some("rescue-screen"),
                )?;
                idx += 2;
            }
            "--min-kmer-hits" => {
                request.min_kmer_hits = gentle_cli_args::parse_required_value(
                    args,
                    idx,
                    "--min-kmer-hits",
                    "N",
                    "rescue-screen",
                    Some("rescue-screen"),
                )?;
                idx += 2;
            }
            "--gene-symbol-tag" => {
                request
                    .gene_symbol_tags
                    .push(gentle_cli_args::required_value(
                        args,
                        idx,
                        "--gene-symbol-tag",
                        "KEY",
                        "rescue-screen",
                    )?);
                idx += 2;
            }
            "--output-prefix" => {
                request.output_prefix = gentle_cli_args::required_value(
                    args,
                    idx,
                    "--output-prefix",
                    "PREFIX",
                    "rescue-screen",
                )?;
                idx += 2;
            }
            other => return Err(format!("Unknown rescue-screen option '{other}'")),
        }
    }
    let summary = run_target_rescue_screen(request)?;
    print_json(&summary)
}
