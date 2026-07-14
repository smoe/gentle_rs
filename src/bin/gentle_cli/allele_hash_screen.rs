//! Direct `gentle_cli allele-hash-screen` adapter.

use super::*;
use gentle::allele_hash_screen::{
    AlleleHashScreenRequest, AlleleReadPairInput, run_allele_hash_screen,
};

pub(super) fn handle_allele_hash_screen(args: &[String], cmd_idx: usize) -> Result<(), String> {
    let mut idx = cmd_idx + 1;
    let mut request = AlleleHashScreenRequest::default();
    while idx < args.len() {
        match args[idx].as_str() {
            "--gene" => {
                request.gene = gentle_cli_args::required_value(
                    args,
                    idx,
                    "--gene",
                    "GENE",
                    "allele-hash-screen",
                )?;
                idx += 2;
            }
            "--transcript-fasta" => {
                request.transcript_fasta = gentle_cli_args::required_value(
                    args,
                    idx,
                    "--transcript-fasta",
                    "PATH",
                    "allele-hash-screen",
                )?;
                idx += 2;
            }
            "--variant-table" => {
                request.variant_table = Some(gentle_cli_args::required_value(
                    args,
                    idx,
                    "--variant-table",
                    "PATH",
                    "allele-hash-screen",
                )?);
                idx += 2;
            }
            "--vcf" => {
                request.vcf = Some(gentle_cli_args::required_value(
                    args,
                    idx,
                    "--vcf",
                    "PATH",
                    "allele-hash-screen",
                )?);
                idx += 2;
            }
            "--transcript-map" => {
                request.transcript_map = Some(gentle_cli_args::required_value(
                    args,
                    idx,
                    "--transcript-map",
                    "PATH",
                    "allele-hash-screen",
                )?);
                idx += 2;
            }
            "--vcf-sample" => {
                request.vcf_sample = Some(gentle_cli_args::required_value(
                    args,
                    idx,
                    "--vcf-sample",
                    "SAMPLE",
                    "allele-hash-screen",
                )?);
                idx += 2;
            }
            "--read-file" => {
                request.read_files.push(gentle_cli_args::required_value(
                    args,
                    idx,
                    "--read-file",
                    "PATH",
                    "allele-hash-screen",
                )?);
                idx += 2;
            }
            "--read-pair" => {
                let raw = gentle_cli_args::required_value(
                    args,
                    idx,
                    "--read-pair",
                    "R1,R2",
                    "allele-hash-screen",
                )?;
                let (read1, read2) = raw
                    .split_once(',')
                    .ok_or_else(|| format!("Invalid --read-pair value '{raw}': expected R1,R2"))?;
                if read1.trim().is_empty() || read2.trim().is_empty() {
                    return Err(format!(
                        "Invalid --read-pair value '{raw}': both paths are required"
                    ));
                }
                request.read_pairs.push(AlleleReadPairInput {
                    read1: read1.trim().to_string(),
                    read2: read2.trim().to_string(),
                });
                idx += 2;
            }
            "--read-id-allowlist" => {
                request.read_id_allowlist = Some(gentle_cli_args::required_value(
                    args,
                    idx,
                    "--read-id-allowlist",
                    "PATH",
                    "allele-hash-screen",
                )?);
                idx += 2;
            }
            "--out" => {
                request.out_dir = gentle_cli_args::required_value(
                    args,
                    idx,
                    "--out",
                    "OUT_DIR",
                    "allele-hash-screen",
                )?;
                idx += 2;
            }
            "--kmer-len" => {
                request.kmer_len = gentle_cli_args::parse_required_value(
                    args,
                    idx,
                    "--kmer-len",
                    "N",
                    "allele-hash-screen",
                    Some("allele-hash-screen"),
                )?;
                idx += 2;
            }
            "--min-unique-kmer-hits" => {
                request.min_unique_kmer_hits = gentle_cli_args::parse_required_value(
                    args,
                    idx,
                    "--min-unique-kmer-hits",
                    "N",
                    "allele-hash-screen",
                    Some("allele-hash-screen"),
                )?;
                idx += 2;
            }
            "--max-inline-read-calls" => {
                request.max_inline_read_calls = gentle_cli_args::parse_required_value(
                    args,
                    idx,
                    "--max-inline-read-calls",
                    "N",
                    "allele-hash-screen",
                    Some("allele-hash-screen"),
                )?;
                idx += 2;
            }
            other => return Err(format!("Unknown allele-hash-screen option '{other}'")),
        }
    }
    let report = run_allele_hash_screen(request)?;
    print_json(&report)
}
