//! Direct `gentle_cli allele-hash-screen` adapter.

use super::*;
use gentle::allele_hash_screen::{
    AlleleHashScreenRequest, AlleleReadPairInput, run_allele_hash_screen,
};

fn parse_allele_hash_screen_request(
    args: &[String],
    cmd_idx: usize,
) -> Result<AlleleHashScreenRequest, String> {
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
            "--from-rna-report" => {
                request.from_rna_report = Some(gentle_cli_args::required_value(
                    args,
                    idx,
                    "--from-rna-report",
                    "REPORT_ID",
                    "allele-hash-screen",
                )?);
                idx += 2;
            }
            "--salmon-unmapped-names" => {
                request.salmon_unmapped_names = Some(gentle_cli_args::required_value(
                    args,
                    idx,
                    "--salmon-unmapped-names",
                    "PATH",
                    "allele-hash-screen",
                )?);
                idx += 2;
            }
            "--salmon-mappings-sam" => {
                request.salmon_mappings_sam = Some(gentle_cli_args::required_value(
                    args,
                    idx,
                    "--salmon-mappings-sam",
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
            "--min-informative-reads" => {
                request.min_informative_reads = gentle_cli_args::parse_required_value(
                    args,
                    idx,
                    "--min-informative-reads",
                    "N",
                    "allele-hash-screen",
                    Some("allele-hash-screen"),
                )?;
                idx += 2;
            }
            "--balanced-band-lo" => {
                request.balanced_band_lo = gentle_cli_args::parse_required_value(
                    args,
                    idx,
                    "--balanced-band-lo",
                    "F",
                    "allele-hash-screen",
                    Some("allele-hash-screen"),
                )?;
                idx += 2;
            }
            "--balanced-band-hi" => {
                request.balanced_band_hi = gentle_cli_args::parse_required_value(
                    args,
                    idx,
                    "--balanced-band-hi",
                    "F",
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
    Ok(request)
}

pub(super) fn handle_allele_hash_screen(
    args: &[String],
    cmd_idx: usize,
    state_path: &str,
) -> Result<(), String> {
    let request = parse_allele_hash_screen_request(args, cmd_idx)?;
    let report = if request.from_rna_report.is_some() {
        let engine = GentleEngine::from_state(load_state(state_path)?);
        engine.run_allele_hash_screen_with_project_sources(request)?
    } else {
        run_allele_hash_screen(request)?
    };
    print_json(&report)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parses_report_and_salmon_source_flags() {
        let args = [
            "gentle_cli",
            "allele-hash-screen",
            "--gene",
            "FUS",
            "--transcript-fasta",
            "transcripts.fa",
            "--variant-table",
            "variants.tsv",
            "--from-rna-report",
            "fus_reads",
            "--read-file",
            "reads.fastq",
            "--salmon-unmapped-names",
            "unmapped_names.txt",
            "--salmon-mappings-sam",
            "mappings.sam",
            "--min-informative-reads",
            "12",
            "--balanced-band-lo",
            "0.45",
            "--balanced-band-hi",
            "0.55",
            "--out",
            "out",
        ]
        .map(str::to_string);
        let request =
            parse_allele_hash_screen_request(&args, 1).expect("parse allele source flags");

        assert_eq!(request.from_rna_report.as_deref(), Some("fus_reads"));
        assert_eq!(request.read_files, vec!["reads.fastq".to_string()]);
        assert_eq!(
            request.salmon_unmapped_names.as_deref(),
            Some("unmapped_names.txt")
        );
        assert_eq!(request.salmon_mappings_sam.as_deref(), Some("mappings.sam"));
        assert_eq!(request.min_informative_reads, 12);
        assert!((request.balanced_band_lo - 0.45).abs() < f64::EPSILON);
        assert!((request.balanced_band_hi - 0.55).abs() < f64::EPSILON);
    }
}
