//! Direct reporter-recommender CLI family.

use super::*;

pub(super) fn handle_reporters_family(args: &[String], cmd_idx: usize) -> Result<(), String> {
    if args.len() <= cmd_idx + 1 {
        usage();
        return Err(
            "reporters requires a subcommand: list, recommend, or export-corpus".to_string(),
        );
    }
    match args[cmd_idx + 1].as_str() {
        "list" => handle_reporters_list(args, cmd_idx + 2),
        "recommend" => handle_reporters_recommend(args, cmd_idx + 2),
        "export-corpus" => handle_reporters_export_corpus(args, cmd_idx + 2),
        other => Err(format!("Unknown reporters subcommand '{other}'")),
    }
}

fn handle_reporters_list(args: &[String], mut idx: usize) -> Result<(), String> {
    let mut catalog_path: Option<String> = None;
    let mut filter: Option<String> = None;
    let mut limit: Option<usize> = None;
    let mut output: Option<String> = None;
    while idx < args.len() {
        match args[idx].as_str() {
            "--catalog" => {
                catalog_path = Some(required_value(args, idx, "--catalog")?);
                idx += 2;
            }
            "--filter" => {
                filter = Some(required_value(args, idx, "--filter")?);
                idx += 2;
            }
            "--limit" => {
                limit = Some(parse_usize_value(args, idx, "--limit")?);
                idx += 2;
            }
            "--output" => {
                output = Some(required_value(args, idx, "--output")?);
                idx += 2;
            }
            other => return Err(format!("Unknown reporters list option '{other}'")),
        }
    }
    let mut engine = GentleEngine::new();
    let result = engine
        .apply(Operation::ListReporterCatalog {
            catalog_path,
            filter,
            limit,
            path: output,
        })
        .map_err(|e| e.to_string())?;
    print_json(&result)
}

fn handle_reporters_recommend(args: &[String], mut idx: usize) -> Result<(), String> {
    let mut catalog_path: Option<String> = None;
    let mut limit: Option<usize> = None;
    let mut output: Option<String> = None;
    let mut constraints = ReporterConstraints::default();
    while idx < args.len() {
        match args[idx].as_str() {
            "--catalog" => {
                catalog_path = Some(required_value(args, idx, "--catalog")?);
                idx += 2;
            }
            "--limit" => {
                limit = Some(parse_usize_value(args, idx, "--limit")?);
                idx += 2;
            }
            "--output" => {
                output = Some(required_value(args, idx, "--output")?);
                idx += 2;
            }
            "--assay" => {
                constraints.intended_assay = Some(required_value(args, idx, "--assay")?);
                idx += 2;
            }
            "--chassis" => {
                constraints.chassis = Some(required_value(args, idx, "--chassis")?);
                idx += 2;
            }
            "--live" => {
                let raw = required_value(args, idx, "--live")?;
                constraints.live_assay = parse_bool_flag(&raw);
                if constraints.live_assay.is_none() {
                    return Err(format!("Invalid --live value '{raw}'; expected true/false"));
                }
                idx += 2;
            }
            "--color" => {
                constraints.desired_color = Some(required_value(args, idx, "--color")?);
                idx += 2;
            }
            "--class" => {
                constraints
                    .allowed_reporter_classes
                    .push(required_value(args, idx, "--class")?);
                idx += 2;
            }
            "--excitation-nm" => {
                constraints.available_excitation_nm.push(parse_u16_value(
                    args,
                    idx,
                    "--excitation-nm",
                )?);
                idx += 2;
            }
            "--emission-nm" => {
                constraints.available_emission_nm.push(parse_u16_value(
                    args,
                    idx,
                    "--emission-nm",
                )?);
                idx += 2;
            }
            "--fusion" => {
                constraints.fusion_mode = Some(required_value(args, idx, "--fusion")?);
                idx += 2;
            }
            "--max-length-bp" => {
                constraints.max_coding_length_bp =
                    Some(parse_usize_value(args, idx, "--max-length-bp")?);
                idx += 2;
            }
            "--forbid-motif" => {
                constraints
                    .forbidden_motifs
                    .push(required_value(args, idx, "--forbid-motif")?);
                idx += 2;
            }
            "--substrate-allowed" => {
                let raw = required_value(args, idx, "--substrate-allowed")?;
                constraints.substrate_allowed = parse_bool_flag(&raw);
                if constraints.substrate_allowed.is_none() {
                    return Err(format!(
                        "Invalid --substrate-allowed value '{raw}'; expected true/false"
                    ));
                }
                idx += 2;
            }
            other => return Err(format!("Unknown reporters recommend option '{other}'")),
        }
    }
    let mut engine = GentleEngine::new();
    let result = engine
        .apply(Operation::RecommendReporters {
            constraints,
            catalog_path,
            limit,
            path: output,
        })
        .map_err(|e| e.to_string())?;
    print_json(&result)
}

fn handle_reporters_export_corpus(args: &[String], mut idx: usize) -> Result<(), String> {
    if idx >= args.len() {
        usage();
        return Err("reporters export-corpus requires OUTPUT path".to_string());
    }
    let output = args[idx].clone();
    idx += 1;
    let mut catalog_path: Option<String> = None;
    let mut format = ReporterCorpusExportFormat::Json;
    while idx < args.len() {
        match args[idx].as_str() {
            "--catalog" => {
                catalog_path = Some(required_value(args, idx, "--catalog")?);
                idx += 2;
            }
            "--format" => {
                let raw = required_value(args, idx, "--format")?;
                format = match raw.as_str() {
                    "json" => ReporterCorpusExportFormat::Json,
                    "jsonl" => ReporterCorpusExportFormat::Jsonl,
                    _ => return Err(format!("Invalid reporters export format '{raw}'")),
                };
                idx += 2;
            }
            other => return Err(format!("Unknown reporters export-corpus option '{other}'")),
        }
    }
    let mut engine = GentleEngine::new();
    let result = engine
        .apply(Operation::ExportReporterCorpus {
            catalog_path,
            path: output,
            format,
        })
        .map_err(|e| e.to_string())?;
    print_json(&result)
}

fn required_value(args: &[String], idx: usize, flag: &str) -> Result<String, String> {
    args.get(idx + 1)
        .cloned()
        .ok_or_else(|| format!("Missing value after {flag}"))
}

fn parse_usize_value(args: &[String], idx: usize, flag: &str) -> Result<usize, String> {
    let raw = required_value(args, idx, flag)?;
    raw.parse::<usize>()
        .map_err(|e| format!("Invalid {flag} value '{raw}': {e}"))
}

fn parse_u16_value(args: &[String], idx: usize, flag: &str) -> Result<u16, String> {
    let raw = required_value(args, idx, flag)?;
    raw.parse::<u16>()
        .map_err(|e| format!("Invalid {flag} value '{raw}': {e}"))
}
