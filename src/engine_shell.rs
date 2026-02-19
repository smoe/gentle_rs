use crate::{
    dna_ladder::LadderMolecule,
    dna_sequence::DNAsequence,
    engine::{Engine, GentleEngine, Operation, ProjectState, RenderSvgMode, Workflow},
    genomes::{GenomeGeneRecord, DEFAULT_GENOME_CATALOG_PATH, DEFAULT_HELPER_GENOME_CATALOG_PATH},
    resource_sync,
};
#[cfg(all(target_os = "macos", feature = "screenshot-capture"))]
use objc2_app_kit::NSApplication;
#[cfg(all(target_os = "macos", feature = "screenshot-capture"))]
use objc2_foundation::MainThreadMarker;
use regex::{Regex, RegexBuilder};
use serde::{Deserialize, Serialize};
use serde_json::{json, Value};
#[cfg(all(target_os = "macos", feature = "screenshot-capture"))]
use std::path::Path;
#[cfg(all(target_os = "macos", feature = "screenshot-capture"))]
use std::process::Command;
use std::{collections::BTreeSet, fs};

#[derive(Debug, Clone)]
pub enum ShellCommand {
    Help,
    Capabilities,
    StateSummary,
    LoadProject {
        path: String,
    },
    SaveProject {
        path: String,
    },
    ScreenshotWindow {
        output: String,
    },
    RenderSvg {
        seq_id: String,
        mode: RenderSvgMode,
        output: String,
    },
    RenderRnaSvg {
        seq_id: String,
        output: String,
    },
    RnaInfo {
        seq_id: String,
    },
    RenderLineageSvg {
        output: String,
    },
    RenderPoolGelSvg {
        inputs: Vec<String>,
        output: String,
        ladders: Option<Vec<String>>,
    },
    LaddersList {
        molecule: LadderMolecule,
        name_filter: Option<String>,
    },
    LaddersExport {
        molecule: LadderMolecule,
        output: String,
        name_filter: Option<String>,
    },
    ExportPool {
        inputs: Vec<String>,
        output: String,
        human_id: Option<String>,
    },
    ImportPool {
        input: String,
        prefix: String,
    },
    ResourcesSyncRebase {
        input: String,
        output: Option<String>,
        commercial_only: bool,
    },
    ResourcesSyncJaspar {
        input: String,
        output: Option<String>,
    },
    ReferenceList {
        helper_mode: bool,
        catalog_path: Option<String>,
    },
    ReferenceValidateCatalog {
        helper_mode: bool,
        catalog_path: Option<String>,
    },
    ReferenceStatus {
        helper_mode: bool,
        genome_id: String,
        catalog_path: Option<String>,
        cache_dir: Option<String>,
    },
    ReferenceGenes {
        helper_mode: bool,
        genome_id: String,
        catalog_path: Option<String>,
        cache_dir: Option<String>,
        filter: String,
        biotypes: Vec<String>,
        limit: usize,
        offset: usize,
    },
    ReferencePrepare {
        helper_mode: bool,
        genome_id: String,
        catalog_path: Option<String>,
        cache_dir: Option<String>,
    },
    ReferenceExtractRegion {
        helper_mode: bool,
        genome_id: String,
        chromosome: String,
        start_1based: usize,
        end_1based: usize,
        output_id: Option<String>,
        catalog_path: Option<String>,
        cache_dir: Option<String>,
    },
    ReferenceExtractGene {
        helper_mode: bool,
        genome_id: String,
        gene_query: String,
        occurrence: Option<usize>,
        output_id: Option<String>,
        catalog_path: Option<String>,
        cache_dir: Option<String>,
    },
    ReferenceBlast {
        helper_mode: bool,
        genome_id: String,
        query_sequence: String,
        max_hits: usize,
        task: Option<String>,
        catalog_path: Option<String>,
        cache_dir: Option<String>,
    },
    TracksImportBed {
        seq_id: String,
        path: String,
        track_name: Option<String>,
        min_score: Option<f64>,
        max_score: Option<f64>,
        clear_existing: bool,
    },
    TracksImportBigWig {
        seq_id: String,
        path: String,
        track_name: Option<String>,
        min_score: Option<f64>,
        max_score: Option<f64>,
        clear_existing: bool,
    },
    Op {
        payload: String,
    },
    Workflow {
        payload: String,
    },
}

#[derive(Debug, Clone)]
pub struct ShellRunResult {
    pub state_changed: bool,
    pub output: Value,
}

#[derive(Debug, Clone, Copy, Default)]
pub struct ShellExecutionOptions {
    pub allow_screenshots: bool,
}

impl ShellExecutionOptions {
    pub fn from_env() -> Self {
        let raw = std::env::var("GENTLE_ALLOW_SCREENSHOTS").unwrap_or_default();
        let allow_screenshots = cfg!(feature = "screenshot-capture")
            && matches!(
                raw.trim().to_ascii_lowercase().as_str(),
                "1" | "true" | "yes" | "on"
            );
        Self { allow_screenshots }
    }
}

#[derive(Debug, Clone, Serialize, Deserialize)]
struct ScreenshotReport {
    schema: String,
    path: String,
    window_title: String,
    captured_at_unix_ms: u128,
    pixel_width: u32,
    pixel_height: u32,
    backend: String,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
struct PoolEnd {
    end_type: String,
    forward_5: String,
    forward_3: String,
    reverse_5: String,
    reverse_3: String,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
struct PoolMember {
    seq_id: String,
    #[serde(default)]
    human_id: String,
    name: Option<String>,
    sequence: String,
    length_bp: usize,
    topology: String,
    ends: PoolEnd,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
struct PoolExport {
    schema: String,
    pool_id: String,
    #[serde(default)]
    human_id: String,
    member_count: usize,
    members: Vec<PoolMember>,
}

impl ShellCommand {
    pub fn preview(&self) -> String {
        match self {
            Self::Help => "show shell command help".to_string(),
            Self::Capabilities => "inspect engine capabilities".to_string(),
            Self::StateSummary => "show sequence/container state summary".to_string(),
            Self::LoadProject { path } => format!("load project state from '{path}'"),
            Self::SaveProject { path } => format!("save current project state to '{path}'"),
            Self::ScreenshotWindow { output } => {
                format!("capture active GENtle window screenshot to '{output}'")
            }
            Self::RenderSvg {
                seq_id,
                mode,
                output,
            } => format!("render {mode:?} SVG for '{seq_id}' to '{output}'"),
            Self::RenderRnaSvg { seq_id, output } => {
                format!("render RNA structure SVG for '{seq_id}' to '{output}'")
            }
            Self::RnaInfo { seq_id } => {
                format!("inspect rnapkin textual RNA report for '{seq_id}'")
            }
            Self::RenderLineageSvg { output } => format!("render lineage SVG to '{output}'"),
            Self::RenderPoolGelSvg {
                inputs,
                output,
                ladders,
            } => {
                let ladders = ladders
                    .as_ref()
                    .map(|v| v.join(","))
                    .unwrap_or_else(|| "auto".to_string());
                format!(
                    "render pool gel SVG from {} input(s) to '{output}' with ladders {ladders}",
                    inputs.len()
                )
            }
            Self::LaddersList {
                molecule,
                name_filter,
            } => {
                let filter = name_filter
                    .as_deref()
                    .map(str::trim)
                    .filter(|v| !v.is_empty())
                    .unwrap_or("-");
                format!(
                    "inspect {} ladders (filter={filter})",
                    molecule.display_name()
                )
            }
            Self::LaddersExport {
                molecule,
                output,
                name_filter,
            } => {
                let filter = name_filter
                    .as_deref()
                    .map(str::trim)
                    .filter(|v| !v.is_empty())
                    .unwrap_or("-");
                format!(
                    "export {} ladders to '{output}' (filter={filter})",
                    molecule.display_name()
                )
            }
            Self::ExportPool {
                inputs,
                output,
                human_id,
            } => {
                let human = human_id.clone().unwrap_or_else(|| "-".to_string());
                format!(
                    "export pool with {} input(s) to '{output}' (human_id={human})",
                    inputs.len()
                )
            }
            Self::ImportPool { input, prefix } => {
                format!("import pool from '{input}' with prefix '{prefix}'")
            }
            Self::ResourcesSyncRebase {
                input,
                output,
                commercial_only,
            } => {
                let output = output
                    .clone()
                    .unwrap_or_else(|| resource_sync::DEFAULT_REBASE_RESOURCE_PATH.to_string());
                format!(
                    "sync REBASE from '{input}' to '{output}'{}",
                    if *commercial_only {
                        " (commercial-only)"
                    } else {
                        ""
                    }
                )
            }
            Self::ResourcesSyncJaspar { input, output } => {
                let output = output
                    .clone()
                    .unwrap_or_else(|| resource_sync::DEFAULT_JASPAR_RESOURCE_PATH.to_string());
                format!("sync JASPAR from '{input}' to '{output}'")
            }
            Self::ReferenceList {
                helper_mode,
                catalog_path,
            } => {
                let label = if *helper_mode { "helpers" } else { "genomes" };
                let catalog = catalog_path
                    .clone()
                    .unwrap_or_else(|| default_catalog_path(*helper_mode).to_string());
                format!("list {label} from catalog '{catalog}'")
            }
            Self::ReferenceValidateCatalog {
                helper_mode,
                catalog_path,
            } => {
                let label = if *helper_mode { "helpers" } else { "genomes" };
                let catalog = catalog_path
                    .clone()
                    .unwrap_or_else(|| default_catalog_path(*helper_mode).to_string());
                format!("validate {label} catalog '{catalog}'")
            }
            Self::ReferenceStatus {
                helper_mode,
                genome_id,
                catalog_path,
                cache_dir,
            } => {
                let label = if *helper_mode { "helper" } else { "genome" };
                let catalog = catalog_path
                    .clone()
                    .unwrap_or_else(|| default_catalog_path(*helper_mode).to_string());
                let cache = cache_dir.clone().unwrap_or_else(|| "-".to_string());
                format!("check {label} '{genome_id}' status (catalog='{catalog}', cache='{cache}')")
            }
            Self::ReferenceGenes {
                helper_mode,
                genome_id,
                catalog_path,
                cache_dir,
                filter,
                biotypes,
                limit,
                offset,
            } => {
                let label = if *helper_mode { "helper" } else { "genome" };
                let catalog = catalog_path
                    .clone()
                    .unwrap_or_else(|| default_catalog_path(*helper_mode).to_string());
                let cache = cache_dir.clone().unwrap_or_else(|| "-".to_string());
                let biotypes = if biotypes.is_empty() {
                    "-".to_string()
                } else {
                    biotypes.join(",")
                };
                format!(
                    "list {label} genes for '{genome_id}' (catalog='{catalog}', cache='{cache}', filter='{filter}', biotypes='{biotypes}', limit={limit}, offset={offset})"
                )
            }
            Self::ReferencePrepare {
                helper_mode,
                genome_id,
                catalog_path,
                cache_dir,
            } => {
                let label = if *helper_mode { "helper" } else { "genome" };
                let catalog = catalog_path
                    .clone()
                    .unwrap_or_else(|| default_catalog_path(*helper_mode).to_string());
                let cache = cache_dir.clone().unwrap_or_else(|| "-".to_string());
                format!("prepare {label} '{genome_id}' (catalog='{catalog}', cache='{cache}')")
            }
            Self::ReferenceExtractRegion {
                helper_mode,
                genome_id,
                chromosome,
                start_1based,
                end_1based,
                output_id,
                catalog_path,
                cache_dir,
            } => {
                let label = if *helper_mode { "helper" } else { "genome" };
                let catalog = catalog_path
                    .clone()
                    .unwrap_or_else(|| default_catalog_path(*helper_mode).to_string());
                let cache = cache_dir.clone().unwrap_or_else(|| "-".to_string());
                let output = output_id.clone().unwrap_or_else(|| "-".to_string());
                format!(
                    "extract {label} region {genome_id}:{chromosome}:{start_1based}-{end_1based} (output='{output}', catalog='{catalog}', cache='{cache}')"
                )
            }
            Self::ReferenceExtractGene {
                helper_mode,
                genome_id,
                gene_query,
                occurrence,
                output_id,
                catalog_path,
                cache_dir,
            } => {
                let label = if *helper_mode { "helper" } else { "genome" };
                let catalog = catalog_path
                    .clone()
                    .unwrap_or_else(|| default_catalog_path(*helper_mode).to_string());
                let cache = cache_dir.clone().unwrap_or_else(|| "-".to_string());
                let occ = occurrence
                    .map(|v| v.to_string())
                    .unwrap_or_else(|| "-".to_string());
                let output = output_id.clone().unwrap_or_else(|| "-".to_string());
                format!(
                    "extract {label} gene '{gene_query}' from '{genome_id}' (occurrence={occ}, output='{output}', catalog='{catalog}', cache='{cache}')"
                )
            }
            Self::ReferenceBlast {
                helper_mode,
                genome_id,
                query_sequence,
                max_hits,
                task,
                catalog_path,
                cache_dir,
            } => {
                let label = if *helper_mode { "helper" } else { "genome" };
                let catalog = catalog_path
                    .clone()
                    .unwrap_or_else(|| default_catalog_path(*helper_mode).to_string());
                let cache = cache_dir.clone().unwrap_or_else(|| "-".to_string());
                let task = task.clone().unwrap_or_else(|| "blastn-short".to_string());
                format!(
                    "blast query (len={}) against {label} '{genome_id}' (max_hits={max_hits}, task='{task}', catalog='{catalog}', cache='{cache}')",
                    query_sequence.len()
                )
            }
            Self::TracksImportBed {
                seq_id,
                path,
                track_name,
                min_score,
                max_score,
                clear_existing,
            } => format!(
                "import BED track for '{seq_id}' from '{}' (track_name='{}', min_score={}, max_score={}, clear_existing={})",
                path,
                track_name
                    .clone()
                    .unwrap_or_else(|| "-".to_string()),
                min_score
                    .map(|v| v.to_string())
                    .unwrap_or_else(|| "-".to_string()),
                max_score
                    .map(|v| v.to_string())
                    .unwrap_or_else(|| "-".to_string()),
                clear_existing
            ),
            Self::TracksImportBigWig {
                seq_id,
                path,
                track_name,
                min_score,
                max_score,
                clear_existing,
            } => format!(
                "import BigWig track for '{seq_id}' from '{}' (track_name='{}', min_score={}, max_score={}, clear_existing={})",
                path,
                track_name
                    .clone()
                    .unwrap_or_else(|| "-".to_string()),
                min_score
                    .map(|v| v.to_string())
                    .unwrap_or_else(|| "-".to_string()),
                max_score
                    .map(|v| v.to_string())
                    .unwrap_or_else(|| "-".to_string()),
                clear_existing
            ),
            Self::Op { .. } => "apply one engine operation from JSON".to_string(),
            Self::Workflow { .. } => "apply engine workflow from JSON".to_string(),
        }
    }

    pub fn is_state_mutating(&self) -> bool {
        matches!(
            self,
            Self::LoadProject { .. }
                | Self::ImportPool { .. }
                | Self::ReferencePrepare { .. }
                | Self::ReferenceExtractRegion { .. }
                | Self::ReferenceExtractGene { .. }
                | Self::TracksImportBed { .. }
                | Self::TracksImportBigWig { .. }
                | Self::Op { .. }
                | Self::Workflow { .. }
        )
    }
}

pub fn shell_help_text() -> String {
    let screenshot_line = if cfg!(feature = "screenshot-capture") {
        "screenshot-window OUTPUT.png\n"
    } else {
        "screenshot-window OUTPUT.png (disabled in this build; enable feature 'screenshot-capture')\n"
    };
    format!(
        "GENtle Shell commands:\n\
help\n\
capabilities\n\
state-summary\n\
load-project PATH\n\
save-project PATH\n\
{screenshot_line}\
render-svg SEQ_ID linear|circular OUTPUT.svg\n\
render-rna-svg SEQ_ID OUTPUT.svg\n\
rna-info SEQ_ID\n\
render-lineage-svg OUTPUT.svg\n\
render-pool-gel-svg IDS OUTPUT.svg [--ladders NAME[,NAME]]\n\
ladders list [--molecule dna|rna] [--filter TEXT]\n\
ladders export OUTPUT.json [--molecule dna|rna] [--filter TEXT]\n\
export-pool IDS OUTPUT.pool.gentle.json [HUMAN_ID]\n\
import-pool INPUT.pool.gentle.json [PREFIX]\n\
resources sync-rebase INPUT.withrefm_or_URL [OUTPUT.rebase.json] [--commercial-only]\n\
resources sync-jaspar INPUT.jaspar_or_URL [OUTPUT.motifs.json]\n\
genomes list [--catalog PATH]\n\
genomes validate-catalog [--catalog PATH]\n\
genomes status GENOME_ID [--catalog PATH] [--cache-dir PATH]\n\
genomes genes GENOME_ID [--catalog PATH] [--cache-dir PATH] [--filter REGEX] [--biotype NAME] [--limit N] [--offset N]\n\
genomes prepare GENOME_ID [--catalog PATH] [--cache-dir PATH]\n\
genomes blast GENOME_ID QUERY_SEQUENCE [--max-hits N] [--task blastn-short|blastn] [--catalog PATH] [--cache-dir PATH]\n\
genomes extract-region GENOME_ID CHR START END [--output-id ID] [--catalog PATH] [--cache-dir PATH]\n\
genomes extract-gene GENOME_ID QUERY [--occurrence N] [--output-id ID] [--catalog PATH] [--cache-dir PATH]\n\
helpers list [--catalog PATH]\n\
helpers validate-catalog [--catalog PATH]\n\
helpers status HELPER_ID [--catalog PATH] [--cache-dir PATH]\n\
helpers genes HELPER_ID [--catalog PATH] [--cache-dir PATH] [--filter REGEX] [--biotype NAME] [--limit N] [--offset N]\n\
helpers prepare HELPER_ID [--catalog PATH] [--cache-dir PATH]\n\
helpers blast HELPER_ID QUERY_SEQUENCE [--max-hits N] [--task blastn-short|blastn] [--catalog PATH] [--cache-dir PATH]\n\
helpers extract-region HELPER_ID CHR START END [--output-id ID] [--catalog PATH] [--cache-dir PATH]\n\
helpers extract-gene HELPER_ID QUERY [--occurrence N] [--output-id ID] [--catalog PATH] [--cache-dir PATH]\n\
tracks import-bed SEQ_ID PATH [--name NAME] [--min-score N] [--max-score N] [--clear-existing]\n\
tracks import-bigwig SEQ_ID PATH [--name NAME] [--min-score N] [--max-score N] [--clear-existing]\n\
op <operation-json-or-@file>\n\
workflow <workflow-json-or-@file>\n\
IDS is comma-separated sequence IDs"
    )
}

fn split_ids(input: &str) -> Vec<String> {
    input
        .split(',')
        .map(|s| s.trim())
        .filter(|s| !s.is_empty())
        .map(|s| s.to_string())
        .collect()
}

fn unique_id(existing: &std::collections::HashMap<String, DNAsequence>, base: &str) -> String {
    if !existing.contains_key(base) {
        return base.to_string();
    }
    let mut i = 2usize;
    loop {
        let candidate = format!("{base}_{i}");
        if !existing.contains_key(&candidate) {
            return candidate;
        }
        i += 1;
    }
}

fn apply_member_overhang(member: &PoolMember, dna: &mut DNAsequence) -> Result<(), String> {
    let mut value =
        serde_json::to_value(&*dna).map_err(|e| format!("Could not serialize sequence: {e}"))?;
    let obj = value
        .as_object_mut()
        .ok_or_else(|| "Internal serialization shape error".to_string())?;
    obj.insert(
        "overhang".to_string(),
        json!({
            "forward_3": member.ends.forward_3.as_bytes(),
            "forward_5": member.ends.forward_5.as_bytes(),
            "reverse_3": member.ends.reverse_3.as_bytes(),
            "reverse_5": member.ends.reverse_5.as_bytes(),
        }),
    );
    let patched: DNAsequence =
        serde_json::from_value(value).map_err(|e| format!("Could not restore overhang: {e}"))?;
    *dna = patched;
    Ok(())
}

fn default_catalog_path(helper_mode: bool) -> &'static str {
    if helper_mode {
        DEFAULT_HELPER_GENOME_CATALOG_PATH
    } else {
        DEFAULT_GENOME_CATALOG_PATH
    }
}

fn operation_catalog_path(catalog_path: &Option<String>, helper_mode: bool) -> Option<String> {
    catalog_path
        .as_deref()
        .map(str::trim)
        .filter(|v| !v.is_empty())
        .map(|v| v.to_string())
        .or_else(|| helper_mode.then(|| default_catalog_path(helper_mode).to_string()))
}

fn resolved_catalog_path<'a>(
    catalog_path: &'a Option<String>,
    helper_mode: bool,
) -> Option<&'a str> {
    if let Some(value) = catalog_path
        .as_deref()
        .map(str::trim)
        .filter(|v| !v.is_empty())
    {
        Some(value)
    } else if helper_mode {
        Some(default_catalog_path(helper_mode))
    } else {
        None
    }
}

fn effective_catalog_path(catalog_path: &Option<String>, helper_mode: bool) -> String {
    catalog_path
        .as_deref()
        .map(str::trim)
        .filter(|v| !v.is_empty())
        .map(|v| v.to_string())
        .unwrap_or_else(|| default_catalog_path(helper_mode).to_string())
}

fn compile_gene_filter_regex(filter: &str) -> Result<Option<Regex>, String> {
    let pattern = filter.trim();
    if pattern.is_empty() {
        return Ok(None);
    }
    RegexBuilder::new(pattern)
        .case_insensitive(true)
        .build()
        .map(Some)
        .map_err(|e| format!("Invalid --filter regex '{}': {}", pattern, e))
}

fn genome_gene_matches_regex_filter(gene: &GenomeGeneRecord, regex: &Regex) -> bool {
    gene.gene_name
        .as_ref()
        .map(|name| regex.is_match(name))
        .unwrap_or(false)
        || gene
            .gene_id
            .as_ref()
            .map(|id| regex.is_match(id))
            .unwrap_or(false)
        || regex.is_match(&gene.chromosome)
}

fn collect_biotypes(genes: &[GenomeGeneRecord]) -> Vec<String> {
    let mut biotypes: BTreeSet<String> = BTreeSet::new();
    for gene in genes {
        let Some(biotype) = gene.biotype.as_ref() else {
            continue;
        };
        let trimmed = biotype.trim();
        if !trimmed.is_empty() {
            biotypes.insert(trimmed.to_string());
        }
    }
    biotypes.into_iter().collect()
}

fn genome_gene_matches_filter(
    gene: &GenomeGeneRecord,
    filter_regex: Option<&Regex>,
    allowed_biotypes_lower: &[String],
) -> bool {
    let regex_ok = filter_regex
        .map(|re| genome_gene_matches_regex_filter(gene, re))
        .unwrap_or(true);
    if !regex_ok {
        return false;
    }
    if allowed_biotypes_lower.is_empty() {
        return true;
    }
    gene.biotype
        .as_ref()
        .map(|b| b.trim().to_ascii_lowercase())
        .map(|probe| allowed_biotypes_lower.iter().any(|b| b == &probe))
        .unwrap_or(false)
}

fn parse_option_path(
    tokens: &[String],
    idx: &mut usize,
    option_name: &str,
    context: &str,
) -> Result<String, String> {
    if *idx + 1 >= tokens.len() {
        return Err(format!("Missing value after {option_name} for {context}"));
    }
    let value = tokens[*idx + 1].clone();
    *idx += 2;
    Ok(value)
}

fn parse_mode(mode: &str) -> Result<RenderSvgMode, String> {
    match mode {
        "linear" => Ok(RenderSvgMode::Linear),
        "circular" => Ok(RenderSvgMode::Circular),
        other => Err(format!(
            "Unknown render mode '{other}', expected 'linear' or 'circular'"
        )),
    }
}

fn parse_ladder_molecule(value: &str) -> Result<LadderMolecule, String> {
    LadderMolecule::parse(value)
        .ok_or_else(|| format!("Unknown ladder molecule '{value}', expected 'dna' or 'rna'"))
}

#[cfg(all(target_os = "macos", feature = "screenshot-capture"))]
fn now_unix_ms() -> u128 {
    std::time::SystemTime::now()
        .duration_since(std::time::UNIX_EPOCH)
        .map(|d| d.as_millis())
        .unwrap_or(0)
}

#[cfg(all(target_os = "macos", feature = "screenshot-capture"))]
fn ensure_parent_dir(path: &str) -> Result<(), String> {
    let parent = Path::new(path)
        .parent()
        .map(|p| p.to_path_buf())
        .unwrap_or_else(|| Path::new(".").to_path_buf());
    fs::create_dir_all(parent)
        .map_err(|e| format!("Could not create output directory for '{path}': {e}"))
}

#[cfg(all(target_os = "macos", feature = "screenshot-capture"))]
fn active_window_info_from_appkit() -> Result<(u64, String), String> {
    let Some(mtm) = MainThreadMarker::new() else {
        return Err("Screenshot capture requires the macOS main thread".to_string());
    };
    let app = NSApplication::sharedApplication(mtm);
    let Some(window) = app.keyWindow().or_else(|| unsafe { app.mainWindow() }) else {
        return Err(
            "No active GENtle window in this process; run from GUI shell with the target window focused"
                .to_string(),
        );
    };
    let raw_window_id = unsafe { window.windowNumber() };
    if raw_window_id <= 0 {
        return Err(format!(
            "Active GENtle window reported invalid window number {raw_window_id}"
        ));
    }
    let window_id = u64::try_from(raw_window_id)
        .map_err(|_| format!("Could not convert window number {raw_window_id} to u64"))?;
    let title = window.title().to_string();
    let window_title = if title.trim().is_empty() {
        "<untitled>".to_string()
    } else {
        title.trim().to_string()
    };
    Ok((window_id, window_title))
}

#[cfg(all(target_os = "macos", feature = "screenshot-capture"))]
fn capture_active_window_screenshot(output: &str) -> Result<ScreenshotReport, String> {
    ensure_parent_dir(output)?;
    let (window_id, window_title) = active_window_info_from_appkit()?;

    let status = Command::new("screencapture")
        .arg("-x")
        .arg("-l")
        .arg(window_id.to_string())
        .arg(output)
        .status()
        .map_err(|e| format!("Could not run macOS screencapture command: {e}"))?;
    if !status.success() {
        return Err(format!(
            "macOS screencapture failed with status {:?}",
            status.code()
        ));
    }
    let (pixel_width, pixel_height) = image::image_dimensions(output)
        .map_err(|e| format!("Screenshot written but dimensions could not be read: {e}"))?;

    Ok(ScreenshotReport {
        schema: "gentle.screenshot.v1".to_string(),
        path: output.to_string(),
        window_title,
        captured_at_unix_ms: now_unix_ms(),
        pixel_width,
        pixel_height,
        backend: "macos.screencapture".to_string(),
    })
}

#[cfg(any(not(target_os = "macos"), not(feature = "screenshot-capture")))]
fn capture_active_window_screenshot(_output: &str) -> Result<ScreenshotReport, String> {
    if !cfg!(feature = "screenshot-capture") {
        Err(
            "screenshot-window is unavailable in this build; enable feature 'screenshot-capture'"
                .to_string(),
        )
    } else {
        Err("screenshot-window is currently supported only on macOS".to_string())
    }
}

fn parse_json_payload(raw: &str) -> Result<String, String> {
    if let Some(path) = raw.strip_prefix('@') {
        fs::read_to_string(path).map_err(|e| format!("Could not read JSON file '{path}': {e}"))
    } else {
        Ok(raw.to_string())
    }
}

fn token_error(command: &str) -> String {
    format!("Invalid '{command}' usage. Try: help")
}

fn parse_reference_command(tokens: &[String], helper_mode: bool) -> Result<ShellCommand, String> {
    let label = if helper_mode { "helpers" } else { "genomes" };
    if tokens.len() < 2 {
        return Err(format!("{label} requires a subcommand"));
    }
    match tokens[1].as_str() {
        "list" => {
            let mut catalog_path: Option<String> = None;
            let mut idx = 2usize;
            while idx < tokens.len() {
                match tokens[idx].as_str() {
                    "--catalog" => {
                        catalog_path = Some(parse_option_path(tokens, &mut idx, "--catalog", label)?)
                    }
                    other => {
                        return Err(format!("Unknown option '{other}' for {label} list"));
                    }
                }
            }
            Ok(ShellCommand::ReferenceList {
                helper_mode,
                catalog_path,
            })
        }
        "validate-catalog" => {
            let mut catalog_path: Option<String> = None;
            let mut idx = 2usize;
            while idx < tokens.len() {
                match tokens[idx].as_str() {
                    "--catalog" => {
                        catalog_path = Some(parse_option_path(tokens, &mut idx, "--catalog", label)?)
                    }
                    other => {
                        return Err(format!(
                            "Unknown option '{other}' for {label} validate-catalog"
                        ));
                    }
                }
            }
            Ok(ShellCommand::ReferenceValidateCatalog {
                helper_mode,
                catalog_path,
            })
        }
        "status" => {
            if tokens.len() < 3 {
                return Err(format!(
                    "{label} status requires GENOME_ID [--catalog PATH] [--cache-dir PATH]"
                ));
            }
            let genome_id = tokens[2].clone();
            let mut catalog_path: Option<String> = None;
            let mut cache_dir: Option<String> = None;
            let mut idx = 3usize;
            while idx < tokens.len() {
                match tokens[idx].as_str() {
                    "--catalog" => {
                        catalog_path =
                            Some(parse_option_path(tokens, &mut idx, "--catalog", label)?)
                    }
                    "--cache-dir" => {
                        cache_dir = Some(parse_option_path(
                            tokens,
                            &mut idx,
                            "--cache-dir",
                            label,
                        )?)
                    }
                    other => {
                        return Err(format!("Unknown option '{other}' for {label} status"));
                    }
                }
            }
            Ok(ShellCommand::ReferenceStatus {
                helper_mode,
                genome_id,
                catalog_path,
                cache_dir,
            })
        }
        "genes" => {
            if tokens.len() < 3 {
                return Err(format!(
                    "{label} genes requires GENOME_ID [--catalog PATH] [--cache-dir PATH] [--filter REGEX] [--biotype NAME] [--limit N] [--offset N]"
                ));
            }
            let genome_id = tokens[2].clone();
            let mut catalog_path: Option<String> = None;
            let mut cache_dir: Option<String> = None;
            let mut filter = String::new();
            let mut biotypes: Vec<String> = vec![];
            let mut limit: usize = 200;
            let mut offset: usize = 0;
            let mut idx = 3usize;
            while idx < tokens.len() {
                match tokens[idx].as_str() {
                    "--catalog" => {
                        catalog_path =
                            Some(parse_option_path(tokens, &mut idx, "--catalog", label)?)
                    }
                    "--cache-dir" => {
                        cache_dir = Some(parse_option_path(
                            tokens,
                            &mut idx,
                            "--cache-dir",
                            label,
                        )?)
                    }
                    "--filter" => {
                        filter = parse_option_path(tokens, &mut idx, "--filter", label)?
                    }
                    "--biotype" => {
                        let biotype = parse_option_path(tokens, &mut idx, "--biotype", label)?;
                        let trimmed = biotype.trim();
                        if !trimmed.is_empty() {
                            biotypes.push(trimmed.to_string());
                        }
                    }
                    "--limit" => {
                        let raw = parse_option_path(tokens, &mut idx, "--limit", label)?;
                        limit = raw
                            .parse::<usize>()
                            .map_err(|e| format!("Invalid --limit value '{raw}': {e}"))?;
                        if limit == 0 {
                            return Err("--limit must be >= 1".to_string());
                        }
                    }
                    "--offset" => {
                        let raw = parse_option_path(tokens, &mut idx, "--offset", label)?;
                        offset = raw
                            .parse::<usize>()
                            .map_err(|e| format!("Invalid --offset value '{raw}': {e}"))?;
                    }
                    other => {
                        return Err(format!("Unknown option '{other}' for {label} genes"));
                    }
                }
            }
            Ok(ShellCommand::ReferenceGenes {
                helper_mode,
                genome_id,
                catalog_path,
                cache_dir,
                filter,
                biotypes,
                limit,
                offset,
            })
        }
        "prepare" => {
            if tokens.len() < 3 {
                return Err(format!(
                    "{label} prepare requires GENOME_ID [--catalog PATH] [--cache-dir PATH]"
                ));
            }
            let genome_id = tokens[2].clone();
            let mut catalog_path: Option<String> = None;
            let mut cache_dir: Option<String> = None;
            let mut idx = 3usize;
            while idx < tokens.len() {
                match tokens[idx].as_str() {
                    "--catalog" => {
                        catalog_path =
                            Some(parse_option_path(tokens, &mut idx, "--catalog", label)?)
                    }
                    "--cache-dir" => {
                        cache_dir = Some(parse_option_path(
                            tokens,
                            &mut idx,
                            "--cache-dir",
                            label,
                        )?)
                    }
                    other => {
                        return Err(format!("Unknown option '{other}' for {label} prepare"));
                    }
                }
            }
            Ok(ShellCommand::ReferencePrepare {
                helper_mode,
                genome_id,
                catalog_path,
                cache_dir,
            })
        }
        "blast" => {
            if tokens.len() < 4 {
                return Err(format!(
                    "{label} blast requires GENOME_ID QUERY_SEQUENCE [--max-hits N] [--task blastn-short|blastn] [--catalog PATH] [--cache-dir PATH]"
                ));
            }
            let genome_id = tokens[2].clone();
            let query_sequence = tokens[3].clone();
            let mut max_hits: usize = 25;
            let mut task: Option<String> = None;
            let mut catalog_path: Option<String> = None;
            let mut cache_dir: Option<String> = None;
            let mut idx = 4usize;
            while idx < tokens.len() {
                match tokens[idx].as_str() {
                    "--max-hits" => {
                        let raw = parse_option_path(tokens, &mut idx, "--max-hits", label)?;
                        max_hits = raw
                            .parse::<usize>()
                            .map_err(|e| format!("Invalid --max-hits value '{raw}': {e}"))?;
                        if max_hits == 0 {
                            return Err("--max-hits must be >= 1".to_string());
                        }
                    }
                    "--task" => {
                        let raw = parse_option_path(tokens, &mut idx, "--task", label)?;
                        let normalized = raw.trim().to_ascii_lowercase();
                        match normalized.as_str() {
                            "blastn-short" | "blastn" => task = Some(normalized),
                            _ => {
                                return Err(format!(
                                    "Unsupported --task value '{}'; expected blastn-short or blastn",
                                    raw
                                ))
                            }
                        }
                    }
                    "--catalog" => {
                        catalog_path =
                            Some(parse_option_path(tokens, &mut idx, "--catalog", label)?)
                    }
                    "--cache-dir" => {
                        cache_dir = Some(parse_option_path(
                            tokens,
                            &mut idx,
                            "--cache-dir",
                            label,
                        )?)
                    }
                    other => {
                        return Err(format!("Unknown option '{other}' for {label} blast"));
                    }
                }
            }
            Ok(ShellCommand::ReferenceBlast {
                helper_mode,
                genome_id,
                query_sequence,
                max_hits,
                task,
                catalog_path,
                cache_dir,
            })
        }
        "extract-region" => {
            if tokens.len() < 6 {
                return Err(format!(
                    "{label} extract-region requires GENOME_ID CHR START END [--output-id ID] [--catalog PATH] [--cache-dir PATH]"
                ));
            }
            let genome_id = tokens[2].clone();
            let chromosome = tokens[3].clone();
            let start_1based = tokens[4]
                .parse::<usize>()
                .map_err(|e| format!("Invalid START coordinate '{}': {e}", tokens[4]))?;
            let end_1based = tokens[5]
                .parse::<usize>()
                .map_err(|e| format!("Invalid END coordinate '{}': {e}", tokens[5]))?;
            let mut output_id: Option<String> = None;
            let mut catalog_path: Option<String> = None;
            let mut cache_dir: Option<String> = None;
            let mut idx = 6usize;
            while idx < tokens.len() {
                match tokens[idx].as_str() {
                    "--output-id" => {
                        output_id = Some(parse_option_path(
                            tokens,
                            &mut idx,
                            "--output-id",
                            label,
                        )?)
                    }
                    "--catalog" => {
                        catalog_path =
                            Some(parse_option_path(tokens, &mut idx, "--catalog", label)?)
                    }
                    "--cache-dir" => {
                        cache_dir = Some(parse_option_path(
                            tokens,
                            &mut idx,
                            "--cache-dir",
                            label,
                        )?)
                    }
                    other => {
                        return Err(format!(
                            "Unknown option '{other}' for {label} extract-region"
                        ));
                    }
                }
            }
            Ok(ShellCommand::ReferenceExtractRegion {
                helper_mode,
                genome_id,
                chromosome,
                start_1based,
                end_1based,
                output_id,
                catalog_path,
                cache_dir,
            })
        }
        "extract-gene" => {
            if tokens.len() < 4 {
                return Err(format!(
                    "{label} extract-gene requires GENOME_ID QUERY [--occurrence N] [--output-id ID] [--catalog PATH] [--cache-dir PATH]"
                ));
            }
            let genome_id = tokens[2].clone();
            let gene_query = tokens[3].clone();
            let mut occurrence: Option<usize> = None;
            let mut output_id: Option<String> = None;
            let mut catalog_path: Option<String> = None;
            let mut cache_dir: Option<String> = None;
            let mut idx = 4usize;
            while idx < tokens.len() {
                match tokens[idx].as_str() {
                    "--occurrence" => {
                        let raw = parse_option_path(tokens, &mut idx, "--occurrence", label)?;
                        let value = raw
                            .parse::<usize>()
                            .map_err(|e| format!("Invalid --occurrence value '{raw}': {e}"))?;
                        if value == 0 {
                            return Err("--occurrence must be >= 1".to_string());
                        }
                        occurrence = Some(value);
                    }
                    "--output-id" => {
                        output_id = Some(parse_option_path(
                            tokens,
                            &mut idx,
                            "--output-id",
                            label,
                        )?)
                    }
                    "--catalog" => {
                        catalog_path =
                            Some(parse_option_path(tokens, &mut idx, "--catalog", label)?)
                    }
                    "--cache-dir" => {
                        cache_dir = Some(parse_option_path(
                            tokens,
                            &mut idx,
                            "--cache-dir",
                            label,
                        )?)
                    }
                    other => {
                        return Err(format!("Unknown option '{other}' for {label} extract-gene"));
                    }
                }
            }
            Ok(ShellCommand::ReferenceExtractGene {
                helper_mode,
                genome_id,
                gene_query,
                occurrence,
                output_id,
                catalog_path,
                cache_dir,
            })
        }
        other => Err(format!(
            "Unknown {label} subcommand '{other}' (expected list, validate-catalog, status, genes, prepare, blast, extract-region, extract-gene)"
        )),
    }
}

pub fn parse_shell_tokens(tokens: &[String]) -> Result<ShellCommand, String> {
    if tokens.is_empty() {
        return Err("Missing shell command".to_string());
    }
    let cmd = tokens[0].as_str();
    match cmd {
        "help" | "-h" | "--help" => Ok(ShellCommand::Help),
        "capabilities" => {
            if tokens.len() == 1 {
                Ok(ShellCommand::Capabilities)
            } else {
                Err(token_error(cmd))
            }
        }
        "state-summary" => {
            if tokens.len() == 1 {
                Ok(ShellCommand::StateSummary)
            } else {
                Err(token_error(cmd))
            }
        }
        "load-project" | "import-state" => {
            if tokens.len() == 2 {
                Ok(ShellCommand::LoadProject {
                    path: tokens[1].clone(),
                })
            } else {
                Err(token_error(cmd))
            }
        }
        "save-project" | "export-state" => {
            if tokens.len() == 2 {
                Ok(ShellCommand::SaveProject {
                    path: tokens[1].clone(),
                })
            } else {
                Err(token_error(cmd))
            }
        }
        "screenshot-window" => {
            if tokens.len() != 2 {
                return Err(token_error(cmd));
            }
            if !cfg!(feature = "screenshot-capture") {
                return Err(
                    "screenshot-window is unavailable in this build; enable feature 'screenshot-capture'"
                        .to_string(),
                );
            }
            let output = tokens[1].trim();
            if output.is_empty() {
                return Err("screenshot-window requires OUTPUT path".to_string());
            }
            Ok(ShellCommand::ScreenshotWindow {
                output: output.to_string(),
            })
        }
        "render-svg" => {
            if tokens.len() != 4 {
                return Err(token_error(cmd));
            }
            Ok(ShellCommand::RenderSvg {
                seq_id: tokens[1].clone(),
                mode: parse_mode(&tokens[2])?,
                output: tokens[3].clone(),
            })
        }
        "render-rna-svg" => {
            if tokens.len() != 3 {
                return Err(token_error(cmd));
            }
            Ok(ShellCommand::RenderRnaSvg {
                seq_id: tokens[1].clone(),
                output: tokens[2].clone(),
            })
        }
        "rna-info" => {
            if tokens.len() != 2 {
                return Err(token_error(cmd));
            }
            Ok(ShellCommand::RnaInfo {
                seq_id: tokens[1].clone(),
            })
        }
        "render-lineage-svg" => {
            if tokens.len() == 2 {
                Ok(ShellCommand::RenderLineageSvg {
                    output: tokens[1].clone(),
                })
            } else {
                Err(token_error(cmd))
            }
        }
        "render-pool-gel-svg" => {
            if tokens.len() < 3 {
                return Err(token_error(cmd));
            }
            let inputs = split_ids(&tokens[1]);
            if inputs.is_empty() {
                return Err("render-pool-gel-svg requires at least one sequence id".to_string());
            }
            let output = tokens[2].clone();
            let mut ladders: Option<Vec<String>> = None;
            let mut idx = 3usize;
            while idx < tokens.len() {
                match tokens[idx].as_str() {
                    "--ladders" => {
                        if idx + 1 >= tokens.len() {
                            return Err("Missing value after --ladders".to_string());
                        }
                        let parsed = split_ids(&tokens[idx + 1]);
                        if !parsed.is_empty() {
                            ladders = Some(parsed);
                        }
                        idx += 2;
                    }
                    other => {
                        return Err(format!(
                            "Unknown argument '{other}' for render-pool-gel-svg"
                        ));
                    }
                }
            }
            Ok(ShellCommand::RenderPoolGelSvg {
                inputs,
                output,
                ladders,
            })
        }
        "ladders" => {
            if tokens.len() < 2 {
                return Err("ladders requires a subcommand: list or export".to_string());
            }
            match tokens[1].as_str() {
                "list" => {
                    let mut molecule = LadderMolecule::Dna;
                    let mut name_filter: Option<String> = None;
                    let mut idx = 2usize;
                    while idx < tokens.len() {
                        match tokens[idx].as_str() {
                            "--molecule" => {
                                if idx + 1 >= tokens.len() {
                                    return Err("Missing value after --molecule".to_string());
                                }
                                molecule = parse_ladder_molecule(&tokens[idx + 1])?;
                                idx += 2;
                            }
                            "--filter" => {
                                if idx + 1 >= tokens.len() {
                                    return Err("Missing value after --filter".to_string());
                                }
                                name_filter = Some(tokens[idx + 1].clone());
                                idx += 2;
                            }
                            other => {
                                return Err(format!("Unknown argument '{other}' for ladders list"));
                            }
                        }
                    }
                    Ok(ShellCommand::LaddersList {
                        molecule,
                        name_filter,
                    })
                }
                "export" => {
                    if tokens.len() < 3 {
                        return Err(
                            "ladders export requires: OUTPUT.json [--molecule dna|rna] [--filter TEXT]".to_string()
                        );
                    }
                    let output = tokens[2].clone();
                    let mut molecule = LadderMolecule::Dna;
                    let mut name_filter: Option<String> = None;
                    let mut idx = 3usize;
                    while idx < tokens.len() {
                        match tokens[idx].as_str() {
                            "--molecule" => {
                                if idx + 1 >= tokens.len() {
                                    return Err("Missing value after --molecule".to_string());
                                }
                                molecule = parse_ladder_molecule(&tokens[idx + 1])?;
                                idx += 2;
                            }
                            "--filter" => {
                                if idx + 1 >= tokens.len() {
                                    return Err("Missing value after --filter".to_string());
                                }
                                name_filter = Some(tokens[idx + 1].clone());
                                idx += 2;
                            }
                            other => {
                                return Err(format!(
                                    "Unknown argument '{other}' for ladders export"
                                ));
                            }
                        }
                    }
                    Ok(ShellCommand::LaddersExport {
                        molecule,
                        output,
                        name_filter,
                    })
                }
                other => Err(format!(
                    "Unknown ladders subcommand '{other}' (expected list or export)"
                )),
            }
        }
        "export-pool" => {
            if tokens.len() < 3 {
                return Err(token_error(cmd));
            }
            let inputs = split_ids(&tokens[1]);
            if inputs.is_empty() {
                return Err("export-pool requires at least one sequence id".to_string());
            }
            Ok(ShellCommand::ExportPool {
                inputs,
                output: tokens[2].clone(),
                human_id: tokens.get(3).cloned(),
            })
        }
        "import-pool" => {
            if tokens.len() > 3 {
                return Err(token_error(cmd));
            }
            if tokens.len() < 2 {
                return Err("import-pool requires: INPUT.pool.gentle.json [PREFIX]".to_string());
            }
            Ok(ShellCommand::ImportPool {
                input: tokens[1].clone(),
                prefix: tokens.get(2).cloned().unwrap_or_else(|| "pool".to_string()),
            })
        }
        "resources" => {
            if tokens.len() < 2 {
                return Err(
                    "resources requires a subcommand: sync-rebase or sync-jaspar".to_string(),
                );
            }
            match tokens[1].as_str() {
                "sync-rebase" => {
                    if tokens.len() < 3 {
                        return Err(
                            "resources sync-rebase requires INPUT.withrefm_or_URL".to_string()
                        );
                    }
                    let input = tokens[2].clone();
                    let mut output: Option<String> = None;
                    let mut commercial_only = false;
                    let mut idx = 3usize;
                    while idx < tokens.len() {
                        match tokens[idx].as_str() {
                            "--commercial-only" => {
                                commercial_only = true;
                                idx += 1;
                            }
                            value if value.starts_with("--") => {
                                return Err(format!(
                                    "Unknown option '{value}' for resources sync-rebase"
                                ));
                            }
                            value => {
                                if output.is_some() {
                                    return Err(format!(
                                        "Unexpected extra positional argument '{value}' for resources sync-rebase"
                                    ));
                                }
                                output = Some(value.to_string());
                                idx += 1;
                            }
                        }
                    }
                    Ok(ShellCommand::ResourcesSyncRebase {
                        input,
                        output,
                        commercial_only,
                    })
                }
                "sync-jaspar" => {
                    if tokens.len() < 3 {
                        return Err(
                            "resources sync-jaspar requires INPUT.jaspar_or_URL".to_string()
                        );
                    }
                    let input = tokens[2].clone();
                    let mut output: Option<String> = None;
                    let mut idx = 3usize;
                    while idx < tokens.len() {
                        let value = tokens[idx].as_str();
                        if value.starts_with("--") {
                            return Err(format!(
                                "Unknown option '{value}' for resources sync-jaspar"
                            ));
                        }
                        if output.is_some() {
                            return Err(format!(
                                "Unexpected extra positional argument '{value}' for resources sync-jaspar"
                            ));
                        }
                        output = Some(value.to_string());
                        idx += 1;
                    }
                    Ok(ShellCommand::ResourcesSyncJaspar { input, output })
                }
                other => Err(format!(
                    "Unknown resources subcommand '{other}' (expected sync-rebase or sync-jaspar)"
                )),
            }
        }
        "tracks" => {
            if tokens.len() < 2 {
                return Err("tracks requires a subcommand: import-bed or import-bigwig".to_string());
            }
            match tokens[1].as_str() {
                "import-bed" => {
                    if tokens.len() < 4 {
                        return Err(
                            "tracks import-bed requires SEQ_ID PATH [--name NAME] [--min-score N] [--max-score N] [--clear-existing]".to_string()
                        );
                    }
                    let seq_id = tokens[2].clone();
                    let path = tokens[3].clone();
                    let mut track_name: Option<String> = None;
                    let mut min_score: Option<f64> = None;
                    let mut max_score: Option<f64> = None;
                    let mut clear_existing = false;
                    let mut idx = 4usize;
                    while idx < tokens.len() {
                        match tokens[idx].as_str() {
                            "--name" => {
                                track_name = Some(parse_option_path(
                                    tokens,
                                    &mut idx,
                                    "--name",
                                    "tracks import-bed",
                                )?);
                            }
                            "--min-score" => {
                                let raw = parse_option_path(
                                    tokens,
                                    &mut idx,
                                    "--min-score",
                                    "tracks import-bed",
                                )?;
                                min_score = Some(raw.parse::<f64>().map_err(|e| {
                                    format!("Invalid --min-score value '{raw}': {e}")
                                })?);
                            }
                            "--max-score" => {
                                let raw = parse_option_path(
                                    tokens,
                                    &mut idx,
                                    "--max-score",
                                    "tracks import-bed",
                                )?;
                                max_score = Some(raw.parse::<f64>().map_err(|e| {
                                    format!("Invalid --max-score value '{raw}': {e}")
                                })?);
                            }
                            "--clear-existing" => {
                                clear_existing = true;
                                idx += 1;
                            }
                            other => {
                                return Err(format!(
                                    "Unknown option '{other}' for tracks import-bed"
                                ));
                            }
                        }
                    }
                    if min_score
                        .zip(max_score)
                        .map(|(min, max)| min > max)
                        .unwrap_or(false)
                    {
                        return Err("--min-score must be <= --max-score".to_string());
                    }
                    Ok(ShellCommand::TracksImportBed {
                        seq_id,
                        path,
                        track_name,
                        min_score,
                        max_score,
                        clear_existing,
                    })
                }
                "import-bigwig" => {
                    if tokens.len() < 4 {
                        return Err(
                            "tracks import-bigwig requires SEQ_ID PATH [--name NAME] [--min-score N] [--max-score N] [--clear-existing]".to_string()
                        );
                    }
                    let seq_id = tokens[2].clone();
                    let path = tokens[3].clone();
                    let mut track_name: Option<String> = None;
                    let mut min_score: Option<f64> = None;
                    let mut max_score: Option<f64> = None;
                    let mut clear_existing = false;
                    let mut idx = 4usize;
                    while idx < tokens.len() {
                        match tokens[idx].as_str() {
                            "--name" => {
                                track_name = Some(parse_option_path(
                                    tokens,
                                    &mut idx,
                                    "--name",
                                    "tracks import-bigwig",
                                )?);
                            }
                            "--min-score" => {
                                let raw = parse_option_path(
                                    tokens,
                                    &mut idx,
                                    "--min-score",
                                    "tracks import-bigwig",
                                )?;
                                min_score = Some(raw.parse::<f64>().map_err(|e| {
                                    format!("Invalid --min-score value '{raw}': {e}")
                                })?);
                            }
                            "--max-score" => {
                                let raw = parse_option_path(
                                    tokens,
                                    &mut idx,
                                    "--max-score",
                                    "tracks import-bigwig",
                                )?;
                                max_score = Some(raw.parse::<f64>().map_err(|e| {
                                    format!("Invalid --max-score value '{raw}': {e}")
                                })?);
                            }
                            "--clear-existing" => {
                                clear_existing = true;
                                idx += 1;
                            }
                            other => {
                                return Err(format!(
                                    "Unknown option '{other}' for tracks import-bigwig"
                                ));
                            }
                        }
                    }
                    if min_score
                        .zip(max_score)
                        .map(|(min, max)| min > max)
                        .unwrap_or(false)
                    {
                        return Err("--min-score must be <= --max-score".to_string());
                    }
                    Ok(ShellCommand::TracksImportBigWig {
                        seq_id,
                        path,
                        track_name,
                        min_score,
                        max_score,
                        clear_existing,
                    })
                }
                other => Err(format!(
                    "Unknown tracks subcommand '{other}' (expected import-bed or import-bigwig)"
                )),
            }
        }
        "genomes" => parse_reference_command(tokens, false),
        "helpers" => parse_reference_command(tokens, true),
        "op" => {
            let payload = tokens[1..].join(" ");
            if payload.trim().is_empty() {
                return Err("Missing operation JSON".to_string());
            }
            Ok(ShellCommand::Op { payload })
        }
        "workflow" => {
            let payload = tokens[1..].join(" ");
            if payload.trim().is_empty() {
                return Err("Missing workflow JSON".to_string());
            }
            Ok(ShellCommand::Workflow { payload })
        }
        other => Err(format!("Unknown shell command '{other}'. Try: help")),
    }
}

pub fn parse_shell_line(line: &str) -> Result<ShellCommand, String> {
    let tokens = split_shell_words(line)?;
    parse_shell_tokens(&tokens)
}

pub fn split_shell_words(line: &str) -> Result<Vec<String>, String> {
    #[derive(Clone, Copy, PartialEq, Eq)]
    enum Mode {
        Normal,
        SingleQuoted,
        DoubleQuoted,
    }

    let mut out = Vec::new();
    let mut current = String::new();
    let mut mode = Mode::Normal;
    let mut chars = line.chars().peekable();

    while let Some(ch) = chars.next() {
        match mode {
            Mode::Normal => match ch {
                '\'' => mode = Mode::SingleQuoted,
                '"' => mode = Mode::DoubleQuoted,
                '\\' => {
                    if let Some(next) = chars.next() {
                        current.push(next);
                    }
                }
                c if c.is_whitespace() => {
                    if !current.is_empty() {
                        out.push(current.clone());
                        current.clear();
                    }
                }
                _ => current.push(ch),
            },
            Mode::SingleQuoted => {
                if ch == '\'' {
                    mode = Mode::Normal;
                } else {
                    current.push(ch);
                }
            }
            Mode::DoubleQuoted => {
                if ch == '"' {
                    mode = Mode::Normal;
                } else if ch == '\\' {
                    if let Some(next) = chars.next() {
                        current.push(next);
                    }
                } else {
                    current.push(ch);
                }
            }
        }
    }

    if mode != Mode::Normal {
        return Err("Unterminated quoted string in shell command".to_string());
    }
    if !current.is_empty() {
        out.push(current);
    }
    if out.is_empty() {
        return Err("Empty shell command".to_string());
    }
    Ok(out)
}

pub fn execute_shell_command(
    engine: &mut GentleEngine,
    command: &ShellCommand,
) -> Result<ShellRunResult, String> {
    execute_shell_command_with_options(engine, command, &ShellExecutionOptions::default())
}

pub fn execute_shell_command_with_options(
    engine: &mut GentleEngine,
    command: &ShellCommand,
    options: &ShellExecutionOptions,
) -> Result<ShellRunResult, String> {
    let result = match command {
        ShellCommand::Help => ShellRunResult {
            state_changed: false,
            output: json!({ "help": shell_help_text() }),
        },
        ShellCommand::Capabilities => ShellRunResult {
            state_changed: false,
            output: serde_json::to_value(GentleEngine::capabilities())
                .map_err(|e| format!("Could not serialize capabilities: {e}"))?,
        },
        ShellCommand::StateSummary => ShellRunResult {
            state_changed: false,
            output: serde_json::to_value(engine.summarize_state())
                .map_err(|e| format!("Could not serialize state summary: {e}"))?,
        },
        ShellCommand::LoadProject { path } => {
            let state = ProjectState::load_from_path(path).map_err(|e| e.to_string())?;
            *engine = GentleEngine::from_state(state);
            ShellRunResult {
                state_changed: true,
                output: json!({
                    "message": format!("Loaded project from '{path}'"),
                    "summary": engine.summarize_state()
                }),
            }
        }
        ShellCommand::SaveProject { path } => {
            engine
                .state()
                .save_to_path(path)
                .map_err(|e| e.to_string())?;
            ShellRunResult {
                state_changed: false,
                output: json!({ "message": format!("Saved project to '{path}'") }),
            }
        }
        ShellCommand::ScreenshotWindow { output } => {
            if !cfg!(feature = "screenshot-capture") {
                return Err(
                    "screenshot capture is unavailable in this build; enable feature 'screenshot-capture'"
                        .to_string(),
                );
            }
            if !options.allow_screenshots {
                return Err(
                    "screenshot capture is disabled; restart with --allow-screenshots".to_string(),
                );
            }
            let report = capture_active_window_screenshot(output)?;
            ShellRunResult {
                state_changed: false,
                output: serde_json::to_value(report)
                    .map_err(|e| format!("Could not serialize screenshot report: {e}"))?,
            }
        }
        ShellCommand::RenderSvg {
            seq_id,
            mode,
            output,
        } => {
            let op_result = engine
                .apply(Operation::RenderSequenceSvg {
                    seq_id: seq_id.clone(),
                    mode: mode.clone(),
                    path: output.clone(),
                })
                .map_err(|e| e.to_string())?;
            ShellRunResult {
                state_changed: false,
                output: json!({ "result": op_result }),
            }
        }
        ShellCommand::RenderRnaSvg { seq_id, output } => {
            let op_result = engine
                .apply(Operation::RenderRnaStructureSvg {
                    seq_id: seq_id.clone(),
                    path: output.clone(),
                })
                .map_err(|e| e.to_string())?;
            ShellRunResult {
                state_changed: false,
                output: json!({ "result": op_result }),
            }
        }
        ShellCommand::RnaInfo { seq_id } => {
            let report = engine
                .inspect_rna_structure(seq_id)
                .map_err(|e| e.to_string())?;
            ShellRunResult {
                state_changed: false,
                output: serde_json::to_value(report)
                    .map_err(|e| format!("Could not serialize RNA report: {e}"))?,
            }
        }
        ShellCommand::RenderLineageSvg { output } => {
            let op_result = engine
                .apply(Operation::RenderLineageSvg {
                    path: output.clone(),
                })
                .map_err(|e| e.to_string())?;
            ShellRunResult {
                state_changed: false,
                output: json!({ "result": op_result }),
            }
        }
        ShellCommand::RenderPoolGelSvg {
            inputs,
            output,
            ladders,
        } => {
            let op_result = engine
                .apply(Operation::RenderPoolGelSvg {
                    inputs: inputs.clone(),
                    path: output.clone(),
                    ladders: ladders.clone(),
                })
                .map_err(|e| e.to_string())?;
            ShellRunResult {
                state_changed: false,
                output: json!({ "result": op_result }),
            }
        }
        ShellCommand::LaddersList {
            molecule,
            name_filter,
        } => ShellRunResult {
            state_changed: false,
            output: match molecule {
                LadderMolecule::Dna => {
                    serde_json::to_value(GentleEngine::inspect_dna_ladders(name_filter.as_deref()))
                        .map_err(|e| format!("Could not serialize DNA ladders catalog: {e}"))?
                }
                LadderMolecule::Rna => {
                    serde_json::to_value(GentleEngine::inspect_rna_ladders(name_filter.as_deref()))
                        .map_err(|e| format!("Could not serialize RNA ladders catalog: {e}"))?
                }
            },
        },
        ShellCommand::LaddersExport {
            molecule,
            output,
            name_filter,
        } => {
            let op_result = match molecule {
                LadderMolecule::Dna => engine
                    .apply(Operation::ExportDnaLadders {
                        path: output.clone(),
                        name_filter: name_filter.clone(),
                    })
                    .map_err(|e| e.to_string())?,
                LadderMolecule::Rna => engine
                    .apply(Operation::ExportRnaLadders {
                        path: output.clone(),
                        name_filter: name_filter.clone(),
                    })
                    .map_err(|e| e.to_string())?,
            };
            ShellRunResult {
                state_changed: false,
                output: json!({ "result": op_result }),
            }
        }
        ShellCommand::ExportPool {
            inputs,
            output,
            human_id,
        } => {
            let op_result = engine
                .apply(Operation::ExportPool {
                    inputs: inputs.clone(),
                    path: output.clone(),
                    pool_id: Some("pool_export".to_string()),
                    human_id: human_id.clone(),
                })
                .map_err(|e| e.to_string())?;
            ShellRunResult {
                state_changed: false,
                output: json!({ "result": op_result }),
            }
        }
        ShellCommand::ImportPool { input, prefix } => {
            let text = fs::read_to_string(input)
                .map_err(|e| format!("Could not read pool file '{input}': {e}"))?;
            let pool: PoolExport = serde_json::from_str(&text)
                .map_err(|e| format!("Invalid pool JSON '{input}': {e}"))?;
            if pool.schema != "gentle.pool.v1" {
                return Err(format!(
                    "Unsupported pool schema '{}', expected 'gentle.pool.v1'",
                    pool.schema
                ));
            }

            let mut state = engine.state().clone();
            let mut imported_ids = Vec::new();
            for (idx, member) in pool.members.iter().enumerate() {
                let mut dna = DNAsequence::from_sequence(&member.sequence)
                    .map_err(|e| format!("Invalid DNA in pool member '{}': {e}", member.seq_id))?;
                if let Some(name) = &member.name {
                    let mut value = serde_json::to_value(&dna)
                        .map_err(|e| format!("Could not serialize sequence: {e}"))?;
                    if let Some(obj) = value.as_object_mut() {
                        if let Some(seq_obj) = obj.get_mut("seq").and_then(|v| v.as_object_mut()) {
                            seq_obj.insert("name".to_string(), json!(name));
                        }
                    }
                    dna = serde_json::from_value(value)
                        .map_err(|e| format!("Could not set sequence name: {e}"))?;
                }
                if member.topology.eq_ignore_ascii_case("circular") {
                    dna.set_circular(true);
                }
                apply_member_overhang(member, &mut dna)?;
                dna.update_computed_features();
                let base = format!("{prefix}_{}", idx + 1);
                let id = unique_id(&state.sequences, &base);
                state.sequences.insert(id.clone(), dna);
                imported_ids.push(id);
            }
            *engine = GentleEngine::from_state(state);
            ShellRunResult {
                state_changed: true,
                output: json!({
                    "message": format!("Imported pool '{}' ({} members)", pool.pool_id, pool.member_count),
                    "input": input,
                    "pool_id": pool.pool_id,
                    "member_count": pool.member_count,
                    "imported_ids": imported_ids,
                }),
            }
        }
        ShellCommand::ResourcesSyncRebase {
            input,
            output,
            commercial_only,
        } => {
            let report = resource_sync::sync_rebase(input, output.as_deref(), *commercial_only)?;
            ShellRunResult {
                state_changed: false,
                output: json!({
                    "message": format!("Synced {} {} entries to '{}'", report.item_count, report.resource, report.output),
                    "report": report,
                }),
            }
        }
        ShellCommand::ResourcesSyncJaspar { input, output } => {
            let report = resource_sync::sync_jaspar(input, output.as_deref())?;
            ShellRunResult {
                state_changed: false,
                output: json!({
                    "message": format!("Synced {} {} entries to '{}'", report.item_count, report.resource, report.output),
                    "report": report,
                }),
            }
        }
        ShellCommand::ReferenceList {
            helper_mode,
            catalog_path,
        } => {
            let resolved_catalog = resolved_catalog_path(catalog_path, *helper_mode);
            let genomes = GentleEngine::list_reference_genomes(resolved_catalog)
                .map_err(|e| e.to_string())?;
            let effective_catalog = effective_catalog_path(catalog_path, *helper_mode);
            ShellRunResult {
                state_changed: false,
                output: json!({
                    "catalog_path": effective_catalog,
                    "genome_count": genomes.len(),
                    "genomes": genomes,
                }),
            }
        }
        ShellCommand::ReferenceValidateCatalog {
            helper_mode,
            catalog_path,
        } => {
            let resolved_catalog = resolved_catalog_path(catalog_path, *helper_mode);
            let effective_catalog = effective_catalog_path(catalog_path, *helper_mode);
            let genomes = GentleEngine::list_reference_genomes(resolved_catalog)
                .map_err(|e| e.to_string())?;
            for genome_id in &genomes {
                GentleEngine::describe_reference_genome_sources(resolved_catalog, genome_id, None)
                    .map_err(|e| e.to_string())?;
            }
            ShellRunResult {
                state_changed: false,
                output: json!({
                    "catalog_path": effective_catalog,
                    "valid": true,
                    "genome_count": genomes.len(),
                    "validated_sources": genomes.len(),
                    "genomes": genomes,
                }),
            }
        }
        ShellCommand::ReferenceStatus {
            helper_mode,
            genome_id,
            catalog_path,
            cache_dir,
        } => {
            let resolved_catalog = resolved_catalog_path(catalog_path, *helper_mode);
            let prepared = GentleEngine::is_reference_genome_prepared(
                resolved_catalog,
                genome_id,
                cache_dir.as_deref(),
            )
            .map_err(|e| e.to_string())?;
            let source_plan = GentleEngine::describe_reference_genome_sources(
                resolved_catalog,
                genome_id,
                cache_dir.as_deref(),
            )
            .map_err(|e| e.to_string())?;
            let effective_catalog = effective_catalog_path(catalog_path, *helper_mode);
            ShellRunResult {
                state_changed: false,
                output: json!({
                    "genome_id": genome_id,
                    "catalog_path": effective_catalog,
                    "cache_dir": cache_dir,
                    "prepared": prepared,
                    "sequence_source_type": source_plan.sequence_source_type,
                    "annotation_source_type": source_plan.annotation_source_type,
                    "sequence_source": source_plan.sequence_source,
                    "annotation_source": source_plan.annotation_source,
                }),
            }
        }
        ShellCommand::ReferenceGenes {
            helper_mode,
            genome_id,
            catalog_path,
            cache_dir,
            filter,
            biotypes,
            limit,
            offset,
        } => {
            let resolved_catalog = resolved_catalog_path(catalog_path, *helper_mode);
            let genes = GentleEngine::list_reference_genome_genes(
                resolved_catalog,
                genome_id,
                cache_dir.as_deref(),
            )
            .map_err(|e| e.to_string())?;
            let filter_regex = compile_gene_filter_regex(filter)?;
            let biotype_filter: Vec<String> = biotypes
                .iter()
                .map(|v| v.trim().to_ascii_lowercase())
                .filter(|v| !v.is_empty())
                .collect();
            let available_biotypes = collect_biotypes(&genes);
            let filtered: Vec<GenomeGeneRecord> = genes
                .into_iter()
                .filter(|g| genome_gene_matches_filter(g, filter_regex.as_ref(), &biotype_filter))
                .collect();
            let total = filtered.len();
            let clamped_offset = (*offset).min(total);
            let returned: Vec<GenomeGeneRecord> = filtered
                .into_iter()
                .skip(clamped_offset)
                .take(*limit)
                .collect();
            let effective_catalog = effective_catalog_path(catalog_path, *helper_mode);
            ShellRunResult {
                state_changed: false,
                output: json!({
                    "genome_id": genome_id,
                    "catalog_path": effective_catalog,
                    "cache_dir": cache_dir,
                    "filter": filter,
                    "biotype_filter": biotypes,
                    "available_biotypes": available_biotypes,
                    "offset": clamped_offset,
                    "limit": limit,
                    "total_matches": total,
                    "returned": returned.len(),
                    "genes": returned,
                }),
            }
        }
        ShellCommand::ReferencePrepare {
            helper_mode,
            genome_id,
            catalog_path,
            cache_dir,
        } => {
            let op_result = engine
                .apply(Operation::PrepareGenome {
                    genome_id: genome_id.clone(),
                    catalog_path: operation_catalog_path(catalog_path, *helper_mode),
                    cache_dir: cache_dir.clone(),
                })
                .map_err(|e| e.to_string())?;
            ShellRunResult {
                state_changed: true,
                output: json!({ "result": op_result }),
            }
        }
        ShellCommand::ReferenceBlast {
            helper_mode,
            genome_id,
            query_sequence,
            max_hits,
            task,
            catalog_path,
            cache_dir,
        } => {
            let resolved_catalog = resolved_catalog_path(catalog_path, *helper_mode);
            let report = if *helper_mode {
                GentleEngine::blast_helper_genome(
                    genome_id,
                    query_sequence,
                    *max_hits,
                    task.as_deref(),
                    resolved_catalog,
                    cache_dir.as_deref(),
                )
            } else {
                GentleEngine::blast_reference_genome(
                    resolved_catalog,
                    genome_id,
                    query_sequence,
                    *max_hits,
                    task.as_deref(),
                    cache_dir.as_deref(),
                )
            }
            .map_err(|e| e.to_string())?;
            let effective_catalog = effective_catalog_path(catalog_path, *helper_mode);
            ShellRunResult {
                state_changed: false,
                output: json!({
                    "catalog_path": effective_catalog,
                    "cache_dir": cache_dir,
                    "report": report,
                }),
            }
        }
        ShellCommand::ReferenceExtractRegion {
            helper_mode,
            genome_id,
            chromosome,
            start_1based,
            end_1based,
            output_id,
            catalog_path,
            cache_dir,
        } => {
            let op_result = engine
                .apply(Operation::ExtractGenomeRegion {
                    genome_id: genome_id.clone(),
                    chromosome: chromosome.clone(),
                    start_1based: *start_1based,
                    end_1based: *end_1based,
                    output_id: output_id.clone(),
                    catalog_path: operation_catalog_path(catalog_path, *helper_mode),
                    cache_dir: cache_dir.clone(),
                })
                .map_err(|e| e.to_string())?;
            let state_changed =
                !op_result.created_seq_ids.is_empty() || !op_result.changed_seq_ids.is_empty();
            ShellRunResult {
                state_changed,
                output: json!({ "result": op_result }),
            }
        }
        ShellCommand::ReferenceExtractGene {
            helper_mode,
            genome_id,
            gene_query,
            occurrence,
            output_id,
            catalog_path,
            cache_dir,
        } => {
            let op_result = engine
                .apply(Operation::ExtractGenomeGene {
                    genome_id: genome_id.clone(),
                    gene_query: gene_query.clone(),
                    occurrence: *occurrence,
                    output_id: output_id.clone(),
                    catalog_path: operation_catalog_path(catalog_path, *helper_mode),
                    cache_dir: cache_dir.clone(),
                })
                .map_err(|e| e.to_string())?;
            let state_changed =
                !op_result.created_seq_ids.is_empty() || !op_result.changed_seq_ids.is_empty();
            ShellRunResult {
                state_changed,
                output: json!({ "result": op_result }),
            }
        }
        ShellCommand::TracksImportBed {
            seq_id,
            path,
            track_name,
            min_score,
            max_score,
            clear_existing,
        } => {
            let op_result = engine
                .apply(Operation::ImportGenomeBedTrack {
                    seq_id: seq_id.clone(),
                    path: path.clone(),
                    track_name: track_name.clone(),
                    min_score: *min_score,
                    max_score: *max_score,
                    clear_existing: Some(*clear_existing),
                })
                .map_err(|e| e.to_string())?;
            let state_changed =
                !op_result.created_seq_ids.is_empty() || !op_result.changed_seq_ids.is_empty();
            ShellRunResult {
                state_changed,
                output: json!({ "result": op_result }),
            }
        }
        ShellCommand::TracksImportBigWig {
            seq_id,
            path,
            track_name,
            min_score,
            max_score,
            clear_existing,
        } => {
            let op_result = engine
                .apply(Operation::ImportGenomeBigWigTrack {
                    seq_id: seq_id.clone(),
                    path: path.clone(),
                    track_name: track_name.clone(),
                    min_score: *min_score,
                    max_score: *max_score,
                    clear_existing: Some(*clear_existing),
                })
                .map_err(|e| e.to_string())?;
            let state_changed =
                !op_result.created_seq_ids.is_empty() || !op_result.changed_seq_ids.is_empty();
            ShellRunResult {
                state_changed,
                output: json!({ "result": op_result }),
            }
        }
        ShellCommand::Op { payload } => {
            let json_text = parse_json_payload(payload)?;
            let op: Operation = serde_json::from_str(&json_text)
                .map_err(|e| format!("Invalid operation JSON: {e}"))?;
            let op_result = engine.apply(op).map_err(|e| e.to_string())?;
            let state_changed =
                !op_result.created_seq_ids.is_empty() || !op_result.changed_seq_ids.is_empty();
            ShellRunResult {
                state_changed,
                output: json!({ "result": op_result }),
            }
        }
        ShellCommand::Workflow { payload } => {
            let json_text = parse_json_payload(payload)?;
            let workflow: Workflow = serde_json::from_str(&json_text)
                .map_err(|e| format!("Invalid workflow JSON: {e}"))?;
            let results = engine.apply_workflow(workflow).map_err(|e| e.to_string())?;
            let state_changed = results
                .iter()
                .any(|r| !r.created_seq_ids.is_empty() || !r.changed_seq_ids.is_empty());
            ShellRunResult {
                state_changed,
                output: json!({ "results": results }),
            }
        }
    };
    Ok(result)
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::fs;
    use tempfile::tempdir;

    #[test]
    fn parse_workflow_payload_keeps_whitespace() {
        let cmd = parse_shell_line("workflow '{ \"run_id\": \"x\", \"ops\": [] }'")
            .expect("workflow command parse");
        match cmd {
            ShellCommand::Workflow { payload } => {
                assert!(payload.contains("\"run_id\""));
                assert!(payload.contains("\"ops\""));
            }
            other => panic!("unexpected command: {other:?}"),
        }
    }

    #[test]
    fn parse_render_pool_gel_with_ladders() {
        let cmd = parse_shell_line("render-pool-gel-svg a,b out.svg --ladders 1kb,100bp")
            .expect("parse command");
        match cmd {
            ShellCommand::RenderPoolGelSvg {
                inputs,
                output,
                ladders,
            } => {
                assert_eq!(inputs, vec!["a".to_string(), "b".to_string()]);
                assert_eq!(output, "out.svg".to_string());
                assert_eq!(ladders, Some(vec!["1kb".to_string(), "100bp".to_string()]));
            }
            other => panic!("unexpected command: {other:?}"),
        }
    }

    #[test]
    fn parse_render_rna_svg() {
        let cmd = parse_shell_line("render-rna-svg rna_seq rna.svg").expect("parse command");
        match cmd {
            ShellCommand::RenderRnaSvg { seq_id, output } => {
                assert_eq!(seq_id, "rna_seq".to_string());
                assert_eq!(output, "rna.svg".to_string());
            }
            other => panic!("unexpected command: {other:?}"),
        }
    }

    #[test]
    fn parse_screenshot_window() {
        #[cfg(feature = "screenshot-capture")]
        {
            let cmd =
                parse_shell_line("screenshot-window docs/images/main.png").expect("parse command");
            match cmd {
                ShellCommand::ScreenshotWindow { output } => {
                    assert_eq!(output, "docs/images/main.png".to_string());
                }
                other => panic!("unexpected command: {other:?}"),
            }
        }
        #[cfg(not(feature = "screenshot-capture"))]
        {
            let err = parse_shell_line("screenshot-window docs/images/main.png").unwrap_err();
            assert!(err.contains("unavailable in this build"));
        }
    }

    #[test]
    fn parse_rna_info() {
        let cmd = parse_shell_line("rna-info rna_seq").expect("parse command");
        match cmd {
            ShellCommand::RnaInfo { seq_id } => {
                assert_eq!(seq_id, "rna_seq".to_string());
            }
            other => panic!("unexpected command: {other:?}"),
        }
    }

    #[test]
    fn parse_ladders_list_with_filter() {
        let cmd = parse_shell_line("ladders list --filter NEB").expect("parse command");
        match cmd {
            ShellCommand::LaddersList {
                molecule,
                name_filter,
            } => {
                assert_eq!(molecule, LadderMolecule::Dna);
                assert_eq!(name_filter, Some("NEB".to_string()));
            }
            other => panic!("unexpected command: {other:?}"),
        }
    }

    #[test]
    fn parse_ladders_export_with_filter() {
        let cmd =
            parse_shell_line("ladders export ladders.json --filter ruler").expect("parse command");
        match cmd {
            ShellCommand::LaddersExport {
                molecule,
                output,
                name_filter,
            } => {
                assert_eq!(molecule, LadderMolecule::Dna);
                assert_eq!(output, "ladders.json".to_string());
                assert_eq!(name_filter, Some("ruler".to_string()));
            }
            other => panic!("unexpected command: {other:?}"),
        }
    }

    #[test]
    fn parse_ladders_list_rna() {
        let cmd = parse_shell_line("ladders list --molecule rna --filter ss").expect("parse");
        match cmd {
            ShellCommand::LaddersList {
                molecule,
                name_filter,
            } => {
                assert_eq!(molecule, LadderMolecule::Rna);
                assert_eq!(name_filter, Some("ss".to_string()));
            }
            other => panic!("unexpected command: {other:?}"),
        }
    }

    #[test]
    fn parse_reference_genes_with_regex_and_biotypes() {
        let cmd = parse_shell_line(
            "helpers genes Helper --filter '^bla$' --biotype promoter --biotype cds --limit 10 --offset 3",
        )
        .expect("parse command");
        match cmd {
            ShellCommand::ReferenceGenes {
                helper_mode,
                genome_id,
                filter,
                biotypes,
                limit,
                offset,
                ..
            } => {
                assert!(helper_mode);
                assert_eq!(genome_id, "Helper");
                assert_eq!(filter, "^bla$");
                assert_eq!(biotypes, vec!["promoter".to_string(), "cds".to_string()]);
                assert_eq!(limit, 10);
                assert_eq!(offset, 3);
            }
            other => panic!("unexpected command: {other:?}"),
        }
    }

    #[test]
    fn parse_genomes_blast_with_options() {
        let cmd = parse_shell_line(
            "genomes blast ToyGenome ACGTACGT --max-hits 12 --task blastn --catalog c.json --cache-dir cache",
        )
        .expect("parse command");
        match cmd {
            ShellCommand::ReferenceBlast {
                helper_mode,
                genome_id,
                query_sequence,
                max_hits,
                task,
                catalog_path,
                cache_dir,
            } => {
                assert!(!helper_mode);
                assert_eq!(genome_id, "ToyGenome");
                assert_eq!(query_sequence, "ACGTACGT");
                assert_eq!(max_hits, 12);
                assert_eq!(task.as_deref(), Some("blastn"));
                assert_eq!(catalog_path.as_deref(), Some("c.json"));
                assert_eq!(cache_dir.as_deref(), Some("cache"));
            }
            other => panic!("unexpected command: {other:?}"),
        }
    }

    #[test]
    fn parse_helpers_blast_defaults() {
        let cmd = parse_shell_line("helpers blast pUC19 ACGTAG").expect("parse command");
        match cmd {
            ShellCommand::ReferenceBlast {
                helper_mode,
                genome_id,
                query_sequence,
                max_hits,
                task,
                ..
            } => {
                assert!(helper_mode);
                assert_eq!(genome_id, "pUC19");
                assert_eq!(query_sequence, "ACGTAG");
                assert_eq!(max_hits, 25);
                assert!(task.is_none());
            }
            other => panic!("unexpected command: {other:?}"),
        }
    }

    #[test]
    fn parse_tracks_import_bed() {
        let cmd = parse_shell_line(
            "tracks import-bed toy_slice test_files/data/peaks.bed.gz --name ChIP --min-score 5 --max-score 50 --clear-existing",
        )
        .expect("parse command");
        match cmd {
            ShellCommand::TracksImportBed {
                seq_id,
                path,
                track_name,
                min_score,
                max_score,
                clear_existing,
            } => {
                assert_eq!(seq_id, "toy_slice".to_string());
                assert_eq!(path, "test_files/data/peaks.bed.gz".to_string());
                assert_eq!(track_name, Some("ChIP".to_string()));
                assert_eq!(min_score, Some(5.0));
                assert_eq!(max_score, Some(50.0));
                assert!(clear_existing);
            }
            other => panic!("unexpected command: {other:?}"),
        }
    }

    #[test]
    fn parse_tracks_import_bigwig() {
        let cmd = parse_shell_line(
            "tracks import-bigwig toy_slice test_files/data/signal.bw --name RNA --min-score 0.5 --max-score 2.5 --clear-existing",
        )
        .expect("parse command");
        match cmd {
            ShellCommand::TracksImportBigWig {
                seq_id,
                path,
                track_name,
                min_score,
                max_score,
                clear_existing,
            } => {
                assert_eq!(seq_id, "toy_slice".to_string());
                assert_eq!(path, "test_files/data/signal.bw".to_string());
                assert_eq!(track_name, Some("RNA".to_string()));
                assert_eq!(min_score, Some(0.5));
                assert_eq!(max_score, Some(2.5));
                assert!(clear_existing);
            }
            other => panic!("unexpected command: {other:?}"),
        }
    }

    #[test]
    fn execute_state_summary_returns_json() {
        let mut engine = GentleEngine::from_state(ProjectState::default());
        let out = execute_shell_command(&mut engine, &ShellCommand::StateSummary)
            .expect("execute state summary");
        assert!(!out.state_changed);
        assert!(out.output.get("sequence_count").is_some());
    }

    #[test]
    fn execute_screenshot_requires_allow_flag() {
        let mut engine = GentleEngine::from_state(ProjectState::default());
        let err = execute_shell_command(
            &mut engine,
            &ShellCommand::ScreenshotWindow {
                output: "out.png".to_string(),
            },
        )
        .unwrap_err();
        if cfg!(feature = "screenshot-capture") {
            assert!(err.contains("--allow-screenshots"));
        } else {
            assert!(err.contains("unavailable in this build"));
        }
    }

    #[test]
    fn parse_resources_sync_rebase_with_flag() {
        let cmd = parse_shell_line(
            "resources sync-rebase https://example.org/rebase.withrefm out.json --commercial-only",
        )
        .expect("parse resources sync-rebase");
        match cmd {
            ShellCommand::ResourcesSyncRebase {
                input,
                output,
                commercial_only,
            } => {
                assert_eq!(input, "https://example.org/rebase.withrefm".to_string());
                assert_eq!(output, Some("out.json".to_string()));
                assert!(commercial_only);
            }
            other => panic!("unexpected command: {other:?}"),
        }
    }

    #[test]
    fn parse_helpers_extract_gene() {
        let cmd = parse_shell_line(
            "helpers extract-gene pUC19 bla --occurrence 2 --output-id out --cache-dir cache",
        )
        .expect("parse helpers extract-gene");
        match cmd {
            ShellCommand::ReferenceExtractGene {
                helper_mode,
                genome_id,
                gene_query,
                occurrence,
                output_id,
                cache_dir,
                ..
            } => {
                assert!(helper_mode);
                assert_eq!(genome_id, "pUC19".to_string());
                assert_eq!(gene_query, "bla".to_string());
                assert_eq!(occurrence, Some(2));
                assert_eq!(output_id, Some("out".to_string()));
                assert_eq!(cache_dir, Some("cache".to_string()));
            }
            other => panic!("unexpected command: {other:?}"),
        }
    }

    #[test]
    fn parse_genomes_validate_catalog() {
        let cmd = parse_shell_line("genomes validate-catalog --catalog assets/genomes.json")
            .expect("parse command");
        match cmd {
            ShellCommand::ReferenceValidateCatalog {
                helper_mode,
                catalog_path,
            } => {
                assert!(!helper_mode);
                assert_eq!(catalog_path, Some("assets/genomes.json".to_string()));
            }
            other => panic!("unexpected command: {other:?}"),
        }
    }

    #[test]
    fn execute_genomes_validate_catalog_reports_valid() {
        let td = tempdir().expect("tempdir");
        let fasta = td.path().join("toy.fa");
        let gtf = td.path().join("toy.gtf");
        let cache = td.path().join("cache");
        fs::write(&fasta, ">chr1\nACGT\n").expect("write fasta");
        fs::write(
            &gtf,
            "chr1\tsrc\tgene\t1\t4\t.\t+\t.\tgene_id \"GENE1\"; gene_name \"GENE1\";\n",
        )
        .expect("write gtf");
        let catalog = td.path().join("catalog.json");
        let catalog_json = format!(
            r#"{{
  "ToyGenome": {{
    "sequence_local": "{}",
    "annotations_local": "{}",
    "cache_dir": "{}"
  }}
}}"#,
            fasta.display(),
            gtf.display(),
            cache.display()
        );
        fs::write(&catalog, catalog_json).expect("write catalog");
        let mut engine = GentleEngine::new();
        let out = execute_shell_command(
            &mut engine,
            &ShellCommand::ReferenceValidateCatalog {
                helper_mode: false,
                catalog_path: Some(catalog.to_string_lossy().to_string()),
            },
        )
        .expect("execute validate-catalog");
        assert!(!out.state_changed);
        assert_eq!(out.output["valid"].as_bool(), Some(true));
        assert_eq!(out.output["genome_count"].as_u64(), Some(1));
    }
}
