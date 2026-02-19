use gentle::{
    about,
    engine::{
        Engine, EngineStateSummary, GentleEngine, Operation, OperationProgress, ProjectState,
        RenderSvgMode, TfbsProgress, Workflow,
    },
    engine_shell::{
        execute_shell_command_with_options, parse_shell_line, parse_shell_tokens, shell_help_text,
        ShellExecutionOptions,
    },
    genomes::{
        GenomeGeneRecord, PrepareGenomeProgress, DEFAULT_GENOME_CATALOG_PATH,
        DEFAULT_HELPER_GENOME_CATALOG_PATH,
    },
};
use regex::{Regex, RegexBuilder};
use serde::{Deserialize, Serialize};
use serde_json::json;
use std::{
    collections::{BTreeSet, HashMap},
    env, fs,
    path::Path,
};

const DEFAULT_STATE_PATH: &str = ".gentle_state.json";
const DEFAULT_REBASE_RESOURCE_PATH: &str = "data/resources/rebase.enzymes.json";
const DEFAULT_JASPAR_RESOURCE_PATH: &str = "data/resources/jaspar.motifs.json";

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

#[derive(Debug, Clone, Serialize, Deserialize)]
struct RebaseEnzymeRecord {
    id: usize,
    name: String,
    sequence: String,
    cut: isize,
    overlap: isize,
    #[serde(rename = "type")]
    enzyme_type: String,
    #[serde(skip_serializing_if = "Option::is_none")]
    note: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    prototype: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    suppliers: Option<Vec<String>>,
    #[serde(skip_serializing_if = "Option::is_none")]
    commercial: Option<bool>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
struct JasparMotifRecord {
    id: String,
    #[serde(skip_serializing_if = "Option::is_none")]
    name: Option<String>,
    consensus_iupac: String,
    length: usize,
    pfm: JasparPfmRows,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
struct JasparPfmRows {
    a: Vec<f64>,
    c: Vec<f64>,
    g: Vec<f64>,
    t: Vec<f64>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
struct JasparMotifSnapshot {
    schema: String,
    source: String,
    fetched_at_unix_ms: u128,
    motif_count: usize,
    motifs: Vec<JasparMotifRecord>,
}

fn now_unix_ms() -> u128 {
    std::time::SystemTime::now()
        .duration_since(std::time::UNIX_EPOCH)
        .map(|d| d.as_millis())
        .unwrap_or(0)
}

fn is_iupac_base(c: char) -> bool {
    matches!(
        c.to_ascii_uppercase(),
        'A' | 'C'
            | 'G'
            | 'T'
            | 'U'
            | 'W'
            | 'S'
            | 'M'
            | 'K'
            | 'R'
            | 'Y'
            | 'B'
            | 'D'
            | 'H'
            | 'V'
            | 'N'
    )
}

fn parse_int_pair_in_parens(spec: &str) -> Option<(isize, isize)> {
    let start = spec.find('(')?;
    let end = spec[start + 1..].find(')')? + start + 1;
    let inside = &spec[start + 1..end];
    let slash = inside.find('/')?;
    let left = inside[..slash].trim().parse::<isize>().ok()?;
    let right = inside[slash + 1..].trim().parse::<isize>().ok()?;
    Some((left, right))
}

fn parse_rebase_site(spec: &str) -> Option<(String, isize, isize)> {
    let sequence: String = spec
        .chars()
        .filter(|c| is_iupac_base(*c))
        .map(|c| c.to_ascii_uppercase())
        .collect();
    if sequence.is_empty() {
        return None;
    }

    if let Some((top, bottom)) = parse_int_pair_in_parens(spec) {
        let cut = if top == 0 { 1 } else { top };
        return Some((sequence, cut, bottom - top));
    }

    let mut pos = 0isize;
    let mut top_cut: Option<isize> = None;
    let mut bottom_cut: Option<isize> = None;
    for ch in spec.chars() {
        if is_iupac_base(ch) {
            pos += 1;
            continue;
        }
        if ch == '^' {
            top_cut = Some(pos);
        } else if ch == '_' {
            bottom_cut = Some(pos);
        }
    }

    let cut = top_cut.unwrap_or(1);
    let overlap = match (top_cut, bottom_cut) {
        (Some(t), Some(b)) => b - t,
        _ => 0,
    };
    Some((sequence, cut, overlap))
}

fn parse_suppliers(field: Option<&String>) -> Vec<String> {
    let Some(value) = field else {
        return vec![];
    };
    value
        .chars()
        .filter(|c| c.is_ascii_alphanumeric())
        .map(|c| c.to_string())
        .collect()
}

fn parse_rebase_withrefm(text: &str, commercial_only: bool) -> Vec<RebaseEnzymeRecord> {
    let mut fields: HashMap<String, String> = HashMap::new();
    let mut last_tag: Option<String> = None;
    let mut out: Vec<RebaseEnzymeRecord> = vec![];

    let push_entry = |fields: &HashMap<String, String>, out: &mut Vec<RebaseEnzymeRecord>| {
        let Some(name) = fields.get("1").cloned() else {
            return;
        };
        let Some(spec) = fields.get("3").cloned() else {
            return;
        };
        let Some((sequence, cut, overlap)) = parse_rebase_site(&spec) else {
            return;
        };
        let suppliers = parse_suppliers(fields.get("7"));
        let commercial = !suppliers.is_empty();
        if commercial_only && !commercial {
            return;
        }
        let id = out.len() + 1;
        out.push(RebaseEnzymeRecord {
            id,
            name: name.trim().to_string(),
            sequence,
            cut,
            overlap,
            enzyme_type: "restriction".to_string(),
            note: fields.get("5").cloned().filter(|s| !s.trim().is_empty()),
            prototype: fields.get("2").cloned().filter(|s| !s.trim().is_empty()),
            suppliers: if suppliers.is_empty() {
                None
            } else {
                Some(suppliers)
            },
            commercial: Some(commercial),
        });
    };

    for raw in text.lines() {
        let line = raw.trim_end();
        if line.trim() == "//" {
            push_entry(&fields, &mut out);
            fields.clear();
            last_tag = None;
            continue;
        }
        if let Some(rest) = line.strip_prefix('<') {
            if let Some(pos) = rest.find('>') {
                let tag = rest[..pos].trim().to_string();
                let value = rest[pos + 1..].trim().to_string();
                fields.insert(tag.clone(), value);
                last_tag = Some(tag);
                continue;
            }
        }
        if let Some(tag) = &last_tag {
            let extra = line.trim();
            if !extra.is_empty() {
                let entry = fields.entry(tag.clone()).or_default();
                if !entry.is_empty() {
                    entry.push(' ');
                }
                entry.push_str(extra);
            }
        }
    }
    if !fields.is_empty() {
        push_entry(&fields, &mut out);
    }
    out
}

fn parse_jaspar_row(line: &str) -> Result<(char, Vec<f64>), String> {
    let trimmed = line.trim();
    if trimmed.is_empty() {
        return Err("Empty JASPAR row".to_string());
    }
    let base = trimmed
        .chars()
        .next()
        .ok_or_else(|| "Missing JASPAR base row label".to_string())?
        .to_ascii_uppercase();
    if !matches!(base, 'A' | 'C' | 'G' | 'T') {
        return Err(format!("Unsupported JASPAR row label '{base}'"));
    }
    let payload = trimmed[1..].replace(['[', ']', ','], " ");
    let values = payload
        .split_whitespace()
        .map(|v| {
            v.parse::<f64>()
                .map_err(|e| format!("Invalid JASPAR matrix number '{v}': {e}"))
        })
        .collect::<Result<Vec<_>, _>>()?;
    if values.is_empty() {
        return Err("JASPAR row has no values".to_string());
    }
    Ok((base, values))
}

fn bases_to_iupac(mut bases: Vec<char>) -> char {
    bases.sort_unstable();
    bases.dedup();
    match bases.as_slice() {
        ['A'] => 'A',
        ['C'] => 'C',
        ['G'] => 'G',
        ['T'] => 'T',
        ['A', 'C'] => 'M',
        ['A', 'G'] => 'R',
        ['A', 'T'] => 'W',
        ['C', 'G'] => 'S',
        ['C', 'T'] => 'Y',
        ['G', 'T'] => 'K',
        ['A', 'C', 'G'] => 'V',
        ['A', 'C', 'T'] => 'H',
        ['A', 'G', 'T'] => 'D',
        ['C', 'G', 'T'] => 'B',
        _ => 'N',
    }
}

fn matrix_consensus_iupac(a: &[f64], c: &[f64], g: &[f64], t: &[f64]) -> Result<String, String> {
    if !(a.len() == c.len() && c.len() == g.len() && g.len() == t.len()) {
        return Err("JASPAR rows do not have equal length".to_string());
    }
    let mut out = String::with_capacity(a.len());
    for i in 0..a.len() {
        let max_val = a[i].max(c[i]).max(g[i]).max(t[i]);
        if max_val <= 0.0 {
            out.push('N');
            continue;
        }
        let mut winners = Vec::new();
        if (a[i] - max_val).abs() < f64::EPSILON {
            winners.push('A');
        }
        if (c[i] - max_val).abs() < f64::EPSILON {
            winners.push('C');
        }
        if (g[i] - max_val).abs() < f64::EPSILON {
            winners.push('G');
        }
        if (t[i] - max_val).abs() < f64::EPSILON {
            winners.push('T');
        }
        out.push(bases_to_iupac(winners));
    }
    Ok(out)
}

fn parse_jaspar_motifs(text: &str) -> Result<Vec<JasparMotifRecord>, String> {
    let lines = text.lines().collect::<Vec<_>>();
    let mut i = 0usize;
    let mut motifs = Vec::new();

    while i < lines.len() {
        let line = lines[i].trim();
        i += 1;
        if line.is_empty() {
            continue;
        }
        if !line.starts_with('>') {
            continue;
        }
        let header = line.trim_start_matches('>').trim();
        let mut parts = header.splitn(2, char::is_whitespace);
        let id = parts
            .next()
            .ok_or_else(|| "Malformed JASPAR header".to_string())?
            .trim()
            .to_string();
        let name = parts
            .next()
            .map(|s| s.trim().to_string())
            .filter(|s| !s.is_empty());

        let mut rows: HashMap<char, Vec<f64>> = HashMap::new();
        while i < lines.len() {
            let row_line = lines[i].trim();
            if row_line.is_empty() {
                i += 1;
                continue;
            }
            if row_line.starts_with('>') {
                break;
            }
            let (base, vals) = parse_jaspar_row(row_line)?;
            rows.insert(base, vals);
            i += 1;
            if rows.len() == 4 {
                break;
            }
        }

        let (Some(a), Some(c), Some(g), Some(t)) = (
            rows.get(&'A'),
            rows.get(&'C'),
            rows.get(&'G'),
            rows.get(&'T'),
        ) else {
            return Err(format!(
                "JASPAR motif '{id}' is missing one or more A/C/G/T rows"
            ));
        };
        let consensus = matrix_consensus_iupac(a, c, g, t)?;
        motifs.push(JasparMotifRecord {
            id,
            name,
            length: consensus.len(),
            consensus_iupac: consensus,
            pfm: JasparPfmRows {
                a: a.clone(),
                c: c.clone(),
                g: g.clone(),
                t: t.clone(),
            },
        });
    }

    Ok(motifs)
}

fn ensure_parent_dir(path: &str) -> Result<(), String> {
    let parent = Path::new(path)
        .parent()
        .map(|p| p.to_path_buf())
        .unwrap_or_else(|| Path::new(".").to_path_buf());
    fs::create_dir_all(parent)
        .map_err(|e| format!("Could not create output directory for '{path}': {e}"))
}

fn usage() {
    eprintln!(
        "Usage:\n  \
  gentle_cli --help\n  \
  gentle_cli --version\n  \
  gentle_cli [--state PATH|--project PATH] [--progress|--progress-stderr|--progress-stdout] [--allow-screenshots] COMMAND ...\n\n  \
  gentle_cli [--state PATH|--project PATH] capabilities\n  \
  gentle_cli [--state PATH|--project PATH] op '<operation-json>'\n  \
  gentle_cli [--state PATH|--project PATH] workflow '<workflow-json>'\n  \
  gentle_cli [--state PATH|--project PATH] state-summary\n  \
  gentle_cli [--state PATH|--project PATH] export-state PATH\n  \
  gentle_cli [--state PATH|--project PATH] import-state PATH\n  \
  gentle_cli [--state PATH|--project PATH] save-project PATH\n  \
  gentle_cli [--state PATH|--project PATH] load-project PATH\n  \
  gentle_cli [--state PATH|--project PATH] render-svg SEQ_ID linear|circular OUTPUT.svg\n  \
  gentle_cli [--state PATH|--project PATH] render-rna-svg SEQ_ID OUTPUT.svg\n  \
  gentle_cli [--state PATH|--project PATH] rna-info SEQ_ID\n  \
  gentle_cli [--state PATH|--project PATH] render-lineage-svg OUTPUT.svg\n\n  \
  gentle_cli [--state PATH|--project PATH] shell 'state-summary'\n  \
  gentle_cli [--state PATH|--project PATH] shell 'op <operation-json>'\n\n  \
  gentle_cli [--state PATH|--project PATH] render-pool-gel-svg IDS OUTPUT.svg [--ladders NAME[,NAME]]\n\n  \
  gentle_cli [--state PATH|--project PATH] [--allow-screenshots] screenshot-window OUTPUT.png\n\n  \
  gentle_cli [--state PATH|--project PATH] ladders list [--molecule dna|rna] [--filter TEXT]\n  \
  gentle_cli [--state PATH|--project PATH] ladders export OUTPUT.json [--molecule dna|rna] [--filter TEXT]\n\n  \
  gentle_cli [--state PATH|--project PATH] export-pool IDS OUTPUT.pool.gentle.json [HUMAN_ID]\n  \
  gentle_cli [--state PATH|--project PATH] import-pool INPUT.pool.gentle.json [PREFIX]\n\n  \
  gentle_cli genomes list [--catalog PATH]\n  \
  gentle_cli genomes validate-catalog [--catalog PATH]\n  \
  gentle_cli genomes status GENOME_ID [--catalog PATH] [--cache-dir PATH]\n  \
  gentle_cli genomes genes GENOME_ID [--catalog PATH] [--cache-dir PATH] [--filter TEXT] [--limit N] [--offset N]\n  \
  gentle_cli [--state PATH|--project PATH] genomes prepare GENOME_ID [--catalog PATH] [--cache-dir PATH]\n  \
  gentle_cli [--state PATH|--project PATH] genomes extract-region GENOME_ID CHR START END [--output-id ID] [--catalog PATH] [--cache-dir PATH]\n  \
  gentle_cli [--state PATH|--project PATH] genomes extract-gene GENOME_ID QUERY [--occurrence N] [--output-id ID] [--catalog PATH] [--cache-dir PATH]\n\n  \
  gentle_cli helpers list [--catalog PATH]\n  \
  gentle_cli helpers validate-catalog [--catalog PATH]\n  \
  gentle_cli helpers status HELPER_ID [--catalog PATH] [--cache-dir PATH]\n  \
  gentle_cli helpers genes HELPER_ID [--catalog PATH] [--cache-dir PATH] [--filter TEXT] [--limit N] [--offset N]\n  \
  gentle_cli [--state PATH|--project PATH] helpers prepare HELPER_ID [--catalog PATH] [--cache-dir PATH]\n  \
  gentle_cli [--state PATH|--project PATH] helpers extract-region HELPER_ID CHR START END [--output-id ID] [--catalog PATH] [--cache-dir PATH]\n  \
  gentle_cli [--state PATH|--project PATH] helpers extract-gene HELPER_ID QUERY [--occurrence N] [--output-id ID] [--catalog PATH] [--cache-dir PATH]\n\n  \
  gentle_cli [--state PATH|--project PATH] tracks import-bed SEQ_ID PATH [--name NAME] [--min-score N] [--max-score N] [--clear-existing]\n\n  \
  gentle_cli resources sync-rebase INPUT.withrefm [OUTPUT.rebase.json] [--commercial-only]\n  \
  gentle_cli resources sync-jaspar INPUT.jaspar.txt [OUTPUT.motifs.json]\n\n  \
  Tip: pass @file.json instead of inline JSON\n  \
  --project is an alias of --state for project.gentle.json files\n\n  \
  Shell help:\n  \
  {shell_help}"
        ,
        shell_help = shell_help_text()
    );
}

fn load_json_arg(value: &str) -> Result<String, String> {
    if let Some(path) = value.strip_prefix('@') {
        fs::read_to_string(path).map_err(|e| format!("Could not read JSON file '{path}': {e}"))
    } else {
        Ok(value.to_string())
    }
}

fn read_text_input(path_or_url: &str) -> Result<String, String> {
    if path_or_url.starts_with("http://") || path_or_url.starts_with("https://") {
        let response = std::panic::catch_unwind(|| reqwest::blocking::get(path_or_url))
            .map_err(|_| {
                format!("Could not fetch URL '{path_or_url}': networking backend panicked")
            })?
            .map_err(|e| format!("Could not fetch URL '{path_or_url}': {e}"))?;
        if !response.status().is_success() {
            return Err(format!(
                "Could not fetch URL '{path_or_url}': HTTP {}",
                response.status()
            ));
        }
        response
            .text()
            .map_err(|e| format!("Could not read URL response '{path_or_url}': {e}"))
    } else {
        fs::read_to_string(path_or_url)
            .map_err(|e| format!("Could not read file '{path_or_url}': {e}"))
    }
}

fn load_state(path: &str) -> Result<ProjectState, String> {
    if std::path::Path::new(path).exists() {
        ProjectState::load_from_path(path).map_err(|e| e.to_string())
    } else {
        Ok(ProjectState::default())
    }
}

fn print_json<T: Serialize>(value: &T) -> Result<(), String> {
    let text = serde_json::to_string_pretty(value)
        .map_err(|e| format!("Could not serialize JSON output: {e}"))?;
    println!("{text}");
    Ok(())
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum ProgressSink {
    Stdout,
    Stderr,
}

#[derive(Clone, Debug)]
struct GlobalCliArgs {
    state_path: String,
    cmd_idx: usize,
    progress_sink: Option<ProgressSink>,
    allow_screenshots: bool,
}

fn parse_global_args(args: &[String]) -> Result<GlobalCliArgs, String> {
    let mut state_path = DEFAULT_STATE_PATH.to_string();
    let mut progress_sink: Option<ProgressSink> = None;
    let mut allow_screenshots = false;
    let mut idx = 1usize;

    while idx < args.len() {
        match args[idx].as_str() {
            "--state" | "--project" => {
                if idx + 1 >= args.len() {
                    return Err(format!("Missing PATH after {}", args[idx]));
                }
                state_path = args[idx + 1].clone();
                idx += 2;
            }
            "--progress" | "--progress-stderr" => {
                progress_sink = Some(ProgressSink::Stderr);
                idx += 1;
            }
            "--progress-stdout" => {
                progress_sink = Some(ProgressSink::Stdout);
                idx += 1;
            }
            "--allow-screenshots" => {
                allow_screenshots = true;
                idx += 1;
            }
            _ => break,
        }
    }

    Ok(GlobalCliArgs {
        state_path,
        cmd_idx: idx,
        progress_sink,
        allow_screenshots,
    })
}

struct ProgressPrinter {
    sink: ProgressSink,
    last_total_tenths: Option<i64>,
    last_motif_index: Option<usize>,
    last_genome_phase: Option<String>,
    last_genome_percent_tenths: Option<i64>,
    last_genome_bytes_bucket: Option<u64>,
}

impl ProgressPrinter {
    fn new(sink: ProgressSink) -> Self {
        Self {
            sink,
            last_total_tenths: None,
            last_motif_index: None,
            last_genome_phase: None,
            last_genome_percent_tenths: None,
            last_genome_bytes_bucket: None,
        }
    }

    fn print_line(&self, line: &str) {
        match self.sink {
            ProgressSink::Stdout => println!("{line}"),
            ProgressSink::Stderr => eprintln!("{line}"),
        }
    }

    fn on_tfbs_progress(&mut self, p: TfbsProgress) {
        let total_tenths = (p.total_percent * 10.0).floor() as i64;
        let motif_changed = self.last_motif_index != Some(p.motif_index);
        let progressed = self
            .last_total_tenths
            .map(|prev| total_tenths > prev)
            .unwrap_or(true);
        let done = (p.total_percent - 100.0).abs() < f64::EPSILON;
        if motif_changed || progressed || done {
            self.last_motif_index = Some(p.motif_index);
            self.last_total_tenths = Some(total_tenths);
            self.print_line(&format!(
                "progress tfbs seq={} motif={}/{} id={} motif_pct={:.1} total_pct={:.1} steps={}/{}",
                p.seq_id,
                p.motif_index,
                p.motif_count,
                p.motif_id,
                p.motif_percent,
                p.total_percent,
                p.scanned_steps,
                p.total_steps
            ));
        }
    }

    fn on_genome_prepare_progress(&mut self, p: PrepareGenomeProgress) {
        let phase_changed = self
            .last_genome_phase
            .as_ref()
            .map(|prev| prev != &p.phase)
            .unwrap_or(true);
        let percent_tenths = p.percent.map(|v| (v * 10.0).floor() as i64);
        let percent_progressed = percent_tenths
            .map(|value| {
                self.last_genome_percent_tenths
                    .map(|prev| value > prev)
                    .unwrap_or(true)
            })
            .unwrap_or(false);
        let bytes_bucket = p.bytes_done / (8 * 1024 * 1024);
        let bytes_progressed = self
            .last_genome_bytes_bucket
            .map(|prev| bytes_bucket > prev)
            .unwrap_or(true);
        let done = p.phase == "ready"
            || p.percent
                .map(|v| (v - 100.0).abs() < f64::EPSILON)
                .unwrap_or(false);
        if phase_changed || percent_progressed || bytes_progressed || done {
            self.last_genome_phase = Some(p.phase.clone());
            if let Some(v) = percent_tenths {
                self.last_genome_percent_tenths = Some(v);
            }
            self.last_genome_bytes_bucket = Some(bytes_bucket);
            let total = p
                .bytes_total
                .map(|v| v.to_string())
                .unwrap_or_else(|| "?".to_string());
            let pct = p
                .percent
                .map(|v| format!("{v:.1}"))
                .unwrap_or_else(|| "?".to_string());
            self.print_line(&format!(
                "progress genome id={} phase={} pct={} bytes={}/{} item={}",
                p.genome_id, p.phase, pct, p.bytes_done, total, p.item
            ));
        }
    }

    fn on_progress(&mut self, progress: OperationProgress) {
        match progress {
            OperationProgress::Tfbs(p) => self.on_tfbs_progress(p),
            OperationProgress::GenomePrepare(p) => self.on_genome_prepare_progress(p),
        }
    }
}

fn summarize_state(engine: &GentleEngine) -> EngineStateSummary {
    engine.summarize_state()
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

fn unique_id(
    existing: &std::collections::HashMap<String, gentle::dna_sequence::DNAsequence>,
    base: &str,
) -> String {
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

fn apply_member_overhang(
    member: &PoolMember,
    dna: &mut gentle::dna_sequence::DNAsequence,
) -> Result<(), String> {
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
    let patched: gentle::dna_sequence::DNAsequence =
        serde_json::from_value(value).map_err(|e| format!("Could not restore overhang: {e}"))?;
    *dna = patched;
    Ok(())
}

fn main() {
    if let Err(e) = run() {
        eprintln!("{e}");
        std::process::exit(1);
    }
}

fn run() -> Result<(), String> {
    let args: Vec<String> = env::args().collect();
    if args.len() <= 1 {
        usage();
        return Err("Missing command".to_string());
    }
    if args.iter().any(|a| a == "--help" || a == "-h") {
        usage();
        return Ok(());
    }
    if args.iter().any(|a| a == "--version" || a == "-V") {
        println!("{}", about::version_cli_text());
        return Ok(());
    }

    let global = parse_global_args(&args)?;
    let state_path = global.state_path.clone();
    let cmd_idx = global.cmd_idx;
    let shell_options = ShellExecutionOptions {
        allow_screenshots: global.allow_screenshots,
    };
    if args.len() <= cmd_idx {
        usage();
        return Err("Missing command".to_string());
    }

    let command = &args[cmd_idx];

    if matches!(
        command.as_str(),
        "genomes"
            | "helpers"
            | "resources"
            | "import-pool"
            | "ladders"
            | "tracks"
            | "screenshot-window"
    ) {
        let tokens = args[cmd_idx..].to_vec();
        let shell_command = parse_shell_tokens(&tokens)?;
        let mut engine = GentleEngine::from_state(load_state(&state_path)?);
        let run = execute_shell_command_with_options(&mut engine, &shell_command, &shell_options)?;
        if run.state_changed {
            engine
                .state()
                .save_to_path(&state_path)
                .map_err(|e| e.to_string())?;
        }
        return print_json(&run.output);
    }

    match command.as_str() {
        "genomes" | "helpers" => {
            let helper_mode = command == "helpers";
            let default_catalog = if helper_mode {
                DEFAULT_HELPER_GENOME_CATALOG_PATH
            } else {
                DEFAULT_GENOME_CATALOG_PATH
            };
            let label = if helper_mode { "helpers" } else { "genomes" };
            if args.len() <= cmd_idx + 1 {
                usage();
                return Err(format!("{label} requires a subcommand"));
            }
            match args[cmd_idx + 1].as_str() {
                "list" => {
                    let mut catalog_path: Option<String> = None;
                    let mut idx = cmd_idx + 2;
                    while idx < args.len() {
                        match args[idx].as_str() {
                            "--catalog" => {
                                if idx + 1 >= args.len() {
                                    return Err(format!(
                                        "Missing PATH after --catalog for {label} list"
                                    ));
                                }
                                catalog_path = Some(args[idx + 1].clone());
                                idx += 2;
                            }
                            other => {
                                return Err(format!("Unknown option '{}' for {label} list", other))
                            }
                        }
                    }
                    let resolved_catalog = catalog_path
                        .as_deref()
                        .map(str::trim)
                        .filter(|v| !v.is_empty())
                        .or_else(|| helper_mode.then_some(default_catalog));
                    let genomes = GentleEngine::list_reference_genomes(resolved_catalog)
                        .map_err(|e| e.to_string())?;
                    let effective_catalog = catalog_path
                        .clone()
                        .filter(|v| !v.trim().is_empty())
                        .unwrap_or_else(|| default_catalog.to_string());
                    print_json(&json!({
                        "catalog_path": effective_catalog,
                        "genome_count": genomes.len(),
                        "genomes": genomes,
                    }))
                }
                "status" => {
                    if args.len() <= cmd_idx + 2 {
                        usage();
                        return Err(format!(
                            "{label} status requires GENOME_ID [--catalog PATH] [--cache-dir PATH]"
                        ));
                    }
                    let genome_id = args[cmd_idx + 2].clone();
                    let mut catalog_path: Option<String> = None;
                    let mut cache_dir: Option<String> = None;
                    let mut idx = cmd_idx + 3;
                    while idx < args.len() {
                        match args[idx].as_str() {
                            "--catalog" => {
                                if idx + 1 >= args.len() {
                                    return Err(format!(
                                        "Missing PATH after --catalog for {label} status"
                                    ));
                                }
                                catalog_path = Some(args[idx + 1].clone());
                                idx += 2;
                            }
                            "--cache-dir" => {
                                if idx + 1 >= args.len() {
                                    return Err(format!(
                                        "Missing PATH after --cache-dir for {label} status"
                                    ));
                                }
                                cache_dir = Some(args[idx + 1].clone());
                                idx += 2;
                            }
                            other => {
                                return Err(format!(
                                    "Unknown option '{}' for {label} status",
                                    other
                                ))
                            }
                        }
                    }
                    let resolved_catalog = catalog_path
                        .as_deref()
                        .map(str::trim)
                        .filter(|v| !v.is_empty())
                        .or_else(|| helper_mode.then_some(default_catalog));
                    let prepared = GentleEngine::is_reference_genome_prepared(
                        resolved_catalog,
                        &genome_id,
                        cache_dir.as_deref(),
                    )
                    .map_err(|e| e.to_string())?;
                    let source_plan = GentleEngine::describe_reference_genome_sources(
                        resolved_catalog,
                        &genome_id,
                        cache_dir.as_deref(),
                    )
                    .map_err(|e| e.to_string())?;
                    let effective_catalog = catalog_path
                        .clone()
                        .filter(|v| !v.trim().is_empty())
                        .unwrap_or_else(|| default_catalog.to_string());
                    print_json(&json!({
                        "genome_id": genome_id,
                        "catalog_path": effective_catalog,
                        "cache_dir": cache_dir,
                        "prepared": prepared,
                        "sequence_source_type": source_plan.sequence_source_type,
                        "annotation_source_type": source_plan.annotation_source_type,
                        "sequence_source": source_plan.sequence_source,
                        "annotation_source": source_plan.annotation_source,
                    }))
                }
                "genes" => {
                    if args.len() <= cmd_idx + 2 {
                        usage();
                        return Err(format!(
                            "{label} genes requires GENOME_ID [--catalog PATH] [--cache-dir PATH] [--filter REGEX] [--biotype NAME] [--limit N] [--offset N]"
                        ));
                    }
                    let genome_id = args[cmd_idx + 2].clone();
                    let mut catalog_path: Option<String> = None;
                    let mut cache_dir: Option<String> = None;
                    let mut filter = String::new();
                    let mut biotype_filters: Vec<String> = vec![];
                    let mut limit: usize = 200;
                    let mut offset: usize = 0;
                    let mut idx = cmd_idx + 3;
                    while idx < args.len() {
                        match args[idx].as_str() {
                            "--catalog" => {
                                if idx + 1 >= args.len() {
                                    return Err(format!(
                                        "Missing PATH after --catalog for {label} genes"
                                    ));
                                }
                                catalog_path = Some(args[idx + 1].clone());
                                idx += 2;
                            }
                            "--cache-dir" => {
                                if idx + 1 >= args.len() {
                                    return Err(format!(
                                        "Missing PATH after --cache-dir for {label} genes"
                                    ));
                                }
                                cache_dir = Some(args[idx + 1].clone());
                                idx += 2;
                            }
                            "--filter" => {
                                if idx + 1 >= args.len() {
                                    return Err(format!(
                                        "Missing TEXT after --filter for {label} genes"
                                    ));
                                }
                                filter = args[idx + 1].clone();
                                idx += 2;
                            }
                            "--biotype" => {
                                if idx + 1 >= args.len() {
                                    return Err(format!(
                                        "Missing NAME after --biotype for {label} genes"
                                    ));
                                }
                                biotype_filters.push(args[idx + 1].clone());
                                idx += 2;
                            }
                            "--limit" => {
                                if idx + 1 >= args.len() {
                                    return Err(format!("Missing N after --limit for {label} genes"));
                                }
                                limit = args[idx + 1].parse::<usize>().map_err(|e| {
                                    format!("Invalid --limit value '{}': {}", args[idx + 1], e)
                                })?;
                                if limit == 0 {
                                    return Err("--limit must be >= 1".to_string());
                                }
                                idx += 2;
                            }
                            "--offset" => {
                                if idx + 1 >= args.len() {
                                    return Err(format!(
                                        "Missing N after --offset for {label} genes"
                                    ));
                                }
                                offset = args[idx + 1].parse::<usize>().map_err(|e| {
                                    format!("Invalid --offset value '{}': {}", args[idx + 1], e)
                                })?;
                                idx += 2;
                            }
                            other => {
                                return Err(format!(
                                    "Unknown option '{}' for {label} genes",
                                    other
                                ))
                            }
                        }
                    }

                    let resolved_catalog = catalog_path
                        .as_deref()
                        .map(str::trim)
                        .filter(|v| !v.is_empty())
                        .or_else(|| helper_mode.then_some(default_catalog));
                    let genes = GentleEngine::list_reference_genome_genes(
                        resolved_catalog,
                        &genome_id,
                        cache_dir.as_deref(),
                    )
                    .map_err(|e| e.to_string())?;
                    let filter_regex = compile_gene_filter_regex(&filter)?;
                    let available_biotypes = collect_biotypes(&genes);
                    let allowed_biotypes_lower: Vec<String> = biotype_filters
                        .iter()
                        .map(|v| v.trim().to_ascii_lowercase())
                        .filter(|v| !v.is_empty())
                        .collect();
                    let filtered: Vec<GenomeGeneRecord> = genes
                        .into_iter()
                        .filter(|g| {
                            genome_gene_matches_filter(
                                g,
                                filter_regex.as_ref(),
                                &allowed_biotypes_lower,
                            )
                        })
                        .collect();
                    let total = filtered.len();
                    let offset = offset.min(total);
                    let returned: Vec<GenomeGeneRecord> =
                        filtered.into_iter().skip(offset).take(limit).collect();
                    let effective_catalog = catalog_path
                        .clone()
                        .filter(|v| !v.trim().is_empty())
                        .unwrap_or_else(|| default_catalog.to_string());
                    print_json(&json!({
                        "genome_id": genome_id,
                        "catalog_path": effective_catalog,
                        "cache_dir": cache_dir,
                        "filter": filter,
                        "biotype_filter": biotype_filters,
                        "available_biotypes": available_biotypes,
                        "offset": offset,
                        "limit": limit,
                        "total_matches": total,
                        "returned": returned.len(),
                        "genes": returned,
                    }))
                }
                "prepare" => {
                    if args.len() <= cmd_idx + 2 {
                        usage();
                        return Err(format!(
                            "{label} prepare requires GENOME_ID [--catalog PATH] [--cache-dir PATH]"
                        ));
                    }
                    let genome_id = args[cmd_idx + 2].clone();
                    let mut catalog_path: Option<String> = None;
                    let mut cache_dir: Option<String> = None;
                    let mut idx = cmd_idx + 3;
                    while idx < args.len() {
                        match args[idx].as_str() {
                            "--catalog" => {
                                if idx + 1 >= args.len() {
                                    return Err(format!(
                                        "Missing PATH after --catalog for {label} prepare"
                                    ));
                                }
                                catalog_path = Some(args[idx + 1].clone());
                                idx += 2;
                            }
                            "--cache-dir" => {
                                if idx + 1 >= args.len() {
                                    return Err(format!(
                                        "Missing PATH after --cache-dir for {label} prepare"
                                    ));
                                }
                                cache_dir = Some(args[idx + 1].clone());
                                idx += 2;
                            }
                            other => {
                                return Err(format!(
                                    "Unknown option '{}' for {label} prepare",
                                    other
                                ))
                            }
                        }
                    }
                    let op_catalog_path = catalog_path
                        .filter(|v| !v.trim().is_empty())
                        .or_else(|| helper_mode.then_some(default_catalog.to_string()));
                    let mut engine = GentleEngine::from_state(load_state(&state_path)?);
                    let op = Operation::PrepareGenome {
                        genome_id,
                        catalog_path: op_catalog_path,
                        cache_dir,
                    };
                    let result = if let Some(sink) = global.progress_sink {
                        let mut printer = ProgressPrinter::new(sink);
                        engine
                            .apply_with_progress(op, |p| printer.on_progress(p))
                            .map_err(|e| e.to_string())?
                    } else {
                        engine.apply(op).map_err(|e| e.to_string())?
                    };
                    engine
                        .state()
                        .save_to_path(&state_path)
                        .map_err(|e| e.to_string())?;
                    print_json(&result)
                }
                "extract-region" => {
                    if args.len() <= cmd_idx + 6 {
                        usage();
                        return Err(format!(
                            "{label} extract-region requires GENOME_ID CHR START END [--output-id ID] [--catalog PATH] [--cache-dir PATH]"
                        ));
                    }
                    let genome_id = args[cmd_idx + 2].clone();
                    let chromosome = args[cmd_idx + 3].clone();
                    let start_1based = args[cmd_idx + 4].parse::<usize>().map_err(|e| {
                        format!("Invalid START coordinate '{}': {}", args[cmd_idx + 4], e)
                    })?;
                    let end_1based = args[cmd_idx + 5]
                        .parse::<usize>()
                        .map_err(|e| format!("Invalid END coordinate '{}': {}", args[cmd_idx + 5], e))?;
                    let mut output_id: Option<String> = None;
                    let mut catalog_path: Option<String> = None;
                    let mut cache_dir: Option<String> = None;
                    let mut idx = cmd_idx + 6;
                    while idx < args.len() {
                        match args[idx].as_str() {
                            "--output-id" => {
                                if idx + 1 >= args.len() {
                                    return Err(format!(
                                        "Missing ID after --output-id for {label} extract-region"
                                    ));
                                }
                                output_id = Some(args[idx + 1].clone());
                                idx += 2;
                            }
                            "--catalog" => {
                                if idx + 1 >= args.len() {
                                    return Err(format!(
                                        "Missing PATH after --catalog for {label} extract-region"
                                    ));
                                }
                                catalog_path = Some(args[idx + 1].clone());
                                idx += 2;
                            }
                            "--cache-dir" => {
                                if idx + 1 >= args.len() {
                                    return Err(format!(
                                        "Missing PATH after --cache-dir for {label} extract-region"
                                    ));
                                }
                                cache_dir = Some(args[idx + 1].clone());
                                idx += 2;
                            }
                            other => {
                                return Err(format!(
                                    "Unknown option '{}' for {label} extract-region",
                                    other
                                ))
                            }
                        }
                    }
                    let op_catalog_path = catalog_path
                        .filter(|v| !v.trim().is_empty())
                        .or_else(|| helper_mode.then_some(default_catalog.to_string()));
                    let mut engine = GentleEngine::from_state(load_state(&state_path)?);
                    let result = engine
                        .apply(Operation::ExtractGenomeRegion {
                            genome_id,
                            chromosome,
                            start_1based,
                            end_1based,
                            output_id,
                            catalog_path: op_catalog_path,
                            cache_dir,
                        })
                        .map_err(|e| e.to_string())?;
                    engine
                        .state()
                        .save_to_path(&state_path)
                        .map_err(|e| e.to_string())?;
                    print_json(&result)
                }
                "extract-gene" => {
                    if args.len() <= cmd_idx + 3 {
                        usage();
                        return Err(format!(
                            "{label} extract-gene requires GENOME_ID QUERY [--occurrence N] [--output-id ID] [--catalog PATH] [--cache-dir PATH]"
                        ));
                    }
                    let genome_id = args[cmd_idx + 2].clone();
                    let gene_query = args[cmd_idx + 3].clone();
                    let mut occurrence: Option<usize> = None;
                    let mut output_id: Option<String> = None;
                    let mut catalog_path: Option<String> = None;
                    let mut cache_dir: Option<String> = None;
                    let mut idx = cmd_idx + 4;
                    while idx < args.len() {
                        match args[idx].as_str() {
                            "--occurrence" => {
                                if idx + 1 >= args.len() {
                                    return Err(format!(
                                        "Missing N after --occurrence for {label} extract-gene"
                                    ));
                                }
                                let occ = args[idx + 1].parse::<usize>().map_err(|e| {
                                    format!(
                                        "Invalid --occurrence value '{}': {}",
                                        args[idx + 1],
                                        e
                                    )
                                })?;
                                if occ == 0 {
                                    return Err("--occurrence must be >= 1".to_string());
                                }
                                occurrence = Some(occ);
                                idx += 2;
                            }
                            "--output-id" => {
                                if idx + 1 >= args.len() {
                                    return Err(format!(
                                        "Missing ID after --output-id for {label} extract-gene"
                                    ));
                                }
                                output_id = Some(args[idx + 1].clone());
                                idx += 2;
                            }
                            "--catalog" => {
                                if idx + 1 >= args.len() {
                                    return Err(format!(
                                        "Missing PATH after --catalog for {label} extract-gene"
                                    ));
                                }
                                catalog_path = Some(args[idx + 1].clone());
                                idx += 2;
                            }
                            "--cache-dir" => {
                                if idx + 1 >= args.len() {
                                    return Err(format!(
                                        "Missing PATH after --cache-dir for {label} extract-gene"
                                    ));
                                }
                                cache_dir = Some(args[idx + 1].clone());
                                idx += 2;
                            }
                            other => {
                                return Err(format!(
                                    "Unknown option '{}' for {label} extract-gene",
                                    other
                                ))
                            }
                        }
                    }
                    let op_catalog_path = catalog_path
                        .filter(|v| !v.trim().is_empty())
                        .or_else(|| helper_mode.then_some(default_catalog.to_string()));
                    let mut engine = GentleEngine::from_state(load_state(&state_path)?);
                    let result = engine
                        .apply(Operation::ExtractGenomeGene {
                            genome_id,
                            gene_query,
                            occurrence,
                            output_id,
                            catalog_path: op_catalog_path,
                            cache_dir,
                        })
                        .map_err(|e| e.to_string())?;
                    engine
                        .state()
                        .save_to_path(&state_path)
                        .map_err(|e| e.to_string())?;
                    print_json(&result)
                }
                other => Err(format!(
                    "Unknown {label} subcommand '{}' (expected list, status, genes, prepare, extract-region, extract-gene)",
                    other
                )),
            }
        }
        "resources" => {
            if args.len() <= cmd_idx + 1 {
                usage();
                return Err(
                    "resources requires a subcommand: sync-rebase or sync-jaspar".to_string(),
                );
            }
            match args[cmd_idx + 1].as_str() {
                "sync-rebase" => {
                    if args.len() <= cmd_idx + 2 {
                        usage();
                        return Err(
                            "resources sync-rebase requires INPUT.withrefm path".to_string()
                        );
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
                    let json = serde_json::to_string_pretty(&enzymes).map_err(|e| {
                        format!("Could not serialize REBASE resource snapshot: {e}")
                    })?;
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
                        return Err(
                            "resources sync-jaspar requires INPUT.jaspar.txt path".to_string()
                        );
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
                    let json = serde_json::to_string_pretty(&snapshot).map_err(|e| {
                        format!("Could not serialize JASPAR resource snapshot: {e}")
                    })?;
                    fs::write(&output, json)
                        .map_err(|e| format!("Could not write JASPAR output '{output}': {e}"))?;
                    println!(
                        "Synced {} JASPAR motifs to '{}'",
                        snapshot.motif_count, output
                    );
                    Ok(())
                }
                other => Err(format!(
                    "Unknown resources subcommand '{other}' (expected sync-rebase or sync-jaspar)"
                )),
            }
        }
        "shell" => {
            if args.len() <= cmd_idx + 1 {
                usage();
                return Err(
                    "shell requires a command string, e.g.: shell 'state-summary'".to_string(),
                );
            }
            let line = args[cmd_idx + 1..].join(" ");
            let command = parse_shell_line(&line)?;
            let mut engine = GentleEngine::from_state(load_state(&state_path)?);
            let run = execute_shell_command_with_options(&mut engine, &command, &shell_options)?;
            if run.state_changed {
                engine
                    .state()
                    .save_to_path(&state_path)
                    .map_err(|e| e.to_string())?;
            }
            print_json(&run.output)
        }
        "capabilities" => {
            print_json(&GentleEngine::capabilities())?;
            Ok(())
        }
        "import-state" | "load-project" => {
            if args.len() <= cmd_idx + 1 {
                usage();
                return Err(format!("Missing path for {command}"));
            }
            let source = &args[cmd_idx + 1];
            let state = ProjectState::load_from_path(source).map_err(|e| e.to_string())?;
            state.save_to_path(&state_path).map_err(|e| e.to_string())?;
            println!("Loaded project from '{source}' into '{state_path}'");
            Ok(())
        }
        "export-state" | "save-project" => {
            if args.len() <= cmd_idx + 1 {
                usage();
                return Err(format!("Missing path for {command}"));
            }
            let target = &args[cmd_idx + 1];
            let state = load_state(&state_path)?;
            state.save_to_path(target).map_err(|e| e.to_string())?;
            println!("Saved project from '{state_path}' to '{target}'");
            Ok(())
        }
        "state-summary" => {
            let state = load_state(&state_path)?;
            let engine = GentleEngine::from_state(state);
            print_json(&summarize_state(&engine))
        }
        "render-svg" => {
            if args.len() <= cmd_idx + 3 {
                usage();
                return Err("render-svg requires: SEQ_ID linear|circular OUTPUT.svg".to_string());
            }
            let seq_id = &args[cmd_idx + 1];
            let mode = &args[cmd_idx + 2];
            let output = &args[cmd_idx + 3];
            let mode = match mode.as_str() {
                "linear" => RenderSvgMode::Linear,
                "circular" => RenderSvgMode::Circular,
                _ => {
                    return Err(format!(
                        "Unknown render mode '{mode}', expected 'linear' or 'circular'"
                    ))
                }
            };
            let mut engine = GentleEngine::from_state(load_state(&state_path)?);
            let result = engine
                .apply(Operation::RenderSequenceSvg {
                    seq_id: seq_id.to_string(),
                    mode,
                    path: output.to_string(),
                })
                .map_err(|e| e.to_string())?;
            engine
                .state()
                .save_to_path(&state_path)
                .map_err(|e| e.to_string())?;
            if let Some(msg) = result.messages.first() {
                println!("{msg}");
            }
            Ok(())
        }
        "render-rna-svg" => {
            if args.len() <= cmd_idx + 2 {
                usage();
                return Err("render-rna-svg requires: SEQ_ID OUTPUT.svg".to_string());
            }
            let seq_id = &args[cmd_idx + 1];
            let output = &args[cmd_idx + 2];
            let mut engine = GentleEngine::from_state(load_state(&state_path)?);
            let result = engine
                .apply(Operation::RenderRnaStructureSvg {
                    seq_id: seq_id.to_string(),
                    path: output.to_string(),
                })
                .map_err(|e| e.to_string())?;
            engine
                .state()
                .save_to_path(&state_path)
                .map_err(|e| e.to_string())?;
            if let Some(msg) = result.messages.first() {
                println!("{msg}");
            }
            Ok(())
        }
        "rna-info" => {
            if args.len() <= cmd_idx + 1 {
                usage();
                return Err("rna-info requires: SEQ_ID".to_string());
            }
            let seq_id = &args[cmd_idx + 1];
            let engine = GentleEngine::from_state(load_state(&state_path)?);
            let report = engine
                .inspect_rna_structure(seq_id)
                .map_err(|e| e.to_string())?;
            print_json(&report)
        }
        "render-lineage-svg" => {
            if args.len() <= cmd_idx + 1 {
                usage();
                return Err("render-lineage-svg requires: OUTPUT.svg".to_string());
            }
            let output = &args[cmd_idx + 1];
            let mut engine = GentleEngine::from_state(load_state(&state_path)?);
            let result = engine
                .apply(Operation::RenderLineageSvg {
                    path: output.to_string(),
                })
                .map_err(|e| e.to_string())?;
            engine
                .state()
                .save_to_path(&state_path)
                .map_err(|e| e.to_string())?;
            if let Some(msg) = result.messages.first() {
                println!("{msg}");
            }
            Ok(())
        }
        "render-pool-gel-svg" => {
            if args.len() <= cmd_idx + 2 {
                usage();
                return Err(
                    "render-pool-gel-svg requires: IDS OUTPUT.svg [--ladders NAME[,NAME]]"
                        .to_string(),
                );
            }
            let ids = args[cmd_idx + 1]
                .split(',')
                .map(|s| s.trim())
                .filter(|s| !s.is_empty())
                .map(|s| s.to_string())
                .collect::<Vec<_>>();
            if ids.is_empty() {
                return Err("render-pool-gel-svg requires at least one sequence id".to_string());
            }
            let output = &args[cmd_idx + 2];
            let mut ladders: Option<Vec<String>> = None;
            let mut idx = cmd_idx + 3;
            while idx < args.len() {
                match args[idx].as_str() {
                    "--ladders" => {
                        if idx + 1 >= args.len() {
                            return Err("Missing value after --ladders".to_string());
                        }
                        let parsed = args[idx + 1]
                            .split(',')
                            .map(|s| s.trim())
                            .filter(|s| !s.is_empty())
                            .map(|s| s.to_string())
                            .collect::<Vec<_>>();
                        if !parsed.is_empty() {
                            ladders = Some(parsed);
                        }
                        idx += 2;
                    }
                    other => {
                        return Err(format!(
                        "Unknown argument '{other}' for render-pool-gel-svg (expected --ladders)"
                    ))
                    }
                }
            }
            let mut engine = GentleEngine::from_state(load_state(&state_path)?);
            let result = engine
                .apply(Operation::RenderPoolGelSvg {
                    inputs: ids,
                    path: output.to_string(),
                    ladders,
                })
                .map_err(|e| e.to_string())?;
            engine
                .state()
                .save_to_path(&state_path)
                .map_err(|e| e.to_string())?;
            if let Some(msg) = result.messages.first() {
                println!("{msg}");
            }
            Ok(())
        }
        "export-pool" => {
            if args.len() <= cmd_idx + 2 {
                usage();
                return Err(
                    "export-pool requires: IDS OUTPUT.pool.gentle.json [HUMAN_ID]".to_string(),
                );
            }
            let ids = args[cmd_idx + 1]
                .split(',')
                .map(|s| s.trim())
                .filter(|s| !s.is_empty())
                .map(|s| s.to_string())
                .collect::<Vec<_>>();
            if ids.is_empty() {
                return Err("export-pool requires at least one sequence id".to_string());
            }
            let output = &args[cmd_idx + 2];
            let human_id = args.get(cmd_idx + 3).cloned();
            let mut engine = GentleEngine::from_state(load_state(&state_path)?);
            let result = engine
                .apply(Operation::ExportPool {
                    inputs: ids,
                    path: output.to_string(),
                    pool_id: Some("pool_export".to_string()),
                    human_id,
                })
                .map_err(|e| e.to_string())?;
            engine
                .state()
                .save_to_path(&state_path)
                .map_err(|e| e.to_string())?;
            print_json(&result)
        }
        "import-pool" => {
            if args.len() <= cmd_idx + 1 {
                usage();
                return Err("import-pool requires: INPUT.pool.gentle.json [PREFIX]".to_string());
            }
            let input = &args[cmd_idx + 1];
            let prefix = args
                .get(cmd_idx + 2)
                .cloned()
                .unwrap_or_else(|| "pool".to_string());
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

            let mut state = load_state(&state_path)?;
            for (idx, member) in pool.members.iter().enumerate() {
                let mut dna = gentle::dna_sequence::DNAsequence::from_sequence(&member.sequence)
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
                state.sequences.insert(id, dna);
            }
            state.save_to_path(&state_path).map_err(|e| e.to_string())?;
            println!(
                "Imported pool '{}' ({} members) into '{}'",
                pool.pool_id, pool.member_count, state_path
            );
            Ok(())
        }
        "op" => {
            if args.len() <= cmd_idx + 1 {
                usage();
                return Err("Missing operation JSON".to_string());
            }
            let json = load_json_arg(&args[cmd_idx + 1])?;
            let op: Operation =
                serde_json::from_str(&json).map_err(|e| format!("Invalid operation JSON: {e}"))?;

            let mut engine = GentleEngine::from_state(load_state(&state_path)?);
            let result = if let Some(sink) = global.progress_sink {
                let mut printer = ProgressPrinter::new(sink);
                engine
                    .apply_with_progress(op, |p| printer.on_progress(p))
                    .map_err(|e| e.to_string())?
            } else {
                engine.apply(op).map_err(|e| e.to_string())?
            };
            engine
                .state()
                .save_to_path(&state_path)
                .map_err(|e| e.to_string())?;
            print_json(&result)
        }
        "workflow" => {
            if args.len() <= cmd_idx + 1 {
                usage();
                return Err("Missing workflow JSON".to_string());
            }
            let json = load_json_arg(&args[cmd_idx + 1])?;
            let workflow: Workflow =
                serde_json::from_str(&json).map_err(|e| format!("Invalid workflow JSON: {e}"))?;

            let mut engine = GentleEngine::from_state(load_state(&state_path)?);
            let results = if let Some(sink) = global.progress_sink {
                let mut printer = ProgressPrinter::new(sink);
                engine
                    .apply_workflow_with_progress(workflow, |p| printer.on_progress(p))
                    .map_err(|e| e.to_string())?
            } else {
                engine.apply_workflow(workflow).map_err(|e| e.to_string())?
            };
            engine
                .state()
                .save_to_path(&state_path)
                .map_err(|e| e.to_string())?;
            print_json(&results)
        }
        _ => {
            usage();
            Err(format!("Unknown command '{command}'"))
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_parse_rebase_site_with_slash_notation() {
        let (seq, cut, overlap) = parse_rebase_site("GAATTC (1/5)").unwrap();
        assert_eq!(seq, "GAATTC");
        assert_eq!(cut, 1);
        assert_eq!(overlap, 4);
    }

    #[test]
    fn test_parse_rebase_withrefm_minimal() {
        let text = r#"
<1>EcoRI
<2>EcoRI
<3>GAATTC (1/5)
<7>N
//
"#;
        let items = parse_rebase_withrefm(text, true);
        assert_eq!(items.len(), 1);
        assert_eq!(items[0].name, "EcoRI");
        assert_eq!(items[0].sequence, "GAATTC");
        assert_eq!(items[0].overlap, 4);
    }

    #[test]
    fn test_parse_jaspar_motifs_consensus() {
        let text = r#"
>MA0001.1 TEST
A [ 10 0 0 0 ]
C [ 0 10 0 0 ]
G [ 0 0 10 0 ]
T [ 0 0 0 10 ]
"#;
        let motifs = parse_jaspar_motifs(text).unwrap();
        assert_eq!(motifs.len(), 1);
        assert_eq!(motifs[0].id, "MA0001.1");
        assert_eq!(motifs[0].consensus_iupac, "ACGT");
    }

    #[test]
    fn test_parse_global_args_progress_stdout_and_state() {
        let args = vec![
            "gentle_cli".to_string(),
            "--state".to_string(),
            "x.json".to_string(),
            "--progress-stdout".to_string(),
            "op".to_string(),
            "{}".to_string(),
        ];
        let parsed = parse_global_args(&args).unwrap();
        assert_eq!(parsed.state_path, "x.json");
        assert_eq!(parsed.progress_sink, Some(ProgressSink::Stdout));
        assert!(!parsed.allow_screenshots);
        assert_eq!(parsed.cmd_idx, 4);
    }

    #[test]
    fn test_parse_global_args_progress_stderr_default_flag() {
        let args = vec![
            "gentle_cli".to_string(),
            "--progress".to_string(),
            "capabilities".to_string(),
        ];
        let parsed = parse_global_args(&args).unwrap();
        assert_eq!(parsed.state_path, DEFAULT_STATE_PATH);
        assert_eq!(parsed.progress_sink, Some(ProgressSink::Stderr));
        assert!(!parsed.allow_screenshots);
        assert_eq!(parsed.cmd_idx, 2);
    }

    #[test]
    fn test_parse_global_args_project_alias() {
        let args = vec![
            "gentle_cli".to_string(),
            "--project".to_string(),
            "project.gentle.json".to_string(),
            "state-summary".to_string(),
        ];
        let parsed = parse_global_args(&args).unwrap();
        assert_eq!(parsed.state_path, "project.gentle.json");
        assert!(!parsed.allow_screenshots);
        assert_eq!(parsed.cmd_idx, 3);
    }

    #[test]
    fn test_parse_global_args_allow_screenshots() {
        let args = vec![
            "gentle_cli".to_string(),
            "--allow-screenshots".to_string(),
            "screenshot-window".to_string(),
            "out.png".to_string(),
        ];
        let parsed = parse_global_args(&args).unwrap();
        assert_eq!(parsed.state_path, DEFAULT_STATE_PATH);
        assert_eq!(parsed.progress_sink, None);
        assert!(parsed.allow_screenshots);
        assert_eq!(parsed.cmd_idx, 2);
    }
}
