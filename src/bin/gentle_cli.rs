//! CLI binary entry point exposing direct commands and shared shell routing.

use gentle::{
    about,
    engine::{
        Engine, EngineStateSummary, GenomeAnnotationScope, GenomeTrackImportProgress, GentleEngine,
        Operation, OperationProgress, ProjectState, RenderSvgMode, RnaReadInterpretProgress,
        TfbsProgress,
    },
    engine_shell::{
        ShellCommand, ShellExecutionOptions, execute_shell_command_with_options, parse_shell_line,
        parse_shell_tokens, parse_workflow_json_payload, shell_help_text,
    },
    genomes::{
        DEFAULT_GENOME_CATALOG_PATH, DEFAULT_HELPER_GENOME_CATALOG_PATH, GenomeGeneRecord,
        PrepareGenomeProgress,
    },
    protocol_cartoon::{ProtocolCartoonKind, protocol_cartoon_catalog_rows},
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
  gentle_cli help [COMMAND ...] [--format text|json|markdown] [--interface all|cli-direct|cli-shell|gui-shell|js|lua|mcp]\n  \
  gentle_cli [--state PATH|--project PATH] [--progress|--progress-stderr|--progress-stdout] COMMAND ...\n\n  \
  gentle_cli [--state PATH|--project PATH] capabilities\n  \
  gentle_cli [--state PATH|--project PATH] op '<operation-json>'|JSON_FILE\n  \
  gentle_cli [--state PATH|--project PATH] workflow '<workflow-json>'|JSON_FILE\n  \
  gentle_cli [--state PATH|--project PATH] state-summary\n  \
  gentle_cli [--state PATH|--project PATH] export-state PATH\n  \
  gentle_cli [--state PATH|--project PATH] import-state PATH\n  \
  gentle_cli [--state PATH|--project PATH] save-project PATH\n  \
  gentle_cli [--state PATH|--project PATH] load-project PATH\n  \
  gentle_cli [--state PATH|--project PATH] render-svg SEQ_ID linear|circular OUTPUT.svg\n  \
  gentle_cli [--state PATH|--project PATH] render-dotplot-svg SEQ_ID DOTPLOT_ID OUTPUT.svg [--flex-track ID] [--display-threshold N] [--intensity-gain N]\n  \
  gentle_cli [--state PATH|--project PATH] inspect-feature-expert SEQ_ID tfbs FEATURE_ID\n  \
  gentle_cli [--state PATH|--project PATH] inspect-feature-expert SEQ_ID restriction CUT_POS_1BASED [--enzyme NAME] [--start START_1BASED] [--end END_1BASED]\n  \
  gentle_cli [--state PATH|--project PATH] inspect-feature-expert SEQ_ID splicing FEATURE_ID\n  \
  gentle_cli [--state PATH|--project PATH] inspect-feature-expert SEQ_ID isoform PANEL_ID\n  \
  gentle_cli [--state PATH|--project PATH] render-feature-expert-svg SEQ_ID tfbs FEATURE_ID OUTPUT.svg\n  \
  gentle_cli [--state PATH|--project PATH] render-feature-expert-svg SEQ_ID restriction CUT_POS_1BASED [--enzyme NAME] [--start START_1BASED] [--end END_1BASED] OUTPUT.svg\n  \
  gentle_cli [--state PATH|--project PATH] render-feature-expert-svg SEQ_ID splicing FEATURE_ID OUTPUT.svg\n  \
  gentle_cli [--state PATH|--project PATH] render-feature-expert-svg SEQ_ID isoform PANEL_ID OUTPUT.svg\n  \
  gentle_cli [--state PATH|--project PATH] panels import-isoform SEQ_ID PANEL_PATH [--panel-id ID] [--strict]\n  \
  gentle_cli [--state PATH|--project PATH] panels inspect-isoform SEQ_ID PANEL_ID\n  \
  gentle_cli [--state PATH|--project PATH] panels render-isoform-svg SEQ_ID PANEL_ID OUTPUT.svg\n  \
  gentle_cli [--state PATH|--project PATH] panels validate-isoform PANEL_PATH [--panel-id ID]\n  \
  gentle_cli [--state PATH|--project PATH] render-rna-svg SEQ_ID OUTPUT.svg\n  \
  gentle_cli [--state PATH|--project PATH] rna-info SEQ_ID\n  \
  gentle_cli [--state PATH|--project PATH] render-lineage-svg OUTPUT.svg\n\n  \
  gentle_cli [--state PATH|--project PATH] protocol-cartoon list\n  \
  gentle_cli [--state PATH|--project PATH] protocol-cartoon render-svg PROTOCOL_ID OUTPUT.svg\n  \
  gentle_cli [--state PATH|--project PATH] protocol-cartoon render-template-svg TEMPLATE.json OUTPUT.svg\n  \
  gentle_cli [--state PATH|--project PATH] protocol-cartoon template-validate TEMPLATE.json\n  \
  gentle_cli [--state PATH|--project PATH] protocol-cartoon render-with-bindings TEMPLATE.json BINDINGS.json OUTPUT.svg\n  \
  gentle_cli [--state PATH|--project PATH] protocol-cartoon template-export PROTOCOL_ID OUTPUT.json\n\n  \
  gentle_cli [--state PATH|--project PATH] gibson preview PLAN_JSON_OR_@FILE [--output OUTPUT.json]\n  \
  gentle_cli [--state PATH|--project PATH] gibson apply PLAN_JSON_OR_@FILE\n\n  \
  gentle_cli [--state PATH|--project PATH] shell 'state-summary'\n  \
  gentle_cli [--state PATH|--project PATH] shell 'op <operation-json>'\n\n  \
  gentle_cli [--state PATH|--project PATH] render-pool-gel-svg IDS|'-' OUTPUT.svg [--ladders NAME[,NAME]] [--containers ID[,ID]] [--arrangement ARR_ID]\n  \
  gentle_cli [--state PATH|--project PATH] render-gel-svg IDS|'-' OUTPUT.svg [--ladders NAME[,NAME]] [--containers ID[,ID]] [--arrangement ARR_ID]\n  \
  gentle_cli [--state PATH|--project PATH] arrange-serial CONTAINER_IDS [--id ARR_ID] [--name TEXT] [--ladders NAME[,NAME]]\n\n  \
  gentle_cli [--state PATH|--project PATH] screenshot-window OUTPUT.png (disabled by security policy)\n\n  \
  gentle_cli [--state PATH|--project PATH] ladders list [--molecule dna|rna] [--filter TEXT]\n  \
  gentle_cli [--state PATH|--project PATH] ladders export OUTPUT.json [--molecule dna|rna] [--filter TEXT]\n\n  \
  gentle_cli [--state PATH|--project PATH] export-pool IDS OUTPUT.pool.gentle.json [HUMAN_ID]\n  \
  gentle_cli [--state PATH|--project PATH] export-run-bundle OUTPUT.run_bundle.json [--run-id RUN_ID]\n  \
  gentle_cli [--state PATH|--project PATH] import-pool INPUT.pool.gentle.json [PREFIX]\n\n  \
  gentle_cli genomes list [--catalog PATH]\n  \
  gentle_cli genomes validate-catalog [--catalog PATH]\n  \
  gentle_cli genomes update-ensembl-specs [--catalog PATH] [--output-catalog PATH]\n  \
  gentle_cli genomes status GENOME_ID [--catalog PATH] [--cache-dir PATH]\n  \
  gentle_cli genomes genes GENOME_ID [--catalog PATH] [--cache-dir PATH] [--filter TEXT] [--limit N] [--offset N]\n  \
  gentle_cli [--state PATH|--project PATH] genomes prepare GENOME_ID [--catalog PATH] [--cache-dir PATH] [--timeout-secs N]\n  \
  gentle_cli genomes remove-prepared GENOME_ID [--catalog PATH] [--cache-dir PATH]\n  \
  gentle_cli genomes remove-catalog-entry GENOME_ID [--catalog PATH] [--output-catalog PATH]\n  \
  gentle_cli genomes blast GENOME_ID QUERY_SEQUENCE [--max-hits N] [--task blastn-short|blastn] [--options-json JSON_OR_@FILE|--options-file PATH] [--catalog PATH] [--cache-dir PATH]\n  \
  gentle_cli [--state PATH|--project PATH] genomes extract-region GENOME_ID CHR START END [--output-id ID] [--annotation-scope none|core|full] [--max-annotation-features N] [--include-genomic-annotation|--no-include-genomic-annotation] [--catalog PATH] [--cache-dir PATH]\n  \
  gentle_cli [--state PATH|--project PATH] genomes extract-gene GENOME_ID QUERY [--occurrence N] [--output-id ID] [--annotation-scope none|core|full] [--max-annotation-features N] [--include-genomic-annotation|--no-include-genomic-annotation] [--catalog PATH] [--cache-dir PATH]\n\n  \
  gentle_cli [--state PATH|--project PATH] genomes extend-anchor SEQ_ID 5p|3p LENGTH_BP [--output-id ID] [--catalog PATH] [--cache-dir PATH] [--prepared-genome GENOME_ID]\n  \
  gentle_cli [--state PATH|--project PATH] genomes verify-anchor SEQ_ID [--catalog PATH] [--cache-dir PATH] [--prepared-genome GENOME_ID]\n\n  \
  gentle_cli helpers list [--catalog PATH]\n  \
  gentle_cli helpers validate-catalog [--catalog PATH]\n  \
  gentle_cli helpers update-ensembl-specs [--catalog PATH] [--output-catalog PATH]\n  \
  gentle_cli helpers status HELPER_ID [--catalog PATH] [--cache-dir PATH]\n  \
  gentle_cli helpers genes HELPER_ID [--catalog PATH] [--cache-dir PATH] [--filter TEXT] [--limit N] [--offset N]\n  \
  gentle_cli [--state PATH|--project PATH] helpers prepare HELPER_ID [--catalog PATH] [--cache-dir PATH] [--timeout-secs N]\n  \
  gentle_cli helpers remove-prepared HELPER_ID [--catalog PATH] [--cache-dir PATH]\n  \
  gentle_cli helpers remove-catalog-entry HELPER_ID [--catalog PATH] [--output-catalog PATH]\n  \
  gentle_cli helpers blast HELPER_ID QUERY_SEQUENCE [--max-hits N] [--task blastn-short|blastn] [--options-json JSON_OR_@FILE|--options-file PATH] [--catalog PATH] [--cache-dir PATH]\n  \
  gentle_cli [--state PATH|--project PATH] helpers extract-region HELPER_ID CHR START END [--output-id ID] [--annotation-scope none|core|full] [--max-annotation-features N] [--include-genomic-annotation|--no-include-genomic-annotation] [--catalog PATH] [--cache-dir PATH]\n  \
  gentle_cli [--state PATH|--project PATH] helpers extract-gene HELPER_ID QUERY [--occurrence N] [--output-id ID] [--annotation-scope none|core|full] [--max-annotation-features N] [--include-genomic-annotation|--no-include-genomic-annotation] [--catalog PATH] [--cache-dir PATH]\n\n  \
  gentle_cli [--state PATH|--project PATH] helpers extend-anchor SEQ_ID 5p|3p LENGTH_BP [--output-id ID] [--catalog PATH] [--cache-dir PATH] [--prepared-genome GENOME_ID]\n  \
  gentle_cli [--state PATH|--project PATH] helpers verify-anchor SEQ_ID [--catalog PATH] [--cache-dir PATH] [--prepared-genome GENOME_ID]\n\n  \
  gentle_cli [--state PATH|--project PATH] tracks import-bed SEQ_ID PATH [--name NAME] [--min-score N] [--max-score N] [--clear-existing]\n\n  \
  gentle_cli [--state PATH|--project PATH] tracks import-bigwig SEQ_ID PATH [--name NAME] [--min-score N] [--max-score N] [--clear-existing]\n\n  \
  gentle_cli [--state PATH|--project PATH] tracks import-vcf SEQ_ID PATH [--name NAME] [--min-score N] [--max-score N] [--clear-existing]\n\n  \
  gentle_cli agents list [--catalog PATH]\n  \
  gentle_cli [--state PATH|--project PATH] agents ask SYSTEM_ID --prompt TEXT [--catalog PATH] [--base-url URL] [--model MODEL] [--timeout-secs N] [--connect-timeout-secs N] [--read-timeout-secs N] [--max-retries N] [--max-response-bytes N] [--allow-auto-exec] [--execute-all] [--execute-index N ...] [--no-state-summary]\n\n  \
  gentle_cli ui intents\n  \
  gentle_cli ui open TARGET [--genome-id GENOME_ID] [--helpers] [--catalog PATH] [--cache-dir PATH] [--filter TEXT] [--species TEXT] [--latest]\n  \
  gentle_cli ui focus TARGET [--genome-id GENOME_ID] [--helpers] [--catalog PATH] [--cache-dir PATH] [--filter TEXT] [--species TEXT] [--latest]\n  \
  gentle_cli ui prepared-genomes [--helpers] [--catalog PATH] [--cache-dir PATH] [--filter TEXT] [--species TEXT] [--latest]\n  \
  gentle_cli ui latest-prepared SPECIES [--helpers] [--catalog PATH] [--cache-dir PATH]\n\n  \
  gentle_cli [--state PATH|--project PATH] macros run SCRIPT_OR_@FILE [--transactional]\n  \
  gentle_cli [--state PATH|--project PATH] macros instance-list\n  \
  gentle_cli [--state PATH|--project PATH] macros instance-show MACRO_INSTANCE_ID\n  \
  gentle_cli [--state PATH|--project PATH] macros template-list\n  \
  gentle_cli [--state PATH|--project PATH] macros template-show TEMPLATE_NAME\n  \
  gentle_cli [--state PATH|--project PATH] macros template-put TEMPLATE_NAME (--script SCRIPT_OR_@FILE|--file PATH) [--description TEXT] [--details-url URL] [--param NAME|NAME=DEFAULT ...] [--input-port PORT_ID:KIND[:one|many][:required|optional][:description]] [--output-port PORT_ID:KIND[:one|many][:required|optional][:description]]\n  \
  gentle_cli [--state PATH|--project PATH] macros template-delete TEMPLATE_NAME\n  \
  gentle_cli [--state PATH|--project PATH] macros template-import PATH\n  \
  gentle_cli [--state PATH|--project PATH] macros template-run TEMPLATE_NAME [--bind KEY=VALUE ...] [--transactional] [--validate-only]\n\n  \
  gentle_cli [--state PATH|--project PATH] candidates list\n  \
  gentle_cli [--state PATH|--project PATH] candidates delete SET_NAME\n  \
  gentle_cli [--state PATH|--project PATH] candidates generate SET_NAME SEQ_ID --length N [--step N] [--feature-kind KIND] [--feature-label-regex REGEX] [--max-distance N] [--feature-geometry feature_span|feature_parts|feature_boundaries] [--feature-boundary any|five_prime|three_prime|start|end] [--strand-relation any|same|opposite] [--limit N]\n  \
  gentle_cli [--state PATH|--project PATH] candidates generate-between-anchors SET_NAME SEQ_ID --length N (--anchor-a-pos N|--anchor-a-json JSON) (--anchor-b-pos N|--anchor-b-json JSON) [--step N] [--limit N]\n  \
  gentle_cli [--state PATH|--project PATH] candidates show SET_NAME [--limit N] [--offset N]\n  \
  gentle_cli [--state PATH|--project PATH] candidates metrics SET_NAME\n  \
  gentle_cli [--state PATH|--project PATH] candidates score SET_NAME METRIC_NAME EXPRESSION\n  \
  gentle_cli [--state PATH|--project PATH] candidates score-distance SET_NAME METRIC_NAME [--feature-kind KIND] [--feature-label-regex REGEX] [--feature-geometry feature_span|feature_parts|feature_boundaries] [--feature-boundary any|five_prime|three_prime|start|end] [--strand-relation any|same|opposite]\n  \
  gentle_cli [--state PATH|--project PATH] candidates score-weighted SET_NAME METRIC_NAME --term METRIC:WEIGHT[:max|min] [--term ...] [--normalize|--no-normalize]\n  \
  gentle_cli [--state PATH|--project PATH] candidates top-k INPUT_SET OUTPUT_SET --metric METRIC_NAME --k N [--direction max|min] [--tie-break seq_start_end|seq_end_start|length_ascending|length_descending|sequence_lexicographic]\n  \
  gentle_cli [--state PATH|--project PATH] candidates pareto INPUT_SET OUTPUT_SET --objective METRIC[:max|min] [--objective ...] [--max-candidates N] [--tie-break seq_start_end|seq_end_start|length_ascending|length_descending|sequence_lexicographic]\n  \
  gentle_cli [--state PATH|--project PATH] candidates filter INPUT_SET OUTPUT_SET --metric METRIC_NAME [--min N] [--max N] [--min-quantile Q] [--max-quantile Q]\n  \
  gentle_cli [--state PATH|--project PATH] candidates set-op union|intersect|subtract LEFT_SET RIGHT_SET OUTPUT_SET\n  \
  gentle_cli [--state PATH|--project PATH] candidates macro SCRIPT_OR_@FILE\n  \
  gentle_cli [--state PATH|--project PATH] candidates template-list\n  \
  gentle_cli [--state PATH|--project PATH] candidates template-show TEMPLATE_NAME\n  \
  gentle_cli [--state PATH|--project PATH] candidates template-put TEMPLATE_NAME (--script SCRIPT_OR_@FILE|--file PATH) [--description TEXT] [--details-url URL] [--param NAME|NAME=DEFAULT ...]\n  \
  gentle_cli [--state PATH|--project PATH] candidates template-delete TEMPLATE_NAME\n  \
  gentle_cli [--state PATH|--project PATH] candidates template-run TEMPLATE_NAME [--bind KEY=VALUE ...] [--transactional]\n\n  \
  gentle_cli [--state PATH|--project PATH] guides list\n  \
  gentle_cli [--state PATH|--project PATH] guides show GUIDE_SET_ID [--limit N] [--offset N]\n  \
  gentle_cli [--state PATH|--project PATH] guides put GUIDE_SET_ID (--json JSON|@FILE|--file PATH)\n  \
  gentle_cli [--state PATH|--project PATH] guides delete GUIDE_SET_ID\n  \
  gentle_cli [--state PATH|--project PATH] guides filter GUIDE_SET_ID [--config JSON|@FILE] [--config-file PATH] [--output-set GUIDE_SET_ID]\n  \
  gentle_cli [--state PATH|--project PATH] guides filter-show GUIDE_SET_ID\n  \
  gentle_cli [--state PATH|--project PATH] guides oligos-generate GUIDE_SET_ID TEMPLATE_ID [--apply-5prime-g-extension] [--output-oligo-set ID] [--passed-only]\n  \
  gentle_cli [--state PATH|--project PATH] guides oligos-list [--guide-set GUIDE_SET_ID]\n  \
  gentle_cli [--state PATH|--project PATH] guides oligos-show OLIGO_SET_ID\n  \
  gentle_cli [--state PATH|--project PATH] guides oligos-export GUIDE_SET_ID OUTPUT_PATH [--format csv_table|plate_csv|fasta] [--plate 96|384] [--oligo-set ID]\n  \
  gentle_cli [--state PATH|--project PATH] guides protocol-export GUIDE_SET_ID OUTPUT_PATH [--oligo-set ID] [--no-qc]\n\n  \
  gentle_cli [--state PATH|--project PATH] primers design REQUEST_JSON_OR_@FILE [--backend auto|internal|primer3] [--primer3-exec PATH]\n  \
  gentle_cli [--state PATH|--project PATH] primers design-qpcr REQUEST_JSON_OR_@FILE [--backend auto|internal|primer3] [--primer3-exec PATH]\n  \
  gentle_cli [--state PATH|--project PATH] primers preflight [--backend auto|internal|primer3] [--primer3-exec PATH]\n  \
  gentle_cli [--state PATH|--project PATH] primers seed-from-feature SEQ_ID FEATURE_ID\n  \
  gentle_cli [--state PATH|--project PATH] primers seed-from-splicing SEQ_ID FEATURE_ID\n  \
  gentle_cli [--state PATH|--project PATH] primers list-reports\n  \
  gentle_cli [--state PATH|--project PATH] primers show-report REPORT_ID\n  \
  gentle_cli [--state PATH|--project PATH] primers export-report REPORT_ID OUTPUT.json\n  \
  gentle_cli [--state PATH|--project PATH] primers list-qpcr-reports\n  \
  gentle_cli [--state PATH|--project PATH] primers show-qpcr-report REPORT_ID\n  \
  gentle_cli [--state PATH|--project PATH] primers export-qpcr-report REPORT_ID OUTPUT.json\n\n  \
  gentle_cli [--state PATH|--project PATH] dotplot compute SEQ_ID [--reference-seq REF_SEQ_ID] [--start N] [--end N] [--ref-start N] [--ref-end N] [--mode self_forward|self_reverse_complement|pair_forward|pair_reverse_complement] [--word-size N] [--step N] [--max-mismatches N] [--tile-bp N] [--id DOTPLOT_ID]\n  \
  gentle_cli [--state PATH|--project PATH] dotplot list [SEQ_ID]\n  \
  gentle_cli [--state PATH|--project PATH] dotplot show DOTPLOT_ID\n  \
  gentle_cli [--state PATH|--project PATH] transcripts derive SEQ_ID [--feature-id N ...] [--scope all_overlapping_both_strands|target_group_any_strand|all_overlapping_target_strand|target_group_target_strand] [--output-prefix PREFIX]\n  \
  gentle_cli [--state PATH|--project PATH] flex compute SEQ_ID [--start N] [--end N] [--model at_richness|at_skew] [--bin-bp N] [--smoothing-bp N] [--id TRACK_ID]\n  \
  gentle_cli [--state PATH|--project PATH] flex list [SEQ_ID]\n  \
  gentle_cli [--state PATH|--project PATH] flex show TRACK_ID\n\n  \
  gentle_cli [--state PATH|--project PATH] splicing-refs derive SEQ_ID START_0BASED END_0BASED [--seed-feature-id N] [--scope all_overlapping_both_strands|target_group_any_strand|all_overlapping_target_strand|target_group_target_strand] [--output-prefix PREFIX]\n  \
  gentle_cli [--state PATH|--project PATH] align compute QUERY_SEQ_ID TARGET_SEQ_ID [--query-start N] [--query-end N] [--target-start N] [--target-end N] [--mode global|local] [--match N] [--mismatch N] [--gap-open N] [--gap-extend N]\n\n  \
  gentle_cli routines list [--catalog PATH] [--family NAME] [--status NAME] [--tag TAG] [--query TEXT]\n  \
  gentle_cli routines explain ROUTINE_ID [--catalog PATH]\n  \
  gentle_cli routines compare ROUTINE_A ROUTINE_B [--catalog PATH]\n\n  \
  gentle_cli planning profile show [--scope global|project_override|confirmed_agent_overlay|effective]\n  \
  gentle_cli planning profile set JSON_OR_@FILE [--scope global|project_override|confirmed_agent_overlay]\n  \
  gentle_cli planning profile clear [--scope global|project_override|confirmed_agent_overlay]\n  \
  gentle_cli planning objective show\n  \
  gentle_cli planning objective set JSON_OR_@FILE\n  \
  gentle_cli planning objective clear\n  \
  gentle_cli planning suggestions list [--status pending|accepted|rejected]\n  \
  gentle_cli planning suggestions accept SUGGESTION_ID\n  \
  gentle_cli planning suggestions reject SUGGESTION_ID [--reason TEXT]\n  \
  gentle_cli planning sync status\n  \
  gentle_cli planning sync pull JSON_OR_@FILE [--source ID] [--confidence N] [--snapshot-id ID]\n  \
  gentle_cli planning sync push JSON_OR_@FILE [--source ID] [--confidence N] [--snapshot-id ID]\n\n  \
  gentle_cli resources sync-rebase INPUT.withrefm [OUTPUT.rebase.json] [--commercial-only]\n  \
  gentle_cli resources sync-jaspar INPUT.jaspar.txt [OUTPUT.motifs.json]\n\n  \
  Tip: pass @file.json instead of inline JSON\n  \
  --project is an alias of --state for project.gentle.json files\n\n  \
  Shell help:\n  \
  {shell_help}",
        shell_help = shell_help_text()
    );
}

const SHELL_FORWARDED_COMMANDS: &[&str] = &[
    "genomes",
    "helpers",
    "agents",
    "ui",
    "routines",
    "planning",
    "gibson",
    "panels",
    "macros",
    "resources",
    "import-pool",
    "ladders",
    "guides",
    "primers",
    "dotplot",
    "transcripts",
    "flex",
    "splicing-refs",
    "align",
    "rna-reads",
    "tracks",
    "genbank",
    "uniprot",
    "screenshot-window",
    "inspect-feature-expert",
    "render-feature-expert-svg",
];

fn is_shell_forwarded_command(command: &str) -> bool {
    SHELL_FORWARDED_COMMANDS.contains(&command)
}

fn parse_forwarded_shell_command(
    args: &[String],
    cmd_idx: usize,
) -> Result<Option<ShellCommand>, String> {
    if args.len() <= cmd_idx {
        return Ok(None);
    }
    let command = &args[cmd_idx];
    if !is_shell_forwarded_command(command.as_str()) {
        return Ok(None);
    }
    let tokens = args[cmd_idx..].to_vec();
    let shell_command = parse_shell_tokens(&tokens)?;
    Ok(Some(shell_command))
}

fn load_json_arg(value: &str) -> Result<String, String> {
    fn strip_shebang_line(raw: &str) -> String {
        if !raw.starts_with("#!") {
            return raw.to_string();
        }
        match raw.split_once('\n') {
            Some((_, rest)) => rest.to_string(),
            None => String::new(),
        }
    }

    let trimmed = value.trim();
    if let Some(path) = trimmed.strip_prefix('@') {
        let path = path.trim();
        if path.is_empty() {
            return Err("Could not read JSON file '': empty path".to_string());
        }
        let text = fs::read_to_string(path)
            .map_err(|e| format!("Could not read JSON file '{path}': {e}"))?;
        return Ok(strip_shebang_line(&text));
    }
    if trimmed.starts_with('{') || trimmed.starts_with('[') {
        return Ok(value.to_string());
    }
    let candidate = Path::new(trimmed);
    if candidate.is_file() {
        let text = fs::read_to_string(candidate).map_err(|e| {
            format!(
                "Could not read JSON file '{}': {e}",
                candidate.to_string_lossy()
            )
        })?;
        return Ok(strip_shebang_line(&text));
    }
    Ok(value.to_string())
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
    if !std::path::Path::new(path).exists() {
        return Ok(ProjectState::default());
    }
    let raw =
        fs::read_to_string(path).map_err(|e| format!("Could not read state file '{path}': {e}"))?;
    if raw.trim().is_empty() {
        return Ok(ProjectState::default());
    }
    ProjectState::load_from_path(path).map_err(|e| e.to_string())
}

fn print_json<T: Serialize>(value: &T) -> Result<(), String> {
    let text = serde_json::to_string_pretty(value)
        .map_err(|e| format!("Could not serialize JSON output: {e}"))?;
    println!("{text}");
    Ok(())
}

#[allow(dead_code)]
fn print_help_output(output: &serde_json::Value) -> Result<(), String> {
    if let Some(help) = output.get("help").and_then(|v| v.as_str()) {
        println!("{help}");
        return Ok(());
    }
    if let Some(help_markdown) = output.get("help_markdown").and_then(|v| v.as_str()) {
        println!("{help_markdown}");
        return Ok(());
    }
    print_json(output)
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
    let allow_screenshots = false;
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
                return Err("--allow-screenshots is disabled by security policy".to_string());
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

    fn on_genome_track_import_progress(&mut self, p: GenomeTrackImportProgress) {
        if p.done || p.parsed_records % 1_000 == 0 {
            self.print_line(&format!(
                "progress track-import seq={} source={} parsed={} imported={} skipped={} done={} path={}",
                p.seq_id,
                p.source,
                p.parsed_records,
                p.imported_features,
                p.skipped_records,
                p.done,
                p.path
            ));
        }
    }

    fn on_rna_read_interpret_progress(&mut self, p: RnaReadInterpretProgress) {
        let stride = p.update_every_reads.max(1);
        if p.done || p.reads_processed % stride == 0 {
            let percent = if p.reads_total == 0 {
                100.0
            } else {
                (p.reads_processed as f64 / p.reads_total as f64) * 100.0
            };
            self.print_line(&format!(
                "progress rna-reads seq={} reads={}/{} ({:.1}%) seed_passed={} aligned={} matched_kmers={} tested_kmers={} done={}",
                p.seq_id,
                p.reads_processed,
                p.reads_total,
                percent,
                p.seed_passed,
                p.aligned,
                p.matched_kmers,
                p.tested_kmers,
                p.done
            ));
        }
    }

    fn on_progress(&mut self, progress: OperationProgress) {
        match progress {
            OperationProgress::Tfbs(p) => self.on_tfbs_progress(p),
            OperationProgress::GenomePrepare(p) => self.on_genome_prepare_progress(p),
            OperationProgress::GenomeTrackImport(p) => self.on_genome_track_import_progress(p),
            OperationProgress::RnaReadInterpret(p) => self.on_rna_read_interpret_progress(p),
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
        allow_agent_commands: true,
    };
    if args.len() <= cmd_idx {
        usage();
        return Err("Missing command".to_string());
    }

    let command = &args[cmd_idx];

    if let Some(shell_command) = parse_forwarded_shell_command(&args, cmd_idx)? {
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
        "help" => {
            let tokens = args[cmd_idx..].to_vec();
            let shell_command = parse_shell_tokens(&tokens)?;
            if !matches!(shell_command, ShellCommand::Help { .. }) {
                return Err("Expected help command".to_string());
            }
            let mut engine = GentleEngine::from_state(ProjectState::default());
            let run =
                execute_shell_command_with_options(&mut engine, &shell_command, &shell_options)?;
            print_help_output(&run.output)
        }
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
                                return Err(format!("Unknown option '{}' for {label} list", other));
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
                "validate-catalog" => {
                    let mut catalog_path: Option<String> = None;
                    let mut idx = cmd_idx + 2;
                    while idx < args.len() {
                        match args[idx].as_str() {
                            "--catalog" => {
                                if idx + 1 >= args.len() {
                                    return Err(format!(
                                        "Missing PATH after --catalog for {label} validate-catalog"
                                    ));
                                }
                                catalog_path = Some(args[idx + 1].clone());
                                idx += 2;
                            }
                            other => {
                                return Err(format!(
                                    "Unknown option '{}' for {label} validate-catalog",
                                    other
                                ));
                            }
                        }
                    }
                    let resolved_catalog = catalog_path
                        .as_deref()
                        .map(str::trim)
                        .filter(|v| !v.is_empty())
                        .or_else(|| helper_mode.then_some(default_catalog))
                        .ok_or_else(|| format!("Missing catalog path for {label} validate-catalog"))?;
                    let genomes = GentleEngine::list_reference_genomes(Some(resolved_catalog))
                        .map_err(|e| e.to_string())?;
                    print_json(&json!({
                        "catalog_path": resolved_catalog,
                        "valid": true,
                        "genome_count": genomes.len(),
                    }))
                }
                "update-ensembl-specs" => {
                    let mut catalog_path: Option<String> = None;
                    let mut output_catalog_path: Option<String> = None;
                    let mut idx = cmd_idx + 2;
                    while idx < args.len() {
                        match args[idx].as_str() {
                            "--catalog" => {
                                if idx + 1 >= args.len() {
                                    return Err(format!(
                                        "Missing PATH after --catalog for {label} update-ensembl-specs"
                                    ));
                                }
                                catalog_path = Some(args[idx + 1].clone());
                                idx += 2;
                            }
                            "--output-catalog" => {
                                if idx + 1 >= args.len() {
                                    return Err(format!(
                                        "Missing PATH after --output-catalog for {label} update-ensembl-specs"
                                    ));
                                }
                                output_catalog_path = Some(args[idx + 1].clone());
                                idx += 2;
                            }
                            other => {
                                return Err(format!(
                                    "Unknown option '{}' for {label} update-ensembl-specs",
                                    other
                                ));
                            }
                        }
                    }
                    let resolved_catalog = catalog_path
                        .as_deref()
                        .map(str::trim)
                        .filter(|v| !v.is_empty())
                        .or_else(|| helper_mode.then_some(default_catalog));
                    let report = GentleEngine::apply_reference_genome_ensembl_catalog_updates(
                        resolved_catalog,
                        output_catalog_path.as_deref(),
                    )
                    .map_err(|e| e.to_string())?;
                    print_json(&report)
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
                                ));
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
                        "nucleotide_length_bp": source_plan.nucleotide_length_bp,
                        "molecular_mass_da": source_plan.molecular_mass_da,
                        "molecular_mass_source": source_plan.molecular_mass_source,
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
                                    return Err(format!(
                                        "Missing N after --limit for {label} genes"
                                    ));
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
                                ));
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
                            "{label} prepare requires GENOME_ID [--catalog PATH] [--cache-dir PATH] [--timeout-secs N]"
                        ));
                    }
                    let genome_id = args[cmd_idx + 2].clone();
                    let mut catalog_path: Option<String> = None;
                    let mut cache_dir: Option<String> = None;
                    let mut timeout_seconds: Option<u64> = None;
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
                            "--timeout-secs" => {
                                if idx + 1 >= args.len() {
                                    return Err(format!(
                                        "Missing N after --timeout-secs for {label} prepare"
                                    ));
                                }
                                let raw = args[idx + 1].trim().to_string();
                                let parsed = raw.parse::<u64>().map_err(|e| {
                                    format!(
                                        "Invalid --timeout-secs value '{}' for {label} prepare: {}",
                                        raw, e
                                    )
                                })?;
                                if parsed == 0 {
                                    timeout_seconds = None;
                                } else {
                                    timeout_seconds = Some(parsed);
                                }
                                idx += 2;
                            }
                            other => {
                                return Err(format!(
                                    "Unknown option '{}' for {label} prepare",
                                    other
                                ));
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
                        timeout_seconds,
                    };
                    let result = if let Some(sink) = global.progress_sink {
                        let mut printer = ProgressPrinter::new(sink);
                        engine
                            .apply_with_progress(op, |p| {
                                printer.on_progress(p);
                                true
                            })
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
                "remove-prepared" => {
                    if args.len() <= cmd_idx + 2 {
                        usage();
                        return Err(format!(
                            "{label} remove-prepared requires GENOME_ID [--catalog PATH] [--cache-dir PATH]"
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
                                        "Missing PATH after --catalog for {label} remove-prepared"
                                    ));
                                }
                                catalog_path = Some(args[idx + 1].clone());
                                idx += 2;
                            }
                            "--cache-dir" => {
                                if idx + 1 >= args.len() {
                                    return Err(format!(
                                        "Missing PATH after --cache-dir for {label} remove-prepared"
                                    ));
                                }
                                cache_dir = Some(args[idx + 1].clone());
                                idx += 2;
                            }
                            other => {
                                return Err(format!(
                                    "Unknown option '{}' for {label} remove-prepared",
                                    other
                                ));
                            }
                        }
                    }
                    let resolved_catalog = catalog_path
                        .as_deref()
                        .map(str::trim)
                        .filter(|v| !v.is_empty())
                        .or_else(|| helper_mode.then_some(default_catalog));
                    let report = GentleEngine::remove_prepared_reference_genome(
                        resolved_catalog,
                        &genome_id,
                        cache_dir.as_deref(),
                    )
                    .map_err(|e| e.to_string())?;
                    print_json(&report)
                }
                "remove-catalog-entry" => {
                    if args.len() <= cmd_idx + 2 {
                        usage();
                        return Err(format!(
                            "{label} remove-catalog-entry requires GENOME_ID [--catalog PATH] [--output-catalog PATH]"
                        ));
                    }
                    let genome_id = args[cmd_idx + 2].clone();
                    let mut catalog_path: Option<String> = None;
                    let mut output_catalog_path: Option<String> = None;
                    let mut idx = cmd_idx + 3;
                    while idx < args.len() {
                        match args[idx].as_str() {
                            "--catalog" => {
                                if idx + 1 >= args.len() {
                                    return Err(format!(
                                        "Missing PATH after --catalog for {label} remove-catalog-entry"
                                    ));
                                }
                                catalog_path = Some(args[idx + 1].clone());
                                idx += 2;
                            }
                            "--output-catalog" => {
                                if idx + 1 >= args.len() {
                                    return Err(format!(
                                        "Missing PATH after --output-catalog for {label} remove-catalog-entry"
                                    ));
                                }
                                output_catalog_path = Some(args[idx + 1].clone());
                                idx += 2;
                            }
                            other => {
                                return Err(format!(
                                    "Unknown option '{}' for {label} remove-catalog-entry",
                                    other
                                ));
                            }
                        }
                    }
                    let resolved_catalog = catalog_path
                        .as_deref()
                        .map(str::trim)
                        .filter(|v| !v.is_empty())
                        .or_else(|| helper_mode.then_some(default_catalog));
                    let report = GentleEngine::remove_reference_genome_catalog_entry(
                        resolved_catalog,
                        &genome_id,
                        output_catalog_path.as_deref(),
                    )
                    .map_err(|e| e.to_string())?;
                    print_json(&report)
                }
                "extract-region" => {
                    if args.len() <= cmd_idx + 6 {
                        usage();
                        return Err(format!(
                            "{label} extract-region requires GENOME_ID CHR START END [--output-id ID] [--annotation-scope none|core|full] [--max-annotation-features N] [--include-genomic-annotation|--no-include-genomic-annotation] [--catalog PATH] [--cache-dir PATH]"
                        ));
                    }
                    let genome_id = args[cmd_idx + 2].clone();
                    let chromosome = args[cmd_idx + 3].clone();
                    let start_1based = args[cmd_idx + 4].parse::<usize>().map_err(|e| {
                        format!("Invalid START coordinate '{}': {}", args[cmd_idx + 4], e)
                    })?;
                    let end_1based = args[cmd_idx + 5].parse::<usize>().map_err(|e| {
                        format!("Invalid END coordinate '{}': {}", args[cmd_idx + 5], e)
                    })?;
                    let mut output_id: Option<String> = None;
                    let mut annotation_scope: Option<GenomeAnnotationScope> = None;
                    let mut max_annotation_features: Option<usize> = None;
                    let mut include_genomic_annotation: Option<bool> = None;
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
                            "--annotation-scope" => {
                                if idx + 1 >= args.len() {
                                    return Err(format!(
                                        "Missing VALUE after --annotation-scope for {label} extract-region"
                                    ));
                                }
                                let value = args[idx + 1].trim().to_ascii_lowercase();
                                let parsed = match value.as_str() {
                                    "none" => GenomeAnnotationScope::None,
                                    "core" => GenomeAnnotationScope::Core,
                                    "full" => GenomeAnnotationScope::Full,
                                    other => {
                                        return Err(format!(
                                            "Invalid --annotation-scope value '{}' for {label} extract-region (expected none|core|full)",
                                            other
                                        ));
                                    }
                                };
                                annotation_scope = Some(parsed);
                                idx += 2;
                            }
                            "--max-annotation-features" => {
                                if idx + 1 >= args.len() {
                                    return Err(format!(
                                        "Missing N after --max-annotation-features for {label} extract-region"
                                    ));
                                }
                                let raw = args[idx + 1].trim().to_string();
                                let parsed = raw.parse::<usize>().map_err(|e| {
                                    format!(
                                        "Invalid --max-annotation-features value '{}' for {label} extract-region: {}",
                                        raw, e
                                    )
                                })?;
                                max_annotation_features = Some(parsed);
                                idx += 2;
                            }
                            "--include-genomic-annotation" => {
                                include_genomic_annotation = Some(true);
                                idx += 1;
                            }
                            "--no-include-genomic-annotation" => {
                                include_genomic_annotation = Some(false);
                                idx += 1;
                            }
                            other => {
                                return Err(format!(
                                    "Unknown option '{}' for {label} extract-region",
                                    other
                                ));
                            }
                        }
                    }
                    if let Some(include) = include_genomic_annotation {
                        let mapped_scope = if include {
                            GenomeAnnotationScope::Core
                        } else {
                            GenomeAnnotationScope::None
                        };
                        if let Some(explicit_scope) = annotation_scope {
                            if explicit_scope != mapped_scope {
                                return Err(format!(
                                    "Conflicting annotation options for {label} extract-region: --annotation-scope={} with legacy include/no-include flag",
                                    match explicit_scope {
                                        GenomeAnnotationScope::None => "none",
                                        GenomeAnnotationScope::Core => "core",
                                        GenomeAnnotationScope::Full => "full",
                                    }
                                ));
                            }
                        } else {
                            annotation_scope = Some(mapped_scope);
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
                            annotation_scope,
                            max_annotation_features,
                            include_genomic_annotation,
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
                            "{label} extract-gene requires GENOME_ID QUERY [--occurrence N] [--output-id ID] [--annotation-scope none|core|full] [--max-annotation-features N] [--include-genomic-annotation|--no-include-genomic-annotation] [--catalog PATH] [--cache-dir PATH]"
                        ));
                    }
                    let genome_id = args[cmd_idx + 2].clone();
                    let gene_query = args[cmd_idx + 3].clone();
                    let mut occurrence: Option<usize> = None;
                    let mut output_id: Option<String> = None;
                    let mut annotation_scope: Option<GenomeAnnotationScope> = None;
                    let mut max_annotation_features: Option<usize> = None;
                    let mut include_genomic_annotation: Option<bool> = None;
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
                                    format!("Invalid --occurrence value '{}': {}", args[idx + 1], e)
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
                            "--annotation-scope" => {
                                if idx + 1 >= args.len() {
                                    return Err(format!(
                                        "Missing VALUE after --annotation-scope for {label} extract-gene"
                                    ));
                                }
                                let parsed = match args[idx + 1]
                                    .trim()
                                    .to_ascii_lowercase()
                                    .as_str()
                                {
                                    "none" => GenomeAnnotationScope::None,
                                    "core" => GenomeAnnotationScope::Core,
                                    "full" => GenomeAnnotationScope::Full,
                                    other => {
                                        return Err(format!(
                                            "Invalid --annotation-scope value '{}' for {label} extract-gene (expected none|core|full)",
                                            other
                                        ));
                                    }
                                };
                                annotation_scope = Some(parsed);
                                idx += 2;
                            }
                            "--max-annotation-features" => {
                                if idx + 1 >= args.len() {
                                    return Err(format!(
                                        "Missing N after --max-annotation-features for {label} extract-gene"
                                    ));
                                }
                                let parsed = args[idx + 1].parse::<usize>().map_err(|e| {
                                    format!(
                                        "Invalid --max-annotation-features value '{}' for {label} extract-gene: {}",
                                        args[idx + 1],
                                        e
                                    )
                                })?;
                                max_annotation_features = Some(parsed);
                                idx += 2;
                            }
                            "--include-genomic-annotation" => {
                                include_genomic_annotation = Some(true);
                                idx += 1;
                            }
                            "--no-include-genomic-annotation" => {
                                include_genomic_annotation = Some(false);
                                idx += 1;
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
                                ));
                            }
                        }
                    }
                    if let Some(include) = include_genomic_annotation {
                        let mapped_scope = if include {
                            GenomeAnnotationScope::Core
                        } else {
                            GenomeAnnotationScope::None
                        };
                        if let Some(explicit_scope) = annotation_scope {
                            if explicit_scope != mapped_scope {
                                return Err(format!(
                                    "Conflicting annotation options for {label} extract-gene: --annotation-scope={} with legacy include/no-include flag",
                                    explicit_scope.as_str()
                                ));
                            }
                        } else {
                            annotation_scope = Some(mapped_scope);
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
                            annotation_scope,
                            max_annotation_features,
                            include_genomic_annotation,
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
                    "Unknown {label} subcommand '{}' (expected list, validate-catalog, update-ensembl-specs, status, genes, prepare, remove-prepared, remove-catalog-entry, extract-region, extract-gene)",
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
        "candidates" => {
            if args.len() <= cmd_idx + 1 {
                usage();
                return Err("candidates requires a subcommand".to_string());
            }
            let tokens = args[cmd_idx..].to_vec();
            let shell_command = parse_shell_tokens(&tokens)?;
            let is_candidates_command = matches!(
                shell_command,
                ShellCommand::CandidatesList
                    | ShellCommand::CandidatesDelete { .. }
                    | ShellCommand::CandidatesGenerate { .. }
                    | ShellCommand::CandidatesGenerateBetweenAnchors { .. }
                    | ShellCommand::CandidatesShow { .. }
                    | ShellCommand::CandidatesMetrics { .. }
                    | ShellCommand::CandidatesScoreExpression { .. }
                    | ShellCommand::CandidatesScoreDistance { .. }
                    | ShellCommand::CandidatesScoreWeightedObjective { .. }
                    | ShellCommand::CandidatesTopK { .. }
                    | ShellCommand::CandidatesParetoFrontier { .. }
                    | ShellCommand::CandidatesFilter { .. }
                    | ShellCommand::CandidatesSetOp { .. }
                    | ShellCommand::CandidatesMacro { .. }
                    | ShellCommand::CandidatesTemplateList
                    | ShellCommand::CandidatesTemplateShow { .. }
                    | ShellCommand::CandidatesTemplateUpsert { .. }
                    | ShellCommand::CandidatesTemplateDelete { .. }
                    | ShellCommand::CandidatesTemplateRun { .. }
            );
            if !is_candidates_command {
                return Err("Expected a candidates subcommand".to_string());
            }
            let mut engine = GentleEngine::from_state(load_state(&state_path)?);
            let run =
                execute_shell_command_with_options(&mut engine, &shell_command, &shell_options)?;
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
                    ));
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
        "render-dotplot-svg" => {
            if args.len() <= cmd_idx + 3 {
                usage();
                return Err(
                    "render-dotplot-svg requires: SEQ_ID DOTPLOT_ID OUTPUT.svg [--flex-track ID] [--display-threshold N] [--intensity-gain N]".to_string(),
                );
            }
            let seq_id = args[cmd_idx + 1].trim();
            let dotplot_id = args[cmd_idx + 2].trim();
            if seq_id.is_empty() {
                return Err("render-dotplot-svg requires non-empty SEQ_ID".to_string());
            }
            if dotplot_id.is_empty() {
                return Err("render-dotplot-svg requires non-empty DOTPLOT_ID".to_string());
            }
            let output = &args[cmd_idx + 3];
            let mut flex_track_id: Option<String> = None;
            let mut display_density_threshold: Option<f32> = None;
            let mut display_intensity_gain: Option<f32> = None;
            let mut idx = cmd_idx + 4;
            while idx < args.len() {
                match args[idx].as_str() {
                    "--flex-track" => {
                        if idx + 1 >= args.len() {
                            return Err("Missing value after --flex-track".to_string());
                        }
                        let value = args[idx + 1].trim();
                        if !value.is_empty() {
                            flex_track_id = Some(value.to_string());
                        }
                        idx += 2;
                    }
                    "--display-threshold" => {
                        if idx + 1 >= args.len() {
                            return Err("Missing value after --display-threshold".to_string());
                        }
                        let raw = args[idx + 1].trim();
                        let parsed = raw.parse::<f32>().map_err(|_| {
                            format!(
                                "Invalid --display-threshold '{}': expected decimal number",
                                raw
                            )
                        })?;
                        display_density_threshold = Some(parsed);
                        idx += 2;
                    }
                    "--intensity-gain" => {
                        if idx + 1 >= args.len() {
                            return Err("Missing value after --intensity-gain".to_string());
                        }
                        let raw = args[idx + 1].trim();
                        let parsed = raw.parse::<f32>().map_err(|_| {
                            format!(
                                "Invalid --intensity-gain '{}': expected decimal number",
                                raw
                            )
                        })?;
                        display_intensity_gain = Some(parsed);
                        idx += 2;
                    }
                    other => {
                        return Err(format!(
                            "Unknown argument '{other}' for render-dotplot-svg (expected --flex-track/--display-threshold/--intensity-gain)"
                        ));
                    }
                }
            }
            let mut engine = GentleEngine::from_state(load_state(&state_path)?);
            let result = engine
                .apply(Operation::RenderDotplotSvg {
                    seq_id: seq_id.to_string(),
                    dotplot_id: dotplot_id.to_string(),
                    path: output.to_string(),
                    flex_track_id,
                    display_density_threshold,
                    display_intensity_gain,
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
        "protocol-cartoon" => {
            if args.len() <= cmd_idx + 1 {
                usage();
                return Err(
                    "protocol-cartoon requires a subcommand: list, render-svg, render-template-svg, template-validate, render-with-bindings, template-export"
                        .to_string(),
                );
            }
            match args[cmd_idx + 1].as_str() {
                "list" => {
                    if args.len() != cmd_idx + 2 {
                        return Err(
                            "protocol-cartoon list takes no additional arguments".to_string()
                        );
                    }
                    print_json(&protocol_cartoon_catalog_rows())
                }
                "render-svg" => {
                    if args.len() != cmd_idx + 4 {
                        return Err(
                            "protocol-cartoon render-svg requires: PROTOCOL_ID OUTPUT.svg"
                                .to_string(),
                        );
                    }
                    let protocol_id = args[cmd_idx + 2].trim();
                    if protocol_id.is_empty() {
                        return Err("protocol-cartoon render-svg requires non-empty PROTOCOL_ID"
                            .to_string());
                    }
                    let output = &args[cmd_idx + 3];
                    let protocol = ProtocolCartoonKind::parse_id(protocol_id).ok_or_else(|| {
                        format!(
                            "Unknown protocol cartoon '{protocol_id}' (run: protocol-cartoon list)"
                        )
                    })?;
                    let mut engine = GentleEngine::from_state(load_state(&state_path)?);
                    let result = engine
                        .apply(Operation::RenderProtocolCartoonSvg {
                            protocol,
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
                "render-template-svg" => {
                    if args.len() != cmd_idx + 4 {
                        return Err(
                            "protocol-cartoon render-template-svg requires: TEMPLATE.json OUTPUT.svg"
                                .to_string(),
                        );
                    }
                    let template_path = args[cmd_idx + 2].trim();
                    if template_path.is_empty() {
                        return Err(
                            "protocol-cartoon render-template-svg requires non-empty TEMPLATE.json"
                                .to_string(),
                        );
                    }
                    let output = &args[cmd_idx + 3];
                    let mut engine = GentleEngine::from_state(load_state(&state_path)?);
                    let result = engine
                        .apply(Operation::RenderProtocolCartoonTemplateSvg {
                            template_path: args[cmd_idx + 2].clone(),
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
                "template-validate" => {
                    if args.len() != cmd_idx + 3 {
                        return Err("protocol-cartoon template-validate requires: TEMPLATE.json"
                            .to_string());
                    }
                    let template_path = args[cmd_idx + 2].trim();
                    if template_path.is_empty() {
                        return Err(
                            "protocol-cartoon template-validate requires non-empty TEMPLATE.json"
                                .to_string(),
                        );
                    }
                    let mut engine = GentleEngine::from_state(load_state(&state_path)?);
                    let result = engine
                        .apply(Operation::ValidateProtocolCartoonTemplate {
                            template_path: args[cmd_idx + 2].clone(),
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
                "render-with-bindings" => {
                    if args.len() != cmd_idx + 5 {
                        return Err(
                            "protocol-cartoon render-with-bindings requires: TEMPLATE.json BINDINGS.json OUTPUT.svg"
                                .to_string(),
                        );
                    }
                    let template_path = args[cmd_idx + 2].trim();
                    if template_path.is_empty() {
                        return Err(
                            "protocol-cartoon render-with-bindings requires non-empty TEMPLATE.json"
                                .to_string(),
                        );
                    }
                    let bindings_path = args[cmd_idx + 3].trim();
                    if bindings_path.is_empty() {
                        return Err(
                            "protocol-cartoon render-with-bindings requires non-empty BINDINGS.json"
                                .to_string(),
                        );
                    }
                    let output = &args[cmd_idx + 4];
                    let mut engine = GentleEngine::from_state(load_state(&state_path)?);
                    let result = engine
                        .apply(Operation::RenderProtocolCartoonTemplateWithBindingsSvg {
                            template_path: args[cmd_idx + 2].clone(),
                            bindings_path: args[cmd_idx + 3].clone(),
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
                "template-export" => {
                    if args.len() != cmd_idx + 4 {
                        return Err(
                            "protocol-cartoon template-export requires: PROTOCOL_ID OUTPUT.json"
                                .to_string(),
                        );
                    }
                    let protocol_id = args[cmd_idx + 2].trim();
                    if protocol_id.is_empty() {
                        return Err(
                            "protocol-cartoon template-export requires non-empty PROTOCOL_ID"
                                .to_string(),
                        );
                    }
                    let protocol = ProtocolCartoonKind::parse_id(protocol_id).ok_or_else(|| {
                        format!(
                            "Unknown protocol cartoon '{protocol_id}' (run: protocol-cartoon list)"
                        )
                    })?;
                    let output = &args[cmd_idx + 3];
                    let mut engine = GentleEngine::from_state(load_state(&state_path)?);
                    let result = engine
                        .apply(Operation::ExportProtocolCartoonTemplateJson {
                            protocol,
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
                other => Err(format!(
                    "Unknown protocol-cartoon subcommand '{other}' (expected list, render-svg, render-template-svg, template-validate, render-with-bindings, template-export)"
                )),
            }
        }
        "render-pool-gel-svg" | "render-gel-svg" => {
            let cmd_name = args[cmd_idx].as_str();
            if args.len() <= cmd_idx + 2 {
                usage();
                return Err(format!(
                    "{cmd_name} requires: IDS|'-' OUTPUT.svg [--ladders NAME[,NAME]] [--containers ID[,ID]] [--arrangement ARR_ID]"
                ));
            }
            let ids = match args[cmd_idx + 1].trim() {
                "-" | "_" => vec![],
                raw => raw
                    .split(',')
                    .map(|s| s.trim())
                    .filter(|s| !s.is_empty())
                    .map(|s| s.to_string())
                    .collect::<Vec<_>>(),
            };
            let output = &args[cmd_idx + 2];
            let mut ladders: Option<Vec<String>> = None;
            let mut container_ids: Option<Vec<String>> = None;
            let mut arrangement_id: Option<String> = None;
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
                    "--containers" => {
                        if idx + 1 >= args.len() {
                            return Err("Missing value after --containers".to_string());
                        }
                        let parsed = args[idx + 1]
                            .split(',')
                            .map(|s| s.trim())
                            .filter(|s| !s.is_empty())
                            .map(|s| s.to_string())
                            .collect::<Vec<_>>();
                        if !parsed.is_empty() {
                            container_ids = Some(parsed);
                        }
                        idx += 2;
                    }
                    "--arrangement" => {
                        if idx + 1 >= args.len() {
                            return Err("Missing value after --arrangement".to_string());
                        }
                        let value = args[idx + 1].trim();
                        if !value.is_empty() {
                            arrangement_id = Some(value.to_string());
                        }
                        idx += 2;
                    }
                    other => {
                        return Err(format!(
                            "Unknown argument '{other}' for {cmd_name} (expected --ladders/--containers/--arrangement)"
                        ));
                    }
                }
            }
            if ids.is_empty()
                && container_ids.as_ref().map_or(true, |v| v.is_empty())
                && arrangement_id.is_none()
            {
                return Err(format!(
                    "{cmd_name} requires sequence ids, --containers, or --arrangement"
                ));
            }
            let mut engine = GentleEngine::from_state(load_state(&state_path)?);
            let result = engine
                .apply(Operation::RenderPoolGelSvg {
                    inputs: ids,
                    path: output.to_string(),
                    ladders,
                    container_ids,
                    arrangement_id,
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
        "arrange-serial" => {
            if args.len() <= cmd_idx + 1 {
                usage();
                return Err(
                    "arrange-serial requires: CONTAINER_IDS [--id ARR_ID] [--name TEXT] [--ladders NAME[,NAME]]"
                        .to_string(),
                );
            }
            let container_ids = args[cmd_idx + 1]
                .split(',')
                .map(|s| s.trim())
                .filter(|s| !s.is_empty())
                .map(|s| s.to_string())
                .collect::<Vec<_>>();
            if container_ids.is_empty() {
                return Err("arrange-serial requires at least one container id".to_string());
            }
            let mut arrangement_id: Option<String> = None;
            let mut name: Option<String> = None;
            let mut ladders: Option<Vec<String>> = None;
            let mut idx = cmd_idx + 2;
            while idx < args.len() {
                match args[idx].as_str() {
                    "--id" => {
                        if idx + 1 >= args.len() {
                            return Err("Missing value after --id".to_string());
                        }
                        let value = args[idx + 1].trim();
                        if !value.is_empty() {
                            arrangement_id = Some(value.to_string());
                        }
                        idx += 2;
                    }
                    "--name" => {
                        if idx + 1 >= args.len() {
                            return Err("Missing value after --name".to_string());
                        }
                        let value = args[idx + 1].trim();
                        if !value.is_empty() {
                            name = Some(value.to_string());
                        }
                        idx += 2;
                    }
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
                            "Unknown argument '{other}' for arrange-serial (expected --id/--name/--ladders)"
                        ));
                    }
                }
            }
            let mut engine = GentleEngine::from_state(load_state(&state_path)?);
            let result = engine
                .apply(Operation::CreateArrangementSerial {
                    container_ids,
                    arrangement_id,
                    name,
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
        "export-run-bundle" => {
            if args.len() <= cmd_idx + 1 {
                usage();
                return Err(
                    "export-run-bundle requires: OUTPUT.run_bundle.json [--run-id RUN_ID]"
                        .to_string(),
                );
            }
            let output = args[cmd_idx + 1].clone();
            let mut run_id: Option<String> = None;
            let mut idx = cmd_idx + 2;
            while idx < args.len() {
                match args[idx].as_str() {
                    "--run-id" => {
                        if idx + 1 >= args.len() {
                            return Err("Missing value after --run-id".to_string());
                        }
                        let value = args[idx + 1].trim();
                        if !value.is_empty() {
                            run_id = Some(value.to_string());
                        }
                        idx += 2;
                    }
                    other => {
                        return Err(format!(
                            "Unknown argument '{other}' for export-run-bundle (expected --run-id)"
                        ));
                    }
                }
            }
            let mut engine = GentleEngine::from_state(load_state(&state_path)?);
            let result = engine
                .apply(Operation::ExportProcessRunBundle {
                    path: output,
                    run_id,
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
                    .apply_with_progress(op, |p| {
                        printer.on_progress(p);
                        true
                    })
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
            let workflow = parse_workflow_json_payload(&json)?;

            let mut engine = GentleEngine::from_state(load_state(&state_path)?);
            let results = if let Some(sink) = global.progress_sink {
                let mut printer = ProgressPrinter::new(sink);
                engine
                    .apply_workflow_with_progress(workflow, |p| {
                        printer.on_progress(p);
                        true
                    })
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
    use gentle::test_support::{
        demo_jaspar_pfm_text, demo_rebase_withrefm_text, write_demo_jaspar_pfm,
        write_demo_pool_json, write_demo_rebase_withrefm, write_demo_workflow_json,
        write_demo_workflow_with_shebang,
    };
    use std::fs;
    use std::sync::Mutex;
    use tempfile::tempdir;

    static TEST_ENV_LOCK: Mutex<()> = Mutex::new(());

    struct EnvVarGuard {
        key: String,
        previous: Option<String>,
    }

    impl EnvVarGuard {
        fn set(key: &str, value: &str) -> Self {
            let previous = std::env::var(key).ok();
            // SAFETY: Test-only mutation is serialized by TEST_ENV_LOCK.
            unsafe { std::env::set_var(key, value) };
            Self {
                key: key.to_string(),
                previous,
            }
        }
    }

    impl Drop for EnvVarGuard {
        fn drop(&mut self) {
            if let Some(value) = &self.previous {
                // SAFETY: Test-only mutation is serialized by TEST_ENV_LOCK.
                unsafe { std::env::set_var(&self.key, value) };
            } else {
                // SAFETY: Test-only mutation is serialized by TEST_ENV_LOCK.
                unsafe { std::env::remove_var(&self.key) };
            }
        }
    }

    fn execute_forwarded_like_cli(
        state: ProjectState,
        args: Vec<String>,
    ) -> (bool, serde_json::Value, ProjectState) {
        let parsed = parse_forwarded_shell_command(&args, 1)
            .expect("parse forwarded command")
            .expect("expected forwarded shell command");
        let mut engine = GentleEngine::from_state(state);
        let run = execute_shell_command_with_options(
            &mut engine,
            &parsed,
            &ShellExecutionOptions {
                allow_screenshots: false,
                allow_agent_commands: true,
            },
        )
        .expect("execute forwarded command");
        (run.state_changed, run.output, engine.state().clone())
    }

    fn execute_shared_shell_tokens(
        state: ProjectState,
        tokens: Vec<String>,
    ) -> (bool, serde_json::Value, ProjectState) {
        let parsed = parse_shell_tokens(&tokens).expect("parse shared shell tokens");
        let mut engine = GentleEngine::from_state(state);
        let run = execute_shell_command_with_options(
            &mut engine,
            &parsed,
            &ShellExecutionOptions {
                allow_screenshots: false,
                allow_agent_commands: true,
            },
        )
        .expect("execute shared shell command");
        (run.state_changed, run.output, engine.state().clone())
    }

    #[test]
    fn test_parse_rebase_site_with_slash_notation() {
        let (seq, cut, overlap) = parse_rebase_site("GAATTC (1/5)").unwrap();
        assert_eq!(seq, "GAATTC");
        assert_eq!(cut, 1);
        assert_eq!(overlap, 4);
    }

    #[test]
    fn test_parse_rebase_withrefm_minimal() {
        let text = demo_rebase_withrefm_text();
        let items = parse_rebase_withrefm(text, true);
        assert_eq!(items.len(), 1);
        assert_eq!(items[0].name, "EcoRI");
        assert_eq!(items[0].sequence, "GAATTC");
        assert_eq!(items[0].overlap, 4);
    }

    #[test]
    fn test_parse_jaspar_motifs_consensus() {
        let text = demo_jaspar_pfm_text();
        let motifs = parse_jaspar_motifs(text).unwrap();
        assert_eq!(motifs.len(), 1);
        assert_eq!(motifs[0].id, "MA0001.1");
        assert_eq!(motifs[0].consensus_iupac, "ACGT");
    }

    #[test]
    fn test_load_json_arg_reads_existing_file_without_at_prefix() {
        let dir = tempdir().expect("temp dir");
        let path = write_demo_workflow_json(dir.path(), "workflow.json", "from_file");
        let loaded = load_json_arg(path.to_str().expect("utf-8 path")).expect("load json");
        assert_eq!(loaded, r#"{"run_id":"from_file","ops":[]}"#);
    }

    #[test]
    fn test_load_json_arg_strips_shebang_from_file() {
        let dir = tempdir().expect("temp dir");
        let path = write_demo_workflow_with_shebang(dir.path(), "workflow.gsh", "with_shebang");
        let loaded = load_json_arg(path.to_str().expect("utf-8 path")).expect("load json");
        assert_eq!(loaded.trim(), r#"{"run_id":"with_shebang","ops":[]}"#);
    }

    #[test]
    fn test_load_json_arg_keeps_inline_json_payload() {
        let payload = r#"{"run_id":"inline","ops":[]}"#;
        let loaded = load_json_arg(payload).expect("load inline json");
        assert_eq!(loaded, payload);
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
    fn test_load_state_returns_default_for_empty_or_whitespace_file() {
        let td = tempdir().expect("tempdir");
        let empty_path = td.path().join("empty_state.json");
        fs::write(&empty_path, "").expect("write empty state");
        let empty_state = load_state(empty_path.to_string_lossy().as_ref())
            .expect("empty file should initialize default state");
        assert!(empty_state.sequences.is_empty());
        assert!(empty_state.metadata.is_empty());
        assert!(empty_state.container_state.containers.is_empty());
        assert!(empty_state.container_state.arrangements.is_empty());

        let whitespace_path = td.path().join("whitespace_state.json");
        fs::write(&whitespace_path, " \n\t\r ").expect("write whitespace state");
        let whitespace_state = load_state(whitespace_path.to_string_lossy().as_ref())
            .expect("whitespace file should initialize default state");
        assert!(whitespace_state.sequences.is_empty());
        assert!(whitespace_state.metadata.is_empty());
        assert!(whitespace_state.container_state.containers.is_empty());
        assert!(whitespace_state.container_state.arrangements.is_empty());
    }

    #[test]
    fn test_shell_forwarded_command_allowlist_contains_shared_shell_commands() {
        for command in SHELL_FORWARDED_COMMANDS {
            assert!(
                is_shell_forwarded_command(command),
                "expected '{command}' to be shell-forwarded"
            );
        }
        for command in ["capabilities", "op", "workflow", "state-summary"] {
            assert!(
                !is_shell_forwarded_command(command),
                "expected '{command}' to be handled by direct CLI branch"
            );
        }
    }

    #[test]
    fn test_shell_forwarded_commands_are_documented_in_glossary() {
        let glossary: serde_json::Value =
            serde_json::from_str(include_str!("../../docs/glossary.json"))
                .expect("parse docs/glossary.json");
        let commands = glossary
            .get("commands")
            .and_then(|value| value.as_array())
            .expect("glossary commands array");
        let declared_count = glossary
            .get("command_count")
            .and_then(|value| value.as_u64())
            .expect("glossary command_count") as usize;
        assert_eq!(
            declared_count,
            commands.len(),
            "docs/glossary.json command_count does not match commands length"
        );
        let paths: Vec<&str> = commands
            .iter()
            .filter_map(|entry| entry.get("path").and_then(|value| value.as_str()))
            .collect();
        for command in SHELL_FORWARDED_COMMANDS {
            let expected_prefix = format!("{command} ");
            assert!(
                paths
                    .iter()
                    .any(|path| *path == *command || path.starts_with(&expected_prefix)),
                "shell-forwarded command '{command}' is missing from docs/glossary.json"
            );
        }
    }

    #[test]
    fn test_forwarded_commands_use_shared_shell_parser_shapes() {
        let resources = parse_shell_tokens(&[
            "resources".to_string(),
            "sync-jaspar".to_string(),
            "motifs.pfm".to_string(),
            "out.json".to_string(),
        ])
        .expect("parse resources sync-jaspar");
        assert!(matches!(
            resources,
            ShellCommand::ResourcesSyncJaspar { .. }
        ));

        let import_pool = parse_shell_tokens(&[
            "import-pool".to_string(),
            "demo.pool.gentle.json".to_string(),
            "pref".to_string(),
        ])
        .expect("parse import-pool");
        assert!(matches!(import_pool, ShellCommand::ImportPool { .. }));

        let extend = parse_shell_tokens(&[
            "genomes".to_string(),
            "extend-anchor".to_string(),
            "anch".to_string(),
            "5p".to_string(),
            "25".to_string(),
        ])
        .expect("parse genomes extend-anchor");
        assert!(matches!(extend, ShellCommand::ReferenceExtendAnchor { .. }));

        let verify = parse_shell_tokens(&[
            "genomes".to_string(),
            "verify-anchor".to_string(),
            "anch".to_string(),
        ])
        .expect("parse genomes verify-anchor");
        assert!(matches!(verify, ShellCommand::ReferenceVerifyAnchor { .. }));

        let agents = parse_shell_tokens(&[
            "agents".to_string(),
            "ask".to_string(),
            "builtin_echo".to_string(),
            "--prompt".to_string(),
            "hello".to_string(),
            "--allow-auto-exec".to_string(),
        ])
        .expect("parse agents ask");
        assert!(matches!(agents, ShellCommand::AgentsAsk { .. }));

        let macros_template_run = parse_shell_tokens(&[
            "macros".to_string(),
            "template-run".to_string(),
            "clone".to_string(),
            "--bind".to_string(),
            "seq_id=seqA".to_string(),
        ])
        .expect("parse macros template-run");
        assert!(matches!(
            macros_template_run,
            ShellCommand::MacrosTemplateRun { .. }
        ));

        let macros_template_import = parse_shell_tokens(&[
            "macros".to_string(),
            "template-import".to_string(),
            "assets/cloning_patterns.json".to_string(),
        ])
        .expect("parse macros template-import");
        assert!(matches!(
            macros_template_import,
            ShellCommand::MacrosTemplateImport { .. }
        ));

        let routines_list = parse_shell_tokens(&[
            "routines".to_string(),
            "list".to_string(),
            "--family".to_string(),
            "crispr".to_string(),
            "--query".to_string(),
            "scan".to_string(),
        ])
        .expect("parse routines list");
        assert!(matches!(routines_list, ShellCommand::RoutinesList { .. }));

        let routines_explain = parse_shell_tokens(&[
            "routines".to_string(),
            "explain".to_string(),
            "golden_gate.type_iis_single_insert".to_string(),
        ])
        .expect("parse routines explain");
        assert!(matches!(
            routines_explain,
            ShellCommand::RoutinesExplain { .. }
        ));

        let routines_compare = parse_shell_tokens(&[
            "routines".to_string(),
            "compare".to_string(),
            "golden_gate.type_iis_single_insert".to_string(),
            "gibson.two_fragment_overlap_preview".to_string(),
        ])
        .expect("parse routines compare");
        assert!(matches!(
            routines_compare,
            ShellCommand::RoutinesCompare { .. }
        ));

        let ui_open = parse_shell_tokens(&[
            "ui".to_string(),
            "open".to_string(),
            "prepared-references".to_string(),
        ])
        .expect("parse ui open");
        assert!(matches!(ui_open, ShellCommand::UiIntent { .. }));

        let panels_import = parse_shell_tokens(&[
            "panels".to_string(),
            "import-isoform".to_string(),
            "seqA".to_string(),
            "assets/panels/tp53_isoforms_v1.json".to_string(),
            "--panel-id".to_string(),
            "tp53_isoforms_v1".to_string(),
        ])
        .expect("parse panels import-isoform");
        assert!(matches!(
            panels_import,
            ShellCommand::PanelsImportIsoform { .. }
        ));

        let panels_validate = parse_shell_tokens(&[
            "panels".to_string(),
            "validate-isoform".to_string(),
            "assets/panels/tp53_isoforms_v1.json".to_string(),
        ])
        .expect("parse panels validate-isoform");
        assert!(matches!(
            panels_validate,
            ShellCommand::PanelsValidateIsoform { .. }
        ));

        let primers_design_qpcr = parse_shell_tokens(&[
            "primers".to_string(),
            "design-qpcr".to_string(),
            "{\"DesignQpcrAssays\":{}}".to_string(),
            "--backend".to_string(),
            "internal".to_string(),
            "--primer3-exec".to_string(),
            "primer3_core".to_string(),
        ])
        .expect("parse primers design-qpcr");
        assert!(matches!(
            primers_design_qpcr,
            ShellCommand::PrimersDesignQpcr {
                backend: Some(_),
                primer3_executable: Some(_),
                ..
            }
        ));

        let primers_preflight = parse_shell_tokens(&[
            "primers".to_string(),
            "preflight".to_string(),
            "--backend".to_string(),
            "primer3".to_string(),
            "--primer3-exec".to_string(),
            "/opt/primer3/primer3_core".to_string(),
        ])
        .expect("parse primers preflight");
        assert!(matches!(
            primers_preflight,
            ShellCommand::PrimersPreflight {
                backend: Some(_),
                primer3_executable: Some(_),
            }
        ));

        let primers_seed_feature = parse_shell_tokens(&[
            "primers".to_string(),
            "seed-from-feature".to_string(),
            "seqA".to_string(),
            "5".to_string(),
        ])
        .expect("parse primers seed-from-feature");
        assert!(matches!(
            primers_seed_feature,
            ShellCommand::PrimersSeedFromFeature { .. }
        ));

        let primers_seed_splicing = parse_shell_tokens(&[
            "primers".to_string(),
            "seed-from-splicing".to_string(),
            "seqA".to_string(),
            "7".to_string(),
        ])
        .expect("parse primers seed-from-splicing");
        assert!(matches!(
            primers_seed_splicing,
            ShellCommand::PrimersSeedFromSplicing { .. }
        ));

        let primers_export_qpcr = parse_shell_tokens(&[
            "primers".to_string(),
            "export-qpcr-report".to_string(),
            "report_a".to_string(),
            "out.json".to_string(),
        ])
        .expect("parse primers export-qpcr-report");
        assert!(matches!(
            primers_export_qpcr,
            ShellCommand::PrimersExportQpcrReport { .. }
        ));

        let rna_reads_interpret = parse_shell_tokens(&[
            "rna-reads".to_string(),
            "interpret".to_string(),
            "seqA".to_string(),
            "7".to_string(),
            "reads.fa".to_string(),
            "--scope".to_string(),
            "all_overlapping_both_strands".to_string(),
        ])
        .expect("parse rna-reads interpret");
        assert!(matches!(
            rna_reads_interpret,
            ShellCommand::RnaReadsInterpret {
                scope: gentle::engine::SplicingScopePreset::AllOverlappingBothStrands,
                ..
            }
        ));

        let splicing_refs = parse_shell_tokens(&[
            "splicing-refs".to_string(),
            "derive".to_string(),
            "seqA".to_string(),
            "10".to_string(),
            "210".to_string(),
            "--scope".to_string(),
            "target_group_any_strand".to_string(),
        ])
        .expect("parse splicing-refs derive");
        assert!(matches!(
            splicing_refs,
            ShellCommand::SplicingRefsDerive {
                scope: gentle::engine::SplicingScopePreset::TargetGroupAnyStrand,
                ..
            }
        ));

        let align = parse_shell_tokens(&[
            "align".to_string(),
            "compute".to_string(),
            "query".to_string(),
            "target".to_string(),
            "--mode".to_string(),
            "local".to_string(),
            "--match".to_string(),
            "3".to_string(),
            "--mismatch".to_string(),
            "-4".to_string(),
        ])
        .expect("parse align compute");
        assert!(matches!(
            align,
            ShellCommand::AlignCompute {
                mode: gentle::engine::PairwiseAlignmentMode::Local,
                match_score: 3,
                mismatch_score: -4,
                ..
            }
        ));
    }

    #[test]
    fn test_parse_forwarded_shell_command_routes_resources_through_shared_parser() {
        let args = vec![
            "gentle_cli".to_string(),
            "resources".to_string(),
            "sync-rebase".to_string(),
            "input.withrefm".to_string(),
            "--bad-flag".to_string(),
        ];
        let err = parse_forwarded_shell_command(&args, 1).expect_err("expected parse failure");
        assert!(
            err.contains("Unknown option '--bad-flag' for resources sync-rebase"),
            "unexpected error: {err}"
        );
    }

    #[test]
    fn test_parse_forwarded_shell_command_routes_import_pool_through_shared_parser() {
        let args = vec![
            "gentle_cli".to_string(),
            "import-pool".to_string(),
            "demo.pool.gentle.json".to_string(),
            "pref".to_string(),
        ];
        let parsed = parse_forwarded_shell_command(&args, 1).expect("parse forwarded");
        assert!(matches!(parsed, Some(ShellCommand::ImportPool { .. })));
    }

    #[test]
    fn test_parse_forwarded_shell_command_routes_genbank_fetch() {
        let args = vec![
            "gentle_cli".to_string(),
            "genbank".to_string(),
            "fetch".to_string(),
            "AY738222".to_string(),
            "--as-id".to_string(),
            "luc_promega".to_string(),
        ];
        let parsed = parse_forwarded_shell_command(&args, 1).expect("parse forwarded");
        match parsed {
            Some(ShellCommand::GenbankFetch { accession, as_id }) => {
                assert_eq!(accession, "AY738222");
                assert_eq!(as_id.as_deref(), Some("luc_promega"));
            }
            other => panic!("unexpected parsed shell command: {other:?}"),
        }
    }

    #[test]
    fn test_parse_forwarded_shell_command_routes_uniprot_fetch() {
        let args = vec![
            "gentle_cli".to_string(),
            "uniprot".to_string(),
            "fetch".to_string(),
            "P04637".to_string(),
            "--entry-id".to_string(),
            "tp53_human".to_string(),
        ];
        let parsed = parse_forwarded_shell_command(&args, 1).expect("parse forwarded");
        match parsed {
            Some(ShellCommand::UniprotFetch { query, entry_id }) => {
                assert_eq!(query, "P04637");
                assert_eq!(entry_id.as_deref(), Some("tp53_human"));
            }
            other => panic!("unexpected parsed shell command: {other:?}"),
        }
    }

    #[test]
    fn test_forwarded_resources_sync_rebase_dispatch_matches_shared_shell_execution() {
        let td = tempdir().expect("tempdir");
        let input_path = write_demo_rebase_withrefm(td.path());
        let output_path = td.path().join("rebase.json");

        let forwarded_args = vec![
            "gentle_cli".to_string(),
            "resources".to_string(),
            "sync-rebase".to_string(),
            input_path.to_string_lossy().to_string(),
            output_path.to_string_lossy().to_string(),
            "--commercial-only".to_string(),
        ];
        let shared_tokens = forwarded_args[1..].to_vec();

        let (forwarded_changed, forwarded_output, forwarded_state) =
            execute_forwarded_like_cli(ProjectState::default(), forwarded_args);
        let (shared_changed, shared_output, shared_state) =
            execute_shared_shell_tokens(ProjectState::default(), shared_tokens);

        assert_eq!(forwarded_changed, shared_changed);
        assert_eq!(forwarded_output, shared_output);
        assert_eq!(
            forwarded_state
                .sequences
                .keys()
                .cloned()
                .collect::<std::collections::BTreeSet<_>>(),
            shared_state
                .sequences
                .keys()
                .cloned()
                .collect::<std::collections::BTreeSet<_>>()
        );
    }

    #[test]
    fn test_forwarded_resources_sync_jaspar_dispatch_matches_shared_shell_execution() {
        let td = tempdir().expect("tempdir");
        let input_path = write_demo_jaspar_pfm(td.path());
        let output_path = td.path().join("motifs.json");

        let forwarded_args = vec![
            "gentle_cli".to_string(),
            "resources".to_string(),
            "sync-jaspar".to_string(),
            input_path.to_string_lossy().to_string(),
            output_path.to_string_lossy().to_string(),
        ];
        let shared_tokens = forwarded_args[1..].to_vec();

        let (forwarded_changed, forwarded_output, forwarded_state) =
            execute_forwarded_like_cli(ProjectState::default(), forwarded_args);
        let (shared_changed, shared_output, shared_state) =
            execute_shared_shell_tokens(ProjectState::default(), shared_tokens);

        assert_eq!(forwarded_changed, shared_changed);
        assert_eq!(forwarded_output, shared_output);
        assert_eq!(
            forwarded_state
                .sequences
                .keys()
                .cloned()
                .collect::<std::collections::BTreeSet<_>>(),
            shared_state
                .sequences
                .keys()
                .cloned()
                .collect::<std::collections::BTreeSet<_>>()
        );
    }

    #[test]
    fn test_forwarded_genbank_fetch_dispatch_matches_shared_shell_execution() {
        let _env_lock = TEST_ENV_LOCK.lock().expect("env lock");
        let td = tempdir().expect("tempdir");
        let mock_dir = td.path().join("mock_genbank");
        fs::create_dir_all(&mock_dir).expect("create mock dir");
        fs::copy(
            "test_files/tp73.ncbi.gb",
            mock_dir.join("NC_000001.gbwithparts"),
        )
        .expect("copy mocked genbank payload");
        let efetch_template = format!("file://{}/{{accession}}.{{rettype}}", mock_dir.display());
        let _efetch_env = EnvVarGuard::set("GENTLE_NCBI_EFETCH_URL", &efetch_template);

        let forwarded_args = vec![
            "gentle_cli".to_string(),
            "genbank".to_string(),
            "fetch".to_string(),
            "NC_000001".to_string(),
            "--as-id".to_string(),
            "tp73_fetch".to_string(),
        ];
        let shared_tokens = forwarded_args[1..].to_vec();

        let (forwarded_changed, forwarded_output, forwarded_state) =
            execute_forwarded_like_cli(ProjectState::default(), forwarded_args);
        let (shared_changed, shared_output, shared_state) =
            execute_shared_shell_tokens(ProjectState::default(), shared_tokens);

        assert_eq!(forwarded_changed, shared_changed);
        assert_eq!(forwarded_output, shared_output);
        assert_eq!(
            forwarded_state
                .sequences
                .keys()
                .cloned()
                .collect::<std::collections::BTreeSet<_>>(),
            shared_state
                .sequences
                .keys()
                .cloned()
                .collect::<std::collections::BTreeSet<_>>()
        );
    }

    #[test]
    fn test_forwarded_import_pool_dispatch_matches_shared_shell_execution() {
        let td = tempdir().expect("tempdir");
        let pool_path = write_demo_pool_json(td.path());

        let forwarded_args = vec![
            "gentle_cli".to_string(),
            "import-pool".to_string(),
            pool_path.to_string_lossy().to_string(),
            "pref".to_string(),
        ];
        let shared_tokens = forwarded_args[1..].to_vec();

        let (forwarded_changed, forwarded_output, forwarded_state) =
            execute_forwarded_like_cli(ProjectState::default(), forwarded_args);
        let (shared_changed, shared_output, shared_state) =
            execute_shared_shell_tokens(ProjectState::default(), shared_tokens);

        assert_eq!(forwarded_changed, shared_changed);
        assert_eq!(forwarded_output, shared_output);
        assert_eq!(
            forwarded_state
                .sequences
                .keys()
                .cloned()
                .collect::<std::collections::BTreeSet<_>>(),
            shared_state
                .sequences
                .keys()
                .cloned()
                .collect::<std::collections::BTreeSet<_>>()
        );
    }

    #[test]
    fn test_forwarded_primers_list_dispatch_matches_shared_shell_execution() {
        let forwarded_args = vec![
            "gentle_cli".to_string(),
            "primers".to_string(),
            "list-reports".to_string(),
        ];
        let shared_tokens = forwarded_args[1..].to_vec();

        let (forwarded_changed, forwarded_output, forwarded_state) =
            execute_forwarded_like_cli(ProjectState::default(), forwarded_args);
        let (shared_changed, shared_output, shared_state) =
            execute_shared_shell_tokens(ProjectState::default(), shared_tokens);

        assert_eq!(forwarded_changed, shared_changed);
        assert_eq!(forwarded_output, shared_output);
        assert_eq!(
            forwarded_state
                .sequences
                .keys()
                .cloned()
                .collect::<std::collections::BTreeSet<_>>(),
            shared_state
                .sequences
                .keys()
                .cloned()
                .collect::<std::collections::BTreeSet<_>>()
        );
    }

    #[test]
    fn test_forwarded_primers_list_qpcr_dispatch_matches_shared_shell_execution() {
        let forwarded_args = vec![
            "gentle_cli".to_string(),
            "primers".to_string(),
            "list-qpcr-reports".to_string(),
        ];
        let shared_tokens = forwarded_args[1..].to_vec();

        let (forwarded_changed, forwarded_output, forwarded_state) =
            execute_forwarded_like_cli(ProjectState::default(), forwarded_args);
        let (shared_changed, shared_output, shared_state) =
            execute_shared_shell_tokens(ProjectState::default(), shared_tokens);

        assert_eq!(forwarded_changed, shared_changed);
        assert_eq!(forwarded_output, shared_output);
        assert_eq!(
            forwarded_state
                .sequences
                .keys()
                .cloned()
                .collect::<std::collections::BTreeSet<_>>(),
            shared_state
                .sequences
                .keys()
                .cloned()
                .collect::<std::collections::BTreeSet<_>>()
        );
    }

    #[test]
    fn test_forwarded_primers_preflight_dispatch_matches_shared_shell_execution() {
        let forwarded_args = vec![
            "gentle_cli".to_string(),
            "primers".to_string(),
            "preflight".to_string(),
            "--backend".to_string(),
            "primer3".to_string(),
            "--primer3-exec".to_string(),
            "__gentle_missing_primer3_for_forwarded_cli_parity_test__".to_string(),
        ];
        let shared_tokens = forwarded_args[1..].to_vec();

        let (forwarded_changed, forwarded_output, forwarded_state) =
            execute_forwarded_like_cli(ProjectState::default(), forwarded_args);
        let (shared_changed, shared_output, shared_state) =
            execute_shared_shell_tokens(ProjectState::default(), shared_tokens);

        assert_eq!(forwarded_changed, shared_changed);
        let mut forwarded_output_normalized = forwarded_output.clone();
        let mut shared_output_normalized = shared_output.clone();
        if let Some(preflight) = forwarded_output_normalized
            .get_mut("preflight")
            .and_then(|value| value.as_object_mut())
        {
            preflight.remove("probe_time_ms");
        }
        if let Some(preflight) = shared_output_normalized
            .get_mut("preflight")
            .and_then(|value| value.as_object_mut())
        {
            preflight.remove("probe_time_ms");
        }
        assert_eq!(forwarded_output_normalized, shared_output_normalized);
        assert_eq!(
            forwarded_output["preflight"]["reachable"].as_bool(),
            Some(false)
        );
        assert_eq!(
            forwarded_state
                .sequences
                .keys()
                .cloned()
                .collect::<std::collections::BTreeSet<_>>(),
            shared_state
                .sequences
                .keys()
                .cloned()
                .collect::<std::collections::BTreeSet<_>>()
        );
    }

    #[test]
    fn test_forwarded_primers_design_roundtrip_preserves_non_annealing_5prime_tail() {
        let mut state = ProjectState::default();
        state.sequences.insert(
            "tpl".to_string(),
            gentle::dna_sequence::DNAsequence::from_sequence(
                "ACGTTGCATGTCAGTACGATCGTACGTAGCTAGTCGATCGTACGATCGTAGCTAGCATCGATGCTAGCTAGTACGTAGCATCGATCGTAGCTAGCATGCTAGCTAGTCGATCGATCGTACGATCG",
            )
            .expect("sequence"),
        );

        let request = serde_json::to_string(&Operation::DesignPrimerPairs {
            template: "tpl".to_string(),
            roi_start_0based: 40,
            roi_end_0based: 80,
            forward: gentle::engine::PrimerDesignSideConstraint {
                min_length: 20,
                max_length: 20,
                location_0based: Some(5),
                start_0based: None,
                end_0based: None,
                min_tm_c: 40.0,
                max_tm_c: 90.0,
                min_gc_fraction: 0.0,
                max_gc_fraction: 1.0,
                max_anneal_hits: 1000,
                non_annealing_5prime_tail: Some("GAATTC".to_string()),
                ..Default::default()
            },
            reverse: gentle::engine::PrimerDesignSideConstraint {
                min_length: 20,
                max_length: 20,
                location_0based: Some(90),
                start_0based: None,
                end_0based: None,
                min_tm_c: 40.0,
                max_tm_c: 90.0,
                min_gc_fraction: 0.0,
                max_gc_fraction: 1.0,
                max_anneal_hits: 1000,
                non_annealing_5prime_tail: Some("CTCGAG".to_string()),
                ..Default::default()
            },
            min_amplicon_bp: 40,
            max_amplicon_bp: 150,
            pair_constraints: gentle::engine::PrimerDesignPairConstraint::default(),
            max_tm_delta_c: Some(100.0),
            max_pairs: Some(10),
            report_id: Some("tail_roundtrip".to_string()),
        })
        .expect("serialize request");

        let (design_changed, design_output, state) = execute_forwarded_like_cli(
            state,
            vec![
                "gentle_cli".to_string(),
                "primers".to_string(),
                "design".to_string(),
                request,
                "--backend".to_string(),
                "internal".to_string(),
            ],
        );
        assert!(design_changed);
        assert_eq!(
            design_output["report"]["report_id"].as_str(),
            Some("tail_roundtrip")
        );
        assert!(
            design_output["report"]["pair_count"]
                .as_u64()
                .unwrap_or_default()
                >= 1
        );

        let design_first_pair = &design_output["report"]["pairs"][0];
        assert_eq!(
            design_first_pair["forward"]["anneal_length_bp"].as_u64(),
            Some(20)
        );
        assert_eq!(
            design_first_pair["reverse"]["anneal_length_bp"].as_u64(),
            Some(20)
        );
        assert_eq!(
            design_first_pair["forward"]["non_annealing_5prime_tail_bp"].as_u64(),
            Some(6)
        );
        assert_eq!(
            design_first_pair["reverse"]["non_annealing_5prime_tail_bp"].as_u64(),
            Some(6)
        );
        assert!(
            design_first_pair["forward"]["sequence"]
                .as_str()
                .unwrap_or_default()
                .starts_with("GAATTC")
        );
        assert!(
            design_first_pair["reverse"]["sequence"]
                .as_str()
                .unwrap_or_default()
                .starts_with("CTCGAG")
        );

        let (show_changed, shown_output, _state) = execute_forwarded_like_cli(
            state,
            vec![
                "gentle_cli".to_string(),
                "primers".to_string(),
                "show-report".to_string(),
                "tail_roundtrip".to_string(),
            ],
        );
        assert!(!show_changed);
        let shown_first_pair = &shown_output["report"]["pairs"][0];
        assert_eq!(
            shown_first_pair["forward"]["non_annealing_5prime_tail_bp"].as_u64(),
            Some(6)
        );
        assert_eq!(
            shown_first_pair["reverse"]["non_annealing_5prime_tail_bp"].as_u64(),
            Some(6)
        );
        assert!(
            shown_first_pair["forward"]["sequence"]
                .as_str()
                .unwrap_or_default()
                .starts_with("GAATTC")
        );
        assert!(
            shown_first_pair["reverse"]["sequence"]
                .as_str()
                .unwrap_or_default()
                .starts_with("CTCGAG")
        );
    }

    #[test]
    fn test_parse_forwarded_shell_command_routes_macros_template_run() {
        let args = vec![
            "gentle_cli".to_string(),
            "macros".to_string(),
            "template-run".to_string(),
            "clone".to_string(),
            "--bind".to_string(),
            "seq_id=seqA".to_string(),
            "--transactional".to_string(),
        ];
        let parsed = parse_forwarded_shell_command(&args, 1).expect("parse forwarded");
        match parsed {
            Some(ShellCommand::MacrosTemplateRun {
                name,
                bindings,
                transactional,
                validate_only,
            }) => {
                assert_eq!(name, "clone");
                assert_eq!(bindings.get("seq_id").map(String::as_str), Some("seqA"));
                assert!(transactional);
                assert!(!validate_only);
            }
            other => panic!("unexpected parsed shell command: {other:?}"),
        }
    }

    #[test]
    fn test_parse_forwarded_shell_command_routes_macros_template_run_validate_only() {
        let args = vec![
            "gentle_cli".to_string(),
            "macros".to_string(),
            "template-run".to_string(),
            "clone".to_string(),
            "--validate-only".to_string(),
        ];
        let parsed = parse_forwarded_shell_command(&args, 1).expect("parse forwarded");
        match parsed {
            Some(ShellCommand::MacrosTemplateRun { validate_only, .. }) => {
                assert!(validate_only);
            }
            other => panic!("unexpected parsed shell command: {other:?}"),
        }
    }

    #[test]
    fn test_parse_forwarded_shell_command_routes_macros_template_import() {
        let args = vec![
            "gentle_cli".to_string(),
            "macros".to_string(),
            "template-import".to_string(),
            "assets/cloning_patterns.json".to_string(),
        ];
        let parsed = parse_forwarded_shell_command(&args, 1).expect("parse forwarded");
        match parsed {
            Some(ShellCommand::MacrosTemplateImport { path }) => {
                assert_eq!(path, "assets/cloning_patterns.json");
            }
            other => panic!("unexpected parsed shell command: {other:?}"),
        }
    }

    #[test]
    fn test_parse_forwarded_shell_command_routes_routines_list() {
        let args = vec![
            "gentle_cli".to_string(),
            "routines".to_string(),
            "list".to_string(),
            "--family".to_string(),
            "crispr".to_string(),
        ];
        let parsed = parse_forwarded_shell_command(&args, 1).expect("parse forwarded");
        match parsed {
            Some(ShellCommand::RoutinesList { family, .. }) => {
                assert_eq!(family.as_deref(), Some("crispr"));
            }
            other => panic!("unexpected parsed shell command: {other:?}"),
        }
    }

    #[test]
    fn test_parse_forwarded_shell_command_routes_routines_compare() {
        let args = vec![
            "gentle_cli".to_string(),
            "routines".to_string(),
            "compare".to_string(),
            "golden_gate.type_iis_single_insert".to_string(),
            "gibson.two_fragment_overlap_preview".to_string(),
        ];
        let parsed = parse_forwarded_shell_command(&args, 1).expect("parse forwarded");
        match parsed {
            Some(ShellCommand::RoutinesCompare {
                left_routine_id,
                right_routine_id,
                ..
            }) => {
                assert_eq!(left_routine_id, "golden_gate.type_iis_single_insert");
                assert_eq!(right_routine_id, "gibson.two_fragment_overlap_preview");
            }
            other => panic!("unexpected parsed shell command: {other:?}"),
        }
    }

    #[test]
    fn test_parse_forwarded_shell_command_routes_gibson_preview() {
        let args = vec![
            "gentle_cli".to_string(),
            "gibson".to_string(),
            "preview".to_string(),
            "@plan.json".to_string(),
            "--output".to_string(),
            "preview.json".to_string(),
        ];
        let parsed = parse_forwarded_shell_command(&args, 1).expect("parse forwarded");
        match parsed {
            Some(ShellCommand::GibsonPreview {
                request_json,
                output_path,
            }) => {
                assert_eq!(request_json, "@plan.json");
                assert_eq!(output_path.as_deref(), Some("preview.json"));
            }
            other => panic!("unexpected parsed shell command: {other:?}"),
        }
    }

    #[test]
    fn test_parse_forwarded_shell_command_routes_gibson_apply() {
        let args = vec![
            "gentle_cli".to_string(),
            "gibson".to_string(),
            "apply".to_string(),
            "@plan.json".to_string(),
        ];
        let parsed = parse_forwarded_shell_command(&args, 1).expect("parse forwarded");
        match parsed {
            Some(ShellCommand::GibsonApply { request_json }) => {
                assert_eq!(request_json, "@plan.json");
            }
            other => panic!("unexpected parsed shell command: {other:?}"),
        }
    }

    #[test]
    fn test_parse_forwarded_shell_command_non_forwarded_returns_none() {
        let args = vec!["gentle_cli".to_string(), "state-summary".to_string()];
        let parsed = parse_forwarded_shell_command(&args, 1).expect("parse forwarded");
        assert!(parsed.is_none());
    }
}
