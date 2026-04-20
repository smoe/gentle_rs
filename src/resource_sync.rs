//! Resource synchronization (REBASE/JASPAR/ATtRACT) parsing and snapshot
//! writing.

use crate::attract_motifs::{
    ATTRACT_MOTIF_SNAPSHOT_SCHEMA, AttractMotifRecord, AttractMotifSnapshot, AttractPfmRows,
    DEFAULT_ATTRACT_RESOURCE_PATH,
};
use serde::{Deserialize, Serialize};
use std::{
    collections::{BTreeSet, HashMap},
    fs,
    path::Path,
    process::Command,
};
use tempfile::NamedTempFile;

pub const DEFAULT_REBASE_RESOURCE_PATH: &str = "data/resources/rebase.enzymes.json";
pub const DEFAULT_JASPAR_RESOURCE_PATH: &str = "data/resources/jaspar.motifs.json";
pub const DEFAULT_JASPAR_REMOTE_METADATA_PATH: &str = "data/resources/jaspar.remote_metadata.json";
pub const DEFAULT_JASPAR_BENCHMARK_PATH: &str = "data/resources/jaspar.registry_benchmark.json";

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

#[derive(Debug, Clone, Default)]
struct AttractDbHeaderMap {
    gene_name: Option<usize>,
    gene_id: Option<usize>,
    organism: Option<usize>,
    matrix_id: Option<usize>,
    motif: Option<usize>,
    len: Option<usize>,
    experiment: Option<usize>,
    family: Option<usize>,
    domain: Option<usize>,
    pubmed: Option<usize>,
    quality_score: Option<usize>,
    source_database: Option<usize>,
}

fn archive_has_named_member(archive_members: &[String], file_name: &str) -> bool {
    archive_members.iter().any(|member| {
        Path::new(member)
            .file_name()
            .and_then(|value| value.to_str())
            .map(|value| value.eq_ignore_ascii_case(file_name))
            .unwrap_or(false)
    })
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SyncReport {
    pub source: String,
    pub output: String,
    pub item_count: usize,
    pub resource: String,
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

fn read_binary_input(path_or_url: &str) -> Result<Vec<u8>, String> {
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
            .bytes()
            .map(|bytes| bytes.to_vec())
            .map_err(|e| format!("Could not read URL response '{path_or_url}': {e}"))
    } else {
        fs::read(path_or_url).map_err(|e| format!("Could not read file '{path_or_url}': {e}"))
    }
}

fn with_local_binary_input_path<T>(
    path_or_url: &str,
    mut f: impl FnMut(&Path) -> Result<T, String>,
) -> Result<T, String> {
    if path_or_url.starts_with("http://") || path_or_url.starts_with("https://") {
        let bytes = read_binary_input(path_or_url)?;
        let mut temp = NamedTempFile::new()
            .map_err(|e| format!("Could not create temporary file for '{path_or_url}': {e}"))?;
        std::io::Write::write_all(&mut temp, &bytes)
            .map_err(|e| format!("Could not write temporary file for '{path_or_url}': {e}"))?;
        f(temp.path())
    } else {
        f(Path::new(path_or_url))
    }
}

fn run_unzip(args: &[&str], archive_path: &Path, member: Option<&str>) -> Result<Vec<u8>, String> {
    let mut command = Command::new("unzip");
    for arg in args {
        command.arg(arg);
    }
    command.arg(archive_path);
    if let Some(member) = member {
        command.arg(member);
    }
    let output = command.output().map_err(|e| {
        format!(
            "Could not run unzip for archive '{}': {e}",
            archive_path.display()
        )
    })?;
    if !output.status.success() {
        let stderr = String::from_utf8_lossy(&output.stderr).trim().to_string();
        let detail = if stderr.is_empty() {
            format!("exit status {}", output.status)
        } else {
            stderr
        };
        return Err(format!(
            "Could not read ATtRACT archive '{}': {detail}",
            archive_path.display()
        ));
    }
    Ok(output.stdout)
}

fn list_zip_members(archive_path: &Path) -> Result<Vec<String>, String> {
    let bytes = run_unzip(&["-Z1"], archive_path, None)?;
    Ok(String::from_utf8_lossy(&bytes)
        .lines()
        .map(str::trim)
        .filter(|line| !line.is_empty())
        .map(|line| line.to_string())
        .collect())
}

fn read_zip_member_text(archive_path: &Path, member: &str) -> Result<String, String> {
    let bytes = run_unzip(&["-p"], archive_path, Some(member))?;
    Ok(String::from_utf8_lossy(&bytes).to_string())
}

fn normalized_column_name(raw: &str) -> String {
    raw.chars()
        .filter(|ch| ch.is_ascii_alphanumeric())
        .map(|ch| ch.to_ascii_lowercase())
        .collect()
}

fn column_value(record: &csv::StringRecord, idx: Option<usize>) -> Option<String> {
    idx.and_then(|idx| record.get(idx))
        .map(str::trim)
        .filter(|value| !value.is_empty())
        .map(|value| value.to_string())
}

fn normalize_iupac_motif(raw: &str) -> Option<String> {
    let cleaned = raw
        .chars()
        .filter(|ch| !ch.is_ascii_whitespace())
        .map(|ch| ch.to_ascii_uppercase())
        .filter(|ch| crate::iupac_code::IupacCode::is_valid_letter(*ch as u8))
        .map(|ch| if ch == 'U' { 'T' } else { ch })
        .collect::<String>();
    (!cleaned.is_empty()).then_some(cleaned)
}

fn sanitize_snapshot_token(raw: &str) -> String {
    raw.chars()
        .map(|ch| {
            if ch.is_ascii_alphanumeric() {
                ch.to_ascii_lowercase()
            } else {
                '_'
            }
        })
        .collect::<String>()
        .trim_matches('_')
        .to_string()
}

fn attract_header_map(headers: &csv::StringRecord) -> AttractDbHeaderMap {
    let mut map = AttractDbHeaderMap::default();
    for (idx, header) in headers.iter().enumerate() {
        match normalized_column_name(header).as_str() {
            "genename" | "gene" | "officialgenename" => map.gene_name = Some(idx),
            "geneid" => map.gene_id = Some(idx),
            "organism" | "species" => map.organism = Some(idx),
            "matrixid" | "motifid" => map.matrix_id = Some(idx),
            "motif" | "sequence" | "reportedsequence" => map.motif = Some(idx),
            "len" | "length" | "sequencelength" => map.len = Some(idx),
            "experiment" | "experiments" | "experimentdescription" => map.experiment = Some(idx),
            "family" => map.family = Some(idx),
            "domain" => map.domain = Some(idx),
            "pubmed" | "pubmedid" => map.pubmed = Some(idx),
            "qualityscore" | "score" | "qscore" => map.quality_score = Some(idx),
            "database" | "sourcedatabase" => map.source_database = Some(idx),
            _ => {}
        }
    }
    map
}

fn normalized_matrix_key(raw: &str) -> String {
    raw.trim().to_ascii_uppercase()
}

fn parse_pwm_numeric_values(raw: &str) -> Vec<f64> {
    raw.chars()
        .map(|ch| match ch {
            '[' | ']' | ',' | ';' | ':' | '|' => ' ',
            _ => ch,
        })
        .collect::<String>()
        .split_whitespace()
        .filter_map(|token| token.parse::<f64>().ok())
        .collect()
}

fn parse_attract_pwm_row(raw: &str) -> Option<(char, Vec<f64>)> {
    let trimmed = raw.trim();
    if trimmed.is_empty() {
        return None;
    }
    let normalized = trimmed
        .chars()
        .map(|ch| match ch {
            '[' | ']' | ',' | ';' | ':' | '|' => ' ',
            _ => ch,
        })
        .collect::<String>();
    let tokens = normalized.split_whitespace().collect::<Vec<_>>();
    if tokens.len() < 2 {
        return None;
    }
    let base = tokens[0]
        .trim_matches(|ch: char| !ch.is_ascii_alphabetic())
        .chars()
        .next()?
        .to_ascii_uppercase();
    let base = match base {
        'A' | 'C' | 'G' => base,
        'T' | 'U' => 'T',
        _ => return None,
    };
    let values = parse_pwm_numeric_values(&tokens[1..].join(" "));
    (!values.is_empty()).then_some((base, values))
}

fn parse_attract_pwm_position_row(raw: &str) -> Option<[f64; 4]> {
    let values = parse_pwm_numeric_values(raw);
    if values.len() == 4 {
        Some([values[0], values[1], values[2], values[3]])
    } else {
        None
    }
}

fn parse_attract_pwm_header_id(raw: &str) -> Option<String> {
    let trimmed = raw.trim();
    if trimmed.is_empty() {
        return None;
    }
    let without_prefix = trimmed.trim_start_matches('>');
    let token = without_prefix.split_whitespace().next()?.trim();
    (!token.is_empty()).then_some(normalized_matrix_key(token))
}

fn build_attract_pfm_rows(
    matrix_id: &str,
    rows: &HashMap<char, Vec<f64>>,
) -> Result<AttractPfmRows, String> {
    let a = rows
        .get(&'A')
        .cloned()
        .ok_or_else(|| format!("ATtRACT PWM '{matrix_id}' is missing an A row"))?;
    let c = rows
        .get(&'C')
        .cloned()
        .ok_or_else(|| format!("ATtRACT PWM '{matrix_id}' is missing a C row"))?;
    let g = rows
        .get(&'G')
        .cloned()
        .ok_or_else(|| format!("ATtRACT PWM '{matrix_id}' is missing a G row"))?;
    let t = rows
        .get(&'T')
        .cloned()
        .ok_or_else(|| format!("ATtRACT PWM '{matrix_id}' is missing a T/U row"))?;
    let len = a.len();
    if len == 0 || c.len() != len || g.len() != len || t.len() != len {
        return Err(format!(
            "ATtRACT PWM '{matrix_id}' has unequal row lengths (A={}, C={}, G={}, T={})",
            a.len(),
            c.len(),
            g.len(),
            t.len()
        ));
    }
    Ok(AttractPfmRows { a, c, g, t })
}

fn parse_attract_pwm_text(
    text: &str,
    source_label: &str,
) -> Result<(HashMap<String, AttractPfmRows>, Vec<String>), String> {
    let mut matrices = HashMap::<String, AttractPfmRows>::new();
    let mut warnings = Vec::<String>::new();
    let mut current_id: Option<String> = None;
    let mut current_rows = HashMap::<char, Vec<f64>>::new();
    let mut current_position_rows = Vec::<[f64; 4]>::new();

    let finalize_current = |matrices: &mut HashMap<String, AttractPfmRows>,
                            warnings: &mut Vec<String>,
                            current_id: &mut Option<String>,
                            current_rows: &mut HashMap<char, Vec<f64>>,
                            current_position_rows: &mut Vec<[f64; 4]>| {
        let Some(matrix_id) = current_id.take() else {
            current_rows.clear();
            current_position_rows.clear();
            return;
        };
        if current_rows.is_empty() && current_position_rows.is_empty() {
            return;
        }
        let parsed = if !current_position_rows.is_empty() {
            let len = current_position_rows.len();
            let mut a = Vec::with_capacity(len);
            let mut c = Vec::with_capacity(len);
            let mut g = Vec::with_capacity(len);
            let mut t = Vec::with_capacity(len);
            for col in current_position_rows.iter() {
                a.push(col[0]);
                c.push(col[1]);
                g.push(col[2]);
                t.push(col[3]);
            }
            Ok(AttractPfmRows { a, c, g, t })
        } else {
            build_attract_pfm_rows(&matrix_id, current_rows)
        };
        match parsed {
            Ok(pfm) => {
                matrices.entry(matrix_id.clone()).or_insert(pfm);
            }
            Err(err) => warnings.push(format!(
                "Skipped ATtRACT PWM block '{}' from '{}': {err}",
                matrix_id, source_label
            )),
        }
        current_rows.clear();
        current_position_rows.clear();
    };

    for raw in text.lines() {
        let trimmed = raw.trim();
        if trimmed.is_empty() {
            finalize_current(
                &mut matrices,
                &mut warnings,
                &mut current_id,
                &mut current_rows,
                &mut current_position_rows,
            );
            continue;
        }
        if let Some(position_row) = parse_attract_pwm_position_row(trimmed) {
            current_position_rows.push(position_row);
            continue;
        }
        if let Some((base, values)) = parse_attract_pwm_row(trimmed) {
            current_rows.insert(base, values);
            continue;
        }
        finalize_current(
            &mut matrices,
            &mut warnings,
            &mut current_id,
            &mut current_rows,
            &mut current_position_rows,
        );
        current_id = parse_attract_pwm_header_id(trimmed);
    }
    finalize_current(
        &mut matrices,
        &mut warnings,
        &mut current_id,
        &mut current_rows,
        &mut current_position_rows,
    );

    if matrices.is_empty() {
        warnings.push(format!(
            "No ATtRACT PWM matrices could be parsed from '{}'; consensus/IUPAC fallback will remain active.",
            source_label
        ));
    }
    Ok((matrices, warnings))
}

fn parse_attract_db_text(
    text: &str,
    source_label: &str,
    archive_members: &[String],
    pwm_by_matrix: &HashMap<String, AttractPfmRows>,
    pwm_warnings: &[String],
) -> Result<AttractMotifSnapshot, String> {
    let mut reader = csv::ReaderBuilder::new()
        .delimiter(b'\t')
        .flexible(true)
        .has_headers(true)
        .from_reader(text.as_bytes());
    let headers = reader
        .headers()
        .map_err(|e| format!("Could not parse ATtRACT_db.txt header from '{source_label}': {e}"))?
        .clone();
    let map = attract_header_map(&headers);
    if map.gene_name.is_none()
        || map.organism.is_none()
        || map.matrix_id.is_none()
        || map.motif.is_none()
    {
        return Err(format!(
            "ATtRACT_db.txt from '{source_label}' is missing required columns (need gene name, organism, matrix id, motif)"
        ));
    }

    let mut motifs = Vec::<AttractMotifRecord>::new();
    let mut warnings = pwm_warnings.to_vec();
    let mut seen_entry_ids = BTreeSet::<String>::new();
    let pwm_file_detected = archive_has_named_member(archive_members, "pwm.txt");
    for (row_idx, row) in reader.records().enumerate() {
        let record = row.map_err(|e| {
            format!(
                "Could not parse ATtRACT_db.txt row {} from '{}': {e}",
                row_idx + 2,
                source_label
            )
        })?;
        let gene_name = match column_value(&record, map.gene_name) {
            Some(value) => value,
            None => continue,
        };
        let organism = match column_value(&record, map.organism) {
            Some(value) => value,
            None => continue,
        };
        let matrix_id = match column_value(&record, map.matrix_id) {
            Some(value) => value,
            None => continue,
        };
        let Some(motif_iupac) =
            column_value(&record, map.motif).and_then(|value| normalize_iupac_motif(&value))
        else {
            warnings.push(format!(
                "Skipped ATtRACT row {} ('{}' / '{}') because the motif field could not be normalized to IUPAC DNA.",
                row_idx + 2,
                gene_name,
                matrix_id
            ));
            continue;
        };
        let len = column_value(&record, map.len)
            .and_then(|value| value.parse::<usize>().ok())
            .unwrap_or_else(|| motif_iupac.len());
        let gene_id = column_value(&record, map.gene_id);
        let experiment = column_value(&record, map.experiment);
        let family = column_value(&record, map.family);
        let domain = column_value(&record, map.domain);
        let pubmed_id = column_value(&record, map.pubmed);
        let quality_score = column_value(&record, map.quality_score)
            .and_then(|value| value.replace(',', ".").parse::<f64>().ok());
        let source_database = column_value(&record, map.source_database);
        let normalized_matrix_id = normalized_matrix_key(&matrix_id);
        let motif_length = motif_iupac.len();
        let candidate_pfm = pwm_by_matrix.get(&normalized_matrix_id).cloned();
        let pfm_match_status = match candidate_pfm.as_ref() {
            Some(pfm) if !pfm.a.is_empty() && pfm.a.len() == motif_length => "exact_length",
            Some(pfm) if !pfm.a.is_empty() => "matrix_id_length_mismatch",
            _ => "none",
        };
        let base_entry_id = format!(
            "{}__{}__{}",
            sanitize_snapshot_token(&gene_name),
            sanitize_snapshot_token(&organism),
            sanitize_snapshot_token(&matrix_id)
        );
        let mut entry_id = base_entry_id.clone();
        let mut duplicate = 1usize;
        while !seen_entry_ids.insert(entry_id.clone()) {
            duplicate += 1;
            entry_id = format!("{base_entry_id}__{duplicate}");
        }
        motifs.push(AttractMotifRecord {
            entry_id,
            matrix_id: matrix_id.clone(),
            gene_name: gene_name.clone(),
            gene_id,
            organism,
            motif_iupac,
            length: motif_length,
            experiment,
            family,
            domain,
            pubmed_id,
            quality_score,
            source_database,
            model_kind: if pfm_match_status == "exact_length" {
                "pwm_counts".to_string()
            } else {
                "consensus_iupac".to_string()
            },
            pwm_present: pwm_file_detected,
            pfm_match_status: pfm_match_status.to_string(),
            pfm: candidate_pfm.clone(),
        });
        if len != motif_length {
            warnings.push(format!(
                "ATtRACT row {} ('{}' / '{}') reported length {} but the normalized motif length is {}.",
                row_idx + 2,
                gene_name,
                matrix_id,
                len,
                motif_length
            ));
        }
        if let Some(candidate_pfm) = candidate_pfm {
            if candidate_pfm.a.len() != motif_length {
                warnings.push(format!(
                    "ATtRACT row {} ('{}' / '{}') keeps consensus-only matching because motif length {} does not match PWM length {}.",
                    row_idx + 2,
                    gene_name,
                    matrix_id,
                    motif_length,
                    candidate_pfm.a.len()
                ));
            }
        }
    }

    if motifs.is_empty() {
        return Err(format!(
            "No ATtRACT motifs were parsed from '{}'; check the archive contents and db columns",
            source_label
        ));
    }
    if pwm_file_detected {
        let mapped_pwm_count = motifs
            .iter()
            .filter(|row| row.pfm_match_status == "exact_length")
            .count();
        if mapped_pwm_count == 0 {
            warnings.push(
                "ATtRACT PWM members were detected in the ZIP archive, but no PWM blocks could be mapped to the normalized Matrix_id entries; consensus/IUPAC scanning remains active for all rows."
                    .to_string(),
            );
        } else if mapped_pwm_count < motifs.len() {
            warnings.push(format!(
                "ATtRACT PWM members were detected in the ZIP archive and mapped for {mapped_pwm_count}/{} normalized motif rows; remaining rows fall back to consensus/IUPAC scanning.",
                motifs.len()
            ));
        } else {
            warnings.push(format!(
                "ATtRACT PWM members were detected in the ZIP archive and mapped for all {} normalized motif rows.",
                motifs.len()
            ));
        }
    }

    Ok(AttractMotifSnapshot {
        schema: ATTRACT_MOTIF_SNAPSHOT_SCHEMA.to_string(),
        source: source_label.to_string(),
        fetched_at_unix_ms: now_unix_ms(),
        motif_count: motifs.len(),
        pwm_row_count: motifs
            .iter()
            .filter(|row| row.pfm_match_status == "exact_length")
            .count(),
        consensus_only_row_count: motifs
            .iter()
            .filter(|row| row.pfm_match_status != "exact_length")
            .count(),
        snapshot_fingerprint: None,
        archive_members: archive_members.to_vec(),
        warnings,
        motifs,
    })
}

pub fn sync_rebase(
    input: &str,
    output: Option<&str>,
    commercial_only: bool,
) -> Result<SyncReport, String> {
    let output = output.unwrap_or(DEFAULT_REBASE_RESOURCE_PATH).to_string();
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
    Ok(SyncReport {
        source: input.to_string(),
        output,
        item_count: enzymes.len(),
        resource: if commercial_only {
            "rebase-commercial".to_string()
        } else {
            "rebase".to_string()
        },
    })
}

pub fn sync_jaspar(input: &str, output: Option<&str>) -> Result<SyncReport, String> {
    let output = output.unwrap_or(DEFAULT_JASPAR_RESOURCE_PATH).to_string();
    let text = read_text_input(input)?;
    let motifs = parse_jaspar_motifs(&text)?;
    if motifs.is_empty() {
        return Err(format!("No JASPAR motifs were parsed from '{input}'"));
    }
    let snapshot = JasparMotifSnapshot {
        schema: "gentle.tf_motifs.v2".to_string(),
        source: input.to_string(),
        fetched_at_unix_ms: now_unix_ms(),
        motif_count: motifs.len(),
        motifs,
    };
    ensure_parent_dir(&output)?;
    let json = serde_json::to_string_pretty(&snapshot)
        .map_err(|e| format!("Could not serialize JASPAR resource snapshot: {e}"))?;
    fs::write(&output, json)
        .map_err(|e| format!("Could not write JASPAR output '{output}': {e}"))?;
    Ok(SyncReport {
        source: input.to_string(),
        output,
        item_count: snapshot.motif_count,
        resource: "jaspar".to_string(),
    })
}

pub fn sync_attract(input: &str, output: Option<&str>) -> Result<SyncReport, String> {
    let output = output.unwrap_or(DEFAULT_ATTRACT_RESOURCE_PATH).to_string();
    let snapshot = with_local_binary_input_path(input, |archive_path| {
        let archive_members = list_zip_members(archive_path)?;
        let db_member = archive_members
            .iter()
            .find(|member| {
                Path::new(member.as_str())
                    .file_name()
                    .and_then(|value| value.to_str())
                    .map(|value| value.eq_ignore_ascii_case("ATtRACT_db.txt"))
                    .unwrap_or(false)
            })
            .cloned()
            .ok_or_else(|| {
                format!(
                    "ATtRACT archive '{}' does not contain ATtRACT_db.txt",
                    archive_path.display()
                )
            })?;
        let db_text = read_zip_member_text(archive_path, &db_member)?;
        let (pwm_by_matrix, pwm_warnings) = if let Some(pwm_member) =
            archive_members.iter().find(|member| {
                Path::new(member.as_str())
                    .file_name()
                    .and_then(|value| value.to_str())
                    .map(|value| value.eq_ignore_ascii_case("pwm.txt"))
                    .unwrap_or(false)
            }) {
            let pwm_text = read_zip_member_text(archive_path, pwm_member)?;
            parse_attract_pwm_text(&pwm_text, input)?
        } else {
            (HashMap::new(), vec![])
        };
        parse_attract_db_text(
            &db_text,
            input,
            &archive_members,
            &pwm_by_matrix,
            &pwm_warnings,
        )
    })?;
    ensure_parent_dir(&output)?;
    let json = serde_json::to_string_pretty(&snapshot)
        .map_err(|e| format!("Could not serialize ATtRACT resource snapshot: {e}"))?;
    fs::write(&output, json)
        .map_err(|e| format!("Could not write ATtRACT output '{output}': {e}"))?;
    Ok(SyncReport {
        source: input.to_string(),
        output,
        item_count: snapshot.motif_count,
        resource: "attract".to_string(),
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use tempfile::tempdir;

    fn write_synthetic_attract_zip() -> std::path::PathBuf {
        let temp = tempdir().expect("tempdir");
        let root = temp.keep();
        let db_path = root.join("ATtRACT_db.txt");
        let pwm_path = root.join("pwm.txt");
        fs::write(
            &db_path,
            "Gene_name\tGene_id\tOrganism\tMatrix_id\tMotif\tLen\tExperiment\tDomain\tPubmed\tQuality_score\nSRSF1\tENSG00000136450\tHomo sapiens\tM001\tGAAGAA\t6\tSELEX\tRRM\t12345\t4.2\nPTBP1\tENSG00000011304\tHomo sapiens\tM002\tUCUU\t4\tCLIP\tRRM\t23456\t3.1\n",
        )
        .expect("write attract db");
        fs::write(
            &pwm_path,
            "M001\nA\t5\t0\t0\t0\t5\t5\nC\t0\t5\t0\t0\t0\t0\nG\t0\t0\t5\t0\t0\t0\nT\t0\t0\t0\t5\t0\t0\n\nM003\nA\t1\t1\t1\t1\nC\t1\t1\t1\t1\nG\t1\t1\t1\t1\nT\t1\t1\t1\t1\n",
        )
        .expect("write pwm placeholder");
        let archive_path = root.join("attract.zip");
        let output = Command::new("zip")
            .arg("-jq")
            .arg(&archive_path)
            .arg(&db_path)
            .arg(&pwm_path)
            .output()
            .expect("run zip");
        assert!(
            output.status.success(),
            "zip command failed: {}",
            String::from_utf8_lossy(&output.stderr)
        );
        archive_path
    }

    #[test]
    fn parses_rebase_edge_fixture_and_extracts_cut_geometry() {
        let text = include_str!("../test_files/fixtures/resources/rebase.edge.withrefm");
        let enzymes = parse_rebase_withrefm(text, false);
        assert_eq!(enzymes.len(), 3);

        let bsa_i = enzymes.iter().find(|e| e.name == "BsaI").expect("BsaI");
        assert_eq!(bsa_i.sequence, "GGTCTC");
        assert_eq!(bsa_i.cut, 1);
        assert_eq!(bsa_i.overlap, 4);
        assert_eq!(bsa_i.commercial, Some(true));

        let eco_ri = enzymes.iter().find(|e| e.name == "EcoRI").expect("EcoRI");
        assert_eq!(eco_ri.sequence, "GAATTC");
        assert_eq!(eco_ri.cut, 1);
        assert_eq!(eco_ri.overlap, 0);
    }

    #[test]
    fn rebase_commercial_filter_keeps_only_supplier_enzymes() {
        let text = include_str!("../test_files/fixtures/resources/rebase.edge.withrefm");
        let enzymes = parse_rebase_withrefm(text, true);
        let names = enzymes.iter().map(|e| e.name.as_str()).collect::<Vec<_>>();
        assert_eq!(names, vec!["BsaI", "EcoRI"]);
    }

    #[test]
    fn parses_jaspar_edge_fixture_with_consensus() {
        let text = include_str!("../test_files/fixtures/resources/jaspar.edge.pfm");
        let motifs = parse_jaspar_motifs(text).expect("parse motifs");
        assert_eq!(motifs.len(), 2);
        assert_eq!(motifs[0].id, "MA0001.1");
        assert_eq!(motifs[0].consensus_iupac, "ACGT");
        assert_eq!(motifs[1].id, "MA0002.1");
        assert_eq!(motifs[1].consensus_iupac, "AACC");
    }

    #[test]
    fn jaspar_rejects_unequal_rows() {
        let broken = "\
>MA9999.1 BROKEN
A [1 0 0]
C [0 1]
G [0 0 1]
T [0 0 0]
";
        let err = parse_jaspar_motifs(broken).expect_err("should fail");
        assert!(err.contains("equal length"), "unexpected error: {err}");
    }

    #[test]
    fn parses_attract_edge_archive_with_consensus_records() {
        let archive_path = write_synthetic_attract_zip();
        let snapshot =
            with_local_binary_input_path(archive_path.to_string_lossy().as_ref(), |path| {
                let members = list_zip_members(path)?;
                let db_member = members
                    .iter()
                    .find(|member| member.ends_with("ATtRACT_db.txt"))
                    .cloned()
                    .expect("db member");
                let db_text = read_zip_member_text(path, &db_member)?;
                let pwm_member = members
                    .iter()
                    .find(|member| member.ends_with("pwm.txt"))
                    .cloned()
                    .expect("pwm member");
                let pwm_text = read_zip_member_text(path, &pwm_member)?;
                let (pwm_by_matrix, pwm_warnings) =
                    parse_attract_pwm_text(&pwm_text, archive_path.to_string_lossy().as_ref())?;
                parse_attract_db_text(
                    &db_text,
                    archive_path.to_string_lossy().as_ref(),
                    &members,
                    &pwm_by_matrix,
                    &pwm_warnings,
                )
            })
            .expect("parse attract archive");
        assert_eq!(snapshot.schema, ATTRACT_MOTIF_SNAPSHOT_SCHEMA);
        assert_eq!(snapshot.motif_count, 2);
        assert_eq!(snapshot.motifs[0].gene_name, "SRSF1");
        assert_eq!(snapshot.motifs[0].motif_iupac, "GAAGAA");
        assert_eq!(snapshot.motifs[0].model_kind, "pwm_counts");
        assert!(snapshot.motifs[0].pwm_present);
        assert!(snapshot.motifs[0].pfm.is_some());
        assert_eq!(snapshot.motifs[1].motif_iupac, "TCTT");
        assert_eq!(snapshot.motifs[1].model_kind, "consensus_iupac");
        assert!(snapshot.motifs[1].pwm_present);
        assert_eq!(snapshot.motifs[1].pfm_match_status, "none");
        assert!(snapshot.motifs[1].pfm.is_none());
        assert!(
            snapshot
                .warnings
                .iter()
                .any(|warning| warning.contains("mapped for 1/2")),
            "warnings: {:?}",
            snapshot.warnings
        );
    }

    #[test]
    fn parses_attract_pwm_position_row_format() {
        let text = "\
>M001_0.6\t3
0.97\t0.01\t0.01\t0.01
0.01\t0.97\t0.01\t0.01
0.01\t0.01\t0.97\t0.01
";
        let (parsed, warnings) =
            parse_attract_pwm_text(text, "synthetic pwm").expect("parse pwm text");
        assert!(warnings.is_empty(), "warnings: {:?}", warnings);
        let pfm = parsed.get("M001_0.6").expect("matrix M001_0.6");
        assert_eq!(pfm.a, vec![0.97, 0.01, 0.01]);
        assert_eq!(pfm.c, vec![0.01, 0.97, 0.01]);
        assert_eq!(pfm.g, vec![0.01, 0.01, 0.97]);
        assert_eq!(pfm.t, vec![0.01, 0.01, 0.01]);
    }
}
