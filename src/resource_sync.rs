use serde::{Deserialize, Serialize};
use std::{collections::HashMap, fs, path::Path};

pub const DEFAULT_REBASE_RESOURCE_PATH: &str = "data/resources/rebase.enzymes.json";
pub const DEFAULT_JASPAR_RESOURCE_PATH: &str = "data/resources/jaspar.motifs.json";

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

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parses_rebase_edge_fixture_and_extracts_cut_geometry() {
        let text = include_str!("../test_files/data/rebase.edge.withrefm");
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
        let text = include_str!("../test_files/data/rebase.edge.withrefm");
        let enzymes = parse_rebase_withrefm(text, true);
        let names = enzymes.iter().map(|e| e.name.as_str()).collect::<Vec<_>>();
        assert_eq!(names, vec!["BsaI", "EcoRI"]);
    }

    #[test]
    fn parses_jaspar_edge_fixture_with_consensus() {
        let text = include_str!("../test_files/data/jaspar.edge.pfm");
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
}
