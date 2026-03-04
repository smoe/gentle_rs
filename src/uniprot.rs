//! UniProt/SWISS-PROT parsing contracts shared by engine/adapters.
//!
//! This module keeps UniProt text ingestion deterministic:
//! - parse complete SWISS-PROT text (`.txt`) into structured records,
//! - preserve raw source text for offline reproducibility,
//! - expose normalized transcript/protein cross-references used by genome mapping.

use serde::{Deserialize, Serialize};
use std::collections::BTreeMap;

/// Parsed UniProt feature location with optional uncertainty flags.
#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
pub struct UniprotFeatureInterval {
    pub raw_location: String,
    pub start_aa: Option<usize>,
    pub end_aa: Option<usize>,
    pub start_uncertain: bool,
    pub end_uncertain: bool,
}

/// Parsed UniProt `FT` feature row with qualifiers preserved.
#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
pub struct UniprotFeature {
    pub key: String,
    pub interval: UniprotFeatureInterval,
    pub note: Option<String>,
    pub qualifiers: BTreeMap<String, String>,
}

/// Parsed `DR   Ensembl; ...` UniProt cross-reference.
#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
pub struct UniprotEnsemblXref {
    pub transcript_id: Option<String>,
    pub protein_id: Option<String>,
    pub gene_id: Option<String>,
    pub isoform_id: Option<String>,
}

/// Canonical parsed UniProt entry persisted in project metadata.
#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
pub struct UniprotEntry {
    pub schema: String,
    pub entry_id: String,
    pub accession: String,
    pub primary_id: String,
    pub reviewed: Option<bool>,
    pub protein_name: Option<String>,
    pub organism: Option<String>,
    pub gene_names: Vec<String>,
    pub sequence: String,
    pub sequence_length: usize,
    pub features: Vec<UniprotFeature>,
    pub ensembl_xrefs: Vec<UniprotEnsemblXref>,
    pub aliases: Vec<String>,
    pub source: String,
    pub source_query: Option<String>,
    pub imported_at_unix_ms: u128,
    pub raw_swiss_prot_text: String,
}

/// Summary row for UniProt entry listing.
#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
pub struct UniprotEntrySummary {
    pub entry_id: String,
    pub accession: String,
    pub primary_id: String,
    pub sequence_length: usize,
    pub feature_count: usize,
    pub ensembl_xref_count: usize,
    pub imported_at_unix_ms: u128,
    pub source: String,
}

/// Projected amino-acid interval onto genomic coordinates.
#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
pub struct UniprotAaGenomicSegment {
    pub aa_start: usize,
    pub aa_end: usize,
    pub genomic_start_1based: usize,
    pub genomic_end_1based: usize,
    pub strand: String,
}

/// One mapped UniProt feature against a transcript/genome projection.
#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
pub struct UniprotFeatureProjection {
    pub feature_key: String,
    pub feature_note: Option<String>,
    pub aa_start: usize,
    pub aa_end: usize,
    pub genomic_segments: Vec<UniprotAaGenomicSegment>,
}

/// UniProt-to-transcript projection payload.
#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
pub struct UniprotTranscriptProjection {
    pub transcript_id: String,
    pub transcript_feature_id: Option<usize>,
    pub strand: String,
    pub aa_segments: Vec<UniprotAaGenomicSegment>,
    pub feature_projections: Vec<UniprotFeatureProjection>,
    pub warnings: Vec<String>,
}

/// Persisted UniProt projection against one sequence context.
#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
pub struct UniprotGenomeProjection {
    pub schema: String,
    pub projection_id: String,
    pub entry_id: String,
    pub seq_id: String,
    pub created_at_unix_ms: u128,
    pub transcript_id_filter: Option<String>,
    pub transcript_projections: Vec<UniprotTranscriptProjection>,
    pub warnings: Vec<String>,
}

/// Summary row for stored UniProt genome projections.
#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
pub struct UniprotGenomeProjectionSummary {
    pub projection_id: String,
    pub entry_id: String,
    pub seq_id: String,
    pub created_at_unix_ms: u128,
    pub transcript_projection_count: usize,
}

fn normalize_organism_text(raw: &str) -> String {
    raw.trim().trim_end_matches('.').trim().to_string()
}

fn parse_numbers(raw: &str) -> Vec<usize> {
    let mut out = vec![];
    let mut current = String::new();
    for ch in raw.chars() {
        if ch.is_ascii_digit() {
            current.push(ch);
            continue;
        }
        if !current.is_empty() {
            if let Ok(value) = current.parse::<usize>() {
                out.push(value);
            }
            current.clear();
        }
    }
    if !current.is_empty()
        && let Ok(value) = current.parse::<usize>()
    {
        out.push(value);
    }
    out
}

fn parse_feature_interval(raw_location: &str) -> UniprotFeatureInterval {
    let location = raw_location.trim().to_string();
    let nums = parse_numbers(&location);
    let (start_aa, end_aa) = match nums.as_slice() {
        [] => (None, None),
        [single] => (Some(*single), Some(*single)),
        [first, second, ..] => (Some(*first), Some(*second)),
    };
    let start_uncertain = location.contains('<') || location.contains('?');
    let end_uncertain = location.contains('>') || location.contains('?');
    UniprotFeatureInterval {
        raw_location: location,
        start_aa,
        end_aa,
        start_uncertain,
        end_uncertain,
    }
}

fn parse_ensembl_xref(line: &str) -> Option<UniprotEnsemblXref> {
    if !line.starts_with("DR   Ensembl;") {
        return None;
    }
    let body = line.trim_start_matches("DR   Ensembl;").trim();
    let parts = body
        .split(';')
        .map(|value| value.trim().trim_end_matches('.').to_string())
        .filter(|value| !value.is_empty())
        .collect::<Vec<_>>();
    let transcript_id = parts.first().cloned();
    let protein_id = parts.get(1).cloned();
    let gene_id = parts.get(2).cloned();
    let mut isoform_id: Option<String> = None;
    for token in &parts {
        let Some(open) = token.find('[') else {
            continue;
        };
        let Some(close) = token[open + 1..].find(']') else {
            continue;
        };
        let candidate = token[open + 1..open + 1 + close].trim();
        if !candidate.is_empty() {
            isoform_id = Some(candidate.to_string());
            break;
        }
    }
    Some(UniprotEnsemblXref {
        transcript_id,
        protein_id,
        gene_id,
        isoform_id,
    })
}

fn strip_wrapped_quotes(raw: &str) -> String {
    let trimmed = raw.trim();
    if trimmed.len() >= 2 && trimmed.starts_with('\"') && trimmed.ends_with('\"') {
        trimmed[1..trimmed.len() - 1].to_string()
    } else {
        trimmed.to_string()
    }
}

#[derive(Debug, Default)]
struct FeatureBuilder {
    key: String,
    raw_location: String,
    qualifiers: BTreeMap<String, String>,
    pending_qualifier: Option<String>,
}

impl FeatureBuilder {
    fn commit_pending(&mut self) {
        if let Some(name) = self.pending_qualifier.take()
            && let Some(value) = self.qualifiers.get_mut(&name)
        {
            *value = strip_wrapped_quotes(value);
        }
    }

    fn into_feature(mut self) -> UniprotFeature {
        self.commit_pending();
        let interval = parse_feature_interval(&self.raw_location);
        let note = self.qualifiers.get("note").cloned();
        UniprotFeature {
            key: self.key,
            interval,
            note,
            qualifiers: self.qualifiers,
        }
    }
}

/// Normalize UniProt record IDs used as project keys.
pub fn normalize_entry_id(raw: &str) -> String {
    let mut out = String::with_capacity(raw.len());
    for ch in raw.trim().chars() {
        if ch.is_ascii_alphanumeric() || matches!(ch, '_' | '-' | '.') {
            out.push(ch.to_ascii_uppercase());
        }
    }
    out
}

/// Parse complete SWISS-PROT text into a structured UniProt entry.
pub fn parse_swiss_prot_text(
    text: &str,
    source: &str,
    source_query: Option<&str>,
    entry_id_override: Option<&str>,
) -> Result<UniprotEntry, String> {
    let mut accession: Option<String> = None;
    let mut primary_id: Option<String> = None;
    let mut reviewed: Option<bool> = None;
    let mut protein_name: Option<String> = None;
    let mut organism = String::new();
    let mut gene_names: Vec<String> = vec![];
    let mut sequence = String::new();
    let mut ensembl_xrefs: Vec<UniprotEnsemblXref> = vec![];
    let mut features: Vec<UniprotFeature> = vec![];
    let mut current_feature: Option<FeatureBuilder> = None;
    let mut in_sequence_block = false;
    let mut saw_terminator = false;

    for line in text.lines() {
        if line.starts_with("//") {
            saw_terminator = true;
            break;
        }
        if line.starts_with("SQ   ") {
            in_sequence_block = true;
            continue;
        }
        if in_sequence_block {
            for ch in line.chars() {
                if ch.is_ascii_alphabetic() {
                    sequence.push(ch.to_ascii_uppercase());
                }
            }
            continue;
        }

        if line.starts_with("ID   ") {
            let body = line.trim_start_matches("ID   ").trim();
            if let Some(token) = body.split_whitespace().next()
                && !token.trim().is_empty()
            {
                primary_id = Some(token.trim().to_string());
            }
            if body.contains("Reviewed;") {
                reviewed = Some(true);
            } else if body.contains("Unreviewed;") {
                reviewed = Some(false);
            }
            continue;
        }
        if line.starts_with("AC   ") {
            if accession.is_none() {
                let body = line.trim_start_matches("AC   ").trim();
                for token in body.split(';') {
                    let trimmed = token.trim();
                    if !trimmed.is_empty() {
                        accession = Some(trimmed.to_string());
                        break;
                    }
                }
            }
            continue;
        }
        if line.starts_with("DE   ") {
            let body = line.trim_start_matches("DE   ").trim();
            if protein_name.is_none()
                && (body.starts_with("RecName: Full=") || body.starts_with("SubName: Full="))
                && let Some((_, value)) = body.split_once('=')
            {
                let candidate = value.trim().trim_end_matches(';').trim().to_string();
                if !candidate.is_empty() {
                    protein_name = Some(candidate);
                }
            }
            continue;
        }
        if line.starts_with("OS   ") {
            let body = line.trim_start_matches("OS   ").trim();
            if !body.is_empty() {
                if !organism.is_empty() {
                    organism.push(' ');
                }
                organism.push_str(body);
            }
            continue;
        }
        if line.starts_with("GN   ") {
            let body = line.trim_start_matches("GN   ").trim();
            for segment in body.split(';') {
                let trimmed = segment.trim();
                if let Some(value) = trimmed.strip_prefix("Name=") {
                    let token = value.trim().trim_end_matches(';').trim();
                    if !token.is_empty() {
                        gene_names.push(token.to_string());
                    }
                } else if let Some(value) = trimmed.strip_prefix("Synonyms=") {
                    for synonym in value.split(',') {
                        let token = synonym.trim().trim_end_matches(';').trim();
                        if !token.is_empty() {
                            gene_names.push(token.to_string());
                        }
                    }
                }
            }
            continue;
        }
        if line.starts_with("DR   Ensembl;") {
            if let Some(xref) = parse_ensembl_xref(line) {
                ensembl_xrefs.push(xref);
            }
            continue;
        }
        if line.starts_with("FT   ") {
            let body = line.get(5..).unwrap_or_default();
            let trimmed = body.trim_end();
            let trimmed_start = trimmed.trim_start();
            if trimmed_start.is_empty() {
                continue;
            }
            if trimmed_start.starts_with('/') {
                let Some(builder) = current_feature.as_mut() else {
                    continue;
                };
                if let Some((name_raw, value_raw)) = trimmed_start[1..].split_once('=') {
                    let name = name_raw.trim().to_string();
                    if name.is_empty() {
                        continue;
                    }
                    let value = value_raw.trim().to_string();
                    builder.qualifiers.insert(name.clone(), value.clone());
                    let quoted_complete = value.starts_with('\"') && value.ends_with('\"');
                    if value.starts_with('\"') && !quoted_complete {
                        builder.pending_qualifier = Some(name);
                    } else if let Some(stored) = builder.qualifiers.get_mut(&name) {
                        *stored = strip_wrapped_quotes(stored);
                    }
                } else if let Some(name) = builder.pending_qualifier.clone()
                    && let Some(stored) = builder.qualifiers.get_mut(&name)
                {
                    if !stored.is_empty() {
                        stored.push(' ');
                    }
                    stored.push_str(trimmed_start);
                    if stored.ends_with('\"') {
                        let cleaned = strip_wrapped_quotes(stored);
                        *stored = cleaned;
                        builder.pending_qualifier = None;
                    }
                }
                continue;
            }

            if let Some(builder) = current_feature.take() {
                features.push(builder.into_feature());
            }
            let mut parts = trimmed_start.split_whitespace();
            let key = parts.next().unwrap_or_default().trim().to_string();
            let location = parts.collect::<Vec<_>>().join(" ");
            if key.is_empty() || location.trim().is_empty() {
                continue;
            }
            current_feature = Some(FeatureBuilder {
                key,
                raw_location: location.trim().to_string(),
                qualifiers: BTreeMap::new(),
                pending_qualifier: None,
            });
            continue;
        }
    }

    if let Some(builder) = current_feature.take() {
        features.push(builder.into_feature());
    }

    if !saw_terminator {
        return Err("SWISS-PROT text is incomplete: missing terminal '//' line".to_string());
    }
    let accession =
        accession.ok_or_else(|| "SWISS-PROT text is missing AC accession".to_string())?;
    let primary_id = primary_id.ok_or_else(|| "SWISS-PROT text is missing ID line".to_string())?;
    if sequence.is_empty() {
        return Err("SWISS-PROT text is missing sequence payload after SQ".to_string());
    }
    let normalized_organism = normalize_organism_text(&organism);
    let organism = (!normalized_organism.is_empty()).then_some(normalized_organism);

    gene_names.sort();
    gene_names.dedup();

    let chosen_entry_id = entry_id_override
        .and_then(|value| {
            let normalized = normalize_entry_id(value);
            (!normalized.is_empty()).then_some(normalized)
        })
        .or_else(|| {
            let normalized = normalize_entry_id(&accession);
            (!normalized.is_empty()).then_some(normalized)
        })
        .or_else(|| {
            let normalized = normalize_entry_id(&primary_id);
            (!normalized.is_empty()).then_some(normalized)
        })
        .ok_or_else(|| "Could not derive non-empty UniProt entry_id".to_string())?;

    let mut aliases = vec![];
    for raw in [
        chosen_entry_id.as_str(),
        accession.as_str(),
        primary_id.as_str(),
    ] {
        let normalized = normalize_entry_id(raw);
        if !normalized.is_empty() && !aliases.iter().any(|existing| existing == &normalized) {
            aliases.push(normalized);
        }
    }

    Ok(UniprotEntry {
        schema: "gentle.uniprot_entry.v1".to_string(),
        entry_id: chosen_entry_id,
        accession,
        primary_id,
        reviewed,
        protein_name,
        organism,
        gene_names,
        sequence_length: sequence.len(),
        sequence,
        features,
        ensembl_xrefs,
        aliases,
        source: source.to_string(),
        source_query: source_query
            .map(str::trim)
            .filter(|value| !value.is_empty())
            .map(|value| value.to_string()),
        imported_at_unix_ms: 0,
        raw_swiss_prot_text: text.to_string(),
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    const TOY_SWISS: &str = r#"ID   TEST_HUMAN              Reviewed;         12 AA.
AC   PTEST1;
DE   RecName: Full=Toy protein;
GN   Name=TP53;
OS   Homo sapiens (Human).
DR   Ensembl; ENST000001.2; ENSP000001.2; ENSG000001.
FT   DOMAIN          2..6
FT                   /note="DNA-binding"
SQ   SEQUENCE   12 AA;  1200 MW;  0000000000000000 CRC64;
     MEEPQSDPSV EP
//
"#;

    #[test]
    fn parses_swiss_prot_text_minimal() {
        let entry =
            parse_swiss_prot_text(TOY_SWISS, "unit-test", Some("PTEST1"), None).expect("parse");
        assert_eq!(entry.accession, "PTEST1");
        assert_eq!(entry.primary_id, "TEST_HUMAN");
        assert_eq!(entry.sequence, "MEEPQSDPSVEP");
        assert_eq!(entry.sequence_length, 12);
        assert_eq!(entry.features.len(), 1);
        assert_eq!(entry.features[0].key, "DOMAIN");
        assert_eq!(entry.features[0].interval.start_aa, Some(2));
        assert_eq!(entry.features[0].interval.end_aa, Some(6));
        assert_eq!(entry.features[0].note.as_deref(), Some("DNA-binding"));
        assert_eq!(entry.ensembl_xrefs.len(), 1);
        assert_eq!(
            entry.ensembl_xrefs[0].transcript_id.as_deref(),
            Some("ENST000001.2")
        );
    }

    #[test]
    fn rejects_incomplete_swiss_prot_text() {
        let err = parse_swiss_prot_text(
            "ID   BAD_HUMAN Reviewed;\nAC   P0BAD;\nSQ   SEQUENCE\n",
            "unit-test",
            None,
            None,
        )
        .expect_err("expected parse error");
        assert!(err.contains("missing terminal '//'"));
    }
}
