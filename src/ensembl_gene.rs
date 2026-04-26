//! Ensembl gene-entry parsing contracts shared by engine and adapters.
//!
//! This module keeps one-off Ensembl REST gene retrieval deterministic:
//! - normalize stable IDs and generated entry keys used as project metadata,
//! - parse the specific REST JSON payloads needed for one gene record
//!   (`lookup/id`, `lookup/symbol`, `sequence/id`),
//! - preserve raw response payloads for offline reproducibility.

use serde::{Deserialize, Serialize};
use std::collections::{BTreeMap, BTreeSet};

/// Summary row for one transcript connected to a fetched Ensembl gene.
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq, Default)]
#[serde(default)]
pub struct EnsemblGeneTranscriptSummary {
    pub transcript_id: String,
    pub transcript_version: Option<usize>,
    pub display_name: Option<String>,
    pub biotype: Option<String>,
    pub start_1based: Option<usize>,
    pub end_1based: Option<usize>,
    pub strand: Option<i8>,
    pub is_canonical: Option<bool>,
    pub gencode_primary: Option<bool>,
    pub translation: Option<EnsemblGeneTranslationSummary>,
    #[serde(default)]
    pub exons: Vec<EnsemblGeneExonSummary>,
}

/// Exon geometry for one transcript in an expanded Ensembl gene lookup.
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq, Default)]
#[serde(default)]
pub struct EnsemblGeneExonSummary {
    pub exon_id: String,
    pub exon_version: Option<usize>,
    pub start_1based: usize,
    pub end_1based: usize,
    pub strand: Option<i8>,
    pub seq_region_name: Option<String>,
}

/// Translation geometry for one protein-coding transcript.
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq, Default)]
#[serde(default)]
pub struct EnsemblGeneTranslationSummary {
    pub translation_id: String,
    pub translation_version: Option<usize>,
    pub length_aa: Option<usize>,
    pub genomic_start_1based: Option<usize>,
    pub genomic_end_1based: Option<usize>,
}

/// Canonical parsed Ensembl gene entry persisted in project metadata.
#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
pub struct EnsemblGeneEntry {
    pub schema: String,
    pub entry_id: String,
    pub gene_id: String,
    pub gene_version: Option<usize>,
    pub gene_symbol: Option<String>,
    pub gene_display_name: Option<String>,
    pub species: Option<String>,
    pub assembly_name: Option<String>,
    pub biotype: Option<String>,
    pub strand: Option<i8>,
    pub seq_region_name: Option<String>,
    pub genomic_start_1based: Option<usize>,
    pub genomic_end_1based: Option<usize>,
    pub sequence: String,
    pub sequence_length: usize,
    #[serde(default)]
    pub transcripts: Vec<EnsemblGeneTranscriptSummary>,
    #[serde(default)]
    pub aliases: Vec<String>,
    pub source: String,
    pub source_query: Option<String>,
    pub imported_at_unix_ms: u128,
    pub lookup_source_url: String,
    pub sequence_source_url: String,
    pub raw_lookup_json: String,
    pub raw_sequence_json: String,
}

/// Summary row for Ensembl gene-entry listing.
#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
pub struct EnsemblGeneEntrySummary {
    pub entry_id: String,
    pub gene_id: String,
    pub gene_symbol: Option<String>,
    pub species: Option<String>,
    pub seq_region_name: Option<String>,
    pub genomic_start_1based: Option<usize>,
    pub genomic_end_1based: Option<usize>,
    pub transcript_count: usize,
    pub sequence_length: usize,
    pub imported_at_unix_ms: u128,
    pub source_query: Option<String>,
}

/// Stored metadata container for fetched Ensembl gene entries.
#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
pub struct EnsemblGeneEntryStore {
    pub schema: String,
    pub updated_at_unix_ms: u128,
    pub entries: BTreeMap<String, EnsemblGeneEntry>,
}

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq)]
pub enum EnsemblGeneQueryKind {
    StableId,
    GeneSymbol,
}

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq)]
pub struct EnsemblGeneResolvedQuery {
    pub normalized_query: String,
    pub kind: EnsemblGeneQueryKind,
}

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq, Default)]
pub struct EnsemblGeneLookupRecord {
    pub gene_id: String,
    pub gene_version: Option<usize>,
    pub gene_symbol: Option<String>,
    pub gene_display_name: Option<String>,
    pub species: Option<String>,
    pub assembly_name: Option<String>,
    pub biotype: Option<String>,
    pub strand: Option<i8>,
    pub seq_region_name: Option<String>,
    pub genomic_start_1based: Option<usize>,
    pub genomic_end_1based: Option<usize>,
    #[serde(default)]
    pub transcripts: Vec<EnsemblGeneTranscriptSummary>,
}

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq)]
pub struct EnsemblGeneSequencePayload {
    pub id: String,
    pub sequence: String,
    pub description: Option<String>,
    pub molecule: Option<String>,
    pub version: Option<usize>,
}

pub fn normalize_entry_id(raw: &str) -> String {
    raw.trim()
        .replace([' ', '\t', '\n', '\r'], "_")
        .replace('/', "_")
        .to_ascii_lowercase()
}

fn looks_like_stable_gene_id(query: &str) -> bool {
    let upper = query.trim().to_ascii_uppercase();
    upper.starts_with("ENSG")
        || upper.starts_with("ENSMUSG")
        || upper.starts_with("ENSRNOG")
        || upper.starts_with("ENSDARG")
        || upper.starts_with("ENSGALG")
        || upper.starts_with("ENSCAFG")
        || upper.starts_with("ENSBTAG")
        || upper.starts_with("ENSXETG")
}

pub fn resolve_query(query: &str) -> Result<EnsemblGeneResolvedQuery, String> {
    let trimmed = query.trim();
    if trimmed.is_empty() {
        return Err("Ensembl gene query cannot be empty".to_string());
    }
    if looks_like_stable_gene_id(trimmed) {
        return Ok(EnsemblGeneResolvedQuery {
            normalized_query: trimmed.to_ascii_uppercase(),
            kind: EnsemblGeneQueryKind::StableId,
        });
    }
    Ok(EnsemblGeneResolvedQuery {
        normalized_query: trimmed.to_string(),
        kind: EnsemblGeneQueryKind::GeneSymbol,
    })
}

pub fn parse_gene_lookup_json(text: &str) -> Result<EnsemblGeneLookupRecord, String> {
    let value: serde_json::Value =
        serde_json::from_str(text).map_err(|e| format!("invalid JSON: {e}"))?;
    let gene_id = value
        .get("id")
        .and_then(|v| v.as_str())
        .unwrap_or_default()
        .trim()
        .to_string();
    if gene_id.is_empty() {
        return Err("lookup payload is missing gene id".to_string());
    }
    let transcripts = value
        .get("Transcript")
        .and_then(|v| v.as_array())
        .map(|rows| {
            rows.iter()
                .filter_map(|row| {
                    let transcript_id = row
                        .get("id")
                        .and_then(|v| v.as_str())
                        .unwrap_or_default()
                        .trim()
                        .to_string();
                    if transcript_id.is_empty() {
                        return None;
                    }
                    let translation = row.get("Translation").and_then(|translation| {
                        let translation_id = translation
                            .get("id")
                            .and_then(|v| v.as_str())
                            .unwrap_or_default()
                            .trim()
                            .to_string();
                        if translation_id.is_empty() {
                            return None;
                        }
                        Some(EnsemblGeneTranslationSummary {
                            translation_id,
                            translation_version: translation
                                .get("version")
                                .and_then(|v| v.as_u64())
                                .and_then(|v| usize::try_from(v).ok()),
                            length_aa: translation
                                .get("length")
                                .and_then(|v| v.as_u64())
                                .and_then(|v| usize::try_from(v).ok()),
                            genomic_start_1based: translation
                                .get("start")
                                .and_then(|v| v.as_u64())
                                .and_then(|v| usize::try_from(v).ok()),
                            genomic_end_1based: translation
                                .get("end")
                                .and_then(|v| v.as_u64())
                                .and_then(|v| usize::try_from(v).ok()),
                        })
                    });
                    let exons = row
                        .get("Exon")
                        .and_then(|v| v.as_array())
                        .map(|exon_rows| {
                            exon_rows
                                .iter()
                                .filter_map(|exon| {
                                    let exon_id = exon
                                        .get("id")
                                        .and_then(|v| v.as_str())
                                        .unwrap_or_default()
                                        .trim()
                                        .to_string();
                                    let start_1based = exon
                                        .get("start")
                                        .and_then(|v| v.as_u64())
                                        .and_then(|v| usize::try_from(v).ok())?;
                                    let end_1based = exon
                                        .get("end")
                                        .and_then(|v| v.as_u64())
                                        .and_then(|v| usize::try_from(v).ok())?;
                                    if exon_id.is_empty() || end_1based < start_1based {
                                        return None;
                                    }
                                    Some(EnsemblGeneExonSummary {
                                        exon_id,
                                        exon_version: exon
                                            .get("version")
                                            .and_then(|v| v.as_u64())
                                            .and_then(|v| usize::try_from(v).ok()),
                                        start_1based,
                                        end_1based,
                                        strand: exon
                                            .get("strand")
                                            .and_then(|v| v.as_i64())
                                            .and_then(|v| i8::try_from(v).ok()),
                                        seq_region_name: exon
                                            .get("seq_region_name")
                                            .and_then(|v| v.as_str())
                                            .map(|v| v.to_string()),
                                    })
                                })
                                .collect::<Vec<_>>()
                        })
                        .unwrap_or_default();
                    Some(EnsemblGeneTranscriptSummary {
                        transcript_id,
                        transcript_version: row
                            .get("version")
                            .and_then(|v| v.as_u64())
                            .and_then(|v| usize::try_from(v).ok()),
                        display_name: row
                            .get("display_name")
                            .and_then(|v| v.as_str())
                            .map(|v| v.to_string()),
                        biotype: row
                            .get("biotype")
                            .and_then(|v| v.as_str())
                            .map(|v| v.to_string()),
                        start_1based: row
                            .get("start")
                            .and_then(|v| v.as_u64())
                            .and_then(|v| usize::try_from(v).ok()),
                        end_1based: row
                            .get("end")
                            .and_then(|v| v.as_u64())
                            .and_then(|v| usize::try_from(v).ok()),
                        strand: row
                            .get("strand")
                            .and_then(|v| v.as_i64())
                            .and_then(|v| i8::try_from(v).ok()),
                        is_canonical: row
                            .get("is_canonical")
                            .and_then(|v| v.as_u64())
                            .map(|v| v != 0),
                        gencode_primary: row
                            .get("gencode_primary")
                            .and_then(|v| v.as_u64())
                            .map(|v| v != 0),
                        translation,
                        exons,
                    })
                })
                .collect::<Vec<_>>()
        })
        .unwrap_or_default();
    Ok(EnsemblGeneLookupRecord {
        gene_id,
        gene_version: value
            .get("version")
            .and_then(|v| v.as_u64())
            .and_then(|v| usize::try_from(v).ok()),
        gene_symbol: value
            .get("display_name")
            .and_then(|v| v.as_str())
            .map(|v| v.to_string()),
        gene_display_name: value
            .get("description")
            .and_then(|v| v.as_str())
            .map(|v| v.to_string())
            .or_else(|| {
                value
                    .get("display_name")
                    .and_then(|v| v.as_str())
                    .map(|v| v.to_string())
            }),
        species: value
            .get("species")
            .and_then(|v| v.as_str())
            .map(|v| v.to_string()),
        assembly_name: value
            .get("assembly_name")
            .and_then(|v| v.as_str())
            .map(|v| v.to_string()),
        biotype: value
            .get("biotype")
            .and_then(|v| v.as_str())
            .map(|v| v.to_string()),
        strand: value
            .get("strand")
            .and_then(|v| v.as_i64())
            .and_then(|v| i8::try_from(v).ok()),
        seq_region_name: value
            .get("seq_region_name")
            .and_then(|v| v.as_str())
            .map(|v| v.to_string()),
        genomic_start_1based: value
            .get("start")
            .and_then(|v| v.as_u64())
            .and_then(|v| usize::try_from(v).ok()),
        genomic_end_1based: value
            .get("end")
            .and_then(|v| v.as_u64())
            .and_then(|v| usize::try_from(v).ok()),
        transcripts,
    })
}

pub fn parse_gene_sequence_json(text: &str) -> Result<EnsemblGeneSequencePayload, String> {
    let value: serde_json::Value =
        serde_json::from_str(text).map_err(|e| format!("invalid JSON: {e}"))?;
    let id = value
        .get("id")
        .and_then(|v| v.as_str())
        .unwrap_or_default()
        .trim()
        .to_string();
    if id.is_empty() {
        return Err("sequence payload is missing gene id".to_string());
    }
    let sequence = value
        .get("seq")
        .and_then(|v| v.as_str())
        .unwrap_or_default()
        .trim()
        .to_string();
    if sequence.is_empty() {
        return Err("sequence payload is missing sequence letters".to_string());
    }
    Ok(EnsemblGeneSequencePayload {
        id,
        sequence,
        description: value
            .get("desc")
            .and_then(|v| v.as_str())
            .map(|v| v.to_string()),
        molecule: value
            .get("molecule")
            .and_then(|v| v.as_str())
            .map(|v| v.to_string()),
        version: value
            .get("version")
            .and_then(|v| v.as_u64())
            .and_then(|v| usize::try_from(v).ok()),
    })
}

pub fn build_entry_from_rest_payloads(
    source_query: &str,
    lookup_source_url: &str,
    lookup_json: &str,
    sequence_source_url: &str,
    sequence_json: &str,
    entry_id_override: Option<&str>,
) -> Result<EnsemblGeneEntry, String> {
    let lookup = parse_gene_lookup_json(lookup_json)?;
    let sequence = parse_gene_sequence_json(sequence_json)?;
    let default_entry_id = entry_id_override
        .map(str::trim)
        .filter(|value| !value.is_empty())
        .map(|value| value.to_string())
        .unwrap_or_else(|| lookup.gene_id.clone());
    let entry_id = normalize_entry_id(&default_entry_id);
    if entry_id.is_empty() {
        return Err("resolved Ensembl gene entry_id is empty".to_string());
    }
    let mut aliases = BTreeSet::new();
    for value in std::iter::once(entry_id.clone())
        .chain(std::iter::once(lookup.gene_id.clone()))
        .chain(lookup.gene_symbol.clone().into_iter())
        .chain(lookup.gene_display_name.clone().into_iter())
        .chain(
            lookup
                .transcripts
                .iter()
                .map(|row| row.transcript_id.clone()),
        )
    {
        let normalized = normalize_entry_id(&value);
        if !normalized.is_empty() {
            aliases.insert(normalized);
        }
    }
    Ok(EnsemblGeneEntry {
        schema: "gentle.ensembl_gene_entry.v1".to_string(),
        entry_id,
        gene_id: lookup.gene_id,
        gene_version: lookup.gene_version,
        gene_symbol: lookup.gene_symbol,
        gene_display_name: lookup.gene_display_name,
        species: lookup.species,
        assembly_name: lookup.assembly_name,
        biotype: lookup.biotype,
        strand: lookup.strand,
        seq_region_name: lookup.seq_region_name,
        genomic_start_1based: lookup.genomic_start_1based,
        genomic_end_1based: lookup.genomic_end_1based,
        sequence: sequence.sequence.clone(),
        sequence_length: sequence.sequence.len(),
        transcripts: lookup.transcripts,
        aliases: aliases.into_iter().collect(),
        source: "ensembl_rest".to_string(),
        source_query: Some(source_query.trim().to_string()),
        imported_at_unix_ms: 0,
        lookup_source_url: lookup_source_url.to_string(),
        sequence_source_url: sequence_source_url.to_string(),
        raw_lookup_json: lookup_json.to_string(),
        raw_sequence_json: sequence_json.to_string(),
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn resolve_query_supports_stable_ids_and_symbols() {
        assert_eq!(
            resolve_query("ensg00000141510").expect("resolve stable id"),
            EnsemblGeneResolvedQuery {
                normalized_query: "ENSG00000141510".to_string(),
                kind: EnsemblGeneQueryKind::StableId,
            }
        );
        assert_eq!(
            resolve_query("TP53").expect("resolve symbol"),
            EnsemblGeneResolvedQuery {
                normalized_query: "TP53".to_string(),
                kind: EnsemblGeneQueryKind::GeneSymbol,
            }
        );
        assert!(resolve_query("").is_err());
    }

    #[test]
    fn parse_rest_payloads_into_entry() {
        let lookup_json = r#"{
          "id":"ENSG00000141510",
          "version":19,
          "display_name":"TP53",
          "description":"tumor protein p53",
          "species":"homo_sapiens",
          "assembly_name":"GRCh38",
          "biotype":"protein_coding",
          "strand":-1,
          "seq_region_name":"17",
          "start":7668402,
          "end":7687550,
          "Transcript":[
            {
              "id":"ENST00000269305",
              "version":8,
              "display_name":"TP53-201",
              "biotype":"protein_coding",
              "start":7668402,
              "end":7687490,
              "strand":-1
            }
          ]
        }"#;
        let sequence_json = r#"{
          "id":"ENSG00000141510",
          "desc":"chromosome:GRCh38:17:7668402:7687550:-1",
          "seq":"ACGTACGTACGT",
          "molecule":"dna",
          "version":19
        }"#;
        let entry = build_entry_from_rest_payloads(
            "TP53",
            "https://rest.ensembl.org/lookup/symbol/homo_sapiens/TP53?content-type=application/json;expand=1",
            lookup_json,
            "https://rest.ensembl.org/sequence/id/ENSG00000141510?content-type=application/json;type=genomic",
            sequence_json,
            None,
        )
        .expect("build entry");
        assert_eq!(entry.entry_id, "ensg00000141510");
        assert_eq!(entry.gene_id, "ENSG00000141510");
        assert_eq!(entry.gene_symbol.as_deref(), Some("TP53"));
        assert_eq!(entry.seq_region_name.as_deref(), Some("17"));
        assert_eq!(entry.genomic_start_1based, Some(7668402));
        assert_eq!(entry.genomic_end_1based, Some(7687550));
        assert_eq!(entry.transcripts.len(), 1);
        assert_eq!(entry.transcripts[0].transcript_id, "ENST00000269305");
        assert_eq!(entry.sequence_length, 12);
        assert!(entry.aliases.contains(&"tp53".to_string()));
    }
}
