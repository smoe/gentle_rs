//! Ensembl protein-entry parsing contracts shared by engine/adapters.
//!
//! This module keeps Ensembl REST ingestion deterministic:
//! - normalize transcript/protein stable IDs used as project keys,
//! - parse the specific REST JSON payloads needed for transcript-first protein
//!   comparison (`lookup/id`, `sequence/id`, `overlap/translation`),
//! - preserve raw response payloads for offline reproducibility.

use serde::{Deserialize, Serialize};
use std::collections::BTreeMap;

/// Normalized Ensembl protein-feature record derived from REST
/// `overlap/translation`.
#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
pub struct EnsemblProteinFeature {
    pub feature_key: String,
    pub feature_type: String,
    pub description: Option<String>,
    pub start_aa: Option<usize>,
    pub end_aa: Option<usize>,
    pub interpro_id: Option<String>,
    pub qualifiers: BTreeMap<String, String>,
}

/// Canonical parsed Ensembl protein entry persisted in project metadata.
#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
pub struct EnsemblProteinEntry {
    pub schema: String,
    pub entry_id: String,
    pub protein_id: String,
    pub protein_version: Option<usize>,
    pub transcript_id: String,
    pub transcript_version: Option<usize>,
    pub gene_id: Option<String>,
    pub gene_symbol: Option<String>,
    pub transcript_display_name: Option<String>,
    pub species: Option<String>,
    #[serde(default)]
    pub transcript_exons: Vec<EnsemblTranscriptExon>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub transcript_translation: Option<EnsemblTranscriptTranslation>,
    pub sequence: String,
    pub sequence_length: usize,
    pub features: Vec<EnsemblProteinFeature>,
    pub aliases: Vec<String>,
    pub source: String,
    pub source_query: Option<String>,
    pub imported_at_unix_ms: u128,
    pub transcript_lookup_source_url: String,
    pub protein_lookup_source_url: Option<String>,
    pub sequence_source_url: String,
    pub feature_source_url: String,
    pub raw_transcript_lookup_json: String,
    pub raw_protein_lookup_json: Option<String>,
    pub raw_sequence_json: String,
    pub raw_feature_json: String,
}

/// Summary row for Ensembl protein entry listing.
#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
pub struct EnsemblProteinEntrySummary {
    pub entry_id: String,
    pub protein_id: String,
    pub transcript_id: String,
    pub gene_symbol: Option<String>,
    pub sequence_length: usize,
    pub feature_count: usize,
    pub imported_at_unix_ms: u128,
    pub source_query: Option<String>,
}

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq)]
pub enum EnsemblProteinQueryKind {
    Transcript,
    Protein,
}

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq)]
pub struct EnsemblProteinResolvedQuery {
    pub normalized_query: String,
    pub kind: EnsemblProteinQueryKind,
}

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq)]
pub struct EnsemblTranslationLookup {
    pub protein_id: String,
    pub protein_version: Option<usize>,
    pub transcript_id: String,
    pub species: Option<String>,
    pub length_aa: Option<usize>,
}

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq)]
pub struct EnsemblTranscriptLookup {
    pub transcript_id: String,
    pub transcript_version: Option<usize>,
    pub translation_id: Option<String>,
    pub translation_version: Option<usize>,
    pub translation_length_aa: Option<usize>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub translation_genomic_start_1based: Option<usize>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub translation_genomic_end_1based: Option<usize>,
    pub gene_id: Option<String>,
    pub transcript_display_name: Option<String>,
    pub species: Option<String>,
    pub strand: Option<i8>,
    pub seq_region_name: Option<String>,
    #[serde(default)]
    pub exons: Vec<EnsemblTranscriptExon>,
}

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq, Default)]
pub struct EnsemblTranscriptExon {
    pub exon_id: String,
    pub start_1based: usize,
    pub end_1based: usize,
    pub strand: i8,
    pub seq_region_name: Option<String>,
    pub version: Option<usize>,
}

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq, Default)]
pub struct EnsemblTranscriptTranslation {
    pub protein_id: String,
    pub protein_version: Option<usize>,
    pub length_aa: Option<usize>,
    pub genomic_start_1based: Option<usize>,
    pub genomic_end_1based: Option<usize>,
    pub species: Option<String>,
}

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq)]
pub struct EnsemblProteinSequencePayload {
    pub protein_id: String,
    pub sequence: String,
    pub molecule: Option<String>,
    pub version: Option<usize>,
    pub description: Option<String>,
}

#[derive(Debug, Clone, Deserialize)]
struct RawEnsemblTranslationLookup {
    id: String,
    #[serde(default)]
    version: Option<usize>,
    #[serde(rename = "Parent")]
    #[serde(default)]
    parent: Option<String>,
    #[serde(default)]
    species: Option<String>,
    #[serde(default)]
    length: Option<usize>,
}

#[derive(Debug, Clone, Deserialize)]
struct RawEnsemblTranscriptLookup {
    id: String,
    #[serde(default)]
    version: Option<usize>,
    #[serde(rename = "Parent")]
    #[serde(default)]
    parent: Option<String>,
    #[serde(default)]
    display_name: Option<String>,
    #[serde(default)]
    species: Option<String>,
    #[serde(default)]
    strand: Option<i8>,
    #[serde(default)]
    seq_region_name: Option<String>,
    #[serde(rename = "Exon")]
    #[serde(default)]
    exons: Vec<RawEnsemblTranscriptExon>,
    #[serde(rename = "Translation")]
    #[serde(default)]
    translation: Option<RawEnsemblEmbeddedTranslation>,
}

#[derive(Debug, Clone, Deserialize)]
struct RawEnsemblEmbeddedTranslation {
    id: String,
    #[serde(default)]
    version: Option<usize>,
    #[serde(default)]
    length: Option<usize>,
    #[serde(default)]
    start: Option<usize>,
    #[serde(default)]
    end: Option<usize>,
    #[serde(default)]
    _species: Option<String>,
}

#[derive(Debug, Clone, Deserialize)]
struct RawEnsemblTranscriptExon {
    id: String,
    start: usize,
    end: usize,
    #[serde(default)]
    strand: Option<i8>,
    #[serde(default)]
    version: Option<usize>,
    #[serde(default)]
    seq_region_name: Option<String>,
}

#[derive(Debug, Clone, Deserialize)]
struct RawEnsemblSequencePayload {
    id: String,
    seq: String,
    #[serde(default)]
    molecule: Option<String>,
    #[serde(default)]
    version: Option<usize>,
    #[serde(default)]
    desc: Option<String>,
}

#[derive(Debug, Clone, Deserialize)]
struct RawEnsemblProteinFeature {
    #[serde(default)]
    id: Option<String>,
    #[serde(default)]
    hseqname: Option<String>,
    #[serde(default)]
    description: Option<String>,
    #[serde(default)]
    start: Option<usize>,
    #[serde(default)]
    end: Option<usize>,
    #[serde(default)]
    interpro: Option<String>,
    #[serde(rename = "type")]
    #[serde(default)]
    feature_type: Option<String>,
    #[serde(default)]
    seq_region_name: Option<String>,
    #[serde(rename = "Parent")]
    #[serde(default)]
    parent: Option<String>,
    #[serde(default)]
    translation_id: Option<usize>,
    #[serde(default)]
    cigar_string: Option<String>,
    #[serde(default)]
    align_type: Option<String>,
    #[serde(default)]
    hit_start: Option<usize>,
    #[serde(default)]
    hit_end: Option<usize>,
}

fn normalize_optional_text(raw: Option<String>) -> Option<String> {
    raw.map(|value| value.trim().to_string())
        .filter(|value| !value.is_empty())
}

fn gene_symbol_from_display_name(display_name: &str) -> Option<String> {
    let trimmed = display_name.trim();
    if trimmed.is_empty() {
        return None;
    }
    let token = trimmed
        .split_once('-')
        .map(|(left, _)| left)
        .unwrap_or(trimmed)
        .trim();
    if token.is_empty() {
        None
    } else {
        Some(token.to_string())
    }
}

/// Normalize Ensembl transcript/protein record IDs used as project keys.
pub fn normalize_entry_id(raw: &str) -> String {
    let mut out = String::with_capacity(raw.len());
    for ch in raw.trim().chars() {
        if ch.is_ascii_alphanumeric() || matches!(ch, '_' | '-' | '.') {
            out.push(ch.to_ascii_uppercase());
        }
    }
    out
}

/// Resolve whether a user-supplied Ensembl query is transcript- or
/// protein-shaped.
pub fn resolve_query(raw: &str) -> Result<EnsemblProteinResolvedQuery, String> {
    let normalized_query = normalize_entry_id(raw);
    if normalized_query.is_empty() {
        return Err("Ensembl protein query cannot be empty".to_string());
    }
    let kind = if normalized_query.starts_with("ENST") {
        EnsemblProteinQueryKind::Transcript
    } else if normalized_query.starts_with("ENSP") {
        EnsemblProteinQueryKind::Protein
    } else {
        return Err(format!(
            "Unsupported Ensembl protein query '{}' (expected ENST... transcript or ENSP... protein stable ID)",
            raw.trim()
        ));
    };
    Ok(EnsemblProteinResolvedQuery {
        normalized_query,
        kind,
    })
}

/// Parse Ensembl `lookup/id/ENSP...` JSON into a translation lookup record.
pub fn parse_translation_lookup_json(text: &str) -> Result<EnsemblTranslationLookup, String> {
    let raw = serde_json::from_str::<RawEnsemblTranslationLookup>(text)
        .map_err(|e| format!("Could not parse Ensembl translation lookup JSON: {e}"))?;
    let protein_id = normalize_entry_id(&raw.id);
    if protein_id.is_empty() {
        return Err("Ensembl translation lookup is missing a protein stable ID".to_string());
    }
    let transcript_id = normalize_optional_text(raw.parent)
        .map(|value| normalize_entry_id(&value))
        .filter(|value| !value.is_empty())
        .ok_or_else(|| {
            format!(
                "Ensembl translation lookup for '{}' did not expose a parent transcript ID",
                protein_id
            )
        })?;
    Ok(EnsemblTranslationLookup {
        protein_id,
        protein_version: raw.version,
        transcript_id,
        species: normalize_optional_text(raw.species),
        length_aa: raw.length,
    })
}

/// Parse Ensembl `lookup/id/ENST...?expand=1` JSON into a transcript lookup
/// record.
pub fn parse_transcript_lookup_json(text: &str) -> Result<EnsemblTranscriptLookup, String> {
    let raw = serde_json::from_str::<RawEnsemblTranscriptLookup>(text)
        .map_err(|e| format!("Could not parse Ensembl transcript lookup JSON: {e}"))?;
    let transcript_id = normalize_entry_id(&raw.id);
    if transcript_id.is_empty() {
        return Err("Ensembl transcript lookup is missing a transcript stable ID".to_string());
    }
    let exons = raw
        .exons
        .into_iter()
        .filter_map(|exon| {
            let exon_id = normalize_entry_id(&exon.id);
            if exon_id.is_empty() || exon.end < exon.start || exon.start == 0 {
                return None;
            }
            Some(EnsemblTranscriptExon {
                exon_id,
                start_1based: exon.start,
                end_1based: exon.end,
                strand: exon.strand.unwrap_or(0),
                seq_region_name: normalize_optional_text(exon.seq_region_name),
                version: exon.version,
            })
        })
        .collect::<Vec<_>>();
    Ok(EnsemblTranscriptLookup {
        transcript_id,
        transcript_version: raw.version,
        translation_id: raw
            .translation
            .as_ref()
            .map(|translation| normalize_entry_id(&translation.id))
            .filter(|value| !value.is_empty()),
        translation_version: raw
            .translation
            .as_ref()
            .and_then(|translation| translation.version),
        translation_length_aa: raw
            .translation
            .as_ref()
            .and_then(|translation| translation.length),
        translation_genomic_start_1based: raw
            .translation
            .as_ref()
            .and_then(|translation| translation.start),
        translation_genomic_end_1based: raw
            .translation
            .as_ref()
            .and_then(|translation| translation.end),
        gene_id: raw
            .parent
            .map(|value| normalize_entry_id(&value))
            .filter(|value| !value.is_empty()),
        transcript_display_name: normalize_optional_text(raw.display_name),
        species: normalize_optional_text(raw.species),
        strand: raw.strand,
        seq_region_name: normalize_optional_text(raw.seq_region_name),
        exons,
    })
}

/// Parse Ensembl `sequence/id/ENSP...?type=protein` JSON into a protein-sequence
/// payload.
pub fn parse_sequence_json(text: &str) -> Result<EnsemblProteinSequencePayload, String> {
    let raw = serde_json::from_str::<RawEnsemblSequencePayload>(text)
        .map_err(|e| format!("Could not parse Ensembl protein sequence JSON: {e}"))?;
    let protein_id = normalize_entry_id(&raw.id);
    if protein_id.is_empty() {
        return Err("Ensembl sequence response is missing a protein stable ID".to_string());
    }
    let sequence = raw.seq.trim().to_ascii_uppercase();
    if sequence.is_empty() {
        return Err(format!(
            "Ensembl sequence response for '{}' contained no protein sequence",
            protein_id
        ));
    }
    Ok(EnsemblProteinSequencePayload {
        protein_id,
        sequence,
        molecule: normalize_optional_text(raw.molecule),
        version: raw.version,
        description: normalize_optional_text(raw.desc),
    })
}

/// Parse Ensembl `overlap/translation` JSON into normalized protein features.
pub fn parse_protein_features_json(text: &str) -> Result<Vec<EnsemblProteinFeature>, String> {
    let raw = serde_json::from_str::<Vec<RawEnsemblProteinFeature>>(text)
        .map_err(|e| format!("Could not parse Ensembl protein-feature JSON: {e}"))?;
    let mut features = vec![];
    for feature in raw {
        let feature_type = feature
            .feature_type
            .as_deref()
            .map(str::trim)
            .filter(|value| !value.is_empty())
            .unwrap_or("protein_feature")
            .to_string();
        let feature_key = normalize_optional_text(feature.id.clone())
            .or_else(|| normalize_optional_text(feature.hseqname.clone()))
            .or_else(|| Some(feature_type.clone()))
            .unwrap_or_else(|| "protein_feature".to_string());
        let mut qualifiers = BTreeMap::new();
        if let Some(value) = normalize_optional_text(feature.id) {
            qualifiers.insert("feature_id".to_string(), value);
        }
        if let Some(value) = normalize_optional_text(feature.hseqname) {
            qualifiers.insert("hit_name".to_string(), value);
        }
        if let Some(value) = normalize_optional_text(feature.seq_region_name) {
            qualifiers.insert("seq_region_name".to_string(), value);
        }
        if let Some(value) = normalize_optional_text(feature.parent) {
            qualifiers.insert(
                "parent_transcript_id".to_string(),
                normalize_entry_id(&value),
            );
        }
        if let Some(value) = feature.translation_id {
            qualifiers.insert("translation_internal_id".to_string(), value.to_string());
        }
        if let Some(value) = normalize_optional_text(feature.cigar_string) {
            qualifiers.insert("cigar_string".to_string(), value);
        }
        if let Some(value) = normalize_optional_text(feature.align_type) {
            qualifiers.insert("align_type".to_string(), value);
        }
        if let Some(value) = feature.hit_start {
            qualifiers.insert("hit_start".to_string(), value.to_string());
        }
        if let Some(value) = feature.hit_end {
            qualifiers.insert("hit_end".to_string(), value.to_string());
        }
        features.push(EnsemblProteinFeature {
            feature_key,
            feature_type,
            description: normalize_optional_text(feature.description),
            start_aa: feature.start,
            end_aa: feature.end,
            interpro_id: normalize_optional_text(feature.interpro),
            qualifiers,
        });
    }
    Ok(features)
}

/// Assemble one canonical Ensembl protein entry from deterministic REST payloads.
pub fn build_entry_from_rest_payloads(
    source_query: &str,
    transcript_lookup_source_url: &str,
    transcript_lookup_json: &str,
    protein_lookup_source_url: Option<&str>,
    protein_lookup_json: Option<&str>,
    sequence_source_url: &str,
    sequence_json: &str,
    feature_source_url: &str,
    feature_json: &str,
    entry_id_override: Option<&str>,
) -> Result<EnsemblProteinEntry, String> {
    let transcript_lookup = parse_transcript_lookup_json(transcript_lookup_json)?;
    let protein_lookup = protein_lookup_json
        .map(parse_translation_lookup_json)
        .transpose()?;
    let sequence = parse_sequence_json(sequence_json)?;
    let features = parse_protein_features_json(feature_json)?;

    let protein_id = transcript_lookup
        .translation_id
        .clone()
        .or_else(|| {
            protein_lookup
                .as_ref()
                .map(|lookup| lookup.protein_id.clone())
        })
        .ok_or_else(|| {
            format!(
                "Ensembl transcript '{}' did not expose a translation/protein stable ID",
                transcript_lookup.transcript_id
            )
        })?;
    if sequence.protein_id != protein_id {
        return Err(format!(
            "Ensembl sequence payload id '{}' does not match resolved protein '{}'",
            sequence.protein_id, protein_id
        ));
    }
    let entry_id = entry_id_override
        .map(normalize_entry_id)
        .filter(|value| !value.is_empty())
        .unwrap_or_else(|| protein_id.clone());
    let gene_symbol = transcript_lookup
        .transcript_display_name
        .as_deref()
        .and_then(gene_symbol_from_display_name);
    let mut aliases = std::collections::BTreeSet::new();
    for candidate in [
        Some(entry_id.clone()),
        Some(protein_id.clone()),
        Some(transcript_lookup.transcript_id.clone()),
        transcript_lookup.gene_id.clone(),
        gene_symbol.clone(),
        transcript_lookup.transcript_display_name.clone(),
    ]
    .into_iter()
    .flatten()
    {
        let normalized = normalize_entry_id(&candidate);
        if !normalized.is_empty() {
            aliases.insert(normalized);
        }
    }
    let transcript_translation =
        transcript_lookup
            .translation_id
            .as_ref()
            .map(|protein_id| EnsemblTranscriptTranslation {
                protein_id: protein_id.clone(),
                protein_version: transcript_lookup.translation_version,
                length_aa: transcript_lookup.translation_length_aa,
                genomic_start_1based: transcript_lookup.translation_genomic_start_1based,
                genomic_end_1based: transcript_lookup.translation_genomic_end_1based,
                species: transcript_lookup.species.clone(),
            });
    let species = transcript_lookup.species.clone().or_else(|| {
        protein_lookup
            .as_ref()
            .and_then(|lookup| lookup.species.clone())
    });
    Ok(EnsemblProteinEntry {
        schema: "gentle.ensembl_protein_entry.v1".to_string(),
        entry_id,
        protein_id,
        protein_version: transcript_lookup
            .translation_version
            .or_else(|| {
                protein_lookup
                    .as_ref()
                    .and_then(|lookup| lookup.protein_version)
            })
            .or(sequence.version),
        transcript_id: transcript_lookup.transcript_id,
        transcript_version: transcript_lookup.transcript_version,
        gene_id: transcript_lookup.gene_id,
        gene_symbol,
        transcript_display_name: transcript_lookup.transcript_display_name,
        species,
        transcript_exons: transcript_lookup.exons,
        transcript_translation,
        sequence: sequence.sequence.clone(),
        sequence_length: sequence.sequence.len(),
        features,
        aliases: aliases.into_iter().collect(),
        source: "ensembl_rest".to_string(),
        source_query: Some(source_query.trim().to_string()),
        imported_at_unix_ms: 0,
        transcript_lookup_source_url: transcript_lookup_source_url.to_string(),
        protein_lookup_source_url: protein_lookup_source_url.map(|value| value.to_string()),
        sequence_source_url: sequence_source_url.to_string(),
        feature_source_url: feature_source_url.to_string(),
        raw_transcript_lookup_json: transcript_lookup_json.to_string(),
        raw_protein_lookup_json: protein_lookup_json.map(|value| value.to_string()),
        raw_sequence_json: sequence_json.to_string(),
        raw_feature_json: feature_json.to_string(),
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn resolve_query_requires_enst_or_ensp() {
        assert_eq!(
            resolve_query("ensp00000288602.7").expect("resolve protein"),
            EnsemblProteinResolvedQuery {
                normalized_query: "ENSP00000288602.7".to_string(),
                kind: EnsemblProteinQueryKind::Protein,
            }
        );
        assert_eq!(
            resolve_query("enst00000288602").expect("resolve transcript"),
            EnsemblProteinResolvedQuery {
                normalized_query: "ENST00000288602".to_string(),
                kind: EnsemblProteinQueryKind::Transcript,
            }
        );
        assert!(resolve_query("BRAF").is_err());
    }

    #[test]
    fn parse_rest_payloads_into_entry() {
        let transcript_lookup_json = r#"{
          "id":"ENST00000288602",
          "version":11,
          "Parent":"ENSG00000157764",
          "display_name":"BRAF-201",
          "species":"homo_sapiens",
          "strand":-1,
          "seq_region_name":"7",
          "Translation":{
            "id":"ENSP00000288602",
            "version":7,
            "length":806,
            "start":140734597,
            "end":140924703,
            "species":"homo_sapiens"
          },
          "Exon":[
            {
              "id":"ENSE00001034876",
              "start":140734486,
              "end":140734770,
              "strand":-1,
              "version":11,
              "seq_region_name":"7"
            },
            {
              "id":"ENSE00003826734",
              "start":140924566,
              "end":140924732,
              "strand":-1,
              "version":1,
              "seq_region_name":"7"
            }
          ]
        }"#;
        let protein_lookup_json = r#"{
          "id":"ENSP00000288602",
          "version":7,
          "Parent":"ENST00000288602",
          "species":"homo_sapiens",
          "length":806
        }"#;
        let sequence_json = r#"{
          "id":"ENSP00000288602",
          "seq":"MEEPQSDPSVEPPLSQETFSDLWKLLPEN",
          "molecule":"protein",
          "version":7
        }"#;
        let feature_json = r#"[
          {
            "Parent":"ENST00000288602",
            "seq_region_name":"ENSP00000288602",
            "start":157,
            "end":225,
            "type":"Pfam",
            "hseqname":"PF02196",
            "id":"PF02196",
            "interpro":"IPR003116",
            "description":"Raf-like Ras-binding"
          }
        ]"#;

        let entry = build_entry_from_rest_payloads(
            "ENSP00000288602",
            "https://rest.ensembl.org/lookup/id/ENST00000288602?content-type=application/json;expand=1",
            transcript_lookup_json,
            Some("https://rest.ensembl.org/lookup/id/ENSP00000288602?content-type=application/json"),
            Some(protein_lookup_json),
            "https://rest.ensembl.org/sequence/id/ENSP00000288602?type=protein;content-type=application/json",
            sequence_json,
            "https://rest.ensembl.org/overlap/translation/ENSP00000288602?feature=protein_feature;content-type=application/json",
            feature_json,
            None,
        )
        .expect("build entry");

        assert_eq!(entry.entry_id, "ENSP00000288602");
        assert_eq!(entry.protein_id, "ENSP00000288602");
        assert_eq!(entry.transcript_id, "ENST00000288602");
        assert_eq!(entry.gene_id.as_deref(), Some("ENSG00000157764"));
        assert_eq!(entry.gene_symbol.as_deref(), Some("BRAF"));
        assert_eq!(entry.sequence_length, 29);
        assert_eq!(entry.features.len(), 1);
        assert_eq!(entry.transcript_exons.len(), 2);
        assert_eq!(entry.transcript_exons[0].exon_id, "ENSE00001034876");
        assert_eq!(
            entry
                .transcript_translation
                .as_ref()
                .and_then(|translation| translation.genomic_start_1based),
            Some(140734597)
        );
        assert_eq!(
            entry
                .transcript_translation
                .as_ref()
                .and_then(|translation| translation.genomic_end_1based),
            Some(140924703)
        );
        assert_eq!(entry.features[0].feature_key, "PF02196");
        assert_eq!(entry.features[0].feature_type, "Pfam");
        assert_eq!(
            entry.features[0].description.as_deref(),
            Some("Raf-like Ras-binding")
        );
        assert!(
            entry.aliases.contains(&"ENST00000288602".to_string()),
            "transcript id should be searchable as alias"
        );
    }
}
