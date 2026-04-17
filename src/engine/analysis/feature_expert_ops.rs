//! Feature-expert, splicing, and isoform helper routines used by engine
//! operations.
//!
//! This module keeps feature-oriented interpretation logic close together so
//! inspection/render/export operations reuse the same qualifier parsing rules.
//!
//! Look here for:
//! - feature qualifier normalization/parsing helpers
//! - splicing-reference derivation and isoform-path interpretation logic
//! - feature-centric summaries that feed both engine reports and expert views

use super::*;
use crate::ensembl_protein::{
    EnsemblProteinEntry, EnsemblProteinEntrySummary, build_entry_from_rest_payloads,
    normalize_entry_id as normalize_ensembl_protein_entry_id,
    resolve_query as resolve_ensembl_query,
};
use crate::uniprot::{UniprotFeature, UniprotFeatureProjection};
use crate::{AMINO_ACIDS, amino_acids::STOP_CODON};
use gentle_protocol::SplicingIntronSignal;

const DEFAULT_DBSNP_REFSNP_ENDPOINT: &str =
    "https://api.ncbi.nlm.nih.gov/variation/v0/refsnp/{refsnp_id}";
const DBSNP_REFSNP_ENV_VAR: &str = "GENTLE_NCBI_DBSNP_REFSNP_URL";
const DEFAULT_ENSEMBL_REST_ENDPOINT: &str = "https://rest.ensembl.org";
const ENSEMBL_REST_BASE_ENV_VAR: &str = "GENTLE_ENSEMBL_REST_BASE_URL";

#[derive(Debug, Clone, PartialEq, Eq)]
pub(super) struct DbsnpResolvedPlacement {
    pub rs_id: String,
    pub chromosome: String,
    pub chromosome_display: String,
    pub position_1based: usize,
    pub assembly_name: Option<String>,
    pub gene_symbols: Vec<String>,
}

#[derive(Clone)]
struct DerivedProteinExpertTranscript {
    transcript_id: String,
    transcript_label: String,
    transcript_feature_id: usize,
    is_reverse: bool,
    transcript_exons_1based: Vec<(usize, usize)>,
    genomic_cds_ranges_1based: Vec<(usize, usize)>,
    intron_ranges_1based: Vec<(usize, usize)>,
    cds_to_protein_segments: Vec<IsoformArchitectureCdsAaSegment>,
    derivation: Option<TranscriptProteinDerivation>,
}

#[derive(Clone)]
pub(crate) struct ExternalProteinOpinionAaSegment {
    pub(crate) aa_start: usize,
    pub(crate) aa_end: usize,
    pub(crate) genomic_start_1based: usize,
    pub(crate) genomic_end_1based: usize,
}

#[derive(Clone)]
pub(crate) struct ExternalProteinOpinionFeatureProjection {
    pub(crate) feature_key: String,
    pub(crate) feature_note: Option<String>,
    pub(crate) genomic_segments: Vec<ExternalProteinOpinionAaSegment>,
}

#[derive(Clone)]
pub(crate) struct ExternalProteinOpinionTranscript {
    pub(crate) transcript_id: String,
    pub(crate) transcript_feature_id: Option<usize>,
    pub(crate) strand: String,
    pub(crate) aa_segments: Vec<ExternalProteinOpinionAaSegment>,
    pub(crate) feature_projections: Vec<ExternalProteinOpinionFeatureProjection>,
}

#[derive(Clone)]
pub(crate) struct ExternalProteinOpinionViewSource {
    pub(crate) source: gentle_protocol::ProteinExternalOpinionSource,
    pub(crate) source_id: String,
    pub(crate) source_label: String,
    pub(crate) panel_source_label: String,
    pub(crate) transcript_id_filter: Option<String>,
    pub(crate) gene_symbol_hint: Option<String>,
    pub(crate) expected_length_aa: Option<usize>,
    pub(crate) feature_inventory_keys: Vec<String>,
    pub(crate) warnings: Vec<String>,
    pub(crate) transcripts: Vec<ExternalProteinOpinionTranscript>,
}

#[derive(Clone)]
enum ExternalProteinOpinionProvider<'a> {
    #[allow(dead_code)]
    Materialized(ExternalProteinOpinionViewSource),
    UniprotProjection {
        projection: &'a UniprotGenomeProjection,
        entry: &'a UniprotEntry,
    },
}

impl ExternalProteinOpinionProvider<'_> {
    fn transcript_id_filter(&self) -> Option<&str> {
        match self {
            Self::Materialized(source) => source.transcript_id_filter.as_deref(),
            Self::UniprotProjection { projection, .. } => {
                projection.transcript_id_filter.as_deref()
            }
        }
    }
}

impl GentleEngine {
    pub(super) fn normalize_dbsnp_rs_id(raw: &str) -> Result<(String, String), EngineError> {
        let trimmed = raw.trim();
        if trimmed.is_empty() {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: "dbSNP rs_id cannot be empty".to_string(),
            });
        }
        let numeric = trimmed
            .strip_prefix("rs")
            .or_else(|| trimmed.strip_prefix("RS"))
            .or_else(|| trimmed.strip_prefix("Rs"))
            .or_else(|| trimmed.strip_prefix("rS"))
            .unwrap_or(trimmed)
            .trim();
        if numeric.is_empty() || !numeric.chars().all(|ch| ch.is_ascii_digit()) {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: format!(
                    "Invalid dbSNP rs_id '{}' (expected digits or an rs-prefixed identifier such as rs9923231)",
                    raw
                ),
            });
        }
        Ok((format!("rs{numeric}"), numeric.to_string()))
    }

    pub(super) fn dbsnp_refsnp_url(refsnp_id: &str) -> String {
        let template = std::env::var(DBSNP_REFSNP_ENV_VAR)
            .ok()
            .filter(|value| !value.trim().is_empty())
            .unwrap_or_else(|| DEFAULT_DBSNP_REFSNP_ENDPOINT.to_string());
        if template.contains("{refsnp_id}") || template.contains("{rs_id}") {
            return template
                .replace("{refsnp_id}", refsnp_id)
                .replace("{rs_id}", refsnp_id);
        }
        if template.ends_with('/') {
            return format!("{template}{refsnp_id}");
        }
        template
    }

    fn json_scalar_to_string(value: &serde_json::Value) -> Option<String> {
        value
            .as_str()
            .map(|raw| raw.to_string())
            .or_else(|| value.as_u64().map(|raw| raw.to_string()))
            .or_else(|| value.as_i64().map(|raw| raw.to_string()))
    }

    fn dbsnp_assembly_family_token(raw: &str) -> Option<String> {
        let token = raw.trim();
        if token.is_empty() {
            return None;
        }
        for chunk in token.split(|c: char| !c.is_ascii_alphanumeric() && c != '.') {
            let lower = chunk.to_ascii_lowercase();
            if lower.is_empty() {
                continue;
            }
            let core = lower
                .split_once(".p")
                .map(|(prefix, _)| prefix)
                .unwrap_or(lower.as_str())
                .trim();
            if core.is_empty() {
                continue;
            }
            let has_alpha = core.chars().any(|c| c.is_ascii_alphabetic());
            let has_digit = core.chars().any(|c| c.is_ascii_digit());
            if has_alpha && has_digit {
                return Some(core.to_string());
            }
        }
        None
    }

    fn dbsnp_accession_chromosome_alias(raw: &str) -> Option<String> {
        let trimmed = raw.trim();
        if trimmed.is_empty() {
            return None;
        }
        let upper = trimmed.to_ascii_uppercase();
        let without_version = upper.split('.').next().unwrap_or(&upper);
        let suffix = without_version
            .strip_prefix("NC_")
            .or_else(|| without_version.strip_prefix("CM"))
            .unwrap_or(without_version)
            .trim_start_matches('0');
        if suffix.is_empty() || !suffix.chars().all(|ch| ch.is_ascii_digit()) {
            return None;
        }
        Some(match suffix {
            "23" => "X".to_string(),
            "24" => "Y".to_string(),
            "12920" => "MT".to_string(),
            other => other.to_string(),
        })
    }

    fn dbsnp_collect_gene_symbols(document: &serde_json::Value) -> Vec<String> {
        let mut symbols = std::collections::BTreeSet::new();
        let allele_annotations = document
            .get("primary_snapshot_data")
            .and_then(|value| value.get("allele_annotations"))
            .and_then(|value| value.as_array())
            .cloned()
            .unwrap_or_default();
        for allele_annotation in allele_annotations {
            let assembly_annotations = allele_annotation
                .get("assembly_annotation")
                .and_then(|value| value.as_array())
                .cloned()
                .unwrap_or_default();
            for assembly_annotation in assembly_annotations {
                let genes = assembly_annotation
                    .get("genes")
                    .and_then(|value| value.as_array())
                    .cloned()
                    .unwrap_or_default();
                for gene in genes {
                    for key in ["locus", "name", "symbol"] {
                        if let Some(symbol) = gene
                            .get(key)
                            .and_then(Self::json_scalar_to_string)
                            .map(|value| value.trim().to_string())
                            .filter(|value| !value.is_empty())
                        {
                            symbols.insert(symbol);
                        }
                    }
                }
            }
        }
        symbols.into_iter().collect()
    }

    pub(super) fn fetch_dbsnp_refsnp_text(
        refsnp_id: &str,
    ) -> Result<(String, String), EngineError> {
        let source_url = Self::dbsnp_refsnp_url(refsnp_id);
        let text = if let Some(path) = source_url.strip_prefix("file://") {
            std::fs::read_to_string(path).map_err(|e| EngineError {
                code: ErrorCode::Io,
                message: format!(
                    "Could not read dbSNP source file '{}' for rs{}: {e}",
                    path, refsnp_id
                ),
            })?
        } else {
            let client = reqwest::blocking::Client::builder()
                .timeout(Duration::from_secs(45))
                .build()
                .map_err(|e| EngineError {
                    code: ErrorCode::Io,
                    message: format!("Could not create dbSNP HTTP client: {e}"),
                })?;
            let response = client.get(&source_url).send().map_err(|e| EngineError {
                code: ErrorCode::Io,
                message: format!(
                    "Could not fetch dbSNP refSNP '{}' from '{}': {e}",
                    refsnp_id, source_url
                ),
            })?;
            if !response.status().is_success() {
                let status = response.status();
                let detail = if status == reqwest::StatusCode::NOT_FOUND {
                    format!(
                        "dbSNP refSNP '{}' was not found (HTTP {}) at '{}'",
                        refsnp_id, status, source_url
                    )
                } else {
                    format!(
                        "dbSNP refSNP '{}' returned HTTP status {} at '{}'",
                        refsnp_id, status, source_url
                    )
                };
                return Err(EngineError {
                    code: ErrorCode::NotFound,
                    message: detail,
                });
            }
            response.text().map_err(|e| EngineError {
                code: ErrorCode::Io,
                message: format!(
                    "Could not read dbSNP response body for rs{}: {e}",
                    refsnp_id
                ),
            })?
        };
        Ok((source_url, text))
    }

    pub(super) fn parse_dbsnp_refsnp_json(
        refsnp_id: &str,
        text: &str,
    ) -> Result<serde_json::Value, EngineError> {
        let document =
            serde_json::from_str::<serde_json::Value>(&text).map_err(|e| EngineError {
                code: ErrorCode::InvalidInput,
                message: format!(
                    "Could not parse dbSNP response JSON for rs{}: {e}",
                    refsnp_id
                ),
            })?;
        if let Some(error_text) = document
            .get("error")
            .and_then(Self::json_scalar_to_string)
            .or_else(|| {
                document
                    .get("error")
                    .and_then(|value| value.get("message"))
                    .and_then(Self::json_scalar_to_string)
            })
        {
            return Err(EngineError {
                code: ErrorCode::NotFound,
                message: format!(
                    "dbSNP refSNP '{}' returned error: {}",
                    refsnp_id, error_text
                ),
            });
        }
        Ok(document)
    }

    pub(super) fn resolve_dbsnp_primary_placement(
        document: &serde_json::Value,
        requested_rs_id: &str,
        requested_assembly_family: Option<&str>,
    ) -> Result<DbsnpResolvedPlacement, EngineError> {
        let resolved_rs_id = document
            .get("refsnp_id")
            .and_then(Self::json_scalar_to_string)
            .map(|value| format!("rs{}", value.trim().trim_start_matches("rs")))
            .unwrap_or_else(|| requested_rs_id.to_string());
        let requested_family = requested_assembly_family
            .map(str::trim)
            .filter(|value| !value.is_empty())
            .map(|value| value.to_ascii_lowercase());
        let placements = document
            .get("primary_snapshot_data")
            .and_then(|value| value.get("placements_with_allele"))
            .and_then(|value| value.as_array())
            .ok_or_else(|| EngineError {
                code: ErrorCode::NotFound,
                message: format!(
                    "dbSNP record '{}' does not include any primary genomic placements",
                    resolved_rs_id
                ),
            })?;
        let mut available_assemblies = std::collections::BTreeSet::new();
        let mut available_families = std::collections::BTreeSet::new();
        let mut best: Option<(usize, String, String, Option<String>)> = None;
        for placement in placements {
            let placement_seq_id = placement
                .get("seq_id")
                .and_then(Self::json_scalar_to_string)
                .unwrap_or_default();
            let traits = placement
                .get("placement_annot")
                .and_then(|value| value.get("seq_id_traits_by_assembly"))
                .and_then(|value| value.as_array())
                .cloned()
                .unwrap_or_default();
            let mut best_trait_score = 0usize;
            let mut matched_assembly_name: Option<String> = None;
            for assembly_trait in traits {
                let is_chromosome = assembly_trait
                    .get("is_chromosome")
                    .and_then(|value| value.as_bool())
                    .unwrap_or(false);
                if !is_chromosome {
                    continue;
                }
                let assembly_name = assembly_trait
                    .get("assembly_name")
                    .and_then(Self::json_scalar_to_string);
                if let Some(name) = assembly_name.as_ref() {
                    available_assemblies.insert(name.clone());
                    if let Some(family) = Self::dbsnp_assembly_family_token(name) {
                        available_families.insert(family);
                    }
                }
                if let Some(requested_family) = requested_family.as_deref()
                    && assembly_name
                        .as_deref()
                        .and_then(Self::dbsnp_assembly_family_token)
                        .as_deref()
                        != Some(requested_family)
                {
                    continue;
                }
                let mut score = 0usize;
                if placement
                    .get("is_ptlp")
                    .and_then(|value| value.as_bool())
                    .unwrap_or(false)
                {
                    score += 8;
                }
                if assembly_trait
                    .get("is_top_level")
                    .and_then(|value| value.as_bool())
                    .unwrap_or(false)
                {
                    score += 4;
                }
                if !assembly_trait
                    .get("is_alt")
                    .and_then(|value| value.as_bool())
                    .unwrap_or(false)
                {
                    score += 2;
                }
                if !assembly_trait
                    .get("is_patch")
                    .and_then(|value| value.as_bool())
                    .unwrap_or(false)
                {
                    score += 1;
                }
                if matched_assembly_name.is_none() || score > best_trait_score {
                    best_trait_score = score;
                    matched_assembly_name = assembly_name;
                }
            }
            if requested_family.is_some() && matched_assembly_name.is_none() {
                continue;
            }
            if requested_family.is_none()
                && matched_assembly_name.is_none()
                && placement_seq_id.is_empty()
            {
                continue;
            }
            let (_position_0based, spdi_seq_id) = placement
                .get("alleles")
                .and_then(|value| value.as_array())
                .and_then(|alleles| {
                    alleles.iter().find_map(|allele| {
                        let spdi = allele.get("allele")?.get("spdi")?;
                        let position = spdi.get("position")?.as_u64()?;
                        let seq_id = spdi.get("seq_id").and_then(Self::json_scalar_to_string);
                        Some((position as usize, seq_id.unwrap_or_default()))
                    })
                })
                .ok_or_else(|| EngineError {
                    code: ErrorCode::NotFound,
                    message: format!(
                        "dbSNP record '{}' does not include chromosome-position allele coordinates",
                        resolved_rs_id
                    ),
                })?;
            let chromosome = if spdi_seq_id.trim().is_empty() {
                placement_seq_id.clone()
            } else {
                spdi_seq_id
            };
            let overall_score = best_trait_score;
            let is_better = best
                .as_ref()
                .map(|(score, _, _, _)| overall_score > *score)
                .unwrap_or(true);
            if is_better {
                best = Some((
                    overall_score,
                    chromosome,
                    Self::dbsnp_accession_chromosome_alias(&placement_seq_id)
                        .unwrap_or_else(|| placement_seq_id.clone()),
                    matched_assembly_name,
                ));
            }
            if best_trait_score >= 15 {
                // Prefer the first strong top-level chromosome placement that matches the requested assembly.
                break;
            }
        }
        let Some((_, chromosome, chromosome_display, assembly_name)) = best else {
            if let Some(requested_family) = requested_family.as_deref()
                && !available_families.contains(requested_family)
            {
                return Self::resolve_dbsnp_primary_placement(document, requested_rs_id, None);
            }
            let available = if available_assemblies.is_empty() {
                "none reported".to_string()
            } else {
                available_assemblies
                    .into_iter()
                    .collect::<Vec<_>>()
                    .join(", ")
            };
            let detail = if let Some(requested_family) = requested_family.as_deref() {
                format!(
                    "dbSNP record '{}' has no top-level chromosome placement matching assembly family '{}' (available assemblies: {})",
                    resolved_rs_id, requested_family, available
                )
            } else {
                format!(
                    "dbSNP record '{}' has no usable top-level chromosome placement",
                    resolved_rs_id
                )
            };
            return Err(EngineError {
                code: ErrorCode::NotFound,
                message: detail,
            });
        };
        let gene_symbols = Self::dbsnp_collect_gene_symbols(document);
        let position_1based = placements
            .iter()
            .find_map(|placement| {
                let placement_seq_id = placement
                    .get("seq_id")
                    .and_then(Self::json_scalar_to_string)
                    .unwrap_or_default();
                let matches_chromosome = placement_seq_id == chromosome
                    || placement
                        .get("alleles")
                        .and_then(|value| value.as_array())
                        .map(|alleles| {
                            alleles.iter().any(|allele| {
                                allele
                                    .get("allele")
                                    .and_then(|value| value.get("spdi"))
                                    .and_then(|value| value.get("seq_id"))
                                    .and_then(Self::json_scalar_to_string)
                                    .map(|seq_id| seq_id == chromosome)
                                    .unwrap_or(false)
                            })
                        })
                        .unwrap_or(false);
                if !matches_chromosome {
                    return None;
                }
                placement
                    .get("alleles")
                    .and_then(|value| value.as_array())
                    .and_then(|alleles| {
                        alleles.iter().find_map(|allele| {
                            allele
                                .get("allele")
                                .and_then(|value| value.get("spdi"))
                                .and_then(|value| value.get("position"))
                                .and_then(|value| value.as_u64())
                                .map(|value| value as usize + 1)
                        })
                    })
            })
            .ok_or_else(|| EngineError {
                code: ErrorCode::NotFound,
                message: format!(
                    "dbSNP record '{}' did not expose a usable 1-based genomic position",
                    resolved_rs_id
                ),
            })?;
        Ok(DbsnpResolvedPlacement {
            rs_id: resolved_rs_id,
            chromosome,
            chromosome_display,
            position_1based,
            assembly_name,
            gene_symbols,
        })
    }

    pub(super) fn feature_qualifier_text(
        feature: &gb_io::seq::Feature,
        key: &str,
    ) -> Option<String> {
        feature
            .qualifier_values(key)
            .map(|value| value.split_whitespace().collect::<Vec<_>>().join(" "))
            .map(|value| value.trim().to_string())
            .find(|value| !value.is_empty())
    }

    pub(super) fn feature_qualifier_f64(feature: &gb_io::seq::Feature, key: &str) -> Option<f64> {
        feature
            .qualifier_values(key)
            .next()
            .and_then(|v| v.trim().parse::<f64>().ok())
    }

    pub(super) fn first_nonempty_feature_qualifier(
        feature: &gb_io::seq::Feature,
        keys: &[&str],
    ) -> Option<String> {
        for key in keys {
            if let Some(value) = Self::feature_qualifier_text(feature, key) {
                return Some(value);
            }
        }
        None
    }

    pub(super) fn serialize_ranges_1based(ranges: &[(usize, usize)]) -> Option<String> {
        let parts = ranges
            .iter()
            .filter_map(|(start, end)| {
                (*start > 0 && *end >= *start).then_some(format!("{start}-{end}"))
            })
            .collect::<Vec<_>>();
        (!parts.is_empty()).then_some(parts.join(","))
    }

    pub(super) fn parse_ranges_1based(raw: &str) -> Vec<(usize, usize)> {
        let mut out = vec![];
        for token in raw.split(',') {
            let trimmed = token.trim();
            if trimmed.is_empty() {
                continue;
            }
            let mut pieces = trimmed.splitn(2, '-');
            let Some(start_raw) = pieces.next() else {
                continue;
            };
            let Some(end_raw) = pieces.next() else {
                continue;
            };
            let Ok(start) = start_raw.trim().parse::<usize>() else {
                continue;
            };
            let Ok(end) = end_raw.trim().parse::<usize>() else {
                continue;
            };
            if start == 0 || end < start {
                continue;
            }
            out.push((start, end));
        }
        out.sort_unstable_by(|a, b| a.0.cmp(&b.0).then(a.1.cmp(&b.1)));
        out.dedup();
        out
    }

    pub(super) fn feature_qualifier_ranges_0based(
        feature: &gb_io::seq::Feature,
        key: &str,
    ) -> Vec<(usize, usize)> {
        Self::feature_qualifier_text(feature, key)
            .map(|raw| {
                Self::parse_ranges_1based(&raw)
                    .into_iter()
                    .map(|(start_1based, end_1based)| (start_1based - 1, end_1based))
                    .collect::<Vec<_>>()
            })
            .unwrap_or_default()
    }

    pub(super) fn is_tfbs_feature(feature: &gb_io::seq::Feature) -> bool {
        matches!(
            feature.kind.to_string().to_ascii_uppercase().as_str(),
            "TFBS" | "TF_BINDING_SITE" | "PROTEIN_BIND"
        )
    }

    pub(super) fn feature_display_label(
        feature: &gb_io::seq::Feature,
        feature_id: usize,
    ) -> String {
        let fallback = format!("{} #{}", feature.kind.to_string(), feature_id + 1);
        let label = Self::first_nonempty_feature_qualifier(
            feature,
            &[
                "label",
                "name",
                "standard_name",
                "gene",
                "protein_id",
                "product",
                "region_name",
                "bound_moiety",
            ],
        )
        .unwrap_or(fallback);
        let trimmed = label.trim();
        if trimmed.is_empty() {
            format!("{} #{}", feature.kind.to_string(), feature_id + 1)
        } else {
            trimmed.to_string()
        }
    }

    pub(super) fn normalize_column_frequencies(counts: [f64; 4]) -> [f64; 4] {
        let total = counts
            .iter()
            .copied()
            .map(|v| if v.is_finite() && v > 0.0 { v } else { 0.0 })
            .sum::<f64>();
        if total <= 0.0 {
            [0.25_f64, 0.25_f64, 0.25_f64, 0.25_f64]
        } else {
            [
                (counts[0].max(0.0) / total).clamp(0.0, 1.0),
                (counts[1].max(0.0) / total).clamp(0.0, 1.0),
                (counts[2].max(0.0) / total).clamp(0.0, 1.0),
                (counts[3].max(0.0) / total).clamp(0.0, 1.0),
            ]
        }
    }

    pub(super) fn information_content_bits(frequencies: &[f64; 4]) -> f64 {
        let entropy = frequencies
            .iter()
            .copied()
            .filter(|p| *p > 0.0)
            .map(|p| -p * p.log2())
            .sum::<f64>();
        (2.0 - entropy).clamp(0.0, 2.0)
    }

    pub(super) fn resolve_tfbs_scoring_motif(
        feature: &gb_io::seq::Feature,
    ) -> Result<(String, Option<String>, Vec<[f64; 4]>), EngineError> {
        let label = Self::feature_qualifier_text(feature, "label");
        let mut tokens = vec![];
        if let Some(tf_id) = Self::feature_qualifier_text(feature, "tf_id") {
            tokens.push(tf_id);
        }
        if let Some(bound) = Self::first_nonempty_feature_qualifier(
            feature,
            &["bound_moiety", "standard_name", "gene", "name"],
        ) {
            tokens.push(bound);
        }
        if let Some(raw_label) = label {
            let trimmed = raw_label.trim().to_string();
            if !trimmed.is_empty() {
                tokens.push(trimmed.clone());
                let label_upper = trimmed.to_ascii_uppercase();
                if let Some(stripped) = label_upper.strip_prefix("TFBS ") {
                    tokens.push(stripped.trim().to_string());
                }
            }
        }

        let mut last_error: Option<EngineError> = None;
        for token in tokens {
            if token.trim().is_empty() {
                continue;
            }
            match Self::resolve_tf_motif_for_scoring(&token) {
                Ok((tf_id, tf_name, _, matrix_counts)) => {
                    return Ok((tf_id, tf_name, matrix_counts));
                }
                Err(err) => {
                    last_error = Some(err);
                }
            }
        }
        Err(last_error.unwrap_or(EngineError {
            code: ErrorCode::InvalidInput,
            message:
                "Could not resolve TF motif for feature (missing tf_id/label resolvable in motif registry)"
                    .to_string(),
        }))
    }

    pub(super) fn build_tfbs_expert_view(
        &self,
        seq_id: &str,
        feature_id: usize,
    ) -> Result<TfbsExpertView, EngineError> {
        let dna = self
            .state
            .sequences
            .get(seq_id)
            .ok_or_else(|| EngineError {
                code: ErrorCode::NotFound,
                message: format!("Sequence '{seq_id}' not found"),
            })?;
        let feature = dna.features().get(feature_id).ok_or_else(|| EngineError {
            code: ErrorCode::NotFound,
            message: format!(
                "Feature id '{}' was not found in sequence '{}'",
                feature_id, seq_id
            ),
        })?;
        if !Self::is_tfbs_feature(feature) {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: format!(
                    "Feature '{}' in '{}' is not a TFBS/protein-bind feature",
                    feature_id, seq_id
                ),
            });
        }
        let (from, _to) = feature.location.find_bounds().map_err(|e| EngineError {
            code: ErrorCode::InvalidInput,
            message: format!("Could not read TFBS feature range: {e}"),
        })?;
        if from < 0 {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: "TFBS feature has negative range bounds".to_string(),
            });
        }

        let (tf_id, tf_name, matrix_counts) = Self::resolve_tfbs_scoring_motif(feature)?;
        if matrix_counts.is_empty() {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: format!("Resolved motif '{tf_id}' has empty matrix"),
            });
        }
        let motif_length = matrix_counts.len();
        let start = from as usize;
        let end = start.saturating_add(motif_length);
        let mut matched_bytes = dna.get_range_safe(start..end).ok_or_else(|| EngineError {
            code: ErrorCode::InvalidInput,
            message: format!(
                "TFBS window {}..{} for feature '{}' exceeds sequence bounds",
                start + 1,
                end,
                feature_id
            ),
        })?;
        let is_reverse = feature_is_reverse(feature);
        if is_reverse {
            matched_bytes = Self::reverse_complement_bytes(&matched_bytes);
        }
        let matched_sequence = String::from_utf8_lossy(&matched_bytes).to_ascii_uppercase();

        let (llr_matrix, true_log_odds_matrix) = Self::prepare_scoring_matrices(&matrix_counts);
        let mut columns: Vec<TfbsExpertColumn> = Vec::with_capacity(motif_length);
        for idx in 0..motif_length {
            let counts = matrix_counts[idx];
            let frequencies = Self::normalize_column_frequencies(counts);
            let information_content_bits = Self::information_content_bits(&frequencies);
            let match_base = matched_bytes
                .get(idx)
                .map(|b| (*b as char).to_ascii_uppercase());
            let match_idx = matched_bytes.get(idx).and_then(|b| Self::base_to_idx(*b));
            let (match_frequency, llr_bits, true_log_odds_bits) = if let Some(base_idx) = match_idx
            {
                (
                    Some(frequencies[base_idx]),
                    Some(llr_matrix[idx][base_idx]),
                    Some(true_log_odds_matrix[idx][base_idx]),
                )
            } else {
                (None, None, None)
            };
            columns.push(TfbsExpertColumn {
                index_1based: idx + 1,
                counts,
                frequencies,
                information_content_bits,
                match_base,
                match_frequency,
                llr_bits,
                true_log_odds_bits,
                llr_rank_desc: None,
            });
        }

        let mut ranked = columns
            .iter()
            .enumerate()
            .filter_map(|(idx, col)| col.llr_bits.map(|score| (idx, score)))
            .collect::<Vec<_>>();
        ranked.sort_by(|a, b| b.1.total_cmp(&a.1).then_with(|| a.0.cmp(&b.0)));
        for (rank_idx, (col_idx, _)) in ranked.iter().enumerate() {
            if let Some(col) = columns.get_mut(*col_idx) {
                col.llr_rank_desc = Some(rank_idx + 1);
            }
        }

        let llr_from_columns = columns.iter().filter_map(|col| col.llr_bits).sum::<f64>();
        let true_log_odds_from_columns = columns
            .iter()
            .filter_map(|col| col.true_log_odds_bits)
            .sum::<f64>();
        let llr_total_bits =
            Self::feature_qualifier_f64(feature, "llr_bits").or(Some(llr_from_columns));
        let llr_quantile = Self::feature_qualifier_f64(feature, "llr_quantile");
        let true_log_odds_total_bits = Self::feature_qualifier_f64(feature, "true_log_odds_bits")
            .or_else(|| Self::feature_qualifier_f64(feature, "log_odds_ratio_bits"))
            .or(Some(true_log_odds_from_columns));
        let true_log_odds_quantile = Self::feature_qualifier_f64(feature, "true_log_odds_quantile")
            .or_else(|| Self::feature_qualifier_f64(feature, "log_odds_ratio_quantile"));

        Ok(TfbsExpertView {
            seq_id: seq_id.to_string(),
            feature_id,
            feature_label: Self::feature_display_label(feature, feature_id),
            tf_id,
            tf_name,
            strand: if is_reverse {
                "-".to_string()
            } else {
                "+".to_string()
            },
            start_1based: start + 1,
            end_1based: end,
            motif_length,
            matched_sequence,
            llr_total_bits,
            llr_quantile,
            true_log_odds_total_bits,
            true_log_odds_quantile,
            instruction: TFBS_EXPERT_INSTRUCTION.to_string(),
            columns,
        })
    }

    pub(super) fn resolve_restriction_site_target(
        &self,
        seq_id: &str,
        cut_pos_1based: usize,
        enzyme: Option<&str>,
        recognition_start_1based: Option<usize>,
        recognition_end_1based: Option<usize>,
    ) -> Result<(RestrictionEnzymeKey, Vec<String>, Option<String>), EngineError> {
        if cut_pos_1based == 0 {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: "Restriction-site cut_pos_1based must be >= 1".to_string(),
            });
        }
        let dna = self
            .state
            .sequences
            .get(seq_id)
            .ok_or_else(|| EngineError {
                code: ErrorCode::NotFound,
                message: format!("Sequence '{seq_id}' not found"),
            })?;
        let enzyme_filter = enzyme
            .map(str::trim)
            .filter(|v| !v.is_empty())
            .map(|v| v.to_ascii_uppercase());
        let mut candidates: Vec<(RestrictionEnzymeKey, Vec<String>)> = Vec::new();
        for (key, names) in dna.restriction_enzyme_groups() {
            if key.pos() < 0 || (key.pos() as usize + 1) != cut_pos_1based {
                continue;
            }
            if let Some(start) = recognition_start_1based {
                if key.from() < 0 || (key.from() as usize + 1) != start {
                    continue;
                }
            }
            if let Some(end) = recognition_end_1based {
                if key.to() < 0 || key.to() as usize != end {
                    continue;
                }
            }
            if let Some(filter) = &enzyme_filter {
                if !names
                    .iter()
                    .any(|name| name.to_ascii_uppercase() == *filter)
                {
                    continue;
                }
            }
            candidates.push((key.clone(), names.clone()));
        }
        if candidates.is_empty() {
            return Err(EngineError {
                code: ErrorCode::NotFound,
                message: format!(
                    "No restriction-site match found for cut position {} in '{}'",
                    cut_pos_1based, seq_id
                ),
            });
        }
        if candidates.len() > 1 {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: format!(
                    "Restriction-site target is ambiguous for cut position {} in '{}'; provide enzyme and/or recognition start/end",
                    cut_pos_1based, seq_id
                ),
            });
        }
        let (key, mut names) = candidates.into_iter().next().unwrap_or_default();
        names.sort_by(|a, b| a.to_ascii_uppercase().cmp(&b.to_ascii_uppercase()));
        names.dedup_by(|left, right| left.eq_ignore_ascii_case(right));
        let selected_enzyme = if let Some(filter) = enzyme_filter {
            names
                .iter()
                .find(|name| name.to_ascii_uppercase() == filter)
                .cloned()
        } else {
            names.first().cloned()
        };
        Ok((key, names, selected_enzyme))
    }

    pub(super) fn complement_iupac_text(seq: &str) -> String {
        seq.as_bytes()
            .iter()
            .map(|b| {
                Self::iupac_letter_complement(*b)
                    .unwrap_or(b'N')
                    .to_ascii_uppercase() as char
            })
            .collect()
    }

    pub(super) fn build_restriction_site_expert_view(
        &self,
        seq_id: &str,
        cut_pos_1based: usize,
        enzyme: Option<&str>,
        recognition_start_1based: Option<usize>,
        recognition_end_1based: Option<usize>,
    ) -> Result<RestrictionSiteExpertView, EngineError> {
        let (key, enzyme_names, selected_enzyme) = self.resolve_restriction_site_target(
            seq_id,
            cut_pos_1based,
            enzyme,
            recognition_start_1based,
            recognition_end_1based,
        )?;
        let dna = self
            .state
            .sequences
            .get(seq_id)
            .ok_or_else(|| EngineError {
                code: ErrorCode::NotFound,
                message: format!("Sequence '{seq_id}' not found"),
            })?;
        let start_1based = key.from().max(0) as usize + 1;
        let end_1based = key.to().max(0) as usize;
        let start0 = start_1based.saturating_sub(1);
        let end0 = end_1based;
        let mut site_sequence = dna
            .get_range_safe(start0..end0)
            .map(|bytes| String::from_utf8_lossy(&bytes).to_ascii_uppercase())
            .unwrap_or_default();
        let selected_enzyme_def = selected_enzyme.as_ref().and_then(|name| {
            dna.restriction_enzymes()
                .iter()
                .find(|enzyme| enzyme.name.eq_ignore_ascii_case(name))
        });
        let recognition_iupac = selected_enzyme_def.map(|enzyme| enzyme.sequence.clone());
        let enzyme_cut_offset_0based = selected_enzyme_def.map(|enzyme| enzyme.cut);
        let overlap_bp = selected_enzyme_def.map(|enzyme| enzyme.overlap);
        let enzyme_note = selected_enzyme_def.and_then(|enzyme| {
            enzyme
                .note
                .as_deref()
                .map(str::trim)
                .filter(|note| !note.is_empty())
                .map(str::to_string)
        });
        let rebase_url = selected_enzyme_def
            .map(|enzyme| format!("https://rebase.neb.com/rebase/enz/{}.html", enzyme.name));
        if site_sequence.is_empty() {
            if let Some(iupac) = &recognition_iupac {
                site_sequence = iupac.to_ascii_uppercase();
            }
        }
        let site_sequence_complement = Self::complement_iupac_text(&site_sequence);
        let max_cut = site_sequence.chars().count();
        let cut_index_0based = key.pos().saturating_sub(key.from()).max(0) as usize;
        let cut_index_0based = cut_index_0based.min(max_cut);
        let paired_cut_pos_1based = key.mate_pos().max(0) as usize + 1;
        let paired_cut_index_0based = key.mate_pos().saturating_sub(key.from()).max(0) as usize;
        let paired_cut_index_0based = paired_cut_index_0based.min(max_cut);

        Ok(RestrictionSiteExpertView {
            seq_id: seq_id.to_string(),
            cut_pos_1based,
            paired_cut_pos_1based,
            recognition_start_1based: start_1based,
            recognition_end_1based: end_1based,
            cut_index_0based,
            paired_cut_index_0based,
            end_geometry: key.cut_geometry().kind_label().to_string(),
            number_of_cuts_for_enzyme: key.number_of_cuts(),
            selected_enzyme,
            enzyme_names,
            recognition_iupac,
            site_sequence,
            site_sequence_complement,
            enzyme_cut_offset_0based,
            overlap_bp,
            enzyme_note,
            rebase_url,
            instruction: RESTRICTION_EXPERT_INSTRUCTION.to_string(),
        })
    }

    pub(super) fn is_mrna_feature(feature: &gb_io::seq::Feature) -> bool {
        let kind = feature.kind.to_string();
        kind.eq_ignore_ascii_case("mRNA") || kind.eq_ignore_ascii_case("transcript")
    }

    pub(super) fn is_splicing_noncoding_transcript_feature(feature: &gb_io::seq::Feature) -> bool {
        let kind = feature.kind.to_string();
        kind.eq_ignore_ascii_case("ncRNA") || kind.eq_ignore_ascii_case("misc_RNA")
    }

    pub(super) fn is_splicing_transcript_feature(feature: &gb_io::seq::Feature) -> bool {
        Self::is_mrna_feature(feature) || Self::is_splicing_noncoding_transcript_feature(feature)
    }

    pub(super) fn is_splicing_seed_feature(feature: &gb_io::seq::Feature) -> bool {
        Self::is_splicing_transcript_feature(feature)
            || Self::is_exon_feature(feature)
            || feature.kind.to_string().eq_ignore_ascii_case("gene")
            || feature.kind.to_string().eq_ignore_ascii_case("cds")
    }

    pub(super) fn is_exon_feature(feature: &gb_io::seq::Feature) -> bool {
        feature.kind.to_string().eq_ignore_ascii_case("exon")
    }

    pub(super) fn splicing_group_label(feature: &gb_io::seq::Feature, fallback: usize) -> String {
        Self::first_nonempty_feature_qualifier(
            feature,
            &[
                "gene",
                "gene_id",
                "locus_tag",
                "standard_name",
                "gene_synonym",
            ],
        )
        .unwrap_or_else(|| format!("mRNA-group-{}", fallback + 1))
    }

    pub(super) fn feature_transcript_id(feature: &gb_io::seq::Feature, fallback: usize) -> String {
        Self::first_nonempty_feature_qualifier(
            feature,
            &[
                "transcript_id",
                "standard_name",
                "product",
                "label",
                "name",
                "gene",
            ],
        )
        .unwrap_or_else(|| format!("transcript-{}", fallback + 1))
    }

    pub(super) fn range_vec_to_splicing(mut ranges: Vec<(usize, usize)>) -> Vec<SplicingRange> {
        ranges.sort_unstable_by(|a, b| a.0.cmp(&b.0).then(a.1.cmp(&b.1)));
        ranges
            .into_iter()
            .filter(|(start, end)| end > start)
            .map(|(start, end)| SplicingRange {
                start_1based: start + 1,
                end_1based: end,
            })
            .collect()
    }

    pub(super) fn range_vec_1based_to_splicing(
        mut ranges: Vec<(usize, usize)>,
    ) -> Vec<SplicingRange> {
        ranges.sort_unstable_by(|a, b| a.0.cmp(&b.0).then(a.1.cmp(&b.1)));
        ranges
            .into_iter()
            .filter(|(start, end)| *start > 0 && *end >= *start)
            .map(|(start, end)| SplicingRange {
                start_1based: start,
                end_1based: end,
            })
            .collect()
    }

    pub(super) fn order_ranges_for_transcript(
        mut ranges: Vec<(usize, usize)>,
        is_reverse: bool,
    ) -> Vec<(usize, usize)> {
        ranges.sort_unstable_by(|a, b| a.0.cmp(&b.0).then(a.1.cmp(&b.1)));
        if is_reverse {
            ranges.reverse();
        }
        ranges
    }

    pub(super) fn cds_ranges_to_reference_aa_segments(
        cds_ranges: Vec<(usize, usize)>,
        is_reverse: bool,
        reference_start_aa: Option<usize>,
    ) -> Vec<IsoformArchitectureCdsAaSegment> {
        let ordered = Self::order_ranges_for_transcript(cds_ranges, is_reverse);
        if ordered.is_empty() {
            return vec![];
        }
        let reference_start = reference_start_aa.unwrap_or(1).max(1);
        let mut cumulative_nt = 0usize;
        let mut out = vec![];
        for (start_0based, end_0based_exclusive) in ordered {
            if end_0based_exclusive <= start_0based {
                continue;
            }
            let segment_nt = end_0based_exclusive - start_0based;
            let local_aa_start = cumulative_nt / 3 + 1;
            let local_aa_end = (cumulative_nt + segment_nt) / 3;
            cumulative_nt = cumulative_nt.saturating_add(segment_nt);
            if local_aa_end < local_aa_start {
                continue;
            }
            let aa_start = reference_start.saturating_add(local_aa_start.saturating_sub(1));
            let aa_end = reference_start.saturating_add(local_aa_end.saturating_sub(1));
            out.push(IsoformArchitectureCdsAaSegment {
                genomic_start_1based: start_0based.saturating_add(1),
                genomic_end_1based: end_0based_exclusive,
                aa_start,
                aa_end,
            });
        }
        out
    }

    pub(super) fn sequence_slice_upper(dna: &DNAsequence, start: usize, end: usize) -> String {
        if end <= start {
            return String::new();
        }
        dna.get_range_safe(start..end)
            .map(|bytes| String::from_utf8_lossy(&bytes).to_ascii_uppercase())
            .unwrap_or_default()
    }

    pub(super) fn normalize_transcript_probe(raw: &str) -> String {
        let trimmed = raw.trim();
        if trimmed.is_empty() {
            return String::new();
        }
        let mut out = trimmed.to_ascii_uppercase();
        if let Some((prefix, suffix)) = out.rsplit_once('.')
            && !prefix.is_empty()
            && !suffix.is_empty()
            && suffix.chars().all(|ch| ch.is_ascii_digit())
        {
            out = prefix.to_string();
        }
        out
    }

    pub(super) fn feature_transcript_match_keys(
        feature: &gb_io::seq::Feature,
        fallback: usize,
    ) -> Vec<String> {
        let mut keys: BTreeSet<String> = BTreeSet::new();
        for key in [
            "transcript_id",
            "standard_name",
            "label",
            "name",
            "product",
            "protein_id",
        ] {
            if let Some(value) = Self::feature_qualifier_text(feature, key) {
                let normalized = Self::normalize_transcript_probe(&value);
                if !normalized.is_empty() {
                    keys.insert(normalized);
                }
            }
        }
        let fallback_id =
            Self::normalize_transcript_probe(&Self::feature_transcript_id(feature, fallback));
        if !fallback_id.is_empty() {
            keys.insert(fallback_id);
        }
        keys.into_iter().collect()
    }

    pub(super) fn normalize_uniprot_projection_id(raw: &str) -> Result<String, EngineError> {
        let trimmed = raw.trim();
        if trimmed.is_empty() {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: "projection_id cannot be empty".to_string(),
            });
        }
        let mut out = String::with_capacity(trimmed.len());
        for ch in trimmed.chars() {
            if ch.is_ascii_alphanumeric() || matches!(ch, '_' | '-' | '.' | '@') {
                out.push(ch);
            } else {
                out.push('_');
            }
        }
        let normalized = out.trim_matches('_').to_string();
        if normalized.is_empty() {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message:
                    "projection_id must contain at least one ASCII letter, digit, '.', '-', '_' or '@'"
                        .to_string(),
            });
        }
        Ok(normalized)
    }

    pub(super) fn load_ensembl_protein_entry_store_from_metadata(
        value: Option<&serde_json::Value>,
    ) -> EnsemblProteinEntryStore {
        let mut store = value
            .cloned()
            .and_then(|v| serde_json::from_value::<EnsemblProteinEntryStore>(v).ok())
            .unwrap_or_default();
        if store.schema.trim().is_empty() {
            store.schema = ENSEMBL_PROTEIN_ENTRIES_SCHEMA.to_string();
        }
        store
    }

    pub(super) fn read_ensembl_protein_entry_store(&self) -> EnsemblProteinEntryStore {
        Self::load_ensembl_protein_entry_store_from_metadata(
            self.state
                .metadata
                .get(ENSEMBL_PROTEIN_ENTRIES_METADATA_KEY),
        )
    }

    pub(super) fn write_ensembl_protein_entry_store(
        &mut self,
        mut store: EnsemblProteinEntryStore,
    ) -> Result<(), EngineError> {
        if store.entries.is_empty() {
            self.state
                .metadata
                .remove(ENSEMBL_PROTEIN_ENTRIES_METADATA_KEY);
            return Ok(());
        }
        store.schema = ENSEMBL_PROTEIN_ENTRIES_SCHEMA.to_string();
        store.updated_at_unix_ms = Self::now_unix_ms();
        let value = serde_json::to_value(store).map_err(|e| EngineError {
            code: ErrorCode::Internal,
            message: format!("Could not serialize Ensembl protein entry metadata: {e}"),
        })?;
        self.state
            .metadata
            .insert(ENSEMBL_PROTEIN_ENTRIES_METADATA_KEY.to_string(), value);
        Ok(())
    }

    pub(crate) fn upsert_ensembl_protein_entry(
        &mut self,
        mut entry: EnsemblProteinEntry,
    ) -> Result<(), EngineError> {
        let normalized = normalize_ensembl_protein_entry_id(&entry.entry_id);
        if normalized.is_empty() {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: "Ensembl protein entry_id cannot be empty".to_string(),
            });
        }
        entry.entry_id = normalized.clone();
        if entry.schema.trim().is_empty() {
            entry.schema = "gentle.ensembl_protein_entry.v1".to_string();
        }
        entry.imported_at_unix_ms = Self::now_unix_ms();
        let mut aliases = BTreeSet::new();
        for value in std::iter::once(entry.entry_id.clone())
            .chain(std::iter::once(entry.protein_id.clone()))
            .chain(std::iter::once(entry.transcript_id.clone()))
            .chain(entry.aliases.clone().into_iter())
            .chain(entry.gene_id.clone().into_iter())
            .chain(entry.gene_symbol.clone().into_iter())
            .chain(entry.transcript_display_name.clone().into_iter())
        {
            let normalized_alias = normalize_ensembl_protein_entry_id(&value);
            if !normalized_alias.is_empty() {
                aliases.insert(normalized_alias);
            }
        }
        entry.aliases = aliases.into_iter().collect::<Vec<_>>();

        let mut store = self.read_ensembl_protein_entry_store();
        store.entries.insert(entry.entry_id.clone(), entry);
        self.write_ensembl_protein_entry_store(store)
    }

    pub(super) fn get_ensembl_protein_entry_from_store(
        store: &EnsemblProteinEntryStore,
        entry_id: &str,
    ) -> Option<EnsemblProteinEntry> {
        let probe = normalize_ensembl_protein_entry_id(entry_id);
        if probe.is_empty() {
            return None;
        }
        if let Some(entry) = store.entries.get(&probe) {
            return Some(entry.clone());
        }
        store
            .entries
            .values()
            .find(|entry| {
                entry
                    .aliases
                    .iter()
                    .any(|alias| normalize_ensembl_protein_entry_id(alias) == probe)
            })
            .cloned()
    }

    pub fn get_ensembl_protein_entry(
        &self,
        entry_id: &str,
    ) -> Result<EnsemblProteinEntry, EngineError> {
        let store = self.read_ensembl_protein_entry_store();
        Self::get_ensembl_protein_entry_from_store(&store, entry_id).ok_or_else(|| EngineError {
            code: ErrorCode::NotFound,
            message: format!(
                "Ensembl protein entry '{}' not found in project metadata",
                entry_id.trim()
            ),
        })
    }

    pub fn list_ensembl_protein_entries(&self) -> Vec<EnsemblProteinEntrySummary> {
        let store = self.read_ensembl_protein_entry_store();
        let mut rows = store
            .entries
            .values()
            .map(|entry| EnsemblProteinEntrySummary {
                entry_id: entry.entry_id.clone(),
                protein_id: entry.protein_id.clone(),
                transcript_id: entry.transcript_id.clone(),
                gene_symbol: entry.gene_symbol.clone(),
                sequence_length: entry.sequence_length,
                feature_count: entry.features.len(),
                imported_at_unix_ms: entry.imported_at_unix_ms,
                source_query: entry.source_query.clone(),
            })
            .collect::<Vec<_>>();
        rows.sort_by(|a, b| a.entry_id.cmp(&b.entry_id));
        rows
    }

    pub(super) fn load_uniprot_entry_store_from_metadata(
        value: Option<&serde_json::Value>,
    ) -> UniprotEntryStore {
        let mut store = value
            .cloned()
            .and_then(|v| serde_json::from_value::<UniprotEntryStore>(v).ok())
            .unwrap_or_default();
        if store.schema.trim().is_empty() {
            store.schema = UNIPROT_ENTRIES_SCHEMA.to_string();
        }
        store
    }

    pub(super) fn load_uniprot_projection_store_from_metadata(
        value: Option<&serde_json::Value>,
    ) -> UniprotGenomeProjectionStore {
        let mut store = value
            .cloned()
            .and_then(|v| serde_json::from_value::<UniprotGenomeProjectionStore>(v).ok())
            .unwrap_or_default();
        if store.schema.trim().is_empty() {
            store.schema = UNIPROT_GENOME_PROJECTIONS_SCHEMA.to_string();
        }
        store
    }

    pub(super) fn read_uniprot_entry_store(&self) -> UniprotEntryStore {
        Self::load_uniprot_entry_store_from_metadata(
            self.state.metadata.get(UNIPROT_ENTRIES_METADATA_KEY),
        )
    }

    pub(super) fn write_uniprot_entry_store(
        &mut self,
        mut store: UniprotEntryStore,
    ) -> Result<(), EngineError> {
        if store.entries.is_empty() {
            self.state.metadata.remove(UNIPROT_ENTRIES_METADATA_KEY);
            return Ok(());
        }
        store.schema = UNIPROT_ENTRIES_SCHEMA.to_string();
        store.updated_at_unix_ms = Self::now_unix_ms();
        let value = serde_json::to_value(store).map_err(|e| EngineError {
            code: ErrorCode::Internal,
            message: format!("Could not serialize UniProt entry metadata: {e}"),
        })?;
        self.state
            .metadata
            .insert(UNIPROT_ENTRIES_METADATA_KEY.to_string(), value);
        Ok(())
    }

    pub(super) fn read_uniprot_projection_store(&self) -> UniprotGenomeProjectionStore {
        Self::load_uniprot_projection_store_from_metadata(
            self.state
                .metadata
                .get(UNIPROT_GENOME_PROJECTIONS_METADATA_KEY),
        )
    }

    pub(super) fn write_uniprot_projection_store(
        &mut self,
        mut store: UniprotGenomeProjectionStore,
    ) -> Result<(), EngineError> {
        if store.projections.is_empty() {
            self.state
                .metadata
                .remove(UNIPROT_GENOME_PROJECTIONS_METADATA_KEY);
            return Ok(());
        }
        store.schema = UNIPROT_GENOME_PROJECTIONS_SCHEMA.to_string();
        store.updated_at_unix_ms = Self::now_unix_ms();
        let value = serde_json::to_value(store).map_err(|e| EngineError {
            code: ErrorCode::Internal,
            message: format!("Could not serialize UniProt projection metadata: {e}"),
        })?;
        self.state
            .metadata
            .insert(UNIPROT_GENOME_PROJECTIONS_METADATA_KEY.to_string(), value);
        Ok(())
    }

    pub(super) fn upsert_uniprot_entry(
        &mut self,
        mut entry: UniprotEntry,
    ) -> Result<(), EngineError> {
        let normalized = normalize_uniprot_entry_id(&entry.entry_id);
        if normalized.is_empty() {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: "UniProt entry_id cannot be empty".to_string(),
            });
        }
        entry.entry_id = normalized.clone();
        if entry.schema.trim().is_empty() {
            entry.schema = "gentle.uniprot_entry.v1".to_string();
        }
        entry.imported_at_unix_ms = Self::now_unix_ms();
        let mut aliases = BTreeSet::new();
        for value in std::iter::once(entry.entry_id.clone())
            .chain(entry.aliases.clone().into_iter())
            .chain(std::iter::once(entry.accession.clone()))
            .chain(std::iter::once(entry.primary_id.clone()))
        {
            let normalized_alias = normalize_uniprot_entry_id(&value);
            if !normalized_alias.is_empty() {
                aliases.insert(normalized_alias);
            }
        }
        entry.aliases = aliases.into_iter().collect::<Vec<_>>();

        let mut store = self.read_uniprot_entry_store();
        store.entries.insert(entry.entry_id.clone(), entry);
        self.write_uniprot_entry_store(store)
    }

    pub(super) fn get_uniprot_entry_from_store(
        store: &UniprotEntryStore,
        entry_id: &str,
    ) -> Option<UniprotEntry> {
        let probe = normalize_uniprot_entry_id(entry_id);
        if probe.is_empty() {
            return None;
        }
        if let Some(entry) = store.entries.get(&probe) {
            return Some(entry.clone());
        }
        store
            .entries
            .values()
            .find(|entry| {
                entry
                    .aliases
                    .iter()
                    .any(|alias| normalize_uniprot_entry_id(alias) == probe)
            })
            .cloned()
    }

    pub fn get_uniprot_entry(&self, entry_id: &str) -> Result<UniprotEntry, EngineError> {
        let store = self.read_uniprot_entry_store();
        Self::get_uniprot_entry_from_store(&store, entry_id).ok_or_else(|| EngineError {
            code: ErrorCode::NotFound,
            message: format!(
                "UniProt entry '{}' not found in project metadata",
                entry_id.trim()
            ),
        })
    }

    pub fn list_uniprot_entries(&self) -> Vec<UniprotEntrySummary> {
        let store = self.read_uniprot_entry_store();
        let mut rows = store
            .entries
            .values()
            .map(|entry| UniprotEntrySummary {
                entry_id: entry.entry_id.clone(),
                accession: entry.accession.clone(),
                primary_id: entry.primary_id.clone(),
                sequence_length: entry.sequence_length,
                feature_count: entry.features.len(),
                ensembl_xref_count: entry.ensembl_xrefs.len(),
                nucleotide_xref_count: entry.nucleotide_xrefs.len(),
                imported_at_unix_ms: entry.imported_at_unix_ms,
                source: entry.source.clone(),
            })
            .collect::<Vec<_>>();
        rows.sort_by(|a, b| a.entry_id.cmp(&b.entry_id));
        rows
    }

    pub fn get_uniprot_genome_projection(
        &self,
        projection_id: &str,
    ) -> Result<UniprotGenomeProjection, EngineError> {
        let projection_id = Self::normalize_uniprot_projection_id(projection_id)?;
        let store = self.read_uniprot_projection_store();
        store
            .projections
            .get(&projection_id)
            .cloned()
            .ok_or_else(|| EngineError {
                code: ErrorCode::NotFound,
                message: format!(
                    "UniProt genome projection '{}' was not found",
                    projection_id
                ),
            })
    }

    pub fn list_uniprot_genome_projections(
        &self,
        seq_filter: Option<&str>,
    ) -> Vec<UniprotGenomeProjectionSummary> {
        let seq_filter = seq_filter
            .map(str::trim)
            .filter(|v| !v.is_empty())
            .map(|v| v.to_string());
        let store = self.read_uniprot_projection_store();
        let mut rows = store
            .projections
            .values()
            .filter(|record| {
                seq_filter
                    .as_deref()
                    .map(|seq_id| record.seq_id == seq_id)
                    .unwrap_or(true)
            })
            .map(|record| UniprotGenomeProjectionSummary {
                projection_id: record.projection_id.clone(),
                entry_id: record.entry_id.clone(),
                seq_id: record.seq_id.clone(),
                created_at_unix_ms: record.created_at_unix_ms,
                op_id: record.op_id.clone(),
                run_id: record.run_id.clone(),
                transcript_id_filter: record.transcript_id_filter.clone(),
                transcript_projection_count: record.transcript_projections.len(),
            })
            .collect::<Vec<_>>();
        rows.sort_by(|a, b| a.projection_id.cmp(&b.projection_id));
        rows
    }

    fn normalize_uniprot_variant_feature_key(feature: &UniprotFeature) -> Option<String> {
        let key = feature.key.trim().to_ascii_uppercase();
        matches!(key.as_str(), "VARIANT" | "VAR_SEQ" | "CONFLICT").then_some(key)
    }

    fn uniprot_init_met_declared(entry: &UniprotEntry) -> bool {
        entry.features.iter().any(|feature| {
            feature.key.eq_ignore_ascii_case("INIT_MET")
                && feature.interval.start_aa == Some(1)
                && feature.interval.end_aa.unwrap_or(1) >= 1
        })
    }

    fn derive_uniprot_projection_transcript_context(
        &self,
        seq_id: &str,
        transcript_projection: &UniprotTranscriptProjection,
    ) -> Result<DerivedProteinExpertTranscript, EngineError> {
        let dna = self
            .state
            .sequences
            .get(seq_id)
            .ok_or_else(|| EngineError {
                code: ErrorCode::NotFound,
                message: format!("Sequence '{seq_id}' not found"),
            })?;
        let features = dna.features();
        let source_sequence_upper = dna.get_forward_string().to_ascii_uppercase().into_bytes();
        if let Some(feature_id) = transcript_projection.transcript_feature_id
            && let Some(feature) = features.get(feature_id)
            && Self::is_transcript_feature_for_derivation(feature)
        {
            return Self::build_derived_protein_expert_transcript(
                &source_sequence_upper,
                feature,
                features,
                feature_id,
                seq_id,
            );
        }
        let normalized_probe =
            Self::normalize_transcript_probe(&transcript_projection.transcript_id);
        for (feature_id, feature) in features.iter().enumerate() {
            if !Self::is_transcript_feature_for_derivation(feature) {
                continue;
            }
            let keys = Self::feature_transcript_match_keys(feature, feature_id);
            if !keys.iter().any(|key| key == &normalized_probe) {
                continue;
            }
            return Self::build_derived_protein_expert_transcript(
                &source_sequence_upper,
                feature,
                features,
                feature_id,
                seq_id,
            );
        }
        Err(EngineError {
            code: ErrorCode::NotFound,
            message: format!(
                "Projected transcript '{}' was not found as a derivable transcript in sequence '{}'",
                transcript_projection.transcript_id, seq_id
            ),
        })
    }

    fn overlap_1based(left: (usize, usize), right: (usize, usize)) -> Option<(usize, usize)> {
        let start = left.0.max(right.0);
        let end = left.1.min(right.1);
        (start > 0 && end >= start).then_some((start, end))
    }

    fn build_transcript_exon_contributions(
        transcript_exons_1based: &[(usize, usize)],
        cds_ranges_1based: &[(usize, usize)],
        is_reverse: bool,
    ) -> Vec<UniprotTranscriptExonContribution> {
        let ordered_exons =
            Self::order_ranges_for_transcript(transcript_exons_1based.to_vec(), is_reverse);
        ordered_exons
            .into_iter()
            .enumerate()
            .map(|(idx, exon)| {
                let overlaps = cds_ranges_1based
                    .iter()
                    .filter_map(|coding| Self::overlap_1based(exon, *coding))
                    .collect::<Vec<_>>();
                let coding_start_1based = overlaps.iter().map(|(start, _)| *start).min();
                let coding_end_1based = overlaps.iter().map(|(_, end)| *end).max();
                let coding_nt = overlaps
                    .iter()
                    .map(|(start, end)| end.saturating_sub(*start).saturating_add(1))
                    .sum::<usize>();
                UniprotTranscriptExonContribution {
                    ordinal: idx + 1,
                    exon_start_1based: exon.0,
                    exon_end_1based: exon.1,
                    exon_nt: exon.1.saturating_sub(exon.0).saturating_add(1),
                    coding_start_1based,
                    coding_end_1based,
                    coding_nt,
                }
            })
            .collect()
    }

    fn summarize_translation_accounting(
        transcript_id: &str,
        transcript_label: &str,
        transcript_feature_id: Option<usize>,
        is_reverse: bool,
        derivation: Option<&TranscriptProteinDerivation>,
        transcript_exons_1based: &[(usize, usize)],
        entry: &UniprotEntry,
    ) -> UniprotProjectionTranscriptAccountingRow {
        let strand = if is_reverse { "-" } else { "+" }.to_string();
        let cds_ranges_1based = derivation
            .map(|value| value.cds_ranges_1based.clone())
            .unwrap_or_default();
        let contributions = Self::build_transcript_exon_contributions(
            transcript_exons_1based,
            &cds_ranges_1based,
            is_reverse,
        );
        let contributing = contributions
            .iter()
            .filter(|exon| exon.coding_nt > 0)
            .cloned()
            .collect::<Vec<_>>();
        let total_exon_nt = contributions.iter().map(|exon| exon.exon_nt).sum::<usize>();
        let translated_nt = contributing
            .iter()
            .map(|exon| exon.coding_nt)
            .sum::<usize>();
        let untranslated_5prime_nt = if let Some(first) = contributing.first() {
            let nt_before = contributions
                .iter()
                .take_while(|exon| exon.ordinal < first.ordinal)
                .map(|exon| exon.exon_nt)
                .sum::<usize>();
            let within = match (first.coding_start_1based, first.coding_end_1based) {
                (Some(_coding_start), Some(coding_end)) if is_reverse => {
                    first.exon_end_1based.saturating_sub(coding_end)
                }
                (Some(coding_start), Some(_)) => {
                    coding_start.saturating_sub(first.exon_start_1based)
                }
                _ => 0,
            };
            nt_before.saturating_add(within)
        } else {
            total_exon_nt
        };
        let untranslated_3prime_nt = total_exon_nt
            .saturating_sub(untranslated_5prime_nt)
            .saturating_sub(translated_nt);
        let expected_aa_count = translated_nt / 3;
        let init_met_declared = Self::uniprot_init_met_declared(entry);
        UniprotProjectionTranscriptAccountingRow {
            transcript_id: transcript_id.to_string(),
            transcript_label: transcript_label.to_string(),
            transcript_feature_id,
            strand,
            contributing_exons: contributing,
            contributing_exon_nt_sum: contributions
                .iter()
                .filter(|exon| exon.coding_nt > 0)
                .map(|exon| exon.exon_nt)
                .sum(),
            untranslated_5prime_nt,
            untranslated_3prime_nt,
            translated_nt,
            translated_nt_divisible_by_3: translated_nt % 3 == 0,
            expected_aa_count,
            derived_protein_length_aa: derivation.map(|value| value.protein_length_aa).unwrap_or(0),
            derived_protein_sequence: derivation
                .map(|value| value.protein_sequence.clone())
                .unwrap_or_default(),
            uniprot_aa_count: entry.sequence_length,
            init_met_declared,
            warnings: derivation
                .map(|value| value.warnings.clone())
                .unwrap_or_default(),
        }
    }

    fn build_uniprot_link_resolution_row(
        entry: &UniprotEntry,
        transcript_projection: &UniprotTranscriptProjection,
    ) -> UniprotEnsemblLinkResolutionRow {
        let normalized_transcript =
            Self::normalize_transcript_probe(&transcript_projection.transcript_id);
        let matched = entry
            .ensembl_xrefs
            .iter()
            .filter(|xref| {
                xref.transcript_id
                    .as_deref()
                    .map(Self::normalize_transcript_probe)
                    .is_some_and(|candidate| candidate == normalized_transcript)
            })
            .cloned()
            .collect::<Vec<_>>();
        let normalized_xref_transcript_ids = matched
            .iter()
            .filter_map(|xref| xref.transcript_id.as_deref())
            .map(Self::normalize_transcript_probe)
            .filter(|value| !value.is_empty())
            .collect::<BTreeSet<_>>()
            .into_iter()
            .collect::<Vec<_>>();
        let normalized_xref_protein_ids = matched
            .iter()
            .filter_map(|xref| xref.protein_id.as_deref())
            .map(normalize_ensembl_protein_entry_id)
            .filter(|value| !value.is_empty())
            .collect::<BTreeSet<_>>()
            .into_iter()
            .collect::<Vec<_>>();
        let normalized_xref_gene_ids = matched
            .iter()
            .filter_map(|xref| xref.gene_id.as_deref())
            .map(normalize_ensembl_protein_entry_id)
            .filter(|value| !value.is_empty())
            .collect::<BTreeSet<_>>()
            .into_iter()
            .collect::<Vec<_>>();
        let matched_xrefs = matched
            .iter()
            .map(|xref| UniprotEnsemblLinkedXref {
                transcript_id: xref.transcript_id.clone(),
                protein_id: xref.protein_id.clone(),
                gene_id: xref.gene_id.clone(),
                isoform_id: xref.isoform_id.clone(),
            })
            .collect::<Vec<_>>();
        let mut diagnostics = vec![];
        let status = if matched_xrefs.is_empty() {
            let available = entry
                .ensembl_xrefs
                .iter()
                .filter_map(|xref| xref.transcript_id.clone())
                .collect::<Vec<_>>();
            if !available.is_empty() {
                diagnostics.push(format!(
                    "No UniProt Ensembl xref matched transcript '{}'; available transcript xrefs: {}",
                    transcript_projection.transcript_id,
                    available.join(", ")
                ));
            }
            "missing".to_string()
        } else if normalized_xref_protein_ids.len() > 1 || normalized_xref_gene_ids.len() > 1 {
            diagnostics.push(format!(
                "Matched UniProt Ensembl links for '{}' are ambiguous across protein/gene IDs.",
                transcript_projection.transcript_id
            ));
            "ambiguous".to_string()
        } else {
            "matched".to_string()
        };
        UniprotEnsemblLinkResolutionRow {
            transcript_id: transcript_projection.transcript_id.clone(),
            transcript_feature_id: transcript_projection.transcript_feature_id,
            matched_xrefs,
            normalized_xref_transcript_ids,
            normalized_xref_protein_ids,
            normalized_xref_gene_ids,
            status,
            diagnostics,
        }
    }

    fn resolve_ensembl_entry_for_projection_transcript(
        &self,
        explicit_entry_id: Option<&str>,
        link_row: &UniprotEnsemblLinkResolutionRow,
    ) -> Result<Option<EnsemblProteinEntry>, EngineError> {
        let store = self.read_ensembl_protein_entry_store();
        if let Some(entry_id) = explicit_entry_id
            .map(str::trim)
            .filter(|value| !value.is_empty())
        {
            return Self::get_ensembl_protein_entry_from_store(&store, entry_id)
                .map(Some)
                .ok_or_else(|| EngineError {
                    code: ErrorCode::NotFound,
                    message: format!(
                        "Ensembl protein entry '{}' not found in project metadata",
                        entry_id
                    ),
                });
        }
        let mut candidates = store
            .entries
            .values()
            .filter(|entry| {
                let transcript_match =
                    link_row
                        .normalized_xref_transcript_ids
                        .iter()
                        .any(|candidate| {
                            Self::normalize_transcript_probe(&entry.transcript_id) == *candidate
                        });
                let protein_match = link_row
                    .normalized_xref_protein_ids
                    .iter()
                    .any(|candidate| {
                        normalize_ensembl_protein_entry_id(&entry.protein_id) == *candidate
                    });
                transcript_match || protein_match
            })
            .cloned()
            .collect::<Vec<_>>();
        candidates.sort_by(|left, right| left.entry_id.cmp(&right.entry_id));
        if candidates.len() <= 1 {
            Ok(candidates.into_iter().next())
        } else {
            Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: format!(
                    "Multiple imported Ensembl protein entries match transcript '{}': {}",
                    link_row.transcript_id,
                    candidates
                        .iter()
                        .map(|entry| entry.entry_id.clone())
                        .collect::<Vec<_>>()
                        .join(", ")
                ),
            })
        }
    }

    fn build_ensembl_accounting_row(
        transcript_id: &str,
        entry: &EnsemblProteinEntry,
        uniprot_entry: &UniprotEntry,
    ) -> Option<UniprotProjectionTranscriptAccountingRow> {
        let translation = entry.transcript_translation.as_ref()?;
        let transcript_exons = entry
            .transcript_exons
            .iter()
            .map(|exon| (exon.start_1based, exon.end_1based))
            .collect::<Vec<_>>();
        let cds_range = (
            translation.genomic_start_1based?,
            translation.genomic_end_1based?,
        );
        let is_reverse = entry
            .transcript_exons
            .first()
            .map(|exon| exon.strand < 0)
            .unwrap_or(false);
        Some(Self::summarize_translation_accounting(
            transcript_id,
            entry
                .transcript_display_name
                .as_deref()
                .unwrap_or(transcript_id),
            None,
            is_reverse,
            Some(&TranscriptProteinDerivation {
                transcript_id: transcript_id.to_string(),
                transcript_label: entry
                    .transcript_display_name
                    .clone()
                    .unwrap_or_else(|| transcript_id.to_string()),
                source_seq_id: entry.entry_id.clone(),
                source_feature_id: 0,
                derivation_mode: TranscriptProteinDerivationMode::AnnotatedCds,
                cds_ranges_1based: vec![cds_range],
                cds_length_bp: if cds_range.1 >= cds_range.0 {
                    cds_range.1 - cds_range.0 + 1
                } else {
                    0
                },
                protein_sequence: entry.sequence.clone(),
                protein_length_aa: entry.sequence_length,
                translation_table: 1,
                translation_table_label: "Standard".to_string(),
                translation_table_source: TranscriptProteinTranslationTableSource::StandardDefault,
                codon_start: 1,
                organism: entry.species.clone(),
                organelle: None,
                translation_speed_profile_hint: None,
                translation_speed_profile_source: None,
                translation_speed_reference_species: None,
                terminal_stop_trimmed: false,
                warnings: vec![],
            }),
            &transcript_exons,
            uniprot_entry,
        ))
    }

    fn build_exon_compare_row(
        current: &UniprotProjectionTranscriptAccountingRow,
        ensembl_entry: Option<&EnsemblProteinEntry>,
        ensembl_accounting: Option<&UniprotProjectionTranscriptAccountingRow>,
    ) -> UniprotEnsemblExonCompareRow {
        if let (Some(entry), Some(ensembl)) = (ensembl_entry, ensembl_accounting) {
            let current_by_ordinal = current
                .contributing_exons
                .iter()
                .map(|row| (row.ordinal, row))
                .collect::<BTreeMap<_, _>>();
            let ensembl_by_ordinal = ensembl
                .contributing_exons
                .iter()
                .map(|row| (row.ordinal, row))
                .collect::<BTreeMap<_, _>>();
            let matched_exon_ordinals = current_by_ordinal
                .iter()
                .filter_map(|(ordinal, current_row)| {
                    ensembl_by_ordinal.get(ordinal).and_then(|ensembl_row| {
                        (current_row.coding_nt == ensembl_row.coding_nt).then_some(*ordinal)
                    })
                })
                .collect::<Vec<_>>();
            let missing_in_ensembl_ordinals = current_by_ordinal
                .keys()
                .filter(|ordinal| !ensembl_by_ordinal.contains_key(ordinal))
                .copied()
                .collect::<Vec<_>>();
            let excess_in_ensembl_ordinals = ensembl_by_ordinal
                .keys()
                .filter(|ordinal| !current_by_ordinal.contains_key(ordinal))
                .copied()
                .collect::<Vec<_>>();
            let mut boundary_differences = vec![];
            for ordinal in current_by_ordinal
                .keys()
                .filter(|ordinal| ensembl_by_ordinal.contains_key(ordinal))
            {
                let current_row = current_by_ordinal[ordinal];
                let ensembl_row = ensembl_by_ordinal[ordinal];
                if current_row.coding_start_1based != ensembl_row.coding_start_1based {
                    boundary_differences.push(UniprotEnsemblExonBoundaryDifference {
                        ordinal: *ordinal,
                        side: "start".to_string(),
                        current_coordinate_1based: current_row.coding_start_1based.unwrap_or(0),
                        ensembl_coordinate_1based: ensembl_row.coding_start_1based.unwrap_or(0),
                        note: "coding start differs".to_string(),
                    });
                }
                if current_row.coding_end_1based != ensembl_row.coding_end_1based {
                    boundary_differences.push(UniprotEnsemblExonBoundaryDifference {
                        ordinal: *ordinal,
                        side: "end".to_string(),
                        current_coordinate_1based: current_row.coding_end_1based.unwrap_or(0),
                        ensembl_coordinate_1based: ensembl_row.coding_end_1based.unwrap_or(0),
                        note: "coding end differs".to_string(),
                    });
                }
            }
            let status = if missing_in_ensembl_ordinals.is_empty()
                && excess_in_ensembl_ordinals.is_empty()
                && boundary_differences.is_empty()
            {
                "matched".to_string()
            } else {
                "mismatch".to_string()
            };
            UniprotEnsemblExonCompareRow {
                transcript_id: current.transcript_id.clone(),
                transcript_label: current.transcript_label.clone(),
                transcript_feature_id: current.transcript_feature_id,
                ensembl_entry_id: Some(entry.entry_id.clone()),
                ensembl_transcript_id: Some(entry.transcript_id.clone()),
                current_contributing_exons: current.contributing_exons.clone(),
                ensembl_contributing_exons: ensembl.contributing_exons.clone(),
                matched_exon_ordinals,
                missing_in_ensembl_ordinals,
                excess_in_ensembl_ordinals,
                boundary_differences,
                status,
                warnings: vec![],
            }
        } else {
            UniprotEnsemblExonCompareRow {
                transcript_id: current.transcript_id.clone(),
                transcript_label: current.transcript_label.clone(),
                transcript_feature_id: current.transcript_feature_id,
                ensembl_entry_id: ensembl_entry.map(|entry| entry.entry_id.clone()),
                ensembl_transcript_id: ensembl_entry.map(|entry| entry.transcript_id.clone()),
                current_contributing_exons: current.contributing_exons.clone(),
                ensembl_contributing_exons: vec![],
                matched_exon_ordinals: vec![],
                missing_in_ensembl_ordinals: current
                    .contributing_exons
                    .iter()
                    .map(|exon| exon.ordinal)
                    .collect(),
                excess_in_ensembl_ordinals: vec![],
                boundary_differences: vec![],
                status: "missing_ensembl_evidence".to_string(),
                warnings: vec![
                    "No imported Ensembl exon/CDS evidence was available for this transcript."
                        .to_string(),
                ],
            }
        }
    }

    fn build_peptide_compare_row(
        &self,
        current: &UniprotProjectionTranscriptAccountingRow,
        ensembl_entry: Option<&EnsemblProteinEntry>,
        entry: &UniprotEntry,
    ) -> Result<UniprotEnsemblPeptideCompareRow, EngineError> {
        let derived_sequence = current.derived_protein_sequence.clone();
        let uniprot_sequence = entry.sequence.clone();
        let ensembl_sequence = ensembl_entry.map(|value| value.sequence.clone());
        let init_met_declared = current.init_met_declared;
        let comparison_mode = if derived_sequence.len() == uniprot_sequence.len() {
            UniprotPeptideComparisonMode::DirectCompare
        } else {
            UniprotPeptideComparisonMode::GlobalAlignment
        };
        let mut direct_mismatches = vec![];
        let mut global_alignment = None;
        let mut status = "matched".to_string();
        if matches!(comparison_mode, UniprotPeptideComparisonMode::DirectCompare) {
            let ensembl_bytes = ensembl_sequence.as_ref().map(|value| value.as_bytes());
            for (idx, (derived, uniprot)) in derived_sequence
                .as_bytes()
                .iter()
                .zip(uniprot_sequence.as_bytes().iter())
                .enumerate()
            {
                if derived != uniprot {
                    direct_mismatches.push(UniprotPeptideDirectMismatch {
                        position_1based: idx + 1,
                        gentle_residue: (*derived as char).to_string(),
                        uniprot_residue: (*uniprot as char).to_string(),
                        ensembl_residue: ensembl_bytes
                            .and_then(|bytes| bytes.get(idx).copied())
                            .map(|value| (value as char).to_string()),
                    });
                }
            }
            if !direct_mismatches.is_empty() {
                status = "mismatch".to_string();
            }
        } else if !derived_sequence.is_empty() && !uniprot_sequence.is_empty() {
            let alignment = Self::compute_pairwise_alignment_report(
                &format!("{}_derived", current.transcript_id),
                &derived_sequence,
                None,
                None,
                &entry.entry_id,
                &uniprot_sequence,
                None,
                None,
                PairwiseAlignmentMode::Global,
                1,
                -1,
                -5,
                -1,
            )?
            .report;
            global_alignment = Some(alignment);
            status = "length_mismatch".to_string();
        }
        let variant_feature_evidence = entry
            .features
            .iter()
            .filter_map(|feature| {
                let feature_key = Self::normalize_uniprot_variant_feature_key(feature)?;
                let start_aa = feature.interval.start_aa;
                let end_aa = feature.interval.end_aa;
                let status = match (comparison_mode, start_aa, end_aa) {
                    (_, Some(start), Some(end)) if start == 0 || end < start => {
                        UniprotVariantFeatureEvidenceStatus::OutsideComparableCoverage
                    }
                    (UniprotPeptideComparisonMode::DirectCompare, Some(start), Some(end)) => {
                        let start_idx = start.saturating_sub(1);
                        let end_idx = end.min(derived_sequence.len());
                        if end_idx > derived_sequence.len()
                            || end_idx > uniprot_sequence.len()
                            || start_idx >= end_idx
                        {
                            UniprotVariantFeatureEvidenceStatus::OutsideComparableCoverage
                        } else {
                            let derived_slice = &derived_sequence[start_idx..end_idx];
                            let uniprot_slice = &uniprot_sequence[start_idx..end_idx];
                            let ensembl_slice = ensembl_sequence
                                .as_deref()
                                .and_then(|value| value.get(start_idx..end_idx));
                            if derived_slice != uniprot_slice
                                || ensembl_slice
                                    .map(|value| value != uniprot_slice)
                                    .unwrap_or(false)
                            {
                                UniprotVariantFeatureEvidenceStatus::Supported
                            } else {
                                UniprotVariantFeatureEvidenceStatus::Contradicted
                            }
                        }
                    }
                    (UniprotPeptideComparisonMode::GlobalAlignment, Some(start), Some(end))
                        if start <= current.derived_protein_length_aa
                            && end <= entry.sequence_length =>
                    {
                        UniprotVariantFeatureEvidenceStatus::Unaligned
                    }
                    _ => UniprotVariantFeatureEvidenceStatus::OutsideComparableCoverage,
                };
                let detail = match status {
                    UniprotVariantFeatureEvidenceStatus::Supported => {
                        "GENtle-derived or Ensembl peptide sequence differs from the UniProt reference across this feature interval.".to_string()
                    }
                    UniprotVariantFeatureEvidenceStatus::Contradicted => {
                        "GENtle-derived and available Ensembl peptide sequence match the UniProt reference across this feature interval.".to_string()
                    }
                    UniprotVariantFeatureEvidenceStatus::Unaligned => {
                        "Protein lengths differ, so this sequence-changing feature remains alignment-dependent in v1.".to_string()
                    }
                    UniprotVariantFeatureEvidenceStatus::OutsideComparableCoverage => {
                        "The annotated feature interval falls outside currently comparable peptide coverage.".to_string()
                    }
                };
                Some(UniprotVariantFeatureEvidenceRow {
                    feature_key,
                    feature_note: feature.note.clone(),
                    start_aa,
                    end_aa,
                    status,
                    detail,
                })
            })
            .collect::<Vec<_>>();
        Ok(UniprotEnsemblPeptideCompareRow {
            transcript_id: current.transcript_id.clone(),
            transcript_label: current.transcript_label.clone(),
            transcript_feature_id: current.transcript_feature_id,
            ensembl_entry_id: ensembl_entry.map(|value| value.entry_id.clone()),
            derived_protein_length_aa: current.derived_protein_length_aa,
            uniprot_protein_length_aa: entry.sequence_length,
            ensembl_protein_length_aa: ensembl_entry.map(|value| value.sequence_length),
            comparison_mode,
            init_met_declared,
            direct_mismatches,
            global_alignment,
            variant_feature_evidence,
            status,
            warnings: vec![],
        })
    }

    fn audit_row_status_from_parts(
        link_row: &UniprotEnsemblLinkResolutionRow,
        accounting: &UniprotProjectionTranscriptAccountingRow,
        exon_compare: &UniprotEnsemblExonCompareRow,
        peptide_compare: &UniprotEnsemblPeptideCompareRow,
    ) -> (UniprotProjectionAuditRowStatus, Vec<String>) {
        let mut reasons = vec![];
        if link_row.status != "matched" {
            reasons.push(format!("link_{}", link_row.status));
        }
        if !accounting.translated_nt_divisible_by_3 {
            reasons.push("translated_nt_not_divisible_by_3".to_string());
        }
        let count_matches = accounting.expected_aa_count == accounting.uniprot_aa_count
            || (accounting.init_met_declared
                && (accounting.expected_aa_count + 1 == accounting.uniprot_aa_count
                    || accounting.uniprot_aa_count + 1 == accounting.expected_aa_count));
        if !count_matches {
            reasons.push("aa_count_mismatch".to_string());
        }
        if exon_compare.status != "matched" {
            reasons.push(format!("ensembl_exon_{}", exon_compare.status));
        }
        if peptide_compare.status != "matched" {
            reasons.push(format!("peptide_{}", peptide_compare.status));
        }
        let status = if reasons.is_empty() {
            UniprotProjectionAuditRowStatus::Consistent
        } else if reasons.iter().all(|reason| {
            reason.starts_with("link_missing")
                || reason.starts_with("ensembl_exon_missing_ensembl_evidence")
        }) {
            UniprotProjectionAuditRowStatus::MissingEvidence
        } else if reasons.iter().any(|reason| {
            matches!(
                reason.as_str(),
                "translated_nt_not_divisible_by_3" | "aa_count_mismatch"
            ) || reason.starts_with("peptide_")
                || reason.starts_with("ensembl_exon_mismatch")
        }) {
            UniprotProjectionAuditRowStatus::Mismatch
        } else {
            UniprotProjectionAuditRowStatus::Warning
        };
        (status, reasons)
    }

    fn compose_uniprot_projection_audit_report(
        report_id: &str,
        projection: &UniprotGenomeProjection,
        entry: &UniprotEntry,
        ensembl_entry_id: Option<&str>,
        op_id: Option<String>,
        run_id: Option<String>,
        link_report: &UniprotEnsemblLinkReport,
        accounting_report: &UniprotProjectionTranscriptAccountingReport,
        exon_compare_report: &UniprotEnsemblExonCompareReport,
        peptide_compare_report: &UniprotEnsemblPeptideCompareReport,
    ) -> UniprotProjectionAuditReport {
        let rows = accounting_report
            .rows
            .iter()
            .map(|accounting| {
                let link_row = link_report
                    .rows
                    .iter()
                    .find(|row| row.transcript_id == accounting.transcript_id)
                    .cloned()
                    .unwrap_or_default();
                let exon_compare = exon_compare_report
                    .rows
                    .iter()
                    .find(|row| row.transcript_id == accounting.transcript_id)
                    .cloned();
                let peptide_compare = peptide_compare_report
                    .rows
                    .iter()
                    .find(|row| row.transcript_id == accounting.transcript_id)
                    .cloned();
                let fallback_exon_compare =
                    exon_compare.unwrap_or_else(|| UniprotEnsemblExonCompareRow {
                        transcript_id: accounting.transcript_id.clone(),
                        transcript_label: accounting.transcript_label.clone(),
                        transcript_feature_id: accounting.transcript_feature_id,
                        ensembl_entry_id: ensembl_entry_id.map(|value| value.to_string()),
                        ensembl_transcript_id: None,
                        current_contributing_exons: accounting.contributing_exons.clone(),
                        ensembl_contributing_exons: vec![],
                        matched_exon_ordinals: vec![],
                        missing_in_ensembl_ordinals: accounting
                            .contributing_exons
                            .iter()
                            .map(|exon| exon.ordinal)
                            .collect(),
                        excess_in_ensembl_ordinals: vec![],
                        boundary_differences: vec![],
                        status: "missing_ensembl_evidence".to_string(),
                        warnings: vec![],
                    });
                let fallback_peptide_compare =
                    peptide_compare.unwrap_or_else(|| UniprotEnsemblPeptideCompareRow {
                        transcript_id: accounting.transcript_id.clone(),
                        transcript_label: accounting.transcript_label.clone(),
                        transcript_feature_id: accounting.transcript_feature_id,
                        ensembl_entry_id: ensembl_entry_id.map(|value| value.to_string()),
                        derived_protein_length_aa: accounting.derived_protein_length_aa,
                        uniprot_protein_length_aa: entry.sequence_length,
                        ensembl_protein_length_aa: None,
                        comparison_mode: UniprotPeptideComparisonMode::Unavailable,
                        init_met_declared: accounting.init_met_declared,
                        direct_mismatches: vec![],
                        global_alignment: None,
                        variant_feature_evidence: vec![],
                        status: "missing_ensembl_evidence".to_string(),
                        warnings: vec![],
                    });
                let (status, mismatch_reasons) = Self::audit_row_status_from_parts(
                    &link_row,
                    accounting,
                    &fallback_exon_compare,
                    &fallback_peptide_compare,
                );
                UniprotProjectionAuditRow {
                    transcript_id: accounting.transcript_id.clone(),
                    transcript_label: accounting.transcript_label.clone(),
                    transcript_feature_id: accounting.transcript_feature_id,
                    status,
                    mismatch_reasons,
                    link_resolution: link_row,
                    accounting: accounting.clone(),
                    ensembl_exon_compare: Some(fallback_exon_compare),
                    peptide_compare: Some(fallback_peptide_compare),
                }
            })
            .collect::<Vec<_>>();
        let email_draft = Self::build_uniprot_projection_audit_email_draft(entry, &rows);
        UniprotProjectionAuditReport {
            schema: UNIPROT_PROJECTION_AUDIT_REPORT_SCHEMA.to_string(),
            report_id: report_id.to_string(),
            projection_id: projection.projection_id.clone(),
            entry_id: projection.entry_id.clone(),
            seq_id: projection.seq_id.clone(),
            transcript_id_filter: projection.transcript_id_filter.clone(),
            ensembl_entry_id: ensembl_entry_id
                .map(str::trim)
                .filter(|value| !value.is_empty())
                .map(|value| value.to_string()),
            generated_at_unix_ms: Self::now_unix_ms(),
            op_id,
            run_id,
            rows,
            warnings: link_report
                .warnings
                .iter()
                .chain(accounting_report.warnings.iter())
                .chain(exon_compare_report.warnings.iter())
                .chain(peptide_compare_report.warnings.iter())
                .cloned()
                .collect(),
            maintainer_email_draft: email_draft,
        }
    }

    fn build_uniprot_projection_audit_email_draft(
        entry: &UniprotEntry,
        rows: &[UniprotProjectionAuditRow],
    ) -> Option<UniprotProjectionAuditEmailDraft> {
        let failing = rows
            .iter()
            .filter(|row| row.status != UniprotProjectionAuditRowStatus::Consistent)
            .collect::<Vec<_>>();
        if failing.is_empty() {
            return None;
        }
        let mut body = String::new();
        body.push_str("Dear UniProt maintainers,\n\n");
        body.push_str(&format!(
            "GENtle flagged a presumed inconsistency in the transcript-to-protein mapping for UniProt entry {} ({}) during a local cross-check against projected transcript structure and imported Ensembl evidence. This is an unsent local draft for human review.\n\n",
            entry.entry_id, entry.accession
        ));
        for row in &failing {
            let peptide_mode = row
                .peptide_compare
                .as_ref()
                .map(|value| value.comparison_mode.as_str())
                .unwrap_or("unavailable");
            body.push_str(&format!(
                "- Transcript {} ({}): exon_nt_sum={}, untranslated_5prime_nt={}, untranslated_3prime_nt={}, translated_nt={}, divisible_by_3={}, expected_aa_count={}, uniprot_aa_count={}, comparison_mode={}. Reasons: {}.\n",
                row.transcript_id,
                row.transcript_label,
                row.accounting.contributing_exon_nt_sum,
                row.accounting.untranslated_5prime_nt,
                row.accounting.untranslated_3prime_nt,
                row.accounting.translated_nt,
                row.accounting.translated_nt_divisible_by_3,
                row.accounting.expected_aa_count,
                row.accounting.uniprot_aa_count,
                peptide_mode,
                if row.mismatch_reasons.is_empty() {
                    "none".to_string()
                } else {
                    row.mismatch_reasons.join(", ")
                }
            ));
        }
        body.push_str(
            "\nThis may reflect a presumed inconsistency in transcript projection, exon/CDS selection, or associated protein annotation rather than a definite UniProt error. We are only surfacing the local comparison numbers here for review.\n",
        );
        Some(UniprotProjectionAuditEmailDraft {
            subject: format!(
                "Draft: presumed transcript/protein inconsistency for UniProt entry {}",
                entry.entry_id
            ),
            body,
            failing_transcript_ids: failing
                .iter()
                .map(|row| row.transcript_id.clone())
                .collect(),
        })
    }

    fn build_uniprot_projection_audit_parity_report(
        report_id: &str,
        direct_report: UniprotProjectionAuditReport,
        composed_report: UniprotProjectionAuditReport,
        op_id: Option<String>,
        run_id: Option<String>,
    ) -> UniprotProjectionAuditParityReport {
        let direct_rows = direct_report
            .rows
            .iter()
            .map(|row| (row.transcript_id.clone(), row))
            .collect::<BTreeMap<_, _>>();
        let composed_rows = composed_report
            .rows
            .iter()
            .map(|row| (row.transcript_id.clone(), row))
            .collect::<BTreeMap<_, _>>();
        let transcript_ids = direct_rows
            .keys()
            .chain(composed_rows.keys())
            .cloned()
            .collect::<BTreeSet<_>>();
        let rows = transcript_ids
            .into_iter()
            .map(|transcript_id| {
                let direct = direct_rows.get(&transcript_id);
                let composed = composed_rows.get(&transcript_id);
                let direct_status = direct.map(|row| row.status).unwrap_or_default();
                let composed_status = composed.map(|row| row.status).unwrap_or_default();
                let statuses_match = direct_status == composed_status;
                let accounting_match = match (direct, composed) {
                    (Some(left), Some(right)) => {
                        left.accounting.contributing_exon_nt_sum
                            == right.accounting.contributing_exon_nt_sum
                            && left.accounting.untranslated_5prime_nt
                                == right.accounting.untranslated_5prime_nt
                            && left.accounting.untranslated_3prime_nt
                                == right.accounting.untranslated_3prime_nt
                            && left.accounting.translated_nt == right.accounting.translated_nt
                            && left.accounting.expected_aa_count
                                == right.accounting.expected_aa_count
                    }
                    _ => false,
                };
                let mismatch_reason_match = match (direct, composed) {
                    (Some(left), Some(right)) => {
                        let mut left_reasons = left.mismatch_reasons.clone();
                        let mut right_reasons = right.mismatch_reasons.clone();
                        left_reasons.sort();
                        right_reasons.sort();
                        left_reasons == right_reasons
                    }
                    _ => false,
                };
                let comparison_mode_match = match (direct, composed) {
                    (Some(left), Some(right)) => {
                        left.peptide_compare
                            .as_ref()
                            .map(|value| value.comparison_mode)
                            == right
                                .peptide_compare
                                .as_ref()
                                .map(|value| value.comparison_mode)
                    }
                    _ => false,
                };
                let mut differences = vec![];
                if !statuses_match {
                    differences.push("status".to_string());
                }
                if !accounting_match {
                    differences.push("accounting".to_string());
                }
                if !mismatch_reason_match {
                    differences.push("mismatch_reasons".to_string());
                }
                if !comparison_mode_match {
                    differences.push("comparison_mode".to_string());
                }
                UniprotProjectionAuditParityRow {
                    transcript_id,
                    direct_status,
                    composed_status,
                    statuses_match,
                    accounting_match,
                    mismatch_reason_match,
                    comparison_mode_match,
                    differences,
                }
            })
            .collect::<Vec<_>>();
        let direct_failing = direct_report
            .maintainer_email_draft
            .as_ref()
            .map(|draft| {
                draft
                    .failing_transcript_ids
                    .iter()
                    .cloned()
                    .collect::<BTreeSet<_>>()
            })
            .unwrap_or_default();
        let composed_failing = composed_report
            .maintainer_email_draft
            .as_ref()
            .map(|draft| {
                draft
                    .failing_transcript_ids
                    .iter()
                    .cloned()
                    .collect::<BTreeSet<_>>()
            })
            .unwrap_or_default();
        UniprotProjectionAuditParityReport {
            schema: UNIPROT_PROJECTION_AUDIT_PARITY_REPORT_SCHEMA.to_string(),
            report_id: report_id.to_string(),
            projection_id: direct_report.projection_id.clone(),
            entry_id: direct_report.entry_id.clone(),
            seq_id: direct_report.seq_id.clone(),
            transcript_id_filter: direct_report.transcript_id_filter.clone(),
            ensembl_entry_id: direct_report.ensembl_entry_id.clone(),
            generated_at_unix_ms: Self::now_unix_ms(),
            op_id,
            run_id,
            direct_report,
            composed_report,
            rows,
            email_draft_transcripts_match: direct_failing == composed_failing,
            warnings: vec![],
        }
    }

    pub fn resolve_uniprot_ensembl_links(
        &self,
        projection_id: &str,
        transcript_id: Option<&str>,
    ) -> Result<UniprotEnsemblLinkReport, EngineError> {
        let projection = self.get_uniprot_genome_projection(projection_id)?;
        let entry = self.get_uniprot_entry(&projection.entry_id)?;
        let normalized_filter = transcript_id
            .map(Self::normalize_transcript_probe)
            .filter(|value| !value.is_empty());
        let rows = projection
            .transcript_projections
            .iter()
            .filter(|transcript| {
                normalized_filter.as_ref().is_none_or(|needle| {
                    Self::normalize_transcript_probe(&transcript.transcript_id) == *needle
                })
            })
            .map(|transcript| Self::build_uniprot_link_resolution_row(&entry, transcript))
            .collect::<Vec<_>>();
        Ok(UniprotEnsemblLinkReport {
            schema: UNIPROT_ENSEMBL_LINKS_SCHEMA.to_string(),
            projection_id: projection.projection_id,
            entry_id: projection.entry_id,
            seq_id: projection.seq_id,
            transcript_id_filter: transcript_id
                .map(str::trim)
                .filter(|value| !value.is_empty())
                .map(|value| value.to_string()),
            rows,
            warnings: projection.warnings,
        })
    }

    pub fn summarize_uniprot_projection_transcript_accounting(
        &self,
        projection_id: &str,
        transcript_id: Option<&str>,
    ) -> Result<UniprotProjectionTranscriptAccountingReport, EngineError> {
        let projection = self.get_uniprot_genome_projection(projection_id)?;
        let entry = self.get_uniprot_entry(&projection.entry_id)?;
        let normalized_filter = transcript_id
            .map(Self::normalize_transcript_probe)
            .filter(|value| !value.is_empty());
        let mut rows = vec![];
        let mut warnings = projection.warnings.clone();
        for transcript in projection
            .transcript_projections
            .iter()
            .filter(|transcript| {
                normalized_filter.as_ref().is_none_or(|needle| {
                    Self::normalize_transcript_probe(&transcript.transcript_id) == *needle
                })
            })
        {
            match self.derive_uniprot_projection_transcript_context(&projection.seq_id, transcript)
            {
                Ok(derived) => rows.push(Self::summarize_translation_accounting(
                    &derived.transcript_id,
                    &derived.transcript_label,
                    Some(derived.transcript_feature_id),
                    derived.is_reverse,
                    derived.derivation.as_ref(),
                    &derived.transcript_exons_1based,
                    &entry,
                )),
                Err(err) => warnings.push(format!(
                    "Transcript '{}' could not be prepared for accounting: {}",
                    transcript.transcript_id, err.message
                )),
            }
        }
        Ok(UniprotProjectionTranscriptAccountingReport {
            schema: UNIPROT_PROJECTION_TRANSCRIPT_ACCOUNTING_SCHEMA.to_string(),
            projection_id: projection.projection_id,
            entry_id: projection.entry_id,
            seq_id: projection.seq_id,
            transcript_id_filter: transcript_id
                .map(str::trim)
                .filter(|value| !value.is_empty())
                .map(|value| value.to_string()),
            rows,
            warnings,
        })
    }

    pub fn compare_uniprot_projection_to_ensembl_exons(
        &self,
        projection_id: &str,
        transcript_id: Option<&str>,
        ensembl_entry_id: Option<&str>,
    ) -> Result<UniprotEnsemblExonCompareReport, EngineError> {
        let projection = self.get_uniprot_genome_projection(projection_id)?;
        let entry = self.get_uniprot_entry(&projection.entry_id)?;
        let accounting =
            self.summarize_uniprot_projection_transcript_accounting(projection_id, transcript_id)?;
        let links = self.resolve_uniprot_ensembl_links(projection_id, transcript_id)?;
        let mut rows = vec![];
        let mut warnings = vec![];
        for current in &accounting.rows {
            let link_row = links
                .rows
                .iter()
                .find(|row| row.transcript_id == current.transcript_id)
                .cloned()
                .unwrap_or_default();
            let resolved_entry = match self
                .resolve_ensembl_entry_for_projection_transcript(ensembl_entry_id, &link_row)
            {
                Ok(entry) => entry,
                Err(err) => {
                    warnings.push(err.message.clone());
                    None
                }
            };
            let ensembl_accounting = resolved_entry.as_ref().and_then(|resolved| {
                Self::build_ensembl_accounting_row(&current.transcript_id, resolved, &entry)
            });
            rows.push(Self::build_exon_compare_row(
                current,
                resolved_entry.as_ref(),
                ensembl_accounting.as_ref(),
            ));
        }
        Ok(UniprotEnsemblExonCompareReport {
            schema: UNIPROT_ENSEMBL_EXON_COMPARE_SCHEMA.to_string(),
            projection_id: projection.projection_id,
            entry_id: projection.entry_id,
            seq_id: projection.seq_id,
            transcript_id_filter: transcript_id
                .map(str::trim)
                .filter(|value| !value.is_empty())
                .map(|value| value.to_string()),
            ensembl_entry_id: ensembl_entry_id
                .map(str::trim)
                .filter(|value| !value.is_empty())
                .map(|value| value.to_string()),
            rows,
            warnings,
        })
    }

    pub fn compare_uniprot_variant_features_to_ensembl_peptide(
        &self,
        projection_id: &str,
        transcript_id: Option<&str>,
        ensembl_entry_id: Option<&str>,
    ) -> Result<UniprotEnsemblPeptideCompareReport, EngineError> {
        let projection = self.get_uniprot_genome_projection(projection_id)?;
        let entry = self.get_uniprot_entry(&projection.entry_id)?;
        let accounting =
            self.summarize_uniprot_projection_transcript_accounting(projection_id, transcript_id)?;
        let links = self.resolve_uniprot_ensembl_links(projection_id, transcript_id)?;
        let mut rows = vec![];
        let mut warnings = vec![];
        for current in &accounting.rows {
            let link_row = links
                .rows
                .iter()
                .find(|row| row.transcript_id == current.transcript_id)
                .cloned()
                .unwrap_or_default();
            let resolved_entry = match self
                .resolve_ensembl_entry_for_projection_transcript(ensembl_entry_id, &link_row)
            {
                Ok(entry) => entry,
                Err(err) => {
                    warnings.push(err.message.clone());
                    None
                }
            };
            rows.push(self.build_peptide_compare_row(current, resolved_entry.as_ref(), &entry)?);
        }
        Ok(UniprotEnsemblPeptideCompareReport {
            schema: UNIPROT_ENSEMBL_PEPTIDE_COMPARE_SCHEMA.to_string(),
            projection_id: projection.projection_id,
            entry_id: projection.entry_id,
            seq_id: projection.seq_id,
            transcript_id_filter: transcript_id
                .map(str::trim)
                .filter(|value| !value.is_empty())
                .map(|value| value.to_string()),
            ensembl_entry_id: ensembl_entry_id
                .map(str::trim)
                .filter(|value| !value.is_empty())
                .map(|value| value.to_string()),
            rows,
            warnings,
        })
    }

    pub fn audit_uniprot_projection_consistency_from_primitives(
        &self,
        projection_id: &str,
        transcript_id: Option<&str>,
        report_id: &str,
        ensembl_entry_id: Option<&str>,
        op_id: Option<String>,
        run_id: Option<String>,
    ) -> Result<UniprotProjectionAuditReport, EngineError> {
        let projection = self.get_uniprot_genome_projection(projection_id)?;
        let entry = self.get_uniprot_entry(&projection.entry_id)?;
        let links = self.resolve_uniprot_ensembl_links(projection_id, transcript_id)?;
        let accounting =
            self.summarize_uniprot_projection_transcript_accounting(projection_id, transcript_id)?;
        let exon_compare = self.compare_uniprot_projection_to_ensembl_exons(
            projection_id,
            transcript_id,
            ensembl_entry_id,
        )?;
        let peptide_compare = self.compare_uniprot_variant_features_to_ensembl_peptide(
            projection_id,
            transcript_id,
            ensembl_entry_id,
        )?;
        Ok(Self::compose_uniprot_projection_audit_report(
            report_id,
            &projection,
            &entry,
            ensembl_entry_id,
            op_id,
            run_id,
            &links,
            &accounting,
            &exon_compare,
            &peptide_compare,
        ))
    }

    pub fn audit_uniprot_projection_parity_from_primitives(
        &self,
        projection_id: &str,
        transcript_id: Option<&str>,
        report_id: &str,
        ensembl_entry_id: Option<&str>,
        op_id: Option<String>,
        run_id: Option<String>,
    ) -> Result<UniprotProjectionAuditParityReport, EngineError> {
        let direct_report = self.audit_uniprot_projection_consistency_from_primitives(
            projection_id,
            transcript_id,
            report_id,
            ensembl_entry_id,
            op_id.clone(),
            run_id.clone(),
        )?;
        let projection = self.get_uniprot_genome_projection(projection_id)?;
        let entry = self.get_uniprot_entry(&projection.entry_id)?;
        let links = self.resolve_uniprot_ensembl_links(projection_id, transcript_id)?;
        let accounting =
            self.summarize_uniprot_projection_transcript_accounting(projection_id, transcript_id)?;
        let exon_compare = self.compare_uniprot_projection_to_ensembl_exons(
            projection_id,
            transcript_id,
            ensembl_entry_id,
        )?;
        let peptide_compare = self.compare_uniprot_variant_features_to_ensembl_peptide(
            projection_id,
            transcript_id,
            ensembl_entry_id,
        )?;
        let composed_report = Self::compose_uniprot_projection_audit_report(
            report_id,
            &projection,
            &entry,
            ensembl_entry_id,
            op_id.clone(),
            run_id.clone(),
            &links,
            &accounting,
            &exon_compare,
            &peptide_compare,
        );
        Ok(Self::build_uniprot_projection_audit_parity_report(
            report_id,
            direct_report,
            composed_report,
            op_id,
            run_id,
        ))
    }

    fn uniprot_feature_query_profile_species_label(
        profile: TranslationSpeedProfile,
    ) -> (&'static str, Option<&'static str>) {
        Self::translation_speed_profile_reference_species(profile)
    }

    fn build_uniprot_feature_query_optimized_dna(
        amino_acid_sequence: &str,
        profile: TranslationSpeedProfile,
    ) -> (String, Vec<String>) {
        let (species_label, profile_warning) =
            Self::uniprot_feature_query_profile_species_label(profile);
        let mut warnings = vec![];
        if let Some(warning) = profile_warning {
            warnings.push(warning.to_string());
        }
        let mut out = String::new();
        for (idx, residue) in amino_acid_sequence.chars().enumerate() {
            let aa = if residue == '*' {
                STOP_CODON
            } else {
                residue.to_ascii_uppercase()
            };
            let codon = AMINO_ACIDS
                .preferred_species_codon(aa, species_label)
                .or_else(|| AMINO_ACIDS.preferred_species_codon(aa, "Default"))
                .or_else(|| {
                    let mut candidates = AMINO_ACIDS.aa2codons(aa, Some(1));
                    candidates.sort();
                    candidates
                        .into_iter()
                        .next()
                        .map(|codon| String::from_utf8_lossy(&codon).to_string())
                })
                .unwrap_or_else(|| {
                    warnings.push(format!(
                        "Protein residue {} ('{}') has no bundled preferred codon; used 'NNN'.",
                        idx + 1,
                        residue
                    ));
                    "NNN".to_string()
                });
            out.push_str(&codon);
        }
        (out, warnings)
    }

    fn uniprot_feature_projection_match_fields(
        feature: &UniprotFeatureProjection,
        query_lower: &str,
    ) -> Vec<String> {
        let mut fields = vec![];
        if feature
            .feature_key
            .to_ascii_lowercase()
            .contains(query_lower)
        {
            fields.push("feature_key".to_string());
        }
        if feature
            .feature_note
            .as_deref()
            .map(|note| note.to_ascii_lowercase().contains(query_lower))
            .unwrap_or(false)
        {
            fields.push("feature_note".to_string());
        }
        fields
    }

    fn uniprot_feature_projection_label(feature: &UniprotFeatureProjection) -> String {
        match feature
            .feature_note
            .as_deref()
            .map(str::trim)
            .filter(|value| !value.is_empty())
        {
            Some(note) => format!(
                "{} '{}' (aa {}..{})",
                feature.feature_key, note, feature.aa_start, feature.aa_end
            ),
            None => format!(
                "{} (aa {}..{})",
                feature.feature_key, feature.aa_start, feature.aa_end
            ),
        }
    }

    fn extract_uniprot_feature_amino_acid_sequence(
        protein_sequence: &str,
        aa_start: usize,
        aa_end: usize,
    ) -> String {
        if aa_start == 0 || aa_end < aa_start {
            return String::new();
        }
        protein_sequence
            .chars()
            .skip(aa_start.saturating_sub(1))
            .take(aa_end.saturating_sub(aa_start).saturating_add(1))
            .collect()
    }

    fn resolve_uniprot_projection_transcript_feature<'a>(
        seq: &'a DNAsequence,
        projection: &UniprotTranscriptProjection,
    ) -> Option<(usize, &'a gb_io::seq::Feature)> {
        let normalized = Self::normalize_transcript_probe(&projection.transcript_id);
        if let Some(feature_id) = projection.transcript_feature_id
            && let Some(feature) = seq.features().get(feature_id)
            && Self::is_mrna_feature(feature)
            && Self::feature_transcript_match_keys(feature, feature_id)
                .iter()
                .any(|candidate| candidate == &normalized)
        {
            return Some((feature_id, feature));
        }
        seq.features()
            .iter()
            .enumerate()
            .find(|(feature_id, feature)| {
                Self::is_mrna_feature(feature)
                    && Self::feature_transcript_match_keys(feature, *feature_id)
                        .iter()
                        .any(|candidate| candidate == &normalized)
            })
    }

    fn transcript_feature_exon_ranges_1based(feature: &gb_io::seq::Feature) -> Vec<(usize, usize)> {
        let mut ranges = vec![];
        collect_location_ranges_usize(&feature.location, &mut ranges);
        ranges.sort_unstable_by(|a, b| a.0.cmp(&b.0).then(a.1.cmp(&b.1)));
        ranges
            .into_iter()
            .filter(|(start, end)| end > start)
            .map(|(start, end)| (start.saturating_add(1), end))
            .collect()
    }

    fn transcript_feature_cds_ranges_1based(feature: &gb_io::seq::Feature) -> Vec<(usize, usize)> {
        Self::feature_qualifier_ranges_0based(feature, "cds_ranges_1based")
            .into_iter()
            .map(|(start, end)| (start.saturating_add(1), end))
            .collect()
    }

    fn transcript_projection_coding_exon_ranges_1based(
        projection: &UniprotTranscriptProjection,
    ) -> Vec<(usize, usize)> {
        let mut ranges = projection
            .aa_segments
            .iter()
            .map(|segment| {
                (
                    segment.genomic_start_1based.min(segment.genomic_end_1based),
                    segment.genomic_start_1based.max(segment.genomic_end_1based),
                )
            })
            .filter(|(start, end)| *start > 0 && *end >= *start)
            .collect::<Vec<_>>();
        ranges.sort_unstable_by(|a, b| a.0.cmp(&b.0).then(a.1.cmp(&b.1)));
        ranges.dedup();
        ranges
    }

    fn order_ranges_1based_for_transcript(
        mut ranges: Vec<(usize, usize)>,
        is_reverse: bool,
    ) -> Vec<(usize, usize)> {
        ranges.sort_unstable_by(|a, b| a.0.cmp(&b.0).then(a.1.cmp(&b.1)));
        if is_reverse {
            ranges.reverse();
        }
        ranges
    }

    fn build_uniprot_feature_exon_spans(
        exon_ranges_1based: &[(usize, usize)],
        genomic_segments: &[UniprotFeatureCodingDnaSegment],
    ) -> Vec<UniprotFeatureCodingDnaExonSpan> {
        let mut by_exon: BTreeMap<usize, UniprotFeatureCodingDnaExonSpan> = BTreeMap::new();
        for segment in genomic_segments {
            let seg_start = segment.genomic_start_1based.min(segment.genomic_end_1based);
            let seg_end = segment.genomic_start_1based.max(segment.genomic_end_1based);
            for (idx, (exon_start, exon_end)) in exon_ranges_1based.iter().enumerate() {
                let overlap_start = seg_start.max(*exon_start);
                let overlap_end = seg_end.min(*exon_end);
                if overlap_end < overlap_start {
                    continue;
                }
                by_exon
                    .entry(idx + 1)
                    .and_modify(|row| {
                        row.coding_start_1based = row.coding_start_1based.min(overlap_start);
                        row.coding_end_1based = row.coding_end_1based.max(overlap_end);
                    })
                    .or_insert_with(|| UniprotFeatureCodingDnaExonSpan {
                        exon_ordinal: idx + 1,
                        exon_start_1based: *exon_start,
                        exon_end_1based: *exon_end,
                        coding_start_1based: overlap_start,
                        coding_end_1based: overlap_end,
                    });
            }
        }
        by_exon.into_values().collect()
    }

    fn build_uniprot_feature_exon_pairs(
        exon_spans: &[UniprotFeatureCodingDnaExonSpan],
    ) -> Vec<UniprotFeatureCodingDnaExonPair> {
        exon_spans
            .windows(2)
            .filter_map(|window| {
                let from = window.first()?;
                let to = window.get(1)?;
                Some(UniprotFeatureCodingDnaExonPair {
                    from_exon_ordinal: from.exon_ordinal,
                    to_exon_ordinal: to.exon_ordinal,
                })
            })
            .collect()
    }

    pub fn query_uniprot_feature_coding_dna(
        &self,
        projection_id: &str,
        feature_query: &str,
        transcript_filter: Option<&str>,
        query_mode: UniprotFeatureCodingDnaQueryMode,
        translation_speed_profile: Option<TranslationSpeedProfile>,
    ) -> Result<UniprotFeatureCodingDnaQueryReport, EngineError> {
        let projection = self.get_uniprot_genome_projection(projection_id)?;
        let entry = self.get_uniprot_entry(&projection.entry_id)?;
        let seq = self
            .state
            .sequences
            .get(&projection.seq_id)
            .ok_or_else(|| EngineError {
                code: ErrorCode::NotFound,
                message: format!("Sequence '{}' not found", projection.seq_id),
            })?;
        let feature_query = feature_query.trim();
        if feature_query.is_empty() {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: "Protein feature query cannot be empty".to_string(),
            });
        }
        let transcript_filter = transcript_filter
            .map(str::trim)
            .filter(|value| !value.is_empty())
            .map(|value| value.to_string());
        let normalized_transcript_filter = transcript_filter
            .as_deref()
            .map(Self::normalize_transcript_probe)
            .filter(|value| !value.is_empty());
        let query_lower = feature_query.to_ascii_lowercase();
        let resolved_profile = match query_mode {
            UniprotFeatureCodingDnaQueryMode::GenomicAsEncoded => None,
            _ => translation_speed_profile
                .map(|profile| {
                    Self::translation_speed_profile_resolution_from_profile(
                        profile,
                        TranslationSpeedProfileSource::ExplicitRequest,
                        None,
                    )
                })
                .or_else(|| {
                    Self::infer_translation_speed_profile_from_organism(entry.organism.as_deref())
                }),
        };

        let mut warnings = projection.warnings.clone();
        if let Some(resolved_profile) = resolved_profile.as_ref() {
            warnings.extend(resolved_profile.warnings.iter().cloned());
        }
        if !matches!(
            query_mode,
            UniprotFeatureCodingDnaQueryMode::GenomicAsEncoded
        ) && resolved_profile.is_none()
        {
            warnings.push(
                "Could not resolve a translation-speed profile automatically; translation_speed_optimized_dna will be omitted unless an explicit profile is supplied."
                    .to_string(),
            );
        }

        let mut matches: Vec<UniprotFeatureCodingDnaMatch> = vec![];
        let mut available_labels = BTreeSet::new();
        for transcript in &projection.transcript_projections {
            if let Some(filter) = normalized_transcript_filter.as_deref()
                && Self::normalize_transcript_probe(&transcript.transcript_id) != filter
            {
                continue;
            }

            let transcript_feature =
                Self::resolve_uniprot_projection_transcript_feature(seq, transcript);
            let is_reverse = transcript.strand.trim() == "-";
            let ordered_exon_ranges_1based = if let Some((_, feature)) = transcript_feature {
                let transcript_exon_ranges = Self::transcript_feature_exon_ranges_1based(feature);
                let transcript_cds_ranges = Self::transcript_feature_cds_ranges_1based(feature);
                let preferred_ranges =
                    if transcript_exon_ranges.len() <= 1 && transcript_cds_ranges.len() > 1 {
                        transcript_cds_ranges
                    } else {
                        transcript_exon_ranges
                    };
                Self::order_ranges_1based_for_transcript(preferred_ranges, is_reverse)
            } else {
                Self::order_ranges_1based_for_transcript(
                    Self::transcript_projection_coding_exon_ranges_1based(transcript),
                    is_reverse,
                )
            };

            for feature in &transcript.feature_projections {
                available_labels.insert(Self::uniprot_feature_projection_label(feature));
                let matched_feature_fields =
                    Self::uniprot_feature_projection_match_fields(feature, &query_lower);
                if matched_feature_fields.is_empty() {
                    continue;
                }

                let amino_acid_sequence = Self::extract_uniprot_feature_amino_acid_sequence(
                    &entry.sequence,
                    feature.aa_start,
                    feature.aa_end,
                );
                let mut match_warnings = transcript.warnings.clone();
                if amino_acid_sequence.is_empty() {
                    match_warnings.push(format!(
                        "Feature aa span {}..{} could not be sliced from UniProt entry '{}'.",
                        feature.aa_start, feature.aa_end, entry.entry_id
                    ));
                }

                let mut genomic_segments = vec![];
                let mut genomic_coding_dna = String::new();
                for segment in &feature.genomic_segments {
                    let seg_start = segment.genomic_start_1based.min(segment.genomic_end_1based);
                    let seg_end = segment.genomic_start_1based.max(segment.genomic_end_1based);
                    let raw = Self::sequence_slice_upper(seq, seg_start.saturating_sub(1), seg_end);
                    if raw.is_empty() {
                        match_warnings.push(format!(
                            "Could not extract genomic DNA for transcript '{}' segment {}..{}.",
                            transcript.transcript_id, seg_start, seg_end
                        ));
                        continue;
                    }
                    let coding_text = if is_reverse {
                        String::from_utf8_lossy(&Self::reverse_complement_bytes(raw.as_bytes()))
                            .to_string()
                    } else {
                        raw
                    };
                    genomic_coding_dna.push_str(&coding_text);
                    genomic_segments.push(UniprotFeatureCodingDnaSegment {
                        aa_start: segment.aa_start,
                        aa_end: segment.aa_end,
                        genomic_start_1based: seg_start,
                        genomic_end_1based: seg_end,
                        strand: segment.strand.clone(),
                    });
                }

                let exon_spans = Self::build_uniprot_feature_exon_spans(
                    &ordered_exon_ranges_1based,
                    &genomic_segments,
                );
                if exon_spans.is_empty() {
                    match_warnings.push(
                        "Could not determine exon attribution for this feature span.".to_string(),
                    );
                } else if transcript_feature.is_none() {
                    match_warnings.push(
                        "Transcript feature could not be resolved back in the current sequence; exon attribution falls back to CDS segments from the stored projection."
                            .to_string(),
                    );
                }
                let exon_pairs = Self::build_uniprot_feature_exon_pairs(&exon_spans);
                let primary_exon_ordinal =
                    (exon_spans.len() == 1).then_some(exon_spans[0].exon_ordinal);
                let primary_exon_pair = (exon_spans.len() == 2)
                    .then(|| exon_pairs.first().cloned())
                    .flatten();

                let mut optimized_warnings = vec![];
                let translation_speed_optimized_dna = if !matches!(
                    query_mode,
                    UniprotFeatureCodingDnaQueryMode::GenomicAsEncoded
                ) {
                    resolved_profile.as_ref().map(|resolution| {
                        let (optimized, local_warnings) =
                            Self::build_uniprot_feature_query_optimized_dna(
                                &amino_acid_sequence,
                                resolution.profile,
                            );
                        optimized_warnings = local_warnings;
                        optimized
                    })
                } else {
                    None
                };
                match_warnings.extend(optimized_warnings);

                matches.push(UniprotFeatureCodingDnaMatch {
                    transcript_id: transcript.transcript_id.clone(),
                    transcript_feature_id: transcript.transcript_feature_id,
                    strand: transcript.strand.clone(),
                    feature_key: feature.feature_key.clone(),
                    feature_note: feature.feature_note.clone(),
                    matched_feature_fields,
                    aa_start: feature.aa_start,
                    aa_end: feature.aa_end,
                    amino_acid_sequence,
                    genomic_segments,
                    genomic_coding_dna,
                    translation_speed_optimized_dna,
                    exon_spans,
                    exon_pairs,
                    primary_exon_ordinal,
                    primary_exon_pair,
                    crosses_exon_junction: feature.genomic_segments.len() > 1,
                    warnings: match_warnings,
                });
            }
        }

        if matches.is_empty() {
            let available = available_labels
                .into_iter()
                .take(8)
                .collect::<Vec<_>>()
                .join(", ");
            let transcript_hint = transcript_filter
                .as_deref()
                .map(|value| format!(" in transcript '{}'", value))
                .unwrap_or_default();
            let suffix = if available.is_empty() {
                "Projection currently exposes no mapped UniProt feature spans.".to_string()
            } else {
                format!("Available mapped features include: {available}.")
            };
            return Err(EngineError {
                code: ErrorCode::NotFound,
                message: format!(
                    "No mapped UniProt feature matched query '{}' in projection '{}'{}; {}",
                    feature_query, projection.projection_id, transcript_hint, suffix
                ),
            });
        }

        Ok(UniprotFeatureCodingDnaQueryReport {
            schema: UNIPROT_FEATURE_CODING_DNA_QUERY_SCHEMA.to_string(),
            projection_id: projection.projection_id,
            entry_id: projection.entry_id,
            seq_id: projection.seq_id,
            feature_query: feature_query.to_string(),
            transcript_filter,
            query_mode,
            requested_translation_speed_profile: translation_speed_profile,
            resolved_translation_speed_profile: resolved_profile.as_ref().map(|r| r.profile),
            resolved_translation_speed_profile_source: resolved_profile.as_ref().map(|r| r.source),
            resolved_translation_speed_reference_species: resolved_profile
                .as_ref()
                .map(|r| r.reference_species.clone()),
            match_count: matches.len(),
            matches,
            warnings,
        })
    }

    fn ensembl_rest_base_url() -> String {
        std::env::var(ENSEMBL_REST_BASE_ENV_VAR)
            .ok()
            .map(|value| value.trim().trim_end_matches('/').to_string())
            .filter(|value| !value.is_empty())
            .unwrap_or_else(|| DEFAULT_ENSEMBL_REST_ENDPOINT.to_string())
    }

    fn fetch_ensembl_rest_text(url: &str, label: &str, query: &str) -> Result<String, EngineError> {
        let client = reqwest::blocking::Client::builder()
            .timeout(Duration::from_secs(45))
            .build()
            .map_err(|e| EngineError {
                code: ErrorCode::Io,
                message: format!("Could not create Ensembl HTTP client: {e}"),
            })?;
        let response = client.get(url).send().map_err(|e| EngineError {
            code: ErrorCode::Io,
            message: format!("Could not fetch Ensembl {label} for '{query}': {e}"),
        })?;
        if !response.status().is_success() {
            let status = response.status();
            let detail = if status == reqwest::StatusCode::NOT_FOUND {
                format!(
                    "Ensembl {label} for '{}' was not found (HTTP {}) at '{}'",
                    query, status, url
                )
            } else {
                format!(
                    "Ensembl {label} for '{}' returned HTTP status {} at '{}'",
                    query, status, url
                )
            };
            return Err(EngineError {
                code: ErrorCode::NotFound,
                message: detail,
            });
        }
        response.text().map_err(|e| EngineError {
            code: ErrorCode::Io,
            message: format!("Could not read Ensembl {label} response for '{query}': {e}"),
        })
    }

    pub(super) fn fetch_ensembl_protein_entry_from_rest(
        query: &str,
        entry_id_override: Option<&str>,
    ) -> Result<EnsemblProteinEntry, EngineError> {
        let query = query.trim();
        if query.is_empty() {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: "Ensembl protein query cannot be empty".to_string(),
            });
        }
        let resolved_query = resolve_ensembl_query(query).map_err(|e| EngineError {
            code: ErrorCode::InvalidInput,
            message: e,
        })?;
        let base_url = Self::ensembl_rest_base_url();

        let (protein_lookup_source_url, protein_lookup_json, transcript_id) = match resolved_query
            .kind
        {
            crate::ensembl_protein::EnsemblProteinQueryKind::Protein => {
                let lookup_url = format!(
                    "{base_url}/lookup/id/{}?content-type=application/json",
                    resolved_query.normalized_query
                );
                let lookup_json =
                    Self::fetch_ensembl_rest_text(&lookup_url, "protein lookup", query)?;
                let lookup = crate::ensembl_protein::parse_translation_lookup_json(&lookup_json)
                    .map_err(|e| EngineError {
                        code: ErrorCode::InvalidInput,
                        message: format!(
                            "Could not parse Ensembl protein lookup for '{query}': {e}"
                        ),
                    })?;
                (
                    Some(lookup_url),
                    Some(lookup_json),
                    lookup.transcript_id.clone(),
                )
            }
            crate::ensembl_protein::EnsemblProteinQueryKind::Transcript => {
                (None, None, resolved_query.normalized_query.clone())
            }
        };

        let transcript_lookup_source_url =
            format!("{base_url}/lookup/id/{transcript_id}?content-type=application/json;expand=1");
        let transcript_lookup_json = Self::fetch_ensembl_rest_text(
            &transcript_lookup_source_url,
            "transcript lookup",
            query,
        )?;
        let transcript_lookup = crate::ensembl_protein::parse_transcript_lookup_json(
            &transcript_lookup_json,
        )
        .map_err(|e| EngineError {
            code: ErrorCode::InvalidInput,
            message: format!("Could not parse Ensembl transcript lookup for '{query}': {e}"),
        })?;
        let protein_id = transcript_lookup
            .translation_id
            .clone()
            .ok_or_else(|| EngineError {
                code: ErrorCode::InvalidInput,
                message: format!(
                    "Transcript '{}' from Ensembl query '{}' has no translated protein stable ID",
                    transcript_lookup.transcript_id, query
                ),
            })?;

        let sequence_source_url = format!(
            "{base_url}/sequence/id/{protein_id}?type=protein;content-type=application/json"
        );
        let sequence_json =
            Self::fetch_ensembl_rest_text(&sequence_source_url, "protein sequence", query)?;
        let feature_source_url = format!(
            "{base_url}/overlap/translation/{protein_id}?feature=protein_feature;content-type=application/json"
        );
        let feature_json =
            Self::fetch_ensembl_rest_text(&feature_source_url, "protein features", query)?;

        build_entry_from_rest_payloads(
            query,
            &transcript_lookup_source_url,
            &transcript_lookup_json,
            protein_lookup_source_url.as_deref(),
            protein_lookup_json.as_deref(),
            &sequence_source_url,
            &sequence_json,
            &feature_source_url,
            &feature_json,
            entry_id_override,
        )
        .map_err(|e| EngineError {
            code: ErrorCode::InvalidInput,
            message: format!("Could not assemble Ensembl protein entry for '{query}': {e}"),
        })
    }

    pub(super) fn import_ensembl_protein_entry_sequence(
        &mut self,
        result: &mut OpResult,
        entry_id: &str,
        output_id: Option<&str>,
    ) -> Result<SeqId, EngineError> {
        let entry = self.get_ensembl_protein_entry(entry_id)?;
        let base_seq_id = output_id
            .map(str::trim)
            .filter(|value| !value.is_empty())
            .map(|value| value.to_string())
            .unwrap_or_else(|| {
                let normalized = normalize_ensembl_protein_entry_id(&entry.entry_id);
                if normalized.is_empty() {
                    format!("ensembl_{}", entry.protein_id)
                } else {
                    normalized
                }
            });
        let seq_id = self.unique_seq_id(&base_seq_id);
        let mut protein = DNAsequence::from_sequence(&entry.sequence).map_err(|e| EngineError {
            code: ErrorCode::Internal,
            message: format!(
                "Could not construct protein sequence for Ensembl protein entry '{}': {e}",
                entry.entry_id
            ),
        })?;
        protein.set_name(
            entry
                .transcript_display_name
                .clone()
                .or_else(|| entry.gene_symbol.clone())
                .unwrap_or_else(|| entry.entry_id.clone()),
        );
        protein.set_molecule_type("protein");
        if protein.len() > 0 {
            let mut protein_qualifiers = vec![
                ("entry_id".into(), Some(entry.entry_id.clone())),
                ("protein_id".into(), Some(entry.protein_id.clone())),
                ("transcript_id".into(), Some(entry.transcript_id.clone())),
                (
                    "sequence_length_aa".into(),
                    Some(entry.sequence_length.to_string()),
                ),
                (
                    "synthetic_origin".into(),
                    Some("ensembl_protein_import".to_string()),
                ),
            ];
            if let Some(gene_id) = entry.gene_id.as_ref() {
                protein_qualifiers.push(("gene_id".into(), Some(gene_id.clone())));
            }
            if let Some(gene_symbol) = entry.gene_symbol.as_ref() {
                protein_qualifiers.push(("gene".into(), Some(gene_symbol.clone())));
                protein_qualifiers.push(("label".into(), Some(gene_symbol.clone())));
            }
            if let Some(display_name) = entry.transcript_display_name.as_ref() {
                protein_qualifiers.push(("product".into(), Some(display_name.clone())));
            }
            if let Some(species) = entry.species.as_ref() {
                protein_qualifiers.push(("organism".into(), Some(species.clone())));
            }
            let protein_len = protein.len();
            protein.features_mut().push(gb_io::seq::Feature {
                kind: "Protein".into(),
                location: gb_io::seq::Location::simple_range(0, protein_len as i64),
                qualifiers: protein_qualifiers,
            });
        }
        for feature in &entry.features {
            let (Some(start_aa), Some(end_aa)) = (feature.start_aa, feature.end_aa) else {
                continue;
            };
            if start_aa == 0 || end_aa < start_aa || end_aa > protein.len() {
                continue;
            }
            let mut qualifiers = vec![
                ("label".into(), Some(feature.feature_key.clone())),
                ("feature_type".into(), Some(feature.feature_type.clone())),
                ("entry_id".into(), Some(entry.entry_id.clone())),
                ("protein_id".into(), Some(entry.protein_id.clone())),
                ("transcript_id".into(), Some(entry.transcript_id.clone())),
            ];
            if let Some(description) = feature.description.as_ref() {
                qualifiers.push(("note".into(), Some(description.clone())));
            }
            if let Some(interpro_id) = feature.interpro_id.as_ref() {
                qualifiers.push(("interpro".into(), Some(interpro_id.clone())));
            }
            for (key, value) in &feature.qualifiers {
                qualifiers.push((key.clone().into(), Some(value.clone())));
            }
            protein.features_mut().push(gb_io::seq::Feature {
                kind: feature.feature_type.clone().into(),
                location: gb_io::seq::Location::simple_range(
                    start_aa.saturating_sub(1) as i64,
                    end_aa as i64,
                ),
                qualifiers,
            });
        }
        Self::prepare_sequence(&mut protein);
        self.state.sequences.insert(seq_id.clone(), protein);
        self.add_lineage_node(&seq_id, SequenceOrigin::Derived, Some(&result.op_id));
        result.created_seq_ids.push(seq_id.clone());
        result.messages.push(format!(
            "Imported Ensembl protein '{}' (transcript '{}') as '{}' ({} aa, features={}).",
            entry.protein_id,
            entry.transcript_id,
            seq_id,
            entry.sequence_length,
            entry.features.len()
        ));
        Ok(seq_id)
    }

    pub(super) fn fetch_uniprot_swiss_prot_text(
        query: &str,
    ) -> Result<(String, String), EngineError> {
        let query = query.trim();
        if query.is_empty() {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: "UniProt query cannot be empty".to_string(),
            });
        }
        let query_url = query.replace(' ', "%20");
        let url = format!("https://rest.uniprot.org/uniprotkb/{query_url}.txt");
        let client = reqwest::blocking::Client::builder()
            .timeout(Duration::from_secs(45))
            .build()
            .map_err(|e| EngineError {
                code: ErrorCode::Io,
                message: format!("Could not create UniProt HTTP client: {e}"),
            })?;
        let response = client.get(&url).send().map_err(|e| EngineError {
            code: ErrorCode::Io,
            message: format!("Could not fetch UniProt entry '{query}': {e}"),
        })?;
        if !response.status().is_success() {
            let status = response.status();
            let detail = if status == reqwest::StatusCode::NOT_FOUND {
                format!(
                    "UniProt query '{}' was not found (HTTP {}) at '{}'; verify accession/entry ID",
                    query, status, url
                )
            } else {
                format!(
                    "UniProt query '{}' returned HTTP status {} at '{}'",
                    query, status, url
                )
            };
            return Err(EngineError {
                code: ErrorCode::NotFound,
                message: detail,
            });
        }
        let text = response.text().map_err(|e| EngineError {
            code: ErrorCode::Io,
            message: format!("Could not read UniProt response body for '{query}': {e}"),
        })?;
        Ok((url, text))
    }

    pub(super) fn fetch_genbank_accession_text(
        accession: &str,
    ) -> Result<(String, String), EngineError> {
        let accession = validate_genbank_accession(accession).map_err(|e| EngineError {
            code: ErrorCode::InvalidInput,
            message: e,
        })?;
        let source_url = build_genbank_efetch_url(&accession, "gbwithparts");
        if let Some(path) = source_url.strip_prefix("file://") {
            let text = std::fs::read_to_string(path).map_err(|e| EngineError {
                code: ErrorCode::Io,
                message: format!(
                    "Could not read GenBank source file '{}' for accession '{}': {e}",
                    path, accession
                ),
            })?;
            return Ok((source_url, text));
        }
        let client = reqwest::blocking::Client::builder()
            .timeout(Duration::from_secs(45))
            .build()
            .map_err(|e| EngineError {
                code: ErrorCode::Io,
                message: format!("Could not create GenBank HTTP client: {e}"),
            })?;
        let response = client.get(&source_url).send().map_err(|e| EngineError {
            code: ErrorCode::Io,
            message: format!(
                "Could not fetch GenBank accession '{}' from '{}': {e}",
                accession, source_url
            ),
        })?;
        if !response.status().is_success() {
            let status = response.status();
            let detail = if status == reqwest::StatusCode::NOT_FOUND {
                format!(
                    "GenBank accession '{}' was not found (HTTP {}) at '{}'",
                    accession, status, source_url
                )
            } else {
                format!(
                    "GenBank accession '{}' returned HTTP status {} at '{}'",
                    accession, status, source_url
                )
            };
            return Err(EngineError {
                code: ErrorCode::NotFound,
                message: detail,
            });
        }
        let text = response.text().map_err(|e| EngineError {
            code: ErrorCode::Io,
            message: format!(
                "Could not read GenBank response body for accession '{}': {e}",
                accession
            ),
        })?;
        Ok((source_url, text))
    }

    pub(super) fn parse_uniprot_entry_text(
        text: &str,
        source: &str,
        source_query: Option<&str>,
        entry_id_override: Option<&str>,
    ) -> Result<UniprotEntry, EngineError> {
        parse_swiss_prot_text(text, source, source_query, entry_id_override).map_err(|e| {
            EngineError {
                code: ErrorCode::InvalidInput,
                message: format!("Could not parse SWISS-PROT text: {e}"),
            }
        })
    }

    pub(super) fn import_uniprot_entry_sequence(
        &mut self,
        result: &mut OpResult,
        entry_id: &str,
        output_id: Option<&str>,
    ) -> Result<SeqId, EngineError> {
        let entry = self.get_uniprot_entry(entry_id)?;
        let base_seq_id = output_id
            .map(str::trim)
            .filter(|value| !value.is_empty())
            .map(|value| value.to_string())
            .unwrap_or_else(|| {
                let normalized = crate::uniprot::normalize_entry_id(&entry.entry_id);
                if normalized.is_empty() {
                    format!("uniprot_{}", entry.primary_id)
                } else {
                    normalized
                }
            });
        let seq_id = self.unique_seq_id(&base_seq_id);
        let mut protein = DNAsequence::from_sequence(&entry.sequence).map_err(|e| EngineError {
            code: ErrorCode::Internal,
            message: format!(
                "Could not construct protein sequence for UniProt entry '{}': {e}",
                entry.entry_id
            ),
        })?;
        protein.set_name(
            entry
                .protein_name
                .clone()
                .unwrap_or_else(|| entry.entry_id.clone()),
        );
        protein.set_molecule_type("protein");
        if protein.len() > 0 {
            let mut protein_qualifiers = vec![
                ("entry_id".into(), Some(entry.entry_id.clone())),
                ("accession".into(), Some(entry.accession.clone())),
                ("primary_id".into(), Some(entry.primary_id.clone())),
                (
                    "sequence_length_aa".into(),
                    Some(entry.sequence_length.to_string()),
                ),
                (
                    "synthetic_origin".into(),
                    Some("uniprot_entry_import".to_string()),
                ),
            ];
            if let Some(reviewed) = entry.reviewed {
                protein_qualifiers.push(("reviewed".into(), Some(reviewed.to_string())));
            }
            if let Some(name) = entry.protein_name.as_ref() {
                protein_qualifiers.push(("product".into(), Some(name.clone())));
                protein_qualifiers.push(("label".into(), Some(name.clone())));
            }
            if let Some(organism) = entry.organism.as_ref() {
                protein_qualifiers.push(("organism".into(), Some(organism.clone())));
            }
            for gene_name in &entry.gene_names {
                if !gene_name.trim().is_empty() {
                    protein_qualifiers.push(("gene".into(), Some(gene_name.clone())));
                }
            }
            let protein_len = protein.len();
            protein.features_mut().push(gb_io::seq::Feature {
                kind: "Protein".into(),
                location: gb_io::seq::Location::simple_range(0, protein_len as i64),
                qualifiers: protein_qualifiers,
            });
        }
        for feature in &entry.features {
            let (Some(start_aa), Some(end_aa)) =
                (feature.interval.start_aa, feature.interval.end_aa)
            else {
                continue;
            };
            if start_aa == 0 || end_aa < start_aa || end_aa > protein.len() {
                continue;
            }
            let mut qualifiers = vec![
                ("label".into(), Some(feature.key.clone())),
                (
                    "raw_location".into(),
                    Some(feature.interval.raw_location.clone()),
                ),
                ("entry_id".into(), Some(entry.entry_id.clone())),
            ];
            if let Some(note) = feature.note.as_ref() {
                qualifiers.push(("note".into(), Some(note.clone())));
            }
            for (key, value) in &feature.qualifiers {
                qualifiers.push((key.clone().into(), Some(value.clone())));
            }
            protein.features_mut().push(gb_io::seq::Feature {
                kind: feature.key.clone().into(),
                location: gb_io::seq::Location::simple_range(
                    start_aa.saturating_sub(1) as i64,
                    end_aa as i64,
                ),
                qualifiers,
            });
        }
        Self::prepare_sequence(&mut protein);
        self.state.sequences.insert(seq_id.clone(), protein);
        self.add_lineage_node(&seq_id, SequenceOrigin::Derived, Some(&result.op_id));
        result.created_seq_ids.push(seq_id.clone());
        result.messages.push(format!(
            "Imported UniProt protein sequence '{}' as '{}' ({} aa, features={}).",
            entry.entry_id,
            seq_id,
            entry.sequence_length,
            entry.features.len()
        ));
        Ok(seq_id)
    }

    pub(super) fn aa_interval_to_genomic_segments(
        aa_start: usize,
        aa_end: usize,
        aa_segments: &[UniprotAaGenomicSegment],
        is_reverse: bool,
    ) -> Vec<UniprotAaGenomicSegment> {
        let mut out = vec![];
        for segment in aa_segments {
            if segment.aa_end < aa_start || segment.aa_start > aa_end {
                continue;
            }
            let overlap_start = aa_start.max(segment.aa_start);
            let overlap_end = aa_end.min(segment.aa_end);
            if overlap_end < overlap_start {
                continue;
            }
            let seg_nt_start = segment.genomic_start_1based.min(segment.genomic_end_1based);
            let seg_nt_end = segment.genomic_start_1based.max(segment.genomic_end_1based);
            let seg_nt_len = seg_nt_end.saturating_sub(seg_nt_start).saturating_add(1);
            let seg_aa_len = segment
                .aa_end
                .saturating_sub(segment.aa_start)
                .saturating_add(1);
            if seg_nt_len == 0 || seg_aa_len == 0 {
                continue;
            }
            let nt_start_offset = overlap_start
                .saturating_sub(segment.aa_start)
                .saturating_mul(3);
            let nt_end_offset_exclusive = overlap_end
                .saturating_sub(segment.aa_start)
                .saturating_add(1)
                .saturating_mul(3);
            let nt_start_offset = nt_start_offset.min(seg_nt_len.saturating_sub(1));
            let nt_end_offset_exclusive = nt_end_offset_exclusive.min(seg_nt_len);
            if nt_end_offset_exclusive <= nt_start_offset {
                continue;
            }
            let (genomic_start_1based, genomic_end_1based) = if is_reverse {
                let end = seg_nt_end.saturating_sub(nt_start_offset);
                let start = seg_nt_end.saturating_sub(nt_end_offset_exclusive.saturating_sub(1));
                (start.min(end), start.max(end))
            } else {
                let start = seg_nt_start.saturating_add(nt_start_offset);
                let end = seg_nt_start.saturating_add(nt_end_offset_exclusive.saturating_sub(1));
                (start.min(end), start.max(end))
            };
            out.push(UniprotAaGenomicSegment {
                aa_start: overlap_start,
                aa_end: overlap_end,
                genomic_start_1based,
                genomic_end_1based,
                strand: segment.strand.clone(),
            });
        }
        out
    }

    pub fn project_uniprot_to_genome(
        &self,
        seq_id: &str,
        entry_id: &str,
        projection_id: Option<&str>,
        transcript_id_filter: Option<&str>,
    ) -> Result<UniprotGenomeProjection, EngineError> {
        let seq = self
            .state
            .sequences
            .get(seq_id)
            .ok_or_else(|| EngineError {
                code: ErrorCode::NotFound,
                message: format!("Sequence '{}' not found", seq_id),
            })?;
        let entry = self.get_uniprot_entry(entry_id)?;
        let mut warnings: Vec<String> = vec![];
        let mut transcripts_requested: Vec<String> = vec![];
        if let Some(filter) = transcript_id_filter
            .map(str::trim)
            .filter(|value| !value.is_empty())
        {
            transcripts_requested.push(filter.to_string());
        } else {
            for xref in &entry.ensembl_xrefs {
                if let Some(transcript_id) = xref
                    .transcript_id
                    .as_deref()
                    .map(str::trim)
                    .filter(|value| !value.is_empty())
                {
                    transcripts_requested.push(transcript_id.to_string());
                }
            }
        }

        let mut transcript_feature_keys: HashMap<String, usize> = HashMap::new();
        let features = seq.features();
        for (idx, feature) in features.iter().enumerate() {
            if !Self::is_mrna_feature(feature) {
                continue;
            }
            for key in Self::feature_transcript_match_keys(feature, idx) {
                transcript_feature_keys.entry(key).or_insert(idx);
            }
        }

        if transcripts_requested.is_empty() {
            warnings.push(
                "UniProt entry has no transcript xref; falling back to all mRNA features in sequence"
                    .to_string(),
            );
            for (idx, feature) in features.iter().enumerate() {
                if !Self::is_mrna_feature(feature) {
                    continue;
                }
                transcripts_requested.push(Self::feature_transcript_id(feature, idx));
            }
        }

        let mut transcript_projections: Vec<UniprotTranscriptProjection> = vec![];
        let mut seen_transcripts = BTreeSet::new();
        for transcript_probe in transcripts_requested {
            let normalized_probe = Self::normalize_transcript_probe(&transcript_probe);
            if normalized_probe.is_empty() || !seen_transcripts.insert(normalized_probe.clone()) {
                continue;
            }
            let Some(feature_index) = transcript_feature_keys.get(&normalized_probe).copied()
            else {
                warnings.push(format!(
                    "Transcript '{}' from UniProt entry '{}' was not found in sequence '{}'",
                    transcript_probe, entry.entry_id, seq_id
                ));
                continue;
            };
            let feature = &features[feature_index];
            let transcript_id = Self::feature_transcript_id(feature, feature_index);
            let is_reverse = feature_is_reverse(feature);
            let mut transcript_warnings: Vec<String> = vec![];
            let mut cds_ranges =
                Self::feature_qualifier_ranges_0based(feature, "cds_ranges_1based");
            if cds_ranges.is_empty() {
                if let Some(inferred_cds_ranges) =
                    Self::infer_cds_ranges_from_compatible_cds_features(features, feature)
                {
                    cds_ranges = inferred_cds_ranges;
                } else {
                    let mut fallback = vec![];
                    collect_location_ranges_usize(&feature.location, &mut fallback);
                    fallback.sort_unstable_by(|a, b| a.0.cmp(&b.0).then(a.1.cmp(&b.1)));
                    fallback.retain(|(start, end)| end > start);
                    cds_ranges = fallback;
                    transcript_warnings.push(format!(
                        "Transcript '{}' has no cds_ranges_1based qualifier or compatible CDS feature; fell back to feature exon spans",
                        transcript_id
                    ));
                }
            }
            let aa_reference_segments =
                Self::cds_ranges_to_reference_aa_segments(cds_ranges.clone(), is_reverse, Some(1))
                    .into_iter()
                    .map(|segment| UniprotAaGenomicSegment {
                        aa_start: segment.aa_start,
                        aa_end: segment.aa_end,
                        genomic_start_1based: segment.genomic_start_1based,
                        genomic_end_1based: segment.genomic_end_1based,
                        strand: if is_reverse {
                            "-".to_string()
                        } else {
                            "+".to_string()
                        },
                    })
                    .collect::<Vec<_>>();
            if aa_reference_segments.is_empty() {
                transcript_warnings.push(format!(
                    "Transcript '{}' has no mappable CDS amino-acid segments",
                    transcript_id
                ));
            }

            let mut feature_projections: Vec<UniprotFeatureProjection> = vec![];
            for feature in &entry.features {
                let (Some(aa_start), Some(aa_end)) =
                    (feature.interval.start_aa, feature.interval.end_aa)
                else {
                    continue;
                };
                if aa_end < aa_start {
                    continue;
                }
                let genomic_segments = Self::aa_interval_to_genomic_segments(
                    aa_start,
                    aa_end,
                    &aa_reference_segments,
                    is_reverse,
                );
                if genomic_segments.is_empty() {
                    continue;
                }
                feature_projections.push(UniprotFeatureProjection {
                    feature_key: feature.key.clone(),
                    feature_note: feature.note.clone(),
                    aa_start,
                    aa_end,
                    genomic_segments,
                });
            }
            transcript_projections.push(UniprotTranscriptProjection {
                transcript_id,
                transcript_feature_id: Some(feature_index),
                strand: if is_reverse {
                    "-".to_string()
                } else {
                    "+".to_string()
                },
                aa_segments: aa_reference_segments,
                feature_projections,
                warnings: transcript_warnings,
            });
        }

        let projection_id = projection_id
            .map(str::trim)
            .filter(|value| !value.is_empty())
            .map(|value| value.to_string())
            .unwrap_or_else(|| format!("{}@{}", entry.entry_id, seq_id));
        let projection_id = Self::normalize_uniprot_projection_id(&projection_id)?;
        Ok(UniprotGenomeProjection {
            schema: UNIPROT_GENOME_PROJECTION_SCHEMA.to_string(),
            projection_id,
            entry_id: entry.entry_id.clone(),
            seq_id: seq_id.to_string(),
            created_at_unix_ms: Self::now_unix_ms(),
            op_id: None,
            run_id: None,
            transcript_id_filter: transcript_id_filter
                .map(str::trim)
                .filter(|value| !value.is_empty())
                .map(|value| value.to_string()),
            transcript_projections,
            warnings,
        })
    }

    fn uniprot_feature_domain_name_from_parts(
        feature_key: &str,
        feature_note: Option<&str>,
    ) -> String {
        let key = feature_key.trim();
        let note = feature_note
            .map(str::trim)
            .filter(|value| !value.is_empty());
        match note {
            Some(note) if !note.eq_ignore_ascii_case(key) => format!("{key}: {note}"),
            _ => key.to_string(),
        }
    }

    fn projected_external_protein_feature_domain(
        feature_projection: &ExternalProteinOpinionFeatureProjection,
    ) -> Option<IsoformArchitectureProteinDomain> {
        let projected_start = feature_projection
            .genomic_segments
            .iter()
            .filter_map(|segment| {
                if segment.aa_start == 0 || segment.aa_end < segment.aa_start {
                    return None;
                }
                Some(segment.aa_start)
            })
            .min()?;
        let projected_end = feature_projection
            .genomic_segments
            .iter()
            .filter_map(|segment| {
                if segment.aa_start == 0 || segment.aa_end < segment.aa_start {
                    return None;
                }
                Some(segment.aa_end.max(segment.aa_start))
            })
            .max()?;
        if projected_end < projected_start {
            return None;
        }
        Some(IsoformArchitectureProteinDomain {
            name: Self::uniprot_feature_domain_name_from_parts(
                &feature_projection.feature_key,
                feature_projection.feature_note.as_deref(),
            ),
            start_aa: projected_start,
            end_aa: projected_end,
            color_hex: Some(Self::uniprot_feature_color_hex(
                &feature_projection.feature_key,
            )),
        })
    }

    fn uniprot_feature_hidden_by_default(feature_key: &str) -> bool {
        matches!(feature_key.trim().to_ascii_uppercase().as_str(), "CONFLICT")
    }

    fn external_protein_feature_hidden_by_default(
        source: gentle_protocol::ProteinExternalOpinionSource,
        feature_key: &str,
    ) -> bool {
        match source {
            gentle_protocol::ProteinExternalOpinionSource::Uniprot => {
                Self::uniprot_feature_hidden_by_default(feature_key)
            }
            gentle_protocol::ProteinExternalOpinionSource::Ensembl
            | gentle_protocol::ProteinExternalOpinionSource::Other => false,
        }
    }

    fn normalized_protein_feature_key(feature_key: &str) -> String {
        feature_key.trim().to_ascii_uppercase()
    }

    fn uniprot_feature_color_hex(feature_key: &str) -> String {
        let palette = [
            "#2563eb", "#9333ea", "#0f766e", "#b45309", "#dc2626", "#0891b2", "#7c3aed", "#4f46e5",
            "#16a34a", "#be123c",
        ];
        let mut hasher = std::collections::hash_map::DefaultHasher::new();
        use std::hash::{Hash, Hasher};
        feature_key.trim().to_ascii_uppercase().hash(&mut hasher);
        let idx = (hasher.finish() as usize) % palette.len();
        palette[idx].to_string()
    }

    fn local_ranges_0based_from_derivation(
        derivation: Option<&TranscriptProteinDerivation>,
    ) -> Vec<(usize, usize)> {
        derivation
            .into_iter()
            .flat_map(|derivation| derivation.cds_ranges_1based.iter().copied())
            .filter_map(|(start_1based, end_1based)| {
                if start_1based == 0 || end_1based < start_1based {
                    return None;
                }
                Some((start_1based.saturating_sub(1), end_1based))
            })
            .collect()
    }

    fn build_transcript_exon_segments_forward(
        exon_ranges_0based: &[(usize, usize)],
    ) -> Vec<(usize, usize, usize, usize)> {
        let mut segments = vec![];
        let mut local_cursor = 0usize;
        for (start_0based, end_0based_exclusive) in exon_ranges_0based {
            if *end_0based_exclusive <= *start_0based {
                continue;
            }
            let local_start = local_cursor;
            local_cursor = local_cursor.saturating_add(end_0based_exclusive - start_0based);
            segments.push((
                *start_0based,
                *end_0based_exclusive,
                local_start,
                local_cursor,
            ));
        }
        segments
    }

    fn build_transcript_introns_from_ranges_1based(
        ranges_1based: &[(usize, usize)],
    ) -> Vec<(usize, usize)> {
        let mut ordered = ranges_1based.to_vec();
        ordered.sort_unstable_by(|a, b| a.0.cmp(&b.0).then(a.1.cmp(&b.1)));
        let mut introns = vec![];
        for pair in ordered.windows(2) {
            if pair[1].0 > pair[0].1.saturating_add(1) {
                introns.push((pair[0].1.saturating_add(1), pair[1].0.saturating_sub(1)));
            }
        }
        introns
    }

    fn build_derived_protein_expert_transcript(
        source_sequence_upper: &[u8],
        source_feature: &gb_io::seq::Feature,
        source_features: &[gb_io::seq::Feature],
        source_feature_id: usize,
        source_seq_id: &str,
    ) -> Result<DerivedProteinExpertTranscript, EngineError> {
        let exon_ranges_0based = Self::feature_location_ranges_0based(source_feature);
        let exon_segments_forward =
            Self::build_transcript_exon_segments_forward(&exon_ranges_0based);
        let transcript_exons_1based = Self::feature_location_ranges_1based(source_feature);
        let (
            derived_transcript,
            transcript_id,
            transcript_label,
            is_reverse,
            _exon_count,
            annotated_derivation,
        ) = Self::derive_transcript_sequence_from_feature(
            source_sequence_upper,
            source_feature,
            source_features,
            source_feature_id,
            source_seq_id,
        )?;
        let derivation = match annotated_derivation {
            Some(derivation) => Some(derivation),
            None => Self::infer_transcript_protein_derivation_without_annotation(
                &derived_transcript.get_forward_string(),
                source_feature,
                source_feature_id,
                source_seq_id,
                source_features,
                &transcript_id,
                &transcript_label,
            )?,
        };
        let local_cds_ranges_0based =
            Self::local_ranges_0based_from_derivation(derivation.as_ref());
        let genomic_cds_ranges_0based = if !local_cds_ranges_0based.is_empty() {
            Self::map_transcript_local_ranges_to_source_ranges_0based(
                &local_cds_ranges_0based,
                &exon_segments_forward,
                is_reverse,
                derived_transcript.len(),
            )
        } else {
            vec![]
        };
        let genomic_cds_ranges_1based = genomic_cds_ranges_0based
            .iter()
            .map(|(start_0based, end_0based_exclusive)| {
                (start_0based.saturating_add(1), *end_0based_exclusive)
            })
            .collect::<Vec<_>>();
        let cds_to_protein_segments = if genomic_cds_ranges_0based.is_empty() {
            vec![]
        } else {
            Self::cds_ranges_to_reference_aa_segments(
                genomic_cds_ranges_0based,
                is_reverse,
                Some(1),
            )
        };
        let intron_ranges_1based =
            Self::build_transcript_introns_from_ranges_1based(&genomic_cds_ranges_1based);
        Ok(DerivedProteinExpertTranscript {
            transcript_id,
            transcript_label,
            transcript_feature_id: source_feature_id,
            is_reverse,
            transcript_exons_1based,
            genomic_cds_ranges_1based,
            intron_ranges_1based,
            cds_to_protein_segments,
            derivation,
        })
    }

    fn ranges_overlap_1based(left: (usize, usize), right: (usize, usize)) -> bool {
        left.0 <= right.1 && right.0 <= left.1
    }

    fn internal_nonoverlapping_ranges_1based(
        ranges_1based: &[(usize, usize)],
        other_ranges_1based: &[(usize, usize)],
        is_reverse: bool,
    ) -> Vec<(usize, usize)> {
        let ordered = Self::order_ranges_for_transcript(
            ranges_1based
                .iter()
                .map(|(start_1based, end_1based)| (start_1based.saturating_sub(1), *end_1based))
                .collect(),
            is_reverse,
        )
        .into_iter()
        .map(|(start_0based, end_0based_exclusive)| {
            (start_0based.saturating_add(1), end_0based_exclusive)
        })
        .collect::<Vec<_>>();
        ordered
            .iter()
            .enumerate()
            .filter_map(|(idx, range)| {
                let internal = idx > 0 && idx + 1 < ordered.len();
                let overlaps_other = other_ranges_1based
                    .iter()
                    .copied()
                    .any(|other| Self::ranges_overlap_1based(*range, other));
                (internal && !overlaps_other).then_some(*range)
            })
            .collect()
    }

    fn map_genomic_range_to_local_aa_ranges(
        genomic_start_1based: usize,
        genomic_end_1based: usize,
        transcript_segments: &[IsoformArchitectureCdsAaSegment],
        is_reverse: bool,
    ) -> Vec<(usize, usize)> {
        let interval_start = genomic_start_1based.min(genomic_end_1based);
        let interval_end = genomic_start_1based.max(genomic_end_1based);
        let mut out = vec![];
        for segment in transcript_segments {
            let segment_start = segment.genomic_start_1based.min(segment.genomic_end_1based);
            let segment_end = segment.genomic_start_1based.max(segment.genomic_end_1based);
            let overlap_start = interval_start.max(segment_start);
            let overlap_end = interval_end.min(segment_end);
            if overlap_end < overlap_start {
                continue;
            }
            let local_start = if is_reverse {
                segment
                    .aa_start
                    .saturating_add(segment_end.saturating_sub(overlap_end) / 3)
            } else {
                segment
                    .aa_start
                    .saturating_add(overlap_start.saturating_sub(segment_start) / 3)
            };
            let local_end = if is_reverse {
                segment
                    .aa_start
                    .saturating_add(segment_end.saturating_sub(overlap_start).saturating_add(1) / 3)
            } else {
                segment
                    .aa_start
                    .saturating_add(overlap_end.saturating_sub(segment_start).saturating_add(1) / 3)
            }
            .saturating_sub(1);
            if local_end >= local_start {
                out.push((local_start, local_end));
            }
        }
        out
    }

    fn projected_external_protein_feature_domain_on_transcript(
        feature_projection: &ExternalProteinOpinionFeatureProjection,
        transcript_segments: &[IsoformArchitectureCdsAaSegment],
        is_reverse: bool,
    ) -> Option<IsoformArchitectureProteinDomain> {
        if transcript_segments.is_empty() {
            return Self::projected_external_protein_feature_domain(feature_projection);
        }
        let mut mapped_local_ranges = vec![];
        for genomic_segment in &feature_projection.genomic_segments {
            mapped_local_ranges.extend(Self::map_genomic_range_to_local_aa_ranges(
                genomic_segment.genomic_start_1based,
                genomic_segment.genomic_end_1based,
                transcript_segments,
                is_reverse,
            ));
        }
        let projected_start = mapped_local_ranges.iter().map(|(start, _)| *start).min()?;
        let projected_end = mapped_local_ranges.iter().map(|(_, end)| *end).max()?;
        (projected_end >= projected_start).then(|| IsoformArchitectureProteinDomain {
            name: Self::uniprot_feature_domain_name_from_parts(
                &feature_projection.feature_key,
                feature_projection.feature_note.as_deref(),
            ),
            start_aa: projected_start,
            end_aa: projected_end,
            color_hex: Some(Self::uniprot_feature_color_hex(
                &feature_projection.feature_key,
            )),
        })
    }

    fn external_protein_transcript_coding_ranges_1based(
        transcript: &ExternalProteinOpinionTranscript,
    ) -> Vec<(usize, usize)> {
        let mut ranges = transcript
            .aa_segments
            .iter()
            .map(|segment| {
                (
                    segment.genomic_start_1based.min(segment.genomic_end_1based),
                    segment.genomic_start_1based.max(segment.genomic_end_1based),
                )
            })
            .filter(|(start, end)| *start > 0 && *end >= *start)
            .collect::<Vec<_>>();
        ranges.sort_unstable_by(|a, b| a.0.cmp(&b.0).then(a.1.cmp(&b.1)));
        ranges.dedup();
        ranges
    }

    fn build_external_protein_opinion_summary(
        source: &ExternalProteinOpinionViewSource,
        transcript: &ExternalProteinOpinionTranscript,
    ) -> gentle_protocol::TranscriptProteinExternalOpinion {
        let reference_start_aa = transcript
            .aa_segments
            .iter()
            .map(|segment| segment.aa_start)
            .min();
        let reference_end_aa = transcript
            .aa_segments
            .iter()
            .map(|segment| segment.aa_end)
            .max();
        gentle_protocol::TranscriptProteinExternalOpinion {
            source: source.source,
            source_id: source.source_id.clone(),
            source_label: source.source_label.clone(),
            expected_length_aa: source.expected_length_aa,
            reference_start_aa,
            reference_end_aa,
            genomic_coding_ranges_1based: Self::external_protein_transcript_coding_ranges_1based(
                transcript,
            ),
        }
    }

    fn build_transcript_protein_comparison(
        transcript_id: &str,
        transcript_label: &str,
        is_reverse: bool,
        derivation: Option<&TranscriptProteinDerivation>,
        derived_coding_ranges_1based: &[(usize, usize)],
        external_opinion: Option<gentle_protocol::TranscriptProteinExternalOpinion>,
    ) -> gentle_protocol::TranscriptProteinComparison {
        let external_ranges_1based = external_opinion
            .as_ref()
            .map(|opinion| opinion.genomic_coding_ranges_1based.clone())
            .unwrap_or_default();
        let derived_only_internal = Self::internal_nonoverlapping_ranges_1based(
            derived_coding_ranges_1based,
            &external_ranges_1based,
            is_reverse,
        );
        let external_only_internal = Self::internal_nonoverlapping_ranges_1based(
            &external_ranges_1based,
            derived_coding_ranges_1based,
            is_reverse,
        );
        let mut mismatch_reasons = vec![];
        if !derived_only_internal.is_empty() {
            mismatch_reasons.push(format!(
                "Derived transcript translation contains internal coding exon(s) absent from the external protein opinion: {}",
                derived_only_internal
                    .iter()
                    .map(|(start, end)| format!("{start}..{end}"))
                    .collect::<Vec<_>>()
                    .join(", ")
            ));
        }
        if !external_only_internal.is_empty() {
            mismatch_reasons.push(format!(
                "External protein opinion contains internal coding exon(s) absent from the transcript-native translation: {}",
                external_only_internal
                    .iter()
                    .map(|(start, end)| format!("{start}..{end}"))
                    .collect::<Vec<_>>()
                    .join(", ")
            ));
        }
        if let (Some(derivation), Some(external_opinion)) = (derivation, external_opinion.as_ref())
        {
            if let Some(external_length_aa) = external_opinion.expected_length_aa
                && derivation.protein_length_aa > 0
                && derivation.protein_length_aa != external_length_aa
            {
                mismatch_reasons.push(format!(
                    "Transcript-native product length is {} aa, while the external protein opinion expects {} aa.",
                    derivation.protein_length_aa, external_length_aa
                ));
            }
        }

        let has_derived = derivation
            .map(|derivation| {
                derivation.protein_length_aa > 0 || !derivation.cds_ranges_1based.is_empty()
            })
            .unwrap_or(false);
        let has_external = external_opinion.is_some();
        let status = match (has_derived, has_external) {
            (true, false) => gentle_protocol::TranscriptProteinComparisonStatus::DerivedOnly,
            (true, true) if mismatch_reasons.is_empty() => {
                gentle_protocol::TranscriptProteinComparisonStatus::ConsistentWithExternalOpinion
            }
            (true, true) => {
                gentle_protocol::TranscriptProteinComparisonStatus::LowConfidenceExternalOpinion
            }
            (false, true) => {
                gentle_protocol::TranscriptProteinComparisonStatus::ExternalOpinionOnly
            }
            (false, false) => gentle_protocol::TranscriptProteinComparisonStatus::NoTranscriptCds,
        };

        gentle_protocol::TranscriptProteinComparison {
            transcript_id: transcript_id.to_string(),
            transcript_label: transcript_label.to_string(),
            status,
            derived: derivation.cloned(),
            external_opinion,
            mismatch_reasons,
            derived_only_exon_ranges_1based: derived_only_internal,
            external_only_exon_ranges_1based: external_only_internal,
        }
    }

    fn build_uniprot_external_protein_opinion_source(
        projection: &UniprotGenomeProjection,
        entry: &UniprotEntry,
    ) -> ExternalProteinOpinionViewSource {
        let source_label = format!("UniProt {} ({})", entry.entry_id, entry.accession);
        let gene_symbol_hint = entry
            .gene_names
            .iter()
            .find(|value| !value.trim().is_empty())
            .cloned()
            .or_else(|| {
                entry
                    .protein_name
                    .clone()
                    .filter(|value| !value.trim().is_empty())
            });
        let feature_inventory_keys = entry
            .features
            .iter()
            .filter_map(|feature| {
                let (Some(start_aa), Some(end_aa)) =
                    (feature.interval.start_aa, feature.interval.end_aa)
                else {
                    return None;
                };
                if start_aa == 0 || end_aa < start_aa || feature.key.eq_ignore_ascii_case("CHAIN") {
                    return None;
                }
                Some(Self::normalized_protein_feature_key(&feature.key))
            })
            .collect::<Vec<_>>();
        let transcripts = projection
            .transcript_projections
            .iter()
            .map(|transcript| ExternalProteinOpinionTranscript {
                transcript_id: transcript.transcript_id.clone(),
                transcript_feature_id: transcript.transcript_feature_id,
                strand: transcript.strand.clone(),
                aa_segments: transcript
                    .aa_segments
                    .iter()
                    .map(|segment| ExternalProteinOpinionAaSegment {
                        aa_start: segment.aa_start,
                        aa_end: segment.aa_end,
                        genomic_start_1based: segment.genomic_start_1based,
                        genomic_end_1based: segment.genomic_end_1based,
                    })
                    .collect(),
                feature_projections: transcript
                    .feature_projections
                    .iter()
                    .map(|feature| ExternalProteinOpinionFeatureProjection {
                        feature_key: feature.feature_key.clone(),
                        feature_note: feature.feature_note.clone(),
                        genomic_segments: feature
                            .genomic_segments
                            .iter()
                            .map(|segment| ExternalProteinOpinionAaSegment {
                                aa_start: segment.aa_start,
                                aa_end: segment.aa_end,
                                genomic_start_1based: segment.genomic_start_1based,
                                genomic_end_1based: segment.genomic_end_1based,
                            })
                            .collect(),
                    })
                    .collect(),
            })
            .collect::<Vec<_>>();

        ExternalProteinOpinionViewSource {
            source: gentle_protocol::ProteinExternalOpinionSource::Uniprot,
            source_id: projection.projection_id.clone(),
            source_label: source_label.clone(),
            panel_source_label: format!(
                "Transcript-native proteins with optional UniProt opinion {} ({})",
                entry.entry_id, entry.accession
            ),
            transcript_id_filter: projection.transcript_id_filter.clone(),
            gene_symbol_hint,
            expected_length_aa: (entry.sequence_length > 0).then_some(entry.sequence_length),
            feature_inventory_keys,
            warnings: projection.warnings.clone(),
            transcripts,
        }
    }

    fn build_ensembl_external_protein_opinion_source(
        &self,
        seq_id: &str,
        transcript_id_filter: Option<&str>,
        entry: &EnsemblProteinEntry,
    ) -> Result<ExternalProteinOpinionViewSource, EngineError> {
        let dna = self
            .state
            .sequences
            .get(seq_id)
            .ok_or_else(|| EngineError {
                code: ErrorCode::NotFound,
                message: format!("Sequence '{seq_id}' not found"),
            })?;
        let features = dna.features();
        let source_sequence_upper = dna.get_forward_string().to_ascii_uppercase().into_bytes();
        let normalized_entry_transcript = Self::normalize_transcript_probe(&entry.transcript_id);
        let normalized_filter = transcript_id_filter
            .map(Self::normalize_transcript_probe)
            .filter(|value| !value.is_empty());
        if let Some(filter) = normalized_filter.as_deref()
            && filter != normalized_entry_transcript
        {
            return Err(EngineError {
                code: ErrorCode::NotFound,
                message: format!(
                    "Ensembl protein entry '{}' belongs to transcript '{}' rather than '{}'",
                    entry.entry_id, entry.transcript_id, filter
                ),
            });
        }

        let feature_inventory_keys = entry
            .features
            .iter()
            .filter_map(|feature| {
                (feature.start_aa.is_some() && feature.end_aa.is_some())
                    .then_some(Self::normalized_protein_feature_key(&feature.feature_key))
            })
            .collect::<Vec<_>>();
        let mut warnings = vec![];
        let mut transcripts = vec![];

        for (feature_id, feature) in features.iter().enumerate() {
            if !Self::is_transcript_feature_for_derivation(feature) {
                continue;
            }
            let match_keys = Self::feature_transcript_match_keys(feature, feature_id);
            if !match_keys
                .iter()
                .any(|key| key == &normalized_entry_transcript)
            {
                continue;
            }
            match Self::build_derived_protein_expert_transcript(
                &source_sequence_upper,
                feature,
                features,
                feature_id,
                seq_id,
            ) {
                Ok(derived) => {
                    let aa_segments = derived
                        .cds_to_protein_segments
                        .iter()
                        .map(|segment| ExternalProteinOpinionAaSegment {
                            aa_start: segment.aa_start,
                            aa_end: segment.aa_end,
                            genomic_start_1based: segment.genomic_start_1based,
                            genomic_end_1based: segment.genomic_end_1based,
                        })
                        .collect::<Vec<_>>();
                    let reference_segments = derived
                        .cds_to_protein_segments
                        .iter()
                        .map(|segment| UniprotAaGenomicSegment {
                            aa_start: segment.aa_start,
                            aa_end: segment.aa_end,
                            genomic_start_1based: segment.genomic_start_1based,
                            genomic_end_1based: segment.genomic_end_1based,
                            strand: if derived.is_reverse {
                                "-".to_string()
                            } else {
                                "+".to_string()
                            },
                        })
                        .collect::<Vec<_>>();
                    let feature_projections = entry
                        .features
                        .iter()
                        .filter_map(|ensembl_feature| {
                            let (Some(aa_start), Some(aa_end)) =
                                (ensembl_feature.start_aa, ensembl_feature.end_aa)
                            else {
                                return None;
                            };
                            if aa_start == 0 || aa_end < aa_start {
                                return None;
                            }
                            let genomic_segments = Self::aa_interval_to_genomic_segments(
                                aa_start,
                                aa_end,
                                &reference_segments,
                                derived.is_reverse,
                            )
                            .into_iter()
                            .map(|segment| ExternalProteinOpinionAaSegment {
                                aa_start: segment.aa_start,
                                aa_end: segment.aa_end,
                                genomic_start_1based: segment.genomic_start_1based,
                                genomic_end_1based: segment.genomic_end_1based,
                            })
                            .collect::<Vec<_>>();
                            if genomic_segments.is_empty() {
                                return None;
                            }
                            Some(ExternalProteinOpinionFeatureProjection {
                                feature_key: ensembl_feature.feature_key.clone(),
                                feature_note: ensembl_feature
                                    .description
                                    .clone()
                                    .or_else(|| Some(ensembl_feature.feature_type.clone())),
                                genomic_segments,
                            })
                        })
                        .collect::<Vec<_>>();
                    transcripts.push(ExternalProteinOpinionTranscript {
                        transcript_id: derived.transcript_id.clone(),
                        transcript_feature_id: Some(derived.transcript_feature_id),
                        strand: if derived.is_reverse {
                            "-".to_string()
                        } else {
                            "+".to_string()
                        },
                        aa_segments,
                        feature_projections,
                    });
                }
                Err(err) => warnings.push(format!(
                    "Transcript '{}' could not be prepared for Ensembl protein comparison: {}",
                    entry.transcript_id, err.message
                )),
            }
        }

        if transcripts.is_empty() {
            return Err(EngineError {
                code: ErrorCode::NotFound,
                message: format!(
                    "Ensembl protein entry '{}' transcript '{}' was not found as a derivable transcript in sequence '{}'",
                    entry.entry_id, entry.transcript_id, seq_id
                ),
            });
        }

        Ok(ExternalProteinOpinionViewSource {
            source: gentle_protocol::ProteinExternalOpinionSource::Ensembl,
            source_id: entry.entry_id.clone(),
            source_label: format!("Ensembl {} ({})", entry.protein_id, entry.transcript_id),
            panel_source_label: format!(
                "Transcript-native proteins with optional Ensembl opinion {} ({})",
                entry.protein_id, entry.transcript_id
            ),
            transcript_id_filter: transcript_id_filter
                .map(str::trim)
                .filter(|value| !value.is_empty())
                .map(|value| value.to_string()),
            gene_symbol_hint: entry
                .gene_symbol
                .clone()
                .or_else(|| entry.transcript_display_name.clone()),
            expected_length_aa: (entry.sequence_length > 0).then_some(entry.sequence_length),
            feature_inventory_keys,
            warnings,
            transcripts,
        })
    }

    fn materialize_external_protein_opinion_provider(
        provider: Option<&ExternalProteinOpinionProvider<'_>>,
    ) -> Option<ExternalProteinOpinionViewSource> {
        match provider {
            Some(ExternalProteinOpinionProvider::Materialized(source)) => Some(source.clone()),
            Some(ExternalProteinOpinionProvider::UniprotProjection { projection, entry }) => Some(
                Self::build_uniprot_external_protein_opinion_source(projection, entry),
            ),
            None => None,
        }
    }

    fn build_external_protein_expert_view_from_source(
        &self,
        seq_id: &str,
        transcript_id_filter: Option<&str>,
        external_source: ExternalProteinOpinionViewSource,
        protein_feature_filter: &ProteinFeatureFilter,
    ) -> Result<IsoformArchitectureExpertView, EngineError> {
        let provider = ExternalProteinOpinionProvider::Materialized(external_source);
        self.build_optional_external_protein_expert_view(
            seq_id,
            transcript_id_filter,
            Some(&provider),
            protein_feature_filter,
        )
    }

    #[cfg(test)]
    pub(crate) fn build_external_protein_expert_view_for_test(
        &self,
        seq_id: &str,
        transcript_id_filter: Option<&str>,
        external_source: ExternalProteinOpinionViewSource,
        protein_feature_filter: &ProteinFeatureFilter,
    ) -> Result<IsoformArchitectureExpertView, EngineError> {
        self.build_external_protein_expert_view_from_source(
            seq_id,
            transcript_id_filter,
            external_source,
            protein_feature_filter,
        )
    }

    fn feature_location_ranges_1based(feature: &gb_io::seq::Feature) -> Vec<(usize, usize)> {
        let mut ranges = vec![];
        let mut raw_ranges = vec![];
        collect_location_ranges_usize(&feature.location, &mut raw_ranges);
        ranges.extend(
            raw_ranges
                .into_iter()
                .filter_map(|(start_0based, end_0based_exclusive)| {
                    if end_0based_exclusive <= start_0based {
                        return None;
                    }
                    Some((start_0based.saturating_add(1), end_0based_exclusive))
                }),
        );
        if ranges.is_empty()
            && let Ok((from, to)) = feature.location.find_bounds()
            && from >= 0
            && to >= 0
        {
            ranges.push((from as usize + 1, to as usize));
        }
        ranges.sort_unstable_by(|a, b| a.0.cmp(&b.0).then(a.1.cmp(&b.1)));
        ranges.retain(|(start, end)| *end >= *start);
        ranges
    }

    fn feature_location_ranges_0based(feature: &gb_io::seq::Feature) -> Vec<(usize, usize)> {
        let mut ranges = vec![];
        collect_location_ranges_usize(&feature.location, &mut ranges);
        ranges.sort_unstable_by(|a, b| a.0.cmp(&b.0).then(a.1.cmp(&b.1)));
        ranges.retain(|(start, end)| *end > *start);
        ranges
    }

    fn is_cds_feature(feature: &gb_io::seq::Feature) -> bool {
        feature.kind.eq_ignore_ascii_case("CDS")
    }

    fn ranges_fully_covered_0based(
        candidate_ranges: &[(usize, usize)],
        container_ranges: &[(usize, usize)],
    ) -> bool {
        candidate_ranges
            .iter()
            .all(|(candidate_start, candidate_end)| {
                container_ranges
                    .iter()
                    .any(|(container_start, container_end)| {
                        candidate_start >= container_start && candidate_end <= container_end
                    })
            })
    }

    fn infer_cds_ranges_from_compatible_cds_features(
        features: &[gb_io::seq::Feature],
        transcript_feature: &gb_io::seq::Feature,
    ) -> Option<Vec<(usize, usize)>> {
        let transcript_ranges = Self::feature_location_ranges_0based(transcript_feature);
        if transcript_ranges.is_empty() {
            return None;
        }
        let transcript_is_reverse = feature_is_reverse(transcript_feature);
        let mut candidates: BTreeSet<Vec<(usize, usize)>> = BTreeSet::new();
        for feature in features {
            if !Self::is_cds_feature(feature)
                || feature_is_reverse(feature) != transcript_is_reverse
            {
                continue;
            }
            let candidate_ranges = Self::feature_location_ranges_0based(feature);
            if candidate_ranges.is_empty() {
                continue;
            }
            if Self::ranges_fully_covered_0based(&candidate_ranges, &transcript_ranges) {
                candidates.insert(candidate_ranges);
            }
        }
        if candidates.is_empty() {
            return None;
        }
        let mut ranked = candidates.into_iter().collect::<Vec<_>>();
        ranked.sort_by(|left, right| {
            let left_bp = left
                .iter()
                .map(|(start, end)| end.saturating_sub(*start))
                .sum::<usize>();
            let right_bp = right
                .iter()
                .map(|(start, end)| end.saturating_sub(*start))
                .sum::<usize>();
            right_bp
                .cmp(&left_bp)
                .then(left.len().cmp(&right.len()))
                .then(left.cmp(right))
        });
        let best = ranked.first()?.clone();
        let best_bp = best
            .iter()
            .map(|(start, end)| end.saturating_sub(*start))
            .sum::<usize>();
        let competing_best = ranked.iter().skip(1).any(|ranges| {
            let bp = ranges
                .iter()
                .map(|(start, end)| end.saturating_sub(*start))
                .sum::<usize>();
            bp == best_bp && ranges != &best
        });
        if competing_best {
            return None;
        }
        Some(best)
    }

    pub(super) fn build_transcript_protein_expert_view(
        &self,
        seq_id: &str,
        transcript_id_filter: Option<&str>,
        protein_feature_filter: &ProteinFeatureFilter,
    ) -> Result<IsoformArchitectureExpertView, EngineError> {
        self.build_optional_external_protein_expert_view(
            seq_id,
            transcript_id_filter,
            None,
            protein_feature_filter,
        )
    }

    pub(super) fn build_ensembl_entry_protein_expert_view(
        &self,
        seq_id: &str,
        transcript_id_filter: Option<&str>,
        entry_id: &str,
        protein_feature_filter: &ProteinFeatureFilter,
    ) -> Result<IsoformArchitectureExpertView, EngineError> {
        let entry = self.get_ensembl_protein_entry(entry_id)?;
        let external_source = self.build_ensembl_external_protein_opinion_source(
            seq_id,
            transcript_id_filter,
            &entry,
        )?;
        self.build_external_protein_expert_view_from_source(
            seq_id,
            transcript_id_filter,
            external_source,
            protein_feature_filter,
        )
    }

    fn build_optional_external_protein_expert_view(
        &self,
        seq_id: &str,
        transcript_id_filter: Option<&str>,
        external_provider: Option<&ExternalProteinOpinionProvider<'_>>,
        protein_feature_filter: &ProteinFeatureFilter,
    ) -> Result<IsoformArchitectureExpertView, EngineError> {
        let materialized_external_source =
            Self::materialize_external_protein_opinion_provider(external_provider);
        let external_source = materialized_external_source.as_ref();
        let dna = self
            .state
            .sequences
            .get(seq_id)
            .ok_or_else(|| EngineError {
                code: ErrorCode::NotFound,
                message: format!("Sequence '{seq_id}' not found"),
            })?;
        let features = dna.features();
        let source_sequence_upper = dna.get_forward_string().to_ascii_uppercase().into_bytes();
        let normalized_filter = transcript_id_filter
            .map(Self::normalize_transcript_probe)
            .filter(|value| !value.is_empty());
        let include_feature_keys = protein_feature_filter
            .include_feature_keys
            .iter()
            .map(|value| Self::normalized_protein_feature_key(value))
            .filter(|value| !value.is_empty())
            .collect::<BTreeSet<_>>();
        let exclude_feature_keys = protein_feature_filter
            .exclude_feature_keys
            .iter()
            .map(|value| Self::normalized_protein_feature_key(value))
            .filter(|value| !value.is_empty())
            .collect::<BTreeSet<_>>();

        let mut warnings = external_source
            .map(|external_source| external_source.warnings.clone())
            .unwrap_or_default();
        if let Some(external_source) = external_source {
            let mut hidden_feature_counts: BTreeMap<String, usize> = BTreeMap::new();
            let mut explicit_filter_hidden_counts: BTreeMap<String, usize> = BTreeMap::new();
            for normalized_key in &external_source.feature_inventory_keys {
                if Self::external_protein_feature_hidden_by_default(
                    external_source.source,
                    normalized_key,
                ) && !include_feature_keys.contains(normalized_key.as_str())
                {
                    *hidden_feature_counts
                        .entry(normalized_key.clone())
                        .or_insert(0) += 1;
                    continue;
                }
                if (!include_feature_keys.is_empty()
                    && !include_feature_keys.contains(normalized_key.as_str()))
                    || exclude_feature_keys.contains(normalized_key.as_str())
                {
                    *explicit_filter_hidden_counts
                        .entry(normalized_key.clone())
                        .or_insert(0) += 1;
                }
            }
            if !hidden_feature_counts.is_empty() {
                let summary = hidden_feature_counts
                    .iter()
                    .map(|(feature_key, count)| format!("{feature_key} x{count}"))
                    .collect::<Vec<_>>()
                    .join(", ");
                warnings.push(format!(
                    "Default external protein-feature filter hid {summary}; hidden feature classes can be added back later through explicit filter controls."
                ));
            }
            if !explicit_filter_hidden_counts.is_empty() {
                let summary = explicit_filter_hidden_counts
                    .iter()
                    .map(|(feature_key, count)| format!("{feature_key} x{count}"))
                    .collect::<Vec<_>>()
                    .join(", ");
                warnings.push(format!(
                    "Protein-feature filter hid {summary} for this protein expert view."
                ));
            }
        }

        let mut derived_transcripts: Vec<DerivedProteinExpertTranscript> = vec![];
        for (feature_id, feature) in features.iter().enumerate() {
            if !Self::is_transcript_feature_for_derivation(feature) {
                continue;
            }
            let match_keys = Self::feature_transcript_match_keys(feature, feature_id);
            if let Some(filter) = normalized_filter.as_ref()
                && !match_keys.iter().any(|key| key == filter)
            {
                continue;
            }
            match Self::build_derived_protein_expert_transcript(
                &source_sequence_upper,
                feature,
                features,
                feature_id,
                seq_id,
            ) {
                Ok(derived) => derived_transcripts.push(derived),
                Err(err) => warnings.push(format!(
                    "Transcript feature n-{} could not be prepared for transcript-native protein comparison: {}",
                    feature_id + 1,
                    err.message
                )),
            }
        }

        let external_transcripts = external_source
            .map(|external_source| external_source.transcripts.clone())
            .unwrap_or_default();
        let mut external_key_index: HashMap<String, usize> = HashMap::new();
        for (idx, transcript) in external_transcripts.iter().enumerate() {
            let key = Self::normalize_transcript_probe(&transcript.transcript_id);
            if !key.is_empty() {
                external_key_index.entry(key).or_insert(idx);
            }
        }

        let mut transcript_lanes = vec![];
        let mut protein_lanes = vec![];
        let mut region_ranges_1based: Vec<(usize, usize)> = vec![];
        let mut matched_external_indices = BTreeSet::new();

        for derived in &derived_transcripts {
            let transcript_key = Self::normalize_transcript_probe(&derived.transcript_id);
            let external_transcript = external_key_index
                .get(&transcript_key)
                .and_then(|idx| external_transcripts.get(*idx));
            if let Some(external_idx) = external_key_index.get(&transcript_key) {
                matched_external_indices.insert(*external_idx);
            }

            let external_opinion = if let (Some(external_source), Some(external_transcript)) =
                (external_source, external_transcript)
            {
                Some(Self::build_external_protein_opinion_summary(
                    external_source,
                    external_transcript,
                ))
            } else {
                None
            };
            if let Some(derivation) = derived.derivation.as_ref() {
                for warning in &derivation.warnings {
                    warnings.push(format!(
                        "Transcript '{}': {}",
                        derived.transcript_id, warning
                    ));
                }
            }
            let comparison = Self::build_transcript_protein_comparison(
                &derived.transcript_id,
                &derived.transcript_label,
                derived.is_reverse,
                derived.derivation.as_ref(),
                &derived.genomic_cds_ranges_1based,
                external_opinion,
            );

            let transcript_note = match comparison.status {
                gentle_protocol::TranscriptProteinComparisonStatus::LowConfidenceExternalOpinion => {
                    Some(format!(
                        "low_confidence_external_opinion: {}",
                        comparison.mismatch_reasons.join(" | ")
                    ))
                }
                gentle_protocol::TranscriptProteinComparisonStatus::ExternalOpinionOnly => Some(
                    "Only an external protein opinion was available; transcript-native CDS translation could not be derived."
                        .to_string(),
                ),
                gentle_protocol::TranscriptProteinComparisonStatus::NoTranscriptCds => Some(
                    "Transcript-native CDS translation could not be derived for this transcript."
                        .to_string(),
                ),
                _ => None,
            };

            let domains = external_transcript
                .into_iter()
                .flat_map(|external_transcript| external_transcript.feature_projections.iter())
                .filter_map(|feature_projection| {
                    let normalized_key =
                        Self::normalized_protein_feature_key(&feature_projection.feature_key);
                    if external_source
                        .map(|external_source| {
                            Self::external_protein_feature_hidden_by_default(
                                external_source.source,
                                &normalized_key,
                            )
                        })
                        .unwrap_or(false)
                        && !include_feature_keys.contains(&normalized_key)
                    {
                        return None;
                    }
                    if (!include_feature_keys.is_empty()
                        && !include_feature_keys.contains(&normalized_key))
                        || exclude_feature_keys.contains(&normalized_key)
                    {
                        return None;
                    }
                    Self::projected_external_protein_feature_domain_on_transcript(
                        feature_projection,
                        &derived.cds_to_protein_segments,
                        derived.is_reverse,
                    )
                })
                .collect::<Vec<_>>();

            region_ranges_1based.extend(derived.transcript_exons_1based.iter().copied());
            region_ranges_1based.extend(derived.genomic_cds_ranges_1based.iter().copied());

            transcript_lanes.push(IsoformArchitectureTranscriptLane {
                isoform_id: derived.transcript_id.clone(),
                label: derived.transcript_label.clone(),
                transcript_id: Some(derived.transcript_id.clone()),
                transcript_feature_id: Some(derived.transcript_feature_id),
                strand: if derived.is_reverse {
                    "-".to_string()
                } else {
                    "+".to_string()
                },
                transcript_exons: Self::range_vec_1based_to_splicing(
                    derived.transcript_exons_1based.clone(),
                ),
                exons: Self::range_vec_1based_to_splicing(
                    derived.genomic_cds_ranges_1based.clone(),
                ),
                introns: Self::range_vec_1based_to_splicing(derived.intron_ranges_1based.clone()),
                mapped: derived
                    .derivation
                    .as_ref()
                    .map(|derivation| derivation.protein_length_aa > 0)
                    .unwrap_or(false),
                transactivation_class: None,
                cds_to_protein_segments: derived.cds_to_protein_segments.clone(),
                note: transcript_note,
            });

            protein_lanes.push(IsoformArchitectureProteinLane {
                isoform_id: derived.transcript_id.clone(),
                label: derived.transcript_label.clone(),
                transcript_id: Some(derived.transcript_id.clone()),
                expected_length_aa: derived
                    .derivation
                    .as_ref()
                    .map(|derivation| derivation.protein_length_aa)
                    .filter(|value| *value > 0)
                    .or_else(|| {
                        comparison
                            .external_opinion
                            .as_ref()
                            .and_then(|opinion| opinion.expected_length_aa)
                    }),
                reference_start_aa: comparison
                    .external_opinion
                    .as_ref()
                    .and_then(|opinion| opinion.reference_start_aa),
                reference_end_aa: comparison
                    .external_opinion
                    .as_ref()
                    .and_then(|opinion| opinion.reference_end_aa),
                domains,
                transactivation_class: None,
                comparison: Some(comparison),
            });
        }

        for (external_idx, transcript) in external_transcripts.iter().enumerate() {
            if matched_external_indices.contains(&external_idx) {
                continue;
            }
            let transcript_key = Self::normalize_transcript_probe(&transcript.transcript_id);
            if let Some(filter) = normalized_filter.as_ref()
                && &transcript_key != filter
            {
                continue;
            }
            let transcript_feature = transcript.transcript_feature_id.and_then(|feature_id| {
                features
                    .get(feature_id)
                    .map(|feature| (feature_id, feature))
            });
            let transcript_exons = transcript_feature
                .as_ref()
                .map(|(_, feature)| Self::feature_location_ranges_1based(feature))
                .unwrap_or_default();
            let transcript_label = transcript_feature
                .as_ref()
                .and_then(|(_, feature)| {
                    Self::first_nonempty_feature_qualifier(
                        feature,
                        &[
                            "label",
                            "name",
                            "standard_name",
                            "product",
                            "transcript_id",
                            "gene",
                        ],
                    )
                })
                .unwrap_or_else(|| transcript.transcript_id.clone());
            let is_reverse = transcript_feature
                .as_ref()
                .map(|(_, feature)| feature_is_reverse(feature))
                .unwrap_or_else(|| transcript.strand.trim() == "-");
            let external_ranges_1based =
                Self::external_protein_transcript_coding_ranges_1based(transcript);
            let comparison = if let Some(external_source) = external_source {
                Self::build_transcript_protein_comparison(
                    &transcript.transcript_id,
                    &transcript_label,
                    is_reverse,
                    None,
                    &[],
                    Some(Self::build_external_protein_opinion_summary(
                        external_source,
                        transcript,
                    )),
                )
            } else {
                Self::build_transcript_protein_comparison(
                    &transcript.transcript_id,
                    &transcript_label,
                    is_reverse,
                    None,
                    &[],
                    None,
                )
            };
            let mut domains = transcript
                .feature_projections
                .iter()
                .filter_map(|feature_projection| {
                    let normalized_key =
                        Self::normalized_protein_feature_key(&feature_projection.feature_key);
                    if external_source
                        .map(|external_source| {
                            Self::external_protein_feature_hidden_by_default(
                                external_source.source,
                                &normalized_key,
                            )
                        })
                        .unwrap_or(false)
                        && !include_feature_keys.contains(&normalized_key)
                    {
                        return None;
                    }
                    if (!include_feature_keys.is_empty()
                        && !include_feature_keys.contains(&normalized_key))
                        || exclude_feature_keys.contains(&normalized_key)
                    {
                        return None;
                    }
                    Self::projected_external_protein_feature_domain(feature_projection)
                })
                .collect::<Vec<_>>();
            domains.sort_by(|left, right| {
                left.start_aa
                    .cmp(&right.start_aa)
                    .then(left.end_aa.cmp(&right.end_aa))
                    .then(left.name.cmp(&right.name))
            });

            region_ranges_1based.extend(transcript_exons.iter().copied());
            region_ranges_1based.extend(external_ranges_1based.iter().copied());

            transcript_lanes.push(IsoformArchitectureTranscriptLane {
                isoform_id: transcript.transcript_id.clone(),
                label: transcript_label.clone(),
                transcript_id: Some(transcript.transcript_id.clone()),
                transcript_feature_id: transcript.transcript_feature_id,
                strand: if is_reverse {
                    "-".to_string()
                } else {
                    "+".to_string()
                },
                transcript_exons: Self::range_vec_1based_to_splicing(transcript_exons),
                exons: Self::range_vec_1based_to_splicing(external_ranges_1based.clone()),
                introns: Self::range_vec_1based_to_splicing(
                    Self::build_transcript_introns_from_ranges_1based(&external_ranges_1based),
                ),
                mapped: !transcript.aa_segments.is_empty(),
                transactivation_class: None,
                cds_to_protein_segments: transcript
                    .aa_segments
                    .iter()
                    .map(|segment| IsoformArchitectureCdsAaSegment {
                        genomic_start_1based: segment.genomic_start_1based,
                        genomic_end_1based: segment.genomic_end_1based,
                        aa_start: segment.aa_start,
                        aa_end: segment.aa_end,
                    })
                    .collect(),
                note: Some(
                    "Only an external protein opinion was available; transcript-native CDS translation could not be derived."
                        .to_string(),
                ),
            });
            protein_lanes.push(IsoformArchitectureProteinLane {
                isoform_id: transcript.transcript_id.clone(),
                label: transcript_label,
                transcript_id: Some(transcript.transcript_id.clone()),
                expected_length_aa: comparison
                    .external_opinion
                    .as_ref()
                    .and_then(|opinion| opinion.expected_length_aa),
                reference_start_aa: comparison
                    .external_opinion
                    .as_ref()
                    .and_then(|opinion| opinion.reference_start_aa),
                reference_end_aa: comparison
                    .external_opinion
                    .as_ref()
                    .and_then(|opinion| opinion.reference_end_aa),
                domains,
                transactivation_class: None,
                comparison: Some(comparison),
            });
        }

        if transcript_lanes.is_empty() {
            let descriptor = normalized_filter
                .as_deref()
                .map(|filter| format!(" matching transcript '{filter}'"))
                .unwrap_or_default();
            return Err(EngineError {
                code: ErrorCode::NotFound,
                message: format!(
                    "Sequence '{}' has no transcript/protein expert rows{}",
                    seq_id, descriptor
                ),
            });
        }

        let region_start_1based = region_ranges_1based
            .iter()
            .map(|(start, _)| *start)
            .min()
            .unwrap_or(1);
        let region_end_1based = region_ranges_1based
            .iter()
            .map(|(_, end)| *end)
            .max()
            .unwrap_or_else(|| dna.len().max(1));
        let gene_symbol = external_source
            .and_then(|external_source| external_source.gene_symbol_hint.clone())
            .or_else(|| {
                features.iter().find_map(|feature| {
                    Self::first_nonempty_feature_qualifier(
                        feature,
                        &["gene", "gene_id", "locus_tag", "standard_name", "name"],
                    )
                })
            })
            .unwrap_or_else(|| seq_id.to_string());
        let panel_id = external_source
            .map(|external_source| external_source.source_id.clone())
            .unwrap_or_else(|| match normalized_filter.as_deref() {
                Some(filter) => format!("protein_compare:{filter}@{seq_id}"),
                None => format!("protein_compare@{seq_id}"),
            });
        let panel_source = external_source
            .map(|external_source| external_source.panel_source_label.clone())
            .or_else(|| {
                Some(
                    "Transcript-native protein derivation (external protein opinions optional)"
                        .to_string(),
                )
            });

        Ok(IsoformArchitectureExpertView {
            seq_id: seq_id.to_string(),
            panel_id,
            gene_symbol,
            transcript_geometry_mode: IsoformTranscriptGeometryMode::Cds.as_str().to_string(),
            panel_source,
            region_start_1based: region_start_1based.max(1),
            region_end_1based: region_end_1based.max(region_start_1based.max(1)),
            instruction: ISOFORM_ARCHITECTURE_EXPERT_INSTRUCTION.to_string(),
            transcript_lanes,
            protein_lanes,
            warnings,
        })
    }

    pub(super) fn build_uniprot_projection_expert_view(
        &self,
        seq_id: &str,
        projection_id: &str,
        protein_feature_filter: &ProteinFeatureFilter,
    ) -> Result<IsoformArchitectureExpertView, EngineError> {
        let projection = self.get_uniprot_genome_projection(projection_id)?;
        if projection.seq_id != seq_id {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: format!(
                    "UniProt projection '{}' belongs to '{}' rather than '{}'",
                    projection_id, projection.seq_id, seq_id
                ),
            });
        }
        let entry = self.get_uniprot_entry(&projection.entry_id)?;
        let provider = ExternalProteinOpinionProvider::UniprotProjection {
            projection: &projection,
            entry: &entry,
        };
        self.build_optional_external_protein_expert_view(
            seq_id,
            provider.transcript_id_filter(),
            Some(&provider),
            protein_feature_filter,
        )
    }

    pub(super) fn upsert_uniprot_projection(
        &mut self,
        projection: UniprotGenomeProjection,
    ) -> Result<(), EngineError> {
        let mut store = self.read_uniprot_projection_store();
        store
            .projections
            .insert(projection.projection_id.clone(), projection);
        self.write_uniprot_projection_store(store)
    }

    pub(super) fn load_isoform_panel_resource(
        path: &str,
        panel_id_override: Option<&str>,
    ) -> Result<IsoformPanelResource, EngineError> {
        let text = std::fs::read_to_string(path).map_err(|e| EngineError {
            code: ErrorCode::Io,
            message: format!("Could not read isoform panel file '{path}': {e}"),
        })?;
        let mut resource =
            serde_json::from_str::<IsoformPanelResource>(&text).map_err(|e| EngineError {
                code: ErrorCode::InvalidInput,
                message: format!("Could not parse isoform panel JSON in '{path}': {e}"),
            })?;
        if resource.schema.trim().is_empty() {
            resource.schema = ISOFORM_PANEL_RESOURCE_SCHEMA.to_string();
        }
        if resource.schema != ISOFORM_PANEL_RESOURCE_SCHEMA {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: format!(
                    "Unsupported isoform-panel schema '{}' in '{}'; expected '{}'",
                    resource.schema, path, ISOFORM_PANEL_RESOURCE_SCHEMA
                ),
            });
        }
        if let Some(override_id) = panel_id_override.map(str::trim).filter(|v| !v.is_empty()) {
            resource.panel_id = override_id.to_string();
        }
        resource.panel_id = resource.panel_id.trim().to_string();
        if resource.panel_id.is_empty() {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: format!("Isoform panel file '{}' has empty panel_id", path),
            });
        }
        resource.gene_symbol = resource.gene_symbol.trim().to_string();
        if resource.gene_symbol.is_empty() {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: format!("Isoform panel file '{}' has empty gene_symbol", path),
            });
        }
        if resource.isoforms.is_empty() {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: format!("Isoform panel file '{}' has no isoforms", path),
            });
        }
        for (idx, isoform) in resource.isoforms.iter_mut().enumerate() {
            isoform.isoform_id = isoform.isoform_id.trim().to_string();
            if isoform.isoform_id.is_empty() {
                return Err(EngineError {
                    code: ErrorCode::InvalidInput,
                    message: format!(
                        "Isoform panel '{}' contains isoform with empty isoform_id at row {}",
                        resource.panel_id,
                        idx + 1
                    ),
                });
            }
            isoform.label = isoform
                .label
                .take()
                .map(|v| v.trim().to_string())
                .filter(|v| !v.is_empty());
            let mut transcript_ids: Vec<String> = vec![];
            let mut seen_transcripts: BTreeSet<String> = BTreeSet::new();
            for entry in &isoform.transcript_ids {
                let normalized = entry.trim().to_string();
                if normalized.is_empty() {
                    continue;
                }
                let key = normalized.to_ascii_uppercase();
                if seen_transcripts.insert(key) {
                    transcript_ids.push(normalized);
                }
            }
            isoform.transcript_ids = transcript_ids;
            for (domain_idx, domain) in isoform.domains.iter_mut().enumerate() {
                domain.name = domain.name.trim().to_string();
                if domain.name.is_empty() {
                    domain.name = format!("domain_{}", domain_idx + 1);
                }
                if domain.start_aa == 0 || domain.end_aa == 0 || domain.end_aa < domain.start_aa {
                    return Err(EngineError {
                        code: ErrorCode::InvalidInput,
                        message: format!(
                            "Isoform panel '{}' has invalid domain range for isoform '{}' domain '{}': {}..{}",
                            resource.panel_id,
                            isoform.isoform_id,
                            domain.name,
                            domain.start_aa,
                            domain.end_aa
                        ),
                    });
                }
                domain.color_hex = domain
                    .color_hex
                    .take()
                    .map(|v| v.trim().to_string())
                    .filter(|v| !v.is_empty());
            }
            isoform.domains.sort_by(|left, right| {
                left.start_aa
                    .cmp(&right.start_aa)
                    .then(left.end_aa.cmp(&right.end_aa))
            });
        }
        Ok(resource)
    }

    pub(super) fn isoform_validation_issue(
        code: &str,
        message: String,
        isoform_id: Option<&str>,
        transcript_probe: Option<&str>,
        domain_name: Option<&str>,
    ) -> IsoformPanelValidationIssue {
        IsoformPanelValidationIssue {
            severity: "warning".to_string(),
            code: code.to_string(),
            message,
            isoform_id: isoform_id.map(|v| v.to_string()),
            transcript_probe: transcript_probe.map(|v| v.to_string()),
            domain_name: domain_name.map(|v| v.to_string()),
        }
    }

    pub fn validate_isoform_panel_resource(
        path: &str,
        panel_id_override: Option<&str>,
    ) -> Result<IsoformPanelValidationReport, EngineError> {
        let resource = Self::load_isoform_panel_resource(path, panel_id_override)?;
        let color_hex_re = Regex::new(r"^#[0-9A-Fa-f]{6}$").map_err(|e| EngineError {
            code: ErrorCode::Internal,
            message: format!("Could not compile isoform color-hex regex: {e}"),
        })?;

        let mut issues: Vec<IsoformPanelValidationIssue> = vec![];
        let mut isoforms: Vec<IsoformPanelValidationIsoformSummary> = vec![];
        let mut transcript_probe_count = 0usize;
        let mut unique_transcript_probes: BTreeSet<String> = BTreeSet::new();
        let mut domain_count = 0usize;
        let mut isoform_id_buckets: HashMap<String, Vec<String>> = HashMap::new();
        let mut transcript_probe_buckets: HashMap<String, BTreeSet<String>> = HashMap::new();

        for isoform in &resource.isoforms {
            let isoform_key = isoform.isoform_id.to_ascii_uppercase();
            isoform_id_buckets
                .entry(isoform_key)
                .or_default()
                .push(isoform.isoform_id.clone());

            if isoform.transcript_ids.is_empty() {
                issues.push(Self::isoform_validation_issue(
                    "missing_transcript_probe",
                    format!(
                        "Isoform '{}' has no transcript probes; mapping can only rely on label heuristics",
                        isoform.isoform_id
                    ),
                    Some(&isoform.isoform_id),
                    None,
                    None,
                ));
            }
            for transcript_probe in &isoform.transcript_ids {
                let normalized = Self::normalize_transcript_probe(transcript_probe);
                if normalized.is_empty() {
                    continue;
                }
                transcript_probe_count += 1;
                unique_transcript_probes.insert(normalized.clone());
                transcript_probe_buckets
                    .entry(normalized)
                    .or_default()
                    .insert(isoform.isoform_id.clone());
            }

            if isoform.domains.is_empty() {
                issues.push(Self::isoform_validation_issue(
                    "missing_domains",
                    format!(
                        "Isoform '{}' has no domains; protein lane may render without functional landmarks",
                        isoform.isoform_id
                    ),
                    Some(&isoform.isoform_id),
                    None,
                    None,
                ));
            }

            domain_count += isoform.domains.len();
            let mut sorted_domains = isoform.domains.clone();
            sorted_domains.sort_by(|left, right| {
                left.start_aa
                    .cmp(&right.start_aa)
                    .then(left.end_aa.cmp(&right.end_aa))
                    .then(left.name.cmp(&right.name))
            });

            for domain in &sorted_domains {
                if let Some(color_hex) = domain.color_hex.as_deref()
                    && !color_hex_re.is_match(color_hex)
                {
                    issues.push(Self::isoform_validation_issue(
                        "invalid_color_hex",
                        format!(
                            "Isoform '{}' domain '{}' uses non-hex color '{}'; expected '#RRGGBB'",
                            isoform.isoform_id, domain.name, color_hex
                        ),
                        Some(&isoform.isoform_id),
                        None,
                        Some(&domain.name),
                    ));
                }
            }

            for pair in sorted_domains.windows(2) {
                if pair[1].start_aa <= pair[0].end_aa {
                    issues.push(Self::isoform_validation_issue(
                        "overlapping_domains",
                        format!(
                            "Isoform '{}' has overlapping domains '{}' ({}..{}) and '{}' ({}..{})",
                            isoform.isoform_id,
                            pair[0].name,
                            pair[0].start_aa,
                            pair[0].end_aa,
                            pair[1].name,
                            pair[1].start_aa,
                            pair[1].end_aa
                        ),
                        Some(&isoform.isoform_id),
                        None,
                        Some(&pair[1].name),
                    ));
                }
            }

            let max_domain_end_aa = sorted_domains.iter().map(|domain| domain.end_aa).max();
            if let (Some(reference_start), Some(reference_end)) =
                (isoform.reference_start_aa, isoform.reference_end_aa)
                && reference_end < reference_start
            {
                issues.push(Self::isoform_validation_issue(
                    "invalid_reference_range",
                    format!(
                        "Isoform '{}' reference range is invalid (start={} > end={})",
                        isoform.isoform_id, reference_start, reference_end
                    ),
                    Some(&isoform.isoform_id),
                    None,
                    None,
                ));
            }
            if let Some(reference_end) = isoform.reference_end_aa {
                if let Some(max_domain_end) = max_domain_end_aa
                    && reference_end < max_domain_end
                {
                    issues.push(Self::isoform_validation_issue(
                        "reference_end_below_domain_end",
                        format!(
                            "Isoform '{}' reference_end_aa={} is smaller than max domain end {}",
                            isoform.isoform_id, reference_end, max_domain_end
                        ),
                        Some(&isoform.isoform_id),
                        None,
                        None,
                    ));
                }
            } else if let (Some(expected), Some(max_domain_end)) =
                (isoform.expected_length_aa, max_domain_end_aa)
                && expected < max_domain_end
            {
                issues.push(Self::isoform_validation_issue(
                    "expected_length_below_domain_end",
                    format!(
                        "Isoform '{}' expected_length_aa={} is smaller than max domain end {}",
                        isoform.isoform_id, expected, max_domain_end
                    ),
                    Some(&isoform.isoform_id),
                    None,
                    None,
                ));
            }

            isoforms.push(IsoformPanelValidationIsoformSummary {
                isoform_id: isoform.isoform_id.clone(),
                label: isoform
                    .label
                    .as_deref()
                    .map(str::trim)
                    .filter(|v| !v.is_empty())
                    .unwrap_or(isoform.isoform_id.as_str())
                    .to_string(),
                transcript_probe_count: isoform.transcript_ids.len(),
                domain_count: isoform.domains.len(),
                expected_length_aa: isoform.expected_length_aa,
                max_domain_end_aa,
            });
        }

        for isoform_ids in isoform_id_buckets.values() {
            if isoform_ids.len() < 2 {
                continue;
            }
            let mut ids = isoform_ids.clone();
            ids.sort();
            ids.dedup();
            issues.push(Self::isoform_validation_issue(
                "duplicate_isoform_id",
                format!(
                    "Panel '{}' contains duplicate isoform_id entries: {}",
                    resource.panel_id,
                    ids.join(", ")
                ),
                ids.first().map(String::as_str),
                None,
                None,
            ));
        }

        for (normalized_probe, isoform_ids) in &transcript_probe_buckets {
            if isoform_ids.len() < 2 {
                continue;
            }
            let ids = isoform_ids.iter().cloned().collect::<Vec<_>>();
            issues.push(Self::isoform_validation_issue(
                "shared_transcript_probe",
                format!(
                    "Transcript probe '{}' is shared by multiple isoforms: {}",
                    normalized_probe,
                    ids.join(", ")
                ),
                ids.first().map(String::as_str),
                Some(normalized_probe.as_str()),
                None,
            ));
        }

        issues.sort_by(|left, right| {
            left.code
                .cmp(&right.code)
                .then(left.isoform_id.cmp(&right.isoform_id))
                .then(left.transcript_probe.cmp(&right.transcript_probe))
                .then(left.domain_name.cmp(&right.domain_name))
                .then(left.message.cmp(&right.message))
        });

        let issue_count = issues.len();
        let status = if issue_count == 0 {
            "ok".to_string()
        } else {
            "warning".to_string()
        };

        Ok(IsoformPanelValidationReport {
            schema: ISOFORM_PANEL_VALIDATION_REPORT_SCHEMA.to_string(),
            path: path.to_string(),
            panel_id: resource.panel_id,
            gene_symbol: resource.gene_symbol,
            assembly: resource.assembly,
            isoform_count: isoforms.len(),
            transcript_probe_count,
            unique_transcript_probe_count: unique_transcript_probes.len(),
            domain_count,
            issue_count,
            status,
            isoforms,
            issues,
        })
    }

    pub(super) fn build_isoform_architecture_expert_view_from_resource(
        &self,
        seq_id: &str,
        panel_id: &str,
        resource: &IsoformPanelResource,
        strict: bool,
    ) -> Result<IsoformArchitectureExpertView, EngineError> {
        let dna = self
            .state
            .sequences
            .get(seq_id)
            .ok_or_else(|| EngineError {
                code: ErrorCode::NotFound,
                message: format!("Sequence '{seq_id}' not found"),
            })?;
        let features = dna.features();

        #[derive(Clone)]
        struct TranscriptProjection {
            feature_id: usize,
            transcript_id: String,
            is_reverse: bool,
            full_exon_ranges: Vec<(usize, usize)>,
            geometry_ranges: Vec<(usize, usize)>,
            cds_ranges: Vec<(usize, usize)>,
            intron_ranges: Vec<(usize, usize)>,
            used_cds_geometry: bool,
        }

        let mut projections: Vec<TranscriptProjection> = vec![];
        let mut projection_key_index: HashMap<String, Vec<usize>> = HashMap::new();
        let use_cds_geometry =
            resource.transcript_geometry_mode == IsoformTranscriptGeometryMode::Cds;

        let gene_probe = resource.gene_symbol.trim().to_ascii_uppercase();
        let mut has_gene_probe_match = gene_probe.is_empty();
        for (feature_id, feature) in features.iter().enumerate() {
            if !has_gene_probe_match
                && let Some(gene_name) = Self::first_nonempty_feature_qualifier(
                    feature,
                    &[
                        "gene",
                        "gene_id",
                        "locus_tag",
                        "standard_name",
                        "name",
                        "product",
                    ],
                )
            {
                has_gene_probe_match = gene_name.to_ascii_uppercase().contains(&gene_probe);
            }

            if !Self::is_mrna_feature(feature) {
                continue;
            }

            let mut exon_ranges = vec![];
            collect_location_ranges_usize(&feature.location, &mut exon_ranges);
            if exon_ranges.is_empty() {
                let (from, to) = feature.location.find_bounds().map_err(|e| EngineError {
                    code: ErrorCode::InvalidInput,
                    message: format!("Could not parse transcript range: {e}"),
                })?;
                if from >= 0 && to >= 0 {
                    exon_ranges.push((from as usize, to as usize));
                }
            }
            exon_ranges.sort_unstable_by(|a, b| a.0.cmp(&b.0).then(a.1.cmp(&b.1)));
            exon_ranges.retain(|(start, end)| end > start);
            if exon_ranges.is_empty() {
                continue;
            }
            let cds_ranges = Self::feature_qualifier_ranges_0based(feature, "cds_ranges_1based");
            let full_exon_ranges = exon_ranges.clone();
            let (geometry_ranges, used_cds_geometry) = if use_cds_geometry && !cds_ranges.is_empty()
            {
                (cds_ranges.clone(), true)
            } else {
                (exon_ranges, false)
            };
            let mut intron_ranges = vec![];
            for pair in geometry_ranges.windows(2) {
                if pair[1].0 > pair[0].1 {
                    intron_ranges.push((pair[0].1, pair[1].0));
                }
            }
            let projection_idx = projections.len();
            let transcript_id = Self::feature_transcript_id(feature, feature_id);
            projections.push(TranscriptProjection {
                feature_id,
                transcript_id: transcript_id.clone(),
                is_reverse: feature_is_reverse(feature),
                full_exon_ranges,
                geometry_ranges,
                cds_ranges,
                intron_ranges,
                used_cds_geometry,
            });
            for key in Self::feature_transcript_match_keys(feature, feature_id) {
                projection_key_index
                    .entry(key)
                    .or_default()
                    .push(projection_idx);
            }
        }

        let mut warnings: Vec<String> = vec![];
        if projections.is_empty() {
            let message = format!(
                "No mRNA/transcript features available in sequence '{}' for isoform panel '{}'",
                seq_id, panel_id
            );
            if strict {
                return Err(EngineError {
                    code: ErrorCode::NotFound,
                    message,
                });
            }
            warnings.push(message);
        }

        if !has_gene_probe_match {
            let message = format!(
                "No feature-gene qualifier in '{}' matched panel gene_symbol '{}'",
                seq_id, resource.gene_symbol
            );
            if strict {
                return Err(EngineError {
                    code: ErrorCode::InvalidInput,
                    message,
                });
            }
            warnings.push(message);
        }

        let mut mapped_ranges: Vec<(usize, usize)> = vec![];
        let mut transcript_lanes: Vec<IsoformArchitectureTranscriptLane> = vec![];
        let mut protein_lanes: Vec<IsoformArchitectureProteinLane> = vec![];

        for isoform in &resource.isoforms {
            let label = isoform
                .label
                .as_deref()
                .map(str::trim)
                .filter(|v| !v.is_empty())
                .unwrap_or(isoform.isoform_id.as_str())
                .to_string();
            let probes: Vec<String> = isoform
                .transcript_ids
                .iter()
                .map(|entry| Self::normalize_transcript_probe(entry))
                .filter(|entry| !entry.is_empty())
                .collect();
            let mut matched_projection: Option<&TranscriptProjection> = None;
            for probe in &probes {
                if let Some(indices) = projection_key_index.get(probe)
                    && let Some(first) = indices.first()
                {
                    matched_projection = projections.get(*first);
                    if matched_projection.is_some() {
                        break;
                    }
                }
            }

            let (
                mapped,
                transcript_id,
                transcript_feature_id,
                strand,
                transcript_exons,
                exons,
                introns,
                note,
                cds_to_protein_segments,
            ) = if let Some(projection) = matched_projection {
                mapped_ranges.extend(projection.geometry_ranges.iter().copied());
                let geometry_note = (use_cds_geometry && !projection.used_cds_geometry)
                    .then_some("CDS geometry unavailable; fell back to exon geometry".to_string());
                let cds_to_protein_segments = Self::cds_ranges_to_reference_aa_segments(
                    projection.cds_ranges.clone(),
                    projection.is_reverse,
                    isoform.reference_start_aa,
                );
                (
                    true,
                    Some(projection.transcript_id.clone()),
                    Some(projection.feature_id),
                    if projection.is_reverse {
                        "-".to_string()
                    } else {
                        "+".to_string()
                    },
                    Self::range_vec_to_splicing(projection.full_exon_ranges.clone()),
                    Self::range_vec_to_splicing(projection.geometry_ranges.clone()),
                    Self::range_vec_to_splicing(projection.intron_ranges.clone()),
                    geometry_note,
                    cds_to_protein_segments,
                )
            } else {
                let probe_text = if isoform.transcript_ids.is_empty() {
                    "(none)".to_string()
                } else {
                    isoform.transcript_ids.join(", ")
                };
                let message = format!(
                    "Isoform '{}' has no mapped transcript in '{}' (probes: {})",
                    isoform.isoform_id, seq_id, probe_text
                );
                if strict {
                    return Err(EngineError {
                        code: ErrorCode::NotFound,
                        message,
                    });
                }
                warnings.push(message.clone());
                (
                    false,
                    isoform.transcript_ids.first().cloned(),
                    None,
                    "?".to_string(),
                    vec![],
                    vec![],
                    vec![],
                    Some("No transcript mapping found".to_string()),
                    vec![],
                )
            };

            transcript_lanes.push(IsoformArchitectureTranscriptLane {
                isoform_id: isoform.isoform_id.clone(),
                label: label.clone(),
                transcript_id: transcript_id.clone(),
                transcript_feature_id,
                strand,
                transcript_exons,
                exons,
                introns,
                mapped,
                transactivation_class: isoform.transactivation_class.clone(),
                cds_to_protein_segments,
                note,
            });

            let mut domains = isoform
                .domains
                .iter()
                .map(|domain| IsoformArchitectureProteinDomain {
                    name: domain.name.clone(),
                    start_aa: domain.start_aa,
                    end_aa: domain.end_aa,
                    color_hex: domain.color_hex.clone(),
                })
                .collect::<Vec<_>>();
            domains.sort_by(|left, right| {
                left.start_aa
                    .cmp(&right.start_aa)
                    .then(left.end_aa.cmp(&right.end_aa))
            });
            let expected_length_aa = isoform.expected_length_aa.or_else(|| {
                domains
                    .iter()
                    .map(|domain| domain.end_aa)
                    .max()
                    .filter(|max_end| *max_end > 0)
            });
            protein_lanes.push(IsoformArchitectureProteinLane {
                isoform_id: isoform.isoform_id.clone(),
                label,
                transcript_id,
                expected_length_aa,
                reference_start_aa: isoform.reference_start_aa,
                reference_end_aa: isoform.reference_end_aa,
                domains,
                transactivation_class: isoform.transactivation_class.clone(),
                comparison: None,
            });
        }

        let (region_start_1based, region_end_1based) = if mapped_ranges.is_empty() {
            (1, dna.len().max(1))
        } else {
            let start0 = mapped_ranges
                .iter()
                .map(|(start, _)| *start)
                .min()
                .unwrap_or(0);
            let end0 = mapped_ranges
                .iter()
                .map(|(_, end)| *end)
                .max()
                .unwrap_or(start0 + 1);
            (start0 + 1, end0.max(start0 + 1))
        };

        Ok(IsoformArchitectureExpertView {
            seq_id: seq_id.to_string(),
            panel_id: panel_id.to_string(),
            gene_symbol: resource.gene_symbol.clone(),
            transcript_geometry_mode: resource.transcript_geometry_mode.as_str().to_string(),
            panel_source: resource.source.clone(),
            region_start_1based,
            region_end_1based,
            instruction: ISOFORM_ARCHITECTURE_EXPERT_INSTRUCTION.to_string(),
            transcript_lanes,
            protein_lanes,
            warnings,
        })
    }

    pub(super) fn build_isoform_architecture_expert_view(
        &self,
        seq_id: &str,
        panel_id: &str,
    ) -> Result<IsoformArchitectureExpertView, EngineError> {
        let record = self.get_isoform_panel_record(seq_id, panel_id)?;
        self.build_isoform_architecture_expert_view_from_resource(
            seq_id,
            &record.panel_id,
            &record.resource,
            false,
        )
    }

    pub(super) fn splice_boundary_markers_for_introns(
        dna: &DNAsequence,
        transcript_feature_id: usize,
        transcript_id: &str,
        is_reverse: bool,
        introns: &[(usize, usize)],
    ) -> Vec<SplicingBoundaryMarker> {
        let mut out = Vec::new();
        for (start, end) in introns {
            if *end <= *start || end - start < 2 {
                continue;
            }
            if is_reverse {
                let donor_raw = Self::sequence_slice_upper(dna, end.saturating_sub(2), *end);
                let acceptor_raw = Self::sequence_slice_upper(dna, *start, (start + 2).min(*end));
                let donor_rc = Self::reverse_complement_bytes(donor_raw.as_bytes());
                let acceptor_rc = Self::reverse_complement_bytes(acceptor_raw.as_bytes());
                let donor_motif = String::from_utf8_lossy(&donor_rc).to_ascii_uppercase();
                let acceptor_motif = String::from_utf8_lossy(&acceptor_rc).to_ascii_uppercase();
                let (canonical_pair, paired_motif_signature, motif_class, annotation) =
                    Self::classify_splice_boundary_pair(&donor_motif, &acceptor_motif);
                out.push(SplicingBoundaryMarker {
                    transcript_feature_id,
                    transcript_id: transcript_id.to_string(),
                    side: "donor".to_string(),
                    position_1based: *end,
                    canonical: donor_motif == "GT",
                    motif_2bp: donor_motif,
                    canonical_pair,
                    partner_position_1based: start + 1,
                    paired_motif_signature: paired_motif_signature.clone(),
                    motif_class: motif_class.to_string(),
                    annotation: annotation.to_string(),
                });
                out.push(SplicingBoundaryMarker {
                    transcript_feature_id,
                    transcript_id: transcript_id.to_string(),
                    side: "acceptor".to_string(),
                    position_1based: start + 1,
                    canonical: acceptor_motif == "AG",
                    motif_2bp: acceptor_motif,
                    canonical_pair,
                    partner_position_1based: *end,
                    paired_motif_signature,
                    motif_class: motif_class.to_string(),
                    annotation: annotation.to_string(),
                });
            } else {
                let donor_motif = Self::sequence_slice_upper(dna, *start, (start + 2).min(*end));
                let acceptor_motif = Self::sequence_slice_upper(dna, end.saturating_sub(2), *end);
                let (canonical_pair, paired_motif_signature, motif_class, annotation) =
                    Self::classify_splice_boundary_pair(&donor_motif, &acceptor_motif);
                out.push(SplicingBoundaryMarker {
                    transcript_feature_id,
                    transcript_id: transcript_id.to_string(),
                    side: "donor".to_string(),
                    position_1based: start + 1,
                    canonical: donor_motif == "GT",
                    motif_2bp: donor_motif,
                    canonical_pair,
                    partner_position_1based: *end,
                    paired_motif_signature: paired_motif_signature.clone(),
                    motif_class: motif_class.to_string(),
                    annotation: annotation.to_string(),
                });
                out.push(SplicingBoundaryMarker {
                    transcript_feature_id,
                    transcript_id: transcript_id.to_string(),
                    side: "acceptor".to_string(),
                    position_1based: *end,
                    canonical: acceptor_motif == "AG",
                    motif_2bp: acceptor_motif,
                    canonical_pair,
                    partner_position_1based: start + 1,
                    paired_motif_signature,
                    motif_class: motif_class.to_string(),
                    annotation: annotation.to_string(),
                });
            }
        }
        out
    }

    pub(crate) fn splice_intron_signals_for_introns(
        dna: &DNAsequence,
        transcript_feature_id: usize,
        transcript_id: &str,
        is_reverse: bool,
        introns: &[(usize, usize)],
    ) -> Vec<SplicingIntronSignal> {
        let mut out = Vec::new();
        for (start, end) in introns {
            if *end <= *start || end - start < 8 {
                continue;
            }
            let oriented_intron = if is_reverse {
                let raw = Self::sequence_slice_upper(dna, *start, *end);
                let rc = Self::reverse_complement_bytes(raw.as_bytes());
                String::from_utf8_lossy(&rc).to_ascii_uppercase()
            } else {
                Self::sequence_slice_upper(dna, *start, *end)
            };
            let intron_len = oriented_intron.len();
            let donor_position_1based = if is_reverse { *end } else { start + 1 };
            let acceptor_position_1based = if is_reverse { start + 1 } else { *end };
            let (
                branchpoint_position_1based,
                branchpoint_motif,
                branchpoint_score,
                branchpoint_annotation,
            ) = match Self::detect_branchpoint_like_site(&oriented_intron) {
                Some((offset_0based, motif, score, annotation)) => (
                    Some(Self::oriented_offset_to_genomic_1based(
                        *start, *end, is_reverse, offset_0based,
                    )),
                    motif,
                    score,
                    annotation.to_string(),
                ),
                None => (
                    None,
                    String::new(),
                    0.0,
                    "No strong branchpoint-like adenine was detected in the usual acceptor-proximal heuristic window.".to_string(),
                ),
            };
            let (
                polypyrimidine_start_1based,
                polypyrimidine_end_1based,
                polypyrimidine_fraction,
                polypyrimidine_annotation,
            ) = match Self::detect_polypyrimidine_tract(&oriented_intron) {
                Some((tract_start_0based, tract_end_0based, fraction, annotation)) => {
                    let start_1based = Self::oriented_offset_to_genomic_1based(
                        *start,
                        *end,
                        is_reverse,
                        tract_start_0based,
                    );
                    let end_1based = Self::oriented_range_end_to_genomic_1based(
                        *start,
                        *end,
                        is_reverse,
                        tract_end_0based,
                    );
                    (
                        Some(start_1based.min(end_1based)),
                        Some(start_1based.max(end_1based)),
                        fraction,
                        annotation.to_string(),
                    )
                }
                None => (
                    None,
                    None,
                    0.0,
                    "No clear acceptor-proximal polypyrimidine-rich tract was detected by the current heuristic.".to_string(),
                ),
            };
            out.push(SplicingIntronSignal {
                transcript_feature_id,
                transcript_id: transcript_id.to_string(),
                donor_position_1based,
                acceptor_position_1based,
                intron_length_bp: intron_len,
                branchpoint_position_1based,
                branchpoint_motif,
                branchpoint_score,
                branchpoint_annotation,
                polypyrimidine_start_1based,
                polypyrimidine_end_1based,
                polypyrimidine_fraction,
                polypyrimidine_annotation,
            });
        }
        out
    }

    pub(super) fn classify_splice_boundary_pair(
        donor_motif: &str,
        acceptor_motif: &str,
    ) -> (bool, String, &'static str, &'static str) {
        let donor = donor_motif.trim().to_ascii_uppercase();
        let acceptor = acceptor_motif.trim().to_ascii_uppercase();
        let signature = format!("{donor}-{acceptor}");
        match (donor.as_str(), acceptor.as_str()) {
            ("GT", "AG") => (
                true,
                signature,
                "gt_ag_major_canonical",
                "Canonical GT-AG intron; the dominant major/U2-type splice-site motif class.",
            ),
            ("GC", "AG") => (
                false,
                signature,
                "gc_ag_major_noncanonical",
                "GC-AG intron; a known non-canonical major/U2-type splice-site motif class.",
            ),
            ("AT", "AC") => (
                false,
                signature,
                "at_ac_minor_u12_like",
                "AT-AC intron; the classic minor/U12-type splice-site motif class.",
            ),
            ("AT", "AG") => (
                false,
                signature,
                "at_ag_minor_u12_like",
                "AT-AG intron; a known minor/U12-like splice-site motif class.",
            ),
            _ => (
                false,
                signature,
                "other_noncanonical",
                "Non-canonical splice-site motif pair; may reflect rarer biology, incomplete annotation, or sequence/assembly issues.",
            ),
        }
    }

    fn oriented_offset_to_genomic_1based(
        intron_start_0based: usize,
        intron_end_0based: usize,
        is_reverse: bool,
        offset_0based: usize,
    ) -> usize {
        if is_reverse {
            intron_end_0based.saturating_sub(offset_0based)
        } else {
            intron_start_0based + offset_0based + 1
        }
    }

    fn oriented_range_end_to_genomic_1based(
        intron_start_0based: usize,
        intron_end_0based: usize,
        is_reverse: bool,
        end_offset_exclusive_0based: usize,
    ) -> usize {
        if is_reverse {
            intron_end_0based.saturating_sub(end_offset_exclusive_0based.saturating_sub(1))
        } else {
            intron_start_0based + end_offset_exclusive_0based
        }
    }

    fn detect_branchpoint_like_site(
        oriented_intron: &str,
    ) -> Option<(usize, String, f32, &'static str)> {
        let bases = oriented_intron.as_bytes();
        if bases.len() < 20 {
            return None;
        }
        let candidate_start = bases.len().saturating_sub(40).max(2);
        let candidate_end_exclusive = bases.len().saturating_sub(18).min(bases.len());
        if candidate_start >= candidate_end_exclusive {
            return None;
        }
        let mut best: Option<(usize, String, f32, &'static str)> = None;
        for idx in candidate_start..candidate_end_exclusive {
            let base = bases[idx].to_ascii_uppercase();
            if base != b'A' {
                continue;
            }
            let mut score = 1.5_f32;
            if idx >= 1 && Self::is_pyrimidine(bases[idx - 1]) {
                score += 1.2;
            }
            if idx >= 2 && Self::is_pyrimidine(bases[idx - 2]) {
                score += 1.0;
            }
            if idx >= 3 && Self::is_pyrimidine(bases[idx - 3]) {
                score += 0.7;
            }
            if idx + 1 < bases.len() && Self::is_purine(bases[idx + 1]) {
                score += 0.6;
            }
            if idx + 2 < bases.len() && Self::is_pyrimidine(bases[idx + 2]) {
                score += 0.4;
            }
            let motif_start = idx.saturating_sub(2);
            let motif_end = (idx + 3).min(bases.len());
            let motif = oriented_intron[motif_start..motif_end].to_string();
            let annotation = if score >= 4.3 {
                "Strong branchpoint-like adenine in the usual 18-40 nt acceptor-proximal window (heuristic, not a splice predictor)."
            } else if score >= 3.3 {
                "Moderate branchpoint-like adenine in the usual acceptor-proximal window (heuristic only)."
            } else {
                "Weak branchpoint-like adenine in the usual acceptor-proximal window; inspect manually before over-interpreting."
            };
            let replace = best
                .as_ref()
                .map(|(_, _, best_score, _)| score > *best_score)
                .unwrap_or(true);
            if replace {
                best = Some((idx, motif, score, annotation));
            }
        }
        best
    }

    fn detect_polypyrimidine_tract(
        oriented_intron: &str,
    ) -> Option<(usize, usize, f32, &'static str)> {
        let bases = oriented_intron.as_bytes();
        if bases.len() < 12 {
            return None;
        }
        let window_lengths = [14usize, 12, 10, 8];
        let region_start = bases.len().saturating_sub(35);
        let region_end_exclusive = bases.len().saturating_sub(2);
        if region_end_exclusive <= region_start + 6 {
            return None;
        }
        let mut best: Option<(usize, usize, f32, usize, f32)> = None;
        for window_len in window_lengths {
            if window_len > region_end_exclusive.saturating_sub(region_start) {
                continue;
            }
            for start in region_start..=region_end_exclusive.saturating_sub(window_len) {
                let end = start + window_len;
                let pyrimidine_count = bases[start..end]
                    .iter()
                    .filter(|base| Self::is_pyrimidine(**base))
                    .count();
                let fraction = pyrimidine_count as f32 / window_len as f32;
                let distance_bonus =
                    (start.saturating_sub(region_start) as f32) / (bases.len().max(1) as f32);
                let effective_score = fraction + distance_bonus * 0.08;
                let replace = best
                    .as_ref()
                    .map(
                        |(_, _, best_fraction, best_window_len, best_effective_score)| {
                            effective_score > *best_effective_score
                                || (effective_score == *best_effective_score
                                    && (fraction > *best_fraction
                                        || (fraction == *best_fraction
                                            && window_len > *best_window_len)))
                        },
                    )
                    .unwrap_or(true);
                if replace {
                    best = Some((start, end, fraction, window_len, effective_score));
                }
            }
        }
        let (start, end, fraction, _window_len, _effective_score) = best?;
        let annotation = if fraction >= 0.85 {
            "Strong acceptor-proximal polypyrimidine tract by simple pyrimidine-density heuristic."
        } else if fraction >= 0.75 {
            "Moderate acceptor-proximal polypyrimidine tract by simple pyrimidine-density heuristic."
        } else if fraction >= 0.65 {
            "Weak pyrimidine-rich tract near the acceptor; review manually before treating it as strong evidence."
        } else {
            return None;
        };
        Some((start, end, fraction, annotation))
    }

    fn is_pyrimidine(base: u8) -> bool {
        matches!(base.to_ascii_uppercase(), b'C' | b'T' | b'U')
    }

    fn is_purine(base: u8) -> bool {
        matches!(base.to_ascii_uppercase(), b'A' | b'G')
    }

    pub(super) fn range_intersection_0based(
        left: (usize, usize),
        right: (usize, usize),
    ) -> Option<(usize, usize)> {
        let start = left.0.max(right.0);
        let end = left.1.min(right.1);
        (end > start).then_some((start, end))
    }

    pub(super) fn exon_cds_phases_for_transcript(
        feature: &gb_io::seq::Feature,
        exon_ranges: &[(usize, usize)],
        is_reverse: bool,
    ) -> Vec<SplicingExonCdsPhase> {
        let mut phases = exon_ranges
            .iter()
            .map(|(start, end)| SplicingExonCdsPhase {
                start_1based: start + 1,
                end_1based: *end,
                left_cds_phase: None,
                right_cds_phase: None,
            })
            .collect::<Vec<_>>();
        if exon_ranges.is_empty() {
            return phases;
        }
        let mut cds_ranges = Self::feature_qualifier_ranges_0based(feature, "cds_ranges_1based");
        if cds_ranges.is_empty() {
            return phases;
        }
        cds_ranges.sort_unstable_by(|a, b| a.0.cmp(&b.0).then(a.1.cmp(&b.1)));
        cds_ranges.dedup();

        let mut consumed_cds_bp = 0usize;
        let exon_indices = if is_reverse {
            (0..exon_ranges.len()).rev().collect::<Vec<_>>()
        } else {
            (0..exon_ranges.len()).collect::<Vec<_>>()
        };

        for exon_idx in exon_indices {
            let exon = exon_ranges[exon_idx];
            let mut coding_segments = cds_ranges
                .iter()
                .filter_map(|cds| Self::range_intersection_0based(exon, *cds))
                .collect::<Vec<_>>();
            if coding_segments.is_empty() {
                continue;
            }
            if is_reverse {
                coding_segments.sort_unstable_by(|a, b| b.0.cmp(&a.0).then(b.1.cmp(&a.1)));
            } else {
                coding_segments.sort_unstable_by(|a, b| a.0.cmp(&b.0).then(a.1.cmp(&b.1)));
            }

            for (seg_start, seg_end) in coding_segments {
                let seg_len = seg_end.saturating_sub(seg_start);
                if seg_len == 0 {
                    continue;
                }
                let entry_phase = (consumed_cds_bp % 3) as u8;
                let exit_phase = ((consumed_cds_bp + seg_len - 1) % 3) as u8;
                if is_reverse {
                    if seg_end == exon.1 && phases[exon_idx].right_cds_phase.is_none() {
                        phases[exon_idx].right_cds_phase = Some(entry_phase);
                    }
                    if seg_start == exon.0 {
                        phases[exon_idx].left_cds_phase = Some(exit_phase);
                    }
                } else {
                    if seg_start == exon.0 && phases[exon_idx].left_cds_phase.is_none() {
                        phases[exon_idx].left_cds_phase = Some(entry_phase);
                    }
                    if seg_end == exon.1 {
                        phases[exon_idx].right_cds_phase = Some(exit_phase);
                    }
                }
                consumed_cds_bp += seg_len;
            }
        }
        phases
    }

    pub(super) fn build_splicing_expert_view(
        &self,
        seq_id: &str,
        feature_id: usize,
        scope: SplicingScopePreset,
    ) -> Result<SplicingExpertView, EngineError> {
        let dna = self
            .state
            .sequences
            .get(seq_id)
            .ok_or_else(|| EngineError {
                code: ErrorCode::NotFound,
                message: format!("Sequence '{seq_id}' not found"),
            })?;
        let features = dna.features();
        let target_feature = features.get(feature_id).ok_or_else(|| EngineError {
            code: ErrorCode::NotFound,
            message: format!(
                "Feature id '{}' was not found in sequence '{}'",
                feature_id, seq_id
            ),
        })?;
        if !Self::is_splicing_seed_feature(target_feature) {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: format!(
                    "Feature '{}' in '{}' is not an mRNA/transcript/ncRNA/misc_RNA/exon/gene/CDS feature and cannot seed a splicing view",
                    feature_id, seq_id
                ),
            });
        }
        let target_group = Self::splicing_group_label(target_feature, feature_id);
        let target_is_reverse = feature_is_reverse(target_feature);
        let restrict_to_target_group = scope.restrict_to_target_group();
        let restrict_to_target_strand = scope.restrict_to_target_strand();

        let mut roi_ranges_0based: Vec<(usize, usize)> = Vec::new();
        for (idx, feature) in features.iter().enumerate() {
            if Self::splicing_group_label(feature, idx) != target_group {
                continue;
            }
            if !Self::is_splicing_transcript_feature(feature) && !Self::is_exon_feature(feature) {
                continue;
            }
            collect_location_ranges_usize(&feature.location, &mut roi_ranges_0based);
        }
        if roi_ranges_0based.is_empty() {
            let (from, to) = target_feature
                .location
                .find_bounds()
                .map_err(|e| EngineError {
                    code: ErrorCode::InvalidInput,
                    message: format!("Could not parse target feature range: {e}"),
                })?;
            if from >= 0 && to >= 0 {
                roi_ranges_0based.push((from as usize, to as usize));
            }
        }
        let roi_start_0based = roi_ranges_0based
            .iter()
            .map(|(start, _)| *start)
            .min()
            .unwrap_or(0);
        let roi_end_0based_exclusive = roi_ranges_0based
            .iter()
            .map(|(_, end)| *end)
            .max()
            .unwrap_or(roi_start_0based.saturating_add(1));

        #[derive(Clone)]
        struct TranscriptWork {
            feature_id: usize,
            transcript_id: String,
            label: String,
            is_reverse: bool,
            exon_ranges: Vec<(usize, usize)>,
            exon_cds_phases: Vec<SplicingExonCdsPhase>,
            intron_ranges: Vec<(usize, usize)>,
        }

        let mut transcripts: Vec<TranscriptWork> = Vec::new();
        for (idx, feature) in features.iter().enumerate() {
            if !Self::is_splicing_transcript_feature(feature) {
                continue;
            }
            if restrict_to_target_group && Self::splicing_group_label(feature, idx) != target_group
            {
                continue;
            }
            let is_reverse = feature_is_reverse(feature);
            if restrict_to_target_strand && is_reverse != target_is_reverse {
                continue;
            }
            if !restrict_to_target_group {
                let mut feature_ranges = Vec::<(usize, usize)>::new();
                collect_location_ranges_usize(&feature.location, &mut feature_ranges);
                let overlaps_roi = feature_ranges.iter().any(|(start, end)| {
                    *end > roi_start_0based && *start < roi_end_0based_exclusive
                });
                if !overlaps_roi {
                    continue;
                }
            }
            let mut exon_ranges = vec![];
            collect_location_ranges_usize(&feature.location, &mut exon_ranges);
            if exon_ranges.is_empty() {
                let (from, to) = feature.location.find_bounds().map_err(|e| EngineError {
                    code: ErrorCode::InvalidInput,
                    message: format!("Could not parse transcript range: {e}"),
                })?;
                if from >= 0 && to >= 0 {
                    exon_ranges.push((from as usize, to as usize));
                }
            }
            exon_ranges.sort_unstable_by(|a, b| a.0.cmp(&b.0).then(a.1.cmp(&b.1)));
            exon_ranges.retain(|(start, end)| end > start);
            if exon_ranges.is_empty() {
                continue;
            }
            let mut introns: Vec<(usize, usize)> = vec![];
            for pair in exon_ranges.windows(2) {
                let left = pair[0];
                let right = pair[1];
                if right.0 > left.1 {
                    introns.push((left.1, right.0));
                }
            }
            let exon_cds_phases =
                Self::exon_cds_phases_for_transcript(feature, &exon_ranges, is_reverse);
            transcripts.push(TranscriptWork {
                feature_id: idx,
                transcript_id: Self::feature_transcript_id(feature, idx),
                label: Self::feature_display_label(feature, idx),
                is_reverse,
                exon_ranges,
                exon_cds_phases,
                intron_ranges: introns,
            });
        }

        if transcripts.is_empty() {
            return Err(EngineError {
                code: ErrorCode::NotFound,
                message: format!(
                    "No transcript-like RNA features found for splicing group '{}' in '{}'",
                    target_group, seq_id
                ),
            });
        }

        transcripts.sort_by(|a, b| a.transcript_id.cmp(&b.transcript_id));

        let mut exon_support: HashMap<(usize, usize), HashSet<usize>> = HashMap::new();
        let mut junction_support: HashMap<(usize, usize), HashSet<usize>> = HashMap::new();
        for transcript in &transcripts {
            for exon in &transcript.exon_ranges {
                exon_support
                    .entry(*exon)
                    .or_default()
                    .insert(transcript.feature_id);
            }
            for intron in &transcript.intron_ranges {
                junction_support
                    .entry(*intron)
                    .or_default()
                    .insert(transcript.feature_id);
            }
        }

        let mut unique_exons: Vec<(usize, usize)> = exon_support.keys().copied().collect();
        unique_exons.sort_unstable_by(|a, b| a.0.cmp(&b.0).then(a.1.cmp(&b.1)));
        let region_start = unique_exons
            .iter()
            .map(|(start, _)| *start)
            .min()
            .unwrap_or(0);
        let region_end = unique_exons
            .iter()
            .map(|(_, end)| *end)
            .max()
            .unwrap_or(region_start + 1);

        let transcript_count = transcripts.len();
        let unique_exon_summaries: Vec<SplicingExonSummary> = unique_exons
            .iter()
            .map(|(start, end)| {
                let support = exon_support
                    .get(&(*start, *end))
                    .map(|set| set.len())
                    .unwrap_or(0);
                SplicingExonSummary {
                    start_1based: start + 1,
                    end_1based: *end,
                    support_transcript_count: support,
                    constitutive: support == transcript_count,
                }
            })
            .collect();

        let transcript_lanes: Vec<SplicingTranscriptLane> = transcripts
            .iter()
            .map(|transcript| SplicingTranscriptLane {
                transcript_feature_id: transcript.feature_id,
                transcript_id: transcript.transcript_id.clone(),
                label: transcript.label.clone(),
                strand: if transcript.is_reverse {
                    "-".to_string()
                } else {
                    "+".to_string()
                },
                exons: Self::range_vec_to_splicing(transcript.exon_ranges.clone()),
                exon_cds_phases: transcript.exon_cds_phases.clone(),
                introns: Self::range_vec_to_splicing(transcript.intron_ranges.clone()),
                has_target_feature: transcript.feature_id == feature_id,
            })
            .collect();

        let matrix_rows: Vec<SplicingMatrixRow> = transcripts
            .iter()
            .map(|transcript| {
                let set = transcript
                    .exon_ranges
                    .iter()
                    .copied()
                    .collect::<HashSet<(usize, usize)>>();
                SplicingMatrixRow {
                    transcript_feature_id: transcript.feature_id,
                    transcript_id: transcript.transcript_id.clone(),
                    label: transcript.label.clone(),
                    exon_presence: unique_exons.iter().map(|exon| set.contains(exon)).collect(),
                }
            })
            .collect();

        let mut boundaries = Vec::new();
        let mut intron_signals = Vec::new();
        for transcript in &transcripts {
            boundaries.extend(Self::splice_boundary_markers_for_introns(
                dna,
                transcript.feature_id,
                &transcript.transcript_id,
                transcript.is_reverse,
                &transcript.intron_ranges,
            ));
            intron_signals.extend(Self::splice_intron_signals_for_introns(
                dna,
                transcript.feature_id,
                &transcript.transcript_id,
                transcript.is_reverse,
                &transcript.intron_ranges,
            ));
        }
        boundaries.sort_by(|a, b| {
            a.position_1based
                .cmp(&b.position_1based)
                .then_with(|| a.side.cmp(&b.side))
                .then_with(|| a.transcript_feature_id.cmp(&b.transcript_feature_id))
        });
        intron_signals.sort_by(|a, b| {
            a.donor_position_1based
                .cmp(&b.donor_position_1based)
                .then_with(|| a.acceptor_position_1based.cmp(&b.acceptor_position_1based))
                .then_with(|| a.transcript_feature_id.cmp(&b.transcript_feature_id))
        });

        let mut junctions: Vec<SplicingJunctionArc> = junction_support
            .into_iter()
            .map(|((donor0, acceptor0), ids)| {
                let mut transcript_feature_ids = ids.into_iter().collect::<Vec<_>>();
                transcript_feature_ids.sort_unstable();
                SplicingJunctionArc {
                    donor_1based: donor0 + 1,
                    acceptor_1based: acceptor0,
                    support_transcript_count: transcript_feature_ids.len(),
                    transcript_feature_ids,
                }
            })
            .collect();
        junctions.sort_by(|a, b| {
            a.donor_1based
                .cmp(&b.donor_1based)
                .then(a.acceptor_1based.cmp(&b.acceptor_1based))
        });

        let mut exon_skip_details = Vec::new();
        let mut intron_retention_details = Vec::new();
        let mut alt5_details = Vec::new();
        let mut alt3_details = Vec::new();
        let mut mex_details = Vec::new();

        let exon_set_by_transcript: HashMap<usize, HashSet<(usize, usize)>> = transcripts
            .iter()
            .map(|transcript| {
                (
                    transcript.feature_id,
                    transcript.exon_ranges.iter().copied().collect(),
                )
            })
            .collect();

        let mut alt5_groups = 0usize;
        let mut alt3_groups = 0usize;
        let mut grouped_by_start: HashMap<usize, Vec<(usize, usize)>> = HashMap::new();
        let mut grouped_by_end: HashMap<usize, Vec<(usize, usize)>> = HashMap::new();
        for exon in &unique_exons {
            grouped_by_start.entry(exon.0).or_default().push(*exon);
            grouped_by_end.entry(exon.1).or_default().push(*exon);
        }
        if target_is_reverse {
            for (start, group) in grouped_by_start {
                let distinct = group
                    .iter()
                    .map(|(_, end)| *end)
                    .collect::<HashSet<_>>()
                    .len();
                if distinct > 1 {
                    alt5_groups += 1;
                    alt5_details.push(format!(
                        "shared start={} with {} alternative 5' exon end(s)",
                        start + 1,
                        distinct
                    ));
                }
            }
            for (end, group) in grouped_by_end {
                let distinct = group
                    .iter()
                    .map(|(start, _)| *start)
                    .collect::<HashSet<_>>()
                    .len();
                if distinct > 1 {
                    alt3_groups += 1;
                    alt3_details.push(format!(
                        "shared end={} with {} alternative 3' exon start(s)",
                        end, distinct
                    ));
                }
            }
        } else {
            for (start, group) in grouped_by_start {
                let distinct = group
                    .iter()
                    .map(|(_, end)| *end)
                    .collect::<HashSet<_>>()
                    .len();
                if distinct > 1 {
                    alt3_groups += 1;
                    alt3_details.push(format!(
                        "shared start={} with {} alternative 3' exon end(s)",
                        start + 1,
                        distinct
                    ));
                }
            }
            for (end, group) in grouped_by_end {
                let distinct = group
                    .iter()
                    .map(|(start, _)| *start)
                    .collect::<HashSet<_>>()
                    .len();
                if distinct > 1 {
                    alt5_groups += 1;
                    alt5_details.push(format!(
                        "shared end={} with {} alternative 5' exon start(s)",
                        end, distinct
                    ));
                }
            }
        }

        let mut exon_skip_count = 0usize;
        for transcript in &transcripts {
            for (intron_start, intron_end) in &transcript.intron_ranges {
                let skipped = unique_exons
                    .iter()
                    .filter(|(exon_start, exon_end)| {
                        *exon_start >= *intron_start
                            && *exon_end <= *intron_end
                            && !transcript
                                .exon_ranges
                                .iter()
                                .any(|own| own.0 == *exon_start && own.1 == *exon_end)
                    })
                    .count();
                if skipped > 0 {
                    exon_skip_count += skipped;
                    exon_skip_details.push(format!(
                        "{} junction {}..{} skips {} exon(s)",
                        transcript.transcript_id,
                        intron_start + 1,
                        intron_end,
                        skipped
                    ));
                }
            }
        }

        let mut intron_retention_count = 0usize;
        for transcript in &transcripts {
            for (intron_start, intron_end) in &transcript.intron_ranges {
                let retained_by = transcripts.iter().find_map(|other| {
                    if other.feature_id == transcript.feature_id {
                        return None;
                    }
                    if other.exon_ranges.iter().any(|(exon_start, exon_end)| {
                        exon_start <= intron_start && exon_end >= intron_end
                    }) {
                        Some(other.transcript_id.clone())
                    } else {
                        None
                    }
                });
                if let Some(other_id) = retained_by {
                    intron_retention_count += 1;
                    intron_retention_details.push(format!(
                        "{} intron {}..{} retained in {}",
                        transcript.transcript_id,
                        intron_start + 1,
                        intron_end,
                        other_id
                    ));
                }
            }
        }

        let mut mex_count = 0usize;
        for left_idx in 0..unique_exons.len() {
            for right_idx in (left_idx + 1)..unique_exons.len() {
                let left = unique_exons[left_idx];
                let right = unique_exons[right_idx];
                if left.1 <= right.0 || right.1 <= left.0 {
                    continue;
                }
                let left_support = exon_support.get(&left).cloned().unwrap_or_default();
                let right_support = exon_support.get(&right).cloned().unwrap_or_default();
                if left_support.is_empty()
                    || right_support.is_empty()
                    || !left_support.is_disjoint(&right_support)
                {
                    continue;
                }
                let co_occurs = transcripts.iter().any(|transcript| {
                    let set = exon_set_by_transcript
                        .get(&transcript.feature_id)
                        .cloned()
                        .unwrap_or_default();
                    set.contains(&left) && set.contains(&right)
                });
                if !co_occurs {
                    mex_count += 1;
                    mex_details.push(format!(
                        "{}..{} and {}..{} appear mutually exclusive",
                        left.0 + 1,
                        left.1,
                        right.0 + 1,
                        right.1
                    ));
                }
            }
        }

        let alt_exon_count = unique_exon_summaries
            .iter()
            .filter(|exon| !exon.constitutive)
            .count();

        let events = vec![
            SplicingEventSummary {
                event_type: "alternative_exon".to_string(),
                count: alt_exon_count,
                details: unique_exon_summaries
                    .iter()
                    .filter(|exon| !exon.constitutive)
                    .take(6)
                    .map(|exon| {
                        format!(
                            "{}..{} support={}/{}",
                            exon.start_1based,
                            exon.end_1based,
                            exon.support_transcript_count,
                            transcript_count
                        )
                    })
                    .collect(),
            },
            SplicingEventSummary {
                event_type: "alternative_5prime".to_string(),
                count: alt5_groups,
                details: alt5_details.into_iter().take(6).collect(),
            },
            SplicingEventSummary {
                event_type: "alternative_3prime".to_string(),
                count: alt3_groups,
                details: alt3_details.into_iter().take(6).collect(),
            },
            SplicingEventSummary {
                event_type: "exon_skipping".to_string(),
                count: exon_skip_count,
                details: exon_skip_details.into_iter().take(6).collect(),
            },
            SplicingEventSummary {
                event_type: "intron_retention".to_string(),
                count: intron_retention_count,
                details: intron_retention_details.into_iter().take(6).collect(),
            },
            SplicingEventSummary {
                event_type: "mutually_exclusive_exons".to_string(),
                count: mex_count,
                details: mex_details.into_iter().take(6).collect(),
            },
        ];

        Ok(SplicingExpertView {
            seq_id: seq_id.to_string(),
            target_feature_id: feature_id,
            group_label: target_group,
            strand: if target_is_reverse {
                "-".to_string()
            } else {
                "+".to_string()
            },
            region_start_1based: region_start + 1,
            region_end_1based: region_end,
            transcript_count,
            unique_exon_count: unique_exons.len(),
            instruction: SPLICING_EXPERT_INSTRUCTION.to_string(),
            transcripts: transcript_lanes,
            unique_exons: unique_exon_summaries,
            matrix_rows,
            boundaries,
            intron_signals,
            junctions,
            events,
        })
    }

    pub fn inspect_feature_expert(
        &self,
        seq_id: &str,
        target: &FeatureExpertTarget,
    ) -> Result<FeatureExpertView, EngineError> {
        match target {
            FeatureExpertTarget::TfbsFeature { feature_id } => self
                .build_tfbs_expert_view(seq_id, *feature_id)
                .map(FeatureExpertView::Tfbs),
            FeatureExpertTarget::RestrictionSite {
                cut_pos_1based,
                enzyme,
                recognition_start_1based,
                recognition_end_1based,
            } => self
                .build_restriction_site_expert_view(
                    seq_id,
                    *cut_pos_1based,
                    enzyme.as_deref(),
                    *recognition_start_1based,
                    *recognition_end_1based,
                )
                .map(FeatureExpertView::RestrictionSite),
            FeatureExpertTarget::SplicingFeature { feature_id, scope } => self
                .build_splicing_expert_view(seq_id, *feature_id, *scope)
                .map(FeatureExpertView::Splicing),
            FeatureExpertTarget::IsoformArchitecture { panel_id } => self
                .build_isoform_architecture_expert_view(seq_id, panel_id)
                .map(FeatureExpertView::IsoformArchitecture),
            FeatureExpertTarget::ProteinComparison {
                transcript_id_filter,
                protein_feature_filter,
                external_source,
                external_entry_id,
            } => {
                match (external_source, external_entry_id.as_deref()) {
                    (None, None) => self
                        .build_transcript_protein_expert_view(
                            seq_id,
                            transcript_id_filter.as_deref(),
                            protein_feature_filter,
                        )
                        .map(FeatureExpertView::IsoformArchitecture),
                    (
                        Some(gentle_protocol::ProteinExternalOpinionSource::Ensembl),
                        Some(entry_id),
                    ) => self
                        .build_ensembl_entry_protein_expert_view(
                            seq_id,
                            transcript_id_filter.as_deref(),
                            entry_id,
                            protein_feature_filter,
                        )
                        .map(FeatureExpertView::IsoformArchitecture),
                    (Some(source), Some(_)) => Err(EngineError {
                        code: ErrorCode::InvalidInput,
                        message: format!(
                            "Protein comparison external source '{}' is not implemented yet",
                            source.as_str()
                        ),
                    }),
                    _ => Err(EngineError {
                        code: ErrorCode::InvalidInput,
                        message:
                            "Protein comparison external_source and external_entry_id must either both be set or both be omitted"
                                .to_string(),
                    }),
                }
            }
            FeatureExpertTarget::UniprotProjection {
                projection_id,
                protein_feature_filter,
            } => self
                .build_uniprot_projection_expert_view(seq_id, projection_id, protein_feature_filter)
                .map(FeatureExpertView::IsoformArchitecture),
        }
    }

    pub fn render_feature_expert_svg_to_path(
        &self,
        seq_id: &str,
        target: &FeatureExpertTarget,
        path: &str,
    ) -> Result<FeatureExpertView, EngineError> {
        let view = self.inspect_feature_expert(seq_id, target)?;
        let svg = render_feature_expert_svg(&view);
        std::fs::write(path, svg).map_err(|e| EngineError {
            code: ErrorCode::Io,
            message: format!("Could not write feature-expert SVG to '{path}': {e}"),
        })?;
        Ok(view)
    }

    pub fn render_isoform_architecture_svg_to_path(
        &self,
        seq_id: &str,
        panel_id: &str,
        path: &str,
    ) -> Result<IsoformArchitectureExpertView, EngineError> {
        let target = FeatureExpertTarget::IsoformArchitecture {
            panel_id: panel_id.to_string(),
        };
        let view = self.render_feature_expert_svg_to_path(seq_id, &target, path)?;
        match view {
            FeatureExpertView::IsoformArchitecture(isoform) => Ok(isoform),
            _ => Err(EngineError {
                code: ErrorCode::Internal,
                message: "Unexpected expert-view payload while rendering isoform architecture SVG"
                    .to_string(),
            }),
        }
    }
}
