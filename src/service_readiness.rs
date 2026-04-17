//! Service-readiness summary for common GENtle-backed analysis flows.
//!
//! This module combines prepared-reference/helper status with active resource
//! snapshots so callers such as ClawBio can answer "what can this GENtle
//! instance work with right now?" from one deterministic record.

use crate::{
    engine::GentleEngine,
    resource_status::{ResourceCatalogReport, resource_catalog_status},
};
use serde::Serialize;
use serde_json::Value;
use std::time::{SystemTime, UNIX_EPOCH};

pub const SERVICE_READINESS_SCHEMA: &str = "gentle.service_readiness.v1";
pub const DEFAULT_REFERENCE_GENOME_IDS: &[&str] = &["Human GRCh38 Ensembl 116"];
pub const DEFAULT_HELPER_IDS: &[&str] = &["Plasmid pUC19 (online)"];

#[derive(Debug, Clone, Serialize)]
pub struct ServiceReadinessReport {
    pub schema: String,
    pub generated_at_unix_ms: u128,
    pub references: Vec<ServiceDependencyStatus>,
    pub helpers: Vec<ServiceDependencyStatus>,
    pub resources: ResourceCatalogReport,
    pub summary_lines: Vec<String>,
}

#[derive(Debug, Clone, Serialize)]
pub struct ServiceDependencyStatus {
    pub dependency_kind: String,
    pub genome_id: String,
    pub prepared: bool,
    pub sequence_source_type: String,
    pub annotation_source_type: String,
    pub sequence_source: Option<String>,
    pub annotation_source: Option<String>,
    pub nucleotide_length_bp: Option<usize>,
    pub molecular_mass_da: Option<f64>,
    pub molecular_mass_source: Option<String>,
    pub cache_dir: Option<String>,
    pub interpretation: Option<Value>,
}

fn now_unix_ms() -> u128 {
    SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .map(|d| d.as_millis())
        .unwrap_or(0)
}

fn inspect_reference_status(genome_id: &str) -> Result<ServiceDependencyStatus, String> {
    let prepared = GentleEngine::is_reference_genome_prepared(None, genome_id, None)
        .map_err(|e| e.to_string())?;
    let source_plan = GentleEngine::describe_reference_genome_sources(None, genome_id, None)
        .map_err(|e| e.to_string())?;
    let cache_dir = GentleEngine::resolve_reference_genome_cache_dir(None, genome_id, None)
        .map_err(|e| e.to_string())?;
    Ok(ServiceDependencyStatus {
        dependency_kind: "reference".to_string(),
        genome_id: genome_id.to_string(),
        prepared,
        sequence_source_type: source_plan.sequence_source_type,
        annotation_source_type: source_plan.annotation_source_type,
        sequence_source: Some(source_plan.sequence_source),
        annotation_source: Some(source_plan.annotation_source),
        nucleotide_length_bp: source_plan.nucleotide_length_bp,
        molecular_mass_da: source_plan.molecular_mass_da,
        molecular_mass_source: source_plan.molecular_mass_source,
        cache_dir: Some(cache_dir),
        interpretation: None,
    })
}

fn inspect_helper_status(genome_id: &str) -> Result<ServiceDependencyStatus, String> {
    let prepared = GentleEngine::is_helper_genome_prepared(genome_id, None, None)
        .map_err(|e| e.to_string())?;
    let source_plan = GentleEngine::describe_helper_genome_sources(genome_id, None, None)
        .map_err(|e| e.to_string())?;
    let cache_dir = GentleEngine::resolve_helper_genome_cache_dir(genome_id, None, None)
        .map_err(|e| e.to_string())?;
    let interpretation = GentleEngine::interpret_helper_genome(genome_id, None)
        .map_err(|e| e.to_string())?
        .map(|record| serde_json::to_value(record))
        .transpose()
        .map_err(|e| format!("Could not serialize helper interpretation: {e}"))?;
    Ok(ServiceDependencyStatus {
        dependency_kind: "helper".to_string(),
        genome_id: genome_id.to_string(),
        prepared,
        sequence_source_type: source_plan.sequence_source_type,
        annotation_source_type: source_plan.annotation_source_type,
        sequence_source: Some(source_plan.sequence_source),
        annotation_source: Some(source_plan.annotation_source),
        nucleotide_length_bp: source_plan.nucleotide_length_bp,
        molecular_mass_da: source_plan.molecular_mass_da,
        molecular_mass_source: source_plan.molecular_mass_source,
        cache_dir: Some(cache_dir),
        interpretation,
    })
}

fn build_summary_lines(
    references: &[ServiceDependencyStatus],
    helpers: &[ServiceDependencyStatus],
    resources: &ResourceCatalogReport,
) -> Vec<String> {
    let mut lines = vec![];
    for reference in references {
        if reference.prepared {
            lines.push(format!(
                "Reference '{}' is prepared locally and can support genome-backed analysis.",
                reference.genome_id
            ));
        } else {
            lines.push(format!(
                "Reference '{}' is known but not yet prepared locally; genome-backed analysis depending on it will first need `genomes prepare`.",
                reference.genome_id
            ));
        }
    }
    for helper in helpers {
        if helper.prepared {
            lines.push(format!(
                "Helper '{}' is prepared locally and can support helper-backed vector/plasmid workflows.",
                helper.genome_id
            ));
        } else {
            lines.push(format!(
                "Helper '{}' is known but not yet prepared locally; helper-backed vector/plasmid workflows may first need `helpers prepare`.",
                helper.genome_id
            ));
        }
    }
    lines.push(format!(
        "JASPAR is active from the {} snapshot ({} motifs).",
        resources.jaspar.active_source, resources.jaspar.active_item_count
    ));
    lines.push(format!(
        "REBASE is active from the {} snapshot ({} enzymes).",
        resources.rebase.active_source, resources.rebase.active_item_count
    ));
    lines.push(
        "ATtRACT is known to GENtle but not yet integrated, so RNA-binding motif services are not available from it yet."
            .to_string(),
    );
    lines
}

pub fn service_readiness_status() -> Result<ServiceReadinessReport, String> {
    let references = DEFAULT_REFERENCE_GENOME_IDS
        .iter()
        .map(|id| inspect_reference_status(id))
        .collect::<Result<Vec<_>, _>>()?;
    let helpers = DEFAULT_HELPER_IDS
        .iter()
        .map(|id| inspect_helper_status(id))
        .collect::<Result<Vec<_>, _>>()?;
    let resources = resource_catalog_status();
    let summary_lines = build_summary_lines(&references, &helpers, &resources);
    Ok(ServiceReadinessReport {
        schema: SERVICE_READINESS_SCHEMA.to_string(),
        generated_at_unix_ms: now_unix_ms(),
        references,
        helpers,
        resources,
        summary_lines,
    })
}
