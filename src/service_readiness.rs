//! Service-readiness summary for common GENtle-backed analysis flows.
//!
//! This module combines prepared-reference/helper status with active resource
//! snapshots so callers such as ClawBio can answer "what can this GENtle
//! instance work with right now?" from one deterministic record.

use crate::{
    engine::GentleEngine,
    genomes::{GenomeCatalog, PrepareGenomeActivityStatus},
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
    pub resource_key: String,
    pub display_name: String,
    pub dependency_kind: String,
    pub genome_id: String,
    pub prepared: bool,
    pub lifecycle_status: String,
    pub availability_status: String,
    pub sequence_source_type: String,
    pub annotation_source_type: String,
    pub sequence_source: Option<String>,
    pub annotation_source: Option<String>,
    pub nucleotide_length_bp: Option<usize>,
    pub molecular_mass_da: Option<f64>,
    pub molecular_mass_source: Option<String>,
    pub cache_dir: Option<String>,
    pub current_activity: Option<PrepareGenomeActivityStatus>,
    pub interpretation: Option<Value>,
}

fn now_unix_ms() -> u128 {
    SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .map(|d| d.as_millis())
        .unwrap_or(0)
}

fn summarize_availability_status(lifecycle_status: &str) -> String {
    match lifecycle_status {
        "running" => "preparing".to_string(),
        "ready" => "prepared".to_string(),
        "missing" => "not_prepared".to_string(),
        other => other.to_string(),
    }
}

fn format_activity_brief(activity: &PrepareGenomeActivityStatus) -> String {
    let mut parts = vec![format!("mode {}", activity.prepare_mode)];
    if let Some(phase) = activity.phase.as_ref() {
        parts.push(format!("phase {}", phase));
    }
    if let Some(percent) = activity.percent {
        parts.push(format!("{percent:.1}%"));
    }
    format!(" ({})", parts.join(", "))
}

fn inspect_reference_status(genome_id: &str) -> Result<ServiceDependencyStatus, String> {
    let prepared = GentleEngine::is_reference_genome_prepared(None, genome_id, None)
        .map_err(|e| e.to_string())?;
    let source_plan = GentleEngine::describe_reference_genome_sources(None, genome_id, None)
        .map_err(|e| e.to_string())?;
    let cache_dir = GentleEngine::resolve_reference_genome_cache_dir(None, genome_id, None)
        .map_err(|e| e.to_string())?;
    let current_activity =
        GentleEngine::inspect_reference_genome_prepare_activity(None, genome_id, None)
            .map_err(|e| e.to_string())?;
    let lifecycle_status =
        GenomeCatalog::derive_prepare_lifecycle_status(prepared, current_activity.as_ref());
    Ok(ServiceDependencyStatus {
        resource_key: format!("reference_genome:{genome_id}"),
        display_name: genome_id.to_string(),
        dependency_kind: "reference".to_string(),
        genome_id: genome_id.to_string(),
        prepared,
        lifecycle_status: lifecycle_status.clone(),
        availability_status: summarize_availability_status(&lifecycle_status),
        sequence_source_type: source_plan.sequence_source_type,
        annotation_source_type: source_plan.annotation_source_type,
        sequence_source: Some(source_plan.sequence_source),
        annotation_source: Some(source_plan.annotation_source),
        nucleotide_length_bp: source_plan.nucleotide_length_bp,
        molecular_mass_da: source_plan.molecular_mass_da,
        molecular_mass_source: source_plan.molecular_mass_source,
        cache_dir: Some(cache_dir),
        current_activity,
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
    let current_activity =
        GentleEngine::inspect_helper_genome_prepare_activity(genome_id, None, None)
            .map_err(|e| e.to_string())?;
    let interpretation = GentleEngine::interpret_helper_genome(genome_id, None)
        .map_err(|e| e.to_string())?
        .map(|record| serde_json::to_value(record))
        .transpose()
        .map_err(|e| format!("Could not serialize helper interpretation: {e}"))?;
    let lifecycle_status =
        GenomeCatalog::derive_prepare_lifecycle_status(prepared, current_activity.as_ref());
    Ok(ServiceDependencyStatus {
        resource_key: format!("helper_genome:{genome_id}"),
        display_name: genome_id.to_string(),
        dependency_kind: "helper".to_string(),
        genome_id: genome_id.to_string(),
        prepared,
        lifecycle_status: lifecycle_status.clone(),
        availability_status: summarize_availability_status(&lifecycle_status),
        sequence_source_type: source_plan.sequence_source_type,
        annotation_source_type: source_plan.annotation_source_type,
        sequence_source: Some(source_plan.sequence_source),
        annotation_source: Some(source_plan.annotation_source),
        nucleotide_length_bp: source_plan.nucleotide_length_bp,
        molecular_mass_da: source_plan.molecular_mass_da,
        molecular_mass_source: source_plan.molecular_mass_source,
        cache_dir: Some(cache_dir),
        current_activity,
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
        if reference.lifecycle_status == "running" {
            let activity_suffix = reference
                .current_activity
                .as_ref()
                .map(format_activity_brief)
                .unwrap_or_default();
            if reference.prepared {
                lines.push(format!(
                    "Reference '{}' is prepared locally; a prepare/reindex run is currently active{}.",
                    reference.genome_id, activity_suffix
                ));
            } else {
                lines.push(format!(
                    "Reference '{}' is currently being prepared locally{}.",
                    reference.genome_id, activity_suffix
                ));
            }
        } else if matches!(
            reference.lifecycle_status.as_str(),
            "failed" | "cancelled" | "stale"
        ) {
            let activity = reference.current_activity.as_ref();
            lines.push(format!(
                "Reference '{}' is {} locally; the most recent prepare run ended as '{}'{}{} Retry with `genomes prepare` when ready.",
                reference.genome_id,
                if reference.prepared { "prepared" } else { "not ready" },
                reference.lifecycle_status,
                activity
                    .map(format_activity_brief)
                    .unwrap_or_default(),
                activity
                    .and_then(|status| status.last_error.as_ref())
                    .map(|err| format!(" ({err})"))
                    .unwrap_or_default()
            ));
        } else if reference.prepared {
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
        if helper.lifecycle_status == "running" {
            let activity_suffix = helper
                .current_activity
                .as_ref()
                .map(format_activity_brief)
                .unwrap_or_default();
            if helper.prepared {
                lines.push(format!(
                    "Helper '{}' is prepared locally; a prepare/reindex run is currently active{}.",
                    helper.genome_id, activity_suffix
                ));
            } else {
                lines.push(format!(
                    "Helper '{}' is currently being prepared locally{}.",
                    helper.genome_id, activity_suffix
                ));
            }
        } else if matches!(
            helper.lifecycle_status.as_str(),
            "failed" | "cancelled" | "stale"
        ) {
            let activity = helper.current_activity.as_ref();
            lines.push(format!(
                "Helper '{}' is {} locally; the most recent prepare run ended as '{}'{}{} Retry with `helpers prepare` when ready.",
                helper.genome_id,
                if helper.prepared { "prepared" } else { "not ready" },
                helper.lifecycle_status,
                activity
                    .map(format_activity_brief)
                    .unwrap_or_default(),
                activity
                    .and_then(|status| status.last_error.as_ref())
                    .map(|err| format!(" ({err})"))
                    .unwrap_or_default()
            ));
        } else if helper.prepared {
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
    if resources.attract.runtime_valid {
        lines.push(format!(
            "ATtRACT is active from the {} snapshot ({} motifs); splice-aware RBP evidence can now draw from the normalized motif set, including PWM-backed rows when Matrix_id-mapped PWM blocks are present.",
            resources.attract.active_source, resources.attract.active_item_count
        ));
    } else {
        lines.push(
            "ATtRACT is known to GENtle, but no valid runtime snapshot is active yet; run `resources sync-attract ATtRACT.zip` before requesting splice-aware RBP evidence."
                .to_string(),
        );
    }
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

#[cfg(test)]
mod tests {
    use super::*;
    use crate::resource_status::resource_catalog_status;

    fn fake_activity(status: &str, phase: &str, percent: f64) -> PrepareGenomeActivityStatus {
        PrepareGenomeActivityStatus {
            genome_id: "Human GRCh38 Ensembl 116".to_string(),
            status_path: "/tmp/.prepare_activity.json".to_string(),
            lock_path: Some("/tmp/.prepare_activity.lock".to_string()),
            lifecycle_status: status.to_string(),
            prepare_mode: "prepare_or_reuse".to_string(),
            phase: Some(phase.to_string()),
            item: Some("demo".to_string()),
            bytes_done: 50,
            bytes_total: Some(100),
            percent: Some(percent),
            step_id: None,
            step_label: None,
            started_at_unix_ms: 1,
            updated_at_unix_ms: 2,
            finished_at_unix_ms: None,
            last_error: None,
            owner_pid: Some(12345),
        }
    }

    fn fake_dependency(
        prepared: bool,
        lifecycle_status: &str,
        current_activity: Option<PrepareGenomeActivityStatus>,
    ) -> ServiceDependencyStatus {
        ServiceDependencyStatus {
            resource_key: "reference_genome:Human GRCh38 Ensembl 116".to_string(),
            display_name: "Human GRCh38 Ensembl 116".to_string(),
            dependency_kind: "reference".to_string(),
            genome_id: "Human GRCh38 Ensembl 116".to_string(),
            prepared,
            lifecycle_status: lifecycle_status.to_string(),
            availability_status: summarize_availability_status(lifecycle_status),
            sequence_source_type: "remote_url".to_string(),
            annotation_source_type: "remote_url".to_string(),
            sequence_source: Some("https://example.invalid/sequence.fa.gz".to_string()),
            annotation_source: Some("https://example.invalid/annotation.gtf.gz".to_string()),
            nucleotide_length_bp: None,
            molecular_mass_da: None,
            molecular_mass_source: None,
            cache_dir: Some("/tmp/genomes".to_string()),
            current_activity,
            interpretation: None,
        }
    }

    #[test]
    fn build_summary_lines_reports_active_prepare_for_unprepared_reference() {
        let resources = resource_catalog_status();
        let lines = build_summary_lines(
            &[fake_dependency(
                false,
                "running",
                Some(fake_activity("running", "download_sequence", 42.0)),
            )],
            &[],
            &resources,
        );
        assert!(
            lines
                .iter()
                .any(|line| line.contains("currently being prepared locally")),
            "expected active-prepare summary, got: {lines:?}"
        );
        assert!(
            lines
                .iter()
                .any(|line| line.contains("phase download_sequence")),
            "expected phase mention, got: {lines:?}"
        );
    }

    #[test]
    fn build_summary_lines_reports_failed_prepare_attempt() {
        let resources = resource_catalog_status();
        let mut activity = fake_activity("failed", "index_blast", 80.0);
        activity.last_error = Some("makeblastdb missing".to_string());
        let lines = build_summary_lines(
            &[fake_dependency(false, "failed", Some(activity))],
            &[],
            &resources,
        );
        assert!(
            lines.iter().any(|line| line.contains("ended as 'failed'")),
            "expected failed-prepare summary, got: {lines:?}"
        );
        assert!(
            lines
                .iter()
                .any(|line| line.contains("makeblastdb missing")),
            "expected error mention, got: {lines:?}"
        );
    }
}
