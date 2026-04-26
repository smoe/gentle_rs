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
use std::{
    env,
    time::{SystemTime, UNIX_EPOCH},
};

pub const SERVICE_READINESS_SCHEMA: &str = "gentle.service_readiness.v1";
pub const SERVICE_HANDOFF_SCHEMA: &str = "gentle.service_handoff.v1";
pub const TELEGRAM_GUIDE_SCHEMA: &str = "gentle.telegram_guide.v1";
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

#[derive(Debug, Clone, Serialize)]
pub struct ServiceHandoffReport {
    pub schema: String,
    pub generated_at_unix_ms: u128,
    pub scope: String,
    pub service_readiness: ServiceReadinessReport,
    pub readiness: Vec<ServiceHandoffReadinessRow>,
    pub summary_lines: Vec<String>,
    pub suggested_actions: Vec<ServiceHandoffAction>,
    pub running_actions: Vec<ServiceHandoffAction>,
    pub blocked_actions: Vec<ServiceHandoffBlockedAction>,
    pub preferred_demo_actions: Vec<ServiceHandoffAction>,
    pub preferred_artifacts: Vec<ServiceHandoffArtifact>,
    pub environment_hints: Vec<ServiceHandoffEnvironmentHint>,
    pub warnings: Vec<String>,
}

#[derive(Debug, Clone, Serialize)]
pub struct ServiceHandoffReadinessRow {
    pub resource_key: String,
    pub display_name: String,
    pub resource_kind: String,
    pub prepared: bool,
    pub lifecycle_status: String,
    pub status_summary: String,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub cache_dir: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub runtime_path: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub source: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub last_error: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub current_activity: Option<PrepareGenomeActivityStatus>,
}

#[derive(Debug, Clone, Serialize)]
pub struct ServiceHandoffAction {
    pub action_id: String,
    pub label: String,
    pub kind: String,
    pub shell_line: String,
    pub timeout_secs: u64,
    pub rationale: String,
    pub requires_confirmation: bool,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub resource_key: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub lifecycle_status: Option<String>,
    #[serde(skip_serializing_if = "Vec::is_empty")]
    pub expected_artifacts: Vec<String>,
}

#[derive(Debug, Clone, Serialize)]
pub struct ServiceHandoffBlockedAction {
    pub action: ServiceHandoffAction,
    pub blocked_reason: String,
    pub unblock_hint: String,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub download_url: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub local_path_hint: Option<String>,
}

#[derive(Debug, Clone, Serialize)]
pub struct ServiceHandoffArtifact {
    pub artifact_id: String,
    pub path: String,
    pub caption: String,
    pub recommended_use: String,
    pub presentation_rank: usize,
    pub is_best_first_artifact: bool,
}

#[derive(Debug, Clone, Serialize)]
pub struct ServiceHandoffEnvironmentHint {
    pub name: String,
    pub is_set: bool,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub current_value: Option<String>,
    pub purpose: String,
    pub recommended_when: String,
}

#[derive(Debug, Clone, Serialize)]
pub struct TelegramGuideReport {
    pub schema: String,
    pub generated_at_unix_ms: u128,
    pub channel: String,
    pub section: String,
    pub gene: Option<String>,
    pub gene_supplied: bool,
    pub summary_lines: Vec<String>,
    pub readiness_summary_lines: Vec<String>,
    pub menu_sections: Vec<TelegramGuideSection>,
    pub suggested_actions: Vec<ServiceHandoffAction>,
    pub blocked_actions: Vec<ServiceHandoffBlockedAction>,
    pub warnings: Vec<String>,
}

#[derive(Debug, Clone, Serialize)]
pub struct TelegramGuideSection {
    pub section_id: String,
    pub title: String,
    pub summary: String,
    pub example_prompts: Vec<String>,
    pub default_genes: Vec<String>,
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

fn shell_quote_arg(value: &str) -> String {
    format!("\"{}\"", value.replace('\\', "\\\\").replace('"', "\\\""))
}

fn action_slug(value: &str) -> String {
    let mut slug = String::new();
    let mut previous_was_separator = false;
    for ch in value.chars() {
        if ch.is_ascii_alphanumeric() {
            slug.push(ch.to_ascii_lowercase());
            previous_was_separator = false;
        } else if !previous_was_separator && !slug.is_empty() {
            slug.push('_');
            previous_was_separator = true;
        }
    }
    while slug.ends_with('_') {
        slug.pop();
    }
    if slug.is_empty() {
        "action".to_string()
    } else {
        slug
    }
}

fn handoff_action(
    label: impl Into<String>,
    kind: impl Into<String>,
    shell_line: impl Into<String>,
    timeout_secs: u64,
    rationale: impl Into<String>,
    requires_confirmation: bool,
    resource_key: Option<String>,
    lifecycle_status: Option<String>,
    expected_artifacts: Vec<String>,
) -> ServiceHandoffAction {
    let label = label.into();
    ServiceHandoffAction {
        action_id: action_slug(&label),
        label,
        kind: kind.into(),
        shell_line: shell_line.into(),
        timeout_secs,
        rationale: rationale.into(),
        requires_confirmation,
        resource_key,
        lifecycle_status,
        expected_artifacts,
    }
}

fn reference_prepare_action(
    reference: &ServiceDependencyStatus,
    retry: bool,
) -> ServiceHandoffAction {
    let shell_line = format!(
        "genomes prepare {} --timeout-secs 7200",
        shell_quote_arg(&reference.genome_id)
    );
    handoff_action(
        if retry {
            format!("Retry prepare for {}", reference.genome_id)
        } else {
            format!("Prepare {}", reference.genome_id)
        },
        "prepare_reference",
        shell_line,
        7500,
        if retry {
            format!(
                "Reference '{}' last ended as {} and is safe to retry when genome-backed analysis is needed.",
                reference.genome_id, reference.lifecycle_status
            )
        } else {
            format!(
                "Reference '{}' is known but not prepared locally; genome-backed analysis needs a prepared cache.",
                reference.genome_id
            )
        },
        true,
        Some(reference.resource_key.clone()),
        Some(reference.lifecycle_status.clone()),
        vec![],
    )
}

fn helper_prepare_action(helper: &ServiceDependencyStatus, retry: bool) -> ServiceHandoffAction {
    let shell_line = format!(
        "helpers prepare {} --timeout-secs 1800",
        shell_quote_arg(&helper.genome_id)
    );
    handoff_action(
        if retry {
            format!("Retry prepare for {}", helper.genome_id)
        } else {
            format!("Prepare {}", helper.genome_id)
        },
        "prepare_helper",
        shell_line,
        2100,
        if retry {
            format!(
                "Helper '{}' last ended as {} and is safe to retry when helper-backed workflows are needed.",
                helper.genome_id, helper.lifecycle_status
            )
        } else {
            format!(
                "Helper '{}' is known but not prepared locally; helper-backed vector/plasmid workflows need a prepared cache.",
                helper.genome_id
            )
        },
        true,
        Some(helper.resource_key.clone()),
        Some(helper.lifecycle_status.clone()),
        vec![],
    )
}

fn status_refresh_action() -> ServiceHandoffAction {
    handoff_action(
        "Re-check services status",
        "refresh_status",
        "services status",
        180,
        "A shared prepare action is already running, so refresh the combined readiness view instead of starting duplicate long-running work.",
        false,
        None,
        Some("running".to_string()),
        vec![],
    )
}

fn dependency_readiness_row(dependency: &ServiceDependencyStatus) -> ServiceHandoffReadinessRow {
    let activity = dependency.current_activity.as_ref();
    ServiceHandoffReadinessRow {
        resource_key: dependency.resource_key.clone(),
        display_name: dependency.display_name.clone(),
        resource_kind: dependency.dependency_kind.clone(),
        prepared: dependency.prepared,
        lifecycle_status: dependency.lifecycle_status.clone(),
        status_summary: match dependency.lifecycle_status.as_str() {
            "running" => format!(
                "{} '{}' is currently preparing{}.",
                dependency.dependency_kind,
                dependency.display_name,
                activity.map(format_activity_brief).unwrap_or_default()
            ),
            "ready" => format!(
                "{} '{}' is ready for reuse.",
                dependency.dependency_kind, dependency.display_name
            ),
            "failed" | "cancelled" | "stale" => format!(
                "{} '{}' is not ready; latest activity is {}{}.",
                dependency.dependency_kind,
                dependency.display_name,
                dependency.lifecycle_status,
                activity
                    .and_then(|status| status.last_error.as_ref())
                    .map(|err| format!(" ({err})"))
                    .unwrap_or_default()
            ),
            _ => format!(
                "{} '{}' is known but not prepared.",
                dependency.dependency_kind, dependency.display_name
            ),
        },
        cache_dir: dependency.cache_dir.clone(),
        runtime_path: None,
        source: dependency.sequence_source.clone(),
        last_error: activity.and_then(|status| status.last_error.clone()),
        current_activity: dependency.current_activity.clone(),
    }
}

fn resource_readiness_rows(resources: &ResourceCatalogReport) -> Vec<ServiceHandoffReadinessRow> {
    let mut rows = vec![];
    rows.push(ServiceHandoffReadinessRow {
        resource_key: "resource:jaspar".to_string(),
        display_name: "JASPAR".to_string(),
        resource_kind: "resource".to_string(),
        prepared: resources.jaspar.active_item_count > 0,
        lifecycle_status: if resources.jaspar.active_item_count > 0 {
            "ready"
        } else {
            "missing"
        }
        .to_string(),
        status_summary: format!(
            "JASPAR is active from the {} snapshot with {} motifs.",
            resources.jaspar.active_source, resources.jaspar.active_item_count
        ),
        cache_dir: None,
        runtime_path: Some(resources.jaspar.runtime_path.clone()),
        source: Some(resources.jaspar.active_source.clone()),
        last_error: resources.jaspar.runtime_error.clone(),
        current_activity: None,
    });
    rows.push(ServiceHandoffReadinessRow {
        resource_key: "resource:rebase".to_string(),
        display_name: "REBASE".to_string(),
        resource_kind: "resource".to_string(),
        prepared: resources.rebase.active_item_count > 0,
        lifecycle_status: if resources.rebase.active_item_count > 0 {
            "ready"
        } else {
            "missing"
        }
        .to_string(),
        status_summary: format!(
            "REBASE is active from the {} snapshot with {} enzymes.",
            resources.rebase.active_source, resources.rebase.active_item_count
        ),
        cache_dir: None,
        runtime_path: Some(resources.rebase.runtime_path.clone()),
        source: Some(resources.rebase.active_source.clone()),
        last_error: resources.rebase.runtime_error.clone(),
        current_activity: None,
    });
    rows.push(ServiceHandoffReadinessRow {
        resource_key: "resource:attract".to_string(),
        display_name: resources.attract.display_name.clone(),
        resource_kind: "resource".to_string(),
        prepared: resources.attract.runtime_valid,
        lifecycle_status: if resources.attract.runtime_valid {
            "ready".to_string()
        } else if resources.attract.runtime_exists {
            "failed".to_string()
        } else {
            "missing".to_string()
        },
        status_summary: if resources.attract.runtime_valid {
            format!(
                "ATtRACT is active from the {} snapshot with {} motifs.",
                resources.attract.active_source, resources.attract.active_item_count
            )
        } else {
            "ATtRACT is known, but no valid runtime snapshot is active yet.".to_string()
        },
        cache_dir: None,
        runtime_path: Some(resources.attract.runtime_path.clone()),
        source: Some(resources.attract.active_source.clone()),
        last_error: resources.attract.runtime_error.clone(),
        current_activity: None,
    });
    rows
}

fn environment_hint(
    name: &str,
    purpose: &str,
    recommended_when: &str,
) -> ServiceHandoffEnvironmentHint {
    let current_value = env::var(name).ok().filter(|value| !value.trim().is_empty());
    ServiceHandoffEnvironmentHint {
        name: name.to_string(),
        is_set: current_value.is_some(),
        current_value,
        purpose: purpose.to_string(),
        recommended_when: recommended_when.to_string(),
    }
}

fn build_environment_hints() -> Vec<ServiceHandoffEnvironmentHint> {
    vec![
        environment_hint(
            "GENTLE_CLI_CMD",
            "Command used by wrappers such as ClawBio to invoke the intended gentle_cli binary or container route.",
            "Set this when the chat host should use a specific checkout/container instead of whatever gentle_cli appears first on PATH.",
        ),
        environment_hint(
            "GENTLE_REPO_ROOT",
            "Local editable GENtle checkout used by gentle_local_checkout_cli.sh.",
            "Set this for ClawBio deployments that pull GENtle from GitHub and execute the local checkout directly.",
        ),
        environment_hint(
            "GENTLE_REFERENCE_CACHE_DIR",
            "Shared prepared reference-genome cache root.",
            "Set this to keep Ensembl/reference downloads outside transient worktrees and reusable across chat sessions.",
        ),
        environment_hint(
            "GENTLE_HELPER_CACHE_DIR",
            "Shared helper-genome/vector cache root.",
            "Set this to keep helper plasmid/vector preparations reusable across chat sessions.",
        ),
        environment_hint(
            "GENTLE_CUTRUN_CACHE_DIR",
            "Shared CUT&RUN dataset cache root.",
            "Set this when dataset-backed regulatory support routes should reuse prepared evidence sets.",
        ),
    ]
}

fn build_preferred_demo_actions(reference_ready: bool) -> Vec<ServiceHandoffAction> {
    let mut actions = vec![
        handoff_action(
            "Render Gibson protocol cartoon",
            "demo_graphic",
            "protocol-cartoon render-svg gibson.two_fragment artifacts/gibson.two_fragment.protocol.svg",
            180,
            "This is the lowest-friction graphical GENtle demo because it does not require prepared Ensembl assets.",
            false,
            None,
            None,
            vec!["artifacts/gibson.two_fragment.protocol.svg".to_string()],
        ),
        handoff_action(
            "Resolve TF query examples",
            "demo_tf_query",
            "resources resolve-tf-query stemness OCT4 \"KLF family\" --output artifacts/tf_query_resolution.stemness_oct4_klf.json",
            180,
            "This demonstrates user-extensible TF group resolution against the active local motif registry.",
            false,
            Some("resource:jaspar".to_string()),
            Some("ready".to_string()),
            vec!["artifacts/tf_query_resolution.stemness_oct4_klf.json".to_string()],
        ),
    ];
    if reference_ready {
        actions.push(handoff_action(
            "Render TERT/TP73 promoter TFBS SVG",
            "demo_promoter_tfbs_svg",
            "genomes promoter-tfbs-svg \"Human GRCh38 Ensembl 116\" --gene TERT --gene TP73 --motif stemness --motif SP1 --upstream-bp 1000 --downstream-bp 200 --score-kind llr_background_tail_log10 artifacts/grch38_tert_tp73_promoters.stemness_sp1.svg",
            1800,
            "This demonstrates the multi-gene promoter-design display once the human Ensembl reference is prepared.",
            false,
            Some("reference_genome:Human GRCh38 Ensembl 116".to_string()),
            Some("ready".to_string()),
            vec!["artifacts/grch38_tert_tp73_promoters.stemness_sp1.svg".to_string()],
        ));
    }
    actions
}

fn normalize_telegram_section(section: Option<&str>) -> Result<String, String> {
    let normalized = section
        .map(str::trim)
        .filter(|value| !value.is_empty())
        .unwrap_or("overview")
        .to_ascii_lowercase()
        .replace('_', "-");
    match normalized.as_str() {
        "overview" | "readiness" | "gene-context" | "tfbs" | "inline-dna" | "cloning"
        | "isoforms" | "follow-up" => Ok(normalized),
        other => Err(format!(
            "Unknown services guide section '{other}' (expected overview, readiness, gene-context, tfbs, inline-dna, cloning, isoforms or follow-up)"
        )),
    }
}

fn normalize_gene_symbol(gene: Option<&str>) -> Option<String> {
    gene.map(str::trim)
        .filter(|value| !value.is_empty())
        .map(|value| value.to_ascii_uppercase())
}

fn guide_section_action(
    label: impl Into<String>,
    section: &str,
    gene: Option<&str>,
    rationale: impl Into<String>,
) -> ServiceHandoffAction {
    let mut shell_line = format!("services guide --channel telegram --section {section}");
    if let Some(gene) = gene {
        shell_line.push_str(" --gene ");
        shell_line.push_str(&shell_quote_arg(gene));
    }
    handoff_action(
        label,
        "guide_section",
        shell_line,
        180,
        rationale,
        false,
        None,
        None,
        vec![],
    )
}

fn default_guide_sections() -> Vec<TelegramGuideSection> {
    vec![
        TelegramGuideSection {
            section_id: "readiness".to_string(),
            title: "Data/readiness".to_string(),
            summary: "Check GENtle version, prepared references, helper vectors, and resource snapshots.".to_string(),
            example_prompts: vec![
                "What GENtle data are ready?".to_string(),
                "Can you prepare Ensembl-backed human data?".to_string(),
            ],
            default_genes: vec![],
        },
        TelegramGuideSection {
            section_id: "gene-context".to_string(),
            title: "Gene context".to_string(),
            summary: "Extract a gene, promoter, or anchored locus from prepared Ensembl-backed references.".to_string(),
            example_prompts: vec![
                "Show me the TP73 gene context.".to_string(),
                "Extract the TERT promoter.".to_string(),
            ],
            default_genes: vec!["TP73".to_string()],
        },
        TelegramGuideSection {
            section_id: "tfbs".to_string(),
            title: "Promoter and TFBS".to_string(),
            summary: "Summarize or render TFBS/PSSM score tracks for promoter windows and factor groups.".to_string(),
            example_prompts: vec![
                "Show stemness and SP1 TFBS upstream of TERT.".to_string(),
                "Compare promoter TFBS score tracks for TP73.".to_string(),
            ],
            default_genes: vec!["TERT".to_string(), "TP73".to_string()],
        },
        TelegramGuideSection {
            section_id: "inline-dna".to_string(),
            title: "Pasted DNA inspection".to_string(),
            summary: "Scan pasted DNA without creating project state when the request is read-only.".to_string(),
            example_prompts: vec![
                "Scan this DNA for EcoRI and SmaI.".to_string(),
                "Find SP1-like TFBS hits in this sequence.".to_string(),
            ],
            default_genes: vec![],
        },
        TelegramGuideSection {
            section_id: "cloning".to_string(),
            title: "Cloning and vectors".to_string(),
            summary: "Render cloning cartoons, inspect helper vectors, and prepare vector/helper caches.".to_string(),
            example_prompts: vec![
                "Show me a Gibson assembly cartoon.".to_string(),
                "Is pUC19 prepared for helper-vector workflows?".to_string(),
            ],
            default_genes: vec![],
        },
        TelegramGuideSection {
            section_id: "isoforms".to_string(),
            title: "Isoforms and protein gels".to_string(),
            summary: "Compare transcript-native isoforms with gel, 2D-gel, or digest-style figures.".to_string(),
            example_prompts: vec![
                "Show TP73 isoforms as a protein gel.".to_string(),
                "Can you compare TP53 splicing or isoform architecture?".to_string(),
            ],
            default_genes: vec!["TP73".to_string(), "TP53".to_string()],
        },
        TelegramGuideSection {
            section_id: "follow-up".to_string(),
            title: "Experimental follow-up".to_string(),
            summary: "Turn a SNP, expression, or splicing observation into a reproducible validation-planning handoff.".to_string(),
            example_prompts: vec![
                "Plan a promoter-reporter follow-up for this SNP.".to_string(),
                "What should we validate for this differentially expressed gene?".to_string(),
            ],
            default_genes: vec!["VKORC1".to_string(), "TERT".to_string()],
        },
    ]
}

fn gene_or_default<'a>(gene: Option<&'a str>, default_gene: &'a str) -> String {
    gene.unwrap_or(default_gene).to_string()
}

fn promoter_tfbs_genes(gene: Option<&str>) -> Vec<String> {
    match gene {
        Some(gene) => vec![gene.to_string()],
        None => vec!["TERT".to_string(), "TP73".to_string()],
    }
}

fn gene_shell_args(genes: &[String]) -> String {
    genes
        .iter()
        .map(|gene| format!("--gene {}", shell_quote_arg(gene)))
        .collect::<Vec<_>>()
        .join(" ")
}

fn first_expected_artifact(path: impl Into<String>) -> Vec<String> {
    vec![path.into()]
}

fn guide_actions_for_section(
    section: &str,
    gene: Option<&str>,
    handoff: &ServiceHandoffReport,
) -> Vec<ServiceHandoffAction> {
    match section {
        "readiness" => {
            let mut actions = vec![handoff_action(
                "Refresh GENtle readiness",
                "refresh_status",
                "services status",
                180,
                "Refresh the combined reference/helper/resource readiness view.",
                false,
                None,
                None,
                vec![],
            )];
            actions.extend(handoff.running_actions.iter().take(2).cloned());
            actions.extend(handoff.suggested_actions.iter().take(2).cloned());
            actions.truncate(4);
            actions
        }
        "gene-context" => {
            let gene = gene_or_default(gene, "TP73");
            let gene_slug = action_slug(&gene);
            vec![
                handoff_action(
                    format!("Extract {gene} gene context"),
                    "gene_context",
                    format!(
                        "genomes extract-gene \"Human GRCh38 Ensembl 116\" {} --occurrence 1 --output-id grch38_{}",
                        shell_quote_arg(&gene),
                        gene_slug
                    ),
                    1200,
                    format!("Extract {gene} from the prepared human Ensembl-backed reference."),
                    true,
                    Some("reference_genome:Human GRCh38 Ensembl 116".to_string()),
                    None,
                    vec![],
                ),
                handoff_action(
                    format!("Extract {gene} promoter"),
                    "gene_context",
                    format!(
                        "genomes extract-promoter \"Human GRCh38 Ensembl 116\" {} --output-id grch38_{}_promoter --upstream-bp 1000 --downstream-bp 200",
                        shell_quote_arg(&gene),
                        gene_slug
                    ),
                    1200,
                    format!("Extract a default 1000 bp upstream / 200 bp downstream promoter window for {gene}."),
                    true,
                    Some("reference_genome:Human GRCh38 Ensembl 116".to_string()),
                    None,
                    vec![],
                ),
                guide_section_action(
                    format!("Show promoter/TFBS options for {gene}"),
                    "tfbs",
                    Some(&gene),
                    "Jump to promoter and transcription-factor analyses for the same gene.",
                ),
            ]
        }
        "tfbs" => {
            let genes = promoter_tfbs_genes(gene);
            let genes_label = genes.join("/");
            let genes_slug = genes
                .iter()
                .map(|gene| action_slug(gene))
                .collect::<Vec<_>>()
                .join("_");
            let gene_args = gene_shell_args(&genes);
            let summary_path =
                format!("artifacts/grch38_{genes_slug}_promoters.stemness_sp1.summary.json");
            let svg_path = format!("artifacts/grch38_{genes_slug}_promoters.stemness_sp1.svg");
            vec![
                handoff_action(
                    format!("Summarize {genes_label} promoter TFBS"),
                    "promoter_tfbs_summary",
                    format!(
                        "genomes promoter-tfbs-summary \"Human GRCh38 Ensembl 116\" {gene_args} --motif stemness --motif SP1 --upstream-bp 1000 --downstream-bp 200 --score-kind llr_background_tail_log10 --path {summary_path}",
                    ),
                    1800,
                    format!("Summarize stemness/Yamanaka and SP1 motif evidence across the {genes_label} promoter window."),
                    true,
                    Some("reference_genome:Human GRCh38 Ensembl 116".to_string()),
                    None,
                    first_expected_artifact(summary_path),
                ),
                handoff_action(
                    format!("Render {genes_label} promoter TFBS figure"),
                    "promoter_tfbs_svg",
                    format!(
                        "genomes promoter-tfbs-svg \"Human GRCh38 Ensembl 116\" {gene_args} --motif stemness --motif SP1 --upstream-bp 1000 --downstream-bp 200 --score-kind llr_background_tail_log10 {svg_path}",
                    ),
                    1800,
                    "Render a Telegram-friendly promoter/TFBS SVG that the wrapper can rasterize into PNG.",
                    true,
                    Some("reference_genome:Human GRCh38 Ensembl 116".to_string()),
                    None,
                    first_expected_artifact(svg_path),
                ),
                handoff_action(
                    "Resolve stemness and KLF TF names",
                    "tf_query",
                    "resources resolve-tf-query stemness OCT4 \"KLF family\" --output artifacts/tf_query_resolution.stemness_oct4_klf.json",
                    180,
                    "Show how GENtle expands user-supplied TF groups or family names before scanning.",
                    false,
                    Some("resource:jaspar".to_string()),
                    Some("ready".to_string()),
                    first_expected_artifact("artifacts/tf_query_resolution.stemness_oct4_klf.json"),
                ),
            ]
        }
        "inline-dna" => vec![
            handoff_action(
                "Run stateless pasted-DNA demo",
                "inline_dna_demo",
                "workflow @docs/examples/workflows/inline_sequence_inspection_stateless_offline.json",
                300,
                "Demonstrate read-only restriction and TFBS inspection without creating project state.",
                false,
                None,
                None,
                vec![],
            ),
            handoff_action(
                "Check restriction and motif resources",
                "refresh_status",
                "resources status",
                180,
                "Confirm REBASE and JASPAR are active before scanning pasted DNA.",
                false,
                None,
                None,
                vec![],
            ),
        ],
        "cloning" => {
            let mut actions = handoff
                .preferred_demo_actions
                .iter()
                .filter(|action| action.kind == "demo_graphic")
                .take(1)
                .cloned()
                .collect::<Vec<_>>();
            if actions.is_empty() {
                actions.push(handoff_action(
                    "Render Gibson protocol cartoon",
                    "demo_graphic",
                    "protocol-cartoon render-svg gibson.two_fragment artifacts/gibson.two_fragment.protocol.svg",
                    180,
                    "Fast graphical cloning demo that does not require Ensembl preparation.",
                    false,
                    None,
                    None,
                    first_expected_artifact("artifacts/gibson.two_fragment.protocol.svg"),
                ));
            }
            actions.push(handoff_action(
                "Check pUC19 helper status",
                "helper_status",
                "helpers status \"Plasmid pUC19 (online)\"",
                180,
                "Check whether the canonical pUC19 helper vector is prepared for vector-backed workflows.",
                false,
                Some("helper_genome:Plasmid pUC19 (online)".to_string()),
                None,
                vec![],
            ));
            actions.truncate(3);
            actions
        }
        "isoforms" => {
            let gene = gene_or_default(gene, "TP73");
            let mut actions = vec![
                handoff_action(
                    "Show TP73 protein gel demo",
                    "isoform_protein_gel",
                    "workflow @docs/examples/workflows/tp73_isoform_protein_gel_offline.json",
                    300,
                    "Render the offline curated TP73 isoform molecular-weight gel.",
                    false,
                    None,
                    None,
                    first_expected_artifact("exports/tp73_isoform_protein_gel.svg"),
                ),
                handoff_action(
                    "Show TP73 2D protein gel demo",
                    "isoform_protein_2d_gel",
                    "workflow @docs/examples/workflows/tp73_isoform_protein_2d_gel_offline.json",
                    300,
                    "Render the offline curated TP73 pI-vs-molecular-weight gel.",
                    false,
                    None,
                    None,
                    first_expected_artifact("exports/tp73_isoform_protein_2d_gel.svg"),
                ),
            ];
            if gene != "TP73" {
                actions.push(guide_section_action(
                    format!("Start with {gene} gene context"),
                    "gene-context",
                    Some(&gene),
                    "Use the selected gene as context before choosing an isoform-specific route.",
                ));
            }
            actions
        }
        "follow-up" => {
            let gene = gene_or_default(gene, "TERT");
            vec![
                guide_section_action(
                    format!("Inspect promoter/TFBS for {gene} first"),
                    "tfbs",
                    Some(&gene),
                    "Jump to promoter evidence before proposing wet-lab validation.",
                ),
                handoff_action(
                    "Run VKORC1 promoter-reporter planning demo",
                    "experimental_followup_demo",
                    "workflow @docs/examples/workflows/vkorc1_rs9923231_promoter_luciferase_assay_planning.json",
                    1800,
                    "Show a deterministic SNP-to-promoter-reporter validation planning example.",
                    false,
                    None,
                    None,
                    vec![
                        "docs/tutorial/reproducibility/vkorc1_rs9923231_promoter_reporter/vkorc1_rs9923231_promoter_context.svg".to_string(),
                        "docs/tutorial/reproducibility/vkorc1_rs9923231_promoter_reporter/vkorc1_rs9923231_reporter_reference.svg".to_string(),
                        "docs/tutorial/reproducibility/vkorc1_rs9923231_promoter_reporter/vkorc1_rs9923231_reporter_alternate.svg".to_string(),
                    ],
                ),
            ]
        }
        _ => {
            let mut actions = vec![guide_section_action(
                "Continue with default genes",
                "overview",
                None,
                "Show the default Telegram guide without gene personalization.",
            )];
            let gene_for_links = gene;
            for (section, label) in [
                ("readiness", "Data/readiness"),
                ("gene-context", "Gene context"),
                ("tfbs", "Promoter and TFBS"),
                ("inline-dna", "Pasted DNA inspection"),
                ("cloning", "Cloning and vectors"),
                ("isoforms", "Isoforms and protein gels"),
                ("follow-up", "Experimental follow-up"),
            ] {
                actions.push(guide_section_action(
                    label,
                    section,
                    gene_for_links,
                    format!("Open the {label} guide section."),
                ));
            }
            actions
        }
    }
}

fn build_telegram_guide_from_handoff(
    handoff: ServiceHandoffReport,
    channel: &str,
    section: &str,
    gene: Option<String>,
) -> TelegramGuideReport {
    let gene_supplied = gene.is_some();
    let gene_ref = gene.as_deref();
    let mut summary_lines = vec![
        "GENtle can guide reproducible sequence, promoter, cloning, isoform, and follow-up work from Telegram.".to_string(),
        "If you have a gene of interest, tell me its symbol. Otherwise I will use defaults for each section.".to_string(),
    ];
    if let Some(gene) = gene_ref {
        summary_lines.push(format!("Personalized guide context: using gene {gene}."));
    } else {
        summary_lines.push("Default guide context: TERT/TP73 for promoter-TFBS, TP73/TP53 for isoforms, TP73 for gene context.".to_string());
    }
    if section != "overview" {
        summary_lines.push(format!("Opened GENtle guide section: {section}."));
    }

    let mut warnings = handoff.warnings.clone();
    if channel != "telegram" {
        warnings.push(format!(
            "Guide channel '{channel}' was accepted, but the current compact presentation is optimized for Telegram."
        ));
    }
    let readiness_summary_lines = handoff
        .service_readiness
        .summary_lines
        .iter()
        .take(5)
        .cloned()
        .collect();
    let suggested_actions = guide_actions_for_section(section, gene_ref, &handoff);

    TelegramGuideReport {
        schema: TELEGRAM_GUIDE_SCHEMA.to_string(),
        generated_at_unix_ms: now_unix_ms(),
        channel: channel.to_string(),
        section: section.to_string(),
        gene,
        gene_supplied,
        summary_lines,
        readiness_summary_lines,
        menu_sections: default_guide_sections(),
        suggested_actions,
        blocked_actions: handoff.blocked_actions,
        warnings,
    }
}

pub fn telegram_guide_report(
    channel: Option<&str>,
    section: Option<&str>,
    gene: Option<&str>,
) -> Result<TelegramGuideReport, String> {
    let channel = channel
        .map(str::trim)
        .filter(|value| !value.is_empty())
        .unwrap_or("telegram")
        .to_ascii_lowercase();
    if channel != "telegram" {
        return Err(format!(
            "Unknown services guide channel '{channel}' (expected telegram)"
        ));
    }
    let section = normalize_telegram_section(section)?;
    let gene = normalize_gene_symbol(gene);
    let handoff = service_handoff_report(Some("telegram"), None)?;
    Ok(build_telegram_guide_from_handoff(
        handoff, &channel, &section, gene,
    ))
}

pub fn service_handoff_report(
    scope: Option<&str>,
    export_path: Option<String>,
) -> Result<ServiceHandoffReport, String> {
    let scope = scope
        .map(str::trim)
        .filter(|value| !value.is_empty())
        .unwrap_or("clawbio")
        .to_string();
    let service_readiness = service_readiness_status()?;
    let mut readiness = vec![];
    let mut suggested_actions = vec![];
    let mut running_actions = vec![];
    let mut blocked_actions = vec![];
    let mut warnings = vec![];
    let mut has_running_prepare = false;

    for reference in &service_readiness.references {
        readiness.push(dependency_readiness_row(reference));
        match reference.lifecycle_status.as_str() {
            "missing" | "not_prepared" => {
                suggested_actions.push(reference_prepare_action(reference, false));
            }
            "failed" | "cancelled" | "stale" => {
                suggested_actions.push(reference_prepare_action(reference, true));
            }
            "running" => {
                has_running_prepare = true;
                running_actions.push(handoff_action(
                    format!("Re-check {} status", reference.genome_id),
                    "refresh_status",
                    format!("genomes status {}", shell_quote_arg(&reference.genome_id)),
                    180,
                    format!(
                        "Reference '{}' is already being prepared; refresh status rather than launching another prepare.",
                        reference.genome_id
                    ),
                    false,
                    Some(reference.resource_key.clone()),
                    Some(reference.lifecycle_status.clone()),
                    vec![],
                ));
            }
            _ => {}
        }
    }
    for helper in &service_readiness.helpers {
        readiness.push(dependency_readiness_row(helper));
        match helper.lifecycle_status.as_str() {
            "missing" | "not_prepared" => {
                suggested_actions.push(helper_prepare_action(helper, false));
            }
            "failed" | "cancelled" | "stale" => {
                suggested_actions.push(helper_prepare_action(helper, true));
            }
            "running" => {
                has_running_prepare = true;
                running_actions.push(handoff_action(
                    format!("Re-check {} status", helper.genome_id),
                    "refresh_status",
                    format!("helpers status {}", shell_quote_arg(&helper.genome_id)),
                    180,
                    format!(
                        "Helper '{}' is already being prepared; refresh status rather than launching another prepare.",
                        helper.genome_id
                    ),
                    false,
                    Some(helper.resource_key.clone()),
                    Some(helper.lifecycle_status.clone()),
                    vec![],
                ));
            }
            _ => {}
        }
    }
    readiness.extend(resource_readiness_rows(&service_readiness.resources));

    if has_running_prepare {
        suggested_actions.push(status_refresh_action());
    }

    if !service_readiness.resources.attract.runtime_valid {
        let action = handoff_action(
            "Sync ATtRACT runtime snapshot",
            "sync_resource",
            "resources sync-attract /path/to/ATtRACT.zip",
            900,
            "ATtRACT is known to GENtle, but the published ZIP must be available locally before GENtle can normalize it into a runtime snapshot.",
            true,
            Some("resource:attract".to_string()),
            Some(if service_readiness.resources.attract.runtime_exists {
                "failed".to_string()
            } else {
                "missing".to_string()
            }),
            vec!["data/resources/attract.motifs.json".to_string()],
        );
        blocked_actions.push(ServiceHandoffBlockedAction {
            action,
            blocked_reason: "requires_local_archive_path".to_string(),
            unblock_hint: "Download ATtRACT.zip first, then replace /path/to/ATtRACT.zip with that local file path.".to_string(),
            download_url: service_readiness.resources.attract.download_url.clone(),
            local_path_hint: Some("/path/to/ATtRACT.zip".to_string()),
        });
        warnings.push(
            "ATtRACT-backed splice/RBP evidence is not ready until the local ZIP has been normalized with `resources sync-attract`."
                .to_string(),
        );
    }

    if service_readiness
        .references
        .iter()
        .any(|reference| reference.lifecycle_status != "ready")
    {
        warnings.push(
            "Genome-backed promoter and locus extraction demos should first prepare the human Ensembl reference."
                .to_string(),
        );
    }

    let reference_ready = service_readiness.references.iter().any(|reference| {
        reference.resource_key == "reference_genome:Human GRCh38 Ensembl 116"
            && reference.lifecycle_status == "ready"
    });
    let preferred_demo_actions = build_preferred_demo_actions(reference_ready);
    let mut preferred_artifacts = vec![];
    if let Some(path) = export_path {
        preferred_artifacts.push(ServiceHandoffArtifact {
            artifact_id: "service_handoff_json".to_string(),
            path,
            caption: "Machine-readable GENtle service handoff report for chat gateways."
                .to_string(),
            recommended_use: "installation_doctor_payload".to_string(),
            presentation_rank: 0,
            is_best_first_artifact: true,
        });
    }

    let mut summary_lines = vec![format!(
        "Built GENtle service handoff for scope '{scope}' with {} readiness rows, {} suggested action(s), {} running action(s), and {} blocked action(s).",
        readiness.len(),
        suggested_actions.len(),
        running_actions.len(),
        blocked_actions.len()
    )];
    summary_lines.extend(service_readiness.summary_lines.clone());
    if !blocked_actions.is_empty() {
        summary_lines.push(
            "Some setup actions are intentionally blocked until the required local input file or path is supplied."
                .to_string(),
        );
    }

    Ok(ServiceHandoffReport {
        schema: SERVICE_HANDOFF_SCHEMA.to_string(),
        generated_at_unix_ms: now_unix_ms(),
        scope,
        service_readiness,
        readiness,
        summary_lines,
        suggested_actions,
        running_actions,
        blocked_actions,
        preferred_demo_actions,
        preferred_artifacts,
        environment_hints: build_environment_hints(),
        warnings,
    })
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

    fn fake_handoff_report(
        reference: ServiceDependencyStatus,
        suggested_actions: Vec<ServiceHandoffAction>,
        running_actions: Vec<ServiceHandoffAction>,
    ) -> ServiceHandoffReport {
        let resources = resource_catalog_status();
        let summary_lines = build_summary_lines(std::slice::from_ref(&reference), &[], &resources);
        ServiceHandoffReport {
            schema: SERVICE_HANDOFF_SCHEMA.to_string(),
            generated_at_unix_ms: 42,
            scope: "telegram".to_string(),
            service_readiness: ServiceReadinessReport {
                schema: SERVICE_READINESS_SCHEMA.to_string(),
                generated_at_unix_ms: 41,
                references: vec![reference.clone()],
                helpers: vec![],
                resources,
                summary_lines,
            },
            readiness: vec![dependency_readiness_row(&reference)],
            summary_lines: vec!["GENtle handoff ready".to_string()],
            suggested_actions,
            running_actions,
            blocked_actions: vec![],
            preferred_demo_actions: build_preferred_demo_actions(false),
            preferred_artifacts: vec![],
            environment_hints: vec![],
            warnings: vec![],
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

    #[test]
    fn reference_prepare_action_uses_canonical_resource_key() {
        let reference = fake_dependency(false, "missing", None);
        let action = reference_prepare_action(&reference, false);
        assert_eq!(action.kind, "prepare_reference");
        assert_eq!(
            action.shell_line,
            "genomes prepare \"Human GRCh38 Ensembl 116\" --timeout-secs 7200"
        );
        assert_eq!(
            action.resource_key.as_deref(),
            Some("reference_genome:Human GRCh38 Ensembl 116")
        );
        assert_eq!(action.lifecycle_status.as_deref(), Some("missing"));
        assert!(action.requires_confirmation);
    }

    #[test]
    fn telegram_guide_overview_uses_suggested_actions_as_section_links() {
        let reference = fake_dependency(false, "missing", None);
        let handoff = fake_handoff_report(reference, vec![], vec![]);
        let guide = build_telegram_guide_from_handoff(handoff, "telegram", "overview", None);

        assert_eq!(guide.schema, TELEGRAM_GUIDE_SCHEMA);
        assert_eq!(guide.section, "overview");
        assert_eq!(guide.gene, None);
        assert!(!guide.gene_supplied);
        assert!(guide
            .summary_lines
            .iter()
            .any(|line| line.contains("If you have a gene of interest")));
        assert!(guide
            .menu_sections
            .iter()
            .any(|section| section.section_id == "tfbs"));
        assert!(guide.suggested_actions.iter().any(|action| {
            action.kind == "guide_section"
                && action.label == "Continue with default genes"
                && action.shell_line == "services guide --channel telegram --section overview"
                && !action.requires_confirmation
        }));
        assert!(guide.suggested_actions.iter().any(|action| {
            action.kind == "guide_section"
                && action.shell_line == "services guide --channel telegram --section tfbs"
                && !action.requires_confirmation
        }));
    }

    #[test]
    fn telegram_guide_tfbs_defaults_to_tert_and_tp73_without_gene() {
        let reference = fake_dependency(false, "missing", None);
        let handoff = fake_handoff_report(reference, vec![], vec![]);
        let guide = build_telegram_guide_from_handoff(handoff, "telegram", "tfbs", None);

        assert_eq!(guide.section, "tfbs");
        let shell_lines = guide
            .suggested_actions
            .iter()
            .map(|action| action.shell_line.as_str())
            .collect::<Vec<_>>();
        assert!(shell_lines.iter().any(|line| {
            line.contains("promoter-tfbs-summary")
                && line.contains("--gene \"TERT\"")
                && line.contains("--gene \"TP73\"")
        }));
        assert!(shell_lines.iter().any(|line| {
            line.contains("promoter-tfbs-svg")
                && line.contains("--gene \"TERT\"")
                && line.contains("--gene \"TP73\"")
        }));
    }

    #[test]
    fn telegram_guide_tfbs_personalizes_gene_actions() {
        let reference = fake_dependency(false, "missing", None);
        let handoff = fake_handoff_report(reference, vec![], vec![]);
        let guide = build_telegram_guide_from_handoff(
            handoff,
            "telegram",
            "tfbs",
            Some("BACH2".to_string()),
        );

        assert_eq!(guide.gene.as_deref(), Some("BACH2"));
        assert!(guide.gene_supplied);
        assert!(guide.suggested_actions.iter().any(|action| {
            action.shell_line.contains("--gene \"BACH2\"")
                && !action.shell_line.contains("--gene \"TERT\"")
                && !action.shell_line.contains("--gene \"TP73\"")
        }));
    }

    #[test]
    fn telegram_guide_readiness_suppresses_prepare_when_target_is_running() {
        let reference = fake_dependency(
            false,
            "running",
            Some(fake_activity("running", "download_sequence", 37.0)),
        );
        let running_action = handoff_action(
            "Re-check Human GRCh38 Ensembl 116 status",
            "refresh_status",
            "genomes status \"Human GRCh38 Ensembl 116\"",
            180,
            "Reference is already being prepared.",
            false,
            Some("reference_genome:Human GRCh38 Ensembl 116".to_string()),
            Some("running".to_string()),
            vec![],
        );
        let handoff = fake_handoff_report(reference, vec![status_refresh_action()], vec![
            running_action,
        ]);
        let guide = build_telegram_guide_from_handoff(
            handoff,
            "telegram",
            "readiness",
            Some("TERT".to_string()),
        );

        assert!(guide
            .suggested_actions
            .iter()
            .any(|action| action.kind == "refresh_status"));
        assert!(!guide
            .suggested_actions
            .iter()
            .any(|action| action.kind == "prepare_reference"));
    }
}
