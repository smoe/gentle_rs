//! Runtime status helpers for built-in and overrideable external resources.
//!
//! This module answers "what resource set is GENtle actually using right now?"
//! for shipped registries such as REBASE and JASPAR, and records clearly when a
//! future/planned external source is not yet integrated.

use crate::{
    attract_motifs::{
        DEFAULT_ATTRACT_RESOURCE_PATH, active_model_counts as active_attract_model_counts,
        active_snapshot_fingerprint as active_attract_snapshot_fingerprint,
        list_motif_summaries as list_attract_motif_summaries, snapshot_fingerprint_from_text,
    },
    enzymes::load_restriction_enzymes_from_json_text,
    resource_sync::DEFAULT_JASPAR_REMOTE_METADATA_PATH,
    tf_motifs::list_motif_summaries,
};
use serde::Serialize;
use serde_json::Value;
use std::{
    fs,
    time::{SystemTime, UNIX_EPOCH},
};

const BUILTIN_REBASE_LABEL: &str = "assets/enzymes.json";
const BUILTIN_REBASE_JSON: &str = include_str!("../assets/enzymes.json");
const BUILTIN_JASPAR_LABEL: &str = "assets/jaspar.motifs.json";
const BUILTIN_JASPAR_JSON: &str = include_str!("../assets/jaspar.motifs.json");
const RUNTIME_REBASE_PATH: &str = "data/resources/rebase.enzymes.json";
const RUNTIME_JASPAR_PATH: &str = "data/resources/jaspar.motifs.json";
const ATTRACT_INDEX_URL: &str = "https://attract.cnic.es/index";
const ATTRACT_DOWNLOAD_URL: &str = "https://attract.cnic.es/attract/static/ATtRACT.zip";
const RESOURCE_STATUS_SCHEMA: &str = "gentle.resource_status.v1";

#[derive(Debug, Clone, Serialize)]
pub struct ResourceCatalogReport {
    pub schema: String,
    pub generated_at_unix_ms: u128,
    pub rebase: ResourceSnapshotStatus,
    pub jaspar: ResourceSnapshotStatus,
    pub attract: AttractResourceStatus,
}

#[derive(Debug, Clone, Serialize)]
pub struct ResourceSnapshotStatus {
    pub resource_id: String,
    pub builtin_label: String,
    pub builtin_available: bool,
    pub builtin_item_count: usize,
    pub runtime_path: String,
    pub runtime_exists: bool,
    pub runtime_valid: bool,
    pub runtime_item_count: Option<usize>,
    pub runtime_error: Option<String>,
    pub active_source: String,
    pub active_item_count: usize,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub remote_metadata_snapshot: Option<DerivedSnapshotStatus>,
}

#[derive(Debug, Clone, Serialize)]
pub struct DerivedSnapshotStatus {
    pub path: String,
    pub exists: bool,
    pub valid: bool,
    pub item_count: Option<usize>,
    pub error: Option<String>,
}

#[derive(Debug, Clone, Serialize)]
pub struct ExternalDatabaseStatus {
    pub resource_id: String,
    pub display_name: String,
    pub support_status: String,
    pub homepage: String,
    pub download_url: Option<String>,
    pub notes: Vec<String>,
}

#[derive(Debug, Clone, Serialize)]
pub struct AttractResourceStatus {
    pub resource_id: String,
    pub display_name: String,
    pub support_status: String,
    pub homepage: String,
    pub download_url: Option<String>,
    pub runtime_path: String,
    pub runtime_exists: bool,
    pub runtime_valid: bool,
    pub runtime_item_count: Option<usize>,
    pub runtime_fingerprint: Option<String>,
    pub runtime_error: Option<String>,
    pub active_source: String,
    pub active_item_count: usize,
    pub active_pwm_row_count: usize,
    pub active_consensus_only_row_count: usize,
    pub active_fingerprint: Option<String>,
    pub notes: Vec<String>,
}

fn now_unix_ms() -> u128 {
    SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .map(|d| d.as_millis())
        .unwrap_or(0)
}

fn count_jaspar_snapshot_items(text: &str) -> Result<usize, String> {
    let parsed: Value =
        serde_json::from_str(text).map_err(|e| format!("Could not parse motif JSON: {e}"))?;
    let schema = parsed
        .get("schema")
        .and_then(|v| v.as_str())
        .unwrap_or_default()
        .to_string();
    if !schema.starts_with("gentle.tf_motifs.v") {
        return Err(format!("Unexpected motif schema '{schema}'"));
    }
    if let Some(count) = parsed.get("motif_count").and_then(|v| v.as_u64()) {
        return Ok(count as usize);
    }
    Ok(parsed
        .get("motifs")
        .and_then(|v| v.as_array())
        .map(|rows| rows.len())
        .unwrap_or(0))
}

fn count_jaspar_remote_metadata_snapshot_items(text: &str) -> Result<usize, String> {
    let parsed: Value = serde_json::from_str(text)
        .map_err(|e| format!("Could not parse JASPAR remote-metadata JSON: {e}"))?;
    let schema = parsed
        .get("schema")
        .and_then(|v| v.as_str())
        .unwrap_or_default()
        .to_string();
    if !schema.starts_with("gentle.jaspar_remote_metadata_snapshot.v") {
        return Err(format!(
            "Unexpected JASPAR remote-metadata schema '{schema}'"
        ));
    }
    if let Some(count) = parsed.get("persisted_entry_count").and_then(|v| v.as_u64()) {
        return Ok(count as usize);
    }
    Ok(parsed
        .get("rows")
        .and_then(|v| v.as_array())
        .map(|rows| rows.len())
        .unwrap_or(0))
}

fn rebase_status() -> ResourceSnapshotStatus {
    let builtin_item_count =
        load_restriction_enzymes_from_json_text(BUILTIN_REBASE_JSON).map_or(0, |v| v.len());
    let runtime_exists = std::path::Path::new(RUNTIME_REBASE_PATH).exists();
    let (runtime_valid, runtime_item_count, runtime_error, active_source, active_item_count) =
        if runtime_exists {
            match fs::read_to_string(RUNTIME_REBASE_PATH)
                .map_err(|e| format!("Could not read runtime REBASE snapshot: {e}"))
                .and_then(|text| {
                    load_restriction_enzymes_from_json_text(&text)
                        .map(|enzymes| enzymes.len())
                        .map_err(|e| format!("Could not parse runtime REBASE snapshot: {e}"))
                }) {
                Ok(count) if count > 0 => (true, Some(count), None, "runtime".to_string(), count),
                Ok(count) => (
                    false,
                    Some(count),
                    Some("Runtime REBASE snapshot is empty".to_string()),
                    "builtin".to_string(),
                    builtin_item_count,
                ),
                Err(error) => (
                    false,
                    None,
                    Some(error),
                    "builtin".to_string(),
                    builtin_item_count,
                ),
            }
        } else {
            (false, None, None, "builtin".to_string(), builtin_item_count)
        };
    ResourceSnapshotStatus {
        resource_id: "rebase".to_string(),
        builtin_label: BUILTIN_REBASE_LABEL.to_string(),
        builtin_available: builtin_item_count > 0,
        builtin_item_count,
        runtime_path: RUNTIME_REBASE_PATH.to_string(),
        runtime_exists,
        runtime_valid,
        runtime_item_count,
        runtime_error,
        active_source,
        active_item_count,
        remote_metadata_snapshot: None,
    }
}

fn jaspar_status() -> ResourceSnapshotStatus {
    let builtin_item_count = count_jaspar_snapshot_items(BUILTIN_JASPAR_JSON).unwrap_or(0);
    let runtime_exists = std::path::Path::new(RUNTIME_JASPAR_PATH).exists();
    let (runtime_valid, runtime_item_count, runtime_error) = if runtime_exists {
        match fs::read_to_string(RUNTIME_JASPAR_PATH)
            .map_err(|e| format!("Could not read runtime JASPAR snapshot: {e}"))
            .and_then(|text| count_jaspar_snapshot_items(&text))
        {
            Ok(count) if count > 0 => (true, Some(count), None),
            Ok(count) => (
                false,
                Some(count),
                Some("Runtime JASPAR snapshot is empty".to_string()),
            ),
            Err(error) => (false, None, Some(error)),
        }
    } else {
        (false, None, None)
    };
    let active_item_count = list_motif_summaries().len();
    let active_source = if runtime_valid {
        "runtime".to_string()
    } else {
        "builtin".to_string()
    };
    let remote_metadata_exists = std::path::Path::new(DEFAULT_JASPAR_REMOTE_METADATA_PATH).exists();
    let (remote_metadata_valid, remote_metadata_item_count, remote_metadata_error) =
        if remote_metadata_exists {
            match fs::read_to_string(DEFAULT_JASPAR_REMOTE_METADATA_PATH)
                .map_err(|e| format!("Could not read JASPAR remote-metadata snapshot: {e}"))
                .and_then(|text| count_jaspar_remote_metadata_snapshot_items(&text))
            {
                Ok(count) if count > 0 => (true, Some(count), None),
                Ok(count) => (
                    false,
                    Some(count),
                    Some("JASPAR remote-metadata snapshot is empty".to_string()),
                ),
                Err(error) => (false, None, Some(error)),
            }
        } else {
            (false, None, None)
        };
    ResourceSnapshotStatus {
        resource_id: "jaspar".to_string(),
        builtin_label: BUILTIN_JASPAR_LABEL.to_string(),
        builtin_available: builtin_item_count > 0,
        builtin_item_count,
        runtime_path: RUNTIME_JASPAR_PATH.to_string(),
        runtime_exists,
        runtime_valid,
        runtime_item_count,
        runtime_error,
        active_source,
        active_item_count,
        remote_metadata_snapshot: Some(DerivedSnapshotStatus {
            path: DEFAULT_JASPAR_REMOTE_METADATA_PATH.to_string(),
            exists: remote_metadata_exists,
            valid: remote_metadata_valid,
            item_count: remote_metadata_item_count,
            error: remote_metadata_error,
        }),
    }
}

fn count_attract_snapshot_items(text: &str) -> Result<usize, String> {
    let parsed: Value =
        serde_json::from_str(text).map_err(|e| format!("Could not parse ATtRACT JSON: {e}"))?;
    let schema = parsed
        .get("schema")
        .and_then(|v| v.as_str())
        .unwrap_or_default()
        .to_string();
    if !schema.starts_with("gentle.attract_motifs.v") {
        return Err(format!("Unexpected ATtRACT schema '{schema}'"));
    }
    if let Some(count) = parsed.get("motif_count").and_then(|v| v.as_u64()) {
        return Ok(count as usize);
    }
    Ok(parsed
        .get("motifs")
        .and_then(|v| v.as_array())
        .map(|rows| rows.len())
        .unwrap_or(0))
}

fn attract_status() -> AttractResourceStatus {
    let runtime_exists = std::path::Path::new(DEFAULT_ATTRACT_RESOURCE_PATH).exists();
    let (runtime_valid, runtime_item_count, runtime_fingerprint, runtime_error) = if runtime_exists
    {
        match fs::read_to_string(DEFAULT_ATTRACT_RESOURCE_PATH)
            .map_err(|e| format!("Could not read runtime ATtRACT snapshot: {e}"))
            .and_then(|text| {
                let count = count_attract_snapshot_items(&text)?;
                Ok((count, snapshot_fingerprint_from_text(&text)))
            }) {
            Ok((count, fingerprint)) if count > 0 => (true, Some(count), fingerprint, None),
            Ok((count, fingerprint)) => (
                false,
                Some(count),
                fingerprint,
                Some("Runtime ATtRACT snapshot is empty".to_string()),
            ),
            Err(error) => (false, None, None, Some(error)),
        }
    } else {
        (false, None, None, None)
    };
    let active_item_count = list_attract_motif_summaries().len();
    let (active_pwm_row_count, active_consensus_only_row_count) = active_attract_model_counts();
    let active_fingerprint = active_attract_snapshot_fingerprint();
    let (support_status, active_source, notes) = if runtime_valid {
        (
            "runtime_snapshot".to_string(),
            "runtime".to_string(),
            vec![
                "ATtRACT motifs are available from a normalized runtime snapshot.".to_string(),
                "Current v1 integration retains normalized consensus/IUPAC motifs from ATtRACT_db.txt and, when PWM blocks can be mapped by Matrix_id, upgrades those rows to PWM-backed splice-aware scoring.".to_string(),
                format!(
                    "Active model mix: {} PWM-backed rows, {} consensus-only rows.",
                    active_pwm_row_count, active_consensus_only_row_count
                ),
                format!("Published ZIP download: {ATTRACT_DOWNLOAD_URL}"),
            ],
        )
    } else if active_item_count > 0 {
        (
            "session_override".to_string(),
            "session_override".to_string(),
            vec![
                "ATtRACT motifs are currently loaded in this GENtle session from a non-default snapshot path.".to_string(),
                format!("Published ZIP download: {ATTRACT_DOWNLOAD_URL}"),
                format!(
                    "Persist the normalized snapshot at '{}' if you want `resources status` / `services status` to pick it up as the default runtime source too.",
                    DEFAULT_ATTRACT_RESOURCE_PATH
                ),
            ],
        )
    } else {
        (
            "known_external_only".to_string(),
            "unavailable".to_string(),
            vec![
                "ATtRACT is an RNA-binding protein and motif database.".to_string(),
                "GENtle now knows how to normalize the published ZIP into a runtime snapshot, but no valid snapshot is active yet.".to_string(),
                format!("Published ZIP download: {ATTRACT_DOWNLOAD_URL}"),
                "Use `resources sync-attract ATtRACT.zip` to turn the downloaded archive into GENtle's runtime snapshot.".to_string(),
            ],
        )
    };
    AttractResourceStatus {
        resource_id: "attract".to_string(),
        display_name: "ATtRACT".to_string(),
        support_status,
        homepage: ATTRACT_INDEX_URL.to_string(),
        download_url: Some(ATTRACT_DOWNLOAD_URL.to_string()),
        runtime_path: DEFAULT_ATTRACT_RESOURCE_PATH.to_string(),
        runtime_exists,
        runtime_valid,
        runtime_item_count,
        runtime_fingerprint,
        runtime_error,
        active_source,
        active_item_count,
        active_pwm_row_count,
        active_consensus_only_row_count,
        active_fingerprint,
        notes,
    }
}

pub fn resource_catalog_status() -> ResourceCatalogReport {
    ResourceCatalogReport {
        schema: RESOURCE_STATUS_SCHEMA.to_string(),
        generated_at_unix_ms: now_unix_ms(),
        rebase: rebase_status(),
        jaspar: jaspar_status(),
        attract: attract_status(),
    }
}
