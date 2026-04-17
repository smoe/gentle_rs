//! Runtime status helpers for built-in and overrideable external resources.
//!
//! This module answers "what resource set is GENtle actually using right now?"
//! for shipped registries such as REBASE and JASPAR, and records clearly when a
//! future/planned external source is not yet integrated.

use crate::{enzymes::load_restriction_enzymes_from_json_text, tf_motifs::list_motif_summaries};
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
    pub attract: ExternalDatabaseStatus,
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
    }
}

fn attract_status() -> ExternalDatabaseStatus {
    ExternalDatabaseStatus {
        resource_id: "attract".to_string(),
        display_name: "ATtRACT".to_string(),
        support_status: "not_yet_integrated".to_string(),
        homepage: ATTRACT_INDEX_URL.to_string(),
        download_url: Some(ATTRACT_DOWNLOAD_URL.to_string()),
        notes: vec![
            "ATtRACT is an RNA-binding protein and motif database.".to_string(),
            "GENtle does not yet import or score against ATtRACT snapshots.".to_string(),
            format!("Published ZIP download: {ATTRACT_DOWNLOAD_URL}"),
            "Current service-readiness reporting covers integrated resources only; ATtRACT is listed here so callers can stay explicit about that gap.".to_string(),
        ],
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
