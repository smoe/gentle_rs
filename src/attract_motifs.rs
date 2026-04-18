//! ATtRACT RNA-binding motif registry and runtime snapshot helpers.
//!
//! This module owns the normalized runtime representation of motifs imported
//! from the published ATtRACT ZIP download. The first integration milestone
//! loaded exact/IUPAC consensus motifs with provenance intact; the current
//! milestone optionally retains per-matrix PWM/PFM rows too so splice-aware
//! interpretation can upgrade to PWM-backed scoring without inventing a second
//! resource schema.

use serde::{Deserialize, Serialize};
use sha1::{Digest, Sha1};
use std::{
    fs,
    sync::{LazyLock, RwLock},
};

pub const ATTRACT_MOTIF_SNAPSHOT_SCHEMA: &str = "gentle.attract_motifs.v1";
pub const DEFAULT_ATTRACT_RESOURCE_PATH: &str = "data/resources/attract.motifs.json";

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
pub struct AttractPfmRows {
    pub a: Vec<f64>,
    pub c: Vec<f64>,
    pub g: Vec<f64>,
    pub t: Vec<f64>,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
pub struct AttractMotifRecord {
    pub entry_id: String,
    pub matrix_id: String,
    pub gene_name: String,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub gene_id: Option<String>,
    pub organism: String,
    pub motif_iupac: String,
    pub length: usize,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub experiment: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub family: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub domain: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub pubmed_id: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub quality_score: Option<f64>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub source_database: Option<String>,
    #[serde(default = "default_model_kind")]
    pub model_kind: String,
    #[serde(default)]
    pub pwm_present: bool,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub pfm: Option<AttractPfmRows>,
}

fn default_model_kind() -> String {
    "consensus_iupac".to_string()
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
pub struct AttractMotifSnapshot {
    pub schema: String,
    pub source: String,
    pub fetched_at_unix_ms: u128,
    pub motif_count: usize,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub snapshot_fingerprint: Option<String>,
    #[serde(default)]
    pub archive_members: Vec<String>,
    #[serde(default)]
    pub warnings: Vec<String>,
    #[serde(default)]
    pub motifs: Vec<AttractMotifRecord>,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
pub struct AttractMotifSummary {
    pub entry_id: String,
    pub matrix_id: String,
    pub gene_name: String,
    pub organism: String,
    pub motif_iupac: String,
    pub model_kind: String,
    pub pwm_present: bool,
}

#[derive(Debug, Clone, Default)]
pub struct AttractMotifDb {
    snapshot: Option<AttractMotifSnapshot>,
    snapshot_fingerprint: Option<String>,
}

impl AttractMotifDb {
    fn from_json(text: &str) -> Option<Self> {
        let mut snapshot = serde_json::from_str::<AttractMotifSnapshot>(text).ok()?;
        if !snapshot.schema.starts_with("gentle.attract_motifs.v") {
            return None;
        }
        let fingerprint = snapshot_fingerprint_from_text(text);
        if snapshot.snapshot_fingerprint.is_none() {
            snapshot.snapshot_fingerprint = fingerprint.clone();
        }
        Some(Self {
            snapshot: Some(snapshot),
            snapshot_fingerprint: fingerprint,
        })
    }

    fn try_load_path(path: &str) -> Option<Self> {
        let text = fs::read_to_string(path).ok()?;
        Self::from_json(&text)
    }

    fn load_from_path(path: Option<&str>) -> Self {
        path.and_then(Self::try_load_path)
            .or_else(|| Self::try_load_path(DEFAULT_ATTRACT_RESOURCE_PATH))
            .unwrap_or_default()
    }

    fn all_motifs(&self) -> Vec<AttractMotifRecord> {
        self.snapshot
            .as_ref()
            .map(|snapshot| snapshot.motifs.clone())
            .unwrap_or_default()
    }

    fn active_snapshot_fingerprint(&self) -> Option<String> {
        self.snapshot_fingerprint.clone()
    }

    fn list_summaries(&self) -> Vec<AttractMotifSummary> {
        let mut rows = self
            .snapshot
            .as_ref()
            .map(|snapshot| {
                snapshot
                    .motifs
                    .iter()
                    .map(|motif| AttractMotifSummary {
                        entry_id: motif.entry_id.clone(),
                        matrix_id: motif.matrix_id.clone(),
                        gene_name: motif.gene_name.clone(),
                        organism: motif.organism.clone(),
                        motif_iupac: motif.motif_iupac.clone(),
                        model_kind: motif.model_kind.clone(),
                        pwm_present: motif.pwm_present,
                    })
                    .collect::<Vec<_>>()
            })
            .unwrap_or_default();
        rows.sort_by(|a, b| {
            a.gene_name
                .to_ascii_uppercase()
                .cmp(&b.gene_name.to_ascii_uppercase())
                .then_with(|| {
                    a.organism
                        .to_ascii_uppercase()
                        .cmp(&b.organism.to_ascii_uppercase())
                })
                .then_with(|| {
                    a.matrix_id
                        .to_ascii_uppercase()
                        .cmp(&b.matrix_id.to_ascii_uppercase())
                })
        });
        rows
    }
}

static ATTRACT_MOTIFS: LazyLock<RwLock<AttractMotifDb>> =
    LazyLock::new(|| RwLock::new(AttractMotifDb::load_from_path(None)));

pub fn reload() {
    reload_from_path(None);
}

pub fn reload_from_path(path: Option<&str>) {
    if let Ok(mut db) = ATTRACT_MOTIFS.write() {
        *db = AttractMotifDb::load_from_path(path);
    }
}

pub fn list_motif_summaries() -> Vec<AttractMotifSummary> {
    ATTRACT_MOTIFS
        .read()
        .ok()
        .map(|db| db.list_summaries())
        .unwrap_or_default()
}

pub fn all_motifs() -> Vec<AttractMotifRecord> {
    ATTRACT_MOTIFS
        .read()
        .ok()
        .map(|db| db.all_motifs())
        .unwrap_or_default()
}

pub fn active_snapshot_fingerprint() -> Option<String> {
    ATTRACT_MOTIFS
        .read()
        .ok()
        .and_then(|db| db.active_snapshot_fingerprint())
}

pub fn snapshot_fingerprint_from_text(text: &str) -> Option<String> {
    let trimmed = text.trim();
    if trimmed.is_empty() {
        return None;
    }
    let mut hasher = Sha1::new();
    hasher.update(trimmed.as_bytes());
    Some(format!("sha1:{:x}", hasher.finalize()))
}
