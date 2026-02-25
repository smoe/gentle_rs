//! TF-motif registry and matching support.

use lazy_static::lazy_static;
use serde::Deserialize;
use std::{collections::HashMap, fs, sync::RwLock};

const RUNTIME_TF_MOTIF_PATH: &str = "data/resources/jaspar.motifs.json";
const BUILTIN_TF_MOTIFS_JSON: &str = include_str!("../assets/jaspar.motifs.json");

#[derive(Debug, Clone, Deserialize)]
struct TfMotifSnapshot {
    schema: String,
    motifs: Vec<TfMotifRecord>,
}

#[derive(Debug, Clone, Deserialize)]
struct TfPfmRows {
    a: Vec<f64>,
    c: Vec<f64>,
    g: Vec<f64>,
    t: Vec<f64>,
}

#[derive(Debug, Clone, Deserialize)]
struct TfMotifRecord {
    id: String,
    name: Option<String>,
    consensus_iupac: String,
    #[serde(default)]
    pfm: Option<TfPfmRows>,
}

#[derive(Debug, Clone)]
pub struct TfMotif {
    pub id: String,
    pub name: Option<String>,
    pub consensus_iupac: String,
    pub matrix_counts: Vec<[f64; 4]>,
}

#[derive(Debug, Clone)]
pub struct TfMotifSummary {
    pub id: String,
    pub name: Option<String>,
}

#[derive(Debug, Clone, Default)]
pub struct TfMotifDb {
    motifs: Vec<TfMotif>,
    by_key: HashMap<String, usize>,
}

fn iupac_counts(letter: char) -> [f64; 4] {
    match letter.to_ascii_uppercase() {
        'A' => [1.0, 0.0, 0.0, 0.0],
        'C' => [0.0, 1.0, 0.0, 0.0],
        'G' => [0.0, 0.0, 1.0, 0.0],
        'T' | 'U' => [0.0, 0.0, 0.0, 1.0],
        'M' => [1.0, 1.0, 0.0, 0.0], // A/C
        'R' => [1.0, 0.0, 1.0, 0.0], // A/G
        'W' => [1.0, 0.0, 0.0, 1.0], // A/T
        'S' => [0.0, 1.0, 1.0, 0.0], // C/G
        'Y' => [0.0, 1.0, 0.0, 1.0], // C/T
        'K' => [0.0, 0.0, 1.0, 1.0], // G/T
        'V' => [1.0, 1.0, 1.0, 0.0], // A/C/G
        'H' => [1.0, 1.0, 0.0, 1.0], // A/C/T
        'D' => [1.0, 0.0, 1.0, 1.0], // A/G/T
        'B' => [0.0, 1.0, 1.0, 1.0], // C/G/T
        _ => [1.0, 1.0, 1.0, 1.0],   // N or unknown
    }
}

impl TfMotifDb {
    fn matrix_from_pfm(pfm: &TfPfmRows) -> Option<Vec<[f64; 4]>> {
        let len = pfm.a.len();
        if len == 0 || pfm.c.len() != len || pfm.g.len() != len || pfm.t.len() != len {
            return None;
        }
        let mut matrix = Vec::with_capacity(len);
        for i in 0..len {
            matrix.push([pfm.a[i], pfm.c[i], pfm.g[i], pfm.t[i]]);
        }
        Some(matrix)
    }

    fn matrix_from_consensus(consensus: &str) -> Vec<[f64; 4]> {
        consensus.chars().map(iupac_counts).collect()
    }

    fn from_json(text: &str) -> Option<Self> {
        let snapshot = serde_json::from_str::<TfMotifSnapshot>(text).ok()?;
        if !snapshot.schema.starts_with("gentle.tf_motifs.v") {
            return None;
        }

        let mut motifs = Vec::new();
        let mut by_key = HashMap::new();
        for m in snapshot.motifs {
            let consensus = m.consensus_iupac.trim().to_ascii_uppercase();
            if consensus.is_empty() {
                continue;
            }
            let matrix_counts = match m.pfm.as_ref().and_then(Self::matrix_from_pfm) {
                Some(matrix) => matrix,
                None => Self::matrix_from_consensus(&consensus),
            };
            if matrix_counts.is_empty() {
                continue;
            }
            let idx = motifs.len();
            let id_key = m.id.trim().to_ascii_uppercase();
            let name_key = m.name.as_ref().map(|n| n.trim().to_ascii_uppercase());
            let motif = TfMotif {
                id: m.id.trim().to_string(),
                name: m.name.as_ref().map(|n| n.trim().to_string()),
                consensus_iupac: consensus,
                matrix_counts,
            };
            motifs.push(motif);
            if !id_key.is_empty() {
                by_key.insert(id_key, idx);
            }
            if let Some(name_key) = name_key {
                if !name_key.is_empty() {
                    by_key.entry(name_key).or_insert(idx);
                }
            }
        }
        Some(Self { motifs, by_key })
    }

    fn try_load_path(path: &str) -> Option<Self> {
        let text = fs::read_to_string(path).ok()?;
        Self::from_json(&text)
    }

    fn load_from_path(path: Option<&str>) -> Self {
        if let Some(db) = path.and_then(Self::try_load_path) {
            return db;
        }
        if let Some(db) = Self::try_load_path(RUNTIME_TF_MOTIF_PATH) {
            return db;
        }
        Self::from_json(BUILTIN_TF_MOTIFS_JSON).unwrap_or_default()
    }

    fn load() -> Self {
        Self::load_from_path(None)
    }

    pub fn resolve(&self, token: &str) -> Option<&TfMotif> {
        let key = token.trim().to_ascii_uppercase();
        let idx = self.by_key.get(&key)?;
        self.motifs.get(*idx)
    }
}

lazy_static! {
    static ref TF_MOTIFS: RwLock<TfMotifDb> = RwLock::new(TfMotifDb::load());
}

pub fn resolve_motif(token: &str) -> Option<String> {
    TF_MOTIFS
        .read()
        .ok()
        .and_then(|db| db.resolve(token).map(|m| m.consensus_iupac.clone()))
}

pub fn resolve_motif_definition(token: &str) -> Option<TfMotif> {
    TF_MOTIFS
        .read()
        .ok()
        .and_then(|db| db.resolve(token).cloned())
}

pub fn reload() {
    reload_from_path(None);
}

pub fn reload_from_path(path: Option<&str>) {
    if let Ok(mut db) = TF_MOTIFS.write() {
        *db = TfMotifDb::load_from_path(path);
    }
}

pub fn list_motif_summaries() -> Vec<TfMotifSummary> {
    let mut out = TF_MOTIFS
        .read()
        .ok()
        .map(|db| {
            db.motifs
                .iter()
                .map(|m| TfMotifSummary {
                    id: m.id.clone(),
                    name: m.name.clone(),
                })
                .collect::<Vec<_>>()
        })
        .unwrap_or_default();
    out.sort_by(|a, b| {
        a.id.to_ascii_uppercase()
            .cmp(&b.id.to_ascii_uppercase())
            .then_with(|| {
                a.name
                    .as_deref()
                    .unwrap_or("")
                    .to_ascii_uppercase()
                    .cmp(&b.name.as_deref().unwrap_or("").to_ascii_uppercase())
            })
    });
    out
}

pub fn all_motif_ids() -> Vec<String> {
    list_motif_summaries().into_iter().map(|m| m.id).collect()
}
