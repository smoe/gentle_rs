//! TF-motif registry and matching support.

use serde::Deserialize;
use std::{
    collections::{BTreeSet, HashMap},
    fs,
    sync::{LazyLock, RwLock},
};

const RUNTIME_TF_MOTIF_PATH: &str = "data/resources/jaspar.motifs.json";
const BUILTIN_TF_MOTIFS_JSON: &str = include_str!("../assets/jaspar.motifs.json");
const BUNDLED_JASPAR_PFM_SUPPLEMENT_JSON: &str = include_str!("../assets/jaspar_2022.json");
const TF_QUERY_STOP_WORDS: &[&str] = &[
    "tf",
    "tfs",
    "factor",
    "factors",
    "family",
    "families",
    "group",
    "groups",
    "motif",
    "motifs",
    "binding",
    "association",
    "associations",
];

const COMMON_MOTIF_ALIASES: &[(&str, &str)] = &[
    ("OCT4", "POU5F1"),
    ("OCT-4", "POU5F1"),
    ("OCT 4", "POU5F1"),
    ("C-MYC", "MYC"),
    ("C MYC", "MYC"),
    ("CMYC", "MYC"),
    ("NRSF", "REST"),
];

#[derive(Debug, Clone)]
pub struct TfQueryGroupDefinition {
    pub id: String,
    pub label: String,
    pub description: String,
    pub aliases: Vec<String>,
    pub member_queries: Vec<String>,
}

#[derive(Debug, Clone)]
pub struct TfQueryResolvedMotif {
    pub motif_id: String,
    pub motif_name: Option<String>,
    pub consensus_iupac: String,
    pub motif_length_bp: usize,
}

#[derive(Debug, Clone)]
pub struct TfQueryResolution {
    pub query: String,
    pub normalized_query: String,
    pub resolution_kind: String,
    pub label: Option<String>,
    pub description: Option<String>,
    pub aliases: Vec<String>,
    pub matches: Vec<TfQueryResolvedMotif>,
    pub unresolved_reason: Option<String>,
}

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

#[derive(Debug, Clone, Deserialize)]
struct BundledJasparPfmRows {
    #[serde(rename = "A")]
    a: Vec<String>,
    #[serde(rename = "C")]
    c: Vec<String>,
    #[serde(rename = "G")]
    g: Vec<String>,
    #[serde(rename = "T")]
    t: Vec<String>,
}

#[derive(Debug, Clone, Deserialize)]
struct BundledJasparPfmRecord {
    matrix_id: String,
    name: Option<String>,
    pfm: BundledJasparPfmRows,
}

#[derive(Debug, Clone)]
struct TfPfmSupplement {
    id: String,
    name: Option<String>,
    consensus_iupac: String,
    matrix_counts: Vec<[f64; 4]>,
}

#[derive(Debug, Clone, Default)]
struct TfPfmSupplements {
    by_id: HashMap<String, TfPfmSupplement>,
    by_unique_name: HashMap<String, TfPfmSupplement>,
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
    pub consensus_iupac: String,
    pub motif_length_bp: usize,
}

fn normalize_lookup_key(raw: &str) -> String {
    raw.chars()
        .map(|ch| {
            if ch.is_ascii_alphanumeric() {
                ch.to_ascii_uppercase()
            } else {
                ' '
            }
        })
        .collect::<String>()
        .split_whitespace()
        .collect::<Vec<_>>()
        .join(" ")
}

fn normalize_query_lookup(raw: &str) -> String {
    let normalized = normalize_lookup_key(raw).to_ascii_lowercase();
    let words = normalized
        .split_whitespace()
        .filter(|word| !TF_QUERY_STOP_WORDS.contains(word))
        .collect::<Vec<_>>();
    if words.is_empty() {
        normalized.trim().to_string()
    } else {
        words.join(" ")
    }
}

fn normalize_dense_lookup(raw: &str) -> String {
    normalize_query_lookup(raw)
        .chars()
        .filter(|ch| ch.is_ascii_alphanumeric())
        .collect()
}

fn common_alias_targets() -> HashMap<String, String> {
    COMMON_MOTIF_ALIASES
        .iter()
        .map(|(alias, target)| (normalize_lookup_key(alias), normalize_lookup_key(target)))
        .collect()
}

fn catalog_group_by_query(query: &str) -> Option<crate::gene_groups::LoadedGeneGroupRecord> {
    let normalized = normalize_query_lookup(query);
    if normalized.is_empty() {
        return None;
    }
    crate::gene_groups::list_tf_query_gene_groups()
        .into_iter()
        .find(|group| {
            normalize_query_lookup(&group.record.id) == normalized
                || normalize_query_lookup(&group.record.label) == normalized
                || group
                    .record
                    .aliases
                    .iter()
                    .any(|alias| normalize_query_lookup(alias) == normalized)
        })
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

fn consensus_from_matrix_counts(matrix_counts: &[[f64; 4]]) -> String {
    const BASES: [char; 4] = ['A', 'C', 'G', 'T'];
    matrix_counts
        .iter()
        .map(|counts| {
            counts
                .iter()
                .copied()
                .enumerate()
                .max_by(|left, right| left.1.total_cmp(&right.1))
                .map(|(idx, _)| BASES[idx])
                .unwrap_or('N')
        })
        .collect()
}

impl TfMotifDb {
    fn bundled_jaspar_pfm_matrix(pfm: &BundledJasparPfmRows) -> Option<Vec<[f64; 4]>> {
        let len = pfm.a.len();
        if len == 0 || pfm.c.len() != len || pfm.g.len() != len || pfm.t.len() != len {
            return None;
        }
        let mut matrix = Vec::with_capacity(len);
        for i in 0..len {
            matrix.push([
                pfm.a[i].parse::<f64>().ok()?,
                pfm.c[i].parse::<f64>().ok()?,
                pfm.g[i].parse::<f64>().ok()?,
                pfm.t[i].parse::<f64>().ok()?,
            ]);
        }
        Some(matrix)
    }

    fn bundled_jaspar_pfm_supplements() -> TfPfmSupplements {
        let records =
            serde_json::from_str::<Vec<BundledJasparPfmRecord>>(BUNDLED_JASPAR_PFM_SUPPLEMENT_JSON)
                .unwrap_or_default();
        let mut by_id = HashMap::new();
        let mut name_candidates: HashMap<String, Option<TfPfmSupplement>> = HashMap::new();
        for record in records {
            let matrix_counts = match Self::bundled_jaspar_pfm_matrix(&record.pfm) {
                Some(matrix_counts) => matrix_counts,
                None => continue,
            };
            let id = record.matrix_id.trim().to_string();
            if id.is_empty() {
                continue;
            }
            let supplement = TfPfmSupplement {
                id: id.clone(),
                name: record.name.as_ref().map(|name| name.trim().to_string()),
                consensus_iupac: consensus_from_matrix_counts(&matrix_counts),
                matrix_counts,
            };
            by_id
                .entry(normalize_lookup_key(&id))
                .or_insert_with(|| supplement.clone());
            if let Some(name) = supplement.name.as_deref() {
                let name_key = normalize_lookup_key(name);
                if !name_key.is_empty() {
                    name_candidates
                        .entry(name_key)
                        .and_modify(|slot| *slot = None)
                        .or_insert_with(|| Some(supplement));
                }
            }
        }
        let by_unique_name = name_candidates
            .into_iter()
            .filter_map(|(name, supplement)| supplement.map(|supplement| (name, supplement)))
            .collect();
        TfPfmSupplements {
            by_id,
            by_unique_name,
        }
    }

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

        static PFM_SUPPLEMENTS: LazyLock<TfPfmSupplements> =
            LazyLock::new(TfMotifDb::bundled_jaspar_pfm_supplements);
        let mut motifs = Vec::new();
        let mut by_key = HashMap::new();
        let alias_targets = common_alias_targets();
        for m in snapshot.motifs {
            let compact_consensus = m.consensus_iupac.trim().to_ascii_uppercase();
            if compact_consensus.is_empty() {
                continue;
            }
            let original_id = m.id.trim().to_string();
            let original_id_key = normalize_lookup_key(&original_id);
            let name_key = m.name.as_ref().map(|name| normalize_lookup_key(name));
            let pfm_matrix = m.pfm.as_ref().and_then(Self::matrix_from_pfm);
            let supplement = if pfm_matrix.is_none() {
                PFM_SUPPLEMENTS.by_id.get(&original_id_key).or_else(|| {
                    name_key
                        .as_ref()
                        .and_then(|name_key| PFM_SUPPLEMENTS.by_unique_name.get(name_key))
                })
            } else {
                None
            };
            let (consensus, matrix_counts, supplement_id_key) = match (pfm_matrix, supplement) {
                (Some(matrix), _) => (compact_consensus, matrix, None),
                (None, Some(supplement)) => (
                    supplement.consensus_iupac.clone(),
                    supplement.matrix_counts.clone(),
                    Some(normalize_lookup_key(&supplement.id)),
                ),
                (None, None) => (
                    compact_consensus.clone(),
                    Self::matrix_from_consensus(&compact_consensus),
                    None,
                ),
            };
            if matrix_counts.is_empty() {
                continue;
            }
            let idx = motifs.len();
            let id_key = normalize_lookup_key(&original_id);
            let motif = TfMotif {
                id: original_id.clone(),
                name: m.name.as_ref().map(|n| n.trim().to_string()),
                consensus_iupac: consensus,
                matrix_counts,
            };
            motifs.push(motif);
            if !id_key.is_empty() {
                by_key.insert(id_key.clone(), idx);
            }
            if let Some(supplement_id_key) = supplement_id_key
                && supplement_id_key != id_key
                && !supplement_id_key.is_empty()
            {
                by_key.entry(supplement_id_key).or_insert(idx);
            }
            if let Some(name_key) = name_key
                && !name_key.is_empty()
            {
                by_key.entry(name_key).or_insert(idx);
            }
            for (alias, target) in &alias_targets {
                if *target == normalize_lookup_key(&m.id)
                    || m.name
                        .as_deref()
                        .map(normalize_lookup_key)
                        .map(|value| value == *target)
                        .unwrap_or(false)
                {
                    by_key.entry(alias.clone()).or_insert(idx);
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
        Self::load_builtin()
    }

    fn load_builtin() -> Self {
        Self::from_json(BUILTIN_TF_MOTIFS_JSON).unwrap_or_default()
    }

    fn load() -> Self {
        Self::load_from_path(None)
    }

    pub fn resolve(&self, token: &str) -> Option<&TfMotif> {
        let key = normalize_lookup_key(token);
        let idx = self.by_key.get(&key)?;
        self.motifs.get(*idx)
    }

    pub fn resolve_cloned(&self, token: &str) -> Option<TfMotif> {
        self.resolve(token).cloned()
    }

    pub fn motif_summaries(&self) -> Vec<TfMotifSummary> {
        let mut out = self
            .motifs
            .iter()
            .map(|m| TfMotifSummary {
                id: m.id.clone(),
                name: m.name.clone(),
                consensus_iupac: m.consensus_iupac.clone(),
                motif_length_bp: m.matrix_counts.len(),
            })
            .collect::<Vec<_>>();
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

    pub fn motif_ids(&self) -> Vec<String> {
        self.motif_summaries().into_iter().map(|m| m.id).collect()
    }
}

static TF_MOTIFS: LazyLock<RwLock<TfMotifDb>> = LazyLock::new(|| RwLock::new(TfMotifDb::load()));

#[cfg(test)]
pub fn test_registry_lock() -> &'static std::sync::Mutex<()> {
    static LOCK: std::sync::OnceLock<std::sync::Mutex<()>> = std::sync::OnceLock::new();
    LOCK.get_or_init(|| std::sync::Mutex::new(()))
}

pub fn snapshot_db() -> TfMotifDb {
    TF_MOTIFS.read().map(|db| db.clone()).unwrap_or_default()
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

#[cfg(test)]
pub fn reload_builtin_for_test() {
    if let Ok(mut db) = TF_MOTIFS.write() {
        *db = TfMotifDb::load_builtin();
    }
}

pub fn list_motif_summaries() -> Vec<TfMotifSummary> {
    snapshot_db().motif_summaries()
}

pub fn all_motif_ids() -> Vec<String> {
    snapshot_db().motif_ids()
}

fn fuzzy_family_like_matches(query: &str) -> Vec<TfQueryResolvedMotif> {
    let normalized_dense = normalize_dense_lookup(query);
    if normalized_dense.len() < 2 {
        return vec![];
    }
    list_motif_summaries()
        .into_iter()
        .filter(|summary| {
            let id_dense = normalize_dense_lookup(&summary.id);
            let name_dense = summary
                .name
                .as_deref()
                .map(normalize_dense_lookup)
                .unwrap_or_default();
            id_dense.starts_with(&normalized_dense)
                || name_dense.starts_with(&normalized_dense)
                || (normalized_dense.len() >= 3
                    && (id_dense.contains(&normalized_dense)
                        || name_dense.contains(&normalized_dense)))
        })
        .map(|summary| TfQueryResolvedMotif {
            motif_id: summary.id,
            motif_name: summary.name,
            consensus_iupac: summary.consensus_iupac,
            motif_length_bp: summary.motif_length_bp,
        })
        .collect()
}

pub fn resolve_tf_query(query: &str) -> TfQueryResolution {
    let trimmed = query.trim().to_string();
    let normalized_query = normalize_query_lookup(&trimmed);
    if trimmed.is_empty() {
        return TfQueryResolution {
            query: trimmed,
            normalized_query,
            resolution_kind: "unresolved".to_string(),
            label: None,
            description: None,
            aliases: vec![],
            matches: vec![],
            unresolved_reason: Some("Empty TF query".to_string()),
        };
    }

    if let Some(motif) = resolve_motif_definition(&trimmed) {
        return TfQueryResolution {
            query: trimmed,
            normalized_query,
            resolution_kind: "exact_motif".to_string(),
            label: Some(
                motif.name
                    .clone()
                    .unwrap_or_else(|| motif.id.clone()),
            ),
            description: Some(
                "Resolved directly against the local motif registry by motif id, TF name, or common alias.".to_string(),
            ),
            aliases: vec![],
            matches: vec![TfQueryResolvedMotif {
                motif_id: motif.id,
                motif_name: motif.name,
                consensus_iupac: motif.consensus_iupac,
                motif_length_bp: motif.matrix_counts.len(),
            }],
            unresolved_reason: None,
        };
    }

    if let Some(group) = catalog_group_by_query(&trimmed) {
        let mut matches = vec![];
        let mut seen = BTreeSet::new();
        for member in &group.record.members {
            if let Some(motif) = resolve_motif_definition(&member.symbol)
                && seen.insert(motif.id.clone())
            {
                matches.push(TfQueryResolvedMotif {
                    motif_id: motif.id,
                    motif_name: motif.name,
                    consensus_iupac: motif.consensus_iupac,
                    motif_length_bp: motif.matrix_counts.len(),
                });
            }
        }
        return TfQueryResolution {
            query: trimmed,
            normalized_query,
            resolution_kind: if group.source_scope == "built-in" {
                "builtin_group".to_string()
            } else {
                "gene_group_catalog".to_string()
            },
            label: Some(group.record.label),
            description: Some(group.record.long_definition),
            aliases: group.record.aliases,
            matches,
            unresolved_reason: None,
        };
    }

    let family_matches = fuzzy_family_like_matches(&trimmed);
    if !family_matches.is_empty() {
        return TfQueryResolution {
            query: trimmed,
            normalized_query,
            resolution_kind: "family_like_match".to_string(),
            label: Some(query.trim().to_string()),
            description: Some(
                "Resolved by local motif-id / TF-name family-style matching against the active motif registry.".to_string(),
            ),
            aliases: vec![],
            matches: family_matches,
            unresolved_reason: None,
        };
    }

    TfQueryResolution {
        query: trimmed,
        normalized_query,
        resolution_kind: "unresolved".to_string(),
        label: None,
        description: None,
        aliases: vec![],
        matches: vec![],
        unresolved_reason: Some(
            "Query did not match a local motif, a built-in TF group, or a family-like local motif-name prefix.".to_string(),
        ),
    }
}

pub fn expand_tf_query_to_motif_ids(query: &str) -> Vec<String> {
    let resolution = resolve_tf_query(query);
    let mut seen = BTreeSet::new();
    resolution
        .matches
        .into_iter()
        .filter_map(|row| seen.insert(row.motif_id.clone()).then_some(row.motif_id))
        .collect()
}

pub fn list_tf_query_groups() -> Vec<TfQueryGroupDefinition> {
    crate::gene_groups::list_tf_query_gene_groups()
        .into_iter()
        .map(|group| TfQueryGroupDefinition {
            id: group.record.id,
            label: group.record.label,
            description: group.record.long_definition,
            aliases: group.record.aliases,
            member_queries: group
                .record
                .members
                .into_iter()
                .map(|member| member.symbol)
                .collect(),
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn resolve_tf_query_supports_exact_aliases() {
        let resolution = resolve_tf_query("OCT4");
        assert_eq!(resolution.resolution_kind, "exact_motif");
        assert_eq!(resolution.matches.len(), 1);
        assert_eq!(resolution.matches[0].motif_name.as_deref(), Some("POU5F1"));
    }

    #[test]
    fn resolve_tf_query_supports_yamanaka_group() {
        let resolution = resolve_tf_query("Yamanaka factors");
        assert_eq!(resolution.resolution_kind, "builtin_group");
        let ids = resolution
            .matches
            .iter()
            .map(|row| {
                row.motif_name
                    .clone()
                    .unwrap_or_else(|| row.motif_id.clone())
            })
            .collect::<Vec<_>>();
        assert!(ids.contains(&"POU5F1".to_string()));
        assert!(ids.contains(&"SOX2".to_string()));
        assert!(ids.contains(&"KLF4".to_string()));
        assert!(ids.contains(&"MYC".to_string()));
    }

    #[test]
    fn resolve_tf_query_supports_family_like_prefixes() {
        let resolution = resolve_tf_query("KLF family");
        assert_eq!(resolution.resolution_kind, "family_like_match");
        assert!(
            resolution
                .matches
                .iter()
                .any(|row| row.motif_name.as_deref() == Some("KLF4"))
        );
    }

    #[test]
    fn resolve_motif_supports_versioned_jaspar_ids() {
        assert!(resolve_motif("MA0001.3").is_some());
    }

    #[test]
    fn builtin_tp73_uses_bundled_pfm_not_consensus_fallback() {
        let db = TfMotifDb::from_json(BUILTIN_TF_MOTIFS_JSON).expect("motif db");
        let resolved = db.resolve("TP73").expect("resolve TP73");
        assert_eq!(resolved.id, "MA0861.2");
        assert_eq!(resolved.matrix_counts.len(), 18);
        assert_eq!(resolved.consensus_iupac, "GACATGTCTGGACATGTC");

        let ambiguous_column = resolved.matrix_counts[8];
        let total = ambiguous_column.iter().sum::<f64>();
        let max_fraction = ambiguous_column
            .iter()
            .copied()
            .map(|value| value / total)
            .fold(0.0_f64, f64::max);
        assert!(
            max_fraction < 0.5,
            "position 9 should remain mixed, not collapsed to a consensus base"
        );

        let compact_id_alias = db.resolve("MA0861.2").expect("resolve compact alias");
        assert_eq!(compact_id_alias.id, "MA0861.2");
        let supplement_id_alias = db.resolve("MA0861.1").expect("resolve supplement alias");
        assert_eq!(supplement_id_alias.id, "MA0861.2");
    }

    #[test]
    fn bundled_pfm_supplement_does_not_replace_canonical_catalog_ids() {
        let db = TfMotifDb::from_json(BUILTIN_TF_MOTIFS_JSON).expect("motif db");
        let myc = db.resolve("MYC").expect("resolve MYC");
        assert_eq!(myc.id, "MA0147.4");
        assert_eq!(
            db.resolve("MA0147.3").map(|motif| motif.id.as_str()),
            Some("MA0147.4")
        );
        let sp1_ids = db
            .motif_summaries()
            .into_iter()
            .filter(|motif| motif.id == "MA0079.5")
            .collect::<Vec<_>>();
        assert_eq!(sp1_ids.len(), 1);
        assert_eq!(sp1_ids[0].name.as_deref(), Some("SP1"));
    }

    #[test]
    fn resolve_motif_supports_names_with_underscores() {
        let db = TfMotifDb::from_json(
            r#"{
  "schema":"gentle.tf_motifs.v1",
  "motifs":[
    {
      "id":"MTEST1.1",
      "name":"CODEX_TEST_MOTIF_1",
      "consensus_iupac":"ACGT"
    }
  ]
}"#,
        )
        .expect("motif db");
        let resolved = db
            .resolve("CODEX_TEST_MOTIF_1")
            .expect("resolve custom motif");
        assert_eq!(resolved.consensus_iupac, "ACGT");
    }
}
