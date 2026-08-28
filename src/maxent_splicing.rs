//! User-supplied MaxEntScan splice-site model snapshots and native scoring.
//!
//! GENtle does not bundle or download the upstream probability tables. The
//! resource sync path normalizes a user-provided MaxEntScan directory/archive,
//! records its provenance, and this module reproduces the original 9 nt donor
//! and 23 nt acceptor scoring formula without invoking Perl at analysis time.

use crate::digest_utils::sha256_prefixed_bytes;
use gentle_protocol::CrypticSplicingModelTableDigest;
use serde::{Deserialize, Serialize};
use std::{
    collections::HashMap,
    fs,
    path::Path,
    sync::{Arc, LazyLock, RwLock},
};

pub const MAXENT_SPLICE_MODEL_SCHEMA: &str = "gentle.maxent_splice_model.v1";
pub const DEFAULT_MAXENT_SPLICE_MODEL_PATH: &str = "data/resources/maxent_splice_model.json";

#[derive(Debug, Clone, Default, PartialEq, Serialize, Deserialize)]
#[serde(default)]
pub struct MaxEntSpliceModelSnapshot {
    pub schema: String,
    pub source: String,
    pub source_url: Option<String>,
    pub retrieved_on: Option<String>,
    pub redistribution_status: String,
    pub model_fingerprint_sha256: String,
    pub table_sha256: Vec<CrypticSplicingModelTableDigest>,
    pub donor_sequences: Vec<String>,
    pub donor_scores: Vec<f64>,
    pub acceptor_tables: Vec<Vec<f64>>,
}

#[derive(Debug)]
pub struct MaxEntSpliceModel {
    snapshot: MaxEntSpliceModelSnapshot,
    donor_index: HashMap<String, usize>,
}

#[derive(Debug, Clone)]
pub struct ActiveMaxEntModelStatus {
    pub path: String,
    pub exists: bool,
    pub model: Option<Arc<MaxEntSpliceModel>>,
    pub error: Option<String>,
}

static ACTIVE_MODEL: LazyLock<RwLock<Option<ActiveMaxEntModelStatus>>> =
    LazyLock::new(|| RwLock::new(None));

fn base4_index(sequence: &str) -> Result<usize, String> {
    let mut index = 0usize;
    for base in sequence.bytes() {
        let digit = match base.to_ascii_uppercase() {
            b'A' => 0,
            b'C' => 1,
            b'G' => 2,
            b'T' => 3,
            other => {
                return Err(format!(
                    "MaxEnt window contains unsupported base '{}'",
                    char::from(other)
                ));
            }
        };
        index = index.saturating_mul(4).saturating_add(digit);
    }
    Ok(index)
}

fn normalized_sequence(sequence: &str, expected_len: usize) -> Result<String, String> {
    let normalized = sequence.trim().to_ascii_uppercase().replace('U', "T");
    if normalized.len() != expected_len {
        return Err(format!(
            "MaxEnt window requires {expected_len} nt, received {}",
            normalized.len()
        ));
    }
    base4_index(&normalized)?;
    Ok(normalized)
}

fn consensus_ratio(base: u8, probabilities: [f64; 4]) -> Result<f64, String> {
    let (index, background) = match base.to_ascii_uppercase() {
        b'A' => (0, 0.27),
        b'C' => (1, 0.23),
        b'G' => (2, 0.23),
        b'T' => (3, 0.27),
        other => {
            return Err(format!(
                "MaxEnt consensus contains unsupported base '{}'",
                char::from(other)
            ));
        }
    };
    Ok(probabilities[index] / background)
}

fn model_content_sha256(snapshot: &MaxEntSpliceModelSnapshot) -> Result<String, String> {
    #[derive(Serialize)]
    struct ModelContent<'a> {
        donor_sequences: &'a [String],
        donor_scores: &'a [f64],
        acceptor_tables: &'a [Vec<f64>],
    }
    let bytes = serde_json::to_vec(&ModelContent {
        donor_sequences: &snapshot.donor_sequences,
        donor_scores: &snapshot.donor_scores,
        acceptor_tables: &snapshot.acceptor_tables,
    })
    .map_err(|error| format!("Could not fingerprint normalized MaxEnt tables: {error}"))?;
    Ok(sha256_prefixed_bytes(&bytes))
}

impl MaxEntSpliceModelSnapshot {
    pub fn finalize_fingerprint(&mut self) -> Result<(), String> {
        self.model_fingerprint_sha256 = model_content_sha256(self)?;
        Ok(())
    }

    pub fn validate(&self) -> Result<(), String> {
        if self.schema != MAXENT_SPLICE_MODEL_SCHEMA {
            return Err(format!(
                "Unexpected MaxEnt snapshot schema '{}'; expected '{}'",
                self.schema, MAXENT_SPLICE_MODEL_SCHEMA
            ));
        }
        if self.donor_sequences.len() != 16_384 || self.donor_scores.len() != 16_384 {
            return Err(format!(
                "MaxEnt donor tables require 16384 sequence/score rows, received {}/{}",
                self.donor_sequences.len(),
                self.donor_scores.len()
            ));
        }
        let expected_acceptor_lengths = [16_384, 16_384, 16_384, 16_384, 16_384, 64, 256, 64, 256];
        if self.acceptor_tables.len() != expected_acceptor_lengths.len() {
            return Err(format!(
                "MaxEnt acceptor model requires 9 tables, received {}",
                self.acceptor_tables.len()
            ));
        }
        for (index, (table, expected_len)) in self
            .acceptor_tables
            .iter()
            .zip(expected_acceptor_lengths)
            .enumerate()
        {
            if table.len() != expected_len {
                return Err(format!(
                    "MaxEnt acceptor table {} requires {} rows, received {}",
                    index + 1,
                    expected_len,
                    table.len()
                ));
            }
        }
        let mut seen = HashMap::<&str, usize>::new();
        for (index, sequence) in self.donor_sequences.iter().enumerate() {
            normalized_sequence(sequence, 7)?;
            if let Some(previous) = seen.insert(sequence.as_str(), index) {
                return Err(format!(
                    "MaxEnt donor sequence '{}' is duplicated at rows {} and {}",
                    sequence,
                    previous + 1,
                    index + 1
                ));
            }
        }
        for value in self
            .donor_scores
            .iter()
            .chain(self.acceptor_tables.iter().flatten())
        {
            if !value.is_finite() || *value <= 0.0 {
                return Err("MaxEnt tables must contain finite positive values".to_string());
            }
        }
        let computed = model_content_sha256(self)?;
        if computed != self.model_fingerprint_sha256 {
            return Err(format!(
                "MaxEnt snapshot fingerprint mismatch: recorded '{}' but normalized tables hash to '{}'",
                self.model_fingerprint_sha256, computed
            ));
        }
        Ok(())
    }
}

impl MaxEntSpliceModel {
    pub fn from_snapshot(snapshot: MaxEntSpliceModelSnapshot) -> Result<Self, String> {
        snapshot.validate()?;
        let donor_index = snapshot
            .donor_sequences
            .iter()
            .enumerate()
            .map(|(index, sequence)| (sequence.clone(), index))
            .collect();
        Ok(Self {
            snapshot,
            donor_index,
        })
    }

    pub fn snapshot(&self) -> &MaxEntSpliceModelSnapshot {
        &self.snapshot
    }

    /// Reproduce `score5.pl` for one 9 nt `-3..+6` donor window.
    pub fn score_donor(&self, sequence: &str) -> Result<f64, String> {
        let sequence = normalized_sequence(sequence, 9)?;
        let bytes = sequence.as_bytes();
        let rest = format!(
            "{}{}{}{}{}{}{}",
            char::from(bytes[0]),
            char::from(bytes[1]),
            char::from(bytes[2]),
            char::from(bytes[5]),
            char::from(bytes[6]),
            char::from(bytes[7]),
            char::from(bytes[8])
        );
        let index = self
            .donor_index
            .get(&rest)
            .copied()
            .ok_or_else(|| format!("MaxEnt donor sequence table does not contain '{rest}'"))?;
        let consensus = consensus_ratio(bytes[3], [0.004, 0.0032, 0.9896, 0.0032])?
            * consensus_ratio(bytes[4], [0.0034, 0.0039, 0.0042, 0.9884])?;
        let score = (consensus * self.snapshot.donor_scores[index]).log2();
        if !score.is_finite() {
            return Err("MaxEnt donor score is not finite".to_string());
        }
        Ok(score)
    }

    /// Reproduce `score3.pl` for one 23 nt `-20..+3` acceptor window.
    pub fn score_acceptor(&self, sequence: &str) -> Result<f64, String> {
        let sequence = normalized_sequence(sequence, 23)?;
        let bytes = sequence.as_bytes();
        let rest = format!("{}{}", &sequence[..18], &sequence[20..]);
        let slices = [
            (0, 7),
            (7, 14),
            (14, 21),
            (4, 11),
            (11, 18),
            (4, 7),
            (7, 11),
            (11, 14),
            (14, 18),
        ];
        let mut scores = [0.0f64; 9];
        for (index, (start, end)) in slices.into_iter().enumerate() {
            let table_index = base4_index(&rest[start..end])?;
            scores[index] = self.snapshot.acceptor_tables[index][table_index];
        }
        let consensus = consensus_ratio(bytes[18], [0.9903, 0.0032, 0.0034, 0.0030])?
            * consensus_ratio(bytes[19], [0.0027, 0.0037, 0.9905, 0.0030])?;
        let model = scores[0] * scores[1] * scores[2] * scores[3] * scores[4]
            / (scores[5] * scores[6] * scores[7] * scores[8]);
        let score = (consensus * model).log2();
        if !score.is_finite() {
            return Err("MaxEnt acceptor score is not finite".to_string());
        }
        Ok(score)
    }
}

fn load_path(path: &str) -> ActiveMaxEntModelStatus {
    if !Path::new(path).is_file() {
        return ActiveMaxEntModelStatus {
            path: path.to_string(),
            exists: false,
            model: None,
            error: None,
        };
    }
    let loaded = fs::read_to_string(path)
        .map_err(|error| format!("Could not read MaxEnt snapshot '{path}': {error}"))
        .and_then(|text| {
            serde_json::from_str::<MaxEntSpliceModelSnapshot>(&text)
                .map_err(|error| format!("Could not parse MaxEnt snapshot '{path}': {error}"))
        })
        .and_then(MaxEntSpliceModel::from_snapshot)
        .map(Arc::new);
    match loaded {
        Ok(model) => ActiveMaxEntModelStatus {
            path: path.to_string(),
            exists: true,
            model: Some(model),
            error: None,
        },
        Err(error) => ActiveMaxEntModelStatus {
            path: path.to_string(),
            exists: true,
            model: None,
            error: Some(error),
        },
    }
}

pub fn reload_from_path(path: Option<&str>) -> ActiveMaxEntModelStatus {
    let status = load_path(path.unwrap_or(DEFAULT_MAXENT_SPLICE_MODEL_PATH));
    if let Ok(mut active) = ACTIVE_MODEL.write() {
        *active = Some(status.clone());
    }
    status
}

pub fn active_model_status() -> ActiveMaxEntModelStatus {
    if let Ok(active) = ACTIVE_MODEL.read()
        && let Some(status) = active.as_ref()
    {
        return status.clone();
    }
    reload_from_path(None)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn base4_sequence(mut index: usize, len: usize) -> String {
        let mut bases = vec![b'A'; len];
        for position in (0..len).rev() {
            bases[position] = match index % 4 {
                0 => b'A',
                1 => b'C',
                2 => b'G',
                _ => b'T',
            };
            index /= 4;
        }
        String::from_utf8(bases).expect("base4 DNA")
    }

    pub(super) fn synthetic_snapshot(scale: f64) -> MaxEntSpliceModelSnapshot {
        let mut snapshot = MaxEntSpliceModelSnapshot {
            schema: MAXENT_SPLICE_MODEL_SCHEMA.to_string(),
            source: "synthetic unit-test tables".to_string(),
            redistribution_status: "synthetic_no_restriction".to_string(),
            donor_sequences: (0..16_384).map(|index| base4_sequence(index, 7)).collect(),
            donor_scores: vec![scale; 16_384],
            acceptor_tables: [16_384, 16_384, 16_384, 16_384, 16_384, 64, 256, 64, 256]
                .into_iter()
                .map(|len| vec![scale; len])
                .collect(),
            ..MaxEntSpliceModelSnapshot::default()
        };
        snapshot.finalize_fingerprint().expect("fingerprint");
        snapshot
    }

    #[test]
    fn base4_index_matches_original_maxent_order() {
        assert_eq!(base4_index("AAAA").expect("index"), 0);
        assert_eq!(base4_index("AAAC").expect("index"), 1);
        assert_eq!(base4_index("TTTT").expect("index"), 255);
        assert!(base4_index("NAAA").is_err());
    }

    #[test]
    fn native_maxent_scoring_reproduces_consensus_and_table_formulae() {
        let model =
            MaxEntSpliceModel::from_snapshot(synthetic_snapshot(1.0)).expect("synthetic model");
        let donor = model.score_donor("AAAGTAAAA").expect("donor score");
        let expected_donor = ((0.9896 / 0.23_f64) * (0.9884 / 0.27_f64)).log2();
        assert!((donor - expected_donor).abs() < 1e-12);

        let acceptor = model
            .score_acceptor("AAAAAAAAAAAAAAAAAAAGAAA")
            .expect("acceptor score");
        let expected_acceptor = ((0.9903 / 0.27_f64) * (0.9905 / 0.23_f64)).log2();
        assert!((acceptor - expected_acceptor).abs() < 1e-12);
    }

    #[test]
    fn model_fingerprint_changes_with_tables_and_rejects_corruption() {
        let first = synthetic_snapshot(1.0);
        let second = synthetic_snapshot(2.0);
        assert_ne!(
            first.model_fingerprint_sha256,
            second.model_fingerprint_sha256
        );

        let mut corrupted = first;
        corrupted.donor_scores[0] = 2.0;
        let error = MaxEntSpliceModel::from_snapshot(corrupted).expect_err("stale fingerprint");
        assert!(error.contains("fingerprint mismatch"), "{error}");
    }
}
