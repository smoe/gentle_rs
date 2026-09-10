//! Typed TFBS track-panel data with a fail-closed adapter for TSS profiles.
//!
//! The raw model permits heterogeneous policies for future faceted consumers;
//! TSS profiles explicitly require homogeneous score kind and display clipping.
//! No names, accessions, row order, or labels are inferred from a runtime preset.
//! The canonical adapter promotes identical per-track policies to panel scope
//! and rejects settings that cannot be represented without loss. A zero display
//! threshold denotes the supported unfiltered profile, not a score mutation.
//! Callers must retain the original bytes for the resolver's panel digest.
//!
//! Test provenance: the 30-track fixture is copied verbatim from the documented
//! Git object in `test_files/fixtures/tfbs_track_panel/README.md`. All other
//! examples are hand-crafted in-memory synthetic policies and PFMs; recreate
//! them by running this module's tests. No global registry override is used.

use crate::digest_utils::sha256_hex_bytes;
use gentle_engine::tss_profiles::validate_panel;
use gentle_protocol::EngineError;
use gentle_protocol::tss_profiles::{
    JasparPanelTrack, JasparTargetPanel, PANEL_SCHEMA, TssCalibrationState, TssScaleMode,
    TssStrandPolicy,
};
use serde::de::{self, MapAccess, SeqAccess, Visitor};
use serde::{Deserialize, Deserializer, Serialize};
use std::collections::BTreeSet;
use std::fmt;

const MAX_PANEL_BYTES: usize = 2 * 1024 * 1024;

/// Raw ordered panel data. Deserialization deliberately does not enforce TSS
/// score homogeneity; consumers must opt into the corresponding policy gate.
#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct TfbsTrackPanel {
    pub schema: String,
    pub tracks: Vec<TfbsPanelTrack>,
}

/// A raw track retains explicit source/factor memberships and display policies.
/// Policy strings remain data even when the current TSS renderer cannot use them.
#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct TfbsPanelTrack {
    pub track_id: String,
    pub label: String,
    pub provider_kind: String,
    pub source_ids: Vec<String>,
    pub factors: Vec<TfbsPanelFactor>,
    pub score_kind: String,
    pub calibration_state: String,
    pub calibration_statement: String,
    #[serde(default)]
    pub calibration_id: Option<String>,
    #[serde(default)]
    pub calibration_sha256: Option<String>,
    pub strand_policy: String,
    pub clip_negative: bool,
    pub display_threshold: f64,
    pub top_hit_count: usize,
    pub scale_mode: String,
    #[serde(default)]
    pub color_hint: Option<String>,
    pub display_order: usize,
}

/// Registry identity and display label are distinct, case-preserving fields.
#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct TfbsPanelFactor {
    pub factor_id: String,
    pub factor_label: String,
}

impl TfbsTrackPanel {
    /// Require one exact score kind and clipping policy, naming conflicting
    /// source accessions instead of silently converting or regrouping tracks.
    pub fn require_homogeneous_score_grammar(&self) -> Result<(), EngineError> {
        let first = self
            .tracks
            .first()
            .ok_or_else(|| invalid("TSS panel has no tracks"))?;
        for track in &self.tracks[1..] {
            if track.score_kind != first.score_kind || track.clip_negative != first.clip_negative {
                return Err(invalid(format!(
                    "Mixed TSS score grammar: {} ({}, clip_negative={}) differs from {} ({}, clip_negative={})",
                    track.source_ids.join(", "),
                    track.score_kind,
                    track.clip_negative,
                    first.source_ids.join(", "),
                    first.score_kind,
                    first.clip_negative,
                )));
            }
        }
        Ok(())
    }

    fn to_tss_panel(&self, input_sha256: &str) -> Result<JasparTargetPanel, EngineError> {
        self.require_homogeneous_score_grammar()?;
        let first = &self.tracks[0];
        let mut track_ids = BTreeSet::new();
        let mut factors = Vec::with_capacity(self.tracks.len());
        for track in &self.tracks {
            metadata_text(&track.track_id, "track_id")?;
            if !track_ids.insert(&track.track_id) {
                return Err(invalid(format!("Duplicate track_id {}", track.track_id)));
            }
            if track.source_ids.len() != 1 || track.factors.len() != 1 {
                return Err(invalid(format!(
                    "Track {} requires exactly one source_id and one factor for TSS profiles",
                    track.track_id,
                )));
            }
            let source = &track.source_ids[0];
            if track.provider_kind != "jaspar_pwm" {
                return Err(invalid(format!(
                    "Track {source}: unsupported provider_kind {}",
                    track.provider_kind
                )));
            }
            if !track.display_threshold.is_finite() || track.display_threshold != 0.0 {
                return Err(invalid(format!(
                    "Track {source}: display_threshold must be zero; nonzero filtering cannot be preserved",
                )));
            }
            for (field, equal) in [
                ("top_hit_count", track.top_hit_count == first.top_hit_count),
                ("scale_mode", track.scale_mode == first.scale_mode),
                ("strand_policy", track.strand_policy == first.strand_policy),
                (
                    "calibration_state",
                    track.calibration_state == first.calibration_state,
                ),
                (
                    "calibration_statement",
                    track.calibration_statement == first.calibration_statement,
                ),
                (
                    "calibration_id",
                    track.calibration_id == first.calibration_id,
                ),
                (
                    "calibration_sha256",
                    track.calibration_sha256 == first.calibration_sha256,
                ),
            ] {
                if !equal {
                    return Err(invalid(format!(
                        "Tracks {source} and {} differ in {field}; the TSS panel cannot preserve mixed values",
                        first.source_ids.join(", "),
                    )));
                }
            }
            let factor = &track.factors[0];
            metadata_text(&factor.factor_label, "factor_label")?;
            factors.push(JasparPanelTrack {
                source_id: source.clone(),
                factor_id: factor.factor_id.clone(),
                label: track.label.clone(),
                display_order: track.display_order,
                color_hint: track.color_hint.clone(),
                score_kind: Some(track.score_kind.clone()),
                track_id: Some(track.track_id.clone()),
                provider_kind: Some(track.provider_kind.clone()),
                factor_label: Some(factor.factor_label.clone()),
            });
        }
        let unsupported = |field: &str, value: &str| {
            invalid(format!(
                "Track {}: unsupported {field} {value}",
                first.source_ids.join(", ")
            ))
        };
        let panel = JasparTargetPanel {
            schema: self.schema.clone(),
            // The raw format declares no panel ID/title. Use a content identity
            // and a neutral title, never a name inferred from its factor list.
            panel_id: format!("sha256:{input_sha256}"),
            label: "JASPAR TFBS track panel".into(),
            score_kind: first.score_kind.clone(),
            clip_negative: first.clip_negative,
            scale_mode: match first.scale_mode.as_str() {
                "independent" => TssScaleMode::Independent,
                "shared" => TssScaleMode::Shared,
                value => return Err(unsupported("scale_mode", value)),
            },
            strand_policy: match first.strand_policy.as_str() {
                "both" => TssStrandPolicy::Both,
                value => return Err(unsupported("strand_policy", value)),
            },
            calibration_state: match first.calibration_state.as_str() {
                "matrix_specific" => TssCalibrationState::MatrixSpecific,
                "cross_source_calibrated" => TssCalibrationState::CrossSourceCalibrated,
                value => return Err(unsupported("calibration_state", value)),
            },
            calibration_statement: first.calibration_statement.clone(),
            calibration_id: first.calibration_id.clone(),
            calibration_sha256: first.calibration_sha256.clone(),
            top_hit_count: first.top_hit_count,
            factors,
        };
        validate_panel(&panel)?;
        Ok(panel)
    }
}

/// Parse either the real `tracks` schema or the provisional canonical `factors`
/// schema. Duplicate keys, ambiguous shapes and unsupported settings fail closed.
/// This performs no registry lookup; pass these same bytes to the strict resolver.
pub fn parse_tss_panel(bytes: &[u8]) -> Result<JasparTargetPanel, EngineError> {
    if bytes.len() > MAX_PANEL_BYTES {
        return Err(invalid("TSS panel JSON exceeds the 2 MiB input limit"));
    }
    serde_json::from_slice::<UniqueJsonKeys>(bytes)
        .map_err(|error| invalid(format!("Invalid TSS panel JSON: {error}")))?;
    let shape: serde_json::Value = serde_json::from_slice(bytes)
        .map_err(|error| invalid(format!("Invalid TSS panel JSON: {error}")))?;
    if shape.get("schema").and_then(serde_json::Value::as_str) != Some(PANEL_SCHEMA) {
        return Err(invalid(format!("Expected explicit {PANEL_SCHEMA} panel")));
    }
    match (
        shape.get("tracks").is_some(),
        shape.get("factors").is_some(),
    ) {
        (true, false) => {
            let raw: TfbsTrackPanel = serde_json::from_slice(bytes)
                .map_err(|error| invalid(format!("Invalid raw track panel: {error}")))?;
            raw.to_tss_panel(&sha256_hex_bytes(bytes))
        }
        (false, true) => {
            let panel: JasparTargetPanel = serde_json::from_slice(bytes)
                .map_err(|error| invalid(format!("Invalid canonical TSS panel: {error}")))?;
            validate_panel(&panel)?;
            Ok(panel)
        }
        _ => Err(invalid(
            "TSS panel must declare exactly one of tracks or factors",
        )),
    }
}

fn invalid(message: impl Into<String>) -> EngineError {
    EngineError::invalid_input(message)
}

fn metadata_text(value: &str, field: &str) -> Result<(), EngineError> {
    if value.trim().is_empty() || value.len() > 4096 || value.chars().any(char::is_control) {
        return Err(invalid(format!(
            "Invalid {field}: expected nonempty printable text within 4096 bytes"
        )));
    }
    Ok(())
}

// Inspect decoded keys recursively before shape detection can collapse objects.
struct UniqueJsonKeys;

impl<'de> Deserialize<'de> for UniqueJsonKeys {
    fn deserialize<D: Deserializer<'de>>(deserializer: D) -> Result<Self, D::Error> {
        struct UniqueVisitor;
        impl<'de> Visitor<'de> for UniqueVisitor {
            type Value = UniqueJsonKeys;

            fn expecting(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
                formatter.write_str("JSON with unique object keys")
            }

            fn visit_map<A: MapAccess<'de>>(self, mut map: A) -> Result<Self::Value, A::Error> {
                let mut keys = BTreeSet::new();
                while let Some(key) = map.next_key::<String>()? {
                    if !keys.insert(key.clone()) {
                        return Err(de::Error::custom(format!(
                            "duplicate JSON object key {key}"
                        )));
                    }
                    map.next_value::<UniqueJsonKeys>()?;
                }
                Ok(UniqueJsonKeys)
            }

            fn visit_seq<A: SeqAccess<'de>>(self, mut seq: A) -> Result<Self::Value, A::Error> {
                while seq.next_element::<UniqueJsonKeys>()?.is_some() {}
                Ok(UniqueJsonKeys)
            }

            fn visit_bool<E: de::Error>(self, _: bool) -> Result<Self::Value, E> {
                Ok(UniqueJsonKeys)
            }
            fn visit_i64<E: de::Error>(self, _: i64) -> Result<Self::Value, E> {
                Ok(UniqueJsonKeys)
            }
            fn visit_u64<E: de::Error>(self, _: u64) -> Result<Self::Value, E> {
                Ok(UniqueJsonKeys)
            }
            fn visit_f64<E: de::Error>(self, _: f64) -> Result<Self::Value, E> {
                Ok(UniqueJsonKeys)
            }
            fn visit_str<E: de::Error>(self, _: &str) -> Result<Self::Value, E> {
                Ok(UniqueJsonKeys)
            }
            fn visit_unit<E: de::Error>(self) -> Result<Self::Value, E> {
                Ok(UniqueJsonKeys)
            }
        }
        deserializer.deserialize_any(UniqueVisitor)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::engine::GentleEngine;
    use crate::tf_motifs::TfMotifDb;
    use serde_json::{Value, json};
    use std::collections::BTreeMap;

    const REAL_PANEL: &[u8] =
        include_bytes!("../test_files/fixtures/tfbs_track_panel/jaspar_30_track_panel.json");

    fn synthetic_panel() -> Value {
        json!({
            "schema": PANEL_SCHEMA,
            "tracks": (["MA9993.7", "MA9991.2"].into_iter().enumerate().map(|(i, id)| json!({
                "track_id": format!("synthetic-track-{i}"),
                "label": format!("Synthetic matrix row {i}"),
                "provider_kind": "jaspar_pwm",
                "source_ids": [id],
                "factors": [{"factor_id": "SyntheticFactor", "factor_label": "Synthetic factor label"}],
                "score_kind": "llr_background_tail_log10",
                "calibration_state": "matrix_specific",
                "calibration_statement": "Synthetic model scores, not biological evidence",
                "strand_policy": "both",
                "clip_negative": true,
                "display_threshold": 0,
                "top_hit_count": 5,
                "scale_mode": "independent",
                "color_hint": "#125a7C",
                "display_order": i * 10,
            })).collect::<Vec<_>>()),
        })
    }

    fn parse(value: &Value) -> Result<JasparTargetPanel, EngineError> {
        parse_tss_panel(&serde_json::to_vec(value).unwrap())
    }

    fn synthetic_registry() -> TfMotifDb {
        let text = serde_json::to_string(&json!({
            "schema": "gentle.tf_motifs.v1",
            "motifs": (["MA9993.7", "MA9991.2"].into_iter().map(|id| json!({
                "id": id, "name": "SyntheticFactor", "consensus_iupac": "AC",
                "pfm": {"a": [6,1], "c": [2,6], "g": [1,2], "t": [1,1]},
            })).collect::<Vec<_>>()),
        }))
        .unwrap();
        TfMotifDb::from_json_for_test(&text).unwrap()
    }

    #[test]
    fn real_panel_preserves_thirty_matrices_twenty_eight_factors_and_adjacent_triplet() {
        let raw: TfbsTrackPanel = serde_json::from_slice(REAL_PANEL).unwrap();
        raw.require_homogeneous_score_grammar().unwrap();
        let panel = parse_tss_panel(REAL_PANEL).unwrap();
        let db =
            TfMotifDb::from_json_for_test(include_str!("../assets/jaspar.motifs.json")).unwrap();
        let resolution = GentleEngine::resolve_tss_panel(panel, REAL_PANEL, &db).unwrap();
        assert_eq!(resolution.panel_sha256, sha256_hex_bytes(REAL_PANEL));
        assert_eq!(resolution.matrices.len(), 30);
        let mut factor_indices = BTreeMap::<&str, Vec<usize>>::new();
        for (index, (matrix, raw)) in resolution.matrices.iter().zip(&raw.tracks).enumerate() {
            let track = &matrix.specification;
            assert_eq!(track.source_id, raw.source_ids[0]);
            assert_eq!(track.factor_id, raw.factors[0].factor_id);
            assert_eq!(
                track.factor_label.as_deref(),
                Some(raw.factors[0].factor_label.as_str())
            );
            assert_eq!(track.track_id.as_deref(), Some(raw.track_id.as_str()));
            assert_eq!(
                track.provider_kind.as_deref(),
                Some(raw.provider_kind.as_str())
            );
            assert_eq!(track.label, raw.label);
            assert_eq!(track.display_order, raw.display_order);
            assert_eq!(track.color_hint, raw.color_hint);
            assert_eq!(track.score_kind.as_deref(), Some(raw.score_kind.as_str()));
            assert!(!matrix.matrix_counts.is_empty());
            assert_eq!(matrix.version, raw.source_ids[0].split_once('.').unwrap().1);
            assert_eq!(resolution.panel.score_kind, raw.score_kind);
            assert_eq!(resolution.panel.clip_negative, raw.clip_negative);
            assert_eq!(resolution.panel.top_hit_count, raw.top_hit_count);
            assert_eq!(
                resolution.panel.calibration_statement,
                raw.calibration_statement
            );
            assert_eq!(
                serde_json::to_value(resolution.panel.calibration_state).unwrap(),
                raw.calibration_state
            );
            assert_eq!(
                serde_json::to_value(resolution.panel.scale_mode).unwrap(),
                raw.scale_mode
            );
            assert_eq!(
                serde_json::to_value(resolution.panel.strand_policy).unwrap(),
                raw.strand_policy
            );
            factor_indices
                .entry(track.factor_id.as_str())
                .or_default()
                .push(index);
        }
        assert_eq!(factor_indices.len(), 28);
        let repeated = factor_indices
            .values()
            .filter(|indices| indices.len() > 1)
            .collect::<Vec<_>>();
        assert_eq!(repeated, [&vec![17, 18, 19]]);
        let triplet = &resolution.matrices[17..20];
        assert_eq!(
            triplet
                .iter()
                .map(|m| &m.specification.source_id)
                .collect::<BTreeSet<_>>()
                .len(),
            3
        );
        assert_eq!(
            triplet
                .iter()
                .map(|m| &m.matrix_sha256)
                .collect::<BTreeSet<_>>()
                .len(),
            3
        );
    }

    #[test]
    fn raw_model_permits_mixed_grammar_but_tss_gate_names_both_conflicting_accessions() {
        for (field, value) in [
            ("score_kind", json!("llr_bits")),
            ("clip_negative", json!(false)),
        ] {
            let mut value_panel = synthetic_panel();
            value_panel["tracks"][1][field] = value;
            let raw: TfbsTrackPanel = serde_json::from_value(value_panel.clone()).unwrap();
            assert_eq!(
                raw.tracks.len(),
                2,
                "generic mixed panels remain representable"
            );
            let error = raw.require_homogeneous_score_grammar().unwrap_err();
            assert!(error.message.contains("MA9993.7"));
            assert!(error.message.contains("MA9991.2"));
            assert!(error.message.contains("Mixed TSS score grammar"));
            assert!(parse(&value_panel).is_err());
        }
        let raw: TfbsTrackPanel = serde_json::from_value(synthetic_panel()).unwrap();
        raw.require_homogeneous_score_grammar().unwrap();
    }

    #[test]
    fn adaptation_preserves_distinct_labels_metadata_and_supported_shared_calibration() {
        let mut value = synthetic_panel();
        for track in value["tracks"].as_array_mut().unwrap() {
            track["score_kind"] = json!("true_log_odds_bits");
            track["clip_negative"] = json!(false);
            track["top_hit_count"] = json!(9);
            track["scale_mode"] = json!("shared");
            track["calibration_state"] = json!("cross_source_calibrated");
            track["calibration_id"] = json!("synthetic-calibration");
            track["calibration_sha256"] = json!(sha256_hex_bytes(b"synthetic calibration bytes"));
        }
        let panel = parse(&value).unwrap();
        assert_eq!(panel.score_kind, "true_log_odds_bits");
        assert!(!panel.clip_negative);
        assert_eq!(panel.top_hit_count, 9);
        assert_eq!(panel.scale_mode, TssScaleMode::Shared);
        assert_eq!(
            panel.calibration_state,
            TssCalibrationState::CrossSourceCalibrated
        );
        assert_eq!(
            panel.calibration_id.as_deref(),
            Some("synthetic-calibration")
        );
        assert_eq!(
            panel.calibration_sha256,
            Some(sha256_hex_bytes(b"synthetic calibration bytes"))
        );
        assert_eq!(panel.factors[0].factor_id, "SyntheticFactor");
        assert_eq!(
            panel.factors[0].factor_label.as_deref(),
            Some("Synthetic factor label")
        );
        assert_eq!(panel.factors[0].label, "Synthetic matrix row 0");
        assert_eq!(
            panel.factors[0].track_id.as_deref(),
            Some("synthetic-track-0")
        );
        assert_eq!(panel.factors[0].color_hint.as_deref(), Some("#125a7C"));
        assert_eq!(
            panel
                .factors
                .iter()
                .map(|t| t.display_order)
                .collect::<Vec<_>>(),
            [0, 10]
        );
    }

    #[test]
    fn unrepresentable_per_track_policies_fail_instead_of_taking_the_first_value() {
        for (field, value) in [
            ("top_hit_count", json!(7)),
            ("scale_mode", json!("shared")),
            ("strand_policy", json!("forward")),
            ("calibration_state", json!("cross_source_calibrated")),
            (
                "calibration_statement",
                json!("Different descriptive statement"),
            ),
            ("calibration_id", json!("different-calibration")),
            ("calibration_sha256", json!("a".repeat(64))),
            ("display_threshold", json!(0.01)),
            ("display_threshold", json!(-1)),
            ("provider_kind", json!("unsupported_provider")),
        ] {
            let mut panel = synthetic_panel();
            panel["tracks"][1][field] = value;
            let raw: TfbsTrackPanel = serde_json::from_value(panel.clone()).unwrap();
            raw.require_homogeneous_score_grammar().unwrap();
            let error = parse(&panel).unwrap_err();
            assert!(error.message.contains(field), "{}", error.message);
            assert!(error.message.contains("MA9991.2"), "{}", error.message);
        }
        for (field, value) in [
            ("score_kind", json!("future_score")),
            ("scale_mode", json!("future_scale")),
            ("strand_policy", json!("forward")),
            ("calibration_state", json!("future_calibration")),
            ("top_hit_count", json!(0)),
        ] {
            let mut panel = synthetic_panel();
            for track in panel["tracks"].as_array_mut().unwrap() {
                track[field] = value.clone();
            }
            assert!(parse(&panel).is_err(), "unsupported homogeneous {field}");
        }
    }

    #[test]
    fn raw_adapter_keeps_strict_accession_identity_checks_at_the_resolver_boundary() {
        let db = synthetic_registry();
        for source in [
            "SyntheticFactor",
            "MA9993",
            "ma9993.7",
            " MA9993.7",
            "MA9993.8",
            "MA9998.1",
            "MA9991.2",
        ] {
            let mut value = synthetic_panel();
            value["tracks"][0]["source_ids"] = json!([source]);
            let bytes = serde_json::to_vec(&value).unwrap();
            let result = parse_tss_panel(&bytes)
                .and_then(|panel| GentleEngine::resolve_tss_panel(panel, &bytes, &db));
            assert!(
                result.is_err(),
                "must not substitute or deduplicate {source}"
            );
        }
        let mut value = synthetic_panel();
        value["tracks"][0]["factors"][0]["factor_id"] = json!("syntheticfactor");
        let bytes = serde_json::to_vec(&value).unwrap();
        let error = GentleEngine::resolve_tss_panel(parse_tss_panel(&bytes).unwrap(), &bytes, &db)
            .unwrap_err();
        assert!(error.message.contains("SyntheticFactor"));
    }

    #[test]
    fn canonical_compatibility_and_original_raw_byte_bindings_are_preserved() {
        let bytes = include_bytes!("../test_files/fixtures/tss_profiles/panel.json");
        let canonical: JasparTargetPanel = serde_json::from_slice(bytes).unwrap();
        assert_eq!(
            serde_json::to_value(parse_tss_panel(bytes).unwrap()).unwrap(),
            serde_json::to_value(canonical).unwrap()
        );
        let raw = synthetic_panel();
        let compact = serde_json::to_vec(&raw).unwrap();
        let pretty = serde_json::to_vec_pretty(&raw).unwrap();
        let db = synthetic_registry();
        for bytes in [&compact, &pretty] {
            let panel = parse_tss_panel(bytes).unwrap();
            assert_eq!(
                panel.panel_id,
                format!("sha256:{}", sha256_hex_bytes(bytes))
            );
            let resolution = GentleEngine::resolve_tss_panel(panel, bytes, &db).unwrap();
            assert_eq!(resolution.panel_sha256, sha256_hex_bytes(bytes));
        }
        assert_ne!(sha256_hex_bytes(&compact), sha256_hex_bytes(&pretty));
    }

    #[test]
    fn duplicate_keys_are_rejected_at_every_depth_including_escaped_keys() {
        let raw = serde_json::to_string(&synthetic_panel()).unwrap();
        let canonical = include_str!("../test_files/fixtures/tss_profiles/panel.json");
        for text in [
            raw.replacen("\"schema\":", "\"schema\":\"ignored\",\"schema\":", 1),
            raw.replacen("\"schema\":", "\"\\u0073chema\":\"ignored\",\"schema\":", 1),
            raw.replacen(
                "\"clip_negative\":",
                "\"clip_negative\":false,\"clip_negative\":",
                1,
            ),
            raw.replacen(
                "\"factor_id\":",
                "\"factor_id\":\"wrong\",\"factor_id\":",
                1,
            ),
            canonical.replacen(
                "\"factor_id\":",
                "\"factor_id\":\"wrong\",\"factor_id\":",
                1,
            ),
        ] {
            let error = parse_tss_panel(text.as_bytes()).unwrap_err();
            assert!(
                error.message.contains("duplicate JSON object key"),
                "{}",
                error.message
            );
        }
    }

    #[test]
    fn unknown_ambiguous_and_lossy_shapes_fail_without_reordering_or_merging() {
        let base = synthetic_panel();
        for (field, value) in [
            ("source_ids", json!(["MA9993.7", "MA9992.1"])),
            ("source_ids", json!([])),
            ("factors", json!([])),
            (
                "factors",
                json!([
                    {"factor_id":"SyntheticFactor", "factor_label":"one"},
                    {"factor_id":"AnotherFactor", "factor_label":"two"},
                ]),
            ),
            ("smoothing", json!(3)),
            ("factor_label", json!("unexpected track field")),
            ("track_id", json!("")),
        ] {
            let mut panel = base.clone();
            panel["tracks"][0][field] = value;
            assert!(parse(&panel).is_err(), "unsupported {field}");
        }
        let mut duplicate = base.clone();
        duplicate["tracks"][1]["track_id"] = duplicate["tracks"][0]["track_id"].clone();
        assert!(
            parse(&duplicate)
                .unwrap_err()
                .message
                .contains("Duplicate track_id")
        );
        let mut reversed = base.clone();
        reversed["tracks"].as_array_mut().unwrap().reverse();
        assert!(
            parse(&reversed)
                .unwrap_err()
                .message
                .contains("display_order")
        );
        let mut ambiguous = base.clone();
        ambiguous["factors"] = json!([]);
        assert!(parse(&ambiguous).is_err());
        let mut unknown = base;
        unknown["future_setting"] = json!(true);
        assert!(parse(&unknown).is_err());
        for bytes in [
            b"null".as_slice(),
            b"{}".as_slice(),
            b"[]".as_slice(),
            b"{}{}".as_slice(),
        ] {
            assert!(parse_tss_panel(bytes).is_err());
        }
        assert!(parse_tss_panel(&vec![b' '; MAX_PANEL_BYTES + 1]).is_err());
    }
}
