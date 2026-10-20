//! Portable measured-gel profiles and descriptive comparisons. Pixel intensity
//! is never relabelled as concentration, molecule counts or independent samples.

use super::gel_image::{GelImageLane, GelMigrationDirection, GelSampleBand, GelSizeKind};
use serde::{Deserialize, Serialize};
use std::{collections::BTreeMap, sync::Arc};

#[derive(Debug, Clone, Copy, Serialize, Deserialize, PartialEq, Eq)]
#[serde(rename_all = "snake_case")]
pub enum GelEvidenceKind {
    ImageDisplay,
    VendorTrace,
    SignedResidual,
    /// Explicitly invented teaching/test integrals, never detector evidence.
    Synthetic,
}

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
#[serde(deny_unknown_fields)]
pub struct GelProfile {
    pub id: String,
    pub label: String,
    pub run_id: String,
    pub lane_id: String,
    pub source_sha256: String,
    pub evidence_kind: GelEvidenceKind,
    pub size_kind: GelSizeKind,
    pub signal_units: String,
    /// Integral of signal over each [edge[i], edge[i+1]] size bin. Values
    /// already carry their integration measure; they are not densities per bp.
    pub edges: Vec<f64>,
    pub weights: Vec<Option<f64>>,
    pub measure: String,
    pub processing: Vec<String>,
    #[serde(default)]
    pub warnings: Vec<String>,
}

/// Explicit integrated-bin CSV import. Point-sampled vendor CSV is not silently
/// integrated against a different axis or assigned an unknown calibration.
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
#[serde(deny_unknown_fields)]
pub struct GelProfileFileRequest {
    pub path: String,
    pub profile_id: String,
    pub label: String,
    pub run_id: String,
    pub lane_id: String,
    pub size_kind: GelSizeKind,
    pub signal_units: String,
    pub measure: String,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct GelNumericSource {
    pub source_name: String,
    pub sha256: String,
    pub original_base64: String,
    /// Replay metadata; path is reduced to its basename before persistence.
    pub request: GelProfileFileRequest,
}

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
pub struct GelNumericSourceSummary {
    pub profile_id: String,
    pub source_name: String,
    pub sha256: String,
    pub signal_units: String,
    pub measure: String,
}

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
#[serde(deny_unknown_fields)]
pub struct GelSourceImportRequest {
    pub image_id: String,
    pub path: String,
    /// PDF raster object number, required if several supported images exist.
    #[serde(default)]
    pub pdf_object: Option<u32>,
    /// Native archive thumbnails require explicit display-only opt-in.
    #[serde(default)]
    pub allow_native_thumbnail: bool,
}

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
pub struct GelSourceRecord {
    pub image_id: String,
    pub source_name: String,
    pub source_sha256: String,
    pub original_base64: String,
    pub extracted_image_sha256: String,
    pub extraction: String,
    pub warnings: Vec<String>,
    /// Safe replay of the exact extraction, with basename rather than local path.
    pub request: GelSourceImportRequest,
}

/// Adapter-facing identity only. Native archives can contain operational logs;
/// their original bytes stay in local project storage, not operation feedback.
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
pub struct GelSourceSummary {
    pub image_id: String,
    pub source_name: String,
    pub source_sha256: String,
    pub extracted_image_sha256: String,
    pub extraction: String,
    pub warnings: Vec<String>,
}

impl From<&GelSourceRecord> for GelSourceSummary {
    fn from(source: &GelSourceRecord) -> Self {
        Self {
            image_id: source.image_id.clone(),
            source_name: source.source_name.clone(),
            source_sha256: source.source_sha256.clone(),
            extracted_image_sha256: source.extracted_image_sha256.clone(),
            extraction: source.extraction.clone(),
            warnings: source.warnings.clone(),
        }
    }
}

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
#[serde(deny_unknown_fields)]
pub struct GelDetectRequest {
    pub image_id: String,
    pub lanes: Vec<GelImageLane>,
    pub migration: GelMigrationDirection,
    pub dark_bands: bool,
    pub minimum_prominence: f64,
    pub minimum_distance_px: u32,
    pub smoothing_radius_px: u32,
}

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
pub struct GelBandSuggestions {
    pub image_sha256: String,
    pub request: GelDetectRequest,
    pub bands: Vec<GelSampleBand>,
    pub warnings: Vec<String>,
}

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
#[serde(deny_unknown_fields)]
pub struct GelExtractProfileRequest {
    pub profile_id: String,
    pub analysis_id: String,
    pub lane_id: String,
    pub run_id: String,
    pub dark_bands: bool,
    /// Explicit background in normalized 0..1 pixel intensity. No fitted or
    /// cross-lane background is silently estimated.
    pub background: f64,
    /// User-reviewed size intervals (e.g. alignment markers). Any intersecting
    /// pixel bin is unavailable rather than split into an invented intensity.
    #[serde(default)]
    pub masked_ranges: Vec<[f64; 2]>,
}

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
#[serde(deny_unknown_fields)]
pub struct GelCompareRequest {
    pub comparison_id: String,
    pub left_profile: String,
    pub right_profile: String,
    pub range: [f64; 2],
    #[serde(default)]
    pub short_threshold: Option<f64>,
}

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
pub struct GelDistributionSummary {
    pub total_weight: f64,
    pub q10: f64,
    pub median: f64,
    pub q90: f64,
    pub fraction_below_threshold: Option<f64>,
}

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
pub struct GelQqPoint {
    pub probability: f64,
    pub left: f64,
    pub right: f64,
}

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
pub struct GelComparisonReport {
    pub schema: String,
    pub request: GelCompareRequest,
    /// Full frozen inputs: no re-query, re-extraction or current form settings
    /// participate when this immutable result is inspected or exported.
    pub left: GelProfile,
    pub right: GelProfile,
    pub common_edges: Vec<f64>,
    pub left_weights: Vec<Option<f64>>,
    pub right_weights: Vec<Option<f64>>,
    pub excluded_intervals: Vec<[f64; 2]>,
    pub left_summary: GelDistributionSummary,
    pub right_summary: GelDistributionSummary,
    pub qq: Vec<GelQqPoint>,
    pub ks_distance: f64,
    pub wasserstein_distance: f64,
    /// No inferential p-value: trace bins/pixels are correlated measurements,
    /// not independent sampled fragment sizes or biological replicates.
    pub p_value: Option<f64>,
    pub warnings: Vec<String>,
}

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
#[serde(deny_unknown_fields)]
pub struct GelAlgebraTerm {
    pub profile_id: String,
    pub coefficient: f64,
}

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
#[serde(deny_unknown_fields)]
pub struct GelAlgebraRequest {
    pub output_id: String,
    pub terms: Vec<GelAlgebraTerm>,
    pub range: [f64; 2],
}

#[derive(Debug, Clone, Copy, Serialize, Deserialize, PartialEq, Eq)]
#[serde(rename_all = "snake_case")]
pub enum GelAlgebraVerification {
    Verified,
    Mismatch,
    Unverifiable,
}

/// Display transform only: each density is the recorded bin integral divided
/// by its size width. It is neither normalized nor calibrated fluorescence.
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
pub struct GelAlgebraSeries {
    pub size_kind: GelSizeKind,
    pub signal_units: String,
    pub measure: String,
    pub edges: Vec<f64>,
    pub weights: Vec<Option<f64>>,
    pub density: Vec<Option<f64>>,
}

/// One nonzero request term projected onto `series.edges`, in request order.
/// Weights are coefficient times rebinned operand integral; density divides
/// that result by bin width. Each term retains its own missing coverage, even
/// where another term makes the residual unavailable. Repeated IDs stay distinct.
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
#[serde(deny_unknown_fields)]
pub struct GelAlgebraContribution {
    pub term_index: usize,
    pub profile_id: String,
    pub coefficient: f64,
    pub weights: Vec<Option<f64>>,
    pub density: Vec<Option<f64>>,
}

/// Read-only, bounded replay of algebra v1. Verification means internal
/// consistency, not source authenticity or scientific acceptance. Failed
/// verification carries no plot series. Inputs include nested dependencies;
/// legacy residual run/lane fields are not identities of the derived result.
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
pub struct GelAlgebraPresentation {
    pub schema: String,
    pub profile_id: String,
    pub label: String,
    pub verification: GelAlgebraVerification,
    pub reason: Option<String>,
    pub input_digest: Option<String>,
    pub request: Option<GelAlgebraRequest>,
    pub series: Option<GelAlgebraSeries>,
    /// Absent in legacy presentations; re-inspect to obtain overlays. Zero
    /// coefficients remain in the request/provenance but have no plotted term.
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub contributions: Vec<GelAlgebraContribution>,
    pub inputs: Vec<GelProfile>,
    pub warnings: Vec<String>,
}

/// Presentation only: neither mode changes frozen integrals or statistics.
#[derive(Debug, Clone, Copy, Serialize, Deserialize, PartialEq, Eq)]
#[serde(rename_all = "snake_case")]
pub enum GelComparisonDisplay {
    BinIntegral,
    SignalDensityPerSize,
}

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
#[serde(deny_unknown_fields)]
pub struct GelWorkbenchExportRequest {
    pub path: String,
    /// Explicit ordering gives the PDF its before/after pair. No pairing is
    /// inferred from adjacent wells or repeated lane names.
    pub analysis_ids: Vec<String>,
    #[serde(default)]
    pub comparison_ids: Vec<String>,
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub algebra_ids: Vec<String>,
    /// PDF comparisons only. Omission preserves the legacy integral display.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub comparison_display: Option<GelComparisonDisplay>,
    pub title: String,
}

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
#[serde(tag = "action", rename_all = "snake_case", deny_unknown_fields)]
pub enum GelWorkbenchRequest {
    ImportSource(GelSourceImportRequest),
    DetectBands(GelDetectRequest),
    ImportProfile { profile: GelProfile },
    ImportProfileFile(GelProfileFileRequest),
    ExtractProfile(GelExtractProfileRequest),
    Compare(GelCompareRequest),
    Algebra(GelAlgebraRequest),
    // An empty struct (not a unit variant) makes Serde enforce unknown fields
    // on this internally tagged action too. Its JSON remains {"action":"inspect"}.
    Inspect {},
    InspectProfile { profile_id: String },
    InspectComparison { comparison_id: String },
    InspectAlgebra { profile_id: String },
    Export(GelWorkbenchExportRequest),
}

impl GelWorkbenchRequest {
    pub fn mutates(&self) -> bool {
        !matches!(
            self,
            Self::Inspect {}
                | Self::InspectProfile { .. }
                | Self::InspectComparison { .. }
                | Self::InspectAlgebra { .. }
                | Self::Export(_)
                | Self::DetectBands(_)
        )
    }
}

#[derive(Debug, Clone, Default, Serialize, Deserialize)]
pub struct GelWorkbenchResult {
    #[serde(skip_serializing_if = "Option::is_none")]
    pub source: Option<GelSourceSummary>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub suggestions: Option<GelBandSuggestions>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub profile: Option<GelProfile>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub comparison: Option<GelComparisonReport>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub algebra: Option<GelAlgebraPresentation>,
    pub profile_ids: Vec<String>,
    pub comparison_ids: Vec<String>,
    pub exported_sha256: Option<String>,
}

#[derive(Debug, Clone, Default, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct GelWorkbenchStore {
    #[serde(default)]
    pub numeric_sources: BTreeMap<String, Arc<GelNumericSource>>,
    #[serde(default)]
    pub sources: BTreeMap<String, Arc<GelSourceRecord>>,
    #[serde(default)]
    pub profiles: BTreeMap<String, Arc<GelProfile>>,
    #[serde(default)]
    pub comparisons: BTreeMap<String, Arc<GelComparisonReport>>,
}

impl GelWorkbenchStore {
    pub fn is_empty(&self) -> bool {
        self.sources.is_empty()
            && self.numeric_sources.is_empty()
            && self.profiles.is_empty()
            && self.comparisons.is_empty()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use serde_json::json;

    #[test]
    fn legacy_gel_store_keeps_its_bytes_without_a_workbench() {
        let legacy = json!({"images": {}, "analyses": {}});
        let store: super::super::gel_image::GelImageStore =
            serde_json::from_value(legacy.clone()).unwrap();
        assert!(store.workbench.is_empty());
        assert_eq!(serde_json::to_value(store).unwrap(), legacy);
    }

    #[test]
    fn shared_workbench_route_keeps_canonical_parity_projection_current() {
        let generated = crate::render_gui_cli_mcp_parity_matrix_markdown();
        assert!(generated.contains("| gel-image workbench |"));
        assert!(generated.contains("| GelWorkbench |"));
        // This generated file already has scoped LF/CRLF checkout coverage in
        // scripts.test_tutorial_checkouts; compare the canonical producer too.
        assert_eq!(
            generated,
            include_str!("../../../docs/gui_cli_mcp_parity.md")
        );
    }

    #[test]
    fn requests_preserve_explicit_actions_defaults_and_mutation_classification() {
        for (value, mutates) in [
            (json!({"action":"inspect"}), false),
            (json!({"action":"inspect_profile", "profile_id":"a"}), false),
            (json!({"action":"inspect_algebra", "profile_id":"a"}), false),
            (
                json!({"action":"inspect_comparison", "comparison_id":"a"}),
                false,
            ),
            (
                json!({"action":"import_source", "image_id":"pdf", "path":"run.pdf"}),
                true,
            ),
            (
                json!({"action":"extract_profile", "profile_id":"p", "analysis_id":"a",
                "lane_id":"B1", "run_id":"run", "dark_bands":true, "background":0}),
                true,
            ),
            (
                json!({"action":"export", "path":"report.pdf", "title":"Saved report",
                "analysis_ids":["before","after"]}),
                false,
            ),
        ] {
            let request: GelWorkbenchRequest = serde_json::from_value(value.clone()).unwrap();
            assert_eq!(request.mutates(), mutates, "{value}");
            let roundtrip = serde_json::to_value(&request).unwrap();
            assert_eq!(
                serde_json::from_value::<GelWorkbenchRequest>(roundtrip).unwrap(),
                request
            );
            let mut unknown = value;
            unknown["unreviewed_option"] = true.into();
            let result = serde_json::from_value::<GelWorkbenchRequest>(unknown.clone());
            assert!(
                result.is_err(),
                "Unknown fields accepted: {unknown}: {result:?}"
            );
        }
        let request: GelWorkbenchRequest = serde_json::from_value(json!({
            "action":"import_source", "image_id":"native", "path":"run.D5000"
        }))
        .unwrap();
        let GelWorkbenchRequest::ImportSource(request) = request else {
            panic!("source request")
        };
        assert!(!request.allow_native_thumbnail);
        assert_eq!(request.pdf_object, None);
    }

    #[test]
    fn profile_roundtrip_distinguishes_zero_missing_and_signed_residuals() {
        let profile = GelProfile {
            id: "residual".into(),
            label: "Synthetic signed residual".into(),
            run_id: "synthetic".into(),
            lane_id: "A1".into(),
            source_sha256: "a".repeat(64),
            evidence_kind: GelEvidenceKind::SignedResidual,
            size_kind: GelSizeKind::LinearDnaBp,
            signal_units: "synthetic_bin_area".into(),
            measure: "explicit_integrated_bins".into(),
            edges: vec![100.0, 200.0, 300.0, 400.0],
            weights: vec![Some(0.0), None, Some(-2.0)],
            processing: vec![],
            warnings: vec![],
        };
        let json = serde_json::to_value(&profile).unwrap();
        assert_eq!(json["weights"], json!([0.0, null, -2.0]));
        assert_eq!(serde_json::from_value::<GelProfile>(json).unwrap(), profile);
    }

    #[test]
    fn algebra_contributions_are_additive_to_legacy_presentations_and_strict() {
        let legacy = json!({
            "schema":"gentle.gel_algebra_presentation.v1", "profile_id":"r", "label":"Synthetic",
            "verification":"unverifiable", "reason":null, "input_digest":null,
            "request":null, "series":null, "inputs":[], "warnings":[]
        });
        let report: GelAlgebraPresentation = serde_json::from_value(legacy.clone()).unwrap();
        assert!(report.contributions.is_empty());
        assert_eq!(serde_json::to_value(report).unwrap(), legacy);
        let mut extended = legacy;
        extended["contributions"] = json!([{
            "term_index":1, "profile_id":"synthetic", "coefficient":-0.5,
            "weights":[-2.0,null], "density":[-0.02,null]
        }]);
        let report: GelAlgebraPresentation = serde_json::from_value(extended.clone()).unwrap();
        assert_eq!(serde_json::to_value(report).unwrap(), extended);
        extended["contributions"][0]["fitted_scale"] = json!(2);
        assert!(serde_json::from_value::<GelAlgebraPresentation>(extended).is_err());
    }

    #[test]
    fn algebra_export_defaults_preserve_legacy_requests_and_result_json() {
        let legacy = json!({"action":"export","path":"saved.json","title":"Saved","analysis_ids":[],"comparison_ids":["c"]});
        let request: GelWorkbenchRequest = serde_json::from_value(legacy.clone()).unwrap();
        assert_eq!(serde_json::to_value(&request).unwrap(), legacy);
        let mut extended = legacy;
        extended["algebra_ids"] = json!(["r"]);
        let request: GelWorkbenchRequest = serde_json::from_value(extended.clone()).unwrap();
        assert!(!request.mutates());
        assert_eq!(serde_json::to_value(request).unwrap(), extended);
        let mut unknown = extended;
        unknown["allow_unverified"] = true.into();
        assert!(serde_json::from_value::<GelWorkbenchRequest>(unknown).is_err());
        let result: GelWorkbenchResult = serde_json::from_value(
            json!({"profile_ids":[],"comparison_ids":[],"exported_sha256":null}),
        )
        .unwrap();
        assert!(result.algebra.is_none());
        assert!(
            serde_json::to_value(result)
                .unwrap()
                .get("algebra")
                .is_none()
        );
    }

    #[test]
    fn comparison_display_is_explicit_optional_and_strict() {
        let legacy = json!({"action":"export","path":"saved.pdf","title":"Saved","analysis_ids":[],"comparison_ids":["c"]});
        let request: GelWorkbenchRequest = serde_json::from_value(legacy.clone()).unwrap();
        assert_eq!(serde_json::to_value(request).unwrap(), legacy);
        for mode in ["bin_integral", "signal_density_per_size"] {
            let mut explicit = legacy.clone();
            explicit["comparison_display"] = mode.into();
            let request: GelWorkbenchRequest = serde_json::from_value(explicit.clone()).unwrap();
            assert!(!request.mutates());
            assert_eq!(serde_json::to_value(request).unwrap(), explicit);
        }
        let mut bad = legacy.clone();
        bad["comparison_display"] = "normalize_each_maximum".into();
        assert!(serde_json::from_value::<GelWorkbenchRequest>(bad).is_err());
        let mut bad = legacy;
        bad["display"] = "bin_integral".into();
        assert!(serde_json::from_value::<GelWorkbenchRequest>(bad).is_err());
    }
}
