//! Portable, accession-pinned TSS profile inputs, computed reports and exports.
//!
//! Window coordinates are transcript-oriented base indices, never motif centers.
//! Computed scores are retained before display clipping; unavailable windows are null.

use serde::{Deserialize, Serialize};
use std::collections::BTreeMap;

pub const PANEL_SCHEMA: &str = "gentle.jaspar_target_panel.v1";
pub const BUNDLE_SCHEMA: &str = "gentle.tss_fasta_bundle.v1";
pub const SELECTION_SCHEMA: &str = "gentle.tss_profile_selection.v1";
pub const REPORT_SCHEMA: &str = "gentle.tss_tfbs_profiles.v1";
pub const RECEIPT_SCHEMA: &str = "gentle.tss_tfbs_profile_receipt.v1";
pub const NON_CLAIMS: &str = "JASPAR tracks are sequence-model predictions, not measured binding, affinity, occupancy, cofactor interaction, promoter activity or functional regulation. Cross-matrix magnitudes are not comparable without a documented calibration. Annotated transcript starts are TSS candidates, not experimentally established initiation sites. Correlation compares model outputs, not biological correctness or co-regulation. A factor also appearing as a target gene does not establish autoregulation.";

#[derive(Debug, Clone, Copy, Serialize, Deserialize, PartialEq, Eq)]
pub enum TssStrand {
    #[serde(rename = "+")]
    Plus,
    #[serde(rename = "-")]
    Minus,
}

impl TssStrand {
    pub fn as_str(self) -> &'static str {
        match self {
            Self::Plus => "+",
            Self::Minus => "-",
        }
    }
    pub fn opposite(self) -> Self {
        match self {
            Self::Plus => Self::Minus,
            Self::Minus => Self::Plus,
        }
    }
}

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq)]
#[serde(deny_unknown_fields)]
pub struct TssReference {
    pub genome_id: String,
    pub assembly: String,
    /// Only a separately declared release, never guessed from a genome label.
    #[serde(default)]
    pub annotation_release: Option<String>,
}

/// Exact input-manifest identity, distinct from the profile producer revision.
#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct TssBundleSource {
    pub schema: String,
    pub manifest_sha256: String,
    pub source_revision: Option<String>,
    pub dataset_id: Option<String>,
    pub producer_sha256: Option<String>,
}

/// A source-bound descriptive selection criterion, not a new activity inference.
#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct TssSelectionEvidence {
    pub label: String,
    pub legend: String,
    pub criterion: String,
    pub factor: Option<String>,
}

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq)]
#[serde(deny_unknown_fields)]
pub struct TssGeometry {
    pub chromosome: String,
    pub strand: TssStrand,
    pub tss_1based: u64,
    pub start_1based: u64,
    pub end_1based: u64,
    pub upstream_bp: usize,
    pub downstream_bp: usize,
}

impl TssGeometry {
    /// Number of bases including the TSS base; None on arithmetic overflow.
    pub fn length(&self) -> Option<usize> {
        self.upstream_bp
            .checked_add(1)?
            .checked_add(self.downstream_bp)
    }

    /// Absolute coordinate of a transcript-oriented base index (not a motif end).
    pub fn genomic_at(&self, index: usize) -> Option<u64> {
        if index >= self.length()? {
            return None;
        }
        let index = u64::try_from(index).ok()?;
        match self.strand {
            TssStrand::Plus => self.start_1based.checked_add(index),
            TssStrand::Minus => self.end_1based.checked_sub(index),
        }
    }

    pub fn relative_at(&self, index: usize) -> Option<i64> {
        if index >= self.length()? {
            return None;
        }
        i64::try_from(index)
            .ok()?
            .checked_sub(i64::try_from(self.upstream_bp).ok()?)
    }
}

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq)]
#[serde(deny_unknown_fields)]
pub struct TssRecord {
    pub promoter_id: String,
    pub gene_id: String,
    pub gene_symbol: String,
    pub geometry: TssGeometry,
    pub transcripts: Vec<String>,
    pub sequence_sha256: String,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct TssBundleManifest {
    pub schema: String,
    pub reference: TssReference,
    /// Relative FASTA filename -> SHA-256 of the exact file bytes.
    pub fasta_files: BTreeMap<String, String>,
    pub records: Vec<TssRecord>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct TssSelection {
    pub schema: String,
    pub reference: TssReference,
    pub selected: Vec<TssSelectedRecord>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct TssSelectedRecord {
    pub promoter_id: String,
    pub gene_id: String,
}

#[derive(Debug, Clone, Copy, Serialize, Deserialize, PartialEq, Eq, Default)]
#[serde(rename_all = "snake_case")]
pub enum TssScaleMode {
    #[default]
    Independent,
    Shared,
}

#[derive(Debug, Clone, Copy, Serialize, Deserialize, PartialEq, Eq)]
#[serde(rename_all = "snake_case")]
pub enum TssCalibrationState {
    MatrixSpecific,
    CrossSourceCalibrated,
}

#[derive(Debug, Clone, Copy, Serialize, Deserialize, PartialEq, Eq, Default)]
#[serde(rename_all = "snake_case")]
pub enum TssStrandPolicy {
    #[default]
    Both,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct JasparTargetPanel {
    pub schema: String,
    pub panel_id: String,
    pub label: String,
    /// Exact existing score-kind spelling; validated against the shared scorer.
    pub score_kind: String,
    pub clip_negative: bool,
    #[serde(default)]
    pub scale_mode: TssScaleMode,
    #[serde(default)]
    pub strand_policy: TssStrandPolicy,
    pub calibration_state: TssCalibrationState,
    pub calibration_statement: String,
    #[serde(default)]
    pub calibration_id: Option<String>,
    #[serde(default)]
    pub calibration_sha256: Option<String>,
    pub top_hit_count: usize,
    /// One entry per exact matrix, not one entry per factor name.
    pub factors: Vec<JasparPanelTrack>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct JasparPanelTrack {
    pub source_id: String,
    pub factor_id: String,
    pub label: String,
    pub display_order: usize,
    #[serde(default)]
    pub color_hint: Option<String>,
    /// Optional repeated policy must equal the panel's score kind in v1.
    #[serde(default)]
    pub score_kind: Option<String>,
    /// Preserve original track identity and labels when normalizing tracks[].
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub track_id: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub provider_kind: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub factor_label: Option<String>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct TssInputBinding {
    pub role: String,
    /// Portable filename, never an absolute host path.
    pub name: String,
    pub sha256: String,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ResolvedTssMatrix {
    pub specification: JasparPanelTrack,
    pub version: String,
    pub consensus: String,
    /// Actual A/C/G/T counts used; also permits report-only logo rendering.
    pub matrix_counts: Vec<[f64; 4]>,
    pub matrix_sha256: String,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct TssPanelResolution {
    pub panel: JasparTargetPanel,
    pub panel_sha256: String,
    pub registry_sources: Vec<TssInputBinding>,
    pub registry_source_url: Option<String>,
    pub matrices: Vec<ResolvedTssMatrix>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct TssPeak {
    pub local_start_0based: usize,
    pub score: f64,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct TssProfileTrack {
    pub accession: String,
    pub motif_length_bp: usize,
    pub forward_scores: Vec<Option<f64>>,
    pub reverse_scores: Vec<Option<f64>>,
    pub forward_maximum: Option<TssPeak>,
    pub reverse_maximum: Option<TssPeak>,
    pub forward_peaks: Vec<TssPeak>,
    pub reverse_peaks: Vec<TssPeak>,
    /// Shared scorer normalization metadata; schema is unchanged from score tracks.
    pub normalization_reference: serde_json::Value,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct TssMatrixComparison {
    pub factor_id: String,
    pub left_accession: String,
    pub right_accession: String,
    pub local_strand: TssStrand,
    pub paired_window_count: usize,
    pub excluded_window_count: usize,
    pub pearson: Option<f64>,
    pub spearman: Option<f64>,
    pub undefined_reason: Option<String>,
    pub method: String,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct TssProfileWindow {
    pub record: TssRecord,
    pub selected: bool,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub selection_evidence: Option<TssSelectionEvidence>,
    pub tracks: Vec<TssProfileTrack>,
    pub comparisons: Vec<TssMatrixComparison>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct TssProfileReport {
    pub schema: String,
    pub reference: TssReference,
    pub panel_resolution: TssPanelResolution,
    pub inputs: Vec<TssInputBinding>,
    pub windows: Vec<TssProfileWindow>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub source: Option<TssBundleSource>,
    pub producer_revision: String,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub producer_executable_sha256: Option<String>,
    pub lockfile_sha256: String,
    pub score_policy: BTreeMap<String, String>,
    pub verification: String,
    pub warnings: Vec<String>,
    pub non_claims: String,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct ComputeTssProfilesRequest {
    pub manifest: String,
    pub panel: String,
    #[serde(default)]
    pub fasta: Vec<String>,
    #[serde(default)]
    pub selection: Option<String>,
    pub expected_genome_id: String,
    #[serde(default)]
    pub expected_assembly: Option<String>,
    #[serde(default)]
    pub expected_annotation_release: Option<String>,
    #[serde(default)]
    pub expected_dataset_id: Option<String>,
}

#[derive(Debug, Clone, Copy, Serialize, Deserialize, PartialEq, Eq)]
#[serde(rename_all = "snake_case")]
pub enum TssExportFormat {
    Svg,
    Png,
    Pdf,
}

fn one_panel() -> usize {
    1
}
fn svg_format() -> Vec<TssExportFormat> {
    vec![TssExportFormat::Svg]
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct TssProfileRenderOptions {
    #[serde(default)]
    pub scale_mode: Option<TssScaleMode>,
    #[serde(default = "one_panel")]
    pub panels_per_page: usize,
}

impl Default for TssProfileRenderOptions {
    fn default() -> Self {
        Self {
            scale_mode: None,
            panels_per_page: 1,
        }
    }
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct ExportTssProfilesRequest {
    pub output_dir: String,
    #[serde(default)]
    pub rendering: TssProfileRenderOptions,
    #[serde(default = "svg_format")]
    pub formats: Vec<TssExportFormat>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct TssProfileReceipt {
    pub schema: String,
    pub report_sha256: String,
    pub input_manifest_sha256: String,
    pub source_revision: Option<String>,
    pub producer_revision: String,
    pub exporter_revision: String,
    pub executable_sha256: String,
    pub lockfile_sha256: String,
    pub inputs: Vec<TssInputBinding>,
    pub rendering: TssProfileRenderOptions,
    pub renderer: String,
    pub render_metadata: Vec<serde_json::Value>,
    pub outputs: BTreeMap<String, String>,
    pub tss_count: usize,
    pub page_count: usize,
    pub non_claims: String,
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn tss_profiles_json_preserves_score_bits_for_report_only_replay() {
        let scores = vec![
            Some(-109.91133953680091_f64),
            None,
            Some(0.0),
            Some(f64::MIN_POSITIVE),
        ];
        let bytes = serde_json::to_vec(&scores).unwrap();
        let restored: Vec<Option<f64>> = serde_json::from_slice(&bytes).unwrap();
        let bits = |values: Vec<Option<f64>>| {
            values
                .into_iter()
                .map(|v| v.map(f64::to_bits))
                .collect::<Vec<_>>()
        };
        assert_eq!(bits(scores), bits(restored));
    }

    #[test]
    fn tss_profiles_geometry_roundtrip_keeps_decreasing_genomic_axis() {
        let geometry = TssGeometry {
            chromosome: "synthetic".into(),
            strand: TssStrand::Minus,
            tss_1based: 1000,
            start_1based: 800,
            end_1based: 1500,
            upstream_bp: 500,
            downstream_bp: 200,
        };
        let restored: TssGeometry =
            serde_json::from_value(serde_json::to_value(&geometry).unwrap()).unwrap();
        assert_eq!(restored, geometry);
        for (i, relative, genomic) in [(0, -500, 1500), (500, 0, 1000), (700, 200, 800)] {
            assert_eq!(restored.relative_at(i), Some(relative));
            assert_eq!(restored.genomic_at(i), Some(genomic));
        }
        assert_eq!(restored.genomic_at(701), None);
    }

    #[test]
    fn tss_profiles_capabilities_classify_export_effects_conservatively() {
        for name in ["ComputeTssTfbsProfiles", "ExportTssTfbsProfiles"] {
            assert_eq!(
                crate::infer_engine_operation_mutation(name),
                crate::CapabilityMutation::External
            );
        }
    }

    #[test]
    fn tss_profiles_requests_fail_on_unknown_settings() {
        let value = serde_json::json!({"manifest":"m.json","panel":"p.json","expected_genome_id":"g","ignore_errors":true});
        assert!(serde_json::from_value::<ComputeTssProfilesRequest>(value).is_err());
        assert!(
            serde_json::from_str::<ExportTssProfilesRequest>(
                r#"{"output_dir":"out","overwrite":true}"#
            )
            .is_err()
        );
    }
}
