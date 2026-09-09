//! Measured-gel image contracts, distinct from simulated gel migration.
//!
//! All coordinates are zero-based pixel centers in the original decoded image,
//! with x increasing right and y down. No EXIF orientation or display transform
//! is implicitly applied. Reported protein sizes are apparent masses, not identities.

use serde::{Deserialize, Serialize};
use std::{collections::BTreeMap, sync::Arc};

pub const GEL_IMAGE_SCHEMA: &str = "gentle.gel_image.v1";
pub const GEL_IMAGE_ANALYSIS_SCHEMA: &str = "gentle.gel_image_analysis.v1";
pub const GEL_CALIBRATION_METHOD: &str = "piecewise_log10_size_v1";

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
#[serde(deny_unknown_fields)]
pub struct GelImageImportRequest {
    pub image_id: String,
    pub path: String,
    #[serde(default)]
    pub tiff_page: Option<u32>,
}

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
#[serde(deny_unknown_fields)]
pub struct GelImageExportRequest {
    pub report_id: String,
    pub path: String,
    pub format: GelImageExportFormat,
}

#[derive(Debug, Clone, Copy, Serialize, Deserialize, PartialEq, Eq)]
#[serde(rename_all = "snake_case")]
pub enum GelSizeKind {
    LinearDnaBp,
    SdsProteinKda,
}

impl GelSizeKind {
    pub fn unit(self) -> &'static str {
        match self {
            Self::LinearDnaBp => "bp",
            Self::SdsProteinKda => "kDa",
        }
    }
}

#[derive(Debug, Clone, Copy, Serialize, Deserialize, PartialEq, Eq)]
#[serde(rename_all = "snake_case")]
pub enum GelMigrationDirection {
    Down,
    Up,
    Right,
    Left,
}

impl GelMigrationDirection {
    /// Signed coordinate increasing along migration, independent of screen orientation.
    pub fn coordinate(self, point: GelImagePoint) -> f64 {
        match self {
            Self::Down => point.y,
            Self::Up => -point.y,
            Self::Right => point.x,
            Self::Left => -point.x,
        }
    }
}

#[derive(Debug, Clone, Copy, Serialize, Deserialize, PartialEq)]
#[serde(deny_unknown_fields)]
pub struct GelImagePoint {
    pub x: f64,
    pub y: f64,
}

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
#[serde(deny_unknown_fields)]
pub struct GelImageLane {
    pub id: String,
    pub label: String,
    /// Inclusive pixel-center bounds in the original decoded image.
    pub min: GelImagePoint,
    pub max: GelImagePoint,
}

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
#[serde(deny_unknown_fields)]
pub struct GelLadderBand {
    pub id: String,
    pub center: GelImagePoint,
    /// bp or kDa, as specified by the request's size_kind. Must be positive.
    pub size: f64,
}

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
#[serde(deny_unknown_fields)]
pub struct GelImageLadder {
    pub lane_id: String,
    pub label: String,
    /// Exact vendor document/catalog or explicit custom/synthetic origin.
    pub source: String,
    #[serde(default)]
    pub gel_system: Option<String>,
    #[serde(default)]
    pub prestained: bool,
    /// Only experimentally visible, user-confirmed assignments, never guessed missing bands.
    pub bands: Vec<GelLadderBand>,
}

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
#[serde(deny_unknown_fields)]
pub struct GelSampleBand {
    pub id: String,
    pub lane_id: String,
    pub label: String,
    pub center: GelImagePoint,
    /// Optional user-supplied localization half-width along migration, not a confidence interval.
    #[serde(default)]
    pub position_half_width_px: Option<f64>,
}

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
#[serde(deny_unknown_fields)]
pub struct GelImageAnalysisRequest {
    pub report_id: String,
    pub image_id: String,
    /// Binds manually supplied coordinates to the exact imported original bytes.
    pub image_sha256: String,
    pub size_kind: GelSizeKind,
    pub migration: GelMigrationDirection,
    pub lanes: Vec<GelImageLane>,
    pub ladder: GelImageLadder,
    pub sample_bands: Vec<GelSampleBand>,
}

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq)]
pub struct GelImageDescriptor {
    pub schema: String,
    pub image_id: String,
    /// File basename only. Absolute workstation paths are not included in portable reports.
    pub source_name: String,
    pub sha256: String,
    pub byte_count: usize,
    pub width: u32,
    pub height: u32,
    pub format: String,
    pub color_type: String,
    #[serde(default)]
    pub tiff_page: Option<u32>,
    pub warnings: Vec<String>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct GelImageRecord {
    pub descriptor: GelImageDescriptor,
    /// Immutable original bytes, retaining native bit depth and acquisition metadata.
    pub original_base64: String,
    /// Bounded, display-only PNG. Never used to determine calibration coordinates.
    pub preview_png_base64: String,
}

#[derive(Debug, Clone, Copy, Serialize, Deserialize, PartialEq, Eq)]
#[serde(rename_all = "snake_case")]
pub enum GelBandSizingStatus {
    Interpolated,
    OutsideCalibratedRange,
}

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
pub struct GelBandSizeEstimate {
    pub band_id: String,
    pub lane_id: String,
    pub label: String,
    pub center: GelImagePoint,
    pub status: GelBandSizingStatus,
    pub estimated_size: Option<f64>,
    pub reference_band_ids: Vec<String>,
    /// Min/max size from the specified localization half-width only. Not total uncertainty.
    pub localization_size_bounds: Option<[f64; 2]>,
    pub warnings: Vec<String>,
}

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
pub struct GelImageAnalysisReport {
    pub schema: String,
    pub algorithm: String,
    pub gentle_version: String,
    pub image: GelImageDescriptor,
    pub request: GelImageAnalysisRequest,
    /// Confirmed reference bands sorted along migration, larger sizes first.
    pub calibration_bands: Vec<GelLadderBand>,
    pub estimates: Vec<GelBandSizeEstimate>,
    pub warnings: Vec<String>,
}

#[derive(Debug, Clone, Copy, Serialize, Deserialize, PartialEq, Eq)]
#[serde(rename_all = "snake_case")]
pub enum GelImageExportFormat {
    Json,
    Tsv,
    Svg,
}

/// Immutable images/reports share allocations across project undo checkpoints.
/// Arc is an in-memory optimization only; project JSON remains self-contained.
#[derive(Debug, Clone, Default, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct GelImageStore {
    #[serde(default)]
    pub images: BTreeMap<String, Arc<GelImageRecord>>,
    #[serde(default)]
    pub analyses: BTreeMap<String, Arc<GelImageAnalysisReport>>,
}

impl GelImageStore {
    pub fn is_empty(&self) -> bool {
        self.images.is_empty() && self.analyses.is_empty()
    }
}
