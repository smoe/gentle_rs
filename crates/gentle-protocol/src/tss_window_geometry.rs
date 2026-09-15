//! Stateless fixed TSS-window geometry for similarity queries, not reporter insert selection.

use crate::tss_profiles::{TssGeometry, TssStrand};
use serde::{Deserialize, Serialize};

pub const TSS_WINDOW_GEOMETRY_REQUEST_SCHEMA: &str = "gentle.tss_window_geometry_request.v1";
pub const TSS_WINDOW_GEOMETRY_REPORT_SCHEMA: &str = "gentle.tss_window_geometry.v1";

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq)]
#[serde(deny_unknown_fields)]
pub struct TssWindowAnchor {
    pub anchor_id: String,
    /// Exact prepared window; requested geometry must fit inside this unclipped source.
    pub source: TssGeometry,
}

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq)]
#[serde(deny_unknown_fields)]
pub struct TssWindowFeature {
    pub feature_id: String,
    pub start_1based: u64,
    pub end_1based: u64,
}

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq)]
#[serde(deny_unknown_fields)]
pub struct TssWindowGeometryGroup {
    pub group_id: String,
    pub chromosome: String,
    pub strand: TssStrand,
    pub anchors: Vec<TssWindowAnchor>,
    pub features: Vec<TssWindowFeature>,
}

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq)]
#[serde(deny_unknown_fields)]
pub struct TssWindowGeometryRequest {
    pub schema: String,
    /// Caller-bound reference label; this operation does not authenticate a reference.
    pub assembly: String,
    pub upstream_bp: usize,
    pub downstream_bp: usize,
    pub groups: Vec<TssWindowGeometryGroup>,
}

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq)]
pub struct TssWindowGeometryWindow {
    pub anchor_id: String,
    pub start_1based: u64,
    pub end_1based: u64,
}

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq)]
pub struct TssConnectedWindowStretch {
    pub start_1based: u64,
    pub end_1based: u64,
    pub anchor_ids: Vec<String>,
}

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq)]
pub struct TssWindowFeatureIntersection {
    pub stretch_index_1based: usize,
    pub feature_id: String,
    pub start_1based: u64,
    pub end_1based: u64,
    /// Prepared windows containing the entire intersection, for bound sequence extraction.
    pub containing_anchor_ids: Vec<String>,
}

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq)]
pub struct TssWindowGeometryGroupResult {
    pub group_id: String,
    pub windows: Vec<TssWindowGeometryWindow>,
    pub stretches: Vec<TssConnectedWindowStretch>,
    pub intersections: Vec<TssWindowFeatureIntersection>,
}

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq)]
pub struct TssWindowGeometryReport {
    pub schema: String,
    pub request: TssWindowGeometryRequest,
    pub request_sha256: String,
    pub groups: Vec<TssWindowGeometryGroupResult>,
}
