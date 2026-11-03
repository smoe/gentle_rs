//! Simulator-neutral rack geometry descriptions, never robot execution authority.

use crate::{RackPhysicalTemplateSpec, RackProfileSnapshot, labware::TubeProfile};
use serde::{Deserialize, Serialize};

/// Independent of the historical rack simulation export and any simulator schema.
pub const RACK_EQUIPMENT_SCHEMA: &str = "gentle.rack_equipment_description.v1";
/// Bump when projection, frames, units or interpretation changes.
pub const RACK_EQUIPMENT_PROJECTION: &str = "gentle.rack_equipment_projection.v1";

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
/// Exact embedded source values; profiles are drawings, not fit approvals.
pub struct RackEquipmentSource {
    pub projection: String,
    #[serde(serialize_with = "serialize_template")]
    pub physical_template: RackPhysicalTemplateSpec,
    pub rack_profile: RackProfileSnapshot,
    pub tube_profiles: Vec<TubeProfile>,
}

// Keep typed API serialization identical to shell/MCP JSON-value projection.
// serde_json otherwise writes f32 text with less precision than to_value uses.
fn serialize_template<S: serde::Serializer>(
    template: &RackPhysicalTemplateSpec,
    serializer: S,
) -> Result<S::Ok, S::Error> {
    let value = serde_json::to_value(template).map_err(serde::ser::Error::custom)?;
    value.serialize(serializer)
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
/// Rack-local frame only: translation in metres, quaternion in x/y/z/w order.
pub struct RackEquipmentFrame {
    pub id: String,
    pub parent: Option<String>,
    pub translation_m: [f64; 3],
    pub rotation_xyzw: [f64; 4],
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
/// Every physical position remains present, including blocked positions.
pub struct RackEquipmentSlot {
    pub coordinate: String,
    pub row: usize,
    pub column: usize,
    pub blocked: bool,
    /// Centre at the carrier top, not an insertion pose or liquid surface.
    pub frame: RackEquipmentFrame,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
/// A descriptive box envelope only, not collision geometry or a fitted lid.
pub struct RackEquipmentEnvelope {
    pub frame_id: String,
    pub center_m: [f64; 3],
    pub size_m: [f64; 3],
    pub basis: String,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
/// Portable equipment foundation. Hashes bind content, not authenticity or safety.
pub struct RackEquipmentDescription {
    pub schema: String,
    pub producer_version: String,
    pub rack_id: String,
    pub asset_id: String,
    /// Template (including grid size) + drawings + projection; no instance/block/fill policy.
    pub asset_revision_sha256: String,
    pub source_sha256: String,
    /// Canonical complete record with this field empty during hashing.
    pub content_sha256: String,
    pub source: RackEquipmentSource,
    pub length_unit: String,
    pub coordinate_convention: String,
    pub base_frame: RackEquipmentFrame,
    pub envelope: RackEquipmentEnvelope,
    pub slots: Vec<RackEquipmentSlot>,
    /// Explicit unknowns, not default physical parameters.
    pub mass_kg: Option<f64>,
    pub inertia_kg_m2: Option<[f64; 6]>,
    pub permitted_actions: Vec<String>,
    pub limitations: Vec<String>,
}
