//! Read-only projection of existing rack geometry; no world or robot execution.

use super::*;
use crate::digest_utils::{canonical_json_string, sha256_prefixed_bytes};
use gentle_protocol::rack_equipment::*;
use std::io::Write;

#[cfg(test)]
#[path = "equipment/tests.rs"]
mod tests;

fn invalid(message: impl Into<String>) -> EngineError {
    EngineError {
        code: ErrorCode::InvalidInput,
        message: message.into(),
        cause_chain: vec![],
    }
}

fn digest(value: &impl serde::Serialize) -> Result<String, EngineError> {
    canonical_json_string(value)
        .map(|json| sha256_prefixed_bytes(json.as_bytes()))
        .map_err(|err| invalid(format!("rack_equipment_serialization: {err}")))
}

impl GentleEngine {
    /// Project one saved rack without changing it or inspecting any external files.
    pub fn rack_equipment_description(
        &self,
        rack_id: &str,
        template: RackPhysicalTemplateKind,
    ) -> Result<RackEquipmentDescription, EngineError> {
        let rack = self
            .state
            .container_state
            .racks
            .get(rack_id)
            .ok_or_else(|| EngineError {
                code: ErrorCode::NotFound,
                message: format!("Rack '{rack_id}' not found"),
                cause_chain: vec![],
            })?;
        let spec = Self::rack_physical_template_spec(
            template,
            &rack.profile,
            rack.carrier_pitch_override_um,
        )?;
        let source = RackEquipmentSource {
            projection: RACK_EQUIPMENT_PROJECTION.into(),
            tube_profiles: spec
                .tube_profiles
                .iter()
                .copied()
                .map(gentle_protocol::labware::tube_profile)
                .collect(),
            physical_template: spec,
            rack_profile: rack.profile.clone(),
        };
        Self::build_rack_equipment_description(rack_id, source, env!("CARGO_PKG_VERSION"))
    }

    fn build_rack_equipment_description(
        rack_id: &str,
        mut source: RackEquipmentSource,
        producer_version: &str,
    ) -> Result<RackEquipmentDescription, EngineError> {
        let count = source
            .rack_profile
            .processing_slot_count()
            .map_err(invalid)?;
        let spec = &source.physical_template;
        if source.projection != RACK_EQUIPMENT_PROJECTION
            || rack_id.trim().is_empty()
            || producer_version.trim().is_empty()
            || spec.rows != source.rack_profile.rows
            || spec.columns != source.rack_profile.columns
            || source.rack_profile.coordinate_scheme != "a1"
        {
            return Err(invalid(
                "rack_equipment_source_invalid: projection, grid or identity mismatch",
            ));
        }
        if source
            .tube_profiles
            .iter()
            .map(|p| p.id)
            .collect::<Vec<_>>()
            != spec.tube_profiles
        {
            return Err(invalid("rack_equipment_profile_binding_mismatch"));
        }
        let positive = [
            spec.overall_width_mm,
            spec.overall_depth_mm,
            spec.rack_height_mm,
            spec.pitch_x_mm,
            spec.pitch_y_mm,
            spec.opening_diameter_mm,
        ];
        let nonnegative = [spec.edge_margin_mm, spec.front_top_clearance_mm];
        if positive.iter().any(|v| !v.is_finite() || *v <= 0.0)
            || nonnegative.iter().any(|v| !v.is_finite() || *v < 0.0)
        {
            return Err(invalid(
                "rack_equipment_dimensions_invalid: finite positive geometry required",
            ));
        }
        // A policy set has one canonical order; physical row order never changes.
        source.rack_profile.blocked_coordinates.sort();
        source.rack_profile.blocked_coordinates.dedup();
        let spec = &source.physical_template;
        let blocked = &source.rack_profile.blocked_coordinates;
        let mut observed_blocked = 0;
        let mut slots = Vec::with_capacity(count);
        for row in 0..spec.rows {
            for column in 0..spec.columns {
                let coordinate =
                    Self::rack_coordinate_from_row_column(&source.rack_profile, row, column)?;
                let (x, y) = Self::rack_physical_hole_center_mm(spec, row, column);
                if !x.is_finite()
                    || !y.is_finite()
                    || x < 0.0
                    || x > spec.overall_width_mm
                    || y < 0.0
                    || y > spec.overall_depth_mm
                {
                    return Err(invalid("rack_equipment_slot_outside_envelope"));
                }
                let is_blocked = blocked.binary_search(&coordinate).is_ok();
                observed_blocked += usize::from(is_blocked);
                slots.push(RackEquipmentSlot {
                    row,
                    column,
                    blocked: is_blocked,
                    frame: RackEquipmentFrame {
                        id: format!("slot_{coordinate}"),
                        parent: Some("rack_base".into()),
                        translation_m: [
                            f64::from(x) / 1000.0,
                            f64::from(y) / 1000.0,
                            f64::from(spec.rack_height_mm) / 1000.0,
                        ],
                        rotation_xyzw: [0.0, 0.0, 0.0, 1.0],
                    },
                    coordinate,
                });
            }
        }
        if observed_blocked != blocked.len() {
            return Err(invalid("rack_equipment_blocked_coordinate_invalid"));
        }
        let size_m = [
            spec.overall_width_mm,
            spec.overall_depth_mm,
            spec.rack_height_mm,
        ]
        .map(|v| f64::from(v) / 1000.0);
        let mut description = RackEquipmentDescription {
            schema: RACK_EQUIPMENT_SCHEMA.into(), producer_version: producer_version.into(),
            rack_id: rack_id.into(), asset_id: format!("gentle_builtin_rack:{}", spec.kind.as_str()),
            asset_revision_sha256: digest(&(&source.projection, spec, &source.tube_profiles))?,
            source_sha256: digest(&source)?, content_sha256: String::new(),
            length_unit: "metre".into(),
            coordinate_convention: "right_handed; origin=front_left_base; +x=increasing_columns; +y=increasing_rows; +z=up; slot_origin=carrier_top; quaternion=xyzw".into(),
            base_frame: RackEquipmentFrame {
                id: "rack_base".into(), parent: None, translation_m: [0.0; 3],
                rotation_xyzw: [0.0, 0.0, 0.0, 1.0],
            },
            envelope: RackEquipmentEnvelope {
                frame_id: "rack_base".into(), center_m: size_m.map(|v| v / 2.0), size_m,
                basis: "unverified_existing_template_body_envelope".into(),
            },
            slots, mass_kg: None, inertia_kg_m2: None, permitted_actions: vec![],
            limitations: vec![
                "Geometry description only; not a simulation run, execution plan or robot-motion approval.".into(),
                "No collision/inertial model, grasp pose, liquid frame, tube insertion depth or world placement is inferred.".into(),
                "Carrier-top slot frames are not validated interaction points. Blocked positions are still physical slots.".into(),
                "Envelope describes the rack body, not detached lids, caps, hinges or contained tubes; no fit or collision-clearance claim.".into(),
                "No calibration, physical safety, manufacturing, biological performance or equipment capability is established.".into(),
                "No sample names, sequences, occupancy or private room geometry are exported. Rack ID and blocked-position policy are explicit instance metadata.".into(),
                "Source dimensions remain millimetres; all derived frames and envelope dimensions are metres. Hashes bind canonical content, not authenticity.".into(),
            ],
            source,
        };
        description.content_sha256 = digest(&description)?;
        Ok(description)
    }

    /// Rebuild geometry and all bindings from the embedded source, without I/O.
    /// This establishes internal consistency, never authenticity of supplied drawings.
    pub fn verify_rack_equipment_description(
        description: &RackEquipmentDescription,
    ) -> Result<(), EngineError> {
        let rebuilt = Self::build_rack_equipment_description(
            &description.rack_id,
            description.source.clone(),
            &description.producer_version,
        )?;
        if canonical_json_string(&rebuilt).map_err(|e| invalid(e.to_string()))?
            != canonical_json_string(description).map_err(|e| invalid(e.to_string()))?
        {
            return Err(invalid(
                "rack_equipment_binding_mismatch: description differs from source projection",
            ));
        }
        Ok(())
    }

    pub(in crate::engine) fn export_rack_equipment_json(
        &self,
        rack_id: &str,
        template: RackPhysicalTemplateKind,
        path: &str,
    ) -> Result<RackEquipmentDescription, EngineError> {
        let description = self.rack_equipment_description(rack_id, template)?;
        // Match OpResult's JSON-value projection, including embedded f32 values.
        let value = serde_json::to_value(&description).map_err(|e| invalid(e.to_string()))?;
        let bytes = serde_json::to_vec_pretty(&value).map_err(|e| invalid(e.to_string()))?;
        let destination = std::path::Path::new(path);
        let parent = destination
            .parent()
            .filter(|p| !p.as_os_str().is_empty())
            .unwrap_or_else(|| std::path::Path::new("."));
        let io_error = |err: std::io::Error| EngineError {
            code: ErrorCode::Io,
            message: format!("Could not publish rack equipment JSON '{path}': {err}"),
            cause_chain: vec![],
        };
        let mut file = tempfile::NamedTempFile::new_in(parent).map_err(io_error)?;
        file.write_all(&bytes).map_err(io_error)?;
        file.as_file().sync_all().map_err(io_error)?;
        // Publish without replacing another file; tempfile owns failed-write cleanup.
        let staged = file.into_temp_path();
        staged
            .persist_noclobber(destination)
            .map_err(|e| io_error(e.error))?;
        Ok(description)
    }
}
