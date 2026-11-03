//! Read-only projection of existing rack geometry; no world or robot execution.

use super::*;
use crate::digest_utils::{canonical_json_string, sha256_prefixed_bytes};
use gentle_protocol::rack_equipment::*;
use std::io::Read;

const VALIDATION_INPUT_LIMIT: u64 = 16 * 1024 * 1024;

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

fn same_json(left: &impl serde::Serialize, right: &impl serde::Serialize) -> bool {
    match (canonical_json_string(left), canonical_json_string(right)) {
        (Ok(left), Ok(right)) => left == right,
        _ => false,
    }
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
            carrier_pitch_override_um: rack.carrier_pitch_override_um,
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
        // Typed callers must not turn forbidden Some(NaN) mass into JSON null.
        if description.mass_kg.is_some() || description.inertia_kg_m2.is_some() {
            return Err(invalid(
                "rack_equipment_binding_mismatch: v1 mass and inertia must be null",
            ));
        }
        let bytes = serde_json::to_vec(description).map_err(|e| invalid(e.to_string()))?;
        let validation = Self::validate_rack_equipment_bytes(&bytes);
        if !validation.internally_consistent {
            return Err(invalid(format!(
                "rack_equipment_binding_mismatch: {}",
                validation
                    .diagnostics
                    .first()
                    .map(String::as_str)
                    .unwrap_or("rack_equipment_source_invalid")
            )));
        }
        Ok(())
    }

    /// Validate the complete JSON value and its embedded-source projection.
    /// No parser diagnostics or source values are copied into finding codes.
    pub fn validate_rack_equipment_bytes(bytes: &[u8]) -> RackEquipmentValidation {
        let mut result = RackEquipmentValidation {
            schema: RACK_EQUIPMENT_VALIDATION_SCHEMA.into(),
            input_sha256: sha256_prefixed_bytes(bytes),
            description_schema: None, claimed_content_sha256: None,
            internally_consistent: false, matches_builtin_catalog: false,
            diagnostics: vec![],
            nonclaims: vec![
                "Consistency and built-in catalog match only; not drawing authenticity, fit, safety or robot authority.".into(),
                "Catalog matching refers to this producer's built-in templates and drawings, not installed equipment or calibration.".into(),
            ],
        };
        let raw: serde_json::Value = match serde_json::from_slice(bytes) {
            Ok(value) => value,
            Err(_) => {
                result
                    .diagnostics
                    .push("rack_equipment_json_invalid".into());
                return result;
            }
        };
        result.description_schema = raw
            .get("schema")
            .and_then(|v| v.as_str())
            .map(str::to_owned);
        result.claimed_content_sha256 = raw
            .get("content_sha256")
            .and_then(|v| v.as_str())
            .map(str::to_owned);
        if result.description_schema.as_deref() != Some(RACK_EQUIPMENT_SCHEMA) {
            result
                .diagnostics
                .push("rack_equipment_schema_unsupported".into());
            return result;
        }
        let typed: RackEquipmentDescription = match serde_json::from_value(raw.clone()) {
            Ok(value) => value,
            Err(_) => {
                result
                    .diagnostics
                    .push("rack_equipment_record_invalid".into());
                return result;
            }
        };
        if !same_json(&typed, &raw) {
            result
                .diagnostics
                .push("rack_equipment_unknown_or_noncanonical_fields".into());
            return result;
        }
        let rebuilt = match Self::build_rack_equipment_description(
            &typed.rack_id,
            typed.source.clone(),
            &typed.producer_version,
        ) {
            Ok(value) => value,
            Err(_) => {
                result
                    .diagnostics
                    .push("rack_equipment_source_invalid".into());
                return result;
            }
        };
        if !same_json(&rebuilt, &raw) {
            result
                .diagnostics
                .push("rack_equipment_binding_mismatch".into());
            return result;
        }
        result.internally_consistent = true;
        let expected = Self::rack_physical_template_spec(
            typed.source.physical_template.kind,
            &typed.source.rack_profile,
            typed.source.carrier_pitch_override_um,
        );
        if let Ok(spec) = expected {
            let profiles = spec
                .tube_profiles
                .iter()
                .copied()
                .map(gentle_protocol::labware::tube_profile)
                .collect::<Vec<_>>();
            result.matches_builtin_catalog = same_json(&spec, &typed.source.physical_template)
                && same_json(&profiles, &typed.source.tube_profiles);
        }
        if !result.matches_builtin_catalog {
            result
                .diagnostics
                .push("rack_equipment_catalog_mismatch".into());
        }
        result
    }

    pub(in crate::engine) fn validate_rack_equipment_json(
        &self,
        path: &str,
    ) -> Result<RackEquipmentValidation, EngineError> {
        let io_error = |err: std::io::Error| EngineError {
            code: ErrorCode::Io,
            message: format!("rack_equipment_validate_input: {err}"),
            cause_chain: vec![],
        };
        if !std::fs::metadata(path).map_err(io_error)?.is_file() {
            return Err(invalid("rack_equipment_validate_input_not_regular_file"));
        }
        let file = std::fs::File::open(path).map_err(io_error)?;
        let mut bytes = Vec::new();
        file.take(VALIDATION_INPUT_LIMIT + 1)
            .read_to_end(&mut bytes)
            .map_err(io_error)?;
        if bytes.len() as u64 > VALIDATION_INPUT_LIMIT {
            return Err(invalid("rack_equipment_validate_input_too_large"));
        }
        Ok(Self::validate_rack_equipment_bytes(&bytes))
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
        self.json_report_destination(path, "rack_equipment_export")?
            .write(&value)?;
        Ok(description)
    }
}
