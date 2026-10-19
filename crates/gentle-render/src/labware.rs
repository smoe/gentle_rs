//! Simplified drawing-backed outer-body visualization, never printable vessels.

use gentle_protocol::labware::{
    DimensionBasis, TubeBottomShape, TubeProfile, TubeProfileId, tube_profile,
};

/// Deterministic OpenSCAD for one reference tube body and collar.
///
/// Cap/hinge outline and wall/bottom radii are insufficiently specified for a
/// manufacturing model. The solid outer silhouette omits the cap, hinge and
/// cavity and explicitly discloses its illustrative tip/collar assumptions.
pub fn render_tube_openscad(id: TubeProfileId) -> Result<String, String> {
    let profile = tube_profile(id);
    render_tube_profile_openscad(&profile)
}

/// Validate required drawing dimensions before producing any model bytes.
pub fn render_tube_profile_openscad(profile: &TubeProfile) -> Result<String, String> {
    let dimension = |value: &gentle_protocol::labware::DrawingDimension, name: &str| {
        value
            .value_mm
            .filter(|value| value.is_finite() && *value > 0.0)
            .filter(|_| value.basis != DimensionBasis::NotDocumented)
            .ok_or_else(|| format!("tube_model_dimension_unavailable: {name}"))
    };
    let length = dimension(&profile.body_length, "body length")?;
    let outer = dimension(&profile.body_outer_diameter, "body outer diameter")?;
    let collar = dimension(&profile.collar_outer_diameter, "collar diameter")?;
    let cylinder = dimension(&profile.cylindrical_section_length, "cylindrical section")?;
    if length <= 1.0 || length - cylinder <= 0.6 || collar < outer {
        return Err(
            "tube_model_geometry_invalid: required silhouette invariants are not satisfied".into(),
        );
    }
    let rounded = profile.bottom_shape == TubeBottomShape::Rounded;
    let mut source = format!(
        "// GENtle tube visualization only; NOT a printable vessel or fit proof.\n\
         // Cap, hinge and inner cavity omitted; cap plan outline is unknown.\n\
         // Illustrative assumptions: collar thickness 1 mm; conical tip radius 0.6 mm; rounded bottom is an ellipsoid.\n\
         // profile_json={}\n\
         $fn = 64;\n\
         body_length_mm = {:.3};\n\
         body_diameter_mm = {:.3};\n\
         collar_diameter_mm = {:.3};\n\
         cylindrical_length_mm = {:.3};\n\
         rounded_bottom = {};\n\
         visual_collar_thickness_mm = 1.0;\n\
         visual_tip_radius_mm = 0.6;\n\n",
        serde_json::to_string(profile)
            .map_err(|error| format!("tube_model_profile_invalid: {error}"))?,
        length,
        outer,
        collar,
        cylinder,
        rounded,
    );
    source.push_str(TUBE_BODY_SCAD);
    Ok(source)
}

const TUBE_BODY_SCAD: &str = r#"module gentle_tube_outer_body() {
    radius = body_diameter_mm / 2;
    lower_height = body_length_mm - cylindrical_length_mm;
    assert(lower_height > visual_tip_radius_mm);
    assert(collar_diameter_mm >= body_diameter_mm);
    assert(body_length_mm > visual_collar_thickness_mm);
    union() {
        if (rounded_bottom) {
            translate([0, 0, lower_height])
                scale([radius, radius, lower_height]) sphere(r=1);
        } else {
            translate([0, 0, visual_tip_radius_mm]) sphere(r=visual_tip_radius_mm);
            translate([0, 0, visual_tip_radius_mm])
                cylinder(h=lower_height - visual_tip_radius_mm,
                         r1=visual_tip_radius_mm, r2=radius);
        }
        translate([0, 0, lower_height])
            cylinder(h=cylindrical_length_mm, r=radius);
        translate([0, 0, body_length_mm - visual_collar_thickness_mm])
            cylinder(h=visual_collar_thickness_mm, d=collar_diameter_mm);
    }
}

gentle_tube_outer_body();
"#;

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn tube_models_bind_their_source_dimensions_and_disclose_omissions() {
        for (id, length, diameter, collar, cylindrical, rounded) in [
            (
                TubeProfileId::EppendorfSafeLock05ml,
                "30.000",
                "6.700",
                "7.800",
                "16.500",
                false,
            ),
            (
                TubeProfileId::EppendorfSafeLock15ml,
                "38.900",
                "9.850",
                "10.800",
                "20.000",
                false,
            ),
            (
                TubeProfileId::EppendorfSafeLock20ml,
                "40.000",
                "9.850",
                "10.800",
                "32.300",
                true,
            ),
        ] {
            let source = render_tube_openscad(id).unwrap();
            assert_eq!(source, render_tube_openscad(id).unwrap());
            for declaration in [
                format!("body_length_mm = {length};"),
                format!("body_diameter_mm = {diameter};"),
                format!("collar_diameter_mm = {collar};"),
                format!("cylindrical_length_mm = {cylindrical};"),
                format!("rounded_bottom = {rounded};"),
            ] {
                assert!(source.contains(&declaration));
            }
            let encoded = source
                .lines()
                .find_map(|line| line.strip_prefix("// profile_json="))
                .unwrap();
            let decoded: gentle_protocol::labware::TubeProfile =
                serde_json::from_str(encoded).unwrap();
            assert_eq!(decoded, tube_profile(id));
            assert!(source.contains("NOT a printable vessel or fit proof"));
            assert!(source.contains("Cap, hinge and inner cavity omitted"));
            assert!(source.contains("Illustrative assumptions:"));
            assert!(source.ends_with("gentle_tube_outer_body();\n"));
        }
    }

    #[test]
    fn labware_tube_model_refuses_missing_or_unproven_dimensions_without_panicking() {
        let original = tube_profile(TubeProfileId::EppendorfSafeLock15ml);
        for field in 0..4 {
            for value in [
                None,
                Some(f64::NAN),
                Some(f64::INFINITY),
                Some(0.0),
                Some(-1.0),
            ] {
                let mut profile = original.clone();
                match field {
                    0 => profile.body_length.value_mm = value,
                    1 => profile.body_outer_diameter.value_mm = value,
                    2 => profile.collar_outer_diameter.value_mm = value,
                    _ => profile.cylindrical_section_length.value_mm = value,
                }
                assert!(
                    render_tube_profile_openscad(&profile)
                        .unwrap_err()
                        .starts_with("tube_model_dimension_unavailable:")
                );
            }
        }
        assert!(
            render_tube_openscad(TubeProfileId::EppendorfPcrTube02ml)
                .unwrap_err()
                .contains("collar diameter")
        );
        let mut profile = original;
        profile.body_length.basis = DimensionBasis::NotDocumented;
        assert!(render_tube_profile_openscad(&profile).is_err());
    }
}
