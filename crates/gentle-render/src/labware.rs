//! Simplified drawing-backed outer-body visualization, never printable vessels.

use gentle_protocol::labware::{
    DimensionBasis, TubeBottomShape, TubeProfile, TubeProfileId, tube_profile,
};

/// Deterministic OpenSCAD for one reference tube or a full PCR strip.
///
/// Cap/hinge outline and wall/bottom radii are insufficiently specified for a
/// manufacturing model. The solid outer silhouette omits the cap, hinge and
/// cavity and explicitly discloses its illustrative tip/collar assumptions.
pub fn render_tube_openscad(id: TubeProfileId) -> Result<String, String> {
    render_tube_openscad_with_count(id, None)
}

/// A strip defaults to all eight tubes; individual-tube exports retain their bytes.
pub fn render_tube_openscad_with_count(
    id: TubeProfileId,
    tube_count: Option<usize>,
) -> Result<String, String> {
    let profile = tube_profile(id);
    let count = profile.model_tube_count(tube_count)?;
    if profile.strip.is_some() {
        return render_strip_profile_openscad(&profile, count);
    }
    render_tube_profile_openscad(&profile)
}

/// Validate required drawing dimensions before producing any model bytes.
pub fn render_tube_profile_openscad(profile: &TubeProfile) -> Result<String, String> {
    if profile.strip.is_some() {
        return render_strip_profile_openscad(profile, profile.model_tube_count(None)?);
    }
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

fn dimension(
    value: &gentle_protocol::labware::DrawingDimension,
    name: &str,
) -> Result<f64, String> {
    value
        .value_mm
        .filter(|value| value.is_finite() && *value > 0.0)
        .filter(|_| value.basis != DimensionBasis::NotDocumented)
        .ok_or_else(|| format!("tube_model_dimension_unavailable: {name}"))
}

fn render_strip_profile_openscad(
    profile: &TubeProfile,
    tube_count: usize,
) -> Result<String, String> {
    profile.model_tube_count(Some(tube_count))?;
    let strip = profile
        .strip
        .as_ref()
        .ok_or_else(|| "tube_model_strip_invalid: strip layout is required".to_string())?;
    let height = dimension(&profile.body_length, "body length")?;
    let body = dimension(&profile.body_outer_diameter, "body outer diameter")?;
    let collar = dimension(&profile.collar_outer_diameter, "collar diameter")?;
    let pitch = dimension(&strip.center_pitch, "strip center pitch")?;
    let rim = dimension(&strip.rim_outer_diameter, "strip rim diameter")?;
    let tip = dimension(&strip.tip_outer_diameter, "strip tip diameter")?;
    if height <= 6.0
        || tip >= body
        || body > rim
        || rim > collar
        || collar >= pitch
        || tip / 2.0 >= height * 0.55
        || height * 0.55 >= height - 4.5
    {
        return Err(
            "tube_model_geometry_invalid: strip silhouette invariants are not satisfied".into(),
        );
    }
    let mut source = format!(
        "// GENtle PCR strip visualization only; NOT a printable vessel or fit proof.\n\
         // Nominal volume 0.2 ml PER TUBE. Separate caps, cavity, numbering and end tab/cut stubs omitted.\n\
         // Contiguous retained tubes; only internal bridges remain. No physical cutting tolerance is modeled.\n\
         // Illustrative assumptions, not drawing dimensions: taper transition at 55% of height; spherical tip;\n\
         // collar 4.5 mm below rim, thickness 0.9 mm; rim thickness 0.6 mm; bridges 1.0 x 0.7 mm.\n\
         // Unknown cylindrical length and cap geometry remain unknown in the profile.\n\
         // profile_json={}\n\
         $fn = 64;\n\
         tube_count = {tube_count};\n\
         max_tubes = {};\n\
         center_pitch_mm = {pitch:.3};\n\
         body_length_mm = {height:.3};\n\
         body_diameter_mm = {body:.3};\n\
         collar_diameter_mm = {collar:.3};\n\
         rim_diameter_mm = {rim:.3};\n\
         tip_diameter_mm = {tip:.3};\n\
         visual_shoulder_fraction = 0.55;\n\
         visual_collar_drop_mm = 4.5;\n\
         visual_collar_thickness_mm = 0.9;\n\
         visual_rim_thickness_mm = 0.6;\n\
         visual_bridge_width_mm = 1.0;\n\
         visual_bridge_height_mm = 0.7;\n\n",
        serde_json::to_string(profile)
            .map_err(|error| format!("tube_model_profile_invalid: {error}"))?,
        strip.max_tubes,
    );
    source.push_str(PCR_STRIP_SCAD);
    Ok(source)
}

const PCR_STRIP_SCAD: &str = r#"assert(tube_count >= 1 && tube_count <= max_tubes && tube_count == floor(tube_count));
assert(max_tubes <= 8);
module gentle_pcr_strip_tube() {
    shoulder = body_length_mm * visual_shoulder_fraction;
    tip_radius = tip_diameter_mm / 2;
    union() {
        translate([0, 0, tip_radius]) sphere(r=tip_radius);
        translate([0, 0, tip_radius])
            cylinder(h=shoulder - tip_radius, r1=tip_radius, r2=body_diameter_mm / 2);
        translate([0, 0, shoulder])
            cylinder(h=body_length_mm - shoulder, d=body_diameter_mm);
        translate([0, 0, body_length_mm - visual_rim_thickness_mm])
            cylinder(h=visual_rim_thickness_mm, d=rim_diameter_mm);
        translate([0, 0, body_length_mm - visual_collar_drop_mm])
            cylinder(h=visual_collar_thickness_mm, d=collar_diameter_mm);
    }
}

union() {
    for (i = [0 : tube_count - 1])
        translate([i * center_pitch_mm, 0, 0]) gentle_pcr_strip_tube();
    if (tube_count > 1)
        for (i = [0 : tube_count - 2])
            translate([i * center_pitch_mm, -visual_bridge_width_mm / 2,
                       body_length_mm - visual_collar_drop_mm])
                cube([center_pitch_mm, visual_bridge_width_mm, visual_bridge_height_mm]);
}
"#;

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
    fn labware_pcr_strip_models_cover_all_cut_lengths_deterministically() {
        let id = TubeProfileId::BrandPcrStrip02ml781320;
        for count in 1..=8 {
            let source = render_tube_openscad_with_count(id, Some(count)).unwrap();
            assert_eq!(
                source,
                render_tube_openscad_with_count(id, Some(count)).unwrap()
            );
            assert!(source.contains(&format!("tube_count = {count};")));
            assert!(source.contains("center_pitch_mm = 9.000;"));
            assert!(source.contains("body_length_mm = 20.750;"));
            assert!(source.contains("if (tube_count > 1)"));
            assert!(source.contains("[0 : tube_count - 2]"));
            assert!(source.contains("NOT a printable vessel or fit proof"));
            assert!(source.contains("Illustrative assumptions, not drawing dimensions"));
            let encoded = source
                .lines()
                .find_map(|line| line.strip_prefix("// profile_json="))
                .unwrap();
            assert_eq!(
                serde_json::from_str::<TubeProfile>(encoded).unwrap(),
                tube_profile(id)
            );
        }
        assert_eq!(
            render_tube_openscad(id).unwrap(),
            render_tube_openscad_with_count(id, Some(8)).unwrap()
        );
        for count in [0, 9, usize::MAX] {
            assert!(render_tube_openscad_with_count(id, Some(count)).is_err());
        }
        for id in [
            TubeProfileId::EppendorfSafeLock05ml,
            TubeProfileId::EppendorfSafeLock15ml,
            TubeProfileId::EppendorfSafeLock20ml,
        ] {
            assert_eq!(
                render_tube_openscad(id).unwrap(),
                render_tube_openscad_with_count(id, Some(1)).unwrap()
            );
            assert!(render_tube_openscad_with_count(id, Some(2)).is_err());
        }
    }

    #[test]
    fn labware_pcr_strip_refuses_missing_geometry_and_keeps_unknowns_unknown() {
        for field in 0..6 {
            for value in [
                None,
                Some(f64::NAN),
                Some(f64::INFINITY),
                Some(0.0),
                Some(-1.0),
            ] {
                let mut profile = tube_profile(TubeProfileId::BrandPcrStrip02ml781320);
                match field {
                    0 => profile.body_length.value_mm = value,
                    1 => profile.body_outer_diameter.value_mm = value,
                    2 => profile.collar_outer_diameter.value_mm = value,
                    3 => profile.strip.as_mut().unwrap().center_pitch.value_mm = value,
                    4 => profile.strip.as_mut().unwrap().rim_outer_diameter.value_mm = value,
                    _ => profile.strip.as_mut().unwrap().tip_outer_diameter.value_mm = value,
                }
                assert!(render_tube_profile_openscad(&profile).is_err());
            }
        }
        let mut profile = tube_profile(TubeProfileId::BrandPcrStrip02ml781320);
        profile.strip.as_mut().unwrap().center_pitch.basis = DimensionBasis::NotDocumented;
        assert!(render_tube_profile_openscad(&profile).is_err());
        profile.strip.as_mut().unwrap().center_pitch.basis = DimensionBasis::Drawing;
        profile.strip.as_mut().unwrap().center_pitch.value_mm = Some(7.0);
        assert!(render_tube_profile_openscad(&profile).is_err());
    }

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
