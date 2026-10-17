//! Source-backed visualization profiles and bounded rack-grid admission.
//!
//! Profiles describe drawing dimensions, not physical inventory, fit or material
//! qualification. Admission is explicit at processing boundaries, not deserialization.

use serde::{Deserialize, Serialize};

/// Resource budget for processing a rectangular grid, not a carrier capacity claim.
pub const MAX_RACK_PROCESSING_POSITIONS: usize = 1_536;

/// Admit a grid before multiplying dimensions, allocating coordinates or painting.
pub fn admit_rack_grid(rows: usize, columns: usize) -> Result<usize, String> {
    if rows == 0 || columns == 0 {
        return Err("rack_grid_empty: rows and columns must both be positive".into());
    }
    let count = rows
        .checked_mul(columns)
        .ok_or_else(|| "rack_grid_overflow: rows times columns exceeds usize".to_string())?;
    if count > MAX_RACK_PROCESSING_POSITIONS {
        return Err(format!(
            "rack_grid_limit: {rows}x{columns} exceeds the {MAX_RACK_PROCESSING_POSITIONS}-position processing budget; saved data is retained"
        ));
    }
    Ok(count)
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
/// Manufacturer-specific identity, not a nominal-volume substitute.
pub enum TubeProfileId {
    #[serde(rename = "eppendorf_safe_lock_0_5ml")]
    EppendorfSafeLock05ml,
    #[serde(rename = "eppendorf_safe_lock_1_5ml")]
    EppendorfSafeLock15ml,
    #[serde(rename = "eppendorf_safe_lock_2_0ml")]
    EppendorfSafeLock20ml,
}

impl TubeProfileId {
    /// Stable inspection order: 0.5, 1.5, then 2.0 ml.
    pub const ALL: [Self; 3] = [
        Self::EppendorfSafeLock05ml,
        Self::EppendorfSafeLock15ml,
        Self::EppendorfSafeLock20ml,
    ];

    /// The same stable identity used by the JSON contract.
    pub fn as_str(self) -> &'static str {
        match self {
            Self::EppendorfSafeLock05ml => "eppendorf_safe_lock_0_5ml",
            Self::EppendorfSafeLock15ml => "eppendorf_safe_lock_1_5ml",
            Self::EppendorfSafeLock20ml => "eppendorf_safe_lock_2_0ml",
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
/// Source status; approximate and absent dimensions must not become exact.
pub enum DimensionBasis {
    Drawing,
    ApproximateDrawing,
    NotDocumented,
}

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
/// Millimetre dimension with its source status preserved.
pub struct DrawingDimension {
    pub value_mm: Option<f64>,
    pub basis: DimensionBasis,
    /// Callout on page 1 of this profile's source drawing; None means unknown.
    pub source_callout: Option<String>,
}

impl DrawingDimension {
    fn measured(value: f64, callout: &str) -> Self {
        Self {
            value_mm: Some(value),
            basis: DimensionBasis::Drawing,
            source_callout: Some(callout.into()),
        }
    }

    fn approximate(value: f64, callout: &str) -> Self {
        Self {
            value_mm: Some(value),
            basis: DimensionBasis::ApproximateDrawing,
            source_callout: Some(callout.into()),
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
/// Qualitative silhouette only; exact bottom radii remain model assumptions.
pub enum TubeBottomShape {
    ConicalWithRoundedTip,
    Rounded,
}

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
/// Portable drawing reference, never a saved inventory or compatibility record.
pub struct TubeProfile {
    pub schema: String,
    pub id: TubeProfileId,
    pub name: String,
    pub nominal_volume_ml: f64,
    pub source_url: String,
    pub source_sha256: String,
    pub source_document: String,
    pub source_page_1based: usize,
    pub body_outer_diameter: DrawingDimension,
    pub body_inner_diameter: DrawingDimension,
    pub collar_outer_diameter: DrawingDimension,
    pub body_length: DrawingDimension,
    pub cylindrical_section_length: DrawingDimension,
    pub closed_height: DrawingDimension,
    /// Side-view width, not a cap plan outline or a verified spacing envelope.
    pub closed_side_view_width: DrawingDimension,
    pub open_90_height: DrawingDimension,
    pub open_90_side_view_width: DrawingDimension,
    pub cap_plan_outline_width: DrawingDimension,
    pub bottom_shape: TubeBottomShape,
    pub nonclaims: Vec<String>,
}

/// Static, source-cited data; never inferred from a project's container contents.
pub fn tube_profile(id: TubeProfileId) -> TubeProfile {
    let (
        volume,
        file,
        digest,
        document,
        outer,
        inner,
        collar,
        length,
        cylindrical,
        closed,
        width,
        opened,
        open_width,
        shape,
    ) = match id {
        TubeProfileId::EppendorfSafeLock05ml => (
            0.5,
            "140020_Technical-Data/Eppendorf_Consumables_Technical-data_Safe-Lock-Tube-05-mL_Safe-Lock-05-mL-technical-drawing.pdf",
            "256d42cddac90d2276f0552cffa7723ef607288461e05469914c270759cb928f",
            "0017010409-04/042016",
            6.7,
            6.1,
            7.8,
            30.0,
            16.5,
            31.0,
            18.0,
            48.0,
            16.0,
            TubeBottomShape::ConicalWithRoundedTip,
        ),
        TubeProfileId::EppendorfSafeLock15ml => (
            1.5,
            "140027_Technical-Data/Eppendorf_Consumables_Technical-data_Safe-Lock-Tube-15-mL_Safe-Lock-15-mL-technical-drawing.pdf",
            "f62ad456464047ea8108fb5117638e03036be7c08880748bf9060ebdecb7fa58",
            "0017010417-05/032024",
            9.85,
            8.7,
            10.8,
            38.9,
            20.0,
            41.0,
            20.0,
            59.0,
            19.0,
            TubeBottomShape::ConicalWithRoundedTip,
        ),
        TubeProfileId::EppendorfSafeLock20ml => (
            2.0,
            "140034_Technical-Data/Eppendorf_Consumables_Technical-data_Safe-Lock-Tube-20-mL_Safe-Lock-20-mL-technical-drawing.pdf",
            "03a3adcfe777f85662550b0d85c1a64b715ba812afc087848c33588dc09ff5c3",
            "0017010425-05/032024",
            9.85,
            8.8,
            10.8,
            40.0,
            32.3,
            42.0,
            20.0,
            60.0,
            19.0,
            TubeBottomShape::Rounded,
        ),
    };
    TubeProfile {
        schema: "gentle.tube_profile.v1".into(),
        id,
        name: format!("Eppendorf Safe-Lock {volume:.1} ml"),
        nominal_volume_ml: volume,
        source_url: format!("https://www.eppendorf.com/product-media/doc/en/{file}"),
        source_sha256: format!("sha256:{digest}"),
        source_document: document.into(),
        source_page_1based: 1,
        body_outer_diameter: DrawingDimension::measured(outer, "body outer diameter"),
        body_inner_diameter: DrawingDimension::measured(inner, "body inner diameter"),
        collar_outer_diameter: DrawingDimension::measured(collar, "collar diameter"),
        body_length: DrawingDimension::measured(length, "overall body length, cap open"),
        cylindrical_section_length: DrawingDimension::measured(
            cylindrical,
            "cylindrical section length",
        ),
        closed_height: DrawingDimension::approximate(closed, "ca. closed height"),
        closed_side_view_width: DrawingDimension::approximate(width, "ca. closed side-view width"),
        open_90_height: DrawingDimension::approximate(opened, "ca. height with cap at 90 degrees"),
        open_90_side_view_width: DrawingDimension::approximate(
            open_width,
            "ca. side-view width with cap at 90 degrees",
        ),
        cap_plan_outline_width: DrawingDimension {
            value_mm: None,
            basis: DimensionBasis::NotDocumented,
            source_callout: None,
        },
        bottom_shape: shape,
        nonclaims: vec![
            "Visualization reference, not a manufacturing drawing or a printable vessel".into(),
            "No verified fit, cap clearance, material, sterility or temperature qualification"
                .into(),
            "Not a PCR tube profile or a dimensional claim for another manufacturer's tube".into(),
        ],
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn rack_grid_admission_checks_zeros_overflow_and_resource_budget() {
        for (rows, columns) in [(0, 1), (1, 0), (0, 0)] {
            assert!(
                admit_rack_grid(rows, columns)
                    .unwrap_err()
                    .starts_with("rack_grid_empty:")
            );
        }
        assert!(
            admit_rack_grid(usize::MAX, 2)
                .unwrap_err()
                .starts_with("rack_grid_overflow:")
        );
        assert!(admit_rack_grid(1_000_000, 1_000_000).is_err());
        assert!(
            admit_rack_grid(1, 1_537)
                .unwrap_err()
                .starts_with("rack_grid_limit:")
        );
        for (rows, columns) in [(1, 1), (8, 12), (16, 24), (32, 48), (1, 1_536)] {
            assert_eq!(admit_rack_grid(rows, columns).unwrap(), rows * columns);
        }
    }

    #[test]
    fn tube_profiles_preserve_drawing_values_unknowns_and_json_identity() {
        for (id, values, shape) in [
            (
                TubeProfileId::EppendorfSafeLock05ml,
                [6.7, 6.1, 7.8, 30.0, 16.5, 31.0, 18.0, 48.0, 16.0],
                TubeBottomShape::ConicalWithRoundedTip,
            ),
            (
                TubeProfileId::EppendorfSafeLock15ml,
                [9.85, 8.7, 10.8, 38.9, 20.0, 41.0, 20.0, 59.0, 19.0],
                TubeBottomShape::ConicalWithRoundedTip,
            ),
            (
                TubeProfileId::EppendorfSafeLock20ml,
                [9.85, 8.8, 10.8, 40.0, 32.3, 42.0, 20.0, 60.0, 19.0],
                TubeBottomShape::Rounded,
            ),
        ] {
            let profile = tube_profile(id);
            let fields = [
                &profile.body_outer_diameter,
                &profile.body_inner_diameter,
                &profile.collar_outer_diameter,
                &profile.body_length,
                &profile.cylindrical_section_length,
                &profile.closed_height,
                &profile.closed_side_view_width,
                &profile.open_90_height,
                &profile.open_90_side_view_width,
            ];
            for (index, (field, value)) in fields.into_iter().zip(values).enumerate() {
                assert_eq!(
                    field.value_mm,
                    Some(value),
                    "{} dimension {index}",
                    id.as_str()
                );
                assert_eq!(
                    field.basis,
                    if index < 5 {
                        DimensionBasis::Drawing
                    } else {
                        DimensionBasis::ApproximateDrawing
                    }
                );
                assert!(
                    field
                        .source_callout
                        .as_ref()
                        .is_some_and(|callout| !callout.is_empty())
                );
            }
            assert_eq!(profile.bottom_shape, shape);
            assert_eq!(profile.cap_plan_outline_width.value_mm, None);
            assert_eq!(
                profile.cap_plan_outline_width.basis,
                DimensionBasis::NotDocumented
            );
            assert_eq!(profile.source_sha256.len(), 71);
            assert_eq!(serde_json::to_value(id).unwrap(), id.as_str());
            let encoded = serde_json::to_string(&profile).unwrap();
            assert_eq!(
                serde_json::from_str::<TubeProfile>(&encoded).unwrap(),
                profile
            );
        }
    }
}
