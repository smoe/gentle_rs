//! Source-backed visualization profiles and bounded rack-grid admission.
//!
//! Profiles describe drawing dimensions, not physical inventory, fit or material
//! qualification. Admission is explicit at processing boundaries, not deserialization.

use serde::{Deserialize, Serialize};

/// Resource budget for processing a rectangular grid, not a carrier capacity claim.
pub const MAX_RACK_PROCESSING_POSITIONS: usize = 1_536;

/// Engineering assumptions in millimetres, not manufacturing tolerances.
pub const RADIAL_CLEARANCE_COLLAR: f64 = 0.15;
pub const RADIAL_CLEARANCE_FLOOR: f64 = 0.30;
pub const MIN_COLLAR_LEDGE: f64 = 0.20;
pub const TIP_CLEARANCE: f64 = 2.0;
pub const CYLINDER_ENGAGEMENT: f64 = 3.0;
pub const COLLAR_TOP_MARGIN: f64 = 5.0;

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
/// Intended support only; neither mode establishes physical fit.
pub enum RackTubeSupportMode {
    CollarSeated,
    FloorSeated,
}

impl RackTubeSupportMode {
    pub fn as_str(self) -> &'static str {
        match self {
            Self::CollarSeated => "collar_seated",
            Self::FloorSeated => "floor_seated",
        }
    }
}

#[derive(Debug, Clone, PartialEq)]
/// Pure drawing-derived straight-bore parameters; dimensions are millimetres.
pub struct TubeCarrierGeometry {
    pub opening_diameter_mm: f64,
    pub bore_depth_mm: f64,
    pub pitch_mm: f64,
    pub edge_margin_mm: f64,
}

fn required_dimension(value: &DrawingDimension, name: &str) -> Result<f64, String> {
    value
        .value_mm
        .filter(|value| value.is_finite() && *value > 0.0)
        .filter(|_| value.basis != DimensionBasis::NotDocumented)
        .ok_or_else(|| format!("tube_dimension_unavailable: {name}"))
}

impl TubeCarrierGeometry {
    /// Recheck a proposed bore against every profile, including synthetic inputs.
    pub fn validate(
        &self,
        profiles: &[TubeProfile],
        mode: RackTubeSupportMode,
    ) -> Result<(), String> {
        if profiles.is_empty() {
            return Err("tube_profiles_empty: carrier needs drawing profiles".into());
        }
        for value in [
            self.opening_diameter_mm,
            self.bore_depth_mm,
            self.pitch_mm,
            self.edge_margin_mm,
        ] {
            if !value.is_finite() || value <= 0.0 {
                return Err(
                    "carrier_geometry_invalid: dimensions must be finite and positive".into(),
                );
            }
        }
        for profile in profiles {
            let body = required_dimension(&profile.body_outer_diameter, "body outer diameter")?;
            let length = required_dimension(&profile.body_length, "body length")?;
            if self.opening_diameter_mm <= body || self.opening_diameter_mm >= self.pitch_mm {
                return Err(
                    "carrier_opening_invalid: body clearance or pitch is insufficient".into(),
                );
            }
            match mode {
                RackTubeSupportMode::CollarSeated => {
                    let collar =
                        required_dimension(&profile.collar_outer_diameter, "collar diameter")?;
                    if self.opening_diameter_mm > collar - 2.0 * MIN_COLLAR_LEDGE {
                        return Err("collar_ledge_insufficient: opening leaves less than 0.20 mm radial ledge".into());
                    }
                    if self.bore_depth_mm < length + TIP_CLEARANCE {
                        return Err(
                            "tip_clearance_insufficient: suspended tip needs 2.0 mm clearance"
                                .into(),
                        );
                    }
                }
                RackTubeSupportMode::FloorSeated => {
                    let cylinder = required_dimension(
                        &profile.cylindrical_section_length,
                        "cylindrical section",
                    )?;
                    if cylinder >= length {
                        return Err(
                            "tube_dimension_invalid: cylinder must be shorter than body".into()
                        );
                    }
                    if self.bore_depth_mm <= length - cylinder + CYLINDER_ENGAGEMENT {
                        return Err("bore_engages_cone_only: require more than 3.0 mm cylindrical engagement".into());
                    }
                    if self.bore_depth_mm >= length - COLLAR_TOP_MARGIN {
                        return Err(
                            "collar_not_above_top: require more than 5.0 mm above the bore".into(),
                        );
                    }
                }
            }
        }
        Ok(())
    }
}

/// Conservative Safe-Lock carrier; side-view width is not a measured cap outline.
pub fn derive_tube_carrier_geometry(
    profiles: &[TubeProfile],
    mode: RackTubeSupportMode,
) -> Result<TubeCarrierGeometry, String> {
    derive_carrier_geometry(profiles, mode, None)
}

/// Explicit 8x12 layout; 9-mm cap compatibility is assumed, never verified.
pub fn derive_pcr96_carrier_geometry(
    profile: &TubeProfile,
    rows: usize,
    columns: usize,
) -> Result<TubeCarrierGeometry, String> {
    if rows != 8 || columns != 12 {
        return Err(
            "carrier_requires_8x12: PCR carrier needs exactly 8 rows and 12 columns".into(),
        );
    }
    derive_carrier_geometry(
        std::slice::from_ref(profile),
        RackTubeSupportMode::FloorSeated,
        Some(9.0),
    )
}

fn derive_carrier_geometry(
    profiles: &[TubeProfile],
    mode: RackTubeSupportMode,
    fixed_pitch: Option<f64>,
) -> Result<TubeCarrierGeometry, String> {
    if profiles.is_empty() {
        return Err("tube_profiles_empty: carrier needs drawing profiles".into());
    }
    let mut body: f64 = 0.0;
    let mut length: f64 = 0.0;
    let mut lower: f64 = 0.0;
    let mut pitch: f64 = 0.0;
    for profile in profiles {
        body = body.max(required_dimension(
            &profile.body_outer_diameter,
            "body outer diameter",
        )?);
        let profile_length = required_dimension(&profile.body_length, "body length")?;
        length = length.max(profile_length);
        if mode == RackTubeSupportMode::FloorSeated {
            let cylinder =
                required_dimension(&profile.cylindrical_section_length, "cylindrical section")?;
            if cylinder >= profile_length {
                return Err("tube_dimension_invalid: cylinder must be shorter than body".into());
            }
            lower = lower.max(profile_length - cylinder);
        }
        if let Some(fixed) = fixed_pitch {
            pitch = fixed;
        } else {
            pitch = pitch.max(
                required_dimension(&profile.closed_side_view_width, "closed side-view width")
                    .map_err(|_| {
                        "cap_envelope_unknown: closed side-view width is unavailable".to_string()
                    })?,
            );
        }
    }
    let (clearance, depth) = match mode {
        RackTubeSupportMode::CollarSeated => (RADIAL_CLEARANCE_COLLAR, length + TIP_CLEARANCE),
        // The strict lower bound is rounded to the next whole millimetre.
        RackTubeSupportMode::FloorSeated => (
            RADIAL_CLEARANCE_FLOOR,
            (lower + CYLINDER_ENGAGEMENT).floor() + 1.0,
        ),
    };
    let opening = body + 2.0 * clearance;
    let geometry = TubeCarrierGeometry {
        opening_diameter_mm: opening,
        bore_depth_mm: depth,
        pitch_mm: pitch,
        edge_margin_mm: (pitch - opening) / 2.0,
    };
    geometry.validate(profiles, mode)?;
    Ok(geometry)
}

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
    #[serde(rename = "eppendorf_pcr_tube_0_2ml")]
    EppendorfPcrTube02ml,
}

impl TubeProfileId {
    /// Stable inspection order: Safe-Lock 0.5/1.5/2.0 ml, then PCR 0.2 ml.
    pub const ALL: [Self; 4] = [
        Self::EppendorfSafeLock05ml,
        Self::EppendorfSafeLock15ml,
        Self::EppendorfSafeLock20ml,
        Self::EppendorfPcrTube02ml,
    ];

    /// The same stable identity used by the JSON contract.
    pub fn as_str(self) -> &'static str {
        match self {
            Self::EppendorfSafeLock05ml => "eppendorf_safe_lock_0_5ml",
            Self::EppendorfSafeLock15ml => "eppendorf_safe_lock_1_5ml",
            Self::EppendorfSafeLock20ml => "eppendorf_safe_lock_2_0ml",
            Self::EppendorfPcrTube02ml => "eppendorf_pcr_tube_0_2ml",
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
    /// Callout on the profile's cited drawing page; None means unknown.
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
        TubeProfileId::EppendorfPcrTube02ml => return pcr_tube_profile(),
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

fn pcr_tube_profile() -> TubeProfile {
    let unknown = || DrawingDimension {
        value_mm: None,
        basis: DimensionBasis::NotDocumented,
        source_callout: None,
    };
    TubeProfile {
        schema: "gentle.tube_profile.v1".into(),
        id: TubeProfileId::EppendorfPcrTube02ml,
        name: "Eppendorf PCR Tube 0.2 ml".into(), nominal_volume_ml: 0.2,
        source_url: "https://www.eppendorf.com/product-media/doc/en/17413036/Consumables_Technical-data_PCR-Tubes_PCR-Tube-Strips_Cap-Strips.pdf".into(),
        source_sha256: "sha256:67b84dea22bcb2cbea26907b8f1116c506210f1917744feb64b2fafb6745b40d".into(),
        source_document: "0012 595 670-00/2026-05".into(), source_page_1based: 4,
        body_outer_diameter: DrawingDimension::measured(6.0, "maximum upper body/rim diameter 6; main cylindrical body diameter 5.8"),
        body_inner_diameter: DrawingDimension::measured(5.2, "upper inner diameter 5.2"),
        collar_outer_diameter: unknown(),
        body_length: DrawingDimension::measured(20.8, "overall body length with cap open"),
        cylindrical_section_length: DrawingDimension::measured(11.1, "top to outer conical transition"),
        closed_height: unknown(), closed_side_view_width: unknown(),
        open_90_height: unknown(), open_90_side_view_width: unknown(),
        cap_plan_outline_width: unknown(), bottom_shape: TubeBottomShape::ConicalWithRoundedTip,
        nonclaims: vec![
            "Drawing reference only; no physical fit, cap clearance or thermal/material qualification".into(),
            "No sufficient seating collar documented: maximum rim 6.0 versus main body 5.8 mm".into(),
            "22 mm is open-cap side-view width, not closed-cap or 90-degree geometry".into(),
        ],
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn labware_pcr_profile_and_96_carrier_preserve_unknowns_and_require_8x12() {
        let profile = tube_profile(TubeProfileId::EppendorfPcrTube02ml);
        assert_eq!(profile.nominal_volume_ml, 0.2);
        assert_eq!(profile.source_page_1based, 4);
        assert_eq!(profile.collar_outer_diameter.value_mm, None);
        assert_eq!(profile.closed_side_view_width.value_mm, None);
        assert_eq!(profile.body_outer_diameter.value_mm, Some(6.0));
        assert_eq!(profile.body_length.value_mm, Some(20.8));
        assert_eq!(profile.cylindrical_section_length.value_mm, Some(11.1));
        let geometry = derive_pcr96_carrier_geometry(&profile, 8, 12).unwrap();
        assert!((geometry.opening_diameter_mm - 6.6).abs() < 1e-10);
        assert_eq!(geometry.bore_depth_mm, 13.0);
        assert_eq!(geometry.pitch_mm, 9.0);
        for (rows, columns) in [(8, 11), (12, 8), (16, 24)] {
            assert!(
                derive_pcr96_carrier_geometry(&profile, rows, columns)
                    .unwrap_err()
                    .starts_with("carrier_requires_8x12:")
            );
        }
        assert!(
            derive_tube_carrier_geometry(&[profile.clone()], RackTubeSupportMode::FloorSeated)
                .unwrap_err()
                .starts_with("cap_envelope_unknown:")
        );
        let encoded = serde_json::to_string(&profile).unwrap();
        assert_eq!(
            serde_json::from_str::<TubeProfile>(&encoded).unwrap(),
            profile
        );
    }

    fn safe_lock_profiles() -> Vec<TubeProfile> {
        [
            TubeProfileId::EppendorfSafeLock15ml,
            TubeProfileId::EppendorfSafeLock20ml,
        ]
        .into_iter()
        .map(tube_profile)
        .collect()
    }

    #[test]
    fn labware_carriers_derive_both_support_modes_from_both_profiles() {
        let profiles = safe_lock_profiles();
        for (mode, opening, depth, edge) in [
            (RackTubeSupportMode::CollarSeated, 10.15, 42.0, 4.925),
            (RackTubeSupportMode::FloorSeated, 10.45, 22.0, 4.775),
        ] {
            let geometry = derive_tube_carrier_geometry(&profiles, mode).unwrap();
            assert!((geometry.opening_diameter_mm - opening).abs() < 1e-10);
            assert_eq!(geometry.bore_depth_mm, depth);
            assert_eq!(geometry.pitch_mm, 20.0);
            assert!((geometry.edge_margin_mm - edge).abs() < 1e-10);
            geometry.validate(&profiles, mode).unwrap();
            let mut reversed = profiles.clone();
            reversed.reverse();
            assert_eq!(
                derive_tube_carrier_geometry(&reversed, mode).unwrap(),
                geometry
            );
        }
    }

    #[test]
    fn labware_carrier_refuses_insufficient_collar_ledge() {
        let mut profiles = safe_lock_profiles();
        profiles[0].collar_outer_diameter.value_mm = Some(10.15);
        assert!(
            derive_tube_carrier_geometry(&profiles, RackTubeSupportMode::CollarSeated)
                .unwrap_err()
                .starts_with("collar_ledge_insufficient:")
        );
    }

    #[test]
    fn labware_carrier_refuses_cone_only_engagement() {
        let mut profiles = safe_lock_profiles();
        let geometry =
            derive_tube_carrier_geometry(&profiles, RackTubeSupportMode::FloorSeated).unwrap();
        // Synthetic longer cone makes this proposed bore unsafe.
        profiles[0].cylindrical_section_length.value_mm = Some(19.0);
        assert!(
            geometry
                .validate(&profiles, RackTubeSupportMode::FloorSeated)
                .unwrap_err()
                .starts_with("bore_engages_cone_only:")
        );
    }

    #[test]
    fn labware_carrier_refuses_collar_at_or_below_top() {
        let mut profiles = safe_lock_profiles();
        profiles[0].body_length.value_mm = Some(25.0);
        profiles[1].cylindrical_section_length.value_mm = Some(20.0);
        // The other synthetic profile requires a bore deeper than this collar permits.
        assert!(
            derive_tube_carrier_geometry(&profiles, RackTubeSupportMode::FloorSeated)
                .unwrap_err()
                .starts_with("collar_not_above_top:")
        );
    }

    #[test]
    fn labware_carrier_refuses_unknown_cap_envelope() {
        let mut profiles = safe_lock_profiles();
        profiles[0].closed_side_view_width.value_mm = None;
        assert!(
            derive_tube_carrier_geometry(&profiles, RackTubeSupportMode::FloorSeated)
                .unwrap_err()
                .starts_with("cap_envelope_unknown:")
        );
    }

    #[test]
    fn labware_carrier_refuses_invalid_and_unproven_dimensions() {
        assert!(derive_tube_carrier_geometry(&[], RackTubeSupportMode::FloorSeated).is_err());
        for value in [
            None,
            Some(f64::NAN),
            Some(f64::INFINITY),
            Some(0.0),
            Some(-1.0),
        ] {
            let mut profiles = safe_lock_profiles();
            profiles[0].body_outer_diameter.value_mm = value;
            assert!(
                derive_tube_carrier_geometry(&profiles, RackTubeSupportMode::FloorSeated).is_err()
            );
        }
        let mut profiles = safe_lock_profiles();
        profiles[0].body_outer_diameter.basis = DimensionBasis::NotDocumented;
        assert!(derive_tube_carrier_geometry(&profiles, RackTubeSupportMode::FloorSeated).is_err());
    }

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
