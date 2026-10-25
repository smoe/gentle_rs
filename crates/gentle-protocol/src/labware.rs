//! Source-backed visualization profiles and bounded rack-grid admission.
//!
//! Profiles describe drawing dimensions, not physical inventory, fit or material
//! qualification. Admission is explicit at processing boundaries, not deserialization.

use serde::{Deserialize, Serialize};

/// Resource budget for processing a rectangular grid, not a carrier capacity claim.
pub const MAX_RACK_PROCESSING_POSITIONS: usize = 1_536;
/// A cut PCR strip retains one through eight contiguous tubes.
pub const MAX_PCR_STRIP_TUBES: usize = 8;

/// Engineering assumptions in millimetres, not manufacturing tolerances.
pub const RADIAL_CLEARANCE_COLLAR: f64 = 0.15;
pub const RADIAL_CLEARANCE_FLOOR: f64 = 0.30;
pub const MIN_COLLAR_LEDGE: f64 = 0.20;
pub const TIP_CLEARANCE: f64 = 2.0;
pub const CYLINDER_ENGAGEMENT: f64 = 3.0;
pub const COLLAR_TOP_MARGIN: f64 = 5.0;
/// Compact-layout engineering allowances, not measured strength/cap clearance.
pub const MIN_WEB: f64 = 0.8;
pub const COLLAR_GAP: f64 = 0.5;
/// Derived storage-box engineering assumptions, not measured fit or tolerances.
pub const RADIAL_CLEARANCE_SNUG: f64 = 0.15;
pub const BOX_DIVIDER_THICKNESS: f64 = 1.2;
pub const BOX_WALL_THICKNESS: f64 = 2.0;
pub const BOX_LID_CLEARANCE: f64 = 3.0;
pub const LID_FIT_CLEARANCE: f64 = 0.3;
pub const LID_TOP_THICKNESS: f64 = 2.0;
pub const LID_SKIRT_HEIGHT: f64 = 8.0;

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
/// GENtle-derived square-cell box; millimetres, never commercial-box or physical-fit evidence.
pub struct StorageBoxGeometry {
    pub rows: usize,
    pub columns: usize,
    pub cell_width_mm: f64,
    pub divider_thickness_mm: f64,
    pub pitch_mm: f64,
    pub divider_height_mm: f64,
    pub wall_thickness_mm: f64,
    pub interior_width_mm: f64,
    pub interior_depth_mm: f64,
    pub interior_height_mm: f64,
    pub exterior_width_mm: f64,
    pub exterior_depth_mm: f64,
    pub exterior_height_mm: f64,
    pub lid_fit_clearance_mm: f64,
    pub lid_top_thickness_mm: f64,
    pub lid_skirt_height_mm: f64,
    pub lid_inner_width_mm: f64,
    pub lid_inner_depth_mm: f64,
    pub lid_exterior_width_mm: f64,
    pub lid_exterior_depth_mm: f64,
    pub lid_exterior_height_mm: f64,
}

/// Admit the grid first, then derive a floor-supported box from the worst profile dimensions.
pub fn derive_storage_box_geometry(
    profiles: &[TubeProfile],
    rows: usize,
    columns: usize,
) -> Result<StorageBoxGeometry, String> {
    admit_rack_grid(rows, columns)?;
    let carrier = derive_tube_carrier_geometry(profiles, RackTubeSupportMode::FloorSeated)?;
    let mut body: f64 = 0.0;
    let mut closed: f64 = 0.0;
    let mut collar: f64 = 0.0;
    for profile in profiles {
        body = body.max(required_dimension(
            &profile.body_outer_diameter,
            "body outer diameter",
        )?);
        collar = collar.max(required_dimension(
            &profile.collar_outer_diameter,
            "collar diameter",
        )?);
        closed = closed.max(
            required_dimension(&profile.closed_height, "closed height").map_err(|_| {
                "cap_height_unknown: storage box requires a documented closed height".to_string()
            })?,
        );
    }
    let cell = body + 2.0 * RADIAL_CLEARANCE_SNUG;
    let pitch = cell + BOX_DIVIDER_THICKNESS;
    if pitch < collar + COLLAR_GAP {
        return Err("box_collar_overlap: adjacent collars need a 0.5 mm gap".into());
    }
    let width = columns as f64 * cell + (columns - 1) as f64 * BOX_DIVIDER_THICKNESS;
    let depth = rows as f64 * cell + (rows - 1) as f64 * BOX_DIVIDER_THICKNESS;
    let exterior_width = width + 2.0 * BOX_WALL_THICKNESS;
    let exterior_depth = depth + 2.0 * BOX_WALL_THICKNESS;
    let lid_inner_width = exterior_width + 2.0 * LID_FIT_CLEARANCE;
    let lid_inner_depth = exterior_depth + 2.0 * LID_FIT_CLEARANCE;
    Ok(StorageBoxGeometry {
        rows,
        columns,
        cell_width_mm: cell,
        divider_thickness_mm: BOX_DIVIDER_THICKNESS,
        pitch_mm: pitch,
        divider_height_mm: carrier.bore_depth_mm,
        wall_thickness_mm: BOX_WALL_THICKNESS,
        interior_width_mm: width,
        interior_depth_mm: depth,
        interior_height_mm: closed + BOX_LID_CLEARANCE,
        exterior_width_mm: exterior_width,
        exterior_depth_mm: exterior_depth,
        exterior_height_mm: BOX_WALL_THICKNESS + closed + BOX_LID_CLEARANCE,
        lid_fit_clearance_mm: LID_FIT_CLEARANCE,
        lid_top_thickness_mm: LID_TOP_THICKNESS,
        lid_skirt_height_mm: LID_SKIRT_HEIGHT,
        lid_inner_width_mm: lid_inner_width,
        lid_inner_depth_mm: lid_inner_depth,
        lid_exterior_width_mm: lid_inner_width + 2.0 * BOX_WALL_THICKNESS,
        lid_exterior_depth_mm: lid_inner_depth + 2.0 * BOX_WALL_THICKNESS,
        lid_exterior_height_mm: LID_TOP_THICKNESS + LID_SKIRT_HEIGHT,
    })
}

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

/// Explicit compact layout only. Preserve conservative outer margins; caps and
/// hinges can collide because their plan-view envelope has not been measured.
pub fn derive_compact_tube_carrier_geometry(
    profiles: &[TubeProfile],
    mode: RackTubeSupportMode,
    pitch_mm: f64,
) -> Result<TubeCarrierGeometry, String> {
    let mut geometry = derive_tube_carrier_geometry(profiles, mode)?;
    if !pitch_mm.is_finite() || pitch_mm <= 0.0 {
        return Err("compact_pitch_invalid: pitch must be finite and positive".into());
    }
    if pitch_mm > geometry.pitch_mm {
        return Err("compact_pitch_not_compact: pitch exceeds conservative spacing".into());
    }
    if pitch_mm < geometry.opening_diameter_mm + 2.0 * MIN_WEB {
        return Err(
            "compact_pitch_below_web: pitch needs two 0.8 mm web allowances beyond the opening"
                .into(),
        );
    }
    for profile in profiles {
        let collar = required_dimension(&profile.collar_outer_diameter, "collar diameter")?;
        if pitch_mm < collar + COLLAR_GAP {
            return Err("compact_pitch_collar_overlap: adjacent collars need a 0.5 mm gap".into());
        }
    }
    geometry.pitch_mm = pitch_mm;
    geometry.validate(profiles, mode)?;
    Ok(geometry)
}

/// Exact positive decimal millimetres to whole micrometres. No exponent,
/// rounding, signs or more than three fractional places are accepted.
pub fn parse_carrier_pitch_um(value: &str) -> Result<u32, String> {
    let invalid = || {
        "compact_pitch_invalid: use positive decimal mm with at most three fractional digits"
            .to_string()
    };
    let value = value.trim();
    let (whole, fraction) = value.split_once('.').unwrap_or((value, ""));
    if whole.is_empty()
        || !whole.bytes().all(|b| b.is_ascii_digit())
        || fraction.len() > 3
        || !fraction.bytes().all(|b| b.is_ascii_digit())
        || (value.contains('.') && fraction.is_empty())
    {
        return Err(invalid());
    }
    let whole = whole.parse::<u32>().map_err(|_| invalid())?;
    let fractional = if fraction.is_empty() {
        0
    } else {
        fraction.parse::<u32>().map_err(|_| invalid())? * 10_u32.pow(3 - fraction.len() as u32)
    };
    whole
        .checked_mul(1000)
        .and_then(|v| v.checked_add(fractional))
        .filter(|v| *v > 0)
        .ok_or_else(invalid)
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
    #[serde(rename = "brand_pcr_strip_0_2ml_781320")]
    BrandPcrStrip02ml781320,
}

impl TubeProfileId {
    /// Stable inspection order; new profiles are appended.
    pub const ALL: [Self; 5] = [
        Self::EppendorfSafeLock05ml,
        Self::EppendorfSafeLock15ml,
        Self::EppendorfSafeLock20ml,
        Self::EppendorfPcrTube02ml,
        Self::BrandPcrStrip02ml781320,
    ];

    /// The same stable identity used by the JSON contract.
    pub fn as_str(self) -> &'static str {
        match self {
            Self::EppendorfSafeLock05ml => "eppendorf_safe_lock_0_5ml",
            Self::EppendorfSafeLock15ml => "eppendorf_safe_lock_1_5ml",
            Self::EppendorfSafeLock20ml => "eppendorf_safe_lock_2_0ml",
            Self::EppendorfPcrTube02ml => "eppendorf_pcr_tube_0_2ml",
            Self::BrandPcrStrip02ml781320 => "brand_pcr_strip_0_2ml_781320",
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
/// Drawing-backed strip layout, separate from illustrative bridge/transition geometry.
pub struct TubeStripLayout {
    pub max_tubes: usize,
    pub center_pitch: DrawingDimension,
    pub rim_outer_diameter: DrawingDimension,
    pub tip_outer_diameter: DrawingDimension,
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
    /// Absent for individual tubes; old profile payloads/bytes stay unchanged.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub strip: Option<TubeStripLayout>,
    pub nonclaims: Vec<String>,
}

impl TubeProfile {
    /// Omission selects the full strip or one individual tube. Never clamp input.
    pub fn model_tube_count(&self, requested: Option<usize>) -> Result<usize, String> {
        let maximum = self.strip.as_ref().map_or(1, |strip| strip.max_tubes);
        if maximum == 0 || maximum > MAX_PCR_STRIP_TUBES {
            return Err("tube_model_strip_invalid: maximum must be in 1..=8".into());
        }
        let count = requested.unwrap_or(maximum);
        if count == 0 || count > maximum {
            return Err(format!(
                "tube_model_count_invalid: {} requires 1..={maximum} tubes, got {count}",
                self.id.as_str()
            ));
        }
        Ok(count)
    }
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
        TubeProfileId::BrandPcrStrip02ml781320 => return brand_pcr_strip_profile(),
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
        strip: None,
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
        strip: None,
        nonclaims: vec![
            "Drawing reference only; no physical fit, cap clearance or thermal/material qualification".into(),
            "No sufficient seating collar documented: maximum rim 6.0 versus main body 5.8 mm".into(),
            "22 mm is open-cap side-view width, not closed-cap or 90-degree geometry".into(),
        ],
    }
}

fn brand_pcr_strip_profile() -> TubeProfile {
    let unknown = || DrawingDimension {
        value_mm: None,
        basis: DimensionBasis::NotDocumented,
        source_callout: None,
    };
    TubeProfile {
        schema: "gentle.tube_profile.v1".into(),
        id: TubeProfileId::BrandPcrStrip02ml781320,
        name: "BRAND PCR strip 0.2 ml - 1-8 tubes (781320)".into(),
        nominal_volume_ml: 0.2,
        source_url: "https://shop.brand.de/media/import/1/27/32406/42485/53752/56329/55264/Data_sheet_PCR_strips_DE.pdf".into(),
        source_sha256: "sha256:8ba466736960b6933c642568d03ac5901c09010db176b7278d6e2989ffbd3c13".into(),
        source_document: "BRAND PCR strips datasheet, article 781320; undated, retrieved 2026-10-22".into(),
        source_page_1based: 6,
        body_outer_diameter: DrawingDimension::measured(6.0, "section A-A: main body diameter 6"),
        body_inner_diameter: unknown(),
        collar_outer_diameter: DrawingDimension::measured(8.0, "section A-A: outer collar diameter 8"),
        body_length: DrawingDimension::measured(20.75, "section A-A: overall uncapped height 20.75"),
        cylindrical_section_length: unknown(),
        closed_height: unknown(),
        closed_side_view_width: unknown(),
        open_90_height: unknown(),
        open_90_side_view_width: unknown(),
        cap_plan_outline_width: unknown(),
        bottom_shape: TubeBottomShape::ConicalWithRoundedTip,
        strip: Some(TubeStripLayout {
            max_tubes: MAX_PCR_STRIP_TUBES,
            center_pitch: DrawingDimension::measured(9.0, "adjacent centers 9; first-to-last centers 63 for eight tubes"),
            rim_outer_diameter: DrawingDimension::measured(7.05, "section A-A: upper rim diameter 7.05"),
            tip_outer_diameter: DrawingDimension::measured(2.75, "section A-A: lower outer diameter 2.75 before rounded tip"),
        }),
        nonclaims: vec![
            "0.2 ml is nominal capacity per tube, not measured model cavity volume".into(),
            "Visualization only; no physical fit, manufacturing, material or temperature qualification".into(),
            "Separate caps, inner cavities, numbering and end tab/cut stubs are omitted".into(),
            "Bridge cross-section, collar position/thickness, taper transition and tip rounding are illustrative, not drawing dimensions".into(),
            "Counts 1-8 represent user-selected contiguous cut lengths, not separate manufacturer articles".into(),
        ],
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn labware_storage_box_geometry_has_square_cells_variable_grids_and_separate_lid() {
        let profiles = [
            tube_profile(TubeProfileId::EppendorfSafeLock15ml),
            tube_profile(TubeProfileId::EppendorfSafeLock20ml),
        ];
        for (rows, columns, width, depth) in [
            (9, 9, 100.95, 100.95),
            (1, 1, 10.15, 10.15),
            (4, 12, 135.0, 44.2),
        ] {
            let g = derive_storage_box_geometry(&profiles, rows, columns).unwrap();
            for (actual, expected) in [
                (g.cell_width_mm, 10.15),
                (g.divider_thickness_mm, 1.2),
                (g.pitch_mm, 11.35),
                (g.divider_height_mm, 22.0),
                (g.interior_height_mm, 45.0),
                (g.interior_width_mm, width),
                (g.interior_depth_mm, depth),
                (g.exterior_width_mm, width + 4.0),
                (g.exterior_depth_mm, depth + 4.0),
                (g.exterior_height_mm, 47.0),
                (g.lid_inner_width_mm, width + 4.6),
                (g.lid_inner_depth_mm, depth + 4.6),
                (g.lid_exterior_width_mm, width + 8.6),
                (g.lid_exterior_depth_mm, depth + 8.6),
                (g.lid_exterior_height_mm, 10.0),
            ] {
                assert!((actual - expected).abs() < 1e-9, "{actual} != {expected}");
            }
            assert_eq!(g.rows, rows);
            assert_eq!(g.columns, columns);
            let roundtrip: StorageBoxGeometry =
                serde_json::from_str(&serde_json::to_string(&g).unwrap()).unwrap();
            assert_eq!(roundtrip, g);
        }
    }

    #[test]
    fn labware_storage_box_refuses_invalid_grids_unknown_caps_and_colliding_collars() {
        // Hand-mutated drawing profiles are synthetic invalid geometry, not vendor data.
        let base = tube_profile(TubeProfileId::EppendorfSafeLock20ml);
        let mut profile = base.clone();
        profile.collar_outer_diameter.value_mm = Some(11.0);
        assert!(
            derive_storage_box_geometry(&[profile], 9, 9)
                .unwrap_err()
                .starts_with("box_collar_overlap:")
        );
        let mut profile = base.clone();
        profile.closed_height.value_mm = None;
        assert!(
            derive_storage_box_geometry(&[profile], 9, 9)
                .unwrap_err()
                .starts_with("cap_height_unknown:")
        );
        let mut profile = base;
        profile.cylindrical_section_length.value_mm = Some(1.0);
        let expected =
            derive_tube_carrier_geometry(&[profile.clone()], RackTubeSupportMode::FloorSeated)
                .unwrap_err();
        assert!(expected.starts_with("collar_not_above_top:"));
        assert_eq!(
            derive_storage_box_geometry(&[profile], 9, 9).unwrap_err(),
            expected
        );
        assert!(
            derive_storage_box_geometry(&[], 0, 5)
                .unwrap_err()
                .starts_with("rack_grid_empty:")
        );
        assert!(
            derive_storage_box_geometry(&[], 1, 1537)
                .unwrap_err()
                .starts_with("rack_grid_limit:")
        );
    }

    #[test]
    fn labware_pcr_strip_count_bounds_and_source_dimensions() {
        let profile = tube_profile(TubeProfileId::BrandPcrStrip02ml781320);
        let strip = profile.strip.as_ref().unwrap();
        assert_eq!(profile.nominal_volume_ml, 0.2);
        assert_eq!(profile.source_page_1based, 6);
        assert_eq!(strip.max_tubes, 8);
        assert_eq!(strip.center_pitch.value_mm, Some(9.0));
        assert_eq!(profile.body_length.value_mm, Some(20.75));
        assert_eq!(profile.body_outer_diameter.value_mm, Some(6.0));
        assert_eq!(profile.collar_outer_diameter.value_mm, Some(8.0));
        assert_eq!(strip.rim_outer_diameter.value_mm, Some(7.05));
        assert_eq!(strip.tip_outer_diameter.value_mm, Some(2.75));
        assert_eq!(profile.cylindrical_section_length.value_mm, None);
        assert_eq!(profile.closed_height.value_mm, None);
        assert_eq!(profile.model_tube_count(None).unwrap(), 8);
        for count in 1..=8 {
            assert_eq!(profile.model_tube_count(Some(count)).unwrap(), count);
        }
        for count in [0, 9, usize::MAX] {
            assert!(profile.model_tube_count(Some(count)).is_err());
        }
        assert_eq!(
            serde_json::from_value::<TubeProfile>(serde_json::to_value(&profile).unwrap()).unwrap(),
            profile
        );
        let single = tube_profile(TubeProfileId::EppendorfSafeLock15ml);
        assert_eq!(single.model_tube_count(None).unwrap(), 1);
        assert!(single.model_tube_count(Some(2)).is_err());
        let old = serde_json::to_value(&single).unwrap();
        assert!(old.get("strip").is_none());
        assert_eq!(serde_json::from_value::<TubeProfile>(old).unwrap(), single);
    }

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
    fn labware_compact_pitch_bounds_preserve_conservative_edges() {
        let profiles = safe_lock_profiles();
        for (mode, minimum, edge) in [
            (RackTubeSupportMode::CollarSeated, 11.75, 4.925),
            (RackTubeSupportMode::FloorSeated, 12.05, 4.775),
        ] {
            for pitch in [minimum, 14.0, 20.0] {
                let g = derive_compact_tube_carrier_geometry(&profiles, mode, pitch).unwrap();
                assert_eq!(g.pitch_mm, pitch);
                assert!((g.edge_margin_mm - edge).abs() < 1e-10);
            }
            for (pitch, code) in [
                (minimum - 0.001, "compact_pitch_below_web:"),
                (20.001, "compact_pitch_not_compact:"),
                (f64::NAN, "compact_pitch_invalid:"),
                (f64::INFINITY, "compact_pitch_invalid:"),
                (0.0, "compact_pitch_invalid:"),
                (-1.0, "compact_pitch_invalid:"),
            ] {
                assert!(
                    derive_compact_tube_carrier_geometry(&profiles, mode, pitch)
                        .unwrap_err()
                        .starts_with(code)
                );
            }
        }
        let mut synthetic = profiles;
        synthetic[0].collar_outer_diameter.value_mm = Some(12.0);
        assert!(
            derive_compact_tube_carrier_geometry(
                &synthetic,
                RackTubeSupportMode::CollarSeated,
                12.1
            )
            .unwrap_err()
            .starts_with("compact_pitch_collar_overlap:")
        );
    }

    #[test]
    fn labware_compact_pitch_decimal_and_legacy_rack_bytes_are_exact() {
        for (mm, um) in [
            ("14", 14000),
            ("14.000", 14000),
            ("11.751", 11751),
            ("0.001", 1),
        ] {
            assert_eq!(parse_carrier_pitch_um(mm).unwrap(), um);
        }
        for mm in [
            "", "0", "-1", "+14", "14.0005", "NaN", "inf", "1e1", "14.", ".5", "1.2.3", "4294968",
        ] {
            assert!(parse_carrier_pitch_um(mm).is_err(), "{mm}");
        }
        let mut rack = crate::Rack::default();
        let legacy = serde_json::to_string(&rack).unwrap();
        assert!(!legacy.contains("carrier_pitch_override_um"));
        assert_eq!(
            serde_json::from_str::<crate::Rack>(&legacy)
                .unwrap()
                .carrier_pitch_override_um,
            None
        );
        rack.carrier_pitch_override_um = Some(14_000);
        assert_eq!(
            serde_json::from_str::<crate::Rack>(&serde_json::to_string(&rack).unwrap()).unwrap(),
            rack
        );
        rack.carrier_pitch_override_um = None;
        assert_eq!(serde_json::to_string(&rack).unwrap(), legacy);
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
