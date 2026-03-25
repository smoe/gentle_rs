//! Protocol-cartoon catalog and deterministic SVG rendering.
//!
//! This module defines a reusable abstraction for protocol cartoons:
//! DNA molecules (linear or circular) are composed of colored feature fragments,
//! optional linear-end styles (hidden, continuation, blunt, or sticky),
//! optional per-strand feature colors, and protocol events arranged as a
//! sequence of snapshots.

use serde::{Deserialize, Serialize};

/// Stable protocol-cartoon identifiers exposed through engine operations.
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq, Hash)]
pub enum ProtocolCartoonKind {
    #[serde(rename = "gibson.two_fragment")]
    GibsonTwoFragment,
    #[serde(rename = "gibson.single_insert_dual_junction")]
    GibsonSingleInsertDualJunction,
    #[serde(rename = "pcr.assay.pair")]
    PcrAssayPair,
    #[serde(rename = "pcr.assay.pair.no_product")]
    PcrAssayPairNoProduct,
    #[serde(rename = "pcr.assay.qpcr")]
    PcrAssayQpcr,
}

impl ProtocolCartoonKind {
    /// Canonical string identifier used in shell/CLI input and output rows.
    pub fn id(&self) -> &'static str {
        match self {
            Self::GibsonTwoFragment => "gibson.two_fragment",
            Self::GibsonSingleInsertDualJunction => "gibson.single_insert_dual_junction",
            Self::PcrAssayPair => "pcr.assay.pair",
            Self::PcrAssayPairNoProduct => "pcr.assay.pair.no_product",
            Self::PcrAssayQpcr => "pcr.assay.qpcr",
        }
    }

    /// Human-readable title for UI/help surfaces.
    pub fn title(&self) -> &'static str {
        match self {
            Self::GibsonTwoFragment => "Gibson Assembly (two-fragment conceptual)",
            Self::GibsonSingleInsertDualJunction => {
                "Gibson Assembly (single-insert dual-junction mechanism)"
            }
            Self::PcrAssayPair => "PCR Assay (pair-primer baseline)",
            Self::PcrAssayPairNoProduct => "PCR Assay (report-only no-product)",
            Self::PcrAssayQpcr => "qPCR Assay (probe-bearing baseline)",
        }
    }

    /// Short semantics summary used in list outputs.
    pub fn summary(&self) -> &'static str {
        match self {
            Self::GibsonTwoFragment => {
                "Event-sequence cartoon with continuation/sticky/blunt ends and strand-separated DNA glyphs"
            }
            Self::GibsonSingleInsertDualJunction => {
                "Single-insert Gibson cartoon showing both destination-insert junctions explicitly"
            }
            Self::PcrAssayPair => {
                "Mechanism-first pair-PCR strip: template context, ROI, assay setup, amplification, and amplicon/report outcome"
            }
            Self::PcrAssayPairNoProduct => {
                "Pair-PCR report-only strip showing a selected ROI and failed/no-product assay outcome without literal primer glyphs"
            }
            Self::PcrAssayQpcr => {
                "Mechanism-first qPCR strip: template context, ROI, probe-bearing assay setup, amplification, and quantitative readout"
            }
        }
    }

    /// Parse supported shell/CLI aliases into the canonical enum.
    pub fn parse_id(raw: &str) -> Option<Self> {
        match raw.trim().to_ascii_lowercase().as_str() {
            "gibson.two_fragment" | "gibson.two-fragment" | "gibson_two_fragment" | "gibson" => {
                Some(Self::GibsonTwoFragment)
            }
            "gibson.single_insert_dual_junction"
            | "gibson.single-insert-dual-junction"
            | "gibson_single_insert_dual_junction"
            | "gibson.single_insert"
            | "gibson.destination_first_single_insert" => {
                Some(Self::GibsonSingleInsertDualJunction)
            }
            "pcr.assay.pair" | "pcr.assay" | "pcr.pair" | "pcr_pair" | "pcr" => {
                Some(Self::PcrAssayPair)
            }
            "pcr.assay.pair.no_product"
            | "pcr.assay.pair.report_only"
            | "pcr.assay.report_only"
            | "pcr_pair_no_product" => Some(Self::PcrAssayPairNoProduct),
            "pcr.assay.qpcr" | "pcr.qpcr" | "qpcr" | "q-pcr" => Some(Self::PcrAssayQpcr),
            _ => None,
        }
    }

    /// Deterministic ordered catalog for list commands.
    pub fn catalog() -> Vec<Self> {
        vec![
            Self::GibsonTwoFragment,
            Self::GibsonSingleInsertDualJunction,
            Self::PcrAssayPair,
            Self::PcrAssayPairNoProduct,
            Self::PcrAssayQpcr,
        ]
    }
}

/// Machine-readable row used by list commands.
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq)]
pub struct ProtocolCartoonCatalogRow {
    pub id: String,
    pub title: String,
    pub summary: String,
}

/// JSON-facing protocol-cartoon template schema (v1).
///
/// A template can be sparse: omitted event/molecule/feature fields are
/// expanded deterministically through `defaults`.
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq)]
pub struct ProtocolCartoonTemplate {
    #[serde(default = "protocol_cartoon_template_schema_id")]
    pub schema: String,
    pub id: String,
    #[serde(default)]
    pub title: String,
    #[serde(default)]
    pub summary: String,
    #[serde(default)]
    pub defaults: ProtocolCartoonTemplateDefaults,
    #[serde(default)]
    pub events: Vec<ProtocolCartoonTemplateEvent>,
}

/// Default values used when template events/molecules/features omit fields.
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq)]
pub struct ProtocolCartoonTemplateDefaults {
    #[serde(default = "default_template_event_caption")]
    pub event_caption: String,
    #[serde(default = "default_template_event_action")]
    pub event_action: ProtocolCartoonAction,
    #[serde(default = "default_template_molecule_topology")]
    pub molecule_topology: DnaTopologyCartoon,
    #[serde(default = "default_template_linear_end_style")]
    pub linear_end_style: DnaEndStyle,
    #[serde(default = "default_template_feature_length_bp")]
    pub feature_length_bp: usize,
    #[serde(default = "default_template_palette")]
    pub palette: Vec<String>,
}

impl Default for ProtocolCartoonTemplateDefaults {
    fn default() -> Self {
        Self {
            event_caption: default_template_event_caption(),
            event_action: default_template_event_action(),
            molecule_topology: default_template_molecule_topology(),
            linear_end_style: default_template_linear_end_style(),
            feature_length_bp: default_template_feature_length_bp(),
            palette: default_template_palette(),
        }
    }
}

/// One event row in a protocol-cartoon template.
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq, Default)]
pub struct ProtocolCartoonTemplateEvent {
    #[serde(default)]
    pub id: String,
    #[serde(default)]
    pub title: String,
    #[serde(default)]
    pub caption: String,
    #[serde(default)]
    pub action: Option<ProtocolCartoonAction>,
    #[serde(default)]
    pub molecules: Vec<ProtocolCartoonTemplateMolecule>,
}

/// One molecule row in a protocol-cartoon template.
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq, Default)]
pub struct ProtocolCartoonTemplateMolecule {
    #[serde(default)]
    pub id: String,
    #[serde(default)]
    pub label: String,
    #[serde(default)]
    pub topology: Option<DnaTopologyCartoon>,
    #[serde(default)]
    pub features: Vec<ProtocolCartoonTemplateFeature>,
    #[serde(default)]
    pub left_end: Option<DnaEndStyle>,
    #[serde(default)]
    pub right_end: Option<DnaEndStyle>,
}

/// One feature-fragment row in a protocol-cartoon template.
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq, Default)]
pub struct ProtocolCartoonTemplateFeature {
    #[serde(default)]
    pub id: String,
    #[serde(default)]
    pub label: String,
    #[serde(default)]
    pub length_bp: Option<usize>,
    #[serde(default)]
    pub top_length_bp: Option<usize>,
    #[serde(default)]
    pub bottom_length_bp: Option<usize>,
    #[serde(default)]
    pub color_hex: Option<String>,
    #[serde(default)]
    pub bottom_color_hex: Option<String>,
    #[serde(default)]
    pub top_nick_after: Option<bool>,
    #[serde(default)]
    pub bottom_nick_after: Option<bool>,
}

/// JSON-facing protocol-cartoon template binding schema (v1).
///
/// Bindings provide deterministic override values (for example protocol-derived
/// overlap lengths/end styles) without modifying the source template file.
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq)]
pub struct ProtocolCartoonTemplateBindings {
    #[serde(default = "protocol_cartoon_template_bindings_schema_id")]
    pub schema: String,
    #[serde(default)]
    pub template_id: Option<String>,
    #[serde(default)]
    pub defaults: ProtocolCartoonTemplateBindingsDefaults,
    #[serde(default)]
    pub event_overrides: Vec<ProtocolCartoonTemplateEventBinding>,
    #[serde(default)]
    pub molecule_overrides: Vec<ProtocolCartoonTemplateMoleculeBinding>,
    #[serde(default)]
    pub feature_overrides: Vec<ProtocolCartoonTemplateFeatureBinding>,
}

/// Optional default overrides for template resolution.
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq, Default)]
pub struct ProtocolCartoonTemplateBindingsDefaults {
    #[serde(default)]
    pub event_caption: Option<String>,
    #[serde(default)]
    pub event_action: Option<ProtocolCartoonAction>,
    #[serde(default)]
    pub molecule_topology: Option<DnaTopologyCartoon>,
    #[serde(default)]
    pub linear_end_style: Option<DnaEndStyle>,
    #[serde(default)]
    pub feature_length_bp: Option<usize>,
    #[serde(default)]
    pub palette: Option<Vec<String>>,
}

/// Event-level binding override.
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq)]
pub struct ProtocolCartoonTemplateEventBinding {
    pub event_id: String,
    #[serde(default)]
    pub title: Option<String>,
    #[serde(default)]
    pub caption: Option<String>,
    #[serde(default)]
    pub action: Option<ProtocolCartoonAction>,
}

/// Molecule-level binding override.
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq)]
pub struct ProtocolCartoonTemplateMoleculeBinding {
    pub event_id: String,
    pub molecule_id: String,
    #[serde(default)]
    pub label: Option<String>,
    #[serde(default)]
    pub topology: Option<DnaTopologyCartoon>,
    #[serde(default)]
    pub left_end: Option<DnaEndStyle>,
    #[serde(default)]
    pub right_end: Option<DnaEndStyle>,
}

/// Feature-level binding override.
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq)]
pub struct ProtocolCartoonTemplateFeatureBinding {
    pub event_id: String,
    pub molecule_id: String,
    pub feature_id: String,
    #[serde(default)]
    pub label: Option<String>,
    #[serde(default)]
    pub length_bp: Option<usize>,
    #[serde(default)]
    pub top_length_bp: Option<usize>,
    #[serde(default)]
    pub bottom_length_bp: Option<usize>,
    #[serde(default)]
    pub color_hex: Option<String>,
    #[serde(default)]
    pub bottom_color_hex: Option<String>,
    #[serde(default)]
    pub top_nick_after: Option<bool>,
    #[serde(default)]
    pub bottom_nick_after: Option<bool>,
}

/// One feature fragment shown as a colored span in a molecule cartoon.
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq)]
pub struct DnaFeatureCartoon {
    pub id: String,
    pub label: String,
    pub length_bp: usize,
    pub top_length_bp: usize,
    pub bottom_length_bp: usize,
    pub color_hex: String,
    pub bottom_color_hex: String,
    pub top_nick_after: bool,
    pub bottom_nick_after: bool,
}

/// Cartoon topology for a DNA molecule glyph.
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq)]
pub enum DnaTopologyCartoon {
    Linear,
    Circular,
}

/// End overhang polarity relative to the top strand rendered left-to-right.
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq)]
pub enum OverhangPolarity {
    FivePrime,
    ThreePrime,
}

impl OverhangPolarity {
    fn label(&self) -> &'static str {
        match self {
            Self::FivePrime => "5'",
            Self::ThreePrime => "3'",
        }
    }
}

/// End style for a linear DNA cartoon.
///
/// `Continuation` means the molecule continues off-panel and should render as a
/// truncated/ruptured edge rather than a biologically meaningful terminus.
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq)]
pub enum DnaEndStyle {
    NotShown,
    Continuation,
    Blunt,
    Sticky {
        polarity: OverhangPolarity,
        nt: usize,
    },
}

/// One DNA molecule drawn in an event panel.
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq)]
pub struct DnaMoleculeCartoon {
    pub id: String,
    pub label: String,
    pub topology: DnaTopologyCartoon,
    pub features: Vec<DnaFeatureCartoon>,
    pub left_end: Option<DnaEndStyle>,
    pub right_end: Option<DnaEndStyle>,
}

/// Semantic action type of an event panel.
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq)]
pub enum ProtocolCartoonAction {
    Context,
    Cut,
    Separate,
    Anneal,
    Extend,
    Seal,
    Ligate,
    Custom { label: String },
}

impl ProtocolCartoonAction {
    fn label(&self) -> &str {
        match self {
            Self::Context => "Context",
            Self::Cut => "Cut",
            Self::Separate => "Separate",
            Self::Anneal => "Anneal",
            Self::Extend => "Extend",
            Self::Seal => "Seal",
            Self::Ligate => "Ligate",
            Self::Custom { label } => label,
        }
    }
}

/// One protocol event snapshot in the strip.
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq)]
pub struct ProtocolCartoonEvent {
    pub id: String,
    pub title: String,
    pub caption: String,
    pub action: ProtocolCartoonAction,
    pub molecules: Vec<DnaMoleculeCartoon>,
}

/// Full protocol cartoon specification.
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq)]
pub struct ProtocolCartoonSpec {
    pub id: String,
    pub title: String,
    pub summary: String,
    pub events: Vec<ProtocolCartoonEvent>,
}

impl ProtocolCartoonSpec {
    fn validate(&self) -> Result<(), String> {
        if self.id.trim().is_empty() {
            return Err("Protocol cartoon id must be non-empty".to_string());
        }
        if self.events.is_empty() {
            return Err("Protocol cartoon must contain at least one event".to_string());
        }
        for (event_idx, event) in self.events.iter().enumerate() {
            if event.id.trim().is_empty() {
                return Err(format!(
                    "Event {} in protocol '{}' has empty id",
                    event_idx + 1,
                    self.id
                ));
            }
            if event.molecules.is_empty() {
                return Err(format!(
                    "Event '{}' in protocol '{}' must contain at least one molecule",
                    event.id, self.id
                ));
            }
            for molecule in &event.molecules {
                if molecule.id.trim().is_empty() {
                    return Err(format!(
                        "Event '{}' in protocol '{}' contains a molecule with empty id",
                        event.id, self.id
                    ));
                }
                if molecule.features.is_empty() {
                    return Err(format!(
                        "Molecule '{}' in event '{}' has no feature fragments",
                        molecule.id, event.id
                    ));
                }
                if molecule.topology == DnaTopologyCartoon::Circular
                    && (molecule.left_end.is_some() || molecule.right_end.is_some())
                {
                    return Err(format!(
                        "Circular molecule '{}' in event '{}' cannot declare linear end styles",
                        molecule.id, event.id
                    ));
                }
                for (side, end) in [
                    ("left", molecule.left_end.as_ref()),
                    ("right", molecule.right_end.as_ref()),
                ] {
                    if let Some(DnaEndStyle::Sticky { nt, .. }) = end {
                        if *nt == 0 {
                            return Err(format!(
                                "Molecule '{}' in event '{}' has a {} sticky end with zero nt overhang",
                                molecule.id, event.id, side
                            ));
                        }
                    }
                }
                for feature in &molecule.features {
                    if feature.top_length_bp == 0 && feature.bottom_length_bp == 0 {
                        return Err(format!(
                            "Feature '{}' in molecule '{}' has zero top/bottom strand length",
                            feature.id, molecule.id
                        ));
                    }
                    if feature.top_length_bp > feature.length_bp
                        || feature.bottom_length_bp > feature.length_bp
                    {
                        return Err(format!(
                            "Feature '{}' in molecule '{}' exceeds canonical feature length on one strand",
                            feature.id, molecule.id
                        ));
                    }
                }
                let top_total_bp: usize = molecule.features.iter().map(|f| f.top_length_bp).sum();
                let bottom_total_bp: usize =
                    molecule.features.iter().map(|f| f.bottom_length_bp).sum();
                if top_total_bp == 0 || bottom_total_bp == 0 {
                    return Err(format!(
                        "Molecule '{}' in event '{}' must retain non-zero top and bottom strand length",
                        molecule.id, event.id
                    ));
                }
            }
        }
        Ok(())
    }
}

fn protocol_cartoon_template_schema_id() -> String {
    "gentle.protocol_cartoon_template.v1".to_string()
}

fn protocol_cartoon_template_bindings_schema_id() -> String {
    "gentle.protocol_cartoon_template_bindings.v1".to_string()
}

fn default_template_event_caption() -> String {
    "Event snapshot".to_string()
}

fn default_template_event_action() -> ProtocolCartoonAction {
    ProtocolCartoonAction::Context
}

fn default_template_molecule_topology() -> DnaTopologyCartoon {
    DnaTopologyCartoon::Linear
}

fn default_template_linear_end_style() -> DnaEndStyle {
    DnaEndStyle::Blunt
}

fn default_template_feature_length_bp() -> usize {
    120
}

fn default_template_palette() -> Vec<String> {
    vec![
        "#0f5d75".to_string(),
        "#1aa3a1".to_string(),
        "#e3b230".to_string(),
        "#4b7f52".to_string(),
    ]
}

fn default_template_summary() -> String {
    "Event-sequence protocol cartoon (feature fragments + end semantics)".to_string()
}

fn default_template_title(protocol_id: &str) -> String {
    format!("GENtle Protocol Cartoon: {}", protocol_id)
}

fn ensure_non_empty(value: &str, fallback: String) -> String {
    let trimmed = value.trim();
    if trimmed.is_empty() {
        fallback
    } else {
        trimmed.to_string()
    }
}

fn normalize_template_defaults(
    defaults: &ProtocolCartoonTemplateDefaults,
) -> ProtocolCartoonTemplateDefaults {
    let mut normalized = defaults.clone();
    normalized.event_caption =
        ensure_non_empty(&normalized.event_caption, default_template_event_caption());
    normalized.feature_length_bp = normalized.feature_length_bp.max(1);
    if normalized.palette.is_empty() {
        normalized.palette = default_template_palette();
    }
    normalized
}

fn update_if_some_string(target: &mut String, value: &Option<String>) {
    if let Some(value) = value {
        *target = value.clone();
    }
}

fn apply_bindings_defaults(
    template_defaults: &mut ProtocolCartoonTemplateDefaults,
    bindings_defaults: &ProtocolCartoonTemplateBindingsDefaults,
) {
    if let Some(value) = bindings_defaults.event_caption.as_ref() {
        template_defaults.event_caption = value.clone();
    }
    if let Some(value) = bindings_defaults.event_action.as_ref() {
        template_defaults.event_action = value.clone();
    }
    if let Some(value) = bindings_defaults.molecule_topology.as_ref() {
        template_defaults.molecule_topology = value.clone();
    }
    if let Some(value) = bindings_defaults.linear_end_style.as_ref() {
        template_defaults.linear_end_style = value.clone();
    }
    if let Some(value) = bindings_defaults.feature_length_bp {
        template_defaults.feature_length_bp = value.max(1);
    }
    if let Some(value) = bindings_defaults.palette.as_ref() {
        template_defaults.palette = value.clone();
    }
}

fn find_event_mut<'a>(
    template: &'a mut ProtocolCartoonTemplate,
    event_id: &str,
) -> Option<&'a mut ProtocolCartoonTemplateEvent> {
    template
        .events
        .iter_mut()
        .find(|event| event.id.trim() == event_id.trim())
}

fn find_molecule_mut<'a>(
    event: &'a mut ProtocolCartoonTemplateEvent,
    molecule_id: &str,
) -> Option<&'a mut ProtocolCartoonTemplateMolecule> {
    event
        .molecules
        .iter_mut()
        .find(|molecule| molecule.id.trim() == molecule_id.trim())
}

fn find_feature_mut<'a>(
    molecule: &'a mut ProtocolCartoonTemplateMolecule,
    feature_id: &str,
) -> Option<&'a mut ProtocolCartoonTemplateFeature> {
    molecule
        .features
        .iter_mut()
        .find(|feature| feature.id.trim() == feature_id.trim())
}

/// Apply deterministic bindings to a protocol cartoon template.
pub fn apply_protocol_cartoon_template_bindings(
    template: &ProtocolCartoonTemplate,
    bindings: &ProtocolCartoonTemplateBindings,
) -> Result<ProtocolCartoonTemplate, String> {
    let expected_schema = protocol_cartoon_template_bindings_schema_id();
    let schema = bindings.schema.trim();
    if !schema.is_empty() && schema != expected_schema {
        return Err(format!(
            "Unsupported protocol cartoon template bindings schema '{}' (expected '{}')",
            schema, expected_schema
        ));
    }

    if let Some(template_id) = bindings
        .template_id
        .as_deref()
        .map(str::trim)
        .filter(|value| !value.is_empty())
    {
        if template.id.trim() != template_id {
            return Err(format!(
                "Template id mismatch: bindings expect '{}' but template id is '{}'",
                template_id, template.id
            ));
        }
    }

    let mut bound = template.clone();
    let bound_template_id = bound.id.clone();
    apply_bindings_defaults(&mut bound.defaults, &bindings.defaults);

    for event_binding in &bindings.event_overrides {
        let event_id = event_binding.event_id.trim();
        if event_id.is_empty() {
            return Err("Event binding has empty event_id".to_string());
        }
        let event = find_event_mut(&mut bound, event_id).ok_or_else(|| {
            format!(
                "Event binding references unknown event_id '{}' in template '{}'",
                event_id, bound_template_id
            )
        })?;
        update_if_some_string(&mut event.title, &event_binding.title);
        update_if_some_string(&mut event.caption, &event_binding.caption);
        if let Some(action) = event_binding.action.as_ref() {
            event.action = Some(action.clone());
        }
    }

    for molecule_binding in &bindings.molecule_overrides {
        let event_id = molecule_binding.event_id.trim();
        let molecule_id = molecule_binding.molecule_id.trim();
        if event_id.is_empty() || molecule_id.is_empty() {
            return Err("Molecule binding requires non-empty event_id and molecule_id".to_string());
        }
        let event = find_event_mut(&mut bound, event_id).ok_or_else(|| {
            format!(
                "Molecule binding references unknown event_id '{}' in template '{}'",
                event_id, bound_template_id
            )
        })?;
        let molecule = find_molecule_mut(event, molecule_id).ok_or_else(|| {
            format!(
                "Molecule binding references unknown molecule_id '{}' in event '{}'",
                molecule_id, event_id
            )
        })?;
        update_if_some_string(&mut molecule.label, &molecule_binding.label);
        if let Some(topology) = molecule_binding.topology.as_ref() {
            molecule.topology = Some(topology.clone());
        }
        if let Some(left_end) = molecule_binding.left_end.as_ref() {
            molecule.left_end = Some(left_end.clone());
        }
        if let Some(right_end) = molecule_binding.right_end.as_ref() {
            molecule.right_end = Some(right_end.clone());
        }
    }

    for feature_binding in &bindings.feature_overrides {
        let event_id = feature_binding.event_id.trim();
        let molecule_id = feature_binding.molecule_id.trim();
        let feature_id = feature_binding.feature_id.trim();
        if event_id.is_empty() || molecule_id.is_empty() || feature_id.is_empty() {
            return Err(
                "Feature binding requires non-empty event_id, molecule_id, and feature_id"
                    .to_string(),
            );
        }
        let event = find_event_mut(&mut bound, event_id).ok_or_else(|| {
            format!(
                "Feature binding references unknown event_id '{}' in template '{}'",
                event_id, bound_template_id
            )
        })?;
        let molecule = find_molecule_mut(event, molecule_id).ok_or_else(|| {
            format!(
                "Feature binding references unknown molecule_id '{}' in event '{}'",
                molecule_id, event_id
            )
        })?;
        let feature = find_feature_mut(molecule, feature_id).ok_or_else(|| {
            format!(
                "Feature binding references unknown feature_id '{}' in molecule '{}'",
                feature_id, molecule_id
            )
        })?;
        update_if_some_string(&mut feature.label, &feature_binding.label);
        if let Some(length_bp) = feature_binding.length_bp {
            feature.length_bp = Some(length_bp.max(1));
        }
        if let Some(top_length_bp) = feature_binding.top_length_bp {
            feature.top_length_bp = Some(top_length_bp);
        }
        if let Some(bottom_length_bp) = feature_binding.bottom_length_bp {
            feature.bottom_length_bp = Some(bottom_length_bp);
        }
        if let Some(color_hex) = feature_binding.color_hex.as_ref() {
            feature.color_hex = Some(color_hex.clone());
        }
        if let Some(bottom_color_hex) = feature_binding.bottom_color_hex.as_ref() {
            feature.bottom_color_hex = Some(bottom_color_hex.clone());
        }
        if let Some(top_nick_after) = feature_binding.top_nick_after {
            feature.top_nick_after = Some(top_nick_after);
        }
        if let Some(bottom_nick_after) = feature_binding.bottom_nick_after {
            feature.bottom_nick_after = Some(bottom_nick_after);
        }
    }

    Ok(bound)
}

/// Resolve a sparse template into a concrete, render-ready protocol cartoon
/// spec with deterministic defaults.
pub fn resolve_protocol_cartoon_template(
    template: &ProtocolCartoonTemplate,
) -> Result<ProtocolCartoonSpec, String> {
    let expected_schema = protocol_cartoon_template_schema_id();
    let schema = template.schema.trim();
    if !schema.is_empty() && schema != expected_schema {
        return Err(format!(
            "Unsupported protocol cartoon template schema '{}' (expected '{}')",
            schema, expected_schema
        ));
    }

    let id = template.id.trim().to_string();
    if id.is_empty() {
        return Err("Protocol cartoon template id must be non-empty".to_string());
    }
    if template.events.is_empty() {
        return Err(format!(
            "Protocol cartoon template '{}' must contain at least one event",
            id
        ));
    }

    let defaults = normalize_template_defaults(&template.defaults);
    let title = ensure_non_empty(&template.title, default_template_title(&id));
    let summary = ensure_non_empty(&template.summary, default_template_summary());

    let mut events: Vec<ProtocolCartoonEvent> = vec![];
    for (event_idx, event) in template.events.iter().enumerate() {
        let event_id = ensure_non_empty(&event.id, format!("event_{}", event_idx + 1));
        let event_title = ensure_non_empty(&event.title, format!("Step {}", event_idx + 1));
        let event_caption = ensure_non_empty(&event.caption, defaults.event_caption.clone());
        let action = event
            .action
            .clone()
            .unwrap_or_else(|| defaults.event_action.clone());

        let template_molecules: Vec<ProtocolCartoonTemplateMolecule> = if event.molecules.is_empty()
        {
            vec![ProtocolCartoonTemplateMolecule::default()]
        } else {
            event.molecules.clone()
        };

        let mut molecules: Vec<DnaMoleculeCartoon> = vec![];
        for (molecule_idx, molecule) in template_molecules.iter().enumerate() {
            let topology = molecule
                .topology
                .clone()
                .unwrap_or_else(|| defaults.molecule_topology.clone());
            let molecule_id = ensure_non_empty(
                &molecule.id,
                format!("{}_molecule_{}", event_id, molecule_idx + 1),
            );
            let molecule_label =
                ensure_non_empty(&molecule.label, format!("Molecule {}", molecule_idx + 1));

            let template_features: Vec<ProtocolCartoonTemplateFeature> =
                if molecule.features.is_empty() {
                    vec![ProtocolCartoonTemplateFeature::default()]
                } else {
                    molecule.features.clone()
                };

            let mut features: Vec<DnaFeatureCartoon> = vec![];
            for (feature_idx, feature) in template_features.iter().enumerate() {
                let feature_id = ensure_non_empty(
                    &feature.id,
                    format!("{}_feature_{}", molecule_id, feature_idx + 1),
                );
                let feature_label =
                    ensure_non_empty(&feature.label, format!("Segment {}", feature_idx + 1));
                let length_bp = feature
                    .length_bp
                    .unwrap_or(defaults.feature_length_bp)
                    .max(1);
                let top_length_bp = feature.top_length_bp.unwrap_or(length_bp);
                let bottom_length_bp = feature.bottom_length_bp.unwrap_or(length_bp);
                let color_hex = feature.color_hex.clone().unwrap_or_else(|| {
                    defaults.palette[feature_idx % defaults.palette.len()].clone()
                });
                let bottom_color_hex = feature
                    .bottom_color_hex
                    .clone()
                    .unwrap_or_else(|| color_hex.clone());
                let top_nick_after = feature.top_nick_after.unwrap_or(false);
                let bottom_nick_after = feature.bottom_nick_after.unwrap_or(false);
                features.push(DnaFeatureCartoon {
                    id: feature_id,
                    label: feature_label,
                    length_bp,
                    top_length_bp,
                    bottom_length_bp,
                    color_hex,
                    bottom_color_hex,
                    top_nick_after,
                    bottom_nick_after,
                });
            }

            let (left_end, right_end) = match topology {
                DnaTopologyCartoon::Circular => (None, None),
                DnaTopologyCartoon::Linear => (
                    Some(
                        molecule
                            .left_end
                            .clone()
                            .unwrap_or_else(|| defaults.linear_end_style.clone()),
                    ),
                    Some(
                        molecule
                            .right_end
                            .clone()
                            .unwrap_or_else(|| defaults.linear_end_style.clone()),
                    ),
                ),
            };

            molecules.push(DnaMoleculeCartoon {
                id: molecule_id,
                label: molecule_label,
                topology,
                features,
                left_end,
                right_end,
            });
        }

        events.push(ProtocolCartoonEvent {
            id: event_id,
            title: event_title,
            caption: event_caption,
            action,
            molecules,
        });
    }

    let spec = ProtocolCartoonSpec {
        id,
        title,
        summary,
        events,
    };
    spec.validate()?;
    Ok(spec)
}

/// Resolve a sparse template into a concrete, render-ready protocol cartoon
/// spec after applying deterministic bindings.
pub fn resolve_protocol_cartoon_template_with_bindings(
    template: &ProtocolCartoonTemplate,
    bindings: &ProtocolCartoonTemplateBindings,
) -> Result<ProtocolCartoonSpec, String> {
    let bound_template = apply_protocol_cartoon_template_bindings(template, bindings)?;
    resolve_protocol_cartoon_template(&bound_template)
}

/// Render one protocol cartoon template as deterministic SVG.
pub fn render_protocol_cartoon_template_svg(template: &ProtocolCartoonTemplate) -> String {
    match resolve_protocol_cartoon_template(template) {
        Ok(spec) => render_protocol_cartoon_spec_svg(&spec),
        Err(message) => {
            let id = ensure_non_empty(&template.id, "unknown".to_string());
            render_invalid_protocol_svg(&id, &message)
        }
    }
}

/// Render one protocol cartoon template as deterministic SVG after applying
/// bindings.
pub fn render_protocol_cartoon_template_with_bindings_svg(
    template: &ProtocolCartoonTemplate,
    bindings: &ProtocolCartoonTemplateBindings,
) -> String {
    match resolve_protocol_cartoon_template_with_bindings(template, bindings) {
        Ok(spec) => render_protocol_cartoon_spec_svg(&spec),
        Err(message) => {
            let id = ensure_non_empty(&template.id, "unknown".to_string());
            render_invalid_protocol_svg(&id, &message)
        }
    }
}

pub fn protocol_cartoon_catalog_rows() -> Vec<ProtocolCartoonCatalogRow> {
    ProtocolCartoonKind::catalog()
        .into_iter()
        .map(|kind| ProtocolCartoonCatalogRow {
            id: kind.id().to_string(),
            title: kind.title().to_string(),
            summary: kind.summary().to_string(),
        })
        .collect()
}

/// Build the reusable protocol-cartoon template for one built-in protocol.
pub fn protocol_cartoon_template_for_kind(kind: &ProtocolCartoonKind) -> ProtocolCartoonTemplate {
    match kind {
        ProtocolCartoonKind::GibsonTwoFragment => gibson_two_fragment_template(),
        ProtocolCartoonKind::GibsonSingleInsertDualJunction => {
            gibson_single_insert_dual_junction_template()
        }
        ProtocolCartoonKind::PcrAssayPair => pcr_assay_pair_template(),
        ProtocolCartoonKind::PcrAssayPairNoProduct => pcr_assay_pair_no_product_template(),
        ProtocolCartoonKind::PcrAssayQpcr => pcr_assay_qpcr_template(),
    }
}

/// Build the reusable event/molecule abstraction for one built-in protocol.
pub fn protocol_cartoon_spec_for_kind(kind: &ProtocolCartoonKind) -> ProtocolCartoonSpec {
    let template = protocol_cartoon_template_for_kind(kind);
    resolve_protocol_cartoon_template(&template).unwrap_or_else(|message| ProtocolCartoonSpec {
        id: template.id.clone(),
        title: "Invalid protocol cartoon template".to_string(),
        summary: message,
        events: vec![ProtocolCartoonEvent {
            id: "invalid_template".to_string(),
            title: "Invalid template".to_string(),
            caption: "Template resolution failed".to_string(),
            action: ProtocolCartoonAction::Custom {
                label: "Invalid".to_string(),
            },
            molecules: vec![DnaMoleculeCartoon {
                id: "invalid_template_molecule".to_string(),
                label: "Invalid template".to_string(),
                topology: DnaTopologyCartoon::Linear,
                features: vec![DnaFeatureCartoon {
                    id: "invalid_feature".to_string(),
                    label: "Invalid".to_string(),
                    length_bp: 1,
                    top_length_bp: 1,
                    bottom_length_bp: 1,
                    color_hex: "#8ea7b1".to_string(),
                    bottom_color_hex: "#8ea7b1".to_string(),
                    top_nick_after: false,
                    bottom_nick_after: false,
                }],
                left_end: Some(DnaEndStyle::Blunt),
                right_end: Some(DnaEndStyle::Blunt),
            }],
        }],
    })
}

/// Render one protocol cartoon as deterministic SVG.
pub fn render_protocol_cartoon_svg(kind: &ProtocolCartoonKind) -> String {
    let spec = protocol_cartoon_spec_for_kind(kind);
    render_protocol_cartoon_spec_svg(&spec)
}

/// Render a protocol cartoon from a full event/molecule specification.
pub fn render_protocol_cartoon_spec_svg(spec: &ProtocolCartoonSpec) -> String {
    if let Err(message) = spec.validate() {
        return render_invalid_protocol_svg(&spec.id, &message);
    }

    let panel_w = 308.0_f32;
    let max_molecules = spec
        .events
        .iter()
        .map(|event| event.molecules.len())
        .max()
        .unwrap_or(1);
    let panel_h = if max_molecules <= 1 { 286.0 } else { 324.0 };
    let panel_gap = 24.0_f32;
    let margin_x = 36.0_f32;
    let margin_y = 28.0_f32;
    let header_h = 86.0_f32;
    let footer_h = 20.0_f32;

    let event_count = spec.events.len() as f32;
    let body_w = event_count * panel_w + (event_count - 1.0).max(0.0) * panel_gap;
    let canvas_w = margin_x * 2.0 + body_w;
    let canvas_h = margin_y * 2.0 + header_h + panel_h + footer_h;

    let mut svg = String::new();
    svg.push_str(&format!(
        "<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"{:.0}\" height=\"{:.0}\" viewBox=\"0 0 {:.0} {:.0}\" role=\"img\" aria-labelledby=\"title desc\" data-protocol-id=\"{}\">",
        canvas_w,
        canvas_h,
        canvas_w,
        canvas_h,
        escape_xml(&spec.id)
    ));
    svg.push_str(&format!(
        "<title id=\"title\">{}</title>",
        escape_xml(&spec.title)
    ));
    svg.push_str(&format!(
        "<desc id=\"desc\">{}</desc>",
        escape_xml(&spec.summary)
    ));
    svg.push_str("<defs>");
    svg.push_str(
        "<linearGradient id=\"pc_bg\" x1=\"0\" y1=\"0\" x2=\"1\" y2=\"1\"><stop offset=\"0%\" stop-color=\"#f5fbff\"/><stop offset=\"100%\" stop-color=\"#edf8f5\"/></linearGradient>",
    );
    svg.push_str(
        "<linearGradient id=\"pc_panel\" x1=\"0\" y1=\"0\" x2=\"0\" y2=\"1\"><stop offset=\"0%\" stop-color=\"#ffffff\"/><stop offset=\"100%\" stop-color=\"#f6fafb\"/></linearGradient>",
    );
    svg.push_str(
        "<marker id=\"pc_arrow\" markerWidth=\"12\" markerHeight=\"12\" refX=\"9\" refY=\"6\" orient=\"auto\"><path d=\"M1,1 L11,6 L1,11 Z\" fill=\"#4e6d78\"/></marker>",
    );
    svg.push_str("<style>");
    svg.push_str(
        ".pc_title { font: 700 36px 'Avenir Next', 'Trebuchet MS', 'Segoe UI', sans-serif; fill: #0f2a36; }",
    );
    svg.push_str(
        ".pc_sub { font: 500 18px 'Avenir Next', 'Trebuchet MS', 'Segoe UI', sans-serif; fill: #36535d; }",
    );
    svg.push_str(
        ".pc_panel_title { font: 700 21px 'Avenir Next', 'Trebuchet MS', 'Segoe UI', sans-serif; fill: #123746; }",
    );
    svg.push_str(
        ".pc_panel_sub { font: 600 13px 'Avenir Next', 'Trebuchet MS', 'Segoe UI', sans-serif; fill: #3f5f68; }",
    );
    svg.push_str(
        ".pc_caption { font: 500 13px 'Avenir Next', 'Trebuchet MS', 'Segoe UI', sans-serif; fill: #466570; }",
    );
    svg.push_str(
        ".pc_mol_label { font: 600 13px 'Avenir Next', 'Trebuchet MS', 'Segoe UI', sans-serif; fill: #2d4f5b; }",
    );
    svg.push_str(
        ".pc_end_label { font: 600 10px 'Avenir Next', 'Trebuchet MS', 'Segoe UI', sans-serif; fill: #375965; }",
    );
    svg.push_str(
        ".pc_feature_marker_text { font: 700 10px 'Avenir Next', 'Trebuchet MS', 'Segoe UI', sans-serif; fill: #ffffff; }",
    );
    svg.push_str(
        ".pc_primer_label { font: 700 10px 'Avenir Next', 'Trebuchet MS', 'Segoe UI', sans-serif; fill: #ffffff; }",
    );
    svg.push_str(
        ".pc_primer_end_label { font: 600 9px 'Avenir Next', 'Trebuchet MS', 'Segoe UI', sans-serif; fill: #375965; }",
    );
    svg.push_str(
        ".pc_step_num { font: 700 13px 'Avenir Next', 'Trebuchet MS', 'Segoe UI', sans-serif; fill: #ffffff; }",
    );
    svg.push_str(".pc_box { fill: url(#pc_panel); stroke: #c8dde3; stroke-width: 2; }");
    svg.push_str(
        ".pc_link { stroke: #4e6d78; stroke-width: 3.2; marker-end: url(#pc_arrow); fill: none; }",
    );
    svg.push_str("</style>");
    svg.push_str("</defs>");

    svg.push_str(&format!(
        "<rect x=\"0\" y=\"0\" width=\"{:.0}\" height=\"{:.0}\" fill=\"url(#pc_bg)\"/>",
        canvas_w, canvas_h
    ));
    svg.push_str(&format!(
        "<rect x=\"14\" y=\"14\" width=\"{:.0}\" height=\"{:.0}\" rx=\"22\" fill=\"none\" stroke=\"#d4e4e8\" stroke-width=\"2\"/>",
        canvas_w - 28.0,
        canvas_h - 28.0
    ));

    svg.push_str(&format!(
        "<text x=\"{:.1}\" y=\"{:.1}\" class=\"pc_title\">{}</text>",
        margin_x,
        margin_y + 38.0,
        escape_xml(&spec.title)
    ));
    svg.push_str(&format!(
        "<text x=\"{:.1}\" y=\"{:.1}\" class=\"pc_sub\">{}</text>",
        margin_x,
        margin_y + 66.0,
        escape_xml(&spec.summary)
    ));

    let panel_top = margin_y + header_h;
    for (idx, event) in spec.events.iter().enumerate() {
        let panel_x = margin_x + idx as f32 * (panel_w + panel_gap);
        svg.push_str(&format!(
            "<rect class=\"pc_box\" x=\"{:.1}\" y=\"{:.1}\" width=\"{:.1}\" height=\"{:.1}\" rx=\"18\" data-event-id=\"{}\" data-action=\"{}\"/>",
            panel_x,
            panel_top,
            panel_w,
            panel_h,
            escape_xml(&event.id),
            escape_xml(event.action.label())
        ));

        let badge_cx = panel_x + 24.0;
        let badge_cy = panel_top + 24.0;
        svg.push_str(&format!(
            "<circle cx=\"{:.1}\" cy=\"{:.1}\" r=\"12\" fill=\"#2b8ea8\"/>",
            badge_cx, badge_cy
        ));
        svg.push_str(&format!(
            "<text x=\"{:.1}\" y=\"{:.1}\" class=\"pc_step_num\" text-anchor=\"middle\" dominant-baseline=\"middle\">{}</text>",
            badge_cx,
            badge_cy + 0.5,
            idx + 1
        ));

        svg.push_str(&format!(
            "<text x=\"{:.1}\" y=\"{:.1}\" class=\"pc_panel_title\">{}</text>",
            panel_x + 46.0,
            panel_top + 30.0,
            escape_xml(&event.title)
        ));
        svg.push_str(&format!(
            "<text x=\"{:.1}\" y=\"{:.1}\" class=\"pc_panel_sub\">{}</text>",
            panel_x + 18.0,
            panel_top + 52.0,
            escape_xml(event.action.label())
        ));

        append_wrapped_text(
            &mut svg,
            panel_x + 18.0,
            panel_top + 72.0,
            46,
            16.0,
            &event.caption,
            "pc_caption",
        );

        let molecule_count = event.molecules.len();
        let molecule_y0 = panel_top + 126.0;
        let molecule_y1 = panel_top + panel_h - 40.0;
        let molecule_step = if molecule_count <= 1 {
            0.0
        } else {
            ((molecule_y1 - molecule_y0) / (molecule_count - 1) as f32).clamp(52.0, 72.0)
        };
        for (molecule_idx, molecule) in event.molecules.iter().enumerate() {
            let molecule_y = molecule_y0 + molecule_idx as f32 * molecule_step;
            render_molecule(
                &mut svg,
                molecule,
                panel_x + 16.0,
                molecule_y,
                panel_w - 32.0,
            );
        }

        if idx + 1 < spec.events.len() {
            let x0 = panel_x + panel_w + 4.0;
            let x1 = panel_x + panel_w + panel_gap - 6.0;
            let y = panel_top + panel_h * 0.54;
            svg.push_str(&format!(
                "<path d=\"M {:.1} {:.1} L {:.1} {:.1}\" class=\"pc_link\"/>",
                x0, y, x1, y
            ));
        }
    }

    svg.push_str(&format!(
        "<text x=\"{:.1}\" y=\"{:.1}\" class=\"pc_caption\">Model: event-sequence abstraction (features + overhangs + topology) | deterministic SVG export.</text>",
        margin_x,
        canvas_h - 16.0
    ));
    svg.push_str("</svg>");
    svg
}

fn render_molecule(svg: &mut String, molecule: &DnaMoleculeCartoon, x: f32, y: f32, width: f32) {
    svg.push_str(&format!(
        "<g data-molecule-id=\"{}\" data-topology=\"{}\">",
        escape_xml(&molecule.id),
        match molecule.topology {
            DnaTopologyCartoon::Linear => "linear",
            DnaTopologyCartoon::Circular => "circular",
        }
    ));

    svg.push_str(&format!(
        "<text x=\"{:.1}\" y=\"{:.1}\" class=\"pc_mol_label\">{}</text>",
        x,
        y,
        escape_xml(&molecule.label)
    ));

    match molecule.topology {
        DnaTopologyCartoon::Linear => {
            render_linear_molecule(svg, molecule, x, y + 18.0, width);
        }
        DnaTopologyCartoon::Circular => {
            render_circular_molecule(svg, molecule, x, y + 18.0, width);
        }
    }

    svg.push_str("</g>");
}

fn render_linear_molecule(
    svg: &mut String,
    molecule: &DnaMoleculeCartoon,
    x: f32,
    y: f32,
    width: f32,
) {
    let canonical_total_bp: usize = molecule
        .features
        .iter()
        .map(|f| f.length_bp)
        .sum::<usize>()
        .max(1);
    let body_w = (width * 0.62).clamp(132.0, 188.0);
    let body_x = x + (width - body_w) * 0.5;

    let mut left_top = body_x;
    let mut left_bottom = body_x;
    let mut right_top = body_x + body_w;
    let mut right_bottom = body_x + body_w;

    if let Some(end) = molecule.left_end.as_ref() {
        apply_left_end_style(&mut left_top, &mut left_bottom, end);
        render_end_annotation(svg, left_top, left_bottom, y + 9.5, true, end);
    }
    if let Some(end) = molecule.right_end.as_ref() {
        apply_right_end_style(&mut right_top, &mut right_bottom, end);
        render_end_annotation(svg, right_top, right_bottom, y + 9.5, false, end);
    }

    let shared_left = left_top.min(left_bottom);
    let shared_right = right_top.max(right_bottom);
    let shared_w = (shared_right - shared_left).max(2.0);
    let mut slot_cursor = shared_left;
    let mut top_spans: Vec<(f32, f32, String)> = vec![];
    let mut bottom_spans: Vec<(f32, f32, String)> = vec![];
    let mut top_nicks: Vec<f32> = vec![];
    let mut bottom_nicks: Vec<f32> = vec![];
    let mut primer_glyphs: Vec<(f32, f32, PrimerGlyphKind, String)> = vec![];
    for (feature_idx, feature) in molecule.features.iter().enumerate() {
        let slot_w = if feature_idx + 1 == molecule.features.len() {
            (shared_right - slot_cursor).max(0.5)
        } else {
            (shared_w * (feature.length_bp as f32 / canonical_total_bp as f32)).max(0.5)
        };
        let top_seg_w = if feature.top_length_bp == 0 {
            0.0
        } else {
            (slot_w * (feature.top_length_bp as f32 / feature.length_bp as f32)).max(0.5)
        };
        let bottom_seg_w = if feature.bottom_length_bp == 0 {
            0.0
        } else {
            (slot_w * (feature.bottom_length_bp as f32 / feature.length_bp as f32)).max(0.5)
        };
        let top_color = normalize_hex_color(&feature.color_hex);
        let bottom_color = normalize_hex_color(&feature.bottom_color_hex);
        if top_seg_w > 0.0 {
            push_linear_span(&mut top_spans, slot_cursor, top_seg_w, top_color);
        }
        if bottom_seg_w > 0.0 {
            push_linear_span(&mut bottom_spans, slot_cursor, bottom_seg_w, bottom_color);
        }
        if let Some(primer_kind) = primer_glyph_kind(feature) {
            primer_glyphs.push((
                slot_cursor + (slot_w * 0.5),
                slot_w,
                primer_kind,
                normalize_hex_color(&feature.color_hex),
            ));
        }
        let boundary_x = slot_cursor + slot_w;
        if feature.top_nick_after
            && feature_idx + 1 < molecule.features.len()
            && feature.top_length_bp > 0
            && molecule.features[feature_idx + 1].top_length_bp > 0
        {
            top_nicks.push(boundary_x);
        }
        if feature.bottom_nick_after
            && feature_idx + 1 < molecule.features.len()
            && feature.bottom_length_bp > 0
            && molecule.features[feature_idx + 1].bottom_length_bp > 0
        {
            bottom_nicks.push(boundary_x);
        }
        slot_cursor += slot_w;
    }

    for (span_x, span_w, color) in top_spans {
        svg.push_str(&format!(
            "<rect x=\"{:.1}\" y=\"{:.1}\" width=\"{:.1}\" height=\"8\" fill=\"{}\"/>",
            span_x, y, span_w, color
        ));
    }
    for (span_x, span_w, color) in bottom_spans {
        svg.push_str(&format!(
            "<rect x=\"{:.1}\" y=\"{:.1}\" width=\"{:.1}\" height=\"8\" fill=\"{}\"/>",
            span_x,
            y + 11.0,
            span_w,
            color
        ));
    }
    for nick_x in top_nicks {
        render_strand_nick(svg, nick_x, y);
    }
    for nick_x in bottom_nicks {
        render_strand_nick(svg, nick_x, y + 11.0);
    }
    for (center_x, primer_w, primer_kind, primer_fill) in primer_glyphs {
        match primer_kind {
            PrimerGlyphKind::Forward => {
                render_primer_glyph(svg, center_x, y, primer_w, "F", &primer_fill, "top", true);
            }
            PrimerGlyphKind::Reverse => {
                render_primer_glyph(
                    svg,
                    center_x,
                    y + 11.0,
                    primer_w,
                    "R",
                    &primer_fill,
                    "bottom",
                    true,
                );
            }
            PrimerGlyphKind::Probe => {
                render_primer_glyph(svg, center_x, y, primer_w, "P", &primer_fill, "top", false);
            }
        }
    }

    if let Some(DnaEndStyle::Continuation) = molecule.left_end.as_ref() {
        render_continuation_marker(svg, left_top.min(left_bottom), y, true);
    }
    if let Some(DnaEndStyle::Continuation) = molecule.right_end.as_ref() {
        render_continuation_marker(svg, right_top.max(right_bottom), y, false);
    }
}

fn push_linear_span(spans: &mut Vec<(f32, f32, String)>, x: f32, width: f32, color: String) {
    if let Some((last_x, last_w, last_color)) = spans.last_mut() {
        let last_end = *last_x + *last_w;
        if *last_color == color && (last_end - x).abs() <= 0.2 {
            *last_w = (x + width) - *last_x;
            return;
        }
    }
    spans.push((x, width, color));
}

fn render_strand_nick(svg: &mut String, x: f32, y: f32) {
    svg.push_str(&format!(
        "<rect x=\"{:.1}\" y=\"{:.1}\" width=\"2.4\" height=\"8\" rx=\"1.0\" fill=\"url(#pc_bg)\"/>",
        x - 1.2,
        y
    ));
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum PrimerGlyphKind {
    Forward,
    Reverse,
    Probe,
}

fn primer_glyph_kind(feature: &DnaFeatureCartoon) -> Option<PrimerGlyphKind> {
    let feature_id = feature.id.to_ascii_lowercase();
    if feature_id.contains("forward_primer") {
        Some(PrimerGlyphKind::Forward)
    } else if feature_id.contains("reverse_primer") {
        Some(PrimerGlyphKind::Reverse)
    } else if feature_id.contains("probe_site") {
        Some(PrimerGlyphKind::Probe)
    } else {
        None
    }
}

fn render_primer_glyph(
    svg: &mut String,
    center_x: f32,
    anchor_y: f32,
    nominal_w: f32,
    label: &str,
    fill: &str,
    anchor: &str,
    with_end_labels: bool,
) {
    let glyph_w = nominal_w.max(20.0).min(40.0);
    let glyph_h = 5.0;
    let anchor_lower = anchor.to_ascii_lowercase();
    let on_bottom = anchor_lower == "bottom";
    let connector_y = if on_bottom {
        anchor_y + 5.0
    } else {
        anchor_y - 2.0
    };
    let bar_y = if on_bottom {
        anchor_y + 10.0
    } else {
        anchor_y - 8.0
    };
    let end_label_y = if on_bottom {
        anchor_y + 21.5
    } else {
        anchor_y - 11.5
    };
    let label_y = bar_y + (glyph_h * 0.5) + 0.4;
    let left_x = center_x - (glyph_w * 0.5);
    let (left_end_label, right_end_label) = if label == "R" {
        ("3'", "5'")
    } else {
        ("5'", "3'")
    };

    svg.push_str(&format!(
        "<g data-primer-glyph=\"{}\" data-primer-anchor=\"{}\" data-primer-rotation=\"0\">",
        escape_xml(label),
        escape_xml(anchor)
    ));
    svg.push_str(&format!(
        "<line x1=\"{:.1}\" y1=\"{:.1}\" x2=\"{:.1}\" y2=\"{:.1}\" stroke=\"{}\" stroke-width=\"1.4\"/>",
        center_x,
        anchor_y,
        center_x,
        connector_y,
        fill
    ));
    svg.push_str(&format!(
        "<rect x=\"{:.1}\" y=\"{:.1}\" width=\"{:.1}\" height=\"{:.1}\" rx=\"2.2\" fill=\"{}\" stroke=\"{}\" stroke-width=\"0.8\"/>",
        left_x,
        bar_y,
        glyph_w,
        glyph_h,
        fill,
        fill
    ));
    svg.push_str(&format!(
        "<text x=\"{:.1}\" y=\"{:.1}\" class=\"pc_primer_label\" text-anchor=\"middle\" dominant-baseline=\"middle\">{}</text>",
        center_x,
        label_y,
        escape_xml(label)
    ));
    if with_end_labels {
        svg.push_str(&format!(
            "<text x=\"{:.1}\" y=\"{:.1}\" class=\"pc_primer_end_label\" text-anchor=\"end\" dominant-baseline=\"middle\">{}</text>",
            left_x - 3.5,
            end_label_y,
            escape_xml(left_end_label)
        ));
        svg.push_str(&format!(
            "<text x=\"{:.1}\" y=\"{:.1}\" class=\"pc_primer_end_label\" text-anchor=\"start\" dominant-baseline=\"middle\">{}</text>",
            left_x + glyph_w + 3.5,
            end_label_y,
            escape_xml(right_end_label)
        ));
    }
    svg.push_str("</g>");
}

fn render_circular_molecule(
    svg: &mut String,
    molecule: &DnaMoleculeCartoon,
    x: f32,
    y: f32,
    width: f32,
) {
    let total_top_bp: usize = molecule
        .features
        .iter()
        .map(|f| f.top_length_bp)
        .sum::<usize>()
        .max(1);
    let total_bottom_bp: usize = molecule
        .features
        .iter()
        .map(|f| f.bottom_length_bp)
        .sum::<usize>()
        .max(1);
    let cx = x + width * 0.5;
    let cy = y + 18.0;
    let top_r = (width * 0.22).clamp(14.0, 28.0);
    let bottom_r = top_r + 7.0;

    svg.push_str(&format!(
        "<circle cx=\"{:.1}\" cy=\"{:.1}\" r=\"{:.1}\" fill=\"none\" stroke=\"#c8d9de\" stroke-width=\"4\"/>",
        cx, cy, top_r
    ));
    svg.push_str(&format!(
        "<circle cx=\"{:.1}\" cy=\"{:.1}\" r=\"{:.1}\" fill=\"none\" stroke=\"#c8d9de\" stroke-width=\"4\"/>",
        cx, cy, bottom_r
    ));

    let mut start_top = -std::f32::consts::FRAC_PI_2;
    let mut start_bottom = -std::f32::consts::FRAC_PI_2;
    for feature in &molecule.features {
        let top_color = normalize_hex_color(&feature.color_hex);
        let bottom_color = normalize_hex_color(&feature.bottom_color_hex);
        if feature.top_length_bp > 0 {
            let top_span =
                std::f32::consts::TAU * (feature.top_length_bp as f32 / total_top_bp as f32);
            let end_top = start_top + top_span;
            let top_path = arc_path(cx, cy, top_r, start_top, end_top);
            svg.push_str(&format!(
                "<path d=\"{}\" fill=\"none\" stroke=\"{}\" stroke-width=\"4\" stroke-linecap=\"butt\"/>",
                top_path, top_color
            ));
            start_top = end_top;
        }
        if feature.bottom_length_bp > 0 {
            let bottom_span =
                std::f32::consts::TAU * (feature.bottom_length_bp as f32 / total_bottom_bp as f32);
            let end_bottom = start_bottom + bottom_span;
            let bottom_path = arc_path(cx, cy, bottom_r, start_bottom, end_bottom);
            svg.push_str(&format!(
                "<path d=\"{}\" fill=\"none\" stroke=\"{}\" stroke-width=\"4\" stroke-linecap=\"butt\"/>",
                bottom_path, bottom_color
            ));
            start_bottom = end_bottom;
        }
    }
}

fn arc_path(cx: f32, cy: f32, r: f32, start: f32, end: f32) -> String {
    let sx = cx + r * start.cos();
    let sy = cy + r * start.sin();
    let ex = cx + r * end.cos();
    let ey = cy + r * end.sin();
    let large_arc = if end - start > std::f32::consts::PI {
        1
    } else {
        0
    };
    format!(
        "M {:.2} {:.2} A {:.2} {:.2} 0 {} 1 {:.2} {:.2}",
        sx, sy, r, r, large_arc, ex, ey
    )
}

fn overhang_len_px(nt: usize) -> f32 {
    if nt == 0 {
        return 0.0;
    }
    (nt as f32 * 1.7).clamp(7.0, 26.0)
}

fn apply_left_end_style(left_top: &mut f32, left_bottom: &mut f32, end: &DnaEndStyle) {
    match end {
        DnaEndStyle::NotShown => {}
        DnaEndStyle::Continuation => {}
        DnaEndStyle::Blunt => {}
        DnaEndStyle::Sticky { polarity, nt } => {
            let delta = overhang_len_px(*nt);
            match polarity {
                OverhangPolarity::FivePrime => *left_top -= delta,
                OverhangPolarity::ThreePrime => *left_bottom -= delta,
            }
        }
    }
}

fn apply_right_end_style(right_top: &mut f32, right_bottom: &mut f32, end: &DnaEndStyle) {
    match end {
        DnaEndStyle::NotShown => {}
        DnaEndStyle::Continuation => {}
        DnaEndStyle::Blunt => {}
        DnaEndStyle::Sticky { polarity, nt } => {
            let delta = overhang_len_px(*nt);
            match polarity {
                OverhangPolarity::FivePrime => *right_bottom += delta,
                OverhangPolarity::ThreePrime => *right_top += delta,
            }
        }
    }
}

fn render_end_annotation(
    svg: &mut String,
    top_edge: f32,
    bottom_edge: f32,
    y: f32,
    is_left: bool,
    end: &DnaEndStyle,
) {
    let text = match end {
        DnaEndStyle::NotShown => return,
        DnaEndStyle::Continuation => return,
        DnaEndStyle::Blunt => "blunt".to_string(),
        DnaEndStyle::Sticky { polarity, nt } => format!("{} {}nt", polarity.label(), nt),
    };
    let (x, anchor) = match end {
        DnaEndStyle::Blunt => {
            if is_left {
                (top_edge - 6.0, "end")
            } else {
                (top_edge + 6.0, "start")
            }
        }
        DnaEndStyle::Sticky { .. } => (((top_edge + bottom_edge) * 0.5), "middle"),
        DnaEndStyle::NotShown | DnaEndStyle::Continuation => unreachable!(),
    };
    let y = match end {
        DnaEndStyle::Blunt => y,
        DnaEndStyle::Sticky { .. } => {
            if is_left {
                y - 1.0
            } else {
                y + 4.0
            }
        }
        DnaEndStyle::NotShown | DnaEndStyle::Continuation => unreachable!(),
    };
    svg.push_str(&format!(
        "<text x=\"{:.1}\" y=\"{:.1}\" class=\"pc_end_label\" text-anchor=\"{}\" dominant-baseline=\"middle\">{}</text>",
        x,
        y,
        anchor,
        escape_xml(&text)
    ));
}

fn render_continuation_marker(svg: &mut String, edge_x: f32, y: f32, is_left: bool) {
    let zig = 7.0_f32;
    let gap = 6.0_f32;
    let dot_base_x = if is_left {
        edge_x - 20.0
    } else {
        edge_x + 20.0
    };
    let dot_step = if is_left { 6.0 } else { -6.0 };
    if is_left {
        svg.push_str(&format!(
            "<path d=\"M {:.1} {:.1} L {:.1} {:.1} L {:.1} {:.1} L {:.1} {:.1} L {:.1} {:.1} Z\" fill=\"url(#pc_bg)\"/>",
            edge_x + 0.4,
            y - 0.8,
            edge_x - zig,
            y + 1.4,
            edge_x + 0.4,
            y + 3.8,
            edge_x - zig,
            y + 6.2,
            edge_x + 0.4,
            y + 8.8
        ));
        svg.push_str(&format!(
            "<path d=\"M {:.1} {:.1} L {:.1} {:.1} L {:.1} {:.1} L {:.1} {:.1} L {:.1} {:.1} Z\" fill=\"url(#pc_bg)\"/>",
            edge_x + 0.4,
            y + 10.2,
            edge_x - zig,
            y + 12.4,
            edge_x + 0.4,
            y + 14.8,
            edge_x - zig,
            y + 17.2,
            edge_x + 0.4,
            y + 19.8
        ));
    } else {
        svg.push_str(&format!(
            "<path d=\"M {:.1} {:.1} L {:.1} {:.1} L {:.1} {:.1} L {:.1} {:.1} L {:.1} {:.1} Z\" fill=\"url(#pc_bg)\"/>",
            edge_x - 0.4,
            y - 0.8,
            edge_x + zig,
            y + 1.4,
            edge_x - 0.4,
            y + 3.8,
            edge_x + zig,
            y + 6.2,
            edge_x - 0.4,
            y + 8.8
        ));
        svg.push_str(&format!(
            "<path d=\"M {:.1} {:.1} L {:.1} {:.1} L {:.1} {:.1} L {:.1} {:.1} L {:.1} {:.1} Z\" fill=\"url(#pc_bg)\"/>",
            edge_x - 0.4,
            y + 10.2,
            edge_x + zig,
            y + 12.4,
            edge_x - 0.4,
            y + 14.8,
            edge_x + zig,
            y + 17.2,
            edge_x - 0.4,
            y + 19.8
        ));
    }
    for idx in 0..3 {
        svg.push_str(&format!(
            "<circle cx=\"{:.1}\" cy=\"{:.1}\" r=\"1.5\" fill=\"#88a3ad\"/>",
            dot_base_x + idx as f32 * dot_step,
            y + gap + idx as f32 * 2.6
        ));
    }
}

fn append_wrapped_text(
    svg: &mut String,
    x: f32,
    y: f32,
    max_chars: usize,
    line_h: f32,
    text: &str,
    class: &str,
) {
    let mut lines: Vec<String> = vec![];
    for paragraph in text.split('\n') {
        let words = paragraph.split_whitespace().collect::<Vec<_>>();
        if words.is_empty() {
            continue;
        }
        let mut current = String::new();
        for word in words {
            let candidate = if current.is_empty() {
                word.to_string()
            } else {
                format!("{} {}", current, word)
            };
            if candidate.chars().count() > max_chars && !current.is_empty() {
                lines.push(current);
                current = word.to_string();
            } else {
                current = candidate;
            }
        }
        if !current.is_empty() {
            lines.push(current);
        }
    }

    for (idx, line) in lines.into_iter().take(3).enumerate() {
        svg.push_str(&format!(
            "<text x=\"{:.1}\" y=\"{:.1}\" class=\"{}\">{}</text>",
            x,
            y + idx as f32 * line_h,
            class,
            escape_xml(&line)
        ));
    }
}

fn normalize_hex_color(raw: &str) -> String {
    let trimmed = raw.trim();
    if trimmed.len() == 7
        && trimmed.starts_with('#')
        && trimmed[1..].chars().all(|c| c.is_ascii_hexdigit())
    {
        return trimmed.to_string();
    }
    "#8ea7b1".to_string()
}

fn escape_xml(raw: &str) -> String {
    raw.replace('&', "&amp;")
        .replace('<', "&lt;")
        .replace('>', "&gt;")
        .replace('"', "&quot;")
        .replace('\'', "&apos;")
}

#[derive(Clone, Copy)]
enum CartoonBlockOccupancy {
    Duplex,
    TopOnly,
    BottomOnly,
}

fn strand_lengths(length_bp: usize, occupancy: CartoonBlockOccupancy) -> (usize, usize) {
    match occupancy {
        CartoonBlockOccupancy::Duplex => (length_bp, length_bp),
        CartoonBlockOccupancy::TopOnly => (length_bp, 0),
        CartoonBlockOccupancy::BottomOnly => (0, length_bp),
    }
}

fn cartoon_feature_block(
    id: &str,
    label: &str,
    length_bp: usize,
    occupancy: CartoonBlockOccupancy,
    top_color_hex: &str,
    bottom_color_hex: &str,
    top_nick_after: bool,
    bottom_nick_after: bool,
) -> DnaFeatureCartoon {
    let (top_length_bp, bottom_length_bp) = strand_lengths(length_bp, occupancy);
    DnaFeatureCartoon {
        id: id.to_string(),
        label: label.to_string(),
        length_bp,
        top_length_bp,
        bottom_length_bp,
        color_hex: top_color_hex.to_string(),
        bottom_color_hex: bottom_color_hex.to_string(),
        top_nick_after,
        bottom_nick_after,
    }
}

fn duplex_feature_block(
    id: &str,
    label: &str,
    length_bp: usize,
    color_hex: &str,
) -> DnaFeatureCartoon {
    cartoon_feature_block(
        id,
        label,
        length_bp,
        CartoonBlockOccupancy::Duplex,
        color_hex,
        color_hex,
        false,
        false,
    )
}

fn duplex_feature_block_with_nicks(
    id: &str,
    label: &str,
    length_bp: usize,
    color_hex: &str,
    top_nick_after: bool,
    bottom_nick_after: bool,
) -> DnaFeatureCartoon {
    cartoon_feature_block(
        id,
        label,
        length_bp,
        CartoonBlockOccupancy::Duplex,
        color_hex,
        color_hex,
        top_nick_after,
        bottom_nick_after,
    )
}

fn top_only_feature_block(
    id: &str,
    label: &str,
    length_bp: usize,
    color_hex: &str,
) -> DnaFeatureCartoon {
    cartoon_feature_block(
        id,
        label,
        length_bp,
        CartoonBlockOccupancy::TopOnly,
        color_hex,
        color_hex,
        false,
        false,
    )
}

fn bottom_only_feature_block(
    id: &str,
    label: &str,
    length_bp: usize,
    color_hex: &str,
) -> DnaFeatureCartoon {
    cartoon_feature_block(
        id,
        label,
        length_bp,
        CartoonBlockOccupancy::BottomOnly,
        color_hex,
        color_hex,
        false,
        false,
    )
}

#[cfg(test)]
fn hybrid_feature_block(
    id: &str,
    label: &str,
    length_bp: usize,
    top_color_hex: &str,
    bottom_color_hex: &str,
) -> DnaFeatureCartoon {
    cartoon_feature_block(
        id,
        label,
        length_bp,
        CartoonBlockOccupancy::Duplex,
        top_color_hex,
        bottom_color_hex,
        false,
        false,
    )
}

fn linear_molecule_block(
    id: &str,
    label: &str,
    features: Vec<DnaFeatureCartoon>,
    left_end: DnaEndStyle,
    right_end: DnaEndStyle,
) -> DnaMoleculeCartoon {
    DnaMoleculeCartoon {
        id: id.to_string(),
        label: label.to_string(),
        topology: DnaTopologyCartoon::Linear,
        features,
        left_end: Some(left_end),
        right_end: Some(right_end),
    }
}

fn event_block(
    id: &str,
    title: &str,
    caption: &str,
    action: ProtocolCartoonAction,
    molecules: Vec<DnaMoleculeCartoon>,
) -> ProtocolCartoonEvent {
    ProtocolCartoonEvent {
        id: id.to_string(),
        title: title.to_string(),
        caption: caption.to_string(),
        action,
        molecules,
    }
}

fn render_invalid_protocol_svg(id: &str, message: &str) -> String {
    format!(
        "<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"900\" height=\"220\" viewBox=\"0 0 900 220\"><rect x=\"0\" y=\"0\" width=\"900\" height=\"220\" fill=\"#fff7ed\"/><text x=\"24\" y=\"56\" font-family=\"Trebuchet MS, Segoe UI, sans-serif\" font-size=\"30\" font-weight=\"700\" fill=\"#8a2f1f\">Invalid protocol cartoon: {}</text><text x=\"24\" y=\"96\" font-family=\"Trebuchet MS, Segoe UI, sans-serif\" font-size=\"18\" fill=\"#7a4a3b\">{}</text></svg>",
        escape_xml(id),
        escape_xml(message)
    )
}

fn protocol_cartoon_template_from_spec(spec: ProtocolCartoonSpec) -> ProtocolCartoonTemplate {
    ProtocolCartoonTemplate {
        schema: protocol_cartoon_template_schema_id(),
        id: spec.id,
        title: spec.title,
        summary: spec.summary,
        defaults: ProtocolCartoonTemplateDefaults::default(),
        events: spec
            .events
            .into_iter()
            .map(|event| ProtocolCartoonTemplateEvent {
                id: event.id,
                title: event.title,
                caption: event.caption,
                action: Some(event.action),
                molecules: event
                    .molecules
                    .into_iter()
                    .map(|molecule| ProtocolCartoonTemplateMolecule {
                        id: molecule.id,
                        label: molecule.label,
                        topology: Some(molecule.topology),
                        features: molecule
                            .features
                            .into_iter()
                            .map(|feature| ProtocolCartoonTemplateFeature {
                                id: feature.id,
                                label: feature.label,
                                length_bp: Some(feature.length_bp),
                                top_length_bp: Some(feature.top_length_bp),
                                bottom_length_bp: Some(feature.bottom_length_bp),
                                color_hex: Some(feature.color_hex),
                                bottom_color_hex: Some(feature.bottom_color_hex),
                                top_nick_after: Some(feature.top_nick_after),
                                bottom_nick_after: Some(feature.bottom_nick_after),
                            })
                            .collect(),
                        left_end: molecule.left_end,
                        right_end: molecule.right_end,
                    })
                    .collect(),
            })
            .collect(),
    }
}

fn gibson_two_fragment_template() -> ProtocolCartoonTemplate {
    protocol_cartoon_template_from_spec(gibson_two_fragment_spec())
}

fn gibson_single_insert_dual_junction_template() -> ProtocolCartoonTemplate {
    protocol_cartoon_template_from_spec(gibson_single_insert_dual_junction_spec())
}

fn pcr_assay_pair_template() -> ProtocolCartoonTemplate {
    protocol_cartoon_template_from_spec(pcr_assay_pair_spec())
}

fn pcr_assay_pair_no_product_template() -> ProtocolCartoonTemplate {
    protocol_cartoon_template_from_spec(pcr_assay_pair_no_product_spec())
}

fn pcr_assay_qpcr_template() -> ProtocolCartoonTemplate {
    protocol_cartoon_template_from_spec(pcr_assay_qpcr_spec())
}

fn pcr_template_context_molecule(
    id: &str,
    label: &str,
    left_context_bp: usize,
    roi_bp: usize,
    right_context_bp: usize,
    context_color: &str,
    roi_color: &str,
) -> DnaMoleculeCartoon {
    linear_molecule_block(
        id,
        label,
        vec![
            duplex_feature_block(
                "template_left_context",
                "Template context",
                left_context_bp,
                context_color,
            ),
            duplex_feature_block("roi", "Selected ROI", roi_bp, roi_color),
            duplex_feature_block(
                "template_right_context",
                "Template context",
                right_context_bp,
                context_color,
            ),
        ],
        DnaEndStyle::Continuation,
        DnaEndStyle::Continuation,
    )
}

fn pcr_constraint_windows_molecule(
    id: &str,
    label: &str,
    left_context_bp: usize,
    forward_window_bp: usize,
    roi_bp: usize,
    reverse_window_bp: usize,
    right_context_bp: usize,
    context_color: &str,
    forward_window_color: &str,
    roi_color: &str,
    reverse_window_color: &str,
) -> DnaMoleculeCartoon {
    linear_molecule_block(
        id,
        label,
        vec![
            duplex_feature_block(
                "template_left_context",
                "Template context",
                left_context_bp,
                context_color,
            ),
            duplex_feature_block(
                "forward_window",
                "Forward primer window",
                forward_window_bp,
                forward_window_color,
            ),
            duplex_feature_block("roi", "Selected ROI", roi_bp, roi_color),
            duplex_feature_block(
                "reverse_window",
                "Reverse primer window",
                reverse_window_bp,
                reverse_window_color,
            ),
            duplex_feature_block(
                "template_right_context",
                "Template context",
                right_context_bp,
                context_color,
            ),
        ],
        DnaEndStyle::Continuation,
        DnaEndStyle::Continuation,
    )
}

fn pcr_primers_within_windows_molecule(
    id: &str,
    label: &str,
    left_context_bp: usize,
    forward_window_margin_bp: usize,
    primer_site_bp: usize,
    roi_bp: usize,
    reverse_window_margin_bp: usize,
    right_context_bp: usize,
    context_color: &str,
    forward_window_color: &str,
    roi_color: &str,
    reverse_window_color: &str,
) -> DnaMoleculeCartoon {
    linear_molecule_block(
        id,
        label,
        vec![
            duplex_feature_block(
                "template_left_context",
                "Template context",
                left_context_bp,
                context_color,
            ),
            duplex_feature_block(
                "forward_window_margin",
                "Forward window margin",
                forward_window_margin_bp,
                forward_window_color,
            ),
            duplex_feature_block(
                "forward_primer_site",
                "Forward primer",
                primer_site_bp,
                forward_window_color,
            ),
            duplex_feature_block("amplicon_roi", "Selected ROI", roi_bp, roi_color),
            duplex_feature_block(
                "reverse_primer_site",
                "Reverse primer",
                primer_site_bp,
                reverse_window_color,
            ),
            duplex_feature_block(
                "reverse_window_margin",
                "Reverse window margin",
                reverse_window_margin_bp,
                reverse_window_color,
            ),
            duplex_feature_block(
                "template_right_context",
                "Template context",
                right_context_bp,
                context_color,
            ),
        ],
        DnaEndStyle::Continuation,
        DnaEndStyle::Continuation,
    )
}

fn pcr_roi_bounded_amplicon_molecule(
    id: &str,
    label: &str,
    primer_site_bp: usize,
    roi_bp: usize,
    forward_primer_color: &str,
    roi_color: &str,
    reverse_primer_color: &str,
) -> DnaMoleculeCartoon {
    linear_molecule_block(
        id,
        label,
        vec![
            duplex_feature_block(
                "forward_primer_site",
                "Forward primer",
                primer_site_bp,
                forward_primer_color,
            ),
            duplex_feature_block("amplicon_roi", "Selected ROI", roi_bp, roi_color),
            duplex_feature_block(
                "reverse_primer_site",
                "Reverse primer",
                primer_site_bp,
                reverse_primer_color,
            ),
        ],
        DnaEndStyle::Blunt,
        DnaEndStyle::Blunt,
    )
}

fn qpcr_constraint_windows_molecule(
    id: &str,
    label: &str,
    left_context_bp: usize,
    forward_window_bp: usize,
    roi_left_bp: usize,
    probe_window_bp: usize,
    roi_right_bp: usize,
    reverse_window_bp: usize,
    right_context_bp: usize,
    context_color: &str,
    forward_window_color: &str,
    roi_color: &str,
    probe_window_color: &str,
    reverse_window_color: &str,
) -> DnaMoleculeCartoon {
    linear_molecule_block(
        id,
        label,
        vec![
            duplex_feature_block(
                "template_left_context",
                "Template context",
                left_context_bp,
                context_color,
            ),
            duplex_feature_block(
                "forward_window",
                "Upstream primer window",
                forward_window_bp,
                forward_window_color,
            ),
            duplex_feature_block("roi_left", "ROI segment", roi_left_bp, roi_color),
            duplex_feature_block(
                "probe_window",
                "Probe window",
                probe_window_bp,
                probe_window_color,
            ),
            duplex_feature_block("roi_right", "ROI segment", roi_right_bp, roi_color),
            duplex_feature_block(
                "reverse_window",
                "Downstream primer window",
                reverse_window_bp,
                reverse_window_color,
            ),
            duplex_feature_block(
                "template_right_context",
                "Template context",
                right_context_bp,
                context_color,
            ),
        ],
        DnaEndStyle::Continuation,
        DnaEndStyle::Continuation,
    )
}

fn qpcr_assay_layout_molecule(
    id: &str,
    label: &str,
    primer_site_bp: usize,
    forward_window_margin_bp: usize,
    roi_left_bp: usize,
    probe_window_left_margin_bp: usize,
    probe_site_bp: usize,
    probe_window_right_margin_bp: usize,
    roi_right_bp: usize,
    reverse_window_margin_bp: usize,
    forward_window_color: &str,
    roi_color: &str,
    probe_window_color: &str,
    reverse_window_color: &str,
) -> DnaMoleculeCartoon {
    linear_molecule_block(
        id,
        label,
        vec![
            duplex_feature_block(
                "forward_primer_site",
                "Forward primer",
                primer_site_bp,
                forward_window_color,
            ),
            duplex_feature_block(
                "forward_window_margin",
                "Forward window margin",
                forward_window_margin_bp,
                forward_window_color,
            ),
            duplex_feature_block("roi_left", "ROI segment", roi_left_bp, roi_color),
            duplex_feature_block(
                "probe_window_left",
                "Probe window",
                probe_window_left_margin_bp,
                probe_window_color,
            ),
            duplex_feature_block(
                "probe_site",
                "Probe oligo",
                probe_site_bp,
                probe_window_color,
            ),
            duplex_feature_block(
                "probe_window_right",
                "Probe window",
                probe_window_right_margin_bp,
                probe_window_color,
            ),
            duplex_feature_block("roi_right", "ROI segment", roi_right_bp, roi_color),
            duplex_feature_block(
                "reverse_window_margin",
                "Reverse window margin",
                reverse_window_margin_bp,
                reverse_window_color,
            ),
            duplex_feature_block(
                "reverse_primer_site",
                "Reverse primer",
                primer_site_bp,
                reverse_window_color,
            ),
        ],
        DnaEndStyle::Blunt,
        DnaEndStyle::Blunt,
    )
}

fn qpcr_quantified_amplicon_molecule(
    id: &str,
    label: &str,
    primer_site_bp: usize,
    forward_window_margin_bp: usize,
    roi_left_bp: usize,
    probe_window_left_margin_bp: usize,
    probe_site_bp: usize,
    probe_window_right_margin_bp: usize,
    roi_right_bp: usize,
    reverse_window_margin_bp: usize,
    forward_window_color: &str,
    roi_color: &str,
    probe_window_color: &str,
    reverse_window_color: &str,
) -> DnaMoleculeCartoon {
    linear_molecule_block(
        id,
        label,
        vec![
            duplex_feature_block(
                "forward_primer_site",
                "Forward primer",
                primer_site_bp,
                forward_window_color,
            ),
            duplex_feature_block(
                "forward_window_margin",
                "Forward window margin",
                forward_window_margin_bp,
                forward_window_color,
            ),
            duplex_feature_block("roi_left", "ROI segment", roi_left_bp, roi_color),
            duplex_feature_block(
                "probe_window_left",
                "Probe window",
                probe_window_left_margin_bp,
                probe_window_color,
            ),
            duplex_feature_block(
                "probe_site",
                "Probe oligo",
                probe_site_bp,
                probe_window_color,
            ),
            duplex_feature_block(
                "probe_window_right",
                "Probe window",
                probe_window_right_margin_bp,
                probe_window_color,
            ),
            duplex_feature_block("roi_right", "ROI segment", roi_right_bp, roi_color),
            duplex_feature_block(
                "reverse_window_margin",
                "Reverse window margin",
                reverse_window_margin_bp,
                reverse_window_color,
            ),
            duplex_feature_block(
                "reverse_primer_site",
                "Reverse primer",
                primer_site_bp,
                reverse_window_color,
            ),
        ],
        DnaEndStyle::Blunt,
        DnaEndStyle::Blunt,
    )
}

fn pcr_assay_pair_spec() -> ProtocolCartoonSpec {
    const TEMPLATE_LEFT_BP: usize = 74;
    const FORWARD_WINDOW_BP: usize = 30;
    const PRIMER_SITE_BP: usize = 18;
    const ROI_BP: usize = 120;
    const REVERSE_WINDOW_BP: usize = 30;
    const TEMPLATE_RIGHT_BP: usize = 74;
    const WINDOW_MARGIN_BP: usize = 12;

    let template_color = "#c9b37e";
    let roi_color = "#2a9d8f";
    let forward_window_color = "#d1495b";
    let reverse_window_color = "#3d6fb6";

    ProtocolCartoonSpec {
        id: ProtocolCartoonKind::PcrAssayPair.id().to_string(),
        title: "GENtle Protocol Cartoon: PCR Assay (pair-primer baseline)".to_string(),
        summary: "Four-step PCR strip: context+ROI, primer constraint windows, primer placement, and ROI-bounded amplicon outcome.".to_string(),
        events: vec![
            event_block(
                "context_roi",
                "Context + ROI",
                "Template context selected.\nROI selected.\nGreen = ROI.",
                ProtocolCartoonAction::Context,
                vec![pcr_template_context_molecule(
                    "template_context",
                    "Source template",
                    TEMPLATE_LEFT_BP,
                    ROI_BP,
                    TEMPLATE_RIGHT_BP,
                    template_color,
                    roi_color,
                )],
            ),
            event_block(
                "primer_constraints",
                "Primer Constraints",
                "Red = Upstream primer window for primer.\nBlue = Downstream primer window for primer.\nPrimer termini must stay inside each window.",
                ProtocolCartoonAction::Custom {
                    label: "Constraints".to_string(),
                },
                vec![pcr_constraint_windows_molecule(
                    "template_with_windows",
                    "Template with ROI + primer windows",
                    TEMPLATE_LEFT_BP,
                    FORWARD_WINDOW_BP,
                    ROI_BP,
                    REVERSE_WINDOW_BP,
                    TEMPLATE_RIGHT_BP,
                    template_color,
                    forward_window_color,
                    roi_color,
                    reverse_window_color,
                )],
            ),
            event_block(
                "primers",
                "Primers",
                "Chosen forward/reverse primers lie within their windows and terminate at the ROI boundaries while genomic context remains visible.",
                ProtocolCartoonAction::Custom {
                    label: "Primers".to_string(),
                },
                vec![pcr_primers_within_windows_molecule(
                    "primer_layout",
                    "Primer placement in windows",
                    TEMPLATE_LEFT_BP,
                    WINDOW_MARGIN_BP,
                    PRIMER_SITE_BP,
                    ROI_BP,
                    WINDOW_MARGIN_BP,
                    TEMPLATE_RIGHT_BP,
                    template_color,
                    forward_window_color,
                    roi_color,
                    reverse_window_color,
                )],
            ),
            event_block(
                "product",
                "Product",
                "PCR yields an ROI-focused amplicon bounded by the chosen primer termini.",
                ProtocolCartoonAction::Custom {
                    label: "Product".to_string(),
                },
                vec![pcr_roi_bounded_amplicon_molecule(
                    "accepted_amplicon",
                    "ROI-bounded amplicon",
                    PRIMER_SITE_BP,
                    ROI_BP,
                    forward_window_color,
                    roi_color,
                    reverse_window_color,
                )],
            ),
        ],
    }
}

fn pcr_assay_pair_no_product_spec() -> ProtocolCartoonSpec {
    const TEMPLATE_LEFT_BP: usize = 74;
    const FORWARD_WINDOW_BP: usize = 30;
    const PRIMER_SITE_BP: usize = 18;
    const ROI_BP: usize = 120;
    const REVERSE_WINDOW_BP: usize = 30;
    const TEMPLATE_RIGHT_BP: usize = 74;
    const WINDOW_MARGIN_BP: usize = 12;

    let template_color = "#c9b37e";
    let roi_color = "#2a9d8f";
    let forward_window_color = "#d1495b";
    let reverse_window_color = "#3d6fb6";

    ProtocolCartoonSpec {
        id: ProtocolCartoonKind::PcrAssayPairNoProduct.id().to_string(),
        title: "GENtle Protocol Cartoon: PCR Assay (report-only no-product)".to_string(),
        summary: "Four-step PCR strip matching the baseline flow, but ending with report-only output when no primer pair is accepted.".to_string(),
        events: vec![
            event_block(
                "context_roi",
                "Context + ROI",
                "Template context selected.\nROI selected.\nGreen = ROI.",
                ProtocolCartoonAction::Context,
                vec![pcr_template_context_molecule(
                    "template_context",
                    "Source template",
                    TEMPLATE_LEFT_BP,
                    ROI_BP,
                    TEMPLATE_RIGHT_BP,
                    template_color,
                    roi_color,
                )],
            ),
            event_block(
                "primer_constraints",
                "Primer Constraints",
                "Red = Upstream primer window for primer.\nBlue = Downstream primer window for primer.\nConstraints may still reject all candidates.",
                ProtocolCartoonAction::Custom {
                    label: "Constraints".to_string(),
                },
                vec![pcr_constraint_windows_molecule(
                    "template_with_windows",
                    "Template with ROI + primer windows",
                    TEMPLATE_LEFT_BP,
                    FORWARD_WINDOW_BP,
                    ROI_BP,
                    REVERSE_WINDOW_BP,
                    TEMPLATE_RIGHT_BP,
                    template_color,
                    forward_window_color,
                    roi_color,
                    reverse_window_color,
                )],
            ),
            event_block(
                "primers",
                "Primers",
                "Candidate forward/reverse primers lie in-window and terminate at ROI boundaries, but no candidate pair survives all constraints.",
                ProtocolCartoonAction::Custom {
                    label: "Primers".to_string(),
                },
                vec![pcr_primers_within_windows_molecule(
                    "primer_layout",
                    "Candidate primer placement",
                    TEMPLATE_LEFT_BP,
                    WINDOW_MARGIN_BP,
                    PRIMER_SITE_BP,
                    ROI_BP,
                    WINDOW_MARGIN_BP,
                    TEMPLATE_RIGHT_BP,
                    template_color,
                    forward_window_color,
                    roi_color,
                    reverse_window_color,
                )],
            ),
            event_block(
                "report_only",
                "Report Only",
                "No amplicon product is emitted; users inspect retained ROI context and the failure report.",
                ProtocolCartoonAction::Custom {
                    label: "Report".to_string(),
                },
                vec![pcr_template_context_molecule(
                    "template_context_report_only",
                    "Template + ROI retained",
                    TEMPLATE_LEFT_BP,
                    ROI_BP,
                    TEMPLATE_RIGHT_BP,
                    template_color,
                    roi_color,
                )],
            ),
        ],
    }
}

fn pcr_assay_qpcr_spec() -> ProtocolCartoonSpec {
    const TEMPLATE_LEFT_BP: usize = 62;
    const FORWARD_WINDOW_BP: usize = 26;
    const PRIMER_SITE_BP: usize = 18;
    const ROI_LEFT_BP: usize = 40;
    const PROBE_WINDOW_BP: usize = 40;
    const PROBE_SITE_BP: usize = 12;
    const ROI_RIGHT_BP: usize = 40;
    const REVERSE_WINDOW_BP: usize = 26;
    const TEMPLATE_RIGHT_BP: usize = 62;
    const FORWARD_WINDOW_MARGIN_BP: usize = FORWARD_WINDOW_BP - PRIMER_SITE_BP;
    const REVERSE_WINDOW_MARGIN_BP: usize = REVERSE_WINDOW_BP - PRIMER_SITE_BP;
    const PROBE_WINDOW_LEFT_MARGIN_BP: usize = (PROBE_WINDOW_BP - PROBE_SITE_BP) / 2;
    const PROBE_WINDOW_RIGHT_MARGIN_BP: usize =
        PROBE_WINDOW_BP - PROBE_SITE_BP - PROBE_WINDOW_LEFT_MARGIN_BP;
    const ROI_TOTAL_BP: usize = ROI_LEFT_BP + PROBE_WINDOW_BP + ROI_RIGHT_BP;

    let template_color = "#c9b37e";
    let roi_color = "#2a9d8f";
    let probe_window_color = "#577590";
    let forward_window_color = "#d1495b";
    let reverse_window_color = "#3d6fb6";

    ProtocolCartoonSpec {
        id: ProtocolCartoonKind::PcrAssayQpcr.id().to_string(),
        title: "GENtle Protocol Cartoon: qPCR Assay (probe-bearing baseline)".to_string(),
        summary: "Four-step qPCR strip aligned to the baseline PCR flow, with an explicit probe window independent of ROI.".to_string(),
        events: vec![
            event_block(
                "context_roi",
                "Context + ROI",
                "Template context selected.\nROI selected.\nGreen = ROI.",
                ProtocolCartoonAction::Context,
                vec![pcr_template_context_molecule(
                    "template_context",
                    "Source template",
                    TEMPLATE_LEFT_BP,
                    ROI_TOTAL_BP,
                    TEMPLATE_RIGHT_BP,
                    template_color,
                    roi_color,
                )],
            ),
            event_block(
                "primer_constraints",
                "Primer Constraints",
                "Red = Upstream primer window.\nBlue = Downstream primer window.\nSlate = Internal probe window.",
                ProtocolCartoonAction::Custom {
                    label: "Constraints".to_string(),
                },
                vec![qpcr_constraint_windows_molecule(
                    "template_with_qpcr_windows",
                    "Template with ROI + primer/probe windows",
                    TEMPLATE_LEFT_BP,
                    FORWARD_WINDOW_BP,
                    ROI_LEFT_BP,
                    PROBE_WINDOW_BP,
                    ROI_RIGHT_BP,
                    REVERSE_WINDOW_BP,
                    TEMPLATE_RIGHT_BP,
                    template_color,
                    forward_window_color,
                    roi_color,
                    probe_window_color,
                    reverse_window_color,
                )],
            ),
            event_block(
                "primers_probe",
                "Primers + Probe",
                "Amplified amplicon only.\nNo genomic flanks shown.\nOligos remain within their windows.",
                ProtocolCartoonAction::Custom {
                    label: "Assay".to_string(),
                },
                vec![qpcr_assay_layout_molecule(
                    "qpcr_layout",
                    "Amplified amplicon with assay oligos",
                    PRIMER_SITE_BP,
                    FORWARD_WINDOW_MARGIN_BP,
                    ROI_LEFT_BP,
                    PROBE_WINDOW_LEFT_MARGIN_BP,
                    PROBE_SITE_BP,
                    PROBE_WINDOW_RIGHT_MARGIN_BP,
                    ROI_RIGHT_BP,
                    REVERSE_WINDOW_MARGIN_BP,
                    forward_window_color,
                    roi_color,
                    probe_window_color,
                    reverse_window_color,
                )],
            ),
            event_block(
                "quantify",
                "Quantify",
                "Fluorescence is read per cycle.\nSame probe-bearing amplicon geometry.",
                ProtocolCartoonAction::Custom {
                    label: "Quantify".to_string(),
                },
                vec![qpcr_quantified_amplicon_molecule(
                    "quantified_amplicon",
                    "Quantified probe-bearing amplicon",
                    PRIMER_SITE_BP,
                    FORWARD_WINDOW_MARGIN_BP,
                    ROI_LEFT_BP,
                    PROBE_WINDOW_LEFT_MARGIN_BP,
                    PROBE_SITE_BP,
                    PROBE_WINDOW_RIGHT_MARGIN_BP,
                    ROI_RIGHT_BP,
                    REVERSE_WINDOW_MARGIN_BP,
                    forward_window_color,
                    roi_color,
                    probe_window_color,
                    reverse_window_color,
                )],
            ),
        ],
    }
}

fn gibson_single_insert_dual_junction_spec() -> ProtocolCartoonSpec {
    const DEST_BODY_BP: usize = 54;
    const INSERT_BODY_BP: usize = 52;
    const LEFT_OVERLAP_BP: usize = 24;
    const RIGHT_OVERLAP_BP: usize = 24;
    const PENDING_FILL_BP: usize = 10;

    let destination_color = "#f2c94c";
    let insert_color = "#3f80e0";
    let left_junction_color = "#55a84f";
    let right_junction_color = "#f2994a";

    ProtocolCartoonSpec {
        id: ProtocolCartoonKind::GibsonSingleInsertDualJunction
            .id()
            .to_string(),
        title: "GENtle Protocol Cartoon: Gibson Single-Insert Assembly".to_string(),
        summary: "Five-step Gibson mechanism with two explicit destination-insert junctions."
            .to_string(),
        events: vec![
            event_block(
                "context",
                "Context",
                "The opened destination contributes one left arm and one right arm, while the insert carries one matching terminal homology at each end.",
                ProtocolCartoonAction::Context,
                vec![
                    linear_molecule_block(
                        "destination_left_context",
                        "Destination left arm",
                        vec![
                            duplex_feature_block(
                                "dest_left_body",
                                "Destination left flank",
                                DEST_BODY_BP,
                                destination_color,
                            ),
                            duplex_feature_block(
                                "left_overlap",
                                "Left junction homology",
                                LEFT_OVERLAP_BP,
                                left_junction_color,
                            ),
                        ],
                        DnaEndStyle::Continuation,
                        DnaEndStyle::Blunt,
                    ),
                    linear_molecule_block(
                        "insert_context",
                        "Insert",
                        vec![
                            duplex_feature_block(
                                "insert_left_overlap",
                                "Left junction homology",
                                LEFT_OVERLAP_BP,
                                left_junction_color,
                            ),
                            duplex_feature_block(
                                "insert_body",
                                "Insert body",
                                INSERT_BODY_BP,
                                insert_color,
                            ),
                            duplex_feature_block(
                                "insert_right_overlap",
                                "Right junction homology",
                                RIGHT_OVERLAP_BP,
                                right_junction_color,
                            ),
                        ],
                        DnaEndStyle::Blunt,
                        DnaEndStyle::Blunt,
                    ),
                    linear_molecule_block(
                        "destination_right_context",
                        "Destination right arm",
                        vec![
                            duplex_feature_block(
                                "right_overlap",
                                "Right junction homology",
                                RIGHT_OVERLAP_BP,
                                right_junction_color,
                            ),
                            duplex_feature_block(
                                "dest_right_body",
                                "Destination right flank",
                                DEST_BODY_BP,
                                destination_color,
                            ),
                        ],
                        DnaEndStyle::Blunt,
                        DnaEndStyle::Continuation,
                    ),
                ],
            ),
            event_block(
                "chew_back",
                "Chew-back",
                "A 5' exonuclease exposes one single-stranded 3' overhang at each destination arm and both ends of the insert, continuing slightly past each homology region so polymerase still has work to do.",
                ProtocolCartoonAction::Custom {
                    label: "5' Exonuclease".to_string(),
                },
                vec![
                    linear_molecule_block(
                        "destination_left_chewed",
                        "Left arm (chewed)",
                        vec![
                            duplex_feature_block(
                                "dest_left_body",
                                "Destination left flank",
                                DEST_BODY_BP,
                                destination_color,
                            ),
                            top_only_feature_block(
                                "dest_left_tail_exposed",
                                "Left-side gap to be filled",
                                PENDING_FILL_BP,
                                destination_color,
                            ),
                            top_only_feature_block(
                                "left_overlap_exposed",
                                "Left junction homology",
                                LEFT_OVERLAP_BP,
                                left_junction_color,
                            ),
                        ],
                        DnaEndStyle::Continuation,
                        DnaEndStyle::Sticky {
                            polarity: OverhangPolarity::ThreePrime,
                            nt: LEFT_OVERLAP_BP + PENDING_FILL_BP,
                        },
                    ),
                    linear_molecule_block(
                        "insert_chewed",
                        "Insert (chewed)",
                        vec![
                            bottom_only_feature_block(
                                "insert_left_overlap_exposed",
                                "Left junction homology",
                                LEFT_OVERLAP_BP,
                                left_junction_color,
                            ),
                            bottom_only_feature_block(
                                "insert_left_tail_exposed",
                                "Insert left-side gap",
                                PENDING_FILL_BP,
                                insert_color,
                            ),
                            duplex_feature_block(
                                "insert_body",
                                "Insert body",
                                INSERT_BODY_BP,
                                insert_color,
                            ),
                            top_only_feature_block(
                                "insert_right_tail_exposed",
                                "Insert right-side gap",
                                PENDING_FILL_BP,
                                insert_color,
                            ),
                            top_only_feature_block(
                                "insert_right_overlap_exposed",
                                "Right junction homology",
                                RIGHT_OVERLAP_BP,
                                right_junction_color,
                            ),
                        ],
                        DnaEndStyle::Sticky {
                            polarity: OverhangPolarity::ThreePrime,
                            nt: LEFT_OVERLAP_BP + PENDING_FILL_BP,
                        },
                        DnaEndStyle::Sticky {
                            polarity: OverhangPolarity::ThreePrime,
                            nt: RIGHT_OVERLAP_BP + PENDING_FILL_BP,
                        },
                    ),
                    linear_molecule_block(
                        "destination_right_chewed",
                        "Right arm (chewed)",
                        vec![
                            bottom_only_feature_block(
                                "right_overlap_exposed",
                                "Right junction homology",
                                RIGHT_OVERLAP_BP,
                                right_junction_color,
                            ),
                            bottom_only_feature_block(
                                "dest_right_tail_exposed",
                                "Right-side gap to be filled",
                                PENDING_FILL_BP,
                                destination_color,
                            ),
                            duplex_feature_block(
                                "dest_right_body",
                                "Destination right flank",
                                DEST_BODY_BP,
                                destination_color,
                            ),
                        ],
                        DnaEndStyle::Sticky {
                            polarity: OverhangPolarity::ThreePrime,
                            nt: RIGHT_OVERLAP_BP + PENDING_FILL_BP,
                        },
                        DnaEndStyle::Continuation,
                    ),
                ],
            ),
            event_block(
                "anneal",
                "Anneal",
                "The insert anneals to the left and right destination arms at both junctions while short single-stranded unique tails remain on both sides of the paired overlaps.",
                ProtocolCartoonAction::Anneal,
                vec![linear_molecule_block(
                    "annealed_intermediate",
                    "Annealed intermediate",
                    vec![
                        duplex_feature_block(
                            "dest_left_body",
                            "Destination left flank",
                            DEST_BODY_BP,
                            destination_color,
                        ),
                        top_only_feature_block(
                            "dest_left_tail",
                            "Left-side gap to be filled",
                            PENDING_FILL_BP,
                            destination_color,
                        ),
                        duplex_feature_block(
                            "left_overlap",
                            "Left junction homology",
                            LEFT_OVERLAP_BP,
                            left_junction_color,
                        ),
                        bottom_only_feature_block(
                            "insert_left_tail",
                            "Insert left-side gap",
                            PENDING_FILL_BP,
                            insert_color,
                        ),
                        duplex_feature_block(
                            "insert_body",
                            "Insert body",
                            INSERT_BODY_BP,
                            insert_color,
                        ),
                        top_only_feature_block(
                            "insert_right_tail",
                            "Insert right-side gap",
                            PENDING_FILL_BP,
                            insert_color,
                        ),
                        duplex_feature_block(
                            "right_overlap",
                            "Right junction homology",
                            RIGHT_OVERLAP_BP,
                            right_junction_color,
                        ),
                        bottom_only_feature_block(
                            "dest_right_tail",
                            "Right-side gap to be filled",
                            PENDING_FILL_BP,
                            destination_color,
                        ),
                        duplex_feature_block(
                            "dest_right_body",
                            "Destination right flank",
                            DEST_BODY_BP,
                            destination_color,
                        ),
                    ],
                    DnaEndStyle::Continuation,
                    DnaEndStyle::Continuation,
                )],
            ),
            event_block(
                "extend",
                "Extend",
                "DNA polymerase fills both remaining gaps, leaving one nick on each strand at the two repaired junctions.",
                ProtocolCartoonAction::Extend,
                vec![linear_molecule_block(
                    "extended_intermediate",
                    "Extended intermediate",
                    vec![
                        duplex_feature_block_with_nicks(
                            "dest_left_body",
                            "Destination left flank",
                            DEST_BODY_BP,
                            destination_color,
                            false,
                            true,
                        ),
                        duplex_feature_block(
                            "dest_left_tail",
                            "Left repaired gap",
                            PENDING_FILL_BP,
                            destination_color,
                        ),
                        duplex_feature_block(
                            "left_overlap",
                            "Left junction homology",
                            LEFT_OVERLAP_BP,
                            left_junction_color,
                        ),
                        duplex_feature_block_with_nicks(
                            "insert_left_tail",
                            "Insert left repaired gap",
                            PENDING_FILL_BP,
                            insert_color,
                            true,
                            false,
                        ),
                        duplex_feature_block_with_nicks(
                            "insert_body",
                            "Insert body",
                            INSERT_BODY_BP,
                            insert_color,
                            false,
                            true,
                        ),
                        duplex_feature_block(
                            "insert_right_tail",
                            "Insert right repaired gap",
                            PENDING_FILL_BP,
                            insert_color,
                        ),
                        duplex_feature_block(
                            "right_overlap",
                            "Right junction homology",
                            RIGHT_OVERLAP_BP,
                            right_junction_color,
                        ),
                        duplex_feature_block_with_nicks(
                            "dest_right_tail",
                            "Right repaired gap",
                            PENDING_FILL_BP,
                            destination_color,
                            true,
                            false,
                        ),
                        duplex_feature_block(
                            "dest_right_body",
                            "Destination right flank",
                            DEST_BODY_BP,
                            destination_color,
                        ),
                    ],
                    DnaEndStyle::Continuation,
                    DnaEndStyle::Continuation,
                )],
            ),
            event_block(
                "seal",
                "Seal",
                "DNA ligase seals both remaining nicks, leaving one continuous destination-insert-destination duplex.",
                ProtocolCartoonAction::Seal,
                vec![linear_molecule_block(
                    "sealed_product",
                    "Sealed product",
                    vec![
                        duplex_feature_block(
                            "dest_left_body",
                            "Destination left flank",
                            DEST_BODY_BP,
                            destination_color,
                        ),
                        duplex_feature_block(
                            "dest_left_tail",
                            "Left repaired gap",
                            PENDING_FILL_BP,
                            destination_color,
                        ),
                        duplex_feature_block(
                            "left_overlap",
                            "Left junction homology",
                            LEFT_OVERLAP_BP,
                            left_junction_color,
                        ),
                        duplex_feature_block(
                            "insert_left_tail",
                            "Insert left repaired gap",
                            PENDING_FILL_BP,
                            insert_color,
                        ),
                        duplex_feature_block(
                            "insert_body",
                            "Insert body",
                            INSERT_BODY_BP,
                            insert_color,
                        ),
                        duplex_feature_block(
                            "insert_right_tail",
                            "Insert right repaired gap",
                            PENDING_FILL_BP,
                            insert_color,
                        ),
                        duplex_feature_block(
                            "right_overlap",
                            "Right junction homology",
                            RIGHT_OVERLAP_BP,
                            right_junction_color,
                        ),
                        duplex_feature_block(
                            "dest_right_tail",
                            "Right repaired gap",
                            PENDING_FILL_BP,
                            destination_color,
                        ),
                        duplex_feature_block(
                            "dest_right_body",
                            "Destination right flank",
                            DEST_BODY_BP,
                            destination_color,
                        ),
                    ],
                    DnaEndStyle::Continuation,
                    DnaEndStyle::Continuation,
                )],
            ),
        ],
    }
}

fn gibson_two_fragment_spec() -> ProtocolCartoonSpec {
    // Keep the homologous overlap visually stable across all panels by using
    // the same displayed canonical total (150 bp) for separate fragments and
    // joined intermediates, while compressing the non-overlap sequence bodies
    // in the later panels.
    const DISPLAY_OVERLAP_BP: usize = 30;
    const DISPLAY_CHEW_TAIL_BP: usize = 12;
    const DISPLAY_JOINED_TAIL_BP: usize = 12;
    const DISPLAY_FRAGMENT_TOTAL_BP: usize = 150;
    const DISPLAY_PRE_JOIN_BODY_BP: usize = DISPLAY_FRAGMENT_TOTAL_BP - DISPLAY_OVERLAP_BP;
    const DISPLAY_CHEW_BODY_BP: usize =
        DISPLAY_FRAGMENT_TOTAL_BP - DISPLAY_OVERLAP_BP - DISPLAY_CHEW_TAIL_BP;
    const DISPLAY_JOINED_BODY_BP: usize =
        (DISPLAY_FRAGMENT_TOTAL_BP - DISPLAY_OVERLAP_BP - (DISPLAY_JOINED_TAIL_BP * 2)) / 2;
    const DISPLAY_SEALED_BODY_BP: usize = (DISPLAY_FRAGMENT_TOTAL_BP - DISPLAY_OVERLAP_BP) / 2;

    let fragment_a_color = "#f2c94c";
    let fragment_b_color = "#2f80ed";
    let overlap_color = "#55a84f";

    ProtocolCartoonSpec {
        id: "gibson.two_fragment".to_string(),
        title: "GENtle Protocol Cartoon: Gibson Two-Fragment (A+B)".to_string(),
        summary:
            "Five-step Gibson mechanism with origin colors (A=yellow, B=blue) and a 30 bp homologous overlap"
                .to_string(),
        events: vec![
            event_block(
                "context",
                "Context",
                "Fragment A ends with a 30 bp homologous sequence at its 3' end, and fragment B starts with the same 30 bp sequence at its 5' end.",
                ProtocolCartoonAction::Context,
                vec![
                    linear_molecule_block(
                        "fragment_a_context",
                        "A",
                        vec![
                            duplex_feature_block(
                                "a_body",
                                "A body",
                                DISPLAY_PRE_JOIN_BODY_BP,
                                fragment_a_color,
                            ),
                            duplex_feature_block(
                                "a_overlap",
                                "A 30 bp overlap",
                                DISPLAY_OVERLAP_BP,
                                overlap_color,
                            ),
                        ],
                        DnaEndStyle::Continuation,
                        DnaEndStyle::Blunt,
                    ),
                    linear_molecule_block(
                        "fragment_b_context",
                        "B",
                        vec![
                            duplex_feature_block(
                                "b_overlap",
                                "B 30 bp overlap",
                                DISPLAY_OVERLAP_BP,
                                overlap_color,
                            ),
                            duplex_feature_block(
                                "b_body",
                                "B body",
                                DISPLAY_PRE_JOIN_BODY_BP,
                                fragment_b_color,
                            ),
                        ],
                        DnaEndStyle::Blunt,
                        DnaEndStyle::Continuation,
                    ),
                ],
            ),
            event_block(
                "chew_back",
                "Chew-back",
                "A 5' exonuclease chews back both fragments past the homologous ends, exposing the full 30 bp overlap plus a short unique single-stranded tail on each 3' end.",
                ProtocolCartoonAction::Custom {
                    label: "5' Exonuclease".to_string(),
                },
                vec![
                    linear_molecule_block(
                        "fragment_a_chewed",
                        "A (chewed)",
                        vec![
                            duplex_feature_block(
                                "a_body_ds",
                                "A body (duplex)",
                                DISPLAY_CHEW_BODY_BP,
                                fragment_a_color,
                            ),
                            top_only_feature_block(
                                "a_body_ss",
                                "A unique 3' tail",
                                DISPLAY_CHEW_TAIL_BP,
                                fragment_a_color,
                            ),
                            top_only_feature_block(
                                "a_overlap_ss",
                                "A exposed 30 bp overlap",
                                DISPLAY_OVERLAP_BP,
                                overlap_color,
                            ),
                        ],
                        DnaEndStyle::Continuation,
                        DnaEndStyle::Sticky {
                            polarity: OverhangPolarity::ThreePrime,
                            nt: 30,
                        },
                    ),
                    linear_molecule_block(
                        "fragment_b_chewed",
                        "B (chewed)",
                        vec![
                            bottom_only_feature_block(
                                "b_overlap_ss",
                                "B exposed 30 bp overlap",
                                DISPLAY_OVERLAP_BP,
                                overlap_color,
                            ),
                            bottom_only_feature_block(
                                "b_body_ss",
                                "B unique 3' tail",
                                DISPLAY_CHEW_TAIL_BP,
                                fragment_b_color,
                            ),
                            duplex_feature_block(
                                "b_body_ds",
                                "B body (duplex)",
                                DISPLAY_CHEW_BODY_BP,
                                fragment_b_color,
                            ),
                        ],
                        DnaEndStyle::Sticky {
                            polarity: OverhangPolarity::ThreePrime,
                            nt: 30,
                        },
                        DnaEndStyle::Continuation,
                    ),
                ],
            ),
            event_block(
                "anneal",
                "Anneal",
                "The exposed 3' overhangs anneal through the 30 bp homologous region, while the extra yellow and blue sequence-unique tails remain single-stranded.",
                ProtocolCartoonAction::Anneal,
                vec![linear_molecule_block(
                    "annealed_intermediate",
                    "A+B annealed",
                    vec![
                        duplex_feature_block(
                            "a_body_ds",
                            "A body (duplex)",
                            DISPLAY_JOINED_BODY_BP,
                            fragment_a_color,
                        ),
                        top_only_feature_block(
                            "a_body_ss",
                            "A unique single-stranded tail",
                            DISPLAY_JOINED_TAIL_BP,
                            fragment_a_color,
                        ),
                        duplex_feature_block(
                            "hybrid_overlap",
                            "30 bp annealed overlap",
                            DISPLAY_OVERLAP_BP,
                            overlap_color,
                        ),
                        bottom_only_feature_block(
                            "b_body_ss",
                            "B unique single-stranded tail",
                            DISPLAY_JOINED_TAIL_BP,
                            fragment_b_color,
                        ),
                        duplex_feature_block(
                            "b_body_ds",
                            "B body (duplex)",
                            DISPLAY_JOINED_BODY_BP,
                            fragment_b_color,
                        ),
                    ],
                    DnaEndStyle::Continuation,
                    DnaEndStyle::Continuation,
                )],
            ),
            event_block(
                "extend",
                "Extend",
                "DNA polymerase fills the single-stranded tails, restoring duplex DNA but leaving one nick on each strand for ligase to seal.",
                ProtocolCartoonAction::Extend,
                vec![linear_molecule_block(
                    "extended_intermediate",
                    "A+B extended duplex",
                    vec![
                        duplex_feature_block_with_nicks(
                            "a_body_ds",
                            "A body (duplex)",
                            DISPLAY_JOINED_BODY_BP,
                            fragment_a_color,
                            false,
                            true,
                        ),
                        duplex_feature_block(
                            "a_tail_filled",
                            "A tail filled by polymerase",
                            DISPLAY_JOINED_TAIL_BP,
                            fragment_a_color,
                        ),
                        duplex_feature_block(
                            "hybrid_overlap",
                            "Annealed overlap",
                            DISPLAY_OVERLAP_BP,
                            overlap_color,
                        ),
                        duplex_feature_block_with_nicks(
                            "b_tail_filled",
                            "B tail filled by polymerase",
                            DISPLAY_JOINED_TAIL_BP,
                            fragment_b_color,
                            true,
                            false,
                        ),
                        duplex_feature_block(
                            "b_body_ds",
                            "B body (duplex)",
                            DISPLAY_JOINED_BODY_BP,
                            fragment_b_color,
                        ),
                    ],
                    DnaEndStyle::Continuation,
                    DnaEndStyle::Continuation,
                )],
            ),
            event_block(
                "seal",
                "Seal",
                "DNA ligase seals nicks; the assembled duplex is covalently continuous while retaining A/B color origin.",
                ProtocolCartoonAction::Seal,
                vec![linear_molecule_block(
                    "sealed_product",
                    "A+B sealed product",
                    vec![
                        duplex_feature_block(
                            "a_body",
                            "A body",
                            DISPLAY_SEALED_BODY_BP,
                            fragment_a_color,
                        ),
                        duplex_feature_block(
                            "hybrid_overlap",
                            "Scarless overlap junction",
                            DISPLAY_OVERLAP_BP,
                            overlap_color,
                        ),
                        duplex_feature_block(
                            "b_body",
                            "B body",
                            DISPLAY_SEALED_BODY_BP,
                            fragment_b_color,
                        ),
                    ],
                    DnaEndStyle::Continuation,
                    DnaEndStyle::Continuation,
                )],
            ),
        ],
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn custom_test_spec() -> ProtocolCartoonSpec {
        ProtocolCartoonSpec {
            id: "custom.test".to_string(),
            title: "Custom Test".to_string(),
            summary: "Linear/circular plus sticky/blunt coverage".to_string(),
            events: vec![ProtocolCartoonEvent {
                id: "evt1".to_string(),
                title: "Any event".to_string(),
                caption: "cut and separate".to_string(),
                action: ProtocolCartoonAction::Separate,
                molecules: vec![
                    DnaMoleculeCartoon {
                        id: "lin".to_string(),
                        label: "Linear sticky".to_string(),
                        topology: DnaTopologyCartoon::Linear,
                        features: vec![
                            DnaFeatureCartoon {
                                id: "f1".to_string(),
                                label: "f1".to_string(),
                                length_bp: 120,
                                top_length_bp: 120,
                                bottom_length_bp: 120,
                                color_hex: "#004488".to_string(),
                                bottom_color_hex: "#004488".to_string(),
                                top_nick_after: false,
                                bottom_nick_after: false,
                            },
                            DnaFeatureCartoon {
                                id: "f2".to_string(),
                                label: "f2".to_string(),
                                length_bp: 80,
                                top_length_bp: 80,
                                bottom_length_bp: 80,
                                color_hex: "#ee9933".to_string(),
                                bottom_color_hex: "#ee9933".to_string(),
                                top_nick_after: false,
                                bottom_nick_after: false,
                            },
                        ],
                        left_end: Some(DnaEndStyle::Sticky {
                            polarity: OverhangPolarity::FivePrime,
                            nt: 6,
                        }),
                        right_end: Some(DnaEndStyle::Blunt),
                    },
                    DnaMoleculeCartoon {
                        id: "cir".to_string(),
                        label: "Circular".to_string(),
                        topology: DnaTopologyCartoon::Circular,
                        features: vec![DnaFeatureCartoon {
                            id: "c1".to_string(),
                            label: "c1".to_string(),
                            length_bp: 300,
                            top_length_bp: 300,
                            bottom_length_bp: 300,
                            color_hex: "#118833".to_string(),
                            bottom_color_hex: "#118833".to_string(),
                            top_nick_after: false,
                            bottom_nick_after: false,
                        }],
                        left_end: None,
                        right_end: None,
                    },
                ],
            }],
        }
    }

    #[test]
    fn shared_building_blocks_preserve_strand_occupancy_and_nicks() {
        let duplex =
            duplex_feature_block_with_nicks("duplex", "Duplex", 20, "#112233", true, false);
        assert_eq!((duplex.top_length_bp, duplex.bottom_length_bp), (20, 20));
        assert!(duplex.top_nick_after);
        assert!(!duplex.bottom_nick_after);

        let top_only = top_only_feature_block("top", "Top only", 11, "#445566");
        assert_eq!((top_only.top_length_bp, top_only.bottom_length_bp), (11, 0));

        let bottom_only = bottom_only_feature_block("bottom", "Bottom only", 9, "#778899");
        assert_eq!(
            (bottom_only.top_length_bp, bottom_only.bottom_length_bp),
            (0, 9)
        );

        let hybrid = hybrid_feature_block("hybrid", "Hybrid", 15, "#abcdef", "#fedcba");
        assert_eq!((hybrid.top_length_bp, hybrid.bottom_length_bp), (15, 15));
        assert_eq!(hybrid.color_hex, "#abcdef");
        assert_eq!(hybrid.bottom_color_hex, "#fedcba");
    }

    #[test]
    fn template_deserialize_applies_schema_and_defaults() {
        let template: ProtocolCartoonTemplate = serde_json::from_str(
            r#"{
                "id": "demo.protocol",
                "events": [
                    { "id": "e1" }
                ]
            }"#,
        )
        .expect("deserialize template");

        assert_eq!(template.schema, "gentle.protocol_cartoon_template.v1");
        assert_eq!(template.defaults.feature_length_bp, 120);
        assert_eq!(template.defaults.palette.len(), 4);
    }

    #[test]
    fn resolve_sparse_template_uses_reasonable_defaults() {
        let template: ProtocolCartoonTemplate = serde_json::from_str(
            r#"{
                "id": "demo.protocol",
                "events": [
                    { "id": "e1" }
                ]
            }"#,
        )
        .expect("deserialize template");

        let spec = resolve_protocol_cartoon_template(&template).expect("resolve template");
        assert_eq!(spec.id, "demo.protocol");
        assert_eq!(spec.title, "GENtle Protocol Cartoon: demo.protocol");
        assert!(spec.summary.contains("Event-sequence"));
        assert_eq!(spec.events.len(), 1);
        assert_eq!(spec.events[0].action, ProtocolCartoonAction::Context);
        assert_eq!(spec.events[0].molecules.len(), 1);
        assert_eq!(
            spec.events[0].molecules[0].topology,
            DnaTopologyCartoon::Linear
        );
        assert_eq!(
            spec.events[0].molecules[0].left_end,
            Some(DnaEndStyle::Blunt)
        );
        assert_eq!(
            spec.events[0].molecules[0].right_end,
            Some(DnaEndStyle::Blunt)
        );
        assert_eq!(spec.events[0].molecules[0].features.len(), 1);
        assert_eq!(spec.events[0].molecules[0].features[0].length_bp, 120);
        assert_eq!(spec.events[0].molecules[0].features[0].color_hex, "#0f5d75");
    }

    #[test]
    fn resolve_template_rejects_unknown_schema() {
        let template: ProtocolCartoonTemplate = serde_json::from_str(
            r#"{
                "schema": "gentle.protocol_cartoon_template.v0",
                "id": "demo.protocol",
                "events": [
                    { "id": "e1" }
                ]
            }"#,
        )
        .expect("deserialize template");
        let err = resolve_protocol_cartoon_template(&template).expect_err("schema mismatch");
        assert!(err.contains("Unsupported protocol cartoon template schema"));
    }

    #[test]
    fn resolve_template_handles_linear_and_circular_end_defaults() {
        let template = ProtocolCartoonTemplate {
            schema: protocol_cartoon_template_schema_id(),
            id: "template.ends".to_string(),
            title: "".to_string(),
            summary: "".to_string(),
            defaults: ProtocolCartoonTemplateDefaults::default(),
            events: vec![ProtocolCartoonTemplateEvent {
                id: "evt".to_string(),
                title: "evt".to_string(),
                caption: "".to_string(),
                action: Some(ProtocolCartoonAction::Separate),
                molecules: vec![
                    ProtocolCartoonTemplateMolecule {
                        id: "linear".to_string(),
                        label: "linear".to_string(),
                        topology: Some(DnaTopologyCartoon::Linear),
                        features: vec![ProtocolCartoonTemplateFeature::default()],
                        left_end: None,
                        right_end: Some(DnaEndStyle::Sticky {
                            polarity: OverhangPolarity::ThreePrime,
                            nt: 7,
                        }),
                    },
                    ProtocolCartoonTemplateMolecule {
                        id: "circular".to_string(),
                        label: "circular".to_string(),
                        topology: Some(DnaTopologyCartoon::Circular),
                        features: vec![ProtocolCartoonTemplateFeature::default()],
                        left_end: Some(DnaEndStyle::Blunt),
                        right_end: Some(DnaEndStyle::Blunt),
                    },
                ],
            }],
        };

        let spec = resolve_protocol_cartoon_template(&template).expect("resolve template");
        assert_eq!(
            spec.events[0].molecules[0].left_end,
            Some(DnaEndStyle::Blunt)
        );
        assert_eq!(
            spec.events[0].molecules[0].right_end,
            Some(DnaEndStyle::Sticky {
                polarity: OverhangPolarity::ThreePrime,
                nt: 7
            })
        );
        assert_eq!(spec.events[0].molecules[1].left_end, None);
        assert_eq!(spec.events[0].molecules[1].right_end, None);
    }

    #[test]
    fn render_template_svg_resolves_sparse_template() {
        let template: ProtocolCartoonTemplate = serde_json::from_str(
            r#"{
                "id": "demo.protocol",
                "events": [
                    { "id": "e1" }
                ]
            }"#,
        )
        .expect("deserialize template");
        let svg = render_protocol_cartoon_template_svg(&template);
        assert!(svg.contains("data-protocol-id=\"demo.protocol\""));
        assert!(svg.contains("Step 1"));
    }

    #[test]
    fn resolve_template_with_bindings_applies_overrides() {
        let template: ProtocolCartoonTemplate = serde_json::from_str(
            r#"{
                "id": "demo.protocol",
                "events": [
                    {
                        "id": "e1",
                        "title": "Original",
                        "molecules": [
                            {
                                "id": "m1",
                                "label": "Mol",
                                "features": [
                                    { "id": "f1", "label": "Feat", "length_bp": 100 }
                                ]
                            }
                        ]
                    }
                ]
            }"#,
        )
        .expect("deserialize template");
        let bindings: ProtocolCartoonTemplateBindings = serde_json::from_str(
            r##"{
                "template_id": "demo.protocol",
                "defaults": { "feature_length_bp": 75 },
                "event_overrides": [
                    { "event_id": "e1", "title": "Bound Title", "action": "Cut" }
                ],
                "molecule_overrides": [
                    { "event_id": "e1", "molecule_id": "m1", "right_end": { "Sticky": { "polarity": "FivePrime", "nt": 9 } } }
                ],
                "feature_overrides": [
                    { "event_id": "e1", "molecule_id": "m1", "feature_id": "f1", "length_bp": 333, "color_hex": "#112233" }
                ]
            }"##,
        )
        .expect("deserialize bindings");

        let spec = resolve_protocol_cartoon_template_with_bindings(&template, &bindings)
            .expect("resolve template with bindings");
        assert_eq!(spec.events[0].title, "Bound Title");
        assert_eq!(spec.events[0].action, ProtocolCartoonAction::Cut);
        assert_eq!(
            spec.events[0].molecules[0].right_end,
            Some(DnaEndStyle::Sticky {
                polarity: OverhangPolarity::FivePrime,
                nt: 9
            })
        );
        assert_eq!(spec.events[0].molecules[0].features[0].length_bp, 333);
        assert_eq!(spec.events[0].molecules[0].features[0].color_hex, "#112233");
    }

    #[test]
    fn resolve_template_with_bindings_rejects_unknown_schema() {
        let template: ProtocolCartoonTemplate = serde_json::from_str(
            r#"{
                "id": "demo.protocol",
                "events": [{ "id": "e1" }]
            }"#,
        )
        .expect("deserialize template");
        let bindings: ProtocolCartoonTemplateBindings = serde_json::from_str(
            r#"{
                "schema": "gentle.protocol_cartoon_template_bindings.v0"
            }"#,
        )
        .expect("deserialize bindings");
        let err = resolve_protocol_cartoon_template_with_bindings(&template, &bindings)
            .expect_err("schema mismatch");
        assert!(err.contains("Unsupported protocol cartoon template bindings schema"));
    }

    #[test]
    fn render_template_with_bindings_svg_uses_bound_values() {
        let template: ProtocolCartoonTemplate = serde_json::from_str(
            r#"{
                "id": "demo.protocol",
                "events": [{ "id": "e1", "title": "T" }]
            }"#,
        )
        .expect("deserialize template");
        let bindings: ProtocolCartoonTemplateBindings = serde_json::from_str(
            r#"{
                "event_overrides": [{ "event_id": "e1", "title": "Bound Event" }]
            }"#,
        )
        .expect("deserialize bindings");
        let svg = render_protocol_cartoon_template_with_bindings_svg(&template, &bindings);
        assert!(svg.contains("Bound Event"));
    }

    #[test]
    fn parse_protocol_cartoon_aliases() {
        assert_eq!(
            ProtocolCartoonKind::parse_id("gibson.two_fragment"),
            Some(ProtocolCartoonKind::GibsonTwoFragment)
        );
        assert_eq!(
            ProtocolCartoonKind::parse_id("gibson"),
            Some(ProtocolCartoonKind::GibsonTwoFragment)
        );
        assert_eq!(
            ProtocolCartoonKind::parse_id("gibson.single_insert_dual_junction"),
            Some(ProtocolCartoonKind::GibsonSingleInsertDualJunction)
        );
        assert_eq!(
            ProtocolCartoonKind::parse_id("pcr"),
            Some(ProtocolCartoonKind::PcrAssayPair)
        );
        assert_eq!(
            ProtocolCartoonKind::parse_id("pcr.assay.pair.no_product"),
            Some(ProtocolCartoonKind::PcrAssayPairNoProduct)
        );
        assert_eq!(
            ProtocolCartoonKind::parse_id("qpcr"),
            Some(ProtocolCartoonKind::PcrAssayQpcr)
        );
        assert!(ProtocolCartoonKind::parse_id("unknown.protocol").is_none());
    }

    #[test]
    fn catalog_rows_are_deterministic() {
        let rows = protocol_cartoon_catalog_rows();
        assert_eq!(rows.len(), 5);
        assert_eq!(rows[0].id, "gibson.two_fragment");
        assert_eq!(rows[1].id, "gibson.single_insert_dual_junction");
        assert_eq!(rows[2].id, "pcr.assay.pair");
        assert_eq!(rows[3].id, "pcr.assay.pair.no_product");
        assert_eq!(rows[4].id, "pcr.assay.qpcr");
        assert!(rows[0].title.contains("Gibson"));
        assert!(rows[2].title.contains("PCR"));
        assert!(rows[4].title.contains("qPCR"));
    }

    #[test]
    fn built_in_kind_exposes_template_schema() {
        let template = protocol_cartoon_template_for_kind(&ProtocolCartoonKind::GibsonTwoFragment);
        assert_eq!(template.schema, "gentle.protocol_cartoon_template.v1");
        assert_eq!(template.id, "gibson.two_fragment");
        assert!(!template.events.is_empty());

        let pcr_template = protocol_cartoon_template_for_kind(&ProtocolCartoonKind::PcrAssayPair);
        assert_eq!(pcr_template.schema, "gentle.protocol_cartoon_template.v1");
        assert_eq!(pcr_template.id, "pcr.assay.pair");
        assert!(!pcr_template.events.is_empty());

        let qpcr_template = protocol_cartoon_template_for_kind(&ProtocolCartoonKind::PcrAssayQpcr);
        assert_eq!(qpcr_template.schema, "gentle.protocol_cartoon_template.v1");
        assert_eq!(qpcr_template.id, "pcr.assay.qpcr");
        assert!(!qpcr_template.events.is_empty());
    }

    #[test]
    fn gibson_spec_is_event_sequence() {
        let spec = protocol_cartoon_spec_for_kind(&ProtocolCartoonKind::GibsonTwoFragment);
        assert_eq!(spec.events.len(), 5);
        assert_eq!(spec.events[0].action, ProtocolCartoonAction::Context);
    }

    #[test]
    fn dual_junction_gibson_spec_shows_three_context_molecules() {
        let spec =
            protocol_cartoon_spec_for_kind(&ProtocolCartoonKind::GibsonSingleInsertDualJunction);
        assert_eq!(spec.events.len(), 5);
        assert_eq!(spec.events[0].molecules.len(), 3);
        assert_eq!(spec.events[1].molecules.len(), 3);
        assert_eq!(spec.events[2].molecules.len(), 1);
    }

    #[test]
    fn dual_junction_gibson_spec_shows_two_explicit_junctions() {
        let spec =
            protocol_cartoon_spec_for_kind(&ProtocolCartoonKind::GibsonSingleInsertDualJunction);
        let seal = &spec.events[4].molecules[0];
        assert!(
            seal.features
                .iter()
                .any(|feature| feature.id == "left_overlap")
        );
        assert!(
            seal.features
                .iter()
                .any(|feature| feature.id == "right_overlap")
        );
    }

    fn event_by_id<'a>(spec: &'a ProtocolCartoonSpec, id: &str) -> &'a ProtocolCartoonEvent {
        spec.events
            .iter()
            .find(|event| event.id == id)
            .unwrap_or_else(|| panic!("missing event '{id}'"))
    }

    fn molecule_by_id<'a>(event: &'a ProtocolCartoonEvent, id: &str) -> &'a DnaMoleculeCartoon {
        event
            .molecules
            .iter()
            .find(|molecule| molecule.id == id)
            .unwrap_or_else(|| panic!("missing molecule '{id}'"))
    }

    fn feature_by_id<'a>(molecule: &'a DnaMoleculeCartoon, id: &str) -> &'a DnaFeatureCartoon {
        molecule
            .features
            .iter()
            .find(|feature| feature.id == id)
            .unwrap_or_else(|| panic!("missing feature '{id}'"))
    }

    #[test]
    fn dual_junction_chew_back_uses_correct_5prime_exonuclease_strands() {
        let spec =
            protocol_cartoon_spec_for_kind(&ProtocolCartoonKind::GibsonSingleInsertDualJunction);
        let chew = &spec.events[1];
        let left_arm = &chew.molecules[0];
        let insert = &chew.molecules[1];
        let right_arm = &chew.molecules[2];

        assert_eq!(
            left_arm
                .features
                .last()
                .map(|f| (f.top_length_bp, f.bottom_length_bp)),
            Some((24, 0))
        );
        assert_eq!(
            insert
                .features
                .first()
                .map(|f| (f.top_length_bp, f.bottom_length_bp)),
            Some((0, 24))
        );
        assert_eq!(
            insert
                .features
                .last()
                .map(|f| (f.top_length_bp, f.bottom_length_bp)),
            Some((24, 0))
        );
        assert_eq!(
            right_arm
                .features
                .first()
                .map(|f| (f.top_length_bp, f.bottom_length_bp)),
            Some((0, 24))
        );
    }

    #[test]
    fn dual_junction_chew_back_extends_past_overlap_for_fill_in() {
        let spec =
            protocol_cartoon_spec_for_kind(&ProtocolCartoonKind::GibsonSingleInsertDualJunction);
        let chew = event_by_id(&spec, "chew_back");
        let extend = event_by_id(&spec, "extend");
        let left_arm = molecule_by_id(chew, "destination_left_chewed");
        let insert = molecule_by_id(chew, "insert_chewed");
        let right_arm = molecule_by_id(chew, "destination_right_chewed");
        let extended = molecule_by_id(extend, "extended_intermediate");

        let expected_dest_left_tail_bp = feature_by_id(extended, "dest_left_tail").length_bp;
        let expected_insert_left_tail_bp = feature_by_id(extended, "insert_left_tail").length_bp;
        let expected_insert_right_tail_bp = feature_by_id(extended, "insert_right_tail").length_bp;
        let expected_dest_right_tail_bp = feature_by_id(extended, "dest_right_tail").length_bp;

        let dest_left_tail = feature_by_id(left_arm, "dest_left_tail_exposed");
        assert_eq!(
            (
                dest_left_tail.top_length_bp,
                dest_left_tail.bottom_length_bp
            ),
            (expected_dest_left_tail_bp, 0)
        );
        let insert_left_tail = feature_by_id(insert, "insert_left_tail_exposed");
        assert_eq!(
            (
                insert_left_tail.top_length_bp,
                insert_left_tail.bottom_length_bp
            ),
            (0, expected_insert_left_tail_bp)
        );
        let insert_right_tail = feature_by_id(insert, "insert_right_tail_exposed");
        assert_eq!(
            (
                insert_right_tail.top_length_bp,
                insert_right_tail.bottom_length_bp
            ),
            (expected_insert_right_tail_bp, 0)
        );
        let dest_right_tail = feature_by_id(right_arm, "dest_right_tail_exposed");
        assert_eq!(
            (
                dest_right_tail.top_length_bp,
                dest_right_tail.bottom_length_bp
            ),
            (0, expected_dest_right_tail_bp)
        );
    }

    #[test]
    fn dual_junction_anneal_leaves_gaps_on_both_sides_of_each_overlap() {
        let spec =
            protocol_cartoon_spec_for_kind(&ProtocolCartoonKind::GibsonSingleInsertDualJunction);
        let anneal = molecule_by_id(event_by_id(&spec, "anneal"), "annealed_intermediate");
        let extend = molecule_by_id(event_by_id(&spec, "extend"), "extended_intermediate");

        let expected_dest_left_tail_bp = feature_by_id(extend, "dest_left_tail").length_bp;
        let expected_insert_left_tail_bp = feature_by_id(extend, "insert_left_tail").length_bp;
        let expected_insert_right_tail_bp = feature_by_id(extend, "insert_right_tail").length_bp;
        let expected_dest_right_tail_bp = feature_by_id(extend, "dest_right_tail").length_bp;

        let dest_left_tail = feature_by_id(anneal, "dest_left_tail");
        assert_eq!(
            (
                dest_left_tail.top_length_bp,
                dest_left_tail.bottom_length_bp
            ),
            (expected_dest_left_tail_bp, 0)
        );
        let insert_left_tail = feature_by_id(anneal, "insert_left_tail");
        assert_eq!(
            (
                insert_left_tail.top_length_bp,
                insert_left_tail.bottom_length_bp
            ),
            (0, expected_insert_left_tail_bp)
        );
        let insert_right_tail = feature_by_id(anneal, "insert_right_tail");
        assert_eq!(
            (
                insert_right_tail.top_length_bp,
                insert_right_tail.bottom_length_bp
            ),
            (expected_insert_right_tail_bp, 0)
        );
        let dest_right_tail = feature_by_id(anneal, "dest_right_tail");
        assert_eq!(
            (
                dest_right_tail.top_length_bp,
                dest_right_tail.bottom_length_bp
            ),
            (0, expected_dest_right_tail_bp)
        );
    }

    #[test]
    fn pcr_pair_spec_is_event_sequence() {
        let spec = protocol_cartoon_spec_for_kind(&ProtocolCartoonKind::PcrAssayPair);
        assert_eq!(spec.events.len(), 4);
        assert_eq!(spec.events[0].action, ProtocolCartoonAction::Context);
        assert_eq!(spec.events[1].title, "Primer Constraints");
        assert_eq!(spec.events[2].title, "Primers");
        assert_eq!(spec.events[3].title, "Product");
    }

    #[test]
    fn pcr_pair_spec_highlights_roi_and_final_amplicon() {
        let spec = protocol_cartoon_spec_for_kind(&ProtocolCartoonKind::PcrAssayPair);
        let context = &spec.events[0].molecules[0];
        let constraints = &spec.events[1].molecules[0];
        let product = &spec.events[3].molecules[0];

        assert!(context.features.iter().any(|feature| feature.id == "roi"));
        assert!(
            !context
                .features
                .iter()
                .any(|feature| feature.id == "forward_primer_site")
        );
        assert!(
            !context
                .features
                .iter()
                .any(|feature| feature.id == "reverse_primer_site")
        );
        assert!(
            constraints
                .features
                .iter()
                .any(|feature| feature.id == "roi")
        );
        assert!(
            !constraints
                .features
                .iter()
                .any(|feature| feature.id == "forward_primer_site")
        );
        assert!(
            !constraints
                .features
                .iter()
                .any(|feature| feature.id == "reverse_primer_site")
        );
        assert!(
            constraints
                .features
                .iter()
                .any(|feature| feature.id == "forward_window")
        );
        assert!(
            constraints
                .features
                .iter()
                .any(|feature| feature.id == "reverse_window")
        );
        assert!(
            product
                .features
                .iter()
                .any(|feature| feature.id == "amplicon_roi")
        );
        assert!(
            product
                .features
                .iter()
                .any(|feature| feature.id == "forward_primer_site")
        );
        assert!(
            product
                .features
                .iter()
                .any(|feature| feature.id == "reverse_primer_site")
        );
    }

    #[test]
    fn pcr_pair_no_product_finishes_with_report_only_context() {
        let spec = protocol_cartoon_spec_for_kind(&ProtocolCartoonKind::PcrAssayPairNoProduct);
        let final_event = &spec.events[3];

        assert_eq!(final_event.title, "Report Only");
        assert!(final_event.caption.contains("No amplicon product"));
        assert!(
            !final_event.molecules[0]
                .features
                .iter()
                .any(|feature| feature.id.contains("amplicon"))
        );
    }

    #[test]
    fn pcr_context_and_roi_panels_remain_primer_free() {
        let pair = protocol_cartoon_spec_for_kind(&ProtocolCartoonKind::PcrAssayPair);
        let no_product =
            protocol_cartoon_spec_for_kind(&ProtocolCartoonKind::PcrAssayPairNoProduct);
        let qpcr = protocol_cartoon_spec_for_kind(&ProtocolCartoonKind::PcrAssayQpcr);

        for spec in [&pair, &no_product, &qpcr] {
            for event in [&spec.events[0], &spec.events[1]] {
                for molecule in &event.molecules {
                    assert!(
                        !molecule
                            .features
                            .iter()
                            .any(|feature| feature.id == "forward_primer_site")
                    );
                    assert!(
                        !molecule
                            .features
                            .iter()
                            .any(|feature| feature.id == "reverse_primer_site")
                    );
                }
            }
        }
    }

    #[test]
    fn qpcr_spec_is_event_sequence_with_probe_window() {
        let spec = protocol_cartoon_spec_for_kind(&ProtocolCartoonKind::PcrAssayQpcr);
        assert_eq!(spec.events.len(), 4);
        assert_eq!(spec.events[1].title, "Primer Constraints");
        assert_eq!(spec.events[2].title, "Primers + Probe");
        assert_eq!(spec.events[3].title, "Quantify");
        assert!(
            spec.events[1].molecules[0]
                .features
                .iter()
                .any(|feature| feature.id == "probe_window")
        );
        assert!(
            !spec.events[1].molecules[0]
                .features
                .iter()
                .any(|feature| feature.id == "probe_site")
        );
    }

    #[test]
    fn qpcr_final_panel_retains_probe_window() {
        let spec = protocol_cartoon_spec_for_kind(&ProtocolCartoonKind::PcrAssayQpcr);
        let assay = &spec.events[2].molecules[0];
        let quantified = &spec.events[3].molecules[0];
        assert!(
            assay
                .features
                .iter()
                .any(|feature| feature.id == "probe_window_left")
        );
        assert!(
            assay
                .features
                .iter()
                .any(|feature| feature.id == "probe_window_right")
        );
        assert!(
            quantified
                .features
                .iter()
                .any(|feature| feature.id == "probe_window_left")
        );
        assert!(
            quantified
                .features
                .iter()
                .any(|feature| feature.id == "probe_window_right")
        );
        assert!(
            assay
                .features
                .iter()
                .all(|feature| !feature.id.contains("template"))
        );
        assert!(
            quantified
                .features
                .iter()
                .all(|feature| !feature.id.contains("template"))
        );
        assert_eq!(quantified.label, "Quantified probe-bearing amplicon");
    }

    #[test]
    fn gibson_chew_back_uses_3prime_overhangs() {
        let spec = protocol_cartoon_spec_for_kind(&ProtocolCartoonKind::GibsonTwoFragment);
        assert_eq!(
            spec.events[1].molecules[0].right_end,
            Some(DnaEndStyle::Sticky {
                polarity: OverhangPolarity::ThreePrime,
                nt: 30
            })
        );
        assert_eq!(
            spec.events[1].molecules[1].left_end,
            Some(DnaEndStyle::Sticky {
                polarity: OverhangPolarity::ThreePrime,
                nt: 30
            })
        );
    }

    #[test]
    fn gibson_anneal_retains_single_stranded_unique_tails() {
        let spec = protocol_cartoon_spec_for_kind(&ProtocolCartoonKind::GibsonTwoFragment);
        let anneal = &spec.events[2].molecules[0];
        assert_eq!(anneal.features[1].top_length_bp, 12);
        assert_eq!(anneal.features[1].bottom_length_bp, 0);
        assert_eq!(anneal.features[2].top_length_bp, 30);
        assert_eq!(anneal.features[2].bottom_length_bp, 30);
        assert_eq!(anneal.features[3].top_length_bp, 0);
        assert_eq!(anneal.features[3].bottom_length_bp, 12);
    }

    #[test]
    fn gibson_extend_leaves_two_strand_specific_nicks() {
        let spec = protocol_cartoon_spec_for_kind(&ProtocolCartoonKind::GibsonTwoFragment);
        let extend = &spec.events[3].molecules[0];
        assert!(extend.features[0].bottom_nick_after);
        assert!(!extend.features[0].top_nick_after);
        assert!(extend.features[3].top_nick_after);
        assert!(!extend.features[3].bottom_nick_after);
    }

    #[test]
    fn gibson_overlap_keeps_stable_display_total_across_panels() {
        let spec = protocol_cartoon_spec_for_kind(&ProtocolCartoonKind::GibsonTwoFragment);
        let display_total = |molecule: &DnaMoleculeCartoon| -> usize {
            molecule
                .features
                .iter()
                .map(|feature| feature.length_bp)
                .sum()
        };

        assert_eq!(display_total(&spec.events[0].molecules[0]), 150);
        assert_eq!(display_total(&spec.events[0].molecules[1]), 150);
        assert_eq!(display_total(&spec.events[1].molecules[0]), 150);
        assert_eq!(display_total(&spec.events[1].molecules[1]), 150);
        assert_eq!(display_total(&spec.events[2].molecules[0]), 150);
        assert_eq!(display_total(&spec.events[3].molecules[0]), 150);
        assert_eq!(display_total(&spec.events[4].molecules[0]), 150);
    }

    #[test]
    fn validate_rejects_empty_protocol() {
        let invalid = ProtocolCartoonSpec {
            id: "".to_string(),
            title: "x".to_string(),
            summary: "x".to_string(),
            events: vec![],
        };
        assert!(invalid.validate().is_err());
    }

    #[test]
    fn validate_rejects_circular_molecule_with_end_styles() {
        let mut spec = custom_test_spec();
        spec.events[0].molecules[1].left_end = Some(DnaEndStyle::Blunt);
        let err = spec.validate().expect_err("spec should be invalid");
        assert!(err.contains("Circular molecule"));
    }

    #[test]
    fn validate_rejects_zero_length_sticky_end() {
        let mut spec = custom_test_spec();
        spec.events[0].molecules[0].left_end = Some(DnaEndStyle::Sticky {
            polarity: OverhangPolarity::ThreePrime,
            nt: 0,
        });
        let err = spec.validate().expect_err("spec should be invalid");
        assert!(err.contains("zero nt overhang"));
    }

    #[test]
    fn render_custom_spec_covers_linear_and_circular() {
        let svg = render_protocol_cartoon_spec_svg(&custom_test_spec());
        assert!(svg.contains("data-topology=\"linear\""));
        assert!(svg.contains("data-topology=\"circular\""));
        assert!(svg.contains("5&apos; 6nt"));
        assert!(svg.contains("blunt"));
    }

    #[test]
    fn render_gibson_svg_contains_expected_labels() {
        let svg = render_protocol_cartoon_svg(&ProtocolCartoonKind::GibsonTwoFragment);
        assert!(svg.contains("<svg"));
        assert!(svg.contains("width=\""));
        assert!(svg.contains("height=\""));
        assert!(svg.contains("Context"));
        assert!(svg.contains("Chew-back"));
        assert!(svg.contains("Extend"));
        assert!(svg.contains("Seal"));
        assert!(svg.contains("30 bp"));
    }

    #[test]
    fn render_dual_junction_gibson_svg_contains_expected_labels() {
        let svg = render_protocol_cartoon_svg(&ProtocolCartoonKind::GibsonSingleInsertDualJunction);
        assert!(svg.contains("<svg"));
        assert!(svg.contains("Destination left arm"));
        assert!(svg.contains("Destination right arm"));
        assert!(svg.contains("Annealed intermediate"));
        assert!(svg.contains("Sealed product"));
    }

    #[test]
    fn render_pcr_pair_svg_contains_expected_labels() {
        let svg = render_protocol_cartoon_svg(&ProtocolCartoonKind::PcrAssayPair);
        assert!(svg.contains("<svg"));
        assert!(svg.contains("PCR Assay"));
        assert!(svg.contains("Context + ROI"));
        assert!(svg.contains("Primer Constraints"));
        assert!(svg.contains("Green ="));
        assert!(svg.contains("Red = Upstream primer window for"));
        assert!(svg.contains("Blue = Downstream primer window for"));
        assert!(svg.contains("ROI-bounded amplicon"));
        assert!(svg.contains("data-primer-glyph=\"F\""));
        assert!(svg.contains("data-primer-glyph=\"R\""));
        assert!(svg.contains("data-primer-anchor=\"bottom\""));
        assert!(!svg.contains("data-primer-rotation=\"180\""));
        assert!(svg.contains("5&apos;"));
        assert!(svg.contains("3&apos;"));
        assert!(!svg.contains("data-feature-marker=\"F\""));
        assert!(!svg.contains("data-feature-marker=\"R\""));
    }

    #[test]
    fn render_pcr_report_only_svg_contains_expected_labels() {
        let svg = render_protocol_cartoon_svg(&ProtocolCartoonKind::PcrAssayPairNoProduct);
        assert!(svg.contains("<svg"));
        assert!(svg.contains("Report Only"));
        assert!(svg.contains("Template + ROI retained"));
    }

    #[test]
    fn render_qpcr_svg_contains_expected_labels() {
        let svg = render_protocol_cartoon_svg(&ProtocolCartoonKind::PcrAssayQpcr);
        assert!(svg.contains("<svg"));
        assert!(svg.contains("qPCR Assay"));
        assert!(svg.contains("Primer Constraints"));
        assert!(svg.contains("Primers + Probe"));
        assert!(svg.contains("Slate = Internal probe window"));
        assert!(svg.contains("Quantify"));
        assert!(svg.contains("data-primer-glyph=\"F\""));
        assert!(svg.contains("data-primer-glyph=\"R\""));
        assert!(svg.contains("data-primer-glyph=\"P\""));
        assert!(svg.contains("data-primer-anchor=\"bottom\""));
        assert!(!svg.contains("data-primer-rotation=\"180\""));
        assert!(!svg.contains("data-feature-marker=\"P\""));
        assert!(!svg.contains("data-feature-marker=\"F\""));
        assert!(!svg.contains("data-feature-marker=\"R\""));
    }

    #[test]
    fn render_not_shown_end_omits_end_label() {
        let mut spec = custom_test_spec();
        spec.events[0].molecules[0].left_end = Some(DnaEndStyle::NotShown);
        let svg = render_protocol_cartoon_spec_svg(&spec);
        assert!(!svg.contains("5&apos; 6nt"));
    }

    #[test]
    fn render_end_annotations_use_compact_mid_strand_layout() {
        let svg = render_protocol_cartoon_spec_svg(&custom_test_spec());
        assert!(svg.contains(".pc_end_label { font: 600 10px"));
        assert!(svg.contains("dominant-baseline=\"middle\""));
    }

    #[test]
    fn render_continuation_end_draws_dots_without_label() {
        let spec = ProtocolCartoonSpec {
            id: "continuation.test".to_string(),
            title: "Continuation".to_string(),
            summary: "continuation markers".to_string(),
            events: vec![ProtocolCartoonEvent {
                id: "evt".to_string(),
                title: "Continuation".to_string(),
                caption: "DNA continues beyond the panel".to_string(),
                action: ProtocolCartoonAction::Context,
                molecules: vec![DnaMoleculeCartoon {
                    id: "mol".to_string(),
                    label: "Fragment".to_string(),
                    topology: DnaTopologyCartoon::Linear,
                    features: vec![DnaFeatureCartoon {
                        id: "body".to_string(),
                        label: "Body".to_string(),
                        length_bp: 60,
                        top_length_bp: 60,
                        bottom_length_bp: 60,
                        color_hex: "#55a84f".to_string(),
                        bottom_color_hex: "#55a84f".to_string(),
                        top_nick_after: false,
                        bottom_nick_after: false,
                    }],
                    left_end: Some(DnaEndStyle::Continuation),
                    right_end: Some(DnaEndStyle::NotShown),
                }],
            }],
        };
        let svg = render_protocol_cartoon_spec_svg(&spec);
        assert_eq!(svg.matches("fill=\"#88a3ad\"").count(), 3);
        assert!(!svg.contains(">blunt<"));
        assert!(!svg.contains("5&apos;"));
        assert!(!svg.contains("3&apos;"));
    }

    #[test]
    fn render_split_overlap_uses_distinct_top_and_bottom_colors() {
        let spec = ProtocolCartoonSpec {
            id: "split.overlap".to_string(),
            title: "Split Overlap".to_string(),
            summary: "top/bottom strand colors differ in overlap".to_string(),
            events: vec![ProtocolCartoonEvent {
                id: "anneal".to_string(),
                title: "Anneal".to_string(),
                caption: "Hybridized overlap".to_string(),
                action: ProtocolCartoonAction::Anneal,
                molecules: vec![DnaMoleculeCartoon {
                    id: "hybrid".to_string(),
                    label: "Hybrid".to_string(),
                    topology: DnaTopologyCartoon::Linear,
                    features: vec![hybrid_feature_block(
                        "overlap", "overlap", 30, "#d7ae2f", "#66aefb",
                    )],
                    left_end: Some(DnaEndStyle::NotShown),
                    right_end: Some(DnaEndStyle::NotShown),
                }],
            }],
        };
        let svg = render_protocol_cartoon_spec_svg(&spec);
        assert!(svg.contains("#d7ae2f"));
        assert!(svg.contains("#66aefb"));
    }

    #[test]
    fn render_merges_adjacent_same_color_spans_on_one_strand() {
        let spec = ProtocolCartoonSpec {
            id: "merged.spans".to_string(),
            title: "Merged Spans".to_string(),
            summary: "adjacent same-color spans should render continuously".to_string(),
            events: vec![ProtocolCartoonEvent {
                id: "evt".to_string(),
                title: "Merged".to_string(),
                caption: "same-color adjacency".to_string(),
                action: ProtocolCartoonAction::Context,
                molecules: vec![DnaMoleculeCartoon {
                    id: "mol".to_string(),
                    label: "Merged".to_string(),
                    topology: DnaTopologyCartoon::Linear,
                    features: vec![
                        DnaFeatureCartoon {
                            id: "body_ds".to_string(),
                            label: "duplex".to_string(),
                            length_bp: 10,
                            top_length_bp: 10,
                            bottom_length_bp: 10,
                            color_hex: "#111111".to_string(),
                            bottom_color_hex: "#111111".to_string(),
                            top_nick_after: false,
                            bottom_nick_after: false,
                        },
                        DnaFeatureCartoon {
                            id: "body_ss".to_string(),
                            label: "top tail".to_string(),
                            length_bp: 5,
                            top_length_bp: 5,
                            bottom_length_bp: 0,
                            color_hex: "#111111".to_string(),
                            bottom_color_hex: "#111111".to_string(),
                            top_nick_after: false,
                            bottom_nick_after: false,
                        },
                        DnaFeatureCartoon {
                            id: "other".to_string(),
                            label: "other".to_string(),
                            length_bp: 10,
                            top_length_bp: 10,
                            bottom_length_bp: 10,
                            color_hex: "#222222".to_string(),
                            bottom_color_hex: "#222222".to_string(),
                            top_nick_after: false,
                            bottom_nick_after: false,
                        },
                    ],
                    left_end: Some(DnaEndStyle::NotShown),
                    right_end: Some(DnaEndStyle::NotShown),
                }],
            }],
        };
        let svg = render_protocol_cartoon_spec_svg(&spec);
        assert_eq!(svg.matches("fill=\"#111111\"").count(), 2);
    }

    #[test]
    fn normalize_hex_color_rejects_invalid_values() {
        assert_eq!(normalize_hex_color("#aabbcc"), "#aabbcc");
        assert_eq!(normalize_hex_color("navy"), "#8ea7b1");
    }

    #[test]
    fn invalid_spec_renders_error_svg() {
        let invalid = ProtocolCartoonSpec {
            id: "broken".to_string(),
            title: "Broken".to_string(),
            summary: "Broken".to_string(),
            events: vec![],
        };
        let svg = render_protocol_cartoon_spec_svg(&invalid);
        assert!(svg.contains("Invalid protocol cartoon"));
        assert!(svg.contains("broken"));
        assert!(svg.contains("width=\"900\""));
        assert!(svg.contains("height=\"220\""));
    }
}
