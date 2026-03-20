//! Protocol-cartoon catalog and deterministic SVG rendering.
//!
//! This module defines a reusable abstraction for protocol cartoons:
//! DNA molecules (linear or circular) are composed of colored feature fragments,
//! optional sticky/blunt end styles, and protocol events arranged as a sequence
//! of snapshots.

use serde::{Deserialize, Serialize};

/// Stable protocol-cartoon identifiers exposed through engine operations.
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq, Hash)]
pub enum ProtocolCartoonKind {
    #[serde(rename = "gibson.two_fragment")]
    GibsonTwoFragment,
}

impl ProtocolCartoonKind {
    /// Canonical string identifier used in shell/CLI input and output rows.
    pub fn id(&self) -> &'static str {
        match self {
            Self::GibsonTwoFragment => "gibson.two_fragment",
        }
    }

    /// Human-readable title for UI/help surfaces.
    pub fn title(&self) -> &'static str {
        match self {
            Self::GibsonTwoFragment => "Gibson Assembly (two-fragment conceptual)",
        }
    }

    /// Short semantics summary used in list outputs.
    pub fn summary(&self) -> &'static str {
        match self {
            Self::GibsonTwoFragment => {
                "Event-sequence cartoon with sticky/blunt ends and linear/circular molecule glyphs"
            }
        }
    }

    /// Parse supported shell/CLI aliases into the canonical enum.
    pub fn parse_id(raw: &str) -> Option<Self> {
        match raw.trim().to_ascii_lowercase().as_str() {
            "gibson.two_fragment" | "gibson.two-fragment" | "gibson_two_fragment" | "gibson" => {
                Some(Self::GibsonTwoFragment)
            }
            _ => None,
        }
    }

    /// Deterministic ordered catalog for list commands.
    pub fn catalog() -> Vec<Self> {
        vec![Self::GibsonTwoFragment]
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
    pub color_hex: Option<String>,
}

/// One feature fragment shown as a colored span in a molecule cartoon.
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq)]
pub struct DnaFeatureCartoon {
    pub id: String,
    pub label: String,
    pub length_bp: usize,
    pub color_hex: String,
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
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq)]
pub enum DnaEndStyle {
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
                    if feature.length_bp == 0 {
                        return Err(format!(
                            "Feature '{}' in molecule '{}' has zero length",
                            feature.id, molecule.id
                        ));
                    }
                }
            }
        }
        Ok(())
    }
}

fn protocol_cartoon_template_schema_id() -> String {
    "gentle.protocol_cartoon_template.v1".to_string()
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
                let color_hex = feature.color_hex.clone().unwrap_or_else(|| {
                    defaults.palette[feature_idx % defaults.palette.len()].clone()
                });
                features.push(DnaFeatureCartoon {
                    id: feature_id,
                    label: feature_label,
                    length_bp,
                    color_hex,
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
                    color_hex: "#8ea7b1".to_string(),
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
    let panel_h = 510.0_f32;
    let panel_gap = 24.0_f32;
    let margin_x = 36.0_f32;
    let margin_y = 28.0_f32;
    let header_h = 92.0_f32;
    let footer_h = 26.0_f32;

    let event_count = spec.events.len() as f32;
    let body_w = event_count * panel_w + (event_count - 1.0).max(0.0) * panel_gap;
    let canvas_w = margin_x * 2.0 + body_w;
    let canvas_h = margin_y * 2.0 + header_h + panel_h + footer_h;

    let mut svg = String::new();
    svg.push_str(&format!(
        "<svg xmlns=\"http://www.w3.org/2000/svg\" viewBox=\"0 0 {:.0} {:.0}\" role=\"img\" aria-labelledby=\"title desc\" data-protocol-id=\"{}\">",
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
        ".pc_end_label { font: 600 11px 'Avenir Next', 'Trebuchet MS', 'Segoe UI', sans-serif; fill: #375965; }",
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

        let mut molecule_y = panel_top + 132.0;
        for molecule in &event.molecules {
            render_molecule(
                &mut svg,
                molecule,
                panel_x + 16.0,
                molecule_y,
                panel_w - 32.0,
            );
            molecule_y += 84.0;
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
            render_linear_molecule(svg, molecule, x, y + 14.0, width);
        }
        DnaTopologyCartoon::Circular => {
            render_circular_molecule(svg, molecule, x, y + 14.0, width);
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
    let total_bp: usize = molecule
        .features
        .iter()
        .map(|f| f.length_bp)
        .sum::<usize>()
        .max(1);
    let body_w = (width - 8.0).max(80.0);
    let body_x = x + 4.0;

    let mut left_top = body_x;
    let mut left_bottom = body_x;
    let mut right_top = body_x + body_w;
    let mut right_bottom = body_x + body_w;

    if let Some(end) = molecule.left_end.as_ref() {
        apply_left_end_style(&mut left_top, &mut left_bottom, end);
        render_end_annotation(svg, x, y + 34.0, true, end);
    }
    if let Some(end) = molecule.right_end.as_ref() {
        apply_right_end_style(&mut right_top, &mut right_bottom, end);
        render_end_annotation(svg, x + width, y + 34.0, false, end);
    }

    svg.push_str(&format!(
        "<rect x=\"{:.1}\" y=\"{:.1}\" width=\"{:.1}\" height=\"8\" rx=\"4\" fill=\"#c8d9de\"/>",
        left_top,
        y,
        (right_top - left_top).max(2.0)
    ));
    svg.push_str(&format!(
        "<rect x=\"{:.1}\" y=\"{:.1}\" width=\"{:.1}\" height=\"8\" rx=\"4\" fill=\"#c8d9de\"/>",
        left_bottom,
        y + 11.0,
        (right_bottom - left_bottom).max(2.0)
    ));

    let mut cursor = body_x;
    for (feature_idx, feature) in molecule.features.iter().enumerate() {
        let frac = feature.length_bp as f32 / total_bp as f32;
        let seg_w = if feature_idx + 1 == molecule.features.len() {
            (body_x + body_w - cursor).max(0.5)
        } else {
            (body_w * frac).max(0.5)
        };
        let color = normalize_hex_color(&feature.color_hex);
        svg.push_str(&format!(
            "<rect x=\"{:.1}\" y=\"{:.1}\" width=\"{:.1}\" height=\"8\" fill=\"{}\"/>",
            cursor, y, seg_w, color
        ));
        svg.push_str(&format!(
            "<rect x=\"{:.1}\" y=\"{:.1}\" width=\"{:.1}\" height=\"8\" fill=\"{}\"/>",
            cursor,
            y + 11.0,
            seg_w,
            color
        ));
        cursor += seg_w;
    }

    svg.push_str(&format!(
        "<line x1=\"{:.1}\" y1=\"{:.1}\" x2=\"{:.1}\" y2=\"{:.1}\" stroke=\"#86a5af\" stroke-width=\"1\"/>",
        body_x,
        y + 24.0,
        body_x + body_w,
        y + 24.0
    ));
}

fn render_circular_molecule(
    svg: &mut String,
    molecule: &DnaMoleculeCartoon,
    x: f32,
    y: f32,
    width: f32,
) {
    let total_bp: usize = molecule
        .features
        .iter()
        .map(|f| f.length_bp)
        .sum::<usize>()
        .max(1);
    let cx = x + width * 0.5;
    let cy = y + 16.0;
    let r = (width * 0.24).clamp(16.0, 30.0);

    svg.push_str(&format!(
        "<circle cx=\"{:.1}\" cy=\"{:.1}\" r=\"{:.1}\" fill=\"none\" stroke=\"#c8d9de\" stroke-width=\"10\"/>",
        cx, cy, r
    ));

    let mut start = -std::f32::consts::FRAC_PI_2;
    for feature in &molecule.features {
        let span = std::f32::consts::TAU * (feature.length_bp as f32 / total_bp as f32);
        let end = start + span;
        let color = normalize_hex_color(&feature.color_hex);
        let path = arc_path(cx, cy, r, start, end);
        svg.push_str(&format!(
            "<path d=\"{}\" fill=\"none\" stroke=\"{}\" stroke-width=\"10\" stroke-linecap=\"butt\"/>",
            path, color
        ));
        start = end;
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

fn render_end_annotation(svg: &mut String, x: f32, y: f32, is_left: bool, end: &DnaEndStyle) {
    let text = match end {
        DnaEndStyle::Blunt => "blunt".to_string(),
        DnaEndStyle::Sticky { polarity, nt } => format!("{} {}nt", polarity.label(), nt),
    };
    let anchor = if is_left { "start" } else { "end" };
    let dx = if is_left { 0.0 } else { -2.0 };
    svg.push_str(&format!(
        "<text x=\"{:.1}\" y=\"{:.1}\" class=\"pc_end_label\" text-anchor=\"{}\">{}</text>",
        x + dx,
        y,
        anchor,
        escape_xml(&text)
    ));
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
    let words = text.split_whitespace().collect::<Vec<_>>();
    if words.is_empty() {
        return;
    }

    let mut lines: Vec<String> = vec![];
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

fn render_invalid_protocol_svg(id: &str, message: &str) -> String {
    format!(
        "<svg xmlns=\"http://www.w3.org/2000/svg\" viewBox=\"0 0 900 220\"><rect x=\"0\" y=\"0\" width=\"900\" height=\"220\" fill=\"#fff7ed\"/><text x=\"24\" y=\"56\" font-family=\"Trebuchet MS, Segoe UI, sans-serif\" font-size=\"30\" font-weight=\"700\" fill=\"#8a2f1f\">Invalid protocol cartoon: {}</text><text x=\"24\" y=\"96\" font-family=\"Trebuchet MS, Segoe UI, sans-serif\" font-size=\"18\" fill=\"#7a4a3b\">{}</text></svg>",
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
                                color_hex: Some(feature.color_hex),
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

fn gibson_two_fragment_spec() -> ProtocolCartoonSpec {
    ProtocolCartoonSpec {
        id: "gibson.two_fragment".to_string(),
        title: "GENtle Protocol Cartoon: Gibson Two-Fragment Assembly".to_string(),
        summary:
            "Feature-colored DNA cartoons with explicit sticky/blunt ends across event sequence"
                .to_string(),
        events: vec![
            ProtocolCartoonEvent {
                id: "context".to_string(),
                title: "Reaction Context".to_string(),
                caption:
                    "Two linear dsDNA molecules with designed overlap enter one isothermal mix."
                        .to_string(),
                action: ProtocolCartoonAction::Context,
                molecules: vec![
                    DnaMoleculeCartoon {
                        id: "fragment_a".to_string(),
                        label: "Fragment A".to_string(),
                        topology: DnaTopologyCartoon::Linear,
                        features: vec![
                            DnaFeatureCartoon {
                                id: "a_backbone".to_string(),
                                label: "Backbone".to_string(),
                                length_bp: 980,
                                color_hex: "#0f5d75".to_string(),
                            },
                            DnaFeatureCartoon {
                                id: "a_overlap".to_string(),
                                label: "Overlap".to_string(),
                                length_bp: 26,
                                color_hex: "#1aa3a1".to_string(),
                            },
                        ],
                        left_end: Some(DnaEndStyle::Blunt),
                        right_end: Some(DnaEndStyle::Blunt),
                    },
                    DnaMoleculeCartoon {
                        id: "fragment_b".to_string(),
                        label: "Fragment B".to_string(),
                        topology: DnaTopologyCartoon::Linear,
                        features: vec![
                            DnaFeatureCartoon {
                                id: "b_overlap".to_string(),
                                label: "Overlap".to_string(),
                                length_bp: 26,
                                color_hex: "#1aa3a1".to_string(),
                            },
                            DnaFeatureCartoon {
                                id: "b_insert".to_string(),
                                label: "Insert".to_string(),
                                length_bp: 760,
                                color_hex: "#e3b230".to_string(),
                            },
                        ],
                        left_end: Some(DnaEndStyle::Blunt),
                        right_end: Some(DnaEndStyle::Blunt),
                    },
                ],
            },
            ProtocolCartoonEvent {
                id: "chew_back".to_string(),
                title: "Chew-back".to_string(),
                caption:
                    "Exonuclease exposes complementary overhangs; sticky ends become explicit."
                        .to_string(),
                action: ProtocolCartoonAction::Cut,
                molecules: vec![
                    DnaMoleculeCartoon {
                        id: "fragment_a_chewed".to_string(),
                        label: "A after exonuclease".to_string(),
                        topology: DnaTopologyCartoon::Linear,
                        features: vec![
                            DnaFeatureCartoon {
                                id: "a_backbone".to_string(),
                                label: "Backbone".to_string(),
                                length_bp: 980,
                                color_hex: "#0f5d75".to_string(),
                            },
                            DnaFeatureCartoon {
                                id: "a_overlap".to_string(),
                                label: "Overlap".to_string(),
                                length_bp: 26,
                                color_hex: "#1aa3a1".to_string(),
                            },
                        ],
                        left_end: Some(DnaEndStyle::Blunt),
                        right_end: Some(DnaEndStyle::Sticky {
                            polarity: OverhangPolarity::ThreePrime,
                            nt: 20,
                        }),
                    },
                    DnaMoleculeCartoon {
                        id: "fragment_b_chewed".to_string(),
                        label: "B after exonuclease".to_string(),
                        topology: DnaTopologyCartoon::Linear,
                        features: vec![
                            DnaFeatureCartoon {
                                id: "b_overlap".to_string(),
                                label: "Overlap".to_string(),
                                length_bp: 26,
                                color_hex: "#1aa3a1".to_string(),
                            },
                            DnaFeatureCartoon {
                                id: "b_insert".to_string(),
                                label: "Insert".to_string(),
                                length_bp: 760,
                                color_hex: "#e3b230".to_string(),
                            },
                        ],
                        left_end: Some(DnaEndStyle::Sticky {
                            polarity: OverhangPolarity::FivePrime,
                            nt: 20,
                        }),
                        right_end: Some(DnaEndStyle::Blunt),
                    },
                ],
            },
            ProtocolCartoonEvent {
                id: "anneal".to_string(),
                title: "Anneal".to_string(),
                caption:
                    "Complementary overhangs pair in register; fragments align before covalent sealing."
                        .to_string(),
                action: ProtocolCartoonAction::Anneal,
                molecules: vec![DnaMoleculeCartoon {
                    id: "annealed_intermediate".to_string(),
                    label: "Annealed intermediate".to_string(),
                    topology: DnaTopologyCartoon::Linear,
                    features: vec![
                        DnaFeatureCartoon {
                            id: "a_backbone".to_string(),
                            label: "Backbone".to_string(),
                            length_bp: 980,
                            color_hex: "#0f5d75".to_string(),
                        },
                        DnaFeatureCartoon {
                            id: "overlap".to_string(),
                            label: "Overlap duplex".to_string(),
                            length_bp: 26,
                            color_hex: "#1aa3a1".to_string(),
                        },
                        DnaFeatureCartoon {
                            id: "b_insert".to_string(),
                            label: "Insert".to_string(),
                            length_bp: 760,
                            color_hex: "#e3b230".to_string(),
                        },
                    ],
                    left_end: Some(DnaEndStyle::Blunt),
                    right_end: Some(DnaEndStyle::Blunt),
                }],
            },
            ProtocolCartoonEvent {
                id: "extend_and_seal".to_string(),
                title: "Extend + Seal".to_string(),
                caption:
                    "Polymerase fills gaps and ligase restores covalent continuity across both strands."
                        .to_string(),
                action: ProtocolCartoonAction::Seal,
                molecules: vec![
                    DnaMoleculeCartoon {
                        id: "assembled_linear".to_string(),
                        label: "Assembled linear product".to_string(),
                        topology: DnaTopologyCartoon::Linear,
                        features: vec![
                            DnaFeatureCartoon {
                                id: "a_backbone".to_string(),
                                label: "Backbone".to_string(),
                                length_bp: 980,
                                color_hex: "#0f5d75".to_string(),
                            },
                            DnaFeatureCartoon {
                                id: "overlap".to_string(),
                                label: "Joined overlap".to_string(),
                                length_bp: 26,
                                color_hex: "#1aa3a1".to_string(),
                            },
                            DnaFeatureCartoon {
                                id: "b_insert".to_string(),
                                label: "Insert".to_string(),
                                length_bp: 760,
                                color_hex: "#e3b230".to_string(),
                            },
                        ],
                        left_end: Some(DnaEndStyle::Blunt),
                        right_end: Some(DnaEndStyle::Blunt),
                    },
                    DnaMoleculeCartoon {
                        id: "assembled_circular".to_string(),
                        label: "Circularized representation".to_string(),
                        topology: DnaTopologyCartoon::Circular,
                        features: vec![
                            DnaFeatureCartoon {
                                id: "a_backbone".to_string(),
                                label: "Backbone".to_string(),
                                length_bp: 980,
                                color_hex: "#0f5d75".to_string(),
                            },
                            DnaFeatureCartoon {
                                id: "overlap".to_string(),
                                label: "Junction".to_string(),
                                length_bp: 26,
                                color_hex: "#1aa3a1".to_string(),
                            },
                            DnaFeatureCartoon {
                                id: "b_insert".to_string(),
                                label: "Insert".to_string(),
                                length_bp: 760,
                                color_hex: "#e3b230".to_string(),
                            },
                        ],
                        left_end: None,
                        right_end: None,
                    },
                ],
            },
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
                                color_hex: "#004488".to_string(),
                            },
                            DnaFeatureCartoon {
                                id: "f2".to_string(),
                                label: "f2".to_string(),
                                length_bp: 80,
                                color_hex: "#ee9933".to_string(),
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
                            color_hex: "#118833".to_string(),
                        }],
                        left_end: None,
                        right_end: None,
                    },
                ],
            }],
        }
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
    fn parse_protocol_cartoon_aliases() {
        assert_eq!(
            ProtocolCartoonKind::parse_id("gibson.two_fragment"),
            Some(ProtocolCartoonKind::GibsonTwoFragment)
        );
        assert_eq!(
            ProtocolCartoonKind::parse_id("gibson"),
            Some(ProtocolCartoonKind::GibsonTwoFragment)
        );
        assert!(ProtocolCartoonKind::parse_id("unknown.protocol").is_none());
    }

    #[test]
    fn catalog_rows_are_deterministic() {
        let rows = protocol_cartoon_catalog_rows();
        assert_eq!(rows.len(), 1);
        assert_eq!(rows[0].id, "gibson.two_fragment");
        assert!(rows[0].title.contains("Gibson"));
    }

    #[test]
    fn built_in_kind_exposes_template_schema() {
        let template = protocol_cartoon_template_for_kind(&ProtocolCartoonKind::GibsonTwoFragment);
        assert_eq!(template.schema, "gentle.protocol_cartoon_template.v1");
        assert_eq!(template.id, "gibson.two_fragment");
        assert!(!template.events.is_empty());
    }

    #[test]
    fn gibson_spec_is_event_sequence() {
        let spec = protocol_cartoon_spec_for_kind(&ProtocolCartoonKind::GibsonTwoFragment);
        assert!(spec.events.len() >= 4);
        assert_eq!(spec.events[0].action, ProtocolCartoonAction::Context);
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
        assert!(svg.contains("Reaction Context"));
        assert!(svg.contains("Chew-back"));
        assert!(svg.contains("Extend + Seal"));
        assert!(svg.contains("Circularized representation"));
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
    }
}
