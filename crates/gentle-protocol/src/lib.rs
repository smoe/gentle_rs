//! Shared machine-readable GENtle contracts.
//!
//! This crate is intentionally scaffolded first in the workspace split so
//! stable cross-adapter payloads can move here before execution or GUI code.
//! The first extracted slice is intentionally small: stable identifier aliases,
//! shared analysis enums, and the portable engine error payload.

pub mod construct_reasoning;
pub mod dna_ladder;

use serde::{Deserialize, Serialize};
use std::{
    collections::{BTreeMap, HashMap},
    error::Error,
    fmt,
};

pub use construct_reasoning::{
    CONSTRUCT_CANDIDATE_SCHEMA, CONSTRUCT_OBJECTIVE_SCHEMA, CONSTRUCT_REASONING_GRAPH_SCHEMA,
    CONSTRUCT_REASONING_STORE_SCHEMA, ConstructCandidate, ConstructObjective,
    ConstructReasoningGraph, ConstructReasoningStore, ConstructRole, DESIGN_DECISION_NODE_SCHEMA,
    DESIGN_EVIDENCE_SCHEMA, DESIGN_FACT_SCHEMA, DecisionMethod, DesignDecisionNode, DesignEvidence,
    DesignFact, EditableStatus, EvidenceClass, EvidenceScope, HelperConstructProfile,
    HostLifecycleRole, HostProfileRecord, HostRouteStep,
};
pub use dna_ladder::{
    DNALadder, DNALadderBand, DNALadders, Ladder, LadderBand, LadderCatalog, LadderMolecule,
    RNALadder, RNALadderBand, RNALadders, default_dna_ladders, default_rna_ladders,
};

/// Stable identifier for one sequence entry stored in project state.
pub type SeqId = String;
/// Stable identifier for one executed operation journal row.
pub type OpId = String;
/// Caller-supplied identifier that groups operations into one workflow/run.
pub type RunId = String;
/// Stable identifier for one lineage graph node.
pub type NodeId = String;
/// Stable identifier for one wet-lab-style container record.
pub type ContainerId = String;

#[derive(Debug, Clone, Serialize, Deserialize)]
/// High-level provenance class for how a sequence entered or was derived in
/// the project graph.
pub enum SequenceOrigin {
    ImportedGenomic,
    ImportedCdna,
    ImportedSynthetic,
    ImportedUnknown,
    Derived,
    InSilicoSelection,
    Branch,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
/// One sequence node in the persisted lineage DAG.
pub struct LineageNode {
    pub node_id: NodeId,
    pub seq_id: SeqId,
    pub created_by_op: Option<OpId>,
    pub origin: SequenceOrigin,
    pub created_at_unix_ms: u128,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
/// Directed lineage edge linking one parent node to one derived node.
pub struct LineageEdge {
    pub from_node_id: NodeId,
    pub to_node_id: NodeId,
    pub op_id: OpId,
    pub run_id: RunId,
}

#[derive(Debug, Clone, Copy, Serialize, Deserialize, PartialEq, Eq, Default)]
#[serde(rename_all = "snake_case")]
/// Terminal status for one recorded macro instance expansion.
pub enum MacroInstanceStatus {
    #[default]
    Ok,
    Failed,
    Cancelled,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
/// Captured binding payload for one typed macro input/output port.
pub struct LineageMacroPortBinding {
    pub port_id: String,
    pub kind: String,
    pub required: bool,
    pub cardinality: String,
    pub values: Vec<String>,
    pub description: Option<String>,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
/// Persistent lineage record for one macro execution attempt.
///
/// Stored in project state for audit/debug and graph visualization across GUI,
/// CLI, and agent adapters.
pub struct LineageMacroInstance {
    pub macro_instance_id: String,
    pub routine_id: Option<String>,
    pub routine_title: Option<String>,
    pub template_name: Option<String>,
    pub run_id: String,
    pub created_at_unix_ms: u128,
    pub bound_inputs: Vec<LineageMacroPortBinding>,
    pub bound_outputs: Vec<LineageMacroPortBinding>,
    pub expanded_op_ids: Vec<String>,
    pub status: MacroInstanceStatus,
    pub status_message: Option<String>,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
/// Persisted lineage DAG plus macro-instance audit trail.
pub struct LineageGraph {
    pub nodes: HashMap<NodeId, LineageNode>,
    pub seq_to_node: HashMap<SeqId, NodeId>,
    pub edges: Vec<LineageEdge>,
    pub macro_instances: Vec<LineageMacroInstance>,
    pub next_node_counter: u64,
    pub next_macro_instance_counter: u64,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
/// Container semantic class used by container-aware operations.
pub enum ContainerKind {
    Singleton,
    Pool,
    Selection,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
/// Group of sequence ids participating in container-aware workflows.
pub struct Container {
    pub container_id: ContainerId,
    pub kind: ContainerKind,
    pub name: Option<String>,
    pub members: Vec<SeqId>,
    #[serde(default = "default_true")]
    pub declared_contents_exclusive: bool,
    pub created_by_op: Option<OpId>,
    pub created_at_unix_ms: u128,
}

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq)]
#[serde(rename_all = "snake_case")]
/// Arrangement layout mode for gel-style or plate-style workflows.
pub enum ArrangementMode {
    Serial,
    Plate,
}

#[derive(Debug, Clone, Copy, Serialize, Deserialize, PartialEq, Eq, Default)]
#[serde(rename_all = "snake_case")]
/// Built-in physical carrier shapes for rack/plate placement.
pub enum RackProfileKind {
    #[default]
    SmallTube4x6,
    Plate96,
    Plate384,
    Custom,
}

impl RackProfileKind {
    pub fn as_str(self) -> &'static str {
        match self {
            Self::SmallTube4x6 => "small_tube_4x6",
            Self::Plate96 => "plate_96",
            Self::Plate384 => "plate_384",
            Self::Custom => "custom",
        }
    }

    pub fn dimensions(self) -> (usize, usize) {
        match self {
            Self::SmallTube4x6 => (4, 6),
            Self::Plate96 => (8, 12),
            Self::Plate384 => (16, 24),
            Self::Custom => (0, 0),
        }
    }

    pub fn capacity(self) -> usize {
        let (rows, columns) = self.dimensions();
        rows * columns
    }
}

#[derive(Debug, Clone, Copy, Serialize, Deserialize, PartialEq, Eq, Default)]
#[serde(rename_all = "snake_case")]
/// Deterministic SVG sheet layouts for rack/arrangement label export.
pub enum RackLabelSheetPreset {
    #[default]
    CompactCards,
    PrintA4,
    WideCards,
}

impl RackLabelSheetPreset {
    pub fn as_str(self) -> &'static str {
        match self {
            Self::CompactCards => "compact_cards",
            Self::PrintA4 => "print_a4",
            Self::WideCards => "wide_cards",
        }
    }
}

#[derive(Debug, Clone, Copy, Serialize, Deserialize, PartialEq, Eq, Default)]
#[serde(rename_all = "snake_case")]
/// Deterministic SVG layouts for carrier-matched front-strip/module label export.
pub enum RackCarrierLabelPreset {
    #[default]
    FrontStripAndCards,
    FrontStripOnly,
    ModuleCardsOnly,
}

impl RackCarrierLabelPreset {
    pub fn as_str(self) -> &'static str {
        match self {
            Self::FrontStripAndCards => "front_strip_and_cards",
            Self::FrontStripOnly => "front_strip_only",
            Self::ModuleCardsOnly => "module_cards_only",
        }
    }
}

#[derive(Debug, Clone, Copy, Serialize, Deserialize, PartialEq, Eq, Default)]
#[serde(rename_all = "snake_case")]
/// Built-in printable physical carrier families layered on top of rack placement.
pub enum RackPhysicalTemplateKind {
    #[default]
    StoragePcrTubeRack,
    PipettingPcrTubeRack,
}

impl RackPhysicalTemplateKind {
    pub fn as_str(self) -> &'static str {
        match self {
            Self::StoragePcrTubeRack => "storage_pcr_tube_rack",
            Self::PipettingPcrTubeRack => "pipetting_pcr_tube_rack",
        }
    }
}

#[derive(Debug, Clone, Copy, Serialize, Deserialize, PartialEq, Eq)]
#[serde(rename_all = "snake_case")]
/// Physical handling intent for printable carrier exports.
pub enum RackPhysicalTemplateFamily {
    Storage,
    Pipetting,
}

impl RackPhysicalTemplateFamily {
    pub fn as_str(self) -> &'static str {
        match self {
            Self::Storage => "storage",
            Self::Pipetting => "pipetting",
        }
    }
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(default)]
/// Deterministic printable/fabrication geometry derived from one rack snapshot.
pub struct RackPhysicalTemplateSpec {
    pub kind: RackPhysicalTemplateKind,
    pub family: RackPhysicalTemplateFamily,
    pub container_format: String,
    pub rows: usize,
    pub columns: usize,
    pub pitch_x_mm: f32,
    pub pitch_y_mm: f32,
    pub opening_diameter_mm: f32,
    pub inner_wall_mm: f32,
    pub outer_wall_mm: f32,
    pub floor_thickness_mm: f32,
    pub rack_height_mm: f32,
    pub edge_margin_mm: f32,
    pub corner_radius_mm: f32,
    pub front_top_clearance_mm: f32,
    pub front_label_strip_depth_mm: f32,
    pub front_label_strip_recess_mm: f32,
    pub overall_width_mm: f32,
    pub overall_depth_mm: f32,
}

impl Default for RackPhysicalTemplateSpec {
    fn default() -> Self {
        Self {
            kind: RackPhysicalTemplateKind::StoragePcrTubeRack,
            family: RackPhysicalTemplateFamily::Storage,
            container_format: "pcr_tube_0_2ml".to_string(),
            rows: 0,
            columns: 0,
            pitch_x_mm: 0.0,
            pitch_y_mm: 0.0,
            opening_diameter_mm: 0.0,
            inner_wall_mm: 0.0,
            outer_wall_mm: 0.0,
            floor_thickness_mm: 0.0,
            rack_height_mm: 0.0,
            edge_margin_mm: 0.0,
            corner_radius_mm: 0.0,
            front_top_clearance_mm: 0.0,
            front_label_strip_depth_mm: 0.0,
            front_label_strip_recess_mm: 0.0,
            overall_width_mm: 0.0,
            overall_depth_mm: 0.0,
        }
    }
}

#[derive(Debug, Clone, Copy, Serialize, Deserialize, PartialEq, Eq, Default)]
#[serde(rename_all = "snake_case")]
/// Engine-owned quick authoring templates for common rack/plate setup styles.
pub enum RackAuthoringTemplate {
    #[default]
    BenchRows,
    PlateColumns,
    PlateEdgeAvoidance,
}

impl RackAuthoringTemplate {
    pub fn as_str(self) -> &'static str {
        match self {
            Self::BenchRows => "bench_rows",
            Self::PlateColumns => "plate_columns",
            Self::PlateEdgeAvoidance => "plate_edge_avoidance",
        }
    }
}

#[derive(Debug, Clone, Copy, Serialize, Deserialize, PartialEq, Eq, Default)]
#[serde(rename_all = "snake_case")]
/// Deterministic fill policy for physical rack/plate placement.
pub enum RackFillDirection {
    #[default]
    RowMajor,
    ColumnMajor,
}

impl RackFillDirection {
    pub fn as_str(self) -> &'static str {
        match self {
            Self::RowMajor => "row_major",
            Self::ColumnMajor => "column_major",
        }
    }
}

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq)]
#[serde(default)]
/// Snapshot of the physical carrier profile used when a rack was created.
pub struct RackProfileSnapshot {
    pub kind: RackProfileKind,
    pub rows: usize,
    pub columns: usize,
    pub coordinate_scheme: String,
    pub fill_direction: RackFillDirection,
    pub blocked_coordinates: Vec<String>,
}

impl Default for RackProfileSnapshot {
    fn default() -> Self {
        Self::from_kind(RackProfileKind::SmallTube4x6)
    }
}

impl RackProfileSnapshot {
    pub fn from_kind(kind: RackProfileKind) -> Self {
        let (rows, columns) = kind.dimensions();
        Self {
            kind,
            rows,
            columns,
            coordinate_scheme: "a1".to_string(),
            fill_direction: RackFillDirection::RowMajor,
            blocked_coordinates: vec![],
        }
    }

    pub fn custom(rows: usize, columns: usize) -> Self {
        Self {
            kind: RackProfileKind::Custom,
            rows,
            columns,
            coordinate_scheme: "a1".to_string(),
            fill_direction: RackFillDirection::RowMajor,
            blocked_coordinates: vec![],
        }
    }

    pub fn capacity(&self) -> usize {
        self.rows.saturating_mul(self.columns)
    }
}

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq)]
#[serde(tag = "kind", rename_all = "snake_case")]
/// Physical occupant placed in one rack coordinate.
pub enum RackOccupant {
    Container { container_id: ContainerId },
    LadderReference { ladder_name: String },
}

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq, Default)]
#[serde(default)]
/// One occupied physical coordinate in a saved rack.
pub struct RackPlacementEntry {
    pub coordinate: String,
    pub occupant: Option<RackOccupant>,
    pub arrangement_id: String,
    pub order_index: usize,
    pub role_label: String,
}

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq, Default)]
#[serde(default)]
/// Saved physical rack/plate draft that can host one or more arrangements.
pub struct Rack {
    pub rack_id: String,
    pub name: String,
    pub profile: RackProfileSnapshot,
    pub placements: Vec<RackPlacementEntry>,
    pub created_by_op: Option<OpId>,
    pub created_at_unix_ms: u128,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
/// Named arrangement definition referencing container lanes.
pub struct Arrangement {
    pub arrangement_id: String,
    pub mode: ArrangementMode,
    pub name: Option<String>,
    pub lane_container_ids: Vec<ContainerId>,
    #[serde(default)]
    pub ladders: Vec<String>,
    #[serde(default)]
    pub lane_role_labels: Vec<String>,
    #[serde(default)]
    pub default_rack_id: Option<String>,
    pub created_by_op: Option<OpId>,
    pub created_at_unix_ms: u128,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
/// Container/arrangement/rack state persisted with the project.
pub struct ContainerState {
    pub containers: HashMap<ContainerId, Container>,
    pub arrangements: HashMap<String, Arrangement>,
    pub racks: HashMap<String, Rack>,
    pub seq_to_latest_container: HashMap<SeqId, ContainerId>,
    pub next_container_counter: u64,
    pub next_arrangement_counter: u64,
    pub next_rack_counter: u64,
}

pub const TFBS_EXPERT_INSTRUCTION: &str = "TFBS expert view: each column is one PSSM position. Bar height is information content (2 - entropy in bits) from column base frequencies. Colored segments show A/C/G/T relative frequencies; the black polyline marks the matched base across positions.";

pub const RESTRICTION_EXPERT_INSTRUCTION: &str = "Restriction-site expert view: top strand is 5'->3', bottom strand is complementary 3'->5'. Aligned cut markers indicate a blunt cut; offset top/bottom markers indicate staggered sticky-end cleavage for the selected enzyme/site.";

pub const SPLICING_EXPERT_INSTRUCTION: &str = "Splicing expert view: one lane per transcript on a shared genomic axis. Exon geometry is coordinate-true (labels never resize exon/intron footprints). Donor/acceptor splice boundaries are marked, junction arcs summarize support across transcripts, and the transcript-vs-exon matrix shows isoform differences.";

pub const ISOFORM_ARCHITECTURE_EXPERT_INSTRUCTION: &str = "Isoform architecture view: top panel shows transcript/exon or transcript/CDS structure on genomic coordinates with 5'->3' orientation left-to-right (strand-aware axis), bottom panel shows per-isoform protein-domain architecture on amino-acid coordinates. Row order is shared across both panels; CDS-to-protein guide lines indicate which coding segments contribute to which amino-acid spans.";

#[derive(Debug, Clone, Copy, Serialize, Deserialize, PartialEq, Eq, Default)]
#[serde(rename_all = "snake_case")]
/// Dotplot comparison mode for dotplot computation/export.
pub enum DotplotMode {
    #[default]
    SelfForward,
    SelfReverseComplement,
    PairForward,
    PairReverseComplement,
}

impl DotplotMode {
    pub fn as_str(self) -> &'static str {
        match self {
            Self::SelfForward => "self_forward",
            Self::SelfReverseComplement => "self_reverse_complement",
            Self::PairForward => "pair_forward",
            Self::PairReverseComplement => "pair_reverse_complement",
        }
    }
}

#[derive(Debug, Clone, Copy, Serialize, Deserialize, PartialEq, Eq, Default)]
#[serde(rename_all = "snake_case")]
/// Pairwise alignment mode for sequence and confirmation alignments.
pub enum PairwiseAlignmentMode {
    #[default]
    Global,
    Local,
}

impl PairwiseAlignmentMode {
    pub fn as_str(self) -> &'static str {
        match self {
            Self::Global => "global",
            Self::Local => "local",
        }
    }
}

#[derive(Debug, Clone, Copy, Serialize, Deserialize, PartialEq, Eq, Default)]
#[serde(rename_all = "snake_case")]
/// Flexibility score model for flexibility-track computation.
pub enum FlexibilityModel {
    #[default]
    AtRichness,
    AtSkew,
}

impl FlexibilityModel {
    pub fn as_str(self) -> &'static str {
        match self {
            Self::AtRichness => "at_richness",
            Self::AtSkew => "at_skew",
        }
    }
}

#[derive(Debug, Clone, Copy, Serialize, Deserialize, PartialEq, Eq, Default)]
#[serde(rename_all = "snake_case")]
/// Shared electrophoresis buffer preset for virtual gel rendering.
pub enum GelBufferModel {
    #[default]
    Tae,
    Tbe,
}

impl GelBufferModel {
    pub fn as_str(self) -> &'static str {
        match self {
            Self::Tae => "tae",
            Self::Tbe => "tbe",
        }
    }
}

#[derive(Debug, Clone, Copy, Serialize, Deserialize, PartialEq, Eq, Default)]
#[serde(rename_all = "snake_case")]
/// Shared gel-topology form hint for electrophoresis rendering.
///
/// This is deliberately narrower than the full sequence topology model. It
/// exists so gel previews/exports can distinguish common circular DNA forms
/// when that information is explicitly known or inferable, while older code
/// paths can still degrade to plain `linear` / generic `circular`.
pub enum GelTopologyForm {
    #[default]
    Linear,
    Circular,
    Supercoiled,
    RelaxedCircular,
    NickedCircular,
}

impl GelTopologyForm {
    pub fn from_hint(raw: &str) -> Option<Self> {
        let lowered = raw.trim().to_ascii_lowercase();
        match lowered.as_str() {
            "" => None,
            "linear" | "linearized" | "linearised" => Some(Self::Linear),
            "circular" => Some(Self::Circular),
            "supercoiled" | "superhelical" | "ccc" | "ccc dna" => Some(Self::Supercoiled),
            "relaxed" | "relaxed_circular" | "relaxed circular" => Some(Self::RelaxedCircular),
            "nicked" | "nicked_circular" | "nicked circular" | "open_circular"
            | "open circular" | "open-circle" => Some(Self::NickedCircular),
            _ => None,
        }
    }

    pub fn as_str(self) -> &'static str {
        match self {
            Self::Linear => "linear",
            Self::Circular => "circular",
            Self::Supercoiled => "supercoiled",
            Self::RelaxedCircular => "relaxed_circular",
            Self::NickedCircular => "nicked_circular",
        }
    }

    pub fn display_label(self) -> &'static str {
        match self {
            Self::Linear => "linear",
            Self::Circular => "circular",
            Self::Supercoiled => "supercoiled",
            Self::RelaxedCircular => "relaxed circular",
            Self::NickedCircular => "nicked circular",
        }
    }

    pub fn is_circular(self) -> bool {
        !matches!(self, Self::Linear)
    }
}

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
#[serde(default)]
/// Shared gel-run conditions used by GUI/CLI/render export.
///
/// This is intentionally a conservative heuristic contract, not a full
/// electrophoresis simulator. Adapters should pass this bundle unchanged to the
/// shared engine/render path so one deterministic model is reused everywhere.
pub struct GelRunConditions {
    pub agarose_percent: f32,
    pub buffer_model: GelBufferModel,
    pub topology_aware: bool,
}

impl Default for GelRunConditions {
    fn default() -> Self {
        Self {
            agarose_percent: 1.0,
            buffer_model: GelBufferModel::Tae,
            topology_aware: true,
        }
    }
}

impl GelRunConditions {
    pub fn normalized(&self) -> Self {
        let agarose_percent = if self.agarose_percent.is_finite() {
            self.agarose_percent.clamp(0.5, 3.0)
        } else {
            1.0
        };
        Self {
            agarose_percent,
            buffer_model: self.buffer_model,
            topology_aware: self.topology_aware,
        }
    }

    pub fn describe(&self) -> String {
        let normalized = self.normalized();
        format!(
            "{:.1}% agarose | {} | topology-aware {}",
            normalized.agarose_percent,
            normalized.buffer_model.as_str().to_ascii_uppercase(),
            if normalized.topology_aware {
                "on"
            } else {
                "off"
            }
        )
    }
}

#[derive(Debug, Clone, Copy, Serialize, Deserialize, PartialEq, Eq, Default)]
#[serde(rename_all = "snake_case")]
/// Shared scope preset for splicing/exon-context views and RNA-read mapping.
pub enum SplicingScopePreset {
    #[default]
    AllOverlappingBothStrands,
    TargetGroupAnyStrand,
    AllOverlappingTargetStrand,
    TargetGroupTargetStrand,
}

impl SplicingScopePreset {
    pub fn as_str(self) -> &'static str {
        match self {
            Self::AllOverlappingBothStrands => "all_overlapping_both_strands",
            Self::TargetGroupAnyStrand => "target_group_any_strand",
            Self::AllOverlappingTargetStrand => "all_overlapping_target_strand",
            Self::TargetGroupTargetStrand => "target_group_target_strand",
        }
    }

    pub fn restrict_to_target_group(self) -> bool {
        matches!(
            self,
            Self::TargetGroupAnyStrand | Self::TargetGroupTargetStrand
        )
    }

    pub fn restrict_to_target_strand(self) -> bool {
        matches!(
            self,
            Self::AllOverlappingTargetStrand | Self::TargetGroupTargetStrand
        )
    }
}

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq)]
#[serde(rename_all = "snake_case")]
pub enum FeatureExpertTarget {
    TfbsFeature {
        feature_id: usize,
    },
    RestrictionSite {
        cut_pos_1based: usize,
        #[serde(default)]
        enzyme: Option<String>,
        #[serde(default)]
        recognition_start_1based: Option<usize>,
        #[serde(default)]
        recognition_end_1based: Option<usize>,
    },
    SplicingFeature {
        feature_id: usize,
        #[serde(default)]
        scope: SplicingScopePreset,
    },
    IsoformArchitecture {
        panel_id: String,
    },
    UniprotProjection {
        projection_id: String,
    },
}

impl FeatureExpertTarget {
    pub fn describe(&self) -> String {
        match self {
            Self::TfbsFeature { feature_id } => format!("tfbs feature #{feature_id}"),
            Self::RestrictionSite {
                cut_pos_1based,
                enzyme,
                recognition_start_1based,
                recognition_end_1based,
            } => {
                let enzyme = enzyme
                    .as_deref()
                    .map(str::trim)
                    .filter(|v| !v.is_empty())
                    .unwrap_or("*");
                let mut out = format!("restriction cut@{cut_pos_1based} enzyme={enzyme}");
                if let (Some(start), Some(end)) = (recognition_start_1based, recognition_end_1based)
                {
                    out.push_str(&format!(" range={start}..{end}"));
                }
                out
            }
            Self::SplicingFeature { feature_id, scope } => {
                format!("splicing feature #{feature_id} scope={}", scope.as_str())
            }
            Self::IsoformArchitecture { panel_id } => {
                format!("isoform architecture panel '{panel_id}'")
            }
            Self::UniprotProjection { projection_id } => {
                format!("UniProt projection '{projection_id}'")
            }
        }
    }
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SplicingRange {
    pub start_1based: usize,
    pub end_1based: usize,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SplicingExonCdsPhase {
    pub start_1based: usize,
    pub end_1based: usize,
    #[serde(default)]
    pub left_cds_phase: Option<u8>,
    #[serde(default)]
    pub right_cds_phase: Option<u8>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SplicingExonSummary {
    pub start_1based: usize,
    pub end_1based: usize,
    pub support_transcript_count: usize,
    pub constitutive: bool,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SplicingBoundaryMarker {
    pub transcript_feature_id: usize,
    pub transcript_id: String,
    pub side: String,
    pub position_1based: usize,
    pub motif_2bp: String,
    pub canonical: bool,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SplicingJunctionArc {
    pub donor_1based: usize,
    pub acceptor_1based: usize,
    pub support_transcript_count: usize,
    pub transcript_feature_ids: Vec<usize>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SplicingTranscriptLane {
    pub transcript_feature_id: usize,
    pub transcript_id: String,
    pub label: String,
    pub strand: String,
    pub exons: Vec<SplicingRange>,
    #[serde(default)]
    pub exon_cds_phases: Vec<SplicingExonCdsPhase>,
    pub introns: Vec<SplicingRange>,
    pub has_target_feature: bool,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SplicingMatrixRow {
    pub transcript_feature_id: usize,
    pub transcript_id: String,
    pub label: String,
    pub exon_presence: Vec<bool>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SplicingEventSummary {
    pub event_type: String,
    pub count: usize,
    pub details: Vec<String>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SplicingExpertView {
    pub seq_id: String,
    pub target_feature_id: usize,
    pub group_label: String,
    pub strand: String,
    pub region_start_1based: usize,
    pub region_end_1based: usize,
    pub transcript_count: usize,
    pub unique_exon_count: usize,
    pub instruction: String,
    pub transcripts: Vec<SplicingTranscriptLane>,
    pub unique_exons: Vec<SplicingExonSummary>,
    pub matrix_rows: Vec<SplicingMatrixRow>,
    pub boundaries: Vec<SplicingBoundaryMarker>,
    pub junctions: Vec<SplicingJunctionArc>,
    pub events: Vec<SplicingEventSummary>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct IsoformArchitectureTranscriptLane {
    pub isoform_id: String,
    pub label: String,
    #[serde(default)]
    pub transcript_id: Option<String>,
    #[serde(default)]
    pub transcript_feature_id: Option<usize>,
    pub strand: String,
    #[serde(default)]
    pub transcript_exons: Vec<SplicingRange>,
    pub exons: Vec<SplicingRange>,
    pub introns: Vec<SplicingRange>,
    pub mapped: bool,
    #[serde(default)]
    pub transactivation_class: Option<String>,
    #[serde(default)]
    pub cds_to_protein_segments: Vec<IsoformArchitectureCdsAaSegment>,
    #[serde(default)]
    pub note: Option<String>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct IsoformArchitectureCdsAaSegment {
    pub genomic_start_1based: usize,
    pub genomic_end_1based: usize,
    pub aa_start: usize,
    pub aa_end: usize,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct IsoformArchitectureProteinDomain {
    pub name: String,
    pub start_aa: usize,
    pub end_aa: usize,
    #[serde(default)]
    pub color_hex: Option<String>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct IsoformArchitectureProteinLane {
    pub isoform_id: String,
    pub label: String,
    #[serde(default)]
    pub transcript_id: Option<String>,
    #[serde(default)]
    pub expected_length_aa: Option<usize>,
    #[serde(default)]
    pub reference_start_aa: Option<usize>,
    #[serde(default)]
    pub reference_end_aa: Option<usize>,
    pub domains: Vec<IsoformArchitectureProteinDomain>,
    #[serde(default)]
    pub transactivation_class: Option<String>,
}

fn default_isoform_transcript_geometry_mode() -> String {
    "exon".to_string()
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct IsoformArchitectureExpertView {
    pub seq_id: String,
    pub panel_id: String,
    pub gene_symbol: String,
    #[serde(default = "default_isoform_transcript_geometry_mode")]
    pub transcript_geometry_mode: String,
    #[serde(default)]
    pub panel_source: Option<String>,
    pub region_start_1based: usize,
    pub region_end_1based: usize,
    pub instruction: String,
    pub transcript_lanes: Vec<IsoformArchitectureTranscriptLane>,
    pub protein_lanes: Vec<IsoformArchitectureProteinLane>,
    #[serde(default)]
    pub warnings: Vec<String>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct TfbsExpertColumn {
    pub index_1based: usize,
    pub counts: [f64; 4],
    pub frequencies: [f64; 4],
    pub information_content_bits: f64,
    pub match_base: Option<char>,
    pub match_frequency: Option<f64>,
    pub llr_bits: Option<f64>,
    pub true_log_odds_bits: Option<f64>,
    pub llr_rank_desc: Option<usize>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct TfbsExpertView {
    pub seq_id: String,
    pub feature_id: usize,
    pub feature_label: String,
    pub tf_id: String,
    #[serde(default)]
    pub tf_name: Option<String>,
    pub strand: String,
    pub start_1based: usize,
    pub end_1based: usize,
    pub motif_length: usize,
    pub matched_sequence: String,
    #[serde(default)]
    pub llr_total_bits: Option<f64>,
    #[serde(default)]
    pub llr_quantile: Option<f64>,
    #[serde(default)]
    pub true_log_odds_total_bits: Option<f64>,
    #[serde(default)]
    pub true_log_odds_quantile: Option<f64>,
    pub instruction: String,
    pub columns: Vec<TfbsExpertColumn>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct RestrictionSiteExpertView {
    pub seq_id: String,
    pub cut_pos_1based: usize,
    #[serde(default)]
    pub paired_cut_pos_1based: usize,
    pub recognition_start_1based: usize,
    pub recognition_end_1based: usize,
    pub cut_index_0based: usize,
    #[serde(default)]
    pub paired_cut_index_0based: usize,
    #[serde(default)]
    pub end_geometry: String,
    pub number_of_cuts_for_enzyme: usize,
    #[serde(default)]
    pub selected_enzyme: Option<String>,
    pub enzyme_names: Vec<String>,
    #[serde(default)]
    pub recognition_iupac: Option<String>,
    pub site_sequence: String,
    pub site_sequence_complement: String,
    #[serde(default)]
    pub enzyme_cut_offset_0based: Option<isize>,
    #[serde(default)]
    pub overlap_bp: Option<isize>,
    #[serde(default)]
    pub enzyme_note: Option<String>,
    #[serde(default)]
    pub rebase_url: Option<String>,
    pub instruction: String,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(tag = "kind", content = "data", rename_all = "snake_case")]
pub enum FeatureExpertView {
    Tfbs(TfbsExpertView),
    RestrictionSite(RestrictionSiteExpertView),
    Splicing(SplicingExpertView),
    IsoformArchitecture(IsoformArchitectureExpertView),
}

impl FeatureExpertView {
    pub fn instruction(&self) -> &str {
        match self {
            Self::Tfbs(_) => TFBS_EXPERT_INSTRUCTION,
            Self::RestrictionSite(_) => RESTRICTION_EXPERT_INSTRUCTION,
            Self::Splicing(_) => SPLICING_EXPERT_INSTRUCTION,
            Self::IsoformArchitecture(_) => ISOFORM_ARCHITECTURE_EXPERT_INSTRUCTION,
        }
    }
}

#[derive(Debug, Clone, Copy, Serialize, Deserialize, PartialEq, Eq, Default)]
#[serde(rename_all = "snake_case")]
pub enum RnaReadInputFormat {
    #[default]
    Fasta,
}

impl RnaReadInputFormat {
    pub fn as_str(self) -> &'static str {
        match self {
            Self::Fasta => "fasta",
        }
    }
}

#[derive(Debug, Clone, Copy, Serialize, Deserialize, PartialEq, Eq, Default)]
#[serde(rename_all = "snake_case")]
pub enum RnaReadInterpretationProfile {
    #[default]
    NanoporeCdnaV1,
    ShortReadV1,
    TransposonV1,
}

impl RnaReadInterpretationProfile {
    pub fn as_str(self) -> &'static str {
        match self {
            Self::NanoporeCdnaV1 => "nanopore_cdna_v1",
            Self::ShortReadV1 => "short_read_v1",
            Self::TransposonV1 => "transposon_v1",
        }
    }
}

#[derive(Debug, Clone, Copy, Serialize, Deserialize, PartialEq, Eq, Default)]
#[serde(rename_all = "snake_case")]
pub enum RnaReadReportMode {
    #[default]
    Full,
    SeedPassedOnly,
}

impl RnaReadReportMode {
    pub fn as_str(self) -> &'static str {
        match self {
            Self::Full => "full",
            Self::SeedPassedOnly => "seed_passed_only",
        }
    }
}

#[derive(Debug, Clone, Copy, Serialize, Deserialize, PartialEq, Eq, Default)]
#[serde(rename_all = "snake_case")]
pub enum RnaReadOriginMode {
    #[default]
    SingleGene,
    MultiGeneSparse,
}

impl RnaReadOriginMode {
    pub fn as_str(self) -> &'static str {
        match self {
            Self::SingleGene => "single_gene",
            Self::MultiGeneSparse => "multi_gene_sparse",
        }
    }
}

#[derive(Debug, Clone, Copy, Serialize, Deserialize, PartialEq, Eq, Default)]
#[serde(rename_all = "snake_case")]
pub enum RnaReadHitSelection {
    All,
    SeedPassed,
    #[default]
    Aligned,
}

impl RnaReadHitSelection {
    pub fn as_str(self) -> &'static str {
        match self {
            Self::All => "all",
            Self::SeedPassed => "seed_passed",
            Self::Aligned => "aligned",
        }
    }
}

#[derive(Debug, Clone, Copy, Serialize, Deserialize, PartialEq, Eq, Default)]
#[serde(rename_all = "snake_case")]
pub enum RnaReadGeneSupportCompleteRule {
    #[default]
    Near,
    Strict,
    Exact,
}

impl RnaReadGeneSupportCompleteRule {
    pub fn as_str(self) -> &'static str {
        match self {
            Self::Near => "near",
            Self::Strict => "strict",
            Self::Exact => "exact",
        }
    }
}

#[derive(Debug, Clone, Copy, Serialize, Deserialize, PartialEq, Eq, Default)]
#[serde(rename_all = "snake_case")]
pub enum RnaReadGeneSupportAuditStatus {
    #[default]
    Unaligned,
    AlignedOtherGene,
    AcceptedFragment,
    AcceptedComplete,
}

impl RnaReadGeneSupportAuditStatus {
    pub fn as_str(self) -> &'static str {
        match self {
            Self::Unaligned => "unaligned",
            Self::AlignedOtherGene => "aligned_other_gene",
            Self::AcceptedFragment => "accepted_fragment",
            Self::AcceptedComplete => "accepted_complete",
        }
    }
}

#[derive(Debug, Clone, Copy, Serialize, Deserialize, PartialEq, Eq, Default)]
#[serde(rename_all = "snake_case")]
pub enum RnaReadGeneSupportAuditCohortFilter {
    #[default]
    All,
    Accepted,
    Fragment,
    Complete,
    Rejected,
}

impl RnaReadGeneSupportAuditCohortFilter {
    pub fn as_str(self) -> &'static str {
        match self {
            Self::All => "all",
            Self::Accepted => "accepted",
            Self::Fragment => "fragment",
            Self::Complete => "complete",
            Self::Rejected => "rejected",
        }
    }
}

#[derive(Debug, Clone, Copy, Serialize, Deserialize, PartialEq, Eq, Default)]
#[serde(rename_all = "snake_case")]
pub enum RnaReadScoreDensityScale {
    Linear,
    #[default]
    Log,
}

impl RnaReadScoreDensityScale {
    pub fn as_str(self) -> &'static str {
        match self {
            Self::Linear => "linear",
            Self::Log => "log",
        }
    }
}

#[derive(Debug, Clone, Copy, Serialize, Deserialize, PartialEq, Eq, Default)]
#[serde(rename_all = "snake_case")]
pub enum RnaReadScoreDensityVariant {
    #[default]
    AllScored,
    CompositeSeedGate,
    RetainedReplayCurrentControls,
}

impl RnaReadScoreDensityVariant {
    pub fn as_str(self) -> &'static str {
        match self {
            Self::AllScored => "all_scored",
            Self::CompositeSeedGate => "composite_seed_gate",
            Self::RetainedReplayCurrentControls => "retained_replay_current_controls",
        }
    }
}

#[derive(Debug, Clone, Copy, Serialize, Deserialize, PartialEq, Eq, Default)]
#[serde(rename_all = "snake_case")]
pub enum RnaReadAlignmentMode {
    #[default]
    Local,
    Semiglobal,
}

impl RnaReadAlignmentMode {
    pub fn as_str(&self) -> &'static str {
        match self {
            Self::Local => "local",
            Self::Semiglobal => "semiglobal",
        }
    }
}

#[derive(Debug, Clone, Copy, Serialize, Deserialize, PartialEq, Eq, Default)]
#[serde(rename_all = "snake_case")]
pub enum RnaReadOriginClass {
    TargetCoherent,
    TargetPartialLocalBlock,
    RoiSameStrandLocalBlock,
    RoiReverseStrandLocalBlock,
    TpFamilyAmbiguous,
    #[default]
    BackgroundLikely,
}

impl RnaReadOriginClass {
    pub fn as_str(&self) -> &'static str {
        match self {
            Self::TargetCoherent => "target_coherent",
            Self::TargetPartialLocalBlock => "target_partial_local_block",
            Self::RoiSameStrandLocalBlock => "roi_same_strand_local_block",
            Self::RoiReverseStrandLocalBlock => "roi_reverse_strand_local_block",
            Self::TpFamilyAmbiguous => "tp_family_ambiguous",
            Self::BackgroundLikely => "background_likely",
        }
    }
}

#[derive(Debug, Clone, Copy, Serialize, Deserialize, PartialEq, Eq, Default)]
#[serde(rename_all = "snake_case")]
pub enum RnaReadAlignmentBackend {
    #[default]
    Banded,
    DenseFallback,
}

impl RnaReadAlignmentBackend {
    pub fn as_str(self) -> &'static str {
        match self {
            Self::Banded => "banded",
            Self::DenseFallback => "dense_fallback",
        }
    }
}

#[derive(Debug, Clone, Copy, Serialize, Deserialize, PartialEq, Eq, Default)]
#[serde(rename_all = "snake_case")]
pub enum RnaReadAlignmentEffect {
    #[default]
    ConfirmedAssignment,
    ReassignedTranscript,
    AlignedWithoutPhase1Assignment,
}

impl RnaReadAlignmentEffect {
    pub fn as_str(self) -> &'static str {
        match self {
            Self::ConfirmedAssignment => "confirmed_assignment",
            Self::ReassignedTranscript => "reassigned_transcript",
            Self::AlignedWithoutPhase1Assignment => "aligned_without_phase1_assignment",
        }
    }
}

#[derive(Debug, Clone, Copy, Serialize, Deserialize, PartialEq, Eq, Default)]
#[serde(rename_all = "snake_case")]
pub enum RnaReadAlignmentInspectionEffectFilter {
    #[default]
    AllAligned,
    ConfirmedOnly,
    DisagreementOnly,
    ReassignedOnly,
    NoPhase1Only,
    SelectedOnly,
}

impl RnaReadAlignmentInspectionEffectFilter {
    pub fn as_str(self) -> &'static str {
        match self {
            Self::AllAligned => "all_aligned",
            Self::ConfirmedOnly => "confirmed_only",
            Self::DisagreementOnly => "disagreement_only",
            Self::ReassignedOnly => "reassigned_only",
            Self::NoPhase1Only => "no_phase1_only",
            Self::SelectedOnly => "selected_only",
        }
    }
}

#[derive(Debug, Clone, Copy, Serialize, Deserialize, PartialEq, Eq, Default)]
#[serde(rename_all = "snake_case")]
pub enum RnaReadAlignmentInspectionSortKey {
    #[default]
    Rank,
    Identity,
    Coverage,
    Score,
}

impl RnaReadAlignmentInspectionSortKey {
    pub fn as_str(self) -> &'static str {
        match self {
            Self::Rank => "rank",
            Self::Identity => "identity",
            Self::Coverage => "coverage",
            Self::Score => "score",
        }
    }
}

#[derive(Debug, Clone, Copy, Serialize, Deserialize, PartialEq, Eq)]
/// Stable high-level error class for engine operation failures.
pub enum ErrorCode {
    InvalidInput,
    NotFound,
    Unsupported,
    Io,
    Internal,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
/// Shared structured error payload returned by engine operations/adapters.
pub struct EngineError {
    pub code: ErrorCode,
    pub message: String,
}

impl fmt::Display for EngineError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{:?}: {}", self.code, self.message)
    }
}

impl Error for EngineError {}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
pub struct DotplotMatchPoint {
    pub x_0based: usize,
    pub y_0based: usize,
    pub mismatches: usize,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
pub struct DotplotBoxplotBin {
    pub query_start_0based: usize,
    pub query_end_0based_exclusive: usize,
    pub hit_count: usize,
    pub min_reference_0based: Option<usize>,
    pub q1_reference_0based: Option<usize>,
    pub median_reference_0based: Option<usize>,
    pub q3_reference_0based: Option<usize>,
    pub max_reference_0based: Option<usize>,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
pub struct DotplotQuerySeries {
    pub series_id: String,
    pub seq_id: String,
    pub label: String,
    pub color_rgb: [u8; 3],
    #[serde(default)]
    pub mode: DotplotMode,
    pub span_start_0based: usize,
    pub span_end_0based: usize,
    pub point_count: usize,
    pub points: Vec<DotplotMatchPoint>,
    pub boxplot_bin_count: usize,
    pub boxplot_bins: Vec<DotplotBoxplotBin>,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
pub struct DotplotReferenceAnnotationInterval {
    pub start_0based: usize,
    pub end_0based_exclusive: usize,
    pub label: String,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
pub struct DotplotReferenceAnnotationTrack {
    pub seq_id: String,
    pub label: String,
    pub interval_count: usize,
    pub intervals: Vec<DotplotReferenceAnnotationInterval>,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
pub struct DotplotOverlayQuerySpec {
    pub seq_id: String,
    pub label: String,
    #[serde(default)]
    pub span_start_0based: Option<usize>,
    #[serde(default)]
    pub span_end_0based: Option<usize>,
    #[serde(default)]
    pub mode: DotplotMode,
    #[serde(default)]
    pub color_rgb: Option<[u8; 3]>,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
pub struct DotplotView {
    pub schema: String,
    pub dotplot_id: String,
    pub owner_seq_id: String,
    pub seq_id: String,
    pub reference_seq_id: Option<String>,
    pub generated_at_unix_ms: u128,
    pub span_start_0based: usize,
    pub span_end_0based: usize,
    pub reference_span_start_0based: usize,
    pub reference_span_end_0based: usize,
    pub mode: DotplotMode,
    pub word_size: usize,
    pub step_bp: usize,
    pub max_mismatches: usize,
    pub tile_bp: Option<usize>,
    pub point_count: usize,
    pub points: Vec<DotplotMatchPoint>,
    pub boxplot_bin_count: usize,
    pub boxplot_bins: Vec<DotplotBoxplotBin>,
    pub series_count: usize,
    pub query_series: Vec<DotplotQuerySeries>,
    #[serde(default)]
    pub reference_annotation: Option<DotplotReferenceAnnotationTrack>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct DotplotViewSummary {
    pub dotplot_id: String,
    pub owner_seq_id: String,
    pub seq_id: String,
    pub reference_seq_id: Option<String>,
    pub generated_at_unix_ms: u128,
    pub span_start_0based: usize,
    pub span_end_0based: usize,
    pub reference_span_start_0based: usize,
    pub reference_span_end_0based: usize,
    pub mode: DotplotMode,
    pub word_size: usize,
    pub step_bp: usize,
    pub max_mismatches: usize,
    pub point_count: usize,
    pub series_count: usize,
}

impl DotplotView {
    pub fn normalize_v3_defaults(&mut self) {
        if self.owner_seq_id.trim().is_empty() {
            self.owner_seq_id = self.seq_id.clone();
        }
        if self.query_series.is_empty() {
            let label = if self.seq_id.trim().is_empty() {
                "<query>".to_string()
            } else {
                self.seq_id.clone()
            };
            self.query_series.push(DotplotQuerySeries {
                series_id: if self.dotplot_id.trim().is_empty() {
                    "series_1".to_string()
                } else {
                    format!("{}_series_1", self.dotplot_id)
                },
                seq_id: self.seq_id.clone(),
                label,
                color_rgb: [29, 78, 216],
                mode: self.mode,
                span_start_0based: self.span_start_0based,
                span_end_0based: self.span_end_0based,
                point_count: if self.point_count == 0 {
                    self.points.len()
                } else {
                    self.point_count
                },
                points: self.points.clone(),
                boxplot_bin_count: if self.boxplot_bin_count == 0 {
                    self.boxplot_bins.len()
                } else {
                    self.boxplot_bin_count
                },
                boxplot_bins: self.boxplot_bins.clone(),
            });
        }
        self.series_count = self.query_series.len();
        if let Some(primary) = self.query_series.first() {
            if self.seq_id.trim().is_empty() {
                self.seq_id = primary.seq_id.clone();
            }
            self.mode = primary.mode;
            self.span_start_0based = primary.span_start_0based;
            self.span_end_0based = primary.span_end_0based;
            self.point_count = primary.point_count.max(primary.points.len());
            self.points = primary.points.clone();
            self.boxplot_bin_count = primary.boxplot_bin_count.max(primary.boxplot_bins.len());
            self.boxplot_bins = primary.boxplot_bins.clone();
        }
        if let Some(annotation) = self.reference_annotation.as_mut() {
            annotation.interval_count = annotation.intervals.len();
        }
    }

    pub fn primary_series(&self) -> Option<&DotplotQuerySeries> {
        self.query_series.first()
    }
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
pub struct SequenceAlignmentReport {
    pub schema: String,
    pub mode: PairwiseAlignmentMode,
    pub query_seq_id: String,
    pub target_seq_id: String,
    pub query_span_start_0based: usize,
    pub query_span_end_0based: usize,
    pub target_span_start_0based: usize,
    pub target_span_end_0based: usize,
    pub aligned_query_start_0based: usize,
    pub aligned_query_end_0based_exclusive: usize,
    pub aligned_target_start_0based: usize,
    pub aligned_target_end_0based_exclusive: usize,
    pub score: i32,
    pub match_score: i32,
    pub mismatch_score: i32,
    pub gap_open: i32,
    pub gap_extend: i32,
    pub aligned_columns: usize,
    pub matches: usize,
    pub mismatches: usize,
    pub insertions: usize,
    pub deletions: usize,
    pub identity_fraction: f64,
    pub query_coverage_fraction: f64,
    pub target_coverage_fraction: f64,
    pub cigar: String,
}

#[derive(Debug, Clone, Copy, Serialize, Deserialize, PartialEq, Eq, Default)]
#[serde(rename_all = "snake_case")]
/// Raw sequencing-trace file format supported by sequencing-evidence intake.
pub enum SequencingTraceFormat {
    #[default]
    AbiAb1,
    Scf,
}

impl SequencingTraceFormat {
    pub fn as_str(self) -> &'static str {
        match self {
            Self::AbiAb1 => "abi_ab1",
            Self::Scf => "scf",
        }
    }
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
/// One raw-channel availability summary for an imported sequencing trace.
pub struct SequencingTraceChannelSummary {
    pub channel: String,
    pub trace_set: String,
    pub point_count: usize,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
/// One stored chromatogram/intensity curve for an imported sequencing trace.
pub struct SequencingTraceChannelData {
    pub channel: String,
    pub trace_set: String,
    pub points: Vec<u32>,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
/// Persisted raw sequencing-trace evidence record.
///
/// This record stores the called bases and the key per-base evidence arrays
/// needed for later trace-aware confirmation without mutating project
/// sequences. Newer schema revisions may also carry raw chromatogram curves
/// and clip-window metadata for GUI inspection.
pub struct SequencingTraceRecord {
    pub schema: String,
    pub trace_id: String,
    pub format: SequencingTraceFormat,
    pub source_path: String,
    pub imported_at_unix_ms: u128,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub seq_id: Option<SeqId>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub sample_name: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub sample_well: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub run_name: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub machine_name: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub machine_model: Option<String>,
    pub called_bases: String,
    #[serde(default)]
    pub called_base_confidence_values: Vec<u8>,
    #[serde(default)]
    pub peak_locations: Vec<u32>,
    #[serde(default)]
    pub channel_data: Vec<SequencingTraceChannelData>,
    #[serde(default)]
    pub channel_summaries: Vec<SequencingTraceChannelSummary>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub clip_start_base_index: Option<usize>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub clip_end_base_index_exclusive: Option<usize>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub comments_text: Option<String>,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
/// Compact row used by adapters to list imported sequencing traces.
pub struct SequencingTraceSummary {
    pub trace_id: String,
    pub format: SequencingTraceFormat,
    pub source_path: String,
    pub imported_at_unix_ms: u128,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub seq_id: Option<SeqId>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub sample_name: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub run_name: Option<String>,
    pub called_base_count: usize,
    pub confidence_value_count: usize,
    pub peak_location_count: usize,
    pub has_curve_data: bool,
    pub channel_count: usize,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
/// Structured result emitted when one raw trace file is imported.
pub struct SequencingTraceImportReport {
    pub schema: String,
    pub trace_id: String,
    pub format: SequencingTraceFormat,
    pub source_path: String,
    pub imported_at_unix_ms: u128,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub seq_id: Option<SeqId>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub sample_name: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub run_name: Option<String>,
    pub called_base_count: usize,
    pub confidence_value_count: usize,
    pub peak_location_count: usize,
    pub has_curve_data: bool,
    pub channel_count: usize,
    #[serde(default)]
    pub warnings: Vec<String>,
}

#[derive(Debug, Clone, Copy, Serialize, Deserialize, PartialEq, Eq, Default)]
#[serde(rename_all = "snake_case")]
/// Which strand-direction a suggested sequencing primer is expected to read on
/// the expected construct.
pub enum SequencingPrimerOrientation {
    #[default]
    ForwardRead,
    ReverseRead,
}

impl SequencingPrimerOrientation {
    pub fn as_str(self) -> &'static str {
        match self {
            Self::ForwardRead => "forward_read",
            Self::ReverseRead => "reverse_read",
        }
    }
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
/// One candidate sequencing-primer hit on an expected construct.
pub struct SequencingPrimerOverlaySuggestion {
    pub primer_seq_id: SeqId,
    pub primer_label: String,
    pub primer_sequence: String,
    pub orientation: SequencingPrimerOrientation,
    pub anneal_sequence: String,
    pub anneal_start_0based: usize,
    pub anneal_end_0based_exclusive: usize,
    pub three_prime_position_0based: usize,
    pub predicted_read_span_start_0based: usize,
    pub predicted_read_span_end_0based_exclusive: usize,
    #[serde(default)]
    pub covered_target_ids: Vec<String>,
    #[serde(default)]
    pub covered_problem_target_ids: Vec<String>,
    #[serde(default)]
    pub covered_variant_ids: Vec<String>,
    #[serde(default)]
    pub covered_problem_variant_ids: Vec<String>,
}

#[derive(Debug, Clone, Copy, Serialize, Deserialize, PartialEq, Eq, Default)]
#[serde(rename_all = "snake_case")]
/// One unresolved confirmation locus category that can receive primer guidance.
pub enum SequencingPrimerProblemKind {
    #[default]
    Target,
    Variant,
}

impl SequencingPrimerProblemKind {
    pub fn as_str(self) -> &'static str {
        match self {
            Self::Target => "target",
            Self::Variant => "variant",
        }
    }
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
/// Guidance row recommending the best existing primer hit for one unresolved locus.
pub struct SequencingPrimerProblemGuidanceRow {
    pub problem_id: String,
    pub problem_kind: SequencingPrimerProblemKind,
    pub problem_label: String,
    pub problem_summary: String,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub recommended_primer_seq_id: Option<SeqId>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub recommended_primer_label: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub recommended_orientation: Option<SequencingPrimerOrientation>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub recommended_read_span_start_0based: Option<usize>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub recommended_read_span_end_0based_exclusive: Option<usize>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub recommended_three_prime_distance_bp: Option<usize>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub recommended_in_read_direction: Option<bool>,
    #[serde(default)]
    pub recommended_problem_target_count: usize,
    #[serde(default)]
    pub recommended_problem_variant_count: usize,
    #[serde(default)]
    pub candidate_count: usize,
    pub reason: String,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
/// One fresh sequencing-primer proposal for an unresolved confirmation locus.
pub struct SequencingPrimerProposalRow {
    pub proposal_id: String,
    pub problem_id: String,
    pub problem_kind: SequencingPrimerProblemKind,
    pub problem_label: String,
    pub problem_summary: String,
    pub orientation: SequencingPrimerOrientation,
    pub primer_sequence: String,
    pub anneal_sequence: String,
    pub anneal_start_0based: usize,
    pub anneal_end_0based_exclusive: usize,
    pub three_prime_position_0based: usize,
    pub predicted_read_span_start_0based: usize,
    pub predicted_read_span_end_0based_exclusive: usize,
    pub three_prime_distance_bp: usize,
    pub tm_c: f64,
    pub gc_fraction: f64,
    pub anneal_hits: usize,
    pub three_prime_gc_clamp: bool,
    pub longest_homopolymer_run_bp: usize,
    pub self_complementary_run_bp: usize,
    pub reason: String,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
/// Transient machine-readable report for sequencing-primer coverage hints.
pub struct SequencingPrimerOverlayReport {
    pub schema: String,
    pub expected_seq_id: SeqId,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub confirmation_report_id: Option<String>,
    pub min_3prime_anneal_bp: usize,
    pub predicted_read_length_bp: usize,
    #[serde(default)]
    pub primer_seq_ids: Vec<SeqId>,
    pub suggestion_count: usize,
    #[serde(default)]
    pub suggestions: Vec<SequencingPrimerOverlaySuggestion>,
    #[serde(default)]
    pub problem_guidance_count: usize,
    #[serde(default)]
    pub problem_guidance: Vec<SequencingPrimerProblemGuidanceRow>,
    #[serde(default)]
    pub proposal_count: usize,
    #[serde(default)]
    pub proposals: Vec<SequencingPrimerProposalRow>,
}

#[derive(Debug, Clone, Copy, Serialize, Deserialize, PartialEq, Eq, Default)]
#[serde(rename_all = "snake_case")]
pub enum SequenceFeatureRangeRelation {
    #[default]
    Overlap,
    Within,
    Contains,
}

impl SequenceFeatureRangeRelation {
    pub fn as_str(self) -> &'static str {
        match self {
            Self::Overlap => "overlap",
            Self::Within => "within",
            Self::Contains => "contains",
        }
    }
}

#[derive(Debug, Clone, Copy, Serialize, Deserialize, PartialEq, Eq, Default)]
#[serde(rename_all = "snake_case")]
pub enum SequenceFeatureStrandFilter {
    #[default]
    Any,
    Forward,
    Reverse,
}

impl SequenceFeatureStrandFilter {
    pub fn as_str(self) -> &'static str {
        match self {
            Self::Any => "any",
            Self::Forward => "forward",
            Self::Reverse => "reverse",
        }
    }
}

#[derive(Debug, Clone, Copy, Serialize, Deserialize, PartialEq, Eq, Default)]
#[serde(rename_all = "snake_case")]
pub enum SequenceFeatureSortBy {
    FeatureId,
    #[default]
    Start,
    End,
    Kind,
    Length,
}

impl SequenceFeatureSortBy {
    pub fn as_str(self) -> &'static str {
        match self {
            Self::FeatureId => "feature_id",
            Self::Start => "start",
            Self::End => "end",
            Self::Kind => "kind",
            Self::Length => "length",
        }
    }
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
pub struct SequenceFeatureQualifierFilter {
    pub key: String,
    pub value_contains: Option<String>,
    pub value_regex: Option<String>,
    pub case_sensitive: bool,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
pub struct SequenceFeatureQuery {
    pub seq_id: SeqId,
    pub include_source: bool,
    pub include_qualifiers: bool,
    pub kind_in: Vec<String>,
    pub kind_not_in: Vec<String>,
    pub start_0based: Option<usize>,
    pub end_0based_exclusive: Option<usize>,
    pub range_relation: SequenceFeatureRangeRelation,
    pub strand: SequenceFeatureStrandFilter,
    pub label_contains: Option<String>,
    pub label_regex: Option<String>,
    pub qualifier_filters: Vec<SequenceFeatureQualifierFilter>,
    pub min_len_bp: Option<usize>,
    pub max_len_bp: Option<usize>,
    pub limit: Option<usize>,
    pub offset: usize,
    pub sort_by: SequenceFeatureSortBy,
    pub descending: bool,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
pub struct SequenceFeatureQueryRow {
    pub feature_id: usize,
    pub kind: String,
    pub start_0based: usize,
    pub end_0based_exclusive: usize,
    pub length_bp: usize,
    pub strand: String,
    pub label: String,
    pub labels: Vec<String>,
    pub qualifiers: BTreeMap<String, Vec<String>>,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
pub struct SequenceFeatureQueryResult {
    pub schema: String,
    pub seq_id: SeqId,
    pub sequence_length_bp: usize,
    pub total_feature_count: usize,
    pub matched_count: usize,
    pub returned_count: usize,
    pub offset: usize,
    pub limit: usize,
    pub range_relation: String,
    pub strand_filter: String,
    pub sort_by: String,
    pub descending: bool,
    pub query: SequenceFeatureQuery,
    pub rows: Vec<SequenceFeatureQueryRow>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
/// Progress payload emitted by TFBS annotation operations.
pub struct TfbsProgress {
    pub seq_id: String,
    pub motif_id: String,
    pub motif_index: usize,
    pub motif_count: usize,
    pub scanned_steps: usize,
    pub total_steps: usize,
    pub motif_percent: f64,
    pub total_percent: f64,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
/// Progress payload emitted by genome track import operations.
pub struct GenomeTrackImportProgress {
    pub seq_id: String,
    pub source: String,
    pub path: String,
    pub parsed_records: usize,
    pub imported_features: usize,
    pub skipped_records: usize,
    pub done: bool,
}

/// Engine capability snapshot used by adapters for discovery/negotiation.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Capabilities {
    pub protocol_version: String,
    pub supported_operations: Vec<String>,
    pub supported_export_formats: Vec<String>,
    pub deterministic_operation_log: bool,
}

fn default_true() -> bool {
    true
}

fn default_poly_t_prefix_min_bp() -> usize {
    18
}

fn default_rna_seed_stride_bp() -> usize {
    1
}

fn default_min_weighted_seed_hit_fraction() -> f64 {
    0.05
}

fn default_min_unique_matched_kmers() -> usize {
    12
}

fn default_max_median_transcript_gap() -> f64 {
    4.0
}

fn default_min_chain_consistency_fraction() -> f64 {
    0.40
}

fn default_min_confirmed_exon_transitions() -> usize {
    1
}

fn default_min_transition_support_fraction() -> f64 {
    0.05
}

fn default_rna_read_checkpoint_every_reads() -> usize {
    10_000
}

/// Composite seed-gate thresholds reused by RNA-read interpretation reports,
/// progress payloads, and adapter-side inspection tools.
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
#[serde(default)]
pub struct RnaReadSeedFilterConfig {
    pub kmer_len: usize,
    #[serde(default = "default_rna_seed_stride_bp")]
    pub seed_stride_bp: usize,
    pub min_seed_hit_fraction: f64,
    #[serde(default = "default_min_weighted_seed_hit_fraction")]
    pub min_weighted_seed_hit_fraction: f64,
    #[serde(default = "default_min_unique_matched_kmers")]
    pub min_unique_matched_kmers: usize,
    #[serde(default = "default_max_median_transcript_gap")]
    pub max_median_transcript_gap: f64,
    #[serde(default = "default_min_chain_consistency_fraction")]
    pub min_chain_consistency_fraction: f64,
    #[serde(default = "default_min_confirmed_exon_transitions")]
    pub min_confirmed_exon_transitions: usize,
    #[serde(default = "default_min_transition_support_fraction")]
    pub min_transition_support_fraction: f64,
    #[serde(default = "default_true")]
    pub cdna_poly_t_flip_enabled: bool,
    #[serde(default = "default_poly_t_prefix_min_bp")]
    pub poly_t_prefix_min_bp: usize,
}

impl Default for RnaReadSeedFilterConfig {
    fn default() -> Self {
        Self {
            kmer_len: 10,
            seed_stride_bp: 1,
            min_seed_hit_fraction: 0.30,
            min_weighted_seed_hit_fraction: 0.05,
            min_unique_matched_kmers: 12,
            max_median_transcript_gap: 4.0,
            min_chain_consistency_fraction: 0.40,
            min_confirmed_exon_transitions: 1,
            min_transition_support_fraction: 0.05,
            cdna_poly_t_flip_enabled: true,
            poly_t_prefix_min_bp: 18,
        }
    }
}

/// Pairwise phase-2 alignment parameters shared by RNA-read mapping reports
/// and inspection/export adapters.
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
#[serde(default)]
pub struct RnaReadAlignConfig {
    pub band_width_bp: usize,
    pub min_identity_fraction: f64,
    pub max_secondary_mappings: usize,
}

impl Default for RnaReadAlignConfig {
    fn default() -> Self {
        Self {
            band_width_bp: 24,
            min_identity_fraction: 0.60,
            max_secondary_mappings: 3,
        }
    }
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
pub struct RnaReadMappingHit {
    #[serde(default)]
    pub alignment_mode: RnaReadAlignmentMode,
    pub transcript_feature_id: usize,
    pub transcript_id: String,
    pub transcript_label: String,
    pub strand: String,
    pub query_start_0based: usize,
    pub query_end_0based_exclusive: usize,
    #[serde(default)]
    pub query_reverse_complemented: bool,
    pub target_start_1based: usize,
    pub target_end_1based: usize,
    #[serde(default)]
    pub target_start_offset_0based: usize,
    #[serde(default)]
    pub target_end_offset_0based_exclusive: usize,
    pub matches: usize,
    pub mismatches: usize,
    pub score: isize,
    pub identity_fraction: f64,
    pub query_coverage_fraction: f64,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
pub struct RnaReadAlignmentDisplay {
    pub transcript_feature_id: usize,
    pub transcript_id: String,
    pub transcript_label: String,
    pub strand: String,
    #[serde(default)]
    pub alignment_mode: RnaReadAlignmentMode,
    #[serde(default)]
    pub query_reverse_complemented: bool,
    #[serde(default)]
    pub query_start_0based: usize,
    #[serde(default)]
    pub query_end_0based_exclusive: usize,
    #[serde(default)]
    pub target_start_1based: usize,
    #[serde(default)]
    pub target_end_1based: usize,
    #[serde(default)]
    pub target_start_offset_0based: usize,
    #[serde(default)]
    pub target_end_offset_0based_exclusive: usize,
    #[serde(default)]
    pub target_length_bp: usize,
    pub score: isize,
    pub identity_fraction: f64,
    pub query_coverage_fraction: f64,
    #[serde(default)]
    pub target_coverage_fraction: f64,
    #[serde(default)]
    pub matches: usize,
    #[serde(default)]
    pub mismatches: usize,
    #[serde(default)]
    pub insertions: usize,
    #[serde(default)]
    pub deletions: usize,
    #[serde(default)]
    pub aligned_columns: usize,
    #[serde(default)]
    pub aligned_query: String,
    #[serde(default)]
    pub aligned_midline: String,
    #[serde(default)]
    pub aligned_target: String,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
pub struct RnaReadStrandAssignmentDiagnostics {
    pub selected_strand: String,
    pub selected_reason: String,
    pub selected_transition_hits: usize,
    pub selected_exon_hits: usize,
    pub plus_best_transcript_id: String,
    pub plus_best_transition_hits: usize,
    pub plus_best_exon_hits: usize,
    pub minus_best_transcript_id: String,
    pub minus_best_transition_hits: usize,
    pub minus_best_exon_hits: usize,
    pub competing_opposite_strand: bool,
    pub ambiguous_near_tie: bool,
    pub chain_preferred_strand: String,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
pub struct RnaReadOriginCandidateContribution {
    pub candidate_role: String,
    pub transcript_id: String,
    pub strand: String,
    pub transition_hits: usize,
    pub exon_hits: usize,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
pub struct RnaReadInterpretationHit {
    pub record_index: usize,
    pub source_byte_offset: usize,
    pub header_id: String,
    pub sequence: String,
    pub read_length_bp: usize,
    pub tested_kmers: usize,
    pub matched_kmers: usize,
    pub seed_hit_fraction: f64,
    #[serde(default)]
    pub weighted_seed_hit_fraction: f64,
    #[serde(default)]
    pub weighted_matched_kmers: f64,
    #[serde(default)]
    pub seed_chain_transcript_id: String,
    #[serde(default)]
    pub seed_chain_support_kmers: usize,
    #[serde(default)]
    pub seed_chain_support_fraction: f64,
    #[serde(default)]
    pub seed_median_transcript_gap: f64,
    #[serde(default)]
    pub seed_transcript_gap_count: usize,
    #[serde(default)]
    pub exon_path_transcript_id: String,
    #[serde(default)]
    pub exon_path: String,
    #[serde(default)]
    pub exon_transitions_confirmed: usize,
    #[serde(default)]
    pub exon_transitions_total: usize,
    #[serde(default)]
    pub reverse_complement_applied: bool,
    #[serde(default)]
    pub strand_diagnostics: RnaReadStrandAssignmentDiagnostics,
    #[serde(default)]
    pub origin_class: RnaReadOriginClass,
    #[serde(default)]
    pub origin_reason: String,
    #[serde(default)]
    pub origin_confidence: f64,
    #[serde(default)]
    pub strand_confidence: f64,
    #[serde(default)]
    pub origin_candidates: Vec<RnaReadOriginCandidateContribution>,
    pub perfect_seed_match: bool,
    pub passed_seed_filter: bool,
    #[serde(default)]
    pub msa_eligible: bool,
    #[serde(default)]
    pub msa_eligibility_reason: String,
    pub best_mapping: Option<RnaReadMappingHit>,
    pub secondary_mappings: Vec<RnaReadMappingHit>,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
pub struct RnaReadExonSupportFrequency {
    pub start_1based: usize,
    pub end_1based: usize,
    pub support_read_count: usize,
    pub support_fraction: f64,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
pub struct RnaReadJunctionSupportFrequency {
    pub donor_1based: usize,
    pub acceptor_1based: usize,
    pub support_read_count: usize,
    pub support_fraction: f64,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
pub struct RnaReadTransitionSupportRow {
    pub from_exon_ordinal: usize,
    pub to_exon_ordinal: usize,
    pub from_start_1based: usize,
    pub from_end_1based: usize,
    pub to_start_1based: usize,
    pub to_end_1based: usize,
    pub support_read_count: usize,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
pub struct RnaReadIsoformSupportRow {
    pub transcript_feature_id: usize,
    pub transcript_id: String,
    pub transcript_label: String,
    pub strand: String,
    pub exon_count: usize,
    pub expected_transition_count: usize,
    pub reads_assigned: usize,
    pub reads_seed_passed: usize,
    pub transition_rows_supported: usize,
    pub transition_rows_supported_fraction: f64,
    pub mean_seed_median_gap: f64,
    pub mean_confirmed_transition_fraction: f64,
    pub best_seed_hit_fraction: f64,
    pub best_weighted_seed_hit_fraction: f64,
    #[serde(default)]
    pub reads_chain_same_strand: usize,
    #[serde(default)]
    pub reads_with_opposite_strand_competition: usize,
    #[serde(default)]
    pub reads_ambiguous_strand_ties: usize,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
pub struct RnaReadMappedIsoformSupportRow {
    pub transcript_feature_id: usize,
    pub transcript_id: String,
    pub transcript_label: String,
    pub strand: String,
    #[serde(default)]
    pub aligned_read_count: usize,
    #[serde(default)]
    pub msa_eligible_read_count: usize,
    #[serde(default)]
    pub mean_identity_fraction: f64,
    #[serde(default)]
    pub mean_query_coverage_fraction: f64,
    #[serde(default)]
    pub best_alignment_score: isize,
    #[serde(default)]
    pub secondary_mapping_total: usize,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
pub struct RnaReadSampleSheetExport {
    pub schema: String,
    pub path: String,
    pub report_count: usize,
    pub appended: bool,
    #[serde(default)]
    pub gene_ids: Vec<String>,
    #[serde(default)]
    pub complete_rule: RnaReadGeneSupportCompleteRule,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
pub struct RnaReadExonPathsExport {
    pub schema: String,
    pub path: String,
    pub report_id: String,
    pub selection: RnaReadHitSelection,
    pub row_count: usize,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
pub struct RnaReadExonAbundanceExport {
    pub schema: String,
    pub path: String,
    pub report_id: String,
    pub selection: RnaReadHitSelection,
    pub selected_read_count: usize,
    pub exon_row_count: usize,
    pub transition_row_count: usize,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
pub struct RnaReadGeneExonSupportRow {
    pub gene_id: String,
    pub exon_ordinal: usize,
    pub start_1based: usize,
    pub end_1based: usize,
    pub support_read_count: usize,
    pub support_fraction: f64,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
pub struct RnaReadGeneExonPairSupportRow {
    pub gene_id: String,
    pub from_exon_ordinal: usize,
    pub from_start_1based: usize,
    pub from_end_1based: usize,
    pub to_exon_ordinal: usize,
    pub to_start_1based: usize,
    pub to_end_1based: usize,
    pub support_read_count: usize,
    pub support_fraction: f64,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
pub struct RnaReadGeneSupportCohortSummary {
    pub read_count: usize,
    pub exon_support: Vec<RnaReadGeneExonSupportRow>,
    pub exon_pair_support: Vec<RnaReadGeneExonPairSupportRow>,
    pub direct_transition_support: Vec<RnaReadGeneExonPairSupportRow>,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
pub struct RnaReadGeneSupportSummary {
    pub schema: String,
    pub report_id: String,
    pub seq_id: String,
    pub requested_gene_ids: Vec<String>,
    pub matched_gene_ids: Vec<String>,
    pub missing_gene_ids: Vec<String>,
    pub selected_record_indices: Vec<usize>,
    #[serde(default)]
    pub complete_rule: RnaReadGeneSupportCompleteRule,
    pub aligned_base_count: usize,
    pub accepted_target_count: usize,
    pub fragment_count: usize,
    pub complete_count: usize,
    pub complete_strict_count: usize,
    pub complete_exact_count: usize,
    pub all_target: RnaReadGeneSupportCohortSummary,
    pub fragments: RnaReadGeneSupportCohortSummary,
    pub complete: RnaReadGeneSupportCohortSummary,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
pub struct RnaReadGeneSupportAuditPair {
    pub from_exon_ordinal: usize,
    pub from_start_1based: usize,
    pub from_end_1based: usize,
    pub to_exon_ordinal: usize,
    pub to_start_1based: usize,
    pub to_end_1based: usize,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
pub struct RnaReadGeneSupportAuditRow {
    pub record_index: usize,
    pub header_id: String,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub gene_id: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub transcript_feature_id: Option<usize>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub transcript_id: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub transcript_label: Option<String>,
    #[serde(default)]
    pub status: RnaReadGeneSupportAuditStatus,
    pub status_reason: String,
    #[serde(default)]
    pub full_length_exact: bool,
    #[serde(default)]
    pub full_length_near: bool,
    #[serde(default)]
    pub full_length_strict: bool,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub full_length_class: Option<String>,
    #[serde(default)]
    pub mapped_exon_ordinals: Vec<usize>,
    #[serde(default)]
    pub exon_pairs: Vec<RnaReadGeneSupportAuditPair>,
    #[serde(default)]
    pub direct_transition_pairs: Vec<RnaReadGeneSupportAuditPair>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub score: Option<isize>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub identity_fraction: Option<f64>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub query_coverage_fraction: Option<f64>,
    #[serde(default)]
    pub passed_seed_filter: bool,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
pub struct RnaReadGeneSupportAudit {
    pub schema: String,
    pub report_id: String,
    pub seq_id: String,
    pub requested_gene_ids: Vec<String>,
    pub matched_gene_ids: Vec<String>,
    pub missing_gene_ids: Vec<String>,
    pub selected_record_indices: Vec<usize>,
    #[serde(default)]
    pub complete_rule: RnaReadGeneSupportCompleteRule,
    #[serde(default)]
    pub cohort_filter: RnaReadGeneSupportAuditCohortFilter,
    pub evaluated_row_count: usize,
    pub row_count: usize,
    pub accepted_target_record_indices: Vec<usize>,
    pub fragment_record_indices: Vec<usize>,
    pub complete_record_indices: Vec<usize>,
    pub complete_strict_record_indices: Vec<usize>,
    pub complete_exact_record_indices: Vec<usize>,
    pub rows: Vec<RnaReadGeneSupportAuditRow>,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
pub struct RnaReadScoreDensitySvgExport {
    pub schema: String,
    pub path: String,
    pub report_id: String,
    pub scale: RnaReadScoreDensityScale,
    #[serde(default)]
    pub variant: RnaReadScoreDensityVariant,
    pub bin_count: usize,
    pub max_bin_count: u64,
    pub total_scored_reads: u64,
    pub derived_from_report_hits_only: bool,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
pub struct RnaReadAlignmentDotplotSvgExport {
    pub schema: String,
    pub path: String,
    pub report_id: String,
    pub selection: RnaReadHitSelection,
    pub point_count: usize,
    pub rendered_point_count: usize,
    pub max_points: usize,
    pub min_score: isize,
    pub max_score: isize,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
pub struct RnaReadAlignmentTsvExport {
    pub schema: String,
    pub path: String,
    pub report_id: String,
    pub selection: RnaReadHitSelection,
    pub row_count: usize,
    pub aligned_count: usize,
    pub limit: Option<usize>,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
pub struct RnaReadAlignmentInspectionSubsetSpec {
    pub effect_filter: RnaReadAlignmentInspectionEffectFilter,
    pub sort_key: RnaReadAlignmentInspectionSortKey,
    pub search: String,
    pub selected_record_indices: Vec<usize>,
    #[serde(default)]
    pub score_density_variant: RnaReadScoreDensityVariant,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub score_density_seed_filter_override: Option<RnaReadSeedFilterConfig>,
    pub score_bin_index: Option<usize>,
    pub score_bin_count: usize,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
pub struct RnaReadMappedSupportExonAttribution {
    pub start_1based: usize,
    pub end_1based: usize,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
pub struct RnaReadMappedSupportJunctionAttribution {
    pub donor_1based: usize,
    pub acceptor_1based: usize,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
/// Portable per-read alignment-inspection row used by GUI/CLI report tables.
pub struct RnaReadAlignmentInspectionRow {
    pub rank: usize,
    pub record_index: usize,
    pub header_id: String,
    #[serde(default)]
    pub phase1_primary_transcript_id: String,
    #[serde(default)]
    pub seed_chain_transcript_id: String,
    #[serde(default)]
    pub exon_path_transcript_id: String,
    #[serde(default)]
    pub exon_path: String,
    #[serde(default)]
    pub exon_transitions_confirmed: usize,
    #[serde(default)]
    pub exon_transitions_total: usize,
    #[serde(default)]
    pub selected_strand: String,
    #[serde(default)]
    pub reverse_complement_applied: bool,
    #[serde(default)]
    pub alignment_effect: RnaReadAlignmentEffect,
    pub transcript_id: String,
    pub transcript_label: String,
    pub strand: String,
    pub alignment_mode: RnaReadAlignmentMode,
    #[serde(default)]
    pub target_start_1based: usize,
    #[serde(default)]
    pub target_end_1based: usize,
    #[serde(default)]
    pub target_length_bp: usize,
    pub score: isize,
    pub identity_fraction: f64,
    pub query_coverage_fraction: f64,
    #[serde(default)]
    pub target_coverage_fraction: f64,
    #[serde(default)]
    pub full_length_exact: bool,
    #[serde(default)]
    pub full_length_near: bool,
    #[serde(default)]
    pub full_length_strict: bool,
    #[serde(default)]
    pub secondary_mapping_count: usize,
    pub seed_hit_fraction: f64,
    pub weighted_seed_hit_fraction: f64,
    pub passed_seed_filter: bool,
    pub msa_eligible: bool,
    pub origin_class: RnaReadOriginClass,
    #[serde(default)]
    pub mapped_exon_support: Vec<RnaReadMappedSupportExonAttribution>,
    #[serde(default)]
    pub mapped_junction_support: Vec<RnaReadMappedSupportJunctionAttribution>,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
pub struct RnaReadAlignmentInspection {
    pub schema: String,
    pub report_id: String,
    pub seq_id: String,
    pub selection: RnaReadHitSelection,
    pub row_count: usize,
    pub aligned_count: usize,
    pub subset_match_count: usize,
    pub limit: usize,
    pub subset_spec: RnaReadAlignmentInspectionSubsetSpec,
    pub align_min_identity_fraction: f64,
    pub max_secondary_mappings: usize,
    pub rows: Vec<RnaReadAlignmentInspectionRow>,
}

#[derive(Debug, Clone, Copy, Serialize, Deserialize, PartialEq, Eq, Default)]
#[serde(rename_all = "snake_case")]
/// How one transcript/CDS protein derivation resolved its genetic code table.
pub enum TranscriptProteinTranslationTableSource {
    #[default]
    StandardDefault,
    ExplicitCdsQualifier,
    ExplicitTranscriptQualifier,
    ExplicitSourceQualifier,
    OrganellePlastidDefault,
    AmbiguousMitochondrialDefault,
}

impl TranscriptProteinTranslationTableSource {
    pub fn as_str(self) -> &'static str {
        match self {
            Self::StandardDefault => "standard_default",
            Self::ExplicitCdsQualifier => "explicit_cds_qualifier",
            Self::ExplicitTranscriptQualifier => "explicit_transcript_qualifier",
            Self::ExplicitSourceQualifier => "explicit_source_qualifier",
            Self::OrganellePlastidDefault => "organelle_plastid_default",
            Self::AmbiguousMitochondrialDefault => "ambiguous_mitochondrial_default",
        }
    }
}

#[derive(Debug, Clone, Copy, Serialize, Deserialize, PartialEq, Eq, Default)]
#[serde(rename_all = "snake_case")]
/// How a transcript-to-protein derivation chose its coding span.
pub enum TranscriptProteinDerivationMode {
    #[default]
    AnnotatedCds,
    InferredOrf,
    HeuristicLongestFrame,
}

impl TranscriptProteinDerivationMode {
    pub fn as_str(self) -> &'static str {
        match self {
            Self::AnnotatedCds => "annotated_cds",
            Self::InferredOrf => "inferred_orf",
            Self::HeuristicLongestFrame => "heuristic_longest_frame",
        }
    }
}

#[derive(Debug, Clone, Copy, Serialize, Deserialize, PartialEq, Eq)]
#[serde(rename_all = "snake_case")]
/// Named codon-bias profile used for protein back-translation.
pub enum TranslationSpeedProfile {
    Human,
    Mouse,
    Yeast,
    Ecoli,
}

impl TranslationSpeedProfile {
    pub fn as_str(self) -> &'static str {
        match self {
            Self::Human => "human",
            Self::Mouse => "mouse",
            Self::Yeast => "yeast",
            Self::Ecoli => "ecoli",
        }
    }
}

#[derive(Debug, Clone, Copy, Serialize, Deserialize, PartialEq, Eq)]
#[serde(rename_all = "snake_case")]
/// Qualitative codon-speed bias for reverse translation.
pub enum TranslationSpeedMark {
    Fast,
    Slow,
}

impl TranslationSpeedMark {
    pub fn as_str(self) -> &'static str {
        match self {
            Self::Fast => "fast",
            Self::Slow => "slow",
        }
    }
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
/// Portable transcript/CDS-to-protein derivation summary.
///
/// This record is intentionally narrower than a first-class protein sequence
/// window. It captures the resolved coding span, genetic code selection, and
/// translated amino-acid sequence so GUI/CLI/shell code can inspect the same
/// deterministic transcript-translation decision without re-implementing the
/// biology locally.
pub struct TranscriptProteinDerivation {
    pub transcript_id: String,
    pub transcript_label: String,
    pub source_seq_id: SeqId,
    pub source_feature_id: usize,
    #[serde(default)]
    pub derivation_mode: TranscriptProteinDerivationMode,
    #[serde(default)]
    pub cds_ranges_1based: Vec<(usize, usize)>,
    pub cds_length_bp: usize,
    pub protein_sequence: String,
    pub protein_length_aa: usize,
    pub translation_table: usize,
    pub translation_table_label: String,
    #[serde(default)]
    pub translation_table_source: TranscriptProteinTranslationTableSource,
    pub codon_start: usize,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub organism: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub organelle: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub translation_speed_profile_hint: Option<String>,
    #[serde(default)]
    pub terminal_stop_trimmed: bool,
    #[serde(default)]
    pub warnings: Vec<String>,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
pub struct RnaReadPairwiseAlignmentDetail {
    pub schema: String,
    pub report_id: String,
    pub seq_id: String,
    pub record_index: usize,
    pub header_id: String,
    pub transcript_id: String,
    pub transcript_label: String,
    pub strand: String,
    #[serde(default)]
    pub alignment_mode: RnaReadAlignmentMode,
    #[serde(default)]
    pub backend: RnaReadAlignmentBackend,
    pub query_length_bp: usize,
    pub target_length_bp: usize,
    pub aligned_query_start_0based: usize,
    pub aligned_query_end_0based_exclusive: usize,
    pub aligned_target_start_offset_0based: usize,
    pub aligned_target_end_offset_0based_exclusive: usize,
    pub target_start_1based: usize,
    pub target_end_1based: usize,
    pub aligned_columns: usize,
    pub matches: usize,
    pub mismatches: usize,
    pub insertions: usize,
    pub deletions: usize,
    pub score: isize,
    pub identity_fraction: f64,
    pub query_coverage_fraction: f64,
    #[serde(default)]
    pub target_coverage_fraction: f64,
    pub cigar: String,
    pub aligned_query: String,
    pub aligned_relation: String,
    pub aligned_target: String,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
pub struct RnaSeedHashCatalogEntry {
    pub seed_bits: u32,
    pub kmer_sequence: String,
    pub transcript_feature_id: usize,
    pub transcript_id: String,
    pub transcript_label: String,
    pub strand: String,
    pub template_offset_0based: usize,
    pub genomic_pos_1based: usize,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
pub struct RnaSeedHashTemplateAuditEntry {
    pub transcript_feature_id: usize,
    pub transcript_id: String,
    pub transcript_label: String,
    pub strand: String,
    pub template_sequence: String,
    pub template_length_bp: usize,
    pub template_first_genomic_pos_1based: usize,
    pub template_last_genomic_pos_1based: usize,
    pub reverse_complemented_from_genome: bool,
}

/// One genome-position bin used for running RNA-read seed-confirmation
/// statistics.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct RnaReadSeedHistogramBin {
    pub start_1based: usize,
    pub end_1based: usize,
    pub confirmed_plus: u64,
    pub confirmed_minus: u64,
}

/// Lightweight top-hit row included in running RNA-read progress updates.
#[derive(Debug, Clone, Serialize, Deserialize, Default)]
pub struct RnaReadTopHitPreview {
    pub record_index: usize,
    pub header_id: String,
    pub seed_hit_fraction: f64,
    pub weighted_seed_hit_fraction: f64,
    #[serde(default)]
    pub weighted_matched_kmers: f64,
    #[serde(default)]
    pub seed_chain_transcript_id: String,
    #[serde(default)]
    pub seed_chain_support_kmers: usize,
    #[serde(default)]
    pub seed_chain_support_fraction: f64,
    #[serde(default)]
    pub seed_median_transcript_gap: f64,
    #[serde(default)]
    pub seed_transcript_gap_count: usize,
    pub matched_kmers: usize,
    pub tested_kmers: usize,
    pub passed_seed_filter: bool,
    #[serde(default)]
    pub reverse_complement_applied: bool,
    #[serde(default)]
    pub selected_strand: String,
    #[serde(default)]
    pub competing_opposite_strand: bool,
    #[serde(default)]
    pub ambiguous_strand_tie: bool,
    #[serde(default)]
    pub origin_class: RnaReadOriginClass,
    #[serde(default)]
    pub origin_reason: String,
    #[serde(default)]
    pub origin_confidence: f64,
    #[serde(default)]
    pub strand_confidence: f64,
    #[serde(default)]
    pub origin_candidates: Vec<RnaReadOriginCandidateContribution>,
    #[serde(default)]
    pub msa_eligible: bool,
    #[serde(default)]
    pub msa_eligibility_reason: String,
    #[serde(default)]
    pub aligned: bool,
    #[serde(default)]
    pub best_alignment_mode: String,
    #[serde(default)]
    pub best_alignment_transcript_id: String,
    #[serde(default)]
    pub best_alignment_transcript_label: String,
    #[serde(default)]
    pub best_alignment_strand: String,
    #[serde(default)]
    pub best_alignment_target_start_1based: usize,
    #[serde(default)]
    pub best_alignment_target_end_1based: usize,
    #[serde(default)]
    pub best_alignment_identity_fraction: f64,
    #[serde(default)]
    pub best_alignment_query_coverage_fraction: f64,
    #[serde(default)]
    pub best_alignment_score: isize,
    #[serde(default)]
    pub secondary_mapping_count: usize,
    pub read_length_bp: usize,
    pub sequence: String,
    pub sequence_preview: String,
}

/// Progress payload emitted by RNA-read interpretation operations.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct RnaReadInterpretProgress {
    pub seq_id: String,
    pub reads_processed: usize,
    pub reads_total: usize,
    #[serde(default)]
    pub read_bases_processed: u64,
    #[serde(default)]
    pub mean_read_length_bp: f64,
    #[serde(default)]
    pub median_read_length_bp: usize,
    #[serde(default)]
    pub p95_read_length_bp: usize,
    #[serde(default)]
    pub input_bytes_processed: u64,
    #[serde(default)]
    pub input_bytes_total: u64,
    pub seed_passed: usize,
    pub aligned: usize,
    pub tested_kmers: usize,
    pub matched_kmers: usize,
    #[serde(default)]
    pub seed_compute_ms: f64,
    #[serde(default)]
    pub align_compute_ms: f64,
    #[serde(default)]
    pub io_read_ms: f64,
    #[serde(default)]
    pub fasta_parse_ms: f64,
    #[serde(default)]
    pub normalize_compute_ms: f64,
    #[serde(default)]
    pub inference_compute_ms: f64,
    #[serde(default)]
    pub progress_emit_ms: f64,
    pub update_every_reads: usize,
    pub done: bool,
    pub bins: Vec<RnaReadSeedHistogramBin>,
    pub score_density_bins: Vec<u64>,
    #[serde(default)]
    pub seed_pass_score_density_bins: Vec<u64>,
    #[serde(default)]
    pub top_hits_preview: Vec<RnaReadTopHitPreview>,
    #[serde(default)]
    pub transition_support_rows: Vec<RnaReadTransitionSupportRow>,
    #[serde(default)]
    pub isoform_support_rows: Vec<RnaReadIsoformSupportRow>,
    #[serde(default)]
    pub mapped_exon_support_frequencies: Vec<RnaReadExonSupportFrequency>,
    #[serde(default)]
    pub mapped_junction_support_frequencies: Vec<RnaReadJunctionSupportFrequency>,
    #[serde(default)]
    pub mapped_isoform_support_rows: Vec<RnaReadMappedIsoformSupportRow>,
    #[serde(default)]
    pub reads_with_transition_support: usize,
    #[serde(default)]
    pub transition_confirmations: usize,
    #[serde(default)]
    pub junction_crossing_seed_bits_indexed: usize,
    #[serde(default)]
    pub origin_class_counts: BTreeMap<String, usize>,
}

/// Persisted RNA-read interpretation report shared across GUI, CLI, and shell
/// inspection/export flows.
#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
pub struct RnaReadInterpretationReport {
    pub schema: String,
    pub report_id: String,
    #[serde(default)]
    pub report_mode: RnaReadReportMode,
    pub seq_id: String,
    pub seed_feature_id: usize,
    pub generated_at_unix_ms: u128,
    pub profile: RnaReadInterpretationProfile,
    pub input_path: String,
    pub input_format: RnaReadInputFormat,
    pub scope: SplicingScopePreset,
    #[serde(default)]
    pub origin_mode: RnaReadOriginMode,
    #[serde(default)]
    pub target_gene_ids: Vec<String>,
    #[serde(default)]
    pub roi_seed_capture_enabled: bool,
    #[serde(default)]
    pub checkpoint_path: Option<String>,
    #[serde(default = "default_rna_read_checkpoint_every_reads")]
    pub checkpoint_every_reads: usize,
    #[serde(default)]
    pub resumed_from_checkpoint: bool,
    pub seed_filter: RnaReadSeedFilterConfig,
    pub align_config: RnaReadAlignConfig,
    pub read_count_total: usize,
    pub read_count_seed_passed: usize,
    pub read_count_aligned: usize,
    #[serde(default)]
    pub retained_count_msa_eligible: usize,
    #[serde(default)]
    pub warnings: Vec<String>,
    #[serde(default)]
    pub hits: Vec<RnaReadInterpretationHit>,
    #[serde(default)]
    pub exon_support_frequencies: Vec<RnaReadExonSupportFrequency>,
    #[serde(default)]
    pub junction_support_frequencies: Vec<RnaReadJunctionSupportFrequency>,
    #[serde(default)]
    pub transition_support_rows: Vec<RnaReadTransitionSupportRow>,
    #[serde(default)]
    pub isoform_support_rows: Vec<RnaReadIsoformSupportRow>,
    #[serde(default)]
    pub mapped_isoform_support_rows: Vec<RnaReadMappedIsoformSupportRow>,
    #[serde(default)]
    pub origin_class_counts: BTreeMap<String, usize>,
    #[serde(default)]
    pub read_length_counts_all: Vec<u64>,
    #[serde(default)]
    pub read_length_counts_seed_passed: Vec<u64>,
    #[serde(default)]
    pub read_length_counts_aligned: Vec<u64>,
    #[serde(default)]
    pub read_length_counts_full_length_exact: Vec<u64>,
    #[serde(default)]
    pub read_length_counts_full_length_near: Vec<u64>,
    #[serde(default)]
    pub read_length_counts_full_length_strict: Vec<u64>,
    #[serde(default)]
    pub score_density_bins: Vec<u64>,
    #[serde(default)]
    pub seed_pass_score_density_bins: Vec<u64>,
}

/// Lightweight listing row for RNA-read reports, kept portable so every
/// adapter can render the same report inventory.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct RnaReadInterpretationReportSummary {
    pub report_id: String,
    #[serde(default)]
    pub report_mode: RnaReadReportMode,
    pub seq_id: String,
    pub generated_at_unix_ms: u128,
    pub profile: RnaReadInterpretationProfile,
    pub input_path: String,
    pub input_format: RnaReadInputFormat,
    pub seed_feature_id: usize,
    pub scope: SplicingScopePreset,
    #[serde(default)]
    pub origin_mode: RnaReadOriginMode,
    #[serde(default)]
    pub target_gene_count: usize,
    #[serde(default)]
    pub roi_seed_capture_enabled: bool,
    pub read_count_total: usize,
    pub read_count_seed_passed: usize,
    pub read_count_aligned: usize,
    #[serde(default)]
    pub retained_count_msa_eligible: usize,
}
