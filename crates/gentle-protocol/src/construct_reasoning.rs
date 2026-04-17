//! Portable construct-reasoning graph records.
//!
//! These types describe inspectable, editable design reasoning that can be
//! shared across GUI/CLI/MCP/JS/Lua adapters without coupling to one frontend.

use crate::{SeqId, TranslationSpeedMark, TranslationSpeedProfile, TranslationSpeedProfileSource};
use serde::{Deserialize, Serialize};
use std::collections::BTreeMap;

pub const CONSTRUCT_OBJECTIVE_SCHEMA: &str = "gentle.construct_objective.v1";
pub const DESIGN_EVIDENCE_SCHEMA: &str = "gentle.design_evidence.v1";
pub const DESIGN_FACT_SCHEMA: &str = "gentle.design_fact.v1";
pub const DESIGN_DECISION_NODE_SCHEMA: &str = "gentle.design_decision_node.v1";
pub const CONSTRUCT_CANDIDATE_SCHEMA: &str = "gentle.construct_candidate.v1";
pub const ANNOTATION_CANDIDATE_SCHEMA: &str = "gentle.annotation_candidate.v1";
pub const ANNOTATION_CANDIDATE_WRITEBACK_SCHEMA: &str = "gentle.annotation_candidate_writeback.v1";
pub const CONSTRUCT_REASONING_GRAPH_SCHEMA: &str = "gentle.construct_reasoning_graph.v1";
pub const CONSTRUCT_REASONING_STORE_SCHEMA: &str = "gentle.construct_reasoning_store.v1";
pub const HOST_PROFILE_CATALOG_SCHEMA: &str = "gentle.host_profile_catalog.v1";

fn default_construct_objective_schema() -> String {
    CONSTRUCT_OBJECTIVE_SCHEMA.to_string()
}

fn default_design_evidence_schema() -> String {
    DESIGN_EVIDENCE_SCHEMA.to_string()
}

fn default_design_fact_schema() -> String {
    DESIGN_FACT_SCHEMA.to_string()
}

fn default_design_decision_node_schema() -> String {
    DESIGN_DECISION_NODE_SCHEMA.to_string()
}

fn default_construct_candidate_schema() -> String {
    CONSTRUCT_CANDIDATE_SCHEMA.to_string()
}

fn default_annotation_candidate_schema() -> String {
    ANNOTATION_CANDIDATE_SCHEMA.to_string()
}

fn default_annotation_candidate_writeback_schema() -> String {
    ANNOTATION_CANDIDATE_WRITEBACK_SCHEMA.to_string()
}

fn default_construct_reasoning_graph_schema() -> String {
    CONSTRUCT_REASONING_GRAPH_SCHEMA.to_string()
}

fn default_construct_reasoning_store_schema() -> String {
    CONSTRUCT_REASONING_STORE_SCHEMA.to_string()
}

#[derive(Debug, Clone, Copy, Serialize, Deserialize, PartialEq, Eq, PartialOrd, Ord, Default)]
#[serde(rename_all = "snake_case")]
/// Ranking intent for protein-to-DNA handoff candidates.
pub enum ProteinToDnaHandoffRankingGoal {
    #[default]
    BalancedProvenance,
    NativeFidelity,
    ExpressionOptimized,
}

impl ProteinToDnaHandoffRankingGoal {
    pub fn as_str(self) -> &'static str {
        match self {
            Self::BalancedProvenance => "balanced_provenance",
            Self::NativeFidelity => "native_fidelity",
            Self::ExpressionOptimized => "expression_optimized",
        }
    }
}

#[derive(Debug, Clone, Copy, Serialize, Deserialize, PartialEq, Eq, PartialOrd, Ord, Default)]
#[serde(rename_all = "snake_case")]
/// Candidate family for a protein-to-DNA handoff.
pub enum ProteinToDnaHandoffStrategy {
    #[default]
    TranscriptNativeReuse,
    FeatureCodingDna,
    ReverseTranslatedSynthetic,
}

impl ProteinToDnaHandoffStrategy {
    pub fn as_str(self) -> &'static str {
        match self {
            Self::TranscriptNativeReuse => "transcript_native_reuse",
            Self::FeatureCodingDna => "feature_coding_dna",
            Self::ReverseTranslatedSynthetic => "reverse_translated_synthetic",
        }
    }
}

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Default)]
#[serde(default, deny_unknown_fields)]
/// Summary of how much of the requested protein evidence one handoff candidate preserves.
pub struct ProteinToDnaHandoffCoverage {
    pub amino_acid_start_1based: Option<usize>,
    pub amino_acid_end_1based: Option<usize>,
    pub covered_amino_acids: usize,
    pub requested_amino_acids: usize,
    pub preserved_feature_labels: Vec<String>,
    pub relaxed_feature_labels: Vec<String>,
}

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Default)]
#[serde(default, deny_unknown_fields)]
/// Typed candidate detail for construct-graph-backed protein-to-DNA handoff reasoning.
pub struct ProteinToDnaHandoffCandidate {
    pub strategy: ProteinToDnaHandoffStrategy,
    pub ranking_goal: ProteinToDnaHandoffRankingGoal,
    pub source_protein_seq_id: SeqId,
    pub source_artifact_refs: Vec<String>,
    pub transcript_id: Option<String>,
    pub projection_id: Option<String>,
    pub ensembl_entry_id: Option<String>,
    pub feature_query: Option<String>,
    pub matched_feature_label: Option<String>,
    pub coverage: ProteinToDnaHandoffCoverage,
    pub preserved_constraints: Vec<String>,
    pub relaxed_constraints: Vec<String>,
    pub translation_table: Option<usize>,
    pub speed_profile: Option<TranslationSpeedProfile>,
    pub speed_profile_source: Option<TranslationSpeedProfileSource>,
    pub speed_mark: Option<TranslationSpeedMark>,
    pub provenance_score: Option<f64>,
    pub codon_policy_summary: String,
    pub next_step_recommendations: Vec<String>,
    pub future_materialization_seq_id: Option<SeqId>,
}

fn default_host_profile_catalog_schema() -> String {
    HOST_PROFILE_CATALOG_SCHEMA.to_string()
}

#[derive(Debug, Clone, Copy, Serialize, Deserialize, PartialEq, Eq, PartialOrd, Ord, Default)]
#[serde(rename_all = "snake_case")]
/// Biological role or design role attributed to one sequence interval.
pub enum ConstructRole {
    Promoter,
    Enhancer,
    Gene,
    Transcript,
    Exon,
    Utr5Prime,
    Cds,
    Utr3Prime,
    Terminator,
    Variant,
    Linker,
    Tag,
    SignalPeptide,
    LocalizationSignal,
    HomologyArm,
    FusionBoundary,
    RestrictionSite,
    SpliceBoundary,
    Tfbs,
    ContextBaggage,
    #[default]
    Other,
}

impl ConstructRole {
    pub fn as_str(self) -> &'static str {
        match self {
            Self::Promoter => "promoter",
            Self::Enhancer => "enhancer",
            Self::Gene => "gene",
            Self::Transcript => "transcript",
            Self::Exon => "exon",
            Self::Utr5Prime => "utr_5prime",
            Self::Cds => "cds",
            Self::Utr3Prime => "utr_3prime",
            Self::Terminator => "terminator",
            Self::Variant => "variant",
            Self::Linker => "linker",
            Self::Tag => "tag",
            Self::SignalPeptide => "signal_peptide",
            Self::LocalizationSignal => "localization_signal",
            Self::HomologyArm => "homology_arm",
            Self::FusionBoundary => "fusion_boundary",
            Self::RestrictionSite => "restriction_site",
            Self::SpliceBoundary => "splice_boundary",
            Self::Tfbs => "tfbs",
            Self::ContextBaggage => "context_baggage",
            Self::Other => "other",
        }
    }
}

#[derive(Debug, Clone, Copy, Serialize, Deserialize, PartialEq, Eq, PartialOrd, Ord, Default)]
#[serde(rename_all = "snake_case")]
/// Confidence class for one evidence node.
pub enum EvidenceClass {
    HardFact,
    #[default]
    ReliableAnnotation,
    ContextEvidence,
    SoftHypothesis,
    UserOverride,
}

impl EvidenceClass {
    pub fn as_str(self) -> &'static str {
        match self {
            Self::HardFact => "hard_fact",
            Self::ReliableAnnotation => "reliable_annotation",
            Self::ContextEvidence => "context_evidence",
            Self::SoftHypothesis => "soft_hypothesis",
            Self::UserOverride => "user_override",
        }
    }
}

#[derive(Debug, Clone, Copy, Serialize, Deserialize, PartialEq, Eq, PartialOrd, Ord, Default)]
#[serde(rename_all = "snake_case")]
/// Reasoning style used by one decision node.
pub enum DecisionMethod {
    #[default]
    HardRule,
    WeightedRule,
    FuzzyRule,
    ModelSummary,
    ParetoRank,
    UserOverride,
}

impl DecisionMethod {
    pub fn as_str(self) -> &'static str {
        match self {
            Self::HardRule => "hard_rule",
            Self::WeightedRule => "weighted_rule",
            Self::FuzzyRule => "fuzzy_rule",
            Self::ModelSummary => "model_summary",
            Self::ParetoRank => "pareto_rank",
            Self::UserOverride => "user_override",
        }
    }
}

#[derive(Debug, Clone, Copy, Serialize, Deserialize, PartialEq, Eq, PartialOrd, Ord, Default)]
#[serde(rename_all = "snake_case")]
/// Shared edit state for evidence, facts, and decisions.
pub enum EditableStatus {
    #[default]
    Draft,
    Accepted,
    Rejected,
    Locked,
}

impl EditableStatus {
    pub fn as_str(self) -> &'static str {
        match self {
            Self::Draft => "draft",
            Self::Accepted => "accepted",
            Self::Rejected => "rejected",
            Self::Locked => "locked",
        }
    }
}

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq, PartialOrd, Ord, Default)]
#[serde(rename_all = "snake_case")]
/// Where one evidence node primarily applies.
pub enum EvidenceScope {
    #[default]
    SequenceSpan,
    WholeConstruct,
    HostProfile,
    HostTransition,
    MediumCondition,
    HelperProfile,
}

impl EvidenceScope {
    pub fn as_str(&self) -> &'static str {
        match self {
            Self::SequenceSpan => "sequence_span",
            Self::WholeConstruct => "whole_construct",
            Self::HostProfile => "host_profile",
            Self::HostTransition => "host_transition",
            Self::MediumCondition => "medium_condition",
            Self::HelperProfile => "helper_profile",
        }
    }
}

#[derive(Debug, Clone, Copy, Serialize, Deserialize, PartialEq, Eq, PartialOrd, Ord, Default)]
#[serde(rename_all = "snake_case")]
/// Lifecycle role for one host step in a construct workflow.
pub enum HostLifecycleRole {
    #[default]
    Propagation,
    Expression,
    Intermediate,
    Storage,
}

impl HostLifecycleRole {
    pub fn as_str(self) -> &'static str {
        match self {
            Self::Propagation => "propagation",
            Self::Expression => "expression",
            Self::Intermediate => "intermediate",
            Self::Storage => "storage",
        }
    }
}

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
#[serde(default, deny_unknown_fields)]
/// One planned host step for propagation, transfer, or expression.
pub struct HostRouteStep {
    pub step_id: String,
    pub host_profile_id: String,
    pub role: HostLifecycleRole,
    pub rationale: String,
    pub notes: Vec<String>,
}

impl Default for HostRouteStep {
    fn default() -> Self {
        Self {
            step_id: String::new(),
            host_profile_id: String::new(),
            role: HostLifecycleRole::Propagation,
            rationale: String::new(),
            notes: vec![],
        }
    }
}

#[derive(Debug, Clone, Copy, Serialize, Deserialize, PartialEq, Eq, PartialOrd, Ord, Default)]
#[serde(rename_all = "snake_case")]
/// Style of blunt-end linker/adapter capture oligos chosen for one cloning intent.
pub enum AdapterCaptureStyle {
    #[default]
    Minimal,
    McsLike,
    Custom,
}

impl AdapterCaptureStyle {
    pub fn as_str(self) -> &'static str {
        match self {
            Self::Minimal => "minimal",
            Self::McsLike => "mcs_like",
            Self::Custom => "custom",
        }
    }
}

#[derive(Debug, Clone, Copy, Serialize, Deserialize, PartialEq, Eq, PartialOrd, Ord, Default)]
#[serde(rename_all = "snake_case")]
/// Intended protection mode when the insert already contains the adapter capture site.
pub enum AdapterCaptureProtectionMode {
    #[default]
    None,
    InsertMethylation,
}

impl AdapterCaptureProtectionMode {
    pub fn as_str(self) -> &'static str {
        match self {
            Self::None => "none",
            Self::InsertMethylation => "insert_methylation",
        }
    }
}

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
#[serde(default, deny_unknown_fields)]
/// One planned blunt-end adapter/linker ligation capture strategy.
pub struct AdapterRestrictionCapturePlan {
    pub capture_id: String,
    pub restriction_enzyme_name: String,
    pub blunt_insert_required: bool,
    pub adapter_style: AdapterCaptureStyle,
    pub protection_mode: AdapterCaptureProtectionMode,
    pub extra_retrieval_enzyme_names: Vec<String>,
    pub notes: Vec<String>,
}

impl Default for AdapterRestrictionCapturePlan {
    fn default() -> Self {
        Self {
            capture_id: String::new(),
            restriction_enzyme_name: String::new(),
            blunt_insert_required: false,
            adapter_style: AdapterCaptureStyle::Minimal,
            protection_mode: AdapterCaptureProtectionMode::None,
            extra_retrieval_enzyme_names: vec![],
            notes: vec![],
        }
    }
}

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
#[serde(default, deny_unknown_fields)]
/// Target outcome used to interpret evidence and derive construct candidates.
pub struct ConstructObjective {
    #[serde(default = "default_construct_objective_schema")]
    pub schema: String,
    pub objective_id: String,
    pub title: String,
    pub goal: String,
    pub host_species: Option<String>,
    pub cell_type: Option<String>,
    pub tissue: Option<String>,
    pub organelle: Option<String>,
    pub expression_intent: Option<String>,
    pub propagation_host_profile_id: Option<String>,
    pub expression_host_profile_id: Option<String>,
    pub host_route: Vec<HostRouteStep>,
    pub medium_conditions: Vec<String>,
    pub helper_profile_id: Option<String>,
    pub adapter_restriction_capture_plans: Vec<AdapterRestrictionCapturePlan>,
    pub required_host_traits: Vec<String>,
    pub forbidden_host_traits: Vec<String>,
    pub required_roles: Vec<ConstructRole>,
    pub forbidden_roles: Vec<ConstructRole>,
    pub preferred_routine_families: Vec<String>,
    pub notes: Vec<String>,
}

impl Default for ConstructObjective {
    fn default() -> Self {
        Self {
            schema: default_construct_objective_schema(),
            objective_id: String::new(),
            title: String::new(),
            goal: String::new(),
            host_species: None,
            cell_type: None,
            tissue: None,
            organelle: None,
            expression_intent: None,
            propagation_host_profile_id: None,
            expression_host_profile_id: None,
            host_route: vec![],
            medium_conditions: vec![],
            helper_profile_id: None,
            adapter_restriction_capture_plans: vec![],
            required_host_traits: vec![],
            forbidden_host_traits: vec![],
            required_roles: vec![],
            forbidden_roles: vec![],
            preferred_routine_families: vec![],
            notes: vec![],
        }
    }
}

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
#[serde(default, deny_unknown_fields)]
/// Sequence-linked evidence node in the construct reasoning graph.
pub struct DesignEvidence {
    #[serde(default = "default_design_evidence_schema")]
    pub schema: String,
    pub evidence_id: String,
    pub seq_id: SeqId,
    pub scope: EvidenceScope,
    pub start_0based: usize,
    pub end_0based_exclusive: usize,
    pub strand: Option<String>,
    pub host_profile_id: Option<String>,
    pub host_route_step_id: Option<String>,
    pub helper_profile_id: Option<String>,
    pub medium_condition_id: Option<String>,
    pub role: ConstructRole,
    pub evidence_class: EvidenceClass,
    pub label: String,
    pub rationale: String,
    pub score: Option<f64>,
    pub confidence: Option<f64>,
    pub specificity_bias: Option<f64>,
    pub sensitivity_bias: Option<f64>,
    pub context_tags: Vec<String>,
    pub provenance_kind: String,
    pub provenance_refs: Vec<String>,
    pub editable_status: EditableStatus,
    pub warnings: Vec<String>,
    pub notes: Vec<String>,
}

impl Default for DesignEvidence {
    fn default() -> Self {
        Self {
            schema: default_design_evidence_schema(),
            evidence_id: String::new(),
            seq_id: String::new(),
            scope: EvidenceScope::SequenceSpan,
            start_0based: 0,
            end_0based_exclusive: 0,
            strand: None,
            host_profile_id: None,
            host_route_step_id: None,
            helper_profile_id: None,
            medium_condition_id: None,
            role: ConstructRole::Other,
            evidence_class: EvidenceClass::ReliableAnnotation,
            label: String::new(),
            rationale: String::new(),
            score: None,
            confidence: None,
            specificity_bias: None,
            sensitivity_bias: None,
            context_tags: vec![],
            provenance_kind: String::new(),
            provenance_refs: vec![],
            editable_status: EditableStatus::Draft,
            warnings: vec![],
            notes: vec![],
        }
    }
}

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
#[serde(default, deny_unknown_fields)]
/// Curated host/strain record used by construct reasoning.
pub struct HostProfileRecord {
    pub profile_id: String,
    pub species: String,
    pub strain: String,
    pub aliases: Vec<String>,
    pub genotype_tags: Vec<String>,
    pub phenotype_tags: Vec<String>,
    pub notes: Vec<String>,
    pub source_notes: Vec<String>,
}

impl Default for HostProfileRecord {
    fn default() -> Self {
        Self {
            profile_id: String::new(),
            species: String::new(),
            strain: String::new(),
            aliases: vec![],
            genotype_tags: vec![],
            phenotype_tags: vec![],
            notes: vec![],
            source_notes: vec![],
        }
    }
}

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
#[serde(default, deny_unknown_fields)]
/// Human-editable catalog of host/strain records used by construct reasoning.
pub struct HostProfileCatalog {
    #[serde(default = "default_host_profile_catalog_schema")]
    pub schema: String,
    pub profiles: Vec<HostProfileRecord>,
}

impl Default for HostProfileCatalog {
    fn default() -> Self {
        Self {
            schema: default_host_profile_catalog_schema(),
            profiles: vec![],
        }
    }
}

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
#[serde(default, deny_unknown_fields)]
/// Curated helper/vector profile used by construct reasoning.
pub struct HelperConstructProfile {
    pub profile_id: String,
    pub helper_seq_id: Option<String>,
    pub helper_genome_id: Option<String>,
    pub vector_family: Option<String>,
    pub backbone_roles: Vec<String>,
    pub host_compatibility_tags: Vec<String>,
    pub notes: Vec<String>,
    pub source_notes: Vec<String>,
}

impl Default for HelperConstructProfile {
    fn default() -> Self {
        Self {
            profile_id: String::new(),
            helper_seq_id: None,
            helper_genome_id: None,
            vector_family: None,
            backbone_roles: vec![],
            host_compatibility_tags: vec![],
            notes: vec![],
            source_notes: vec![],
        }
    }
}

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
#[serde(default, deny_unknown_fields)]
/// Derived reusable fact in the construct reasoning graph.
pub struct DesignFact {
    #[serde(default = "default_design_fact_schema")]
    pub schema: String,
    pub fact_id: String,
    pub fact_type: String,
    pub label: String,
    pub rationale: String,
    pub based_on_evidence_ids: Vec<String>,
    pub value_json: serde_json::Value,
    pub confidence: Option<f64>,
    pub editable_status: EditableStatus,
}

impl Default for DesignFact {
    fn default() -> Self {
        Self {
            schema: default_design_fact_schema(),
            fact_id: String::new(),
            fact_type: String::new(),
            label: String::new(),
            rationale: String::new(),
            based_on_evidence_ids: vec![],
            value_json: serde_json::Value::Null,
            confidence: None,
            editable_status: EditableStatus::Draft,
        }
    }
}

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
#[serde(default, deny_unknown_fields)]
/// Explicit reasoning step that consumes evidence/facts and emits facts/candidates.
pub struct DesignDecisionNode {
    #[serde(default = "default_design_decision_node_schema")]
    pub schema: String,
    pub decision_id: String,
    pub decision_type: String,
    pub method: DecisionMethod,
    pub title: String,
    pub rationale: String,
    pub input_evidence_ids: Vec<String>,
    pub input_fact_ids: Vec<String>,
    pub output_fact_ids: Vec<String>,
    pub output_candidate_ids: Vec<String>,
    pub parameters_json: serde_json::Value,
    pub editable_status: EditableStatus,
}

impl Default for DesignDecisionNode {
    fn default() -> Self {
        Self {
            schema: default_design_decision_node_schema(),
            decision_id: String::new(),
            decision_type: String::new(),
            method: DecisionMethod::HardRule,
            title: String::new(),
            rationale: String::new(),
            input_evidence_ids: vec![],
            input_fact_ids: vec![],
            output_fact_ids: vec![],
            output_candidate_ids: vec![],
            parameters_json: serde_json::Value::Null,
            editable_status: EditableStatus::Draft,
        }
    }
}

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
#[serde(default, deny_unknown_fields)]
/// One candidate construct derived from the reasoning graph.
pub struct ConstructCandidate {
    #[serde(default = "default_construct_candidate_schema")]
    pub schema: String,
    pub candidate_id: String,
    pub objective_id: String,
    pub title: String,
    pub component_ids: Vec<String>,
    pub derived_from_fact_ids: Vec<String>,
    pub suggested_routine_ids: Vec<String>,
    pub compactness_score: Option<f64>,
    pub confidence_score: Option<f64>,
    pub cloning_complexity_score: Option<f64>,
    pub host_fit_score: Option<f64>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub protein_to_dna_handoff: Option<ProteinToDnaHandoffCandidate>,
    pub warnings: Vec<String>,
    pub notes: Vec<String>,
}

impl Default for ConstructCandidate {
    fn default() -> Self {
        Self {
            schema: default_construct_candidate_schema(),
            candidate_id: String::new(),
            objective_id: String::new(),
            title: String::new(),
            component_ids: vec![],
            derived_from_fact_ids: vec![],
            suggested_routine_ids: vec![],
            compactness_score: None,
            confidence_score: None,
            cloning_complexity_score: None,
            host_fit_score: None,
            protein_to_dna_handoff: None,
            warnings: vec![],
            notes: vec![],
        }
    }
}

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
#[serde(default, deny_unknown_fields)]
/// One inspectable annotation-grade interval derived from construct reasoning.
pub struct AnnotationCandidate {
    #[serde(default = "default_annotation_candidate_schema")]
    pub schema: String,
    pub annotation_id: String,
    pub evidence_id: String,
    pub seq_id: SeqId,
    pub start_0based: usize,
    pub end_0based_exclusive: usize,
    pub strand: Option<String>,
    pub role: ConstructRole,
    pub label: String,
    pub rationale: String,
    pub source_kind: String,
    pub supporting_fact_ids: Vec<String>,
    pub supporting_fact_labels: Vec<String>,
    pub supporting_decision_ids: Vec<String>,
    pub supporting_decision_titles: Vec<String>,
    pub transcript_context_status: Option<String>,
    pub effect_tags: Vec<String>,
    pub editable_status: EditableStatus,
    pub warnings: Vec<String>,
    pub notes: Vec<String>,
}

impl Default for AnnotationCandidate {
    fn default() -> Self {
        Self {
            schema: default_annotation_candidate_schema(),
            annotation_id: String::new(),
            evidence_id: String::new(),
            seq_id: String::new(),
            start_0based: 0,
            end_0based_exclusive: 0,
            strand: None,
            role: ConstructRole::Other,
            label: String::new(),
            rationale: String::new(),
            source_kind: String::new(),
            supporting_fact_ids: vec![],
            supporting_fact_labels: vec![],
            supporting_decision_ids: vec![],
            supporting_decision_titles: vec![],
            transcript_context_status: None,
            effect_tags: vec![],
            editable_status: EditableStatus::Draft,
            warnings: vec![],
            notes: vec![],
        }
    }
}

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
#[serde(default, deny_unknown_fields)]
/// Result of writing one accepted construct-reasoning annotation candidate back
/// into the sequence as an ordinary feature.
pub struct AnnotationCandidateWriteback {
    #[serde(default = "default_annotation_candidate_writeback_schema")]
    pub schema: String,
    pub graph_id: String,
    pub annotation_id: String,
    pub evidence_id: String,
    pub seq_id: SeqId,
    pub feature_kind: String,
    pub feature_label: String,
    pub start_0based: usize,
    pub end_0based_exclusive: usize,
    pub strand: Option<String>,
    pub source_kind: String,
    pub editable_status: EditableStatus,
    pub created: bool,
    pub already_present: bool,
    pub warnings: Vec<String>,
    pub notes: Vec<String>,
}

impl Default for AnnotationCandidateWriteback {
    fn default() -> Self {
        Self {
            schema: default_annotation_candidate_writeback_schema(),
            graph_id: String::new(),
            annotation_id: String::new(),
            evidence_id: String::new(),
            seq_id: String::new(),
            feature_kind: String::new(),
            feature_label: String::new(),
            start_0based: 0,
            end_0based_exclusive: 0,
            strand: None,
            source_kind: String::new(),
            editable_status: EditableStatus::Draft,
            created: false,
            already_present: false,
            warnings: vec![],
            notes: vec![],
        }
    }
}

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
#[serde(default, deny_unknown_fields)]
/// Complete sequence-linked reasoning graph for one design objective.
pub struct ConstructReasoningGraph {
    #[serde(default = "default_construct_reasoning_graph_schema")]
    pub schema: String,
    pub graph_id: String,
    pub seq_id: SeqId,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub op_id: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub run_id: Option<String>,
    pub objective: ConstructObjective,
    pub generated_at_unix_ms: u128,
    pub evidence: Vec<DesignEvidence>,
    pub facts: Vec<DesignFact>,
    pub decisions: Vec<DesignDecisionNode>,
    pub candidates: Vec<ConstructCandidate>,
    pub annotation_candidates: Vec<AnnotationCandidate>,
    pub notes: Vec<String>,
}

impl Default for ConstructReasoningGraph {
    fn default() -> Self {
        Self {
            schema: default_construct_reasoning_graph_schema(),
            graph_id: String::new(),
            seq_id: String::new(),
            op_id: None,
            run_id: None,
            objective: ConstructObjective::default(),
            generated_at_unix_ms: 0,
            evidence: vec![],
            facts: vec![],
            decisions: vec![],
            candidates: vec![],
            annotation_candidates: vec![],
            notes: vec![],
        }
    }
}

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
#[serde(default, deny_unknown_fields)]
/// Persisted collection of construct objectives and reasoning graphs.
pub struct ConstructReasoningStore {
    #[serde(default = "default_construct_reasoning_store_schema")]
    pub schema: String,
    pub updated_at_unix_ms: u128,
    pub objectives: BTreeMap<String, ConstructObjective>,
    pub graphs: BTreeMap<String, ConstructReasoningGraph>,
    pub preferred_graph_id: Option<String>,
}

impl Default for ConstructReasoningStore {
    fn default() -> Self {
        Self {
            schema: default_construct_reasoning_store_schema(),
            updated_at_unix_ms: 0,
            objectives: BTreeMap::new(),
            graphs: BTreeMap::new(),
            preferred_graph_id: None,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn construct_reasoning_graph_round_trip_preserves_schema_defaults() {
        let mut graph = ConstructReasoningGraph::default();
        graph.graph_id = "graph_demo".to_string();
        graph.seq_id = "seq1".to_string();
        graph.objective.objective_id = "obj1".to_string();
        graph.objective.goal = "demo".to_string();
        graph.objective.propagation_host_profile_id = Some("ecoli_k12".to_string());
        graph.objective.host_route.push(HostRouteStep {
            step_id: "step1".to_string(),
            host_profile_id: "ecoli_k12".to_string(),
            role: HostLifecycleRole::Propagation,
            rationale: "initial propagation".to_string(),
            ..HostRouteStep::default()
        });
        graph.evidence.push(DesignEvidence {
            evidence_id: "ev1".to_string(),
            seq_id: "seq1".to_string(),
            scope: EvidenceScope::HostTransition,
            host_profile_id: Some("ecoli_k12".to_string()),
            host_route_step_id: Some("step1".to_string()),
            role: ConstructRole::RestrictionSite,
            evidence_class: EvidenceClass::HardFact,
            label: "EcoRI".to_string(),
            rationale: "Restriction site".to_string(),
            ..DesignEvidence::default()
        });

        let value = serde_json::to_value(&graph).expect("serialize");
        let round_trip: ConstructReasoningGraph =
            serde_json::from_value(value).expect("deserialize");
        assert_eq!(round_trip.schema, CONSTRUCT_REASONING_GRAPH_SCHEMA);
        assert_eq!(round_trip.objective.schema, CONSTRUCT_OBJECTIVE_SCHEMA);
        assert_eq!(round_trip.evidence[0].schema, DESIGN_EVIDENCE_SCHEMA);
        assert_eq!(
            round_trip.objective.propagation_host_profile_id.as_deref(),
            Some("ecoli_k12")
        );
        assert_eq!(round_trip.objective.host_route.len(), 1);
        assert_eq!(round_trip.evidence[0].scope, EvidenceScope::HostTransition);
        assert_eq!(
            round_trip.evidence[0].host_profile_id.as_deref(),
            Some("ecoli_k12")
        );
    }
}
