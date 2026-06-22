//! Shared machine-readable GENtle contracts.
//!
//! This crate is intentionally scaffolded first in the workspace split so
//! stable cross-adapter payloads can move here before execution or GUI code.
//! The first extracted slice is intentionally small: stable identifier aliases,
//! shared analysis enums, and the portable engine error payload.

pub mod construct_reasoning;
pub mod dna_ladder;
pub mod gene_groups;
pub mod gene_sets;
pub mod orthologs;
pub mod reporter;

use serde::{Deserialize, Deserializer, Serialize};
use serde_json::{Value, json};
use std::{
    collections::{BTreeMap, BTreeSet, HashMap},
    error::Error,
    fmt,
    sync::LazyLock,
};

pub use construct_reasoning::{
    ANNOTATION_CANDIDATE_SCHEMA, ANNOTATION_CANDIDATE_SUMMARY_SCHEMA,
    ANNOTATION_CANDIDATE_WRITEBACK_SCHEMA, AdapterCaptureProtectionMode, AdapterCaptureStyle,
    AdapterRestrictionCapturePlan, AnnotationCandidate, AnnotationCandidateSummary,
    AnnotationCandidateWriteback, CONSTRUCT_CANDIDATE_SCHEMA, CONSTRUCT_OBJECTIVE_SCHEMA,
    CONSTRUCT_REASONING_GRAPH_SCHEMA, CONSTRUCT_REASONING_INSPECTION_ACTION_SCHEMA,
    CONSTRUCT_REASONING_STORE_SCHEMA, ConstructCandidate, ConstructObjective,
    ConstructReasoningGraph, ConstructReasoningInspectionAction,
    ConstructReasoningInspectionActionKind, ConstructReasoningRepeatFamilyProvenance,
    ConstructReasoningRiskTask, ConstructReasoningSeverity, ConstructReasoningStore,
    ConstructReasoningTaskSeverity, ConstructRole, DESIGN_DECISION_NODE_SCHEMA,
    DESIGN_EVIDENCE_SCHEMA, DESIGN_FACT_SCHEMA, DecisionMethod, DesignDecisionNode, DesignEvidence,
    DesignFact, EditableStatus, EvidenceClass, EvidenceScope, HOST_PROFILE_CATALOG_SCHEMA,
    HelperConstructProfile, HostLifecycleRole, HostProfileCatalog, HostProfileRecord,
    HostRouteStep, ProteinToDnaHandoffCandidate, ProteinToDnaHandoffCoverage,
    ProteinToDnaHandoffRankingGoal, ProteinToDnaHandoffStrategy,
};
pub use dna_ladder::{
    DNALadder, DNALadderBand, DNALadders, Ladder, LadderBand, LadderCatalog, LadderMolecule,
    RNALadder, RNALadderBand, RNALadders, default_dna_ladders, default_rna_ladders,
};
pub use gene_groups::{
    GENE_GROUP_CATALOG_SCHEMA, GENE_GROUP_DOCTOR_REPORT_SCHEMA, GENE_GROUP_DRAFT_REPORT_SCHEMA,
    GENE_GROUP_LIST_REPORT_SCHEMA, GENE_GROUP_RESOLVE_REPORT_SCHEMA, GENE_GROUP_SHOW_REPORT_SCHEMA,
    GeneGroupCatalog, GeneGroupCatalogSourceReport, GeneGroupDoctorReport, GeneGroupDraftReport,
    GeneGroupExternalMapping, GeneGroupExternalResource, GeneGroupListEntry, GeneGroupListReport,
    GeneGroupMember, GeneGroupRecord, GeneGroupResolveReport, GeneGroupShowReport,
};
pub use gene_sets::{
    GENE_SET_CUTRUN_REGULATORY_SUPPORT_SCHEMA, GENE_SET_PROMOTER_COHORT_SCHEMA,
    GENE_SET_RESOLUTION_SCHEMA, GeneSetCohortRelationship, GeneSetCohortRelationshipFlag,
    GeneSetCutRunEvaluationState, GeneSetCutRunMemberSupport, GeneSetCutRunRegulatorySupportReport,
    GeneSetCutRunSupportAggregate, GeneSetPromoterCohortReport, GeneSetPromoterWindow,
    GeneSetProvenanceRow, GeneSetRandomProvenance, GeneSetRequest, GeneSetResolutionReport,
    GeneSetResolvedMember, GeneSetUnresolvedMember,
};
pub use orthologs::{
    ORTHOLOG_PROMOTER_COHORT_SCHEMA, ORTHOLOG_PROMOTER_COMPARISON_SCHEMA, ORTHOLOG_RESOURCE_SCHEMA,
    OrthologAmbiguityPolicy, OrthologCutRunSupportRow, OrthologCutRunSupportStatus,
    OrthologExpressionAssignment, OrthologMappingRow, OrthologPairwiseTfbsSimilarity,
    OrthologPromoterCohortReport, OrthologPromoterCohortRequest, OrthologPromoterComparisonReport,
    OrthologPromoterRole, OrthologPromoterRow, OrthologResource, OrthologSequenceSimilarityRow,
    OrthologSpeciesAlias, OrthologTfbsPeakSummary, OrthologTfbsSummaryRow, OrthologUnresolvedRow,
};
pub use reporter::{
    PortBindingStatus, REPORTER_CATALOG_REPORT_SCHEMA, REPORTER_CATALOG_SCHEMA,
    REPORTER_CONSTRUCT_HANDOFF_SCHEMA, REPORTER_CORPUS_EXPORT_SCHEMA,
    REPORTER_RECOMMENDATION_SCHEMA, ReporterAnnotatedRecord, ReporterBackboneResolution,
    ReporterBackboneResolutionStatus, ReporterCatalog, ReporterCatalogReport,
    ReporterComputedAnnotation, ReporterConstraints, ReporterConstructHandoffCommand,
    ReporterConstructHandoffPlan, ReporterConstructHandoffProvenance, ReporterConstructPortBinding,
    ReporterConstructSelectedFragment, ReporterConstructSelectedReporter, ReporterCorpusExport,
    ReporterCorpusExportFormat, ReporterPreferenceWeights, ReporterQuarantinedRecord,
    ReporterRecommendation, ReporterRecommendationResult, ReporterRecord,
    ReporterRejectedCandidate, ReporterSourceRef, ReporterSpectralProfile,
};

/// Stable identifier for one sequence entry stored in project state.
pub type SeqId = String;
/// Stable identifier for one executed operation journal row.
pub type OpId = String;
/// Caller-supplied identifier that groups operations into one workflow/run.
pub type RunId = String;

/// Provider capability catalog for vendor/CRO/eProcurement integrations.
pub const EXTERNAL_SERVICE_PROVIDER_CATALOG_SCHEMA: &str =
    "gentle.external_service_provider_catalog.v1";
/// Overlay-discoverable provider behavior/configuration catalog.
pub const EXTERNAL_SERVICE_PROVIDER_CONFIG_SCHEMA: &str =
    "gentle.external_service_provider_config.v1";
/// Doctor report for external-service provider configuration overlays.
pub const EXTERNAL_SERVICE_PROVIDER_CONFIG_DOCTOR_SCHEMA: &str =
    "gentle.external_service_provider_config_doctor.v1";
/// Portable request contract for external-service preflight and handoff.
pub const EXTERNAL_SERVICE_REQUEST_SCHEMA: &str = "gentle.external_service_request.v1";
/// Local capability/eligibility report for an external-service request.
pub const EXTERNAL_SERVICE_PREFLIGHT_SCHEMA: &str = "gentle.external_service_preflight.v1";
/// Quote/handoff report for external services that have not been submitted.
pub const EXTERNAL_SERVICE_QUOTE_SCHEMA: &str = "gentle.external_service_quote.v1";
/// Provider/service-kind routing report for ambiguous sequence-delivery requests.
pub const EXTERNAL_SERVICE_DELIVERY_ROUTE_SCHEMA: &str =
    "gentle.external_service_delivery_route.v1";
/// Input contract for selecting an external-service route from sequence context.
pub const EXTERNAL_SERVICE_DELIVERY_ROUTE_REQUEST_SCHEMA: &str =
    "gentle.external_service_delivery_route_request.v1";

/// Engine-owned list of external providers and their supported service kinds.
#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
pub struct ExternalServiceProviderCatalog {
    pub schema: String,
    pub generated_at_unix_ms: u128,
    pub providers: Vec<ExternalServiceProviderRecord>,
    pub summary_lines: Vec<String>,
}

/// One external provider entry, such as the initial GeneArt capability row.
#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
pub struct ExternalServiceProviderRecord {
    pub provider: String,
    pub display_name: String,
    pub support_status: String,
    pub website_url: String,
    pub dashboard_url: String,
    pub api_documentation_url: Option<String>,
    pub capabilities: Vec<ExternalServiceCapability>,
    pub account_enablement_notes: Vec<String>,
    pub warnings: Vec<String>,
}

/// Overlay-discoverable provider configuration used to build provider records.
#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
pub struct ExternalServiceProviderConfigCatalog {
    pub schema: String,
    pub providers: Vec<ExternalServiceProviderConfigRecord>,
    pub summary_lines: Vec<String>,
}

/// One provider configuration row, including provider-specific handoff mapping.
#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
pub struct ExternalServiceProviderConfigRecord {
    pub provider: String,
    pub display_name: String,
    pub support_status: String,
    pub website_url: String,
    pub dashboard_url: String,
    pub api_documentation_url: Option<String>,
    pub capabilities: Vec<ExternalServiceCapability>,
    pub channels: Vec<ExternalServiceChannelConfig>,
    pub product_templates: Vec<ExternalServiceProductTemplateMapping>,
    pub validation_rules: Vec<ExternalServiceValidationRule>,
    pub default_delivery_hints: Vec<String>,
    pub default_purification_hints: Vec<String>,
    pub default_qc_hints: Vec<String>,
    pub required_followup: Vec<String>,
    pub account_enablement_notes: Vec<String>,
    pub warnings: Vec<String>,
}

/// Vendor contact or handoff channel such as WOP or email+Excel.
#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
pub struct ExternalServiceChannelConfig {
    pub channel: String,
    pub display_name: String,
    pub url: Option<String>,
    pub email: Option<String>,
    pub notes: Vec<String>,
}

/// Mapping from a provider-neutral service kind to a vendor order template.
#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
pub struct ExternalServiceProductTemplateMapping {
    pub service_kind: String,
    pub product_name: String,
    pub channel: String,
    pub template_url: Option<String>,
    pub template_format: Option<String>,
    pub local_template_path: Option<String>,
    pub field_map: BTreeMap<String, String>,
    pub notes: Vec<String>,
}

/// Deterministic provider-specific validation rule for request preflight.
#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
pub struct ExternalServiceValidationRule {
    pub service_kind: String,
    pub required_source_fields: Vec<String>,
    pub required_delivery_fields: Vec<String>,
    pub warnings: Vec<String>,
    pub required_followup: Vec<String>,
}

/// One provider-configuration source inspected by the doctor route.
#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
pub struct ExternalServiceProviderConfigSourceReport {
    pub scope: String,
    pub path: String,
    pub status: String,
    pub provider_count: usize,
    pub sha1: Option<String>,
    pub warnings: Vec<String>,
    pub errors: Vec<String>,
}

/// Doctor report for external-service provider configuration overlays.
#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
pub struct ExternalServiceProviderConfigDoctorReport {
    pub schema: String,
    pub catalog_label: String,
    pub source_count: usize,
    pub parsed_source_count: usize,
    pub provider_count: usize,
    pub warning_count: usize,
    pub error_count: usize,
    pub sources: Vec<ExternalServiceProviderConfigSourceReport>,
    pub warnings: Vec<String>,
    pub errors: Vec<String>,
}

/// One service kind offered by an external provider.
#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
pub struct ExternalServiceCapability {
    pub service_kind: String,
    pub track: String,
    pub display_name: String,
    pub quote_handoff_supported: bool,
    pub direct_api_documented: bool,
    pub direct_api_implemented: bool,
    pub supported_submission_modes: Vec<String>,
    pub status_tracking: String,
    pub artifact_kinds: Vec<String>,
    pub notes: Vec<String>,
}

/// Vendor-neutral project request accepted by external-service preflight/quote routes.
#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(default)]
pub struct ExternalServiceRequest {
    pub schema: String,
    pub provider: String,
    pub service_kind: String,
    pub source_target: Value,
    pub optimization_target: Option<Value>,
    pub vector_spec: Option<Value>,
    pub delivery_options: Option<Value>,
    pub commercial_context_ref: Option<String>,
    pub return_spec: ExternalServiceReturnSpec,
    pub request_metadata: Option<Value>,
}

impl Default for ExternalServiceRequest {
    fn default() -> Self {
        Self {
            schema: EXTERNAL_SERVICE_REQUEST_SCHEMA.to_string(),
            provider: String::new(),
            service_kind: String::new(),
            source_target: Value::Null,
            optimization_target: None,
            vector_spec: None,
            delivery_options: None,
            commercial_context_ref: None,
            return_spec: ExternalServiceReturnSpec::default(),
            request_metadata: None,
        }
    }
}

/// Caller preferences for the payloads GENtle should return to automation.
#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(default)]
pub struct ExternalServiceReturnSpec {
    pub requested_payloads: Vec<String>,
    pub inline_max_bytes: Option<usize>,
    pub redact_commercial_fields: bool,
    pub prefer_artifact_bundle: bool,
}

impl Default for ExternalServiceReturnSpec {
    fn default() -> Self {
        Self {
            requested_payloads: vec!["quote_metadata".to_string(), "handoff_bundle".to_string()],
            inline_max_bytes: Some(32 * 1024),
            redact_commercial_fields: true,
            prefer_artifact_bundle: true,
        }
    }
}

/// Provider-neutral request for classifying "deliver this sequence" intent.
#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(default)]
pub struct ExternalServiceDeliveryRouteRequest {
    pub schema: String,
    pub source_target: Value,
    pub optimization_target: Option<Value>,
    pub vector_spec: Option<Value>,
    pub delivery_options: Option<Value>,
    pub commercial_context_ref: Option<String>,
    pub return_spec: ExternalServiceReturnSpec,
    pub request_metadata: Option<Value>,
}

impl Default for ExternalServiceDeliveryRouteRequest {
    fn default() -> Self {
        Self {
            schema: EXTERNAL_SERVICE_DELIVERY_ROUTE_REQUEST_SCHEMA.to_string(),
            source_target: Value::Null,
            optimization_target: None,
            vector_spec: None,
            delivery_options: None,
            commercial_context_ref: None,
            return_spec: ExternalServiceReturnSpec::default(),
            request_metadata: None,
        }
    }
}

/// Candidate provider route for a classified sequence-delivery request.
#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
pub struct ExternalServiceDeliveryRouteCandidate {
    pub provider: String,
    pub provider_display_name: String,
    pub service_kind: String,
    pub service_display_name: String,
    pub confidence: String,
    pub rationale: Vec<String>,
    pub request: ExternalServiceRequest,
}

/// Deterministic routing report for generic sequence-delivery wording.
#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
pub struct ExternalServiceDeliveryRouteReport {
    pub schema: String,
    pub generated_at_unix_ms: u128,
    pub status: String,
    pub molecule_type: String,
    pub sequence_kind: String,
    pub sequence_count: usize,
    pub sequence_length: Option<usize>,
    pub max_sequence_length: Option<usize>,
    pub length_unit: String,
    pub recommended_provider: Option<String>,
    pub recommended_service_kind: Option<String>,
    pub candidates: Vec<ExternalServiceDeliveryRouteCandidate>,
    pub summary_lines: Vec<String>,
    pub rationale: Vec<String>,
    pub clarification_questions: Vec<String>,
    pub warnings: Vec<String>,
}

/// Deterministic local eligibility/capability report for one service request.
#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
pub struct ExternalServicePreflightReport {
    pub schema: String,
    pub generated_at_unix_ms: u128,
    pub provider: String,
    pub provider_display_name: String,
    pub service_kind: String,
    pub capability_status: String,
    pub eligible: bool,
    pub quote_handoff_available: bool,
    pub direct_submission_available: bool,
    pub supported_submission_modes: Vec<String>,
    pub blocking_issues: Vec<String>,
    pub warnings: Vec<String>,
    pub estimated_turnaround: Option<String>,
    pub estimated_cost_hint: Option<String>,
    pub required_followup: Vec<String>,
    pub dashboard_links: Vec<ExternalServiceLink>,
    pub request_summary: Vec<String>,
}

/// Service-ready quote/handoff packet; V1 never implies vendor submission.
#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
pub struct ExternalServiceQuoteReport {
    pub schema: String,
    pub generated_at_unix_ms: u128,
    pub provider: String,
    pub service_kind: String,
    pub quote_status: String,
    pub quote_mode: String,
    pub preflight: ExternalServicePreflightReport,
    pub dashboard_links: Vec<ExternalServiceLink>,
    pub required_followup: Vec<String>,
    pub service_ready_bundle: ExternalServiceArtifactBundle,
    pub return_spec: ExternalServiceReturnSpec,
    pub warnings: Vec<String>,
}

/// Provider documentation/dashboard link attached to an external-service report.
#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
pub struct ExternalServiceLink {
    pub label: String,
    pub url: String,
    pub purpose: String,
}

/// Bundle of external-service artifacts selected by a request return spec.
#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
pub struct ExternalServiceArtifactBundle {
    pub schema: String,
    pub provider: String,
    pub service_kind: String,
    pub artifact_id: String,
    pub local_files: Vec<ExternalServiceArtifactRef>,
    pub inline_payloads: Vec<ExternalServiceInlinePayload>,
    pub notes: Vec<String>,
}

/// File-backed artifact reference owned by a GENtle report.
#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
pub struct ExternalServiceArtifactRef {
    pub artifact_kind: String,
    pub path: String,
    pub checksum_sha256: Option<String>,
    pub description: String,
}

/// Small inline artifact payload for handoff summaries and redacted request JSON.
#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
pub struct ExternalServiceInlinePayload {
    pub payload_kind: String,
    pub content_type: String,
    pub text: String,
    pub description: String,
}
/// Stable identifier for one lineage graph node.
pub type NodeId = String;
/// Stable identifier for one wet-lab-style container record.
pub type ContainerId = String;

#[derive(Debug, Clone, Serialize, Deserialize)]
/// Sequence SVG render layout requested by shell/CLI/export adapters.
pub enum RenderSvgMode {
    Linear,
    Circular,
}

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
    #[serde(rename = "pcr.assay.pair.with_tail")]
    PcrAssayPairWithTail,
    #[serde(rename = "pcr.oe.substitution")]
    PcrOeSubstitution,
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
            Self::PcrAssayPairWithTail => "pcr.assay.pair.with_tail",
            Self::PcrOeSubstitution => "pcr.oe.substitution",
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
            Self::PcrAssayPairWithTail => "PCR Assay (insertion-first tailed pair-PCR)",
            Self::PcrOeSubstitution => "PCR Mutagenesis (overlap-extension substitution)",
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
            Self::PcrAssayPairWithTail => {
                "Insertion-first pair-PCR strip with anchored 5' extensions and carried-in product tails"
            }
            Self::PcrOeSubstitution => {
                "Six-step overlap-extension substitution strip with primer set a-f and strand-specific anneal/fill states"
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
            "pcr.assay.pair.with_tail"
            | "pcr.assay.pair.tailed"
            | "pcr.assay.tailed"
            | "pcr.tailed"
            | "tailed_pcr" => Some(Self::PcrAssayPairWithTail),
            "pcr.oe.substitution"
            | "pcr.oe"
            | "oe_pcr"
            | "oe-pcr"
            | "overlap_extension_pcr"
            | "overlap-extension-pcr" => Some(Self::PcrOeSubstitution),
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
            Self::PcrAssayPairWithTail,
            Self::PcrOeSubstitution,
            Self::PcrAssayQpcr,
        ]
    }
}

#[derive(Debug, Clone, Copy, Serialize, Deserialize, PartialEq, Eq)]
#[serde(rename_all = "snake_case")]
/// Strand-contextual anchor extension side.
///
/// Interpretation is biological 5'/3' relative to anchor strand, not absolute
/// genomic coordinate direction.
pub enum GenomeAnchorSide {
    FivePrime,
    ThreePrime,
}

impl GenomeAnchorSide {
    pub fn as_str(self) -> &'static str {
        match self {
            Self::FivePrime => "5prime",
            Self::ThreePrime => "3prime",
        }
    }
}

#[derive(Debug, Clone, Copy, Serialize, Deserialize, PartialEq, Eq)]
#[serde(rename_all = "snake_case")]
/// Annotation projection policy for prepared/remote genome extraction.
pub enum GenomeAnnotationScope {
    None,
    Core,
    Full,
}

impl GenomeAnnotationScope {
    pub fn as_str(self) -> &'static str {
        match self {
            Self::None => "none",
            Self::Core => "core",
            Self::Full => "full",
        }
    }
}

#[derive(Debug, Clone, Copy, Serialize, Deserialize, PartialEq, Eq, Default)]
#[serde(rename_all = "snake_case")]
/// Interval policy for prepared-reference gene extraction.
pub enum GenomeGeneExtractMode {
    #[default]
    Gene,
    CodingWithPromoter,
}

impl GenomeGeneExtractMode {
    pub fn as_str(self) -> &'static str {
        match self {
            Self::Gene => "gene",
            Self::CodingWithPromoter => "coding_with_promoter",
        }
    }
}

#[derive(Debug, Clone, Copy, Serialize, Deserialize, PartialEq, Eq, Default)]
#[serde(rename_all = "snake_case")]
/// External track source format used by tracked-import subscriptions.
pub enum GenomeTrackSource {
    #[default]
    Bed,
    BigWig,
    Vcf,
}

impl GenomeTrackSource {
    pub fn from_path(path: &str) -> Self {
        let lower = path.trim().to_ascii_lowercase();
        if lower.ends_with(".bw") || lower.ends_with(".bigwig") {
            Self::BigWig
        } else if lower.ends_with(".vcf") || lower.ends_with(".vcf.gz") {
            Self::Vcf
        } else {
            Self::Bed
        }
    }

    pub fn label(self) -> &'static str {
        match self {
            Self::Bed => "BED",
            Self::BigWig => "BigWig",
            Self::Vcf => "VCF",
        }
    }
}

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
#[serde(default)]
/// Persisted genome-track import subscription shared by GUI and shell routes.
#[derive(Default)]
pub struct GenomeTrackSubscription {
    pub source: GenomeTrackSource,
    pub path: String,
    pub track_name: Option<String>,
    pub min_score: Option<f64>,
    pub max_score: Option<f64>,
    pub clear_existing: bool,
}

#[derive(Debug, Clone, Copy, Serialize, Deserialize, PartialEq, Eq, Default)]
#[serde(rename_all = "snake_case")]
/// Primer-design backend requested by adapter-facing primer/qPCR commands.
pub enum PrimerDesignBackend {
    #[default]
    Auto,
    Internal,
    Primer3,
}

impl PrimerDesignBackend {
    pub fn as_str(self) -> &'static str {
        match self {
            Self::Auto => "auto",
            Self::Internal => "internal",
            Self::Primer3 => "primer3",
        }
    }
}

#[derive(Debug, Clone, Copy, Serialize, Deserialize, PartialEq, Eq, Default)]
#[serde(rename_all = "snake_case")]
/// Whether primer specificity should merely be reported or enforced by a
/// design operation.
pub enum PrimerSpecificityCheckMode {
    #[default]
    None,
    ReportOnly,
    RequirePass,
}

impl PrimerSpecificityCheckMode {
    pub fn as_str(self) -> &'static str {
        match self {
            Self::None => "none",
            Self::ReportOnly => "report_only",
            Self::RequirePass => "require_pass",
        }
    }
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(default, deny_unknown_fields)]
/// Local BLAST specificity policy shared by standalone confirmation and future
/// design-time filtering.
pub struct PrimerSpecificityPolicy {
    pub specificity_check: PrimerSpecificityCheckMode,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub specificity_target_genome_id: Option<String>,
    pub max_target_amplicon_bp: usize,
    pub min_primer_coverage_fraction: f64,
    pub max_3prime_mismatches: usize,
    pub three_prime_window_bp: usize,
    pub min_total_mismatches_to_unintended_target: usize,
    pub allow_same_gene_splice_variants: bool,
    pub max_hits_per_primer: usize,
    pub avoid_known_variants: bool,
    pub avoid_rmsk_repeats: bool,
    pub avoid_low_complexity: bool,
}

impl Default for PrimerSpecificityPolicy {
    fn default() -> Self {
        Self {
            specificity_check: PrimerSpecificityCheckMode::ReportOnly,
            specificity_target_genome_id: None,
            max_target_amplicon_bp: 4_000,
            min_primer_coverage_fraction: 0.80,
            max_3prime_mismatches: 0,
            three_prime_window_bp: 5,
            min_total_mismatches_to_unintended_target: 2,
            allow_same_gene_splice_variants: false,
            max_hits_per_primer: 500,
            avoid_known_variants: false,
            avoid_rmsk_repeats: false,
            avoid_low_complexity: false,
        }
    }
}

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq)]
#[serde(rename_all = "snake_case")]
/// Transcript-aware qPCR design intent for splicing-driven assays.
#[derive(Default)]
pub enum QpcrTranscriptTargetingMode {
    #[default]
    SharedGene,
    DistinguishTranscript,
}

impl QpcrTranscriptTargetingMode {
    pub fn as_str(&self) -> &'static str {
        match self {
            Self::SharedGene => "shared_gene",
            Self::DistinguishTranscript => "distinguish_transcript",
        }
    }
}

#[derive(Debug, Clone, Copy, Serialize, Deserialize, PartialEq, Eq, Default)]
#[serde(rename_all = "snake_case")]
/// Which transcript-specific evidence kind a transcript-distinguishing qPCR
/// assay must satisfy.
pub enum QpcrTranscriptSpecificityEvidence {
    #[default]
    JunctionOnly,
    UniqueExonOrChain,
    EitherPreferJunction,
}

impl QpcrTranscriptSpecificityEvidence {
    pub fn as_str(self) -> &'static str {
        match self {
            Self::JunctionOnly => "junction_only",
            Self::UniqueExonOrChain => "unique_exon_or_chain",
            Self::EitherPreferJunction => "either_prefer_junction",
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize, Default)]
#[serde(rename_all = "snake_case")]
/// Transcript-row ordering policy for cDNA PCR/qPCR transcript maps.
pub enum CdnaAssayTranscriptOrder {
    #[default]
    TranscriptId,
    GenomicFirstExon,
    GenomicLastExon,
    AntisenseFirstExon,
}

impl CdnaAssayTranscriptOrder {
    pub fn as_str(self) -> &'static str {
        match self {
            Self::TranscriptId => "transcript_id",
            Self::GenomicFirstExon => "genomic_first_exon",
            Self::GenomicLastExon => "genomic_last_exon",
            Self::AntisenseFirstExon => "antisense_first_exon",
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize, Default)]
#[serde(rename_all = "snake_case")]
/// Coordinate system used by portable cDNA PCR/qPCR transcript maps.
pub enum CdnaAssayTranscriptMapCoordinateMode {
    #[default]
    Cdna,
    GenomicAligned,
}

impl CdnaAssayTranscriptMapCoordinateMode {
    pub fn as_str(self) -> &'static str {
        match self {
            Self::Cdna => "cdna",
            Self::GenomicAligned => "genomic_aligned",
        }
    }
}

#[derive(Debug, Clone, Copy, Serialize, Deserialize, PartialEq, Eq, Default)]
#[serde(rename_all = "snake_case")]
/// Restriction-cloning handoff mode for tailed PCR primer workflows.
pub enum RestrictionCloningPcrHandoffMode {
    #[default]
    SingleSite,
    DirectedPair,
}

impl RestrictionCloningPcrHandoffMode {
    pub fn as_str(self) -> &'static str {
        match self {
            Self::SingleSite => "single_site",
            Self::DirectedPair => "directed_pair",
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
/// Cleanup mode for conservative prepared-cache cleanup.
pub enum PreparedCacheCleanupMode {
    BlastDbOnly,
    DerivedIndexesOnly,
    SelectedPreparedInstalls,
    AllPreparedInCache,
}

impl PreparedCacheCleanupMode {
    pub fn label(self) -> &'static str {
        match self {
            Self::BlastDbOnly => "blast_db_only",
            Self::DerivedIndexesOnly => "derived_indexes_only",
            Self::SelectedPreparedInstalls => "selected_prepared_installs",
            Self::AllPreparedInCache => "all_prepared_in_cache",
        }
    }

    pub fn allows_orphaned_remnants(self) -> bool {
        matches!(
            self,
            Self::SelectedPreparedInstalls | Self::AllPreparedInCache
        )
    }
}

#[derive(Debug, Clone, Serialize, Deserialize)]
/// Request payload for deterministic prepared-cache cleanup.
pub struct PreparedCacheCleanupRequest {
    pub mode: PreparedCacheCleanupMode,
    #[serde(default)]
    pub cache_roots: Vec<String>,
    #[serde(default)]
    pub prepared_ids: Vec<String>,
    #[serde(default)]
    pub prepared_paths: Vec<String>,
    #[serde(default)]
    pub include_orphaned_remnants: bool,
}

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq)]
/// Named boundary selector for local sequence anchors.
pub enum AnchorBoundary {
    Start,
    End,
    Middle,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
/// Directional expansion selector for local sequence anchors.
pub enum AnchorDirection {
    Upstream,
    Downstream,
}

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq)]
/// Root-independent sequence anchor used by candidate and extraction commands.
pub enum SequenceAnchor {
    Position {
        zero_based: usize,
    },
    FeatureBoundary {
        feature_kind: Option<String>,
        feature_label: Option<String>,
        boundary: AnchorBoundary,
        occurrence: Option<usize>,
    },
}

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
/// Per-transcription-factor threshold override used by TFBS shell/operation
/// contracts.
pub struct TfThresholdOverride {
    pub tf: String,
    pub min_llr_bits: Option<f64>,
    pub min_llr_quantile: Option<f64>,
}

#[derive(Debug, Clone, Copy, Serialize, Deserialize, PartialEq, Eq, Default)]
#[serde(rename_all = "snake_case")]
/// How feature-based candidate queries turn matching annotations into geometry.
///
/// Use this when finding intervals around features or feature boundaries.
pub enum CandidateFeatureGeometryMode {
    #[default]
    FeatureSpan,
    FeatureParts,
    FeatureBoundaries,
}

impl CandidateFeatureGeometryMode {
    pub fn as_str(self) -> &'static str {
        match self {
            Self::FeatureSpan => "feature_span",
            Self::FeatureParts => "feature_parts",
            Self::FeatureBoundaries => "feature_boundaries",
        }
    }
}

#[derive(Debug, Clone, Copy, Serialize, Deserialize, PartialEq, Eq, Default)]
#[serde(rename_all = "snake_case")]
/// Which boundary of a matched feature is eligible when boundary mode is used.
pub enum CandidateFeatureBoundaryMode {
    #[default]
    Any,
    FivePrime,
    ThreePrime,
    Start,
    End,
}

impl CandidateFeatureBoundaryMode {
    pub fn as_str(self) -> &'static str {
        match self {
            Self::Any => "any",
            Self::FivePrime => "five_prime",
            Self::ThreePrime => "three_prime",
            Self::Start => "start",
            Self::End => "end",
        }
    }
}

#[derive(Debug, Clone, Copy, Serialize, Deserialize, PartialEq, Eq, Default)]
#[serde(rename_all = "snake_case")]
/// Strand relation required between a candidate query and matched feature.
pub enum CandidateFeatureStrandRelation {
    #[default]
    Any,
    Same,
    Opposite,
}

impl CandidateFeatureStrandRelation {
    pub fn as_str(self) -> &'static str {
        match self {
            Self::Any => "any",
            Self::Same => "same",
            Self::Opposite => "opposite",
        }
    }
}

#[derive(Debug, Clone, Copy, Serialize, Deserialize, PartialEq, Eq)]
#[serde(rename_all = "snake_case")]
/// Deterministic set algebra supported by candidate-set combination commands.
pub enum CandidateSetOperator {
    Union,
    Intersect,
    Subtract,
}

impl CandidateSetOperator {
    pub fn parse(raw: &str) -> Option<Self> {
        match raw.trim().to_ascii_lowercase().as_str() {
            "union" => Some(Self::Union),
            "intersect" | "intersection" => Some(Self::Intersect),
            "subtract" | "difference" => Some(Self::Subtract),
            _ => None,
        }
    }

    pub fn as_str(self) -> &'static str {
        match self {
            Self::Union => "union",
            Self::Intersect => "intersect",
            Self::Subtract => "subtract",
        }
    }
}

#[derive(Debug, Clone, Copy, Serialize, Deserialize, PartialEq, Eq, Default)]
#[serde(rename_all = "snake_case")]
pub enum CandidateObjectiveDirection {
    #[default]
    Maximize,
    Minimize,
}

impl CandidateObjectiveDirection {
    pub fn as_str(self) -> &'static str {
        match self {
            Self::Maximize => "maximize",
            Self::Minimize => "minimize",
        }
    }
}

#[derive(Debug, Clone, Serialize, Deserialize)]
/// One objective dimension for Pareto-frontier ranking.
pub struct CandidateObjectiveSpec {
    pub metric: String,
    #[serde(default)]
    pub direction: CandidateObjectiveDirection,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
/// Weighted scalar objective term used by candidate objective scoring.
pub struct CandidateWeightedObjectiveTerm {
    pub metric: String,
    pub weight: f64,
    #[serde(default)]
    pub direction: CandidateObjectiveDirection,
}

#[derive(Debug, Clone, Copy, Serialize, Deserialize, PartialEq, Eq, Default)]
#[serde(rename_all = "snake_case")]
/// Stable tie-breaker used after objective scores compare equal.
pub enum CandidateTieBreakPolicy {
    #[default]
    SeqStartEnd,
    SeqEndStart,
    LengthAscending,
    LengthDescending,
    SequenceLexicographic,
}

impl CandidateTieBreakPolicy {
    pub fn as_str(self) -> &'static str {
        match self {
            Self::SeqStartEnd => "seq_start_end",
            Self::SeqEndStart => "seq_end_start",
            Self::LengthAscending => "length_ascending",
            Self::LengthDescending => "length_descending",
            Self::SequenceLexicographic => "sequence_lexicographic",
        }
    }
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
/// One typed parameter exposed by a candidate macro template.
pub struct CandidateMacroTemplateParam {
    pub name: String,
    pub default_value: Option<String>,
    pub required: bool,
}

#[cfg(test)]
mod shell_contract_extraction_tests {
    use super::{
        CandidateFeatureGeometryMode, CandidateSetOperator, CandidateTieBreakPolicy,
        GenomeAnnotationScope, GenomeTrackSource, PreparedCacheCleanupMode,
        PrimerSpecificityCheckMode, ProtocolCartoonKind, QpcrTranscriptSpecificityEvidence,
        QpcrTranscriptTargetingMode,
    };

    #[test]
    fn moved_shell_payload_methods_keep_stable_spellings() {
        assert_eq!(GenomeAnnotationScope::Core.as_str(), "core");
        assert_eq!(
            QpcrTranscriptTargetingMode::SharedGene.as_str(),
            "shared_gene"
        );
        assert_eq!(
            QpcrTranscriptSpecificityEvidence::EitherPreferJunction.as_str(),
            "either_prefer_junction"
        );
        assert_eq!(
            PreparedCacheCleanupMode::SelectedPreparedInstalls.label(),
            "selected_prepared_installs"
        );
        assert_eq!(
            PrimerSpecificityCheckMode::RequirePass.as_str(),
            "require_pass"
        );
        assert_eq!(ProtocolCartoonKind::PcrAssayQpcr.id(), "pcr.assay.qpcr");
        assert_eq!(
            CandidateFeatureGeometryMode::FeatureBoundaries.as_str(),
            "feature_boundaries"
        );
        assert_eq!(
            CandidateTieBreakPolicy::LengthDescending.as_str(),
            "length_descending"
        );
    }

    #[test]
    fn genome_track_source_detects_common_extensions() {
        assert_eq!(
            GenomeTrackSource::from_path("signal.bigWig").label(),
            "BigWig"
        );
        assert_eq!(
            GenomeTrackSource::from_path("variants.vcf.gz").label(),
            "VCF"
        );
        assert_eq!(GenomeTrackSource::from_path("peaks.bed").label(), "BED");
    }

    #[test]
    fn protocol_cartoon_kind_accepts_shell_aliases() {
        assert_eq!(
            ProtocolCartoonKind::parse_id("oe-pcr"),
            Some(ProtocolCartoonKind::PcrOeSubstitution)
        );
        assert!(ProtocolCartoonKind::parse_id("unknown.protocol").is_none());
    }

    #[test]
    fn candidate_set_operator_accepts_shell_aliases() {
        assert_eq!(
            CandidateSetOperator::parse("intersection"),
            Some(CandidateSetOperator::Intersect)
        );
        assert_eq!(
            CandidateSetOperator::parse("difference"),
            Some(CandidateSetOperator::Subtract)
        );
    }
}

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
    #[serde(
        rename = "small_tube_4x6",
        alias = "small_tube4x6",
        alias = "SmallTube4x6"
    )]
    SmallTube4x6,
    #[serde(rename = "plate_6", alias = "plate6", alias = "Plate6")]
    Plate6,
    #[serde(rename = "plate_12")]
    Plate12,
    #[serde(rename = "plate_24")]
    Plate24,
    #[serde(rename = "plate_48")]
    Plate48,
    #[serde(rename = "plate_96", alias = "plate96", alias = "Plate96")]
    Plate96,
    #[serde(rename = "plate_384", alias = "plate384", alias = "Plate384")]
    Plate384,
    Custom,
}

impl RackProfileKind {
    pub fn as_str(self) -> &'static str {
        match self {
            Self::SmallTube4x6 => "small_tube_4x6",
            Self::Plate6 => "plate_6",
            Self::Plate12 => "plate_12",
            Self::Plate24 => "plate_24",
            Self::Plate48 => "plate_48",
            Self::Plate96 => "plate_96",
            Self::Plate384 => "plate_384",
            Self::Custom => "custom",
        }
    }

    pub fn dimensions(self) -> (usize, usize) {
        match self {
            Self::SmallTube4x6 => (4, 6),
            Self::Plate6 => (2, 3),
            Self::Plate12 => (3, 4),
            Self::Plate24 => (4, 6),
            Self::Plate48 => (6, 8),
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
    #[serde(rename = "cell_culture_plate")]
    CellCulturePlate,
}

impl RackPhysicalTemplateKind {
    pub fn as_str(self) -> &'static str {
        match self {
            Self::StoragePcrTubeRack => "storage_pcr_tube_rack",
            Self::PipettingPcrTubeRack => "pipetting_pcr_tube_rack",
            Self::CellCulturePlate => "cell_culture_plate",
        }
    }
}

#[derive(Debug, Clone, Copy, Serialize, Deserialize, PartialEq, Eq)]
#[serde(rename_all = "snake_case")]
/// Physical handling intent for printable carrier exports.
pub enum RackPhysicalTemplateFamily {
    Storage,
    Pipetting,
    CellCulture,
}

impl RackPhysicalTemplateFamily {
    pub fn as_str(self) -> &'static str {
        match self {
            Self::Storage => "storage",
            Self::Pipetting => "pipetting",
            Self::CellCulture => "cell_culture",
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
/// Overlay x-axis layout for multi-query reference-centered dotplots.
pub enum DotplotOverlayXAxisMode {
    #[default]
    PercentLength,
    LeftAlignedBp,
    RightAlignedBp,
    SharedExonAnchor,
    QueryAnchorBp,
}

impl DotplotOverlayXAxisMode {
    pub fn as_str(self) -> &'static str {
        match self {
            Self::PercentLength => "percent_length",
            Self::LeftAlignedBp => "left_aligned_bp",
            Self::RightAlignedBp => "right_aligned_bp",
            Self::SharedExonAnchor => "shared_exon_anchor",
            Self::QueryAnchorBp => "query_anchor_bp",
        }
    }

    pub fn ui_label(self) -> &'static str {
        match self {
            Self::PercentLength => "% transcript",
            Self::LeftAlignedBp => "bp left-aligned",
            Self::RightAlignedBp => "bp right-aligned",
            Self::SharedExonAnchor => "shared exon anchored",
            Self::QueryAnchorBp => "manual/domain anchor",
        }
    }

    pub fn axis_label(self) -> &'static str {
        match self {
            Self::PercentLength => "x: transcript length (%)",
            Self::LeftAlignedBp => "x: isoform query bp (left-aligned)",
            Self::RightAlignedBp => "x: isoform query bp (right-aligned)",
            Self::SharedExonAnchor => "x: isoform query bp (shared exon anchored)",
            Self::QueryAnchorBp => "x: isoform query bp (manual/domain anchor)",
        }
    }

    pub fn plot_query_span_bp(
        self,
        max_query_span_bp: usize,
        average_query_span_bp: usize,
    ) -> usize {
        match self {
            Self::PercentLength => average_query_span_bp.max(1),
            Self::LeftAlignedBp
            | Self::RightAlignedBp
            | Self::SharedExonAnchor
            | Self::QueryAnchorBp => max_query_span_bp.max(1),
        }
    }

    pub fn axis_edge_labels(self, max_query_span_bp: usize) -> (String, String) {
        match self {
            Self::PercentLength => ("0%".to_string(), "100%".to_string()),
            Self::LeftAlignedBp
            | Self::RightAlignedBp
            | Self::SharedExonAnchor
            | Self::QueryAnchorBp => ("1".to_string(), max_query_span_bp.max(1).to_string()),
        }
    }

    pub fn point_fraction(
        self,
        point_x_0based: usize,
        span_start_0based: usize,
        span_end_0based: usize,
        max_query_span_bp: usize,
    ) -> f32 {
        let query_span_bp = span_end_0based.saturating_sub(span_start_0based).max(1);
        let query_span_max = query_span_bp.saturating_sub(1).max(1);
        let max_query_span_max = max_query_span_bp.max(1).saturating_sub(1).max(1);
        let query_local_bp = point_x_0based
            .saturating_sub(span_start_0based)
            .min(query_span_max);
        match self {
            Self::PercentLength => (query_local_bp as f32 / query_span_max as f32).clamp(0.0, 1.0),
            Self::LeftAlignedBp => {
                (query_local_bp as f32 / max_query_span_max as f32).clamp(0.0, 1.0)
            }
            Self::RightAlignedBp => {
                let offset_bp = max_query_span_bp.saturating_sub(query_span_bp);
                ((offset_bp.saturating_add(query_local_bp)) as f32 / max_query_span_max as f32)
                    .clamp(0.0, 1.0)
            }
            Self::SharedExonAnchor => {
                (query_local_bp as f32 / max_query_span_max as f32).clamp(0.0, 1.0)
            }
            Self::QueryAnchorBp => {
                (query_local_bp as f32 / max_query_span_max as f32).clamp(0.0, 1.0)
            }
        }
    }

    pub fn query_coordinate_at_fraction(
        self,
        axis_fraction: f32,
        span_start_0based: usize,
        span_end_0based: usize,
        max_query_span_bp: usize,
    ) -> Option<usize> {
        let query_span_bp = span_end_0based.saturating_sub(span_start_0based).max(1);
        let query_span_max = query_span_bp.saturating_sub(1).max(1);
        let max_query_span_max = max_query_span_bp.max(1).saturating_sub(1).max(1);
        let axis_fraction = axis_fraction.clamp(0.0, 1.0);
        match self {
            Self::PercentLength => Some(
                span_start_0based
                    .saturating_add((axis_fraction * query_span_max as f32).round() as usize),
            ),
            Self::LeftAlignedBp => {
                let global_bp = (axis_fraction * max_query_span_max as f32).round() as usize;
                if global_bp > query_span_max {
                    None
                } else {
                    Some(span_start_0based.saturating_add(global_bp))
                }
            }
            Self::RightAlignedBp => {
                let offset_bp = max_query_span_bp.saturating_sub(query_span_bp);
                let global_bp = (axis_fraction * max_query_span_max as f32).round() as usize;
                if global_bp < offset_bp || global_bp > offset_bp.saturating_add(query_span_max) {
                    None
                } else {
                    Some(span_start_0based.saturating_add(global_bp - offset_bp))
                }
            }
            Self::SharedExonAnchor => Some(
                span_start_0based
                    .saturating_add((axis_fraction * max_query_span_max as f32).round() as usize),
            ),
            Self::QueryAnchorBp => Some(
                span_start_0based
                    .saturating_add((axis_fraction * max_query_span_max as f32).round() as usize),
            ),
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
    AllOverlappingAnyStrand,
    TargetGroupAnyStrand,
    AllOverlappingTargetStrand,
    TargetGroupTargetStrand,
}

impl SplicingScopePreset {
    pub fn as_str(self) -> &'static str {
        match self {
            Self::AllOverlappingAnyStrand => "all_overlapping_any_strand",
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
#[serde(default)]
#[derive(Default)]
pub struct ProteinFeatureFilter {
    pub include_feature_keys: Vec<String>,
    pub exclude_feature_keys: Vec<String>,
}

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq)]
#[serde(rename_all = "snake_case")]
pub enum FeatureExpertTarget {
    #[serde(alias = "TfbsFeature")]
    TfbsFeature { feature_id: usize },
    #[serde(alias = "RestrictionSite")]
    RestrictionSite {
        cut_pos_1based: usize,
        #[serde(default)]
        enzyme: Option<String>,
        #[serde(default)]
        recognition_start_1based: Option<usize>,
        #[serde(default)]
        recognition_end_1based: Option<usize>,
    },
    #[serde(alias = "SplicingFeature")]
    SplicingFeature {
        feature_id: usize,
        #[serde(default)]
        scope: SplicingScopePreset,
    },
    #[serde(alias = "IsoformArchitecture")]
    IsoformArchitecture { panel_id: String },
    #[serde(alias = "ProteinComparison")]
    ProteinComparison {
        #[serde(default)]
        transcript_id_filter: Option<String>,
        #[serde(default)]
        protein_feature_filter: ProteinFeatureFilter,
        #[serde(default, skip_serializing_if = "Option::is_none")]
        external_source: Option<ProteinExternalOpinionSource>,
        #[serde(default, skip_serializing_if = "Option::is_none")]
        external_entry_id: Option<String>,
    },
    #[serde(alias = "UniprotProjection")]
    UniprotProjection {
        projection_id: String,
        #[serde(default)]
        protein_feature_filter: ProteinFeatureFilter,
    },
}

impl FeatureExpertTarget {
    pub fn uniprot_projection(projection_id: impl Into<String>) -> Self {
        Self::UniprotProjection {
            projection_id: projection_id.into(),
            protein_feature_filter: ProteinFeatureFilter::default(),
        }
    }

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
            Self::ProteinComparison {
                transcript_id_filter,
                protein_feature_filter,
                external_source,
                external_entry_id,
            } => {
                let mut out = "protein comparison".to_string();
                if let Some(transcript_id_filter) = transcript_id_filter
                    .as_deref()
                    .map(str::trim)
                    .filter(|v| !v.is_empty())
                {
                    out.push_str(&format!(" transcript={transcript_id_filter}"));
                }
                if !protein_feature_filter.include_feature_keys.is_empty() {
                    out.push_str(&format!(
                        " include={}",
                        protein_feature_filter.include_feature_keys.join(",")
                    ));
                }
                if !protein_feature_filter.exclude_feature_keys.is_empty() {
                    out.push_str(&format!(
                        " exclude={}",
                        protein_feature_filter.exclude_feature_keys.join(",")
                    ));
                }
                if let (Some(external_source), Some(external_entry_id)) =
                    (external_source, external_entry_id.as_deref())
                {
                    out.push_str(&format!(
                        " external={}:'{}'",
                        external_source.as_str(),
                        external_entry_id
                    ));
                }
                out
            }
            Self::UniprotProjection {
                projection_id,
                protein_feature_filter,
            } => {
                let mut out = format!("UniProt projection '{projection_id}'");
                if !protein_feature_filter.include_feature_keys.is_empty() {
                    out.push_str(&format!(
                        " include={}",
                        protein_feature_filter.include_feature_keys.join(",")
                    ));
                }
                if !protein_feature_filter.exclude_feature_keys.is_empty() {
                    out.push_str(&format!(
                        " exclude={}",
                        protein_feature_filter.exclude_feature_keys.join(",")
                    ));
                }
                out
            }
        }
    }
}

#[cfg(test)]
mod feature_expert_target_tests {
    use super::{FeatureExpertTarget, SplicingScopePreset};

    #[test]
    fn feature_expert_target_accepts_legacy_pascal_case_tags() {
        let target: FeatureExpertTarget = serde_json::from_str(
            r#"{"SplicingFeature":{"feature_id":2,"scope":"all_overlapping_any_strand"}}"#,
        )
        .expect("deserialize legacy splicing feature target");
        assert_eq!(
            target,
            FeatureExpertTarget::SplicingFeature {
                feature_id: 2,
                scope: SplicingScopePreset::AllOverlappingAnyStrand,
            }
        );

        let target: FeatureExpertTarget = serde_json::from_str(
            r#"{"UniprotProjection":{"projection_id":"tp53_uniprot_p04637"}}"#,
        )
        .expect("deserialize legacy UniProt projection target");
        assert_eq!(
            target,
            FeatureExpertTarget::uniprot_projection("tp53_uniprot_p04637")
        );
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

pub const EXON_SKIP_SELECTION_PLAN_SCHEMA: &str = "gentle.exon_skip_selection_plan.v1";
pub const EXON_SKIP_MATERIALIZATION_SCHEMA: &str = "gentle.exon_skip_materialization.v1";

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum ExonSkipSelectionCriterion {
    ManualExonIds {
        candidate_ids: Vec<String>,
    },
    ExplicitIntervals {
        intervals_1based: Vec<SplicingRange>,
    },
    CurrentMapSelection {
        start_1based: usize,
        end_1based: usize,
    },
    FeatureOverlap {
        query: SequenceFeatureQuery,
    },
    LengthMod3 {
        values: Vec<u8>,
    },
    CodingMod3 {
        values: Vec<u8>,
    },
    CodingContext {
        contexts: Vec<String>,
    },
    CdsPhaseEntryKind {
        kinds: Vec<String>,
    },
    ReasoningCandidateIds {
        source_id: String,
        candidate_ids: Vec<String>,
    },
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
pub struct ExonSkipCandidateExon {
    pub candidate_id: String,
    pub ordinal: usize,
    pub start_1based: usize,
    pub end_1based: usize,
    pub length_bp: usize,
    pub length_mod3: usize,
    pub frame_neutral_length: bool,
    pub coding_skip_bp: usize,
    pub coding_skip_mod3: usize,
    pub frame_neutral_coding_skip: bool,
    #[serde(default = "default_exon_skip_coding_context")]
    pub coding_context: String,
    pub support_transcript_count: usize,
    pub support_transcript_total: usize,
    pub support_fraction: f64,
    pub constitutive: bool,
    pub transcript_exon_count: usize,
    #[serde(default = "default_exon_skip_transcript_position")]
    pub transcript_position: String,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub upstream_intron_bp: Option<usize>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub downstream_intron_bp: Option<usize>,
    pub present_in_base_transcript: bool,
    pub cds_overlap: bool,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub left_cds_phase: Option<u8>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub right_cds_phase: Option<u8>,
    #[serde(default = "default_exon_skip_phase_entry_kind")]
    pub cds_phase_entry_kind: String,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub cds_phase_warning: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub coding_frame_note: Option<String>,
    pub selected: bool,
    pub selection_sources: Vec<String>,
    pub matched_feature_ids: Vec<usize>,
    pub rationale: Vec<String>,
}

fn default_exon_skip_phase_entry_kind() -> String {
    "unavailable".to_string()
}

fn default_exon_skip_coding_context() -> String {
    "unknown".to_string()
}

fn default_exon_skip_transcript_position() -> String {
    "unknown".to_string()
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
pub struct ExonSkipSelectionPlan {
    pub schema: String,
    pub plan_id: String,
    pub seq_id: SeqId,
    pub transcript_feature_id: usize,
    pub transcript_id: String,
    pub transcript_label: String,
    pub strand: String,
    pub region_start_1based: usize,
    pub region_end_1based: usize,
    pub criteria: Vec<ExonSkipSelectionCriterion>,
    pub candidate_exons: Vec<ExonSkipCandidateExon>,
    pub selected_candidate_ids: Vec<String>,
    pub warnings: Vec<String>,
    pub messages: Vec<String>,
}

#[derive(Debug, Clone, Copy, Serialize, Deserialize, PartialEq, Eq, Hash)]
#[serde(rename_all = "snake_case")]
#[derive(Default)]
pub enum ExonSkipReturnKind {
    #[default]
    Genbank,
    CdnaFasta,
    AminoAcidSequence,
    AminoAcidFasta,
}

impl ExonSkipReturnKind {
    pub fn as_str(self) -> &'static str {
        match self {
            Self::Genbank => "genbank",
            Self::CdnaFasta => "cdna_fasta",
            Self::AminoAcidSequence => "amino_acid_sequence",
            Self::AminoAcidFasta => "amino_acid_fasta",
        }
    }
}

#[derive(Debug, Clone, Serialize, Deserialize, Default, PartialEq, Eq)]
#[serde(default)]
pub struct ExonSkipReturnPayload {
    pub kind: ExonSkipReturnKind,
    pub available: bool,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub seq_id: Option<SeqId>,
    pub label: String,
    pub mime_type: String,
    pub text: String,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub message: Option<String>,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default, PartialEq, Eq)]
#[serde(default)]
pub struct ExonSkipMaterializationReport {
    pub schema: String,
    pub plan_id: String,
    pub source_seq_id: SeqId,
    pub transcript_feature_id: usize,
    pub skipped_candidate_ids: Vec<String>,
    pub retained_exon_count: usize,
    pub skipped_exon_count: usize,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub genomic_seq_id: Option<SeqId>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub cdna_seq_id: Option<SeqId>,
    #[serde(skip_serializing_if = "Vec::is_empty")]
    pub requested_returns: Vec<ExonSkipReturnKind>,
    #[serde(skip_serializing_if = "Vec::is_empty")]
    pub return_payloads: Vec<ExonSkipReturnPayload>,
    pub warnings: Vec<String>,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
/// One transcript-derived PCR/qPCR product materialized as a reusable sequence.
pub struct CdnaAssayMaterializedProductRow {
    pub product_seq_id: SeqId,
    pub transcript_id: String,
    pub transcript_feature_id: usize,
    pub product_index: usize,
    pub amplicon_length_bp: usize,
    pub amplicon_start_0based: usize,
    pub amplicon_end_0based_exclusive: usize,
    pub probe_supported: bool,
    #[serde(skip_serializing_if = "Vec::is_empty")]
    pub covered_junction_labels: Vec<String>,
    #[serde(skip_serializing_if = "String::is_empty")]
    pub genomic_carryover_risk: String,
    #[serde(skip_serializing_if = "String::is_empty")]
    pub genomic_carryover_rationale: String,
    pub created: bool,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
/// Text-friendly row describing one rendered PCR/qPCR product-gel band.
pub struct CdnaAssayProductGelBandRow {
    pub lane_name: String,
    pub band_index: usize,
    pub apparent_bp: usize,
    pub min_bp: usize,
    pub max_bp: usize,
    pub product_count: usize,
    #[serde(skip_serializing_if = "Vec::is_empty")]
    pub labels: Vec<String>,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
/// Summary of optional cDNA PCR/qPCR product materialization into concrete
/// GENtle sequence entries, a reusable product container, and an optional gel.
pub struct CdnaAssayProductMaterialization {
    pub schema: String,
    pub assay_kind: String,
    pub source_seq_id: String,
    pub source_feature_id: usize,
    pub group_label: String,
    pub product_count: usize,
    #[serde(skip_serializing_if = "Vec::is_empty")]
    pub product_seq_ids: Vec<SeqId>,
    #[serde(skip_serializing_if = "Vec::is_empty")]
    pub created_product_seq_ids: Vec<SeqId>,
    #[serde(skip_serializing_if = "Vec::is_empty")]
    pub reused_product_seq_ids: Vec<SeqId>,
    #[serde(skip_serializing_if = "Vec::is_empty")]
    pub product_rows: Vec<CdnaAssayMaterializedProductRow>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub container_id: Option<ContainerId>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub container_kind: Option<ContainerKind>,
    pub container_created: bool,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub product_gel_svg_path: Option<String>,
    #[serde(skip_serializing_if = "Vec::is_empty")]
    pub gel_band_rows: Vec<CdnaAssayProductGelBandRow>,
    #[serde(skip_serializing_if = "Vec::is_empty")]
    pub gel_summary_lines: Vec<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub product_output_prefix: Option<String>,
    pub idempotent_reuse: bool,
    #[serde(skip_serializing_if = "Vec::is_empty")]
    pub warnings: Vec<String>,
}

#[cfg(test)]
mod exon_skip_return_contract_tests {
    use super::{
        CdnaAssayMaterializedProductRow, CdnaAssayProductGelBandRow,
        CdnaAssayProductMaterialization, ContainerKind, ExonSkipReturnKind, ExonSkipReturnPayload,
    };

    #[test]
    fn exon_skip_return_payloads_keep_stable_public_spellings() {
        assert_eq!(ExonSkipReturnKind::Genbank.as_str(), "genbank");
        assert_eq!(ExonSkipReturnKind::CdnaFasta.as_str(), "cdna_fasta");
        assert_eq!(
            ExonSkipReturnKind::AminoAcidSequence.as_str(),
            "amino_acid_sequence"
        );
        assert_eq!(
            ExonSkipReturnKind::AminoAcidFasta.as_str(),
            "amino_acid_fasta"
        );

        let payload = ExonSkipReturnPayload {
            kind: ExonSkipReturnKind::AminoAcidSequence,
            available: true,
            seq_id: None,
            label: "skipped isoform protein".to_string(),
            mime_type: "text/plain".to_string(),
            text: "MK".to_string(),
            message: None,
        };
        let json = serde_json::to_string(&payload).expect("serialize return payload");
        assert!(json.contains(r#""kind":"amino_acid_sequence""#));
    }

    #[test]
    fn cdna_assay_product_materialization_contract_round_trips_idempotency_fields() {
        let payload = CdnaAssayProductMaterialization {
            schema: "gentle.cdna_assay_product_materialization.v1".to_string(),
            assay_kind: "pcr".to_string(),
            source_seq_id: "source".to_string(),
            source_feature_id: 7,
            group_label: "GENE1".to_string(),
            product_count: 1,
            product_seq_ids: vec!["gene1_pcr_TX1_p1_42bp".to_string()],
            reused_product_seq_ids: vec!["gene1_pcr_TX1_p1_42bp".to_string()],
            product_rows: vec![CdnaAssayMaterializedProductRow {
                product_seq_id: "gene1_pcr_TX1_p1_42bp".to_string(),
                transcript_id: "TX1".to_string(),
                transcript_feature_id: 3,
                product_index: 1,
                amplicon_length_bp: 42,
                amplicon_start_0based: 10,
                amplicon_end_0based_exclusive: 52,
                probe_supported: false,
                genomic_carryover_risk: "low".to_string(),
                created: false,
                ..CdnaAssayMaterializedProductRow::default()
            }],
            container_id: Some("container-1".to_string()),
            container_kind: Some(ContainerKind::Singleton),
            gel_band_rows: vec![CdnaAssayProductGelBandRow {
                lane_name: "cDNA PCR products (GENE1)".to_string(),
                band_index: 1,
                apparent_bp: 42,
                min_bp: 42,
                max_bp: 42,
                product_count: 1,
                labels: vec!["gene1_pcr_TX1_p1_42bp (42 bp)".to_string()],
            }],
            gel_summary_lines: vec![
                "Product gel lane 'cDNA PCR products (GENE1)' has 1 band(s).".to_string(),
            ],
            product_output_prefix: Some("gene1_pcr".to_string()),
            idempotent_reuse: true,
            ..CdnaAssayProductMaterialization::default()
        };
        let json =
            serde_json::to_string(&payload).expect("serialize product materialization contract");
        assert!(json.contains(r#""reused_product_seq_ids":["gene1_pcr_TX1_p1_42bp"]"#));
        assert!(json.contains(r#""gel_summary_lines":["#));
        let round_tripped: CdnaAssayProductMaterialization =
            serde_json::from_str(&json).expect("round-trip product materialization contract");
        assert!(round_tripped.idempotent_reuse);
        assert!(!round_tripped.product_rows[0].created);
        assert_eq!(round_tripped.gel_band_rows[0].apparent_bp, 42);
    }
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SplicingBoundaryMarker {
    pub transcript_feature_id: usize,
    pub transcript_id: String,
    pub side: String,
    pub position_1based: usize,
    pub motif_2bp: String,
    pub canonical: bool,
    #[serde(default)]
    pub canonical_pair: bool,
    #[serde(default)]
    pub partner_position_1based: usize,
    #[serde(default)]
    pub paired_motif_signature: String,
    #[serde(default)]
    pub motif_class: String,
    #[serde(default)]
    pub annotation: String,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SplicingIntronSignal {
    pub transcript_feature_id: usize,
    pub transcript_id: String,
    pub donor_position_1based: usize,
    pub acceptor_position_1based: usize,
    pub intron_length_bp: usize,
    #[serde(default)]
    pub branchpoint_position_1based: Option<usize>,
    #[serde(default)]
    pub branchpoint_motif: String,
    #[serde(default)]
    pub branchpoint_score: f32,
    #[serde(default)]
    pub branchpoint_annotation: String,
    #[serde(default)]
    pub polypyrimidine_start_1based: Option<usize>,
    #[serde(default)]
    pub polypyrimidine_end_1based: Option<usize>,
    #[serde(default)]
    pub polypyrimidine_fraction: f32,
    #[serde(default)]
    pub polypyrimidine_annotation: String,
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
    #[serde(default)]
    pub scope: SplicingScopePreset,
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
    #[serde(default)]
    pub intron_signals: Vec<SplicingIntronSignal>,
    pub junctions: Vec<SplicingJunctionArc>,
    pub events: Vec<SplicingEventSummary>,
}

#[derive(Debug, Clone, Copy, Serialize, Deserialize, PartialEq, Eq, Default)]
#[serde(rename_all = "snake_case")]
pub enum AttractSpeciesMatchMode {
    #[default]
    ExactOrganism,
    FallbackAllCompatible,
}

impl AttractSpeciesMatchMode {
    pub fn as_str(self) -> &'static str {
        match self {
            Self::ExactOrganism => "exact_organism",
            Self::FallbackAllCompatible => "fallback_all_compatible",
        }
    }
}

#[derive(Debug, Clone, Copy, Serialize, Deserialize, PartialEq, Eq, Default)]
#[serde(rename_all = "snake_case")]
pub enum AttractRegionClass {
    #[default]
    ExonBody,
    DonorFlank,
    AcceptorFlank,
    IntronBody,
}

impl AttractRegionClass {
    pub fn as_str(self) -> &'static str {
        match self {
            Self::ExonBody => "exon_body",
            Self::DonorFlank => "donor_flank",
            Self::AcceptorFlank => "acceptor_flank",
            Self::IntronBody => "intron_body",
        }
    }
}

#[derive(Debug, Clone, Copy, Serialize, Deserialize, PartialEq, Eq, Default)]
#[serde(rename_all = "snake_case")]
pub enum AttractPwmMappingPolicy {
    #[default]
    StrictSameLength,
    WindowedSubmatrix,
}

impl AttractPwmMappingPolicy {
    pub fn as_str(self) -> &'static str {
        match self {
            Self::StrictSameLength => "strict_same_length",
            Self::WindowedSubmatrix => "windowed_submatrix",
        }
    }
}

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
#[serde(default)]
pub struct AttractSplicingEvidenceSettings {
    pub scope: SplicingScopePreset,
    pub transcript_strand_only: bool,
    pub boundary_flank_bp: usize,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub requested_organism: Option<String>,
    pub allow_species_fallback: bool,
    pub minimum_quality_score: f64,
    pub minimum_match_quantile: f64,
    pub pwm_mapping_policy: AttractPwmMappingPolicy,
    pub compare_alternate_policy: bool,
}

impl Default for AttractSplicingEvidenceSettings {
    fn default() -> Self {
        Self {
            scope: SplicingScopePreset::TargetGroupTargetStrand,
            transcript_strand_only: true,
            boundary_flank_bp: 25,
            requested_organism: None,
            allow_species_fallback: true,
            minimum_quality_score: 0.0,
            minimum_match_quantile: 0.99,
            pwm_mapping_policy: AttractPwmMappingPolicy::StrictSameLength,
            compare_alternate_policy: false,
        }
    }
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
pub struct AttractSplicingEvidencePolicySummary {
    pub pwm_mapping_policy: AttractPwmMappingPolicy,
    pub unique_rbp_count: usize,
    pub hit_count: usize,
    pub pwm_scored_hit_count: usize,
    pub exact_length_pwm_hit_count: usize,
    pub windowed_pwm_hit_count: usize,
    pub consensus_hit_count: usize,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
pub struct AttractSplicingEvidenceSummaryRow {
    pub gene_name: String,
    pub organism: String,
    pub matrix_id: String,
    pub motif_iupac: String,
    pub model_kind: String,
    pub hit_count: usize,
    pub strongest_score: f64,
    #[serde(default)]
    pub strongest_score_kind: String,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub strongest_score_quantile: Option<f64>,
    pub max_quality_score: f64,
    pub exon_body_hits: usize,
    pub donor_flank_hits: usize,
    pub acceptor_flank_hits: usize,
    pub intron_body_hits: usize,
    #[serde(default)]
    pub exact_length_pwm_hits: usize,
    #[serde(default)]
    pub windowed_pwm_hits: usize,
    #[serde(default)]
    pub consensus_only_hits: usize,
    #[serde(default)]
    pub supporting_transcript_ids: Vec<String>,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
pub struct AttractSplicingEvidenceHitRow {
    pub transcript_feature_id: usize,
    pub transcript_id: String,
    pub transcript_label: String,
    pub transcript_strand: String,
    pub gene_name: String,
    pub organism: String,
    pub matrix_id: String,
    pub motif_iupac: String,
    pub model_kind: String,
    pub region_class: AttractRegionClass,
    pub region_start_1based: usize,
    pub region_end_1based: usize,
    pub region_local_start_1based: usize,
    pub region_local_end_1based: usize,
    pub matched_sequence: String,
    pub match_score: f64,
    #[serde(default)]
    pub match_score_kind: String,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub match_score_quantile: Option<f64>,
    pub quality_score: f64,
    #[serde(default)]
    pub exact_species_match: bool,
    #[serde(default)]
    pub pwm_mapping_status: String,
    #[serde(default)]
    pub mapping_policy_used: String,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub pfm_subwindow_start_1based: Option<usize>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub pfm_subwindow_end_1based: Option<usize>,
    #[serde(default)]
    pub warnings: Vec<String>,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
pub struct AttractSplicingEvidenceView {
    pub schema: String,
    pub seq_id: String,
    pub target_feature_id: usize,
    pub scope: SplicingScopePreset,
    pub group_label: String,
    pub target_strand: String,
    pub settings: AttractSplicingEvidenceSettings,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub requested_organism: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub resolved_organism: Option<String>,
    pub species_match_mode: AttractSpeciesMatchMode,
    pub scanned_transcript_count: usize,
    pub scanned_window_count: usize,
    pub unique_rbp_count: usize,
    pub hit_count: usize,
    pub pwm_scored_hit_count: usize,
    #[serde(default)]
    pub exact_length_pwm_hit_count: usize,
    #[serde(default)]
    pub windowed_pwm_hit_count: usize,
    pub consensus_hit_count: usize,
    pub active_resource_source: String,
    pub active_resource_item_count: usize,
    pub active_resource_pwm_row_count: usize,
    pub active_resource_consensus_only_row_count: usize,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub active_resource_fingerprint: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub alternate_policy_summary: Option<AttractSplicingEvidencePolicySummary>,
    #[serde(default)]
    pub summary_rows: Vec<AttractSplicingEvidenceSummaryRow>,
    #[serde(default)]
    pub hit_rows: Vec<AttractSplicingEvidenceHitRow>,
    #[serde(default)]
    pub warnings: Vec<String>,
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
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub comparison: Option<TranscriptProteinComparison>,
}

/// Optional expression/readout matrix rendered beside isoform architecture lanes.
///
/// `rows[].values` are aligned to `sample_labels`; `None` represents a missing
/// isoform/sample cell rather than a zero value.
#[derive(Debug, Clone, Serialize, Deserialize, Default)]
pub struct IsoformExpressionMatrix {
    pub sample_labels: Vec<String>,
    pub rows: Vec<IsoformExpressionRow>,
    #[serde(default)]
    pub warnings: Vec<String>,
}

/// One isoform row in an expression matrix.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct IsoformExpressionRow {
    pub isoform_id: String,
    pub values: Vec<Option<f64>>,
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
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub expression_matrix: Option<IsoformExpressionMatrix>,
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
    #[serde(default)]
    pub tooltip_lines: Vec<String>,
    pub instruction: String,
}

impl RestrictionSiteExpertView {
    pub fn enzyme_display_label(&self) -> String {
        self.selected_enzyme
            .clone()
            .or_else(|| (!self.enzyme_names.is_empty()).then(|| self.enzyme_names.join(", ")))
            .unwrap_or_else(|| "Restriction site".to_string())
    }

    pub fn paired_cut_pos_1based_effective(&self) -> usize {
        if self.paired_cut_pos_1based == 0 {
            self.cut_pos_1based
        } else {
            self.paired_cut_pos_1based
        }
    }

    pub fn paired_cut_index_0based_effective(&self) -> usize {
        if self.paired_cut_pos_1based == 0 && self.paired_cut_index_0based == 0 {
            self.cut_index_0based
        } else {
            self.paired_cut_index_0based
        }
    }

    pub fn site_count_label(&self) -> String {
        match self.number_of_cuts_for_enzyme {
            1 => "1 site".to_string(),
            count => format!("{count} sites"),
        }
    }

    pub fn cut_position_label(&self) -> String {
        let paired = self.paired_cut_pos_1based_effective();
        if paired == self.cut_pos_1based {
            format!("cut at {} bp", self.cut_pos_1based)
        } else {
            format!("cuts at {}|{} bp", self.cut_pos_1based, paired)
        }
    }

    pub fn geometry_display_label(&self) -> String {
        match self.end_geometry.as_str() {
            "5prime_overhang" => format!(
                "5' overhang ({} bp)",
                self.overlap_bp.unwrap_or(0).unsigned_abs()
            ),
            "3prime_overhang" => format!(
                "3' overhang ({} bp)",
                self.overlap_bp.unwrap_or(0).unsigned_abs()
            ),
            _ => match self.overlap_bp {
                Some(value) if value > 0 => format!("5' overhang ({} bp)", value as usize),
                Some(value) if value < 0 => format!("3' overhang ({} bp)", value.unsigned_abs()),
                _ => "blunt".to_string(),
            },
        }
    }

    pub fn marked_top_sequence(&self) -> String {
        sequence_with_cut_marker(&self.site_sequence, self.cut_index_0based)
    }

    pub fn marked_bottom_sequence(&self) -> String {
        sequence_with_cut_marker(
            &self.site_sequence_complement,
            self.paired_cut_index_0based_effective(),
        )
    }

    pub fn tooltip_summary_lines(&self) -> Vec<String> {
        if !self.tooltip_lines.is_empty() {
            return self.tooltip_lines.clone();
        }
        self.computed_tooltip_summary_lines()
    }

    pub fn computed_tooltip_summary_lines(&self) -> Vec<String> {
        let mut lines = vec![
            format!(
                "{} | {} | {}",
                self.enzyme_display_label(),
                self.site_count_label(),
                self.geometry_display_label()
            ),
            format!(
                "{} | recognition {}..{}",
                self.cut_position_label(),
                self.recognition_start_1based,
                self.recognition_end_1based
            ),
        ];
        if let Some(iupac) = self.recognition_iupac.as_deref() {
            lines.push(format!("recognition_iupac={iupac}"));
        }
        lines.push(format!("5' {} 3'", self.marked_top_sequence()));
        lines.push(format!("3' {} 5'", self.marked_bottom_sequence()));
        if let Some(note) = self.enzyme_note.as_deref() {
            lines.push(format!("note={note}"));
        }
        if let Some(url) = self.rebase_url.as_deref() {
            lines.push(format!("REBASE: {url}"));
        }
        lines
    }
}

fn sequence_with_cut_marker(sequence: &str, index_0based: usize) -> String {
    let chars = sequence.chars().collect::<Vec<_>>();
    let split = index_0based.min(chars.len());
    let left = chars[..split].iter().collect::<String>();
    let right = chars[split..].iter().collect::<String>();
    format!("{left}^{right}")
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

pub const READ_ACQUISITION_REPORT_SCHEMA: &str = "gentle.read_acquisition_report.v1";

#[derive(Debug, Clone, Copy, Serialize, Deserialize, PartialEq, Eq, Hash, Default)]
#[serde(rename_all = "snake_case")]
pub enum ReadAcquisitionAnalysisFormat {
    #[default]
    Fasta,
    Fastq,
}

impl ReadAcquisitionAnalysisFormat {
    pub fn as_str(self) -> &'static str {
        match self {
            Self::Fasta => "fasta",
            Self::Fastq => "fastq",
        }
    }
}

#[derive(Debug, Clone, Copy, Serialize, Deserialize, PartialEq, Eq, Hash, Default)]
#[serde(rename_all = "snake_case")]
pub enum ReadAcquisitionReadLayout {
    #[default]
    SingleEnd,
    PairedEnd,
    SplitSpot,
}

impl ReadAcquisitionReadLayout {
    pub fn as_str(self) -> &'static str {
        match self {
            Self::SingleEnd => "single_end",
            Self::PairedEnd => "paired_end",
            Self::SplitSpot => "split_spot",
        }
    }
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
pub struct ReadAcquisitionManifestRow {
    pub row_number: usize,
    pub sample_id: String,
    pub sra_accession: String,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub sample_name: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub assay_kind: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub read_layout: Option<ReadAcquisitionReadLayout>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub analysis_format: Option<ReadAcquisitionAnalysisFormat>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub note: Option<String>,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
pub struct ReadAcquisitionReadLengthStats {
    pub read_count: usize,
    pub total_bases: usize,
    pub min_length_bp: usize,
    pub max_length_bp: usize,
    pub mean_length_bp: f64,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
pub struct ReadAcquisitionOutputPath {
    pub role: String,
    pub path: String,
    pub analysis_format: ReadAcquisitionAnalysisFormat,
    pub file_size_bytes: u64,
    pub checksum_sha1: String,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub read_length_stats: Option<ReadAcquisitionReadLengthStats>,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
pub struct ReadAcquisitionCommandProvenance {
    pub phase: String,
    pub command: String,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub exit_code: Option<i32>,
    pub stdout_log_path: String,
    pub stderr_log_path: String,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
pub struct ReadAcquisitionRunReport {
    pub sample_id: String,
    pub sra_accession: String,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub sample_name: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub assay_kind: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub note: Option<String>,
    pub resource_key: String,
    pub lifecycle_status: String,
    pub read_layout: ReadAcquisitionReadLayout,
    pub analysis_format: ReadAcquisitionAnalysisFormat,
    pub sra_path: String,
    pub run_dir: String,
    #[serde(default)]
    pub output_paths: Vec<ReadAcquisitionOutputPath>,
    #[serde(default)]
    pub command_provenance: Vec<ReadAcquisitionCommandProvenance>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub current_activity: Option<SharedAssetActivityStatus>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub error: Option<String>,
    #[serde(default)]
    pub warnings: Vec<String>,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
pub struct ReadAcquisitionReport {
    pub schema: String,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub manifest_path: Option<String>,
    pub cache_dir: String,
    pub work_dir: String,
    pub generated_at_unix_ms: u128,
    pub lifecycle_status: String,
    pub sample_count: usize,
    pub ready_count: usize,
    pub running_count: usize,
    pub failed_count: usize,
    pub cancelled_count: usize,
    pub stale_count: usize,
    pub missing_count: usize,
    #[serde(default)]
    pub rows: Vec<ReadAcquisitionRunReport>,
    #[serde(default)]
    pub warnings: Vec<String>,
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
pub enum CutRunReadLayout {
    #[default]
    SingleEnd,
    PairedEnd,
}

impl CutRunReadLayout {
    pub fn as_str(self) -> &'static str {
        match self {
            Self::SingleEnd => "single_end",
            Self::PairedEnd => "paired_end",
        }
    }
}

#[derive(Debug, Clone, Copy, Serialize, Deserialize, PartialEq, Eq, Default)]
#[serde(rename_all = "snake_case")]
pub enum CutRunInputFormat {
    #[default]
    Fasta,
    Fastq,
}

impl CutRunInputFormat {
    pub fn as_str(self) -> &'static str {
        match self {
            Self::Fasta => "fasta",
            Self::Fastq => "fastq",
        }
    }
}

#[derive(Debug, Clone, Copy, Serialize, Deserialize, PartialEq, Eq, Default)]
#[serde(rename_all = "snake_case")]
pub enum CutRunCoverageKind {
    #[default]
    Coverage,
    CutSites,
    Fragments,
}

impl CutRunCoverageKind {
    pub fn as_str(self) -> &'static str {
        match self {
            Self::Coverage => "coverage",
            Self::CutSites => "cut_sites",
            Self::Fragments => "fragments",
        }
    }
}

#[derive(Debug, Clone, Copy, Serialize, Deserialize, PartialEq, Eq, Default)]
#[serde(rename_all = "snake_case")]
pub enum CutRunReadOrientation {
    #[default]
    Forward,
    Reverse,
}

impl CutRunReadOrientation {
    pub fn as_str(self) -> &'static str {
        match self {
            Self::Forward => "forward",
            Self::Reverse => "reverse",
        }
    }
}

#[derive(Debug, Clone, Copy, Serialize, Deserialize, PartialEq, Eq, PartialOrd, Ord, Default)]
#[serde(rename_all = "snake_case")]
pub enum CutRunReadUnitStatus {
    #[default]
    SingleMapped,
    SingleUnmapped,
    ConcordantPair,
    OrphanR1,
    OrphanR2,
    R1OnlyMapped,
    R2OnlyMapped,
    DiscordantPair,
    PairUnmapped,
}

impl CutRunReadUnitStatus {
    pub fn as_str(self) -> &'static str {
        match self {
            Self::SingleMapped => "single_mapped",
            Self::SingleUnmapped => "single_unmapped",
            Self::ConcordantPair => "concordant_pair",
            Self::OrphanR1 => "orphan_r1",
            Self::OrphanR2 => "orphan_r2",
            Self::R1OnlyMapped => "r1_only_mapped",
            Self::R2OnlyMapped => "r2_only_mapped",
            Self::DiscordantPair => "discordant_pair",
            Self::PairUnmapped => "pair_unmapped",
        }
    }
}

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
#[serde(default)]
pub struct CutRunSeedFilterConfig {
    pub kmer_len: usize,
    pub min_seed_matches: usize,
}

impl Default for CutRunSeedFilterConfig {
    fn default() -> Self {
        Self {
            kmer_len: 12,
            min_seed_matches: 1,
        }
    }
}

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
#[serde(default)]
pub struct CutRunAlignConfig {
    pub max_mismatches: usize,
    pub min_identity_fraction: f64,
    pub max_fragment_span_bp: usize,
}

impl Default for CutRunAlignConfig {
    fn default() -> Self {
        Self {
            max_mismatches: 4,
            min_identity_fraction: 0.9,
            max_fragment_span_bp: 800,
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
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub cause_chain: Vec<String>,
}

impl EngineError {
    pub fn new(code: ErrorCode, message: impl Into<String>) -> Self {
        Self {
            code,
            message: message.into(),
            cause_chain: vec![],
        }
    }

    pub fn invalid_input(message: impl Into<String>) -> Self {
        Self::new(ErrorCode::InvalidInput, message)
    }

    pub fn internal(message: impl Into<String>) -> Self {
        Self::new(ErrorCode::Internal, message)
    }

    pub fn with_cause(mut self, source: impl fmt::Display) -> Self {
        self.cause_chain.push(source.to_string());
        self
    }

    pub fn portable_payload(&self, extra_cause_chain: Vec<String>) -> Value {
        let mut cause_chain = self.cause_chain.clone();
        cause_chain.extend(extra_cause_chain);
        json!({
            "code": self.code,
            "message": self.message,
            "cause_chain": cause_chain
        })
    }
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
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub transcript_feature_id: Option<usize>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub query_anchor_0based: Option<usize>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub query_anchor_label: Option<String>,
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
    pub kind: String,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub strand: Option<String>,
    pub lane: usize,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub color_rgb: Option<[u8; 3]>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub detail: Option<String>,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
pub struct DotplotReferenceAnnotationTrack {
    pub seq_id: String,
    pub label: String,
    pub interval_count: usize,
    pub intervals: Vec<DotplotReferenceAnnotationInterval>,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default, PartialEq, Eq)]
#[serde(default)]
pub struct DotplotOverlayAnchorExonRef {
    pub start_1based: usize,
    pub end_1based: usize,
}

impl DotplotOverlayAnchorExonRef {
    pub fn token(&self) -> String {
        format!("{}..{}", self.start_1based, self.end_1based)
    }

    pub fn parse(raw: &str) -> Result<Self, String> {
        let trimmed = raw.trim();
        let (left, right) = trimmed.split_once("..").ok_or_else(|| {
            format!(
                "Invalid shared-exon anchor '{raw}'; expected START..END with 1-based inclusive coordinates"
            )
        })?;
        let start_1based = left.trim().parse::<usize>().map_err(|e| {
            format!(
                "Invalid shared-exon anchor start '{}' in '{}': {e}",
                left.trim(),
                raw
            )
        })?;
        let end_1based = right.trim().parse::<usize>().map_err(|e| {
            format!(
                "Invalid shared-exon anchor end '{}' in '{}': {e}",
                right.trim(),
                raw
            )
        })?;
        if start_1based == 0 || end_1based == 0 {
            return Err(format!(
                "Invalid shared-exon anchor '{raw}'; coordinates must be >= 1"
            ));
        }
        if end_1based < start_1based {
            return Err(format!(
                "Invalid shared-exon anchor '{raw}'; end must be >= start"
            ));
        }
        Ok(Self {
            start_1based,
            end_1based,
        })
    }
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
pub struct DotplotOverlayAnchorSeriesSupport {
    pub series_id: String,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub transcript_feature_id: Option<usize>,
    pub query_start_0based: usize,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
pub struct DotplotOverlayAnchorExon {
    pub exon: DotplotOverlayAnchorExonRef,
    pub support_series_count: usize,
    pub max_query_start_0based: usize,
    pub supporting_series: Vec<DotplotOverlayAnchorSeriesSupport>,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
pub struct DotplotOverlayQuerySpec {
    pub seq_id: String,
    pub label: String,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub transcript_feature_id: Option<usize>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub query_anchor_0based: Option<usize>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub query_anchor_label: Option<String>,
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
    #[serde(default)]
    pub overlay_anchor_exons: Vec<DotplotOverlayAnchorExon>,
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

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct DotplotOverlayResolvedAnchorSeries {
    pub series_index: usize,
    pub shift_bp: usize,
    pub plotted_span_end_0based: usize,
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
                transcript_feature_id: None,
                query_anchor_0based: None,
                query_anchor_label: None,
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
            for interval in &mut annotation.intervals {
                if interval.kind.trim().is_empty() {
                    interval.kind = "exon".to_string();
                }
                if interval.label.trim().is_empty() {
                    interval.label = interval.kind.clone();
                }
            }
        }
        for anchor in &mut self.overlay_anchor_exons {
            anchor.support_series_count = anchor.supporting_series.len();
        }
    }

    pub fn primary_series(&self) -> Option<&DotplotQuerySeries> {
        self.query_series.first()
    }

    pub fn overlay_anchor_by_ref(
        &self,
        exon: &DotplotOverlayAnchorExonRef,
    ) -> Option<&DotplotOverlayAnchorExon> {
        self.overlay_anchor_exons
            .iter()
            .find(|anchor| anchor.exon == *exon)
    }

    pub fn resolve_overlay_anchor_series(
        &self,
        exon: &DotplotOverlayAnchorExonRef,
    ) -> Vec<DotplotOverlayResolvedAnchorSeries> {
        let Some(anchor) = self.overlay_anchor_by_ref(exon) else {
            return vec![];
        };
        let mut resolved = vec![];
        for support in &anchor.supporting_series {
            let Some(series_index) = self
                .query_series
                .iter()
                .position(|series| series.series_id == support.series_id)
            else {
                continue;
            };
            let Some(series) = self.query_series.get(series_index) else {
                continue;
            };
            let shift_bp = anchor
                .max_query_start_0based
                .saturating_sub(support.query_start_0based);
            let series_span_bp = series
                .span_end_0based
                .saturating_sub(series.span_start_0based)
                .max(1);
            resolved.push(DotplotOverlayResolvedAnchorSeries {
                series_index,
                shift_bp,
                plotted_span_end_0based: shift_bp.saturating_add(series_span_bp),
            });
        }
        resolved
    }

    pub fn resolve_query_anchor_series(&self) -> Vec<DotplotOverlayResolvedAnchorSeries> {
        let max_query_anchor_0based = self
            .query_series
            .iter()
            .filter_map(|series| {
                series
                    .query_anchor_0based
                    .map(|anchor| anchor.saturating_sub(series.span_start_0based))
            })
            .max();
        let Some(max_query_anchor_0based) = max_query_anchor_0based else {
            return vec![];
        };
        let mut resolved = vec![];
        for (series_index, series) in self.query_series.iter().enumerate() {
            let Some(anchor_0based) = series.query_anchor_0based else {
                continue;
            };
            let query_anchor_local_0based = anchor_0based.saturating_sub(series.span_start_0based);
            let shift_bp = max_query_anchor_0based.saturating_sub(query_anchor_local_0based);
            let series_span_bp = series
                .span_end_0based
                .saturating_sub(series.span_start_0based)
                .max(1);
            resolved.push(DotplotOverlayResolvedAnchorSeries {
                series_index,
                shift_bp,
                plotted_span_end_0based: shift_bp.saturating_add(series_span_bp),
            });
        }
        resolved
    }

    pub fn query_anchor_label(&self) -> Option<&str> {
        let mut labels = self
            .query_series
            .iter()
            .filter_map(|series| series.query_anchor_label.as_deref())
            .map(str::trim)
            .filter(|label| !label.is_empty());
        let first = labels.next()?;
        if labels.all(|label| label.eq_ignore_ascii_case(first)) {
            Some(first)
        } else {
            None
        }
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

#[derive(Debug, Clone, Copy, Serialize, Deserialize, PartialEq, Eq, Default)]
#[serde(rename_all = "snake_case")]
pub enum FeatureBedCoordinateMode {
    #[default]
    Auto,
    Local,
    Genomic,
}

impl FeatureBedCoordinateMode {
    pub fn as_str(self) -> &'static str {
        match self {
            Self::Auto => "auto",
            Self::Local => "local",
            Self::Genomic => "genomic",
        }
    }
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
pub struct SequenceFeatureBedExportReport {
    pub schema: String,
    pub seq_id: SeqId,
    pub path: String,
    pub coordinate_mode: String,
    pub include_restriction_sites: bool,
    pub restriction_enzyme_filters: Vec<String>,
    pub bed_columns: Vec<String>,
    pub matched_sequence_feature_count: usize,
    pub matched_restriction_site_count: usize,
    pub matched_row_count: usize,
    pub exportable_row_count: usize,
    pub exported_row_count: usize,
    pub offset: usize,
    pub limit: Option<usize>,
    pub local_coordinate_row_count: usize,
    pub genomic_coordinate_row_count: usize,
    pub skipped_missing_genomic_coordinates: usize,
    pub query: SequenceFeatureQuery,
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
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub task_kind: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub stage_label: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub detail: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub stage_percent: Option<f64>,
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
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub capability_registry: Vec<CapabilityDescriptor>,
}

/// Adapter families that can project a shared capability descriptor.
#[derive(Debug, Clone, Copy, Serialize, Deserialize, PartialEq, Eq, PartialOrd, Ord, Hash)]
#[serde(rename_all = "snake_case")]
pub enum CapabilityAdapter {
    Gui,
    Cli,
    Mcp,
    Js,
    Lua,
    Clawbio,
}

impl CapabilityAdapter {
    pub fn as_str(self) -> &'static str {
        match self {
            Self::Gui => "gui",
            Self::Cli => "cli",
            Self::Mcp => "mcp",
            Self::Js => "js",
            Self::Lua => "lua",
            Self::Clawbio => "clawbio",
        }
    }

    pub fn label(self) -> &'static str {
        match self {
            Self::Gui => "GUI",
            Self::Cli => "gentle_cli",
            Self::Mcp => "MCP",
            Self::Js => "JS",
            Self::Lua => "Lua",
            Self::Clawbio => "ClawBio",
        }
    }
}

/// Mutation class consumed by adapter safety boundaries.
#[derive(Debug, Clone, Copy, Serialize, Deserialize, PartialEq, Eq)]
pub enum CapabilityMutation {
    #[serde(rename = "false")]
    ReadOnly,
    #[serde(rename = "true")]
    Mutating,
    #[serde(rename = "external")]
    External,
}

impl CapabilityMutation {
    pub fn as_str(self) -> &'static str {
        match self {
            Self::ReadOnly => "false",
            Self::Mutating => "true",
            Self::External => "external",
        }
    }
}

/// Source class for a capability row in the shared registry.
#[derive(Debug, Clone, Copy, Serialize, Deserialize, PartialEq, Eq, PartialOrd, Ord, Hash)]
#[serde(rename_all = "snake_case")]
pub enum CapabilitySource {
    GlossaryCommand,
    EngineOperation,
    McpTool,
}

impl CapabilitySource {
    pub fn as_str(self) -> &'static str {
        match self {
            Self::GlossaryCommand => "glossary-command",
            Self::EngineOperation => "engine-operation",
            Self::McpTool => "mcp-tool",
        }
    }
}

/// Per-adapter surfacing state for one shared capability descriptor.
#[derive(Debug, Clone, Copy, Serialize, Deserialize, PartialEq, Eq)]
#[serde(rename_all = "snake_case")]
pub enum AdapterSurfacing {
    Prominent,
    ShellPassthrough,
    Gap,
    NotApplicable,
}

impl AdapterSurfacing {
    pub fn as_str(self) -> &'static str {
        match self {
            Self::Prominent => "prominent",
            Self::ShellPassthrough => "shell-only",
            Self::Gap => "gap",
            Self::NotApplicable => "n/a",
        }
    }

    pub fn is_reachable(self) -> bool {
        matches!(self, Self::Prominent | Self::ShellPassthrough)
    }
}

/// Machine-readable capability row shared by GUI, CLI, MCP, JS, Lua, and ClawBio adapters.
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
pub struct CapabilityDescriptor {
    pub name: String,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub title: Option<String>,
    pub description: String,
    pub input_schema: Value,
    pub output_schema: Value,
    pub mutating: CapabilityMutation,
    pub gui: AdapterSurfacing,
    pub cli: AdapterSurfacing,
    pub mcp: AdapterSurfacing,
    pub js: AdapterSurfacing,
    pub lua: AdapterSurfacing,
    pub clawbio: AdapterSurfacing,
    #[serde(default, skip_serializing_if = "BTreeMap::is_empty")]
    pub surfacing_justifications: BTreeMap<String, String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub inline_operand_ok: Option<bool>,
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub engine_operations: Vec<String>,
    pub source: CapabilitySource,
}

/// Machine-readable descriptor for one GENtle-local slash alias.
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
pub struct ShellAliasDescriptor {
    /// Stable alias key, such as `/fetch genbank`.
    pub alias: String,
    /// Human/agent-facing syntax for this alias.
    pub surface_form: String,
    /// Canonical shared-shell command or operation template that the alias
    /// resolves to.
    pub canonical_command: String,
    /// Short description used by docs, prompt contracts, and rejection help.
    pub description: String,
    /// Mutation/network safety class consumed by agent confirmation policy.
    pub mutating: CapabilityMutation,
    pub gui: AdapterSurfacing,
    pub cli: AdapterSurfacing,
    pub mcp: AdapterSurfacing,
    pub js: AdapterSurfacing,
    pub lua: AdapterSurfacing,
    pub clawbio: AdapterSurfacing,
    /// Engine operation variants reached by this alias, if any.
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub engine_operations: Vec<String>,
    /// Registry capability names this alias delegates to.
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub capability_names: Vec<String>,
    /// Nearby valid alternatives returned when the alias family is misused.
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub supported_alternatives: Vec<String>,
}

impl ShellAliasDescriptor {
    pub fn surfacing_for_adapter(&self, adapter: CapabilityAdapter) -> AdapterSurfacing {
        match adapter {
            CapabilityAdapter::Gui => self.gui,
            CapabilityAdapter::Cli => self.cli,
            CapabilityAdapter::Mcp => self.mcp,
            CapabilityAdapter::Js => self.js,
            CapabilityAdapter::Lua => self.lua,
            CapabilityAdapter::Clawbio => self.clawbio,
        }
    }
}

impl CapabilityDescriptor {
    pub fn surfacing_for_adapter(&self, adapter: CapabilityAdapter) -> AdapterSurfacing {
        match adapter {
            CapabilityAdapter::Gui => self.gui,
            CapabilityAdapter::Cli => self.cli,
            CapabilityAdapter::Mcp => self.mcp,
            CapabilityAdapter::Js => self.js,
            CapabilityAdapter::Lua => self.lua,
            CapabilityAdapter::Clawbio => self.clawbio,
        }
    }

    pub fn is_reachable_from_adapter(&self, adapter: CapabilityAdapter) -> bool {
        self.surfacing_for_adapter(adapter).is_reachable()
    }
}

#[derive(Debug, Deserialize)]
struct GlossaryRegistry {
    commands: Vec<GlossaryRegistryCommand>,
}

#[derive(Debug, Deserialize)]
struct GlossaryRegistryCommand {
    path: String,
    summary: String,
    #[serde(default)]
    interfaces: Vec<String>,
    #[serde(default)]
    engine_operations: Vec<String>,
    usage: String,
}

#[derive(Debug, Deserialize)]
struct ParityMatrixOverrides {
    overrides: Vec<ParityMatrixOverride>,
}

#[derive(Debug, Deserialize)]
struct ParityMatrixOverride {
    source: CapabilitySource,
    name: String,
    adapter: CapabilityAdapter,
    surfacing: AdapterSurfacing,
    reason: String,
}

const PUBLIC_ENGINE_OPERATION_NAMES: &[&str] = &[
    "LoadFile",
    "SaveFile",
    "RenderSequenceSvg",
    "RenderDotplotSvg",
    "RenderTfbsScoreTracksSvg",
    "RenderTfbsScoreTrackCorrelationSvg",
    "RenderFeatureExpertSvg",
    "RenderIsoformArchitectureSvg",
    "RenderRnaStructureSvg",
    "RenderLineageSvg",
    "RenderPoolGelSvg",
    "RenderProteinGelSvg",
    "RenderProteinGelReportsSvg",
    "RenderProteaseDigestGelSvg",
    "RenderProtein2dGelSvg",
    "RenderProtocolCartoonSvg",
    "RenderProtocolCartoonTemplateSvg",
    "ValidateProtocolCartoonTemplate",
    "RenderProtocolCartoonTemplateWithBindingsSvg",
    "ExportProtocolCartoonTemplateJson",
    "ApplyGibsonAssemblyPlan",
    "CreateArrangementSerial",
    "SetArrangementLadders",
    "SetContainerDeclaredContentsExclusive",
    "CreateRackFromArrangement",
    "PlaceArrangementOnRack",
    "MoveRackPlacement",
    "MoveRackSamples",
    "MoveRackArrangementBlocks",
    "SetRackProfile",
    "ApplyRackTemplate",
    "SetRackFillDirection",
    "SetRackProfileCustom",
    "SetRackBlockedCoordinates",
    "ExportRackLabelsSvg",
    "ExportRackFabricationSvg",
    "ExportRackIsometricSvg",
    "ExportRackHeroSvg",
    "ExportRackOpenScad",
    "ExportRackCarrierLabelsSvg",
    "ExportRackSimulationJson",
    "ExportDnaLadders",
    "ExportRnaLadders",
    "ExportPool",
    "ExportProcessRunBundle",
    "ExportLabAssistantInstructions",
    "PrepareGenome",
    "ExtractGenomeRegion",
    "ExtractGenomeGene",
    "ExtractGenomePromoterSlice",
    "ExtendGenomeAnchor",
    "ImportGenomeBedTrack",
    "ImportGenomeBigWigTrack",
    "ImportGenomeVcfTrack",
    "ProjectMicroarrayTrack",
    "ProjectGenomeInterval",
    "ListCutRunDatasets",
    "ShowCutRunDatasetStatus",
    "PrepareCutRunDataset",
    "ProjectCutRunDataset",
    "InterpretCutRunReads",
    "ListCutRunReadReports",
    "ShowCutRunReadReport",
    "ExportCutRunReadCoverage",
    "InspectCutRunRegulatorySupport",
    "ImportIsoformPanel",
    "ImportUniprotSwissProt",
    "FetchUniprotSwissProt",
    "FetchEnsemblGene",
    "FetchEnsemblRegion",
    "FetchEnsemblProtein",
    "FetchGenBankAccession",
    "FetchDbSnpRegion",
    "FetchUniprotLinkedGenBank",
    "ImportUniprotEntrySequence",
    "ImportEnsemblGeneSequence",
    "ImportEnsemblProteinSequence",
    "ProjectUniprotToGenome",
    "QueryProteinResidueGenomicCoordinates",
    "AuditUniprotProjectionConsistency",
    "AuditUniprotProjectionParity",
    "ImportBlastHitsTrack",
    "DigestContainer",
    "MergeContainersById",
    "LigationContainer",
    "FilterContainerByMolecularWeight",
    "Digest",
    "FindRestrictionSites",
    "QueryRepeatAnnotations",
    "QueryRepeatOverlaps",
    "MaterializeRepeatFeatures",
    "BuildRepeatEnvironmentCohort",
    "MergeContainers",
    "Ligation",
    "Pcr",
    "PcrAdvanced",
    "PcrMutagenesis",
    "DesignPrimerPairs",
    "DesignInsertionPrimerPairs",
    "ExportPrimerDesignReport",
    "AssessPrimerPairSpecificity",
    "PrepareRestrictionCloningPcrHandoff",
    "PcrOverlapExtensionMutagenesis",
    "DesignQpcrAssays",
    "TestCdnaPcr",
    "TestCdnaQpcr",
    "TestCdnaQpcrFasta",
    "DeriveTranscriptSequences",
    "PlanExonSkippedIsoform",
    "MaterializeExonSkippedIsoform",
    "DeriveProteinSequences",
    "ReverseTranslateProteinSequence",
    "ProteaseDigestProteinSequence",
    "BuildProteinToDnaHandoffReasoning",
    "ComputeDotplot",
    "ComputeDotplotOverlay",
    "ComputeFlexibilityTrack",
    "DeriveSplicingReferences",
    "AlignSequences",
    "ConfirmConstructReads",
    "SuggestSequencingPrimers",
    "ListSequencingConfirmationReports",
    "ShowSequencingConfirmationReport",
    "ExportSequencingConfirmationReport",
    "ExportSequencingConfirmationSupportTsv",
    "ReadAcquireStatus",
    "ReadAcquirePrepare",
    "ReadAcquireInspect",
    "ReadAcquireCancel",
    "InterpretRnaReads",
    "AlignRnaReadReport",
    "PreflightRnaReadIsoforms",
    "ListRnaReadReports",
    "ShowRnaReadReport",
    "SummarizeRnaReadGeneSupport",
    "InspectRnaReadGeneSupport",
    "ExportRnaReadIsoformTriageTsv",
    "RunRnaReadBatchMap",
    "SummarizeTfbsRegion",
    "SummarizeTfbsScoreTracks",
    "SummarizeTfbsTrackSimilarity",
    "SummarizeAlternativePromoterComparison",
    "SummarizePromoterEvidenceMatrix",
    "SummarizeIsoformPromoterComparison",
    "SummarizePromoterExpressionEvidence",
    "ExportPromoterArtifactManifest",
    "SummarizeMultiGenePromoterTfbs",
    "RenderMultiGenePromoterTfbsSvg",
    "ListReporterCatalog",
    "RecommendReporters",
    "ExportReporterCorpus",
    "PlanReporterConstructHandoff",
    "ScanTfbsHits",
    "InspectJasparEntry",
    "SummarizeJasparEntries",
    "ResolveTfQueries",
    "BenchmarkJasparRegistry",
    "ListJasparCatalog",
    "SyncJasparRemoteMetadata",
    "AnnotatePromoterWindows",
    "SummarizeVariantPromoterContext",
    "SuggestPromoterReporterFragments",
    "MaterializeVariantAllele",
    "ExportRnaReadReport",
    "ExportRnaReadHitsFasta",
    "ExportRnaReadSampleSheet",
    "ExportRnaReadTargetQuality",
    "ExportRnaReadExonPathsTsv",
    "ExportRnaReadExonAbundanceTsv",
    "ExportRnaReadScoreDensitySvg",
    "ExportRnaReadAlignmentsTsv",
    "ExportRnaReadAlignmentDotplotSvg",
    "MaterializeRnaReadHitSequences",
    "ExtractRegion",
    "ExtractAnchoredRegion",
    "SelectCandidate",
    "FilterByMolecularWeight",
    "FilterByDesignConstraints",
    "GenerateCandidateSet",
    "GenerateCandidateSetBetweenAnchors",
    "DeleteCandidateSet",
    "UpsertGuideSet",
    "DeleteGuideSet",
    "FilterGuidesPractical",
    "GenerateGuideOligos",
    "ExportGuideOligos",
    "ExportGuideProtocolText",
    "ExportFeaturesBed",
    "InspectSequenceContextView",
    "ExportSequenceContextBundle",
    "ScoreCandidateSetExpression",
    "ScoreCandidateSetDistance",
    "FilterCandidateSet",
    "CandidateSetOp",
    "ScoreCandidateSetWeightedObjective",
    "TopKCandidateSet",
    "ParetoFrontierCandidateSet",
    "UpsertWorkflowMacroTemplate",
    "DeleteWorkflowMacroTemplate",
    "UpsertCandidateMacroTemplate",
    "DeleteCandidateMacroTemplate",
    "Reverse",
    "Complement",
    "ReverseComplement",
    "Branch",
    "SetDisplayVisibility",
    "SetLinearViewport",
    "SetTopology",
    "RecomputeFeatures",
    "SetParameter",
    "AnnotateTfbs",
];

const MCP_TOOL_NAMES: &[(&str, &str, &str, CapabilityMutation)] = &[
    (
        "capabilities",
        "Capabilities",
        "Return shared GENtle engine capabilities.",
        CapabilityMutation::ReadOnly,
    ),
    (
        "state_summary",
        "State Summary",
        "Return a deterministic summary of the current project state.",
        CapabilityMutation::ReadOnly,
    ),
    (
        "exon_skip_plan",
        "Exon Skip Plan",
        "Build and persist an inspectable exon-skip selection plan through the shared `transcripts exon-skip-plan` shell contract.",
        CapabilityMutation::Mutating,
    ),
    (
        "exon_skip_materialize",
        "Exon Skip Materialize",
        "Materialize one stored exon-skip plan through the shared `transcripts exon-skip-materialize` shell contract and optionally return caller-requested payloads.",
        CapabilityMutation::Mutating,
    ),
    (
        "restriction_site_detail",
        "Restriction Site Detail",
        "Return one restriction-site expert record through the shared `inspect-feature-expert ... restriction` shell contract.",
        CapabilityMutation::ReadOnly,
    ),
    (
        "agent_systems",
        "Agent Systems",
        "Return configured GENtle agent systems via the shared `agents list` contract.",
        CapabilityMutation::ReadOnly,
    ),
    (
        "agent_preflight",
        "Agent Preflight",
        "Return transport/runtime preflight details for one GENtle agent system.",
        CapabilityMutation::External,
    ),
    (
        "agent_models",
        "Agent Models",
        "Discover model ids for one GENtle OpenAI-compatible/native agent system.",
        CapabilityMutation::External,
    ),
    (
        "agent_plan",
        "Agent Plan",
        "Compile free prose into a typed GENtle execution plan via the shared `agents plan` contract.",
        CapabilityMutation::External,
    ),
    (
        "agent_execute_plan",
        "Agent Execute Plan",
        "Execute one stored planner candidate through the shared `agents execute-plan` contract.",
        CapabilityMutation::Mutating,
    ),
    (
        "op",
        "Apply Operation",
        "Apply one operation via shared engine contract and persist state (requires confirm=true).",
        CapabilityMutation::Mutating,
    ),
    (
        "workflow",
        "Apply Workflow",
        "Apply a workflow via shared engine contract and persist state (requires confirm=true).",
        CapabilityMutation::Mutating,
    ),
    (
        "help",
        "Help",
        "Return shell command reference content from docs/glossary.json.",
        CapabilityMutation::ReadOnly,
    ),
    (
        "reference_catalog_entries",
        "Reference Catalog Entries",
        "Return structured reference catalog entries via the shared `genomes list` contract.",
        CapabilityMutation::ReadOnly,
    ),
    (
        "helper_catalog_entries",
        "Helper Catalog Entries",
        "Return structured helper catalog entries, including normalized helper interpretations when available.",
        CapabilityMutation::ReadOnly,
    ),
    (
        "helper_semantics_vocabulary",
        "Helper Semantics Vocabulary",
        "Return resolved helper semantics vocabulary terms, aliases, descriptions, sources, and routine hints.",
        CapabilityMutation::ReadOnly,
    ),
    (
        "host_profile_catalog_entries",
        "Host Profile Catalog Entries",
        "Return structured host-profile catalog entries used by construct reasoning.",
        CapabilityMutation::ReadOnly,
    ),
    (
        "ensembl_installable_genomes",
        "Ensembl Installable Genomes",
        "Return Ensembl species directories that currently appear installable because both sequence and annotation listings are present.",
        CapabilityMutation::External,
    ),
    (
        "construct_reasoning_graphs",
        "Construct Reasoning Graphs",
        "Return stored construct-reasoning graph summaries through the shared `construct-reasoning list-graphs` shell contract.",
        CapabilityMutation::ReadOnly,
    ),
    (
        "construct_reasoning_graph",
        "Construct Reasoning Graph",
        "Return one stored construct-reasoning graph plus compact summary details through the shared `construct-reasoning show-graph` shell contract.",
        CapabilityMutation::ReadOnly,
    ),
    (
        "construct_reasoning_inspection_actions",
        "Construct Reasoning Inspection Actions",
        "Return recommended inspection actions from one stored construct-reasoning graph through the shared `construct-reasoning list-inspection-actions` shell contract.",
        CapabilityMutation::ReadOnly,
    ),
    (
        "construct_reasoning_run_inspection_action",
        "Run Construct Reasoning Inspection Action",
        "Compute the dotplot recommended by one construct-reasoning inspection action through the shared `construct-reasoning run-inspection-action` shell contract.",
        CapabilityMutation::Mutating,
    ),
    (
        "construct_reasoning_set_annotation_status",
        "Set Construct Reasoning Annotation Status",
        "Update one stored construct-reasoning annotation candidate to draft, accepted, rejected, or locked through the shared `construct-reasoning set-annotation-status` shell contract.",
        CapabilityMutation::Mutating,
    ),
    (
        "construct_reasoning_write_annotation",
        "Write Back Construct Reasoning Annotation",
        "Materialize one accepted or locked generated construct-reasoning annotation candidate as an ordinary sequence feature through the shared `construct-reasoning write-annotation` shell contract.",
        CapabilityMutation::Mutating,
    ),
    (
        "helper_interpretation",
        "Helper Interpretation",
        "Return the normalized helper-construct interpretation for one helper id or alias.",
        CapabilityMutation::ReadOnly,
    ),
    (
        "ui_intents",
        "UI Intents Catalog",
        "Return deterministic UI-intent target/command catalog (shared `ui intents` contract).",
        CapabilityMutation::ReadOnly,
    ),
    (
        "ui_intent",
        "UI Intent",
        "Resolve/record one UI intent through shared `ui open|focus` parser/executor path.",
        CapabilityMutation::Mutating,
    ),
    (
        "ui_prepared_genomes",
        "UI Prepared Genomes Query",
        "Return deterministic prepared-genome rows via shared `ui prepared-genomes` contract.",
        CapabilityMutation::ReadOnly,
    ),
    (
        "ui_latest_prepared",
        "UI Latest Prepared",
        "Resolve latest prepared genome for a species via shared `ui latest-prepared` contract.",
        CapabilityMutation::ReadOnly,
    ),
    (
        "blast_async_start",
        "BLAST Async Start",
        "Start one async BLAST job through the shared shell contract (`genomes/helpers blast-start`).",
        CapabilityMutation::External,
    ),
    (
        "blast_async_status",
        "BLAST Async Status",
        "Inspect one async BLAST job status and optional report payload.",
        CapabilityMutation::ReadOnly,
    ),
    (
        "blast_async_cancel",
        "BLAST Async Cancel",
        "Request cancellation for one async BLAST job.",
        CapabilityMutation::External,
    ),
    (
        "blast_async_list",
        "BLAST Async List",
        "List known async BLAST jobs for genome/helper scope.",
        CapabilityMutation::ReadOnly,
    ),
];

const MCP_PROMINENT_GLOSSARY_COMMAND_PATHS: &[&str] = &[
    "capabilities",
    "state-summary",
    "agents list",
    "agents preflight",
    "agents discover-models",
    "agents plan",
    "agents execute-plan",
    "op",
    "workflow",
    "help",
    "genomes list",
    "helpers list",
    "helpers vocabulary list",
    "hosts list",
    "genomes ensembl-available",
    "construct-reasoning list-graphs",
    "construct-reasoning show-graph",
    "construct-reasoning set-annotation-status",
    "construct-reasoning write-annotation",
    "ui intents",
    "ui open",
    "ui focus",
    "ui prepared-genomes",
    "ui latest-prepared",
    "genomes blast-start",
    "helpers blast-start",
    "genomes blast-status",
    "helpers blast-status",
    "genomes blast-cancel",
    "helpers blast-cancel",
    "genomes blast-list",
    "helpers blast-list",
    "transcripts exon-skip-plan",
    "transcripts exon-skip-materialize",
];

fn is_mcp_prominent_glossary_command(path: &str) -> bool {
    MCP_PROMINENT_GLOSSARY_COMMAND_PATHS.contains(&path)
}

static CAPABILITY_REGISTRY: LazyLock<Vec<CapabilityDescriptor>> =
    LazyLock::new(build_capability_registry);
static SHELL_ALIAS_REGISTRY: LazyLock<Vec<ShellAliasDescriptor>> =
    LazyLock::new(build_shell_alias_registry);
static PARITY_MATRIX_OVERRIDES: LazyLock<Vec<ParityMatrixOverride>> =
    LazyLock::new(load_parity_matrix_overrides);

const CAPABILITY_PARITY_ADAPTERS: [CapabilityAdapter; 6] = [
    CapabilityAdapter::Gui,
    CapabilityAdapter::Cli,
    CapabilityAdapter::Mcp,
    CapabilityAdapter::Js,
    CapabilityAdapter::Lua,
    CapabilityAdapter::Clawbio,
];

/// Protocol-owned names for public engine operation variants.
pub fn public_engine_operation_names() -> &'static [&'static str] {
    PUBLIC_ENGINE_OPERATION_NAMES
}

/// Shared capability registry projected by adapters.
pub fn capability_registry() -> &'static [CapabilityDescriptor] {
    &CAPABILITY_REGISTRY
}

/// GENtle-local slash aliases accepted by the shared shell parser.
pub fn shell_alias_registry() -> &'static [ShellAliasDescriptor] {
    &SHELL_ALIAS_REGISTRY
}

/// Return one GENtle-local slash alias descriptor by stable alias key.
pub fn shell_alias_descriptor(alias: &str) -> Option<&'static ShellAliasDescriptor> {
    shell_alias_registry()
        .iter()
        .find(|descriptor| descriptor.alias == alias)
}

/// Suggested alternatives for invalid or unsupported local slash commands.
pub fn shell_alias_supported_alternatives() -> Vec<String> {
    shell_alias_registry()
        .iter()
        .map(|descriptor| descriptor.surface_form.clone())
        .collect()
}

/// Adapter order used by the generated GUI/CLI/MCP parity matrix.
pub fn capability_parity_adapters() -> &'static [CapabilityAdapter] {
    &CAPABILITY_PARITY_ADAPTERS
}

/// Return one capability descriptor by stable name and source.
pub fn capability_descriptor(
    source: CapabilitySource,
    name: &str,
) -> Option<&'static CapabilityDescriptor> {
    capability_registry()
        .iter()
        .find(|descriptor| descriptor.source == source && descriptor.name == name)
}

/// Return descriptors exposed by one adapter.
pub fn capability_registry_for_adapter(adapter: CapabilityAdapter) -> Vec<CapabilityDescriptor> {
    capability_registry()
        .iter()
        .filter(|descriptor| descriptor.is_reachable_from_adapter(adapter))
        .cloned()
        .collect()
}

fn load_parity_matrix_overrides() -> Vec<ParityMatrixOverride> {
    let overrides: ParityMatrixOverrides =
        serde_json::from_str(include_str!("../../../docs/parity_matrix_overrides.json"))
            .expect("docs/parity_matrix_overrides.json must parse");
    for entry in &overrides.overrides {
        assert_eq!(
            entry.surfacing,
            AdapterSurfacing::NotApplicable,
            "parity matrix override for {:?} `{}` on {:?} must declare surfacing=not_applicable",
            entry.source,
            entry.name,
            entry.adapter
        );
        assert!(
            !entry.reason.trim().is_empty(),
            "parity matrix override for {:?} `{}` on {:?} must have a reason",
            entry.source,
            entry.name,
            entry.adapter
        );
    }
    overrides.overrides
}

fn not_applicable_override_reason(
    source: CapabilitySource,
    name: &str,
    adapter: CapabilityAdapter,
) -> Option<String> {
    PARITY_MATRIX_OVERRIDES
        .iter()
        .find(|entry| entry.source == source && entry.name == name && entry.adapter == adapter)
        .or_else(|| {
            PARITY_MATRIX_OVERRIDES.iter().find(|entry| {
                entry.source == source && entry.name == "*" && entry.adapter == adapter
            })
        })
        .map(|entry| entry.reason.clone())
}

/// Render the GUI/CLI/MCP parity matrix from the protocol capability registry.
pub fn render_gui_cli_mcp_parity_matrix_markdown() -> String {
    render_capability_parity_matrix_markdown(capability_registry())
}

fn render_capability_parity_matrix_markdown(registry: &[CapabilityDescriptor]) -> String {
    let mut out = String::new();
    out.push_str("# GUI / CLI / MCP Parity Matrix\n\n");
    out.push_str("Generated automatically from the protocol capability registry.\n\n");
    out.push_str("## Method\n\n");
    out.push_str(
        "This file is a projection of `gentle_protocol::capability_registry()`. \
The cell vocabulary is the architecture-doc parity policy:\n\n",
    );
    out.push_str("- `prominent`: `AdapterSurfacing::Prominent` for this adapter\n");
    out.push_str("- `shell-only`: `AdapterSurfacing::ShellPassthrough`\n");
    out.push_str(
        "- `n/a`: `AdapterSurfacing::NotApplicable`, with a one-line Notes justification\n",
    );
    out.push_str("- `gap`: operation should be reachable but currently is not\n\n");
    out.push_str(
        "Only `gap` signals implementation work. Human-readable Notes are populated \
from each descriptor's `surfacing_justifications` map.\n\n",
    );
    out.push_str("## Findings\n\n");
    out.push_str("| Adapter | prominent | shell-only | gap |\n");
    out.push_str("|---|---:|---:|---:|\n");
    for adapter in capability_parity_adapters() {
        let mut prominent = 0usize;
        let mut shell_only = 0usize;
        let mut gap = 0usize;
        for descriptor in registry {
            match descriptor.surfacing_for_adapter(*adapter) {
                AdapterSurfacing::Prominent => prominent += 1,
                AdapterSurfacing::ShellPassthrough => shell_only += 1,
                AdapterSurfacing::Gap => gap += 1,
                AdapterSurfacing::NotApplicable => {}
            }
        }
        out.push_str(&format!(
            "| {} | {} | {} | {} |\n",
            adapter.label(),
            prominent,
            shell_only,
            gap
        ));
    }
    out.push('\n');
    push_parity_section(
        &mut out,
        "Glossary Commands",
        registry
            .iter()
            .filter(|descriptor| descriptor.source == CapabilitySource::GlossaryCommand),
    );
    push_parity_section(
        &mut out,
        "Engine Operations",
        registry
            .iter()
            .filter(|descriptor| descriptor.source == CapabilitySource::EngineOperation),
    );
    push_parity_section(
        &mut out,
        "MCP Tools",
        registry
            .iter()
            .filter(|descriptor| descriptor.source == CapabilitySource::McpTool),
    );
    push_open_gaps_section(&mut out, registry);
    while out.ends_with("\n\n") {
        out.pop();
    }
    out
}

fn push_open_gaps_section(out: &mut String, registry: &[CapabilityDescriptor]) {
    out.push_str("## Open Gaps\n\n");
    out.push_str("| Adapter | Capability | Source | Triage | Note |\n");
    out.push_str("|---|---|---|---|---|\n");
    let mut gaps = Vec::new();
    for descriptor in registry {
        for adapter in capability_parity_adapters() {
            if descriptor.surfacing_for_adapter(*adapter) != AdapterSurfacing::Gap {
                continue;
            }
            let (triage, note) = gap_triage_note(descriptor);
            gaps.push((
                adapter.label().to_string(),
                descriptor.name.clone(),
                descriptor.source.as_str().to_string(),
                triage.to_string(),
                note,
            ));
        }
    }
    gaps.sort_by(|left, right| {
        left.0
            .cmp(&right.0)
            .then_with(|| left.1.cmp(&right.1))
            .then_with(|| left.2.cmp(&right.2))
    });
    if gaps.is_empty() {
        out.push_str("| _None_ |  |  |  |\n");
    } else {
        for (adapter, name, source, triage, note) in gaps {
            out.push_str("| ");
            out.push_str(&markdown_cell(&adapter));
            out.push_str(" | ");
            out.push_str(&markdown_cell(&name));
            out.push_str(" | ");
            out.push_str(&markdown_cell(&source));
            out.push_str(" | ");
            out.push_str(&markdown_cell(&triage));
            out.push_str(" | ");
            out.push_str(&markdown_cell(&note));
            out.push_str(" |\n");
        }
    }
    out.push('\n');
}

fn gap_triage_note(descriptor: &CapabilityDescriptor) -> (&'static str, String) {
    if descriptor.engine_operations.is_empty() {
        (
            "needs operation work",
            "No portable engine operation or curated local-only override is declared.".to_string(),
        )
    } else {
        (
            "ready to wire",
            format!(
                "Engine operation route is declared: {}.",
                descriptor.engine_operations.join(", ")
            ),
        )
    }
}

fn push_parity_section<'a>(
    out: &mut String,
    title: &str,
    descriptors: impl Iterator<Item = &'a CapabilityDescriptor>,
) {
    out.push_str(&format!("## {title}\n\n"));
    out.push_str("| Capability | Source | GUI | gentle_cli | MCP | JS | Lua | ClawBio | Notes |\n");
    out.push_str("|---|---|---|---|---|---|---|---|---|\n");
    for descriptor in descriptors {
        out.push_str("| ");
        out.push_str(&markdown_cell(&descriptor.name));
        out.push_str(" | ");
        out.push_str(descriptor.source.as_str());
        for adapter in capability_parity_adapters() {
            let value = descriptor.surfacing_for_adapter(*adapter).as_str();
            out.push_str(" | ");
            out.push_str(value);
        }
        out.push_str(" | ");
        out.push_str(&descriptor_notes(descriptor));
        out.push_str(" |\n");
    }
    out.push('\n');
}

fn descriptor_notes(descriptor: &CapabilityDescriptor) -> String {
    let notes = capability_parity_adapters()
        .iter()
        .filter_map(|adapter| {
            descriptor
                .surfacing_justifications
                .get(adapter.as_str())
                .map(|note| format!("{}: {}", adapter.label(), note))
        })
        .collect::<Vec<_>>();
    markdown_cell(&notes.join("<br>"))
}

fn markdown_cell(raw: &str) -> String {
    raw.replace('\n', " ").replace('|', "\\|")
}

fn build_capability_registry() -> Vec<CapabilityDescriptor> {
    let glossary: GlossaryRegistry =
        serde_json::from_str(include_str!("../../../docs/glossary.json"))
            .expect("docs/glossary.json must parse for capability registry");
    let interfaces_by_path = glossary_interfaces_by_path(&glossary.commands);
    let mut descriptors = vec![];
    for command in &glossary.commands {
        descriptors.push(glossary_command_descriptor(command, &interfaces_by_path));
    }
    for op_name in PUBLIC_ENGINE_OPERATION_NAMES {
        descriptors.push(engine_operation_descriptor(op_name, &glossary.commands));
    }
    for (name, title, description, mutating) in MCP_TOOL_NAMES {
        descriptors.push(mcp_tool_descriptor(name, title, description, *mutating));
    }
    descriptors
}

fn build_shell_alias_registry() -> Vec<ShellAliasDescriptor> {
    let prominent = AdapterSurfacing::Prominent;
    let shell = AdapterSurfacing::ShellPassthrough;
    let na = AdapterSurfacing::NotApplicable;
    let read_only_alternatives = &["/help", "/list"];
    let local_alternatives = &[
        "/help",
        "/list",
        "/open",
        "/import",
        "/open file PATH [--id ID]",
        "/import file PATH [--id ID]",
        "/paste sequence --sequence-text DNA [--id ID]",
        "/features restriction-scan SEQ_ID [--enzyme NAME]",
        "/fetch genbank ACCESSION [--id ID]",
        "/fetch ensembl QUERY [--species NAME] [--id ID]",
        "/fetch uniprot QUERY [--id ID]",
    ];
    vec![
        shell_alias_descriptor_row(
            "/help",
            "/help [TOPIC]",
            "help [TOPIC]",
            "Show GENtle shared-shell help.",
            CapabilityMutation::ReadOnly,
            &[],
            &["help"],
            read_only_alternatives,
            prominent,
            shell,
            shell,
            na,
            na,
            na,
        ),
        shell_alias_descriptor_row(
            "/list",
            "/list",
            "state-summary",
            "List the current GENtle workspace state.",
            CapabilityMutation::ReadOnly,
            &[],
            &["state-summary"],
            read_only_alternatives,
            prominent,
            shell,
            shell,
            na,
            na,
            na,
        ),
        shell_alias_descriptor_row(
            "/open",
            "/open",
            "ui open open-sequence",
            "Open the GUI sequence-file picker.",
            CapabilityMutation::Mutating,
            &[],
            &[],
            local_alternatives,
            prominent,
            shell,
            shell,
            na,
            na,
            na,
        ),
        shell_alias_descriptor_row(
            "/import",
            "/import",
            "ui open open-sequence",
            "Open the GUI sequence-file picker.",
            CapabilityMutation::Mutating,
            &[],
            &[],
            local_alternatives,
            prominent,
            shell,
            shell,
            na,
            na,
            na,
        ),
        shell_alias_descriptor_row(
            "/open file",
            "/open file PATH [--id ID]",
            "op '{\"LoadFile\":{\"path\":\"PATH\",\"as_id\":\"ID\"}}'",
            "Load a user-specified local sequence file into the project.",
            CapabilityMutation::Mutating,
            &["LoadFile"],
            &["LoadFile"],
            local_alternatives,
            prominent,
            shell,
            shell,
            na,
            na,
            na,
        ),
        shell_alias_descriptor_row(
            "/import file",
            "/import file PATH [--id ID]",
            "op '{\"LoadFile\":{\"path\":\"PATH\",\"as_id\":\"ID\"}}'",
            "Load a user-specified local sequence file into the project.",
            CapabilityMutation::Mutating,
            &["LoadFile"],
            &["LoadFile"],
            local_alternatives,
            prominent,
            shell,
            shell,
            na,
            na,
            na,
        ),
        shell_alias_descriptor_row(
            "/paste sequence",
            "/paste sequence --sequence-text DNA [--id ID]",
            "sequence create --sequence-text DNA [--id ID]",
            "Create a project sequence from explicit pasted IUPAC sequence text.",
            CapabilityMutation::Mutating,
            &["CreateSequenceFromText"],
            &["sequence create"],
            local_alternatives,
            prominent,
            shell,
            shell,
            na,
            na,
            na,
        ),
        shell_alias_descriptor_row(
            "/features restriction-scan",
            "/features restriction-scan SEQ_ID [--enzyme NAME]",
            "features restriction-scan SEQ_ID [--enzyme NAME]",
            "Scan restriction-enzyme sites on one loaded sequence.",
            CapabilityMutation::ReadOnly,
            &["FindRestrictionSites"],
            &["features restriction-scan"],
            local_alternatives,
            prominent,
            shell,
            shell,
            na,
            na,
            na,
        ),
        shell_alias_descriptor_row(
            "/fetch genbank",
            "/fetch genbank ACCESSION [--id ID]",
            "genbank fetch ACCESSION [--as-id ID]",
            "Fetch an NCBI GenBank nucleotide accession; --id names the local GENtle sequence id.",
            CapabilityMutation::External,
            &["FetchGenBankAccession"],
            &["genbank fetch"],
            local_alternatives,
            prominent,
            shell,
            shell,
            na,
            na,
            na,
        ),
        shell_alias_descriptor_row(
            "/fetch ncbi",
            "/fetch ncbi ACCESSION [--id ID]",
            "genbank fetch ACCESSION [--as-id ID]",
            "Synonym for `/fetch genbank`; ACCESSION is an NCBI GenBank nucleotide accession and --id names the local GENtle sequence id.",
            CapabilityMutation::External,
            &["FetchGenBankAccession"],
            &["genbank fetch"],
            local_alternatives,
            prominent,
            shell,
            shell,
            na,
            na,
            na,
        ),
        shell_alias_descriptor_row(
            "/fetch uniprot",
            "/fetch uniprot QUERY [--id ID]",
            "uniprot fetch QUERY [--entry-id ID]",
            "Fetch a UniProtKB/Swiss-Prot entry by accession or entry name; --id names the local GENtle metadata entry.",
            CapabilityMutation::External,
            &["FetchUniprotSwissProt"],
            &["uniprot fetch"],
            local_alternatives,
            prominent,
            shell,
            shell,
            na,
            na,
            na,
        ),
        shell_alias_descriptor_row(
            "/fetch ensembl",
            "/fetch ensembl QUERY [--species NAME] [--id ID]",
            "ensembl-gene fetch QUERY [--species NAME] [--entry-id ID]",
            "Fetch an Ensembl gene lookup record by HGNC-approved symbol or stable Ensembl gene id; --species uses names such as homo_sapiens and --id names the local GENtle metadata entry.",
            CapabilityMutation::External,
            &["FetchEnsemblGene"],
            &["ensembl-gene fetch"],
            local_alternatives,
            prominent,
            shell,
            shell,
            na,
            na,
            na,
        ),
        shell_alias_descriptor_row(
            "/fetch ensembl-gene",
            "/fetch ensembl-gene QUERY [--species NAME] [--id ID]",
            "ensembl-gene fetch QUERY [--species NAME] [--entry-id ID]",
            "Fetch an Ensembl gene lookup record by HGNC-approved symbol or stable Ensembl gene id; --species uses names such as homo_sapiens and --id names the local GENtle metadata entry.",
            CapabilityMutation::External,
            &["FetchEnsemblGene"],
            &["ensembl-gene fetch"],
            local_alternatives,
            prominent,
            shell,
            shell,
            na,
            na,
            na,
        ),
        shell_alias_descriptor_row(
            "/fetch ensembl-protein",
            "/fetch ensembl-protein QUERY [--id ID]",
            "ensembl-protein fetch QUERY [--entry-id ID]",
            "Fetch an Ensembl protein/translation lookup record by stable id or accepted lookup query; --id names the local GENtle metadata entry.",
            CapabilityMutation::External,
            &["FetchEnsemblProtein"],
            &["ensembl-protein fetch"],
            local_alternatives,
            prominent,
            shell,
            shell,
            na,
            na,
            na,
        ),
        shell_alias_descriptor_row(
            "/fetch ensembl-region",
            "/fetch ensembl-region SPECIES CHR START END [--strand +|-] [--id ID]",
            "ensembl-region fetch SPECIES CHR START END [--strand +|-] [--output-id ID]",
            "Fetch an Ensembl REST genomic region; SPECIES is an Ensembl species name, coordinates are assembly positions, and --id names the local GENtle sequence id.",
            CapabilityMutation::External,
            &["FetchEnsemblRegion"],
            &["ensembl-region fetch"],
            local_alternatives,
            prominent,
            shell,
            shell,
            na,
            na,
            na,
        ),
        shell_alias_descriptor_row(
            "/fetch dbsnp",
            "/fetch dbsnp RS_ID GENOME_ID [--id ID]",
            "dbsnp fetch RS_ID GENOME_ID [--output-id ID]",
            "Fetch a dbSNP rsID-centered region from a prepared GENOME_ID reference; --id names the local GENtle sequence id.",
            CapabilityMutation::External,
            &["FetchDbSnpRegion"],
            &["dbsnp fetch"],
            local_alternatives,
            prominent,
            shell,
            shell,
            na,
            na,
            na,
        ),
    ]
}

#[allow(clippy::too_many_arguments)]
fn shell_alias_descriptor_row(
    alias: &str,
    surface_form: &str,
    canonical_command: &str,
    description: &str,
    mutating: CapabilityMutation,
    engine_operations: &[&str],
    capability_names: &[&str],
    supported_alternatives: &[&str],
    gui: AdapterSurfacing,
    cli: AdapterSurfacing,
    mcp: AdapterSurfacing,
    js: AdapterSurfacing,
    lua: AdapterSurfacing,
    clawbio: AdapterSurfacing,
) -> ShellAliasDescriptor {
    ShellAliasDescriptor {
        alias: alias.to_string(),
        surface_form: surface_form.to_string(),
        canonical_command: canonical_command.to_string(),
        description: description.to_string(),
        mutating,
        gui,
        cli,
        mcp,
        js,
        lua,
        clawbio,
        engine_operations: engine_operations
            .iter()
            .map(|value| (*value).to_string())
            .collect(),
        capability_names: capability_names
            .iter()
            .map(|value| (*value).to_string())
            .collect(),
        supported_alternatives: supported_alternatives
            .iter()
            .map(|value| (*value).to_string())
            .collect(),
    }
}

fn glossary_interfaces_by_path(
    commands: &[GlossaryRegistryCommand],
) -> BTreeMap<String, BTreeSet<String>> {
    let mut interfaces_by_path = BTreeMap::new();
    for command in commands {
        let interfaces = interfaces_by_path
            .entry(command.path.clone())
            .or_insert_with(BTreeSet::new);
        interfaces.extend(command.interfaces.iter().cloned());
    }
    interfaces_by_path
}

fn glossary_command_descriptor(
    command: &GlossaryRegistryCommand,
    interfaces_by_path: &BTreeMap<String, BTreeSet<String>>,
) -> CapabilityDescriptor {
    let surfacing = glossary_command_surfacing(command, interfaces_by_path);
    CapabilityDescriptor {
        name: command.path.clone(),
        title: Some(title_from_stable_name(&command.path)),
        description: command.summary.clone(),
        input_schema: json!({
            "type": "object",
            "description": "Command arguments are documented by the glossary usage string.",
            "properties": {
                "usage": { "const": command.usage }
            },
            "additionalProperties": true
        }),
        output_schema: generic_output_schema(),
        mutating: infer_command_mutation(&command.path, &command.engine_operations),
        gui: surfacing.gui,
        cli: surfacing.cli,
        mcp: surfacing.mcp,
        js: surfacing.js,
        lua: surfacing.lua,
        clawbio: surfacing.clawbio,
        surfacing_justifications: surfacing_justifications_for_glossary_command(
            &surfacing, command,
        ),
        inline_operand_ok: inline_operand_ok_for_operations(&command.engine_operations),
        engine_operations: command.engine_operations.clone(),
        source: CapabilitySource::GlossaryCommand,
    }
}

fn engine_operation_descriptor(
    op_name: &str,
    commands: &[GlossaryRegistryCommand],
) -> CapabilityDescriptor {
    let mut description = None;
    let surfacing = engine_operation_surfacing(op_name);
    for command in commands {
        if command
            .engine_operations
            .iter()
            .any(|operation| operation == op_name)
        {
            description.get_or_insert_with(|| command.summary.clone());
        }
    }
    CapabilityDescriptor {
        name: op_name.to_string(),
        title: Some(title_from_stable_name(op_name)),
        description: description.unwrap_or_else(|| format!("GENtle engine operation `{op_name}`.")),
        input_schema: json!({
            "type": "object",
            "description": format!("Serialized Operation::{op_name} payload.")
        }),
        output_schema: json!({
            "type": "object",
            "description": "Serialized OpResult payload."
        }),
        mutating: infer_engine_operation_mutation(op_name),
        gui: surfacing.gui,
        cli: surfacing.cli,
        mcp: surfacing.mcp,
        js: surfacing.js,
        lua: surfacing.lua,
        clawbio: surfacing.clawbio,
        surfacing_justifications: surfacing_justifications_for_engine_operation(
            &surfacing, op_name,
        ),
        inline_operand_ok: inline_operand_ok_for_operation(op_name),
        engine_operations: vec![op_name.to_string()],
        source: CapabilitySource::EngineOperation,
    }
}

fn mcp_tool_descriptor(
    name: &str,
    title: &str,
    description: &str,
    mutating: CapabilityMutation,
) -> CapabilityDescriptor {
    let surfacing = mcp_tool_surfacing(name);
    CapabilityDescriptor {
        name: name.to_string(),
        title: Some(title.to_string()),
        description: description.to_string(),
        input_schema: json!({
            "type": "object",
            "description": "See the MCP `tools/list` inputSchema for field-level details."
        }),
        output_schema: generic_output_schema(),
        mutating,
        gui: surfacing.gui,
        cli: surfacing.cli,
        mcp: surfacing.mcp,
        js: surfacing.js,
        lua: surfacing.lua,
        clawbio: surfacing.clawbio,
        surfacing_justifications: surfacing_justifications_for_mcp_tool(&surfacing, name),
        inline_operand_ok: None,
        engine_operations: vec![],
        source: CapabilitySource::McpTool,
    }
}

#[derive(Debug, Clone, Copy)]
struct CapabilitySurfacingSet {
    gui: AdapterSurfacing,
    cli: AdapterSurfacing,
    mcp: AdapterSurfacing,
    js: AdapterSurfacing,
    lua: AdapterSurfacing,
    clawbio: AdapterSurfacing,
}

impl CapabilitySurfacingSet {
    fn for_adapter(self, adapter: CapabilityAdapter) -> AdapterSurfacing {
        match adapter {
            CapabilityAdapter::Gui => self.gui,
            CapabilityAdapter::Cli => self.cli,
            CapabilityAdapter::Mcp => self.mcp,
            CapabilityAdapter::Js => self.js,
            CapabilityAdapter::Lua => self.lua,
            CapabilityAdapter::Clawbio => self.clawbio,
        }
    }

    fn set_for_adapter(&mut self, adapter: CapabilityAdapter, surfacing: AdapterSurfacing) {
        match adapter {
            CapabilityAdapter::Gui => self.gui = surfacing,
            CapabilityAdapter::Cli => self.cli = surfacing,
            CapabilityAdapter::Mcp => self.mcp = surfacing,
            CapabilityAdapter::Js => self.js = surfacing,
            CapabilityAdapter::Lua => self.lua = surfacing,
            CapabilityAdapter::Clawbio => self.clawbio = surfacing,
        }
    }
}

fn glossary_command_surfacing(
    command: &GlossaryRegistryCommand,
    interfaces_by_path: &BTreeMap<String, BTreeSet<String>>,
) -> CapabilitySurfacingSet {
    let interfaces = interfaces_by_path
        .get(&command.path)
        .expect("glossary command path exists in interface index");
    let mut surfacing = CapabilitySurfacingSet {
        gui: if interfaces.contains("gui-menu") {
            AdapterSurfacing::Prominent
        } else if interfaces.contains("gui-shell") {
            AdapterSurfacing::ShellPassthrough
        } else {
            AdapterSurfacing::Gap
        },
        cli: if interfaces.contains("cli-direct") {
            AdapterSurfacing::Prominent
        } else if interfaces.contains("cli-shell")
            || command_is_reachable_via_generic_operation_route(command, CapabilityAdapter::Cli)
        {
            AdapterSurfacing::ShellPassthrough
        } else {
            AdapterSurfacing::Gap
        },
        mcp: if interfaces.contains("mcp") || is_mcp_prominent_glossary_command(&command.path) {
            AdapterSurfacing::Prominent
        } else if command_is_reachable_via_generic_operation_route(command, CapabilityAdapter::Mcp)
        {
            AdapterSurfacing::ShellPassthrough
        } else {
            AdapterSurfacing::Gap
        },
        js: if interfaces.contains("js") {
            AdapterSurfacing::Prominent
        } else if command_is_reachable_via_generic_operation_route(command, CapabilityAdapter::Js) {
            AdapterSurfacing::ShellPassthrough
        } else {
            AdapterSurfacing::Gap
        },
        lua: if interfaces.contains("lua") {
            AdapterSurfacing::Prominent
        } else if command_is_reachable_via_generic_operation_route(command, CapabilityAdapter::Lua)
        {
            AdapterSurfacing::ShellPassthrough
        } else {
            AdapterSurfacing::Gap
        },
        clawbio: AdapterSurfacing::Gap,
    };
    apply_parity_matrix_overrides(
        CapabilitySource::GlossaryCommand,
        &command.path,
        &mut surfacing,
    );
    surfacing
}

fn command_is_reachable_via_generic_operation_route(
    command: &GlossaryRegistryCommand,
    adapter: CapabilityAdapter,
) -> bool {
    !command.engine_operations.is_empty()
        && matches!(
            adapter,
            CapabilityAdapter::Cli
                | CapabilityAdapter::Mcp
                | CapabilityAdapter::Js
                | CapabilityAdapter::Lua
        )
}

fn engine_operation_surfacing(op_name: &str) -> CapabilitySurfacingSet {
    let mut surfacing = CapabilitySurfacingSet {
        gui: AdapterSurfacing::ShellPassthrough,
        cli: AdapterSurfacing::ShellPassthrough,
        mcp: AdapterSurfacing::ShellPassthrough,
        js: AdapterSurfacing::ShellPassthrough,
        lua: AdapterSurfacing::ShellPassthrough,
        clawbio: AdapterSurfacing::NotApplicable,
    };
    apply_parity_matrix_overrides(CapabilitySource::EngineOperation, op_name, &mut surfacing);
    surfacing
}

fn mcp_tool_surfacing(tool_name: &str) -> CapabilitySurfacingSet {
    let mut surfacing = CapabilitySurfacingSet {
        gui: AdapterSurfacing::NotApplicable,
        cli: AdapterSurfacing::NotApplicable,
        mcp: AdapterSurfacing::Prominent,
        js: AdapterSurfacing::NotApplicable,
        lua: AdapterSurfacing::NotApplicable,
        clawbio: AdapterSurfacing::NotApplicable,
    };
    apply_parity_matrix_overrides(CapabilitySource::McpTool, tool_name, &mut surfacing);
    surfacing
}

fn apply_parity_matrix_overrides(
    source: CapabilitySource,
    name: &str,
    surfacing: &mut CapabilitySurfacingSet,
) {
    for adapter in capability_parity_adapters() {
        if not_applicable_override_reason(source, name, *adapter).is_none() {
            continue;
        }
        let current = surfacing.for_adapter(*adapter);
        assert!(
            !current.is_reachable(),
            "parity matrix override for {:?} `{}` on {:?} would demote reachable surfacing {:?}",
            source,
            name,
            adapter,
            current
        );
        surfacing.set_for_adapter(*adapter, AdapterSurfacing::NotApplicable);
    }
}

fn surfacing_justifications_for_glossary_command(
    surfacing: &CapabilitySurfacingSet,
    command: &GlossaryRegistryCommand,
) -> BTreeMap<String, String> {
    let mut justifications = BTreeMap::new();
    for adapter in capability_parity_adapters() {
        if surfacing.for_adapter(*adapter) == AdapterSurfacing::NotApplicable {
            let note = not_applicable_override_reason(
                CapabilitySource::GlossaryCommand,
                &command.path,
                *adapter,
            )
            .unwrap_or_else(|| {
                panic!(
                    "missing curated n/a override for glossary command `{}` on {:?}",
                    command.path, adapter
                )
            });
            justifications.insert(adapter.as_str().to_string(), note);
        }
    }
    justifications
}

fn surfacing_justifications_for_engine_operation(
    surfacing: &CapabilitySurfacingSet,
    op_name: &str,
) -> BTreeMap<String, String> {
    capability_parity_adapters()
        .iter()
        .filter(|&adapter| surfacing.for_adapter(*adapter) == AdapterSurfacing::NotApplicable)
        .map(|adapter| {
            (
                adapter.as_str().to_string(),
                not_applicable_override_reason(
                    CapabilitySource::EngineOperation,
                    op_name,
                    *adapter,
                )
                .unwrap_or_else(|| {
                    panic!(
                        "missing curated n/a override for engine operation `{op_name}` on {:?}",
                        adapter
                    )
                }),
            )
        })
        .collect()
}

fn surfacing_justifications_for_mcp_tool(
    surfacing: &CapabilitySurfacingSet,
    tool_name: &str,
) -> BTreeMap<String, String> {
    capability_parity_adapters()
        .iter()
        .filter(|&adapter| surfacing.for_adapter(*adapter) == AdapterSurfacing::NotApplicable)
        .map(|adapter| {
            (
                adapter.as_str().to_string(),
                not_applicable_override_reason(CapabilitySource::McpTool, tool_name, *adapter)
                    .unwrap_or_else(|| {
                        panic!(
                            "missing curated n/a override for MCP tool `{tool_name}` on {:?}",
                            adapter
                        )
                    }),
            )
        })
        .collect()
}

fn generic_output_schema() -> Value {
    json!({
        "type": "object",
        "description": "Adapter-specific structured output payload."
    })
}

fn title_from_stable_name(name: &str) -> String {
    name.replace(['-', '_'], " ")
        .split_whitespace()
        .map(|part| {
            let mut chars = part.chars();
            match chars.next() {
                Some(first) => first.to_ascii_uppercase().to_string() + chars.as_str(),
                None => String::new(),
            }
        })
        .collect::<Vec<_>>()
        .join(" ")
}

fn infer_command_mutation(path: &str, operations: &[String]) -> CapabilityMutation {
    if operations
        .iter()
        .any(|operation| infer_engine_operation_mutation(operation) == CapabilityMutation::External)
        || path.starts_with("agents ask")
        || path.starts_with("agents models")
        || path.starts_with("agents plan")
        || path.contains(" blast-start")
        || path.contains(" fetch")
    {
        CapabilityMutation::External
    } else if operations
        .iter()
        .any(|operation| infer_engine_operation_mutation(operation) == CapabilityMutation::Mutating)
        || path.starts_with("save")
        || path.starts_with("load")
        || path.starts_with("import")
        || path.starts_with("export")
        || path.contains(" materialize")
        || path.contains(" write")
        || path.contains(" set-")
        || path.contains(" sync")
    {
        CapabilityMutation::Mutating
    } else {
        CapabilityMutation::ReadOnly
    }
}

fn infer_engine_operation_mutation(operation: &str) -> CapabilityMutation {
    if operation.starts_with("Fetch")
        || operation.starts_with("PrepareGenome")
        || operation.starts_with("PrepareCutRun")
        || operation.starts_with("ReadAcquire")
        || operation.starts_with("Sync")
        || operation.starts_with("Benchmark")
        || operation.starts_with("InterpretCutRun")
        || operation.starts_with("InterpretRna")
        || operation.starts_with("RunRna")
        || operation.starts_with("ImportGenomeBigWig")
    {
        CapabilityMutation::External
    } else if operation.starts_with("Render")
        || operation.starts_with("Export")
        || operation.starts_with("List")
        || operation.starts_with("Show")
        || operation.starts_with("Inspect")
        || operation.starts_with("Query")
        || operation.starts_with("Summarize")
        || operation.starts_with("Compute")
        || operation.starts_with("Validate")
        || operation.starts_with("Audit")
        || operation.starts_with("Resolve")
        || operation.starts_with("Recommend")
        || operation == "SaveFile"
        || operation == "FindRestrictionSites"
        || operation == "AlignSequences"
        || operation == "AssessPrimerPairSpecificity"
        || operation == "TestCdnaPcr"
        || operation == "TestCdnaQpcr"
        || operation == "TestCdnaQpcrFasta"
        || operation == "SuggestSequencingPrimers"
        || operation == "SuggestPromoterReporterFragments"
        || operation == "BuildRepeatEnvironmentCohort"
        || operation == "BuildProteinToDnaHandoffReasoning"
    {
        CapabilityMutation::ReadOnly
    } else {
        CapabilityMutation::Mutating
    }
}

fn inline_operand_ok_for_operations(operations: &[String]) -> Option<bool> {
    let mut saw_known = false;
    for operation in operations {
        if let Some(value) = inline_operand_ok_for_operation(operation) {
            saw_known = true;
            if value {
                return Some(true);
            }
        }
    }
    saw_known.then_some(false)
}

fn inline_operand_ok_for_operation(operation: &str) -> Option<bool> {
    match operation {
        "RenderTfbsScoreTracksSvg"
        | "FindRestrictionSites"
        | "SummarizeTfbsScoreTracks"
        | "SummarizeTfbsTrackSimilarity"
        | "ScanTfbsHits"
        | "AlignSequences" => Some(true),
        "RenderSequenceSvg"
        | "RenderRnaStructureSvg"
        | "RenderTfbsScoreTrackCorrelationSvg"
        | "ComputeDotplot"
        | "ComputeFlexibilityTrack" => Some(false),
        _ => None,
    }
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

impl RnaReadSeedFilterConfig {
    /// Stable machine token for how input reads are interpreted before scoring.
    pub fn input_orientation_mode(&self) -> &'static str {
        if self.cdna_poly_t_flip_enabled {
            "cdna_oriented"
        } else {
            "direct_rna"
        }
    }

    /// Short human label matching [`Self::input_orientation_mode`].
    pub fn input_orientation_label(&self) -> &'static str {
        if self.cdna_poly_t_flip_enabled {
            "cDNA-oriented"
        } else {
            "direct-RNA"
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

pub const RNA_READ_ALIGNMENT_DISPLAY_BATCH_SCHEMA: &str =
    "gentle.rna_read_alignment_display_batch.v1";

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
pub struct RnaReadAlignmentDisplayBatchEntry {
    pub record_index: usize,
    pub header_id: String,
    pub alignment: RnaReadAlignmentDisplay,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
pub struct RnaReadAlignmentDisplayBatchSkippedRecord {
    pub record_index: usize,
    pub header_id: String,
    pub reason: String,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
pub struct RnaReadAlignmentDisplayBatch {
    pub schema: String,
    pub report_id: String,
    pub seq_id: String,
    pub seed_feature_id: usize,
    pub requested_gene_ids: Vec<String>,
    pub matched_gene_ids: Vec<String>,
    pub missing_gene_ids: Vec<String>,
    pub selection_mode: String,
    #[serde(default)]
    pub cohort_filter: RnaReadGeneSupportAuditCohortFilter,
    #[serde(default)]
    pub complete_rule: RnaReadGeneSupportCompleteRule,
    pub selected_record_indices: Vec<usize>,
    pub limit: Option<usize>,
    pub entry_count: usize,
    pub skipped_record_indices: Vec<usize>,
    pub skipped_records: Vec<RnaReadAlignmentDisplayBatchSkippedRecord>,
    pub entries: Vec<RnaReadAlignmentDisplayBatchEntry>,
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
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub dominant_triage_bin: Option<RnaReadIsoformTriageBin>,
    #[serde(default, skip_serializing_if = "BTreeMap::is_empty")]
    pub triage_bin_counts: BTreeMap<String, usize>,
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

pub const RNA_READ_BATCH_MAP_REPORT_SCHEMA: &str = "gentle.rna_read_batch_map_report.v1";
pub const RNA_READ_GENE_SCREEN_SUMMARY_SCHEMA: &str = "gentle.rna_read_gene_screen_summary.v1";

#[derive(Debug, Clone, Copy, Serialize, Deserialize, PartialEq, Eq, Default)]
#[serde(rename_all = "snake_case")]
pub enum RnaReadBatchMapSampleStatus {
    #[default]
    Ok,
    Failed,
    NeedsPreparation,
}

impl RnaReadBatchMapSampleStatus {
    pub fn as_str(self) -> &'static str {
        match self {
            Self::Ok => "ok",
            Self::Failed => "failed",
            Self::NeedsPreparation => "needs_preparation",
        }
    }
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
pub struct RnaReadBatchMapSampleRow {
    pub sample_id: String,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub sample_name: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub sample_description: Option<String>,
    pub status: RnaReadBatchMapSampleStatus,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub input_path: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub sra_accession: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub report_id: Option<String>,
    pub seq_id: String,
    pub seed_feature_id: usize,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub elapsed_ms: Option<u128>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub error: Option<String>,
    #[serde(default)]
    pub warnings: Vec<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub gene_support_summary_json_path: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub gene_support_audit_json_path: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub concatemer_json_path: Option<String>,
    pub read_count_total: usize,
    pub read_count_seed_passed: usize,
    pub read_count_aligned: usize,
    pub seed_pass_fraction: f64,
    pub aligned_fraction: f64,
    pub mean_read_length_bp: f64,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub all_q0_bp: Option<usize>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub all_q25_bp: Option<usize>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub all_q50_bp: Option<usize>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub all_q75_bp: Option<usize>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub all_q90_bp: Option<usize>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub all_q95_bp: Option<usize>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub all_q99_bp: Option<usize>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub all_q100_bp: Option<usize>,
    #[serde(default)]
    pub origin_class_counts: BTreeMap<String, usize>,
    #[serde(default)]
    pub requested_gene_ids: Vec<String>,
    #[serde(default)]
    pub matched_gene_ids: Vec<String>,
    #[serde(default)]
    pub missing_gene_ids: Vec<String>,
    pub aligned_base_count: usize,
    pub accepted_target_count: usize,
    pub accepted_target_fraction_total: f64,
    pub accepted_target_fraction_aligned: f64,
    pub aligned_other_gene_count: usize,
    pub aligned_other_gene_fraction_aligned: f64,
    pub fragment_count: usize,
    pub complete_count: usize,
    pub complete_strict_count: usize,
    pub complete_exact_count: usize,
    pub mean_assigned_read_length_bp: f64,
    pub isoform_support_count: usize,
    pub concatemer_inspected_count: usize,
    pub concatemer_suspicious_count: usize,
    pub concatemer_strong_count: usize,
    pub concatemer_multi_gene_fragment_count: usize,
    pub target_partner_gene_fragment_count: usize,
    pub internal_adapter_match_count: usize,
    pub disjoint_secondary_mapping_count: usize,
    pub low_primary_coverage_count: usize,
    pub internal_poly_a_count: usize,
    pub internal_poly_t_count: usize,
    pub phase1_partial_origin_count: usize,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub seed_passed_q90_bp: Option<usize>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub seed_passed_q95_bp: Option<usize>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub seed_passed_q99_bp: Option<usize>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub seed_passed_max_bp: Option<usize>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub seed_passed_mean_bp: Option<f64>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub accepted_target_max_bp: Option<usize>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub accepted_target_mean_bp: Option<f64>,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
pub struct RnaReadBatchIsoformSupportRow {
    pub sample_id: String,
    pub report_id: String,
    pub seq_id: String,
    pub gene_id: String,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub transcript_feature_id: Option<usize>,
    pub transcript_id: String,
    pub transcript_label: String,
    pub aligned_count: usize,
    pub fragment_count: usize,
    pub complete_count: usize,
    pub complete_strict_count: usize,
    pub complete_exact_count: usize,
    pub mean_read_length_bp: f64,
    pub mean_identity_fraction: f64,
    pub mean_query_coverage_fraction: f64,
    pub exon_support_json: String,
    pub exon_pair_support_json: String,
    pub direct_transition_support_json: String,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
pub struct RnaReadBatchConcatemerPartnerRow {
    pub sample_id: String,
    pub report_id: String,
    pub seq_id: String,
    pub partner_kind: String,
    pub gene_id: String,
    pub transcript_id: String,
    pub transcript_label: String,
    pub suspicious_read_count: usize,
    pub fragment_count: usize,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
pub struct RnaReadBatchMapSraPreparationRow {
    pub sample_id: String,
    pub sra_accession: String,
    pub planned_fasta_path: String,
    pub preparation_command: String,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
pub struct RnaReadBatchMapReport {
    pub schema: String,
    pub manifest_path: String,
    pub out_dir: String,
    pub generated_at_unix_ms: u128,
    pub seq_id: String,
    pub seed_feature_id: usize,
    #[serde(default)]
    pub requested_gene_ids: Vec<String>,
    #[serde(default)]
    pub target_gene_ids: Vec<String>,
    pub profile: String,
    pub input_format: String,
    pub scope: String,
    pub origin_mode: String,
    pub report_mode: String,
    pub align_selection: String,
    pub complete_rule: String,
    pub max_secondary_mappings: usize,
    pub continue_on_error: bool,
    pub batch_report_json_path: String,
    pub batch_summary_tsv_path: String,
    pub gene_screen_summary_tsv_path: String,
    pub sample_sheet_tsv_path: String,
    pub isoform_support_tsv_path: String,
    pub concatemer_partner_summary_tsv_path: String,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub sra_preparation_plan_tsv_path: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub sra_preparation_commands_sh_path: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub read_acquisition_report: Option<ReadAcquisitionReport>,
    pub sample_count: usize,
    pub ok_count: usize,
    pub failed_count: usize,
    pub needs_preparation_count: usize,
    #[serde(default)]
    pub rows: Vec<RnaReadBatchMapSampleRow>,
    #[serde(default)]
    pub isoform_support_rows: Vec<RnaReadBatchIsoformSupportRow>,
    #[serde(default)]
    pub concatemer_partner_rows: Vec<RnaReadBatchConcatemerPartnerRow>,
    #[serde(default)]
    pub sra_preparation_rows: Vec<RnaReadBatchMapSraPreparationRow>,
    #[serde(default)]
    pub warnings: Vec<String>,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
pub struct RnaReadTargetQualityComparisonEntry {
    pub entry_id: String,
    pub comparison_label: String,
    pub gentle_version: String,
    pub report_id: String,
    pub seq_id: String,
    pub report_generated_at_unix_ms: u128,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub report_op_id: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub report_run_id: Option<String>,
    #[serde(default)]
    pub requested_gene_ids: Vec<String>,
    #[serde(default)]
    pub matched_gene_ids: Vec<String>,
    #[serde(default)]
    pub missing_gene_ids: Vec<String>,
    #[serde(default)]
    pub report_target_gene_ids: Vec<String>,
    #[serde(default)]
    pub complete_rule: RnaReadGeneSupportCompleteRule,
    #[serde(default)]
    pub profile: RnaReadInterpretationProfile,
    #[serde(default)]
    pub scope: SplicingScopePreset,
    #[serde(default)]
    pub origin_mode: RnaReadOriginMode,
    #[serde(default)]
    pub report_mode: RnaReadReportMode,
    #[serde(default)]
    pub seed_filter: RnaReadSeedFilterConfig,
    #[serde(default)]
    pub align_config: RnaReadAlignConfig,
    #[serde(default)]
    pub all_read_lengths: RnaReadLengthDistributionSummary,
    pub summary: RnaReadGeneSupportSummary,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
pub struct RnaReadTargetQualityComparisonBundle {
    pub schema: String,
    pub generated_at_unix_ms: u128,
    #[serde(default)]
    pub entries: Vec<RnaReadTargetQualityComparisonEntry>,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
pub struct RnaReadTargetQualityExport {
    pub schema: String,
    pub requested_path: String,
    pub written_path: String,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub bundle_path: Option<String>,
    pub format: String,
    pub report_id: String,
    #[serde(default)]
    pub requested_gene_ids: Vec<String>,
    #[serde(default)]
    pub complete_rule: RnaReadGeneSupportCompleteRule,
    pub entry_count: usize,
    #[serde(default)]
    pub appended_to_existing_bundle: bool,
    #[serde(default)]
    pub reused_existing_entry_slot: bool,
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
pub struct RnaReadLengthDistributionSummary {
    pub sample_count: usize,
    pub mean_length_bp: f64,
    pub min_length_bp: usize,
    pub q25_length_bp: usize,
    pub median_length_bp: usize,
    pub q75_length_bp: usize,
    pub max_length_bp: usize,
    pub p95_length_bp: usize,
    pub length_counts: Vec<u64>,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
pub struct RnaReadFractionDistributionSummary {
    pub sample_count: usize,
    pub mean_fraction: f64,
    pub median_fraction: f64,
    pub p95_fraction: f64,
    pub bin_counts: Vec<u64>,
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
    pub evaluated_read_lengths: RnaReadLengthDistributionSummary,
    pub accepted_target_read_lengths: RnaReadLengthDistributionSummary,
    pub accepted_target_fragment_lengths: RnaReadLengthDistributionSummary,
    pub accepted_target_query_coverage: RnaReadFractionDistributionSummary,
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
pub struct CutRunCatalogEntry {
    pub description: Option<String>,
    pub summary: Option<String>,
    #[serde(default)]
    pub aliases: Vec<String>,
    #[serde(default)]
    pub tags: Vec<String>,
    #[serde(default)]
    pub search_terms: Vec<String>,
    pub species: Option<String>,
    pub assembly_label: Option<String>,
    pub target_factor: Option<String>,
    pub sample_label: Option<String>,
    pub tissue_or_cell_type: Option<String>,
    pub condition: Option<String>,
    pub replicate: Option<String>,
    pub assay_kind: Option<String>,
    #[serde(default)]
    pub supported_reference_genome_ids: Vec<String>,
    pub provider: Option<String>,
    pub source_accession: Option<String>,
    pub reference_url: Option<String>,
    pub peaks_remote: Option<String>,
    pub peaks_local: Option<String>,
    pub signal_remote: Option<String>,
    pub signal_local: Option<String>,
    pub reads_r1_remote: Option<String>,
    pub reads_r1_local: Option<String>,
    pub reads_r2_remote: Option<String>,
    pub reads_r2_local: Option<String>,
    pub reads_sra_accession: Option<String>,
    #[serde(default)]
    pub read_layout: CutRunReadLayout,
    pub cache_dir: Option<String>,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
pub struct CutRunCatalogListEntry {
    pub dataset_id: String,
    pub description: Option<String>,
    pub summary: Option<String>,
    #[serde(default)]
    pub aliases: Vec<String>,
    #[serde(default)]
    pub tags: Vec<String>,
    #[serde(default)]
    pub search_terms: Vec<String>,
    pub species: Option<String>,
    pub assembly_label: Option<String>,
    pub target_factor: Option<String>,
    pub sample_label: Option<String>,
    pub tissue_or_cell_type: Option<String>,
    pub condition: Option<String>,
    pub replicate: Option<String>,
    pub assay_kind: Option<String>,
    #[serde(default)]
    pub supported_reference_genome_ids: Vec<String>,
    pub provider: Option<String>,
    pub source_accession: Option<String>,
    pub reference_url: Option<String>,
    pub has_peaks_asset: bool,
    pub has_signal_asset: bool,
    pub has_raw_reads: bool,
    #[serde(default)]
    pub read_layout: CutRunReadLayout,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
pub struct CutRunDatasetListReport {
    pub schema: String,
    pub catalog_origin_label: String,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub requested_catalog_path: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub filter: Option<String>,
    pub dataset_count: usize,
    pub datasets: Vec<CutRunCatalogListEntry>,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
pub struct CutRunPreparedAssetManifest {
    pub source: String,
    pub local_path: String,
    pub file_name: String,
    pub file_size_bytes: u64,
    pub checksum_sha1: String,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
pub struct CutRunPreparedManifest {
    pub schema: String,
    pub dataset_id: String,
    pub prepared_at_unix_ms: u128,
    pub catalog_origin_label: String,
    pub install_dir: String,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub peaks: Option<CutRunPreparedAssetManifest>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub signal: Option<CutRunPreparedAssetManifest>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub reads_r1: Option<CutRunPreparedAssetManifest>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub reads_r2: Option<CutRunPreparedAssetManifest>,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
pub struct CutRunPreparedAssetStatus {
    pub configured: bool,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub source: Option<String>,
    pub prepared: bool,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub local_path: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub file_name: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub file_size_bytes: Option<u64>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub checksum_sha1: Option<String>,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
/// Shared lease/heartbeat status for one prepared shared asset install.
pub struct SharedAssetActivityStatus {
    pub resource_key: String,
    pub display_name: String,
    pub status_path: String,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub lock_path: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub cancel_path: Option<String>,
    pub lifecycle_status: String,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub phase: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub item: Option<String>,
    #[serde(default)]
    pub bytes_done: u64,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub bytes_total: Option<u64>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub percent: Option<f64>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub monitored_free_bytes: Option<u64>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub minimum_free_bytes: Option<u64>,
    pub started_at_unix_ms: u128,
    pub updated_at_unix_ms: u128,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub finished_at_unix_ms: Option<u128>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub last_error: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub owner_pid: Option<u32>,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
/// Portable readiness/status snapshot for one CUT&RUN prepared dataset.
pub struct CutRunDatasetStatus {
    pub schema: String,
    pub dataset_id: String,
    pub resource_key: String,
    pub display_name: String,
    pub catalog_origin_label: String,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub requested_catalog_path: Option<String>,
    pub effective_cache_dir: String,
    pub install_dir: String,
    pub prepared: bool,
    pub lifecycle_status: String,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub current_activity: Option<SharedAssetActivityStatus>,
    pub description: Option<String>,
    pub summary: Option<String>,
    pub species: Option<String>,
    pub assembly_label: Option<String>,
    pub target_factor: Option<String>,
    pub sample_label: Option<String>,
    pub tissue_or_cell_type: Option<String>,
    pub condition: Option<String>,
    pub replicate: Option<String>,
    pub assay_kind: Option<String>,
    #[serde(default)]
    pub supported_reference_genome_ids: Vec<String>,
    pub provider: Option<String>,
    pub source_accession: Option<String>,
    pub reference_url: Option<String>,
    #[serde(default)]
    pub read_layout: CutRunReadLayout,
    pub peaks: CutRunPreparedAssetStatus,
    pub signal: CutRunPreparedAssetStatus,
    pub reads_r1: CutRunPreparedAssetStatus,
    pub reads_r2: CutRunPreparedAssetStatus,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub manifest_path: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub manifest: Option<CutRunPreparedManifest>,
    #[serde(default)]
    pub warnings: Vec<String>,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
pub struct CutRunDatasetProjectionReport {
    pub schema: String,
    pub seq_id: String,
    pub dataset_id: String,
    pub include_peaks: bool,
    pub include_signal: bool,
    pub clear_existing: bool,
    pub projected_peak_features: usize,
    pub projected_signal_features: usize,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub peak_track_name: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub signal_track_name: Option<String>,
    #[serde(default)]
    pub warnings: Vec<String>,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
pub struct CutRunReadPlacement {
    pub orientation: CutRunReadOrientation,
    pub local_start_1based: usize,
    pub local_end_1based: usize,
    pub genomic_start_1based: usize,
    pub genomic_end_1based: usize,
    pub mismatches: usize,
    pub identity_fraction: f64,
    pub seed_matches: usize,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
pub struct CutRunReadUnitRow {
    pub unit_index: usize,
    pub normalized_read_id: String,
    pub status: CutRunReadUnitStatus,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub r1_record_index: Option<usize>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub r1_header_id: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub r1_read_length_bp: Option<usize>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub r1_placement: Option<CutRunReadPlacement>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub r2_record_index: Option<usize>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub r2_header_id: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub r2_read_length_bp: Option<usize>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub r2_placement: Option<CutRunReadPlacement>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub fragment_local_start_1based: Option<usize>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub fragment_local_end_1based: Option<usize>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub fragment_genomic_start_1based: Option<usize>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub fragment_genomic_end_1based: Option<usize>,
    #[serde(default)]
    pub deduplicated: bool,
    #[serde(default)]
    pub duplicate_count: usize,
    #[serde(default)]
    pub notes: Vec<String>,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
pub struct CutRunFragmentSpan {
    pub fragment_id: String,
    pub normalized_read_id: String,
    pub status: CutRunReadUnitStatus,
    pub local_start_1based: usize,
    pub local_end_1based: usize,
    pub genomic_start_1based: usize,
    pub genomic_end_1based: usize,
    pub length_bp: usize,
    pub left_cut_site_local_1based: usize,
    pub right_cut_site_local_1based: usize,
    #[serde(default)]
    pub duplicate_count: usize,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
pub struct CutRunSupportCluster {
    pub cluster_index: usize,
    pub local_start_1based: usize,
    pub local_end_1based: usize,
    pub genomic_start_1based: usize,
    pub genomic_end_1based: usize,
    pub peak_coverage: u32,
    pub total_cut_sites: u32,
    pub fragment_count: usize,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
pub struct CutRunReadReportSummary {
    pub report_id: String,
    pub seq_id: String,
    pub generated_at_unix_ms: u128,
    pub input_format: CutRunInputFormat,
    pub read_layout: CutRunReadLayout,
    pub roi_flank_bp: usize,
    pub reference_window_length_bp: usize,
    pub total_units: usize,
    pub mapped_units: usize,
    pub fragment_count: usize,
    pub concordant_pair_count: usize,
    pub orphan_unit_count: usize,
    pub mean_read_length_bp: f64,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
pub struct CutRunReadReport {
    pub schema: String,
    pub report_id: String,
    pub seq_id: String,
    pub generated_at_unix_ms: u128,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub op_id: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub run_id: Option<String>,
    pub input_r1_path: String,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub input_r2_path: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub dataset_id: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub target_factor: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub species: Option<String>,
    pub input_format: CutRunInputFormat,
    pub read_layout: CutRunReadLayout,
    pub roi_flank_bp: usize,
    pub deduplicate_fragments: bool,
    pub seed_filter: CutRunSeedFilterConfig,
    pub align_config: CutRunAlignConfig,
    pub genome_id: String,
    pub chromosome: String,
    pub reference_window_start_1based: usize,
    pub reference_window_end_1based: usize,
    pub reference_window_orientation: String,
    pub roi_local_start_1based: usize,
    pub roi_local_end_1based: usize,
    pub reference_window_length_bp: usize,
    pub total_units: usize,
    pub mapped_units: usize,
    pub fragment_count: usize,
    pub concordant_pair_count: usize,
    pub orphan_r1_count: usize,
    pub orphan_r2_count: usize,
    pub unmatched_pair_count: usize,
    pub mean_read_length_bp: f64,
    pub mean_mapped_read_length_bp: f64,
    #[serde(default)]
    pub coverage: Vec<u32>,
    #[serde(default)]
    pub cut_site_counts: Vec<u32>,
    #[serde(default)]
    pub units: Vec<CutRunReadUnitRow>,
    #[serde(default)]
    pub fragments: Vec<CutRunFragmentSpan>,
    #[serde(default)]
    pub support_clusters: Vec<CutRunSupportCluster>,
    #[serde(default)]
    pub warnings: Vec<String>,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
pub struct CutRunReadReportStore {
    pub schema: String,
    pub updated_at_unix_ms: u128,
    #[serde(default)]
    pub reports: HashMap<String, CutRunReadReport>,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
pub struct CutRunReadCoverageExport {
    pub schema: String,
    pub path: String,
    pub report_id: String,
    pub kind: CutRunCoverageKind,
    pub row_count: usize,
}

#[derive(Debug, Clone, Copy, Serialize, Deserialize, PartialEq, Eq, Default)]
#[serde(rename_all = "snake_case")]
pub enum CutRunRegulatoryEvidenceSourceKind {
    #[default]
    Dataset,
    ReadReport,
}

#[derive(Debug, Clone, Copy, Serialize, Deserialize, PartialEq, Eq, Default)]
#[serde(rename_all = "snake_case")]
pub enum CutRunSupportStrength {
    #[default]
    Weak,
    Moderate,
    Strong,
}

#[derive(Debug, Clone, Copy, Serialize, Deserialize, PartialEq, Eq, Default)]
#[serde(rename_all = "snake_case")]
pub enum CutRunRegulatoryTfbsConfirmationStatus {
    #[default]
    Unconfirmed,
    Confirmed,
    Nearby,
    Absent,
    MotifPoor,
}

#[derive(Debug, Clone, Copy, Serialize, Deserialize, PartialEq, Eq, Default)]
#[serde(rename_all = "snake_case")]
pub enum CutRunMotifContextScope {
    #[default]
    InsideWindow,
    NeighborWindow,
}

#[derive(Debug, Clone, Copy, Serialize, Deserialize, PartialEq, Eq, Default)]
#[serde(rename_all = "snake_case")]
pub enum CutRunMotifAbsentOccupancyInterpretation {
    #[default]
    ContextSupportedByOtherMotifs,
    MotifPoorSupported,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
pub struct CutRunRegulatoryEvidenceSourceRef {
    pub source_kind: CutRunRegulatoryEvidenceSourceKind,
    pub source_id: String,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub description: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub dataset_id: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub report_id: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub target_factor: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub species: Option<String>,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
pub struct CutRunSupportWindowRecord {
    pub window_id: String,
    pub local_start_0based: usize,
    pub local_end_0based_exclusive: usize,
    pub genomic_start_1based: usize,
    pub genomic_end_1based: usize,
    pub support_strength: CutRunSupportStrength,
    pub overlapping_peak_count: usize,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub max_signal_value: Option<f64>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub mean_signal_value: Option<f64>,
    pub supporting_fragment_count: usize,
    pub cut_site_count: u32,
    #[serde(default)]
    pub contributing_dataset_ids: Vec<String>,
    #[serde(default)]
    pub contributing_read_report_ids: Vec<String>,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
pub struct CutRunRegulatoryTfbsRow {
    pub feature_id: usize,
    pub feature_label: String,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub motif_id: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub motif_label: Option<String>,
    pub local_start_0based: usize,
    pub local_end_0based_exclusive: usize,
    pub genomic_start_1based: usize,
    pub genomic_end_1based: usize,
    pub strand: String,
    pub confirmation_status: CutRunRegulatoryTfbsConfirmationStatus,
    pub support_status: CutRunRegulatoryTfbsConfirmationStatus,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub strongest_support_window_id: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub strongest_support_strength: Option<CutRunSupportStrength>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub support_distance_bp: Option<usize>,
    pub overlapping_peak_count: usize,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub max_signal_value: Option<f64>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub mean_signal_value: Option<f64>,
    pub supporting_fragment_count: usize,
    pub cut_site_count: u32,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
pub struct CutRunMotifContextHit {
    pub motif_id: String,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub motif_label: Option<String>,
    pub hit_count: usize,
    pub best_llr_bits: f64,
    pub best_true_log_odds_bits: f64,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
pub struct CutRunMotifAbsentSupportWindow {
    pub window_id: String,
    pub local_start_0based: usize,
    pub local_end_0based_exclusive: usize,
    pub genomic_start_1based: usize,
    pub genomic_end_1based: usize,
    pub support_strength: CutRunSupportStrength,
    pub overlapping_peak_count: usize,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub max_signal_value: Option<f64>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub mean_signal_value: Option<f64>,
    pub supporting_fragment_count: usize,
    pub cut_site_count: u32,
    pub target_motif_resolved: bool,
    pub target_motif_present: bool,
    #[serde(default)]
    pub motifs_inside_window: Vec<CutRunMotifContextHit>,
    #[serde(default)]
    pub motifs_in_neighbor_window: Vec<CutRunMotifContextHit>,
    pub occupancy_interpretation: CutRunMotifAbsentOccupancyInterpretation,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
pub struct CutRunMotifContextSummaryRow {
    pub motif_id: String,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub motif_label: Option<String>,
    pub context_scope: CutRunMotifContextScope,
    pub window_count: usize,
    pub window_fraction: f64,
    pub mean_best_score: f64,
    pub max_best_score: f64,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
pub struct CutRunRegulatorySupportReport {
    pub schema: String,
    pub seq_id: String,
    pub generated_at_unix_ms: u128,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub op_id: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub run_id: Option<String>,
    #[serde(default)]
    pub evidence_sources: Vec<CutRunRegulatoryEvidenceSourceRef>,
    pub promoter_search_start_0based: usize,
    pub promoter_search_end_0based_exclusive: usize,
    pub neighbor_window_bp: usize,
    pub motif_context_min_llr_quantile: f64,
    #[serde(default)]
    pub species_filters: Vec<String>,
    #[serde(default)]
    pub support_windows: Vec<CutRunSupportWindowRecord>,
    #[serde(default)]
    pub confirmed_tfbs_rows: Vec<CutRunRegulatoryTfbsRow>,
    #[serde(default)]
    pub unconfirmed_tfbs_rows: Vec<CutRunRegulatoryTfbsRow>,
    #[serde(default)]
    pub motif_absent_supported_windows: Vec<CutRunMotifAbsentSupportWindow>,
    #[serde(default)]
    pub common_motifs_inside_supported_windows: Vec<CutRunMotifContextSummaryRow>,
    #[serde(default)]
    pub common_motifs_near_supported_windows: Vec<CutRunMotifContextSummaryRow>,
    #[serde(default)]
    pub warnings: Vec<String>,
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

/// Conservative per-read isoform triage bins for aligned RNA-read reports.
#[derive(Debug, Clone, Copy, Serialize, Deserialize, PartialEq, Eq, PartialOrd, Ord, Default)]
#[serde(rename_all = "snake_case")]
pub enum RnaReadIsoformTriageBin {
    KnownIsoformConfirmed,
    KnownIsoformAmbiguous,
    GeneSupportedNoIsoformCall,
    #[default]
    OffTargetOrBadSeed,
}

impl RnaReadIsoformTriageBin {
    pub fn as_str(self) -> &'static str {
        match self {
            Self::KnownIsoformConfirmed => "known_isoform_confirmed",
            Self::KnownIsoformAmbiguous => "known_isoform_ambiguous",
            Self::GeneSupportedNoIsoformCall => "gene_supported_no_isoform_call",
            Self::OffTargetOrBadSeed => "off_target_or_bad_seed",
        }
    }
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
pub struct RnaReadIsoformTriageTsvExport {
    pub schema: String,
    pub path: String,
    pub report_id: String,
    pub selection: RnaReadHitSelection,
    pub row_count: usize,
    pub limit: Option<usize>,
    pub min_identity_fraction: f64,
    pub min_query_coverage_fraction: f64,
    pub min_confirmed_transition_fraction: f64,
    pub max_secondary_mappings: usize,
    pub bin_counts: BTreeMap<String, usize>,
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

fn default_rna_read_concatemer_internal_homopolymer_min_bp() -> usize {
    18
}

fn default_rna_read_concatemer_end_margin_bp() -> usize {
    30
}

fn default_rna_read_concatemer_max_primary_query_coverage_fraction() -> f64 {
    0.85
}

fn default_rna_read_concatemer_min_secondary_identity_fraction() -> f64 {
    0.85
}

fn default_rna_read_concatemer_max_secondary_query_overlap_fraction() -> f64 {
    0.20
}

fn default_rna_read_concatemer_adapter_min_match_bp() -> usize {
    16
}

fn default_rna_read_concatemer_fragment_min_bp() -> usize {
    60
}

fn default_rna_read_concatemer_fragment_max_parts() -> usize {
    4
}

fn default_rna_read_concatemer_fragment_min_identity_fraction() -> f64 {
    0.80
}

fn default_rna_read_concatemer_fragment_min_query_coverage_fraction() -> f64 {
    0.35
}

fn deserialize_rna_read_concatemer_transcript_fasta_paths<'de, D>(
    deserializer: D,
) -> Result<Vec<String>, D::Error>
where
    D: Deserializer<'de>,
{
    #[derive(Deserialize)]
    #[serde(untagged)]
    enum OneOrMany {
        One(String),
        Many(Vec<String>),
    }

    let value = Option::<OneOrMany>::deserialize(deserializer)?;
    Ok(match value {
        Some(OneOrMany::One(path)) => vec![path],
        Some(OneOrMany::Many(paths)) => paths,
        None => Vec::new(),
    })
}

#[derive(Debug, Clone, Copy, Serialize, Deserialize, PartialEq, Eq, Default)]
#[serde(rename_all = "snake_case")]
/// Conservative severity label for fragment/concatemer-like RNA-read rows.
pub enum RnaReadConcatemerSuspicionLevel {
    #[default]
    Background,
    Weak,
    Moderate,
    Strong,
}

impl RnaReadConcatemerSuspicionLevel {
    pub fn as_str(self) -> &'static str {
        match self {
            Self::Background => "background",
            Self::Weak => "weak",
            Self::Moderate => "moderate",
            Self::Strong => "strong",
        }
    }
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(default)]
/// Shared thresholds for the RNA-read fragment/concatemer suspicion audit.
pub struct RnaReadConcatemerInspectionSettings {
    #[serde(default = "default_rna_read_concatemer_internal_homopolymer_min_bp")]
    pub internal_homopolymer_min_bp: usize,
    #[serde(default = "default_rna_read_concatemer_end_margin_bp")]
    pub end_margin_bp: usize,
    #[serde(default = "default_rna_read_concatemer_max_primary_query_coverage_fraction")]
    pub max_primary_query_coverage_fraction: f64,
    #[serde(default = "default_rna_read_concatemer_min_secondary_identity_fraction")]
    pub min_secondary_identity_fraction: f64,
    #[serde(default = "default_rna_read_concatemer_max_secondary_query_overlap_fraction")]
    pub max_secondary_query_overlap_fraction: f64,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub adapter_fasta_path: Option<String>,
    #[serde(default = "default_rna_read_concatemer_adapter_min_match_bp")]
    pub adapter_min_match_bp: usize,
    #[serde(default = "default_rna_read_concatemer_fragment_min_bp")]
    pub fragment_min_bp: usize,
    #[serde(default = "default_rna_read_concatemer_fragment_max_parts")]
    pub fragment_max_parts: usize,
    #[serde(default = "default_rna_read_concatemer_fragment_min_identity_fraction")]
    pub fragment_min_identity_fraction: f64,
    #[serde(default = "default_rna_read_concatemer_fragment_min_query_coverage_fraction")]
    pub fragment_min_query_coverage_fraction: f64,
    #[serde(
        default,
        alias = "transcript_fasta_path",
        deserialize_with = "deserialize_rna_read_concatemer_transcript_fasta_paths",
        skip_serializing_if = "Vec::is_empty"
    )]
    pub transcript_fasta_paths: Vec<String>,
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub transcript_index_paths: Vec<String>,
}

impl Default for RnaReadConcatemerInspectionSettings {
    fn default() -> Self {
        Self {
            internal_homopolymer_min_bp: default_rna_read_concatemer_internal_homopolymer_min_bp(),
            end_margin_bp: default_rna_read_concatemer_end_margin_bp(),
            max_primary_query_coverage_fraction:
                default_rna_read_concatemer_max_primary_query_coverage_fraction(),
            min_secondary_identity_fraction:
                default_rna_read_concatemer_min_secondary_identity_fraction(),
            max_secondary_query_overlap_fraction:
                default_rna_read_concatemer_max_secondary_query_overlap_fraction(),
            adapter_fasta_path: None,
            adapter_min_match_bp: default_rna_read_concatemer_adapter_min_match_bp(),
            fragment_min_bp: default_rna_read_concatemer_fragment_min_bp(),
            fragment_max_parts: default_rna_read_concatemer_fragment_max_parts(),
            fragment_min_identity_fraction:
                default_rna_read_concatemer_fragment_min_identity_fraction(),
            fragment_min_query_coverage_fraction:
                default_rna_read_concatemer_fragment_min_query_coverage_fraction(),
            transcript_fasta_paths: Vec::new(),
            transcript_index_paths: Vec::new(),
        }
    }
}

pub const RNA_READ_TRANSCRIPT_CATALOG_INDEX_SCHEMA: &str =
    "gentle.rna_read_transcript_catalog_index.v1";

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
/// One normalized transcript template row stored in a reusable RNA-read
/// transcript catalog index.
pub struct RnaReadTranscriptCatalogTemplateRecord {
    pub transcript_id: String,
    pub transcript_label: String,
    pub gene_id: String,
    pub strand: String,
    pub sequence: String,
    pub kmer_positions: BTreeMap<u32, Vec<usize>>,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
/// Reusable transcript catalog prepared for concatemer/fragment-origin audits.
pub struct RnaReadTranscriptCatalogIndex {
    pub schema: String,
    pub generated_at_unix_ms: u128,
    pub seed_kmer_len: usize,
    pub source_paths: Vec<String>,
    pub transcript_count: usize,
    pub gene_count: usize,
    pub warnings: Vec<String>,
    pub templates: Vec<RnaReadTranscriptCatalogTemplateRecord>,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
/// Seed-gate score for one pseudo-read used by the RNA-read isoform preflight.
pub struct RnaReadIsoformPreflightScore {
    pub transcript_id: String,
    pub transcript_label: String,
    pub sequence_length_bp: usize,
    pub passed_seed_filter: bool,
    pub raw_hit_fraction: f64,
    pub weighted_hit_fraction: f64,
    pub unique_matched_kmers: usize,
    pub chain_consistency_fraction: f64,
    pub seed_median_transcript_gap: f64,
    pub confirmed_transitions: usize,
    pub total_transitions: usize,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
/// Per-gene/per-family control summary emitted by RNA-read isoform preflight.
pub struct RnaReadIsoformPreflightControlSummary {
    pub control_id: String,
    pub source_paths: Vec<String>,
    pub transcript_count: usize,
    pub passed_transcript_count: usize,
    pub weighted_pass_probability: f64,
    pub best_transcript_id: Option<String>,
    pub best_transcript_label: Option<String>,
    pub best_raw_hit_fraction: f64,
    pub best_weighted_hit_fraction: f64,
    pub best_unique_matched_kmers: usize,
    pub best_chain_consistency_fraction: f64,
    pub best_seed_median_transcript_gap: f64,
    pub best_confirmed_transitions: usize,
    pub best_total_transitions: usize,
    pub worst_case_ambiguity: f64,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
/// Paste-ready seed-filter thresholds derived from target-vs-control isoform
/// preflight.
pub struct RnaReadIsoformPreflightThresholdRecommendation {
    pub basis: String,
    pub target_pass_probability: f64,
    pub positive_pass_probability: f64,
    pub max_control_pass_probability: f64,
    pub limiting_control_id: Option<String>,
    pub limiting_control_transcript_id: Option<String>,
    pub limiting_control_transcript_label: Option<String>,
    pub control_raw_hit_margin: f64,
    pub control_weighted_hit_margin: f64,
    pub control_unique_kmer_margin: isize,
    pub seed_filter_cli_fragment: String,
    pub interpret_command_fragment: String,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
/// Upfront RNA-read isoform seed-gate preflight over target transcripts and
/// explicit external control transcript variants.
pub struct RnaReadIsoformPreflightReport {
    pub schema: String,
    pub seq_id: String,
    pub seed_feature_id: usize,
    pub scope: SplicingScopePreset,
    pub seed_filter: RnaReadSeedFilterConfig,
    pub optimize_parameters: bool,
    pub positive_transcript_fasta_paths: Vec<String>,
    pub control_transcript_fasta_paths: Vec<String>,
    pub max_control_match_probability: f64,
    pub target_transcript_count: usize,
    pub target_passed_transcript_count: usize,
    pub target_pass_probability: f64,
    pub positive_control_transcript_count: usize,
    pub positive_control_passed_transcript_count: usize,
    pub positive_control_pass_probability: f64,
    pub target_transcripts: Vec<RnaReadIsoformPreflightScore>,
    pub positive_control_transcripts: Vec<RnaReadIsoformPreflightScore>,
    pub control_summaries: Vec<RnaReadIsoformPreflightControlSummary>,
    pub recommended_seed_filter: RnaReadSeedFilterConfig,
    pub threshold_recommendation: RnaReadIsoformPreflightThresholdRecommendation,
    pub recommended_command_fragment: String,
    pub warnings: Vec<String>,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
/// One internal adapter-signature match found away from the read ends.
pub struct RnaReadConcatemerAdapterHit {
    pub label: String,
    pub orientation: String,
    pub start_0based: usize,
    pub end_0based_exclusive: usize,
    pub matched_bp: usize,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
/// One iteratively recovered transcript-origin fragment from a suspicious read.
pub struct RnaReadConcatemerFragmentOrigin {
    pub fragment_rank: usize,
    pub query_start_0based: usize,
    pub query_end_0based_exclusive: usize,
    pub query_length_bp: usize,
    pub transcript_feature_id: usize,
    pub transcript_id: String,
    pub transcript_label: String,
    pub gene_id: String,
    pub strand: String,
    pub identity_fraction: f64,
    pub query_coverage_fraction: f64,
    pub score: isize,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
/// Aggregated non-primary gene/group partners recurring across suspicious reads.
pub struct RnaReadConcatemerPartnerGeneSummary {
    pub rank: usize,
    pub gene_id: String,
    pub suspicious_read_count: usize,
    pub fragment_count: usize,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
/// Aggregated non-primary transcript partners recurring across suspicious reads.
pub struct RnaReadConcatemerPartnerTranscriptSummary {
    pub rank: usize,
    pub transcript_id: String,
    pub transcript_label: String,
    pub gene_id: String,
    pub suspicious_read_count: usize,
    pub fragment_count: usize,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
/// One retained RNA-read row ranked by fragment/concatemer suspicion.
pub struct RnaReadConcatemerSuspicionRow {
    pub rank: usize,
    pub record_index: usize,
    pub header_id: String,
    pub read_length_bp: usize,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub phase1_primary_transcript_id: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub best_transcript_id: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub best_transcript_label: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub best_gene_id: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub best_identity_fraction: Option<f64>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub best_query_coverage_fraction: Option<f64>,
    #[serde(default)]
    pub secondary_mapping_count: usize,
    #[serde(default)]
    pub internal_poly_a_run_bp: usize,
    #[serde(default)]
    pub internal_poly_t_run_bp: usize,
    #[serde(default)]
    pub internal_adapter_hit_count: usize,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub top_internal_adapter_label: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub top_internal_adapter_match_bp: Option<usize>,
    #[serde(default)]
    pub disjoint_secondary_mapping_count: usize,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub top_disjoint_secondary_transcript_id: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub top_disjoint_secondary_identity_fraction: Option<f64>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub top_disjoint_secondary_query_coverage_fraction: Option<f64>,
    #[serde(default)]
    pub origin_class: RnaReadOriginClass,
    #[serde(default)]
    pub suspicion_score: usize,
    #[serde(default)]
    pub suspicion_level: RnaReadConcatemerSuspicionLevel,
    #[serde(default)]
    pub suspicion_signals: Vec<String>,
    #[serde(default)]
    pub suspicion_summary: String,
    #[serde(default)]
    pub fragment_origin_gene_count: usize,
    #[serde(default)]
    pub fragment_origin_gene_ids: Vec<String>,
    #[serde(default)]
    pub partner_gene_count: usize,
    #[serde(default)]
    pub partner_gene_ids: Vec<String>,
    #[serde(default)]
    pub partner_transcript_count: usize,
    #[serde(default)]
    pub partner_transcript_ids: Vec<String>,
    #[serde(default)]
    pub adapter_hits: Vec<RnaReadConcatemerAdapterHit>,
    #[serde(default)]
    pub fragment_origins: Vec<RnaReadConcatemerFragmentOrigin>,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
/// Shared fragment/concatemer suspicion audit over one RNA-read report.
pub struct RnaReadConcatemerInspection {
    pub schema: String,
    pub report_id: String,
    pub seq_id: String,
    pub selection: RnaReadHitSelection,
    #[serde(default)]
    pub selected_record_indices: Vec<usize>,
    pub subset_match_count: usize,
    pub inspected_count: usize,
    pub suspicious_count: usize,
    pub strong_count: usize,
    pub low_query_coverage_count: usize,
    pub internal_poly_a_count: usize,
    pub internal_poly_t_count: usize,
    pub internal_adapter_match_count: usize,
    pub disjoint_secondary_mapping_count: usize,
    pub phase1_partial_origin_count: usize,
    pub multi_gene_fragment_count: usize,
    pub limit: usize,
    pub max_secondary_mappings: usize,
    pub settings: RnaReadConcatemerInspectionSettings,
    #[serde(default)]
    pub warnings: Vec<String>,
    #[serde(default)]
    pub partner_gene_summaries: Vec<RnaReadConcatemerPartnerGeneSummary>,
    #[serde(default)]
    pub partner_transcript_summaries: Vec<RnaReadConcatemerPartnerTranscriptSummary>,
    #[serde(default)]
    pub rows: Vec<RnaReadConcatemerSuspicionRow>,
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
    OrganismEcoliDefault,
    OrganellePlastidDefault,
    OrganelleVertebrateMitochondrialDefault,
    OrganelleInvertebrateMitochondrialDefault,
    OrganelleYeastMitochondrialDefault,
    AmbiguousMitochondrialDefault,
}

impl TranscriptProteinTranslationTableSource {
    pub fn as_str(self) -> &'static str {
        match self {
            Self::StandardDefault => "standard_default",
            Self::ExplicitCdsQualifier => "explicit_cds_qualifier",
            Self::ExplicitTranscriptQualifier => "explicit_transcript_qualifier",
            Self::ExplicitSourceQualifier => "explicit_source_qualifier",
            Self::OrganismEcoliDefault => "organism_ecoli_default",
            Self::OrganellePlastidDefault => "organelle_plastid_default",
            Self::OrganelleVertebrateMitochondrialDefault => {
                "organelle_vertebrate_mitochondrial_default"
            }
            Self::OrganelleInvertebrateMitochondrialDefault => {
                "organelle_invertebrate_mitochondrial_default"
            }
            Self::OrganelleYeastMitochondrialDefault => "organelle_yeast_mitochondrial_default",
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

#[derive(Debug, Clone, Copy, Serialize, Deserialize, PartialEq, Eq, Default)]
#[serde(rename_all = "snake_case")]
/// How a translation-speed profile was resolved.
pub enum TranslationSpeedProfileSource {
    #[default]
    SourceOrganismScientificName,
    SourceOrganismCommonAlias,
    FeatureQualifierHint,
    ExplicitRequest,
}

impl TranslationSpeedProfileSource {
    pub fn as_str(self) -> &'static str {
        match self {
            Self::SourceOrganismScientificName => "source_organism_scientific_name",
            Self::SourceOrganismCommonAlias => "source_organism_common_alias",
            Self::FeatureQualifierHint => "feature_qualifier_hint",
            Self::ExplicitRequest => "explicit_request",
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

#[derive(Debug, Clone, Copy, Serialize, Deserialize, PartialEq, Eq, Default)]
#[serde(rename_all = "snake_case")]
/// Output selection for protein-feature coding DNA queries.
pub enum UniprotFeatureCodingDnaQueryMode {
    GenomicAsEncoded,
    TranslationSpeedOptimized,
    #[default]
    Both,
}

impl UniprotFeatureCodingDnaQueryMode {
    pub fn as_str(self) -> &'static str {
        match self {
            Self::GenomicAsEncoded => "genomic_as_encoded",
            Self::TranslationSpeedOptimized => "translation_speed_optimized",
            Self::Both => "both",
        }
    }
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
/// One coding segment contributing to a queried protein feature.
pub struct UniprotFeatureCodingDnaSegment {
    pub aa_start: usize,
    pub aa_end: usize,
    pub genomic_start_1based: usize,
    pub genomic_end_1based: usize,
    pub strand: String,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
/// One exon overlap contributing coding DNA to a queried protein feature.
pub struct UniprotFeatureCodingDnaExonSpan {
    pub exon_ordinal: usize,
    pub exon_start_1based: usize,
    pub exon_end_1based: usize,
    pub coding_start_1based: usize,
    pub coding_end_1based: usize,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
/// Ordered exon transition crossed by a queried protein feature.
pub struct UniprotFeatureCodingDnaExonPair {
    pub from_exon_ordinal: usize,
    pub to_exon_ordinal: usize,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
/// One transcript-specific answer for a queried UniProt feature.
pub struct UniprotFeatureCodingDnaMatch {
    pub transcript_id: String,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub transcript_feature_id: Option<usize>,
    pub strand: String,
    pub feature_key: String,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub feature_note: Option<String>,
    #[serde(default)]
    pub matched_feature_fields: Vec<String>,
    pub aa_start: usize,
    pub aa_end: usize,
    pub amino_acid_sequence: String,
    #[serde(default)]
    pub genomic_segments: Vec<UniprotFeatureCodingDnaSegment>,
    pub genomic_coding_dna: String,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub translation_speed_optimized_dna: Option<String>,
    #[serde(default)]
    pub exon_spans: Vec<UniprotFeatureCodingDnaExonSpan>,
    #[serde(default)]
    pub exon_pairs: Vec<UniprotFeatureCodingDnaExonPair>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub primary_exon_ordinal: Option<usize>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub primary_exon_pair: Option<UniprotFeatureCodingDnaExonPair>,
    #[serde(default)]
    pub crosses_exon_junction: bool,
    #[serde(default)]
    pub warnings: Vec<String>,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
/// Portable response for querying coding DNA behind one projected UniProt feature.
pub struct UniprotFeatureCodingDnaQueryReport {
    pub schema: String,
    pub projection_id: String,
    pub entry_id: String,
    pub seq_id: String,
    pub feature_query: String,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub transcript_filter: Option<String>,
    #[serde(default)]
    pub query_mode: UniprotFeatureCodingDnaQueryMode,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub requested_translation_speed_profile: Option<TranslationSpeedProfile>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub resolved_translation_speed_profile: Option<TranslationSpeedProfile>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub resolved_translation_speed_profile_source: Option<TranslationSpeedProfileSource>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub resolved_translation_speed_reference_species: Option<String>,
    pub match_count: usize,
    #[serde(default)]
    pub matches: Vec<UniprotFeatureCodingDnaMatch>,
    #[serde(default)]
    pub warnings: Vec<String>,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
/// One nucleotide of a queried protein-residue codon, reported in coding
/// transcript order with its corresponding genomic coordinate.
pub struct ProteinResidueGenomicCoordinateBase {
    /// Zero-based codon position within the residue: `0`, `1`, or `2`.
    pub codon_offset_0based: usize,
    /// One-based genomic coordinate on the source sequence.
    pub genomic_pos_1based: usize,
    pub base: String,
    /// One-based exon ordinal in transcript order.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub exon_ordinal: Option<usize>,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
/// One transcript-specific genomic coordinate answer for a protein residue.
pub struct ProteinResidueGenomicCoordinateMatch {
    pub transcript_id: String,
    pub transcript_label: String,
    /// Zero-based source feature index for the matched transcript/mRNA feature.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub transcript_feature_id: Option<usize>,
    pub strand: String,
    /// One-based residue index within the derived protein sequence.
    pub residue_index_1based: usize,
    pub amino_acid: String,
    pub codon: String,
    /// Minimum one-based genomic coordinate across the three codon bases.
    pub genomic_codon_start_1based: usize,
    /// Maximum one-based genomic coordinate across the three codon bases.
    pub genomic_codon_end_1based: usize,
    pub spans_exon_junction: bool,
    #[serde(default)]
    pub genomic_bases: Vec<ProteinResidueGenomicCoordinateBase>,
    #[serde(default)]
    pub warnings: Vec<String>,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
/// Portable response for mapping transcript-native protein residues back to
/// genomic codon nucleotide coordinates.
pub struct ProteinResidueGenomicCoordinateReport {
    pub schema: String,
    pub seq_id: SeqId,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub transcript_filter: Option<String>,
    /// One-based inclusive residue-range start requested by the caller.
    pub residue_start_1based: usize,
    /// One-based inclusive residue-range end requested by the caller.
    pub residue_end_1based: usize,
    pub match_count: usize,
    #[serde(default)]
    pub matches: Vec<ProteinResidueGenomicCoordinateMatch>,
    #[serde(default)]
    pub warnings: Vec<String>,
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
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub translation_speed_profile_source: Option<TranslationSpeedProfileSource>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub translation_speed_reference_species: Option<String>,
    #[serde(default)]
    pub terminal_stop_trimmed: bool,
    #[serde(default)]
    pub warnings: Vec<String>,
}

#[derive(Debug, Clone, Copy, Serialize, Deserialize, PartialEq, Eq, Default)]
#[serde(rename_all = "snake_case")]
/// Optional external protein-evidence sources compared against
/// transcript-native translation.
pub enum ProteinExternalOpinionSource {
    #[default]
    Uniprot,
    Ensembl,
    Other,
}

impl ProteinExternalOpinionSource {
    pub fn as_str(self) -> &'static str {
        match self {
            Self::Uniprot => "uniprot",
            Self::Ensembl => "ensembl",
            Self::Other => "other",
        }
    }
}

#[derive(Debug, Clone, Copy, Serialize, Deserialize, PartialEq, Eq, Default)]
#[serde(rename_all = "snake_case")]
/// Comparison status between transcript-native translation and one optional
/// external protein opinion.
pub enum TranscriptProteinComparisonStatus {
    #[default]
    DerivedOnly,
    ConsistentWithExternalOpinion,
    LowConfidenceExternalOpinion,
    NoTranscriptCds,
    ExternalOpinionOnly,
}

impl TranscriptProteinComparisonStatus {
    pub fn as_str(self) -> &'static str {
        match self {
            Self::DerivedOnly => "derived_only",
            Self::ConsistentWithExternalOpinion => "consistent_with_external_opinion",
            Self::LowConfidenceExternalOpinion => "low_confidence_external_opinion",
            Self::NoTranscriptCds => "no_transcript_cds",
            Self::ExternalOpinionOnly => "external_opinion_only",
        }
    }
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
/// One optional external protein interpretation layered onto a
/// transcript-native product.
pub struct TranscriptProteinExternalOpinion {
    #[serde(default)]
    pub source: ProteinExternalOpinionSource,
    pub source_id: String,
    pub source_label: String,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub expected_length_aa: Option<usize>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub reference_start_aa: Option<usize>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub reference_end_aa: Option<usize>,
    #[serde(default)]
    pub genomic_coding_ranges_1based: Vec<(usize, usize)>,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
/// Source-neutral transcript protein comparison record used by shared expert
/// views.
pub struct TranscriptProteinComparison {
    pub transcript_id: String,
    pub transcript_label: String,
    #[serde(default)]
    pub status: TranscriptProteinComparisonStatus,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub derived: Option<TranscriptProteinDerivation>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub external_opinion: Option<TranscriptProteinExternalOpinion>,
    #[serde(default)]
    pub mismatch_reasons: Vec<String>,
    #[serde(default)]
    pub derived_only_exon_ranges_1based: Vec<(usize, usize)>,
    #[serde(default)]
    pub external_only_exon_ranges_1based: Vec<(usize, usize)>,
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
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub op_id: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub run_id: Option<String>,
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
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub op_id: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub run_id: Option<String>,
    pub profile: RnaReadInterpretationProfile,
    pub input_path: String,
    pub input_format: RnaReadInputFormat,
    pub seed_feature_id: usize,
    pub scope: SplicingScopePreset,
    #[serde(default)]
    pub origin_mode: RnaReadOriginMode,
    #[serde(default)]
    pub input_orientation_mode: String,
    #[serde(default)]
    pub input_orientation_label: String,
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

#[cfg(test)]
mod dotplot_and_concatemer_setting_tests {
    use super::{
        CutRunRegulatoryTfbsConfirmationStatus, CutRunRegulatoryTfbsRow, DotplotMode,
        DotplotOverlayAnchorExonRef, DotplotOverlayResolvedAnchorSeries, DotplotOverlayXAxisMode,
        DotplotQuerySeries, DotplotView, RestrictionSiteExpertView,
        RnaReadConcatemerInspectionSettings, RnaReadMappedIsoformSupportRow,
    };

    #[test]
    fn dotplot_overlay_x_axis_bp_alignment_projects_left_and_right_variants() {
        let left_fraction = DotplotOverlayXAxisMode::LeftAlignedBp.point_fraction(4, 0, 8, 12);
        let right_fraction = DotplotOverlayXAxisMode::RightAlignedBp.point_fraction(4, 0, 8, 12);
        assert!(left_fraction < right_fraction);

        assert_eq!(
            DotplotOverlayXAxisMode::LeftAlignedBp.query_coordinate_at_fraction(0.2, 0, 8, 12),
            Some(2)
        );
        assert_eq!(
            DotplotOverlayXAxisMode::LeftAlignedBp.query_coordinate_at_fraction(0.95, 0, 8, 12),
            None
        );
        assert_eq!(
            DotplotOverlayXAxisMode::RightAlignedBp.query_coordinate_at_fraction(0.05, 0, 8, 12),
            None
        );
        assert_eq!(
            DotplotOverlayXAxisMode::RightAlignedBp.query_coordinate_at_fraction(0.82, 0, 8, 12),
            Some(5)
        );
    }

    #[test]
    fn cutrun_tfbs_row_deserializes_old_two_state_payload_without_support_status() {
        let row: CutRunRegulatoryTfbsRow = serde_json::from_value(serde_json::json!({
            "feature_id": 7,
            "feature_label": "SP1 site",
            "local_start_0based": 10,
            "local_end_0based_exclusive": 20,
            "genomic_start_1based": 110,
            "genomic_end_1based": 119,
            "strand": "+",
            "confirmation_status": "confirmed",
            "overlapping_peak_count": 1,
            "supporting_fragment_count": 2,
            "cut_site_count": 4
        }))
        .expect("deserialize old CUT&RUN TFBS row");
        assert_eq!(
            row.confirmation_status,
            CutRunRegulatoryTfbsConfirmationStatus::Confirmed
        );
        assert_eq!(
            row.support_status,
            CutRunRegulatoryTfbsConfirmationStatus::Unconfirmed
        );
    }

    #[test]
    fn dotplot_overlay_anchor_exon_ref_parses_range_token() {
        let exon = DotplotOverlayAnchorExonRef::parse("27..34").expect("parse anchor exon");
        assert_eq!(exon.start_1based, 27);
        assert_eq!(exon.end_1based, 34);
        assert_eq!(exon.token(), "27..34");
    }

    #[test]
    fn dotplot_view_resolves_manual_query_anchor_series() {
        let view = DotplotView {
            query_series: vec![
                DotplotQuerySeries {
                    series_id: "tp73".to_string(),
                    seq_id: "tp73".to_string(),
                    label: "TP73".to_string(),
                    color_rgb: [29, 78, 216],
                    transcript_feature_id: None,
                    query_anchor_0based: Some(926),
                    query_anchor_label: Some("shared core motif".to_string()),
                    mode: DotplotMode::PairForward,
                    span_start_0based: 0,
                    span_end_0based: 1200,
                    point_count: 0,
                    points: vec![],
                    boxplot_bin_count: 0,
                    boxplot_bins: vec![],
                },
                DotplotQuerySeries {
                    series_id: "tp63".to_string(),
                    seq_id: "tp63".to_string(),
                    label: "TP63".to_string(),
                    color_rgb: [220, 38, 38],
                    transcript_feature_id: None,
                    query_anchor_0based: Some(1044),
                    query_anchor_label: Some("shared core motif".to_string()),
                    mode: DotplotMode::PairForward,
                    span_start_0based: 0,
                    span_end_0based: 1300,
                    point_count: 0,
                    points: vec![],
                    boxplot_bin_count: 0,
                    boxplot_bins: vec![],
                },
                DotplotQuerySeries {
                    series_id: "tp53".to_string(),
                    seq_id: "tp53".to_string(),
                    label: "TP53".to_string(),
                    color_rgb: [5, 150, 105],
                    transcript_feature_id: None,
                    query_anchor_0based: Some(707),
                    query_anchor_label: Some("shared core motif".to_string()),
                    mode: DotplotMode::PairForward,
                    span_start_0based: 0,
                    span_end_0based: 950,
                    point_count: 0,
                    points: vec![],
                    boxplot_bin_count: 0,
                    boxplot_bins: vec![],
                },
            ],
            ..DotplotView::default()
        };
        let resolved = view.resolve_query_anchor_series();
        assert_eq!(resolved.len(), 3);
        assert_eq!(
            resolved,
            vec![
                DotplotOverlayResolvedAnchorSeries {
                    series_index: 0,
                    shift_bp: 118,
                    plotted_span_end_0based: 1318,
                },
                DotplotOverlayResolvedAnchorSeries {
                    series_index: 1,
                    shift_bp: 0,
                    plotted_span_end_0based: 1300,
                },
                DotplotOverlayResolvedAnchorSeries {
                    series_index: 2,
                    shift_bp: 337,
                    plotted_span_end_0based: 1287,
                },
            ]
        );
        assert_eq!(view.query_anchor_label(), Some("shared core motif"));
    }

    #[test]
    fn rna_read_concatemer_settings_accept_legacy_single_transcript_fasta_path() {
        let settings: RnaReadConcatemerInspectionSettings =
            serde_json::from_str("{\"transcript_fasta_path\":\"data/transcriptome.fa.gz\"}")
                .expect("deserialize legacy concatemer settings");
        assert_eq!(
            settings.transcript_fasta_paths,
            vec!["data/transcriptome.fa.gz".to_string()]
        );
    }

    #[test]
    fn rna_read_concatemer_settings_accept_multiple_transcript_fasta_paths() {
        let settings: RnaReadConcatemerInspectionSettings = serde_json::from_str(
            "{\"transcript_fasta_paths\":[\"data/cdna.fa.gz\",\"data/ncrna.fa.gz\"]}",
        )
        .expect("deserialize plural concatemer settings");
        assert_eq!(
            settings.transcript_fasta_paths,
            vec![
                "data/cdna.fa.gz".to_string(),
                "data/ncrna.fa.gz".to_string()
            ]
        );
    }

    #[test]
    fn mapped_isoform_support_row_accepts_legacy_json_without_triage_fields() {
        let row: RnaReadMappedIsoformSupportRow = serde_json::from_str(
            r#"{
                "transcript_feature_id": 7,
                "transcript_id": "TX1",
                "transcript_label": "isoform 1",
                "strand": "+",
                "aligned_read_count": 3
            }"#,
        )
        .expect("legacy mapped isoform row should deserialize");

        assert_eq!(row.transcript_id, "TX1");
        assert_eq!(row.aligned_read_count, 3);
        assert!(row.dominant_triage_bin.is_none());
        assert!(row.triage_bin_counts.is_empty());
    }

    #[test]
    fn restriction_site_expert_view_builds_tooltip_lines_with_cut_markers() {
        let view = RestrictionSiteExpertView {
            seq_id: "p".to_string(),
            cut_pos_1based: 4,
            paired_cut_pos_1based: 8,
            recognition_start_1based: 3,
            recognition_end_1based: 8,
            cut_index_0based: 1,
            paired_cut_index_0based: 5,
            end_geometry: "5prime_overhang".to_string(),
            number_of_cuts_for_enzyme: 1,
            selected_enzyme: Some("EcoRI".to_string()),
            enzyme_names: vec!["EcoRI".to_string()],
            recognition_iupac: Some("GAATTC".to_string()),
            site_sequence: "GAATTC".to_string(),
            site_sequence_complement: "CTTAAG".to_string(),
            enzyme_cut_offset_0based: Some(1),
            overlap_bp: Some(4),
            enzyme_note: None,
            rebase_url: Some("https://rebase.neb.com/rebase/enz/EcoRI.html".to_string()),
            tooltip_lines: vec![],
            instruction: "inspect".to_string(),
        };

        assert_eq!(view.enzyme_display_label(), "EcoRI");
        assert_eq!(view.site_count_label(), "1 site");
        assert_eq!(view.cut_position_label(), "cuts at 4|8 bp");
        assert_eq!(view.geometry_display_label(), "5' overhang (4 bp)");
        assert_eq!(view.marked_top_sequence(), "G^AATTC");
        assert_eq!(view.marked_bottom_sequence(), "CTTAA^G");
        let lines = view.tooltip_summary_lines();
        assert!(lines.iter().any(|line| line.contains("EcoRI | 1 site")));
        assert!(lines.iter().any(|line| line.contains("5' G^AATTC 3'")));
        assert!(lines.iter().any(|line| line.contains("3' CTTAA^G 5'")));
    }
}
