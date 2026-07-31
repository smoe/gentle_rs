//! Portable ortholog and cross-species promoter comparison contracts.
//!
//! Ortholog resources are local, reviewable mapping tables. The engine resolves
//! them into promoter windows through prepared genome catalogs, then emits
//! evidence-separated comparison reports for sequence, TFBS, expression, and
//! CUT&RUN/occupancy signals.

use crate::{
    BiologicalContextRegistry, BiologicalContextResolutionError, GeneSetCohortRelationship,
    GeneSetCohortRelationshipFlag,
};
use serde::{Deserialize, Deserializer, Serialize, Serializer};
use std::{collections::BTreeMap, fmt, ops::Deref};

/// Local offline ortholog mapping resource schema.
pub const ORTHOLOG_RESOURCE_SCHEMA: &str = "gentle.ortholog_resource.v1";
/// Resolved cross-species promoter cohort report schema.
pub const ORTHOLOG_PROMOTER_COHORT_SCHEMA: &str = "gentle.ortholog_promoter_cohort.v1";
/// Cross-species promoter comparison report schema.
pub const ORTHOLOG_PROMOTER_COMPARISON_SCHEMA: &str = "gentle.ortholog_promoter_comparison.v1";

/// Canonical cardinality represented by a recognized orthology-type value.
#[derive(Debug, Clone, Copy, Serialize, Deserialize, PartialEq, Eq)]
#[serde(rename_all = "snake_case")]
pub enum OrthologyCardinality {
    OneToOne,
    OneToMany,
    ManyToOne,
    ManyToMany,
}

impl OrthologyCardinality {
    pub fn as_str(self) -> &'static str {
        match self {
            Self::OneToOne => "one_to_one",
            Self::OneToMany => "one_to_many",
            Self::ManyToOne => "many_to_one",
            Self::ManyToMany => "many_to_many",
        }
    }
}

/// Open orthology-type vocabulary that preserves provider-specific legacy text.
///
/// Known spellings expose a typed cardinality while unknown values round-trip
/// unchanged. This keeps existing local resources such as `synthetic_ortholog`
/// valid instead of treating a new vocabulary as a migration gate.
#[derive(Debug, Clone, Default, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct OrthologyType(String);

impl OrthologyType {
    pub const ONE_TO_ONE: &'static str = "one_to_one";
    pub const ONE_TO_MANY: &'static str = "one_to_many";
    pub const MANY_TO_ONE: &'static str = "many_to_one";
    pub const MANY_TO_MANY: &'static str = "many_to_many";

    pub fn new(value: impl Into<String>) -> Self {
        Self(value.into())
    }

    pub fn as_str(&self) -> &str {
        &self.0
    }

    pub fn cardinality(&self) -> Option<OrthologyCardinality> {
        match normalized_orthology_token(&self.0).as_str() {
            "one_to_one" | "one2one" | "ortholog_one_to_one" | "ortholog_one2one" => {
                Some(OrthologyCardinality::OneToOne)
            }
            "one_to_many" | "one2many" | "ortholog_one_to_many" | "ortholog_one2many" => {
                Some(OrthologyCardinality::OneToMany)
            }
            "many_to_one" | "many2one" | "ortholog_many_to_one" | "ortholog_many2one" => {
                Some(OrthologyCardinality::ManyToOne)
            }
            "many_to_many" | "many2many" | "ortholog_many_to_many" | "ortholog_many2many" => {
                Some(OrthologyCardinality::ManyToMany)
            }
            _ => None,
        }
    }
}

impl From<String> for OrthologyType {
    fn from(value: String) -> Self {
        Self::new(value)
    }
}

impl From<&str> for OrthologyType {
    fn from(value: &str) -> Self {
        Self::new(value)
    }
}

impl Deref for OrthologyType {
    type Target = str;

    fn deref(&self) -> &Self::Target {
        self.as_str()
    }
}

impl fmt::Display for OrthologyType {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        formatter.write_str(self.as_str())
    }
}

impl Serialize for OrthologyType {
    fn serialize<S>(&self, serializer: S) -> Result<S::Ok, S::Error>
    where
        S: Serializer,
    {
        serializer.serialize_str(self.as_str())
    }
}

impl<'de> Deserialize<'de> for OrthologyType {
    fn deserialize<D>(deserializer: D) -> Result<Self, D::Error>
    where
        D: Deserializer<'de>,
    {
        String::deserialize(deserializer).map(Self)
    }
}

/// Canonical confidence tier represented by a recognized confidence value.
#[derive(Debug, Clone, Copy, Serialize, Deserialize, PartialEq, Eq)]
#[serde(rename_all = "snake_case")]
pub enum OrthologConfidenceLevel {
    High,
    Medium,
    Low,
}

impl OrthologConfidenceLevel {
    pub fn as_str(self) -> &'static str {
        match self {
            Self::High => "high",
            Self::Medium => "medium",
            Self::Low => "low",
        }
    }
}

/// Open ortholog-confidence vocabulary preserving provider-specific values.
#[derive(Debug, Clone, Default, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct OrthologConfidence(String);

impl OrthologConfidence {
    pub const HIGH: &'static str = "high";
    pub const MEDIUM: &'static str = "medium";
    pub const LOW: &'static str = "low";

    pub fn new(value: impl Into<String>) -> Self {
        Self(value.into())
    }

    pub fn as_str(&self) -> &str {
        &self.0
    }

    pub fn level(&self) -> Option<OrthologConfidenceLevel> {
        match normalized_orthology_token(&self.0).as_str() {
            "high" => Some(OrthologConfidenceLevel::High),
            "medium" | "moderate" => Some(OrthologConfidenceLevel::Medium),
            "low" => Some(OrthologConfidenceLevel::Low),
            _ => None,
        }
    }
}

impl From<String> for OrthologConfidence {
    fn from(value: String) -> Self {
        Self::new(value)
    }
}

impl From<&str> for OrthologConfidence {
    fn from(value: &str) -> Self {
        Self::new(value)
    }
}

impl Deref for OrthologConfidence {
    type Target = str;

    fn deref(&self) -> &Self::Target {
        self.as_str()
    }
}

impl fmt::Display for OrthologConfidence {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        formatter.write_str(self.as_str())
    }
}

impl Serialize for OrthologConfidence {
    fn serialize<S>(&self, serializer: S) -> Result<S::Ok, S::Error>
    where
        S: Serializer,
    {
        serializer.serialize_str(self.as_str())
    }
}

impl<'de> Deserialize<'de> for OrthologConfidence {
    fn deserialize<D>(deserializer: D) -> Result<Self, D::Error>
    where
        D: Deserializer<'de>,
    {
        String::deserialize(deserializer).map(Self)
    }
}

fn normalized_orthology_token(value: &str) -> String {
    value.trim().to_ascii_lowercase().replace(['-', ' '], "_")
}

/// How to handle multiple local ortholog rows for one anchor/target pair.
#[derive(Debug, Clone, Copy, Serialize, Deserialize, PartialEq, Eq, Default)]
#[serde(rename_all = "snake_case")]
pub enum OrthologAmbiguityPolicy {
    #[default]
    Reject,
    First,
    Preserve,
}

impl OrthologAmbiguityPolicy {
    pub fn as_str(self) -> &'static str {
        match self {
            Self::Reject => "reject",
            Self::First => "first",
            Self::Preserve => "preserve",
        }
    }
}

/// Anchor or target row role in a resolved ortholog promoter cohort.
#[derive(Debug, Clone, Copy, Serialize, Deserialize, PartialEq, Eq, Default)]
#[serde(rename_all = "snake_case")]
pub enum OrthologPromoterRole {
    #[default]
    Anchor,
    Target,
}

/// Cross-species CUT&RUN/occupancy support status.
#[derive(Debug, Clone, Copy, Serialize, Deserialize, PartialEq, Eq, Default)]
#[serde(rename_all = "snake_case")]
pub enum OrthologCutRunSupportStatus {
    Confirmed,
    Nearby,
    MotifOnly,
    OccupancyOnly,
    #[default]
    NoData,
    NotComparable,
}

impl OrthologCutRunSupportStatus {
    pub fn as_str(self) -> &'static str {
        match self {
            Self::Confirmed => "confirmed",
            Self::Nearby => "nearby",
            Self::MotifOnly => "motif_only",
            Self::OccupancyOnly => "occupancy_only",
            Self::NoData => "no_data",
            Self::NotComparable => "not_comparable",
        }
    }
}

/// Whether normalized cross-species CUT&RUN values were compared.
#[derive(Debug, Clone, Copy, Serialize, Deserialize, PartialEq, Eq, Default)]
#[serde(rename_all = "snake_case")]
pub enum OrthologCutRunQuantitativeComparisonStatus {
    #[default]
    NotComparable,
    Comparable,
}

impl OrthologCutRunQuantitativeComparisonStatus {
    pub fn as_str(self) -> &'static str {
        match self {
            Self::NotComparable => "not_comparable",
            Self::Comparable => "comparable",
        }
    }
}

/// One caller-supplied normalized CUT&RUN value for a resolved promoter.
///
/// Values are never synthesized from prepared peak scores or read counts.
/// Source ids and provenance bind each value to the evidence from which the
/// caller derived it.
#[derive(Debug, Clone, Serialize, Deserialize, Default, PartialEq)]
#[serde(default)]
pub struct OrthologCutRunNormalizedValueInput {
    pub species: String,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub gene_label: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub transcript_id: Option<String>,
    pub normalized_value: f64,
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub contributing_dataset_ids: Vec<String>,
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub contributing_read_report_ids: Vec<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub provenance: Option<String>,
}

/// Explicit normalization contract for quantitative cross-promoter CUT&RUN.
///
/// A method, unit, shared comparison reference, and provenance statement are
/// all required before the engine will compare the supplied values.
#[derive(Debug, Clone, Serialize, Deserialize, Default, PartialEq)]
#[serde(default)]
pub struct OrthologCutRunNormalizationInput {
    pub normalization_method: String,
    pub unit: String,
    pub comparison_reference: String,
    pub provenance: String,
    #[serde(default)]
    pub values: Vec<OrthologCutRunNormalizedValueInput>,
}

/// One normalized value assigned to a resolved ortholog promoter row.
#[derive(Debug, Clone, Serialize, Deserialize, Default, PartialEq)]
#[serde(default)]
pub struct OrthologCutRunNormalizedAssignment {
    pub species: String,
    pub gene_label: String,
    pub transcript_id: String,
    pub normalized_value: f64,
    #[serde(default)]
    pub contributing_dataset_ids: Vec<String>,
    #[serde(default)]
    pub contributing_read_report_ids: Vec<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub provenance: Option<String>,
}

/// Pairwise difference between two explicitly normalized CUT&RUN values.
#[derive(Debug, Clone, Serialize, Deserialize, Default, PartialEq)]
#[serde(default)]
pub struct OrthologCutRunPairwiseQuantitativeComparison {
    pub left_species: String,
    pub right_species: String,
    pub left_gene_label: String,
    pub right_gene_label: String,
    pub left_normalized_value: f64,
    pub right_normalized_value: f64,
    pub delta_right_minus_left: f64,
    pub absolute_delta: f64,
}

/// Fail-closed quantitative CUT&RUN comparison carried beside qualitative
/// motif/occupancy states.
#[derive(Debug, Clone, Serialize, Deserialize, Default, PartialEq)]
#[serde(default)]
pub struct OrthologCutRunQuantitativeComparison {
    pub status: OrthologCutRunQuantitativeComparisonStatus,
    pub normalization: OrthologCutRunNormalizationInput,
    #[serde(default)]
    pub assignments: Vec<OrthologCutRunNormalizedAssignment>,
    #[serde(default)]
    pub pairwise_comparisons: Vec<OrthologCutRunPairwiseQuantitativeComparison>,
    pub detail: String,
}

/// Species alias mapping used by a local ortholog resource.
#[derive(Debug, Clone, Serialize, Deserialize, Default, PartialEq, Eq)]
#[serde(default)]
pub struct OrthologSpeciesAlias {
    pub species: String,
    pub aliases: Vec<String>,
}

/// One directional row in a local ortholog mapping resource.
#[derive(Debug, Clone, Serialize, Deserialize, Default, PartialEq)]
#[serde(default)]
pub struct OrthologMappingRow {
    pub source_species: String,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub source_context_id: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub source_gene_id: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub source_gene_symbol: Option<String>,
    pub target_species: String,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub target_context_id: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub target_gene_id: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub target_gene_symbol: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub orthology_type: Option<OrthologyType>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub confidence: Option<OrthologConfidence>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub source: Option<String>,
    #[serde(default)]
    pub evidence: Vec<String>,
}

/// Reviewable local ortholog mapping resource.
#[derive(Debug, Clone, Serialize, Deserialize, Default, PartialEq)]
#[serde(default)]
pub struct OrthologResource {
    pub schema: String,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub id: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub label: Option<String>,
    #[serde(flatten)]
    pub biological_contexts: BiologicalContextRegistry,
    #[serde(default)]
    pub species_aliases: Vec<OrthologSpeciesAlias>,
    #[serde(default)]
    pub rows: Vec<OrthologMappingRow>,
    #[serde(default)]
    pub warnings: Vec<String>,
}

impl OrthologResource {
    /// Validate the registry and every optional row-level context reference.
    pub fn validate_context_references(&self) -> Result<(), BiologicalContextResolutionError> {
        self.biological_contexts.validate()?;
        for row in &self.rows {
            for context_id in [
                row.source_context_id.as_deref(),
                row.target_context_id.as_deref(),
            ]
            .into_iter()
            .flatten()
            {
                self.biological_contexts.context(context_id)?;
            }
        }
        Ok(())
    }
}

/// Request echoed into a resolved ortholog promoter cohort.
#[derive(Debug, Clone, Serialize, Deserialize, Default, PartialEq, Eq)]
#[serde(default)]
pub struct OrthologPromoterCohortRequest {
    pub anchor_species: String,
    pub anchor_genome_id: String,
    pub anchor_gene_query: String,
    pub target_species: Vec<String>,
    #[serde(default)]
    pub target_genome_ids: BTreeMap<String, String>,
    #[serde(default)]
    pub transcript_ids: BTreeMap<String, String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub ortholog_resource_path: Option<String>,
    pub upstream_bp: usize,
    pub downstream_bp: usize,
    pub ambiguity_policy: OrthologAmbiguityPolicy,
    #[serde(default)]
    pub relationship: GeneSetCohortRelationship,
}

/// One ordered mapping candidate retained for an unresolved ambiguous target.
///
/// Candidate rows preserve declared identity and evidence without treating the
/// mapping as selected or resolving it into a promoter-cohort member.
#[derive(Debug, Clone, Serialize, Deserialize, Default, PartialEq, Eq)]
#[serde(default)]
pub struct OrthologAmbiguityCandidate {
    pub candidate_rank: usize,
    pub candidate_label: String,
    pub target_species: String,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub target_genome_id: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub source_context_id: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub target_context_id: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub target_gene_id: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub target_gene_symbol: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub orthology_type: Option<OrthologyType>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub confidence: Option<OrthologConfidence>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub source: Option<String>,
    #[serde(default)]
    pub evidence: Vec<String>,
}

/// One unresolved species/gene row from ortholog promoter resolution.
#[derive(Debug, Clone, Serialize, Deserialize, Default, PartialEq, Eq)]
#[serde(default)]
pub struct OrthologUnresolvedRow {
    pub species: String,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub context_id: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub genome_id: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub gene_query: Option<String>,
    pub reason: String,
    #[serde(default)]
    pub candidates: Vec<String>,
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub candidate_mappings: Vec<OrthologAmbiguityCandidate>,
}

/// One promoter window resolved for an anchor or target ortholog.
#[derive(Debug, Clone, Serialize, Deserialize, Default, PartialEq)]
#[serde(default)]
pub struct OrthologPromoterRow {
    pub species: String,
    pub genome_id: String,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub context_id: Option<String>,
    pub role: OrthologPromoterRole,
    pub gene_query: String,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub gene_symbol: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub gene_id: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub transcript_id_requested: Option<String>,
    pub transcript_id: String,
    pub display_label: String,
    pub chromosome: String,
    pub strand: String,
    pub promoter_start_1based: usize,
    pub promoter_end_1based: usize,
    pub promoter_length_bp: usize,
    pub tss_1based: usize,
    pub tss_position_0based: usize,
    pub sequence_orientation: String,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub promoter_sequence: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub orthology_type: Option<OrthologyType>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub confidence: Option<OrthologConfidence>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub orthology_source_context_id: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub orthology_target_context_id: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub orthology_source: Option<String>,
    #[serde(default)]
    pub orthology_evidence: Vec<String>,
    #[serde(default)]
    pub warnings: Vec<String>,
}

/// Portable resolved ortholog promoter cohort.
#[derive(Debug, Clone, Serialize, Deserialize, Default, PartialEq)]
#[serde(default)]
pub struct OrthologPromoterCohortReport {
    pub schema: String,
    pub generated_at_unix_ms: u128,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub op_id: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub run_id: Option<String>,
    pub request: OrthologPromoterCohortRequest,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub ortholog_resource_label: Option<String>,
    #[serde(flatten)]
    pub biological_contexts: BiologicalContextRegistry,
    pub resolved_promoter_count: usize,
    pub unresolved_count: usize,
    #[serde(default)]
    pub rows: Vec<OrthologPromoterRow>,
    #[serde(default)]
    pub unresolved_rows: Vec<OrthologUnresolvedRow>,
    #[serde(default)]
    pub relationship: GeneSetCohortRelationship,
    #[serde(default)]
    pub relationship_flags: Vec<GeneSetCohortRelationshipFlag>,
    #[serde(default)]
    pub warnings: Vec<String>,
}

/// Per-species motif evidence summary in an ortholog promoter comparison.
#[derive(Debug, Clone, Serialize, Deserialize, Default, PartialEq)]
#[serde(default)]
pub struct OrthologTfbsSummaryRow {
    pub species: String,
    pub gene_label: String,
    pub transcript_id: String,
    pub tf_id: String,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub tf_name: Option<String>,
    pub max_score: f64,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub peak_position_0based: Option<usize>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub peak_position_promoter_relative_bp: Option<i64>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub peak_genomic_position_1based: Option<usize>,
    pub positive_fraction: f64,
}

/// Pairwise TFBS score-track similarity across two ortholog promoters.
#[derive(Debug, Clone, Serialize, Deserialize, Default, PartialEq)]
#[serde(default)]
pub struct OrthologPairwiseTfbsSimilarity {
    pub left_species: String,
    pub right_species: String,
    pub left_gene_label: String,
    pub right_gene_label: String,
    pub shared_motif_count: usize,
    pub mean_raw_pearson: f64,
    pub mean_smoothed_spearman: f64,
    #[serde(default)]
    pub motif_ids: Vec<String>,
}

/// Motif peak summary shared across, or specific to, species in a cohort.
#[derive(Debug, Clone, Serialize, Deserialize, Default, PartialEq)]
#[serde(default)]
pub struct OrthologTfbsPeakSummary {
    pub tf_id: String,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub tf_name: Option<String>,
    pub promoter_count: usize,
    #[serde(default)]
    pub species: Vec<String>,
    #[serde(default)]
    pub gene_labels: Vec<String>,
    pub max_score: f64,
    #[serde(default)]
    pub peak_positions_promoter_relative_bp: Vec<i64>,
}

/// Simple promoter-sequence similarity row; kept distinct from TFBS evidence.
#[derive(Debug, Clone, Serialize, Deserialize, Default, PartialEq)]
#[serde(default)]
pub struct OrthologSequenceSimilarityRow {
    pub left_species: String,
    pub right_species: String,
    pub left_gene_label: String,
    pub right_gene_label: String,
    pub alignment_mode: String,
    pub compared_length_bp: usize,
    pub identical_bp: usize,
    pub identity_fraction: f64,
}

/// Cross-species CUT&RUN/occupancy assignment row.
#[derive(Debug, Clone, Serialize, Deserialize, Default, PartialEq, Eq)]
#[serde(default)]
pub struct OrthologCutRunSupportRow {
    pub species: String,
    pub gene_label: String,
    pub status: OrthologCutRunSupportStatus,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub nearest_peak_distance_bp: Option<i64>,
    #[serde(default)]
    pub contributing_dataset_ids: Vec<String>,
    #[serde(default)]
    pub contributing_read_report_ids: Vec<String>,
    pub detail: String,
}

/// Expression evidence row carried into an ortholog promoter comparison.
#[derive(Debug, Clone, Serialize, Deserialize, Default, PartialEq)]
#[serde(default)]
pub struct OrthologExpressionAssignment {
    pub species: String,
    pub gene_label: String,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub sample_id: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub condition: Option<String>,
    pub value: f64,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub unit: Option<String>,
    pub source: String,
    pub assignment_note: String,
}

/// Portable cross-species promoter comparison report.
#[derive(Debug, Clone, Serialize, Deserialize, Default, PartialEq)]
#[serde(default)]
pub struct OrthologPromoterComparisonReport {
    pub schema: String,
    pub generated_at_unix_ms: u128,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub op_id: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub run_id: Option<String>,
    pub cohort: OrthologPromoterCohortReport,
    #[serde(default)]
    pub motifs_requested: Vec<String>,
    pub score_kind: String,
    pub clip_negative: bool,
    #[serde(default)]
    pub promoter_summaries: Vec<OrthologTfbsSummaryRow>,
    #[serde(default)]
    pub pairwise_tfbs_similarity: Vec<OrthologPairwiseTfbsSimilarity>,
    #[serde(default)]
    pub conserved_tfbs_peaks: Vec<OrthologTfbsPeakSummary>,
    #[serde(default)]
    pub species_specific_tfbs_peaks: Vec<OrthologTfbsPeakSummary>,
    #[serde(default)]
    pub sequence_similarity: Vec<OrthologSequenceSimilarityRow>,
    #[serde(default)]
    pub cutrun_support: Vec<OrthologCutRunSupportRow>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub cutrun_quantitative_comparison: Option<OrthologCutRunQuantitativeComparison>,
    #[serde(default)]
    pub expression_assignments: Vec<OrthologExpressionAssignment>,
    #[serde(default)]
    pub relationship: GeneSetCohortRelationship,
    #[serde(default)]
    pub relationship_flags: Vec<GeneSetCohortRelationshipFlag>,
    #[serde(default)]
    pub warnings: Vec<String>,
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn ortholog_resource_defaults_to_offline_schema_shape() {
        let resource: OrthologResource = serde_json::from_value(serde_json::json!({
            "schema": ORTHOLOG_RESOURCE_SCHEMA,
            "species_aliases": [{"species": "Homo sapiens", "aliases": ["human"]}],
            "rows": [{
                "source_species": "Homo sapiens",
                "source_gene_symbol": "TP73",
                "target_species": "Mus musculus",
                "target_gene_symbol": "Trp73",
                "orthology_type": "one_to_one",
                "confidence": "high"
            }]
        }))
        .expect("deserialize resource");
        assert_eq!(resource.rows.len(), 1);
        assert_eq!(resource.rows[0].confidence.as_deref(), Some("high"));
        assert_eq!(
            resource.rows[0]
                .confidence
                .as_ref()
                .and_then(OrthologConfidence::level),
            Some(OrthologConfidenceLevel::High)
        );
        assert_eq!(
            resource.rows[0]
                .orthology_type
                .as_ref()
                .and_then(OrthologyType::cardinality),
            Some(OrthologyCardinality::OneToOne)
        );
        assert_eq!(resource.species_aliases[0].aliases, vec!["human"]);
    }

    #[test]
    fn ortholog_open_vocabularies_preserve_provider_specific_legacy_text() {
        let row: OrthologMappingRow = serde_json::from_value(serde_json::json!({
            "source_species": "Homo sapiens",
            "source_gene_symbol": "TP73",
            "target_species": "Mus musculus",
            "target_gene_symbol": "Trp73",
            "orthology_type": "synthetic_ortholog",
            "confidence": "tutorial"
        }))
        .expect("deserialize provider-specific values");

        assert_eq!(row.orthology_type.as_deref(), Some("synthetic_ortholog"));
        assert_eq!(row.confidence.as_deref(), Some("tutorial"));
        assert_eq!(
            row.orthology_type
                .as_ref()
                .and_then(OrthologyType::cardinality),
            None
        );
        assert_eq!(
            row.confidence.as_ref().and_then(OrthologConfidence::level),
            None
        );
        let serialized = serde_json::to_value(row).expect("serialize provider-specific values");
        assert_eq!(serialized["orthology_type"], "synthetic_ortholog");
        assert_eq!(serialized["confidence"], "tutorial");
    }

    #[test]
    fn ortholog_resource_context_references_are_validated() {
        let mut resource: OrthologResource = serde_json::from_value(serde_json::json!({
            "schema": ORTHOLOG_RESOURCE_SCHEMA,
            "contexts": [
                {
                    "context_id": "human",
                    "organism": "Homo sapiens",
                    "genome_id": "GRCh38"
                },
                {
                    "context_id": "mouse",
                    "organism": "Mus musculus",
                    "genome_id": "GRCm39"
                }
            ],
            "rows": [{
                "source_species": "Homo sapiens",
                "source_context_id": "human",
                "source_gene_symbol": "TP73",
                "target_species": "Mus musculus",
                "target_context_id": "mouse",
                "target_gene_symbol": "Trp73"
            }]
        }))
        .expect("deserialize context-bound resource");

        resource
            .validate_context_references()
            .expect("known context references");
        resource.rows[0].target_context_id = Some("unknown".to_string());
        assert!(matches!(
            resource.validate_context_references(),
            Err(BiologicalContextResolutionError::UnknownContextId {
                context_id
            }) if context_id == "unknown"
        ));
    }

    #[test]
    fn preserve_ambiguity_policy_and_candidate_mapping_round_trip() {
        let unresolved: OrthologUnresolvedRow = serde_json::from_value(serde_json::json!({
            "species": "Mus musculus",
            "reason": "ambiguous",
            "candidate_mappings": [{
                "candidate_rank": 1,
                "candidate_label": "Mus musculus:Trp73:ENSMUSG_TRP73:one_to_one",
                "target_species": "Mus musculus",
                "target_genome_id": "GRCm39",
                "source_context_id": "human",
                "target_context_id": "mouse",
                "target_gene_id": "ENSMUSG_TRP73",
                "target_gene_symbol": "Trp73",
                "orthology_type": "one_to_one",
                "confidence": "provider_reviewed",
                "source": "curated provider",
                "evidence": ["provider_assertion"]
            }]
        }))
        .expect("deserialize preserved ambiguity");
        assert_eq!(unresolved.candidate_mappings.len(), 1);
        assert_eq!(
            unresolved.candidate_mappings[0]
                .orthology_type
                .as_ref()
                .and_then(OrthologyType::cardinality),
            Some(OrthologyCardinality::OneToOne)
        );
        assert_eq!(
            unresolved.candidate_mappings[0].confidence.as_deref(),
            Some("provider_reviewed")
        );

        let policy = OrthologAmbiguityPolicy::Preserve;
        assert_eq!(policy.as_str(), "preserve");
        assert_eq!(
            serde_json::to_value(policy).expect("serialize policy"),
            "preserve"
        );
        assert_eq!(
            serde_json::from_value::<OrthologAmbiguityPolicy>(serde_json::json!("preserve"))
                .expect("deserialize policy"),
            OrthologAmbiguityPolicy::Preserve
        );
        assert!(
            serde_json::from_value::<OrthologAmbiguityPolicy>(serde_json::json!("future_policy"))
                .is_err(),
            "operation-control policies are closed and must not silently default"
        );
        let serialized = serde_json::to_value(unresolved).expect("serialize preserved ambiguity");
        assert_eq!(
            serialized["candidate_mappings"][0]["confidence"],
            "provider_reviewed"
        );

        let legacy: OrthologUnresolvedRow = serde_json::from_value(serde_json::json!({
            "species": "Mus musculus",
            "reason": "ambiguous",
            "candidates": ["Trp73", "Trp73b"]
        }))
        .expect("deserialize legacy unresolved row");
        assert!(legacy.candidate_mappings.is_empty());
        assert!(
            serde_json::to_value(legacy)
                .expect("serialize legacy unresolved row")
                .get("candidate_mappings")
                .is_none()
        );
    }

    #[test]
    fn old_promoter_cohort_defaults_new_optional_rows() {
        let report: OrthologPromoterCohortReport = serde_json::from_value(serde_json::json!({
            "schema": ORTHOLOG_PROMOTER_COHORT_SCHEMA,
            "generated_at_unix_ms": 1,
            "request": {
                "anchor_species": "Homo sapiens",
                "anchor_genome_id": "HumanToy",
                "anchor_gene_query": "TP73"
            }
        }))
        .expect("deserialize old-shaped cohort");
        assert_eq!(
            report.request.ambiguity_policy,
            OrthologAmbiguityPolicy::Reject
        );
        assert_eq!(
            report.request.relationship,
            GeneSetCohortRelationship::Unspecified
        );
        assert_eq!(report.relationship, GeneSetCohortRelationship::Unspecified);
        assert!(report.relationship_flags.is_empty());
        assert!(report.rows.is_empty());
        assert!(report.unresolved_rows.is_empty());
        assert!(report.warnings.is_empty());
        assert!(report.biological_contexts.contexts.is_empty());
    }

    #[test]
    fn cutrun_status_serializes_snake_case() {
        let row = OrthologCutRunSupportRow {
            species: "Homo sapiens".to_string(),
            gene_label: "TP73".to_string(),
            status: OrthologCutRunSupportStatus::MotifOnly,
            detail: "motif present; no occupancy data".to_string(),
            ..OrthologCutRunSupportRow::default()
        };
        let value = serde_json::to_value(row).expect("serialize row");
        assert_eq!(value["status"], "motif_only");
    }

    #[test]
    fn normalized_cutrun_contract_round_trips_without_changing_old_reports() {
        let old_report: OrthologPromoterComparisonReport =
            serde_json::from_value(serde_json::json!({
                "schema": ORTHOLOG_PROMOTER_COMPARISON_SCHEMA,
                "cohort": {"schema": ORTHOLOG_PROMOTER_COHORT_SCHEMA}
            }))
            .expect("deserialize comparison without normalized CUT&RUN");
        assert!(old_report.cutrun_quantitative_comparison.is_none());

        let comparison = OrthologCutRunQuantitativeComparison {
            status: OrthologCutRunQuantitativeComparisonStatus::Comparable,
            normalization: OrthologCutRunNormalizationInput {
                normalization_method: "spike_in_scaled_cpm".to_string(),
                unit: "normalized_fragments_per_million".to_string(),
                comparison_reference: "shared_batch_1".to_string(),
                provenance: "synthetic reviewed normalization table".to_string(),
                values: vec![OrthologCutRunNormalizedValueInput {
                    species: "Homo sapiens".to_string(),
                    gene_label: Some("TP73".to_string()),
                    normalized_value: 4.5,
                    contributing_read_report_ids: vec!["human_reads".to_string()],
                    ..OrthologCutRunNormalizedValueInput::default()
                }],
            },
            detail: "Explicit normalized values were comparable.".to_string(),
            ..OrthologCutRunQuantitativeComparison::default()
        };
        let value = serde_json::to_value(comparison).expect("serialize normalized comparison");
        assert_eq!(value["status"], "comparable");
        assert_eq!(
            value["normalization"]["normalization_method"],
            "spike_in_scaled_cpm"
        );
    }
}
