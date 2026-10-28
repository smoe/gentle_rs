//! Pure observed-junction kernel; no file access, project mutation or isoform/protein assertion.
//!
//! Callers bind and verify input bytes. This S0 API validates in-memory geometry,
//! retains alternatives, and compares explicit stable gene IDs in existing payload
//! annotations. S1's engine boundary verifies files and publishes a separate report.

use crate::digest_utils::{canonical_json_string, sha256_prefixed_str};
use crate::engine::{
    EngineError, PrimerPayloadAnnotation, PrimerPayloadReferenceScope, PrimerPayloadSegment,
};
use serde::{Deserialize, Serialize};
use std::collections::{BTreeMap, BTreeSet};

pub(crate) mod annotation;
mod comparison;
mod core;
pub(crate) mod dag;
pub(crate) mod geometry;
mod kernel;
// S2.0b is deliberately private until a separately reviewed ingress contract exists.
#[allow(dead_code)]
mod sam;
// S2.2a decodes supplied raw lines only; verified file/unit ingress stays separate.
#[allow(dead_code)]
mod sam_text;
#[cfg(test)]
mod tests;

pub use geometry::Block as SupportedBlock;
pub use kernel::summarize;
pub use kernel::summarize_with_control;

pub type Interval = (usize, usize);

/// Coordinate authority, not an annotation/provider label or a selected gene.
#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct ReferenceBinding {
    pub assembly: String,
    pub fasta_sha256: String,
}

/// A gene ID in one explicitly identified annotation; no symbol guessing.
#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct ReferenceGene {
    pub source_id: String,
    pub gene_id: String,
}

/// Reuses Phase 2's annotation geometry without introducing a third transcript model.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct AnnotationSource {
    pub source_id: String,
    pub provider: Option<String>,
    pub annotation_sha256: String,
    pub annotation: PrimerPayloadAnnotation,
}

/// Original record identity: aliases and repeated QNAMEs cannot create extra support.
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq, PartialOrd, Ord)]
pub struct ReadRecord {
    pub source_sha256: String,
    pub record_index: u64,
    pub read_id: String,
}

#[derive(Debug, Clone, Copy, Serialize, Deserialize, PartialEq, Eq, PartialOrd, Ord)]
#[serde(rename_all = "snake_case")]
pub enum EndConstraint {
    PrimerBounded,
    Masked,
    ReadEndUnverified,
    Clipped,
    Unknown,
}

/// One query/anchor view of a record. Missing masks do not mean an empty mask.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct QueryEvidence {
    pub record: ReadRecord,
    pub alignment_sha256: String,
    pub query_id: String,
    pub query_length: usize,
    pub scope: PrimerPayloadReferenceScope,
    pub target_lengths: BTreeMap<String, usize>,
    pub masks: Option<Vec<Interval>>,
    pub ends: [Vec<EndConstraint>; 2],
    pub search_complete: bool,
    pub segments: Vec<PrimerPayloadSegment>,
}

/// Resource ceilings, not truncation policies or calibrated biological defaults.
#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Budgets {
    pub max_queries: usize,
    pub max_records: usize,
    pub max_segments_per_query: usize,
    pub max_operations: usize,
    pub max_skips_per_chain: usize,
    pub max_chain_families: usize,
    pub max_alternatives_per_record: usize,
    pub max_annotation_rows: usize,
    pub max_model_comparisons_per_family: usize,
    pub max_total_work: usize,
}

impl Budgets {
    /// Shared admission before either in-memory interpretation or file/control effects.
    pub(crate) fn validate(&self, flank: usize) -> Result<(), EngineError> {
        if !(1..=1_000_000).contains(&flank)
            || [
                self.max_queries,
                self.max_records,
                self.max_segments_per_query,
                self.max_operations,
                self.max_skips_per_chain,
                self.max_chain_families,
                self.max_alternatives_per_record,
                self.max_annotation_rows,
                self.max_model_comparisons_per_family,
            ]
            .iter()
            .any(|n| !(1..=1_000_000).contains(n))
            || !(1..=100_000_000).contains(&self.max_total_work)
        {
            return Err(bad("observed_structure_admission"));
        }
        Ok(())
    }
}

/// Pure S0 input. The selector is mandatory; other providers need explicit ID bindings.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Input {
    pub reference: ReferenceBinding,
    pub reference_gene: ReferenceGene,
    pub comparison_genes: Vec<ReferenceGene>,
    pub annotation_sources: Vec<AnnotationSource>,
    pub queries: Vec<QueryEvidence>,
    /// Excluded/missing views must not strengthen support from another anchor of the same record.
    #[serde(default)]
    pub unassessed_records: Vec<ReadRecord>,
    pub required_flank_bases: usize,
    pub budgets: Budgets,
}

#[derive(Debug, Clone, Copy, Serialize, Deserialize, PartialEq, Eq)]
#[serde(rename_all = "snake_case")]
pub enum GeneResolution {
    Resolved,
    MissingGene,
    AmbiguousLocus,
    UnmappedSource,
}

#[derive(Debug, Clone, Copy, Serialize, Deserialize, PartialEq, Eq)]
#[serde(rename_all = "snake_case")]
pub enum DetailStatus {
    Detailed,
    CoarseUnassessed,
    TranscriptomeUnassessed,
}

#[derive(Debug, Clone, Copy, Serialize, Deserialize, PartialEq, Eq, PartialOrd, Ord)]
#[serde(rename_all = "snake_case")]
pub enum Topology {
    Spliced,
    Unspliced,
    Segmented,
    Unassessed,
}

#[derive(Debug, Clone, Copy, Serialize, Deserialize, PartialEq, Eq)]
#[serde(rename_all = "snake_case")]
pub enum Qualification {
    Qualified,
    InsufficientSupport,
    Unassessed,
}

#[derive(Debug, Clone, Copy, Serialize, Deserialize, PartialEq, Eq)]
#[serde(rename_all = "snake_case")]
pub enum ModelRelation {
    ExactFullChain,
    ExactPartialChain,
    TerminalCoverageDivergent,
    IncompatibleCoverage,
    WithinExon,
    ExonIntronBoundary,
    Intronic,
    OutsideModel,
}

#[derive(Debug, Clone, Copy, Serialize, Deserialize, PartialEq, Eq)]
#[serde(rename_all = "snake_case")]
pub enum JunctionRelation {
    ContainsUnannotatedJunction,
    AllJunctionsAnnotatedNoSingleModelContainsRun,
    ConsecutiveModelRun,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct TranscriptMatch {
    pub transcript_id: String,
    pub gene_ids: Vec<String>,
    /// Annotation evidence, not an observed RNA-strand assignment.
    pub model_strand: char,
    pub relation: ModelRelation,
    pub retained_intron_supported_bases: Option<usize>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
/// Exact per-source model membership; absent comparison is represented by the enclosing Option.
pub struct JunctionMatch {
    pub skipped_interval: Interval,
    pub annotated_in_models: Vec<String>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SourceComparison {
    pub source_id: String,
    pub gene_id: Option<String>,
    pub resolution: GeneResolution,
    pub junction_relation: Option<JunctionRelation>,
    /// None means unavailable comparison, not an unannotated/absent junction.
    pub junctions: Option<Vec<JunctionMatch>>,
    pub models: Vec<TranscriptMatch>,
}

/// One segment observation; continuous alignment skips are never invented between segments.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Observation {
    pub observation_id: String,
    pub record: ReadRecord,
    pub query_id: String,
    pub contig: String,
    pub observed_span: Interval,
    pub aligned_query_interval: Interval,
    pub payload_query_length: usize,
    pub alignment_clipped: [bool; 2],
    pub query_to_reference_strand: char,
    pub reported_mapq: Option<u8>,
    /// Constraints on payload ends, not a claim that the alignment reaches them.
    pub ends: [Vec<EndConstraint>; 2],
    pub detail: DetailStatus,
    pub topology: Topology,
    pub blocks: Vec<SupportedBlock>,
    pub excluded_query_intervals: Option<Vec<Interval>>,
    pub skips: Vec<Interval>,
    pub flank_qualification: Vec<Qualification>,
    pub observed_junction_chain: Option<String>,
    pub comparisons: Vec<SourceComparison>,
}

#[derive(Debug, Clone, Copy, Serialize, Deserialize, PartialEq, Eq)]
#[serde(rename_all = "snake_case")]
pub enum Reconciliation {
    Agreement,
    ContainedPartial,
    Alternatives,
    Unassessed,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct RecordSupport {
    pub record: ReadRecord,
    /// Each path contains observation IDs, not an unqualified union of DAG nodes.
    pub hypotheses: Vec<Vec<String>>,
    pub reconciliation: Reconciliation,
    pub topologies: Vec<Topology>,
    pub search_declared_complete: bool,
    pub all_query_views_assessed: bool,
    pub definite_chains: Vec<String>,
    pub possible_chains: Vec<String>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ChainFamily {
    pub chain_id: String,
    pub contig: String,
    pub skips: Vec<Interval>,
    pub supporting_records: usize,
    pub definite_records: usize,
}

/// Source receipts without duplicating the annotation's transcript geometry.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct AnnotationBinding {
    pub source_id: String,
    pub provider: Option<String>,
    pub annotation_sha256: String,
    pub reference_id: String,
    pub assembly: String,
    pub release: String,
    pub reference_fasta_sha256: String,
}

/// Rust-level result; S1 v1 wraps this serialized shape in a versioned artifact.
/// Stored-v1 semantics require an explicit compatibility/version decision before changing.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Summary {
    pub context_sha256: String,
    pub reference: ReferenceBinding,
    pub annotation_bindings: Vec<AnnotationBinding>,
    pub reference_gene: ReferenceGene,
    pub gene_resolutions: Vec<SourceComparison>,
    pub observations: Vec<Observation>,
    pub records: Vec<RecordSupport>,
    pub families: Vec<ChainFamily>,
    pub comparison_work: usize,
}

fn bad(message: &str) -> EngineError {
    EngineError::invalid_input(message)
}

fn digest<T: Serialize>(value: &T) -> Result<String, EngineError> {
    Ok(sha256_prefixed_str(
        &canonical_json_string(value).map_err(|_| bad("observed_structure_digest"))?,
    ))
}

fn is_hash(value: &str) -> bool {
    value.strip_prefix("sha256:").is_some_and(|s| {
        s.len() == 64
            && s.bytes()
                .all(|x| x.is_ascii_digit() || (b'a'..=b'f').contains(&x))
    })
}

fn union(mut values: Vec<Interval>) -> Vec<Interval> {
    values.sort_unstable();
    let mut result: Vec<Interval> = Vec::new();
    for (a, z) in values {
        if let Some(last) = result.last_mut()
            && a <= last.1
        {
            last.1 = last.1.max(z);
        } else {
            result.push((a, z));
        }
    }
    result
}

fn unmasked(range: Interval, masks: &[Interval]) -> Vec<Interval> {
    let mut start = range.0;
    let mut result = Vec::new();
    for &(a, z) in &masks[masks.partition_point(|&(_, z)| z <= range.0)..] {
        if a >= range.1 {
            break;
        }
        if z <= start {
            continue;
        }
        if a > start {
            result.push((start, a));
        }
        start = start.max(z);
    }
    if start < range.1 {
        result.push((start, range.1));
    }
    result
}
