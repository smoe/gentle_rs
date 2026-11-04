//! Prospective, unweighted RT-primer terminal compatibility against explicit FASTA inputs.

use super::{InlineSequenceTopology, PrimerPayloadFile, RnaLibraryPrimer};
use serde::{Deserialize, Serialize};

pub const RT_PRIMER_ASSESSMENT_REQUEST_SCHEMA: &str = "gentle.rt_primer_assessment_request.v1";
pub const RT_PRIMER_ASSESSMENT_SCHEMA: &str = "gentle.rt_primer_assessment.v1";

#[derive(Debug, Clone, Copy, Serialize, Deserialize, PartialEq, Eq)]
#[serde(rename_all = "snake_case")]
pub enum RtPrimerBackgroundScope {
    Genome,
    Transcriptome,
    MitochondrialRna,
    RibosomalRna,
    CustomRna,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct RtPrimerBackground {
    pub reference_id: String,
    pub scope: RtPrimerBackgroundScope,
    /// Assembly/release or explicit synthetic provenance. Labels do not verify identity.
    pub source: String,
    /// Raw, uncompressed FASTA bytes; RNA is represented using T, in biological 5′–3′ order.
    pub fasta: Option<PrimerPayloadFile>,
    /// Required exactly when fasta is absent; an absent input is never a negative screen.
    pub unavailable_reason: Option<String>,
    /// Topology is explicit per record; unspecified records are linear.
    #[serde(default)]
    pub circular_record_ids: Vec<String>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct RtPrimerIntendedSite {
    pub primer_id: String,
    pub reference_id: String,
    pub record_id: String,
    pub tail_start_0based: usize,
    pub tail_end_0based_exclusive: usize,
    /// Strand of the template bound by the primer, relative to the FASTA record.
    pub template_strand: char,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(default, deny_unknown_fields)]
pub struct RtPrimerAssessmentBudgets {
    pub max_input_bytes: u64,
    pub max_reference_bases: usize,
    pub max_record_bases: usize,
    pub max_records: usize,
    /// Counts individual base comparisons, including both genomic template orientations.
    pub max_base_comparisons: u64,
    /// Per matrix cell; counts remain exhaustive when the example list fills.
    pub max_examples_per_cell: usize,
}

impl Default for RtPrimerAssessmentBudgets {
    fn default() -> Self {
        Self {
            max_input_bytes: 16_000_000,
            max_reference_bases: 1_000_000,
            max_record_bases: 1_000_000,
            max_records: 10_000,
            max_base_comparisons: 100_000_000,
            max_examples_per_cell: 8,
        }
    }
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct AssessRtPrimerCandidatesRequest {
    pub schema: String,
    /// Existing complete oligos and tiled common_handle/terminal target_core components.
    pub primers: Vec<RnaLibraryPrimer>,
    pub references: Vec<RtPrimerBackground>,
    #[serde(default)]
    pub intended_sites: Vec<RtPrimerIntendedSite>,
    #[serde(default)]
    pub budgets: RtPrimerAssessmentBudgets,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct RtPrimerSeedCounts {
    /// Physical start positions, not orientation comparisons or molecules.
    pub assessed_windows: u64,
    pub ambiguous_windows: u64,
    /// Index is positional mismatches; best remaining template orientation per physical window.
    pub mismatch_histogram: Vec<u64>,
    /// Index is the exact contiguous run from the primer's 3′ end; longest remaining orientation.
    pub terminal_run_histogram: Vec<u64>,
    pub exact_matching_records: usize,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct RtPrimerSeedExample {
    pub record_id: String,
    pub start_0based: usize,
    /// May exceed record length for a circular-origin window.
    pub end_0based_exclusive: usize,
    pub topology: InlineSequenceTopology,
    /// All exact template orientations for this physical site, sorted '+' then '-'.
    pub template_strands: Vec<char>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct RtPrimerSeedAssessment {
    pub primer_id: String,
    pub reference_id: String,
    pub seed_length: usize,
    pub status: String,
    pub reason: Option<String>,
    pub records_shorter_than_seed: usize,
    pub all_sites: RtPrimerSeedCounts,
    pub excluding_intended_orientations: RtPrimerSeedCounts,
    pub exact_examples: Vec<RtPrimerSeedExample>,
    pub exact_examples_truncated: bool,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct RtPrimerBackgroundAssessment {
    pub reference_id: String,
    pub status: String,
    pub reason: Option<String>,
    pub records: usize,
    pub bases: usize,
    pub template_strands: Vec<char>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct RtPrimerAssessmentReport {
    pub schema: String,
    pub content_sha256: String,
    pub request_sha256: String,
    pub request: AssessRtPrimerCandidatesRequest,
    pub references: Vec<RtPrimerBackgroundAssessment>,
    pub cells: Vec<RtPrimerSeedAssessment>,
    pub base_comparisons: u64,
    pub abundance_status: String,
    pub retrospective_status: String,
    pub warnings: Vec<String>,
}
