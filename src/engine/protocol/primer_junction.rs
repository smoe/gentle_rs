//! Read-only joined inspection of saved trace/query/payload evidence; no read or reference replay.

use super::*;

pub const PRIMER_JUNCTION_INSPECTION_REQUEST_SCHEMA: &str =
    "gentle.primer_junction_inspection_request.v1";
pub const PRIMER_JUNCTION_INSPECTION_SCHEMA: &str = "gentle.primer_junction_inspection.v1";

/// Conservative saved-evidence inspection defaults, independent of the producer's admission ceilings.
pub fn primer_junction_inspection_default_budgets() -> PrimerPayloadBudgets {
    PrimerPayloadBudgets {
        max_input_bytes: 64_000_000,
        max_output_bytes: 8_000_000,
        max_records: 4096,
        max_queries: 4096,
        max_query_bases: 2_000_000,
        max_line_bytes: 1_000_000,
        max_alignments_per_query: 128,
        max_cigar_operations: 4096,
        max_annotation_rows: 100_000,
        max_reference_context_bases: 1_000_000,
        max_chain_comparisons: 1_000_000,
    }
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct InspectPrimerPayloadJunctionsRequest {
    pub schema: String,
    pub source: PrimerPayloadResultBinding,
    pub source_report: PrimerPayloadFile,
    /// Row filters and pagination never filter the verified cohort/profile denominators.
    pub file_id: Option<String>,
    pub primer_id: Option<String>,
    #[serde(default)]
    pub offset: usize,
    #[serde(default = "primer_junction_default_rows")]
    pub rows: usize,
    #[serde(default = "primer_junction_inspection_default_budgets")]
    pub budgets: PrimerPayloadBudgets,
}
fn primer_junction_default_rows() -> usize {
    100
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PrimerPayloadJunctionRow {
    pub file_id: String,
    pub source_sha256: String,
    pub sample_id: String,
    pub barcode: String,
    pub priming: String,
    /// Original-read coordinates and all hit alternatives, not reoriented into a guessed isoform.
    pub trace: RnaPrimerTraceRow,
    /// Saved primer-oriented geometry. The representative anchor does not classify every alternative primer.
    pub queries: Vec<PrimerPayloadQuery>,
    /// First at most 64 already-oriented, already-masked saved query bases; keyed by query ID.
    pub payload_previews: BTreeMap<String, String>,
    pub evidence: Vec<PrimerPayloadReadEvidence>,
    pub unassessed: Vec<PrimerPayloadUnassessed>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PrimerPayloadJunctionInspection {
    pub schema: String,
    pub content_sha256: String,
    pub request_sha256: String,
    pub request: InspectPrimerPayloadJunctionsRequest,
    pub source_status: String,
    pub source_control_error: Option<String>,
    pub trace_status: String,
    pub trace_control_error: Option<RnaReadPanelControlError>,
    pub query_status: String,
    pub query_control_error: Option<String>,
    pub cohort: PrimerPayloadCohort,
    pub primer_panel: RnaPrimerPanel,
    pub trace_inputs_not_completed: Vec<RnaReadPanelInput>,
    pub query_files_not_completed: Vec<String>,
    pub batches_not_completed: Vec<String>,
    /// Selected, verified query-source units: one physical source record, including handle-only/excluded reads.
    pub cohort_distinct_reads: usize,
    pub matched_trace_rows: usize,
    pub rows_truncated: bool,
    pub rows: Vec<PrimerPayloadJunctionRow>,
    pub profiles: Vec<PrimerPayloadPartnerProfile>,
    pub breakpoint_clusters: Vec<PrimerPayloadBreakpointCluster>,
    pub saved_artifacts_verified: bool,
    pub original_reads_reverified: bool,
    pub reference_bases_reverified: bool,
    pub warnings: Vec<String>,
}
