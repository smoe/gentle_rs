//! Shared deterministic operation engine and project-state contract.
//!
//! This module is the single execution layer for GENtle biology/workflow
//! behavior. GUI, CLI, JavaScript, and Lua adapters are expected to translate
//! user intent into operations defined here instead of re-implementing domain
//! logic in frontend code.
//!
//! Core responsibilities:
//! - Define persistent state contracts (`ProjectState`, lineage, containers,
//!   display settings, and metadata schemas).
//! - Define typed operation/workflow payloads and structured results/warnings.
//! - Execute cloning/editing/analysis/render-export operations through one
//!   deterministic path.
//! - Maintain provenance and operation identifiers used by replay, audit, and
//!   adapter-equivalent automation.
//!
//! Invariants:
//! - Behavioral parity across adapters: same operation input yields equivalent
//!   state/result output regardless of entry point.
//! - Operation semantics are explicit and machine-readable; hidden adapter state
//!   must not change core biological outcomes.
//! - Long-running operations use cooperative progress callbacks and explicit
//!   cancellation signaling where supported.

use crate::{
    DNA_LADDERS, RNA_LADDERS,
    app::GENtleApp,
    dna_sequence::DNAsequence,
    enzymes::active_restriction_enzymes,
    feature_location::{collect_location_ranges_usize, feature_is_reverse},
    genomes::{
        BlastExternalBinaryPreflightReport, DEFAULT_GENOME_CATALOG_PATH,
        DEFAULT_HELPER_GENOME_CATALOG_PATH, EnsemblCatalogUpdatePreview,
        EnsemblCatalogUpdateReport, GenomeBlastReport, GenomeCatalog,
        GenomeCatalogEntryRemovalReport, GenomeGeneRecord, GenomeSourcePlan,
        GenomeTranscriptRecord, PrepareGenomePlan, PrepareGenomeProgress, PrepareGenomeReport,
        PreparedCacheCleanupReport, PreparedCacheCleanupRequest, PreparedCacheInspectionReport,
        PreparedGenomeFallbackPolicy, PreparedGenomeInspection, PreparedGenomeRemovalReport,
        blast_external_binary_preflight_report, build_genbank_efetch_url,
        clear_prepared_cache_roots, inspect_prepared_cache_roots, is_prepare_cancelled_error,
        validate_genbank_accession,
    },
    iupac_code::IupacCode,
    lineage_export::export_lineage_svg,
    methylation_sites::MethylationMode,
    pool_gel::{GelSampleInput, export_pool_gel_svg},
    protocol_cartoon::{ProtocolCartoonKind, ProtocolCartoonTemplateBindings},
    render_export::{export_circular_svg, export_linear_svg},
    render_feature_expert::render_feature_expert_svg,
    restriction_enzyme::{RestrictionEnzyme, RestrictionEnzymeKey},
    rna_structure::{self, RnaStructureError, RnaStructureSvgReport, RnaStructureTextReport},
    tf_motifs,
    uniprot::{
        UniprotAaGenomicSegment, UniprotEntry, UniprotEntrySummary, UniprotFeatureProjection,
        UniprotGenomeProjection, UniprotGenomeProjectionSummary, UniprotTranscriptProjection,
        normalize_entry_id as normalize_uniprot_entry_id, parse_swiss_prot_text,
    },
};
use flate2::read::MultiGzDecoder;
use rayon::join;
use regex::{Regex, RegexBuilder};
use serde::{Deserialize, Serialize};
use serde_json::json;
use std::{
    cell::Cell,
    cmp::Ordering,
    collections::hash_map::DefaultHasher,
    collections::{BTreeMap, BTreeSet, BinaryHeap, HashMap, HashSet},
    env,
    fs::{File, OpenOptions},
    hash::{Hash, Hasher},
    io::{BufRead, BufReader, BufWriter, Read, Write},
    path::{Path, PathBuf},
    process::{Command, Stdio},
    sync::{
        Arc,
        atomic::{AtomicU64, Ordering as AtomicOrdering},
    },
    time::{Duration, Instant},
};
use tempfile::NamedTempFile;

pub use gentle_protocol::{ContainerId, NodeId, OpId, RunId, SeqId};

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum PrepareReferenceGenomeMode {
    PrepareOrReuse,
    ReindexCachedFiles,
    RefreshFromSources,
}
pub use crate::feature_expert::{
    FeatureExpertTarget, FeatureExpertView, ISOFORM_ARCHITECTURE_EXPERT_INSTRUCTION,
    IsoformArchitectureCdsAaSegment, IsoformArchitectureExpertView,
    IsoformArchitectureProteinDomain, IsoformArchitectureProteinLane,
    IsoformArchitectureTranscriptLane, RESTRICTION_EXPERT_INSTRUCTION, RestrictionSiteExpertView,
    SPLICING_EXPERT_INSTRUCTION, SplicingBoundaryMarker, SplicingEventSummary,
    SplicingExonCdsPhase, SplicingExonSummary, SplicingExpertView, SplicingJunctionArc,
    SplicingMatrixRow, SplicingRange, SplicingScopePreset, SplicingTranscriptLane,
    TFBS_EXPERT_INSTRUCTION, TfbsExpertColumn, TfbsExpertView,
};
const PROVENANCE_METADATA_KEY: &str = "provenance";
const GENOME_EXTRACTIONS_METADATA_KEY: &str = "genome_extractions";
pub const GENOME_TRACK_SUBSCRIPTIONS_METADATA_KEY: &str = "genome_bed_track_subscriptions";
const GENOME_TRACK_KNOWN_ANCHORS_METADATA_KEY: &str = "genome_track_autosync_known_anchors";
pub const CANDIDATE_SETS_METADATA_KEY: &str = "candidate_sets";
const CANDIDATE_SETS_SCHEMA: &str = "gentle.candidate_sets.v1";
const CANDIDATE_SETS_REF_SCHEMA: &str = "gentle.candidate_sets.ref.v1";
const CANDIDATE_SETS_DISK_INDEX_SCHEMA: &str = "gentle.candidate_sets.disk_index.v1";
const CANDIDATE_SETS_LOAD_WARNING_METADATA_KEY: &str = "candidate_sets_load_warning";
const CANDIDATE_STORE_STRICT_LOAD_ENV: &str = "GENTLE_CANDIDATE_STORE_STRICT_LOAD";
pub const GUIDE_DESIGN_METADATA_KEY: &str = "guide_design";
const GUIDE_DESIGN_SCHEMA: &str = "gentle.guide_design.v1";
pub const PRIMER_DESIGN_REPORTS_METADATA_KEY: &str = "primer_design_reports";
const PRIMER_DESIGN_REPORTS_SCHEMA: &str = "gentle.primer_design_reports.v1";
const PRIMER_DESIGN_REPORT_SCHEMA: &str = "gentle.primer_design_report.v1";
const QPCR_DESIGN_REPORT_SCHEMA: &str = "gentle.qpcr_design_report.v1";
pub const SEQUENCING_TRACES_METADATA_KEY: &str = "sequencing_traces";
const SEQUENCING_TRACES_SCHEMA: &str = "gentle.sequencing_traces.v1";
pub const SEQUENCING_TRACE_RECORD_SCHEMA: &str = "gentle.sequencing_trace_record.v1";
pub const SEQUENCING_TRACE_IMPORT_REPORT_SCHEMA: &str = "gentle.sequencing_trace_import_report.v1";
pub const SEQUENCING_CONFIRMATION_REPORTS_METADATA_KEY: &str = "sequencing_confirmation_reports";
const SEQUENCING_CONFIRMATION_REPORTS_SCHEMA: &str = "gentle.sequencing_confirmation_reports.v1";
pub const SEQUENCING_CONFIRMATION_REPORT_SCHEMA: &str = "gentle.sequencing_confirmation_report.v1";
pub const SEQUENCING_CONFIRMATION_SUPPORT_TSV_SCHEMA: &str =
    "gentle.sequencing_confirmation_support_tsv.v1";
pub const PLANNING_METADATA_KEY: &str = "planning";
pub const PLANNING_PROFILE_SCHEMA: &str = "gentle.planning_profile.v1";
pub const PLANNING_OBJECTIVE_SCHEMA: &str = "gentle.planning_objective.v1";
pub const PLANNING_ESTIMATE_SCHEMA: &str = "gentle.planning_estimate.v1";
pub const PLANNING_SUGGESTION_SCHEMA: &str = "gentle.planning_suggestion.v1";
pub const PLANNING_SYNC_STATUS_SCHEMA: &str = "gentle.planning_sync_status.v1";
const PLANNING_STORE_SCHEMA: &str = "gentle.planning_store.v1";
pub const WORKFLOW_MACRO_TEMPLATES_METADATA_KEY: &str = "workflow_macro_templates";
const WORKFLOW_MACRO_TEMPLATES_SCHEMA: &str = "gentle.workflow_macro_templates.v1";
pub const CLONING_MACRO_TEMPLATE_SCHEMA: &str = "gentle.cloning_macro_template.v1";
pub const CANDIDATE_MACRO_TEMPLATES_METADATA_KEY: &str = "candidate_macro_templates";
const CANDIDATE_MACRO_TEMPLATES_SCHEMA: &str = "gentle.candidate_macro_templates.v1";
const GENOME_BED_TRACK_GENERATED_TAG: &str = "genome_bed_track";
const GENOME_BIGWIG_TRACK_GENERATED_TAG: &str = "genome_bigwig_track";
const GENOME_VCF_TRACK_GENERATED_TAG: &str = "genome_vcf_track";
const DBSNP_VARIANT_MARKER_GENERATED_TAG: &str = "dbsnp_variant_marker";
const BLAST_HIT_TRACK_GENERATED_TAG: &str = "blast_hit_track";
const GENOME_ANNOTATION_PROJECTION_GENERATED_TAG: &str = "genome_annotation_projection";
const CONTEXT_LAYER_QUALIFIER_KEY: &str = "gentle_context_layer";
const CONTEXT_LAYER_TRANSCRIPT: &str = "contextual_transcript";
const CONTEXT_LAYER_GENE: &str = "contextual_gene";
const BLAST_OPTIONS_OVERRIDE_METADATA_KEY: &str = "blast_options_override";
const BLAST_OPTIONS_DEFAULTS_PATH_METADATA_KEY: &str = "blast_options_defaults_path";
const DEFAULT_BLAST_OPTIONS_PATH: &str = "assets/blast_defaults.json";
const ISOFORM_PANELS_METADATA_KEY: &str = "isoform_panels";
const ISOFORM_PANELS_SCHEMA: &str = "gentle.isoform_panels.v1";
const ISOFORM_PANEL_RESOURCE_SCHEMA: &str = "gentle.isoform_panel_resource.v1";
const ISOFORM_PANEL_VALIDATION_REPORT_SCHEMA: &str = "gentle.isoform_panel_validation_report.v1";
const UNIPROT_ENTRIES_METADATA_KEY: &str = "uniprot_entries";
const UNIPROT_ENTRIES_SCHEMA: &str = "gentle.uniprot_entries.v1";
const UNIPROT_GENOME_PROJECTIONS_METADATA_KEY: &str = "uniprot_genome_projections";
const UNIPROT_GENOME_PROJECTIONS_SCHEMA: &str = "gentle.uniprot_genome_projections.v1";
const UNIPROT_GENOME_PROJECTION_SCHEMA: &str = "gentle.uniprot_genome_projection.v1";
const PROCESS_RUN_BUNDLE_SCHEMA: &str = "gentle.process_run_bundle.v1";
pub const ROUTINE_DECISION_TRACES_METADATA_KEY: &str = "routine_decision_traces";
pub const ROUTINE_DECISION_TRACE_SCHEMA: &str = "gentle.routine_decision_trace.v1";
pub const ROUTINE_DECISION_TRACE_STORE_SCHEMA: &str = "gentle.routine_decision_trace_store.v1";
pub const DOTPLOT_ANALYSIS_METADATA_KEY: &str = "dotplot_analysis";
const DOTPLOT_ANALYSIS_SCHEMA: &str = "gentle.dotplot_analysis_store.v1";
const DOTPLOT_VIEW_SCHEMA: &str = "gentle.dotplot_view.v3";
const FLEXIBILITY_TRACK_SCHEMA: &str = "gentle.flexibility_track.v1";
const SEQUENCE_ALIGNMENT_REPORT_SCHEMA: &str = "gentle.sequence_alignment_report.v1";
pub const RNA_READ_REPORTS_METADATA_KEY: &str = "rna_read_reports";
const RNA_READ_REPORTS_SCHEMA: &str = "gentle.rna_read_reports.v1";
const RNA_READ_REPORT_SCHEMA: &str = "gentle.rna_read_report.v1";
const RNA_READ_CHECKPOINT_SCHEMA: &str = "gentle.rna_read_interpret_checkpoint.v1";
const RNA_READ_SAMPLE_SHEET_EXPORT_SCHEMA: &str = "gentle.rna_read_sample_sheet_export.v1";
const RNA_READ_EXON_PATHS_EXPORT_SCHEMA: &str = "gentle.rna_read_exon_paths_export.v1";
const RNA_READ_EXON_ABUNDANCE_EXPORT_SCHEMA: &str = "gentle.rna_read_exon_abundance_export.v1";
const RNA_READ_GENE_SUPPORT_SUMMARY_SCHEMA: &str = "gentle.rna_read_gene_support_summary.v1";
const RNA_READ_SCORE_DENSITY_SVG_EXPORT_SCHEMA: &str =
    "gentle.rna_read_score_density_svg_export.v1";
const RNA_READ_ALIGNMENT_TSV_EXPORT_SCHEMA: &str = "gentle.rna_read_alignment_tsv_export.v1";
const RNA_READ_ALIGNMENT_DOTPLOT_SVG_EXPORT_SCHEMA: &str =
    "gentle.rna_read_alignment_dotplot_svg_export.v1";
const RNA_READ_ALIGNMENT_INSPECTION_SCHEMA: &str = "gentle.rna_read_alignment_inspection.v1";
const RNA_READ_ALIGNMENT_DETAIL_SCHEMA: &str = "gentle.rna_read_alignment_detail.v1";
#[cfg(debug_assertions)]
const RNA_READ_PROGRESS_UPDATE_EVERY_READS: usize = 1_000;
#[cfg(not(debug_assertions))]
const RNA_READ_PROGRESS_UPDATE_EVERY_READS: usize = 10_000;
const RNA_READ_ALIGNMENT_PROGRESS_UPDATE_EVERY_READS: usize = 1;
const RNA_READ_PROGRESS_UPDATE_MAX_INTERVAL: Duration = Duration::from_secs(2);
const RNA_READ_PROGRESS_MAX_HISTOGRAM_BINS: usize = 200;
const RNA_READ_SCORE_DENSITY_BIN_COUNT: usize = 40;
const RNA_READ_RETAINED_HITS_MAX: usize = 5_000;
const RNA_READ_RETAINED_TOP_SCORE_GUARANTEE_COUNT: usize = 2_000;
const RNA_READ_RETAINED_HIGH_SCORE_BIN_GUARANTEE_COUNT: usize = 20;
const RNA_READ_CHECKPOINT_DEFAULT_EVERY_READS: usize = 10_000;
const RNA_READ_COOPERATIVE_YIELD_EVERY_READS: usize = 512;
const RNA_READ_PROGRESS_TOP_HITS_PREVIEW_MAX: usize = 20;
const RNA_READ_SEED_CHAIN_MAX_CANDIDATES_PER_BIT: usize = 64;
const RNA_READ_INFER_PARALLEL_MIN_MATCHED_BITS: usize = 96;
const RNA_READ_INFER_PARALLEL_MIN_MATCHED_OBSERVATIONS: usize = 256;
const MAX_DOTPLOT_POINTS: usize = 250_000;
const DOTPLOT_BOXPLOT_DEFAULT_BINS: usize = 96;
pub const MAX_DOTPLOT_PAIR_EVALUATIONS: usize = 100_000_000;
pub const DEFAULT_BIGWIG_TO_BEDGRAPH_BIN: &str = "bigWigToBedGraph";
pub const BIGWIG_TO_BEDGRAPH_ENV_BIN: &str = "GENTLE_BIGWIG_TO_BEDGRAPH_BIN";
const MAX_IMPORTED_SIGNAL_FEATURES: usize = 25_000;
const DEFAULT_EXTRACT_REGION_ANNOTATION_FEATURE_CAP: usize = 5_000;
const EXON_CONCAT_SPACER_BP: usize = 24;
const HELPER_MCS_GENERATED_TAG: &str = "helper_mcs";
const PUC19_MCS_SEQUENCE: &str = "GAATTCGAGCTCGGTACCCGGGGATCCTCTAGAGTCGACCTGCAGGCATGCAAGCTT";
const PUC18_MCS_SEQUENCE: &str = "AAGCTTGCATGCCTGCAGGTCGACTCTAGAGGATCCCCGGGTACCGAGCTCGAATTC";
const PUC_MCS_EXPECTED_SITES: &str = "EcoRI,SacI,KpnI,SmaI,BamHI,XbaI,SalI,PstI,SphI,HindIII";
const PRIMER_PREFERRED_MIN_LENGTH_BP: usize = 20;
const PRIMER_PREFERRED_MAX_LENGTH_BP: usize = 30;
const PRIMER_RECOMMENDED_MAX_HOMOPOLYMER_RUN_BP: usize = 4;
const PRIMER_RECOMMENDED_MAX_SELF_COMPLEMENTARY_RUN_BP: usize = 6;
const PRIMER_RECOMMENDED_MAX_PAIR_DIMER_RUN_BP: usize = 6;
const PRIMER_RECOMMENDED_MAX_PAIR_3PRIME_DIMER_RUN_BP: usize = 3;
const PRIMER_INTERNAL_MAX_PAIR_EVALUATIONS: usize = 1_000_000;
const FEATURE_QUERY_RESULT_SCHEMA: &str = "gentle.sequence_feature_query_result.v1";
const FEATURE_QUERY_DEFAULT_LIMIT: usize = 200;
const FEATURE_QUERY_MAX_LIMIT: usize = 10_000;

// Private decomposition slices of the engine implementation. Shared public
// contracts stay in this file; heavy helpers and operation families live in the
// corresponding `src/engine/*` module so future edits can land in one focused
// area without changing adapter-visible APIs.
#[path = "engine/ops/candidate_guides.rs"]
mod candidate_guides;
#[path = "engine/analysis/candidate_metrics.rs"]
mod candidate_metrics;
#[path = "engine/state/feature_coordinate_formulas.rs"]
mod feature_coordinate_formulas;
#[path = "engine/analysis/feature_expert_ops.rs"]
mod feature_expert_ops;
#[path = "engine/io/genome_tracks.rs"]
mod genome_tracks;
#[path = "engine/io/import_anchors.rs"]
mod import_anchors;
#[path = "engine/state/lineage_containers.rs"]
mod lineage_containers;
#[path = "engine/ops/operation_handlers.rs"]
mod operation_handlers;
#[path = "engine/analysis/rna_reads.rs"]
mod rna_reads;
#[path = "engine/state/sequence_ops.rs"]
mod sequence_ops;
#[path = "engine/analysis/sequencing_confirmation.rs"]
mod sequencing_confirmation;
#[path = "engine/io/sequencing_traces.rs"]
mod sequencing_traces;

#[path = "engine/protocol.rs"]
mod protocol;
pub use feature_coordinate_formulas::{
    parse_feature_coordinate_term_on_sequence, parse_required_usize_or_formula_text_on_sequence,
    resolve_formula_roi_range_inputs_0based_on_sequence,
    resolve_selection_formula_range_0based_on_sequence, split_feature_formula_range_expression,
};
pub use protocol::*;

// Shared default helpers used by both engine operations and the extracted
// public protocol layer.
fn default_rna_align_report_selection() -> RnaReadHitSelection {
    RnaReadHitSelection::All
}

fn default_true() -> bool {
    true
}

fn default_rna_read_checkpoint_every_reads() -> usize {
    RNA_READ_CHECKPOINT_DEFAULT_EVERY_READS
}

fn default_rna_read_alignment_dotplot_max_points() -> usize {
    2_500
}

fn default_splicing_reference_scope() -> SplicingScopePreset {
    SplicingScopePreset::TargetGroupTargetStrand
}

fn default_pairwise_match_score() -> i32 {
    2
}

fn default_pairwise_mismatch_score() -> i32 {
    -3
}

fn default_pairwise_gap_open() -> i32 {
    -5
}

fn default_pairwise_gap_extend() -> i32 {
    -1
}

fn default_sequencing_confirmation_alignment_mode() -> PairwiseAlignmentMode {
    PairwiseAlignmentMode::Local
}

fn default_sequencing_confirmation_min_identity_fraction() -> f64 {
    0.80
}

fn default_sequencing_confirmation_min_target_coverage_fraction() -> f64 {
    1.0
}

// Private RNA-read execution state remains in the engine implementation.
// Only the stable adapter-facing report/progress contracts live in
// `engine/protocol.rs`.
#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
struct RnaReadReportStore {
    schema: String,
    updated_at_unix_ms: u128,
    reports: BTreeMap<String, RnaReadInterpretationReport>,
}

#[derive(Debug, Clone, Default)]
struct ParsedFastaReadRecord {
    record_index: usize,
    source_byte_offset: usize,
    header_id: String,
    sequence: Vec<u8>,
}

#[derive(Debug, Clone, Copy, Default)]
struct FastaVisitProgress {
    records_processed: usize,
    input_bytes_processed: u64,
    input_bytes_total: u64,
    io_read_ms: f64,
    record_parse_ms: f64,
}

#[derive(Debug)]
struct CountingReader<R> {
    inner: R,
    bytes_read: Arc<AtomicU64>,
}

impl<R> CountingReader<R> {
    fn new(inner: R, bytes_read: Arc<AtomicU64>) -> Self {
        Self { inner, bytes_read }
    }
}

impl<R: Read> Read for CountingReader<R> {
    fn read(&mut self, buf: &mut [u8]) -> std::io::Result<usize> {
        let read = self.inner.read(buf)?;
        self.bytes_read
            .fetch_add(read as u64, AtomicOrdering::Relaxed);
        Ok(read)
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
struct RnaReadRetentionRank {
    passed_seed_filter: bool,
    has_alignment: bool,
    alignment_identity_ppm: u32,
    alignment_query_coverage_ppm: u32,
    alignment_score: i64,
    weighted_support_milli: u64,
    weighted_seed_hit_ppm: u32,
    seed_hit_ppm: u32,
    matched_kmers: usize,
    tested_kmers: usize,
    read_length_bp: usize,
    record_index: usize,
}

impl Ord for RnaReadRetentionRank {
    fn cmp(&self, other: &Self) -> Ordering {
        self.passed_seed_filter
            .cmp(&other.passed_seed_filter)
            .then(self.has_alignment.cmp(&other.has_alignment))
            .then(
                self.alignment_query_coverage_ppm
                    .cmp(&other.alignment_query_coverage_ppm),
            )
            .then(
                self.alignment_identity_ppm
                    .cmp(&other.alignment_identity_ppm),
            )
            .then(self.alignment_score.cmp(&other.alignment_score))
            .then(
                self.weighted_support_milli
                    .cmp(&other.weighted_support_milli),
            )
            .then(self.weighted_seed_hit_ppm.cmp(&other.weighted_seed_hit_ppm))
            .then(self.seed_hit_ppm.cmp(&other.seed_hit_ppm))
            .then(self.matched_kmers.cmp(&other.matched_kmers))
            .then(self.tested_kmers.cmp(&other.tested_kmers))
            .then(self.read_length_bp.cmp(&other.read_length_bp))
            .then(other.record_index.cmp(&self.record_index))
    }
}

impl PartialOrd for RnaReadRetentionRank {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
struct RnaReadPhase1ScoreRank {
    seed_hit_ppm: u32,
    weighted_seed_hit_ppm: u32,
    weighted_support_milli: u64,
    matched_kmers: usize,
    tested_kmers: usize,
    read_length_bp: usize,
    record_index: usize,
}

impl Ord for RnaReadPhase1ScoreRank {
    fn cmp(&self, other: &Self) -> Ordering {
        self.seed_hit_ppm
            .cmp(&other.seed_hit_ppm)
            .then(self.weighted_seed_hit_ppm.cmp(&other.weighted_seed_hit_ppm))
            .then(
                self.weighted_support_milli
                    .cmp(&other.weighted_support_milli),
            )
            .then(self.matched_kmers.cmp(&other.matched_kmers))
            .then(self.tested_kmers.cmp(&other.tested_kmers))
            .then(self.read_length_bp.cmp(&other.read_length_bp))
            .then(other.record_index.cmp(&self.record_index))
    }
}

impl PartialOrd for RnaReadPhase1ScoreRank {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

#[derive(Debug, Clone)]
struct RetainedRnaReadHit {
    rank: RnaReadRetentionRank,
    hit: RnaReadInterpretationHit,
}

impl Ord for RetainedRnaReadHit {
    fn cmp(&self, other: &Self) -> Ordering {
        // Keep the weakest retained hit at heap top for O(log N) replacement.
        other.rank.cmp(&self.rank)
    }
}

impl PartialOrd for RetainedRnaReadHit {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl PartialEq for RetainedRnaReadHit {
    fn eq(&self, other: &Self) -> bool {
        self.rank == other.rank
    }
}

impl Eq for RetainedRnaReadHit {}

#[derive(Debug, Clone)]
struct RetainedRnaReadScoreHit {
    rank: RnaReadPhase1ScoreRank,
    hit: RnaReadInterpretationHit,
}

impl Ord for RetainedRnaReadScoreHit {
    fn cmp(&self, other: &Self) -> Ordering {
        // Keep the weakest guaranteed-score hit at heap top for O(log N)
        // replacement, mirroring the main retained-hit heap behavior.
        other.rank.cmp(&self.rank)
    }
}

impl PartialOrd for RetainedRnaReadScoreHit {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl PartialEq for RetainedRnaReadScoreHit {
    fn eq(&self, other: &Self) -> bool {
        self.rank == other.rank
    }
}

impl Eq for RetainedRnaReadScoreHit {}

#[derive(Debug, Clone)]
struct RetainedRnaReadPreviewHit {
    rank: RnaReadRetentionRank,
    preview: RnaReadTopHitPreview,
}

impl Ord for RetainedRnaReadPreviewHit {
    fn cmp(&self, other: &Self) -> Ordering {
        other.rank.cmp(&self.rank)
    }
}

impl PartialOrd for RetainedRnaReadPreviewHit {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl PartialEq for RetainedRnaReadPreviewHit {
    fn eq(&self, other: &Self) -> bool {
        self.rank == other.rank
    }
}

impl Eq for RetainedRnaReadPreviewHit {}

#[derive(Debug, Clone, Default)]
struct SplicingTranscriptTemplate {
    transcript_feature_id: usize,
    transcript_id: String,
    transcript_label: String,
    strand: String,
    sequence: Vec<u8>,
    genomic_positions_1based: Vec<usize>,
    kmer_positions: HashMap<u32, Vec<usize>>,
}

#[derive(Debug, Clone, Default)]
struct TranscriptExonPathModel {
    transcript_feature_id: usize,
    transcript_id: String,
    transcript_label: String,
    strand: String,
    exon_ordinals: Vec<usize>,
    transitions: Vec<(usize, usize)>,
}

#[derive(Debug, Clone, Copy, Default)]
struct SeedHistogramWeight {
    bin_index: usize,
    strand_minus: bool,
}

#[derive(Debug, Clone, Copy, Default)]
struct SeedTemplatePosition {
    template_idx: usize,
    template_pos: usize,
}

#[derive(Debug, Clone, Copy, Default)]
struct SeedMatchObservation {
    read_start: usize,
    bits: u32,
}

#[derive(Debug, Clone, Default)]
struct SeedChainSpacingMetrics {
    transcript_id: String,
    support_kmers: usize,
    support_fraction: f64,
    median_transcript_gap: f64,
    transcript_gap_count: usize,
}

#[derive(Debug, Clone, Default)]
struct RnaReadOriginClassification {
    origin_class: RnaReadOriginClass,
    reason: String,
    origin_confidence: f64,
    strand_confidence: f64,
}

#[derive(Debug, Clone, Default)]
struct ReadExonPathInference {
    path: String,
    confirmed_transitions: usize,
    total_transitions: usize,
    transcript_id: String,
    strand: String,
    strand_diagnostics: RnaReadStrandAssignmentDiagnostics,
}

#[derive(Debug, Clone, Copy)]
struct TranscriptSupportScore<'a> {
    model: &'a TranscriptExonPathModel,
    exon_hits: usize,
    transition_hits: usize,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
struct IsoformSupportAccumulator {
    transcript_feature_id: usize,
    transcript_id: String,
    transcript_label: String,
    strand: String,
    exon_count: usize,
    expected_transition_count: usize,
    reads_assigned: usize,
    reads_seed_passed: usize,
    transition_rows_supported: HashSet<(usize, usize)>,
    seed_gap_sum: f64,
    seed_gap_count: usize,
    confirmed_transition_fraction_sum: f64,
    confirmed_transition_fraction_count: usize,
    best_seed_hit_fraction: f64,
    best_weighted_seed_hit_fraction: f64,
    reads_chain_same_strand: usize,
    reads_with_opposite_strand_competition: usize,
    reads_ambiguous_strand_ties: usize,
}

#[derive(Debug, Clone)]
struct RnaReadInterpretOptions {
    report_mode: RnaReadReportMode,
    checkpoint_path: Option<String>,
    checkpoint_every_reads: usize,
    resume_from_checkpoint: bool,
}

impl Default for RnaReadInterpretOptions {
    fn default() -> Self {
        Self {
            report_mode: RnaReadReportMode::Full,
            checkpoint_path: None,
            checkpoint_every_reads: default_rna_read_checkpoint_every_reads(),
            resume_from_checkpoint: false,
        }
    }
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
struct RnaReadInterpretCheckpoint {
    schema: String,
    created_at_unix_ms: u128,
    report_id: String,
    report_mode: RnaReadReportMode,
    seq_id: String,
    seed_feature_id: usize,
    profile: RnaReadInterpretationProfile,
    input_path: String,
    input_format: RnaReadInputFormat,
    scope: SplicingScopePreset,
    origin_mode: RnaReadOriginMode,
    target_gene_ids: Vec<String>,
    roi_seed_capture_enabled: bool,
    seed_filter: RnaReadSeedFilterConfig,
    align_config: RnaReadAlignConfig,
    checkpoint_path: Option<String>,
    checkpoint_every_reads: usize,
    reads_processed: usize,
    read_count_seed_passed: usize,
    read_count_aligned: usize,
    input_bytes_processed: u64,
    input_bytes_total: u64,
    cumulative_tested_kmers: usize,
    cumulative_matched_kmers: usize,
    cumulative_seed_compute_ms: f64,
    cumulative_align_compute_ms: f64,
    cumulative_io_read_ms: f64,
    cumulative_fasta_parse_ms: f64,
    cumulative_normalize_compute_ms: f64,
    cumulative_inference_compute_ms: f64,
    cumulative_progress_emit_ms: f64,
    cumulative_read_bases_processed: u64,
    read_length_counts_all: Vec<u64>,
    read_length_counts_seed_passed: Vec<u64>,
    read_length_counts_aligned: Vec<u64>,
    read_length_counts_full_length_exact: Vec<u64>,
    read_length_counts_full_length_near: Vec<u64>,
    read_length_counts_full_length_strict: Vec<u64>,
    support_aligned_reads: usize,
    support_exon_counts: Vec<usize>,
    support_junction_counts: Vec<usize>,
    reads_with_transition_support: usize,
    transition_confirmations: usize,
    transition_support_rows: Vec<RnaReadTransitionSupportRow>,
    isoform_support_accumulators: BTreeMap<String, IsoformSupportAccumulator>,
    origin_class_counts: BTreeMap<String, usize>,
    bins: Vec<RnaReadSeedHistogramBin>,
    score_density_bins: Vec<u64>,
    seed_pass_score_density_bins: Vec<u64>,
    retained_hits: Vec<RnaReadInterpretationHit>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
/// High-level provenance class for how a sequence entered/was derived in the
/// project graph.
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
/// One sequence node in lineage DAG state.
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
/// Stored in `ProjectState.lineage.macro_instances` for audit/debug and graph
/// visualization.
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
/// Lineage DAG plus macro-instance records.
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

#[derive(Debug, Clone, Serialize, Deserialize)]
/// Named arrangement definition referencing container lanes.
pub struct Arrangement {
    pub arrangement_id: String,
    pub mode: ArrangementMode,
    pub name: Option<String>,
    pub lane_container_ids: Vec<ContainerId>,
    #[serde(default)]
    pub ladders: Vec<String>,
    pub created_by_op: Option<OpId>,
    pub created_at_unix_ms: u128,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
/// Container/arrangement state persisted with the project.
pub struct ContainerState {
    pub containers: HashMap<ContainerId, Container>,
    pub arrangements: HashMap<String, Arrangement>,
    pub seq_to_latest_container: HashMap<SeqId, ContainerId>,
    pub next_container_counter: u64,
    pub next_arrangement_counter: u64,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(default)]
/// Engine-level behavioral parameters (non-rendering).
pub struct EngineParameters {
    pub max_fragments_per_container: usize,
    pub require_verified_genome_anchor_for_extension: bool,
    pub genome_anchor_prepared_fallback_policy: GenomeAnchorPreparedFallbackPolicy,
    pub primer_design_backend: PrimerDesignBackend,
    pub primer3_executable: String,
}

impl Default for EngineParameters {
    fn default() -> Self {
        Self {
            max_fragments_per_container: 80_000,
            require_verified_genome_anchor_for_extension: false,
            genome_anchor_prepared_fallback_policy:
                GenomeAnchorPreparedFallbackPolicy::SingleCompatible,
            primer_design_backend: PrimerDesignBackend::Auto,
            primer3_executable: "primer3_core".to_string(),
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
/// Prepared-genome fallback policy when an anchor references a non-prepared id.
pub enum GenomeAnchorPreparedFallbackPolicy {
    Off,
    SingleCompatible,
    AlwaysExplicit,
}

impl GenomeAnchorPreparedFallbackPolicy {
    pub fn as_str(self) -> &'static str {
        match self {
            Self::Off => "off",
            Self::SingleCompatible => "single_compatible",
            Self::AlwaysExplicit => "always_explicit",
        }
    }
}

#[derive(Debug, Clone, Default, Serialize, Deserialize)]
/// Complete persisted project state consumed by `GentleEngine`.
///
/// `metadata` stores versioned extension payloads; critical shared keys are
/// defined as constants in this module.
pub struct ProjectState {
    pub sequences: HashMap<SeqId, DNAsequence>,
    pub metadata: HashMap<String, serde_json::Value>,
    #[serde(default)]
    pub display: DisplaySettings,
    #[serde(default)]
    pub lineage: LineageGraph,
    #[serde(default)]
    pub parameters: EngineParameters,
    #[serde(default)]
    pub container_state: ContainerState,
}

#[derive(Debug)]
enum CandidateSidecarTransaction {
    Noop,
    Replace {
        final_dir: PathBuf,
        staging_dir: PathBuf,
    },
    Remove {
        final_dir: PathBuf,
    },
}

#[derive(Debug)]
enum CandidateSidecarCommitted {
    Noop,
    Replaced {
        final_dir: PathBuf,
        backup_dir: Option<PathBuf>,
    },
    Removed {
        final_dir: PathBuf,
        backup_dir: Option<PathBuf>,
    },
}

impl ProjectState {
    /// Load project state JSON from disk.
    ///
    /// This also hydrates candidate-set sidecar references when present.
    pub fn load_from_path(path: &str) -> Result<Self, EngineError> {
        let text = std::fs::read_to_string(path).map_err(|e| EngineError {
            code: ErrorCode::Io,
            message: format!("Could not read state file '{path}': {e}"),
        })?;
        let mut state: Self = serde_json::from_str(&text).map_err(|e| EngineError {
            code: ErrorCode::InvalidInput,
            message: format!("Could not parse state JSON '{path}': {e}"),
        })?;
        if let Err(err) = state.hydrate_candidate_store_from_external_ref(path) {
            if Self::strict_candidate_store_load_enabled() {
                return Err(err);
            }
            let warning = err.message;
            state.metadata.remove(CANDIDATE_SETS_METADATA_KEY);
            state.metadata.insert(
                CANDIDATE_SETS_LOAD_WARNING_METADATA_KEY.to_string(),
                json!({
                    "message": warning,
                    "strict_env": CANDIDATE_STORE_STRICT_LOAD_ENV,
                }),
            );
        } else {
            state
                .metadata
                .remove(CANDIDATE_SETS_LOAD_WARNING_METADATA_KEY);
        }
        Ok(state)
    }

    /// Save project state JSON to disk using atomic replacement semantics.
    ///
    /// Candidate-set sidecar data is staged/committed in lockstep to avoid
    /// partially-written project metadata.
    pub fn save_to_path(&self, path: &str) -> Result<(), EngineError> {
        let mut state_for_disk = self.clone();
        let project_path = Path::new(path);
        let sidecar_tx =
            state_for_disk.prepare_candidate_store_sidecar_transaction(project_path)?;
        let text = serde_json::to_string_pretty(&state_for_disk).map_err(|e| EngineError {
            code: ErrorCode::Internal,
            message: format!("Could not serialize state: {e}"),
        })?;
        let committed = Self::commit_candidate_store_transaction(sidecar_tx)?;
        if let Err(write_err) = Self::write_text_file_atomically(project_path, &text) {
            if let Err(rollback_err) = Self::rollback_candidate_store_transaction(committed) {
                return Err(EngineError {
                    code: ErrorCode::Io,
                    message: format!(
                        "{}; candidate-sidecar rollback also failed: {}",
                        write_err.message, rollback_err.message
                    ),
                });
            }
            return Err(write_err);
        }
        Self::finalize_candidate_store_transaction(committed)
    }

    fn candidate_store_sidecar_basename(project_path: &Path) -> String {
        let fallback = "project.gentle.json".to_string();
        let base = project_path
            .file_name()
            .and_then(|name| name.to_str())
            .unwrap_or(&fallback);
        format!("{base}.candidates")
    }

    fn candidate_store_sidecar_dir(project_path: &Path) -> PathBuf {
        let parent = project_path
            .parent()
            .map(Path::to_path_buf)
            .unwrap_or_else(|| PathBuf::from("."));
        parent.join(Self::candidate_store_sidecar_basename(project_path))
    }

    fn candidate_store_sidecar_index_rel(project_path: &Path) -> PathBuf {
        PathBuf::from(Self::candidate_store_sidecar_basename(project_path)).join("index.json")
    }

    fn strict_candidate_store_load_enabled() -> bool {
        matches!(
            env::var(CANDIDATE_STORE_STRICT_LOAD_ENV)
                .unwrap_or_default()
                .trim()
                .to_ascii_lowercase()
                .as_str(),
            "1" | "true" | "yes" | "on"
        )
    }

    fn write_text_file_atomically(path: &Path, text: &str) -> Result<(), EngineError> {
        let parent = path
            .parent()
            .map(Path::to_path_buf)
            .unwrap_or_else(|| PathBuf::from("."));
        std::fs::create_dir_all(&parent).map_err(|e| EngineError {
            code: ErrorCode::Io,
            message: format!(
                "Could not create parent directory '{}' for state save: {e}",
                parent.display()
            ),
        })?;
        let mut tmp = NamedTempFile::new_in(&parent).map_err(|e| EngineError {
            code: ErrorCode::Io,
            message: format!(
                "Could not create temporary state file in '{}': {e}",
                parent.display()
            ),
        })?;
        tmp.write_all(text.as_bytes()).map_err(|e| EngineError {
            code: ErrorCode::Io,
            message: format!(
                "Could not write temporary state file for '{}': {e}",
                path.display()
            ),
        })?;
        tmp.flush().map_err(|e| EngineError {
            code: ErrorCode::Io,
            message: format!(
                "Could not flush temporary state file for '{}': {e}",
                path.display()
            ),
        })?;
        tmp.as_file().sync_all().map_err(|e| EngineError {
            code: ErrorCode::Io,
            message: format!(
                "Could not sync temporary state file for '{}': {e}",
                path.display()
            ),
        })?;
        tmp.persist(path).map_err(|e| EngineError {
            code: ErrorCode::Io,
            message: format!(
                "Could not replace state file '{}': {}",
                path.display(),
                e.error
            ),
        })?;
        Ok(())
    }

    fn candidate_sidecar_staging_dir(final_dir: &Path) -> Result<PathBuf, EngineError> {
        let parent = final_dir
            .parent()
            .map(Path::to_path_buf)
            .unwrap_or_else(|| PathBuf::from("."));
        let base = final_dir
            .file_name()
            .and_then(|name| name.to_str())
            .unwrap_or("candidate_sidecar");
        let nonce = format!("{}-{}", std::process::id(), GentleEngine::now_unix_ms());
        for attempt in 0..64 {
            let suffix = if attempt == 0 {
                nonce.clone()
            } else {
                format!("{nonce}-{attempt}")
            };
            let candidate = parent.join(format!(".{base}.staging.{suffix}"));
            match std::fs::create_dir(&candidate) {
                Ok(_) => return Ok(candidate),
                Err(e) if e.kind() == std::io::ErrorKind::AlreadyExists => continue,
                Err(e) => {
                    return Err(EngineError {
                        code: ErrorCode::Io,
                        message: format!(
                            "Could not create candidate-sidecar staging directory '{}': {e}",
                            candidate.display()
                        ),
                    });
                }
            }
        }
        Err(EngineError {
            code: ErrorCode::Io,
            message: format!(
                "Could not create unique candidate-sidecar staging directory under '{}'",
                parent.display()
            ),
        })
    }

    fn candidate_sidecar_backup_dir(final_dir: &Path) -> Result<PathBuf, EngineError> {
        let parent = final_dir
            .parent()
            .map(Path::to_path_buf)
            .unwrap_or_else(|| PathBuf::from("."));
        let base = final_dir
            .file_name()
            .and_then(|name| name.to_str())
            .unwrap_or("candidate_sidecar");
        let nonce = format!("{}-{}", std::process::id(), GentleEngine::now_unix_ms());
        for attempt in 0..64 {
            let suffix = if attempt == 0 {
                nonce.clone()
            } else {
                format!("{nonce}-{attempt}")
            };
            let candidate = parent.join(format!(".{base}.backup.{suffix}"));
            if !candidate.exists() {
                return Ok(candidate);
            }
        }
        Err(EngineError {
            code: ErrorCode::Io,
            message: format!(
                "Could not allocate candidate-sidecar backup path under '{}'",
                parent.display()
            ),
        })
    }

    fn sanitize_candidate_set_file_stem(name: &str) -> String {
        let mut out = String::with_capacity(name.len());
        for ch in name.chars() {
            if ch.is_ascii_alphanumeric() {
                out.push(ch.to_ascii_lowercase());
            } else if ch.is_whitespace() || matches!(ch, '-' | '_' | '.') {
                out.push('_');
            }
        }
        let trimmed = out.trim_matches('_');
        if trimmed.is_empty() {
            "set".to_string()
        } else {
            trimmed.chars().take(48).collect()
        }
    }

    fn parse_inline_candidate_store_value(value: &serde_json::Value) -> Option<CandidateStore> {
        let schema = value
            .get("schema")
            .and_then(|v| v.as_str())
            .unwrap_or_default();
        let has_sets = value.get("sets").map(|v| v.is_object()).unwrap_or(false);
        if schema == CANDIDATE_SETS_REF_SCHEMA {
            return None;
        }
        if !has_sets && schema != CANDIDATE_SETS_SCHEMA {
            return None;
        }
        let mut store: CandidateStore = serde_json::from_value(value.clone()).ok()?;
        if store.schema.trim().is_empty() {
            store.schema = CANDIDATE_SETS_SCHEMA.to_string();
        }
        Some(store)
    }

    fn candidate_set_metrics_for_disk(set: &CandidateSet) -> Vec<String> {
        let mut names = BTreeSet::new();
        for candidate in &set.candidates {
            for metric in candidate.metrics.keys() {
                names.insert(metric.to_string());
            }
        }
        names.into_iter().collect()
    }

    fn load_candidate_store_from_ref(
        project_path: &Path,
        reference: &CandidateStoreReference,
    ) -> Result<CandidateStore, EngineError> {
        let index_path = if Path::new(&reference.index_path).is_absolute() {
            PathBuf::from(&reference.index_path)
        } else {
            project_path
                .parent()
                .map(Path::to_path_buf)
                .unwrap_or_else(|| PathBuf::from("."))
                .join(&reference.index_path)
        };
        let index_text = std::fs::read_to_string(&index_path).map_err(|e| EngineError {
            code: ErrorCode::Io,
            message: format!(
                "Could not read candidate-store index '{}': {e}",
                index_path.display()
            ),
        })?;
        let index: CandidateStoreDiskIndex =
            serde_json::from_str(&index_text).map_err(|e| EngineError {
                code: ErrorCode::InvalidInput,
                message: format!(
                    "Could not parse candidate-store index '{}': {e}",
                    index_path.display()
                ),
            })?;
        if !index.schema.trim().is_empty() && index.schema != CANDIDATE_SETS_DISK_INDEX_SCHEMA {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: format!(
                    "Unsupported candidate-store index schema '{}' in '{}'",
                    index.schema,
                    index_path.display()
                ),
            });
        }
        let index_dir = index_path
            .parent()
            .map(Path::to_path_buf)
            .unwrap_or_else(|| PathBuf::from("."));
        let mut sets: HashMap<String, CandidateSet> = HashMap::new();
        for entry in index.sets {
            let records_path = if Path::new(&entry.records_path).is_absolute() {
                PathBuf::from(&entry.records_path)
            } else {
                index_dir.join(&entry.records_path)
            };
            let file = File::open(&records_path).map_err(|e| EngineError {
                code: ErrorCode::Io,
                message: format!(
                    "Could not read candidate records '{}' for set '{}': {e}",
                    records_path.display(),
                    entry.name
                ),
            })?;
            let reader = BufReader::new(file);
            let mut candidates: Vec<CandidateRecord> = vec![];
            for (line_no, line_result) in reader.lines().enumerate() {
                let line = line_result.map_err(|e| EngineError {
                    code: ErrorCode::Io,
                    message: format!(
                        "Could not read candidate record line {} from '{}': {e}",
                        line_no + 1,
                        records_path.display()
                    ),
                })?;
                let trimmed = line.trim();
                if trimmed.is_empty() {
                    continue;
                }
                let candidate =
                    serde_json::from_str::<CandidateRecord>(trimmed).map_err(|e| EngineError {
                        code: ErrorCode::InvalidInput,
                        message: format!(
                            "Invalid candidate record JSON at '{}':{}: {}",
                            records_path.display(),
                            line_no + 1,
                            e
                        ),
                    })?;
                candidates.push(candidate);
            }
            sets.insert(
                entry.name.clone(),
                CandidateSet {
                    name: entry.name,
                    created_at_unix_ms: entry.created_at_unix_ms,
                    source_seq_ids: entry.source_seq_ids,
                    candidates,
                },
            );
        }
        Ok(CandidateStore {
            schema: CANDIDATE_SETS_SCHEMA.to_string(),
            updated_at_unix_ms: index.updated_at_unix_ms,
            sets,
        })
    }

    fn hydrate_candidate_store_from_external_ref(
        &mut self,
        project_path: &str,
    ) -> Result<(), EngineError> {
        let Some(value) = self.metadata.get(CANDIDATE_SETS_METADATA_KEY).cloned() else {
            return Ok(());
        };
        if Self::parse_inline_candidate_store_value(&value).is_some() {
            return Ok(());
        }
        let reference: CandidateStoreReference =
            serde_json::from_value(value.clone()).map_err(|e| EngineError {
                code: ErrorCode::InvalidInput,
                message: format!(
                    "Could not parse candidate-store reference metadata '{}': {e}",
                    CANDIDATE_SETS_METADATA_KEY
                ),
            })?;
        if reference.schema != CANDIDATE_SETS_REF_SCHEMA {
            return Ok(());
        }
        let store = Self::load_candidate_store_from_ref(Path::new(project_path), &reference)?;
        let store_value = serde_json::to_value(store).map_err(|e| EngineError {
            code: ErrorCode::Internal,
            message: format!("Could not hydrate candidate-store metadata: {e}"),
        })?;
        self.metadata
            .insert(CANDIDATE_SETS_METADATA_KEY.to_string(), store_value);
        Ok(())
    }

    fn load_candidate_store_for_externalization(
        &self,
        project_path: &Path,
        value: serde_json::Value,
    ) -> Result<Option<CandidateStore>, EngineError> {
        if let Some(store) = Self::parse_inline_candidate_store_value(&value) {
            return Ok(Some(store));
        }
        let reference: CandidateStoreReference =
            serde_json::from_value(value).map_err(|e| EngineError {
                code: ErrorCode::InvalidInput,
                message: format!(
                    "Could not parse candidate-store reference metadata '{}': {e}",
                    CANDIDATE_SETS_METADATA_KEY
                ),
            })?;
        if reference.schema != CANDIDATE_SETS_REF_SCHEMA {
            return Ok(None);
        }
        let store = Self::load_candidate_store_from_ref(project_path, &reference)?;
        Ok(Some(store))
    }

    fn prepare_candidate_store_sidecar_transaction(
        &mut self,
        project_path: &Path,
    ) -> Result<CandidateSidecarTransaction, EngineError> {
        let sidecar_dir = Self::candidate_store_sidecar_dir(project_path);
        let Some(value) = self.metadata.get(CANDIDATE_SETS_METADATA_KEY).cloned() else {
            if sidecar_dir.exists() {
                return Ok(CandidateSidecarTransaction::Remove {
                    final_dir: sidecar_dir,
                });
            }
            return Ok(CandidateSidecarTransaction::Noop);
        };
        let Some(mut store) = self.load_candidate_store_for_externalization(project_path, value)?
        else {
            return Ok(CandidateSidecarTransaction::Noop);
        };
        if store.sets.is_empty() {
            self.metadata.remove(CANDIDATE_SETS_METADATA_KEY);
            if sidecar_dir.exists() {
                return Ok(CandidateSidecarTransaction::Remove {
                    final_dir: sidecar_dir,
                });
            }
            return Ok(CandidateSidecarTransaction::Noop);
        }
        if store.schema.trim().is_empty() {
            store.schema = CANDIDATE_SETS_SCHEMA.to_string();
        }
        let staging_dir = Self::candidate_sidecar_staging_dir(&sidecar_dir)?;

        let mut set_names: Vec<String> = store.sets.keys().cloned().collect();
        set_names.sort_unstable();
        let mut index_entries: Vec<CandidateStoreDiskSetIndexEntry> = vec![];
        for (idx, set_name) in set_names.iter().enumerate() {
            let Some(set) = store.sets.get(set_name) else {
                continue;
            };
            let mut hasher = DefaultHasher::new();
            set_name.hash(&mut hasher);
            let stem = Self::sanitize_candidate_set_file_stem(set_name);
            let records_filename = format!("{:03}_{}_{}.jsonl", idx + 1, stem, hasher.finish());
            let records_path = staging_dir.join(&records_filename);
            let records_file = File::create(&records_path).map_err(|e| EngineError {
                code: ErrorCode::Io,
                message: format!(
                    "Could not create candidate-set records file '{}': {e}",
                    records_path.display()
                ),
            })?;
            let mut writer = BufWriter::new(records_file);
            for candidate in &set.candidates {
                serde_json::to_writer(&mut writer, candidate).map_err(|e| EngineError {
                    code: ErrorCode::Internal,
                    message: format!(
                        "Could not serialize candidate record for set '{}': {e}",
                        set_name
                    ),
                })?;
                writer.write_all(b"\n").map_err(|e| EngineError {
                    code: ErrorCode::Io,
                    message: format!(
                        "Could not write candidate record for set '{}': {e}",
                        set_name
                    ),
                })?;
            }
            writer.flush().map_err(|e| EngineError {
                code: ErrorCode::Io,
                message: format!(
                    "Could not flush candidate-set records file '{}': {e}",
                    records_path.display()
                ),
            })?;

            index_entries.push(CandidateStoreDiskSetIndexEntry {
                name: set.name.clone(),
                created_at_unix_ms: set.created_at_unix_ms,
                source_seq_ids: set.source_seq_ids.clone(),
                candidate_count: set.candidates.len(),
                metrics: Self::candidate_set_metrics_for_disk(set),
                records_path: records_filename,
            });
        }

        let index = CandidateStoreDiskIndex {
            schema: CANDIDATE_SETS_DISK_INDEX_SCHEMA.to_string(),
            updated_at_unix_ms: store.updated_at_unix_ms,
            set_count: index_entries.len(),
            sets: index_entries,
        };
        let index_path = staging_dir.join("index.json");
        let index_text = serde_json::to_string_pretty(&index).map_err(|e| EngineError {
            code: ErrorCode::Internal,
            message: format!("Could not serialize candidate-store index: {e}"),
        })?;
        std::fs::write(&index_path, index_text).map_err(|e| EngineError {
            code: ErrorCode::Io,
            message: format!(
                "Could not write candidate-store index '{}': {e}",
                index_path.display()
            ),
        })?;

        let reference = CandidateStoreReference {
            schema: CANDIDATE_SETS_REF_SCHEMA.to_string(),
            storage: "jsonl_indexed".to_string(),
            index_path: Self::candidate_store_sidecar_index_rel(project_path)
                .to_string_lossy()
                .to_string(),
            set_count: index.set_count,
            updated_at_unix_ms: index.updated_at_unix_ms,
        };
        let reference_value = serde_json::to_value(reference).map_err(|e| EngineError {
            code: ErrorCode::Internal,
            message: format!("Could not serialize candidate-store reference metadata: {e}"),
        })?;
        self.metadata
            .insert(CANDIDATE_SETS_METADATA_KEY.to_string(), reference_value);
        Ok(CandidateSidecarTransaction::Replace {
            final_dir: sidecar_dir,
            staging_dir,
        })
    }

    fn commit_candidate_store_transaction(
        tx: CandidateSidecarTransaction,
    ) -> Result<CandidateSidecarCommitted, EngineError> {
        match tx {
            CandidateSidecarTransaction::Noop => Ok(CandidateSidecarCommitted::Noop),
            CandidateSidecarTransaction::Replace {
                final_dir,
                staging_dir,
            } => {
                let mut backup_dir: Option<PathBuf> = None;
                if final_dir.exists() {
                    let backup = Self::candidate_sidecar_backup_dir(&final_dir)?;
                    std::fs::rename(&final_dir, &backup).map_err(|e| EngineError {
                        code: ErrorCode::Io,
                        message: format!(
                            "Could not move existing candidate-sidecar directory '{}' to backup '{}': {e}",
                            final_dir.display(),
                            backup.display()
                        ),
                    })?;
                    backup_dir = Some(backup);
                }
                if let Err(e) = std::fs::rename(&staging_dir, &final_dir) {
                    if let Some(backup) = &backup_dir {
                        let _ = std::fs::rename(backup, &final_dir);
                    }
                    let _ = std::fs::remove_dir_all(&staging_dir);
                    return Err(EngineError {
                        code: ErrorCode::Io,
                        message: format!(
                            "Could not activate candidate-sidecar directory '{}' from staging '{}': {e}",
                            final_dir.display(),
                            staging_dir.display()
                        ),
                    });
                }
                Ok(CandidateSidecarCommitted::Replaced {
                    final_dir,
                    backup_dir,
                })
            }
            CandidateSidecarTransaction::Remove { final_dir } => {
                if !final_dir.exists() {
                    return Ok(CandidateSidecarCommitted::Noop);
                }
                let backup = Self::candidate_sidecar_backup_dir(&final_dir)?;
                std::fs::rename(&final_dir, &backup).map_err(|e| EngineError {
                    code: ErrorCode::Io,
                    message: format!(
                        "Could not stage candidate-sidecar removal '{}' -> '{}': {e}",
                        final_dir.display(),
                        backup.display()
                    ),
                })?;
                Ok(CandidateSidecarCommitted::Removed {
                    final_dir,
                    backup_dir: Some(backup),
                })
            }
        }
    }

    fn rollback_candidate_store_transaction(
        committed: CandidateSidecarCommitted,
    ) -> Result<(), EngineError> {
        match committed {
            CandidateSidecarCommitted::Noop => Ok(()),
            CandidateSidecarCommitted::Replaced {
                final_dir,
                backup_dir,
            } => {
                if let Some(backup) = backup_dir {
                    if final_dir.exists() {
                        std::fs::remove_dir_all(&final_dir).map_err(|e| EngineError {
                            code: ErrorCode::Io,
                            message: format!(
                                "Could not remove partially committed candidate-sidecar directory '{}': {e}",
                                final_dir.display()
                            ),
                        })?;
                    }
                    std::fs::rename(&backup, &final_dir).map_err(|e| EngineError {
                        code: ErrorCode::Io,
                        message: format!(
                            "Could not restore candidate-sidecar backup '{}' -> '{}': {e}",
                            backup.display(),
                            final_dir.display()
                        ),
                    })?;
                } else if final_dir.exists() {
                    std::fs::remove_dir_all(&final_dir).map_err(|e| EngineError {
                        code: ErrorCode::Io,
                        message: format!(
                            "Could not remove new candidate-sidecar directory '{}' during rollback: {e}",
                            final_dir.display()
                        ),
                    })?;
                }
                Ok(())
            }
            CandidateSidecarCommitted::Removed {
                final_dir,
                backup_dir,
            } => {
                if let Some(backup) = backup_dir {
                    std::fs::rename(&backup, &final_dir).map_err(|e| EngineError {
                        code: ErrorCode::Io,
                        message: format!(
                            "Could not restore candidate-sidecar directory '{}' -> '{}': {e}",
                            backup.display(),
                            final_dir.display()
                        ),
                    })?;
                }
                Ok(())
            }
        }
    }

    fn finalize_candidate_store_transaction(
        committed: CandidateSidecarCommitted,
    ) -> Result<(), EngineError> {
        match committed {
            CandidateSidecarCommitted::Noop => Ok(()),
            CandidateSidecarCommitted::Replaced { backup_dir, .. }
            | CandidateSidecarCommitted::Removed { backup_dir, .. } => {
                if let Some(backup) = backup_dir {
                    std::fs::remove_dir_all(&backup).map_err(|e| EngineError {
                        code: ErrorCode::Io,
                        message: format!(
                            "Could not remove candidate-sidecar backup directory '{}': {e}",
                            backup.display()
                        ),
                    })?;
                }
                Ok(())
            }
        }
    }
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub enum ExportFormat {
    GenBank,
    Fasta,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub enum RenderSvgMode {
    Linear,
    Circular,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub enum PrimerLibraryMode {
    Enumerate,
    Sample,
}

#[derive(Debug, Clone, Copy, Serialize, Deserialize, PartialEq, Eq, Default)]
#[serde(rename_all = "snake_case")]
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

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PcrPrimerSpec {
    pub sequence: String,
    pub anneal_len: Option<usize>,
    pub max_mismatches: Option<usize>,
    pub require_3prime_exact_bases: Option<usize>,
    pub library_mode: Option<PrimerLibraryMode>,
    pub max_variants: Option<usize>,
    pub sample_seed: Option<u64>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SnpMutationSpec {
    pub zero_based_position: usize,
    pub reference: String,
    pub alternate: String,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default, deny_unknown_fields)]
pub struct PrimerDesignBaseLock {
    pub offset_0based: usize,
    pub base: String,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(default, deny_unknown_fields)]
pub struct PrimerDesignPairConstraint {
    pub require_roi_flanking: bool,
    #[serde(default)]
    pub required_amplicon_motifs: Vec<String>,
    #[serde(default)]
    pub forbidden_amplicon_motifs: Vec<String>,
    pub fixed_amplicon_start_0based: Option<usize>,
    pub fixed_amplicon_end_0based_exclusive: Option<usize>,
}

impl Default for PrimerDesignPairConstraint {
    fn default() -> Self {
        Self {
            require_roi_flanking: false,
            required_amplicon_motifs: vec![],
            forbidden_amplicon_motifs: vec![],
            fixed_amplicon_start_0based: None,
            fixed_amplicon_end_0based_exclusive: None,
        }
    }
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(default)]
pub struct PrimerDesignSideConstraint {
    pub min_length: usize,
    pub max_length: usize,
    pub location_0based: Option<usize>,
    pub start_0based: Option<usize>,
    pub end_0based: Option<usize>,
    pub min_tm_c: f64,
    pub max_tm_c: f64,
    pub min_gc_fraction: f64,
    pub max_gc_fraction: f64,
    pub max_anneal_hits: usize,
    pub non_annealing_5prime_tail: Option<String>,
    pub fixed_5prime: Option<String>,
    pub fixed_3prime: Option<String>,
    #[serde(default)]
    pub required_motifs: Vec<String>,
    #[serde(default)]
    pub forbidden_motifs: Vec<String>,
    #[serde(default)]
    pub locked_positions: Vec<PrimerDesignBaseLock>,
}

impl Default for PrimerDesignSideConstraint {
    fn default() -> Self {
        Self {
            min_length: PRIMER_PREFERRED_MIN_LENGTH_BP,
            max_length: PRIMER_PREFERRED_MAX_LENGTH_BP,
            location_0based: None,
            start_0based: None,
            end_0based: None,
            min_tm_c: 55.0,
            max_tm_c: 68.0,
            min_gc_fraction: 0.35,
            max_gc_fraction: 0.70,
            max_anneal_hits: 1,
            non_annealing_5prime_tail: None,
            fixed_5prime: None,
            fixed_3prime: None,
            required_motifs: vec![],
            forbidden_motifs: vec![],
            locked_positions: vec![],
        }
    }
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
pub struct PrimerDesignPrimerRecord {
    pub sequence: String,
    pub start_0based: usize,
    pub end_0based_exclusive: usize,
    pub length_bp: usize,
    pub anneal_length_bp: usize,
    pub non_annealing_5prime_tail_bp: usize,
    pub tm_c: f64,
    pub gc_fraction: f64,
    pub anneal_hits: usize,
    pub three_prime_base: String,
    pub three_prime_gc_clamp: bool,
    pub longest_homopolymer_run_bp: usize,
    pub self_complementary_run_bp: usize,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
pub struct PrimerDesignPairRuleFlags {
    pub roi_covered: bool,
    pub amplicon_size_in_range: bool,
    pub tm_delta_in_range: bool,
    pub forward_secondary_structure_ok: bool,
    pub reverse_secondary_structure_ok: bool,
    pub primer_pair_dimer_risk_low: bool,
    pub forward_three_prime_gc_clamp: bool,
    pub reverse_three_prime_gc_clamp: bool,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
pub struct PrimerDesignPairRecord {
    pub rank: usize,
    pub score: f64,
    pub forward: PrimerDesignPrimerRecord,
    pub reverse: PrimerDesignPrimerRecord,
    pub amplicon_start_0based: usize,
    pub amplicon_end_0based_exclusive: usize,
    pub amplicon_length_bp: usize,
    pub tm_delta_c: f64,
    pub primer_pair_complementary_run_bp: usize,
    pub primer_pair_3prime_complementary_run_bp: usize,
    pub rule_flags: PrimerDesignPairRuleFlags,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default, deny_unknown_fields)]
pub struct PrimerInsertionIntent {
    pub requested_forward_3prime_end_0based_exclusive: usize,
    pub requested_reverse_3prime_start_0based: usize,
    pub forward_extension_5prime: String,
    pub reverse_extension_5prime: String,
    pub forward_window_start_0based: usize,
    pub forward_window_end_0based_exclusive: usize,
    pub reverse_window_start_0based: usize,
    pub reverse_window_end_0based_exclusive: usize,
    pub max_anchor_shift_bp: Option<usize>,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
pub struct PrimerInsertionPairCompensation {
    pub rank: usize,
    pub forward_anchor_shift_bp: isize,
    pub reverse_anchor_shift_bp: isize,
    pub within_shift_budget: bool,
    pub compensable: bool,
    pub forward_compensation_5prime: String,
    pub reverse_compensation_5prime: String,
    pub compensated_forward_5prime_tail: String,
    pub compensated_reverse_5prime_tail: String,
    pub compensated_forward_sequence: String,
    pub compensated_reverse_sequence: String,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
pub struct PrimerInsertionContextReport {
    pub requested_forward_3prime_end_0based_exclusive: usize,
    pub requested_reverse_3prime_start_0based: usize,
    pub forward_extension_5prime: String,
    pub reverse_extension_5prime: String,
    pub forward_window_start_0based: usize,
    pub forward_window_end_0based_exclusive: usize,
    pub reverse_window_start_0based: usize,
    pub reverse_window_end_0based_exclusive: usize,
    pub max_anchor_shift_bp: usize,
    pub uncompensable_pair_count: usize,
    pub out_of_shift_budget_pair_count: usize,
    #[serde(default)]
    pub pairs: Vec<PrimerInsertionPairCompensation>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(default, deny_unknown_fields)]
pub struct OverlapExtensionMutagenesisConstraints {
    pub overlap_bp: usize,
    pub outer_forward: PrimerDesignSideConstraint,
    pub outer_reverse: PrimerDesignSideConstraint,
    pub inner_forward: PrimerDesignSideConstraint,
    pub inner_reverse: PrimerDesignSideConstraint,
}

impl Default for OverlapExtensionMutagenesisConstraints {
    fn default() -> Self {
        let side = PrimerDesignSideConstraint {
            min_length: 18,
            max_length: 30,
            min_tm_c: 50.0,
            max_tm_c: 72.0,
            min_gc_fraction: 0.30,
            max_gc_fraction: 0.75,
            max_anneal_hits: 1,
            ..PrimerDesignSideConstraint::default()
        };
        Self {
            overlap_bp: 24,
            outer_forward: side.clone(),
            outer_reverse: side.clone(),
            inner_forward: side.clone(),
            inner_reverse: side,
        }
    }
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
pub struct PrimerDesignRejectionSummary {
    pub out_of_window: usize,
    pub gc_or_tm_out_of_bounds: usize,
    pub non_unique_anneal: usize,
    pub amplicon_or_roi_failure: usize,
    pub primer_constraint_failure: usize,
    pub pair_constraint_failure: usize,
    pub pair_evaluation_limit_skipped: usize,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
pub struct PrimerDesignReport {
    pub schema: String,
    pub report_id: String,
    pub template: String,
    pub generated_at_unix_ms: u128,
    pub roi_start_0based: usize,
    pub roi_end_0based: usize,
    pub forward: PrimerDesignSideConstraint,
    pub reverse: PrimerDesignSideConstraint,
    #[serde(default)]
    pub pair_constraints: PrimerDesignPairConstraint,
    pub min_amplicon_bp: usize,
    pub max_amplicon_bp: usize,
    pub max_tm_delta_c: f64,
    pub max_pairs: usize,
    pub pair_count: usize,
    #[serde(default)]
    pub pairs: Vec<PrimerDesignPairRecord>,
    #[serde(default)]
    pub rejection_summary: PrimerDesignRejectionSummary,
    #[serde(default)]
    pub backend: PrimerDesignBackendInfo,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub insertion_context: Option<PrimerInsertionContextReport>,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
pub struct PrimerDesignBackendInfo {
    pub requested: String,
    pub used: String,
    pub fallback_reason: Option<String>,
    pub primer3_executable: Option<String>,
    pub primer3_version: Option<String>,
    pub primer3_explain: Option<String>,
    pub primer3_request_boulder_io: Option<String>,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
pub struct Primer3PreflightReport {
    pub backend: String,
    pub executable: String,
    pub reachable: bool,
    pub version_probe_ok: bool,
    pub status_code: Option<i32>,
    pub version: Option<String>,
    pub detail: Option<String>,
    pub error: Option<String>,
    pub probe_time_ms: u128,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
pub struct PrimerDesignReportSummary {
    pub report_id: String,
    pub template: String,
    pub generated_at_unix_ms: u128,
    pub roi_start_0based: usize,
    pub roi_end_0based: usize,
    pub pair_count: usize,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
pub struct QpcrAssayRuleFlags {
    pub roi_covered: bool,
    pub amplicon_size_in_range: bool,
    pub primer_tm_delta_in_range: bool,
    pub probe_inside_amplicon: bool,
    pub probe_tm_delta_in_range: bool,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
pub struct QpcrAssayRecord {
    pub rank: usize,
    pub score: f64,
    pub forward: PrimerDesignPrimerRecord,
    pub reverse: PrimerDesignPrimerRecord,
    pub probe: PrimerDesignPrimerRecord,
    pub amplicon_start_0based: usize,
    pub amplicon_end_0based_exclusive: usize,
    pub amplicon_length_bp: usize,
    pub primer_tm_delta_c: f64,
    pub probe_tm_delta_c: f64,
    pub rule_flags: QpcrAssayRuleFlags,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
pub struct QpcrDesignRejectionSummary {
    pub primer_pair: PrimerDesignRejectionSummary,
    pub probe_out_of_window: usize,
    pub probe_gc_or_tm_out_of_bounds: usize,
    pub probe_non_unique_anneal: usize,
    pub probe_or_assay_failure: usize,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
pub struct QpcrDesignReport {
    pub schema: String,
    pub report_id: String,
    pub template: String,
    pub generated_at_unix_ms: u128,
    pub roi_start_0based: usize,
    pub roi_end_0based: usize,
    pub forward: PrimerDesignSideConstraint,
    pub reverse: PrimerDesignSideConstraint,
    pub probe: PrimerDesignSideConstraint,
    #[serde(default)]
    pub pair_constraints: PrimerDesignPairConstraint,
    pub min_amplicon_bp: usize,
    pub max_amplicon_bp: usize,
    pub max_tm_delta_c: f64,
    pub max_probe_tm_delta_c: f64,
    pub max_assays: usize,
    pub assay_count: usize,
    #[serde(default)]
    pub assays: Vec<QpcrAssayRecord>,
    #[serde(default)]
    pub rejection_summary: QpcrDesignRejectionSummary,
    #[serde(default)]
    pub backend: PrimerDesignBackendInfo,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
pub struct QpcrDesignReportSummary {
    pub report_id: String,
    pub template: String,
    pub generated_at_unix_ms: u128,
    pub roi_start_0based: usize,
    pub roi_end_0based: usize,
    pub assay_count: usize,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct TfThresholdOverride {
    pub tf: String,
    pub min_llr_bits: Option<f64>,
    pub min_llr_quantile: Option<f64>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub enum LigationProtocol {
    Sticky,
    Blunt,
}

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq)]
pub enum AnchorBoundary {
    Start,
    End,
    Middle,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub enum AnchorDirection {
    Upstream,
    Downstream,
}

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq)]
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

// Backward-compatible alias kept while adapter/docs migrate to SequenceAnchor.
pub type AnchoredRegionAnchor = SequenceAnchor;

#[derive(Debug, Clone, Copy, Serialize, Deserialize, PartialEq, Eq, Default)]
#[serde(rename_all = "snake_case")]
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
pub struct GenomeTrackSubscription {
    pub source: GenomeTrackSource,
    pub path: String,
    pub track_name: Option<String>,
    pub min_score: Option<f64>,
    pub max_score: Option<f64>,
    pub clear_existing: bool,
}

impl Default for GenomeTrackSubscription {
    fn default() -> Self {
        Self {
            source: GenomeTrackSource::Bed,
            path: String::new(),
            track_name: None,
            min_score: None,
            max_score: None,
            clear_existing: false,
        }
    }
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
pub struct GenomeTrackSyncReport {
    pub subscriptions_considered: usize,
    pub target_sequences: usize,
    pub applied_imports: usize,
    pub failed_imports: usize,
    pub warnings_count: usize,
    pub errors: Vec<String>,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
struct CandidateStore {
    schema: String,
    updated_at_unix_ms: u128,
    sets: HashMap<String, CandidateSet>,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
struct CandidateStoreReference {
    schema: String,
    storage: String,
    index_path: String,
    set_count: usize,
    updated_at_unix_ms: u128,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
struct CandidateStoreDiskIndex {
    schema: String,
    updated_at_unix_ms: u128,
    set_count: usize,
    sets: Vec<CandidateStoreDiskSetIndexEntry>,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
struct CandidateStoreDiskSetIndexEntry {
    name: String,
    created_at_unix_ms: u128,
    source_seq_ids: Vec<String>,
    candidate_count: usize,
    metrics: Vec<String>,
    records_path: String,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
/// Persisted named collection of candidate intervals plus their scored metrics.
///
/// Candidate-set operations in GUI/CLI/macros exchange this record rather than
/// adapter-local tables. Start here when tracing candidate generation,
/// filtering, or top-k/pareto workflows.
pub struct CandidateSet {
    pub name: String,
    pub created_at_unix_ms: u128,
    pub source_seq_ids: Vec<String>,
    pub candidates: Vec<CandidateRecord>,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
/// One concrete candidate interval extracted from a source sequence.
///
/// Coordinates are 0-based half-open on `seq_id`; `metrics` stores arbitrary
/// numeric scores produced by candidate analysis operations.
pub struct CandidateRecord {
    pub seq_id: String,
    pub start_0based: usize,
    pub end_0based: usize,
    pub sequence: String,
    pub metrics: HashMap<String, f64>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
/// Lightweight listing row for one candidate set without materializing records.
pub struct CandidateSetSummary {
    pub name: String,
    pub created_at_unix_ms: u128,
    pub source_seq_ids: Vec<String>,
    pub candidate_count: usize,
    pub metrics: Vec<String>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
/// Coverage summary for one metric across all rows in a candidate set.
pub struct CandidateMetricSummary {
    pub metric: String,
    pub present_in_candidates: usize,
    pub missing_in_candidates: usize,
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

#[derive(Debug, Clone, Copy, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
/// Deterministic set algebra supported by candidate-set combination commands.
pub enum CandidateSetOperator {
    Union,
    Intersect,
    Subtract,
}

impl CandidateSetOperator {
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
/// Weighted scalar objective term used by `CandidatesScoreWeightedObjective`.
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

#[derive(Debug, Clone, Copy, Serialize, Deserialize, PartialEq, Eq, Default)]
#[serde(rename_all = "snake_case")]
pub enum GuideU6TerminatorWindow {
    SpacerOnly,
    #[default]
    SpacerPlusTail,
}

impl GuideU6TerminatorWindow {
    pub fn as_str(self) -> &'static str {
        match self {
            Self::SpacerOnly => "spacer_only",
            Self::SpacerPlusTail => "spacer_plus_tail",
        }
    }
}

#[derive(Debug, Clone, Copy, Serialize, Deserialize, PartialEq, Eq)]
#[serde(rename_all = "snake_case")]
pub enum GuideOligoExportFormat {
    CsvTable,
    PlateCsv,
    Fasta,
}

impl GuideOligoExportFormat {
    pub fn as_str(self) -> &'static str {
        match self {
            Self::CsvTable => "csv_table",
            Self::PlateCsv => "plate_csv",
            Self::Fasta => "fasta",
        }
    }
}

#[derive(Debug, Clone, Copy, Serialize, Deserialize, PartialEq, Eq, Default)]
#[serde(rename_all = "snake_case")]
pub enum GuideOligoPlateFormat {
    #[default]
    Plate96,
    Plate384,
}

impl GuideOligoPlateFormat {
    pub fn as_str(self) -> &'static str {
        match self {
            Self::Plate96 => "96",
            Self::Plate384 => "384",
        }
    }

    fn dimensions(self) -> (usize, usize) {
        match self {
            Self::Plate96 => (8, 12),
            Self::Plate384 => (16, 24),
        }
    }
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
struct GuideDesignStore {
    schema: String,
    updated_at_unix_ms: u128,
    guide_sets: HashMap<String, GuideSet>,
    practical_filter_reports: HashMap<String, GuidePracticalFilterReport>,
    oligo_sets: HashMap<String, GuideOligoSet>,
    latest_oligo_set_by_guide_set: HashMap<String, String>,
    audit_log: Vec<GuideDesignAuditEntry>,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
struct PrimerDesignStore {
    schema: String,
    updated_at_unix_ms: u128,
    reports: HashMap<String, PrimerDesignReport>,
    qpcr_reports: HashMap<String, QpcrDesignReport>,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
struct SequencingTraceStore {
    schema: String,
    updated_at_unix_ms: u128,
    traces: BTreeMap<String, SequencingTraceRecord>,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
struct SequencingConfirmationReportStore {
    schema: String,
    updated_at_unix_ms: u128,
    reports: BTreeMap<String, SequencingConfirmationReport>,
}

#[derive(Debug, Clone, Copy, Serialize, Deserialize, PartialEq, Eq, Default)]
#[serde(rename_all = "snake_case")]
/// Which stored planning profile layer a read/write command refers to.
///
/// The effective merge order is `global -> confirmed_agent_overlay ->
/// project_override`.
pub enum PlanningProfileScope {
    #[default]
    Global,
    ProjectOverride,
    ConfirmedAgentOverlay,
    Effective,
}

impl PlanningProfileScope {
    pub fn as_str(self) -> &'static str {
        match self {
            Self::Global => "global",
            Self::ProjectOverride => "project_override",
            Self::ConfirmedAgentOverlay => "confirmed_agent_overlay",
            Self::Effective => "effective",
        }
    }
}

#[derive(Debug, Clone, Copy, Serialize, Deserialize, PartialEq, Eq, Default)]
#[serde(rename_all = "snake_case")]
/// Review state for an advisory planning suggestion.
pub enum PlanningSuggestionStatus {
    #[default]
    Pending,
    Accepted,
    Rejected,
}

impl PlanningSuggestionStatus {
    pub fn as_str(self) -> &'static str {
        match self {
            Self::Pending => "pending",
            Self::Accepted => "accepted",
            Self::Rejected => "rejected",
        }
    }
}

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
#[serde(default, deny_unknown_fields)]
/// Availability/cost note for one inventory item referenced by planning.
pub struct PlanningInventoryItem {
    pub available: bool,
    pub unit_cost: Option<f64>,
    pub procurement_business_days: Option<f64>,
    pub note: Option<String>,
}

impl Default for PlanningInventoryItem {
    fn default() -> Self {
        Self {
            available: true,
            unit_cost: None,
            procurement_business_days: None,
            note: None,
        }
    }
}

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
#[serde(default, deny_unknown_fields)]
/// Availability/queue note for one machine or platform used in planning.
pub struct PlanningMachineAvailability {
    pub available: bool,
    pub queue_business_days: f64,
    pub run_cost_per_hour: Option<f64>,
    pub note: Option<String>,
}

impl Default for PlanningMachineAvailability {
    fn default() -> Self {
        Self {
            available: true,
            queue_business_days: 0.0,
            run_cost_per_hour: None,
            note: None,
        }
    }
}

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
#[serde(default, deny_unknown_fields)]
/// User- or agent-supplied planning profile describing local capabilities.
///
/// This is the main durable contract for inventory, procurement latency, and
/// machine availability. Start here when tracing `planning profile ...`.
pub struct PlanningProfile {
    pub schema: String,
    pub profile_id: Option<String>,
    pub currency: Option<String>,
    pub procurement_business_days_default: f64,
    pub capabilities: Vec<String>,
    pub inventory: HashMap<String, PlanningInventoryItem>,
    pub machine_availability: HashMap<String, PlanningMachineAvailability>,
    pub notes: Option<String>,
}

impl Default for PlanningProfile {
    fn default() -> Self {
        Self {
            schema: PLANNING_PROFILE_SCHEMA.to_string(),
            profile_id: None,
            currency: None,
            procurement_business_days_default: 10.0,
            capabilities: vec![],
            inventory: HashMap::new(),
            machine_availability: HashMap::new(),
            notes: None,
        }
    }
}

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
#[serde(default, deny_unknown_fields)]
/// Optimization weights and hard guardrails used by planning estimates.
pub struct PlanningObjective {
    pub schema: String,
    pub weight_time: f64,
    pub weight_cost: f64,
    pub weight_local_fit: f64,
    pub max_cost: Option<f64>,
    pub max_time_hours: Option<f64>,
    pub required_capabilities: Vec<String>,
    pub enforce_guardrails: bool,
}

impl Default for PlanningObjective {
    fn default() -> Self {
        Self {
            schema: PLANNING_OBJECTIVE_SCHEMA.to_string(),
            weight_time: 1.0,
            weight_cost: 1.0,
            weight_local_fit: 1.0,
            max_cost: None,
            max_time_hours: None,
            required_capabilities: vec![],
            enforce_guardrails: true,
        }
    }
}

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
#[serde(default, deny_unknown_fields)]
/// Deterministic estimate payload produced by planning evaluation.
pub struct PlanningEstimate {
    pub schema: String,
    pub estimated_time_hours: f64,
    pub estimated_cost: f64,
    pub local_fit_score: f64,
    pub composite_meta_score: f64,
    pub passes_guardrails: bool,
    pub guardrail_failures: Vec<String>,
    pub explanation: serde_json::Value,
}

impl Default for PlanningEstimate {
    fn default() -> Self {
        Self {
            schema: PLANNING_ESTIMATE_SCHEMA.to_string(),
            estimated_time_hours: 0.0,
            estimated_cost: 0.0,
            local_fit_score: 1.0,
            composite_meta_score: 0.0,
            passes_guardrails: true,
            guardrail_failures: vec![],
            explanation: json!({}),
        }
    }
}

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
#[serde(default, deny_unknown_fields)]
/// Advisory pull/push suggestion awaiting explicit user acceptance or rejection.
pub struct PlanningSuggestion {
    pub schema: String,
    pub suggestion_id: String,
    pub status: PlanningSuggestionStatus,
    pub direction: String,
    pub source: String,
    pub confidence: Option<f64>,
    pub snapshot_id: Option<String>,
    pub message: Option<String>,
    pub profile_patch: Option<PlanningProfile>,
    pub objective_patch: Option<PlanningObjective>,
    pub diff: serde_json::Value,
    pub created_at_unix_ms: u128,
    pub resolved_at_unix_ms: Option<u128>,
    pub rejection_reason: Option<String>,
}

impl Default for PlanningSuggestion {
    fn default() -> Self {
        Self {
            schema: PLANNING_SUGGESTION_SCHEMA.to_string(),
            suggestion_id: String::new(),
            status: PlanningSuggestionStatus::Pending,
            direction: "pull".to_string(),
            source: String::new(),
            confidence: None,
            snapshot_id: None,
            message: None,
            profile_patch: None,
            objective_patch: None,
            diff: json!({}),
            created_at_unix_ms: 0,
            resolved_at_unix_ms: None,
            rejection_reason: None,
        }
    }
}

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
#[serde(default, deny_unknown_fields)]
/// High-level sync status for the planning suggestion channel.
pub struct PlanningSyncStatus {
    pub schema: String,
    pub pending_suggestion_count: usize,
    pub last_pull_at_unix_ms: Option<u128>,
    pub last_push_at_unix_ms: Option<u128>,
    pub last_source: Option<String>,
    pub last_snapshot_id: Option<String>,
    pub last_error: Option<String>,
}

impl Default for PlanningSyncStatus {
    fn default() -> Self {
        Self {
            schema: PLANNING_SYNC_STATUS_SCHEMA.to_string(),
            pending_suggestion_count: 0,
            last_pull_at_unix_ms: None,
            last_push_at_unix_ms: None,
            last_source: None,
            last_snapshot_id: None,
            last_error: None,
        }
    }
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
struct PlanningStore {
    schema: String,
    updated_at_unix_ms: u128,
    global_profile: Option<PlanningProfile>,
    project_override_profile: Option<PlanningProfile>,
    confirmed_agent_overlay_profile: Option<PlanningProfile>,
    objective: Option<PlanningObjective>,
    suggestions: HashMap<String, PlanningSuggestion>,
    sync_status: PlanningSyncStatus,
    next_suggestion_counter: u64,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
struct DotplotAnalysisStore {
    schema: String,
    updated_at_unix_ms: u128,
    dotplots: HashMap<String, DotplotView>,
    flexibility_tracks: HashMap<String, FlexibilityTrack>,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default, deny_unknown_fields)]
pub struct IsoformPanelDomainSpec {
    pub name: String,
    pub start_aa: usize,
    pub end_aa: usize,
    #[serde(default)]
    pub color_hex: Option<String>,
}

#[derive(Debug, Clone, Copy, Serialize, Deserialize, PartialEq, Eq, Default)]
#[serde(rename_all = "snake_case")]
pub enum IsoformTranscriptGeometryMode {
    #[default]
    Exon,
    Cds,
}

impl IsoformTranscriptGeometryMode {
    pub fn as_str(self) -> &'static str {
        match self {
            Self::Exon => "exon",
            Self::Cds => "cds",
        }
    }
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default, deny_unknown_fields)]
pub struct IsoformPanelIsoformSpec {
    pub isoform_id: String,
    #[serde(default)]
    pub label: Option<String>,
    #[serde(default)]
    pub transcript_ids: Vec<String>,
    #[serde(default)]
    pub transactivation_class: Option<String>,
    #[serde(default)]
    pub expected_length_aa: Option<usize>,
    #[serde(default)]
    pub reference_start_aa: Option<usize>,
    #[serde(default)]
    pub reference_end_aa: Option<usize>,
    #[serde(default)]
    pub domains: Vec<IsoformPanelDomainSpec>,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default, deny_unknown_fields)]
pub struct IsoformPanelResource {
    pub schema: String,
    pub panel_id: String,
    pub gene_symbol: String,
    #[serde(default)]
    pub transcript_geometry_mode: IsoformTranscriptGeometryMode,
    #[serde(default)]
    pub assembly: Option<String>,
    #[serde(default)]
    pub source: Option<String>,
    #[serde(default)]
    pub notes: Option<String>,
    #[serde(default)]
    pub isoforms: Vec<IsoformPanelIsoformSpec>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(default, deny_unknown_fields)]
pub struct IsoformPanelValidationIssue {
    pub severity: String,
    pub code: String,
    pub message: String,
    pub isoform_id: Option<String>,
    pub transcript_probe: Option<String>,
    pub domain_name: Option<String>,
}

impl Default for IsoformPanelValidationIssue {
    fn default() -> Self {
        Self {
            severity: "warning".to_string(),
            code: String::new(),
            message: String::new(),
            isoform_id: None,
            transcript_probe: None,
            domain_name: None,
        }
    }
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default, deny_unknown_fields)]
pub struct IsoformPanelValidationIsoformSummary {
    pub isoform_id: String,
    pub label: String,
    pub transcript_probe_count: usize,
    pub domain_count: usize,
    pub expected_length_aa: Option<usize>,
    pub max_domain_end_aa: Option<usize>,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default, deny_unknown_fields)]
pub struct IsoformPanelValidationReport {
    pub schema: String,
    pub path: String,
    pub panel_id: String,
    pub gene_symbol: String,
    pub assembly: Option<String>,
    pub isoform_count: usize,
    pub transcript_probe_count: usize,
    pub unique_transcript_probe_count: usize,
    pub domain_count: usize,
    pub issue_count: usize,
    pub status: String,
    pub isoforms: Vec<IsoformPanelValidationIsoformSummary>,
    pub issues: Vec<IsoformPanelValidationIssue>,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default, deny_unknown_fields)]
struct IsoformPanelRecord {
    seq_id: String,
    panel_id: String,
    imported_at_unix_ms: u128,
    source_path: String,
    strict: bool,
    resource: IsoformPanelResource,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default, deny_unknown_fields)]
struct IsoformPanelStore {
    schema: String,
    updated_at_unix_ms: u128,
    records: Vec<IsoformPanelRecord>,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default, deny_unknown_fields)]
struct UniprotEntryStore {
    schema: String,
    updated_at_unix_ms: u128,
    entries: HashMap<String, UniprotEntry>,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default, deny_unknown_fields)]
struct UniprotGenomeProjectionStore {
    schema: String,
    updated_at_unix_ms: u128,
    projections: HashMap<String, UniprotGenomeProjection>,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
pub struct GuideSet {
    pub guide_set_id: String,
    pub guides: Vec<GuideCandidate>,
    pub created_at_unix_ms: u128,
    pub updated_at_unix_ms: u128,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
pub struct GuideCandidate {
    pub guide_id: String,
    pub seq_id: String,
    pub start_0based: usize,
    pub end_0based_exclusive: usize,
    pub strand: String,
    pub protospacer: String,
    pub pam: String,
    pub nuclease: String,
    pub cut_offset_from_protospacer_start: usize,
    pub rank: Option<usize>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct GuideSetSummary {
    pub guide_set_id: String,
    pub guide_count: usize,
    pub created_at_unix_ms: u128,
    pub updated_at_unix_ms: u128,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(default)]
pub struct GuidePracticalFilterConfig {
    pub gc_min: Option<f64>,
    pub gc_max: Option<f64>,
    pub max_homopolymer_run: Option<usize>,
    #[serde(default)]
    pub max_homopolymer_run_per_base: HashMap<String, usize>,
    pub reject_ambiguous_bases: bool,
    pub avoid_u6_terminator_tttt: bool,
    pub u6_terminator_window: GuideU6TerminatorWindow,
    pub max_dinucleotide_repeat_units: Option<usize>,
    #[serde(default)]
    pub forbidden_motifs: Vec<String>,
    pub required_5prime_base: Option<String>,
    pub allow_5prime_g_extension: bool,
}

impl Default for GuidePracticalFilterConfig {
    fn default() -> Self {
        Self {
            gc_min: None,
            gc_max: None,
            max_homopolymer_run: None,
            max_homopolymer_run_per_base: HashMap::new(),
            reject_ambiguous_bases: true,
            avoid_u6_terminator_tttt: true,
            u6_terminator_window: GuideU6TerminatorWindow::SpacerPlusTail,
            max_dinucleotide_repeat_units: None,
            forbidden_motifs: vec![],
            required_5prime_base: None,
            allow_5prime_g_extension: false,
        }
    }
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
pub struct GuideFilterReason {
    pub code: String,
    pub message: String,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
pub struct GuidePracticalFilterResult {
    pub guide_id: String,
    pub passed: bool,
    #[serde(default)]
    pub reasons: Vec<GuideFilterReason>,
    #[serde(default)]
    pub warnings: Vec<String>,
    #[serde(default)]
    pub metrics: HashMap<String, f64>,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
pub struct GuidePracticalFilterReport {
    pub guide_set_id: String,
    pub generated_at_unix_ms: u128,
    pub config: GuidePracticalFilterConfig,
    pub passed_count: usize,
    pub rejected_count: usize,
    #[serde(default)]
    pub results: Vec<GuidePracticalFilterResult>,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
pub struct GuideOligoTemplateSpec {
    pub template_id: String,
    pub description: String,
    pub forward_prefix: String,
    pub forward_suffix: String,
    pub reverse_prefix: String,
    pub reverse_suffix: String,
    pub reverse_uses_reverse_complement_of_spacer: bool,
    pub uppercase_output: bool,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
pub struct GuideOligoRecord {
    pub guide_id: String,
    pub rank: Option<usize>,
    pub forward_oligo: String,
    pub reverse_oligo: String,
    #[serde(default)]
    pub notes: Vec<String>,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
pub struct GuideOligoSet {
    pub oligo_set_id: String,
    pub guide_set_id: String,
    pub generated_at_unix_ms: u128,
    pub template: GuideOligoTemplateSpec,
    pub apply_5prime_g_extension: bool,
    #[serde(default)]
    pub records: Vec<GuideOligoRecord>,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
pub struct GuideOligoExportReport {
    pub guide_set_id: String,
    pub oligo_set_id: String,
    pub format: String,
    pub path: String,
    pub exported_records: usize,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
pub struct GuideProtocolExportReport {
    pub guide_set_id: String,
    pub oligo_set_id: String,
    pub path: String,
    pub guide_count: usize,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
struct GuideDesignAuditEntry {
    unix_ms: u128,
    operation: String,
    guide_set_id: String,
    payload: serde_json::Value,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
struct WorkflowMacroTemplateStore {
    schema: String,
    updated_at_unix_ms: u128,
    templates: HashMap<String, WorkflowMacroTemplate>,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
/// Stored workflow-macro template shared by shell, GUI helpers, and future MCP.
///
/// A workflow macro expands to one workflow script plus typed parameters and
/// declared input/output ports for lineage/protocol reporting.
pub struct WorkflowMacroTemplate {
    pub name: String,
    pub description: Option<String>,
    pub details_url: Option<String>,
    pub parameters: Vec<WorkflowMacroTemplateParam>,
    #[serde(default)]
    pub input_ports: Vec<WorkflowMacroTemplatePort>,
    #[serde(default)]
    pub output_ports: Vec<WorkflowMacroTemplatePort>,
    #[serde(default = "default_cloning_macro_template_schema")]
    pub template_schema: String,
    pub script: String,
    pub created_at_unix_ms: u128,
    pub updated_at_unix_ms: u128,
}

fn default_cloning_macro_template_schema() -> String {
    CLONING_MACRO_TEMPLATE_SCHEMA.to_string()
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
/// One typed parameter exposed by a workflow macro template.
pub struct WorkflowMacroTemplateParam {
    pub name: String,
    pub default_value: Option<String>,
    pub required: bool,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
/// Declared input or output port on a workflow macro template.
pub struct WorkflowMacroTemplatePort {
    pub port_id: String,
    pub kind: String,
    pub required: bool,
    pub cardinality: String,
    pub description: Option<String>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
/// Listing row for one workflow macro template.
pub struct WorkflowMacroTemplateSummary {
    pub name: String,
    pub description: Option<String>,
    pub details_url: Option<String>,
    pub parameter_count: usize,
    pub created_at_unix_ms: u128,
    pub updated_at_unix_ms: u128,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
struct CandidateMacroTemplateStore {
    schema: String,
    updated_at_unix_ms: u128,
    templates: HashMap<String, CandidateMacroTemplate>,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
/// Stored candidate-macro template that expands to candidate-shell commands.
pub struct CandidateMacroTemplate {
    pub name: String,
    pub description: Option<String>,
    pub details_url: Option<String>,
    pub parameters: Vec<CandidateMacroTemplateParam>,
    pub script: String,
    pub created_at_unix_ms: u128,
    pub updated_at_unix_ms: u128,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
/// One typed parameter exposed by a candidate macro template.
pub struct CandidateMacroTemplateParam {
    pub name: String,
    pub default_value: Option<String>,
    pub required: bool,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
/// Listing row for one candidate macro template.
pub struct CandidateMacroTemplateSummary {
    pub name: String,
    pub description: Option<String>,
    pub details_url: Option<String>,
    pub parameter_count: usize,
    pub created_at_unix_ms: u128,
    pub updated_at_unix_ms: u128,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
/// BLAST hit reduced to the fields needed for feature import/overlay pipelines.
pub struct BlastHitFeatureInput {
    pub subject_id: String,
    pub query_start_1based: usize,
    pub query_end_1based: usize,
    pub subject_start_1based: usize,
    pub subject_end_1based: usize,
    pub identity_percent: f64,
    pub bit_score: f64,
    pub evalue: f64,
    pub query_coverage_percent: Option<f64>,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
/// Provenance bundle describing exactly how a BLAST search was invoked.
///
/// This is the record to inspect when reproducing a search or explaining which
/// executable/options/catalog paths produced a report.
pub struct BlastInvocationProvenance {
    pub genome_id: String,
    pub query_label: String,
    pub query_length: usize,
    pub max_hits: usize,
    pub task: String,
    pub blastn_executable: String,
    pub blast_db_prefix: String,
    pub command: Vec<String>,
    pub command_line: String,
    pub catalog_path: Option<String>,
    pub cache_dir: Option<String>,
    #[serde(default)]
    pub options_override_json: Option<serde_json::Value>,
    #[serde(default)]
    pub effective_options_json: Option<serde_json::Value>,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default, deny_unknown_fields)]
/// Optional threshold overrides for genome BLAST post-filtering.
pub struct BlastThresholdOptions {
    pub max_evalue: Option<f64>,
    pub min_identity_percent: Option<f64>,
    pub min_query_coverage_percent: Option<f64>,
    pub min_alignment_length_bp: Option<usize>,
    pub min_bit_score: Option<f64>,
    pub unique_best_hit: Option<bool>,
}

impl BlastThresholdOptions {
    fn merge_from(&mut self, other: &Self) {
        if other.max_evalue.is_some() {
            self.max_evalue = other.max_evalue;
        }
        if other.min_identity_percent.is_some() {
            self.min_identity_percent = other.min_identity_percent;
        }
        if other.min_query_coverage_percent.is_some() {
            self.min_query_coverage_percent = other.min_query_coverage_percent;
        }
        if other.min_alignment_length_bp.is_some() {
            self.min_alignment_length_bp = other.min_alignment_length_bp;
        }
        if other.min_bit_score.is_some() {
            self.min_bit_score = other.min_bit_score;
        }
        if other.unique_best_hit.is_some() {
            self.unique_best_hit = other.unique_best_hit;
        }
    }
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default, deny_unknown_fields)]
/// Request-time BLAST execution options before defaults/project overrides merge.
pub struct BlastRunOptions {
    pub task: Option<String>,
    pub max_hits: Option<usize>,
    #[serde(default)]
    pub thresholds: BlastThresholdOptions,
}

impl BlastRunOptions {
    fn merge_from(&mut self, other: &Self) {
        if let Some(task) = other.task.as_ref() {
            self.task = Some(task.clone());
        }
        if other.max_hits.is_some() {
            self.max_hits = other.max_hits;
        }
        self.thresholds.merge_from(&other.thresholds);
    }
}

#[derive(Debug, Clone, Serialize, Deserialize)]
/// Resolved BLAST options after defaults/project/request layering.
pub struct BlastResolvedOptions {
    pub task: String,
    pub max_hits: usize,
    #[serde(default)]
    pub thresholds: BlastThresholdOptions,
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
    fn as_str(self) -> &'static str {
        match self {
            Self::FivePrime => "5prime",
            Self::ThreePrime => "3prime",
        }
    }
}

#[derive(Debug, Clone, Copy, Serialize, Deserialize, PartialEq, Eq)]
#[serde(rename_all = "snake_case")]
/// Annotation projection policy for `ExtractGenomeRegion`.
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
/// Interval policy for `ExtractGenomeGene`.
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
/// Verdict assigned to one confirmation target, read, or whole report.
pub enum SequencingConfirmationStatus {
    Confirmed,
    Contradicted,
    #[default]
    InsufficientEvidence,
}

impl SequencingConfirmationStatus {
    pub fn as_str(self) -> &'static str {
        match self {
            Self::Confirmed => "confirmed",
            Self::Contradicted => "contradicted",
            Self::InsufficientEvidence => "insufficient_evidence",
        }
    }
}

#[derive(Debug, Clone, Copy, Serialize, Deserialize, PartialEq, Eq, Default)]
#[serde(rename_all = "snake_case")]
/// Orientation chosen for the best alignment of one confirmation read.
pub enum SequencingReadOrientation {
    #[default]
    Forward,
    ReverseComplement,
}

impl SequencingReadOrientation {
    pub fn as_str(self) -> &'static str {
        match self {
            Self::Forward => "forward",
            Self::ReverseComplement => "reverse_complement",
        }
    }
}

#[derive(Debug, Clone, Copy, Serialize, Deserialize, PartialEq, Eq, Default)]
#[serde(rename_all = "snake_case")]
/// Origin of one sequencing-confirmation evidence row.
pub enum SequencingConfirmationEvidenceKind {
    #[default]
    Sequence,
    Trace,
}

impl SequencingConfirmationEvidenceKind {
    pub fn as_str(self) -> &'static str {
        match self {
            Self::Sequence => "sequence",
            Self::Trace => "trace",
        }
    }
}

#[derive(Debug, Clone, Copy, Serialize, Deserialize, PartialEq, Eq, Default)]
#[serde(rename_all = "snake_case")]
/// Biology-facing target classes supported by sequencing confirmation v1.
pub enum SequencingConfirmationTargetKind {
    #[default]
    FullSpan,
    Junction,
    FeaturePresence,
    ExpectedEdit,
    RestrictionSite,
}

impl SequencingConfirmationTargetKind {
    pub fn as_str(self) -> &'static str {
        match self {
            Self::FullSpan => "full_span",
            Self::Junction => "junction",
            Self::FeaturePresence => "feature_presence",
            Self::ExpectedEdit => "expected_edit",
            Self::RestrictionSite => "restriction_site",
        }
    }
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(default)]
/// One requested construct-confirmation target on the expected sequence.
///
/// Coordinates are 0-based half-open on the expected construct. Junction
/// targets also set `junction_left_end_0based` to mark the break between the
/// two halves that must both be supported.
pub struct SequencingConfirmationTargetSpec {
    pub target_id: String,
    pub label: String,
    pub kind: SequencingConfirmationTargetKind,
    pub start_0based: usize,
    pub end_0based_exclusive: usize,
    pub junction_left_end_0based: Option<usize>,
    pub required: bool,
}

impl Default for SequencingConfirmationTargetSpec {
    fn default() -> Self {
        Self {
            target_id: String::new(),
            label: String::new(),
            kind: SequencingConfirmationTargetKind::FullSpan,
            start_0based: 0,
            end_0based_exclusive: 0,
            junction_left_end_0based: None,
            required: true,
        }
    }
}

#[derive(Debug, Clone, Copy, Serialize, Deserialize, PartialEq, Eq, Default)]
#[serde(rename_all = "snake_case")]
/// Alignment-level discrepancy class recorded for one confirmation read.
pub enum SequencingConfirmationDiscrepancyKind {
    #[default]
    Mismatch,
    Insertion,
    Deletion,
}

impl SequencingConfirmationDiscrepancyKind {
    pub fn as_str(self) -> &'static str {
        match self {
            Self::Mismatch => "mismatch",
            Self::Insertion => "insertion",
            Self::Deletion => "deletion",
        }
    }
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
/// One mismatch/indel segment extracted from the best read alignment.
pub struct SequencingConfirmationDiscrepancy {
    pub kind: SequencingConfirmationDiscrepancyKind,
    pub query_start_0based: usize,
    pub query_end_0based_exclusive: usize,
    pub target_start_0based: usize,
    pub target_end_0based_exclusive: usize,
    pub query_bases: String,
    pub target_bases: String,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
/// Construct-level verdict for one requested confirmation target.
pub struct SequencingConfirmationTargetResult {
    pub target_id: String,
    pub label: String,
    pub kind: SequencingConfirmationTargetKind,
    pub start_0based: usize,
    pub end_0based_exclusive: usize,
    pub junction_left_end_0based: Option<usize>,
    pub required: bool,
    pub status: SequencingConfirmationStatus,
    pub covered_bp: usize,
    pub target_length_bp: usize,
    pub support_read_ids: Vec<String>,
    pub contradicting_read_ids: Vec<String>,
    pub reason: String,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
/// Best-alignment and discrepancy summary for one input sequencing evidence row.
pub struct SequencingConfirmationReadResult {
    pub evidence_kind: SequencingConfirmationEvidenceKind,
    pub evidence_id: String,
    pub read_seq_id: String,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub trace_id: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub linked_seq_id: Option<String>,
    pub orientation: SequencingReadOrientation,
    pub usable: bool,
    pub best_alignment: SequenceAlignmentReport,
    pub covered_target_ids: Vec<String>,
    pub confirmed_target_ids: Vec<String>,
    pub contradicted_target_ids: Vec<String>,
    pub discrepancies: Vec<SequencingConfirmationDiscrepancy>,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
/// Stored construct-confirmation report shared by shell, GUI, and exports.
///
/// Start here when following the sequencing-confirmation workflow end to end:
/// it ties together requested targets, per-evidence alignments, and the overall
/// confirmed/contradicted/insufficient-evidence verdict.
pub struct SequencingConfirmationReport {
    pub schema: String,
    pub report_id: String,
    pub expected_seq_id: String,
    pub generated_at_unix_ms: u128,
    pub overall_status: SequencingConfirmationStatus,
    pub alignment_mode: PairwiseAlignmentMode,
    pub match_score: i32,
    pub mismatch_score: i32,
    pub gap_open: i32,
    pub gap_extend: i32,
    pub min_identity_fraction: f64,
    pub min_target_coverage_fraction: f64,
    pub allow_reverse_complement: bool,
    pub read_seq_ids: Vec<String>,
    #[serde(default)]
    pub trace_ids: Vec<String>,
    pub target_count: usize,
    pub reads: Vec<SequencingConfirmationReadResult>,
    pub targets: Vec<SequencingConfirmationTargetResult>,
    pub warnings: Vec<String>,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
/// Lightweight listing row for one sequencing-confirmation report.
pub struct SequencingConfirmationReportSummary {
    pub report_id: String,
    pub expected_seq_id: String,
    pub generated_at_unix_ms: u128,
    pub overall_status: SequencingConfirmationStatus,
    pub read_count: usize,
    pub target_count: usize,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
pub struct FlexibilityBinScore {
    pub start_0based: usize,
    pub end_0based_exclusive: usize,
    pub score: f64,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
pub struct FlexibilityTrack {
    pub schema: String,
    pub track_id: String,
    pub seq_id: String,
    pub generated_at_unix_ms: u128,
    pub span_start_0based: usize,
    pub span_end_0based: usize,
    pub model: FlexibilityModel,
    pub bin_bp: usize,
    pub smoothing_bp: Option<usize>,
    pub min_score: f64,
    pub max_score: f64,
    pub bins: Vec<FlexibilityBinScore>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct FlexibilityTrackSummary {
    pub track_id: String,
    pub seq_id: String,
    pub generated_at_unix_ms: u128,
    pub span_start_0based: usize,
    pub span_end_0based: usize,
    pub model: FlexibilityModel,
    pub bin_bp: usize,
    pub smoothing_bp: Option<usize>,
    pub bin_count: usize,
    pub min_score: f64,
    pub max_score: f64,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
/// Canonical engine operation contract.
///
/// All adapters (GUI/CLI/JS/Lua/MCP) should map user intent to this enum and
/// rely on `GentleEngine::apply` for execution. This preserves one deterministic
/// behavior surface and avoids adapter-specific biology logic branches.
pub enum Operation {
    LoadFile {
        path: String,
        as_id: Option<SeqId>,
    },
    SaveFile {
        seq_id: SeqId,
        path: String,
        format: ExportFormat,
    },
    RenderSequenceSvg {
        seq_id: SeqId,
        mode: RenderSvgMode,
        path: String,
    },
    RenderDotplotSvg {
        seq_id: SeqId,
        dotplot_id: String,
        path: String,
        #[serde(default)]
        flex_track_id: Option<String>,
        #[serde(default)]
        display_density_threshold: Option<f32>,
        #[serde(default)]
        display_intensity_gain: Option<f32>,
    },
    RenderFeatureExpertSvg {
        seq_id: SeqId,
        target: FeatureExpertTarget,
        path: String,
    },
    RenderIsoformArchitectureSvg {
        seq_id: SeqId,
        panel_id: String,
        path: String,
    },
    RenderRnaStructureSvg {
        seq_id: SeqId,
        path: String,
    },
    RenderLineageSvg {
        path: String,
    },
    RenderPoolGelSvg {
        inputs: Vec<SeqId>,
        path: String,
        ladders: Option<Vec<String>>,
        #[serde(default)]
        container_ids: Option<Vec<ContainerId>>,
        #[serde(default)]
        arrangement_id: Option<String>,
    },
    RenderProtocolCartoonSvg {
        protocol: ProtocolCartoonKind,
        path: String,
    },
    RenderProtocolCartoonTemplateSvg {
        template_path: String,
        path: String,
    },
    ValidateProtocolCartoonTemplate {
        template_path: String,
    },
    RenderProtocolCartoonTemplateWithBindingsSvg {
        template_path: String,
        bindings_path: String,
        path: String,
    },
    ExportProtocolCartoonTemplateJson {
        protocol: ProtocolCartoonKind,
        path: String,
    },
    ApplyGibsonAssemblyPlan {
        plan_json: String,
    },
    CreateArrangementSerial {
        container_ids: Vec<ContainerId>,
        arrangement_id: Option<String>,
        name: Option<String>,
        #[serde(default)]
        ladders: Option<Vec<String>>,
    },
    SetArrangementLadders {
        arrangement_id: String,
        #[serde(default)]
        ladders: Option<Vec<String>>,
    },
    ExportDnaLadders {
        path: String,
        #[serde(default)]
        name_filter: Option<String>,
    },
    ExportRnaLadders {
        path: String,
        #[serde(default)]
        name_filter: Option<String>,
    },
    ExportPool {
        inputs: Vec<SeqId>,
        path: String,
        pool_id: Option<String>,
        human_id: Option<String>,
    },
    ExportProcessRunBundle {
        path: String,
        #[serde(default)]
        run_id: Option<RunId>,
    },
    PrepareGenome {
        genome_id: String,
        catalog_path: Option<String>,
        cache_dir: Option<String>,
        #[serde(default)]
        timeout_seconds: Option<u64>,
    },
    ExtractGenomeRegion {
        genome_id: String,
        chromosome: String,
        start_1based: usize,
        end_1based: usize,
        output_id: Option<SeqId>,
        #[serde(default, skip_serializing_if = "Option::is_none")]
        annotation_scope: Option<GenomeAnnotationScope>,
        #[serde(default, skip_serializing_if = "Option::is_none")]
        max_annotation_features: Option<usize>,
        /// Legacy compatibility flag; prefer `annotation_scope`.
        #[serde(default, skip_serializing_if = "Option::is_none")]
        include_genomic_annotation: Option<bool>,
        catalog_path: Option<String>,
        cache_dir: Option<String>,
    },
    ExtractGenomeGene {
        genome_id: String,
        gene_query: String,
        occurrence: Option<usize>,
        output_id: Option<SeqId>,
        #[serde(default, skip_serializing_if = "Option::is_none")]
        extract_mode: Option<GenomeGeneExtractMode>,
        #[serde(default, skip_serializing_if = "Option::is_none")]
        promoter_upstream_bp: Option<usize>,
        #[serde(default, skip_serializing_if = "Option::is_none")]
        annotation_scope: Option<GenomeAnnotationScope>,
        #[serde(default, skip_serializing_if = "Option::is_none")]
        max_annotation_features: Option<usize>,
        /// Legacy compatibility flag; prefer `annotation_scope`.
        #[serde(default, skip_serializing_if = "Option::is_none")]
        include_genomic_annotation: Option<bool>,
        catalog_path: Option<String>,
        cache_dir: Option<String>,
    },
    ExtendGenomeAnchor {
        seq_id: SeqId,
        side: GenomeAnchorSide,
        length_bp: usize,
        output_id: Option<SeqId>,
        catalog_path: Option<String>,
        cache_dir: Option<String>,
        #[serde(default)]
        prepared_genome_id: Option<String>,
    },
    VerifyGenomeAnchor {
        seq_id: SeqId,
        catalog_path: Option<String>,
        cache_dir: Option<String>,
        #[serde(default)]
        prepared_genome_id: Option<String>,
    },
    ImportGenomeBedTrack {
        seq_id: SeqId,
        path: String,
        track_name: Option<String>,
        min_score: Option<f64>,
        max_score: Option<f64>,
        clear_existing: Option<bool>,
    },
    ImportGenomeBigWigTrack {
        seq_id: SeqId,
        path: String,
        track_name: Option<String>,
        min_score: Option<f64>,
        max_score: Option<f64>,
        clear_existing: Option<bool>,
    },
    ImportGenomeVcfTrack {
        seq_id: SeqId,
        path: String,
        track_name: Option<String>,
        min_score: Option<f64>,
        max_score: Option<f64>,
        clear_existing: Option<bool>,
    },
    ImportIsoformPanel {
        seq_id: SeqId,
        panel_path: String,
        panel_id: Option<String>,
        #[serde(default)]
        strict: bool,
    },
    ImportUniprotSwissProt {
        path: String,
        entry_id: Option<String>,
    },
    FetchUniprotSwissProt {
        query: String,
        entry_id: Option<String>,
    },
    FetchGenBankAccession {
        accession: String,
        as_id: Option<SeqId>,
    },
    FetchDbSnpRegion {
        rs_id: String,
        genome_id: String,
        flank_bp: Option<usize>,
        output_id: Option<SeqId>,
        annotation_scope: Option<GenomeAnnotationScope>,
        max_annotation_features: Option<usize>,
        catalog_path: Option<String>,
        cache_dir: Option<String>,
    },
    FetchUniprotLinkedGenBank {
        entry_id: String,
        #[serde(default, skip_serializing_if = "Option::is_none")]
        accession: Option<String>,
        #[serde(default, skip_serializing_if = "Option::is_none")]
        as_id: Option<SeqId>,
    },
    ImportUniprotEntrySequence {
        entry_id: String,
        output_id: Option<SeqId>,
    },
    ProjectUniprotToGenome {
        seq_id: SeqId,
        entry_id: String,
        projection_id: Option<String>,
        transcript_id: Option<String>,
    },
    ImportBlastHitsTrack {
        seq_id: SeqId,
        hits: Vec<BlastHitFeatureInput>,
        track_name: Option<String>,
        clear_existing: Option<bool>,
        #[serde(default)]
        blast_provenance: Option<BlastInvocationProvenance>,
    },
    DigestContainer {
        container_id: ContainerId,
        enzymes: Vec<String>,
        output_prefix: Option<String>,
    },
    MergeContainersById {
        container_ids: Vec<ContainerId>,
        output_prefix: Option<String>,
    },
    LigationContainer {
        container_id: ContainerId,
        circularize_if_possible: bool,
        output_id: Option<SeqId>,
        protocol: LigationProtocol,
        output_prefix: Option<String>,
        unique: Option<bool>,
    },
    FilterContainerByMolecularWeight {
        container_id: ContainerId,
        min_bp: usize,
        max_bp: usize,
        error: f64,
        unique: bool,
        output_prefix: Option<String>,
    },
    Digest {
        input: SeqId,
        enzymes: Vec<String>,
        output_prefix: Option<String>,
    },
    Ligation {
        inputs: Vec<SeqId>,
        circularize_if_possible: bool,
        output_id: Option<SeqId>,
        protocol: LigationProtocol,
        output_prefix: Option<String>,
        unique: Option<bool>,
    },
    MergeContainers {
        inputs: Vec<SeqId>,
        output_prefix: Option<String>,
    },
    Pcr {
        template: SeqId,
        forward_primer: String,
        reverse_primer: String,
        output_id: Option<SeqId>,
        unique: Option<bool>,
    },
    PcrAdvanced {
        template: SeqId,
        forward_primer: PcrPrimerSpec,
        reverse_primer: PcrPrimerSpec,
        output_id: Option<SeqId>,
        unique: Option<bool>,
    },
    PcrMutagenesis {
        template: SeqId,
        forward_primer: PcrPrimerSpec,
        reverse_primer: PcrPrimerSpec,
        mutations: Vec<SnpMutationSpec>,
        output_id: Option<SeqId>,
        unique: Option<bool>,
        require_all_mutations: Option<bool>,
    },
    DesignPrimerPairs {
        template: SeqId,
        roi_start_0based: usize,
        roi_end_0based: usize,
        #[serde(default)]
        forward: PrimerDesignSideConstraint,
        #[serde(default)]
        reverse: PrimerDesignSideConstraint,
        #[serde(default)]
        pair_constraints: PrimerDesignPairConstraint,
        min_amplicon_bp: usize,
        max_amplicon_bp: usize,
        max_tm_delta_c: Option<f64>,
        max_pairs: Option<usize>,
        report_id: Option<String>,
    },
    DesignInsertionPrimerPairs {
        template: SeqId,
        insertion: PrimerInsertionIntent,
        #[serde(default)]
        forward: PrimerDesignSideConstraint,
        #[serde(default)]
        reverse: PrimerDesignSideConstraint,
        #[serde(default)]
        pair_constraints: PrimerDesignPairConstraint,
        min_amplicon_bp: usize,
        max_amplicon_bp: usize,
        max_tm_delta_c: Option<f64>,
        max_pairs: Option<usize>,
        report_id: Option<String>,
    },
    PcrOverlapExtensionMutagenesis {
        template: SeqId,
        edit_start_0based: usize,
        edit_end_0based_exclusive: usize,
        #[serde(default)]
        insert_sequence: String,
        #[serde(default)]
        constraints: OverlapExtensionMutagenesisConstraints,
        output_prefix: Option<String>,
    },
    DesignQpcrAssays {
        template: SeqId,
        roi_start_0based: usize,
        roi_end_0based: usize,
        #[serde(default)]
        forward: PrimerDesignSideConstraint,
        #[serde(default)]
        reverse: PrimerDesignSideConstraint,
        #[serde(default)]
        probe: PrimerDesignSideConstraint,
        #[serde(default)]
        pair_constraints: PrimerDesignPairConstraint,
        min_amplicon_bp: usize,
        max_amplicon_bp: usize,
        max_tm_delta_c: Option<f64>,
        max_probe_tm_delta_c: Option<f64>,
        max_assays: Option<usize>,
        report_id: Option<String>,
    },
    DeriveTranscriptSequences {
        seq_id: SeqId,
        #[serde(default)]
        feature_ids: Vec<usize>,
        #[serde(default)]
        scope: Option<SplicingScopePreset>,
        #[serde(default)]
        output_prefix: Option<String>,
    },
    ComputeDotplot {
        seq_id: SeqId,
        #[serde(default)]
        reference_seq_id: Option<SeqId>,
        #[serde(default)]
        span_start_0based: Option<usize>,
        #[serde(default)]
        span_end_0based: Option<usize>,
        #[serde(default)]
        reference_span_start_0based: Option<usize>,
        #[serde(default)]
        reference_span_end_0based: Option<usize>,
        #[serde(default)]
        mode: DotplotMode,
        word_size: usize,
        step_bp: usize,
        #[serde(default)]
        max_mismatches: usize,
        #[serde(default)]
        tile_bp: Option<usize>,
        #[serde(default)]
        store_as: Option<String>,
    },
    ComputeDotplotOverlay {
        owner_seq_id: SeqId,
        reference_seq_id: SeqId,
        #[serde(default)]
        reference_span_start_0based: Option<usize>,
        #[serde(default)]
        reference_span_end_0based: Option<usize>,
        #[serde(default)]
        queries: Vec<DotplotOverlayQuerySpec>,
        word_size: usize,
        step_bp: usize,
        #[serde(default)]
        max_mismatches: usize,
        #[serde(default)]
        tile_bp: Option<usize>,
        #[serde(default)]
        store_as: Option<String>,
    },
    ComputeFlexibilityTrack {
        seq_id: SeqId,
        #[serde(default)]
        span_start_0based: Option<usize>,
        #[serde(default)]
        span_end_0based: Option<usize>,
        #[serde(default)]
        model: FlexibilityModel,
        bin_bp: usize,
        #[serde(default)]
        smoothing_bp: Option<usize>,
        #[serde(default)]
        store_as: Option<String>,
    },
    DeriveSplicingReferences {
        seq_id: SeqId,
        span_start_0based: usize,
        span_end_0based: usize,
        #[serde(default)]
        seed_feature_id: Option<usize>,
        #[serde(default = "default_splicing_reference_scope")]
        scope: SplicingScopePreset,
        #[serde(default)]
        output_prefix: Option<SeqId>,
    },
    AlignSequences {
        query_seq_id: SeqId,
        target_seq_id: SeqId,
        #[serde(default)]
        query_span_start_0based: Option<usize>,
        #[serde(default)]
        query_span_end_0based: Option<usize>,
        #[serde(default)]
        target_span_start_0based: Option<usize>,
        #[serde(default)]
        target_span_end_0based: Option<usize>,
        #[serde(default)]
        mode: PairwiseAlignmentMode,
        #[serde(default = "default_pairwise_match_score")]
        match_score: i32,
        #[serde(default = "default_pairwise_mismatch_score")]
        mismatch_score: i32,
        #[serde(default = "default_pairwise_gap_open")]
        gap_open: i32,
        #[serde(default = "default_pairwise_gap_extend")]
        gap_extend: i32,
    },
    ImportSequencingTrace {
        path: String,
        #[serde(default)]
        trace_id: Option<String>,
        #[serde(default)]
        seq_id: Option<SeqId>,
    },
    ListSequencingTraces {
        #[serde(default)]
        seq_id: Option<SeqId>,
    },
    ShowSequencingTrace {
        trace_id: String,
    },
    ConfirmConstructReads {
        expected_seq_id: SeqId,
        #[serde(default)]
        read_seq_ids: Vec<SeqId>,
        #[serde(default)]
        trace_ids: Vec<String>,
        #[serde(default)]
        targets: Vec<SequencingConfirmationTargetSpec>,
        #[serde(default = "default_sequencing_confirmation_alignment_mode")]
        alignment_mode: PairwiseAlignmentMode,
        #[serde(default = "default_pairwise_match_score")]
        match_score: i32,
        #[serde(default = "default_pairwise_mismatch_score")]
        mismatch_score: i32,
        #[serde(default = "default_pairwise_gap_open")]
        gap_open: i32,
        #[serde(default = "default_pairwise_gap_extend")]
        gap_extend: i32,
        #[serde(default = "default_sequencing_confirmation_min_identity_fraction")]
        min_identity_fraction: f64,
        #[serde(default = "default_sequencing_confirmation_min_target_coverage_fraction")]
        min_target_coverage_fraction: f64,
        #[serde(default = "default_true")]
        allow_reverse_complement: bool,
        #[serde(default)]
        report_id: Option<String>,
    },
    ListSequencingConfirmationReports {
        #[serde(default)]
        expected_seq_id: Option<SeqId>,
    },
    ShowSequencingConfirmationReport {
        report_id: String,
    },
    ExportSequencingConfirmationReport {
        report_id: String,
        path: String,
    },
    ExportSequencingConfirmationSupportTsv {
        report_id: String,
        path: String,
    },
    InterpretRnaReads {
        seq_id: SeqId,
        seed_feature_id: usize,
        #[serde(default)]
        profile: RnaReadInterpretationProfile,
        input_path: String,
        #[serde(default)]
        input_format: RnaReadInputFormat,
        #[serde(default)]
        scope: SplicingScopePreset,
        #[serde(default)]
        origin_mode: RnaReadOriginMode,
        #[serde(default)]
        target_gene_ids: Vec<String>,
        #[serde(default)]
        roi_seed_capture_enabled: bool,
        #[serde(default)]
        seed_filter: RnaReadSeedFilterConfig,
        #[serde(default)]
        align_config: RnaReadAlignConfig,
        #[serde(default)]
        report_id: Option<String>,
        #[serde(default)]
        report_mode: RnaReadReportMode,
        #[serde(default)]
        checkpoint_path: Option<String>,
        #[serde(default = "default_rna_read_checkpoint_every_reads")]
        checkpoint_every_reads: usize,
        #[serde(default)]
        resume_from_checkpoint: bool,
    },
    AlignRnaReadReport {
        report_id: String,
        #[serde(default = "default_rna_align_report_selection")]
        selection: RnaReadHitSelection,
        #[serde(default)]
        align_config_override: Option<RnaReadAlignConfig>,
        #[serde(default, skip_serializing_if = "Vec::is_empty")]
        selected_record_indices: Vec<usize>,
    },
    ListRnaReadReports {
        #[serde(default)]
        seq_id: Option<SeqId>,
    },
    ShowRnaReadReport {
        report_id: String,
    },
    SummarizeRnaReadGeneSupport {
        report_id: String,
        gene_ids: Vec<String>,
        #[serde(default, skip_serializing_if = "Vec::is_empty")]
        selected_record_indices: Vec<usize>,
        #[serde(default)]
        complete_rule: RnaReadGeneSupportCompleteRule,
        #[serde(default, skip_serializing_if = "Option::is_none")]
        path: Option<String>,
    },
    ExportRnaReadReport {
        report_id: String,
        path: String,
    },
    ExportRnaReadHitsFasta {
        report_id: String,
        path: String,
        #[serde(default)]
        selection: RnaReadHitSelection,
        #[serde(default, skip_serializing_if = "Vec::is_empty")]
        selected_record_indices: Vec<usize>,
        #[serde(default, skip_serializing_if = "Option::is_none")]
        subset_spec: Option<String>,
    },
    ExportRnaReadSampleSheet {
        path: String,
        #[serde(default)]
        seq_id: Option<SeqId>,
        #[serde(default)]
        report_ids: Vec<String>,
        #[serde(default)]
        append: bool,
    },
    ExportRnaReadExonPathsTsv {
        report_id: String,
        path: String,
        #[serde(default)]
        selection: RnaReadHitSelection,
        #[serde(default, skip_serializing_if = "Vec::is_empty")]
        selected_record_indices: Vec<usize>,
        #[serde(default, skip_serializing_if = "Option::is_none")]
        subset_spec: Option<String>,
    },
    ExportRnaReadExonAbundanceTsv {
        report_id: String,
        path: String,
        #[serde(default)]
        selection: RnaReadHitSelection,
        #[serde(default, skip_serializing_if = "Vec::is_empty")]
        selected_record_indices: Vec<usize>,
        #[serde(default, skip_serializing_if = "Option::is_none")]
        subset_spec: Option<String>,
    },
    ExportRnaReadScoreDensitySvg {
        report_id: String,
        path: String,
        #[serde(default)]
        scale: RnaReadScoreDensityScale,
        #[serde(default)]
        variant: RnaReadScoreDensityVariant,
        #[serde(default, skip_serializing_if = "Option::is_none")]
        seed_filter_override: Option<RnaReadSeedFilterConfig>,
    },
    ExportRnaReadAlignmentsTsv {
        report_id: String,
        path: String,
        #[serde(default)]
        selection: RnaReadHitSelection,
        #[serde(default)]
        limit: Option<usize>,
        #[serde(default, skip_serializing_if = "Vec::is_empty")]
        selected_record_indices: Vec<usize>,
        #[serde(default, skip_serializing_if = "Option::is_none")]
        subset_spec: Option<String>,
    },
    ExportRnaReadAlignmentDotplotSvg {
        report_id: String,
        path: String,
        #[serde(default)]
        selection: RnaReadHitSelection,
        #[serde(default = "default_rna_read_alignment_dotplot_max_points")]
        max_points: usize,
    },
    MaterializeRnaReadHitSequences {
        report_id: String,
        #[serde(default = "default_rna_align_report_selection")]
        selection: RnaReadHitSelection,
        #[serde(default, skip_serializing_if = "Vec::is_empty")]
        selected_record_indices: Vec<usize>,
        #[serde(default)]
        output_prefix: Option<SeqId>,
    },
    ExtractRegion {
        input: SeqId,
        from: usize,
        to: usize,
        output_id: Option<SeqId>,
    },
    ExtractAnchoredRegion {
        input: SeqId,
        anchor: SequenceAnchor,
        direction: AnchorDirection,
        target_length_bp: usize,
        length_tolerance_bp: usize,
        required_re_sites: Vec<String>,
        required_tf_motifs: Vec<String>,
        forward_primer: Option<String>,
        reverse_primer: Option<String>,
        output_prefix: Option<String>,
        unique: Option<bool>,
        max_candidates: Option<usize>,
    },
    GenerateCandidateSetBetweenAnchors {
        set_name: String,
        seq_id: SeqId,
        anchor_a: SequenceAnchor,
        anchor_b: SequenceAnchor,
        length_bp: usize,
        step_bp: usize,
        limit: Option<usize>,
    },
    SelectCandidate {
        input: SeqId,
        criterion: String,
        output_id: Option<SeqId>,
    },
    FilterByMolecularWeight {
        inputs: Vec<SeqId>,
        min_bp: usize,
        max_bp: usize,
        error: f64,
        unique: bool,
        output_prefix: Option<String>,
    },
    #[serde(alias = "FilterBySequenceQuality")]
    FilterByDesignConstraints {
        inputs: Vec<SeqId>,
        gc_min: Option<f64>,
        gc_max: Option<f64>,
        max_homopolymer_run: Option<usize>,
        reject_ambiguous_bases: Option<bool>,
        avoid_u6_terminator_tttt: Option<bool>,
        #[serde(default)]
        forbidden_motifs: Vec<String>,
        unique: bool,
        output_prefix: Option<String>,
    },
    GenerateCandidateSet {
        set_name: String,
        seq_id: SeqId,
        length_bp: usize,
        step_bp: usize,
        #[serde(default)]
        feature_kinds: Vec<String>,
        feature_label_regex: Option<String>,
        max_distance_bp: Option<usize>,
        #[serde(default)]
        feature_geometry_mode: Option<CandidateFeatureGeometryMode>,
        #[serde(default)]
        feature_boundary_mode: Option<CandidateFeatureBoundaryMode>,
        #[serde(default)]
        feature_strand_relation: Option<CandidateFeatureStrandRelation>,
        limit: Option<usize>,
    },
    DeleteCandidateSet {
        set_name: String,
    },
    UpsertGuideSet {
        guide_set_id: String,
        guides: Vec<GuideCandidate>,
    },
    DeleteGuideSet {
        guide_set_id: String,
    },
    FilterGuidesPractical {
        guide_set_id: String,
        #[serde(default)]
        config: GuidePracticalFilterConfig,
        #[serde(default)]
        output_guide_set_id: Option<String>,
    },
    GenerateGuideOligos {
        guide_set_id: String,
        template_id: String,
        #[serde(default)]
        apply_5prime_g_extension: Option<bool>,
        #[serde(default)]
        output_oligo_set_id: Option<String>,
        #[serde(default)]
        passed_only: Option<bool>,
    },
    ExportGuideOligos {
        guide_set_id: String,
        #[serde(default)]
        oligo_set_id: Option<String>,
        format: GuideOligoExportFormat,
        path: String,
        #[serde(default)]
        plate_format: Option<GuideOligoPlateFormat>,
    },
    ExportGuideProtocolText {
        guide_set_id: String,
        #[serde(default)]
        oligo_set_id: Option<String>,
        path: String,
        #[serde(default)]
        include_qc_checklist: Option<bool>,
    },
    ScoreCandidateSetExpression {
        set_name: String,
        metric: String,
        expression: String,
    },
    ScoreCandidateSetDistance {
        set_name: String,
        metric: String,
        #[serde(default)]
        feature_kinds: Vec<String>,
        feature_label_regex: Option<String>,
        #[serde(default)]
        feature_geometry_mode: Option<CandidateFeatureGeometryMode>,
        #[serde(default)]
        feature_boundary_mode: Option<CandidateFeatureBoundaryMode>,
        #[serde(default)]
        feature_strand_relation: Option<CandidateFeatureStrandRelation>,
    },
    FilterCandidateSet {
        input_set: String,
        output_set: String,
        metric: String,
        min: Option<f64>,
        max: Option<f64>,
        min_quantile: Option<f64>,
        max_quantile: Option<f64>,
    },
    CandidateSetOp {
        op: CandidateSetOperator,
        left_set: String,
        right_set: String,
        output_set: String,
    },
    ScoreCandidateSetWeightedObjective {
        set_name: String,
        metric: String,
        #[serde(default)]
        objectives: Vec<CandidateWeightedObjectiveTerm>,
        #[serde(default)]
        normalize_metrics: Option<bool>,
    },
    TopKCandidateSet {
        input_set: String,
        output_set: String,
        metric: String,
        k: usize,
        #[serde(default)]
        direction: Option<CandidateObjectiveDirection>,
        #[serde(default)]
        tie_break: Option<CandidateTieBreakPolicy>,
    },
    ParetoFrontierCandidateSet {
        input_set: String,
        output_set: String,
        #[serde(default)]
        objectives: Vec<CandidateObjectiveSpec>,
        #[serde(default)]
        max_candidates: Option<usize>,
        #[serde(default)]
        tie_break: Option<CandidateTieBreakPolicy>,
    },
    UpsertWorkflowMacroTemplate {
        name: String,
        description: Option<String>,
        #[serde(default)]
        details_url: Option<String>,
        #[serde(default)]
        parameters: Vec<WorkflowMacroTemplateParam>,
        #[serde(default)]
        input_ports: Vec<WorkflowMacroTemplatePort>,
        #[serde(default)]
        output_ports: Vec<WorkflowMacroTemplatePort>,
        script: String,
    },
    DeleteWorkflowMacroTemplate {
        name: String,
    },
    UpsertCandidateMacroTemplate {
        name: String,
        description: Option<String>,
        #[serde(default)]
        details_url: Option<String>,
        #[serde(default)]
        parameters: Vec<CandidateMacroTemplateParam>,
        script: String,
    },
    DeleteCandidateMacroTemplate {
        name: String,
    },
    Reverse {
        input: SeqId,
        output_id: Option<SeqId>,
    },
    Complement {
        input: SeqId,
        output_id: Option<SeqId>,
    },
    ReverseComplement {
        input: SeqId,
        output_id: Option<SeqId>,
    },
    Branch {
        input: SeqId,
        output_id: Option<SeqId>,
    },
    SetDisplayVisibility {
        target: DisplayTarget,
        visible: bool,
    },
    SetLinearViewport {
        start_bp: usize,
        span_bp: usize,
    },
    SetTopology {
        seq_id: SeqId,
        circular: bool,
    },
    RecomputeFeatures {
        seq_id: SeqId,
    },
    SetParameter {
        name: String,
        value: serde_json::Value,
    },
    AnnotateTfbs {
        seq_id: SeqId,
        motifs: Vec<String>,
        min_llr_bits: Option<f64>,
        min_llr_quantile: Option<f64>,
        #[serde(default)]
        per_tf_thresholds: Vec<TfThresholdOverride>,
        clear_existing: Option<bool>,
        #[serde(default)]
        max_hits: Option<usize>,
    },
}

#[derive(Debug, Clone, Serialize, Deserialize)]
struct PoolEnd {
    end_type: String,
    forward_5: String,
    forward_3: String,
    reverse_5: String,
    reverse_3: String,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
struct PoolMember {
    seq_id: SeqId,
    human_id: String,
    name: Option<String>,
    sequence: String,
    length_bp: usize,
    topology: String,
    ends: PoolEnd,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
struct PoolExport {
    schema: String,
    pool_id: String,
    human_id: String,
    member_count: usize,
    members: Vec<PoolMember>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ProcessRunBundleOperationInputSummary {
    pub op_id: String,
    pub run_id: String,
    pub operation: String,
    pub record_index: usize,
    #[serde(default)]
    pub sequence_ids: Vec<String>,
    #[serde(default)]
    pub container_ids: Vec<String>,
    #[serde(default)]
    pub arrangement_ids: Vec<String>,
    #[serde(default)]
    pub candidate_set_ids: Vec<String>,
    #[serde(default)]
    pub guide_set_ids: Vec<String>,
    #[serde(default)]
    pub genome_ids: Vec<String>,
    #[serde(default)]
    pub file_paths: Vec<String>,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
pub struct ProcessRunBundleInputs {
    pub root_sequence_ids: Vec<String>,
    pub referenced_sequence_ids: Vec<String>,
    pub referenced_container_ids: Vec<String>,
    pub referenced_arrangement_ids: Vec<String>,
    pub referenced_candidate_set_ids: Vec<String>,
    pub referenced_guide_set_ids: Vec<String>,
    pub referenced_genome_ids: Vec<String>,
    pub file_inputs: Vec<String>,
    pub operation_inputs: Vec<ProcessRunBundleOperationInputSummary>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ProcessRunBundleParameterOverride {
    pub op_id: String,
    pub run_id: String,
    pub record_index: usize,
    pub name: String,
    pub value: serde_json::Value,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
pub struct ProcessRunBundleOutputs {
    pub created_seq_ids: Vec<String>,
    pub changed_seq_ids: Vec<String>,
    pub final_sequences: Vec<EngineSequenceSummary>,
    pub created_container_ids: Vec<String>,
    pub created_arrangement_ids: Vec<String>,
    pub exported_paths: Vec<String>,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
pub struct RoutineDecisionTraceComparison {
    pub left_routine_id: String,
    pub right_routine_id: String,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
pub struct RoutineDecisionTraceDisambiguationQuestion {
    pub question_id: String,
    pub question_text: String,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
pub struct RoutineDecisionTraceDisambiguationAnswer {
    pub question_id: String,
    pub answer_text: String,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
pub struct RoutineDecisionTracePreflightSnapshot {
    pub can_execute: bool,
    pub warnings: Vec<String>,
    pub errors: Vec<String>,
    pub contract_source: Option<String>,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
pub struct RoutineDecisionTraceExportEvent {
    pub run_bundle_path: String,
    pub exported_at_unix_ms: u128,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
pub struct RoutineDecisionTrace {
    pub schema: String,
    pub trace_id: String,
    pub source: String,
    pub status: String,
    pub created_at_unix_ms: u128,
    pub updated_at_unix_ms: u128,
    pub goal_text: String,
    pub query_text: String,
    pub candidate_routine_ids: Vec<String>,
    pub selected_routine_id: Option<String>,
    pub selected_routine_title: Option<String>,
    pub selected_routine_family: Option<String>,
    pub alternatives_presented: Vec<String>,
    pub comparisons: Vec<RoutineDecisionTraceComparison>,
    pub disambiguation_questions_presented: Vec<RoutineDecisionTraceDisambiguationQuestion>,
    pub disambiguation_answers: Vec<RoutineDecisionTraceDisambiguationAnswer>,
    pub bindings_snapshot: BTreeMap<String, String>,
    pub preflight_history: Vec<RoutineDecisionTracePreflightSnapshot>,
    pub preflight_snapshot: Option<RoutineDecisionTracePreflightSnapshot>,
    pub execution_attempted: bool,
    pub execution_success: Option<bool>,
    pub transactional: Option<bool>,
    pub macro_instance_id: Option<String>,
    pub emitted_operation_ids: Vec<String>,
    pub execution_error: Option<String>,
    pub export_events: Vec<RoutineDecisionTraceExportEvent>,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
pub struct RoutineDecisionTraceStore {
    pub schema: String,
    pub traces: Vec<RoutineDecisionTrace>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ProcessRunBundleExport {
    pub schema: String,
    pub generated_at_unix_ms: u128,
    #[serde(default)]
    pub run_id_filter: Option<String>,
    pub selected_record_count: usize,
    pub inputs: ProcessRunBundleInputs,
    #[serde(default)]
    pub parameter_overrides: Vec<ProcessRunBundleParameterOverride>,
    #[serde(default)]
    pub decision_traces: Vec<RoutineDecisionTrace>,
    pub operation_log: Vec<OperationRecord>,
    pub outputs: ProcessRunBundleOutputs,
    pub parameter_snapshot: serde_json::Value,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct DnaLadderBandInfo {
    pub length_bp: f64,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub relative_strength: Option<f64>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct DnaLadderInfo {
    pub name: String,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub loading_hint: Option<f64>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub min_bp: Option<usize>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub max_bp: Option<usize>,
    pub band_count: usize,
    pub bands: Vec<DnaLadderBandInfo>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct DnaLadderCatalog {
    pub schema: String,
    pub ladder_count: usize,
    pub ladders: Vec<DnaLadderInfo>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct DnaLadderExportReport {
    pub path: String,
    pub ladder_count: usize,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct RnaLadderBandInfo {
    pub length_nt: f64,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub relative_strength: Option<f64>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct RnaLadderInfo {
    pub name: String,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub loading_hint: Option<f64>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub min_nt: Option<usize>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub max_nt: Option<usize>,
    pub band_count: usize,
    pub bands: Vec<RnaLadderBandInfo>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct RnaLadderCatalog {
    pub schema: String,
    pub ladder_count: usize,
    pub ladders: Vec<RnaLadderInfo>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct RnaLadderExportReport {
    pub path: String,
    pub ladder_count: usize,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct GenomeExtractionProvenance {
    pub seq_id: SeqId,
    pub recorded_at_unix_ms: u128,
    pub operation: String,
    pub genome_id: String,
    pub catalog_path: String,
    pub cache_dir: Option<String>,
    pub chromosome: Option<String>,
    pub start_1based: Option<usize>,
    pub end_1based: Option<usize>,
    pub gene_query: Option<String>,
    pub occurrence: Option<usize>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub gene_extract_mode: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub promoter_upstream_bp: Option<usize>,
    pub gene_id: Option<String>,
    pub gene_name: Option<String>,
    pub strand: Option<char>,
    pub anchor_strand: Option<char>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub anchor_verified: Option<bool>,
    pub sequence_source_type: Option<String>,
    pub annotation_source_type: Option<String>,
    pub sequence_source: Option<String>,
    pub annotation_source: Option<String>,
    pub sequence_sha1: Option<String>,
    pub annotation_sha1: Option<String>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SequenceGenomeAnchorSummary {
    pub seq_id: String,
    pub genome_id: String,
    pub chromosome: String,
    pub start_1based: usize,
    pub end_1based: usize,
    pub strand: Option<char>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub anchor_verified: Option<bool>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SequenceAnchorPreparedGenomeOptionsSummary {
    pub seq_id: String,
    pub requested_genome_id: String,
    pub requested_catalog_key: String,
    pub requested_family: Option<String>,
    pub exact_prepared: bool,
    pub compatible_prepared_options: Vec<String>,
}

#[derive(Debug, Clone)]
struct GenomeSequenceAnchor {
    genome_id: String,
    chromosome: String,
    start_1based: usize,
    end_1based: usize,
    strand: Option<char>,
    anchor_verified: Option<bool>,
    catalog_path: Option<String>,
    cache_dir: Option<String>,
}

#[derive(Debug, Clone, Default)]
struct EngineHistoryCheckpoint {
    state: ProjectState,
    journal: Vec<OperationRecord>,
    op_counter: u64,
}

#[derive(Debug, Clone)]
struct BedRecord {
    chromosome: String,
    start_0based: usize,
    end_0based: usize,
    name: Option<String>,
    score: Option<f64>,
    strand: Option<char>,
}

#[derive(Debug, Clone)]
struct VcfRecord {
    chromosome: String,
    pos_1based: usize,
    id: Option<String>,
    reference: String,
    alternates: Vec<String>,
    qual: Option<f64>,
    filter: Option<String>,
    info: Option<String>,
    format: Option<String>,
    sample_columns: Vec<String>,
}

#[derive(Debug, Clone, Default)]
struct GenomeBedTrackImportReport {
    track_name: String,
    parsed_records: usize,
    imported_features: usize,
    skipped_records: usize,
    skipped_invalid: usize,
    skipped_wrong_chromosome: usize,
    skipped_non_overlap: usize,
    skipped_missing_score: usize,
    skipped_outside_score_range: usize,
    truncated_at_limit: bool,
    cancelled: bool,
    warnings: Vec<String>,
    skipped_wrong_chromosome_examples: Vec<String>,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum VcfVariantClass {
    Snp,
    Ins,
    Del,
    Sv,
    Other,
}

impl VcfVariantClass {
    fn as_str(self) -> &'static str {
        match self {
            Self::Snp => "SNP",
            Self::Ins => "INS",
            Self::Del => "DEL",
            Self::Sv => "SV",
            Self::Other => "OTHER",
        }
    }
}

#[derive(Debug, Clone, Default)]
struct VcfAltGenotypeSummary {
    carriers: usize,
    het: usize,
    hom_alt: usize,
    mixed_alt: usize,
    haploid_alt: usize,
    phased_carriers: usize,
    unphased_carriers: usize,
    carrier_samples: Vec<String>,
}

#[derive(Debug, Clone)]
struct FeatureDistanceTarget {
    feature_index: usize,
    kind_upper: String,
    labels_upper: Vec<String>,
    start_0based: usize,
    end_0based: usize,
    strand: Option<char>,
}

#[derive(Debug, Clone, Copy)]
struct FeatureLocationSegment {
    start_0based: usize,
    end_0based: usize,
    strand: Option<char>,
}

#[derive(Debug, Clone)]
struct PrimerDesignCandidate {
    sequence: String,
    start_0based: usize,
    end_0based_exclusive: usize,
    tm_c: f64,
    gc_fraction: f64,
    anneal_hits: usize,
}

#[derive(Debug, Clone, Copy, Default)]
struct PrimerHeuristicMetrics {
    length_bp: usize,
    three_prime_gc_clamp: bool,
    three_prime_base: u8,
    longest_homopolymer_run_bp: usize,
    self_complementary_run_bp: usize,
}

#[derive(Debug, Clone, Copy, Default)]
struct PrimerPairDimerMetrics {
    max_complementary_run_bp: usize,
    max_3prime_complementary_run_bp: usize,
}

#[derive(Debug, Clone, Default)]
struct NormalizedPrimerSideSequenceConstraints {
    non_annealing_5prime_tail: Option<String>,
    fixed_5prime: Option<String>,
    fixed_3prime: Option<String>,
    required_motifs: Vec<String>,
    forbidden_motifs: Vec<String>,
    locked_positions: Vec<(usize, String)>,
}

#[derive(Debug, Clone, Default)]
struct NormalizedPrimerPairConstraints {
    require_roi_flanking: bool,
    required_amplicon_motifs: Vec<String>,
    forbidden_amplicon_motifs: Vec<String>,
    fixed_amplicon_start_0based: Option<usize>,
    fixed_amplicon_end_0based_exclusive: Option<usize>,
}

#[derive(Debug, Clone, PartialEq)]
enum ExpressionToken {
    Number(f64),
    Ident(String),
    Plus,
    Minus,
    Star,
    Slash,
    LParen,
    RParen,
    Comma,
}

#[derive(Debug, Clone)]
enum MetricExpr {
    Number(f64),
    Variable(String),
    UnaryMinus(Box<MetricExpr>),
    Binary {
        op: ExpressionBinaryOp,
        left: Box<MetricExpr>,
        right: Box<MetricExpr>,
    },
    Function {
        name: String,
        args: Vec<MetricExpr>,
    },
}

#[derive(Debug, Clone, Copy)]
enum ExpressionBinaryOp {
    Add,
    Subtract,
    Multiply,
    Divide,
}

struct MetricExpressionParser {
    tokens: Vec<ExpressionToken>,
    index: usize,
}

impl MetricExpressionParser {
    fn new(tokens: Vec<ExpressionToken>) -> Self {
        Self { tokens, index: 0 }
    }

    fn parse(mut self) -> Result<MetricExpr, String> {
        let expr = self.parse_add_sub()?;
        if self.index != self.tokens.len() {
            return Err("Unexpected trailing expression tokens".to_string());
        }
        Ok(expr)
    }

    fn peek(&self) -> Option<&ExpressionToken> {
        self.tokens.get(self.index)
    }

    fn consume(&mut self) -> Option<ExpressionToken> {
        let token = self.tokens.get(self.index).cloned();
        if token.is_some() {
            self.index += 1;
        }
        token
    }

    fn parse_add_sub(&mut self) -> Result<MetricExpr, String> {
        let mut expr = self.parse_mul_div()?;
        loop {
            let op = match self.peek() {
                Some(ExpressionToken::Plus) => ExpressionBinaryOp::Add,
                Some(ExpressionToken::Minus) => ExpressionBinaryOp::Subtract,
                _ => break,
            };
            let _ = self.consume();
            let rhs = self.parse_mul_div()?;
            expr = MetricExpr::Binary {
                op,
                left: Box::new(expr),
                right: Box::new(rhs),
            };
        }
        Ok(expr)
    }

    fn parse_mul_div(&mut self) -> Result<MetricExpr, String> {
        let mut expr = self.parse_unary()?;
        loop {
            let op = match self.peek() {
                Some(ExpressionToken::Star) => ExpressionBinaryOp::Multiply,
                Some(ExpressionToken::Slash) => ExpressionBinaryOp::Divide,
                _ => break,
            };
            let _ = self.consume();
            let rhs = self.parse_unary()?;
            expr = MetricExpr::Binary {
                op,
                left: Box::new(expr),
                right: Box::new(rhs),
            };
        }
        Ok(expr)
    }

    fn parse_unary(&mut self) -> Result<MetricExpr, String> {
        if matches!(self.peek(), Some(ExpressionToken::Minus)) {
            let _ = self.consume();
            return Ok(MetricExpr::UnaryMinus(Box::new(self.parse_unary()?)));
        }
        self.parse_primary()
    }

    fn parse_primary(&mut self) -> Result<MetricExpr, String> {
        match self.consume() {
            Some(ExpressionToken::Number(value)) => Ok(MetricExpr::Number(value)),
            Some(ExpressionToken::Ident(name)) => {
                if matches!(self.peek(), Some(ExpressionToken::LParen)) {
                    let _ = self.consume();
                    let mut args = vec![];
                    if !matches!(self.peek(), Some(ExpressionToken::RParen)) {
                        loop {
                            args.push(self.parse_add_sub()?);
                            if matches!(self.peek(), Some(ExpressionToken::Comma)) {
                                let _ = self.consume();
                                continue;
                            }
                            break;
                        }
                    }
                    if !matches!(self.consume(), Some(ExpressionToken::RParen)) {
                        return Err(format!("Function '{}' is missing closing ')'", name));
                    }
                    Ok(MetricExpr::Function { name, args })
                } else {
                    Ok(MetricExpr::Variable(name))
                }
            }
            Some(ExpressionToken::LParen) => {
                let expr = self.parse_add_sub()?;
                if !matches!(self.consume(), Some(ExpressionToken::RParen)) {
                    return Err("Missing ')' in expression".to_string());
                }
                Ok(expr)
            }
            _ => Err("Unexpected token while parsing expression".to_string()),
        }
    }
}

#[derive(Debug, Clone, Default)]
struct ExtractRegionAnnotationProjectionBatch {
    features: Vec<gb_io::seq::Feature>,
    gene_count: usize,
    transcript_count: usize,
    exon_count: usize,
    cds_count: usize,
}

impl ExtractRegionAnnotationProjectionBatch {
    fn feature_count(&self) -> usize {
        self.features.len()
    }
}

#[derive(Debug, Clone)]
struct ExonConcatenatedBlock {
    genomic_start_1based: usize,
    genomic_end_1based: usize,
    local_start_0based: usize,
    local_end_0based_exclusive: usize,
}

#[derive(Debug, Clone, Default)]
struct ExonConcatenatedProjection {
    sequence: String,
    blocks: Vec<ExonConcatenatedBlock>,
}

/// Minimal execution contract shared by concrete engine implementations.
///
/// Adapters should prefer this trait boundary when they only need to submit
/// operations/workflows and inspect the resulting snapshot.
pub trait Engine {
    /// Apply one operation and append it to the deterministic operation log.
    fn apply(&mut self, op: Operation) -> Result<OpResult, EngineError>;
    /// Apply a workflow in order using the workflow's caller-supplied `run_id`.
    fn apply_workflow(&mut self, wf: Workflow) -> Result<Vec<OpResult>, EngineError>;
    /// Borrow the current canonical project snapshot.
    fn snapshot(&self) -> &ProjectState;
}

/// Default in-process engine implementation used by GUI, CLI, and scripting
/// adapters.
///
/// The engine owns canonical project state plus a deterministic operation
/// journal and undo/redo history. Helper modules under `src/engine/` extend
/// this type rather than introducing parallel execution layers.
#[derive(Debug, Clone, Default, Serialize, Deserialize)]
pub struct GentleEngine {
    state: ProjectState,
    journal: Vec<OperationRecord>,
    op_counter: u64,
    #[serde(skip, default)]
    undo_stack: Vec<EngineHistoryCheckpoint>,
    #[serde(skip, default)]
    redo_stack: Vec<EngineHistoryCheckpoint>,
    #[serde(skip, default = "GentleEngine::default_history_limit")]
    history_limit: usize,
}

impl GentleEngine {
    fn default_history_limit() -> usize {
        256
    }

    /// Construct a new empty engine with default history settings.
    pub fn new() -> Self {
        let mut ret = Self::default();
        if ret.history_limit == 0 {
            ret.history_limit = Self::default_history_limit();
        }
        ret
    }

    /// Construct an engine from persisted state and reconcile derived indexes.
    ///
    /// This is the preferred rehydration path when loading project JSON because
    /// it restores lineage/container helper structures that may lag behind older
    /// on-disk snapshots.
    pub fn from_state(state: ProjectState) -> Self {
        let mut ret = Self {
            state,
            ..Self::default()
        };
        if ret.history_limit == 0 {
            ret.history_limit = Self::default_history_limit();
        }
        ret.reconcile_lineage_nodes();
        ret.reconcile_containers();
        ret.reseed_op_counter_from_state();
        ret
    }

    fn parse_op_counter_from_id(op_id: &str) -> Option<u64> {
        let trimmed = op_id.trim();
        let (prefix, suffix) = trimmed.split_once('-')?;
        if !prefix.eq_ignore_ascii_case("op") {
            return None;
        }
        suffix.parse::<u64>().ok()
    }

    fn reseed_op_counter_from_state(&mut self) {
        let mut max_seen = self.op_counter;
        let mut consider = |raw: &str| {
            if let Some(value) = Self::parse_op_counter_from_id(raw) {
                max_seen = max_seen.max(value);
            }
        };

        for node in self.state.lineage.nodes.values() {
            if let Some(op_id) = node.created_by_op.as_deref() {
                consider(op_id);
            }
        }
        for edge in &self.state.lineage.edges {
            consider(&edge.op_id);
        }
        for instance in &self.state.lineage.macro_instances {
            for op_id in &instance.expanded_op_ids {
                consider(op_id);
            }
        }
        for container in self.state.container_state.containers.values() {
            if let Some(op_id) = container.created_by_op.as_deref() {
                consider(op_id);
            }
        }
        for arrangement in self.state.container_state.arrangements.values() {
            if let Some(op_id) = arrangement.created_by_op.as_deref() {
                consider(op_id);
            }
        }
        self.op_counter = max_seen;
    }

    /// Borrow the canonical mutable-independent project snapshot.
    pub fn state(&self) -> &ProjectState {
        &self.state
    }

    /// Mutably borrow the canonical project snapshot.
    ///
    /// Direct mutation is intended for tightly controlled internal call sites;
    /// adapter code should normally prefer `apply`/`apply_workflow` so lineage,
    /// journaling, and parity guarantees remain intact.
    pub fn state_mut(&mut self) -> &mut ProjectState {
        &mut self.state
    }

    pub fn list_sequences_with_genome_anchor(&self) -> Vec<String> {
        let mut seq_ids: Vec<String> = self.state.sequences.keys().cloned().collect();
        seq_ids.sort_unstable();
        seq_ids
            .into_iter()
            .filter(|seq_id| self.latest_genome_anchor_for_seq(seq_id).is_ok())
            .collect()
    }

    pub fn describe_sequence_genome_anchor(&self, seq_id: &str) -> Result<String, EngineError> {
        let anchor = self.latest_genome_anchor_for_seq(seq_id)?;
        let strand = anchor.strand.unwrap_or('+');
        let verification = match anchor.anchor_verified {
            Some(true) => "verified",
            Some(false) => "unverified",
            None => "verification n/a",
        };
        Ok(format!(
            "{}:{}-{} ({}, strand {}, {})",
            anchor.chromosome,
            anchor.start_1based,
            anchor.end_1based,
            anchor.genome_id,
            strand,
            verification
        ))
    }

    pub fn sequence_genome_anchor_summary(
        &self,
        seq_id: &str,
    ) -> Result<SequenceGenomeAnchorSummary, EngineError> {
        let anchor = self.latest_genome_anchor_for_seq(seq_id)?;
        Ok(SequenceGenomeAnchorSummary {
            seq_id: seq_id.to_string(),
            genome_id: anchor.genome_id,
            chromosome: anchor.chromosome,
            start_1based: anchor.start_1based,
            end_1based: anchor.end_1based,
            strand: anchor.strand,
            anchor_verified: anchor.anchor_verified,
        })
    }

    pub fn sequence_anchor_prepared_genome_options(
        &self,
        seq_id: &str,
        catalog_path: Option<&str>,
        cache_dir: Option<&str>,
    ) -> Result<SequenceAnchorPreparedGenomeOptionsSummary, EngineError> {
        let anchor = self.latest_genome_anchor_for_seq(seq_id)?;
        let requested_catalog_path = catalog_path
            .map(str::trim)
            .filter(|v| !v.is_empty())
            .map(|v| v.to_string())
            .or_else(|| anchor.catalog_path.clone())
            .unwrap_or_else(|| DEFAULT_GENOME_CATALOG_PATH.to_string());
        let resolved_cache_dir = cache_dir
            .map(str::trim)
            .filter(|v| !v.is_empty())
            .map(|v| v.to_string())
            .or_else(|| anchor.cache_dir.clone());
        let catalog =
            GenomeCatalog::from_json_file(&requested_catalog_path).map_err(|e| EngineError {
                code: ErrorCode::InvalidInput,
                message: format!(
                    "Could not open genome catalog '{}': {}",
                    requested_catalog_path, e
                ),
            })?;
        let inspection = catalog
            .inspect_prepared_genome_compatibility(&anchor.genome_id, resolved_cache_dir.as_deref())
            .map_err(|e| EngineError {
                code: ErrorCode::InvalidInput,
                message: format!(
                    "Could not inspect prepared genome compatibility for '{}': {}",
                    anchor.genome_id, e
                ),
            })?;
        Ok(SequenceAnchorPreparedGenomeOptionsSummary {
            seq_id: seq_id.to_string(),
            requested_genome_id: anchor.genome_id,
            requested_catalog_key: inspection.requested_catalog_key,
            requested_family: inspection.requested_family,
            exact_prepared: inspection.exact_prepared,
            compatible_prepared_options: inspection.compatible_prepared_options,
        })
    }

    pub fn list_sequence_genome_anchor_summaries(&self) -> Vec<SequenceGenomeAnchorSummary> {
        let mut seq_ids: Vec<String> = self.state.sequences.keys().cloned().collect();
        seq_ids.sort_unstable();
        seq_ids
            .into_iter()
            .filter_map(|seq_id| self.sequence_genome_anchor_summary(&seq_id).ok())
            .collect()
    }

    pub fn list_genome_track_subscriptions(&self) -> Vec<GenomeTrackSubscription> {
        Self::read_track_subscriptions_from_metadata(
            self.state
                .metadata
                .get(GENOME_TRACK_SUBSCRIPTIONS_METADATA_KEY),
        )
    }

    pub fn add_genome_track_subscription(
        &mut self,
        subscription: GenomeTrackSubscription,
    ) -> Result<bool, EngineError> {
        let normalized = Self::normalize_genome_track_subscription(subscription)?;
        let mut subscriptions = self.list_genome_track_subscriptions();
        if subscriptions.iter().any(|existing| existing == &normalized) {
            return Ok(false);
        }
        subscriptions.push(normalized);
        Self::sort_track_subscriptions(&mut subscriptions);
        self.write_track_subscriptions_to_metadata(&subscriptions)?;
        Ok(true)
    }

    pub fn remove_genome_track_subscription(
        &mut self,
        index: usize,
    ) -> Result<GenomeTrackSubscription, EngineError> {
        let mut subscriptions = self.list_genome_track_subscriptions();
        if index >= subscriptions.len() {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: format!(
                    "Track subscription index {} is out of bounds (len={})",
                    index,
                    subscriptions.len()
                ),
            });
        }
        let removed = subscriptions.remove(index);
        self.write_track_subscriptions_to_metadata(&subscriptions)?;
        Ok(removed)
    }

    pub fn clear_genome_track_subscriptions(&mut self) {
        self.state
            .metadata
            .remove(GENOME_TRACK_SUBSCRIPTIONS_METADATA_KEY);
        self.state
            .metadata
            .remove(GENOME_TRACK_KNOWN_ANCHORS_METADATA_KEY);
    }

    pub fn import_genome_track_to_all_anchored(
        &mut self,
        subscription: GenomeTrackSubscription,
        track_subscription: bool,
    ) -> Result<GenomeTrackSyncReport, EngineError> {
        let normalized = Self::normalize_genome_track_subscription(subscription)?;
        let seq_ids = self.list_sequences_with_genome_anchor();
        if seq_ids.is_empty() {
            return Err(EngineError {
                code: ErrorCode::NotFound,
                message:
                    "No genome-anchored sequence is available. Extract a genome region or gene first."
                        .to_string(),
            });
        }
        let mut report = self.apply_track_subscription_to_seq_ids(&seq_ids, &normalized);
        report.subscriptions_considered = 1;
        report.target_sequences = seq_ids.len();
        if track_subscription {
            let _ = self.add_genome_track_subscription(normalized)?;
            let known: BTreeSet<String> = seq_ids.iter().cloned().collect();
            self.write_known_track_anchor_ids(&known);
        }
        Ok(report)
    }

    pub fn apply_tracked_genome_track_subscription(
        &mut self,
        index: usize,
    ) -> Result<GenomeTrackSyncReport, EngineError> {
        let subscriptions = self.list_genome_track_subscriptions();
        let Some(subscription) = subscriptions.get(index).cloned() else {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: format!(
                    "Track subscription index {} is out of bounds (len={})",
                    index,
                    subscriptions.len()
                ),
            });
        };
        let seq_ids = self.list_sequences_with_genome_anchor();
        if seq_ids.is_empty() {
            return Err(EngineError {
                code: ErrorCode::NotFound,
                message:
                    "No genome-anchored sequence is available. Extract a genome region or gene first."
                        .to_string(),
            });
        }
        let mut report = self.apply_track_subscription_to_seq_ids(&seq_ids, &subscription);
        report.subscriptions_considered = 1;
        report.target_sequences = seq_ids.len();
        let known: BTreeSet<String> = seq_ids.iter().cloned().collect();
        self.write_known_track_anchor_ids(&known);
        Ok(report)
    }

    pub fn sync_tracked_genome_track_subscriptions(
        &mut self,
        only_new_anchors: bool,
    ) -> Result<GenomeTrackSyncReport, EngineError> {
        let subscriptions = self.list_genome_track_subscriptions();
        let current_anchors = self.list_sequences_with_genome_anchor();
        let current_set: BTreeSet<String> = current_anchors.iter().cloned().collect();
        let target_seq_ids: Vec<String> = if only_new_anchors {
            let known = self.read_known_track_anchor_ids();
            current_set.difference(&known).cloned().collect()
        } else {
            current_anchors.clone()
        };

        if subscriptions.is_empty() || target_seq_ids.is_empty() {
            self.write_known_track_anchor_ids(&current_set);
            return Ok(GenomeTrackSyncReport {
                subscriptions_considered: subscriptions.len(),
                target_sequences: target_seq_ids.len(),
                ..GenomeTrackSyncReport::default()
            });
        }

        let mut report = GenomeTrackSyncReport {
            subscriptions_considered: subscriptions.len(),
            target_sequences: target_seq_ids.len(),
            ..GenomeTrackSyncReport::default()
        };
        for seq_id in &target_seq_ids {
            for subscription in &subscriptions {
                let op = Self::track_subscription_to_operation(seq_id, subscription);
                match self.apply(op) {
                    Ok(op_result) => {
                        report.applied_imports += 1;
                        report.warnings_count += op_result.warnings.len();
                    }
                    Err(e) => {
                        report.failed_imports += 1;
                        if report.errors.len() < 20 {
                            report.errors.push(format!(
                                "{} '{}' @ {}: {}",
                                subscription.source.label(),
                                subscription.path,
                                seq_id,
                                e.message
                            ));
                        }
                    }
                }
            }
        }
        self.write_known_track_anchor_ids(&current_set);
        Ok(report)
    }

    pub fn capabilities() -> Capabilities {
        Capabilities {
            protocol_version: "v1".to_string(),
            supported_operations: vec![
                "LoadFile".to_string(),
                "SaveFile".to_string(),
                "RenderSequenceSvg".to_string(),
                "RenderDotplotSvg".to_string(),
                "RenderFeatureExpertSvg".to_string(),
                "RenderIsoformArchitectureSvg".to_string(),
                "RenderRnaStructureSvg".to_string(),
                "RenderLineageSvg".to_string(),
                "RenderPoolGelSvg".to_string(),
                "RenderProtocolCartoonSvg".to_string(),
                "RenderProtocolCartoonTemplateSvg".to_string(),
                "ValidateProtocolCartoonTemplate".to_string(),
                "RenderProtocolCartoonTemplateWithBindingsSvg".to_string(),
                "ExportProtocolCartoonTemplateJson".to_string(),
                "ApplyGibsonAssemblyPlan".to_string(),
                "CreateArrangementSerial".to_string(),
                "SetArrangementLadders".to_string(),
                "ExportDnaLadders".to_string(),
                "ExportRnaLadders".to_string(),
                "ExportPool".to_string(),
                "ExportProcessRunBundle".to_string(),
                "PrepareGenome".to_string(),
                "ExtractGenomeRegion".to_string(),
                "ExtractGenomeGene".to_string(),
                "ExtendGenomeAnchor".to_string(),
                "ImportGenomeBedTrack".to_string(),
                "ImportGenomeBigWigTrack".to_string(),
                "ImportGenomeVcfTrack".to_string(),
                "ImportIsoformPanel".to_string(),
                "ImportUniprotSwissProt".to_string(),
                "FetchUniprotSwissProt".to_string(),
                "FetchGenBankAccession".to_string(),
                "FetchDbSnpRegion".to_string(),
                "FetchUniprotLinkedGenBank".to_string(),
                "ImportUniprotEntrySequence".to_string(),
                "ProjectUniprotToGenome".to_string(),
                "ImportBlastHitsTrack".to_string(),
                "DigestContainer".to_string(),
                "MergeContainersById".to_string(),
                "LigationContainer".to_string(),
                "FilterContainerByMolecularWeight".to_string(),
                "Digest".to_string(),
                "MergeContainers".to_string(),
                "Ligation".to_string(),
                "Pcr".to_string(),
                "PcrAdvanced".to_string(),
                "PcrMutagenesis".to_string(),
                "DesignPrimerPairs".to_string(),
                "DesignInsertionPrimerPairs".to_string(),
                "PcrOverlapExtensionMutagenesis".to_string(),
                "DesignQpcrAssays".to_string(),
                "DeriveTranscriptSequences".to_string(),
                "ComputeDotplot".to_string(),
                "ComputeDotplotOverlay".to_string(),
                "ComputeFlexibilityTrack".to_string(),
                "DeriveSplicingReferences".to_string(),
                "AlignSequences".to_string(),
                "ConfirmConstructReads".to_string(),
                "ListSequencingConfirmationReports".to_string(),
                "ShowSequencingConfirmationReport".to_string(),
                "ExportSequencingConfirmationReport".to_string(),
                "ExportSequencingConfirmationSupportTsv".to_string(),
                "InterpretRnaReads".to_string(),
                "AlignRnaReadReport".to_string(),
                "ListRnaReadReports".to_string(),
                "ShowRnaReadReport".to_string(),
                "SummarizeRnaReadGeneSupport".to_string(),
                "ExportRnaReadReport".to_string(),
                "ExportRnaReadHitsFasta".to_string(),
                "ExportRnaReadSampleSheet".to_string(),
                "ExportRnaReadExonPathsTsv".to_string(),
                "ExportRnaReadExonAbundanceTsv".to_string(),
                "ExportRnaReadScoreDensitySvg".to_string(),
                "ExportRnaReadAlignmentsTsv".to_string(),
                "ExportRnaReadAlignmentDotplotSvg".to_string(),
                "MaterializeRnaReadHitSequences".to_string(),
                "ExtractRegion".to_string(),
                "ExtractAnchoredRegion".to_string(),
                "SelectCandidate".to_string(),
                "FilterByMolecularWeight".to_string(),
                "FilterByDesignConstraints".to_string(),
                "GenerateCandidateSet".to_string(),
                "GenerateCandidateSetBetweenAnchors".to_string(),
                "DeleteCandidateSet".to_string(),
                "UpsertGuideSet".to_string(),
                "DeleteGuideSet".to_string(),
                "FilterGuidesPractical".to_string(),
                "GenerateGuideOligos".to_string(),
                "ExportGuideOligos".to_string(),
                "ExportGuideProtocolText".to_string(),
                "ScoreCandidateSetExpression".to_string(),
                "ScoreCandidateSetDistance".to_string(),
                "FilterCandidateSet".to_string(),
                "CandidateSetOp".to_string(),
                "ScoreCandidateSetWeightedObjective".to_string(),
                "TopKCandidateSet".to_string(),
                "ParetoFrontierCandidateSet".to_string(),
                "UpsertWorkflowMacroTemplate".to_string(),
                "DeleteWorkflowMacroTemplate".to_string(),
                "UpsertCandidateMacroTemplate".to_string(),
                "DeleteCandidateMacroTemplate".to_string(),
                "Reverse".to_string(),
                "Complement".to_string(),
                "ReverseComplement".to_string(),
                "Branch".to_string(),
                "SetDisplayVisibility".to_string(),
                "SetLinearViewport".to_string(),
                "SetTopology".to_string(),
                "RecomputeFeatures".to_string(),
                "SetParameter".to_string(),
                "AnnotateTfbs".to_string(),
            ],
            supported_export_formats: vec!["GenBank".to_string(), "Fasta".to_string()],
            deterministic_operation_log: true,
        }
    }

    pub fn inspect_dna_ladders(name_filter: Option<&str>) -> DnaLadderCatalog {
        let filter = name_filter
            .map(str::trim)
            .filter(|v| !v.is_empty())
            .map(|v| v.to_ascii_lowercase());

        let mut ladders: Vec<DnaLadderInfo> = vec![];
        for name in DNA_LADDERS.names_sorted() {
            if let Some(filter_text) = &filter {
                if !name.to_ascii_lowercase().contains(filter_text) {
                    continue;
                }
            }
            let Some(ladder) = DNA_LADDERS.get(&name) else {
                continue;
            };
            let bands = ladder
                .bands()
                .iter()
                .map(|band| DnaLadderBandInfo {
                    length_bp: band.length_bp(),
                    relative_strength: band.relative_strength,
                })
                .collect::<Vec<_>>();
            ladders.push(DnaLadderInfo {
                name,
                loading_hint: ladder.loading_hint(),
                min_bp: ladder.min_bp(),
                max_bp: ladder.max_bp(),
                band_count: bands.len(),
                bands,
            });
        }

        DnaLadderCatalog {
            schema: "gentle.dna_ladders.v1".to_string(),
            ladder_count: ladders.len(),
            ladders,
        }
    }

    pub fn export_dna_ladders(
        path: &str,
        name_filter: Option<&str>,
    ) -> Result<DnaLadderExportReport, EngineError> {
        let catalog = Self::inspect_dna_ladders(name_filter);
        let text = serde_json::to_string_pretty(&catalog).map_err(|e| EngineError {
            code: ErrorCode::Internal,
            message: format!("Could not serialize DNA ladders JSON: {e}"),
        })?;
        std::fs::write(path, text).map_err(|e| EngineError {
            code: ErrorCode::Io,
            message: format!("Could not write DNA ladders file '{path}': {e}"),
        })?;
        Ok(DnaLadderExportReport {
            path: path.to_string(),
            ladder_count: catalog.ladder_count,
        })
    }

    pub fn inspect_rna_ladders(name_filter: Option<&str>) -> RnaLadderCatalog {
        let filter = name_filter
            .map(str::trim)
            .filter(|v| !v.is_empty())
            .map(|v| v.to_ascii_lowercase());

        let mut ladders: Vec<RnaLadderInfo> = vec![];
        for name in RNA_LADDERS.names_sorted() {
            if let Some(filter_text) = &filter {
                if !name.to_ascii_lowercase().contains(filter_text) {
                    continue;
                }
            }
            let Some(ladder) = RNA_LADDERS.get(&name) else {
                continue;
            };
            let bands = ladder
                .bands()
                .iter()
                .map(|band| RnaLadderBandInfo {
                    length_nt: band.length_nt(),
                    relative_strength: band.relative_strength,
                })
                .collect::<Vec<_>>();
            ladders.push(RnaLadderInfo {
                name,
                loading_hint: ladder.loading_hint(),
                min_nt: ladder.min_nt(),
                max_nt: ladder.max_nt(),
                band_count: bands.len(),
                bands,
            });
        }

        RnaLadderCatalog {
            schema: "gentle.rna_ladders.v1".to_string(),
            ladder_count: ladders.len(),
            ladders,
        }
    }

    pub fn export_rna_ladders(
        path: &str,
        name_filter: Option<&str>,
    ) -> Result<RnaLadderExportReport, EngineError> {
        let catalog = Self::inspect_rna_ladders(name_filter);
        let text = serde_json::to_string_pretty(&catalog).map_err(|e| EngineError {
            code: ErrorCode::Internal,
            message: format!("Could not serialize RNA ladders JSON: {e}"),
        })?;
        std::fs::write(path, text).map_err(|e| EngineError {
            code: ErrorCode::Io,
            message: format!("Could not write RNA ladders file '{path}': {e}"),
        })?;
        Ok(RnaLadderExportReport {
            path: path.to_string(),
            ladder_count: catalog.ladder_count,
        })
    }

    fn map_rna_structure_error(err: RnaStructureError) -> EngineError {
        match err {
            RnaStructureError::UnsupportedBiotype { .. } | RnaStructureError::EmptySequence => {
                EngineError {
                    code: ErrorCode::InvalidInput,
                    message: err.to_string(),
                }
            }
            RnaStructureError::ToolNotFound { .. } => EngineError {
                code: ErrorCode::Unsupported,
                message: err.to_string(),
            },
            RnaStructureError::ToolFailed { .. } => EngineError {
                code: ErrorCode::Internal,
                message: err.to_string(),
            },
            RnaStructureError::Io { .. } => EngineError {
                code: ErrorCode::Io,
                message: err.to_string(),
            },
        }
    }

    pub fn inspect_rna_structure(
        &self,
        seq_id: &str,
    ) -> Result<RnaStructureTextReport, EngineError> {
        let dna = self
            .state
            .sequences
            .get(seq_id)
            .ok_or_else(|| EngineError {
                code: ErrorCode::NotFound,
                message: format!("Sequence '{seq_id}' not found"),
            })?;
        rna_structure::inspect_text(dna).map_err(Self::map_rna_structure_error)
    }

    pub fn render_rna_structure_svg_to_path(
        &self,
        seq_id: &str,
        path: &str,
    ) -> Result<RnaStructureSvgReport, EngineError> {
        let dna = self
            .state
            .sequences
            .get(seq_id)
            .ok_or_else(|| EngineError {
                code: ErrorCode::NotFound,
                message: format!("Sequence '{seq_id}' not found"),
            })?;
        rna_structure::render_svg(dna, path).map_err(Self::map_rna_structure_error)
    }

    fn open_reference_genome_catalog(
        catalog_path: Option<&str>,
    ) -> Result<(GenomeCatalog, String), EngineError> {
        let catalog_path = catalog_path
            .map(str::trim)
            .filter(|v| !v.is_empty())
            .unwrap_or(DEFAULT_GENOME_CATALOG_PATH)
            .to_string();
        let catalog = GenomeCatalog::from_json_file(&catalog_path).map_err(|e| EngineError {
            code: ErrorCode::InvalidInput,
            message: format!("Could not open genome catalog '{catalog_path}': {e}"),
        })?;
        Ok((catalog, catalog_path))
    }

    pub fn list_reference_genomes(catalog_path: Option<&str>) -> Result<Vec<String>, EngineError> {
        let (catalog, _) = Self::open_reference_genome_catalog(catalog_path)?;
        Ok(catalog.list_genomes())
    }

    pub fn describe_reference_genome_sources(
        catalog_path: Option<&str>,
        genome_id: &str,
        cache_dir: Option<&str>,
    ) -> Result<GenomeSourcePlan, EngineError> {
        let (catalog, catalog_path) = Self::open_reference_genome_catalog(catalog_path)?;
        catalog
            .source_plan(
                genome_id,
                cache_dir.map(str::trim).filter(|v| !v.is_empty()),
            )
            .map_err(|e| EngineError {
                code: ErrorCode::InvalidInput,
                message: format!(
                    "Could not resolve source plan for genome '{}' in catalog '{}': {}",
                    genome_id, catalog_path, e
                ),
            })
    }

    pub fn preview_reference_genome_ensembl_catalog_updates(
        catalog_path: Option<&str>,
    ) -> Result<EnsemblCatalogUpdatePreview, EngineError> {
        let (catalog, catalog_path) = Self::open_reference_genome_catalog(catalog_path)?;
        catalog
            .preview_ensembl_catalog_updates()
            .map_err(|e| EngineError {
                code: ErrorCode::Io,
                message: format!(
                    "Could not preview Ensembl catalog updates for '{}': {}",
                    catalog_path, e
                ),
            })
    }

    pub fn apply_reference_genome_ensembl_catalog_updates(
        catalog_path: Option<&str>,
        output_catalog_path: Option<&str>,
    ) -> Result<EnsemblCatalogUpdateReport, EngineError> {
        let (catalog, catalog_path) = Self::open_reference_genome_catalog(catalog_path)?;
        catalog
            .apply_ensembl_catalog_updates(output_catalog_path)
            .map_err(|e| EngineError {
                code: ErrorCode::Io,
                message: format!(
                    "Could not apply Ensembl catalog updates for '{}': {}",
                    catalog_path, e
                ),
            })
    }

    pub fn remove_reference_genome_catalog_entry(
        catalog_path: Option<&str>,
        genome_id: &str,
        output_catalog_path: Option<&str>,
    ) -> Result<GenomeCatalogEntryRemovalReport, EngineError> {
        let (catalog, catalog_path) = Self::open_reference_genome_catalog(catalog_path)?;
        catalog
            .remove_catalog_entry(genome_id, output_catalog_path)
            .map_err(|e| EngineError {
                code: ErrorCode::Io,
                message: format!(
                    "Could not remove genome catalog entry '{}' from '{}': {}",
                    genome_id, catalog_path, e
                ),
            })
    }

    pub fn remove_prepared_reference_genome(
        catalog_path: Option<&str>,
        genome_id: &str,
        cache_dir: Option<&str>,
    ) -> Result<PreparedGenomeRemovalReport, EngineError> {
        let (catalog, catalog_path) = Self::open_reference_genome_catalog(catalog_path)?;
        catalog
            .remove_prepared_genome_install(
                genome_id,
                cache_dir.map(str::trim).filter(|v| !v.is_empty()),
            )
            .map_err(|e| EngineError {
                code: ErrorCode::Io,
                message: format!(
                    "Could not remove prepared genome '{}' from '{}': {}",
                    genome_id, catalog_path, e
                ),
            })
    }

    pub fn inspect_prepared_cache_roots(
        cache_roots: &[String],
    ) -> Result<PreparedCacheInspectionReport, EngineError> {
        inspect_prepared_cache_roots(cache_roots).map_err(|e| EngineError {
            code: ErrorCode::Io,
            message: format!("Could not inspect prepared cache roots: {e}"),
        })
    }

    pub fn clear_prepared_cache_roots(
        request: &PreparedCacheCleanupRequest,
    ) -> Result<PreparedCacheCleanupReport, EngineError> {
        clear_prepared_cache_roots(request).map_err(|e| EngineError {
            code: ErrorCode::Io,
            message: format!("Could not clear prepared cache roots: {e}"),
        })
    }

    pub fn list_helper_genomes(catalog_path: Option<&str>) -> Result<Vec<String>, EngineError> {
        let chosen = catalog_path
            .map(str::trim)
            .filter(|v| !v.is_empty())
            .unwrap_or(DEFAULT_HELPER_GENOME_CATALOG_PATH);
        Self::list_reference_genomes(Some(chosen))
    }

    pub fn describe_helper_genome_sources(
        genome_id: &str,
        catalog_path: Option<&str>,
        cache_dir: Option<&str>,
    ) -> Result<GenomeSourcePlan, EngineError> {
        let chosen = catalog_path
            .map(str::trim)
            .filter(|v| !v.is_empty())
            .unwrap_or(DEFAULT_HELPER_GENOME_CATALOG_PATH);
        Self::describe_reference_genome_sources(Some(chosen), genome_id, cache_dir)
    }

    pub fn is_reference_genome_prepared(
        catalog_path: Option<&str>,
        genome_id: &str,
        cache_dir: Option<&str>,
    ) -> Result<bool, EngineError> {
        let (catalog, _) = Self::open_reference_genome_catalog(catalog_path)?;
        catalog
            .is_prepared(
                genome_id,
                cache_dir.map(str::trim).filter(|v| !v.is_empty()),
            )
            .map_err(|e| EngineError {
                code: ErrorCode::InvalidInput,
                message: format!(
                    "Could not check prepared status for genome '{}': {}",
                    genome_id, e
                ),
            })
    }

    pub fn is_helper_genome_prepared(
        genome_id: &str,
        catalog_path: Option<&str>,
        cache_dir: Option<&str>,
    ) -> Result<bool, EngineError> {
        let chosen = catalog_path
            .map(str::trim)
            .filter(|v| !v.is_empty())
            .unwrap_or(DEFAULT_HELPER_GENOME_CATALOG_PATH);
        Self::is_reference_genome_prepared(Some(chosen), genome_id, cache_dir)
    }

    pub fn list_reference_genome_genes(
        catalog_path: Option<&str>,
        genome_id: &str,
        cache_dir: Option<&str>,
    ) -> Result<Vec<GenomeGeneRecord>, EngineError> {
        let (catalog, catalog_path) = Self::open_reference_genome_catalog(catalog_path)?;
        catalog
            .list_gene_regions(
                genome_id,
                cache_dir.map(str::trim).filter(|v| !v.is_empty()),
            )
            .map_err(|e| EngineError {
                code: ErrorCode::InvalidInput,
                message: format!(
                    "Could not list genes for genome '{}' in catalog '{}': {}",
                    genome_id, catalog_path, e
                ),
            })
    }

    pub fn list_helper_genome_features(
        genome_id: &str,
        catalog_path: Option<&str>,
        cache_dir: Option<&str>,
    ) -> Result<Vec<GenomeGeneRecord>, EngineError> {
        let chosen = catalog_path
            .map(str::trim)
            .filter(|v| !v.is_empty())
            .unwrap_or(DEFAULT_HELPER_GENOME_CATALOG_PATH);
        Self::list_reference_genome_genes(Some(chosen), genome_id, cache_dir)
    }

    fn blast_reference_genome_raw(
        catalog_path: Option<&str>,
        genome_id: &str,
        query_sequence: &str,
        max_hits: usize,
        task: Option<&str>,
        cache_dir: Option<&str>,
        should_cancel: &mut dyn FnMut() -> bool,
    ) -> Result<GenomeBlastReport, EngineError> {
        let (catalog, catalog_path) = Self::open_reference_genome_catalog(catalog_path)?;
        catalog
            .blast_sequence_with_cache_and_cancel(
                genome_id,
                query_sequence,
                max_hits,
                task,
                cache_dir.map(str::trim).filter(|v| !v.is_empty()),
                should_cancel,
            )
            .map_err(|e| EngineError {
                code: ErrorCode::Io,
                message: format!(
                    "Could not run BLAST search against genome '{}' in catalog '{}': {}",
                    genome_id, catalog_path, e
                ),
            })
    }

    fn blast_helper_genome_raw(
        genome_id: &str,
        query_sequence: &str,
        max_hits: usize,
        task: Option<&str>,
        catalog_path: Option<&str>,
        cache_dir: Option<&str>,
        should_cancel: &mut dyn FnMut() -> bool,
    ) -> Result<GenomeBlastReport, EngineError> {
        let chosen = catalog_path
            .map(str::trim)
            .filter(|v| !v.is_empty())
            .unwrap_or(DEFAULT_HELPER_GENOME_CATALOG_PATH);
        Self::blast_reference_genome_raw(
            Some(chosen),
            genome_id,
            query_sequence,
            max_hits,
            task,
            cache_dir,
            should_cancel,
        )
    }

    fn normalize_blast_task_value(raw: &str) -> Result<String, EngineError> {
        let normalized = raw.trim().to_ascii_lowercase();
        match normalized.as_str() {
            "blastn-short" | "blastn" => Ok(normalized),
            _ => Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: format!(
                    "Unsupported BLAST task '{}'. Expected one of: blastn-short, blastn",
                    raw
                ),
            }),
        }
    }

    fn validate_blast_thresholds(thresholds: &BlastThresholdOptions) -> Result<(), EngineError> {
        if let Some(max_evalue) = thresholds.max_evalue {
            if !max_evalue.is_finite() || max_evalue < 0.0 {
                return Err(EngineError {
                    code: ErrorCode::InvalidInput,
                    message: "BLAST max_evalue must be a finite number >= 0.0".to_string(),
                });
            }
        }
        if let Some(identity) = thresholds.min_identity_percent {
            if !identity.is_finite() || !(0.0..=100.0).contains(&identity) {
                return Err(EngineError {
                    code: ErrorCode::InvalidInput,
                    message: "BLAST min_identity_percent must be within [0, 100]".to_string(),
                });
            }
        }
        if let Some(qcov) = thresholds.min_query_coverage_percent {
            if !qcov.is_finite() || !(0.0..=100.0).contains(&qcov) {
                return Err(EngineError {
                    code: ErrorCode::InvalidInput,
                    message: "BLAST min_query_coverage_percent must be within [0, 100]".to_string(),
                });
            }
        }
        if let Some(min_len) = thresholds.min_alignment_length_bp {
            if min_len == 0 {
                return Err(EngineError {
                    code: ErrorCode::InvalidInput,
                    message: "BLAST min_alignment_length_bp must be >= 1".to_string(),
                });
            }
        }
        if let Some(min_bitscore) = thresholds.min_bit_score {
            if !min_bitscore.is_finite() || min_bitscore < 0.0 {
                return Err(EngineError {
                    code: ErrorCode::InvalidInput,
                    message: "BLAST min_bit_score must be a finite number >= 0.0".to_string(),
                });
            }
        }
        Ok(())
    }

    fn parse_blast_run_options_from_value(
        value: &serde_json::Value,
        source: &str,
    ) -> Result<BlastRunOptions, EngineError> {
        if value.is_null() {
            return Ok(BlastRunOptions::default());
        }
        if !value.is_object() {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: format!("BLAST options in {source} must be a JSON object"),
            });
        }
        serde_json::from_value::<BlastRunOptions>(value.clone()).map_err(|e| EngineError {
            code: ErrorCode::InvalidInput,
            message: format!("Invalid BLAST options in {source}: {e}"),
        })
    }

    fn built_in_blast_run_options() -> BlastRunOptions {
        BlastRunOptions {
            task: Some("blastn-short".to_string()),
            max_hits: Some(25),
            thresholds: BlastThresholdOptions::default(),
        }
    }

    fn load_blast_options_from_file(
        path: &str,
        strict_missing: bool,
    ) -> Result<Option<BlastRunOptions>, EngineError> {
        let trimmed = path.trim();
        if trimmed.is_empty() {
            return Ok(None);
        }
        let p = Path::new(trimmed);
        if !p.exists() {
            if strict_missing {
                return Err(EngineError {
                    code: ErrorCode::InvalidInput,
                    message: format!("BLAST defaults file '{}' does not exist", trimmed),
                });
            }
            return Ok(None);
        }
        let text = std::fs::read_to_string(p).map_err(|e| EngineError {
            code: ErrorCode::Io,
            message: format!("Could not read BLAST defaults file '{}': {}", trimmed, e),
        })?;
        let value = serde_json::from_str::<serde_json::Value>(&text).map_err(|e| EngineError {
            code: ErrorCode::InvalidInput,
            message: format!("Could not parse BLAST defaults file '{}': {}", trimmed, e),
        })?;
        let parsed = Self::parse_blast_run_options_from_value(
            &value,
            &format!("defaults file '{}'", trimmed),
        )?;
        Ok(Some(parsed))
    }

    pub fn resolve_blast_options_with_layers(
        defaults_path: Option<&str>,
        project_override_json: Option<&serde_json::Value>,
        request_override_json: Option<&serde_json::Value>,
        legacy_task: Option<&str>,
        legacy_max_hits: Option<usize>,
    ) -> Result<BlastResolvedOptions, EngineError> {
        let mut merged = Self::built_in_blast_run_options();

        let explicit_defaults_path = defaults_path.map(str::trim).filter(|v| !v.is_empty());
        let defaults_path = explicit_defaults_path.unwrap_or(DEFAULT_BLAST_OPTIONS_PATH);
        if let Some(from_file) =
            Self::load_blast_options_from_file(defaults_path, explicit_defaults_path.is_some())?
        {
            merged.merge_from(&from_file);
        }

        if let Some(project_json) = project_override_json {
            let from_project =
                Self::parse_blast_run_options_from_value(project_json, "project override")?;
            merged.merge_from(&from_project);
        }

        if legacy_task.is_some() || legacy_max_hits.is_some() {
            let legacy = BlastRunOptions {
                task: legacy_task.map(|v| v.to_string()),
                max_hits: legacy_max_hits,
                thresholds: BlastThresholdOptions::default(),
            };
            merged.merge_from(&legacy);
        }

        if let Some(request_json) = request_override_json {
            let from_request =
                Self::parse_blast_run_options_from_value(request_json, "request override")?;
            merged.merge_from(&from_request);
        }

        let task = Self::normalize_blast_task_value(
            merged.task.as_deref().unwrap_or("blastn-short").trim(),
        )?;
        let max_hits = merged.max_hits.unwrap_or(25);
        if max_hits == 0 {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: "BLAST max_hits must be >= 1".to_string(),
            });
        }
        Self::validate_blast_thresholds(&merged.thresholds)?;
        Ok(BlastResolvedOptions {
            task,
            max_hits,
            thresholds: merged.thresholds,
        })
    }

    fn blast_report_passes_thresholds(
        hit: &crate::genomes::BlastHit,
        t: &BlastThresholdOptions,
    ) -> bool {
        if let Some(max_evalue) = t.max_evalue {
            if hit.evalue > max_evalue {
                return false;
            }
        }
        if let Some(min_identity_percent) = t.min_identity_percent {
            if hit.identity_percent < min_identity_percent {
                return false;
            }
        }
        if let Some(min_query_coverage_percent) = t.min_query_coverage_percent {
            match hit.query_coverage_percent {
                Some(v) if v >= min_query_coverage_percent => {}
                _ => return false,
            }
        }
        if let Some(min_alignment_length_bp) = t.min_alignment_length_bp {
            if hit.alignment_length < min_alignment_length_bp {
                return false;
            }
        }
        if let Some(min_bit_score) = t.min_bit_score {
            if hit.bit_score < min_bit_score {
                return false;
            }
        }
        true
    }

    fn apply_blast_thresholds_to_report(
        report: &mut GenomeBlastReport,
        thresholds: &BlastThresholdOptions,
    ) -> Result<(), EngineError> {
        let before = report.hits.len();
        report
            .hits
            .retain(|hit| Self::blast_report_passes_thresholds(hit, thresholds));
        let removed = before.saturating_sub(report.hits.len());
        report.hit_count = report.hits.len();
        if removed > 0 {
            report.warnings.push(format!(
                "BLAST thresholds removed {} hit(s) (remaining={})",
                removed, report.hit_count
            ));
        }
        if thresholds.unique_best_hit.unwrap_or(false) && report.hit_count != 1 {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: format!(
                    "BLAST unique_best_hit requires exactly one remaining hit, found {}",
                    report.hit_count
                ),
            });
        }
        Ok(())
    }

    pub fn resolve_blast_options_for_request(
        &self,
        request_override_json: Option<&serde_json::Value>,
        legacy_task: Option<&str>,
        legacy_max_hits: Option<usize>,
    ) -> Result<BlastResolvedOptions, EngineError> {
        let defaults_path = self
            .state
            .metadata
            .get(BLAST_OPTIONS_DEFAULTS_PATH_METADATA_KEY)
            .and_then(|v| v.as_str())
            .map(|v| v.to_string());
        let project_override = self.state.metadata.get(BLAST_OPTIONS_OVERRIDE_METADATA_KEY);
        Self::resolve_blast_options_with_layers(
            defaults_path.as_deref(),
            project_override,
            request_override_json,
            legacy_task,
            legacy_max_hits,
        )
    }

    pub fn blast_external_binary_preflight_report(&self) -> BlastExternalBinaryPreflightReport {
        blast_external_binary_preflight_report()
    }

    pub fn blast_reference_genome_with_project_and_request_options(
        &self,
        catalog_path: Option<&str>,
        genome_id: &str,
        query_sequence: &str,
        request_override_json: Option<&serde_json::Value>,
        legacy_task: Option<&str>,
        legacy_max_hits: Option<usize>,
        cache_dir: Option<&str>,
    ) -> Result<GenomeBlastReport, EngineError> {
        let mut never_cancel = || false;
        self.blast_reference_genome_with_project_and_request_options_and_cancel(
            catalog_path,
            genome_id,
            query_sequence,
            request_override_json,
            legacy_task,
            legacy_max_hits,
            cache_dir,
            &mut never_cancel,
        )
    }

    pub fn blast_reference_genome_with_project_and_request_options_and_cancel(
        &self,
        catalog_path: Option<&str>,
        genome_id: &str,
        query_sequence: &str,
        request_override_json: Option<&serde_json::Value>,
        legacy_task: Option<&str>,
        legacy_max_hits: Option<usize>,
        cache_dir: Option<&str>,
        should_cancel: &mut dyn FnMut() -> bool,
    ) -> Result<GenomeBlastReport, EngineError> {
        let resolved = self.resolve_blast_options_for_request(
            request_override_json,
            legacy_task,
            legacy_max_hits,
        )?;
        let mut report = Self::blast_reference_genome_raw(
            catalog_path,
            genome_id,
            query_sequence,
            resolved.max_hits,
            Some(&resolved.task),
            cache_dir,
            should_cancel,
        )?;
        Self::apply_blast_thresholds_to_report(&mut report, &resolved.thresholds)?;
        report.options_override_json = request_override_json.cloned();
        report.effective_options_json = serde_json::to_value(&resolved).ok();
        Ok(report)
    }

    pub fn blast_helper_genome_with_project_and_request_options(
        &self,
        genome_id: &str,
        query_sequence: &str,
        request_override_json: Option<&serde_json::Value>,
        legacy_task: Option<&str>,
        legacy_max_hits: Option<usize>,
        catalog_path: Option<&str>,
        cache_dir: Option<&str>,
    ) -> Result<GenomeBlastReport, EngineError> {
        let mut never_cancel = || false;
        self.blast_helper_genome_with_project_and_request_options_and_cancel(
            genome_id,
            query_sequence,
            request_override_json,
            legacy_task,
            legacy_max_hits,
            catalog_path,
            cache_dir,
            &mut never_cancel,
        )
    }

    pub fn blast_helper_genome_with_project_and_request_options_and_cancel(
        &self,
        genome_id: &str,
        query_sequence: &str,
        request_override_json: Option<&serde_json::Value>,
        legacy_task: Option<&str>,
        legacy_max_hits: Option<usize>,
        catalog_path: Option<&str>,
        cache_dir: Option<&str>,
        should_cancel: &mut dyn FnMut() -> bool,
    ) -> Result<GenomeBlastReport, EngineError> {
        let resolved = self.resolve_blast_options_for_request(
            request_override_json,
            legacy_task,
            legacy_max_hits,
        )?;
        let mut report = Self::blast_helper_genome_raw(
            genome_id,
            query_sequence,
            resolved.max_hits,
            Some(&resolved.task),
            catalog_path,
            cache_dir,
            should_cancel,
        )?;
        Self::apply_blast_thresholds_to_report(&mut report, &resolved.thresholds)?;
        report.options_override_json = request_override_json.cloned();
        report.effective_options_json = serde_json::to_value(&resolved).ok();
        Ok(report)
    }

    pub fn blast_reference_genome_with_request_options(
        catalog_path: Option<&str>,
        genome_id: &str,
        query_sequence: &str,
        request_override_json: Option<&serde_json::Value>,
        cache_dir: Option<&str>,
    ) -> Result<GenomeBlastReport, EngineError> {
        let engine = Self::default();
        engine.blast_reference_genome_with_project_and_request_options(
            catalog_path,
            genome_id,
            query_sequence,
            request_override_json,
            None,
            None,
            cache_dir,
        )
    }

    pub fn blast_helper_genome_with_request_options(
        genome_id: &str,
        query_sequence: &str,
        request_override_json: Option<&serde_json::Value>,
        catalog_path: Option<&str>,
        cache_dir: Option<&str>,
    ) -> Result<GenomeBlastReport, EngineError> {
        let engine = Self::default();
        engine.blast_helper_genome_with_project_and_request_options(
            genome_id,
            query_sequence,
            request_override_json,
            None,
            None,
            catalog_path,
            cache_dir,
        )
    }

    pub fn blast_reference_genome(
        catalog_path: Option<&str>,
        genome_id: &str,
        query_sequence: &str,
        max_hits: usize,
        task: Option<&str>,
        cache_dir: Option<&str>,
    ) -> Result<GenomeBlastReport, EngineError> {
        let mut request = serde_json::Map::new();
        request.insert("max_hits".to_string(), serde_json::Value::from(max_hits));
        if let Some(task) = task.map(str::trim).filter(|v| !v.is_empty()) {
            request.insert("task".to_string(), serde_json::Value::from(task));
        }
        let request_json = serde_json::Value::Object(request);
        Self::blast_reference_genome_with_request_options(
            catalog_path,
            genome_id,
            query_sequence,
            Some(&request_json),
            cache_dir,
        )
    }

    pub fn blast_helper_genome(
        genome_id: &str,
        query_sequence: &str,
        max_hits: usize,
        task: Option<&str>,
        catalog_path: Option<&str>,
        cache_dir: Option<&str>,
    ) -> Result<GenomeBlastReport, EngineError> {
        let mut request = serde_json::Map::new();
        request.insert("max_hits".to_string(), serde_json::Value::from(max_hits));
        if let Some(task) = task.map(str::trim).filter(|v| !v.is_empty()) {
            request.insert("task".to_string(), serde_json::Value::from(task));
        }
        let request_json = serde_json::Value::Object(request);
        Self::blast_helper_genome_with_request_options(
            genome_id,
            query_sequence,
            Some(&request_json),
            catalog_path,
            cache_dir,
        )
    }

    fn prepare_reference_genome_plan_with_options(
        genome_id: &str,
        catalog_path: Option<&str>,
        cache_dir: Option<&str>,
        mode: PrepareReferenceGenomeMode,
    ) -> Result<PrepareGenomePlan, EngineError> {
        let (catalog, _) = Self::open_reference_genome_catalog(catalog_path)?;
        let result = match mode {
            PrepareReferenceGenomeMode::PrepareOrReuse => catalog.prepare_genome_plan(
                genome_id,
                cache_dir.map(str::trim).filter(|v| !v.is_empty()),
            ),
            PrepareReferenceGenomeMode::ReindexCachedFiles => catalog.reindex_genome_plan(
                genome_id,
                cache_dir.map(str::trim).filter(|v| !v.is_empty()),
            ),
            PrepareReferenceGenomeMode::RefreshFromSources => catalog.reinstall_genome_plan(
                genome_id,
                cache_dir.map(str::trim).filter(|v| !v.is_empty()),
            ),
        };
        result.map_err(|e| EngineError {
            code: ErrorCode::Io,
            message: match mode {
                PrepareReferenceGenomeMode::PrepareOrReuse => {
                    format!("Could not plan genome preparation for '{genome_id}': {e}")
                }
                PrepareReferenceGenomeMode::ReindexCachedFiles => {
                    format!("Could not plan genome reindex for '{genome_id}': {e}")
                }
                PrepareReferenceGenomeMode::RefreshFromSources => {
                    format!("Could not plan genome refresh for '{genome_id}': {e}")
                }
            },
        })
    }

    fn prepare_reference_genome_once_with_options(
        genome_id: &str,
        catalog_path: Option<&str>,
        cache_dir: Option<&str>,
        timeout_seconds: Option<u64>,
        mode: PrepareReferenceGenomeMode,
        on_progress: &mut dyn FnMut(PrepareGenomeProgress) -> bool,
    ) -> Result<PrepareGenomeReport, EngineError> {
        let (catalog, _) = Self::open_reference_genome_catalog(catalog_path)?;
        let timeout = timeout_seconds.map(Duration::from_secs);
        let started = Instant::now();
        let mut timed_out = false;
        let mut guarded_progress = |progress: PrepareGenomeProgress| -> bool {
            if let Some(limit) = timeout {
                if started.elapsed() >= limit {
                    timed_out = true;
                    return false;
                }
            }
            on_progress(progress)
        };
        let action_label = match mode {
            PrepareReferenceGenomeMode::PrepareOrReuse => "preparation",
            PrepareReferenceGenomeMode::ReindexCachedFiles => "reindex",
            PrepareReferenceGenomeMode::RefreshFromSources => "refresh",
        };
        let result = match mode {
            PrepareReferenceGenomeMode::PrepareOrReuse => catalog
                .prepare_genome_once_with_progress(
                    genome_id,
                    cache_dir.map(str::trim).filter(|v| !v.is_empty()),
                    &mut guarded_progress,
                ),
            PrepareReferenceGenomeMode::ReindexCachedFiles => catalog
                .reindex_genome_once_with_progress(
                    genome_id,
                    cache_dir.map(str::trim).filter(|v| !v.is_empty()),
                    &mut guarded_progress,
                ),
            PrepareReferenceGenomeMode::RefreshFromSources => catalog
                .reinstall_genome_once_with_progress(
                    genome_id,
                    cache_dir.map(str::trim).filter(|v| !v.is_empty()),
                    &mut guarded_progress,
                ),
        };
        result.map_err(|e| EngineError {
            code: ErrorCode::Io,
            message: if is_prepare_cancelled_error(&e) {
                if timed_out {
                    format!(
                        "Genome {action_label} timed out for '{genome_id}' after {} second(s)",
                        timeout_seconds.unwrap_or(0)
                    )
                } else {
                    format!("Genome {action_label} cancelled for '{genome_id}'")
                }
            } else {
                format!("Could not {action_label} genome '{genome_id}': {e}")
            },
        })
    }

    pub fn prepare_reference_genome_plan(
        genome_id: &str,
        catalog_path: Option<&str>,
        cache_dir: Option<&str>,
    ) -> Result<PrepareGenomePlan, EngineError> {
        Self::prepare_reference_genome_plan_with_options(
            genome_id,
            catalog_path,
            cache_dir,
            PrepareReferenceGenomeMode::PrepareOrReuse,
        )
    }

    pub fn reindex_reference_genome_plan(
        genome_id: &str,
        catalog_path: Option<&str>,
        cache_dir: Option<&str>,
    ) -> Result<PrepareGenomePlan, EngineError> {
        Self::prepare_reference_genome_plan_with_options(
            genome_id,
            catalog_path,
            cache_dir,
            PrepareReferenceGenomeMode::ReindexCachedFiles,
        )
    }

    pub fn reinstall_reference_genome_plan(
        genome_id: &str,
        catalog_path: Option<&str>,
        cache_dir: Option<&str>,
    ) -> Result<PrepareGenomePlan, EngineError> {
        Self::prepare_reference_genome_plan_with_options(
            genome_id,
            catalog_path,
            cache_dir,
            PrepareReferenceGenomeMode::RefreshFromSources,
        )
    }

    pub fn prepare_reference_genome_once(
        genome_id: &str,
        catalog_path: Option<&str>,
        cache_dir: Option<&str>,
        timeout_seconds: Option<u64>,
        on_progress: &mut dyn FnMut(PrepareGenomeProgress) -> bool,
    ) -> Result<PrepareGenomeReport, EngineError> {
        Self::prepare_reference_genome_once_with_options(
            genome_id,
            catalog_path,
            cache_dir,
            timeout_seconds,
            PrepareReferenceGenomeMode::PrepareOrReuse,
            on_progress,
        )
    }

    pub fn reindex_reference_genome_once(
        genome_id: &str,
        catalog_path: Option<&str>,
        cache_dir: Option<&str>,
        timeout_seconds: Option<u64>,
        on_progress: &mut dyn FnMut(PrepareGenomeProgress) -> bool,
    ) -> Result<PrepareGenomeReport, EngineError> {
        Self::prepare_reference_genome_once_with_options(
            genome_id,
            catalog_path,
            cache_dir,
            timeout_seconds,
            PrepareReferenceGenomeMode::ReindexCachedFiles,
            on_progress,
        )
    }

    pub fn reinstall_reference_genome_once(
        genome_id: &str,
        catalog_path: Option<&str>,
        cache_dir: Option<&str>,
        timeout_seconds: Option<u64>,
        on_progress: &mut dyn FnMut(PrepareGenomeProgress) -> bool,
    ) -> Result<PrepareGenomeReport, EngineError> {
        Self::prepare_reference_genome_once_with_options(
            genome_id,
            catalog_path,
            cache_dir,
            timeout_seconds,
            PrepareReferenceGenomeMode::RefreshFromSources,
            on_progress,
        )
    }

    pub fn prepare_helper_genome_once(
        genome_id: &str,
        catalog_path: Option<&str>,
        cache_dir: Option<&str>,
        timeout_seconds: Option<u64>,
        on_progress: &mut dyn FnMut(PrepareGenomeProgress) -> bool,
    ) -> Result<PrepareGenomeReport, EngineError> {
        let chosen = catalog_path
            .map(str::trim)
            .filter(|v| !v.is_empty())
            .unwrap_or(DEFAULT_HELPER_GENOME_CATALOG_PATH);
        Self::prepare_reference_genome_once(
            genome_id,
            Some(chosen),
            cache_dir,
            timeout_seconds,
            on_progress,
        )
    }

    pub fn format_prepare_genome_message(
        genome_id: &str,
        cache_dir: Option<&str>,
        report: &PrepareGenomeReport,
    ) -> String {
        let status = if report.reused_existing {
            "reused existing local cache"
        } else {
            "downloaded and installed"
        };
        let sequence_type = report.sequence_source_type.as_deref().unwrap_or("unknown");
        let annotation_type = report
            .annotation_source_type
            .as_deref()
            .unwrap_or("unknown");
        let blast_status = if report.blast_index_ready {
            format!(
                "ready ({})",
                report
                    .blast_db_prefix
                    .as_deref()
                    .unwrap_or("unknown BLAST DB prefix")
            )
        } else {
            "not available".to_string()
        };
        let cached_sequence_status = if report.cached_contig_count == 0 {
            "unavailable".to_string()
        } else {
            let longest = report.cached_longest_contig.as_deref().unwrap_or("unknown");
            let longest_bp = report.cached_longest_contig_bp.unwrap_or(0);
            let preview = if report.cached_contig_preview.is_empty() {
                String::new()
            } else {
                format!(" [{}]", report.cached_contig_preview.join(", "))
            };
            format!(
                "{} contigs, total_span={} bp, longest={} ({} bp){}",
                report.cached_contig_count,
                report.cached_total_span_bp,
                longest,
                longest_bp,
                preview
            )
        };
        format!(
            "Prepared genome '{}' ({status}). cache='{}' sequence='{}' [{}], annotation='{}' [{}], blast_index={}, cached_sequence={}",
            genome_id,
            cache_dir
                .map(str::trim)
                .filter(|v| !v.is_empty())
                .unwrap_or("catalog/default"),
            report.sequence_path,
            sequence_type,
            report.annotation_path,
            annotation_type,
            blast_status,
            cached_sequence_status
        )
    }

    /// Return the immutable deterministic operation journal for this engine.
    pub fn operation_log(&self) -> &[OperationRecord] {
        &self.journal
    }

    fn history_limit_or_default(&self) -> usize {
        if self.history_limit == 0 {
            Self::default_history_limit()
        } else {
            self.history_limit
        }
    }

    fn capture_history_checkpoint(&self) -> EngineHistoryCheckpoint {
        EngineHistoryCheckpoint {
            state: self.state.clone(),
            journal: self.journal.clone(),
            op_counter: self.op_counter,
        }
    }

    fn restore_history_checkpoint(&mut self, checkpoint: EngineHistoryCheckpoint) {
        self.state = checkpoint.state;
        self.journal = checkpoint.journal;
        self.op_counter = checkpoint.op_counter;
        self.reconcile_lineage_nodes();
        self.reconcile_containers();
    }

    fn op_records_history_checkpoint(op: &Operation) -> bool {
        !matches!(
            op,
            Operation::SaveFile { .. }
                | Operation::RenderSequenceSvg { .. }
                | Operation::RenderDotplotSvg { .. }
                | Operation::RenderFeatureExpertSvg { .. }
                | Operation::RenderIsoformArchitectureSvg { .. }
                | Operation::RenderRnaStructureSvg { .. }
                | Operation::RenderLineageSvg { .. }
                | Operation::RenderPoolGelSvg { .. }
                | Operation::RenderProtocolCartoonSvg { .. }
                | Operation::RenderProtocolCartoonTemplateSvg { .. }
                | Operation::ValidateProtocolCartoonTemplate { .. }
                | Operation::RenderProtocolCartoonTemplateWithBindingsSvg { .. }
                | Operation::ExportProtocolCartoonTemplateJson { .. }
                | Operation::ExportDnaLadders { .. }
                | Operation::ExportRnaLadders { .. }
                | Operation::ExportPool { .. }
                | Operation::ExportProcessRunBundle { .. }
                | Operation::ListRnaReadReports { .. }
                | Operation::ShowRnaReadReport { .. }
                | Operation::SummarizeRnaReadGeneSupport { .. }
                | Operation::ExportRnaReadReport { .. }
                | Operation::ExportRnaReadHitsFasta { .. }
                | Operation::ExportRnaReadSampleSheet { .. }
                | Operation::ExportRnaReadExonPathsTsv { .. }
                | Operation::ExportRnaReadExonAbundanceTsv { .. }
                | Operation::ExportRnaReadScoreDensitySvg { .. }
                | Operation::ExportRnaReadAlignmentsTsv { .. }
                | Operation::ExportRnaReadAlignmentDotplotSvg { .. }
                | Operation::ListSequencingConfirmationReports { .. }
                | Operation::ShowSequencingConfirmationReport { .. }
                | Operation::ExportSequencingConfirmationReport { .. }
                | Operation::ExportSequencingConfirmationSupportTsv { .. }
                | Operation::AlignSequences { .. }
        )
    }

    fn maybe_capture_checkpoint(&self, op: &Operation) -> Option<EngineHistoryCheckpoint> {
        Self::op_records_history_checkpoint(op).then(|| self.capture_history_checkpoint())
    }

    fn push_undo_checkpoint(&mut self, checkpoint: EngineHistoryCheckpoint) {
        self.undo_stack.push(checkpoint);
        let limit = self.history_limit_or_default();
        if self.undo_stack.len() > limit {
            let drain_len = self.undo_stack.len() - limit;
            self.undo_stack.drain(0..drain_len);
        }
        self.redo_stack.clear();
    }

    pub fn undo_available(&self) -> usize {
        self.undo_stack.len()
    }

    pub fn redo_available(&self) -> usize {
        self.redo_stack.len()
    }

    pub fn undo_last_operation(&mut self) -> Result<(), EngineError> {
        let Some(previous) = self.undo_stack.pop() else {
            return Err(EngineError {
                code: ErrorCode::NotFound,
                message: "No operation to undo".to_string(),
            });
        };
        let current = self.capture_history_checkpoint();
        self.redo_stack.push(current);
        let limit = self.history_limit_or_default();
        if self.redo_stack.len() > limit {
            let drain_len = self.redo_stack.len() - limit;
            self.redo_stack.drain(0..drain_len);
        }
        self.restore_history_checkpoint(previous);
        Ok(())
    }

    pub fn redo_last_operation(&mut self) -> Result<(), EngineError> {
        let Some(next) = self.redo_stack.pop() else {
            return Err(EngineError {
                code: ErrorCode::NotFound,
                message: "No operation to redo".to_string(),
            });
        };
        let current = self.capture_history_checkpoint();
        self.undo_stack.push(current);
        let limit = self.history_limit_or_default();
        if self.undo_stack.len() > limit {
            let drain_len = self.undo_stack.len() - limit;
            self.undo_stack.drain(0..drain_len);
        }
        self.restore_history_checkpoint(next);
        Ok(())
    }

    pub fn apply_with_progress<F>(
        &mut self,
        op: Operation,
        mut on_progress: F,
    ) -> Result<OpResult, EngineError>
    where
        F: FnMut(OperationProgress) -> bool,
    {
        let run_id = "interactive".to_string();
        let checkpoint = self.maybe_capture_checkpoint(&op);
        let result = self.apply_internal(op.clone(), &run_id, &mut on_progress)?;
        self.journal.push(OperationRecord {
            run_id,
            op,
            result: result.clone(),
        });
        if let Some(checkpoint) = checkpoint {
            self.push_undo_checkpoint(checkpoint);
        }
        Ok(result)
    }

    /// Apply a workflow while streaming long-running progress payloads.
    ///
    /// Progress callbacks are cooperative: returning `false` requests
    /// cancellation for operations that support it.
    pub fn apply_workflow_with_progress<F>(
        &mut self,
        wf: Workflow,
        mut on_progress: F,
    ) -> Result<Vec<OpResult>, EngineError>
    where
        F: FnMut(OperationProgress) -> bool,
    {
        let mut results = Vec::new();
        for op in &wf.ops {
            let checkpoint = self.maybe_capture_checkpoint(op);
            let result = self.apply_internal(op.clone(), &wf.run_id, &mut on_progress)?;
            self.journal.push(OperationRecord {
                run_id: wf.run_id.clone(),
                op: op.clone(),
                result: result.clone(),
            });
            if let Some(checkpoint) = checkpoint {
                self.push_undo_checkpoint(checkpoint);
            }
            results.push(result);
        }
        Ok(results)
    }

    fn next_op_id(&mut self) -> OpId {
        self.op_counter += 1;
        format!("op-{}", self.op_counter)
    }

    fn derive_seq_id(path: &str) -> SeqId {
        Path::new(path)
            .file_stem()
            .map(|s| s.to_string_lossy().to_string())
            .filter(|s| !s.is_empty())
            .unwrap_or_else(|| "sequence".to_string())
    }

    fn normalize_id_token(raw: &str) -> String {
        let mut out = String::new();
        for c in raw.chars() {
            if c.is_ascii_alphanumeric() {
                out.push(c.to_ascii_lowercase());
            } else if matches!(c, '_' | '-' | '.') {
                out.push('_');
            }
        }
        let trimmed = out.trim_matches('_');
        if trimmed.is_empty() {
            "region".to_string()
        } else {
            trimmed.to_string()
        }
    }

    fn genome_gene_matches_exact(record: &GenomeGeneRecord, query: &str) -> bool {
        record
            .gene_name
            .as_ref()
            .map(|v| v.eq_ignore_ascii_case(query))
            .unwrap_or(false)
            || record
                .gene_id
                .as_ref()
                .map(|v| v.eq_ignore_ascii_case(query))
                .unwrap_or(false)
    }

    fn genome_gene_matches_contains(record: &GenomeGeneRecord, query_lower: &str) -> bool {
        record
            .gene_name
            .as_ref()
            .map(|v| v.to_ascii_lowercase().contains(query_lower))
            .unwrap_or(false)
            || record
                .gene_id
                .as_ref()
                .map(|v| v.to_ascii_lowercase().contains(query_lower))
                .unwrap_or(false)
    }

    fn genome_gene_display_label(record: &GenomeGeneRecord) -> String {
        let label = record
            .gene_name
            .as_ref()
            .or(record.gene_id.as_ref())
            .cloned()
            .unwrap_or_else(|| "unnamed_gene".to_string());
        format!(
            "{}:{}-{} ({})",
            record.chromosome, record.start_1based, record.end_1based, label
        )
    }

    fn default_extract_genome_gene_output_id(
        genome_id: &str,
        gene: &GenomeGeneRecord,
        extract_mode: GenomeGeneExtractMode,
        promoter_upstream_bp: usize,
    ) -> String {
        let label = gene
            .gene_name
            .as_deref()
            .or(gene.gene_id.as_deref())
            .unwrap_or("gene");
        let genome_token = Self::normalize_id_token(genome_id);
        let label_token = Self::normalize_id_token(label);
        match extract_mode {
            GenomeGeneExtractMode::Gene => format!(
                "{}_{}_{}_{}",
                genome_token, label_token, gene.start_1based, gene.end_1based
            ),
            GenomeGeneExtractMode::CodingWithPromoter => {
                if promoter_upstream_bp == 0 {
                    format!("{genome_token}_{label_token}_coding")
                } else {
                    format!(
                        "{genome_token}_{label_token}_coding_promoter_{}bp",
                        promoter_upstream_bp
                    )
                }
            }
        }
    }

    fn resolve_extract_genome_gene_interval(
        selected_gene: &GenomeGeneRecord,
        transcript_records: &[GenomeTranscriptRecord],
        extract_mode: GenomeGeneExtractMode,
        promoter_upstream_bp: usize,
    ) -> Result<(usize, usize), String> {
        match extract_mode {
            GenomeGeneExtractMode::Gene => {
                Ok((selected_gene.start_1based, selected_gene.end_1based))
            }
            GenomeGeneExtractMode::CodingWithPromoter => {
                let strand = selected_gene
                    .strand
                    .or_else(|| transcript_records.iter().filter_map(|record| record.strand).next())
                    .ok_or_else(|| {
                        format!(
                            "Gene '{}' has no strand annotation; cannot resolve CDS/promoter interval",
                            Self::genome_gene_display_label(selected_gene)
                        )
                    })?;
                let mut cds_start_1based: Option<usize> = None;
                let mut cds_end_1based: Option<usize> = None;
                for record in transcript_records {
                    for (start_1based, end_1based) in &record.cds_1based {
                        cds_start_1based = Some(
                            cds_start_1based
                                .map(|current| current.min(*start_1based))
                                .unwrap_or(*start_1based),
                        );
                        cds_end_1based = Some(
                            cds_end_1based
                                .map(|current| current.max(*end_1based))
                                .unwrap_or(*end_1based),
                        );
                    }
                }
                let Some(cds_start_1based) = cds_start_1based else {
                    return Err(format!(
                        "Gene '{}' has no CDS annotation; extract_mode=coding_with_promoter requires CDS-bearing transcripts",
                        Self::genome_gene_display_label(selected_gene)
                    ));
                };
                let Some(cds_end_1based) = cds_end_1based else {
                    return Err(format!(
                        "Gene '{}' has no CDS annotation; extract_mode=coding_with_promoter requires CDS-bearing transcripts",
                        Self::genome_gene_display_label(selected_gene)
                    ));
                };
                if strand == '-' {
                    Ok((
                        cds_start_1based,
                        cds_end_1based.saturating_add(promoter_upstream_bp),
                    ))
                } else {
                    Ok((
                        cds_start_1based.saturating_sub(promoter_upstream_bp).max(1),
                        cds_end_1based,
                    ))
                }
            }
        }
    }

    fn normalize_genome_chromosome_token(raw: &str) -> String {
        let canonical_source =
            Self::genome_accession_chromosome_alias(raw).unwrap_or_else(|| raw.trim().to_string());
        let trimmed = canonical_source.trim();
        let without_chr = trimmed
            .strip_prefix("chr")
            .or_else(|| trimmed.strip_prefix("Chr"))
            .or_else(|| trimmed.strip_prefix("CHR"))
            .unwrap_or(trimmed);
        match without_chr.to_ascii_lowercase().as_str() {
            "m" | "mt" => "mt".to_string(),
            other => other.to_string(),
        }
    }

    fn genome_accession_chromosome_alias(raw: &str) -> Option<String> {
        let trimmed = raw.trim();
        if trimmed.is_empty() {
            return None;
        }
        let upper = trimmed.to_ascii_uppercase();
        let without_version = upper.split('.').next().unwrap_or(&upper);
        let suffix = without_version
            .strip_prefix("NC_")
            .or_else(|| without_version.strip_prefix("CM"))?;
        if suffix.is_empty() || !suffix.chars().all(|ch| ch.is_ascii_digit()) {
            return None;
        }
        let normalized = suffix.trim_start_matches('0');
        Some(match normalized {
            "" => "0".to_string(),
            "23" => "x".to_string(),
            "24" => "y".to_string(),
            "12920" => "mt".to_string(),
            other => other.to_string(),
        })
    }

    fn genome_chromosome_match_tokens(raw: &str) -> std::collections::BTreeSet<String> {
        let mut tokens = std::collections::BTreeSet::new();
        let normalized = Self::normalize_genome_chromosome_token(raw);
        if !normalized.is_empty() {
            tokens.insert(normalized);
        }
        if let Some(alias) = Self::genome_accession_chromosome_alias(raw) {
            let alias_normalized = Self::normalize_genome_chromosome_token(&alias);
            if !alias_normalized.is_empty() {
                tokens.insert(alias_normalized);
            }
        }
        tokens
    }

    fn genome_chromosome_matches(left: &str, right: &str) -> bool {
        let left_trimmed = left.trim();
        let right_trimmed = right.trim();
        if left_trimmed.eq_ignore_ascii_case(right_trimmed) {
            return true;
        }
        let left_tokens = Self::genome_chromosome_match_tokens(left_trimmed);
        let right_tokens = Self::genome_chromosome_match_tokens(right_trimmed);
        left_tokens.iter().any(|token| right_tokens.contains(token))
    }

    fn resolve_extract_region_annotation_scope(
        annotation_scope: Option<GenomeAnnotationScope>,
        include_genomic_annotation: Option<bool>,
    ) -> GenomeAnnotationScope {
        if let Some(scope) = annotation_scope {
            return scope;
        }
        match include_genomic_annotation {
            Some(false) => GenomeAnnotationScope::None,
            Some(true) => GenomeAnnotationScope::Core,
            None => GenomeAnnotationScope::Core,
        }
    }

    fn helper_mcs_preset_hint(genome_id: &str) -> Option<&'static str> {
        let lowered = genome_id.to_ascii_lowercase();
        if lowered.contains("puc19") {
            Some("pUC19")
        } else if lowered.contains("puc18") {
            Some("pUC18")
        } else {
            None
        }
    }

    fn helper_mcs_sequence_by_preset(preset: &str) -> Option<&'static [u8]> {
        if preset.eq_ignore_ascii_case("pUC19") {
            Some(PUC19_MCS_SEQUENCE.as_bytes())
        } else if preset.eq_ignore_ascii_case("pUC18") {
            Some(PUC18_MCS_SEQUENCE.as_bytes())
        } else {
            None
        }
    }

    fn detect_helper_mcs(
        sequence: &[u8],
        preferred_preset: &str,
    ) -> Option<(&'static str, usize, usize)> {
        let mut candidates: Vec<(&'static str, &'static [u8])> = vec![];
        if preferred_preset.eq_ignore_ascii_case("pUC19") {
            candidates.push(("pUC19", PUC19_MCS_SEQUENCE.as_bytes()));
            candidates.push(("pUC18", PUC18_MCS_SEQUENCE.as_bytes()));
        } else if preferred_preset.eq_ignore_ascii_case("pUC18") {
            candidates.push(("pUC18", PUC18_MCS_SEQUENCE.as_bytes()));
            candidates.push(("pUC19", PUC19_MCS_SEQUENCE.as_bytes()));
        } else {
            return None;
        }
        for (preset, motif) in candidates {
            let hits = Self::find_all_subsequences(sequence, motif);
            if let Some(start_0based) = hits.first().copied() {
                return Some((preset, start_0based, hits.len()));
            }
        }
        None
    }

    fn is_generated_helper_mcs_feature(feature: &gb_io::seq::Feature) -> bool {
        feature
            .qualifier_values("gentle_generated")
            .any(|v| v.eq_ignore_ascii_case(HELPER_MCS_GENERATED_TAG))
    }

    fn feature_has_qualifier_value(feature: &gb_io::seq::Feature, key: &str, value: &str) -> bool {
        feature
            .qualifier_values(key)
            .any(|entry| entry.trim().eq_ignore_ascii_case(value))
    }

    fn mark_feature_as_genome_annotation_context(
        feature: &mut gb_io::seq::Feature,
        context_layer: &'static str,
    ) {
        if !Self::feature_has_qualifier_value(
            feature,
            "gentle_generated",
            GENOME_ANNOTATION_PROJECTION_GENERATED_TAG,
        ) {
            feature.qualifiers.push((
                "gentle_generated".into(),
                Some(GENOME_ANNOTATION_PROJECTION_GENERATED_TAG.to_string()),
            ));
        }
        if !Self::feature_has_qualifier_value(feature, CONTEXT_LAYER_QUALIFIER_KEY, context_layer) {
            feature.qualifiers.push((
                CONTEXT_LAYER_QUALIFIER_KEY.into(),
                Some(context_layer.to_string()),
            ));
        }
    }

    fn text_mentions_mcs(text: &str) -> bool {
        let lower = text.to_ascii_lowercase();
        lower.contains("multiple cloning site")
            || lower == "mcs"
            || lower.contains(" mcs ")
            || lower.starts_with("mcs ")
            || lower.ends_with(" mcs")
            || lower.contains("(mcs)")
    }

    fn normalize_enzyme_match_token(text: &str) -> String {
        text.chars()
            .filter(|c| c.is_ascii_alphanumeric())
            .map(|c| c.to_ascii_uppercase())
            .collect()
    }

    fn extract_rebase_enzyme_names_from_text(text: &str) -> Vec<String> {
        let tokens = text
            .split(|c: char| !c.is_ascii_alphanumeric())
            .map(str::trim)
            .filter(|value| !value.is_empty())
            .map(str::to_string)
            .collect::<Vec<_>>();
        if tokens.is_empty() {
            return vec![];
        }
        let mut lookup: HashMap<String, String> = HashMap::new();
        for enzyme in active_restriction_enzymes() {
            let normalized_name = Self::normalize_enzyme_match_token(&enzyme.name);
            if !normalized_name.is_empty() {
                lookup.entry(normalized_name).or_insert(enzyme.name);
            }
        }
        let mut out: Vec<String> = vec![];
        let mut seen: HashSet<String> = HashSet::new();
        let mut push_candidate = |candidate: String| {
            let normalized = Self::normalize_enzyme_match_token(&candidate);
            if let Some(name) = lookup.get(&normalized) {
                if seen.insert(name.clone()) {
                    out.push(name.clone());
                }
            }
        };
        for idx in 0..tokens.len() {
            push_candidate(tokens[idx].clone());
            if idx + 1 < tokens.len() {
                push_candidate(format!("{}{}", tokens[idx], tokens[idx + 1]));
            }
            if idx + 2 < tokens.len() {
                push_candidate(format!(
                    "{}{}{}",
                    tokens[idx],
                    tokens[idx + 1],
                    tokens[idx + 2]
                ));
            }
        }
        out
    }

    fn feature_looks_like_mcs(feature: &gb_io::seq::Feature) -> bool {
        if Self::is_generated_helper_mcs_feature(feature) {
            return true;
        }
        let text = Self::first_nonempty_feature_qualifier(
            feature,
            &["label", "note", "gene", "name", "standard_name"],
        )
        .unwrap_or_default();
        Self::text_mentions_mcs(&text)
    }

    fn feature_overlaps_span(
        feature: &gb_io::seq::Feature,
        start_0based: usize,
        end_0based_exclusive: usize,
    ) -> bool {
        if end_0based_exclusive <= start_0based {
            return false;
        }
        let mut ranges = vec![];
        collect_location_ranges_usize(&feature.location, &mut ranges);
        ranges.into_iter().any(|(part_start, part_end_exclusive)| {
            part_end_exclusive > start_0based && part_start < end_0based_exclusive
        })
    }

    fn build_helper_mcs_feature(
        preset: &str,
        start_0based: usize,
        end_0based_exclusive: usize,
        overlap_labels: &[String],
    ) -> gb_io::seq::Feature {
        let location =
            gb_io::seq::Location::simple_range(start_0based as i64, end_0based_exclusive as i64);
        let mut note = format!(
            "Multiple cloning site (MCS), {preset} orientation. Expected unique restriction sites: {}.",
            PUC_MCS_EXPECTED_SITES
        );
        if !overlap_labels.is_empty() {
            note.push_str(&format!(
                " Overlaps existing feature(s): {}. Insertion at this MCS may disrupt these features.",
                overlap_labels.join(", ")
            ));
        }
        gb_io::seq::Feature {
            kind: "misc_feature".into(),
            location,
            qualifiers: vec![
                ("label".into(), Some(format!("MCS ({preset})"))),
                ("note".into(), Some(note)),
                ("mcs_preset".into(), Some(preset.to_string())),
                (
                    "mcs_expected_sites".into(),
                    Some(PUC_MCS_EXPECTED_SITES.to_string()),
                ),
                (
                    "gentle_generated".into(),
                    Some(HELPER_MCS_GENERATED_TAG.to_string()),
                ),
            ],
        }
    }

    fn maybe_attach_known_helper_mcs_annotation(
        &mut self,
        seq_id: &str,
        genome_id: &str,
        result: &mut OpResult,
    ) {
        let Some(preferred_preset) = Self::helper_mcs_preset_hint(genome_id) else {
            return;
        };
        let Some(dna) = self.state.sequences.get_mut(seq_id) else {
            return;
        };
        if dna.features().iter().any(Self::feature_looks_like_mcs) {
            result.messages.push(format!(
                "Detected existing MCS annotation on '{}'; skipped fallback helper MCS annotation.",
                seq_id
            ));
            return;
        }
        let sequence = dna.get_forward_string().into_bytes();
        let Some((detected_preset, start_0based, match_count)) =
            Self::detect_helper_mcs(&sequence, preferred_preset)
        else {
            result.warnings.push(format!(
                "Could not auto-detect canonical MCS motif for helper genome '{}' on extracted sequence '{}'",
                genome_id, seq_id
            ));
            return;
        };
        let motif_len = Self::helper_mcs_sequence_by_preset(detected_preset)
            .map(|motif| motif.len())
            .unwrap_or(0usize);
        if motif_len == 0 || start_0based.saturating_add(motif_len) > sequence.len() {
            return;
        }
        if match_count != 1 {
            result.warnings.push(format!(
                "Found {match_count} canonical MCS motif matches on '{}'; expected exactly one, so no MCS fallback annotation was applied.",
                seq_id
            ));
            return;
        }
        let end_0based_exclusive = start_0based.saturating_add(motif_len);
        let overlap_labels: Vec<String> = dna
            .features()
            .iter()
            .enumerate()
            .filter_map(|(feature_id, feature)| {
                let kind = feature.kind.to_string().to_ascii_uppercase();
                if !matches!(kind.as_str(), "GENE" | "CDS" | "MRNA") {
                    return None;
                }
                if !Self::feature_overlaps_span(feature, start_0based, end_0based_exclusive) {
                    return None;
                }
                Some(Self::feature_display_label(feature, feature_id))
            })
            .collect();
        dna.features_mut().push(Self::build_helper_mcs_feature(
            detected_preset,
            start_0based,
            end_0based_exclusive,
            &overlap_labels,
        ));
        Self::prepare_sequence(dna);
        result.messages.push(format!(
            "Annotated helper MCS on '{}' using {} preset ({}..{}).",
            seq_id,
            detected_preset,
            start_0based.saturating_add(1),
            end_0based_exclusive
        ));
        if !detected_preset.eq_ignore_ascii_case(preferred_preset) {
            result.warnings.push(format!(
                "Helper genome '{}' suggests {} but matched {} MCS orientation; using detected orientation.",
                genome_id, preferred_preset, detected_preset
            ));
        }
    }

    fn genomic_interval_to_local_location(
        extracted_start_1based: usize,
        clipped_start_1based: usize,
        clipped_end_1based: usize,
        strand: Option<char>,
    ) -> Option<gb_io::seq::Location> {
        if clipped_end_1based < clipped_start_1based
            || clipped_start_1based < extracted_start_1based
        {
            return None;
        }
        let local_start_0based = clipped_start_1based.saturating_sub(extracted_start_1based);
        let local_end_exclusive = clipped_end_1based
            .saturating_sub(extracted_start_1based)
            .saturating_add(1);
        if local_end_exclusive <= local_start_0based {
            return None;
        }
        let base = gb_io::seq::Location::simple_range(
            local_start_0based as i64,
            local_end_exclusive as i64,
        );
        Some(if strand == Some('-') {
            gb_io::seq::Location::Complement(Box::new(base))
        } else {
            base
        })
    }

    fn gene_feature_from_genome_record(
        record: &GenomeGeneRecord,
        extracted_start_1based: usize,
        extracted_end_1based: usize,
    ) -> Option<gb_io::seq::Feature> {
        let clipped_start = record.start_1based.max(extracted_start_1based);
        let clipped_end = record.end_1based.min(extracted_end_1based);
        if clipped_end < clipped_start {
            return None;
        }
        let location = Self::genomic_interval_to_local_location(
            extracted_start_1based,
            clipped_start,
            clipped_end,
            record.strand,
        )?;
        let mut qualifiers = vec![
            ("chromosome".into(), Some(record.chromosome.clone())),
            (
                "genomic_start_1based".into(),
                Some(record.start_1based.to_string()),
            ),
            (
                "genomic_end_1based".into(),
                Some(record.end_1based.to_string()),
            ),
        ];
        if let Some(gene_name) = &record.gene_name {
            qualifiers.push(("gene".into(), Some(gene_name.clone())));
            qualifiers.push(("label".into(), Some(gene_name.clone())));
        }
        if let Some(gene_id) = &record.gene_id {
            qualifiers.push(("gene_id".into(), Some(gene_id.clone())));
            if record.gene_name.is_none() {
                qualifiers.push(("label".into(), Some(gene_id.clone())));
            }
        }
        if let Some(strand) = record.strand {
            qualifiers.push(("strand".into(), Some(strand.to_string())));
        }
        if let Some(biotype) = record
            .biotype
            .as_deref()
            .map(str::trim)
            .filter(|value| !value.is_empty())
        {
            qualifiers.push(("biotype".into(), Some(biotype.to_string())));
        }
        if let Some(mcs_text) = record
            .gene_name
            .as_ref()
            .or(record.gene_id.as_ref())
            .filter(|value| Self::text_mentions_mcs(value))
        {
            qualifiers.push(("note".into(), Some(mcs_text.clone())));
            let linked_enzymes = Self::extract_rebase_enzyme_names_from_text(mcs_text);
            if !linked_enzymes.is_empty() {
                qualifiers.push(("mcs_expected_sites".into(), Some(linked_enzymes.join(","))));
            }
        }
        let mut feature = gb_io::seq::Feature {
            kind: "gene".into(),
            location,
            qualifiers,
        };
        Self::mark_feature_as_genome_annotation_context(&mut feature, CONTEXT_LAYER_GENE);
        Some(feature)
    }

    fn transcript_subfeatures_from_genome_record(
        record: &GenomeTranscriptRecord,
        extracted_start_1based: usize,
        extracted_end_1based: usize,
    ) -> (Vec<gb_io::seq::Feature>, usize, usize) {
        let mut features: Vec<gb_io::seq::Feature> = vec![];
        let mut exons_attached = 0usize;
        let mut cds_attached = 0usize;

        let mut exon_ranges = record
            .exons_1based
            .iter()
            .filter_map(|(start, end)| {
                let clipped_start = (*start).max(extracted_start_1based);
                let clipped_end = (*end).min(extracted_end_1based);
                (clipped_end >= clipped_start).then_some((clipped_start, clipped_end))
            })
            .collect::<Vec<_>>();
        exon_ranges.sort_unstable_by(|a, b| a.0.cmp(&b.0).then(a.1.cmp(&b.1)));
        exon_ranges.dedup();
        for (idx, (start, end)) in exon_ranges.iter().enumerate() {
            let Some(location) = Self::genomic_interval_to_local_location(
                extracted_start_1based,
                *start,
                *end,
                record.strand,
            ) else {
                continue;
            };
            let mut qualifiers = vec![
                ("transcript_id".into(), Some(record.transcript_id.clone())),
                ("chromosome".into(), Some(record.chromosome.clone())),
                ("genomic_start_1based".into(), Some(start.to_string())),
                ("genomic_end_1based".into(), Some(end.to_string())),
                ("exon_number".into(), Some((idx + 1).to_string())),
            ];
            if let Some(gene_name) = &record.gene_name {
                qualifiers.push(("gene".into(), Some(gene_name.clone())));
            }
            if let Some(gene_id) = &record.gene_id {
                qualifiers.push(("gene_id".into(), Some(gene_id.clone())));
            }
            if let Some(strand) = record.strand {
                qualifiers.push(("strand".into(), Some(strand.to_string())));
            }
            let mut feature = gb_io::seq::Feature {
                kind: "exon".into(),
                location,
                qualifiers,
            };
            Self::mark_feature_as_genome_annotation_context(&mut feature, CONTEXT_LAYER_TRANSCRIPT);
            features.push(feature);
            exons_attached = exons_attached.saturating_add(1);
        }

        let mut cds_ranges = record
            .cds_1based
            .iter()
            .filter_map(|(start, end)| {
                let clipped_start = (*start).max(extracted_start_1based);
                let clipped_end = (*end).min(extracted_end_1based);
                (clipped_end >= clipped_start).then_some((clipped_start, clipped_end))
            })
            .collect::<Vec<_>>();
        cds_ranges.sort_unstable_by(|a, b| a.0.cmp(&b.0).then(a.1.cmp(&b.1)));
        cds_ranges.dedup();
        for (start, end) in cds_ranges {
            let Some(location) = Self::genomic_interval_to_local_location(
                extracted_start_1based,
                start,
                end,
                record.strand,
            ) else {
                continue;
            };
            let mut qualifiers = vec![
                ("transcript_id".into(), Some(record.transcript_id.clone())),
                ("chromosome".into(), Some(record.chromosome.clone())),
                ("genomic_start_1based".into(), Some(start.to_string())),
                ("genomic_end_1based".into(), Some(end.to_string())),
            ];
            if let Some(gene_name) = &record.gene_name {
                qualifiers.push(("gene".into(), Some(gene_name.clone())));
            }
            if let Some(gene_id) = &record.gene_id {
                qualifiers.push(("gene_id".into(), Some(gene_id.clone())));
            }
            if let Some(strand) = record.strand {
                qualifiers.push(("strand".into(), Some(strand.to_string())));
            }
            let mut feature = gb_io::seq::Feature {
                kind: "CDS".into(),
                location,
                qualifiers,
            };
            Self::mark_feature_as_genome_annotation_context(&mut feature, CONTEXT_LAYER_TRANSCRIPT);
            features.push(feature);
            cds_attached = cds_attached.saturating_add(1);
        }
        (features, exons_attached, cds_attached)
    }

    fn build_extract_region_annotation_projection(
        genes: &[GenomeGeneRecord],
        transcripts: &[GenomeTranscriptRecord],
        extracted_start_1based: usize,
        extracted_end_1based: usize,
        scope: GenomeAnnotationScope,
    ) -> ExtractRegionAnnotationProjectionBatch {
        let mut batch = ExtractRegionAnnotationProjectionBatch::default();
        if matches!(scope, GenomeAnnotationScope::None) {
            return batch;
        }
        for record in genes {
            if let Some(feature) = Self::gene_feature_from_genome_record(
                record,
                extracted_start_1based,
                extracted_end_1based,
            ) {
                batch.features.push(feature);
                batch.gene_count = batch.gene_count.saturating_add(1);
            }
        }
        for record in transcripts {
            if let Some(feature) = Self::transcript_feature_from_genome_record(
                record,
                extracted_start_1based,
                extracted_end_1based,
            ) {
                batch.features.push(feature);
                batch.transcript_count = batch.transcript_count.saturating_add(1);
            }
            if matches!(scope, GenomeAnnotationScope::Full) {
                let (subfeatures, exons_attached, cds_attached) =
                    Self::transcript_subfeatures_from_genome_record(
                        record,
                        extracted_start_1based,
                        extracted_end_1based,
                    );
                batch.exon_count = batch.exon_count.saturating_add(exons_attached);
                batch.cds_count = batch.cds_count.saturating_add(cds_attached);
                batch.features.extend(subfeatures);
            }
        }
        batch
    }

    fn build_exon_concatenated_projection(
        extracted_sequence: &str,
        extracted_start_1based: usize,
        extracted_end_1based: usize,
        strand: Option<char>,
        transcripts: &[GenomeTranscriptRecord],
        spacer_bp: usize,
    ) -> Option<ExonConcatenatedProjection> {
        if extracted_end_1based < extracted_start_1based || extracted_sequence.is_empty() {
            return None;
        }
        let mut clipped_exons_1based: Vec<(usize, usize)> = transcripts
            .iter()
            .flat_map(|record| record.exons_1based.iter())
            .filter_map(|(start, end)| {
                let clipped_start = (*start).max(extracted_start_1based);
                let clipped_end = (*end).min(extracted_end_1based);
                (clipped_end >= clipped_start).then_some((clipped_start, clipped_end))
            })
            .collect();
        if clipped_exons_1based.is_empty() {
            return None;
        }
        clipped_exons_1based.sort_unstable_by(|a, b| a.0.cmp(&b.0).then(a.1.cmp(&b.1)));
        let mut merged_exons_1based: Vec<(usize, usize)> = vec![];
        for (start, end) in clipped_exons_1based {
            if let Some(last) = merged_exons_1based.last_mut()
                && start <= last.1
            {
                last.1 = last.1.max(end);
            } else {
                merged_exons_1based.push((start, end));
            }
        }
        if merged_exons_1based.is_empty() {
            return None;
        }

        let sequence_bytes = extracted_sequence.as_bytes();
        let mut assembled: Vec<u8> = Vec::with_capacity(sequence_bytes.len());
        let mut blocks: Vec<ExonConcatenatedBlock> = vec![];
        for (index, (start, end)) in merged_exons_1based.iter().enumerate() {
            let local_start = start.saturating_sub(extracted_start_1based);
            let local_end_exclusive = end.saturating_sub(extracted_start_1based).saturating_add(1);
            if local_end_exclusive <= local_start || local_end_exclusive > sequence_bytes.len() {
                continue;
            }
            if index > 0 && spacer_bp > 0 {
                for _ in 0..spacer_bp {
                    assembled.push(b'N');
                }
            }
            let out_start = assembled.len();
            assembled.extend_from_slice(&sequence_bytes[local_start..local_end_exclusive]);
            let out_end = assembled.len();
            blocks.push(ExonConcatenatedBlock {
                genomic_start_1based: *start,
                genomic_end_1based: *end,
                local_start_0based: out_start,
                local_end_0based_exclusive: out_end,
            });
        }
        if blocks.is_empty() {
            return None;
        }

        let mut sequence = String::from_utf8(assembled).ok()?;
        if strand == Some('-') {
            let total_len = sequence.len();
            sequence = Self::reverse_complement(&sequence);
            for block in &mut blocks {
                let new_start = total_len.saturating_sub(block.local_end_0based_exclusive);
                let new_end = total_len.saturating_sub(block.local_start_0based);
                block.local_start_0based = new_start;
                block.local_end_0based_exclusive = new_end;
            }
            blocks.sort_unstable_by(|a, b| {
                a.local_start_0based.cmp(&b.local_start_0based).then(
                    a.local_end_0based_exclusive
                        .cmp(&b.local_end_0based_exclusive),
                )
            });
        }

        Some(ExonConcatenatedProjection { sequence, blocks })
    }

    fn transcript_feature_from_genome_record(
        record: &GenomeTranscriptRecord,
        extracted_start_1based: usize,
        extracted_end_1based: usize,
    ) -> Option<gb_io::seq::Feature> {
        let mut exons = record
            .exons_1based
            .iter()
            .filter_map(|(start, end)| {
                let clipped_start = (*start).max(extracted_start_1based);
                let clipped_end = (*end).min(extracted_end_1based);
                (clipped_end >= clipped_start).then_some((clipped_start, clipped_end))
            })
            .collect::<Vec<_>>();
        exons.sort_unstable_by(|a, b| a.0.cmp(&b.0).then(a.1.cmp(&b.1)));
        exons.dedup();
        if exons.is_empty() {
            return None;
        }
        let mut cds_local_1based = record
            .cds_1based
            .iter()
            .filter_map(|(start, end)| {
                let clipped_start = (*start).max(extracted_start_1based);
                let clipped_end = (*end).min(extracted_end_1based);
                if clipped_end < clipped_start {
                    return None;
                }
                let local_start = clipped_start
                    .saturating_sub(extracted_start_1based)
                    .saturating_add(1);
                let local_end = clipped_end
                    .saturating_sub(extracted_start_1based)
                    .saturating_add(1);
                (local_start > 0 && local_end >= local_start).then_some((local_start, local_end))
            })
            .collect::<Vec<_>>();
        cds_local_1based.sort_unstable_by(|a, b| a.0.cmp(&b.0).then(a.1.cmp(&b.1)));
        cds_local_1based.dedup();
        let mut parts = exons
            .iter()
            .filter_map(|(start, end)| {
                if *start < extracted_start_1based {
                    return None;
                }
                let local_start_0based = start.saturating_sub(extracted_start_1based);
                let local_end_exclusive =
                    end.saturating_sub(extracted_start_1based).saturating_add(1);
                (local_end_exclusive > local_start_0based).then_some(
                    gb_io::seq::Location::simple_range(
                        local_start_0based as i64,
                        local_end_exclusive as i64,
                    ),
                )
            })
            .collect::<Vec<_>>();
        if parts.is_empty() {
            return None;
        }
        let base_location = if parts.len() == 1 {
            parts.remove(0)
        } else {
            gb_io::seq::Location::Join(parts)
        };
        let location = if record.strand == Some('-') {
            gb_io::seq::Location::Complement(Box::new(base_location))
        } else {
            base_location
        };

        let mut qualifiers = vec![
            ("transcript_id".into(), Some(record.transcript_id.clone())),
            ("label".into(), Some(record.transcript_id.clone())),
            ("chromosome".into(), Some(record.chromosome.clone())),
            (
                "genomic_start_1based".into(),
                Some(record.transcript_start_1based.to_string()),
            ),
            (
                "genomic_end_1based".into(),
                Some(record.transcript_end_1based.to_string()),
            ),
        ];
        if let Some(gene_name) = &record.gene_name {
            qualifiers.push(("gene".into(), Some(gene_name.clone())));
        }
        if let Some(gene_id) = &record.gene_id {
            qualifiers.push(("gene_id".into(), Some(gene_id.clone())));
        }
        if let Some(strand) = record.strand {
            qualifiers.push(("strand".into(), Some(strand.to_string())));
        }
        if let Some(cds_encoded) = Self::serialize_ranges_1based(&cds_local_1based) {
            qualifiers.push(("cds_ranges_1based".into(), Some(cds_encoded)));
        }

        let mut feature = gb_io::seq::Feature {
            kind: "mRNA".into(),
            location,
            qualifiers,
        };
        Self::mark_feature_as_genome_annotation_context(&mut feature, CONTEXT_LAYER_TRANSCRIPT);
        Some(feature)
    }

    fn import_genome_slice_sequence(
        &mut self,
        result: &mut OpResult,
        sequence: String,
        default_id: String,
    ) -> Result<SeqId, EngineError> {
        let mut dna = DNAsequence::from_sequence(&sequence).map_err(|e| EngineError {
            code: ErrorCode::Internal,
            message: format!("Could not construct DNA sequence from genome slice: {e}"),
        })?;
        let seq_id = self.unique_seq_id(&default_id);
        dna.set_name(seq_id.clone());
        Self::prepare_sequence(&mut dna);
        self.state.sequences.insert(seq_id.clone(), dna);
        self.add_lineage_node(
            &seq_id,
            SequenceOrigin::ImportedGenomic,
            Some(&result.op_id),
        );
        result.created_seq_ids.push(seq_id.clone());
        Ok(seq_id)
    }

    fn now_unix_ms() -> u128 {
        std::time::SystemTime::now()
            .duration_since(std::time::UNIX_EPOCH)
            .map(|d| d.as_millis())
            .unwrap_or(0)
    }

    fn append_genome_extraction_provenance(&mut self, entry: GenomeExtractionProvenance) {
        let mut provenance = self
            .state
            .metadata
            .get(PROVENANCE_METADATA_KEY)
            .cloned()
            .unwrap_or_else(|| json!({}));
        if !provenance.is_object() {
            provenance = json!({});
        }
        if let Some(obj) = provenance.as_object_mut() {
            let mut records: Vec<GenomeExtractionProvenance> = obj
                .get(GENOME_EXTRACTIONS_METADATA_KEY)
                .cloned()
                .and_then(|v| serde_json::from_value(v).ok())
                .unwrap_or_default();
            records.push(entry);
            if let Ok(records_value) = serde_json::to_value(records) {
                obj.insert(GENOME_EXTRACTIONS_METADATA_KEY.to_string(), records_value);
            }
            obj.insert("updated_at_unix_ms".to_string(), json!(Self::now_unix_ms()));
            self.state
                .metadata
                .insert(PROVENANCE_METADATA_KEY.to_string(), provenance);
        }
    }

    fn genome_source_snapshot(
        source_plan: Option<&GenomeSourcePlan>,
        inspection: Option<&PreparedGenomeInspection>,
    ) -> (
        Option<String>,
        Option<String>,
        Option<String>,
        Option<String>,
        Option<String>,
        Option<String>,
    ) {
        let sequence_source_type = source_plan
            .map(|p| p.sequence_source_type.clone())
            .or_else(|| inspection.map(|i| i.sequence_source_type.clone()));
        let annotation_source_type = source_plan
            .map(|p| p.annotation_source_type.clone())
            .or_else(|| inspection.map(|i| i.annotation_source_type.clone()));
        let sequence_source = source_plan
            .map(|p| p.sequence_source.clone())
            .or_else(|| inspection.map(|i| i.sequence_source.clone()));
        let annotation_source = source_plan
            .map(|p| p.annotation_source.clone())
            .or_else(|| inspection.map(|i| i.annotation_source.clone()));
        let sequence_sha1 = inspection.and_then(|i| i.sequence_sha1.clone());
        let annotation_sha1 = inspection.and_then(|i| i.annotation_sha1.clone());
        (
            sequence_source_type,
            annotation_source_type,
            sequence_source,
            annotation_source,
            sequence_sha1,
            annotation_sha1,
        )
    }

    fn normalize_genome_track_subscription(
        mut subscription: GenomeTrackSubscription,
    ) -> Result<GenomeTrackSubscription, EngineError> {
        subscription.path = subscription.path.trim().to_string();
        if subscription.path.is_empty() {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: "Track path cannot be empty".to_string(),
            });
        }
        if subscription
            .min_score
            .zip(subscription.max_score)
            .map(|(min, max)| min > max)
            .unwrap_or(false)
        {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: "Track subscription requires min_score <= max_score".to_string(),
            });
        }
        subscription.track_name = subscription
            .track_name
            .map(|v| v.trim().to_string())
            .filter(|v| !v.is_empty());
        Ok(subscription)
    }

    fn sort_track_subscriptions(subscriptions: &mut Vec<GenomeTrackSubscription>) {
        subscriptions.sort_by(|a, b| {
            a.source
                .label()
                .cmp(b.source.label())
                .then(a.path.cmp(&b.path))
                .then(
                    a.track_name
                        .clone()
                        .unwrap_or_default()
                        .cmp(&b.track_name.clone().unwrap_or_default()),
                )
        });
    }

    fn read_track_subscriptions_from_metadata(
        value: Option<&serde_json::Value>,
    ) -> Vec<GenomeTrackSubscription> {
        let mut subscriptions = value
            .cloned()
            .and_then(|v| serde_json::from_value::<Vec<GenomeTrackSubscription>>(v).ok())
            .unwrap_or_default();
        subscriptions.retain(|subscription| !subscription.path.trim().is_empty());
        Self::sort_track_subscriptions(&mut subscriptions);
        subscriptions
    }

    fn write_track_subscriptions_to_metadata(
        &mut self,
        subscriptions: &[GenomeTrackSubscription],
    ) -> Result<(), EngineError> {
        if subscriptions.is_empty() {
            self.state
                .metadata
                .remove(GENOME_TRACK_SUBSCRIPTIONS_METADATA_KEY);
            return Ok(());
        }
        let value = serde_json::to_value(subscriptions).map_err(|e| EngineError {
            code: ErrorCode::Internal,
            message: format!("Could not serialize genome track subscriptions: {e}"),
        })?;
        self.state
            .metadata
            .insert(GENOME_TRACK_SUBSCRIPTIONS_METADATA_KEY.to_string(), value);
        Ok(())
    }

    fn read_known_track_anchor_ids(&self) -> BTreeSet<String> {
        self.state
            .metadata
            .get(GENOME_TRACK_KNOWN_ANCHORS_METADATA_KEY)
            .and_then(|v| v.as_array())
            .map(|values| {
                values
                    .iter()
                    .filter_map(|v| v.as_str())
                    .map(str::trim)
                    .filter(|v| !v.is_empty())
                    .map(str::to_string)
                    .collect::<BTreeSet<_>>()
            })
            .unwrap_or_default()
    }

    fn write_known_track_anchor_ids(&mut self, anchors: &BTreeSet<String>) {
        if anchors.is_empty() {
            self.state
                .metadata
                .remove(GENOME_TRACK_KNOWN_ANCHORS_METADATA_KEY);
            return;
        }
        let values = anchors
            .iter()
            .map(|v| serde_json::Value::String(v.clone()))
            .collect::<Vec<_>>();
        let new_value = serde_json::Value::Array(values);
        if self
            .state
            .metadata
            .get(GENOME_TRACK_KNOWN_ANCHORS_METADATA_KEY)
            == Some(&new_value)
        {
            return;
        }
        self.state.metadata.insert(
            GENOME_TRACK_KNOWN_ANCHORS_METADATA_KEY.to_string(),
            new_value,
        );
    }

    fn track_subscription_to_operation(
        seq_id: &str,
        subscription: &GenomeTrackSubscription,
    ) -> Operation {
        match subscription.source {
            GenomeTrackSource::Bed => Operation::ImportGenomeBedTrack {
                seq_id: seq_id.to_string(),
                path: subscription.path.clone(),
                track_name: subscription.track_name.clone(),
                min_score: subscription.min_score,
                max_score: subscription.max_score,
                clear_existing: Some(subscription.clear_existing),
            },
            GenomeTrackSource::BigWig => Operation::ImportGenomeBigWigTrack {
                seq_id: seq_id.to_string(),
                path: subscription.path.clone(),
                track_name: subscription.track_name.clone(),
                min_score: subscription.min_score,
                max_score: subscription.max_score,
                clear_existing: Some(subscription.clear_existing),
            },
            GenomeTrackSource::Vcf => Operation::ImportGenomeVcfTrack {
                seq_id: seq_id.to_string(),
                path: subscription.path.clone(),
                track_name: subscription.track_name.clone(),
                min_score: subscription.min_score,
                max_score: subscription.max_score,
                clear_existing: Some(subscription.clear_existing),
            },
        }
    }

    fn apply_track_subscription_to_seq_ids(
        &mut self,
        seq_ids: &[String],
        subscription: &GenomeTrackSubscription,
    ) -> GenomeTrackSyncReport {
        let mut report = GenomeTrackSyncReport {
            subscriptions_considered: 1,
            target_sequences: seq_ids.len(),
            ..GenomeTrackSyncReport::default()
        };
        for seq_id in seq_ids {
            let op = Self::track_subscription_to_operation(seq_id, subscription);
            match self.apply(op) {
                Ok(op_result) => {
                    report.applied_imports += 1;
                    report.warnings_count += op_result.warnings.len();
                }
                Err(e) => {
                    report.failed_imports += 1;
                    if report.errors.len() < 20 {
                        report.errors.push(format!(
                            "{} '{}' @ {}: {}",
                            subscription.source.label(),
                            subscription.path,
                            seq_id,
                            e.message
                        ));
                    }
                }
            }
        }
        report
    }

    fn normalize_candidate_set_name(raw: &str) -> Result<String, EngineError> {
        let trimmed = raw.trim();
        if trimmed.is_empty() {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: "Candidate set name cannot be empty".to_string(),
            });
        }
        Ok(trimmed.to_string())
    }

    fn normalize_metric_name(raw: &str) -> String {
        let mut out = String::with_capacity(raw.len());
        for ch in raw.chars() {
            if ch.is_ascii_alphanumeric() {
                out.push(ch.to_ascii_lowercase());
            } else {
                out.push('_');
            }
        }
        let trimmed = out.trim_matches('_');
        if trimmed.is_empty() {
            "metric".to_string()
        } else {
            trimmed.to_string()
        }
    }

    fn read_candidate_store_from_metadata(value: Option<&serde_json::Value>) -> CandidateStore {
        let mut store = value
            .cloned()
            .and_then(|v| serde_json::from_value::<CandidateStore>(v).ok())
            .unwrap_or_default();
        if store.schema.trim().is_empty() {
            store.schema = CANDIDATE_SETS_SCHEMA.to_string();
        }
        store
    }

    fn read_candidate_store(&self) -> CandidateStore {
        Self::read_candidate_store_from_metadata(
            self.state.metadata.get(CANDIDATE_SETS_METADATA_KEY),
        )
    }

    fn write_candidate_store(&mut self, mut store: CandidateStore) -> Result<(), EngineError> {
        if store.sets.is_empty() {
            self.state.metadata.remove(CANDIDATE_SETS_METADATA_KEY);
            return Ok(());
        }
        store.schema = CANDIDATE_SETS_SCHEMA.to_string();
        store.updated_at_unix_ms = Self::now_unix_ms();
        let value = serde_json::to_value(store).map_err(|e| EngineError {
            code: ErrorCode::Internal,
            message: format!("Could not serialize candidate-set metadata: {e}"),
        })?;
        self.state
            .metadata
            .insert(CANDIDATE_SETS_METADATA_KEY.to_string(), value);
        Ok(())
    }

    fn read_isoform_panel_store_from_metadata(
        value: Option<&serde_json::Value>,
    ) -> IsoformPanelStore {
        let mut store = value
            .cloned()
            .and_then(|v| serde_json::from_value::<IsoformPanelStore>(v).ok())
            .unwrap_or_default();
        if store.schema.trim().is_empty() {
            store.schema = ISOFORM_PANELS_SCHEMA.to_string();
        }
        store
    }

    fn read_isoform_panel_store(&self) -> IsoformPanelStore {
        Self::read_isoform_panel_store_from_metadata(
            self.state.metadata.get(ISOFORM_PANELS_METADATA_KEY),
        )
    }

    fn write_isoform_panel_store(
        &mut self,
        mut store: IsoformPanelStore,
    ) -> Result<(), EngineError> {
        if store.records.is_empty() {
            self.state.metadata.remove(ISOFORM_PANELS_METADATA_KEY);
            return Ok(());
        }
        store.schema = ISOFORM_PANELS_SCHEMA.to_string();
        store.updated_at_unix_ms = Self::now_unix_ms();
        let value = serde_json::to_value(store).map_err(|e| EngineError {
            code: ErrorCode::Internal,
            message: format!("Could not serialize isoform-panel metadata: {e}"),
        })?;
        self.state
            .metadata
            .insert(ISOFORM_PANELS_METADATA_KEY.to_string(), value);
        Ok(())
    }

    fn upsert_isoform_panel_record(
        &mut self,
        record: IsoformPanelRecord,
    ) -> Result<(), EngineError> {
        let mut store = self.read_isoform_panel_store();
        store.records.retain(|existing| {
            !(existing.seq_id == record.seq_id && existing.panel_id == record.panel_id)
        });
        store.records.push(record);
        self.write_isoform_panel_store(store)
    }

    fn get_isoform_panel_record(
        &self,
        seq_id: &str,
        panel_id: &str,
    ) -> Result<IsoformPanelRecord, EngineError> {
        let probe = panel_id.trim();
        if probe.is_empty() {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: "panel_id cannot be empty".to_string(),
            });
        }
        let store = self.read_isoform_panel_store();
        store
            .records
            .iter()
            .rev()
            .find(|record| record.seq_id == seq_id && record.panel_id == probe)
            .cloned()
            .ok_or_else(|| EngineError {
                code: ErrorCode::NotFound,
                message: format!(
                    "Isoform panel '{}' was not imported for sequence '{}'",
                    probe, seq_id
                ),
            })
    }

    fn read_primer_design_store_from_metadata(
        value: Option<&serde_json::Value>,
    ) -> PrimerDesignStore {
        let mut store = value
            .cloned()
            .and_then(|v| serde_json::from_value::<PrimerDesignStore>(v).ok())
            .unwrap_or_default();
        if store.schema.trim().is_empty() {
            store.schema = PRIMER_DESIGN_REPORTS_SCHEMA.to_string();
        }
        store
    }

    fn read_primer_design_store(&self) -> PrimerDesignStore {
        Self::read_primer_design_store_from_metadata(
            self.state.metadata.get(PRIMER_DESIGN_REPORTS_METADATA_KEY),
        )
    }

    fn write_primer_design_store(
        &mut self,
        mut store: PrimerDesignStore,
    ) -> Result<(), EngineError> {
        if store.reports.is_empty() && store.qpcr_reports.is_empty() {
            self.state
                .metadata
                .remove(PRIMER_DESIGN_REPORTS_METADATA_KEY);
            return Ok(());
        }
        store.schema = PRIMER_DESIGN_REPORTS_SCHEMA.to_string();
        store.updated_at_unix_ms = Self::now_unix_ms();
        let value = serde_json::to_value(store).map_err(|e| EngineError {
            code: ErrorCode::Internal,
            message: format!("Could not serialize primer-design metadata: {e}"),
        })?;
        self.state
            .metadata
            .insert(PRIMER_DESIGN_REPORTS_METADATA_KEY.to_string(), value);
        Ok(())
    }

    fn normalize_primer_design_report_id(raw: &str) -> Result<String, EngineError> {
        let trimmed = raw.trim();
        if trimmed.is_empty() {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: "report_id cannot be empty".to_string(),
            });
        }
        let mut out = String::with_capacity(trimmed.len());
        for ch in trimmed.chars() {
            if ch.is_ascii_alphanumeric() || matches!(ch, '_' | '-' | '.') {
                out.push(ch);
            } else {
                out.push('_');
            }
        }
        if out.is_empty() {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: "report_id must contain at least one ASCII letter, digit, '-', '_' or '.'"
                    .to_string(),
            });
        }
        Ok(out)
    }

    pub fn list_primer_design_reports(&self) -> Vec<PrimerDesignReportSummary> {
        let store = self.read_primer_design_store();
        let mut ids = store.reports.keys().cloned().collect::<Vec<_>>();
        ids.sort();
        ids.into_iter()
            .filter_map(|id| store.reports.get(&id))
            .map(|report| PrimerDesignReportSummary {
                report_id: report.report_id.clone(),
                template: report.template.clone(),
                generated_at_unix_ms: report.generated_at_unix_ms,
                roi_start_0based: report.roi_start_0based,
                roi_end_0based: report.roi_end_0based,
                pair_count: report.pair_count,
            })
            .collect()
    }

    pub fn get_primer_design_report(
        &self,
        report_id: &str,
    ) -> Result<PrimerDesignReport, EngineError> {
        let report_id = Self::normalize_primer_design_report_id(report_id)?;
        let store = self.read_primer_design_store();
        store
            .reports
            .get(&report_id)
            .cloned()
            .ok_or_else(|| EngineError {
                code: ErrorCode::NotFound,
                message: format!("Primer-design report '{}' not found", report_id),
            })
    }

    pub fn export_primer_design_report(
        &self,
        report_id: &str,
        path: &str,
    ) -> Result<PrimerDesignReport, EngineError> {
        let report = self.get_primer_design_report(report_id)?;
        let text = serde_json::to_string_pretty(&report).map_err(|e| EngineError {
            code: ErrorCode::Internal,
            message: format!(
                "Could not serialize primer-design report '{}': {e}",
                report.report_id
            ),
        })?;
        std::fs::write(path, text).map_err(|e| EngineError {
            code: ErrorCode::Io,
            message: format!("Could not write primer-design report to '{path}': {e}"),
        })?;
        Ok(report)
    }

    pub fn list_qpcr_design_reports(&self) -> Vec<QpcrDesignReportSummary> {
        let store = self.read_primer_design_store();
        let mut ids = store.qpcr_reports.keys().cloned().collect::<Vec<_>>();
        ids.sort();
        ids.into_iter()
            .filter_map(|id| store.qpcr_reports.get(&id))
            .map(|report| QpcrDesignReportSummary {
                report_id: report.report_id.clone(),
                template: report.template.clone(),
                generated_at_unix_ms: report.generated_at_unix_ms,
                roi_start_0based: report.roi_start_0based,
                roi_end_0based: report.roi_end_0based,
                assay_count: report.assay_count,
            })
            .collect()
    }

    pub fn get_qpcr_design_report(&self, report_id: &str) -> Result<QpcrDesignReport, EngineError> {
        let report_id = Self::normalize_primer_design_report_id(report_id)?;
        let store = self.read_primer_design_store();
        store
            .qpcr_reports
            .get(&report_id)
            .cloned()
            .ok_or_else(|| EngineError {
                code: ErrorCode::NotFound,
                message: format!("qPCR design report '{}' not found", report_id),
            })
    }

    pub fn export_qpcr_design_report(
        &self,
        report_id: &str,
        path: &str,
    ) -> Result<QpcrDesignReport, EngineError> {
        let report = self.get_qpcr_design_report(report_id)?;
        let text = serde_json::to_string_pretty(&report).map_err(|e| EngineError {
            code: ErrorCode::Internal,
            message: format!(
                "Could not serialize qPCR design report '{}': {e}",
                report.report_id
            ),
        })?;
        std::fs::write(path, text).map_err(|e| EngineError {
            code: ErrorCode::Io,
            message: format!("Could not write qPCR design report to '{path}': {e}"),
        })?;
        Ok(report)
    }

    fn feature_query_label_values(feature: &gb_io::seq::Feature) -> Vec<String> {
        let mut labels: Vec<String> = vec![];
        let mut seen: HashSet<String> = HashSet::new();
        for key in [
            "label",
            "name",
            "standard_name",
            "gene",
            "locus_tag",
            "product",
            "note",
        ] {
            for value in feature.qualifier_values(key) {
                let trimmed = value.trim();
                if trimmed.is_empty() {
                    continue;
                }
                let normalized = trimmed.to_ascii_uppercase();
                if seen.insert(normalized) {
                    labels.push(trimmed.to_string());
                }
            }
        }
        labels
    }

    fn normalize_feature_kind_token(raw: &str) -> Option<String> {
        let trimmed = raw.trim();
        if trimmed.is_empty() {
            None
        } else {
            Some(trimmed.to_ascii_uppercase())
        }
    }

    pub fn query_sequence_features(
        &self,
        mut query: SequenceFeatureQuery,
    ) -> Result<SequenceFeatureQueryResult, EngineError> {
        query.seq_id = query.seq_id.trim().to_string();
        if query.seq_id.is_empty() {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: "query_sequence_features requires non-empty seq_id".to_string(),
            });
        }
        query.kind_in = query
            .kind_in
            .iter()
            .filter_map(|value| Self::normalize_feature_kind_token(value))
            .collect::<BTreeSet<_>>()
            .into_iter()
            .collect();
        query.kind_not_in = query
            .kind_not_in
            .iter()
            .filter_map(|value| Self::normalize_feature_kind_token(value))
            .collect::<BTreeSet<_>>()
            .into_iter()
            .collect();
        query.label_contains = query
            .label_contains
            .as_deref()
            .map(str::trim)
            .filter(|value| !value.is_empty())
            .map(|value| value.to_string());
        query.label_regex = query
            .label_regex
            .as_deref()
            .map(str::trim)
            .filter(|value| !value.is_empty())
            .map(|value| value.to_string());

        let label_contains_upper = query
            .label_contains
            .as_ref()
            .map(|value| value.to_ascii_uppercase());
        let label_regex_compiled = query
            .label_regex
            .as_ref()
            .map(|pattern| {
                RegexBuilder::new(pattern)
                    .case_insensitive(true)
                    .build()
                    .map_err(|e| EngineError {
                        code: ErrorCode::InvalidInput,
                        message: format!("Invalid label_regex '{}': {e}", pattern),
                    })
            })
            .transpose()?;

        if let (Some(min_len), Some(max_len)) = (query.min_len_bp, query.max_len_bp)
            && min_len > max_len
        {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: format!(
                    "Invalid feature length filter: min_len_bp ({min_len}) must be <= max_len_bp ({max_len})"
                ),
            });
        }

        if query.start_0based.is_some() != query.end_0based_exclusive.is_some() {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: "Feature range filter requires both start_0based and end_0based_exclusive"
                    .to_string(),
            });
        }

        let dna = self
            .state
            .sequences
            .get(&query.seq_id)
            .ok_or_else(|| EngineError {
                code: ErrorCode::NotFound,
                message: format!("Sequence '{}' not found", query.seq_id),
            })?;
        let sequence_length_bp = dna.len();

        let range_filter = if let (Some(start), Some(end_exclusive)) =
            (query.start_0based, query.end_0based_exclusive)
        {
            if sequence_length_bp == 0 {
                return Err(EngineError {
                    code: ErrorCode::InvalidInput,
                    message: "Feature range filter cannot be applied on empty sequence".to_string(),
                });
            }
            if end_exclusive <= start {
                return Err(EngineError {
                    code: ErrorCode::InvalidInput,
                    message: format!(
                        "Invalid feature range filter: start ({start}) must be < end ({end_exclusive})"
                    ),
                });
            }
            if start >= sequence_length_bp || end_exclusive > sequence_length_bp {
                return Err(EngineError {
                    code: ErrorCode::InvalidInput,
                    message: format!(
                        "Feature range filter {}..{} is outside sequence length {}",
                        start, end_exclusive, sequence_length_bp
                    ),
                });
            }
            Some((start, end_exclusive))
        } else {
            None
        };

        let mut qualifier_filters: Vec<(SequenceFeatureQualifierFilter, Option<Regex>)> = vec![];
        for mut filter in query.qualifier_filters.clone() {
            filter.key = filter.key.trim().to_string();
            if filter.key.is_empty() {
                return Err(EngineError {
                    code: ErrorCode::InvalidInput,
                    message: "Qualifier filter key must not be empty".to_string(),
                });
            }
            filter.value_contains = filter
                .value_contains
                .as_deref()
                .map(str::trim)
                .filter(|value| !value.is_empty())
                .map(|value| value.to_string());
            filter.value_regex = filter
                .value_regex
                .as_deref()
                .map(str::trim)
                .filter(|value| !value.is_empty())
                .map(|value| value.to_string());
            let compiled = filter
                .value_regex
                .as_ref()
                .map(|pattern| {
                    RegexBuilder::new(pattern)
                        .case_insensitive(!filter.case_sensitive)
                        .build()
                        .map_err(|e| EngineError {
                            code: ErrorCode::InvalidInput,
                            message: format!(
                                "Invalid qualifier regex for key '{}': {e}",
                                filter.key
                            ),
                        })
                })
                .transpose()?;
            qualifier_filters.push((filter, compiled));
        }
        query.qualifier_filters = qualifier_filters
            .iter()
            .map(|(filter, _)| filter.clone())
            .collect();

        let mut rows: Vec<SequenceFeatureQueryRow> = vec![];
        let total_feature_count = dna.features().len();
        let kind_in = query.kind_in.iter().collect::<HashSet<_>>();
        let kind_not_in = query.kind_not_in.iter().collect::<HashSet<_>>();

        for (feature_id, feature) in dna.features().iter().enumerate() {
            let kind = feature.kind.to_string();
            let kind_upper = kind.to_ascii_uppercase();
            if !query.include_source && kind_upper == "SOURCE" {
                continue;
            }
            if !kind_in.is_empty() && !kind_in.contains(&kind_upper) {
                continue;
            }
            if kind_not_in.contains(&kind_upper) {
                continue;
            }

            let mut ranges: Vec<(usize, usize)> = vec![];
            collect_location_ranges_usize(&feature.location, &mut ranges);
            if ranges.is_empty() {
                let Ok((from, to)) = feature.location.find_bounds() else {
                    continue;
                };
                if from < 0 || to < 0 {
                    continue;
                }
                ranges.push((from as usize, to as usize));
            }
            let Some(start_0based) = ranges.iter().map(|(start, _)| *start).min() else {
                continue;
            };
            let Some(end_0based_exclusive) = ranges
                .iter()
                .map(|(_, end)| *end)
                .max()
                .map(|end| end.min(sequence_length_bp))
            else {
                continue;
            };
            if end_0based_exclusive <= start_0based || start_0based >= sequence_length_bp {
                continue;
            }
            let length_bp = end_0based_exclusive.saturating_sub(start_0based);

            if let Some(min_len_bp) = query.min_len_bp
                && length_bp < min_len_bp
            {
                continue;
            }
            if let Some(max_len_bp) = query.max_len_bp
                && length_bp > max_len_bp
            {
                continue;
            }

            if let Some((range_start, range_end)) = range_filter {
                let range_ok = match query.range_relation {
                    SequenceFeatureRangeRelation::Overlap => {
                        end_0based_exclusive > range_start && start_0based < range_end
                    }
                    SequenceFeatureRangeRelation::Within => {
                        start_0based >= range_start && end_0based_exclusive <= range_end
                    }
                    SequenceFeatureRangeRelation::Contains => {
                        start_0based <= range_start && end_0based_exclusive >= range_end
                    }
                };
                if !range_ok {
                    continue;
                }
            }

            let is_reverse = feature_is_reverse(feature);
            let strand_text = if is_reverse { "reverse" } else { "forward" };
            let strand_ok = match query.strand {
                SequenceFeatureStrandFilter::Any => true,
                SequenceFeatureStrandFilter::Forward => !is_reverse,
                SequenceFeatureStrandFilter::Reverse => is_reverse,
            };
            if !strand_ok {
                continue;
            }

            let labels = Self::feature_query_label_values(feature);
            if let Some(needle_upper) = label_contains_upper.as_ref() {
                let label_match = labels.iter().any(|value| {
                    let upper = value.to_ascii_uppercase();
                    upper == *needle_upper || upper.contains(needle_upper)
                });
                if !label_match {
                    continue;
                }
            }
            if let Some(label_regex) = label_regex_compiled.as_ref() {
                let label_match = labels.iter().any(|value| label_regex.is_match(value));
                if !label_match {
                    continue;
                }
            }

            let mut qualifiers_map: BTreeMap<String, Vec<String>> = BTreeMap::new();
            if query.include_qualifiers {
                for (key, value) in &feature.qualifiers {
                    let Some(value) = value.as_deref().map(str::trim).filter(|v| !v.is_empty())
                    else {
                        continue;
                    };
                    qualifiers_map
                        .entry(key.to_string())
                        .or_default()
                        .push(value.to_string());
                }
            }

            let qualifiers_ok = qualifier_filters.iter().all(|(filter, regex)| {
                let values = feature
                    .qualifier_values(filter.key.as_str())
                    .map(|value| value.trim().to_string())
                    .filter(|value| !value.is_empty())
                    .collect::<Vec<_>>();
                if values.is_empty() {
                    return false;
                }
                if let Some(needle) = filter.value_contains.as_ref() {
                    if filter.case_sensitive {
                        if !values.iter().any(|value| value.contains(needle)) {
                            return false;
                        }
                    } else {
                        let needle_upper = needle.to_ascii_uppercase();
                        if !values
                            .iter()
                            .any(|value| value.to_ascii_uppercase().contains(&needle_upper))
                        {
                            return false;
                        }
                    }
                }
                if let Some(regex) = regex {
                    if !values.iter().any(|value| regex.is_match(value)) {
                        return false;
                    }
                }
                true
            });
            if !qualifiers_ok {
                continue;
            }

            rows.push(SequenceFeatureQueryRow {
                feature_id,
                kind,
                start_0based,
                end_0based_exclusive,
                length_bp,
                strand: strand_text.to_string(),
                label: Self::feature_display_label(feature, feature_id),
                labels,
                qualifiers: qualifiers_map,
            });
        }

        rows.sort_by(|a, b| match query.sort_by {
            SequenceFeatureSortBy::FeatureId => a.feature_id.cmp(&b.feature_id),
            SequenceFeatureSortBy::Start => a
                .start_0based
                .cmp(&b.start_0based)
                .then(a.end_0based_exclusive.cmp(&b.end_0based_exclusive))
                .then(a.feature_id.cmp(&b.feature_id)),
            SequenceFeatureSortBy::End => a
                .end_0based_exclusive
                .cmp(&b.end_0based_exclusive)
                .then(a.start_0based.cmp(&b.start_0based))
                .then(a.feature_id.cmp(&b.feature_id)),
            SequenceFeatureSortBy::Kind => a
                .kind
                .to_ascii_uppercase()
                .cmp(&b.kind.to_ascii_uppercase())
                .then(a.start_0based.cmp(&b.start_0based))
                .then(a.end_0based_exclusive.cmp(&b.end_0based_exclusive))
                .then(a.feature_id.cmp(&b.feature_id)),
            SequenceFeatureSortBy::Length => a
                .length_bp
                .cmp(&b.length_bp)
                .then(a.start_0based.cmp(&b.start_0based))
                .then(a.end_0based_exclusive.cmp(&b.end_0based_exclusive))
                .then(a.feature_id.cmp(&b.feature_id)),
        });
        if query.descending {
            rows.reverse();
        }

        let matched_count = rows.len();
        let offset = query.offset.min(matched_count);
        let limit = query
            .limit
            .unwrap_or(FEATURE_QUERY_DEFAULT_LIMIT)
            .clamp(1, FEATURE_QUERY_MAX_LIMIT);
        let rows = rows
            .into_iter()
            .skip(offset)
            .take(limit)
            .collect::<Vec<_>>();
        let returned_count = rows.len();
        query.limit = Some(limit);
        query.offset = offset;

        Ok(SequenceFeatureQueryResult {
            schema: FEATURE_QUERY_RESULT_SCHEMA.to_string(),
            seq_id: query.seq_id.clone(),
            sequence_length_bp,
            total_feature_count,
            matched_count,
            returned_count,
            offset,
            limit,
            range_relation: query.range_relation.as_str().to_string(),
            strand_filter: query.strand.as_str().to_string(),
            sort_by: query.sort_by.as_str().to_string(),
            descending: query.descending,
            query,
            rows,
        })
    }

    pub(crate) fn normalize_planning_class_key(raw: &str) -> String {
        let mut out = String::with_capacity(raw.len());
        for ch in raw.trim().chars() {
            if ch.is_ascii_alphanumeric() {
                out.push(ch.to_ascii_lowercase());
            } else if matches!(ch, '_' | '-' | '.' | ':') {
                out.push('_');
            }
        }
        let trimmed = out.trim_matches('_');
        if trimmed.is_empty() {
            "item".to_string()
        } else {
            trimmed.to_string()
        }
    }

    fn validate_planning_profile_schema(profile: &PlanningProfile) -> Result<(), EngineError> {
        let schema = profile.schema.trim();
        if !schema.is_empty() && schema != PLANNING_PROFILE_SCHEMA {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: format!(
                    "Unsupported planning profile schema '{}' (expected '{}')",
                    schema, PLANNING_PROFILE_SCHEMA
                ),
            });
        }
        Ok(())
    }

    fn validate_planning_objective_schema(
        objective: &PlanningObjective,
    ) -> Result<(), EngineError> {
        let schema = objective.schema.trim();
        if !schema.is_empty() && schema != PLANNING_OBJECTIVE_SCHEMA {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: format!(
                    "Unsupported planning objective schema '{}' (expected '{}')",
                    schema, PLANNING_OBJECTIVE_SCHEMA
                ),
            });
        }
        Ok(())
    }

    fn normalize_planning_profile(mut profile: PlanningProfile) -> PlanningProfile {
        profile.schema = PLANNING_PROFILE_SCHEMA.to_string();
        profile.profile_id = profile
            .profile_id
            .and_then(|v| (!v.trim().is_empty()).then(|| v.trim().to_string()));
        profile.currency = profile
            .currency
            .and_then(|v| (!v.trim().is_empty()).then(|| v.trim().to_ascii_uppercase()));
        if !profile.procurement_business_days_default.is_finite()
            || profile.procurement_business_days_default <= 0.0
        {
            profile.procurement_business_days_default = 10.0;
        }
        profile.notes = profile
            .notes
            .and_then(|v| (!v.trim().is_empty()).then(|| v.trim().to_string()));

        let mut normalized_caps = profile
            .capabilities
            .iter()
            .map(|value| Self::normalize_planning_class_key(value))
            .filter(|value| !value.is_empty())
            .collect::<BTreeSet<_>>()
            .into_iter()
            .collect::<Vec<_>>();
        normalized_caps.sort_by_key(|value| value.to_ascii_lowercase());
        profile.capabilities = normalized_caps;

        let mut inventory: HashMap<String, PlanningInventoryItem> = HashMap::new();
        for (key, mut item) in profile.inventory {
            let normalized = Self::normalize_planning_class_key(key.as_str());
            item.procurement_business_days = item
                .procurement_business_days
                .filter(|value| value.is_finite() && *value > 0.0);
            item.unit_cost = item
                .unit_cost
                .filter(|value| value.is_finite() && *value >= 0.0);
            item.note = item
                .note
                .and_then(|v| (!v.trim().is_empty()).then(|| v.trim().to_string()));
            inventory.insert(normalized, item);
        }
        profile.inventory = inventory;

        let mut machines: HashMap<String, PlanningMachineAvailability> = HashMap::new();
        for (key, mut machine) in profile.machine_availability {
            let normalized = Self::normalize_planning_class_key(key.as_str());
            if !machine.queue_business_days.is_finite() || machine.queue_business_days < 0.0 {
                machine.queue_business_days = 0.0;
            }
            machine.run_cost_per_hour = machine
                .run_cost_per_hour
                .filter(|value| value.is_finite() && *value >= 0.0);
            machine.note = machine
                .note
                .and_then(|v| (!v.trim().is_empty()).then(|| v.trim().to_string()));
            machines.insert(normalized, machine);
        }
        profile.machine_availability = machines;
        profile
    }

    fn normalize_planning_objective(mut objective: PlanningObjective) -> PlanningObjective {
        objective.schema = PLANNING_OBJECTIVE_SCHEMA.to_string();
        if !objective.weight_time.is_finite() || objective.weight_time < 0.0 {
            objective.weight_time = 0.0;
        }
        if !objective.weight_cost.is_finite() || objective.weight_cost < 0.0 {
            objective.weight_cost = 0.0;
        }
        if !objective.weight_local_fit.is_finite() || objective.weight_local_fit < 0.0 {
            objective.weight_local_fit = 0.0;
        }
        if objective.weight_time == 0.0
            && objective.weight_cost == 0.0
            && objective.weight_local_fit == 0.0
        {
            objective.weight_time = 1.0;
            objective.weight_cost = 1.0;
            objective.weight_local_fit = 1.0;
        }
        objective.max_cost = objective
            .max_cost
            .filter(|value| value.is_finite() && *value >= 0.0);
        objective.max_time_hours = objective
            .max_time_hours
            .filter(|value| value.is_finite() && *value >= 0.0);
        let mut required_caps = objective
            .required_capabilities
            .iter()
            .map(|value| Self::normalize_planning_class_key(value))
            .filter(|value| !value.is_empty())
            .collect::<BTreeSet<_>>()
            .into_iter()
            .collect::<Vec<_>>();
        required_caps.sort_by_key(|value| value.to_ascii_lowercase());
        objective.required_capabilities = required_caps;
        objective
    }

    fn merge_planning_profile(base: &mut PlanningProfile, overlay: &PlanningProfile) {
        if let Some(profile_id) = overlay
            .profile_id
            .as_deref()
            .map(str::trim)
            .filter(|value| !value.is_empty())
        {
            base.profile_id = Some(profile_id.to_string());
        }
        if let Some(currency) = overlay
            .currency
            .as_deref()
            .map(str::trim)
            .filter(|value| !value.is_empty())
        {
            base.currency = Some(currency.to_ascii_uppercase());
        }
        if overlay.procurement_business_days_default.is_finite()
            && overlay.procurement_business_days_default > 0.0
        {
            base.procurement_business_days_default = overlay.procurement_business_days_default;
        }
        for capability in &overlay.capabilities {
            let key = Self::normalize_planning_class_key(capability);
            if !key.is_empty() && !base.capabilities.iter().any(|existing| existing == &key) {
                base.capabilities.push(key);
            }
        }
        for (key, item) in &overlay.inventory {
            let normalized = Self::normalize_planning_class_key(key.as_str());
            base.inventory.insert(normalized, item.clone());
        }
        for (key, machine) in &overlay.machine_availability {
            let normalized = Self::normalize_planning_class_key(key.as_str());
            base.machine_availability
                .insert(normalized, machine.clone());
        }
        if let Some(notes) = overlay
            .notes
            .as_deref()
            .map(str::trim)
            .filter(|value| !value.is_empty())
        {
            base.notes = Some(notes.to_string());
        }
        base.capabilities
            .sort_by_key(|value| value.to_ascii_lowercase());
        base.capabilities.dedup();
        base.schema = PLANNING_PROFILE_SCHEMA.to_string();
    }

    fn read_planning_store_from_metadata(value: Option<&serde_json::Value>) -> PlanningStore {
        let mut store = value
            .cloned()
            .and_then(|v| serde_json::from_value::<PlanningStore>(v).ok())
            .unwrap_or_default();
        if store.schema.trim().is_empty() {
            store.schema = PLANNING_STORE_SCHEMA.to_string();
        }
        store.global_profile = store
            .global_profile
            .take()
            .map(Self::normalize_planning_profile);
        store.project_override_profile = store
            .project_override_profile
            .take()
            .map(Self::normalize_planning_profile);
        store.confirmed_agent_overlay_profile = store
            .confirmed_agent_overlay_profile
            .take()
            .map(Self::normalize_planning_profile);
        store.objective = store
            .objective
            .take()
            .map(Self::normalize_planning_objective);

        for suggestion in store.suggestions.values_mut() {
            suggestion.schema = PLANNING_SUGGESTION_SCHEMA.to_string();
            suggestion.source = suggestion.source.trim().to_string();
            suggestion.direction = suggestion.direction.trim().to_ascii_lowercase();
            if !matches!(suggestion.direction.as_str(), "pull" | "push") {
                suggestion.direction = "pull".to_string();
            }
            suggestion.profile_patch = suggestion
                .profile_patch
                .take()
                .map(Self::normalize_planning_profile);
            suggestion.objective_patch = suggestion
                .objective_patch
                .take()
                .map(Self::normalize_planning_objective);
            suggestion.confidence = suggestion
                .confidence
                .filter(|value| value.is_finite() && *value >= 0.0 && *value <= 1.0);
        }
        store.sync_status.schema = PLANNING_SYNC_STATUS_SCHEMA.to_string();
        store.sync_status.pending_suggestion_count = store
            .suggestions
            .values()
            .filter(|row| row.status == PlanningSuggestionStatus::Pending)
            .count();
        if store.next_suggestion_counter == 0 {
            store.next_suggestion_counter = 1;
        }
        store
    }

    fn read_planning_store(&self) -> PlanningStore {
        Self::read_planning_store_from_metadata(self.state.metadata.get(PLANNING_METADATA_KEY))
    }

    fn write_planning_store(&mut self, mut store: PlanningStore) -> Result<(), EngineError> {
        store.schema = PLANNING_STORE_SCHEMA.to_string();
        store.updated_at_unix_ms = Self::now_unix_ms();
        store.sync_status.schema = PLANNING_SYNC_STATUS_SCHEMA.to_string();
        store.sync_status.pending_suggestion_count = store
            .suggestions
            .values()
            .filter(|row| row.status == PlanningSuggestionStatus::Pending)
            .count();
        if store.next_suggestion_counter == 0 {
            store.next_suggestion_counter = 1;
        }

        let is_empty = store.global_profile.is_none()
            && store.project_override_profile.is_none()
            && store.confirmed_agent_overlay_profile.is_none()
            && store.objective.is_none()
            && store.suggestions.is_empty()
            && store.sync_status.pending_suggestion_count == 0
            && store.sync_status.last_pull_at_unix_ms.is_none()
            && store.sync_status.last_push_at_unix_ms.is_none()
            && store.sync_status.last_source.is_none()
            && store.sync_status.last_snapshot_id.is_none()
            && store.sync_status.last_error.is_none();
        if is_empty {
            self.state.metadata.remove(PLANNING_METADATA_KEY);
            return Ok(());
        }

        let value = serde_json::to_value(store).map_err(|e| EngineError {
            code: ErrorCode::Internal,
            message: format!("Could not serialize planning metadata: {e}"),
        })?;
        self.state
            .metadata
            .insert(PLANNING_METADATA_KEY.to_string(), value);
        Ok(())
    }

    pub fn planning_profile(&self, scope: PlanningProfileScope) -> Option<PlanningProfile> {
        let store = self.read_planning_store();
        match scope {
            PlanningProfileScope::Global => store.global_profile,
            PlanningProfileScope::ProjectOverride => store.project_override_profile,
            PlanningProfileScope::ConfirmedAgentOverlay => store.confirmed_agent_overlay_profile,
            PlanningProfileScope::Effective => Some(self.planning_effective_profile()),
        }
    }

    pub fn planning_effective_profile(&self) -> PlanningProfile {
        let store = self.read_planning_store();
        let mut effective = PlanningProfile::default();
        if let Some(global) = store.global_profile {
            Self::merge_planning_profile(&mut effective, &global);
        }
        if let Some(overlay) = store.confirmed_agent_overlay_profile {
            Self::merge_planning_profile(&mut effective, &overlay);
        }
        if let Some(project_override) = store.project_override_profile {
            Self::merge_planning_profile(&mut effective, &project_override);
        }
        effective.schema = PLANNING_PROFILE_SCHEMA.to_string();
        effective
    }

    pub fn set_planning_profile(
        &mut self,
        scope: PlanningProfileScope,
        profile: Option<PlanningProfile>,
    ) -> Result<(), EngineError> {
        if scope == PlanningProfileScope::Effective {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: "Cannot set profile scope 'effective' directly".to_string(),
            });
        }
        let mut store = self.read_planning_store();
        if let Some(profile) = profile.as_ref() {
            Self::validate_planning_profile_schema(profile)?;
        }
        let normalized = profile.map(Self::normalize_planning_profile);
        match scope {
            PlanningProfileScope::Global => store.global_profile = normalized,
            PlanningProfileScope::ProjectOverride => store.project_override_profile = normalized,
            PlanningProfileScope::ConfirmedAgentOverlay => {
                store.confirmed_agent_overlay_profile = normalized
            }
            PlanningProfileScope::Effective => {}
        }
        self.write_planning_store(store)
    }

    pub fn planning_objective(&self) -> PlanningObjective {
        self.read_planning_store()
            .objective
            .map(Self::normalize_planning_objective)
            .unwrap_or_default()
    }

    pub fn set_planning_objective(
        &mut self,
        objective: Option<PlanningObjective>,
    ) -> Result<(), EngineError> {
        let mut store = self.read_planning_store();
        if let Some(objective) = objective.as_ref() {
            Self::validate_planning_objective_schema(objective)?;
        }
        store.objective = objective.map(Self::normalize_planning_objective);
        self.write_planning_store(store)
    }

    pub fn planning_meta_enabled(&self) -> bool {
        let store = self.read_planning_store();
        store.global_profile.is_some()
            || store.project_override_profile.is_some()
            || store.confirmed_agent_overlay_profile.is_some()
            || store.objective.is_some()
    }

    fn build_planning_diff(
        current_profile: Option<&PlanningProfile>,
        profile_patch: Option<&PlanningProfile>,
        current_objective: Option<&PlanningObjective>,
        objective_patch: Option<&PlanningObjective>,
    ) -> serde_json::Value {
        let mut out = serde_json::Map::<String, serde_json::Value>::new();
        if let Some(patch) = profile_patch {
            let current_value = current_profile
                .and_then(|value| serde_json::to_value(value).ok())
                .unwrap_or_else(|| json!({}));
            let patch_value = serde_json::to_value(patch).unwrap_or_else(|_| json!({}));
            out.insert("profile_current".to_string(), current_value);
            out.insert("profile_patch".to_string(), patch_value);
        }
        if let Some(patch) = objective_patch {
            let current_value = current_objective
                .and_then(|value| serde_json::to_value(value).ok())
                .unwrap_or_else(|| json!({}));
            let patch_value = serde_json::to_value(patch).unwrap_or_else(|_| json!({}));
            out.insert("objective_current".to_string(), current_value);
            out.insert("objective_patch".to_string(), patch_value);
        }
        serde_json::Value::Object(out)
    }

    pub fn propose_planning_suggestion(
        &mut self,
        direction: &str,
        source: &str,
        confidence: Option<f64>,
        snapshot_id: Option<&str>,
        profile_patch: Option<PlanningProfile>,
        objective_patch: Option<PlanningObjective>,
        message: Option<&str>,
    ) -> Result<PlanningSuggestion, EngineError> {
        let direction = direction.trim().to_ascii_lowercase();
        if !matches!(direction.as_str(), "pull" | "push") {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: format!("Unsupported planning sync direction '{direction}'"),
            });
        }
        if profile_patch.is_none() && objective_patch.is_none() {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message:
                    "Planning suggestion requires at least one of profile_patch/objective_patch"
                        .to_string(),
            });
        }
        let source = source.trim();
        if source.is_empty() {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: "Planning suggestion source cannot be empty".to_string(),
            });
        }

        let mut store = self.read_planning_store();
        let counter = store.next_suggestion_counter.max(1);
        store.next_suggestion_counter = counter.saturating_add(1);
        let suggestion_id = format!("planning_suggestion_{counter:06}");

        if let Some(profile_patch) = profile_patch.as_ref() {
            Self::validate_planning_profile_schema(profile_patch)?;
        }
        if let Some(objective_patch) = objective_patch.as_ref() {
            Self::validate_planning_objective_schema(objective_patch)?;
        }
        let normalized_profile_patch = profile_patch.map(Self::normalize_planning_profile);
        let normalized_objective_patch = objective_patch.map(Self::normalize_planning_objective);
        let diff = Self::build_planning_diff(
            store.confirmed_agent_overlay_profile.as_ref(),
            normalized_profile_patch.as_ref(),
            store.objective.as_ref(),
            normalized_objective_patch.as_ref(),
        );
        let now = Self::now_unix_ms();
        let clean_message = message
            .map(str::trim)
            .filter(|value| !value.is_empty())
            .map(|value| value.to_string());
        let confidence =
            confidence.filter(|value| value.is_finite() && *value >= 0.0 && *value <= 1.0);
        let snapshot_id = snapshot_id
            .map(str::trim)
            .filter(|value| !value.is_empty())
            .map(|value| value.to_string());

        let suggestion = PlanningSuggestion {
            schema: PLANNING_SUGGESTION_SCHEMA.to_string(),
            suggestion_id: suggestion_id.clone(),
            status: PlanningSuggestionStatus::Pending,
            direction: direction.clone(),
            source: source.to_string(),
            confidence,
            snapshot_id: snapshot_id.clone(),
            message: clean_message,
            profile_patch: normalized_profile_patch,
            objective_patch: normalized_objective_patch,
            diff,
            created_at_unix_ms: now,
            resolved_at_unix_ms: None,
            rejection_reason: None,
        };
        if direction == "pull" {
            store.sync_status.last_pull_at_unix_ms = Some(now);
        } else {
            store.sync_status.last_push_at_unix_ms = Some(now);
        }
        store.sync_status.last_source = Some(source.to_string());
        store.sync_status.last_snapshot_id = snapshot_id;
        store.sync_status.last_error = None;
        store.suggestions.insert(suggestion_id, suggestion.clone());
        self.write_planning_store(store)?;
        Ok(suggestion)
    }

    pub fn list_planning_suggestions(
        &self,
        status: Option<PlanningSuggestionStatus>,
    ) -> Vec<PlanningSuggestion> {
        let store = self.read_planning_store();
        let mut rows = store
            .suggestions
            .values()
            .filter(|row| status.is_none_or(|probe| row.status == probe))
            .cloned()
            .collect::<Vec<_>>();
        rows.sort_by(|left, right| {
            left.created_at_unix_ms
                .cmp(&right.created_at_unix_ms)
                .then(left.suggestion_id.cmp(&right.suggestion_id))
        });
        rows
    }

    pub fn accept_planning_suggestion(
        &mut self,
        suggestion_id: &str,
    ) -> Result<PlanningSuggestion, EngineError> {
        let target = suggestion_id.trim();
        if target.is_empty() {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: "suggestion_id cannot be empty".to_string(),
            });
        }
        let mut store = self.read_planning_store();
        let Some(mut suggestion) = store.suggestions.get(target).cloned() else {
            return Err(EngineError {
                code: ErrorCode::NotFound,
                message: format!("Planning suggestion '{}' not found", target),
            });
        };
        if suggestion.status != PlanningSuggestionStatus::Pending {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: format!(
                    "Planning suggestion '{}' is already {}",
                    target,
                    suggestion.status.as_str()
                ),
            });
        }
        if let Some(profile_patch) = suggestion.profile_patch.clone() {
            let mut overlay = store
                .confirmed_agent_overlay_profile
                .clone()
                .unwrap_or_default();
            Self::merge_planning_profile(&mut overlay, &profile_patch);
            store.confirmed_agent_overlay_profile = Some(Self::normalize_planning_profile(overlay));
        }
        if let Some(objective_patch) = suggestion.objective_patch.clone() {
            store.objective = Some(Self::normalize_planning_objective(objective_patch));
        }
        suggestion.status = PlanningSuggestionStatus::Accepted;
        suggestion.resolved_at_unix_ms = Some(Self::now_unix_ms());
        suggestion.rejection_reason = None;
        store
            .suggestions
            .insert(target.to_string(), suggestion.clone());
        self.write_planning_store(store)?;
        Ok(suggestion)
    }

    pub fn reject_planning_suggestion(
        &mut self,
        suggestion_id: &str,
        reason: Option<&str>,
    ) -> Result<PlanningSuggestion, EngineError> {
        let target = suggestion_id.trim();
        if target.is_empty() {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: "suggestion_id cannot be empty".to_string(),
            });
        }
        let mut store = self.read_planning_store();
        let Some(mut suggestion) = store.suggestions.get(target).cloned() else {
            return Err(EngineError {
                code: ErrorCode::NotFound,
                message: format!("Planning suggestion '{}' not found", target),
            });
        };
        if suggestion.status != PlanningSuggestionStatus::Pending {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: format!(
                    "Planning suggestion '{}' is already {}",
                    target,
                    suggestion.status.as_str()
                ),
            });
        }
        suggestion.status = PlanningSuggestionStatus::Rejected;
        suggestion.resolved_at_unix_ms = Some(Self::now_unix_ms());
        suggestion.rejection_reason = reason
            .map(str::trim)
            .filter(|value| !value.is_empty())
            .map(|value| value.to_string());
        store
            .suggestions
            .insert(target.to_string(), suggestion.clone());
        self.write_planning_store(store)?;
        Ok(suggestion)
    }

    pub fn planning_sync_status(&self) -> PlanningSyncStatus {
        self.read_planning_store().sync_status
    }

    pub fn mark_planning_sync_error(
        &mut self,
        direction: Option<&str>,
        source: Option<&str>,
        snapshot_id: Option<&str>,
        error: &str,
    ) -> Result<PlanningSyncStatus, EngineError> {
        let clean_error = error.trim();
        if clean_error.is_empty() {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: "Planning sync error message cannot be empty".to_string(),
            });
        }
        let mut store = self.read_planning_store();
        let now = Self::now_unix_ms();
        if let Some(direction) = direction
            .map(str::trim)
            .filter(|value| !value.is_empty())
            .map(|value| value.to_ascii_lowercase())
        {
            match direction.as_str() {
                "pull" => store.sync_status.last_pull_at_unix_ms = Some(now),
                "push" => store.sync_status.last_push_at_unix_ms = Some(now),
                _ => {}
            }
        }
        if let Some(source) = source
            .map(str::trim)
            .filter(|value| !value.is_empty())
            .map(|value| value.to_string())
        {
            store.sync_status.last_source = Some(source);
        }
        if let Some(snapshot_id) = snapshot_id
            .map(str::trim)
            .filter(|value| !value.is_empty())
            .map(|value| value.to_string())
        {
            store.sync_status.last_snapshot_id = Some(snapshot_id);
        }
        store.sync_status.last_error = Some(clean_error.to_string());
        let status = store.sync_status.clone();
        self.write_planning_store(store)?;
        Ok(status)
    }

    fn read_dotplot_analysis_store_from_metadata(
        value: Option<&serde_json::Value>,
    ) -> DotplotAnalysisStore {
        let mut store = value
            .cloned()
            .and_then(|v| serde_json::from_value::<DotplotAnalysisStore>(v).ok())
            .unwrap_or_default();
        if store.schema.trim().is_empty() {
            store.schema = DOTPLOT_ANALYSIS_SCHEMA.to_string();
        }
        for view in store.dotplots.values_mut() {
            view.normalize_v3_defaults();
            view.schema = DOTPLOT_VIEW_SCHEMA.to_string();
        }
        store
    }

    fn read_dotplot_analysis_store(&self) -> DotplotAnalysisStore {
        Self::read_dotplot_analysis_store_from_metadata(
            self.state.metadata.get(DOTPLOT_ANALYSIS_METADATA_KEY),
        )
    }

    fn write_dotplot_analysis_store(
        &mut self,
        mut store: DotplotAnalysisStore,
    ) -> Result<(), EngineError> {
        if store.dotplots.is_empty() && store.flexibility_tracks.is_empty() {
            self.state.metadata.remove(DOTPLOT_ANALYSIS_METADATA_KEY);
            return Ok(());
        }
        store.schema = DOTPLOT_ANALYSIS_SCHEMA.to_string();
        store.updated_at_unix_ms = Self::now_unix_ms();
        let value = serde_json::to_value(store).map_err(|e| EngineError {
            code: ErrorCode::Internal,
            message: format!("Could not serialize dotplot analysis metadata: {e}"),
        })?;
        self.state
            .metadata
            .insert(DOTPLOT_ANALYSIS_METADATA_KEY.to_string(), value);
        Ok(())
    }

    fn normalize_analysis_id(raw: &str, kind: &str) -> Result<String, EngineError> {
        let trimmed = raw.trim();
        if trimmed.is_empty() {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: format!("{kind}_id cannot be empty"),
            });
        }
        let mut out = String::with_capacity(trimmed.len());
        for ch in trimmed.chars() {
            if ch.is_ascii_alphanumeric() || matches!(ch, '_' | '-' | '.') {
                out.push(ch);
            } else {
                out.push('_');
            }
        }
        if out.is_empty() {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: format!(
                    "{kind}_id must contain at least one ASCII letter, digit, '-', '_' or '.'"
                ),
            });
        }
        Ok(out)
    }

    fn resolve_analysis_span(
        seq_len: usize,
        span_start_0based: Option<usize>,
        span_end_0based: Option<usize>,
    ) -> Result<(usize, usize), EngineError> {
        let start = span_start_0based.unwrap_or(0);
        let end = span_end_0based.unwrap_or(seq_len);
        if start >= seq_len {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: format!(
                    "span_start_0based ({start}) must be within sequence length ({seq_len})"
                ),
            });
        }
        if end > seq_len {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: format!("span_end_0based ({end}) must be <= sequence length ({seq_len})"),
            });
        }
        if end <= start {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: format!(
                    "Invalid span {start}..{end}; span_end_0based must be > span_start_0based"
                ),
            });
        }
        Ok((start, end))
    }

    fn count_mismatches_capped(left: &[u8], right: &[u8], cap: usize) -> usize {
        let mut mismatches = 0usize;
        for (a, b) in left.iter().zip(right.iter()) {
            if !a.eq_ignore_ascii_case(b) {
                mismatches += 1;
                if mismatches > cap {
                    break;
                }
            }
        }
        mismatches
    }

    fn default_dotplot_series_color(index: usize) -> [u8; 3] {
        const PALETTE: [[u8; 3]; 8] = [
            [29, 78, 216],
            [220, 38, 38],
            [5, 150, 105],
            [217, 119, 6],
            [124, 58, 237],
            [190, 24, 93],
            [8, 145, 178],
            [71, 85, 105],
        ];
        PALETTE[index % PALETTE.len()]
    }

    fn build_dotplot_query_series(
        series_id: String,
        seq_id: String,
        label: String,
        color_rgb: [u8; 3],
        mode: DotplotMode,
        span_start_0based: usize,
        span_end_0based: usize,
        points: Vec<DotplotMatchPoint>,
        boxplot_bins: Vec<DotplotBoxplotBin>,
    ) -> DotplotQuerySeries {
        DotplotQuerySeries {
            series_id,
            seq_id,
            label,
            color_rgb,
            mode,
            span_start_0based,
            span_end_0based,
            point_count: points.len(),
            points,
            boxplot_bin_count: boxplot_bins.len(),
            boxplot_bins,
        }
    }

    fn build_dotplot_reference_annotation_track(
        &self,
        reference_seq_id: &str,
        reference_span_start_0based: usize,
        reference_span_end_0based: usize,
    ) -> Option<DotplotReferenceAnnotationTrack> {
        let dna = self.state.sequences.get(reference_seq_id)?;
        let mut intervals: Vec<(usize, usize)> = vec![];
        for feature in dna.features() {
            if !feature.kind.to_string().eq_ignore_ascii_case("exon") {
                continue;
            }
            let mut ranges: Vec<(usize, usize)> = vec![];
            collect_location_ranges_usize(&feature.location, &mut ranges);
            for (start_0based, end_0based_exclusive) in ranges {
                let clipped_start = start_0based.max(reference_span_start_0based);
                let clipped_end = end_0based_exclusive.min(reference_span_end_0based);
                if clipped_end > clipped_start {
                    intervals.push((clipped_start, clipped_end));
                }
            }
        }
        if intervals.is_empty() {
            return None;
        }
        intervals.sort_unstable_by(|a, b| a.0.cmp(&b.0).then(a.1.cmp(&b.1)));
        let mut merged: Vec<(usize, usize)> = vec![];
        for (start, end) in intervals {
            if let Some(last) = merged.last_mut()
                && start <= last.1
            {
                last.1 = last.1.max(end);
            } else {
                merged.push((start, end));
            }
        }
        let intervals = merged
            .into_iter()
            .map(
                |(start_0based, end_0based_exclusive)| DotplotReferenceAnnotationInterval {
                    start_0based,
                    end_0based_exclusive,
                    label: "exon".to_string(),
                },
            )
            .collect::<Vec<_>>();
        Some(DotplotReferenceAnnotationTrack {
            seq_id: reference_seq_id.to_string(),
            label: "merged exons".to_string(),
            interval_count: intervals.len(),
            intervals,
        })
    }

    fn dotplot_boxplot_quantile(sorted: &[usize], q: f64) -> Option<usize> {
        if sorted.is_empty() {
            return None;
        }
        let q = q.clamp(0.0, 1.0);
        let idx = ((sorted.len().saturating_sub(1)) as f64 * q).round() as usize;
        sorted.get(idx.min(sorted.len().saturating_sub(1))).copied()
    }

    fn compute_dotplot_boxplot_bins(
        points: &[DotplotMatchPoint],
        query_span_start_0based: usize,
        query_span_end_0based: usize,
        requested_bins: usize,
    ) -> Vec<DotplotBoxplotBin> {
        if points.is_empty() || query_span_end_0based <= query_span_start_0based {
            return vec![];
        }
        let query_span = query_span_end_0based.saturating_sub(query_span_start_0based);
        let bin_count = requested_bins.clamp(1, 512).min(query_span.max(1));
        let bin_width = query_span.div_ceil(bin_count);
        if bin_width == 0 {
            return vec![];
        }
        let mut bins: Vec<Vec<usize>> = vec![vec![]; bin_count];
        for point in points {
            if point.x_0based < query_span_start_0based || point.x_0based >= query_span_end_0based {
                continue;
            }
            let local = point.x_0based.saturating_sub(query_span_start_0based);
            let idx = local
                .checked_div(bin_width)
                .unwrap_or(0)
                .min(bin_count.saturating_sub(1));
            bins[idx].push(point.y_0based);
        }
        let mut out = Vec::with_capacity(bin_count);
        for (idx, mut values) in bins.into_iter().enumerate() {
            let start = query_span_start_0based.saturating_add(idx.saturating_mul(bin_width));
            let end = (start + bin_width).min(query_span_end_0based);
            if values.is_empty() {
                out.push(DotplotBoxplotBin {
                    query_start_0based: start,
                    query_end_0based_exclusive: end,
                    hit_count: 0,
                    min_reference_0based: None,
                    q1_reference_0based: None,
                    median_reference_0based: None,
                    q3_reference_0based: None,
                    max_reference_0based: None,
                });
                continue;
            }
            values.sort_unstable();
            let hit_count = values.len();
            out.push(DotplotBoxplotBin {
                query_start_0based: start,
                query_end_0based_exclusive: end,
                hit_count,
                min_reference_0based: values.first().copied(),
                q1_reference_0based: Self::dotplot_boxplot_quantile(&values, 0.25),
                median_reference_0based: Self::dotplot_boxplot_quantile(&values, 0.50),
                q3_reference_0based: Self::dotplot_boxplot_quantile(&values, 0.75),
                max_reference_0based: values.last().copied(),
            });
        }
        out
    }

    fn compute_dotplot_points(
        query_sequence: &[u8],
        reference_sequence: &[u8],
        query_span_start_0based: usize,
        reference_span_start_0based: usize,
        mode: DotplotMode,
        word_size: usize,
        step_bp: usize,
        max_mismatches: usize,
        max_points: usize,
    ) -> Result<(Vec<DotplotMatchPoint>, bool), EngineError> {
        if word_size == 0 {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: "ComputeDotplot requires word_size >= 1".to_string(),
            });
        }
        if step_bp == 0 {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: "ComputeDotplot requires step_bp >= 1".to_string(),
            });
        }
        if query_sequence.len() < word_size {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: format!(
                    "ComputeDotplot word_size ({word_size}) exceeds selected query span length ({})",
                    query_sequence.len()
                ),
            });
        }
        if reference_sequence.len() < word_size {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: format!(
                    "ComputeDotplot word_size ({word_size}) exceeds selected reference span length ({})",
                    reference_sequence.len()
                ),
            });
        }
        let query_positions: Vec<usize> = (0..=query_sequence.len() - word_size)
            .step_by(step_bp)
            .collect();
        let reference_positions: Vec<usize> = (0..=reference_sequence.len() - word_size)
            .step_by(step_bp)
            .collect();

        if max_mismatches == 0 {
            let mut reference_index: HashMap<Vec<u8>, Vec<usize>> = HashMap::new();
            for y_start in &reference_positions {
                let key: Vec<u8> = match mode {
                    DotplotMode::SelfForward | DotplotMode::PairForward => reference_sequence
                        [*y_start..*y_start + word_size]
                        .iter()
                        .map(|b| b.to_ascii_uppercase())
                        .collect(),
                    DotplotMode::SelfReverseComplement | DotplotMode::PairReverseComplement => {
                        Self::reverse_complement_bytes(
                            &reference_sequence[*y_start..*y_start + word_size],
                        )
                        .into_iter()
                        .map(|b| b.to_ascii_uppercase())
                        .collect()
                    }
                };
                reference_index.entry(key).or_default().push(*y_start);
            }
            let mut points: Vec<DotplotMatchPoint> = vec![];
            let mut truncated = false;
            for x_start in &query_positions {
                let left_key: Vec<u8> = query_sequence[*x_start..*x_start + word_size]
                    .iter()
                    .map(|b| b.to_ascii_uppercase())
                    .collect();
                if let Some(y_hits) = reference_index.get(&left_key) {
                    for y_start in y_hits {
                        points.push(DotplotMatchPoint {
                            x_0based: query_span_start_0based + *x_start,
                            y_0based: reference_span_start_0based + *y_start,
                            mismatches: 0,
                        });
                        if points.len() >= max_points {
                            truncated = true;
                            return Ok((points, truncated));
                        }
                    }
                }
            }
            return Ok((points, truncated));
        }

        let pair_evaluations = query_positions
            .len()
            .saturating_mul(reference_positions.len())
            .max(query_positions.len());
        if pair_evaluations > MAX_DOTPLOT_PAIR_EVALUATIONS {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: format!(
                    "ComputeDotplot requested {} pair evaluations (limit {}); increase step_bp or reduce span",
                    pair_evaluations, MAX_DOTPLOT_PAIR_EVALUATIONS
                ),
            });
        }

        let mut reverse_words: Vec<Vec<u8>> = vec![];
        if matches!(
            mode,
            DotplotMode::SelfReverseComplement | DotplotMode::PairReverseComplement
        ) {
            reverse_words = reference_positions
                .iter()
                .map(|start| {
                    Self::reverse_complement_bytes(&reference_sequence[*start..*start + word_size])
                        .into_iter()
                        .map(|b| b.to_ascii_uppercase())
                        .collect::<Vec<_>>()
                })
                .collect();
        }

        let mut points: Vec<DotplotMatchPoint> = vec![];
        let mut truncated = false;
        for (x_idx, x_start) in query_positions.iter().enumerate() {
            let left = &query_sequence[*x_start..*x_start + word_size];
            for (y_idx, y_start) in reference_positions.iter().enumerate() {
                let mismatches = match mode {
                    DotplotMode::SelfForward | DotplotMode::PairForward => {
                        Self::count_mismatches_capped(
                            left,
                            &reference_sequence[*y_start..*y_start + word_size],
                            max_mismatches,
                        )
                    }
                    DotplotMode::SelfReverseComplement | DotplotMode::PairReverseComplement => {
                        Self::count_mismatches_capped(left, &reverse_words[y_idx], max_mismatches)
                    }
                };
                if mismatches <= max_mismatches {
                    points.push(DotplotMatchPoint {
                        x_0based: query_span_start_0based + query_positions[x_idx],
                        y_0based: reference_span_start_0based + reference_positions[y_idx],
                        mismatches,
                    });
                    if points.len() >= max_points {
                        truncated = true;
                        return Ok((points, truncated));
                    }
                }
            }
        }
        Ok((points, truncated))
    }

    /// Build a transient pairwise dotplot view without storing it in engine state.
    ///
    /// This reuses the same deterministic point-generation logic as
    /// `Operation::ComputeDotplot`, but returns the resulting view directly so
    /// UI surfaces can inspect a read/reference comparison inline before
    /// deciding whether to export or persist it.
    pub fn preview_pair_dotplot_view(
        query_seq_id: &str,
        query_text: &str,
        reference_seq_id: &str,
        reference_text: &str,
        reference_span_start_0based: usize,
        reference_span_end_0based: usize,
        mode: DotplotMode,
        word_size: usize,
        step_bp: usize,
        max_mismatches: usize,
        tile_bp: Option<usize>,
    ) -> Result<DotplotView, EngineError> {
        let query_text = query_text.trim().to_ascii_uppercase();
        let reference_text = reference_text.trim().to_ascii_uppercase();
        let query_bytes = query_text.as_bytes();
        let reference_bytes = reference_text.as_bytes();
        let (reference_span_start_0based, reference_span_end_0based) = Self::resolve_analysis_span(
            reference_bytes.len(),
            Some(reference_span_start_0based),
            Some(reference_span_end_0based),
        )?;
        let (points, _truncated) = Self::compute_dotplot_points(
            query_bytes,
            &reference_bytes[reference_span_start_0based..reference_span_end_0based],
            0,
            reference_span_start_0based,
            mode,
            word_size,
            step_bp,
            max_mismatches,
            MAX_DOTPLOT_POINTS,
        )?;
        let boxplot_bins = Self::compute_dotplot_boxplot_bins(
            &points,
            0,
            query_bytes.len(),
            DOTPLOT_BOXPLOT_DEFAULT_BINS,
        );
        let query_label = query_seq_id.trim().to_string();
        let primary_series = Self::build_dotplot_query_series(
            "preview_series_1".to_string(),
            query_seq_id.trim().to_string(),
            query_label.clone(),
            Self::default_dotplot_series_color(0),
            mode,
            0,
            query_bytes.len(),
            points,
            boxplot_bins,
        );
        Ok(DotplotView {
            schema: DOTPLOT_VIEW_SCHEMA.to_string(),
            dotplot_id: String::new(),
            owner_seq_id: query_seq_id.trim().to_string(),
            seq_id: query_seq_id.trim().to_string(),
            reference_seq_id: Some(reference_seq_id.trim().to_string()),
            generated_at_unix_ms: 0,
            span_start_0based: primary_series.span_start_0based,
            span_end_0based: primary_series.span_end_0based,
            reference_span_start_0based,
            reference_span_end_0based,
            mode,
            word_size,
            step_bp,
            max_mismatches,
            tile_bp,
            point_count: primary_series.point_count,
            points: primary_series.points.clone(),
            boxplot_bin_count: primary_series.boxplot_bin_count,
            boxplot_bins: primary_series.boxplot_bins.clone(),
            series_count: 1,
            query_series: vec![primary_series],
            reference_annotation: None,
        })
    }

    fn compute_flexibility_track_bins(
        sequence: &[u8],
        span_start_0based: usize,
        model: FlexibilityModel,
        bin_bp: usize,
        smoothing_bp: Option<usize>,
    ) -> Result<Vec<FlexibilityBinScore>, EngineError> {
        if bin_bp == 0 {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: "ComputeFlexibilityTrack requires bin_bp >= 1".to_string(),
            });
        }
        if sequence.is_empty() {
            return Ok(vec![]);
        }
        let mut bins: Vec<FlexibilityBinScore> = vec![];
        let mut idx = 0usize;
        while idx < sequence.len() {
            let end = (idx + bin_bp).min(sequence.len());
            let slice = &sequence[idx..end];
            let mut a_count = 0usize;
            let mut t_count = 0usize;
            for b in slice {
                match b.to_ascii_uppercase() {
                    b'A' => a_count += 1,
                    b'T' => t_count += 1,
                    _ => {}
                }
            }
            let at_total = a_count + t_count;
            let score = match model {
                FlexibilityModel::AtRichness => {
                    if slice.is_empty() {
                        0.0
                    } else {
                        at_total as f64 / slice.len() as f64
                    }
                }
                FlexibilityModel::AtSkew => {
                    if at_total == 0 {
                        0.0
                    } else {
                        (a_count as f64 - t_count as f64) / at_total as f64
                    }
                }
            };
            bins.push(FlexibilityBinScore {
                start_0based: span_start_0based + idx,
                end_0based_exclusive: span_start_0based + end,
                score,
            });
            idx = end;
        }
        if bins.is_empty() {
            return Ok(bins);
        }
        let smoothing = smoothing_bp
            .and_then(|value| (value > 0).then_some(value))
            .map(|value| value.div_ceil(bin_bp))
            .unwrap_or(0);
        if smoothing > 1 {
            let radius = smoothing / 2;
            let raw = bins.iter().map(|bin| bin.score).collect::<Vec<_>>();
            for (index, bin) in bins.iter_mut().enumerate() {
                let start = index.saturating_sub(radius);
                let end = (index + radius + 1).min(raw.len());
                if end > start {
                    let sum = raw[start..end].iter().sum::<f64>();
                    bin.score = sum / (end - start) as f64;
                }
            }
        }
        Ok(bins)
    }

    fn upsert_dotplot_view(&mut self, view: DotplotView) -> Result<(), EngineError> {
        let mut store = self.read_dotplot_analysis_store();
        let mut view = view;
        view.schema = DOTPLOT_VIEW_SCHEMA.to_string();
        view.normalize_v3_defaults();
        store.dotplots.insert(view.dotplot_id.clone(), view);
        self.write_dotplot_analysis_store(store)
    }

    fn upsert_flexibility_track(&mut self, track: FlexibilityTrack) -> Result<(), EngineError> {
        let mut store = self.read_dotplot_analysis_store();
        store
            .flexibility_tracks
            .insert(track.track_id.clone(), track);
        self.write_dotplot_analysis_store(store)
    }

    pub fn list_dotplot_views(&self, seq_id_filter: Option<&str>) -> Vec<DotplotViewSummary> {
        let filter = seq_id_filter
            .map(str::trim)
            .filter(|value| !value.is_empty())
            .map(|value| value.to_ascii_lowercase());
        let mut rows = self
            .read_dotplot_analysis_store()
            .dotplots
            .values()
            .filter(|view| {
                filter.as_ref().is_none_or(|needle| {
                    view.owner_seq_id
                        .to_ascii_lowercase()
                        .eq_ignore_ascii_case(needle)
                })
            })
            .map(|view| DotplotViewSummary {
                dotplot_id: view.dotplot_id.clone(),
                owner_seq_id: view.owner_seq_id.clone(),
                seq_id: view.seq_id.clone(),
                reference_seq_id: view.reference_seq_id.clone(),
                generated_at_unix_ms: view.generated_at_unix_ms,
                span_start_0based: view.span_start_0based,
                span_end_0based: view.span_end_0based,
                reference_span_start_0based: view.reference_span_start_0based,
                reference_span_end_0based: view.reference_span_end_0based,
                mode: view.mode,
                word_size: view.word_size,
                step_bp: view.step_bp,
                max_mismatches: view.max_mismatches,
                point_count: view.point_count,
                series_count: view.series_count.max(view.query_series.len()),
            })
            .collect::<Vec<_>>();
        rows.sort_by(|left, right| {
            left.owner_seq_id
                .to_ascii_lowercase()
                .cmp(&right.owner_seq_id.to_ascii_lowercase())
                .then(left.generated_at_unix_ms.cmp(&right.generated_at_unix_ms))
                .then(
                    left.dotplot_id
                        .to_ascii_lowercase()
                        .cmp(&right.dotplot_id.to_ascii_lowercase()),
                )
        });
        rows
    }

    pub fn get_dotplot_view(&self, dotplot_id: &str) -> Result<DotplotView, EngineError> {
        let dotplot_id = Self::normalize_analysis_id(dotplot_id, "dotplot")?;
        self.read_dotplot_analysis_store()
            .dotplots
            .get(dotplot_id.as_str())
            .cloned()
            .ok_or_else(|| EngineError {
                code: ErrorCode::NotFound,
                message: format!("Dotplot view '{}' not found", dotplot_id),
            })
    }

    pub fn list_flexibility_tracks(
        &self,
        seq_id_filter: Option<&str>,
    ) -> Vec<FlexibilityTrackSummary> {
        let filter = seq_id_filter
            .map(str::trim)
            .filter(|value| !value.is_empty())
            .map(|value| value.to_ascii_lowercase());
        let mut rows = self
            .read_dotplot_analysis_store()
            .flexibility_tracks
            .values()
            .filter(|track| {
                filter.as_ref().is_none_or(|needle| {
                    track
                        .seq_id
                        .to_ascii_lowercase()
                        .eq_ignore_ascii_case(needle)
                })
            })
            .map(|track| FlexibilityTrackSummary {
                track_id: track.track_id.clone(),
                seq_id: track.seq_id.clone(),
                generated_at_unix_ms: track.generated_at_unix_ms,
                span_start_0based: track.span_start_0based,
                span_end_0based: track.span_end_0based,
                model: track.model,
                bin_bp: track.bin_bp,
                smoothing_bp: track.smoothing_bp,
                bin_count: track.bins.len(),
                min_score: track.min_score,
                max_score: track.max_score,
            })
            .collect::<Vec<_>>();
        rows.sort_by(|left, right| {
            left.seq_id
                .to_ascii_lowercase()
                .cmp(&right.seq_id.to_ascii_lowercase())
                .then(left.generated_at_unix_ms.cmp(&right.generated_at_unix_ms))
                .then(
                    left.track_id
                        .to_ascii_lowercase()
                        .cmp(&right.track_id.to_ascii_lowercase()),
                )
        });
        rows
    }

    pub fn get_flexibility_track(&self, track_id: &str) -> Result<FlexibilityTrack, EngineError> {
        let track_id = Self::normalize_analysis_id(track_id, "track")?;
        self.read_dotplot_analysis_store()
            .flexibility_tracks
            .get(track_id.as_str())
            .cloned()
            .ok_or_else(|| EngineError {
                code: ErrorCode::NotFound,
                message: format!("Flexibility track '{}' not found", track_id),
            })
    }

    fn read_guide_design_store_from_metadata(
        value: Option<&serde_json::Value>,
    ) -> GuideDesignStore {
        let mut store = value
            .cloned()
            .and_then(|v| serde_json::from_value::<GuideDesignStore>(v).ok())
            .unwrap_or_default();
        if store.schema.trim().is_empty() {
            store.schema = GUIDE_DESIGN_SCHEMA.to_string();
        }
        store
    }

    fn read_guide_design_store(&self) -> GuideDesignStore {
        Self::read_guide_design_store_from_metadata(
            self.state.metadata.get(GUIDE_DESIGN_METADATA_KEY),
        )
    }

    fn write_guide_design_store(&mut self, mut store: GuideDesignStore) -> Result<(), EngineError> {
        if store.guide_sets.is_empty()
            && store.practical_filter_reports.is_empty()
            && store.oligo_sets.is_empty()
            && store.latest_oligo_set_by_guide_set.is_empty()
            && store.audit_log.is_empty()
        {
            self.state.metadata.remove(GUIDE_DESIGN_METADATA_KEY);
            return Ok(());
        }
        store.schema = GUIDE_DESIGN_SCHEMA.to_string();
        store.updated_at_unix_ms = Self::now_unix_ms();
        let value = serde_json::to_value(store).map_err(|e| EngineError {
            code: ErrorCode::Internal,
            message: format!("Could not serialize guide-design metadata: {e}"),
        })?;
        self.state
            .metadata
            .insert(GUIDE_DESIGN_METADATA_KEY.to_string(), value);
        Ok(())
    }

    fn normalize_guide_set_id(raw: &str) -> Result<String, EngineError> {
        let trimmed = raw.trim();
        if trimmed.is_empty() {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: "guide_set_id cannot be empty".to_string(),
            });
        }
        Ok(trimmed.to_string())
    }

    fn normalize_guide_id(raw: &str, index: usize) -> String {
        let trimmed = raw.trim();
        if trimmed.is_empty() {
            format!("g_{:04}", index + 1)
        } else {
            trimmed.to_string()
        }
    }

    fn normalize_guide_strand(raw: &str) -> Result<String, EngineError> {
        match raw.trim().to_ascii_lowercase().as_str() {
            "+" | "plus" | "forward" | "fwd" => Ok("+".to_string()),
            "-" | "minus" | "reverse" | "rev" => Ok("-".to_string()),
            other => Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: format!("Unsupported guide strand '{}'; expected '+' or '-'", other),
            }),
        }
    }

    fn normalize_guide_candidate(
        guide: GuideCandidate,
        index: usize,
    ) -> Result<GuideCandidate, EngineError> {
        if guide.end_0based_exclusive <= guide.start_0based {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: format!(
                    "Guide {} has invalid genomic interval {}..{}",
                    index + 1,
                    guide.start_0based,
                    guide.end_0based_exclusive
                ),
            });
        }
        if guide.seq_id.trim().is_empty() {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: format!("Guide {} has empty seq_id", index + 1),
            });
        }
        if matches!(guide.rank, Some(0)) {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: format!("Guide {} rank must be >= 1", index + 1),
            });
        }
        let protospacer = Self::normalize_iupac_text(&guide.protospacer)?;
        if protospacer.is_empty() {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: format!("Guide {} has empty protospacer", index + 1),
            });
        }
        let pam = Self::normalize_iupac_text(&guide.pam)?;
        if pam.is_empty() {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: format!("Guide {} has empty PAM", index + 1),
            });
        }
        if guide.cut_offset_from_protospacer_start >= protospacer.len() {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: format!(
                    "Guide {} cut_offset_from_protospacer_start ({}) must be < protospacer length ({})",
                    index + 1,
                    guide.cut_offset_from_protospacer_start,
                    protospacer.len()
                ),
            });
        }
        Ok(GuideCandidate {
            guide_id: Self::normalize_guide_id(&guide.guide_id, index),
            seq_id: guide.seq_id.trim().to_string(),
            start_0based: guide.start_0based,
            end_0based_exclusive: guide.end_0based_exclusive,
            strand: Self::normalize_guide_strand(&guide.strand)?,
            protospacer,
            pam,
            nuclease: if guide.nuclease.trim().is_empty() {
                "SpCas9".to_string()
            } else {
                guide.nuclease.trim().to_string()
            },
            cut_offset_from_protospacer_start: guide.cut_offset_from_protospacer_start,
            rank: guide.rank,
        })
    }

    fn normalize_guide_candidates(
        guides: Vec<GuideCandidate>,
    ) -> Result<Vec<GuideCandidate>, EngineError> {
        if guides.is_empty() {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: "guide set requires at least one guide".to_string(),
            });
        }
        let mut normalized = Vec::with_capacity(guides.len());
        let mut seen_ids = HashSet::new();
        for (idx, guide) in guides.into_iter().enumerate() {
            let guide = Self::normalize_guide_candidate(guide, idx)?;
            if !seen_ids.insert(guide.guide_id.clone()) {
                return Err(EngineError {
                    code: ErrorCode::InvalidInput,
                    message: format!("Duplicate guide_id '{}' in guide set", guide.guide_id),
                });
            }
            normalized.push(guide);
        }
        normalized.sort_by(|a, b| {
            a.rank
                .unwrap_or(usize::MAX)
                .cmp(&b.rank.unwrap_or(usize::MAX))
                .then(a.guide_id.cmp(&b.guide_id))
        });
        Ok(normalized)
    }

    fn max_homopolymer_run_for_base(sequence: &[u8], base: u8) -> usize {
        let base = base.to_ascii_uppercase();
        let mut best = 0usize;
        let mut current = 0usize;
        for b in sequence {
            if b.to_ascii_uppercase() == base {
                current += 1;
                best = best.max(current);
            } else {
                current = 0;
            }
        }
        best
    }

    fn max_dinucleotide_repeat_units(sequence: &[u8]) -> usize {
        if sequence.len() < 2 {
            return 0;
        }
        let canonical = |b: u8| matches!(b.to_ascii_uppercase(), b'A' | b'C' | b'G' | b'T');
        let mut best = 1usize;
        for start in 0..(sequence.len() - 1) {
            let b0 = sequence[start].to_ascii_uppercase();
            let b1 = sequence[start + 1].to_ascii_uppercase();
            if !canonical(b0) || !canonical(b1) {
                continue;
            }
            let mut units = 1usize;
            let mut idx = start + 2;
            while idx + 1 < sequence.len()
                && sequence[idx].to_ascii_uppercase() == b0
                && sequence[idx + 1].to_ascii_uppercase() == b1
            {
                units += 1;
                idx += 2;
            }
            best = best.max(units);
        }
        best
    }

    fn normalize_practical_filter_config(
        mut config: GuidePracticalFilterConfig,
    ) -> Result<GuidePracticalFilterConfig, EngineError> {
        if let Some(min) = config.gc_min {
            if !(0.0..=1.0).contains(&min) {
                return Err(EngineError {
                    code: ErrorCode::InvalidInput,
                    message: format!("gc_min ({min}) must be between 0.0 and 1.0"),
                });
            }
        }
        if let Some(max) = config.gc_max {
            if !(0.0..=1.0).contains(&max) {
                return Err(EngineError {
                    code: ErrorCode::InvalidInput,
                    message: format!("gc_max ({max}) must be between 0.0 and 1.0"),
                });
            }
        }
        if let (Some(min), Some(max)) = (config.gc_min, config.gc_max) {
            if min > max {
                return Err(EngineError {
                    code: ErrorCode::InvalidInput,
                    message: format!("gc_min ({min}) must be <= gc_max ({max})"),
                });
            }
        }
        if let Some(max_run) = config.max_homopolymer_run {
            if max_run == 0 {
                return Err(EngineError {
                    code: ErrorCode::InvalidInput,
                    message: "max_homopolymer_run must be >= 1".to_string(),
                });
            }
        }
        if let Some(max_repeat) = config.max_dinucleotide_repeat_units {
            if max_repeat == 0 {
                return Err(EngineError {
                    code: ErrorCode::InvalidInput,
                    message: "max_dinucleotide_repeat_units must be >= 1".to_string(),
                });
            }
        }

        let mut normalized_per_base = HashMap::new();
        for (base, value) in &config.max_homopolymer_run_per_base {
            let key = base.trim().to_ascii_uppercase();
            if !matches!(key.as_str(), "A" | "C" | "G" | "T" | "U") {
                return Err(EngineError {
                    code: ErrorCode::InvalidInput,
                    message: format!(
                        "Unsupported max_homopolymer_run_per_base key '{}'; expected A/C/G/T",
                        base
                    ),
                });
            }
            if *value == 0 {
                return Err(EngineError {
                    code: ErrorCode::InvalidInput,
                    message: format!("max_homopolymer_run_per_base for '{}' must be >= 1", base),
                });
            }
            let key = if key == "U" {
                "T".to_string()
            } else {
                key.to_string()
            };
            normalized_per_base.insert(key, *value);
        }
        config.max_homopolymer_run_per_base = normalized_per_base;

        if let Some(required) = config.required_5prime_base.as_ref() {
            let normalized = Self::normalize_iupac_text(required)?;
            if normalized.len() != 1
                || !matches!(normalized.as_bytes()[0], b'A' | b'C' | b'G' | b'T' | b'U')
            {
                return Err(EngineError {
                    code: ErrorCode::InvalidInput,
                    message: "required_5prime_base must be one canonical nucleotide (A/C/G/T)"
                        .to_string(),
                });
            }
            config.required_5prime_base = Some(if normalized == "U" {
                "T".to_string()
            } else {
                normalized
            });
        }

        let mut motifs = vec![];
        for motif in &config.forbidden_motifs {
            let normalized = Self::normalize_iupac_text(motif)?;
            if !normalized.is_empty() {
                motifs.push(normalized);
            }
        }
        motifs.sort();
        motifs.dedup();
        config.forbidden_motifs = motifs;
        Ok(config)
    }

    pub fn list_guide_sets(&self) -> Vec<GuideSetSummary> {
        let store = self.read_guide_design_store();
        let mut names = store.guide_sets.keys().cloned().collect::<Vec<_>>();
        names.sort();
        names
            .into_iter()
            .filter_map(|name| store.guide_sets.get(&name).cloned())
            .map(|set| GuideSetSummary {
                guide_set_id: set.guide_set_id,
                guide_count: set.guides.len(),
                created_at_unix_ms: set.created_at_unix_ms,
                updated_at_unix_ms: set.updated_at_unix_ms,
            })
            .collect()
    }

    pub fn inspect_guide_set_page(
        &self,
        guide_set_id: &str,
        limit: usize,
        offset: usize,
    ) -> Result<(GuideSet, usize, usize), EngineError> {
        if limit == 0 {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: "Guide page limit must be >= 1".to_string(),
            });
        }
        let guide_set_id = Self::normalize_guide_set_id(guide_set_id)?;
        let store = self.read_guide_design_store();
        let mut set = store
            .guide_sets
            .get(&guide_set_id)
            .cloned()
            .ok_or_else(|| EngineError {
                code: ErrorCode::NotFound,
                message: format!("Guide set '{}' not found", guide_set_id),
            })?;
        let total = set.guides.len();
        let clamped_offset = offset.min(total);
        set.guides = set
            .guides
            .into_iter()
            .skip(clamped_offset)
            .take(limit)
            .collect();
        Ok((set, total, clamped_offset))
    }

    pub fn get_guide_practical_filter_report(
        &self,
        guide_set_id: &str,
    ) -> Result<GuidePracticalFilterReport, EngineError> {
        let guide_set_id = Self::normalize_guide_set_id(guide_set_id)?;
        let store = self.read_guide_design_store();
        store
            .practical_filter_reports
            .get(&guide_set_id)
            .cloned()
            .ok_or_else(|| EngineError {
                code: ErrorCode::NotFound,
                message: format!(
                    "No practical guide filter report found for '{}'",
                    guide_set_id
                ),
            })
    }

    pub fn list_guide_oligo_sets(&self, guide_set_id: Option<&str>) -> Vec<GuideOligoSet> {
        let store = self.read_guide_design_store();
        let filter = guide_set_id.map(|v| v.trim().to_string());
        let mut ids = store.oligo_sets.keys().cloned().collect::<Vec<_>>();
        ids.sort();
        ids.into_iter()
            .filter_map(|id| store.oligo_sets.get(&id).cloned())
            .filter(|set| {
                filter
                    .as_ref()
                    .map(|f| set.guide_set_id == *f)
                    .unwrap_or(true)
            })
            .collect()
    }

    pub fn get_guide_oligo_set(&self, oligo_set_id: &str) -> Result<GuideOligoSet, EngineError> {
        let id = oligo_set_id.trim();
        if id.is_empty() {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: "oligo_set_id cannot be empty".to_string(),
            });
        }
        let store = self.read_guide_design_store();
        store
            .oligo_sets
            .get(id)
            .cloned()
            .ok_or_else(|| EngineError {
                code: ErrorCode::NotFound,
                message: format!("Guide oligo set '{}' not found", id),
            })
    }

    fn built_in_guide_oligo_template(template_id: &str) -> Option<GuideOligoTemplateSpec> {
        match template_id.trim() {
            "lenti_bsmbi_u6_default" => Some(GuideOligoTemplateSpec {
                template_id: "lenti_bsmbi_u6_default".to_string(),
                description: "U6 sgRNA cloning oligos with BsmBI overhangs".to_string(),
                forward_prefix: "CACC".to_string(),
                forward_suffix: "".to_string(),
                reverse_prefix: "AAAC".to_string(),
                reverse_suffix: "C".to_string(),
                reverse_uses_reverse_complement_of_spacer: true,
                uppercase_output: true,
            }),
            "plain_forward_reverse" => Some(GuideOligoTemplateSpec {
                template_id: "plain_forward_reverse".to_string(),
                description: "Raw spacer and reverse-complement spacer".to_string(),
                forward_prefix: "".to_string(),
                forward_suffix: "".to_string(),
                reverse_prefix: "".to_string(),
                reverse_suffix: "".to_string(),
                reverse_uses_reverse_complement_of_spacer: true,
                uppercase_output: true,
            }),
            _ => None,
        }
    }

    fn normalize_oligo_set_id(raw: &str) -> Result<String, EngineError> {
        let trimmed = raw.trim();
        if trimmed.is_empty() {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: "oligo_set_id cannot be empty".to_string(),
            });
        }
        Ok(trimmed.to_string())
    }

    fn unique_oligo_set_id(store: &GuideDesignStore, base: &str) -> String {
        if !store.oligo_sets.contains_key(base) {
            return base.to_string();
        }
        let mut idx = 2usize;
        loop {
            let candidate = format!("{base}_{idx}");
            if !store.oligo_sets.contains_key(&candidate) {
                return candidate;
            }
            idx += 1;
        }
    }

    fn resolve_oligo_set_for_export(
        store: &GuideDesignStore,
        guide_set_id: &str,
        requested_oligo_set_id: Option<&str>,
    ) -> Result<GuideOligoSet, EngineError> {
        if let Some(id) = requested_oligo_set_id {
            let id = id.trim();
            if id.is_empty() {
                return Err(EngineError {
                    code: ErrorCode::InvalidInput,
                    message: "oligo_set_id cannot be empty".to_string(),
                });
            }
            let set = store.oligo_sets.get(id).ok_or_else(|| EngineError {
                code: ErrorCode::NotFound,
                message: format!("Guide oligo set '{}' not found", id),
            })?;
            if set.guide_set_id != guide_set_id {
                return Err(EngineError {
                    code: ErrorCode::InvalidInput,
                    message: format!(
                        "Guide oligo set '{}' belongs to '{}' (expected '{}')",
                        id, set.guide_set_id, guide_set_id
                    ),
                });
            }
            return Ok(set.clone());
        }

        let latest = store
            .latest_oligo_set_by_guide_set
            .get(guide_set_id)
            .ok_or_else(|| EngineError {
                code: ErrorCode::NotFound,
                message: format!(
                    "No oligo set is registered for guide set '{}'; run GenerateGuideOligos first",
                    guide_set_id
                ),
            })?;
        store
            .oligo_sets
            .get(latest)
            .cloned()
            .ok_or_else(|| EngineError {
                code: ErrorCode::NotFound,
                message: format!(
                    "Guide oligo set '{}' referenced by '{}' is missing",
                    latest, guide_set_id
                ),
            })
    }

    fn normalize_workflow_macro_template_name(raw: &str) -> Result<String, EngineError> {
        let trimmed = raw.trim();
        if trimmed.is_empty() {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: "Workflow macro template name cannot be empty".to_string(),
            });
        }
        Ok(trimmed.to_string())
    }

    fn normalize_workflow_macro_param_name(raw: &str) -> Result<String, EngineError> {
        let trimmed = raw.trim();
        if trimmed.is_empty() {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: "Workflow macro parameter name cannot be empty".to_string(),
            });
        }
        let valid = trimmed.chars().enumerate().all(|(idx, ch)| match idx {
            0 => ch.is_ascii_alphabetic() || ch == '_',
            _ => ch.is_ascii_alphanumeric() || ch == '_',
        });
        if !valid {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: format!(
                    "Invalid workflow macro parameter name '{}' (expected [A-Za-z_][A-Za-z0-9_]*)",
                    trimmed
                ),
            });
        }
        Ok(trimmed.to_string())
    }

    fn normalize_workflow_macro_template_details_url(
        raw: Option<String>,
    ) -> Result<Option<String>, EngineError> {
        let Some(raw) = raw else {
            return Ok(None);
        };
        let trimmed = raw.trim();
        if trimmed.is_empty() {
            return Ok(None);
        }
        let lower = trimmed.to_ascii_lowercase();
        if !(lower.starts_with("https://") || lower.starts_with("http://")) {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: format!(
                    "Workflow macro template details_url '{}' must start with http:// or https://",
                    trimmed
                ),
            });
        }
        Ok(Some(trimmed.to_string()))
    }

    fn read_workflow_macro_template_store_from_metadata(
        value: Option<&serde_json::Value>,
    ) -> WorkflowMacroTemplateStore {
        let mut store = value
            .cloned()
            .and_then(|v| serde_json::from_value::<WorkflowMacroTemplateStore>(v).ok())
            .unwrap_or_default();
        if store.schema.trim().is_empty() {
            store.schema = WORKFLOW_MACRO_TEMPLATES_SCHEMA.to_string();
        }
        store
    }

    fn read_workflow_macro_template_store(&self) -> WorkflowMacroTemplateStore {
        Self::read_workflow_macro_template_store_from_metadata(
            self.state
                .metadata
                .get(WORKFLOW_MACRO_TEMPLATES_METADATA_KEY),
        )
    }

    fn write_workflow_macro_template_store(
        &mut self,
        mut store: WorkflowMacroTemplateStore,
    ) -> Result<(), EngineError> {
        if store.templates.is_empty() {
            self.state
                .metadata
                .remove(WORKFLOW_MACRO_TEMPLATES_METADATA_KEY);
            return Ok(());
        }
        store.schema = WORKFLOW_MACRO_TEMPLATES_SCHEMA.to_string();
        store.updated_at_unix_ms = Self::now_unix_ms();
        let value = serde_json::to_value(store).map_err(|e| EngineError {
            code: ErrorCode::Internal,
            message: format!("Could not serialize workflow macro template metadata: {e}"),
        })?;
        self.state
            .metadata
            .insert(WORKFLOW_MACRO_TEMPLATES_METADATA_KEY.to_string(), value);
        Ok(())
    }

    pub fn list_workflow_macro_templates(&self) -> Vec<WorkflowMacroTemplateSummary> {
        let store = self.read_workflow_macro_template_store();
        let mut names = store.templates.keys().cloned().collect::<Vec<_>>();
        names.sort();
        names
            .into_iter()
            .filter_map(|name| store.templates.get(&name))
            .map(|template| WorkflowMacroTemplateSummary {
                name: template.name.clone(),
                description: template.description.clone(),
                details_url: template.details_url.clone(),
                parameter_count: template.parameters.len(),
                created_at_unix_ms: template.created_at_unix_ms,
                updated_at_unix_ms: template.updated_at_unix_ms,
            })
            .collect()
    }

    pub fn get_workflow_macro_template(
        &self,
        name: &str,
    ) -> Result<WorkflowMacroTemplate, EngineError> {
        let name = Self::normalize_workflow_macro_template_name(name)?;
        let store = self.read_workflow_macro_template_store();
        store
            .templates
            .get(&name)
            .cloned()
            .ok_or_else(|| EngineError {
                code: ErrorCode::NotFound,
                message: format!("Workflow macro template '{}' not found", name),
            })
    }

    pub fn render_workflow_macro_template_script(
        &self,
        name: &str,
        bindings: &HashMap<String, String>,
    ) -> Result<String, EngineError> {
        let template = self.get_workflow_macro_template(name)?;
        let declared = template
            .parameters
            .iter()
            .map(|p| p.name.clone())
            .collect::<HashSet<_>>();
        for key in bindings.keys() {
            if !declared.contains(key) {
                return Err(EngineError {
                    code: ErrorCode::InvalidInput,
                    message: format!(
                        "Workflow macro template '{}' does not define parameter '{}'",
                        template.name, key
                    ),
                });
            }
        }

        let mut resolved: HashMap<String, String> = HashMap::new();
        for param in &template.parameters {
            if let Some(value) = bindings.get(&param.name) {
                resolved.insert(param.name.clone(), value.clone());
            } else if let Some(default_value) = &param.default_value {
                resolved.insert(param.name.clone(), default_value.clone());
            } else if param.required {
                return Err(EngineError {
                    code: ErrorCode::InvalidInput,
                    message: format!(
                        "Workflow macro template '{}' is missing required parameter '{}'",
                        template.name, param.name
                    ),
                });
            }
        }

        let placeholder_regex =
            Regex::new(r"\$\{([A-Za-z_][A-Za-z0-9_]*)\}").map_err(|e| EngineError {
                code: ErrorCode::Internal,
                message: format!("Could not compile workflow macro placeholder regex: {e}"),
            })?;
        let mut missing: Vec<String> = vec![];
        for captures in placeholder_regex.captures_iter(&template.script) {
            if let Some(name) = captures.get(1).map(|m| m.as_str()) {
                if !declared.contains(name) {
                    return Err(EngineError {
                        code: ErrorCode::InvalidInput,
                        message: format!(
                            "Workflow macro template '{}' references undeclared parameter '{}'",
                            template.name, name
                        ),
                    });
                }
                if !resolved.contains_key(name) {
                    missing.push(name.to_string());
                }
            }
        }
        if !missing.is_empty() {
            missing.sort();
            missing.dedup();
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: format!(
                    "Workflow macro template '{}' is missing parameter bindings for: {}",
                    template.name,
                    missing.join(", ")
                ),
            });
        }

        let rendered = placeholder_regex
            .replace_all(&template.script, |captures: &regex::Captures<'_>| {
                let key = captures.get(1).map(|m| m.as_str()).unwrap_or_default();
                resolved.get(key).cloned().unwrap_or_default()
            })
            .to_string();
        Ok(rendered)
    }

    fn normalize_candidate_macro_template_name(raw: &str) -> Result<String, EngineError> {
        let trimmed = raw.trim();
        if trimmed.is_empty() {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: "Candidate macro template name cannot be empty".to_string(),
            });
        }
        Ok(trimmed.to_string())
    }

    fn normalize_candidate_macro_param_name(raw: &str) -> Result<String, EngineError> {
        let trimmed = raw.trim();
        if trimmed.is_empty() {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: "Candidate macro parameter name cannot be empty".to_string(),
            });
        }
        let valid = trimmed.chars().enumerate().all(|(idx, ch)| match idx {
            0 => ch.is_ascii_alphabetic() || ch == '_',
            _ => ch.is_ascii_alphanumeric() || ch == '_',
        });
        if !valid {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: format!(
                    "Invalid candidate macro parameter name '{}' (expected [A-Za-z_][A-Za-z0-9_]*)",
                    trimmed
                ),
            });
        }
        Ok(trimmed.to_string())
    }

    fn normalize_candidate_macro_template_details_url(
        raw: Option<String>,
    ) -> Result<Option<String>, EngineError> {
        let Some(raw) = raw else {
            return Ok(None);
        };
        let trimmed = raw.trim();
        if trimmed.is_empty() {
            return Ok(None);
        }
        let lower = trimmed.to_ascii_lowercase();
        if !(lower.starts_with("https://") || lower.starts_with("http://")) {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: format!(
                    "Candidate macro template details_url '{}' must start with http:// or https://",
                    trimmed
                ),
            });
        }
        Ok(Some(trimmed.to_string()))
    }

    fn read_candidate_macro_template_store_from_metadata(
        value: Option<&serde_json::Value>,
    ) -> CandidateMacroTemplateStore {
        let mut store = value
            .cloned()
            .and_then(|v| serde_json::from_value::<CandidateMacroTemplateStore>(v).ok())
            .unwrap_or_default();
        if store.schema.trim().is_empty() {
            store.schema = CANDIDATE_MACRO_TEMPLATES_SCHEMA.to_string();
        }
        store
    }

    fn read_candidate_macro_template_store(&self) -> CandidateMacroTemplateStore {
        Self::read_candidate_macro_template_store_from_metadata(
            self.state
                .metadata
                .get(CANDIDATE_MACRO_TEMPLATES_METADATA_KEY),
        )
    }

    fn write_candidate_macro_template_store(
        &mut self,
        mut store: CandidateMacroTemplateStore,
    ) -> Result<(), EngineError> {
        if store.templates.is_empty() {
            self.state
                .metadata
                .remove(CANDIDATE_MACRO_TEMPLATES_METADATA_KEY);
            return Ok(());
        }
        store.schema = CANDIDATE_MACRO_TEMPLATES_SCHEMA.to_string();
        store.updated_at_unix_ms = Self::now_unix_ms();
        let value = serde_json::to_value(store).map_err(|e| EngineError {
            code: ErrorCode::Internal,
            message: format!("Could not serialize candidate macro template metadata: {e}"),
        })?;
        self.state
            .metadata
            .insert(CANDIDATE_MACRO_TEMPLATES_METADATA_KEY.to_string(), value);
        Ok(())
    }

    pub fn list_candidate_macro_templates(&self) -> Vec<CandidateMacroTemplateSummary> {
        let store = self.read_candidate_macro_template_store();
        let mut names = store.templates.keys().cloned().collect::<Vec<_>>();
        names.sort();
        names
            .into_iter()
            .filter_map(|name| store.templates.get(&name))
            .map(|template| CandidateMacroTemplateSummary {
                name: template.name.clone(),
                description: template.description.clone(),
                details_url: template.details_url.clone(),
                parameter_count: template.parameters.len(),
                created_at_unix_ms: template.created_at_unix_ms,
                updated_at_unix_ms: template.updated_at_unix_ms,
            })
            .collect()
    }

    pub fn get_candidate_macro_template(
        &self,
        name: &str,
    ) -> Result<CandidateMacroTemplate, EngineError> {
        let name = Self::normalize_candidate_macro_template_name(name)?;
        let store = self.read_candidate_macro_template_store();
        store
            .templates
            .get(&name)
            .cloned()
            .ok_or_else(|| EngineError {
                code: ErrorCode::NotFound,
                message: format!("Candidate macro template '{}' not found", name),
            })
    }

    pub fn render_candidate_macro_template_script(
        &self,
        name: &str,
        bindings: &HashMap<String, String>,
    ) -> Result<String, EngineError> {
        let template = self.get_candidate_macro_template(name)?;
        let declared = template
            .parameters
            .iter()
            .map(|p| p.name.clone())
            .collect::<HashSet<_>>();
        for key in bindings.keys() {
            if !declared.contains(key) {
                return Err(EngineError {
                    code: ErrorCode::InvalidInput,
                    message: format!(
                        "Candidate macro template '{}' does not define parameter '{}'",
                        template.name, key
                    ),
                });
            }
        }

        let mut resolved: HashMap<String, String> = HashMap::new();
        for param in &template.parameters {
            if let Some(value) = bindings.get(&param.name) {
                resolved.insert(param.name.clone(), value.clone());
            } else if let Some(default_value) = &param.default_value {
                resolved.insert(param.name.clone(), default_value.clone());
            } else if param.required {
                return Err(EngineError {
                    code: ErrorCode::InvalidInput,
                    message: format!(
                        "Candidate macro template '{}' is missing required parameter '{}'",
                        template.name, param.name
                    ),
                });
            }
        }

        let placeholder_regex =
            Regex::new(r"\$\{([A-Za-z_][A-Za-z0-9_]*)\}").map_err(|e| EngineError {
                code: ErrorCode::Internal,
                message: format!("Could not compile candidate macro placeholder regex: {e}"),
            })?;
        let mut missing: Vec<String> = vec![];
        for captures in placeholder_regex.captures_iter(&template.script) {
            if let Some(name) = captures.get(1).map(|m| m.as_str()) {
                if !declared.contains(name) {
                    return Err(EngineError {
                        code: ErrorCode::InvalidInput,
                        message: format!(
                            "Candidate macro template '{}' references undeclared parameter '{}'",
                            template.name, name
                        ),
                    });
                }
                if !resolved.contains_key(name) {
                    missing.push(name.to_string());
                }
            }
        }
        if !missing.is_empty() {
            missing.sort();
            missing.dedup();
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: format!(
                    "Candidate macro template '{}' is missing parameter bindings for: {}",
                    template.name,
                    missing.join(", ")
                ),
            });
        }

        let rendered = placeholder_regex
            .replace_all(&template.script, |captures: &regex::Captures<'_>| {
                let key = captures.get(1).map(|m| m.as_str()).unwrap_or_default();
                resolved.get(key).cloned().unwrap_or_default()
            })
            .to_string();
        Ok(rendered)
    }

    fn latest_genome_anchor_for_seq(
        &self,
        seq_id: &str,
    ) -> Result<GenomeSequenceAnchor, EngineError> {
        let Some(provenance) = self.state.metadata.get(PROVENANCE_METADATA_KEY) else {
            return Err(EngineError {
                code: ErrorCode::NotFound,
                message: format!(
                    "Sequence '{seq_id}' has no genome anchor provenance (missing metadata '{}')",
                    PROVENANCE_METADATA_KEY
                ),
            });
        };
        let Some(entries) = provenance
            .get(GENOME_EXTRACTIONS_METADATA_KEY)
            .and_then(|v| v.as_array())
        else {
            return Err(EngineError {
                code: ErrorCode::NotFound,
                message: format!(
                    "Sequence '{seq_id}' has no genome anchor provenance (missing '{}')",
                    GENOME_EXTRACTIONS_METADATA_KEY
                ),
            });
        };

        let mut latest: Option<GenomeSequenceAnchor> = None;
        let mut latest_recorded_at = 0u128;

        for entry in entries {
            let Some(entry_seq_id) = entry.get("seq_id").and_then(|v| v.as_str()) else {
                continue;
            };
            if entry_seq_id != seq_id {
                continue;
            }
            let Some(chromosome) = entry.get("chromosome").and_then(|v| v.as_str()) else {
                continue;
            };
            let Some(start_1based) = entry.get("start_1based").and_then(|v| v.as_u64()) else {
                continue;
            };
            let Some(end_1based) = entry.get("end_1based").and_then(|v| v.as_u64()) else {
                continue;
            };
            if start_1based == 0 || end_1based < start_1based {
                continue;
            }
            let anchor_strand = entry
                .get("anchor_strand")
                .and_then(|v| v.as_str())
                .and_then(|v| v.trim().chars().next())
                .filter(|c| matches!(c, '+' | '-'));
            let anchor_verified = entry.get("anchor_verified").and_then(|v| v.as_bool());
            let catalog_path = entry
                .get("catalog_path")
                .and_then(|v| v.as_str())
                .map(str::trim)
                .filter(|v| !v.is_empty())
                .map(|v| v.to_string());
            let cache_dir = entry
                .get("cache_dir")
                .and_then(|v| v.as_str())
                .map(str::trim)
                .filter(|v| !v.is_empty())
                .map(|v| v.to_string());
            let recorded_at = entry
                .get("recorded_at_unix_ms")
                .and_then(|v| v.as_u64())
                .map(|v| v as u128)
                .unwrap_or(0);
            if latest.is_none() || recorded_at >= latest_recorded_at {
                latest_recorded_at = recorded_at;
                latest = Some(GenomeSequenceAnchor {
                    genome_id: entry
                        .get("genome_id")
                        .and_then(|v| v.as_str())
                        .unwrap_or_default()
                        .to_string(),
                    chromosome: chromosome.to_string(),
                    start_1based: start_1based as usize,
                    end_1based: end_1based as usize,
                    strand: anchor_strand,
                    anchor_verified,
                    catalog_path,
                    cache_dir,
                });
            }
        }

        latest.ok_or_else(|| EngineError {
            code: ErrorCode::NotFound,
            message: format!(
                "Sequence '{seq_id}' is not anchored to a genome interval; run ExtractGenomeRegion/ExtractGenomeGene/ExtendGenomeAnchor first"
            ),
        })
    }

    fn normalize_chromosome_alias(raw: &str) -> String {
        let trimmed = raw.trim();
        if trimmed.is_empty() {
            return String::new();
        }
        let lower = trimmed.to_ascii_lowercase();
        let core = lower.strip_prefix("chr").unwrap_or(&lower);
        if core.is_empty() {
            return String::new();
        }
        if matches!(core, "m" | "mt" | "mitochondria" | "mitochondrion") {
            return "mt".to_string();
        }
        if core.chars().all(|ch| ch.is_ascii_digit()) {
            let normalized = core.trim_start_matches('0');
            if normalized.is_empty() {
                return "0".to_string();
            }
            return normalized.to_string();
        }
        if core.len() == 1 && matches!(core, "x" | "y" | "w" | "z") {
            return core.to_ascii_uppercase();
        }
        core.to_string()
    }

    fn chromosomes_match(left: &str, right: &str) -> bool {
        Self::normalize_chromosome_alias(left) == Self::normalize_chromosome_alias(right)
    }

    fn append_chromosome_mismatch_warning(
        report: &mut GenomeBedTrackImportReport,
        anchor_chromosome: &str,
        source_label: &str,
        mismatch_counts: &HashMap<String, usize>,
    ) {
        if mismatch_counts.is_empty() {
            return;
        }
        let mut sorted = mismatch_counts
            .iter()
            .map(|(chrom, count)| (chrom.clone(), *count))
            .collect::<Vec<_>>();
        sorted.sort_by(|a, b| b.1.cmp(&a.1).then(a.0.cmp(&b.0)));
        report.skipped_wrong_chromosome_examples = sorted
            .iter()
            .take(6)
            .map(|(chrom, count)| format!("{chrom} ({count})"))
            .collect();
        let seen = sorted
            .iter()
            .take(3)
            .map(|(chrom, count)| format!("{chrom} ({count})"))
            .collect::<Vec<_>>()
            .join(", ");
        report.warnings.push(format!(
            "{} record(s) in {} input did not match anchor chromosome '{}' (examples: {})",
            report.skipped_wrong_chromosome, source_label, anchor_chromosome, seen
        ));
    }
}

impl Engine for GentleEngine {
    fn apply(&mut self, op: Operation) -> Result<OpResult, EngineError> {
        let run_id = "interactive".to_string();
        let mut noop = |_p: OperationProgress| true;
        let checkpoint = self.maybe_capture_checkpoint(&op);
        let result = self.apply_internal(op.clone(), &run_id, &mut noop)?;
        self.journal.push(OperationRecord {
            run_id,
            op,
            result: result.clone(),
        });
        if let Some(checkpoint) = checkpoint {
            self.push_undo_checkpoint(checkpoint);
        }
        Ok(result)
    }

    fn apply_workflow(&mut self, wf: Workflow) -> Result<Vec<OpResult>, EngineError> {
        let mut results = Vec::new();
        for op in &wf.ops {
            let mut noop = |_p: OperationProgress| true;
            let checkpoint = self.maybe_capture_checkpoint(op);
            let result = self.apply_internal(op.clone(), &wf.run_id, &mut noop)?;
            self.journal.push(OperationRecord {
                run_id: wf.run_id.clone(),
                op: op.clone(),
                result: result.clone(),
            });
            if let Some(checkpoint) = checkpoint {
                self.push_undo_checkpoint(checkpoint);
            }
            results.push(result);
        }
        Ok(results)
    }

    fn snapshot(&self) -> &ProjectState {
        &self.state
    }
}

#[cfg(test)]
mod tests;
