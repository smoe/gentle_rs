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
    amino_acids::{STOP_CODON, UNKNOWN_CODON},
    app::GENtleApp,
    dna_sequence::DNAsequence,
    ensembl_protein::EnsemblProteinEntry,
    enzymes::{active_restriction_enzymes, default_preferred_restriction_enzyme_names},
    feature_location::{collect_location_ranges_usize, feature_is_reverse},
    genomes::{
        BlastExternalBinaryPreflightReport, DEFAULT_HELPER_CATALOG_DISCOVERY_TOKEN,
        DEFAULT_REFERENCE_CATALOG_DISCOVERY_TOKEN, EnsemblCatalogUpdatePreview,
        EnsemblCatalogUpdateReport, EnsemblInstallableGenomeCatalog,
        EnsemblQuickInstallCatalogWriteReport, EnsemblQuickInstallPreview,
        EnsemblQuickInstallReport, GenomeBlastReport, GenomeCatalog,
        GenomeCatalogEntryRemovalReport, GenomeCatalogListEntry, GenomeGeneRecord,
        GenomeSourcePlan, GenomeTranscriptRecord, HelperConstructInterpretation,
        PrepareGenomeActivityStatus, PrepareGenomePlan, PrepareGenomeProgress, PrepareGenomeReport,
        PreparedCacheCleanupReport, PreparedCacheCleanupRequest, PreparedCacheInspectionReport,
        PreparedGenomeCompatibilityInspection, PreparedGenomeFallbackPolicy,
        PreparedGenomeInspection, PreparedGenomeRemovalReport,
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
        UniprotAaGenomicSegment, UniprotEntry, UniprotEntrySummary, UniprotGenomeProjection,
        UniprotGenomeProjectionSummary, UniprotTranscriptProjection,
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

pub use gentle_protocol::{
    Arrangement, ArrangementMode, Container, ContainerId, ContainerKind, ContainerState,
    GelBufferModel, GelRunConditions, GelTopologyForm, LineageEdge, LineageGraph,
    LineageMacroInstance, LineageMacroPortBinding, LineageNode, MacroInstanceStatus, NodeId, OpId,
    ProteinExternalOpinionSource, ProteinFeatureFilter, Rack, RackAuthoringTemplate,
    RackCarrierLabelPreset, RackFillDirection, RackLabelSheetPreset, RackOccupant,
    RackPhysicalTemplateFamily, RackPhysicalTemplateKind, RackPhysicalTemplateSpec,
    RackPlacementEntry, RackProfileKind, RackProfileSnapshot, RunId, SeqId, SequenceOrigin,
};

pub const DEFAULT_HOST_PROFILE_CATALOG_PATH: &str = "assets/host_profiles.json";

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum PrepareReferenceGenomeMode {
    PrepareOrReuse,
    ReindexCachedFiles,
    RefreshFromSources,
}

#[derive(Debug, Clone, Copy)]
struct ConstructSelectionRule {
    rule_id: &'static str,
    label: &'static str,
    category: &'static str,
    feature_match_tokens: &'static [&'static str],
    condition_signal_tokens: &'static [&'static str],
    helper_offered_function_tokens: &'static [&'static str],
    helper_component_kind_tokens: &'static [&'static str],
    helper_component_tag_tokens: &'static [&'static str],
    helper_component_attribute_terms: &'static [&'static str],
}

#[derive(Debug, Clone, Serialize, Default)]
struct ConstructConditionSignal {
    category: String,
    token: String,
    label: String,
    source_condition: String,
    matched_terms: Vec<String>,
    temperature_celsius: Option<f64>,
    notes: Vec<String>,
}

#[derive(Debug, Clone, Serialize, Default, PartialEq, Eq)]
pub(crate) struct RoutinePreferenceContext {
    pub helper_profile_id: Option<String>,
    pub construct_reasoning_seq_id: Option<String>,
    pub helper_resolution_status: String,
    pub explicit_preferred_routine_families: Vec<String>,
    pub helper_derived_preferred_routine_families: Vec<String>,
    pub construct_strategy_derived_preferred_routine_families: Vec<String>,
    pub variant_derived_preferred_routine_families: Vec<String>,
    pub effective_preferred_routine_families: Vec<String>,
    pub helper_offered_functions: Vec<String>,
    pub helper_component_labels: Vec<String>,
    pub variant_effect_tags: Vec<String>,
    pub variant_suggested_assay_ids: Vec<String>,
    pub rationale: Vec<String>,
}

#[derive(Debug, Clone, Serialize, Default)]
struct ConstructSelectionRuleMatch {
    rule_id: String,
    label: String,
    category: String,
    status: String,
    matched_medium_conditions: Vec<String>,
    matched_medium_terms: Vec<String>,
    matched_condition_signal_tokens: Vec<String>,
    candidate_labels: Vec<String>,
    candidate_evidence_ids: Vec<String>,
    candidate_match_tokens: Vec<String>,
    helper_profile_id: Option<String>,
    helper_offered_function_matches: Vec<String>,
    helper_component_ids: Vec<String>,
    helper_component_labels: Vec<String>,
    helper_component_kind_matches: Vec<String>,
    helper_component_tag_matches: Vec<String>,
    helper_component_attribute_matches: Vec<String>,
}

#[derive(Debug, Clone, Serialize, Default)]
struct ConstructSequenceRestrictionMethylationPatterns {
    dam_site_count: usize,
    dcm_site_count: usize,
    eco_k_target_site_count: usize,
    present_pattern_tokens: Vec<String>,
}

#[derive(Debug, Clone, Serialize, Default)]
struct ConstructHostRestrictionMethylationStep {
    step_id: String,
    host_profile_id: String,
    role: String,
    trait_tokens: Vec<String>,
    notes: Vec<String>,
}

#[derive(Debug, Clone, Serialize, Default)]
struct ConstructRestrictionMethylationConflict {
    conflict_id: String,
    label: String,
    upstream_step_id: String,
    downstream_step_id: String,
    rationale: String,
}

#[derive(Debug, Clone, Serialize, Default)]
struct ConstructAdapterRestrictionCaptureMotifSummary {
    enzyme_name: String,
    motif_role: String,
    resolution_status: String,
    site_geometry: Option<String>,
    internal_site_count: usize,
    internal_site_status: String,
    internal_site_ranges_0based: Vec<[usize; 2]>,
}

#[derive(Debug, Clone, Serialize, Default)]
struct ConstructAdapterRestrictionCapturePlanSummary {
    capture_id: String,
    restriction_enzyme_name: String,
    enzyme_resolution_status: String,
    adapter_style: String,
    blunt_insert_required: bool,
    protection_mode: String,
    extra_retrieval_enzyme_names: Vec<String>,
    capture_site_geometry: Option<String>,
    internal_site_count: usize,
    internal_site_status: String,
    internal_site_ranges_0based: Vec<[usize; 2]>,
    named_motif_count: usize,
    resolved_motif_count: usize,
    all_named_motifs_present_on_insert: bool,
    motif_summaries: Vec<ConstructAdapterRestrictionCaptureMotifSummary>,
    notes: Vec<String>,
}

#[derive(Debug, Clone, Serialize, Default)]
struct ConstructAdapterRestrictionCaptureReviewItem {
    capture_id: String,
    issue_id: String,
    label: String,
    rationale: String,
}

#[derive(Debug, Clone, Serialize, Default)]
struct ConstructVariantCodingPrediction {
    cds_label: String,
    transcript_label: Option<String>,
    effect_tag: String,
    codon_index_1based: usize,
    ref_codon: String,
    alt_codon: String,
    ref_amino_acid: String,
    alt_amino_acid: String,
    translation_table: usize,
}

#[derive(Debug, Clone, Serialize, Default)]
struct ConstructVariantReasoningSummary {
    evidence_id: String,
    label: String,
    start_0based: usize,
    end_0based_exclusive: usize,
    strand: Option<String>,
    variant_class: Option<String>,
    genomic_ref: Option<String>,
    genomic_alt: Option<String>,
    overlapping_roles: Vec<String>,
    overlapping_labels: Vec<String>,
    gene_labels: Vec<String>,
    transcript_labels: Vec<String>,
    effect_tags: Vec<String>,
    suggested_assay_ids: Vec<String>,
    transcript_context_status: String,
    transcript_effect_summaries: Vec<ConstructTranscriptVariantEffectSummary>,
    coding_predictions: Vec<ConstructVariantCodingPrediction>,
}

#[derive(Debug, Clone, Serialize, Default)]
struct ConstructTranscriptVariantEffectSummary {
    transcript_label: String,
    effect_tags: Vec<String>,
    suggested_assay_ids: Vec<String>,
    source: String,
}

pub use crate::feature_expert::{
    FeatureExpertTarget, FeatureExpertView, ISOFORM_ARCHITECTURE_EXPERT_INSTRUCTION,
    IsoformArchitectureCdsAaSegment, IsoformArchitectureExpertView,
    IsoformArchitectureProteinDomain, IsoformArchitectureProteinLane,
    IsoformArchitectureTranscriptLane, RESTRICTION_EXPERT_INSTRUCTION, RestrictionSiteExpertView,
    SPLICING_EXPERT_INSTRUCTION, SplicingBoundaryMarker, SplicingEventSummary,
    SplicingExonCdsPhase, SplicingExonSummary, SplicingExpertView, SplicingJunctionArc,
    SplicingMatrixRow, SplicingRange, SplicingScopePreset, SplicingTranscriptLane,
    TFBS_EXPERT_INSTRUCTION, TfbsExpertColumn, TfbsExpertView, TranscriptProteinComparison,
    TranscriptProteinComparisonStatus, TranscriptProteinDerivation,
    TranscriptProteinExternalOpinion,
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
const RESTRICTION_CLONING_PCR_HANDOFF_REPORT_SCHEMA: &str =
    "gentle.restriction_cloning_pcr_handoff.v1";
pub const PROTEIN_DERIVATION_REPORTS_METADATA_KEY: &str = "protein_derivation_reports";
const PROTEIN_DERIVATION_REPORTS_SCHEMA: &str = "gentle.protein_derivation_reports.v1";
pub const PROTEIN_DERIVATION_REPORT_SCHEMA: &str = "gentle.protein_derivation_report.v1";
pub const REVERSE_TRANSLATION_REPORTS_METADATA_KEY: &str = "reverse_translation_reports";
const REVERSE_TRANSLATION_REPORTS_SCHEMA: &str = "gentle.reverse_translation_reports.v1";
pub const REVERSE_TRANSLATION_REPORT_SCHEMA: &str = "gentle.reverse_translation_report.v1";
pub const SEQUENCING_TRACES_METADATA_KEY: &str = "sequencing_traces";
const SEQUENCING_TRACES_SCHEMA: &str = "gentle.sequencing_traces.v1";
pub const SEQUENCING_TRACE_RECORD_SCHEMA: &str = "gentle.sequencing_trace_record.v2";
pub const SEQUENCING_TRACE_IMPORT_REPORT_SCHEMA: &str = "gentle.sequencing_trace_import_report.v2";
pub const SEQUENCING_CONFIRMATION_REPORTS_METADATA_KEY: &str = "sequencing_confirmation_reports";
const SEQUENCING_CONFIRMATION_REPORTS_SCHEMA: &str = "gentle.sequencing_confirmation_reports.v1";
pub const SEQUENCING_CONFIRMATION_REPORT_SCHEMA: &str = "gentle.sequencing_confirmation_report.v1";
pub const SEQUENCING_PRIMER_OVERLAY_REPORT_SCHEMA: &str =
    "gentle.sequencing_primer_overlay_report.v1";
pub const SEQUENCING_CONFIRMATION_SUPPORT_TSV_SCHEMA: &str =
    "gentle.sequencing_confirmation_support_tsv.v1";
pub const PLANNING_METADATA_KEY: &str = "planning";
pub const PLANNING_PROFILE_SCHEMA: &str = "gentle.planning_profile.v1";
pub const PLANNING_OBJECTIVE_SCHEMA: &str = "gentle.planning_objective.v1";
pub const PLANNING_ESTIMATE_SCHEMA: &str = "gentle.planning_estimate.v1";
pub const PLANNING_SUGGESTION_SCHEMA: &str = "gentle.planning_suggestion.v1";
pub const PLANNING_SYNC_STATUS_SCHEMA: &str = "gentle.planning_sync_status.v1";
const PLANNING_STORE_SCHEMA: &str = "gentle.planning_store.v1";
pub const CONSTRUCT_REASONING_METADATA_KEY: &str = "construct_reasoning";
pub const CONSTRUCT_OBJECTIVE_SCHEMA: &str = gentle_protocol::CONSTRUCT_OBJECTIVE_SCHEMA;
pub const DESIGN_EVIDENCE_SCHEMA: &str = gentle_protocol::DESIGN_EVIDENCE_SCHEMA;
pub const DESIGN_FACT_SCHEMA: &str = gentle_protocol::DESIGN_FACT_SCHEMA;
pub const DESIGN_DECISION_NODE_SCHEMA: &str = gentle_protocol::DESIGN_DECISION_NODE_SCHEMA;
pub const CONSTRUCT_CANDIDATE_SCHEMA: &str = gentle_protocol::CONSTRUCT_CANDIDATE_SCHEMA;
pub const ANNOTATION_CANDIDATE_SCHEMA: &str = gentle_protocol::ANNOTATION_CANDIDATE_SCHEMA;
pub const ANNOTATION_CANDIDATE_SUMMARY_SCHEMA: &str =
    gentle_protocol::ANNOTATION_CANDIDATE_SUMMARY_SCHEMA;
pub const CONSTRUCT_REASONING_GRAPH_SCHEMA: &str =
    gentle_protocol::CONSTRUCT_REASONING_GRAPH_SCHEMA;
pub const CONSTRUCT_REASONING_STORE_SCHEMA: &str =
    gentle_protocol::CONSTRUCT_REASONING_STORE_SCHEMA;
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
const ENSEMBL_PROTEIN_ENTRIES_METADATA_KEY: &str = "ensembl_protein_entries";
const ENSEMBL_PROTEIN_ENTRIES_SCHEMA: &str = "gentle.ensembl_protein_entries.v1";
const UNIPROT_GENOME_PROJECTIONS_METADATA_KEY: &str = "uniprot_genome_projections";
const UNIPROT_GENOME_PROJECTIONS_SCHEMA: &str = "gentle.uniprot_genome_projections.v1";
const UNIPROT_GENOME_PROJECTION_SCHEMA: &str = "gentle.uniprot_genome_projection.v1";
const UNIPROT_FEATURE_CODING_DNA_QUERY_SCHEMA: &str = "gentle.uniprot_feature_coding_dna_query.v1";
pub const UNIPROT_PROJECTION_AUDIT_REPORTS_METADATA_KEY: &str = "uniprot_projection_audit_reports";
const UNIPROT_PROJECTION_AUDIT_REPORTS_SCHEMA: &str = "gentle.uniprot_projection_audit_reports.v1";
pub const UNIPROT_PROJECTION_AUDIT_REPORT_SCHEMA: &str = "gentle.uniprot_projection_audit.v1";
pub const UNIPROT_PROJECTION_AUDIT_PARITY_REPORTS_METADATA_KEY: &str =
    "uniprot_projection_audit_parity_reports";
const UNIPROT_PROJECTION_AUDIT_PARITY_REPORTS_SCHEMA: &str =
    "gentle.uniprot_projection_audit_parity_reports.v1";
pub const UNIPROT_PROJECTION_AUDIT_PARITY_REPORT_SCHEMA: &str =
    "gentle.uniprot_projection_audit_parity.v1";
pub const UNIPROT_ENSEMBL_LINKS_SCHEMA: &str = "gentle.uniprot_ensembl_links.v1";
pub const UNIPROT_PROJECTION_TRANSCRIPT_ACCOUNTING_SCHEMA: &str =
    "gentle.uniprot_projection_transcript_accounting.v1";
pub const UNIPROT_ENSEMBL_EXON_COMPARE_SCHEMA: &str = "gentle.uniprot_ensembl_exon_compare.v1";
pub const UNIPROT_ENSEMBL_PEPTIDE_COMPARE_SCHEMA: &str =
    "gentle.uniprot_ensembl_peptide_compare.v1";
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
const RNA_READ_GENE_SUPPORT_AUDIT_SCHEMA: &str = "gentle.rna_read_gene_support_audit.v1";
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
const FEATURE_BED_EXPORT_REPORT_SCHEMA: &str = "gentle.sequence_feature_bed_export.v1";
const RESTRICTION_SITE_SCAN_REPORT_SCHEMA: &str = "gentle.restriction_site_scan.v1";
const SEQUENCE_CONTEXT_VIEW_SCHEMA: &str = "gentle.sequence_context_view.v1";
const SEQUENCE_CONTEXT_BUNDLE_SCHEMA: &str = "gentle.sequence_context_bundle.v1";
const FEATURE_QUERY_DEFAULT_LIMIT: usize = 200;
const FEATURE_QUERY_MAX_LIMIT: usize = 10_000;
const TFBS_REGION_SUMMARY_SCHEMA: &str = "gentle.tfbs_region_summary.v1";
const TFBS_SCORE_TRACK_REPORT_SCHEMA: &str = "gentle.tfbs_score_tracks.v1";
const JASPAR_ENTRY_EXPERT_VIEW_SCHEMA: &str = "gentle.jaspar_entry_expert.v1";
const JASPAR_ENTRY_PRESENTATION_REPORT_SCHEMA: &str = "gentle.jaspar_entry_presentation.v1";
const JASPAR_REGISTRY_BENCHMARK_REPORT_SCHEMA: &str = "gentle.jaspar_registry_benchmark.v1";
const JASPAR_CATALOG_REPORT_SCHEMA: &str = "gentle.jaspar_catalog.v1";
const JASPAR_REMOTE_METADATA_SNAPSHOT_SCHEMA: &str = "gentle.jaspar_remote_metadata_snapshot.v1";
const VARIANT_PROMOTER_CONTEXT_SCHEMA: &str = "gentle.variant_promoter_context.v1";
const PROMOTER_REPORTER_CANDIDATES_SCHEMA: &str = "gentle.promoter_reporter_candidates.v1";
const TFBS_REGION_SUMMARY_DEFAULT_LIMIT: usize = 200;
const TFBS_REGION_SUMMARY_MAX_LIMIT: usize = 10_000;
pub(crate) const DEFAULT_PROMOTER_WINDOW_UPSTREAM_BP: usize = 1000;
pub(crate) const DEFAULT_PROMOTER_WINDOW_DOWNSTREAM_BP: usize = 200;
pub const DEFAULT_JASPAR_PRESENTATION_RANDOM_SEQUENCE_LENGTH_BP: usize = 10_000;
pub const DEFAULT_JASPAR_PRESENTATION_RANDOM_SEED: u64 = 0x4A_41_53_50_41_52_5F_31;
const DEFAULT_VARIANT_PROMOTER_TFBS_FOCUS_HALF_WINDOW_BP: usize = 100;
const DEFAULT_PROMOTER_REPORTER_RETAIN_DOWNSTREAM_FROM_TSS_BP: usize = 200;
const DEFAULT_PROMOTER_REPORTER_RETAIN_UPSTREAM_BEYOND_VARIANT_BP: usize = 500;
const DEFAULT_PROMOTER_REPORTER_MAX_CANDIDATES: usize = 5;
const ANNOTATE_PROMOTER_WINDOWS_GENERATED_TAG: &str = "annotate_promoter_windows";
const CONSTRUCT_REASONING_ANNOTATION_WRITEBACK_GENERATED_TAG: &str = "generated";
const CONSTRUCT_REASONING_ANNOTATION_WRITEBACK_GENERATED_BY: &str =
    "ConstructReasoningWriteAnnotation";
const LOW_COMPLEXITY_WINDOW_BP: usize = 96;
const LOW_COMPLEXITY_STEP_BP: usize = 24;
const LOW_COMPLEXITY_SCORE_THRESHOLD: f64 = 0.45;
const LOW_COMPLEXITY_HOMOPOLYMER_THRESHOLD_BP: usize = 9;
const HOMOPOLYMER_ANNOTATION_THRESHOLD_BP: usize = 10;
const TANDEM_REPEAT_MIN_TOTAL_BP: usize = 18;
const DIRECT_REPEAT_SEED_BP: usize = 12;
const DIRECT_REPEAT_MAX_CLUSTER_GAP_BP: usize = 512;
const INVERTED_REPEAT_SEED_BP: usize = 10;
const INVERTED_REPEAT_MAX_SPAN_BP: usize = 256;
const ALU_LIKE_MIN_BP: usize = 240;
const ALU_LIKE_MAX_BP: usize = 340;
const ALU_LIKE_POLYA_MIN_BP: usize = 12;

fn default_tfbs_region_summary_min_focus_occurrences() -> usize {
    1
}

fn default_tfbs_score_track_clip_negative() -> bool {
    true
}

fn default_tfbs_score_track_value_kind() -> TfbsScoreTrackValueKind {
    TfbsScoreTrackValueKind::LlrBits
}

fn default_promoter_window_upstream_bp() -> usize {
    DEFAULT_PROMOTER_WINDOW_UPSTREAM_BP
}

fn default_promoter_window_downstream_bp() -> usize {
    DEFAULT_PROMOTER_WINDOW_DOWNSTREAM_BP
}

fn default_jaspar_presentation_random_sequence_length_bp() -> usize {
    DEFAULT_JASPAR_PRESENTATION_RANDOM_SEQUENCE_LENGTH_BP
}

fn default_jaspar_presentation_random_seed() -> u64 {
    DEFAULT_JASPAR_PRESENTATION_RANDOM_SEED
}

fn default_variant_promoter_tfbs_focus_half_window_bp() -> usize {
    DEFAULT_VARIANT_PROMOTER_TFBS_FOCUS_HALF_WINDOW_BP
}

fn default_promoter_reporter_retain_downstream_from_tss_bp() -> usize {
    DEFAULT_PROMOTER_REPORTER_RETAIN_DOWNSTREAM_FROM_TSS_BP
}

fn default_promoter_reporter_retain_upstream_beyond_variant_bp() -> usize {
    DEFAULT_PROMOTER_REPORTER_RETAIN_UPSTREAM_BEYOND_VARIANT_BP
}

fn default_promoter_reporter_max_candidates() -> usize {
    DEFAULT_PROMOTER_REPORTER_MAX_CANDIDATES
}

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
#[path = "engine/analysis/jaspar.rs"]
mod jaspar;
#[path = "engine/state/lineage_containers.rs"]
mod lineage_containers;
#[path = "engine/analysis/motif_statistics.rs"]
mod motif_statistics;
#[path = "engine/ops/operation_handlers.rs"]
mod operation_handlers;
#[path = "engine/analysis/promoter_design.rs"]
mod promoter_design;
#[path = "engine/analysis/protein_handoff.rs"]
mod protein_handoff;
#[path = "engine/analysis/rna_reads.rs"]
mod rna_reads;
#[path = "engine/state/sequence_ops.rs"]
mod sequence_ops;
#[path = "engine/analysis/sequencing_confirmation.rs"]
mod sequencing_confirmation;
#[path = "engine/io/sequencing_traces.rs"]
mod sequencing_traces;
#[path = "engine/analysis/variant_promoter.rs"]
mod variant_promoter;

#[path = "engine/protocol.rs"]
pub mod protocol;
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

fn default_sequencing_primer_min_3prime_anneal_bp() -> usize {
    18
}

fn default_sequencing_primer_predicted_read_length_bp() -> usize {
    800
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

impl PrimerDesignBackend {
    pub fn as_str(self) -> &'static str {
        match self {
            Self::Auto => "auto",
            Self::Internal => "internal",
            Self::Primer3 => "primer3",
        }
    }
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

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct TfThresholdOverride {
    pub tf: String,
    pub min_llr_bits: Option<f64>,
    pub min_llr_quantile: Option<f64>,
}

// Backward-compatible alias kept while adapter/docs migrate to SequenceAnchor.
pub type AnchoredRegionAnchor = SequenceAnchor;

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

impl GuideU6TerminatorWindow {
    pub fn as_str(self) -> &'static str {
        match self {
            Self::SpacerOnly => "spacer_only",
            Self::SpacerPlusTail => "spacer_plus_tail",
        }
    }
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
    restriction_cloning_handoffs: HashMap<String, RestrictionCloningPcrHandoffReport>,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
struct ProteinDerivationReportStore {
    schema: String,
    updated_at_unix_ms: u128,
    reports: BTreeMap<String, ProteinDerivationReport>,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
struct UniprotProjectionAuditReportStore {
    schema: String,
    updated_at_unix_ms: u128,
    reports: BTreeMap<String, UniprotProjectionAuditReport>,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
struct UniprotProjectionAuditParityReportStore {
    schema: String,
    updated_at_unix_ms: u128,
    reports: BTreeMap<String, UniprotProjectionAuditParityReport>,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default)]
struct ReverseTranslationReportStore {
    schema: String,
    updated_at_unix_ms: u128,
    reports: BTreeMap<String, ReverseTranslationReport>,
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

impl PlanningSuggestionStatus {
    pub fn as_str(self) -> &'static str {
        match self {
            Self::Pending => "pending",
            Self::Accepted => "accepted",
            Self::Rejected => "rejected",
        }
    }
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
            helper_profile_id: None,
            preferred_routine_families: vec![],
            enforce_guardrails: true,
        }
    }
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

impl IsoformTranscriptGeometryMode {
    pub fn as_str(self) -> &'static str {
        match self {
            Self::Exon => "exon",
            Self::Cds => "cds",
        }
    }
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
struct EnsemblProteinEntryStore {
    schema: String,
    updated_at_unix_ms: u128,
    entries: HashMap<String, EnsemblProteinEntry>,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
#[serde(default, deny_unknown_fields)]
struct UniprotGenomeProjectionStore {
    schema: String,
    updated_at_unix_ms: u128,
    projections: HashMap<String, UniprotGenomeProjection>,
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
struct CandidateMacroTemplateStore {
    schema: String,
    updated_at_unix_ms: u128,
    templates: HashMap<String, CandidateMacroTemplate>,
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

impl GenomeAnchorSide {
    fn as_str(self) -> &'static str {
        match self {
            Self::FivePrime => "5prime",
            Self::ThreePrime => "3prime",
        }
    }
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

impl GenomeGeneExtractMode {
    pub fn as_str(self) -> &'static str {
        match self {
            Self::Gene => "gene",
            Self::CodingWithPromoter => "coding_with_promoter",
        }
    }
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
        #[serde(default)]
        overlay_x_axis_mode: DotplotOverlayXAxisMode,
        #[serde(default)]
        overlay_anchor_exon: Option<DotplotOverlayAnchorExonRef>,
    },
    RenderTfbsScoreTracksSvg {
        seq_id: SeqId,
        motifs: Vec<String>,
        start_0based: usize,
        end_0based_exclusive: usize,
        #[serde(default = "default_tfbs_score_track_value_kind")]
        score_kind: TfbsScoreTrackValueKind,
        clip_negative: bool,
        path: String,
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
        #[serde(default)]
        conditions: Option<GelRunConditions>,
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
    SetContainerDeclaredContentsExclusive {
        container_id: ContainerId,
        exclusive: bool,
    },
    CreateRackFromArrangement {
        arrangement_id: String,
        #[serde(default)]
        rack_id: Option<String>,
        #[serde(default)]
        name: Option<String>,
        #[serde(default)]
        profile: Option<RackProfileKind>,
    },
    PlaceArrangementOnRack {
        arrangement_id: String,
        rack_id: String,
    },
    MoveRackPlacement {
        rack_id: String,
        from_coordinate: String,
        to_coordinate: String,
        #[serde(default)]
        move_block: bool,
    },
    MoveRackSamples {
        rack_id: String,
        from_coordinates: Vec<String>,
        to_coordinate: String,
    },
    MoveRackArrangementBlocks {
        rack_id: String,
        arrangement_ids: Vec<String>,
        to_coordinate: String,
    },
    SetRackProfile {
        rack_id: String,
        profile: RackProfileKind,
    },
    ApplyRackTemplate {
        rack_id: String,
        template: RackAuthoringTemplate,
    },
    SetRackFillDirection {
        rack_id: String,
        fill_direction: RackFillDirection,
    },
    SetRackProfileCustom {
        rack_id: String,
        rows: usize,
        columns: usize,
    },
    SetRackBlockedCoordinates {
        rack_id: String,
        #[serde(default)]
        blocked_coordinates: Vec<String>,
    },
    ExportRackLabelsSvg {
        rack_id: String,
        path: String,
        #[serde(default)]
        arrangement_id: Option<String>,
        #[serde(default)]
        preset: RackLabelSheetPreset,
    },
    ExportRackFabricationSvg {
        rack_id: String,
        path: String,
        #[serde(default)]
        template: RackPhysicalTemplateKind,
    },
    ExportRackIsometricSvg {
        rack_id: String,
        path: String,
        #[serde(default)]
        template: RackPhysicalTemplateKind,
    },
    ExportRackOpenScad {
        rack_id: String,
        path: String,
        #[serde(default)]
        template: RackPhysicalTemplateKind,
    },
    ExportRackCarrierLabelsSvg {
        rack_id: String,
        path: String,
        #[serde(default)]
        arrangement_id: Option<String>,
        #[serde(default)]
        template: RackPhysicalTemplateKind,
        #[serde(default)]
        preset: RackCarrierLabelPreset,
    },
    ExportRackSimulationJson {
        rack_id: String,
        path: String,
        #[serde(default)]
        template: RackPhysicalTemplateKind,
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
    ExtractGenomePromoterSlice {
        genome_id: String,
        gene_query: String,
        occurrence: Option<usize>,
        #[serde(default, skip_serializing_if = "Option::is_none")]
        transcript_id: Option<String>,
        output_id: Option<SeqId>,
        #[serde(default = "default_promoter_window_upstream_bp")]
        upstream_bp: usize,
        #[serde(default = "default_promoter_window_downstream_bp")]
        downstream_bp: usize,
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
    FetchEnsemblProtein {
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
    ImportEnsemblProteinSequence {
        entry_id: String,
        output_id: Option<SeqId>,
    },
    ProjectUniprotToGenome {
        seq_id: SeqId,
        entry_id: String,
        projection_id: Option<String>,
        transcript_id: Option<String>,
    },
    AuditUniprotProjectionConsistency {
        projection_id: String,
        #[serde(default)]
        transcript_id: Option<String>,
        #[serde(default)]
        report_id: Option<String>,
        #[serde(default)]
        ensembl_entry_id: Option<String>,
    },
    AuditUniprotProjectionParity {
        projection_id: String,
        #[serde(default)]
        transcript_id: Option<String>,
        #[serde(default)]
        report_id: Option<String>,
        #[serde(default)]
        ensembl_entry_id: Option<String>,
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
    FindRestrictionSites {
        target: SequenceScanTarget,
        #[serde(default)]
        enzymes: Vec<String>,
        #[serde(default, skip_serializing_if = "Option::is_none")]
        max_sites_per_enzyme: Option<usize>,
        #[serde(default = "default_true")]
        include_cut_geometry: bool,
        #[serde(default, skip_serializing_if = "Option::is_none")]
        path: Option<String>,
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
    PrepareRestrictionCloningPcrHandoff {
        template: SeqId,
        primer_report_id: String,
        pair_index: usize,
        destination_vector_seq_id: SeqId,
        mode: RestrictionCloningPcrHandoffMode,
        forward_enzyme: String,
        reverse_enzyme: Option<String>,
        forward_leader_5prime: Option<String>,
        reverse_leader_5prime: Option<String>,
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
    DeriveProteinSequences {
        seq_id: SeqId,
        #[serde(default)]
        feature_ids: Vec<usize>,
        #[serde(default)]
        scope: Option<SplicingScopePreset>,
        #[serde(default)]
        output_prefix: Option<String>,
    },
    ReverseTranslateProteinSequence {
        seq_id: SeqId,
        #[serde(default)]
        output_id: Option<SeqId>,
        #[serde(default)]
        speed_profile: Option<TranslationSpeedProfile>,
        #[serde(default)]
        speed_mark: Option<TranslationSpeedMark>,
        #[serde(default)]
        translation_table: Option<usize>,
        #[serde(default)]
        target_anneal_tm_c: Option<f64>,
        #[serde(default)]
        anneal_window_bp: Option<usize>,
    },
    BuildProteinToDnaHandoffReasoning {
        seq_id: SeqId,
        protein_seq_id: SeqId,
        #[serde(default)]
        transcript_filter: Option<String>,
        #[serde(default)]
        projection_id: Option<String>,
        #[serde(default)]
        ensembl_entry_id: Option<String>,
        #[serde(default)]
        feature_query: Option<String>,
        #[serde(default)]
        ranking_goal: ProteinToDnaHandoffRankingGoal,
        #[serde(default)]
        speed_profile: Option<TranslationSpeedProfile>,
        #[serde(default)]
        speed_mark: Option<TranslationSpeedMark>,
        #[serde(default)]
        translation_table: Option<usize>,
        #[serde(default)]
        target_anneal_tm_c: Option<f64>,
        #[serde(default)]
        anneal_window_bp: Option<usize>,
        #[serde(default)]
        objective_id: Option<String>,
        #[serde(default)]
        graph_id: Option<String>,
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
    SuggestSequencingPrimers {
        expected_seq_id: SeqId,
        #[serde(default)]
        primer_seq_ids: Vec<SeqId>,
        #[serde(default)]
        confirmation_report_id: Option<String>,
        #[serde(default = "default_sequencing_primer_min_3prime_anneal_bp")]
        min_3prime_anneal_bp: usize,
        #[serde(default = "default_sequencing_primer_predicted_read_length_bp")]
        predicted_read_length_bp: usize,
    },
    ConfirmConstructReads {
        expected_seq_id: SeqId,
        #[serde(default)]
        baseline_seq_id: Option<SeqId>,
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
    InspectRnaReadGeneSupport {
        report_id: String,
        gene_ids: Vec<String>,
        #[serde(default, skip_serializing_if = "Vec::is_empty")]
        selected_record_indices: Vec<usize>,
        #[serde(default)]
        complete_rule: RnaReadGeneSupportCompleteRule,
        #[serde(default)]
        cohort_filter: RnaReadGeneSupportAuditCohortFilter,
        #[serde(default, skip_serializing_if = "Option::is_none")]
        path: Option<String>,
    },
    SummarizeTfbsRegion {
        seq_id: String,
        focus_start_0based: usize,
        focus_end_0based_exclusive: usize,
        #[serde(default, skip_serializing_if = "Option::is_none")]
        context_start_0based: Option<usize>,
        #[serde(default, skip_serializing_if = "Option::is_none")]
        context_end_0based_exclusive: Option<usize>,
        #[serde(default = "default_tfbs_region_summary_min_focus_occurrences")]
        min_focus_occurrences: usize,
        #[serde(default)]
        min_context_occurrences: usize,
        #[serde(default, skip_serializing_if = "Option::is_none")]
        limit: Option<usize>,
        #[serde(default, skip_serializing_if = "Option::is_none")]
        path: Option<String>,
    },
    SummarizeTfbsScoreTracks {
        seq_id: String,
        motifs: Vec<String>,
        start_0based: usize,
        end_0based_exclusive: usize,
        #[serde(default = "default_tfbs_score_track_value_kind")]
        score_kind: TfbsScoreTrackValueKind,
        #[serde(default = "default_tfbs_score_track_clip_negative")]
        clip_negative: bool,
        #[serde(default, skip_serializing_if = "Option::is_none")]
        path: Option<String>,
    },
    SummarizeJasparEntries {
        #[serde(default, skip_serializing_if = "Vec::is_empty")]
        motifs: Vec<String>,
        #[serde(default = "default_jaspar_presentation_random_sequence_length_bp")]
        random_sequence_length_bp: usize,
        #[serde(default = "default_jaspar_presentation_random_seed")]
        random_seed: u64,
        #[serde(default, skip_serializing_if = "Option::is_none")]
        path: Option<String>,
    },
    BenchmarkJasparRegistry {
        #[serde(default = "default_jaspar_presentation_random_sequence_length_bp")]
        random_sequence_length_bp: usize,
        #[serde(default = "default_jaspar_presentation_random_seed")]
        random_seed: u64,
        #[serde(default, skip_serializing_if = "Option::is_none")]
        path: Option<String>,
    },
    ListJasparCatalog {
        #[serde(default, skip_serializing_if = "Option::is_none")]
        filter: Option<String>,
        #[serde(default, skip_serializing_if = "Option::is_none")]
        limit: Option<usize>,
        #[serde(default)]
        include_remote_metadata: bool,
        #[serde(default)]
        refresh_remote_metadata: bool,
        #[serde(default, skip_serializing_if = "Option::is_none")]
        path: Option<String>,
    },
    SyncJasparRemoteMetadata {
        #[serde(default, skip_serializing_if = "Vec::is_empty")]
        motifs: Vec<String>,
        #[serde(default, skip_serializing_if = "Option::is_none")]
        filter: Option<String>,
        #[serde(default, skip_serializing_if = "Option::is_none")]
        limit: Option<usize>,
        #[serde(default, skip_serializing_if = "Option::is_none")]
        path: Option<String>,
    },
    InspectJasparEntry {
        motif: String,
        #[serde(default = "default_jaspar_presentation_random_sequence_length_bp")]
        random_sequence_length_bp: usize,
        #[serde(default = "default_jaspar_presentation_random_seed")]
        random_seed: u64,
        #[serde(default)]
        include_remote_metadata: bool,
        #[serde(default)]
        refresh_remote_metadata: bool,
        #[serde(default, skip_serializing_if = "Option::is_none")]
        path: Option<String>,
    },
    AnnotatePromoterWindows {
        input: SeqId,
        #[serde(default, skip_serializing_if = "Option::is_none")]
        gene_label: Option<String>,
        #[serde(default, skip_serializing_if = "Option::is_none")]
        transcript_id: Option<String>,
        #[serde(default = "default_promoter_window_upstream_bp")]
        upstream_bp: usize,
        #[serde(default = "default_promoter_window_downstream_bp")]
        downstream_bp: usize,
        #[serde(default)]
        collapse_mode: PromoterWindowCollapseMode,
    },
    SummarizeVariantPromoterContext {
        input: SeqId,
        #[serde(default, skip_serializing_if = "Option::is_none")]
        variant_label_or_id: Option<String>,
        #[serde(default, skip_serializing_if = "Option::is_none")]
        gene_label: Option<String>,
        #[serde(default, skip_serializing_if = "Option::is_none")]
        transcript_id: Option<String>,
        #[serde(default = "default_promoter_window_upstream_bp")]
        promoter_upstream_bp: usize,
        #[serde(default = "default_promoter_window_downstream_bp")]
        promoter_downstream_bp: usize,
        #[serde(default = "default_variant_promoter_tfbs_focus_half_window_bp")]
        tfbs_focus_half_window_bp: usize,
        #[serde(default, skip_serializing_if = "Option::is_none")]
        path: Option<String>,
    },
    SuggestPromoterReporterFragments {
        input: SeqId,
        #[serde(default, skip_serializing_if = "Option::is_none")]
        variant_label_or_id: Option<String>,
        #[serde(default, skip_serializing_if = "Option::is_none")]
        gene_label: Option<String>,
        #[serde(default, skip_serializing_if = "Option::is_none")]
        transcript_id: Option<String>,
        #[serde(default = "default_promoter_reporter_retain_downstream_from_tss_bp")]
        retain_downstream_from_tss_bp: usize,
        #[serde(default = "default_promoter_reporter_retain_upstream_beyond_variant_bp")]
        retain_upstream_beyond_variant_bp: usize,
        #[serde(default = "default_promoter_reporter_max_candidates")]
        max_candidates: usize,
        #[serde(default, skip_serializing_if = "Option::is_none")]
        path: Option<String>,
    },
    MaterializeVariantAllele {
        input: SeqId,
        #[serde(default, skip_serializing_if = "Option::is_none")]
        variant_label_or_id: Option<String>,
        #[serde(default)]
        allele: VariantAlleleChoice,
        #[serde(default, skip_serializing_if = "Option::is_none")]
        output_id: Option<SeqId>,
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
        #[serde(default, skip_serializing_if = "Vec::is_empty")]
        gene_ids: Vec<String>,
        #[serde(default)]
        complete_rule: RnaReadGeneSupportCompleteRule,
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
    ExportFeaturesBed {
        query: SequenceFeatureQuery,
        path: String,
        #[serde(default)]
        coordinate_mode: Option<FeatureBedCoordinateMode>,
        #[serde(default)]
        include_restriction_sites: Option<bool>,
        #[serde(default)]
        restriction_enzymes: Vec<String>,
    },
    InspectSequenceContextView {
        seq_id: String,
        #[serde(default, skip_serializing_if = "Option::is_none")]
        mode: Option<RenderSvgMode>,
        #[serde(default, skip_serializing_if = "Option::is_none")]
        viewport_start_0based: Option<usize>,
        #[serde(default, skip_serializing_if = "Option::is_none")]
        viewport_end_0based_exclusive: Option<usize>,
        #[serde(default, skip_serializing_if = "Vec::is_empty")]
        include_visible_classes: Vec<String>,
        #[serde(default, skip_serializing_if = "Option::is_none")]
        coordinate_mode: Option<FeatureBedCoordinateMode>,
        #[serde(default, skip_serializing_if = "Option::is_none")]
        limit: Option<usize>,
    },
    ExportSequenceContextBundle {
        seq_id: String,
        #[serde(default, skip_serializing_if = "Option::is_none")]
        mode: Option<RenderSvgMode>,
        #[serde(default, skip_serializing_if = "Option::is_none")]
        viewport_start_0based: Option<usize>,
        #[serde(default, skip_serializing_if = "Option::is_none")]
        viewport_end_0based_exclusive: Option<usize>,
        #[serde(default, skip_serializing_if = "Option::is_none")]
        coordinate_mode: Option<FeatureBedCoordinateMode>,
        #[serde(default, skip_serializing_if = "Option::is_none")]
        include_feature_bed: Option<bool>,
        #[serde(default, skip_serializing_if = "Option::is_none")]
        include_text_summary: Option<bool>,
        #[serde(default, skip_serializing_if = "Option::is_none")]
        include_restriction_sites: Option<bool>,
        #[serde(default, skip_serializing_if = "Vec::is_empty")]
        restriction_enzymes: Vec<String>,
        output_dir: String,
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
struct FeatureBedExportRow {
    chrom: String,
    chrom_start_0based: usize,
    chrom_end_0based_exclusive: usize,
    name: String,
    score_0_to_1000: usize,
    strand: char,
    kind: String,
    row_id: String,
    coordinate_source: &'static str,
    qualifiers_json: String,
    sort_feature_id: usize,
    length_bp: usize,
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
            .or_else(|| anchor.catalog_path.clone());
        let resolved_cache_dir = cache_dir
            .map(str::trim)
            .filter(|v| !v.is_empty())
            .map(|v| v.to_string())
            .or_else(|| anchor.cache_dir.clone());
        let (catalog, _) = Self::open_reference_genome_catalog(requested_catalog_path.as_deref())?;
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
                "RenderTfbsScoreTracksSvg".to_string(),
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
                "SetContainerDeclaredContentsExclusive".to_string(),
                "CreateRackFromArrangement".to_string(),
                "PlaceArrangementOnRack".to_string(),
                "MoveRackPlacement".to_string(),
                "MoveRackSamples".to_string(),
                "MoveRackArrangementBlocks".to_string(),
                "SetRackProfile".to_string(),
                "ApplyRackTemplate".to_string(),
                "SetRackFillDirection".to_string(),
                "SetRackProfileCustom".to_string(),
                "SetRackBlockedCoordinates".to_string(),
                "ExportRackLabelsSvg".to_string(),
                "ExportRackFabricationSvg".to_string(),
                "ExportRackIsometricSvg".to_string(),
                "ExportRackOpenScad".to_string(),
                "ExportRackCarrierLabelsSvg".to_string(),
                "ExportRackSimulationJson".to_string(),
                "ExportDnaLadders".to_string(),
                "ExportRnaLadders".to_string(),
                "ExportPool".to_string(),
                "ExportProcessRunBundle".to_string(),
                "PrepareGenome".to_string(),
                "ExtractGenomeRegion".to_string(),
                "ExtractGenomeGene".to_string(),
                "ExtractGenomePromoterSlice".to_string(),
                "ExtendGenomeAnchor".to_string(),
                "ImportGenomeBedTrack".to_string(),
                "ImportGenomeBigWigTrack".to_string(),
                "ImportGenomeVcfTrack".to_string(),
                "ImportIsoformPanel".to_string(),
                "ImportUniprotSwissProt".to_string(),
                "FetchUniprotSwissProt".to_string(),
                "FetchEnsemblProtein".to_string(),
                "FetchGenBankAccession".to_string(),
                "FetchDbSnpRegion".to_string(),
                "FetchUniprotLinkedGenBank".to_string(),
                "ImportUniprotEntrySequence".to_string(),
                "ImportEnsemblProteinSequence".to_string(),
                "ProjectUniprotToGenome".to_string(),
                "AuditUniprotProjectionConsistency".to_string(),
                "AuditUniprotProjectionParity".to_string(),
                "ImportBlastHitsTrack".to_string(),
                "DigestContainer".to_string(),
                "MergeContainersById".to_string(),
                "LigationContainer".to_string(),
                "FilterContainerByMolecularWeight".to_string(),
                "Digest".to_string(),
                "FindRestrictionSites".to_string(),
                "MergeContainers".to_string(),
                "Ligation".to_string(),
                "Pcr".to_string(),
                "PcrAdvanced".to_string(),
                "PcrMutagenesis".to_string(),
                "DesignPrimerPairs".to_string(),
                "DesignInsertionPrimerPairs".to_string(),
                "PrepareRestrictionCloningPcrHandoff".to_string(),
                "PcrOverlapExtensionMutagenesis".to_string(),
                "DesignQpcrAssays".to_string(),
                "DeriveTranscriptSequences".to_string(),
                "DeriveProteinSequences".to_string(),
                "ReverseTranslateProteinSequence".to_string(),
                "BuildProteinToDnaHandoffReasoning".to_string(),
                "ComputeDotplot".to_string(),
                "ComputeDotplotOverlay".to_string(),
                "ComputeFlexibilityTrack".to_string(),
                "DeriveSplicingReferences".to_string(),
                "AlignSequences".to_string(),
                "ConfirmConstructReads".to_string(),
                "SuggestSequencingPrimers".to_string(),
                "ListSequencingConfirmationReports".to_string(),
                "ShowSequencingConfirmationReport".to_string(),
                "ExportSequencingConfirmationReport".to_string(),
                "ExportSequencingConfirmationSupportTsv".to_string(),
                "InterpretRnaReads".to_string(),
                "AlignRnaReadReport".to_string(),
                "ListRnaReadReports".to_string(),
                "ShowRnaReadReport".to_string(),
                "SummarizeRnaReadGeneSupport".to_string(),
                "InspectRnaReadGeneSupport".to_string(),
                "SummarizeTfbsRegion".to_string(),
                "SummarizeTfbsScoreTracks".to_string(),
                "InspectJasparEntry".to_string(),
                "SummarizeJasparEntries".to_string(),
                "BenchmarkJasparRegistry".to_string(),
                "ListJasparCatalog".to_string(),
                "SyncJasparRemoteMetadata".to_string(),
                "AnnotatePromoterWindows".to_string(),
                "SummarizeVariantPromoterContext".to_string(),
                "SuggestPromoterReporterFragments".to_string(),
                "MaterializeVariantAllele".to_string(),
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
                "ExportFeaturesBed".to_string(),
                "InspectSequenceContextView".to_string(),
                "ExportSequenceContextBundle".to_string(),
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

    fn open_catalog_with_default_mode(
        catalog_path: Option<&str>,
        helper_mode: bool,
    ) -> Result<(GenomeCatalog, String), EngineError> {
        let requested_catalog = catalog_path
            .map(str::trim)
            .filter(|v| !v.is_empty())
            .map(|v| v.to_string())
            .unwrap_or_default();
        let catalog = match requested_catalog.as_str() {
            "" => GenomeCatalog::from_default_discovery(helper_mode),
            DEFAULT_REFERENCE_CATALOG_DISCOVERY_TOKEN => {
                GenomeCatalog::from_default_discovery(false)
            }
            DEFAULT_HELPER_CATALOG_DISCOVERY_TOKEN => GenomeCatalog::from_default_discovery(true),
            path => GenomeCatalog::from_json_file(path),
        }
        .map_err(|e| EngineError {
            code: ErrorCode::InvalidInput,
            message: if requested_catalog.is_empty() {
                format!("Could not open default genome catalog discovery: {e}")
            } else {
                format!("Could not open genome catalog '{}': {e}", requested_catalog)
            },
        })?;
        let origin = catalog.catalog_origin_label().to_string();
        Ok((catalog, origin))
    }

    fn open_reference_genome_catalog(
        catalog_path: Option<&str>,
    ) -> Result<(GenomeCatalog, String), EngineError> {
        Self::open_catalog_with_default_mode(catalog_path, false)
    }

    fn open_helper_genome_catalog(
        catalog_path: Option<&str>,
    ) -> Result<(GenomeCatalog, String), EngineError> {
        Self::open_catalog_with_default_mode(catalog_path, true)
    }

    fn open_host_profile_catalog(
        catalog_path: Option<&str>,
    ) -> Result<(gentle_protocol::HostProfileCatalog, String), EngineError> {
        let resolved_catalog_path = catalog_path
            .map(str::trim)
            .filter(|value| !value.is_empty())
            .unwrap_or(DEFAULT_HOST_PROFILE_CATALOG_PATH)
            .to_string();
        let raw = std::fs::read_to_string(&resolved_catalog_path).map_err(|e| EngineError {
            code: ErrorCode::Io,
            message: format!(
                "Could not read host profile catalog '{}': {}",
                resolved_catalog_path, e
            ),
        })?;
        let catalog =
            serde_json::from_str::<gentle_protocol::HostProfileCatalog>(&raw).map_err(|e| {
                EngineError {
                    code: ErrorCode::InvalidInput,
                    message: format!(
                        "Could not parse host profile catalog '{}': {}",
                        resolved_catalog_path, e
                    ),
                }
            })?;
        Ok((catalog, resolved_catalog_path))
    }

    fn host_profile_matches_filter(profile: &HostProfileRecord, filter: Option<&str>) -> bool {
        let tokens = filter
            .map(str::trim)
            .filter(|value| !value.is_empty())
            .map(|value| {
                value
                    .split_whitespace()
                    .map(|token| token.trim().to_ascii_lowercase())
                    .filter(|token| !token.is_empty())
                    .collect::<Vec<_>>()
            })
            .unwrap_or_default();
        if tokens.is_empty() {
            return true;
        }
        let mut haystack = vec![
            profile.profile_id.as_str(),
            profile.species.as_str(),
            profile.strain.as_str(),
        ]
        .join("\n");
        for value in profile.aliases.iter().chain(profile.genotype_tags.iter()) {
            haystack.push('\n');
            haystack.push_str(value);
        }
        for value in profile
            .phenotype_tags
            .iter()
            .chain(profile.notes.iter())
            .chain(profile.source_notes.iter())
        {
            haystack.push('\n');
            haystack.push_str(value);
        }
        let haystack = haystack.to_ascii_lowercase();
        tokens.iter().all(|token| haystack.contains(token))
    }

    pub fn list_host_profile_catalog_entries(
        catalog_path: Option<&str>,
        filter: Option<&str>,
    ) -> Result<Vec<HostProfileRecord>, EngineError> {
        let (catalog, _) = Self::open_host_profile_catalog(catalog_path)?;
        let mut profiles = catalog
            .profiles
            .into_iter()
            .filter(|profile| Self::host_profile_matches_filter(profile, filter))
            .collect::<Vec<_>>();
        profiles.sort_by(|left, right| {
            left.profile_id
                .cmp(&right.profile_id)
                .then_with(|| left.species.cmp(&right.species))
                .then_with(|| left.strain.cmp(&right.strain))
        });
        Ok(profiles)
    }

    pub fn list_reference_genomes(catalog_path: Option<&str>) -> Result<Vec<String>, EngineError> {
        let (catalog, _) = Self::open_reference_genome_catalog(catalog_path)?;
        Ok(catalog.list_genomes())
    }

    pub fn list_reference_catalog_entries(
        catalog_path: Option<&str>,
        filter: Option<&str>,
    ) -> Result<Vec<GenomeCatalogListEntry>, EngineError> {
        let (catalog, _) = Self::open_reference_genome_catalog(catalog_path)?;
        Ok(catalog.list_entries(filter))
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

    pub fn preview_helper_genome_ensembl_catalog_updates(
        catalog_path: Option<&str>,
    ) -> Result<EnsemblCatalogUpdatePreview, EngineError> {
        let (catalog, catalog_path) = Self::open_helper_genome_catalog(catalog_path)?;
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

    pub fn discover_ensembl_installable_genomes(
        collection_filter: Option<&str>,
        filter: Option<&str>,
    ) -> Result<EnsemblInstallableGenomeCatalog, EngineError> {
        GenomeCatalog::discover_ensembl_installable_genomes(collection_filter, filter).map_err(
            |e| EngineError {
                code: ErrorCode::Io,
                message: format!("Could not discover installable Ensembl genomes: {e}"),
            },
        )
    }

    fn preview_ensembl_quick_install_with_scope(
        helper_mode: bool,
        catalog_path: Option<&str>,
        collection: &str,
        species_dir: &str,
        output_catalog_path: Option<&str>,
        genome_id: Option<&str>,
    ) -> Result<EnsemblQuickInstallPreview, EngineError> {
        let (catalog, catalog_origin) =
            Self::open_catalog_with_default_mode(catalog_path, helper_mode)?;
        catalog
            .preview_ensembl_quick_install_for_domain(
                helper_mode,
                collection,
                species_dir,
                output_catalog_path,
                genome_id,
            )
            .map_err(|e| EngineError {
                code: ErrorCode::Io,
                message: format!(
                    "Could not preview Ensembl quick install for '{}' from '{}': {}",
                    species_dir, catalog_origin, e
                ),
            })
    }

    fn apply_ensembl_quick_install_with_scope(
        helper_mode: bool,
        catalog_path: Option<&str>,
        collection: &str,
        species_dir: &str,
        output_catalog_path: Option<&str>,
        genome_id: Option<&str>,
    ) -> Result<EnsemblQuickInstallCatalogWriteReport, EngineError> {
        let (catalog, catalog_origin) =
            Self::open_catalog_with_default_mode(catalog_path, helper_mode)?;
        catalog
            .apply_ensembl_quick_install_for_domain(
                helper_mode,
                collection,
                species_dir,
                output_catalog_path,
                genome_id,
            )
            .map_err(|e| EngineError {
                code: ErrorCode::Io,
                message: format!(
                    "Could not apply Ensembl quick install for '{}' from '{}': {}",
                    species_dir, catalog_origin, e
                ),
            })
    }

    fn quick_install_ensembl_with_scope(
        helper_mode: bool,
        catalog_path: Option<&str>,
        collection: &str,
        species_dir: &str,
        output_catalog_path: Option<&str>,
        genome_id: Option<&str>,
        cache_dir: Option<&str>,
        timeout_seconds: Option<u64>,
        on_progress: &mut dyn FnMut(PrepareGenomeProgress) -> bool,
    ) -> Result<EnsemblQuickInstallReport, EngineError> {
        let (catalog, catalog_origin) =
            Self::open_catalog_with_default_mode(catalog_path, helper_mode)?;
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
        catalog
            .quick_install_ensembl_genome_once_with_progress_for_domain(
                helper_mode,
                collection,
                species_dir,
                output_catalog_path,
                genome_id,
                cache_dir,
                &mut guarded_progress,
            )
            .map_err(|e| EngineError {
                code: ErrorCode::Io,
                message: if is_prepare_cancelled_error(&e) {
                    if timed_out {
                        format!(
                            "Ensembl quick install timed out for '{}' after {} second(s)",
                            species_dir,
                            timeout_seconds.unwrap_or(0)
                        )
                    } else {
                        format!("Ensembl quick install cancelled for '{}'", species_dir)
                    }
                } else {
                    format!(
                        "Could not quick-install Ensembl species '{}' from '{}': {}",
                        species_dir, catalog_origin, e
                    )
                },
            })
    }

    pub fn preview_reference_genome_ensembl_quick_install(
        catalog_path: Option<&str>,
        collection: &str,
        species_dir: &str,
        output_catalog_path: Option<&str>,
        genome_id: Option<&str>,
    ) -> Result<EnsemblQuickInstallPreview, EngineError> {
        Self::preview_ensembl_quick_install_with_scope(
            false,
            catalog_path,
            collection,
            species_dir,
            output_catalog_path,
            genome_id,
        )
    }

    pub fn preview_helper_genome_ensembl_quick_install(
        catalog_path: Option<&str>,
        collection: &str,
        species_dir: &str,
        output_catalog_path: Option<&str>,
        genome_id: Option<&str>,
    ) -> Result<EnsemblQuickInstallPreview, EngineError> {
        Self::preview_ensembl_quick_install_with_scope(
            true,
            catalog_path,
            collection,
            species_dir,
            output_catalog_path,
            genome_id,
        )
    }

    pub fn apply_reference_genome_ensembl_quick_install(
        catalog_path: Option<&str>,
        collection: &str,
        species_dir: &str,
        output_catalog_path: Option<&str>,
        genome_id: Option<&str>,
    ) -> Result<EnsemblQuickInstallCatalogWriteReport, EngineError> {
        Self::apply_ensembl_quick_install_with_scope(
            false,
            catalog_path,
            collection,
            species_dir,
            output_catalog_path,
            genome_id,
        )
    }

    pub fn apply_helper_genome_ensembl_quick_install(
        catalog_path: Option<&str>,
        collection: &str,
        species_dir: &str,
        output_catalog_path: Option<&str>,
        genome_id: Option<&str>,
    ) -> Result<EnsemblQuickInstallCatalogWriteReport, EngineError> {
        Self::apply_ensembl_quick_install_with_scope(
            true,
            catalog_path,
            collection,
            species_dir,
            output_catalog_path,
            genome_id,
        )
    }

    pub fn quick_install_reference_genome_from_ensembl(
        catalog_path: Option<&str>,
        collection: &str,
        species_dir: &str,
        output_catalog_path: Option<&str>,
        genome_id: Option<&str>,
        cache_dir: Option<&str>,
        timeout_seconds: Option<u64>,
        on_progress: &mut dyn FnMut(PrepareGenomeProgress) -> bool,
    ) -> Result<EnsemblQuickInstallReport, EngineError> {
        Self::quick_install_ensembl_with_scope(
            false,
            catalog_path,
            collection,
            species_dir,
            output_catalog_path,
            genome_id,
            cache_dir,
            timeout_seconds,
            on_progress,
        )
    }

    pub fn quick_install_helper_genome_from_ensembl(
        catalog_path: Option<&str>,
        collection: &str,
        species_dir: &str,
        output_catalog_path: Option<&str>,
        genome_id: Option<&str>,
        cache_dir: Option<&str>,
        timeout_seconds: Option<u64>,
        on_progress: &mut dyn FnMut(PrepareGenomeProgress) -> bool,
    ) -> Result<EnsemblQuickInstallReport, EngineError> {
        Self::quick_install_ensembl_with_scope(
            true,
            catalog_path,
            collection,
            species_dir,
            output_catalog_path,
            genome_id,
            cache_dir,
            timeout_seconds,
            on_progress,
        )
    }

    pub fn apply_helper_genome_ensembl_catalog_updates(
        catalog_path: Option<&str>,
        output_catalog_path: Option<&str>,
    ) -> Result<EnsemblCatalogUpdateReport, EngineError> {
        let (catalog, catalog_path) = Self::open_helper_genome_catalog(catalog_path)?;
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

    pub fn remove_helper_genome_catalog_entry(
        catalog_path: Option<&str>,
        genome_id: &str,
        output_catalog_path: Option<&str>,
    ) -> Result<GenomeCatalogEntryRemovalReport, EngineError> {
        let (catalog, catalog_path) = Self::open_helper_genome_catalog(catalog_path)?;
        catalog
            .remove_catalog_entry(genome_id, output_catalog_path)
            .map_err(|e| EngineError {
                code: ErrorCode::Io,
                message: format!(
                    "Could not remove helper catalog entry '{}' from '{}': {}",
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

    pub fn remove_prepared_helper_genome(
        catalog_path: Option<&str>,
        genome_id: &str,
        cache_dir: Option<&str>,
    ) -> Result<PreparedGenomeRemovalReport, EngineError> {
        let (catalog, catalog_path) = Self::open_helper_genome_catalog(catalog_path)?;
        catalog
            .remove_prepared_genome_install(
                genome_id,
                cache_dir.map(str::trim).filter(|v| !v.is_empty()),
            )
            .map_err(|e| EngineError {
                code: ErrorCode::Io,
                message: format!(
                    "Could not remove prepared helper '{}' from '{}': {}",
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
        let (catalog, _) = Self::open_helper_genome_catalog(catalog_path)?;
        Ok(catalog.list_genomes())
    }

    pub fn list_helper_catalog_entries(
        catalog_path: Option<&str>,
        filter: Option<&str>,
    ) -> Result<Vec<GenomeCatalogListEntry>, EngineError> {
        let (catalog, _) = Self::open_helper_genome_catalog(catalog_path)?;
        Ok(catalog.list_entries(filter))
    }

    pub fn describe_helper_genome_sources(
        genome_id: &str,
        catalog_path: Option<&str>,
        cache_dir: Option<&str>,
    ) -> Result<GenomeSourcePlan, EngineError> {
        let (catalog, catalog_path) = Self::open_helper_genome_catalog(catalog_path)?;
        catalog
            .source_plan(
                genome_id,
                cache_dir.map(str::trim).filter(|v| !v.is_empty()),
            )
            .map_err(|e| EngineError {
                code: ErrorCode::InvalidInput,
                message: format!(
                    "Could not resolve source plan for helper '{}' in catalog '{}': {}",
                    genome_id, catalog_path, e
                ),
            })
    }

    pub fn interpret_helper_genome(
        genome_id: &str,
        catalog_path: Option<&str>,
    ) -> Result<Option<HelperConstructInterpretation>, EngineError> {
        let (catalog, catalog_path) = Self::open_helper_genome_catalog(catalog_path)?;
        catalog
            .helper_construct_interpretation(genome_id)
            .map_err(|e| EngineError {
                code: ErrorCode::InvalidInput,
                message: format!(
                    "Could not interpret helper '{}' in catalog '{}': {}",
                    genome_id, catalog_path, e
                ),
            })
    }

    pub fn resolve_reference_genome_cache_dir(
        catalog_path: Option<&str>,
        genome_id: &str,
        cache_dir: Option<&str>,
    ) -> Result<String, EngineError> {
        let (catalog, catalog_path) = Self::open_reference_genome_catalog(catalog_path)?;
        catalog
            .effective_cache_dir(
                genome_id,
                cache_dir.map(str::trim).filter(|v| !v.is_empty()),
            )
            .map_err(|e| EngineError {
                code: ErrorCode::InvalidInput,
                message: format!(
                    "Could not resolve cache dir for genome '{}' in catalog '{}': {}",
                    genome_id, catalog_path, e
                ),
            })
    }

    pub fn resolve_helper_genome_cache_dir(
        genome_id: &str,
        catalog_path: Option<&str>,
        cache_dir: Option<&str>,
    ) -> Result<String, EngineError> {
        let (catalog, catalog_path) = Self::open_helper_genome_catalog(catalog_path)?;
        catalog
            .effective_cache_dir(
                genome_id,
                cache_dir.map(str::trim).filter(|v| !v.is_empty()),
            )
            .map_err(|e| EngineError {
                code: ErrorCode::InvalidInput,
                message: format!(
                    "Could not resolve cache dir for helper '{}' in catalog '{}': {}",
                    genome_id, catalog_path, e
                ),
            })
    }

    pub fn inspect_reference_genome_prepared_compatibility(
        catalog_path: Option<&str>,
        genome_id: &str,
        cache_dir: Option<&str>,
    ) -> Result<PreparedGenomeCompatibilityInspection, EngineError> {
        let (catalog, catalog_path) = Self::open_reference_genome_catalog(catalog_path)?;
        catalog
            .inspect_prepared_genome_compatibility(
                genome_id,
                cache_dir.map(str::trim).filter(|v| !v.is_empty()),
            )
            .map_err(|e| EngineError {
                code: ErrorCode::InvalidInput,
                message: format!(
                    "Could not inspect prepared compatibility for genome '{}' in catalog '{}': {}",
                    genome_id, catalog_path, e
                ),
            })
    }

    pub fn inspect_helper_genome_prepared_compatibility(
        genome_id: &str,
        catalog_path: Option<&str>,
        cache_dir: Option<&str>,
    ) -> Result<PreparedGenomeCompatibilityInspection, EngineError> {
        let (catalog, catalog_path) = Self::open_helper_genome_catalog(catalog_path)?;
        catalog
            .inspect_prepared_genome_compatibility(
                genome_id,
                cache_dir.map(str::trim).filter(|v| !v.is_empty()),
            )
            .map_err(|e| EngineError {
                code: ErrorCode::InvalidInput,
                message: format!(
                    "Could not inspect prepared compatibility for helper '{}' in catalog '{}': {}",
                    genome_id, catalog_path, e
                ),
            })
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
        let (catalog, _) = Self::open_helper_genome_catalog(catalog_path)?;
        catalog
            .is_prepared(
                genome_id,
                cache_dir.map(str::trim).filter(|v| !v.is_empty()),
            )
            .map_err(|e| EngineError {
                code: ErrorCode::InvalidInput,
                message: format!(
                    "Could not check prepared status for helper '{}': {}",
                    genome_id, e
                ),
            })
    }

    pub fn inspect_reference_genome_prepare_activity(
        catalog_path: Option<&str>,
        genome_id: &str,
        cache_dir: Option<&str>,
    ) -> Result<Option<PrepareGenomeActivityStatus>, EngineError> {
        let (catalog, _) = Self::open_reference_genome_catalog(catalog_path)?;
        catalog
            .inspect_prepare_activity_status(
                genome_id,
                cache_dir.map(str::trim).filter(|v| !v.is_empty()),
            )
            .map_err(|e| EngineError {
                code: ErrorCode::InvalidInput,
                message: format!(
                    "Could not inspect prepare activity for genome '{}': {}",
                    genome_id, e
                ),
            })
    }

    pub fn inspect_helper_genome_prepare_activity(
        genome_id: &str,
        catalog_path: Option<&str>,
        cache_dir: Option<&str>,
    ) -> Result<Option<PrepareGenomeActivityStatus>, EngineError> {
        let (catalog, _) = Self::open_helper_genome_catalog(catalog_path)?;
        catalog
            .inspect_prepare_activity_status(
                genome_id,
                cache_dir.map(str::trim).filter(|v| !v.is_empty()),
            )
            .map_err(|e| EngineError {
                code: ErrorCode::InvalidInput,
                message: format!(
                    "Could not inspect prepare activity for helper '{}': {}",
                    genome_id, e
                ),
            })
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
        let (catalog, catalog_path) = Self::open_helper_genome_catalog(catalog_path)?;
        catalog
            .list_gene_regions(
                genome_id,
                cache_dir.map(str::trim).filter(|v| !v.is_empty()),
            )
            .map_err(|e| EngineError {
                code: ErrorCode::InvalidInput,
                message: format!(
                    "Could not list features for helper '{}' in catalog '{}': {}",
                    genome_id, catalog_path, e
                ),
            })
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
        let (catalog, catalog_path) = Self::open_helper_genome_catalog(catalog_path)?;
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
                    "Could not run BLAST search against helper '{}' in catalog '{}': {}",
                    genome_id, catalog_path, e
                ),
            })
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
        let (catalog, _) = Self::open_helper_genome_catalog(catalog_path)?;
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
        catalog
            .prepare_genome_once_with_progress(
                genome_id,
                cache_dir.map(str::trim).filter(|v| !v.is_empty()),
                &mut guarded_progress,
            )
            .map_err(|e| EngineError {
                code: ErrorCode::Io,
                message: if is_prepare_cancelled_error(&e) {
                    if timed_out {
                        format!(
                            "Genome preparation timed out for '{}' after {} second(s)",
                            genome_id,
                            timeout_seconds.unwrap_or(0)
                        )
                    } else {
                        format!("Genome preparation cancelled for '{genome_id}'")
                    }
                } else {
                    format!("Could not prepare helper genome '{}': {e}", genome_id)
                },
            })
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
                | Operation::RenderTfbsScoreTracksSvg { .. }
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
                | Operation::ExportRackLabelsSvg { .. }
                | Operation::ExportRackFabricationSvg { .. }
                | Operation::ExportRackIsometricSvg { .. }
                | Operation::ExportRackOpenScad { .. }
                | Operation::ExportRackCarrierLabelsSvg { .. }
                | Operation::ExportRackSimulationJson { .. }
                | Operation::ExportDnaLadders { .. }
                | Operation::ExportRnaLadders { .. }
                | Operation::ExportPool { .. }
                | Operation::ExportProcessRunBundle { .. }
                | Operation::ExportFeaturesBed { .. }
                | Operation::InspectSequenceContextView { .. }
                | Operation::ExportSequenceContextBundle { .. }
                | Operation::ListRnaReadReports { .. }
                | Operation::ShowRnaReadReport { .. }
                | Operation::SummarizeRnaReadGeneSupport { .. }
                | Operation::InspectRnaReadGeneSupport { .. }
                | Operation::FindRestrictionSites { .. }
                | Operation::SummarizeTfbsRegion { .. }
                | Operation::SummarizeTfbsScoreTracks { .. }
                | Operation::InspectJasparEntry { .. }
                | Operation::SummarizeJasparEntries { .. }
                | Operation::BenchmarkJasparRegistry { .. }
                | Operation::ListJasparCatalog { .. }
                | Operation::SyncJasparRemoteMetadata { .. }
                | Operation::SummarizeVariantPromoterContext { .. }
                | Operation::SuggestPromoterReporterFragments { .. }
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
                | Operation::SuggestSequencingPrimers { .. }
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

    fn default_extract_genome_promoter_slice_output_id(
        genome_id: &str,
        gene: &GenomeGeneRecord,
        transcript_id: &str,
        upstream_bp: usize,
        downstream_bp: usize,
    ) -> String {
        let label = gene
            .gene_name
            .as_deref()
            .or(gene.gene_id.as_deref())
            .unwrap_or("gene");
        let genome_token = Self::normalize_id_token(genome_id);
        let label_token = Self::normalize_id_token(label);
        let transcript_token = Self::normalize_id_token(transcript_id);
        format!(
            "{}_{}_{}_promoter_{}up_{}down",
            genome_token,
            if transcript_token.is_empty() {
                label_token
            } else {
                format!("{label_token}_{transcript_token}")
            },
            "slice",
            upstream_bp,
            downstream_bp
        )
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

    fn transcript_record_tss_1based(
        selected_gene: &GenomeGeneRecord,
        transcript_record: &GenomeTranscriptRecord,
    ) -> Result<(char, usize), String> {
        let strand = transcript_record
            .strand
            .or(selected_gene.strand)
            .ok_or_else(|| {
                format!(
                    "Transcript '{}' for gene '{}' has no strand annotation",
                    transcript_record.transcript_id,
                    Self::genome_gene_display_label(selected_gene)
                )
            })?;
        let tss_1based = if strand == '-' {
            transcript_record.transcript_end_1based
        } else {
            transcript_record.transcript_start_1based
        };
        Ok((strand, tss_1based))
    }

    fn select_promoter_transcript_record(
        selected_gene: &GenomeGeneRecord,
        transcript_records: &[GenomeTranscriptRecord],
        transcript_id: Option<&str>,
    ) -> Result<GenomeTranscriptRecord, String> {
        if transcript_records.is_empty() {
            return Err(format!(
                "Gene '{}' has no transcript annotation; cannot derive promoter slice from transcript TSS",
                Self::genome_gene_display_label(selected_gene)
            ));
        }
        if let Some(requested) = transcript_id
            .map(str::trim)
            .filter(|value| !value.is_empty())
        {
            let requested_norm = Self::normalize_id_token(requested);
            return transcript_records
                .iter()
                .find(|record| {
                    Self::normalize_id_token(&record.transcript_id) == requested_norm
                        || record
                            .gene_name
                            .as_deref()
                            .map(Self::normalize_id_token)
                            .map(|value| value == requested_norm)
                            .unwrap_or(false)
                })
                .cloned()
                .ok_or_else(|| {
                    format!(
                        "Gene '{}' has no transcript matching '{}'",
                        Self::genome_gene_display_label(selected_gene),
                        requested
                    )
                });
        }

        let mut ranked = transcript_records
            .iter()
            .filter_map(|record| {
                Self::transcript_record_tss_1based(selected_gene, record)
                    .ok()
                    .map(|(strand, tss_1based)| (record.clone(), strand, tss_1based))
            })
            .collect::<Vec<_>>();
        if ranked.is_empty() {
            return Err(format!(
                "Gene '{}' has transcript records but none had usable TSS geometry",
                Self::genome_gene_display_label(selected_gene)
            ));
        }
        ranked.sort_by(
            |(left_record, left_strand, left_tss), (right_record, right_strand, right_tss)| match (
                left_strand,
                right_strand,
            ) {
                ('-', '-') => right_tss
                    .cmp(left_tss)
                    .then(left_record.transcript_id.cmp(&right_record.transcript_id)),
                _ => left_tss
                    .cmp(right_tss)
                    .then(left_record.transcript_id.cmp(&right_record.transcript_id)),
            },
        );
        Ok(ranked.remove(0).0)
    }

    fn resolve_extract_genome_promoter_slice_interval(
        selected_gene: &GenomeGeneRecord,
        transcript_records: &[GenomeTranscriptRecord],
        transcript_id: Option<&str>,
        upstream_bp: usize,
        downstream_bp: usize,
    ) -> Result<(GenomeTranscriptRecord, usize, usize, usize), String> {
        let transcript_record = Self::select_promoter_transcript_record(
            selected_gene,
            transcript_records,
            transcript_id,
        )?;
        let (strand, tss_1based) =
            Self::transcript_record_tss_1based(selected_gene, &transcript_record)?;
        let (start_1based, end_1based) = if strand == '-' {
            (
                tss_1based.saturating_sub(downstream_bp).max(1),
                tss_1based.saturating_add(upstream_bp),
            )
        } else {
            (
                tss_1based.saturating_sub(upstream_bp).max(1),
                tss_1based.saturating_add(downstream_bp),
            )
        };
        Ok((transcript_record, tss_1based, start_1based, end_1based))
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

    fn rebase_name_lookup_by_normalized() -> HashMap<String, String> {
        let mut lookup = HashMap::new();
        for enzyme in active_restriction_enzymes() {
            let normalized_name = Self::normalize_enzyme_match_token(&enzyme.name);
            if !normalized_name.is_empty() {
                lookup.entry(normalized_name).or_insert(enzyme.name);
            }
        }
        lookup
    }

    fn canonicalize_rebase_enzyme_name(
        raw: &str,
        lookup: &HashMap<String, String>,
    ) -> Option<String> {
        let normalized = Self::normalize_enzyme_match_token(raw);
        lookup.get(&normalized).cloned()
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
        let lookup = Self::rebase_name_lookup_by_normalized();
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

    fn restriction_cloning_cut_positions_0based(dna: &DNAsequence) -> BTreeMap<String, Vec<usize>> {
        let seq_len = dna.len();
        let mut by_enzyme: BTreeMap<String, Vec<usize>> = BTreeMap::new();
        for site in dna
            .restriction_enzyme_sites()
            .iter()
            .filter(|site| site.forward_strand)
        {
            let Some((cut_start, _)) = site.recessed_opening_window_0based(seq_len) else {
                continue;
            };
            by_enzyme
                .entry(site.enzyme.name.clone())
                .or_default()
                .push(cut_start);
        }
        by_enzyme
    }

    fn restriction_cloning_mcs_expected_site_order(dna: &DNAsequence) -> Vec<String> {
        let lookup = Self::rebase_name_lookup_by_normalized();
        let mut ordered = vec![];
        let mut seen = HashSet::new();
        for feature in dna
            .features()
            .iter()
            .filter(|feature| Self::feature_looks_like_mcs(feature))
        {
            for raw in feature.qualifier_values("mcs_expected_sites") {
                for token in raw
                    .split(',')
                    .map(str::trim)
                    .filter(|value| !value.is_empty())
                {
                    if let Some(name) = Self::canonicalize_rebase_enzyme_name(token, &lookup) {
                        if seen.insert(name.clone()) {
                            ordered.push(name);
                        }
                    }
                }
            }
        }
        ordered
    }

    pub fn restriction_cloning_vector_enzyme_suggestions(
        &self,
        seq_id: &str,
    ) -> Result<RestrictionCloningVectorEnzymeSuggestions, EngineError> {
        let dna = self
            .state
            .sequences
            .get(seq_id)
            .ok_or_else(|| EngineError {
                code: ErrorCode::NotFound,
                message: format!("Sequence '{}' not found", seq_id),
            })?;
        let cut_positions = Self::restriction_cloning_cut_positions_0based(dna);
        let unique_cutters = cut_positions
            .iter()
            .filter_map(|(name, cuts)| (cuts.len() == 1).then_some(name.clone()))
            .collect::<HashSet<_>>();
        let lookup = Self::rebase_name_lookup_by_normalized();
        let mut selected_mcs = vec![];
        let mut missing_mcs = vec![];
        let mut seen_selected = HashSet::new();
        let mut seen_missing = HashSet::new();
        for feature in dna
            .features()
            .iter()
            .filter(|feature| Self::feature_looks_like_mcs(feature))
        {
            let mut candidates = vec![];
            for raw in feature.qualifier_values("mcs_expected_sites") {
                for token in raw
                    .split(',')
                    .map(str::trim)
                    .filter(|value| !value.is_empty())
                {
                    if let Some(name) = Self::canonicalize_rebase_enzyme_name(token, &lookup) {
                        candidates.push(name);
                    } else {
                        candidates.extend(Self::extract_rebase_enzyme_names_from_text(token));
                    }
                }
            }
            for key in ["note", "label", "gene", "name", "standard_name"] {
                for raw in feature.qualifier_values(key) {
                    candidates.extend(Self::extract_rebase_enzyme_names_from_text(raw));
                }
            }
            for name in candidates {
                if unique_cutters.contains(&name) {
                    if seen_selected.insert(name.clone()) {
                        selected_mcs.push(name);
                    }
                } else if seen_missing.insert(name.clone()) {
                    missing_mcs.push(name);
                }
            }
        }
        let mut other_unique = unique_cutters
            .into_iter()
            .filter(|name| !seen_selected.contains(name))
            .collect::<Vec<_>>();
        other_unique.sort_by(|left, right| {
            cut_positions
                .get(left)
                .and_then(|cuts| cuts.first().copied())
                .unwrap_or(usize::MAX)
                .cmp(
                    &cut_positions
                        .get(right)
                        .and_then(|cuts| cuts.first().copied())
                        .unwrap_or(usize::MAX),
                )
                .then(left.cmp(right))
        });
        let recommended_single_site = selected_mcs
            .iter()
            .chain(other_unique.iter())
            .filter_map(|name| {
                cut_positions
                    .get(name)
                    .and_then(|cuts| cuts.first().copied())
                    .map(
                        |cut_position_0based| RestrictionCloningSingleSiteSuggestion {
                            enzyme: name.clone(),
                            cut_position_0based,
                        },
                    )
            })
            .collect::<Vec<_>>();
        let recommended_directed_pairs = if selected_mcs.len() >= 2 {
            let mut pairs = vec![];
            for left_idx in 0..selected_mcs.len() {
                for right_idx in (left_idx + 1)..selected_mcs.len() {
                    let forward = &selected_mcs[left_idx];
                    let reverse = &selected_mcs[right_idx];
                    let Some(forward_cut_position_0based) = cut_positions
                        .get(forward)
                        .and_then(|cuts| cuts.first().copied())
                    else {
                        continue;
                    };
                    let Some(reverse_cut_position_0based) = cut_positions
                        .get(reverse)
                        .and_then(|cuts| cuts.first().copied())
                    else {
                        continue;
                    };
                    pairs.push(RestrictionCloningDirectedPairSuggestion {
                        order_source: "mcs_expected_sites".to_string(),
                        forward_enzyme: forward.clone(),
                        reverse_enzyme: reverse.clone(),
                        forward_cut_position_0based,
                        reverse_cut_position_0based,
                    });
                }
            }
            pairs
        } else {
            let mut ordered_unique = other_unique
                .iter()
                .filter_map(|name| {
                    cut_positions
                        .get(name)
                        .and_then(|cuts| cuts.first().copied())
                        .map(|cut| (name.clone(), cut))
                })
                .collect::<Vec<_>>();
            ordered_unique.sort_by(|left, right| left.1.cmp(&right.1).then(left.0.cmp(&right.0)));
            let mut pairs = vec![];
            for left_idx in 0..ordered_unique.len() {
                for right_idx in (left_idx + 1)..ordered_unique.len() {
                    let (forward_enzyme, forward_cut_position_0based) = &ordered_unique[left_idx];
                    let (reverse_enzyme, reverse_cut_position_0based) = &ordered_unique[right_idx];
                    pairs.push(RestrictionCloningDirectedPairSuggestion {
                        order_source: "vector_unique_cut_position_order".to_string(),
                        forward_enzyme: forward_enzyme.clone(),
                        reverse_enzyme: reverse_enzyme.clone(),
                        forward_cut_position_0based: *forward_cut_position_0based,
                        reverse_cut_position_0based: *reverse_cut_position_0based,
                    });
                }
            }
            pairs
        };
        Ok(RestrictionCloningVectorEnzymeSuggestions {
            seq_id: seq_id.to_string(),
            selected_mcs,
            other_unique,
            missing_mcs,
            recommended_single_site,
            recommended_directed_pairs,
        })
    }

    fn restriction_cloning_suggested_enzyme_lookup(
        suggestions: &RestrictionCloningVectorEnzymeSuggestions,
    ) -> HashMap<String, String> {
        let mut lookup = HashMap::new();
        for enzyme in suggestions
            .selected_mcs
            .iter()
            .chain(suggestions.other_unique.iter())
        {
            let normalized = Self::normalize_enzyme_match_token(enzyme);
            if !normalized.is_empty() {
                lookup.entry(normalized).or_insert_with(|| enzyme.clone());
            }
        }
        lookup
    }

    fn restriction_cloning_available_enzyme_names(
        suggestions: &RestrictionCloningVectorEnzymeSuggestions,
    ) -> Vec<String> {
        suggestions
            .selected_mcs
            .iter()
            .chain(suggestions.other_unique.iter())
            .cloned()
            .collect::<Vec<_>>()
    }

    pub fn seed_restriction_cloning_pcr_handoff_request(
        &self,
        primer_report_id: &str,
        destination_vector_seq_id: &str,
        pair_rank_1based: Option<usize>,
        mode: RestrictionCloningPcrHandoffMode,
        forward_enzyme: Option<&str>,
        reverse_enzyme: Option<&str>,
        forward_leader_5prime: Option<&str>,
        reverse_leader_5prime: Option<&str>,
    ) -> Result<RestrictionCloningPcrHandoffSeedRequest, EngineError> {
        let report = self.get_primer_design_report(primer_report_id)?;
        if report.pairs.is_empty() {
            return Err(EngineError {
                code: ErrorCode::NotFound,
                message: format!(
                    "Primer-design report '{}' has no accepted primer pairs",
                    report.report_id
                ),
            });
        }
        let pair_rank = pair_rank_1based.unwrap_or(1);
        if pair_rank == 0 {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: "Restriction-cloning pair rank must be >= 1".to_string(),
            });
        }
        let pair_index = pair_rank - 1;
        let selected_pair = report
            .pairs
            .get(pair_index)
            .cloned()
            .ok_or_else(|| EngineError {
                code: ErrorCode::NotFound,
                message: format!(
                    "Primer-design report '{}' does not have accepted pair rank {}",
                    report.report_id, pair_rank
                ),
            })?;
        let pair_geometry = report.pair_core_geometry(&selected_pair);
        let destination_vector_seq_id = destination_vector_seq_id.trim();
        if destination_vector_seq_id.is_empty() {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: "Restriction-cloning destination vector is empty".to_string(),
            });
        }
        let suggestions =
            self.restriction_cloning_vector_enzyme_suggestions(destination_vector_seq_id)?;
        let enzyme_lookup = Self::restriction_cloning_suggested_enzyme_lookup(&suggestions);
        let available_names = Self::restriction_cloning_available_enzyme_names(&suggestions);
        let canonicalize_enzyme = |raw: &str, side_label: &str| -> Result<String, EngineError> {
            let trimmed = raw.trim();
            let normalized = Self::normalize_enzyme_match_token(trimmed);
            enzyme_lookup
                    .get(&normalized)
                    .cloned()
                    .ok_or_else(|| EngineError {
                        code: ErrorCode::InvalidInput,
                        message: if available_names.is_empty() {
                            format!(
                                "Restriction-cloning {} enzyme '{}' is not a unique usable cutter on vector '{}'",
                                side_label, trimmed, destination_vector_seq_id
                            )
                        } else {
                            format!(
                                "Restriction-cloning {} enzyme '{}' is not a unique usable cutter on vector '{}'; available suggestions: {}",
                                side_label,
                                trimmed,
                                destination_vector_seq_id,
                                available_names.join(", ")
                            )
                        },
                    })
        };
        let explicit_forward = forward_enzyme
            .map(str::trim)
            .filter(|value| !value.is_empty());
        let explicit_reverse = reverse_enzyme
            .map(str::trim)
            .filter(|value| !value.is_empty());
        let (resolved_forward, resolved_reverse, selection_source, suggestion_order_source) =
            match mode {
                RestrictionCloningPcrHandoffMode::SingleSite => {
                    let explicit = explicit_forward.or(explicit_reverse);
                    if let (Some(fwd), Some(rev)) = (explicit_forward, explicit_reverse)
                        && !fwd.eq_ignore_ascii_case(rev)
                    {
                        return Err(EngineError {
                            code: ErrorCode::InvalidInput,
                            message: format!(
                                "Restriction-cloning single_site seed requires matching forward/reverse enzymes, got '{}' and '{}'",
                                fwd, rev
                            ),
                        });
                    }
                    if let Some(raw) = explicit {
                        let canonical = canonicalize_enzyme(raw, "single-site")?;
                        (
                            canonical.clone(),
                            canonical,
                            "explicit_single_site".to_string(),
                            None,
                        )
                    } else {
                        let suggestion = suggestions
                            .recommended_single_site
                            .first()
                            .ok_or_else(|| EngineError {
                                code: ErrorCode::NotFound,
                                message: format!(
                                    "Vector '{}' has no recommended single-site restriction cutters for restriction-cloning seed generation",
                                    destination_vector_seq_id
                                ),
                            })?;
                        (
                            suggestion.enzyme.clone(),
                            suggestion.enzyme.clone(),
                            "recommended_single_site".to_string(),
                            None,
                        )
                    }
                }
                RestrictionCloningPcrHandoffMode::DirectedPair => {
                    if explicit_forward.is_some() ^ explicit_reverse.is_some() {
                        return Err(EngineError {
                            code: ErrorCode::InvalidInput,
                            message: "Restriction-cloning directed_pair seed requires both forward and reverse enzymes, or neither to use the top recommendation".to_string(),
                        });
                    }
                    if let (Some(fwd), Some(rev)) = (explicit_forward, explicit_reverse) {
                        let canonical_forward = canonicalize_enzyme(fwd, "forward")?;
                        let canonical_reverse = canonicalize_enzyme(rev, "reverse")?;
                        let matching = suggestions.recommended_directed_pairs.iter().find(|pair| {
                            pair.forward_enzyme.eq_ignore_ascii_case(&canonical_forward)
                                && pair.reverse_enzyme.eq_ignore_ascii_case(&canonical_reverse)
                        });
                        let Some(matching) = matching else {
                            let recommended = suggestions
                                .recommended_directed_pairs
                                .iter()
                                .map(|pair| {
                                    format!("{} -> {}", pair.forward_enzyme, pair.reverse_enzyme)
                                })
                                .collect::<Vec<_>>();
                            return Err(EngineError {
                                code: ErrorCode::InvalidInput,
                                message: if recommended.is_empty() {
                                    format!(
                                        "Restriction-cloning directed pair '{} -> {}' is not available on vector '{}'",
                                        canonical_forward,
                                        canonical_reverse,
                                        destination_vector_seq_id
                                    )
                                } else {
                                    format!(
                                        "Restriction-cloning directed pair '{} -> {}' does not match the valid order on vector '{}'; recommended pairs: {}",
                                        canonical_forward,
                                        canonical_reverse,
                                        destination_vector_seq_id,
                                        recommended.join(", ")
                                    )
                                },
                            });
                        };
                        (
                            canonical_forward,
                            canonical_reverse,
                            "explicit_directed_pair".to_string(),
                            Some(matching.order_source.clone()),
                        )
                    } else {
                        let suggestion = suggestions
                            .recommended_directed_pairs
                            .first()
                            .ok_or_else(|| EngineError {
                                code: ErrorCode::NotFound,
                                message: format!(
                                    "Vector '{}' has no recommended directed restriction-site pairs for restriction-cloning seed generation",
                                    destination_vector_seq_id
                                ),
                            })?;
                        (
                            suggestion.forward_enzyme.clone(),
                            suggestion.reverse_enzyme.clone(),
                            "recommended_directed_pair".to_string(),
                            Some(suggestion.order_source.clone()),
                        )
                    }
                }
            };
        let forward_leader_5prime = forward_leader_5prime
            .map(str::trim)
            .filter(|value| !value.is_empty())
            .unwrap_or_default()
            .to_string();
        let reverse_leader_5prime = reverse_leader_5prime
            .map(str::trim)
            .filter(|value| !value.is_empty())
            .unwrap_or_default()
            .to_string();
        let operation = Operation::PrepareRestrictionCloningPcrHandoff {
            template: report.template.clone(),
            primer_report_id: report.report_id.clone(),
            pair_index,
            destination_vector_seq_id: destination_vector_seq_id.to_string(),
            mode,
            forward_enzyme: resolved_forward.clone(),
            reverse_enzyme: Some(resolved_reverse.clone()),
            forward_leader_5prime: (!forward_leader_5prime.is_empty())
                .then_some(forward_leader_5prime.clone()),
            reverse_leader_5prime: (!reverse_leader_5prime.is_empty())
                .then_some(reverse_leader_5prime.clone()),
        };
        Ok(RestrictionCloningPcrHandoffSeedRequest {
            schema: "gentle.restriction_cloning_pcr_handoff_seed.v1".to_string(),
            primer_report_id: report.report_id.clone(),
            template: report.template.clone(),
            destination_vector_seq_id: destination_vector_seq_id.to_string(),
            pair_index,
            pair_rank,
            selected_pair,
            selected_pair_core_geometry: pair_geometry,
            mode,
            forward_enzyme: resolved_forward,
            reverse_enzyme: resolved_reverse,
            forward_leader_5prime,
            reverse_leader_5prime,
            selection_source,
            suggestion_order_source,
            vector_suggestions: suggestions,
            operation,
        })
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
        if store.reports.is_empty()
            && store.qpcr_reports.is_empty()
            && store.restriction_cloning_handoffs.is_empty()
        {
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

    fn read_protein_derivation_report_store_from_metadata(
        value: Option<&serde_json::Value>,
    ) -> ProteinDerivationReportStore {
        let mut store = value
            .cloned()
            .and_then(|v| serde_json::from_value::<ProteinDerivationReportStore>(v).ok())
            .unwrap_or_default();
        if store.schema.trim().is_empty() {
            store.schema = PROTEIN_DERIVATION_REPORTS_SCHEMA.to_string();
        }
        store
    }

    fn read_protein_derivation_report_store(&self) -> ProteinDerivationReportStore {
        Self::read_protein_derivation_report_store_from_metadata(
            self.state
                .metadata
                .get(PROTEIN_DERIVATION_REPORTS_METADATA_KEY),
        )
    }

    fn write_protein_derivation_report_store(
        &mut self,
        mut store: ProteinDerivationReportStore,
    ) -> Result<(), EngineError> {
        if store.reports.is_empty() {
            self.state
                .metadata
                .remove(PROTEIN_DERIVATION_REPORTS_METADATA_KEY);
            return Ok(());
        }
        store.schema = PROTEIN_DERIVATION_REPORTS_SCHEMA.to_string();
        store.updated_at_unix_ms = Self::now_unix_ms();
        let value = serde_json::to_value(store).map_err(|e| EngineError {
            code: ErrorCode::Internal,
            message: format!("Could not serialize protein-derivation metadata: {e}"),
        })?;
        self.state
            .metadata
            .insert(PROTEIN_DERIVATION_REPORTS_METADATA_KEY.to_string(), value);
        Ok(())
    }

    fn upsert_protein_derivation_report(
        &mut self,
        report: ProteinDerivationReport,
    ) -> Result<(), EngineError> {
        let mut store = self.read_protein_derivation_report_store();
        store.reports.insert(report.report_id.clone(), report);
        self.write_protein_derivation_report_store(store)
    }

    fn normalize_protein_derivation_report_id(raw: &str) -> Result<String, EngineError> {
        let trimmed = raw.trim();
        if trimmed.is_empty() {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: "report_id cannot be empty".to_string(),
            });
        }
        Ok(trimmed.to_string())
    }

    fn read_uniprot_projection_audit_report_store_from_metadata(
        value: Option<&serde_json::Value>,
    ) -> UniprotProjectionAuditReportStore {
        let mut store = value
            .cloned()
            .and_then(|v| serde_json::from_value::<UniprotProjectionAuditReportStore>(v).ok())
            .unwrap_or_default();
        if store.schema.trim().is_empty() {
            store.schema = UNIPROT_PROJECTION_AUDIT_REPORTS_SCHEMA.to_string();
        }
        store
    }

    fn read_uniprot_projection_audit_report_store(&self) -> UniprotProjectionAuditReportStore {
        Self::read_uniprot_projection_audit_report_store_from_metadata(
            self.state
                .metadata
                .get(UNIPROT_PROJECTION_AUDIT_REPORTS_METADATA_KEY),
        )
    }

    fn write_uniprot_projection_audit_report_store(
        &mut self,
        mut store: UniprotProjectionAuditReportStore,
    ) -> Result<(), EngineError> {
        if store.reports.is_empty() {
            self.state
                .metadata
                .remove(UNIPROT_PROJECTION_AUDIT_REPORTS_METADATA_KEY);
            return Ok(());
        }
        store.schema = UNIPROT_PROJECTION_AUDIT_REPORTS_SCHEMA.to_string();
        store.updated_at_unix_ms = Self::now_unix_ms();
        let value = serde_json::to_value(store).map_err(|e| EngineError {
            code: ErrorCode::Internal,
            message: format!("Could not serialize UniProt projection audit metadata: {e}"),
        })?;
        self.state.metadata.insert(
            UNIPROT_PROJECTION_AUDIT_REPORTS_METADATA_KEY.to_string(),
            value,
        );
        Ok(())
    }

    fn read_uniprot_projection_audit_parity_report_store_from_metadata(
        value: Option<&serde_json::Value>,
    ) -> UniprotProjectionAuditParityReportStore {
        let mut store = value
            .cloned()
            .and_then(|v| serde_json::from_value::<UniprotProjectionAuditParityReportStore>(v).ok())
            .unwrap_or_default();
        if store.schema.trim().is_empty() {
            store.schema = UNIPROT_PROJECTION_AUDIT_PARITY_REPORTS_SCHEMA.to_string();
        }
        store
    }

    fn read_uniprot_projection_audit_parity_report_store(
        &self,
    ) -> UniprotProjectionAuditParityReportStore {
        Self::read_uniprot_projection_audit_parity_report_store_from_metadata(
            self.state
                .metadata
                .get(UNIPROT_PROJECTION_AUDIT_PARITY_REPORTS_METADATA_KEY),
        )
    }

    fn write_uniprot_projection_audit_parity_report_store(
        &mut self,
        mut store: UniprotProjectionAuditParityReportStore,
    ) -> Result<(), EngineError> {
        if store.reports.is_empty() {
            self.state
                .metadata
                .remove(UNIPROT_PROJECTION_AUDIT_PARITY_REPORTS_METADATA_KEY);
            return Ok(());
        }
        store.schema = UNIPROT_PROJECTION_AUDIT_PARITY_REPORTS_SCHEMA.to_string();
        store.updated_at_unix_ms = Self::now_unix_ms();
        let value = serde_json::to_value(store).map_err(|e| EngineError {
            code: ErrorCode::Internal,
            message: format!("Could not serialize UniProt projection audit parity metadata: {e}"),
        })?;
        self.state.metadata.insert(
            UNIPROT_PROJECTION_AUDIT_PARITY_REPORTS_METADATA_KEY.to_string(),
            value,
        );
        Ok(())
    }

    fn normalize_uniprot_projection_audit_report_id(raw: &str) -> Result<String, EngineError> {
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

    fn upsert_uniprot_projection_audit_report(
        &mut self,
        report: UniprotProjectionAuditReport,
    ) -> Result<(), EngineError> {
        let mut store = self.read_uniprot_projection_audit_report_store();
        store.reports.insert(report.report_id.clone(), report);
        self.write_uniprot_projection_audit_report_store(store)
    }

    fn upsert_uniprot_projection_audit_parity_report(
        &mut self,
        report: UniprotProjectionAuditParityReport,
    ) -> Result<(), EngineError> {
        let mut store = self.read_uniprot_projection_audit_parity_report_store();
        store.reports.insert(report.report_id.clone(), report);
        self.write_uniprot_projection_audit_parity_report_store(store)
    }

    pub fn list_uniprot_projection_audit_reports(
        &self,
        seq_id_filter: Option<&str>,
    ) -> Vec<UniprotProjectionAuditReportSummary> {
        let filter = seq_id_filter
            .map(str::trim)
            .filter(|value| !value.is_empty())
            .map(|value| value.to_ascii_lowercase());
        let mut rows = self
            .read_uniprot_projection_audit_report_store()
            .reports
            .values()
            .filter(|report| {
                filter
                    .as_ref()
                    .is_none_or(|needle| report.seq_id.to_ascii_lowercase() == *needle)
            })
            .map(|report| UniprotProjectionAuditReportSummary {
                report_id: report.report_id.clone(),
                projection_id: report.projection_id.clone(),
                entry_id: report.entry_id.clone(),
                seq_id: report.seq_id.clone(),
                ensembl_entry_id: report.ensembl_entry_id.clone(),
                generated_at_unix_ms: report.generated_at_unix_ms,
                op_id: report.op_id.clone(),
                run_id: report.run_id.clone(),
                transcript_count: report.rows.len(),
                failing_transcript_count: report
                    .rows
                    .iter()
                    .filter(|row| row.status != UniprotProjectionAuditRowStatus::Consistent)
                    .count(),
            })
            .collect::<Vec<_>>();
        rows.sort_by(|left, right| {
            left.seq_id
                .to_ascii_lowercase()
                .cmp(&right.seq_id.to_ascii_lowercase())
                .then(left.generated_at_unix_ms.cmp(&right.generated_at_unix_ms))
                .then(
                    left.report_id
                        .to_ascii_lowercase()
                        .cmp(&right.report_id.to_ascii_lowercase()),
                )
        });
        rows
    }

    pub fn get_uniprot_projection_audit_report(
        &self,
        report_id: &str,
    ) -> Result<UniprotProjectionAuditReport, EngineError> {
        let report_id = Self::normalize_uniprot_projection_audit_report_id(report_id)?;
        self.read_uniprot_projection_audit_report_store()
            .reports
            .get(&report_id)
            .cloned()
            .ok_or_else(|| EngineError {
                code: ErrorCode::NotFound,
                message: format!("UniProt projection audit report '{}' not found", report_id),
            })
    }

    pub fn export_uniprot_projection_audit_report(
        &self,
        report_id: &str,
        path: &str,
    ) -> Result<UniprotProjectionAuditReport, EngineError> {
        let report = self.get_uniprot_projection_audit_report(report_id)?;
        let text = serde_json::to_string_pretty(&report).map_err(|e| EngineError {
            code: ErrorCode::Internal,
            message: format!(
                "Could not serialize UniProt projection audit report '{}': {e}",
                report.report_id
            ),
        })?;
        std::fs::write(path, text).map_err(|e| EngineError {
            code: ErrorCode::Io,
            message: format!("Could not write UniProt projection audit report to '{path}': {e}"),
        })?;
        Ok(report)
    }

    pub fn list_uniprot_projection_audit_parity_reports(
        &self,
        seq_id_filter: Option<&str>,
    ) -> Vec<UniprotProjectionAuditParityReportSummary> {
        let filter = seq_id_filter
            .map(str::trim)
            .filter(|value| !value.is_empty())
            .map(|value| value.to_ascii_lowercase());
        let mut rows = self
            .read_uniprot_projection_audit_parity_report_store()
            .reports
            .values()
            .filter(|report| {
                filter
                    .as_ref()
                    .is_none_or(|needle| report.seq_id.to_ascii_lowercase() == *needle)
            })
            .map(|report| UniprotProjectionAuditParityReportSummary {
                report_id: report.report_id.clone(),
                projection_id: report.projection_id.clone(),
                entry_id: report.entry_id.clone(),
                seq_id: report.seq_id.clone(),
                ensembl_entry_id: report.ensembl_entry_id.clone(),
                generated_at_unix_ms: report.generated_at_unix_ms,
                op_id: report.op_id.clone(),
                run_id: report.run_id.clone(),
                transcript_count: report.rows.len(),
                divergent_transcript_count: report
                    .rows
                    .iter()
                    .filter(|row| {
                        !row.statuses_match
                            || !row.accounting_match
                            || !row.mismatch_reason_match
                            || !row.comparison_mode_match
                    })
                    .count(),
            })
            .collect::<Vec<_>>();
        rows.sort_by(|left, right| {
            left.seq_id
                .to_ascii_lowercase()
                .cmp(&right.seq_id.to_ascii_lowercase())
                .then(left.generated_at_unix_ms.cmp(&right.generated_at_unix_ms))
                .then(
                    left.report_id
                        .to_ascii_lowercase()
                        .cmp(&right.report_id.to_ascii_lowercase()),
                )
        });
        rows
    }

    pub fn get_uniprot_projection_audit_parity_report(
        &self,
        report_id: &str,
    ) -> Result<UniprotProjectionAuditParityReport, EngineError> {
        let report_id = Self::normalize_uniprot_projection_audit_report_id(report_id)?;
        self.read_uniprot_projection_audit_parity_report_store()
            .reports
            .get(&report_id)
            .cloned()
            .ok_or_else(|| EngineError {
                code: ErrorCode::NotFound,
                message: format!(
                    "UniProt projection audit parity report '{}' not found",
                    report_id
                ),
            })
    }

    pub fn export_uniprot_projection_audit_parity_report(
        &self,
        report_id: &str,
        path: &str,
    ) -> Result<UniprotProjectionAuditParityReport, EngineError> {
        let report = self.get_uniprot_projection_audit_parity_report(report_id)?;
        let text = serde_json::to_string_pretty(&report).map_err(|e| EngineError {
            code: ErrorCode::Internal,
            message: format!(
                "Could not serialize UniProt projection audit parity report '{}': {e}",
                report.report_id
            ),
        })?;
        std::fs::write(path, text).map_err(|e| EngineError {
            code: ErrorCode::Io,
            message: format!(
                "Could not write UniProt projection audit parity report to '{path}': {e}"
            ),
        })?;
        Ok(report)
    }

    fn summarize_protein_derivation_modes(rows: &[ProteinDerivationReportRow]) -> String {
        let mut labels = rows
            .iter()
            .map(|row| row.derivation.derivation_mode.as_str().to_string())
            .collect::<BTreeSet<_>>();
        match labels.len() {
            0 => "-".to_string(),
            1 => labels.pop_first().unwrap_or_else(|| "-".to_string()),
            _ => "mixed".to_string(),
        }
    }

    fn read_reverse_translation_report_store_from_metadata(
        value: Option<&serde_json::Value>,
    ) -> ReverseTranslationReportStore {
        let mut store = value
            .cloned()
            .and_then(|v| serde_json::from_value::<ReverseTranslationReportStore>(v).ok())
            .unwrap_or_default();
        if store.schema.trim().is_empty() {
            store.schema = REVERSE_TRANSLATION_REPORTS_SCHEMA.to_string();
        }
        store
    }

    fn read_reverse_translation_report_store(&self) -> ReverseTranslationReportStore {
        Self::read_reverse_translation_report_store_from_metadata(
            self.state
                .metadata
                .get(REVERSE_TRANSLATION_REPORTS_METADATA_KEY),
        )
    }

    fn write_reverse_translation_report_store(
        &mut self,
        mut store: ReverseTranslationReportStore,
    ) -> Result<(), EngineError> {
        if store.reports.is_empty() {
            self.state
                .metadata
                .remove(REVERSE_TRANSLATION_REPORTS_METADATA_KEY);
            return Ok(());
        }
        store.schema = REVERSE_TRANSLATION_REPORTS_SCHEMA.to_string();
        store.updated_at_unix_ms = Self::now_unix_ms();
        let value = serde_json::to_value(store).map_err(|e| EngineError {
            code: ErrorCode::Internal,
            message: format!("Could not serialize reverse-translation metadata: {e}"),
        })?;
        self.state
            .metadata
            .insert(REVERSE_TRANSLATION_REPORTS_METADATA_KEY.to_string(), value);
        Ok(())
    }

    fn upsert_reverse_translation_report(
        &mut self,
        report: ReverseTranslationReport,
    ) -> Result<(), EngineError> {
        let mut store = self.read_reverse_translation_report_store();
        store.reports.insert(report.report_id.clone(), report);
        self.write_reverse_translation_report_store(store)
    }

    fn summarize_reverse_translation_speed_profile(report: &ReverseTranslationReport) -> String {
        let mut summary = report
            .resolved_speed_profile
            .map(|profile| profile.as_str().to_string())
            .unwrap_or_else(|| "auto".to_string());
        if let Some(mark) = report.speed_mark {
            summary.push(':');
            summary.push_str(mark.as_str());
        }
        summary
    }

    fn summarize_reverse_translation_diagnostics(report: &ReverseTranslationReport) -> String {
        let mut parts = Vec::new();
        let synonym_total =
            report.preferred_synonymous_choice_count + report.alternative_synonymous_choice_count;
        if synonym_total > 0 {
            parts.push(format!(
                "preferred={}/{}",
                report.preferred_synonymous_choice_count, synonym_total
            ));
            if report.alternative_synonymous_choice_count > 0 {
                parts.push(format!(
                    "alt={}",
                    report.alternative_synonymous_choice_count
                ));
            }
        }
        if report.fallback_unknown_codon_count > 0 {
            parts.push(format!("fallback={}", report.fallback_unknown_codon_count));
        }
        if let Some(gc_fraction) = report.gc_fraction {
            parts.push(format!("gc={:.0}%", gc_fraction * 100.0));
        }
        if let Some(realized_tm) = report.realized_anneal_tm_c {
            parts.push(format!("tm={realized_tm:.1}C"));
        }
        if parts.is_empty() {
            "-".to_string()
        } else {
            parts.join(", ")
        }
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
                op_id: report.op_id.clone(),
                run_id: report.run_id.clone(),
                roi_start_0based: report.roi_start_0based,
                roi_end_0based: report.roi_end_0based,
                pair_count: report.pair_count,
                backend_used: report.backend.used.clone(),
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
                op_id: report.op_id.clone(),
                run_id: report.run_id.clone(),
                roi_start_0based: report.roi_start_0based,
                roi_end_0based: report.roi_end_0based,
                assay_count: report.assay_count,
                best_assay_probe_placement: report.best_assay_probe_placement.clone(),
                best_assay_summary: report.best_assay_summary.clone(),
                backend_used: report.backend.used.clone(),
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

    pub fn list_restriction_cloning_pcr_handoffs(
        &self,
    ) -> Vec<RestrictionCloningPcrHandoffReportSummary> {
        let store = self.read_primer_design_store();
        let mut ids = store
            .restriction_cloning_handoffs
            .keys()
            .cloned()
            .collect::<Vec<_>>();
        ids.sort();
        ids.into_iter()
            .filter_map(|id| store.restriction_cloning_handoffs.get(&id))
            .map(|report| RestrictionCloningPcrHandoffReportSummary {
                report_id: report.report_id.clone(),
                template: report.template.clone(),
                primer_report_id: report.primer_report_id.clone(),
                pair_index: report.pair_index,
                pair_rank: report.pair_rank,
                destination_vector_seq_id: report.destination_vector_seq_id.clone(),
                generated_at_unix_ms: report.generated_at_unix_ms,
                op_id: report.op_id.clone(),
                run_id: report.run_id.clone(),
                mode: report.mode.as_str().to_string(),
                forward_enzyme: report.forward_enzyme.clone(),
                reverse_enzyme: report.reverse_enzyme.clone(),
                compatibility_status: report.compatibility.status.clone(),
            })
            .collect()
    }

    pub fn get_restriction_cloning_pcr_handoff(
        &self,
        report_id: &str,
    ) -> Result<RestrictionCloningPcrHandoffReport, EngineError> {
        let report_id = Self::normalize_primer_design_report_id(report_id)?;
        let store = self.read_primer_design_store();
        store
            .restriction_cloning_handoffs
            .get(&report_id)
            .cloned()
            .ok_or_else(|| EngineError {
                code: ErrorCode::NotFound,
                message: format!(
                    "Restriction-cloning PCR handoff report '{}' not found",
                    report_id
                ),
            })
    }

    pub fn export_restriction_cloning_pcr_handoff(
        &self,
        report_id: &str,
        path: &str,
    ) -> Result<RestrictionCloningPcrHandoffReport, EngineError> {
        let report = self.get_restriction_cloning_pcr_handoff(report_id)?;
        let text = serde_json::to_string_pretty(&report).map_err(|e| EngineError {
            code: ErrorCode::Internal,
            message: format!(
                "Could not serialize restriction-cloning PCR handoff report '{}': {e}",
                report.report_id
            ),
        })?;
        std::fs::write(path, text).map_err(|e| EngineError {
            code: ErrorCode::Io,
            message: format!(
                "Could not write restriction-cloning PCR handoff report to '{path}': {e}"
            ),
        })?;
        Ok(report)
    }

    pub fn list_protein_derivation_reports(
        &self,
        seq_id_filter: Option<&str>,
    ) -> Vec<ProteinDerivationReportSummary> {
        let filter = seq_id_filter
            .map(str::trim)
            .filter(|value| !value.is_empty())
            .map(|value| value.to_ascii_lowercase());
        let mut rows = self
            .read_protein_derivation_report_store()
            .reports
            .values()
            .filter(|report| {
                filter
                    .as_ref()
                    .is_none_or(|needle| report.seq_id.to_ascii_lowercase() == *needle)
            })
            .map(|report| ProteinDerivationReportSummary {
                report_id: report.report_id.clone(),
                seq_id: report.seq_id.clone(),
                generated_at_unix_ms: report.generated_at_unix_ms,
                op_id: report.op_id.clone(),
                run_id: report.run_id.clone(),
                effective_output_prefix: report.effective_output_prefix.clone(),
                derived_count: report.derived_count,
                derivation_mode_summary: Self::summarize_protein_derivation_modes(&report.rows),
            })
            .collect::<Vec<_>>();
        rows.sort_by(|left, right| {
            left.seq_id
                .to_ascii_lowercase()
                .cmp(&right.seq_id.to_ascii_lowercase())
                .then(left.generated_at_unix_ms.cmp(&right.generated_at_unix_ms))
                .then(
                    left.report_id
                        .to_ascii_lowercase()
                        .cmp(&right.report_id.to_ascii_lowercase()),
                )
        });
        rows
    }

    pub fn get_protein_derivation_report(
        &self,
        report_id: &str,
    ) -> Result<ProteinDerivationReport, EngineError> {
        let report_id = Self::normalize_protein_derivation_report_id(report_id)?;
        self.read_protein_derivation_report_store()
            .reports
            .get(&report_id)
            .cloned()
            .ok_or_else(|| EngineError {
                code: ErrorCode::NotFound,
                message: format!("Protein-derivation report '{}' not found", report_id),
            })
    }

    pub fn export_protein_derivation_report(
        &self,
        report_id: &str,
        path: &str,
    ) -> Result<ProteinDerivationReport, EngineError> {
        let report = self.get_protein_derivation_report(report_id)?;
        let text = serde_json::to_string_pretty(&report).map_err(|e| EngineError {
            code: ErrorCode::Internal,
            message: format!(
                "Could not serialize protein-derivation report '{}': {e}",
                report.report_id
            ),
        })?;
        std::fs::write(path, text).map_err(|e| EngineError {
            code: ErrorCode::Io,
            message: format!("Could not write protein-derivation report to '{path}': {e}"),
        })?;
        Ok(report)
    }

    pub fn list_reverse_translation_reports(
        &self,
        protein_seq_id_filter: Option<&str>,
    ) -> Vec<ReverseTranslationReportSummary> {
        let filter = protein_seq_id_filter
            .map(str::trim)
            .filter(|value| !value.is_empty())
            .map(|value| value.to_ascii_lowercase());
        let mut rows = self
            .read_reverse_translation_report_store()
            .reports
            .values()
            .filter(|report| {
                filter
                    .as_ref()
                    .is_none_or(|needle| report.protein_seq_id.to_ascii_lowercase() == *needle)
            })
            .map(|report| ReverseTranslationReportSummary {
                report_id: report.report_id.clone(),
                protein_seq_id: report.protein_seq_id.clone(),
                coding_seq_id: report.coding_seq_id.clone(),
                generated_at_unix_ms: report.generated_at_unix_ms,
                op_id: report.op_id.clone(),
                run_id: report.run_id.clone(),
                protein_length_aa: report.protein_length_aa,
                coding_length_bp: report.coding_length_bp,
                translation_table: report.translation_table,
                speed_profile_summary: Self::summarize_reverse_translation_speed_profile(report),
                diagnostics_summary: Self::summarize_reverse_translation_diagnostics(report),
            })
            .collect::<Vec<_>>();
        rows.sort_by(|left, right| {
            left.protein_seq_id
                .to_ascii_lowercase()
                .cmp(&right.protein_seq_id.to_ascii_lowercase())
                .then(left.generated_at_unix_ms.cmp(&right.generated_at_unix_ms))
                .then(
                    left.report_id
                        .to_ascii_lowercase()
                        .cmp(&right.report_id.to_ascii_lowercase()),
                )
        });
        rows
    }

    pub fn get_reverse_translation_report(
        &self,
        report_id: &str,
    ) -> Result<ReverseTranslationReport, EngineError> {
        let report_id = Self::normalize_protein_derivation_report_id(report_id)?;
        self.read_reverse_translation_report_store()
            .reports
            .get(&report_id)
            .cloned()
            .ok_or_else(|| EngineError {
                code: ErrorCode::NotFound,
                message: format!("Reverse-translation report '{}' not found", report_id),
            })
    }

    pub fn export_reverse_translation_report(
        &self,
        report_id: &str,
        path: &str,
    ) -> Result<ReverseTranslationReport, EngineError> {
        let report = self.get_reverse_translation_report(report_id)?;
        let text = serde_json::to_string_pretty(&report).map_err(|e| EngineError {
            code: ErrorCode::Internal,
            message: format!(
                "Could not serialize reverse-translation report '{}': {e}",
                report.report_id
            ),
        })?;
        std::fs::write(path, text).map_err(|e| EngineError {
            code: ErrorCode::Io,
            message: format!("Could not write reverse-translation report to '{path}': {e}"),
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

    fn normalize_sequence_context_visible_class_token(raw: &str) -> Option<String> {
        let normalized = raw.trim().to_ascii_lowercase();
        match normalized.as_str() {
            "" => None,
            "gene" | "genes" => Some("gene".to_string()),
            "mrna" | "transcript" | "transcripts" => Some("mrna".to_string()),
            "exon" | "exons" => Some("exon".to_string()),
            "cds" => Some("cds".to_string()),
            "promoter" | "promoters" => Some("promoter".to_string()),
            "variation" | "variant" | "variants" | "snp" | "snps" => Some("variation".to_string()),
            "tfbs" | "motif" | "motifs" => Some("tfbs".to_string()),
            _ => None,
        }
    }

    fn sequence_context_feature_kinds_for_class(class_id: &str) -> &'static [&'static str] {
        match class_id {
            "gene" => &["GENE"],
            "mrna" => &["MRNA"],
            "exon" => &["EXON"],
            "cds" => &["CDS"],
            "promoter" => &["PROMOTER"],
            "variation" => &["VARIATION"],
            "tfbs" => &["TFBS"],
            _ => &[],
        }
    }

    fn default_sequence_context_visible_classes(&self) -> Vec<String> {
        let display = &self.state.display;
        let mut classes: Vec<String> = vec![];
        if display.show_features {
            if display.show_gene_features {
                classes.push("gene".to_string());
            }
            if display.show_mrna_features {
                classes.push("mrna".to_string());
            }
            classes.push("exon".to_string());
            if display.show_cds_features {
                classes.push("cds".to_string());
            }
            classes.push("promoter".to_string());
            classes.push("variation".to_string());
            if display.show_tfbs {
                classes.push("tfbs".to_string());
            }
        }
        classes
    }

    fn genomic_interval_from_anchor(
        anchor: &SequenceGenomeAnchorSummary,
        start_0based: usize,
        end_0based_exclusive: usize,
    ) -> Option<(usize, usize)> {
        if end_0based_exclusive <= start_0based {
            return None;
        }
        let span = end_0based_exclusive - start_0based;
        if span == 0 {
            return None;
        }
        match anchor.strand.unwrap_or('+') {
            '-' => {
                let genomic_end = anchor.end_1based.checked_sub(start_0based)?;
                let genomic_start = anchor.end_1based.checked_sub(end_0based_exclusive - 1)?;
                Some((
                    genomic_start.min(genomic_end),
                    genomic_start.max(genomic_end),
                ))
            }
            _ => {
                let genomic_start = anchor.start_1based.checked_add(start_0based)?;
                let genomic_end = anchor.start_1based.checked_add(end_0based_exclusive - 1)?;
                Some((
                    genomic_start.min(genomic_end),
                    genomic_start.max(genomic_end),
                ))
            }
        }
    }

    fn sequence_context_row_genomic_coordinates(
        row: &SequenceFeatureQueryRow,
        genome_anchor: Option<&SequenceGenomeAnchorSummary>,
    ) -> (Option<String>, Option<usize>, Option<usize>) {
        let chromosome = row
            .qualifiers
            .get("chromosome")
            .and_then(|values| values.first())
            .map(|value| value.trim().to_string())
            .filter(|value| !value.is_empty());
        let qualifier_start = row
            .qualifiers
            .get("genomic_start_1based")
            .and_then(|values| values.first())
            .and_then(|value| value.trim().parse::<usize>().ok());
        let qualifier_end = row
            .qualifiers
            .get("genomic_end_1based")
            .and_then(|values| values.first())
            .and_then(|value| value.trim().parse::<usize>().ok());
        if chromosome.is_some() && qualifier_start.is_some() && qualifier_end.is_some() {
            return (chromosome, qualifier_start, qualifier_end);
        }
        let Some(anchor) = genome_anchor else {
            return (None, None, None);
        };
        let Some((genomic_start, genomic_end)) =
            Self::genomic_interval_from_anchor(anchor, row.start_0based, row.end_0based_exclusive)
        else {
            return (None, None, None);
        };
        (
            Some(anchor.chromosome.clone()),
            Some(genomic_start),
            Some(genomic_end),
        )
    }

    pub fn inspect_sequence_context_view(
        &self,
        seq_id: &str,
        mode: Option<RenderSvgMode>,
        viewport_start_0based: Option<usize>,
        viewport_end_0based_exclusive: Option<usize>,
        include_visible_classes: &[String],
        coordinate_mode: Option<FeatureBedCoordinateMode>,
        limit: Option<usize>,
    ) -> Result<SequenceContextViewReport, EngineError> {
        let seq_id = seq_id.trim();
        if seq_id.is_empty() {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: "inspect_sequence_context_view requires non-empty seq_id".to_string(),
            });
        }
        let dna = self
            .state
            .sequences
            .get(seq_id)
            .ok_or_else(|| EngineError {
                code: ErrorCode::NotFound,
                message: format!("Sequence '{}' not found", seq_id),
            })?;
        let sequence_length_bp = dna.len();
        let mode = mode.unwrap_or(RenderSvgMode::Linear);
        if viewport_start_0based.is_some() != viewport_end_0based_exclusive.is_some() {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: "Sequence-context inspection requires both viewport_start_0based and viewport_end_0based_exclusive"
                    .to_string(),
            });
        }
        let (viewport_start_0based, viewport_end_0based_exclusive) = if let (
            Some(start),
            Some(end),
        ) =
            (viewport_start_0based, viewport_end_0based_exclusive)
        {
            if end <= start {
                return Err(EngineError {
                    code: ErrorCode::InvalidInput,
                    message: format!(
                        "Invalid sequence-context viewport: start ({start}) must be < end ({end})"
                    ),
                });
            }
            if end > sequence_length_bp {
                return Err(EngineError {
                    code: ErrorCode::InvalidInput,
                    message: format!(
                        "Sequence-context viewport {}..{} is outside sequence length {}",
                        start, end, sequence_length_bp
                    ),
                });
            }
            (start, end)
        } else if matches!(mode, RenderSvgMode::Linear)
            && self.state.display.linear_view_span_bp > 0
            && sequence_length_bp > 0
        {
            let max_start = sequence_length_bp.saturating_sub(1);
            let start = self.state.display.linear_view_start_bp.min(max_start);
            let span = self
                .state
                .display
                .linear_view_span_bp
                .clamp(1, sequence_length_bp.saturating_sub(start).max(1));
            (start, start + span)
        } else {
            (0, sequence_length_bp)
        };
        let viewport_span_bp = viewport_end_0based_exclusive.saturating_sub(viewport_start_0based);
        let coordinate_mode = coordinate_mode.unwrap_or(FeatureBedCoordinateMode::Auto);
        let coordinate_mode_label = match coordinate_mode {
            FeatureBedCoordinateMode::Auto => "auto",
            FeatureBedCoordinateMode::Local => "local",
            FeatureBedCoordinateMode::Genomic => "genomic",
        };
        let mode_label = match mode {
            RenderSvgMode::Linear => "linear",
            RenderSvgMode::Circular => "circular",
        };

        let visible_classes = if include_visible_classes.is_empty() {
            self.default_sequence_context_visible_classes()
        } else {
            include_visible_classes
                .iter()
                .filter_map(|value| Self::normalize_sequence_context_visible_class_token(value))
                .collect::<BTreeSet<_>>()
                .into_iter()
                .collect::<Vec<_>>()
        };

        let genome_anchor = self.sequence_genome_anchor_summary(seq_id).ok();
        let (genomic_view_start_1based, genomic_view_end_1based) = genome_anchor
            .as_ref()
            .and_then(|anchor| {
                Self::genomic_interval_from_anchor(
                    anchor,
                    viewport_start_0based,
                    viewport_end_0based_exclusive,
                )
            })
            .map(|(start, end)| (Some(start), Some(end)))
            .unwrap_or((None, None));

        let limit = limit
            .unwrap_or(FEATURE_QUERY_DEFAULT_LIMIT.min(25))
            .clamp(1, FEATURE_QUERY_MAX_LIMIT);
        let union_kinds = visible_classes
            .iter()
            .flat_map(|class_id| Self::sequence_context_feature_kinds_for_class(class_id))
            .map(|kind| kind.to_string())
            .collect::<BTreeSet<_>>()
            .into_iter()
            .collect::<Vec<_>>();

        let mut report = SequenceContextViewReport {
            schema: SEQUENCE_CONTEXT_VIEW_SCHEMA.to_string(),
            seq_id: seq_id.to_string(),
            sequence_length_bp,
            generated_at_unix_ms: Self::now_unix_ms(),
            op_id: None,
            run_id: None,
            mode: mode_label.to_string(),
            coordinate_mode: coordinate_mode_label.to_string(),
            viewport_start_0based,
            viewport_end_0based_exclusive,
            viewport_span_bp,
            genome_anchor,
            genomic_view_start_1based,
            genomic_view_end_1based,
            visible_classes: vec![],
            matched_feature_count: 0,
            returned_feature_count: 0,
            limit,
            rows: vec![],
            summary_lines: vec![],
        };

        for class_id in &visible_classes {
            let class_kinds = Self::sequence_context_feature_kinds_for_class(class_id);
            let class_result = self.query_sequence_features(SequenceFeatureQuery {
                seq_id: seq_id.to_string(),
                kind_in: class_kinds.iter().map(|kind| kind.to_string()).collect(),
                start_0based: (viewport_span_bp > 0).then_some(viewport_start_0based),
                end_0based_exclusive: (viewport_span_bp > 0)
                    .then_some(viewport_end_0based_exclusive),
                range_relation: SequenceFeatureRangeRelation::Overlap,
                sort_by: SequenceFeatureSortBy::Start,
                include_qualifiers: true,
                limit: Some(limit.min(5)),
                ..SequenceFeatureQuery::default()
            })?;
            let prominent_labels = class_result
                .rows
                .iter()
                .map(|row| row.label.trim().to_string())
                .filter(|label| !label.is_empty())
                .collect::<BTreeSet<_>>()
                .into_iter()
                .take(3)
                .collect::<Vec<_>>();
            report.visible_classes.push(SequenceContextVisibleClass {
                class_id: class_id.clone(),
                feature_kinds: class_kinds.iter().map(|kind| kind.to_string()).collect(),
                matched_count: class_result.matched_count,
                returned_count: class_result.returned_count,
                prominent_labels,
            });
        }

        if !union_kinds.is_empty() {
            let query_result = self.query_sequence_features(SequenceFeatureQuery {
                seq_id: seq_id.to_string(),
                kind_in: union_kinds,
                start_0based: (viewport_span_bp > 0).then_some(viewport_start_0based),
                end_0based_exclusive: (viewport_span_bp > 0)
                    .then_some(viewport_end_0based_exclusive),
                range_relation: SequenceFeatureRangeRelation::Overlap,
                sort_by: SequenceFeatureSortBy::Start,
                include_qualifiers: true,
                limit: Some(limit),
                ..SequenceFeatureQuery::default()
            })?;
            report.matched_feature_count = query_result.matched_count;
            report.returned_feature_count = query_result.returned_count;
            report.rows = query_result
                .rows
                .into_iter()
                .map(|row| {
                    let (chromosome, genomic_start_1based, genomic_end_1based) =
                        Self::sequence_context_row_genomic_coordinates(
                            &row,
                            report.genome_anchor.as_ref(),
                        );
                    SequenceContextFeatureRow {
                        feature_id: row.feature_id,
                        kind: row.kind,
                        start_0based: row.start_0based,
                        end_0based_exclusive: row.end_0based_exclusive,
                        length_bp: row.length_bp,
                        strand: row.strand,
                        label: row.label,
                        labels: row.labels,
                        chromosome,
                        genomic_start_1based,
                        genomic_end_1based,
                    }
                })
                .collect();
        }

        report.summary_lines.push(format!(
            "{} view {}..{} ({} bp) of {} bp",
            report.mode,
            report.viewport_start_0based + 1,
            report.viewport_end_0based_exclusive,
            report.viewport_span_bp,
            report.sequence_length_bp
        ));
        if let (Some(anchor), Some(genomic_start), Some(genomic_end)) = (
            report.genome_anchor.as_ref(),
            report.genomic_view_start_1based,
            report.genomic_view_end_1based,
        ) {
            report.summary_lines.push(format!(
                "genome anchor: {}:{}-{} ({}, strand {}, {})",
                anchor.chromosome,
                genomic_start,
                genomic_end,
                anchor.genome_id,
                anchor.strand.unwrap_or('+'),
                match anchor.anchor_verified {
                    Some(true) => "verified",
                    Some(false) => "unverified",
                    None => "verification n/a",
                }
            ));
        }
        if !report.visible_classes.is_empty() {
            report.summary_lines.push(format!(
                "visible classes: {}",
                report
                    .visible_classes
                    .iter()
                    .map(|class| format!("{}({})", class.class_id, class.matched_count))
                    .collect::<Vec<_>>()
                    .join(", ")
            ));
        } else {
            report
                .summary_lines
                .push("visible classes: none currently enabled".to_string());
        }
        let top_labels = report
            .rows
            .iter()
            .map(|row| row.label.trim().to_string())
            .filter(|label| !label.is_empty())
            .collect::<BTreeSet<_>>()
            .into_iter()
            .take(5)
            .collect::<Vec<_>>();
        if !top_labels.is_empty() {
            report
                .summary_lines
                .push(format!("top labels: {}", top_labels.join(", ")));
        }

        Ok(report)
    }

    fn write_pretty_json_file<T: Serialize>(
        &self,
        value: &T,
        path: &str,
        artifact_label: &str,
    ) -> Result<(), EngineError> {
        let text = serde_json::to_string_pretty(value).map_err(|e| EngineError {
            code: ErrorCode::Internal,
            message: format!("Could not serialize {artifact_label} for '{path}': {e}"),
        })?;
        std::fs::write(path, text).map_err(|e| EngineError {
            code: ErrorCode::Io,
            message: format!("Could not write {artifact_label} to '{path}': {e}"),
        })
    }

    fn write_text_file(
        &self,
        path: &str,
        contents: &str,
        artifact_label: &str,
    ) -> Result<(), EngineError> {
        std::fs::write(path, contents).map_err(|e| EngineError {
            code: ErrorCode::Io,
            message: format!("Could not write {artifact_label} to '{path}': {e}"),
        })
    }

    fn sequence_feature_bed_columns() -> Vec<String> {
        [
            "chrom",
            "chrom_start_0based",
            "chrom_end_0based_exclusive",
            "name",
            "score_0_to_1000",
            "strand",
            "kind",
            "row_id",
            "coordinate_source",
            "qualifiers_json",
        ]
        .into_iter()
        .map(str::to_string)
        .collect()
    }

    fn empty_sequence_feature_bed_export_report(
        seq_id: &str,
        path: &str,
        coordinate_mode: FeatureBedCoordinateMode,
        include_restriction_sites: bool,
        restriction_enzymes: Vec<String>,
        query: SequenceFeatureQuery,
    ) -> SequenceFeatureBedExportReport {
        SequenceFeatureBedExportReport {
            schema: FEATURE_BED_EXPORT_REPORT_SCHEMA.to_string(),
            seq_id: seq_id.to_string(),
            path: path.to_string(),
            coordinate_mode: coordinate_mode.as_str().to_string(),
            include_restriction_sites,
            restriction_enzyme_filters: restriction_enzymes,
            bed_columns: Self::sequence_feature_bed_columns(),
            matched_sequence_feature_count: 0,
            matched_restriction_site_count: 0,
            matched_row_count: 0,
            exportable_row_count: 0,
            exported_row_count: 0,
            offset: 0,
            limit: query.limit,
            local_coordinate_row_count: 0,
            genomic_coordinate_row_count: 0,
            skipped_missing_genomic_coordinates: 0,
            query,
        }
    }

    fn sequence_context_artifact_caption(
        artifact_kind: &str,
        seq_id: &str,
        sequence_context_view: &SequenceContextViewReport,
    ) -> String {
        let span_label = format!(
            "{}..{} ({} bp)",
            sequence_context_view.viewport_start_0based + 1,
            sequence_context_view.viewport_end_0based_exclusive,
            sequence_context_view.viewport_span_bp
        );
        match artifact_kind {
            "context_svg" => format!(
                "Sequence-context {} figure for '{}' covering {}",
                sequence_context_view.mode, seq_id, span_label
            ),
            "context_summary_json" => format!(
                "Machine-readable sequence-context summary for '{}' covering {}",
                seq_id, span_label
            ),
            "context_summary_text" => format!(
                "Compact text summary for '{}' covering {}",
                seq_id, span_label
            ),
            "context_feature_bed" => format!(
                "Coordinate-bearing feature table for '{}' covering {}",
                seq_id, span_label
            ),
            "bundle_manifest" => {
                format!("Bundle manifest for '{}' sequence-context export", seq_id)
            }
            _ => format!("Sequence-context artifact '{artifact_kind}' for '{seq_id}'"),
        }
    }

    fn sequence_context_bundle_artifact_manifest(
        seq_id: &str,
        sequence_context_view: &SequenceContextViewReport,
        svg_path: &Path,
        summary_json_path: &Path,
        summary_text_path: Option<&Path>,
        feature_bed_path: Option<&Path>,
        bundle_json_path: &Path,
    ) -> Vec<SequenceContextBundleArtifact> {
        let mut artifacts = vec![SequenceContextBundleArtifact {
            artifact_id: "context_svg".to_string(),
            path: svg_path.to_string_lossy().into_owned(),
            media_type: "image/svg+xml".to_string(),
            artifact_kind: "figure".to_string(),
            caption: Self::sequence_context_artifact_caption(
                "context_svg",
                seq_id,
                sequence_context_view,
            ),
            recommended_use: "best_first_figure".to_string(),
            presentation_rank: 0,
            is_best_first_artifact: true,
        }];
        artifacts.push(SequenceContextBundleArtifact {
            artifact_id: "context_summary_json".to_string(),
            path: summary_json_path.to_string_lossy().into_owned(),
            media_type: "application/json".to_string(),
            artifact_kind: "summary_json".to_string(),
            caption: Self::sequence_context_artifact_caption(
                "context_summary_json",
                seq_id,
                sequence_context_view,
            ),
            recommended_use: "machine_readable_context".to_string(),
            presentation_rank: 1,
            is_best_first_artifact: false,
        });
        if let Some(summary_text_path) = summary_text_path {
            artifacts.push(SequenceContextBundleArtifact {
                artifact_id: "context_summary_text".to_string(),
                path: summary_text_path.to_string_lossy().into_owned(),
                media_type: "text/plain".to_string(),
                artifact_kind: "summary_text".to_string(),
                caption: Self::sequence_context_artifact_caption(
                    "context_summary_text",
                    seq_id,
                    sequence_context_view,
                ),
                recommended_use: "compact_chat_summary".to_string(),
                presentation_rank: 2,
                is_best_first_artifact: false,
            });
        }
        if let Some(feature_bed_path) = feature_bed_path {
            artifacts.push(SequenceContextBundleArtifact {
                artifact_id: "context_feature_bed".to_string(),
                path: feature_bed_path.to_string_lossy().into_owned(),
                media_type: "text/plain".to_string(),
                artifact_kind: "feature_bed".to_string(),
                caption: Self::sequence_context_artifact_caption(
                    "context_feature_bed",
                    seq_id,
                    sequence_context_view,
                ),
                recommended_use: "coordinate_table".to_string(),
                presentation_rank: 3,
                is_best_first_artifact: false,
            });
        }
        artifacts.push(SequenceContextBundleArtifact {
            artifact_id: "bundle_manifest".to_string(),
            path: bundle_json_path.to_string_lossy().into_owned(),
            media_type: "application/json".to_string(),
            artifact_kind: "bundle_manifest".to_string(),
            caption: Self::sequence_context_artifact_caption(
                "bundle_manifest",
                seq_id,
                sequence_context_view,
            ),
            recommended_use: "replay_manifest".to_string(),
            presentation_rank: 4,
            is_best_first_artifact: false,
        });
        artifacts
    }

    pub fn export_sequence_context_bundle(
        &self,
        seq_id: &str,
        mode: Option<RenderSvgMode>,
        viewport_start_0based: Option<usize>,
        viewport_end_0based_exclusive: Option<usize>,
        coordinate_mode: Option<FeatureBedCoordinateMode>,
        include_feature_bed: Option<bool>,
        include_text_summary: Option<bool>,
        include_restriction_sites: Option<bool>,
        restriction_enzymes: &[String],
        output_dir: &str,
        op_id: Option<&str>,
        run_id: Option<&str>,
    ) -> Result<SequenceContextBundleExport, EngineError> {
        let output_dir = output_dir.trim();
        if output_dir.is_empty() {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: "export_sequence_context_bundle requires non-empty output_dir".to_string(),
            });
        }
        let seq_id = seq_id.trim();
        if seq_id.is_empty() {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: "export_sequence_context_bundle requires non-empty seq_id".to_string(),
            });
        }

        let mode = mode.unwrap_or(RenderSvgMode::Linear);
        let coordinate_mode = coordinate_mode.unwrap_or(FeatureBedCoordinateMode::Auto);
        let include_feature_bed = include_feature_bed.unwrap_or(true);
        let include_text_summary = include_text_summary.unwrap_or(true);
        let include_restriction_sites = include_restriction_sites.unwrap_or(false);
        let restriction_enzymes = restriction_enzymes
            .iter()
            .map(|value| value.trim())
            .filter(|value| !value.is_empty())
            .map(|value| value.to_string())
            .collect::<BTreeSet<_>>()
            .into_iter()
            .collect::<Vec<_>>();

        let mut sequence_context_view = self.inspect_sequence_context_view(
            seq_id,
            Some(mode.clone()),
            viewport_start_0based,
            viewport_end_0based_exclusive,
            &[],
            Some(coordinate_mode),
            None,
        )?;
        sequence_context_view.op_id = op_id.map(str::to_string);
        sequence_context_view.run_id = run_id.map(str::to_string);

        let output_dir_path = PathBuf::from(output_dir);
        std::fs::create_dir_all(&output_dir_path).map_err(|e| EngineError {
            code: ErrorCode::Io,
            message: format!(
                "Could not create sequence-context bundle directory '{}': {e}",
                output_dir
            ),
        })?;

        let svg_path = output_dir_path.join("context.svg");
        let summary_json_path = output_dir_path.join("context_summary.json");
        let summary_text_path =
            include_text_summary.then(|| output_dir_path.join("context_summary.txt"));
        let feature_bed_path =
            include_feature_bed.then(|| output_dir_path.join("context_features.bed"));
        let bundle_json_path = output_dir_path.join("bundle.json");

        let dna = self
            .state
            .sequences
            .get(seq_id)
            .ok_or_else(|| EngineError {
                code: ErrorCode::NotFound,
                message: format!("Sequence '{seq_id}' not found"),
            })?;
        let mut display = self.state.display.clone();
        if matches!(mode, RenderSvgMode::Linear) && sequence_context_view.sequence_length_bp > 0 {
            display.linear_view_start_bp = sequence_context_view.viewport_start_0based;
            display.linear_view_span_bp = sequence_context_view.viewport_span_bp.max(1);
        }
        let svg = match mode {
            RenderSvgMode::Linear => export_linear_svg(dna, &display),
            RenderSvgMode::Circular => export_circular_svg(dna, &display),
        };
        self.write_text_file(
            svg_path.to_string_lossy().as_ref(),
            &svg,
            "sequence-context SVG",
        )?;
        self.write_pretty_json_file(
            &sequence_context_view,
            summary_json_path.to_string_lossy().as_ref(),
            "sequence-context summary JSON",
        )?;
        if let Some(summary_text_path) = summary_text_path.as_ref() {
            let summary_text = if sequence_context_view.summary_lines.is_empty() {
                String::new()
            } else {
                format!("{}\n", sequence_context_view.summary_lines.join("\n"))
            };
            self.write_text_file(
                summary_text_path.to_string_lossy().as_ref(),
                &summary_text,
                "sequence-context summary text",
            )?;
        }

        let feature_bed_export = if let Some(feature_bed_path) = feature_bed_path.as_ref() {
            let visible_kinds = sequence_context_view
                .visible_classes
                .iter()
                .flat_map(|class| Self::sequence_context_feature_kinds_for_class(&class.class_id))
                .copied()
                .map(str::to_string)
                .collect::<BTreeSet<_>>()
                .into_iter()
                .collect::<Vec<_>>();
            let query = SequenceFeatureQuery {
                seq_id: seq_id.to_string(),
                include_qualifiers: true,
                kind_in: visible_kinds.clone(),
                start_0based: Some(sequence_context_view.viewport_start_0based),
                end_0based_exclusive: Some(sequence_context_view.viewport_end_0based_exclusive),
                range_relation: SequenceFeatureRangeRelation::Overlap,
                sort_by: SequenceFeatureSortBy::Start,
                ..SequenceFeatureQuery::default()
            };
            if visible_kinds.is_empty() && !include_restriction_sites {
                File::create(feature_bed_path).map_err(|e| EngineError {
                    code: ErrorCode::Io,
                    message: format!(
                        "Could not create empty sequence-context BED '{}': {e}",
                        feature_bed_path.display()
                    ),
                })?;
                Some(Self::empty_sequence_feature_bed_export_report(
                    seq_id,
                    feature_bed_path.to_string_lossy().as_ref(),
                    coordinate_mode,
                    include_restriction_sites,
                    restriction_enzymes.clone(),
                    query,
                ))
            } else {
                Some(self.export_sequence_features_bed(
                    query,
                    feature_bed_path.to_string_lossy().as_ref(),
                    coordinate_mode,
                    include_restriction_sites,
                    &restriction_enzymes,
                )?)
            }
        } else {
            None
        };

        let artifacts = Self::sequence_context_bundle_artifact_manifest(
            seq_id,
            &sequence_context_view,
            &svg_path,
            &summary_json_path,
            summary_text_path.as_deref(),
            feature_bed_path.as_deref(),
            &bundle_json_path,
        );
        let best_first_artifact = artifacts
            .iter()
            .find(|artifact| artifact.is_best_first_artifact)
            .or_else(|| artifacts.first());
        let best_first_artifact_id =
            best_first_artifact.map(|artifact| artifact.artifact_id.clone());
        let best_first_artifact_path = best_first_artifact.map(|artifact| artifact.path.clone());
        let bundle = SequenceContextBundleExport {
            schema: SEQUENCE_CONTEXT_BUNDLE_SCHEMA.to_string(),
            seq_id: seq_id.to_string(),
            generated_at_unix_ms: Self::now_unix_ms(),
            op_id: op_id.map(str::to_string),
            run_id: run_id.map(str::to_string),
            output_dir: output_dir.to_string(),
            svg_path: svg_path.to_string_lossy().into_owned(),
            summary_json_path: summary_json_path.to_string_lossy().into_owned(),
            summary_text_path: summary_text_path
                .as_ref()
                .map(|path| path.to_string_lossy().into_owned()),
            feature_bed_path: feature_bed_path
                .as_ref()
                .map(|path| path.to_string_lossy().into_owned()),
            bundle_json_path: bundle_json_path.to_string_lossy().into_owned(),
            include_text_summary,
            include_feature_bed,
            include_restriction_sites,
            restriction_enzymes,
            artifacts,
            best_first_artifact_id,
            best_first_artifact_path,
            sequence_context_view,
            feature_bed_export,
        };
        self.write_pretty_json_file(
            &bundle,
            bundle_json_path.to_string_lossy().as_ref(),
            "sequence-context bundle manifest",
        )?;
        Ok(bundle)
    }

    fn restriction_end_geometry_label(enzyme: &RestrictionEnzyme) -> &'static str {
        match enzyme.end_geometry() {
            crate::restriction_enzyme::RestrictionEndGeometry::Blunt => "blunt",
            crate::restriction_enzyme::RestrictionEndGeometry::FivePrimeOverhang(_) => {
                "5prime_overhang"
            }
            crate::restriction_enzyme::RestrictionEndGeometry::ThreePrimeOverhang(_) => {
                "3prime_overhang"
            }
        }
    }

    fn validate_sequence_scan_span(
        sequence_length_bp: usize,
        span_start_0based: Option<usize>,
        span_end_0based_exclusive: Option<usize>,
        context: &str,
    ) -> Result<(usize, usize), EngineError> {
        if sequence_length_bp == 0 {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: format!("{context} cannot be computed on an empty sequence"),
            });
        }
        if span_start_0based.is_some() != span_end_0based_exclusive.is_some() {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: format!(
                    "{context} span requires both span_start_0based and span_end_0based_exclusive"
                ),
            });
        }
        let (start, end_exclusive) = match (span_start_0based, span_end_0based_exclusive) {
            (Some(start), Some(end_exclusive)) => (start, end_exclusive),
            (None, None) => (0, sequence_length_bp),
            _ => unreachable!("span presence already validated"),
        };
        if end_exclusive <= start {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: format!(
                    "Invalid {context} span {start}..{end_exclusive}: end must be > start"
                ),
            });
        }
        if start >= sequence_length_bp || end_exclusive > sequence_length_bp {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: format!(
                    "{context} span {start}..{end_exclusive} is outside sequence length {sequence_length_bp}"
                ),
            });
        }
        Ok((start, end_exclusive))
    }

    fn resolve_sequence_scan_target(
        &self,
        target: &SequenceScanTarget,
    ) -> Result<
        (
            String,
            String,
            usize,
            usize,
            usize,
            InlineSequenceTopology,
            DNAsequence,
        ),
        EngineError,
    > {
        match target {
            SequenceScanTarget::SeqId {
                seq_id,
                span_start_0based,
                span_end_0based_exclusive,
            } => {
                let seq_id = seq_id.trim();
                if seq_id.is_empty() {
                    return Err(EngineError {
                        code: ErrorCode::InvalidInput,
                        message: "FindRestrictionSites requires non-empty seq_id".to_string(),
                    });
                }
                let dna = self
                    .state
                    .sequences
                    .get(seq_id)
                    .ok_or_else(|| EngineError {
                        code: ErrorCode::NotFound,
                        message: format!("Sequence '{seq_id}' not found"),
                    })?;
                let source_sequence_length_bp = dna.len();
                let (scan_start_0based, scan_end_0based_exclusive) =
                    Self::validate_sequence_scan_span(
                        source_sequence_length_bp,
                        *span_start_0based,
                        *span_end_0based_exclusive,
                        "FindRestrictionSites",
                    )?;
                let full_span = scan_start_0based == 0
                    && scan_end_0based_exclusive == source_sequence_length_bp;
                let scan_topology = if full_span && dna.is_circular() {
                    InlineSequenceTopology::Circular
                } else {
                    InlineSequenceTopology::Linear
                };
                let scan_dna = if full_span {
                    dna.clone()
                } else {
                    dna.extract_region_preserving_features(
                        scan_start_0based,
                        scan_end_0based_exclusive,
                    )
                    .ok_or_else(|| EngineError {
                        code: ErrorCode::InvalidInput,
                        message: format!(
                            "Could not extract scan span {}..{} from '{}'",
                            scan_start_0based, scan_end_0based_exclusive, seq_id
                        ),
                    })?
                };
                Ok((
                    "seq_id".to_string(),
                    seq_id.to_string(),
                    source_sequence_length_bp,
                    scan_start_0based,
                    scan_end_0based_exclusive,
                    scan_topology,
                    scan_dna,
                ))
            }
            SequenceScanTarget::InlineSequence {
                sequence_text,
                topology,
                id_hint,
                span_start_0based,
                span_end_0based_exclusive,
            } => {
                let mut dna =
                    DNAsequence::from_sequence(sequence_text).map_err(|e| EngineError {
                        code: ErrorCode::InvalidInput,
                        message: format!("Could not parse inline DNA sequence: {e}"),
                    })?;
                dna.set_circular(matches!(topology, InlineSequenceTopology::Circular));
                let source_sequence_length_bp = dna.len();
                let (scan_start_0based, scan_end_0based_exclusive) =
                    Self::validate_sequence_scan_span(
                        source_sequence_length_bp,
                        *span_start_0based,
                        *span_end_0based_exclusive,
                        "FindRestrictionSites",
                    )?;
                let full_span = scan_start_0based == 0
                    && scan_end_0based_exclusive == source_sequence_length_bp;
                let scan_topology = if full_span {
                    *topology
                } else {
                    InlineSequenceTopology::Linear
                };
                let scan_dna = if full_span {
                    dna
                } else {
                    dna.extract_region_preserving_features(
                        scan_start_0based,
                        scan_end_0based_exclusive,
                    )
                    .ok_or_else(|| EngineError {
                        code: ErrorCode::InvalidInput,
                        message: format!(
                            "Could not extract inline scan span {}..{}",
                            scan_start_0based, scan_end_0based_exclusive
                        ),
                    })?
                };
                let target_label = id_hint
                    .as_deref()
                    .map(str::trim)
                    .filter(|value| !value.is_empty())
                    .unwrap_or("inline_sequence")
                    .to_string();
                Ok((
                    "inline_sequence".to_string(),
                    target_label,
                    source_sequence_length_bp,
                    scan_start_0based,
                    scan_end_0based_exclusive,
                    scan_topology,
                    scan_dna,
                ))
            }
        }
    }

    fn resolve_restriction_scan_enzymes(
        &self,
        enzymes: &[String],
    ) -> Result<(Vec<String>, Vec<RestrictionEnzyme>), EngineError> {
        let effective_requested = if enzymes.is_empty() {
            let preferred = self
                .state
                .display
                .preferred_restriction_enzymes
                .iter()
                .map(|value| value.trim())
                .filter(|value| !value.is_empty())
                .map(|value| value.to_string())
                .collect::<Vec<_>>();
            if preferred.is_empty() {
                default_preferred_restriction_enzyme_names()
            } else {
                preferred
            }
        } else {
            enzymes
                .iter()
                .map(|value| value.trim())
                .filter(|value| !value.is_empty())
                .map(|value| value.to_string())
                .collect::<Vec<_>>()
        };

        if effective_requested.is_empty() {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message:
                    "FindRestrictionSites requires at least one enzyme or configured preferred restriction enzymes"
                        .to_string(),
            });
        }

        let catalog = active_restriction_enzymes();
        let mut seen: HashSet<String> = HashSet::new();
        let mut filters: Vec<String> = vec![];
        let mut resolved: Vec<RestrictionEnzyme> = vec![];
        let mut missing: Vec<String> = vec![];
        for raw in effective_requested {
            let normalized = raw.to_ascii_uppercase();
            if !seen.insert(normalized) {
                continue;
            }
            if let Some(enzyme) = catalog
                .iter()
                .find(|enzyme| enzyme.name.eq_ignore_ascii_case(&raw))
            {
                filters.push(enzyme.name.clone());
                resolved.push(enzyme.clone());
            } else {
                missing.push(raw);
            }
        }
        if !missing.is_empty() {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: format!(
                    "Unknown restriction enzyme(s) for FindRestrictionSites: {}",
                    missing.join(", ")
                ),
            });
        }
        Ok((filters, resolved))
    }

    pub fn find_restriction_sites(
        &self,
        target: SequenceScanTarget,
        enzymes: &[String],
        max_sites_per_enzyme: Option<usize>,
        include_cut_geometry: bool,
        op_id: Option<&str>,
        run_id: Option<&str>,
    ) -> Result<RestrictionSiteScanReport, EngineError> {
        if max_sites_per_enzyme.is_some_and(|value| value == 0) {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: "FindRestrictionSites max_sites_per_enzyme must be >= 1".to_string(),
            });
        }
        let (
            target_kind,
            target_label,
            source_sequence_length_bp,
            scan_start_0based,
            scan_end_0based_exclusive,
            scan_topology,
            scan_dna,
        ) = self.resolve_sequence_scan_target(&target)?;
        let (enzyme_filters, resolved_enzymes) = self.resolve_restriction_scan_enzymes(enzymes)?;

        let mut report = RestrictionSiteScanReport {
            schema: RESTRICTION_SITE_SCAN_REPORT_SCHEMA.to_string(),
            target_kind,
            target_label,
            source_sequence_length_bp,
            scan_start_0based,
            scan_end_0based_exclusive,
            scan_length_bp: scan_dna.len(),
            scan_topology,
            generated_at_unix_ms: Self::now_unix_ms(),
            op_id: op_id.map(str::to_string),
            run_id: run_id.map(str::to_string),
            enzyme_filters: enzyme_filters.clone(),
            enzymes_scanned: enzyme_filters,
            max_sites_per_enzyme,
            include_cut_geometry,
            matched_site_count: 0,
            skipped_enzyme_names_due_to_max_sites: vec![],
            path: None,
            rows: vec![],
        };

        for enzyme in resolved_enzymes {
            let sites = enzyme.get_sites(&scan_dna, None);
            if let Some(max) = max_sites_per_enzyme
                && sites.len() > max
            {
                report
                    .skipped_enzyme_names_due_to_max_sites
                    .push(enzyme.name.clone());
                continue;
            }
            for site in sites {
                let Some((recognition_start_0based, recognition_end_0based_exclusive)) =
                    site.recognition_bounds_0based(scan_dna.len())
                else {
                    continue;
                };
                let (
                    forward_cut_0based,
                    reverse_cut_0based,
                    opening_start_0based,
                    opening_end_0based_exclusive,
                ) = if include_cut_geometry {
                    let (forward_cut_0based, reverse_cut_0based) = site
                        .strand_cut_positions_0based(scan_dna.len())
                        .unwrap_or((recognition_start_0based, recognition_start_0based));
                    let (opening_start_0based, opening_end_0based_exclusive) = site
                        .recessed_opening_window_0based(scan_dna.len())
                        .unwrap_or((forward_cut_0based, reverse_cut_0based));
                    (
                        Some(forward_cut_0based),
                        Some(reverse_cut_0based),
                        Some(opening_start_0based),
                        Some(opening_end_0based_exclusive),
                    )
                } else {
                    (None, None, None, None)
                };
                report.rows.push(RestrictionSiteScanHit {
                    enzyme_name: site.enzyme.name.clone(),
                    recognition_sequence: site.enzyme.sequence.clone(),
                    recognition_start_0based,
                    recognition_end_0based_exclusive,
                    source_recognition_start_0based: scan_start_0based
                        .saturating_add(recognition_start_0based),
                    source_recognition_end_0based_exclusive: scan_start_0based
                        .saturating_add(recognition_end_0based_exclusive),
                    recognition_length_bp: recognition_end_0based_exclusive
                        .saturating_sub(recognition_start_0based),
                    forward_strand: site.forward_strand,
                    end_geometry: Self::restriction_end_geometry_label(&site.enzyme).to_string(),
                    note: site
                        .enzyme
                        .note
                        .as_deref()
                        .map(str::trim)
                        .filter(|value| !value.is_empty())
                        .map(str::to_string),
                    forward_cut_0based,
                    reverse_cut_0based,
                    opening_start_0based,
                    opening_end_0based_exclusive,
                    source_forward_cut_0based: forward_cut_0based
                        .map(|value| scan_start_0based.saturating_add(value)),
                    source_reverse_cut_0based: reverse_cut_0based
                        .map(|value| scan_start_0based.saturating_add(value)),
                    source_opening_start_0based: opening_start_0based
                        .map(|value| scan_start_0based.saturating_add(value)),
                    source_opening_end_0based_exclusive: opening_end_0based_exclusive
                        .map(|value| scan_start_0based.saturating_add(value)),
                });
            }
        }

        report.skipped_enzyme_names_due_to_max_sites.sort();
        report.rows.sort_by(|a, b| {
            a.recognition_start_0based
                .cmp(&b.recognition_start_0based)
                .then(a.enzyme_name.cmp(&b.enzyme_name))
                .then(
                    a.recognition_end_0based_exclusive
                        .cmp(&b.recognition_end_0based_exclusive),
                )
                .then(
                    a.source_recognition_start_0based
                        .cmp(&b.source_recognition_start_0based),
                )
        });
        report.matched_site_count = report.rows.len();
        Ok(report)
    }

    fn tfbs_summary_group_name(feature: &gb_io::seq::Feature, feature_id: usize) -> String {
        Self::first_nonempty_feature_qualifier(
            feature,
            &["bound_moiety", "standard_name", "gene", "name"],
        )
        .or_else(|| Self::feature_qualifier_text(feature, "tf_id"))
        .unwrap_or_else(|| Self::feature_display_label(feature, feature_id))
    }

    pub fn summarize_tfbs_region(
        &self,
        mut request: TfbsRegionSummaryRequest,
    ) -> Result<TfbsRegionSummary, EngineError> {
        #[derive(Default)]
        struct TfbsSummaryAccumulator {
            tf_name: String,
            motif_ids: BTreeSet<String>,
            focus_occurrences: usize,
            context_occurrences: usize,
        }

        request.seq_id = request.seq_id.trim().to_string();
        if request.seq_id.is_empty() {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: "summarize_tfbs_region requires non-empty seq_id".to_string(),
            });
        }

        let dna = self
            .state
            .sequences
            .get(&request.seq_id)
            .ok_or_else(|| EngineError {
                code: ErrorCode::NotFound,
                message: format!("Sequence '{}' not found", request.seq_id),
            })?;
        let sequence_length_bp = dna.len();
        if sequence_length_bp == 0 {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: "TFBS region summary cannot be computed on an empty sequence".to_string(),
            });
        }

        let focus_start = request.focus_start_0based;
        let focus_end = request.focus_end_0based_exclusive;
        if focus_end <= focus_start {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: format!(
                    "Invalid TFBS focus range {}..{}: end must be > start",
                    focus_start, focus_end
                ),
            });
        }
        if focus_start >= sequence_length_bp || focus_end > sequence_length_bp {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: format!(
                    "TFBS focus range {}..{} is outside sequence length {}",
                    focus_start, focus_end, sequence_length_bp
                ),
            });
        }

        if request.context_start_0based.is_some() != request.context_end_0based_exclusive.is_some()
        {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message:
                    "TFBS context range requires both context_start_0based and context_end_0based_exclusive"
                        .to_string(),
            });
        }
        let (context_start, context_end) = match (
            request.context_start_0based,
            request.context_end_0based_exclusive,
        ) {
            (Some(start), Some(end)) => (start, end),
            (None, None) => (0, sequence_length_bp),
            _ => unreachable!("context range presence already validated"),
        };
        if context_end <= context_start {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: format!(
                    "Invalid TFBS context range {}..{}: end must be > start",
                    context_start, context_end
                ),
            });
        }
        if context_start >= sequence_length_bp || context_end > sequence_length_bp {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: format!(
                    "TFBS context range {}..{} is outside sequence length {}",
                    context_start, context_end, sequence_length_bp
                ),
            });
        }
        if context_start > focus_start || context_end < focus_end {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: format!(
                    "TFBS context range {}..{} must contain focus range {}..{}",
                    context_start, context_end, focus_start, focus_end
                ),
            });
        }

        if request.limit.is_some_and(|limit| limit == 0) {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: "TFBS region summary limit must be >= 1".to_string(),
            });
        }
        let limit = request
            .limit
            .unwrap_or(TFBS_REGION_SUMMARY_DEFAULT_LIMIT)
            .clamp(1, TFBS_REGION_SUMMARY_MAX_LIMIT);
        request.limit = Some(limit);
        request.context_start_0based = Some(context_start);
        request.context_end_0based_exclusive = Some(context_end);

        let focus_width_bp = focus_end.saturating_sub(focus_start);
        let context_width_bp = context_end.saturating_sub(context_start);
        let outside_focus_width_bp = context_width_bp.saturating_sub(focus_width_bp);
        let total_tfbs_feature_count = dna
            .features()
            .iter()
            .filter(|feature| Self::is_tfbs_feature(feature))
            .count();

        let mut grouped: BTreeMap<String, TfbsSummaryAccumulator> = BTreeMap::new();
        let mut focus_hit_count = 0usize;
        let mut context_hit_count = 0usize;
        for (feature_id, feature) in dna.features().iter().enumerate() {
            if !Self::is_tfbs_feature(feature) {
                continue;
            }
            let overlaps_context = Self::feature_overlaps_span(feature, context_start, context_end);
            if !overlaps_context {
                continue;
            }
            let overlaps_focus = Self::feature_overlaps_span(feature, focus_start, focus_end);
            context_hit_count += 1;
            if overlaps_focus {
                focus_hit_count += 1;
            }

            let tf_name = Self::tfbs_summary_group_name(feature, feature_id);
            let group_key = tf_name.to_ascii_uppercase();
            let entry = grouped.entry(group_key).or_default();
            if entry.tf_name.is_empty() {
                entry.tf_name = tf_name;
            }
            entry.context_occurrences += 1;
            if overlaps_focus {
                entry.focus_occurrences += 1;
            }
            for tf_id in feature.qualifier_values("tf_id") {
                let trimmed = tf_id.trim();
                if !trimmed.is_empty() {
                    entry.motif_ids.insert(trimmed.to_string());
                }
            }
        }

        let mut rows = grouped
            .into_values()
            .filter(|row| {
                row.focus_occurrences >= request.min_focus_occurrences
                    && row.context_occurrences >= request.min_context_occurrences
            })
            .map(|row| {
                let outside_focus_occurrences = row
                    .context_occurrences
                    .saturating_sub(row.focus_occurrences);
                let focus_density_per_kb =
                    row.focus_occurrences as f64 * 1000.0 / focus_width_bp as f64;
                let context_density_per_kb =
                    row.context_occurrences as f64 * 1000.0 / context_width_bp as f64;
                let outside_focus_density_per_kb = if outside_focus_width_bp > 0 {
                    outside_focus_occurrences as f64 * 1000.0 / outside_focus_width_bp as f64
                } else {
                    0.0
                };
                let focus_share_of_context_occurrences = if row.context_occurrences > 0 {
                    row.focus_occurrences as f64 / row.context_occurrences as f64
                } else {
                    0.0
                };
                TfbsRegionSummaryRow {
                    tf_name: row.tf_name,
                    motif_ids: row.motif_ids.into_iter().collect(),
                    focus_occurrences: row.focus_occurrences,
                    context_occurrences: row.context_occurrences,
                    outside_focus_occurrences,
                    focus_density_per_kb,
                    context_density_per_kb,
                    outside_focus_density_per_kb,
                    focus_share_of_context_occurrences,
                    focus_vs_context_density_ratio: (context_density_per_kb > 0.0)
                        .then_some(focus_density_per_kb / context_density_per_kb),
                    focus_vs_outside_density_ratio: (outside_focus_width_bp > 0
                        && outside_focus_density_per_kb > 0.0)
                        .then_some(focus_density_per_kb / outside_focus_density_per_kb),
                }
            })
            .collect::<Vec<_>>();

        rows.sort_by(|a, b| {
            b.focus_occurrences
                .cmp(&a.focus_occurrences)
                .then(b.context_occurrences.cmp(&a.context_occurrences))
                .then(
                    b.outside_focus_occurrences
                        .cmp(&a.outside_focus_occurrences),
                )
                .then(
                    a.tf_name
                        .to_ascii_uppercase()
                        .cmp(&b.tf_name.to_ascii_uppercase()),
                )
        });

        let matched_tf_count = rows.len();
        rows.truncate(limit);
        let returned_tf_count = rows.len();

        Ok(TfbsRegionSummary {
            schema: TFBS_REGION_SUMMARY_SCHEMA.to_string(),
            seq_id: request.seq_id,
            sequence_length_bp,
            focus_start_0based: focus_start,
            focus_end_0based_exclusive: focus_end,
            context_start_0based: context_start,
            context_end_0based_exclusive: context_end,
            focus_width_bp,
            context_width_bp,
            outside_focus_width_bp,
            total_tfbs_feature_count,
            focus_hit_count,
            context_hit_count,
            matched_tf_count,
            returned_tf_count,
            min_focus_occurrences: request.min_focus_occurrences,
            min_context_occurrences: request.min_context_occurrences,
            limit,
            rows,
        })
    }

    pub(crate) fn write_tfbs_region_summary_json(
        &self,
        summary: &TfbsRegionSummary,
        path: &str,
    ) -> Result<(), EngineError> {
        let text = serde_json::to_string_pretty(summary).map_err(|e| EngineError {
            code: ErrorCode::Internal,
            message: format!(
                "Could not serialize TFBS region summary '{}' for '{}': {e}",
                summary.seq_id, path
            ),
        })?;
        std::fs::write(path, text).map_err(|e| EngineError {
            code: ErrorCode::Io,
            message: format!("Could not write TFBS region summary to '{path}': {e}"),
        })
    }

    pub fn export_sequence_features_bed(
        &self,
        mut query: SequenceFeatureQuery,
        path: &str,
        coordinate_mode: FeatureBedCoordinateMode,
        include_restriction_sites: bool,
        restriction_enzymes: &[String],
    ) -> Result<SequenceFeatureBedExportReport, EngineError> {
        fn bed_safe_field(raw: &str) -> String {
            let sanitized = raw.replace(['\t', '\n', '\r'], " ");
            let trimmed = sanitized.trim();
            if trimmed.is_empty() {
                ".".to_string()
            } else {
                trimmed.to_string()
            }
        }

        fn collect_qualifiers(feature: &gb_io::seq::Feature) -> BTreeMap<String, Vec<String>> {
            let mut qualifiers: BTreeMap<String, Vec<String>> = BTreeMap::new();
            for (key, value) in &feature.qualifiers {
                let Some(value) = value.as_deref().map(str::trim).filter(|v| !v.is_empty()) else {
                    continue;
                };
                qualifiers
                    .entry(key.to_string())
                    .or_default()
                    .push(value.to_string());
            }
            qualifiers
        }

        fn compile_qualifier_filters(
            query: &SequenceFeatureQuery,
        ) -> Result<Vec<(SequenceFeatureQualifierFilter, Option<Regex>)>, EngineError> {
            let mut qualifier_filters: Vec<(SequenceFeatureQualifierFilter, Option<Regex>)> =
                vec![];
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
            Ok(qualifier_filters)
        }

        fn qualifiers_match(
            qualifiers: &BTreeMap<String, Vec<String>>,
            qualifier_filters: &[(SequenceFeatureQualifierFilter, Option<Regex>)],
        ) -> bool {
            qualifier_filters.iter().all(|(filter, regex)| {
                let values = qualifiers
                    .get(filter.key.as_str())
                    .cloned()
                    .unwrap_or_default()
                    .into_iter()
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
                if let Some(regex) = regex
                    && !values.iter().any(|value| regex.is_match(value))
                {
                    return false;
                }
                true
            })
        }

        fn label_matches(
            labels: &[String],
            label_contains_upper: Option<&String>,
            label_regex: Option<&Regex>,
        ) -> bool {
            if let Some(needle_upper) = label_contains_upper {
                let matched = labels.iter().any(|value| {
                    let upper = value.to_ascii_uppercase();
                    upper == *needle_upper || upper.contains(needle_upper)
                });
                if !matched {
                    return false;
                }
            }
            if let Some(regex) = label_regex
                && !labels.iter().any(|value| regex.is_match(value))
            {
                return false;
            }
            true
        }

        fn range_matches(
            start_0based: usize,
            end_0based_exclusive: usize,
            query: &SequenceFeatureQuery,
        ) -> bool {
            let Some(range_start) = query.start_0based else {
                return true;
            };
            let Some(range_end) = query.end_0based_exclusive else {
                return true;
            };
            match query.range_relation {
                SequenceFeatureRangeRelation::Overlap => {
                    end_0based_exclusive > range_start && start_0based < range_end
                }
                SequenceFeatureRangeRelation::Within => {
                    start_0based >= range_start && end_0based_exclusive <= range_end
                }
                SequenceFeatureRangeRelation::Contains => {
                    start_0based <= range_start && end_0based_exclusive >= range_end
                }
            }
        }

        fn strand_matches(strand: char, filter: SequenceFeatureStrandFilter) -> bool {
            match filter {
                SequenceFeatureStrandFilter::Any => true,
                SequenceFeatureStrandFilter::Forward => strand == '+',
                SequenceFeatureStrandFilter::Reverse => strand == '-',
            }
        }

        fn parse_bed_score(raw: &str) -> Option<usize> {
            raw.trim()
                .parse::<f64>()
                .ok()
                .map(|value| value.round().clamp(0.0, 1000.0) as usize)
        }

        fn parse_quantile_score(raw: &str) -> Option<usize> {
            raw.trim()
                .parse::<f64>()
                .ok()
                .map(|value| (value.clamp(0.0, 1.0) * 1000.0).round() as usize)
        }

        fn feature_bed_score(feature: &gb_io::seq::Feature) -> usize {
            for key in ["score", "bed_score"] {
                if let Some(score) = feature.qualifier_values(key).find_map(parse_bed_score) {
                    return score;
                }
            }
            for key in ["llr_quantile", "true_log_odds_quantile"] {
                if let Some(score) = feature.qualifier_values(key).find_map(parse_quantile_score) {
                    return score;
                }
            }
            0
        }

        fn genomic_bed_coordinates(
            feature: &gb_io::seq::Feature,
        ) -> Option<(String, usize, usize)> {
            let chromosome = feature.qualifier_values("chromosome").find_map(|value| {
                let trimmed = value.trim();
                (!trimmed.is_empty()).then_some(trimmed.to_string())
            })?;
            let start_1based = feature
                .qualifier_values("genomic_start_1based")
                .find_map(|value| value.trim().parse::<usize>().ok())?;
            let end_1based = feature
                .qualifier_values("genomic_end_1based")
                .find_map(|value| value.trim().parse::<usize>().ok())?;
            if start_1based == 0 || end_1based < start_1based {
                return None;
            }
            Some((chromosome, start_1based.saturating_sub(1), end_1based))
        }

        fn restriction_end_geometry_label(enzyme: &RestrictionEnzyme) -> &'static str {
            match enzyme.end_geometry() {
                crate::restriction_enzyme::RestrictionEndGeometry::Blunt => "blunt",
                crate::restriction_enzyme::RestrictionEndGeometry::FivePrimeOverhang(_) => {
                    "5prime_overhang"
                }
                crate::restriction_enzyme::RestrictionEndGeometry::ThreePrimeOverhang(_) => {
                    "3prime_overhang"
                }
            }
        }

        query.seq_id = query.seq_id.trim().to_string();
        if query.seq_id.is_empty() {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: "export_sequence_features_bed requires non-empty seq_id".to_string(),
            });
        }
        if query.limit.is_some_and(|limit| limit == 0) {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: "Feature BED export limit must be >= 1".to_string(),
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

        let mut all_feature_rows: Vec<SequenceFeatureQueryRow> = vec![];
        let mut page_offset = 0usize;
        loop {
            let mut page_query = query.clone();
            page_query.offset = page_offset;
            page_query.limit = Some(FEATURE_QUERY_MAX_LIMIT);
            let page = self.query_sequence_features(page_query)?;
            if page.rows.is_empty() {
                break;
            }
            page_offset = page_offset.saturating_add(page.rows.len());
            all_feature_rows.extend(page.rows);
            if page_offset >= page.matched_count {
                break;
            }
        }
        let matched_sequence_feature_count = all_feature_rows.len();

        let normalized_kind_in = query
            .kind_in
            .iter()
            .filter_map(|value| Self::normalize_feature_kind_token(value))
            .collect::<HashSet<_>>();
        let normalized_kind_not_in = query
            .kind_not_in
            .iter()
            .filter_map(|value| Self::normalize_feature_kind_token(value))
            .collect::<HashSet<_>>();
        let restriction_kind_allowed = (normalized_kind_in.is_empty()
            || normalized_kind_in.contains("RESTRICTION_SITE"))
            && !normalized_kind_not_in.contains("RESTRICTION_SITE");
        let label_contains_upper = query
            .label_contains
            .as_deref()
            .map(str::trim)
            .filter(|value| !value.is_empty())
            .map(|value| value.to_ascii_uppercase());
        let label_regex = query
            .label_regex
            .as_deref()
            .map(str::trim)
            .filter(|value| !value.is_empty())
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
        let qualifier_filters = compile_qualifier_filters(&query)?;

        let mut normalized_restriction_filters: Vec<String> = vec![];
        let mut seen_restriction_filters: HashSet<String> = HashSet::new();
        for raw in restriction_enzymes {
            let trimmed = raw.trim();
            if trimmed.is_empty() {
                continue;
            }
            let normalized = trimmed.to_ascii_uppercase();
            if seen_restriction_filters.insert(normalized) {
                normalized_restriction_filters.push(trimmed.to_string());
            }
        }

        let mut skipped_missing_genomic_coordinates = 0usize;
        let mut export_rows: Vec<FeatureBedExportRow> = vec![];
        for row in all_feature_rows {
            let feature = dna
                .features()
                .get(row.feature_id)
                .ok_or_else(|| EngineError {
                    code: ErrorCode::Internal,
                    message: format!(
                        "Feature query row referenced missing feature_id {} on '{}'",
                        row.feature_id, query.seq_id
                    ),
                })?;
            let qualifiers = collect_qualifiers(feature);
            let (chrom, chrom_start_0based, chrom_end_0based_exclusive, coordinate_source) =
                match coordinate_mode {
                    FeatureBedCoordinateMode::Auto => {
                        if let Some((chrom, start, end)) = genomic_bed_coordinates(feature) {
                            (chrom, start, end, "genomic")
                        } else {
                            (
                                query.seq_id.clone(),
                                row.start_0based,
                                row.end_0based_exclusive,
                                "local",
                            )
                        }
                    }
                    FeatureBedCoordinateMode::Local => (
                        query.seq_id.clone(),
                        row.start_0based,
                        row.end_0based_exclusive,
                        "local",
                    ),
                    FeatureBedCoordinateMode::Genomic => {
                        let Some((chrom, start, end)) = genomic_bed_coordinates(feature) else {
                            skipped_missing_genomic_coordinates =
                                skipped_missing_genomic_coordinates.saturating_add(1);
                            continue;
                        };
                        (chrom, start, end, "genomic")
                    }
                };
            let qualifiers_json = serde_json::to_string(&qualifiers).map_err(|e| EngineError {
                code: ErrorCode::Internal,
                message: format!(
                    "Could not serialize qualifiers for feature {} on '{}': {e}",
                    row.feature_id, query.seq_id
                ),
            })?;
            export_rows.push(FeatureBedExportRow {
                chrom,
                chrom_start_0based,
                chrom_end_0based_exclusive,
                name: row.label.clone(),
                score_0_to_1000: feature_bed_score(feature),
                strand: if row.strand.eq_ignore_ascii_case("reverse") {
                    '-'
                } else {
                    '+'
                },
                kind: row.kind.clone(),
                row_id: format!("feature:{}", row.feature_id),
                coordinate_source,
                qualifiers_json,
                sort_feature_id: row.feature_id,
                length_bp: row.length_bp,
            });
        }

        let mut matched_restriction_site_count = 0usize;
        if include_restriction_sites && restriction_kind_allowed {
            let restriction_filter_set = normalized_restriction_filters
                .iter()
                .map(|value| value.to_ascii_uppercase())
                .collect::<HashSet<_>>();
            for (site_idx, site) in dna.restriction_enzyme_sites().iter().enumerate() {
                if !restriction_filter_set.is_empty()
                    && !restriction_filter_set.contains(&site.enzyme.name.to_ascii_uppercase())
                {
                    continue;
                }
                let Some((start_0based, end_0based_exclusive)) =
                    site.recognition_bounds_0based(dna.len())
                else {
                    continue;
                };
                if end_0based_exclusive <= start_0based || start_0based >= dna.len() {
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
                if !range_matches(start_0based, end_0based_exclusive, &query) {
                    continue;
                }
                let strand = if site.forward_strand { '+' } else { '-' };
                if !strand_matches(strand, query.strand) {
                    continue;
                }
                let labels = vec![site.enzyme.name.clone(), site.enzyme.sequence.clone()];
                if !label_matches(&labels, label_contains_upper.as_ref(), label_regex.as_ref()) {
                    continue;
                }
                let (forward_cut, reverse_cut) = site
                    .strand_cut_positions_0based(dna.len())
                    .unwrap_or((start_0based, start_0based));
                let mut qualifiers: BTreeMap<String, Vec<String>> = BTreeMap::new();
                for (key, value) in [
                    ("enzyme", site.enzyme.name.clone()),
                    ("recognition_sequence", site.enzyme.sequence.clone()),
                    (
                        "end_geometry",
                        restriction_end_geometry_label(&site.enzyme).to_string(),
                    ),
                    ("forward_cut_0based", forward_cut.to_string()),
                    ("reverse_cut_0based", reverse_cut.to_string()),
                    ("gentle_generated", "restriction_site".to_string()),
                ] {
                    qualifiers.entry(key.to_string()).or_default().push(value);
                }
                if let Some(note) = site
                    .enzyme
                    .note
                    .as_deref()
                    .map(str::trim)
                    .filter(|value| !value.is_empty())
                {
                    qualifiers
                        .entry("note".to_string())
                        .or_default()
                        .push(note.to_string());
                }
                if !qualifiers_match(&qualifiers, &qualifier_filters) {
                    continue;
                }
                matched_restriction_site_count = matched_restriction_site_count.saturating_add(1);
                let (chrom, chrom_start_0based, chrom_end_0based_exclusive, coordinate_source) =
                    match coordinate_mode {
                        FeatureBedCoordinateMode::Auto | FeatureBedCoordinateMode::Local => (
                            query.seq_id.clone(),
                            start_0based,
                            end_0based_exclusive,
                            "local",
                        ),
                        FeatureBedCoordinateMode::Genomic => {
                            skipped_missing_genomic_coordinates =
                                skipped_missing_genomic_coordinates.saturating_add(1);
                            continue;
                        }
                    };
                let qualifiers_json =
                    serde_json::to_string(&qualifiers).map_err(|e| EngineError {
                        code: ErrorCode::Internal,
                        message: format!(
                            "Could not serialize restriction-site qualifiers for '{}' on '{}': {e}",
                            site.enzyme.name, query.seq_id
                        ),
                    })?;
                export_rows.push(FeatureBedExportRow {
                    chrom,
                    chrom_start_0based,
                    chrom_end_0based_exclusive,
                    name: site.enzyme.name.clone(),
                    score_0_to_1000: 0,
                    strand,
                    kind: "restriction_site".to_string(),
                    row_id: format!(
                        "restriction_site:{}:{}:{}:{}",
                        Self::normalize_id_token(&site.enzyme.name),
                        start_0based,
                        end_0based_exclusive,
                        site_idx
                    ),
                    coordinate_source,
                    qualifiers_json,
                    sort_feature_id: usize::MAX / 2usize + site_idx,
                    length_bp,
                });
            }
        }

        let matched_row_count =
            matched_sequence_feature_count.saturating_add(matched_restriction_site_count);
        export_rows.sort_by(|a, b| match query.sort_by {
            SequenceFeatureSortBy::FeatureId => a
                .sort_feature_id
                .cmp(&b.sort_feature_id)
                .then(a.chrom_start_0based.cmp(&b.chrom_start_0based))
                .then(
                    a.chrom_end_0based_exclusive
                        .cmp(&b.chrom_end_0based_exclusive),
                )
                .then(a.row_id.cmp(&b.row_id)),
            SequenceFeatureSortBy::Start => a
                .chrom_start_0based
                .cmp(&b.chrom_start_0based)
                .then(
                    a.chrom_end_0based_exclusive
                        .cmp(&b.chrom_end_0based_exclusive),
                )
                .then(a.sort_feature_id.cmp(&b.sort_feature_id))
                .then(a.row_id.cmp(&b.row_id)),
            SequenceFeatureSortBy::End => a
                .chrom_end_0based_exclusive
                .cmp(&b.chrom_end_0based_exclusive)
                .then(a.chrom_start_0based.cmp(&b.chrom_start_0based))
                .then(a.sort_feature_id.cmp(&b.sort_feature_id))
                .then(a.row_id.cmp(&b.row_id)),
            SequenceFeatureSortBy::Kind => a
                .kind
                .to_ascii_uppercase()
                .cmp(&b.kind.to_ascii_uppercase())
                .then(a.chrom_start_0based.cmp(&b.chrom_start_0based))
                .then(
                    a.chrom_end_0based_exclusive
                        .cmp(&b.chrom_end_0based_exclusive),
                )
                .then(a.sort_feature_id.cmp(&b.sort_feature_id))
                .then(a.row_id.cmp(&b.row_id)),
            SequenceFeatureSortBy::Length => a
                .length_bp
                .cmp(&b.length_bp)
                .then(a.chrom_start_0based.cmp(&b.chrom_start_0based))
                .then(
                    a.chrom_end_0based_exclusive
                        .cmp(&b.chrom_end_0based_exclusive),
                )
                .then(a.sort_feature_id.cmp(&b.sort_feature_id))
                .then(a.row_id.cmp(&b.row_id)),
        });
        if query.descending {
            export_rows.reverse();
        }

        let exportable_row_count = export_rows.len();
        let offset = query.offset.min(exportable_row_count);
        let selected_rows = match query.limit {
            Some(limit) => export_rows
                .into_iter()
                .skip(offset)
                .take(limit)
                .collect::<Vec<_>>(),
            None => export_rows.into_iter().skip(offset).collect::<Vec<_>>(),
        };
        let local_coordinate_row_count = selected_rows
            .iter()
            .filter(|row| row.coordinate_source == "local")
            .count();
        let genomic_coordinate_row_count = selected_rows
            .iter()
            .filter(|row| row.coordinate_source == "genomic")
            .count();

        let file = File::create(path).map_err(|e| EngineError {
            code: ErrorCode::Io,
            message: format!("Could not create BED export '{}': {e}", path),
        })?;
        let mut writer = BufWriter::new(file);
        for row in &selected_rows {
            writeln!(
                writer,
                "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
                bed_safe_field(&row.chrom),
                row.chrom_start_0based,
                row.chrom_end_0based_exclusive,
                bed_safe_field(&row.name),
                row.score_0_to_1000,
                row.strand,
                bed_safe_field(&row.kind),
                bed_safe_field(&row.row_id),
                row.coordinate_source,
                bed_safe_field(&row.qualifiers_json),
            )
            .map_err(|e| EngineError {
                code: ErrorCode::Io,
                message: format!("Could not write BED export '{}': {e}", path),
            })?;
        }
        writer.flush().map_err(|e| EngineError {
            code: ErrorCode::Io,
            message: format!("Could not flush BED export '{}': {e}", path),
        })?;

        Ok(SequenceFeatureBedExportReport {
            schema: FEATURE_BED_EXPORT_REPORT_SCHEMA.to_string(),
            seq_id: query.seq_id.clone(),
            path: path.to_string(),
            coordinate_mode: coordinate_mode.as_str().to_string(),
            include_restriction_sites,
            restriction_enzyme_filters: normalized_restriction_filters,
            bed_columns: Self::sequence_feature_bed_columns(),
            matched_sequence_feature_count,
            matched_restriction_site_count,
            matched_row_count,
            exportable_row_count,
            exported_row_count: selected_rows.len(),
            offset,
            limit: query.limit,
            local_coordinate_row_count,
            genomic_coordinate_row_count,
            skipped_missing_genomic_coordinates,
            query,
        })
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

    fn normalize_routine_family_token(raw: &str) -> String {
        let mut out = String::with_capacity(raw.len());
        let mut last_was_separator = false;
        for ch in raw.trim().chars() {
            if ch.is_ascii_alphanumeric() {
                out.push(ch.to_ascii_lowercase());
                last_was_separator = false;
            } else if !last_was_separator {
                out.push('_');
                last_was_separator = true;
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
        objective.helper_profile_id =
            Self::normalize_optional_id_token_text(objective.helper_profile_id.take());
        let mut preferred_families = objective
            .preferred_routine_families
            .iter()
            .map(|value| Self::normalize_routine_family_token(value))
            .filter(|value| !value.is_empty())
            .collect::<BTreeSet<_>>()
            .into_iter()
            .collect::<Vec<_>>();
        preferred_families.sort_by_key(|value| value.to_ascii_lowercase());
        objective.preferred_routine_families = preferred_families;
        objective
    }

    fn normalize_routine_family_preferences(values: &[String]) -> Vec<String> {
        let mut families = values
            .iter()
            .map(|value| Self::normalize_routine_family_token(value))
            .filter(|value| !value.is_empty())
            .collect::<BTreeSet<_>>()
            .into_iter()
            .collect::<Vec<_>>();
        families.sort_by_key(|value| value.to_ascii_lowercase());
        families
    }

    fn helper_construct_routine_preference_haystack(
        helper_interpretation: &HelperConstructInterpretation,
    ) -> String {
        let mut fields = vec![
            helper_interpretation.helper_id.as_str(),
            helper_interpretation
                .description
                .as_deref()
                .unwrap_or_default(),
            helper_interpretation.summary.as_deref().unwrap_or_default(),
        ];
        for alias in &helper_interpretation.aliases {
            fields.push(alias.as_str());
        }
        for kind in &helper_interpretation.helper_kinds {
            fields.push(kind.as_str());
        }
        for host in &helper_interpretation.host_systems {
            fields.push(host.as_str());
        }
        for function in &helper_interpretation.offered_functions {
            fields.push(function.as_str());
        }
        for component in &helper_interpretation.components {
            fields.push(component.id.as_str());
            fields.push(component.kind.as_str());
            fields.push(component.label.as_deref().unwrap_or_default());
            fields.push(component.description.as_deref().unwrap_or_default());
            for tag in &component.tags {
                fields.push(tag.as_str());
            }
            for value in component.attributes.values() {
                fields.push(value.as_str());
            }
        }
        for relationship in &helper_interpretation.relationships {
            fields.push(relationship.subject.as_str());
            fields.push(relationship.predicate.as_str());
            fields.push(relationship.object.as_str());
            fields.push(relationship.note.as_deref().unwrap_or_default());
        }
        fields.join(" ").to_ascii_lowercase()
    }

    fn helper_construct_haystack_contains_any(haystack: &str, terms: &[&str]) -> bool {
        let lower = haystack.to_ascii_lowercase();
        terms
            .iter()
            .any(|term| lower.contains(&term.to_ascii_lowercase()))
    }

    fn helper_construct_derived_routine_preferences(
        helper_interpretation: &HelperConstructInterpretation,
    ) -> (Vec<String>, Vec<String>) {
        let haystack = Self::helper_construct_routine_preference_haystack(helper_interpretation);
        let offered_functions = helper_interpretation
            .offered_functions
            .iter()
            .map(|value| Self::normalize_planning_class_key(value))
            .collect::<BTreeSet<_>>();
        let component_kinds = helper_interpretation
            .components
            .iter()
            .map(|component| Self::normalize_planning_class_key(&component.kind))
            .collect::<BTreeSet<_>>();
        let component_tags = helper_interpretation
            .components
            .iter()
            .flat_map(|component| component.tags.iter())
            .map(|value| Self::normalize_planning_class_key(value))
            .collect::<BTreeSet<_>>();

        let has_insert_cloning = offered_functions.contains("insert_cloning")
            || component_kinds.contains("cloning_site")
            || component_tags.contains("insert_cloning");
        let has_in_frame_fusion_context = offered_functions.contains("in_frame_fusion_design")
            || offered_functions.contains("fusion_tagging")
            || offered_functions.contains("protease_cleavage")
            || offered_functions.contains("protease_tag_removal");

        let mut families = BTreeSet::new();
        let mut rationale = BTreeSet::new();

        if Self::helper_construct_haystack_contains_any(
            &haystack,
            &["gateway", "clonase", "attl", "attr", "attb", "attp"],
        ) {
            families.insert("gateway".to_string());
            rationale.insert(
                "Helper interpretation mentions Gateway-specific cloning tokens, so Gateway routines stay near the top of routine ranking.".to_string(),
            );
        }
        if Self::helper_construct_haystack_contains_any(
            &haystack,
            &["golden gate", "type_iis", "type iis", "type-iis"],
        ) {
            families.insert("golden_gate".to_string());
            rationale.insert(
                "Helper interpretation mentions Type IIS / Golden Gate grammar, so Golden Gate routines are treated as a compatible fit.".to_string(),
            );
        }
        if Self::helper_construct_haystack_contains_any(&haystack, &["gibson", "gibson assembly"]) {
            families.insert("gibson".to_string());
            rationale.insert(
                "Helper interpretation explicitly mentions Gibson-style overlap assembly."
                    .to_string(),
            );
        }
        if Self::helper_construct_haystack_contains_any(&haystack, &["in-fusion", "infusion"]) {
            families.insert("infusion".to_string());
            rationale.insert(
                "Helper interpretation explicitly mentions In-Fusion-style overlap assembly."
                    .to_string(),
            );
        }
        if Self::helper_construct_haystack_contains_any(
            &haystack,
            &["nebuilder", "hi-fi assembly", "hifi assembly"],
        ) {
            families.insert("nebuilder_hifi".to_string());
            rationale.insert(
                "Helper interpretation explicitly mentions NEBuilder HiFi-style overlap assembly."
                    .to_string(),
            );
        }
        if Self::helper_construct_haystack_contains_any(&haystack, &["topo cloning", "topo vector"])
        {
            families.insert("topo".to_string());
            rationale
                .insert("Helper interpretation mentions TOPO-style cloning semantics.".to_string());
        }
        if Self::helper_construct_haystack_contains_any(
            &haystack,
            &[
                "ta cloning",
                "t-overhang",
                "t overhang",
                "a-overhang",
                "a overhang",
            ],
        ) {
            families.insert("ta_gc".to_string());
            rationale.insert(
                "Helper interpretation mentions TA-style overhang cloning semantics.".to_string(),
            );
        }
        if has_insert_cloning
            || Self::helper_construct_haystack_contains_any(
                &haystack,
                &[
                    "multiple cloning site",
                    "mcs",
                    "restriction",
                    "digest",
                    "ligation",
                ],
            )
        {
            families.insert("restriction".to_string());
            rationale.insert(
                "Helper interpretation exposes a cloning-site / insert-cloning context, so restriction-family subcloning stays explicitly favored.".to_string(),
            );
        }
        if has_in_frame_fusion_context {
            families.insert("gibson".to_string());
            families.insert("infusion".to_string());
            families.insert("nebuilder_hifi".to_string());
            rationale.insert(
                "Helper interpretation exposes in-frame fusion / protease-tag context, so overlap-assembly families are favored for scar-aware CDS insertion planning.".to_string(),
            );
        }

        (
            families.into_iter().collect(),
            rationale.into_iter().collect(),
        )
    }

    fn variant_assay_routine_family_preferences(assay_id: &str) -> Vec<String> {
        match assay_id.trim() {
            "allele_paired_promoter_luciferase_reporter"
            | "allele_paired_regulatory_reporter"
            | "allele_paired_expression_compare"
            | "utr_reporter_translation_compare" => vec![
                "gibson".to_string(),
                "infusion".to_string(),
                "nebuilder_hifi".to_string(),
                "restriction".to_string(),
                "golden_gate".to_string(),
            ],
            "minigene_splicing_assay" => vec![
                "gibson".to_string(),
                "infusion".to_string(),
                "nebuilder_hifi".to_string(),
                "restriction".to_string(),
                "golden_gate".to_string(),
                "pcr".to_string(),
            ],
            _ => vec![],
        }
    }

    fn variant_assay_ids_to_routine_family_preferences(assay_ids: &[String]) -> Vec<String> {
        assay_ids
            .iter()
            .flat_map(|assay_id| Self::variant_assay_routine_family_preferences(assay_id))
            .collect::<BTreeSet<_>>()
            .into_iter()
            .collect::<Vec<_>>()
    }

    fn construct_reasoning_variant_preference_summary(
        graph: Option<&ConstructReasoningGraph>,
    ) -> (
        Option<String>,
        Vec<String>,
        Vec<String>,
        Vec<String>,
        Option<String>,
    ) {
        let Some(graph) = graph else {
            return (None, vec![], vec![], vec![], None);
        };
        let effect_tags = graph
            .facts
            .iter()
            .find(|fact| fact.fact_type == "variant_effect_context")
            .and_then(|fact| {
                fact.value_json
                    .get("effect_tags")
                    .and_then(serde_json::Value::as_array)
                    .map(|rows| {
                        rows.iter()
                            .filter_map(serde_json::Value::as_str)
                            .map(|value| value.to_string())
                            .collect::<BTreeSet<_>>()
                            .into_iter()
                            .collect::<Vec<_>>()
                    })
            })
            .unwrap_or_default();
        let assay_ids = graph
            .facts
            .iter()
            .find(|fact| fact.fact_type == "variant_assay_context")
            .and_then(|fact| {
                fact.value_json
                    .get("suggested_assay_ids")
                    .and_then(serde_json::Value::as_array)
                    .map(|rows| {
                        rows.iter()
                            .filter_map(serde_json::Value::as_str)
                            .map(|value| value.to_string())
                            .collect::<BTreeSet<_>>()
                            .into_iter()
                            .collect::<Vec<_>>()
                    })
            })
            .unwrap_or_default();
        let routine_families = Self::variant_assay_ids_to_routine_family_preferences(&assay_ids);
        let rationale = if assay_ids.is_empty() {
            None
        } else {
            Some(format!(
                "Construct reasoning for '{}' suggests assay family/ies {} from variant context (effect tags: {}).",
                graph.seq_id,
                assay_ids.join(", "),
                if effect_tags.is_empty() {
                    "none".to_string()
                } else {
                    effect_tags.join(", ")
                }
            ))
        };
        (
            Some(graph.seq_id.clone()),
            effect_tags,
            assay_ids,
            routine_families,
            rationale,
        )
    }

    fn construct_reasoning_strategy_preference_summary(
        graph: Option<&ConstructReasoningGraph>,
    ) -> (Option<String>, Vec<String>, Option<String>) {
        let Some(graph) = graph else {
            return (None, vec![], None);
        };
        let strategy_families = graph
            .facts
            .iter()
            .find(|fact| fact.fact_type == "adapter_restriction_capture_context")
            .and_then(|fact| {
                fact.value_json
                    .get("derived_preferred_routine_families")
                    .and_then(serde_json::Value::as_array)
                    .map(|rows| {
                        rows.iter()
                            .filter_map(serde_json::Value::as_str)
                            .map(|value| value.to_string())
                            .collect::<BTreeSet<_>>()
                            .into_iter()
                            .collect::<Vec<_>>()
                    })
            })
            .unwrap_or_default();
        let rationale = if strategy_families.is_empty() {
            None
        } else {
            Some(format!(
                "Construct reasoning for '{}' suggests routine family/ies {} from adapter/linker capture context.",
                graph.seq_id,
                strategy_families.join(", ")
            ))
        };
        (Some(graph.seq_id.clone()), strategy_families, rationale)
    }

    fn build_routine_preference_context(
        helper_profile_id: Option<&str>,
        explicit_preferred_routine_families: &[String],
        helper_interpretation: Option<&HelperConstructInterpretation>,
        helper_resolution_status: &str,
        helper_resolution_note: Option<String>,
        construct_reasoning_seq_id: Option<&str>,
        construct_strategy_derived_preferred_routine_families: &[String],
        variant_effect_tags: &[String],
        variant_suggested_assay_ids: &[String],
        variant_derived_preferred_routine_families: &[String],
        construct_strategy_reasoning_note: Option<String>,
        variant_reasoning_note: Option<String>,
    ) -> RoutinePreferenceContext {
        let explicit_preferred_routine_families =
            Self::normalize_routine_family_preferences(explicit_preferred_routine_families);
        let construct_strategy_derived_preferred_routine_families =
            Self::normalize_routine_family_preferences(
                construct_strategy_derived_preferred_routine_families,
            );
        let variant_derived_preferred_routine_families =
            Self::normalize_routine_family_preferences(variant_derived_preferred_routine_families);
        let helper_profile_id = helper_profile_id
            .map(str::trim)
            .filter(|value| !value.is_empty())
            .map(|value| value.to_string());
        let construct_reasoning_seq_id = construct_reasoning_seq_id
            .map(str::trim)
            .filter(|value| !value.is_empty())
            .map(|value| value.to_string());
        let helper_offered_functions = helper_interpretation
            .map(|row| row.offered_functions.clone())
            .unwrap_or_default();
        let helper_component_labels = helper_interpretation
            .map(|row| {
                row.components
                    .iter()
                    .map(|component| {
                        component
                            .label
                            .clone()
                            .unwrap_or_else(|| component.id.clone())
                    })
                    .collect::<Vec<_>>()
            })
            .unwrap_or_default();
        let (helper_derived_preferred_routine_families, mut rationale) = helper_interpretation
            .map(Self::helper_construct_derived_routine_preferences)
            .unwrap_or_else(|| (vec![], vec![]));
        if !explicit_preferred_routine_families.is_empty() {
            rationale.push(format!(
                "Objective explicitly prefers routine families: {}.",
                explicit_preferred_routine_families.join(", ")
            ));
        }
        if let Some(note) = helper_resolution_note {
            rationale.push(note);
        }
        if let Some(note) = construct_strategy_reasoning_note {
            rationale.push(note);
        }
        if let Some(note) = variant_reasoning_note {
            rationale.push(note);
        }
        rationale.sort();
        rationale.dedup();
        let mut effective_preferred_routine_families = explicit_preferred_routine_families
            .iter()
            .chain(helper_derived_preferred_routine_families.iter())
            .chain(construct_strategy_derived_preferred_routine_families.iter())
            .chain(variant_derived_preferred_routine_families.iter())
            .cloned()
            .collect::<BTreeSet<_>>()
            .into_iter()
            .collect::<Vec<_>>();
        effective_preferred_routine_families.sort_by_key(|value| value.to_ascii_lowercase());

        RoutinePreferenceContext {
            helper_profile_id,
            construct_reasoning_seq_id,
            helper_resolution_status: helper_resolution_status.to_string(),
            explicit_preferred_routine_families,
            helper_derived_preferred_routine_families,
            construct_strategy_derived_preferred_routine_families,
            variant_derived_preferred_routine_families,
            effective_preferred_routine_families,
            helper_offered_functions,
            helper_component_labels,
            variant_effect_tags: variant_effect_tags.to_vec(),
            variant_suggested_assay_ids: variant_suggested_assay_ids.to_vec(),
            rationale,
        }
    }

    fn routine_preference_context_record(
        context: &RoutinePreferenceContext,
    ) -> RoutinePreferenceContextRecord {
        RoutinePreferenceContextRecord {
            helper_profile_id: context.helper_profile_id.clone(),
            construct_reasoning_seq_id: context.construct_reasoning_seq_id.clone(),
            helper_resolution_status: context.helper_resolution_status.clone(),
            explicit_preferred_routine_families: context
                .explicit_preferred_routine_families
                .clone(),
            helper_derived_preferred_routine_families: context
                .helper_derived_preferred_routine_families
                .clone(),
            construct_strategy_derived_preferred_routine_families: context
                .construct_strategy_derived_preferred_routine_families
                .clone(),
            variant_derived_preferred_routine_families: context
                .variant_derived_preferred_routine_families
                .clone(),
            effective_preferred_routine_families: context
                .effective_preferred_routine_families
                .clone(),
            helper_offered_functions: context.helper_offered_functions.clone(),
            helper_component_labels: context.helper_component_labels.clone(),
            variant_effect_tags: context.variant_effect_tags.clone(),
            variant_suggested_assay_ids: context.variant_suggested_assay_ids.clone(),
            rationale: context.rationale.clone(),
        }
    }

    pub(crate) fn planning_objective_routine_preference_context(
        objective: &PlanningObjective,
    ) -> RoutinePreferenceContext {
        let explicit_preferred_routine_families =
            Self::normalize_routine_family_preferences(&objective.preferred_routine_families);
        let Some(helper_profile_id) = objective.helper_profile_id.as_deref() else {
            return Self::build_routine_preference_context(
                None,
                &explicit_preferred_routine_families,
                None,
                "not_requested",
                None,
                None,
                &[],
                &[],
                &[],
                &[],
                None,
                None,
            );
        };
        match Self::interpret_helper_genome(helper_profile_id, None) {
            Ok(Some(helper_interpretation)) => Self::build_routine_preference_context(
                Some(helper_profile_id),
                &explicit_preferred_routine_families,
                Some(&helper_interpretation),
                "resolved",
                None,
                None,
                &[],
                &[],
                &[],
                &[],
                None,
                None,
            ),
            Ok(None) => Self::build_routine_preference_context(
                Some(helper_profile_id),
                &explicit_preferred_routine_families,
                None,
                "not_found",
                Some(format!(
                    "Planning objective selected helper profile '{helper_profile_id}', but no helper interpretation was found in the active helper catalog."
                )),
                None,
                &[],
                &[],
                &[],
                &[],
                None,
                None,
            ),
            Err(error) => Self::build_routine_preference_context(
                Some(helper_profile_id),
                &explicit_preferred_routine_families,
                None,
                "error",
                Some(format!(
                    "Planning objective could not resolve helper profile '{helper_profile_id}' for routine-preference synthesis: {}.",
                    error.message
                )),
                None,
                &[],
                &[],
                &[],
                &[],
                None,
                None,
            ),
        }
    }

    fn construct_objective_routine_preference_context(
        objective: &ConstructObjective,
        helper_interpretation: Option<&HelperConstructInterpretation>,
        construct_reasoning_seq_id: Option<&str>,
        construct_strategy_derived_preferred_routine_families: &[String],
        variant_effect_tags: &[String],
        variant_suggested_assay_ids: &[String],
        variant_derived_preferred_routine_families: &[String],
    ) -> RoutinePreferenceContext {
        let explicit_preferred_routine_families =
            Self::normalize_routine_family_preferences(&objective.preferred_routine_families);
        let helper_resolution_status = match (
            objective.helper_profile_id.as_ref(),
            helper_interpretation.is_some(),
        ) {
            (None, _) => "not_requested",
            (Some(_), true) => "resolved",
            (Some(_), false) => "not_found",
        };
        let helper_resolution_note = match (
            objective.helper_profile_id.as_ref(),
            helper_interpretation.is_some(),
        ) {
            (Some(helper_profile_id), false) => Some(format!(
                "Construct objective selected helper profile '{helper_profile_id}', but no helper interpretation was available for routine-preference synthesis."
            )),
            _ => None,
        };
        Self::build_routine_preference_context(
            objective.helper_profile_id.as_deref(),
            &explicit_preferred_routine_families,
            helper_interpretation,
            helper_resolution_status,
            helper_resolution_note,
            construct_reasoning_seq_id,
            construct_strategy_derived_preferred_routine_families,
            variant_effect_tags,
            variant_suggested_assay_ids,
            variant_derived_preferred_routine_families,
            (!construct_strategy_derived_preferred_routine_families.is_empty()).then(|| {
                format!(
                    "Construct reasoning records adapter/linker capture strategy-derived routine family/ies {} for sequence '{}'.",
                    construct_strategy_derived_preferred_routine_families.join(", "),
                    construct_reasoning_seq_id.unwrap_or_default()
                )
            }),
            (!variant_suggested_assay_ids.is_empty()).then(|| {
                format!(
                    "Variant-aware construct reasoning suggests assay family/ies {} for sequence '{}'.",
                    variant_suggested_assay_ids.join(", "),
                    construct_reasoning_seq_id.unwrap_or_default()
                )
            }),
        )
    }

    fn macro_template_family_match_terms(family: &str) -> Vec<String> {
        let normalized = Self::normalize_routine_family_token(family);
        let mut terms = BTreeSet::new();
        if !normalized.is_empty() {
            terms.insert(normalized.clone());
            terms.insert(normalized.replace('_', " "));
        }
        match normalized.as_str() {
            "golden_gate" => {
                terms.insert("golden gate".to_string());
                terms.insert("type_iis".to_string());
                terms.insert("type iis".to_string());
            }
            "gibson" => {
                terms.insert("gibson".to_string());
                terms.insert("overlap assembly".to_string());
            }
            "restriction" => {
                terms.insert("restriction".to_string());
                terms.insert("digest".to_string());
                terms.insert("ligation".to_string());
            }
            "gateway" => {
                terms.insert("gateway".to_string());
                terms.insert("clonase".to_string());
                terms.insert("attl".to_string());
                terms.insert("attr".to_string());
            }
            "infusion" => {
                terms.insert("in-fusion".to_string());
                terms.insert("infusion".to_string());
            }
            "nebuilder_hifi" => {
                terms.insert("nebuilder".to_string());
                terms.insert("hifi".to_string());
                terms.insert("hi fi".to_string());
            }
            "topo" => {
                terms.insert("topo".to_string());
                terms.insert("topo cloning".to_string());
            }
            "ta_gc" => {
                terms.insert("ta cloning".to_string());
                terms.insert("t overhang".to_string());
                terms.insert("a overhang".to_string());
            }
            _ => {}
        }
        terms.into_iter().collect()
    }

    fn macro_template_suggestion_haystack(
        name: &str,
        description: Option<&str>,
        script: &str,
    ) -> String {
        format!(
            "{} {} {}",
            name.trim(),
            description.unwrap_or_default().trim(),
            script.trim()
        )
        .to_ascii_lowercase()
    }

    fn routine_macro_id_match_terms(routine_id: &str) -> Vec<String> {
        let compact = routine_id.trim();
        if compact.is_empty() {
            return vec![];
        }
        let mut terms = BTreeSet::new();
        terms.insert(compact.to_ascii_lowercase());
        terms.insert(compact.replace('.', "_").to_ascii_lowercase());
        terms.insert(
            compact
                .replace('.', " ")
                .replace('_', " ")
                .to_ascii_lowercase(),
        );
        terms.into_iter().collect()
    }

    fn score_macro_template_suggestion(
        macro_kind: &str,
        template_name: &str,
        description: Option<&str>,
        details_url: Option<&str>,
        script: &str,
        selected_routine_id: Option<&str>,
        selected_routine_family: Option<&str>,
        context: &RoutinePreferenceContext,
    ) -> Option<MacroTemplateSuggestion> {
        let haystack = Self::macro_template_suggestion_haystack(template_name, description, script);
        let mut score = 0.0;
        let mut matched_routine_families = BTreeSet::new();
        let mut matched_terms = BTreeSet::new();
        let mut rationale = BTreeSet::new();

        if let Some(selected_family) = selected_routine_family {
            let normalized_selected_family = Self::normalize_routine_family_token(selected_family);
            let family_terms = Self::macro_template_family_match_terms(&normalized_selected_family);
            let family_hits = family_terms
                .iter()
                .filter(|term| haystack.contains(term.as_str()))
                .cloned()
                .collect::<Vec<_>>();
            if !family_hits.is_empty() {
                score += 0.60;
                matched_routine_families.insert(normalized_selected_family.clone());
                matched_terms.extend(family_hits.clone());
                rationale.insert(format!(
                    "{} macro template matches selected routine family '{}'.",
                    macro_kind, normalized_selected_family
                ));
            }
        }

        for family in &context.effective_preferred_routine_families {
            let family_terms = Self::macro_template_family_match_terms(family);
            let family_hits = family_terms
                .iter()
                .filter(|term| haystack.contains(term.as_str()))
                .cloned()
                .collect::<Vec<_>>();
            if family_hits.is_empty() {
                continue;
            }
            score += 0.25;
            matched_routine_families.insert(family.clone());
            matched_terms.extend(family_hits);
            rationale.insert(format!(
                "{} macro template aligns with synthesized preferred routine family '{}'.",
                macro_kind, family
            ));
        }

        if let Some(routine_id) = selected_routine_id {
            let matched_id_terms = Self::routine_macro_id_match_terms(routine_id)
                .into_iter()
                .filter(|term| haystack.contains(term))
                .collect::<Vec<_>>();
            if !matched_id_terms.is_empty() {
                score += 0.15;
                matched_terms.extend(matched_id_terms);
                rationale.insert(format!(
                    "{} macro template references selected routine '{}'.",
                    macro_kind, routine_id
                ));
            }
        }

        if score <= 0.0 {
            return None;
        }

        Some(MacroTemplateSuggestion {
            macro_kind: macro_kind.to_string(),
            template_name: template_name.to_string(),
            description: description.map(|value| value.to_string()),
            details_url: details_url.map(|value| value.to_string()),
            score,
            matched_routine_families: matched_routine_families.into_iter().collect(),
            matched_terms: matched_terms.into_iter().collect(),
            rationale: rationale.into_iter().collect(),
        })
    }

    pub fn planning_routine_preference_context_record(&self) -> RoutinePreferenceContextRecord {
        let objective = self.planning_objective();
        let context = Self::planning_objective_routine_preference_context(&objective);
        Self::routine_preference_context_record(&context)
    }

    fn planning_objective_routine_preference_context_for_sequence(
        &mut self,
        objective: &PlanningObjective,
        seq_id: Option<&str>,
    ) -> RoutinePreferenceContext {
        let mut context = Self::planning_objective_routine_preference_context(objective);
        let graph = seq_id
            .map(str::trim)
            .filter(|value| !value.is_empty())
            .and_then(|seq_id| {
                self.refresh_construct_reasoning_graph_for_seq_id(seq_id)
                    .ok()
            });
        let (
            construct_strategy_seq_id,
            construct_strategy_derived_preferred_routine_families,
            construct_strategy_reasoning_note,
        ) = Self::construct_reasoning_strategy_preference_summary(graph.as_ref());
        let (
            construct_reasoning_seq_id,
            variant_effect_tags,
            variant_suggested_assay_ids,
            variant_derived_preferred_routine_families,
            variant_reasoning_note,
        ) = Self::construct_reasoning_variant_preference_summary(graph.as_ref());
        if construct_reasoning_seq_id.is_some()
            || construct_strategy_seq_id.is_some()
            || !construct_strategy_derived_preferred_routine_families.is_empty()
            || !variant_effect_tags.is_empty()
            || !variant_suggested_assay_ids.is_empty()
            || !variant_derived_preferred_routine_families.is_empty()
        {
            context.construct_reasoning_seq_id =
                construct_reasoning_seq_id.or(construct_strategy_seq_id);
            context.construct_strategy_derived_preferred_routine_families =
                construct_strategy_derived_preferred_routine_families;
            context.variant_effect_tags = variant_effect_tags;
            context.variant_suggested_assay_ids = variant_suggested_assay_ids;
            context.variant_derived_preferred_routine_families =
                variant_derived_preferred_routine_families;
            if let Some(note) = construct_strategy_reasoning_note {
                context.rationale.push(note);
            }
            if let Some(note) = variant_reasoning_note {
                context.rationale.push(note);
            }
            context.rationale.sort();
            context.rationale.dedup();
            context.effective_preferred_routine_families = context
                .explicit_preferred_routine_families
                .iter()
                .chain(context.helper_derived_preferred_routine_families.iter())
                .chain(
                    context
                        .construct_strategy_derived_preferred_routine_families
                        .iter(),
                )
                .chain(context.variant_derived_preferred_routine_families.iter())
                .cloned()
                .collect::<BTreeSet<_>>()
                .into_iter()
                .collect::<Vec<_>>();
            context
                .effective_preferred_routine_families
                .sort_by_key(|value| value.to_ascii_lowercase());
        }
        context
    }

    pub fn planning_routine_preference_context_record_for_sequence(
        &mut self,
        seq_id: Option<&str>,
    ) -> RoutinePreferenceContextRecord {
        let objective = self.planning_objective();
        let context =
            self.planning_objective_routine_preference_context_for_sequence(&objective, seq_id);
        Self::routine_preference_context_record(&context)
    }

    pub(crate) fn planning_routine_preference_context_for_sequence(
        &mut self,
        seq_id: Option<&str>,
    ) -> RoutinePreferenceContext {
        let objective = self.planning_objective();
        self.planning_objective_routine_preference_context_for_sequence(&objective, seq_id)
    }

    pub fn suggest_macro_templates_for_routine(
        &self,
        selected_routine_id: Option<&str>,
        selected_routine_family: Option<&str>,
        limit: usize,
    ) -> Vec<MacroTemplateSuggestion> {
        let objective = self.planning_objective();
        let context = Self::planning_objective_routine_preference_context(&objective);
        Self::suggest_macro_templates_for_routine_with_context(
            self,
            &context,
            selected_routine_id,
            selected_routine_family,
            limit,
        )
    }

    pub fn suggest_macro_templates_for_routine_for_sequence(
        &mut self,
        selected_routine_id: Option<&str>,
        selected_routine_family: Option<&str>,
        seq_id: Option<&str>,
        limit: usize,
    ) -> Vec<MacroTemplateSuggestion> {
        let objective = self.planning_objective();
        let context =
            self.planning_objective_routine_preference_context_for_sequence(&objective, seq_id);
        Self::suggest_macro_templates_for_routine_with_context(
            self,
            &context,
            selected_routine_id,
            selected_routine_family,
            limit,
        )
    }

    fn suggest_macro_templates_for_routine_with_context(
        &self,
        context: &RoutinePreferenceContext,
        selected_routine_id: Option<&str>,
        selected_routine_family: Option<&str>,
        limit: usize,
    ) -> Vec<MacroTemplateSuggestion> {
        let mut suggestions = vec![];

        let workflow_store = self.read_workflow_macro_template_store();
        for template in workflow_store.templates.values() {
            if let Some(suggestion) = Self::score_macro_template_suggestion(
                "workflow",
                &template.name,
                template.description.as_deref(),
                template.details_url.as_deref(),
                &template.script,
                selected_routine_id,
                selected_routine_family,
                context,
            ) {
                suggestions.push(suggestion);
            }
        }

        let candidate_store = self.read_candidate_macro_template_store();
        for template in candidate_store.templates.values() {
            if let Some(suggestion) = Self::score_macro_template_suggestion(
                "candidate",
                &template.name,
                template.description.as_deref(),
                template.details_url.as_deref(),
                &template.script,
                selected_routine_id,
                selected_routine_family,
                context,
            ) {
                suggestions.push(suggestion);
            }
        }

        suggestions.sort_by(|left, right| {
            right
                .score
                .total_cmp(&left.score)
                .then_with(|| left.macro_kind.cmp(&right.macro_kind))
                .then_with(|| left.template_name.cmp(&right.template_name))
        });
        if limit > 0 && suggestions.len() > limit {
            suggestions.truncate(limit);
        }
        suggestions
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

    fn validate_construct_objective_schema(
        objective: &ConstructObjective,
    ) -> Result<(), EngineError> {
        let schema = objective.schema.trim();
        if !schema.is_empty() && schema != CONSTRUCT_OBJECTIVE_SCHEMA {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: format!(
                    "Unsupported construct objective schema '{}' (expected '{}')",
                    schema, CONSTRUCT_OBJECTIVE_SCHEMA
                ),
            });
        }
        Ok(())
    }

    fn normalize_role_list(roles: &mut Vec<ConstructRole>) {
        roles.sort();
        roles.dedup();
    }

    fn normalize_optional_note_text(values: &mut Vec<String>) {
        let mut seen: HashSet<String> = HashSet::new();
        *values = values
            .drain(..)
            .map(|value| value.trim().to_string())
            .filter(|value| !value.is_empty())
            .filter(|value| seen.insert(value.to_ascii_uppercase()))
            .collect();
    }

    fn normalize_optional_id_token_text(value: Option<String>) -> Option<String> {
        value
            .map(|value| value.trim().to_string())
            .filter(|value| !value.is_empty())
            .map(|value| Self::normalize_id_token(&value))
    }

    fn normalize_tag_like_text(values: &mut Vec<String>) {
        let mut seen: HashSet<String> = HashSet::new();
        *values = values
            .drain(..)
            .map(|value| value.trim().to_ascii_lowercase())
            .filter(|value| !value.is_empty())
            .filter(|value| seen.insert(value.clone()))
            .collect();
    }

    fn normalize_host_route_steps(steps: &mut Vec<HostRouteStep>) {
        for (idx, step) in steps.iter_mut().enumerate() {
            step.step_id = if step.step_id.trim().is_empty() {
                format!("host_step_{}", idx.saturating_add(1))
            } else {
                Self::normalize_id_token(step.step_id.trim())
            };
            step.host_profile_id = if step.host_profile_id.trim().is_empty() {
                String::new()
            } else {
                Self::normalize_id_token(step.host_profile_id.trim())
            };
            step.rationale = step.rationale.trim().to_string();
            Self::normalize_optional_note_text(&mut step.notes);
        }
        steps.retain(|step| !step.host_profile_id.is_empty() || !step.rationale.is_empty());
        let mut seen: HashSet<String> = HashSet::new();
        steps.retain(|step| seen.insert(step.step_id.clone()));
    }

    fn normalize_restriction_enzyme_name(raw: &str) -> String {
        let trimmed = raw.trim();
        if trimmed.is_empty() {
            return String::new();
        }
        active_restriction_enzymes()
            .into_iter()
            .find(|enzyme| enzyme.name.eq_ignore_ascii_case(trimmed))
            .map(|enzyme| enzyme.name)
            .unwrap_or_else(|| trimmed.to_string())
    }

    fn normalize_adapter_restriction_capture_plans(plans: &mut Vec<AdapterRestrictionCapturePlan>) {
        for (idx, plan) in plans.iter_mut().enumerate() {
            plan.capture_id = if plan.capture_id.trim().is_empty() {
                format!("adapter_capture_{}", idx.saturating_add(1))
            } else {
                Self::normalize_id_token(plan.capture_id.trim())
            };
            plan.restriction_enzyme_name =
                Self::normalize_restriction_enzyme_name(&plan.restriction_enzyme_name);
            plan.extra_retrieval_enzyme_names = plan
                .extra_retrieval_enzyme_names
                .iter()
                .map(|value| Self::normalize_restriction_enzyme_name(value))
                .filter(|value| !value.is_empty())
                .collect();
            plan.extra_retrieval_enzyme_names
                .sort_by_key(|value| value.to_ascii_lowercase());
            plan.extra_retrieval_enzyme_names.dedup();
            Self::normalize_optional_note_text(&mut plan.notes);
        }
        plans.retain(|plan| {
            !plan.restriction_enzyme_name.is_empty()
                || !plan.extra_retrieval_enzyme_names.is_empty()
                || !plan.notes.is_empty()
        });
        let mut seen: HashSet<String> = HashSet::new();
        plans.retain(|plan| seen.insert(plan.capture_id.clone()));
    }

    fn normalize_construct_objective(mut objective: ConstructObjective) -> ConstructObjective {
        objective.schema = CONSTRUCT_OBJECTIVE_SCHEMA.to_string();
        objective.objective_id = if !objective.objective_id.trim().is_empty() {
            Self::normalize_id_token(objective.objective_id.trim())
        } else {
            let preferred = if !objective.title.trim().is_empty() {
                objective.title.trim()
            } else if !objective.goal.trim().is_empty() {
                objective.goal.trim()
            } else {
                "construct_objective"
            };
            format!(
                "construct_objective_{}",
                Self::normalize_id_token(preferred)
            )
        };
        objective.title = objective.title.trim().to_string();
        objective.goal = objective.goal.trim().to_string();
        objective.host_species = objective
            .host_species
            .take()
            .map(|value| value.trim().to_string())
            .filter(|value| !value.is_empty());
        objective.cell_type = objective
            .cell_type
            .take()
            .map(|value| value.trim().to_string())
            .filter(|value| !value.is_empty());
        objective.tissue = objective
            .tissue
            .take()
            .map(|value| value.trim().to_string())
            .filter(|value| !value.is_empty());
        objective.organelle = objective
            .organelle
            .take()
            .map(|value| value.trim().to_string())
            .filter(|value| !value.is_empty());
        objective.expression_intent = objective
            .expression_intent
            .take()
            .map(|value| value.trim().to_string())
            .filter(|value| !value.is_empty());
        objective.propagation_host_profile_id =
            Self::normalize_optional_id_token_text(objective.propagation_host_profile_id.take());
        objective.expression_host_profile_id =
            Self::normalize_optional_id_token_text(objective.expression_host_profile_id.take());
        objective.helper_profile_id =
            Self::normalize_optional_id_token_text(objective.helper_profile_id.take());
        Self::normalize_host_route_steps(&mut objective.host_route);
        Self::normalize_tag_like_text(&mut objective.medium_conditions);
        Self::normalize_adapter_restriction_capture_plans(
            &mut objective.adapter_restriction_capture_plans,
        );
        Self::normalize_tag_like_text(&mut objective.required_host_traits);
        Self::normalize_tag_like_text(&mut objective.forbidden_host_traits);
        Self::normalize_role_list(&mut objective.required_roles);
        Self::normalize_role_list(&mut objective.forbidden_roles);
        objective.preferred_routine_families = objective
            .preferred_routine_families
            .iter()
            .map(|value| value.trim().to_ascii_lowercase())
            .filter(|value| !value.is_empty())
            .collect();
        objective.preferred_routine_families.sort();
        objective.preferred_routine_families.dedup();
        Self::normalize_optional_note_text(&mut objective.notes);
        objective
    }

    fn normalize_confidence_score(value: Option<f64>) -> Option<f64> {
        value.filter(|value| value.is_finite() && *value >= 0.0 && *value <= 1.0)
    }

    fn normalize_optional_score(value: Option<f64>) -> Option<f64> {
        value.filter(|value| value.is_finite())
    }

    fn normalize_design_evidence(mut evidence: DesignEvidence, idx: usize) -> DesignEvidence {
        evidence.schema = DESIGN_EVIDENCE_SCHEMA.to_string();
        evidence.evidence_id = if evidence.evidence_id.trim().is_empty() {
            match evidence.scope {
                EvidenceScope::SequenceSpan => format!(
                    "evidence_{}_{}_{}_{}",
                    evidence.role.as_str(),
                    evidence.start_0based,
                    evidence.end_0based_exclusive,
                    idx
                ),
                _ => format!(
                    "evidence_{}_{}_{}",
                    evidence.role.as_str(),
                    evidence.scope.as_str(),
                    idx
                ),
            }
        } else {
            evidence.evidence_id.trim().to_string()
        };
        evidence.seq_id = evidence.seq_id.trim().to_string();
        if evidence.end_0based_exclusive < evidence.start_0based {
            std::mem::swap(
                &mut evidence.start_0based,
                &mut evidence.end_0based_exclusive,
            );
        }
        evidence.strand = evidence
            .strand
            .take()
            .map(|value| value.trim().to_string())
            .filter(|value| !value.is_empty());
        evidence.host_profile_id =
            Self::normalize_optional_id_token_text(evidence.host_profile_id.take());
        evidence.host_route_step_id =
            Self::normalize_optional_id_token_text(evidence.host_route_step_id.take());
        evidence.helper_profile_id =
            Self::normalize_optional_id_token_text(evidence.helper_profile_id.take());
        evidence.medium_condition_id =
            Self::normalize_optional_id_token_text(evidence.medium_condition_id.take());
        evidence.label = evidence.label.trim().to_string();
        evidence.rationale = evidence.rationale.trim().to_string();
        evidence.score = Self::normalize_optional_score(evidence.score);
        evidence.confidence = Self::normalize_confidence_score(evidence.confidence);
        evidence.specificity_bias = Self::normalize_confidence_score(evidence.specificity_bias);
        evidence.sensitivity_bias = Self::normalize_confidence_score(evidence.sensitivity_bias);
        evidence.context_tags = evidence
            .context_tags
            .iter()
            .map(|value| value.trim().to_ascii_lowercase())
            .filter(|value| !value.is_empty())
            .collect();
        evidence.context_tags.sort();
        evidence.context_tags.dedup();
        evidence.provenance_kind = evidence.provenance_kind.trim().to_string();
        evidence.provenance_refs = evidence
            .provenance_refs
            .iter()
            .map(|value| value.trim().to_string())
            .filter(|value| !value.is_empty())
            .collect();
        evidence.provenance_refs.sort();
        evidence.provenance_refs.dedup();
        Self::normalize_optional_note_text(&mut evidence.warnings);
        Self::normalize_optional_note_text(&mut evidence.notes);
        evidence
    }

    fn normalize_design_fact(mut fact: DesignFact, idx: usize) -> DesignFact {
        fact.schema = DESIGN_FACT_SCHEMA.to_string();
        fact.fact_id = if fact.fact_id.trim().is_empty() {
            format!("fact_{}_{}", Self::normalize_id_token(&fact.fact_type), idx)
        } else {
            fact.fact_id.trim().to_string()
        };
        fact.fact_type = fact.fact_type.trim().to_string();
        fact.label = fact.label.trim().to_string();
        fact.rationale = fact.rationale.trim().to_string();
        fact.based_on_evidence_ids = fact
            .based_on_evidence_ids
            .iter()
            .map(|value| value.trim().to_string())
            .filter(|value| !value.is_empty())
            .collect();
        fact.based_on_evidence_ids.sort();
        fact.based_on_evidence_ids.dedup();
        fact.confidence = Self::normalize_confidence_score(fact.confidence);
        fact
    }

    fn normalize_design_decision_node(
        mut node: DesignDecisionNode,
        idx: usize,
    ) -> DesignDecisionNode {
        node.schema = DESIGN_DECISION_NODE_SCHEMA.to_string();
        node.decision_id = if node.decision_id.trim().is_empty() {
            format!(
                "decision_{}_{}",
                Self::normalize_id_token(&node.decision_type),
                idx
            )
        } else {
            node.decision_id.trim().to_string()
        };
        node.decision_type = node.decision_type.trim().to_string();
        node.title = node.title.trim().to_string();
        node.rationale = node.rationale.trim().to_string();
        for ids in [
            &mut node.input_evidence_ids,
            &mut node.input_fact_ids,
            &mut node.output_fact_ids,
            &mut node.output_candidate_ids,
        ] {
            *ids = ids
                .iter()
                .map(|value| value.trim().to_string())
                .filter(|value| !value.is_empty())
                .collect();
            ids.sort();
            ids.dedup();
        }
        node
    }

    fn normalize_construct_candidate(
        mut candidate: ConstructCandidate,
        idx: usize,
    ) -> ConstructCandidate {
        candidate.schema = CONSTRUCT_CANDIDATE_SCHEMA.to_string();
        candidate.candidate_id = if candidate.candidate_id.trim().is_empty() {
            format!("candidate_{}", idx)
        } else {
            candidate.candidate_id.trim().to_string()
        };
        candidate.objective_id = candidate.objective_id.trim().to_string();
        candidate.title = candidate.title.trim().to_string();
        for ids in [
            &mut candidate.component_ids,
            &mut candidate.derived_from_fact_ids,
            &mut candidate.suggested_routine_ids,
        ] {
            *ids = ids
                .iter()
                .map(|value| value.trim().to_string())
                .filter(|value| !value.is_empty())
                .collect();
            ids.sort();
            ids.dedup();
        }
        candidate.compactness_score = Self::normalize_optional_score(candidate.compactness_score);
        candidate.confidence_score = Self::normalize_confidence_score(candidate.confidence_score);
        candidate.cloning_complexity_score =
            Self::normalize_optional_score(candidate.cloning_complexity_score);
        candidate.host_fit_score = Self::normalize_confidence_score(candidate.host_fit_score);
        Self::normalize_optional_note_text(&mut candidate.warnings);
        Self::normalize_optional_note_text(&mut candidate.notes);
        candidate
    }

    fn normalize_annotation_candidate(
        mut candidate: AnnotationCandidate,
        idx: usize,
    ) -> AnnotationCandidate {
        candidate.schema = ANNOTATION_CANDIDATE_SCHEMA.to_string();
        candidate.annotation_id = if candidate.annotation_id.trim().is_empty() {
            let evidence_token = if candidate.evidence_id.trim().is_empty() {
                idx.to_string()
            } else {
                Self::normalize_id_token(&candidate.evidence_id)
            };
            format!(
                "annotation_{}_{}",
                Self::normalize_id_token(candidate.role.as_str()),
                evidence_token
            )
        } else {
            candidate.annotation_id.trim().to_string()
        };
        candidate.evidence_id = candidate.evidence_id.trim().to_string();
        candidate.seq_id = candidate.seq_id.trim().to_string();
        if candidate.end_0based_exclusive < candidate.start_0based {
            std::mem::swap(
                &mut candidate.start_0based,
                &mut candidate.end_0based_exclusive,
            );
        }
        candidate.strand = candidate
            .strand
            .take()
            .map(|value| value.trim().to_string())
            .filter(|value| !value.is_empty());
        candidate.label = candidate.label.trim().to_string();
        candidate.rationale = candidate.rationale.trim().to_string();
        candidate.source_kind = candidate.source_kind.trim().to_string();
        for values in [
            &mut candidate.supporting_fact_ids,
            &mut candidate.supporting_fact_labels,
            &mut candidate.supporting_decision_ids,
            &mut candidate.supporting_decision_titles,
            &mut candidate.effect_tags,
        ] {
            *values = values
                .iter()
                .map(|value| value.trim().to_string())
                .filter(|value| !value.is_empty())
                .collect();
            values.sort();
            values.dedup();
        }
        candidate.transcript_context_status = candidate
            .transcript_context_status
            .take()
            .map(|value| value.trim().to_string())
            .filter(|value| !value.is_empty());
        Self::normalize_optional_note_text(&mut candidate.warnings);
        Self::normalize_optional_note_text(&mut candidate.notes);
        candidate
    }

    fn construct_reasoning_annotation_candidate_summary_group_key(
        candidate: &AnnotationCandidate,
    ) -> (
        String,
        ConstructRole,
        Option<String>,
        String,
        String,
        String,
    ) {
        (
            candidate.seq_id.clone(),
            candidate.role,
            candidate.strand.clone(),
            candidate.source_kind.clone(),
            candidate.supporting_fact_ids.join("\u{1f}"),
            candidate.supporting_decision_ids.join("\u{1f}"),
        )
    }

    fn construct_reasoning_role_display_label(role: ConstructRole) -> &'static str {
        match role {
            ConstructRole::Promoter => "Promoter",
            ConstructRole::Enhancer => "Enhancer",
            ConstructRole::Gene => "Gene",
            ConstructRole::Transcript => "Transcript",
            ConstructRole::Exon => "Exon",
            ConstructRole::Utr5Prime => "5' UTR",
            ConstructRole::Cds => "CDS",
            ConstructRole::Utr3Prime => "3' UTR",
            ConstructRole::Terminator => "Terminator",
            ConstructRole::Variant => "Variant",
            ConstructRole::Linker => "Linker",
            ConstructRole::Tag => "Tag",
            ConstructRole::SignalPeptide => "Signal peptide",
            ConstructRole::LocalizationSignal => "Localization signal",
            ConstructRole::HomologyArm => "Homology arm",
            ConstructRole::FusionBoundary => "Fusion boundary",
            ConstructRole::RepeatRegion => "Repeat/similarity",
            ConstructRole::MobileElement => "Mobile element",
            ConstructRole::RestrictionSite => "Restriction site",
            ConstructRole::SpliceBoundary => "Splice boundary",
            ConstructRole::Tfbs => "TFBS",
            ConstructRole::ContextBaggage => "Context baggage",
            ConstructRole::Other => "Other",
        }
    }

    fn construct_reasoning_annotation_summary_review_status(
        candidates: &[&AnnotationCandidate],
    ) -> String {
        let mut statuses = candidates
            .iter()
            .map(|candidate| candidate.editable_status.as_str().to_string())
            .collect::<Vec<_>>();
        statuses.sort();
        statuses.dedup();
        match statuses.len() {
            0 => "draft".to_string(),
            1 => statuses[0].clone(),
            _ => format!("mixed: {}", statuses.join(", ")),
        }
    }

    fn construct_reasoning_annotation_summary_title(
        role: ConstructRole,
        candidates: &[&AnnotationCandidate],
    ) -> String {
        let mut labels = candidates
            .iter()
            .filter_map(|candidate| {
                let value = candidate.label.trim();
                (!value.is_empty()).then(|| value.to_string())
            })
            .collect::<Vec<_>>();
        labels.sort();
        labels.dedup();
        let count = candidates.len();
        if labels.len() == 1 {
            if count > 1 {
                format!("{} ({count})", labels[0])
            } else {
                labels[0].clone()
            }
        } else if count > 1 {
            format!(
                "{} ({count})",
                Self::construct_reasoning_role_display_label(role)
            )
        } else {
            Self::construct_reasoning_role_display_label(role).to_string()
        }
    }

    fn construct_reasoning_annotation_summary_subtitle(
        role: ConstructRole,
        candidate_count: usize,
        transcript_context_statuses: &[String],
    ) -> String {
        if transcript_context_statuses
            .iter()
            .any(|status| status == "multi_transcript_ambiguous")
        {
            return "Transcript interpretations disagree".to_string();
        }
        if candidate_count > 1 {
            return match role {
                ConstructRole::Promoter => {
                    format!("{candidate_count} overlapping promoter-linked candidates")
                }
                ConstructRole::Cds => format!("{candidate_count} overlapping coding candidates"),
                ConstructRole::Exon => format!("{candidate_count} overlapping exon candidates"),
                ConstructRole::RepeatRegion => {
                    format!("{candidate_count} overlapping repeat/similarity candidates")
                }
                ConstructRole::MobileElement => {
                    format!("{candidate_count} overlapping mobile-element candidates")
                }
                ConstructRole::SpliceBoundary => {
                    format!("{candidate_count} overlapping splice-relevant candidates")
                }
                ConstructRole::Variant => {
                    format!("{candidate_count} overlapping variant-linked candidates")
                }
                _ => format!("{candidate_count} overlapping candidates collapsed"),
            };
        }
        match role {
            ConstructRole::Promoter => "Promoter candidate".to_string(),
            ConstructRole::Enhancer => "Enhancer candidate".to_string(),
            ConstructRole::Cds => "Coding candidate".to_string(),
            ConstructRole::Exon => "Exon candidate".to_string(),
            ConstructRole::RepeatRegion => "Repeat/similarity candidate".to_string(),
            ConstructRole::MobileElement => "Mobile-element candidate".to_string(),
            ConstructRole::SpliceBoundary => "Splice-relevant candidate".to_string(),
            ConstructRole::Variant => "Variant-linked candidate".to_string(),
            ConstructRole::Transcript => "Transcript-linked candidate".to_string(),
            _ => format!(
                "{} candidate",
                Self::construct_reasoning_role_display_label(role).to_ascii_lowercase()
            ),
        }
    }

    fn construct_reasoning_build_annotation_candidate_summary(
        seq_id: &str,
        role: ConstructRole,
        strand: Option<&str>,
        candidates: &[&AnnotationCandidate],
    ) -> AnnotationCandidateSummary {
        let start_0based = candidates
            .iter()
            .map(|candidate| candidate.start_0based)
            .min()
            .unwrap_or_default();
        let end_0based_exclusive = candidates
            .iter()
            .map(|candidate| candidate.end_0based_exclusive)
            .max()
            .unwrap_or(start_0based);
        let mut annotation_ids = candidates
            .iter()
            .map(|candidate| candidate.annotation_id.clone())
            .collect::<Vec<_>>();
        annotation_ids.sort();
        annotation_ids.dedup();
        let mut source_kinds = candidates
            .iter()
            .filter_map(|candidate| {
                let value = candidate.source_kind.trim();
                (!value.is_empty()).then(|| value.to_string())
            })
            .collect::<Vec<_>>();
        source_kinds.sort();
        source_kinds.dedup();
        let mut transcript_context_statuses = candidates
            .iter()
            .filter_map(|candidate| candidate.transcript_context_status.clone())
            .collect::<Vec<_>>();
        transcript_context_statuses.sort();
        transcript_context_statuses.dedup();
        let mut effect_tags = candidates
            .iter()
            .flat_map(|candidate| candidate.effect_tags.iter().cloned())
            .collect::<Vec<_>>();
        effect_tags.sort();
        effect_tags.dedup();
        let mut supporting_fact_labels = candidates
            .iter()
            .flat_map(|candidate| candidate.supporting_fact_labels.iter().cloned())
            .collect::<Vec<_>>();
        supporting_fact_labels.sort();
        supporting_fact_labels.dedup();
        let mut supporting_decision_titles = candidates
            .iter()
            .flat_map(|candidate| candidate.supporting_decision_titles.iter().cloned())
            .collect::<Vec<_>>();
        supporting_decision_titles.sort();
        supporting_decision_titles.dedup();
        let mut warnings = candidates
            .iter()
            .flat_map(|candidate| candidate.warnings.iter().cloned())
            .collect::<Vec<_>>();
        warnings.sort();
        warnings.dedup();
        let mut notes = candidates
            .iter()
            .flat_map(|candidate| candidate.notes.iter().cloned())
            .collect::<Vec<_>>();
        notes.sort();
        notes.dedup();
        AnnotationCandidateSummary {
            summary_id: format!(
                "summary_{}_{}_{}_{}",
                Self::normalize_id_token(role.as_str()),
                start_0based,
                end_0based_exclusive,
                annotation_ids
                    .first()
                    .map(|value| Self::normalize_id_token(value))
                    .unwrap_or_else(|| "annotation".to_string())
            ),
            seq_id: seq_id.to_string(),
            start_0based,
            end_0based_exclusive,
            strand: strand.map(str::to_string),
            role,
            title: Self::construct_reasoning_annotation_summary_title(role, candidates),
            subtitle: Self::construct_reasoning_annotation_summary_subtitle(
                role,
                candidates.len(),
                &transcript_context_statuses,
            ),
            annotation_ids,
            source_kinds,
            transcript_context_statuses,
            effect_tags,
            candidate_count: candidates.len(),
            review_status_summary: Self::construct_reasoning_annotation_summary_review_status(
                candidates,
            ),
            supporting_fact_labels,
            supporting_decision_titles,
            warnings,
            notes,
            ..AnnotationCandidateSummary::default()
        }
    }

    fn construct_reasoning_build_annotation_candidate_summaries(
        annotation_candidates: &[AnnotationCandidate],
    ) -> Vec<AnnotationCandidateSummary> {
        let mut by_signature: BTreeMap<
            (
                String,
                ConstructRole,
                Option<String>,
                String,
                String,
                String,
            ),
            Vec<&AnnotationCandidate>,
        > = BTreeMap::new();
        for candidate in annotation_candidates {
            by_signature
                .entry(Self::construct_reasoning_annotation_candidate_summary_group_key(candidate))
                .or_default()
                .push(candidate);
        }

        let mut summaries = Vec::new();
        for ((seq_id, role, strand, _, _, _), mut candidates) in by_signature {
            candidates.sort_by(|left, right| {
                left.start_0based
                    .cmp(&right.start_0based)
                    .then(left.end_0based_exclusive.cmp(&right.end_0based_exclusive))
                    .then(left.annotation_id.cmp(&right.annotation_id))
            });
            let mut cluster: Vec<&AnnotationCandidate> = Vec::new();
            let mut cluster_end = 0usize;
            for candidate in candidates {
                if cluster.is_empty() {
                    cluster_end = candidate.end_0based_exclusive;
                    cluster.push(candidate);
                } else if candidate.start_0based <= cluster_end {
                    cluster_end = cluster_end.max(candidate.end_0based_exclusive);
                    cluster.push(candidate);
                } else {
                    summaries.push(
                        Self::construct_reasoning_build_annotation_candidate_summary(
                            &seq_id,
                            role,
                            strand.as_deref(),
                            &cluster,
                        ),
                    );
                    cluster = vec![candidate];
                    cluster_end = candidate.end_0based_exclusive;
                }
            }
            if !cluster.is_empty() {
                summaries.push(
                    Self::construct_reasoning_build_annotation_candidate_summary(
                        &seq_id,
                        role,
                        strand.as_deref(),
                        &cluster,
                    ),
                );
            }
        }
        summaries
    }

    fn normalize_annotation_candidate_summary(
        mut summary: AnnotationCandidateSummary,
        idx: usize,
    ) -> AnnotationCandidateSummary {
        summary.schema = ANNOTATION_CANDIDATE_SUMMARY_SCHEMA.to_string();
        summary.summary_id = if summary.summary_id.trim().is_empty() {
            format!(
                "summary_{}_{}_{}_{}",
                Self::normalize_id_token(summary.role.as_str()),
                summary.start_0based,
                summary.end_0based_exclusive,
                idx
            )
        } else {
            summary.summary_id.trim().to_string()
        };
        summary.seq_id = summary.seq_id.trim().to_string();
        if summary.end_0based_exclusive < summary.start_0based {
            std::mem::swap(&mut summary.start_0based, &mut summary.end_0based_exclusive);
        }
        summary.strand = summary
            .strand
            .take()
            .map(|value| value.trim().to_string())
            .filter(|value| !value.is_empty());
        summary.title = summary.title.trim().to_string();
        summary.subtitle = summary.subtitle.trim().to_string();
        summary.review_status_summary = summary.review_status_summary.trim().to_string();
        for values in [
            &mut summary.annotation_ids,
            &mut summary.source_kinds,
            &mut summary.transcript_context_statuses,
            &mut summary.effect_tags,
            &mut summary.supporting_fact_labels,
            &mut summary.supporting_decision_titles,
        ] {
            *values = values
                .iter()
                .map(|value| value.trim().to_string())
                .filter(|value| !value.is_empty())
                .collect();
            values.sort();
            values.dedup();
        }
        summary.candidate_count = summary.candidate_count.max(summary.annotation_ids.len());
        Self::normalize_optional_note_text(&mut summary.warnings);
        Self::normalize_optional_note_text(&mut summary.notes);
        summary
    }

    fn construct_reasoning_annotation_candidate_match_key(
        candidate: &AnnotationCandidate,
    ) -> (String, usize, usize, ConstructRole) {
        (
            candidate.evidence_id.clone(),
            candidate.start_0based,
            candidate.end_0based_exclusive,
            candidate.role,
        )
    }

    fn construct_reasoning_annotation_candidate_span_key(
        candidate: &AnnotationCandidate,
    ) -> (usize, usize, ConstructRole, String) {
        (
            candidate.start_0based,
            candidate.end_0based_exclusive,
            candidate.role,
            candidate.label.trim().to_ascii_lowercase(),
        )
    }

    fn preserve_construct_reasoning_annotation_candidate_statuses(
        annotation_candidates: &mut [AnnotationCandidate],
        previous_graph: Option<&ConstructReasoningGraph>,
    ) {
        let Some(previous_graph) = previous_graph else {
            return;
        };
        let by_annotation_id = previous_graph
            .annotation_candidates
            .iter()
            .map(|candidate| (candidate.annotation_id.as_str(), candidate.editable_status))
            .collect::<HashMap<_, _>>();
        let by_match_key = previous_graph
            .annotation_candidates
            .iter()
            .map(|candidate| {
                (
                    Self::construct_reasoning_annotation_candidate_match_key(candidate),
                    candidate.editable_status,
                )
            })
            .collect::<BTreeMap<_, _>>();
        let by_span_key = previous_graph
            .annotation_candidates
            .iter()
            .map(|candidate| {
                (
                    Self::construct_reasoning_annotation_candidate_span_key(candidate),
                    candidate.editable_status,
                )
            })
            .collect::<BTreeMap<_, _>>();
        for candidate in annotation_candidates.iter_mut() {
            let preserved = by_annotation_id
                .get(candidate.annotation_id.as_str())
                .copied()
                .or_else(|| {
                    by_match_key
                        .get(&Self::construct_reasoning_annotation_candidate_match_key(
                            candidate,
                        ))
                        .copied()
                })
                .or_else(|| {
                    by_span_key
                        .get(&Self::construct_reasoning_annotation_candidate_span_key(
                            candidate,
                        ))
                        .copied()
                });
            if let Some(status) = preserved {
                candidate.editable_status = status;
            }
        }
    }

    fn normalize_construct_reasoning_graph(
        mut graph: ConstructReasoningGraph,
    ) -> ConstructReasoningGraph {
        graph.schema = CONSTRUCT_REASONING_GRAPH_SCHEMA.to_string();
        graph.seq_id = graph.seq_id.trim().to_string();
        graph.op_id = graph
            .op_id
            .take()
            .map(|value| value.trim().to_string())
            .filter(|value| !value.is_empty());
        graph.run_id = graph
            .run_id
            .take()
            .map(|value| value.trim().to_string())
            .filter(|value| !value.is_empty());
        graph.objective = Self::normalize_construct_objective(graph.objective);
        graph.graph_id = if graph.graph_id.trim().is_empty() {
            format!(
                "construct_reasoning_{}_{}",
                Self::normalize_id_token(&graph.seq_id),
                Self::normalize_id_token(&graph.objective.objective_id)
            )
        } else {
            graph.graph_id.trim().to_string()
        };
        if graph.generated_at_unix_ms == 0 {
            graph.generated_at_unix_ms = Self::now_unix_ms();
        }
        for (idx, evidence) in graph.evidence.iter_mut().enumerate() {
            let mut normalized = Self::normalize_design_evidence(evidence.clone(), idx);
            if normalized.seq_id.is_empty() {
                normalized.seq_id = graph.seq_id.clone();
            }
            *evidence = normalized;
        }
        graph.evidence.sort_by(|a, b| {
            a.start_0based
                .cmp(&b.start_0based)
                .then(a.end_0based_exclusive.cmp(&b.end_0based_exclusive))
                .then(a.role.cmp(&b.role))
                .then(a.evidence_class.cmp(&b.evidence_class))
                .then(a.label.cmp(&b.label))
                .then(a.evidence_id.cmp(&b.evidence_id))
        });
        graph.facts = graph
            .facts
            .into_iter()
            .enumerate()
            .map(|(idx, fact)| Self::normalize_design_fact(fact, idx))
            .collect();
        graph.facts.sort_by(|a, b| a.fact_id.cmp(&b.fact_id));
        graph.decisions = graph
            .decisions
            .into_iter()
            .enumerate()
            .map(|(idx, node)| Self::normalize_design_decision_node(node, idx))
            .collect();
        graph
            .decisions
            .sort_by(|a, b| a.decision_id.cmp(&b.decision_id));
        graph.candidates = graph
            .candidates
            .into_iter()
            .enumerate()
            .map(|(idx, candidate)| Self::normalize_construct_candidate(candidate, idx))
            .collect();
        graph
            .candidates
            .sort_by(|a, b| a.candidate_id.cmp(&b.candidate_id));
        if graph.annotation_candidates.is_empty() {
            graph.annotation_candidates = Self::construct_reasoning_build_annotation_candidates(
                &graph.seq_id,
                &graph.evidence,
                &graph.facts,
                &graph.decisions,
            );
        }
        graph.annotation_candidates = graph
            .annotation_candidates
            .into_iter()
            .enumerate()
            .map(|(idx, candidate)| Self::normalize_annotation_candidate(candidate, idx))
            .collect();
        graph.annotation_candidates.sort_by(|a, b| {
            a.start_0based
                .cmp(&b.start_0based)
                .then(a.end_0based_exclusive.cmp(&b.end_0based_exclusive))
                .then(a.role.cmp(&b.role))
                .then(a.label.cmp(&b.label))
                .then(a.annotation_id.cmp(&b.annotation_id))
        });
        if graph.annotation_candidate_summaries.is_empty() {
            graph.annotation_candidate_summaries =
                Self::construct_reasoning_build_annotation_candidate_summaries(
                    &graph.annotation_candidates,
                );
        }
        graph.annotation_candidate_summaries = graph
            .annotation_candidate_summaries
            .into_iter()
            .enumerate()
            .map(|(idx, summary)| Self::normalize_annotation_candidate_summary(summary, idx))
            .collect();
        graph.annotation_candidate_summaries.sort_by(|a, b| {
            a.start_0based
                .cmp(&b.start_0based)
                .then(a.end_0based_exclusive.cmp(&b.end_0based_exclusive))
                .then(a.role.cmp(&b.role))
                .then(a.title.cmp(&b.title))
                .then(a.summary_id.cmp(&b.summary_id))
        });
        Self::normalize_optional_note_text(&mut graph.notes);
        graph
    }

    fn read_construct_reasoning_store_from_metadata(
        value: Option<&serde_json::Value>,
    ) -> ConstructReasoningStore {
        let mut store = value
            .cloned()
            .and_then(|v| serde_json::from_value::<ConstructReasoningStore>(v).ok())
            .unwrap_or_default();
        if store.schema.trim().is_empty() {
            store.schema = CONSTRUCT_REASONING_STORE_SCHEMA.to_string();
        }
        store.objectives = store
            .objectives
            .into_values()
            .map(Self::normalize_construct_objective)
            .map(|objective| (objective.objective_id.clone(), objective))
            .collect();
        store.graphs = store
            .graphs
            .into_values()
            .map(Self::normalize_construct_reasoning_graph)
            .map(|graph| (graph.graph_id.clone(), graph))
            .collect();
        store.preferred_graph_id = store
            .preferred_graph_id
            .take()
            .map(|value| value.trim().to_string())
            .filter(|value| !value.is_empty() && store.graphs.contains_key(value));
        store
    }

    fn read_construct_reasoning_store(&self) -> ConstructReasoningStore {
        Self::read_construct_reasoning_store_from_metadata(
            self.state.metadata.get(CONSTRUCT_REASONING_METADATA_KEY),
        )
    }

    fn write_construct_reasoning_store(
        &mut self,
        mut store: ConstructReasoningStore,
    ) -> Result<(), EngineError> {
        store.schema = CONSTRUCT_REASONING_STORE_SCHEMA.to_string();
        store.updated_at_unix_ms = Self::now_unix_ms();
        store.objectives = store
            .objectives
            .into_values()
            .map(Self::normalize_construct_objective)
            .map(|objective| (objective.objective_id.clone(), objective))
            .collect();
        store.graphs = store
            .graphs
            .into_values()
            .map(Self::normalize_construct_reasoning_graph)
            .map(|graph| (graph.graph_id.clone(), graph))
            .collect();
        if let Some(preferred_graph_id) = store.preferred_graph_id.clone() {
            if preferred_graph_id.trim().is_empty()
                || !store.graphs.contains_key(&preferred_graph_id)
            {
                store.preferred_graph_id = None;
            }
        }
        if store.objectives.is_empty()
            && store.graphs.is_empty()
            && store.preferred_graph_id.is_none()
        {
            self.state.metadata.remove(CONSTRUCT_REASONING_METADATA_KEY);
            return Ok(());
        }
        let value = serde_json::to_value(store).map_err(|e| EngineError {
            code: ErrorCode::Internal,
            message: format!("Could not serialize construct-reasoning metadata: {e}"),
        })?;
        self.state
            .metadata
            .insert(CONSTRUCT_REASONING_METADATA_KEY.to_string(), value);
        Ok(())
    }

    pub fn construct_reasoning_store(&self) -> ConstructReasoningStore {
        self.read_construct_reasoning_store()
    }

    pub fn list_construct_reasoning_graph_summaries(
        &self,
        seq_id_filter: Option<&str>,
    ) -> Vec<ConstructReasoningGraphSummary> {
        let filter = seq_id_filter
            .map(str::trim)
            .filter(|value| !value.is_empty())
            .map(|value| value.to_ascii_lowercase());
        let mut rows =
            self.read_construct_reasoning_store()
                .graphs
                .values()
                .filter(|graph| {
                    filter
                        .as_ref()
                        .is_none_or(|needle| graph.seq_id.to_ascii_lowercase() == *needle)
                })
                .map(|graph| {
                    let summary =
                        Self::summarize_process_run_bundle_construct_reasoning_graph(graph);
                    let (
                        protein_to_dna_handoff_candidate_count,
                        protein_to_dna_source_protein_seq_ids,
                    ) = Self::construct_reasoning_graph_protein_to_dna_handoff_summary(graph);
                    ConstructReasoningGraphSummary {
                        graph_id: graph.graph_id.clone(),
                        seq_id: graph.seq_id.clone(),
                        generated_at_unix_ms: graph.generated_at_unix_ms,
                        op_id: graph.op_id.clone(),
                        run_id: graph.run_id.clone(),
                        objective_id: graph.objective.objective_id.clone(),
                        objective_title: graph.objective.title.clone(),
                        objective_goal: graph.objective.goal.clone(),
                        evidence_count: graph.evidence.len(),
                        decision_count: graph.decisions.len(),
                        candidate_count: graph.candidates.len(),
                        contains_protein_to_dna_handoff: protein_to_dna_handoff_candidate_count > 0,
                        protein_to_dna_handoff_candidate_count,
                        protein_to_dna_source_protein_seq_ids,
                        summary_lines: summary.summary_lines,
                        warning_lines: summary.warning_lines,
                    }
                })
                .collect::<Vec<_>>();
        rows.sort_by(|left, right| {
            left.seq_id
                .to_ascii_lowercase()
                .cmp(&right.seq_id.to_ascii_lowercase())
                .then(left.generated_at_unix_ms.cmp(&right.generated_at_unix_ms))
                .then(
                    left.graph_id
                        .to_ascii_lowercase()
                        .cmp(&right.graph_id.to_ascii_lowercase()),
                )
        });
        rows
    }

    pub fn upsert_construct_objective(
        &mut self,
        objective: ConstructObjective,
    ) -> Result<ConstructObjective, EngineError> {
        Self::validate_construct_objective_schema(&objective)?;
        let normalized = Self::normalize_construct_objective(objective);
        let mut store = self.read_construct_reasoning_store();
        store
            .objectives
            .insert(normalized.objective_id.clone(), normalized.clone());
        self.write_construct_reasoning_store(store)?;
        Ok(normalized)
    }

    pub fn construct_reasoning_graph(
        &self,
        graph_id: &str,
    ) -> Result<ConstructReasoningGraph, EngineError> {
        let target = graph_id.trim();
        if target.is_empty() {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: "graph_id cannot be empty".to_string(),
            });
        }
        self.read_construct_reasoning_store()
            .graphs
            .get(target)
            .cloned()
            .ok_or_else(|| EngineError {
                code: ErrorCode::NotFound,
                message: format!("Construct reasoning graph '{}' not found", target),
            })
    }

    fn select_construct_reasoning_graph_for_seq_id(
        store: &ConstructReasoningStore,
        seq_id: &str,
    ) -> Option<ConstructReasoningGraph> {
        if let Some(graph) = store
            .preferred_graph_id
            .as_deref()
            .and_then(|graph_id| store.graphs.get(graph_id))
            .filter(|graph| graph.seq_id == seq_id)
        {
            return Some(graph.clone());
        }
        store
            .graphs
            .values()
            .filter(|graph| graph.seq_id == seq_id)
            .max_by(|left, right| {
                left.generated_at_unix_ms
                    .cmp(&right.generated_at_unix_ms)
                    .then_with(|| left.graph_id.cmp(&right.graph_id))
            })
            .cloned()
    }

    pub fn construct_reasoning_graph_for_seq_id(
        &self,
        seq_id: &str,
    ) -> Result<ConstructReasoningGraph, EngineError> {
        let target = seq_id.trim();
        if target.is_empty() {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: "seq_id cannot be empty".to_string(),
            });
        }
        Self::select_construct_reasoning_graph_for_seq_id(
            &self.read_construct_reasoning_store(),
            target,
        )
        .ok_or_else(|| EngineError {
            code: ErrorCode::NotFound,
            message: format!(
                "No construct reasoning graph stored for sequence '{}'",
                target
            ),
        })
    }

    pub fn refresh_construct_reasoning_graph_for_seq_id(
        &mut self,
        seq_id: &str,
    ) -> Result<ConstructReasoningGraph, EngineError> {
        let target = seq_id.trim();
        if target.is_empty() {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: "seq_id cannot be empty".to_string(),
            });
        }
        let existing = Self::select_construct_reasoning_graph_for_seq_id(
            &self.read_construct_reasoning_store(),
            target,
        );
        let objective_id = existing
            .as_ref()
            .map(|graph| graph.objective.objective_id.as_str());
        let graph_id = existing.as_ref().map(|graph| graph.graph_id.as_str());
        self.build_construct_reasoning_graph(target, objective_id, graph_id)
    }

    pub fn upsert_construct_reasoning_graph(
        &mut self,
        graph: ConstructReasoningGraph,
    ) -> Result<ConstructReasoningGraph, EngineError> {
        let normalized = Self::normalize_construct_reasoning_graph(graph);
        let mut store = self.read_construct_reasoning_store();
        store.objectives.insert(
            normalized.objective.objective_id.clone(),
            normalized.objective.clone(),
        );
        store
            .graphs
            .insert(normalized.graph_id.clone(), normalized.clone());
        if store.preferred_graph_id.is_none() {
            store.preferred_graph_id = Some(normalized.graph_id.clone());
        }
        self.write_construct_reasoning_store(store)?;
        Ok(normalized)
    }

    pub fn set_construct_reasoning_annotation_candidate_status(
        &mut self,
        graph_id: &str,
        annotation_id: &str,
        editable_status: EditableStatus,
    ) -> Result<ConstructReasoningGraph, EngineError> {
        let graph_id = graph_id.trim();
        if graph_id.is_empty() {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: "graph_id cannot be empty".to_string(),
            });
        }
        let annotation_id = annotation_id.trim();
        if annotation_id.is_empty() {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: "annotation_id cannot be empty".to_string(),
            });
        }
        let mut store = self.read_construct_reasoning_store();
        let Some(mut graph) = store.graphs.get(graph_id).cloned() else {
            return Err(EngineError {
                code: ErrorCode::NotFound,
                message: format!("Construct reasoning graph '{}' not found", graph_id),
            });
        };
        let Some(candidate) = graph
            .annotation_candidates
            .iter_mut()
            .find(|candidate| candidate.annotation_id == annotation_id)
        else {
            return Err(EngineError {
                code: ErrorCode::NotFound,
                message: format!(
                    "Annotation candidate '{}' not found in construct reasoning graph '{}'",
                    annotation_id, graph_id
                ),
            });
        };
        if candidate.editable_status == EditableStatus::Locked
            && editable_status != EditableStatus::Locked
        {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: format!(
                    "Annotation candidate '{}' in graph '{}' is locked",
                    annotation_id, graph_id
                ),
            });
        }
        candidate.editable_status = editable_status;
        graph.generated_at_unix_ms = Self::now_unix_ms();
        let graph = Self::normalize_construct_reasoning_graph(graph);
        store.graphs.insert(graph.graph_id.clone(), graph.clone());
        self.write_construct_reasoning_store(store)?;
        Ok(graph)
    }

    fn construct_reasoning_annotation_candidate_is_writeback_eligible(
        candidate: &AnnotationCandidate,
        evidence: &DesignEvidence,
    ) -> bool {
        candidate.source_kind == "generated_annotation"
            && evidence.scope == EvidenceScope::SequenceSpan
    }

    fn construct_reasoning_annotation_candidate_feature_kind(role: ConstructRole) -> &'static str {
        match role {
            ConstructRole::Gene => "gene",
            ConstructRole::Transcript => "mRNA",
            ConstructRole::Exon => "exon",
            ConstructRole::Cds => "CDS",
            ConstructRole::Utr5Prime => "5'UTR",
            ConstructRole::Utr3Prime => "3'UTR",
            ConstructRole::Promoter => "promoter",
            ConstructRole::Enhancer => "enhancer",
            ConstructRole::Terminator => "terminator",
            ConstructRole::Variant => "variation",
            ConstructRole::Linker => "misc_feature",
            ConstructRole::Tag => "misc_feature",
            ConstructRole::SignalPeptide => "sig_peptide",
            ConstructRole::LocalizationSignal => "misc_feature",
            ConstructRole::HomologyArm => "misc_feature",
            ConstructRole::FusionBoundary => "misc_feature",
            ConstructRole::RepeatRegion => "repeat_region",
            ConstructRole::MobileElement => "mobile_element",
            ConstructRole::Tfbs | ConstructRole::SpliceBoundary => "regulatory",
            ConstructRole::RestrictionSite
            | ConstructRole::ContextBaggage
            | ConstructRole::Other => "misc_feature",
        }
    }

    fn construct_reasoning_annotation_candidate_feature_location(
        candidate: &AnnotationCandidate,
    ) -> gb_io::seq::Location {
        let base = gb_io::seq::Location::simple_range(
            candidate.start_0based as i64,
            candidate.end_0based_exclusive as i64,
        );
        if candidate.strand.as_deref() == Some("-") {
            gb_io::seq::Location::Complement(Box::new(base))
        } else {
            base
        }
    }

    fn construct_reasoning_feature_matches_annotation_candidate(
        feature: &gb_io::seq::Feature,
        annotation_id: &str,
    ) -> bool {
        Self::feature_has_qualifier_value(
            feature,
            "construct_reasoning_annotation_id",
            annotation_id,
        )
    }

    fn build_construct_reasoning_annotation_candidate_feature(
        graph_id: &str,
        candidate: &AnnotationCandidate,
        evidence: &DesignEvidence,
    ) -> gb_io::seq::Feature {
        let label = candidate
            .label
            .trim()
            .to_string()
            .chars()
            .collect::<String>();
        let feature_kind =
            Self::construct_reasoning_annotation_candidate_feature_kind(candidate.role);
        let display_label = if label.trim().is_empty() {
            candidate.role.as_str().to_string()
        } else {
            label
        };
        let mut note = format!(
            "Written back from construct reasoning graph '{}' candidate '{}' (source_kind='{}').",
            graph_id, candidate.annotation_id, candidate.source_kind
        );
        if !candidate.rationale.trim().is_empty() {
            note.push(' ');
            note.push_str(candidate.rationale.trim());
        }
        let mut qualifiers = vec![
            ("label".into(), Some(display_label)),
            ("note".into(), Some(note)),
            (
                "gentle_generated".into(),
                Some(CONSTRUCT_REASONING_ANNOTATION_WRITEBACK_GENERATED_TAG.to_string()),
            ),
            (
                "generated_by".into(),
                Some(CONSTRUCT_REASONING_ANNOTATION_WRITEBACK_GENERATED_BY.to_string()),
            ),
            (
                "construct_reasoning_graph_id".into(),
                Some(graph_id.to_string()),
            ),
            (
                "construct_reasoning_annotation_id".into(),
                Some(candidate.annotation_id.clone()),
            ),
            (
                "construct_reasoning_evidence_id".into(),
                Some(candidate.evidence_id.clone()),
            ),
            (
                "construct_reasoning_source_kind".into(),
                Some(candidate.source_kind.clone()),
            ),
            (
                "construct_reasoning_editable_status".into(),
                Some(candidate.editable_status.as_str().to_string()),
            ),
            (
                "construct_reasoning_role".into(),
                Some(candidate.role.as_str().to_string()),
            ),
            (
                "construct_reasoning_evidence_provenance_kind".into(),
                Some(evidence.provenance_kind.clone()),
            ),
        ];
        if let Some(strand) = candidate.strand.as_deref() {
            qualifiers.push(("strand".into(), Some(strand.to_string())));
        }
        if let Some(status) = candidate.transcript_context_status.as_deref() {
            qualifiers.push(("transcript_context_status".into(), Some(status.to_string())));
        }
        if !candidate.effect_tags.is_empty() {
            qualifiers.push((
                "construct_reasoning_effect_tags".into(),
                Some(candidate.effect_tags.join(",")),
            ));
        }
        if !candidate.supporting_fact_ids.is_empty() {
            qualifiers.push((
                "construct_reasoning_supporting_fact_ids".into(),
                Some(candidate.supporting_fact_ids.join(",")),
            ));
        }
        if !candidate.supporting_decision_ids.is_empty() {
            qualifiers.push((
                "construct_reasoning_supporting_decision_ids".into(),
                Some(candidate.supporting_decision_ids.join(",")),
            ));
        }
        match candidate.role {
            ConstructRole::Promoter => {
                qualifiers.push(("regulatory_class".into(), Some("promoter".to_string())));
            }
            ConstructRole::Enhancer => {
                qualifiers.push(("regulatory_class".into(), Some("enhancer".to_string())));
            }
            ConstructRole::Terminator => {
                qualifiers.push(("regulatory_class".into(), Some("terminator".to_string())));
            }
            ConstructRole::Tfbs => {
                qualifiers.push(("regulatory_class".into(), Some("tfbs".to_string())));
            }
            ConstructRole::SpliceBoundary => {
                qualifiers.push((
                    "regulatory_class".into(),
                    Some("splice_boundary".to_string()),
                ));
            }
            _ => {}
        }
        gb_io::seq::Feature {
            kind: feature_kind.into(),
            location: Self::construct_reasoning_annotation_candidate_feature_location(candidate),
            qualifiers,
        }
    }

    pub fn write_back_construct_reasoning_annotation_candidate(
        &mut self,
        graph_id: &str,
        annotation_id: &str,
    ) -> Result<(ConstructReasoningGraph, AnnotationCandidateWriteback), EngineError> {
        let graph_id = graph_id.trim();
        if graph_id.is_empty() {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: "graph_id cannot be empty".to_string(),
            });
        }
        let annotation_id = annotation_id.trim();
        if annotation_id.is_empty() {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: "annotation_id cannot be empty".to_string(),
            });
        }
        let graph = self.construct_reasoning_graph(graph_id)?;
        let candidate = graph
            .annotation_candidates
            .iter()
            .find(|candidate| candidate.annotation_id == annotation_id)
            .cloned()
            .ok_or_else(|| EngineError {
                code: ErrorCode::NotFound,
                message: format!(
                    "Annotation candidate '{}' not found in construct reasoning graph '{}'",
                    annotation_id, graph_id
                ),
            })?;
        if !matches!(
            candidate.editable_status,
            EditableStatus::Accepted | EditableStatus::Locked
        ) {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: format!(
                    "Annotation candidate '{}' must be accepted before write-back",
                    annotation_id
                ),
            });
        }
        let evidence = graph
            .evidence
            .iter()
            .find(|row| row.evidence_id == candidate.evidence_id)
            .cloned()
            .ok_or_else(|| EngineError {
                code: ErrorCode::NotFound,
                message: format!(
                    "Evidence '{}' referenced by annotation candidate '{}' was not found in graph '{}'",
                    candidate.evidence_id, annotation_id, graph_id
                ),
            })?;
        let mut writeback = AnnotationCandidateWriteback {
            graph_id: graph.graph_id.clone(),
            annotation_id: candidate.annotation_id.clone(),
            evidence_id: candidate.evidence_id.clone(),
            seq_id: graph.seq_id.clone(),
            feature_kind: Self::construct_reasoning_annotation_candidate_feature_kind(
                candidate.role,
            )
            .to_string(),
            feature_label: if candidate.label.trim().is_empty() {
                candidate.role.as_str().to_string()
            } else {
                candidate.label.clone()
            },
            start_0based: candidate.start_0based,
            end_0based_exclusive: candidate.end_0based_exclusive,
            strand: candidate.strand.clone(),
            source_kind: candidate.source_kind.clone(),
            editable_status: candidate.editable_status,
            ..AnnotationCandidateWriteback::default()
        };
        if evidence.provenance_kind == "sequence_feature_annotation" {
            writeback.already_present = true;
            writeback.notes.push(
                "Annotation candidate is already backed by an ordinary sequence feature."
                    .to_string(),
            );
            return Ok((graph, writeback));
        }
        if !Self::construct_reasoning_annotation_candidate_is_writeback_eligible(
            &candidate, &evidence,
        ) {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: format!(
                    "Annotation candidate '{}' is not currently eligible for automatic write-back",
                    annotation_id
                ),
            });
        }
        let seq_id = graph.seq_id.clone();
        let dna = self
            .state
            .sequences
            .get_mut(&seq_id)
            .ok_or_else(|| EngineError {
                code: ErrorCode::NotFound,
                message: format!("Sequence '{}' not found", seq_id),
            })?;
        if dna.features().iter().any(|feature| {
            Self::construct_reasoning_feature_matches_annotation_candidate(
                feature,
                &candidate.annotation_id,
            )
        }) {
            writeback.already_present = true;
            writeback.notes.push(
                "Annotation candidate was already written back as a sequence feature.".to_string(),
            );
            return Ok((graph, writeback));
        }
        dna.features_mut().push(
            Self::build_construct_reasoning_annotation_candidate_feature(
                graph_id, &candidate, &evidence,
            ),
        );
        Self::prepare_sequence(dna);
        writeback.created = true;
        writeback.notes.push(
            "Accepted annotation candidate was written back as an ordinary sequence feature."
                .to_string(),
        );
        let refreshed_graph = self.refresh_construct_reasoning_graph_for_seq_id(&seq_id)?;
        if let Some(updated_candidate) = refreshed_graph.annotation_candidates.iter().find(|row| {
            Self::construct_reasoning_annotation_candidate_span_key(row)
                == Self::construct_reasoning_annotation_candidate_span_key(&candidate)
        }) {
            writeback.editable_status = updated_candidate.editable_status;
        }
        Ok((refreshed_graph, writeback))
    }

    fn construct_reasoning_default_objective(
        seq_id: &str,
        dna: &DNAsequence,
    ) -> ConstructObjective {
        let title = dna
            .name()
            .as_ref()
            .map(|value| value.trim().to_string())
            .filter(|value| !value.is_empty())
            .unwrap_or_else(|| format!("Construct reasoning for {seq_id}"));
        Self::normalize_construct_objective(ConstructObjective {
            objective_id: format!("construct_objective_{}", Self::normalize_id_token(seq_id)),
            title,
            goal: format!(
                "Inspect sequence-linked construct evidence and candidate boundaries on '{seq_id}'"
            ),
            ..ConstructObjective::default()
        })
    }

    fn construct_reasoning_feature_has_cdna_confirmation_hint(
        feature: &gb_io::seq::Feature,
    ) -> bool {
        for key in [
            "experiment",
            "evidence",
            "inference",
            "note",
            "support",
            "transcript_support_level",
        ] {
            for value in feature.qualifier_values(key) {
                let normalized = value.trim().to_ascii_lowercase();
                if normalized.contains("cdna")
                    || normalized.contains("c dna")
                    || normalized.contains("mrna-supported")
                    || normalized.contains("mrna supported")
                    || normalized.contains("validated by mrna")
                {
                    return true;
                }
            }
        }
        false
    }

    fn construct_reasoning_role_from_feature(
        feature: &gb_io::seq::Feature,
    ) -> Option<ConstructRole> {
        let kind = feature.kind.to_string().to_ascii_uppercase();
        match kind.as_str() {
            "GENE" => Some(ConstructRole::Gene),
            "MRNA" | "TRANSCRIPT" | "NCRNA" | "RRNA" | "TRNA" => Some(ConstructRole::Transcript),
            "EXON" => Some(ConstructRole::Exon),
            "CDS" => Some(ConstructRole::Cds),
            "5'UTR" | "5UTR" | "FIVE_PRIME_UTR" => Some(ConstructRole::Utr5Prime),
            "3'UTR" | "3UTR" | "THREE_PRIME_UTR" => Some(ConstructRole::Utr3Prime),
            "PROMOTER" => Some(ConstructRole::Promoter),
            "ENHANCER" => Some(ConstructRole::Enhancer),
            "TERMINATOR" => Some(ConstructRole::Terminator),
            "REPEAT_REGION" | "REPEAT" => Some(ConstructRole::RepeatRegion),
            "MOBILE_ELEMENT" | "SINE" | "LINE" | "LTR" | "RETROTRANSPOSON" => {
                Some(ConstructRole::MobileElement)
            }
            "VARIATION" | "VARIANT" | "SNP" | "SNV" | "MUTATION" => Some(ConstructRole::Variant),
            "SIG_PEPTIDE" | "SIGNAL_PEPTIDE" => Some(ConstructRole::SignalPeptide),
            "TFBS" | "TF_BINDING_SITE" | "PROTEIN_BIND" => Some(ConstructRole::Tfbs),
            "REGULATORY" | "MISC_FEATURE" => {
                let hint = Self::first_nonempty_feature_qualifier(
                    feature,
                    &["regulatory_class", "function", "note", "label", "name"],
                )
                .unwrap_or_default()
                .to_ascii_lowercase();
                if hint.contains("promoter") {
                    Some(ConstructRole::Promoter)
                } else if hint.contains("enhancer") {
                    Some(ConstructRole::Enhancer)
                } else if hint.contains("terminator") {
                    Some(ConstructRole::Terminator)
                } else if hint.contains("alu")
                    || hint.contains("sine")
                    || hint.contains("retrotransposon")
                {
                    Some(ConstructRole::MobileElement)
                } else if hint.contains("repeat") || hint.contains("low complexity") {
                    Some(ConstructRole::RepeatRegion)
                } else if hint.contains("splice") {
                    Some(ConstructRole::SpliceBoundary)
                } else if hint.contains("tfbs") || hint.contains("binding") {
                    Some(ConstructRole::Tfbs)
                } else {
                    None
                }
            }
            _ => None,
        }
    }

    fn construct_reasoning_class_for_feature(
        feature: &gb_io::seq::Feature,
        role: ConstructRole,
    ) -> EvidenceClass {
        if role == ConstructRole::Variant {
            let generated = Self::feature_qualifier_text(feature, "gentle_generated")
                .unwrap_or_default()
                .trim()
                .to_ascii_lowercase();
            if generated == GENOME_VCF_TRACK_GENERATED_TAG
                || generated == DBSNP_VARIANT_MARKER_GENERATED_TAG
            {
                return EvidenceClass::HardFact;
            }
            return EvidenceClass::ReliableAnnotation;
        }
        if role == ConstructRole::Tfbs {
            return EvidenceClass::SoftHypothesis;
        }
        if role == ConstructRole::MobileElement {
            return EvidenceClass::SoftHypothesis;
        }
        if matches!(role, ConstructRole::Exon | ConstructRole::SpliceBoundary)
            && Self::construct_reasoning_feature_has_cdna_confirmation_hint(feature)
        {
            return EvidenceClass::HardFact;
        }
        EvidenceClass::ReliableAnnotation
    }

    fn construct_reasoning_feature_context_tags(feature: &gb_io::seq::Feature) -> Vec<String> {
        let mut tags = vec![feature.kind.to_string().to_ascii_lowercase()];
        if Self::construct_reasoning_feature_has_cdna_confirmation_hint(feature) {
            tags.push("cdna_confirmed".to_string());
        }
        if let Some(generated) = Self::feature_qualifier_text(feature, "gentle_generated")
            .map(|value| value.trim().to_ascii_lowercase())
            .filter(|value| !value.is_empty())
        {
            tags.push(generated);
        }
        if let Some(regulatory_class) = Self::feature_qualifier_text(feature, "regulatory_class") {
            tags.push(regulatory_class.trim().to_ascii_lowercase());
        }
        if let Some(variant_class) = Self::feature_qualifier_text(feature, "vcf_variant_class")
            .map(|value| value.trim().to_ascii_lowercase())
            .filter(|value| !value.is_empty())
        {
            tags.push(variant_class);
        }
        tags.sort();
        tags.dedup();
        tags
    }

    fn construct_reasoning_feature_rationale(
        feature: &gb_io::seq::Feature,
        role: ConstructRole,
        evidence_class: EvidenceClass,
    ) -> String {
        if role == ConstructRole::Variant {
            let generated = Self::feature_qualifier_text(feature, "gentle_generated")
                .unwrap_or_default()
                .trim()
                .to_ascii_lowercase();
            if generated == GENOME_VCF_TRACK_GENERATED_TAG {
                return "Variant marker imported from a VCF track; locus and observed allele annotation are treated as an inspectable hard fact.".to_string();
            }
            if generated == DBSNP_VARIANT_MARKER_GENERATED_TAG {
                return "Variant marker imported from dbSNP / NCBI Variation; locus identity is treated as an inspectable hard fact while downstream functional impact remains to be reasoned from context.".to_string();
            }
            return "Variant annotation imported from the sequence record as inspectable construct-planning context.".to_string();
        }
        if role == ConstructRole::Tfbs {
            return "TFBS-style annotation is kept as a soft hypothesis because binding evidence is context-sensitive.".to_string();
        }
        if role == ConstructRole::MobileElement {
            return "Mobile-element style annotation is kept as a soft hypothesis until a curated repeat-family catalog confirms the family assignment.".to_string();
        }
        if role == ConstructRole::RepeatRegion {
            return "Repeat/similarity annotation imported from the sequence record as inspectable structural context.".to_string();
        }
        if evidence_class == EvidenceClass::HardFact {
            return format!(
                "{} annotation is treated as a hard fact because its qualifiers explicitly hint at cDNA/mRNA confirmation.",
                role.as_str()
            );
        }
        format!(
            "{} annotation imported from the sequence record as reliable construct-planning context.",
            role.as_str()
        )
    }

    fn construct_reasoning_nucleotide_index(base: u8) -> Option<usize> {
        match base.to_ascii_uppercase() {
            b'A' => Some(0),
            b'C' => Some(1),
            b'G' => Some(2),
            b'T' => Some(3),
            _ => None,
        }
    }

    fn construct_reasoning_longest_homopolymer_run(bytes: &[u8]) -> usize {
        let mut best = 0usize;
        let mut current = 0usize;
        let mut current_base = 0u8;
        for base in bytes.iter().map(|value| value.to_ascii_uppercase()) {
            if Self::construct_reasoning_nucleotide_index(base).is_none() {
                current = 0;
                current_base = 0;
                continue;
            }
            if current > 0 && base == current_base {
                current += 1;
            } else {
                current = 1;
                current_base = base;
            }
            best = best.max(current);
        }
        best
    }

    fn construct_reasoning_normalized_entropy(bytes: &[u8]) -> f64 {
        let mut counts = [0usize; 4];
        let mut total = 0usize;
        for base in bytes {
            if let Some(idx) = Self::construct_reasoning_nucleotide_index(*base) {
                counts[idx] += 1;
                total += 1;
            }
        }
        if total == 0 {
            return 0.0;
        }
        let entropy = counts
            .into_iter()
            .filter(|count| *count > 0)
            .map(|count| {
                let p = count as f64 / total as f64;
                -(p * p.log2())
            })
            .sum::<f64>();
        (entropy / 2.0).clamp(0.0, 1.0)
    }

    fn construct_reasoning_unique_kmer_ratio(bytes: &[u8], k: usize) -> f64 {
        if k == 0 {
            return 1.0;
        }
        let Some(total) = bytes.len().checked_sub(k).map(|value| value + 1) else {
            return 1.0;
        };
        let mut unique = HashSet::new();
        let mut valid = 0usize;
        for window in bytes.windows(k) {
            if window
                .iter()
                .all(|base| Self::construct_reasoning_nucleotide_index(*base).is_some())
            {
                unique.insert(
                    window
                        .iter()
                        .map(|value| value.to_ascii_uppercase())
                        .collect::<Vec<_>>(),
                );
                valid += 1;
            }
        }
        if valid == 0 {
            return 0.0;
        }
        (unique.len() as f64 / total as f64).clamp(0.0, 1.0)
    }

    fn construct_reasoning_most_common_base_run(bytes: &[u8]) -> Option<(char, usize)> {
        let mut best = 0usize;
        let mut best_base = None;
        let mut current = 0usize;
        let mut current_base = 0u8;
        for base in bytes.iter().map(|value| value.to_ascii_uppercase()) {
            if Self::construct_reasoning_nucleotide_index(base).is_none() {
                current = 0;
                current_base = 0;
                continue;
            }
            if current > 0 && base == current_base {
                current += 1;
            } else {
                current = 1;
                current_base = base;
            }
            if current > best {
                best = current;
                best_base = Some(base as char);
            }
        }
        best_base.map(|base| (base, best))
    }

    fn construct_reasoning_kmer_is_low_value(kmer: &[u8]) -> bool {
        if kmer
            .iter()
            .any(|base| Self::construct_reasoning_nucleotide_index(*base).is_none())
        {
            return true;
        }
        let distinct = kmer
            .iter()
            .map(|base| base.to_ascii_uppercase())
            .collect::<BTreeSet<_>>()
            .len();
        distinct < 3
            || Self::construct_reasoning_longest_homopolymer_run(kmer)
                >= kmer.len().saturating_sub(2)
    }

    fn construct_reasoning_shared_kmer_fraction(left: &[u8], right: &[u8], k: usize) -> f64 {
        if k == 0 || left.len() < k || right.len() < k {
            return 0.0;
        }
        let left_set = left
            .windows(k)
            .filter(|window| {
                window
                    .iter()
                    .all(|base| Self::construct_reasoning_nucleotide_index(*base).is_some())
            })
            .map(|window| {
                window
                    .iter()
                    .map(|value| value.to_ascii_uppercase())
                    .collect::<Vec<_>>()
            })
            .collect::<BTreeSet<_>>();
        let right_set = right
            .windows(k)
            .filter(|window| {
                window
                    .iter()
                    .all(|base| Self::construct_reasoning_nucleotide_index(*base).is_some())
            })
            .map(|window| {
                window
                    .iter()
                    .map(|value| value.to_ascii_uppercase())
                    .collect::<Vec<_>>()
            })
            .collect::<BTreeSet<_>>();
        let denom = left_set.len().min(right_set.len());
        if denom == 0 {
            return 0.0;
        }
        let shared = left_set.intersection(&right_set).count();
        (shared as f64 / denom as f64).clamp(0.0, 1.0)
    }

    fn construct_reasoning_push_generated_sequence_evidence(
        evidence: &mut Vec<DesignEvidence>,
        seq_id: &str,
        evidence_id: String,
        start_0based: usize,
        end_0based_exclusive: usize,
        strand: Option<String>,
        role: ConstructRole,
        evidence_class: EvidenceClass,
        label: String,
        rationale: String,
        score: Option<f64>,
        confidence: Option<f64>,
        specificity_bias: Option<f64>,
        sensitivity_bias: Option<f64>,
        mut context_tags: Vec<String>,
        mut provenance_refs: Vec<String>,
        notes: Vec<String>,
    ) {
        context_tags.push("generated".to_string());
        context_tags.sort();
        context_tags.dedup();
        provenance_refs.sort();
        provenance_refs.dedup();
        evidence.push(DesignEvidence {
            evidence_id,
            seq_id: seq_id.to_string(),
            start_0based,
            end_0based_exclusive,
            strand,
            role,
            evidence_class,
            label,
            rationale,
            score,
            confidence,
            specificity_bias,
            sensitivity_bias,
            context_tags,
            provenance_kind: "derived_sequence_similarity".to_string(),
            provenance_refs,
            notes,
            ..DesignEvidence::default()
        });
    }

    fn build_construct_reasoning_sequence_similarity_evidence(
        seq_id: &str,
        dna: &DNAsequence,
    ) -> Vec<DesignEvidence> {
        #[derive(Clone, Debug)]
        struct LowComplexitySpan {
            start_0based: usize,
            end_0based_exclusive: usize,
            complexity_score: f64,
            normalized_entropy: f64,
            unique_kmer_ratio: f64,
            longest_homopolymer: usize,
        }

        #[derive(Clone, Debug)]
        struct MotifRepeatSpan {
            start_0based: usize,
            end_0based_exclusive: usize,
            label: String,
            risk_score: f64,
            context_tags: Vec<String>,
            notes: Vec<String>,
        }

        #[derive(Clone, Debug)]
        struct RepeatClusterSpan {
            start_0based: usize,
            end_0based_exclusive: usize,
            label: String,
            seed_examples: Vec<String>,
            pair_count: usize,
            risk_score: f64,
            context_tags: Vec<String>,
            notes: Vec<String>,
        }

        #[derive(Clone, Debug)]
        struct AluLikeSpan {
            start_0based: usize,
            end_0based_exclusive: usize,
            strand: String,
            tail_length_bp: usize,
            arm_similarity: f64,
            candidate_length_bp: usize,
            score: f64,
        }

        fn merge_low_complexity_spans(spans: Vec<LowComplexitySpan>) -> Vec<LowComplexitySpan> {
            let mut spans = spans;
            spans.sort_by_key(|row| (row.start_0based, row.end_0based_exclusive));
            let mut merged: Vec<LowComplexitySpan> = vec![];
            for span in spans {
                if let Some(last) = merged.last_mut() {
                    if span.start_0based
                        <= last
                            .end_0based_exclusive
                            .saturating_add(LOW_COMPLEXITY_STEP_BP)
                    {
                        last.end_0based_exclusive =
                            last.end_0based_exclusive.max(span.end_0based_exclusive);
                        last.complexity_score = last.complexity_score.min(span.complexity_score);
                        last.normalized_entropy =
                            last.normalized_entropy.min(span.normalized_entropy);
                        last.unique_kmer_ratio = last.unique_kmer_ratio.min(span.unique_kmer_ratio);
                        last.longest_homopolymer =
                            last.longest_homopolymer.max(span.longest_homopolymer);
                        continue;
                    }
                }
                merged.push(span);
            }
            merged
        }

        fn merge_repeat_cluster_spans(
            mut spans: Vec<RepeatClusterSpan>,
            max_gap_bp: usize,
        ) -> Vec<RepeatClusterSpan> {
            spans.sort_by_key(|row| (row.start_0based, row.end_0based_exclusive));
            let mut merged: Vec<RepeatClusterSpan> = vec![];
            for span in spans {
                if let Some(last) = merged.last_mut() {
                    if last.label == span.label
                        && span.start_0based <= last.end_0based_exclusive.saturating_add(max_gap_bp)
                    {
                        last.end_0based_exclusive =
                            last.end_0based_exclusive.max(span.end_0based_exclusive);
                        last.pair_count = last.pair_count.max(span.pair_count);
                        last.risk_score = last.risk_score.max(span.risk_score);
                        last.seed_examples.extend(span.seed_examples);
                        last.seed_examples.sort();
                        last.seed_examples.dedup();
                        last.context_tags.extend(span.context_tags);
                        last.context_tags.sort();
                        last.context_tags.dedup();
                        last.notes.extend(span.notes);
                        last.notes.sort();
                        last.notes.dedup();
                        continue;
                    }
                }
                merged.push(span);
            }
            merged
        }

        fn merge_alu_like_spans(mut spans: Vec<AluLikeSpan>) -> Vec<AluLikeSpan> {
            spans.sort_by_key(|row| (row.start_0based, row.end_0based_exclusive));
            let mut merged: Vec<AluLikeSpan> = vec![];
            for span in spans {
                if let Some(last) = merged.last_mut() {
                    if last.strand == span.strand
                        && span.start_0based <= last.end_0based_exclusive.saturating_add(24)
                    {
                        if span.score > last.score {
                            last.start_0based = last.start_0based.min(span.start_0based);
                            last.end_0based_exclusive =
                                last.end_0based_exclusive.max(span.end_0based_exclusive);
                            last.tail_length_bp = last.tail_length_bp.max(span.tail_length_bp);
                            last.arm_similarity = last.arm_similarity.max(span.arm_similarity);
                            last.candidate_length_bp =
                                last.end_0based_exclusive.saturating_sub(last.start_0based);
                            last.score = span.score.max(last.score);
                        }
                        continue;
                    }
                }
                merged.push(span);
            }
            merged
        }

        let bytes = dna
            .forward_bytes()
            .iter()
            .map(|base| base.to_ascii_uppercase())
            .collect::<Vec<_>>();
        if bytes.is_empty() {
            return vec![];
        }

        let mut evidence = vec![];

        let low_complexity_window = LOW_COMPLEXITY_WINDOW_BP.min(bytes.len()).max(1);
        let low_complexity_step = LOW_COMPLEXITY_STEP_BP.min(low_complexity_window).max(1);
        let mut low_complexity_raw = vec![];
        let mut scan_starts = (0..=bytes.len().saturating_sub(low_complexity_window))
            .step_by(low_complexity_step)
            .collect::<Vec<_>>();
        if scan_starts.is_empty() {
            scan_starts.push(0);
        } else {
            let last_start = bytes.len().saturating_sub(low_complexity_window);
            if scan_starts.last().copied() != Some(last_start) {
                scan_starts.push(last_start);
            }
        }
        for start_0based in scan_starts {
            let end_0based_exclusive = start_0based + low_complexity_window;
            let window = &bytes[start_0based..end_0based_exclusive];
            let normalized_entropy = Self::construct_reasoning_normalized_entropy(window);
            let unique_kmer_ratio = Self::construct_reasoning_unique_kmer_ratio(window, 3);
            let longest_homopolymer = Self::construct_reasoning_longest_homopolymer_run(window);
            let complexity_score =
                (normalized_entropy * 0.55 + unique_kmer_ratio * 0.45).clamp(0.0, 1.0);
            if complexity_score <= LOW_COMPLEXITY_SCORE_THRESHOLD
                || longest_homopolymer >= LOW_COMPLEXITY_HOMOPOLYMER_THRESHOLD_BP
            {
                low_complexity_raw.push(LowComplexitySpan {
                    start_0based,
                    end_0based_exclusive,
                    complexity_score,
                    normalized_entropy,
                    unique_kmer_ratio,
                    longest_homopolymer,
                });
            }
        }
        for (idx, span) in merge_low_complexity_spans(low_complexity_raw)
            .into_iter()
            .enumerate()
        {
            let slice = &bytes[span.start_0based..span.end_0based_exclusive];
            let normalized_entropy = Self::construct_reasoning_normalized_entropy(slice);
            let unique_kmer_ratio = Self::construct_reasoning_unique_kmer_ratio(slice, 3);
            let longest_homopolymer = Self::construct_reasoning_longest_homopolymer_run(slice);
            let complexity_score =
                (normalized_entropy * 0.55 + unique_kmer_ratio * 0.45).clamp(0.0, 1.0);
            let risk_score = (1.0 - complexity_score).max(
                ((longest_homopolymer.saturating_sub(LOW_COMPLEXITY_HOMOPOLYMER_THRESHOLD_BP))
                    as f64
                    / 12.0)
                    .clamp(0.0, 1.0),
            );
            let mut tags = vec![
                "low_complexity".to_string(),
                "repeat_region".to_string(),
                "mapping_ambiguity_risk".to_string(),
                "nanopore_signal_risk".to_string(),
            ];
            if longest_homopolymer >= HOMOPOLYMER_ANNOTATION_THRESHOLD_BP {
                tags.push("polymerase_slippage_risk".to_string());
                tags.push("nanopore_homopolymer_risk".to_string());
            }
            let dominant_run = Self::construct_reasoning_most_common_base_run(slice)
                .map(|(base, run)| format!("dominant_run={base}x{run}"))
                .unwrap_or_else(|| "dominant_run=unknown".to_string());
            Self::construct_reasoning_push_generated_sequence_evidence(
                &mut evidence,
                seq_id,
                format!(
                    "similarity_low_complexity_{}_{}_{}",
                    span.start_0based, span.end_0based_exclusive, idx
                ),
                span.start_0based,
                span.end_0based_exclusive,
                Some("+".to_string()),
                ConstructRole::RepeatRegion,
                EvidenceClass::ContextEvidence,
                "Low-complexity region".to_string(),
                "Merged sequence window shows reduced composition entropy and/or k-mer diversity, so it is kept as generated repeat/similarity context rather than a hard biological claim.".to_string(),
                Some(risk_score),
                Some(0.72),
                Some(0.52),
                Some(0.83),
                tags,
                vec!["low_complexity_scan".to_string()],
                vec![
                    format!("complexity_score={complexity_score:.3}"),
                    format!("normalized_entropy={normalized_entropy:.3}"),
                    format!("unique_trimer_ratio={unique_kmer_ratio:.3}"),
                    format!("longest_homopolymer={longest_homopolymer}"),
                    dominant_run,
                ],
            );
        }

        let mut tandem_repeat_spans = vec![];
        let mut pos = 0usize;
        while pos < bytes.len() {
            let mut best: Option<MotifRepeatSpan> = None;
            for motif_len in 1..=6usize {
                if pos + motif_len > bytes.len() {
                    break;
                }
                let motif = &bytes[pos..pos + motif_len];
                if motif
                    .iter()
                    .any(|base| Self::construct_reasoning_nucleotide_index(*base).is_none())
                {
                    continue;
                }
                let mut repeat_count = 1usize;
                while pos + (repeat_count + 1) * motif_len <= bytes.len()
                    && bytes[pos + repeat_count * motif_len..pos + (repeat_count + 1) * motif_len]
                        == *motif
                {
                    repeat_count += 1;
                }
                let total_len = repeat_count * motif_len;
                let meets_threshold = if motif_len == 1 {
                    repeat_count >= HOMOPOLYMER_ANNOTATION_THRESHOLD_BP
                } else {
                    repeat_count >= 3 && total_len >= TANDEM_REPEAT_MIN_TOTAL_BP
                };
                if !meets_threshold {
                    continue;
                }
                let label = if motif_len == 1 {
                    format!("Homopolymer run ({} x{repeat_count})", motif[0] as char)
                } else {
                    format!(
                        "Tandem repeat (({}) x{repeat_count})",
                        String::from_utf8_lossy(motif)
                    )
                };
                let mut tags = vec![
                    "repeat_region".to_string(),
                    if motif_len == 1 {
                        "homopolymer_run".to_string()
                    } else {
                        "tandem_repeat".to_string()
                    },
                    "polymerase_slippage_risk".to_string(),
                    "nanopore_signal_risk".to_string(),
                ];
                if motif_len == 1 {
                    tags.push("nanopore_homopolymer_risk".to_string());
                }
                if motif_len == 1 || motif_len <= 3 {
                    tags.push("low_complexity".to_string());
                }
                let risk_score = ((total_len as f64 / 32.0).min(1.0)
                    + (repeat_count as f64 / 8.0).min(1.0))
                    / 2.0;
                let span = MotifRepeatSpan {
                    start_0based: pos,
                    end_0based_exclusive: pos + total_len,
                    label,
                    risk_score: risk_score.clamp(0.0, 1.0),
                    context_tags: tags,
                    notes: vec![
                        format!("motif_len={motif_len}"),
                        format!("repeat_count={repeat_count}"),
                        format!("total_bp={total_len}"),
                    ],
                };
                if best
                    .as_ref()
                    .map(|candidate| {
                        candidate
                            .end_0based_exclusive
                            .saturating_sub(candidate.start_0based)
                    })
                    .unwrap_or_default()
                    < total_len
                {
                    best = Some(span);
                }
            }
            if let Some(span) = best {
                pos = span.end_0based_exclusive;
                tandem_repeat_spans.push(span);
            } else {
                pos += 1;
            }
        }
        for (idx, span) in tandem_repeat_spans.into_iter().enumerate() {
            Self::construct_reasoning_push_generated_sequence_evidence(
                &mut evidence,
                seq_id,
                format!(
                    "similarity_repeat_{}_{}_{}",
                    span.start_0based, span.end_0based_exclusive, idx
                ),
                span.start_0based,
                span.end_0based_exclusive,
                Some("+".to_string()),
                ConstructRole::RepeatRegion,
                EvidenceClass::ContextEvidence,
                span.label,
                "Exact tandem motif reuse is a sequence-determined repeat pattern. GENtle keeps it as generated repeat context because downstream operational risk depends on the experiment.".to_string(),
                Some(span.risk_score),
                Some(0.88),
                Some(0.81),
                Some(0.63),
                span.context_tags,
                vec!["exact_tandem_repeat_scan".to_string()],
                span.notes,
            );
        }

        let mut direct_repeat_raw = vec![];
        if bytes.len() >= DIRECT_REPEAT_SEED_BP {
            let mut seed_positions: HashMap<Vec<u8>, Vec<usize>> = HashMap::new();
            for (start_0based, window) in bytes.windows(DIRECT_REPEAT_SEED_BP).enumerate() {
                if Self::construct_reasoning_kmer_is_low_value(window) {
                    continue;
                }
                seed_positions
                    .entry(window.to_vec())
                    .or_default()
                    .push(start_0based);
            }
            for (seed, positions) in seed_positions {
                if positions.len() < 2 {
                    continue;
                }
                let mut cluster_start = positions[0];
                let mut cluster_end = positions[0] + DIRECT_REPEAT_SEED_BP;
                let mut cluster_count = 1usize;
                for position in positions.iter().copied().skip(1) {
                    if position <= cluster_end.saturating_add(DIRECT_REPEAT_MAX_CLUSTER_GAP_BP) {
                        cluster_end = position + DIRECT_REPEAT_SEED_BP;
                        cluster_count += 1;
                    } else {
                        if cluster_count >= 2 {
                            let risk_score = ((cluster_count as f64 / 4.0).min(1.0)
                                + ((cluster_end.saturating_sub(cluster_start)) as f64 / 600.0)
                                    .min(1.0))
                                / 2.0;
                            let mut tags = vec![
                                "repeat_region".to_string(),
                                "direct_repeat".to_string(),
                                "repeat_dense".to_string(),
                                "mapping_ambiguity_risk".to_string(),
                            ];
                            if cluster_count >= 3 {
                                tags.push("cloning_stability_risk".to_string());
                            }
                            direct_repeat_raw.push(RepeatClusterSpan {
                                start_0based: cluster_start,
                                end_0based_exclusive: cluster_end,
                                label: "Direct repeat cluster".to_string(),
                                seed_examples: vec![String::from_utf8_lossy(&seed).to_string()],
                                pair_count: cluster_count,
                                risk_score: risk_score.clamp(0.0, 1.0),
                                context_tags: tags,
                                notes: vec![format!("repeat_seed_count={cluster_count}")],
                            });
                        }
                        cluster_start = position;
                        cluster_end = position + DIRECT_REPEAT_SEED_BP;
                        cluster_count = 1;
                    }
                }
                if cluster_count >= 2 {
                    let risk_score = ((cluster_count as f64 / 4.0).min(1.0)
                        + ((cluster_end.saturating_sub(cluster_start)) as f64 / 600.0).min(1.0))
                        / 2.0;
                    let mut tags = vec![
                        "repeat_region".to_string(),
                        "direct_repeat".to_string(),
                        "repeat_dense".to_string(),
                        "mapping_ambiguity_risk".to_string(),
                    ];
                    if cluster_count >= 3 {
                        tags.push("cloning_stability_risk".to_string());
                    }
                    direct_repeat_raw.push(RepeatClusterSpan {
                        start_0based: cluster_start,
                        end_0based_exclusive: cluster_end,
                        label: "Direct repeat cluster".to_string(),
                        seed_examples: vec![String::from_utf8_lossy(&seed).to_string()],
                        pair_count: cluster_count,
                        risk_score: risk_score.clamp(0.0, 1.0),
                        context_tags: tags,
                        notes: vec![format!("repeat_seed_count={cluster_count}")],
                    });
                }
            }
        }
        for (idx, span) in merge_repeat_cluster_spans(direct_repeat_raw, DIRECT_REPEAT_SEED_BP)
            .into_iter()
            .enumerate()
        {
            let mut notes = span.notes;
            notes.push(format!("seed_examples={}", span.seed_examples.join(",")));
            Self::construct_reasoning_push_generated_sequence_evidence(
                &mut evidence,
                seq_id,
                format!(
                    "similarity_direct_repeat_{}_{}_{}",
                    span.start_0based, span.end_0based_exclusive, idx
                ),
                span.start_0based,
                span.end_0based_exclusive,
                Some("+".to_string()),
                ConstructRole::RepeatRegion,
                EvidenceClass::ContextEvidence,
                span.label,
                "Exact local k-mer reuse suggests a direct-repeat cluster. GENtle keeps this as generated similarity context because downstream risk depends on whether PCR, mapping, or cloning stability matters for the current task.".to_string(),
                Some(span.risk_score),
                Some(0.84),
                Some(0.79),
                Some(0.67),
                span.context_tags,
                vec!["direct_repeat_cluster_scan".to_string()],
                notes,
            );
        }

        let mut inverted_repeat_raw = vec![];
        if bytes.len() >= INVERTED_REPEAT_SEED_BP {
            let mut seed_positions: HashMap<Vec<u8>, Vec<usize>> = HashMap::new();
            for (start_0based, window) in bytes.windows(INVERTED_REPEAT_SEED_BP).enumerate() {
                if Self::construct_reasoning_kmer_is_low_value(window) {
                    continue;
                }
                seed_positions
                    .entry(window.to_vec())
                    .or_default()
                    .push(start_0based);
            }
            let mut seen_pairs = HashSet::new();
            for (seed, positions) in &seed_positions {
                let rc_seed = Self::reverse_complement_bytes(seed);
                let Some(partners) = seed_positions.get(&rc_seed) else {
                    continue;
                };
                for left in positions.iter().copied() {
                    for right in partners.iter().copied() {
                        if right <= left {
                            continue;
                        }
                        if right.saturating_sub(left) > INVERTED_REPEAT_MAX_SPAN_BP {
                            break;
                        }
                        if !seen_pairs.insert((left, right)) {
                            continue;
                        }
                        let span_len = right + INVERTED_REPEAT_SEED_BP - left;
                        let risk_score = ((span_len as f64 / 220.0).min(1.0) + 0.6) / 1.6;
                        inverted_repeat_raw.push(RepeatClusterSpan {
                            start_0based: left,
                            end_0based_exclusive: right + INVERTED_REPEAT_SEED_BP,
                            label: "Inverted repeat cluster".to_string(),
                            seed_examples: vec![String::from_utf8_lossy(seed).to_string()],
                            pair_count: 1,
                            risk_score: risk_score.clamp(0.0, 1.0),
                            context_tags: vec![
                                "repeat_region".to_string(),
                                "inverted_repeat".to_string(),
                                "palindromic_repeat_risk".to_string(),
                                "inversion_risk".to_string(),
                                "cloning_stability_risk".to_string(),
                            ],
                            notes: vec![format!("span_bp={span_len}")],
                        });
                    }
                }
            }
        }
        for (idx, span) in merge_repeat_cluster_spans(inverted_repeat_raw, INVERTED_REPEAT_SEED_BP)
            .into_iter()
            .enumerate()
        {
            let mut notes = span.notes;
            notes.push(format!("seed_examples={}", span.seed_examples.join(",")));
            notes.push("interpretation=palindromic_or_inverted_seed_reuse".to_string());
            Self::construct_reasoning_push_generated_sequence_evidence(
                &mut evidence,
                seq_id,
                format!(
                    "similarity_inverted_repeat_{}_{}_{}",
                    span.start_0based, span.end_0based_exclusive, idx
                ),
                span.start_0based,
                span.end_0based_exclusive,
                Some("+".to_string()),
                ConstructRole::RepeatRegion,
                EvidenceClass::ContextEvidence,
                span.label,
                "Local reverse-complement seed reuse suggests an inverted-repeat/palindromic cluster, which can matter for inversion or cloning-instability review without being a standalone biological fact.".to_string(),
                Some(span.risk_score),
                Some(0.82),
                Some(0.76),
                Some(0.61),
                span.context_tags,
                vec!["inverted_repeat_cluster_scan".to_string()],
                notes,
            );
        }

        let mut alu_like_raw = vec![];
        let mut poly_a_runs = vec![];
        let mut run_start = None;
        for (idx, base) in bytes.iter().copied().enumerate() {
            if base == b'A' {
                run_start.get_or_insert(idx);
            } else if let Some(start) = run_start.take() {
                if idx.saturating_sub(start) >= ALU_LIKE_POLYA_MIN_BP {
                    poly_a_runs.push((start, idx));
                }
            }
        }
        if let Some(start) = run_start.take() {
            if bytes.len().saturating_sub(start) >= ALU_LIKE_POLYA_MIN_BP {
                poly_a_runs.push((start, bytes.len()));
            }
        }
        for (tail_start, tail_end) in poly_a_runs {
            for candidate_length_bp in (ALU_LIKE_MIN_BP..=ALU_LIKE_MAX_BP).step_by(20) {
                if tail_start < candidate_length_bp {
                    continue;
                }
                let start_0based = tail_start - candidate_length_bp;
                let end_0based_exclusive = tail_end;
                let body = &bytes[start_0based..tail_start];
                let body_mid = body.len() / 2;
                if body_mid < 60 {
                    continue;
                }
                let left = &body[..body_mid];
                let right = &body[body_mid..];
                let arm_similarity = Self::construct_reasoning_shared_kmer_fraction(left, right, 8);
                let linker_start = body_mid.saturating_sub(12);
                let linker_end = (body_mid + 12).min(body.len());
                let linker = &body[linker_start..linker_end];
                let at_fraction = if linker.is_empty() {
                    0.0
                } else {
                    linker
                        .iter()
                        .filter(|base| matches!(base.to_ascii_uppercase(), b'A' | b'T'))
                        .count() as f64
                        / linker.len() as f64
                };
                let score = (arm_similarity * 0.55
                    + at_fraction * 0.15
                    + ((tail_end - tail_start) as f64 / 24.0).min(1.0) * 0.30)
                    .clamp(0.0, 1.0);
                if arm_similarity >= 0.18 && score >= 0.55 {
                    alu_like_raw.push(AluLikeSpan {
                        start_0based,
                        end_0based_exclusive,
                        strand: "+".to_string(),
                        tail_length_bp: tail_end - tail_start,
                        arm_similarity,
                        candidate_length_bp: end_0based_exclusive - start_0based,
                        score,
                    });
                }
            }
        }

        let mut poly_t_runs = vec![];
        let mut run_start = None;
        for (idx, base) in bytes.iter().copied().enumerate() {
            if base == b'T' {
                run_start.get_or_insert(idx);
            } else if let Some(start) = run_start.take() {
                if idx.saturating_sub(start) >= ALU_LIKE_POLYA_MIN_BP {
                    poly_t_runs.push((start, idx));
                }
            }
        }
        if let Some(start) = run_start.take() {
            if bytes.len().saturating_sub(start) >= ALU_LIKE_POLYA_MIN_BP {
                poly_t_runs.push((start, bytes.len()));
            }
        }
        for (tail_start, tail_end) in poly_t_runs {
            for candidate_length_bp in (ALU_LIKE_MIN_BP..=ALU_LIKE_MAX_BP).step_by(20) {
                if tail_end + candidate_length_bp > bytes.len() {
                    continue;
                }
                let start_0based = tail_start;
                let end_0based_exclusive = tail_end + candidate_length_bp;
                let body = Self::reverse_complement_bytes(&bytes[tail_end..end_0based_exclusive]);
                let body_mid = body.len() / 2;
                if body_mid < 60 {
                    continue;
                }
                let left = &body[..body_mid];
                let right = &body[body_mid..];
                let arm_similarity = Self::construct_reasoning_shared_kmer_fraction(left, right, 8);
                let linker_start = body_mid.saturating_sub(12);
                let linker_end = (body_mid + 12).min(body.len());
                let linker = &body[linker_start..linker_end];
                let at_fraction = if linker.is_empty() {
                    0.0
                } else {
                    linker
                        .iter()
                        .filter(|base| matches!(base.to_ascii_uppercase(), b'A' | b'T'))
                        .count() as f64
                        / linker.len() as f64
                };
                let score = (arm_similarity * 0.55
                    + at_fraction * 0.15
                    + ((tail_end - tail_start) as f64 / 24.0).min(1.0) * 0.30)
                    .clamp(0.0, 1.0);
                if arm_similarity >= 0.18 && score >= 0.55 {
                    alu_like_raw.push(AluLikeSpan {
                        start_0based,
                        end_0based_exclusive,
                        strand: "-".to_string(),
                        tail_length_bp: tail_end - tail_start,
                        arm_similarity,
                        candidate_length_bp: end_0based_exclusive - start_0based,
                        score,
                    });
                }
            }
        }
        for (idx, span) in merge_alu_like_spans(alu_like_raw).into_iter().enumerate() {
            Self::construct_reasoning_push_generated_sequence_evidence(
                &mut evidence,
                seq_id,
                format!(
                    "similarity_alu_like_{}_{}_{}",
                    span.start_0based, span.end_0based_exclusive, idx
                ),
                span.start_0based,
                span.end_0based_exclusive,
                Some(span.strand.clone()),
                ConstructRole::MobileElement,
                EvidenceClass::SoftHypothesis,
                "Alu-like SINE candidate".to_string(),
                "Region shows a two-arm repeat structure with an A-rich tail consistent with a cautious Alu-like SINE heuristic. This remains a soft hypothesis until a curated repeat-family catalog is integrated.".to_string(),
                Some(span.score),
                Some(0.46),
                Some(0.71),
                Some(0.37),
                vec![
                    "alu_like".to_string(),
                    "mobile_element".to_string(),
                    "sine".to_string(),
                    "retrotransposon".to_string(),
                    "repeat_region".to_string(),
                    "mapping_ambiguity_risk".to_string(),
                ],
                vec!["alu_like_sine_heuristic".to_string()],
                vec![
                    format!("candidate_length_bp={}", span.candidate_length_bp),
                    format!("tail_length_bp={}", span.tail_length_bp),
                    format!("arm_similarity={:.3}", span.arm_similarity),
                    "family_assignment=soft_heuristic".to_string(),
                ],
            );
        }

        evidence
    }

    fn push_construct_reasoning_objective_context_evidence(
        evidence: &mut Vec<DesignEvidence>,
        seq_id: &str,
        objective: &ConstructObjective,
        scope: EvidenceScope,
        label: String,
        rationale: String,
        host_profile_id: Option<String>,
        host_route_step_id: Option<String>,
        helper_profile_id: Option<String>,
        medium_condition_id: Option<String>,
        mut context_tags: Vec<String>,
        mut provenance_refs: Vec<String>,
        notes: Vec<String>,
    ) {
        context_tags.push("construct_objective".to_string());
        context_tags.sort();
        context_tags.dedup();
        provenance_refs.insert(0, objective.objective_id.clone());
        provenance_refs.sort();
        provenance_refs.dedup();
        evidence.push(DesignEvidence {
            seq_id: seq_id.to_string(),
            scope,
            start_0based: 0,
            end_0based_exclusive: 0,
            host_profile_id,
            host_route_step_id,
            helper_profile_id,
            medium_condition_id,
            role: ConstructRole::ContextBaggage,
            evidence_class: EvidenceClass::UserOverride,
            label,
            rationale,
            confidence: Some(1.0),
            context_tags,
            provenance_kind: "construct_objective".to_string(),
            provenance_refs,
            notes,
            ..DesignEvidence::default()
        });
    }

    fn build_construct_reasoning_objective_context_evidence(
        seq_id: &str,
        objective: &ConstructObjective,
    ) -> Vec<DesignEvidence> {
        let mut evidence: Vec<DesignEvidence> = vec![];

        if let Some(profile_id) = objective.propagation_host_profile_id.as_ref() {
            Self::push_construct_reasoning_objective_context_evidence(
                &mut evidence,
                seq_id,
                objective,
                EvidenceScope::HostProfile,
                format!("Propagation host: {profile_id}"),
                format!(
                    "Construct objective explicitly selects host profile '{profile_id}' for construct propagation."
                ),
                Some(profile_id.clone()),
                None,
                None,
                None,
                vec!["host_profile".to_string(), "propagation_host".to_string()],
                vec![profile_id.clone()],
                vec![],
            );
        }
        if let Some(profile_id) = objective.expression_host_profile_id.as_ref() {
            Self::push_construct_reasoning_objective_context_evidence(
                &mut evidence,
                seq_id,
                objective,
                EvidenceScope::HostProfile,
                format!("Expression host: {profile_id}"),
                format!(
                    "Construct objective explicitly selects host profile '{profile_id}' for downstream expression."
                ),
                Some(profile_id.clone()),
                None,
                None,
                None,
                vec!["expression_host".to_string(), "host_profile".to_string()],
                vec![profile_id.clone()],
                vec![],
            );
        }
        if let Some(helper_profile_id) = objective.helper_profile_id.as_ref() {
            Self::push_construct_reasoning_objective_context_evidence(
                &mut evidence,
                seq_id,
                objective,
                EvidenceScope::HelperProfile,
                format!("Helper profile: {helper_profile_id}"),
                format!(
                    "Construct objective explicitly selects helper/vector profile '{helper_profile_id}'."
                ),
                None,
                None,
                Some(helper_profile_id.clone()),
                None,
                vec!["helper_profile".to_string()],
                vec![helper_profile_id.clone()],
                vec![],
            );
        }
        for capture_plan in &objective.adapter_restriction_capture_plans {
            let mut context_tags = vec![
                "adapter_restriction_capture".to_string(),
                "linker".to_string(),
                capture_plan.adapter_style.as_str().to_string(),
                capture_plan.protection_mode.as_str().to_string(),
            ];
            if capture_plan.blunt_insert_required {
                context_tags.push("blunt_insert_required".to_string());
            }
            let mut provenance_refs = vec![
                capture_plan.capture_id.clone(),
                capture_plan.restriction_enzyme_name.clone(),
            ];
            provenance_refs.extend(capture_plan.extra_retrieval_enzyme_names.clone());
            Self::push_construct_reasoning_objective_context_evidence(
                &mut evidence,
                seq_id,
                objective,
                EvidenceScope::WholeConstruct,
                format!(
                    "Adapter capture: {} via {}",
                    capture_plan.capture_id, capture_plan.restriction_enzyme_name
                ),
                format!(
                    "Construct objective records a {} adapter/linker capture plan using restriction site '{}'.",
                    capture_plan.adapter_style.as_str(),
                    capture_plan.restriction_enzyme_name
                ),
                None,
                None,
                None,
                None,
                context_tags,
                provenance_refs,
                capture_plan.notes.clone(),
            );
            if !capture_plan.extra_retrieval_enzyme_names.is_empty() {
                Self::push_construct_reasoning_objective_context_evidence(
                    &mut evidence,
                    seq_id,
                    objective,
                    EvidenceScope::WholeConstruct,
                    format!(
                        "Adapter retrieval sites: {}",
                        capture_plan.extra_retrieval_enzyme_names.join(", ")
                    ),
                    format!(
                        "Construct objective records extra retrieval site(s) on adapter capture plan '{}'.",
                        capture_plan.capture_id
                    ),
                    None,
                    None,
                    None,
                    None,
                    vec![
                        "adapter_restriction_capture".to_string(),
                        "retrieval_sites".to_string(),
                        capture_plan.adapter_style.as_str().to_string(),
                    ],
                    capture_plan.extra_retrieval_enzyme_names.clone(),
                    capture_plan.notes.clone(),
                );
            }
        }
        for step in &objective.host_route {
            let rationale = if step.rationale.is_empty() {
                format!(
                    "Construct objective explicitly records a {} host transition step through '{}'.",
                    step.role.as_str(),
                    step.host_profile_id
                )
            } else {
                step.rationale.clone()
            };
            Self::push_construct_reasoning_objective_context_evidence(
                &mut evidence,
                seq_id,
                objective,
                EvidenceScope::HostTransition,
                format!(
                    "Host route ({}): {}",
                    step.role.as_str(),
                    step.host_profile_id
                ),
                rationale,
                Some(step.host_profile_id.clone()),
                Some(step.step_id.clone()),
                None,
                None,
                vec![
                    "host_transition".to_string(),
                    step.role.as_str().to_string(),
                ],
                vec![step.host_profile_id.clone(), step.step_id.clone()],
                step.notes.clone(),
            );
        }
        for condition in &objective.medium_conditions {
            Self::push_construct_reasoning_objective_context_evidence(
                &mut evidence,
                seq_id,
                objective,
                EvidenceScope::MediumCondition,
                format!("Medium condition: {condition}"),
                format!(
                    "Construct objective explicitly records medium/selection condition '{condition}'."
                ),
                None,
                None,
                None,
                Some(condition.clone()),
                vec!["medium_condition".to_string()],
                vec![condition.clone()],
                vec![],
            );
        }
        for trait_tag in &objective.required_host_traits {
            Self::push_construct_reasoning_objective_context_evidence(
                &mut evidence,
                seq_id,
                objective,
                EvidenceScope::WholeConstruct,
                format!("Required host trait: {trait_tag}"),
                format!("Construct objective explicitly requires host trait/tag '{trait_tag}'."),
                None,
                None,
                None,
                None,
                vec!["host_trait".to_string(), "required_host_trait".to_string()],
                vec![trait_tag.clone()],
                vec![],
            );
        }
        for trait_tag in &objective.forbidden_host_traits {
            Self::push_construct_reasoning_objective_context_evidence(
                &mut evidence,
                seq_id,
                objective,
                EvidenceScope::WholeConstruct,
                format!("Forbidden host trait: {trait_tag}"),
                format!("Construct objective explicitly forbids host trait/tag '{trait_tag}'."),
                None,
                None,
                None,
                None,
                vec!["forbidden_host_trait".to_string(), "host_trait".to_string()],
                vec![trait_tag.clone()],
                vec![],
            );
        }
        for (label, value, tag) in [
            (
                "Host species",
                objective.host_species.as_ref(),
                "host_species",
            ),
            ("Cell type", objective.cell_type.as_ref(), "cell_type"),
            ("Tissue", objective.tissue.as_ref(), "tissue"),
            ("Organelle", objective.organelle.as_ref(), "organelle"),
            (
                "Expression intent",
                objective.expression_intent.as_ref(),
                "expression_intent",
            ),
        ] {
            if let Some(value) = value {
                Self::push_construct_reasoning_objective_context_evidence(
                    &mut evidence,
                    seq_id,
                    objective,
                    EvidenceScope::WholeConstruct,
                    format!("{label}: {value}"),
                    format!("Construct objective explicitly records {tag} context '{value}'."),
                    None,
                    None,
                    None,
                    None,
                    vec![tag.to_string()],
                    vec![value.clone()],
                    vec![],
                );
            }
        }

        evidence
    }

    fn build_construct_reasoning_sequence_evidence_with_promoter_params(
        &self,
        seq_id: &str,
        dna: &DNAsequence,
        promoter_upstream_bp: usize,
        promoter_downstream_bp: usize,
    ) -> Vec<DesignEvidence> {
        let seq_len = dna.len();
        let mut evidence: Vec<DesignEvidence> = vec![];

        for (site_idx, site) in dna.restriction_enzyme_sites().iter().enumerate() {
            let Some((start_0based, end_0based_exclusive)) =
                site.recognition_bounds_0based(seq_len)
            else {
                continue;
            };
            let (forward_cut, reverse_cut) = site
                .strand_cut_positions_0based(seq_len)
                .unwrap_or((start_0based, start_0based));
            let geometry_tag = match site.enzyme.end_geometry() {
                crate::restriction_enzyme::RestrictionEndGeometry::Blunt => "blunt",
                crate::restriction_enzyme::RestrictionEndGeometry::FivePrimeOverhang(_) => {
                    "five_prime_overhang"
                }
                crate::restriction_enzyme::RestrictionEndGeometry::ThreePrimeOverhang(_) => {
                    "three_prime_overhang"
                }
            };
            evidence.push(DesignEvidence {
                evidence_id: format!(
                    "restriction_site_{}_{}_{}_{}",
                    Self::normalize_id_token(&site.enzyme.name),
                    start_0based,
                    end_0based_exclusive,
                    site_idx
                ),
                seq_id: seq_id.to_string(),
                start_0based,
                end_0based_exclusive,
                strand: Some(if site.forward_strand { "+" } else { "-" }.to_string()),
                role: ConstructRole::RestrictionSite,
                evidence_class: EvidenceClass::HardFact,
                label: site.enzyme.name.clone(),
                rationale: format!(
                    "Restriction site '{}' is a sequence-determined hard fact with {} cleavage geometry.",
                    site.enzyme.name, geometry_tag
                ),
                score: Some(1.0),
                confidence: Some(1.0),
                context_tags: vec![
                    "restriction_site".to_string(),
                    geometry_tag.to_string(),
                    site.enzyme.sequence.to_ascii_lowercase(),
                ],
                provenance_kind: "deterministic_restriction_site".to_string(),
                provenance_refs: vec![site.enzyme.name.clone()],
                notes: vec![format!(
                    "recognition={} forward_cut_0based={} reverse_cut_0based={}",
                    site.enzyme.sequence, forward_cut, reverse_cut
                )],
                ..DesignEvidence::default()
            });
        }

        for (feature_id, feature) in dna.features().iter().enumerate() {
            let Some(role) = Self::construct_reasoning_role_from_feature(feature) else {
                continue;
            };
            let evidence_class = Self::construct_reasoning_class_for_feature(feature, role);
            let mut ranges = vec![];
            collect_location_ranges_usize(&feature.location, &mut ranges);
            if ranges.is_empty() {
                continue;
            }
            let label = Self::feature_display_label(feature, feature_id);
            let rationale =
                Self::construct_reasoning_feature_rationale(feature, role, evidence_class);
            let context_tags = Self::construct_reasoning_feature_context_tags(feature);
            let strand = Some(
                if feature_is_reverse(feature) {
                    "-"
                } else {
                    "+"
                }
                .to_string(),
            );
            let total_parts = ranges.len();
            for (part_idx, (start_0based, end_0based_exclusive)) in ranges.into_iter().enumerate() {
                let mut part_label = label.clone();
                if total_parts > 1 {
                    part_label = format!(
                        "{} (part {}/{})",
                        part_label,
                        part_idx.saturating_add(1),
                        total_parts
                    );
                }
                evidence.push(DesignEvidence {
                    evidence_id: format!(
                        "feature_{}_{}_{}_{}_{}",
                        Self::normalize_id_token(&feature.kind.to_string()),
                        feature_id,
                        start_0based,
                        end_0based_exclusive,
                        part_idx
                    ),
                    seq_id: seq_id.to_string(),
                    start_0based,
                    end_0based_exclusive,
                    strand: strand.clone(),
                    role,
                    evidence_class,
                    label: part_label,
                    rationale: rationale.clone(),
                    confidence: match evidence_class {
                        EvidenceClass::HardFact => Some(1.0),
                        EvidenceClass::ReliableAnnotation => Some(0.9),
                        EvidenceClass::ContextEvidence => Some(0.7),
                        EvidenceClass::SoftHypothesis => Some(0.5),
                        EvidenceClass::UserOverride => Some(1.0),
                    },
                    context_tags: context_tags.clone(),
                    provenance_kind: "sequence_feature_annotation".to_string(),
                    provenance_refs: vec![feature.kind.to_string()],
                    ..DesignEvidence::default()
                });
            }
        }

        evidence.extend(self.build_construct_reasoning_generated_promoter_evidence(
            seq_id,
            dna,
            promoter_upstream_bp,
            promoter_downstream_bp,
        ));
        evidence.extend(Self::build_construct_reasoning_sequence_similarity_evidence(seq_id, dna));
        evidence
    }

    fn build_construct_reasoning_sequence_evidence(
        &self,
        seq_id: &str,
        dna: &DNAsequence,
    ) -> Vec<DesignEvidence> {
        self.build_construct_reasoning_sequence_evidence_with_promoter_params(
            seq_id,
            dna,
            DEFAULT_PROMOTER_WINDOW_UPSTREAM_BP,
            DEFAULT_PROMOTER_WINDOW_DOWNSTREAM_BP,
        )
    }

    fn construct_reasoning_ranges_for_feature(
        feature: &gb_io::seq::Feature,
    ) -> Vec<(usize, usize)> {
        let mut ranges = vec![];
        collect_location_ranges_usize(&feature.location, &mut ranges);
        ranges
    }

    fn construct_reasoning_first_range_for_feature(
        feature: &gb_io::seq::Feature,
    ) -> Option<(usize, usize)> {
        Self::construct_reasoning_ranges_for_feature(feature)
            .into_iter()
            .find(|(start, end)| end > start)
    }

    fn construct_reasoning_ranges_overlap(
        left_start: usize,
        left_end_exclusive: usize,
        right_start: usize,
        right_end_exclusive: usize,
    ) -> bool {
        left_start < right_end_exclusive && right_start < left_end_exclusive
    }

    fn construct_reasoning_normalized_dna_qualifier(
        feature: &gb_io::seq::Feature,
        key: &str,
    ) -> Option<String> {
        let value = Self::feature_qualifier_text(feature, key)?;
        let normalized = Self::normalize_dna_text(&value).to_ascii_uppercase();
        (!normalized.is_empty()).then_some(normalized)
    }

    fn construct_reasoning_feature_translation_table(feature: &gb_io::seq::Feature) -> usize {
        for key in ["transl_table", "translation_table"] {
            if let Some(value) = Self::feature_qualifier_text(feature, key)
                .and_then(|raw| raw.trim().parse::<usize>().ok())
                .filter(|value| *value > 0)
            {
                return value;
            }
        }
        1
    }

    fn construct_reasoning_feature_codon_start_offset(feature: &gb_io::seq::Feature) -> usize {
        if let Some(value) = Self::feature_qualifier_text(feature, "codon_start")
            .and_then(|raw| raw.trim().parse::<usize>().ok())
            .filter(|value| (1..=3).contains(value))
        {
            return value.saturating_sub(1);
        }
        if let Some(value) = Self::feature_qualifier_text(feature, "phase")
            .and_then(|raw| raw.trim().parse::<usize>().ok())
            .filter(|value| *value <= 2)
        {
            return value;
        }
        0
    }

    fn construct_reasoning_aa_display(aa: char) -> String {
        if aa == STOP_CODON {
            "*".to_string()
        } else if aa == UNKNOWN_CODON {
            "X".to_string()
        } else {
            aa.to_string()
        }
    }

    fn construct_reasoning_translate_codon(codon: &[u8; 3], translation_table: usize) -> char {
        crate::AMINO_ACIDS.codon2aa(*codon, Some(translation_table))
    }

    fn construct_reasoning_suggest_variant_assays_from_effect_tags(
        effect_tags: &[String],
    ) -> Vec<String> {
        let tags = effect_tags
            .iter()
            .map(String::as_str)
            .collect::<HashSet<_>>();
        let mut suggested_assay_ids = vec![];
        if tags.contains("promoter_variant_candidate") {
            suggested_assay_ids.push("allele_paired_promoter_luciferase_reporter".to_string());
        } else if tags.contains("enhancer_variant_candidate")
            || tags.contains("tfbs_overlap_candidate")
            || tags.contains("regulatory_variant_candidate")
        {
            suggested_assay_ids.push("allele_paired_regulatory_reporter".to_string());
        }
        if tags.contains("splice_variant_candidate") {
            suggested_assay_ids.push("minigene_splicing_assay".to_string());
        }
        if tags.contains("coding_variant_candidate")
            || tags.contains("missense_variant")
            || tags.contains("synonymous_variant")
            || tags.contains("nonsense_variant")
            || tags.contains("stop_lost_variant")
        {
            suggested_assay_ids.push("allele_paired_expression_compare".to_string());
        }
        if tags.contains("utr_variant_candidate")
            && !suggested_assay_ids
                .iter()
                .any(|value| value == "allele_paired_promoter_luciferase_reporter")
        {
            suggested_assay_ids.push("utr_reporter_translation_compare".to_string());
        }
        suggested_assay_ids.sort();
        suggested_assay_ids.dedup();
        suggested_assay_ids
    }

    fn construct_reasoning_predict_variant_coding_effects(
        dna: &DNAsequence,
        variant_feature: &gb_io::seq::Feature,
    ) -> Vec<ConstructVariantCodingPrediction> {
        let Some((variant_start_0based, variant_end_0based_exclusive)) =
            Self::construct_reasoning_first_range_for_feature(variant_feature)
        else {
            return vec![];
        };
        if variant_end_0based_exclusive.saturating_sub(variant_start_0based) != 1 {
            return vec![];
        }
        let Some(reference_allele) =
            Self::construct_reasoning_normalized_dna_qualifier(variant_feature, "vcf_ref")
        else {
            return vec![];
        };
        let Some(alternate_allele) =
            Self::construct_reasoning_normalized_dna_qualifier(variant_feature, "vcf_alt")
        else {
            return vec![];
        };
        if reference_allele.len() != 1 || alternate_allele.len() != 1 {
            return vec![];
        }

        let sequence = dna.get_forward_string().to_ascii_uppercase();
        let bytes = sequence.as_bytes();
        if variant_start_0based >= bytes.len() {
            return vec![];
        }
        let sequence_ref = bytes[variant_start_0based] as char;
        let expected_ref = reference_allele.as_bytes()[0] as char;
        if sequence_ref.to_ascii_uppercase() != expected_ref.to_ascii_uppercase() {
            return vec![];
        }

        let mut predictions = vec![];
        for (feature_id, feature) in dna.features().iter().enumerate() {
            if Self::construct_reasoning_role_from_feature(feature) != Some(ConstructRole::Cds) {
                continue;
            }
            let ranges = Self::construct_reasoning_ranges_for_feature(feature);
            if ranges.is_empty()
                || !ranges.iter().any(|(start, end)| {
                    Self::construct_reasoning_ranges_overlap(
                        variant_start_0based,
                        variant_end_0based_exclusive,
                        *start,
                        *end,
                    )
                })
            {
                continue;
            }

            let is_reverse = feature_is_reverse(feature);
            let mut ordered_ranges = ranges
                .iter()
                .copied()
                .filter(|(start, end)| end > start && *end <= bytes.len())
                .collect::<Vec<_>>();
            ordered_ranges.sort();
            if is_reverse {
                ordered_ranges.reverse();
            }
            if ordered_ranges.is_empty() {
                continue;
            }

            let mut oriented_bases: Vec<u8> = vec![];
            let mut oriented_positions: Vec<usize> = vec![];
            for (start, end) in &ordered_ranges {
                if is_reverse {
                    for pos in (*start..*end).rev() {
                        let rc = Self::reverse_complement_bytes(&[bytes[pos]]);
                        if let Some(base) = rc.first() {
                            oriented_bases.push(*base);
                            oriented_positions.push(pos);
                        }
                    }
                } else {
                    oriented_bases.extend_from_slice(&bytes[*start..*end]);
                    oriented_positions.extend(*start..*end);
                }
            }
            let codon_offset = Self::construct_reasoning_feature_codon_start_offset(feature);
            if codon_offset >= oriented_bases.len() {
                continue;
            }
            let trimmed_bases = &oriented_bases[codon_offset..];
            let trimmed_positions = &oriented_positions[codon_offset..];
            let Some(variant_index_in_cds) = trimmed_positions
                .iter()
                .position(|pos| *pos == variant_start_0based)
            else {
                continue;
            };
            let codon_start_0based = (variant_index_in_cds / 3) * 3;
            if codon_start_0based + 3 > trimmed_bases.len() {
                continue;
            }

            let ref_codon = [
                trimmed_bases[codon_start_0based],
                trimmed_bases[codon_start_0based + 1],
                trimmed_bases[codon_start_0based + 2],
            ];
            let mut alt_codon = ref_codon;
            let codon_pos = variant_index_in_cds % 3;
            let oriented_ref = if is_reverse {
                Self::reverse_complement(&reference_allele)
            } else {
                reference_allele.clone()
            };
            let oriented_alt = if is_reverse {
                Self::reverse_complement(&alternate_allele)
            } else {
                alternate_allele.clone()
            };
            if oriented_ref.len() != 1 || oriented_alt.len() != 1 {
                continue;
            }
            if ref_codon[codon_pos].to_ascii_uppercase()
                != oriented_ref.as_bytes()[0].to_ascii_uppercase()
            {
                continue;
            }
            alt_codon[codon_pos] = oriented_alt.as_bytes()[0].to_ascii_uppercase();

            let translation_table = Self::construct_reasoning_feature_translation_table(feature);
            let ref_aa = Self::construct_reasoning_translate_codon(&ref_codon, translation_table);
            let alt_aa = Self::construct_reasoning_translate_codon(&alt_codon, translation_table);
            let effect_tag = if ref_aa == alt_aa {
                "synonymous_variant"
            } else if alt_aa == STOP_CODON {
                "nonsense_variant"
            } else if ref_aa == STOP_CODON {
                "stop_lost_variant"
            } else if ref_aa == UNKNOWN_CODON || alt_aa == UNKNOWN_CODON {
                "coding_variant_candidate"
            } else {
                "missense_variant"
            };
            let transcript_label = Self::first_nonempty_feature_qualifier(
                feature,
                &["transcript_id", "transcript_name", "product"],
            );
            predictions.push(ConstructVariantCodingPrediction {
                cds_label: Self::feature_display_label(feature, feature_id),
                transcript_label,
                effect_tag: effect_tag.to_string(),
                codon_index_1based: codon_start_0based / 3 + 1,
                ref_codon: String::from_utf8_lossy(&ref_codon).to_string(),
                alt_codon: String::from_utf8_lossy(&alt_codon).to_string(),
                ref_amino_acid: Self::construct_reasoning_aa_display(ref_aa),
                alt_amino_acid: Self::construct_reasoning_aa_display(alt_aa),
                translation_table,
            });
        }

        predictions.sort_by(|left, right| {
            left.cds_label
                .to_ascii_lowercase()
                .cmp(&right.cds_label.to_ascii_lowercase())
                .then_with(|| left.codon_index_1based.cmp(&right.codon_index_1based))
        });
        predictions.dedup_by(|left, right| {
            left.cds_label == right.cds_label
                && left.codon_index_1based == right.codon_index_1based
                && left.effect_tag == right.effect_tag
        });
        predictions
    }

    fn construct_reasoning_variant_summary(
        dna: &DNAsequence,
        variant_feature: &gb_io::seq::Feature,
        variant_evidence: &DesignEvidence,
        evidence: &[DesignEvidence],
    ) -> ConstructVariantReasoningSummary {
        let overlapping = evidence
            .iter()
            .filter(|row| {
                row.scope == EvidenceScope::SequenceSpan
                    && row.evidence_id != variant_evidence.evidence_id
                    && row.role != ConstructRole::Variant
                    && Self::construct_reasoning_ranges_overlap(
                        variant_evidence.start_0based,
                        variant_evidence.end_0based_exclusive,
                        row.start_0based,
                        row.end_0based_exclusive,
                    )
            })
            .collect::<Vec<_>>();
        let has_role = |role| overlapping.iter().any(|row| row.role == role);
        let labels_for_role = |role| {
            overlapping
                .iter()
                .filter(|row| row.role == role)
                .map(|row| row.label.clone())
                .collect::<BTreeSet<_>>()
                .into_iter()
                .collect::<Vec<_>>()
        };

        let coding_predictions =
            Self::construct_reasoning_predict_variant_coding_effects(dna, variant_feature);

        let mut effect_tags = vec![];
        if has_role(ConstructRole::Promoter) {
            effect_tags.push("promoter_variant_candidate".to_string());
        }
        if has_role(ConstructRole::Enhancer) {
            effect_tags.push("enhancer_variant_candidate".to_string());
        }
        if has_role(ConstructRole::Tfbs) {
            effect_tags.push("tfbs_overlap_candidate".to_string());
        }
        if has_role(ConstructRole::Promoter) || has_role(ConstructRole::Enhancer) {
            effect_tags.push("regulatory_variant_candidate".to_string());
        }
        if has_role(ConstructRole::Promoter) && has_role(ConstructRole::Tfbs) {
            effect_tags.push("promoter_tfbs_regulatory_candidate".to_string());
        }
        if has_role(ConstructRole::Cds) || !coding_predictions.is_empty() {
            effect_tags.push("coding_variant_candidate".to_string());
        }
        if has_role(ConstructRole::Exon) && !has_role(ConstructRole::Cds) {
            effect_tags.push("exonic_variant_candidate".to_string());
        }
        if has_role(ConstructRole::SpliceBoundary) {
            effect_tags.push("splice_variant_candidate".to_string());
        }
        if has_role(ConstructRole::Utr5Prime) || has_role(ConstructRole::Utr3Prime) {
            effect_tags.push("utr_variant_candidate".to_string());
        }
        if has_role(ConstructRole::Gene) || has_role(ConstructRole::Transcript) {
            effect_tags.push("genic_variant_candidate".to_string());
        }
        effect_tags.extend(
            coding_predictions
                .iter()
                .map(|row| row.effect_tag.clone())
                .collect::<Vec<_>>(),
        );
        effect_tags.sort();
        effect_tags.dedup();

        let suggested_assay_ids =
            Self::construct_reasoning_suggest_variant_assays_from_effect_tags(&effect_tags);

        let overlapping_roles = overlapping
            .iter()
            .map(|row| row.role.as_str().to_string())
            .collect::<BTreeSet<_>>()
            .into_iter()
            .collect::<Vec<_>>();
        let overlapping_labels = overlapping
            .iter()
            .map(|row| row.label.clone())
            .collect::<BTreeSet<_>>()
            .into_iter()
            .collect::<Vec<_>>();
        let transcript_labels = labels_for_role(ConstructRole::Transcript);
        let general_transcript_tags = effect_tags
            .iter()
            .filter(|tag| {
                matches!(
                    tag.as_str(),
                    "promoter_variant_candidate"
                        | "enhancer_variant_candidate"
                        | "tfbs_overlap_candidate"
                        | "regulatory_variant_candidate"
                        | "promoter_tfbs_regulatory_candidate"
                        | "splice_variant_candidate"
                        | "utr_variant_candidate"
                        | "genic_variant_candidate"
                        | "exonic_variant_candidate"
                )
            })
            .cloned()
            .collect::<Vec<_>>();
        let mut transcript_effect_map: BTreeMap<
            String,
            (BTreeSet<String>, BTreeSet<String>, String),
        > = BTreeMap::new();
        for transcript_label in &transcript_labels {
            let assays = Self::construct_reasoning_suggest_variant_assays_from_effect_tags(
                &general_transcript_tags,
            );
            transcript_effect_map.insert(
                transcript_label.clone(),
                (
                    general_transcript_tags.iter().cloned().collect(),
                    assays.into_iter().collect(),
                    "shared_overlap_context".to_string(),
                ),
            );
        }
        for prediction in &coding_predictions {
            let transcript_label = prediction
                .transcript_label
                .clone()
                .unwrap_or_else(|| prediction.cds_label.clone());
            let entry = transcript_effect_map
                .entry(transcript_label)
                .or_insert_with(|| {
                    (
                        BTreeSet::new(),
                        BTreeSet::new(),
                        "coding_prediction".to_string(),
                    )
                });
            entry.0.insert("coding_variant_candidate".to_string());
            entry.0.insert(prediction.effect_tag.clone());
            for assay in Self::construct_reasoning_suggest_variant_assays_from_effect_tags(
                &entry.0.iter().cloned().collect::<Vec<_>>(),
            ) {
                entry.1.insert(assay);
            }
            if entry.2 != "coding_prediction" {
                entry.2 = "shared_overlap_and_coding_prediction".to_string();
            }
        }
        let transcript_effect_summaries = transcript_effect_map
            .into_iter()
            .map(
                |(transcript_label, (effect_tags, suggested_assay_ids, source))| {
                    ConstructTranscriptVariantEffectSummary {
                        transcript_label,
                        effect_tags: effect_tags.into_iter().collect(),
                        suggested_assay_ids: suggested_assay_ids.into_iter().collect(),
                        source,
                    }
                },
            )
            .collect::<Vec<_>>();
        let transcript_context_status = if transcript_effect_summaries.len() <= 1 {
            "single_or_unspecified"
        } else {
            let unique_effect_sets = transcript_effect_summaries
                .iter()
                .map(|row| row.effect_tags.join("|"))
                .collect::<BTreeSet<_>>();
            let unique_assay_sets = transcript_effect_summaries
                .iter()
                .map(|row| row.suggested_assay_ids.join("|"))
                .collect::<BTreeSet<_>>();
            if unique_effect_sets.len() > 1 || unique_assay_sets.len() > 1 {
                "multi_transcript_ambiguous"
            } else {
                "multi_transcript_consistent"
            }
        }
        .to_string();
        let mut effect_tags = effect_tags;
        if transcript_context_status == "multi_transcript_ambiguous" {
            effect_tags.push("transcript_context_ambiguous".to_string());
            effect_tags.sort();
            effect_tags.dedup();
        }

        ConstructVariantReasoningSummary {
            evidence_id: variant_evidence.evidence_id.clone(),
            label: variant_evidence.label.clone(),
            start_0based: variant_evidence.start_0based,
            end_0based_exclusive: variant_evidence.end_0based_exclusive,
            strand: variant_evidence.strand.clone(),
            variant_class: Self::feature_qualifier_text(variant_feature, "vcf_variant_class"),
            genomic_ref: Self::construct_reasoning_normalized_dna_qualifier(
                variant_feature,
                "vcf_ref",
            ),
            genomic_alt: Self::construct_reasoning_normalized_dna_qualifier(
                variant_feature,
                "vcf_alt",
            ),
            overlapping_roles,
            overlapping_labels,
            gene_labels: labels_for_role(ConstructRole::Gene),
            transcript_labels,
            effect_tags,
            suggested_assay_ids,
            transcript_context_status,
            transcript_effect_summaries,
            coding_predictions,
        }
    }

    fn construct_reasoning_collect_variant_summaries(
        dna: &DNAsequence,
        evidence: &[DesignEvidence],
        variant_evidence_rows: &[&DesignEvidence],
    ) -> Vec<ConstructVariantReasoningSummary> {
        if variant_evidence_rows.is_empty() {
            return vec![];
        }
        let mut variant_summaries = vec![];
        for (feature_id, feature) in dna.features().iter().enumerate() {
            if Self::construct_reasoning_role_from_feature(feature) != Some(ConstructRole::Variant)
            {
                continue;
            }
            let Some((start_0based, end_0based_exclusive)) =
                Self::construct_reasoning_first_range_for_feature(feature)
            else {
                continue;
            };
            let label = Self::feature_display_label(feature, feature_id);
            let label_token = Self::construct_reasoning_label_match_token(&label);
            let Some(variant_evidence) = variant_evidence_rows.iter().find(|row| {
                row.start_0based == start_0based
                    && row.end_0based_exclusive == end_0based_exclusive
                    && Self::construct_reasoning_label_match_token(&row.label) == label_token
            }) else {
                continue;
            };
            variant_summaries.push(Self::construct_reasoning_variant_summary(
                dna,
                feature,
                variant_evidence,
                evidence,
            ));
        }

        variant_summaries.sort_by(|left, right| {
            left.start_0based.cmp(&right.start_0based).then_with(|| {
                left.label
                    .to_ascii_lowercase()
                    .cmp(&right.label.to_ascii_lowercase())
            })
        });
        variant_summaries
    }

    fn construct_reasoning_label_match_token(label: &str) -> String {
        let base = label.split(" (part ").next().unwrap_or(label).trim();
        Self::normalize_id_token(base)
    }

    fn construct_reasoning_selection_rules() -> &'static [ConstructSelectionRule] {
        const EMPTY: &[&str] = &[];
        const PROLINE_MATCH_TOKENS: &[&str] = &["proa", "prob"];
        const PROLINE_CONDITION_SIGNALS: &[&str] = &["proline_omission"];
        const AMPICILLIN_FEATURE_TOKENS: &[&str] = &["ampr", "bla", "betalactamase"];
        const AMPICILLIN_CONDITION_SIGNALS: &[&str] =
            &["ampicillin_selection", "carbenicillin_selection"];
        const AMPICILLIN_HELPER_FUNCTIONS: &[&str] = &["ampicillin_selection"];
        const AMPICILLIN_HELPER_COMPONENT_KINDS: &[&str] = &["selectable_marker"];
        const AMPICILLIN_HELPER_COMPONENT_TAGS: &[&str] = &["selection"];
        const AMPICILLIN_HELPER_COMPONENT_ATTRIBUTES: &[&str] = &["ampicillin", "carbenicillin"];
        const RULES: &[ConstructSelectionRule] = &[
            ConstructSelectionRule {
                rule_id: "proline_auxotrophy_rescue",
                label: "Proline complementation",
                category: "complementation",
                feature_match_tokens: PROLINE_MATCH_TOKENS,
                condition_signal_tokens: PROLINE_CONDITION_SIGNALS,
                helper_offered_function_tokens: EMPTY,
                helper_component_kind_tokens: EMPTY,
                helper_component_tag_tokens: EMPTY,
                helper_component_attribute_terms: EMPTY,
            },
            ConstructSelectionRule {
                rule_id: "ampicillin_vector_selection",
                label: "Ampicillin vector selection",
                category: "selection",
                feature_match_tokens: AMPICILLIN_FEATURE_TOKENS,
                condition_signal_tokens: AMPICILLIN_CONDITION_SIGNALS,
                helper_offered_function_tokens: AMPICILLIN_HELPER_FUNCTIONS,
                helper_component_kind_tokens: AMPICILLIN_HELPER_COMPONENT_KINDS,
                helper_component_tag_tokens: AMPICILLIN_HELPER_COMPONENT_TAGS,
                helper_component_attribute_terms: AMPICILLIN_HELPER_COMPONENT_ATTRIBUTES,
            },
        ];
        RULES
    }

    fn construct_reasoning_condition_terms_matching(
        condition: &str,
        terms: &[&str],
    ) -> Vec<String> {
        let lower = condition.to_ascii_lowercase();
        let mut matches = terms
            .iter()
            .filter(|term| lower.contains(**term))
            .map(|term| (*term).to_string())
            .collect::<Vec<_>>();
        matches.sort();
        matches.dedup();
        matches
    }

    fn construct_reasoning_selection_status_from_match(
        has_medium_match: bool,
        has_candidates: bool,
    ) -> &'static str {
        if has_medium_match && has_candidates {
            "supported"
        } else if has_medium_match {
            "review_needed"
        } else if has_candidates {
            "candidates_detected"
        } else {
            "context_recorded"
        }
    }

    fn construct_reasoning_push_condition_signal(
        signals: &mut Vec<ConstructConditionSignal>,
        seen: &mut HashSet<String>,
        category: &str,
        token: &str,
        label: &str,
        source_condition: &str,
        matched_terms: Vec<String>,
        temperature_celsius: Option<f64>,
        notes: Vec<String>,
    ) {
        let key = format!(
            "{}|{}|{}|{}",
            category,
            token,
            source_condition,
            temperature_celsius
                .map(|value| format!("{value:.2}"))
                .unwrap_or_default()
        );
        if !seen.insert(key) {
            return;
        }
        signals.push(ConstructConditionSignal {
            category: category.to_string(),
            token: token.to_string(),
            label: label.to_string(),
            source_condition: source_condition.to_string(),
            matched_terms,
            temperature_celsius,
            notes,
        });
    }

    fn construct_reasoning_extract_temperature_celsius(condition: &str) -> Option<f64> {
        let regex = Regex::new(r"(?i)(-?\d+(?:\.\d+)?)\s*°?\s*c\b").ok()?;
        let captures = regex.captures(condition)?;
        captures.get(1)?.as_str().parse::<f64>().ok()
    }

    fn construct_reasoning_temperature_signal_token(temp_celsius: f64) -> String {
        if (temp_celsius.fract()).abs() < 0.01 {
            format!("temperature_{:.0}c", temp_celsius)
        } else {
            format!(
                "temperature_{}c",
                temp_celsius.to_string().replace('.', "_")
            )
        }
    }

    fn construct_reasoning_interpret_medium_conditions(
        objective: &ConstructObjective,
    ) -> Vec<ConstructConditionSignal> {
        let mut signals = vec![];
        let mut seen: HashSet<String> = HashSet::new();
        for condition in &objective.medium_conditions {
            let proline_terms = Self::construct_reasoning_condition_terms_matching(
                condition,
                &[
                    "proline-free",
                    "without proline",
                    "no proline",
                    "-pro",
                    "dropout proline",
                ],
            );
            if !proline_terms.is_empty() {
                Self::construct_reasoning_push_condition_signal(
                    &mut signals,
                    &mut seen,
                    "nutrient_condition",
                    "proline_omission",
                    "Proline omitted",
                    condition,
                    proline_terms,
                    None,
                    vec![],
                );
            }

            for (token, label, terms) in [
                (
                    "ampicillin_selection",
                    "Ampicillin selection",
                    vec!["ampicillin"],
                ),
                (
                    "carbenicillin_selection",
                    "Carbenicillin selection",
                    vec!["carbenicillin"],
                ),
            ] {
                let matched_terms =
                    Self::construct_reasoning_condition_terms_matching(condition, &terms);
                if matched_terms.is_empty() {
                    continue;
                }
                Self::construct_reasoning_push_condition_signal(
                    &mut signals,
                    &mut seen,
                    "selection_agent",
                    token,
                    label,
                    condition,
                    matched_terms,
                    None,
                    vec![],
                );
            }

            for (token, label, terms) in [
                ("heat_shock", "Heat shock", vec!["heat shock"]),
                ("cold_shock", "Cold shock", vec!["cold shock"]),
            ] {
                let matched_terms =
                    Self::construct_reasoning_condition_terms_matching(condition, &terms);
                if matched_terms.is_empty() {
                    continue;
                }
                Self::construct_reasoning_push_condition_signal(
                    &mut signals,
                    &mut seen,
                    "temperature_process",
                    token,
                    label,
                    condition,
                    matched_terms,
                    None,
                    vec![],
                );
            }

            if let Some(temp_celsius) =
                Self::construct_reasoning_extract_temperature_celsius(condition)
            {
                let exact_token = Self::construct_reasoning_temperature_signal_token(temp_celsius);
                Self::construct_reasoning_push_condition_signal(
                    &mut signals,
                    &mut seen,
                    "temperature",
                    &exact_token,
                    &format!("{temp_celsius:.0} C"),
                    condition,
                    vec![format!("{temp_celsius:.0}c")],
                    Some(temp_celsius),
                    vec!["temperature parsed from free-text condition".to_string()],
                );
                if temp_celsius >= 40.0 {
                    Self::construct_reasoning_push_condition_signal(
                        &mut signals,
                        &mut seen,
                        "temperature",
                        "high_temperature",
                        "High temperature",
                        condition,
                        vec![format!("{temp_celsius:.0}c")],
                        Some(temp_celsius),
                        vec![],
                    );
                }
                if temp_celsius <= 20.0 {
                    Self::construct_reasoning_push_condition_signal(
                        &mut signals,
                        &mut seen,
                        "temperature",
                        "low_temperature",
                        "Low temperature",
                        condition,
                        vec![format!("{temp_celsius:.0}c")],
                        Some(temp_celsius),
                        vec![],
                    );
                }
            }
        }

        signals.sort_by(|left, right| {
            (
                left.category.as_str(),
                left.token.as_str(),
                left.source_condition.as_str(),
            )
                .cmp(&(
                    right.category.as_str(),
                    right.token.as_str(),
                    right.source_condition.as_str(),
                ))
        });
        signals
    }

    fn construct_reasoning_value_matches_tokens(value: &str, tokens: &[&str]) -> bool {
        if tokens.is_empty() {
            return false;
        }
        let token = Self::normalize_id_token(value);
        tokens.contains(&token.as_str())
    }

    fn construct_reasoning_value_contains_any_term(value: &str, terms: &[&str]) -> Vec<String> {
        if terms.is_empty() {
            return vec![];
        }
        let lower = value.to_ascii_lowercase();
        let mut matches = terms
            .iter()
            .filter(|term| lower.contains(**term))
            .map(|term| (*term).to_string())
            .collect::<Vec<_>>();
        matches.sort();
        matches.dedup();
        matches
    }

    fn construct_reasoning_build_selection_rule_match(
        rule: &ConstructSelectionRule,
        condition_signals: &[ConstructConditionSignal],
        evidence: &[DesignEvidence],
        helper_interpretation: Option<&HelperConstructInterpretation>,
    ) -> Option<ConstructSelectionRuleMatch> {
        let mut matched_medium_conditions = vec![];
        let mut matched_medium_terms = vec![];
        let mut matched_condition_signal_tokens = vec![];
        for signal in condition_signals.iter().filter(|signal| {
            rule.condition_signal_tokens
                .contains(&signal.token.as_str())
        }) {
            matched_medium_conditions.push(signal.source_condition.clone());
            matched_medium_terms.extend(signal.matched_terms.clone());
            matched_condition_signal_tokens.push(signal.token.clone());
        }
        matched_medium_conditions.sort();
        matched_medium_conditions.dedup();
        matched_medium_terms.sort();
        matched_medium_terms.dedup();
        matched_condition_signal_tokens.sort();
        matched_condition_signal_tokens.dedup();

        let mut candidate_labels = vec![];
        let mut candidate_evidence_ids = vec![];
        let mut candidate_match_tokens = vec![];
        for row in evidence.iter().filter(|row| {
            row.scope == EvidenceScope::SequenceSpan
                && matches!(row.role, ConstructRole::Gene | ConstructRole::Cds)
        }) {
            let token = Self::construct_reasoning_label_match_token(&row.label);
            if !rule.feature_match_tokens.contains(&token.as_str()) {
                continue;
            }
            candidate_labels.push(row.label.clone());
            candidate_evidence_ids.push(row.evidence_id.clone());
            candidate_match_tokens.push(token);
        }
        candidate_labels.sort();
        candidate_labels.dedup();
        candidate_evidence_ids.sort();
        candidate_evidence_ids.dedup();
        candidate_match_tokens.sort();
        candidate_match_tokens.dedup();

        let mut helper_offered_function_matches = vec![];
        let mut helper_component_ids = vec![];
        let mut helper_component_labels = vec![];
        let mut helper_component_kind_matches = vec![];
        let mut helper_component_tag_matches = vec![];
        let mut helper_component_attribute_matches = vec![];
        if let Some(helper_interpretation) = helper_interpretation {
            for function in &helper_interpretation.offered_functions {
                if Self::construct_reasoning_value_matches_tokens(
                    function,
                    rule.helper_offered_function_tokens,
                ) {
                    helper_offered_function_matches.push(function.clone());
                }
            }
            for component in &helper_interpretation.components {
                let kind_matched = Self::construct_reasoning_value_matches_tokens(
                    &component.kind,
                    rule.helper_component_kind_tokens,
                );
                let tag_matches = component
                    .tags
                    .iter()
                    .filter(|tag| {
                        Self::construct_reasoning_value_matches_tokens(
                            tag,
                            rule.helper_component_tag_tokens,
                        )
                    })
                    .cloned()
                    .collect::<Vec<_>>();
                let mut attribute_matches = component
                    .attributes
                    .iter()
                    .flat_map(|(_, value)| {
                        Self::construct_reasoning_value_contains_any_term(
                            value,
                            rule.helper_component_attribute_terms,
                        )
                    })
                    .collect::<Vec<_>>();
                attribute_matches.sort();
                attribute_matches.dedup();

                if !(kind_matched || !tag_matches.is_empty() || !attribute_matches.is_empty()) {
                    continue;
                }

                helper_component_ids.push(component.id.clone());
                helper_component_labels.push(
                    component
                        .label
                        .clone()
                        .unwrap_or_else(|| component.id.clone()),
                );
                if kind_matched {
                    helper_component_kind_matches.push(component.kind.clone());
                }
                helper_component_tag_matches.extend(tag_matches);
                helper_component_attribute_matches.extend(attribute_matches);
            }
        }
        helper_offered_function_matches.sort();
        helper_offered_function_matches.dedup();
        helper_component_ids.sort();
        helper_component_ids.dedup();
        helper_component_labels.sort();
        helper_component_labels.dedup();
        helper_component_kind_matches.sort();
        helper_component_kind_matches.dedup();
        helper_component_tag_matches.sort();
        helper_component_tag_matches.dedup();
        helper_component_attribute_matches.sort();
        helper_component_attribute_matches.dedup();

        if let Some(helper_interpretation) = helper_interpretation {
            for label in &helper_component_labels {
                candidate_labels.push(format!("{}: {}", helper_interpretation.helper_id, label));
            }
        }
        candidate_labels.sort();
        candidate_labels.dedup();

        if matched_medium_conditions.is_empty()
            && candidate_labels.is_empty()
            && helper_offered_function_matches.is_empty()
        {
            return None;
        }

        Some(ConstructSelectionRuleMatch {
            rule_id: rule.rule_id.to_string(),
            label: rule.label.to_string(),
            category: rule.category.to_string(),
            status: Self::construct_reasoning_selection_status_from_match(
                !matched_medium_conditions.is_empty(),
                !candidate_labels.is_empty() || !helper_offered_function_matches.is_empty(),
            )
            .to_string(),
            matched_medium_conditions,
            matched_medium_terms,
            matched_condition_signal_tokens,
            candidate_labels,
            candidate_evidence_ids,
            candidate_match_tokens,
            helper_profile_id: helper_interpretation.map(|row| row.helper_id.clone()),
            helper_offered_function_matches,
            helper_component_ids,
            helper_component_labels,
            helper_component_kind_matches,
            helper_component_tag_matches,
            helper_component_attribute_matches,
        })
    }

    fn construct_reasoning_count_pattern_with_n(sequence: &[u8], pattern: &[u8]) -> usize {
        if pattern.is_empty() || sequence.len() < pattern.len() {
            return 0;
        }
        let mut count = 0usize;
        for start in 0..=sequence.len().saturating_sub(pattern.len()) {
            let matched = pattern.iter().enumerate().all(|(offset, expected)| {
                let observed = sequence[start + offset].to_ascii_uppercase();
                *expected == b'N' || observed == expected.to_ascii_uppercase()
            });
            if matched {
                count = count.saturating_add(1);
            }
        }
        count
    }

    fn construct_reasoning_scan_sequence_restriction_methylation_patterns(
        dna: &DNAsequence,
    ) -> ConstructSequenceRestrictionMethylationPatterns {
        let sequence = dna.forward_bytes();
        let mut dam_mode = crate::methylation_sites::MethylationMode::default();
        dam_mode.set_dam(true);
        let dam_site_count =
            crate::methylation_sites::MethylationSites::new_from_sequence(sequence, dam_mode)
                .sites()
                .len();
        let mut dcm_mode = crate::methylation_sites::MethylationMode::default();
        dcm_mode.set_dcm(true);
        let dcm_site_count =
            crate::methylation_sites::MethylationSites::new_from_sequence(sequence, dcm_mode)
                .sites()
                .len();
        let eco_k_target_site_count =
            Self::construct_reasoning_count_pattern_with_n(sequence, b"AACNNNNNNGTGC")
                .saturating_add(Self::construct_reasoning_count_pattern_with_n(
                    sequence,
                    b"GCACNNNNNNGTT",
                ));
        let mut present_pattern_tokens = vec![];
        if dam_site_count > 0 {
            present_pattern_tokens.push("dam_site".to_string());
        }
        if dcm_site_count > 0 {
            present_pattern_tokens.push("dcm_site".to_string());
        }
        if eco_k_target_site_count > 0 {
            present_pattern_tokens.push("eco_k_target_like_site".to_string());
        }
        ConstructSequenceRestrictionMethylationPatterns {
            dam_site_count,
            dcm_site_count,
            eco_k_target_site_count,
            present_pattern_tokens,
        }
    }

    fn construct_reasoning_host_text_mentions_positive(text: &str, token: &str) -> bool {
        let compact = text
            .chars()
            .filter(|ch| !ch.is_whitespace())
            .collect::<String>()
            .to_ascii_lowercase();
        compact.contains(&format!("{token}+"))
            || compact.contains(&format!("{token}plus"))
            || compact.contains(&format!("{token}positive"))
    }

    fn construct_reasoning_host_text_mentions_negative(text: &str, token: &str) -> bool {
        let compact = text
            .chars()
            .filter(|ch| !ch.is_whitespace())
            .collect::<String>()
            .to_ascii_lowercase();
        compact.contains(&format!("{token}-"))
            || compact.contains(&format!("{token}minus"))
            || compact.contains(&format!("{token}negative"))
    }

    fn construct_reasoning_collect_host_route_trait_tokens(step: &HostRouteStep) -> Vec<String> {
        let text = std::iter::once(step.host_profile_id.as_str())
            .chain(std::iter::once(step.rationale.as_str()))
            .chain(step.notes.iter().map(String::as_str))
            .collect::<Vec<_>>()
            .join(" ");
        let compact = text
            .chars()
            .filter(|ch| !ch.is_whitespace())
            .collect::<String>()
            .to_ascii_lowercase();
        let mut tokens = BTreeSet::new();
        for token in [
            "hsdr", "hsdm", "hsds", "dam", "dcm", "mdrs", "mcr", "mcra", "mcrbc", "mrr",
        ] {
            if Self::construct_reasoning_host_text_mentions_positive(&text, token) {
                tokens.insert(format!("{token}_present"));
            }
            if Self::construct_reasoning_host_text_mentions_negative(&text, token) {
                tokens.insert(format!("{token}_absent"));
            }
        }
        if compact.contains("hsdr-") || compact.contains("hsdrminus") {
            if compact.contains("m+") || compact.contains("mplus") {
                tokens.insert("hsdm_present".to_string());
            }
            if compact.contains("m-") || compact.contains("mminus") {
                tokens.insert("hsdm_absent".to_string());
            }
        }
        tokens.into_iter().collect()
    }

    fn construct_reasoning_build_host_route_restriction_methylation_steps(
        objective: &ConstructObjective,
    ) -> Vec<ConstructHostRestrictionMethylationStep> {
        objective
            .host_route
            .iter()
            .map(|step| ConstructHostRestrictionMethylationStep {
                step_id: step.step_id.clone(),
                host_profile_id: step.host_profile_id.clone(),
                role: step.role.as_str().to_string(),
                trait_tokens: Self::construct_reasoning_collect_host_route_trait_tokens(step),
                notes: step.notes.clone(),
            })
            .collect()
    }

    fn construct_reasoning_route_step_has_any_trait(
        step: &ConstructHostRestrictionMethylationStep,
        tokens: &[&str],
    ) -> bool {
        tokens
            .iter()
            .any(|token| step.trait_tokens.iter().any(|value| value == token))
    }

    fn construct_reasoning_detect_route_restriction_methylation_conflicts(
        steps: &[ConstructHostRestrictionMethylationStep],
        patterns: &ConstructSequenceRestrictionMethylationPatterns,
    ) -> Vec<ConstructRestrictionMethylationConflict> {
        let mut conflicts = vec![];
        for window in steps.windows(2) {
            let [upstream, downstream] = match window {
                [upstream, downstream] => [upstream, downstream],
                _ => continue,
            };
            let upstream_has_hsd_methylase =
                Self::construct_reasoning_route_step_has_any_trait(upstream, &["hsdm_present"]);
            let upstream_has_dam_dcm = Self::construct_reasoning_route_step_has_any_trait(
                upstream,
                &["dam_present", "dcm_present"],
            );
            let downstream_has_hsd_restriction =
                Self::construct_reasoning_route_step_has_any_trait(downstream, &["hsdr_present"]);
            let downstream_has_methylated_dna_restriction =
                Self::construct_reasoning_route_step_has_any_trait(
                    downstream,
                    &[
                        "mdrs_present",
                        "mcr_present",
                        "mcra_present",
                        "mcrbc_present",
                        "mrr_present",
                    ],
                );

            if upstream_has_hsd_methylase
                && downstream_has_hsd_restriction
                && patterns.eco_k_target_site_count > 0
            {
                conflicts.push(ConstructRestrictionMethylationConflict {
                    conflict_id: "ecoki_hsdr_transition".to_string(),
                    label: "HsdR/HsdM host transition review".to_string(),
                    upstream_step_id: upstream.step_id.clone(),
                    downstream_step_id: downstream.step_id.clone(),
                    rationale: format!(
                        "Upstream step '{}' carries HsdM-like methylase context while downstream step '{}' carries HsdR-like restriction context, and the sequence contains {} EcoK target-like site(s). Direct transfer across that route should stay inspectable.",
                        upstream.step_id,
                        downstream.step_id,
                        patterns.eco_k_target_site_count
                    ),
                });
            }

            if upstream_has_dam_dcm
                && downstream_has_methylated_dna_restriction
                && (patterns.dam_site_count > 0 || patterns.dcm_site_count > 0)
            {
                conflicts.push(ConstructRestrictionMethylationConflict {
                    conflict_id: "dam_dcm_methylated_dna_transition".to_string(),
                    label: "Dam/Dcm methylation-dependent restriction review".to_string(),
                    upstream_step_id: upstream.step_id.clone(),
                    downstream_step_id: downstream.step_id.clone(),
                    rationale: format!(
                        "Upstream step '{}' carries Dam/Dcm-like methylation context while downstream step '{}' carries methylation-dependent restriction context, and the sequence contains {} Dam site(s) plus {} Dcm site(s).",
                        upstream.step_id,
                        downstream.step_id,
                        patterns.dam_site_count,
                        patterns.dcm_site_count
                    ),
                });
            }
        }
        conflicts
    }

    fn construct_reasoning_find_restriction_enzyme_by_name(
        enzyme_name: &str,
    ) -> Option<RestrictionEnzyme> {
        let target = enzyme_name.trim();
        if target.is_empty() {
            return None;
        }
        active_restriction_enzymes()
            .into_iter()
            .find(|enzyme| enzyme.name.eq_ignore_ascii_case(target))
    }

    fn construct_reasoning_internal_restriction_site_ranges(
        dna: &DNAsequence,
        enzyme: &RestrictionEnzyme,
    ) -> Vec<[usize; 2]> {
        enzyme
            .get_sites(dna, None)
            .into_iter()
            .filter_map(|site| {
                site.recognition_bounds_0based(dna.len())
                    .map(|(start, end)| [start, end])
            })
            .collect()
    }

    fn construct_reasoning_adapter_capture_methylation_review_note() -> &'static str {
        "A possible rescue is methylation-based protection if a compatible methylation system blocks cleavage at that motif, but GENtle's enzyme-specific methylation knowledge base is still incomplete."
    }

    fn construct_reasoning_summarize_adapter_restriction_capture_plans(
        dna: &DNAsequence,
        objective: &ConstructObjective,
    ) -> (
        Vec<ConstructAdapterRestrictionCapturePlanSummary>,
        Vec<ConstructAdapterRestrictionCaptureReviewItem>,
        Vec<String>,
    ) {
        let mut summaries = vec![];
        let mut review_items = vec![];
        let mut derived_preferred_routine_families = BTreeSet::new();

        for plan in &objective.adapter_restriction_capture_plans {
            derived_preferred_routine_families.insert("restriction".to_string());
            let mut seen_named_motifs = BTreeSet::new();
            let named_motifs =
                std::iter::once((plan.restriction_enzyme_name.clone(), "capture".to_string()))
                    .chain(
                        plan.extra_retrieval_enzyme_names
                            .iter()
                            .cloned()
                            .map(|enzyme_name| (enzyme_name, "retrieval".to_string())),
                    )
                    .filter_map(|(enzyme_name, motif_role)| {
                        let trimmed = enzyme_name.trim();
                        if trimmed.is_empty() {
                            return None;
                        }
                        let dedupe_key = trimmed.to_ascii_lowercase();
                        seen_named_motifs
                            .insert(dedupe_key)
                            .then(|| (trimmed.to_string(), motif_role))
                    })
                    .collect::<Vec<_>>();

            let mut capture_restriction_enzyme_name = plan.restriction_enzyme_name.clone();
            let mut capture_enzyme_resolution_status = "not_found".to_string();
            let mut capture_site_geometry = None;
            let mut capture_internal_site_count = 0usize;
            let mut capture_internal_site_status = "enzyme_not_found".to_string();
            let mut capture_internal_site_ranges_0based = vec![];
            let mut motif_summaries = vec![];
            let mut resolved_motif_count = 0usize;
            let mut every_named_motif_present_on_insert = !named_motifs.is_empty();

            for (enzyme_name, motif_role) in &named_motifs {
                let Some(enzyme) =
                    Self::construct_reasoning_find_restriction_enzyme_by_name(enzyme_name)
                else {
                    motif_summaries.push(ConstructAdapterRestrictionCaptureMotifSummary {
                        enzyme_name: enzyme_name.clone(),
                        motif_role: motif_role.clone(),
                        resolution_status: "not_found".to_string(),
                        site_geometry: None,
                        internal_site_count: 0,
                        internal_site_status: "enzyme_not_found".to_string(),
                        internal_site_ranges_0based: vec![],
                    });
                    every_named_motif_present_on_insert = false;
                    let (issue_id, label) = if motif_role == "capture" {
                        (
                            "capture_enzyme_not_found",
                            "Adapter capture enzyme needs review",
                        )
                    } else {
                        (
                            "retrieval_enzyme_not_found",
                            "Adapter retrieval enzyme needs review",
                        )
                    };
                    review_items.push(ConstructAdapterRestrictionCaptureReviewItem {
                        capture_id: plan.capture_id.clone(),
                        issue_id: issue_id.to_string(),
                        label: label.to_string(),
                        rationale: format!(
                            "Adapter/linker capture plan '{}' references {} restriction enzyme '{}', but the active restriction catalog could not resolve that enzyme name.",
                            plan.capture_id,
                            motif_role,
                            enzyme_name
                        ),
                    });
                    continue;
                };

                resolved_motif_count = resolved_motif_count.saturating_add(1);
                let internal_site_ranges =
                    Self::construct_reasoning_internal_restriction_site_ranges(dna, &enzyme);
                let internal_site_count = internal_site_ranges.len();
                let internal_site_status = if internal_site_count == 0 {
                    "no_internal_site_conflict"
                } else if plan.protection_mode == AdapterCaptureProtectionMode::InsertMethylation {
                    "methylation_protection_requested"
                } else {
                    "internal_site_conflict"
                };
                if internal_site_count == 0 {
                    every_named_motif_present_on_insert = false;
                }

                if motif_role == "capture" {
                    capture_restriction_enzyme_name = enzyme.name.clone();
                    capture_enzyme_resolution_status = "resolved".to_string();
                    capture_site_geometry = Some(enzyme.end_geometry().kind_label().to_string());
                    capture_internal_site_count = internal_site_count;
                    capture_internal_site_status = internal_site_status.to_string();
                    capture_internal_site_ranges_0based = internal_site_ranges.clone();
                }

                motif_summaries.push(ConstructAdapterRestrictionCaptureMotifSummary {
                    enzyme_name: enzyme.name.clone(),
                    motif_role: motif_role.clone(),
                    resolution_status: "resolved".to_string(),
                    site_geometry: Some(enzyme.end_geometry().kind_label().to_string()),
                    internal_site_count,
                    internal_site_status: internal_site_status.to_string(),
                    internal_site_ranges_0based: internal_site_ranges.clone(),
                });

                if internal_site_count == 0 {
                    continue;
                }

                let methylation_note =
                    Self::construct_reasoning_adapter_capture_methylation_review_note();
                if plan.protection_mode == AdapterCaptureProtectionMode::InsertMethylation {
                    let (issue_id, label) = if motif_role == "capture" {
                        (
                            "methylation_protection_requires_enzyme_specific_review",
                            "Insert methylation protection requires review",
                        )
                    } else {
                        (
                            "retrieval_motif_methylation_protection_requires_review",
                            "Adapter retrieval motif protection requires review",
                        )
                    };
                    review_items.push(ConstructAdapterRestrictionCaptureReviewItem {
                        capture_id: plan.capture_id.clone(),
                        issue_id: issue_id.to_string(),
                        label: label.to_string(),
                        rationale: format!(
                            "Adapter/linker capture plan '{}' reuses {} restriction site '{}' that already occurs {} time(s) on the insert. Planned insert methylation keeps that protection strategy explicit, but enzyme-specific methylation sensitivity still needs review before assuming adapter-only cleavage. {methylation_note}",
                            plan.capture_id,
                            motif_role,
                            enzyme.name,
                            internal_site_count
                        ),
                    });
                } else {
                    let (issue_id, label) = if motif_role == "capture" {
                        (
                            "internal_capture_site_conflict",
                            "Internal adapter capture site conflict",
                        )
                    } else {
                        (
                            "internal_retrieval_site_conflict",
                            "Internal adapter retrieval site conflict",
                        )
                    };
                    review_items.push(ConstructAdapterRestrictionCaptureReviewItem {
                        capture_id: plan.capture_id.clone(),
                        issue_id: issue_id.to_string(),
                        label: label.to_string(),
                        rationale: format!(
                            "Adapter/linker capture plan '{}' uses {} restriction site '{}', and that site already occurs {} time(s) on the insert without an explicit protection mode. {methylation_note}",
                            plan.capture_id,
                            motif_role,
                            enzyme.name,
                            internal_site_count
                        ),
                    });
                }
            }

            let all_named_motifs_present_on_insert = named_motifs.len() > 1
                && resolved_motif_count == named_motifs.len()
                && every_named_motif_present_on_insert;
            if all_named_motifs_present_on_insert {
                review_items.push(ConstructAdapterRestrictionCaptureReviewItem {
                    capture_id: plan.capture_id.clone(),
                    issue_id: "all_adapter_motifs_already_present_on_insert".to_string(),
                    label: "All adapter motifs already occur on the insert".to_string(),
                    rationale: format!(
                        "Adapter/linker capture plan '{}' names {} motif(s), and every resolved adapter motif already occurs internally on the insert. That means the adapter may contribute little retrieval discrimination after ligation/digest. {}",
                        plan.capture_id,
                        named_motifs.len(),
                        Self::construct_reasoning_adapter_capture_methylation_review_note()
                    ),
                });
            }

            summaries.push(ConstructAdapterRestrictionCapturePlanSummary {
                capture_id: plan.capture_id.clone(),
                restriction_enzyme_name: capture_restriction_enzyme_name,
                enzyme_resolution_status: capture_enzyme_resolution_status,
                adapter_style: plan.adapter_style.as_str().to_string(),
                blunt_insert_required: plan.blunt_insert_required,
                protection_mode: plan.protection_mode.as_str().to_string(),
                extra_retrieval_enzyme_names: plan.extra_retrieval_enzyme_names.clone(),
                capture_site_geometry,
                internal_site_count: capture_internal_site_count,
                internal_site_status: capture_internal_site_status,
                internal_site_ranges_0based: capture_internal_site_ranges_0based,
                named_motif_count: named_motifs.len(),
                resolved_motif_count,
                all_named_motifs_present_on_insert,
                motif_summaries,
                notes: plan.notes.clone(),
            });
        }

        (
            summaries,
            review_items,
            derived_preferred_routine_families.into_iter().collect(),
        )
    }

    fn construct_reasoning_evidence_ids_matching<F>(
        evidence: &[DesignEvidence],
        mut predicate: F,
    ) -> Vec<String>
    where
        F: FnMut(&DesignEvidence) -> bool,
    {
        let mut ids = evidence
            .iter()
            .filter(|row| predicate(row))
            .map(|row| row.evidence_id.clone())
            .collect::<Vec<_>>();
        ids.sort();
        ids.dedup();
        ids
    }

    fn construct_reasoning_build_fact(
        fact_id: &str,
        fact_type: &str,
        label: String,
        rationale: String,
        based_on_evidence_ids: Vec<String>,
        value_json: serde_json::Value,
    ) -> DesignFact {
        DesignFact {
            fact_id: fact_id.to_string(),
            fact_type: fact_type.to_string(),
            label,
            rationale,
            based_on_evidence_ids,
            value_json,
            confidence: Some(1.0),
            editable_status: EditableStatus::Draft,
            ..DesignFact::default()
        }
    }

    fn construct_reasoning_build_decision(
        decision_id: &str,
        decision_type: &str,
        title: String,
        rationale: String,
        input_evidence_ids: Vec<String>,
        output_fact_ids: Vec<String>,
        parameters_json: serde_json::Value,
    ) -> DesignDecisionNode {
        DesignDecisionNode {
            decision_id: decision_id.to_string(),
            decision_type: decision_type.to_string(),
            method: DecisionMethod::HardRule,
            title,
            rationale,
            input_evidence_ids,
            output_fact_ids,
            parameters_json,
            editable_status: EditableStatus::Draft,
            ..DesignDecisionNode::default()
        }
    }

    fn construct_reasoning_sequence_evidence_is_annotation_grade(
        evidence: &DesignEvidence,
    ) -> bool {
        if evidence.scope != EvidenceScope::SequenceSpan {
            return false;
        }
        if evidence
            .context_tags
            .iter()
            .any(|tag| tag == "cdna_confirmed")
            && matches!(
                evidence.role,
                ConstructRole::Exon
                    | ConstructRole::SpliceBoundary
                    | ConstructRole::Cds
                    | ConstructRole::Transcript
            )
        {
            return true;
        }
        let looks_generated = evidence.context_tags.iter().any(|tag| {
            matches!(
                tag.as_str(),
                "generated" | ANNOTATE_PROMOTER_WINDOWS_GENERATED_TAG
            )
        }) || evidence.provenance_kind.starts_with("derived_");
        looks_generated
            && !matches!(
                evidence.role,
                ConstructRole::RestrictionSite
                    | ConstructRole::Tfbs
                    | ConstructRole::ContextBaggage
            )
    }

    fn construct_reasoning_annotation_source_kind(
        evidence: &DesignEvidence,
        has_support: bool,
    ) -> String {
        if evidence.context_tags.iter().any(|tag| {
            matches!(
                tag.as_str(),
                "generated" | ANNOTATE_PROMOTER_WINDOWS_GENERATED_TAG
            )
        }) || evidence.provenance_kind.starts_with("derived_")
        {
            return "generated_annotation".to_string();
        }
        if evidence
            .context_tags
            .iter()
            .any(|tag| tag == "cdna_confirmed")
        {
            return "confirmed_annotation".to_string();
        }
        if has_support {
            return "supporting_annotation".to_string();
        }
        "imported_annotation".to_string()
    }

    fn construct_reasoning_variant_annotation_context(
        facts: &[DesignFact],
        evidence_id: &str,
    ) -> (Option<String>, Vec<String>) {
        for fact in facts {
            if fact.fact_type != "variant_effect_context"
                || !fact
                    .based_on_evidence_ids
                    .iter()
                    .any(|id| id == evidence_id)
            {
                continue;
            }
            let Some(variants) = fact
                .value_json
                .get("variants")
                .and_then(serde_json::Value::as_array)
            else {
                continue;
            };
            for variant in variants {
                let matches_evidence = variant
                    .get("evidence_id")
                    .and_then(serde_json::Value::as_str)
                    .map(|id| id == evidence_id)
                    .unwrap_or(false);
                if !matches_evidence {
                    continue;
                }
                let transcript_context_status = variant
                    .get("transcript_context_status")
                    .and_then(serde_json::Value::as_str)
                    .map(|value| value.to_string());
                let mut effect_tags = variant
                    .get("effect_tags")
                    .and_then(serde_json::Value::as_array)
                    .map(|rows| {
                        rows.iter()
                            .filter_map(serde_json::Value::as_str)
                            .map(|value| value.to_string())
                            .collect::<Vec<_>>()
                    })
                    .unwrap_or_default();
                effect_tags.sort();
                effect_tags.dedup();
                return (transcript_context_status, effect_tags);
            }
        }
        (None, vec![])
    }

    fn construct_reasoning_build_annotation_candidates(
        seq_id: &str,
        evidence: &[DesignEvidence],
        facts: &[DesignFact],
        decisions: &[DesignDecisionNode],
    ) -> Vec<AnnotationCandidate> {
        let mut fact_support_by_evidence: HashMap<&str, (Vec<String>, Vec<String>)> =
            HashMap::new();
        for fact in facts {
            for evidence_id in &fact.based_on_evidence_ids {
                let entry = fact_support_by_evidence
                    .entry(evidence_id.as_str())
                    .or_insert_with(|| (vec![], vec![]));
                entry.0.push(fact.fact_id.clone());
                if !fact.label.trim().is_empty() {
                    entry.1.push(fact.label.clone());
                }
            }
        }
        let mut decision_support_by_evidence: HashMap<&str, (Vec<String>, Vec<String>)> =
            HashMap::new();
        for decision in decisions {
            for evidence_id in &decision.input_evidence_ids {
                let entry = decision_support_by_evidence
                    .entry(evidence_id.as_str())
                    .or_insert_with(|| (vec![], vec![]));
                entry.0.push(decision.decision_id.clone());
                if !decision.title.trim().is_empty() {
                    entry.1.push(decision.title.clone());
                }
            }
        }

        let mut rows = vec![];
        for evidence_row in evidence {
            if evidence_row.scope != EvidenceScope::SequenceSpan {
                continue;
            }
            let (supporting_fact_ids, supporting_fact_labels) = fact_support_by_evidence
                .get(evidence_row.evidence_id.as_str())
                .cloned()
                .unwrap_or_default();
            let (supporting_decision_ids, supporting_decision_titles) =
                decision_support_by_evidence
                    .get(evidence_row.evidence_id.as_str())
                    .cloned()
                    .unwrap_or_default();
            let has_support =
                !supporting_fact_ids.is_empty() || !supporting_decision_ids.is_empty();
            if !has_support
                && !Self::construct_reasoning_sequence_evidence_is_annotation_grade(evidence_row)
            {
                continue;
            }
            let (transcript_context_status, effect_tags) =
                Self::construct_reasoning_variant_annotation_context(
                    facts,
                    &evidence_row.evidence_id,
                );
            rows.push(AnnotationCandidate {
                annotation_id: format!(
                    "annotation_{}_{}",
                    Self::normalize_id_token(evidence_row.role.as_str()),
                    Self::normalize_id_token(&evidence_row.evidence_id)
                ),
                evidence_id: evidence_row.evidence_id.clone(),
                seq_id: if evidence_row.seq_id.trim().is_empty() {
                    seq_id.to_string()
                } else {
                    evidence_row.seq_id.clone()
                },
                start_0based: evidence_row.start_0based,
                end_0based_exclusive: evidence_row.end_0based_exclusive,
                strand: evidence_row.strand.clone(),
                role: evidence_row.role,
                label: evidence_row.label.clone(),
                rationale: evidence_row.rationale.clone(),
                source_kind: Self::construct_reasoning_annotation_source_kind(
                    evidence_row,
                    has_support,
                ),
                supporting_fact_ids,
                supporting_fact_labels,
                supporting_decision_ids,
                supporting_decision_titles,
                transcript_context_status,
                effect_tags,
                editable_status: evidence_row.editable_status,
                warnings: evidence_row.warnings.clone(),
                notes: evidence_row.notes.clone(),
                ..AnnotationCandidate::default()
            });
        }
        rows
    }

    fn build_construct_reasoning_facts_and_decisions(
        seq_id: &str,
        dna: &DNAsequence,
        objective: &ConstructObjective,
        evidence: &[DesignEvidence],
        helper_interpretation: Option<&HelperConstructInterpretation>,
    ) -> (Vec<DesignFact>, Vec<DesignDecisionNode>) {
        let mut facts: Vec<DesignFact> = vec![];
        let mut decisions: Vec<DesignDecisionNode> = vec![];

        let propagation_evidence_ids =
            Self::construct_reasoning_evidence_ids_matching(evidence, |row| {
                row.scope == EvidenceScope::HostProfile
                    && row.context_tags.iter().any(|tag| tag == "propagation_host")
            });
        let expression_evidence_ids =
            Self::construct_reasoning_evidence_ids_matching(evidence, |row| {
                row.scope == EvidenceScope::HostProfile
                    && row.context_tags.iter().any(|tag| tag == "expression_host")
            });
        let host_route_evidence_ids =
            Self::construct_reasoning_evidence_ids_matching(evidence, |row| {
                row.scope == EvidenceScope::HostTransition
            });
        let helper_evidence_ids =
            Self::construct_reasoning_evidence_ids_matching(evidence, |row| {
                row.scope == EvidenceScope::HelperProfile
            });
        let adapter_capture_evidence_ids =
            Self::construct_reasoning_evidence_ids_matching(evidence, |row| {
                row.scope == EvidenceScope::WholeConstruct
                    && row.role == ConstructRole::ContextBaggage
                    && row
                        .context_tags
                        .iter()
                        .any(|tag| tag == "adapter_restriction_capture")
            });
        let medium_evidence_ids =
            Self::construct_reasoning_evidence_ids_matching(evidence, |row| {
                row.scope == EvidenceScope::MediumCondition
            });
        let variant_evidence_rows = evidence
            .iter()
            .filter(|row| {
                row.scope == EvidenceScope::SequenceSpan && row.role == ConstructRole::Variant
            })
            .collect::<Vec<_>>();
        let condition_signals = Self::construct_reasoning_interpret_medium_conditions(objective);
        let host_constraint_evidence_ids =
            Self::construct_reasoning_evidence_ids_matching(evidence, |row| {
                row.scope == EvidenceScope::WholeConstruct
                    && row.context_tags.iter().any(|tag| {
                        matches!(
                            tag.as_str(),
                            "required_host_trait"
                                | "forbidden_host_trait"
                                | "host_species"
                                | "cell_type"
                                | "tissue"
                                | "organelle"
                                | "expression_intent"
                        )
                    })
            });
        let (
            adapter_capture_plan_summaries,
            adapter_capture_review_items,
            adapter_capture_derived_preferred_routine_families,
        ) = Self::construct_reasoning_summarize_adapter_restriction_capture_plans(dna, objective);
        let adapter_capture_restriction_evidence_ids = if adapter_capture_plan_summaries.is_empty()
        {
            vec![]
        } else {
            let resolved_enzyme_names = adapter_capture_plan_summaries
                .iter()
                .filter(|row| row.enzyme_resolution_status == "resolved")
                .map(|row| row.restriction_enzyme_name.clone())
                .collect::<BTreeSet<_>>();
            Self::construct_reasoning_evidence_ids_matching(evidence, |row| {
                row.scope == EvidenceScope::SequenceSpan
                    && row.role == ConstructRole::RestrictionSite
                    && resolved_enzyme_names.contains(&row.label)
            })
        };
        let host_fit_context_present = objective.propagation_host_profile_id.is_some()
            || objective.expression_host_profile_id.is_some()
            || !objective.host_route.is_empty()
            || !objective.required_host_traits.is_empty()
            || !objective.forbidden_host_traits.is_empty()
            || objective.host_species.is_some()
            || objective.cell_type.is_some()
            || objective.tissue.is_some()
            || objective.organelle.is_some()
            || objective.expression_intent.is_some();

        if host_fit_context_present {
            let propagation_status = if objective.propagation_host_profile_id.is_some() {
                "specified"
            } else {
                "unspecified"
            };
            let propagation_label =
                if let Some(profile_id) = objective.propagation_host_profile_id.as_ref() {
                    format!("Propagation host specified: {profile_id}")
                } else {
                    "Propagation host unspecified".to_string()
                };
            let propagation_rationale = if let Some(profile_id) =
                objective.propagation_host_profile_id.as_ref()
            {
                format!(
                    "Construct objective explicitly selects propagation host profile '{profile_id}'. Current deterministic host-fit reasoning records the choice and the stated host constraints without yet consulting a curated strain catalog."
                )
            } else {
                "Construct objective does not yet select a dedicated propagation host profile. Current deterministic host-fit reasoning treats this as an inspectable gap rather than guessing a host.".to_string()
            };
            let propagation_fact = Self::construct_reasoning_build_fact(
                "fact_propagation_host_context",
                "propagation_host_context",
                propagation_label,
                propagation_rationale.clone(),
                propagation_evidence_ids.clone(),
                json!({
                    "seq_id": seq_id,
                    "status": propagation_status,
                    "profile_id": objective.propagation_host_profile_id.clone(),
                    "required_host_traits": objective.required_host_traits.clone(),
                    "forbidden_host_traits": objective.forbidden_host_traits.clone(),
                }),
            );
            facts.push(propagation_fact);
            let mut propagation_input_ids = propagation_evidence_ids.clone();
            propagation_input_ids.extend(host_constraint_evidence_ids.clone());
            propagation_input_ids.sort();
            propagation_input_ids.dedup();
            decisions.push(Self::construct_reasoning_build_decision(
                "decision_evaluate_propagation_host_fit",
                "evaluate_propagation_host_fit",
                "Evaluate Propagation Host Fit".to_string(),
                propagation_rationale,
                propagation_input_ids,
                vec!["fact_propagation_host_context".to_string()],
                json!({
                    "status": propagation_status,
                    "profile_id": objective.propagation_host_profile_id.clone(),
                }),
            ));

            let expression_status = if objective.expression_host_profile_id.is_some() {
                "specified"
            } else {
                "unspecified"
            };
            let expression_label =
                if let Some(profile_id) = objective.expression_host_profile_id.as_ref() {
                    format!("Expression host specified: {profile_id}")
                } else {
                    "Expression host unspecified".to_string()
                };
            let expression_rationale = if let Some(profile_id) =
                objective.expression_host_profile_id.as_ref()
            {
                format!(
                    "Construct objective explicitly selects expression host profile '{profile_id}'. Current deterministic host-fit reasoning records the downstream host together with the stated species/cell/tissue context without yet applying catalog-based suitability scoring."
                )
            } else {
                "Construct objective does not yet select a dedicated expression host profile. Current deterministic host-fit reasoning keeps that gap explicit so later host/cell-type evaluation can stay inspectable.".to_string()
            };
            let expression_fact = Self::construct_reasoning_build_fact(
                "fact_expression_host_context",
                "expression_host_context",
                expression_label,
                expression_rationale.clone(),
                expression_evidence_ids.clone(),
                json!({
                    "seq_id": seq_id,
                    "status": expression_status,
                    "profile_id": objective.expression_host_profile_id.clone(),
                    "host_species": objective.host_species.clone(),
                    "cell_type": objective.cell_type.clone(),
                    "tissue": objective.tissue.clone(),
                    "organelle": objective.organelle.clone(),
                    "expression_intent": objective.expression_intent.clone(),
                }),
            );
            facts.push(expression_fact);
            let mut expression_input_ids = expression_evidence_ids.clone();
            expression_input_ids.extend(host_constraint_evidence_ids.clone());
            expression_input_ids.sort();
            expression_input_ids.dedup();
            decisions.push(Self::construct_reasoning_build_decision(
                "decision_evaluate_expression_host_fit",
                "evaluate_expression_host_fit",
                "Evaluate Expression Host Fit".to_string(),
                expression_rationale,
                expression_input_ids,
                vec!["fact_expression_host_context".to_string()],
                json!({
                    "status": expression_status,
                    "profile_id": objective.expression_host_profile_id.clone(),
                }),
            ));
        }

        if objective.propagation_host_profile_id.is_some()
            || objective.expression_host_profile_id.is_some()
            || !objective.host_route.is_empty()
        {
            let mut ordered_hosts: Vec<String> = vec![];
            if let Some(profile_id) = objective.propagation_host_profile_id.as_ref() {
                ordered_hosts.push(profile_id.clone());
            }
            ordered_hosts.extend(
                objective
                    .host_route
                    .iter()
                    .map(|step| step.host_profile_id.clone()),
            );
            if let Some(profile_id) = objective.expression_host_profile_id.as_ref() {
                ordered_hosts.push(profile_id.clone());
            }
            let distinct_host_profiles = ordered_hosts.iter().cloned().collect::<BTreeSet<_>>();
            let transition_status = match (
                objective.propagation_host_profile_id.as_ref(),
                objective.expression_host_profile_id.as_ref(),
                objective.host_route.is_empty(),
            ) {
                (Some(left), Some(right), true) if left == right => "same_host",
                (Some(left), Some(right), false) if left == right => "explicit_route",
                (Some(left), Some(right), true) if left != right => "review_needed",
                (Some(_), Some(_), false) => "explicit_route",
                (None, None, false) | (Some(_), None, false) | (None, Some(_), false) => {
                    "partial_route_only"
                }
                (Some(_), None, true) | (None, Some(_), true) => "one_endpoint_only",
                _ => "unspecified",
            };
            let transition_label = match transition_status {
                "same_host" => "Host transition stays on one host".to_string(),
                "explicit_route" => "Host transition route recorded".to_string(),
                "review_needed" => "Host transition requires review".to_string(),
                "partial_route_only" => "Partial host route recorded".to_string(),
                "one_endpoint_only" => "Host transition only partially specified".to_string(),
                _ => "Host transition unspecified".to_string(),
            };
            let transition_rationale = match transition_status {
                "same_host" => {
                    "Propagation and expression currently point to the same host profile, so no inter-host transfer step is required by the recorded objective context.".to_string()
                }
                "explicit_route" => format!(
                    "Construct objective records a host route with {} explicit step(s), so inter-host transfer context is inspectable instead of implicit.",
                    objective.host_route.len()
                ),
                "review_needed" => "Propagation and expression point to different host profiles, but no explicit host route is recorded yet. Current deterministic reasoning flags that missing transition documentation for review rather than inventing one.".to_string(),
                "partial_route_only" => format!(
                    "Construct objective records {} host-route step(s) but does not fully anchor them with both propagation and expression endpoints.",
                    objective.host_route.len()
                ),
                "one_endpoint_only" => "Construct objective records only one host endpoint and no explicit intermediate route, so host transfer intent remains only partially documented.".to_string(),
                _ => "Construct objective does not yet provide host-transition context.".to_string(),
            };
            let mut transition_evidence_ids = propagation_evidence_ids.clone();
            transition_evidence_ids.extend(expression_evidence_ids.clone());
            transition_evidence_ids.extend(host_route_evidence_ids.clone());
            transition_evidence_ids.sort();
            transition_evidence_ids.dedup();
            facts.push(Self::construct_reasoning_build_fact(
                "fact_host_transition_context",
                "host_transition_context",
                transition_label,
                transition_rationale.clone(),
                transition_evidence_ids.clone(),
                json!({
                    "seq_id": seq_id,
                    "status": transition_status,
                    "propagation_host_profile_id": objective.propagation_host_profile_id.clone(),
                    "expression_host_profile_id": objective.expression_host_profile_id.clone(),
                    "host_route_step_count": objective.host_route.len(),
                    "ordered_host_profile_ids": ordered_hosts,
                    "distinct_host_profile_count": distinct_host_profiles.len(),
                }),
            ));
            decisions.push(Self::construct_reasoning_build_decision(
                "decision_evaluate_host_transition_risk",
                "evaluate_host_transition_risk",
                "Evaluate Host Transition Risk".to_string(),
                transition_rationale,
                transition_evidence_ids,
                vec!["fact_host_transition_context".to_string()],
                json!({
                    "status": transition_status,
                    "host_route_step_count": objective.host_route.len(),
                }),
            ));
        }

        let host_route_restriction_steps =
            Self::construct_reasoning_build_host_route_restriction_methylation_steps(objective);
        let sequence_restriction_methylation_patterns =
            Self::construct_reasoning_scan_sequence_restriction_methylation_patterns(dna);
        let route_restriction_methylation_conflicts =
            Self::construct_reasoning_detect_route_restriction_methylation_conflicts(
                &host_route_restriction_steps,
                &sequence_restriction_methylation_patterns,
            );
        if !host_route_restriction_steps.is_empty()
            || !sequence_restriction_methylation_patterns
                .present_pattern_tokens
                .is_empty()
        {
            let route_status = if !route_restriction_methylation_conflicts.is_empty() {
                "review_needed"
            } else if !host_route_restriction_steps.is_empty() {
                "context_recorded"
            } else {
                "sequence_patterns_only"
            };
            let route_label = match route_status {
                "review_needed" => "Restriction/methylation transfer route requires review",
                "context_recorded" => "Restriction/methylation transfer context recorded",
                _ => "Sequence restriction/methylation patterns recorded",
            }
            .to_string();
            let route_rationale = if !route_restriction_methylation_conflicts.is_empty() {
                route_restriction_methylation_conflicts
                    .iter()
                    .map(|conflict| conflict.rationale.clone())
                    .collect::<Vec<_>>()
                    .join(" ")
            } else if !host_route_restriction_steps.is_empty() {
                format!(
                    "Construct objective records host-route restriction/methylation traits, and the current deterministic engine also summarizes Dam/Dcm plus EcoK target-like sequence motifs (Dam={}, Dcm={}, EcoK-like={}).",
                    sequence_restriction_methylation_patterns.dam_site_count,
                    sequence_restriction_methylation_patterns.dcm_site_count,
                    sequence_restriction_methylation_patterns.eco_k_target_site_count
                )
            } else {
                format!(
                    "Sequence scan records methylation-sensitive motifs for downstream host-route review (Dam={}, Dcm={}, EcoK-like={}).",
                    sequence_restriction_methylation_patterns.dam_site_count,
                    sequence_restriction_methylation_patterns.dcm_site_count,
                    sequence_restriction_methylation_patterns.eco_k_target_site_count
                )
            };
            let mut route_evidence_ids = host_route_evidence_ids.clone();
            route_evidence_ids.extend(host_constraint_evidence_ids.clone());
            route_evidence_ids.sort();
            route_evidence_ids.dedup();
            facts.push(Self::construct_reasoning_build_fact(
                "fact_host_restriction_methylation_context",
                "host_restriction_methylation_context",
                route_label,
                route_rationale.clone(),
                route_evidence_ids.clone(),
                json!({
                    "seq_id": seq_id,
                    "status": route_status,
                    "route_step_signals": host_route_restriction_steps,
                    "sequence_pattern_counts": sequence_restriction_methylation_patterns,
                    "detected_conflicts": route_restriction_methylation_conflicts,
                }),
            ));
            decisions.push(Self::construct_reasoning_build_decision(
                "decision_evaluate_methylation_restriction_risk",
                "evaluate_methylation_restriction_risk",
                "Evaluate Restriction/Methylation Transfer Risk".to_string(),
                route_rationale,
                route_evidence_ids,
                vec!["fact_host_restriction_methylation_context".to_string()],
                json!({
                    "status": route_status,
                }),
            ));
        }

        if !adapter_capture_plan_summaries.is_empty() {
            let adapter_capture_status = if !adapter_capture_review_items.is_empty() {
                "review_needed"
            } else {
                "context_recorded"
            };
            let adapter_capture_label = if !adapter_capture_review_items.is_empty() {
                "Adapter/linker restriction capture requires review".to_string()
            } else {
                "Adapter/linker restriction capture context recorded".to_string()
            };
            let adapter_capture_rationale = if !adapter_capture_review_items.is_empty() {
                adapter_capture_review_items
                    .iter()
                    .map(|item| item.rationale.clone())
                    .collect::<Vec<_>>()
                    .join(" ")
            } else {
                format!(
                    "Construct objective records {} adapter/linker capture plan(s), and the current deterministic engine inspects internal reuse of the chosen capture restriction site(s) on the insert together with any requested insert-methylation protection mode.",
                    adapter_capture_plan_summaries.len()
                )
            };
            let mut adapter_fact_evidence_ids = adapter_capture_evidence_ids.clone();
            adapter_fact_evidence_ids.extend(adapter_capture_restriction_evidence_ids.clone());
            adapter_fact_evidence_ids.sort();
            adapter_fact_evidence_ids.dedup();
            facts.push(Self::construct_reasoning_build_fact(
                "fact_adapter_restriction_capture_context",
                "adapter_restriction_capture_context",
                adapter_capture_label,
                adapter_capture_rationale.clone(),
                adapter_fact_evidence_ids.clone(),
                json!({
                    "seq_id": seq_id,
                    "status": adapter_capture_status,
                    "capture_plans": adapter_capture_plan_summaries,
                    "review_items": adapter_capture_review_items,
                    "derived_preferred_routine_families": adapter_capture_derived_preferred_routine_families,
                }),
            ));
            decisions.push(Self::construct_reasoning_build_decision(
                "decision_evaluate_adapter_restriction_capture_strategy",
                "evaluate_adapter_restriction_capture_strategy",
                "Evaluate Adapter/Linker Restriction Capture".to_string(),
                adapter_capture_rationale,
                adapter_fact_evidence_ids,
                vec!["fact_adapter_restriction_capture_context".to_string()],
                json!({
                    "status": adapter_capture_status,
                }),
            ));
        }

        if !objective.medium_conditions.is_empty() {
            let signal_categories = condition_signals
                .iter()
                .map(|signal| signal.category.clone())
                .collect::<BTreeSet<_>>()
                .into_iter()
                .collect::<Vec<_>>();
            let temperature_celsius_values = condition_signals
                .iter()
                .filter_map(|signal| signal.temperature_celsius)
                .map(|value| format!("{value:.2}"))
                .collect::<BTreeSet<_>>()
                .into_iter()
                .filter_map(|value| value.parse::<f64>().ok())
                .collect::<Vec<_>>();
            let growth_label = if condition_signals.is_empty() {
                "Growth/condition context recorded".to_string()
            } else {
                "Growth/condition context interpreted".to_string()
            };
            let growth_rationale = if condition_signals.is_empty() {
                "Construct objective records medium/growth conditions, but the current deterministic interpreter did not yet classify any canonical condition signals from them.".to_string()
            } else {
                format!(
                    "Construct objective records medium/growth conditions that the deterministic condition interpreter normalizes into inspectable signal(s): {}.",
                    condition_signals
                        .iter()
                        .map(|signal| signal.label.clone())
                        .collect::<BTreeSet<_>>()
                        .into_iter()
                        .collect::<Vec<_>>()
                        .join(", ")
                )
            };
            facts.push(Self::construct_reasoning_build_fact(
                "fact_growth_condition_context",
                "growth_condition_context",
                growth_label,
                growth_rationale.clone(),
                medium_evidence_ids.clone(),
                json!({
                    "seq_id": seq_id,
                    "medium_conditions": objective.medium_conditions.clone(),
                    "condition_signals": condition_signals.clone(),
                    "signal_categories": signal_categories.clone(),
                    "temperature_celsius_values": temperature_celsius_values.clone(),
                }),
            ));
            decisions.push(Self::construct_reasoning_build_decision(
                "decision_evaluate_growth_condition_context",
                "evaluate_growth_condition_context",
                "Evaluate Growth/Condition Context".to_string(),
                growth_rationale,
                medium_evidence_ids.clone(),
                vec!["fact_growth_condition_context".to_string()],
                json!({
                    "signal_count": condition_signals.len(),
                    "signal_categories": signal_categories,
                }),
            ));
        }

        let helper_mcs_feature_labels = dna
            .features()
            .iter()
            .enumerate()
            .filter(|(_, feature)| Self::feature_looks_like_mcs(feature))
            .map(|(feature_id, feature)| Self::feature_display_label(feature, feature_id))
            .collect::<Vec<_>>();
        let helper_offered_functions = helper_interpretation
            .map(|row| row.offered_functions.clone())
            .unwrap_or_default();
        let helper_constraints = helper_interpretation
            .map(|row| row.constraints.clone())
            .unwrap_or_default();
        let helper_component_labels = helper_interpretation
            .map(|row| {
                row.components
                    .iter()
                    .map(|component| {
                        component
                            .label
                            .clone()
                            .unwrap_or_else(|| component.id.clone())
                    })
                    .collect::<Vec<_>>()
            })
            .unwrap_or_default();
        let variant_summaries = Self::construct_reasoning_collect_variant_summaries(
            dna,
            evidence,
            &variant_evidence_rows,
        );
        let variant_effect_tags = variant_summaries
            .iter()
            .flat_map(|row| row.effect_tags.iter().cloned())
            .collect::<BTreeSet<_>>()
            .into_iter()
            .collect::<Vec<_>>();
        let suggested_variant_assay_ids = variant_summaries
            .iter()
            .flat_map(|row| row.suggested_assay_ids.iter().cloned())
            .collect::<BTreeSet<_>>()
            .into_iter()
            .collect::<Vec<_>>();
        let variant_derived_preferred_routine_families =
            Self::variant_assay_ids_to_routine_family_preferences(&suggested_variant_assay_ids);
        let routine_preference_context = Self::construct_objective_routine_preference_context(
            objective,
            helper_interpretation,
            (!variant_summaries.is_empty()).then_some(seq_id),
            &adapter_capture_derived_preferred_routine_families,
            &variant_effect_tags,
            &suggested_variant_assay_ids,
            &variant_derived_preferred_routine_families,
        );
        if objective.helper_profile_id.is_some() || !helper_mcs_feature_labels.is_empty() {
            let helper_status = match (
                objective.helper_profile_id.as_ref(),
                helper_mcs_feature_labels.is_empty(),
            ) {
                (Some(_), false) => "helper_profile_with_mcs_hint",
                (Some(_), true) => "helper_profile_selected",
                (None, false) => "mcs_feature_detected",
                (None, true) => "unspecified",
            };
            let helper_label = match helper_status {
                "helper_profile_with_mcs_hint" => {
                    "Helper context documented with MCS hint".to_string()
                }
                "helper_profile_selected" => "Helper context documented".to_string(),
                "mcs_feature_detected" => "Helper-like sequence context detected".to_string(),
                _ => "Helper context unspecified".to_string(),
            };
            let helper_rationale = match (
                objective.helper_profile_id.as_ref(),
                helper_mcs_feature_labels.is_empty(),
            ) {
                (Some(profile_id), false) => format!(
                    "Construct objective selects helper/vector profile '{profile_id}', and sequence feature inspection also sees MCS-style annotation(s): {}.",
                    helper_mcs_feature_labels.join(", ")
                ),
                (Some(profile_id), true) => format!(
                    "Construct objective selects helper/vector profile '{profile_id}'. Current deterministic helper reasoning records that explicit choice even when no MCS-style sequence annotation is visible yet."
                ),
                (None, false) => format!(
                    "Sequence feature inspection sees MCS-style annotation(s): {}. Current deterministic reasoning records this as helper-like context even without an explicit helper profile selection.",
                    helper_mcs_feature_labels.join(", ")
                ),
                (None, true) => "Helper context is not yet specified.".to_string(),
            };
            facts.push(Self::construct_reasoning_build_fact(
                "fact_helper_context",
                "helper_context",
                helper_label,
                helper_rationale.clone(),
                helper_evidence_ids.clone(),
                json!({
                    "seq_id": seq_id,
                    "status": helper_status,
                    "helper_profile_id": objective.helper_profile_id.clone(),
                    "mcs_feature_labels": helper_mcs_feature_labels,
                    "helper_offered_functions": helper_offered_functions,
                    "helper_constraints": helper_constraints,
                    "helper_component_labels": helper_component_labels,
                }),
            ));
            decisions.push(Self::construct_reasoning_build_decision(
                "decision_evaluate_helper_profile_context",
                "evaluate_helper_profile_context",
                "Evaluate Helper Profile Context".to_string(),
                helper_rationale,
                helper_evidence_ids.clone(),
                vec!["fact_helper_context".to_string()],
                json!({
                    "status": helper_status,
                }),
            ));
        }

        let repeat_region_rows = evidence
            .iter()
            .filter(|row| {
                row.scope == EvidenceScope::SequenceSpan && row.role == ConstructRole::RepeatRegion
            })
            .collect::<Vec<_>>();
        let mobile_element_rows = evidence
            .iter()
            .filter(|row| {
                row.scope == EvidenceScope::SequenceSpan && row.role == ConstructRole::MobileElement
            })
            .collect::<Vec<_>>();
        let has_similarity_tag =
            |row: &DesignEvidence, tag: &str| row.context_tags.iter().any(|value| value == tag);
        let similarity_labels = |rows: &[&DesignEvidence]| {
            rows.iter()
                .filter_map(|row| {
                    let value = row.label.trim();
                    (!value.is_empty()).then(|| value.to_string())
                })
                .collect::<BTreeSet<_>>()
                .into_iter()
                .collect::<Vec<_>>()
        };
        let similarity_ids = |rows: &[&DesignEvidence]| {
            rows.iter()
                .map(|row| row.evidence_id.clone())
                .collect::<BTreeSet<_>>()
                .into_iter()
                .collect::<Vec<_>>()
        };
        let similarity_max_score =
            |rows: &[&DesignEvidence]| rows.iter().filter_map(|row| row.score).fold(0.0, f64::max);

        let low_complexity_rows = repeat_region_rows
            .iter()
            .copied()
            .filter(|row| has_similarity_tag(row, "low_complexity"))
            .collect::<Vec<_>>();
        let homopolymer_rows = repeat_region_rows
            .iter()
            .copied()
            .filter(|row| has_similarity_tag(row, "homopolymer_run"))
            .collect::<Vec<_>>();
        let tandem_rows = repeat_region_rows
            .iter()
            .copied()
            .filter(|row| has_similarity_tag(row, "tandem_repeat"))
            .collect::<Vec<_>>();
        let direct_repeat_rows = repeat_region_rows
            .iter()
            .copied()
            .filter(|row| has_similarity_tag(row, "direct_repeat"))
            .collect::<Vec<_>>();
        let inverted_repeat_rows = repeat_region_rows
            .iter()
            .copied()
            .filter(|row| has_similarity_tag(row, "inverted_repeat"))
            .collect::<Vec<_>>();
        let slippage_risk_rows = repeat_region_rows
            .iter()
            .copied()
            .filter(|row| has_similarity_tag(row, "polymerase_slippage_risk"))
            .collect::<Vec<_>>();
        let nanopore_signal_rows = repeat_region_rows
            .iter()
            .copied()
            .filter(|row| has_similarity_tag(row, "nanopore_signal_risk"))
            .collect::<Vec<_>>();
        let nanopore_homopolymer_rows = repeat_region_rows
            .iter()
            .copied()
            .filter(|row| has_similarity_tag(row, "nanopore_homopolymer_risk"))
            .collect::<Vec<_>>();
        let mapping_ambiguity_rows = repeat_region_rows
            .iter()
            .copied()
            .filter(|row| has_similarity_tag(row, "mapping_ambiguity_risk"))
            .collect::<Vec<_>>();
        let inversion_risk_rows = repeat_region_rows
            .iter()
            .copied()
            .filter(|row| {
                has_similarity_tag(row, "inversion_risk")
                    || has_similarity_tag(row, "palindromic_repeat_risk")
            })
            .collect::<Vec<_>>();
        let cloning_risk_rows = repeat_region_rows
            .iter()
            .copied()
            .filter(|row| has_similarity_tag(row, "cloning_stability_risk"))
            .collect::<Vec<_>>();
        let alu_like_rows = mobile_element_rows
            .iter()
            .copied()
            .filter(|row| has_similarity_tag(row, "alu_like"))
            .collect::<Vec<_>>();
        let mobile_mapping_rows = mobile_element_rows
            .iter()
            .copied()
            .filter(|row| has_similarity_tag(row, "mapping_ambiguity_risk"))
            .collect::<Vec<_>>();

        if !low_complexity_rows.is_empty()
            || !homopolymer_rows.is_empty()
            || !tandem_rows.is_empty()
        {
            let complexity_status = if !homopolymer_rows.is_empty() || !tandem_rows.is_empty() {
                "repeat_patterns_detected"
            } else {
                "low_complexity_detected"
            };
            let complexity_label = match complexity_status {
                "repeat_patterns_detected" => {
                    "Low-complexity / tandem-repeat context detected".to_string()
                }
                _ => "Low-complexity context detected".to_string(),
            };
            let mut rationale_bits = vec![];
            if !low_complexity_rows.is_empty() {
                rationale_bits.push(format!(
                    "{} low-complexity region(s)",
                    low_complexity_rows.len()
                ));
            }
            if !homopolymer_rows.is_empty() {
                rationale_bits.push(format!("{} homopolymer run(s)", homopolymer_rows.len()));
            }
            if !tandem_rows.is_empty() {
                rationale_bits.push(format!("{} tandem repeat(s)", tandem_rows.len()));
            }
            let complexity_rationale = format!(
                "Sequence-derived complexity heuristics detected {}. These are kept as generated structural annotations rather than hard biological facts because their practical meaning depends on the downstream task.",
                rationale_bits.join(", ")
            );
            let mut complexity_evidence_ids = similarity_ids(&low_complexity_rows);
            complexity_evidence_ids.extend(similarity_ids(&homopolymer_rows));
            complexity_evidence_ids.extend(similarity_ids(&tandem_rows));
            complexity_evidence_ids.sort();
            complexity_evidence_ids.dedup();
            facts.push(Self::construct_reasoning_build_fact(
                "fact_sequence_complexity_context",
                "sequence_complexity_context",
                complexity_label,
                complexity_rationale.clone(),
                complexity_evidence_ids.clone(),
                json!({
                    "seq_id": seq_id,
                    "status": complexity_status,
                    "low_complexity_region_count": low_complexity_rows.len(),
                    "homopolymer_run_count": homopolymer_rows.len(),
                    "tandem_repeat_count": tandem_rows.len(),
                    "labels": similarity_labels(&low_complexity_rows),
                    "homopolymer_labels": similarity_labels(&homopolymer_rows),
                    "tandem_repeat_labels": similarity_labels(&tandem_rows),
                    "max_low_complexity_risk_score": similarity_max_score(&low_complexity_rows),
                    "max_polymerase_slippage_risk_score": similarity_max_score(&slippage_risk_rows),
                }),
            ));
            decisions.push(Self::construct_reasoning_build_decision(
                "decision_evaluate_sequence_complexity_context",
                "evaluate_sequence_complexity_context",
                "Evaluate Sequence Complexity".to_string(),
                complexity_rationale,
                complexity_evidence_ids,
                vec!["fact_sequence_complexity_context".to_string()],
                json!({
                    "status": complexity_status,
                }),
            ));
        }

        if !direct_repeat_rows.is_empty() || !inverted_repeat_rows.is_empty() {
            let repeat_architecture_status = if !inverted_repeat_rows.is_empty() {
                "direct_and_inverted_repeat_clusters"
            } else {
                "direct_repeat_clusters"
            };
            let repeat_architecture_label = if !inverted_repeat_rows.is_empty() {
                "Repeat/similarity architecture detected".to_string()
            } else {
                "Direct repeat architecture detected".to_string()
            };
            let repeat_architecture_rationale = format!(
                "Sequence-derived repeat scanning detected {} direct-repeat cluster(s) and {} inverted-repeat cluster(s). These exact sequence patterns are kept as generated structural context because inversion or recombination relevance depends on the workflow.",
                direct_repeat_rows.len(),
                inverted_repeat_rows.len()
            );
            let mut repeat_architecture_ids = similarity_ids(&direct_repeat_rows);
            repeat_architecture_ids.extend(similarity_ids(&inverted_repeat_rows));
            repeat_architecture_ids.sort();
            repeat_architecture_ids.dedup();
            facts.push(Self::construct_reasoning_build_fact(
                "fact_repeat_architecture_context",
                "repeat_architecture_context",
                repeat_architecture_label,
                repeat_architecture_rationale.clone(),
                repeat_architecture_ids.clone(),
                json!({
                    "seq_id": seq_id,
                    "status": repeat_architecture_status,
                    "direct_repeat_cluster_count": direct_repeat_rows.len(),
                    "inverted_repeat_cluster_count": inverted_repeat_rows.len(),
                    "direct_repeat_labels": similarity_labels(&direct_repeat_rows),
                    "inverted_repeat_labels": similarity_labels(&inverted_repeat_rows),
                    "max_direct_repeat_risk_score": similarity_max_score(&direct_repeat_rows),
                    "max_inverted_repeat_risk_score": similarity_max_score(&inverted_repeat_rows),
                }),
            ));
            decisions.push(Self::construct_reasoning_build_decision(
                "decision_evaluate_repeat_architecture_context",
                "evaluate_repeat_architecture_context",
                "Evaluate Repeat/Similarity Architecture".to_string(),
                repeat_architecture_rationale,
                repeat_architecture_ids,
                vec!["fact_repeat_architecture_context".to_string()],
                json!({
                    "status": repeat_architecture_status,
                }),
            ));
        }

        if !mobile_element_rows.is_empty() {
            let mobile_element_status = if !alu_like_rows.is_empty() {
                "alu_like_candidates_detected"
            } else {
                "mobile_element_candidates_detected"
            };
            let mobile_element_label = match mobile_element_status {
                "alu_like_candidates_detected" => {
                    "Alu-like mobile-element context detected".to_string()
                }
                _ => "Mobile-element context detected".to_string(),
            };
            let mobile_element_rationale = if !alu_like_rows.is_empty() {
                format!(
                    "Sequence-derived mobile-element heuristics detected {} Alu-like SINE candidate(s). These remain soft hypotheses until a curated repeat-family catalog or external masker confirms the family assignment.",
                    alu_like_rows.len()
                )
            } else {
                format!(
                    "Sequence-derived mobile-element heuristics detected {} candidate region(s), but none currently match the narrower Alu-like rule.",
                    mobile_element_rows.len()
                )
            };
            let mobile_element_ids = similarity_ids(&mobile_element_rows);
            facts.push(Self::construct_reasoning_build_fact(
                "fact_mobile_element_context",
                "mobile_element_context",
                mobile_element_label,
                mobile_element_rationale.clone(),
                mobile_element_ids.clone(),
                json!({
                    "seq_id": seq_id,
                    "status": mobile_element_status,
                    "mobile_element_candidate_count": mobile_element_rows.len(),
                    "alu_like_candidate_count": alu_like_rows.len(),
                    "labels": similarity_labels(&mobile_element_rows),
                    "max_mobile_element_score": similarity_max_score(&mobile_element_rows),
                }),
            ));
            decisions.push(Self::construct_reasoning_build_decision(
                "decision_evaluate_mobile_element_context",
                "evaluate_mobile_element_context",
                "Evaluate Mobile-Element Context".to_string(),
                mobile_element_rationale,
                mobile_element_ids,
                vec!["fact_mobile_element_context".to_string()],
                json!({
                    "status": mobile_element_status,
                }),
            ));
        }

        if !slippage_risk_rows.is_empty()
            || !homopolymer_rows.is_empty()
            || !tandem_rows.is_empty()
            || !low_complexity_rows.is_empty()
        {
            let pcr_risk_status = if !homopolymer_rows.is_empty() || !tandem_rows.is_empty() {
                "review_needed_repeat_slippage_risk"
            } else {
                "review_needed_low_complexity_pcr_risk"
            };
            let pcr_risk_rationale = format!(
                "Sequence-derived repeat/complexity annotations imply {} low-complexity span(s), {} homopolymer run(s), {} tandem-repeat span(s), and {} explicit slippage-risk span(s). GENtle keeps these as PCR/amplification review cues because the practical severity depends on primer placement, amplicon length, and polymerase choice.",
                low_complexity_rows.len(),
                homopolymer_rows.len(),
                tandem_rows.len(),
                slippage_risk_rows.len()
            );
            let mut pcr_risk_ids = similarity_ids(&low_complexity_rows);
            pcr_risk_ids.extend(similarity_ids(&homopolymer_rows));
            pcr_risk_ids.extend(similarity_ids(&tandem_rows));
            pcr_risk_ids.extend(similarity_ids(&slippage_risk_rows));
            pcr_risk_ids.sort();
            pcr_risk_ids.dedup();
            facts.push(Self::construct_reasoning_build_fact(
                "fact_pcr_operational_risk_context",
                "pcr_operational_risk_context",
                "PCR/amplification review suggested".to_string(),
                pcr_risk_rationale.clone(),
                pcr_risk_ids.clone(),
                json!({
                    "seq_id": seq_id,
                    "status": pcr_risk_status,
                    "low_complexity_region_count": low_complexity_rows.len(),
                    "homopolymer_run_count": homopolymer_rows.len(),
                    "tandem_repeat_count": tandem_rows.len(),
                    "polymerase_slippage_risk_count": slippage_risk_rows.len(),
                    "low_complexity_labels": similarity_labels(&low_complexity_rows),
                    "homopolymer_labels": similarity_labels(&homopolymer_rows),
                    "tandem_repeat_labels": similarity_labels(&tandem_rows),
                    "slippage_labels": similarity_labels(&slippage_risk_rows),
                    "max_pcr_risk_score": similarity_max_score(&low_complexity_rows)
                        .max(similarity_max_score(&slippage_risk_rows)),
                }),
            ));
            decisions.push(Self::construct_reasoning_build_decision(
                "decision_evaluate_pcr_operational_risk",
                "evaluate_pcr_operational_risk",
                "Evaluate PCR/Amplification Risk".to_string(),
                pcr_risk_rationale,
                pcr_risk_ids,
                vec!["fact_pcr_operational_risk_context".to_string()],
                json!({
                    "status": pcr_risk_status,
                }),
            ));
        }

        if !nanopore_signal_rows.is_empty() {
            let nanopore_risk_status = if !nanopore_homopolymer_rows.is_empty() {
                "review_needed_homopolymer_signal_risk"
            } else {
                "review_needed_low_complexity_signal_risk"
            };
            let nanopore_risk_rationale = format!(
                "Sequence-derived repeat/complexity annotations imply {} nanopore signal-risk span(s), including {} homopolymer-dominated span(s). GENtle keeps these as direct-sequencing review cues because homopolymers and low-complexity sequence can affect raw-signal interpretation, basecalling, and local mapping confidence.",
                nanopore_signal_rows.len(),
                nanopore_homopolymer_rows.len()
            );
            let mut nanopore_risk_ids = similarity_ids(&nanopore_signal_rows);
            nanopore_risk_ids.extend(similarity_ids(&nanopore_homopolymer_rows));
            nanopore_risk_ids.sort();
            nanopore_risk_ids.dedup();
            facts.push(Self::construct_reasoning_build_fact(
                "fact_nanopore_operational_risk_context",
                "nanopore_operational_risk_context",
                "Nanopore/direct-sequencing review suggested".to_string(),
                nanopore_risk_rationale.clone(),
                nanopore_risk_ids.clone(),
                json!({
                    "seq_id": seq_id,
                    "status": nanopore_risk_status,
                    "nanopore_signal_risk_count": nanopore_signal_rows.len(),
                    "nanopore_homopolymer_risk_count": nanopore_homopolymer_rows.len(),
                    "signal_labels": similarity_labels(&nanopore_signal_rows),
                    "homopolymer_labels": similarity_labels(&nanopore_homopolymer_rows),
                    "max_nanopore_risk_score": similarity_max_score(&nanopore_signal_rows),
                }),
            ));
            decisions.push(Self::construct_reasoning_build_decision(
                "decision_evaluate_nanopore_operational_risk",
                "evaluate_nanopore_operational_risk",
                "Evaluate Nanopore/Direct-Sequencing Risk".to_string(),
                nanopore_risk_rationale,
                nanopore_risk_ids,
                vec!["fact_nanopore_operational_risk_context".to_string()],
                json!({
                    "status": nanopore_risk_status,
                }),
            ));
        }

        if !mapping_ambiguity_rows.is_empty() || !mobile_mapping_rows.is_empty() {
            let mapping_risk_status = if !mobile_mapping_rows.is_empty() {
                "review_needed_repeat_family_mapping_risk"
            } else {
                "review_needed_repeat_mapping_risk"
            };
            let mapping_risk_rationale = format!(
                "Sequence-derived repeat/mobile-element annotations imply {} repeat-driven mapping-ambiguity span(s) and {} mobile-element candidate span(s). GENtle keeps these as mapping-review cues because ambiguity depends on read length, mapper behavior, and whether repeated family members are expected elsewhere in the genome.",
                mapping_ambiguity_rows.len(),
                mobile_mapping_rows.len()
            );
            let mut mapping_risk_ids = similarity_ids(&mapping_ambiguity_rows);
            mapping_risk_ids.extend(similarity_ids(&mobile_mapping_rows));
            mapping_risk_ids.sort();
            mapping_risk_ids.dedup();
            facts.push(Self::construct_reasoning_build_fact(
                "fact_mapping_operational_risk_context",
                "mapping_operational_risk_context",
                "Repeat-driven mapping review suggested".to_string(),
                mapping_risk_rationale.clone(),
                mapping_risk_ids.clone(),
                json!({
                    "seq_id": seq_id,
                    "status": mapping_risk_status,
                    "mapping_ambiguity_risk_count": mapping_ambiguity_rows.len(),
                    "mobile_element_candidate_count": mobile_mapping_rows.len(),
                    "mapping_labels": similarity_labels(&mapping_ambiguity_rows),
                    "mobile_element_labels": similarity_labels(&mobile_mapping_rows),
                    "max_mapping_risk_score": similarity_max_score(&mapping_ambiguity_rows)
                        .max(similarity_max_score(&mobile_mapping_rows)),
                }),
            ));
            decisions.push(Self::construct_reasoning_build_decision(
                "decision_evaluate_mapping_operational_risk",
                "evaluate_mapping_operational_risk",
                "Evaluate Repeat-Driven Mapping Risk".to_string(),
                mapping_risk_rationale,
                mapping_risk_ids,
                vec!["fact_mapping_operational_risk_context".to_string()],
                json!({
                    "status": mapping_risk_status,
                }),
            ));
        }

        if !inversion_risk_rows.is_empty() || !cloning_risk_rows.is_empty() {
            let cloning_risk_status = if !inversion_risk_rows.is_empty() {
                "review_needed_inversion_recombination_risk"
            } else {
                "review_needed_repeat_cloning_stability"
            };
            let cloning_risk_rationale = format!(
                "Sequence-derived repeat annotations imply {} inversion/palindromic-risk span(s) and {} cloning-stability span(s). GENtle keeps these as cloning review cues because repeated or inverted sequence can destabilize recovery, maintenance, or rearrangement-prone constructs.",
                inversion_risk_rows.len(),
                cloning_risk_rows.len()
            );
            let mut cloning_risk_ids = similarity_ids(&inversion_risk_rows);
            cloning_risk_ids.extend(similarity_ids(&cloning_risk_rows));
            cloning_risk_ids.sort();
            cloning_risk_ids.dedup();
            facts.push(Self::construct_reasoning_build_fact(
                "fact_cloning_stability_context",
                "cloning_stability_context",
                "Cloning stability review suggested".to_string(),
                cloning_risk_rationale.clone(),
                cloning_risk_ids.clone(),
                json!({
                    "seq_id": seq_id,
                    "status": cloning_risk_status,
                    "inversion_risk_count": inversion_risk_rows.len(),
                    "cloning_stability_risk_count": cloning_risk_rows.len(),
                    "inversion_labels": similarity_labels(&inversion_risk_rows),
                    "cloning_labels": similarity_labels(&cloning_risk_rows),
                    "max_cloning_risk_score": similarity_max_score(&inversion_risk_rows)
                        .max(similarity_max_score(&cloning_risk_rows)),
                }),
            ));
            decisions.push(Self::construct_reasoning_build_decision(
                "decision_evaluate_cloning_stability_context",
                "evaluate_cloning_stability_context",
                "Evaluate Cloning Stability".to_string(),
                cloning_risk_rationale,
                cloning_risk_ids,
                vec!["fact_cloning_stability_context".to_string()],
                json!({
                    "status": cloning_risk_status,
                }),
            ));
        }

        if !slippage_risk_rows.is_empty()
            || !mapping_ambiguity_rows.is_empty()
            || !inversion_risk_rows.is_empty()
            || !cloning_risk_rows.is_empty()
        {
            let similarity_risk_status = if !inversion_risk_rows.is_empty() {
                "review_needed_inversion_and_repeat_risk"
            } else if !cloning_risk_rows.is_empty() {
                "review_needed_cloning_repeat_risk"
            } else {
                "review_needed_repeat_risk"
            };
            let similarity_risk_rationale = format!(
                "Sequence-derived repeat/complexity annotations imply {} slippage-risk span(s), {} mapping-ambiguity span(s), {} inversion-risk span(s), and {} cloning-stability span(s). GENtle keeps these as operational review cues because the practical severity depends on PCR, mapping, nanopore, or cloning intent.",
                slippage_risk_rows.len(),
                mapping_ambiguity_rows.len(),
                inversion_risk_rows.len(),
                cloning_risk_rows.len()
            );
            let mut similarity_risk_ids = similarity_ids(&slippage_risk_rows);
            similarity_risk_ids.extend(similarity_ids(&mapping_ambiguity_rows));
            similarity_risk_ids.extend(similarity_ids(&inversion_risk_rows));
            similarity_risk_ids.extend(similarity_ids(&cloning_risk_rows));
            similarity_risk_ids.sort();
            similarity_risk_ids.dedup();
            facts.push(Self::construct_reasoning_build_fact(
                "fact_similarity_operational_risk_context",
                "similarity_operational_risk_context",
                "Similarity/low-complexity operational risks detected".to_string(),
                similarity_risk_rationale.clone(),
                similarity_risk_ids.clone(),
                json!({
                    "seq_id": seq_id,
                    "status": similarity_risk_status,
                    "polymerase_slippage_risk_count": slippage_risk_rows.len(),
                    "mapping_ambiguity_risk_count": mapping_ambiguity_rows.len(),
                    "inversion_risk_count": inversion_risk_rows.len(),
                    "cloning_stability_risk_count": cloning_risk_rows.len(),
                    "max_risk_score": similarity_max_score(&slippage_risk_rows)
                        .max(similarity_max_score(&mapping_ambiguity_rows))
                        .max(similarity_max_score(&inversion_risk_rows))
                        .max(similarity_max_score(&cloning_risk_rows)),
                }),
            ));
            decisions.push(Self::construct_reasoning_build_decision(
                "decision_evaluate_similarity_operational_risk",
                "evaluate_similarity_operational_risk",
                "Evaluate Similarity/Repeat Operational Risk".to_string(),
                similarity_risk_rationale,
                similarity_risk_ids,
                vec!["fact_similarity_operational_risk_context".to_string()],
                json!({
                    "status": similarity_risk_status,
                }),
            ));
        }

        if objective.helper_profile_id.is_some()
            || !objective.preferred_routine_families.is_empty()
            || !routine_preference_context
                .helper_derived_preferred_routine_families
                .is_empty()
            || !routine_preference_context
                .construct_strategy_derived_preferred_routine_families
                .is_empty()
            || !routine_preference_context
                .variant_derived_preferred_routine_families
                .is_empty()
        {
            let planning_status = match (
                routine_preference_context
                    .explicit_preferred_routine_families
                    .is_empty(),
                routine_preference_context
                    .helper_derived_preferred_routine_families
                    .is_empty(),
                routine_preference_context
                    .construct_strategy_derived_preferred_routine_families
                    .is_empty(),
                routine_preference_context
                    .variant_derived_preferred_routine_families
                    .is_empty(),
                routine_preference_context.helper_resolution_status.as_str(),
            ) {
                (false, false, false, false, _) => "explicit_helper_strategy_and_variant_derived",
                (false, false, false, true, _) => "explicit_helper_and_strategy_derived",
                (false, false, true, false, _) => "explicit_helper_and_variant_derived",
                (false, true, false, false, _) => "explicit_strategy_and_variant_derived",
                (true, false, false, false, _) => "helper_strategy_and_variant_derived",
                (false, false, true, true, _) => "explicit_and_helper_derived",
                (false, true, false, true, _) => "explicit_and_strategy_derived",
                (false, true, true, false, _) => "explicit_and_variant_derived",
                (true, false, false, true, _) => "helper_and_strategy_derived",
                (true, false, true, false, _) => "helper_and_variant_derived",
                (true, true, false, false, _) => "strategy_and_variant_derived",
                (false, true, true, true, _) => "explicit_only",
                (true, false, true, true, _) => "helper_derived",
                (true, true, false, true, _) => "strategy_derived",
                (true, true, true, false, _) => "variant_derived",
                (true, true, true, true, "error") => "helper_resolution_error",
                (true, true, true, true, "not_found") => "helper_not_found",
                _ => "unspecified",
            };
            let planning_label = match planning_status {
                "explicit_helper_strategy_and_variant_derived" => {
                    "Routine planning preferences synthesized".to_string()
                }
                "explicit_helper_and_strategy_derived" => {
                    "Routine planning preferences synthesized".to_string()
                }
                "explicit_helper_and_variant_derived" => {
                    "Routine planning preferences synthesized".to_string()
                }
                "explicit_strategy_and_variant_derived" => {
                    "Routine planning preferences synthesized".to_string()
                }
                "helper_strategy_and_variant_derived" => {
                    "Routine planning preferences synthesized".to_string()
                }
                "explicit_and_helper_derived" => {
                    "Routine planning preferences synthesized".to_string()
                }
                "explicit_and_strategy_derived" => {
                    "Routine planning preferences synthesized".to_string()
                }
                "explicit_and_variant_derived" => {
                    "Routine planning preferences synthesized".to_string()
                }
                "helper_and_strategy_derived" => {
                    "Routine planning preferences synthesized".to_string()
                }
                "helper_and_variant_derived" => {
                    "Routine planning preferences synthesized".to_string()
                }
                "strategy_and_variant_derived" => {
                    "Routine planning preferences synthesized".to_string()
                }
                "explicit_only" => "Routine planning preferences recorded".to_string(),
                "helper_derived" => {
                    "Routine planning preferences derived from helper profile".to_string()
                }
                "strategy_derived" => {
                    "Routine planning preferences derived from construct strategy".to_string()
                }
                "variant_derived" => {
                    "Routine planning preferences derived from variant assay context".to_string()
                }
                "helper_resolution_error" => {
                    "Helper-derived routine planning context could not be resolved".to_string()
                }
                "helper_not_found" => {
                    "Helper-derived routine planning context missing from catalog".to_string()
                }
                _ => "Routine planning context unspecified".to_string(),
            };
            let planning_rationale = if !routine_preference_context.rationale.is_empty() {
                routine_preference_context.rationale.join(" ")
            } else {
                "Current deterministic construct reasoning did not synthesize any routine-family preferences from the recorded objective context."
                    .to_string()
            };
            facts.push(Self::construct_reasoning_build_fact(
                "fact_routine_planning_context",
                "routine_planning_context",
                planning_label,
                planning_rationale.clone(),
                helper_evidence_ids.clone(),
                json!({
                    "seq_id": seq_id,
                    "status": planning_status,
                    "context": routine_preference_context,
                }),
            ));
            decisions.push(Self::construct_reasoning_build_decision(
                "decision_prioritize_helper_compatible_routines",
                "prioritize_helper_compatible_routines",
                "Prioritize Helper-Compatible Routines".to_string(),
                planning_rationale,
                helper_evidence_ids.clone(),
                vec!["fact_routine_planning_context".to_string()],
                json!({
                    "status": planning_status,
                    "effective_preferred_routine_families": routine_preference_context
                        .effective_preferred_routine_families,
                }),
            ));
        }

        let selection_rule_matches = Self::construct_reasoning_selection_rules()
            .iter()
            .filter_map(|rule| {
                Self::construct_reasoning_build_selection_rule_match(
                    rule,
                    &condition_signals,
                    evidence,
                    helper_interpretation,
                )
            })
            .collect::<Vec<_>>();
        let selection_candidates = selection_rule_matches
            .iter()
            .flat_map(|row| row.candidate_labels.iter().cloned())
            .collect::<BTreeSet<_>>();
        let selection_candidate_ids = selection_rule_matches
            .iter()
            .flat_map(|row| row.candidate_evidence_ids.iter().cloned())
            .collect::<BTreeSet<_>>()
            .into_iter()
            .collect::<Vec<_>>();
        let matching_medium_conditions = selection_rule_matches
            .iter()
            .flat_map(|row| row.matched_medium_conditions.iter().cloned())
            .collect::<BTreeSet<_>>()
            .into_iter()
            .collect::<Vec<_>>();
        let matching_rule_ids = selection_rule_matches
            .iter()
            .map(|row| row.rule_id.clone())
            .collect::<Vec<_>>();
        let helper_backed_selection_present = selection_rule_matches.iter().any(|row| {
            row.helper_profile_id.is_some()
                || !row.helper_offered_function_matches.is_empty()
                || !row.helper_component_ids.is_empty()
        });
        if !objective.medium_conditions.is_empty()
            || !selection_candidates.is_empty()
            || !selection_rule_matches.is_empty()
        {
            let selection_status = if selection_rule_matches
                .iter()
                .any(|row| row.status == "supported")
            {
                "supported"
            } else if selection_rule_matches
                .iter()
                .any(|row| row.status == "review_needed")
            {
                "review_needed"
            } else if selection_rule_matches
                .iter()
                .any(|row| row.status == "candidates_detected")
            {
                "candidates_detected"
            } else {
                "context_recorded"
            };
            let selection_label = match selection_status {
                "supported" => "Selection/complementation context supported".to_string(),
                "review_needed" => "Selection/complementation context requires review".to_string(),
                "candidates_detected" => {
                    "Selection/complementation candidates detected".to_string()
                }
                _ => "Selection context recorded".to_string(),
            };
            let selection_rationale = match selection_status {
                "supported" => format!(
                    "Recorded medium conditions match engine-owned selection/complementation rule(s): {}. Construct or helper context includes candidate(s): {}.",
                    selection_rule_matches
                        .iter()
                        .filter(|row| row.status == "supported")
                        .map(|row| row.label.clone())
                        .collect::<Vec<_>>()
                        .join(", "),
                    selection_candidates.iter().cloned().collect::<Vec<_>>().join(", ")
                ),
                "review_needed" => format!(
                    "Recorded medium conditions match engine-owned selection/complementation rule(s): {}. No matching construct or helper candidates were found yet.",
                    selection_rule_matches
                        .iter()
                        .filter(|row| row.status == "review_needed")
                        .map(|row| row.label.clone())
                        .collect::<Vec<_>>()
                        .join(", ")
                ),
                "candidates_detected" => format!(
                    "Construct or helper context matches engine-owned selection/complementation candidate(s): {}. No matching medium condition is recorded yet, so current deterministic reasoning keeps this as inspectable context rather than claiming a selection strategy.",
                    selection_candidates.iter().cloned().collect::<Vec<_>>().join(", ")
                ),
                _ => "Construct objective records medium/selection conditions, but none currently match the available engine-owned selection/complementation rules.".to_string(),
            };
            let mut selection_evidence_ids = medium_evidence_ids.clone();
            selection_evidence_ids.extend(selection_candidate_ids);
            if helper_backed_selection_present {
                selection_evidence_ids.extend(helper_evidence_ids.clone());
            }
            selection_evidence_ids.sort();
            selection_evidence_ids.dedup();
            facts.push(Self::construct_reasoning_build_fact(
                "fact_selection_context",
                "selection_context",
                selection_label,
                selection_rationale.clone(),
                selection_evidence_ids.clone(),
                json!({
                    "seq_id": seq_id,
                    "status": selection_status,
                    "medium_conditions": objective.medium_conditions.clone(),
                    "matching_medium_conditions": matching_medium_conditions,
                    "selection_candidates": selection_candidates.iter().cloned().collect::<Vec<_>>(),
                    "selection_rule_matches": selection_rule_matches,
                }),
            ));
            decisions.push(Self::construct_reasoning_build_decision(
                "decision_evaluate_selection_or_complementation_fit",
                "evaluate_selection_or_complementation_fit",
                "Evaluate Selection/Complementation Fit".to_string(),
                selection_rationale,
                selection_evidence_ids,
                vec!["fact_selection_context".to_string()],
                json!({
                    "status": selection_status,
                    "matching_rule_ids": matching_rule_ids,
                }),
            ));
        }

        if !variant_summaries.is_empty() {
            let variant_labels = variant_summaries
                .iter()
                .map(|row| row.label.clone())
                .collect::<Vec<_>>();
            let coding_prediction_count = variant_summaries
                .iter()
                .map(|row| row.coding_predictions.len())
                .sum::<usize>();
            let ambiguous_transcript_count = variant_summaries
                .iter()
                .filter(|row| row.transcript_context_status == "multi_transcript_ambiguous")
                .count();
            let mut variant_related_evidence_ids = variant_evidence_rows
                .iter()
                .map(|row| row.evidence_id.clone())
                .collect::<Vec<_>>();
            variant_related_evidence_ids.extend(evidence.iter().filter_map(|row| {
                (row.scope == EvidenceScope::SequenceSpan
                    && row.role != ConstructRole::Variant
                    && variant_evidence_rows.iter().any(|variant| {
                        Self::construct_reasoning_ranges_overlap(
                            variant.start_0based,
                            variant.end_0based_exclusive,
                            row.start_0based,
                            row.end_0based_exclusive,
                        )
                    }))
                .then(|| row.evidence_id.clone())
            }));
            variant_related_evidence_ids.sort();
            variant_related_evidence_ids.dedup();

            let variant_effect_status = if !variant_effect_tags.is_empty() {
                "effect_candidates_detected"
            } else {
                "variant_context_recorded"
            };
            let variant_effect_label = match variant_effect_status {
                "effect_candidates_detected" => "Variant effect candidates derived".to_string(),
                _ => "Variant context recorded".to_string(),
            };
            let variant_effect_rationale = if variant_effect_tags.is_empty() {
                format!(
                    "Detected {} variant marker(s) ({}) but the current deterministic overlap rules did not yet classify promoter/TFBS/CDS/UTR/splice effect candidates for them.",
                    variant_summaries.len(),
                    variant_labels.join(", ")
                )
            } else {
                let coding_text = if coding_prediction_count > 0 {
                    format!(
                        " {} codon-level coding prediction(s) were possible because allele-bearing VCF context was available.",
                        coding_prediction_count
                    )
                } else {
                    " Codon-level refinement remains absent for markers without explicit allele context.".to_string()
                };
                let transcript_text = if ambiguous_transcript_count > 0 {
                    format!(
                        " {} variant(s) remain transcript-ambiguous, so transcript-specific effect summaries are kept inspectable instead of being flattened away.",
                        ambiguous_transcript_count
                    )
                } else {
                    String::new()
                };
                format!(
                    "Variant placement against promoter/TFBS/CDS/UTR/splice annotations yields deterministic effect candidate tag(s): {}. Variant(s): {}.{}{}",
                    variant_effect_tags.join(", "),
                    variant_labels.join(", "),
                    coding_text,
                    transcript_text
                )
            };
            facts.push(Self::construct_reasoning_build_fact(
                "fact_variant_effect_context",
                "variant_effect_context",
                variant_effect_label,
                variant_effect_rationale.clone(),
                variant_related_evidence_ids.clone(),
                json!({
                    "seq_id": seq_id,
                    "status": variant_effect_status,
                    "variant_count": variant_summaries.len(),
                    "variant_labels": variant_labels,
                    "effect_tags": variant_effect_tags,
                    "variants": variant_summaries,
                }),
            ));
            decisions.push(Self::construct_reasoning_build_decision(
                "decision_evaluate_variant_effect_context",
                "evaluate_variant_effect_context",
                "Evaluate Variant Effect Context".to_string(),
                variant_effect_rationale,
                variant_related_evidence_ids.clone(),
                vec!["fact_variant_effect_context".to_string()],
                json!({
                    "variant_count": variant_summaries.len(),
                    "coding_prediction_count": coding_prediction_count,
                }),
            ));

            let variant_assay_status = if !suggested_variant_assay_ids.is_empty() {
                "assay_candidates_suggested"
            } else {
                "no_assay_rule_matched"
            };
            let variant_assay_label = match variant_assay_status {
                "assay_candidates_suggested" => "Promoter-design assays suggested".to_string(),
                _ => "Promoter-design assays require manual review".to_string(),
            };
            let variant_assay_rationale = if suggested_variant_assay_ids.is_empty() {
                "Variant markers are present, but the current deterministic assay suggester did not match a promoter/regulatory, coding, splice, or UTR assay family from the available overlap evidence.".to_string()
            } else {
                format!(
                    "Current deterministic assay-family rules suggest {} from the mapped variant context (for example promoter -> reporter, CDS -> paired expression, splice -> minigene, UTR -> reporter/translation compare).",
                    suggested_variant_assay_ids.join(", ")
                )
            };
            facts.push(Self::construct_reasoning_build_fact(
                "fact_variant_assay_context",
                "variant_assay_context",
                variant_assay_label,
                variant_assay_rationale.clone(),
                variant_related_evidence_ids.clone(),
                json!({
                    "seq_id": seq_id,
                    "status": variant_assay_status,
                    "suggested_assay_ids": suggested_variant_assay_ids,
                    "variant_labels": variant_summaries
                        .iter()
                        .map(|row| row.label.clone())
                        .collect::<Vec<_>>(),
                    "variant_assay_map": variant_summaries
                        .iter()
                        .map(|row| json!({
                            "label": row.label,
                            "suggested_assay_ids": row.suggested_assay_ids,
                            "effect_tags": row.effect_tags,
                        }))
                        .collect::<Vec<_>>(),
                }),
            ));
            decisions.push(Self::construct_reasoning_build_decision(
                "decision_suggest_variant_follow_up_assays",
                "suggest_variant_follow_up_assays",
                "Suggest Promoter-design Assays".to_string(),
                variant_assay_rationale,
                variant_related_evidence_ids,
                vec!["fact_variant_assay_context".to_string()],
                json!({
                    "variant_count": variant_summaries.len(),
                    "suggested_assay_ids": suggested_variant_assay_ids,
                }),
            ));
        }

        (facts, decisions)
    }

    fn assemble_construct_reasoning_graph(
        &mut self,
        seq_id: &str,
        objective_id: Option<&str>,
        graph_id: Option<&str>,
    ) -> Result<ConstructReasoningGraph, EngineError> {
        let seq_id = seq_id.trim();
        if seq_id.is_empty() {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: "seq_id cannot be empty".to_string(),
            });
        }
        let dna = self
            .state
            .sequences
            .get(seq_id)
            .ok_or_else(|| EngineError {
                code: ErrorCode::NotFound,
                message: format!("Sequence '{}' not found", seq_id),
            })?;
        let store = self.read_construct_reasoning_store();
        let existing_graph = graph_id
            .map(str::trim)
            .filter(|value| !value.is_empty())
            .and_then(|graph_id| store.graphs.get(graph_id))
            .filter(|graph| graph.seq_id.eq_ignore_ascii_case(seq_id))
            .cloned()
            .or_else(|| Self::select_construct_reasoning_graph_for_seq_id(&store, seq_id));
        let objective = if let Some(objective_id) = objective_id
            .map(str::trim)
            .filter(|value| !value.is_empty())
        {
            store
                .objectives
                .get(objective_id)
                .cloned()
                .ok_or_else(|| EngineError {
                    code: ErrorCode::NotFound,
                    message: format!(
                        "Construct objective '{}' not found in construct reasoning store",
                        objective_id
                    ),
                })?
        } else {
            Self::construct_reasoning_default_objective(seq_id, dna)
        };

        let mut evidence =
            Self::build_construct_reasoning_objective_context_evidence(seq_id, &objective);
        evidence.extend(self.build_construct_reasoning_sequence_evidence(seq_id, dna));

        let helper_interpretation = objective.helper_profile_id.as_ref().and_then(|profile_id| {
            Self::interpret_helper_genome(profile_id, None)
                .ok()
                .flatten()
        });
        let (facts, decisions) = Self::build_construct_reasoning_facts_and_decisions(
            seq_id,
            dna,
            &objective,
            &evidence,
            helper_interpretation.as_ref(),
        );
        let mut annotation_candidates = Self::construct_reasoning_build_annotation_candidates(
            seq_id, &evidence, &facts, &decisions,
        );
        Self::preserve_construct_reasoning_annotation_candidate_statuses(
            &mut annotation_candidates,
            existing_graph.as_ref(),
        );

        let graph = ConstructReasoningGraph {
            graph_id: graph_id
                .map(str::trim)
                .filter(|value| !value.is_empty())
                .map(|value| value.to_string())
                .unwrap_or_else(|| {
                    format!(
                        "construct_reasoning_{}_{}",
                        Self::normalize_id_token(seq_id),
                        Self::normalize_id_token(&objective.objective_id)
                    )
            }),
            seq_id: seq_id.to_string(),
            objective,
            generated_at_unix_ms: Self::now_unix_ms(),
            evidence,
            facts,
            decisions,
            annotation_candidates,
            notes: vec![
                "v1 deterministic reasoning graph includes construct-objective context evidence, sequence-backed restriction/annotation/variant evidence, interpreted growth-condition signals, host-route plus adapter-capture restriction/methylation review, and first hard-rule host/helper/selection/variant summary decisions."
                    .to_string(),
            ],
            ..ConstructReasoningGraph::default()
        };
        Ok(graph)
    }

    pub fn build_construct_reasoning_graph(
        &mut self,
        seq_id: &str,
        objective_id: Option<&str>,
        graph_id: Option<&str>,
    ) -> Result<ConstructReasoningGraph, EngineError> {
        let graph = self.assemble_construct_reasoning_graph(seq_id, objective_id, graph_id)?;
        self.upsert_construct_reasoning_graph(graph)
    }

    pub fn export_construct_reasoning_graph_json(
        &self,
        graph_id: &str,
        path: &str,
    ) -> Result<ConstructReasoningGraph, EngineError> {
        let graph = self.construct_reasoning_graph(graph_id)?;
        let text = serde_json::to_string_pretty(&graph).map_err(|e| EngineError {
            code: ErrorCode::Internal,
            message: format!(
                "Could not serialize construct reasoning graph '{}': {e}",
                graph.graph_id
            ),
        })?;
        std::fs::write(path, text).map_err(|e| EngineError {
            code: ErrorCode::Io,
            message: format!("Could not write construct reasoning graph to '{path}': {e}"),
        })?;
        Ok(graph)
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
        transcript_feature_id: Option<usize>,
        query_anchor_0based: Option<usize>,
        query_anchor_label: Option<String>,
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
            transcript_feature_id,
            query_anchor_0based,
            query_anchor_label,
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

    fn dotplot_transcript_feature_exon_local_starts(
        feature: &gb_io::seq::Feature,
    ) -> Vec<((usize, usize), usize, usize)> {
        let mut exon_ranges: Vec<(usize, usize)> = vec![];
        collect_location_ranges_usize(&feature.location, &mut exon_ranges);
        exon_ranges.sort_unstable_by(|a, b| a.0.cmp(&b.0).then(a.1.cmp(&b.1)));
        exon_ranges.dedup();
        exon_ranges.retain(|(start, end)| *end > *start);
        let is_reverse = feature_is_reverse(feature);
        let ordered = if is_reverse {
            exon_ranges.iter().copied().rev().collect::<Vec<_>>()
        } else {
            exon_ranges.clone()
        };
        let mut local_cursor = 0usize;
        let mut out = vec![];
        for (start_0based, end_0based_exclusive) in ordered {
            let exon_len = end_0based_exclusive.saturating_sub(start_0based);
            out.push((
                (start_0based.saturating_add(1), end_0based_exclusive),
                local_cursor,
                exon_len,
            ));
            local_cursor = local_cursor.saturating_add(exon_len);
        }
        out
    }

    fn build_dotplot_overlay_anchor_exons(
        &self,
        reference_seq_id: &str,
        query_specs: &[DotplotOverlayQuerySpec],
        query_series: &[DotplotQuerySeries],
    ) -> Vec<DotplotOverlayAnchorExon> {
        let Some(reference_dna) = self.state.sequences.get(reference_seq_id) else {
            return vec![];
        };
        let mut exon_support: BTreeMap<(usize, usize), Vec<DotplotOverlayAnchorSeriesSupport>> =
            BTreeMap::new();
        for (query, series) in query_specs.iter().zip(query_series.iter()) {
            let Some(transcript_feature_id) = query.transcript_feature_id else {
                continue;
            };
            let Some(feature) = reference_dna.features().get(transcript_feature_id) else {
                continue;
            };
            for ((start_1based, end_1based), transcript_local_start_0based, exon_len_bp) in
                Self::dotplot_transcript_feature_exon_local_starts(feature)
            {
                let transcript_local_end_0based =
                    transcript_local_start_0based.saturating_add(exon_len_bp);
                if transcript_local_start_0based < series.span_start_0based
                    || transcript_local_end_0based > series.span_end_0based
                {
                    continue;
                }
                exon_support
                    .entry((start_1based, end_1based))
                    .or_default()
                    .push(DotplotOverlayAnchorSeriesSupport {
                        series_id: series.series_id.clone(),
                        transcript_feature_id: Some(transcript_feature_id),
                        query_start_0based: transcript_local_start_0based
                            .saturating_sub(series.span_start_0based),
                    });
            }
        }
        exon_support
            .into_iter()
            .filter_map(|((start_1based, end_1based), supporting_series)| {
                if supporting_series.len() < 2 {
                    return None;
                }
                let max_query_start_0based = supporting_series
                    .iter()
                    .map(|support| support.query_start_0based)
                    .max()
                    .unwrap_or(0);
                Some(DotplotOverlayAnchorExon {
                    exon: DotplotOverlayAnchorExonRef {
                        start_1based,
                        end_1based,
                    },
                    support_series_count: supporting_series.len(),
                    max_query_start_0based,
                    supporting_series,
                })
            })
            .collect()
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
            None,
            None,
            None,
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
            overlay_anchor_exons: vec![],
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
