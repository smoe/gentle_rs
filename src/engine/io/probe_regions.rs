//! Affymetrix probe/probeset-region planning helpers.
//!
//! This first slice is intentionally plan-only: it checks user-supplied CEL,
//! metadata, annotation/library, platform, and local tool availability without
//! running R, APT, or any summarization backend.

use super::*;
use crate::publication_resources;

const PROBE_REGION_STAGE: &str = "plan_only";
const PROBE_REGION_IMPLEMENTATION_STATUS: &str =
    "stage_1_preflight_only_no_cel_summarization_backend";
const PROBE_REGION_OLIGO_HELPER: &str = "scripts/probe_regions_oligo.R";
const PROBE_REGION_TABLE_FILE: &str = "region_intensity_chrom_order.csv";
const PROBE_REGION_PROBE_TABLE_FILE: &str = "probe_intensity_chrom_order.csv";
const PROBE_REGION_SAMPLE_TABLE_FILE: &str = "sample_table.tsv";
const PROBE_REGION_PLAN_FILE: &str = "plan.json";
const PROBE_REGION_MATRIX_MANIFEST_FILE: &str = "normalized_feature_matrix_manifest.json";
const PROBE_REGION_PROVENANCE_FILE: &str = "provenance.json";
const PROBE_REGION_MATRIX_MANIFEST_SCHEMA: &str =
    "gentle.probe_region_normalized_matrix_manifest.v1";
const PROBE_REGION_BACKEND_PROVENANCE_SCHEMA: &str = "gentle.probe_region_backend_provenance.v1";
const PROBE_REGION_BACKEND_RUN_SCHEMA: &str = "gentle.probe_region_backend_run.v1";
const PROBE_REGION_FINGERPRINT_SHA1_MAX_BYTES: u64 = 10 * 1024 * 1024;
const CLARIOM_D_HUMAN_VENDOR_SUPPORT_DIR: &str =
    "data/resources/affymetrix/clariom_d_human_na36_hg38";
const CLARIOM_D_HUMAN_PROBESET_ZIP: &str = "Clariom_D_Human-na36-hg38-probeset-csv.zip";
const CLARIOM_D_HUMAN_TRANSCRIPT_ZIP: &str = "Clariom_D_Human.r1.na36.hg38.a1.transcript.csv.zip";
const CLARIOM_D_HUMAN_PROBESET_ZIP_ALIASES: &[&str] =
    &["TFS-Assets_LSG_Support-Files_Clariom_D_Human-na36-hg38-probeset-csv.zip"];
const CLARIOM_D_HUMAN_TRANSCRIPT_ZIP_ALIASES: &[&str] =
    &["TFS-Assets_LSG_Support-Files_Clariom_D_Human.r1.na36.hg38.a1.transcript.csv.zip"];

#[derive(Default)]
struct ProbeRegionTableSummary {
    row_count: usize,
    column_count: usize,
    feature_count: usize,
    transcript_cluster_count: usize,
    chromosome_count: usize,
    chromosomes: Vec<String>,
    gene_symbols: Vec<String>,
    sample_columns: Vec<String>,
    condition_summary_columns: Vec<String>,
    logfc_columns: Vec<String>,
    preview_rows: Vec<ProbeRegionOutputPreviewRow>,
    required_columns_missing: Vec<String>,
    warnings: Vec<String>,
}

#[derive(Default)]
struct ProbeRegionProbeTableSummary {
    row_count: usize,
    parent_feature_count: usize,
    required_columns_missing: Vec<String>,
}

#[derive(Clone)]
struct ProbeRegionPlotTrack {
    label: String,
    values: Vec<Option<f64>>,
}

#[derive(Clone)]
struct ProbeRegionPlotRow {
    chromosome: String,
    start_1based: usize,
    stop_1based: Option<usize>,
    probeset_or_region_id: String,
    gene_symbol: String,
}

struct ProbeRegionPlotData {
    rows: Vec<ProbeRegionPlotRow>,
    intensity_tracks: Vec<ProbeRegionPlotTrack>,
    logfc_tracks: Vec<ProbeRegionPlotTrack>,
    chromosome_count: usize,
    warnings: Vec<String>,
}

#[derive(Clone)]
struct ProbeRegionProjectionRow {
    chromosome: String,
    start_1based: usize,
    end_1based: usize,
    strand: Option<char>,
    feature_id: String,
    parent_feature_id: Option<String>,
    intensity_source: Option<String>,
    transcript_cluster_id: Option<String>,
    gene_symbol: Option<String>,
    logfc_values: BTreeMap<String, f64>,
}

#[derive(Clone)]
struct ProbeRegionProjectedEvidence {
    evidence_id: String,
    level: String,
    feature_id: String,
    parent_feature_id: Option<String>,
    intensity_source: Option<String>,
    chromosome: Option<String>,
    start_1based: Option<usize>,
    end_1based: Option<usize>,
    strand: Option<String>,
    logfc: Option<f64>,
    assembly_check: Option<String>,
    coordinate_frame: String,
    anchor_genome_id: Option<String>,
    anchor_start_1based: Option<usize>,
    anchor_end_1based: Option<usize>,
    anchor_strand: Option<char>,
    ranges_0based: Vec<(usize, usize)>,
}

struct ProbeRegionTranscriptModel {
    transcript_id: String,
    gene: Option<String>,
    label: Option<String>,
    strand: Option<String>,
    exons: Vec<ProbeRegionTranscriptExon>,
    exon_ranges_0based: Vec<(usize, usize)>,
    span_0based: Option<(usize, usize)>,
}

#[derive(Clone)]
struct ProbeRegionTranscriptExon {
    ordinal: usize,
    range_0based: (usize, usize),
}

#[derive(Default)]
struct ProbeRegionTranscriptEvidenceCounts {
    compatible: usize,
    constraining: usize,
    shared: usize,
    unique: usize,
    compatible_score: f64,
    constraining_score: f64,
    shared_score: f64,
    unique_score: f64,
}

#[derive(Default)]
struct ProbeRegionAptLibraryPlan {
    pgf_path: Option<String>,
    clf_path: Option<String>,
    mps_path: Option<String>,
    source_detail: Option<String>,
}

struct ProbeRegionAptAnnotationRow {
    chromosome: String,
    start_1based: usize,
    end_1based: usize,
    strand: String,
    feature_id: String,
    transcript_cluster_id: String,
    number_of_probes: String,
    gene_symbol: String,
}

#[derive(Default)]
struct ProbeRegionAptAnnotationTable {
    regions: HashMap<String, ProbeRegionAptAnnotationRow>,
    probes: Vec<ProbeRegionAptProbeAnnotationRow>,
}

struct ProbeRegionAptProbeAnnotationRow {
    chromosome: String,
    start_1based: usize,
    end_1based: usize,
    strand: String,
    parent_feature_id: String,
    transcript_cluster_id: String,
    gene_symbol: String,
    probe_id: String,
    x: String,
    y: String,
}

#[derive(Default)]
struct ProbeRegionAptMetadataSummary {
    path: String,
    sample_column: String,
    condition_column: String,
    conditions: Vec<ProbeRegionAptConditionGroup>,
    contrasts: Vec<ProbeRegionAptContrast>,
    warnings: Vec<String>,
}

struct ProbeRegionAptConditionGroup {
    label: String,
    column_label: String,
    sample_indices: Vec<usize>,
}

struct ProbeRegionAptContrast {
    label: String,
    numerator_index: usize,
    denominator_index: usize,
}

struct ProbeRegionAptProbeIntensityTable {
    path: String,
    probe_id_column: String,
    sample_columns: Vec<String>,
    value_header: Vec<String>,
    rows: BTreeMap<String, Vec<String>>,
    warnings: Vec<String>,
}

#[path = "probe_regions/import.rs"]
mod import;
#[path = "probe_regions/inspection.rs"]
mod inspection;
#[path = "probe_regions/interpretation.rs"]
mod interpretation;
#[path = "probe_regions/planning_backend.rs"]
mod planning_backend;
#[path = "probe_regions/projection.rs"]
mod projection;
#[path = "probe_regions/render_output_svg.rs"]
mod render_output_svg;
