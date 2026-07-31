//! Allele-aware RNA-read hash screen over transcript-coordinate variants.
//!
//! The screen builds allele-aware transcript hashes from a local transcript
//! FASTA plus either an explicit transcript-coordinate variant table or a
//! reviewed local VCF with an explicit transcript-coordinate map. It classifies
//! reads and fragments against phase-block-specific k-mer sets without
//! inventing phase for unphased genotypes. This is a deterministic evidence
//! screen: it reports sequence support for allele imbalance, but it does not
//! call biological significance.

use crate::target_rescue::{
    canonical_kmers_for_each, normalize_read_id, open_maybe_gz_reader, read_id_set_from_path,
    target_mapped_read_ids_from_sam, visit_fasta_records, visit_paired_read_records,
    visit_read_records,
};
use serde::{Deserialize, Serialize};
use std::{
    collections::{BTreeMap, BTreeSet, HashMap, HashSet},
    fs::{self, File},
    io::{BufRead, BufReader, BufWriter, Write},
    path::Path,
};

pub const ALLELE_HASH_SCREEN_SCHEMA: &str = "gentle.rna_allele_hash_screen.v2";
#[cfg(test)]
const ALLELE_HASH_SCREEN_SCHEMA_V1: &str = "gentle.rna_allele_hash_screen.v1";
const DEFAULT_KMER_LEN: usize = 21;
const DEFAULT_MIN_UNIQUE_KMER_HITS: u64 = 1;
const DEFAULT_MAX_INLINE_READ_CALLS: usize = 10_000;
const GENE_SYMBOL_TAGS: &[&str] = &["gene_symbol", "gene_name", "symbol", "gene"];

fn default_max_inline_read_calls() -> usize {
    DEFAULT_MAX_INLINE_READ_CALLS
}

/// Request for the standalone allele-aware RNA-read screen.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct AlleleHashScreenRequest {
    pub gene: String,
    pub transcript_fasta: String,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub variant_table: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub vcf: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub transcript_map: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub vcf_sample: Option<String>,
    pub read_files: Vec<String>,
    #[serde(default)]
    pub read_pairs: Vec<AlleleReadPairInput>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub read_id_allowlist: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub from_rna_report: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub salmon_unmapped_names: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub salmon_mappings_sam: Option<String>,
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub rna_report_reads: Vec<AlleleRnaReportReadInput>,
    pub out_dir: String,
    pub kmer_len: usize,
    pub min_unique_kmer_hits: u64,
    #[serde(default = "default_max_inline_read_calls")]
    pub max_inline_read_calls: usize,
}

impl Default for AlleleHashScreenRequest {
    fn default() -> Self {
        Self {
            gene: String::new(),
            transcript_fasta: String::new(),
            variant_table: None,
            vcf: None,
            transcript_map: None,
            vcf_sample: None,
            read_files: Vec::new(),
            read_pairs: Vec::new(),
            read_id_allowlist: None,
            from_rna_report: None,
            salmon_unmapped_names: None,
            salmon_mappings_sam: None,
            rna_report_reads: Vec::new(),
            out_dir: String::new(),
            kmer_len: DEFAULT_KMER_LEN,
            min_unique_kmer_hits: DEFAULT_MIN_UNIQUE_KMER_HITS,
            max_inline_read_calls: DEFAULT_MAX_INLINE_READ_CALLS,
        }
    }
}

/// Top-level report emitted by the allele-aware RNA-read screen.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct AlleleHashScreenReport {
    pub schema: String,
    pub gene: String,
    pub phase_mode: AlleleHashPhaseMode,
    pub params: AlleleHashScreenParams,
    pub transcript_count: usize,
    pub variant_count: usize,
    pub read_count_total: u64,
    pub read_count_selected: u64,
    #[serde(default)]
    pub fragment_count_total: u64,
    #[serde(default)]
    pub fragment_count_selected: u64,
    #[serde(default)]
    pub evidence_observation_count_selected: u64,
    #[serde(default)]
    pub read_calls_truncated: bool,
    #[serde(default)]
    pub phase_blocks: Vec<AllelePhaseBlockSummary>,
    pub output_files: AlleleHashScreenOutputFiles,
    pub haplotype_fastas: Vec<HaplotypeFastaReport>,
    pub transcript_summaries: Vec<AlleleTranscriptSummary>,
    pub variant_summaries: Vec<AlleleVariantSummary>,
    pub classification_counts: BTreeMap<String, u64>,
    #[serde(default)]
    pub source_provenance: Vec<AlleleReadSourceProvenance>,
    pub reads: Vec<AlleleReadCall>,
    pub warnings: Vec<String>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct AlleleHashScreenParams {
    pub transcript_fasta: String,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub variant_table: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub vcf: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub transcript_map: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub vcf_sample: Option<String>,
    pub read_files: Vec<String>,
    #[serde(default)]
    pub read_pairs: Vec<AlleleReadPairInput>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub read_id_allowlist: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub from_rna_report: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub salmon_unmapped_names: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub salmon_mappings_sam: Option<String>,
    pub kmer_len: usize,
    pub min_unique_kmer_hits: u64,
    #[serde(default = "default_max_inline_read_calls")]
    pub max_inline_read_calls: usize,
}

#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct AlleleReadPairInput {
    pub read1: String,
    pub read2: String,
}

/// One target-gene read resolved from a persisted RNA-read report by the engine.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct AlleleRnaReportReadInput {
    pub read_id: String,
    pub sequence: String,
    pub report_id: String,
    pub record_index: usize,
}

/// Input-basis tags attached to each selected allele-screen observation.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum AlleleReadSourceOrigin {
    ExplicitFile,
    RnaReportTargetMapped,
    SalmonUnassigned,
    SalmonTargetMapped,
}

impl AlleleReadSourceOrigin {
    pub fn as_str(self) -> &'static str {
        match self {
            Self::ExplicitFile => "explicit_file",
            Self::RnaReportTargetMapped => "rna_report_target_mapped",
            Self::SalmonUnassigned => "salmon_unassigned",
            Self::SalmonTargetMapped => "salmon_target_mapped",
        }
    }
}

/// Selected evidence counts by source tag. Counts may overlap across tags.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct AlleleReadSourceProvenance {
    pub origin: AlleleReadSourceOrigin,
    pub evidence_observation_count: u64,
    pub input_read_count: u64,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct AlleleHashScreenOutputFiles {
    pub report_json: String,
    pub reads_tsv: String,
    pub reference_fasta: String,
    pub hap1_fasta: String,
    pub hap2_fasta: String,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct HaplotypeFastaReport {
    pub haplotype: String,
    pub path: String,
    pub transcript_count: usize,
    #[serde(default)]
    pub inferred: bool,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum AlleleHashPhaseMode {
    Phased,
    PhasedBlocks,
    Mixed,
    UnphasedAlleleLevelOnly,
}

impl AlleleHashPhaseMode {
    pub fn as_str(self) -> &'static str {
        match self {
            Self::Phased => "phased",
            Self::PhasedBlocks => "phased_blocks",
            Self::Mixed => "mixed",
            Self::UnphasedAlleleLevelOnly => "unphased_allele_level_only",
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum AlleleReadClassification {
    Hap1,
    Hap2,
    Alternate,
    ReferenceOnly,
    Ambiguous,
    Uninformative,
    OffTarget,
}

impl AlleleReadClassification {
    pub fn as_str(self) -> &'static str {
        match self {
            Self::Hap1 => "hap1",
            Self::Hap2 => "hap2",
            Self::Alternate => "alternate",
            Self::ReferenceOnly => "reference_only",
            Self::Ambiguous => "ambiguous",
            Self::Uninformative => "uninformative",
            Self::OffTarget => "off_target",
        }
    }
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct AlleleReadCall {
    pub read_id: String,
    pub source_file: String,
    #[serde(default)]
    pub source_origins: Vec<AlleleReadSourceOrigin>,
    #[serde(default)]
    pub evidence_unit: AlleleEvidenceUnit,
    #[serde(default = "default_input_read_count")]
    pub input_read_count: u8,
    pub classification: AlleleReadClassification,
    pub read_length: usize,
    pub read_kmer_count: u64,
    pub reference_hits: u64,
    pub hap1_hits: u64,
    pub hap2_hits: u64,
    pub hap1_unique_hits: u64,
    pub hap2_unique_hits: u64,
    #[serde(default)]
    pub alternate_unique_hits: u64,
    pub reference_only_hits: u64,
    pub matched_transcripts: Vec<String>,
    pub supporting_variants: Vec<String>,
    #[serde(default)]
    pub phase_block_calls: Vec<AllelePhaseBlockReadCall>,
}

fn default_input_read_count() -> u8 {
    1
}

#[derive(Debug, Clone, Copy, Default, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum AlleleEvidenceUnit {
    #[default]
    Read,
    Fragment,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum AllelePhaseBlockStatus {
    Phased,
    Unphased,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum AllelePhaseBlockClassification {
    Haplotype1,
    Haplotype2,
    Alternate,
    Reference,
    Ambiguous,
    Uninformative,
}

impl AllelePhaseBlockClassification {
    fn as_str(self) -> &'static str {
        match self {
            Self::Haplotype1 => "haplotype1",
            Self::Haplotype2 => "haplotype2",
            Self::Alternate => "alternate",
            Self::Reference => "reference",
            Self::Ambiguous => "ambiguous",
            Self::Uninformative => "uninformative",
        }
    }
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct AllelePhaseBlockReadCall {
    pub block_id: String,
    pub classification: AllelePhaseBlockClassification,
    pub reference_unique_hits: u64,
    pub hap1_unique_hits: u64,
    pub hap2_unique_hits: u64,
    pub alternate_unique_hits: u64,
    pub variant_ids: Vec<String>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct AllelePhaseBlockSummary {
    pub block_id: String,
    pub transcript_id: String,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub phase_set: Option<String>,
    pub status: AllelePhaseBlockStatus,
    pub variant_ids: Vec<String>,
    pub informative_reads: u64,
    pub hap1_reads: u64,
    pub hap2_reads: u64,
    pub alternate_reads: u64,
    pub reference_reads: u64,
    pub ambiguous_reads: u64,
    pub uninformative_reads: u64,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct AlleleTranscriptSummary {
    pub transcript_id: String,
    pub length_nt: usize,
    pub variant_count: usize,
    pub informative_reads: u64,
    pub hap1_reads: u64,
    pub hap2_reads: u64,
    #[serde(default)]
    pub alternate_reads: u64,
    pub reference_only_reads: u64,
    pub ambiguous_reads: u64,
    pub uninformative_reads: u64,
    pub off_target_reads: u64,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub imbalance_ratio_hap1_over_informative: Option<f64>,
    pub coverage_blind_spots: Vec<String>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct AlleleVariantSummary {
    pub variant_id: String,
    pub transcript_id: String,
    pub cdna_pos_1based: usize,
    pub ref_allele: String,
    pub alt_allele: String,
    pub genotype: String,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub phase_set: Option<String>,
    pub informative_reads: u64,
    pub hap1_reads: u64,
    pub hap2_reads: u64,
    #[serde(default)]
    pub alternate_reads: u64,
    pub reference_only_reads: u64,
    pub ambiguous_reads: u64,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub imbalance_ratio_hap1_over_informative: Option<f64>,
    pub coverage_blind_spots: Vec<String>,
}

#[derive(Debug, Clone)]
struct TranscriptRecord {
    id: String,
    header: String,
    reference: String,
    hap1: Vec<u8>,
    hap2: Vec<u8>,
}

#[derive(Debug, Clone)]
struct VariantRecord {
    variant_id: String,
    transcript_id: String,
    cdna_pos_1based: usize,
    ref_allele: String,
    alt_allele: String,
    genotype: String,
    phase_set: Option<String>,
    hap1_alt: bool,
    hap2_alt: bool,
    phased: bool,
}

#[derive(Debug, Clone, Default)]
struct KmerIndex {
    reference: HashSet<u64>,
    hap1: HashSet<u64>,
    hap2: HashSet<u64>,
    alternate: HashSet<u64>,
    reference_only: HashSet<u64>,
    hap1_unique: HashSet<u64>,
    hap2_unique: HashSet<u64>,
    alternate_unique: HashSet<u64>,
}

#[derive(Debug, Clone)]
struct TranscriptKmerIndex {
    transcript_id: String,
    all: HashSet<u64>,
}

#[derive(Debug, Clone)]
struct PhaseBlockIndex {
    block_id: String,
    transcript_id: String,
    phase_set: Option<String>,
    status: AllelePhaseBlockStatus,
    variant_ids: Vec<String>,
    variant_kmers: BTreeMap<String, HashSet<u64>>,
    reference_only: HashSet<u64>,
    hap1_unique: HashSet<u64>,
    hap2_unique: HashSet<u64>,
    alternate_unique: HashSet<u64>,
    all: HashSet<u64>,
}

#[derive(Debug, Clone, Default)]
struct CountSummary {
    informative_reads: u64,
    hap1_reads: u64,
    hap2_reads: u64,
    alternate_reads: u64,
    reference_only_reads: u64,
    ambiguous_reads: u64,
    uninformative_reads: u64,
    off_target_reads: u64,
}

/// Run the allele-aware RNA-read hash screen and write FASTA/TSV/JSON outputs.
pub fn run_allele_hash_screen(
    request: AlleleHashScreenRequest,
) -> Result<AlleleHashScreenReport, String> {
    validate_request(&request)?;
    let mut warnings = Vec::<String>::new();
    let mut transcripts = load_gene_transcripts(&request.transcript_fasta, &request.gene)?;
    let variants = match (
        request.variant_table.as_deref(),
        request.vcf.as_deref(),
        request.transcript_map.as_deref(),
    ) {
        (Some(variant_table), None, _) => load_variant_table(variant_table, &transcripts)?,
        (None, Some(vcf), Some(transcript_map)) => load_vcf_variants(
            vcf,
            transcript_map,
            request.vcf_sample.as_deref(),
            &transcripts,
            &mut warnings,
        )?,
        _ => {
            return Err(
                "allele-hash-screen requires either --variant-table PATH or --vcf PATH with --transcript-map PATH"
                    .to_string(),
            );
        }
    };
    let phase_block_indexes = build_phase_block_indexes(&transcripts, &variants, request.kmer_len)?;
    let phased_block_count = phase_block_indexes
        .iter()
        .filter(|block| block.status == AllelePhaseBlockStatus::Phased)
        .count();
    let unphased_block_count = phase_block_indexes.len() - phased_block_count;
    let phase_mode = match (phased_block_count, unphased_block_count) {
        (1, 0) => AlleleHashPhaseMode::Phased,
        (n, 0) if n > 1 => AlleleHashPhaseMode::PhasedBlocks,
        (0, _) => AlleleHashPhaseMode::UnphasedAlleleLevelOnly,
        _ => AlleleHashPhaseMode::Mixed,
    };
    let has_global_haplotype = phase_mode == AlleleHashPhaseMode::Phased;
    if unphased_block_count > 0 {
        warnings.push(
            "At least one variant is unphased; those rows are screened as reference-versus-alternate evidence and are not materialized as hap1/hap2."
                .to_string(),
        );
    }
    if phase_block_indexes.len() > 1 {
        warnings.push(format!(
            "Variants resolve to {} disconnected phase block(s); GENtle reports block-local evidence and does not fabricate a gene-wide haplotype.",
            phase_block_indexes.len()
        ));
    }
    if has_global_haplotype {
        apply_variants_to_haplotypes(&mut transcripts, &variants)?;
    } else {
        warnings.push(
            "Legacy hap1/hap2 FASTA paths contain reference copies because no single global phase block is supported; use phase_blocks for allele interpretation."
                .to_string(),
        );
    }
    let out_dir = Path::new(&request.out_dir);
    fs::create_dir_all(out_dir).map_err(|e| {
        format!(
            "Could not create output directory '{}': {e}",
            out_dir.display()
        )
    })?;
    let output_files = output_files(out_dir, &request.gene);
    write_haplotype_fastas(&transcripts, &output_files)?;

    let kmer_index = build_kmer_index(&transcripts, &phase_block_indexes, request.kmer_len);
    let transcript_indexes =
        build_transcript_kmer_indexes(&transcripts, &phase_block_indexes, request.kmer_len);
    let allowlist = request
        .read_id_allowlist
        .as_deref()
        .map(read_id_set_from_path)
        .transpose()?;
    let salmon_unassigned_ids = request
        .salmon_unmapped_names
        .as_deref()
        .map(read_id_set_from_path)
        .transpose()?;
    let target_transcript_ids = transcripts
        .iter()
        .flat_map(|transcript| {
            let unversioned = transcript
                .id
                .split_once('.')
                .map_or(transcript.id.as_str(), |(id, _)| id)
                .to_string();
            [transcript.id.clone(), unversioned]
        })
        .collect::<HashSet<_>>();
    let salmon_target_mapped_ids = request
        .salmon_mappings_sam
        .as_deref()
        .map(|path| target_mapped_read_ids_from_sam(path, &target_transcript_ids))
        .transpose()?;
    let salmon_selector_active =
        salmon_unassigned_ids.is_some() || salmon_target_mapped_ids.is_some();
    let rna_report_read_ids = request
        .rna_report_reads
        .iter()
        .map(|read| normalize_read_id(&read.read_id))
        .collect::<HashSet<_>>();
    let mut read_count_total = 0u64;
    let mut read_count_selected = 0u64;
    let mut fragment_count_total = 0u64;
    let mut fragment_count_selected = 0u64;
    let mut evidence_observation_count_selected = 0u64;
    let mut calls = Vec::<AlleleReadCall>::new();
    let mut transcript_counts = BTreeMap::<String, CountSummary>::new();
    let mut variant_counts = BTreeMap::<String, CountSummary>::new();
    let mut phase_block_counts = BTreeMap::<String, CountSummary>::new();
    let mut classification_counts = empty_classification_counts();
    let mut source_provenance_counts = BTreeMap::<AlleleReadSourceOrigin, (u64, u64)>::new();
    let mut duplicate_explicit_report_records = 0u64;
    let mut reads_writer = BufWriter::new(File::create(&output_files.reads_tsv).map_err(|e| {
        format!(
            "Could not create read TSV '{}': {e}",
            output_files.reads_tsv
        )
    })?);
    write_read_calls_header(&mut reads_writer, &output_files.reads_tsv)?;
    {
        let mut process_record = |record: crate::target_rescue::SequenceRecord,
                                  source_file: &str,
                                  evidence_unit: AlleleEvidenceUnit,
                                  input_read_count: u8,
                                  mut source_origins: Vec<AlleleReadSourceOrigin>,
                                  selected_by_source: bool|
         -> Result<(), String> {
            read_count_total = read_count_total.saturating_add(u64::from(input_read_count));
            fragment_count_total = fragment_count_total.saturating_add(1);
            if !selected_by_source {
                return Ok(());
            }
            if let Some(allowlist) = &allowlist
                && !allowlist.contains(&record.id)
            {
                return Ok(());
            }
            source_origins.sort_unstable();
            source_origins.dedup();
            read_count_selected = read_count_selected.saturating_add(u64::from(input_read_count));
            fragment_count_selected = fragment_count_selected.saturating_add(1);
            evidence_observation_count_selected =
                evidence_observation_count_selected.saturating_add(1);
            for origin in &source_origins {
                let counts = source_provenance_counts.entry(*origin).or_default();
                counts.0 = counts.0.saturating_add(1);
                counts.1 = counts.1.saturating_add(u64::from(input_read_count));
            }
            let call = classify_read(
                ReadClassificationInput {
                    read_id: record.id,
                    source_file,
                    source_origins,
                    evidence_unit,
                    input_read_count,
                    sequence: &record.sequence,
                    source_length: record.source_length,
                },
                ReadClassificationContext {
                    k: request.kmer_len,
                    min_unique_hits: request.min_unique_kmer_hits,
                    index: &kmer_index,
                    transcript_indexes: &transcript_indexes,
                    phase_blocks: &phase_block_indexes,
                },
            );
            write_read_call_tsv(&mut reads_writer, &output_files.reads_tsv, &call)?;
            accumulate_read_call(
                &call,
                &mut transcript_counts,
                &mut variant_counts,
                &mut phase_block_counts,
                &mut classification_counts,
            );
            if calls.len() < request.max_inline_read_calls {
                calls.push(call);
            }
            Ok(())
        };
        for report_read in &request.rna_report_reads {
            let read_id = normalize_read_id(&report_read.read_id);
            let mut source_origins = vec![AlleleReadSourceOrigin::RnaReportTargetMapped];
            if salmon_unassigned_ids
                .as_ref()
                .is_some_and(|ids| ids.contains(&read_id))
            {
                source_origins.push(AlleleReadSourceOrigin::SalmonUnassigned);
            }
            if salmon_target_mapped_ids
                .as_ref()
                .is_some_and(|ids| ids.contains(&read_id))
            {
                source_origins.push(AlleleReadSourceOrigin::SalmonTargetMapped);
            }
            let source_file = format!("rna-report:{}", report_read.report_id);
            process_record(
                crate::target_rescue::SequenceRecord {
                    id: read_id,
                    source_length: report_read.sequence.len(),
                    sequence: report_read.sequence.clone(),
                },
                &source_file,
                AlleleEvidenceUnit::Read,
                1,
                source_origins,
                true,
            )?;
        }
        for read_file in &request.read_files {
            visit_read_records(read_file, |record| {
                if rna_report_read_ids.contains(&record.id) {
                    duplicate_explicit_report_records =
                        duplicate_explicit_report_records.saturating_add(1);
                    return Ok(());
                }
                let mut source_origins = vec![AlleleReadSourceOrigin::ExplicitFile];
                let mut selected_by_source = !salmon_selector_active;
                if salmon_unassigned_ids
                    .as_ref()
                    .is_some_and(|ids| ids.contains(&record.id))
                {
                    source_origins.push(AlleleReadSourceOrigin::SalmonUnassigned);
                    selected_by_source = true;
                }
                if salmon_target_mapped_ids
                    .as_ref()
                    .is_some_and(|ids| ids.contains(&record.id))
                {
                    source_origins.push(AlleleReadSourceOrigin::SalmonTargetMapped);
                    selected_by_source = true;
                }
                process_record(
                    record,
                    read_file,
                    AlleleEvidenceUnit::Read,
                    1,
                    source_origins,
                    selected_by_source,
                )
            })?;
        }
        for read_pair in &request.read_pairs {
            let source_file = format!("{};{}", read_pair.read1, read_pair.read2);
            visit_paired_read_records(&read_pair.read1, &read_pair.read2, |record| {
                if rna_report_read_ids.contains(&record.id) {
                    duplicate_explicit_report_records =
                        duplicate_explicit_report_records.saturating_add(1);
                    return Ok(());
                }
                let mut source_origins = vec![AlleleReadSourceOrigin::ExplicitFile];
                let mut selected_by_source = !salmon_selector_active;
                if salmon_unassigned_ids
                    .as_ref()
                    .is_some_and(|ids| ids.contains(&record.id))
                {
                    source_origins.push(AlleleReadSourceOrigin::SalmonUnassigned);
                    selected_by_source = true;
                }
                if salmon_target_mapped_ids
                    .as_ref()
                    .is_some_and(|ids| ids.contains(&record.id))
                {
                    source_origins.push(AlleleReadSourceOrigin::SalmonTargetMapped);
                    selected_by_source = true;
                }
                process_record(
                    record,
                    &source_file,
                    AlleleEvidenceUnit::Fragment,
                    2,
                    source_origins,
                    selected_by_source,
                )
            })?;
        }
    }
    reads_writer
        .flush()
        .map_err(|e| format!("Could not flush read TSV '{}': {e}", output_files.reads_tsv))?;

    let mut transcript_summaries =
        summarize_transcripts(&transcripts, &variants, &transcript_counts);
    let mut variant_summaries = summarize_variants(&variants, &variant_counts);
    for summary in &mut transcript_summaries {
        if summary.informative_reads == 0 {
            summary
                .coverage_blind_spots
                .push("no haplotype-informative reads overlapped this transcript".to_string());
        }
    }
    for summary in &mut variant_summaries {
        if summary.informative_reads == 0 {
            summary
                .coverage_blind_spots
                .push("no reads carried haplotype-informative k-mers for this variant".to_string());
        }
    }
    let phase_blocks = summarize_phase_blocks(&phase_block_indexes, &phase_block_counts);
    let read_calls_truncated = evidence_observation_count_selected > calls.len() as u64;
    if read_calls_truncated {
        warnings.push(format!(
            "Inline JSON read calls were capped at {}; the complete selected-read table is in '{}'.",
            request.max_inline_read_calls, output_files.reads_tsv
        ));
    }
    if duplicate_explicit_report_records > 0 {
        warnings.push(format!(
            "Skipped {duplicate_explicit_report_records} explicit read record(s) whose normalized ids were already sourced from RNA-read report '{}'.",
            request.from_rna_report.as_deref().unwrap_or("unknown")
        ));
    }
    let source_provenance = source_provenance_counts
        .into_iter()
        .map(|(origin, (evidence_observation_count, input_read_count))| {
            AlleleReadSourceProvenance {
                origin,
                evidence_observation_count,
                input_read_count,
            }
        })
        .collect::<Vec<_>>();
    let report = AlleleHashScreenReport {
        schema: ALLELE_HASH_SCREEN_SCHEMA.to_string(),
        gene: request.gene.clone(),
        phase_mode,
        params: AlleleHashScreenParams {
            transcript_fasta: request.transcript_fasta.clone(),
            variant_table: request.variant_table.clone(),
            vcf: request.vcf.clone(),
            transcript_map: request.transcript_map.clone(),
            vcf_sample: request.vcf_sample.clone(),
            read_files: request.read_files.clone(),
            read_pairs: request.read_pairs.clone(),
            read_id_allowlist: request.read_id_allowlist.clone(),
            from_rna_report: request.from_rna_report.clone(),
            salmon_unmapped_names: request.salmon_unmapped_names.clone(),
            salmon_mappings_sam: request.salmon_mappings_sam.clone(),
            kmer_len: request.kmer_len,
            min_unique_kmer_hits: request.min_unique_kmer_hits,
            max_inline_read_calls: request.max_inline_read_calls,
        },
        transcript_count: transcripts.len(),
        variant_count: variants.len(),
        read_count_total,
        read_count_selected,
        fragment_count_total,
        fragment_count_selected,
        evidence_observation_count_selected,
        read_calls_truncated,
        phase_blocks,
        output_files: output_files.clone(),
        haplotype_fastas: vec![
            HaplotypeFastaReport {
                haplotype: "reference".to_string(),
                path: output_files.reference_fasta.clone(),
                transcript_count: transcripts.len(),
                inferred: false,
            },
            HaplotypeFastaReport {
                haplotype: "hap1".to_string(),
                path: output_files.hap1_fasta.clone(),
                transcript_count: transcripts.len(),
                inferred: has_global_haplotype,
            },
            HaplotypeFastaReport {
                haplotype: "hap2".to_string(),
                path: output_files.hap2_fasta.clone(),
                transcript_count: transcripts.len(),
                inferred: has_global_haplotype,
            },
        ],
        transcript_summaries,
        variant_summaries,
        classification_counts,
        source_provenance,
        reads: calls,
        warnings,
    };
    write_json_file(&output_files.report_json, &report)?;
    Ok(report)
}

fn validate_request(request: &AlleleHashScreenRequest) -> Result<(), String> {
    if request.gene.trim().is_empty() {
        return Err("allele-hash-screen requires --gene GENE".to_string());
    }
    if request.transcript_fasta.trim().is_empty() {
        return Err("allele-hash-screen requires --transcript-fasta PATH".to_string());
    }
    if request.variant_table.is_some() && request.vcf.is_some() {
        return Err("allele-hash-screen accepts --variant-table or --vcf, not both".to_string());
    }
    if request.variant_table.is_none() && request.vcf.is_none() {
        return Err(
            "allele-hash-screen requires --variant-table PATH or --vcf PATH with --transcript-map PATH"
                .to_string(),
        );
    }
    if request.vcf.is_some() && request.transcript_map.is_none() {
        return Err("allele-hash-screen --vcf requires --transcript-map PATH".to_string());
    }
    if request.from_rna_report.is_some() && request.rna_report_reads.is_empty() {
        return Err(
            "allele-hash-screen --from-rna-report must be resolved through a project-backed engine and contain at least one target-mapped read"
                .to_string(),
        );
    }
    if (request.salmon_unmapped_names.is_some() || request.salmon_mappings_sam.is_some())
        && request.read_files.is_empty()
        && request.read_pairs.is_empty()
    {
        return Err(
            "allele-hash-screen Salmon selectors require at least one --read-file PATH or --read-pair R1,R2 because Salmon name/mapping files do not contain every selected read sequence"
                .to_string(),
        );
    }
    if request.read_files.is_empty()
        && request.read_pairs.is_empty()
        && request.rna_report_reads.is_empty()
    {
        return Err(
            "allele-hash-screen requires --from-rna-report REPORT_ID, --read-file PATH, or --read-pair R1,R2"
                .to_string(),
        );
    }
    if request.out_dir.trim().is_empty() {
        return Err("allele-hash-screen requires --out OUT_DIR".to_string());
    }
    if !(1..=31).contains(&request.kmer_len) {
        return Err(format!(
            "Invalid --kmer-len value '{}': expected 1..=31",
            request.kmer_len
        ));
    }
    if request.min_unique_kmer_hits == 0 {
        return Err("Invalid --min-unique-kmer-hits value '0': expected at least 1".to_string());
    }
    Ok(())
}

fn load_gene_transcripts(path: &str, gene: &str) -> Result<Vec<TranscriptRecord>, String> {
    let wanted = normalize_gene(gene);
    let mut records = Vec::new();
    let mut saw_gene_tag = false;
    visit_fasta_records(path, |header, sequence| {
        let gene_symbol = parse_gene_symbol_from_header(header);
        if gene_symbol.is_some() {
            saw_gene_tag = true;
        }
        if gene_symbol
            .as_deref()
            .is_some_and(|symbol| normalize_gene(symbol) != wanted)
        {
            return Ok(());
        }
        if gene_symbol.is_none()
            || gene_symbol
                .as_deref()
                .is_some_and(|symbol| normalize_gene(symbol) == wanted)
        {
            let id = transcript_id_from_header(header);
            if !id.is_empty() {
                records.push(TranscriptRecord {
                    id,
                    header: header.to_string(),
                    reference: normalize_sequence(sequence),
                    hap1: normalize_sequence(sequence).into_bytes(),
                    hap2: normalize_sequence(sequence).into_bytes(),
                });
            }
        }
        Ok(())
    })?;
    if records.is_empty() {
        if saw_gene_tag {
            return Err(format!(
                "No transcripts for gene '{gene}' were found in '{path}'"
            ));
        }
        return Err(format!("No transcript records were found in '{path}'"));
    }
    Ok(records)
}

fn transcript_id_from_header(header: &str) -> String {
    header
        .split_ascii_whitespace()
        .next()
        .unwrap_or("")
        .trim_start_matches('>')
        .to_string()
}

fn normalize_sequence(raw: &str) -> String {
    raw.bytes()
        .filter(|byte| !byte.is_ascii_whitespace())
        .map(|byte| match byte.to_ascii_uppercase() {
            b'U' => 'T',
            b'A' => 'A',
            b'C' => 'C',
            b'G' => 'G',
            b'T' => 'T',
            _ => 'N',
        })
        .collect()
}

fn parse_gene_symbol_from_header(header: &str) -> Option<String> {
    let mut values = HashMap::<String, String>::new();
    for token in header.split_ascii_whitespace() {
        let Some((key, value)) = token.split_once(':') else {
            continue;
        };
        values
            .entry(key.to_ascii_lowercase())
            .or_insert_with(|| value.trim().to_string());
    }
    for tag in GENE_SYMBOL_TAGS {
        if let Some(value) = values.get(*tag)
            && !value.trim().is_empty()
        {
            return Some(value.trim().to_string());
        }
    }
    None
}

fn normalize_gene(raw: &str) -> String {
    raw.trim().to_ascii_uppercase()
}

fn load_variant_table(
    path: &str,
    transcripts: &[TranscriptRecord],
) -> Result<Vec<VariantRecord>, String> {
    let transcript_ids = transcripts
        .iter()
        .map(|record| record.id.as_str())
        .collect::<HashSet<_>>();
    let file =
        File::open(path).map_err(|e| format!("Could not open variant table '{path}': {e}"))?;
    let reader = BufReader::new(file);
    let mut header = Vec::<String>::new();
    let mut variants = Vec::<VariantRecord>::new();
    for (line_idx, line) in reader.lines().enumerate() {
        let line = line.map_err(|e| format!("Could not read variant table '{path}': {e}"))?;
        let trimmed = line.trim();
        if trimmed.is_empty() || trimmed.starts_with('#') {
            continue;
        }
        if header.is_empty() {
            header = trimmed
                .split('\t')
                .map(|value| value.trim().to_ascii_lowercase())
                .collect();
            continue;
        }
        let values = trimmed.split('\t').collect::<Vec<_>>();
        let row = header
            .iter()
            .zip(values.iter().copied())
            .map(|(key, value)| (key.as_str(), value.trim()))
            .collect::<HashMap<_, _>>();
        let transcript_id = required_table_value(&row, "transcript_id", path, line_idx + 1)?;
        if !transcript_ids.contains(transcript_id) {
            continue;
        }
        let pos_raw = required_table_value(&row, "cdna_pos_1based", path, line_idx + 1)?;
        let cdna_pos_1based = pos_raw.parse::<usize>().map_err(|e| {
            format!(
                "Invalid cdna_pos_1based '{pos_raw}' in variant table '{path}' line {}: {e}",
                line_idx + 1
            )
        })?;
        let ref_allele =
            required_table_value(&row, "ref", path, line_idx + 1)?.to_ascii_uppercase();
        let alt_allele =
            required_table_value(&row, "alt", path, line_idx + 1)?.to_ascii_uppercase();
        if ref_allele.len() != 1 || alt_allele.len() != 1 {
            return Err(format!(
                "Variant table '{path}' line {} is not an SNV-like row: ref='{}' alt='{}'",
                line_idx + 1,
                ref_allele,
                alt_allele
            ));
        }
        let genotype = required_table_value(&row, "genotype", path, line_idx + 1)?.to_string();
        let assignment = parse_genotype_assignment(&genotype).map_err(|detail| {
            format!(
                "Invalid genotype '{genotype}' in variant table '{path}' line {}: {detail}",
                line_idx + 1
            )
        })?;
        let variant_id = row
            .get("variant_id")
            .or_else(|| row.get("id"))
            .map(|value| value.trim())
            .filter(|value| !value.is_empty())
            .map(str::to_string)
            .unwrap_or_else(|| {
                format!("{transcript_id}:{cdna_pos_1based}:{ref_allele}>{alt_allele}")
            });
        variants.push(VariantRecord {
            variant_id,
            transcript_id: transcript_id.to_string(),
            cdna_pos_1based,
            ref_allele,
            alt_allele,
            genotype,
            phase_set: row
                .get("phase_set")
                .map(|value| value.trim().to_string())
                .filter(|value| !value.is_empty() && value != "."),
            hap1_alt: assignment.0,
            hap2_alt: assignment.1,
            phased: assignment.2,
        });
    }
    if variants.is_empty() {
        return Err(format!(
            "Variant table '{path}' did not contain transcript-coordinate variants for the selected transcripts"
        ));
    }
    Ok(variants)
}

#[derive(Debug, Clone)]
struct TranscriptCoordinateProjection {
    transcript_id: String,
    cdna_pos_1based: usize,
    strand: char,
}

fn load_vcf_variants(
    vcf_path: &str,
    transcript_map_path: &str,
    requested_sample: Option<&str>,
    transcripts: &[TranscriptRecord],
    warnings: &mut Vec<String>,
) -> Result<Vec<VariantRecord>, String> {
    let projections = load_transcript_coordinate_map(transcript_map_path, transcripts)?;
    let mut reader = open_maybe_gz_reader(vcf_path)?;
    let mut line = String::new();
    let mut sample_names = Vec::<String>::new();
    let mut sample_index: Option<usize> = None;
    let mut variants = Vec::new();
    let mut skipped_filter_rows = 0usize;
    let mut skipped_unprojected_rows = 0usize;
    loop {
        line.clear();
        let bytes = reader
            .read_line(&mut line)
            .map_err(|e| format!("Could not read VCF '{vcf_path}': {e}"))?;
        if bytes == 0 {
            break;
        }
        let trimmed = line.trim_end_matches(['\r', '\n']);
        if trimmed.starts_with("##") || trimmed.is_empty() {
            continue;
        }
        if trimmed.starts_with("#CHROM") {
            let fields = trimmed.split('\t').collect::<Vec<_>>();
            sample_names = fields
                .iter()
                .skip(9)
                .map(|value| (*value).to_string())
                .collect();
            sample_index = Some(resolve_vcf_sample_index(
                &sample_names,
                requested_sample,
                vcf_path,
            )?);
            continue;
        }
        if trimmed.starts_with('#') {
            continue;
        }
        let selected_sample = sample_index.ok_or_else(|| {
            format!("VCF '{vcf_path}' is missing a #CHROM header with sample columns")
        })?;
        let fields = trimmed.split('\t').collect::<Vec<_>>();
        if fields.len() < 10 + selected_sample {
            return Err(format!(
                "VCF '{vcf_path}' row has {} columns but sample '{}' requires column {}",
                fields.len(),
                sample_names[selected_sample],
                10 + selected_sample
            ));
        }
        if fields[6] != "PASS" {
            skipped_filter_rows = skipped_filter_rows.saturating_add(1);
            continue;
        }
        let genomic_pos_1based = fields[1].parse::<usize>().map_err(|e| {
            format!(
                "Invalid VCF position '{}' in '{vcf_path}' for {}: {e}",
                fields[1], fields[0]
            )
        })?;
        let Some(mapped) = projections.get(&(fields[0].to_string(), genomic_pos_1based)) else {
            skipped_unprojected_rows = skipped_unprojected_rows.saturating_add(1);
            continue;
        };
        let alt_values = fields[4].split(',').collect::<Vec<_>>();
        if fields[3].len() != 1 || alt_values.len() != 1 || alt_values[0].len() != 1 {
            return Err(format!(
                "VCF '{vcf_path}' {}:{} is not a biallelic SNV: REF='{}' ALT='{}'",
                fields[0], fields[1], fields[3], fields[4]
            ));
        }
        let format_keys = fields[8].split(':').collect::<Vec<_>>();
        let sample_values = fields[9 + selected_sample].split(':').collect::<Vec<_>>();
        let genotype = vcf_sample_value(&format_keys, &sample_values, "GT").ok_or_else(|| {
            format!(
                "VCF '{vcf_path}' {}:{} sample '{}' is missing GT",
                fields[0], fields[1], sample_names[selected_sample]
            )
        })?;
        let assignment = parse_genotype_assignment(genotype).map_err(|detail| {
            format!(
                "Invalid VCF genotype '{genotype}' at {}:{} for sample '{}': {detail}",
                fields[0], fields[1], sample_names[selected_sample]
            )
        })?;
        let phase_set = vcf_sample_value(&format_keys, &sample_values, "PS")
            .filter(|value| !value.is_empty() && *value != ".")
            .map(str::to_string);
        for projection in mapped {
            let mut ref_allele = fields[3].as_bytes()[0].to_ascii_uppercase();
            let mut alt_allele = alt_values[0].as_bytes()[0].to_ascii_uppercase();
            if projection.strand == '-' {
                ref_allele = complement_base(ref_allele)?;
                alt_allele = complement_base(alt_allele)?;
            }
            let base_id = if fields[2].is_empty() || fields[2] == "." {
                format!("{}:{}:{}>{}", fields[0], fields[1], fields[3], fields[4])
            } else {
                fields[2].to_string()
            };
            variants.push(VariantRecord {
                variant_id: format!("{base_id}@{}", projection.transcript_id),
                transcript_id: projection.transcript_id.clone(),
                cdna_pos_1based: projection.cdna_pos_1based,
                ref_allele: (ref_allele as char).to_string(),
                alt_allele: (alt_allele as char).to_string(),
                genotype: genotype.to_string(),
                phase_set: phase_set.clone(),
                hap1_alt: assignment.0,
                hap2_alt: assignment.1,
                phased: assignment.2,
            });
        }
    }
    if skipped_filter_rows > 0 {
        warnings.push(format!(
            "Skipped {skipped_filter_rows} non-PASS VCF row(s); only explicitly PASS rows enter allele hashes."
        ));
    }
    if skipped_unprojected_rows > 0 {
        warnings.push(format!(
            "Skipped {skipped_unprojected_rows} VCF row(s) without a selected-transcript coordinate mapping."
        ));
    }
    if variants.is_empty() {
        return Err(format!(
            "VCF '{vcf_path}' did not yield PASS biallelic SNVs for the selected transcripts through '{transcript_map_path}'"
        ));
    }
    Ok(variants)
}

fn resolve_vcf_sample_index(
    sample_names: &[String],
    requested_sample: Option<&str>,
    vcf_path: &str,
) -> Result<usize, String> {
    if let Some(requested) = requested_sample {
        return sample_names
            .iter()
            .position(|sample| sample == requested)
            .ok_or_else(|| {
                format!(
                    "VCF '{vcf_path}' does not contain requested sample '{requested}' (available: {})",
                    sample_names.join(",")
                )
            });
    }
    match sample_names.len() {
        1 => Ok(0),
        0 => Err(format!(
            "VCF '{vcf_path}' has no sample genotype column; a sample GT is required"
        )),
        _ => Err(format!(
            "VCF '{vcf_path}' contains multiple samples ({}); select one with --vcf-sample",
            sample_names.join(",")
        )),
    }
}

fn vcf_sample_value<'a>(
    format_keys: &[&str],
    sample_values: &'a [&str],
    key: &str,
) -> Option<&'a str> {
    format_keys
        .iter()
        .position(|candidate| *candidate == key)
        .and_then(|index| sample_values.get(index).copied())
}

fn load_transcript_coordinate_map(
    path: &str,
    transcripts: &[TranscriptRecord],
) -> Result<BTreeMap<(String, usize), Vec<TranscriptCoordinateProjection>>, String> {
    let selected_transcripts = transcripts
        .iter()
        .map(|record| record.id.as_str())
        .collect::<HashSet<_>>();
    let mut reader = open_maybe_gz_reader(path)?;
    let mut line = String::new();
    let mut header = Vec::<String>::new();
    let mut projections = BTreeMap::<(String, usize), Vec<TranscriptCoordinateProjection>>::new();
    loop {
        line.clear();
        let bytes = reader
            .read_line(&mut line)
            .map_err(|e| format!("Could not read transcript map '{path}': {e}"))?;
        if bytes == 0 {
            break;
        }
        let trimmed = line.trim();
        if trimmed.is_empty() || trimmed.starts_with('#') {
            continue;
        }
        if header.is_empty() {
            header = trimmed
                .split('\t')
                .map(|value| value.trim().to_ascii_lowercase())
                .collect();
            continue;
        }
        let values = trimmed.split('\t').collect::<Vec<_>>();
        let chrom = mapped_column(&header, &values, &["chrom", "chromosome"], path)?;
        let genomic_pos = mapped_column(
            &header,
            &values,
            &["genomic_pos_1based", "pos_1based", "pos"],
            path,
        )?
        .parse::<usize>()
        .map_err(|e| format!("Invalid genomic position in transcript map '{path}': {e}"))?;
        let transcript_id = mapped_column(&header, &values, &["transcript_id"], path)?;
        if !selected_transcripts.contains(transcript_id) {
            continue;
        }
        let cdna_pos = mapped_column(
            &header,
            &values,
            &["cdna_pos_1based", "transcript_pos_1based"],
            path,
        )?
        .parse::<usize>()
        .map_err(|e| format!("Invalid cDNA position in transcript map '{path}': {e}"))?;
        let strand = mapped_optional_column(&header, &values, &["strand"])
            .unwrap_or("+")
            .chars()
            .next()
            .unwrap_or('+');
        if !matches!(strand, '+' | '-') {
            return Err(format!(
                "Invalid strand '{strand}' in transcript map '{path}': expected '+' or '-'"
            ));
        }
        projections
            .entry((chrom.to_string(), genomic_pos))
            .or_default()
            .push(TranscriptCoordinateProjection {
                transcript_id: transcript_id.to_string(),
                cdna_pos_1based: cdna_pos,
                strand,
            });
    }
    if projections.is_empty() {
        return Err(format!(
            "Transcript map '{path}' did not contain rows for the selected transcripts"
        ));
    }
    Ok(projections)
}

fn mapped_column<'a>(
    header: &[String],
    values: &'a [&str],
    aliases: &[&str],
    path: &str,
) -> Result<&'a str, String> {
    mapped_optional_column(header, values, aliases).ok_or_else(|| {
        format!(
            "Transcript map '{path}' is missing required column {}",
            aliases.join("/")
        )
    })
}

fn mapped_optional_column<'a>(
    header: &[String],
    values: &'a [&str],
    aliases: &[&str],
) -> Option<&'a str> {
    header
        .iter()
        .position(|column| aliases.iter().any(|alias| column == alias))
        .and_then(|index| values.get(index).copied())
        .map(str::trim)
        .filter(|value| !value.is_empty())
}

fn complement_base(base: u8) -> Result<u8, String> {
    match base.to_ascii_uppercase() {
        b'A' => Ok(b'T'),
        b'C' => Ok(b'G'),
        b'G' => Ok(b'C'),
        b'T' | b'U' => Ok(b'A'),
        other => Err(format!(
            "Cannot complement non-ACGT VCF allele '{}'",
            other as char
        )),
    }
}

fn required_table_value<'a>(
    row: &'a HashMap<&str, &str>,
    key: &str,
    path: &str,
    line_number: usize,
) -> Result<&'a str, String> {
    row.get(key)
        .copied()
        .filter(|value| !value.trim().is_empty())
        .ok_or_else(|| format!("Variant table '{path}' line {line_number} is missing '{key}'"))
}

fn parse_genotype_assignment(genotype: &str) -> Result<(bool, bool, bool), String> {
    let trimmed = genotype.trim();
    if let Some((left, right)) = trimmed.split_once('|') {
        let left = parse_biallelic_genotype_allele(left)?;
        let right = parse_biallelic_genotype_allele(right)?;
        return Ok((left, right, true));
    }
    if let Some((left, right)) = trimmed.split_once('/') {
        let has_alt =
            parse_biallelic_genotype_allele(left)? || parse_biallelic_genotype_allele(right)?;
        return Ok((has_alt, false, false));
    }
    Ok((parse_biallelic_genotype_allele(trimmed)?, false, false))
}

fn parse_biallelic_genotype_allele(raw: &str) -> Result<bool, String> {
    match raw.trim() {
        "0" => Ok(false),
        "1" => Ok(true),
        other => Err(format!(
            "expected a biallelic 0/1 genotype allele, found '{other}'"
        )),
    }
}

fn apply_variants_to_haplotypes(
    transcripts: &mut [TranscriptRecord],
    variants: &[VariantRecord],
) -> Result<(), String> {
    let mut by_id = transcripts
        .iter_mut()
        .map(|record| (record.id.clone(), record))
        .collect::<HashMap<_, _>>();
    for variant in variants {
        let Some(record) = by_id.get_mut(&variant.transcript_id) else {
            continue;
        };
        let idx = variant
            .cdna_pos_1based
            .checked_sub(1)
            .ok_or_else(|| format!("Variant '{}' has position 0", variant.variant_id))?;
        if idx >= record.reference.len() {
            return Err(format!(
                "Variant '{}' position {} is outside transcript '{}' length {}",
                variant.variant_id,
                variant.cdna_pos_1based,
                record.id,
                record.reference.len()
            ));
        }
        let observed = record.reference.as_bytes()[idx].to_ascii_uppercase();
        let expected = variant.ref_allele.as_bytes()[0].to_ascii_uppercase();
        if observed != expected {
            return Err(format!(
                "Variant '{}' expected ref '{}' at {}:{}, but transcript FASTA has '{}'; refusing to build allele hashes from mismatched coordinates",
                variant.variant_id,
                variant.ref_allele,
                record.id,
                variant.cdna_pos_1based,
                observed as char
            ));
        }
        let alt = variant.alt_allele.as_bytes()[0].to_ascii_uppercase();
        if variant.hap1_alt {
            record.hap1[idx] = alt;
        }
        if variant.hap2_alt {
            record.hap2[idx] = alt;
        }
    }
    Ok(())
}

fn phase_block_id(variant: &VariantRecord) -> String {
    if variant.phased {
        if let Some(phase_set) = variant.phase_set.as_deref() {
            return format!("{}:phase_set:{phase_set}", variant.transcript_id);
        }
        return format!(
            "{}:phased_variant:{}",
            variant.transcript_id, variant.variant_id
        );
    }
    format!(
        "{}:unphased_variant:{}",
        variant.transcript_id, variant.variant_id
    )
}

fn build_phase_block_indexes(
    transcripts: &[TranscriptRecord],
    variants: &[VariantRecord],
    k: usize,
) -> Result<Vec<PhaseBlockIndex>, String> {
    let by_transcript = transcripts
        .iter()
        .map(|record| (record.id.as_str(), record))
        .collect::<HashMap<_, _>>();
    let mut grouped = BTreeMap::<String, Vec<&VariantRecord>>::new();
    for variant in variants {
        grouped
            .entry(phase_block_id(variant))
            .or_default()
            .push(variant);
    }
    grouped
        .into_iter()
        .map(|(block_id, block_variants)| {
            let first = block_variants[0];
            let transcript = by_transcript
                .get(first.transcript_id.as_str())
                .copied()
                .ok_or_else(|| {
                    format!(
                        "Phase block '{block_id}' references missing transcript '{}'",
                        first.transcript_id
                    )
                })?;
            let status = if first.phased {
                AllelePhaseBlockStatus::Phased
            } else {
                AllelePhaseBlockStatus::Unphased
            };
            let mut hap1 = transcript.reference.as_bytes().to_vec();
            let mut hap2 = transcript.reference.as_bytes().to_vec();
            let mut alternate = transcript.reference.as_bytes().to_vec();
            for variant in &block_variants {
                let idx = variant.cdna_pos_1based.checked_sub(1).ok_or_else(|| {
                    format!("Variant '{}' has position 0", variant.variant_id)
                })?;
                if idx >= transcript.reference.len() {
                    return Err(format!(
                        "Variant '{}' position {} is outside transcript '{}' length {}",
                        variant.variant_id,
                        variant.cdna_pos_1based,
                        transcript.id,
                        transcript.reference.len()
                    ));
                }
                let observed = transcript.reference.as_bytes()[idx].to_ascii_uppercase();
                let expected = variant.ref_allele.as_bytes()[0].to_ascii_uppercase();
                if observed != expected {
                    return Err(format!(
                        "Variant '{}' expected ref '{}' at {}:{}, but transcript FASTA has '{}'; refusing to build allele hashes from mismatched coordinates",
                        variant.variant_id,
                        variant.ref_allele,
                        transcript.id,
                        variant.cdna_pos_1based,
                        observed as char
                    ));
                }
                let alt = variant.alt_allele.as_bytes()[0].to_ascii_uppercase();
                match status {
                    AllelePhaseBlockStatus::Phased => {
                        if variant.hap1_alt {
                            hap1[idx] = alt;
                        }
                        if variant.hap2_alt {
                            hap2[idx] = alt;
                        }
                    }
                    AllelePhaseBlockStatus::Unphased => {
                        if variant.hap1_alt || variant.hap2_alt {
                            alternate[idx] = alt;
                        }
                    }
                }
            }
            let reference = canonical_kmer_set(transcript.reference.as_bytes(), k);
            let hap1_set = canonical_kmer_set(&hap1, k);
            let hap2_set = canonical_kmer_set(&hap2, k);
            let alternate_set = canonical_kmer_set(&alternate, k);
            let (reference_only, hap1_unique, hap2_unique, alternate_unique) = match status {
                AllelePhaseBlockStatus::Phased => (
                    set_difference_from_two(&reference, &hap1_set, &hap2_set),
                    set_difference_from_two(&hap1_set, &hap2_set, &reference),
                    set_difference_from_two(&hap2_set, &hap1_set, &reference),
                    HashSet::new(),
                ),
                AllelePhaseBlockStatus::Unphased => (
                    reference.difference(&alternate_set).copied().collect(),
                    HashSet::new(),
                    HashSet::new(),
                    alternate_set.difference(&reference).copied().collect(),
                ),
            };
            let mut all = HashSet::new();
            let mut variant_kmers = BTreeMap::new();
            for variant in &block_variants {
                let pos = variant.cdna_pos_1based - 1;
                let mut local = overlapping_kmers(transcript.reference.as_bytes(), pos, k);
                local.extend(overlapping_kmers(&hap1, pos, k));
                local.extend(overlapping_kmers(&hap2, pos, k));
                local.extend(overlapping_kmers(&alternate, pos, k));
                all.extend(local.iter().copied());
                variant_kmers.insert(variant.variant_id.clone(), local);
            }
            Ok(PhaseBlockIndex {
                block_id,
                transcript_id: first.transcript_id.clone(),
                phase_set: first.phase_set.clone(),
                status,
                variant_ids: block_variants
                    .iter()
                    .map(|variant| variant.variant_id.clone())
                    .collect(),
                variant_kmers,
                reference_only,
                hap1_unique,
                hap2_unique,
                alternate_unique,
                all,
            })
        })
        .collect()
}

fn canonical_kmer_set(sequence: &[u8], k: usize) -> HashSet<u64> {
    let mut kmers = HashSet::new();
    canonical_kmers_for_each(sequence, k, |kmer| {
        kmers.insert(kmer);
    });
    kmers
}

fn overlapping_kmers(sequence: &[u8], pos_0based: usize, k: usize) -> HashSet<u64> {
    if k == 0 || sequence.len() < k {
        return HashSet::new();
    }
    let min_start = pos_0based.saturating_add(1).saturating_sub(k);
    let max_start = pos_0based.min(sequence.len().saturating_sub(k));
    let mut kmers = HashSet::new();
    for start in min_start..=max_start {
        canonical_kmers_for_each(&sequence[start..start + k], k, |kmer| {
            kmers.insert(kmer);
        });
    }
    kmers
}

fn set_difference_from_two(
    source: &HashSet<u64>,
    other_a: &HashSet<u64>,
    other_b: &HashSet<u64>,
) -> HashSet<u64> {
    source
        .iter()
        .filter(|kmer| !other_a.contains(kmer) && !other_b.contains(kmer))
        .copied()
        .collect()
}

fn build_kmer_index(
    transcripts: &[TranscriptRecord],
    phase_blocks: &[PhaseBlockIndex],
    k: usize,
) -> KmerIndex {
    let mut index = KmerIndex::default();
    for transcript in transcripts {
        canonical_kmers_for_each(transcript.reference.as_bytes(), k, |kmer| {
            index.reference.insert(kmer);
        });
        canonical_kmers_for_each(&transcript.hap1, k, |kmer| {
            index.hap1.insert(kmer);
        });
        canonical_kmers_for_each(&transcript.hap2, k, |kmer| {
            index.hap2.insert(kmer);
        });
    }
    index.reference_only = index
        .reference
        .difference(&index.hap1)
        .copied()
        .collect::<HashSet<_>>()
        .difference(&index.hap2)
        .copied()
        .collect();
    index.hap1_unique = index
        .hap1
        .difference(&index.hap2)
        .copied()
        .collect::<HashSet<_>>()
        .difference(&index.reference)
        .copied()
        .collect();
    index.hap2_unique = index
        .hap2
        .difference(&index.hap1)
        .copied()
        .collect::<HashSet<_>>()
        .difference(&index.reference)
        .copied()
        .collect();
    for block in phase_blocks {
        if block.status == AllelePhaseBlockStatus::Unphased {
            index
                .alternate_unique
                .extend(block.alternate_unique.iter().copied());
            index.alternate.extend(block.all.iter().copied());
        }
    }
    index
}

fn build_transcript_kmer_indexes(
    transcripts: &[TranscriptRecord],
    phase_blocks: &[PhaseBlockIndex],
    k: usize,
) -> Vec<TranscriptKmerIndex> {
    transcripts
        .iter()
        .map(|transcript| {
            let mut all = HashSet::new();
            canonical_kmers_for_each(transcript.reference.as_bytes(), k, |kmer| {
                all.insert(kmer);
            });
            canonical_kmers_for_each(&transcript.hap1, k, |kmer| {
                all.insert(kmer);
            });
            canonical_kmers_for_each(&transcript.hap2, k, |kmer| {
                all.insert(kmer);
            });
            for block in phase_blocks
                .iter()
                .filter(|block| block.transcript_id == transcript.id)
            {
                all.extend(block.all.iter().copied());
            }
            TranscriptKmerIndex {
                transcript_id: transcript.id.clone(),
                all,
            }
        })
        .collect()
}

struct ReadClassificationInput<'a> {
    read_id: String,
    source_file: &'a str,
    source_origins: Vec<AlleleReadSourceOrigin>,
    evidence_unit: AlleleEvidenceUnit,
    input_read_count: u8,
    sequence: &'a str,
    source_length: usize,
}

struct ReadClassificationContext<'a> {
    k: usize,
    min_unique_hits: u64,
    index: &'a KmerIndex,
    transcript_indexes: &'a [TranscriptKmerIndex],
    phase_blocks: &'a [PhaseBlockIndex],
}

fn classify_read(
    input: ReadClassificationInput<'_>,
    context: ReadClassificationContext<'_>,
) -> AlleleReadCall {
    let normalized = normalize_sequence(input.sequence);
    let mut read_kmers = HashSet::<u64>::new();
    canonical_kmers_for_each(normalized.as_bytes(), context.k, |kmer| {
        read_kmers.insert(kmer);
    });
    let read_kmer_count = read_kmers.len() as u64;
    let reference_hits = intersection_count(&read_kmers, &context.index.reference);
    let hap1_hits = intersection_count(&read_kmers, &context.index.hap1);
    let hap2_hits = intersection_count(&read_kmers, &context.index.hap2);
    let hap1_unique_hits = intersection_count(&read_kmers, &context.index.hap1_unique);
    let hap2_unique_hits = intersection_count(&read_kmers, &context.index.hap2_unique);
    let alternate_unique_hits = intersection_count(&read_kmers, &context.index.alternate_unique);
    let reference_only_hits = intersection_count(&read_kmers, &context.index.reference_only);
    let phase_block_calls = context
        .phase_blocks
        .iter()
        .map(|block| classify_phase_block(&read_kmers, block, context.min_unique_hits))
        .collect::<Vec<_>>();
    let informative_block_calls = phase_block_calls
        .iter()
        .filter(|call| call.classification != AllelePhaseBlockClassification::Uninformative)
        .collect::<Vec<_>>();
    let classification = if read_kmer_count == 0 {
        AlleleReadClassification::Uninformative
    } else if context.phase_blocks.len() == 1 {
        match phase_block_calls[0].classification {
            AllelePhaseBlockClassification::Haplotype1 => AlleleReadClassification::Hap1,
            AllelePhaseBlockClassification::Haplotype2 => AlleleReadClassification::Hap2,
            AllelePhaseBlockClassification::Alternate => AlleleReadClassification::Alternate,
            AllelePhaseBlockClassification::Reference => AlleleReadClassification::ReferenceOnly,
            AllelePhaseBlockClassification::Ambiguous => AlleleReadClassification::Ambiguous,
            AllelePhaseBlockClassification::Uninformative => {
                if reference_hits >= context.min_unique_hits {
                    AlleleReadClassification::Ambiguous
                } else {
                    AlleleReadClassification::OffTarget
                }
            }
        }
    } else if informative_block_calls.is_empty() {
        if reference_hits >= context.min_unique_hits {
            AlleleReadClassification::Ambiguous
        } else {
            AlleleReadClassification::OffTarget
        }
    } else if informative_block_calls
        .iter()
        .all(|call| call.classification == AllelePhaseBlockClassification::Reference)
    {
        AlleleReadClassification::ReferenceOnly
    } else if informative_block_calls
        .iter()
        .all(|call| call.classification == AllelePhaseBlockClassification::Alternate)
    {
        AlleleReadClassification::Alternate
    } else {
        AlleleReadClassification::Ambiguous
    };
    let supporting_variants = informative_block_calls
        .iter()
        .flat_map(|call| call.variant_ids.iter().cloned())
        .collect::<BTreeSet<_>>()
        .into_iter()
        .collect();
    let matched_transcripts = matched_transcripts_for_read(&read_kmers, context.transcript_indexes);
    AlleleReadCall {
        read_id: input.read_id,
        source_file: input.source_file.to_string(),
        source_origins: input.source_origins,
        evidence_unit: input.evidence_unit,
        input_read_count: input.input_read_count,
        classification,
        read_length: input.source_length,
        read_kmer_count,
        reference_hits,
        hap1_hits,
        hap2_hits,
        hap1_unique_hits,
        hap2_unique_hits,
        alternate_unique_hits,
        reference_only_hits,
        matched_transcripts,
        supporting_variants,
        phase_block_calls,
    }
}

fn classify_phase_block(
    read_kmers: &HashSet<u64>,
    block: &PhaseBlockIndex,
    min_unique_hits: u64,
) -> AllelePhaseBlockReadCall {
    let reference_unique_hits = intersection_count(read_kmers, &block.reference_only);
    let hap1_unique_hits = intersection_count(read_kmers, &block.hap1_unique);
    let hap2_unique_hits = intersection_count(read_kmers, &block.hap2_unique);
    let alternate_unique_hits = intersection_count(read_kmers, &block.alternate_unique);
    let all_hits = intersection_count(read_kmers, &block.all);
    let passing = [
        reference_unique_hits >= min_unique_hits,
        hap1_unique_hits >= min_unique_hits,
        hap2_unique_hits >= min_unique_hits,
        alternate_unique_hits >= min_unique_hits,
    ];
    let passing_count = passing.iter().filter(|passes| **passes).count();
    let classification = if passing_count > 1 {
        AllelePhaseBlockClassification::Ambiguous
    } else if hap1_unique_hits >= min_unique_hits {
        AllelePhaseBlockClassification::Haplotype1
    } else if hap2_unique_hits >= min_unique_hits {
        AllelePhaseBlockClassification::Haplotype2
    } else if alternate_unique_hits >= min_unique_hits {
        AllelePhaseBlockClassification::Alternate
    } else if reference_unique_hits >= min_unique_hits {
        AllelePhaseBlockClassification::Reference
    } else if all_hits >= min_unique_hits {
        AllelePhaseBlockClassification::Ambiguous
    } else {
        AllelePhaseBlockClassification::Uninformative
    };
    let variant_ids = block
        .variant_kmers
        .iter()
        .filter(|(_, kmers)| !read_kmers.is_disjoint(kmers))
        .map(|(variant_id, _)| variant_id.clone())
        .collect();
    AllelePhaseBlockReadCall {
        block_id: block.block_id.clone(),
        classification,
        reference_unique_hits,
        hap1_unique_hits,
        hap2_unique_hits,
        alternate_unique_hits,
        variant_ids,
    }
}

fn matched_transcripts_for_read(
    read_kmers: &HashSet<u64>,
    transcript_indexes: &[TranscriptKmerIndex],
) -> Vec<String> {
    if read_kmers.is_empty() {
        return Vec::new();
    }
    transcript_indexes
        .iter()
        .filter(|index| !index.all.is_disjoint(read_kmers))
        .map(|index| index.transcript_id.clone())
        .collect::<BTreeSet<_>>()
        .into_iter()
        .collect()
}

fn intersection_count(read_kmers: &HashSet<u64>, target: &HashSet<u64>) -> u64 {
    read_kmers
        .iter()
        .filter(|kmer| target.contains(kmer))
        .count() as u64
}

fn summarize_transcripts(
    transcripts: &[TranscriptRecord],
    variants: &[VariantRecord],
    accumulated: &BTreeMap<String, CountSummary>,
) -> Vec<AlleleTranscriptSummary> {
    let mut summaries = transcripts
        .iter()
        .map(|transcript| {
            (
                transcript.id.as_str(),
                accumulated.get(&transcript.id).cloned().unwrap_or_default(),
            )
        })
        .collect::<BTreeMap<_, _>>();
    let variant_counts =
        variants
            .iter()
            .fold(BTreeMap::<&str, usize>::new(), |mut acc, variant| {
                *acc.entry(variant.transcript_id.as_str()).or_insert(0) += 1;
                acc
            });
    transcripts
        .iter()
        .map(|transcript| {
            let counts = summaries.remove(transcript.id.as_str()).unwrap_or_default();
            transcript_summary(
                transcript,
                counts,
                variant_counts
                    .get(transcript.id.as_str())
                    .copied()
                    .unwrap_or_default(),
            )
        })
        .collect()
}

fn summarize_variants(
    variants: &[VariantRecord],
    accumulated: &BTreeMap<String, CountSummary>,
) -> Vec<AlleleVariantSummary> {
    let mut counts = variants
        .iter()
        .map(|variant| {
            (
                variant.variant_id.as_str(),
                accumulated
                    .get(&variant.variant_id)
                    .cloned()
                    .unwrap_or_default(),
            )
        })
        .collect::<BTreeMap<_, _>>();
    variants
        .iter()
        .map(|variant| {
            let count = counts
                .remove(variant.variant_id.as_str())
                .unwrap_or_default();
            variant_summary(variant, count)
        })
        .collect()
}

fn accumulate_read_call(
    call: &AlleleReadCall,
    transcript_counts: &mut BTreeMap<String, CountSummary>,
    variant_counts: &mut BTreeMap<String, CountSummary>,
    phase_block_counts: &mut BTreeMap<String, CountSummary>,
    classification_counts: &mut BTreeMap<String, u64>,
) {
    *classification_counts
        .entry(call.classification.as_str().to_string())
        .or_insert(0) += 1;
    for transcript_id in &call.matched_transcripts {
        add_classification_count(
            transcript_counts.entry(transcript_id.clone()).or_default(),
            call.classification,
        );
    }
    for block_call in &call.phase_block_calls {
        add_phase_block_classification_count(
            phase_block_counts
                .entry(block_call.block_id.clone())
                .or_default(),
            block_call.classification,
        );
        for variant_id in &block_call.variant_ids {
            add_phase_block_classification_count(
                variant_counts.entry(variant_id.clone()).or_default(),
                block_call.classification,
            );
        }
    }
}

fn add_phase_block_classification_count(
    count: &mut CountSummary,
    classification: AllelePhaseBlockClassification,
) {
    match classification {
        AllelePhaseBlockClassification::Haplotype1 => {
            count.informative_reads = count.informative_reads.saturating_add(1);
            count.hap1_reads = count.hap1_reads.saturating_add(1);
        }
        AllelePhaseBlockClassification::Haplotype2 => {
            count.informative_reads = count.informative_reads.saturating_add(1);
            count.hap2_reads = count.hap2_reads.saturating_add(1);
        }
        AllelePhaseBlockClassification::Alternate => {
            count.informative_reads = count.informative_reads.saturating_add(1);
            count.alternate_reads = count.alternate_reads.saturating_add(1);
        }
        AllelePhaseBlockClassification::Reference => {
            count.informative_reads = count.informative_reads.saturating_add(1);
            count.reference_only_reads = count.reference_only_reads.saturating_add(1);
        }
        AllelePhaseBlockClassification::Ambiguous => {
            count.ambiguous_reads = count.ambiguous_reads.saturating_add(1);
        }
        AllelePhaseBlockClassification::Uninformative => {
            count.uninformative_reads = count.uninformative_reads.saturating_add(1);
        }
    }
}

fn summarize_phase_blocks(
    blocks: &[PhaseBlockIndex],
    accumulated: &BTreeMap<String, CountSummary>,
) -> Vec<AllelePhaseBlockSummary> {
    blocks
        .iter()
        .map(|block| {
            let counts = accumulated
                .get(&block.block_id)
                .cloned()
                .unwrap_or_default();
            AllelePhaseBlockSummary {
                block_id: block.block_id.clone(),
                transcript_id: block.transcript_id.clone(),
                phase_set: block.phase_set.clone(),
                status: block.status,
                variant_ids: block.variant_ids.clone(),
                informative_reads: counts.informative_reads,
                hap1_reads: counts.hap1_reads,
                hap2_reads: counts.hap2_reads,
                alternate_reads: counts.alternate_reads,
                reference_reads: counts.reference_only_reads,
                ambiguous_reads: counts.ambiguous_reads,
                uninformative_reads: counts.uninformative_reads,
            }
        })
        .collect()
}

fn add_classification_count(count: &mut CountSummary, classification: AlleleReadClassification) {
    match classification {
        AlleleReadClassification::Hap1 => {
            count.informative_reads = count.informative_reads.saturating_add(1);
            count.hap1_reads = count.hap1_reads.saturating_add(1);
        }
        AlleleReadClassification::Hap2 => {
            count.informative_reads = count.informative_reads.saturating_add(1);
            count.hap2_reads = count.hap2_reads.saturating_add(1);
        }
        AlleleReadClassification::Alternate => {
            count.informative_reads = count.informative_reads.saturating_add(1);
            count.alternate_reads = count.alternate_reads.saturating_add(1);
        }
        AlleleReadClassification::ReferenceOnly => {
            count.reference_only_reads = count.reference_only_reads.saturating_add(1);
        }
        AlleleReadClassification::Ambiguous => {
            count.ambiguous_reads = count.ambiguous_reads.saturating_add(1);
        }
        AlleleReadClassification::Uninformative => {
            count.uninformative_reads = count.uninformative_reads.saturating_add(1);
        }
        AlleleReadClassification::OffTarget => {
            count.off_target_reads = count.off_target_reads.saturating_add(1);
        }
    }
}

fn transcript_summary(
    transcript: &TranscriptRecord,
    count: CountSummary,
    variant_count: usize,
) -> AlleleTranscriptSummary {
    AlleleTranscriptSummary {
        transcript_id: transcript.id.clone(),
        length_nt: transcript.reference.len(),
        variant_count,
        informative_reads: count.informative_reads,
        hap1_reads: count.hap1_reads,
        hap2_reads: count.hap2_reads,
        alternate_reads: count.alternate_reads,
        reference_only_reads: count.reference_only_reads,
        ambiguous_reads: count.ambiguous_reads,
        uninformative_reads: count.uninformative_reads,
        off_target_reads: count.off_target_reads,
        imbalance_ratio_hap1_over_informative: imbalance_ratio(count.hap1_reads, count.hap2_reads),
        coverage_blind_spots: Vec::new(),
    }
}

fn variant_summary(variant: &VariantRecord, count: CountSummary) -> AlleleVariantSummary {
    AlleleVariantSummary {
        variant_id: variant.variant_id.clone(),
        transcript_id: variant.transcript_id.clone(),
        cdna_pos_1based: variant.cdna_pos_1based,
        ref_allele: variant.ref_allele.clone(),
        alt_allele: variant.alt_allele.clone(),
        genotype: variant.genotype.clone(),
        phase_set: variant.phase_set.clone(),
        informative_reads: count.informative_reads,
        hap1_reads: count.hap1_reads,
        hap2_reads: count.hap2_reads,
        alternate_reads: count.alternate_reads,
        reference_only_reads: count.reference_only_reads,
        ambiguous_reads: count.ambiguous_reads,
        imbalance_ratio_hap1_over_informative: imbalance_ratio(count.hap1_reads, count.hap2_reads),
        coverage_blind_spots: Vec::new(),
    }
}

fn imbalance_ratio(hap1: u64, hap2: u64) -> Option<f64> {
    let total = hap1.saturating_add(hap2);
    (total > 0).then(|| hap1 as f64 / total as f64)
}

fn empty_classification_counts() -> BTreeMap<String, u64> {
    let mut counts = BTreeMap::<String, u64>::new();
    for classification in [
        AlleleReadClassification::Hap1,
        AlleleReadClassification::Hap2,
        AlleleReadClassification::Alternate,
        AlleleReadClassification::ReferenceOnly,
        AlleleReadClassification::Ambiguous,
        AlleleReadClassification::Uninformative,
        AlleleReadClassification::OffTarget,
    ] {
        counts.insert(classification.as_str().to_string(), 0);
    }
    counts
}

fn output_files(out_dir: &Path, gene: &str) -> AlleleHashScreenOutputFiles {
    let stem = sanitize_file_stem(gene);
    AlleleHashScreenOutputFiles {
        report_json: out_dir
            .join(format!("{stem}.allele_hash_screen.json"))
            .display()
            .to_string(),
        reads_tsv: out_dir
            .join(format!("{stem}.allele_hash_screen.reads.tsv"))
            .display()
            .to_string(),
        reference_fasta: out_dir
            .join(format!("{stem}.reference.transcripts.fa"))
            .display()
            .to_string(),
        hap1_fasta: out_dir
            .join(format!("{stem}.hap1.transcripts.fa"))
            .display()
            .to_string(),
        hap2_fasta: out_dir
            .join(format!("{stem}.hap2.transcripts.fa"))
            .display()
            .to_string(),
    }
}

fn sanitize_file_stem(raw: &str) -> String {
    let stem = raw
        .chars()
        .map(|ch| {
            if ch.is_ascii_alphanumeric() || matches!(ch, '-' | '_') {
                ch
            } else {
                '_'
            }
        })
        .collect::<String>();
    if stem.is_empty() {
        "allele_hash_screen".to_string()
    } else {
        stem
    }
}

fn write_haplotype_fastas(
    transcripts: &[TranscriptRecord],
    outputs: &AlleleHashScreenOutputFiles,
) -> Result<(), String> {
    write_haplotype_fasta(
        &outputs.reference_fasta,
        transcripts,
        "reference",
        |record| record.reference.as_bytes().to_vec(),
    )?;
    write_haplotype_fasta(&outputs.hap1_fasta, transcripts, "hap1", |record| {
        record.hap1.clone()
    })?;
    write_haplotype_fasta(&outputs.hap2_fasta, transcripts, "hap2", |record| {
        record.hap2.clone()
    })
}

fn write_haplotype_fasta(
    path: &str,
    transcripts: &[TranscriptRecord],
    haplotype: &str,
    sequence_for: impl Fn(&TranscriptRecord) -> Vec<u8>,
) -> Result<(), String> {
    let mut writer = BufWriter::new(
        File::create(path).map_err(|e| format!("Could not create FASTA '{path}': {e}"))?,
    );
    for record in transcripts {
        writeln!(
            writer,
            ">{}_{} {} haplotype:{}",
            record.id, haplotype, record.header, haplotype
        )
        .map_err(|e| format!("Could not write FASTA '{path}': {e}"))?;
        let sequence = sequence_for(record);
        for chunk in sequence.chunks(80) {
            writer
                .write_all(chunk)
                .and_then(|_| writer.write_all(b"\n"))
                .map_err(|e| format!("Could not write FASTA '{path}': {e}"))?;
        }
    }
    writer
        .flush()
        .map_err(|e| format!("Could not flush FASTA '{path}': {e}"))
}

fn write_read_calls_header(writer: &mut impl Write, path: &str) -> Result<(), String> {
    writeln!(
        writer,
        "read_id\tsource_file\tsource_origins\tevidence_unit\tinput_read_count\tclassification\tread_length\tread_kmer_count\treference_hits\thap1_hits\thap2_hits\thap1_unique_hits\thap2_unique_hits\talternate_unique_hits\treference_only_hits\tmatched_transcripts\tsupporting_variants\tphase_block_calls"
    )
    .map_err(|e| format!("Could not write read TSV '{path}': {e}"))
}

fn write_read_call_tsv(
    writer: &mut impl Write,
    path: &str,
    call: &AlleleReadCall,
) -> Result<(), String> {
    let block_calls = call
        .phase_block_calls
        .iter()
        .map(|block| {
            format!(
                "{}:{}:ref={}:h1={}:h2={}:alt={}",
                block.block_id,
                block.classification.as_str(),
                block.reference_unique_hits,
                block.hap1_unique_hits,
                block.hap2_unique_hits,
                block.alternate_unique_hits
            )
        })
        .collect::<Vec<_>>()
        .join(";");
    writeln!(
        writer,
        "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
        tsv_cell(&call.read_id),
        tsv_cell(&call.source_file),
        call.source_origins
            .iter()
            .map(|origin| origin.as_str())
            .collect::<Vec<_>>()
            .join(";"),
        match call.evidence_unit {
            AlleleEvidenceUnit::Read => "read",
            AlleleEvidenceUnit::Fragment => "fragment",
        },
        call.input_read_count,
        call.classification.as_str(),
        call.read_length,
        call.read_kmer_count,
        call.reference_hits,
        call.hap1_hits,
        call.hap2_hits,
        call.hap1_unique_hits,
        call.hap2_unique_hits,
        call.alternate_unique_hits,
        call.reference_only_hits,
        tsv_cell(&call.matched_transcripts.join(";")),
        tsv_cell(&call.supporting_variants.join(";")),
        tsv_cell(&block_calls)
    )
    .map_err(|e| format!("Could not write read TSV '{path}': {e}"))
}

fn write_json_file<T: Serialize>(path: &str, value: &T) -> Result<(), String> {
    let file = File::create(path).map_err(|e| format!("Could not create JSON '{path}': {e}"))?;
    serde_json::to_writer_pretty(BufWriter::new(file), value)
        .map_err(|e| format!("Could not write JSON '{path}': {e}"))
}

fn tsv_cell(raw: &str) -> String {
    raw.replace('\\', "\\\\")
        .replace('\t', "\\t")
        .replace('\r', "\\r")
        .replace('\n', "\\n")
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::path::PathBuf;
    use tempfile::tempdir;

    fn fixture_path(name: &str) -> String {
        PathBuf::from("test_files/fixtures/allele_hash_screen")
            .join(name)
            .display()
            .to_string()
    }

    fn provenance_count(
        report: &AlleleHashScreenReport,
        origin: AlleleReadSourceOrigin,
    ) -> (u64, u64) {
        report
            .source_provenance
            .iter()
            .find(|row| row.origin == origin)
            .map(|row| (row.evidence_observation_count, row.input_read_count))
            .unwrap_or_default()
    }

    #[test]
    fn fus_fixture_classifies_haplotype_and_ambiguous_reads() {
        let dir = tempdir().expect("temp dir");
        let report = run_allele_hash_screen(AlleleHashScreenRequest {
            gene: "FUS".to_string(),
            transcript_fasta: fixture_path("fus_transcripts.fa"),
            variant_table: Some(fixture_path("fus_variants.tsv")),
            read_files: vec![fixture_path("fus_reads.fastq")],
            read_id_allowlist: Some(fixture_path("read_allowlist.txt")),
            out_dir: dir.path().display().to_string(),
            kmer_len: 9,
            min_unique_kmer_hits: 1,
            ..AlleleHashScreenRequest::default()
        })
        .expect("allele screen should run");
        assert_eq!(report.schema, ALLELE_HASH_SCREEN_SCHEMA);
        assert_eq!(report.phase_mode, AlleleHashPhaseMode::Phased);
        assert_eq!(report.phase_blocks.len(), 1);
        assert!(report.haplotype_fastas[1].inferred);
        assert_eq!(report.transcript_count, 2);
        assert_eq!(report.variant_count, 2);
        let by_id = report
            .reads
            .iter()
            .map(|call| (call.read_id.as_str(), call))
            .collect::<BTreeMap<_, _>>();
        assert_eq!(
            by_id["fus_hap1_alt_v2"].classification,
            AlleleReadClassification::Hap1
        );
        assert_eq!(
            by_id["fus_hap2_alt_v1"].classification,
            AlleleReadClassification::Hap2
        );
        assert_eq!(by_id["fus_hap2_alt_v1"].reference_hits, 0);
        assert!(by_id["fus_hap2_alt_v1"].hap2_unique_hits > 0);
        assert_eq!(
            by_id["fus_hap1_alt_v2"].supporting_variants,
            vec!["v2".to_string()]
        );
        assert_eq!(
            by_id["fus_hap2_alt_v1"].supporting_variants,
            vec!["v1".to_string()]
        );
        assert_eq!(
            by_id["fus_shared_region"].classification,
            AlleleReadClassification::Ambiguous
        );
        assert_eq!(
            by_id["fus_off_target"].classification,
            AlleleReadClassification::OffTarget
        );
        assert_eq!(
            by_id["fus_short_uninformative"].classification,
            AlleleReadClassification::Uninformative
        );
        assert!(Path::new(&report.output_files.report_json).exists());
        assert!(Path::new(&report.output_files.reads_tsv).exists());
        assert!(Path::new(&report.output_files.hap1_fasta).exists());
        assert!(Path::new(&report.output_files.hap2_fasta).exists());
    }

    #[test]
    fn salmon_selectors_reproduce_fixture_calls_and_provenance() {
        let dir = tempdir().expect("temp dir");
        let report = run_allele_hash_screen(AlleleHashScreenRequest {
            gene: "FUS".to_string(),
            transcript_fasta: fixture_path("fus_transcripts.fa"),
            variant_table: Some(fixture_path("fus_variants.tsv")),
            read_files: vec![fixture_path("fus_reads.fastq")],
            salmon_unmapped_names: Some(fixture_path("salmon_unmapped_names.txt")),
            salmon_mappings_sam: Some(fixture_path("salmon_mappings.sam")),
            out_dir: dir.path().display().to_string(),
            kmer_len: 9,
            min_unique_kmer_hits: 1,
            ..AlleleHashScreenRequest::default()
        })
        .expect("Salmon-sourced allele screen should run");

        assert_eq!(report.read_count_total, 5);
        assert_eq!(report.read_count_selected, 5);
        assert_eq!(
            provenance_count(&report, AlleleReadSourceOrigin::ExplicitFile),
            (5, 5)
        );
        assert_eq!(
            provenance_count(&report, AlleleReadSourceOrigin::SalmonUnassigned),
            (4, 4)
        );
        assert_eq!(
            provenance_count(&report, AlleleReadSourceOrigin::SalmonTargetMapped),
            (1, 1)
        );
        let by_id = report
            .reads
            .iter()
            .map(|call| (call.read_id.as_str(), call))
            .collect::<BTreeMap<_, _>>();
        assert_eq!(
            by_id["fus_hap1_alt_v2"].classification,
            AlleleReadClassification::Hap1
        );
        assert_eq!(
            by_id["fus_hap2_alt_v1"].classification,
            AlleleReadClassification::Hap2
        );
        assert_eq!(by_id["fus_hap2_alt_v1"].reference_hits, 0);
        assert!(
            by_id["fus_hap2_alt_v1"]
                .source_origins
                .contains(&AlleleReadSourceOrigin::SalmonUnassigned)
        );
        assert_eq!(
            by_id["fus_shared_region"].classification,
            AlleleReadClassification::Ambiguous
        );
        assert_eq!(
            by_id["fus_off_target"].classification,
            AlleleReadClassification::OffTarget
        );
        assert_eq!(
            by_id["fus_short_uninformative"].classification,
            AlleleReadClassification::Uninformative
        );
    }

    #[test]
    fn unphased_variant_labels_report_without_inventing_phase() {
        let dir = tempdir().expect("temp dir");
        let variant_table = dir.path().join("unphased.tsv");
        fs::write(
            &variant_table,
            "variant_id\ttranscript_id\tcdna_pos_1based\tref\talt\tgenotype\tphase_set\nv1\tFUS_TX1\t20\tC\tT\t0/1\t.\n",
        )
        .expect("write variant table");
        let report = run_allele_hash_screen(AlleleHashScreenRequest {
            gene: "FUS".to_string(),
            transcript_fasta: fixture_path("fus_transcripts.fa"),
            variant_table: Some(variant_table.display().to_string()),
            read_files: vec![fixture_path("fus_reads.fastq")],
            out_dir: dir.path().join("out").display().to_string(),
            kmer_len: 9,
            min_unique_kmer_hits: 1,
            ..AlleleHashScreenRequest::default()
        })
        .expect("unphased screen should run");
        assert_eq!(
            report.phase_mode,
            AlleleHashPhaseMode::UnphasedAlleleLevelOnly
        );
        assert!(
            report
                .warnings
                .iter()
                .any(|warning| warning.contains("unphased"))
        );
        let alternate = report
            .reads
            .iter()
            .find(|call| call.read_id == "fus_hap2_alt_v1")
            .expect("alternate-bearing read");
        assert_eq!(
            alternate.classification,
            AlleleReadClassification::Alternate
        );
        assert!(alternate.alternate_unique_hits > 0);
        assert_eq!(report.phase_blocks.len(), 1);
        assert_eq!(report.phase_blocks[0].alternate_reads, 1);
        assert!(!report.haplotype_fastas[1].inferred);
        assert!(!report.haplotype_fastas[2].inferred);
    }

    #[test]
    fn disconnected_phase_sets_remain_separate_blocks() {
        let dir = tempdir().expect("temp dir");
        let variant_table = dir.path().join("phase_blocks.tsv");
        fs::write(
            &variant_table,
            concat!(
                "variant_id\ttranscript_id\tcdna_pos_1based\tref\talt\tgenotype\tphase_set\n",
                "v1\tFUS_TX1\t20\tC\tT\t0|1\tps_one\n",
                "v2\tFUS_TX1\t45\tA\tG\t1|0\tps_two\n"
            ),
        )
        .expect("write variant table");
        let report = run_allele_hash_screen(AlleleHashScreenRequest {
            gene: "FUS".to_string(),
            transcript_fasta: fixture_path("fus_transcripts.fa"),
            variant_table: Some(variant_table.display().to_string()),
            read_files: vec![fixture_path("fus_reads.fastq")],
            out_dir: dir.path().join("out").display().to_string(),
            kmer_len: 9,
            min_unique_kmer_hits: 1,
            ..AlleleHashScreenRequest::default()
        })
        .expect("phase-block screen should run");
        assert_eq!(report.phase_mode, AlleleHashPhaseMode::PhasedBlocks);
        assert_eq!(report.phase_blocks.len(), 2);
        assert_eq!(report.phase_blocks[0].phase_set.as_deref(), Some("ps_one"));
        assert_eq!(report.phase_blocks[1].phase_set.as_deref(), Some("ps_two"));
        assert!(
            report
                .warnings
                .iter()
                .any(|warning| warning.contains("does not fabricate a gene-wide haplotype"))
        );
        assert!(
            report.haplotype_fastas[1..]
                .iter()
                .all(|fasta| !fasta.inferred)
        );
    }

    #[test]
    fn reference_mismatch_is_rejected_before_hashing() {
        let dir = tempdir().expect("temp dir");
        let variant_table = dir.path().join("mismatch.tsv");
        fs::write(
            &variant_table,
            "variant_id\ttranscript_id\tcdna_pos_1based\tref\talt\tgenotype\tphase_set\nv1\tFUS_TX1\t20\tA\tT\t0|1\tps1\n",
        )
        .expect("write variant table");
        let error = run_allele_hash_screen(AlleleHashScreenRequest {
            gene: "FUS".to_string(),
            transcript_fasta: fixture_path("fus_transcripts.fa"),
            variant_table: Some(variant_table.display().to_string()),
            read_files: vec![fixture_path("fus_reads.fastq")],
            out_dir: dir.path().join("out").display().to_string(),
            kmer_len: 9,
            min_unique_kmer_hits: 1,
            ..AlleleHashScreenRequest::default()
        })
        .expect_err("reference mismatch must fail");
        assert!(error.contains("refusing to build allele hashes"));
    }

    #[test]
    fn inline_read_calls_are_capped_but_tsv_remains_complete() {
        let dir = tempdir().expect("temp dir");
        let report = run_allele_hash_screen(AlleleHashScreenRequest {
            gene: "FUS".to_string(),
            transcript_fasta: fixture_path("fus_transcripts.fa"),
            variant_table: Some(fixture_path("fus_variants.tsv")),
            read_files: vec![fixture_path("fus_reads.fastq")],
            out_dir: dir.path().display().to_string(),
            kmer_len: 9,
            min_unique_kmer_hits: 1,
            max_inline_read_calls: 1,
            ..AlleleHashScreenRequest::default()
        })
        .expect("bounded report should run");
        assert_eq!(report.read_count_selected, 5);
        assert_eq!(report.reads.len(), 1);
        assert!(report.read_calls_truncated);
        let tsv = fs::read_to_string(&report.output_files.reads_tsv).expect("read TSV");
        assert_eq!(tsv.lines().count(), 6);
    }

    #[test]
    fn paired_files_are_counted_as_one_fragment_observation() {
        let dir = tempdir().expect("temp dir");
        let read1 = dir.path().join("reads_1.fastq");
        let read2 = dir.path().join("reads_2.fastq");
        let mut alternate =
            b"ATGACCAAGTTGACCTGACCGTACCGATGCTAGCTTACGATCGTACGATCTAGCTAGGCTA".to_vec();
        alternate[19] = b'T';
        let informative = std::str::from_utf8(&alternate[15..24]).expect("ASCII sequence");
        fs::write(
            &read1,
            format!("@fragment_1/1\n{informative}\n+\nIIIIIIIII\n"),
        )
        .expect("write read 1");
        fs::write(&read2, "@fragment_1/2\nNNNNNNNNN\n+\nIIIIIIIII\n").expect("write read 2");
        let report = run_allele_hash_screen(AlleleHashScreenRequest {
            gene: "FUS".to_string(),
            transcript_fasta: fixture_path("fus_transcripts.fa"),
            variant_table: Some(fixture_path("fus_variants.tsv")),
            read_pairs: vec![AlleleReadPairInput {
                read1: read1.display().to_string(),
                read2: read2.display().to_string(),
            }],
            out_dir: dir.path().join("out").display().to_string(),
            kmer_len: 9,
            min_unique_kmer_hits: 1,
            ..AlleleHashScreenRequest::default()
        })
        .expect("paired screen should run");
        assert_eq!(report.read_count_total, 2);
        assert_eq!(report.read_count_selected, 2);
        assert_eq!(report.fragment_count_total, 1);
        assert_eq!(report.fragment_count_selected, 1);
        assert_eq!(report.evidence_observation_count_selected, 1);
        assert_eq!(report.reads.len(), 1);
        assert_eq!(report.reads[0].evidence_unit, AlleleEvidenceUnit::Fragment);
        assert_eq!(report.reads[0].input_read_count, 2);
        assert_eq!(
            report.reads[0].classification,
            AlleleReadClassification::Hap2
        );
    }

    #[test]
    fn local_pass_vcf_projects_into_transcript_allele_hashes() {
        let dir = tempdir().expect("temp dir");
        let vcf = dir.path().join("sample.vcf");
        let transcript_map = dir.path().join("transcript_map.tsv");
        fs::write(
            &vcf,
            concat!(
                "##fileformat=VCFv4.2\n",
                "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tALS_SAMPLE\n",
                "chr1\t100\trs_demo\tC\tT\t60\tPASS\t.\tGT:PS\t0|1:ps_demo\n"
            ),
        )
        .expect("write VCF");
        fs::write(
            &transcript_map,
            concat!(
                "chrom\tgenomic_pos_1based\ttranscript_id\tcdna_pos_1based\tstrand\n",
                "chr1\t100\tFUS_TX1\t20\t+\n"
            ),
        )
        .expect("write transcript map");
        let report = run_allele_hash_screen(AlleleHashScreenRequest {
            gene: "FUS".to_string(),
            transcript_fasta: fixture_path("fus_transcripts.fa"),
            vcf: Some(vcf.display().to_string()),
            transcript_map: Some(transcript_map.display().to_string()),
            vcf_sample: Some("ALS_SAMPLE".to_string()),
            read_files: vec![fixture_path("fus_reads.fastq")],
            out_dir: dir.path().join("out").display().to_string(),
            kmer_len: 9,
            min_unique_kmer_hits: 1,
            ..AlleleHashScreenRequest::default()
        })
        .expect("VCF-backed screen should run");
        assert_eq!(report.variant_count, 1);
        assert_eq!(report.phase_mode, AlleleHashPhaseMode::Phased);
        assert_eq!(report.phase_blocks[0].phase_set.as_deref(), Some("ps_demo"));
        assert_eq!(report.variant_summaries[0].variant_id, "rs_demo@FUS_TX1");
        assert!(
            report
                .reads
                .iter()
                .any(|call| call.classification == AlleleReadClassification::Hap2)
        );
    }

    #[test]
    fn old_request_payload_defaults_inline_call_cap() {
        let request: AlleleHashScreenRequest = serde_json::from_value(serde_json::json!({
            "gene": "FUS",
            "transcript_fasta": "tx.fa",
            "variant_table": "variants.tsv",
            "read_files": ["reads.fq"],
            "out_dir": "out",
            "kmer_len": 21,
            "min_unique_kmer_hits": 1
        }))
        .expect("old request payload");
        assert_eq!(request.max_inline_read_calls, DEFAULT_MAX_INLINE_READ_CALLS);
        assert!(request.from_rna_report.is_none());
        assert!(request.salmon_unmapped_names.is_none());
        assert!(request.salmon_mappings_sam.is_none());
        assert!(request.rna_report_reads.is_empty());
    }

    #[test]
    fn v1_report_payload_deserializes_without_source_provenance() {
        let dir = tempdir().expect("temp dir");
        let report = run_allele_hash_screen(AlleleHashScreenRequest {
            gene: "FUS".to_string(),
            transcript_fasta: fixture_path("fus_transcripts.fa"),
            variant_table: Some(fixture_path("fus_variants.tsv")),
            read_files: vec![fixture_path("fus_reads.fastq")],
            out_dir: dir.path().display().to_string(),
            kmer_len: 9,
            min_unique_kmer_hits: 1,
            ..AlleleHashScreenRequest::default()
        })
        .expect("current report");
        let mut payload = serde_json::to_value(report).expect("serialize current report");
        let object = payload.as_object_mut().expect("report object");
        object.insert(
            "schema".to_string(),
            serde_json::Value::String(ALLELE_HASH_SCREEN_SCHEMA_V1.to_string()),
        );
        object.remove("source_provenance");
        object
            .get_mut("params")
            .and_then(serde_json::Value::as_object_mut)
            .expect("params object")
            .retain(|key, _| {
                !matches!(
                    key.as_str(),
                    "from_rna_report" | "salmon_unmapped_names" | "salmon_mappings_sam"
                )
            });
        for read in object
            .get_mut("reads")
            .and_then(serde_json::Value::as_array_mut)
            .expect("read rows")
        {
            read.as_object_mut()
                .expect("read object")
                .remove("source_origins");
        }

        let decoded: AlleleHashScreenReport =
            serde_json::from_value(payload).expect("v1 report remains readable");
        assert_eq!(decoded.schema, ALLELE_HASH_SCREEN_SCHEMA_V1);
        assert!(decoded.source_provenance.is_empty());
        assert!(
            decoded
                .reads
                .iter()
                .all(|read| read.source_origins.is_empty())
        );
    }
}
