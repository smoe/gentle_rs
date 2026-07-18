//! Standalone RNA-seq target-region rescue screen.
//!
//! This module implements a deliberately small, deterministic calibration
//! screen for reads that should be compared against requested transcript
//! targets without entering the shared RNA-read interpretation engine. It
//! builds redundant per-gene k-mer sets, streams single reads or lockstep read
//! pairs once, and writes TSV/JSON reports suitable for later folding into
//! richer RNA reports. Junction and exon-intron-boundary catalogs contribute
//! only k-mers that cross their declared boundary; retained-intron evidence
//! remains a candidate and is reported separately with or without same-fragment
//! annotated-junction support.

use flate2::read::MultiGzDecoder;
use serde::{Deserialize, Serialize};
use std::{
    collections::{BTreeMap, BTreeSet, HashMap, HashSet},
    fs::{self, File},
    io::{BufRead, BufReader, BufWriter, Write},
    path::Path,
};

const DEFAULT_KMER_LEN: usize = 25;
const DEFAULT_MIN_KMER_HITS: u64 = 3;
const DEFAULT_GENE_SYMBOL_TAGS: &[&str] = &["gene_symbol", "gene_name", "symbol", "gene"];
const SUMMARY_SCHEMA: &str = "gentle.rna_target_rescue.summary.v1";
const METADATA_SCHEMA: &str = "gentle.rna_target_rescue.run_metadata.v1";
const COMPARISON_SCHEMA: &str = "gentle.rna_target_rescue.comparison.v1";

/// Request for the standalone RNA target-rescue screen.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct TargetRescueRequest {
    pub transcript_fastas: Vec<String>,
    pub genes: Vec<String>,
    #[serde(default)]
    pub reads: Vec<String>,
    #[serde(default)]
    pub read_pairs: Vec<TargetRescueReadPairInput>,
    #[serde(default)]
    pub structural_evidence: Vec<StructuralEvidenceSource>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub read_id_allowlist: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub salmon_unmapped_names: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub salmon_mappings_sam: Option<String>,
    pub kmer_len: usize,
    pub min_kmer_hits: u64,
    pub gene_symbol_tags: Vec<String>,
    pub output_prefix: String,
}

impl Default for TargetRescueRequest {
    fn default() -> Self {
        Self {
            transcript_fastas: Vec::new(),
            genes: Vec::new(),
            reads: Vec::new(),
            read_pairs: Vec::new(),
            structural_evidence: Vec::new(),
            read_id_allowlist: None,
            salmon_unmapped_names: None,
            salmon_mappings_sam: None,
            kmer_len: DEFAULT_KMER_LEN,
            min_kmer_hits: DEFAULT_MIN_KMER_HITS,
            gene_symbol_tags: default_gene_symbol_tags(),
            output_prefix: String::new(),
        }
    }
}

/// Stable parameter snapshot echoed into JSON reports.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct TargetRescueParams {
    pub transcript_fastas: Vec<String>,
    pub reads: Vec<String>,
    #[serde(default)]
    pub read_pairs: Vec<TargetRescueReadPairInput>,
    #[serde(default)]
    pub structural_evidence: Vec<StructuralEvidenceSource>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub read_id_allowlist: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub salmon_unmapped_names: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub salmon_mappings_sam: Option<String>,
    pub kmer_len: usize,
    pub min_kmer_hits: u64,
    pub gene_symbol_tags: Vec<String>,
}

/// Per-gene transcript and hash coverage row.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct GeneManifestEntry {
    pub gene_symbol: String,
    pub requested: bool,
    pub present: bool,
    pub transcript_count: usize,
    pub distinct_kmers: usize,
    pub total_kmer_positions: u64,
    pub kmer_len: usize,
    pub transcript_ids: Vec<String>,
}

/// Per-universe read-count summary.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct UniverseSummary {
    pub universe: String,
    /// Legacy name retained for v1 compatibility; counts selected evidence units.
    pub selected_reads: u64,
    /// Legacy name retained for v1 compatibility; counts evaluated evidence units.
    pub evaluated_reads: u64,
    /// Legacy name retained for v1 compatibility; counts matching evidence units.
    pub reads_matching_any: u64,
    /// Legacy name retained for v1 compatibility; counts ambiguous evidence units.
    pub ambiguous_reads: u64,
    #[serde(default)]
    pub selected_input_reads: u64,
    #[serde(default)]
    pub selected_evidence_units: u64,
    #[serde(default)]
    pub evaluated_evidence_units: u64,
    #[serde(default)]
    pub evidence_units_matching_any: u64,
    #[serde(default)]
    pub ambiguous_evidence_units: u64,
}

/// Lockstep paired-end input for the target-rescue screen.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct TargetRescueReadPairInput {
    pub read1: String,
    pub read2: String,
}

#[derive(Debug, Clone, Copy, Default, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum TargetRescueEvidenceUnit {
    #[default]
    Read,
    Fragment,
}

impl TargetRescueEvidenceUnit {
    fn as_str(self) -> &'static str {
        match self {
            Self::Read => "read",
            Self::Fragment => "fragment",
        }
    }
}

/// Per-gene rescue hit counts inside one read universe.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct GeneHitSummary {
    pub universe: String,
    pub gene_symbol: String,
    pub present: bool,
    pub reads_hit: u64,
    pub reads_hit_unique: u64,
    pub reads_hit_ambiguous: u64,
    pub total_matching_kmers: u64,
}

/// Kind of local sequence catalog used for a conservative structural screen.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum StructuralEvidenceKind {
    Exon,
    AnnotatedJunction,
    ExonIntronBoundary,
    Intron,
    GenomicRegion,
}

impl StructuralEvidenceKind {
    fn as_str(self) -> &'static str {
        match self {
            Self::Exon => "exon",
            Self::AnnotatedJunction => "annotated_junction",
            Self::ExonIntronBoundary => "exon_intron_boundary",
            Self::Intron => "intron",
            Self::GenomicRegion => "genomic_region",
        }
    }

    fn candidate_label(self) -> &'static str {
        match self {
            Self::Exon => "known_exonic",
            Self::AnnotatedJunction => "known_junction",
            Self::ExonIntronBoundary => "retained_intron_candidate",
            Self::Intron => "intronic_candidate",
            Self::GenomicRegion => "genomic_region_candidate",
        }
    }

    fn requires_boundary_span(self) -> bool {
        matches!(self, Self::AnnotatedJunction | Self::ExonIntronBoundary)
    }
}

/// One local FASTA catalog used to screen reads for structural-region evidence.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct StructuralEvidenceSource {
    pub kind: StructuralEvidenceKind,
    pub path: String,
}

/// Aggregate matches to one gene and structural evidence class.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct StructuralEvidenceHitSummary {
    pub universe: String,
    pub gene_symbol: String,
    pub evidence_kind: StructuralEvidenceKind,
    pub candidate_label: String,
    pub record_count: usize,
    pub distinct_kmers: usize,
    pub reads_hit: u64,
    pub total_matching_kmers: u64,
    #[serde(default)]
    pub boundary_spanning_only: bool,
    #[serde(default)]
    pub rna_anchored_hits: u64,
    #[serde(default)]
    pub unanchored_hits: u64,
}

/// Paths written by a target-rescue run.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct TargetRescueOutputFiles {
    pub gene_manifest_tsv: String,
    pub gene_hits_tsv: String,
    pub read_hits_tsv: String,
    #[serde(default)]
    pub structural_hits_tsv: String,
    pub run_metadata_json: String,
    pub summary_json: String,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub comparison_json: Option<String>,
}

/// Top-level machine-readable summary printed by the CLI.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct TargetRescueSummary {
    pub schema: String,
    pub params: TargetRescueParams,
    pub requested_genes: Vec<String>,
    pub missing_genes: Vec<String>,
    pub total_input_reads: u64,
    #[serde(default)]
    pub total_evidence_units: u64,
    pub universes: Vec<UniverseSummary>,
    pub gene_hits: Vec<GeneHitSummary>,
    #[serde(default)]
    pub structural_hits: Vec<StructuralEvidenceHitSummary>,
    pub manifest: Vec<GeneManifestEntry>,
    pub output_files: TargetRescueOutputFiles,
    #[serde(default)]
    pub warnings: Vec<String>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
struct TargetRescueRunMetadata {
    schema: String,
    tool_version: String,
    params: TargetRescueParams,
    requested_genes: Vec<String>,
    missing_genes: Vec<String>,
    transcript_sources: Vec<TranscriptSourceReport>,
    read_sources: Vec<ReadSourceReport>,
    total_input_reads: u64,
    total_evidence_units: u64,
    universes: Vec<UniverseSummary>,
    gene_hits: Vec<GeneHitSummary>,
    structural_hits: Vec<StructuralEvidenceHitSummary>,
    structural_sources: Vec<StructuralEvidenceSourceReport>,
    manifest: Vec<GeneManifestEntry>,
    output_files: TargetRescueOutputFiles,
    warnings: Vec<String>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
struct TranscriptSourceReport {
    path: String,
    bytes: u64,
    record_count: u64,
    matched_record_count: u64,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
struct ReadSourceReport {
    path: String,
    bytes: u64,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
struct StructuralEvidenceSourceReport {
    kind: StructuralEvidenceKind,
    path: String,
    bytes: u64,
    record_count: u64,
    matched_record_count: u64,
    boundary_spanning_record_count: u64,
    distinct_kmers: usize,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
struct TargetRescueComparison {
    schema: String,
    universes: Vec<UniverseSummary>,
    gene_rows: Vec<ComparisonGeneRow>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
struct ComparisonGeneRow {
    gene_symbol: String,
    salmon_unmapped_reads_hit: u64,
    target_mapped_reads_hit: u64,
    delta_unmapped_minus_target_mapped: i64,
}

#[derive(Debug, Clone)]
struct GeneIndex {
    symbol: String,
    transcript_ids: BTreeSet<String>,
    kmers: HashSet<u64>,
    total_kmer_positions: u64,
}

#[derive(Debug, Clone)]
struct StructuralEvidenceIndex {
    gene_symbol: String,
    kind: StructuralEvidenceKind,
    record_ids: BTreeSet<String>,
    kmers: HashSet<u64>,
    boundary_spanning_only: bool,
}

#[derive(Debug, Clone, Default)]
struct UniverseTally {
    selected_reads: u64,
    evaluated_reads: u64,
    reads_matching_any: u64,
    ambiguous_reads: u64,
    selected_input_reads: u64,
}

#[derive(Debug, Clone, Default)]
struct GeneTally {
    reads_hit: u64,
    reads_hit_unique: u64,
    reads_hit_ambiguous: u64,
    total_matching_kmers: u64,
}

#[derive(Debug, Clone, Default)]
struct StructuralEvidenceTally {
    reads_hit: u64,
    total_matching_kmers: u64,
    rna_anchored_hits: u64,
    unanchored_hits: u64,
}

#[derive(Debug, Clone)]
struct UniverseState {
    name: String,
    selected_ids: Option<HashSet<String>>,
    tally: UniverseTally,
    gene_tallies: Vec<GeneTally>,
    structural_tallies: Vec<StructuralEvidenceTally>,
}

impl UniverseState {
    fn selects(&self, read_id: &str) -> bool {
        match &self.selected_ids {
            Some(ids) => ids.contains(read_id),
            None => true,
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum ReadFormat {
    Fasta,
    Fastq,
}

#[derive(Debug, Clone)]
pub(crate) struct SequenceRecord {
    pub(crate) id: String,
    pub(crate) sequence: String,
    pub(crate) source_length: usize,
}

pub(crate) struct ReadRecordStream {
    path: String,
    reader: Box<dyn BufRead>,
    format: ReadFormat,
    pending_fasta_header: Option<String>,
}

impl ReadRecordStream {
    pub(crate) fn open(path: &str) -> Result<Self, String> {
        let mut reader = open_maybe_gz_reader(path)?;
        let format = detect_read_format(path, reader.as_mut())?;
        Ok(Self {
            path: path.to_string(),
            reader,
            format,
            pending_fasta_header: None,
        })
    }

    pub(crate) fn next_record(&mut self) -> Result<Option<SequenceRecord>, String> {
        match self.format {
            ReadFormat::Fasta => self.next_fasta_record(),
            ReadFormat::Fastq => self.next_fastq_record(),
        }
    }

    fn next_fasta_record(&mut self) -> Result<Option<SequenceRecord>, String> {
        let mut line = String::new();
        if self.pending_fasta_header.is_none() {
            loop {
                line.clear();
                let bytes = self
                    .reader
                    .read_line(&mut line)
                    .map_err(|e| format!("Could not read FASTA reads '{}': {e}", self.path))?;
                if bytes == 0 {
                    return Ok(None);
                }
                let trimmed = line.trim_end_matches(['\r', '\n']);
                if let Some(header) = trimmed.strip_prefix('>') {
                    self.pending_fasta_header = Some(header.to_string());
                    break;
                }
                if !trimmed.trim().is_empty() {
                    return Err(format!(
                        "Invalid FASTA reads '{}': sequence data appeared before the first '>' header",
                        self.path
                    ));
                }
            }
        }
        let header = self.pending_fasta_header.take().ok_or_else(|| {
            format!(
                "Could not establish a FASTA header while reading '{}'",
                self.path
            )
        })?;
        let mut sequence = String::new();
        loop {
            line.clear();
            let bytes = self
                .reader
                .read_line(&mut line)
                .map_err(|e| format!("Could not read FASTA reads '{}': {e}", self.path))?;
            if bytes == 0 {
                break;
            }
            let trimmed = line.trim_end_matches(['\r', '\n']);
            if let Some(next_header) = trimmed.strip_prefix('>') {
                self.pending_fasta_header = Some(next_header.to_string());
                break;
            }
            sequence.push_str(trimmed.trim());
        }
        let source_length = sequence.len();
        Ok(Some(SequenceRecord {
            id: normalize_read_id(&header),
            sequence,
            source_length,
        }))
    }

    fn next_fastq_record(&mut self) -> Result<Option<SequenceRecord>, String> {
        let mut header = String::new();
        loop {
            header.clear();
            let bytes = self
                .reader
                .read_line(&mut header)
                .map_err(|e| format!("Could not read FASTQ reads '{}': {e}", self.path))?;
            if bytes == 0 {
                return Ok(None);
            }
            if !header.trim().is_empty() {
                break;
            }
        }
        if !header.starts_with('@') {
            return Err(format!(
                "Invalid FASTQ record in '{}': expected '@' header, got '{}'",
                self.path,
                header.trim_end()
            ));
        }
        let mut sequence = String::new();
        let mut plus = String::new();
        let mut quality = String::new();
        if self
            .reader
            .read_line(&mut sequence)
            .map_err(|e| format!("Could not read FASTQ sequence from '{}': {e}", self.path))?
            == 0
        {
            return Err(format!(
                "Truncated FASTQ record in '{}' after header",
                self.path
            ));
        }
        if self
            .reader
            .read_line(&mut plus)
            .map_err(|e| format!("Could not read FASTQ '+' line from '{}': {e}", self.path))?
            == 0
        {
            return Err(format!(
                "Truncated FASTQ record in '{}' before '+' line",
                self.path
            ));
        }
        if !plus.starts_with('+') {
            return Err(format!(
                "Invalid FASTQ record in '{}': expected '+' line, got '{}'",
                self.path,
                plus.trim_end()
            ));
        }
        if self
            .reader
            .read_line(&mut quality)
            .map_err(|e| format!("Could not read FASTQ quality from '{}': {e}", self.path))?
            == 0
        {
            return Err(format!(
                "Truncated FASTQ record in '{}' before quality line",
                self.path
            ));
        }
        let sequence = sequence.trim().to_string();
        let source_length = sequence.len();
        Ok(Some(SequenceRecord {
            id: normalize_read_id(&header),
            sequence,
            source_length,
        }))
    }
}

/// Run the target-region rescue screen and write all requested reports.
pub fn run_target_rescue_screen(
    request: TargetRescueRequest,
) -> Result<TargetRescueSummary, String> {
    validate_request(&request)?;
    let requested_genes = normalize_requested_genes(&request.genes)?;
    let gene_symbol_tags = if request.gene_symbol_tags.is_empty() {
        default_gene_symbol_tags()
    } else {
        request.gene_symbol_tags.clone()
    };
    let params = TargetRescueParams {
        transcript_fastas: request.transcript_fastas.clone(),
        reads: request.reads.clone(),
        read_pairs: request.read_pairs.clone(),
        structural_evidence: request.structural_evidence.clone(),
        read_id_allowlist: request.read_id_allowlist.clone(),
        salmon_unmapped_names: request.salmon_unmapped_names.clone(),
        salmon_mappings_sam: request.salmon_mappings_sam.clone(),
        kmer_len: request.kmer_len,
        min_kmer_hits: request.min_kmer_hits,
        gene_symbol_tags,
    };

    let mut gene_indices = requested_genes
        .iter()
        .map(|symbol| GeneIndex {
            symbol: symbol.clone(),
            transcript_ids: BTreeSet::new(),
            kmers: HashSet::new(),
            total_kmer_positions: 0,
        })
        .collect::<Vec<_>>();
    let gene_lookup = requested_genes
        .iter()
        .enumerate()
        .map(|(idx, symbol)| (normalize_gene_symbol(symbol), idx))
        .collect::<HashMap<_, _>>();
    let transcript_sources = build_gene_kmer_indexes(
        &request.transcript_fastas,
        request.kmer_len,
        &params.gene_symbol_tags,
        &gene_lookup,
        &mut gene_indices,
    )?;
    let (structural_indices, structural_sources) = build_structural_evidence_indexes(
        &request.structural_evidence,
        request.kmer_len,
        &params.gene_symbol_tags,
        &gene_lookup,
        &requested_genes,
    )?;
    let manifest = build_manifest(&gene_indices, request.kmer_len);
    let missing_genes = manifest
        .iter()
        .filter(|entry| !entry.present)
        .map(|entry| entry.gene_symbol.clone())
        .collect::<Vec<_>>();

    let allowlist = request
        .read_id_allowlist
        .as_deref()
        .map(read_id_set_from_path)
        .transpose()?;
    let target_transcripts = target_transcript_lookup(&gene_indices);
    let salmon_bridge =
        request.salmon_unmapped_names.is_some() || request.salmon_mappings_sam.is_some();
    let mut universes = build_universes(
        &request,
        &target_transcripts,
        allowlist.as_ref(),
        gene_indices.len(),
        structural_indices.len(),
        salmon_bridge,
    )?;
    let output_files = output_files_for_prefix(
        &request.output_prefix,
        salmon_bridge && universes.len() >= 2,
    );
    ensure_output_parent(&request.output_prefix)?;
    let mut read_hits_writer = BufWriter::new(create_file(&output_files.read_hits_tsv)?);
    write_read_hits_header(&mut read_hits_writer)?;
    let mut structural_hits_writer =
        BufWriter::new(create_file(&output_files.structural_hits_tsv)?);
    write_structural_hits_header(&mut structural_hits_writer)?;

    let mut total_input_reads = 0u64;
    let mut total_evidence_units = 0u64;
    {
        let mut process_record = |record: SequenceRecord,
                                  source_file: &str,
                                  evidence_unit: TargetRescueEvidenceUnit,
                                  input_read_count: u8|
         -> Result<(), String> {
            total_input_reads = total_input_reads.saturating_add(u64::from(input_read_count));
            total_evidence_units = total_evidence_units.saturating_add(1);
            let selected_universes = universes
                .iter()
                .enumerate()
                .filter_map(|(idx, universe)| universe.selects(&record.id).then_some(idx))
                .collect::<Vec<_>>();
            if selected_universes.is_empty() {
                return Ok(());
            }
            let (read_kmer_count, gene_match_counts) =
                count_gene_kmer_matches(&record.sequence, request.kmer_len, &gene_indices);
            let structural_match_counts = count_structural_kmer_matches(
                &record.sequence,
                request.kmer_len,
                &structural_indices,
            );
            let hit_gene_indices = gene_match_counts
                .iter()
                .enumerate()
                .filter_map(|(idx, count)| (*count >= request.min_kmer_hits).then_some(idx))
                .collect::<Vec<_>>();
            let mut hit_genes = hit_gene_indices
                .iter()
                .map(|idx| gene_indices[*idx].symbol.clone())
                .collect::<Vec<_>>();
            hit_genes.sort();
            let hit_genes_joined = hit_genes.join(";");
            let ambiguous = hit_gene_indices.len() > 1;
            let evaluated = read_kmer_count > 0;
            let passing_structural_indices = structural_match_counts
                .iter()
                .enumerate()
                .filter_map(|(idx, count)| (*count >= request.min_kmer_hits).then_some(idx))
                .collect::<Vec<_>>();
            let rna_anchor_genes = passing_structural_indices
                .iter()
                .filter_map(|idx| {
                    let evidence = &structural_indices[*idx];
                    (evidence.kind == StructuralEvidenceKind::AnnotatedJunction)
                        .then_some(evidence.gene_symbol.as_str())
                })
                .collect::<HashSet<_>>();
            for universe_idx in selected_universes {
                let universe = &mut universes[universe_idx];
                universe.tally.selected_reads = universe.tally.selected_reads.saturating_add(1);
                universe.tally.selected_input_reads = universe
                    .tally
                    .selected_input_reads
                    .saturating_add(u64::from(input_read_count));
                if evaluated {
                    universe.tally.evaluated_reads =
                        universe.tally.evaluated_reads.saturating_add(1);
                }
                if !hit_gene_indices.is_empty() {
                    universe.tally.reads_matching_any =
                        universe.tally.reads_matching_any.saturating_add(1);
                }
                if ambiguous {
                    universe.tally.ambiguous_reads =
                        universe.tally.ambiguous_reads.saturating_add(1);
                }
                for gene_idx in &hit_gene_indices {
                    let matching_kmers = gene_match_counts[*gene_idx];
                    let gene_tally = &mut universe.gene_tallies[*gene_idx];
                    gene_tally.reads_hit = gene_tally.reads_hit.saturating_add(1);
                    if ambiguous {
                        gene_tally.reads_hit_ambiguous =
                            gene_tally.reads_hit_ambiguous.saturating_add(1);
                    } else {
                        gene_tally.reads_hit_unique = gene_tally.reads_hit_unique.saturating_add(1);
                    }
                    gene_tally.total_matching_kmers = gene_tally
                        .total_matching_kmers
                        .saturating_add(matching_kmers);
                    write_read_hit_row(
                        &mut read_hits_writer,
                        &universe.name,
                        &record.id,
                        source_file,
                        evidence_unit,
                        input_read_count,
                        &gene_indices[*gene_idx].symbol,
                        matching_kmers,
                        read_kmer_count,
                        record.source_length,
                        ambiguous,
                        &hit_genes_joined,
                    )?;
                }
                for structural_idx in &passing_structural_indices {
                    let matching_kmers = structural_match_counts[*structural_idx];
                    let evidence = &structural_indices[*structural_idx];
                    let rna_anchored = match evidence.kind {
                        StructuralEvidenceKind::AnnotatedJunction => true,
                        StructuralEvidenceKind::ExonIntronBoundary
                        | StructuralEvidenceKind::Intron => {
                            rna_anchor_genes.contains(evidence.gene_symbol.as_str())
                        }
                        StructuralEvidenceKind::Exon | StructuralEvidenceKind::GenomicRegion => {
                            false
                        }
                    };
                    let tally = &mut universe.structural_tallies[*structural_idx];
                    tally.reads_hit = tally.reads_hit.saturating_add(1);
                    tally.total_matching_kmers =
                        tally.total_matching_kmers.saturating_add(matching_kmers);
                    if rna_anchored {
                        tally.rna_anchored_hits = tally.rna_anchored_hits.saturating_add(1);
                    } else {
                        tally.unanchored_hits = tally.unanchored_hits.saturating_add(1);
                    }
                    write_structural_hit_row(
                        &mut structural_hits_writer,
                        &universe.name,
                        &record.id,
                        source_file,
                        evidence_unit,
                        input_read_count,
                        evidence,
                        matching_kmers,
                        read_kmer_count,
                        record.source_length,
                        rna_anchored,
                    )?;
                }
            }
            Ok(())
        };
        for read_path in &request.reads {
            visit_read_records(read_path, |record| {
                process_record(record, read_path, TargetRescueEvidenceUnit::Read, 1)
            })?;
        }
        for read_pair in &request.read_pairs {
            let source_file = format!("{};{}", read_pair.read1, read_pair.read2);
            visit_paired_read_records(&read_pair.read1, &read_pair.read2, |record| {
                process_record(record, &source_file, TargetRescueEvidenceUnit::Fragment, 2)
            })?;
        }
    }
    read_hits_writer
        .flush()
        .map_err(|e| format!("Could not flush read hits TSV: {e}"))?;
    structural_hits_writer
        .flush()
        .map_err(|e| format!("Could not flush structural hits TSV: {e}"))?;

    let universe_summaries = build_universe_summaries(&universes);
    let gene_hits = build_gene_hit_summaries(&universes, &gene_indices);
    let structural_hits = build_structural_hit_summaries(&universes, &structural_indices);
    write_gene_manifest_tsv(&output_files.gene_manifest_tsv, &manifest)?;
    write_gene_hits_tsv(&output_files.gene_hits_tsv, &gene_hits)?;
    let read_source_paths = request.reads.iter().chain(
        request
            .read_pairs
            .iter()
            .flat_map(|pair| [&pair.read1, &pair.read2]),
    );
    let read_sources = read_source_paths
        .map(|path| {
            Ok(ReadSourceReport {
                path: path.clone(),
                bytes: file_size(path)?,
            })
        })
        .collect::<Result<Vec<_>, String>>()?;
    let warnings = structural_interpretation_warnings(&request.structural_evidence);

    let summary = TargetRescueSummary {
        schema: SUMMARY_SCHEMA.to_string(),
        params: params.clone(),
        requested_genes: requested_genes.clone(),
        missing_genes: missing_genes.clone(),
        total_input_reads,
        total_evidence_units,
        universes: universe_summaries.clone(),
        gene_hits: gene_hits.clone(),
        structural_hits: structural_hits.clone(),
        manifest: manifest.clone(),
        output_files: output_files.clone(),
        warnings: warnings.clone(),
    };
    let metadata = TargetRescueRunMetadata {
        schema: METADATA_SCHEMA.to_string(),
        tool_version: env!("CARGO_PKG_VERSION").to_string(),
        params,
        requested_genes,
        missing_genes,
        transcript_sources,
        read_sources,
        total_input_reads,
        total_evidence_units,
        universes: universe_summaries,
        gene_hits,
        structural_hits,
        structural_sources,
        manifest,
        output_files: output_files.clone(),
        warnings,
    };
    write_json_file(&output_files.run_metadata_json, &metadata)?;
    write_json_file(&output_files.summary_json, &summary)?;
    if let Some(comparison_path) = &output_files.comparison_json {
        let comparison = build_comparison_report(&summary);
        write_json_file(comparison_path, &comparison)?;
    }
    Ok(summary)
}

fn validate_request(request: &TargetRescueRequest) -> Result<(), String> {
    if request.transcript_fastas.is_empty() {
        return Err("rescue-screen requires at least one --transcript-fasta PATH".to_string());
    }
    if request.genes.is_empty() {
        return Err("rescue-screen requires at least one --genes or --gene value".to_string());
    }
    if request.reads.is_empty() && request.read_pairs.is_empty() {
        return Err(
            "rescue-screen requires at least one --reads PATH or --read-pair R1,R2".to_string(),
        );
    }
    if request.output_prefix.trim().is_empty() {
        return Err("rescue-screen requires --output-prefix PREFIX".to_string());
    }
    if !(1..=31).contains(&request.kmer_len) {
        return Err(format!(
            "Invalid --kmer-len value '{}': expected 1..=31",
            request.kmer_len
        ));
    }
    if request.min_kmer_hits == 0 {
        return Err("Invalid --min-kmer-hits value '0': expected at least 1".to_string());
    }
    Ok(())
}

fn default_gene_symbol_tags() -> Vec<String> {
    DEFAULT_GENE_SYMBOL_TAGS
        .iter()
        .map(|tag| (*tag).to_string())
        .collect()
}

fn normalize_requested_genes(raw_genes: &[String]) -> Result<Vec<String>, String> {
    let mut seen = HashSet::new();
    let mut genes = Vec::new();
    for raw in raw_genes {
        let trimmed = raw.trim();
        if trimmed.is_empty() {
            continue;
        }
        let key = normalize_gene_symbol(trimmed);
        if seen.insert(key) {
            genes.push(trimmed.to_string());
        }
    }
    if genes.is_empty() {
        return Err("rescue-screen requires at least one non-empty gene symbol".to_string());
    }
    Ok(genes)
}

fn normalize_gene_symbol(raw: &str) -> String {
    raw.trim().to_ascii_uppercase()
}

fn build_gene_kmer_indexes(
    fasta_paths: &[String],
    k: usize,
    gene_symbol_tags: &[String],
    gene_lookup: &HashMap<String, usize>,
    gene_indices: &mut [GeneIndex],
) -> Result<Vec<TranscriptSourceReport>, String> {
    let mut reports = Vec::new();
    for path in fasta_paths {
        let mut record_count = 0u64;
        let mut matched_record_count = 0u64;
        visit_fasta_records(path, |header, sequence| {
            record_count = record_count.saturating_add(1);
            let Some(symbol) = parse_gene_symbol_from_header(header, gene_symbol_tags) else {
                return Ok(());
            };
            let key = normalize_gene_symbol(&symbol);
            let Some(gene_idx) = gene_lookup.get(&key).copied() else {
                return Ok(());
            };
            matched_record_count = matched_record_count.saturating_add(1);
            let transcript_id = header
                .split_ascii_whitespace()
                .next()
                .unwrap_or("")
                .trim_start_matches('>')
                .to_string();
            if !transcript_id.is_empty() {
                gene_indices[gene_idx].transcript_ids.insert(transcript_id);
            }
            canonical_kmers_for_each(sequence.as_bytes(), k, |kmer| {
                gene_indices[gene_idx].total_kmer_positions = gene_indices[gene_idx]
                    .total_kmer_positions
                    .saturating_add(1);
                gene_indices[gene_idx].kmers.insert(kmer);
            });
            Ok(())
        })?;
        reports.push(TranscriptSourceReport {
            path: path.clone(),
            bytes: file_size(path)?,
            record_count,
            matched_record_count,
        });
    }
    Ok(reports)
}

fn build_structural_evidence_indexes(
    sources: &[StructuralEvidenceSource],
    k: usize,
    gene_symbol_tags: &[String],
    gene_lookup: &HashMap<String, usize>,
    requested_genes: &[String],
) -> Result<
    (
        Vec<StructuralEvidenceIndex>,
        Vec<StructuralEvidenceSourceReport>,
    ),
    String,
> {
    let mut indexes = BTreeMap::<(usize, StructuralEvidenceKind), StructuralEvidenceIndex>::new();
    let mut reports = Vec::new();
    for source in sources {
        let mut record_count = 0u64;
        let mut matched_record_count = 0u64;
        let mut boundary_spanning_record_count = 0u64;
        let mut source_kmers = HashSet::new();
        visit_fasta_records(&source.path, |header, sequence| {
            record_count = record_count.saturating_add(1);
            let Some(symbol) = parse_gene_symbol_from_header(header, gene_symbol_tags) else {
                return Ok(());
            };
            let Some(gene_idx) = gene_lookup.get(&normalize_gene_symbol(&symbol)).copied() else {
                return Ok(());
            };
            let record_id = header
                .split_ascii_whitespace()
                .next()
                .unwrap_or("")
                .trim_start_matches('>');
            let record_kmers = if source.kind.requires_boundary_span() {
                let raw_boundary = header_token_value(header, "boundary_after").ok_or_else(|| {
                    format!(
                        "Structural evidence record '{}' in '{}' requires boundary_after:N in its FASTA header",
                        record_id, source.path
                    )
                })?;
                let boundary_after = raw_boundary.parse::<usize>().map_err(|e| {
                    format!(
                        "Invalid boundary_after '{}' for structural evidence record '{}' in '{}': {e}",
                        raw_boundary, record_id, source.path
                    )
                })?;
                let kmers =
                    canonical_kmers_crossing_boundary(sequence.as_bytes(), k, boundary_after)?;
                if kmers.is_empty() {
                    return Err(format!(
                        "Structural evidence record '{}' in '{}' has no valid boundary-spanning {}-mers",
                        record_id, source.path, k
                    ));
                }
                boundary_spanning_record_count = boundary_spanning_record_count.saturating_add(1);
                kmers
            } else {
                canonical_kmer_set(sequence.as_bytes(), k)
            };
            matched_record_count = matched_record_count.saturating_add(1);
            let index =
                indexes
                    .entry((gene_idx, source.kind))
                    .or_insert_with(|| StructuralEvidenceIndex {
                        gene_symbol: requested_genes[gene_idx].clone(),
                        kind: source.kind,
                        record_ids: BTreeSet::new(),
                        kmers: HashSet::new(),
                        boundary_spanning_only: source.kind.requires_boundary_span(),
                    });
            if !record_id.is_empty() {
                index.record_ids.insert(record_id.to_string());
            }
            source_kmers.extend(record_kmers.iter().copied());
            index.kmers.extend(record_kmers);
            Ok(())
        })?;
        reports.push(StructuralEvidenceSourceReport {
            kind: source.kind,
            path: source.path.clone(),
            bytes: file_size(&source.path)?,
            record_count,
            matched_record_count,
            boundary_spanning_record_count,
            distinct_kmers: source_kmers.len(),
        });
    }
    Ok((indexes.into_values().collect(), reports))
}

fn canonical_kmer_set(sequence: &[u8], k: usize) -> HashSet<u64> {
    let mut kmers = HashSet::new();
    canonical_kmers_for_each(sequence, k, |kmer| {
        kmers.insert(kmer);
    });
    kmers
}

fn canonical_kmers_crossing_boundary(
    sequence: &[u8],
    k: usize,
    boundary_after: usize,
) -> Result<HashSet<u64>, String> {
    if boundary_after == 0 || boundary_after >= sequence.len() {
        return Err(format!(
            "boundary_after must be between 1 and sequence_length-1 (got {boundary_after} for length {})",
            sequence.len()
        ));
    }
    if k == 0 || k > 31 || sequence.len() < k {
        return Ok(HashSet::new());
    }
    let min_start = boundary_after.saturating_sub(k.saturating_sub(1));
    let max_start = boundary_after
        .saturating_sub(1)
        .min(sequence.len().saturating_sub(k));
    let mut kmers = HashSet::new();
    for start in min_start..=max_start {
        canonical_kmers_for_each(&sequence[start..start + k], k, |kmer| {
            kmers.insert(kmer);
        });
    }
    Ok(kmers)
}

fn build_manifest(gene_indices: &[GeneIndex], k: usize) -> Vec<GeneManifestEntry> {
    gene_indices
        .iter()
        .map(|gene| GeneManifestEntry {
            gene_symbol: gene.symbol.clone(),
            requested: true,
            present: !gene.transcript_ids.is_empty(),
            transcript_count: gene.transcript_ids.len(),
            distinct_kmers: gene.kmers.len(),
            total_kmer_positions: gene.total_kmer_positions,
            kmer_len: k,
            transcript_ids: gene.transcript_ids.iter().cloned().collect(),
        })
        .collect()
}

fn target_transcript_lookup(gene_indices: &[GeneIndex]) -> HashSet<String> {
    let mut transcripts = HashSet::new();
    for gene in gene_indices {
        for transcript_id in &gene.transcript_ids {
            transcripts.insert(transcript_id.clone());
            transcripts.insert(strip_transcript_version(transcript_id).to_string());
        }
    }
    transcripts
}

fn build_universes(
    request: &TargetRescueRequest,
    target_transcripts: &HashSet<String>,
    allowlist: Option<&HashSet<String>>,
    gene_count: usize,
    structural_count: usize,
    salmon_bridge: bool,
) -> Result<Vec<UniverseState>, String> {
    if salmon_bridge {
        let mut universes = Vec::new();
        if let Some(path) = &request.salmon_unmapped_names {
            let mut ids = read_id_set_from_path(path)?;
            apply_allowlist_filter(&mut ids, allowlist);
            universes.push(UniverseState {
                name: "salmon_unmapped".to_string(),
                selected_ids: Some(ids),
                tally: UniverseTally::default(),
                gene_tallies: vec![GeneTally::default(); gene_count],
                structural_tallies: vec![StructuralEvidenceTally::default(); structural_count],
            });
        }
        if let Some(path) = &request.salmon_mappings_sam {
            let mut ids = target_mapped_read_ids_from_sam(path, target_transcripts)?;
            apply_allowlist_filter(&mut ids, allowlist);
            universes.push(UniverseState {
                name: "target_mapped".to_string(),
                selected_ids: Some(ids),
                tally: UniverseTally::default(),
                gene_tallies: vec![GeneTally::default(); gene_count],
                structural_tallies: vec![StructuralEvidenceTally::default(); structural_count],
            });
        }
        if universes.is_empty() {
            return Err(
                "rescue-screen Salmon bridge requires --salmon-unmapped-names or --salmon-mappings-sam"
                    .to_string(),
            );
        }
        return Ok(universes);
    }
    Ok(vec![UniverseState {
        name: if allowlist.is_some() {
            "filtered".to_string()
        } else {
            "all".to_string()
        },
        selected_ids: allowlist.cloned(),
        tally: UniverseTally::default(),
        gene_tallies: vec![GeneTally::default(); gene_count],
        structural_tallies: vec![StructuralEvidenceTally::default(); structural_count],
    }])
}

fn apply_allowlist_filter(ids: &mut HashSet<String>, allowlist: Option<&HashSet<String>>) {
    if let Some(allowlist) = allowlist {
        ids.retain(|id| allowlist.contains(id));
    }
}

fn output_files_for_prefix(prefix: &str, salmon_bridge: bool) -> TargetRescueOutputFiles {
    TargetRescueOutputFiles {
        gene_manifest_tsv: format!("{prefix}.gene_manifest.tsv"),
        gene_hits_tsv: format!("{prefix}.gene_hits.tsv"),
        read_hits_tsv: format!("{prefix}.read_hits.tsv"),
        structural_hits_tsv: format!("{prefix}.structural_hits.tsv"),
        run_metadata_json: format!("{prefix}.run_metadata.json"),
        summary_json: format!("{prefix}.summary.json"),
        comparison_json: salmon_bridge.then(|| format!("{prefix}.comparison.json")),
    }
}

fn ensure_output_parent(prefix: &str) -> Result<(), String> {
    let parent = Path::new(prefix).parent().unwrap_or_else(|| Path::new("."));
    fs::create_dir_all(parent).map_err(|e| {
        format!(
            "Could not create output directory '{}': {e}",
            parent.display()
        )
    })
}

fn create_file(path: &str) -> Result<File, String> {
    File::create(path).map_err(|e| format!("Could not create output file '{path}': {e}"))
}

fn write_json_file<T: Serialize>(path: &str, value: &T) -> Result<(), String> {
    let file = create_file(path)?;
    serde_json::to_writer_pretty(BufWriter::new(file), value)
        .map_err(|e| format!("Could not write JSON file '{path}': {e}"))
}

fn build_universe_summaries(universes: &[UniverseState]) -> Vec<UniverseSummary> {
    universes
        .iter()
        .map(|universe| UniverseSummary {
            universe: universe.name.clone(),
            selected_reads: universe.tally.selected_reads,
            evaluated_reads: universe.tally.evaluated_reads,
            reads_matching_any: universe.tally.reads_matching_any,
            ambiguous_reads: universe.tally.ambiguous_reads,
            selected_input_reads: universe.tally.selected_input_reads,
            selected_evidence_units: universe.tally.selected_reads,
            evaluated_evidence_units: universe.tally.evaluated_reads,
            evidence_units_matching_any: universe.tally.reads_matching_any,
            ambiguous_evidence_units: universe.tally.ambiguous_reads,
        })
        .collect()
}

fn build_gene_hit_summaries(
    universes: &[UniverseState],
    gene_indices: &[GeneIndex],
) -> Vec<GeneHitSummary> {
    let mut rows = Vec::new();
    for universe in universes {
        for (idx, gene) in gene_indices.iter().enumerate() {
            let tally = &universe.gene_tallies[idx];
            rows.push(GeneHitSummary {
                universe: universe.name.clone(),
                gene_symbol: gene.symbol.clone(),
                present: !gene.transcript_ids.is_empty(),
                reads_hit: tally.reads_hit,
                reads_hit_unique: tally.reads_hit_unique,
                reads_hit_ambiguous: tally.reads_hit_ambiguous,
                total_matching_kmers: tally.total_matching_kmers,
            });
        }
    }
    rows
}

fn build_structural_hit_summaries(
    universes: &[UniverseState],
    indexes: &[StructuralEvidenceIndex],
) -> Vec<StructuralEvidenceHitSummary> {
    let mut rows = Vec::new();
    for universe in universes {
        for (index_idx, index) in indexes.iter().enumerate() {
            let tally = &universe.structural_tallies[index_idx];
            rows.push(StructuralEvidenceHitSummary {
                universe: universe.name.clone(),
                gene_symbol: index.gene_symbol.clone(),
                evidence_kind: index.kind,
                candidate_label: index.kind.candidate_label().to_string(),
                record_count: index.record_ids.len(),
                distinct_kmers: index.kmers.len(),
                reads_hit: tally.reads_hit,
                total_matching_kmers: tally.total_matching_kmers,
                boundary_spanning_only: index.boundary_spanning_only,
                rna_anchored_hits: tally.rna_anchored_hits,
                unanchored_hits: tally.unanchored_hits,
            });
        }
    }
    rows
}

fn build_comparison_report(summary: &TargetRescueSummary) -> TargetRescueComparison {
    let by_universe_gene = summary
        .gene_hits
        .iter()
        .map(|row| {
            (
                (row.universe.as_str(), row.gene_symbol.as_str()),
                row.reads_hit,
            )
        })
        .collect::<HashMap<_, _>>();
    let gene_rows = summary
        .requested_genes
        .iter()
        .map(|gene| {
            let salmon_unmapped = by_universe_gene
                .get(&("salmon_unmapped", gene.as_str()))
                .copied()
                .unwrap_or(0);
            let target_mapped = by_universe_gene
                .get(&("target_mapped", gene.as_str()))
                .copied()
                .unwrap_or(0);
            ComparisonGeneRow {
                gene_symbol: gene.clone(),
                salmon_unmapped_reads_hit: salmon_unmapped,
                target_mapped_reads_hit: target_mapped,
                delta_unmapped_minus_target_mapped: salmon_unmapped as i64 - target_mapped as i64,
            }
        })
        .collect();
    TargetRescueComparison {
        schema: COMPARISON_SCHEMA.to_string(),
        universes: summary.universes.clone(),
        gene_rows,
    }
}

fn write_gene_manifest_tsv(path: &str, manifest: &[GeneManifestEntry]) -> Result<(), String> {
    let mut writer = BufWriter::new(create_file(path)?);
    writeln!(
        writer,
        "gene_symbol\trequested\tpresent\ttranscript_count\tdistinct_kmers\ttotal_kmer_positions\tkmer_len\ttranscript_ids"
    )
    .map_err(|e| format!("Could not write gene manifest TSV '{path}': {e}"))?;
    for row in manifest {
        writeln!(
            writer,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            tsv_cell(&row.gene_symbol),
            row.requested,
            row.present,
            row.transcript_count,
            row.distinct_kmers,
            row.total_kmer_positions,
            row.kmer_len,
            tsv_cell(&row.transcript_ids.join(";"))
        )
        .map_err(|e| format!("Could not write gene manifest TSV '{path}': {e}"))?;
    }
    writer
        .flush()
        .map_err(|e| format!("Could not flush gene manifest TSV '{path}': {e}"))
}

fn write_gene_hits_tsv(path: &str, gene_hits: &[GeneHitSummary]) -> Result<(), String> {
    let mut writer = BufWriter::new(create_file(path)?);
    writeln!(
        writer,
        "universe\tgene_symbol\tpresent\treads_hit\treads_hit_unique\treads_hit_ambiguous\ttotal_matching_kmers"
    )
    .map_err(|e| format!("Could not write gene hits TSV '{path}': {e}"))?;
    for row in gene_hits {
        writeln!(
            writer,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}",
            tsv_cell(&row.universe),
            tsv_cell(&row.gene_symbol),
            row.present,
            row.reads_hit,
            row.reads_hit_unique,
            row.reads_hit_ambiguous,
            row.total_matching_kmers
        )
        .map_err(|e| format!("Could not write gene hits TSV '{path}': {e}"))?;
    }
    writer
        .flush()
        .map_err(|e| format!("Could not flush gene hits TSV '{path}': {e}"))
}

fn write_read_hits_header(writer: &mut dyn Write) -> Result<(), String> {
    writeln!(
        writer,
        "universe\tread_id\tsource_file\tevidence_unit\tinput_read_count\tgene_symbol\tmatching_kmer_count\tread_kmer_count\tread_length\tambiguous\thit_genes"
    )
    .map_err(|e| format!("Could not write read hits TSV header: {e}"))
}

#[allow(clippy::too_many_arguments)]
fn write_read_hit_row(
    writer: &mut dyn Write,
    universe: &str,
    read_id: &str,
    source_file: &str,
    evidence_unit: TargetRescueEvidenceUnit,
    input_read_count: u8,
    gene_symbol: &str,
    matching_kmer_count: u64,
    read_kmer_count: u64,
    read_length: usize,
    ambiguous: bool,
    hit_genes: &str,
) -> Result<(), String> {
    writeln!(
        writer,
        "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
        tsv_cell(universe),
        tsv_cell(read_id),
        tsv_cell(source_file),
        evidence_unit.as_str(),
        input_read_count,
        tsv_cell(gene_symbol),
        matching_kmer_count,
        read_kmer_count,
        read_length,
        ambiguous,
        tsv_cell(hit_genes)
    )
    .map_err(|e| format!("Could not write read hits TSV row: {e}"))
}

fn write_structural_hits_header(writer: &mut dyn Write) -> Result<(), String> {
    writeln!(
        writer,
        "universe\tread_id\tsource_file\tevidence_unit\tinput_read_count\tgene_symbol\tevidence_kind\tcandidate_label\tboundary_spanning_only\trna_anchored\tmatching_kmer_count\tread_kmer_count\tread_length\tevidence_record_ids"
    )
    .map_err(|e| format!("Could not write structural hits TSV header: {e}"))
}

#[allow(clippy::too_many_arguments)]
fn write_structural_hit_row(
    writer: &mut dyn Write,
    universe: &str,
    read_id: &str,
    source_file: &str,
    evidence_unit: TargetRescueEvidenceUnit,
    input_read_count: u8,
    evidence: &StructuralEvidenceIndex,
    matching_kmer_count: u64,
    read_kmer_count: u64,
    read_length: usize,
    rna_anchored: bool,
) -> Result<(), String> {
    writeln!(
        writer,
        "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
        tsv_cell(universe),
        tsv_cell(read_id),
        tsv_cell(source_file),
        evidence_unit.as_str(),
        input_read_count,
        tsv_cell(&evidence.gene_symbol),
        evidence.kind.as_str(),
        evidence.kind.candidate_label(),
        evidence.boundary_spanning_only,
        rna_anchored,
        matching_kmer_count,
        read_kmer_count,
        read_length,
        tsv_cell(
            &evidence
                .record_ids
                .iter()
                .cloned()
                .collect::<Vec<_>>()
                .join(";")
        )
    )
    .map_err(|e| format!("Could not write structural hits TSV row: {e}"))
}

fn tsv_cell(raw: &str) -> String {
    raw.replace('\\', "\\\\")
        .replace('\t', "\\t")
        .replace('\r', "\\r")
        .replace('\n', "\\n")
}

fn count_gene_kmer_matches(
    sequence: &str,
    k: usize,
    gene_indices: &[GeneIndex],
) -> (u64, Vec<u64>) {
    let mut read_kmer_count = 0u64;
    let mut gene_counts = vec![0u64; gene_indices.len()];
    canonical_kmers_for_each(sequence.as_bytes(), k, |kmer| {
        read_kmer_count = read_kmer_count.saturating_add(1);
        for (idx, gene) in gene_indices.iter().enumerate() {
            if gene.kmers.contains(&kmer) {
                gene_counts[idx] = gene_counts[idx].saturating_add(1);
            }
        }
    });
    (read_kmer_count, gene_counts)
}

fn count_structural_kmer_matches(
    sequence: &str,
    k: usize,
    indexes: &[StructuralEvidenceIndex],
) -> Vec<u64> {
    let mut counts = vec![0u64; indexes.len()];
    canonical_kmers_for_each(sequence.as_bytes(), k, |kmer| {
        for (idx, evidence) in indexes.iter().enumerate() {
            if evidence.kmers.contains(&kmer) {
                counts[idx] = counts[idx].saturating_add(1);
            }
        }
    });
    counts
}

pub(crate) fn canonical_kmers_for_each(sequence: &[u8], k: usize, mut on_kmer: impl FnMut(u64)) {
    if k == 0 || k > 31 {
        return;
    }
    let mask = if k == 32 {
        u64::MAX
    } else {
        (1u64 << (2 * k)) - 1
    };
    let shift = 2 * (k - 1);
    let mut forward = 0u64;
    let mut reverse = 0u64;
    let mut valid = 0usize;
    for base in sequence {
        let Some(bits) = base_bits(*base) else {
            forward = 0;
            reverse = 0;
            valid = 0;
            continue;
        };
        forward = ((forward << 2) | u64::from(bits)) & mask;
        let rc_bits = u64::from(3u8.saturating_sub(bits));
        reverse = (reverse >> 2) | (rc_bits << shift);
        valid = valid.saturating_add(1).min(k);
        if valid == k {
            on_kmer(forward.min(reverse));
        }
    }
}

fn base_bits(base: u8) -> Option<u8> {
    match base.to_ascii_uppercase() {
        b'A' => Some(0),
        b'C' => Some(1),
        b'G' => Some(2),
        b'T' | b'U' => Some(3),
        _ => None,
    }
}

fn parse_gene_symbol_from_header(header: &str, tags: &[String]) -> Option<String> {
    let mut values = HashMap::<String, String>::new();
    for token in header.split_ascii_whitespace() {
        let Some((key, value)) = token.split_once(':') else {
            continue;
        };
        values
            .entry(key.to_ascii_lowercase())
            .or_insert_with(|| value.trim().to_string());
    }
    for tag in tags {
        if let Some(value) = values.get(&tag.to_ascii_lowercase())
            && !value.trim().is_empty()
        {
            return Some(value.trim().to_string());
        }
    }
    None
}

fn header_token_value<'a>(header: &'a str, requested_key: &str) -> Option<&'a str> {
    header.split_ascii_whitespace().find_map(|token| {
        let (key, value) = token.split_once(':')?;
        key.eq_ignore_ascii_case(requested_key)
            .then_some(value.trim())
            .filter(|value| !value.is_empty())
    })
}

fn structural_interpretation_warnings(sources: &[StructuralEvidenceSource]) -> Vec<String> {
    let kinds = sources
        .iter()
        .map(|source| source.kind)
        .collect::<BTreeSet<_>>();
    let mut warnings = Vec::new();
    if kinds.contains(&StructuralEvidenceKind::ExonIntronBoundary)
        || kinds.contains(&StructuralEvidenceKind::Intron)
    {
        warnings.push(
            "Retained-intron and intronic hits remain compatible with incompletely processed pre-mRNA or genomic-DNA contamination. rna_anchored_hits require an annotated splice-junction signal for the same gene in the same read or fragment, but are still candidate evidence rather than a retained-intron call."
                .to_string(),
        );
    }
    if kinds.contains(&StructuralEvidenceKind::GenomicRegion) {
        warnings.push(
            "A genomic-region catalog was explicitly enabled; genomic-region hits are not RNA-specific and should not be interpreted as transcript evidence."
                .to_string(),
        );
    }
    warnings
}

fn strip_transcript_version(transcript_id: &str) -> &str {
    transcript_id
        .split_once('.')
        .map_or(transcript_id, |(id, _)| id)
}

fn normalize_read_id(raw: &str) -> String {
    let first = raw
        .trim()
        .trim_start_matches('@')
        .trim_start_matches('>')
        .split_ascii_whitespace()
        .next()
        .unwrap_or("");
    first
        .strip_suffix("/1")
        .or_else(|| first.strip_suffix("/2"))
        .unwrap_or(first)
        .to_string()
}

pub(crate) fn read_id_set_from_path(path: &str) -> Result<HashSet<String>, String> {
    let mut reader = open_maybe_gz_reader(path)?;
    let mut ids = HashSet::new();
    let mut line = String::new();
    loop {
        line.clear();
        let bytes = reader
            .read_line(&mut line)
            .map_err(|e| format!("Could not read ID list '{path}': {e}"))?;
        if bytes == 0 {
            break;
        }
        let trimmed = line.trim();
        if trimmed.is_empty() || trimmed.starts_with('#') {
            continue;
        }
        let id = normalize_read_id(trimmed);
        if !id.is_empty() {
            ids.insert(id);
        }
    }
    Ok(ids)
}

fn target_mapped_read_ids_from_sam(
    path: &str,
    target_transcripts: &HashSet<String>,
) -> Result<HashSet<String>, String> {
    let mut reader = open_maybe_gz_reader(path)?;
    let mut ids = HashSet::new();
    let mut line = String::new();
    loop {
        line.clear();
        let bytes = reader
            .read_line(&mut line)
            .map_err(|e| format!("Could not read SAM mappings '{path}': {e}"))?;
        if bytes == 0 {
            break;
        }
        let trimmed = line.trim();
        if trimmed.is_empty() || trimmed.starts_with('@') {
            continue;
        }
        let fields = trimmed.split('\t').collect::<Vec<_>>();
        if fields.len() < 3 {
            continue;
        }
        let flag = fields[1].parse::<u32>().unwrap_or(0);
        if flag & 0x4 != 0 {
            continue;
        }
        let rname = fields[2];
        if rname == "*" {
            continue;
        }
        if target_transcripts.contains(rname)
            || target_transcripts.contains(strip_transcript_version(rname))
        {
            let id = normalize_read_id(fields[0]);
            if !id.is_empty() {
                ids.insert(id);
            }
        }
    }
    Ok(ids)
}

pub(crate) fn visit_fasta_records(
    path: &str,
    mut on_record: impl FnMut(&str, &str) -> Result<(), String>,
) -> Result<(), String> {
    let mut reader = open_maybe_gz_reader(path)?;
    let mut line = String::new();
    let mut header: Option<String> = None;
    let mut sequence = String::new();
    loop {
        line.clear();
        let bytes = reader
            .read_line(&mut line)
            .map_err(|e| format!("Could not read FASTA '{path}': {e}"))?;
        if bytes == 0 {
            break;
        }
        let trimmed = line.trim_end_matches(['\r', '\n']);
        if let Some(next_header) = trimmed.strip_prefix('>') {
            if let Some(current_header) = header.replace(next_header.to_string()) {
                on_record(&current_header, &sequence)?;
                sequence.clear();
            }
        } else {
            sequence.push_str(trimmed.trim());
        }
    }
    if let Some(current_header) = header {
        on_record(&current_header, &sequence)?;
    }
    Ok(())
}

pub(crate) fn visit_read_records(
    path: &str,
    mut on_record: impl FnMut(SequenceRecord) -> Result<(), String>,
) -> Result<(), String> {
    let mut records = ReadRecordStream::open(path)?;
    while let Some(record) = records.next_record()? {
        on_record(record)?;
    }
    Ok(())
}

pub(crate) fn visit_paired_read_records(
    read1_path: &str,
    read2_path: &str,
    mut on_fragment: impl FnMut(SequenceRecord) -> Result<(), String>,
) -> Result<(), String> {
    let mut read1 = ReadRecordStream::open(read1_path)?;
    let mut read2 = ReadRecordStream::open(read2_path)?;
    loop {
        match (read1.next_record()?, read2.next_record()?) {
            (None, None) => return Ok(()),
            (Some(_), None) | (None, Some(_)) => {
                return Err(format!(
                    "Paired read files '{read1_path}' and '{read2_path}' contain different record counts"
                ));
            }
            (Some(left), Some(right)) => {
                if left.id != right.id {
                    return Err(format!(
                        "Paired read files are out of sync: '{}' in '{read1_path}' does not match '{}' in '{read2_path}'",
                        left.id, right.id
                    ));
                }
                let mut sequence = left.sequence;
                sequence.push('N');
                sequence.push_str(&right.sequence);
                on_fragment(SequenceRecord {
                    id: left.id,
                    sequence,
                    source_length: left.source_length.saturating_add(right.source_length),
                })?;
            }
        }
    }
}

fn detect_read_format(path: &str, reader: &mut dyn BufRead) -> Result<ReadFormat, String> {
    let lower = strip_gz_suffix(path).to_ascii_lowercase();
    if lower.ends_with(".fa")
        || lower.ends_with(".fasta")
        || lower.ends_with(".fna")
        || lower.ends_with(".fas")
    {
        return Ok(ReadFormat::Fasta);
    }
    if lower.ends_with(".fq") || lower.ends_with(".fastq") {
        return Ok(ReadFormat::Fastq);
    }
    let buffer = reader
        .fill_buf()
        .map_err(|e| format!("Could not sniff read format for '{path}': {e}"))?;
    let first = buffer
        .iter()
        .copied()
        .find(|byte| !byte.is_ascii_whitespace())
        .ok_or_else(|| format!("Could not sniff read format for empty file '{path}'"))?;
    match first {
        b'>' => Ok(ReadFormat::Fasta),
        b'@' => Ok(ReadFormat::Fastq),
        _ => Err(format!(
            "Could not detect read format for '{path}': expected FASTA '>' or FASTQ '@'"
        )),
    }
}

fn strip_gz_suffix(path: &str) -> &str {
    path.strip_suffix(".gz")
        .or_else(|| path.strip_suffix(".GZ"))
        .unwrap_or(path)
}

pub(crate) fn open_maybe_gz_reader(path: &str) -> Result<Box<dyn BufRead>, String> {
    let file = File::open(path).map_err(|e| format!("Could not open '{path}': {e}"))?;
    if path.to_ascii_lowercase().ends_with(".gz") {
        Ok(Box::new(BufReader::new(MultiGzDecoder::new(file))))
    } else {
        Ok(Box::new(BufReader::new(file)))
    }
}

fn file_size(path: &str) -> Result<u64, String> {
    fs::metadata(path)
        .map(|meta| meta.len())
        .map_err(|e| format!("Could not inspect '{path}': {e}"))
}

#[cfg(test)]
mod tests {
    use super::*;
    use flate2::{Compression, write::GzEncoder};
    use std::collections::BTreeMap;
    use tempfile::{TempDir, tempdir};

    const GENEA_SEQ: &str =
        "CTCAGATGGCTTAGCCAGCCTACGCGTGCCACTTGCGAAACTTCAAGTATCAAACTTGGTCTGGCGCATCACACGTTGCG";
    const GENEB_SEQ: &str =
        "TAAAGTGTGGCACCTGCAACGCACCCAAATTCTCCCAGACCCCGTGATTAAGAGACGTAGGACGACCTACGGCTGAGCCG";

    fn canonical_set(sequence: &str, k: usize) -> HashSet<u64> {
        let mut kmers = HashSet::new();
        canonical_kmers_for_each(sequence.as_bytes(), k, |kmer| {
            kmers.insert(kmer);
        });
        kmers
    }

    fn valid_kmer_count(sequence: &str, k: usize) -> u64 {
        let mut count = 0u64;
        canonical_kmers_for_each(sequence.as_bytes(), k, |_| {
            count += 1;
        });
        count
    }

    fn reverse_complement(sequence: &str) -> String {
        sequence
            .bytes()
            .rev()
            .map(|base| match base.to_ascii_uppercase() {
                b'A' => 'T',
                b'C' => 'G',
                b'G' => 'C',
                b'T' | b'U' => 'A',
                _ => 'N',
            })
            .collect()
    }

    fn write_text(path: &Path, text: &str) {
        fs::write(path, text).expect("test fixture should write");
    }

    fn gzipped_committed_reads(base: &Path, dir: &TempDir) -> String {
        let source = base.join("mini_reads.fastq");
        assert!(
            source.exists(),
            "committed target_rescue FASTQ fixture is required"
        );
        let reads = fs::read(&source).expect("committed FASTQ fixture should read");
        let gz_path = dir.path().join("mini_reads.fastq.gz");
        let file = File::create(&gz_path).expect("temporary gzipped FASTQ should create");
        let mut encoder = GzEncoder::new(file, Compression::default());
        encoder
            .write_all(&reads)
            .expect("temporary gzipped FASTQ should write");
        encoder
            .finish()
            .expect("temporary gzipped FASTQ should finish");
        gz_path.display().to_string()
    }

    fn synthetic_salmon_bridge_inputs(dir: &TempDir) -> (String, String) {
        let unmapped_names = dir.path().join("salmon_unmapped_names.txt");
        write_text(
            &unmapped_names,
            "readA_geneA u\nreadN_nohit u\nreadAB_ambiguous u\nreadA2_geneA u\n",
        );
        let mappings_sam = dir.path().join("salmon_mappings.sam");
        write_text(
            &mappings_sam,
            concat!(
                "@HD\tVN:1.0\tSO:unsorted\n",
                "@SQ\tSN:ENSTB1.1\tLN:80\n",
                "@SQ\tSN:ENSTC1.1\tLN:80\n",
                "readB_geneB\t0\tENSTB1.1\t9\t255\t48M\t*\t0\t0\t",
                "GGCACCTGCAACGCACCCAAATTCTCCCAGACCCCGTGATTAAGAGAC\t*\n",
                "readA_geneA\t0\tENSTC1.1\t1\t255\t44M\t*\t0\t0\t",
                "TTAGCCAGCCTACGCGTGCCACTTGCGAAACTTCAAGTATCAAA\t*\n",
            ),
        );
        (
            unmapped_names.display().to_string(),
            mappings_sam.display().to_string(),
        )
    }

    fn synthetic_fasta() -> String {
        format!(
            ">ENSTA1.1 cdna gene_symbol:GENEA description:test\n{GENEA_SEQ}\n>ENSTB1.1 cdna gene_symbol:GENEB description:test\n{GENEB_SEQ}\n"
        )
    }

    fn synthetic_fastq() -> String {
        "@readA_geneA\nTTAGCCAGCCTACGCGTGCCACTTGCGAAACTTCAAGTATCAAA\n+\nIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII\n@readAB_ambiguous\nCTCAGATGGCTTAGCCAGCCTACGCGTGCCTAAAGTGTGGCACCTGCAACGCACCCAAAT\n+\nIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII\n@readN_nohit\nGCGGAGGAAGCTGATCTCTGTTCCCTGTAACCGCCGACCTTGCCCACGATGGGGCGCATA\n+\nIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII\n".to_string()
    }

    fn run_synthetic(prefix: &str, reads: String) -> (TempDir, TargetRescueSummary) {
        let dir = tempdir().expect("temp dir");
        let fasta = dir.path().join("transcripts.fasta");
        let fastq = dir.path().join("reads.fastq");
        write_text(&fasta, &synthetic_fasta());
        write_text(&fastq, &reads);
        let summary = run_target_rescue_screen(TargetRescueRequest {
            transcript_fastas: vec![fasta.display().to_string()],
            genes: vec![
                "GENEA".to_string(),
                "GENEB".to_string(),
                "MISSINGX".to_string(),
            ],
            reads: vec![fastq.display().to_string()],
            kmer_len: 11,
            min_kmer_hits: 3,
            output_prefix: dir.path().join(prefix).display().to_string(),
            ..TargetRescueRequest::default()
        })
        .expect("synthetic rescue screen should run");
        (dir, summary)
    }

    #[test]
    fn canonical_kmers_are_strand_agnostic() {
        let sequence = "ACGTTGCAACGTAC";
        let reverse = reverse_complement(sequence);
        assert_eq!(canonical_set(sequence, 7), canonical_set(&reverse, 7));
    }

    #[test]
    fn kmers_spanning_n_are_skipped() {
        assert_eq!(valid_kmer_count("AAAAANAAAAA", 5), 2);
        assert_eq!(valid_kmer_count("NNNNN", 5), 0);
    }

    #[test]
    fn synthetic_gene_hashes_are_disjoint_for_fixture_k() {
        let gene_a = canonical_set(GENEA_SEQ, 11);
        let gene_b = canonical_set(GENEB_SEQ, 11);
        assert!(gene_a.is_disjoint(&gene_b));
    }

    #[test]
    fn missing_gene_is_reported_and_ambiguous_read_hits_both_genes() {
        let (_dir, summary) = run_synthetic("rescue", synthetic_fastq());
        assert_eq!(summary.missing_genes, vec!["MISSINGX"]);
        let present = summary
            .manifest
            .iter()
            .map(|entry| (entry.gene_symbol.as_str(), entry.present))
            .collect::<BTreeMap<_, _>>();
        assert_eq!(present.get("GENEA"), Some(&true));
        assert_eq!(present.get("GENEB"), Some(&true));
        assert_eq!(present.get("MISSINGX"), Some(&false));

        let read_hits = fs::read_to_string(&summary.output_files.read_hits_tsv)
            .expect("read hits TSV should exist");
        let ambiguous_rows = read_hits
            .lines()
            .filter(|line| line.contains("readAB_ambiguous"))
            .collect::<Vec<_>>();
        assert_eq!(ambiguous_rows.len(), 2);
        assert!(ambiguous_rows.iter().all(|line| line.contains("\ttrue\t")));
        assert!(
            ambiguous_rows
                .iter()
                .all(|line| line.ends_with("\tGENEA;GENEB"))
        );
        let pure_rows = read_hits
            .lines()
            .filter(|line| line.contains("readA_geneA"))
            .collect::<Vec<_>>();
        assert_eq!(pure_rows.len(), 1);
        assert!(pure_rows[0].contains("\tfalse\tGENEA"));
    }

    #[test]
    fn structural_catalogs_emit_candidate_evidence_without_novelty_claims() {
        let dir = tempdir().expect("temp dir");
        let transcripts = dir.path().join("transcripts.fasta");
        let reads = dir.path().join("reads.fastq");
        let exons = dir.path().join("exons.fasta");
        let boundaries = dir.path().join("boundaries.fasta");
        write_text(&transcripts, &synthetic_fasta());
        write_text(&reads, &synthetic_fastq());
        write_text(
            &exons,
            &format!(">GENEA_exon_1 gene_symbol:GENEA\n{}\n", &GENEA_SEQ[..36]),
        );
        write_text(
            &boundaries,
            &format!(
                ">GENEB_exon_intron_1 gene_symbol:GENEB boundary_after:18\n{}\n",
                &GENEB_SEQ[..36]
            ),
        );

        let summary = run_target_rescue_screen(TargetRescueRequest {
            transcript_fastas: vec![transcripts.display().to_string()],
            genes: vec!["GENEA".to_string(), "GENEB".to_string()],
            reads: vec![reads.display().to_string()],
            structural_evidence: vec![
                StructuralEvidenceSource {
                    kind: StructuralEvidenceKind::Exon,
                    path: exons.display().to_string(),
                },
                StructuralEvidenceSource {
                    kind: StructuralEvidenceKind::ExonIntronBoundary,
                    path: boundaries.display().to_string(),
                },
            ],
            kmer_len: 11,
            min_kmer_hits: 3,
            output_prefix: dir.path().join("structural").display().to_string(),
            ..TargetRescueRequest::default()
        })
        .expect("structural evidence run should pass");

        assert_eq!(summary.structural_hits.len(), 2);
        let exon = summary
            .structural_hits
            .iter()
            .find(|row| row.evidence_kind == StructuralEvidenceKind::Exon)
            .expect("exon summary");
        assert_eq!(exon.gene_symbol, "GENEA");
        assert_eq!(exon.candidate_label, "known_exonic");
        assert_eq!(exon.reads_hit, 2);
        let boundary = summary
            .structural_hits
            .iter()
            .find(|row| row.evidence_kind == StructuralEvidenceKind::ExonIntronBoundary)
            .expect("boundary summary");
        assert_eq!(boundary.gene_symbol, "GENEB");
        assert_eq!(boundary.candidate_label, "retained_intron_candidate");
        assert_eq!(boundary.reads_hit, 1);
        assert!(boundary.boundary_spanning_only);
        assert_eq!(boundary.rna_anchored_hits, 0);
        assert_eq!(boundary.unanchored_hits, 1);
        assert!(
            summary
                .warnings
                .iter()
                .any(|warning| warning.contains("genomic-DNA contamination"))
        );

        let tsv = fs::read_to_string(&summary.output_files.structural_hits_tsv)
            .expect("structural hit TSV should exist");
        assert!(tsv.contains("\tknown_exonic\t"));
        assert!(tsv.contains("\tretained_intron_candidate\t"));
        assert!(!tsv.contains("novel_exon"));
    }

    #[test]
    fn boundary_catalog_indexes_only_kmers_that_cross_the_declared_boundary() {
        let sequence = b"ACGTCAGTACGATTGGCCTAACGT";
        let boundary_kmers = canonical_kmers_crossing_boundary(sequence, 7, 12)
            .expect("declared boundary should be valid");
        let left_flank = canonical_set("ACGTCAGTACGA", 7);
        let crossing_read = canonical_set("GTACGATTGGCC", 7);

        assert!(!boundary_kmers.is_empty());
        assert!(boundary_kmers.is_disjoint(&left_flank));
        assert!(!boundary_kmers.is_disjoint(&crossing_read));
    }

    #[test]
    fn paired_fragments_distinguish_rna_anchored_and_unanchored_boundary_hits() {
        let dir = tempdir().expect("temp dir");
        let transcripts = dir.path().join("transcripts.fasta");
        let boundaries = dir.path().join("boundaries.fasta");
        let junctions = dir.path().join("junctions.fasta");
        let read1 = dir.path().join("reads_1.fastq");
        let read2 = dir.path().join("reads_2.fastq");
        write_text(
            &transcripts,
            &format!(">ENSTA1.1 cdna gene_symbol:GENEA\n{GENEA_SEQ}\n"),
        );
        write_text(
            &boundaries,
            ">GENEA_retained_boundary gene_symbol:GENEA boundary_after:12\nACGTCAGTACGATTGGCCTAACGT\n",
        );
        write_text(
            &junctions,
            ">GENEA_known_junction gene_symbol:GENEA boundary_after:12\nGATCCGATAGCACCGTTAAGGCTT\n",
        );
        write_text(
            &read1,
            concat!(
                "@rna_anchored/1\nGTACGATTGGCC\n+\nIIIIIIIIIIII\n",
                "@boundary_only/1\nGTACGATTGGCC\n+\nIIIIIIIIIIII\n",
                "@flank_only/1\nACGTCAGTACGA\n+\nIIIIIIIIIIII\n",
            ),
        );
        write_text(
            &read2,
            concat!(
                "@rna_anchored/2\nATAGCACCGTTA\n+\nIIIIIIIIIIII\n",
                "@boundary_only/2\nNNNNNNNNNNNN\n+\nIIIIIIIIIIII\n",
                "@flank_only/2\nNNNNNNNNNNNN\n+\nIIIIIIIIIIII\n",
            ),
        );

        let summary = run_target_rescue_screen(TargetRescueRequest {
            transcript_fastas: vec![transcripts.display().to_string()],
            genes: vec!["GENEA".to_string()],
            reads: Vec::new(),
            read_pairs: vec![TargetRescueReadPairInput {
                read1: read1.display().to_string(),
                read2: read2.display().to_string(),
            }],
            structural_evidence: vec![
                StructuralEvidenceSource {
                    kind: StructuralEvidenceKind::AnnotatedJunction,
                    path: junctions.display().to_string(),
                },
                StructuralEvidenceSource {
                    kind: StructuralEvidenceKind::ExonIntronBoundary,
                    path: boundaries.display().to_string(),
                },
            ],
            kmer_len: 7,
            min_kmer_hits: 1,
            output_prefix: dir.path().join("paired").display().to_string(),
            ..TargetRescueRequest::default()
        })
        .expect("paired structural rescue should run");

        assert_eq!(summary.total_input_reads, 6);
        assert_eq!(summary.total_evidence_units, 3);
        assert_eq!(summary.universes[0].selected_reads, 3);
        assert_eq!(summary.universes[0].selected_input_reads, 6);
        assert_eq!(summary.universes[0].selected_evidence_units, 3);
        let boundary = summary
            .structural_hits
            .iter()
            .find(|row| row.evidence_kind == StructuralEvidenceKind::ExonIntronBoundary)
            .expect("boundary summary");
        assert_eq!(boundary.reads_hit, 2);
        assert_eq!(boundary.rna_anchored_hits, 1);
        assert_eq!(boundary.unanchored_hits, 1);
        let junction = summary
            .structural_hits
            .iter()
            .find(|row| row.evidence_kind == StructuralEvidenceKind::AnnotatedJunction)
            .expect("junction summary");
        assert_eq!(junction.reads_hit, 1);
        assert_eq!(junction.rna_anchored_hits, 1);

        let tsv = fs::read_to_string(&summary.output_files.structural_hits_tsv)
            .expect("structural TSV should exist");
        let retained_rows = tsv
            .lines()
            .filter(|line| line.contains("retained_intron_candidate"))
            .collect::<Vec<_>>();
        assert_eq!(retained_rows.len(), 2);
        assert!(retained_rows.iter().any(|line| {
            line.contains("rna_anchored")
                && line.contains("\tfragment\t2\t")
                && line.contains("\ttrue\t")
        }));
        assert!(retained_rows.iter().any(|line| {
            line.contains("boundary_only")
                && line.contains("\tfragment\t2\t")
                && line.contains("\tfalse\t")
        }));
        assert!(!tsv.contains("flank_only\t"));
    }

    #[test]
    fn legacy_request_defaults_to_no_paired_inputs() {
        let request: TargetRescueRequest = serde_json::from_value(serde_json::json!({
            "transcript_fastas": ["transcripts.fa"],
            "genes": ["GENEA"],
            "reads": ["reads.fq"],
            "structural_evidence": [],
            "kmer_len": 25,
            "min_kmer_hits": 3,
            "gene_symbol_tags": ["gene_symbol"],
            "output_prefix": "out/rescue"
        }))
        .expect("legacy request should deserialize");
        assert!(request.read_pairs.is_empty());

        let paired_only: TargetRescueRequest = serde_json::from_value(serde_json::json!({
            "transcript_fastas": ["transcripts.fa"],
            "genes": ["GENEA"],
            "read_pairs": [{"read1": "reads_1.fq", "read2": "reads_2.fq"}],
            "structural_evidence": [],
            "kmer_len": 25,
            "min_kmer_hits": 3,
            "gene_symbol_tags": ["gene_symbol"],
            "output_prefix": "out/rescue"
        }))
        .expect("paired-only request should deserialize without an empty reads field");
        assert!(paired_only.reads.is_empty());
        assert_eq!(paired_only.read_pairs.len(), 1);
    }

    #[test]
    fn committed_fixture_gzip_allowlist_counts_streamed_and_selected_reads() {
        let base = Path::new("test_files/fixtures/target_rescue");
        let dir = tempdir().expect("temp dir");
        let gz_reads = gzipped_committed_reads(base, &dir);
        let allowlist = dir.path().join("allowlist.txt");
        write_text(&allowlist, "readA_geneA\nreadA2_geneA\n");
        let prefix = dir.path().join("fixture_allowlist");
        let summary = run_target_rescue_screen(TargetRescueRequest {
            transcript_fastas: vec![base.join("mini_transcripts.fasta").display().to_string()],
            genes: vec![
                "GENEA".to_string(),
                "GENEB".to_string(),
                "MISSINGX".to_string(),
            ],
            reads: vec![gz_reads],
            read_id_allowlist: Some(allowlist.display().to_string()),
            kmer_len: 11,
            min_kmer_hits: 3,
            output_prefix: prefix.display().to_string(),
            ..TargetRescueRequest::default()
        })
        .expect("fixture gzip allowlist run should pass");
        assert_eq!(summary.total_input_reads, 5);
        assert_eq!(summary.universes[0].universe, "filtered");
        assert_eq!(summary.universes[0].selected_reads, 2);
        assert_eq!(summary.universes[0].reads_matching_any, 2);
        let gene_hits = summary
            .gene_hits
            .iter()
            .map(|row| (row.gene_symbol.as_str(), row.reads_hit))
            .collect::<BTreeMap<_, _>>();
        assert_eq!(gene_hits.get("GENEA"), Some(&2));
        assert_eq!(gene_hits.get("GENEB"), Some(&0));
    }

    #[test]
    fn committed_fixture_salmon_bridge_builds_two_universes_and_comparison() {
        let base = Path::new("test_files/fixtures/target_rescue");
        let dir = tempdir().expect("temp dir");
        let gz_reads = gzipped_committed_reads(base, &dir);
        let (salmon_unmapped_names, salmon_mappings_sam) = synthetic_salmon_bridge_inputs(&dir);
        let prefix = dir.path().join("fixture_salmon");
        let summary = run_target_rescue_screen(TargetRescueRequest {
            transcript_fastas: vec![base.join("mini_transcripts.fasta").display().to_string()],
            genes: vec![
                "GENEA".to_string(),
                "GENEB".to_string(),
                "MISSINGX".to_string(),
            ],
            reads: vec![gz_reads],
            salmon_unmapped_names: Some(salmon_unmapped_names),
            salmon_mappings_sam: Some(salmon_mappings_sam),
            kmer_len: 11,
            min_kmer_hits: 3,
            output_prefix: prefix.display().to_string(),
            ..TargetRescueRequest::default()
        })
        .expect("fixture Salmon bridge run should pass");
        assert_eq!(
            summary
                .universes
                .iter()
                .map(|row| row.universe.as_str())
                .collect::<Vec<_>>(),
            vec!["salmon_unmapped", "target_mapped"]
        );
        assert_eq!(summary.universes[0].selected_reads, 4);
        assert_eq!(summary.universes[0].ambiguous_reads, 1);
        assert_eq!(summary.universes[1].selected_reads, 1);
        let unmapped_gene_a = summary
            .gene_hits
            .iter()
            .find(|row| row.universe == "salmon_unmapped" && row.gene_symbol == "GENEA")
            .expect("GENEA salmon_unmapped row");
        assert_eq!(unmapped_gene_a.reads_hit, 3);
        let target_gene_b = summary
            .gene_hits
            .iter()
            .find(|row| row.universe == "target_mapped" && row.gene_symbol == "GENEB")
            .expect("GENEB target_mapped row");
        assert_eq!(target_gene_b.reads_hit, 1);
        let comparison_path = summary
            .output_files
            .comparison_json
            .as_ref()
            .expect("comparison report path");
        let comparison = fs::read_to_string(comparison_path).expect("comparison JSON should exist");
        assert!(comparison.contains(COMPARISON_SCHEMA));
    }
}
