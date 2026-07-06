//! Allele-aware RNA-read hash screen over transcript-coordinate variants.
//!
//! The screen materializes reference, hap1, and hap2 transcript FASTAs from a
//! local transcript FASTA plus an explicit transcript-coordinate variant table,
//! then classifies reads against haplotype-specific k-mer sets. It is a small
//! deterministic evidence screen: it reports sequence support for allele
//! imbalance, but it does not call biological significance.

use crate::target_rescue::{
    canonical_kmers_for_each, read_id_set_from_path, visit_fasta_records, visit_read_records,
};
use serde::{Deserialize, Serialize};
use std::{
    collections::{BTreeMap, BTreeSet, HashMap, HashSet},
    fs::{self, File},
    io::{BufRead, BufReader, BufWriter, Write},
    path::Path,
};

pub const ALLELE_HASH_SCREEN_SCHEMA: &str = "gentle.rna_allele_hash_screen.v1";
const DEFAULT_KMER_LEN: usize = 21;
const DEFAULT_MIN_UNIQUE_KMER_HITS: u64 = 1;
const GENE_SYMBOL_TAGS: &[&str] = &["gene_symbol", "gene_name", "symbol", "gene"];

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
    pub read_files: Vec<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub read_id_allowlist: Option<String>,
    pub out_dir: String,
    pub kmer_len: usize,
    pub min_unique_kmer_hits: u64,
}

impl Default for AlleleHashScreenRequest {
    fn default() -> Self {
        Self {
            gene: String::new(),
            transcript_fasta: String::new(),
            variant_table: None,
            vcf: None,
            transcript_map: None,
            read_files: Vec::new(),
            read_id_allowlist: None,
            out_dir: String::new(),
            kmer_len: DEFAULT_KMER_LEN,
            min_unique_kmer_hits: DEFAULT_MIN_UNIQUE_KMER_HITS,
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
    pub output_files: AlleleHashScreenOutputFiles,
    pub haplotype_fastas: Vec<HaplotypeFastaReport>,
    pub transcript_summaries: Vec<AlleleTranscriptSummary>,
    pub variant_summaries: Vec<AlleleVariantSummary>,
    pub classification_counts: BTreeMap<String, u64>,
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
    pub read_files: Vec<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub read_id_allowlist: Option<String>,
    pub kmer_len: usize,
    pub min_unique_kmer_hits: u64,
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
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum AlleleHashPhaseMode {
    Phased,
    UnphasedAlleleLevelOnly,
}

impl AlleleHashPhaseMode {
    pub fn as_str(self) -> &'static str {
        match self {
            Self::Phased => "phased",
            Self::UnphasedAlleleLevelOnly => "unphased_allele_level_only",
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum AlleleReadClassification {
    Hap1,
    Hap2,
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
    pub classification: AlleleReadClassification,
    pub read_length: usize,
    pub read_kmer_count: u64,
    pub reference_hits: u64,
    pub hap1_hits: u64,
    pub hap2_hits: u64,
    pub hap1_unique_hits: u64,
    pub hap2_unique_hits: u64,
    pub reference_only_hits: u64,
    pub matched_transcripts: Vec<String>,
    pub supporting_variants: Vec<String>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct AlleleTranscriptSummary {
    pub transcript_id: String,
    pub length_nt: usize,
    pub variant_count: usize,
    pub informative_reads: u64,
    pub hap1_reads: u64,
    pub hap2_reads: u64,
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
    reference_only: HashSet<u64>,
    hap1_unique: HashSet<u64>,
    hap2_unique: HashSet<u64>,
}

#[derive(Debug, Clone)]
struct TranscriptKmerIndex {
    transcript_id: String,
    all: HashSet<u64>,
}

#[derive(Debug, Clone)]
struct VariantSupportKmers {
    variant_id: String,
    transcript_id: String,
    reference_only: HashSet<u64>,
    hap1_unique: HashSet<u64>,
    hap2_unique: HashSet<u64>,
}

#[derive(Debug, Clone, Default)]
struct CountSummary {
    informative_reads: u64,
    hap1_reads: u64,
    hap2_reads: u64,
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
    let variant_table = request.variant_table.as_deref().ok_or_else(|| {
        "allele-hash-screen requires --variant-table transcript-coordinate TSV in Phase 1"
            .to_string()
    })?;
    let variants = load_variant_table(variant_table, &transcripts)?;
    if variants.iter().any(|variant| !variant.phased) {
        warnings.push(
            "At least one variant is unphased; report is unphased, allele-level only and does not fabricate haplotype phase."
                .to_string(),
        );
    }
    let phase_mode = if variants.iter().all(|variant| variant.phased) {
        AlleleHashPhaseMode::Phased
    } else {
        AlleleHashPhaseMode::UnphasedAlleleLevelOnly
    };
    apply_variants_to_haplotypes(&mut transcripts, &variants, &mut warnings)?;
    let out_dir = Path::new(&request.out_dir);
    fs::create_dir_all(out_dir).map_err(|e| {
        format!(
            "Could not create output directory '{}': {e}",
            out_dir.display()
        )
    })?;
    let output_files = output_files(out_dir, &request.gene);
    write_haplotype_fastas(&transcripts, &output_files)?;

    let kmer_index = build_kmer_index(&transcripts, request.kmer_len);
    let transcript_indexes = build_transcript_kmer_indexes(&transcripts, request.kmer_len);
    let variant_support = build_variant_support_kmers(&transcripts, &variants, request.kmer_len);
    let allowlist = request
        .read_id_allowlist
        .as_deref()
        .map(read_id_set_from_path)
        .transpose()?;
    let mut read_count_total = 0u64;
    let mut read_count_selected = 0u64;
    let mut calls = Vec::<AlleleReadCall>::new();
    for read_file in &request.read_files {
        visit_read_records(read_file, |record| {
            read_count_total = read_count_total.saturating_add(1);
            if let Some(allowlist) = &allowlist
                && !allowlist.contains(&record.id)
            {
                return Ok(());
            }
            read_count_selected = read_count_selected.saturating_add(1);
            let call = classify_read(
                record.id,
                read_file,
                &record.sequence,
                request.kmer_len,
                request.min_unique_kmer_hits,
                &kmer_index,
                &transcript_indexes,
                &variant_support,
            );
            calls.push(call);
            Ok(())
        })?;
    }

    let mut transcript_summaries = summarize_transcripts(&transcripts, &variants, &calls);
    let mut variant_summaries = summarize_variants(&variants, &variant_support, &calls);
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
    let classification_counts = classification_counts(&calls);
    write_read_calls_tsv(&output_files.reads_tsv, &calls)?;
    let report = AlleleHashScreenReport {
        schema: ALLELE_HASH_SCREEN_SCHEMA.to_string(),
        gene: request.gene.clone(),
        phase_mode,
        params: AlleleHashScreenParams {
            transcript_fasta: request.transcript_fasta.clone(),
            variant_table: request.variant_table.clone(),
            vcf: request.vcf.clone(),
            transcript_map: request.transcript_map.clone(),
            read_files: request.read_files.clone(),
            read_id_allowlist: request.read_id_allowlist.clone(),
            kmer_len: request.kmer_len,
            min_unique_kmer_hits: request.min_unique_kmer_hits,
        },
        transcript_count: transcripts.len(),
        variant_count: variants.len(),
        read_count_total,
        read_count_selected,
        output_files: output_files.clone(),
        haplotype_fastas: vec![
            HaplotypeFastaReport {
                haplotype: "reference".to_string(),
                path: output_files.reference_fasta.clone(),
                transcript_count: transcripts.len(),
            },
            HaplotypeFastaReport {
                haplotype: "hap1".to_string(),
                path: output_files.hap1_fasta.clone(),
                transcript_count: transcripts.len(),
            },
            HaplotypeFastaReport {
                haplotype: "hap2".to_string(),
                path: output_files.hap2_fasta.clone(),
                transcript_count: transcripts.len(),
            },
        ],
        transcript_summaries,
        variant_summaries,
        classification_counts,
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
    if request.variant_table.is_none() {
        if request.vcf.is_some() {
            return Err(
                "allele-hash-screen VCF input requires an explicit transcript-coordinate mapping; use --variant-table for the deterministic Phase 1 path"
                    .to_string(),
            );
        }
        return Err("allele-hash-screen requires --variant-table PATH".to_string());
    }
    if request.vcf.is_some() {
        return Err(
            "allele-hash-screen VCF projection is not implemented in Phase 1; provide --variant-table transcript-coordinate TSV"
                .to_string(),
        );
    }
    if request.read_files.is_empty() {
        return Err("allele-hash-screen requires at least one --read-file PATH".to_string());
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
        let assignment = parse_genotype_assignment(&genotype);
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
                .filter(|value| !value.is_empty()),
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

fn parse_genotype_assignment(genotype: &str) -> (bool, bool, bool) {
    let trimmed = genotype.trim();
    if let Some((left, right)) = trimmed.split_once('|') {
        return (left.trim() == "1", right.trim() == "1", true);
    }
    if let Some((left, right)) = trimmed.split_once('/') {
        let has_alt = left.trim() == "1" || right.trim() == "1";
        return (has_alt, has_alt, false);
    }
    let has_alt = trimmed == "1";
    (has_alt, has_alt, false)
}

fn apply_variants_to_haplotypes(
    transcripts: &mut [TranscriptRecord],
    variants: &[VariantRecord],
    warnings: &mut Vec<String>,
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
            warnings.push(format!(
                "Variant '{}' expected ref '{}' at {}:{}, but transcript FASTA has '{}'; using supplied alternate only where requested.",
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

fn build_kmer_index(transcripts: &[TranscriptRecord], k: usize) -> KmerIndex {
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
    index
}

fn build_transcript_kmer_indexes(
    transcripts: &[TranscriptRecord],
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
            TranscriptKmerIndex {
                transcript_id: transcript.id.clone(),
                all,
            }
        })
        .collect()
}

fn build_variant_support_kmers(
    transcripts: &[TranscriptRecord],
    variants: &[VariantRecord],
    k: usize,
) -> Vec<VariantSupportKmers> {
    let by_id = transcripts
        .iter()
        .map(|record| (record.id.as_str(), record))
        .collect::<HashMap<_, _>>();
    variants
        .iter()
        .filter_map(|variant| {
            let transcript = by_id.get(variant.transcript_id.as_str())?;
            let pos = variant.cdna_pos_1based.checked_sub(1)?;
            let reference = overlapping_kmers(transcript.reference.as_bytes(), pos, k);
            let hap1 = overlapping_kmers(&transcript.hap1, pos, k);
            let hap2 = overlapping_kmers(&transcript.hap2, pos, k);
            let reference_only = reference
                .difference(&hap1)
                .copied()
                .collect::<HashSet<_>>()
                .difference(&hap2)
                .copied()
                .collect();
            let hap1_unique = hap1
                .difference(&hap2)
                .copied()
                .collect::<HashSet<_>>()
                .difference(&reference)
                .copied()
                .collect();
            let hap2_unique = hap2
                .difference(&hap1)
                .copied()
                .collect::<HashSet<_>>()
                .difference(&reference)
                .copied()
                .collect();
            Some(VariantSupportKmers {
                variant_id: variant.variant_id.clone(),
                transcript_id: variant.transcript_id.clone(),
                reference_only,
                hap1_unique,
                hap2_unique,
            })
        })
        .collect()
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

fn classify_read(
    read_id: String,
    source_file: &str,
    sequence: &str,
    k: usize,
    min_unique_hits: u64,
    index: &KmerIndex,
    transcript_indexes: &[TranscriptKmerIndex],
    variant_support: &[VariantSupportKmers],
) -> AlleleReadCall {
    let normalized = normalize_sequence(sequence);
    let mut read_kmers = HashSet::<u64>::new();
    canonical_kmers_for_each(normalized.as_bytes(), k, |kmer| {
        read_kmers.insert(kmer);
    });
    let read_kmer_count = read_kmers.len() as u64;
    let reference_hits = intersection_count(&read_kmers, &index.reference);
    let hap1_hits = intersection_count(&read_kmers, &index.hap1);
    let hap2_hits = intersection_count(&read_kmers, &index.hap2);
    let hap1_unique_hits = intersection_count(&read_kmers, &index.hap1_unique);
    let hap2_unique_hits = intersection_count(&read_kmers, &index.hap2_unique);
    let reference_only_hits = intersection_count(&read_kmers, &index.reference_only);
    let classification = if read_kmer_count == 0 {
        AlleleReadClassification::Uninformative
    } else if hap1_unique_hits >= min_unique_hits && hap2_unique_hits < min_unique_hits {
        AlleleReadClassification::Hap1
    } else if hap2_unique_hits >= min_unique_hits && hap1_unique_hits < min_unique_hits {
        AlleleReadClassification::Hap2
    } else if hap1_unique_hits >= min_unique_hits && hap2_unique_hits >= min_unique_hits {
        AlleleReadClassification::Ambiguous
    } else if reference_only_hits >= min_unique_hits {
        AlleleReadClassification::ReferenceOnly
    } else if reference_hits >= min_unique_hits
        || hap1_hits >= min_unique_hits
        || hap2_hits >= min_unique_hits
    {
        AlleleReadClassification::Ambiguous
    } else {
        AlleleReadClassification::OffTarget
    };
    let supporting_variants =
        supporting_variants_for_read(&read_kmers, classification, variant_support);
    let matched_transcripts = matched_transcripts_for_read(&read_kmers, transcript_indexes);
    AlleleReadCall {
        read_id,
        source_file: source_file.to_string(),
        classification,
        read_length: normalized.len(),
        read_kmer_count,
        reference_hits,
        hap1_hits,
        hap2_hits,
        hap1_unique_hits,
        hap2_unique_hits,
        reference_only_hits,
        matched_transcripts,
        supporting_variants,
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

fn supporting_variants_for_read(
    read_kmers: &HashSet<u64>,
    classification: AlleleReadClassification,
    variant_support: &[VariantSupportKmers],
) -> Vec<String> {
    let mut ids = BTreeSet::new();
    for variant in variant_support {
        let supported = match classification {
            AlleleReadClassification::Hap1 => !variant.hap1_unique.is_disjoint(read_kmers),
            AlleleReadClassification::Hap2 => !variant.hap2_unique.is_disjoint(read_kmers),
            AlleleReadClassification::ReferenceOnly => {
                !variant.reference_only.is_disjoint(read_kmers)
            }
            AlleleReadClassification::Ambiguous => {
                !variant.hap1_unique.is_disjoint(read_kmers)
                    || !variant.hap2_unique.is_disjoint(read_kmers)
                    || !variant.reference_only.is_disjoint(read_kmers)
            }
            AlleleReadClassification::Uninformative | AlleleReadClassification::OffTarget => false,
        };
        if supported {
            ids.insert(variant.variant_id.clone());
        }
    }
    ids.into_iter().collect()
}

fn summarize_transcripts(
    transcripts: &[TranscriptRecord],
    variants: &[VariantRecord],
    calls: &[AlleleReadCall],
) -> Vec<AlleleTranscriptSummary> {
    let mut summaries = transcripts
        .iter()
        .map(|transcript| (transcript.id.as_str(), CountSummary::default()))
        .collect::<BTreeMap<_, _>>();
    for call in calls {
        for transcript_id in &call.matched_transcripts {
            if let Some(summary) = summaries.get_mut(transcript_id.as_str()) {
                add_classification_count(summary, call.classification);
            }
        }
    }
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
    support: &[VariantSupportKmers],
    calls: &[AlleleReadCall],
) -> Vec<AlleleVariantSummary> {
    let mut counts = variants
        .iter()
        .map(|variant| (variant.variant_id.as_str(), CountSummary::default()))
        .collect::<BTreeMap<_, _>>();
    let variant_transcripts = support
        .iter()
        .map(|row| (row.variant_id.as_str(), row.transcript_id.as_str()))
        .collect::<BTreeMap<_, _>>();
    for call in calls {
        for variant_id in &call.supporting_variants {
            if let Some(summary) = counts.get_mut(variant_id.as_str()) {
                add_classification_count(summary, call.classification);
            }
        }
    }
    variants
        .iter()
        .map(|variant| {
            let count = counts
                .remove(variant.variant_id.as_str())
                .unwrap_or_default();
            let _transcript_id = variant_transcripts
                .get(variant.variant_id.as_str())
                .copied()
                .unwrap_or(variant.transcript_id.as_str());
            variant_summary(variant, count)
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

fn classification_counts(calls: &[AlleleReadCall]) -> BTreeMap<String, u64> {
    let mut counts = BTreeMap::<String, u64>::new();
    for classification in [
        AlleleReadClassification::Hap1,
        AlleleReadClassification::Hap2,
        AlleleReadClassification::ReferenceOnly,
        AlleleReadClassification::Ambiguous,
        AlleleReadClassification::Uninformative,
        AlleleReadClassification::OffTarget,
    ] {
        counts.insert(classification.as_str().to_string(), 0);
    }
    for call in calls {
        *counts
            .entry(call.classification.as_str().to_string())
            .or_insert(0) += 1;
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

fn write_read_calls_tsv(path: &str, calls: &[AlleleReadCall]) -> Result<(), String> {
    let mut writer = BufWriter::new(
        File::create(path).map_err(|e| format!("Could not create read TSV '{path}': {e}"))?,
    );
    writeln!(
        writer,
        "read_id\tsource_file\tclassification\tread_length\tread_kmer_count\treference_hits\thap1_hits\thap2_hits\thap1_unique_hits\thap2_unique_hits\treference_only_hits\tmatched_transcripts\tsupporting_variants"
    )
    .map_err(|e| format!("Could not write read TSV '{path}': {e}"))?;
    for call in calls {
        writeln!(
            writer,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            tsv_cell(&call.read_id),
            tsv_cell(&call.source_file),
            call.classification.as_str(),
            call.read_length,
            call.read_kmer_count,
            call.reference_hits,
            call.hap1_hits,
            call.hap2_hits,
            call.hap1_unique_hits,
            call.hap2_unique_hits,
            call.reference_only_hits,
            tsv_cell(&call.matched_transcripts.join(";")),
            tsv_cell(&call.supporting_variants.join(";"))
        )
        .map_err(|e| format!("Could not write read TSV '{path}': {e}"))?;
    }
    writer
        .flush()
        .map_err(|e| format!("Could not flush read TSV '{path}': {e}"))
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
    }
}
