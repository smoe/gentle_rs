//! UCSC RepeatMasker (`rmsk`) table resource contracts.
//!
//! This module keeps the row parser, portable resource descriptor, and index
//! recommendations together so GUI, CLI, and future genome-track ingestion
//! paths can agree on the same `rmsk` semantics without re-parsing UCSC table
//! details in adapter code.

use serde::{Deserialize, Serialize};
use std::collections::BTreeMap;
use std::io::{BufRead, Write};
use std::path::Path;

pub const UCSC_RMSK_RESOURCE_SCHEMA: &str = "gentle.ucsc_rmsk_resource.v1";
pub const UCSC_RMSK_DESCRIPTOR_SCHEMA: &str = "gentle.ucsc_rmsk_descriptor.v1";
pub const UCSC_RMSK_INTERVAL_INDEX_SCHEMA: &str = "gentle.ucsc_rmsk_interval_index.v1";
pub const DEFAULT_UCSC_RMSK_ASSEMBLY: &str = "hg38";
pub const DEFAULT_UCSC_RMSK_RESOURCE_PATH: &str = "data/resources/ucsc.rmsk.hg38.json";
pub const DEFAULT_UCSC_RMSK_INDEX_PATH: &str = "data/resources/ucsc.rmsk.hg38.interval-index.json";
pub const UCSC_RMSK_TABLE: &str = "rmsk";

pub const UCSC_RMSK_FIELD_ORDER: &[&str] = &[
    "bin",
    "swScore",
    "milliDiv",
    "milliDel",
    "milliIns",
    "genoName",
    "genoStart",
    "genoEnd",
    "genoLeft",
    "strand",
    "repName",
    "repClass",
    "repFamily",
    "repStart",
    "repEnd",
    "repLeft",
    "id",
];

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq)]
pub struct UcscRmskFieldSpec {
    pub name: String,
    pub sql_type: String,
    pub description: String,
}

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq)]
pub struct UcscRmskIndexRecommendation {
    pub index_id: String,
    pub label: String,
    pub storage: String,
    pub key_fields: Vec<String>,
    pub query_pattern: String,
    pub notes: Vec<String>,
}

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq)]
pub struct UcscRmskResourceDescriptor {
    pub schema: String,
    pub provider: String,
    pub assembly_database: String,
    pub table: String,
    pub download_url: String,
    pub schema_url: String,
    pub track_url: String,
    pub coordinate_system: String,
    pub field_order: Vec<String>,
    pub fields: Vec<UcscRmskFieldSpec>,
    pub index_recommendations: Vec<UcscRmskIndexRecommendation>,
}

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq)]
pub struct UcscRmskRecord {
    pub bin: u32,
    #[serde(rename = "swScore")]
    pub sw_score: u32,
    #[serde(rename = "milliDiv")]
    pub milli_div: u32,
    #[serde(rename = "milliDel")]
    pub milli_del: u32,
    #[serde(rename = "milliIns")]
    pub milli_ins: u32,
    #[serde(rename = "genoName")]
    pub geno_name: String,
    #[serde(rename = "genoStart")]
    pub geno_start: u64,
    #[serde(rename = "genoEnd")]
    pub geno_end: u64,
    #[serde(rename = "genoLeft")]
    pub geno_left: i64,
    pub strand: String,
    #[serde(rename = "repName")]
    pub rep_name: String,
    #[serde(rename = "repClass")]
    pub rep_class: String,
    #[serde(rename = "repFamily")]
    pub rep_family: String,
    #[serde(rename = "repStart")]
    pub rep_start: i64,
    #[serde(rename = "repEnd")]
    pub rep_end: i64,
    #[serde(rename = "repLeft")]
    pub rep_left: i64,
    pub id: String,
}

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq)]
pub struct UcscRmskResourceSnapshot {
    pub schema: String,
    pub source: String,
    pub provider: String,
    pub assembly_database: String,
    pub table: String,
    pub coordinate_system: String,
    #[serde(default)]
    pub download_url: String,
    #[serde(default)]
    pub schema_url: String,
    #[serde(default)]
    pub track_url: String,
    #[serde(default)]
    pub fetched_at_unix_ms: u128,
    #[serde(default)]
    pub field_order: Vec<String>,
    #[serde(default)]
    pub fields: Vec<UcscRmskFieldSpec>,
    #[serde(default)]
    pub index_recommendations: Vec<UcscRmskIndexRecommendation>,
    #[serde(default)]
    pub records: Vec<UcscRmskRecord>,
    #[serde(default)]
    pub row_count: usize,
    #[serde(default)]
    pub truncated: bool,
    #[serde(default)]
    pub warnings: Vec<String>,
}

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq)]
pub struct UcscRmskIndexedInterval {
    pub source_row: usize,
    pub bin: u32,
    pub chromosome: String,
    pub start_0based: usize,
    pub end_0based_exclusive: usize,
    pub strand: String,
    pub rep_name: String,
    pub rep_class: String,
    pub rep_family: String,
    pub normalized_alias: String,
    #[serde(rename = "swScore")]
    pub sw_score: u32,
    #[serde(rename = "milliDiv")]
    pub milli_div: u32,
    #[serde(rename = "milliDel")]
    pub milli_del: u32,
    #[serde(rename = "milliIns")]
    pub milli_ins: u32,
    pub id: String,
}

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq)]
pub struct UcscRmskIntervalIndex {
    pub schema: String,
    pub source_resource: String,
    pub provider: String,
    pub assembly_database: String,
    pub table: String,
    pub coordinate_system: String,
    pub row_count: usize,
    pub chromosome_count: usize,
    pub chromosomes: BTreeMap<String, Vec<UcscRmskIndexedInterval>>,
    pub class_family_index: BTreeMap<String, Vec<usize>>,
    pub repeat_name_index: BTreeMap<String, Vec<usize>>,
}

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq)]
pub struct UcscRmskIndexSummary {
    pub output: String,
    pub source_resource: String,
    pub assembly_database: String,
    pub row_count: usize,
    pub chromosome_count: usize,
}

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq)]
pub struct UcscRmskProjectionAnchor {
    pub chromosome: String,
    pub start_1based: usize,
    pub end_1based: usize,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub strand: Option<char>,
}

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq)]
pub struct UcscRmskProjectedOverlap {
    pub interval: UcscRmskIndexedInterval,
    pub local_start_0based: usize,
    pub local_end_0based_exclusive: usize,
    pub local_strand: String,
    pub genomic_start_0based: usize,
    pub genomic_end_0based_exclusive: usize,
    pub overlap_bp: usize,
    pub clipped: bool,
}

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq)]
pub struct UcscRmskSyncSummary {
    pub output: String,
    pub source: String,
    pub assembly_database: String,
    pub row_count: usize,
    pub truncated: bool,
}

pub fn ucsc_rmsk_download_url(assembly_database: &str) -> String {
    format!("https://hgdownload.soe.ucsc.edu/goldenPath/{assembly_database}/database/rmsk.txt.gz")
}

pub fn ucsc_rmsk_schema_url(assembly_database: &str) -> String {
    format!("https://hgdownload.soe.ucsc.edu/goldenPath/{assembly_database}/database/rmsk.sql")
}

pub fn ucsc_rmsk_track_url(assembly_database: &str) -> String {
    format!("https://genome.ucsc.edu/cgi-bin/hgTrackUi?db={assembly_database}&g=rmsk")
}

pub fn default_ucsc_rmsk_resource_path(assembly_database: &str) -> String {
    if assembly_database == DEFAULT_UCSC_RMSK_ASSEMBLY {
        DEFAULT_UCSC_RMSK_RESOURCE_PATH.to_string()
    } else {
        let token = assembly_database
            .chars()
            .map(|ch| {
                if ch.is_ascii_alphanumeric() {
                    ch.to_ascii_lowercase()
                } else {
                    '_'
                }
            })
            .collect::<String>()
            .trim_matches('_')
            .to_string();
        format!("data/resources/ucsc.rmsk.{token}.json")
    }
}

pub fn default_ucsc_rmsk_index_path(assembly_database: &str) -> String {
    if assembly_database == DEFAULT_UCSC_RMSK_ASSEMBLY {
        DEFAULT_UCSC_RMSK_INDEX_PATH.to_string()
    } else {
        let token = assembly_database
            .chars()
            .map(|ch| {
                if ch.is_ascii_alphanumeric() {
                    ch.to_ascii_lowercase()
                } else {
                    '_'
                }
            })
            .collect::<String>()
            .trim_matches('_')
            .to_string();
        format!("data/resources/ucsc.rmsk.{token}.interval-index.json")
    }
}

pub fn normalized_repeat_alias(rep_name: &str, rep_class: &str, rep_family: &str) -> String {
    let class = rep_class.trim();
    let family = rep_family.trim();
    let name = rep_name.trim();
    let class_upper = class.to_ascii_uppercase();
    let family_upper = family.to_ascii_uppercase();
    let name_upper = name.to_ascii_uppercase();
    if class_upper == "SINE" && (family_upper.contains("ALU") || name_upper.contains("ALU")) {
        return "SINE/Alu".to_string();
    }
    if class_upper == "LINE" && (family_upper.starts_with("L1") || name_upper.starts_with("L1")) {
        return "LINE/L1".to_string();
    }
    if class_upper == "LTR" && (family_upper.contains("ERV") || name_upper.contains("ERV")) {
        return "LTR/ERV".to_string();
    }
    if !class.is_empty() && !family.is_empty() {
        format!("{class}/{family}")
    } else if !class.is_empty() {
        class.to_string()
    } else {
        name.to_string()
    }
}

pub fn ucsc_rmsk_field_specs() -> Vec<UcscRmskFieldSpec> {
    [
        (
            "bin",
            "smallint unsigned",
            "UCSC bin value for chromosome range queries.",
        ),
        ("swScore", "int unsigned", "Smith-Waterman alignment score."),
        (
            "milliDiv",
            "int unsigned",
            "Base mismatches in parts per thousand.",
        ),
        (
            "milliDel",
            "int unsigned",
            "Deleted bases in parts per thousand.",
        ),
        (
            "milliIns",
            "int unsigned",
            "Inserted bases in parts per thousand.",
        ),
        ("genoName", "varchar", "Genomic sequence/chromosome name."),
        (
            "genoStart",
            "int unsigned",
            "0-based start in the genomic sequence.",
        ),
        (
            "genoEnd",
            "int unsigned",
            "0-based half-open end in the genomic sequence.",
        ),
        (
            "genoLeft",
            "int",
            "Negative bases remaining after the genomic match.",
        ),
        ("strand", "char", "Relative repeat orientation, '+' or '-'."),
        ("repName", "varchar", "Repeat name."),
        ("repClass", "varchar", "Repeat class."),
        ("repFamily", "varchar", "Repeat family."),
        (
            "repStart",
            "int",
            "Start in repeat consensus, or negative bases after match on reverse strand.",
        ),
        ("repEnd", "int", "End in repeat consensus."),
        (
            "repLeft",
            "int",
            "Negative bases after match, or consensus start on reverse strand.",
        ),
        (
            "id",
            "char",
            "First digit of RepeatMasker .out id; UCSC advises it is best ignored.",
        ),
    ]
    .into_iter()
    .map(|(name, sql_type, description)| UcscRmskFieldSpec {
        name: name.to_string(),
        sql_type: sql_type.to_string(),
        description: description.to_string(),
    })
    .collect()
}

pub fn ucsc_rmsk_index_recommendations() -> Vec<UcscRmskIndexRecommendation> {
    vec![
        UcscRmskIndexRecommendation {
            index_id: "interval_by_chrom_bin".to_string(),
            label: "Primary interval lookup".to_string(),
            storage: "per-assembly sidecar with chrom -> sorted intervals, optionally retaining UCSC bin levels".to_string(),
            key_fields: vec![
                "genoName".to_string(),
                "bin".to_string(),
                "genoStart".to_string(),
                "genoEnd".to_string(),
            ],
            query_pattern: "Return repeats overlapping one extracted genomic window.".to_string(),
            notes: vec![
                "This should be the first implemented index because GENtle region/gene extraction asks chromosome-range questions.".to_string(),
                "Keep coordinates in UCSC 0-based half-open form internally and convert only at feature materialization boundaries.".to_string(),
                "A compact per-chromosome sorted Vec plus binary search is sufficient for local windows; add the UCSC bin column as a coarse prefilter when chromosome rows become large.".to_string(),
            ],
        },
        UcscRmskIndexRecommendation {
            index_id: "class_family_partition".to_string(),
            label: "Class/family display partition".to_string(),
            storage: "secondary dictionary keyed by normalized repClass and repFamily to row offsets".to_string(),
            key_fields: vec![
                "repClass".to_string(),
                "repFamily".to_string(),
                "repName".to_string(),
            ],
            query_pattern: "Filter or summarize SINE/LINE/LTR/DNA/simple-repeat families without scanning all rows.".to_string(),
            notes: vec![
                "Normalize uncertain suffixes such as 'DNA?' for grouping while preserving original strings in qualifiers.".to_string(),
                "This index feeds feature-tree grouping, legends, and future repeat-density summaries.".to_string(),
            ],
        },
        UcscRmskIndexRecommendation {
            index_id: "repeat_name_lookup".to_string(),
            label: "Repeat-name lookup".to_string(),
            storage: "dictionary keyed by repName to row offsets and class/family metadata".to_string(),
            key_fields: vec!["repName".to_string()],
            query_pattern: "Find all AluY/L1PA2/etc. intervals or build name-specific statistics.".to_string(),
            notes: vec![
                "Useful for expert panels and for matching external repeat-family annotations to local genomic intervals.".to_string(),
                "Do not make this the only index; name lookup is secondary to region-window lookup in GENtle's current workflows.".to_string(),
            ],
        },
        UcscRmskIndexRecommendation {
            index_id: "quality_metric_columns".to_string(),
            label: "Display-quality metric columns".to_string(),
            storage: "columnar sidecar for swScore, milliDiv, milliDel, milliIns keyed by row offset".to_string(),
            key_fields: vec![
                "swScore".to_string(),
                "milliDiv".to_string(),
                "milliDel".to_string(),
                "milliIns".to_string(),
            ],
            query_pattern: "Shade repeat features or summarize repeat-divergence quality after interval lookup.".to_string(),
            notes: vec![
                "Keep metrics separate from geometry so small region imports do not need to deserialize the entire raw row payload.".to_string(),
                "UCSC uses mismatch/deletion/insertion values for display shading; GENtle can later mirror that without changing the primary interval index.".to_string(),
            ],
        },
    ]
}

pub fn ucsc_rmsk_descriptor(assembly_database: &str) -> UcscRmskResourceDescriptor {
    UcscRmskResourceDescriptor {
        schema: UCSC_RMSK_DESCRIPTOR_SCHEMA.to_string(),
        provider: "UCSC Genome Browser".to_string(),
        assembly_database: assembly_database.to_string(),
        table: UCSC_RMSK_TABLE.to_string(),
        download_url: ucsc_rmsk_download_url(assembly_database),
        schema_url: ucsc_rmsk_schema_url(assembly_database),
        track_url: ucsc_rmsk_track_url(assembly_database),
        coordinate_system: "UCSC 0-based half-open genomic intervals".to_string(),
        field_order: UCSC_RMSK_FIELD_ORDER
            .iter()
            .map(|field| (*field).to_string())
            .collect(),
        fields: ucsc_rmsk_field_specs(),
        index_recommendations: ucsc_rmsk_index_recommendations(),
    }
}

fn parse_u32(raw: &str, field: &str, line_number: usize) -> Result<u32, String> {
    raw.parse::<u32>()
        .map_err(|e| format!("Invalid UCSC rmsk {field} '{raw}' on line {line_number}: {e}"))
}

fn parse_u64(raw: &str, field: &str, line_number: usize) -> Result<u64, String> {
    raw.parse::<u64>()
        .map_err(|e| format!("Invalid UCSC rmsk {field} '{raw}' on line {line_number}: {e}"))
}

fn parse_i64(raw: &str, field: &str, line_number: usize) -> Result<i64, String> {
    raw.parse::<i64>()
        .map_err(|e| format!("Invalid UCSC rmsk {field} '{raw}' on line {line_number}: {e}"))
}

pub fn parse_ucsc_rmsk_record(
    line: &str,
    line_number: usize,
) -> Result<Option<UcscRmskRecord>, String> {
    let trimmed = line.trim();
    if trimmed.is_empty()
        || trimmed.starts_with('#')
        || trimmed.starts_with("track ")
        || trimmed.starts_with("browser ")
    {
        return Ok(None);
    }
    let fields = trimmed.split_whitespace().collect::<Vec<_>>();
    if fields
        .first()
        .is_some_and(|value| value.eq_ignore_ascii_case("bin"))
    {
        return Ok(None);
    }
    if fields.len() != UCSC_RMSK_FIELD_ORDER.len() {
        return Err(format!(
            "UCSC rmsk record on line {line_number} has {} fields; expected {}",
            fields.len(),
            UCSC_RMSK_FIELD_ORDER.len()
        ));
    }
    let geno_start = parse_u64(fields[6], "genoStart", line_number)?;
    let geno_end = parse_u64(fields[7], "genoEnd", line_number)?;
    if geno_end <= geno_start {
        return Err(format!(
            "Invalid UCSC rmsk interval on line {line_number}: genoEnd ({geno_end}) must be greater than genoStart ({geno_start})"
        ));
    }
    let strand = fields[9];
    if !matches!(strand, "+" | "-") {
        return Err(format!(
            "Invalid UCSC rmsk strand '{strand}' on line {line_number}; expected '+' or '-'"
        ));
    }
    Ok(Some(UcscRmskRecord {
        bin: parse_u32(fields[0], "bin", line_number)?,
        sw_score: parse_u32(fields[1], "swScore", line_number)?,
        milli_div: parse_u32(fields[2], "milliDiv", line_number)?,
        milli_del: parse_u32(fields[3], "milliDel", line_number)?,
        milli_ins: parse_u32(fields[4], "milliIns", line_number)?,
        geno_name: fields[5].to_string(),
        geno_start,
        geno_end,
        geno_left: parse_i64(fields[8], "genoLeft", line_number)?,
        strand: strand.to_string(),
        rep_name: fields[10].to_string(),
        rep_class: fields[11].to_string(),
        rep_family: fields[12].to_string(),
        rep_start: parse_i64(fields[13], "repStart", line_number)?,
        rep_end: parse_i64(fields[14], "repEnd", line_number)?,
        rep_left: parse_i64(fields[15], "repLeft", line_number)?,
        id: fields[16].to_string(),
    }))
}

pub fn parse_ucsc_rmsk_records(text: &str) -> Result<Vec<UcscRmskRecord>, String> {
    let mut records = vec![];
    for (idx, line) in text.lines().enumerate() {
        if let Some(record) = parse_ucsc_rmsk_record(line, idx + 1)? {
            records.push(record);
        }
    }
    Ok(records)
}

fn record_to_indexed_interval(
    source_row: usize,
    record: &UcscRmskRecord,
) -> Result<UcscRmskIndexedInterval, String> {
    Ok(UcscRmskIndexedInterval {
        source_row,
        bin: record.bin,
        chromosome: record.geno_name.clone(),
        start_0based: usize::try_from(record.geno_start).map_err(|_| {
            format!(
                "UCSC rmsk row {} has genoStart outside this platform's usize range",
                source_row + 1
            )
        })?,
        end_0based_exclusive: usize::try_from(record.geno_end).map_err(|_| {
            format!(
                "UCSC rmsk row {} has genoEnd outside this platform's usize range",
                source_row + 1
            )
        })?,
        strand: normalize_rmsk_strand(&record.strand),
        rep_name: record.rep_name.clone(),
        rep_class: record.rep_class.clone(),
        rep_family: record.rep_family.clone(),
        normalized_alias: normalized_repeat_alias(
            &record.rep_name,
            &record.rep_class,
            &record.rep_family,
        ),
        sw_score: record.sw_score,
        milli_div: record.milli_div,
        milli_del: record.milli_del,
        milli_ins: record.milli_ins,
        id: record.id.clone(),
    })
}

fn normalize_rmsk_strand(raw: &str) -> String {
    match raw.trim() {
        "-" | "C" | "c" => "-".to_string(),
        "+" => "+".to_string(),
        "" => ".".to_string(),
        other => other.to_string(),
    }
}

fn flip_strand(raw: &str) -> String {
    match normalize_rmsk_strand(raw).as_str() {
        "+" => "-".to_string(),
        "-" => "+".to_string(),
        other => other.to_string(),
    }
}

pub fn build_ucsc_rmsk_interval_index(
    source_resource: &str,
    assembly_database: &str,
    records: &[UcscRmskRecord],
) -> Result<UcscRmskIntervalIndex, String> {
    let mut chromosomes = BTreeMap::<String, Vec<UcscRmskIndexedInterval>>::new();
    let mut class_family_index = BTreeMap::<String, Vec<usize>>::new();
    let mut repeat_name_index = BTreeMap::<String, Vec<usize>>::new();
    for (source_row, record) in records.iter().enumerate() {
        let interval = record_to_indexed_interval(source_row, record)?;
        class_family_index
            .entry(format!("{}/{}", interval.rep_class, interval.rep_family))
            .or_default()
            .push(source_row);
        repeat_name_index
            .entry(interval.rep_name.clone())
            .or_default()
            .push(source_row);
        chromosomes
            .entry(interval.chromosome.clone())
            .or_default()
            .push(interval);
    }
    for rows in chromosomes.values_mut() {
        rows.sort_by(|left, right| {
            left.start_0based
                .cmp(&right.start_0based)
                .then(left.end_0based_exclusive.cmp(&right.end_0based_exclusive))
                .then(left.rep_name.cmp(&right.rep_name))
                .then(left.source_row.cmp(&right.source_row))
        });
    }
    Ok(UcscRmskIntervalIndex {
        schema: UCSC_RMSK_INTERVAL_INDEX_SCHEMA.to_string(),
        source_resource: source_resource.to_string(),
        provider: "UCSC Genome Browser".to_string(),
        assembly_database: assembly_database.to_string(),
        table: UCSC_RMSK_TABLE.to_string(),
        coordinate_system: "UCSC 0-based half-open genomic intervals".to_string(),
        row_count: records.len(),
        chromosome_count: chromosomes.len(),
        chromosomes,
        class_family_index,
        repeat_name_index,
    })
}

impl UcscRmskIntervalIndex {
    pub fn chromosome_intervals(&self, chromosome: &str) -> Option<&[UcscRmskIndexedInterval]> {
        if let Some(rows) = self.chromosomes.get(chromosome) {
            return Some(rows);
        }
        self.chromosomes
            .iter()
            .find(|(key, _)| key.eq_ignore_ascii_case(chromosome))
            .map(|(_, rows)| rows.as_slice())
    }

    pub fn overlapping_intervals(
        &self,
        chromosome: &str,
        start_0based: usize,
        end_0based_exclusive: usize,
    ) -> Vec<&UcscRmskIndexedInterval> {
        if end_0based_exclusive <= start_0based {
            return vec![];
        }
        let Some(rows) = self.chromosome_intervals(chromosome) else {
            return vec![];
        };
        rows.iter()
            .filter(|row| row.end_0based_exclusive > start_0based)
            .take_while(|row| row.start_0based < end_0based_exclusive)
            .filter(|row| row.end_0based_exclusive > start_0based)
            .collect()
    }
}

pub fn read_ucsc_rmsk_resource_snapshot(path: &str) -> Result<UcscRmskResourceSnapshot, String> {
    let text = std::fs::read_to_string(path)
        .map_err(|e| format!("Could not read UCSC rmsk resource snapshot '{path}': {e}"))?;
    let snapshot: UcscRmskResourceSnapshot = serde_json::from_str(&text)
        .map_err(|e| format!("Could not parse UCSC rmsk resource snapshot '{path}': {e}"))?;
    if snapshot.schema != UCSC_RMSK_RESOURCE_SCHEMA {
        return Err(format!(
            "Unexpected UCSC rmsk resource schema '{}' in '{path}'",
            snapshot.schema
        ));
    }
    Ok(snapshot)
}

pub fn read_ucsc_rmsk_interval_index(path: &str) -> Result<UcscRmskIntervalIndex, String> {
    let text = std::fs::read_to_string(path)
        .map_err(|e| format!("Could not read UCSC rmsk interval index '{path}': {e}"))?;
    let index: UcscRmskIntervalIndex = serde_json::from_str(&text)
        .map_err(|e| format!("Could not parse UCSC rmsk interval index '{path}': {e}"))?;
    if index.schema != UCSC_RMSK_INTERVAL_INDEX_SCHEMA {
        return Err(format!(
            "Unexpected UCSC rmsk interval index schema '{}' in '{path}'",
            index.schema
        ));
    }
    Ok(index)
}

pub fn write_ucsc_rmsk_interval_index_from_resource(
    resource_path: &str,
    output: &str,
) -> Result<UcscRmskIndexSummary, String> {
    let snapshot = read_ucsc_rmsk_resource_snapshot(resource_path)?;
    if snapshot.truncated {
        return Err(format!(
            "Refusing to build a UCSC rmsk interval index from truncated snapshot '{resource_path}'"
        ));
    }
    let output_path = Path::new(output);
    let parent = output_path.parent().unwrap_or_else(|| Path::new("."));
    std::fs::create_dir_all(parent)
        .map_err(|e| format!("Could not create output directory for '{output}': {e}"))?;
    let index = build_ucsc_rmsk_interval_index(
        resource_path,
        &snapshot.assembly_database,
        &snapshot.records,
    )?;
    let mut text = serde_json::to_string_pretty(&index)
        .map_err(|e| format!("Could not serialize UCSC rmsk interval index: {e}"))?;
    text.push('\n');
    std::fs::write(output_path, text)
        .map_err(|e| format!("Could not write UCSC rmsk interval index '{output}': {e}"))?;
    Ok(UcscRmskIndexSummary {
        output: output.to_string(),
        source_resource: resource_path.to_string(),
        assembly_database: index.assembly_database,
        row_count: index.row_count,
        chromosome_count: index.chromosome_count,
    })
}

pub fn project_interval_to_anchor(
    interval: &UcscRmskIndexedInterval,
    anchor: &UcscRmskProjectionAnchor,
) -> Option<UcscRmskProjectedOverlap> {
    if !interval.chromosome.eq_ignore_ascii_case(&anchor.chromosome)
        || anchor.end_1based < anchor.start_1based
    {
        return None;
    }
    let anchor_start_0based = anchor.start_1based.checked_sub(1)?;
    let anchor_end_0based_exclusive = anchor.end_1based;
    let genomic_start_0based = interval.start_0based.max(anchor_start_0based);
    let genomic_end_0based_exclusive = interval
        .end_0based_exclusive
        .min(anchor_end_0based_exclusive);
    if genomic_end_0based_exclusive <= genomic_start_0based {
        return None;
    }
    let anchor_is_reverse = anchor.strand == Some('-');
    let (local_start_0based, local_end_0based_exclusive, local_strand) = if anchor_is_reverse {
        (
            anchor_end_0based_exclusive.checked_sub(genomic_end_0based_exclusive)?,
            anchor_end_0based_exclusive.checked_sub(genomic_start_0based)?,
            flip_strand(&interval.strand),
        )
    } else {
        (
            genomic_start_0based.checked_sub(anchor_start_0based)?,
            genomic_end_0based_exclusive.checked_sub(anchor_start_0based)?,
            normalize_rmsk_strand(&interval.strand),
        )
    };
    Some(UcscRmskProjectedOverlap {
        interval: interval.clone(),
        local_start_0based,
        local_end_0based_exclusive,
        local_strand,
        genomic_start_0based,
        genomic_end_0based_exclusive,
        overlap_bp: genomic_end_0based_exclusive - genomic_start_0based,
        clipped: genomic_start_0based != interval.start_0based
            || genomic_end_0based_exclusive != interval.end_0based_exclusive,
    })
}

pub fn project_index_overlaps_to_anchor(
    index: &UcscRmskIntervalIndex,
    anchor: &UcscRmskProjectionAnchor,
    start_0based: Option<usize>,
    end_0based_exclusive: Option<usize>,
    limit: Option<usize>,
) -> Vec<UcscRmskProjectedOverlap> {
    let anchor_start_0based = anchor.start_1based.saturating_sub(1);
    let anchor_end_0based_exclusive = anchor.end_1based;
    let (query_start, query_end) = if anchor.strand == Some('-') {
        let local_start = start_0based.unwrap_or(0);
        let local_end = end_0based_exclusive
            .unwrap_or(anchor_end_0based_exclusive.saturating_sub(anchor_start_0based));
        (
            anchor_end_0based_exclusive.saturating_sub(local_end),
            anchor_end_0based_exclusive.saturating_sub(local_start),
        )
    } else {
        (
            anchor_start_0based.saturating_add(start_0based.unwrap_or(0)),
            end_0based_exclusive
                .map(|end| anchor_start_0based.saturating_add(end))
                .unwrap_or(anchor_end_0based_exclusive),
        )
    };
    let effective_limit = limit.unwrap_or(usize::MAX);
    index
        .overlapping_intervals(&anchor.chromosome, query_start, query_end)
        .into_iter()
        .filter_map(|interval| {
            let mut overlap = project_interval_to_anchor(interval, anchor)?;
            let clipped_local_start = start_0based
                .map(|start| overlap.local_start_0based.max(start))
                .unwrap_or(overlap.local_start_0based);
            let clipped_local_end = end_0based_exclusive
                .map(|end| overlap.local_end_0based_exclusive.min(end))
                .unwrap_or(overlap.local_end_0based_exclusive);
            if clipped_local_end <= clipped_local_start {
                return None;
            }
            if clipped_local_start != overlap.local_start_0based
                || clipped_local_end != overlap.local_end_0based_exclusive
            {
                if anchor.strand == Some('-') {
                    overlap.genomic_start_0based =
                        anchor_end_0based_exclusive.saturating_sub(clipped_local_end);
                    overlap.genomic_end_0based_exclusive =
                        anchor_end_0based_exclusive.saturating_sub(clipped_local_start);
                } else {
                    overlap.genomic_start_0based =
                        anchor_start_0based.saturating_add(clipped_local_start);
                    overlap.genomic_end_0based_exclusive =
                        anchor_start_0based.saturating_add(clipped_local_end);
                }
                overlap.local_start_0based = clipped_local_start;
                overlap.local_end_0based_exclusive = clipped_local_end;
                overlap.overlap_bp = clipped_local_end - clipped_local_start;
                overlap.clipped = true;
            }
            Some(overlap)
        })
        .take(effective_limit)
        .collect()
}

pub fn write_ucsc_rmsk_resource_snapshot<R: BufRead>(
    reader: R,
    source: &str,
    assembly_database: &str,
    output: &str,
    max_records: Option<usize>,
    fetched_at_unix_ms: u128,
) -> Result<UcscRmskSyncSummary, String> {
    let output_path = Path::new(output);
    let parent = output_path.parent().unwrap_or_else(|| Path::new("."));
    std::fs::create_dir_all(parent)
        .map_err(|e| format!("Could not create output directory for '{output}': {e}"))?;
    let mut temp = tempfile::NamedTempFile::new_in(parent)
        .map_err(|e| format!("Could not create temporary output for '{output}': {e}"))?;
    {
        let mut writer = std::io::BufWriter::new(&mut temp);
        let descriptor = ucsc_rmsk_descriptor(assembly_database);
        write!(writer, "{{\n  \"schema\": ")
            .map_err(|e| format!("Could not write UCSC rmsk output '{output}': {e}"))?;
        serde_json::to_writer(&mut writer, UCSC_RMSK_RESOURCE_SCHEMA)
            .map_err(|e| format!("Could not serialize UCSC rmsk schema: {e}"))?;
        write!(writer, ",\n  \"source\": ")
            .map_err(|e| format!("Could not write UCSC rmsk output '{output}': {e}"))?;
        serde_json::to_writer(&mut writer, source)
            .map_err(|e| format!("Could not serialize UCSC rmsk source: {e}"))?;
        write!(writer, ",\n  \"provider\": ")
            .map_err(|e| format!("Could not write UCSC rmsk output '{output}': {e}"))?;
        serde_json::to_writer(&mut writer, &descriptor.provider)
            .map_err(|e| format!("Could not serialize UCSC rmsk provider: {e}"))?;
        write!(writer, ",\n  \"assembly_database\": ")
            .map_err(|e| format!("Could not write UCSC rmsk output '{output}': {e}"))?;
        serde_json::to_writer(&mut writer, assembly_database)
            .map_err(|e| format!("Could not serialize UCSC rmsk assembly: {e}"))?;
        write!(
            writer,
            ",\n  \"table\": \"{}\",\n  \"coordinate_system\": ",
            UCSC_RMSK_TABLE
        )
        .map_err(|e| format!("Could not write UCSC rmsk output '{output}': {e}"))?;
        serde_json::to_writer(&mut writer, &descriptor.coordinate_system)
            .map_err(|e| format!("Could not serialize UCSC rmsk coordinates: {e}"))?;
        write!(writer, ",\n  \"download_url\": ")
            .map_err(|e| format!("Could not write UCSC rmsk output '{output}': {e}"))?;
        serde_json::to_writer(&mut writer, &descriptor.download_url)
            .map_err(|e| format!("Could not serialize UCSC rmsk URL: {e}"))?;
        write!(writer, ",\n  \"schema_url\": ")
            .map_err(|e| format!("Could not write UCSC rmsk output '{output}': {e}"))?;
        serde_json::to_writer(&mut writer, &descriptor.schema_url)
            .map_err(|e| format!("Could not serialize UCSC rmsk schema URL: {e}"))?;
        write!(writer, ",\n  \"track_url\": ")
            .map_err(|e| format!("Could not write UCSC rmsk output '{output}': {e}"))?;
        serde_json::to_writer(&mut writer, &descriptor.track_url)
            .map_err(|e| format!("Could not serialize UCSC rmsk track URL: {e}"))?;
        write!(
            writer,
            ",\n  \"fetched_at_unix_ms\": {fetched_at_unix_ms},\n  \"field_order\": "
        )
        .map_err(|e| format!("Could not write UCSC rmsk output '{output}': {e}"))?;
        serde_json::to_writer(&mut writer, &descriptor.field_order)
            .map_err(|e| format!("Could not serialize UCSC rmsk field order: {e}"))?;
        write!(writer, ",\n  \"fields\": ")
            .map_err(|e| format!("Could not write UCSC rmsk output '{output}': {e}"))?;
        serde_json::to_writer(&mut writer, &descriptor.fields)
            .map_err(|e| format!("Could not serialize UCSC rmsk fields: {e}"))?;
        write!(writer, ",\n  \"index_recommendations\": ")
            .map_err(|e| format!("Could not write UCSC rmsk output '{output}': {e}"))?;
        serde_json::to_writer(&mut writer, &descriptor.index_recommendations)
            .map_err(|e| format!("Could not serialize UCSC rmsk index recommendations: {e}"))?;
        write!(writer, ",\n  \"records\": [\n")
            .map_err(|e| format!("Could not write UCSC rmsk output '{output}': {e}"))?;

        let mut row_count = 0usize;
        let mut truncated = false;
        for (idx, line) in reader.lines().enumerate() {
            let line = line.map_err(|e| {
                format!(
                    "Could not read UCSC rmsk input '{}' at line {}: {e}",
                    source,
                    idx + 1
                )
            })?;
            let Some(record) = parse_ucsc_rmsk_record(&line, idx + 1)? else {
                continue;
            };
            if let Some(limit) = max_records
                && row_count >= limit
            {
                truncated = true;
                break;
            }
            if row_count > 0 {
                write!(writer, ",\n")
                    .map_err(|e| format!("Could not write UCSC rmsk output '{output}': {e}"))?;
            }
            write!(writer, "    ")
                .map_err(|e| format!("Could not write UCSC rmsk output '{output}': {e}"))?;
            serde_json::to_writer(&mut writer, &record)
                .map_err(|e| format!("Could not serialize UCSC rmsk record: {e}"))?;
            row_count += 1;
        }
        let mut warnings = Vec::<String>::new();
        if truncated {
            warnings.push(format!(
                "Snapshot truncated to {} UCSC rmsk row(s) by --limit; do not use this as a complete assembly resource.",
                row_count
            ));
        }
        write!(
            writer,
            "\n  ],\n  \"row_count\": {row_count},\n  \"truncated\": {truncated},\n  \"warnings\": "
        )
        .map_err(|e| format!("Could not write UCSC rmsk output '{output}': {e}"))?;
        serde_json::to_writer(&mut writer, &warnings)
            .map_err(|e| format!("Could not serialize UCSC rmsk warnings: {e}"))?;
        write!(writer, "\n}}\n")
            .map_err(|e| format!("Could not write UCSC rmsk output '{output}': {e}"))?;
        writer
            .flush()
            .map_err(|e| format!("Could not flush UCSC rmsk output '{output}': {e}"))?;

        drop(writer);
        temp.persist(output_path).map_err(|e| {
            format!(
                "Could not move temporary UCSC rmsk output into '{}': {}",
                output_path.display(),
                e.error
            )
        })?;
        Ok(UcscRmskSyncSummary {
            output: output.to_string(),
            source: source.to_string(),
            assembly_database: assembly_database.to_string(),
            row_count,
            truncated,
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use tempfile::tempdir;

    #[test]
    fn parses_ucsc_rmsk_records_with_header_skip() {
        let text = include_str!("../test_files/fixtures/resources/ucsc.rmsk.hg38.edge.txt");
        let records = parse_ucsc_rmsk_records(text).expect("parse rmsk fixture");
        assert_eq!(records.len(), 4);
        assert_eq!(records[0].geno_name, "chr1");
        assert_eq!(records[0].geno_start, 10000);
        assert_eq!(records[0].geno_end, 10468);
        assert_eq!(records[0].rep_name, "(TAACCC)n");
        assert_eq!(records[0].rep_class, "Simple_repeat");
        assert_eq!(records[2].strand, "-");
        assert_eq!(records[2].rep_class, "LINE");
        assert_eq!(records[3].rep_family, "hAT-Charlie");
    }

    #[test]
    fn descriptor_carries_ucsc_urls_and_index_recommendations() {
        let descriptor = ucsc_rmsk_descriptor("hg38");
        assert_eq!(descriptor.schema, UCSC_RMSK_DESCRIPTOR_SCHEMA);
        assert!(
            descriptor
                .download_url
                .ends_with("/hg38/database/rmsk.txt.gz")
        );
        assert!(descriptor.schema_url.ends_with("/hg38/database/rmsk.sql"));
        assert!(
            descriptor
                .index_recommendations
                .iter()
                .any(|row| row.index_id == "interval_by_chrom_bin")
        );
    }

    #[test]
    fn writes_bounded_resource_snapshot() {
        let text = include_str!("../test_files/fixtures/resources/ucsc.rmsk.hg38.edge.txt");
        let td = tempdir().expect("tempdir");
        let output = td.path().join("rmsk.json");
        let summary = write_ucsc_rmsk_resource_snapshot(
            std::io::Cursor::new(text),
            "fixture",
            "hg38",
            output.to_string_lossy().as_ref(),
            Some(2),
            123,
        )
        .expect("write snapshot");
        assert_eq!(summary.row_count, 2);
        assert!(summary.truncated);
        let json: serde_json::Value =
            serde_json::from_str(&std::fs::read_to_string(output).expect("read snapshot"))
                .expect("parse json");
        assert_eq!(json["schema"], UCSC_RMSK_RESOURCE_SCHEMA);
        assert_eq!(json["assembly_database"], "hg38");
        assert_eq!(json["row_count"], 2);
        assert_eq!(json["records"].as_array().map(|rows| rows.len()), Some(2));
        assert!(
            json["index_recommendations"]
                .as_array()
                .is_some_and(|rows| !rows.is_empty())
        );
    }

    #[test]
    fn interval_index_queries_overlapping_chromosome_rows() {
        let text = include_str!("../test_files/fixtures/resources/ucsc.rmsk.hg38.edge.txt");
        let records = parse_ucsc_rmsk_records(text).expect("parse rmsk fixture");
        let index =
            build_ucsc_rmsk_interval_index("fixture", "hg38", &records).expect("build index");
        assert_eq!(index.schema, UCSC_RMSK_INTERVAL_INDEX_SCHEMA);
        assert_eq!(index.row_count, 4);
        assert!(index.repeat_name_index.contains_key("(TAACCC)n"));
        let rows = index.overlapping_intervals("chr1", 10_300, 10_800);
        assert_eq!(rows.len(), 2);
        assert_eq!(rows[0].rep_name, "(TAACCC)n");
        assert_eq!(rows[1].rep_class, "Satellite");
    }

    #[test]
    fn projects_plus_anchor_repeat_overlap_to_local_feature_span() {
        let text = include_str!("../test_files/fixtures/resources/ucsc.rmsk.hg38.edge.txt");
        let records = parse_ucsc_rmsk_records(text).expect("parse rmsk fixture");
        let index =
            build_ucsc_rmsk_interval_index("fixture", "hg38", &records).expect("build index");
        let anchor = UcscRmskProjectionAnchor {
            chromosome: "chr1".to_string(),
            start_1based: 10_001,
            end_1based: 11_000,
            strand: Some('+'),
        };
        let overlaps = project_index_overlaps_to_anchor(&index, &anchor, None, None, None);
        assert_eq!(overlaps.len(), 2);
        assert_eq!(overlaps[0].local_start_0based, 0);
        assert_eq!(overlaps[0].local_end_0based_exclusive, 468);
        assert_eq!(overlaps[0].local_strand, "+");
        assert!(!overlaps[0].clipped);
        assert_eq!(overlaps[1].local_start_0based, 468);
        assert_eq!(overlaps[1].local_end_0based_exclusive, 1000);
        assert_eq!(overlaps[1].local_strand, "-");
        assert!(overlaps[1].clipped);
    }

    #[test]
    fn projects_reverse_anchor_repeat_overlap_by_flipping_local_strand() {
        let text = include_str!("../test_files/fixtures/resources/ucsc.rmsk.hg38.edge.txt");
        let records = parse_ucsc_rmsk_records(text).expect("parse rmsk fixture");
        let index =
            build_ucsc_rmsk_interval_index("fixture", "hg38", &records).expect("build index");
        let anchor = UcscRmskProjectionAnchor {
            chromosome: "chr1".to_string(),
            start_1based: 10_001,
            end_1based: 11_000,
            strand: Some('-'),
        };
        let overlaps = project_index_overlaps_to_anchor(&index, &anchor, None, None, None);
        assert_eq!(overlaps.len(), 2);
        assert_eq!(overlaps[0].local_start_0based, 532);
        assert_eq!(overlaps[0].local_end_0based_exclusive, 1000);
        assert_eq!(overlaps[0].local_strand, "-");
        assert_eq!(overlaps[1].local_start_0based, 0);
        assert_eq!(overlaps[1].local_end_0based_exclusive, 532);
        assert_eq!(overlaps[1].local_strand, "+");
        assert!(overlaps[1].clipped);
    }

    #[test]
    fn projects_partial_overlap_and_marks_clipping() {
        let text = include_str!("../test_files/fixtures/resources/ucsc.rmsk.hg38.edge.txt");
        let records = parse_ucsc_rmsk_records(text).expect("parse rmsk fixture");
        let index =
            build_ucsc_rmsk_interval_index("fixture", "hg38", &records).expect("build index");
        let anchor = UcscRmskProjectionAnchor {
            chromosome: "chr1".to_string(),
            start_1based: 10_201,
            end_1based: 10_400,
            strand: Some('+'),
        };
        let overlaps = project_index_overlaps_to_anchor(&index, &anchor, None, None, None);
        assert_eq!(overlaps.len(), 1);
        assert_eq!(overlaps[0].local_start_0based, 0);
        assert_eq!(overlaps[0].local_end_0based_exclusive, 200);
        assert_eq!(overlaps[0].genomic_start_0based, 10_200);
        assert_eq!(overlaps[0].genomic_end_0based_exclusive, 10_400);
        assert_eq!(overlaps[0].overlap_bp, 200);
        assert!(overlaps[0].clipped);
    }

    #[test]
    fn writes_interval_index_from_resource_snapshot() {
        let text = include_str!("../test_files/fixtures/resources/ucsc.rmsk.hg38.edge.txt");
        let td = tempdir().expect("tempdir");
        let resource = td.path().join("rmsk.resource.json");
        let index_path = td.path().join("rmsk.index.json");
        write_ucsc_rmsk_resource_snapshot(
            std::io::Cursor::new(text),
            "fixture",
            "hg38",
            resource.to_string_lossy().as_ref(),
            None,
            123,
        )
        .expect("write snapshot");
        let summary = write_ucsc_rmsk_interval_index_from_resource(
            resource.to_string_lossy().as_ref(),
            index_path.to_string_lossy().as_ref(),
        )
        .expect("write index");
        assert_eq!(summary.row_count, 4);
        let index =
            read_ucsc_rmsk_interval_index(index_path.to_string_lossy().as_ref()).expect("read");
        assert_eq!(index.chromosome_count, 1);
        assert_eq!(index.overlapping_intervals("chr1", 10_000, 10_100).len(), 1);
    }
}
