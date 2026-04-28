//! UCSC RepeatMasker (`rmsk`) table resource contracts.
//!
//! This module keeps the row parser, portable resource descriptor, and index
//! recommendations together so GUI, CLI, and future genome-track ingestion
//! paths can agree on the same `rmsk` semantics without re-parsing UCSC table
//! details in adapter code.

use serde::{Deserialize, Serialize};
use std::io::{BufRead, Write};
use std::path::Path;

pub const UCSC_RMSK_RESOURCE_SCHEMA: &str = "gentle.ucsc_rmsk_resource.v1";
pub const UCSC_RMSK_DESCRIPTOR_SCHEMA: &str = "gentle.ucsc_rmsk_descriptor.v1";
pub const DEFAULT_UCSC_RMSK_ASSEMBLY: &str = "hg38";
pub const DEFAULT_UCSC_RMSK_RESOURCE_PATH: &str = "data/resources/ucsc.rmsk.hg38.json";
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
}
