//! NCBI GenBank XML (`GBSet/GBSeq`) parser and dialect detection helpers.
//!
//! This module intentionally supports only `GBSet/GBSeq` in the current pass.
//! Other XML dialects (notably `INSDSet/INSDSeq`) are detected and rejected
//! with explicit diagnostics so import behavior stays deterministic.

use anyhow::{Result, anyhow};
use gb_io::seq::{Feature, FeatureKind, Location, Seq, Topology};
use serde::Deserialize;

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum NcbiXmlDialect {
    GbSetGbSeq,
    InsdSetInsdSeq,
    Unknown,
}

impl NcbiXmlDialect {
    pub fn label(self) -> &'static str {
        match self {
            Self::GbSetGbSeq => "GBSet/GBSeq",
            Self::InsdSetInsdSeq => "INSDSet/INSDSeq",
            Self::Unknown => "unknown",
        }
    }
}

pub fn detect_ncbi_xml_dialect(input: &str) -> NcbiXmlDialect {
    let lower = input.to_ascii_lowercase();
    if lower.contains("<gbset") {
        NcbiXmlDialect::GbSetGbSeq
    } else if lower.contains("<insdset") {
        NcbiXmlDialect::InsdSetInsdSeq
    } else {
        NcbiXmlDialect::Unknown
    }
}

pub fn parse_gbseq_xml_file(path: &str) -> Result<Vec<Seq>> {
    let text = std::fs::read_to_string(path)
        .map_err(|e| anyhow!("Could not read XML file '{path}': {e}"))?;
    parse_gbseq_xml_text(&text).map_err(|e| anyhow!("Could not parse XML file '{path}': {e}"))
}

pub fn parse_gbseq_xml_text(xml: &str) -> Result<Vec<Seq>> {
    match detect_ncbi_xml_dialect(xml) {
        NcbiXmlDialect::GbSetGbSeq => {}
        NcbiXmlDialect::InsdSetInsdSeq => {
            return Err(anyhow!(
                "Unsupported XML dialect '{}'; only GBSet/GBSeq is currently supported",
                NcbiXmlDialect::InsdSetInsdSeq.label()
            ));
        }
        NcbiXmlDialect::Unknown => {
            return Err(anyhow!(
                "Unsupported XML dialect: expected '{}' root element",
                NcbiXmlDialect::GbSetGbSeq.label()
            ));
        }
    }

    let parsed: GbSetXml =
        quick_xml::de::from_str(xml).map_err(|e| anyhow!("Malformed GBSet XML: {e}"))?;
    if parsed.sequences.is_empty() {
        return Err(anyhow!("Malformed GBSet XML: no GBSeq records found"));
    }

    parsed
        .sequences
        .iter()
        .enumerate()
        .map(|(record_idx, record)| gbseq_record_to_seq(record, record_idx))
        .collect()
}

#[derive(Debug, Deserialize)]
#[serde(rename = "GBSet")]
struct GbSetXml {
    #[serde(rename = "GBSeq", default)]
    sequences: Vec<GbSeqXml>,
}

#[derive(Debug, Deserialize)]
struct GbSeqXml {
    #[serde(rename = "GBSeq_locus")]
    locus: Option<String>,
    #[serde(rename = "GBSeq_moltype")]
    moltype: Option<String>,
    #[serde(rename = "GBSeq_topology")]
    topology: Option<String>,
    #[serde(rename = "GBSeq_division")]
    division: Option<String>,
    #[serde(rename = "GBSeq_definition")]
    definition: Option<String>,
    #[serde(rename = "GBSeq_primary-accession")]
    primary_accession: Option<String>,
    #[serde(rename = "GBSeq_accession-version")]
    accession_version: Option<String>,
    #[serde(rename = "GBSeq_source")]
    _source: Option<String>,
    #[serde(rename = "GBSeq_sequence")]
    sequence: Option<String>,
    #[serde(rename = "GBSeq_feature-table")]
    feature_table: Option<GbFeatureTableXml>,
}

#[derive(Debug, Deserialize)]
struct GbFeatureTableXml {
    #[serde(rename = "GBFeature", default)]
    features: Vec<GbFeatureXml>,
}

#[derive(Debug, Deserialize)]
struct GbFeatureXml {
    #[serde(rename = "GBFeature_key")]
    key: Option<String>,
    #[serde(rename = "GBFeature_location")]
    location: Option<String>,
    #[serde(rename = "GBFeature_intervals")]
    intervals: Option<GbFeatureIntervalsXml>,
    #[serde(rename = "GBFeature_quals")]
    qualifiers: Option<GbFeatureQualsXml>,
}

#[derive(Debug, Deserialize)]
struct GbFeatureIntervalsXml {
    #[serde(rename = "GBInterval", default)]
    intervals: Vec<GbIntervalXml>,
}

#[derive(Debug, Deserialize)]
struct GbIntervalXml {
    #[serde(rename = "GBInterval_from")]
    from: Option<usize>,
    #[serde(rename = "GBInterval_to")]
    to: Option<usize>,
    #[serde(rename = "GBInterval_point")]
    point: Option<usize>,
    #[serde(rename = "GBInterval_iscomp")]
    iscomp: Option<bool>,
}

#[derive(Debug, Deserialize)]
struct GbFeatureQualsXml {
    #[serde(rename = "GBQualifier", default)]
    qualifiers: Vec<GbQualifierXml>,
}

#[derive(Debug, Deserialize)]
struct GbQualifierXml {
    #[serde(rename = "GBQualifier_name")]
    name: Option<String>,
    #[serde(rename = "GBQualifier_value")]
    value: Option<String>,
}

fn gbseq_record_to_seq(record: &GbSeqXml, record_idx: usize) -> Result<Seq> {
    let mut seq = Seq::empty();

    let sequence_text: String = record
        .sequence
        .as_deref()
        .unwrap_or_default()
        .chars()
        .filter(|ch| ch.is_ascii_alphabetic())
        .map(|ch| ch.to_ascii_uppercase())
        .collect();
    if sequence_text.is_empty() {
        return Err(anyhow!(
            "GBSeq record {} is missing sequence data",
            record_idx + 1
        ));
    }

    let accession = nonempty_owned(record.primary_accession.as_deref())
        .or_else(|| accession_from_accession_version(record.accession_version.as_deref()));
    let version = nonempty_owned(record.accession_version.as_deref());

    seq.name = nonempty_owned(record.locus.as_deref())
        .or_else(|| accession.clone())
        .or_else(|| Some(format!("record_{}", record_idx + 1)));
    seq.topology = parse_topology(record.topology.as_deref());
    seq.molecule_type = nonempty_owned(record.moltype.as_deref());
    seq.division = nonempty_owned(record.division.as_deref()).unwrap_or_default();
    seq.definition = nonempty_owned(record.definition.as_deref());
    seq.accession = accession;
    seq.version = version;
    seq.source = None;
    seq.seq = sequence_text.into_bytes();
    seq.len = Some(seq.seq.len());
    seq.features = parse_features(record, seq.name.as_deref().unwrap_or("<unnamed>"))?;

    Ok(seq)
}

fn parse_features(record: &GbSeqXml, seq_label: &str) -> Result<Vec<Feature>> {
    let mut features = Vec::new();
    let Some(feature_table) = &record.feature_table else {
        return Ok(features);
    };
    for (feature_idx, raw_feature) in feature_table.features.iter().enumerate() {
        let feature_key = nonempty_owned(raw_feature.key.as_deref())
            .unwrap_or_else(|| "misc_feature".to_string());
        let feature_kind = FeatureKind::from(feature_key.as_str());
        let location_text = resolve_feature_location_text(raw_feature).ok_or_else(|| {
            anyhow!(
                "GBSeq '{}' feature #{} ('{}') is missing location",
                seq_label,
                feature_idx + 1,
                feature_key
            )
        })?;
        let parsed_location = Location::from_gb_format(location_text.as_str()).map_err(|e| {
            anyhow!(
                "GBSeq '{}' feature #{} ('{}') has invalid location '{}': {}",
                seq_label,
                feature_idx + 1,
                feature_key,
                location_text,
                e
            )
        })?;
        let location = canonicalize_location(parsed_location).map_err(|e| {
            anyhow!(
                "GBSeq '{}' feature #{} ('{}') could not canonicalize location '{}': {}",
                seq_label,
                feature_idx + 1,
                feature_key,
                location_text,
                e
            )
        })?;
        let qualifiers = parse_qualifiers(raw_feature);
        features.push(Feature {
            kind: feature_kind,
            location,
            qualifiers,
        });
    }
    Ok(features)
}

fn resolve_feature_location_text(feature: &GbFeatureXml) -> Option<String> {
    if let Some(location) = nonempty_owned(feature.location.as_deref()) {
        return Some(location);
    }
    location_from_intervals(
        feature
            .intervals
            .as_ref()
            .map(|raw| raw.intervals.as_slice())
            .unwrap_or_default(),
    )
}

fn location_from_intervals(intervals: &[GbIntervalXml]) -> Option<String> {
    if intervals.is_empty() {
        return None;
    }
    let mut parts: Vec<String> = vec![];
    for interval in intervals {
        let from = interval.from.or(interval.point)?;
        let to = interval.to.or(interval.point).unwrap_or(from);
        if from == 0 || to == 0 {
            return None;
        }
        let (start, end) = if from <= to { (from, to) } else { (to, from) };
        let mut piece = format!("{start}..{end}");
        if interval.iscomp.unwrap_or(false) {
            piece = format!("complement({piece})");
        }
        parts.push(piece);
    }
    if parts.len() == 1 {
        parts.into_iter().next()
    } else {
        Some(format!("join({})", parts.join(",")))
    }
}

fn parse_qualifiers(feature: &GbFeatureXml) -> Vec<(gb_io::seq::QualifierKey, Option<String>)> {
    let mut qualifiers = vec![];
    let Some(raw_quals) = feature.qualifiers.as_ref() else {
        return qualifiers;
    };
    for qualifier in &raw_quals.qualifiers {
        let Some(name) = nonempty_owned(qualifier.name.as_deref()) else {
            continue;
        };
        let value = nonempty_owned(qualifier.value.as_deref());
        qualifiers.push((name.into(), value));
    }
    qualifiers
}

fn parse_topology(raw: Option<&str>) -> Topology {
    if raw
        .map(|value| value.trim().eq_ignore_ascii_case("circular"))
        .unwrap_or(false)
    {
        Topology::Circular
    } else {
        Topology::Linear
    }
}

fn accession_from_accession_version(raw: Option<&str>) -> Option<String> {
    let version = nonempty_owned(raw)?;
    Some(
        version
            .split_once('.')
            .map(|(accession, _)| accession.trim().to_string())
            .filter(|accession| !accession.is_empty())
            .unwrap_or(version),
    )
}

fn nonempty_owned(raw: Option<&str>) -> Option<String> {
    let text = raw.unwrap_or_default().trim();
    (!text.is_empty()).then_some(text.to_string())
}

fn canonicalize_location(location: Location) -> Result<Location> {
    let value = serde_json::to_value(&location)?;
    Ok(serde_json::from_value(value)?)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_detect_ncbi_xml_dialect_gbset() {
        let xml = r#"<?xml version="1.0"?><GBSet><GBSeq/></GBSet>"#;
        assert_eq!(detect_ncbi_xml_dialect(xml), NcbiXmlDialect::GbSetGbSeq);
    }

    #[test]
    fn test_detect_ncbi_xml_dialect_insdset() {
        let xml = r#"<?xml version="1.0"?><INSDSet><INSDSeq/></INSDSet>"#;
        assert_eq!(detect_ncbi_xml_dialect(xml), NcbiXmlDialect::InsdSetInsdSeq);
    }

    #[test]
    fn test_parse_gbseq_xml_text_toy_small() {
        let text = std::fs::read_to_string("test_files/fixtures/import_parity/toy.small.gbseq.xml")
            .expect("read toy XML");
        let parsed = parse_gbseq_xml_text(&text).expect("parse GBSet/GBSeq");
        assert_eq!(parsed.len(), 1);
        let seq = &parsed[0];
        assert_eq!(seq.seq.len(), 120);
        assert!(
            seq.features
                .iter()
                .any(|feature| feature.kind.to_string().eq_ignore_ascii_case("gene"))
        );
    }

    #[test]
    fn test_parse_gbseq_xml_text_rejects_insdset() {
        let err = parse_gbseq_xml_text("<INSDSet><INSDSeq/></INSDSet>")
            .expect_err("INSDSet should be rejected");
        assert!(
            err.to_string().contains("Unsupported XML dialect"),
            "expected unsupported-dialect error, got: {err}"
        );
    }
}
