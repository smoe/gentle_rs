//! NCBI sequence XML parser and dialect detection helpers.
//!
//! `GBSet/GBSeq` and `INSDSet/INSDSeq` records are normalized into the same
//! `gb_io::seq::Seq` representation used by GenBank/EMBL import, keeping XML
//! support on the shared sequence/feature semantics path.

use anyhow::{Result, anyhow};
use gb_io::seq::{Feature, Location, Seq, Topology};
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
    Ok(parse_gbseq_xml_file_with_dialect(path)?.0)
}

pub fn parse_gbseq_xml_file_with_dialect(path: &str) -> Result<(Vec<Seq>, NcbiXmlDialect)> {
    let text = std::fs::read_to_string(path)
        .map_err(|e| anyhow!("Could not read XML file '{path}': {e}"))?;
    parse_gbseq_xml_text_with_dialect(&text)
        .map_err(|e| anyhow!("Could not parse XML file '{path}': {e}"))
}

pub fn parse_gbseq_xml_text(xml: &str) -> Result<Vec<Seq>> {
    Ok(parse_gbseq_xml_text_with_dialect(xml)?.0)
}

pub fn parse_gbseq_xml_text_with_dialect(xml: &str) -> Result<(Vec<Seq>, NcbiXmlDialect)> {
    let dialect = detect_ncbi_xml_dialect(xml);
    let parsed = match dialect {
        NcbiXmlDialect::GbSetGbSeq => parse_gbset_xml_text(xml)?,
        NcbiXmlDialect::InsdSetInsdSeq => parse_insdset_xml_text(xml)?,
        NcbiXmlDialect::Unknown => {
            return Err(anyhow!(
                "Unsupported XML dialect: expected '{}' or '{}' root element",
                NcbiXmlDialect::GbSetGbSeq.label(),
                NcbiXmlDialect::InsdSetInsdSeq.label()
            ));
        }
    };
    Ok((parsed, dialect))
}

fn parse_gbset_xml_text(xml: &str) -> Result<Vec<Seq>> {
    let parsed: GbSetXml =
        quick_xml::de::from_str(xml).map_err(|e| anyhow!("Malformed GBSet XML: {e}"))?;
    if parsed.sequences.is_empty() {
        return Err(anyhow!("Malformed GBSet XML: no GBSeq records found"));
    }

    parsed
        .sequences
        .into_iter()
        .map(NcbiXmlRecord::from)
        .enumerate()
        .map(|(record_idx, record)| ncbi_xml_record_to_seq(&record, "GBSeq", record_idx))
        .collect()
}

fn parse_insdset_xml_text(xml: &str) -> Result<Vec<Seq>> {
    let parsed: InsdSetXml =
        quick_xml::de::from_str(xml).map_err(|e| anyhow!("Malformed INSDSet XML: {e}"))?;
    if parsed.sequences.is_empty() {
        return Err(anyhow!("Malformed INSDSet XML: no INSDSeq records found"));
    }

    parsed
        .sequences
        .into_iter()
        .map(NcbiXmlRecord::from)
        .enumerate()
        .map(|(record_idx, record)| ncbi_xml_record_to_seq(&record, "INSDSeq", record_idx))
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

#[derive(Debug, Deserialize)]
#[serde(rename = "INSDSet")]
struct InsdSetXml {
    #[serde(rename = "INSDSeq", default)]
    sequences: Vec<InsdSeqXml>,
}

#[derive(Debug, Deserialize)]
struct InsdSeqXml {
    #[serde(rename = "INSDSeq_locus")]
    locus: Option<String>,
    #[serde(rename = "INSDSeq_moltype")]
    moltype: Option<String>,
    #[serde(rename = "INSDSeq_topology")]
    topology: Option<String>,
    #[serde(rename = "INSDSeq_division")]
    division: Option<String>,
    #[serde(rename = "INSDSeq_definition")]
    definition: Option<String>,
    #[serde(rename = "INSDSeq_primary-accession")]
    primary_accession: Option<String>,
    #[serde(rename = "INSDSeq_accession-version")]
    accession_version: Option<String>,
    #[serde(rename = "INSDSeq_source")]
    _source: Option<String>,
    #[serde(rename = "INSDSeq_sequence")]
    sequence: Option<String>,
    #[serde(rename = "INSDSeq_feature-table")]
    feature_table: Option<InsdFeatureTableXml>,
}

#[derive(Debug, Deserialize)]
struct InsdFeatureTableXml {
    #[serde(rename = "INSDFeature", default)]
    features: Vec<InsdFeatureXml>,
}

#[derive(Debug, Deserialize)]
struct InsdFeatureXml {
    #[serde(rename = "INSDFeature_key")]
    key: Option<String>,
    #[serde(rename = "INSDFeature_location")]
    location: Option<String>,
    #[serde(rename = "INSDFeature_intervals")]
    intervals: Option<InsdFeatureIntervalsXml>,
    #[serde(rename = "INSDFeature_quals")]
    qualifiers: Option<InsdFeatureQualsXml>,
}

#[derive(Debug, Deserialize)]
struct InsdFeatureIntervalsXml {
    #[serde(rename = "INSDInterval", default)]
    intervals: Vec<InsdIntervalXml>,
}

#[derive(Debug, Deserialize)]
struct InsdIntervalXml {
    #[serde(rename = "INSDInterval_from")]
    from: Option<usize>,
    #[serde(rename = "INSDInterval_to")]
    to: Option<usize>,
    #[serde(rename = "INSDInterval_point")]
    point: Option<usize>,
    #[serde(rename = "INSDInterval_iscomp")]
    iscomp: Option<bool>,
}

#[derive(Debug, Deserialize)]
struct InsdFeatureQualsXml {
    #[serde(rename = "INSDQualifier", default)]
    qualifiers: Vec<InsdQualifierXml>,
}

#[derive(Debug, Deserialize)]
struct InsdQualifierXml {
    #[serde(rename = "INSDQualifier_name")]
    name: Option<String>,
    #[serde(rename = "INSDQualifier_value")]
    value: Option<String>,
}

#[derive(Debug)]
struct NcbiXmlRecord {
    locus: Option<String>,
    moltype: Option<String>,
    topology: Option<String>,
    division: Option<String>,
    definition: Option<String>,
    primary_accession: Option<String>,
    accession_version: Option<String>,
    sequence: Option<String>,
    features: Vec<NcbiXmlFeature>,
}

#[derive(Debug)]
struct NcbiXmlFeature {
    key: Option<String>,
    location: Option<String>,
    intervals: Vec<NcbiXmlInterval>,
    qualifiers: Vec<NcbiXmlQualifier>,
}

#[derive(Debug)]
struct NcbiXmlInterval {
    from: Option<usize>,
    to: Option<usize>,
    point: Option<usize>,
    iscomp: Option<bool>,
}

#[derive(Debug)]
struct NcbiXmlQualifier {
    name: Option<String>,
    value: Option<String>,
}

impl From<GbSeqXml> for NcbiXmlRecord {
    fn from(record: GbSeqXml) -> Self {
        Self {
            locus: record.locus,
            moltype: record.moltype,
            topology: record.topology,
            division: record.division,
            definition: record.definition,
            primary_accession: record.primary_accession,
            accession_version: record.accession_version,
            sequence: record.sequence,
            features: record
                .feature_table
                .map(|table| {
                    table
                        .features
                        .into_iter()
                        .map(NcbiXmlFeature::from)
                        .collect()
                })
                .unwrap_or_default(),
        }
    }
}

impl From<GbFeatureXml> for NcbiXmlFeature {
    fn from(feature: GbFeatureXml) -> Self {
        Self {
            key: feature.key,
            location: feature.location,
            intervals: feature
                .intervals
                .map(|intervals| {
                    intervals
                        .intervals
                        .into_iter()
                        .map(NcbiXmlInterval::from)
                        .collect()
                })
                .unwrap_or_default(),
            qualifiers: feature
                .qualifiers
                .map(|qualifiers| {
                    qualifiers
                        .qualifiers
                        .into_iter()
                        .map(NcbiXmlQualifier::from)
                        .collect()
                })
                .unwrap_or_default(),
        }
    }
}

impl From<GbIntervalXml> for NcbiXmlInterval {
    fn from(interval: GbIntervalXml) -> Self {
        Self {
            from: interval.from,
            to: interval.to,
            point: interval.point,
            iscomp: interval.iscomp,
        }
    }
}

impl From<GbQualifierXml> for NcbiXmlQualifier {
    fn from(qualifier: GbQualifierXml) -> Self {
        Self {
            name: qualifier.name,
            value: qualifier.value,
        }
    }
}

impl From<InsdSeqXml> for NcbiXmlRecord {
    fn from(record: InsdSeqXml) -> Self {
        Self {
            locus: record.locus,
            moltype: record.moltype,
            topology: record.topology,
            division: record.division,
            definition: record.definition,
            primary_accession: record.primary_accession,
            accession_version: record.accession_version,
            sequence: record.sequence,
            features: record
                .feature_table
                .map(|table| {
                    table
                        .features
                        .into_iter()
                        .map(NcbiXmlFeature::from)
                        .collect()
                })
                .unwrap_or_default(),
        }
    }
}

impl From<InsdFeatureXml> for NcbiXmlFeature {
    fn from(feature: InsdFeatureXml) -> Self {
        Self {
            key: feature.key,
            location: feature.location,
            intervals: feature
                .intervals
                .map(|intervals| {
                    intervals
                        .intervals
                        .into_iter()
                        .map(NcbiXmlInterval::from)
                        .collect()
                })
                .unwrap_or_default(),
            qualifiers: feature
                .qualifiers
                .map(|qualifiers| {
                    qualifiers
                        .qualifiers
                        .into_iter()
                        .map(NcbiXmlQualifier::from)
                        .collect()
                })
                .unwrap_or_default(),
        }
    }
}

impl From<InsdIntervalXml> for NcbiXmlInterval {
    fn from(interval: InsdIntervalXml) -> Self {
        Self {
            from: interval.from,
            to: interval.to,
            point: interval.point,
            iscomp: interval.iscomp,
        }
    }
}

impl From<InsdQualifierXml> for NcbiXmlQualifier {
    fn from(qualifier: InsdQualifierXml) -> Self {
        Self {
            name: qualifier.name,
            value: qualifier.value,
        }
    }
}

fn ncbi_xml_record_to_seq(
    record: &NcbiXmlRecord,
    record_label: &str,
    record_idx: usize,
) -> Result<Seq> {
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
            "{} record {} is missing sequence data",
            record_label,
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
    seq.features = parse_features(
        record,
        record_label,
        seq.name.as_deref().unwrap_or("<unnamed>"),
    )?;

    Ok(seq)
}

fn parse_features(
    record: &NcbiXmlRecord,
    record_label: &str,
    seq_label: &str,
) -> Result<Vec<Feature>> {
    let mut features = Vec::new();
    for (feature_idx, raw_feature) in record.features.iter().enumerate() {
        let feature_key = nonempty_owned(raw_feature.key.as_deref())
            .unwrap_or_else(|| "misc_feature".to_string());
        let feature_kind = feature_key.clone().into();
        let location_text = resolve_feature_location_text(raw_feature).ok_or_else(|| {
            anyhow!(
                "{} '{}' feature #{} ('{}') is missing location",
                record_label,
                seq_label,
                feature_idx + 1,
                feature_key
            )
        })?;
        let parsed_location = Location::from_gb_format(location_text.as_str()).map_err(|e| {
            anyhow!(
                "{} '{}' feature #{} ('{}') has invalid location '{}': {}",
                record_label,
                seq_label,
                feature_idx + 1,
                feature_key,
                location_text,
                e
            )
        })?;
        let location = canonicalize_location(parsed_location).map_err(|e| {
            anyhow!(
                "{} '{}' feature #{} ('{}') could not canonicalize location '{}': {}",
                record_label,
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

fn resolve_feature_location_text(feature: &NcbiXmlFeature) -> Option<String> {
    if let Some(location) = nonempty_owned(feature.location.as_deref()) {
        return Some(location);
    }
    location_from_intervals(feature.intervals.as_slice())
}

fn location_from_intervals(intervals: &[NcbiXmlInterval]) -> Option<String> {
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

fn parse_qualifiers(feature: &NcbiXmlFeature) -> Vec<gb_io::seq::Qualifier> {
    let mut qualifiers = vec![];
    for qualifier in &feature.qualifiers {
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
        let (parsed, dialect) =
            parse_gbseq_xml_text_with_dialect(&text).expect("parse GBSet/GBSeq");
        assert_eq!(dialect, NcbiXmlDialect::GbSetGbSeq);
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
    fn test_parse_gbseq_xml_text_toy_small_insdset() {
        let text =
            std::fs::read_to_string("test_files/fixtures/import_parity/toy.small.insdseq.xml")
                .expect("read toy INSD XML");
        let (parsed, dialect) =
            parse_gbseq_xml_text_with_dialect(&text).expect("parse INSDSet/INSDSeq");
        assert_eq!(dialect, NcbiXmlDialect::InsdSetInsdSeq);
        assert_eq!(parsed.len(), 1);
        let seq = &parsed[0];
        assert_eq!(seq.seq.len(), 120);
        assert!(seq.features.iter().any(|feature| {
            feature.kind.to_string().eq_ignore_ascii_case("gene")
                && feature
                    .qualifier_values("gene")
                    .any(|value| value == "toyB")
        }));
    }

    #[test]
    fn test_parse_gbseq_xml_text_rejects_unknown_xml() {
        let err = parse_gbseq_xml_text("<OtherSet><OtherSeq/></OtherSet>")
            .expect_err("unknown XML dialect should be rejected");
        assert!(
            err.to_string().contains("Unsupported XML dialect"),
            "expected unsupported-dialect error, got: {err}"
        );
    }
}
