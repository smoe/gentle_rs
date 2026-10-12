//! Shared GenBank/EMBL serialization of one annotated nucleotide record.
//! Features use the same INSDC location/qualifier writer. Format-specific header
//! details remain explicit provenance comments, never a second sequence model.

use anyhow::{Result, bail};
use gb_io::seq::{Date, Reference, Seq, Source};
use std::fmt::Write as _;

const DBLINK: &str = "GENtle original GenBank DBLINK: ";
const REFERENCE: &str = "GENtle original GenBank reference description: ";
const SOURCE: &str = "GENtle original GenBank SOURCE: ";

/// Write GenBank with correctly padded topology, retaining the original body.
/// Extended identifiers retain gb-io's extended-name behavior, without truncation.
pub fn genbank_bytes(seq: &Seq) -> Result<Vec<u8>> {
    let mut bytes = Vec::new();
    gb_io::writer::write(&mut bytes, seq)?;
    let newline = bytes
        .iter()
        .position(|b| *b == b'\n')
        .ok_or_else(|| anyhow::anyhow!("GenBank writer omitted LOCUS terminator"))?;
    let name = seq
        .name
        .as_deref()
        .or(seq.accession.as_deref())
        .unwrap_or("UNTITLED")
        .split_whitespace()
        .collect::<Vec<_>>()
        .join("_");
    let length = seq.len().to_string();
    let padding = 28usize.saturating_sub(name.len() + length.len()).max(1);
    // Match the nucleotide model's default. An empty type lets token-based
    // readers consume "circular" as the molecule type and lose topology.
    let molecule = seq.molecule_type.as_deref().unwrap_or("DNA");
    let topology = if seq.is_circular() {
        "circular"
    } else {
        "linear"
    };
    let date = seq
        .date
        .clone()
        .unwrap_or(Date::from_ymd(1970, 1, 1).unwrap());
    let locus = format!(
        "LOCUS       {name}{}{length} bp    {molecule:<7} {topology:<8} {} {date}\n",
        " ".repeat(padding),
        seq.division
    );
    bytes.splice(..=newline, locus.bytes());
    Ok(bytes)
}

fn field(out: &mut String, tag: &str, value: &str) {
    // Metadata has no significant line wrapping; avoid splitting Unicode or words.
    let mut line = String::new();
    for word in value.split_whitespace() {
        if !line.is_empty() && line.len() + word.len() + 1 > 74 {
            writeln!(out, "{tag}   {line}").unwrap();
            line.clear();
        }
        if !line.is_empty() {
            line.push(' ');
        }
        line.push_str(word);
    }
    writeln!(out, "{tag}   {line}").unwrap();
}

/// Write EMBL nucleotide interchange text, without inventing annotations or scores.
/// This is not a submission validator; missing database metadata stays missing.
pub fn embl_bytes(seq: &Seq) -> Result<Vec<u8>> {
    if seq.seq.is_empty()
        || !seq
            .seq
            .iter()
            .all(|b| b"ACGTURYSWKMBDHVN".contains(&b.to_ascii_uppercase()))
    {
        bail!("EMBL nucleotide export requires nonempty IUPAC nucleotide sequence data");
    }
    let name = seq
        .name
        .as_deref()
        .or(seq.accession.as_deref())
        .unwrap_or("UNTITLED");
    if name.is_empty()
        || !name.is_ascii()
        || name.bytes().any(|b| b.is_ascii_whitespace() || b == b';')
    {
        bail!("EMBL record identifier must be nonempty ASCII without whitespace or semicolons");
    }
    let molecule = seq.molecule_type.as_deref().unwrap_or("DNA");
    if !matches!(
        molecule.to_ascii_lowercase().as_str(),
        "dna"
            | "rna"
            | "mrna"
            | "rrna"
            | "trna"
            | "urna"
            | "snrna"
            | "snorna"
            | "ss-dna"
            | "ds-dna"
            | "ss-rna"
            | "ds-rna"
            | "genomic dna"
            | "genomic rna"
    ) {
        bail!("EMBL nucleotide export does not support molecule type '{molecule}'");
    }
    let version = seq
        .version
        .as_deref()
        .and_then(|v| v.rsplit('.').next())
        .filter(|v| !v.is_empty() && v.bytes().all(|b| b.is_ascii_digit()))
        .unwrap_or("0");
    let topology = if seq.is_circular() {
        "circular"
    } else {
        "linear"
    };
    let division = if seq.division == "UNK" {
        "UNC"
    } else {
        &seq.division
    };
    let mut out = format!(
        "ID   {name}; SV {version}; {topology}; {molecule}; STD; {division}; {} BP.\nXX\n",
        seq.seq.len()
    );
    if let Some(value) = &seq.accession {
        field(&mut out, "AC", &format!("{value};"));
    }
    if let Some(value) = &seq.version {
        field(&mut out, "SV", value);
    }
    if let Some(value) = &seq.date {
        field(&mut out, "DT", &format!("{value} (Last updated)"));
    }
    if let Some(value) = &seq.definition {
        field(&mut out, "DE", value);
    }
    if let Some(value) = &seq.keywords {
        field(&mut out, "KW", value);
    }
    if let Some(value) = &seq.source {
        field(&mut out, "OS", &value.source);
        // gb-io flattens ORGANISM and taxonomy into one string. Do not mislabel
        // that string as EMBL's taxonomy-only OC field or guess its boundary.
        out.push_str("XX\n");
        field(
            &mut out,
            "CC",
            &format!("{SOURCE}{}", serde_json::to_string(value)?),
        );
    }
    for (i, reference) in seq.references.iter().enumerate() {
        writeln!(out, "XX\nRN   [{}]", i + 1).unwrap();
        // The original reference description may include nonstandard ranges.
        field(
            &mut out,
            "RC",
            &format!(
                "{REFERENCE}{}",
                serde_json::to_string(&reference.description)?
            ),
        );
        if let Some(value) = &reference.authors {
            field(&mut out, "RA", &format!("{value};"));
        }
        if let Some(value) = &reference.consortium {
            field(&mut out, "RG", &format!("{value};"));
        }
        field(
            &mut out,
            "RT",
            &format!("\"{}\";", reference.title.replace('"', "\"\"")),
        );
        if let Some(value) = &reference.journal {
            field(&mut out, "RL", value);
        }
        if let Some(value) = &reference.pubmed {
            field(&mut out, "RX", &format!("PUBMED; {value}."));
        }
        if let Some(value) = &reference.remark {
            out.push_str("XX\n");
            field(&mut out, "RC", value);
        }
    }
    if let Some(value) = &seq.dblink {
        out.push_str("XX\n");
        field(
            &mut out,
            "CC",
            &format!("{DBLINK}{}", serde_json::to_string(value)?),
        );
    }
    for comment in &seq.comments {
        out.push_str("XX\n");
        field(&mut out, "CC", comment);
    }
    out.push_str("XX\nFH   Key             Location/Qualifiers\nFH\n");
    let gb = String::from_utf8(genbank_bytes(seq)?)?;
    let mut features = false;
    for line in gb.lines() {
        if line.starts_with("FEATURES") {
            features = true;
            continue;
        }
        if features && !line.starts_with("     ") {
            break;
        }
        if features {
            writeln!(out, "FT   {}", &line[5..]).unwrap();
        }
    }
    if let Some(contig) = &seq.contig {
        field(&mut out, "CO", &contig.to_gb_format());
    }
    let counts = b"ACGT".map(|base| {
        seq.seq
            .iter()
            .filter(|b| b.to_ascii_uppercase() == base)
            .count()
    });
    writeln!(
        out,
        "XX\nSQ   Sequence {} BP; {} A; {} C; {} G; {} T; {} other;",
        seq.seq.len(),
        counts[0],
        counts[1],
        counts[2],
        counts[3],
        seq.seq.len() - counts.iter().sum::<usize>()
    )
    .unwrap();
    for (line, chunk) in seq.seq.chunks(60).enumerate() {
        let bases = chunk
            .chunks(10)
            .map(|c| String::from_utf8_lossy(c).to_ascii_lowercase())
            .collect::<Vec<_>>()
            .join(" ");
        writeln!(out, "     {bases:<65}{:>9}", line * 60 + chunk.len()).unwrap();
    }
    out.push_str("//\n");
    Ok(out.into_bytes())
}

fn append(target: &mut Option<String>, value: &str) {
    if let Some(existing) = target {
        if !existing.is_empty() {
            existing.push(' ');
        }
        existing.push_str(value);
    } else {
        *target = Some(value.to_string());
    }
}

/// Restore the standard EMBL metadata omitted by the legacy feature reader.
pub(crate) fn read_embl_metadata(lines: &[String], seq: &mut Seq) -> Result<()> {
    let mut groups: Vec<(String, String)> = Vec::new();
    for line in lines {
        if line.starts_with("SQ") {
            break;
        }
        let tag = line.get(..2).unwrap_or("");
        if tag == "FT" || tag == "FH" {
            continue;
        }
        let value = line.get(5..).unwrap_or("").trim();
        if let Some((last_tag, last_value)) = groups.last_mut() {
            if last_tag == tag && tag != "RN" && tag != "ID" && tag != "DT" {
                if !last_value.is_empty() {
                    last_value.push(' ');
                }
                last_value.push_str(value);
                continue;
            }
        }
        groups.push((tag.into(), value.into()));
    }
    let mut source = None;
    let mut organism = None;
    let mut original_source = None;
    for (tag, value) in groups {
        match tag.as_str() {
            "ID" => {
                let parts = value.split(';').map(str::trim).collect::<Vec<_>>();
                if parts.len() == 7 && parts[1].starts_with("SV ") {
                    seq.molecule_type = Some(
                        match parts[3] {
                            "genomic DNA" => "DNA",
                            "genomic RNA" => "RNA",
                            other => other,
                        }
                        .into(),
                    );
                    seq.division = if parts[5] == "UNC" { "UNK" } else { parts[5] }.into();
                    if seq.version.is_none() && parts[1] != "SV 0" {
                        if let Some(accession) = &seq.accession {
                            seq.version = Some(format!("{accession}.{}", &parts[1][3..]));
                        }
                    }
                }
            }
            "DT" => {
                if let Some(date) = value.split_whitespace().next() {
                    let pieces = date.split('-').collect::<Vec<_>>();
                    if pieces.len() == 3 {
                        let month = [
                            "JAN", "FEB", "MAR", "APR", "MAY", "JUN", "JUL", "AUG", "SEP", "OCT",
                            "NOV", "DEC",
                        ]
                        .iter()
                        .position(|m| *m == pieces[1]);
                        if let (Ok(day), Some(month), Ok(year)) =
                            (pieces[0].parse(), month, pieces[2].parse())
                        {
                            seq.date = Date::from_ymd(year, month as u32 + 1, day).ok();
                        }
                    }
                }
            }
            "KW" => append(&mut seq.keywords, &value),
            "OS" => append(&mut source, &value),
            "OC" => append(&mut organism, &value),
            "CC" => {
                if let Some(raw) = value.strip_prefix(DBLINK) {
                    seq.dblink = Some(serde_json::from_str(raw)?);
                } else if let Some(raw) = value.strip_prefix(SOURCE) {
                    original_source = Some(serde_json::from_str(raw)?);
                } else {
                    seq.comments.push(value);
                }
            }
            "PR" | "DR" => {
                // Preserve native cross-reference text rather than dropping unknown databases.
                seq.comments.push(format!("EMBL {tag}: {value}"));
            }
            "RN" => seq.references.push(Reference {
                description: value.trim_matches(['[', ']']).into(),
                authors: None,
                consortium: None,
                title: String::new(),
                journal: None,
                pubmed: None,
                remark: None,
            }),
            "RA" | "RG" | "RT" | "RL" | "RX" | "RC" | "RP" => {
                if let Some(reference) = seq.references.last_mut() {
                    match tag.as_str() {
                        "RA" => {
                            reference.authors =
                                Some(value.strip_suffix(';').unwrap_or(&value).into())
                        }
                        "RG" => {
                            reference.consortium =
                                Some(value.strip_suffix(';').unwrap_or(&value).into())
                        }
                        "RT" => {
                            let title = value.strip_suffix(';').unwrap_or(&value);
                            reference.title = title
                                .strip_prefix('"')
                                .and_then(|s| s.strip_suffix('"'))
                                .unwrap_or(title)
                                .replace("\"\"", "\"")
                        }
                        "RL" => reference.journal = Some(value),
                        "RX" if value.starts_with("PUBMED;") => {
                            reference.pubmed = Some(value[7..].trim().trim_end_matches('.').into())
                        }
                        "RC" if value.starts_with(REFERENCE) => {
                            reference.description = serde_json::from_str(&value[REFERENCE.len()..])?
                        }
                        "RP" => {
                            reference.description = format!(
                                "{} (bases {})",
                                reference.description,
                                value.replace('-', " to ")
                            )
                        }
                        _ => append(&mut reference.remark, &value),
                    }
                }
            }
            "CO" => {
                seq.contig = Some(gb_io::seq::Location::from_gb_format(
                    &value.replace(' ', ""),
                )?)
            }
            _ => {}
        }
    }
    if let Some(source) = source {
        let organism = organism.map(|taxonomy| format!("{source} {taxonomy}"));
        seq.source = Some(Source { source, organism });
    }
    if let Some(source) = original_source {
        seq.source = Some(source);
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use gb_io::seq::{Feature, Location, Topology};

    // Synthetic annotated interchange record, recreated here; not biological evidence.
    fn fixture() -> Seq {
        let mut seq = Seq::empty();
        seq.name = Some("SYNTHETIC_LOCAL".into());
        seq.accession = Some("SYNTH001".into());
        seq.version = Some("SYNTH001.3".into());
        seq.seq = b"ACGTN".repeat(24);
        seq.molecule_type = Some("DNA".into());
        seq.date = Some(Date::from_ymd(2026, 9, 13).unwrap());
        seq.division = "SYN".into();
        seq.definition = Some("Hand-crafted serialization fixture".into());
        seq.keywords = Some("synthetic; test.".into());
        seq.source = Some(Source {
            source: "synthetic construct".into(),
            organism: Some("artificial sequences.".into()),
        });
        seq.dblink = Some("BioProject: SYNTHETIC_ONLY".into());
        seq.comments = vec![
            "Evidence class: synthetic; not experimentally validated.".into(),
            "Source-report SHA256: no real report in this hand-crafted fixture.".into(),
        ];
        seq.references.push(Reference {
            description: "1 (bases 1 to 120)".into(),
            authors: Some("Synthetic A.,Synthetic B.".into()),
            consortium: Some("Test only".into()),
            title: "Not a real \"publication\"".into(),
            journal: Some("Unpublished".into()),
            pubmed: None,
            remark: Some("Provenance test, not a citation.".into()),
        });
        seq.features.push(Feature {
            kind: "misc_feature".into(),
            location: Location::from_gb_format("complement(join(<1..10,25..>40))").unwrap(),
            qualifiers: vec![
                (
                    "note".into(),
                    Some("Quoted \"label\"; raw_score=4; source_id=synthetic".into()),
                ),
                ("note".into(), Some("Repeated qualifier preserved".into())),
                ("pseudo".into(), None),
                ("empty".into(), Some(String::new())),
                ("translation".into(), Some("M".repeat(150))),
            ],
        });
        seq
    }

    fn read_embl(bytes: &[u8]) -> Seq {
        crate::dna_sequence::parse_embl_records(std::str::from_utf8(bytes).unwrap())
            .unwrap()
            .remove(0)
    }

    #[test]
    #[ignore = "explicit temporary artifacts for optional independent-reader verification"]
    fn write_annotated_format_parity_interop_example() {
        let output = tempfile::tempdir().unwrap().keep();
        for topology in [Topology::Linear, Topology::Circular] {
            let mut seq = fixture();
            seq.topology = topology;
            let stem = if seq.is_circular() {
                "circular"
            } else {
                "linear"
            };
            std::fs::write(
                output.join(format!("{stem}.gb")),
                genbank_bytes(&seq).unwrap(),
            )
            .unwrap();
            std::fs::write(
                output.join(format!("{stem}.embl")),
                embl_bytes(&seq).unwrap(),
            )
            .unwrap();
        }
        println!("Synthetic interoperability fixtures: {}", output.display());
    }

    #[test]
    fn annotated_format_parity_preserves_features_metadata_and_both_topologies() {
        for topology in [Topology::Linear, Topology::Circular] {
            let mut seq = fixture();
            seq.topology = topology;
            let before = seq.clone();
            let output = embl_bytes(&seq).unwrap();
            assert_eq!(output, embl_bytes(&seq).unwrap());
            let parsed = read_embl(&output);
            assert_eq!(parsed.seq, seq.seq);
            assert_eq!(parsed.topology, seq.topology);
            assert_eq!(parsed.name, seq.name);
            assert_eq!(parsed.accession, seq.accession);
            assert_eq!(parsed.version, seq.version);
            assert_eq!(parsed.molecule_type, seq.molecule_type);
            assert_eq!(parsed.date, seq.date);
            assert_eq!(parsed.division, seq.division);
            assert_eq!(parsed.definition, seq.definition);
            assert_eq!(parsed.source, seq.source);
            assert_eq!(parsed.keywords, seq.keywords);
            assert_eq!(parsed.comments, seq.comments);
            assert_eq!(parsed.dblink, seq.dblink);
            assert_eq!(parsed.references, seq.references);
            let gb = genbank_bytes(&parsed).unwrap();
            let gb = gb_io::reader::SeqReader::new(gb.as_slice())
                .next()
                .unwrap()
                .unwrap();
            assert_eq!(gb.seq.to_ascii_uppercase(), seq.seq);
            for (expected, actual) in seq.features.iter().zip(&parsed.features) {
                assert_eq!(
                    expected.location.to_gb_format(),
                    actual.location.to_gb_format()
                );
                assert_eq!(expected.qualifiers, actual.qualifiers);
            }
            assert_eq!(parsed.features.len(), seq.features.len());
            assert_eq!(gb.features.len(), seq.features.len());
            let unwrap_values = |feature: &Feature| {
                feature
                    .qualifiers
                    .iter()
                    .map(|(key, value)| {
                        (
                            key.clone(),
                            value
                                .as_ref()
                                .map(|text| text.split_whitespace().collect::<Vec<_>>().join(" ")),
                        )
                    })
                    .collect::<Vec<_>>()
            };
            assert_eq!(
                unwrap_values(&gb.features[0]),
                unwrap_values(&seq.features[0])
            );
            assert_eq!(seq, before);
        }
    }

    #[test]
    fn annotated_format_parity_genbank_columns_and_nucleotide_rejections() {
        let mut seq = fixture();
        for topology in [Topology::Linear, Topology::Circular] {
            seq.topology = topology;
            let data = genbank_bytes(&seq).unwrap();
            let locus = std::str::from_utf8(data.split(|b| *b == b'\n').next().unwrap()).unwrap();
            assert_eq!(locus.len(), 79);
            assert_eq!(&locus[41..43], "bp");
            assert_eq!(&locus[47..54], "DNA    ");
            assert_eq!(
                &locus[55..63],
                if seq.is_circular() {
                    "circular"
                } else {
                    "linear  "
                }
            );
            assert_eq!(&locus[64..67], "SYN");
            assert_eq!(&locus[68..79], "13-SEP-2026");
        }
        seq.molecule_type = Some("protein".into());
        assert!(embl_bytes(&seq).is_err());
        seq.molecule_type = Some("RNA".into());
        seq.seq = b"AUGCNN".to_vec();
        assert_eq!(
            read_embl(&embl_bytes(&seq).unwrap())
                .molecule_type
                .as_deref(),
            Some("RNA")
        );
        seq.name = Some("invalid; name".into());
        assert!(embl_bytes(&seq).is_err());
    }

    #[test]
    fn annotated_format_parity_quoted_slash_continuation_is_not_a_new_qualifier() {
        // Hand-crafted wrapped quote, as can be emitted by INSDC feature writers.
        let native = "ID   SYNTH; SV 0; linear; DNA; STD; SYN; 4 BP.\nFT   misc_feature    1..4\nFT                   /note=\"path continues\nFT                   /inside/the/quoted/value\"\nFT                   /pseudo\nSQ   Sequence 4 BP;\n     acgt 4\n//\n";
        let seq = crate::dna_sequence::parse_embl_records(native)
            .unwrap()
            .remove(0);
        assert_eq!(
            seq.features[0].qualifiers,
            vec![
                (
                    "note".into(),
                    Some("path continues /inside/the/quoted/value".into())
                ),
                ("pseudo".into(), None),
            ]
        );
    }

    #[test]
    fn annotated_format_parity_native_embl_metadata_and_multiple_records() {
        let native = "ID   SYNTH; SV 7; circular; DNA; STD; SYN; 4 BP.\nAC   SYNTH;\nDT   01-JAN-2000 (Created)\nDT   02-FEB-2001 (Last updated)\nCC   Evidence stays a comment.\nSQ   Sequence 4 BP;\n     acgt 4\n//\n\n";
        let records = crate::dna_sequence::parse_embl_records(&native.repeat(2)).unwrap();
        assert_eq!(records.len(), 2);
        assert_eq!(records[0].version.as_deref(), Some("SYNTH.7"));
        assert_eq!(records[0].date, Some(Date::from_ymd(2001, 2, 2).unwrap()));
        assert_eq!(records[0].comments, ["Evidence stays a comment."]);
    }
}
