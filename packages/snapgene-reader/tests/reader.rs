use snapgene_reader::{PacketType, PrimerRecord, SnapGeneError, parse_bytes, parse_path};

fn packet(tag: u8, payload: &[u8]) -> Vec<u8> {
    let mut out = Vec::with_capacity(1 + 4 + payload.len());
    out.push(tag);
    out.extend_from_slice(&(payload.len() as u32).to_be_bytes());
    out.extend_from_slice(payload);
    out
}

fn cookie() -> Vec<u8> {
    let mut payload = Vec::new();
    payload.extend_from_slice(b"SnapGene");
    payload.extend_from_slice(&1u16.to_be_bytes());
    payload.extend_from_slice(&1u16.to_be_bytes());
    payload.extend_from_slice(&1u16.to_be_bytes());
    packet(0x09, &payload)
}

fn dna(sequence: &str, circular: bool) -> Vec<u8> {
    let mut payload = Vec::with_capacity(sequence.len() + 1);
    payload.push(if circular { 0x01 } else { 0x00 });
    payload.extend_from_slice(sequence.as_bytes());
    packet(0x00, &payload)
}

fn xml_packet(tag: u8, xml: &str) -> Vec<u8> {
    packet(tag, xml.as_bytes())
}

fn assemble(parts: &[Vec<u8>]) -> Vec<u8> {
    let total = parts.iter().map(Vec::len).sum();
    let mut out = Vec::with_capacity(total);
    for part in parts {
        out.extend_from_slice(part);
    }
    out
}

#[test]
fn test_parse_fixture_from_path() {
    let parsed = parse_path("tests/data/toy.small.dna").expect("parse fixture");
    assert_eq!(parsed.dna.sequence.len(), 120);
    assert_eq!(parsed.features.len(), 4);
    assert!(parsed.notes.is_some());
}

#[test]
fn test_invalid_cookie_is_rejected() {
    let mut payload = Vec::new();
    payload.extend_from_slice(b"NotGene!");
    payload.extend_from_slice(&1u16.to_be_bytes());
    payload.extend_from_slice(&1u16.to_be_bytes());
    payload.extend_from_slice(&1u16.to_be_bytes());
    let bytes = assemble(&[packet(0x09, &payload), dna("ACGT", false)]);
    let err = parse_bytes(&bytes).expect_err("invalid cookie should fail");
    assert!(matches!(err, SnapGeneError::InvalidCookie { .. }));
}

#[test]
fn test_truncated_packet_is_rejected() {
    let mut bytes = cookie();
    bytes.push(0x00);
    bytes.extend_from_slice(&10u32.to_be_bytes());
    bytes.extend_from_slice(b"ACGT");
    let err = parse_bytes(&bytes).expect_err("truncated payload should fail");
    assert!(matches!(err, SnapGeneError::UnexpectedEof { .. }));
}

#[test]
fn test_missing_dna_packet_is_rejected() {
    let err = parse_bytes(&cookie()).expect_err("missing DNA packet should fail");
    assert!(matches!(err, SnapGeneError::MissingDnaPacket));
}

#[test]
fn test_unknown_packets_are_retained() {
    let bytes = assemble(&[
        cookie(),
        dna("ACGTACGT", false),
        packet(0x03, b"enzyme-state"),
        packet(0x0D, b"display-state"),
    ]);
    let parsed = parse_bytes(&bytes).expect("parse unknown packets");
    assert_eq!(parsed.extra_packets.len(), 2);
    assert_eq!(parsed.extra_packets[0].packet_type, PacketType::Unknown(0x03));
    assert_eq!(parsed.extra_packets[0].payload, b"enzyme-state");
    assert_eq!(parsed.extra_packets[1].packet_type, PacketType::Unknown(0x0D));
    assert_eq!(parsed.extra_packets[1].payload, b"display-state");
}

#[cfg(feature = "gb-io")]
#[test]
fn test_joined_feature_segments_roundtrip_to_gb_io() {
    let features_xml = r##"<?xml version="1.0"?><Features><Feature name="joined" type="gene" directionality="1"><Segment range="5-10" color="#ffffff" type="standard"/><Segment range="20-30" color="#ffffff" type="standard"/><Q name="gene"><V text="joinedA"/></Q></Feature></Features>"##;
    let bytes = assemble(&[cookie(), dna(&"A".repeat(60), false), xml_packet(0x0A, features_xml)]);
    let parsed = parse_bytes(&bytes).expect("parse joined feature");
    let seq = gb_io::seq::Seq::try_from(&parsed).expect("adapt to gb_io");
    assert_eq!(seq.features.len(), 1);
    assert_eq!(seq.features[0].location.to_gb_format(), "join(5..10,20..30)");
}

#[cfg(feature = "gb-io")]
#[test]
fn test_reverse_joined_feature_segments_roundtrip_to_gb_io() {
    let features_xml = r##"<?xml version="1.0"?><Features><Feature name="joined_reverse" type="gene" directionality="2"><Segment range="5-10" color="#ffffff" type="standard"/><Segment range="20-30" color="#ffffff" type="standard"/><Q name="gene"><V text="joinedRev"/></Q></Feature></Features>"##;
    let bytes = assemble(&[cookie(), dna(&"A".repeat(60), false), xml_packet(0x0A, features_xml)]);
    let parsed = parse_bytes(&bytes).expect("parse reverse joined feature");
    let seq = gb_io::seq::Seq::try_from(&parsed).expect("adapt to gb_io");
    assert_eq!(
        seq.features[0].location.to_gb_format(),
        "complement(join(5..10,20..30))"
    );
}

#[cfg(feature = "gb-io")]
#[test]
fn test_wraparound_feature_segments_roundtrip_to_gb_io() {
    let features_xml = r##"<?xml version="1.0"?><Features><Feature name="wrap" type="misc_feature" directionality="1"><Segment range="110-10" color="#ffffff" type="standard"/></Feature></Features>"##;
    let bytes = assemble(&[cookie(), dna(&"A".repeat(120), true), xml_packet(0x0A, features_xml)]);
    let parsed = parse_bytes(&bytes).expect("parse wrap feature");
    let seq = gb_io::seq::Seq::try_from(&parsed).expect("adapt to gb_io");
    let text = seq.features[0].location.to_gb_format();
    assert!(text.contains("110..120"));
    assert!(text.contains("1..10"));
}

#[test]
fn test_primer_binding_visibility_helper_filters_hidden_sites() {
    let primers_xml = r#"<?xml version="1.0"?><Primers><HybridizationParams minContinuousMatchLen="10" minMeltingTemperature="40" allowMismatch="1" showAdditionalFivePrimeMatches="1"/><Primer name="p1" sequence="ATGCATGCATGC" description="demo"><BindingSite location="5-20" boundStrand="0" simplified="1" annealedBases="AAAAAAAAAAAA" meltingTemperature="55"/><BindingSite location="5-20" boundStrand="0" simplified="1" annealedBases="AAAAAAAAAAAA" meltingTemperature="55"/><BindingSite location="30-40" boundStrand="1" simplified="0" annealedBases="AAAAA" meltingTemperature="30"/></Primer></Primers>"#;
    let bytes = assemble(&[cookie(), dna(&"A".repeat(80), false), xml_packet(0x05, primers_xml)]);
    let parsed = parse_bytes(&bytes).expect("parse primers");
    let primer: &PrimerRecord = parsed.primers.first().expect("one primer");
    let visible = primer.visible_binding_sites(parsed.primer_hybridization_params.as_ref());
    assert_eq!(visible.len(), 1);
    assert_eq!(visible[0].location_start_1based, 5);
    assert_eq!(visible[0].location_end_1based, 20);
}

#[cfg(feature = "gb-io")]
#[test]
fn test_fixture_adapter_parity_core_fields() {
    let parsed = parse_bytes(include_bytes!("data/toy.small.dna")).expect("parse fixture");
    let seq = gb_io::seq::Seq::try_from(&parsed).expect("adapt fixture");
    assert_eq!(seq.seq.len(), 120);
    assert_eq!(seq.topology, gb_io::seq::Topology::Linear);
    assert!(
        seq.features
            .iter()
            .any(|feature| feature.kind.to_string().eq_ignore_ascii_case("gene"))
    );
    assert!(
        seq.features
            .iter()
            .flat_map(|feature| feature.qualifier_values("gene"))
            .any(|value| value == "toyA")
    );
    assert!(
        seq.features
            .iter()
            .flat_map(|feature| feature.qualifier_values("gene"))
            .any(|value| value == "toyB")
    );
}
