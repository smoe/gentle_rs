//! Inline hand-crafted SAM 1.6 grammar/resource fixtures, not aligner output or real reads.
//! Recreate through the `sam_text::tests` filter; no fixtures, files, tools or networks are opened.

use super::*;

fn dictionary() -> BTreeMap<String, usize> {
    BTreeMap::from([("ctg_a".into(), 1000), ("ctg_b".into(), 1000)])
}

fn limits() -> Limits {
    Limits {
        max_input_bytes: 1_000_000,
        max_line_bytes: 65536,
        max_lines: 100,
        max_rows: 100,
        max_dictionary_entries: 10,
        max_header_ids: 10,
        max_fields_per_line: 64,
        max_operations: 10000,
        max_sa_entries: 100,
        max_retained_bytes: 4 * 1024 * 1024,
        max_work: 1_000_000,
    }
}

const HEADER: &str = "@HD\tVN:1.6\tSO:unsorted\n@SQ\tSN:ctg_a\tLN:1000\n@SQ\tSN:ctg_b\tLN:1000\n";
const MAPPED: &str = "same_name\t0\tctg_a\t101\t0\t10M10N10M10N10M\t*\t0\t0\t*\t*";

fn input() -> String {
    format!(
        "{HEADER}@RG\tID:library_a\tDS:synthetisch\u{e4}\n@PG\tID:producer\tCL:unused \u{3a9} command\n@CO\tUnicode \u{3a9} comment\n{MAPPED}\tRG:Z:library_a\tSA:Z:ctg_b,201,-,10M15N10M,255,0;\tNM:i:3\tAS:i:-12\n\
        same_name\t3856\tctg_b\t201\t255\t10M15N10M\t*\t0\t0\t*\t*\n\
        unmapped\t4\t*\t0\t255\t*\t*\t0\t0\t*\t*\n"
    )
}

fn decode(
    raw: &[u8],
    d: &BTreeMap<String, usize>,
    b: &Limits,
    check: &mut dyn FnMut() -> Result<(), EngineError>,
) -> Result<Document, EngineError> {
    let mut consumer = Consumer::new(d, b, check)?;
    for line in raw.split_inclusive(|&b| b == b'\n') {
        consumer.push_line(line, check)?;
    }
    consumer.finish(check)
}

fn run(raw: &str) -> Document {
    decode(raw.as_bytes(), &dictionary(), &limits(), &mut || Ok(())).unwrap()
}
fn bad_row(row: &str) -> EngineError {
    decode(
        format!("{HEADER}{row}\n").as_bytes(),
        &dictionary(),
        &limits(),
        &mut || Ok(()),
    )
    .unwrap_err()
}

#[test]
fn sam_text_lf_crlf_hashes_ordinals_and_nonexclusive_flag_counts() {
    let raw = input();
    let a = run(&raw);
    let crlf = raw.replace('\n', "\r\n");
    let b = run(&crlf);
    assert_eq!(a.grammar, "sam_text_strict.private.v1");
    assert_eq!(
        a.raw_sha256,
        crate::digest_utils::sha256_prefixed_bytes(raw.as_bytes())
    );
    assert_eq!(
        b.raw_sha256,
        crate::digest_utils::sha256_prefixed_bytes(crlf.as_bytes())
    );
    assert_ne!(a.raw_sha256, b.raw_sha256);
    assert_eq!(a.raw_bytes, raw.len());
    assert_eq!(a.counts, b.counts);
    assert_eq!(
        a.counts,
        Counts {
            lines: 9,
            header_lines: 6,
            alignment_rows: 3,
            unmapped: 1,
            secondary: 1,
            supplementary: 1,
            qc_fail: 1,
            duplicate: 1
        }
    );
    assert_eq!(a.dictionary_order, vec!["ctg_a", "ctg_b"]);
    assert_eq!(
        a.rows.iter().map(|r| r.row.ordinal).collect::<Vec<_>>(),
        vec![1, 2, 3]
    );
    // Equal QNAME is retained twice, not used to group or fabricate read identity.
    assert_eq!(a.rows[0].row.qname, a.rows[1].row.qname);
    assert_eq!(a.rows.len(), 3);
    assert_eq!(a.rows[0].read_group.as_deref(), Some("library_a"));
    assert_eq!(a.rows[1].row.mapq, 255);
    assert_eq!(a.rows[0].row.mapq, 0);
    assert_eq!(
        serde_json::to_value(a.rows.iter().map(|r| &r.row).collect::<Vec<_>>()).unwrap(),
        serde_json::to_value(b.rows.iter().map(|r| &r.row).collect::<Vec<_>>()).unwrap()
    );
    let no_newline = run(raw.trim_end_matches('\n'));
    assert_eq!(no_newline.counts, a.counts);
    assert_ne!(no_newline.raw_sha256, a.raw_sha256);
}

#[test]
fn sam_text_consumes_only_rg_sa_cg_without_retaining_sequence_scores_or_commands() {
    let seq = "ACGTTGCAACGTAGCTTGCATGCAACGTAC";
    assert_eq!(seq.len(), 30);
    let quality = "J".repeat(30);
    let raw = input().replacen("\t*\t*\tRG:", &format!("\t{seq}\t{quality}\tRG:"), 1);
    let a = run(&raw);
    assert_eq!(a.rows[0].row.sequence_length, Some(30));
    assert_eq!(a.rows[0].row.quality_length, Some(30));
    assert_eq!(
        a.rows[0].companions,
        vec![Companion {
            target: "ctg_b".into(),
            position_1based: 201,
            strand: '-',
            cigar: vec![('M', 10), ('N', 15), ('M', 10)],
            mapq: 255
        }]
    );
    let debug = format!("{a:?}");
    for absent in [
        seq,
        quality.as_str(),
        "unused \u{3a9} command",
        "AS:i",
        "NM:i",
        "Unicode \u{3a9} comment",
    ] {
        assert!(!debug.contains(absent), "{absent}");
    }
    assert_eq!(a.rows.len(), 3, "SA is not another observed row");
    let cg = run(&format!("{HEADER}{MAPPED}\tCG:B:I,160,163\n"));
    assert!(cg.rows[0].row.has_cg);
    assert_eq!(cg.operations, 7); // Five actual CIGAR ops plus two encoded array elements.
}

#[test]
fn sam_text_headers_are_ordered_unique_full_dictionary_and_exact_names() {
    for raw in [
        MAPPED.to_string(),
        format!("@SQ\tSN:ctg_a\tLN:1000\n{MAPPED}\n"),
        format!("@SQ\tSN:ctg_a\tLN:1000\n@SQ\tSN:ctg_a\tLN:1000\n"),
        HEADER.replace("LN:1000", "LN:999"),
        HEADER.replace("ctg_b", "alias_b"),
        HEADER.replace("SN:ctg_b", "SN:ctg_b\tSN:ctg_b"),
        HEADER.replace("VN:1.6", "VN:1.5"),
        format!("{HEADER}@HD\tVN:1.6\n"),
        format!("{HEADER}{MAPPED}\n@CO\tlate header\n"),
        format!("{HEADER}@RG\tID:x\n@RG\tID:x\n"),
        format!("{HEADER}@PG\tID:x\n@PG\tID:x\n"),
        format!("{HEADER}@RG\tDS:missing id\n"),
        HEADER.replace("SN:ctg_a", "SN:ctg_a\tTP:circular"),
        format!("{HEADER}@XX\tID:x\n"),
    ] {
        assert!(
            decode(raw.as_bytes(), &dictionary(), &limits(), &mut || Ok(())).is_err(),
            "{raw}"
        );
    }
    let mut reversed = HEADER.replace(
        "@SQ\tSN:ctg_a\tLN:1000\n@SQ\tSN:ctg_b\tLN:1000",
        "@SQ\tSN:ctg_b\tLN:1000\n@SQ\tSN:ctg_a\tLN:1000",
    );
    reversed.push_str(MAPPED);
    assert_eq!(run(&reversed).dictionary_order, vec!["ctg_b", "ctg_a"]);
    let header_only = run(HEADER);
    assert!(header_only.rows.is_empty());
    // Declaration strings/aliases are not a FASTA check or alternate contig resolver.
    assert!(
        bad_row(&MAPPED.replace("ctg_a", "alias_a"))
            .message
            .contains("target")
    );
}

#[test]
fn sam_text_mandatory_fields_ranges_clips_and_lengths_fail_closed() {
    let fields = MAPPED.split('\t').collect::<Vec<_>>();
    for (index, value) in [
        (0, ""),
        (0, "bad name"),
        (0, "read@name"),
        (1, "4096"),
        (1, "-1"),
        (1, "65536"),
        (2, "missing"),
        (3, "0"),
        (3, "2147483648"),
        (3, "+101"),
        (4, "256"),
        (5, "0M"),
        (5, "1M1S1M"),
        (5, "3M2H3M"),
        (5, "3M2"),
        (5, "3Z"),
        (6, "other"),
        (7, "2147483648"),
        (8, "-2147483648"),
        (9, "AC?T"),
        (10, "\u{e9}"),
    ] {
        let mut f = fields.clone();
        f[index] = value;
        bad_row(&f.join("\t"));
    }
    bad_row(&format!("{MAPPED}\textra"));
    bad_row(&fields[..10].join("\t"));
    let mut f = fields.clone();
    f[9] = "ACGT";
    assert!(bad_row(&f.join("\t")).message.contains("sequence_quality"));
    let mut f = fields;
    f[10] = "IIII";
    assert!(bad_row(&f.join("\t")).message.contains("sequence_quality"));
    let unavailable = run(&format!("{HEADER}*\t4\t*\t0\t255\t*\t*\t0\t0\t*\t*\n"));
    assert_eq!(unavailable.rows[0].row.qname, "*");
    assert!(unavailable.rows[0].row.sequence_length.is_none());
}

#[test]
fn sam_text_optional_types_validate_all_values_but_do_not_inventory_ignored_tags() {
    for value in [
        "aa:A:!",
        "ab:i:-2147483648",
        "ac:i:+4294967295",
        "ad:f:-1.25e2",
        "ae:Z: spaces : preserved only raw",
        "af:H:00ABFF",
        "ag:B:c,-128,127",
        "ah:B:C,0,255",
        "ai:B:s,-32768,32767",
        "aj:B:S,0,65535",
        "ak:B:i,-2147483648,2147483647",
        "al:B:I,4294967295",
        "am:B:f,0,-1.5,2e3",
        "an:B:I",
        "ao:Z:",
        "ap:H:",
    ] {
        assert_eq!(run(&format!("{HEADER}{MAPPED}\t{value}\n")).rows.len(), 1);
    }
    for value in [
        "aa:A: ",
        "aa:A:xx",
        "aa:i:4294967296",
        "aa:i:--1",
        "aa:f:1.",
        "aa:f:nan",
        "aa:f:1e99",
        "aa:H:abc",
        "aa:H:ab",
        "aa:B:C,-1",
        "aa:B:c,128",
        "aa:B:I,4294967296",
        "aa:B:x,1",
        "aa:B:f,NaN",
        "aa:B:I,",
        "a:Z:x",
        "1a:Z:x",
        "aa:Q:x",
        "aa:Z:x\taa:Z:y",
        "RG:i:1",
        "RG:Z:undeclared",
        "SA:i:1",
        "SA:Z:",
        "SA:Z:ctg_b,201,+,10M,0,0",
        "SA:Z:unknown,201,+,10M,0,0;",
        "SA:Z:ctg_b,201,?,10M,0,0;",
        "SA:Z:ctg_b,0,+,10M,0,0;",
        "SA:Z:ctg_b,999,+,10M,0,0;",
        "SA:Z:ctg_b,201,+,*,0,0;",
        "SA:Z:ctg_b,201,+,10M,0,-1;",
        "CG:Z:ignored",
        "CG:B:f,1",
        "CG:B:I,-1",
        "CG:B:I,1\tCG:B:I,1",
    ] {
        bad_row(&format!("{MAPPED}\t{value}"));
    }
}

#[test]
fn sam_text_valid_but_unsupported_geometry_rows_remain_accounted_without_inference() {
    for (flag, cigar) in [
        (1, "10M"),
        (0, "*"),
        (0, "10D"),
        (0, "20S10N"),
        (0, "2000001S1M"),
    ] {
        let raw = format!("{HEADER}q\t{flag}\tctg_a\t101\t255\t{cigar}\t*\t0\t0\t*\t*\n");
        let d = run(&raw);
        assert_eq!(d.rows.len(), 1);
        assert_eq!(d.rows[0].row.flag, flag);
        assert!(d.rows[0].row.sequence_length.is_none());
        assert!(d.rows[0].row.quality_length.is_none());
        // Text decoding is not unit admission or positive splice/search evidence.
        assert!(d.rows[0].companions.is_empty());
    }
}

#[test]
fn sam_text_budget_exact_boundaries_and_one_over_are_atomic() {
    let raw = input();
    let d = dictionary();
    let b = limits();
    let a = decode(raw.as_bytes(), &d, &b, &mut || Ok(())).unwrap();
    let max_line = raw.split_inclusive('\n').map(str::len).max().unwrap();
    let max_fields = raw
        .lines()
        .filter(|s| !s.starts_with("@CO"))
        .map(|s| s.split('\t').count())
        .max()
        .unwrap();
    for field in 0..11 {
        let mut exact = b.clone();
        match field {
            0 => exact.max_input_bytes = raw.len(),
            1 => exact.max_line_bytes = max_line,
            2 => exact.max_lines = a.counts.lines,
            3 => exact.max_rows = a.rows.len(),
            4 => exact.max_dictionary_entries = d.len(),
            5 => exact.max_header_ids = 2,
            6 => exact.max_fields_per_line = max_fields,
            7 => exact.max_operations = a.operations,
            8 => exact.max_sa_entries = a.sa_entries,
            9 => exact.max_retained_bytes = a.retained_peak_bytes,
            _ => exact.max_work = a.work,
        }
        assert!(
            decode(raw.as_bytes(), &d, &exact, &mut || Ok(())).is_ok(),
            "exact {field}"
        );
        match field {
            0 => exact.max_input_bytes -= 1,
            1 => exact.max_line_bytes -= 1,
            2 => exact.max_lines -= 1,
            3 => exact.max_rows -= 1,
            4 => exact.max_dictionary_entries -= 1,
            5 => exact.max_header_ids -= 1,
            6 => exact.max_fields_per_line -= 1,
            7 => exact.max_operations -= 1,
            8 => exact.max_sa_entries -= 1,
            9 => exact.max_retained_bytes -= 1,
            _ => exact.max_work -= 1,
        }
        assert!(
            decode(raw.as_bytes(), &d, &exact, &mut || Ok(())).is_err(),
            "one over {field}"
        );
    }
}

#[test]
fn sam_text_limits_admission_raw_framing_utf8_and_poison_are_not_partial_success() {
    let d = dictionary();
    let b = limits();
    for bytes in [
        b"\n".as_slice(),
        b"\xef\xbb\xbf@HD\tVN:1.6\n",
        b"\xff\n",
        b"@CO\tbad\rcontrol\n",
        b"@CO\tfirst\n@CO\tsecond\n",
    ] {
        let mut p = Consumer::new(&d, &b, &mut || Ok(())).unwrap();
        assert!(p.push_line(bytes, &mut || Ok(())).is_err());
        assert!(
            p.finish(&mut || Ok(()))
                .unwrap_err()
                .message
                .contains("poisoned")
        );
    }
    let mut p = Consumer::new(&d, &b, &mut || Ok(())).unwrap();
    p.push_line(b"@CO\tunterminated", &mut || Ok(())).unwrap();
    assert!(p.push_line(b"@CO\tmore\n", &mut || Ok(())).is_err());
    assert!(p.finish(&mut || Ok(())).is_err());
    let mut too_large = b.clone();
    too_large.max_line_bytes = 1024 * 1024 + 1;
    assert!(Consumer::new(&d, &too_large, &mut || Ok(())).is_err());
    let invalid = BTreeMap::from([("ctg,alias".into(), 1000)]);
    assert!(Consumer::new(&invalid, &b, &mut || Ok(())).is_err());
}

#[test]
fn sam_text_cancellation_propagates_at_every_poll_and_finish_has_no_prefix_result() {
    let raw = input();
    let d = dictionary();
    let b = limits();
    let mut polls = 0;
    decode(raw.as_bytes(), &d, &b, &mut || {
        polls += 1;
        Ok(())
    })
    .unwrap();
    assert!(polls > 10);
    for stop in 0..polls {
        let mut seen = 0;
        let e = decode(raw.as_bytes(), &d, &b, &mut || {
            let at = seen;
            seen += 1;
            if at == stop {
                Err(fail("synthetic_cancel"))
            } else {
                Ok(())
            }
        })
        .unwrap_err();
        assert!(e.message.contains("synthetic_cancel"), "{stop}: {e:?}");
    }
    let mut p = Consumer::new(&d, &b, &mut || Ok(())).unwrap();
    assert!(
        p.push_line(HEADER.lines().next().unwrap().as_bytes(), &mut || Err(
            fail("synthetic_cancel")
        ))
        .is_err()
    );
    assert!(
        p.finish(&mut || Ok(()))
            .unwrap_err()
            .message
            .contains("poisoned")
    );
}
