//! Inline synthetic original reads; literal expectations are independent of SAM SEQ orientation.

use super::*;

fn meter() -> Meter<'static> {
    Meter {
        limits: Limits::default(),
        work: 0,
        control: None,
    }
}
fn observation(strand: char) -> pileup::Observation {
    pileup::Observation {
        observation_id: "synthetic_observation".into(),
        row_ordinal: 1,
        blocks: vec![SupportedBlock {
            query: (12, 112),
            target: (100, 200),
        }],
        connectors: vec![],
        query_to_reference_strand: strand,
    }
}
fn original() -> String {
    format!("{}{}{}", "A".repeat(12), "C".repeat(100), "T".repeat(18))
}
fn scan(sequence: &str) -> Vec<Run> {
    literals(sequence, Policy::default(), 1000).unwrap()
}
fn variant() -> Variant {
    Variant {
        read_id: "synthetic_read".into(),
        sequence_sha256: format!("sha256:{}", "a".repeat(64)),
        query_length: 130,
        runs: scan(&original()),
    }
}
fn unit() -> ObservedSamQueryUnitKey {
    ObservedSamQueryUnitKey {
        source_sha256: format!("sha256:{}", "b".repeat(64)),
        unit_id: "synthetic_unit".into(),
    }
}

#[test]
fn read_end_oriented_literals_keep_12a_18t_without_counting_internal_runs_as_terminal() {
    let runs = scan(&original());
    assert_eq!(
        runs.iter()
            .map(|r| (r.base, r.length, r.query_interval))
            .collect::<Vec<_>>(),
        vec![('A', 12, (0, 12)), ('T', 18, (112, 130))]
    );
    assert_eq!(runs[0].read_ends, vec![End::FivePrime]);
    assert_eq!(runs[1].read_ends, vec![End::ThreePrime]);
    let internal = scan(&format!(
        "{}{}{}",
        "C".repeat(70),
        "A".repeat(12),
        "C".repeat(70)
    ));
    assert_eq!(internal.len(), 1);
    assert!(!internal[0].terminal && internal[0].read_ends.is_empty());
    let short = scan(&"T".repeat(18));
    assert_eq!(short.len(), 1);
    assert_eq!(short[0].read_ends, vec![End::FivePrime, End::ThreePrime]);
    assert!(scan(&"A".repeat(9)).is_empty());
    assert!(literals(&original(), Policy::default(), 1).is_err());
}

#[test]
fn read_end_oriented_strand_matrix_never_complements_literal_bases_or_query_coordinates() {
    for transcript in ['+', '-'] {
        for query in ['+', '-'] {
            let result = project(
                "sample",
                &unit(),
                transcript,
                &[observation(query)],
                &[variant()],
                None,
                Policy::default(),
                &mut meter(),
            )
            .unwrap();
            assert!(
                result.direction
                    == if transcript == query {
                        Direction::Same
                    } else {
                        Direction::Opposite
                    }
            );
            let runs = &result.variants[0].runs;
            assert_eq!(
                (runs[0].base, runs[0].length, runs[0].query_interval),
                ('A', 12, (0, 12))
            );
            assert_eq!(
                (runs[1].base, runs[1].length, runs[1].query_interval),
                ('T', 18, (112, 130))
            );
            assert_eq!(
                runs[0].transcript_ends,
                Some(vec![if transcript == query {
                    End::FivePrime
                } else {
                    End::ThreePrime
                }])
            );
            assert!(
                runs.iter()
                    .all(|r| !r.aligned_to_reference && r.reference_anchor.is_some())
            );
        }
    }
}

#[test]
fn read_end_oriented_alignment_keeps_softclips_and_partially_aligned_runs() {
    let mut o = observation('+');
    o.blocks = vec![SupportedBlock {
        query: (6, 130),
        target: (100, 224),
    }];
    let e = project(
        "sample",
        &unit(),
        '+',
        &[o],
        &[variant()],
        None,
        Policy::default(),
        &mut meter(),
    )
    .unwrap();
    assert!(!e.variants[0].runs[0].aligned_to_reference);
    assert_eq!(e.variants[0].runs[0].aligned_query_intervals, vec![(6, 12)]);
    assert_eq!(e.variants[0].runs[0].reference_intervals, vec![(100, 106)]);
    assert!(e.variants[0].runs[1].aligned_to_reference);
    let mut o = observation('-');
    o.blocks = vec![SupportedBlock {
        query: (0, 130),
        target: (100, 230),
    }];
    let e = project(
        "sample",
        &unit(),
        '-',
        &[o],
        &[variant()],
        None,
        Policy::default(),
        &mut meter(),
    )
    .unwrap();
    assert_eq!(e.variants[0].runs[0].reference_intervals, vec![(218, 230)]);
    assert_eq!(e.variants[0].runs[1].reference_intervals, vec![(100, 118)]);
    assert!(e.variants[0].runs.iter().all(|r| r.aligned_to_reference));
}

#[test]
fn read_end_oriented_mixed_mapping_retains_alternatives_and_unavailable_is_not_zero() {
    let e = project(
        "sample",
        &unit(),
        '+',
        &[observation('+'), observation('-')],
        &[variant()],
        None,
        Policy::default(),
        &mut meter(),
    )
    .unwrap();
    assert!(e.direction == Direction::Alternatives);
    assert!(
        e.variants[0]
            .runs
            .iter()
            .all(|r| r.transcript_ends.is_none())
    );
    let empty = project(
        "sample",
        &unit(),
        '+',
        &[observation('+')],
        &[],
        None,
        Policy::default(),
        &mut meter(),
    )
    .unwrap();
    assert_eq!(empty.availability, "no_bound_basecalls");
    assert_eq!(empty.read_key, e.read_key);
}

fn document() -> (ObservedSamAnalysis, Document) {
    let (_, mut analysis) = super::super::tests::fixture("plus");
    analysis.units[0].unclipped_axis_length = Some(130);
    analysis.units[0].rows = vec![serde_json::from_value(serde_json::json!({
        "ordinal":1,"qname":"synthetic_read","read_group":null,"flag":16,"target":"synthetic_contig",
        "position_1based":101,"reported_mapq":60,"cigar":[["S",18],["M",100],["S",12]],
        "sequence_length":130,"quality_length":null,"axis_length":130,"clipping":null,"unavailable":null
    })).unwrap()];
    let doc = Document {
        schema: BASES_SCHEMA.into(),
        analysis_context_sha256: analysis.context_sha256.clone(),
        orientation: "original_sequencing_5_to_3".into(),
        provenance: "Hand-crafted original read, not reverse-complemented SAM SEQ".into(),
        reads: vec![OriginalRead {
            unit: analysis.units[0].unit.clone(),
            read_id: "synthetic_read".into(),
            sequence: original(),
        }],
    };
    (analysis, doc)
}

#[test]
fn read_end_oriented_import_refuses_misbindings_duplicates_and_budget_overflow() {
    let (analysis, doc) = document();
    let hash = format!("sha256:{}", "e".repeat(64));
    let check = |doc: &Document| {
        prepare(
            Some(Input {
                document: doc,
                file_sha256: &hash,
            }),
            &analysis,
            Policy::default(),
            &mut meter(),
        )
    };
    assert_eq!(check(&doc).unwrap().len(), 1);
    for field in [
        "orientation",
        "analysis_context_sha256",
        "sequence",
        "read_id",
        "unit",
    ] {
        let mut bad = doc.clone();
        match field {
            "orientation" => bad.orientation = "sam_seq".into(),
            "analysis_context_sha256" => bad.analysis_context_sha256.clear(),
            "sequence" => {
                bad.reads[0].sequence.pop();
            }
            "read_id" => bad.reads[0].read_id = "another_read".into(),
            _ => bad.reads[0].unit.unit_id = "another_unit".into(),
        }
        assert!(check(&bad).is_err(), "{field}");
    }
    let mut duplicate = doc.clone();
    duplicate.reads.push(doc.reads[0].clone());
    assert!(check(&duplicate).is_err());
    let mut tiny = meter();
    tiny.limits.max_work = 129;
    assert!(
        prepare(
            Some(Input {
                document: &doc,
                file_sha256: &hash
            }),
            &analysis,
            Policy::default(),
            &mut tiny
        )
        .is_err()
    );
    assert!(
        prepare(
            None,
            &analysis,
            Policy {
                window_bp: 0,
                minimum_run_bp: 10
            },
            &mut meter()
        )
        .is_err()
    );
}

#[test]
fn read_end_oriented_conflicting_original_sequences_are_deterministic_alternatives() {
    let (analysis, mut doc) = document();
    let mut second = doc.reads[0].clone();
    second.sequence = original().replace('A', "T");
    doc.reads.push(second);
    let hash = format!("sha256:{}", "e".repeat(64));
    let run = |doc: &Document| {
        prepare(
            Some(Input {
                document: doc,
                file_sha256: &hash,
            }),
            &analysis,
            Policy::default(),
            &mut meter(),
        )
        .unwrap()
    };
    let a = run(&doc);
    doc.reads.reverse();
    let b = run(&doc);
    assert_eq!(
        serde_json::to_vec(&a[&analysis.units[0].unit]).unwrap(),
        serde_json::to_vec(&b[&analysis.units[0].unit]).unwrap()
    );
    let e = project(
        "sample",
        &analysis.units[0].unit,
        '+',
        &[observation('+')],
        &a[&analysis.units[0].unit],
        None,
        Policy::default(),
        &mut meter(),
    )
    .unwrap();
    assert_eq!(e.availability, "alternatives");
    assert_eq!(e.variants.len(), 2);
}

#[test]
fn read_end_oriented_hardclipped_sam_cannot_substitute_for_full_original_bases() {
    let (mut analysis, mut doc) = document();
    let row = &mut analysis.units[0].rows[0];
    row.sequence_length = Some(100);
    row.cigar = Some(vec![('H', 18), ('M', 100), ('H', 12)]);
    let hash = format!("sha256:{}", "e".repeat(64));
    assert!(
        prepare(
            Some(Input {
                document: &doc,
                file_sha256: &hash
            }),
            &analysis,
            Policy::default(),
            &mut meter()
        )
        .is_ok()
    );
    doc.reads[0].sequence = "C".repeat(100);
    assert!(
        prepare(
            Some(Input {
                document: &doc,
                file_sha256: &hash
            }),
            &analysis,
            Policy::default(),
            &mut meter()
        )
        .is_err()
    );
}
