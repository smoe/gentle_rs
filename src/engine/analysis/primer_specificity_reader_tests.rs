//! Synthetic old/new reader equivalence for hit, warning, product and summary data.
//! The fixture is shared with the prepared-reader tests; no external tools run.

use super::*;
use crate::genomes::{PreparedSequenceRegionReader, prepared_sequence_reader_tests::fixture};

fn hits() -> Vec<(&'static str, BlastHit)> {
    let hit = BlastHit {
        subject_id: "1".into(),
        identity_percent: 100.0,
        alignment_length: 8,
        mismatches: 0,
        gap_opens: 0,
        query_start: 1,
        query_end: 8,
        subject_start: 1,
        subject_end: 8,
        evalue: 1e-8,
        bit_score: 20.0,
        query_coverage_percent: Some(100.0),
        query_length: Some(8),
        aligned_query: Some("ACGTACGT".into()),
        aligned_subject: Some("ACGTACGT".into()),
    };
    vec![
        ("ACGTACGT", hit.clone()),
        (
            "GGTTAAAA",
            BlastHit {
                subject_id: "NC_000099.1".into(),
                subject_start: 8,
                subject_end: 1,
                aligned_query: Some("GGTTAAAA".into()),
                aligned_subject: Some("GGTTAAAA".into()),
                ..hit.clone()
            },
        ),
        (
            "ACGTACGT",
            BlastHit {
                query_start: 3,
                query_end: 6,
                subject_start: 11,
                subject_end: 14,
                alignment_length: 4,
                query_length: None,
                aligned_query: None,
                aligned_subject: None,
                ..hit.clone()
            },
        ),
        (
            "ACGTACGT",
            BlastHit {
                subject_end: 9,
                alignment_length: 9,
                gap_opens: 1,
                aligned_query: Some("ACGTACG-T".into()),
                aligned_subject: Some("ACGTACGTA".into()),
                ..hit.clone()
            },
        ),
        (
            "ACGTACGT",
            BlastHit {
                aligned_query: Some("ACGT".into()),
                ..hit.clone()
            },
        ),
        (
            "ACGTACGT",
            BlastHit {
                subject_id: "missing".into(),
                ..hit.clone()
            },
        ),
        (
            "ACGTACGT",
            BlastHit {
                subject_start: 999,
                subject_end: 1006,
                query_length: None,
                aligned_query: None,
                aligned_subject: None,
                ..hit.clone()
            },
        ),
        (
            "ACGTACGT",
            BlastHit {
                query_start: 6,
                query_end: 3,
                subject_start: 14,
                subject_end: 11,
                alignment_length: 4,
                query_length: None,
                aligned_query: None,
                aligned_subject: None,
                ..hit
            },
        ),
    ]
}

fn evaluate(
    reader: &mut PreparedSequenceRegionReader<'_>,
    policy: &PrimerSpecificityPolicy,
) -> (
    Vec<PrimerSpecificityPrimerHit>,
    Vec<PrimerSpecificityPrimerHit>,
    Vec<String>,
) {
    let aliases = BTreeMap::from([("1".into(), "chr1".into()), ("chr1".into(), "chr1".into())]);
    let mut warnings = Vec::new();
    let mut roles = Vec::new();
    for role in [
        PrimerSpecificityPrimerRole::Forward,
        PrimerSpecificityPrimerRole::Reverse,
    ] {
        roles.push(
            hits()
                .iter()
                .enumerate()
                .map(|(index, (query, hit))| {
                    GentleEngine::primer_specificity_hit_from_blast(
                        reader,
                        "ReaderFixture",
                        None,
                        role.clone(),
                        index,
                        query,
                        hit,
                        &aliases,
                        policy,
                        &mut warnings,
                    )
                })
                .collect::<Vec<_>>(),
        );
    }
    let reverse = roles.pop().unwrap();
    let forward = roles.pop().unwrap();
    (forward, reverse, warnings)
}

#[test]
fn primer_specificity_local_reader_preserves_hits_warnings_products_and_summary() {
    for newline in ["\n", "\r\n"] {
        let temp = tempfile::tempdir().unwrap();
        let f = fixture(temp.path(), newline);
        for mode in [
            PrimerSpecificityFullAlignmentMode::Disabled,
            PrimerSpecificityFullAlignmentMode::BestEffort,
            PrimerSpecificityFullAlignmentMode::Required,
        ] {
            let mut policy = PrimerSpecificityPolicy::default();
            policy.full_alignment.mode = mode;
            let mut old = PreparedSequenceRegionReader::uncached(&f.catalog);
            let mut local = PreparedSequenceRegionReader::new(&f.catalog);
            let expected = evaluate(&mut old, &policy);
            let actual = evaluate(&mut local, &policy);
            assert_eq!(
                serde_json::to_value(&actual).unwrap(),
                serde_json::to_value(&expected).unwrap()
            );
            assert_eq!(
                local.index_load_count(),
                1,
                "both directions and fallback windows share the index"
            );
            assert_eq!(actual.0.len(), hits().len(), "no HSP truncation");
            assert_eq!(actual.0[0].subject_id, "chr1");
            assert_eq!(actual.0[0].raw_subject_id.as_deref(), Some("1"));
            assert_eq!(actual.0[1].strand, "-");
            assert_eq!(actual.0[1].three_prime_mismatches, 0);
            let intended = PrimerSpecificityIntendedTarget::default();
            let products = |rows: &(
                Vec<PrimerSpecificityPrimerHit>,
                Vec<PrimerSpecificityPrimerHit>,
                Vec<String>,
            )| {
                let mut products = GentleEngine::primer_specificity_collect_amplicons_for_hits(
                    &rows.0, &rows.1, &policy,
                );
                GentleEngine::primer_specificity_finalize_amplicons(
                    &mut products,
                    &intended,
                    &policy,
                    crate::genomes::BlastDatabaseIndexKind::GenomicDna,
                );
                let summary = GentleEngine::primer_specificity_summary(
                    &rows.0,
                    &rows.1,
                    &products,
                    &intended,
                    crate::genomes::BlastDatabaseIndexKind::GenomicDna,
                    true,
                );
                serde_json::to_value((
                    products,
                    summary,
                    GentleEngine::primer_specificity_aggregate_warnings(rows.2.clone()),
                ))
                .unwrap()
            };
            assert_eq!(products(&actual), products(&expected));
            local.validate_unchanged().unwrap();
        }
    }
}

#[test]
fn primer_specificity_local_reader_retains_reference_warnings_with_complete_blast_strings() {
    let temp = tempfile::tempdir().unwrap();
    let f = fixture(temp.path(), "\n");
    fs::remove_file(&f.annotation).unwrap();
    let policy = PrimerSpecificityPolicy::default();
    let mut old = PreparedSequenceRegionReader::uncached(&f.catalog);
    let mut local = PreparedSequenceRegionReader::new(&f.catalog);
    let expected = evaluate(&mut old, &policy);
    let actual = evaluate(&mut local, &policy);
    assert_eq!(
        serde_json::to_value(&actual).unwrap(),
        serde_json::to_value(&expected).unwrap()
    );
    assert!(
        actual
            .2
            .iter()
            .any(|message| message.contains("Could not fetch subject sequence"))
    );
    assert_eq!(
        actual.0[0].full_alignment.as_ref().unwrap().algorithm,
        "blast_full_query_aligned_strings_v1"
    );
    assert_eq!(local.index_load_count(), 0);
}
