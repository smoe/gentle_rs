//! Hand-crafted report/annotation pairs, not experimental or JASPAR scorer evidence.
//! Recreated by `fixture`; exercised by headless projection and native frame tests.

use super::*;
use gb_io::seq::{Feature, Location, Seq};
use gentle_protocol::genomic_motif_evidence::*;

pub(crate) fn fixture(minus: bool) -> (DNAsequence, TssProfileReport) {
    let mut report = crate::tss_profile_export::tests::synthetic_report();
    let window = &report.windows[usize::from(minus)];
    let g = &window.record.geometry;
    let mut seq = Seq::empty();
    seq.seq = b"ACGTA".to_vec();
    seq.comments = vec![
        format!(
            "GENtle promoter_id={}; sequence_sha256={}",
            window.record.promoter_id, window.record.sequence_sha256
        ),
        format!(
            "Reference={}; assembly={}; annotation_release={}; chromosome={}; genomic={}..{}; genomic_strand={}; local_axis=transcript_5prime_to_3prime; TSS_local_1based={}",
            report.reference.genome_id,
            report.reference.assembly,
            report.reference.annotation_release.as_deref().unwrap(),
            g.chromosome,
            g.start_1based,
            g.end_1based,
            g.strand.as_str(),
            g.upstream_bp + 1
        ),
    ];
    seq.features.push(Feature {
        kind: "misc_feature".into(),
        location: Location::single(g.upstream_bp as i64),
        qualifiers: vec![
            ("label".into(), Some("Annotated TSS candidate".into())),
            (
                "note".into(),
                Some(format!("Genomic {}; annotation-derived", g.tss_1based)),
            ),
        ],
    });
    let accession = "MA0001.1".to_string();
    let mut evidence = GenomicMotifEvidenceReport {
        schema: GENOMIC_MOTIF_EVIDENCE_SCHEMA.into(),
        report_id: "synthetic-sparse-query".into(),
        availability: GenomicMotifEvidenceAvailability::Available,
        request: GenomicMotifEvidenceRequest {
            motif_ids: vec![accession.clone()],
            ..Default::default()
        },
        provider: Some(GenomicMotifEvidenceProviderProvenance {
            genome_id: report.reference.genome_id.clone(),
            assembly_name: Some(report.reference.assembly.clone()),
            coordinate_mode: "bed_0based_half_open".into(),
            score_mode: "synthetic_package_log2".into(),
            provider_kind: "synthetic-test-only".into(),
            run_id: "synthetic-run".into(),
            manifest_sha256: "a".repeat(64),
            ..Default::default()
        }),
        motif_coverage: vec![GenomicMotifEvidenceMotifCoverage {
            motif_id: accession.clone(),
            source_minimum_score: Some(-1.0),
            density_limited: Some(true),
            status: GenomicMotifEvidenceCoverageStatus::TruncatedAtMaxRows,
            returned_hit_count: 2,
            ..Default::default()
        }],
        matched_hit_count: 3,
        returned_hit_count: 2,
        truncated: true,
        query_complete: false,
        ..Default::default()
    };
    for (i, (start, end, hit_start, strand, score)) in
        [(100, 102, 98, "+", -1.0), (202, 204, 202, "+", 8.0)]
            .into_iter()
            .enumerate()
    {
        let interval_id = format!("region-{i}");
        evidence.regions.push(GenomicMotifEvidenceResolvedRegion {
            interval_id: interval_id.clone(),
            requested_chromosome: g.chromosome.clone(),
            resolved_chromosome: Some(g.chromosome.clone()),
            start_0based: start,
            end_0based_exclusive: end,
            compatibility_status:
                GenomicMotifEvidenceCompatibilityStatus::ContigGeometryMatchedOnly,
            ..Default::default()
        });
        evidence.hits.push(GenomicMotifEvidenceHit {
            interval_id,
            chromosome: g.chromosome.clone(),
            start_0based: hit_start,
            end_0based_exclusive: hit_start + 3,
            motif_id: accession.clone(),
            strand: strand.into(),
            score,
            score_mode: "synthetic_package_log2".into(),
            ..Default::default()
        });
    }
    let digest = sha256_hex_bytes(&serde_json::to_vec(&evidence).unwrap());
    report
        .imported_motif_evidence
        .push(TssImportedMotifEvidence {
            source: TssInputBinding {
                role: "genomic_motif_evidence".into(),
                name: "synthetic-query.json".into(),
                sha256: digest.clone(),
            },
            report_sha256: digest,
            report: evidence,
        });
    (DNAsequence::from_genbank_seq(seq), report)
}

#[test]
fn report_arrays_and_imported_footprints_keep_strands_gaps_scales_and_provenance() {
    for minus in [false, true] {
        let (dna, report) = fixture(minus);
        let before = serde_json::to_vec(&report).unwrap();
        let base = TssSequenceView::from_dna(&dna).unwrap();
        let view = base.with_profile(&report).unwrap();
        let curve = view
            .lanes
            .iter()
            .find(|l| l.kind == TssLaneKind::ScoreTrace)
            .unwrap();
        let trace = curve.trace.as_ref().unwrap();
        assert_eq!(trace.forward, vec![Some(-2.0), None, Some(3.0)]);
        assert_eq!(trace.reverse, vec![Some(4.0), None, Some(-1.0)]);
        assert_eq!((curve.scale_min, curve.scale_max), (0.0, 4.0));
        assert_eq!(curve.units, "llr_bits");
        let imported = view
            .lanes
            .iter()
            .find(|l| l.kind == TssLaneKind::ImportedMotif)
            .unwrap();
        assert_eq!(imported.units, "synthetic_package_log2");
        assert_eq!((imported.scale_min, imported.scale_max), (-1.0, 8.0));
        assert!(imported.state.contains("PARTIAL"));
        assert!(imported.state.contains("TRUNCATED"));
        assert!(imported.details.contains("density limited Some(true)"));
        let hit = &imported.features[0];
        assert_eq!((hit.start, hit.end), (0, 2));
        assert_eq!(
            hit.reverse, minus,
            "same genomic + hit flips relative to negative-strand displayed DNA"
        );
        assert!(hit.clipped);
        assert_eq!(hit.feature_id, None);
        assert_eq!(
            view.geometry.genomic_at(0),
            Some(if minus { 204 } else { 100 })
        );
        assert!(
            hit.details
                .contains(&report.imported_motif_evidence[0].report_sha256)
        );
        assert_eq!(serde_json::to_vec(&report).unwrap(), before);
        let again = view.with_profile(&report).unwrap();
        assert_eq!(
            again.lanes.len(),
            view.lanes.len(),
            "reattaching must not duplicate evidence"
        );
    }
}

#[test]
fn same_sequence_at_another_locus_reference_or_tss_is_not_a_match() {
    let (dna, report) = fixture(false);
    let base = TssSequenceView::from_dna(&dna).unwrap();
    for field in [
        "promoter", "assembly", "genome", "release", "strand", "span", "hash",
    ] {
        let mut view = base.clone();
        match field {
            "promoter" => view.promoter_id = "other".into(),
            "assembly" => view.assembly = "other".into(),
            "genome" => view.genome_id = Some("other".into()),
            "release" => view.annotation_release = None,
            "strand" => view.geometry.strand = TssStrand::Minus,
            "span" => view.geometry.start_1based += 1,
            _ => view.sequence_sha256 = "f".repeat(64),
        }
        assert!(
            view.with_profile(&report).is_err(),
            "must reject {field} mismatch"
        );
    }
}

#[test]
fn report_scale_policy_keeps_zero_and_unavailable_distinct_without_division_by_zero() {
    let (_, mut report) = fixture(false);
    let accession = report.windows[0].tracks[0].accession.clone();
    for track in &mut report.windows[0].tracks {
        track.forward_scores.fill(Some(0.0));
        track.reverse_scores.fill(None);
    }
    assert_eq!(
        trace_range(&report, &report.windows[0], &accession),
        (0.0, 1.0, true)
    );
    report.panel_resolution.panel.scale_mode = TssScaleMode::SharedAcrossTss;
    assert_eq!(
        trace_range(&report, &report.windows[0], &accession),
        (0.0, 4.0, false)
    );
    report.panel_resolution.panel.clip_negative = false;
    assert_eq!(
        trace_range(&report, &report.windows[0], &accession),
        (-2.0, 4.0, false)
    );
    report.panel_resolution.panel.scale_mode = TssScaleMode::Shared;
    report.windows[0].tracks[1].forward_scores[0] = Some(7.0);
    assert_eq!(
        trace_range(&report, &report.windows[0], &accession),
        (0.0, 7.0, false)
    );
}

#[test]
fn report_content_validation_is_not_bypassed_by_matching_window_hash() {
    let (dna, mut report) = fixture(false);
    let view = TssSequenceView::from_dna(&dna).unwrap();
    report.imported_motif_evidence[0].report.hits[0].score = 99.0;
    assert!(view.with_profile(&report).unwrap_err().contains("hash"));
    report.imported_motif_evidence.clear();
    report.windows[0].tracks[0].forward_scores.pop();
    assert!(view.with_profile(&report).is_err());
}

#[test]
fn loading_json_binds_exact_file_bytes_and_detaching_restores_annotations() {
    let (dna, report) = fixture(false);
    let base = TssSequenceView::from_dna(&dna).unwrap();
    let temp = tempfile::tempdir().unwrap();
    let path = temp.path().join("report with spaces.json");
    let bytes = serde_json::to_vec_pretty(&report).unwrap();
    std::fs::write(&path, &bytes).unwrap();
    let mut loaded = base.load_profile(&path).unwrap();
    assert_eq!(
        loaded.profile.as_ref().unwrap().file_sha256,
        sha256_hex_bytes(&bytes)
    );
    loaded.clear_profile();
    assert_eq!(
        serde_json::to_vec(&loaded).unwrap(),
        serde_json::to_vec(&base).unwrap()
    );
    std::fs::write(&path, b"{}").unwrap();
    assert!(base.load_profile(&path).is_err());
    assert!(base.load_profile(temp.path()).is_err());
    std::fs::File::create(&path)
        .unwrap()
        .set_len(256 * 1024 * 1024 + 1)
        .unwrap();
    assert!(base.load_profile(&path).unwrap_err().contains("256 MiB"));
}
