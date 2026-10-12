//! Explicit, bounded attachment of saved sparse motif queries to TSS exports.
//! Rendering never queries DuckDB or alters source reports/local score arrays.

use super::*;
use std::io::Read;

pub(super) fn attach(
    report: &mut TssProfileReport,
    paths: &[String],
    should_continue: &mut dyn FnMut() -> bool,
) -> Result<(), EngineError> {
    use gentle_protocol::tss_motif_evidence::MAX_EVIDENCE_REPORTS;
    if paths.len() + report.imported_motif_evidence.len() > MAX_EVIDENCE_REPORTS {
        return Err(invalid(
            "At most 64 imported motif evidence reports are supported",
        ));
    }
    let mut bytes_read = 0;
    for path in paths {
        if !should_continue() {
            return Err(GentleEngine::tfbs_cancelled_error(
                "imported motif evidence",
            ));
        }
        let file = crate::tss_fasta_bundle::open_regular_input(Path::new(path))
            .map_err(|e| invalid(format!("Imported motif evidence input: {e}")))?;
        let mut bytes = Vec::new();
        file.take(32 * 1024 * 1024 + 1)
            .read_to_end(&mut bytes)
            .map_err(|e| invalid(format!("Read imported motif evidence: {e}")))?;
        bytes_read += bytes.len();
        if bytes.len() > 32 * 1024 * 1024 || bytes_read > 64 * 1024 * 1024 {
            return Err(invalid(
                "Imported motif evidence exceeds 32 MiB/file or 64 MiB total",
            ));
        }
        let evidence: gentle_protocol::genomic_motif_evidence::GenomicMotifEvidenceReport =
            serde_json::from_slice(&bytes).map_err(|e| {
                invalid(format!("Invalid saved genomic motif evidence report: {e}"))
            })?;
        let source_hash = sha256_hex_bytes(&bytes);
        let canonical = serde_json::to_vec(&evidence).map_err(|e| invalid(e.to_string()))?;
        report
            .imported_motif_evidence
            .push(TssImportedMotifEvidence {
                source: TssInputBinding {
                    role: "genomic_motif_evidence".into(),
                    name: format!("motif-evidence-{source_hash}.json"),
                    sha256: source_hash,
                },
                report_sha256: sha256_hex_bytes(&canonical),
                report: evidence,
            });
    }
    crate::tss_profile_export::validate_tss_profile_report(report)
}

#[cfg(test)]
mod tests {
    use super::*;
    use gentle_protocol::genomic_motif_evidence::*;

    // Hand-crafted sparse-query results over the existing synthetic export
    // fixture. No live DuckDB, biological evidence or downloaded sequences.
    fn evidence(report: &TssProfileReport) -> GenomicMotifEvidenceReport {
        let mut r = GenomicMotifEvidenceReport {
            schema: GENOMIC_MOTIF_EVIDENCE_SCHEMA.into(),
            report_id: "synthetic-imported-hits".into(),
            availability: GenomicMotifEvidenceAvailability::Available,
            request: GenomicMotifEvidenceRequest {
                motif_ids: vec!["MA0001.1".into()],
                ..Default::default()
            },
            provider: Some(GenomicMotifEvidenceProviderProvenance {
                genome_id: report.reference.genome_id.clone(),
                assembly_name: Some(report.reference.assembly.clone()),
                coordinate_mode: "bed_0based_half_open".into(),
                score_mode: "log2_relative_risk".into(),
                run_id: "synthetic-run".into(),
                manifest_sha256: format!("sha256:{}", "a".repeat(64)),
                ..Default::default()
            }),
            motif_coverage: vec![GenomicMotifEvidenceMotifCoverage {
                motif_id: "MA0001.1".into(),
                source_minimum_score: Some(-1.0),
                status: GenomicMotifEvidenceCoverageStatus::TruncatedAtMaxRows,
                ..Default::default()
            }],
            truncated: true,
            query_complete: false,
            ..Default::default()
        };
        for (i, w) in report.windows.iter().enumerate() {
            let g = &w.record.geometry;
            let id = format!("region-{i}");
            r.regions.push(GenomicMotifEvidenceResolvedRegion {
                interval_id: id.clone(),
                resolved_chromosome: Some(g.chromosome.clone()),
                start_0based: g.start_1based - 1,
                end_0based_exclusive: g.end_1based,
                compatibility_status:
                    GenomicMotifEvidenceCompatibilityStatus::ContigGeometryMatchedOnly,
                ..Default::default()
            });
            for (offset, strand, score) in [(0, "+", 4.0), (1, "-", 2.0), (3, "+", -0.5)] {
                r.hits.push(GenomicMotifEvidenceHit {
                    interval_id: id.clone(),
                    chromosome: g.chromosome.clone(),
                    start_0based: g.start_1based - 1 + offset,
                    end_0based_exclusive: g.start_1based + 2 + offset,
                    motif_id: "MA0001.1".into(),
                    strand: strand.into(),
                    score,
                    score_mode: "log2_relative_risk".into(),
                    ..Default::default()
                });
            }
        }
        r.returned_hit_count = r.hits.len();
        r.matched_hit_count = r.hits.len();
        r.motif_coverage[0].returned_hit_count = r.hits.len();
        r
    }

    #[test]
    fn tss_imported_motif_export_binds_sources_preserves_scores_and_replays_offline() {
        let report = crate::tss_profile_export::tests::synthetic_report();
        let before = serde_json::to_value(&report.windows).unwrap();
        let scratch = tempfile::tempdir().unwrap();
        let root = std::fs::canonicalize(scratch.path()).unwrap();
        let path = root.join("hits.json");
        let evidence = evidence(&report);
        let bytes = serde_json::to_vec_pretty(&evidence).unwrap();
        std::fs::write(&path, &bytes).unwrap();
        let output = root.join("enriched");
        let request = ExportTssProfilesRequest {
            output_dir: output.display().to_string(),
            context_manifest: None,
            genomic_motif_evidence: vec![path.display().to_string()],
            rendering: TssProfileRenderOptions::default(),
            formats: vec![TssExportFormat::Svg],
        };
        let (enriched, receipt) =
            GentleEngine::export_tss_profile_report(&report, &request, &mut |_| true).unwrap();
        let enriched = enriched.unwrap();
        assert_eq!(serde_json::to_value(&enriched.windows).unwrap(), before);
        assert_eq!(
            enriched.imported_motif_evidence[0].source.sha256,
            sha256_hex_bytes(&bytes)
        );
        assert_eq!(
            serde_json::to_value(&enriched.imported_motif_evidence[0].report).unwrap(),
            serde_json::to_value(&evidence).unwrap()
        );
        assert!(
            receipt
                .inputs
                .iter()
                .any(|i| i.role == "genomic_motif_evidence")
        );
        crate::tss_profile_export::verify_tss_profile_receipt(&output, &receipt).unwrap();
        let tsv = std::fs::read_to_string(output.join("imported-motif-hits.tsv")).unwrap();
        assert_eq!(tsv.lines().count(), 7);
        assert!(tsv.contains("\ttrue\t-0.5\tlog2_relative_risk\t"));
        for file in receipt.outputs.keys().filter(|n| n.ends_with(".svg")) {
            let svg = std::fs::read_to_string(output.join(file)).unwrap();
            assert_eq!(svg.matches("data-role=\"imported-motif-hit\"").count(), 3);
            assert!(svg.contains("truncated_at_max_rows"));
            assert!(svg.contains("genomic +"));
        }
        let mut tampered = enriched.clone();
        tampered.imported_motif_evidence[0].report.hits[0].score += 1.0;
        assert!(
            crate::tss_profile_export::validate_tss_profile_report(&tampered)
                .unwrap_err()
                .message
                .contains("hash")
        );
        std::fs::remove_file(path).unwrap();
        let replay = root.join("replay");
        let (unchanged, receipt) = GentleEngine::export_tss_profile_report(
            &enriched,
            &ExportTssProfilesRequest {
                genomic_motif_evidence: vec![],
                output_dir: replay.display().to_string(),
                ..request
            },
            &mut |_| true,
        )
        .unwrap();
        assert!(unchanged.is_none());
        crate::tss_profile_export::verify_tss_profile_receipt(&replay, &receipt).unwrap();
        assert_eq!(
            std::fs::read(replay.join("report.json")).unwrap(),
            std::fs::read(output.join("report.json")).unwrap()
        );
    }

    #[test]
    fn tss_imported_motif_attachment_rejects_wrong_reference_geometry_and_scores() {
        let report = crate::tss_profile_export::tests::synthetic_report();
        let base = evidence(&report);
        let scratch = tempfile::tempdir().unwrap();
        let path = scratch.path().join("hits.json");
        let mutations: Vec<Box<dyn Fn(&mut GenomicMotifEvidenceReport)>> = vec![
            Box::new(|r| r.provider.as_mut().unwrap().genome_id.push_str("-other")),
            Box::new(|r| r.provider.as_mut().unwrap().assembly_name = Some("other".into())),
            Box::new(|r| r.provider.as_mut().unwrap().coordinate_mode = "1based".into()),
            Box::new(|r| r.hits[0].strand = "?".into()),
            Box::new(|r| r.hits[0].motif_id = "MA0001.99".into()),
            Box::new(|r| r.hits[0].end_0based_exclusive += 1),
            Box::new(|r| r.hits[0].score_mode = "tail".into()),
            Box::new(|r| r.hits[0].chromosome = "other".into()),
            Box::new(|r| r.motif_coverage.clear()),
            Box::new(|r| r.returned_hit_count = 0),
        ];
        for mutate in mutations {
            let mut r = base.clone();
            mutate(&mut r);
            std::fs::write(&path, serde_json::to_vec(&r).unwrap()).unwrap();
            assert!(
                attach(
                    &mut report.clone(),
                    &[path.display().to_string()],
                    &mut || true
                )
                .is_err()
            );
        }
        let mut r = report.clone();
        std::fs::write(&path, serde_json::to_vec(&base).unwrap()).unwrap();
        assert!(
            attach(
                &mut r,
                &[path.display().to_string(), path.display().to_string()],
                &mut || true,
            )
            .is_err()
        );
    }

    #[test]
    fn tss_imported_motif_attachment_cancels_before_reading_or_publishing() {
        let report = crate::tss_profile_export::tests::synthetic_report();
        let scratch = tempfile::tempdir().unwrap();
        let output = scratch.path().join("cancelled");
        let result = GentleEngine::export_tss_profile_report(
            &report,
            &ExportTssProfilesRequest {
                output_dir: output.display().to_string(),
                context_manifest: None,
                genomic_motif_evidence: vec!["this-input-does-not-exist.json".into()],
                rendering: Default::default(),
                formats: vec![TssExportFormat::Svg],
            },
            &mut |_| false,
        );
        assert!(
            result
                .unwrap_err()
                .message
                .to_ascii_lowercase()
                .contains("cancel")
        );
        assert!(!output.exists());
        assert!(report.imported_motif_evidence.is_empty());
    }
}
