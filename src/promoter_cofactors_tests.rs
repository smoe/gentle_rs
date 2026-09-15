//! Synthetic package contract tests; real Parquet execution is explicit opt-in.

use super::*;
use crate::engine::Engine;
use serde_json::json;

fn write_json(path: &Path, value: &Value) {
    fs::write(path, serde_json::to_vec_pretty(value).unwrap()).unwrap();
}
fn inventory(root: &Path, names: &[&str]) -> Vec<Value> {
    names.iter().map(|name|{
    let path=root.join(name);json!({"path":name,"bytes":fs::metadata(&path).unwrap().len(),"sha256":sha256_file_hex(&path).unwrap()})
}).collect()
}
fn publish(root: &Path) {
    let names = [
        "anchors.parquet",
        "feature_MA9000.1.parquet",
        "anchor_promoter.parquet",
        "promoter.parquet",
        "promoter_gene.parquet",
        "cofactor_distance_isoform_comparison.parquet",
    ];
    let mut files = inventory(root, &names);
    let mut manifest: Value = serde_json::from_str(include_str!(
        "../test_files/fixtures/promoter_cofactors/manifest.template.json"
    ))
    .unwrap();
    manifest["files"] = json!(files);
    let identity = manifest["identity"].clone();
    write_json(&root.join("manifest.json"), &manifest);
    files.extend(inventory(root, &["manifest.json"]));
    write_json(
        &root.join("complete.json"),
        &json!({"identity":identity,"files":files}),
    );
}
fn package() -> tempfile::TempDir {
    let dir = tempfile::tempdir().unwrap();
    for name in [
        "anchors",
        "feature_MA9000.1",
        "anchor_promoter",
        "promoter",
        "promoter_gene",
        "cofactor_distance_isoform_comparison",
    ] {
        fs::write(
            dir.path().join(format!("{name}.parquet")),
            b"synthetic integrity-only placeholder",
        )
        .unwrap();
    }
    publish(dir.path());
    dir
}
fn request(root: &Path) -> PromoterCofactorRequest {
    PromoterCofactorRequest {
        package_path: root.display().to_string(),
        ..Default::default()
    }
}

/// Hand-crafted, complete report for transport tests; no production package.
pub(crate) fn handoff_report() -> PromoterCofactorReport {
    let p = package();
    let mut report = query(&request(p.path())).unwrap();
    report.anchors = vec![serde_json::from_value(json!({"anchor_id":1,"chrom":"1","anchor_start":100,"anchor_end":116,"anchor_score":0.25,
        "supported_tp73_saos2_TA":true,"depth_tp73_saos2_TA":9,"depth_negative_control_saos2_TA":0})).unwrap()];
    report.details = vec![serde_json::from_value(json!({"anchor_id":1,"motif_id":"MA9000.1","distance_band":"gap_6_20","hit_start":133,"hit_end":148,"best_score":4.25,
        "plus_score":4.25,"minus_score":-0.5,"best_strand":"+","interval_distance_bp":17,"genomic_side":"right","n_source_loci":2,"n_score_zero_loci":1,"present_at_requested_threshold":true})).unwrap()];
    report.promoters = vec![serde_json::from_value(json!({"regulatory_feature_id":"synthetic-promoter","chrom":"1","extended_start":90,"extended_end":140,
        "gene_links":[{"gene_id":"GENE-A","link_source":"synthetic","annotation_release":"toy"},{"gene_id":"GENE-B","link_source":"synthetic","annotation_release":"toy"}]})).unwrap()];
    report.memberships = vec![
        serde_json::from_value(json!({"anchor_id":1,"regulatory_feature_id":"synthetic-promoter"}))
            .unwrap(),
    ];
    report
}

fn hit_target() -> CofactorRegionTarget {
    CofactorRegionTarget::Hit {
        anchor_id: 1,
        motif_id: "MA9000.1".into(),
        distance_band: "gap_6_20".into(),
    }
}

#[test]
fn promoter_cofactors_capture_uses_shared_shell_and_preserves_both_orientations() {
    use crate::engine::{GentleEngine, Operation};
    use gentle_protocol::*;
    let source = handoff_report();
    let original = serde_json::to_value(&source).unwrap();
    for (orientation, start, end, strand) in [
        ("+", 33, 48, GenomicRegionStrand::Plus),
        ("-", 52, 67, GenomicRegionStrand::Minus),
    ] {
        let mut engine = GentleEngine::default();
        engine.state_mut().sequences.insert(
            "demo".into(),
            crate::dna_sequence::DNAsequence::from_sequence(&"ACGT".repeat(25)).unwrap(),
        );
        engine.state_mut().metadata.insert("provenance".into(), json!({"genome_extractions":[{"seq_id":"demo","genome_id":"GRCh38","chromosome":"1","start_1based":101,"end_1based":200,"anchor_strand":orientation,"anchor_verified":true}]}));
        let request = GenomicRegionCaptureRequest {
            set_id: "cofactors".into(),
            source: GenomicRegionCaptureSource::PromoterCofactor {
                report: Box::new(source.clone()),
                target: hit_target(),
                seq_id: Some("demo".into()),
            },
            ..Default::default()
        };
        let line = format!(
            "regions capture '{}'",
            serde_json::to_string(&request).unwrap()
        );
        let command = crate::engine_shell::parse_shell_line(&line).unwrap();
        let output = crate::engine_shell::execute_shell_command(&mut engine, &command).unwrap();
        assert!(output.state_changed);
        let inspected = engine
            .apply(Operation::ListGenomicRegions {
                request: Default::default(),
            })
            .unwrap();
        let json = serde_json::to_value(inspected.genomic_region_operation.unwrap()).unwrap();
        assert!(json.to_string().contains("cofactors"));
        let store = engine.genomic_region_store_snapshot().unwrap();
        let region = &store.sets[0].regions[0];
        assert_eq!(
            (
                region.interval.start_0based,
                region.interval.end_0based_exclusive
            ),
            (133, 148)
        );
        assert_eq!(region.interval.strand, GenomicRegionStrand::Plus);
        let projection = region.local_projection.as_ref().unwrap();
        assert_eq!(
            (
                projection.local_start_0based,
                projection.local_end_0based_exclusive,
                projection.local_strand
            ),
            (start, end, strand)
        );
        assert_eq!(region.evidence[0].source_record.as_ref(), Some(&original));
        assert_eq!(region.evidence[0].associated_gene_ids, ["GENE-A", "GENE-B"]);
        assert!(region.evidence[0].max_signal_value.is_none());
        let encoded = serde_json::to_vec(&original).unwrap();
        assert_eq!(
            region.evidence[0].source_sha256.as_deref(),
            Some(format!("sha256:{}", sha256_hex_bytes(&encoded)).as_str())
        );
        assert_eq!(engine.state().sequences["demo"].features().len(), 0);
        let mut wrong = request.clone();
        if let GenomicRegionCaptureSource::PromoterCofactor { report, .. } = &mut wrong.source {
            report.request.assembly = "mm10".into();
        }
        assert!(
            engine
                .apply(Operation::CaptureGenomicRegion { request: wrong })
                .is_err()
        );
        for genome_id in ["GRCh37", "Human GRCh38 Ensembl 116"] {
            engine.state_mut().metadata.get_mut("provenance").unwrap()["genome_extractions"][0]["genome_id"] =
                json!(genome_id);
            let mut mismatch = request.clone();
            mismatch.set_id = "must_not_create".into();
            assert!(
                engine
                    .apply(Operation::CaptureGenomicRegion { request: mismatch })
                    .is_err()
            );
            assert_eq!(
                engine.genomic_region_store_snapshot().unwrap().sets.len(),
                1
            );
        }
    }
    assert_eq!(serde_json::to_value(source).unwrap(), original);
}

#[test]
fn promoter_cofactors_region_rejects_tampered_source_record() {
    use crate::engine::{GentleEngine, Operation};
    use gentle_protocol::*;
    let (interval, mut evidence) =
        capture_region_evidence(&handoff_report(), &hit_target()).unwrap();
    evidence.source_record.as_mut().unwrap()["details"][0]["best_score"] = json!(300);
    let mut engine = GentleEngine::default();
    let err = engine
        .apply(Operation::CreateGenomicRegion {
            request: GenomicRegionCreateRequest {
                set_id: "tampered".into(),
                interval,
                evidence: vec![evidence],
                ..Default::default()
            },
        })
        .unwrap_err();
    assert!(
        err.to_string().contains("source record digest mismatch"),
        "{err}"
    );
    assert!(
        engine
            .genomic_region_store_snapshot()
            .unwrap()
            .sets
            .is_empty()
    );
}

#[test]
fn promoter_cofactors_capture_rejects_absence_bad_geometry_and_failed_evidence() {
    let report = handoff_report();
    for target in [
        CofactorRegionTarget::Anchor { anchor_id: 1 },
        CofactorRegionTarget::Promoter {
            regulatory_feature_id: "synthetic-promoter".into(),
        },
    ] {
        let (interval, _) = capture_region_evidence(&report, &target).unwrap();
        assert_eq!(
            interval.strand,
            gentle_protocol::GenomicRegionStrand::Unstranded
        );
    }
    let mut absent = report.clone();
    absent.details[0].hit_start = None;
    absent.details[0].hit_end = None;
    assert!(capture_region_evidence(&absent, &hit_target()).is_err());
    let mut bad = report.clone();
    bad.details[0].hit_end = Some(149);
    bad.details[0].interval_distance_bp = Some(18);
    assert!(capture_region_evidence(&bad, &hit_target()).is_err());
    let mut bad = report.clone();
    bad.availability = CofactorAvailability::QueryFailed;
    assert!(capture_region_evidence(&bad, &hit_target()).is_err());
    let mut tie = report.clone();
    tie.details[0].minus_score = Some(4.25);
    tie.details[0].best_strand = Some(".".into());
    assert_eq!(
        capture_region_evidence(&tie, &hit_target())
            .unwrap()
            .0
            .strand,
        gentle_protocol::GenomicRegionStrand::Unstranded
    );
    let mut r = PromoterCofactorRequest::default();
    assert_eq!(r.timeout_seconds, 30);
    r.timeout_seconds = 120;
    assert!(validate_request(&r).is_ok());
    r.timeout_seconds = 121;
    assert!(validate_request(&r).is_err());
}

#[test]
fn promoter_cofactors_inspection_is_relocatable_and_bound() {
    let a = package();
    let b = tempfile::tempdir().unwrap();
    for entry in fs::read_dir(a.path()).unwrap() {
        let entry = entry.unwrap();
        fs::copy(entry.path(), b.path().join(entry.file_name())).unwrap();
    }
    let first = query(&request(a.path())).unwrap();
    let second = query(&request(b.path())).unwrap();
    assert_eq!(first.availability, CofactorAvailability::Available);
    assert_eq!(first.report_id, second.report_id);
    assert_eq!(first.verified_file_sha256.len(), 7);
    assert!(!first.coverage.as_ref().unwrap().complete_genome_scan);
    assert!(
        first.coverage.unwrap().requested_candidates[0]
            .motif_ids
            .is_empty()
    );
}

#[test]
fn promoter_cofactors_integrity_assembly_coverage_and_runtime_fail_closed() {
    let p = package();
    let mut r = request(p.path());
    r.assembly = "mm10".into();
    assert_eq!(
        query(&r).unwrap().availability,
        CofactorAvailability::AssemblyMismatch
    );
    r.assembly = "GRCh38".into();
    r.query = CofactorQuery::AnchorDetail;
    r.anchor_id = Some(1);
    r.motif = Some("MA9001.1".into());
    assert_eq!(
        query(&r).unwrap().availability,
        CofactorAvailability::UnsupportedCoverage
    );
    r.motif = Some("MA9000.1".into());
    r.region = Some(CofactorRegion {
        chromosome: "X".into(),
        start_0based: 0,
        end_0based_exclusive: 100,
    });
    assert_eq!(
        query(&r).unwrap().availability,
        CofactorAvailability::UnsupportedCoverage
    );
    r.region = None;
    r.duckdb_executable = Some(p.path().join("missing-runtime").display().to_string());
    assert_eq!(
        query(&r).unwrap().availability,
        CofactorAvailability::RuntimeUnavailable
    );
    fs::write(p.path().join("feature_MA9000.1.parquet"), b"replacement").unwrap();
    let failed = query(&r).unwrap();
    assert_eq!(failed.availability, CofactorAvailability::InvalidPackage);
    assert!(failed.details.is_empty());
}

#[test]
fn promoter_cofactors_missing_package_and_unsafe_inventory() {
    let absent = tempfile::tempdir().unwrap();
    assert_eq!(
        query(&request(&absent.path().join("absent")))
            .unwrap()
            .availability,
        CofactorAvailability::PackageMissing
    );
    let p = package();
    let path = p.path().join("complete.json");
    let mut value: Value = serde_json::from_slice(&fs::read(&path).unwrap()).unwrap();
    value["files"][0]["path"] = json!("../escape");
    write_json(&path, &value);
    assert_eq!(
        query(&request(p.path())).unwrap().availability,
        CofactorAvailability::InvalidPackage
    );
}

#[test]
fn promoter_cofactors_replacement_at_same_path_changes_identity() {
    let p = package();
    let first = query(&request(p.path())).unwrap();
    fs::write(p.path().join("anchors.parquet"), b"new synthetic content").unwrap();
    publish(p.path());
    let second = query(&request(p.path())).unwrap();
    assert_eq!(second.availability, CofactorAvailability::Available);
    assert_ne!(first.report_id, second.report_id);
    assert_ne!(
        first.package_manifest_sha256,
        second.package_manifest_sha256
    );
}

#[test]
fn promoter_cofactors_rejects_unknown_version_and_duplicate_inventory() {
    let p = package();
    let path = p.path().join("manifest.json");
    let mut manifest: Value = serde_json::from_slice(&fs::read(&path).unwrap()).unwrap();
    manifest["schema_version"] = json!(999);
    write_json(&path, &manifest);
    assert_eq!(
        query(&request(p.path())).unwrap().availability,
        CofactorAvailability::InvalidPackage
    );
    publish(p.path());
    let path = p.path().join("complete.json");
    let mut complete: Value = serde_json::from_slice(&fs::read(&path).unwrap()).unwrap();
    let duplicate = complete["files"][0].clone();
    complete["files"].as_array_mut().unwrap().push(duplicate);
    write_json(&path, &complete);
    assert_eq!(
        query(&request(p.path())).unwrap().availability,
        CofactorAvailability::InvalidPackage
    );
}

#[cfg(unix)]
#[test]
fn promoter_cofactors_symlink_escape_and_stalled_runtime_are_bounded() {
    use std::os::unix::fs::{PermissionsExt, symlink};
    let p = package();
    let outside = tempfile::NamedTempFile::new().unwrap();
    fs::remove_file(p.path().join("anchors.parquet")).unwrap();
    symlink(outside.path(), p.path().join("anchors.parquet")).unwrap();
    assert_eq!(
        query(&request(p.path())).unwrap().availability,
        CofactorAvailability::InvalidPackage
    );
    fs::remove_file(p.path().join("anchors.parquet")).unwrap();
    fs::write(p.path().join("anchors.parquet"), b"restored").unwrap();
    publish(p.path());
    let runtime = p.path().join("stall.sh");
    fs::write(&runtime, "#!/bin/sh\nexec sleep 10\n").unwrap();
    fs::set_permissions(&runtime, fs::Permissions::from_mode(0o700)).unwrap();
    let r = PromoterCofactorRequest {
        query: CofactorQuery::Rankings,
        timeout_seconds: 1,
        duckdb_executable: Some(runtime.display().to_string()),
        ..request(p.path())
    };
    let start = Instant::now();
    let result = query(&r).unwrap();
    assert_eq!(
        result.availability,
        CofactorAvailability::RuntimeUnavailable
    );
    assert!(start.elapsed() < Duration::from_secs(5));
}

#[test]
fn promoter_cofactors_request_rejects_unsupported_inferences() {
    for value in [-2.0, f64::NAN, f64::INFINITY] {
        let r = PromoterCofactorRequest {
            presence_threshold: value,
            ..Default::default()
        };
        assert!(validate_request(&r).is_err());
    }
    let r = PromoterCofactorRequest {
        query: CofactorQuery::Rankings,
        gene_id: Some("GENE-A".into()),
        ..Default::default()
    };
    assert!(validate_request(&r).is_err());
    let r = PromoterCofactorRequest {
        max_rows: 2001,
        ..Default::default()
    };
    assert!(validate_request(&r).is_err());
    let r = PromoterCofactorRequest {
        query: CofactorQuery::AnchorDetail,
        ..Default::default()
    };
    assert!(validate_request(&r).is_err());
}

#[test]
fn promoter_cofactors_geometry_boundaries_and_strand_ties() {
    let a: CofactorAnchor = serde_json::from_value(
        json!({"anchor_id":1,"chrom":"1","anchor_start":100,"anchor_end":116,"anchor_score":0.25}),
    )
    .unwrap();
    for (start, band, gap) in [
        (110, "overlap", -6),
        (116, "adjacent_0_5", 0),
        (121, "adjacent_0_5", 5),
        (122, "gap_6_20", 6),
        (136, "gap_6_20", 20),
        (137, "gap_21_50", 21),
        (166, "gap_21_50", 50),
        (167, "gap_51_100", 51),
        (216, "gap_51_100", 100),
        (217, "gap_101_150", 101),
        (266, "gap_101_150", 150),
    ] {
        let mut d:CofactorDetail=serde_json::from_value(json!({"anchor_id":1,"motif_id":"MA9000.1","distance_band":band,"hit_start":start,"hit_end":start+15,"best_score":-0.5,"plus_score":-0.5,"minus_score":-0.5,"best_strand":".","interval_distance_bp":gap,"genomic_side":"right","n_source_loci":1,"n_score_zero_loci":0,"present_at_requested_threshold":false})).unwrap();
        validate_detail(&a, &d).unwrap();
        d.best_strand = Some("+".into());
        assert!(validate_detail(&a, &d).is_err());
    }
}

#[test]
fn promoter_cofactors_shell_and_operation_leave_state_unchanged() {
    let p = package();
    let r = request(p.path());
    let mut engine = crate::engine::GentleEngine::default();
    let before = serde_json::to_value(engine.state()).unwrap();
    let command = format!(
        "features promoter-cofactors '{}'",
        serde_json::to_string(&r).unwrap()
    );
    let parsed = crate::engine_shell::parse_shell_line(&command).unwrap();
    let result = crate::engine_shell::execute_shell_command(&mut engine, &parsed).unwrap();
    assert!(!result.state_changed);
    assert_eq!(result.output["report"]["availability"], "available");
    assert_eq!(before, serde_json::to_value(engine.state()).unwrap());
}

#[test]
fn promoter_cofactors_real_parquet_queries() {
    let Ok(executable) = std::env::var("GENTLE_TEST_DUCKDB") else {
        eprintln!("SKIP real Parquet: GENTLE_TEST_DUCKDB not set");
        return;
    };
    let dir = tempfile::tempdir().unwrap();
    let mut sql = include_str!("../test_files/fixtures/promoter_cofactors/fixture.sql").to_string();
    for (table, file) in [
        ("anchors", "anchors"),
        ("feature", "feature_MA9000.1"),
        ("anchor_promoter", "anchor_promoter"),
        ("promoter", "promoter"),
        ("promoter_gene", "promoter_gene"),
        (
            "cofactor_distance_isoform_comparison",
            "cofactor_distance_isoform_comparison",
        ),
    ] {
        sql.push_str(&format!(
            "COPY {table} TO {} (FORMAT PARQUET);",
            literal(&dir.path().join(format!("{file}.parquet")).to_string_lossy())
        ));
    }
    let out = std::process::Command::new(&executable)
        .args(["-no-init", ":memory:", "-c", &sql])
        .output()
        .unwrap();
    assert!(
        out.status.success(),
        "{}",
        String::from_utf8_lossy(&out.stderr)
    );
    publish(dir.path());
    let mut r = request(dir.path());
    // Includes integrity verification and several bounded DuckDB subprocesses.
    // Loaded developer/CI hosts get the existing maximum, not a product change.
    r.timeout_seconds = 120;
    r.duckdb_executable = Some(executable);
    r.query = CofactorQuery::AnchorDetail;
    r.anchor_id = Some(1);
    r.motif = Some("MA9000.1".into());
    let result = query(&r).unwrap();
    assert_eq!(
        result.availability,
        CofactorAvailability::Available,
        "{:?}",
        result.diagnostic
    );
    assert_eq!(result.anchors.len(), 1);
    assert_eq!(result.details.len(), 6);
    assert_eq!(result.promoters.len(), 2);
    assert_eq!(result.promoters[0].gene_links.len(), 2);
    assert_eq!(result.memberships.len(), 2);
    assert_eq!(result.details[2].best_score, Some(4.25));
    assert_eq!(result.details[2].interval_distance_bp, Some(17));
    assert!(result.details[0].best_score.is_none());
    assert_eq!(result.details[0].n_source_loci, 0);
    assert!(!result.details[0].present_at_requested_threshold);
    r.presence_threshold = 5.0;
    assert!(!query(&r).unwrap().details[2].present_at_requested_threshold);
    r.max_rows = 3;
    let limited = query(&r).unwrap();
    assert_eq!(limited.availability, CofactorAvailability::RowLimitExceeded);
    assert!(limited.details.is_empty());
    r.max_rows = 200;
    r.query = CofactorQuery::Anchors;
    r.anchor_id = None;
    r.motif = None;
    r.gene_id = Some("GENE-A".into());
    let by_gene = query(&r).unwrap();
    assert_eq!(by_gene.anchors.len(), 1);
    assert_eq!(by_gene.promoters.len(), 2);
    r.query = CofactorQuery::Rankings;
    r.gene_id = None;
    r.max_rows = 1;
    let ranking = query(&r).unwrap();
    assert_eq!(ranking.rankings.len(), 1);
    assert!(ranking.more_rankings_available);
    assert_eq!(ranking.rankings[0].motif_id, "MA9000.1");
    r.motif = Some("MA9001.1".into());
    assert_eq!(query(&r).unwrap().rankings.len(), 1);
    r.motif = Some("nothing' OR true --".into());
    assert!(query(&r).unwrap().rankings.is_empty());
    r.motif = None;
    r.ranking = CofactorRanking::DnEnriched;
    r.max_q_value = Some(0.05);
    assert!(query(&r).unwrap().rankings.is_empty());
    r.max_q_value = None;
    r.source_species = Some("not this species".into());
    assert!(query(&r).unwrap().rankings.is_empty());
    r.source_species = None;
    r.query = CofactorQuery::Anchors;
    r.region = Some(CofactorRegion {
        chromosome: "1".into(),
        start_0based: 116,
        end_0based_exclusive: 133,
    });
    assert_eq!(
        query(&r).unwrap().availability,
        CofactorAvailability::UnsupportedCoverage
    );
}
