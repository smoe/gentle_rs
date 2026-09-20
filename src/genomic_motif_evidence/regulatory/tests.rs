//! Offline reader regressions over synthetic, producer-layout DuckDB/Parquet.
//! No human genome, downloaded matrix, private evidence or remote process.

use super::*;
use crate::digest_utils::sha256_file_hex;

fn fixture() -> (tempfile::TempDir, PathBuf, String) {
    let tmp = tempfile::tempdir().unwrap();
    let root = tmp.path().join("delivered package with spaces and 'quote'");
    let duckdb = env::var(DUCKDB_BIN_ENV).unwrap_or_else(|_| "duckdb".into());
    let python = if cfg!(windows) { "python" } else { "python3" };
    let script = Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("test_files/fixtures/regulatory_motif_subset/make_fixture.py");
    let output = Command::new(python)
        .arg(script)
        .arg(&root)
        .arg("--duckdb")
        .arg(&duckdb)
        .current_dir(tmp.path())
        .output()
        .expect("Python and DuckDB required for this explicit integration test");
    assert!(
        output.status.success(),
        "{}",
        String::from_utf8_lossy(&output.stderr)
    );
    (tmp, root, duckdb)
}

fn request(root: &Path, duckdb: &str) -> GenomicMotifEvidenceRequest {
    GenomicMotifEvidenceRequest {
        package_root: Some(root.display().to_string()),
        duckdb_executable: Some(duckdb.into()),
        expected_genome_id: Some("fixture-genome".into()),
        motif_ids: vec!["MA0861.2".into(), "MA1961.2".into()],
        target: GenomicMotifEvidenceTarget::PackageTssWindows {
            gene_query: Some("TOY".into()),
            tss_id: None,
        },
        ..Default::default()
    }
}

fn run(
    r: &GenomicMotifEvidenceRequest,
    regions: &[GenomicMotifQueryRegion],
) -> GenomicMotifEvidenceReport {
    query_genomic_motif_evidence(r, regions).expect("valid typed request")
}

fn available(report: &GenomicMotifEvidenceReport) {
    assert_eq!(
        report.availability,
        GenomicMotifEvidenceAvailability::Available,
        "{:?}",
        report.warnings
    );
}

fn read(path: &Path) -> Value {
    serde_json::from_slice(&fs::read(path).unwrap()).unwrap()
}
fn write(path: &Path, value: &Value) {
    fs::write(path, serde_json::to_vec_pretty(value).unwrap()).unwrap();
}
fn rebind(root: &Path, name: &str) {
    let mut m = read(&root.join("manifest.json"));
    let f = m["files"]
        .as_array_mut()
        .unwrap()
        .iter_mut()
        .find(|f| f["path"] == name)
        .unwrap();
    f["sha256"] = sha256_file_hex(&root.join(name)).unwrap().into();
    f["bytes"] = fs::metadata(root.join(name)).unwrap().len().into();
    write(&root.join("manifest.json"), &m);
}

fn sql(root: &Path, duckdb: &str, query: &str) {
    let output = Command::new(duckdb)
        .args(["-no-init", "-batch", "-bail"])
        .arg(root.join("regulatory_tfbs.duckdb"))
        .arg("-c")
        .arg(query)
        .output()
        .unwrap();
    assert!(
        output.status.success(),
        "{}",
        String::from_utf8_lossy(&output.stderr)
    );
}

fn snapshot(root: &Path) -> BTreeMap<PathBuf, String> {
    let mut out = BTreeMap::new();
    for entry in fs::read_dir(root).unwrap() {
        let p = entry.unwrap().path();
        if p.is_dir() {
            out.extend(snapshot(&p));
        } else {
            out.insert(p.clone(), sha256_file_hex(&p).unwrap());
        }
    }
    out
}

#[test]
fn regulatory_catalog_targets_are_explicitly_bounded() {
    assert!(
        validate_target(&GenomicMotifEvidenceTarget::PackageCatalog {
            search: String::new(),
            offset: 0,
            limit: 100,
        })
        .is_ok()
    );
    assert!(
        validate_target(&GenomicMotifEvidenceTarget::PackageCatalog {
            search: String::new(),
            offset: 0,
            limit: 0,
        })
        .is_err()
    );
    assert!(
        validate_target(&GenomicMotifEvidenceTarget::PackageTssWindows {
            gene_query: Some("TEST".into()),
            tss_id: Some("tss".into()),
        })
        .is_err()
    );
    let mut r = GenomicMotifEvidenceRequest {
        target: GenomicMotifEvidenceTarget::PackageCatalog {
            search: String::new(),
            offset: 0,
            limit: 100,
        },
        ..Default::default()
    };
    assert!(validate_genomic_motif_request(&r, &[]).is_ok());
    r.minimum_score = Some(0.0);
    assert!(validate_genomic_motif_request(&r, &[]).is_err());
}

#[test]
fn regulatory_unknown_pseudocount_is_not_zero() {
    let mut provider = GenomicMotifEvidenceProviderProvenance::default();
    let mut count = None;
    merge_policy(&mut provider, &mut count, &serde_json::json!({})).unwrap();
    assert_eq!(count, None);
    merge_policy(
        &mut provider,
        &mut count,
        &serde_json::json!({"pseudocount":0.0}),
    )
    .unwrap();
    assert_eq!(count, Some(0.0));
    assert!(
        merge_policy(
            &mut provider,
            &mut count,
            &serde_json::json!({"pseudocount":1.0})
        )
        .is_err()
    );
}

#[test]
fn regulatory_paths_and_total_budget_fail_closed() {
    let tmp = tempfile::tempdir().unwrap();
    for invalid in [
        "../escape",
        "a/../b",
        "/absolute",
        "C:/outside",
        "a\\b",
        "a//b",
        "./manifest.json",
    ] {
        assert!(checked_file(tmp.path(), invalid).is_err(), "{invalid}");
    }
    let runtime = Runtime {
        executable: "not-run".into(),
        start: Instant::now() - Duration::from_secs(2),
        timeout: Duration::from_secs(1),
    };
    assert!(runtime.remaining().unwrap_err().contains("wall-clock"));
}

#[test]
#[ignore = "requires Python 3 and real DuckDB CLI; run explicitly, no mocked Parquets"]
fn regulatory_real_duckdb_query_catalog_strands_owners_and_replay() {
    let (_tmp, root, duckdb) = fixture();
    let before = snapshot(&root);
    let mut r = request(&root, &duckdb);
    let report = run(&r, &[]);
    available(&report);
    assert_eq!(report.hits.len(), 7);
    assert_eq!(
        report.selected_payload_emitted_hit_count, 7,
        "not original atlas counts"
    );
    assert!(report.query_complete);
    assert!(
        report
            .motif_coverage
            .iter()
            .all(|c| c.density_limited.is_none() && c.informative_threshold.is_none()),
        "do not fabricate atlas calibration statistics"
    );
    assert_eq!(
        report.hits.iter().filter(|h| h.start_0based == 295).count(),
        2
    );
    assert!(report.hits.iter().any(|h| h.score == -4.0));
    assert!(report.hits.iter().any(|h| h.score == -0.5));
    assert!(report.hits.iter().any(|h| h.score == 0.0));
    let subset = report.regulatory_subset.as_ref().unwrap();
    assert_eq!(subset.tss_windows.len(), 2);
    assert_eq!(
        subset.transcript_owners.len(),
        4,
        "shared physical TSS retains all owners"
    );
    assert_eq!(subset.regulatory_gene_links[0]["gene_id"], "geneB");
    assert_eq!(subset.tss_windows[0].start_0based, 300);
    assert_eq!(subset.tss_windows[0].end_0based_exclusive, 1301);
    assert_eq!(subset.tss_windows[1].start_0based, 1100);
    assert_eq!(subset.tss_windows[1].end_0based_exclusive, 2101);
    assert_eq!(subset.tss_windows[1].strand, "-");
    assert!(
        subset
            .hit_annotations
            .iter()
            .any(|a| a.query_interval_ids.len() == 2)
    );
    assert_eq!(
        serde_json::from_value::<GenomicMotifEvidenceReport>(
            serde_json::to_value(&report).unwrap()
        )
        .unwrap(),
        report
    );
    r.minimum_score = Some(0.0);
    let positive = run(&r, &[]);
    available(&positive);
    assert_eq!(positive.hits.len(), 4);
    r.minimum_score = Some(-6.0);
    let censored = run(&r, &[]);
    available(&censored);
    assert!(!censored.query_complete);
    assert!(
        censored
            .motif_coverage
            .iter()
            .all(|c| c.status == GenomicMotifEvidenceCoverageStatus::IncompleteBelowStorageFloor)
    );
    r.minimum_score = None;
    r.max_rows = 2;
    let truncated = run(&r, &[]);
    available(&truncated);
    assert!(truncated.truncated && !truncated.query_complete);
    r = request(&root, &duckdb);
    r.target = GenomicMotifEvidenceTarget::GenomicIntervals { intervals: vec![] };
    let mut reverse = region("reverse".into(), "1".into(), 1100, 2101);
    reverse.source_seq_id = Some("negative-locus".into());
    reverse.source_start_0based = Some(0);
    reverse.source_end_0based_exclusive = Some(1001);
    reverse.source_anchor_start_1based = Some(1101);
    reverse.source_anchor_end_1based = Some(2101);
    reverse.source_anchor_reverse = true;
    reverse.source_sequence_length_bp = Some(1001);
    let minus = run(&r, &[reverse]);
    available(&minus);
    let hit = minus.hits.iter().find(|h| h.start_0based == 1600).unwrap();
    assert_eq!(hit.source_start_0based, Some(485));
    assert_eq!(hit.source_end_0based_exclusive, Some(501));
    assert_eq!(hit.source_forward_strand, Some(true));
    let abut = run(&r, &[region("abut".into(), "1".into(), 311, 320)]);
    available(&abut);
    assert!(abut.hits.is_empty(), "half-open abutment is not overlap");
    r.motif_ids.clear();
    r.target = GenomicMotifEvidenceTarget::PackageCatalog {
        search: "patz".into(),
        offset: 0,
        limit: 1,
    };
    let catalog = run(&r, &[]);
    available(&catalog);
    assert_eq!(catalog.selected_payload_file_count, 0);
    let subset = catalog.regulatory_subset.unwrap();
    assert_eq!(subset.catalog_total_motifs, 2);
    assert_eq!(subset.catalog_matched_motifs, 1);
    assert_eq!(subset.motif_metadata[0]["species"], "fixture-other-taxon");
    assert!(
        !subset
            .verified_files
            .iter()
            .any(|f| f.path.starts_with("hits/"))
    );
    assert_eq!(
        snapshot(&root),
        before,
        "read-only package including catalog"
    );
    let relocated = root.with_file_name("relocated package");
    fs::rename(&root, &relocated).unwrap();
    r.package_root = Some(relocated.display().to_string());
    available(&run(&r, &[]));
}

#[test]
#[ignore = "requires Python 3 and real DuckDB CLI"]
fn regulatory_real_duckdb_empty_missing_corrupt_and_incompatible() {
    let (_tmp, root, duckdb) = fixture();
    let mut r = request(&root, &duckdb);
    r.target = GenomicMotifEvidenceTarget::GenomicIntervals { intervals: vec![] };
    let empty = run(&r, &[region("empty".into(), "2".into(), 300, 1301)]);
    available(&empty);
    assert!(empty.query_complete && empty.hits.is_empty());
    assert_eq!(empty.selected_payload_file_count, 2);
    let mt = run(&r, &[region("mt".into(), "MT".into(), 300, 1301)]);
    available(&mt);
    assert!(mt.query_complete && mt.hits.is_empty());
    assert_eq!(mt.selected_payload_file_count, 0);
    assert!(
        mt.regulatory_subset
            .unwrap()
            .coverage
            .iter()
            .all(|c| c.state == RegulatoryMotifCoverageState::KnownEmptyIntersection)
    );
    let unknown = run(&r, &[region("unknown".into(), "chr1".into(), 300, 1301)]);
    available(&unknown);
    assert!(
        !unknown.query_complete && unknown.hits.is_empty(),
        "no contig alias guesses"
    );
    assert!(
        unknown
            .motif_coverage
            .iter()
            .all(|c| c.status == GenomicMotifEvidenceCoverageStatus::NotAssessed)
    );
    assert!(
        unknown
            .warnings
            .iter()
            .all(|w| !w.contains("not represented in the package threshold registry"))
    );
    r.motif_ids = vec!["MA9999.1".into()];
    let missing = run(&r, &[region("missing".into(), "1".into(), 300, 1301)]);
    available(&missing);
    assert!(!missing.query_complete);
    assert_eq!(
        missing.motif_coverage[0].status,
        GenomicMotifEvidenceCoverageStatus::MotifNotInPackage
    );
    r = request(&root, &duckdb);
    r.expected_genome_id = Some("wrong".into());
    assert_eq!(
        run(&r, &[]).availability,
        GenomicMotifEvidenceAvailability::IncompatiblePackage
    );
    r = request(&root, &duckdb);
    r.duckdb_executable = Some(root.join("missing-executable").display().to_string());
    assert_eq!(
        run(&r, &[]).availability,
        GenomicMotifEvidenceAvailability::DuckdbUnavailable
    );
    r = request(&root, &duckdb);
    r.max_payload_files = 1;
    assert_eq!(
        run(&r, &[]).availability,
        GenomicMotifEvidenceAvailability::QueryFailed
    );
    r = request(&root, &duckdb);
    let payload = root.join("hits/chrom=1/MA0861.2.parquet");
    let original = fs::read(&payload).unwrap();
    let mut corrupted = original.clone();
    corrupted[10] ^= 1;
    fs::write(&payload, corrupted).unwrap();
    let failed = run(&r, &[]);
    assert_eq!(
        failed.availability,
        GenomicMotifEvidenceAvailability::InvalidPackage
    );
    assert!(failed.warnings[0].contains("SHA-256"));
    fs::write(&payload, &original).unwrap();
    available(&run(&r, &[]));
    fs::remove_file(&payload).unwrap();
    assert_eq!(
        run(&r, &[]).availability,
        GenomicMotifEvidenceAvailability::InvalidPackage
    );
}

#[test]
#[ignore = "requires Python 3 and real DuckDB CLI"]
fn regulatory_real_duckdb_rejects_views_bad_tags_and_metadata() {
    let (_tmp, root, duckdb) = fixture();
    let r = request(&root, &duckdb);
    let original = read(&root.join("manifest.json"));
    for (key, value) in [
        ("schema_version", serde_json::json!(2)),
        ("scope", serde_json::json!("regulatory_or_tss")),
        ("complete_genome_scan", serde_json::json!(true)),
    ] {
        let mut m = original.clone();
        m[key] = value;
        write(&root.join("manifest.json"), &m);
        assert_eq!(
            run(&r, &[]).availability,
            GenomicMotifEvidenceAvailability::InvalidPackage
        );
    }
    write(&root.join("manifest.json"), &original);
    let mut am = read(&root.join("annotation/manifest.json"));
    am["annotation_release"] = "wrong".into();
    write(&root.join("annotation/manifest.json"), &am);
    rebind(&root, "annotation/manifest.json");
    assert_eq!(
        run(&r, &[]).availability,
        GenomicMotifEvidenceAvailability::InvalidPackage
    );
    am["annotation_release"] = "fixture-gtf".into();
    write(&root.join("annotation/manifest.json"), &am);
    rebind(&root, "annotation/manifest.json");
    available(&run(&r, &[]));
    sql(
        &root,
        &duckdb,
        "UPDATE motif_metadata SET motif_id='MA9999.1' WHERE motif_id='MA0861.2';",
    );
    rebind(&root, "regulatory_tfbs.duckdb");
    let wrong_catalog = run(&r, &[]);
    assert_eq!(
        wrong_catalog.availability,
        GenomicMotifEvidenceAvailability::InvalidPackage
    );
    assert!(wrong_catalog.warnings[0].contains("different accessions"));
    sql(
        &root,
        &duckdb,
        "UPDATE motif_metadata SET motif_id='MA0861.2' WHERE motif_id='MA9999.1';",
    );
    rebind(&root, "regulatory_tfbs.duckdb");
    let inventory_original = read(&root.join("file_inventory.json"));
    let mut duplicate = inventory_original.clone();
    duplicate[1] = duplicate[0].clone();
    write(&root.join("file_inventory.json"), &duplicate);
    rebind(&root, "file_inventory.json");
    assert_eq!(
        run(&r, &[]).availability,
        GenomicMotifEvidenceAvailability::InvalidPackage
    );
    write(&root.join("file_inventory.json"), &inventory_original);
    rebind(&root, "file_inventory.json");

    // Rebind a deliberately nonqualifying bridging hit, so failure must be
    // semantic admission rather than a stale checksum.
    let payload = root.join("hits/chrom=1/MA0861.2.parquet");
    let original_payload = fs::read(&payload).unwrap();
    let replacement = root.join("replacement.parquet");
    sql(
        &root,
        &duckdb,
        &format!(
            "COPY (SELECT * REPLACE(false AS overlaps_regulatory_tss_intersection) FROM {}) TO {} (FORMAT PARQUET);",
            parquet(&payload),
            sql_string(&replacement.to_string_lossy())
        ),
    );
    fs::remove_file(&payload).unwrap();
    fs::rename(&replacement, &payload).unwrap();
    let hash = sha256_file_hex(&payload).unwrap();
    let size = fs::metadata(&payload).unwrap().len();
    let mut inventory = inventory_original.clone();
    let row = inventory
        .as_array_mut()
        .unwrap()
        .iter_mut()
        .find(|r| r["chrom"] == "1" && r["motif_id"] == "MA0861.2")
        .unwrap();
    row["sha256"] = hash.clone().into();
    row["bytes"] = size.into();
    write(&root.join("file_inventory.json"), &inventory);
    sql(
        &root,
        &duckdb,
        &format!(
            "UPDATE file_inventory SET sha256={},bytes={size} WHERE chrom='1' AND motif_id='MA0861.2';",
            sql_string(&hash)
        ),
    );
    let mut manifest = read(&root.join("manifest.json"));
    manifest["parquet_bytes"] = inventory
        .as_array()
        .unwrap()
        .iter()
        .map(|r| r["bytes"].as_u64().unwrap())
        .sum::<u64>()
        .into();
    write(&root.join("manifest.json"), &manifest);
    rebind(&root, "file_inventory.json");
    rebind(&root, "regulatory_tfbs.duckdb");
    let invalid = run(&r, &[]);
    assert_eq!(
        invalid.availability,
        GenomicMotifEvidenceAvailability::InvalidPackage
    );
    assert!(
        invalid.warnings[0].contains("intersection hit"),
        "{:?}",
        invalid.warnings
    );
    fs::write(&payload, original_payload).unwrap();
    write(&root.join("file_inventory.json"), &inventory_original);
    sql(
        &root,
        &duckdb,
        &format!(
            "DELETE FROM file_inventory; INSERT INTO file_inventory SELECT * FROM read_json_auto({});",
            sql_string(&root.join("file_inventory.json").to_string_lossy())
        ),
    );
    manifest["parquet_bytes"] = original["parquet_bytes"].clone();
    write(&root.join("manifest.json"), &manifest);
    rebind(&root, "file_inventory.json");
    rebind(&root, "regulatory_tfbs.duckdb");
    sql(
        &root,
        &duckdb,
        "ALTER TABLE motif_metadata RENAME TO hidden; CREATE VIEW motif_metadata AS SELECT * FROM hidden;",
    );
    rebind(&root, "regulatory_tfbs.duckdb");
    let failed = run(&r, &[]);
    assert_eq!(
        failed.availability,
        GenomicMotifEvidenceAvailability::InvalidPackage
    );
    assert!(failed.warnings[0].contains("materialized"));
}

#[test]
#[ignore = "requires Python 3 and real DuckDB CLI"]
fn regulatory_real_duckdb_large_inventory_and_shell_parity() {
    use crate::engine::GentleEngine;
    use crate::engine_shell::{
        execute_shell_command, parse_shell_line, shell_quote as quote_shell_arg,
    };
    let (_tmp, root, duckdb) = fixture();
    let line = format!(
        "features genomic-motif-evidence --gene TOY --package {} --duckdb {} --motif MA0861.2 --path {}",
        quote_shell_arg(&root.to_string_lossy()),
        quote_shell_arg(&duckdb),
        quote_shell_arg(&root.parent().unwrap().join("saved.json").to_string_lossy())
    );
    let mut engine = GentleEngine::default();
    let result = execute_shell_command(&mut engine, &parse_shell_line(&line).unwrap()).unwrap();
    assert!(result.output["report"]["regulatory_subset"].is_object());
    let saved: GenomicMotifEvidenceReport =
        serde_json::from_value(read(&root.parent().unwrap().join("saved.json"))).unwrap();
    available(&saved);
    assert_eq!(saved.hits.len(), 4);
    assert!(
        engine.state().sequences.is_empty(),
        "no DNA imported or annotated"
    );
    // Production shape: metadata-only 65,825 chromosome/motif rows. Non-selected
    // payload paths deliberately do not exist; no whole-atlas probing is allowed.
    let template = read(&root.join("file_inventory.json"));
    let mut inventory = Vec::new();
    for i in 0..22 {
        for m in 0..2633 {
            let mut row = template[0].clone();
            row["chrom"] = (i + 1).to_string().into();
            row["motif_id"] = format!("M{m}").into();
            row["path"] = format!("hits/chrom={}/M{m}.parquet", i + 1).into();
            row["rows"] = 0.into();
            inventory.push(row);
        }
    }
    for c in ["X", "Y", "MT"] {
        for m in 0..2633 {
            let mut row = template[0].clone();
            row["chrom"] = c.into();
            row["motif_id"] = format!("M{m}").into();
            row["rows"] = 0.into();
            row["path"] = format!("hits/chrom={c}/M{m}.parquet").into();
            if c == "MT" {
                row["path"] = Value::Null;
                row["sha256"] = Value::Null;
                row["bytes"] = 0.into();
                row["state"] = "known_empty_intersection".into();
            }
            inventory.push(row);
        }
    }
    assert_eq!(inventory.len(), 65_825);
    let mut m = read(&root.join("manifest.json"));
    let mut a = read(&root.join("annotation/manifest.json"));
    let chromosomes = (1..=22)
        .map(|n| n.to_string())
        .chain(["X".into(), "Y".into(), "MT".into()])
        .collect::<Vec<_>>();
    a["chromosomes"] = chromosomes
        .iter()
        .map(|c| {
            let mut v = a["chromosomes"][if c == "MT" { 2 } else { 0 }].clone();
            v["chrom"] = c.clone().into();
            v
        })
        .collect::<Vec<_>>()
        .into();
    m["chromosomes"] = serde_json::to_value(chromosomes).unwrap();
    m["motif_count"] = 2633.into();
    m["rows"] = 0.into();
    m["parquet_bytes"] = inventory
        .iter()
        .map(|r| r["bytes"].as_u64().unwrap())
        .sum::<u64>()
        .into();
    write(&root.join("file_inventory.json"), &Value::Array(inventory));
    write(&root.join("annotation/manifest.json"), &a);
    write(&root.join("manifest.json"), &m);
    rebind(&root, "file_inventory.json");
    rebind(&root, "annotation/manifest.json");
    let r = request(&root, &duckdb);
    let paths = resolve_package_paths(&r).unwrap();
    let runtime = Runtime {
        executable: duckdb,
        start: Instant::now(),
        timeout: Duration::from_secs(30),
    };
    let package = Package::open(&paths, &runtime).unwrap();
    assert_eq!(package.entries.len(), 65_825);
    // Exercise the public catalog path too, with all metadata rows but no large
    // payloads. A reader that probes every inventory file would fail here.
    sql(
        &root,
        &runtime.executable,
        &format!(
            "DROP TABLE motif_metadata; CREATE TABLE motif_metadata AS SELECT 'M'||i AS motif_id, 'Synthetic matrix '||i AS motif_name, 16 AS motif_length FROM range(2633) t(i); \
         DELETE FROM file_inventory; INSERT INTO file_inventory SELECT * FROM read_json_auto({});",
            sql_string(&root.join("file_inventory.json").to_string_lossy())
        ),
    );
    rebind(&root, "regulatory_tfbs.duckdb");
    let mut inspect = r;
    inspect.motif_ids.clear();
    inspect.target = GenomicMotifEvidenceTarget::PackageCatalog {
        search: String::new(),
        offset: 2600,
        limit: 20,
    };
    let catalog = run(&inspect, &[]);
    available(&catalog);
    let subset = catalog.regulatory_subset.unwrap();
    assert_eq!(subset.catalog_total_motifs, 2633);
    assert_eq!(subset.motif_metadata.len(), 20);
    assert!(subset.catalog_has_more);
    assert_eq!(catalog.selected_payload_file_count, 0);
}
