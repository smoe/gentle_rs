//! Focused regression tests for strict TSS panel resolution and profile scoring.
//!
//! Fixture provenance: all MA999x accessions, names, counts, sequences, loci,
//! memberships, and comparison signals below are hand-crafted synthetic data,
//! not JASPAR entries or reference-genome evidence. Recreate them by running this
//! module: helpers serialize the literal matrices into an in-memory registry and
//! construct transcript-oriented `VerifiedTssBundle` records with SHA-256 bindings.
//! They exercise the resolver, verified-bundle scorer, and comparison arithmetic;
//! FASTA parsing/reference authenticity belong to the bundle-reader tests.
//! The legacy-scorer parity test uses equivalent literal DNA motifs under the
//! registry test lock, without changing the runtime registry or any fixture file.
//!
//! The typed-operation smoke uses temporary copies of the documented synthetic
//! windows in `test_files/fixtures/tss_profiles`, plus its panel pinning the
//! bundled JASPAR MA0004.1/Arnt PFM. It exports SVG/PNG/raster-backed PDF twice,
//! including a report-only replay after removing those temporary inputs. Set
//! `GENTLE_TEST_TSS_PROFILE_OUTPUT_DIR` to an existing scratch directory to retain
//! the uniquely named output directory on success; otherwise artifacts are temporary.
//! A separate two-TSS/three-matrix smoke uses only the in-memory synthetic registry
//! and verified-bundle helpers, then audits JSON/TSV/SVG projections and receipts.

use super::*;
use crate::tf_motifs::TfMotifDb;
use crate::tss_fasta_bundle::VerifiedTssBundle;

const FACTOR: &str = "SyntheticFactor";
const IDS: [&str; 3] = ["MA9993.7", "MA9991.2", "MA9992.1"];
const AC_COUNTS: [[f64; 4]; 2] = [[6.0, 2.0, 1.0, 1.0], [1.0, 6.0, 2.0, 1.0]];
const GT_COUNTS: [[f64; 4]; 2] = [[1.0, 1.0, 6.0, 2.0], [2.0, 1.0, 1.0, 6.0]];
const ACG_COUNTS: [[f64; 4]; 3] = [
    [6.0, 2.0, 1.0, 1.0],
    [1.0, 6.0, 2.0, 1.0],
    [1.0, 1.0, 6.0, 2.0],
];

fn motif_entry(id: &str, consensus: &str, counts: &[[f64; 4]]) -> serde_json::Value {
    json!({
        "id": id,
        "name": FACTOR,
        "consensus_iupac": consensus,
        "pfm": {
            "a": counts.iter().map(|column| column[0]).collect::<Vec<_>>(),
            "c": counts.iter().map(|column| column[1]).collect::<Vec<_>>(),
            "g": counts.iter().map(|column| column[2]).collect::<Vec<_>>(),
            "t": counts.iter().map(|column| column[3]).collect::<Vec<_>>(),
        },
    })
}

fn registry(entries: Vec<serde_json::Value>) -> (String, TfMotifDb) {
    let text = serde_json::to_string(&json!({
        "schema": "gentle.tf_motifs.v1",
        "motifs": entries,
    }))
    .unwrap();
    let registry = TfMotifDb::from_json_for_test(&text).expect("synthetic registry JSON");
    (text, registry)
}

fn three_matrix_registry() -> (String, TfMotifDb) {
    registry(vec![
        motif_entry(IDS[2], "ACG", &ACG_COUNTS),
        motif_entry(IDS[0], "AC", &AC_COUNTS),
        motif_entry(IDS[1], "GT", &GT_COUNTS),
    ])
}

fn panel(ids: &[&str], score_kind: TfbsScoreTrackValueKind) -> JasparTargetPanel {
    JasparTargetPanel {
        schema: PANEL_SCHEMA.into(),
        panel_id: "synthetic-strict-panel".into(),
        label: "Synthetic ordered matrices, not biological evidence".into(),
        score_kind: score_kind.as_str().into(),
        clip_negative: false,
        scale_mode: TssScaleMode::Independent,
        strand_policy: TssStrandPolicy::Both,
        calibration_state: TssCalibrationState::MatrixSpecific,
        calibration_statement: "Synthetic matrix-specific scores are not cross-calibrated".into(),
        calibration_id: None,
        calibration_sha256: None,
        top_hit_count: 5,
        factors: ids
            .iter()
            .enumerate()
            .map(|(index, id)| JasparPanelTrack {
                source_id: (*id).into(),
                factor_id: FACTOR.into(),
                label: format!("Synthetic row {} ({id})", index + 1),
                display_order: (index + 1) * 10,
                color_hint: Some("#125a7C".into()),
                score_kind: Some(score_kind.as_str().into()),
                track_id: None,
                provider_kind: None,
                factor_label: None,
            })
            .collect(),
    }
}

fn resolve(panel: JasparTargetPanel, registry: &TfMotifDb) -> TssPanelResolution {
    let bytes = serde_json::to_vec(&panel).unwrap();
    GentleEngine::resolve_tss_panel(panel, &bytes, registry)
        .expect("resolve synthetic strict panel")
}

fn assert_rejected(panel: JasparTargetPanel, registry: &TfMotifDb, diagnostic: &str) {
    let bytes = serde_json::to_vec(&panel).unwrap();
    let error = GentleEngine::resolve_tss_panel(panel, &bytes, registry)
        .expect_err("strict panel must be rejected");
    assert_eq!(error.code, ErrorCode::InvalidInput);
    assert!(
        error.message.contains(diagnostic),
        "expected {diagnostic:?} in {:?}",
        error.message
    );
}

fn record(
    id: &str,
    strand: TssStrand,
    tss_1based: u64,
    sequence: &str,
    selected: bool,
) -> (TssRecord, String, bool) {
    assert!(!sequence.is_empty());
    let upstream_bp = 3.min(sequence.len() - 1);
    let downstream_bp = sequence.len() - upstream_bp - 1;
    let (left, right) = match strand {
        TssStrand::Plus => (upstream_bp, downstream_bp),
        TssStrand::Minus => (downstream_bp, upstream_bp),
    };
    let record = TssRecord {
        promoter_id: id.into(),
        gene_id: "synthetic_gene_a".into(),
        gene_symbol: "SyntheticGeneA".into(),
        geometry: TssGeometry {
            chromosome: "synthetic_chr1".into(),
            strand,
            tss_1based,
            start_1based: tss_1based - left as u64,
            end_1based: tss_1based + right as u64,
            upstream_bp,
            downstream_bp,
        },
        transcripts: vec![format!("{id}_tx2"), format!("{id}_tx1")],
        sequence_sha256: sha256_hex_bytes(sequence.as_bytes()),
    };
    gentle_engine::tss_profiles::validate_record(&record).unwrap();
    (record, sequence.into(), selected)
}

fn bundle(records: Vec<(TssRecord, String, bool)>) -> VerifiedTssBundle {
    let reference = TssReference {
        genome_id: "synthetic_genome".into(),
        assembly: "synthetic_assembly_v1".into(),
        annotation_release: Some("synthetic_annotation_v2".into()),
    };
    // This private scoring seam binds a synthetic in-memory manifest, not a
    // reader-verified file bundle. The typed-operation smoke covers that reader.
    let fasta = records
        .iter()
        .map(|(record, sequence, _)| format!(">{}\n{sequence}\n", record.promoter_id))
        .collect::<String>();
    let bytes = serde_json::to_vec(&TssBundleManifest {
        schema: BUNDLE_SCHEMA.into(),
        reference: reference.clone(),
        fasta_files: BTreeMap::from([("in-memory.fa".into(), sha256_hex_bytes(fasta.as_bytes()))]),
        records: records
            .iter()
            .map(|(record, _, _)| record.clone())
            .collect(),
    })
    .unwrap();
    let manifest_sha256 = sha256_hex_bytes(&bytes);
    VerifiedTssBundle {
        reference,
        source: TssBundleSource {
            schema: BUNDLE_SCHEMA.into(),
            manifest_sha256: manifest_sha256.clone(),
            source_revision: None,
            dataset_id: None,
            producer_sha256: None,
        },
        selection_evidence: BTreeMap::new(),
        records,
        inputs: vec![TssInputBinding {
            role: "bundle_manifest".into(),
            name: "in-memory-manifest.json".into(),
            sha256: manifest_sha256,
        }],
        warnings: vec!["Synthetic records; no independent reference verification".into()],
    }
}

fn engine_with_sentinel() -> GentleEngine {
    let mut engine = GentleEngine::new();
    engine.state.metadata.insert(
        "tss_test_sentinel".into(),
        json!({"keep": ["existing", "project", "state"]}),
    );
    engine
}

fn assert_close(actual: f64, expected: f64) {
    assert!(
        actual.is_finite() && (actual - expected).abs() <= 1e-10,
        "expected {expected:.15}, got {actual:.15}"
    );
}

#[test]
fn strict_resolution_freezes_three_same_factor_matrices_in_declared_order() {
    let (registry_bytes, registry) = three_matrix_registry();
    let panel = panel(&IDS, TfbsScoreTrackValueKind::LlrBits);
    let panel_bytes = serde_json::to_vec_pretty(&panel).unwrap();
    let expected_panel = serde_json::to_value(&panel).unwrap();
    let resolution =
        GentleEngine::resolve_tss_panel(panel.clone(), &panel_bytes, &registry).unwrap();

    assert_eq!(
        serde_json::to_value(&resolution.panel).unwrap(),
        expected_panel
    );
    assert_eq!(resolution.panel_sha256, sha256_hex_bytes(&panel_bytes));
    assert!(resolution.registry_source_url.is_none());
    let active = resolution
        .registry_sources
        .iter()
        .find(|binding| binding.role == "active_registry")
        .expect("actual registry byte binding");
    assert_eq!(active.sha256, sha256_hex_bytes(registry_bytes.as_bytes()));
    assert_eq!(resolution.matrices.len(), 3);
    for (index, (expected_counts, version)) in [
        (AC_COUNTS.as_slice(), "7"),
        (GT_COUNTS.as_slice(), "2"),
        (ACG_COUNTS.as_slice(), "1"),
    ]
    .into_iter()
    .enumerate()
    {
        let matrix = &resolution.matrices[index];
        assert_eq!(matrix.specification.source_id, IDS[index]);
        assert_eq!(matrix.specification.factor_id, FACTOR);
        assert_eq!(matrix.specification.display_order, (index + 1) * 10);
        assert_eq!(matrix.version, version);
        assert_eq!(matrix.matrix_counts, expected_counts);
        let expected_identity =
            serde_json::to_vec(&(IDS[index], Some(FACTOR), expected_counts)).unwrap();
        assert_eq!(matrix.matrix_sha256, sha256_hex_bytes(&expected_identity));
    }
    let compact = resolve(panel, &registry);
    assert_ne!(compact.panel_sha256, resolution.panel_sha256);
    assert_eq!(
        serde_json::to_value(compact.matrices).unwrap(),
        serde_json::to_value(resolution.matrices).unwrap()
    );
}

#[test]
fn strict_resolution_requires_exact_case_and_version_without_alias_or_latest_fallback() {
    let (_, registry) = registry(vec![
        motif_entry("MA9993.9", "GT", &GT_COUNTS),
        motif_entry(IDS[0], "AC", &AC_COUNTS),
    ]);
    let valid = panel(&IDS[..1], TfbsScoreTrackValueKind::LlrBits);
    let resolved = resolve(valid.clone(), &registry);
    assert_eq!(resolved.matrices[0].version, "7");
    assert_eq!(resolved.matrices[0].matrix_counts, AC_COUNTS);

    for factor_id in ["syntheticfactor", "SYNTHETICFACTOR", "OtherFactor"] {
        let mut input = valid.clone();
        input.factors[0].factor_id = factor_id.into();
        assert_rejected(input, &registry, IDS[0]);
    }
    for accession in [FACTOR, "MA9993", "ma9993.7", "MA9993.8", "MA9994.1"] {
        let mut input = valid.clone();
        input.factors[0].source_id = accession.into();
        assert_rejected(input, &registry, accession);
    }
}

#[test]
fn strict_resolution_rejects_duplicate_panel_or_registry_ids_and_conflicting_order() {
    let (_, registry) = three_matrix_registry();
    let valid = panel(&IDS[..2], TfbsScoreTrackValueKind::LlrBits);
    let mut duplicate = valid.clone();
    duplicate.factors[1].source_id = IDS[0].into();
    assert_rejected(duplicate, &registry, "Duplicate");

    let mut tied = valid.clone();
    tied.factors[1].display_order = tied.factors[0].display_order;
    assert_rejected(tied, &registry, "display_order");
    let mut reordered = valid;
    reordered.factors.swap(0, 1);
    assert_rejected(reordered, &registry, "display_order");

    for duplicate_counts in [&AC_COUNTS, &GT_COUNTS] {
        let (_, duplicate_registry) = self::registry(vec![
            motif_entry(IDS[0], "AC", &AC_COUNTS),
            motif_entry(IDS[0], "AC", duplicate_counts),
        ]);
        assert_rejected(
            panel(&IDS[..1], TfbsScoreTrackValueKind::LlrBits),
            &duplicate_registry,
            "expected one exact registry entry",
        );
    }
}

#[test]
fn strict_resolution_rejects_duplicates_even_when_one_consensus_is_empty() {
    for consensus in ["", " \t "] {
        let valid = motif_entry(IDS[0], "AC", &AC_COUNTS);
        let skipped = motif_entry(IDS[0], consensus, &AC_COUNTS);
        for entries in [
            vec![valid.clone(), skipped.clone()],
            vec![skipped.clone(), valid.clone()],
        ] {
            let (_, registry) = registry(entries);
            assert_eq!(
                registry
                    .resolve(IDS[0])
                    .expect("legacy lookup retains the valid row")
                    .matrix_counts,
                AC_COUNTS
            );
            assert_rejected(
                panel(&IDS[..1], TfbsScoreTrackValueKind::LlrBits),
                &registry,
                "found 2",
            );
        }
    }
}

#[test]
fn strict_resolution_does_not_canonicalize_whitespace_padded_raw_accessions() {
    for raw_id in [
        format!(" {}", IDS[0]),
        format!("{} ", IDS[0]),
        format!("\t{}\n", IDS[0]),
    ] {
        let (_, registry) = registry(vec![motif_entry(&raw_id, "AC", &AC_COUNTS)]);
        let legacy = registry
            .resolve(IDS[0])
            .expect("legacy lookup still normalizes IDs");
        assert_eq!(legacy.id, IDS[0]);
        assert_eq!(legacy.matrix_counts, AC_COUNTS);
        assert_rejected(
            panel(&IDS[..1], TfbsScoreTrackValueKind::LlrBits),
            &registry,
            "expected one exact registry entry",
        );
    }
}

#[test]
fn structured_source_metadata_preserves_legacy_loading_and_strict_resolution() {
    for source in [
        json!({"provider": "synthetic", "release": 1, "url": "https://example.invalid/matrices.json"}),
        json!([{"provider": "synthetic", "local": true}]),
    ] {
        let registry_bytes = serde_json::to_string(&json!({
            "schema": "gentle.tf_motifs.v1",
            "source": source,
            "motifs": [motif_entry(IDS[0], "AC", &AC_COUNTS)],
        }))
        .unwrap();
        let registry = TfMotifDb::from_json_for_test(&registry_bytes)
            .expect("structured source metadata must not discard a valid legacy registry");
        let legacy = registry
            .resolve(FACTOR)
            .expect("legacy factor lookup remains available");
        assert_eq!(legacy.id, IDS[0]);
        assert_eq!(legacy.matrix_counts, AC_COUNTS);
        let resolution = resolve(
            panel(&IDS[..1], TfbsScoreTrackValueKind::LlrBits),
            &registry,
        );
        assert_eq!(resolution.matrices.len(), 1);
        assert_eq!(resolution.matrices[0].specification.source_id, IDS[0]);
        assert_eq!(resolution.matrices[0].matrix_counts, AC_COUNTS);
        assert!(
            resolution.registry_source_url.is_none(),
            "do not infer a URL from structured metadata"
        );
        let active = resolution
            .registry_sources
            .iter()
            .find(|binding| binding.role == "active_registry")
            .expect("retain the actual registry byte binding");
        assert_eq!(active.sha256, sha256_hex_bytes(registry_bytes.as_bytes()));
    }
}

#[test]
fn strict_resolution_rejects_consensus_fallback_and_malformed_full_pfms() {
    let mut consensus_only = motif_entry(IDS[0], "AC", &AC_COUNTS);
    consensus_only.as_object_mut().unwrap().remove("pfm");
    let mut invalid_entries = vec![consensus_only];
    for pfm in [
        json!({"a": [6.0], "c": [2.0, 6.0], "g": [1.0, 2.0], "t": [1.0, 1.0]}),
        json!({"a": [], "c": [], "g": [], "t": []}),
        json!({"a": [-1.0, 1.0], "c": [2.0, 6.0], "g": [1.0, 2.0], "t": [1.0, 1.0]}),
        json!({"a": [0.0, 1.0], "c": [0.0, 6.0], "g": [0.0, 2.0], "t": [0.0, 1.0]}),
        json!({"a": [1e308, 1.0], "c": [1e308, 6.0], "g": [1.0, 2.0], "t": [1.0, 1.0]}),
    ] {
        let mut entry = motif_entry(IDS[0], "AC", &AC_COUNTS);
        entry["pfm"] = pfm;
        invalid_entries.push(entry);
    }
    for entry in invalid_entries {
        let (_, registry) = registry(vec![entry, motif_entry(IDS[1], "GT", &GT_COUNTS)]);
        assert_rejected(
            panel(&IDS[..1], TfbsScoreTrackValueKind::LlrBits),
            &registry,
            "valid full PFM required",
        );
    }
}

#[test]
fn strict_resolution_rejects_mixed_or_unsupported_scores_at_the_resolver_boundary() {
    let (_, registry) = three_matrix_registry();
    let mut mixed = panel(&IDS[..2], TfbsScoreTrackValueKind::LlrBits);
    mixed.factors[1].score_kind = Some("true_log_odds_bits".into());
    assert_rejected(mixed, &registry, IDS[1]);

    let mut unsupported = panel(&IDS[..1], TfbsScoreTrackValueKind::LlrBits);
    unsupported.score_kind = "binding_affinity".into();
    assert_rejected(unsupported, &registry, "binding_affinity");
}

#[test]
fn plus_and_minus_profiles_match_an_independent_full_pfm_hand_calculation() {
    let (_, registry) = registry(vec![motif_entry(IDS[0], "AC", &AC_COUNTS)]);
    let engine = engine_with_sentinel();
    let before = serde_json::to_value(&engine).unwrap();
    for kind in [
        TfbsScoreTrackValueKind::LlrBits,
        TfbsScoreTrackValueKind::TrueLogOddsBits,
    ] {
        let resolution = resolve(panel(&IDS[..1], kind), &registry);
        let report = engine
            .compute_verified_tss_profiles(
                bundle(vec![
                    record("plus", TssStrand::Plus, 1_000, "ACGTACTAC", false),
                    record("minus", TssStrand::Minus, 2_000, "ACGTACTAC", false),
                ]),
                resolution,
                &mut |_| true,
            )
            .unwrap();
        // Each column totals 10; the documented pseudocount is 10 * 1e-9
        // per base. These expressions do not call any scoring helper.
        let weight = |count: f64| {
            let p = (count + 1e-8) / (10.0 + 4e-8);
            match kind {
                TfbsScoreTrackValueKind::LlrBits => (4.0 * p).log2(),
                TfbsScoreTrackValueKind::TrueLogOddsBits => (3.0 * p / (1.0 - p)).log2(),
                _ => unreachable!(),
            }
        };
        let high = 2.0 * weight(6.0);
        let low = 2.0 * weight(1.0);
        let expected_forward = [
            high,
            2.0 * weight(2.0),
            low,
            low,
            high,
            weight(2.0) + weight(1.0),
            low,
            high,
        ];
        let expected_reverse = [
            low,
            2.0 * weight(2.0),
            high,
            low,
            low,
            weight(6.0) + weight(2.0),
            low,
            low,
        ];
        assert_eq!(report.windows.len(), 2);
        for window in &report.windows {
            assert!(!window.selected);
            assert!(window.comparisons.is_empty());
            let track = &window.tracks[0];
            assert_eq!(track.motif_length_bp, 2);
            assert_eq!(track.forward_scores.len(), 8);
            assert_eq!(track.reverse_scores.len(), 8);
            for (actual, expected) in track.forward_scores.iter().zip(expected_forward) {
                assert_close(actual.expect("valid forward window"), expected);
            }
            for (actual, expected) in track.reverse_scores.iter().zip(expected_reverse) {
                assert_close(actual.expect("valid reverse window"), expected);
            }
            assert_eq!(
                track.forward_maximum.as_ref().unwrap().local_start_0based,
                0
            );
            assert_eq!(
                track.reverse_maximum.as_ref().unwrap().local_start_0based,
                2
            );
            assert_close(track.forward_maximum.as_ref().unwrap().score, high);
            assert_close(track.reverse_maximum.as_ref().unwrap().score, high);
            let geometry = &window.record.geometry;
            assert_eq!(geometry.relative_at(0), Some(-3));
            assert_eq!(geometry.relative_at(3), Some(0));
            assert_eq!(geometry.relative_at(7), Some(4));
            let coordinates = match geometry.strand {
                TssStrand::Plus => [997, 1_000, 1_004],
                TssStrand::Minus => [2_003, 2_000, 1_996],
            };
            for (index, coordinate) in [0, 3, 7].into_iter().zip(coordinates) {
                assert_eq!(geometry.genomic_at(index), Some(coordinate));
            }
            assert_eq!(geometry.genomic_at(9), None);
        }
        assert_eq!(
            report.windows[0].tracks[0].forward_scores,
            report.windows[1].tracks[0].forward_scores
        );
        assert_eq!(
            report.windows[0].tracks[0].reverse_scores,
            report.windows[1].tracks[0].reverse_scores
        );
    }
    assert_eq!(serde_json::to_value(&engine).unwrap(), before);
}

#[test]
fn profiles_keep_three_matrix_order_memberships_and_selected_first_ties_per_gene() {
    let (_, registry) = three_matrix_registry();
    let resolution = resolve(panel(&IDS, TfbsScoreTrackValueKind::LlrBits), &registry);
    let mut other_gene_selected =
        record("gene-b-selected", TssStrand::Plus, 500, "ACGTACTAC", true);
    other_gene_selected.0.gene_id = "synthetic_gene_b".into();
    other_gene_selected.0.gene_symbol = "SyntheticGeneB".into();
    let mut chromosome_two = record("chromosome-two", TssStrand::Plus, 100, "ACGTACTAC", true);
    chromosome_two.0.geometry.chromosome = "synthetic_chr2".into();
    let records = vec![
        other_gene_selected,
        chromosome_two,
        record("unselected", TssStrand::Plus, 100, "ACGTACTAC", false),
        record("tie-b", TssStrand::Plus, 400, "ACGTACTAC", true),
        record("tie-a", TssStrand::Minus, 402, "ACGTACTAC", true),
        record("selected-earlier", TssStrand::Plus, 300, "ACGTACTAC", true),
    ];
    let expected_records = records
        .iter()
        .map(|(record, _, selected)| (record.promoter_id.clone(), (record.clone(), *selected)))
        .collect::<BTreeMap<_, _>>();
    let bundle = bundle(records);
    let expected_reference = bundle.reference.clone();
    let expected_warnings = bundle.warnings.clone();
    let expected_inputs = serde_json::to_value(&bundle.inputs).unwrap();
    let engine = engine_with_sentinel();
    let before = serde_json::to_value(&engine).unwrap();
    let mut background_starts = vec![];
    let mut background_stage = None;
    let mut last_progress = None;
    let report = engine
        .compute_verified_tss_profiles(bundle, resolution, &mut |progress| {
            if let OperationProgress::Tfbs(progress) = progress {
                if progress.stage_label.as_deref() == Some("background calibration") {
                    // Setup and the scan both emit zero; count stage entries,
                    // not callbacks, to detect recalibration between TSSs.
                    if background_stage.as_ref() != Some(&progress.motif_id) {
                        assert_eq!(progress.scanned_steps, 0);
                        background_starts.push(progress.motif_id.clone());
                        background_stage = Some(progress.motif_id.clone());
                    }
                } else {
                    background_stage = None;
                }
                last_progress = Some(progress);
            }
            true
        })
        .unwrap();
    assert_eq!(
        report
            .windows
            .iter()
            .map(|w| w.record.promoter_id.as_str())
            .collect::<Vec<_>>(),
        [
            "selected-earlier",
            "tie-a",
            "tie-b",
            "chromosome-two",
            "unselected",
            "gene-b-selected"
        ]
    );
    assert_eq!(background_starts, IDS);
    assert_close(last_progress.unwrap().total_percent, 100.0);
    for window in &report.windows {
        let (expected_record, expected_selected) = &expected_records[&window.record.promoter_id];
        assert_eq!(&window.record, expected_record);
        assert_eq!(window.selected, *expected_selected);
        assert_eq!(
            window
                .tracks
                .iter()
                .map(|t| t.accession.as_str())
                .collect::<Vec<_>>(),
            IDS
        );
        assert_eq!(
            window
                .tracks
                .iter()
                .map(|track| track.motif_length_bp)
                .collect::<Vec<_>>(),
            [2, 2, 3]
        );
        assert_ne!(
            window.tracks[0].forward_scores,
            window.tracks[1].forward_scores
        );
        assert_eq!(window.comparisons.len(), 6);
        for (pair, (left, right)) in [(IDS[0], IDS[1]), (IDS[0], IDS[2]), (IDS[1], IDS[2])]
            .into_iter()
            .enumerate()
        {
            for (offset, strand) in [TssStrand::Plus, TssStrand::Minus].into_iter().enumerate() {
                let comparison = &window.comparisons[pair * 2 + offset];
                assert_eq!(comparison.factor_id, FACTOR);
                assert_eq!(comparison.left_accession, left);
                assert_eq!(comparison.right_accession, right);
                assert_eq!(comparison.local_strand, strand);
                assert_eq!(
                    comparison.paired_window_count,
                    if pair == 0 { 8 } else { 7 }
                );
                assert_eq!(comparison.excluded_window_count, usize::from(pair != 0));
                assert!(comparison.undefined_reason.is_none());
                assert!(comparison.pearson.is_some_and(f64::is_finite));
                assert!(comparison.spearman.is_some_and(f64::is_finite));
            }
        }
    }
    assert_eq!(report.schema, REPORT_SCHEMA);
    assert_eq!(report.reference, expected_reference);
    assert_eq!(report.warnings, expected_warnings);
    assert_eq!(
        serde_json::to_value(&report.inputs[..1]).unwrap(),
        expected_inputs
    );
    assert_eq!(report.inputs.len(), 2);
    assert_eq!(report.inputs[1].role, "panel");
    assert_eq!(
        report.inputs[1].sha256,
        report.panel_resolution.panel_sha256
    );
    assert_eq!(report.non_claims, NON_CLAIMS);
    assert_eq!(serde_json::to_value(&engine).unwrap(), before);
}

#[test]
fn valid_profile_arrays_match_shared_summary_without_clipping_or_registry_override() {
    let _guard = tf_motifs::test_registry_lock().lock().unwrap();
    // "AC" is a real legacy registry alias; these longer DNA tokens are not.
    let literals = ["ACGT", "ACGTA"];
    let tokens = literals.map(str::to_owned);
    assert_eq!(
        GentleEngine::expand_tf_query_tokens(&tokens).unwrap(),
        tokens
    );
    let (_, registry) = registry(
        literals
            .iter()
            .enumerate()
            .map(|(index, literal)| {
                let counts = literal
                    .bytes()
                    .map(|base| {
                        [b'A', b'C', b'G', b'T']
                            .map(|candidate| f64::from(u8::from(base == candidate)))
                    })
                    .collect::<Vec<_>>();
                motif_entry(IDS[index], literal, &counts)
            })
            .collect(),
    );
    let engine = GentleEngine::new();
    let sequence = "ACGTACGTNACGTACGT";
    for kind in [
        TfbsScoreTrackValueKind::LlrBits,
        TfbsScoreTrackValueKind::LlrBackgroundTailLog10,
    ] {
        let report = engine
            .compute_verified_tss_profiles(
                bundle(vec![
                    record("plus", TssStrand::Plus, 1_000, sequence, false),
                    record("minus", TssStrand::Minus, 2_000, sequence, false),
                ]),
                resolve(panel(&IDS[..2], kind), &registry),
                &mut |_| true,
            )
            .unwrap();
        let summary = engine
            .summarize_tfbs_score_tracks(
                SequenceScanTarget::InlineSequence {
                    sequence_text: sequence.into(),
                    topology: InlineSequenceTopology::Linear,
                    id_hint: Some("synthetic-parity".into()),
                    span_start_0based: None,
                    span_end_0based_exclusive: None,
                },
                &tokens,
                kind,
                false,
            )
            .unwrap();
        assert!(!summary.clip_negative);
        assert_eq!(summary.tracks.len(), 2);
        for window in &report.windows {
            for (index, (track, shared)) in window.tracks.iter().zip(&summary.tracks).enumerate() {
                assert_eq!(track.accession, IDS[index]);
                assert_eq!(shared.tf_id, literals[index]);
                assert_eq!(track.motif_length_bp, shared.motif_length_bp);
                for (scores, expected) in [
                    (&track.forward_scores, &shared.forward_scores),
                    (&track.reverse_scores, &shared.reverse_scores),
                ] {
                    assert_eq!(scores.len(), expected.len());
                    for (offset, (value, shared_value)) in scores.iter().zip(expected).enumerate() {
                        let valid = !sequence.as_bytes()[offset..offset + track.motif_length_bp]
                            .contains(&b'N');
                        assert_eq!(value.is_some(), valid);
                        if let Some(value) = value {
                            assert_close(*value, *shared_value);
                        }
                    }
                    assert!(scores.iter().any(Option::is_some));
                    assert!(scores.iter().any(Option::is_none));
                }
            }
        }
    }
}

#[test]
fn ambiguous_windows_are_null_and_trailing_or_too_short_windows_are_not_zero_padded() {
    let (_, registry) = registry(vec![
        motif_entry(IDS[0], "AC", &AC_COUNTS),
        motif_entry(IDS[1], "ACG", &ACG_COUNTS),
    ]);
    let report = GentleEngine::new()
        .compute_verified_tss_profiles(
            bundle(vec![
                record("partly-ambiguous", TssStrand::Plus, 100, "ACNTA", false),
                record("all-ambiguous", TssStrand::Minus, 200, "NNNNN", false),
                record("too-short", TssStrand::Plus, 300, "A", false),
            ]),
            resolve(
                panel(&IDS[..2], TfbsScoreTrackValueKind::LlrBits),
                &registry,
            ),
            &mut |_| true,
        )
        .unwrap();
    let partial = &report.windows[0];
    assert_eq!(partial.record.promoter_id, "partly-ambiguous");
    for scores in [
        &partial.tracks[0].forward_scores,
        &partial.tracks[0].reverse_scores,
    ] {
        assert_eq!(
            scores.iter().map(Option::is_some).collect::<Vec<_>>(),
            [true, false, false, true]
        );
        let encoded = serde_json::to_value(scores).unwrap();
        assert!(encoded[1].is_null());
        assert!(encoded[2].is_null());
    }
    for window in &report.windows {
        let sequence_length = window.record.geometry.length().unwrap();
        for track in &window.tracks {
            let expected_count = sequence_length
                .checked_sub(track.motif_length_bp)
                .map_or(0, |n| n + 1);
            assert_eq!(track.forward_scores.len(), expected_count);
            assert_eq!(track.reverse_scores.len(), expected_count);
            for offset in expected_count..sequence_length {
                assert_eq!(track.forward_scores.get(offset), None);
                assert_eq!(track.reverse_scores.get(offset), None);
            }
            if window.record.promoter_id != "partly-ambiguous" || track.motif_length_bp == 3 {
                assert!(track.forward_scores.iter().all(Option::is_none));
                assert!(track.reverse_scores.iter().all(Option::is_none));
                assert!(track.forward_maximum.is_none());
                assert!(track.reverse_maximum.is_none());
                assert!(track.forward_peaks.is_empty());
                assert!(track.reverse_peaks.is_empty());
                assert!(track.normalization_reference.is_null());
            }
        }
        for comparison in &window.comparisons {
            assert_eq!(comparison.paired_window_count, 0);
            assert_eq!(
                comparison.excluded_window_count,
                window.tracks[0].forward_scores.len()
            );
            assert_eq!(
                comparison.undefined_reason.as_deref(),
                Some("insufficient_common_valid_windows")
            );
            assert!(comparison.pearson.is_none());
            assert!(comparison.spearman.is_none());
        }
    }
}

#[test]
fn negative_maxima_and_more_than_three_peaks_survive_display_clipping() {
    let (_, registry) = registry(vec![motif_entry(IDS[0], "AC", &AC_COUNTS)]);
    let engine = GentleEngine::new();
    let mut raw_panel = panel(&IDS[..1], TfbsScoreTrackValueKind::LlrBits);
    raw_panel.top_hit_count = 5;
    let mut clipped_panel = raw_panel.clone();
    clipped_panel.clip_negative = true;
    let compute = |panel| {
        engine
            .compute_verified_tss_profiles(
                bundle(vec![record(
                    "negative",
                    TssStrand::Plus,
                    100,
                    "TTTTTTTTTTTT",
                    false,
                )]),
                resolve(panel, &registry),
                &mut |_| true,
            )
            .unwrap()
    };
    let raw = compute(raw_panel);
    let clipped = compute(clipped_panel);
    assert!(!raw.panel_resolution.panel.clip_negative);
    assert!(clipped.panel_resolution.panel.clip_negative);
    assert_eq!(
        serde_json::to_value(&raw.windows).unwrap(),
        serde_json::to_value(&clipped.windows).unwrap()
    );
    let track = &clipped.windows[0].tracks[0];
    let low = ((1.0_f64 + 1e-8) / (10.0 + 4e-8) * 4.0).log2();
    let high = ((6.0_f64 + 1e-8) / (10.0 + 4e-8) * 4.0).log2();
    for (scores, maximum, peaks, expected) in [
        (
            &track.forward_scores,
            &track.forward_maximum,
            &track.forward_peaks,
            2.0 * low,
        ),
        (
            &track.reverse_scores,
            &track.reverse_maximum,
            &track.reverse_peaks,
            high + low,
        ),
    ] {
        assert!(
            scores
                .iter()
                .all(|score| score.is_some_and(|score| score < 0.0))
        );
        let maximum = maximum
            .as_ref()
            .expect("negative maximum is still a maximum");
        assert_eq!(maximum.local_start_0based, 0);
        assert_close(maximum.score, expected);
        assert_eq!(
            peaks
                .iter()
                .map(|peak| peak.local_start_0based)
                .collect::<Vec<_>>(),
            [0, 2, 4, 6, 8]
        );
        for peak in peaks {
            assert_close(peak.score, expected);
        }
    }
}

fn comparison_track(
    accession: &str,
    motif_length_bp: usize,
    forward_scores: Vec<Option<f64>>,
    reverse_scores: Vec<Option<f64>>,
) -> TssProfileTrack {
    TssProfileTrack {
        accession: accession.into(),
        motif_length_bp,
        forward_scores,
        reverse_scores,
        forward_maximum: None,
        reverse_maximum: None,
        forward_peaks: vec![],
        reverse_peaks: vec![],
        normalization_reference: serde_json::Value::Null,
    }
}

fn comparison_window(tracks: Vec<TssProfileTrack>) -> TssProfileWindow {
    TssProfileWindow {
        detail_context: None,
        record: record("comparison", TssStrand::Minus, 100, "ACGTACG", false).0,
        selected: false,
        selection_evidence: None,
        tracks,
        comparisons: vec![],
    }
}

#[test]
fn comparisons_use_common_starts_same_strand_unclipped_values_and_average_tied_ranks() {
    let (_, registry) = registry(vec![
        motif_entry(IDS[0], "AC", &AC_COUNTS),
        motif_entry(IDS[1], "ACG", &ACG_COUNTS),
    ]);
    let resolution = resolve(
        panel(&IDS[..2], TfbsScoreTrackValueKind::LlrBits),
        &registry,
    );
    let window = comparison_window(vec![
        comparison_track(
            IDS[0],
            2,
            vec![
                Some(-3.0),
                Some(-1.0),
                None,
                Some(1.0),
                Some(3.0),
                Some(99.0),
            ],
            vec![Some(1.0), Some(1.0), Some(3.0), Some(3.0), None, Some(99.0)],
        ),
        comparison_track(
            IDS[1],
            3,
            vec![Some(3.0), Some(1.0), Some(42.0), Some(-1.0), Some(-3.0)],
            vec![Some(2.0), Some(4.0), Some(8.0), Some(16.0), Some(55.0)],
        ),
    ]);
    let comparisons = GentleEngine::tss_matrix_comparisons(&window, &resolution);
    assert_eq!(comparisons.len(), 2);
    let forward = &comparisons[0];
    assert_eq!(forward.local_strand, TssStrand::Plus);
    assert_close(forward.pearson.unwrap(), -1.0);
    assert_close(forward.spearman.unwrap(), -1.0);
    let reverse = &comparisons[1];
    assert_eq!(reverse.local_strand, TssStrand::Minus);
    // Pearson: covariance 18, variances 4 and 115. Spearman ranks the
    // tied left signal as [1.5, 1.5, 3.5, 3.5], not ordinal ranks.
    assert_close(reverse.pearson.unwrap(), 9.0 / 115.0_f64.sqrt());
    assert_close(reverse.spearman.unwrap(), 2.0 / 5.0_f64.sqrt());
    for comparison in &comparisons {
        assert_eq!(comparison.paired_window_count, 4);
        assert_eq!(comparison.excluded_window_count, 2);
        assert!(comparison.undefined_reason.is_none());
        assert!(comparison.method.contains("no smoothing"));
        assert!(comparison.method.contains("unclipped"));
    }
    let mut display_only = resolution.clone();
    display_only.panel.clip_negative = true;
    display_only.panel.factors[0].color_hint = Some("#cc1122".into());
    assert_eq!(
        serde_json::to_value(GentleEngine::tss_matrix_comparisons(&window, &display_only)).unwrap(),
        serde_json::to_value(comparisons).unwrap()
    );
}

#[test]
fn constant_and_insufficient_comparisons_are_undefined_not_numeric_zero() {
    let (_, registry) = registry(vec![
        motif_entry(
            IDS[0],
            "ACAC",
            &[AC_COUNTS[0], AC_COUNTS[1], AC_COUNTS[0], AC_COUNTS[1]],
        ),
        motif_entry(
            IDS[1],
            "GTGT",
            &[GT_COUNTS[0], GT_COUNTS[1], GT_COUNTS[0], GT_COUNTS[1]],
        ),
    ]);
    let resolution = resolve(
        panel(&IDS[..2], TfbsScoreTrackValueKind::LlrBits),
        &registry,
    );
    let variable = vec![Some(1.0), Some(2.0), Some(3.0), Some(4.0)];
    for (left, right, paired, reason) in [
        (vec![Some(0.0); 4], variable.clone(), 4, "constant_signal"),
        (variable.clone(), vec![Some(-2.0); 4], 4, "constant_signal"),
        (
            vec![Some(7.0), Some(7.0), None, Some(7.0)],
            variable.clone(),
            3,
            "constant_signal",
        ),
        (
            vec![None, Some(2.0), Some(3.0), None],
            vec![Some(1.0), None, Some(4.0), None],
            1,
            "insufficient_common_valid_windows",
        ),
        (
            vec![None; 4],
            variable,
            0,
            "insufficient_common_valid_windows",
        ),
    ] {
        let window = comparison_window(vec![
            comparison_track(IDS[0], 4, left.clone(), right.clone()),
            comparison_track(IDS[1], 4, right, left),
        ]);
        let comparisons = GentleEngine::tss_matrix_comparisons(&window, &resolution);
        assert_eq!(comparisons.len(), 2);
        for comparison in comparisons {
            assert_eq!(comparison.paired_window_count, paired);
            assert_eq!(comparison.excluded_window_count, 4 - paired);
            assert_eq!(comparison.undefined_reason.as_deref(), Some(reason));
            assert!(comparison.pearson.is_none());
            assert!(comparison.spearman.is_none());
            let encoded = serde_json::to_value(&comparison).unwrap();
            assert!(encoded["pearson"].is_null());
            assert!(encoded["spearman"].is_null());
        }
    }
}

#[test]
fn near_constant_identical_profiles_remain_perfectly_correlated() {
    let (_, registry) = registry(vec![
        motif_entry(IDS[0], "AAAAA", &[AC_COUNTS[0]; 5]),
        motif_entry(IDS[1], "GGGGG", &[GT_COUNTS[0]; 5]),
    ]);
    let resolution = resolve(
        panel(&IDS[..2], TfbsScoreTrackValueKind::LlrBits),
        &registry,
    );
    let tiny = vec![Some(0.0), Some(1e-10), Some(2e-10)];
    let window = comparison_window(vec![
        comparison_track(IDS[0], 5, tiny.clone(), tiny.clone()),
        comparison_track(IDS[1], 5, tiny.clone(), tiny),
    ]);
    let comparisons = GentleEngine::tss_matrix_comparisons(&window, &resolution);
    assert_eq!(comparisons.len(), 2);
    for comparison in comparisons {
        assert_eq!(comparison.paired_window_count, 3);
        assert_eq!(comparison.excluded_window_count, 0);
        assert!(comparison.undefined_reason.is_none());
        assert_close(
            comparison
                .pearson
                .expect("small amplitude is not a constant signal"),
            1.0,
        );
        assert_close(comparison.spearman.unwrap(), 1.0);
    }
}

#[test]
fn extreme_finite_identical_profiles_remain_perfectly_correlated() {
    let (_, registry) = registry(vec![
        motif_entry(IDS[0], "AAAAA", &[AC_COUNTS[0]; 5]),
        motif_entry(IDS[1], "GGGGG", &[GT_COUNTS[0]; 5]),
    ]);
    let resolution = resolve(
        panel(&IDS[..2], TfbsScoreTrackValueKind::LlrBits),
        &registry,
    );
    let extremes = vec![Some(f64::MIN), Some(0.0), Some(f64::MAX)];
    let window = comparison_window(vec![
        comparison_track(IDS[0], 5, extremes.clone(), extremes.clone()),
        comparison_track(IDS[1], 5, extremes.clone(), extremes),
    ]);
    let comparisons = GentleEngine::tss_matrix_comparisons(&window, &resolution);
    assert_eq!(comparisons.len(), 2);
    for comparison in comparisons {
        assert_eq!(comparison.paired_window_count, 3);
        assert_eq!(comparison.excluded_window_count, 0);
        assert!(comparison.undefined_reason.is_none());
        assert_close(
            comparison
                .pearson
                .expect("finite extreme range must not overflow normalization"),
            1.0,
        );
        assert_close(comparison.spearman.unwrap(), 1.0);
    }
}

#[test]
fn cancellation_at_background_or_partial_tss_stages_leaves_engine_unchanged() {
    let (_, registry) = three_matrix_registry();
    let engine = engine_with_sentinel();
    let before = serde_json::to_value(&engine).unwrap();
    let revisions = (
        engine.execution_revision,
        engine.mutation_revision,
        engine.structural_revision,
    );
    for cancellation_stage in [
        "background_setup",
        "background_scan",
        "tss_scan",
        "second_matrix",
    ] {
        let mut cancelled = false;
        let mut finished_first_matrix_scan = false;
        let error = engine
            .compute_verified_tss_profiles(
                bundle(vec![
                    record("first", TssStrand::Plus, 100, "ACGTACTAC", false),
                    record("second", TssStrand::Minus, 200, "ACGTACTAC", false),
                ]),
                resolve(
                    panel(&IDS[..2], TfbsScoreTrackValueKind::LlrBits),
                    &registry,
                ),
                &mut |progress| {
                    assert!(!cancelled, "progress must stop after cancellation");
                    if let OperationProgress::Tfbs(progress) = progress {
                        let background =
                            progress.stage_label.as_deref() == Some("background calibration");
                        let tss_scan = progress.stage_label.as_deref() == Some("TSS scan");
                        if tss_scan
                            && progress.motif_index == 1
                            && progress.scanned_steps == progress.total_steps
                        {
                            finished_first_matrix_scan = true;
                        }
                        cancelled = match cancellation_stage {
                            "background_setup" => background && progress.scanned_steps == 0,
                            "background_scan" => background && progress.scanned_steps > 0,
                            "tss_scan" => tss_scan && progress.scanned_steps > 0,
                            "second_matrix" => {
                                tss_scan && progress.motif_index == 2 && progress.scanned_steps > 0
                            }
                            _ => unreachable!(),
                        };
                    }
                    !cancelled
                },
            )
            .expect_err("cancellation must not return a partially populated report");
        assert!(cancelled, "did not reach {cancellation_stage}");
        assert!(error.message.to_ascii_lowercase().contains("cancel"));
        if cancellation_stage == "second_matrix" {
            assert!(
                finished_first_matrix_scan,
                "exercise cancellation after partial work"
            );
        }
        assert_eq!(serde_json::to_value(&engine).unwrap(), before);
        assert_eq!(
            (
                engine.execution_revision,
                engine.mutation_revision,
                engine.structural_revision
            ),
            revisions
        );
        assert!(engine.undo_stack.is_empty());
        assert!(engine.redo_stack.is_empty());
    }
}

fn exported_tsv_rows(text: &str) -> Vec<BTreeMap<&str, &str>> {
    let mut lines = text
        .lines()
        .filter(|line| !line.is_empty() && !line.starts_with('#'));
    let header = lines
        .next()
        .expect("TSV header after its preamble")
        .split('\t')
        .collect::<Vec<_>>();
    lines
        .map(|line| {
            let cells = line.split('\t').collect::<Vec<_>>();
            assert_eq!(cells.len(), header.len());
            header.iter().copied().zip(cells).collect()
        })
        .collect()
}

#[test]
fn tp73_maximum_tail_survives_shared_engine_and_exported_tss_profile() {
    // Actual accession/PFM from the committed JASPAR registry; maximizing and
    // ambiguous words in synthetic TSS geometry, not a human promoter fixture.
    let registry =
        TfMotifDb::from_json_for_test(include_str!("../../../assets/jaspar.motifs.json")).unwrap();
    let mut input_panel = panel(
        &["MA0861.2"],
        TfbsScoreTrackValueKind::LlrBackgroundTailLog10,
    );
    input_panel.factors[0].factor_id = "TP73".into();
    let sequence = "ACATGTCTGGACATGT";
    let resolution = resolve(input_panel, &registry);
    let (llr, _) = GentleEngine::prepare_scoring_matrices(&resolution.matrices[0].matrix_counts);
    assert!(
        (GentleEngine::score_matrix_window(sequence.as_bytes(), &llr).unwrap() - 19.543680326691)
            .abs()
            < 1e-11
    );
    let report = engine_with_sentinel()
        .compute_verified_tss_profiles(
            bundle(vec![
                record("tp73-maximum", TssStrand::Plus, 1000, sequence, true),
                record(
                    "tp73-ambiguous",
                    TssStrand::Minus,
                    2000,
                    "NCATGTCTGGACATGT",
                    false,
                ),
            ]),
            resolution,
            &mut |_| true,
        )
        .unwrap();
    let track = &report
        .windows
        .iter()
        .find(|w| w.record.promoter_id == "tp73-maximum")
        .unwrap()
        .tracks[0];
    assert!((track.forward_scores[0].unwrap() - 9.632959861247).abs() < 1e-11);
    assert!(
        (track.normalization_reference["observed_peak_modeled_tail_probability"]
            .as_f64()
            .unwrap()
            / 4.0_f64.powi(-16)
            - 1.0)
            .abs()
            < 1e-12
    );
    for window in &report.windows {
        for value in window.tracks[0]
            .forward_scores
            .iter()
            .chain(&window.tracks[0].reverse_scores)
            .flatten()
        {
            assert!(*value <= 16.0 * 4.0_f64.log10() + 1e-11);
        }
        if window.record.promoter_id == "tp73-ambiguous" {
            assert_eq!(window.tracks[0].forward_scores, [None]);
            assert_eq!(window.tracks[0].reverse_scores, [None]);
        }
    }
    assert_eq!(
        report.score_policy["modeled_tail_method"],
        GentleEngine::TFBS_MODELED_TAIL_METHOD
    );
    let scratch = tempfile::tempdir().unwrap();
    let directory = std::fs::canonicalize(scratch.path())
        .unwrap()
        .join("tp73-tail-export");
    let receipt = crate::tss_profile_export::export_tss_profiles(
        &report,
        &ExportTssProfilesRequest {
            context_manifest: None,
            output_dir: directory.to_string_lossy().into_owned(),
            rendering: TssProfileRenderOptions {
                scale_mode: Some(TssScaleMode::SharedAcrossTss),
                panels_per_page: 1,
            },
            formats: vec![TssExportFormat::Svg],
        },
    )
    .unwrap();
    crate::tss_profile_export::verify_tss_profile_receipt(&directory, &receipt).unwrap();
    assert_eq!(
        receipt.rendering.scale_mode,
        Some(TssScaleMode::SharedAcrossTss)
    );
    assert_eq!(
        std::fs::read(directory.join("report.json")).unwrap(),
        serde_json::to_vec(&report).unwrap()
    );
    for file in receipt.outputs.keys().filter(|name| name.ends_with(".svg")) {
        let svg = std::fs::read_to_string(directory.join(file)).unwrap();
        assert!(svg.contains("shared_across_tss"));
        assert!(svg.contains("9.633"));
        assert!(!svg.contains("legacy/unversioned background scoring"));
    }
}

#[test]
fn two_tss_three_matrix_producer_exports_preserve_scores_pairs_and_receipt_bindings() {
    let (registry_bytes, registry) = three_matrix_registry();
    let mut input_panel = panel(&IDS, TfbsScoreTrackValueKind::LlrBits);
    input_panel.clip_negative = true;
    let engine = engine_with_sentinel();
    let before = serde_json::to_value(&engine).unwrap();
    let report = engine
        .compute_verified_tss_profiles(
            bundle(vec![
                record("matrix-plus", TssStrand::Plus, 1_000, "ACGTACTAC", false),
                record("matrix-minus", TssStrand::Minus, 2_000, "GTACTACGT", true),
            ]),
            resolve(input_panel, &registry),
            &mut |_| true,
        )
        .expect("produce the two-TSS/three-matrix report from frozen synthetic PFMs");
    assert_eq!(report.windows.len(), 2);
    assert_eq!(report.panel_resolution.matrices.len(), 3);
    assert!(report.producer_executable_sha256.is_none());
    let source = report.source.as_ref().expect("synthetic manifest source");
    assert_eq!(source.schema, BUNDLE_SCHEMA);
    assert_eq!(source.manifest_sha256, report.inputs[0].sha256);
    assert_eq!(report.inputs[0].role, "bundle_manifest");
    assert_eq!(
        report
            .windows
            .iter()
            .map(|w| w.record.promoter_id.as_str())
            .collect::<Vec<_>>(),
        ["matrix-minus", "matrix-plus"]
    );
    for window in &report.windows {
        assert_eq!(window.record.geometry.length(), Some(9));
        assert_eq!(
            window
                .tracks
                .iter()
                .map(|track| track.accession.as_str())
                .collect::<Vec<_>>(),
            IDS
        );
        assert_eq!(window.comparisons.len(), 6);
    }

    let scratch = tempfile::tempdir().unwrap();
    let directory = std::fs::canonicalize(scratch.path())
        .unwrap()
        .join("three-matrix-export");
    let receipt = crate::tss_profile_export::export_tss_profiles(
        &report,
        &ExportTssProfilesRequest {
            context_manifest: None,
            output_dir: directory.to_string_lossy().into_owned(),
            rendering: TssProfileRenderOptions::default(),
            formats: vec![TssExportFormat::Svg],
        },
    )
    .expect("export the producer report without adapting scores or matrix digests");
    crate::tss_profile_export::verify_tss_profile_receipt(&directory, &receipt).unwrap();
    let loaded_receipt =
        crate::tss_profile_export::read_and_verify_tss_profile_receipt(&directory).unwrap();
    assert_eq!(loaded_receipt.outputs, receipt.outputs);
    assert_eq!(receipt.tss_count, 2);
    assert_eq!(receipt.page_count, 2);
    assert_eq!(receipt.outputs.len(), 10);
    assert!(!receipt.outputs.contains_key("receipt.json"));
    assert_eq!(std::fs::read_dir(&directory).unwrap().count(), 11);
    for matrix in &report.panel_resolution.matrices {
        assert!(receipt.inputs.iter().any(|binding| binding.name
            == format!("{}.json", matrix.specification.source_id)
            && binding.sha256 == matrix.matrix_sha256));
    }
    assert!(
        receipt
            .inputs
            .iter()
            .any(|binding| binding.role == "active_registry"
                && binding.sha256 == sha256_hex_bytes(registry_bytes.as_bytes()))
    );

    let report_bytes = std::fs::read(directory.join("report.json")).unwrap();
    assert_eq!(report_bytes, serde_json::to_vec(&report).unwrap());
    assert_eq!(sha256_hex_bytes(&report_bytes), receipt.report_sha256);
    let exported_json: serde_json::Value = serde_json::from_slice(&report_bytes).unwrap();
    let index: serde_json::Value =
        serde_json::from_slice(&std::fs::read(directory.join("index.json")).unwrap()).unwrap();
    assert_eq!(index["tss_count"], json!(2));
    assert_eq!(index["page_count"], json!(2));
    let genes = index["genes"].as_array().unwrap();
    assert_eq!(genes.len(), 1);
    let gene = &genes[0];
    assert_eq!(gene["promoter_ids"], json!(["matrix-minus", "matrix-plus"]));
    let gene_json: serde_json::Value = serde_json::from_slice(
        &std::fs::read(directory.join(gene["report"].as_str().unwrap())).unwrap(),
    )
    .unwrap();
    assert_eq!(
        gene_json, exported_json,
        "one-gene JSON retains both TSSs, all matrices and all pairs"
    );

    let scores_text =
        std::fs::read_to_string(directory.join(gene["scores"].as_str().unwrap())).unwrap();
    assert!(scores_text.contains(NON_CLAIMS));
    let score_rows = exported_tsv_rows(&scores_text);
    assert_eq!(score_rows.len(), 2 * 3 * 9 * 2);
    let tracks = report
        .windows
        .iter()
        .flat_map(|window| window.tracks.iter().map(move |track| (window, track)));
    for (rows, (window, track)) in score_rows.chunks_exact(18).zip(tracks) {
        for (offset, row) in rows.iter().enumerate() {
            let start = offset / 2;
            let (strand, scores) = if offset % 2 == 0 {
                (TssStrand::Plus, &track.forward_scores)
            } else {
                (TssStrand::Minus, &track.reverse_scores)
            };
            assert_eq!(row["promoter_id"], window.record.promoter_id);
            assert_eq!(row["accession"], track.accession);
            assert_eq!(row["factor_id"], FACTOR);
            assert_eq!(row["score_kind"], "llr_bits");
            assert_eq!(
                row["local_window_start_0based"].parse::<usize>().unwrap(),
                start
            );
            assert_eq!(
                row["motif_length_bp"].parse::<usize>().unwrap(),
                track.motif_length_bp
            );
            assert_eq!(row["local_motif_strand"], strand.as_str());
            let genomic_strand = if strand == TssStrand::Plus {
                window.record.geometry.strand
            } else {
                window.record.geometry.strand.opposite()
            };
            assert_eq!(row["genomic_motif_strand"], genomic_strand.as_str());
            assert_eq!(
                row["genomic_window_start_1based"].parse::<u64>().unwrap(),
                window.record.geometry.genomic_at(start).unwrap()
            );
            if let Some(raw) = scores.get(start).copied().flatten() {
                assert_eq!(row["availability"], "available");
                assert_close(row["raw_score"].parse().unwrap(), raw);
                assert_close(row["display_score"].parse().unwrap(), raw.max(0.0));
            } else {
                assert_eq!(row["availability"], "outside_window");
                assert_eq!(row["raw_score"], "null");
                assert_eq!(row["display_score"], "null");
            }
        }
    }

    let comparisons_text = std::fs::read_to_string(directory.join("comparisons.tsv")).unwrap();
    assert!(comparisons_text.contains(NON_CLAIMS));
    assert_eq!(
        comparisons_text,
        std::fs::read_to_string(directory.join(gene["comparisons"].as_str().unwrap()),).unwrap()
    );
    let comparison_rows = exported_tsv_rows(&comparisons_text);
    assert_eq!(comparison_rows.len(), 12);
    let pairs = [(IDS[0], IDS[1]), (IDS[0], IDS[2]), (IDS[1], IDS[2])];
    for (rows, window) in comparison_rows.chunks_exact(6).zip(&report.windows) {
        for (index, (row, comparison)) in rows.iter().zip(&window.comparisons).enumerate() {
            let (left, right) = pairs[index / 2];
            let strand = [TssStrand::Plus, TssStrand::Minus][index % 2];
            assert_eq!(row["promoter_id"], window.record.promoter_id);
            assert_eq!(row["factor_id"], FACTOR);
            assert_eq!(row["left_accession"], left);
            assert_eq!(row["right_accession"], right);
            assert_eq!(row["local_strand"], strand.as_str());
            assert_eq!(
                row["paired_window_count"].parse::<usize>().unwrap(),
                if index < 2 { 8 } else { 7 }
            );
            assert_eq!(
                row["excluded_window_count"].parse::<usize>().unwrap(),
                usize::from(index >= 2)
            );
            assert_eq!(row["input_values"], "raw_unclipped");
            assert_eq!(row["undefined_reason"], "null");
            assert_close(row["pearson"].parse().unwrap(), comparison.pearson.unwrap());
            assert_close(
                row["spearman"].parse().unwrap(),
                comparison.spearman.unwrap(),
            );
        }
    }

    let pages = gene["pages"].as_array().unwrap();
    assert_eq!(pages.len(), 2);
    for (page, window) in pages.iter().zip(&report.windows) {
        assert_eq!(page["promoter_ids"], json!([window.record.promoter_id]));
        let files = page["files"].as_array().unwrap();
        assert_eq!(files.len(), 1);
        let file = files[0].as_str().unwrap();
        assert!(file.ends_with(".svg"));
        let svg = std::fs::read_to_string(directory.join(file)).unwrap();
        let elements = svg::read(&svg)
            .unwrap()
            .filter_map(|event| match event {
                svg::parser::Event::Error(error) => panic!("invalid SVG: {error}"),
                svg::parser::Event::Tag(_, _, attributes) => Some(attributes),
                _ => None,
            })
            .collect::<Vec<_>>();
        let tagged = |role: &str| {
            elements
                .iter()
                .filter(|attributes| {
                    attributes
                        .get("data-role")
                        .is_some_and(|value| value.to_string() == role)
                })
                .collect::<Vec<_>>()
        };
        let panels = tagged("tss-panel");
        assert_eq!(panels.len(), 1);
        assert_eq!(
            panels[0]["data-promoter-id"].to_string(),
            window.record.promoter_id
        );
        assert_eq!(
            panels[0]["data-selected"].to_string(),
            window.selected.to_string()
        );
        let rows = tagged("matrix-row");
        assert_eq!(rows.len(), 3);
        for (row, matrix) in rows.iter().zip(&report.panel_resolution.matrices) {
            assert_eq!(
                row["data-accession"].to_string(),
                matrix.specification.source_id
            );
            assert_eq!(row["data-matrix-sha256"].to_string(), matrix.matrix_sha256);
            assert_eq!(
                row["data-display-order"].to_string(),
                matrix.specification.display_order.to_string()
            );
        }
        assert_eq!(tagged("pfm-logo").len(), 3);
        assert_eq!(tagged("matrix-comparison").len(), 6);
    }
    assert_eq!(serde_json::to_value(&engine).unwrap(), before);
}

fn assert_smoke_export(
    directory: &std::path::Path,
    receipt: &TssProfileReceipt,
    expected_report: &TssProfileReport,
) {
    crate::tss_profile_export::verify_tss_profile_receipt(directory, receipt)
        .expect("verify the exact published inventory and its report/input bindings");
    assert_eq!(receipt.schema, RECEIPT_SCHEMA);
    assert_eq!(receipt.tss_count, 2);
    assert_eq!(receipt.page_count, 2);
    assert_eq!(receipt.outputs.len(), 17);
    assert!(!receipt.outputs.contains_key("receipt.json"));
    assert_eq!(std::fs::read_dir(directory).unwrap().count(), 18);
    assert_eq!(receipt.non_claims, NON_CLAIMS);
    let report_bytes = std::fs::read(directory.join("report.json")).unwrap();
    assert_eq!(sha256_hex_bytes(&report_bytes), receipt.report_sha256);
    assert_eq!(
        report_bytes,
        serde_json::to_vec(expected_report).unwrap(),
        "export preserves exact report bytes"
    );
    assert_eq!(
        serde_json::from_slice::<serde_json::Value>(&report_bytes).unwrap(),
        serde_json::to_value(expected_report).unwrap(),
        "JSON parsing preserves every scientific value exactly"
    );

    for extension in ["svg", "png", "pdf"] {
        let files = receipt
            .outputs
            .keys()
            .filter(|name| name.ends_with(&format!(".{extension}")))
            .collect::<Vec<_>>();
        assert_eq!(files.len(), 2, "one {extension} page per TSS");
        for name in files {
            let bytes = std::fs::read(directory.join(name)).unwrap();
            assert!(!bytes.is_empty());
            match extension {
                "svg" => assert!(std::str::from_utf8(&bytes).unwrap().contains("<svg")),
                "png" => assert!(bytes.starts_with(b"\x89PNG\r\n\x1a\n")),
                "pdf" => assert!(bytes.starts_with(b"%PDF-1.4")),
                _ => unreachable!(),
            }
        }
    }
    let mut availability = BTreeMap::<String, usize>::new();
    for name in receipt
        .outputs
        .keys()
        .filter(|name| name.ends_with(".scores.tsv"))
    {
        let text = std::fs::read_to_string(directory.join(name)).unwrap();
        for row in exported_tsv_rows(&text) {
            let status = row["availability"];
            *availability.entry(status.into()).or_default() += 1;
            if status != "available" {
                assert_eq!(row["raw_score"], "null");
            }
        }
    }
    assert_eq!(
        availability,
        BTreeMap::from([
            ("available".into(), 14),
            ("unavailable".into(), 10),
            ("outside_window".into(), 20),
        ])
    );
}

#[test]
fn typed_compute_and_report_only_replay_export_all_formats_without_source_changes() {
    let _guard = tf_motifs::test_registry_lock().lock().unwrap();
    let fixture =
        std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("test_files/fixtures/tss_profiles");
    let source_files = [
        "manifest.json",
        "plus.fa",
        "minus.fa",
        "selection.json",
        "SHA256SUMS",
        "panel.json",
    ]
    .into_iter()
    .map(|name| (name, std::fs::read(fixture.join(name)).unwrap()))
    .collect::<BTreeMap<_, _>>();
    let input_copy = tempfile::tempdir().unwrap();
    for name in source_files.keys() {
        std::fs::copy(fixture.join(name), input_copy.path().join(name)).unwrap();
    }
    let input_path = input_copy.path().to_path_buf();
    let retain_parent = std::env::var_os("GENTLE_TEST_TSS_PROFILE_OUTPUT_DIR");
    let scratch = match &retain_parent {
        Some(parent) => tempfile::Builder::new()
            .prefix("gentle-tss-profile-smoke-")
            .tempdir_in(parent),
        None => tempfile::Builder::new()
            .prefix("gentle-tss-profile-smoke-")
            .tempdir(),
    }
    .expect("create unique output scratch directory");
    // Resolve macOS /var or /tmp aliases before the exporter's no-symlink checks.
    let output_root = std::fs::canonicalize(scratch.path()).unwrap();
    let computed_output = output_root.join("computed");
    let replay_output = output_root.join("replayed");
    let export_request = |path: &std::path::Path| ExportTssProfilesRequest {
        context_manifest: None,
        output_dir: path.to_string_lossy().into_owned(),
        rendering: TssProfileRenderOptions::default(),
        formats: vec![
            TssExportFormat::Svg,
            TssExportFormat::Png,
            TssExportFormat::Pdf,
        ],
    };
    let mut state = engine_with_sentinel().state;
    state.sequences.insert(
        "untouched-sequence".into(),
        DNAsequence::from_sequence("ACGTTGCA").unwrap(),
    );
    let mut engine = GentleEngine::from_state(state);
    let before_state = serde_json::to_value(&engine.state).unwrap();
    let before_revisions = (engine.mutation_revision, engine.structural_revision);
    let mut compute_stages = vec![];
    let computed = engine
        .apply_with_progress(
            Operation::ComputeTssTfbsProfiles {
                request: Box::new(ComputeTssProfilesRequest {
                    manifest: input_path
                        .join("manifest.json")
                        .to_string_lossy()
                        .into_owned(),
                    panel: input_path.join("panel.json").to_string_lossy().into_owned(),
                    fasta: vec!["minus.fa".into(), "plus.fa".into()],
                    selection: Some("selection.json".into()),
                    expected_genome_id: "synthetic-genome-v1".into(),
                    expected_assembly: Some("synthetic-assembly-v1".into()),
                    expected_annotation_release: Some("synthetic-annotation-v1".into()),
                    expected_dataset_id: None,
                }),
                export: Some(export_request(&computed_output)),
            },
            |progress| {
                if let OperationProgress::Tfbs(progress) = progress {
                    compute_stages.push(progress.stage_label.unwrap_or_default());
                }
                true
            },
        )
        .expect("typed computation must export its own returned report without digest adaptation");
    let report = computed
        .tss_tfbs_profiles
        .expect("typed computational report");
    let receipt = computed
        .tss_tfbs_profile_receipt
        .expect("typed compute/export receipt");
    for stage in ["background calibration", "TSS scan", "document export"] {
        assert!(
            compute_stages.iter().any(|seen| seen == stage),
            "missing {stage} progress"
        );
    }
    assert_eq!(report.windows.len(), 2);
    assert_eq!(report.panel_resolution.matrices.len(), 1);
    let matrix = &report.panel_resolution.matrices[0];
    assert_eq!(matrix.specification.source_id, "MA0004.1");
    assert_eq!(matrix.specification.factor_id, "Arnt");
    assert_eq!(matrix.matrix_counts.len(), 6);
    assert_eq!(
        matrix.matrix_sha256,
        sha256_hex_bytes(
            &serde_json::to_vec(&("MA0004.1", "Arnt", &matrix.matrix_counts,)).unwrap()
        )
    );
    assert_eq!(report.windows[0].record.promoter_id, "synthetic-minus");
    assert!(report.windows[0].selected);
    assert_eq!(report.windows[1].record.promoter_id, "synthetic-plus");
    assert!(!report.windows[1].selected);
    for window in &report.windows {
        assert_eq!(window.record.geometry.length(), Some(11));
        assert_eq!(window.tracks.len(), 1);
        assert_eq!(window.tracks[0].forward_scores.len(), 6);
        assert_eq!(window.tracks[0].reverse_scores.len(), 6);
    }
    let report_bytes = serde_json::to_vec(report.as_ref()).unwrap();
    assert_smoke_export(&computed_output, &receipt, &report);
    let source = report
        .source
        .as_ref()
        .expect("reader-bound manifest provenance");
    assert_eq!(
        source.manifest_sha256,
        sha256_hex_bytes(&source_files["manifest.json"])
    );
    gentle_engine::tss_profiles::validate_sha256(
        report
            .producer_executable_sha256
            .as_deref()
            .expect("public producer executable digest"),
    )
    .unwrap();
    assert_eq!(serde_json::to_value(&engine.state).unwrap(), before_state);
    for (name, expected) in &source_files {
        assert_eq!(std::fs::read(input_path.join(name)).unwrap(), *expected);
    }

    let replay: TssProfileReport =
        serde_json::from_slice(&std::fs::read(computed_output.join("report.json")).unwrap())
            .unwrap();
    assert_eq!(
        serde_json::to_vec(&replay).unwrap(),
        report_bytes,
        "typed JSON replay must retain every score exactly, not merely within a tolerance"
    );
    drop(input_copy);
    assert!(
        !input_path.exists(),
        "report replay cannot rely on its original input files"
    );
    let mut export_callbacks = 0;
    let replayed = engine
        .apply_with_progress(
            Operation::ExportTssTfbsProfiles {
                report: Box::new(replay),
                request: export_request(&replay_output),
            },
            |progress| {
                if let OperationProgress::Tfbs(progress) = progress {
                    assert_eq!(
                        progress.stage_label.as_deref(),
                        Some("document export"),
                        "report-only replay must not calibrate or score again"
                    );
                    export_callbacks += 1;
                }
                true
            },
        )
        .expect("typed report replay succeeds without source FASTA/panel files");
    assert!(export_callbacks > 0);
    let replay_receipt = replayed
        .tss_tfbs_profile_receipt
        .expect("report-only export receipt");
    assert_smoke_export(&replay_output, &replay_receipt, &report);
    assert_eq!(
        replay_receipt.outputs, receipt.outputs,
        "same report and renderer produce identical artifacts"
    );
    assert_eq!(serde_json::to_value(&engine.state).unwrap(), before_state);
    assert_eq!(
        (engine.mutation_revision, engine.structural_revision),
        before_revisions
    );
    assert!(engine.undo_stack.is_empty());
    assert!(engine.redo_stack.is_empty());
    for (name, expected) in &source_files {
        assert_eq!(std::fs::read(fixture.join(name)).unwrap(), *expected);
    }
    if retain_parent.is_some() {
        eprintln!(
            "Retained TSS profile smoke artifacts: {}",
            scratch.keep().display()
        );
    }
}
