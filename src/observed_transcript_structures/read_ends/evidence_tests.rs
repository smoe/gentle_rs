//! Hand-crafted caller/unit fixtures, not actual tail calls or calibrated kit support.

use super::super::tests::fixture;
use super::*;

fn input() -> (ObservedSamAnalysis, TailImport) {
    let (_, mut analysis) = fixture("plus");
    analysis.units[0].rows.push(ObservedSamRowAssessment {
        ordinal: 1,
        qname: "synthetic_read".into(),
        read_group: None,
        flag: 0,
        target: Some("synthetic_contig".into()),
        position_1based: 101,
        reported_mapq: None,
        cigar: None,
        sequence_length: None,
        quality_length: None,
        axis_length: Some(25),
        clipping: None,
        unavailable: None,
    });
    let parameters = serde_json::json!({"synthetic_parameter":"explicit"});
    let configuration = CallerConfiguration {
        caller: "synthetic_caller".into(),
        version: "synthetic_version".into(),
        model: "synthetic_model".into(),
        kit_chemistry: "synthetic_kit".into(),
        molecule_type: "synthetic_molecule_type".into(),
        basecaller: "synthetic_basecaller".into(),
        basecaller_version: "synthetic_version".into(),
        basecaller_model: "synthetic_model".into(),
        parameter_sha256: crate::digest_utils::sha256_prefixed_str(
            &crate::digest_utils::canonical_json_string(&parameters).unwrap(),
        ),
        parameters,
    };
    let input = TailImport {
        schema: TAIL_SCHEMA.into(),
        raw_signal_manifest_sha256: format!("sha256:{}", "c".repeat(64)),
        configuration,
        rows: vec![TailRow {
            unit: analysis.units[0].unit.clone(),
            read_id: "synthetic_read".into(),
            tail_estimate_nt: Some(123.0),
            uncertainty_nt: Some(4.0),
            quality_score: Some(1.0),
            qc_status: QcStatus::Passed,
            qc_details: "synthetic declared QC".into(),
        }],
    };
    (analysis, input)
}

#[test]
fn read_end_evidence_configuration_and_qc_are_separate_from_estimates() {
    let (analysis, mut input) = input();
    let before = serde_json::to_value(&analysis).unwrap();
    let accepted = import_tail_estimates(
        &input,
        &[input.configuration.clone()],
        &input.raw_signal_manifest_sha256,
        &analysis,
        1,
    )
    .unwrap();
    assert_eq!(accepted[0].accepted_estimate_nt, Some(123.0));
    let unknown =
        import_tail_estimates(&input, &[], &input.raw_signal_manifest_sha256, &analysis, 1)
            .unwrap();
    assert_eq!(unknown[0].availability, TailAvailability::Unsupported);
    assert_eq!(unknown[0].accepted_estimate_nt, None);
    input.rows[0].qc_status = QcStatus::Failed;
    let failed = import_tail_estimates(
        &input,
        &[input.configuration.clone()],
        &input.raw_signal_manifest_sha256,
        &analysis,
        1,
    )
    .unwrap();
    assert_eq!(failed[0].availability, TailAvailability::RejectedQc);
    assert_eq!(failed[0].accepted_estimate_nt, None);
    assert_eq!(serde_json::to_value(&analysis).unwrap(), before);
}

#[test]
fn read_end_evidence_source_unit_read_and_parameter_mismatches_refuse() {
    let (analysis, input) = input();
    for change in 0..5 {
        let mut changed = input.clone();
        match change {
            0 => changed.rows[0].unit.source_sha256 = format!("sha256:{}", "e".repeat(64)),
            1 => changed.rows[0].read_id = "wrong_read".into(),
            2 => changed.configuration.parameter_sha256 = format!("sha256:{}", "e".repeat(64)),
            3 => changed.rows.push(changed.rows[0].clone()),
            _ => changed.rows[0].tail_estimate_nt = Some(f64::NAN),
        }
        assert!(
            import_tail_estimates(
                &changed,
                &[input.configuration.clone()],
                &input.raw_signal_manifest_sha256,
                &analysis,
                10
            )
            .is_err()
        );
    }
    assert!(
        import_tail_estimates(
            &input,
            &[],
            &format!("sha256:{}", "e".repeat(64)),
            &analysis,
            10
        )
        .is_err()
    );
}

#[test]
fn read_end_evidence_basecall_channel_is_explicitly_unavailable_without_scan() {
    let policy = CompositionPolicy {
        orientation: Provenance::Undeclared,
        trimming: Provenance::Undeclared,
        adapter_sequences: vec![],
        primer_sequences: vec![],
        window_bp: 32,
        minimum_run_bp: 7,
        maximum_mismatches: 0,
    };
    let evidence = unavailable_composition(policy.clone()).unwrap();
    assert_eq!(evidence.availability, "unavailable");
    assert_eq!(evidence.reason, "no_bound_basecalls");
    assert_eq!(evidence.diagnostics.len(), 2);
    let mut invalid = policy;
    invalid.window_bp = 0;
    assert!(unavailable_composition(invalid).is_err());
}
