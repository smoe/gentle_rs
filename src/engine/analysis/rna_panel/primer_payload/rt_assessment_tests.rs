//! Hand-crafted oligos and FASTA records recreated in private temporary directories.
//! These synthetic tests exercise prospective counts and real byte intake, not assay performance.

use super::*;
use crate::engine_shell::{execute_shell_command, parse_shell_line, quote_shell_arg};

fn rc(s: &str) -> String {
    String::from_utf8(GentleEngine::reverse_complement_bytes(s.as_bytes())).unwrap()
}

fn request(root: &Path, fasta: &str, tail: &str) -> AssessRtPrimerCandidatesRequest {
    let path = root.join("synthetic background's file.fa");
    fs::write(&path, fasta.as_bytes()).unwrap();
    let handle = "TTTTTT";
    AssessRtPrimerCandidatesRequest {
        schema: RT_PRIMER_ASSESSMENT_REQUEST_SCHEMA.into(),
        primers: vec![RnaLibraryPrimer {
            primer_id: "invented".into(),
            assay_label: "RT synthetic".into(),
            sequence: format!("{handle}{tail}"),
            components: vec![
                RnaPrimerComponent {
                    component_id: "handle".into(),
                    role: "common_handle".into(),
                    start: 0,
                    end: handle.len(),
                },
                RnaPrimerComponent {
                    component_id: "tail".into(),
                    role: "target_core".into(),
                    start: handle.len(),
                    end: handle.len() + tail.len(),
                },
            ],
            core_orientation: PrimerCoreOrientation::Antisense3prime,
            expected_targets: vec![],
            variant: "synthetic-v1".into(),
            source: "Hand-crafted synthetic; no biological reference".into(),
            source_sha256: None,
        }],
        references: vec![RtPrimerBackground {
            reference_id: "bg".into(),
            scope: RtPrimerBackgroundScope::CustomRna,
            source: "Hand-crafted synthetic FASTA, written by this test".into(),
            fasta: Some(PrimerPayloadFile {
                path: path.to_str().unwrap().into(),
                sha256: sha256_prefixed_str(fasta),
                bytes: fasta.len() as u64,
            }),
            unavailable_reason: None,
            circular_record_ids: vec![],
        }],
        intended_sites: vec![],
        budgets: RtPrimerAssessmentBudgets::default(),
    }
}

fn run(req: AssessRtPrimerCandidatesRequest) -> RtPrimerAssessmentReport {
    GentleEngine::new()
        .assess_rt_primer_candidates(req, &mut |_| true)
        .unwrap()
}

fn cell(report: &RtPrimerAssessmentReport, n: usize) -> &RtPrimerSeedAssessment {
    report.cells.iter().find(|c| c.seed_length == n).unwrap()
}

#[test]
fn rt_primer_assessment_positional_matches_are_not_terminal_runs() {
    let td = tempfile::tempdir().unwrap();
    let fasta = format!(
        ">internal\n{}\n>terminal\n{}\n",
        rc("TTGTACGA"),
        rc("ACGTACCC")
    );
    let report = run(request(td.path(), &fasta, "GGGGACGTACGA"));
    let c = cell(&report, 8);
    assert_eq!(c.all_sites.assessed_windows, 2);
    assert_eq!(c.all_sites.mismatch_histogram[2], 2); // ≥6/8 for both.
    assert_eq!(c.all_sites.terminal_run_histogram[6], 1);
    assert_eq!(c.all_sites.terminal_run_histogram[0], 1);
    assert_eq!(report.abundance_status, "not_assessed");
    assert_eq!(report.retrospective_status, "not_assessed");
}

#[test]
fn rt_primer_assessment_rna_and_genomic_orientation_are_explicit() {
    let td = tempfile::tempdir().unwrap();
    let mut req = request(td.path(), ">a\nACGTACGA\n", "GGGGACGTACGA");
    assert_eq!(
        cell(&run(req.clone()), 8).all_sites.mismatch_histogram[0],
        0
    );
    req.references[0].scope = RtPrimerBackgroundScope::Genome;
    let report = run(req);
    assert_eq!(cell(&report, 8).all_sites.mismatch_histogram[0], 1);
    assert_eq!(
        cell(&report, 8).exact_examples[0].template_strands,
        vec!['-']
    );
}

#[test]
fn rt_primer_assessment_palindromes_are_one_physical_site_with_two_orientations() {
    let td = tempfile::tempdir().unwrap();
    let mut req = request(td.path(), ">a\nACGTACGT\n", "GGGGACGTACGT");
    req.references[0].scope = RtPrimerBackgroundScope::Genome;
    let report = run(req);
    let c = cell(&report, 8);
    assert_eq!(c.all_sites.assessed_windows, 1);
    assert_eq!(c.all_sites.mismatch_histogram[0], 1);
    assert_eq!(c.exact_examples[0].template_strands, vec!['+', '-']);
}

#[test]
fn rt_primer_assessment_exclusions_verify_full_tail_and_only_remove_intended_orientation() {
    let td = tempfile::tempdir().unwrap();
    let tail = "ACGTACGTACGT";
    let mut req = request(td.path(), &format!(">a\n{}\n", rc(tail)), tail);
    req.references[0].scope = RtPrimerBackgroundScope::Genome;
    req.intended_sites.push(RtPrimerIntendedSite {
        primer_id: "invented".into(),
        reference_id: "bg".into(),
        record_id: "a".into(),
        tail_start_0based: 0,
        tail_end_0based_exclusive: 12,
        template_strand: '+',
    });
    let report = run(req.clone());
    let c = cell(&report, 12);
    assert_eq!(c.all_sites.mismatch_histogram[0], 1);
    assert_eq!(c.excluding_intended_orientations.mismatch_histogram[0], 1); // other strand remains.
    req.references[0].scope = RtPrimerBackgroundScope::CustomRna;
    assert_eq!(
        cell(&run(req.clone()), 12)
            .excluding_intended_orientations
            .assessed_windows,
        0
    );
    req.intended_sites[0].record_id = "missing".into();
    assert!(
        GentleEngine::new()
            .assess_rt_primer_candidates(req, &mut |_| true)
            .is_err()
    );
    let mut bad = report.request;
    bad.intended_sites[0].template_strand = '-';
    bad.intended_sites[0].tail_end_0based_exclusive = 13;
    assert!(
        GentleEngine::new()
            .assess_rt_primer_candidates(bad, &mut |_| true)
            .is_err()
    );
}

#[test]
fn rt_primer_assessment_nonpalindromic_intended_site_and_reverse_seed_coordinates() {
    let td = tempfile::tempdir().unwrap();
    let tail = "GGGGACGTACGA";
    for strand in ['+', '-'] {
        let seq = if strand == '+' { rc(tail) } else { tail.into() };
        let mut req = request(td.path(), &format!(">a\n{seq}\n"), tail);
        req.references[0].scope = RtPrimerBackgroundScope::Genome;
        req.intended_sites.push(RtPrimerIntendedSite {
            primer_id: "invented".into(),
            reference_id: "bg".into(),
            record_id: "a".into(),
            tail_start_0based: 0,
            tail_end_0based_exclusive: 12,
            template_strand: strand,
        });
        let report = run(req);
        assert_eq!(cell(&report, 12).all_sites.mismatch_histogram[0], 1);
        assert_eq!(
            cell(&report, 12)
                .excluding_intended_orientations
                .mismatch_histogram[0],
            0
        );
        assert_eq!(cell(&report, 8).all_sites.mismatch_histogram[0], 1);
        assert_eq!(
            cell(&report, 8)
                .excluding_intended_orientations
                .mismatch_histogram[0],
            0
        );
    }
}

#[test]
fn rt_primer_assessment_circular_windows_cross_origin_without_multilap_matches() {
    let td = tempfile::tempdir().unwrap();
    // Physical circular RNA is the rotation of the reverse complement of ACGTACGA.
    let mut req = request(td.path(), ">a\nACGTTCGT\n", "GGGGACGTACGA");
    let linear = run(req.clone());
    req.references[0].circular_record_ids.push("a".into());
    let circular = run(req.clone());
    assert_eq!(cell(&linear, 8).all_sites.assessed_windows, 1);
    assert_eq!(cell(&circular, 8).all_sites.assessed_windows, 8);
    assert_eq!(cell(&circular, 8).all_sites.mismatch_histogram[0], 1);
    assert!(cell(&circular, 8).exact_examples[0].end_0based_exclusive > 8);
    assert_eq!(cell(&circular, 12).records_shorter_than_seed, 1);
    req.references[0].circular_record_ids.push("absent".into());
    assert!(
        GentleEngine::new()
            .assess_rt_primer_candidates(req, &mut |_| true)
            .is_err()
    );
}

#[test]
fn rt_primer_assessment_missing_short_ambiguous_and_handle_only_are_not_clean_negatives() {
    let td = tempfile::tempdir().unwrap();
    let mut req = request(
        td.path(),
        ">short\nACG\n>ambiguous\nNNNNNNNN\n>handle_only\nAAAAAAAA\n",
        "GGGGACGTACGA",
    );
    let mut missing = req.references[0].clone();
    missing.reference_id = "missing_transcriptome".into();
    missing.scope = RtPrimerBackgroundScope::Transcriptome;
    missing.fasta = None;
    missing.unavailable_reason = Some("not supplied".into());
    req.references.push(missing);
    let report = run(req);
    let c = cell(&report, 8);
    assert_eq!(c.all_sites.ambiguous_windows, 1);
    assert_eq!(c.records_shorter_than_seed, 1);
    assert_eq!(c.all_sites.mismatch_histogram[0], 0);
    assert_eq!(report.references[1].status, "not_assessed");
    assert_eq!(report.cells[4].reason.as_deref(), Some("not supplied"));
    assert_eq!(cell(&report, 12).status, "not_assessed");
}

#[test]
fn rt_primer_assessment_counts_records_and_positions_without_pooling_reference_scopes() {
    let td = tempfile::tempdir().unwrap();
    let mut req = request(td.path(), ">a\nAAAAAAAAA\n>b\nAAAAAAAA\n", "GGGGTTTTTTTT");
    req.budgets.max_examples_per_cell = 1;
    let mut duplicate = req.references[0].clone();
    duplicate.reference_id = "overlap".into();
    req.references.push(duplicate);
    let report = run(req);
    let c = cell(&report, 8);
    assert_eq!(c.all_sites.mismatch_histogram[0], 3);
    assert_eq!(c.all_sites.exact_matching_records, 2);
    assert_eq!(c.exact_examples.len(), 1);
    assert!(c.exact_examples_truncated);
    assert_eq!(report.cells[5].all_sites.mismatch_histogram[0], 3);
}

#[test]
fn rt_primer_assessment_real_file_boundary_hash_newline_and_duplicate_record_guards() {
    let td = tempfile::tempdir().unwrap();
    for newline in ["\n", "\r\n"] {
        let fasta = format!(">a{newline}TCGTACGT{newline}");
        let req = request(td.path(), &fasta, "GGGGACGTACGA");
        assert_eq!(
            cell(&run(req.clone()), 8).all_sites.mismatch_histogram[0],
            1
        );
        // Same length and valid FASTA: only the byte digest can refuse this edit.
        fs::write(
            &req.references[0].fasta.as_ref().unwrap().path,
            fasta.replace("TCGT", "ACGT"),
        )
        .unwrap();
        assert!(
            GentleEngine::new()
                .assess_rt_primer_candidates(req, &mut |_| true)
                .is_err()
        );
    }
    for fasta in [
        ">a\nTCGTACGT\n>a\nTCGTACGT\n",
        ">a\n",
        "ACGT\n",
        ">a\nACGU\n",
    ] {
        let req = request(td.path(), fasta, "GGGGACGTACGA");
        assert!(
            GentleEngine::new()
                .assess_rt_primer_candidates(req, &mut |_| true)
                .is_err()
        );
    }
}

#[test]
fn rt_primer_assessment_budget_and_cancellation_refuse_without_partial_report() {
    let td = tempfile::tempdir().unwrap();
    let req = request(td.path(), ">a\nTCGTACGT\n", "GGGGACGTACGA");
    let mut engine = GentleEngine::new();
    let before = serde_json::to_value(engine.snapshot()).unwrap();
    assert!(
        engine
            .apply_with_progress(
                Operation::AssessRtPrimerCandidates {
                    request: Box::new(req.clone())
                },
                |_| false
            )
            .is_err()
    );
    assert_eq!(serde_json::to_value(engine.snapshot()).unwrap(), before);
    for which in 0..4 {
        let mut limited = req.clone();
        match which {
            0 => limited.budgets.max_base_comparisons = 1,
            1 => limited.budgets.max_record_bases = 1,
            2 => limited.budgets.max_reference_bases = 1,
            _ => limited.budgets.max_input_bytes = 1,
        }
        assert!(
            engine
                .assess_rt_primer_candidates(limited, &mut |_| true)
                .is_err()
        );
    }
    let mut polls = 0;
    assert!(
        engine
            .assess_rt_primer_candidates(req, &mut |_| {
                polls += 1;
                polls < 2
            })
            .is_err()
    );
    assert_eq!(polls, 2); // Final callback cancellation cannot release a completed report.
}

#[test]
fn rt_primer_assessment_partition_validation_identity_and_legacy_panel_bytes() {
    let td = tempfile::tempdir().unwrap();
    let req = request(td.path(), ">a\nTCGTACGT\n", "GGGGACGTACGA");
    let panel = RnaPrimerPanel {
        schema: RNA_PRIMER_PANEL_SCHEMA.into(),
        content_sha256: String::new(),
        primers: req.primers.clone(),
        policy: RnaPrimerTracePolicy::default(),
    };
    let panel = GentleEngine::validate_rna_primer_panel(
        &RnaPrimerPanelBinding::Inline {
            panel: Box::new(panel),
        },
        true,
    )
    .unwrap();
    let bytes = encode(&panel).unwrap();
    let report = run(req.clone());
    assert_eq!(encode(&panel).unwrap(), bytes);
    let mut basis = report.clone();
    basis.content_sha256.clear();
    assert_eq!(report.content_sha256, hash(&basis).unwrap());
    assert_eq!(report.request_sha256, hash(&req).unwrap());
    assert_eq!(report.content_sha256, run(req.clone()).content_sha256);
    for which in 0..4 {
        let mut invalid = req.clone();
        match which {
            0 => invalid.primers[0].components[1].end -= 1,
            1 => invalid.primers[0].core_orientation = PrimerCoreOrientation::Sense5prime,
            2 => invalid.primers.push(invalid.primers[0].clone()),
            _ => invalid.references[0].unavailable_reason = Some("conflict".into()),
        }
        assert!(
            GentleEngine::new()
                .assess_rt_primer_candidates(invalid, &mut |_| true)
                .is_err()
        );
    }
    let mut changed = req.clone();
    changed.primers[0].sequence.replace_range(..1, "A");
    assert_ne!(report.content_sha256, run(changed).content_sha256);
    let mut json = serde_json::to_value(&req).unwrap();
    json["unexpected"] = serde_json::json!(true);
    assert!(serde_json::from_value::<AssessRtPrimerCandidatesRequest>(json).is_err());
}

#[test]
fn rt_primer_assessment_typed_operation_and_shell_share_report_without_project_changes() {
    let td = tempfile::tempdir().unwrap();
    let req = request(td.path(), ">a\nTCGTACGT\n", "GGGGACGTACGA");
    let mut engine = GentleEngine::new();
    let revision = engine.mutation_revision();
    let undo = engine.undo_available();
    let before = serde_json::to_value(engine.snapshot()).unwrap();
    let result = engine
        .apply(Operation::AssessRtPrimerCandidates {
            request: Box::new(req.clone()),
        })
        .unwrap();
    assert_eq!(serde_json::to_value(engine.snapshot()).unwrap(), before);
    assert_eq!(engine.mutation_revision(), revision);
    assert_eq!(engine.undo_available(), undo);
    let journal = engine.journal_len();
    let expected = serde_json::to_value(result.rt_primer_assessment.unwrap()).unwrap();
    let line = format!(
        "primers assess-rt-candidates {}",
        quote_shell_arg(&serde_json::to_string(&req).unwrap())
    );
    let command = parse_shell_line(&line).unwrap();
    assert!(!command.is_state_mutating());
    let actual = execute_shell_command(&mut engine, &command).unwrap();
    assert!(!actual.state_changed);
    assert_eq!(engine.journal_len(), journal);
    assert_eq!(actual.output["result"]["rt_primer_assessment"], expected);
    assert_eq!(serde_json::to_value(engine.snapshot()).unwrap(), before);
    let path = td.path().join("request's file.json");
    fs::write(&path, serde_json::to_vec(&req).unwrap()).unwrap();
    let command = parse_shell_line(&format!(
        "primers assess-rt-candidates {}",
        quote_shell_arg(&format!("@{}", path.display()))
    ))
    .unwrap();
    assert_eq!(
        execute_shell_command(&mut engine, &command).unwrap().output["result"]["rt_primer_assessment"],
        expected
    );
}
