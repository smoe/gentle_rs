//! Literal synthetic SNV windows and simulated anchor metadata, recreated here.
//! No fixture is a fetched organism, verified haplotype, editor or laboratory result.

use super::*;
use crate::dna_sequence::DNAsequence;
use crate::engine_shell::{execute_shell_command, parse_shell_line, quote_shell_arg};
use gentle_protocol::{BiologicalContext, GenomicRegionReference};
use serde_json::json;

const SOURCE_HASH: &str = "sha256:e76797dcfb9f3a4639fdf013e8a72d0bd943a105bfc397387335c6540c496559";
const REVERSE_HASH: &str =
    "sha256:49c92ef4f8cf2d905248239e333475f727b30d296b1961f9c9a0dda7087da6c6";

fn request(reverse: bool) -> GenomeEditIntentRequest {
    GenomeEditIntentRequest {
        schema: GENOME_EDIT_INTENT_REQUEST_SCHEMA.into(),
        target: GenomeEditTarget::InlineSequence {
            sequence: if reverse { "GTACGT" } else { "ACGTAC" }.into(),
        },
        reference: GenomeEditReference::Synthetic {
            namespace: "synthetic:snv-tests".into(),
            reference_id: "six-base-window".into(),
        },
        anchor: GenomeEditAnchor {
            start_1based: 1001,
            end_1based: 1006,
            strand: if reverse {
                GenomicRegionStrand::Minus
            } else {
                GenomicRegionStrand::Plus
            },
        },
        expected_sequence_sha256: if reverse { REVERSE_HASH } else { SOURCE_HASH }.into(),
        design_haplotype: "reference".into(),
        context: None,
        intent: GenomeEditIntent {
            start_0based: 1001,
            end_0based_exclusive: 1002,
            deleted: "C".into(),
            inserted: "T".into(),
        },
    }
}

fn genomic(request: &mut GenomeEditIntentRequest) {
    request.reference = GenomeEditReference::Genomic {
        reference: GenomicRegionReference {
            taxon_id: Some(9606),
            assembly_name: "TEST_ASSEMBLY.1".into(),
            assembly_accession: Some("TEST_ASSEMBLY.1".into()),
            contig_name: "TEST_CONTIG.2".into(),
            contig_accession: Some("TEST_CONTIG.2".into()),
            ..Default::default()
        },
    };
}

fn loaded(reverse: bool, genomic_mode: bool) -> (GentleEngine, GenomeEditIntentRequest) {
    let mut engine = GentleEngine::new();
    engine.state_mut().sequences.insert(
        "synthetic_loaded".into(),
        DNAsequence::from_sequence(if reverse { "GTACGT" } else { "ACGTAC" }).unwrap(),
    );
    let mut req = request(reverse);
    req.target = GenomeEditTarget::LoadedSequence {
        seq_id: "synthetic_loaded".into(),
    };
    if genomic_mode {
        genomic(&mut req);
        engine.state_mut().metadata.insert(PROVENANCE_METADATA_KEY.into(), json!({
            "genome_extractions":[{
                "seq_id":"synthetic_loaded", "operation":"synthetic-test-anchor", "catalog_path":"not-opened.json",
                "genome_id":"TEST_ASSEMBLY.1", "chromosome":"TEST_CONTIG.2",
                "start_1based":1001,"end_1based":1006,"anchor_strand":if reverse {"-"} else {"+"},
                "anchor_verified":true,"recorded_at_unix_ms":123
            }]
        }));
    }
    (engine, req)
}

fn evaluate_unchanged(
    engine: &mut GentleEngine,
    req: GenomeEditIntentRequest,
) -> Result<GenomeEditIntentReport, EngineError> {
    let before = serde_json::to_vec(engine.state()).unwrap();
    let revision = (
        engine.mutation_revision(),
        engine.structural_revision(),
        engine.undo_available(),
        engine.redo_available(),
    );
    let result = engine.apply(Operation::EvaluateGenomeEditIntent {
        request: Box::new(req),
    });
    assert_eq!(serde_json::to_vec(engine.state()).unwrap(), before);
    assert_eq!(
        (
            engine.mutation_revision(),
            engine.structural_revision(),
            engine.undo_available(),
            engine.redo_available()
        ),
        revision
    );
    result.map(|result| {
        assert!(result.created_seq_ids.is_empty() && result.changed_seq_ids.is_empty());
        *result.genome_edit_intent.unwrap()
    })
}

#[test]
fn genome_edit_intent_preserves_unreconciled_state_on_success_and_refusal() {
    for reverse in [false, true] {
        for inline in [false, true] {
            let (mut engine, loaded_request) = loaded(reverse, false);
            let req = if inline {
                request(reverse)
            } else {
                loaded_request
            };
            assert!(engine.state().lineage.nodes.is_empty());
            assert!(engine.state().lineage.seq_to_node.is_empty());
            assert!(engine.state().container_state.containers.is_empty());
            assert!(
                engine
                    .state()
                    .container_state
                    .seq_to_latest_container
                    .is_empty()
            );
            let original = serde_json::to_vec(engine.state()).unwrap();

            // Independent controls prove neither repair would be a no-op here.
            let mut lineage_repaired = engine.clone();
            lineage_repaired.reconcile_lineage_nodes();
            assert!(
                lineage_repaired
                    .state()
                    .lineage
                    .seq_to_node
                    .contains_key("synthetic_loaded")
            );
            assert!(
                lineage_repaired
                    .state()
                    .container_state
                    .containers
                    .is_empty()
            );
            assert_ne!(
                serde_json::to_vec(lineage_repaired.state()).unwrap(),
                original
            );
            let mut containers_repaired = engine.clone();
            containers_repaired.reconcile_containers();
            assert!(
                containers_repaired
                    .state()
                    .container_state
                    .seq_to_latest_container
                    .contains_key("synthetic_loaded")
            );
            assert!(containers_repaired.state().lineage.nodes.is_empty());
            assert_ne!(
                serde_json::to_vec(containers_repaired.state()).unwrap(),
                original
            );

            let mut refused = req.clone();
            refused.expected_sequence_sha256 = format!("sha256:{}", "0".repeat(64));
            let error = evaluate_unchanged(&mut engine, refused).unwrap_err();
            assert!(error.cause_chain[0].ends_with("sequence_digest_mismatch"));
            let report = evaluate_unchanged(&mut engine, req).unwrap();
            assert_eq!(report.desired_genomic_forward_sequence, "ATGTAC");
            assert_eq!(serde_json::to_vec(engine.state()).unwrap(), original);
        }
    }
}

#[test]
fn genome_edit_intent_both_strands_reconstruct_exactly_once_without_persistence() {
    for reverse in [false, true] {
        let mut engine = GentleEngine::new();
        let report = evaluate_unchanged(&mut engine, request(reverse)).unwrap();
        assert_eq!(
            report.source_sequence_sha256,
            if reverse { REVERSE_HASH } else { SOURCE_HASH }
        );
        assert_eq!(report.genomic_forward_sequence, "ACGTAC");
        assert_eq!(report.desired_genomic_forward_sequence, "ATGTAC");
        assert_eq!(
            report.desired_local_sequence,
            if reverse { "GTACAT" } else { "ATGTAC" }
        );
        assert_eq!(
            (report.local_start_0based, report.local_end_0based_exclusive),
            if reverse { (4, 5) } else { (1, 2) }
        );
        assert_eq!(
            report.annotation_status,
            GenomeEditAnnotationStatus::AnnotationNotSupplied
        );
        assert_eq!(
            report.reference_verification,
            GenomeEditReferenceVerification::SyntheticDeclared
        );
        assert!(report.local_projection.is_none());
        assert_eq!(
            serde_json::from_slice::<GenomeEditIntentReport>(&serde_json::to_vec(&report).unwrap())
                .unwrap(),
            report
        );
    }
}

#[test]
fn genome_edit_intent_anchor_edges_and_all_snv_substitutions_are_not_editor_rules() {
    let cases = [
        (0, "A", "C", "CCGTAC"),
        (0, "A", "G", "GCGTAC"),
        (0, "A", "T", "TCGTAC"),
        (1, "C", "A", "AAGTAC"),
        (1, "C", "G", "AGGTAC"),
        (1, "C", "T", "ATGTAC"),
        (2, "G", "A", "ACATAC"),
        (2, "G", "C", "ACCTAC"),
        (2, "G", "T", "ACTTAC"),
        (3, "T", "A", "ACGAAC"),
        (3, "T", "C", "ACGCAC"),
        (3, "T", "G", "ACGGAC"),
        (5, "C", "G", "ACGTAG"),
    ];
    for reverse in [false, true] {
        for (offset, deleted, inserted, desired) in cases {
            let mut req = request(reverse);
            req.intent = GenomeEditIntent {
                start_0based: 1000 + offset,
                end_0based_exclusive: 1001 + offset,
                deleted: deleted.into(),
                inserted: inserted.into(),
            };
            let report = evaluate_unchanged(&mut GentleEngine::new(), req).unwrap();
            assert_eq!(report.desired_genomic_forward_sequence, desired);
            let local_start = if reverse { 5 - offset } else { offset };
            assert_eq!(
                (report.local_start_0based, report.local_end_0based_exclusive),
                (local_start, local_start + 1)
            );
        }
    }
}

#[test]
fn genome_edit_intent_large_genomic_origin_never_becomes_an_absolute_usize_index() {
    let mut req = request(false);
    req.anchor.start_1based = u64::MAX - 5;
    req.anchor.end_1based = u64::MAX;
    req.intent.start_0based = u64::MAX - 6;
    req.intent.end_0based_exclusive = u64::MAX - 5;
    req.intent.deleted = "A".into();
    req.intent.inserted = "G".into();
    let report = evaluate_unchanged(&mut GentleEngine::new(), req).unwrap();
    assert_eq!(
        (report.local_start_0based, report.local_end_0based_exclusive),
        (0, 1)
    );
    assert_eq!(report.desired_genomic_forward_sequence, "GCGTAC");
}

#[test]
fn genome_edit_intent_documented_operation_is_the_tested_synthetic_request() {
    let section = include_str!(concat!(env!("CARGO_MANIFEST_DIR"), "/docs/cli/snv-intent.md"))
        .split("## Read-Only SNV Intent")
        .nth(1)
        .unwrap();
    let json = section
        .split("```json")
        .nth(1)
        .unwrap()
        .split("```")
        .next()
        .unwrap();
    let operation: Operation = serde_json::from_str(json).unwrap();
    let Operation::EvaluateGenomeEditIntent {
        request: documented,
    } = operation
    else {
        panic!("documented SNV operation")
    };
    assert_eq!(*documented, request(false));
    assert_eq!(
        evaluate_unchanged(&mut GentleEngine::new(), *documented)
            .unwrap()
            .desired_local_sequence,
        "ATGTAC"
    );
}

#[test]
fn genome_edit_intent_bad_required_bindings_fail_closed() {
    let base = serde_json::to_value(request(false)).unwrap();
    let cases = [
        ("/schema", json!("wrong")),
        ("/expected_sequence_sha256", json!("")),
        ("/expected_sequence_sha256", json!("sha256:000")),
        (
            "/expected_sequence_sha256",
            json!(format!("sha256:{}", "0".repeat(64))),
        ),
        ("/design_haplotype", json!("")),
        ("/design_haplotype", json!("unknown")),
        ("/anchor/start_1based", json!(0)),
        ("/anchor/end_1based", json!(999)),
        ("/anchor/end_1based", json!(1007)),
        ("/anchor/strand", json!("unstranded")),
        ("/anchor/end_1based", json!(u64::MAX)),
        ("/intent/start_0based", json!(1003)),
        ("/intent/end_0based_exclusive", json!(1003)),
        ("/intent/deleted", json!("G")),
        ("/intent/inserted", json!("C")),
        ("/reference/namespace", json!("biological")),
    ];
    for (path, value) in cases {
        let mut raw = base.clone();
        *raw.pointer_mut(path).unwrap() = value;
        let error = evaluate_unchanged(
            &mut GentleEngine::new(),
            serde_json::from_value(raw).unwrap(),
        )
        .unwrap_err();
        assert_eq!(
            error.code,
            ErrorCode::InvalidInput,
            "{path}: {}",
            error.message
        );
        assert!(error.cause_chain[0].starts_with("genome_edit_intent."));
    }
    for (start, end, deleted) in [(999, 1000, "A"), (1006, 1007, "C")] {
        let mut req = request(false);
        req.intent.start_0based = start;
        req.intent.end_0based_exclusive = end;
        req.intent.deleted = deleted.into();
        assert!(
            evaluate_unchanged(&mut GentleEngine::new(), req)
                .unwrap_err()
                .cause_chain[0]
                .ends_with("outside_sequence")
        );
    }
}

#[test]
fn genome_edit_intent_ambiguity_is_independent_of_optional_guide_filters() {
    for symbol in b"NRYSWKMBDHVU-?" {
        for field in ["sequence", "deleted", "inserted"] {
            let mut req = request(false);
            let value = (*symbol as char).to_string();
            match field {
                "sequence" => {
                    req.target = GenomeEditTarget::InlineSequence {
                        sequence: format!("A{value}GTAC"),
                    }
                }
                "deleted" => req.intent.deleted = value,
                _ => req.intent.inserted = value,
            }
            let error = evaluate_unchanged(&mut GentleEngine::new(), req).unwrap_err();
            assert_eq!(error.code, ErrorCode::InvalidInput);
            assert!(error.cause_chain[0].ends_with("ambiguous_dna"));
        }
    }
    let mut req = request(false);
    req.target = GenomeEditTarget::InlineSequence {
        sequence: "acgtac".into(),
    };
    req.intent.deleted = "c".into();
    req.intent.inserted = "t".into();
    let report = evaluate_unchanged(&mut GentleEngine::new(), req).unwrap();
    assert_eq!(report.source_sequence_sha256, SOURCE_HASH);
    assert_eq!(report.source_sequence, "ACGTAC");
    assert_eq!(report.request.intent.deleted, "C");
    assert_eq!(report.request.intent.inserted, "T");
    assert_eq!(
        report.request.target,
        GenomeEditTarget::InlineSequence {
            sequence: "acgtac".into(),
        }
    );
}

#[test]
fn genome_edit_intent_admission_limits_and_legacy_defaults_do_not_authorize_a_binding() {
    for oversized in 0..3 {
        let mut req = request(false);
        match oversized {
            0 => {
                req.target = GenomeEditTarget::InlineSequence {
                    sequence: "A".repeat(MAX_GENOME_EDIT_WINDOW_BP + 1),
                }
            }
            1 => req.intent.inserted = "A".repeat(MAX_GENOME_EDIT_WINDOW_BP + 1),
            _ => req.design_haplotype = "x".repeat(MAX_GENOME_EDIT_ID_BYTES + 1),
        }
        assert_eq!(
            evaluate_unchanged(&mut GentleEngine::new(), req)
                .unwrap_err()
                .code,
            ErrorCode::InvalidInput
        );
    }
    // Preserve old serde behavior; the editing boundary must not trust it.
    let legacy: gentle_protocol::GenomicRegionLocalProjection = serde_json::from_str("{}").unwrap();
    assert_eq!(
        legacy.status,
        gentle_protocol::GenomicRegionLocalProjectionStatus::Current
    );
    let mut raw = serde_json::to_value(request(false)).unwrap();
    raw["local_projection"] = serde_json::to_value(legacy).unwrap();
    assert!(
        serde_json::from_value::<Operation>(json!({"EvaluateGenomeEditIntent":{"request":raw}}))
            .is_err()
    );
    let mut req = request(false);
    req.reference = GenomeEditReference::Genomic {
        reference: GenomicRegionReference::default(),
    };
    assert!(evaluate_unchanged(&mut GentleEngine::new(), req).is_err());
    let mut req = request(false);
    req.context = Some(BiologicalContext::default());
    assert!(evaluate_unchanged(&mut GentleEngine::new(), req).is_err());
    for field in ["assembly_accession", "contig_accession"] {
        let (mut engine, req) = loaded(false, true);
        let mut raw = serde_json::to_value(req).unwrap();
        raw["reference"]["reference"][field] = json!("OTHER_VERSION.9");
        assert!(evaluate_unchanged(&mut engine, serde_json::from_value(raw).unwrap()).is_err());
    }
}

#[test]
fn genome_edit_intent_indels_and_circular_dna_are_explicitly_unsupported() {
    for (start, end, deleted, inserted) in [
        (1003, 1003, "", "GG"),
        (1002, 1004, "GT", ""),
        (1002, 1004, "GT", "AA"),
    ] {
        let mut req = request(false);
        req.intent = GenomeEditIntent {
            start_0based: start,
            end_0based_exclusive: end,
            deleted: deleted.into(),
            inserted: inserted.into(),
        };
        assert_eq!(
            evaluate_unchanged(&mut GentleEngine::new(), req)
                .unwrap_err()
                .code,
            ErrorCode::Unsupported
        );
    }
    let (mut engine, req) = loaded(false, false);
    engine
        .state_mut()
        .sequences
        .get_mut("synthetic_loaded")
        .unwrap()
        .set_circular(true);
    assert_eq!(
        evaluate_unchanged(&mut engine, req).unwrap_err().code,
        ErrorCode::Unsupported
    );
    for molecule_type in ["RNA", "protein"] {
        let (mut engine, req) = loaded(false, false);
        engine
            .state_mut()
            .sequences
            .get_mut("synthetic_loaded")
            .unwrap()
            .set_molecule_type(molecule_type);
        assert_eq!(
            evaluate_unchanged(&mut engine, req).unwrap_err().code,
            ErrorCode::Unsupported,
            "{molecule_type}"
        );
    }
}

#[test]
fn genome_edit_intent_genomic_identity_and_annotation_unknowns_are_distinct() {
    let mut req = request(false);
    genomic(&mut req);
    let report = evaluate_unchanged(&mut GentleEngine::new(), req.clone()).unwrap();
    assert_eq!(
        report.reference_verification,
        GenomeEditReferenceVerification::InlineCallerDeclared
    );
    assert!(report.local_projection.is_none());
    req.context = Some(BiologicalContext {
        context_id: "fixture".into(),
        taxon_id: Some("9606".into()),
        assembly_accession: Some("TEST_ASSEMBLY.1".into()),
        annotation_source: Some("synthetic-annotations".into()),
        annotation_release: Some("1".into()),
        ..Default::default()
    });
    let report = evaluate_unchanged(&mut GentleEngine::new(), req.clone()).unwrap();
    assert_eq!(
        report.annotation_status,
        GenomeEditAnnotationStatus::SuppliedNotEvaluated
    );
    for field in ["annotation_release", "taxon_id", "assembly_accession"] {
        let mut raw = serde_json::to_value(&req).unwrap();
        raw["context"][field] = if field == "annotation_release" {
            json!(null)
        } else {
            json!("other")
        };
        assert_eq!(
            evaluate_unchanged(
                &mut GentleEngine::new(),
                serde_json::from_value(raw).unwrap()
            )
            .unwrap_err()
            .code,
            ErrorCode::InvalidInput
        );
    }
    for field in ["assembly_accession", "contig_accession", "taxon_id"] {
        let mut raw = serde_json::to_value(&req).unwrap();
        raw["reference"]["reference"][field] = json!(null);
        assert_eq!(
            evaluate_unchanged(
                &mut GentleEngine::new(),
                serde_json::from_value(raw).unwrap()
            )
            .unwrap_err()
            .code,
            ErrorCode::InvalidInput
        );
    }
    let mut synthetic = request(false);
    synthetic.context = req.context;
    assert!(evaluate_unchanged(&mut GentleEngine::new(), synthetic).is_err());
}

#[test]
fn genome_edit_intent_loaded_anchor_is_fresh_versioned_and_constructed_only_by_engine() {
    for reverse in [false, true] {
        let (mut engine, req) = loaded(reverse, true);
        let report = evaluate_unchanged(&mut engine, req.clone()).unwrap();
        assert_eq!(
            report.reference_verification,
            GenomeEditReferenceVerification::RetainedAnchorMatched
        );
        let projection = report.local_projection.unwrap();
        assert_eq!(projection.local_start_0based, if reverse { 4 } else { 1 });
        assert_eq!(
            projection.sequence_sha256,
            if reverse { REVERSE_HASH } else { SOURCE_HASH }
        );
        assert_eq!(
            projection.status,
            gentle_protocol::GenomicRegionLocalProjectionStatus::Current
        );
        let (mut lowercase, req) = loaded(reverse, true);
        lowercase.state_mut().sequences.insert(
            "synthetic_loaded".into(),
            DNAsequence::from_sequence(if reverse { "gtacgt" } else { "acgtac" }).unwrap(),
        );
        let report = evaluate_unchanged(&mut lowercase, req).unwrap();
        assert_eq!(
            report.source_sequence_sha256,
            if reverse { REVERSE_HASH } else { SOURCE_HASH }
        );
        assert_eq!(report.genomic_forward_sequence, "ACGTAC");
        assert_eq!(
            report.local_projection.unwrap().sequence_sha256,
            if reverse {
                "sha256:7a94e800a9ee2177bb3e8e46c59d03dacf62f2272457bb20400a9ee8ecc80405"
            } else {
                "sha256:f0d40db2ea6593e1297db2d5d96627b9de588913141b800d88eafdfcd657a680"
            }
        );
        for (field, value) in [
            ("anchor_verified", json!(false)),
            ("anchor_strand", json!(if reverse { "+" } else { "-" })),
            ("start_1based", json!(1002)),
            ("genome_id", json!("TEST_ASSEMBLY")),
            ("chromosome", json!("other")),
        ] {
            let (mut changed, req) = loaded(reverse, true);
            changed
                .state_mut()
                .metadata
                .get_mut(PROVENANCE_METADATA_KEY)
                .unwrap()["genome_extractions"][0][field] = value;
            assert_eq!(
                evaluate_unchanged(&mut changed, req).unwrap_err().code,
                ErrorCode::InvalidInput,
                "{field}"
            );
        }
        let (mut changed, req) = loaded(reverse, true);
        changed.state_mut().sequences.insert(
            "synthetic_loaded".into(),
            DNAsequence::from_sequence("AAAAAA").unwrap(),
        );
        assert!(
            evaluate_unchanged(&mut changed, req)
                .unwrap_err()
                .cause_chain[0]
                .ends_with("sequence_digest_mismatch")
        );
        let (mut engine, mut req) = loaded(reverse, true);
        req.reference = request(reverse).reference;
        assert!(
            evaluate_unchanged(&mut engine, req).is_err(),
            "no genomic-to-synthetic fallback"
        );
    }
    let (mut engine, req) = loaded(false, false);
    assert_eq!(
        evaluate_unchanged(&mut engine, req)
            .unwrap()
            .reference_verification,
        GenomeEditReferenceVerification::SyntheticDeclared
    );
    let (_, mut missing) = loaded(false, false);
    missing.target = GenomeEditTarget::LoadedSequence {
        seq_id: "absent".into(),
    };
    assert_eq!(
        evaluate_unchanged(&mut GentleEngine::new(), missing)
            .unwrap_err()
            .code,
        ErrorCode::NotFound
    );
}

#[test]
fn genome_edit_intent_latest_provenance_cannot_fall_back_or_mix_records() {
    let (mut engine, req) = loaded(false, true);
    let entries = engine
        .state_mut()
        .metadata
        .get_mut(PROVENANCE_METADATA_KEY)
        .unwrap()["genome_extractions"]
        .as_array_mut()
        .unwrap();
    let mut newer = entries[0].clone();
    newer["recorded_at_unix_ms"] = json!(124);
    newer["start_1based"] = json!(null);
    entries.push(newer);
    assert!(
        evaluate_unchanged(&mut engine, req).is_err(),
        "a newer incomplete anchor must not silently use an older complete one"
    );

    let (mut engine, mut req) = loaded(false, true);
    let entries = engine
        .state_mut()
        .metadata
        .get_mut(PROVENANCE_METADATA_KEY)
        .unwrap()["genome_extractions"]
        .as_array_mut()
        .unwrap();
    let mut conflicting = entries[0].clone();
    conflicting["start_1based"] = json!(1000);
    conflicting["end_1based"] = json!(1005);
    conflicting["chromosome"] = json!("OTHER_CONTIG.9");
    entries.push(conflicting);
    // The typed helper breaks equal timestamps by coordinates, while the
    // legacy summary uses the last record. A caller alias cannot join them.
    req.anchor.start_1based = 1000;
    req.anchor.end_1based = 1005;
    req.intent.start_0based = 1000;
    req.intent.end_0based_exclusive = 1001;
    let GenomeEditReference::Genomic { reference } = &mut req.reference else {
        unreachable!()
    };
    reference.contig_aliases.push("OTHER_CONTIG.9".into());
    assert!(evaluate_unchanged(&mut engine, req).is_err());
}

#[test]
fn genome_edit_intent_shell_parity_readiness_effects_and_existing_redo_are_preserved() {
    let mut engine = GentleEngine::new();
    engine
        .apply(Operation::CreateSequenceFromText {
            sequence_text: "AAAA".into(),
            output_id: Some("unrelated".into()),
            name: None,
            circular: false,
        })
        .unwrap();
    engine.undo_last_operation().unwrap();
    assert_eq!(engine.redo_available(), 1);
    let req = request(false);
    let expected = engine.evaluate_genome_edit_intent(req.clone()).unwrap();
    let command = parse_shell_line(&format!(
        "op {}",
        quote_shell_arg(&json!({"EvaluateGenomeEditIntent":{"request":req}}).to_string())
    ))
    .unwrap();
    assert!(!command.is_state_mutating());
    let before = serde_json::to_vec(engine.state()).unwrap();
    let revisions = (engine.mutation_revision(), engine.structural_revision());
    let result = execute_shell_command(&mut engine, &command).unwrap();
    assert!(!result.state_changed);
    assert_eq!(
        result.output["result"]["genome_edit_intent"],
        serde_json::to_value(expected).unwrap()
    );
    assert_eq!(serde_json::to_vec(engine.state()).unwrap(), before);
    assert_eq!(
        (engine.mutation_revision(), engine.structural_revision()),
        revisions
    );
    assert_eq!(engine.redo_available(), 1);
    for line in [
        "introspect capabilities",
        "introspect readiness EvaluateGenomeEditIntent",
        "introspect verify-effects EvaluateGenomeEditIntent",
    ] {
        let out = execute_shell_command(&mut engine, &parse_shell_line(line).unwrap()).unwrap();
        assert!(!out.state_changed);
        assert!(out.output.to_string().contains("EvaluateGenomeEditIntent"));
    }
    let descriptor = genome_edit_intent_descriptor();
    assert_eq!(descriptor["mutating"], "false");
    assert_eq!(descriptor["requires_confirmation"], false);
    assert_eq!(descriptor["effects"], json!([]));
    let mut legacy = result.output["result"].clone();
    legacy.as_object_mut().unwrap().remove("genome_edit_intent");
    assert!(
        serde_json::from_value::<OpResult>(legacy)
            .unwrap()
            .genome_edit_intent
            .is_none()
    );
}

fn genome_edit_intent_descriptor() -> serde_json::Value {
    let mut engine = GentleEngine::new();
    let result = execute_shell_command(
        &mut engine,
        &parse_shell_line("introspect capabilities").unwrap(),
    )
    .unwrap();
    // Find the descriptor without depending on unrelated introspection envelope keys.
    fn find(value: &serde_json::Value) -> Option<serde_json::Value> {
        if value.get("id").and_then(|v| v.as_str()) == Some("EvaluateGenomeEditIntent") {
            return Some(value.clone());
        }
        match value {
            serde_json::Value::Array(rows) => rows.iter().find_map(find),
            serde_json::Value::Object(rows) => rows.values().find_map(find),
            _ => None,
        }
    }
    find(&result.output).expect("typed operation descriptor")
}
