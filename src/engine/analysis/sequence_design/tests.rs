//! Inline hand-crafted synthetic coding inserts, recreated from literal strings.
//! Used only for the shared engine/shell approval boundary, never as natural genes.

use super::*;

fn request() -> DnaSequenceDesignRequest {
    DnaSequenceDesignRequest {
        schema: REQUEST_SCHEMA.into(),
        target: DnaDesignTarget::InlineSequence {
            sequence: "ATGGAATTCTAA".into(),
        },
        purpose: "synthetic_coding_insert".into(),
        cds: DesignInterval {
            start_0based: 0,
            end_0based_exclusive: 12,
        },
        protein_sequence: "MEF".into(),
        genetic_code: 1,
        protected_intervals: vec![],
        avoid_motifs: vec![AvoidDesignMotif {
            pattern: "GAATTC".into(),
            strand: DesignStrand::Both,
        }],
        max_evaluations: 4096,
        search_strategy: DnaDesignSearchStrategy::FullEnumeration,
        output_seq_id: "synthetic_without_ecori".into(),
    }
}

#[test]
fn sequence_design_conflict_search_is_explicit_bound_and_separately_approved() {
    let mut engine = GentleEngine::new();
    let original = serde_json::to_value(engine.state()).unwrap();
    let legacy = plan(&mut engine, request());
    assert_eq!(legacy.algorithm, core::ALGORITHM);
    let legacy_json = serde_json::to_value(&legacy).unwrap();
    assert!(legacy_json["request"].get("search_strategy").is_none());
    let old_wire: DnaSequenceDesignReport = serde_json::from_value(legacy_json).unwrap();
    assert_eq!(
        GentleEngine::dna_design_approval(&old_wire).unwrap(),
        legacy.approval_digest.clone().unwrap()
    );

    let mut directed = request();
    directed.search_strategy = DnaDesignSearchStrategy::ConflictDirected;
    let preview = plan(&mut engine, directed.clone());
    assert_eq!(preview.algorithm, core::CONFLICT_ALGORITHM);
    assert_eq!(preview.output_sequence, legacy.output_sequence);
    assert!(preview.minimum_edits_proven && preview.optimization_complete);
    assert_eq!(preview, plan(&mut engine, directed));
    assert_ne!(preview.approval_digest, legacy.approval_digest);
    let roundtrip: DnaSequenceDesignReport =
        serde_json::from_value(serde_json::to_value(&preview).unwrap()).unwrap();
    assert_eq!(preview, roundtrip);

    for field in 0..3 {
        let mut altered = preview.clone();
        match field {
            0 => altered.algorithm = core::ALGORITHM.into(),
            1 => altered.request.search_strategy = DnaDesignSearchStrategy::FullEnumeration,
            2 => altered.nonclaims = legacy.nonclaims.clone(),
            _ => unreachable!(),
        }
        // A content hash is not a signature: semantic algorithm/policy validation
        // must reject an inconsistent client-rehashed preview too.
        altered.approval_digest = Some(GentleEngine::dna_design_approval(&altered).unwrap());
        assert!(apply(&mut engine, altered).is_err());
    }
    assert_eq!(serde_json::to_value(engine.state()).unwrap(), original);
    let result = apply(&mut engine, roundtrip).unwrap();
    assert_eq!(result.created_seq_ids, ["synthetic_without_ecori"]);
    assert_eq!(
        engine.state.sequences["synthetic_without_ecori"].get_forward_string(),
        "ATGGAATTTTAA"
    );
    engine.undo_last_operation().unwrap();
    assert!(
        !engine
            .state
            .sequences
            .contains_key("synthetic_without_ecori")
    );
    // The legacy preview/digest is still applicable after undo, without re-planning.
    apply(&mut engine, old_wire).unwrap();
}

#[test]
fn sequence_design_conflict_budget_and_cancellation_never_supply_approval() {
    let mut engine = GentleEngine::new();
    let mut directed = request();
    directed.search_strategy = DnaDesignSearchStrategy::ConflictDirected;
    directed.max_evaluations = 1;
    let exhausted = plan(&mut engine, directed.clone());
    assert_eq!(exhausted.algorithm, core::CONFLICT_ALGORITHM);
    assert_eq!(exhausted.status, DnaDesignStatus::SearchExhausted);
    assert!(exhausted.approval_digest.is_none());
    assert!(apply(&mut engine, exhausted).is_err());
    directed.max_evaluations = 4096;
    let cancelled = engine
        .plan_dna_sequence_design(directed, |evaluated| evaluated == 2)
        .unwrap();
    assert_eq!(cancelled.status, DnaDesignStatus::Cancelled);
    assert!(cancelled.approval_digest.is_none() && cancelled.output_sequence.is_none());
    assert!(apply(&mut engine, cancelled).is_err());
    assert!(engine.state.sequences.is_empty());
}

fn plan(engine: &mut GentleEngine, request: DnaSequenceDesignRequest) -> DnaSequenceDesignReport {
    *engine
        .apply(Operation::PlanDnaSequenceDesign {
            request: Box::new(request),
            path: None,
        })
        .unwrap()
        .dna_sequence_design
        .unwrap()
}

fn apply(
    engine: &mut GentleEngine,
    proposal: DnaSequenceDesignReport,
) -> Result<OpResult, EngineError> {
    let approval_digest = proposal.approval_digest.clone().unwrap_or_default();
    engine.apply(Operation::ApplyDnaSequenceDesign {
        proposal: Box::new(proposal),
        approval_digest,
    })
}

#[test]
fn sequence_design_engine_preview_exact_apply_and_undo() {
    let mut engine = GentleEngine::new();
    let original = serde_json::to_value(engine.state()).unwrap();
    let preview = plan(&mut engine, request());
    assert_eq!(preview.status, DnaDesignStatus::Feasible);
    assert_eq!(preview.output_sequence.as_deref(), Some("ATGGAATTTTAA"));
    assert_eq!(preview, plan(&mut engine, request()));
    assert_eq!(serde_json::to_value(engine.state()).unwrap(), original);
    let roundtrip: DnaSequenceDesignReport =
        serde_json::from_str(&serde_json::to_string(&preview).unwrap()).unwrap();
    assert_eq!(preview, roundtrip);
    let applied = apply(&mut engine, roundtrip).unwrap();
    let receipt = applied.dna_sequence_design_receipt.unwrap();
    assert_eq!(applied.created_seq_ids, [receipt.created_seq_id.clone()]);
    let dna = &engine.state.sequences[&receipt.created_seq_id];
    assert_eq!(
        dna.get_forward_string(),
        preview.output_sequence.as_deref().unwrap()
    );
    assert_eq!(dna.features().len(), 1);
    assert_eq!(&*dna.features()[0].kind, "CDS");
    assert!(
        engine
            .state
            .metadata
            .contains_key(&format!("dna_sequence_design:{}", receipt.created_seq_id))
    );
    assert!(apply(&mut engine, preview.clone()).is_err());
    engine.undo_last_operation().unwrap();
    assert!(!engine.state.sequences.contains_key(&receipt.created_seq_id));
    assert!(
        !engine
            .state
            .metadata
            .contains_key(&format!("dna_sequence_design:{}", receipt.created_seq_id))
    );
    engine.redo_last_operation().unwrap();
    assert_eq!(
        engine.state.sequences[&receipt.created_seq_id].get_forward_string(),
        "ATGGAATTTTAA"
    );
}

#[test]
fn sequence_design_rejects_tampering_stale_sources_and_preserves_lineage() {
    let mut engine = GentleEngine::new();
    engine
        .apply(Operation::CreateSequenceFromText {
            sequence_text: "ATGGAATTCTAA".into(),
            output_id: Some("synthetic_source".into()),
            name: None,
            circular: false,
        })
        .unwrap();
    let feature = gb_io::seq::Feature {
        kind: "misc_feature".into(),
        location: gb_io::seq::Location::simple_range(3, 9),
        qualifiers: vec![("note".into(), Some("Unverified source function".into()))],
    };
    engine
        .state
        .sequences
        .get_mut("synthetic_source")
        .unwrap()
        .features_mut()
        .push(feature.clone());
    let mut request = request();
    request.target = DnaDesignTarget::LoadedSequence {
        seq_id: "synthetic_source".into(),
    };
    let preview = plan(&mut engine, request);
    assert_eq!(preview.omitted_source_feature_count, 1);
    let original = serde_json::to_value(engine.state()).unwrap();
    let mut altered = preview.clone();
    altered.output_sequence = Some("ATGGAATTTTAG".into());
    assert!(apply(&mut engine, altered).is_err());
    let mut altered = preview.clone();
    altered.request.protected_intervals.push(DesignInterval {
        start_0based: 3,
        end_0based_exclusive: 9,
    });
    assert!(apply(&mut engine, altered).is_err());
    let mut altered = preview.clone();
    altered.genetic_code_mapping[0].residue = "X".into();
    assert!(apply(&mut engine, altered).is_err());
    assert_eq!(serde_json::to_value(engine.state()).unwrap(), original);
    engine
        .state
        .sequences
        .get_mut("synthetic_source")
        .unwrap()
        .features_mut()
        .push(feature);
    assert!(apply(&mut engine, preview.clone()).is_err());
    engine
        .state
        .sequences
        .get_mut("synthetic_source")
        .unwrap()
        .features_mut()
        .pop();
    let source = engine.state.sequences["synthetic_source"].clone();
    engine.state.sequences.insert(
        "synthetic_source".into(),
        DNAsequence::from_sequence("ATGGAATTTTAA").unwrap(),
    );
    assert!(apply(&mut engine, preview.clone()).is_err());
    engine
        .state
        .sequences
        .insert("synthetic_source".into(), source);
    apply(&mut engine, preview.clone()).unwrap();
    assert_eq!(
        engine.state.sequences["synthetic_source"].get_forward_string(),
        "ATGGAATTCTAA"
    );
    assert_eq!(
        engine.state.sequences["synthetic_without_ecori"]
            .features()
            .len(),
        1
    );
    let lineage = serde_json::to_value(engine.state()).unwrap()["lineage"].to_string();
    assert!(lineage.contains("synthetic_source"), "{lineage}");
}

#[test]
fn sequence_design_unresolved_cancelled_invalid_and_unsupported_have_no_approval() {
    let mut engine = GentleEngine::new();
    let mut bounded = request();
    bounded.max_evaluations = 1;
    let preview = plan(&mut engine, bounded);
    assert_eq!(preview.status, DnaDesignStatus::SearchExhausted);
    assert!(preview.approval_digest.is_none());
    assert!(apply(&mut engine, preview).is_err());
    let cancelled = engine
        .apply_with_progress(
            Operation::PlanDnaSequenceDesign {
                request: Box::new(request()),
                path: None,
            },
            |_| false,
        )
        .unwrap();
    let preview = *cancelled.dna_sequence_design.unwrap();
    assert_eq!(preview.status, DnaDesignStatus::Cancelled);
    assert!(preview.approval_digest.is_none() && preview.output_sequence.is_none());
    assert!(apply(&mut engine, preview).is_err());
    let mut invalid = request();
    invalid.protein_sequence = "MKF".into();
    assert_eq!(plan(&mut engine, invalid).status, DnaDesignStatus::Invalid);
    let mut unsupported = request();
    unsupported.genetic_code = 2;
    assert_eq!(
        plan(&mut engine, unsupported).status,
        DnaDesignStatus::Unsupported
    );
    let mut undeclared = request();
    undeclared.purpose = "natural_assay_target".into();
    assert!(
        engine
            .apply(Operation::PlanDnaSequenceDesign {
                request: Box::new(undeclared),
                path: None
            })
            .is_err()
    );
    let mut protected = request();
    protected.protected_intervals = vec![protected.cds.clone()];
    assert_eq!(
        plan(&mut engine, protected).status,
        DnaDesignStatus::ProvenInfeasible
    );
}

#[test]
fn sequence_design_shared_shell_discovery_and_json_file_parity() {
    assert_sequence_design_shell_parity(DnaDesignSearchStrategy::FullEnumeration);
}

#[test]
fn sequence_design_conflict_shared_shell_discovery_and_json_file_parity() {
    assert_sequence_design_shell_parity(DnaDesignSearchStrategy::ConflictDirected);
}

fn assert_sequence_design_shell_parity(strategy: DnaDesignSearchStrategy) {
    use crate::engine_shell::{execute_shell_command, parse_shell_line, quote_shell_arg};
    let mut request = request();
    request.search_strategy = strategy;
    let directory = tempfile::tempdir().unwrap();
    let request_path = directory.path().join("request with apostrophe's.json");
    let output_path = directory.path().join("preview.json");
    std::fs::write(&request_path, serde_json::to_vec_pretty(&request).unwrap()).unwrap();
    let command = parse_shell_line(&format!(
        "sequence-design plan {} --path {}",
        quote_shell_arg(&format!("@{}", request_path.display())),
        quote_shell_arg(&output_path.display().to_string())
    ))
    .unwrap();
    assert!(!command.is_state_mutating());
    let mut shell_engine = GentleEngine::new();
    let mut direct = GentleEngine::new();
    let result = execute_shell_command(&mut shell_engine, &command).unwrap();
    assert!(!result.state_changed);
    let saved: DnaSequenceDesignReport =
        serde_json::from_slice(&std::fs::read(&output_path).unwrap()).unwrap();
    assert_eq!(saved, plan(&mut direct, request));
    assert_eq!(
        result.output["result"]["dna_sequence_design"],
        serde_json::to_value(&saved).unwrap()
    );
    let command = parse_shell_line(&format!(
        "sequence-design apply {} --approve {}",
        quote_shell_arg(&format!("@{}", output_path.display())),
        saved.approval_digest.as_deref().unwrap()
    ))
    .unwrap();
    assert!(command.is_state_mutating());
    assert!(
        execute_shell_command(&mut shell_engine, &command)
            .unwrap()
            .state_changed
    );
    assert!(parse_shell_line("sequence-design apply @preview.json").is_err());
    let capability = execute_shell_command(
        &mut shell_engine,
        &parse_shell_line("introspect capabilities").unwrap(),
    )
    .unwrap();
    let text = capability.output.to_string();
    for id in [
        "sequence-design plan",
        "sequence-design apply",
        "PlanDnaSequenceDesign",
        "ApplyDnaSequenceDesign",
    ] {
        assert!(text.contains(id), "missing {id}");
    }
    assert!(crate::agent_bridge::AGENT_BRIDGE_SYSTEM_PROMPT.contains("sequence-design plan"));
    assert!(
        crate::agent_bridge::AGENT_BRIDGE_SYSTEM_PROMPT
            .contains("search_strategy=conflict_directed")
    );
    assert!(text.contains("conflict_directed"));
}

#[test]
fn sequence_design_admission_rejects_oversized_requests_and_previews_without_mutation() {
    let mut engine = GentleEngine::new();
    let original = serde_json::to_value(engine.state()).unwrap();
    let mut oversized = request();
    oversized
        .avoid_motifs
        .resize(core::MAX_MOTIFS + 1, oversized.avoid_motifs[0].clone());
    assert!(
        engine
            .apply(Operation::PlanDnaSequenceDesign {
                request: Box::new(oversized),
                path: None
            })
            .is_err()
    );
    let mut preview = plan(&mut engine, request());
    preview.output_sequence = Some("A".repeat(core::MAX_SEQUENCE_BP + 1));
    assert!(apply(&mut engine, preview).is_err());
    let preview = plan(&mut engine, request());
    for field in 0..6 {
        let mut oversized = preview.clone();
        let digest = "x".repeat(MAX_DIGEST_BYTES + 1);
        match field {
            0 => {
                oversized.request.target = DnaDesignTarget::InlineSequence {
                    sequence: "A".repeat(core::MAX_SEQUENCE_BP + 1),
                };
            }
            1 => {
                oversized.request.target = DnaDesignTarget::LoadedSequence {
                    seq_id: "x".repeat(MAX_SOURCE_ID_BYTES + 1),
                };
            }
            2 => oversized.source_sha256 = digest,
            3 => oversized.genetic_code_mapping_sha256 = digest,
            4 => oversized.source_features_sha256 = Some(digest),
            5 => oversized.output_sha256 = Some(digest),
            _ => unreachable!(),
        }
        // This is not just a stale-hash test: a client can re-hash its own JSON.
        oversized.approval_digest = Some(GentleEngine::dna_design_approval(&oversized).unwrap());
        let error = apply(&mut engine, oversized).unwrap_err();
        assert!(error.message.contains("bounded"), "{}", error.message);
    }
    assert_eq!(serde_json::to_value(engine.state()).unwrap(), original);
}
