//! Inline hand-crafted synthetic coding inserts, recreated from literal strings.
//! Used only for the shared engine/shell approval boundary, never as natural genes.
//! The two-edit MEF example recreates Glen's 2026-10-05 reporting-integrity audit.

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
        gc_content: None,
        max_evaluations: 4096,
        search_strategy: DnaDesignSearchStrategy::FullEnumeration,
        output_seq_id: "synthetic_without_ecori".into(),
    }
}

fn gc_request() -> DnaSequenceDesignRequest {
    let mut request = request();
    request.target = DnaDesignTarget::InlineSequence {
        sequence: "ATGTTTAAATAA".into(),
    };
    request.protein_sequence = "MFK".into();
    request.avoid_motifs.clear();
    request.gc_content = Some(DesignGcBounds {
        min_basis_points: 1666,
        max_basis_points: 1667,
    });
    request.output_seq_id = "synthetic_gc".into();
    request
}

#[test]
fn sequence_design_gc_facts_exact_apply_and_undo_share_the_validated_contract() {
    for strategy in [
        DnaDesignSearchStrategy::FullEnumeration,
        DnaDesignSearchStrategy::ConflictDirected,
    ] {
        let mut engine = GentleEngine::new();
        let original = serde_json::to_value(engine.state()).unwrap();
        let mut request = gc_request();
        request.search_strategy = strategy;
        let preview = plan(&mut engine, request);
        assert_eq!(
            preview.algorithm,
            search_strategy(strategy).algorithm_with_gc(true)
        );
        assert_eq!(preview.output_sequence.as_deref(), Some("ATGTTCAAATAA"));
        let facts = preview.gc_content.as_ref().unwrap();
        assert_eq!(
            (
                facts.denominator_bases,
                facts.minimum_gc_bases,
                facts.maximum_gc_bases
            ),
            (12, 2, 2)
        );
        assert_eq!(facts.input.gc_bases, 1);
        assert!(!facts.input.satisfies_bounds);
        assert_eq!(facts.output.as_ref().unwrap().gc_bases, 2);
        assert!(facts.output.as_ref().unwrap().satisfies_bounds);
        let decoded: DnaSequenceDesignReport =
            serde_json::from_str(&serde_json::to_string(&preview).unwrap()).unwrap();
        assert_eq!(decoded, preview);
        assert_eq!(serde_json::to_value(engine.state()).unwrap(), original);
        let journal_len = engine.journal.len();
        let undo_count = engine.undo_available();
        let applied = apply(&mut engine, decoded).unwrap();
        let receipt = applied.dna_sequence_design_receipt.unwrap();
        assert!(receipt.output_constraints_verified && !receipt.search_claims_verified);
        assert_single_design_derivation(&engine, "synthetic_gc", &applied.op_id);
        assert_eq!(engine.journal.len(), journal_len + 1);
        assert_eq!(engine.undo_available(), undo_count + 1);
        assert_eq!(
            engine.state.sequences["synthetic_gc"].get_forward_string(),
            "ATGTTCAAATAA"
        );
        engine.undo_last_operation().unwrap();
        assert_eq!(serde_json::to_value(engine.state()).unwrap(), original);
        assert_eq!(engine.journal.len(), journal_len);
        assert_eq!(engine.undo_available(), undo_count);
    }
}

#[test]
fn sequence_design_gc_rehashed_invalid_output_and_false_facts_are_rejected() {
    let mut engine = GentleEngine::new();
    let original = serde_json::to_value(engine.state()).unwrap();
    let preview = plan(&mut engine, gc_request());
    for kind in 0..5 {
        let mut changed = preview.clone();
        match kind {
            0 => {
                // Valid translation/frozen bases and no motifs, but GC is below minimum.
                changed.output_sequence = Some("ATGTTTAAATAA".into());
                changed.output_sha256 = Some(sha256_prefixed_str("ATGTTTAAATAA"));
                changed.edits.clear();
            }
            1 => {
                changed
                    .gc_content
                    .as_mut()
                    .unwrap()
                    .output
                    .as_mut()
                    .unwrap()
                    .gc_bases = 3
            }
            2 => changed.gc_content.as_mut().unwrap().input.gc_bases = 2,
            3 => changed.request.gc_content = None,
            4 => changed.algorithm = core::ALGORITHM.into(),
            _ => unreachable!(),
        }
        changed.approval_digest = Some(GentleEngine::dna_design_approval(&changed).unwrap());
        assert!(apply(&mut engine, changed).is_err(), "tamper {kind}");
        assert_eq!(serde_json::to_value(engine.state()).unwrap(), original);
    }
}

#[test]
fn sequence_design_gc_unknown_output_and_unverified_search_claims_remain_explicit() {
    let mut engine = GentleEngine::new();
    let mut request = gc_request();
    request.max_evaluations = 1;
    let stopped = plan(&mut engine, request.clone());
    assert_eq!(stopped.status, DnaDesignStatus::SearchExhausted);
    assert!(stopped.approval_digest.is_none());
    assert!(stopped.gc_content.unwrap().output.is_none());
    request.max_evaluations = 4096;
    let cancelled = engine
        .plan_dna_sequence_design(request.clone(), |n| n == 2)
        .unwrap();
    assert_eq!(cancelled.status, DnaDesignStatus::Cancelled);
    assert!(cancelled.output_sequence.is_none() && cancelled.approval_digest.is_none());
    assert!(!cancelled.minimum_edits_proven && !cancelled.optimization_complete);
    assert!(cancelled.gc_content.unwrap().output.is_none());
    for (min, max, expected) in [
        (5000, 4999, DnaDesignStatus::Invalid),
        (1667, 1667, DnaDesignStatus::ProvenInfeasible),
    ] {
        request.gc_content = Some(DesignGcBounds {
            min_basis_points: min,
            max_basis_points: max,
        });
        let preview = plan(&mut engine, request.clone());
        assert_eq!(preview.status, expected);
        assert!(preview.approval_digest.is_none());
        if expected == DnaDesignStatus::Invalid {
            assert!(preview.gc_content.is_none());
        } else {
            assert!(preview.gc_content.unwrap().output.is_none());
        }
    }
    // A client can submit a different valid nonminimum output. Constraint validity
    // is verified, but its rehashed completeness/minimum claim is not authenticated.
    request.gc_content = Some(DesignGcBounds {
        min_basis_points: 1666,
        max_basis_points: 2500,
    });
    let mut submitted = plan(&mut engine, request);
    submitted.output_sequence = Some("ATGTTCAAGTAA".into());
    submitted.output_sha256 = Some(sha256_prefixed_str("ATGTTCAAGTAA"));
    submitted.edits = vec![
        DesignNucleotideEdit {
            position_0based: 5,
            before: "T".into(),
            after: "C".into(),
        },
        DesignNucleotideEdit {
            position_0based: 8,
            before: "A".into(),
            after: "G".into(),
        },
    ];
    submitted.gc_content.as_mut().unwrap().output = Some(DesignGcMeasurement {
        gc_bases: 3,
        satisfies_bounds: true,
    });
    submitted.approval_digest = Some(GentleEngine::dna_design_approval(&submitted).unwrap());
    let receipt = apply(&mut engine, submitted)
        .unwrap()
        .dna_sequence_design_receipt
        .unwrap();
    assert!(receipt.output_constraints_verified && !receipt.search_claims_verified);
}

#[test]
fn sequence_design_conflict_search_is_explicit_bound_and_separately_approved() {
    let mut engine = GentleEngine::new();
    let original = serde_json::to_value(engine.state()).unwrap();
    let legacy = plan(&mut engine, request());
    assert_eq!(legacy.algorithm, core::ALGORITHM);
    let legacy_json = serde_json::to_value(&legacy).unwrap();
    assert!(legacy_json["request"].get("search_strategy").is_none());
    assert!(legacy_json["request"].get("gc_content").is_none());
    assert!(legacy_json.get("gc_content").is_none());
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

fn assert_single_design_derivation(engine: &GentleEngine, seq_id: &str, op_id: &str) {
    let containers = engine
        .state
        .container_state
        .containers
        .values()
        .filter(|container| container.members.iter().any(|member| member == seq_id))
        .collect::<Vec<_>>();
    assert_eq!(containers.len(), 1, "one apply must create one container");
    let container = containers[0];
    assert!(matches!(&container.kind, ContainerKind::Singleton));
    assert_eq!(container.members, vec![seq_id.to_string()]);
    assert_eq!(container.created_by_op.as_deref(), Some(op_id));
    assert_eq!(
        engine
            .state
            .container_state
            .seq_to_latest_container
            .get(seq_id),
        Some(&container.container_id)
    );
    let nodes = engine
        .state
        .lineage
        .nodes
        .values()
        .filter(|node| node.seq_id == seq_id)
        .collect::<Vec<_>>();
    assert_eq!(nodes.len(), 1, "one apply must create one lineage node");
    let node = nodes[0];
    assert!(matches!(&node.origin, SequenceOrigin::Derived));
    assert_eq!(node.created_by_op.as_deref(), Some(op_id));
    assert_eq!(
        engine.state.lineage.seq_to_node.get(seq_id),
        Some(&node.node_id)
    );
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
    let journal_len = engine.journal.len();
    let undo_count = engine.undo_available();
    let applied = apply(&mut engine, roundtrip).unwrap();
    let receipt = applied.dna_sequence_design_receipt.unwrap();
    assert_eq!(applied.created_seq_ids, [receipt.created_seq_id.clone()]);
    assert_eq!(engine.journal.len(), journal_len + 1);
    assert_eq!(engine.undo_available(), undo_count + 1);
    assert!(receipt.output_constraints_verified);
    assert!(!receipt.search_claims_verified);
    assert_single_design_derivation(&engine, &receipt.created_seq_id, &applied.op_id);
    assert_eq!(engine.state.container_state.containers.len(), 1);
    assert_eq!(engine.state.lineage.nodes.len(), 1);
    assert!(engine.state.lineage.edges.is_empty());
    let dna = &engine.state.sequences[&receipt.created_seq_id];
    assert_eq!(
        dna.get_forward_string(),
        preview.output_sequence.as_deref().unwrap()
    );
    assert_eq!(dna.features().len(), 1);
    assert_eq!(&*dna.features()[0].kind, "CDS");
    let record = &engine.state.metadata[&format!("dna_sequence_design:{}", receipt.created_seq_id)];
    assert_eq!(
        record["submitted_proposal"],
        serde_json::to_value(&preview).unwrap()
    );
    assert!(record.get("proposal").is_none());
    let applied_state = serde_json::to_value(engine.state()).unwrap();
    assert!(apply(&mut engine, preview.clone()).is_err());
    assert_eq!(serde_json::to_value(engine.state()).unwrap(), applied_state);
    engine.undo_last_operation().unwrap();
    assert!(!engine.state.sequences.contains_key(&receipt.created_seq_id));
    assert!(
        !engine
            .state
            .metadata
            .contains_key(&format!("dna_sequence_design:{}", receipt.created_seq_id))
    );
    assert_eq!(serde_json::to_value(engine.state()).unwrap(), original);
    assert_eq!(engine.journal.len(), journal_len);
    assert_eq!(engine.undo_available(), undo_count);
    engine.redo_last_operation().unwrap();
    assert_eq!(serde_json::to_value(engine.state()).unwrap(), applied_state);
    assert_single_design_derivation(&engine, &receipt.created_seq_id, &applied.op_id);
    assert_eq!(engine.journal.len(), journal_len + 1);
    assert_eq!(engine.undo_available(), undo_count + 1);
}

#[test]
fn sequence_design_direct_apply_method_uses_one_operation_and_undo_boundary() {
    let mut engine = GentleEngine::new();
    let preview = plan(&mut engine, request());
    let original = serde_json::to_value(engine.state()).unwrap();
    let journal_len = engine.journal.len();
    let approval = preview.approval_digest.clone().unwrap();
    let receipt = engine
        .apply_dna_sequence_design(preview, &approval)
        .unwrap();
    assert_eq!(engine.journal.len(), journal_len + 1);
    let operation = engine.journal.last().unwrap();
    assert!(matches!(
        &operation.op,
        Operation::ApplyDnaSequenceDesign { .. }
    ));
    assert_single_design_derivation(&engine, &receipt.created_seq_id, &operation.result.op_id);
    assert_eq!(engine.undo_available(), 1);
    engine.undo_last_operation().unwrap();
    assert_eq!(serde_json::to_value(engine.state()).unwrap(), original);
}

#[test]
fn sequence_design_rehashed_nonminimum_output_does_not_verify_search_claims() {
    for strategy in [
        DnaDesignSearchStrategy::FullEnumeration,
        DnaDesignSearchStrategy::ConflictDirected,
    ] {
        let mut engine = GentleEngine::new();
        let mut request = request();
        request.search_strategy = strategy;
        let mut submitted = plan(&mut engine, request);
        assert_eq!(submitted.output_sequence.as_deref(), Some("ATGGAATTTTAA"));
        assert_eq!(submitted.edits.len(), 1);
        assert!(submitted.minimum_edits_proven && submitted.optimization_complete);
        let original = serde_json::to_value(engine.state()).unwrap();

        // GAG and TTT still encode EF, but need two edits instead of the one above.
        submitted.output_sequence = Some("ATGGAGTTTTAA".into());
        submitted.output_sha256 = Some(sha256_prefixed_str("ATGGAGTTTTAA"));
        submitted.edits = vec![
            DesignNucleotideEdit {
                position_0based: 5,
                before: "A".into(),
                after: "G".into(),
            },
            DesignNucleotideEdit {
                position_0based: 8,
                before: "C".into(),
                after: "T".into(),
            },
        ];
        submitted.reason = "Fabricated completed search proving a two-edit minimum".into();
        submitted.evaluated_candidates = 0;
        submitted.approval_digest = Some(GentleEngine::dna_design_approval(&submitted).unwrap());

        let mut invalid = submitted.clone();
        invalid.output_sequence = Some("ATGGAGTTTTAG".into());
        invalid.output_sha256 = Some(sha256_prefixed_str("ATGGAGTTTTAG"));
        invalid.edits.push(DesignNucleotideEdit {
            position_0based: 11,
            before: "A".into(),
            after: "G".into(),
        });
        invalid.approval_digest = Some(GentleEngine::dna_design_approval(&invalid).unwrap());
        assert!(
            apply(&mut engine, invalid).is_err(),
            "frozen stop must remain enforced"
        );
        assert_eq!(serde_json::to_value(engine.state()).unwrap(), original);

        let applied = apply(&mut engine, submitted.clone()).unwrap();
        let receipt = applied.dna_sequence_design_receipt.unwrap();
        assert!(receipt.output_constraints_verified);
        assert!(!receipt.search_claims_verified);
        assert_eq!(
            receipt.approval_digest,
            submitted.approval_digest.clone().unwrap()
        );
        assert_eq!(
            engine.state.sequences[&receipt.created_seq_id].get_forward_string(),
            "ATGGAGTTTTAA"
        );
        assert!(
            receipt
                .nonclaims
                .iter()
                .any(|statement| statement == APPLY_SEARCH_NONCLAIM)
        );
        let record =
            &engine.state.metadata[&format!("dna_sequence_design:{}", receipt.created_seq_id)];
        assert_eq!(
            record["submitted_proposal"],
            serde_json::to_value(&submitted).unwrap()
        );
        assert_eq!(record["receipt"], serde_json::to_value(&receipt).unwrap());
        assert!(record.get("proposal").is_none());
        assert_single_design_derivation(&engine, &receipt.created_seq_id, &applied.op_id);
        engine.undo_last_operation().unwrap();
        assert_eq!(serde_json::to_value(engine.state()).unwrap(), original);
    }
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
    let applied = apply(&mut engine, preview.clone()).unwrap();
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
    assert_single_design_derivation(&engine, "synthetic_without_ecori", &applied.op_id);
    assert_eq!(engine.state.container_state.containers.len(), 2);
    assert_eq!(engine.state.lineage.nodes.len(), 2);
    assert_eq!(engine.state.lineage.edges.len(), 1);
    let edge = &engine.state.lineage.edges[0];
    assert_eq!(
        edge.from_node_id,
        engine.state.lineage.seq_to_node["synthetic_source"]
    );
    assert_eq!(
        edge.to_node_id,
        engine.state.lineage.seq_to_node["synthetic_without_ecori"]
    );
    assert_eq!(edge.op_id, applied.op_id);
    assert_eq!(edge.run_id, "interactive");
    let source_node = &edge.from_node_id;
    assert_eq!(
        serde_json::to_value(&engine.state.lineage.nodes[source_node]).unwrap(),
        original["lineage"]["nodes"][source_node]
    );
    let source_container =
        &engine.state.container_state.seq_to_latest_container["synthetic_source"];
    assert_eq!(
        serde_json::to_value(&engine.state.container_state.containers[source_container]).unwrap(),
        original["container_state"]["containers"][source_container]
    );
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
    assert_sequence_design_shell_parity(DnaDesignSearchStrategy::FullEnumeration, false);
}

#[test]
fn sequence_design_conflict_shared_shell_discovery_and_json_file_parity() {
    assert_sequence_design_shell_parity(DnaDesignSearchStrategy::ConflictDirected, false);
}

#[test]
fn sequence_design_gc_shared_shell_discovery_and_json_file_parity() {
    for strategy in [
        DnaDesignSearchStrategy::FullEnumeration,
        DnaDesignSearchStrategy::ConflictDirected,
    ] {
        assert_sequence_design_shell_parity(strategy, true);
    }
}

fn assert_sequence_design_shell_parity(strategy: DnaDesignSearchStrategy, gc: bool) {
    use crate::engine_shell::{execute_shell_command, parse_shell_line, quote_shell_arg};
    let mut request = request();
    if gc {
        request = gc_request();
    }
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
    let applied = execute_shell_command(&mut shell_engine, &command).unwrap();
    assert!(applied.state_changed);
    let direct_applied = apply(&mut direct, saved).unwrap();
    assert_eq!(
        applied.output["result"]["dna_sequence_design_receipt"],
        serde_json::to_value(direct_applied.dna_sequence_design_receipt.unwrap()).unwrap()
    );
    assert_single_design_derivation(
        &shell_engine,
        if gc {
            "synthetic_gc"
        } else {
            "synthetic_without_ecori"
        },
        applied.output["result"]["op_id"].as_str().unwrap(),
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
    assert!(text.contains("gc_content"));
    assert!(crate::agent_bridge::AGENT_BRIDGE_SYSTEM_PROMPT.contains("min_basis_points"));
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
