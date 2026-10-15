//! Explicit, provenance-preserving gene-set to serial/storage-rack projection.
//! Reuses existing source bindings and arrangement/rack creation, not pooling.

use super::*;

impl GentleEngine {
    pub(super) fn create_gene_set_arrangement(
        &mut self,
        resolution_id: &str,
        member_bindings: Vec<GeneSetPoolMemberBinding>,
        arrangement_id: &str,
        name: Option<String>,
        confirm_physical_placement: bool,
        op_id: &str,
        run_id: &str,
    ) -> Result<(CollectionOperationReport, String), EngineError> {
        if !confirm_physical_placement {
            return Err(EngineError::invalid_input(
                "Gene-set storage projection requires confirm_physical_placement=true; a logical set never authorizes physical placement",
            ));
        }
        let arrangement_id = arrangement_id.trim();
        if arrangement_id.is_empty()
            || self
                .state
                .container_state
                .arrangements
                .contains_key(arrangement_id)
            || self
                .state
                .metadata
                .contains_key(&format!("gene_set_arrangement:{arrangement_id}"))
        {
            return Err(EngineError::invalid_input(
                "Gene-set storage requires a fresh non-empty arrangement_id",
            ));
        }
        let subject = CollectionSubjectRef::GeneSetResolution {
            report_id: resolution_id.trim().to_string(),
        };
        let resolution = self.get_gene_set_resolution_artifact(resolution_id.trim())?;
        if resolution.unresolved_member_count > 0 || !resolution.unresolved_members.is_empty() {
            return Err(EngineError::invalid_input(
                "Storage projection requires a complete resolved gene set",
            ));
        }
        let resolution_sha256 = sha256_prefixed_str(
            &serde_json::to_string(&resolution)
                .map_err(|e| EngineError::internal(e.to_string()))?,
        );
        let (subject, members, mut provenance, biological_contexts) = self
            .collection_subject_members_for_project_sequences_or_gene_set(
                &subject,
                "storage projection",
            )?;
        if members.is_empty() || members.len() > 24 {
            return Err(EngineError::invalid_input(
                "Minimal storage projection supports 1..24 members on small_tube_4x6",
            ));
        }
        let bindings = Self::gene_set_pool_binding_map(member_bindings)?;
        let member_ids = members
            .iter()
            .map(|m| m.stable_member_id.as_str())
            .collect::<BTreeSet<_>>();
        if member_ids.len() != members.len()
            || bindings.len() != members.len()
            || bindings.keys().any(|id| !member_ids.contains(id.as_str()))
        {
            return Err(EngineError::invalid_input(
                "Storage projection requires exactly one binding for every distinct resolved member, with no extra bindings",
            ));
        }
        let lift_policy = collection_lift_policy(
            CapabilitySource::EngineOperation,
            "CreateGeneSetArrangement",
            subject.kind(),
        )
        .cloned()
        .ok_or_else(|| EngineError::internal("Missing storage projection lift policy"))?;
        if lift_policy.context_requirement != CollectionContextRequirement::ContextAgnostic
            || !matches!(
                &lift_policy.support,
                CollectionLiftSupport::Supported {
                    mode: CollectionLiftingMode::Arrange,
                    ..
                }
            )
        {
            return Err(EngineError::internal(
                "Invalid storage projection lift policy",
            ));
        }
        let fingerprint = sha256_prefixed_str(&canonical_collection_membership_json(
            subject.kind(),
            &members,
        ));
        let mut containers = Vec::with_capacity(members.len());
        let mut rows = Vec::with_capacity(members.len());
        let mut source_bindings = Vec::with_capacity(members.len());
        for mut member in members {
            let container_id = &bindings[&member.stable_member_id];
            let container = self
                .state
                .container_state
                .containers
                .get(container_id)
                .ok_or_else(|| {
                    EngineError::invalid_input(format!("Source tube '{container_id}' not found"))
                })?;
            if !matches!(&container.kind, ContainerKind::Singleton)
                || !container.declared_contents_exclusive
                || container.members.len() != 1
            {
                return Err(EngineError::invalid_input(format!(
                    "Source tube '{container_id}' must be an exclusive singleton"
                )));
            }
            let seq_id = &container.members[0];
            let sequence = self.state.sequences.get(seq_id).ok_or_else(|| {
                EngineError::invalid_input(format!(
                    "Source tube '{container_id}' has a missing sequence"
                ))
            })?;
            let snapshot =
                Self::collection_sequence_snapshot_sha256(sequence, "storage projection source")?;
            source_bindings.push(json!({"stable_member_id":member.stable_member_id,"source_container_id":container_id,"source_seq_id":seq_id,"sequence_snapshot_sha256":snapshot}));
            member.seq_id = Some(seq_id.clone());
            member.source_provenance.push(GeneSetProvenanceRow {
                source_kind: "physical_source_container".into(), source_id: container_id.clone(),
                source_label: container.name.clone(), source_path: None,
                note: Some(format!("Explicit caller binding; sequence {seq_id}; snapshot {snapshot}; not proof of physical stock")),
            });
            containers.push(container_id.clone());
            rows.push(CollectionMemberStatusRow {
                member,
                outcome: CollectionMemberOutcome::Succeeded,
                error: None,
                produced_report_ids: vec![],
            });
        }
        let identity = serde_json::to_string(&json!({"resolution_sha256":resolution_sha256,"arrangement_id":arrangement_id,"members":rows,"name":name})).map_err(|e| EngineError::internal(e.to_string()))?;
        provenance.push(GeneSetProvenanceRow {
            source_kind: "serial_arrangement".into(), source_id: arrangement_id.into(),
            source_label: name.clone(), source_path: None,
            note: Some("Explicit storage projection; logical set and source tubes retained without aliquoting or pooling".into()),
        });
        let report = CollectionOperationReport {
            schema: COLLECTION_OPERATION_REPORT_SCHEMA.into(), report_id: short_sha256_id("gene_set_arrangement", &identity),
            op_id: Some(op_id.into()), run_id: Some(run_id.into()), generated_at_unix_ms: Self::now_unix_ms(),
            capability_source: CapabilitySource::EngineOperation, capability_name: "CreateGeneSetArrangement".into(),
            collection_subject: subject, lifting_mode: CollectionLiftingMode::Arrange, lift_policy,
            fingerprint_algorithm: COLLECTION_MEMBERSHIP_FINGERPRINT_ALGORITHM.into(), collection_membership_fingerprint_sha256: fingerprint,
            biological_contexts, dry_run: false, applied: true, per_member_status: rows,
            aggregate_warnings: vec!["Digital placement draft, not physical stock, freezer tracking or fit/manufacturing validation".into()], provenance,
        };
        // Serialize provenance before creating the layout; retain the set snapshot
        // and ordered tube mapping, not a mutable live-only association.
        let mut record = json!({"schema":"gentle.gene_set_storage_projection.v1","source_resolution_sha256":resolution_sha256,"source_bindings":source_bindings,"report":report,"arrangement_id":arrangement_id,"physical_template":"storage_pcr_tube_rack"});
        self.add_serial_arrangement(
            &containers,
            Some(arrangement_id.into()),
            name,
            None,
            None,
            Some(op_id),
        )?;
        let rack_id = self.state.container_state.arrangements[arrangement_id]
            .default_rack_id
            .clone()
            .expect("serial arrangement creates its rack");
        record["rack_id"] = json!(rack_id);
        self.state
            .metadata
            .insert(format!("gene_set_arrangement:{arrangement_id}"), record);
        Ok((report, rack_id))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn fixture() -> (GentleEngine, Operation) {
        // Literal synthetic tubes: no biological accession or physical stock.
        let mut state = ProjectState::default();
        for id in ["tube_a", "tube_b"] {
            state
                .sequences
                .insert(id.into(), DNAsequence::from_sequence("ACGTACGT").unwrap());
        }
        let mut engine = GentleEngine::from_state(state);
        engine
            .upsert_gene_set_resolution_artifact(GeneSetResolutionReport {
                schema: GENE_SET_RESOLUTION_SCHEMA.into(),
                op_id: Some("synthetic_storage_set".into()),
                resolved_member_count: 2,
                resolved_members: ["GENE_B", "GENE_A"]
                    .into_iter()
                    .map(|id| GeneSetResolvedMember {
                        dedup_key: id.into(),
                        symbol: id.into(),
                        ..Default::default()
                    })
                    .collect(),
                ..Default::default()
            })
            .unwrap();
        let member_bindings = [("GENE_A", "tube_a"), ("GENE_B", "tube_b")]
            .into_iter()
            .map(|(member, seq)| GeneSetPoolMemberBinding {
                stable_member_id: member.into(),
                source_container_id: engine.state().container_state.seq_to_latest_container[seq]
                    .clone(),
            })
            .collect();
        (
            engine,
            Operation::CreateGeneSetArrangement {
                resolution_id: "synthetic_storage_set".into(),
                member_bindings,
                arrangement_id: "storage_review".into(),
                name: Some("Synthetic storage".into()),
                confirm_physical_placement: true,
            },
        )
    }

    #[test]
    fn gene_set_arrangement_preserves_order_sources_provenance_and_undo() {
        let (mut engine, operation) = fixture();
        let sources = engine
            .state()
            .container_state
            .seq_to_latest_container
            .clone();
        let result = engine.apply(operation).unwrap();
        let report = result.collection_operation.unwrap();
        assert_eq!(report.lifting_mode, CollectionLiftingMode::Arrange);
        assert_eq!(
            report
                .per_member_status
                .iter()
                .map(|row| row.member.stable_member_id.as_str())
                .collect::<Vec<_>>(),
            vec!["GENE_B", "GENE_A"]
        );
        assert!(report.per_member_status.iter().all(|row| {
            row.member.source_provenance.iter().any(|p| {
                p.source_kind == "physical_source_container"
                    && p.note.as_ref().unwrap().contains("snapshot sha256:")
            })
        }));
        let arrangement = &engine.state().container_state.arrangements["storage_review"];
        assert_eq!(
            arrangement.lane_container_ids,
            vec![sources["tube_b"].clone(), sources["tube_a"].clone()]
        );
        assert_eq!(
            arrangement.created_by_op.as_deref(),
            Some(result.op_id.as_str())
        );
        let rack =
            &engine.state().container_state.racks[arrangement.default_rack_id.as_ref().unwrap()];
        assert_eq!(rack.profile.kind, RackProfileKind::SmallTube4x6);
        assert_eq!(rack.placements.len(), 2);
        assert_eq!(engine.state().sequences.len(), 2);
        assert_eq!(engine.state().container_state.containers.len(), 2);
        assert_eq!(
            engine.state().container_state.seq_to_latest_container,
            sources
        );
        let record = engine.state().metadata["gene_set_arrangement:storage_review"].clone();
        assert_eq!(record["schema"], "gentle.gene_set_storage_projection.v1");
        assert_eq!(
            record["report"]["collection_subject"]["report_id"],
            "synthetic_storage_set"
        );
        assert!(
            record["source_resolution_sha256"]
                .as_str()
                .unwrap()
                .starts_with("sha256:")
        );
        let reopened: ProjectState =
            serde_json::from_str(&serde_json::to_string(engine.state()).unwrap()).unwrap();
        assert_eq!(
            reopened.metadata["gene_set_arrangement:storage_review"],
            record
        );
        engine.undo_last_operation().unwrap();
        assert!(
            !engine
                .state()
                .container_state
                .arrangements
                .contains_key("storage_review")
        );
        assert!(
            !engine
                .state()
                .metadata
                .contains_key("gene_set_arrangement:storage_review")
        );
        engine.redo_last_operation().unwrap();
        assert_eq!(
            engine.state().metadata["gene_set_arrangement:storage_review"],
            record
        );
    }

    #[test]
    fn gene_set_arrangement_refuses_implicit_missing_duplicate_and_nonexclusive_bindings() {
        for case in [
            "implicit",
            "missing",
            "extra",
            "duplicate",
            "nonexclusive",
            "unresolved",
            "missing_sequence",
        ] {
            let (mut engine, mut operation) = fixture();
            let Operation::CreateGeneSetArrangement {
                member_bindings,
                confirm_physical_placement,
                ..
            } = &mut operation
            else {
                unreachable!()
            };
            match case {
                "implicit" => *confirm_physical_placement = false,
                "missing" => {
                    member_bindings.pop();
                }
                "extra" => member_bindings.push(GeneSetPoolMemberBinding {
                    stable_member_id: "NOT_IN_SET".into(),
                    source_container_id: "unknown".into(),
                }),
                "duplicate" => member_bindings.push(member_bindings[0].clone()),
                "nonexclusive" => {
                    engine
                        .state
                        .container_state
                        .containers
                        .get_mut(&member_bindings[0].source_container_id)
                        .unwrap()
                        .declared_contents_exclusive = false;
                }
                "unresolved" => {
                    let mut report = engine
                        .get_gene_set_resolution_artifact("synthetic_storage_set")
                        .unwrap();
                    report.unresolved_member_count = 1;
                    engine.upsert_gene_set_resolution_artifact(report).unwrap();
                }
                "missing_sequence" => {
                    engine.state.sequences.remove("tube_a");
                }
                _ => unreachable!(),
            }
            let before = serde_json::to_string(engine.state()).unwrap();
            assert!(engine.apply(operation).is_err(), "{case}");
            assert_eq!(
                serde_json::to_string(engine.state()).unwrap(),
                before,
                "{case}"
            );
        }
    }

    #[test]
    fn gene_set_arrangement_generic_shell_and_engine_agree() {
        let (mut direct, operation) = fixture();
        let (mut shell, _) = fixture();
        let expected = direct
            .apply(operation.clone())
            .unwrap()
            .collection_operation
            .unwrap();
        let command = crate::engine_shell::parse_shell_line(&format!(
            "op {}",
            crate::engine_shell::quote_shell_arg(&serde_json::to_string(&operation).unwrap())
        ))
        .unwrap();
        let actual = crate::engine_shell::execute_shell_command(&mut shell, &command).unwrap();
        assert!(actual.state_changed);
        let mut actual = actual.output["result"]["collection_operation"].clone();
        let mut expected = serde_json::to_value(expected).unwrap();
        for value in [&mut actual, &mut expected] {
            value
                .as_object_mut()
                .unwrap()
                .remove("generated_at_unix_ms");
            value.as_object_mut().unwrap().remove("run_id");
        }
        assert_eq!(actual, expected);
    }

    #[test]
    fn gene_set_arrangement_identity_is_stable_and_capacity_collision_refusals_are_atomic() {
        let (mut first, operation) = fixture();
        let (mut second, mut reversed) = fixture();
        if let Operation::CreateGeneSetArrangement {
            member_bindings, ..
        } = &mut reversed
        {
            member_bindings.reverse();
        }
        let left = first
            .apply(operation.clone())
            .unwrap()
            .collection_operation
            .unwrap();
        let right = second
            .apply(reversed)
            .unwrap()
            .collection_operation
            .unwrap();
        assert_eq!(left.report_id, right.report_id);
        assert_eq!(
            first.state().metadata["gene_set_arrangement:storage_review"]["source_bindings"],
            second.state().metadata["gene_set_arrangement:storage_review"]["source_bindings"]
        );
        let before = serde_json::to_string(first.state()).unwrap();
        assert!(first.apply(operation).is_err());
        assert_eq!(serde_json::to_string(first.state()).unwrap(), before);
        let (mut engine, operation) = fixture();
        let mut resolution = engine
            .get_gene_set_resolution_artifact("synthetic_storage_set")
            .unwrap();
        resolution.resolved_members = (0..25)
            .map(|index| GeneSetResolvedMember {
                dedup_key: format!("GENE_{index}"),
                symbol: format!("GENE_{index}"),
                ..Default::default()
            })
            .collect();
        resolution.resolved_member_count = 25;
        engine
            .upsert_gene_set_resolution_artifact(resolution)
            .unwrap();
        let before = serde_json::to_string(engine.state()).unwrap();
        assert!(
            engine
                .apply(operation)
                .unwrap_err()
                .message
                .contains("1..24")
        );
        assert_eq!(serde_json::to_string(engine.state()).unwrap(), before);
    }
}
