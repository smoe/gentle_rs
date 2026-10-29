//! Slice 6 synthetic persistence regressions. All DNA, indexes and JSONL rows
//! are invented below, recreated in temporary folders through the real saver.
//! No user projects, external tools or scientific data are used.

use super::*;

struct Fixture {
    project: PathBuf,
    index: PathBuf,
    records: PathBuf,
    reference: serde_json::Value,
}

fn fixture(root: &Path, ending: &str) -> Fixture {
    fs::create_dir_all(root).unwrap();
    let project = root.join("synthetic's project.json");
    let mut state = ProjectState::default();
    state.sequences.insert("source".into(), seq("ACGTACGT"));
    state.metadata.insert(
        CANDIDATE_SETS_METADATA_KEY.into(),
        json!({
            "schema": CANDIDATE_SETS_SCHEMA, "updated_at_unix_ms": 0,
            "sets": {"windows": {"name":"windows", "created_at_unix_ms":0,
                "source_seq_ids":["source"], "candidates":[{"seq_id":"source",
                    "start_0based":0, "end_0based":4, "sequence":"ACGT", "metrics":{}}]}}
        }),
    );
    state.save_to_path(project.to_str().unwrap()).unwrap();
    let mut disk: serde_json::Value = serde_json::from_slice(&fs::read(&project).unwrap()).unwrap();
    disk["metadata"][CANDIDATE_SETS_METADATA_KEY]["synthetic_extension"] =
        json!({"preserve": "verbatim"});
    fs::write(&project, serde_json::to_string_pretty(&disk).unwrap()).unwrap();
    let reference = disk["metadata"][CANDIDATE_SETS_METADATA_KEY].clone();
    let index = root.join(reference["index_path"].as_str().unwrap());
    let index_value: serde_json::Value =
        serde_json::from_slice(&fs::read(&index).unwrap()).unwrap();
    let records = index
        .parent()
        .unwrap()
        .join(index_value["sets"][0]["records_path"].as_str().unwrap());
    for path in [&project, &index, &records] {
        let text = fs::read_to_string(path).unwrap();
        fs::write(path, text.replace('\n', ending)).unwrap();
    }
    Fixture {
        project,
        index,
        records,
        reference,
    }
}

fn tree(root: &Path) -> BTreeMap<PathBuf, Option<Vec<u8>>> {
    fn walk(root: &Path, dir: &Path, out: &mut BTreeMap<PathBuf, Option<Vec<u8>>>) {
        for entry in fs::read_dir(dir).unwrap() {
            let path = entry.unwrap().path();
            let relative = path.strip_prefix(root).unwrap().to_path_buf();
            if path.is_dir() {
                out.insert(relative, None);
                walk(root, &path, out);
            } else {
                out.insert(relative, Some(fs::read(path).unwrap()));
            }
        }
    }
    let mut out = BTreeMap::new();
    walk(root, root, &mut out);
    out
}

fn unavailable(error: EngineError) {
    assert_eq!(error.code, ErrorCode::InvalidInput);
    assert!(
        error.message.starts_with("candidate_store_unavailable:"),
        "{error}"
    );
    assert!(!error.cause_chain.is_empty());
}

#[test]
fn candidate_store_unresolved_load_preserves_reference_and_refuses_save_and_save_as() {
    let _lock = candidate_store_env_lock().lock().unwrap();
    let _env = EnvVarGuard::set(CANDIDATE_STORE_STRICT_LOAD_ENV, "0");
    for ending in ["\n", "\r\n"] {
        for fault in [
            "missing_index",
            "corrupt_index",
            "unsupported_index",
            "missing_records",
            "corrupt_records",
            "malformed_reference",
            "backslash_path",
        ] {
            let td = tempdir().unwrap();
            let mut f = fixture(td.path(), ending);
            fs::write(
                f.index.parent().unwrap().join("intact-sibling.jsonl"),
                b"synthetic sibling\n",
            )
            .unwrap();
            match fault {
                "missing_index" => fs::remove_file(&f.index).unwrap(),
                "corrupt_index" => fs::write(&f.index, b"not json").unwrap(),
                "unsupported_index" => {
                    let mut value: serde_json::Value =
                        serde_json::from_slice(&fs::read(&f.index).unwrap()).unwrap();
                    value["schema"] = json!("synthetic.unsupported_index.v9");
                    fs::write(&f.index, serde_json::to_vec(&value).unwrap()).unwrap();
                }
                "malformed_reference" | "backslash_path" => {
                    f.reference["index_path"] = if fault == "malformed_reference" {
                        json!(42)
                    } else {
                        json!("absent-synthetic.candidates\\index.json")
                    };
                    let mut value: serde_json::Value =
                        serde_json::from_slice(&fs::read(&f.project).unwrap()).unwrap();
                    value["metadata"][CANDIDATE_SETS_METADATA_KEY] = f.reference.clone();
                    fs::write(&f.project, serde_json::to_vec(&value).unwrap()).unwrap();
                }
                "missing_records" => fs::remove_file(&f.records).unwrap(),
                "corrupt_records" => fs::write(&f.records, b"not json\n").unwrap(),
                _ => unreachable!(),
            }
            let existing = td.path().join("existing's target.json");
            fs::write(&existing, b"preserve destination").unwrap();
            let destination = td.path().join("destination");
            fs::create_dir(&destination).unwrap();
            let decoy = destination.join(f.index.parent().unwrap().file_name().unwrap());
            fs::create_dir(&decoy).unwrap();
            fs::write(
                decoy.join("index.json"),
                b"do not hydrate or replace this foreign index",
            )
            .unwrap();
            let before = tree(td.path());
            let state = ProjectState::load_from_path(f.project.to_str().unwrap()).unwrap();
            assert_eq!(state.metadata[CANDIDATE_SETS_METADATA_KEY], f.reference);
            let warning = &state.metadata[CANDIDATE_SETS_LOAD_WARNING_METADATA_KEY];
            assert_eq!(warning["schema"], CANDIDATE_STORE_LOAD_DIAGNOSTIC_SCHEMA);
            assert_eq!(warning["status"], "unresolved");
            assert_eq!(warning["code"], "candidate_store_unavailable");
            assert!(warning["message"].is_string());
            for target in [
                &f.project,
                &existing,
                &destination.join("save-as.json"),
                &td.path().join("absent parent/new project.json"),
            ] {
                unavailable(state.save_to_path(target.to_str().unwrap()).unwrap_err());
                assert_eq!(
                    tree(td.path()),
                    before,
                    "{fault}: no write, staging or directory creation"
                );
            }
        }
    }
}

#[test]
fn candidate_store_unresolved_reads_and_writes_refuse_without_hiding_other_sequences() {
    let _lock = candidate_store_env_lock().lock().unwrap();
    let _env = EnvVarGuard::set(CANDIDATE_STORE_STRICT_LOAD_ENV, "0");
    let td = tempdir().unwrap();
    let f = fixture(td.path(), "\n");
    fs::remove_file(&f.index).unwrap();
    let mut engine = GentleEngine::from_state(
        ProjectState::load_from_path(f.project.to_str().unwrap()).unwrap(),
    );
    let before = engine.state().metadata.clone();
    unavailable(engine.list_candidate_sets().unwrap_err());
    unavailable(
        engine
            .inspect_candidate_set_page("windows", 10, 0)
            .unwrap_err(),
    );
    unavailable(engine.list_candidate_set_metrics("windows").unwrap_err());
    unavailable(
        engine
            .apply(Operation::DeleteCandidateSet {
                set_name: "windows".into(),
            })
            .unwrap_err(),
    );
    unavailable(
        engine
            .apply(Operation::GenerateCandidateSet {
                set_name: "replacement".into(),
                seq_id: "source".into(),
                length_bp: 4,
                step_bp: 4,
                feature_kinds: vec![],
                feature_label_regex: None,
                max_distance_bp: None,
                feature_geometry_mode: None,
                feature_boundary_mode: None,
                feature_strand_relation: None,
                limit: Some(10),
            })
            .unwrap_err(),
    );
    assert_eq!(engine.state().metadata, before);
    assert_eq!(engine.state().sequences["source"].len(), 8);
    assert_eq!(engine.summarize_state().sequences.len(), 1);
    let graph = engine.project_fact_graph();
    assert_eq!(graph.unavailable_facts, ["candidate_set.exists"]);
    assert!(
        graph
            .facts
            .iter()
            .any(|fact| fact.fact == "sequence.exists")
    );
    assert!(
        !graph
            .facts
            .iter()
            .any(|fact| fact.fact == "candidate_set.exists")
    );
    let expression = FactExpression::Atom(FactAtom {
        fact: "candidate_set.exists".into(),
        id: Some("windows".into()),
        ..FactAtom::default()
    });
    assert_eq!(
        GentleEngine::evaluate_fact_expression_against_graph(&expression, &graph).truth,
        FactTruth::Unknown
    );
    assert_eq!(
        GentleEngine::default()
            .evaluate_fact_expression(&expression, &[])
            .truth,
        FactTruth::Unsatisfied
    );
    let wire: ProjectFactGraph = serde_json::from_value(json!({"facts":[]})).unwrap();
    assert!(
        wire.unavailable_facts.is_empty(),
        "legacy graph reading stays additive"
    );
    let context = crate::agent_bridge::build_agent_introspection_context(&graph);
    assert_eq!(
        serde_json::to_value(&wire).unwrap(),
        json!({"schema":PROJECT_FACT_GRAPH_SCHEMA, "facts":[]})
    );
    assert!(
        !context
            .fact_type_counts
            .contains_key("candidate_set.exists")
    );
    assert!(
        context
            .notes
            .iter()
            .any(|note| note.contains("candidate_set.exists: unavailable"))
    );
    let shell = crate::engine_shell::parse_shell_line("candidates list").unwrap();
    let err = crate::engine_shell::execute_shell_command(&mut engine, &shell).unwrap_err();
    assert!(err.contains("candidate_store_unavailable:"));
}

#[test]
fn candidate_store_unresolved_warning_only_and_unhydrated_json_cannot_be_silently_healed() {
    let _lock = candidate_store_env_lock().lock().unwrap();
    let _env = EnvVarGuard::set(CANDIDATE_STORE_STRICT_LOAD_ENV, "0");
    let td = tempdir().unwrap();
    let f = fixture(td.path(), "\n");
    let raw: ProjectState = serde_json::from_slice(&fs::read(&f.project).unwrap()).unwrap();
    unavailable(
        raw.save_to_path(td.path().join("new.json").to_str().unwrap())
            .unwrap_err(),
    );
    let mut legacy = ProjectState::default();
    let warning =
        json!({"message":"old missing sidecar", "strict_env": CANDIDATE_STORE_STRICT_LOAD_ENV});
    legacy.metadata.insert(
        CANDIDATE_SETS_LOAD_WARNING_METADATA_KEY.into(),
        warning.clone(),
    );
    fs::write(&f.project, serde_json::to_vec(&legacy).unwrap()).unwrap();
    let loaded = ProjectState::load_from_path(f.project.to_str().unwrap()).unwrap();
    assert_eq!(
        loaded.metadata[CANDIDATE_SETS_LOAD_WARNING_METADATA_KEY],
        warning
    );
    assert!(!loaded.metadata.contains_key(CANDIDATE_SETS_METADATA_KEY));
    let round_trip: ProjectState =
        serde_json::from_value(serde_json::to_value(loaded.clone()).unwrap()).unwrap();
    let before = tree(td.path());
    for state in [loaded, round_trip] {
        unavailable(state.save_to_path(f.project.to_str().unwrap()).unwrap_err());
        unavailable(
            GentleEngine::from_state(state)
                .list_candidate_sets()
                .unwrap_err(),
        );
    }
    assert_eq!(tree(td.path()), before);
}

#[test]
fn candidate_store_unresolved_recovery_requires_reopen_and_preserves_healthy_empty_controls() {
    let _lock = candidate_store_env_lock().lock().unwrap();
    let _env = EnvVarGuard::set(CANDIDATE_STORE_STRICT_LOAD_ENV, "0");
    let td = tempdir().unwrap();
    let f = fixture(td.path(), "\r\n");
    // Existing ancillary defaults remain accepted; no count-validation policy change.
    let mut index: serde_json::Value =
        serde_json::from_slice(&fs::read(&f.index).unwrap()).unwrap();
    index.as_object_mut().unwrap().remove("updated_at_unix_ms");
    for key in ["metrics", "created_at_unix_ms", "source_seq_ids"] {
        index["sets"][0].as_object_mut().unwrap().remove(key);
    }
    fs::write(&f.index, serde_json::to_vec(&index).unwrap()).unwrap();
    ProjectState::load_from_path(f.project.to_str().unwrap())
        .unwrap()
        .ensure_candidate_store_available()
        .unwrap();
    let original_index = fs::read(&f.index).unwrap();
    fs::write(&f.index, b"broken").unwrap();
    let degraded = ProjectState::load_from_path(f.project.to_str().unwrap()).unwrap();
    fs::write(&f.index, original_index).unwrap();
    unavailable(
        degraded
            .save_to_path(f.project.to_str().unwrap())
            .unwrap_err(),
    );
    let mut persisted = serde_json::to_value(&degraded).unwrap();
    // Synthetic previously retained warning + original reference: full hydration clears it.
    persisted["metadata"][CANDIDATE_SETS_METADATA_KEY] = f.reference;
    fs::write(&f.project, serde_json::to_vec(&persisted).unwrap()).unwrap();
    let healthy = ProjectState::load_from_path(f.project.to_str().unwrap()).unwrap();
    assert!(
        !healthy
            .metadata
            .contains_key(CANDIDATE_SETS_LOAD_WARNING_METADATA_KEY)
    );
    let target = td.path().join("other's folder/new.json");
    healthy.save_to_path(target.to_str().unwrap()).unwrap();
    let reopened = ProjectState::load_from_path(target.to_str().unwrap()).unwrap();
    assert_eq!(
        GentleEngine::from_state(reopened)
            .list_candidate_sets()
            .unwrap()[0]
            .candidate_count,
        1
    );
    let before_deletion = tree(td.path());
    let mut deleting = GentleEngine::from_state(healthy);
    deleting
        .apply(Operation::DeleteCandidateSet {
            set_name: "windows".into(),
        })
        .unwrap();
    deleting
        .state()
        .save_to_path(f.project.to_str().unwrap())
        .unwrap();
    assert_ne!(
        tree(td.path()),
        before_deletion,
        "negative control detects intentional deletion"
    );
    assert!(!ProjectState::candidate_store_sidecar_dir(&f.project).exists());
    let empty = ProjectState::default();
    assert!(
        GentleEngine::from_state(empty.clone())
            .list_candidate_sets()
            .unwrap()
            .is_empty()
    );
    empty.save_to_path(f.project.to_str().unwrap()).unwrap();
    assert!(
        !ProjectState::candidate_store_sidecar_dir(&f.project).exists(),
        "healthy intentional removal is unchanged"
    );
}

#[test]
fn candidate_store_unresolved_strict_load_still_refuses_missing_but_keeps_legacy_load_compatibility()
 {
    let _lock = candidate_store_env_lock().lock().unwrap();
    let _env = EnvVarGuard::set(CANDIDATE_STORE_STRICT_LOAD_ENV, "1");
    let td = tempdir().unwrap();
    let f = fixture(td.path(), "\n");
    fs::remove_file(&f.index).unwrap();
    assert!(
        ProjectState::load_from_path(f.project.to_str().unwrap())
            .unwrap_err()
            .message
            .contains("candidate-store index")
    );
    let mut legacy = ProjectState::default();
    legacy.metadata.insert(
        CANDIDATE_SETS_LOAD_WARNING_METADATA_KEY.into(),
        json!({"message":"legacy"}),
    );
    fs::write(&f.project, serde_json::to_vec(&legacy).unwrap()).unwrap();
    let loaded = ProjectState::load_from_path(f.project.to_str().unwrap()).unwrap();
    unavailable(
        loaded
            .save_to_path(f.project.to_str().unwrap())
            .unwrap_err(),
    );
}

#[test]
fn candidate_store_unresolved_unknown_schemas_load_but_cannot_be_used_or_saved() {
    let _lock = candidate_store_env_lock().lock().unwrap();
    let td = tempdir().unwrap();
    let project = td.path().join("unknown.json");
    for strict in ["0", "1"] {
        let _env = EnvVarGuard::set(CANDIDATE_STORE_STRICT_LOAD_ENV, strict);
        for value in [
            json!({"schema":"synthetic.future_reference.v9", "index_path":"not-read/index.json"}),
            json!({"schema":"synthetic.future_inline.v9", "sets":{}}),
        ] {
            let mut state = ProjectState::default();
            state
                .metadata
                .insert(CANDIDATE_SETS_METADATA_KEY.into(), value.clone());
            fs::write(&project, serde_json::to_vec(&state).unwrap()).unwrap();
            let before = tree(td.path());
            let loaded = ProjectState::load_from_path(project.to_str().unwrap()).unwrap();
            assert_eq!(loaded.metadata[CANDIDATE_SETS_METADATA_KEY], value);
            unavailable(loaded.save_to_path(project.to_str().unwrap()).unwrap_err());
            unavailable(
                GentleEngine::from_state(loaded)
                    .list_candidate_sets()
                    .unwrap_err(),
            );
            assert_eq!(tree(td.path()), before);
        }
    }
}

#[test]
fn candidate_store_unresolved_status_survives_history_snapshots_and_serialization() {
    let mut state = ProjectState::default();
    state.metadata.insert(
        CANDIDATE_SETS_LOAD_WARNING_METADATA_KEY.into(),
        json!({"message":"synthetic legacy warning"}),
    );
    let mut engine = GentleEngine::from_state(state);
    engine
        .apply(Operation::CreateSequenceFromText {
            sequence_text: "ACGT".into(),
            output_id: Some("unrelated_synthetic".into()),
            name: None,
            circular: false,
        })
        .unwrap();
    unavailable(engine.list_candidate_sets().unwrap_err());
    engine.undo_last_operation().unwrap();
    unavailable(engine.list_candidate_sets().unwrap_err());
    engine.redo_last_operation().unwrap();
    let fork = engine.fork_detached_execution();
    unavailable(fork.engine().list_candidate_sets().unwrap_err());
    let mut imported = engine.fork_detached_execution();
    imported.import_checkpoint_engine(engine.clone());
    unavailable(imported.engine().list_candidate_sets().unwrap_err());
    let round_trip: ProjectState =
        serde_json::from_value(serde_json::to_value(imported.engine().state()).unwrap()).unwrap();
    unavailable(
        GentleEngine::from_state(round_trip)
            .list_candidate_sets()
            .unwrap_err(),
    );
}

#[test]
fn candidate_store_unresolved_does_not_tighten_existing_index_count_defaults() {
    let _lock = candidate_store_env_lock().lock().unwrap();
    let _env = EnvVarGuard::set(CANDIDATE_STORE_STRICT_LOAD_ENV, "1");
    let td = tempdir().unwrap();
    let f = fixture(td.path(), "\n");
    let mut index: serde_json::Value =
        serde_json::from_slice(&fs::read(&f.index).unwrap()).unwrap();
    index.as_object_mut().unwrap().remove("schema");
    index.as_object_mut().unwrap().remove("set_count");
    index["sets"][0]
        .as_object_mut()
        .unwrap()
        .remove("candidate_count");
    fs::write(&f.index, serde_json::to_vec(&index).unwrap()).unwrap();
    let loaded = ProjectState::load_from_path(f.project.to_str().unwrap()).unwrap();
    assert_eq!(
        GentleEngine::from_state(loaded)
            .list_candidate_sets()
            .unwrap()[0]
            .candidate_count,
        1
    );
}
