//! Optimistic detached execution for long-running GUI engine operations.
//!
//! Workers fork a versioned engine execution snapshot under a read lock,
//! execute without holding the shared lock, and commit with a short guarded
//! swap. Inherited undo/redo snapshots stay in the live engine. A concurrent
//! edit makes the result stale and leaves the live project untouched.

use crate::engine::{Engine, EngineError, ErrorCode, GentleEngine, OpResult, Operation};
use std::sync::{Arc, RwLock};

pub(crate) fn execute_on_engine_snapshot<T>(
    shared: &Arc<RwLock<GentleEngine>>,
    work: impl FnOnce(&mut GentleEngine) -> Result<T, EngineError>,
) -> Result<T, EngineError> {
    let mut detached = {
        let guard = shared.read().map_err(|_| EngineError {
            code: ErrorCode::Internal,
            message: "Engine lock poisoned while snapshotting background work".to_string(),
            cause_chain: vec![],
        })?;
        guard.fork_detached_execution()
    };

    let result = work(detached.engine_mut())?;
    let old_snapshot = {
        let mut guard = shared.write().map_err(|_| EngineError {
            code: ErrorCode::Internal,
            message: "Engine lock poisoned while committing background work".to_string(),
            cause_chain: vec![],
        })?;
        guard.commit_detached_execution(&mut detached)
    }?;
    drop(old_snapshot);
    Ok(result)
}

/// Apply a portable inspection operation to a detached snapshot and discard
/// its private journal/state after returning the typed result.
pub(crate) fn execute_read_only_operation_on_engine_snapshot(
    shared: &Arc<RwLock<GentleEngine>>,
    operation: Operation,
) -> Result<OpResult, EngineError> {
    let mut detached = {
        let guard = shared.read().map_err(|_| EngineError {
            code: ErrorCode::Internal,
            message: "Engine lock poisoned while snapshotting a read-only operation".to_string(),
            cause_chain: vec![],
        })?;
        guard.fork_detached_execution()
    };
    detached.engine_mut().apply(operation)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::dna_sequence::DNAsequence;
    use crate::engine::{Engine, GenomeTrackSource, GenomeTrackSubscription, Operation, Workflow};

    #[test]
    fn detached_fork_excludes_inherited_history_but_keeps_execution_baseline() {
        let mut live = GentleEngine::new();
        live.apply(Operation::CreateSequenceFromText {
            sequence_text: "ATGC".to_string(),
            output_id: Some("existing".to_string()),
            name: None,
            circular: false,
        })
        .expect("create existing sequence");
        let undone_result = live
            .apply(Operation::CreateSequenceFromText {
                sequence_text: "GCTA".to_string(),
                output_id: Some("undone".to_string()),
                name: None,
                circular: false,
            })
            .expect("create sequence to undo");
        live.undo_last_operation().expect("create live redo entry");

        assert_eq!(live.undo_available(), 1);
        assert_eq!(live.redo_available(), 1);
        let baseline_execution_revision = live.execution_revision();
        let baseline_mutation_revision = live.mutation_revision();
        let baseline_structural_revision = live.structural_revision();
        let baseline_journal_len = live.operation_log().len();
        let baseline_history_limit = live.history_summary().history_limit;

        let mut detached = live.fork_detached_execution();
        let snapshot = detached.engine();
        assert!(snapshot.state().sequences.contains_key("existing"));
        assert!(!snapshot.state().sequences.contains_key("undone"));
        assert_eq!(snapshot.operation_log().len(), baseline_journal_len);
        assert_eq!(snapshot.execution_revision(), baseline_execution_revision);
        assert_eq!(snapshot.mutation_revision(), baseline_mutation_revision);
        assert_eq!(snapshot.structural_revision(), baseline_structural_revision);
        assert_eq!(
            snapshot.history_summary().history_limit,
            baseline_history_limit
        );
        assert_eq!(snapshot.undo_available(), 0);
        assert_eq!(snapshot.redo_available(), 0);

        let detached_result = detached
            .engine_mut()
            .apply(Operation::CreateSequenceFromText {
                sequence_text: "TTAA".to_string(),
                output_id: Some("detached".to_string()),
                name: None,
                circular: false,
            })
            .expect("execute from inherited operation counter");
        assert_eq!(detached_result.op_id, undone_result.op_id);
    }

    #[test]
    fn detached_read_only_snapshot_excludes_inherited_history() {
        let mut live = GentleEngine::new();
        live.apply(Operation::CreateSequenceFromText {
            sequence_text: "ATGC".to_string(),
            output_id: Some("existing".to_string()),
            name: None,
            circular: false,
        })
        .expect("create existing sequence");
        live.apply(Operation::CreateSequenceFromText {
            sequence_text: "GCTA".to_string(),
            output_id: Some("redoable".to_string()),
            name: None,
            circular: false,
        })
        .expect("create redoable sequence");
        live.undo_last_operation().expect("prepare redo history");

        let snapshot = live.clone_without_history();
        assert_eq!(snapshot.undo_available(), 0);
        assert_eq!(snapshot.redo_available(), 0);
        assert_eq!(
            serde_json::to_value(snapshot.operation_log()).expect("snapshot journal"),
            serde_json::to_value(live.operation_log()).expect("live journal")
        );
        assert_eq!(snapshot.execution_revision(), live.execution_revision());
        assert_eq!(snapshot.mutation_revision(), live.mutation_revision());
        assert_eq!(snapshot.structural_revision(), live.structural_revision());
        assert_eq!(
            serde_json::to_value(snapshot.state()).expect("snapshot state"),
            serde_json::to_value(live.state()).expect("live state")
        );
        assert!(snapshot.state().sequences.contains_key("existing"));
        assert!(!snapshot.state().sequences.contains_key("redoable"));
    }

    #[test]
    fn read_only_operation_returns_data_without_committing_state() {
        let sequence = "AAAGTCCCCCTACTAACCCCCCCCCCCCCCCCCCCAGAAA";
        let mut live = GentleEngine::new();
        live.apply(Operation::CreateSequenceFromText {
            sequence_text: sequence.to_string(),
            output_id: Some("existing".to_string()),
            name: None,
            circular: false,
        })
        .expect("create existing sequence");
        let revision_before = live.mutation_revision();
        let journal_len_before = live.operation_log().len();
        let shared = Arc::new(RwLock::new(live));
        let result = execute_read_only_operation_on_engine_snapshot(
            &shared,
            Operation::InspectCrypticSplicingScreen {
                request: gentle_protocol::CrypticSplicingScreenRequest {
                    seq_id: "existing".to_string(),
                    start_1based: 1,
                    end_1based: sequence.len(),
                    min_pseudo_intron_bp: 20,
                    max_pseudo_intron_bp: 100,
                    max_candidate_pairs: 10,
                    ..gentle_protocol::CrypticSplicingScreenRequest::default()
                },
                path: None,
            },
        )
        .expect("read-only operation");
        assert!(result.cryptic_splicing_screen.is_some());
        let guard = shared.read().expect("live engine");
        assert_eq!(guard.mutation_revision(), revision_before);
        assert_eq!(guard.operation_log().len(), journal_len_before);
    }

    #[test]
    fn read_only_worker_sources_use_history_free_engine_clones() {
        let worker_sources = [
            include_str!("main_area_dna/rna_read_mapping_ui.rs"),
            include_str!("app/routine_and_agent_assistant_ui.rs"),
        ];
        for source in worker_sources {
            assert!(
                !source.contains("guard.clone()"),
                "read-only background workers must not clone complete engine history"
            );
            assert!(
                source.contains("guard.clone_without_history()"),
                "read-only worker should use the explicit history-free clone"
            );
        }
    }

    #[test]
    fn collection_construct_reasoning_inspection_excludes_inherited_history() {
        let source = include_str!("engine/ops/operation_handlers.rs");
        assert!(
            !source.contains("let mut inspection_engine = self.clone();"),
            "collection inspection must not clone complete engine history"
        );
        assert!(
            source.contains("let mut inspection_engine = self.clone_without_history();"),
            "collection inspection should use the explicit history-free clone"
        );
    }

    #[test]
    fn detached_commit_publishes_snapshot_with_revisions() {
        let shared = Arc::new(RwLock::new(GentleEngine::new()));
        execute_on_engine_snapshot(&shared, |snapshot| {
            snapshot
                .state_mut()
                .metadata
                .insert("finished".to_string(), serde_json::Value::Bool(true));
            Ok(())
        })
        .expect("detached commit");

        let guard = shared.read().expect("engine");
        assert_eq!(
            guard.state().metadata.get("finished"),
            Some(&serde_json::Value::Bool(true))
        );
        assert!(guard.execution_revision() > 0);
        assert!(guard.mutation_revision() > 0);
    }

    #[test]
    fn detached_track_metadata_commit_marks_dirty_and_invalidates_live_redo() {
        let shared = Arc::new(RwLock::new(GentleEngine::new()));
        {
            let mut live = shared.write().expect("engine");
            live.apply(Operation::CreateSequenceFromText {
                sequence_text: "ATGC".to_string(),
                output_id: Some("existing".to_string()),
                name: None,
                circular: false,
            })
            .expect("create existing sequence");
            live.apply(Operation::CreateSequenceFromText {
                sequence_text: "GCTA".to_string(),
                output_id: Some("redoable".to_string()),
                name: None,
                circular: false,
            })
            .expect("create redoable sequence");
            live.undo_last_operation().expect("prepare redo history");
            assert_eq!(live.undo_available(), 1);
            assert_eq!(live.redo_available(), 1);
        }
        let baseline_mutation_revision = shared.read().expect("engine").mutation_revision();

        execute_on_engine_snapshot(&shared, |snapshot| {
            assert!(
                snapshot.add_genome_track_subscription(GenomeTrackSubscription {
                    source: GenomeTrackSource::Bed,
                    path: "tracked.bed".to_string(),
                    track_name: Some("tracked".to_string()),
                    min_score: None,
                    max_score: None,
                    clear_existing: false,
                })?
            );
            Ok(())
        })
        .expect("commit metadata-only detached work");

        let live = shared.write().expect("engine");
        assert!(live.mutation_revision() > baseline_mutation_revision);
        assert_eq!(live.undo_available(), 1);
        assert_eq!(live.redo_available(), 0);
        assert_eq!(live.list_genome_track_subscriptions().len(), 1);
        assert!(!live.state().sequences.contains_key("redoable"));
    }

    #[test]
    fn clearing_track_subscriptions_advances_mutation_revision() {
        let mut engine = GentleEngine::new();
        engine
            .add_genome_track_subscription(GenomeTrackSubscription {
                source: GenomeTrackSource::Bed,
                path: "tracked.bed".to_string(),
                track_name: Some("tracked".to_string()),
                min_score: None,
                max_score: None,
                clear_existing: false,
            })
            .expect("add subscription");
        let before_clear = engine.mutation_revision();

        engine.clear_genome_track_subscriptions();

        assert!(engine.mutation_revision() > before_clear);
        assert!(engine.list_genome_track_subscriptions().is_empty());
    }

    #[test]
    fn detached_commit_rejects_intervening_project_mutation() {
        let shared = Arc::new(RwLock::new(GentleEngine::new()));
        let worker_engine = shared.clone();
        let (started_tx, started_rx) = std::sync::mpsc::channel();
        let (continue_tx, continue_rx) = std::sync::mpsc::channel();
        let worker = std::thread::spawn(move || {
            execute_on_engine_snapshot(&worker_engine, |snapshot| {
                started_tx.send(()).expect("signal snapshot");
                continue_rx.recv().expect("continue worker");
                snapshot
                    .state_mut()
                    .metadata
                    .insert("worker".to_string(), serde_json::Value::Bool(true));
                Ok(())
            })
        });

        started_rx.recv().expect("worker started");
        shared
            .write()
            .expect("engine")
            .state_mut()
            .sequences
            .insert(
                "live".to_string(),
                DNAsequence::from_sequence("ATGC").expect("sequence"),
            );
        continue_tx.send(()).expect("continue");

        let error = worker.join().expect("worker join").expect_err("stale");
        assert!(error.message.contains("became stale"));
        let guard = shared.read().expect("engine");
        assert!(guard.state().sequences.contains_key("live"));
        assert!(!guard.state().metadata.contains_key("worker"));
    }

    #[test]
    fn detached_commit_preserves_concurrent_display_changes() {
        let shared = Arc::new(RwLock::new(GentleEngine::new()));
        let worker_engine = shared.clone();
        let (started_tx, started_rx) = std::sync::mpsc::channel();
        let (continue_tx, continue_rx) = std::sync::mpsc::channel();
        let worker = std::thread::spawn(move || {
            execute_on_engine_snapshot(&worker_engine, |snapshot| {
                started_tx.send(()).expect("signal snapshot");
                continue_rx.recv().expect("continue worker");
                snapshot
                    .state_mut()
                    .metadata
                    .insert("worker".to_string(), serde_json::Value::Bool(true));
                Ok(())
            })
        });

        started_rx.recv().expect("worker started");
        shared
            .write()
            .expect("engine")
            .display_state_mut()
            .linear_view_start_bp = 123;
        continue_tx.send(()).expect("continue");
        worker.join().expect("worker join").expect("commit");

        let guard = shared.read().expect("engine");
        assert_eq!(guard.state().display.linear_view_start_bp, 123);
        assert_eq!(
            guard.state().metadata.get("worker"),
            Some(&serde_json::Value::Bool(true))
        );
    }

    #[test]
    fn detached_commit_preserves_concurrent_auxiliary_metadata_through_undo() {
        let shared = Arc::new(RwLock::new(GentleEngine::new()));
        let worker_engine = shared.clone();
        let (started_tx, started_rx) = std::sync::mpsc::channel();
        let (continue_tx, continue_rx) = std::sync::mpsc::channel();
        let worker = std::thread::spawn(move || {
            execute_on_engine_snapshot(&worker_engine, |snapshot| {
                started_tx.send(()).expect("signal snapshot");
                continue_rx.recv().expect("continue worker");
                snapshot.apply(Operation::CreateSequenceFromText {
                    sequence_text: "ATGC".to_string(),
                    output_id: Some("worker".to_string()),
                    name: None,
                    circular: false,
                })?;
                Ok(())
            })
        });

        started_rx.recv().expect("worker started");
        shared
            .write()
            .expect("engine")
            .auxiliary_metadata_mut()
            .insert("workspace".to_string(), serde_json::Value::Bool(true));
        continue_tx.send(()).expect("continue");
        worker.join().expect("worker join").expect("commit");

        let mut guard = shared.write().expect("engine");
        assert!(guard.state().sequences.contains_key("worker"));
        assert_eq!(
            guard.state().metadata.get("workspace"),
            Some(&serde_json::Value::Bool(true))
        );
        guard
            .undo_last_operation()
            .expect("undo detached operation");
        assert!(!guard.state().sequences.contains_key("worker"));
        assert_eq!(
            guard.state().metadata.get("workspace"),
            Some(&serde_json::Value::Bool(true))
        );
    }

    #[test]
    fn detached_commit_rejects_conflicting_auxiliary_metadata() {
        let shared = Arc::new(RwLock::new(GentleEngine::new()));
        let worker_engine = shared.clone();
        let (started_tx, started_rx) = std::sync::mpsc::channel();
        let (continue_tx, continue_rx) = std::sync::mpsc::channel();
        let worker = std::thread::spawn(move || {
            execute_on_engine_snapshot(&worker_engine, |snapshot| {
                started_tx.send(()).expect("signal snapshot");
                continue_rx.recv().expect("continue worker");
                snapshot
                    .state_mut()
                    .metadata
                    .insert("shared".to_string(), serde_json::json!("worker"));
                Ok(())
            })
        });

        started_rx.recv().expect("worker started");
        shared
            .write()
            .expect("engine")
            .auxiliary_metadata_mut()
            .insert("shared".to_string(), serde_json::json!("live"));
        continue_tx.send(()).expect("continue");

        let error = worker.join().expect("worker join").expect_err("stale");
        assert!(error.message.contains("metadata key 'shared'"));
        assert_eq!(
            shared
                .read()
                .expect("engine")
                .state()
                .metadata
                .get("shared"),
            Some(&serde_json::json!("live"))
        );
    }

    #[test]
    fn detached_workflow_preserves_live_display_in_each_new_undo_checkpoint() {
        let shared = Arc::new(RwLock::new(GentleEngine::new()));
        let worker_engine = shared.clone();
        let (started_tx, started_rx) = std::sync::mpsc::channel();
        let (continue_tx, continue_rx) = std::sync::mpsc::channel();
        let worker = std::thread::spawn(move || {
            execute_on_engine_snapshot(&worker_engine, |snapshot| {
                started_tx.send(()).expect("signal snapshot");
                continue_rx.recv().expect("continue worker");
                snapshot.apply_workflow(Workflow {
                    run_id: "detached-batch".to_string(),
                    ops: vec![
                        Operation::CreateSequenceFromText {
                            sequence_text: "ATGC".to_string(),
                            output_id: Some("first".to_string()),
                            name: None,
                            circular: false,
                        },
                        Operation::CreateSequenceFromText {
                            sequence_text: "GCTA".to_string(),
                            output_id: Some("second".to_string()),
                            name: None,
                            circular: false,
                        },
                    ],
                })?;
                Ok(())
            })
        });

        started_rx.recv().expect("worker started");
        shared
            .write()
            .expect("engine")
            .display_state_mut()
            .linear_view_start_bp = 123;
        continue_tx.send(()).expect("continue");
        worker.join().expect("worker join").expect("commit");

        let mut guard = shared.write().expect("engine");
        guard.undo_last_operation().expect("undo second operation");
        assert_eq!(guard.state().display.linear_view_start_bp, 123);
        guard.undo_last_operation().expect("undo first operation");
        assert_eq!(guard.state().display.linear_view_start_bp, 123);
    }

    #[test]
    fn detached_commit_appends_new_checkpoints_after_live_history() {
        let shared = Arc::new(RwLock::new(GentleEngine::new()));
        shared
            .write()
            .expect("engine")
            .apply(Operation::CreateSequenceFromText {
                sequence_text: "ATGC".to_string(),
                output_id: Some("existing".to_string()),
                name: None,
                circular: false,
            })
            .expect("create existing sequence");

        execute_on_engine_snapshot(&shared, |snapshot| {
            snapshot.apply_workflow(Workflow {
                run_id: "detached-history".to_string(),
                ops: vec![
                    Operation::CreateSequenceFromText {
                        sequence_text: "GCTA".to_string(),
                        output_id: Some("detached_first".to_string()),
                        name: None,
                        circular: false,
                    },
                    Operation::CreateSequenceFromText {
                        sequence_text: "TTAA".to_string(),
                        output_id: Some("detached_second".to_string()),
                        name: None,
                        circular: false,
                    },
                ],
            })?;
            Ok(())
        })
        .expect("commit detached workflow");

        let mut guard = shared.write().expect("engine");
        assert_eq!(guard.undo_available(), 3);
        assert_eq!(guard.redo_available(), 0);
        assert_eq!(guard.operation_log().len(), 3);

        guard.undo_last_operation().expect("undo detached second");
        assert!(!guard.state().sequences.contains_key("detached_second"));
        assert!(guard.state().sequences.contains_key("detached_first"));
        assert!(guard.state().sequences.contains_key("existing"));

        guard.undo_last_operation().expect("undo detached first");
        assert!(!guard.state().sequences.contains_key("detached_first"));
        assert!(guard.state().sequences.contains_key("existing"));

        guard
            .undo_last_operation()
            .expect("undo pre-existing operation");
        assert!(!guard.state().sequences.contains_key("existing"));
        assert_eq!(guard.undo_available(), 0);
        assert_eq!(guard.redo_available(), 3);
    }

    #[test]
    fn detached_read_only_work_does_not_rewrite_an_older_undo_checkpoint() {
        let shared = Arc::new(RwLock::new(GentleEngine::new()));
        shared
            .write()
            .expect("engine")
            .apply(Operation::CreateSequenceFromText {
                sequence_text: "ATGC".to_string(),
                output_id: Some("existing".to_string()),
                name: None,
                circular: false,
            })
            .expect("create sequence");
        shared
            .write()
            .expect("engine")
            .display_state_mut()
            .linear_view_start_bp = 123;

        execute_on_engine_snapshot(&shared, |snapshot| {
            snapshot.apply(Operation::ListJasparCatalog {
                filter: None,
                limit: Some(1),
                include_remote_metadata: false,
                refresh_remote_metadata: false,
                path: None,
            })?;
            Ok(())
        })
        .expect("read-only detached operation");

        let mut guard = shared.write().expect("engine");
        guard.undo_last_operation().expect("undo sequence creation");
        assert_eq!(guard.state().display.linear_view_start_bp, 0);
    }

    #[test]
    fn detached_read_only_operation_preserves_undo_and_redo_availability() {
        let shared = Arc::new(RwLock::new(GentleEngine::new()));
        {
            let mut guard = shared.write().expect("engine");
            guard
                .apply(Operation::CreateSequenceFromText {
                    sequence_text: "ATGC".to_string(),
                    output_id: Some("redoable".to_string()),
                    name: None,
                    circular: false,
                })
                .expect("create sequence");
            guard.undo_last_operation().expect("prepare redo history");
            assert_eq!(guard.undo_available(), 0);
            assert_eq!(guard.redo_available(), 1);
        }

        execute_on_engine_snapshot(&shared, |snapshot| {
            snapshot.apply(Operation::ListJasparCatalog {
                filter: None,
                limit: Some(1),
                include_remote_metadata: false,
                refresh_remote_metadata: false,
                path: None,
            })?;
            Ok(())
        })
        .expect("read-only detached operation");

        let mut guard = shared.write().expect("engine");
        assert_eq!(guard.undo_available(), 0);
        assert_eq!(guard.redo_available(), 1);
        guard
            .redo_last_operation()
            .expect("redo pre-existing operation");
        assert!(guard.state().sequences.contains_key("redoable"));
    }
}
