//! Optimistic detached execution for long-running GUI engine operations.
//!
//! Workers clone a versioned engine snapshot under a read lock, execute without
//! holding the shared lock, and commit with a short guarded swap. A concurrent
//! edit makes the result stale and leaves the live project untouched.

use crate::engine::{EngineError, ErrorCode, GentleEngine};
use std::sync::{Arc, RwLock};

pub(crate) fn execute_on_engine_snapshot<T>(
    shared: &Arc<RwLock<GentleEngine>>,
    work: impl FnOnce(&mut GentleEngine) -> Result<T, EngineError>,
) -> Result<T, EngineError> {
    let (base_structural_revision, base_mutation_revision, base_journal_len, mut snapshot) = {
        let guard = shared.read().map_err(|_| EngineError {
            code: ErrorCode::Internal,
            message: "Engine lock poisoned while snapshotting background work".to_string(),
            cause_chain: vec![],
        })?;
        (
            guard.structural_revision(),
            guard.mutation_revision(),
            guard.operation_log().len(),
            guard.clone(),
        )
    };
    let base_metadata = snapshot.state().metadata.clone();

    let result = work(&mut snapshot)?;
    let structural_changed = snapshot.structural_revision() != base_structural_revision;
    let persisted_state_changed = snapshot.mutation_revision() != base_mutation_revision;
    let old_snapshot = {
        let mut guard = shared.write().map_err(|_| EngineError {
            code: ErrorCode::Internal,
            message: "Engine lock poisoned while committing background work".to_string(),
            cause_chain: vec![],
        })?;
        if guard.structural_revision() != base_structural_revision
            || guard.operation_log().len() != base_journal_len
        {
            return Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: "Background result became stale because project data or operation history changed while it was running; rerun the operation"
                    .to_string(),
                cause_chain: vec![],
            });
        }
        snapshot.rebase_detached_commit_revisions(
            &guard,
            &base_metadata,
            base_journal_len,
            structural_changed,
            persisted_state_changed,
        )?;
        std::mem::replace(&mut *guard, snapshot)
    };
    drop(old_snapshot);
    Ok(result)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::dna_sequence::DNAsequence;
    use crate::engine::{Engine, Operation, Workflow};

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
}
