//! Inline synthetic files and tutorial 06.07's invented CDS; no private inputs.
//! Recreated in temporary directories; Rust/native execution belongs to Glen/CI.

use super::*;
use crate::engine::{Engine, Workflow, lifecycle_slice3_tests as fixture};
use serde_json::json;

fn prepare(engine: &GentleEngine, path: &Path) -> Result<PreviewDestination, EngineError> {
    engine.sequence_design_preview_destination(path.to_str().unwrap())
}

fn bound(path: &Path) -> GentleEngine {
    let mut engine = GentleEngine::new();
    engine.set_host_project_path(Some(path)).unwrap();
    engine
}

#[test]
fn lifecycle_slice4_fresh_exact_json_and_all_existing_destinations_refused() {
    let td = tempfile::tempdir().unwrap();
    let engine = GentleEngine::new();
    let value = json!({"synthetic":[1,"two"],"absent":null});
    let path = td.path().join("preview.json");
    prepare(&engine, &path).unwrap().write(&value).unwrap();
    let expected = serde_json::to_string_pretty(&value).unwrap().into_bytes();
    assert_eq!(fs::read(&path).unwrap(), expected);
    assert!(
        prepare(&engine, &path)
            .unwrap_err()
            .message
            .contains("destination_exists")
    );
    assert_eq!(fs::read(&path).unwrap(), expected);
    let empty = td.path().join("empty.json");
    fs::write(&empty, b"").unwrap();
    assert!(prepare(&engine, &empty).is_err());
    assert_eq!(fs::metadata(&empty).unwrap().len(), 0);
    let directory = td.path().join("directory");
    fs::create_dir(&directory).unwrap();
    assert!(prepare(&engine, &directory).is_err());
    let link = td.path().join("hardlink.json");
    fs::hard_link(&path, &link).unwrap();
    let identity = directory_identity(&link).unwrap();
    assert!(prepare(&engine, &link).is_err());
    assert_eq!(directory_identity(&path).unwrap(), identity);
    assert_eq!(fs::read(&link).unwrap(), expected);
}

#[test]
fn lifecycle_slice4_missing_parent_and_invalid_leaf_never_create_directories() {
    let td = tempfile::tempdir().unwrap();
    let engine = GentleEngine::new();
    let missing = td.path().join("missing");
    assert!(prepare(&engine, &missing.join("preview.json")).is_err());
    assert!(!missing.exists());
    for raw in [
        "",
        " ",
        ".",
        "..",
        "invalid\0file",
        "parent/",
        "parent/.",
        "parent/..",
    ] {
        assert!(
            engine.sequence_design_preview_destination(raw).is_err(),
            "{raw:?}"
        );
    }
    let inside = td.path().join("inside");
    fs::create_dir(&inside).unwrap();
    let path = inside.join("..").join("resolvable.json");
    prepare(&engine, &path)
        .unwrap()
        .write(&json!({"synthetic":true}))
        .unwrap();
    assert!(td.path().join("resolvable.json").is_file());
}

#[test]
fn lifecycle_slice4_bound_missing_project_case_and_save_owned_subtree_are_protected() {
    let td = tempfile::tempdir().unwrap();
    let project = td.path().join("project.json");
    let engine = bound(&project);
    for leaf in [
        "project.json",
        "PROJECT.JSON",
        "project.json.candidates",
        "PROJECT.JSON.CANDIDATES",
    ] {
        let path = td.path().join(leaf);
        assert!(
            prepare(&engine, &path)
                .unwrap_err()
                .message
                .contains("protected_path")
        );
        assert!(!path.exists());
    }
    let sidecar = ProjectState::candidate_store_sidecar_dir(&project);
    fs::create_dir(&sidecar).unwrap();
    let nested = sidecar.join("records");
    fs::create_dir(&nested).unwrap();
    let before = b"synthetic candidate records\r\n";
    let existing = nested.join("records.jsonl");
    fs::write(&existing, before).unwrap();
    for path in [sidecar.join("new.json"), nested.join("unused.json")] {
        assert!(
            prepare(&engine, &path)
                .unwrap_err()
                .message
                .contains("protected_path")
        );
        assert!(!path.exists());
    }
    assert_eq!(fs::read(existing).unwrap(), before);
    // Nonstandard external folders are not the save-owned derived subtree.
    let external = td.path().join("external");
    fs::create_dir(&external).unwrap();
    prepare(&engine, &external.join("new.json"))
        .unwrap()
        .write(&json!({"synthetic":true}))
        .unwrap();
    // A genuinely absent project parent is not an unavailable-identity error.
    let elsewhere = bound(&td.path().join("absent").join("project.json"));
    prepare(&elsewhere, &td.path().join("elsewhere.json")).unwrap();
}

#[test]
fn lifecycle_slice4_unicode_parents_work_but_ambiguous_protected_names_fail_closed() {
    let td = tempfile::tempdir().unwrap();
    let parent = td.path().join("parent's \u{e9}");
    fs::create_dir(&parent).unwrap();
    let engine = bound(&parent.join("project.json"));
    prepare(&engine, &parent.join("preview.json"))
        .unwrap()
        .write(&json!({"synthetic":true}))
        .unwrap();
    for name in ["pr\u{e9}view.json", "pre\u{301}view.json"] {
        assert!(
            prepare(&engine, &parent.join(name))
                .unwrap_err()
                .message
                .contains("unsupported_name_comparison")
        );
    }
    let unicode_project = bound(&parent.join("proje\u{301}ct.json"));
    assert!(prepare(&unicode_project, &parent.join("ascii.json")).is_err());
    let outside = td.path().join("outside");
    fs::create_dir(&outside).unwrap();
    prepare(&unicode_project, &outside.join("preview.json")).unwrap();
    // Pure state-only use has no filename to infer or protect.
    prepare(&GentleEngine::new(), &parent.join("pr\u{e9}view.json")).unwrap();
}

#[cfg(unix)]
#[test]
fn lifecycle_slice4_symlinks_dangling_links_and_parent_aliases() {
    use std::os::unix::fs::symlink;
    let td = tempfile::tempdir().unwrap();
    let target = td.path().join("sentinel");
    fs::write(&target, b"unchanged").unwrap();
    let link = td.path().join("link.json");
    symlink(&target, &link).unwrap();
    assert!(prepare(&GentleEngine::new(), &link).is_err());
    let dangling = td.path().join("dangling.json");
    symlink(td.path().join("absent"), &dangling).unwrap();
    assert!(prepare(&GentleEngine::new(), &dangling).is_err());
    assert!(
        fs::symlink_metadata(&dangling)
            .unwrap()
            .file_type()
            .is_symlink()
    );
    assert_eq!(fs::read(target).unwrap(), b"unchanged");
    let socket_path = td.path().join("socket.json");
    let _socket = std::os::unix::net::UnixListener::bind(&socket_path).unwrap();
    assert!(prepare(&GentleEngine::new(), &socket_path).is_err());
    let real = td.path().join("real");
    fs::create_dir(&real).unwrap();
    let alias = td.path().join("alias");
    symlink(&real, &alias).unwrap();
    let project = real.join("project.json");
    let engine = bound(&project);
    assert!(prepare(&engine, &alias.join("PROJECT.JSON")).is_err());
    let sidecar = ProjectState::candidate_store_sidecar_dir(&project);
    fs::create_dir(&sidecar).unwrap();
    let external_alias = td.path().join("sidecar-alias");
    symlink(&sidecar, &external_alias).unwrap();
    assert!(prepare(&engine, &external_alias.join("new.json")).is_err());
    prepare(&engine, &alias.join("fresh.json"))
        .unwrap()
        .write(&json!({"synthetic":true}))
        .unwrap();
    assert!(real.join("fresh.json").is_file());
}

#[cfg(windows)]
#[test]
fn lifecycle_slice4_windows_unsupported_native_names_are_refused() {
    for raw in [
        r"C:preview.json",
        r"\preview.json",
        r"\\.\NUL",
        r"C:\preview.json:stream",
        r"C:\NUL.json",
        r"C:\COM1.json",
        r"C:\preview.json.",
    ] {
        assert!(
            GentleEngine::new()
                .sequence_design_preview_destination(raw)
                .unwrap_err()
                .message
                .contains("unsupported_path"),
            "{raw}"
        );
    }
}

#[test]
fn lifecycle_slice4_final_name_race_does_not_replace_competitor() {
    let td = tempfile::tempdir().unwrap();
    let path = td.path().join("race.json");
    let destination = prepare(&GentleEngine::new(), &path).unwrap();
    let error = destination
        .write_with(&json!({"synthetic":true}), |file, text| {
            stage(file, text)?;
            fs::write(&path, b"competitor").unwrap();
            Ok(())
        })
        .unwrap_err();
    assert!(error.message.contains("sequence_design_report.publish"));
    assert_eq!(fs::read(&path).unwrap(), b"competitor");
    assert_eq!(
        fs::read_dir(td.path()).unwrap().count(),
        1,
        "owned staging file cleaned up"
    );
}

#[test]
fn lifecycle_slice4_staged_faults_cleanup_only_owned_temporary_files() {
    for reason in ["write", "flush", "sync"] {
        let td = tempfile::tempdir().unwrap();
        let path = td.path().join("failed.json");
        let destination = prepare(&GentleEngine::new(), &path).unwrap();
        let sentinel = td.path().join("keep");
        fs::write(&sentinel, b"unchanged").unwrap();
        let error = destination
            .write_with(&json!({"synthetic":true}), |file, text| {
                file.write_all(&text.as_bytes()[..text.len() / 2]).unwrap();
                if reason != "write" {
                    file.flush().unwrap();
                }
                Err(io_error(
                    reason,
                    io::Error::other("synthetic injected staging failure"),
                ))
            })
            .unwrap_err();
        assert!(
            error
                .message
                .starts_with(&format!("sequence_design_report.{reason}:"))
        );
        assert!(!path.exists());
        assert_eq!(fs::read(sentinel).unwrap(), b"unchanged");
        assert_eq!(fs::read_dir(td.path()).unwrap().count(), 1);
    }
}

#[test]
fn lifecycle_slice4_serialization_failure_leaves_no_staging_or_final_file() {
    struct InvalidSynthetic;
    impl Serialize for InvalidSynthetic {
        fn serialize<S: serde::Serializer>(&self, _: S) -> Result<S::Ok, S::Error> {
            Err(serde::ser::Error::custom("synthetic serialization failure"))
        }
    }
    let td = tempfile::tempdir().unwrap();
    let destination = prepare(&GentleEngine::new(), &td.path().join("preview.json")).unwrap();
    assert!(
        destination
            .write(&InvalidSynthetic)
            .unwrap_err()
            .message
            .contains("sequence_design_report.serialize")
    );
    assert_eq!(fs::read_dir(td.path()).unwrap().count(), 0);
}

#[test]
fn lifecycle_slice4_admission_refusal_precedes_search_and_preserves_redo() {
    let td = tempfile::tempdir().unwrap();
    let project = td.path().join("project.json");
    let mut engine = fixture::live_with_redo();
    engine.set_host_project_path(Some(&project)).unwrap();
    let state = serde_json::to_value(engine.state()).unwrap();
    let revision = (engine.mutation_revision(), engine.structural_revision());
    let history = engine.history_summary();
    let mut callbacks = 0;
    let result = engine.apply_with_progress(
        fixture::plan(Some(project.display().to_string()), true),
        |_| {
            callbacks += 1;
            true
        },
    );
    assert!(result.unwrap_err().message.contains("protected_path"));
    assert_eq!(callbacks, 0);
    assert_eq!(serde_json::to_value(engine.state()).unwrap(), state);
    assert_eq!(
        (engine.mutation_revision(), engine.structural_revision()),
        revision
    );
    assert_eq!(
        engine.history_summary().operation_log_count,
        history.operation_log_count
    );
    assert_eq!(
        (engine.undo_available(), engine.redo_available()),
        (history.undo_count, history.redo_count)
    );
    assert!(!project.exists());
    engine.redo_last_operation().unwrap();
    assert!(engine.state().sequences.contains_key("synthetic_temporary"));
}

#[test]
fn lifecycle_slice4_cancelled_and_complete_report_bytes_match_results() {
    use gentle_protocol::sequence_design::DnaDesignStatus;
    for cancel in [false, true] {
        let td = tempfile::tempdir().unwrap();
        let report_path = td.path().join("report.json");
        let mut engine = fixture::live_with_redo();
        let before = serde_json::to_value(engine.state()).unwrap();
        let redo = engine.redo_available();
        let report = engine
            .apply_with_progress(
                fixture::plan(Some(report_path.display().to_string()), true),
                |_| !cancel,
            )
            .unwrap()
            .dna_sequence_design
            .unwrap();
        assert_eq!(
            report.status,
            if cancel {
                DnaDesignStatus::Cancelled
            } else {
                DnaDesignStatus::Feasible
            }
        );
        assert_eq!(
            fs::read(&report_path).unwrap(),
            serde_json::to_string_pretty(&report).unwrap().as_bytes()
        );
        assert_eq!(serde_json::to_value(engine.state()).unwrap(), before);
        assert_eq!(engine.redo_available(), redo);
    }
}

#[test]
fn lifecycle_slice4_engine_publish_race_preserves_competitor_state_and_journal() {
    let td = tempfile::tempdir().unwrap();
    let path = td.path().join("competitor.json");
    let mut engine = fixture::live_with_redo();
    let before = serde_json::to_value(engine.state()).unwrap();
    let revisions = (engine.mutation_revision(), engine.structural_revision());
    let history = engine.history_summary();
    let mut competed = false;
    let error = engine
        .apply_with_progress(
            fixture::plan(Some(path.display().to_string()), true),
            |_| {
                if !competed {
                    fs::write(&path, b"competing writer").unwrap();
                    competed = true;
                }
                true
            },
        )
        .unwrap_err();
    assert!(competed, "must reach search after successful preflight");
    assert!(error.message.contains("sequence_design_report.publish"));
    assert_eq!(fs::read(path).unwrap(), b"competing writer");
    assert_eq!(serde_json::to_value(engine.state()).unwrap(), before);
    assert_eq!(
        (engine.mutation_revision(), engine.structural_revision()),
        revisions
    );
    assert_eq!(
        (
            engine.history_summary().operation_log_count,
            engine.undo_available(),
            engine.redo_available()
        ),
        (
            history.operation_log_count,
            history.undo_count,
            history.redo_count
        )
    );
}

#[test]
fn lifecycle_slice4_workflow_duplicate_destination_keeps_published_prefix() {
    let td = tempfile::tempdir().unwrap();
    let first = td.path().join("first.json");
    let second = td.path().join("second.json");
    let mut fresh = GentleEngine::new();
    fresh
        .apply_workflow(Workflow {
            run_id: "synthetic-distinct".into(),
            ops: vec![
                fixture::plan(Some(first.display().to_string()), false),
                fixture::plan(Some(second.display().to_string()), false),
            ],
        })
        .unwrap();
    assert_eq!(fresh.history_summary().operation_log_count, 2);
    assert_eq!(fs::read(first).unwrap(), fs::read(second).unwrap());
    let path = td.path().join("prefix.json");
    let op = fixture::plan(Some(path.display().to_string()), false);
    let mut engine = GentleEngine::new();
    let error = engine
        .apply_workflow(Workflow {
            run_id: "synthetic-duplicate".into(),
            ops: vec![op.clone(), op],
        })
        .unwrap_err();
    assert!(error.message.contains("destination_exists"));
    assert_eq!(engine.history_summary().operation_log_count, 1);
    let report: serde_json::Value = serde_json::from_slice(&fs::read(path).unwrap()).unwrap();
    assert_eq!(report["output_sequence"], "ATGTTTAAGTAA");
}

#[test]
fn lifecycle_slice4_host_binding_is_not_serialized_or_a_project_edit() {
    let td = tempfile::tempdir().unwrap();
    let old = td.path().join("old.json");
    let new = td.path().join("new.json");
    let mut live = fixture::live_with_redo();
    let state = serde_json::to_value(live.state()).unwrap();
    let revision = (live.mutation_revision(), live.structural_revision());
    let history = serde_json::to_value(live.history_summary()).unwrap();
    live.set_host_project_path(Some(&old)).unwrap();
    assert_eq!(serde_json::to_value(live.state()).unwrap(), state);
    assert_eq!(
        (live.mutation_revision(), live.structural_revision()),
        revision
    );
    assert_eq!(
        serde_json::to_value(live.history_summary()).unwrap(),
        history
    );
    let mut encoded = serde_json::to_value(&live).unwrap();
    assert!(encoded.get("host_project_path").is_none());
    encoded["host_project_path"] = json!("untrusted-project.json");
    let decoded: GentleEngine = serde_json::from_value(encoded).unwrap();
    assert!(decoded.host_project_path().is_none());
    let mut checkpoint_import = live.fork_detached_execution();
    checkpoint_import.import_checkpoint_engine(decoded);
    assert_eq!(
        checkpoint_import.engine().host_project_path(),
        Some(old.as_path())
    );
    let mut detached = live.fork_detached_execution();
    assert_eq!(detached.engine().host_project_path(), Some(old.as_path()));
    detached.engine_mut().apply(fixture::snv()).unwrap();
    live.set_host_project_path(Some(&new)).unwrap();
    drop(live.commit_detached_execution(&mut detached).unwrap());
    assert_eq!(
        live.host_project_path(),
        Some(new.as_path()),
        "commit must not restore the snapshot's old filename"
    );
    let result: Result<(), &str> = live.with_rollback_on_error(|engine| {
        engine.set_host_project_path(None).unwrap();
        Err("synthetic rollback")
    });
    assert!(result.is_err());
    assert_eq!(live.host_project_path(), Some(new.as_path()));
    live.set_host_project_path(None).unwrap();
    assert!(live.host_project_path().is_none());
}

#[test]
fn lifecycle_slice4_shell_load_preserves_host_output_binding() {
    let td = tempfile::tempdir().unwrap();
    let source = td.path().join("source.json");
    fs::write(&source, serde_json::to_vec(&fixture::state()).unwrap()).unwrap();
    let output = td.path().join("host-output.json");
    let mut engine = bound(&output);
    crate::engine_shell::execute_shell_command(
        &mut engine,
        &crate::engine_shell::ShellCommand::LoadProject {
            path: source.display().to_string(),
        },
    )
    .unwrap();
    assert_eq!(engine.host_project_path(), Some(output.as_path()));
    assert!(prepare(&engine, &output).is_err());
}
