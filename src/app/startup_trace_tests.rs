//! Synthetic startup trace regressions through the real project-loading adapter.
//! Records are constructed in temporary directories; no user profile is changed.

use super::*;
use crate::gui_profiler::startup_trace;

#[test]
fn startup_trace_project_decode_and_install_keep_failure_and_success_separate() {
    let temp = tempfile::tempdir().unwrap();
    let missing = temp.path().join("private missing project.gentle");
    let path = temp.path().join("private project.gentle");
    ProjectState::default()
        .save_to_path(path.to_str().unwrap())
        .unwrap();
    let mut app = GENtleApp::default();
    let (_, report) = startup_trace::capture(|| {
        assert!(
            app.load_project_from_file_with_recent(missing.to_str().unwrap(), false)
                .is_err()
        );
        app.load_project_from_file_with_recent(path.to_str().unwrap(), false)
            .unwrap();
    });
    let events = report["events"].as_array().unwrap();
    for phase in ["project_load", "project_read_decode"] {
        assert_eq!(
            events
                .iter()
                .filter(|event| event["phase"] == phase && event["kind"] == "failed")
                .count(),
            1
        );
        assert_eq!(
            events
                .iter()
                .filter(|event| event["phase"] == phase && event["kind"] == "completed")
                .count(),
            1
        );
    }
    assert_eq!(
        events
            .iter()
            .filter(|event| event["phase"] == "project_install" && event["kind"] == "completed")
            .count(),
        1
    );
    assert_eq!(report["dropped_events"], 0);
    assert!(!report.to_string().contains("private"));
    assert!(!report.to_string().contains(path.to_str().unwrap()));
    assert_eq!(
        PathBuf::from(app.current_project_path.as_ref().unwrap()),
        path.canonicalize().unwrap()
    );
}

#[test]
fn startup_trace_open_dispatch_does_not_claim_a_missing_sequence_loaded() {
    let mut app = GENtleApp::default();
    let (_, report) = startup_trace::capture(|| app.open_sequence_window("missing"));
    assert!(app.new_windows.is_empty());
    let events = report["events"].as_array().unwrap();
    assert_eq!(events.len(), 2);
    assert!(
        events
            .iter()
            .all(|event| event["phase"] == "dna_open_dispatch")
    );
    assert_eq!(events[1]["kind"], "completed");
}
