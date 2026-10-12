//! Startup trace regressions through the real help and project-loading adapters.
//! Synthetic records/cache sentinels use temporary paths; help uses the public
//! repository catalog. No user profile is changed and no timing target is asserted.

use super::*;
use crate::gui_profiler::startup_trace;

#[test]
fn startup_trace_help_image_adapter_keeps_cache_and_failure_behavior() {
    let temp = tempfile::tempdir().unwrap();
    let svg = temp.path().join("private missing.svg");
    let png = GENtleApp::help_svg_png_cache_path(&svg).unwrap();
    let ordinary = temp.path().join("ordinary.png");
    let (_, failed) = startup_trace::capture(|| {
        assert_eq!(GENtleApp::help_image_render_path(&ordinary), ordinary);
        for _ in 0..2 {
            assert_eq!(GENtleApp::help_image_render_path(&svg), svg);
        }
    });
    assert_eq!(failed["help_image_work"]["svg_references"], 2);
    assert_eq!(failed["help_image_work"]["rasterization_attempts"], 2);
    assert_eq!(failed["help_image_work"]["rasterization_failures"], 2);
    assert!(!png.exists());

    // The production hit rule is is_file(), not image decoding. Only this
    // temporary source's unique cache entry is created and removed by the test.
    fs::create_dir_all(png.parent().unwrap()).unwrap();
    fs::write(&png, b"synthetic cache-hit sentinel").unwrap();
    let (result, cached) = startup_trace::capture(|| GENtleApp::help_image_render_path(&svg));
    fs::remove_file(&png).unwrap();
    assert_eq!(result, png);
    assert_eq!(cached["help_image_work"]["svg_references"], 1);
    assert_eq!(cached["help_image_work"]["cache_hits"], 1);
    assert_eq!(cached["help_image_work"]["rasterization_attempts"], 0);
    for report in [failed, cached] {
        assert_eq!(report["dropped_events"], 0);
        assert!(!report.to_string().contains("private"));
        assert!(!report.to_string().contains(svg.to_str().unwrap()));
    }
}

#[test]
fn startup_trace_help_refresh_and_first_use_preserve_contents() {
    let mut app = GENtleApp::default();
    let (_, report) = startup_trace::capture(|| {
        startup_trace::measure(startup_trace::Phase::HelpPreparation, || {
            app.refresh_help_docs()
        });
        let manuals = (
            app.help_gui_markdown.clone(),
            app.help_cli_markdown.clone(),
            app.help_agent_interface_markdown.clone(),
            app.help_reviewer_preview_markdown.clone(),
            app.help_shell_markdown.clone(),
            app.help_tutorial_markdown.clone(),
        );
        for _ in 0..2 {
            app.show_help_dialog = false;
            app.open_help_doc(HelpDoc::Gui);
            app.open_help_tutorial_doc(app.help_tutorial_selected);
        }
        assert_eq!(
            manuals,
            (
                app.help_gui_markdown.clone(),
                app.help_cli_markdown.clone(),
                app.help_agent_interface_markdown.clone(),
                app.help_reviewer_preview_markdown.clone(),
                app.help_shell_markdown.clone(),
                app.help_tutorial_markdown.clone(),
            )
        );
        let first = GENtleApp::tutorial_project_guided_walkthrough_entries();
        let second = GENtleApp::tutorial_project_guided_walkthrough_entries();
        let identity = |rows: &[HelpTutorialDocEntry]| {
            rows.iter()
                .map(|row| {
                    (
                        row.title.clone(),
                        row.path.clone(),
                        row.summary.clone(),
                        row.group_label.clone(),
                        row.group_order,
                        row.group_position,
                    )
                })
                .collect::<Vec<_>>()
        };
        assert_eq!(identity(&first), identity(&second));
        assert!(app.help_tutorial_entries.len() > 1);
        let initial = app.help_tutorial_selected;
        let next = (initial + 1) % app.help_tutorial_entries.len();
        for index in [next, initial] {
            assert!(app.set_help_tutorial_selected(index));
            assert_eq!(app.help_tutorial_selected, index);
        }
    });
    let events = report["events"].as_array().unwrap();
    let completed: Vec<_> = events
        .iter()
        .filter(|row| row["kind"] == "completed")
        .map(|row| row["phase"].as_str().unwrap())
        .collect();
    assert_eq!(
        completed,
        [
            "help_manuals",
            "help_shell_reference",
            "help_tutorial_discovery",
            "help_tutorial_selected_load",
            "help_preparation",
            "help_open",
            "help_tutorial_open",
            "help_tutorial_menu_discovery",
            "help_tutorial_switch",
        ]
    );
    assert_eq!(events.first().unwrap()["phase"], "help_preparation");
    assert_eq!(report["dropped_events"], 0);
}

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
