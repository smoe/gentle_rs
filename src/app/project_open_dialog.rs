//! Nonblocking project-file selection with project ownership and dirty guards.
//!
//! The native picker is polled, never synchronously awaited, on the GUI thread.
//! In particular, macOS must not enter `runModal` from a winit event callback.

use super::{GENtleApp, ProjectAction};
use eframe::egui;
use std::{
    future::Future,
    path::PathBuf,
    pin::Pin,
    sync::Arc,
    task::{Context, Poll, Wake, Waker},
};

type ProjectSelection = Pin<Box<dyn Future<Output = Option<PathBuf>> + Send>>;

pub(super) struct PendingProjectOpenDialog {
    owner_instance_id: u64,
    approved_mutation_revision: u64,
    selection: ProjectSelection,
}

struct RepaintProjectDialog(egui::Context);

impl Wake for RepaintProjectDialog {
    fn wake(self: Arc<Self>) {
        self.wake_by_ref();
    }

    fn wake_by_ref(self: &Arc<Self>) {
        self.0.request_repaint_of(egui::ViewportId::ROOT);
    }
}

impl GENtleApp {
    pub(super) fn prompt_open_project(&mut self) {
        self.begin_project_open_dialog_with(|| {
            Box::pin(async {
                rfd::AsyncFileDialog::new()
                    .add_filter("GENtle project", &["json"])
                    .pick_file()
                    .await
                    .map(|file| file.path().to_path_buf())
            })
        });
    }

    fn begin_project_open_dialog_with(&mut self, select: impl FnOnce() -> ProjectSelection) {
        if self.pending_project_open_dialog.is_some() {
            self.app_status = "A project file dialog is already open".into();
            return;
        }
        let (owner_instance_id, approved_mutation_revision) = {
            let engine = self.engine.read().unwrap();
            (engine.instance_id(), engine.mutation_revision())
        };
        self.pending_project_open_dialog = Some(PendingProjectOpenDialog {
            owner_instance_id,
            approved_mutation_revision,
            selection: select(),
        });
        self.app_status = "Choose a project file; Cancel keeps the current project".into();
    }

    pub(super) fn poll_project_open_dialog(&mut self, ctx: &egui::Context) {
        let Some(pending) = self.pending_project_open_dialog.as_mut() else {
            return;
        };
        let waker = Waker::from(Arc::new(RepaintProjectDialog(ctx.clone())));
        let mut context = Context::from_waker(&waker);
        let Poll::Ready(path) = pending.selection.as_mut().poll(&mut context) else {
            return;
        };
        let owner_instance_id = pending.owner_instance_id;
        let approved_mutation_revision = pending.approved_mutation_revision;
        self.pending_project_open_dialog = None;
        ctx.request_repaint_of(egui::ViewportId::ROOT);

        let (current_owner, current_revision) = {
            let engine = self.engine.read().unwrap();
            (engine.instance_id(), engine.mutation_revision())
        };
        if current_owner != owner_instance_id {
            self.app_status = "Ignored project file selection after project replacement".into();
            return;
        }
        if self.pending_project_action.is_some() || self.pending_app_quit {
            self.app_status =
                "Ignored project file selection while another project action is pending".into();
            return;
        }
        if let Some(path) = path {
            let action = ProjectAction::OpenPath(path.display().to_string());
            if current_revision == approved_mutation_revision {
                self.execute_project_action(action);
            } else {
                // Initial discard approval does not cover edits made during selection.
                self.request_project_action(action);
            }
        } else {
            self.app_status = "Project opening cancelled".into();
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{
        egui_compat,
        engine::{Engine, Operation, ProjectState},
    };
    use std::{
        fs,
        sync::{
            Mutex,
            atomic::{AtomicUsize, Ordering},
        },
    };

    // All project bytes are generated from hand-crafted, synthetic state.
    fn fixture(path: &std::path::Path, crlf: bool) -> Vec<u8> {
        let mut state = ProjectState::default();
        state.metadata.insert("picked_fixture".into(), true.into());
        let text = serde_json::to_string_pretty(&state).unwrap();
        let bytes = if crlf {
            text.replace('\n', "\r\n")
        } else {
            text
        }
        .into_bytes();
        fs::write(path, &bytes).unwrap();
        bytes
    }

    fn select(app: &mut GENtleApp, path: Option<PathBuf>) {
        app.begin_project_open_dialog_with(|| Box::pin(std::future::ready(path)));
    }

    fn add_unsaved_sequence(app: &mut GENtleApp) {
        app.engine
            .write()
            .unwrap()
            .apply(Operation::CreateSequenceFromText {
                sequence_text: "ACGT".into(),
                output_id: Some("unsaved".into()),
                name: None,
                circular: false,
            })
            .unwrap();
    }

    #[test]
    fn project_open_dialog_request_defers_native_picker_and_preserves_state() {
        let mut app = GENtleApp::default();
        let before = serde_json::to_value(app.engine.read().unwrap().state()).unwrap();
        app.request_project_action(ProjectAction::Open);
        assert!(app.pending_project_open_dialog.is_some());
        assert!(app.pending_project_action.is_none());
        assert_eq!(
            serde_json::to_value(app.engine.read().unwrap().state()).unwrap(),
            before
        );
        // Drop the unpolled selection: this test must never open an OS dialog.
        app.pending_project_open_dialog = None;
    }

    #[test]
    fn project_open_dialog_pending_future_is_polled_once_and_can_wake_root() {
        let mut app = GENtleApp::default();
        let ctx = egui::Context::default();
        let polls = Arc::new(AtomicUsize::new(0));
        let observed_polls = polls.clone();
        let wake = Arc::new(Mutex::new(None::<Waker>));
        let observed_wake = wake.clone();
        let repaints = Arc::new(Mutex::new(Vec::new()));
        let observed_repaints = repaints.clone();
        ctx.set_repaint_observer(move |info| {
            observed_repaints.lock().unwrap().push(info.viewport_id);
        });
        app.begin_project_open_dialog_with(|| {
            Box::pin(std::future::poll_fn(move |cx| {
                observed_polls.fetch_add(1, Ordering::Relaxed);
                *observed_wake.lock().unwrap() = Some(cx.waker().clone());
                Poll::<Option<PathBuf>>::Pending
            }))
        });
        app.poll_project_open_dialog(&ctx);
        assert_eq!(polls.load(Ordering::Relaxed), 1);
        assert!(app.pending_project_open_dialog.is_some());
        assert!(repaints.lock().unwrap().is_empty());
        wake.lock().unwrap().take().unwrap().wake();
        assert_eq!(*repaints.lock().unwrap(), [egui::ViewportId::ROOT]);
    }

    #[test]
    fn project_open_dialog_duplicate_request_does_not_create_second_picker() {
        let mut app = GENtleApp::default();
        select(&mut app, None);
        app.begin_project_open_dialog_with(|| panic!("second picker must not be constructed"));
        assert!(app.pending_project_open_dialog.is_some());
        app.poll_project_open_dialog(&egui::Context::default());
        assert!(app.pending_project_open_dialog.is_none());
        assert_eq!(app.app_status, "Project opening cancelled");
    }

    #[test]
    fn project_open_dialog_cancel_preserves_project_path_dirty_state_and_history() {
        let mut app = GENtleApp::default();
        add_unsaved_sequence(&mut app);
        {
            let mut engine = app.engine.write().unwrap();
            engine
                .apply(Operation::CreateSequenceFromText {
                    sequence_text: "TGCA".into(),
                    output_id: Some("undo_control".into()),
                    name: None,
                    circular: false,
                })
                .unwrap();
            engine.undo_last_operation().unwrap();
        }
        let owner = app.engine.read().unwrap().instance_id();
        let before = serde_json::to_value(app.engine.read().unwrap().state()).unwrap();
        app.current_project_path = Some("original.gentle.json".into());
        let path_before = app.current_project_path.clone();
        let undo_before = app.engine.read().unwrap().undo_available();
        let redo_before = app.engine.read().unwrap().redo_available();
        assert_eq!(redo_before, 1);
        select(&mut app, None);
        app.poll_project_open_dialog(&egui::Context::default());
        assert_eq!(app.engine.read().unwrap().instance_id(), owner);
        assert_eq!(
            serde_json::to_value(app.engine.read().unwrap().state()).unwrap(),
            before
        );
        assert_eq!(app.engine.read().unwrap().undo_available(), undo_before);
        assert_eq!(app.engine.read().unwrap().redo_available(), redo_before);
        assert_eq!(app.current_project_path, path_before);
        assert!(app.is_project_dirty());
        assert!(app.pending_project_action.is_none());
    }

    #[test]
    fn project_open_dialog_selected_project_loads_once_with_lf_crlf_and_unicode_path() {
        for crlf in [false, true] {
            let temp = tempfile::tempdir().unwrap();
            let path = temp.path().join("Owner's Pr\u{00fc}fung.gentle.json");
            let bytes = fixture(&path, crlf);
            let mut app = GENtleApp::default();
            // Avoid writing user preferences; the real recent-path code is unchanged.
            app.recent_project_paths =
                vec![GENtleApp::normalize_project_path(path.to_str().unwrap())];
            let owner = app.engine.read().unwrap().instance_id();
            select(&mut app, Some(path.clone()));
            app.poll_project_open_dialog(&egui::Context::default());
            let loaded_owner = app.engine.read().unwrap().instance_id();
            assert_ne!(loaded_owner, owner);
            assert_eq!(app.current_project_path.as_deref(), path.to_str());
            assert_eq!(
                app.engine.read().unwrap().state().metadata["picked_fixture"],
                true
            );
            assert!(!app.is_project_dirty());
            assert!(app.pending_project_open_dialog.is_none());
            app.poll_project_open_dialog(&egui::Context::default());
            assert_eq!(app.engine.read().unwrap().instance_id(), loaded_owner);
            assert_eq!(fs::read(&path).unwrap(), bytes);
        }
    }

    #[test]
    fn project_open_dialog_late_edit_requires_review_before_selected_file_io() {
        let temp = tempfile::tempdir().unwrap();
        let path = temp.path().join("not-yet-created.json");
        let mut app = GENtleApp::default();
        select(&mut app, Some(path.clone()));
        add_unsaved_sequence(&mut app);
        app.poll_project_open_dialog(&egui::Context::default());
        assert!(matches!(
            app.pending_project_action.as_ref(),
            Some(ProjectAction::OpenPath(selected)) if selected == path.to_str().unwrap()
        ));
        assert!(
            app.engine
                .read()
                .unwrap()
                .state()
                .sequences
                .contains_key("unsaved")
        );
        assert!(app.current_project_path.is_none());
        assert!(!path.exists());
        assert!(app.is_project_dirty());
    }

    #[test]
    fn project_open_dialog_dirty_initial_project_is_reviewed_before_picker() {
        let mut app = GENtleApp::default();
        add_unsaved_sequence(&mut app);
        app.request_project_action(ProjectAction::Open);
        assert!(matches!(
            app.pending_project_action,
            Some(ProjectAction::Open)
        ));
        assert!(app.pending_project_open_dialog.is_none());
    }

    #[test]
    fn project_open_dialog_unchanged_discard_approval_is_not_requested_twice() {
        let temp = tempfile::tempdir().unwrap();
        let path = temp.path().join("selected.gentle.json");
        let bytes = fixture(&path, false);
        let mut app = GENtleApp::default();
        app.recent_project_paths = vec![GENtleApp::normalize_project_path(path.to_str().unwrap())];
        add_unsaved_sequence(&mut app);
        app.request_project_action(ProjectAction::Open);
        let approved_action = app.pending_project_action.take().unwrap();
        assert!(matches!(approved_action, ProjectAction::Open));
        // This is the existing Unsaved Changes -> Discard action, not a new bypass.
        app.execute_project_action(approved_action);
        app.pending_project_open_dialog.as_mut().unwrap().selection =
            Box::pin(std::future::ready(Some(path.clone())));
        app.poll_project_open_dialog(&egui::Context::default());
        assert!(app.pending_project_action.is_none());
        assert_eq!(app.current_project_path.as_deref(), path.to_str());
        assert!(
            !app.engine
                .read()
                .unwrap()
                .state()
                .sequences
                .contains_key("unsaved")
        );
        assert_eq!(
            app.engine.read().unwrap().state().metadata["picked_fixture"],
            true
        );
        assert_eq!(fs::read(path).unwrap(), bytes);
    }

    #[test]
    fn project_open_dialog_project_replacement_discards_old_selection() {
        let mut app = GENtleApp::default();
        select(&mut app, Some(PathBuf::from("must-not-be-read.json")));
        app.new_project();
        let owner = app.engine.read().unwrap().instance_id();
        app.poll_project_open_dialog(&egui::Context::default());
        assert_eq!(app.engine.read().unwrap().instance_id(), owner);
        assert!(app.current_project_path.is_none());
        assert!(app.pending_project_action.is_none());
        assert!(app.pending_project_open_dialog.is_none());
        assert_eq!(
            app.app_status,
            "Ignored project file selection after project replacement"
        );
    }

    #[test]
    fn project_open_dialog_result_cannot_replace_pending_confirmation_or_quit() {
        for quit in [false, true] {
            let mut app = GENtleApp::default();
            select(&mut app, Some(PathBuf::from("must-not-be-read.json")));
            if quit {
                app.pending_app_quit = true;
            } else {
                app.pending_project_action = Some(ProjectAction::New);
            }
            app.poll_project_open_dialog(&egui::Context::default());
            assert!(app.current_project_path.is_none());
            assert_eq!(app.pending_app_quit, quit);
            if !quit {
                assert!(matches!(
                    app.pending_project_action,
                    Some(ProjectAction::New)
                ));
            }
            assert!(app.pending_project_open_dialog.is_none());
        }
    }

    #[test]
    fn project_open_dialog_polls_selection_during_gui_pass_without_native_dialog() {
        let ctx = egui::Context::default();
        let mut app = GENtleApp::default();
        select(&mut app, None);
        ctx.begin_pass(egui::RawInput::default());
        app.poll_project_open_dialog(&ctx);
        egui_compat::discard_test_pass_output(&ctx);
        assert!(app.pending_project_open_dialog.is_none());
        assert_eq!(app.app_status, "Project opening cancelled");
    }
}
