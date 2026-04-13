//! DNA sequence-window wrapper and per-window controls.
//!
//! `WindowDna` is the viewport/lifecycle shell around
//! [`MainAreaDna`](crate::main_area_dna::MainAreaDna):
//! - it owns deferred sequence loading for lazy-open paths,
//! - it handles viewport-scoped controls such as help/main/configuration,
//! - it forwards ready-to-render sequence-window work into `MainAreaDna`.
//!
//! Keep adapter/window-lifecycle behavior here and route sequence-window
//! editing/inspection logic into `MainAreaDna` so the architectural boundary
//! stays clear for future edits.

use crate::{
    app::{request_open_graphics_settings_from_native_menu, request_open_help_from_native_menu},
    dna_sequence::DNAsequence,
    engine::GentleEngine,
    main_area_dna::MainAreaDna,
    window_backdrop::{
        WindowBackdropKind, current_window_backdrop_settings, paint_window_backdrop,
        with_window_content_inset,
    },
};
use eframe::egui;
use std::panic::{AssertUnwindSafe, catch_unwind};
use std::sync::{
    Arc, Mutex, RwLock,
    mpsc::{self, Receiver, TryRecvError},
};
use std::thread;

#[derive(Clone, Debug)]
enum DeferredAnalysisFocus {
    Dotplot(String),
    FlexibilityTrack(String),
    SequencingConfirmation(String),
    UniprotProjection(String),
    TranscriptProtein(Option<String>),
}

#[derive(Clone, Debug)]
/// Viewport shell around the main DNA-window controller.
///
/// The main rule is:
/// - `WindowDna` manages opening, deferred loading, and viewport-scoped chrome.
/// - `MainAreaDna` owns the actual sequence-window UI state and rendering.
pub struct WindowDna {
    main_area: MainAreaDna,
    pending_dna_load: Option<Arc<Mutex<Receiver<Result<DNAsequence, String>>>>>,
    deferred_load_message: Option<String>,
    deferred_analysis_focus: Option<DeferredAnalysisFocus>,
    close_requested: bool,
}

impl WindowDna {
    /// Construct an eager sequence window when sequence content is already in hand.
    pub fn new(dna: DNAsequence, seq_id: String, engine: Arc<RwLock<GentleEngine>>) -> Self {
        Self {
            main_area: MainAreaDna::new(dna, Some(seq_id), Some(engine)),
            pending_dna_load: None,
            deferred_load_message: None,
            deferred_analysis_focus: None,
            close_requested: false,
        }
    }

    /// Construct a lazy sequence window that resolves sequence payload from the
    /// shared engine in the background before handing off to `MainAreaDna`.
    pub fn new_lazy(seq_id: String, engine: Arc<RwLock<GentleEngine>>) -> Self {
        let (tx, rx) = mpsc::channel::<Result<DNAsequence, String>>();
        let thread_engine = engine.clone();
        let thread_seq_id = seq_id.clone();
        thread::spawn(move || {
            let result = thread_engine
                .read()
                .map_err(|_| "Engine lock poisoned during deferred load".to_string())
                .and_then(|guard| {
                    guard
                        .state()
                        .sequences
                        .get(&thread_seq_id)
                        .cloned()
                        .ok_or_else(|| {
                            format!(
                                "Sequence '{}' no longer exists while opening window",
                                thread_seq_id
                            )
                        })
                });
            let _ = tx.send(result);
        });
        let placeholder = DNAsequence::from_sequence("").expect("valid empty sequence");
        let mut main_area = MainAreaDna::new(placeholder, Some(seq_id), Some(engine));
        main_area.defer_feature_tree_until_interaction();
        Self {
            main_area,
            pending_dna_load: Some(Arc::new(Mutex::new(rx))),
            deferred_load_message: None,
            deferred_analysis_focus: None,
            close_requested: false,
        }
    }

    fn apply_deferred_analysis_focus(&mut self) {
        let Some(request) = self.deferred_analysis_focus.take() else {
            return;
        };
        match request {
            DeferredAnalysisFocus::Dotplot(dotplot_id) => {
                self.main_area.focus_dotplot_analysis(&dotplot_id);
            }
            DeferredAnalysisFocus::FlexibilityTrack(track_id) => {
                self.main_area.focus_flexibility_track_analysis(&track_id);
            }
            DeferredAnalysisFocus::SequencingConfirmation(report_id) => {
                self.main_area
                    .focus_sequencing_confirmation_report(&report_id);
            }
            DeferredAnalysisFocus::UniprotProjection(projection_id) => {
                self.main_area
                    .focus_uniprot_projection_expert(&projection_id);
            }
            DeferredAnalysisFocus::TranscriptProtein(transcript_id_filter) => {
                self.main_area
                    .focus_transcript_protein_expert(transcript_id_filter.as_deref());
            }
        }
    }

    fn poll_deferred_load(&mut self) {
        if let Some(rx) = self.pending_dna_load.as_ref() {
            let recv_result = rx
                .lock()
                .map(|guard| guard.try_recv())
                .map_err(|_| "Deferred sequence load lock poisoned".to_string());
            match recv_result {
                Ok(Ok(Ok(dna))) => {
                    self.main_area.replace_loaded_sequence(dna);
                    self.pending_dna_load = None;
                    self.deferred_load_message = None;
                    self.apply_deferred_analysis_focus();
                }
                Ok(Ok(Err(message))) => {
                    self.pending_dna_load = None;
                    self.deferred_load_message = Some(message);
                }
                Ok(Err(TryRecvError::Empty)) => {}
                Ok(Err(TryRecvError::Disconnected)) => {
                    self.pending_dna_load = None;
                    self.deferred_load_message = Some(
                        "Deferred sequence load channel disconnected unexpectedly".to_string(),
                    );
                }
                Err(message) => {
                    self.pending_dna_load = None;
                    self.deferred_load_message = Some(message);
                }
            }
        }
    }

    fn render_nav_row(&mut self, ui: &mut egui::Ui) {
        ui.horizontal(|ui| {
            if ui
                .button("Help")
                .on_hover_text("Open GUI help (F1 on Windows/Linux, Cmd+Shift+/ on macOS)")
                .clicked()
            {
                request_open_help_from_native_menu();
            }
            if ui
                .button("Main")
                .on_hover_text("Bring the main project workspace to front")
                .clicked()
            {
                ui.ctx().send_viewport_cmd_to(
                    egui::ViewportId::ROOT,
                    egui::ViewportCommand::Visible(true),
                );
                ui.ctx()
                    .send_viewport_cmd_to(egui::ViewportId::ROOT, egui::ViewportCommand::Focus);
            }
            if ui
                .button("Configuration")
                .on_hover_text("Open Configuration window on Graphics settings")
                .clicked()
            {
                request_open_graphics_settings_from_native_menu();
            }
            if ui
                .button("Close")
                .on_hover_text("Close this sequence window (Cmd/Ctrl+W)")
                .clicked()
            {
                self.close_requested = true;
            }
        });
    }

    /// Drive one viewport update.
    ///
    /// This polls deferred load state, renders viewport-scoped controls, and
    /// then delegates the actual sequence-window layout to `MainAreaDna`.
    pub fn update(&mut self, ctx: &egui::Context) {
        self.poll_deferred_load();
        let result = catch_unwind(AssertUnwindSafe(|| {
            let open_help_f1 = egui::KeyboardShortcut::new(egui::Modifiers::NONE, egui::Key::F1);
            let open_help_ctrl_f1 =
                egui::KeyboardShortcut::new(egui::Modifiers::CTRL, egui::Key::F1);
            let open_help_cmd_shift_slash = egui::KeyboardShortcut::new(
                egui::Modifiers::COMMAND | egui::Modifiers::SHIFT,
                egui::Key::Slash,
            );
            if ctx.input_mut(|i| i.consume_shortcut(&open_help_f1))
                || ctx.input_mut(|i| i.consume_shortcut(&open_help_ctrl_f1))
                || ctx.input_mut(|i| i.consume_shortcut(&open_help_cmd_shift_slash))
            {
                request_open_help_from_native_menu();
            }
            let kind = if self.main_area.opened_from_pool_context() {
                WindowBackdropKind::Pool
            } else {
                WindowBackdropKind::Sequence
            };
            let settings = current_window_backdrop_settings();
            let nav_panel_id = egui::Id::new(("window_dna_nav", self.main_area.sequence_id()));
            crate::egui_compat::show_top_panel(
                ctx,
                nav_panel_id,
                egui::Panel::top(nav_panel_id).frame(egui::Frame::NONE),
                |ui| {
                    paint_window_backdrop(ui, kind, &settings);
                    self.render_nav_row(ui);
                },
            );
            if self.pending_dna_load.is_some() {
                crate::egui_compat::show_central_panel(
                    ctx,
                    egui::CentralPanel::default().frame(egui::Frame::NONE),
                    |ui| {
                        paint_window_backdrop(ui, kind, &settings);
                        with_window_content_inset(ui, |ui| {
                            ui.vertical_centered(|ui| {
                                ui.add_space(48.0);
                                ui.add(egui::Spinner::new());
                                ui.add_space(6.0);
                                ui.label(
                                    "Opening sequence window: loading sequence content in background...",
                                );
                            });
                        });
                    },
                );
                ctx.request_repaint();
            } else if let Some(message) = self.deferred_load_message.as_deref() {
                crate::egui_compat::show_central_panel(
                    ctx,
                    egui::CentralPanel::default().frame(egui::Frame::NONE),
                    |ui| {
                        paint_window_backdrop(ui, kind, &settings);
                        with_window_content_inset(ui, |ui| {
                            ui.colored_label(egui::Color32::from_rgb(180, 32, 32), message);
                        });
                    },
                );
            } else {
                // MainAreaDna owns the root panel layout for sequence windows.
                self.main_area.render(ctx);
            }
        }));
        if result.is_err() {
            eprintln!("E WindowDna: recovered from panic while rendering DNA window");
        }
    }

    pub fn update_embedded(&mut self, ui: &mut egui::Ui) {
        self.poll_deferred_load();
        let result = catch_unwind(AssertUnwindSafe(|| {
            let open_help_f1 = egui::KeyboardShortcut::new(egui::Modifiers::NONE, egui::Key::F1);
            let open_help_ctrl_f1 =
                egui::KeyboardShortcut::new(egui::Modifiers::CTRL, egui::Key::F1);
            let open_help_cmd_shift_slash = egui::KeyboardShortcut::new(
                egui::Modifiers::COMMAND | egui::Modifiers::SHIFT,
                egui::Key::Slash,
            );
            if ui.ctx().input_mut(|i| i.consume_shortcut(&open_help_f1))
                || ui
                    .ctx()
                    .input_mut(|i| i.consume_shortcut(&open_help_ctrl_f1))
                || ui
                    .ctx()
                    .input_mut(|i| i.consume_shortcut(&open_help_cmd_shift_slash))
            {
                request_open_help_from_native_menu();
            }
            let kind = if self.main_area.opened_from_pool_context() {
                WindowBackdropKind::Pool
            } else {
                WindowBackdropKind::Sequence
            };
            let settings = current_window_backdrop_settings();
            let nav_panel_id =
                egui::Id::new(("window_dna_nav_embedded", self.main_area.panel_scope_key()));
            crate::egui_compat::show_top_panel_inside(
                ui,
                egui::Panel::top(nav_panel_id).frame(egui::Frame::NONE),
                |ui| {
                    paint_window_backdrop(ui, kind, &settings);
                    self.render_nav_row(ui);
                },
            );
            if self.pending_dna_load.is_some() {
                crate::egui_compat::show_central_panel_inside(
                    ui,
                    egui::CentralPanel::default().frame(egui::Frame::NONE),
                    |ui| {
                        paint_window_backdrop(ui, kind, &settings);
                        with_window_content_inset(ui, |ui| {
                            ui.vertical_centered(|ui| {
                                ui.add_space(48.0);
                                ui.add(egui::Spinner::new());
                                ui.add_space(6.0);
                                ui.label(
                                    "Opening sequence window: loading sequence content in background...",
                                );
                            });
                        });
                    },
                );
                ui.ctx().request_repaint();
            } else if let Some(message) = self.deferred_load_message.as_deref() {
                crate::egui_compat::show_central_panel_inside(
                    ui,
                    egui::CentralPanel::default().frame(egui::Frame::NONE),
                    |ui| {
                        paint_window_backdrop(ui, kind, &settings);
                        with_window_content_inset(ui, |ui| {
                            ui.colored_label(egui::Color32::from_rgb(180, 32, 32), message);
                        });
                    },
                );
            } else {
                self.main_area.render_inside(ui);
            }
        }));
        if result.is_err() {
            eprintln!("E WindowDna: recovered from panic while rendering embedded DNA window");
        }
    }

    pub fn name(&self) -> String {
        self.main_area.window_title()
    }

    pub fn sequence_id(&self) -> Option<String> {
        self.main_area.sequence_id().map(|v| v.to_string())
    }

    pub fn take_close_requested(&mut self) -> bool {
        let requested = self.close_requested;
        self.close_requested = false;
        requested
    }

    pub fn set_window_scope_id(&mut self, scope_id: String) {
        self.main_area.set_window_scope_id(scope_id);
    }

    pub fn selection_range_0based(&self) -> Option<(usize, usize)> {
        self.main_area.selection_range_0based()
    }

    pub fn collect_open_auxiliary_window_entries(&self) -> Vec<(egui::ViewportId, String, String)> {
        self.main_area.collect_open_auxiliary_window_entries()
    }

    pub fn embedded_auxiliary_window_layer_id(
        &self,
        viewport_id: egui::ViewportId,
    ) -> Option<egui::LayerId> {
        self.main_area
            .embedded_auxiliary_window_layer_id(viewport_id)
    }

    pub fn render_pcr_designer_specialist(&mut self, ui: &mut egui::Ui, ctx: &egui::Context) {
        self.main_area.render_pcr_designer_specialist(ui, ctx);
    }

    pub fn render_sequencing_confirmation_specialist(
        &mut self,
        ui: &mut egui::Ui,
        ctx: &egui::Context,
    ) {
        self.main_area
            .render_sequencing_confirmation_specialist(ui, ctx);
    }

    pub fn set_pool_context(&mut self, pool_seq_ids: Vec<String>) {
        self.main_area.set_pool_context(pool_seq_ids);
    }

    pub fn enable_compact_lane_layout(&mut self) {
        self.main_area.enable_compact_lane_layout();
    }

    pub fn focus_dotplot_analysis(&mut self, dotplot_id: &str) {
        if self.pending_dna_load.is_some() {
            self.deferred_analysis_focus =
                Some(DeferredAnalysisFocus::Dotplot(dotplot_id.to_string()));
            return;
        }
        self.main_area.focus_dotplot_analysis(dotplot_id);
    }

    pub fn focus_flexibility_track_analysis(&mut self, track_id: &str) {
        if self.pending_dna_load.is_some() {
            self.deferred_analysis_focus = Some(DeferredAnalysisFocus::FlexibilityTrack(
                track_id.to_string(),
            ));
            return;
        }
        self.main_area.focus_flexibility_track_analysis(track_id);
    }

    pub fn focus_sequencing_confirmation_report(&mut self, report_id: &str) {
        if self.pending_dna_load.is_some() {
            self.deferred_analysis_focus = Some(DeferredAnalysisFocus::SequencingConfirmation(
                report_id.to_string(),
            ));
            return;
        }
        self.main_area
            .focus_sequencing_confirmation_report(report_id);
    }

    pub fn focus_uniprot_projection_expert(&mut self, projection_id: &str) {
        if self.pending_dna_load.is_some() {
            self.deferred_analysis_focus = Some(DeferredAnalysisFocus::UniprotProjection(
                projection_id.to_string(),
            ));
            return;
        }
        self.main_area
            .focus_uniprot_projection_expert(projection_id);
    }

    pub fn focus_transcript_protein_expert(&mut self, transcript_id_filter: Option<&str>) {
        if self.pending_dna_load.is_some() {
            self.deferred_analysis_focus = Some(DeferredAnalysisFocus::TranscriptProtein(
                transcript_id_filter
                    .map(str::trim)
                    .filter(|value| !value.is_empty())
                    .map(|value| value.to_string()),
            ));
            return;
        }
        self.main_area
            .focus_transcript_protein_expert(transcript_id_filter);
    }

    pub fn refresh_from_engine_settings(&mut self) {
        self.main_area.refresh_from_engine_settings();
    }

    pub fn refresh_from_engine_state(&mut self) {
        self.main_area.refresh_from_engine_sequence_state();
        self.main_area.refresh_from_engine_settings();
    }
}
