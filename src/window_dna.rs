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
use std::time::Duration;

const DEFERRED_LOAD_REPAINT_INTERVAL: Duration = Duration::from_millis(100);

#[derive(Clone, Debug)]
enum DeferredAnalysisFocus {
    Dotplot(String),
    FlexibilityTrack(String),
    CrypticSplicing,
    RnaReadReport(String),
    PrimerDesign(String),
    PrimerSpecificity(String),
    QpcrDesign(String),
    RestrictionCloningPcrHandoff(String),
    SequencingConfirmation(String),
    ConstructReasoning(String),
    UniprotProjection {
        projection_id: String,
        protein_feature_filter: gentle_protocol::ProteinFeatureFilter,
    },
    TranscriptProtein(Option<String>),
    EnsemblProtein {
        transcript_id_filter: Option<String>,
        entry_id: String,
        protein_feature_filter: gentle_protocol::ProteinFeatureFilter,
    },
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
    awaiting_first_content_paint: bool,
    deferred_load_message: Option<String>,
    deferred_analysis_focus: Option<DeferredAnalysisFocus>,
    close_requested: bool,
}

impl WindowDna {
    fn render_deferred_load_indicator(ui: &mut egui::Ui) {
        let phase = ((ui.input(|input| input.time) * 10.0) as usize) % 4;
        let marker = match phase {
            0 => "[   ]",
            1 => "[.  ]",
            2 => "[.. ]",
            _ => "[...]",
        };
        ui.monospace(marker);
        ui.add_space(6.0);
        ui.label(crate::i18n::tr("sequence.loading"));
        ui.ctx()
            .request_repaint_after(DEFERRED_LOAD_REPAINT_INTERVAL);
    }

    /// Construct an eager sequence window when sequence content is already in hand.
    pub fn new(dna: DNAsequence, seq_id: String, engine: Arc<RwLock<GentleEngine>>) -> Self {
        Self {
            main_area: MainAreaDna::new(dna, Some(seq_id), Some(engine)),
            pending_dna_load: None,
            awaiting_first_content_paint: true,
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
        crate::gentle_gui_profile_scope!("WindowDna::new_lazy.placeholder");
        let mut main_area = MainAreaDna::new(placeholder, Some(seq_id), Some(engine));
        main_area.defer_feature_tree_until_interaction();
        Self {
            main_area,
            pending_dna_load: Some(Arc::new(Mutex::new(rx))),
            awaiting_first_content_paint: true,
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
            DeferredAnalysisFocus::CrypticSplicing => {
                self.main_area.focus_cryptic_splicing_screen();
            }
            DeferredAnalysisFocus::RnaReadReport(report_id) => {
                self.main_area.focus_rna_read_report(&report_id);
            }
            DeferredAnalysisFocus::PrimerDesign(report_id) => {
                self.main_area.focus_primer_design_report(&report_id);
            }
            DeferredAnalysisFocus::PrimerSpecificity(report_id) => {
                self.main_area.focus_primer_specificity_report(&report_id);
            }
            DeferredAnalysisFocus::QpcrDesign(report_id) => {
                self.main_area.focus_qpcr_design_report(&report_id);
            }
            DeferredAnalysisFocus::RestrictionCloningPcrHandoff(report_id) => {
                self.main_area
                    .focus_restriction_cloning_handoff_report(&report_id);
            }
            DeferredAnalysisFocus::SequencingConfirmation(report_id) => {
                self.main_area
                    .focus_sequencing_confirmation_report(&report_id);
            }
            DeferredAnalysisFocus::ConstructReasoning(graph_id) => {
                self.main_area.focus_construct_reasoning_graph(&graph_id);
            }
            DeferredAnalysisFocus::UniprotProjection {
                projection_id,
                protein_feature_filter,
            } => {
                self.main_area
                    .focus_uniprot_projection_expert(&projection_id, protein_feature_filter);
            }
            DeferredAnalysisFocus::TranscriptProtein(transcript_id_filter) => {
                self.main_area
                    .focus_transcript_protein_expert(transcript_id_filter.as_deref());
            }
            DeferredAnalysisFocus::EnsemblProtein {
                transcript_id_filter,
                entry_id,
                protein_feature_filter,
            } => {
                self.main_area.focus_ensembl_entry_protein_expert(
                    transcript_id_filter.as_deref(),
                    &entry_id,
                    protein_feature_filter,
                );
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
                    crate::gentle_gui_profile_scope!("WindowDna::poll_deferred_load.hydrate");
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
        let agent_help_window_title = self.name();
        ui.horizontal(|ui| {
            if ui
                .button(crate::i18n::tr("button.help"))
                .on_hover_text("Open GUI help (F1 on Windows/Linux, Cmd+Shift+/ on macOS)")
                .clicked()
            {
                request_open_help_from_native_menu();
            }
            if ui
                .button(crate::i18n::tr("button.main"))
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
                .button(crate::i18n::tr("button.configuration"))
                .on_hover_text("Open Configuration window on Graphics settings")
                .clicked()
            {
                request_open_graphics_settings_from_native_menu();
            }
            crate::agent_help::render_agent_help_button(
                ui,
                agent_help_window_title.clone(),
                "window.dna_viewer",
            );
            if ui
                .button(crate::i18n::tr("button.close"))
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
        crate::gentle_gui_profile_scope!("WindowDna::update");
        #[cfg(feature = "gui-test-support")]
        crate::gui_test_support::begin_viewport_frame(ctx);
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
                                Self::render_deferred_load_indicator(ui);
                            });
                        });
                    },
                );
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
                self.awaiting_first_content_paint = false;
            }
        }));
        if result.is_err() {
            eprintln!("E WindowDna: recovered from panic while rendering DNA window");
        }
    }

    pub fn update_embedded(&mut self, ui: &mut egui::Ui) {
        crate::gentle_gui_profile_scope!("WindowDna::update_embedded");
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
                                Self::render_deferred_load_indicator(ui);
                            });
                        });
                    },
                );
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
                self.render_bounded_embedded_main_area(ui);
                self.awaiting_first_content_paint = false;
            }
        }));
        if result.is_err() {
            eprintln!("E WindowDna: recovered from panic while rendering embedded DNA window");
        }
    }

    fn render_bounded_embedded_main_area(&mut self, ui: &mut egui::Ui) {
        crate::gentle_gui_profile_scope!("WindowDna::render_bounded_embedded_main_area");
        let body_rect = ui.available_rect_before_wrap();
        let body_id = egui::Id::new(("window_dna_embedded_body", self.main_area.panel_scope_key()));
        // Do not allocate the child's expanded min_rect in the parent: the DNA
        // panels can overflow while the hosted window remains user-sized.
        let mut body_ui = ui.new_child(
            egui::UiBuilder::new()
                .id_salt(body_id)
                .max_rect(body_rect)
                .layout(*ui.layout()),
        );
        body_ui.set_clip_rect(body_rect.intersect(ui.clip_rect()));
        self.main_area
            .render_inside_without_auxiliary_windows(&mut body_ui);
        ui.advance_cursor_after_rect(body_rect);
    }

    /// Keep detached auxiliary workspaces alive without redrawing the parent
    /// sequence viewport that originally spawned them.
    pub fn update_auxiliary_windows_only(&mut self, ctx: &egui::Context) {
        self.poll_deferred_load();
        if self.pending_dna_load.is_some() || self.deferred_load_message.is_some() {
            return;
        }
        let result = catch_unwind(AssertUnwindSafe(|| {
            self.main_area.render_auxiliary_windows_only(ctx);
        }));
        if result.is_err() {
            eprintln!("E WindowDna: recovered from panic while rendering detached auxiliaries");
        }
    }

    pub fn name(&self) -> String {
        self.main_area.window_title()
    }

    pub(crate) fn open_window_registry_signature(&self) -> Option<u64> {
        self.main_area.open_window_registry_signature()
    }

    pub fn sequence_id(&self) -> Option<String> {
        self.main_area.sequence_id().map(|v| v.to_string())
    }

    pub fn is_opening(&self) -> bool {
        self.deferred_load_message.is_none()
            && (self.pending_dna_load.is_some() || self.awaiting_first_content_paint)
    }

    pub fn opening_phase(&self) -> &'static str {
        if self.pending_dna_load.is_some() {
            "payload-loading"
        } else if self.deferred_load_message.is_some() {
            "failed"
        } else if self.awaiting_first_content_paint {
            "first-paint"
        } else {
            "ready"
        }
    }

    pub fn deferred_load_error(&self) -> Option<&str> {
        self.deferred_load_message.as_deref()
    }

    #[cfg(test)]
    pub(crate) fn mark_first_content_painted_for_tests(&mut self) {
        self.pending_dna_load = None;
        self.deferred_load_message = None;
        self.awaiting_first_content_paint = false;
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

    pub fn focus_feature_location_editor(&mut self, feature_index: Option<usize>) {
        self.main_area.focus_feature_location_editor(feature_index);
    }

    pub fn close_feature_location_editor(&mut self) -> bool {
        self.main_area.close_feature_location_editor()
    }

    pub fn feature_location_editor_is_open(&self) -> bool {
        self.main_area.feature_location_editor_is_open()
    }

    pub fn set_selection_range_0based(
        &mut self,
        start: usize,
        end_exclusive: usize,
    ) -> Result<(), String> {
        self.main_area
            .set_selection_range_0based(start, end_exclusive)
    }

    pub fn collect_open_auxiliary_window_entries(&self) -> Vec<(egui::ViewportId, String, String)> {
        self.main_area.collect_open_auxiliary_window_entries()
    }

    pub fn has_open_auxiliary_windows(&self) -> bool {
        !self
            .main_area
            .collect_open_auxiliary_window_entries()
            .is_empty()
    }

    pub fn embedded_auxiliary_window_layer_id(
        &self,
        viewport_id: egui::ViewportId,
    ) -> Option<egui::LayerId> {
        self.main_area
            .embedded_auxiliary_window_layer_id(viewport_id)
    }

    pub fn request_focus_auxiliary_window(&mut self, viewport_id: egui::ViewportId) -> bool {
        self.main_area.request_focus_auxiliary_window(viewport_id)
    }

    pub fn render_pcr_designer_specialist(&mut self, ui: &mut egui::Ui, ctx: &egui::Context) {
        self.main_area.set_rendering_pcr_specialist(true);
        self.main_area.render_pcr_designer_specialist(ui, ctx);
        self.main_area.set_rendering_pcr_specialist(false);
    }

    pub fn render_sequencing_confirmation_specialist(
        &mut self,
        ui: &mut egui::Ui,
        ctx: &egui::Context,
    ) {
        self.main_area
            .render_sequencing_confirmation_specialist(ui, ctx);
    }

    #[cfg(test)]
    pub(crate) fn seed_rna_read_mapping_window_for_tests(
        &mut self,
        seq_id: &str,
        feature_id: usize,
        group_label: &str,
    ) {
        self.main_area
            .seed_rna_read_mapping_window_for_tests(seq_id, feature_id, group_label);
    }

    #[cfg(test)]
    pub(crate) fn seed_splicing_expert_window_for_tests(
        &mut self,
        seq_id: &str,
        feature_id: usize,
        group_label: &str,
    ) {
        self.main_area
            .seed_splicing_expert_window_for_tests(seq_id, feature_id, group_label);
    }

    #[cfg(test)]
    pub(crate) fn seed_variant_followup_window_for_tests(
        &mut self,
        seq_id: &str,
        feature_id: usize,
        gene_label: &str,
    ) {
        self.main_area
            .seed_variant_followup_window_for_tests(seq_id, feature_id, gene_label);
    }

    #[cfg(test)]
    pub(crate) fn splicing_expert_focus_requested_for_tests(&self) -> bool {
        self.main_area.splicing_expert_focus_requested_for_tests()
    }

    #[cfg(test)]
    pub(crate) fn rna_read_mapping_focus_requested_for_tests(&self) -> bool {
        self.main_area.rna_read_mapping_focus_requested_for_tests()
    }

    #[cfg(test)]
    pub(crate) fn variant_followup_focus_requested_for_tests(&self) -> bool {
        self.main_area.variant_followup_focus_requested_for_tests()
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

    pub fn focus_cryptic_splicing_screen(&mut self) {
        if self.pending_dna_load.is_some() {
            self.deferred_analysis_focus = Some(DeferredAnalysisFocus::CrypticSplicing);
            return;
        }
        self.main_area.focus_cryptic_splicing_screen();
    }

    pub fn focus_primer_design_report(&mut self, report_id: &str) {
        if self.pending_dna_load.is_some() {
            self.deferred_analysis_focus =
                Some(DeferredAnalysisFocus::PrimerDesign(report_id.to_string()));
            return;
        }
        self.main_area.focus_primer_design_report(report_id);
    }

    pub fn focus_primer_specificity_report(&mut self, report_id: &str) {
        if self.pending_dna_load.is_some() {
            self.deferred_analysis_focus = Some(DeferredAnalysisFocus::PrimerSpecificity(
                report_id.to_string(),
            ));
            return;
        }
        self.main_area.focus_primer_specificity_report(report_id);
    }

    pub fn focus_rna_read_report(&mut self, report_id: &str) {
        if self.pending_dna_load.is_some() {
            self.deferred_analysis_focus =
                Some(DeferredAnalysisFocus::RnaReadReport(report_id.to_string()));
            return;
        }
        self.main_area.focus_rna_read_report(report_id);
    }

    pub fn focus_qpcr_design_report(&mut self, report_id: &str) {
        if self.pending_dna_load.is_some() {
            self.deferred_analysis_focus =
                Some(DeferredAnalysisFocus::QpcrDesign(report_id.to_string()));
            return;
        }
        self.main_area.focus_qpcr_design_report(report_id);
    }

    pub fn focus_restriction_cloning_handoff_report(&mut self, report_id: &str) {
        if self.pending_dna_load.is_some() {
            self.deferred_analysis_focus = Some(
                DeferredAnalysisFocus::RestrictionCloningPcrHandoff(report_id.to_string()),
            );
            return;
        }
        self.main_area
            .focus_restriction_cloning_handoff_report(report_id);
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

    pub fn focus_construct_reasoning_graph(&mut self, graph_id: &str) {
        if self.pending_dna_load.is_some() {
            self.deferred_analysis_focus = Some(DeferredAnalysisFocus::ConstructReasoning(
                graph_id.to_string(),
            ));
            return;
        }
        self.main_area.focus_construct_reasoning_graph(graph_id);
    }

    pub fn focus_uniprot_projection_expert(
        &mut self,
        projection_id: &str,
        protein_feature_filter: gentle_protocol::ProteinFeatureFilter,
    ) {
        if self.pending_dna_load.is_some() {
            self.deferred_analysis_focus = Some(DeferredAnalysisFocus::UniprotProjection {
                projection_id: projection_id.to_string(),
                protein_feature_filter,
            });
            return;
        }
        self.main_area
            .focus_uniprot_projection_expert(projection_id, protein_feature_filter);
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

    pub fn focus_ensembl_entry_protein_expert(
        &mut self,
        transcript_id_filter: Option<&str>,
        entry_id: &str,
        protein_feature_filter: gentle_protocol::ProteinFeatureFilter,
    ) {
        let normalized_entry_id = entry_id.trim();
        if normalized_entry_id.is_empty() {
            return;
        }
        if self.pending_dna_load.is_some() {
            self.deferred_analysis_focus = Some(DeferredAnalysisFocus::EnsemblProtein {
                transcript_id_filter: transcript_id_filter
                    .map(str::trim)
                    .filter(|value| !value.is_empty())
                    .map(|value| value.to_string()),
                entry_id: normalized_entry_id.to_string(),
                protein_feature_filter,
            });
            return;
        }
        self.main_area.focus_ensembl_entry_protein_expert(
            transcript_id_filter,
            normalized_entry_id,
            protein_feature_filter,
        );
    }

    pub fn refresh_from_engine_settings(&mut self) {
        self.main_area.refresh_from_engine_settings();
    }

    pub fn refresh_from_engine_state(&mut self) {
        self.main_area.refresh_from_engine_sequence_state();
        self.main_area.refresh_from_engine_settings();
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::engine::Engine;

    #[test]
    fn deferred_sequence_loading_uses_bounded_repaint_polling() {
        assert!(DEFERRED_LOAD_REPAINT_INTERVAL >= Duration::from_millis(50));
        assert!(DEFERRED_LOAD_REPAINT_INTERVAL <= Duration::from_millis(250));
    }

    #[test]
    fn deferred_sequence_loading_preserves_populated_restriction_site_catalog() {
        let mut engine = GentleEngine::default();
        let result = engine
            .apply(crate::engine::Operation::CreateSequenceFromText {
                sequence_text: "GAATTC".to_string(),
                output_id: Some("lazy_restriction_sites".to_string()),
                name: None,
                circular: false,
            })
            .expect("create prepared engine sequence");
        let seq_id = result
            .created_seq_ids
            .first()
            .expect("created sequence id")
            .clone();
        assert!(
            !engine
                .state()
                .sequences
                .get(&seq_id)
                .expect("created engine sequence")
                .restriction_enzyme_sites()
                .is_empty(),
            "the shared create operation must prepare the restriction-site catalog"
        );
        let engine = Arc::new(RwLock::new(engine));
        let mut window = WindowDna::new_lazy(seq_id, engine);
        let deadline = std::time::Instant::now() + Duration::from_secs(2);
        while window.pending_dna_load.is_some() && std::time::Instant::now() < deadline {
            window.poll_deferred_load();
            thread::sleep(Duration::from_millis(1));
        }

        assert!(
            window.pending_dna_load.is_none(),
            "deferred sequence load did not complete before the test deadline"
        );
        assert!(window.deferred_load_message.is_none());
        let loaded = window.main_area.dna().read().expect("loaded DNA");
        assert_eq!(loaded.get_forward_string(), "GAATTC");
        assert!(
            !loaded.restriction_enzyme_sites().is_empty(),
            "replace_loaded_sequence must leave the lazy window with the engine sequence's populated restriction-site catalog"
        );
    }
}
