//! DNA sequence-window wrapper and per-window controls.

use crate::{
    app::{request_open_graphics_settings_from_native_menu, request_open_help_from_native_menu},
    dna_sequence::DNAsequence,
    engine::GentleEngine,
    main_area_dna::MainAreaDna,
    window_backdrop::{
        WindowBackdropKind, current_window_backdrop_settings, paint_window_backdrop,
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
pub struct WindowDna {
    main_area: MainAreaDna,
    pending_dna_load: Option<Arc<Mutex<Receiver<Result<DNAsequence, String>>>>>,
    deferred_load_message: Option<String>,
}

impl WindowDna {
    pub fn new(dna: DNAsequence, seq_id: String, engine: Arc<RwLock<GentleEngine>>) -> Self {
        Self {
            main_area: MainAreaDna::new(dna, Some(seq_id), Some(engine)),
            pending_dna_load: None,
            deferred_load_message: None,
        }
    }

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
        }
    }

    pub fn update(&mut self, ctx: &egui::Context) {
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
            let nav_panel_id = egui::Id::new(("window_dna_nav", self.main_area.sequence_id()));
            egui::TopBottomPanel::top(nav_panel_id).show(ctx, |ui| {
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
                        .on_hover_text("Bring the main project window to front")
                        .clicked()
                    {
                        ctx.send_viewport_cmd_to(
                            egui::ViewportId::ROOT,
                            egui::ViewportCommand::Visible(true),
                        );
                        ctx.send_viewport_cmd_to(
                            egui::ViewportId::ROOT,
                            egui::ViewportCommand::Focus,
                        );
                    }
                    if ui
                        .button("Configuration")
                        .on_hover_text("Open Configuration window on Graphics settings")
                        .clicked()
                    {
                        request_open_graphics_settings_from_native_menu();
                    }
                });
            });
            let kind = if self.main_area.opened_from_pool_context() {
                WindowBackdropKind::Pool
            } else {
                WindowBackdropKind::Sequence
            };
            let settings = current_window_backdrop_settings();
            if self.pending_dna_load.is_some() {
                egui::CentralPanel::default()
                    .frame(egui::Frame::NONE)
                    .show(ctx, |ui| {
                        paint_window_backdrop(ui, kind, &settings);
                        ui.vertical_centered(|ui| {
                            ui.add_space(48.0);
                            ui.add(egui::Spinner::new());
                            ui.add_space(6.0);
                            ui.label(
                                "Opening sequence window: loading sequence content in background...",
                            );
                        });
                    });
                ctx.request_repaint();
            } else if let Some(message) = self.deferred_load_message.as_deref() {
                egui::CentralPanel::default()
                    .frame(egui::Frame::NONE)
                    .show(ctx, |ui| {
                        paint_window_backdrop(ui, kind, &settings);
                        ui.colored_label(egui::Color32::from_rgb(180, 32, 32), message);
                    });
            } else {
                // MainAreaDna owns the root panel layout for sequence windows.
                self.main_area.render(ctx);
            }
        }));
        if result.is_err() {
            eprintln!("E WindowDna: recovered from panic while rendering DNA window");
        }
    }

    pub fn name(&self) -> String {
        self.main_area.window_title()
    }

    pub fn sequence_id(&self) -> Option<String> {
        self.main_area.sequence_id().map(|v| v.to_string())
    }

    pub fn collect_open_auxiliary_window_entries(
        &self,
    ) -> Vec<(egui::ViewportId, String, String)> {
        self.main_area.collect_open_auxiliary_window_entries()
    }

    pub fn set_pool_context(&mut self, pool_seq_ids: Vec<String>) {
        self.main_area.set_pool_context(pool_seq_ids);
    }

    pub fn refresh_from_engine_settings(&mut self) {
        self.main_area.refresh_from_engine_settings();
    }

    pub fn refresh_from_engine_state(&mut self) {
        self.main_area.refresh_from_engine_sequence_state();
        self.main_area.refresh_from_engine_settings();
    }
}
