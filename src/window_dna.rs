use crate::{
    app::request_open_help_from_native_menu,
    dna_sequence::DNAsequence,
    engine::GentleEngine,
    main_area_dna::MainAreaDna,
    window_backdrop::{
        WindowBackdropKind, current_window_backdrop_settings, paint_window_backdrop,
    },
};
use eframe::egui;
use std::panic::{AssertUnwindSafe, catch_unwind};
use std::sync::{Arc, RwLock};

#[derive(Clone, Debug)]
pub struct WindowDna {
    main_area: MainAreaDna,
}

impl WindowDna {
    pub fn new(dna: DNAsequence, seq_id: String, engine: Arc<RwLock<GentleEngine>>) -> Self {
        Self {
            main_area: MainAreaDna::new(dna, Some(seq_id), Some(engine)),
        }
    }

    pub fn update(&mut self, ctx: &egui::Context) {
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
                });
            });
            egui::CentralPanel::default().show(ctx, |ui| {
                let kind = if self.main_area.opened_from_pool_context() {
                    WindowBackdropKind::Pool
                } else {
                    WindowBackdropKind::Sequence
                };
                let settings = current_window_backdrop_settings();
                paint_window_backdrop(ui, kind, &settings);
                self.main_area.render(ctx, ui);
            });
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
