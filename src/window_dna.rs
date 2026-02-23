use crate::{
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
            let nav_panel_id = egui::Id::new(("window_dna_nav", self.main_area.sequence_id()));
            egui::TopBottomPanel::top(nav_panel_id).show(ctx, |ui| {
                ui.horizontal(|ui| {
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
}
