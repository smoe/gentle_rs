//! Compatibility helpers for top-level egui panel entry points.
//!
//! `egui 0.34` prefers `show_inside(...)`, but the top-level panel allocation
//! path still relies on crate-private internals in egui itself. GENtle
//! therefore centralizes the remaining top-level `show(ctx, ...)` calls here so
//! the rest of the codebase stays on the newer surface and builds warning-free.

use eframe::egui;

pub(crate) fn show_central_panel<R>(
    ctx: &egui::Context,
    panel: egui::CentralPanel,
    add_contents: impl FnOnce(&mut egui::Ui) -> R,
) -> egui::InnerResponse<R> {
    #[allow(deprecated)]
    {
        panel.show(ctx, add_contents)
    }
}

pub(crate) fn show_central_panel_inside<R>(
    ui: &mut egui::Ui,
    panel: egui::CentralPanel,
    add_contents: impl FnOnce(&mut egui::Ui) -> R,
) -> egui::InnerResponse<R> {
    panel.show_inside(ui, add_contents)
}

pub(crate) fn show_top_panel<R>(
    ctx: &egui::Context,
    _root_id: egui::Id,
    panel: egui::Panel,
    add_contents: impl FnOnce(&mut egui::Ui) -> R,
) -> egui::InnerResponse<R> {
    #[allow(deprecated)]
    {
        panel.show(ctx, add_contents)
    }
}

pub(crate) fn show_top_panel_inside<R>(
    ui: &mut egui::Ui,
    panel: egui::Panel,
    add_contents: impl FnOnce(&mut egui::Ui) -> R,
) -> egui::InnerResponse<R> {
    panel.show_inside(ui, add_contents)
}

pub(crate) fn show_bottom_panel<R>(
    ctx: &egui::Context,
    _root_id: egui::Id,
    panel: egui::Panel,
    add_contents: impl FnOnce(&mut egui::Ui) -> R,
) -> egui::InnerResponse<R> {
    #[allow(deprecated)]
    {
        panel.show(ctx, add_contents)
    }
}

pub(crate) fn show_bottom_panel_inside<R>(
    ui: &mut egui::Ui,
    panel: egui::Panel,
    add_contents: impl FnOnce(&mut egui::Ui) -> R,
) -> egui::InnerResponse<R> {
    panel.show_inside(ui, add_contents)
}
