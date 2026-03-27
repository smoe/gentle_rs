//! Compatibility helpers for top-level egui panel entry points.
//!
//! `egui 0.34` deprecates `Panel::show(ctx, ...)` and `CentralPanel::show(ctx, ...)`
//! in favor of `show_inside(...)`. GENtle still has several multi-viewport/root
//! call sites that are simplest to keep in top-level form for now, so this
//! module centralizes the deprecated calls behind narrow helpers while the
//! broader UI migration continues.

use eframe::egui;

/// Show a top-level central panel without emitting per-callsite deprecation warnings.
#[allow(deprecated)]
pub(crate) fn show_central_panel<R>(
    ctx: &egui::Context,
    panel: egui::CentralPanel,
    add_contents: impl FnOnce(&mut egui::Ui) -> R,
) -> egui::InnerResponse<R> {
    panel.show(ctx, add_contents)
}

/// Show a top-level side/top/bottom panel without emitting per-callsite deprecation warnings.
#[allow(deprecated)]
pub(crate) fn show_panel<R>(
    ctx: &egui::Context,
    panel: egui::Panel,
    add_contents: impl FnOnce(&mut egui::Ui) -> R,
) -> egui::InnerResponse<R> {
    panel.show(ctx, add_contents)
}
