//! Compatibility helpers for top-level egui panel entry points.
//!
//! `egui 0.34` prefers `show_inside(...)`, but the top-level panel allocation
//! path still relies on crate-private internals in egui itself. GENtle
//! therefore centralizes the remaining top-level `show(ctx, ...)` calls here so
//! the rest of the codebase stays on the newer surface and builds warning-free.

use eframe::egui;

const HOSTED_WINDOW_SAFE_INSET_X_PX: f32 = 28.0;
const HOSTED_WINDOW_SAFE_INSET_TOP_PX: f32 = 36.0;
const HOSTED_WINDOW_SAFE_INSET_BOTTOM_PX: f32 = 42.0;

pub(crate) fn hosted_window_safe_rect_for_rect(rect: egui::Rect) -> egui::Rect {
    let shrunk = egui::Rect::from_min_max(
        egui::pos2(
            rect.min.x + HOSTED_WINDOW_SAFE_INSET_X_PX,
            rect.min.y + HOSTED_WINDOW_SAFE_INSET_TOP_PX,
        ),
        egui::pos2(
            rect.max.x - HOSTED_WINDOW_SAFE_INSET_X_PX,
            rect.max.y - HOSTED_WINDOW_SAFE_INSET_BOTTOM_PX,
        ),
    );
    if shrunk.width() >= 160.0 && shrunk.height() >= 120.0 {
        shrunk
    } else {
        rect
    }
}

pub(crate) fn hosted_window_safe_rect(ctx: &egui::Context) -> egui::Rect {
    hosted_window_safe_rect_for_rect(ctx.content_rect())
}

pub(crate) fn clamp_hosted_window_default_pos(
    desired: Option<egui::Pos2>,
    constrain_rect: egui::Rect,
    default_size: egui::Vec2,
) -> egui::Pos2 {
    let fallback = egui::pos2(constrain_rect.min.x + 16.0, constrain_rect.min.y + 12.0);
    let desired = desired.unwrap_or(fallback);
    let max_x = (constrain_rect.max.x - default_size.x).max(constrain_rect.min.x);
    let max_y = (constrain_rect.max.y - default_size.y).max(constrain_rect.min.y);
    egui::pos2(
        desired.x.clamp(constrain_rect.min.x, max_x),
        desired.y.clamp(constrain_rect.min.y, max_y),
    )
}

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

#[cfg(test)]
mod tests {
    use super::{
        HOSTED_WINDOW_SAFE_INSET_BOTTOM_PX, HOSTED_WINDOW_SAFE_INSET_TOP_PX,
        HOSTED_WINDOW_SAFE_INSET_X_PX, clamp_hosted_window_default_pos,
        hosted_window_safe_rect_for_rect,
    };
    use eframe::egui::{Rect, pos2, vec2};

    #[test]
    fn hosted_window_safe_rect_applies_expected_inset() {
        let rect = Rect::from_min_max(pos2(0.0, 0.0), pos2(1000.0, 800.0));
        let safe = hosted_window_safe_rect_for_rect(rect);
        assert_eq!(safe.min.x, HOSTED_WINDOW_SAFE_INSET_X_PX);
        assert_eq!(safe.min.y, HOSTED_WINDOW_SAFE_INSET_TOP_PX);
        assert_eq!(safe.max.x, 1000.0 - HOSTED_WINDOW_SAFE_INSET_X_PX);
        assert_eq!(safe.max.y, 800.0 - HOSTED_WINDOW_SAFE_INSET_BOTTOM_PX);
    }

    #[test]
    fn clamp_hosted_window_default_pos_keeps_window_inside_safe_rect() {
        let safe = Rect::from_min_max(pos2(28.0, 36.0), pos2(972.0, 758.0));
        let pos =
            clamp_hosted_window_default_pos(Some(pos2(900.0, 700.0)), safe, vec2(320.0, 240.0));
        assert!(pos.x <= safe.max.x - 320.0);
        assert!(pos.y <= safe.max.y - 240.0);
        assert!(pos.x >= safe.min.x);
        assert!(pos.y >= safe.min.y);
    }
}
