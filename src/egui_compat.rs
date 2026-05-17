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
const HOSTED_WINDOW_STALE_OUTER_SIZE_TOLERANCE_PX: f32 = 96.0;
const HOSTED_WINDOW_STALE_ROOT_OCCUPANCY_THRESHOLD: f32 = 0.92;

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

pub(crate) fn clamp_hosted_window_default_size(
    desired: egui::Vec2,
    constrain_rect: egui::Rect,
    min_size: egui::Vec2,
) -> egui::Vec2 {
    let max_width = constrain_rect.width().max(min_size.x);
    let max_height = constrain_rect.height().max(min_size.y);
    egui::vec2(
        desired.x.clamp(min_size.x, max_width),
        desired.y.clamp(min_size.y, max_height),
    )
}

pub(crate) fn hosted_window_max_inner_size(
    constrain_rect: egui::Rect,
    min_size: egui::Vec2,
    drag_margin: egui::Vec2,
) -> egui::Vec2 {
    let available_width = constrain_rect.width() - (drag_margin.x.max(0.0) * 2.0);
    let available_height = constrain_rect.height() - (drag_margin.y.max(0.0) * 2.0);
    egui::vec2(
        available_width.max(min_size.x),
        available_height.max(min_size.y),
    )
}

pub(crate) fn reset_hosted_window_area_if_stale_oversized(
    ctx: &egui::Context,
    stable_id: egui::Id,
    constrain_rect: egui::Rect,
    max_inner_size: egui::Vec2,
) -> bool {
    let Some(rect) = ctx.memory(|mem| mem.area_rect(stable_id)) else {
        return false;
    };
    let tolerance = HOSTED_WINDOW_STALE_OUTER_SIZE_TOLERANCE_PX;
    let over_width = rect.width() > max_inner_size.x + tolerance;
    let over_height = rect.height() > max_inner_size.y + tolerance;
    let fills_root_width =
        rect.width() > constrain_rect.width() * HOSTED_WINDOW_STALE_ROOT_OCCUPANCY_THRESHOLD;
    let fills_root_height =
        rect.height() > constrain_rect.height() * HOSTED_WINDOW_STALE_ROOT_OCCUPANCY_THRESHOLD;
    if !(over_width && fills_root_width || over_height && fills_root_height) {
        return false;
    }
    ctx.memory_mut(|mem| mem.reset_areas());
    ctx.request_repaint();
    true
}

#[derive(Clone, Debug)]
pub(crate) struct HostedWindowSpec {
    pub(crate) title: String,
    pub(crate) stable_id: egui::Id,
    pub(crate) default_size: egui::Vec2,
    pub(crate) min_size: egui::Vec2,
    pub(crate) drag_margin: egui::Vec2,
    pub(crate) initial_pos: Option<egui::Pos2>,
    pub(crate) resizable: bool,
    pub(crate) collapsible: bool,
    pub(crate) foreground: bool,
    pub(crate) cleanup_legacy_title_layer: bool,
    pub(crate) legacy_layer_ids: Vec<egui::LayerId>,
}

impl HostedWindowSpec {
    pub(crate) fn new(
        title: impl Into<String>,
        stable_id: egui::Id,
        default_size: egui::Vec2,
        min_size: egui::Vec2,
    ) -> Self {
        Self {
            title: title.into(),
            stable_id,
            default_size,
            min_size,
            drag_margin: egui::Vec2::ZERO,
            initial_pos: None,
            resizable: true,
            collapsible: false,
            foreground: false,
            cleanup_legacy_title_layer: true,
            legacy_layer_ids: vec![],
        }
    }

    pub(crate) fn initial_pos(mut self, initial_pos: Option<egui::Pos2>) -> Self {
        self.initial_pos = initial_pos;
        self
    }

    pub(crate) fn drag_margin(mut self, drag_margin: egui::Vec2) -> Self {
        self.drag_margin = egui::vec2(drag_margin.x.max(0.0), drag_margin.y.max(0.0));
        self
    }

    pub(crate) fn resizable(mut self, resizable: bool) -> Self {
        self.resizable = resizable;
        self
    }

    pub(crate) fn collapsible(mut self, collapsible: bool) -> Self {
        self.collapsible = collapsible;
        self
    }

    pub(crate) fn foreground(mut self, foreground: bool) -> Self {
        self.foreground = foreground;
        self
    }

    pub(crate) fn cleanup_legacy_title_layer(mut self, cleanup: bool) -> Self {
        self.cleanup_legacy_title_layer = cleanup;
        self
    }

    pub(crate) fn legacy_layer_id(mut self, layer_id: egui::LayerId) -> Self {
        self.legacy_layer_ids.push(layer_id);
        self
    }
}

pub(crate) fn hosted_window_title_layer_id(title: &str) -> egui::LayerId {
    egui::LayerId::new(egui::Order::Middle, egui::Id::new(title.to_string()))
}

pub(crate) fn reset_hosted_window_areas_if_legacy_title_layer_visible(
    ctx: &egui::Context,
    title: &str,
) -> bool {
    let stale_title_layer = hosted_window_title_layer_id(title);
    if ctx.memory(|mem| mem.areas().is_visible(&stale_title_layer)) {
        ctx.memory_mut(|mem| mem.reset_areas());
        true
    } else {
        false
    }
}

pub(crate) fn viewport_builder_for_hosted_window(spec: &HostedWindowSpec) -> egui::ViewportBuilder {
    egui::ViewportBuilder::default()
        .with_title(spec.title.clone())
        .with_inner_size([spec.default_size.x, spec.default_size.y])
        .with_min_inner_size([spec.min_size.x, spec.min_size.y])
}

#[derive(Clone, Debug)]
pub(crate) struct ModalWindowSpec {
    pub(crate) title: String,
    pub(crate) stable_id: egui::Id,
    pub(crate) default_size: Option<egui::Vec2>,
    pub(crate) min_size: Option<egui::Vec2>,
    pub(crate) resizable: bool,
    pub(crate) movable: bool,
    pub(crate) foreground: bool,
}

impl ModalWindowSpec {
    pub(crate) fn new(title: impl Into<String>, stable_id: egui::Id) -> Self {
        Self {
            title: title.into(),
            stable_id,
            default_size: None,
            min_size: None,
            resizable: false,
            movable: false,
            foreground: true,
        }
    }

    pub(crate) fn default_size(mut self, default_size: egui::Vec2) -> Self {
        self.default_size = Some(default_size);
        self
    }

    pub(crate) fn min_size(mut self, min_size: egui::Vec2) -> Self {
        self.min_size = Some(min_size);
        self
    }

    pub(crate) fn resizable(mut self, resizable: bool) -> Self {
        self.resizable = resizable;
        self
    }
}

pub(crate) fn show_modal_window<R>(
    ctx: &egui::Context,
    spec: &ModalWindowSpec,
    open: &mut bool,
    add_contents: impl FnOnce(&mut egui::Ui) -> R,
) -> Option<egui::InnerResponse<Option<R>>> {
    let constrain_rect = hosted_window_safe_rect(ctx);
    let min_size = spec.min_size.unwrap_or_else(|| egui::vec2(180.0, 96.0));
    let default_size = spec
        .default_size
        .map(|size| clamp_hosted_window_default_size(size, constrain_rect, min_size));
    let mut window = egui::Window::new(spec.title.clone())
        .id(spec.stable_id)
        .open(open)
        .collapsible(false)
        .movable(spec.movable)
        .resizable(spec.resizable)
        .anchor(egui::Align2::CENTER_CENTER, egui::Vec2::ZERO)
        .min_size(min_size)
        .max_size(constrain_rect.size())
        .constrain_to(constrain_rect);
    if let Some(default_size) = default_size {
        window = window.default_size(default_size);
    }
    if spec.foreground {
        window = window.order(egui::Order::Foreground);
    }
    window.show(ctx, add_contents)
}

pub(crate) fn show_hosted_window<R>(
    ctx: &egui::Context,
    spec: &HostedWindowSpec,
    open: &mut bool,
    add_contents: impl FnOnce(&mut egui::Ui) -> R,
) -> Option<egui::InnerResponse<Option<R>>> {
    let constrain_rect = hosted_window_safe_rect(ctx);
    let max_size = hosted_window_max_inner_size(constrain_rect, spec.min_size, spec.drag_margin);
    let stale_extra_visible = spec
        .legacy_layer_ids
        .iter()
        .any(|layer_id| ctx.memory(|mem| mem.areas().is_visible(layer_id)));
    let legacy_title_layer_reset = stale_extra_visible
        || (spec.cleanup_legacy_title_layer
            && reset_hosted_window_areas_if_legacy_title_layer_visible(ctx, spec.title.as_str()));
    if legacy_title_layer_reset {
        if stale_extra_visible {
            ctx.memory_mut(|mem| mem.reset_areas());
        }
        ctx.request_repaint();
    }
    if spec.drag_margin != egui::Vec2::ZERO {
        reset_hosted_window_area_if_stale_oversized(ctx, spec.stable_id, constrain_rect, max_size);
    }
    let default_size = egui::vec2(
        spec.default_size.x.clamp(spec.min_size.x, max_size.x),
        spec.default_size.y.clamp(spec.min_size.y, max_size.y),
    );
    let default_pos =
        clamp_hosted_window_default_pos(spec.initial_pos, constrain_rect, default_size);
    let mut window = egui::Window::new(spec.title.clone())
        .id(spec.stable_id)
        .open(open)
        .collapsible(spec.collapsible)
        .resizable(spec.resizable)
        .default_pos(default_pos)
        .default_size(default_size)
        .min_size(spec.min_size)
        .max_size(max_size)
        .constrain_to(constrain_rect);
    if spec.foreground {
        window = window.order(egui::Order::Foreground);
    }
    window.show(ctx, add_contents)
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
        HOSTED_WINDOW_SAFE_INSET_X_PX, HOSTED_WINDOW_STALE_OUTER_SIZE_TOLERANCE_PX,
        HOSTED_WINDOW_STALE_ROOT_OCCUPANCY_THRESHOLD, HostedWindowSpec, ModalWindowSpec,
        clamp_hosted_window_default_pos, clamp_hosted_window_default_size,
        hosted_window_max_inner_size, hosted_window_safe_rect_for_rect,
        reset_hosted_window_area_if_stale_oversized, show_hosted_window, show_modal_window,
    };
    use eframe::egui::{self, Rect, pos2, vec2};

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

    #[test]
    fn clamp_hosted_window_default_size_fits_large_defaults_into_safe_rect() {
        let safe = Rect::from_min_max(pos2(28.0, 36.0), pos2(972.0, 758.0));
        let size = clamp_hosted_window_default_size(vec2(1360.0, 900.0), safe, vec2(720.0, 420.0));
        assert_eq!(size.x, safe.width());
        assert_eq!(size.y, safe.height());
    }

    #[test]
    fn hosted_window_max_inner_size_never_drops_below_min_size() {
        let safe = Rect::from_min_size(pos2(0.0, 0.0), vec2(898.0, 540.0));
        let max_size = hosted_window_max_inner_size(safe, vec2(900.0, 560.0), vec2(0.0, 0.0));

        assert_eq!(max_size, vec2(900.0, 560.0));
    }

    #[test]
    fn show_hosted_window_tolerates_safe_rect_narrower_than_min_size() {
        let ctx = egui::Context::default();
        let mut open = true;
        let stable_id = egui::Id::new("narrow_hosted_window");
        let spec = HostedWindowSpec::new(
            "Narrow Host",
            stable_id,
            vec2(1280.0, 860.0),
            vec2(900.0, 560.0),
        );

        ctx.begin_pass(egui::RawInput {
            screen_rect: Some(Rect::from_min_size(pos2(0.0, 0.0), vec2(954.0, 638.0))),
            ..Default::default()
        });
        show_hosted_window(&ctx, &spec, &mut open, |ui| {
            ui.label("content");
        });
        assert!(ctx.memory(|mem| {
            mem.areas()
                .is_visible(&egui::LayerId::new(egui::Order::Middle, stable_id))
        }));
        let _ = ctx.end_pass();
    }

    #[test]
    fn show_hosted_window_uses_stable_id_not_visible_title_id() {
        let ctx = egui::Context::default();
        let stable_id = egui::Id::new("stable_hosted_window_test_id");
        let spec = HostedWindowSpec::new(
            "Stable Hosted",
            stable_id,
            vec2(300.0, 220.0),
            vec2(160.0, 120.0),
        );
        let stable_layer = egui::LayerId::new(egui::Order::Middle, stable_id);
        let title_layer = super::hosted_window_title_layer_id("Stable Hosted");
        let mut open = true;

        ctx.begin_pass(egui::RawInput::default());
        show_hosted_window(&ctx, &spec, &mut open, |ui| {
            ui.label("content");
        });
        assert!(ctx.memory(|mem| mem.areas().is_visible(&stable_layer)));
        assert!(!ctx.memory(|mem| mem.areas().is_visible(&title_layer)));
        let _ = ctx.end_pass();
    }

    #[test]
    fn show_hosted_window_resets_legacy_title_layer() {
        let ctx = egui::Context::default();
        let stable_id = egui::Id::new("legacy_cleanup_stable_id");
        let spec = HostedWindowSpec::new(
            "Legacy Cleanup",
            stable_id,
            vec2(300.0, 220.0),
            vec2(160.0, 120.0),
        );
        let stale_title_layer = super::hosted_window_title_layer_id("Legacy Cleanup");
        let mut open = true;

        ctx.begin_pass(egui::RawInput::default());
        egui::Window::new("Legacy Cleanup").show(&ctx, |ui| {
            ui.label("legacy shell");
        });
        assert!(ctx.memory(|mem| mem.areas().is_visible(&stale_title_layer)));
        show_hosted_window(&ctx, &spec, &mut open, |ui| {
            ui.label("stable shell");
        });
        assert!(!ctx.memory(|mem| mem.areas().is_visible(&stale_title_layer)));
        let _ = ctx.end_pass();
    }

    #[test]
    fn show_hosted_window_foreground_renders_above_middle_window() {
        let ctx = egui::Context::default();
        let normal_id = egui::Id::new("normal_hosted_layer");
        let foreground_id = egui::Id::new("foreground_hosted_layer");
        let mut normal_open = true;
        let mut foreground_open = true;

        ctx.begin_pass(egui::RawInput {
            screen_rect: Some(Rect::from_min_size(pos2(0.0, 0.0), vec2(800.0, 600.0))),
            ..Default::default()
        });
        show_hosted_window(
            &ctx,
            &HostedWindowSpec::new("Normal", normal_id, vec2(300.0, 220.0), vec2(160.0, 120.0)),
            &mut normal_open,
            |ui| {
                ui.label("normal");
            },
        );
        show_hosted_window(
            &ctx,
            &HostedWindowSpec::new(
                "Foreground",
                foreground_id,
                vec2(300.0, 220.0),
                vec2(160.0, 120.0),
            )
            .foreground(true),
            &mut foreground_open,
            |ui| {
                ui.label("foreground");
            },
        );
        assert!(ctx.memory(|mem| {
            mem.areas()
                .is_visible(&egui::LayerId::new(egui::Order::Middle, normal_id))
        }));
        assert!(ctx.memory(|mem| {
            mem.areas()
                .is_visible(&egui::LayerId::new(egui::Order::Foreground, foreground_id))
        }));
        let _ = ctx.end_pass();
    }

    #[test]
    fn show_hosted_window_drag_margin_resets_stale_oversized_area() {
        let ctx = egui::Context::default();
        let stable_id = egui::Id::new("drag_margin_hosted_layer");
        let spec = HostedWindowSpec::new(
            "Drag Margin",
            stable_id,
            vec2(900.0, 620.0),
            vec2(180.0, 120.0),
        )
        .drag_margin(vec2(96.0, 48.0));
        let mut open = true;
        let screen_rect = Rect::from_min_size(pos2(0.0, 0.0), vec2(1000.0, 720.0));
        let safe = hosted_window_safe_rect_for_rect(screen_rect);
        let max_size = hosted_window_max_inner_size(safe, spec.min_size, spec.drag_margin);

        ctx.begin_pass(egui::RawInput {
            screen_rect: Some(screen_rect),
            ..Default::default()
        });
        show_hosted_window(
            &ctx,
            &HostedWindowSpec::new(
                "Drag Margin",
                stable_id,
                vec2(5000.0, 5000.0),
                vec2(180.0, 120.0),
            ),
            &mut open,
            |ui| {
                ui.set_min_size(vec2(5000.0, 5000.0));
                ui.label("stale full-size shell");
            },
        );
        let stale_rect = ctx
            .memory(|mem| mem.area_rect(stable_id))
            .expect("stale hosted area should be visible");
        assert!(
            stale_rect.width() > max_size.x + HOSTED_WINDOW_STALE_OUTER_SIZE_TOLERANCE_PX
                || stale_rect.height() > max_size.y + HOSTED_WINDOW_STALE_OUTER_SIZE_TOLERANCE_PX,
            "stale_rect={stale_rect:?}, max_size={max_size:?}"
        );
        let _ = ctx.end_pass();

        ctx.begin_pass(egui::RawInput {
            screen_rect: Some(screen_rect),
            ..Default::default()
        });
        show_hosted_window(&ctx, &spec, &mut open, |ui| {
            ui.label("normal content");
        });
        let rect = ctx
            .memory(|mem| mem.area_rect(stable_id))
            .expect("hosted area should be visible");
        assert!(
            rect.width() <= max_size.x + HOSTED_WINDOW_STALE_OUTER_SIZE_TOLERANCE_PX,
            "rect={rect:?}, max_size={max_size:?}"
        );
        assert!(
            rect.height() <= max_size.y + HOSTED_WINDOW_STALE_OUTER_SIZE_TOLERANCE_PX,
            "rect={rect:?}, max_size={max_size:?}"
        );
        let _ = ctx.end_pass();
    }

    #[test]
    fn hosted_window_stale_reset_ignores_non_full_root_natural_width() {
        let ctx = egui::Context::default();
        let stable_id = egui::Id::new("natural_width_hosted_layer");
        let screen_rect = Rect::from_min_size(pos2(0.0, 0.0), vec2(1400.0, 900.0));
        let safe = hosted_window_safe_rect_for_rect(screen_rect);
        let max_size = hosted_window_max_inner_size(safe, vec2(820.0, 520.0), vec2(140.0, 72.0));
        let natural_width = max_size.x + 120.0;

        ctx.begin_pass(egui::RawInput {
            screen_rect: Some(screen_rect),
            ..Default::default()
        });
        egui::Area::new(stable_id)
            .fixed_pos(safe.min)
            .show(&ctx, |ui| {
                ui.set_min_size(vec2(natural_width, 240.0));
            });
        let rect = ctx
            .memory(|mem| mem.area_rect(stable_id))
            .expect("test area should be visible");
        assert!(
            rect.width() > max_size.x + HOSTED_WINDOW_STALE_OUTER_SIZE_TOLERANCE_PX,
            "rect={rect:?}, max_size={max_size:?}"
        );
        assert!(
            rect.width() < safe.width() * HOSTED_WINDOW_STALE_ROOT_OCCUPANCY_THRESHOLD,
            "rect={rect:?}, safe={safe:?}"
        );
        assert!(
            !reset_hosted_window_area_if_stale_oversized(&ctx, stable_id, safe, max_size),
            "natural-width content should not be treated like a stale full-root shell"
        );
        let retained_rect = ctx
            .memory(|mem| mem.area_rect(stable_id))
            .expect("test area should remain visible");
        assert_eq!(retained_rect, rect);
        let _ = ctx.end_pass();
    }

    #[test]
    fn show_hosted_window_with_drag_margin_moves_on_title_bar_drag() {
        let ctx = egui::Context::default();
        let stable_id = egui::Id::new("drag_room_move_layer");
        let spec = HostedWindowSpec::new(
            "Drag Room",
            stable_id,
            vec2(480.0, 320.0),
            vec2(180.0, 120.0),
        )
        .drag_margin(vec2(96.0, 48.0));
        let mut open = true;
        let screen_rect = Rect::from_min_size(pos2(0.0, 0.0), vec2(1200.0, 800.0));

        ctx.begin_pass(egui::RawInput {
            screen_rect: Some(screen_rect),
            ..Default::default()
        });
        show_hosted_window(&ctx, &spec, &mut open, |ui| {
            ui.label("content");
        });
        let initial_rect = ctx
            .memory(|mem| mem.area_rect(stable_id))
            .expect("hosted area should be visible");
        let drag_start = pos2(initial_rect.center().x, initial_rect.top() + 18.0);
        let _ = ctx.end_pass();

        ctx.begin_pass(egui::RawInput {
            screen_rect: Some(screen_rect),
            events: vec![
                egui::Event::PointerMoved(drag_start),
                egui::Event::PointerButton {
                    pos: drag_start,
                    button: egui::PointerButton::Primary,
                    pressed: true,
                    modifiers: egui::Modifiers::default(),
                },
            ],
            ..Default::default()
        });
        show_hosted_window(&ctx, &spec, &mut open, |ui| {
            ui.label("content");
        });
        let _ = ctx.end_pass();

        let drag_end = drag_start + vec2(80.0, 0.0);
        ctx.begin_pass(egui::RawInput {
            screen_rect: Some(screen_rect),
            events: vec![egui::Event::PointerMoved(drag_end)],
            ..Default::default()
        });
        show_hosted_window(&ctx, &spec, &mut open, |ui| {
            ui.label("content");
        });
        let moved_rect = ctx
            .memory(|mem| mem.area_rect(stable_id))
            .expect("hosted area should remain visible");
        assert!(
            moved_rect.min.x > initial_rect.min.x + 40.0,
            "initial_rect={initial_rect:?}, moved_rect={moved_rect:?}"
        );
        let _ = ctx.end_pass();

        ctx.begin_pass(egui::RawInput {
            screen_rect: Some(screen_rect),
            events: vec![egui::Event::PointerButton {
                pos: drag_end,
                button: egui::PointerButton::Primary,
                pressed: false,
                modifiers: egui::Modifiers::default(),
            }],
            ..Default::default()
        });
        show_hosted_window(&ctx, &spec, &mut open, |ui| {
            ui.label("content");
        });
        let _ = ctx.end_pass();
    }

    #[test]
    fn show_modal_window_uses_stable_foreground_layer() {
        let ctx = egui::Context::default();
        let modal_id = egui::Id::new("stable_modal_window_test_id");
        let spec = ModalWindowSpec::new("Stable Modal", modal_id);
        let modal_layer = egui::LayerId::new(egui::Order::Foreground, modal_id);
        let title_layer = super::hosted_window_title_layer_id("Stable Modal");
        let mut open = true;

        ctx.begin_pass(egui::RawInput::default());
        show_modal_window(&ctx, &spec, &mut open, |ui| {
            ui.label("modal content");
        });
        assert!(ctx.memory(|mem| mem.areas().is_visible(&modal_layer)));
        assert!(!ctx.memory(|mem| mem.areas().is_visible(&title_layer)));
        let _ = ctx.end_pass();
    }

    #[test]
    fn show_modal_window_clamps_default_size_to_small_safe_rect() {
        let ctx = egui::Context::default();
        let modal_id = egui::Id::new("small_modal_window_test_id");
        let spec = ModalWindowSpec::new("Small Modal", modal_id)
            .default_size(vec2(600.0, 460.0))
            .min_size(vec2(180.0, 96.0))
            .resizable(true);
        let mut open = true;

        ctx.begin_pass(egui::RawInput {
            screen_rect: Some(Rect::from_min_size(pos2(0.0, 0.0), vec2(320.0, 220.0))),
            ..Default::default()
        });
        show_modal_window(&ctx, &spec, &mut open, |ui| {
            ui.label("modal content");
        });
        let rect = ctx
            .memory(|mem| mem.area_rect(modal_id))
            .expect("modal area should be visible");
        assert!(rect.width() <= 320.0);
        assert!(rect.height() <= 220.0);
        let _ = ctx.end_pass();
    }

    #[test]
    fn show_modal_window_reconciles_open_state() {
        let ctx = egui::Context::default();
        let modal_id = egui::Id::new("closed_modal_window_test_id");
        let spec = ModalWindowSpec::new("Closed Modal", modal_id);
        let mut open = false;

        ctx.begin_pass(egui::RawInput::default());
        let response = show_modal_window(&ctx, &spec, &mut open, |ui| {
            ui.label("modal content");
        });
        assert!(response.is_none());
        assert!(!open);
        assert!(!ctx.memory(|mem| {
            mem.areas()
                .is_visible(&egui::LayerId::new(egui::Order::Foreground, modal_id))
        }));
        let _ = ctx.end_pass();
    }
}
