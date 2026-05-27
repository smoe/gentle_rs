//! Compatibility helpers for top-level egui panel entry points.
//!
//! `egui 0.34` prefers `show_inside(...)`, but the top-level panel allocation
//! path still relies on crate-private internals in egui itself. GENtle
//! therefore centralizes the remaining top-level `show(ctx, ...)` calls here so
//! the rest of the codebase stays on the newer surface and builds warning-free.

use eframe::egui;
use std::cell::{Cell, RefCell};
use std::collections::{HashMap, HashSet};
use std::ptr::NonNull;

thread_local! {
    static CURRENT_ROOT_UI: Cell<Option<NonNull<egui::Ui>>> = const { Cell::new(None) };
    static HOSTED_WINDOW_STALE_REPAINT_REQUESTS: RefCell<HashSet<egui::Id>> = RefCell::new(HashSet::new());
    static HOSTED_WINDOW_FRAME_DRAG_IDS: RefCell<HashMap<egui::Id, egui::Id>> = RefCell::new(HashMap::new());
    static HOSTED_WINDOW_LAYERS: RefCell<HashMap<egui::LayerId, egui::Id>> = RefCell::new(HashMap::new());
    static HOSTED_WINDOW_FRAME_DRAG_OWNER: Cell<Option<egui::Id>> = const { Cell::new(None) };
}

pub(crate) struct RootUiGuard {
    previous: Option<NonNull<egui::Ui>>,
}

impl Drop for RootUiGuard {
    fn drop(&mut self) {
        CURRENT_ROOT_UI.with(|current| {
            current.set(self.previous);
        });
    }
}

pub(crate) fn with_current_root_ui<R>(
    ui: &mut egui::Ui,
    add_contents: impl FnOnce(&egui::Context) -> R,
) -> R {
    let ctx = ui.ctx().clone();
    let _guard = install_current_root_ui(ui);
    add_contents(&ctx)
}

pub(crate) fn install_current_root_ui(ui: &mut egui::Ui) -> RootUiGuard {
    let previous = CURRENT_ROOT_UI.with(|current| current.replace(Some(NonNull::from(ui))));
    RootUiGuard { previous }
}

pub(crate) trait ContextHost {
    fn ctx(&self) -> &egui::Context;
}

impl ContextHost for &egui::Context {
    fn ctx(&self) -> &egui::Context {
        self
    }
}

impl ContextHost for &egui::Ui {
    fn ctx(&self) -> &egui::Context {
        egui::Ui::ctx(self)
    }
}

impl ContextHost for &mut egui::Ui {
    fn ctx(&self) -> &egui::Context {
        egui::Ui::ctx(self)
    }
}

pub(crate) trait PanelHost: ContextHost {
    fn with_panel_ui<R>(
        self,
        add_contents: impl FnOnce(&mut egui::Ui) -> egui::InnerResponse<R>,
    ) -> egui::InnerResponse<R>;
}

impl PanelHost for &mut egui::Ui {
    fn with_panel_ui<R>(
        self,
        add_contents: impl FnOnce(&mut egui::Ui) -> egui::InnerResponse<R>,
    ) -> egui::InnerResponse<R> {
        add_contents(self)
    }
}

impl PanelHost for &egui::Context {
    fn with_panel_ui<R>(
        self,
        add_contents: impl FnOnce(&mut egui::Ui) -> egui::InnerResponse<R>,
    ) -> egui::InnerResponse<R> {
        if let Some(mut current_ui) = CURRENT_ROOT_UI.with(Cell::get) {
            // `with_current_root_ui` installs a pointer only for the dynamic
            // extent of egui's own UI callback. This keeps the egui-main trial
            // compatible with the existing Context-based GENtle call surface
            // without cloning the whole window-rendering stack to pass `&mut Ui`.
            // TODO: replace this with explicit `&mut Ui` threading once the
            // egui-main root-panel API settles.
            return add_contents(unsafe { current_ui.as_mut() });
        }

        panic!(
            "show_central_panel(&Context) requires an active root Ui; call \
             egui_compat::with_current_root_ui or pass &mut Ui"
        )
    }
}

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

pub(crate) fn hosted_window_title_layer_visible(ctx: &egui::Context, title: &str) -> bool {
    let stale_title_layer = hosted_window_title_layer_id(title);
    ctx.memory(|mem| mem.areas().is_visible(&stale_title_layer))
}

#[cfg(test)]
pub(crate) fn show_legacy_layer_for_tests(
    ctx: &egui::Context,
    layer_id: egui::LayerId,
    add_contents: impl FnOnce(&mut egui::Ui),
) {
    egui::Area::new(layer_id.id)
        .order(layer_id.order)
        .show(ctx, add_contents);
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
    host: impl ContextHost,
    spec: &ModalWindowSpec,
    open: &mut bool,
    add_contents: impl FnOnce(&mut egui::Ui) -> R,
) -> Option<egui::InnerResponse<Option<R>>> {
    let ctx = host.ctx();
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
    host: impl ContextHost,
    spec: &HostedWindowSpec,
    open: &mut bool,
    add_contents: impl FnOnce(&mut egui::Ui) -> R,
) -> Option<egui::InnerResponse<Option<R>>> {
    let ctx = host.ctx();
    register_hosted_window_frame_drag_ids(spec);
    let frame_drag_owner = hosted_window_frame_drag_owner(ctx);
    let blocked_by_other_hosted_window_drag =
        frame_drag_owner.is_some_and(|owner| owner != spec.stable_id);
    let use_area_drag_for_title_move = hosted_window_title_area_drag_active(ctx, spec);
    let constrain_rect = hosted_window_safe_rect(ctx);
    let max_size = hosted_window_max_inner_size(constrain_rect, spec.min_size, spec.drag_margin);
    let stale_extra_visible = spec
        .legacy_layer_ids
        .iter()
        .any(|layer_id| ctx.memory(|mem| mem.areas().is_visible(layer_id)));
    let stale_title_visible = spec.cleanup_legacy_title_layer
        && hosted_window_title_layer_visible(ctx, spec.title.as_str());
    if stale_extra_visible || stale_title_visible {
        if should_request_hosted_window_stale_repaint(spec) {
            ctx.request_repaint();
        }
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
    if use_area_drag_for_title_move {
        window = window.drag_area(egui::WindowDrag::Anywhere);
    }
    if blocked_by_other_hosted_window_drag {
        window = window.interactable(false);
    }
    if spec.foreground {
        window = window.order(egui::Order::Foreground);
    }
    window.show(ctx, |ui| {
        if blocked_by_other_hosted_window_drag {
            ui.add_enabled_ui(false, add_contents).inner
        } else {
            add_contents(ui)
        }
    })
}

fn should_request_hosted_window_stale_repaint(spec: &HostedWindowSpec) -> bool {
    let stale_repaint_key = egui::Id::new((
        "hosted_window_stale_layer_repaint_requested",
        spec.stable_id,
    ));
    HOSTED_WINDOW_STALE_REPAINT_REQUESTS
        .with(|requests| requests.borrow_mut().insert(stale_repaint_key))
}

fn register_hosted_window_frame_drag_ids(spec: &HostedWindowSpec) {
    HOSTED_WINDOW_FRAME_DRAG_IDS.with(|registry| {
        let mut registry = registry.borrow_mut();
        for id in hosted_window_frame_drag_ids(spec) {
            registry.insert(id, spec.stable_id);
        }
    });
    HOSTED_WINDOW_LAYERS.with(|registry| {
        let mut registry = registry.borrow_mut();
        for order in [egui::Order::Middle, egui::Order::Foreground] {
            registry.insert(egui::LayerId::new(order, spec.stable_id), spec.stable_id);
        }
    });
}

fn hosted_window_frame_drag_owner(ctx: &egui::Context) -> Option<egui::Id> {
    if !ctx.input(|input| input.pointer.primary_down()) {
        HOSTED_WINDOW_FRAME_DRAG_OWNER.set(None);
        return None;
    }
    if let Some(owner) = HOSTED_WINDOW_FRAME_DRAG_OWNER.with(Cell::get) {
        return Some(owner);
    }
    if let Some(dragged_id) = ctx.dragged_id() {
        let owner = HOSTED_WINDOW_FRAME_DRAG_IDS
            .with(|registry| registry.borrow().get(&dragged_id).copied());
        if let Some(owner) = owner {
            HOSTED_WINDOW_FRAME_DRAG_OWNER.set(Some(owner));
            return Some(owner);
        }
    }
    let owner = hosted_window_press_origin_owner(ctx);
    if let Some(owner) = owner {
        HOSTED_WINDOW_FRAME_DRAG_OWNER.set(Some(owner));
    }
    owner
}

fn hosted_window_press_origin_owner(ctx: &egui::Context) -> Option<egui::Id> {
    let origin = ctx.input(|input| input.pointer.press_origin())?;
    ctx.memory(|memory| {
        let mut layer_ids = memory.layer_ids().collect::<Vec<_>>();
        layer_ids.reverse();
        for layer_id in layer_ids {
            let owner =
                HOSTED_WINDOW_LAYERS.with(|registry| registry.borrow().get(&layer_id).copied());
            let Some(owner) = owner else {
                continue;
            };
            if !memory.areas().is_visible(&layer_id) {
                continue;
            }
            let Some(rect) = memory.area_rect(layer_id.id) else {
                continue;
            };
            if hosted_window_capture_rect(ctx, rect).contains(origin) {
                return Some(owner);
            }
        }
        None
    })
}

fn hosted_window_capture_rect(ctx: &egui::Context, rect: egui::Rect) -> egui::Rect {
    let grab_radius = ctx.global_style().interaction.resize_grab_radius_side;
    rect.expand(grab_radius.max(0.0))
}

fn hosted_window_title_area_drag_active(ctx: &egui::Context, spec: &HostedWindowSpec) -> bool {
    if ctx.dragged_id() != Some(spec.stable_id.with("move")) {
        return false;
    }
    let Some(origin) = ctx.input(|input| input.pointer.press_origin()) else {
        return false;
    };
    let Some(rect) = ctx.memory(|mem| mem.area_rect(spec.stable_id)) else {
        return false;
    };
    let style = ctx.global_style();
    let grab_radius = style.interaction.resize_grab_radius_side.max(0.0);
    let title_height =
        (style.spacing.interact_size.y + 2.0 * style.spacing.item_spacing.y).clamp(24.0, 44.0);
    let title_rect = egui::Rect::from_min_max(
        egui::pos2(rect.left() + grab_radius, rect.top() + grab_radius),
        egui::pos2(rect.right() - grab_radius, rect.top() + title_height),
    );
    title_rect.contains(origin)
}

fn hosted_window_frame_drag_ids(spec: &HostedWindowSpec) -> Vec<egui::Id> {
    let mut ids = vec![
        spec.stable_id.with("move"),
        spec.stable_id.with("__title_click"),
    ];
    for order in [egui::Order::Middle, egui::Order::Foreground] {
        let edge_drag_id =
            egui::Id::new(egui::LayerId::new(order, spec.stable_id)).with("edge_drag");
        ids.push(edge_drag_id);
        for edge in [
            "left",
            "right",
            "top",
            "bottom",
            "right_bottom",
            "right_top",
            "left_bottom",
            "left_top",
        ] {
            ids.push(edge_drag_id.with(edge));
        }
    }
    ids
}

pub(crate) fn show_central_panel<R>(
    host: impl PanelHost,
    panel: egui::CentralPanel,
    add_contents: impl FnOnce(&mut egui::Ui) -> R,
) -> egui::InnerResponse<R> {
    host.with_panel_ui(|ui| panel.show(ui, add_contents))
}

#[cfg(test)]
pub(crate) fn show_central_panel_for_test_context<R>(
    ctx: &egui::Context,
    panel: egui::CentralPanel,
    add_contents: impl FnOnce(&mut egui::Ui) -> R,
) -> egui::InnerResponse<R> {
    let test_layer_id = egui::LayerId::new(
        egui::Order::Middle,
        egui::Id::new((ctx.viewport_id(), "__gentle_test_root_panel")),
    );
    let mut root_ui = egui::Ui::new(
        ctx.clone(),
        test_layer_id.id,
        egui::UiBuilder::new()
            .layer_id(test_layer_id)
            .max_rect(ctx.viewport_rect()),
    );
    show_central_panel(&mut root_ui, panel, add_contents)
}

pub(crate) fn show_central_panel_inside<R>(
    ui: &mut egui::Ui,
    panel: egui::CentralPanel,
    add_contents: impl FnOnce(&mut egui::Ui) -> R,
) -> egui::InnerResponse<R> {
    panel.show(ui, add_contents)
}

pub(crate) fn show_top_panel<R>(
    host: impl PanelHost,
    _root_id: egui::Id,
    panel: egui::Panel,
    add_contents: impl FnOnce(&mut egui::Ui) -> R,
) -> egui::InnerResponse<R> {
    host.with_panel_ui(|ui| panel.show(ui, add_contents))
}

pub(crate) fn show_top_panel_inside<R>(
    ui: &mut egui::Ui,
    panel: egui::Panel,
    add_contents: impl FnOnce(&mut egui::Ui) -> R,
) -> egui::InnerResponse<R> {
    panel.show(ui, add_contents)
}

pub(crate) fn show_bottom_panel<R>(
    host: impl PanelHost,
    _root_id: egui::Id,
    panel: egui::Panel,
    add_contents: impl FnOnce(&mut egui::Ui) -> R,
) -> egui::InnerResponse<R> {
    host.with_panel_ui(|ui| panel.show(ui, add_contents))
}

pub(crate) fn show_bottom_panel_inside<R>(
    ui: &mut egui::Ui,
    panel: egui::Panel,
    add_contents: impl FnOnce(&mut egui::Ui) -> R,
) -> egui::InnerResponse<R> {
    panel.show(ui, add_contents)
}

#[cfg(test)]
mod tests {
    use super::{
        clamp_hosted_window_default_pos, clamp_hosted_window_default_size,
        hosted_window_frame_drag_ids, hosted_window_frame_drag_owner, hosted_window_max_inner_size,
        hosted_window_press_origin_owner, hosted_window_safe_rect_for_rect,
        register_hosted_window_frame_drag_ids, should_request_hosted_window_stale_repaint,
        show_hosted_window, show_modal_window, HostedWindowSpec, ModalWindowSpec,
        HOSTED_WINDOW_SAFE_INSET_BOTTOM_PX, HOSTED_WINDOW_SAFE_INSET_TOP_PX,
        HOSTED_WINDOW_SAFE_INSET_X_PX,
    };
    use eframe::egui::{self, pos2, vec2, Rect};

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
    fn show_hosted_window_legacy_title_layer_expires_without_resetting_sibling() {
        let ctx = egui::Context::default();
        let stable_id = egui::Id::new("legacy_cleanup_stable_id");
        let sibling_id = egui::Id::new("legacy_cleanup_sibling_id");
        let spec = HostedWindowSpec::new(
            "Legacy Cleanup",
            stable_id,
            vec2(300.0, 220.0),
            vec2(160.0, 120.0),
        );
        let sibling_spec = HostedWindowSpec::new(
            "Sibling",
            sibling_id,
            vec2(360.0, 260.0),
            vec2(160.0, 120.0),
        )
        .initial_pos(Some(pos2(180.0, 130.0)));
        let stale_title_layer = super::hosted_window_title_layer_id("Legacy Cleanup");
        let mut open = true;
        let mut sibling_open = true;

        ctx.begin_pass(egui::RawInput::default());
        super::show_legacy_layer_for_tests(&ctx, stale_title_layer, |ui| {
            ui.label("legacy shell");
        });
        assert!(ctx.memory(|mem| mem.areas().is_visible(&stale_title_layer)));
        let _ = ctx.end_pass();

        ctx.begin_pass(egui::RawInput::default());
        show_hosted_window(&ctx, &sibling_spec, &mut sibling_open, |ui| {
            ui.label("sibling shell");
        });
        let sibling_before = ctx
            .memory(|mem| mem.area_rect(sibling_id))
            .expect("sibling area should be visible");
        show_hosted_window(&ctx, &spec, &mut open, |ui| {
            ui.label("stable shell");
        });
        let sibling_after = ctx
            .memory(|mem| mem.area_rect(sibling_id))
            .expect("sibling area should remain visible");
        assert_eq!(sibling_after, sibling_before);
        let _ = ctx.end_pass();
        assert!(!ctx.memory(|mem| mem.areas().is_visible(&stale_title_layer)));
    }

    #[test]
    fn show_hosted_window_legacy_extra_layer_expires_without_resetting_sibling() {
        let ctx = egui::Context::default();
        let stable_id = egui::Id::new("legacy_extra_cleanup_stable_id");
        let sibling_id = egui::Id::new("legacy_extra_cleanup_sibling_id");
        let stale_extra_id = egui::Id::new("legacy_extra_layer");
        let stale_extra_layer = egui::LayerId::new(egui::Order::Middle, stale_extra_id);
        let spec = HostedWindowSpec::new(
            "Legacy Extra Cleanup",
            stable_id,
            vec2(300.0, 220.0),
            vec2(160.0, 120.0),
        )
        .legacy_layer_id(stale_extra_layer);
        let sibling_spec = HostedWindowSpec::new(
            "Extra Sibling",
            sibling_id,
            vec2(360.0, 260.0),
            vec2(160.0, 120.0),
        )
        .initial_pos(Some(pos2(180.0, 130.0)));
        let mut open = true;
        let mut sibling_open = true;

        ctx.begin_pass(egui::RawInput::default());
        egui::Area::new(stale_extra_id).show(&ctx, |ui| {
            ui.label("legacy extra shell");
        });
        assert!(ctx.memory(|mem| mem.areas().is_visible(&stale_extra_layer)));
        let _ = ctx.end_pass();

        ctx.begin_pass(egui::RawInput::default());
        show_hosted_window(&ctx, &sibling_spec, &mut sibling_open, |ui| {
            ui.label("sibling shell");
        });
        let sibling_before = ctx
            .memory(|mem| mem.area_rect(sibling_id))
            .expect("sibling area should be visible");
        show_hosted_window(&ctx, &spec, &mut open, |ui| {
            ui.label("stable shell");
        });
        let sibling_after = ctx
            .memory(|mem| mem.area_rect(sibling_id))
            .expect("sibling area should remain visible");
        assert_eq!(sibling_after, sibling_before);
        let _ = ctx.end_pass();
        assert!(!ctx.memory(|mem| mem.areas().is_visible(&stale_extra_layer)));
    }

    #[test]
    fn show_hosted_window_stale_layer_repaint_request_is_one_shot() {
        let stable_id = egui::Id::new("stale_repaint_once_stable_id");
        let stale_extra_id = egui::Id::new("stale_repaint_once_legacy_layer");
        let stale_extra_layer = egui::LayerId::new(egui::Order::Middle, stale_extra_id);
        let spec = HostedWindowSpec::new(
            "Stale Repaint Once",
            stable_id,
            vec2(300.0, 220.0),
            vec2(160.0, 120.0),
        )
        .legacy_layer_id(stale_extra_layer);

        assert!(should_request_hosted_window_stale_repaint(&spec));
        assert!(!should_request_hosted_window_stale_repaint(&spec));
    }

    #[test]
    fn hosted_window_frame_drag_owner_tracks_registered_resize_edge() {
        let ctx = egui::Context::default();
        let owner_id = egui::Id::new("hosted_frame_drag_owner_window");
        let owner_spec =
            HostedWindowSpec::new("Owner", owner_id, vec2(300.0, 220.0), vec2(160.0, 120.0));
        register_hosted_window_frame_drag_ids(&owner_spec);
        let edge_drag_id = hosted_window_frame_drag_ids(&owner_spec)
            .into_iter()
            .find(|id| *id != owner_id.with("move") && *id != owner_id.with("__title_click"))
            .expect("edge drag id");

        ctx.begin_pass(egui::RawInput {
            events: vec![egui::Event::PointerButton {
                pos: pos2(10.0, 10.0),
                button: egui::PointerButton::Primary,
                pressed: true,
                modifiers: egui::Modifiers::default(),
            }],
            ..Default::default()
        });
        ctx.set_dragged_id(edge_drag_id);

        assert_eq!(hosted_window_frame_drag_owner(&ctx), Some(owner_id));
        let _ = ctx.end_pass();

        ctx.begin_pass(egui::RawInput {
            events: vec![egui::Event::PointerButton {
                pos: pos2(10.0, 10.0),
                button: egui::PointerButton::Primary,
                pressed: false,
                modifiers: egui::Modifiers::default(),
            }],
            ..Default::default()
        });
        assert_eq!(hosted_window_frame_drag_owner(&ctx), None);
        let _ = ctx.end_pass();
    }

    #[test]
    fn hosted_window_press_origin_owner_uses_topmost_registered_window() {
        let ctx = egui::Context::default();
        let lower_id = egui::Id::new("hosted_press_origin_lower_window");
        let upper_id = egui::Id::new("hosted_press_origin_upper_window");
        let lower_spec =
            HostedWindowSpec::new("Lower", lower_id, vec2(320.0, 240.0), vec2(160.0, 120.0))
                .initial_pos(Some(pos2(40.0, 40.0)));
        let upper_spec =
            HostedWindowSpec::new("Upper", upper_id, vec2(320.0, 240.0), vec2(160.0, 120.0))
                .initial_pos(Some(pos2(90.0, 90.0)));
        let mut lower_open = true;
        let mut upper_open = true;

        ctx.begin_pass(egui::RawInput::default());
        show_hosted_window(&ctx, &lower_spec, &mut lower_open, |ui| {
            ui.label("lower");
        });
        show_hosted_window(&ctx, &upper_spec, &mut upper_open, |ui| {
            ui.label("upper");
        });
        let _ = ctx.end_pass();

        ctx.begin_pass(egui::RawInput {
            events: vec![egui::Event::PointerButton {
                pos: pos2(100.0, 100.0),
                button: egui::PointerButton::Primary,
                pressed: true,
                modifiers: egui::Modifiers::default(),
            }],
            ..Default::default()
        });

        assert_eq!(hosted_window_press_origin_owner(&ctx), Some(upper_id));
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
    fn show_hosted_window_drag_margin_clamps_stale_oversized_area() {
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
            stale_rect.width() > max_size.x || stale_rect.height() > max_size.y,
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
            rect.width() <= max_size.x,
            "rect={rect:?}, max_size={max_size:?}"
        );
        assert!(
            rect.height() <= max_size.y,
            "rect={rect:?}, max_size={max_size:?}"
        );
        let _ = ctx.end_pass();
    }

    #[test]
    fn stale_oversized_hosted_window_does_not_reset_sibling_area() {
        let ctx = egui::Context::default();
        let sibling_id = egui::Id::new("sibling_hosted_layer");
        let stale_id = egui::Id::new("stale_hosted_layer");
        let mut sibling_open = true;
        let mut stale_open = true;
        let screen_rect = Rect::from_min_size(pos2(0.0, 0.0), vec2(1000.0, 720.0));
        let sibling_spec = HostedWindowSpec::new(
            "Sibling",
            sibling_id,
            vec2(420.0, 280.0),
            vec2(180.0, 120.0),
        )
        .initial_pos(Some(pos2(240.0, 160.0)));
        let stale_spec =
            HostedWindowSpec::new("Stale", stale_id, vec2(900.0, 620.0), vec2(180.0, 120.0))
                .drag_margin(vec2(96.0, 48.0));

        ctx.begin_pass(egui::RawInput {
            screen_rect: Some(screen_rect),
            ..Default::default()
        });
        show_hosted_window(
            &ctx,
            &HostedWindowSpec::new("Stale", stale_id, vec2(5000.0, 5000.0), vec2(180.0, 120.0)),
            &mut stale_open,
            |ui| {
                ui.set_min_size(vec2(5000.0, 5000.0));
                ui.label("stale full-size shell");
            },
        );
        let _ = ctx.end_pass();

        ctx.begin_pass(egui::RawInput {
            screen_rect: Some(screen_rect),
            ..Default::default()
        });
        show_hosted_window(&ctx, &sibling_spec, &mut sibling_open, |ui| {
            ui.label("sibling content");
        });
        let sibling_before = ctx
            .memory(|mem| mem.area_rect(sibling_id))
            .expect("sibling hosted area should be visible");
        show_hosted_window(&ctx, &stale_spec, &mut stale_open, |ui| {
            ui.label("normal content");
        });
        let sibling_after = ctx
            .memory(|mem| mem.area_rect(sibling_id))
            .expect("sibling hosted area should remain visible");
        assert_eq!(sibling_after, sibling_before);
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
        let drag_start = pos2(initial_rect.center().x, initial_rect.top() + 10.0);
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

        let drag_mid = drag_start + vec2(40.0, 0.0);
        ctx.begin_pass(egui::RawInput {
            screen_rect: Some(screen_rect),
            events: vec![egui::Event::PointerMoved(drag_mid)],
            ..Default::default()
        });
        show_hosted_window(&ctx, &spec, &mut open, |ui| {
            ui.label("content");
        });
        let _ = ctx.end_pass();

        let drag_end = drag_start + vec2(100.0, 0.0);
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
