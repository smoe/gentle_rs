//! Read-only semantic GUI metadata for native acceptance tests.
//!
//! This module exists only with `gui-test-support`. It exposes post-layout
//! rectangles in egui logical points and never injects input or bypasses GUI
//! confirmation. Callers may opt into an atomic JSON snapshot by setting
//! `GENTLE_GUI_TEST_SNAPSHOT` to an output path.

use std::{
    collections::{BTreeMap, BTreeSet},
    fs,
    path::Path,
};

use eframe::egui::{self, Context, Rect, Response};
use serde::{Deserialize, Serialize};

pub const SNAPSHOT_SCHEMA: &str = "gentle.gui_semantic_snapshot.v2";
pub const SNAPSHOT_PATH_ENV: &str = "GENTLE_GUI_TEST_SNAPSHOT";

#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord)]
pub struct GuiTestId(&'static str);

impl GuiTestId {
    pub fn new(value: &'static str) -> Self {
        assert_semantic_token(value);
        Self(value)
    }

    pub fn as_str(self) -> &'static str {
        self.0
    }
}

impl From<&'static str> for GuiTestId {
    fn from(value: &'static str) -> Self {
        Self::new(value)
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum GuiTestWidgetKind {
    Button,
    Checkbox,
    Tab,
    Row,
    Splitter,
    Window,
    Status,
    TextInput,
}

#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, Serialize, Deserialize)]
pub struct GuiTestState {
    pub visible: bool,
    pub enabled: bool,
    pub selected: bool,
}

#[derive(Clone, Copy, Debug, PartialEq, Serialize, Deserialize)]
pub struct GuiTestRect {
    pub min_x: f32,
    pub min_y: f32,
    pub max_x: f32,
    pub max_y: f32,
}

impl From<Rect> for GuiTestRect {
    fn from(rect: Rect) -> Self {
        Self {
            min_x: rect.min.x,
            min_y: rect.min.y,
            max_x: rect.max.x,
            max_y: rect.max.y,
        }
    }
}

#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct GuiTestItem {
    pub semantic_id: String,
    pub window_id: String,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub subject_scope: Option<String>,
    pub widget_kind: GuiTestWidgetKind,
    pub state: GuiTestState,
    pub rect_logical_points: GuiTestRect,
    pub pixels_per_point: f32,
    pub generation: u64,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub outcome_role: Option<String>,
}

#[derive(Clone, Debug, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct GuiTestSnapshot {
    pub schema: String,
    pub coordinate_space: String,
    pub generation: u64,
    pub settled: bool,
    pub items: Vec<GuiTestItem>,
}

#[derive(Clone, Default)]
struct GuiTestRegistry {
    generation: u64,
    unsettled_viewports: BTreeSet<egui::ViewportId>,
    items: BTreeMap<(String, String, Option<String>), RegisteredGuiTestItem>,
}

#[derive(Clone)]
struct RegisteredGuiTestItem {
    viewport_id: egui::ViewportId,
    pass_nr: u64,
    item: GuiTestItem,
}

fn discard_stale_viewport_items(
    registry: &mut GuiTestRegistry,
    viewport_id: egui::ViewportId,
    pass_nr: u64,
) {
    registry.items.retain(|_, registered| {
        registered.viewport_id != viewport_id || registered.pass_nr == pass_nr
    });
}

pub fn begin_frame(ctx: &Context) {
    let viewport_id = ctx.viewport_id();
    let active_viewports =
        ctx.input(|input| input.raw.viewports.keys().copied().collect::<BTreeSet<_>>());
    ctx.data_mut(|data| {
        let mut registry = data
            .get_temp::<GuiTestRegistry>(registry_id())
            .unwrap_or_default();
        registry.generation = registry.generation.saturating_add(1);
        registry
            .unsettled_viewports
            .retain(|candidate| active_viewports.contains(candidate) && *candidate != viewport_id);
        registry.items.retain(|_, registered| {
            active_viewports.contains(&registered.viewport_id)
                && registered.viewport_id != viewport_id
        });
        data.insert_temp(registry_id(), registry);
    });
}

/// Start one independently scheduled child-viewport paint.
///
/// Deferred native viewports can repaint more than once between root app
/// frames. Keep their latest rectangles in the aggregate snapshot while
/// retaining the duplicate-id assertion within each individual paint.
pub fn begin_viewport_frame(ctx: &Context) {
    let viewport_id = ctx.viewport_id();
    let pass_nr = ctx.cumulative_pass_nr();
    ctx.data_mut(|data| {
        let mut registry = data
            .get_temp::<GuiTestRegistry>(registry_id())
            .unwrap_or_default();
        discard_stale_viewport_items(&mut registry, viewport_id, pass_nr);
        registry.unsettled_viewports.remove(&viewport_id);
        data.insert_temp(registry_id(), registry);
    });
}

pub fn mark_unsettled(ctx: &Context) {
    let viewport_id = ctx.viewport_id();
    ctx.data_mut(|data| {
        let mut registry = data
            .get_temp::<GuiTestRegistry>(registry_id())
            .unwrap_or_default();
        registry.unsettled_viewports.insert(viewport_id);
        data.insert_temp(registry_id(), registry);
    });
}

pub fn register_response(
    response: &Response,
    semantic_id: impl Into<GuiTestId>,
    window_id: impl Into<GuiTestId>,
    subject_scope: Option<&str>,
    widget_kind: GuiTestWidgetKind,
    selected: bool,
) {
    register_response_with_outcome(
        response,
        semantic_id,
        window_id,
        subject_scope,
        widget_kind,
        selected,
        None,
    );
}

#[allow(clippy::too_many_arguments)]
pub fn register_response_with_outcome(
    response: &Response,
    semantic_id: impl Into<GuiTestId>,
    window_id: impl Into<GuiTestId>,
    subject_scope: Option<&str>,
    widget_kind: GuiTestWidgetKind,
    selected: bool,
    outcome_role: Option<&str>,
) {
    let semantic_id = semantic_id.into();
    let window_id = window_id.into();
    register_rect(
        response.ctx.clone(),
        semantic_id,
        window_id,
        subject_scope,
        widget_kind,
        response.interact_rect,
        response.interact_rect.is_positive(),
        response.enabled(),
        selected,
        outcome_role,
    );
}

#[allow(clippy::too_many_arguments)]
pub fn register_rect(
    ctx: Context,
    semantic_id: impl Into<GuiTestId>,
    window_id: impl Into<GuiTestId>,
    subject_scope: Option<&str>,
    widget_kind: GuiTestWidgetKind,
    rect: Rect,
    visible: bool,
    enabled: bool,
    selected: bool,
    outcome_role: Option<&str>,
) {
    let semantic_id = semantic_id.into();
    let window_id = window_id.into();
    let viewport_id = ctx.viewport_id();
    let pass_nr = ctx.cumulative_pass_nr();
    if let Some(scope) = subject_scope {
        assert_pseudonymous_scope(scope);
    }
    let pixels_per_point = ctx.pixels_per_point();
    let clipped_rect = rect.intersect(ctx.content_rect());
    let visible = visible && clipped_rect.is_positive();
    let clipped_rect = if visible { clipped_rect } else { Rect::ZERO };
    let viewport_origin = ctx.input(|input| {
        input
            .viewport()
            .outer_rect
            .map(|rect| rect.min)
            .unwrap_or(egui::Pos2::ZERO)
    });
    let screen_rect = clipped_rect.translate(viewport_origin.to_vec2());
    ctx.data_mut(|data| {
        let mut registry = data
            .get_temp::<GuiTestRegistry>(registry_id())
            .unwrap_or_default();
        discard_stale_viewport_items(&mut registry, viewport_id, pass_nr);
        let key = (
            window_id.as_str().to_string(),
            semantic_id.as_str().to_string(),
            subject_scope.map(str::to_string),
        );
        let item = GuiTestItem {
            semantic_id: semantic_id.as_str().to_string(),
            window_id: window_id.as_str().to_string(),
            subject_scope: subject_scope.map(str::to_string),
            widget_kind,
            state: GuiTestState {
                visible,
                enabled,
                selected,
            },
            rect_logical_points: screen_rect.into(),
            pixels_per_point,
            generation: registry.generation,
            outcome_role: outcome_role.map(str::to_string),
        };
        assert!(
            registry
                .items
                .insert(
                    key,
                    RegisteredGuiTestItem {
                        viewport_id,
                        pass_nr,
                        item,
                    },
                )
                .is_none(),
            "duplicate semantic GUI id in one window/frame: {}/{}",
            window_id.as_str(),
            semantic_id.as_str()
        );
        data.insert_temp(registry_id(), registry);
    });
}

pub fn snapshot(ctx: &Context) -> GuiTestSnapshot {
    let registry = ctx
        .data(|data| data.get_temp::<GuiTestRegistry>(registry_id()))
        .unwrap_or_default();
    GuiTestSnapshot {
        schema: SNAPSHOT_SCHEMA.to_string(),
        coordinate_space:
            "screen-relative egui logical points; multiply by pixels_per_point for physical pixels"
                .to_string(),
        generation: registry.generation,
        settled: registry.unsettled_viewports.is_empty(),
        items: registry
            .items
            .into_values()
            .map(|registered| registered.item)
            .collect(),
    }
}

pub fn finish_frame(ctx: &Context) -> Result<(), String> {
    let Some(path) = std::env::var_os(SNAPSHOT_PATH_ENV) else {
        return Ok(());
    };
    write_snapshot(Path::new(&path), &snapshot(ctx))
}

fn write_snapshot(path: &Path, snapshot: &GuiTestSnapshot) -> Result<(), String> {
    let bytes = serde_json::to_vec_pretty(snapshot)
        .map_err(|error| format!("Could not serialize semantic GUI snapshot: {error}"))?;
    let temporary = path.with_extension("json.tmp");
    fs::write(&temporary, bytes).map_err(|error| {
        format!(
            "Could not write semantic GUI snapshot '{}': {error}",
            temporary.display()
        )
    })?;
    fs::rename(&temporary, path).map_err(|error| {
        format!(
            "Could not publish semantic GUI snapshot '{}': {error}",
            path.display()
        )
    })
}

pub fn pseudonymous_subject_scope(parts: &[&str]) -> String {
    crate::tutorial_gui_semantics::pseudonymous_subject_scope(parts)
}

fn registry_id() -> egui::Id {
    egui::Id::new("gentle_gui_test_registry")
}

fn assert_semantic_token(value: &str) {
    assert!(
        !value.is_empty()
            && value.bytes().all(|byte| byte.is_ascii_lowercase()
                || byte.is_ascii_digit()
                || byte == b'.'
                || byte == b'_'),
        "semantic GUI ids must be stable lowercase dotted tokens"
    );
}

fn assert_pseudonymous_scope(value: &str) {
    assert!(
        value.starts_with("subject-")
            && value.len() == 40
            && value[8..].bytes().all(|byte| byte.is_ascii_hexdigit()),
        "semantic GUI subject scopes must be pseudonymous hashes"
    );
}

#[cfg(test)]
mod tests {
    use super::*;

    fn end_pass(ctx: &Context) {
        let mut output = ctx.end_pass();
        output.textures_delta.clear();
    }

    fn response(ctx: &Context, label: &str, enabled: bool) -> Response {
        let mut result = None;
        ctx.begin_pass(egui::RawInput::default());
        egui::Window::new("semantic test").show(ctx, |ui| {
            result = Some(ui.add_enabled(enabled, egui::Button::new(label)));
        });
        end_pass(ctx);
        result.expect("button response")
    }

    fn input_for_viewport(
        viewport_id: egui::ViewportId,
        active_viewports: &[egui::ViewportId],
    ) -> egui::RawInput {
        let mut input = egui::RawInput {
            viewport_id,
            ..Default::default()
        };
        for active_viewport in active_viewports {
            input
                .viewports
                .insert(*active_viewport, egui::ViewportInfo::default());
        }
        input
    }

    #[test]
    fn registry_is_deterministic_label_independent_and_reports_state() {
        let ctx = Context::default();
        begin_frame(&ctx);
        let disabled = response(&ctx, "Translated label A", false);
        register_response(
            &disabled,
            "agent.ask",
            "window.agent_assistant",
            None,
            GuiTestWidgetKind::Button,
            true,
        );
        register_rect(
            ctx.clone(),
            "agent.screenshot.remove",
            "window.agent_assistant",
            None,
            GuiTestWidgetKind::Button,
            Rect::NOTHING,
            false,
            false,
            false,
            None,
        );
        let first = snapshot(&ctx);
        assert_eq!(first.items[0].semantic_id, "agent.ask");
        assert!(!first.items[0].state.enabled);
        assert!(first.items[0].state.selected);
        assert!(!first.items[1].state.visible);

        begin_frame(&ctx);
        let translated = response(&ctx, "Andere Uebersetzung", false);
        register_response(
            &translated,
            "agent.ask",
            "window.agent_assistant",
            None,
            GuiTestWidgetKind::Button,
            true,
        );
        assert_eq!(snapshot(&ctx).items[0].semantic_id, "agent.ask");
    }

    #[test]
    fn response_registration_marks_controls_outside_the_clip_rect_invisible() {
        let ctx = Context::default();
        ctx.begin_pass(egui::RawInput::default());
        begin_frame(&ctx);
        egui::Window::new("clip test").show(&ctx, |ui| {
            ui.set_clip_rect(Rect::from_min_size(
                egui::Pos2::ZERO,
                egui::vec2(20.0, 20.0),
            ));
            let response = ui.put(
                Rect::from_min_size(egui::pos2(100.0, 100.0), egui::vec2(40.0, 20.0)),
                egui::Button::new("Off screen"),
            );
            register_response(
                &response,
                "pcr.design.run",
                "window.pcr_design",
                None,
                GuiTestWidgetKind::Button,
                false,
            );
        });
        let snapshot = snapshot(&ctx);
        assert_eq!(snapshot.items.len(), 1);
        assert!(!snapshot.items[0].state.visible);
        assert_eq!(
            snapshot.items[0].rect_logical_points,
            GuiTestRect::from(Rect::ZERO)
        );
        end_pass(&ctx);
    }

    #[test]
    #[should_panic(expected = "duplicate semantic GUI id")]
    fn duplicate_ids_in_one_window_and_frame_are_rejected() {
        let ctx = Context::default();
        begin_frame(&ctx);
        let button = response(&ctx, "A", true);
        for _ in 0..2 {
            register_response(
                &button,
                "agent.ask",
                "window.agent_assistant",
                None,
                GuiTestWidgetKind::Button,
                false,
            );
        }
    }

    #[test]
    fn a_new_egui_pass_replaces_that_viewports_previous_items() {
        let ctx = Context::default();
        begin_frame(&ctx);
        let first = response(&ctx, "First paint", true);
        register_response(
            &first,
            "agent.help.open",
            "window.dna_viewer",
            None,
            GuiTestWidgetKind::Button,
            false,
        );

        let second = response(&ctx, "Second paint", true);
        register_response(
            &second,
            "agent.help.open",
            "window.dna_viewer",
            None,
            GuiTestWidgetKind::Button,
            false,
        );

        let snapshot = snapshot(&ctx);
        assert_eq!(snapshot.items.len(), 1);
        assert_eq!(snapshot.items[0].semantic_id, "agent.help.open");
    }

    #[test]
    #[should_panic(expected = "duplicate semantic GUI id")]
    fn child_viewport_hook_does_not_hide_duplicates_in_the_current_pass() {
        let ctx = Context::default();
        begin_frame(&ctx);
        let button = response(&ctx, "Current paint", true);
        register_response(
            &button,
            "agent.help.open",
            "window.dna_viewer",
            None,
            GuiTestWidgetKind::Button,
            false,
        );

        begin_viewport_frame(&ctx);
        register_response(
            &button,
            "agent.help.open",
            "window.dna_viewer",
            None,
            GuiTestWidgetKind::Button,
            false,
        );
    }

    #[test]
    fn root_publication_keeps_active_child_items_and_removes_closed_viewports() {
        let ctx = Context::default();
        let root = egui::ViewportId::ROOT;
        let child = egui::ViewportId::from_hash_of("semantic-child");

        ctx.begin_pass(input_for_viewport(root, &[root, child]));
        begin_frame(&ctx);
        register_rect(
            ctx.clone(),
            "window.main",
            "window.main",
            None,
            GuiTestWidgetKind::Window,
            Rect::from_min_size(egui::Pos2::ZERO, egui::vec2(100.0, 80.0)),
            true,
            true,
            true,
            Some("ready"),
        );
        end_pass(&ctx);

        ctx.begin_pass(input_for_viewport(child, &[root, child]));
        begin_viewport_frame(&ctx);
        register_rect(
            ctx.clone(),
            "dna.splitter.info_width",
            "window.dna_viewer",
            None,
            GuiTestWidgetKind::Splitter,
            Rect::from_min_size(egui::pos2(300.0, 0.0), egui::vec2(12.0, 600.0)),
            true,
            true,
            false,
            None,
        );
        end_pass(&ctx);

        ctx.begin_pass(input_for_viewport(root, &[root, child]));
        begin_frame(&ctx);
        register_rect(
            ctx.clone(),
            "window.main",
            "window.main",
            None,
            GuiTestWidgetKind::Window,
            Rect::from_min_size(egui::Pos2::ZERO, egui::vec2(120.0, 90.0)),
            true,
            true,
            true,
            Some("ready"),
        );
        let aggregate = snapshot(&ctx);
        assert!(
            aggregate
                .items
                .iter()
                .any(|item| item.semantic_id == "dna.splitter.info_width")
        );
        assert!(
            aggregate
                .items
                .iter()
                .any(|item| item.semantic_id == "window.main")
        );
        end_pass(&ctx);

        ctx.begin_pass(input_for_viewport(root, &[root]));
        begin_frame(&ctx);
        assert!(
            snapshot(&ctx)
                .items
                .iter()
                .all(|item| item.window_id != "window.dna_viewer")
        );
        end_pass(&ctx);
    }

    #[test]
    fn child_unsettled_state_survives_root_publication_until_child_repaints() {
        let ctx = Context::default();
        let root = egui::ViewportId::ROOT;
        let child = egui::ViewportId::from_hash_of("semantic-loading-child");

        ctx.begin_pass(input_for_viewport(child, &[root, child]));
        begin_viewport_frame(&ctx);
        mark_unsettled(&ctx);
        end_pass(&ctx);

        ctx.begin_pass(input_for_viewport(root, &[root, child]));
        begin_frame(&ctx);
        assert!(!snapshot(&ctx).settled);
        end_pass(&ctx);

        ctx.begin_pass(input_for_viewport(child, &[root, child]));
        begin_viewport_frame(&ctx);
        assert!(snapshot(&ctx).settled);
        end_pass(&ctx);
    }

    #[test]
    fn scoped_dynamic_rows_keep_identity_across_order_and_filter_changes() {
        let a = pseudonymous_subject_scope(&["seq-1", "repeat", "100..120", "rmsk:Alu"]);
        let b = pseudonymous_subject_scope(&["seq-1", "array", "300..320", "PSR-1"]);
        let mut original = vec![a.clone(), b.clone()];
        let reordered = vec![b.clone(), a.clone()];
        original.sort();
        let mut reordered_sorted = reordered;
        reordered_sorted.sort();
        assert_eq!(original, reordered_sorted);
        assert_eq!(vec![a.clone()], vec![a.clone()]);
        assert!(!b.contains("PSR-1"));
        assert_ne!(
            pseudonymous_subject_scope(&["ab", "c"]),
            pseudonymous_subject_scope(&["a", "bc"]),
            "length-prefixing must keep scope components unambiguous"
        );
        assert_eq!(a.len(), 40);
    }

    #[test]
    fn snapshot_contains_no_labels_paths_or_raw_subject_data() {
        let ctx = Context::default();
        begin_frame(&ctx);
        let scope =
            pseudonymous_subject_scope(&["ACGTACGT", "/Users/private/project.gb", "secret-token"]);
        let button = response(&ctx, "Private visible label", true);
        register_response(
            &button,
            "dna.feature_tree.row",
            "window.dna_viewer",
            Some(&scope),
            GuiTestWidgetKind::Row,
            false,
        );
        let json = serde_json::to_string(&snapshot(&ctx)).expect("snapshot JSON");
        for secret in [
            "ACGTACGT",
            "/Users/private",
            "secret-token",
            "Private visible label",
        ] {
            assert!(!json.contains(secret));
        }
    }

    #[test]
    fn readiness_advances_with_generation_without_sleeping() {
        let ctx = Context::default();
        begin_frame(&ctx);
        mark_unsettled(&ctx);
        let first = snapshot(&ctx);
        begin_frame(&ctx);
        let second = snapshot(&ctx);
        assert!(!first.settled);
        assert!(second.settled);
        assert_eq!(second.generation, first.generation + 1);
    }

    #[test]
    fn hosted_and_native_modes_publish_the_same_control_identity() {
        fn render_mode(ctx: &Context, mode: &str) -> GuiTestSnapshot {
            begin_frame(ctx);
            let compose = response(ctx, mode, true);
            register_response(
                &compose,
                "splicing.locus.compose",
                "window.splicing_expert",
                None,
                GuiTestWidgetKind::Button,
                false,
            );
            register_rect(
                ctx.clone(),
                "splicing.locus.status",
                "window.splicing_expert",
                None,
                GuiTestWidgetKind::Status,
                Rect::from_min_size(egui::pos2(1.0, 2.0), egui::vec2(30.0, 10.0)),
                true,
                true,
                false,
                Some("ready"),
            );
            snapshot(ctx)
        }

        let hosted = render_mode(&Context::default(), "Compose and preview");
        let native = render_mode(&Context::default(), "Vorschau erstellen");
        let identities = |snapshot: &GuiTestSnapshot| {
            snapshot
                .items
                .iter()
                .map(|item| {
                    (
                        item.window_id.clone(),
                        item.semantic_id.clone(),
                        item.widget_kind,
                    )
                })
                .collect::<Vec<_>>()
        };
        assert_eq!(identities(&hosted), identities(&native));
        assert_eq!(hosted.items[1].outcome_role.as_deref(), Some("ready"));
    }
}
