//! Thin, read-only protocol graph over the shared routines-explain report.

use crate::planning_protocol::{
    ProtocolOverview, ProtocolOverviewStatus, from_explanation, legend, node_details,
};
use std::sync::Arc;

#[derive(Clone)]
struct ViewState {
    input_sha256: String,
    report: Result<Arc<ProtocolOverview>, String>,
    zoom: f32,
    selected: Option<String>,
}

fn semantic(response: &egui::Response, id: &str, subject: Option<&str>, selected: bool) {
    #[cfg(feature = "gui-test-support")]
    crate::gui_test_support::register_response(
        response,
        id,
        "window.routine_assistant",
        subject,
        crate::gui_test_support::GuiTestWidgetKind::Button,
        selected,
    );
    #[cfg(not(feature = "gui-test-support"))]
    let _ = (response, id, subject, selected);
}

pub(super) fn show(ui: &mut egui::Ui, expected_routine: &str, output: Option<&serde_json::Value>) {
    let Some(output) = output else {
        ui.label("Bench-step overview unavailable: no explanation loaded.");
        return;
    };
    if output["routine"]["routine_id"].as_str() != Some(expected_routine) {
        ui.label("Bench-step overview unavailable: explanation belongs to another routine.");
        return;
    }
    let key = egui::Id::new("routine_protocol_overview_view");
    // Hash the actual source and presentation bytes, not a claimed wire digest
    // or pointer identity. Verify only on change, outside the drawing loop.
    let bytes =
        serde_json::to_vec(&(output.get("routine"), output.get("planning"))).unwrap_or_default();
    let hash = crate::digest_utils::sha256_prefixed_bytes(&bytes);
    let previous = ui.ctx().data_mut(|data| data.get_temp::<ViewState>(key));
    let mut state = match previous {
        Some(previous) if previous.input_sha256 == hash => previous,
        _ => {
            let report = from_explanation(output).map(Arc::new);
            let zoom = report
                .as_ref()
                .map(|report| (ui.available_width() / report.width as f32).clamp(0.15, 1.0))
                .unwrap_or(1.0);
            ViewState {
                input_sha256: hash,
                report,
                zoom,
                selected: None,
            }
        }
    };
    match state.report.clone() {
        Err(error) => {
            ui.label(format!("Bench-step overview unavailable: {error}"));
        }
        Ok(report) => {
            ui.label(&report.reason);
            ui.small(legend(&report));
            if let Ok(estimate) = serde_json::from_value::<crate::engine::PlanningEstimate>(
                output["planning"]["estimate"].clone(),
            ) {
                for line in crate::engine_shell::planning_estimate_display_lines(&estimate)
                    .into_iter()
                    .take(4)
                {
                    ui.small(line);
                }
            }
            let fit_zoom = (ui.available_width() / report.width as f32).clamp(0.15, 1.0);
            ui.horizontal_wrapped(|ui| {
                ui.add(egui::Slider::new(&mut state.zoom, 0.15..=2.0).text("Zoom").suffix("x"));
                let fit = ui.button("Fit").on_hover_text("Fit the recorded dependency graph to the available width");
                semantic(&fit, "routine.protocol.fit", None, false);
                if fit.clicked() {
                    state.zoom = fit_zoom;
                }
                let copy_svg = ui.button("Copy SVG").on_hover_text("Copy the shared source-bound SVG, including step hover details");
                semantic(&copy_svg, "routine.protocol.copy_svg", None, false);
                if copy_svg.clicked() { ui.ctx().copy_text(report.svg.clone()); }
                let copy_json = ui.button("Copy JSON").on_hover_text("Copy the full explanation, including the estimate needed to verify this diagram");
                semantic(&copy_json, "routine.protocol.copy_json", None, false);
                if copy_json.clicked() { ui.ctx().copy_text(serde_json::to_string_pretty(output).unwrap_or_default()); }
            });
            if report.status == ProtocolOverviewStatus::Available {
                egui::ScrollArea::both()
                    .id_salt("routine_protocol_canvas")
                    .max_height(340.0)
                    .scroll_source(egui::scroll_area::ScrollSource::ALL)
                    .auto_shrink([false, false])
                    .show(ui, |ui| {
                        let zoom = state.zoom;
                        let size = egui::vec2(
                            report.width as f32 * zoom,
                            (report.height - 90) as f32 * zoom,
                        );
                        let (canvas, _) = ui.allocate_exact_size(size, egui::Sense::hover());
                        let point = |x: u32, y: u32| {
                            canvas.min + egui::vec2(x as f32, (y - 90) as f32) * zoom
                        };
                        let painter = ui.painter();
                        let grey = egui::Color32::from_rgb(110, 110, 110);
                        let critical = egui::Color32::from_rgb(177, 43, 98);
                        for edge in &report.edges {
                            let points = edge
                                .points
                                .iter()
                                .map(|p| point(p[0], p[1]))
                                .collect::<Vec<_>>();
                            let end = *points.last().unwrap();
                            let stroke = egui::Stroke::new(
                                1.5,
                                if edge.typical_critical_path {
                                    critical
                                } else {
                                    grey
                                },
                            );
                            painter.add(egui::Shape::line(points, stroke));
                            painter.arrow(
                                end - egui::vec2(12.0 * zoom, 0.0),
                                egui::vec2(12.0 * zoom, 0.0),
                                stroke,
                            );
                        }
                        for node in &report.nodes {
                            let rect = egui::Rect::from_min_size(
                                point(node.x, node.y),
                                egui::vec2(node.width as f32, node.height as f32) * zoom,
                            );
                            let selected =
                                state.selected.as_deref() == Some(node.step.step_id.as_str());
                            painter.rect_filled(rect, 4.0, egui::Color32::from_rgb(250, 250, 250));
                            painter.rect_stroke(
                                rect,
                                4.0,
                                egui::Stroke::new(
                                    if selected { 3.0 } else { 1.5 },
                                    if node.typical_critical_path {
                                        critical
                                    } else {
                                        grey
                                    },
                                ),
                                egui::StrokeKind::Inside,
                            );
                            for (index, line) in node.label_lines.iter().enumerate() {
                                painter.text(
                                    point(node.x + 12, node.y + 18 + index as u32 * 18),
                                    egui::Align2::LEFT_CENTER,
                                    line,
                                    egui::FontId::proportional(13.0 * zoom),
                                    egui::Color32::from_rgb(34, 34, 34),
                                );
                            }
                            for (y, text, color) in [
                                (
                                    node.y + 74,
                                    format!("Hands-on {:.2} h", node.step.hands_on_hours.typical),
                                    egui::Color32::from_rgb(36, 115, 75),
                                ),
                                (
                                    node.y + 94,
                                    format!(
                                        "Unattended {:.2} h",
                                        node.step.unattended_hours.typical
                                    ),
                                    egui::Color32::from_rgb(144, 100, 0),
                                ),
                            ] {
                                painter.text(
                                    point(node.x + 12, y),
                                    egui::Align2::LEFT_CENTER,
                                    text,
                                    egui::FontId::proportional(12.0 * zoom),
                                    color,
                                );
                            }
                            let response = ui
                                .interact(
                                    rect,
                                    ui.id().with(&node.step.step_id),
                                    egui::Sense::click(),
                                )
                                .on_hover_text(node_details(node));
                            semantic(
                                &response,
                                &format!("routine.protocol.step.{}", node.step.step_id),
                                Some(&node.step.step_id),
                                selected,
                            );
                            if response.clicked() {
                                state.selected = Some(node.step.step_id.clone());
                            }
                        }
                    });
                if let Some(node) = report
                    .nodes
                    .iter()
                    .find(|node| Some(&node.step.step_id) == state.selected.as_ref())
                {
                    ui.separator();
                    ui.label(node_details(node));
                }
            } else {
                ui.label("Bench-step overview unavailable");
            }
            egui::CollapsingHeader::new("Source and limits").show(ui, |ui| {
                ui.add(
                    egui::Label::new(
                        egui::RichText::new(&report.source_estimate_sha256).monospace(),
                    )
                    .wrap(),
                );
                for nonclaim in &report.nonclaims {
                    ui.label(nonclaim);
                }
            });
        }
    }
    ui.ctx().data_mut(|data| data.insert_temp(key, state));
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn protocol_overview_cache_verifies_changed_bytes_even_when_wire_digest_is_retained() {
        // Synthetic bundled bench defaults; no external protocol or screenshot.
        let mut engine = crate::engine::GentleEngine::default();
        let before = serde_json::to_value(engine.state()).unwrap();
        let command = crate::engine_shell::parse_shell_line(
            "routines explain restriction.digest_ligate_extract_sticky",
        )
        .unwrap();
        let mut output = crate::engine_shell::execute_shell_command(&mut engine, &command)
            .unwrap()
            .output;
        let ctx = egui::Context::default();
        for valid in [true, false] {
            if !valid {
                output["planning"]["protocol_overview"]["svg"] =
                    serde_json::json!("<svg>tampered</svg>");
            }
            ctx.begin_pass(egui::RawInput::default());
            crate::egui_compat::show_central_panel_for_test_context(
                &ctx,
                egui::CentralPanel::default(),
                |ui| {
                    show(
                        ui,
                        "restriction.digest_ligate_extract_sticky",
                        Some(&output),
                    );
                },
            );
            let _ = crate::egui_compat::end_test_pass(&ctx);
            let state = ctx
                .data_mut(|data| {
                    data.get_temp::<ViewState>(egui::Id::new("routine_protocol_overview_view"))
                })
                .unwrap();
            assert_eq!(state.report.is_ok(), valid);
        }
        assert_eq!(serde_json::to_value(engine.state()).unwrap(), before);
    }
}
