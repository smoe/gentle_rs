//! Presentation-only navigation to the current selection, including off-screen spans.

use super::*;

impl MainAreaDna {
    fn navigate_to_selection(&mut self, fit: bool) -> bool {
        if self.is_circular() {
            return false;
        }
        let Some((start, end)) = self.current_selection_range_0based() else {
            return false;
        };
        let (_, current_span, sequence_length) = self.current_linear_viewport();
        if start >= end || end > sequence_length {
            return false;
        }
        let span = if fit { end - start } else { current_span };
        let center = start + (end - start) / 2;
        self.set_linear_viewport(center.saturating_sub(span / 2), span);
        self.sync_linear_view_input_fields_to_viewport();
        self.fit_linear_features_in_view();
        self.map_sequence.request_scroll_to_selection();
        if !self.show_map {
            self.show_map = true;
            self.set_display_visibility(DisplayTarget::MapPanel, true);
        }
        if matches!(self.primary_map_mode, PrimaryMapMode::Dotplot) {
            self.primary_map_mode = PrimaryMapMode::Standard;
            self.save_engine_ops_state();
        }
        true
    }

    /// Both the toolbar and the map context menu navigate the same stored selection.
    pub(super) fn render_selection_navigation_controls(&mut self, ui: &mut egui::Ui) -> bool {
        let enabled = !self.is_circular() && self.current_selection_range_0based().is_some();
        let mut activated = false;
        for (fit, label, hint) in [
            (
                true,
                "sequence.zoom_selection",
                "sequence.zoom_selection_hover",
            ),
            (
                false,
                "sequence.go_selection",
                "sequence.go_selection_hover",
            ),
        ] {
            let response = ui
                .add_enabled(enabled, egui::Button::new(Self::tr(label)))
                .on_hover_text(Self::tr(hint));
            #[cfg(feature = "gui-test-support")]
            crate::gui_test_support::register_response(
                &response,
                if fit {
                    "dna.selection.zoom"
                } else {
                    "dna.selection.go"
                },
                crate::tutorial_gui_semantics::WINDOW_DNA_VIEWER,
                Some(&crate::gui_test_support::pseudonymous_subject_scope(&[
                    self.seq_id.as_deref().unwrap_or("unnamed"),
                ])),
                crate::gui_test_support::GuiTestWidgetKind::Button,
                false,
            );
            if response.clicked() {
                activated |= self.navigate_to_selection(fit);
            }
        }
        activated
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::engine::{GentleEngine, ProjectState};

    // Synthetic sequence and feature geometry; no downloaded biological fixtures.
    fn area() -> MainAreaDna {
        let dna = DNAsequence::from_sequence(&"ACGT".repeat(2500)).expect("sequence");
        MainAreaDna::new(dna, Some("navigation_test".to_owned()), None)
    }

    #[test]
    fn selection_navigation_fits_offscreen_formula_and_preserves_selection() {
        let mut area = area();
        area.set_linear_viewport(7000, 600);
        area.selection_formula_text = "=150 .. 950".into();
        area.apply_selection_formula();
        assert_eq!(area.current_linear_viewport(), (7000, 600, 10000));
        area.show_map = false;
        area.primary_map_mode = PrimaryMapMode::Dotplot;
        assert!(area.navigate_to_selection(true));
        assert_eq!(area.current_linear_viewport(), (150, 800, 10000));
        assert_eq!(area.current_selection_range_0based(), Some((150, 950)));
        assert_eq!(area.linear_view_start_1based_input, "151");
        assert_eq!(area.linear_view_end_1based_input, "950");
        assert!(area.show_map);
        assert!(matches!(area.primary_map_mode, PrimaryMapMode::Standard));
    }

    #[test]
    fn selection_navigation_go_preserves_zoom_and_clamps_at_sequence_edges() {
        let mut area = area();
        for (start, end, expected_start) in [(8000, 8100, 7750), (0, 1, 0), (9999, 10000, 9400)] {
            area.set_linear_viewport(3000, 600);
            area.set_selection_range_0based(start, end)
                .expect("selection");
            assert!(area.navigate_to_selection(false));
            assert_eq!(area.current_linear_viewport(), (expected_start, 600, 10000));
            assert_eq!(area.current_selection_range_0based(), Some((start, end)));
        }
        for (start, end) in [(0, 1), (9999, 10000), (0, 10000)] {
            area.set_selection_range_0based(start, end)
                .expect("selection");
            assert!(area.navigate_to_selection(true));
            assert_eq!(area.current_linear_viewport(), (start, end - start, 10000));
        }
    }

    #[test]
    fn selection_navigation_noops_without_selection_or_in_circular_view() {
        let mut area = area();
        area.set_linear_viewport(4000, 500);
        assert!(!area.navigate_to_selection(true));
        assert!(!area.navigate_to_selection(false));
        area.set_selection_range_0based(100, 200)
            .expect("selection");
        area.dna.write().expect("dna").set_circular(true);
        assert!(!area.navigate_to_selection(true));
        assert!(!area.navigate_to_selection(false));
        assert_eq!(area.current_linear_viewport(), (4000, 500, 10000));
    }

    #[test]
    fn selection_navigation_changes_only_project_display_settings() {
        let mut area = area();
        let mut state = ProjectState::default();
        state.sequences.insert(
            "navigation_test".into(),
            area.dna.read().expect("dna").clone(),
        );
        let engine = Arc::new(RwLock::new(GentleEngine::from_state(state)));
        area.engine = Some(engine.clone());
        let mut before =
            serde_json::to_value(engine.read().expect("engine").state()).expect("state");
        area.set_selection_range_0based(25, 450).expect("selection");
        assert!(area.navigate_to_selection(true));
        assert_eq!(
            engine
                .read()
                .expect("engine")
                .state()
                .display
                .linear_view_span_bp,
            425
        );
        let after = serde_json::to_value(engine.read().expect("engine").state()).expect("state");
        before["display"] = after["display"].clone();
        assert_eq!(before, after);
    }

    fn text_rects(shape: &egui::epaint::Shape, out: &mut Vec<(String, egui::Rect)>) {
        match shape {
            egui::epaint::Shape::Text(text) => {
                out.push((text.galley.job.text.clone(), text.visual_bounding_rect()))
            }
            egui::epaint::Shape::Vec(shapes) => {
                for shape in shapes {
                    text_rects(shape, out);
                }
            }
            _ => {}
        }
    }

    #[test]
    fn selection_navigation_buttons_work_after_manual_formula_application() {
        let mut area = area();
        area.set_linear_viewport(7000, 600);
        area.selection_formula_text = "=150 .. 950".into();
        area.apply_selection_formula();
        for (label, expected) in [
            ("Go to selection", (250, 600, 10000)),
            ("Zoom to selection", (150, 800, 10000)),
        ] {
            let ctx = egui::Context::default();
            let mut render = |events| {
                let mut out = ctx.run_ui(
                    egui::RawInput {
                        screen_rect: Some(egui::Rect::from_min_size(
                            egui::Pos2::ZERO,
                            egui::vec2(650.0, 240.0),
                        )),
                        events,
                        ..Default::default()
                    },
                    |ui| {
                        ui.horizontal_wrapped(|ui| {
                            area.render_selection_navigation_controls(ui);
                        });
                    },
                );
                let mut texts = Vec::new();
                for shape in &out.shapes {
                    text_rects(&shape.shape, &mut texts);
                }
                out.textures_delta.clear();
                texts
            };
            let texts = render(Vec::new());
            let pos = texts
                .iter()
                .find(|(text, _)| text == label)
                .expect("navigation button")
                .1
                .center();
            for pressed in [true, false] {
                render(vec![
                    egui::Event::PointerMoved(pos),
                    egui::Event::PointerButton {
                        pos,
                        button: egui::PointerButton::Primary,
                        pressed,
                        modifiers: egui::Modifiers::default(),
                    },
                ]);
            }
            assert_eq!(area.current_linear_viewport(), expected);
        }
    }

    #[test]
    fn dna_toolbar_wraps_at_window_width_and_leaves_room_for_the_map() {
        for width in [480.0, 900.0, 1300.0] {
            let ctx = egui::Context::default();
            let mut area = area();
            area.selection_formula_text = "=150 .. 950".into();
            area.apply_selection_formula();
            area.op_status = "Invalid selection_formula.start: resolved coordinate -617 is out of bounds for sequence length 10000".into();
            let screen = egui::Rect::from_min_size(egui::Pos2::ZERO, egui::vec2(width, 800.0));
            for frame in 0..4 {
                let mut out = ctx.run_ui(
                    egui::RawInput {
                        screen_rect: Some(screen),
                        ..Default::default()
                    },
                    |ui| {
                        egui::Panel::top("toolbar_test")
                            .frame(egui::Frame::NONE)
                            .show(ui, |ui| {
                                area.render_top_panel(ui);
                                assert!(
                                    ui.min_rect().width() <= width + 1.0,
                                    "width {width}, frame {frame}: {:?}",
                                    ui.min_rect()
                                );
                                assert!(
                                    ui.min_rect().height() <= 361.0,
                                    "toolbar must leave map space: {:?}",
                                    ui.min_rect()
                                );
                            });
                    },
                );
                let mut texts = Vec::new();
                for shape in &out.shapes {
                    text_rects(&shape.shape, &mut texts);
                }
                out.textures_delta.clear();
                for (text, rect) in texts {
                    assert!(
                        rect.right() <= width + 1.0,
                        "width {width}, frame {frame}: clipped {text:?} at {rect:?}"
                    );
                }
            }
        }
    }
}
