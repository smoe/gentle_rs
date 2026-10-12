//! Catalog-backed tutorial submenus, bounded to the current viewport.

use super::{GENtleApp, HelpTutorialDocEntry};
use eframe::egui;
use std::collections::BTreeMap;

fn config() -> egui::containers::menu::MenuConfig {
    egui::containers::menu::MenuConfig::new()
        .close_behavior(egui::PopupCloseBehavior::CloseOnClickOutside)
}

/// Scrollbar clicks must leave tutorial menus open; callers close on selection.
pub(super) fn submenu<'a>(
    ui: &mut egui::Ui,
    title: impl egui::IntoAtoms<'a>,
    contents: impl FnOnce(&mut egui::Ui),
) {
    egui::containers::menu::SubMenuButton::new(title)
        .config(config())
        .ui(ui, contents);
}

fn groups(entries: &[HelpTutorialDocEntry]) -> Vec<(String, Vec<usize>)> {
    let mut groups = BTreeMap::<(usize, String), Vec<usize>>::new();
    for (index, entry) in entries.iter().enumerate() {
        let label = GENtleApp::tutorial_audience_group_label(entry);
        let rank = GENtleApp::tutorial_audience_group_rank(entry, &label);
        groups.entry((rank, label)).or_default().push(index);
    }
    groups
        .into_iter()
        .map(|((_, label), indices)| (label, indices))
        .collect()
}

/// Return the original catalog index, not the index within a submenu.
pub(super) fn help_entries(ui: &mut egui::Ui, entries: &[HelpTutorialDocEntry]) -> Option<usize> {
    let mut selected = None;
    scroll(ui, |ui| {
        for (group, indices) in groups(entries) {
            ui.menu_button(format!("{group} ({})", indices.len()), |ui| {
                scroll(ui, |ui| {
                    for index in indices {
                        let entry = &entries[index];
                        let label = GENtleApp::tutorial_display_label(
                            entry.decimal_id.as_deref(),
                            None,
                            &entry.title,
                        );
                        if ui
                            .add(egui::Button::new(label).wrap())
                            .on_hover_text(format!(
                                "{}\n{}",
                                GENtleApp::help_tutorial_review_label(entry),
                                entry.summary
                            ))
                            .clicked()
                        {
                            selected = Some(index);
                        }
                    }
                });
            });
        }
    });
    selected
}

fn edge_scroll_offset(
    offset: f32,
    maximum: f32,
    rect: egui::Rect,
    pointer: egui::Pos2,
    dt: f32,
) -> f32 {
    if maximum <= 0.0 || !rect.contains(pointer) {
        return offset;
    }
    let band = 22.0_f32.min(rect.height() / 3.0);
    let direction = if pointer.y < rect.top() + band {
        -1.0
    } else if pointer.y > rect.bottom() - band {
        1.0
    } else {
        0.0
    };
    (offset + direction * 300.0 * dt.clamp(0.0, 0.05)).clamp(0.0, maximum)
}

/// Use egui scrolling for input/scrollbars, adding edge-hover movement only to
/// the menu under the pointer. Child submenus and scrollbar drags take priority.
pub(super) fn scroll<R>(
    ui: &mut egui::Ui,
    contents: impl FnOnce(&mut egui::Ui) -> R,
) -> egui::scroll_area::ScrollAreaOutput<R> {
    let viewport = ui.ctx().content_rect();
    let width = (viewport.width() - 48.0).clamp(1.0, 420.0);
    let height = (viewport.height() - 48.0).clamp(1.0, 480.0);
    ui.set_max_width(width);
    ui.style_mut().wrap_mode = Some(egui::TextWrapMode::Wrap);
    // Read before ScrollArea consumes the wheel/trackpad delta.
    let (pointer, busy, dt) = ui.input(|input| {
        (
            input.pointer.hover_pos(),
            input.pointer.any_down() || input.smooth_scroll_delta.y != 0.0,
            input.stable_dt,
        )
    });
    let mut output = egui::ScrollArea::vertical()
        .id_salt("tutorial_menu_scroll")
        .max_width(width)
        .max_height(height)
        .min_scrolled_height(1.0)
        .auto_shrink([true, true])
        .show(ui, contents);

    if !busy
        && output.state.velocity().y == 0.0
        && let Some(pointer) = pointer
        && ui.ctx().layer_id_at(pointer) == Some(ui.layer_id())
    {
        let rect = output.inner_rect.intersect(ui.clip_rect());
        let maximum = (output.content_size.y - output.inner_rect.height()).max(0.0);
        let offset = edge_scroll_offset(output.state.offset.y, maximum, rect, pointer, dt);
        if offset != output.state.offset.y {
            output.state.offset.y = offset;
            output.state.store(ui.ctx(), output.id);
            ui.ctx().request_repaint();
        }
    }
    output
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn tutorial_menu_groups_preserve_every_catalog_index_once() {
        let entries = GENtleApp::discover_help_tutorial_entries();
        assert!(!entries.is_empty());
        let grouped = groups(&entries);
        assert!(grouped.len() > 1 && grouped.len() < entries.len());
        let mut indices = Vec::new();
        let mut previous_rank = 0;
        for (label, members) in grouped {
            assert!(!members.is_empty());
            let rank = GENtleApp::tutorial_audience_group_rank(&entries[members[0]], &label);
            assert!(rank >= previous_rank);
            previous_rank = rank;
            for index in members {
                assert_eq!(
                    GENtleApp::tutorial_audience_group_label(&entries[index]),
                    label
                );
                indices.push(index);
            }
        }
        indices.sort_unstable();
        assert_eq!(indices, (0..entries.len()).collect::<Vec<_>>());
    }

    #[test]
    fn tutorial_menu_submenu_click_returns_original_entry_in_small_viewport() {
        let mut entries = GENtleApp::discover_help_tutorial_entries();
        entries.truncate(3);
        assert_eq!(entries.len(), 3);
        for (index, entry) in entries.iter_mut().enumerate() {
            entry.group_label = Some(if index == 1 { "Second" } else { "First" }.to_string());
            entry.group_order = Some(if index == 1 { 2 } else { 1 });
            entry.decimal_id = None;
            entry.title = format!(
                "Synthetic tutorial {index} with a long title that wraps inside the submenu"
            );
        }
        let ctx = egui::Context::default();
        let mut time = 0.0;
        let mut selected = None;
        let viewport = egui::Rect::from_min_size(egui::Pos2::ZERO, egui::vec2(440.0, 300.0));
        let mut draw = |events: Vec<egui::Event>| {
            time += 1.0 / 60.0;
            let output = ctx.run_ui(
                egui::RawInput {
                    screen_rect: Some(viewport),
                    time: Some(time),
                    events,
                    ..Default::default()
                },
                |ui| {
                    egui::MenuBar::new().ui(ui, |ui| {
                        egui::containers::menu::MenuButton::new("Tutorials")
                            .config(config())
                            .ui(ui, |ui| {
                                if let Some(index) = help_entries(ui, &entries) {
                                    selected = Some(index);
                                    ui.close();
                                }
                            });
                    });
                },
            );
            let mut texts = Vec::new();
            fn collect(
                shape: &egui::epaint::Shape,
                clip: egui::Rect,
                texts: &mut Vec<(String, egui::Rect)>,
            ) {
                match shape {
                    egui::epaint::Shape::Vec(shapes) => {
                        for shape in shapes {
                            collect(shape, clip, texts);
                        }
                    }
                    egui::epaint::Shape::Text(text)
                        if clip.intersects(text.visual_bounding_rect()) =>
                    {
                        texts.push((text.galley.text().to_string(), text.visual_bounding_rect()));
                    }
                    _ => {}
                }
            }
            for shape in &output.shapes {
                collect(&shape.shape, shape.clip_rect, &mut texts);
            }
            output.drop_without_applying_deltas();
            texts
        };
        draw(vec![]);
        let mut texts = draw(vec![]);
        for label in ["Tutorials", "Second (1)", entries[1].title.as_str()] {
            let rect = texts
                .iter()
                .find(|(text, _)| text == label)
                .unwrap_or_else(|| panic!("missing menu label {label:?}: {texts:?}"))
                .1;
            assert!(
                viewport.contains_rect(rect),
                "menu text outside viewport: {rect:?}"
            );
            for pressed in [true, false] {
                draw(vec![
                    egui::Event::PointerMoved(rect.center()),
                    egui::Event::PointerButton {
                        pos: rect.center(),
                        button: egui::PointerButton::Primary,
                        pressed,
                        modifiers: Default::default(),
                    },
                ]);
            }
            draw(vec![]);
            texts = draw(vec![]);
        }
        assert_eq!(selected, Some(1));
        assert!(!texts.iter().any(|(text, _)| text == "Second (1)"));
        assert!(!texts.iter().any(|(text, _)| text == &entries[1].title));
    }

    #[test]
    fn tutorial_menu_edge_scrolling_is_bounded_and_only_at_edges() {
        let rect = egui::Rect::from_min_size(egui::pos2(10.0, 10.0), egui::vec2(200.0, 200.0));
        let top = egui::pos2(50.0, 12.0);
        let bottom = egui::pos2(50.0, 208.0);
        assert_eq!(edge_scroll_offset(50.0, 100.0, rect, top, 0.02), 44.0);
        assert_eq!(edge_scroll_offset(50.0, 100.0, rect, bottom, 0.02), 56.0);
        assert_eq!(edge_scroll_offset(0.0, 100.0, rect, top, 0.02), 0.0);
        assert_eq!(edge_scroll_offset(98.0, 100.0, rect, bottom, 0.02), 100.0);
        assert_eq!(edge_scroll_offset(0.0, 0.0, rect, bottom, 0.02), 0.0);
        assert_eq!(
            edge_scroll_offset(50.0, 100.0, rect, rect.center(), 0.02),
            50.0
        );
        assert_eq!(
            edge_scroll_offset(50.0, 100.0, rect, egui::pos2(250.0, 208.0), 0.02),
            50.0
        );
        assert_eq!(edge_scroll_offset(0.0, 100.0, rect, bottom, 10.0), 15.0);
    }

    struct MenuFrame {
        area: egui::scroll_area::ScrollAreaOutput<Vec<egui::Rect>>,
        selected: Option<usize>,
    }

    // Hand-crafted GUI rows exercise overflow without loading or running a tutorial.
    fn frame(ctx: &egui::Context, time: &mut f64, events: Vec<egui::Event>) -> MenuFrame {
        *time += 1.0 / 60.0;
        let mut result = None;
        ctx.run_ui(
            egui::RawInput {
                screen_rect: Some(egui::Rect::from_min_size(egui::Pos2::ZERO, egui::vec2(360.0, 280.0))),
                time: Some(*time),
                events,
                ..Default::default()
            },
            |ui| {
                egui::Area::new(egui::Id::new("tutorial_menu_test"))
                    .fixed_pos(egui::pos2(8.0, 8.0))
                    .show(ui.ctx(), |ui| {
                        let mut selected = None;
                        let area = scroll(ui, |ui| {
                            (0..40)
                                .map(|index| {
                                    let response = ui.add(egui::Button::new(format!(
                                        "{index:02} A deliberately long synthetic tutorial title that must wrap in a small window"
                                    )).wrap());
                                    if response.clicked() {
                                        selected = Some(index);
                                    }
                                    response.rect
                                })
                                .collect()
                        });
                        result = Some(MenuFrame { area, selected });
                    });
            },
        ).drop_without_applying_deltas();
        result.expect("menu rendered")
    }

    #[test]
    fn tutorial_menu_small_window_reaches_and_selects_last_row_then_scrolls_back() {
        let ctx = egui::Context::default();
        let mut time = 0.0;
        frame(&ctx, &mut time, vec![]);
        let mut rendered = frame(&ctx, &mut time, vec![]);
        assert!(rendered.area.inner_rect.height() <= 232.0);
        assert!(rendered.area.content_size.x <= 312.0);
        assert!(rendered.area.content_size.y > rendered.area.inner_rect.height());
        let bottom = rendered.area.inner_rect.center_bottom() - egui::vec2(0.0, 2.0);
        for _ in 0..900 {
            rendered = frame(&ctx, &mut time, vec![egui::Event::PointerMoved(bottom)]);
            let maximum = rendered.area.content_size.y - rendered.area.inner_rect.height();
            if rendered.area.state.offset.y >= maximum {
                break;
            }
        }
        // The changed offset is painted in the following frame.
        rendered = frame(&ctx, &mut time, vec![]);
        let last = *rendered.area.inner.last().expect("last tutorial");
        assert!(
            rendered.area.inner_rect.contains_rect(last),
            "last row is reachable: {last:?}"
        );
        let click = last.center();
        for pressed in [true, false] {
            rendered = frame(
                &ctx,
                &mut time,
                vec![
                    egui::Event::PointerMoved(click),
                    egui::Event::PointerButton {
                        pos: click,
                        button: egui::PointerButton::Primary,
                        pressed,
                        modifiers: Default::default(),
                    },
                ],
            );
        }
        assert_eq!(rendered.selected, Some(39));
        let offset = rendered.area.state.offset.y;
        let top = rendered.area.inner_rect.center_top() + egui::vec2(0.0, 2.0);
        rendered = frame(&ctx, &mut time, vec![egui::Event::PointerMoved(top)]);
        assert!(rendered.area.state.offset.y < offset);
    }

    #[test]
    fn tutorial_menu_wheel_scrolling_and_drag_do_not_compete_with_edge_hover() {
        let ctx = egui::Context::default();
        let mut time = 0.0;
        frame(&ctx, &mut time, vec![]);
        let mut rendered = frame(&ctx, &mut time, vec![]);
        let center = rendered.area.inner_rect.center();
        rendered = frame(
            &ctx,
            &mut time,
            vec![
                egui::Event::PointerMoved(center),
                egui::Event::MouseWheel {
                    unit: egui::MouseWheelUnit::Point,
                    delta: egui::vec2(0.0, -90.0),
                    phase: egui::TouchPhase::Move,
                    modifiers: Default::default(),
                },
            ],
        );
        assert!(rendered.area.state.offset.y > 0.0);
        // The same wheel input at the upper edge must not be counteracted by
        // hover scrolling after egui consumes its input delta.
        let edge_ctx = egui::Context::default();
        let mut edge_time = 0.0;
        frame(&edge_ctx, &mut edge_time, vec![]);
        let edge = frame(&edge_ctx, &mut edge_time, vec![]);
        let edge = frame(
            &edge_ctx,
            &mut edge_time,
            vec![
                egui::Event::PointerMoved(edge.area.inner_rect.center_top() + egui::vec2(0.0, 2.0)),
                egui::Event::MouseWheel {
                    unit: egui::MouseWheelUnit::Point,
                    delta: egui::vec2(0.0, -90.0),
                    phase: egui::TouchPhase::Move,
                    modifiers: Default::default(),
                },
            ],
        );
        assert_eq!(edge.area.state.offset.y, rendered.area.state.offset.y);
        for _ in 0..90 {
            rendered = frame(&ctx, &mut time, vec![]);
        }
        let offset = rendered.area.state.offset.y;
        let bottom = rendered.area.inner_rect.center_bottom() - egui::vec2(0.0, 2.0);
        rendered = frame(
            &ctx,
            &mut time,
            vec![
                egui::Event::PointerMoved(bottom),
                egui::Event::PointerButton {
                    pos: bottom,
                    button: egui::PointerButton::Primary,
                    pressed: true,
                    modifiers: Default::default(),
                },
            ],
        );
        assert_eq!(rendered.area.state.offset.y, offset);
    }
}
