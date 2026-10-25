//! Shared GUI scroll/pan/zoom intent mapping and keyboard scroll helpers.
//!
//! This module provides a centralized input handling policy for scroll, pan, and zoom
//! interactions across all GENtle GUI canvases and viewports. It translates raw input
//! events (mouse wheel, trackpad, keyboard) into semantic intents that renderers can
//! consume consistently.
//!
//! ## Core Responsibilities
//!
//! - **Wheel/Trackpad Intent Classification**: Interprets mouse wheel and trackpad
//!   scroll events based on modifier keys. By default:
//!   - No modifier: Pan/Scroll in the scroll direction
//!   - Shift key: Zoom (in/out based on scroll direction)
//!
//! - **Keyboard Navigation**: Maps arrow keys to pan or zoom intents for
//!   zoomable canvases, with configurable step sizes.
//!
//! - **Intent Types**: Defines semantic intent types that decouple input handling
//!   from renderer logic:
//!   - `WheelIntent`: Wheel/trackpad scroll classification (Pan, Zoom, None)
//!   - `ArrowIntent`: Keyboard arrow key classification (Pan, Zoom, None)
//!   - `ZoomDirection`: Direction of zoom (In, Out)
//!
//! - **Configuration Constants**: Exposes tunable parameters for input sensitivity:
//!   - `DEFAULT_CANVAS_PAN_STEP`: Arrow key pan distance for zoomable canvases (40px)
//!   - `DEFAULT_SCROLLAREA_KEYBOARD_STEP`: Arrow key scroll distance for list/table panes (48px)
//!
//! ## Design Principles
//!
//! - **Modifier-Based Mode Switching**: Shift key toggles between pan and zoom modes
//!   for wheel input, providing a standard interaction pattern across all viewports.
//!
//! - **Dominant Axis Detection**: For zoom intents, uses the dominant axis (larger
//!   delta magnitude) to determine zoom direction, allowing both vertical and
//!   horizontal trackpad gestures to zoom consistently.
//!
//! - **Epsilon-Based Dead Zone**: Ignores scroll deltas below `DELTA_EPSILON` (0.0001)
//!   to prevent jitter from noise in input devices.
//!
//! - **Explicit Intent Passing**: Renderers receive semantic intents (`WheelIntent`,
//!   `ArrowIntent`) rather than raw input events, making behavior easier to test
//!   and modify.
//!
//! - **Cursor Feedback**: Returns appropriate `CursorIcon` values to provide visual
//!   feedback about the current interaction mode (e.g., grab cursor for pan, zoom
//!   cursor for zoom).
//!
//! ## Invariants
//!
//! - Wheel zoom direction: positive delta = zoom in, negative delta = zoom out
//! - Arrow pan direction: matches the arrow key direction
//! - Shift modifier always triggers zoom mode for wheel input
//! - Small deltas (< epsilon) are always ignored
//!
//! ## Failure Modes
//!
//! - Zero or NaN deltas are treated as no-op (not pan/zoom)
//! - Missing input state defaults to no modifiers (pan mode)
//! - Negative step sizes are treated as positive (absolute value)
//!
//! ## Usage Patterns
//!
//! Typical renderer integration:
//! ```ignore
//! let wheel_intent = scroll_input_policy::wheel_intent(input.scroll_delta, input.modifiers);
//! match wheel_intent {
//!     WheelIntent::Pan { delta } => { /* adjust viewport offset */ }
//!     WheelIntent::Zoom { direction, amount } => { /* adjust zoom level */ }
//!     WheelIntent::None => {}
//! }
//! ```
//!
//! For keyboard navigation in zoomable canvas:
//! ```ignore
//! let arrow_intent = scroll_input_policy::arrow_intent(&input, modifiers);
//! match arrow_intent {
//!     ArrowIntent::Pan { delta } => { /* adjust viewport offset */ }
//!     ArrowIntent::Zoom { direction } => { /* adjust zoom level */ }
//!     ArrowIntent::None => {}
//! }
//! ```

use eframe::egui::{self, CursorIcon, InputState, Key, Modifiers, Vec2};

const DELTA_EPSILON: f32 = 0.0001;

/// Default arrow-key pan step for zoomable canvases.
pub const DEFAULT_CANVAS_PAN_STEP: f32 = 40.0;

/// Default arrow-key scroll step for list/table/text panes.
pub const DEFAULT_SCROLLAREA_KEYBOARD_STEP: f32 = 48.0;

/// Direction of a zoom intent.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum ZoomDirection {
    In,
    Out,
}

/// Wheel/trackpad intent after applying modifier policy.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum WheelIntent {
    None,
    Pan {
        delta: Vec2,
    },
    Zoom {
        direction: ZoomDirection,
        amount: f32,
    },
}

/// Arrow-key intent on zoomable canvases.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum ArrowIntent {
    None,
    Pan { delta: Vec2 },
    Zoom { direction: ZoomDirection },
}

/// Classify wheel/trackpad input:
/// - no Shift => pan/scroll
/// - Shift => zoom
pub fn wheel_intent(delta: Vec2, modifiers: Modifiers) -> WheelIntent {
    if delta.x.abs() <= DELTA_EPSILON && delta.y.abs() <= DELTA_EPSILON {
        return WheelIntent::None;
    }
    if modifiers.shift {
        let dominant = if delta.y.abs() >= delta.x.abs() {
            delta.y
        } else {
            delta.x
        };
        let direction = if dominant >= 0.0 {
            ZoomDirection::In
        } else {
            ZoomDirection::Out
        };
        return WheelIntent::Zoom {
            direction,
            amount: dominant.abs(),
        };
    }
    WheelIntent::Pan { delta }
}

/// Lineage graph wheel classification:
/// - Shift zoom is primary
/// - Cmd/Ctrl wheel remains a compatibility zoom alias
pub fn lineage_graph_wheel_intent(delta: Vec2, modifiers: Modifiers) -> WheelIntent {
    if modifiers.shift || modifiers.command || modifiers.ctrl {
        return wheel_intent(
            delta,
            Modifiers {
                shift: true,
                ..modifiers
            },
        );
    }
    wheel_intent(delta, modifiers)
}

/// Map one arrow key press to pan (no Shift) or zoom (Shift) intent.
pub fn arrow_intent_for_key(key: Key, shift: bool, pan_step: f32) -> ArrowIntent {
    match (shift, key) {
        (true, Key::ArrowUp) | (true, Key::ArrowRight) => ArrowIntent::Zoom {
            direction: ZoomDirection::In,
        },
        (true, Key::ArrowDown) | (true, Key::ArrowLeft) => ArrowIntent::Zoom {
            direction: ZoomDirection::Out,
        },
        (false, Key::ArrowLeft) => ArrowIntent::Pan {
            delta: Vec2::new(-pan_step, 0.0),
        },
        (false, Key::ArrowRight) => ArrowIntent::Pan {
            delta: Vec2::new(pan_step, 0.0),
        },
        (false, Key::ArrowUp) => ArrowIntent::Pan {
            delta: Vec2::new(0.0, -pan_step),
        },
        (false, Key::ArrowDown) => ArrowIntent::Pan {
            delta: Vec2::new(0.0, pan_step),
        },
        _ => ArrowIntent::None,
    }
}

/// Aggregate arrow-key pan delta for this frame on a zoomable canvas.
pub fn canvas_keyboard_pan_delta(input: &InputState, pan_step: f32) -> Vec2 {
    if input.modifiers.shift {
        return Vec2::ZERO;
    }
    let mut delta = Vec2::ZERO;
    if input.key_down(Key::ArrowLeft) {
        delta.x -= pan_step;
    }
    if input.key_down(Key::ArrowRight) {
        delta.x += pan_step;
    }
    if input.key_down(Key::ArrowUp) {
        delta.y -= pan_step;
    }
    if input.key_down(Key::ArrowDown) {
        delta.y += pan_step;
    }
    delta
}

/// Aggregate zoom steps from Shift+arrow presses on a zoomable canvas.
pub fn canvas_keyboard_zoom_steps(input: &InputState) -> i32 {
    if !input.modifiers.shift {
        return 0;
    }
    let mut steps = 0;
    if input.key_pressed(Key::ArrowUp) {
        steps += 1;
    }
    if input.key_pressed(Key::ArrowRight) {
        steps += 1;
    }
    if input.key_pressed(Key::ArrowDown) {
        steps -= 1;
    }
    if input.key_pressed(Key::ArrowLeft) {
        steps -= 1;
    }
    steps
}

/// Keyboard scroll delta for non-zoom `ScrollArea`s.
///
/// Uses arrows for continuous scroll and PageUp/PageDown/Home/End for jumps.
pub fn scrollarea_keyboard_delta(input: &InputState, step: f32) -> Vec2 {
    let mut delta = Vec2::ZERO;
    if input.key_down(Key::ArrowUp) {
        delta.y += step;
    }
    if input.key_down(Key::ArrowDown) {
        delta.y -= step;
    }
    if input.key_down(Key::ArrowLeft) {
        delta.x += step;
    }
    if input.key_down(Key::ArrowRight) {
        delta.x -= step;
    }
    if input.key_pressed(Key::PageUp) {
        delta.y += step * 8.0;
    }
    if input.key_pressed(Key::PageDown) {
        delta.y -= step * 8.0;
    }
    if input.key_pressed(Key::Home) {
        delta.y += step * 16.0;
    }
    if input.key_pressed(Key::End) {
        delta.y -= step * 16.0;
    }
    delta
}

/// Apply keyboard scrolling to a hovered/non-text-edit scroll pane.
pub fn apply_scrollarea_keyboard_navigation(ui: &egui::Ui, step: f32) {
    if !ui.ui_contains_pointer() || ui.ctx().egui_wants_keyboard_input() {
        return;
    }
    let delta = ui.input(|i| scrollarea_keyboard_delta(i, step));
    if delta != Vec2::ZERO {
        ui.scroll_with_delta(delta);
    }
}

/// Primary hand-pan drag modifier (`Option`/`Alt`).
pub fn option_pan_modifier_active(modifiers: Modifiers) -> bool {
    modifiers.alt
}

/// Cursor icon policy for canvases with Option-pan and Shift-zoom affordances.
pub fn canvas_hover_cursor(
    modifiers: Modifiers,
    zoomable: bool,
    active_pan_drag: bool,
    active_node_drag: bool,
    last_wheel_intent: WheelIntent,
) -> Option<CursorIcon> {
    if active_pan_drag || active_node_drag {
        return Some(CursorIcon::Grabbing);
    }
    if option_pan_modifier_active(modifiers) {
        return Some(CursorIcon::Grab);
    }
    if zoomable && modifiers.shift {
        return match last_wheel_intent {
            WheelIntent::Zoom {
                direction: ZoomDirection::Out,
                ..
            } => Some(CursorIcon::ZoomOut),
            _ => Some(CursorIcon::ZoomIn),
        };
    }
    None
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn wheel_without_shift_maps_to_pan() {
        let intent = wheel_intent(Vec2::new(12.0, -6.0), Modifiers::NONE);
        assert_eq!(
            intent,
            WheelIntent::Pan {
                delta: Vec2::new(12.0, -6.0)
            }
        );
    }

    #[test]
    fn wheel_with_shift_maps_to_zoom_direction() {
        let intent_in = wheel_intent(
            Vec2::new(0.0, 8.0),
            Modifiers {
                shift: true,
                ..Modifiers::NONE
            },
        );
        let intent_out = wheel_intent(
            Vec2::new(0.0, -4.0),
            Modifiers {
                shift: true,
                ..Modifiers::NONE
            },
        );
        assert_eq!(
            intent_in,
            WheelIntent::Zoom {
                direction: ZoomDirection::In,
                amount: 8.0
            }
        );
        assert_eq!(
            intent_out,
            WheelIntent::Zoom {
                direction: ZoomDirection::Out,
                amount: 4.0
            }
        );
    }

    #[test]
    fn option_sets_pan_modifier_and_grab_cursor() {
        let modifiers = Modifiers {
            alt: true,
            ..Modifiers::NONE
        };
        assert!(option_pan_modifier_active(modifiers));
        assert_eq!(
            canvas_hover_cursor(modifiers, true, false, false, WheelIntent::None),
            Some(CursorIcon::Grab)
        );
    }

    #[test]
    fn shift_arrows_map_to_zoom_directions() {
        assert_eq!(
            arrow_intent_for_key(Key::ArrowUp, true, 10.0),
            ArrowIntent::Zoom {
                direction: ZoomDirection::In
            }
        );
        assert_eq!(
            arrow_intent_for_key(Key::ArrowRight, true, 10.0),
            ArrowIntent::Zoom {
                direction: ZoomDirection::In
            }
        );
        assert_eq!(
            arrow_intent_for_key(Key::ArrowDown, true, 10.0),
            ArrowIntent::Zoom {
                direction: ZoomDirection::Out
            }
        );
        assert_eq!(
            arrow_intent_for_key(Key::ArrowLeft, true, 10.0),
            ArrowIntent::Zoom {
                direction: ZoomDirection::Out
            }
        );
    }

    #[test]
    fn lineage_graph_wheel_accepts_cmd_ctrl_zoom_alias() {
        let cmd_modifiers = Modifiers {
            command: true,
            ..Modifiers::NONE
        };
        let ctrl_modifiers = Modifiers {
            ctrl: true,
            ..Modifiers::NONE
        };
        assert_eq!(
            lineage_graph_wheel_intent(Vec2::new(0.0, 6.0), cmd_modifiers),
            WheelIntent::Zoom {
                direction: ZoomDirection::In,
                amount: 6.0
            }
        );
        assert_eq!(
            lineage_graph_wheel_intent(Vec2::new(0.0, -3.0), ctrl_modifiers),
            WheelIntent::Zoom {
                direction: ZoomDirection::Out,
                amount: 3.0
            }
        );
    }
}
