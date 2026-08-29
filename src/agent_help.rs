//! User-invoked screenshots for window-bound Agent Assistant diagnostics.
//!
//! A capture request is created inside the viewport where the user clicked
//! `Agent help`. The Agent Assistant is opened only after that request has
//! produced an image, so focusing the assistant cannot change the target.

use eframe::egui;
use std::{
    collections::VecDeque,
    sync::{
        LazyLock, Mutex,
        atomic::{AtomicU64, Ordering},
    },
};

static NEXT_CAPTURE_REQUEST_ID: AtomicU64 = AtomicU64::new(1);
static CAPTURE_EVENTS: LazyLock<Mutex<VecDeque<AgentHelpCaptureEvent>>> =
    LazyLock::new(|| Mutex::new(VecDeque::new()));
static CONSUMED_EGUI_CAPTURE_IDS: LazyLock<Mutex<VecDeque<u64>>> =
    LazyLock::new(|| Mutex::new(VecDeque::new()));

#[derive(Clone, Debug)]
struct EguiCaptureToken {
    request_id: u64,
    window_title: String,
}

#[derive(Clone, Debug)]
pub struct AgentHelpCapturedImage {
    pub request_id: u64,
    pub window_title: String,
    pub backend: String,
    pub image: egui::ColorImage,
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub enum AgentHelpCaptureFailureKind {
    PermissionRequired,
    RestartRequired,
    Unsupported,
    CaptureFailed,
}

#[derive(Clone, Debug)]
pub struct AgentHelpCaptureFailure {
    pub request_id: u64,
    pub window_title: String,
    pub kind: AgentHelpCaptureFailureKind,
    pub message: String,
}

#[derive(Clone, Debug)]
pub enum AgentHelpCaptureEvent {
    Captured(AgentHelpCapturedImage),
    Failed(AgentHelpCaptureFailure),
}

fn next_request_id() -> u64 {
    NEXT_CAPTURE_REQUEST_ID.fetch_add(1, Ordering::Relaxed)
}

fn push_capture_event(event: AgentHelpCaptureEvent) {
    if let Ok(mut events) = CAPTURE_EVENTS.lock() {
        events.push_back(event);
    }
}

fn claim_egui_capture(request_id: u64) -> bool {
    let Ok(mut consumed) = CONSUMED_EGUI_CAPTURE_IDS.lock() else {
        return false;
    };
    if consumed.contains(&request_id) {
        return false;
    }
    consumed.push_back(request_id);
    if consumed.len() > 256 {
        consumed.pop_front();
    }
    true
}

pub fn take_capture_events() -> Vec<AgentHelpCaptureEvent> {
    CAPTURE_EVENTS
        .lock()
        .map(|mut events| events.drain(..).collect())
        .unwrap_or_default()
}

/// Collect screenshot replies for requests issued from this egui viewport.
pub fn collect_egui_capture_events(ctx: &egui::Context) {
    collect_egui_capture_events_for(ctx, ctx.viewport_id());
}

/// Collect screenshot replies for a specific registered egui viewport.
pub fn collect_egui_capture_events_for(ctx: &egui::Context, viewport_id: egui::ViewportId) {
    let captures = ctx.input_for(viewport_id, |input| {
        input
            .events
            .iter()
            .filter_map(|event| {
                let egui::Event::Screenshot {
                    user_data, image, ..
                } = event
                else {
                    return None;
                };
                let token = user_data
                    .data
                    .as_ref()?
                    .downcast_ref::<EguiCaptureToken>()?;
                if !claim_egui_capture(token.request_id) {
                    return None;
                }
                Some((token.clone(), image.as_ref().clone()))
            })
            .collect::<Vec<_>>()
    });
    if !captures.is_empty() {
        ctx.request_repaint_of(viewport_id);
    }
    for (token, image) in captures {
        push_capture_event(AgentHelpCaptureEvent::Captured(AgentHelpCapturedImage {
            request_id: token.request_id,
            window_title: token.window_title,
            backend: "egui.viewport".to_string(),
            image,
        }));
    }
}

pub fn request_egui_viewport_capture(ctx: &egui::Context, window_title: impl Into<String>) -> u64 {
    request_egui_viewport_capture_for(ctx, ctx.viewport_id(), window_title)
}

/// Request one screenshot from an explicitly selected registered egui viewport.
pub fn request_egui_viewport_capture_for(
    ctx: &egui::Context,
    viewport_id: egui::ViewportId,
    window_title: impl Into<String>,
) -> u64 {
    let token = EguiCaptureToken {
        request_id: next_request_id(),
        window_title: window_title.into(),
    };
    let request_id = token.request_id;
    ctx.send_viewport_cmd_to(
        viewport_id,
        egui::ViewportCommand::Screenshot(egui::UserData::new(token)),
    );
    ctx.request_repaint_of(viewport_id);
    request_id
}

/// Render the common help affordance. Left-click captures the current GENtle
/// viewport; on macOS, the context menu offers native full-window capture.
pub fn render_agent_help_button(ui: &mut egui::Ui, window_title: impl Into<String>) {
    collect_egui_capture_events(ui.ctx());
    let window_title = window_title.into();
    let response = ui
        .button(crate::i18n::tr("agent.help_button"))
        .on_hover_text(crate::i18n::tr("agent.help_button.tooltip"));
    #[cfg(feature = "gui-test-support")]
    {
        let semantic_window = if window_title.starts_with("Splicing Expert") {
            "window.splicing_expert"
        } else if window_title.starts_with("Agent Assistant") {
            "window.agent_assistant"
        } else {
            "window.dna_viewer"
        };
        let subject_scope = crate::gui_test_support::pseudonymous_subject_scope(&[&window_title]);
        crate::gui_test_support::register_response(
            &response,
            "agent.help.open",
            semantic_window,
            Some(&subject_scope),
            crate::gui_test_support::GuiTestWidgetKind::Button,
            false,
        );
    }
    if response.clicked() {
        request_egui_viewport_capture(ui.ctx(), window_title.clone());
    }
    response.context_menu(|ui| {
        if ui
            .button(crate::i18n::tr("agent.help_native_window"))
            .on_hover_text(crate::i18n::tr("agent.help_native_window.tooltip"))
            .clicked()
        {
            request_native_window_capture(ui.ctx().clone(), window_title.clone());
            ui.close();
        }
    });
}

#[cfg(all(target_os = "macos", feature = "screenshot-capture"))]
#[link(name = "CoreGraphics", kind = "framework")]
unsafe extern "C" {
    fn CGPreflightScreenCaptureAccess() -> bool;
    fn CGRequestScreenCaptureAccess() -> bool;
}

#[cfg(all(target_os = "macos", feature = "screenshot-capture"))]
fn active_macos_window_target() -> Result<(u32, f64), String> {
    use objc2::MainThreadMarker;
    use objc2_app_kit::NSApplication;

    let Some(mtm) = MainThreadMarker::new() else {
        return Err(
            "Native window capture must be requested from the macOS main thread".to_string(),
        );
    };
    let app = NSApplication::sharedApplication(mtm);
    let Some(window) = app.keyWindow().or_else(|| app.mainWindow()) else {
        return Err("No active GENtle window is available for capture".to_string());
    };
    let raw_window_id = window.windowNumber();
    let window_id = u32::try_from(raw_window_id)
        .map_err(|_| format!("Active GENtle window reported invalid id {raw_window_id}"))?;
    Ok((window_id, window.backingScaleFactor() as f64))
}

#[cfg(all(target_os = "macos", feature = "screenshot-capture"))]
fn native_capture_permission_granted() -> bool {
    // SAFETY: CoreGraphics exposes these process-global TCC queries without
    // pointer arguments or caller-owned memory.
    unsafe { CGPreflightScreenCaptureAccess() }
}

#[cfg(all(target_os = "macos", feature = "screenshot-capture"))]
fn request_native_capture_permission() -> bool {
    // SAFETY: This invokes the standard macOS Screen Recording consent prompt.
    unsafe { CGRequestScreenCaptureAccess() }
}

#[cfg(all(target_os = "macos", feature = "screenshot-capture"))]
fn native_screenshot_manager_available() -> bool {
    objc2::runtime::AnyClass::get(c"SCScreenshotManager").is_some()
}

#[cfg(all(target_os = "macos", feature = "screenshot-capture"))]
pub fn request_native_window_capture(ctx: egui::Context, window_title: String) {
    let request_id = next_request_id();
    if !native_screenshot_manager_available() {
        push_capture_event(AgentHelpCaptureEvent::Failed(AgentHelpCaptureFailure {
            request_id,
            window_title,
            kind: AgentHelpCaptureFailureKind::Unsupported,
            message: "Native full-window capture requires macOS 14 or newer. Use the normal Agent help click to capture GENtle's viewport instead."
                .to_string(),
        }));
        ctx.request_repaint();
        return;
    }
    if !native_capture_permission_granted() {
        let granted_now = request_native_capture_permission();
        push_capture_event(AgentHelpCaptureEvent::Failed(AgentHelpCaptureFailure {
            request_id,
            window_title,
            kind: if granted_now {
                AgentHelpCaptureFailureKind::RestartRequired
            } else {
                AgentHelpCaptureFailureKind::PermissionRequired
            },
            message: if granted_now {
                "macOS Screen Recording access was granted. Restart GENtle before using native window capture."
                    .to_string()
            } else {
                "Native window capture needs Screen Recording access in System Settings > Privacy & Security."
                    .to_string()
            },
        }));
        ctx.request_repaint();
        return;
    }
    let (window_id, scale_factor) = match active_macos_window_target() {
        Ok(target) => target,
        Err(message) => {
            push_capture_event(AgentHelpCaptureEvent::Failed(AgentHelpCaptureFailure {
                request_id,
                window_title,
                kind: AgentHelpCaptureFailureKind::CaptureFailed,
                message,
            }));
            ctx.request_repaint();
            return;
        }
    };
    std::thread::spawn(move || {
        let result = capture_macos_window(window_id, scale_factor).map(|image| {
            AgentHelpCaptureEvent::Captured(AgentHelpCapturedImage {
                request_id,
                window_title: window_title.clone(),
                backend: "macos.screencapturekit".to_string(),
                image,
            })
        });
        push_capture_event(result.unwrap_or_else(|message| {
            AgentHelpCaptureEvent::Failed(AgentHelpCaptureFailure {
                request_id,
                window_title,
                kind: AgentHelpCaptureFailureKind::CaptureFailed,
                message,
            })
        }));
        ctx.request_repaint();
    });
}

#[cfg(all(target_os = "macos", feature = "screenshot-capture"))]
fn capture_macos_window(window_id: u32, scale_factor: f64) -> Result<egui::ColorImage, String> {
    use screencapturekit::{
        prelude::{PixelFormat, SCContentFilter, SCShareableContent, SCStreamConfiguration},
        screenshot_manager::{CGImageExt, SCScreenshotManager},
    };

    let content = SCShareableContent::get()
        .map_err(|error| format!("Could not access shareable macOS windows: {error}"))?;
    let window = content
        .windows()
        .into_iter()
        .find(|window| window.window_id() == window_id)
        .ok_or_else(|| {
            format!("The selected GENtle window ({window_id}) is no longer available")
        })?;
    let owner_pid = window
        .owning_application()
        .map(|app| app.process_id())
        .ok_or_else(|| "The selected window has no owning application".to_string())?;
    if owner_pid != std::process::id() as i32 {
        return Err("Refusing to capture a window that is not owned by GENtle".to_string());
    }
    let frame = window.frame();
    let width = (frame.size.width * scale_factor).round().max(1.0) as u32;
    let height = (frame.size.height * scale_factor).round().max(1.0) as u32;
    let filter = SCContentFilter::create().with_window(&window).build();
    let config = SCStreamConfiguration::new()
        .with_width(width)
        .with_height(height)
        .with_pixel_format(PixelFormat::BGRA);
    let image = SCScreenshotManager::capture_image(&filter, &config)
        .map_err(|error| format!("ScreenCaptureKit could not capture the window: {error}"))?;
    let rgba = image
        .rgba_data()
        .map_err(|error| format!("Could not decode the captured window image: {error}"))?;
    let width = image.width();
    let height = image.height();
    Ok(egui::ColorImage::from_rgba_unmultiplied(
        [width, height],
        &rgba,
    ))
}

#[cfg(not(all(target_os = "macos", feature = "screenshot-capture")))]
pub fn request_native_window_capture(ctx: egui::Context, window_title: String) {
    let request_id = next_request_id();
    push_capture_event(AgentHelpCaptureEvent::Failed(AgentHelpCaptureFailure {
        request_id,
        window_title,
        kind: AgentHelpCaptureFailureKind::Unsupported,
        message: "Native full-window capture is available only in macOS builds with screenshot capture enabled."
            .to_string(),
    }));
    ctx.request_repaint();
}

pub fn open_macos_screen_recording_settings() -> Result<(), String> {
    #[cfg(target_os = "macos")]
    {
        std::process::Command::new("open")
            .arg("x-apple.systempreferences:com.apple.preference.security?Privacy_ScreenCapture")
            .status()
            .map_err(|error| format!("Could not open macOS System Settings: {error}"))
            .and_then(|status| {
                if status.success() {
                    Ok(())
                } else {
                    Err(format!(
                        "Opening macOS System Settings failed with {status}"
                    ))
                }
            })
    }
    #[cfg(not(target_os = "macos"))]
    {
        Err("Screen Recording settings are specific to macOS".to_string())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn egui_capture_event_keeps_explicit_window_identity() {
        let ctx = egui::Context::default();
        let _ = ctx.begin_pass(egui::RawInput::default());
        request_egui_viewport_capture(&ctx, "Splicing Expert");
        let output = ctx.end_pass();
        assert!(output.viewport_output.values().any(|viewport| {
            viewport.commands.iter().any(|command| {
                let egui::ViewportCommand::Screenshot(user_data) = command else {
                    return false;
                };
                user_data
                    .data
                    .as_ref()
                    .and_then(|data| data.downcast_ref::<EguiCaptureToken>())
                    .map(|token| token.window_title == "Splicing Expert")
                    .unwrap_or(false)
            })
        }));
    }

    #[test]
    fn egui_capture_command_targets_selected_viewport_once() {
        let ctx = egui::Context::default();
        let target = egui::ViewportId::from_hash_of("selected-content-window");
        let mut input = egui::RawInput {
            viewport_id: target,
            ..Default::default()
        };
        input.viewports.insert(target, Default::default());
        let _ = ctx.begin_pass(input);
        let request_id = request_egui_viewport_capture_for(&ctx, target, "Sequence map");
        let output = ctx.end_pass();

        let screenshot_tokens = output
            .viewport_output
            .get(&target)
            .into_iter()
            .flat_map(|viewport| &viewport.commands)
            .filter_map(|command| {
                let egui::ViewportCommand::Screenshot(user_data) = command else {
                    return None;
                };
                user_data
                    .data
                    .as_ref()
                    .and_then(|data| data.downcast_ref::<EguiCaptureToken>())
            })
            .collect::<Vec<_>>();
        assert_eq!(screenshot_tokens.len(), 1);
        assert_eq!(screenshot_tokens[0].request_id, request_id);
        assert_eq!(screenshot_tokens[0].window_title, "Sequence map");
    }

    #[test]
    fn egui_capture_reply_is_collected_only_once_per_request() {
        let _ = take_capture_events();
        let ctx = egui::Context::default();
        let request_id = next_request_id();
        let event = egui::Event::Screenshot {
            viewport_id: egui::ViewportId::ROOT,
            user_data: egui::UserData::new(EguiCaptureToken {
                request_id,
                window_title: "Promoter design".to_string(),
            }),
            image: std::sync::Arc::new(egui::ColorImage::filled([2, 2], egui::Color32::WHITE)),
        };
        let _ = ctx.begin_pass(egui::RawInput {
            events: vec![event],
            ..egui::RawInput::default()
        });

        collect_egui_capture_events(&ctx);
        collect_egui_capture_events(&ctx);
        let events = take_capture_events();
        let _ = ctx.end_pass();

        assert_eq!(events.len(), 1);
        let AgentHelpCaptureEvent::Captured(capture) = &events[0] else {
            panic!("expected captured image")
        };
        assert_eq!(capture.request_id, request_id);
        assert_eq!(capture.window_title, "Promoter design");
    }
}
