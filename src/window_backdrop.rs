use eframe::egui::{self, Color32, Ui};
use serde::{Deserialize, Serialize};
use std::{
    path::Path,
    sync::{OnceLock, RwLock},
};

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum WindowBackdropKind {
    Main,
    Sequence,
    Pool,
    Configuration,
    Help,
}

impl WindowBackdropKind {
    fn label(self) -> &'static str {
        match self {
            Self::Main => "GENtle",
            Self::Sequence => "DNA",
            Self::Pool => "POOL",
            Self::Configuration => "SETTINGS",
            Self::Help => "HELP",
        }
    }
}

#[derive(Clone, Debug, Serialize, Deserialize)]
#[serde(default)]
pub struct WindowBackdropSettings {
    pub enabled: bool,
    pub draw_images: bool,
    pub tint_opacity: f32,
    pub image_opacity: f32,
    pub show_text_watermark: bool,
    pub main_image_path: String,
    pub sequence_image_path: String,
    pub pool_image_path: String,
    pub configuration_image_path: String,
    pub help_image_path: String,
}

impl Default for WindowBackdropSettings {
    fn default() -> Self {
        Self {
            enabled: false,
            draw_images: true,
            tint_opacity: 0.06,
            image_opacity: 0.10,
            show_text_watermark: true,
            main_image_path: String::new(),
            sequence_image_path: String::new(),
            pool_image_path: String::new(),
            configuration_image_path: String::new(),
            help_image_path: String::new(),
        }
    }
}

impl WindowBackdropSettings {
    pub fn image_path_for_kind(&self, kind: WindowBackdropKind) -> &str {
        match kind {
            WindowBackdropKind::Main => self.main_image_path.as_str(),
            WindowBackdropKind::Sequence => self.sequence_image_path.as_str(),
            WindowBackdropKind::Pool => self.pool_image_path.as_str(),
            WindowBackdropKind::Configuration => self.configuration_image_path.as_str(),
            WindowBackdropKind::Help => self.help_image_path.as_str(),
        }
    }
}

static ACTIVE_SETTINGS: OnceLock<RwLock<WindowBackdropSettings>> = OnceLock::new();

fn active_settings() -> &'static RwLock<WindowBackdropSettings> {
    ACTIVE_SETTINGS.get_or_init(|| RwLock::new(WindowBackdropSettings::default()))
}

pub fn current_window_backdrop_settings() -> WindowBackdropSettings {
    active_settings()
        .read()
        .map(|v| v.clone())
        .unwrap_or_default()
}

pub fn set_window_backdrop_settings(settings: WindowBackdropSettings) {
    if let Ok(mut slot) = active_settings().write() {
        *slot = settings;
    }
}

fn kind_color(kind: WindowBackdropKind) -> Color32 {
    match kind {
        WindowBackdropKind::Main => Color32::from_rgb(28, 80, 120),
        WindowBackdropKind::Sequence => Color32::from_rgb(26, 118, 92),
        WindowBackdropKind::Pool => Color32::from_rgb(124, 82, 34),
        WindowBackdropKind::Configuration => Color32::from_rgb(90, 74, 118),
        WindowBackdropKind::Help => Color32::from_rgb(44, 96, 112),
    }
}

fn alpha_to_u8(opacity: f32) -> u8 {
    (opacity.clamp(0.0, 1.0) * 255.0).round() as u8
}

pub fn paint_window_backdrop(
    ui: &mut Ui,
    kind: WindowBackdropKind,
    settings: &WindowBackdropSettings,
) {
    if !settings.enabled {
        return;
    }

    let rect = ui.max_rect();
    let accent = kind_color(kind);
    let tint_alpha = alpha_to_u8(settings.tint_opacity);
    let image_alpha = alpha_to_u8(settings.image_opacity);
    let painter = ui.painter_at(rect);

    if tint_alpha > 0 {
        painter.rect_filled(
            rect,
            0.0,
            Color32::from_rgba_unmultiplied(accent.r(), accent.g(), accent.b(), tint_alpha),
        );
    }

    if settings.draw_images {
        let path = settings.image_path_for_kind(kind).trim();
        if !path.is_empty() && Path::new(path).exists() {
            let image = egui::Image::new(path.to_string())
                .fit_to_exact_size(rect.size())
                .tint(Color32::from_rgba_unmultiplied(
                    accent.r(),
                    accent.g(),
                    accent.b(),
                    image_alpha,
                ))
                .sense(egui::Sense::hover());
            let _ = ui.put(rect, image);
        }
    }

    if settings.show_text_watermark {
        let text_color =
            Color32::from_rgba_unmultiplied(accent.r(), accent.g(), accent.b(), alpha_to_u8(0.08));
        painter.text(
            rect.center(),
            egui::Align2::CENTER_CENTER,
            kind.label(),
            egui::FontId::proportional((rect.height() * 0.12).clamp(32.0, 128.0)),
            text_color,
        );
    }
}
