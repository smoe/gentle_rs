use eframe::egui::{self, Color32, Ui};
use serde::{Deserialize, Serialize};
use std::{
    env,
    path::PathBuf,
    sync::{OnceLock, RwLock},
};

const DEFAULT_MAIN_IMAGE_PATH: &str = "assets/backgrounds/Rostock_Neuer_Markt.jpeg";
const DEFAULT_SEQUENCE_IMAGE_PATH: &str = "assets/backgrounds/Lab_PCR_Machines.jpg";
const DEFAULT_POOL_IMAGE_PATH: &str = "assets/backgrounds/Lab_Plates_Piles.jpg";
const DEFAULT_CONFIGURATION_IMAGE_PATH: &str = "assets/backgrounds/Lab_Glas_Dry.jpg";
const DEFAULT_HELP_IMAGE_PATH: &str = "assets/backgrounds/Lab_Tubes.jpg";

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
            enabled: true,
            draw_images: true,
            tint_opacity: 0.16,
            image_opacity: 0.22,
            show_text_watermark: false,
            main_image_path: DEFAULT_MAIN_IMAGE_PATH.to_string(),
            sequence_image_path: DEFAULT_SEQUENCE_IMAGE_PATH.to_string(),
            pool_image_path: DEFAULT_POOL_IMAGE_PATH.to_string(),
            configuration_image_path: DEFAULT_CONFIGURATION_IMAGE_PATH.to_string(),
            help_image_path: DEFAULT_HELP_IMAGE_PATH.to_string(),
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

    pub fn apply_runtime_defaults_if_legacy(&mut self) {
        self.main_image_path = migrate_legacy_image_path(&self.main_image_path);
        self.sequence_image_path = migrate_legacy_image_path(&self.sequence_image_path);
        self.pool_image_path = migrate_legacy_image_path(&self.pool_image_path);
        self.configuration_image_path = migrate_legacy_image_path(&self.configuration_image_path);
        self.help_image_path = migrate_legacy_image_path(&self.help_image_path);

        let all_empty = self.main_image_path.trim().is_empty()
            && self.sequence_image_path.trim().is_empty()
            && self.pool_image_path.trim().is_empty()
            && self.configuration_image_path.trim().is_empty()
            && self.help_image_path.trim().is_empty();
        if all_empty {
            self.main_image_path = DEFAULT_MAIN_IMAGE_PATH.to_string();
            self.sequence_image_path = DEFAULT_SEQUENCE_IMAGE_PATH.to_string();
            self.pool_image_path = DEFAULT_POOL_IMAGE_PATH.to_string();
            self.configuration_image_path = DEFAULT_CONFIGURATION_IMAGE_PATH.to_string();
            self.help_image_path = DEFAULT_HELP_IMAGE_PATH.to_string();
        }

        if all_empty && !self.enabled {
            self.enabled = true;
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
        WindowBackdropKind::Main => Color32::from_rgb(158, 108, 66),
        WindowBackdropKind::Sequence => Color32::from_rgb(176, 116, 72),
        WindowBackdropKind::Pool => Color32::from_rgb(148, 95, 58),
        WindowBackdropKind::Configuration => Color32::from_rgb(169, 118, 79),
        WindowBackdropKind::Help => Color32::from_rgb(141, 92, 56),
    }
}

fn alpha_to_u8(opacity: f32) -> u8 {
    (opacity.clamp(0.0, 1.0) * 255.0).round() as u8
}

fn migrate_legacy_image_path(path: &str) -> String {
    let trimmed = path.trim();
    if trimmed.is_empty() {
        return String::new();
    }
    if PathBuf::from(trimmed).exists() {
        return trimmed.to_string();
    }
    if let Some(suffix) = trimmed.strip_prefix("docs/backgrounds/") {
        let candidate = format!("assets/backgrounds/{suffix}");
        if PathBuf::from(&candidate).exists() {
            return candidate;
        }
    }
    trimmed.to_string()
}

fn resolve_runtime_asset_path(path: &str) -> Option<String> {
    let trimmed = path.trim();
    if trimmed.is_empty() {
        return None;
    }

    let direct = PathBuf::from(trimmed);
    if direct.exists() {
        return Some(direct.to_string_lossy().to_string());
    }

    if let Ok(exe_path) = env::current_exe() {
        if let Some(exe_dir) = exe_path.parent() {
            let bundled = exe_dir.join("../Resources").join(trimmed);
            if bundled.exists() {
                return Some(bundled.to_string_lossy().to_string());
            }
        }
    }
    None
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

    if settings.draw_images {
        let path = settings.image_path_for_kind(kind).trim();
        if let Some(resolved_path) = resolve_runtime_asset_path(path) {
            let image = egui::Image::new(resolved_path)
                .fit_to_exact_size(rect.size())
                .tint(Color32::from_rgba_unmultiplied(
                    196,
                    196,
                    196,
                    image_alpha,
                ))
                .sense(egui::Sense::hover());
            let _ = ui.put(rect, image);
        }
    }

    if tint_alpha > 0 {
        painter.rect_filled(
            rect,
            0.0,
            Color32::from_rgba_unmultiplied(accent.r(), accent.g(), accent.b(), tint_alpha),
        );
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
