//! Window backdrop configuration and rendering helpers.

use egui::{self, Color32, Rect, Ui, Vec2};
use serde::{Deserialize, Serialize};
use std::{
    collections::{HashMap, HashSet},
    env,
    path::PathBuf,
    sync::{OnceLock, RwLock},
};

const DEFAULT_MAIN_IMAGE_PATH: &str = "assets/backgrounds/Rostock_Neuer_Markt.jpeg";
const DEFAULT_SEQUENCE_IMAGE_PATH: &str = "assets/backgrounds/Lab_PCR_Machines.jpg";
const DEFAULT_SPLICING_IMAGE_PATH: &str = "assets/backgrounds/Lab_Hood.jpg";
const DEFAULT_POOL_IMAGE_PATH: &str = "assets/backgrounds/Lab_Plates_Piles.jpg";
const DEFAULT_CONFIGURATION_IMAGE_PATH: &str = "assets/backgrounds/Lab_Glas_Dry.jpg";
const DEFAULT_HELP_IMAGE_PATH: &str = "assets/backgrounds/Lab_Tubes.jpg";

const DEFAULT_MAIN_TINT_RGB: [u8; 3] = [158, 108, 66];
const DEFAULT_SEQUENCE_TINT_RGB: [u8; 3] = [176, 116, 72];
const DEFAULT_SPLICING_TINT_RGB: [u8; 3] = [113, 135, 86];
const DEFAULT_POOL_TINT_RGB: [u8; 3] = [148, 95, 58];
const DEFAULT_CONFIGURATION_TINT_RGB: [u8; 3] = [169, 118, 79];
const DEFAULT_HELP_TINT_RGB: [u8; 3] = [141, 92, 56];
const BACKDROP_TEXTURE_HINT_PX: f32 = 1024.0;
const BACKDROP_IMAGE_TINT_RGB: [u8; 3] = [255, 252, 214];
const BACKDROP_TINT_BRIGHTEN_GAMMA: f32 = 1.18;
const BACKDROP_IMAGE_OPACITY_SCALE: f32 = 0.50;
const BACKDROP_TINT_OPACITY_SCALE: f32 = 0.42;
const BACKDROP_WARM_ACCENT_BLEND: f32 = 0.28;
const BACKDROP_WARM_ACCENT_RGB: [u8; 3] = [255, 230, 156];
const BACKDROP_BASE_FILL_RGB: [u8; 3] = [247, 243, 232];
const BACKDROP_BASE_ACCENT_BLEND: f32 = 0.10;
const BACKDROP_CONTENT_VEIL_RGB: [u8; 3] = [251, 249, 242];
const BACKDROP_OVERDRAW_PX: f32 = 2.0;
pub const WINDOW_CONTENT_INSET_PX: f32 = 8.0;

fn default_main_image_path() -> String {
    DEFAULT_MAIN_IMAGE_PATH.to_string()
}

fn default_sequence_image_path() -> String {
    DEFAULT_SEQUENCE_IMAGE_PATH.to_string()
}

fn default_splicing_image_path() -> String {
    DEFAULT_SPLICING_IMAGE_PATH.to_string()
}

fn default_pool_image_path() -> String {
    DEFAULT_POOL_IMAGE_PATH.to_string()
}

fn default_configuration_image_path() -> String {
    DEFAULT_CONFIGURATION_IMAGE_PATH.to_string()
}

fn default_help_image_path() -> String {
    DEFAULT_HELP_IMAGE_PATH.to_string()
}

fn default_main_tint_rgb() -> [u8; 3] {
    DEFAULT_MAIN_TINT_RGB
}

fn default_sequence_tint_rgb() -> [u8; 3] {
    DEFAULT_SEQUENCE_TINT_RGB
}

fn default_splicing_tint_rgb() -> [u8; 3] {
    DEFAULT_SPLICING_TINT_RGB
}

fn default_pool_tint_rgb() -> [u8; 3] {
    DEFAULT_POOL_TINT_RGB
}

fn default_configuration_tint_rgb() -> [u8; 3] {
    DEFAULT_CONFIGURATION_TINT_RGB
}

fn default_help_tint_rgb() -> [u8; 3] {
    DEFAULT_HELP_TINT_RGB
}

pub fn default_tint_rgb_for_kind(kind: WindowBackdropKind) -> [u8; 3] {
    match kind {
        WindowBackdropKind::Main => DEFAULT_MAIN_TINT_RGB,
        WindowBackdropKind::Sequence => DEFAULT_SEQUENCE_TINT_RGB,
        WindowBackdropKind::Splicing => DEFAULT_SPLICING_TINT_RGB,
        WindowBackdropKind::Pool => DEFAULT_POOL_TINT_RGB,
        WindowBackdropKind::Configuration => DEFAULT_CONFIGURATION_TINT_RGB,
        WindowBackdropKind::Help => DEFAULT_HELP_TINT_RGB,
    }
}

fn sanitize_opacity(value: f32, default: f32) -> f32 {
    if value.is_finite() {
        value.clamp(0.0, 0.25)
    } else {
        default
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum WindowBackdropKind {
    Main,
    Sequence,
    Splicing,
    Pool,
    Configuration,
    Help,
}

impl WindowBackdropKind {
    fn label(self) -> &'static str {
        match self {
            Self::Main => "GENtle",
            Self::Sequence => "DNA",
            Self::Splicing => "SPLICING",
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
    #[serde(default = "default_main_image_path")]
    pub main_image_path: String,
    #[serde(default = "default_sequence_image_path")]
    pub sequence_image_path: String,
    #[serde(default = "default_splicing_image_path")]
    pub splicing_image_path: String,
    #[serde(default = "default_pool_image_path")]
    pub pool_image_path: String,
    #[serde(default = "default_configuration_image_path")]
    pub configuration_image_path: String,
    #[serde(default = "default_help_image_path")]
    pub help_image_path: String,
    #[serde(default = "default_main_tint_rgb")]
    pub main_tint_rgb: [u8; 3],
    #[serde(default = "default_sequence_tint_rgb")]
    pub sequence_tint_rgb: [u8; 3],
    #[serde(default = "default_splicing_tint_rgb")]
    pub splicing_tint_rgb: [u8; 3],
    #[serde(default = "default_pool_tint_rgb")]
    pub pool_tint_rgb: [u8; 3],
    #[serde(default = "default_configuration_tint_rgb")]
    pub configuration_tint_rgb: [u8; 3],
    #[serde(default = "default_help_tint_rgb")]
    pub help_tint_rgb: [u8; 3],
}

impl Default for WindowBackdropSettings {
    fn default() -> Self {
        Self {
            enabled: true,
            draw_images: true,
            tint_opacity: 0.16,
            image_opacity: 0.22,
            show_text_watermark: false,
            main_image_path: default_main_image_path(),
            sequence_image_path: default_sequence_image_path(),
            splicing_image_path: default_splicing_image_path(),
            pool_image_path: default_pool_image_path(),
            configuration_image_path: default_configuration_image_path(),
            help_image_path: default_help_image_path(),
            main_tint_rgb: default_main_tint_rgb(),
            sequence_tint_rgb: default_sequence_tint_rgb(),
            splicing_tint_rgb: default_splicing_tint_rgb(),
            pool_tint_rgb: default_pool_tint_rgb(),
            configuration_tint_rgb: default_configuration_tint_rgb(),
            help_tint_rgb: default_help_tint_rgb(),
        }
    }
}

impl WindowBackdropSettings {
    pub fn image_path_for_kind(&self, kind: WindowBackdropKind) -> &str {
        match kind {
            WindowBackdropKind::Main => self.main_image_path.as_str(),
            WindowBackdropKind::Sequence => self.sequence_image_path.as_str(),
            WindowBackdropKind::Splicing => self.splicing_image_path.as_str(),
            WindowBackdropKind::Pool => self.pool_image_path.as_str(),
            WindowBackdropKind::Configuration => self.configuration_image_path.as_str(),
            WindowBackdropKind::Help => self.help_image_path.as_str(),
        }
    }

    pub fn tint_color_for_kind(&self, kind: WindowBackdropKind) -> Color32 {
        let [r, g, b] = match kind {
            WindowBackdropKind::Main => self.main_tint_rgb,
            WindowBackdropKind::Sequence => self.sequence_tint_rgb,
            WindowBackdropKind::Splicing => self.splicing_tint_rgb,
            WindowBackdropKind::Pool => self.pool_tint_rgb,
            WindowBackdropKind::Configuration => self.configuration_tint_rgb,
            WindowBackdropKind::Help => self.help_tint_rgb,
        };
        Color32::from_rgb(r, g, b)
    }

    pub fn apply_runtime_defaults_if_legacy(&mut self) {
        self.main_image_path = migrate_legacy_image_path(&self.main_image_path);
        self.sequence_image_path = migrate_legacy_image_path(&self.sequence_image_path);
        self.splicing_image_path = migrate_legacy_image_path(&self.splicing_image_path);
        self.pool_image_path = migrate_legacy_image_path(&self.pool_image_path);
        self.configuration_image_path = migrate_legacy_image_path(&self.configuration_image_path);
        self.help_image_path = migrate_legacy_image_path(&self.help_image_path);

        self.tint_opacity = sanitize_opacity(self.tint_opacity, 0.16);
        self.image_opacity = sanitize_opacity(self.image_opacity, 0.22);

        let legacy_rows_empty = self.main_image_path.trim().is_empty()
            && self.sequence_image_path.trim().is_empty()
            && self.pool_image_path.trim().is_empty()
            && self.configuration_image_path.trim().is_empty()
            && self.help_image_path.trim().is_empty();
        if legacy_rows_empty {
            self.main_image_path = default_main_image_path();
            self.sequence_image_path = default_sequence_image_path();
            if self.splicing_image_path.trim().is_empty() {
                self.splicing_image_path = default_splicing_image_path();
            }
            self.pool_image_path = default_pool_image_path();
            self.configuration_image_path = default_configuration_image_path();
            self.help_image_path = default_help_image_path();
        }

        if legacy_rows_empty && !self.enabled {
            self.enabled = true;
        }
    }
}

static ACTIVE_SETTINGS: OnceLock<RwLock<WindowBackdropSettings>> = OnceLock::new();
static RESOLVED_URI_CACHE: OnceLock<RwLock<HashMap<String, Option<String>>>> = OnceLock::new();
static PREWARMED_URI_SET: OnceLock<RwLock<HashSet<String>>> = OnceLock::new();

fn active_settings() -> &'static RwLock<WindowBackdropSettings> {
    ACTIVE_SETTINGS.get_or_init(|| RwLock::new(WindowBackdropSettings::default()))
}

fn resolved_uri_cache() -> &'static RwLock<HashMap<String, Option<String>>> {
    RESOLVED_URI_CACHE.get_or_init(|| RwLock::new(HashMap::new()))
}

fn prewarmed_uri_set() -> &'static RwLock<HashSet<String>> {
    PREWARMED_URI_SET.get_or_init(|| RwLock::new(HashSet::new()))
}

fn clear_runtime_caches() {
    if let Ok(mut cache) = resolved_uri_cache().write() {
        cache.clear();
    }
    if let Ok(mut warmed) = prewarmed_uri_set().write() {
        warmed.clear();
    }
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
    clear_runtime_caches();
}

fn alpha_to_u8(opacity: f32) -> u8 {
    (opacity.clamp(0.0, 1.0) * 255.0).round() as u8
}

fn kind_image_opacity_scale(kind: WindowBackdropKind) -> f32 {
    match kind {
        WindowBackdropKind::Main => 0.86,
        WindowBackdropKind::Sequence => 0.56,
        WindowBackdropKind::Splicing => 0.28,
        WindowBackdropKind::Pool => 0.35,
        WindowBackdropKind::Configuration => 0.48,
        WindowBackdropKind::Help => 0.44,
    }
}

fn kind_tint_opacity_scale(kind: WindowBackdropKind) -> f32 {
    match kind {
        WindowBackdropKind::Main => 0.64,
        WindowBackdropKind::Sequence => 0.48,
        WindowBackdropKind::Splicing => 0.34,
        WindowBackdropKind::Pool => 0.38,
        WindowBackdropKind::Configuration => 0.50,
        WindowBackdropKind::Help => 0.46,
    }
}

fn kind_content_veil_opacity(kind: WindowBackdropKind) -> f32 {
    match kind {
        WindowBackdropKind::Main => 0.24,
        WindowBackdropKind::Sequence => 0.40,
        WindowBackdropKind::Splicing => 0.58,
        WindowBackdropKind::Pool => 0.50,
        WindowBackdropKind::Configuration => 0.40,
        WindowBackdropKind::Help => 0.38,
    }
}

fn backdrop_texture_hint() -> Vec2 {
    Vec2::splat(BACKDROP_TEXTURE_HINT_PX)
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

fn resolve_runtime_asset_uri_uncached(path: &str) -> Option<String> {
    let trimmed = path.trim();
    if trimmed.is_empty() {
        return None;
    }

    // egui file loading expects URI schemes such as file://
    if trimmed.starts_with("file://")
        || trimmed.starts_with("http://")
        || trimmed.starts_with("https://")
        || trimmed.starts_with("bytes://")
    {
        return Some(trimmed.to_string());
    }

    let direct = PathBuf::from(trimmed);
    if direct.exists() {
        return Some(format!("file://{}", direct.to_string_lossy()));
    }

    if let Ok(exe_path) = env::current_exe() {
        if let Some(exe_dir) = exe_path.parent() {
            let bundled = exe_dir.join("../Resources").join(trimmed);
            if bundled.exists() {
                return Some(format!("file://{}", bundled.to_string_lossy()));
            }
        }
    }
    None
}

fn resolve_runtime_asset_uri(path: &str) -> Option<String> {
    let trimmed = path.trim();
    if trimmed.is_empty() {
        return None;
    }

    if trimmed.starts_with("file://")
        || trimmed.starts_with("http://")
        || trimmed.starts_with("https://")
        || trimmed.starts_with("bytes://")
    {
        return Some(trimmed.to_string());
    }

    if let Ok(cache) = resolved_uri_cache().read() {
        if let Some(cached) = cache.get(trimmed) {
            return cached.clone();
        }
    }

    let resolved = resolve_runtime_asset_uri_uncached(trimmed);
    if let Ok(mut cache) = resolved_uri_cache().write() {
        cache.insert(trimmed.to_string(), resolved.clone());
    }
    resolved
}

pub fn preload_window_backdrop_images(ctx: &egui::Context, settings: &WindowBackdropSettings) {
    if !settings.enabled || !settings.draw_images || alpha_to_u8(settings.image_opacity) == 0 {
        return;
    }

    for kind in [
        WindowBackdropKind::Main,
        WindowBackdropKind::Sequence,
        WindowBackdropKind::Splicing,
        WindowBackdropKind::Pool,
        WindowBackdropKind::Configuration,
        WindowBackdropKind::Help,
    ] {
        let Some(resolved_uri) = resolve_runtime_asset_uri(settings.image_path_for_kind(kind))
        else {
            continue;
        };

        let already_warmed = prewarmed_uri_set()
            .read()
            .map(|set| set.contains(&resolved_uri))
            .unwrap_or(false);
        if already_warmed {
            continue;
        }

        match egui::Image::new(resolved_uri.clone())
            .show_loading_spinner(false)
            // Keep one shared texture-size hint across windows to avoid per-window-size
            // cache misses and repeated image decode/resample work.
            .load_for_size(ctx, backdrop_texture_hint())
        {
            Ok(egui::load::TexturePoll::Ready { .. })
            | Ok(egui::load::TexturePoll::Pending { .. }) => {
                if let Ok(mut set) = prewarmed_uri_set().write() {
                    set.insert(resolved_uri);
                }
            }
            Err(_) => {}
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct WindowBackdropPreloadStatus {
    pub active: bool,
    pub configured_images: usize,
    pub queued_images: usize,
}

pub fn window_backdrop_preload_status(
    settings: &WindowBackdropSettings,
) -> WindowBackdropPreloadStatus {
    if !settings.enabled || !settings.draw_images || alpha_to_u8(settings.image_opacity) == 0 {
        return WindowBackdropPreloadStatus {
            active: false,
            configured_images: 0,
            queued_images: 0,
        };
    }

    let mut configured_uris = HashSet::new();
    for kind in [
        WindowBackdropKind::Main,
        WindowBackdropKind::Sequence,
        WindowBackdropKind::Splicing,
        WindowBackdropKind::Pool,
        WindowBackdropKind::Configuration,
        WindowBackdropKind::Help,
    ] {
        if let Some(uri) = resolve_runtime_asset_uri(settings.image_path_for_kind(kind)) {
            configured_uris.insert(uri);
        }
    }

    let queued_images = prewarmed_uri_set()
        .read()
        .map(|set| {
            configured_uris
                .iter()
                .filter(|uri| set.contains(*uri))
                .count()
        })
        .unwrap_or(0);

    WindowBackdropPreloadStatus {
        active: true,
        configured_images: configured_uris.len(),
        queued_images,
    }
}

fn cover_uv_rect(texture_size: Vec2, target_size: Vec2) -> Rect {
    if texture_size.x <= 0.0
        || texture_size.y <= 0.0
        || target_size.x <= 0.0
        || target_size.y <= 0.0
    {
        return Rect::from_min_max(egui::pos2(0.0, 0.0), egui::pos2(1.0, 1.0));
    }

    let image_aspect = texture_size.x / texture_size.y;
    let target_aspect = target_size.x / target_size.y;
    if !image_aspect.is_finite() || !target_aspect.is_finite() {
        return Rect::from_min_max(egui::pos2(0.0, 0.0), egui::pos2(1.0, 1.0));
    }

    if image_aspect > target_aspect {
        // Wider image: crop horizontal edges.
        let visible_fraction = (target_aspect / image_aspect).clamp(0.0, 1.0);
        let margin = (1.0 - visible_fraction) * 0.5;
        Rect::from_min_max(egui::pos2(margin, 0.0), egui::pos2(1.0 - margin, 1.0))
    } else if image_aspect < target_aspect {
        // Taller image: crop vertical edges.
        let visible_fraction = (image_aspect / target_aspect).clamp(0.0, 1.0);
        let margin = (1.0 - visible_fraction) * 0.5;
        Rect::from_min_max(egui::pos2(0.0, margin), egui::pos2(1.0, 1.0 - margin))
    } else {
        Rect::from_min_max(egui::pos2(0.0, 0.0), egui::pos2(1.0, 1.0))
    }
}

pub fn validate_window_backdrop_image_path(path: &str) -> Result<String, String> {
    let trimmed = path.trim();
    if trimmed.is_empty() {
        return Err("No image path configured".to_string());
    }
    resolve_runtime_asset_uri(trimmed).ok_or_else(|| {
        "Path not found in filesystem or app bundle Resources/<path> lookup.".to_string()
    })
}

pub fn paint_window_backdrop(
    ui: &mut Ui,
    kind: WindowBackdropKind,
    settings: &WindowBackdropSettings,
) {
    if !settings.enabled {
        return;
    }

    // Overdraw slightly to avoid thin panel-edge seams from backend sampling/frame rounding.
    let rect = ui.max_rect().expand(BACKDROP_OVERDRAW_PX);
    // Keep a warm, low-contrast palette even when backend color management changes.
    let base_accent = settings.tint_color_for_kind(kind);
    let blend = BACKDROP_WARM_ACCENT_BLEND.clamp(0.0, 1.0);
    let mixed_r = ((base_accent.r() as f32 * (1.0 - blend))
        + (BACKDROP_WARM_ACCENT_RGB[0] as f32 * blend))
        .round()
        .clamp(0.0, 255.0) as u8;
    let mixed_g = ((base_accent.g() as f32 * (1.0 - blend))
        + (BACKDROP_WARM_ACCENT_RGB[1] as f32 * blend))
        .round()
        .clamp(0.0, 255.0) as u8;
    let mixed_b = ((base_accent.b() as f32 * (1.0 - blend))
        + (BACKDROP_WARM_ACCENT_RGB[2] as f32 * blend))
        .round()
        .clamp(0.0, 255.0) as u8;
    let accent =
        Color32::from_rgb(mixed_r, mixed_g, mixed_b).gamma_multiply(BACKDROP_TINT_BRIGHTEN_GAMMA);
    let tint_alpha = alpha_to_u8(
        settings.tint_opacity * BACKDROP_TINT_OPACITY_SCALE * kind_tint_opacity_scale(kind),
    );
    let image_alpha = alpha_to_u8(
        settings.image_opacity * BACKDROP_IMAGE_OPACITY_SCALE * kind_image_opacity_scale(kind),
    );
    let painter = ui.painter_at(rect);

    // Provide a deterministic light base so image/tint do not blend against a dark panel clear color.
    let base_blend = BACKDROP_BASE_ACCENT_BLEND.clamp(0.0, 1.0);
    let base_r = ((BACKDROP_BASE_FILL_RGB[0] as f32 * (1.0 - base_blend))
        + (accent.r() as f32 * base_blend))
        .round()
        .clamp(0.0, 255.0) as u8;
    let base_g = ((BACKDROP_BASE_FILL_RGB[1] as f32 * (1.0 - base_blend))
        + (accent.g() as f32 * base_blend))
        .round()
        .clamp(0.0, 255.0) as u8;
    let base_b = ((BACKDROP_BASE_FILL_RGB[2] as f32 * (1.0 - base_blend))
        + (accent.b() as f32 * base_blend))
        .round()
        .clamp(0.0, 255.0) as u8;
    painter.rect_filled(rect, 0.0, Color32::from_rgb(base_r, base_g, base_b));

    if settings.draw_images && image_alpha > 0 {
        let path = settings.image_path_for_kind(kind).trim();
        if let Some(resolved_uri) = resolve_runtime_asset_uri(path) {
            let image_tint = Color32::from_rgba_unmultiplied(
                BACKDROP_IMAGE_TINT_RGB[0],
                BACKDROP_IMAGE_TINT_RGB[1],
                BACKDROP_IMAGE_TINT_RGB[2],
                image_alpha,
            );
            match egui::Image::new(resolved_uri)
                .show_loading_spinner(false)
                // Use a stable texture-size hint instead of current viewport dimensions,
                // so opening differently sized windows does not trigger extra reload paths.
                .load_for_size(ui.ctx(), backdrop_texture_hint())
            {
                Ok(egui::load::TexturePoll::Ready { texture }) => {
                    let uv = cover_uv_rect(texture.size, rect.size());
                    painter.image(texture.id, rect, uv, image_tint);
                }
                Ok(egui::load::TexturePoll::Pending { .. }) => {
                    // Keep repainting while the backdrop is loading, but do not show
                    // a spinner arc for full-window background images.
                    ui.ctx().request_repaint();
                }
                Err(_) => {
                    // Keep tint/watermark fallback only when image loading fails.
                }
            }
        }
    }

    if tint_alpha > 0 {
        painter.rect_filled(
            rect,
            0.0,
            Color32::from_rgba_unmultiplied(accent.r(), accent.g(), accent.b(), tint_alpha),
        );
    }

    let veil_alpha = alpha_to_u8(kind_content_veil_opacity(kind));
    if veil_alpha > 0 {
        painter.rect_filled(
            rect,
            0.0,
            Color32::from_rgba_unmultiplied(
                BACKDROP_CONTENT_VEIL_RGB[0],
                BACKDROP_CONTENT_VEIL_RGB[1],
                BACKDROP_CONTENT_VEIL_RGB[2],
                veil_alpha,
            ),
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

pub fn with_window_content_inset<R>(ui: &mut Ui, add_contents: impl FnOnce(&mut Ui) -> R) -> R {
    let outer_rect =
        resolve_window_content_outer_rect(ui.available_rect_before_wrap(), ui.max_rect());
    let inset_limit = (outer_rect.width().min(outer_rect.height()) * 0.25).max(0.0);
    let inset = WINDOW_CONTENT_INSET_PX.min(inset_limit);
    let inner_rect = outer_rect.shrink2(Vec2::splat(inset));
    let target_rect = if inner_rect.width() >= 1.0 && inner_rect.height() >= 1.0 {
        inner_rect
    } else {
        outer_rect
    };
    let builder = egui::UiBuilder::new()
        .max_rect(target_rect)
        .layout(*ui.layout());
    let mut child = ui.new_child(builder);
    add_contents(&mut child)
}

fn resolve_window_content_outer_rect(
    available_rect: egui::Rect,
    max_rect: egui::Rect,
) -> egui::Rect {
    let overlap = available_rect.intersect(max_rect);
    if overlap.width() >= 1.0 && overlap.height() >= 1.0 {
        overlap
    } else if available_rect.width() >= 1.0 && available_rect.height() >= 1.0 {
        available_rect
    } else {
        max_rect
    }
}

#[cfg(test)]
mod tests {
    use super::{
        DEFAULT_SEQUENCE_IMAGE_PATH, DEFAULT_SPLICING_IMAGE_PATH, WindowBackdropKind,
        WindowBackdropSettings, cover_uv_rect, default_tint_rgb_for_kind,
        kind_content_veil_opacity, kind_image_opacity_scale, kind_tint_opacity_scale,
        resolve_runtime_asset_uri, resolve_window_content_outer_rect,
    };
    use egui::{Rect, Vec2, pos2};

    #[test]
    fn resolve_runtime_asset_uri_preserves_supported_uri_schemes() {
        assert_eq!(
            resolve_runtime_asset_uri("file://assets/backgrounds/Lab_Tubes.jpg").as_deref(),
            Some("file://assets/backgrounds/Lab_Tubes.jpg")
        );
        assert_eq!(
            resolve_runtime_asset_uri("https://example.org/bg.jpg").as_deref(),
            Some("https://example.org/bg.jpg")
        );
        assert_eq!(
            resolve_runtime_asset_uri("bytes://bg.png").as_deref(),
            Some("bytes://bg.png")
        );
    }

    #[test]
    fn resolve_runtime_asset_uri_converts_existing_file_to_file_uri() {
        let tmp = tempfile::NamedTempFile::new().expect("temp file");
        let path = tmp.path().to_string_lossy().to_string();
        let uri = resolve_runtime_asset_uri(&path).expect("resolved file uri");
        assert!(uri.starts_with("file://"));
        assert!(uri.ends_with(&path), "uri='{uri}' path='{path}'");
    }

    #[test]
    fn cover_uv_rect_crops_wider_images_horizontally() {
        let uv = cover_uv_rect(Vec2::new(2000.0, 1000.0), Vec2::new(1000.0, 1000.0));
        assert!(uv.min.x > 0.0);
        assert!(uv.max.x < 1.0);
        assert_eq!(uv.min.y, 0.0);
        assert_eq!(uv.max.y, 1.0);
    }

    #[test]
    fn cover_uv_rect_crops_taller_images_vertically() {
        let uv = cover_uv_rect(Vec2::new(1000.0, 2000.0), Vec2::new(1000.0, 1000.0));
        assert!(uv.min.y > 0.0);
        assert!(uv.max.y < 1.0);
        assert_eq!(uv.min.x, 0.0);
        assert_eq!(uv.max.x, 1.0);
    }

    #[test]
    fn default_settings_include_splicing_backdrop_path() {
        let settings = WindowBackdropSettings::default();
        assert_eq!(
            settings.image_path_for_kind(WindowBackdropKind::Splicing),
            DEFAULT_SPLICING_IMAGE_PATH
        );
    }

    #[test]
    fn default_settings_include_distinct_splicing_tint() {
        let settings = WindowBackdropSettings::default();
        let sequence = settings.tint_color_for_kind(WindowBackdropKind::Sequence);
        let splicing = settings.tint_color_for_kind(WindowBackdropKind::Splicing);
        assert_ne!(sequence, splicing);
    }

    #[test]
    fn serde_defaults_fill_missing_image_and_color_fields() {
        let legacy_json = r#"{
            "enabled": true,
            "draw_images": true,
            "tint_opacity": 0.16,
            "image_opacity": 0.22,
            "show_text_watermark": false,
            "main_image_path": "assets/backgrounds/legacy_main.jpg"
        }"#;
        let parsed: WindowBackdropSettings =
            serde_json::from_str(legacy_json).expect("deserialize legacy backdrop settings");

        assert_eq!(parsed.main_image_path, "assets/backgrounds/legacy_main.jpg");
        assert_eq!(parsed.sequence_image_path, DEFAULT_SEQUENCE_IMAGE_PATH);
        assert_eq!(parsed.splicing_image_path, DEFAULT_SPLICING_IMAGE_PATH);
        let main_tint = parsed.tint_color_for_kind(WindowBackdropKind::Main);
        let sequence_tint = parsed.tint_color_for_kind(WindowBackdropKind::Sequence);
        assert_ne!(main_tint, egui::Color32::BLACK);
        assert_ne!(sequence_tint, egui::Color32::BLACK);
    }

    #[test]
    fn default_tint_accessor_matches_defaults() {
        let settings = WindowBackdropSettings::default();
        for kind in [
            WindowBackdropKind::Main,
            WindowBackdropKind::Sequence,
            WindowBackdropKind::Splicing,
            WindowBackdropKind::Pool,
            WindowBackdropKind::Configuration,
            WindowBackdropKind::Help,
        ] {
            let expected = default_tint_rgb_for_kind(kind);
            let actual = settings.tint_color_for_kind(kind);
            assert_eq!(
                actual,
                egui::Color32::from_rgb(expected[0], expected[1], expected[2])
            );
        }
    }

    #[test]
    fn apply_runtime_defaults_clamps_invalid_opacity_values() {
        let mut settings = WindowBackdropSettings {
            tint_opacity: f32::NAN,
            image_opacity: 10.0,
            ..WindowBackdropSettings::default()
        };
        settings.apply_runtime_defaults_if_legacy();
        assert_eq!(settings.tint_opacity, 0.16);
        assert_eq!(settings.image_opacity, 0.25);

        settings.tint_opacity = -1.0;
        settings.image_opacity = -0.2;
        settings.apply_runtime_defaults_if_legacy();
        assert_eq!(settings.tint_opacity, 0.0);
        assert_eq!(settings.image_opacity, 0.0);
    }

    #[test]
    fn sequence_backdrop_uses_stronger_image_and_tint_scales() {
        assert!(kind_image_opacity_scale(WindowBackdropKind::Sequence) >= 0.50);
        assert!(kind_tint_opacity_scale(WindowBackdropKind::Sequence) >= 0.45);
        assert!(kind_content_veil_opacity(WindowBackdropKind::Sequence) <= 0.42);
    }

    #[test]
    fn resolve_window_content_outer_rect_prefers_available_rect_inside_window() {
        let max_rect = Rect::from_min_max(pos2(0.0, 0.0), pos2(1000.0, 800.0));
        let available_rect = Rect::from_min_max(pos2(40.0, 60.0), pos2(620.0, 520.0));
        assert_eq!(
            resolve_window_content_outer_rect(available_rect, max_rect),
            available_rect
        );
    }

    #[test]
    fn resolve_window_content_outer_rect_falls_back_to_max_rect_when_available_is_invalid() {
        let max_rect = Rect::from_min_max(pos2(0.0, 0.0), pos2(1000.0, 800.0));
        let available_rect = Rect::from_min_max(pos2(40.0, 60.0), pos2(40.0, 60.0));
        assert_eq!(
            resolve_window_content_outer_rect(available_rect, max_rect),
            max_rect
        );
    }
}
