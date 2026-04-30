//! Embedded icon/resource helpers for GUI rendering.

// Some icons used in the GUI

use egui;
use std::sync::LazyLock;

pub static APP_ICON: LazyLock<egui::Image<'static>> =
    LazyLock::new(|| egui::Image::new(egui::include_image!("../../../assets/icon.png")));
pub static SPLASH_SCREEN: LazyLock<egui::Image<'static>> =
    LazyLock::new(|| egui::Image::new(egui::include_image!("../../../icons/GENtle.png")));
pub static ICON_CIRCULAR_LINEAR: LazyLock<egui::Image<'static>> = LazyLock::new(|| {
    egui::Image::new(egui::include_image!(
        "../../../icons/display_circular_linear.png"
    ))
});
pub static ICON_SHOW_SEQUENCE: LazyLock<egui::Image<'static>> =
    LazyLock::new(|| egui::Image::new(egui::include_image!("../../../icons/show_sequence.png")));
pub static ICON_SHOW_MAP: LazyLock<egui::Image<'static>> =
    LazyLock::new(|| egui::Image::new(egui::include_image!("../../../icons/show_map.png")));
pub static ICON_RESTRICTION_ENZYMES: LazyLock<egui::Image<'static>> = LazyLock::new(|| {
    egui::Image::new(egui::include_image!(
        "../../../icons/restriction_enzymes.png"
    ))
});
pub static ICON_OPEN_READING_FRAMES: LazyLock<egui::Image<'static>> =
    LazyLock::new(|| egui::Image::new(egui::include_image!("../../../icons/display_orfs.png")));
pub static ICON_FEATURES: LazyLock<egui::Image<'static>> =
    LazyLock::new(|| egui::Image::new(egui::include_image!("../../../icons/display_features.png")));
pub static ICON_GC_CONTENT: LazyLock<egui::Image<'static>> =
    LazyLock::new(|| egui::Image::new(egui::include_image!("../../../icons/gc_content.png")));
pub static ICON_METHYLATION_SITES: LazyLock<egui::Image<'static>> = LazyLock::new(|| {
    egui::Image::new(egui::include_image!(
        "../../../icons/accessories-calculator.png"
    ))
});
