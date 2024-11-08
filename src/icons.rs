// Some icons used in the GUI

use eframe::egui;
use lazy_static::lazy_static;

lazy_static! {
    pub static ref SPLASH_SCREEN: egui::Image<'static> =
        egui::Image::new(egui::include_image!("../icons/GENtle.png"));
    pub static ref ICON_CIRCULAR_LINEAR: egui::Image<'static> =
        egui::Image::new(egui::include_image!("../icons/display_circular_linear.png"));
    pub static ref ICON_SHOW_SEQUENCE: egui::Image<'static> =
        egui::Image::new(egui::include_image!("../icons/show_sequence.png"));
    pub static ref ICON_SHOW_MAP: egui::Image<'static> =
        egui::Image::new(egui::include_image!("../icons/show_map.png"));
    pub static ref ICON_RESTRICTION_ENZYMES: egui::Image<'static> =
        egui::Image::new(egui::include_image!("../icons/restriction_enzymes.png"));
    pub static ref ICON_OPEN_READING_FRAMES: egui::Image<'static> =
        egui::Image::new(egui::include_image!("../icons/display_orfs.png"));
    pub static ref ICON_FEATURES: egui::Image<'static> =
        egui::Image::new(egui::include_image!("../icons/display_features.png"));
    pub static ref ICON_GC_CONTENT: egui::Image<'static> =
        egui::Image::new(egui::include_image!("../icons/gc_content.png"));
}
