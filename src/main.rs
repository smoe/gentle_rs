use crate::facility::Facility;
use lazy_static::lazy_static;

pub mod amino_acids;
pub mod app;
pub mod dna_sequence;
pub mod enzymes;
pub mod facility;
pub mod icons;
pub mod left_panel;
pub mod main_area;
pub mod main_area_dna;
pub mod protease;
pub mod render_dna_circular;
pub mod render_dna_linear;
pub mod restriction_enzyme;

lazy_static! {
    // IUPAC codes etc.
    pub static ref FACILITY: Facility = Facility::new();
}

fn main() {
    let native_options = eframe::NativeOptions {
        persist_window: true,
        // resizable: true,
        // initial_window_size: Some(egui::Vec2 { x: 300.0, y: 300.0 }),
        // min_window_size: Some(egui::Vec2 { x: 300.0, y: 300.0 }),
        // icon_data: Some(eframe::IconData {
        //     rgba: ICON.to_vec(),
        //     width: 32,
        //     height: 32,
        // }),
        ..Default::default()
    };

    eframe::run_native(
        "GENtle",
        native_options,
        Box::new(|cc| Ok(Box::new(app::GENtleApp::new(cc)))),
    )
    .unwrap();
}
