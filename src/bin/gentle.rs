use eframe::egui;
use gentle::app;

fn main() {
    let options = eframe::NativeOptions {
        viewport: egui::ViewportBuilder::default().with_inner_size([320.0, 240.0]),
        ..Default::default()
    };
    let _ = eframe::run_native(
        "GENtle",
        options,
        Box::new(|_cc| Ok(Box::new(app::GENtleApp::new()))),
    );
}
