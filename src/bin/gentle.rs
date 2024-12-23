use eframe::{egui,NativeOptions};
use gentle::app;

fn load_icon(path: &str) -> egui::IconData {
  let (icon_rgba, icon_width, icon_height) = {
      let image = image::open(path)
          .expect("Failed to open icon path")
          .into_rgba8();
      let (width, height) = image.dimensions();
      let rgba = image.into_raw();
      (rgba, width, height)
  };

  egui::IconData {
      rgba: icon_rgba,
      width: icon_width,
      height: icon_height,
  }
}

fn main() -> eframe::Result<()> {
    let options = NativeOptions {
        viewport: egui::ViewportBuilder::default()
        .with_inner_size([800.0, 600.0])
        .with_min_inner_size([300.0, 220.0])
        .with_icon(load_icon("assets/icon.png")),

        ..Default::default()
    };

    eframe::run_native(
        "GENtle",
        options,
        Box::new(|_cc| {
            Ok(Box::new(app::GENtleApp::new()))
        }),
    )
}
