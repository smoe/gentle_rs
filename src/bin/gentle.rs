use eframe::{egui, NativeOptions};
use gentle::{about, app};
use std::{env, path::PathBuf};

#[cfg(target_os = "macos")]
fn configure_macos_process_name() {
    use objc2_foundation::{ns_string, NSProcessInfo};
    // Winit builds the macOS app menu title from NSProcessInfo::processName.
    // Set it early so the native menu shows "About GENtle" instead of "About gentle".
    unsafe {
        NSProcessInfo::processInfo().setProcessName(ns_string!("GENtle"));
    }
}

#[cfg(not(target_os = "macos"))]
fn configure_macos_process_name() {}

fn resolve_runtime_asset_path(path: &str) -> PathBuf {
    let direct = PathBuf::from(path);
    if direct.exists() {
        return direct;
    }
    if let Ok(exe_path) = env::current_exe() {
        if let Some(exe_dir) = exe_path.parent() {
            let bundled = exe_dir.join("../Resources").join(path);
            if bundled.exists() {
                return bundled;
            }
        }
    }
    PathBuf::from(path)
}

fn load_icon(path: &str) -> egui::IconData {
    let icon_path = resolve_runtime_asset_path(path);
    let (icon_rgba, icon_width, icon_height) = {
        let image = image::open(&icon_path).unwrap_or_else(|e| {
            panic!(
                "Failed to open icon path '{}': {}",
                icon_path.display(),
                e
            )
        })
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
    configure_macos_process_name();

    let args: Vec<String> = env::args().skip(1).collect();
    if args.iter().any(|a| a == "--version" || a == "-V") {
        println!("{}", about::version_cli_text());
        return Ok(());
    }
    let project_path = args.iter().find(|a| !a.starts_with('-')).cloned();

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
        Box::new(move |_cc| {
            Ok(Box::new(app::GENtleApp::new_with_project(
                project_path.as_deref(),
            )))
        }),
    )
}
