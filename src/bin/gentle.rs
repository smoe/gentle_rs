//! GUI binary entry point that boots the native GENtle desktop application.

use eframe::{NativeOptions, egui};
use gentle::{about, app};
use std::{
    backtrace::Backtrace,
    env,
    fs::OpenOptions,
    io::Write,
    panic,
    path::{Path, PathBuf},
    time::{SystemTime, UNIX_EPOCH},
};

#[cfg(target_os = "macos")]
fn configure_macos_process_name() {
    use objc2_foundation::{NSProcessInfo, ns_string};
    // Winit builds the macOS app menu title from NSProcessInfo::processName.
    // Set it early so the native menu shows "About GENtle" instead of "About gentle".
    unsafe {
        NSProcessInfo::processInfo().setProcessName(ns_string!("GENtle"));
    }
}

#[cfg(not(target_os = "macos"))]
fn configure_macos_process_name() {}

fn resolve_runtime_asset_path_from(
    path: &str,
    current_exe: Option<&Path>,
    manifest_dir: Option<&Path>,
) -> PathBuf {
    let direct = PathBuf::from(path);
    if direct.exists() {
        return direct;
    }
    if let Some(manifest_dir) = manifest_dir {
        let repo_relative = manifest_dir.join(path);
        if repo_relative.exists() {
            return repo_relative;
        }
    }
    if let Some(exe_path) = current_exe {
        if let Some(exe_dir) = exe_path.parent() {
            let mut candidates: Vec<PathBuf> = Vec::new();
            if let Some(parent) = exe_dir.parent() {
                candidates.push(parent.join("Resources").join(path));
                if let Some(grandparent) = parent.parent() {
                    candidates.push(grandparent.join(path));
                }
            }
            for candidate in candidates {
                if candidate.exists() {
                    return candidate;
                }
            }
        }
    }
    PathBuf::from(path)
}

fn resolve_runtime_asset_path(path: &str) -> PathBuf {
    let current_exe = env::current_exe().ok();
    let manifest_dir = option_env!("CARGO_MANIFEST_DIR").map(Path::new);
    resolve_runtime_asset_path_from(path, current_exe.as_deref(), manifest_dir)
}

fn load_icon(path: &str) -> Option<egui::IconData> {
    let icon_path = resolve_runtime_asset_path(path);
    let (icon_rgba, icon_width, icon_height) = {
        let image = match image::open(&icon_path) {
            Ok(image) => image.into_rgba8(),
            Err(err) => {
                eprintln!(
                    "Warning: could not load app icon '{}': {}; continuing without custom icon",
                    icon_path.display(),
                    err
                );
                return None;
            }
        };
        let (width, height) = image.dimensions();
        let rgba = image.into_raw();
        (rgba, width, height)
    };

    Some(egui::IconData {
        rgba: icon_rgba,
        width: icon_width,
        height: icon_height,
    })
}

fn install_panic_logging() {
    let default_hook = panic::take_hook();
    panic::set_hook(Box::new(move |panic_info| {
        let payload = if let Some(s) = panic_info.payload().downcast_ref::<&str>() {
            (*s).to_string()
        } else if let Some(s) = panic_info.payload().downcast_ref::<String>() {
            s.clone()
        } else {
            "non-string panic payload".to_string()
        };
        let location = panic_info
            .location()
            .map(|loc| format!("{}:{}:{}", loc.file(), loc.line(), loc.column()))
            .unwrap_or_else(|| "<unknown>".to_string());
        let timestamp = SystemTime::now()
            .duration_since(UNIX_EPOCH)
            .map(|d| d.as_secs())
            .unwrap_or(0);
        let backtrace = Backtrace::force_capture();
        let message = format!(
            "[{timestamp}] panic at {location}\npayload: {payload}\nbacktrace:\n{backtrace}\n"
        );
        eprintln!("{message}");
        let mut log_path = env::temp_dir();
        log_path.push("gentle_panic.log");
        if let Ok(mut file) = OpenOptions::new().create(true).append(true).open(&log_path) {
            let _ = writeln!(file, "{message}");
        }
        default_hook(panic_info);
    }));
}

#[derive(Debug, Default)]
struct CliArgs {
    show_help: bool,
    show_version: bool,
    project_path: Option<String>,
    allow_screenshots: bool,
}

fn print_help() {
    println!(
        "Usage:\n  \
gentle [--help|-h] [--version|-V]\n  \
gentle [--project PATH] [PATH]\n\n  \
PATH can be a project file such as 'project.gentle.json'."
    );
}

fn parse_cli_args(args: &[String]) -> Result<CliArgs, String> {
    let mut parsed = CliArgs::default();
    let mut idx = 0usize;
    while idx < args.len() {
        match args[idx].as_str() {
            "--help" | "-h" => {
                parsed.show_help = true;
                idx += 1;
            }
            "--version" | "-V" => {
                parsed.show_version = true;
                idx += 1;
            }
            "--project" => {
                if idx + 1 >= args.len() {
                    return Err("Missing PATH after --project".to_string());
                }
                parsed.project_path = Some(args[idx + 1].clone());
                idx += 2;
            }
            "--allow-screenshots" => {
                return Err("--allow-screenshots is disabled by security policy".to_string());
            }
            arg if arg.starts_with('-') => {
                return Err(format!("Unknown option '{arg}'"));
            }
            path => {
                if parsed.project_path.is_some() {
                    return Err(format!(
                        "Multiple project paths provided ('{}' and '{}')",
                        parsed.project_path.as_deref().unwrap_or_default(),
                        path
                    ));
                }
                parsed.project_path = Some(path.to_string());
                idx += 1;
            }
        }
    }
    Ok(parsed)
}

fn main() -> eframe::Result<()> {
    install_panic_logging();
    configure_macos_process_name();

    let args: Vec<String> = env::args().skip(1).collect();
    let cli = match parse_cli_args(&args) {
        Ok(parsed) => parsed,
        Err(e) => {
            eprintln!("{e}");
            print_help();
            return Ok(());
        }
    };
    if cli.show_help {
        print_help();
        return Ok(());
    }
    if cli.show_version {
        println!("{}", about::version_cli_text());
        return Ok(());
    }
    if cli.allow_screenshots {
        eprintln!("--allow-screenshots is not supported in this build");
    }

    let mut viewport = egui::ViewportBuilder::default()
        .with_inner_size([800.0, 600.0])
        .with_min_inner_size([300.0, 220.0]);
    if let Some(icon) = load_icon("assets/icon.png") {
        viewport = viewport.with_icon(icon);
    }
    let options = NativeOptions {
        viewport,
        ..Default::default()
    };

    eframe::run_native(
        "GENtle",
        options,
        Box::new(move |_cc| {
            Ok(Box::new(app::GENtleApp::new_with_project(
                cli.project_path.as_deref(),
            )))
        }),
    )
}

#[cfg(test)]
mod tests {
    use super::{load_icon, resolve_runtime_asset_path_from};
    use std::path::Path;

    #[test]
    fn resolve_runtime_asset_path_uses_manifest_dir_fallback() {
        let temp = tempfile::tempdir().expect("temp dir");
        let manifest_dir = temp.path();
        let relative = "test_runtime_assets/icon.png";
        let assets_dir = manifest_dir.join("test_runtime_assets");
        std::fs::create_dir_all(&assets_dir).expect("create assets dir");
        let icon_path = assets_dir.join("icon.png");
        std::fs::write(&icon_path, b"placeholder").expect("write placeholder icon");

        let resolved = resolve_runtime_asset_path_from(relative, None, Some(manifest_dir));
        assert_eq!(resolved, icon_path);
    }

    #[test]
    fn resolve_runtime_asset_path_uses_executable_relative_dev_fallback() {
        let temp = tempfile::tempdir().expect("temp dir");
        let repo_root = temp.path();
        let relative = "test_runtime_assets/icon.png";
        let assets_dir = repo_root.join("test_runtime_assets");
        std::fs::create_dir_all(&assets_dir).expect("create assets dir");
        let icon_path = assets_dir.join("icon.png");
        std::fs::write(&icon_path, b"placeholder").expect("write placeholder icon");

        let exe_path = repo_root.join("target").join("debug").join("gentle");
        let resolved = resolve_runtime_asset_path_from(relative, Some(&exe_path), None);
        assert_eq!(resolved, icon_path);
    }

    #[test]
    fn load_icon_missing_file_returns_none() {
        assert!(load_icon("assets/does-not-exist.png").is_none());
    }

    #[test]
    fn load_icon_can_read_repo_icon_from_manifest_fallback() {
        let manifest_dir = Path::new(env!("CARGO_MANIFEST_DIR"));
        let repo_icon = manifest_dir.join("assets").join("icon.png");
        assert!(
            repo_icon.exists(),
            "expected repo icon fixture at '{}'",
            repo_icon.display()
        );
        assert!(load_icon("assets/icon.png").is_some());
    }
}
