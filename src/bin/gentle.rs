//! GUI binary entry point that boots the native GENtle desktop application.

use eframe::{NativeOptions, egui};
use gentle::{
    about,
    agent_bridge::{
        AGENT_BASE_URL_ENV, AGENT_MODEL_ENV, ANTHROPIC_API_KEY_ENV, MISTRAL_API_KEY_ENV,
        OPENAI_API_KEY_ENV,
    },
    app,
    cli_support::{SingleProjectCliOptions, parse_single_project_cli_args},
    gui_profiler,
};
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
    if let Some(exe_path) = current_exe
        && let Some(exe_dir) = exe_path.parent()
    {
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

fn help_text() -> String {
    format!(
        "Usage:\n  \
gentle [--help|-h] [--version|-V]\n  \
gentle [--project PATH] [PATH]\n\n  \
PATH can be a project file such as 'project.gentle.json'.\n\n\
Related GENtle executables:\n  \
gentle_cli              automation, JSON operations/workflows, and shared shell commands\n  \
gentle_mcp              MCP stdio server for external agents and tool clients\n  \
gentle_js               JavaScript shell, when built with js-interface\n  \
gentle_lua              Lua shell, when built with lua-interface\n  \
gentle_examples_docs    regenerate examples and tutorial artifacts\n\n\
If you expected a terminal command surface, try: gentle_cli --help\n\n\
Agent Assistant API environment:\n  \
{OPENAI_API_KEY_ENV}       OpenAI Platform API key for native OpenAI transport\n  \
{ANTHROPIC_API_KEY_ENV}    Anthropic Console API key for native Claude transport\n  \
{MISTRAL_API_KEY_ENV}      Mistral La Plateforme API key for native Mistral transport\n  \
{AGENT_BASE_URL_ENV}  optional native/OpenAI-compatible base URL override\n  \
{AGENT_MODEL_ENV}     optional model override\n\n\
GUI diagnostics:\n  \
{gui_profile_env}=1    enable Puffin GUI profiling when built with --features gui-profiler\n  \
{gui_profile_addr_env}  optional profiler bind address (default {gui_profile_addr})\n\n  \
ChatGPT, Claude.ai/Claude Code, and Le Chat subscription/login tokens are not API keys.",
        gui_profile_env = gui_profiler::GUI_PROFILER_ENV,
        gui_profile_addr_env = gui_profiler::GUI_PROFILER_ADDR_ENV,
        gui_profile_addr = gui_profiler::DEFAULT_GUI_PROFILER_ADDR
    )
}

fn print_help() {
    println!("{}", help_text());
}

fn persist_root_window_state() -> bool {
    !cfg!(target_os = "macos")
}

fn main() -> eframe::Result<()> {
    install_panic_logging();
    gui_profiler::init_from_env();
    configure_macos_process_name();

    let args: Vec<String> = env::args().skip(1).collect();
    let cli = match parse_single_project_cli_args(&args, SingleProjectCliOptions::gui()) {
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
    let mut viewport = egui::ViewportBuilder::default()
        .with_inner_size([800.0, 600.0])
        .with_min_inner_size([300.0, 220.0])
        .with_fullscreen(false)
        .with_maximized(false);
    if let Some(icon) = load_icon("assets/icon.png") {
        viewport = viewport.with_icon(icon);
    }
    let options = NativeOptions {
        viewport,
        persist_window: persist_root_window_state(),
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
    use super::{help_text, load_icon, persist_root_window_state, resolve_runtime_asset_path_from};
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

    #[test]
    fn root_window_state_persistence_matches_platform_policy() {
        assert_eq!(persist_root_window_state(), !cfg!(target_os = "macos"));
    }

    #[test]
    fn help_text_points_to_cli_and_agent_api_environment() {
        let help = help_text();
        assert!(help.contains("gentle_cli --help"));
        assert!(help.contains("gentle_mcp"));
        assert!(help.contains("gentle_js"));
        assert!(help.contains("gentle_lua"));
        assert!(help.contains(gentle::agent_bridge::OPENAI_API_KEY_ENV));
        assert!(help.contains(gentle::agent_bridge::ANTHROPIC_API_KEY_ENV));
        assert!(help.contains(gentle::agent_bridge::MISTRAL_API_KEY_ENV));
        assert!(help.contains(gentle::gui_profiler::GUI_PROFILER_ENV));
        assert!(help.contains(gentle::gui_profiler::GUI_PROFILER_ADDR_ENV));
        assert!(help.contains("subscription/login tokens are not API keys"));
    }
}
