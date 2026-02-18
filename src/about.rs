pub const GENTLE_DISPLAY_VERSION: &str = env!("GENTLE_DISPLAY_VERSION");
pub const GENTLE_BUILD_N: &str = env!("GENTLE_BUILD_N");

pub fn version_cli_text() -> String {
    format!(
        "GENtle {}\nBuild {}\nCross-platform DNA cloning workbench",
        GENTLE_DISPLAY_VERSION, GENTLE_BUILD_N
    )
}

#[cfg(target_os = "macos")]
pub fn show_native_about_panel() -> bool {
    use objc2_app_kit::NSApplication;
    use objc2_foundation::MainThreadMarker;

    let Some(mtm) = MainThreadMarker::new() else {
        return false;
    };
    let app = NSApplication::sharedApplication(mtm);
    unsafe {
        app.orderFrontStandardAboutPanel(None);
    }
    true
}

#[cfg(not(target_os = "macos"))]
pub fn show_native_about_panel() -> bool {
    false
}
