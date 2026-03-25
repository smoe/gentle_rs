//! Build script that stamps compile-time version metadata and embeds platform resources.

use std::time::{SystemTime, UNIX_EPOCH};

fn display_version_prefix() -> String {
    let pkg_version = std::env::var("CARGO_PKG_VERSION").unwrap_or_else(|_| "0.1.0".to_string());
    let mut parts = pkg_version.split('.');
    let major = parts.next().unwrap_or("0");
    let minor = parts.next().unwrap_or("1");
    format!("{major}.{minor}")
}

fn emit_build_version() {
    let n = SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .map(|d| d.as_millis())
        .unwrap_or(0);
    let display_prefix = display_version_prefix();
    println!("cargo:rustc-env=GENTLE_BUILD_N={n}");
    println!("cargo:rustc-env=GENTLE_DISPLAY_VERSION={display_prefix}.{n}");
}

fn main() {
    // Keep debug iteration incremental: without an explicit rerun rule Cargo
    // re-executes build.rs on any package file change, which would restamp the
    // compile-time version metadata and invalidate more of the crate than
    // needed on every small edit.
    println!("cargo:rerun-if-changed=build.rs");
    emit_build_version();

    #[cfg(target_os = "windows")]
    {
        let mut res = winres::WindowsResource::new();
        res.set_icon("assets/icon.ico"); // Windows-compatible .ico
        res.compile().expect("Could not compile Windows resources");
    }
}
