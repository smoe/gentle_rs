use std::time::{SystemTime, UNIX_EPOCH};

fn emit_build_version() {
    let n = SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .map(|d| d.as_millis())
        .unwrap_or(0);
    println!("cargo:rustc-env=GENTLE_BUILD_N={n}");
    println!("cargo:rustc-env=GENTLE_DISPLAY_VERSION=0.0.{n}");
}

fn main() {
    emit_build_version();

    #[cfg(target_os = "windows")]
    {
        let mut res = winres::WindowsResource::new();
        res.set_icon("assets/icon.ico"); // Windows-compatible .ico
        res.compile().expect("Could not compile Windows resources");
    }
}
