
#[cfg(target_os = "windows")]
fn main() {
    let mut res = winres::WindowsResource::new();
    res.set_icon("assets/icon.ico"); // Windows-compatible .ico
    res.compile().unwrap();
}

#[cfg(not(target_os = "windows"))]
fn main() {
    // Do nothing for macOS or Linux
}
