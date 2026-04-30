//! Build script that stamps compile-time version metadata and embeds platform resources.

use std::path::PathBuf;
use std::time::{SystemTime, UNIX_EPOCH};

fn verify_embedded_runtime_resources() {
    let manifest_dir = std::env::var("CARGO_MANIFEST_DIR")
        .map(PathBuf::from)
        .unwrap_or_else(|_| PathBuf::from("."));
    let cwd = std::env::current_dir().unwrap_or_else(|_| PathBuf::from("."));

    // These paths are referenced by include_str!/include_image! in non-test code.
    // Verifying them here gives one actionable message with cwd + manifest context
    // instead of a long macro-expansion error cascade.
    let required_files = [
        "assets/amino_acids.json",
        "assets/codon_catalog.csv",
        "assets/codon_tables.json",
        "assets/dna_ladders.json",
        "assets/enzymes.json",
        "assets/icon.png",
        "assets/jaspar.motifs.json",
        "assets/rna_ladders.json",
        "docs/agent_interface.md",
        "docs/cli.md",
        "docs/examples/protocol_cartoon/oe_substitution_figure1_template.json",
        "docs/glossary.json",
        "docs/gui.md",
        "docs/reviewer_preview.md",
        "icons/GENtle.png",
        "icons/accessories-calculator.png",
        "icons/display_circular_linear.png",
        "icons/display_features.png",
        "icons/display_orfs.png",
        "icons/gc_content.png",
        "icons/restriction_enzymes.png",
        "icons/show_map.png",
        "icons/show_sequence.png",
    ];

    let mut missing: Vec<(String, String)> = vec![];
    for rel in required_files {
        println!("cargo:rerun-if-changed={rel}");
        let abs = manifest_dir.join(rel);
        if !abs.exists() {
            missing.push((rel.to_string(), abs.display().to_string()));
        }
    }
    if missing.is_empty() {
        return;
    }

    let mut message = String::new();
    message.push_str(
        "Missing required embedded resource file(s) for compile-time include_* macros.\n",
    );
    message.push_str(&format!("CARGO_MANIFEST_DIR: {}\n", manifest_dir.display()));
    message.push_str(&format!("Current working directory: {}\n", cwd.display()));
    message.push_str("Missing paths (relative -> absolute):\n");
    for (rel, abs) in missing {
        message.push_str(&format!("  - {rel} -> {abs}\n"));
    }
    message.push_str(
        "If this happens in Docker/CI, ensure assets/docs/icons are copied before `cargo build`.",
    );
    panic!("{message}");
}

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
    verify_embedded_runtime_resources();
    emit_build_version();

    #[cfg(target_os = "windows")]
    {
        let mut res = winres::WindowsResource::new();
        res.set_icon("assets/icon.ico"); // Windows-compatible .ico
        res.compile().expect("Could not compile Windows resources");
    }
}
