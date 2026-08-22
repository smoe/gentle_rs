//! Build script that stamps compile-time version metadata and embeds platform resources.

use std::path::{Path, PathBuf};
use std::process::Command;
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

fn git_stdout(manifest_dir: &Path, args: &[&str]) -> Option<String> {
    let output = Command::new("git")
        .arg("-C")
        .arg(manifest_dir)
        .args(args)
        .output()
        .ok()?;
    if !output.status.success() {
        return None;
    }
    let value = String::from_utf8(output.stdout).ok()?;
    let value = value.trim();
    (!value.is_empty()).then(|| value.to_string())
}

fn normalize_git_commit(value: &str, source: &str) -> String {
    let value = value.trim();
    if !(7..=64).contains(&value.len()) || !value.bytes().all(|byte| byte.is_ascii_hexdigit()) {
        panic!("{source} must contain a 7-64 character hexadecimal Git commit, got '{value}'");
    }
    value.to_ascii_lowercase()
}

fn explicit_git_commit(name: &str) -> Option<String> {
    let value = std::env::var(name).ok()?;
    let value = value.trim();
    (!value.is_empty()).then(|| normalize_git_commit(value, name))
}

fn emit_git_rerun_paths(manifest_dir: &Path) {
    for git_path in ["HEAD", "packed-refs"] {
        if let Some(path) = git_stdout(manifest_dir, &["rev-parse", "--git-path", git_path]) {
            println!("cargo:rerun-if-changed={path}");
        }
    }
    if let Some(branch_ref) = git_stdout(manifest_dir, &["symbolic-ref", "--quiet", "HEAD"])
        && let Some(path) = git_stdout(
            manifest_dir,
            &["rev-parse", "--git-path", branch_ref.as_str()],
        )
    {
        println!("cargo:rerun-if-changed={path}");
    }
}

fn emit_source_revision() {
    println!("cargo:rerun-if-env-changed=GENTLE_GIT_COMMIT");
    println!("cargo:rerun-if-env-changed=GIT_COMMIT_HASH");

    let manifest_dir = std::env::var("CARGO_MANIFEST_DIR")
        .map(PathBuf::from)
        .unwrap_or_else(|_| PathBuf::from("."));
    emit_git_rerun_paths(&manifest_dir);

    let git_commit = explicit_git_commit("GENTLE_GIT_COMMIT")
        .or_else(|| explicit_git_commit("GIT_COMMIT_HASH"))
        .or_else(|| {
            git_stdout(&manifest_dir, &["rev-parse", "--verify", "HEAD"])
                .map(|value| normalize_git_commit(&value, "git rev-parse HEAD"))
        })
        .unwrap_or_else(|| "unknown".to_string());
    let package_version =
        std::env::var("CARGO_PKG_VERSION").unwrap_or_else(|_| "unknown".to_string());
    let source_revision = if git_commit == "unknown" {
        package_version
    } else {
        format!("{package_version}+git.{git_commit}")
    };
    println!("cargo:rustc-env=GENTLE_GIT_COMMIT={git_commit}");
    println!("cargo:rustc-env=GENTLE_SOURCE_REVISION={source_revision}");
}

fn main() {
    // Keep debug iteration incremental: without an explicit rerun rule Cargo
    // re-executes build.rs on any package file change, which would restamp the
    // compile-time version metadata and invalidate more of the crate than
    // needed on every small edit.
    println!("cargo:rerun-if-changed=build.rs");
    verify_embedded_runtime_resources();
    emit_build_version();
    emit_source_revision();

    #[cfg(target_os = "windows")]
    {
        let mut res = winres::WindowsResource::new();
        res.set_icon("assets/icon.ico"); // Windows-compatible .ico
        res.compile().expect("Could not compile Windows resources");
    }
}
