//! Guards release-facing version strings against Cargo package-version drift.

use std::{fs, path::Path};

fn read(root: &Path, relative: &str) -> String {
    let path = root.join(relative);
    fs::read_to_string(&path).unwrap_or_else(|error| panic!("read {}: {error}", path.display()))
}

#[test]
fn release_metadata_matches_cargo_package_version() {
    let root = Path::new(env!("CARGO_MANIFEST_DIR"));
    let package_name = env!("CARGO_PKG_NAME");
    let version = env!("CARGO_PKG_VERSION");
    let tag = format!("v{version}");

    let lock = read(root, "Cargo.lock");
    assert!(
        lock.contains(&format!(
            "name = \"{package_name}\"\nversion = \"{version}\""
        )),
        "Cargo.lock must contain the root package at {version}"
    );

    let readme = read(root, "README.md");
    assert!(
        readme.contains(&format!("Current package version: `{version}`.")),
        "README project status must name Cargo package version {version}"
    );
    let badge_version = version.replace('-', "--");
    assert!(
        readme.contains(&format!("version-{badge_version}-")),
        "README version badge must encode Cargo package version {version}"
    );

    if version.contains("-internal.") {
        let release_note_relative = format!("docs/release_notes/release_notes_v{version}.md");
        let release_note = read(root, &release_note_relative);
        assert!(
            release_note
                .lines()
                .next()
                .is_some_and(|line| line.contains(&tag)),
            "{release_note_relative} title must name {tag}"
        );

        let changelog = read(root, "docs/CHANGELOG.md");
        assert!(
            changelog
                .lines()
                .any(|line| line.starts_with("## ") && line.contains(&format!("`{tag}`"))),
            "docs/CHANGELOG.md must contain a version heading for {tag}"
        );

        let release_index = read(root, "docs/release_notes/README.md");
        let release_note_name = format!("release_notes_v{version}.md");
        assert!(
            release_index.contains(&release_note_name),
            "release-note index must link {release_note_name}"
        );
        assert!(
            readme.contains(&format!("docs/release_notes/{release_note_name}")),
            "README must link the current internal release note"
        );
    }
}
