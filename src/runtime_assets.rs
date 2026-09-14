//! Runtime lookup for shipped assets from a checkout or packaged desktop app.

use std::path::{Path, PathBuf};

/// Resolve a built-in asset; explicit user paths should not use this fallback.
pub fn resolve_runtime_asset_path(path: &str) -> PathBuf {
    let current_exe = std::env::current_exe().ok();
    resolve_runtime_asset_path_from(
        path,
        current_exe.as_deref(),
        option_env!("CARGO_MANIFEST_DIR").map(Path::new),
    )
}

/// Inject runtime locations for deterministic asset-discovery tests.
/// A direct path takes precedence; unavailable assets retain the requested path.
pub fn resolve_runtime_asset_path_from(
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
    if let Some(exe_dir) = current_exe.and_then(Path::parent) {
        let mut candidates = vec![exe_dir.join(path)];
        if let Some(parent) = exe_dir.parent() {
            candidates.push(parent.join("Resources").join(path));
            if let Some(grandparent) = parent.parent() {
                candidates.push(grandparent.join(path));
            }
        }
        if let Some(candidate) = candidates.into_iter().find(|candidate| candidate.exists()) {
            return candidate;
        }
    }
    direct
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn packaged_catalog_is_found_without_checkout_or_cwd_assets() {
        let temp = tempfile::tempdir().unwrap();
        let relative = "synthetic-assets/catalog";
        let macos = temp.path().join("GENtle.app/Contents/MacOS/gentle");
        let resources = temp
            .path()
            .join("GENtle.app/Contents/Resources")
            .join(relative);
        std::fs::create_dir_all(&resources).unwrap();
        assert_eq!(
            resolve_runtime_asset_path_from(relative, Some(&macos), None),
            resources
        );
        let linux = temp.path().join("linux/gentle");
        let assets = temp.path().join("linux").join(relative);
        std::fs::create_dir_all(&assets).unwrap();
        assert_eq!(
            resolve_runtime_asset_path_from(relative, Some(&linux), None),
            assets
        );
    }

    #[test]
    fn checkout_fallback_and_explicit_path_are_preserved() {
        let temp = tempfile::tempdir().unwrap();
        let relative = "synthetic-assets/catalog";
        let catalog = temp.path().join(relative);
        std::fs::create_dir_all(&catalog).unwrap();
        assert_eq!(
            resolve_runtime_asset_path_from(relative, None, Some(temp.path())),
            catalog
        );
        assert_eq!(
            resolve_runtime_asset_path_from(catalog.to_str().unwrap(), None, None),
            catalog
        );
        let missing = temp.path().join("missing");
        assert_eq!(
            resolve_runtime_asset_path_from(missing.to_str().unwrap(), None, None),
            missing
        );
    }
}
