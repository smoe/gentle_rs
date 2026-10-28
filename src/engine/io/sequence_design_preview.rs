//! Fresh-only sequence-design report publication and host-project exclusions.
//! Native directory identity is checked in execution, not parsing or admission.
//! No hostile parent-swap, crash-durability or universal atomicity guarantee.

use crate::engine::{EngineError, ErrorCode, GentleEngine, ProjectState};
use same_file::Handle;
use serde::Serialize;
use std::{
    fs,
    io::{self, Write},
    path::{Path, PathBuf},
};
use tempfile::NamedTempFile;

fn invalid(reason: &str, detail: impl std::fmt::Display) -> EngineError {
    EngineError {
        code: ErrorCode::InvalidInput,
        message: format!("sequence_design_report.{reason}: {detail}"),
        cause_chain: vec![format!("sequence_design_report.{reason}")],
    }
}

fn io_error(reason: &str, error: io::Error) -> EngineError {
    EngineError {
        code: ErrorCode::Io,
        message: format!("sequence_design_report.{reason}: {error}"),
        cause_chain: vec![
            format!("sequence_design_report.{reason}"),
            error.to_string(),
        ],
    }
}

fn directory_identity(path: &Path) -> Result<Handle, EngineError> {
    Handle::from_path(path).map_err(|e| io_error("identity", e))
}

#[cfg(windows)]
fn check_native_path(path: &Path) -> Result<(), EngineError> {
    use std::path::{Component, Prefix};
    if !path.is_absolute() && path.has_root() {
        return Err(invalid(
            "unsupported_path",
            "root-relative paths are not supported",
        ));
    }
    for component in path.components() {
        match component {
            Component::Prefix(prefix) => match prefix.kind() {
                Prefix::Disk(_)
                | Prefix::UNC(..)
                | Prefix::VerbatimDisk(_)
                | Prefix::VerbatimUNC(..)
                    if path.has_root() => {}
                _ => return Err(invalid("unsupported_path", "device or drive-relative path")),
            },
            Component::Normal(name) if name != std::ffi::OsStr::new("..") => {
                let name = name
                    .to_str()
                    .ok_or_else(|| invalid("unsupported_path", "non-text path component"))?;
                let stem = name
                    .split('.')
                    .next()
                    .unwrap_or_default()
                    .to_ascii_uppercase();
                let numbered = ["COM", "LPT"].iter().any(|prefix| {
                    stem.strip_prefix(*prefix).is_some_and(|tail| {
                        tail.len() == 1 && matches!(tail.as_bytes()[0], b'1'..=b'9')
                    })
                });
                if name
                    .chars()
                    .any(|c| matches!(c, ':' | '*' | '?' | '"' | '<' | '>' | '|'))
                    || name.ends_with(' ')
                    || name.ends_with('.')
                    || matches!(stem.as_str(), "CON" | "PRN" | "AUX" | "NUL")
                    || numbered
                {
                    return Err(invalid(
                        "unsupported_path",
                        "stream, device or ambiguous Windows component",
                    ));
                }
            }
            _ => {}
        }
    }
    Ok(())
}

#[cfg(not(windows))]
fn check_native_path(_: &Path) -> Result<(), EngineError> {
    Ok(())
}

/// An admitted absolute destination; final no-clobber publication is the race guard.
#[derive(Debug)]
pub(crate) struct PreviewDestination {
    path: PathBuf,
}

impl GentleEngine {
    pub(crate) fn sequence_design_preview_destination(
        &self,
        raw: &str,
    ) -> Result<PreviewDestination, EngineError> {
        if raw.trim().is_empty()
            || raw.chars().any(char::is_control)
            || raw.ends_with('/')
            || (cfg!(windows) && raw.ends_with('\\'))
        {
            return Err(invalid(
                "invalid_path",
                "expected a nonempty path without control characters",
            ));
        }
        let path = Path::new(raw);
        let raw_leaf = raw
            .rsplit(|c| c == '/' || (cfg!(windows) && c == '\\'))
            .next()
            .unwrap_or_default();
        if matches!(raw_leaf, "." | "..") {
            return Err(invalid("invalid_path", "missing file name"));
        }
        check_native_path(path)?;
        let leaf = path
            .file_name()
            .ok_or_else(|| invalid("invalid_path", "missing file name"))?;
        let parent = path
            .parent()
            .filter(|p| !p.as_os_str().is_empty())
            .unwrap_or(Path::new("."));
        let parent = fs::canonicalize(parent).map_err(|e| io_error("parent", e))?;
        if !fs::metadata(&parent)
            .map_err(|e| io_error("parent", e))?
            .is_dir()
        {
            return Err(invalid(
                "parent",
                "output parent must be an existing directory",
            ));
        }
        let destination = parent.join(leaf);
        match fs::symlink_metadata(&destination) {
            Ok(_) => {
                return Err(invalid(
                    "destination_exists",
                    "existing destinations are never overwritten",
                ));
            }
            Err(error) if error.kind() == io::ErrorKind::NotFound => {}
            Err(error) => return Err(io_error("destination", error)),
        }
        if let Some(project) = self.host_project_path() {
            let project_parent = project
                .parent()
                .ok_or_else(|| invalid("host_path", "missing project parent"))?;
            match fs::canonicalize(project_parent) {
                Err(error) if error.kind() == io::ErrorKind::NotFound => {}
                Err(error) => return Err(io_error("project_parent", error)),
                Ok(project_parent) => {
                    let output_identity = directory_identity(&parent)?;
                    if output_identity == directory_identity(&project_parent)? {
                        let project_leaf = project
                            .file_name()
                            .ok_or_else(|| invalid("host_path", "missing project filename"))?;
                        let sidecar = ProjectState::candidate_store_sidecar_dir(project);
                        let sidecar_leaf = sidecar
                            .file_name()
                            .ok_or_else(|| invalid("host_path", "missing sidecar filename"))?;
                        for protected in [project_leaf, sidecar_leaf] {
                            let (Some(leaf), Some(protected)) = (leaf.to_str(), protected.to_str())
                            else {
                                return Err(invalid(
                                    "unsupported_name_comparison",
                                    "non-text protected filename",
                                ));
                            };
                            if !leaf.is_ascii() || !protected.is_ascii() {
                                return Err(invalid(
                                    "unsupported_name_comparison",
                                    "non-ASCII same-parent names require a different output directory",
                                ));
                            }
                            if leaf.eq_ignore_ascii_case(protected) {
                                return Err(invalid(
                                    "protected_path",
                                    "output is a bound project or sidecar path",
                                ));
                            }
                        }
                    }
                    let sidecar = ProjectState::candidate_store_sidecar_dir(project);
                    match fs::metadata(&sidecar) {
                        Ok(metadata) if metadata.is_dir() => {
                            let sidecar_identity = directory_identity(&sidecar)?;
                            for ancestor in parent.ancestors() {
                                if directory_identity(ancestor)? == sidecar_identity {
                                    return Err(invalid(
                                        "protected_path",
                                        "output is inside the project's save-owned sidecar tree",
                                    ));
                                }
                            }
                        }
                        Ok(_) => {}
                        Err(error) if error.kind() == io::ErrorKind::NotFound => {}
                        Err(error) => return Err(io_error("sidecar", error)),
                    }
                }
            }
        }
        Ok(PreviewDestination { path: destination })
    }
}

fn stage(file: &mut NamedTempFile, text: &str) -> Result<(), EngineError> {
    file.write_all(text.as_bytes())
        .map_err(|e| io_error("write", e))?;
    file.flush().map_err(|e| io_error("flush", e))?;
    file.as_file().sync_all().map_err(|e| io_error("sync", e))
}

impl PreviewDestination {
    pub(crate) fn write<T: Serialize>(&self, value: &T) -> Result<(), EngineError> {
        self.write_with(value, stage)
    }

    // Private staging seam permits deterministic failures/races without global
    // fault flags, permission assumptions or sleeps. Production always uses stage.
    fn write_with<T: Serialize>(
        &self,
        value: &T,
        stage: impl FnOnce(&mut NamedTempFile, &str) -> Result<(), EngineError>,
    ) -> Result<(), EngineError> {
        let text = serde_json::to_string_pretty(value).map_err(|e| {
            EngineError::new(
                ErrorCode::Internal,
                format!("sequence_design_report.serialize: {e}"),
            )
        })?;
        let mut temp = NamedTempFile::new_in(self.path.parent().expect("admitted absolute path"))
            .map_err(|e| io_error("stage", e))?;
        stage(&mut temp, &text)?;
        temp.persist_noclobber(&self.path)
            .map_err(|e| io_error("publish", e.error))?;
        Ok(())
    }
}

#[cfg(test)]
#[path = "sequence_design_preview/tests.rs"]
mod tests;
