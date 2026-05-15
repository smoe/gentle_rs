//! Agent-loop diagnostics for the direct `gentle_cli doctor --agent` route.

use super::{now_unix_ms, print_json};
use serde::Serialize;
use std::{
    env,
    fs::{self, OpenOptions},
    path::{Path, PathBuf},
    process::Command,
    time::UNIX_EPOCH,
};

#[derive(Debug, Serialize)]
struct AgentDoctorReport {
    schema: &'static str,
    ok: bool,
    cargo: AgentDoctorCargoCheck,
    binary: AgentDoctorBinaryFreshnessCheck,
    state: AgentDoctorStateWritabilityCheck,
    examples_docs_check: AgentDoctorCommandCheck,
}

#[derive(Debug, Serialize)]
struct AgentDoctorCargoCheck {
    ok: bool,
    program: &'static str,
    version: Option<String>,
    error: Option<String>,
}

#[derive(Debug, Serialize)]
struct AgentDoctorBinaryFreshnessCheck {
    ok: bool,
    binary_path: Option<String>,
    binary_exists: bool,
    binary_mtime_unix: Option<u64>,
    newest_source_path: Option<String>,
    newest_source_mtime_unix: Option<u64>,
    fresh: Option<bool>,
    error: Option<String>,
}

#[derive(Debug, Serialize)]
struct AgentDoctorStateWritabilityCheck {
    ok: bool,
    path: String,
    parent_path: String,
    parent_exists: bool,
    parent_writable: bool,
    file_writable_if_exists: bool,
    error: Option<String>,
}

#[derive(Debug, Serialize)]
struct AgentDoctorCommandCheck {
    ok: bool,
    command: Vec<String>,
    status_code: Option<i32>,
    stdout_tail: Option<String>,
    stderr_tail: Option<String>,
    error: Option<String>,
}

pub(super) fn run_agent_doctor(state_path: &str) -> Result<(), String> {
    let report = build_agent_doctor_report(state_path);
    let ok = report.ok;
    print_json(&report)?;
    if ok {
        Ok(())
    } else {
        std::process::exit(1);
    }
}

fn build_agent_doctor_report(state_path: &str) -> AgentDoctorReport {
    let cargo = agent_doctor_cargo_check();
    let binary = agent_doctor_binary_freshness_check();
    let state = agent_doctor_state_writability_check(state_path);
    let examples_docs_check = agent_doctor_examples_docs_check();
    let ok = cargo.ok && binary.ok && state.ok && examples_docs_check.ok;

    AgentDoctorReport {
        schema: "gentle.agent_dev_doctor.v1",
        ok,
        cargo,
        binary,
        state,
        examples_docs_check,
    }
}

fn agent_doctor_cargo_check() -> AgentDoctorCargoCheck {
    match Command::new("cargo").arg("--version").output() {
        Ok(output) => {
            let version = agent_doctor_tail(&output.stdout, 400);
            let stderr = agent_doctor_tail(&output.stderr, 400);
            AgentDoctorCargoCheck {
                ok: output.status.success(),
                program: "cargo",
                version,
                error: if output.status.success() {
                    None
                } else {
                    Some(stderr.unwrap_or_else(|| {
                        format!("cargo --version exited with {:?}", output.status.code())
                    }))
                },
            }
        }
        Err(e) => AgentDoctorCargoCheck {
            ok: false,
            program: "cargo",
            version: None,
            error: Some(e.to_string()),
        },
    }
}

fn agent_doctor_binary_freshness_check() -> AgentDoctorBinaryFreshnessCheck {
    let binary_path = match env::current_exe() {
        Ok(path) => path,
        Err(e) => {
            return AgentDoctorBinaryFreshnessCheck {
                ok: false,
                binary_path: None,
                binary_exists: false,
                binary_mtime_unix: None,
                newest_source_path: None,
                newest_source_mtime_unix: None,
                fresh: None,
                error: Some(format!("Could not resolve current executable: {e}")),
            };
        }
    };
    let binary_exists = binary_path.exists();
    let binary_mtime_unix = fs::metadata(&binary_path)
        .and_then(|metadata| metadata.modified())
        .ok()
        .and_then(agent_doctor_unix_secs);
    let repo_root = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    let newest_source = agent_doctor_newest_source_mtime(&repo_root);
    let (newest_source_path, newest_source_mtime_unix) = newest_source
        .map(|(path, mtime)| (Some(path.display().to_string()), Some(mtime)))
        .unwrap_or((None, None));
    let fresh = match (binary_mtime_unix, newest_source_mtime_unix) {
        (Some(binary_mtime), Some(source_mtime)) => Some(binary_mtime >= source_mtime),
        _ => None,
    };
    let ok = binary_exists && fresh == Some(true);
    let error = if ok {
        None
    } else if !binary_exists {
        Some("Current executable path does not exist".to_string())
    } else if fresh == Some(false) {
        Some("Current executable is older than at least one source file".to_string())
    } else {
        Some("Could not compare executable and source mtimes".to_string())
    };

    AgentDoctorBinaryFreshnessCheck {
        ok,
        binary_path: Some(binary_path.display().to_string()),
        binary_exists,
        binary_mtime_unix,
        newest_source_path,
        newest_source_mtime_unix,
        fresh,
        error,
    }
}

fn agent_doctor_state_writability_check(state_path: &str) -> AgentDoctorStateWritabilityCheck {
    let path = Path::new(state_path);
    let parent = path
        .parent()
        .filter(|parent| !parent.as_os_str().is_empty())
        .unwrap_or_else(|| Path::new("."));
    let parent_exists = parent.exists();
    let mut errors = Vec::new();
    let parent_writable = if parent_exists {
        let probe_path = parent.join(format!(
            ".gentle_doctor_write_test_{}_{}",
            std::process::id(),
            now_unix_ms()
        ));
        match OpenOptions::new()
            .write(true)
            .create_new(true)
            .open(&probe_path)
        {
            Ok(_) => {
                if let Err(e) = fs::remove_file(&probe_path) {
                    errors.push(format!(
                        "Could not remove state write probe '{}': {e}",
                        probe_path.display()
                    ));
                }
                true
            }
            Err(e) => {
                errors.push(format!(
                    "Could not create state write probe in '{}': {e}",
                    parent.display()
                ));
                false
            }
        }
    } else {
        errors.push(format!(
            "State parent directory '{}' does not exist",
            parent.display()
        ));
        false
    };

    let file_writable_if_exists = if path.exists() {
        match OpenOptions::new().write(true).open(path) {
            Ok(_) => true,
            Err(e) => {
                errors.push(format!(
                    "Existing state path '{}' is not writable: {e}",
                    path.display()
                ));
                false
            }
        }
    } else {
        true
    };
    let ok = parent_writable && file_writable_if_exists;

    AgentDoctorStateWritabilityCheck {
        ok,
        path: state_path.to_string(),
        parent_path: parent.display().to_string(),
        parent_exists,
        parent_writable,
        file_writable_if_exists,
        error: if errors.is_empty() {
            None
        } else {
            Some(errors.join("; "))
        },
    }
}

fn agent_doctor_examples_docs_check() -> AgentDoctorCommandCheck {
    let command = vec![
        "cargo".to_string(),
        "run".to_string(),
        "-q".to_string(),
        "--bin".to_string(),
        "gentle_examples_docs".to_string(),
        "--".to_string(),
        "--check".to_string(),
    ];
    match Command::new("cargo")
        .args([
            "run",
            "-q",
            "--bin",
            "gentle_examples_docs",
            "--",
            "--check",
        ])
        .output()
    {
        Ok(output) => AgentDoctorCommandCheck {
            ok: output.status.success(),
            command,
            status_code: output.status.code(),
            stdout_tail: agent_doctor_tail(&output.stdout, 2000),
            stderr_tail: agent_doctor_tail(&output.stderr, 2000),
            error: if output.status.success() {
                None
            } else {
                Some(format!(
                    "gentle_examples_docs --check exited with {:?}",
                    output.status.code()
                ))
            },
        },
        Err(e) => AgentDoctorCommandCheck {
            ok: false,
            command,
            status_code: None,
            stdout_tail: None,
            stderr_tail: None,
            error: Some(e.to_string()),
        },
    }
}

fn agent_doctor_newest_source_mtime(repo_root: &Path) -> Option<(PathBuf, u64)> {
    let mut newest: Option<(PathBuf, u64)> = None;
    for relative in ["src", "crates", "Cargo.toml", "Cargo.lock", "build.rs"] {
        let path = repo_root.join(relative);
        if path.exists() {
            agent_doctor_visit_source_mtimes(&path, &mut newest);
        }
    }
    newest
}

fn agent_doctor_visit_source_mtimes(path: &Path, newest: &mut Option<(PathBuf, u64)>) {
    let metadata = match fs::metadata(path) {
        Ok(metadata) => metadata,
        Err(_) => return,
    };
    if metadata.is_dir() {
        let entries = match fs::read_dir(path) {
            Ok(entries) => entries,
            Err(_) => return,
        };
        for entry in entries.flatten() {
            agent_doctor_visit_source_mtimes(&entry.path(), newest);
        }
        return;
    }
    if !metadata.is_file() || !agent_doctor_is_source_file(path) {
        return;
    }
    let Some(mtime) = metadata.modified().ok().and_then(agent_doctor_unix_secs) else {
        return;
    };
    if newest
        .as_ref()
        .map(|(_, seen)| mtime > *seen)
        .unwrap_or(true)
    {
        *newest = Some((path.to_path_buf(), mtime));
    }
}

fn agent_doctor_is_source_file(path: &Path) -> bool {
    path.file_name()
        .and_then(|name| name.to_str())
        .map(|name| name == "Cargo.lock" || name == "Cargo.toml" || name == "build.rs")
        .unwrap_or(false)
        || path
            .extension()
            .and_then(|extension| extension.to_str())
            .map(|extension| matches!(extension, "rs" | "toml" | "lock"))
            .unwrap_or(false)
}

fn agent_doctor_unix_secs(time: std::time::SystemTime) -> Option<u64> {
    time.duration_since(UNIX_EPOCH).ok().map(|d| d.as_secs())
}

fn agent_doctor_tail(bytes: &[u8], max_chars: usize) -> Option<String> {
    let text = String::from_utf8_lossy(bytes).trim().to_string();
    if text.is_empty() {
        return None;
    }
    let char_count = text.chars().count();
    if char_count <= max_chars {
        return Some(text);
    }
    let tail = text
        .chars()
        .skip(char_count.saturating_sub(max_chars))
        .collect::<String>();
    Some(tail)
}

#[cfg(test)]
mod tests {
    use super::*;
    use tempfile::tempdir;

    #[test]
    fn state_writability_accepts_temp_state_path() {
        let dir = tempdir().expect("tempdir");
        let state_path = dir.path().join("agent-session.gentle.json");
        let check = agent_doctor_state_writability_check(&state_path.display().to_string());
        assert!(check.ok, "{check:?}");
        assert!(check.parent_writable);
        assert!(check.file_writable_if_exists);
    }

    #[test]
    fn tail_keeps_short_text_and_trims_long_text() {
        assert_eq!(agent_doctor_tail(b"short\n", 20).as_deref(), Some("short"));
        assert_eq!(agent_doctor_tail(b"abcdef", 3).as_deref(), Some("def"));
    }
}
