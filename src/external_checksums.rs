//! External checksum helpers for legacy download verification.
//!
//! SHA-1 is retained only where upstream resources or existing manifests still
//! expose legacy SHA-1 fields. GENtle deliberately delegates such checks to
//! platform tools instead of implementing SHA-1 verification in-process.

use serde::Serialize;
use std::fmt;
use std::io::Write;
use std::path::Path;
use std::process::{Command, Output};

pub const SHA1_TOOL_ENV_BIN: &str = "GENTLE_SHA1_TOOL";
pub const SHA1_DISABLE_ENV: &str = "GENTLE_DISABLE_LEGACY_SHA1";

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ExternalSha1Digest {
    pub value: String,
    pub tool: String,
    pub executable: String,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum ExternalSha1ErrorKind {
    Disabled,
    ToolUnavailable,
    ToolFailed,
    ParseFailed,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ExternalSha1Error {
    pub kind: ExternalSha1ErrorKind,
    pub message: String,
}

impl ExternalSha1Error {
    fn disabled() -> Self {
        Self {
            kind: ExternalSha1ErrorKind::Disabled,
            message: format!("legacy SHA-1 verification disabled by {SHA1_DISABLE_ENV}"),
        }
    }

    fn unavailable(tools: &[&str]) -> Self {
        Self {
            kind: ExternalSha1ErrorKind::ToolUnavailable,
            message: format!(
                "no external SHA-1 tool available; tried {}",
                tools.join(", ")
            ),
        }
    }

    fn failed(tool: &str, output: &Output) -> Self {
        let stderr = String::from_utf8_lossy(&output.stderr).trim().to_string();
        let stdout = String::from_utf8_lossy(&output.stdout).trim().to_string();
        let detail = if !stderr.is_empty() {
            stderr
        } else if !stdout.is_empty() {
            stdout
        } else {
            format!("exit status {}", output.status)
        };
        Self {
            kind: ExternalSha1ErrorKind::ToolFailed,
            message: format!("external SHA-1 tool '{tool}' failed: {detail}"),
        }
    }

    fn parse_failed(tool: &str, stdout: &str) -> Self {
        Self {
            kind: ExternalSha1ErrorKind::ParseFailed,
            message: format!("could not parse SHA-1 digest from '{tool}' output: {stdout:?}"),
        }
    }

    pub fn is_tool_unavailable(&self) -> bool {
        self.kind == ExternalSha1ErrorKind::ToolUnavailable
    }

    pub fn allows_unverified_fallback(&self) -> bool {
        matches!(
            self.kind,
            ExternalSha1ErrorKind::Disabled | ExternalSha1ErrorKind::ToolUnavailable
        )
    }
}

impl fmt::Display for ExternalSha1Error {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}", self.message)
    }
}

impl std::error::Error for ExternalSha1Error {}

#[derive(Debug, Clone, Serialize)]
pub struct LegacySha1ToolStatus {
    pub resource_id: String,
    pub display_name: String,
    pub support_status: String,
    pub env_var: String,
    pub disable_env_var: String,
    pub configured_executable: Option<String>,
    pub resolved_executable: Option<String>,
    pub tool_kind: Option<String>,
    pub available: bool,
    pub probe_digest: Option<String>,
    pub error: Option<String>,
    pub notes: Vec<String>,
}

pub fn compute_sha1_with_external_tool(
    path: &Path,
) -> Result<ExternalSha1Digest, ExternalSha1Error> {
    if legacy_sha1_disabled() {
        return Err(ExternalSha1Error::disabled());
    }
    let candidates = sha1_tool_candidates();
    let tool_names = candidates
        .iter()
        .map(|tool| tool.tool_kind.command_name())
        .collect::<Vec<_>>();

    let mut parse_error: Option<ExternalSha1Error> = None;
    for candidate in candidates {
        let output = match candidate.run(path) {
            Ok(output) => output,
            Err(error) if error.kind() == std::io::ErrorKind::NotFound && !candidate.configured => {
                continue;
            }
            Err(error) => {
                return Err(ExternalSha1Error {
                    kind: ExternalSha1ErrorKind::ToolFailed,
                    message: format!(
                        "could not run external SHA-1 tool '{}': {error}",
                        candidate.executable
                    ),
                });
            }
        };
        if !output.status.success() {
            return Err(ExternalSha1Error::failed(&candidate.executable, &output));
        }
        let stdout = String::from_utf8_lossy(&output.stdout);
        match candidate.parse_output(&stdout) {
            Some(value) => {
                return Ok(ExternalSha1Digest {
                    value,
                    tool: candidate.tool_kind.command_name().to_string(),
                    executable: candidate.executable,
                });
            }
            None => {
                parse_error = Some(ExternalSha1Error::parse_failed(
                    &candidate.executable,
                    &stdout,
                ));
            }
        }
    }
    Err(parse_error.unwrap_or_else(|| ExternalSha1Error::unavailable(&tool_names)))
}

pub fn legacy_sha1_tool_status() -> LegacySha1ToolStatus {
    let configured_executable = configured_sha1_tool();
    let notes = vec![
        "Optional legacy helper used only when upstream resources or compatibility manifests expose SHA-1 fields.".to_string(),
        "GENtle does not implement SHA-1 verification in-process; absence of this tool means SHA-1 verification is skipped.".to_string(),
        format!("Override executable with {SHA1_TOOL_ENV_BIN}; disable with {SHA1_DISABLE_ENV}=1."),
    ];
    if legacy_sha1_disabled() {
        return LegacySha1ToolStatus {
            resource_id: "legacy_sha1".to_string(),
            display_name: "Legacy SHA-1 verification tool".to_string(),
            support_status: "disabled".to_string(),
            env_var: SHA1_TOOL_ENV_BIN.to_string(),
            disable_env_var: SHA1_DISABLE_ENV.to_string(),
            configured_executable,
            resolved_executable: None,
            tool_kind: None,
            available: false,
            probe_digest: None,
            error: Some(format!("Disabled by {SHA1_DISABLE_ENV}")),
            notes,
        };
    }

    let mut probe = match tempfile::NamedTempFile::new() {
        Ok(file) => file,
        Err(error) => {
            return LegacySha1ToolStatus {
                resource_id: "legacy_sha1".to_string(),
                display_name: "Legacy SHA-1 verification tool".to_string(),
                support_status: "probe_failed".to_string(),
                env_var: SHA1_TOOL_ENV_BIN.to_string(),
                disable_env_var: SHA1_DISABLE_ENV.to_string(),
                configured_executable,
                resolved_executable: None,
                tool_kind: None,
                available: false,
                probe_digest: None,
                error: Some(format!(
                    "Could not create temporary SHA-1 probe file: {error}"
                )),
                notes,
            };
        }
    };
    if let Err(error) = probe.write_all(b"abc") {
        return LegacySha1ToolStatus {
            resource_id: "legacy_sha1".to_string(),
            display_name: "Legacy SHA-1 verification tool".to_string(),
            support_status: "probe_failed".to_string(),
            env_var: SHA1_TOOL_ENV_BIN.to_string(),
            disable_env_var: SHA1_DISABLE_ENV.to_string(),
            configured_executable,
            resolved_executable: None,
            tool_kind: None,
            available: false,
            probe_digest: None,
            error: Some(format!(
                "Could not write temporary SHA-1 probe file: {error}"
            )),
            notes,
        };
    }

    match compute_sha1_with_external_tool(probe.path()) {
        Ok(digest) => LegacySha1ToolStatus {
            resource_id: "legacy_sha1".to_string(),
            display_name: "Legacy SHA-1 verification tool".to_string(),
            support_status: "ready_runtime".to_string(),
            env_var: SHA1_TOOL_ENV_BIN.to_string(),
            disable_env_var: SHA1_DISABLE_ENV.to_string(),
            configured_executable,
            resolved_executable: Some(digest.executable),
            tool_kind: Some(digest.tool),
            available: true,
            probe_digest: Some(digest.value),
            error: None,
            notes,
        },
        Err(error) if error.kind == ExternalSha1ErrorKind::ToolUnavailable => {
            LegacySha1ToolStatus {
                resource_id: "legacy_sha1".to_string(),
                display_name: "Legacy SHA-1 verification tool".to_string(),
                support_status: "missing".to_string(),
                env_var: SHA1_TOOL_ENV_BIN.to_string(),
                disable_env_var: SHA1_DISABLE_ENV.to_string(),
                configured_executable,
                resolved_executable: None,
                tool_kind: None,
                available: false,
                probe_digest: None,
                error: Some(error.message),
                notes,
            }
        }
        Err(error) => LegacySha1ToolStatus {
            resource_id: "legacy_sha1".to_string(),
            display_name: "Legacy SHA-1 verification tool".to_string(),
            support_status: "runtime_invalid".to_string(),
            env_var: SHA1_TOOL_ENV_BIN.to_string(),
            disable_env_var: SHA1_DISABLE_ENV.to_string(),
            configured_executable,
            resolved_executable: None,
            tool_kind: None,
            available: false,
            probe_digest: None,
            error: Some(error.message),
            notes,
        },
    }
}

pub fn legacy_sha1_unverified_warning(context: &str) -> Option<String> {
    let status = legacy_sha1_tool_status();
    if status.available {
        None
    } else {
        let reason = status
            .error
            .as_deref()
            .unwrap_or("no external SHA-1 verification tool available");
        Some(format!(
            "Legacy SHA-1 verification skipped for {context}: {reason}. Basic corruption, size, gzip, and parser checks may still apply."
        ))
    }
}

#[derive(Debug, Clone, Copy)]
enum Sha1ToolKind {
    Sha1sum,
    Shasum,
    Certutil,
}

impl Sha1ToolKind {
    fn command_name(self) -> &'static str {
        match self {
            Sha1ToolKind::Sha1sum => "sha1sum",
            Sha1ToolKind::Shasum => "shasum",
            Sha1ToolKind::Certutil => "certutil",
        }
    }
}

#[derive(Debug, Clone)]
struct Sha1Tool {
    tool_kind: Sha1ToolKind,
    executable: String,
    configured: bool,
}

impl Sha1Tool {
    fn run(&self, path: &Path) -> Result<Output, std::io::Error> {
        match self.tool_kind {
            Sha1ToolKind::Sha1sum => Command::new(&self.executable).arg(path).output(),
            Sha1ToolKind::Shasum => Command::new(&self.executable)
                .arg("-a")
                .arg("1")
                .arg(path)
                .output(),
            Sha1ToolKind::Certutil => Command::new(&self.executable)
                .arg("-hashfile")
                .arg(path)
                .arg("SHA1")
                .output(),
        }
    }

    fn parse_output(&self, stdout: &str) -> Option<String> {
        match self.tool_kind {
            Sha1ToolKind::Sha1sum | Sha1ToolKind::Shasum => parse_sha1sum_like_output(stdout),
            Sha1ToolKind::Certutil => parse_certutil_output(stdout),
        }
    }
}

fn configured_sha1_tool() -> Option<String> {
    let value = crate::tool_overrides::configured_or_env(SHA1_TOOL_ENV_BIN);
    let trimmed = value.trim();
    if trimmed.is_empty() {
        None
    } else {
        Some(trimmed.to_string())
    }
}

fn legacy_sha1_disabled() -> bool {
    let value = std::env::var(SHA1_DISABLE_ENV).unwrap_or_default();
    matches!(
        value.trim().to_ascii_lowercase().as_str(),
        "1" | "true" | "yes" | "on"
    )
}

fn infer_sha1_tool_kind(executable: &str) -> Sha1ToolKind {
    let basename = Path::new(executable)
        .file_name()
        .and_then(|value| value.to_str())
        .unwrap_or(executable)
        .to_ascii_lowercase();
    if basename.contains("shasum") {
        Sha1ToolKind::Shasum
    } else if basename.contains("certutil") {
        Sha1ToolKind::Certutil
    } else {
        Sha1ToolKind::Sha1sum
    }
}

fn sha1_tool_candidates() -> Vec<Sha1Tool> {
    if let Some(configured) = configured_sha1_tool() {
        return vec![Sha1Tool {
            tool_kind: infer_sha1_tool_kind(&configured),
            executable: configured,
            configured: true,
        }];
    }
    let mut candidates = Vec::new();
    #[cfg(windows)]
    candidates.push(Sha1Tool {
        tool_kind: Sha1ToolKind::Certutil,
        executable: "certutil".to_string(),
        configured: false,
    });
    candidates.push(Sha1Tool {
        tool_kind: Sha1ToolKind::Sha1sum,
        executable: "sha1sum".to_string(),
        configured: false,
    });
    candidates.push(Sha1Tool {
        tool_kind: Sha1ToolKind::Shasum,
        executable: "shasum".to_string(),
        configured: false,
    });
    candidates
}

fn normalize_sha1_token(token: &str) -> Option<String> {
    let cleaned = token.trim();
    if cleaned.len() == 40 && cleaned.chars().all(|c| c.is_ascii_hexdigit()) {
        Some(cleaned.to_ascii_lowercase())
    } else {
        None
    }
}

fn parse_sha1sum_like_output(stdout: &str) -> Option<String> {
    stdout
        .lines()
        .filter_map(|line| line.split_whitespace().next())
        .find_map(normalize_sha1_token)
}

fn parse_certutil_output(stdout: &str) -> Option<String> {
    stdout.lines().find_map(|line| {
        let compact = line
            .chars()
            .filter(|c| !c.is_ascii_whitespace())
            .collect::<String>();
        normalize_sha1_token(&compact)
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parse_sha1sum_output() {
        let digest = "0123456789abcdef0123456789abcdef01234567";
        assert_eq!(
            parse_sha1sum_like_output(&format!("{digest}  downloaded.fa\n")),
            Some(digest.to_string())
        );
    }

    #[test]
    fn parse_shasum_output() {
        let digest = "abcdef0123456789abcdef0123456789abcdef01";
        assert_eq!(
            parse_sha1sum_like_output(&format!("{digest}  /tmp/downloaded.fa\n")),
            Some(digest.to_string())
        );
    }

    #[test]
    fn parse_certutil_output_with_spaced_digest() {
        let digest = "AB CD EF 01 23 45 67 89 AB CD EF 01 23 45 67 89 AB CD EF 01";
        assert_eq!(
            parse_certutil_output(&format!(
                "SHA1 hash of file downloaded.fa:\n{digest}\nCertUtil: -hashfile command completed successfully.\n"
            )),
            Some("abcdef0123456789abcdef0123456789abcdef01".to_string())
        );
    }

    #[test]
    fn rejects_non_digest_output() {
        assert_eq!(parse_sha1sum_like_output("not a digest\n"), None);
        assert_eq!(parse_certutil_output("CertUtil: command failed\n"), None);
    }

    #[test]
    fn computes_known_digest_when_external_tool_is_available() {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("abc.txt");
        std::fs::write(&path, b"abc").unwrap();
        match compute_sha1_with_external_tool(&path) {
            Ok(digest) => {
                assert_eq!(digest.value, "a9993e364706816aba3e25717850c26c9cd0d89d");
                assert!(!digest.tool.trim().is_empty());
            }
            Err(error) if error.allows_unverified_fallback() => {}
            Err(error) => panic!("unexpected external SHA-1 error: {error}"),
        }
    }
}
