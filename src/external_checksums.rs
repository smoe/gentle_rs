//! External checksum helpers for legacy download verification.
//!
//! SHA-1 is retained only where upstream resources or existing manifests still
//! expose legacy SHA-1 fields. GENtle deliberately delegates such checks to
//! platform tools instead of implementing SHA-1 verification in-process.

use std::fmt;
use std::path::Path;
use std::process::{Command, Output};

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ExternalSha1Digest {
    pub value: String,
    pub tool: String,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum ExternalSha1ErrorKind {
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
}

impl fmt::Display for ExternalSha1Error {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}", self.message)
    }
}

impl std::error::Error for ExternalSha1Error {}

pub fn compute_sha1_with_external_tool(
    path: &Path,
) -> Result<ExternalSha1Digest, ExternalSha1Error> {
    #[cfg(windows)]
    let candidates = [Sha1Tool::Certutil, Sha1Tool::Sha1sum, Sha1Tool::Shasum];
    #[cfg(not(windows))]
    let candidates = [Sha1Tool::Sha1sum, Sha1Tool::Shasum];
    let tool_names = candidates
        .iter()
        .map(|tool| tool.command_name())
        .collect::<Vec<_>>();

    let mut parse_error: Option<ExternalSha1Error> = None;
    for candidate in candidates {
        let output = match candidate.run(path) {
            Ok(output) => output,
            Err(error) if error.kind() == std::io::ErrorKind::NotFound => continue,
            Err(error) => {
                return Err(ExternalSha1Error {
                    kind: ExternalSha1ErrorKind::ToolFailed,
                    message: format!(
                        "could not run external SHA-1 tool '{}': {error}",
                        candidate.command_name()
                    ),
                });
            }
        };
        if !output.status.success() {
            return Err(ExternalSha1Error::failed(candidate.command_name(), &output));
        }
        let stdout = String::from_utf8_lossy(&output.stdout);
        match candidate.parse_output(&stdout) {
            Some(value) => {
                return Ok(ExternalSha1Digest {
                    value,
                    tool: candidate.command_name().to_string(),
                });
            }
            None => {
                parse_error = Some(ExternalSha1Error::parse_failed(
                    candidate.command_name(),
                    &stdout,
                ));
            }
        }
    }
    Err(parse_error.unwrap_or_else(|| ExternalSha1Error::unavailable(&tool_names)))
}

#[derive(Debug, Clone, Copy)]
enum Sha1Tool {
    Sha1sum,
    Shasum,
    #[cfg(windows)]
    Certutil,
}

impl Sha1Tool {
    fn command_name(self) -> &'static str {
        match self {
            Sha1Tool::Sha1sum => "sha1sum",
            Sha1Tool::Shasum => "shasum",
            #[cfg(windows)]
            Sha1Tool::Certutil => "certutil",
        }
    }

    fn run(self, path: &Path) -> Result<Output, std::io::Error> {
        match self {
            Sha1Tool::Sha1sum => Command::new("sha1sum").arg(path).output(),
            Sha1Tool::Shasum => Command::new("shasum").arg("-a").arg("1").arg(path).output(),
            #[cfg(windows)]
            Sha1Tool::Certutil => Command::new("certutil")
                .arg("-hashfile")
                .arg(path)
                .arg("SHA1")
                .output(),
        }
    }

    fn parse_output(self, stdout: &str) -> Option<String> {
        match self {
            Sha1Tool::Sha1sum | Sha1Tool::Shasum => parse_sha1sum_like_output(stdout),
            #[cfg(windows)]
            Sha1Tool::Certutil => parse_certutil_output(stdout),
        }
    }
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

#[cfg(any(windows, test))]
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
            Err(error) if error.is_tool_unavailable() => {}
            Err(error) => panic!("unexpected external SHA-1 error: {error}"),
        }
    }
}
