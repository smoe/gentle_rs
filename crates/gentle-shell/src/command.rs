//! Shell-owned command-adjacent contracts that do not require engine types.

use serde::Serialize;
use serde_json::Value;

pub use gentle_protocol::CandidateSetOperator;

/// Deterministic result envelope returned after one shell command executes.
#[derive(Debug, Clone)]
pub struct ShellRunResult {
    pub state_changed: bool,
    pub output: Value,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize)]
pub enum BatchManifestDelimiter {
    Auto,
    Tsv,
    Csv,
}

impl BatchManifestDelimiter {
    pub fn as_str(self) -> &'static str {
        match self {
            Self::Auto => "auto",
            Self::Tsv => "tsv",
            Self::Csv => "csv",
        }
    }

    pub fn parse(raw: &str) -> Result<Self, String> {
        match raw.trim().to_ascii_lowercase().as_str() {
            "auto" => Ok(Self::Auto),
            "tsv" | "tab" | "tabs" => Ok(Self::Tsv),
            "csv" | "comma" => Ok(Self::Csv),
            other => Err(format!(
                "Unsupported batch manifest delimiter '{other}' (expected auto|tsv|csv)"
            )),
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize)]
pub enum BatchStateMode {
    None,
    Shared,
    CopyPerRow,
}

impl BatchStateMode {
    pub fn as_str(self) -> &'static str {
        match self {
            Self::None => "none",
            Self::Shared => "shared",
            Self::CopyPerRow => "copy-per-row",
        }
    }

    pub fn parse(raw: &str) -> Result<Self, String> {
        match raw.trim().to_ascii_lowercase().as_str() {
            "none" | "no-state" => Ok(Self::None),
            "shared" => Ok(Self::Shared),
            "copy-per-row" | "copy_per_row" | "copy" => Ok(Self::CopyPerRow),
            other => Err(format!(
                "Unsupported batch state mode '{other}' (expected none|shared|copy-per-row)"
            )),
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize)]
pub enum BatchEmitMode {
    Local,
    Slurm,
}

impl BatchEmitMode {
    pub fn as_str(self) -> &'static str {
        match self {
            Self::Local => "local",
            Self::Slurm => "slurm",
        }
    }

    pub fn parse(raw: &str) -> Result<Self, String> {
        match raw.trim().to_ascii_lowercase().as_str() {
            "local" | "shell" | "bash" => Ok(Self::Local),
            "slurm" | "sbatch" => Ok(Self::Slurm),
            other => Err(format!(
                "Unsupported batch emit mode '{other}' (expected local|slurm)"
            )),
        }
    }
}

/// Which prepared-genome cache tree a cleanup command should inspect or mutate.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CacheCleanupScope {
    References,
    Helpers,
    Both,
}

impl CacheCleanupScope {
    pub fn label(self) -> &'static str {
        match self {
            Self::References => "references",
            Self::Helpers => "helpers",
            Self::Both => "both",
        }
    }

    pub fn from_flag(raw: &str) -> Option<Self> {
        match raw {
            "--references" => Some(Self::References),
            "--helpers" => Some(Self::Helpers),
            "--both" => Some(Self::Both),
            _ => None,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::{
        BatchEmitMode, BatchManifestDelimiter, BatchStateMode, CacheCleanupScope,
        CandidateSetOperator,
    };

    #[test]
    fn batch_contract_enums_parse_stable_aliases() {
        assert_eq!(
            BatchManifestDelimiter::parse("tab").unwrap().as_str(),
            "tsv"
        );
        assert_eq!(
            BatchStateMode::parse("copy_per_row").unwrap().as_str(),
            "copy-per-row"
        );
        assert_eq!(BatchEmitMode::parse("sbatch").unwrap().as_str(), "slurm");
    }

    #[test]
    fn candidate_set_operator_accepts_shell_aliases() {
        assert_eq!(
            CandidateSetOperator::parse("intersection")
                .expect("parse operator")
                .as_str(),
            "intersect"
        );
    }

    #[test]
    fn cache_cleanup_scope_flags_are_stable() {
        assert_eq!(
            CacheCleanupScope::from_flag("--helpers")
                .expect("parse scope")
                .label(),
            "helpers"
        );
        assert!(CacheCleanupScope::from_flag("--unknown").is_none());
    }
}
