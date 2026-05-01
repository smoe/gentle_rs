//! Shell-owned command-adjacent contracts that do not require engine types.

use serde::Serialize;
use serde_json::Value;

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

#[derive(Debug, Clone, Copy)]
/// Parser-local set-algebra enum for `candidates set-op ...` shell commands.
///
/// This stays local to the shell layer because it only exists to parse command
/// text before converting into engine operations.
pub enum CandidateSetOperator {
    Union,
    Intersect,
    Subtract,
}

impl CandidateSetOperator {
    pub fn parse(raw: &str) -> Option<Self> {
        match raw.trim().to_ascii_lowercase().as_str() {
            "union" => Some(Self::Union),
            "intersect" | "intersection" => Some(Self::Intersect),
            "subtract" | "difference" => Some(Self::Subtract),
            _ => None,
        }
    }

    pub fn as_str(self) -> &'static str {
        match self {
            Self::Union => "union",
            Self::Intersect => "intersect",
            Self::Subtract => "subtract",
        }
    }
}

#[cfg(test)]
mod tests {
    use super::{BatchEmitMode, BatchManifestDelimiter, BatchStateMode, CandidateSetOperator};

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
}
