use crate::dna_sequence::DNAsequence;
use serde::{Deserialize, Serialize};
use std::{
    fmt,
    io::ErrorKind,
    path::Path,
    process::{Command, Output},
};

const DEFAULT_RNAPKIN_BIN: &str = "rnapkin";
const RNAPKIN_ENV_BIN: &str = "GENTLE_RNAPKIN_BIN";

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct RnaStructureTextReport {
    pub tool: String,
    pub executable: String,
    pub sequence_length: usize,
    pub command: Vec<String>,
    pub stdout: String,
    pub stderr: String,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct RnaStructureSvgReport {
    pub tool: String,
    pub executable: String,
    pub sequence_length: usize,
    pub output_path: String,
    pub command: Vec<String>,
    pub stdout: String,
    pub stderr: String,
}

#[derive(Debug, Clone)]
pub enum RnaStructureError {
    UnsupportedBiotype {
        molecule_type: Option<String>,
    },
    EmptySequence,
    ToolNotFound {
        executable: String,
    },
    ToolFailed {
        executable: String,
        args: Vec<String>,
        status: Option<i32>,
        stdout: String,
        stderr: String,
    },
    Io {
        message: String,
    },
}

impl fmt::Display for RnaStructureError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::UnsupportedBiotype { molecule_type } => write!(
                f,
                "RNA structure analysis requires single-stranded RNA biotype, found {:?}",
                molecule_type
            ),
            Self::EmptySequence => write!(f, "RNA sequence is empty"),
            Self::ToolNotFound { executable } => write!(
                f,
                "Could not find rnapkin executable '{}'. Install rnapkin or set {}",
                executable, RNAPKIN_ENV_BIN
            ),
            Self::ToolFailed {
                executable,
                args,
                status,
                stdout,
                stderr,
            } => write!(
                f,
                "rnapkin command failed: {} {} (status={:?}, stdout='{}', stderr='{}')",
                executable,
                args.join(" "),
                status,
                stdout.trim(),
                stderr.trim()
            ),
            Self::Io { message } => write!(f, "{message}"),
        }
    }
}

impl std::error::Error for RnaStructureError {}

pub fn is_single_stranded_rna(dna: &DNAsequence) -> bool {
    matches!(
        dna.molecule_type(),
        Some(v) if v.eq_ignore_ascii_case("RNA") || v.eq_ignore_ascii_case("ssRNA")
    )
}

fn rnapkin_executable() -> String {
    std::env::var(RNAPKIN_ENV_BIN)
        .ok()
        .map(|v| v.trim().to_string())
        .filter(|v| !v.is_empty())
        .unwrap_or_else(|| DEFAULT_RNAPKIN_BIN.to_string())
}

fn normalized_rna_sequence(dna: &DNAsequence) -> String {
    dna.get_forward_string()
        .chars()
        .filter(|c| !c.is_ascii_whitespace())
        .map(|c| match c.to_ascii_uppercase() {
            'T' => 'U',
            other => other,
        })
        .collect()
}

fn run_rnapkin(executable: &str, args: &[String]) -> Result<Output, RnaStructureError> {
    Command::new(executable).args(args).output().map_err(|e| {
        if e.kind() == ErrorKind::NotFound {
            RnaStructureError::ToolNotFound {
                executable: executable.to_string(),
            }
        } else {
            RnaStructureError::Io {
                message: format!(
                    "Could not run rnapkin executable '{}' with args [{}]: {}",
                    executable,
                    args.join(" "),
                    e
                ),
            }
        }
    })
}

fn ensure_rna_input(dna: &DNAsequence) -> Result<String, RnaStructureError> {
    if !is_single_stranded_rna(dna) {
        return Err(RnaStructureError::UnsupportedBiotype {
            molecule_type: dna.molecule_type().map(ToString::to_string),
        });
    }
    let sequence = normalized_rna_sequence(dna);
    if sequence.is_empty() {
        return Err(RnaStructureError::EmptySequence);
    }
    Ok(sequence)
}

pub fn inspect_text(dna: &DNAsequence) -> Result<RnaStructureTextReport, RnaStructureError> {
    let sequence = ensure_rna_input(dna)?;
    let executable = rnapkin_executable();
    let args = vec!["-v".to_string(), "-p".to_string(), sequence];
    let output = run_rnapkin(&executable, &args)?;
    let stdout = String::from_utf8_lossy(&output.stdout).to_string();
    let stderr = String::from_utf8_lossy(&output.stderr).to_string();

    if !output.status.success() {
        return Err(RnaStructureError::ToolFailed {
            executable,
            args,
            status: output.status.code(),
            stdout,
            stderr,
        });
    }

    Ok(RnaStructureTextReport {
        tool: "rnapkin".to_string(),
        executable,
        sequence_length: dna.len(),
        command: args,
        stdout,
        stderr,
    })
}

pub fn render_svg(
    dna: &DNAsequence,
    output_path: &str,
) -> Result<RnaStructureSvgReport, RnaStructureError> {
    let sequence = ensure_rna_input(dna)?;
    let output_path = output_path.trim();
    if output_path.is_empty() {
        return Err(RnaStructureError::Io {
            message: "RNA structure output path is empty".to_string(),
        });
    }

    if let Some(parent) = Path::new(output_path).parent() {
        if !parent.as_os_str().is_empty() {
            std::fs::create_dir_all(parent).map_err(|e| RnaStructureError::Io {
                message: format!("Could not create RNA structure output directory: {e}"),
            })?;
        }
    }

    let executable = rnapkin_executable();
    let args = vec![sequence, output_path.to_string()];
    let output = run_rnapkin(&executable, &args)?;
    let stdout = String::from_utf8_lossy(&output.stdout).to_string();
    let stderr = String::from_utf8_lossy(&output.stderr).to_string();

    if !output.status.success() {
        return Err(RnaStructureError::ToolFailed {
            executable,
            args,
            status: output.status.code(),
            stdout,
            stderr,
        });
    }

    if !Path::new(output_path).exists() {
        return Err(RnaStructureError::Io {
            message: format!(
                "rnapkin reported success but did not produce output file '{}'",
                output_path
            ),
        });
    }

    Ok(RnaStructureSvgReport {
        tool: "rnapkin".to_string(),
        executable,
        sequence_length: dna.len(),
        output_path: output_path.to_string(),
        command: args,
        stdout,
        stderr,
    })
}
