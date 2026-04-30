//! RNA structure wrappers and tool integration glue.

use crate::dna_sequence::DNAsequence;
use serde::{Deserialize, Serialize};
use std::{
    fmt,
    io::{ErrorKind, Write},
    path::Path,
    process::{Command, Output, Stdio},
};

pub const DEFAULT_RNAFOLD_BIN: &str = "RNAfold";
pub const RNAFOLD_ENV_BIN: &str = "GENTLE_RNAFOLD_BIN";
pub const DEFAULT_RNAPKIN_BIN: &str = "rnapkin";
pub const RNAPKIN_ENV_BIN: &str = "GENTLE_RNAPKIN_BIN";

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct RnaStructureTextReport {
    pub tool: String,
    pub executable: String,
    pub sequence_length: usize,
    pub command: Vec<String>,
    pub stdout: String,
    pub stderr: String,
    #[serde(default)]
    pub structure: String,
    #[serde(default)]
    pub mfe_kcal_per_mol: Option<f64>,
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
    #[serde(default)]
    pub fold_executable: String,
    #[serde(default)]
    pub fold_command: Vec<String>,
    #[serde(default)]
    pub fold_stdout: String,
    #[serde(default)]
    pub fold_stderr: String,
    #[serde(default)]
    pub structure: String,
    #[serde(default)]
    pub mfe_kcal_per_mol: Option<f64>,
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
                "Could not find RNA structure executable '{}'. Install RNAfold/rnapkin or set {} / {}",
                executable, RNAFOLD_ENV_BIN, RNAPKIN_ENV_BIN
            ),
            Self::ToolFailed {
                executable,
                args,
                status,
                stdout,
                stderr,
            } => write!(
                f,
                "RNA structure command failed: {} {} (status={:?}, stdout='{}', stderr='{}')",
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
    crate::tool_overrides::resolve_tool_executable(RNAPKIN_ENV_BIN, DEFAULT_RNAPKIN_BIN)
}

fn rnafold_executable() -> String {
    crate::tool_overrides::resolve_tool_executable(RNAFOLD_ENV_BIN, DEFAULT_RNAFOLD_BIN)
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

fn run_tool(
    executable: &str,
    args: &[String],
    stdin_text: Option<&str>,
) -> Result<Output, RnaStructureError> {
    if let Some(stdin_text) = stdin_text {
        let mut child = Command::new(executable)
            .args(args)
            .stdin(Stdio::piped())
            .stdout(Stdio::piped())
            .stderr(Stdio::piped())
            .spawn()
            .map_err(|e| map_spawn_error(executable, args, e))?;
        if let Some(mut stdin) = child.stdin.take() {
            stdin
                .write_all(stdin_text.as_bytes())
                .map_err(|e| RnaStructureError::Io {
                    message: format!(
                        "Could not write RNA structure input to '{}' with args [{}]: {}",
                        executable,
                        args.join(" "),
                        e
                    ),
                })?;
        }
        return child.wait_with_output().map_err(|e| RnaStructureError::Io {
            message: format!(
                "Could not collect RNA structure output from '{}' with args [{}]: {}",
                executable,
                args.join(" "),
                e
            ),
        });
    }

    Command::new(executable)
        .args(args)
        .output()
        .map_err(|e| map_spawn_error(executable, args, e))
}

fn map_spawn_error(executable: &str, args: &[String], e: std::io::Error) -> RnaStructureError {
    if e.kind() == ErrorKind::NotFound {
        RnaStructureError::ToolNotFound {
            executable: executable.to_string(),
        }
    } else {
        RnaStructureError::Io {
            message: format!(
                "Could not run RNA structure executable '{}' with args [{}]: {}",
                executable,
                args.join(" "),
                e
            ),
        }
    }
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

#[derive(Debug, Clone)]
struct RnaFoldResult {
    executable: String,
    args: Vec<String>,
    stdout: String,
    stderr: String,
    structure: String,
    mfe_kcal_per_mol: Option<f64>,
}

fn parse_rnafold_structure(stdout: &str) -> Result<(String, Option<f64>), RnaStructureError> {
    let mut lines = stdout.lines().filter(|line| !line.trim().is_empty());
    let _sequence_line = lines.next().ok_or_else(|| RnaStructureError::Io {
        message: "RNAfold did not emit a folded sequence line".to_string(),
    })?;
    let structure_line = lines.next().ok_or_else(|| RnaStructureError::Io {
        message: "RNAfold did not emit a dot-bracket structure line".to_string(),
    })?;
    let structure = structure_line
        .split_whitespace()
        .next()
        .unwrap_or_default()
        .trim()
        .to_string();
    if structure.is_empty() {
        return Err(RnaStructureError::Io {
            message: "RNAfold emitted an empty dot-bracket structure".to_string(),
        });
    }
    let rest = structure_line
        .strip_prefix(&structure)
        .unwrap_or_default()
        .trim();
    let mfe_kcal_per_mol = rest
        .strip_prefix('(')
        .and_then(|v| v.strip_suffix(')'))
        .and_then(|v| v.trim().parse::<f64>().ok());
    Ok((structure, mfe_kcal_per_mol))
}

fn fold_rna_sequence(sequence: &str) -> Result<RnaFoldResult, RnaStructureError> {
    let executable = rnafold_executable();
    let args = vec!["--noPS".to_string()];
    let input = format!("{sequence}\n");
    let output = run_tool(&executable, &args, Some(&input))?;
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

    let (structure, mfe_kcal_per_mol) = parse_rnafold_structure(&stdout)?;
    Ok(RnaFoldResult {
        executable,
        args,
        stdout,
        stderr,
        structure,
        mfe_kcal_per_mol,
    })
}

pub fn inspect_text(dna: &DNAsequence) -> Result<RnaStructureTextReport, RnaStructureError> {
    let sequence = ensure_rna_input(dna)?;
    let fold = fold_rna_sequence(&sequence)?;

    Ok(RnaStructureTextReport {
        tool: "RNAfold".to_string(),
        executable: fold.executable,
        sequence_length: dna.len(),
        command: fold.args,
        stdout: fold.stdout,
        stderr: fold.stderr,
        structure: fold.structure,
        mfe_kcal_per_mol: fold.mfe_kcal_per_mol,
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

    let fold = fold_rna_sequence(&sequence)?;
    let mut rnapkin_input = tempfile::NamedTempFile::new().map_err(|e| RnaStructureError::Io {
        message: format!("Could not create temporary rnapkin input: {e}"),
    })?;
    writeln!(rnapkin_input, "{sequence}").map_err(|e| RnaStructureError::Io {
        message: format!("Could not write temporary rnapkin sequence: {e}"),
    })?;
    writeln!(rnapkin_input, "{}", fold.structure).map_err(|e| RnaStructureError::Io {
        message: format!("Could not write temporary rnapkin structure: {e}"),
    })?;
    rnapkin_input.flush().map_err(|e| RnaStructureError::Io {
        message: format!("Could not flush temporary rnapkin input: {e}"),
    })?;

    let executable = rnapkin_executable();
    let input_path = rnapkin_input.path().display().to_string();
    let args = vec!["-o".to_string(), output_path.to_string(), input_path];
    let output = run_tool(&executable, &args, None)?;
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
        fold_executable: fold.executable,
        fold_command: fold.args,
        fold_stdout: fold.stdout,
        fold_stderr: fold.stderr,
        structure: fold.structure,
        mfe_kcal_per_mol: fold.mfe_kcal_per_mol,
    })
}
