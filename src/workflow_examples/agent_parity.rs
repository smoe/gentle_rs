//! Non-executing tutorial command validation against the shared shell parser.

use super::{TutorialAgentExecutionIntent, TutorialAgentParityCase, TutorialSourceUnit};
use crate::engine_shell::{ShellCommand, parse_shell_line, shell_quote, split_shell_words};
use pulldown_cmark::{CodeBlockKind, Event, Parser, Tag, TagEnd};
use serde::{Deserialize, Serialize};
use std::{
    fs,
    path::{Path, PathBuf},
};

/// Author intent and conservative parser classification are counted separately.
#[derive(Debug, Clone, Default, Serialize, Deserialize, PartialEq, Eq)]
pub struct TutorialAgentParitySummary {
    pub tutorials: usize,
    pub cases: usize,
    pub declared_mutating: usize,
    pub parser_state_mutating: usize,
}

/// Stable code plus source identities, suitable for CLI and library consumers.
#[derive(Debug, Clone, Serialize, PartialEq, Eq)]
pub struct TutorialAgentParityFinding {
    pub code: String,
    pub tutorial_id: String,
    pub case_id: String,
    pub detail: String,
}

#[derive(Debug, Clone, Default, Serialize)]
pub struct TutorialAgentParityCheck {
    pub summary: TutorialAgentParitySummary,
    pub findings: Vec<TutorialAgentParityFinding>,
}

impl TutorialAgentParityCheck {
    /// Preserve tutorial-check's stderr/exit-1 convention, with one JSON finding
    /// per line so identifiers and embedded newlines remain unambiguous.
    pub fn failure_message(&self) -> Option<String> {
        if self.findings.is_empty() {
            return None;
        }
        let mut message = format!(
            "tutorial agent parity failed ({} findings):",
            self.findings.len()
        );
        for finding in &self.findings {
            message.push('\n');
            message.push_str(&serde_json::to_string(finding).expect("string-only finding"));
        }
        Some(message)
    }
}

type Diagnostic = (&'static str, String);

// Contracts use portable repository-relative paths. Check the original spelling
// before canonicalization so traversal is rejected even when it resolves inside.
fn validate_relative_path(raw: &str) -> Result<(), Diagnostic> {
    if raw.is_empty()
        || raw.contains(['\\', ':'])
        || Path::new(raw).is_absolute()
        || raw.split('/').any(|part| part == "..")
    {
        return Err((
            "unsafe_path",
            format!("Expected a portable repository-relative path: {raw:?}"),
        ));
    }
    Ok(())
}

fn bound_file(root: &Path, raw: &str) -> Result<PathBuf, Diagnostic> {
    validate_relative_path(raw)?;
    let path = root
        .join(raw)
        .canonicalize()
        .map_err(|error| ("file_unavailable", format!("{raw}: {error}")))?;
    if !path.starts_with(root) {
        return Err((
            "unsafe_path",
            format!("Path escapes the repository via a symlink: {raw}"),
        ));
    }
    if !path.is_file() {
        return Err(("file_unavailable", format!("Not a file: {raw}")));
    }
    Ok(path)
}

fn read_bound_file(root: &Path, raw: &str) -> Result<String, Diagnostic> {
    fs::read_to_string(bound_file(root, raw)?)
        .map_err(|error| ("file_unavailable", format!("{raw}: {error}")))
}

fn guide_json_blocks(guide: &str) -> Vec<serde_json::Value> {
    let mut blocks = Vec::new();
    let mut current = None;
    for event in Parser::new(guide) {
        match event {
            Event::Start(Tag::CodeBlock(CodeBlockKind::Fenced(info))) if info.trim() == "json" => {
                current = Some(String::new());
            }
            Event::Text(text) => {
                if let Some(block) = &mut current {
                    block.push_str(&text);
                }
            }
            Event::End(TagEnd::CodeBlock) => {
                if let Some(block) = current.take()
                    && let Ok(value) = serde_json::from_str(&block)
                {
                    blocks.push(value);
                }
            }
            _ => {}
        }
    }
    blocks
}

fn parser_command(
    case: &TutorialAgentParityCase,
    root: &Path,
    guide: &str,
) -> Result<String, Diagnostic> {
    let mut tokens =
        split_shell_words(&case.command).map_err(|error| ("parse_failed", error.to_string()))?;
    let mapping = case.parser_payload.as_ref();
    if let Some(mapping) = mapping {
        validate_relative_path(&mapping.file)?;
        let expected = format!("@{}", mapping.file);
        if tokens.iter().filter(|token| *token == &expected).count() != 1 {
            return Err((
                "payload_binding_invalid",
                format!("Expected exactly one {expected:?} token"),
            ));
        }
        let template = read_bound_file(root, &mapping.template)?;
        let value: serde_json::Value = serde_json::from_str(&template)
            .map_err(|error| ("payload_invalid_json", error.to_string()))?;
        if !guide_json_blocks(guide).contains(&value) {
            return Err((
                "payload_not_in_guide",
                format!(
                    "Template {} does not match a fenced JSON block in the guide",
                    mapping.template
                ),
            ));
        }
        // Use the exact validated JSON, not an invented execution-ready payload.
        for token in &mut tokens {
            if *token == expected {
                *token = template.clone();
            }
        }
    }
    for token in &mut tokens {
        if let Some(raw) = token.strip_prefix('@') {
            let path = bound_file(root, raw)?;
            let path = path
                .to_str()
                .ok_or_else(|| ("unsafe_path", "Payload path is not UTF-8".to_string()))?;
            *token = format!("@{path}");
        }
    }
    Ok(tokens
        .iter()
        .map(|token| shell_quote(token))
        .collect::<Vec<_>>()
        .join(" "))
}

/// Validate every declared contract, without an engine, model, approval, or
/// command execution. Callers load source units through the shared source loader.
/// Relative guide/template/`@file` paths are bound to `repo_root`, never the CWD.
pub fn check_tutorial_agent_parity(
    units: &[TutorialSourceUnit],
    repo_root: &Path,
) -> TutorialAgentParityCheck {
    let mut result = TutorialAgentParityCheck::default();
    let root = repo_root
        .canonicalize()
        .map_err(|error| ("repo_root_unavailable", error.to_string()));
    for unit in units {
        let Some(contract) = &unit.agent_parity else {
            continue;
        };
        result.summary.tutorials += 1;
        let guide = root.as_ref().map_err(Clone::clone).and_then(|root| {
            let path = unit
                .generated_chapter
                .as_ref()
                .and_then(|chapter| chapter.guide_path.as_deref())
                .or_else(|| unit.catalog.as_ref().map(|catalog| catalog.path.as_str()))
                .ok_or_else(|| ("guide_unavailable", "No guide path declared".to_string()))?;
            read_bound_file(root, path)
        });
        for case in &contract.cases {
            result.summary.cases += 1;
            result.summary.declared_mutating += usize::from(case.mutating);
            let mut record = |(code, detail): Diagnostic| {
                result.findings.push(TutorialAgentParityFinding {
                    code: code.to_string(),
                    tutorial_id: unit.id.clone(),
                    case_id: case.id.clone(),
                    detail,
                })
            };
            let guide = match &guide {
                Ok(guide) => guide,
                Err(error) => {
                    record(error.clone());
                    continue;
                }
            };
            if !guide.contains(&case.command) {
                record((
                    "guide_missing_command",
                    format!("Command not present verbatim: {}", case.command),
                ));
            }
            let command =
                match parser_command(case, root.as_ref().expect("guide needs root"), guide) {
                    Ok(command) => command,
                    Err(error) => {
                        record(error);
                        continue;
                    }
                };
            match parse_shell_line(&command) {
                Err(error) => record(("parse_failed", error.to_string())),
                Ok(parsed) => {
                    let mutating = parsed.is_state_mutating();
                    result.summary.parser_state_mutating += usize::from(mutating);
                    if case.mutating && !mutating {
                        record((
                            "declared_mutation_not_parser_mutating",
                            "Author declares mutation but parser does not".to_string(),
                        ));
                    }
                    if case.execution == TutorialAgentExecutionIntent::Auto
                        && (mutating
                            || matches!(
                                parsed,
                                ShellCommand::HistoryUndo | ShellCommand::HistoryRedo
                            ))
                    {
                        record(("auto_blocked_by_runtime", "Runtime requires explicit review for parser-classified mutations and undo/redo".to_string()));
                    }
                }
            }
        }
    }
    result
}

#[cfg(test)]
mod tests;
