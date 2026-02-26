//! Generated/derived shell-help documentation helpers.

use serde::Deserialize;
use serde_json::{Value, json};
use std::sync::OnceLock;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum HelpOutputFormat {
    Text,
    Json,
    Markdown,
}

impl HelpOutputFormat {
    pub fn parse(raw: &str) -> Result<Self, String> {
        match raw.trim().to_ascii_lowercase().as_str() {
            "text" => Ok(Self::Text),
            "json" => Ok(Self::Json),
            "markdown" | "md" => Ok(Self::Markdown),
            other => Err(format!(
                "Unsupported help format '{other}' (expected text|json|markdown)"
            )),
        }
    }

    pub fn as_str(self) -> &'static str {
        match self {
            Self::Text => "text",
            Self::Json => "json",
            Self::Markdown => "markdown",
        }
    }
}

#[derive(Debug, Clone, Deserialize)]
struct ShellGlossary {
    schema: String,
    #[serde(default)]
    interfaces: Vec<String>,
    #[serde(default)]
    commands: Vec<ShellCommandDoc>,
}

#[derive(Debug, Clone, Deserialize)]
struct ShellCommandDoc {
    path: String,
    usage: String,
    summary: String,
    #[serde(default)]
    interfaces: Vec<String>,
    #[serde(default)]
    engine_operations: Vec<String>,
    #[serde(default)]
    aliases: Vec<String>,
}

static GLOSSARY: OnceLock<Result<ShellGlossary, String>> = OnceLock::new();

fn parse_glossary(raw: &str) -> Result<ShellGlossary, String> {
    serde_json::from_str::<ShellGlossary>(raw)
        .map_err(|e| format!("Could not parse docs/glossary.json: {e}"))
}

fn glossary() -> Result<&'static ShellGlossary, String> {
    match GLOSSARY.get_or_init(|| {
        let raw = include_str!("../docs/glossary.json");
        parse_glossary(raw)
    }) {
        Ok(glossary) => Ok(glossary),
        Err(err) => Err(err.clone()),
    }
}

fn normalize_interface_filter(raw: Option<&str>) -> Result<Option<String>, String> {
    let Some(raw) = raw else {
        return Ok(None);
    };
    let normalized = raw.trim().to_ascii_lowercase();
    if normalized.is_empty() || normalized == "all" {
        return Ok(None);
    }
    if normalized == "mcp" {
        // MCP currently reuses the shared shell command surface.
        return Ok(Some("cli-shell".to_string()));
    }
    let valid = ["cli-direct", "cli-shell", "gui-shell", "js", "lua"];
    if valid.contains(&normalized.as_str()) {
        Ok(Some(normalized))
    } else {
        Err(format!(
            "Unsupported --interface '{}' (expected all|cli-direct|cli-shell|gui-shell|js|lua|mcp)",
            raw
        ))
    }
}

fn path_tokens(path: &str) -> Vec<String> {
    path.split_whitespace()
        .map(|t| t.to_ascii_lowercase())
        .collect()
}

fn topic_tokens(topic: &[String]) -> Vec<String> {
    topic
        .iter()
        .map(|t| t.trim().to_ascii_lowercase())
        .filter(|t| !t.is_empty())
        .collect()
}

fn is_prefix_path(path: &str, topic: &[String]) -> bool {
    let path = path_tokens(path);
    if topic.len() < path.len() {
        return false;
    }
    path.iter().zip(topic.iter()).all(|(a, b)| a == b)
}

fn docs_for_interface(
    interface_filter: Option<&str>,
) -> Result<Vec<&'static ShellCommandDoc>, String> {
    let glossary = glossary()?;
    docs_for_interface_from_glossary(glossary, interface_filter)
}

fn docs_for_interface_from_glossary<'a>(
    glossary: &'a ShellGlossary,
    interface_filter: Option<&str>,
) -> Result<Vec<&'a ShellCommandDoc>, String> {
    let interface_filter = normalize_interface_filter(interface_filter)?;
    Ok(glossary
        .commands
        .iter()
        .filter(|doc| {
            interface_filter
                .as_ref()
                .map_or(true, |f| doc.interfaces.iter().any(|iface| iface == f))
        })
        .collect())
}

fn doc_record(doc: &ShellCommandDoc) -> Value {
    json!({
        "path": doc.path,
        "usage": doc.usage,
        "summary": doc.summary,
        "interfaces": doc.interfaces,
        "engine_operations": doc.engine_operations,
        "aliases": doc.aliases
    })
}

fn find_doc_for_topic<'a>(
    docs: &'a [&'static ShellCommandDoc],
    topic: &[String],
) -> Option<&'a ShellCommandDoc> {
    let topic = topic_tokens(topic);
    docs.iter()
        .filter_map(|doc| {
            let mut longest = None::<usize>;
            if is_prefix_path(&doc.path, &topic) {
                longest = Some(path_tokens(&doc.path).len());
            }
            for alias in &doc.aliases {
                if is_prefix_path(alias, &topic) {
                    let score = path_tokens(alias).len();
                    if longest.map_or(true, |v| score > v) {
                        longest = Some(score);
                    }
                }
            }
            longest.map(|score| (score, *doc))
        })
        .max_by_key(|(score, _)| *score)
        .map(|(_, doc)| doc)
}

fn topic_not_found_message(topic: &[String], docs: &[&ShellCommandDoc]) -> String {
    let topic_joined = topic.join(" ");
    let first = topic
        .first()
        .map(|v| v.to_ascii_lowercase())
        .unwrap_or_default();
    let mut suggestions = docs
        .iter()
        .map(|doc| doc.path.as_str())
        .filter(|path| first.is_empty() || path.starts_with(&first))
        .take(8)
        .collect::<Vec<_>>();
    if suggestions.is_empty() {
        suggestions = docs
            .iter()
            .map(|doc| doc.path.as_str())
            .take(8)
            .collect::<Vec<_>>();
    }
    format!(
        "Unknown help topic '{}'. Try one of: {}",
        topic_joined,
        suggestions.join(", ")
    )
}

fn render_topic_text(doc: &ShellCommandDoc) -> String {
    let mut out = String::new();
    out.push_str("GENtle command help\n");
    out.push_str(&format!("Path: {}\n", doc.path));
    out.push_str(&format!("Usage: {}\n", doc.usage));
    out.push_str(&format!("Summary: {}\n", doc.summary));
    out.push_str(&format!("Interfaces: {}\n", doc.interfaces.join(", ")));
    if !doc.engine_operations.is_empty() {
        out.push_str(&format!(
            "Engine operations: {}\n",
            doc.engine_operations.join(", ")
        ));
    }
    if !doc.aliases.is_empty() {
        out.push_str(&format!("Aliases: {}\n", doc.aliases.join(", ")));
    }
    out
}

fn render_topic_markdown(doc: &ShellCommandDoc) -> String {
    let mut out = String::new();
    out.push_str(&format!("## `{}`\n", doc.path));
    out.push_str(&format!("- Usage: `{}`\n", doc.usage));
    out.push_str(&format!("- Summary: {}\n", doc.summary));
    out.push_str(&format!(
        "- Interfaces: `{}`\n",
        doc.interfaces.join("`, `")
    ));
    if !doc.engine_operations.is_empty() {
        out.push_str(&format!(
            "- Engine operations: `{}`\n",
            doc.engine_operations.join("`, `")
        ));
    }
    if !doc.aliases.is_empty() {
        out.push_str(&format!("- Aliases: `{}`\n", doc.aliases.join("`, `")));
    }
    out
}

pub fn shell_help_text(interface_filter: Option<&str>) -> Result<String, String> {
    let docs = docs_for_interface(interface_filter)?;
    let mut out = String::from("GENtle Shell commands:\n");
    for doc in docs {
        out.push_str("- ");
        out.push_str(&doc.usage);
        out.push_str("\n  ");
        out.push_str(&doc.summary);
        out.push('\n');
    }
    out.push_str(
        "Use `help COMMAND ...` for command-specific help.\n\
Use `help [--format json|markdown]` to export machine-readable docs.\n\
Use `--interface` to filter (`all|cli-direct|cli-shell|gui-shell|js|lua|mcp`).",
    );
    Ok(out)
}

pub fn shell_help_markdown(interface_filter: Option<&str>) -> Result<String, String> {
    let docs = docs_for_interface(interface_filter)?;
    let mut out = String::from("# GENtle Shell Command Reference\n\n");
    for doc in docs {
        out.push_str(&render_topic_markdown(doc));
        out.push('\n');
    }
    Ok(out)
}

pub fn shell_help_markdown_from_glossary_json(
    raw: &str,
    interface_filter: Option<&str>,
) -> Result<String, String> {
    let glossary = parse_glossary(raw)?;
    let docs = docs_for_interface_from_glossary(&glossary, interface_filter)?;
    let mut out = String::from("# GENtle Shell Command Reference\n\n");
    for doc in docs {
        out.push_str(&render_topic_markdown(doc));
        out.push('\n');
    }
    Ok(out)
}

pub fn shell_help_json(interface_filter: Option<&str>) -> Result<Value, String> {
    let glossary = glossary()?;
    let docs = docs_for_interface(interface_filter)?;
    Ok(json!({
        "schema": "gentle.shell_help_catalog.v1",
        "source_schema": glossary.schema,
        "known_interfaces": glossary.interfaces,
        "interface_filter": normalize_interface_filter(interface_filter)?.unwrap_or_else(|| "all".to_string()),
        "command_count": docs.len(),
        "commands": docs.iter().map(|doc| doc_record(doc)).collect::<Vec<_>>()
    }))
}

pub fn shell_topic_help_text(
    topic: &[String],
    interface_filter: Option<&str>,
) -> Result<String, String> {
    let docs = docs_for_interface(interface_filter)?;
    let doc =
        find_doc_for_topic(&docs, topic).ok_or_else(|| topic_not_found_message(topic, &docs))?;
    Ok(render_topic_text(doc))
}

pub fn shell_topic_help_markdown(
    topic: &[String],
    interface_filter: Option<&str>,
) -> Result<String, String> {
    let docs = docs_for_interface(interface_filter)?;
    let doc =
        find_doc_for_topic(&docs, topic).ok_or_else(|| topic_not_found_message(topic, &docs))?;
    Ok(render_topic_markdown(doc))
}

pub fn shell_topic_help_json(
    topic: &[String],
    interface_filter: Option<&str>,
) -> Result<Value, String> {
    let docs = docs_for_interface(interface_filter)?;
    let doc =
        find_doc_for_topic(&docs, topic).ok_or_else(|| topic_not_found_message(topic, &docs))?;
    Ok(json!({
        "schema": "gentle.shell_help_topic.v1",
        "topic": topic.join(" "),
        "doc": doc_record(doc)
    }))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn interface_filter_accepts_mcp_alias() {
        let normalized = normalize_interface_filter(Some("mcp")).expect("normalize filter");
        assert_eq!(normalized.as_deref(), Some("cli-shell"));
    }

    #[test]
    fn shell_help_text_lists_mcp_interface_value() {
        let text = shell_help_text(None).expect("render help");
        assert!(text.contains("|mcp"));
    }
}
