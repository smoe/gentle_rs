//! Shared GENtle shell grammar and executor.
//!
//! This crate now owns the glossary-driven shell help rendering layer, and it
//! remains the intended home for the shared textual command grammar/executor
//! reused by GUI Shell, CLI shell, and shell-like adapters.

mod shell_docs;

pub use shell_docs::{
    HelpOutputFormat, shell_help_json, shell_help_markdown, shell_help_markdown_from_glossary_json,
    shell_help_text, shell_topic_help_json, shell_topic_help_markdown, shell_topic_help_text,
};
