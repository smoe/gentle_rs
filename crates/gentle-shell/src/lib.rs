//! Shared GENtle shell grammar and executor.
//!
//! This crate now owns the glossary-driven shell help rendering layer, and it
//! remains the intended home for the shared textual command grammar/executor
//! reused by GUI Shell, CLI shell, and shell-like adapters.

mod command;
mod parser;
mod shell_docs;
mod ui_intent;

pub use command::{
    BatchEmitMode, BatchManifestDelimiter, BatchStateMode, CandidateSetOperator, ShellRunResult,
};
pub use parser::split_shell_words;
pub use shell_docs::{
    HelpOutputFormat, shell_help_json, shell_help_markdown, shell_help_markdown_from_glossary_json,
    shell_help_text, shell_topic_help_json, shell_topic_help_markdown, shell_topic_help_text,
};
pub use ui_intent::{
    UiIntentAction, UiIntentTarget, UiIntentTargetCatalogRow, ui_intent_target_catalog,
};
