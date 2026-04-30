//! Compatibility shim for glossary-driven shell help rendering.
//!
//! The implementation now lives in `gentle-shell`; this module preserves the
//! existing root-crate import path while the shell crate extraction is staged.

pub use gentle_shell::{
    HelpOutputFormat, shell_help_json, shell_help_markdown, shell_help_markdown_from_glossary_json,
    shell_help_text, shell_topic_help_json, shell_topic_help_markdown, shell_topic_help_text,
};
