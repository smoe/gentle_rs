//! Shared shell tokenization utilities.
//!
//! Full typed command parsing is still staged in the root compatibility layer
//! until the engine crate extraction is complete. The quoting/escaping rules
//! live here first so every adapter can share the same low-level shell syntax.

/// Split one raw shell line using GENtle's shared quoting/escaping rules.
pub fn split_shell_words(line: &str) -> Result<Vec<String>, String> {
    #[derive(Clone, Copy, PartialEq, Eq)]
    enum Mode {
        Normal,
        SingleQuoted,
        DoubleQuoted,
    }

    let mut out = Vec::new();
    let mut current = String::new();
    let mut mode = Mode::Normal;
    let mut chars = line.chars().peekable();

    while let Some(ch) = chars.next() {
        match mode {
            Mode::Normal => match ch {
                '\'' => mode = Mode::SingleQuoted,
                '"' => mode = Mode::DoubleQuoted,
                '\\' => {
                    if let Some(next) = chars.next() {
                        current.push(next);
                    }
                }
                c if c.is_whitespace() => {
                    if !current.is_empty() {
                        out.push(current.clone());
                        current.clear();
                    }
                }
                _ => current.push(ch),
            },
            Mode::SingleQuoted => {
                if ch == '\'' {
                    mode = Mode::Normal;
                } else {
                    current.push(ch);
                }
            }
            Mode::DoubleQuoted => {
                if ch == '"' {
                    mode = Mode::Normal;
                } else if ch == '\\' {
                    if let Some(next) = chars.next() {
                        current.push(next);
                    }
                } else {
                    current.push(ch);
                }
            }
        }
    }

    if mode != Mode::Normal {
        return Err("Unterminated quoted string in shell command".to_string());
    }
    if !current.is_empty() {
        out.push(current);
    }
    if out.is_empty() {
        return Err("Empty shell command".to_string());
    }
    Ok(out)
}

#[cfg(test)]
mod tests {
    use super::split_shell_words;

    #[test]
    fn split_shell_words_preserves_quoted_payloads() {
        let words = split_shell_words("rna-reads inspect-alignments report --search 'TP73 alpha'")
            .expect("split shell words");
        assert_eq!(
            words,
            vec![
                "rna-reads",
                "inspect-alignments",
                "report",
                "--search",
                "TP73 alpha"
            ]
        );
    }

    #[test]
    fn split_shell_words_handles_double_quote_escapes() {
        let words = split_shell_words(r#"ui focus pcr-design --filter "TP73 \"alpha\"""#)
            .expect("split shell words");
        assert_eq!(
            words,
            vec!["ui", "focus", "pcr-design", "--filter", "TP73 \"alpha\""]
        );
    }

    #[test]
    fn split_shell_words_rejects_unclosed_quotes() {
        let err = split_shell_words("primers design 'unterminated")
            .expect_err("unterminated quote should fail");
        assert!(err.contains("Unterminated quoted string"));
    }
}
