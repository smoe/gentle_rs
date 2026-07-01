//! Small argument helpers for the direct `gentle_cli` command surface.

use std::str::FromStr;

pub(super) fn required_value(
    args: &[String],
    idx: usize,
    flag: &str,
    expected: &str,
    context: &str,
) -> Result<String, String> {
    if idx + 1 >= args.len() {
        return Err(format!("Missing {expected} after {flag} for {context}"));
    }
    Ok(args[idx + 1].clone())
}

pub(super) fn parse_value<T>(raw: &str, flag: &str, context: Option<&str>) -> Result<T, String>
where
    T: FromStr,
    T::Err: std::fmt::Display,
{
    raw.parse::<T>().map_err(|e| match context {
        Some(context) => format!("Invalid {flag} value '{raw}' for {context}: {e}"),
        None => format!("Invalid {flag} value '{raw}': {e}"),
    })
}

pub(super) fn parse_required_value<T>(
    args: &[String],
    idx: usize,
    flag: &str,
    expected: &str,
    context: &str,
    invalid_context: Option<&str>,
) -> Result<T, String>
where
    T: FromStr,
    T::Err: std::fmt::Display,
{
    let raw = required_value(args, idx, flag, expected, context)?;
    parse_value(&raw, flag, invalid_context)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn argv(values: &[&str]) -> Vec<String> {
        values.iter().map(|value| value.to_string()).collect()
    }

    #[test]
    fn required_value_reports_missing_with_context() {
        let args = argv(&["gentle_cli", "genomes", "list", "--catalog"]);
        let err = required_value(&args, 3, "--catalog", "PATH", "genomes list")
            .expect_err("missing flag value should fail");
        assert_eq!(err, "Missing PATH after --catalog for genomes list");
    }

    #[test]
    fn required_value_returns_following_token() {
        let args = argv(&["gentle_cli", "genomes", "list", "--filter", "tp73"]);
        let value = required_value(&args, 3, "--filter", "TEXT", "genomes list")
            .expect("flag value should parse");
        assert_eq!(value, "tp73");
    }

    #[test]
    fn parse_value_reports_legacy_invalid_number_without_context() {
        let err =
            parse_value::<usize>("abc", "--limit", None).expect_err("invalid integer should fail");
        assert!(err.starts_with("Invalid --limit value 'abc': "));
    }

    #[test]
    fn parse_required_value_reports_invalid_number_with_context() {
        let args = argv(&["gentle_cli", "genomes", "prepare", "--timeout-secs", "abc"]);
        let err = parse_required_value::<u64>(
            &args,
            3,
            "--timeout-secs",
            "N",
            "genomes prepare",
            Some("genomes prepare"),
        )
        .expect_err("invalid integer should fail");
        assert!(err.starts_with("Invalid --timeout-secs value 'abc' for genomes prepare: "));
    }
}
