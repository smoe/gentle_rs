//! Direct `candidates` CLI command handling kept out of `run()`.
//!
//! `candidates ...` invocations are first-class direct CLI routes that wrap the
//! shared shell parser and executor. This keeps candidate command behavior
//! engine-owned while the CLI adapter is split into smaller families.

use super::*;

pub(super) const USAGE: &str = "\
  gentle_cli [--state PATH|--project PATH] candidates list
  gentle_cli [--state PATH|--project PATH] candidates delete SET_NAME
  gentle_cli [--state PATH|--project PATH] candidates generate SET_NAME SEQ_ID --length N [--step N] [--feature-kind KIND] [--feature-label-regex REGEX] [--max-distance N] [--feature-geometry feature_span|feature_parts|feature_boundaries] [--feature-boundary any|five_prime|three_prime|start|end] [--strand-relation any|same|opposite] [--limit N]
  gentle_cli [--state PATH|--project PATH] candidates generate-between-anchors SET_NAME SEQ_ID --length N (--anchor-a-pos N|--anchor-a-json JSON) (--anchor-b-pos N|--anchor-b-json JSON) [--step N] [--limit N]
  gentle_cli [--state PATH|--project PATH] candidates show SET_NAME [--limit N] [--offset N]
  gentle_cli [--state PATH|--project PATH] candidates metrics SET_NAME
  gentle_cli [--state PATH|--project PATH] candidates score SET_NAME METRIC_NAME EXPRESSION
  gentle_cli [--state PATH|--project PATH] candidates score-distance SET_NAME METRIC_NAME [--feature-kind KIND] [--feature-label-regex REGEX] [--feature-geometry feature_span|feature_parts|feature_boundaries] [--feature-boundary any|five_prime|three_prime|start|end] [--strand-relation any|same|opposite]
  gentle_cli [--state PATH|--project PATH] candidates score-weighted SET_NAME METRIC_NAME --term METRIC:WEIGHT[:max|min] [--term ...] [--normalize|--no-normalize]
  gentle_cli [--state PATH|--project PATH] candidates top-k INPUT_SET OUTPUT_SET --metric METRIC_NAME --k N [--direction max|min] [--tie-break seq_start_end|seq_end_start|length_ascending|length_descending|sequence_lexicographic]
  gentle_cli [--state PATH|--project PATH] candidates pareto INPUT_SET OUTPUT_SET --objective METRIC[:max|min] [--objective ...] [--max-candidates N] [--tie-break seq_start_end|seq_end_start|length_ascending|length_descending|sequence_lexicographic]
  gentle_cli [--state PATH|--project PATH] candidates filter INPUT_SET OUTPUT_SET --metric METRIC_NAME [--min N] [--max N] [--min-quantile Q] [--max-quantile Q]
  gentle_cli [--state PATH|--project PATH] candidates set-op union|intersect|subtract LEFT_SET RIGHT_SET OUTPUT_SET
  gentle_cli [--state PATH|--project PATH] candidates macro SCRIPT_OR_@FILE
  gentle_cli [--state PATH|--project PATH] candidates template-list
  gentle_cli [--state PATH|--project PATH] candidates template-show TEMPLATE_NAME
  gentle_cli [--state PATH|--project PATH] candidates template-put TEMPLATE_NAME (--script SCRIPT_OR_@FILE|--file PATH) [--description TEXT] [--details-url URL] [--param NAME|NAME=DEFAULT ...]
  gentle_cli [--state PATH|--project PATH] candidates template-delete TEMPLATE_NAME
  gentle_cli [--state PATH|--project PATH] candidates template-run TEMPLATE_NAME [--bind KEY=VALUE ...] [--transactional]

";

pub(super) fn is_candidates_shell_command(shell_command: &ShellCommand) -> bool {
    matches!(
        shell_command,
        ShellCommand::CandidatesList
            | ShellCommand::CandidatesDelete { .. }
            | ShellCommand::CandidatesGenerate { .. }
            | ShellCommand::CandidatesGenerateBetweenAnchors { .. }
            | ShellCommand::CandidatesShow { .. }
            | ShellCommand::CandidatesMetrics { .. }
            | ShellCommand::CandidatesScoreExpression { .. }
            | ShellCommand::CandidatesScoreDistance { .. }
            | ShellCommand::CandidatesScoreWeightedObjective { .. }
            | ShellCommand::CandidatesTopK { .. }
            | ShellCommand::CandidatesParetoFrontier { .. }
            | ShellCommand::CandidatesFilter { .. }
            | ShellCommand::CandidatesSetOp { .. }
            | ShellCommand::CandidatesMacro { .. }
            | ShellCommand::CandidatesTemplateList
            | ShellCommand::CandidatesTemplateShow { .. }
            | ShellCommand::CandidatesTemplateUpsert { .. }
            | ShellCommand::CandidatesTemplateDelete { .. }
            | ShellCommand::CandidatesTemplateRun { .. }
    )
}

pub(super) fn handle_candidates_family(
    args: &[String],
    cmd_idx: usize,
    state_path: &str,
    shell_options: &ShellExecutionOptions,
) -> Result<(), String> {
    if args.len() <= cmd_idx + 1 {
        usage();
        return Err("candidates requires a subcommand".to_string());
    }
    let tokens = args[cmd_idx..].to_vec();
    let shell_command = parse_shell_tokens(&tokens)?;
    if !is_candidates_shell_command(&shell_command) {
        return Err("Expected a candidates subcommand".to_string());
    }
    let mut engine = GentleEngine::from_state(load_state(state_path)?);
    let run = execute_shell_command_with_options(&mut engine, &shell_command, shell_options)?;
    if run.state_changed {
        engine
            .state()
            .save_to_path(state_path)
            .map_err(|e| e.to_string())?;
    }
    print_json(&run.output)
}

#[cfg(test)]
mod tests {
    use super::*;
    use tempfile::tempdir;

    fn argv(values: &[&str]) -> Vec<String> {
        values.iter().map(|value| value.to_string()).collect()
    }

    #[test]
    fn usage_block_keeps_candidate_routes() {
        assert!(USAGE.contains("candidates list"));
        assert!(USAGE.contains("candidates generate SET_NAME"));
        assert!(USAGE.contains("candidates template-run"));
    }

    #[test]
    fn candidates_missing_subcommand_reports_legacy_message() {
        let args = argv(&["gentle_cli", "candidates"]);
        let options = ShellExecutionOptions {
            allow_screenshots: false,
            allow_agent_commands: true,
            progress_callback: None,
        };
        let err = handle_candidates_family(&args, 1, ".gentle_state.json", &options)
            .expect_err("missing subcommand should fail");
        assert_eq!(err, "candidates requires a subcommand");
    }

    #[test]
    fn candidates_list_runs_through_direct_wrapper() {
        let td = tempdir().expect("tempdir");
        let state_path = td.path().join("state.json");
        let args = argv(&["gentle_cli", "candidates", "list"]);
        let options = ShellExecutionOptions {
            allow_screenshots: false,
            allow_agent_commands: true,
            progress_callback: None,
        };
        handle_candidates_family(&args, 1, &state_path.to_string_lossy(), &options)
            .expect("candidates list should execute");
        assert!(
            !state_path.exists(),
            "read-only candidates list should not persist state"
        );
    }

    #[test]
    fn candidate_command_filter_accepts_template_list() {
        assert!(is_candidates_shell_command(
            &ShellCommand::CandidatesTemplateList
        ));
    }

    #[test]
    fn candidate_command_filter_rejects_other_shell_commands() {
        assert!(!is_candidates_shell_command(
            &ShellCommand::ServicesProvidersList
        ));
    }
}
