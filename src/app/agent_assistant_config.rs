//! Agent Assistant defaults, prompt templates, and system-selection helpers.
//!
//! Keeping this small policy bundle outside `app.rs` makes the top-level GUI
//! coordinator less dense while preserving the same Agent Assistant behavior.

use std::env;

use crate::agent_bridge::{
    AGENT_CONNECT_TIMEOUT_SECS_ENV, AGENT_MAX_RESPONSE_BYTES_ENV, AGENT_MAX_RETRIES_ENV,
    AGENT_READ_TIMEOUT_SECS_ENV, AGENT_TIMEOUT_SECS_ENV, AgentSystemSpec, AgentSystemTransport,
    OPENAI_COMPAT_UNSPECIFIED_MODEL,
};

pub(super) const AGENT_PROMPT_TEMPLATE_DEFAULT_ID: &str = "structured";

pub(super) fn normalize_agent_model_name(raw: &str) -> Option<String> {
    let trimmed = raw.trim();
    if trimmed.is_empty() || trimmed.eq_ignore_ascii_case(OPENAI_COMPAT_UNSPECIFIED_MODEL) {
        None
    } else {
        Some(trimmed.to_string())
    }
}

pub(super) fn preferred_openai_agent_system_id(systems: &[AgentSystemSpec]) -> Option<String> {
    for preferred_id in ["openai_gpt5_native", "openai_gpt5_stdio"] {
        if systems.iter().any(|system| system.id == preferred_id) {
            return Some(preferred_id.to_string());
        }
    }
    systems
        .iter()
        .find(|system| matches!(system.transport, AgentSystemTransport::NativeOpenai))
        .map(|system| system.id.clone())
}

pub(super) fn preferred_anthropic_agent_system_id(systems: &[AgentSystemSpec]) -> Option<String> {
    for preferred_id in ["anthropic_claude_sonnet_native", "claude_sonnet_native"] {
        if systems.iter().any(|system| system.id == preferred_id) {
            return Some(preferred_id.to_string());
        }
    }
    systems
        .iter()
        .find(|system| matches!(system.transport, AgentSystemTransport::NativeAnthropic))
        .map(|system| system.id.clone())
}

pub(super) fn preferred_mistral_agent_system_id(systems: &[AgentSystemSpec]) -> Option<String> {
    for preferred_id in ["mistral_large_native", "mistral_native"] {
        if systems.iter().any(|system| system.id == preferred_id) {
            return Some(preferred_id.to_string());
        }
    }
    systems
        .iter()
        .find(|system| matches!(system.transport, AgentSystemTransport::NativeMistral))
        .map(|system| system.id.clone())
}

pub(super) fn preferred_local_agent_system_id(systems: &[AgentSystemSpec]) -> Option<String> {
    for preferred_id in [
        "msty_mlx_local_compat_template",
        "msty_local_compat_template",
        "local_llama_compat",
        "jan_local_compat_template",
    ] {
        if systems.iter().any(|system| system.id == preferred_id) {
            return Some(preferred_id.to_string());
        }
    }
    systems
        .iter()
        .find(|system| matches!(system.transport, AgentSystemTransport::NativeOpenaiCompat))
        .map(|system| system.id.clone())
}

pub(super) fn agent_prompt_template_options() -> &'static [(&'static str, &'static str)] {
    &[
        ("structured", "Structured (recommended)"),
        ("compact_intro", "Compact intro (no state)"),
        ("candidate_anchors", "Candidate between anchors"),
        ("blast_specificity", "BLAST specificity check"),
        ("track_intersection", "Track import + prioritization"),
        ("macro_template", "Macro/template authoring"),
    ]
}

pub(super) fn agent_prompt_template_includes_state_summary_by_default(id: &str) -> bool {
    !matches!(id, "compact_intro")
}

pub(super) fn agent_prompt_template_label(id: &str) -> &'static str {
    agent_prompt_template_options()
        .iter()
        .find(|(value, _)| *value == id)
        .map(|(_, label)| *label)
        .unwrap_or("Structured (recommended)")
}

pub(super) fn agent_prompt_template_text(id: &str) -> &'static str {
    match id {
        "compact_intro" => {
            r#"Task:
Introduce yourself briefly as GENtle's internal Agent Assistant.

Context policy:
- Do not assume any loaded project or sequence context.
- Answer in the user's language if it is clear from their wording; otherwise use concise English.
- Keep the reply compact enough for a live demo.
- Use GENtle docs/glossary.json plus docs/cli.md operand conventions when proposing commands.
- If an operand such as QUERY, ID, SEQ_ID, ENTRY_ID, PATH, or SPECIES is unclear, ask instead of guessing.

Output wanted:
- 5-8 bullets about what you can help with inside GENtle.
- 2-4 safe suggested_commands using GENtle shared-shell commands only.
- Each suggested command needs a clear title as the user intent and preconditions[] when it depends on state.
- Each suggested command should include expected_outcomes[] describing what should be observable if the command succeeds; these are expected effects, not guarantees.
- When you know the structured fact-graph logic, include precondition_expr and expected_effects alongside the prose; otherwise omit them rather than guessing.
- For negative requirements such as "no EcoRI site", require or produce a positive proof fact such as restriction_site.absent based on a complete scan; do not infer absence from missing features.
- On an empty or unknown project, prefer orientation/open/retrieve commands first; do not suggest feature scans as runnable first actions.
- Mark runnable suggestions execution="ask"; use execution="chat" only for purely explanatory rows that should not run.
- Mention that external database/network actions require explicit confirmation.

Good first-step demo commands:
- /help
- /list (current GENtle project state and loaded sequences, not filesystem files)
- state-summary
- /open (GUI file picker for local sequence files)
- /open file test_files/pGEX_3X.fa --id pgex
- /paste sequence --sequence-text GAATTCGCGGCCGCTTCTAGA --id demo_seq
- /fetch ensembl FUS --species homo_sapiens --id fus_live
- /fetch genbank NM_001126241.3 --id tp73_cdna

Follow-up demo command after a sequence exists:
- /features restriction-scan demo_seq --enzyme EcoRI
  preconditions[] = ["Sequence demo_seq exists in the current GENtle project."]
  expected_outcomes[] = ["A restriction-site report for demo_seq is available if the scan succeeds."]
  expected_effects[] may include restriction_site.absent only when the scan proves zero matching sites over the stated range.

Continuing earlier work:
- If the user wants to continue an earlier project, suggest GUI path `File -> Open Project...` or `File -> Open Recent Project...`.
- Do not invent a shell/slash command for recent projects until GENtle exposes one.
- If the user supplies an exact saved project path, tell them to open it through `File -> Open Project...` or by launching GENtle with that project path.

Do not describe /list as a directory listing. Do not suggest placeholder commands such as /open file PATH [--id ID] unless the user supplied a real PATH."#
        }
        "candidate_anchors" => {
            r#"Objective:
Generate candidate windows between two local anchors and rank them.

Context:
Project sequence ID: <SEQ_ID>

Inputs:
- seq_id: <SEQ_ID>
- anchor A: <feature boundary or absolute position>
- anchor B: <feature boundary or absolute position>

Constraints:
- length: 20
- step: 1
- GC range: 40-80%
- additional constraints: <motifs/sites/strand>

Output wanted:
- exact `gentle_cli shell "candidates ..."` commands
- scoring + filter + top-k steps
- validation checklist

Execution policy:
ask-before-run"#
        }
        "blast_specificity" => {
            r#"Objective:
Run a specificity check for one sequence with BLAST.

Context:
Target catalog: genomes | helpers

Inputs:
- genome_id/helper_id: <ID>
- query_sequence: <ACGT...>

Constraints:
- max_hits: 20
- task: blastn-short

Output wanted:
- exact BLAST command
- concise interpretation checklist for top hits

Execution policy:
chat-only"#
        }
        "track_intersection" => {
            r#"Objective:
Import external track evidence and prioritize candidates near track features.

Context:
Anchored sequence ID: <SEQ_ID>

Inputs:
- seq_id: <SEQ_ID>
- track file path: <BED/BED.GZ/BIGWIG/VCF path>

Constraints:
- keep imported features in TRACK groups
- do not modify original sequence content

Output wanted:
- exact track-import commands
- candidate generation/scoring/filter commands near imported TRACK features
- validation checklist

Execution policy:
ask-before-run"#
        }
        "macro_template" => {
            r#"Objective:
Create or update a reusable candidate macro template and run it with bindings.

Context:
Template name: <NAME>

Inputs:
- template parameters: <param list>
- script intent: <generate/score/filter/top-k/...>

Constraints:
- transactional run enabled
- deterministic tie-break policy where applicable

Output wanted:
- `candidates template-put` and `candidates template-run` commands
- brief note on expected outputs and rollback behavior

Execution policy:
ask-before-run"#
        }
        _ => {
            r#"Objective:
<one clear goal>

Documentation context:
Use GENtle docs/glossary.json for command paths, docs/cli.md for operand conventions, and docs/protocol.md for request/response semantics. If a placeholder such as QUERY, ID, SEQ_ID, ENTRY_ID, PATH, or SPECIES is unclear, ask instead of guessing.

Context:
<sequence/genome/helper IDs and short background>

Inputs:
- seq_id / genome_id / helper_id: ...
- anchors or coordinates: ...
- feature labels/kinds: ...

Constraints:
- length: ...
- GC range: ...
- motifs/sites to require or avoid: ...
- strand assumptions: ...

Output wanted:
- plan
- exact gentle_cli commands
- validation checklist

Execution policy:
chat-only | ask-before-run | allow-auto-exec"#
        }
    }
}

pub(super) fn default_agent_timeout_secs_string() -> String {
    default_env_string(AGENT_TIMEOUT_SECS_ENV)
}

pub(super) fn default_agent_connect_timeout_secs_string() -> String {
    default_env_string(AGENT_CONNECT_TIMEOUT_SECS_ENV)
}

pub(super) fn default_agent_read_timeout_secs_string() -> String {
    default_env_string(AGENT_READ_TIMEOUT_SECS_ENV)
}

pub(super) fn default_agent_max_retries_string() -> String {
    default_env_string(AGENT_MAX_RETRIES_ENV)
}

pub(super) fn default_agent_max_response_bytes_string() -> String {
    default_env_string(AGENT_MAX_RESPONSE_BYTES_ENV)
}

fn default_env_string(name: &str) -> String {
    env::var(name)
        .ok()
        .map(|v| v.trim().to_string())
        .filter(|v| !v.is_empty())
        .unwrap_or_default()
}

#[cfg(test)]
mod tests {
    use super::*;

    fn test_agent_system(id: &str, transport: AgentSystemTransport) -> AgentSystemSpec {
        AgentSystemSpec {
            id: id.to_string(),
            label: id.to_string(),
            description: None,
            transport,
            command: vec![],
            env: Default::default(),
            base_url: None,
            model: None,
            working_dir: None,
        }
    }

    #[test]
    fn model_name_normalization_treats_blank_and_unspecified_as_absent() {
        assert_eq!(normalize_agent_model_name(""), None);
        assert_eq!(normalize_agent_model_name("   "), None);
        assert_eq!(
            normalize_agent_model_name(OPENAI_COMPAT_UNSPECIFIED_MODEL),
            None
        );
        assert_eq!(
            normalize_agent_model_name(" gpt-test "),
            Some("gpt-test".to_string())
        );
    }

    #[test]
    fn prompt_templates_keep_default_and_fallback_stable() {
        assert!(
            agent_prompt_template_options()
                .iter()
                .any(|(id, _)| *id == AGENT_PROMPT_TEMPLATE_DEFAULT_ID)
        );
        assert_eq!(
            agent_prompt_template_label("unknown"),
            "Structured (recommended)"
        );
        assert!(
            !agent_prompt_template_includes_state_summary_by_default("compact_intro"),
            "compact intro should stay stateless for fast live demos"
        );
        assert!(agent_prompt_template_includes_state_summary_by_default(
            "structured"
        ));
        assert!(agent_prompt_template_text("compact_intro").contains("shared-shell"));
        assert!(agent_prompt_template_text("compact_intro").contains("not filesystem files"));
        assert!(agent_prompt_template_text("compact_intro").contains("Do not describe /list"));
        assert!(
            agent_prompt_template_text("compact_intro")
                .contains("/open file test_files/pGEX_3X.fa --id pgex")
        );
        assert!(
            agent_prompt_template_text("compact_intro")
                .contains("/fetch ensembl FUS --species homo_sapiens --id fus_live")
        );
        assert!(agent_prompt_template_text("compact_intro").contains("preconditions[]"));
        assert!(agent_prompt_template_text("compact_intro").contains("expected_outcomes[]"));
        assert!(agent_prompt_template_text("compact_intro").contains("precondition_expr"));
        assert!(agent_prompt_template_text("compact_intro").contains("expected_effects"));
        assert!(agent_prompt_template_text("compact_intro").contains("restriction_site.absent"));
        assert!(
            agent_prompt_template_text("compact_intro")
                .contains("Sequence demo_seq exists in the current GENtle project.")
        );
        assert!(
            agent_prompt_template_text("compact_intro")
                .contains("restriction-site report for demo_seq is available")
        );
        assert!(
            agent_prompt_template_text("compact_intro")
                .contains("do not suggest feature scans as runnable first actions")
        );
        assert!(
            agent_prompt_template_text("compact_intro").contains("File -> Open Recent Project...")
        );
        assert!(
            agent_prompt_template_text("compact_intro")
                .contains("Do not invent a shell/slash command for recent projects")
        );
        assert!(agent_prompt_template_text("compact_intro").contains("docs/glossary.json"));
        assert!(agent_prompt_template_text("candidate_anchors").contains("candidates"));
        assert!(agent_prompt_template_text("unknown").contains("Objective:"));
        assert!(agent_prompt_template_text("unknown").contains("operand conventions"));
    }

    #[test]
    fn preferred_agent_system_helpers_select_specific_templates_first() {
        let systems = vec![
            test_agent_system(
                "msty_mlx_local_compat_template",
                AgentSystemTransport::NativeOpenaiCompat,
            ),
            test_agent_system(
                "local_llama_compat",
                AgentSystemTransport::NativeOpenaiCompat,
            ),
            test_agent_system("openai_fallback", AgentSystemTransport::NativeOpenai),
            test_agent_system("openai_gpt5_native", AgentSystemTransport::NativeOpenai),
        ];
        assert_eq!(
            preferred_openai_agent_system_id(&systems).as_deref(),
            Some("openai_gpt5_native")
        );
        assert_eq!(
            preferred_local_agent_system_id(&systems).as_deref(),
            Some("msty_mlx_local_compat_template")
        );
    }
}
