//! Agent Assistant defaults, prompt templates, and system-selection helpers.
//!
//! Keeping this small policy bundle outside `app.rs` makes the top-level GUI
//! coordinator less dense while preserving the same Agent Assistant behavior.

use std::{collections::HashMap, env, fs, io::ErrorKind, path::Path};

use crate::agent_bridge::{
    AGENT_CONNECT_TIMEOUT_SECS_ENV, AGENT_MAX_RESPONSE_BYTES_ENV, AGENT_MAX_RETRIES_ENV,
    AGENT_READ_TIMEOUT_SECS_ENV, AGENT_TIMEOUT_SECS_ENV, ANTHROPIC_API_KEY_ENV, AgentSystemSpec,
    AgentSystemTransport, MISTRAL_API_KEY_ENV, OPENAI_API_KEY_ENV, OPENAI_COMPAT_UNSPECIFIED_MODEL,
};

pub(super) const AGENT_PROMPT_TEMPLATE_DEFAULT_ID: &str = "structured";
const AGENT_TOKEN_FILE_MAX_BYTES: u64 = 16 * 1024;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(super) enum AgentTokenFileStatus {
    Loaded,
    Missing,
    Empty,
    Invalid,
    Unreadable,
}

#[derive(Clone)]
pub(super) struct AgentTokenFileCredential {
    pub(super) display_path: &'static str,
    pub(super) status: AgentTokenFileStatus,
    pub(super) detail: String,
    pub(super) broadly_readable: bool,
    value: Option<String>,
}

impl AgentTokenFileCredential {
    pub(super) fn value(&self) -> Option<&str> {
        self.value.as_deref()
    }
}

struct AgentTokenFileSpec {
    env_var: &'static str,
    file_name: &'static str,
    display_path: &'static str,
}

const AGENT_TOKEN_FILE_SPECS: [AgentTokenFileSpec; 3] = [
    AgentTokenFileSpec {
        env_var: OPENAI_API_KEY_ENV,
        file_name: ".codex_token",
        display_path: "~/.codex_token",
    },
    AgentTokenFileSpec {
        env_var: ANTHROPIC_API_KEY_ENV,
        file_name: ".claude_token",
        display_path: "~/.claude_token",
    },
    AgentTokenFileSpec {
        env_var: MISTRAL_API_KEY_ENV,
        file_name: ".mistral_token",
        display_path: "~/.mistral_token",
    },
];

pub(super) fn agent_api_key_source(
    transport: AgentSystemTransport,
) -> Option<(&'static str, &'static str)> {
    match transport {
        AgentSystemTransport::NativeOpenai => Some((OPENAI_API_KEY_ENV, "~/.codex_token")),
        AgentSystemTransport::NativeAnthropic => Some((ANTHROPIC_API_KEY_ENV, "~/.claude_token")),
        AgentSystemTransport::NativeMistral => Some((MISTRAL_API_KEY_ENV, "~/.mistral_token")),
        _ => None,
    }
}

pub(super) fn load_agent_token_file_credentials(
    home: Option<&Path>,
) -> HashMap<String, AgentTokenFileCredential> {
    let mut credentials = HashMap::new();
    for spec in AGENT_TOKEN_FILE_SPECS {
        let Some(home) = home else {
            credentials.insert(
                spec.env_var.to_string(),
                AgentTokenFileCredential {
                    display_path: spec.display_path,
                    status: AgentTokenFileStatus::Unreadable,
                    detail: "home directory is unavailable".to_string(),
                    broadly_readable: false,
                    value: None,
                },
            );
            continue;
        };
        let path = home.join(spec.file_name);
        let metadata = match fs::metadata(&path) {
            Ok(metadata) => metadata,
            Err(err) if err.kind() == ErrorKind::NotFound => {
                credentials.insert(
                    spec.env_var.to_string(),
                    AgentTokenFileCredential {
                        display_path: spec.display_path,
                        status: AgentTokenFileStatus::Missing,
                        detail: String::new(),
                        broadly_readable: false,
                        value: None,
                    },
                );
                continue;
            }
            Err(err) => {
                credentials.insert(
                    spec.env_var.to_string(),
                    AgentTokenFileCredential {
                        display_path: spec.display_path,
                        status: AgentTokenFileStatus::Unreadable,
                        detail: err.to_string(),
                        broadly_readable: false,
                        value: None,
                    },
                );
                continue;
            }
        };
        let broadly_readable = token_file_is_broadly_readable(&metadata);
        if metadata.len() > AGENT_TOKEN_FILE_MAX_BYTES {
            credentials.insert(
                spec.env_var.to_string(),
                AgentTokenFileCredential {
                    display_path: spec.display_path,
                    status: AgentTokenFileStatus::Invalid,
                    detail: format!(
                        "file is larger than the {}-byte safety limit",
                        AGENT_TOKEN_FILE_MAX_BYTES
                    ),
                    broadly_readable,
                    value: None,
                },
            );
            continue;
        }
        let raw = match fs::read_to_string(&path) {
            Ok(raw) => raw,
            Err(err) => {
                credentials.insert(
                    spec.env_var.to_string(),
                    AgentTokenFileCredential {
                        display_path: spec.display_path,
                        status: AgentTokenFileStatus::Unreadable,
                        detail: err.to_string(),
                        broadly_readable,
                        value: None,
                    },
                );
                continue;
            }
        };
        let trimmed = raw.trim();
        let (status, detail, value) = if trimmed.is_empty() {
            (
                AgentTokenFileStatus::Empty,
                "file is empty".to_string(),
                None,
            )
        } else if trimmed.chars().any(char::is_whitespace) {
            (
                AgentTokenFileStatus::Invalid,
                "expected exactly one token without whitespace".to_string(),
                None,
            )
        } else {
            (
                AgentTokenFileStatus::Loaded,
                String::new(),
                Some(trimmed.to_string()),
            )
        };
        credentials.insert(
            spec.env_var.to_string(),
            AgentTokenFileCredential {
                display_path: spec.display_path,
                status,
                detail,
                broadly_readable,
                value,
            },
        );
    }
    credentials
}

#[cfg(unix)]
fn token_file_is_broadly_readable(metadata: &fs::Metadata) -> bool {
    use std::os::unix::fs::PermissionsExt;

    metadata.permissions().mode() & 0o077 != 0
}

#[cfg(not(unix))]
fn token_file_is_broadly_readable(_metadata: &fs::Metadata) -> bool {
    false
}

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
        (
            "regulatory_reporter_study",
            "Evidence-guided regulatory reporter study",
        ),
        ("biological_adjustment", "Adapt a biological workflow"),
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
- Include Configuration among useful starting actions when x_gui_context is available.
- Use x_gui_context as the authoritative list of recent projects, tutorial projects, and Configuration sections; use exact open_command values from its rows.
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
- If x_gui_context contains recent_projects, list those rows and use the selected row's exact `ui open recent-project ITEM_ID` command.
- If x_gui_context contains tutorial_projects, list those rows instead of claiming that no tutorial projects are known.
- Use `ui open configuration SECTION` from configuration_sections to put the user on the relevant settings tab; opening it does not itself change credentials or paths.
- If no recent row matches, suggest GUI path `File -> Open Project...` or `File -> Open Recent Project...` and report any x_gui_context warning.
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
        "regulatory_reporter_study" => {
            r#"Objective:
Plan an evidence-guided study that places native regulatory DNA from responsive genes into a validated reporter backbone.

Context to establish before suggesting execution:
- perturbation factor and form/isoform: <FACTOR_AND_FORM>
- biological system, cell type, time point, and contrast: <CONTEXT>
- responsive target genes and evidence source: <GENES_AND_EVIDENCE>
- species and exact genome assembly: <SPECIES_AND_BUILD>
- exact reporter-vector catalog identity: <VECTOR_OR_UNKNOWN>

Required distinctions:
- expression response, sequence motif, occupancy evidence, and functional reporter evidence are separate evidence classes;
- a motif match alone supports a candidate binding site, not occupancy, causality, or a cofactor claim;
- enumerate transcript-derived TSS/promoter classes before choosing one fragment per gene;
- keep native regulatory fragments as the default; engineered motif controls are optional per member and require an explicit mutation policy;
- when retaining the complete annotated 5' UTR but excluding source translation, use canonical_cds_start_exclusive with one explicit transcript_id;
- treat ready_to_plan, ready_for_approval, and ready_to_materialize as distinct transitions; only the last permits the exact approved proposal to proceed;
- do not materialize constructs, download external data, or run a mutation without explicit confirmation.

Useful shared GENtle surfaces after their operands are known:
- `features promoter-evidence-matrix SEQ_ID --gene-label GENE ...`
- `op '{"SuggestPromoterReporterFragments":{...}}'`
- `promoters panel-readiness request @REQUEST.json --path READINESS.json`
- `promoters panel-plan @REQUEST.json --path PROPOSAL.json`
- `promoters panel-readiness proposal @PROPOSAL.json [--approve DIGEST] --path READINESS.json`
- `promoters panel-materialize @PROPOSAL.json --approve DIGEST`

Output wanted:
1. assumptions and missing inputs;
2. one evidence ledger per gene, keeping response/motif/occupancy provenance separate;
3. TSS/promoter alternatives and fragment-policy rationale;
4. a native-only panel request first, with optional member-specific controls listed separately;
5. exact parser-valid commands only after every operand is known;
6. proposal review checklist covering vector identity, fragment bounds, anchor/evidence provenance, primers, product annotations, and non-claims;
7. an explicit stop point before materialization.

Execution policy:
ask-before-run; chat-only while species/build, evidence, TSS policy, vector identity, or local paths remain unresolved"#
        }
        "biological_adjustment" => {
            r#"Objective:
Determine how GENtle should support this biological adjustment or related workflow request:
<DESCRIBE THE REQUEST>

Context:
- existing GENtle workflow or operation, if known: <WORKFLOW_OR_UNKNOWN>
- project sequence/report IDs, if relevant: <IDS_OR_UNKNOWN>
- biological purpose: <WHY_THIS_ADJUSTMENT_IS_NEEDED>

Fixed requirements:
- sequence orientation and coordinate interpretation: <DETAILS>
- values or motifs that must remain exact: <DETAILS>
- evidence inputs that are authoritative: <DETAILS>

Adjustable policy:
- candidate window or search scope: <DETAILS>
- scoring/ranking preferences: <DETAILS>
- acceptable ambiguity and required user choices: <DETAILS>

Classify the request first:
1. parameter-only use of an existing capability;
2. composition of existing GENtle operations or reports;
3. a missing engine capability that needs implementation.

Output wanted:
- cite capability/readiness evidence before proposing execution;
- if already supported, provide exact parser-valid GENtle commands and validation;
- if unsupported, do not invent a command. Produce a concise implementation brief covering the biological invariant, reusable GENtle engine services, versioned request/report, provenance, explicit non-claims, adapter reachability, and deterministic edge-case tests;
- distinguish what GENtle can execute now from what a coding agent or contributor must add.

Execution policy:
chat-only until the classification and biological contract have been reviewed"#
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
    use tempfile::tempdir;

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
            supports_image_attachments: false,
            supports_web_research: false,
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
            agent_prompt_template_text("compact_intro").contains("ui open recent-project ITEM_ID")
        );
        assert!(agent_prompt_template_text("compact_intro").contains("tutorial_projects"));
        assert!(
            agent_prompt_template_text("compact_intro").contains("ui open configuration SECTION")
        );
        assert!(agent_prompt_template_text("compact_intro").contains("docs/glossary.json"));
        assert!(agent_prompt_template_text("candidate_anchors").contains("candidates"));
        assert!(
            agent_prompt_template_text("regulatory_reporter_study")
                .contains("canonical_cds_start_exclusive")
        );
        assert!(
            agent_prompt_template_text("regulatory_reporter_study").contains("motif match alone")
        );
        assert!(
            agent_prompt_template_text("regulatory_reporter_study")
                .contains("`op '{\"SuggestPromoterReporterFragments\"")
        );
        assert!(
            agent_prompt_template_text("regulatory_reporter_study")
                .contains("promoters panel-readiness request")
        );
        assert!(!agent_prompt_template_text("regulatory_reporter_study").contains("op apply"));
        assert!(
            agent_prompt_template_text("biological_adjustment")
                .contains("parameter-only use of an existing capability")
        );
        assert!(
            agent_prompt_template_text("biological_adjustment")
                .contains("if unsupported, do not invent a command")
        );
        assert!(
            agent_prompt_template_text("biological_adjustment")
                .contains("versioned request/report")
        );
        assert!(agent_prompt_template_text("unknown").contains("Objective:"));
        assert!(agent_prompt_template_text("unknown").contains("operand conventions"));
    }

    #[test]
    fn regulatory_reporter_template_uses_a_parser_valid_operation_command() {
        let line = r#"op '{"SuggestPromoterReporterFragments":{"input":"seq","gene_label":"GENE","retain_downstream_from_tss_bp":200,"retain_upstream_beyond_variant_bp":500,"max_candidates":10,"fragment_policy":{"anchor":{"kind":"motif_hit","label_or_id":"FACTOR","occurrence":1},"collapse_mode":{"tss_cluster":{"tolerance_bp":50}},"promoter_upstream_baseline_bp":500,"anchor_flank_bp":150,"max_fragment_length_bp":5000},"path":"candidates.json"}}'"#;

        assert!(matches!(
            crate::engine_shell::parse_shell_line(line).expect("template operation command"),
            crate::engine_shell::ShellCommand::Op { .. }
        ));
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

    #[test]
    fn native_provider_token_files_have_stable_fallback_paths() {
        assert_eq!(
            agent_api_key_source(AgentSystemTransport::NativeOpenai),
            Some((OPENAI_API_KEY_ENV, "~/.codex_token"))
        );
        assert_eq!(
            agent_api_key_source(AgentSystemTransport::NativeAnthropic),
            Some((ANTHROPIC_API_KEY_ENV, "~/.claude_token"))
        );
        assert_eq!(
            agent_api_key_source(AgentSystemTransport::NativeMistral),
            Some((MISTRAL_API_KEY_ENV, "~/.mistral_token"))
        );
        assert_eq!(
            agent_api_key_source(AgentSystemTransport::ExternalJsonStdio),
            None,
            "Codex/Claude local stdio transports must retain their own login semantics"
        );
    }

    #[test]
    fn token_file_loader_reads_single_tokens_and_reports_missing_files() {
        let temp = tempdir().expect("temp home");
        fs::write(temp.path().join(".codex_token"), "codex-test-token\n")
            .expect("write codex token");
        fs::write(temp.path().join(".claude_token"), "two words")
            .expect("write invalid claude token");

        let loaded = load_agent_token_file_credentials(Some(temp.path()));
        let codex = loaded.get(OPENAI_API_KEY_ENV).expect("codex credential");
        assert_eq!(codex.status, AgentTokenFileStatus::Loaded);
        assert_eq!(codex.value(), Some("codex-test-token"));

        let claude = loaded
            .get(ANTHROPIC_API_KEY_ENV)
            .expect("claude credential");
        assert_eq!(claude.status, AgentTokenFileStatus::Invalid);
        assert!(claude.value().is_none());

        let mistral = loaded.get(MISTRAL_API_KEY_ENV).expect("mistral credential");
        assert_eq!(mistral.status, AgentTokenFileStatus::Missing);
        assert!(mistral.value().is_none());
    }

    #[cfg(unix)]
    #[test]
    fn token_file_loader_reports_permissions_visible_beyond_owner() {
        use std::os::unix::fs::PermissionsExt;

        let temp = tempdir().expect("temp home");
        let path = temp.path().join(".codex_token");
        fs::write(&path, "codex-test-token").expect("write codex token");
        fs::set_permissions(&path, fs::Permissions::from_mode(0o644))
            .expect("set token permissions");

        let loaded = load_agent_token_file_credentials(Some(temp.path()));
        assert!(
            loaded
                .get(OPENAI_API_KEY_ENV)
                .expect("codex credential")
                .broadly_readable
        );
    }
}
