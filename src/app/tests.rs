use super::{
    AGENT_BASE_URL_ENV, AGENT_CONNECT_TIMEOUT_SECS_ENV, AGENT_MAX_RESPONSE_BYTES_ENV,
    AGENT_MAX_RETRIES_ENV, AGENT_MODEL_ENV, AGENT_READ_TIMEOUT_SECS_ENV, AGENT_TIMEOUT_SECS_ENV,
    ANTHROPIC_API_KEY_AUTH_HINT, ANTHROPIC_API_KEY_ENV, APP_CONFIGURATION_SCHEMA_VERSION,
    BACKGROUND_JOB_HISTORY_METADATA_KEY, BACKGROUND_JOB_HISTORY_SCHEMA,
    BACKGROUND_JOBS_RECENT_JOB_EVENTS_SCROLL_ID, BACKGROUND_JOBS_RETRY_CLEANUP_AUDIT_SCROLL_ID,
    BACKGROUND_JOBS_RETRY_SNAPSHOTS_REMOVED_PREVIEW_SCROLL_ID,
    BACKGROUND_JOBS_RETRY_SNAPSHOTS_RETAINED_PREVIEW_SCROLL_ID,
    BACKGROUND_JOBS_RETRY_SNAPSHOTS_SCROLL_ID, BackgroundJobEventPhase, BackgroundJobKind,
    CacheCleanupScope, CloningRoutineCatalogRow, CommandPaletteAction, ConfigurationTab,
    ContainerRow, DEFAULT_DBSNP_TUTORIAL_RS_ID, DEFAULT_HELPER_GENOME_CACHE_DIR,
    DEFAULT_HELPER_GENOME_CATALOG_PATH, DEFAULT_LINEAGE_MAIN_SPLIT_FRACTION, DbSnpFetchTask,
    DbSnpFetchTaskMessage, EngineError, ErrorCode, GENtleApp, GenomeBlastOptionsPreset,
    GenomeBlastTask, GenomeBlastTaskMessage, GenomeDialogScope, GenomePrepareLaunchMode,
    GenomePrepareTask, GenomePrepareTaskMessage, GenomeTrackImportTask,
    GenomeTrackImportTaskMessage, GibsonUiInsertOrientation, GibsonUiInsertRow,
    GibsonUiOpeningMode, HelpDoc, HelpSearchMatch, HelpTutorialDocEntry,
    LINEAGE_GRAPH_WORKSPACE_METADATA_KEY, LINEAGE_MAIN_TOP_PANEL_MIN_HEIGHT, LineageAnalysisKind,
    LineageCopyPayloadKind, LineageNodeKind, LineageRow, MAX_RECENT_PROJECTS,
    MISTRAL_API_KEY_AUTH_HINT, MISTRAL_API_KEY_ENV, OPENAI_API_KEY_ENV,
    OPERATION_HISTORY_SCROLL_ID, PendingEnsemblCatalogUpdateDialog,
    PendingEnsemblInstallableGenomeDialog, PendingEnsemblQuickInstallDialog,
    PersistedConfiguration, PersistedLineageGraphWorkspace, PersistedLineageNodeGroup,
    PersistedRackWorkspace, PrepareGenomeDialogPrimaryAction, PrepareGenomeFailureRecovery,
    PrepareGenomeUiStepState, PrepareGenomeUiStepStatus, PreparedGenomeReinstallDialogHost,
    PreparedGenomeReinstallRequest, ProjectAction, ProjectOverviewMetric, ProjectOverviewTarget,
    RACK_HELP_AUTO_MINIMIZE_MOVE_THRESHOLD, RACK_WORKSPACE_METADATA_KEY,
    ROUTINE_DECISION_TRACE_SCHEMA, ROUTINE_DECISION_TRACE_STORE_SCHEMA,
    ROUTINE_DECISION_TRACES_METADATA_KEY, RackDragState, RetryCleanupAuditActionFilter,
    RetrySnapshotKindFilter, RetrySnapshotPendingCleanupAction, RoutineAssistantStage,
    TutorialProjectOpenOutcome, TutorialProjectTask, TutorialProjectTaskMessage,
    TutorialProjectTaskProgress, gui_prominent_glossary_entries,
    preferred_anthropic_agent_system_id, preferred_local_agent_system_id,
    preferred_mistral_agent_system_id, preferred_openai_agent_system_id,
};
use crate::{
    agent_bridge::{AgentSystemSpec, AgentSystemTransport},
    dna_sequence::DNAsequence,
    engine::{
        Arrangement, ArrangementMode, BlastHitFeatureInput, BlastInvocationProvenance, Container,
        ContainerKind, DbSnpFetchProgress, DbSnpFetchStage, DisplaySettings, DotplotMode, Engine,
        FlexibilityModel, GenomeAnnotationProjectionTelemetry, GenomeGeneExtractMode, GentleEngine,
        LineageEdge, LineageNode, LinearSequenceLetterLayoutMode, OpResult, Operation,
        PLANNING_ESTIMATE_SCHEMA, PairwiseAlignmentMode, PlanningEstimate, PlanningObjective,
        PrimerDesignPairConstraint, PrimerDesignSideConstraint, ProjectState,
        ProteinToDnaHandoffRankingGoal, ProteinToDnaHandoffStrategy, Rack, RackAuthoringTemplate,
        RackFillDirection, RackProfileKind, RackProfileSnapshot, RenderSvgMode,
        RestrictionCloningPcrHandoffMode, RestrictionEnzymeDisplayMode, ReverseTranslationReport,
        RoutineDecisionTraceDisambiguationAnswer, RoutineDecisionTraceDisambiguationQuestion,
        RoutineDecisionTracePreflightSnapshot, RoutineDecisionTraceStore, SequenceOrigin,
        TranslationSpeedMark, TranslationSpeedProfile, TranslationSpeedProfileSource,
        UniprotFeatureCodingDnaQueryMode, UniprotFeatureCodingDnaQueryReport,
    },
    engine_shell::{ShellCommand, UiIntentTarget, parse_shell_line},
    ensembl_protein::{EnsemblProteinEntry, EnsemblProteinEntrySummary, EnsemblProteinFeature},
    genomes::{
        EnsemblCatalogUpdatePreview, EnsemblInstallableGenomeCatalog, EnsemblQuickInstallPreview,
        HelperConstructInterpretation, PrepareGenomePlan, PrepareGenomePlanStep,
        PrepareGenomeProgress, PrepareGenomeStepId, PreparedCacheArtifactGroup,
        PreparedCacheArtifactStat, PreparedCacheCleanupItemReport, PreparedCacheCleanupMode,
        PreparedCacheCleanupReport, PreparedCacheEntryKind, PreparedCacheInspectionEntry,
        PreparedCacheInspectionReport, PreparedGenomeInspection,
    },
    gibson_planning::{
        GibsonAssemblyPlan, GibsonAssemblyPreview, GibsonCartoonPreview,
        GibsonDesignAdjustmentTarget, GibsonDestinationOpeningSuggestion, GibsonPreviewDestination,
        GibsonPreviewInsert, GibsonPrimerSegment, GibsonPrimerSuggestion,
        GibsonResolvedJunctionPreview, GibsonRoutineHandoffPreview,
        GibsonSuggestedDesignAdjustment,
    },
    i18n::UiLanguage,
    uniprot::UniprotEntrySummary,
    window::Window,
};
use eframe::egui;
use std::{
    collections::{HashMap, HashSet},
    env, fs,
    path::{Path, PathBuf},
    sync::{
        Arc, RwLock,
        atomic::{AtomicBool, Ordering},
        mpsc,
    },
    time::{Duration, Instant},
};
use tempfile::tempdir;

#[test]
fn splash_screen_is_visible_for_fresh_app() {
    let app = GENtleApp::default();
    assert!(app.splash_should_render_at(app.splash_started_at));
}

#[test]
fn splash_screen_hides_after_timeout_or_dismissal() {
    let mut app = GENtleApp::default();
    let now = Instant::now();
    app.splash_started_at = now - Duration::from_secs(10);
    assert!(!app.splash_should_render_at(now));

    app.splash_started_at = now;
    app.dismiss_splash_screen();
    assert!(!app.splash_should_render_at(now));
}

#[test]
fn direct_egui_windows_are_limited_to_window_shell_wrappers() {
    let manifest_dir = Path::new(env!("CARGO_MANIFEST_DIR"));
    let mut source_files = Vec::new();
    collect_rust_source_files(&manifest_dir.join("src"), &mut source_files);
    let mut unexpected = Vec::new();

    for path in source_files {
        let relative = path
            .strip_prefix(manifest_dir)
            .expect("source path should be inside manifest dir")
            .to_string_lossy()
            .replace('\\', "/");
        if matches!(
            relative.as_str(),
            "src/egui_compat.rs" | "src/bin/gentle_egui_window_repro/main.rs"
        ) {
            continue;
        }
        if relative.ends_with("/tests.rs") {
            continue;
        }
        let source = fs::read_to_string(&path)
            .unwrap_or_else(|err| panic!("read {relative} for window-shell guard: {err}"));
        collect_direct_window_new_lines(&relative, &source, &mut unexpected);
    }

    assert!(
        unexpected.is_empty(),
        "production windows must use egui_compat::show_hosted_window for hosted/specialist surfaces or egui_compat::show_modal_window for one-shot prompts; direct egui::Window::new belongs only in egui_compat or tests:\n{}",
        unexpected.join("\n")
    );
}

fn collect_rust_source_files(dir: &Path, files: &mut Vec<PathBuf>) {
    for entry in fs::read_dir(dir).unwrap_or_else(|err| panic!("read {}: {err}", dir.display())) {
        let entry = entry.expect("read source directory entry");
        let path = entry.path();
        if path.is_dir() {
            collect_rust_source_files(&path, files);
        } else if path.extension().and_then(|ext| ext.to_str()) == Some("rs") {
            files.push(path);
        }
    }
}

fn collect_direct_window_new_lines(relative: &str, source: &str, unexpected: &mut Vec<String>) {
    let mut pending_cfg_test = false;
    let mut skipping_test_module_depth: Option<i32> = None;
    for (line_index, line) in source.lines().enumerate() {
        if let Some(depth) = skipping_test_module_depth.as_mut() {
            *depth += brace_delta(line);
            if *depth <= 0 {
                skipping_test_module_depth = None;
            }
            continue;
        }

        let trimmed = line.trim();
        if pending_cfg_test && trimmed.starts_with("mod tests") {
            let depth = brace_delta(line);
            if depth > 0 {
                skipping_test_module_depth = Some(depth);
            }
            pending_cfg_test = false;
            continue;
        }
        if pending_cfg_test && !trimmed.starts_with("#[") {
            pending_cfg_test = false;
        }
        if trimmed == "#[cfg(test)]" {
            pending_cfg_test = true;
            continue;
        }

        if line.contains("egui::Window::new") {
            unexpected.push(format!("{relative}:{}: {}", line_index + 1, line.trim()));
        }
    }
}

fn brace_delta(line: &str) -> i32 {
    let opens = line.chars().filter(|ch| *ch == '{').count() as i32;
    let closes = line.chars().filter(|ch| *ch == '}').count() as i32;
    opens - closes
}

fn make_prepare_plan(step_ids: &[PrepareGenomeStepId]) -> PrepareGenomePlan {
    PrepareGenomePlan {
        genome_id: "ToyGenome".to_string(),
        steps: step_ids
            .iter()
            .copied()
            .map(|step_id| PrepareGenomePlanStep {
                step_id,
                label: step_id.label().to_string(),
                operation_summary: format!("Plan {}", step_id.label()),
                determinate_hint: step_id != PrepareGenomeStepId::BlastIndex,
            })
            .collect(),
    }
}

fn test_agent_system(id: &str, transport: AgentSystemTransport) -> AgentSystemSpec {
    AgentSystemSpec {
        id: id.to_string(),
        label: id.to_string(),
        description: None,
        transport,
        command: vec![],
        env: HashMap::new(),
        base_url: None,
        model: None,
        working_dir: None,
    }
}

fn write_toy_prepare_catalog(root: &std::path::Path) -> (String, String) {
    let fasta = root.join("toy.fa");
    let ann = root.join("toy.gtf");
    let cache_dir = root.join("cache");
    let catalog_path = root.join("catalog.json");
    fs::write(&fasta, ">chr1\nACGTACGT\n").expect("write fasta");
    fs::write(
        &ann,
        "chr1\tsrc\tgene\t1\t8\t.\t+\t.\tgene_id \"GENE1\"; gene_name \"ONE\";\n",
    )
    .expect("write annotation");
    fs::write(
        &catalog_path,
        format!(
            r#"{{
  "ToyGenome": {{
    "description": "toy test genome",
    "sequence_local": "{}",
    "annotations_local": "{}",
    "cache_dir": "{}"
  }}
}}"#,
            fasta.display(),
            ann.display(),
            cache_dir.display()
        ),
    )
    .expect("write catalog");
    (
        catalog_path.display().to_string(),
        cache_dir.display().to_string(),
    )
}

fn write_app_test_helper_catalog(root: &std::path::Path) -> String {
    let pgex_fasta = root.join("pgex.fa");
    let pgex_ann = root.join("pgex.gb");
    let neutral_fasta = root.join("neutral.fa");
    let neutral_ann = root.join("neutral.gb");
    let helper_cache_dir = root.join("helper_cache");
    let catalog_path = root.join("helper_catalog.json");
    fs::write(&pgex_fasta, ">pgex\nATGGGATCCGAACTCGAGTAA\n").expect("write pgex fasta");
    fs::write(
        &pgex_ann,
        "LOCUS       pGEXTEST                  22 bp    DNA\n",
    )
    .expect("write pgex annotation");
    fs::write(&neutral_fasta, ">neutral\nATGCGTACTGATCGATCGTA\n").expect("write neutral fasta");
    fs::write(
        &neutral_ann,
        "LOCUS       NEUTRAL                  20 bp    DNA\n",
    )
    .expect("write neutral annotation");
    fs::write(
        &catalog_path,
        format!(
            r#"{{
  "pGEX_like_vector": {{
    "description": "GST affinity fusion vector",
    "summary": "E. coli GST-tag fusion helper with Factor Xa cleavage before the MCS.",
    "aliases": ["pGEX", "PGEX-001"],
    "tags": ["expression", "purification"],
    "search_terms": ["factor xa", "gst", "affinity purification"],
    "helper_kind": "plasmid_vector",
    "host_system": "Escherichia coli",
    "sequence_availability": "source-backed synthetic test fixture",
    "redistribution_status": "synthetic fixture; redistribution unrestricted in tests",
    "biological_safety_note": "synthetic non-biological test data",
    "usable_as_empty_backbone": true,
    "procurement": {{
      "vendor_name": "ExampleBio",
      "catalog_number": "PGEX-001",
      "order_url": "https://example.invalid/pgex"
    }},
    "semantics": {{
      "schema": "gentle.helper_semantics.v1",
      "affordances": [
        "bacterial_expression",
        "affinity_purification",
        "fusion_tagging",
        "protease_tag_removal"
      ],
      "constraints": ["reading_frame_must_be_preserved"],
      "components": [
        {{
          "id": "gst_tag",
          "kind": "fusion_tag",
          "label": "GST",
          "description": "Glutathione affinity tag",
          "tags": ["purification"],
          "attributes": {{"position": "n_terminal"}}
        }},
        {{
          "id": "factor_xa_site",
          "kind": "protease_site",
          "label": "Factor Xa",
          "description": "Protease cleavage site"
        }},
        {{
          "id": "mcs",
          "kind": "cloning_site",
          "label": "MCS",
          "description": "In-frame coding insert window"
        }}
      ],
      "relationships": [
        {{
          "subject": "gst_tag",
          "predicate": "upstream_of",
          "object": "factor_xa_site",
          "note": "Tag removal before cloning window"
        }},
        {{
          "subject": "factor_xa_site",
          "predicate": "upstream_of",
          "object": "mcs"
        }}
      ]
    }},
    "sequence_local": "{}",
    "annotations_local": "{}",
    "cache_dir": "{}"
  }},
  "NeutralBackbone": {{
    "description": "Generic cloning backbone",
    "summary": "Neutral helper without typed semantics.",
    "aliases": ["NEUTRAL-1"],
    "helper_kind": "plasmid_vector",
    "host_system": "Escherichia coli",
    "sequence_availability": "source-backed synthetic test fixture",
    "redistribution_status": "synthetic fixture; redistribution unrestricted in tests",
    "biological_safety_note": "synthetic non-biological test data",
    "usable_as_empty_backbone": true,
    "semantics": {{
      "schema": "gentle.helper_semantics.v1",
      "components": [
        {{
          "id": "neutral_origin",
          "kind": "origin",
          "label": "Neutral origin"
        }}
      ]
    }},
    "sequence_local": "{}",
    "annotations_local": "{}",
    "cache_dir": "{}"
  }}
}}"#,
            pgex_fasta.display(),
            pgex_ann.display(),
            helper_cache_dir.display(),
            neutral_fasta.display(),
            neutral_ann.display(),
            helper_cache_dir.display()
        ),
    )
    .expect("write helper catalog");
    catalog_path.display().to_string()
}

fn make_prepared_genome_inspection() -> PreparedGenomeInspection {
    PreparedGenomeInspection {
        genome_id: "ToyGenome".to_string(),
        install_dir: "/tmp/toy_install".to_string(),
        manifest_path: "/tmp/toy_install/manifest.json".to_string(),
        sequence_source_type: "remote_http".to_string(),
        annotation_source_type: "remote_http".to_string(),
        sequence_source: "https://example.invalid/sequence.fa.gz".to_string(),
        annotation_source: "https://example.invalid/annotation.gtf.gz".to_string(),
        sequence_path: "/tmp/toy_install/sequence.fa".to_string(),
        annotation_path: "/tmp/toy_install/annotation.gtf".to_string(),
        fasta_index_path: "/tmp/toy_install/sequence.fa.fai".to_string(),
        gene_index_path: "/tmp/toy_install/genes.json".to_string(),
        transcript_index_path: Some("/tmp/toy_install/transcripts.json".to_string()),
        blast_db_prefix: Some("/tmp/toy_install/blastdb/genome".to_string()),
        blast_index_ready: true,
        sequence_sha1: None,
        annotation_sha1: None,
        sequence_present: true,
        annotation_present: true,
        fasta_index_ready: true,
        gene_index_ready: true,
        transcript_index_ready: true,
        total_size_bytes: 1234,
        installed_at_unix_ms: 0,
        cached_contig_count: 7,
        cached_total_span_bp: 1000,
        cached_longest_contig: Some("chr1".to_string()),
        cached_longest_contig_bp: Some(1000),
        cached_contig_preview: vec!["chr1".to_string()],
    }
}

fn make_test_app_with_open_windows(seq_ids: &[&str]) -> GENtleApp {
    let mut state = ProjectState::default();
    for seq_id in seq_ids {
        state.sequences.insert(
            (*seq_id).to_string(),
            DNAsequence::from_sequence("ACGTACGT").expect("sequence"),
        );
    }
    let engine = Arc::new(RwLock::new(GentleEngine::from_state(state)));
    let mut app = GENtleApp::default();
    app.engine = engine.clone();
    for (index, seq_id) in seq_ids.iter().enumerate() {
        let window = Window::new_dna(
            DNAsequence::from_sequence("ACGTACGT").expect("sequence"),
            (*seq_id).to_string(),
            engine.clone(),
        );
        app.windows.insert(
            egui::ViewportId::from_hash_of(format!("test-window-{index}-{seq_id}")),
            Arc::new(RwLock::new(window)),
        );
    }
    app
}

fn wait_for_tutorial_project_task(app: &mut GENtleApp) {
    let ctx = egui::Context::default();
    let deadline = Instant::now() + Duration::from_secs(30);
    while app.tutorial_project_task.is_some() && Instant::now() < deadline {
        app.poll_tutorial_project_task(&ctx);
        std::thread::sleep(Duration::from_millis(10));
    }
    assert!(
        app.tutorial_project_task.is_none(),
        "tutorial task did not finish in time: {}",
        app.tutorial_project_status
    );
}

fn make_test_gibson_routine_row() -> CloningRoutineCatalogRow {
    CloningRoutineCatalogRow {
        routine_id: "gibson.two_fragment_overlap_preview".to_string(),
        title: "Gibson test routine".to_string(),
        family: "gibson".to_string(),
        status: "experimental".to_string(),
        template_name: "gibson_two_fragment_overlap_preview".to_string(),
        input_ports: vec![
            serde_json::json!({
                "port_id": "left_seq_id",
                "kind": "sequence",
                "required": true
            }),
            serde_json::json!({
                "port_id": "right_seq_id",
                "kind": "sequence",
                "required": true
            }),
        ],
        ..Default::default()
    }
}

#[test]
fn preferred_agent_quickstart_helpers_pick_expected_templates() {
    let systems = vec![
        test_agent_system("builtin_echo", AgentSystemTransport::BuiltinEcho),
        test_agent_system(
            "msty_local_compat_template",
            AgentSystemTransport::NativeOpenaiCompat,
        ),
        test_agent_system(
            "anthropic_claude_sonnet_native",
            AgentSystemTransport::NativeAnthropic,
        ),
        test_agent_system("mistral_large_native", AgentSystemTransport::NativeMistral),
        test_agent_system("openai_gpt5_native", AgentSystemTransport::NativeOpenai),
    ];
    assert_eq!(
        preferred_openai_agent_system_id(&systems).as_deref(),
        Some("openai_gpt5_native")
    );
    assert_eq!(
        preferred_anthropic_agent_system_id(&systems).as_deref(),
        Some("anthropic_claude_sonnet_native")
    );
    assert_eq!(
        preferred_mistral_agent_system_id(&systems).as_deref(),
        Some("mistral_large_native")
    );
    assert_eq!(
        preferred_local_agent_system_id(&systems).as_deref(),
        Some("msty_local_compat_template")
    );
}

#[test]
fn selected_agent_session_env_overrides_use_anthropic_key_for_claude() {
    let mut app = GENtleApp::default();
    let system = test_agent_system(
        "anthropic_claude_sonnet_native",
        AgentSystemTransport::NativeAnthropic,
    );
    app.agent_openai_api_key = "sk-ant-test".to_string();
    app.agent_base_url_override = "https://api.anthropic.com/v1".to_string();
    app.agent_model_override = "claude-sonnet-4-6".to_string();
    let overrides = app
        .selected_agent_session_env_overrides(&system)
        .expect("agent overrides");
    assert_eq!(
        overrides.get(ANTHROPIC_API_KEY_ENV).map(String::as_str),
        Some("sk-ant-test")
    );
    assert_eq!(overrides.get(OPENAI_API_KEY_ENV), None);
    assert_eq!(
        overrides.get(AGENT_BASE_URL_ENV).map(String::as_str),
        Some("https://api.anthropic.com/v1")
    );
    assert_eq!(
        overrides.get(AGENT_MODEL_ENV).map(String::as_str),
        Some("claude-sonnet-4-6")
    );
}

#[test]
fn selected_agent_session_env_overrides_use_mistral_key_for_mistral() {
    let mut app = GENtleApp::default();
    let system = test_agent_system("mistral_large_native", AgentSystemTransport::NativeMistral);
    app.agent_openai_api_key = "mistral-test-key".to_string();
    app.agent_base_url_override = "https://api.mistral.ai/v1".to_string();
    app.agent_model_override = "mistral-large-latest".to_string();
    let overrides = app
        .selected_agent_session_env_overrides(&system)
        .expect("agent overrides");
    assert_eq!(
        overrides.get(MISTRAL_API_KEY_ENV).map(String::as_str),
        Some("mistral-test-key")
    );
    assert_eq!(overrides.get(OPENAI_API_KEY_ENV), None);
    assert_eq!(overrides.get(ANTHROPIC_API_KEY_ENV), None);
    assert_eq!(
        overrides.get(AGENT_BASE_URL_ENV).map(String::as_str),
        Some("https://api.mistral.ai/v1")
    );
    assert_eq!(
        overrides.get(AGENT_MODEL_ENV).map(String::as_str),
        Some("mistral-large-latest")
    );
}

#[test]
fn selected_agent_session_env_overrides_validate_and_include_gui_key() {
    let mut app = GENtleApp::default();
    let system = test_agent_system("openai_gpt5_native", AgentSystemTransport::NativeOpenai);
    app.agent_openai_api_key = "sk-test".to_string();
    app.agent_base_url_override = "https://api.openai.com/v1".to_string();
    app.agent_model_override = "gpt-5".to_string();
    app.agent_timeout_secs = "120".to_string();
    app.agent_connect_timeout_secs = "11".to_string();
    app.agent_read_timeout_secs = "180".to_string();
    app.agent_max_retries = "3".to_string();
    app.agent_max_response_bytes = "2048".to_string();
    let overrides = app
        .selected_agent_session_env_overrides(&system)
        .expect("agent overrides");
    assert_eq!(
        overrides.get(OPENAI_API_KEY_ENV).map(String::as_str),
        Some("sk-test")
    );
    assert_eq!(
        overrides.get(AGENT_BASE_URL_ENV).map(String::as_str),
        Some("https://api.openai.com/v1")
    );
    assert_eq!(
        overrides.get(AGENT_MODEL_ENV).map(String::as_str),
        Some("gpt-5")
    );
    assert_eq!(
        overrides.get(AGENT_TIMEOUT_SECS_ENV).map(String::as_str),
        Some("120")
    );
    assert_eq!(
        overrides
            .get(AGENT_CONNECT_TIMEOUT_SECS_ENV)
            .map(String::as_str),
        Some("11")
    );
    assert_eq!(
        overrides
            .get(AGENT_READ_TIMEOUT_SECS_ENV)
            .map(String::as_str),
        Some("180")
    );
    assert_eq!(
        overrides.get(AGENT_MAX_RETRIES_ENV).map(String::as_str),
        Some("3")
    );
    assert_eq!(
        overrides
            .get(AGENT_MAX_RESPONSE_BYTES_ENV)
            .map(String::as_str),
        Some("2048")
    );
}

#[test]
fn selected_agent_session_env_overrides_reject_invalid_timeout() {
    let mut app = GENtleApp::default();
    let system = test_agent_system("openai_gpt5_native", AgentSystemTransport::NativeOpenai);
    app.agent_timeout_secs = "twelve".to_string();
    let err = app
        .selected_agent_session_env_overrides(&system)
        .expect_err("invalid timeout should fail");
    assert!(err.contains("Invalid timeout_sec"));
}

#[test]
fn agent_test_setup_uses_live_for_native_transports_only() {
    let native = test_agent_system("openai_gpt5_native", AgentSystemTransport::NativeOpenai);
    let anthropic = test_agent_system(
        "anthropic_claude_sonnet_native",
        AgentSystemTransport::NativeAnthropic,
    );
    let mistral = test_agent_system("mistral_large_native", AgentSystemTransport::NativeMistral);
    let compat = test_agent_system(
        "msty_local_compat_template",
        AgentSystemTransport::NativeOpenaiCompat,
    );
    let echo = test_agent_system("builtin_echo", AgentSystemTransport::BuiltinEcho);
    assert!(GENtleApp::agent_test_setup_uses_live_probe(&native));
    assert!(GENtleApp::agent_test_setup_uses_live_probe(&anthropic));
    assert!(GENtleApp::agent_test_setup_uses_live_probe(&mistral));
    assert!(GENtleApp::agent_test_setup_uses_live_probe(&compat));
    assert!(!GENtleApp::agent_test_setup_uses_live_probe(&echo));
}

#[test]
fn agent_preflight_next_actions_explain_anthropic_auth_tokens() {
    let preflight = crate::agent_transport::AgentSystemPreflight {
        transport: AgentSystemTransport::NativeAnthropic.as_str().to_string(),
        live_probe: Some(crate::agent_transport::AgentSystemLiveProbe {
            enabled: true,
            status_class: crate::agent_transport::AgentLiveProbeStatusClass::AuthFailed,
            ..Default::default()
        }),
        ..Default::default()
    };

    let actions = GENtleApp::agent_preflight_next_actions(&preflight);
    assert_eq!(actions, vec![ANTHROPIC_API_KEY_AUTH_HINT.to_string()]);
}

#[test]
fn agent_preflight_next_actions_explain_mistral_auth_tokens() {
    let preflight = crate::agent_transport::AgentSystemPreflight {
        transport: AgentSystemTransport::NativeMistral.as_str().to_string(),
        live_probe: Some(crate::agent_transport::AgentSystemLiveProbe {
            enabled: true,
            status_class: crate::agent_transport::AgentLiveProbeStatusClass::AuthFailed,
            ..Default::default()
        }),
        ..Default::default()
    };

    let actions = GENtleApp::agent_preflight_next_actions(&preflight);
    assert_eq!(actions, vec![MISTRAL_API_KEY_AUTH_HINT.to_string()]);
}

#[test]
fn agent_setup_input_change_invalidates_stale_preflight() {
    let mut app = GENtleApp::default();
    app.agent_preflight_output = Some(crate::agent_transport::AgentSystemPreflight::default());
    app.invalidate_agent_preflight_after_setup_input_change();
    assert!(app.agent_preflight_output.is_none());
}

#[test]
fn agent_model_discovery_uses_env_key_when_session_key_is_empty() {
    let _lock = crate::genomes::genbank_env_lock()
        .lock()
        .unwrap_or_else(|e| e.into_inner());
    let _guard = EnvVarGuard::set(OPENAI_API_KEY_ENV, "sk-test-from-env");
    let app = GENtleApp::default();
    let system = test_agent_system("openai_gpt5_native", AgentSystemTransport::NativeOpenai);

    assert_eq!(
        app.selected_agent_model_discovery_key_label(&system),
        "env-openai-api-key"
    );
}

#[test]
fn agent_model_discovery_labels_anthropic_env_key() {
    let _lock = crate::genomes::genbank_env_lock()
        .lock()
        .unwrap_or_else(|e| e.into_inner());
    let _guard = EnvVarGuard::set(ANTHROPIC_API_KEY_ENV, "sk-ant-test-from-env");
    let app = GENtleApp::default();
    let system = test_agent_system(
        "anthropic_claude_sonnet_native",
        AgentSystemTransport::NativeAnthropic,
    );

    assert_eq!(
        app.selected_agent_model_discovery_key_label(&system),
        "env-anthropic-api-key"
    );
}

#[test]
fn agent_model_discovery_labels_mistral_env_key() {
    let _lock = crate::genomes::genbank_env_lock()
        .lock()
        .unwrap_or_else(|e| e.into_inner());
    let _guard = EnvVarGuard::set(MISTRAL_API_KEY_ENV, "mistral-test-from-env");
    let app = GENtleApp::default();
    let system = test_agent_system("mistral_large_native", AgentSystemTransport::NativeMistral);

    assert_eq!(
        app.selected_agent_model_discovery_key_label(&system),
        "env-mistral-api-key"
    );
}

#[test]
fn agent_model_discovery_does_not_auto_retry_failed_source() {
    let mut app = GENtleApp::default();
    let system = test_agent_system(
        "msty_local_compat_template",
        AgentSystemTransport::NativeOpenaiCompat,
    );
    app.agent_base_url_override = "http://127.0.0.1:9/v1".to_string();
    let source_key = app
        .selected_agent_model_discovery_source_key(&system)
        .expect("source key");
    app.agent_model_discovery_failed_source_key = source_key;
    app.agent_model_discovery_status =
        "Model discovery failed after 0.1s: unauthorized. Authentication failed.".to_string();

    app.start_agent_model_discovery_task(&system, false);

    assert!(app.agent_model_discovery_task.is_none());
    assert!(
        app.agent_model_discovery_status
            .contains("Authentication failed")
    );
}

#[test]
fn agent_model_discovery_auth_failure_hint_mentions_api_key_type() {
    let hint = GENtleApp::agent_model_discovery_failure_hint(
            "model discovery failed at https://api.openai.com/v1/models (status=401 Unauthorized): invalid_api_key",
        )
        .expect("auth hint");

    assert!(hint.contains("OpenAI Platform API key"));
    assert!(hint.contains("ChatGPT/Codex subscription tokens"));
}

#[test]
fn agent_model_discovery_auth_failure_hint_mentions_anthropic_console_key() {
    let hint = GENtleApp::agent_model_discovery_failure_hint(
            "Anthropic model discovery failed at https://api.anthropic.com/v1/models (status=401 Unauthorized): authentication_error invalid x-api-key",
        )
        .expect("auth hint");

    assert!(hint.contains("Anthropic Console API key"));
    assert!(hint.contains("Claude Code/Claude.ai"));
}

#[test]
fn agent_model_discovery_auth_failure_hint_mentions_mistral_api_key() {
    let hint = GENtleApp::agent_model_discovery_failure_hint(
            "Mistral model discovery failed at https://api.mistral.ai/v1/models (status=401 Unauthorized): authentication_error",
        )
        .expect("auth hint");

    assert!(hint.contains("Mistral La Plateforme API key"));
    assert!(hint.contains("Le Chat"));
}

#[test]
fn agent_assistant_suggestions_block_nested_agent_family_commands() {
    let mut app = GENtleApp::default();
    for (idx, command) in [
        "agents ask builtin_echo --prompt 'ask: capabilities'",
        "agents plan builtin_echo --prompt 'ask: capabilities'",
        "agents execute-plan @plan.json --candidate-id candidate-1",
    ]
    .iter()
    .enumerate()
    {
        app.execute_agent_suggested_command(idx + 1, command, "manual");
        assert!(
            app.agent_status.contains("agent-to-agent"),
            "unexpected status for {command}: {}",
            app.agent_status
        );
        let entry = app
            .agent_execution_log
            .last()
            .expect("blocked command should be logged");
        assert!(!entry.ok);
        assert_eq!(entry.command, *command);
        assert_eq!(entry.summary, "agent-to-agent agents command blocked");
    }
}

#[test]
fn external_agent_mcp_snippet_includes_binary_and_state_path() {
    let mut app = GENtleApp::default();
    assert_eq!(app.external_agent_mcp_state_path(), ".gentle_state.json");
    let default_snippet = app.external_agent_mcp_command_snippet();
    assert!(default_snippet.contains("gentle_mcp"));
    assert!(default_snippet.contains(".gentle_state.json"));

    app.current_project_path = Some("/tmp/GENtle Project/demo.gentle.json".to_string());
    let saved_snippet = app.external_agent_mcp_command_snippet();
    assert!(saved_snippet.contains("gentle_mcp"));
    assert!(saved_snippet.contains("'/tmp/GENtle Project/demo.gentle.json'"));
}

struct EnvVarGuard {
    key: &'static str,
    previous: Option<String>,
}

impl EnvVarGuard {
    fn set(key: &'static str, value: &str) -> Self {
        let previous = env::var(key).ok();
        unsafe {
            env::set_var(key, value);
        }
        Self { key, previous }
    }
}

impl Drop for EnvVarGuard {
    fn drop(&mut self) {
        match &self.previous {
            Some(value) => unsafe {
                env::set_var(self.key, value);
            },
            None => unsafe {
                env::remove_var(self.key);
            },
        }
    }
}

#[test]
fn load_help_doc_returns_fallback_when_missing() {
    let loaded = GENtleApp::load_help_doc(
        "/definitely/missing/gentle/help-doc.md",
        "fallback help markdown",
    );
    assert_eq!(loaded, "fallback help markdown");
}

#[test]
fn suggest_chromosome_names_prefers_alias_equivalent_contigs() {
    let available = vec![
        "NC_000017.11".to_string(),
        "chr17".to_string(),
        "17_random".to_string(),
        "chr1".to_string(),
    ];
    let suggestions = GENtleApp::suggest_chromosome_names("17", &available, 8);
    assert!(!suggestions.is_empty());
    assert_eq!(suggestions[0], "chr17");
    assert!(suggestions.iter().any(|name| name == "NC_000017.11"));
}

#[test]
fn suggest_chromosome_names_does_not_offer_shorter_numeric_prefix_only_match() {
    let available = vec!["1".to_string()];
    let suggestions = GENtleApp::suggest_chromosome_names("17", &available, 8);
    assert!(suggestions.is_empty());
}

#[test]
fn load_lineage_workspace_derives_main_split_fraction_from_legacy_heights() {
    let mut workspace = PersistedLineageGraphWorkspace::default();
    workspace.main_split_fraction = None;
    workspace.graph_area_height = 900.0;
    workspace.container_area_height = 300.0;

    let mut state = ProjectState::default();
    state.metadata.insert(
        LINEAGE_GRAPH_WORKSPACE_METADATA_KEY.to_string(),
        serde_json::to_value(workspace).expect("serialize lineage workspace"),
    );

    let mut app = GENtleApp::default();
    app.engine = Arc::new(RwLock::new(GentleEngine::from_state(state)));
    app.load_lineage_graph_workspace_from_state();

    assert!((app.lineage_main_split_fraction - 0.75).abs() <= 0.0001);
}

#[test]
fn lineage_workspace_defaults_to_graph_view_without_metadata() {
    let mut app = GENtleApp::default();
    app.lineage_graph_view = false;
    app.load_lineage_graph_workspace_from_state();

    assert!(app.lineage_graph_view);
}

#[test]
fn load_lineage_workspace_restores_persisted_table_view() {
    let workspace = PersistedLineageGraphWorkspace {
        graph_view: Some(false),
        ..PersistedLineageGraphWorkspace::default()
    };

    let mut state = ProjectState::default();
    state.metadata.insert(
        LINEAGE_GRAPH_WORKSPACE_METADATA_KEY.to_string(),
        serde_json::to_value(workspace).expect("serialize lineage workspace"),
    );

    let mut app = GENtleApp::default();
    app.engine = Arc::new(RwLock::new(GentleEngine::from_state(state)));
    app.load_lineage_graph_workspace_from_state();

    assert!(!app.lineage_graph_view);
}

#[test]
fn persist_lineage_workspace_stores_main_split_fraction() {
    let mut app = GENtleApp::default();
    app.lineage_main_split_fraction = 0.72;
    app.persist_lineage_graph_workspace_to_state();

    let workspace_value = app
        .engine
        .read()
        .unwrap()
        .state()
        .metadata
        .get(LINEAGE_GRAPH_WORKSPACE_METADATA_KEY)
        .cloned()
        .expect("workspace metadata");
    let workspace: PersistedLineageGraphWorkspace =
        serde_json::from_value(workspace_value).expect("deserialize lineage workspace");

    assert_eq!(workspace.main_split_fraction, Some(0.72));
}

#[test]
fn persist_lineage_workspace_stores_explicit_table_view_choice() {
    let mut app = GENtleApp::default();
    app.lineage_graph_view = false;
    app.persist_lineage_graph_workspace_to_state();

    let workspace_value = app
        .engine
        .read()
        .unwrap()
        .state()
        .metadata
        .get(LINEAGE_GRAPH_WORKSPACE_METADATA_KEY)
        .cloned()
        .expect("workspace metadata");
    let workspace: PersistedLineageGraphWorkspace =
        serde_json::from_value(workspace_value).expect("deserialize lineage workspace");

    assert_eq!(workspace.graph_view, Some(false));
}

#[test]
fn project_overview_metrics_count_and_route_project_areas() {
    let mut sequence = make_lineage_row("n_seq", "seq");
    sequence.pool_size = 1;
    let mut pool = make_lineage_row("n_pool", "pool");
    pool.pool_size = 3;
    let mut analysis = make_lineage_row("n_analysis", "analysis");
    analysis.kind = LineageNodeKind::Analysis;
    let rows = vec![sequence, pool, analysis];

    let metrics = GENtleApp::project_overview_metrics(&rows, 2, 1);

    assert_eq!(
        metrics,
        [
            ProjectOverviewMetric {
                label: "sequences",
                count: 1,
                target: ProjectOverviewTarget::Lineage,
                hover: "Focus the lineage graph/table sequence nodes",
            },
            ProjectOverviewMetric {
                label: "pools",
                count: 1,
                target: ProjectOverviewTarget::Lineage,
                hover: "Focus the lineage graph/table pool nodes",
            },
            ProjectOverviewMetric {
                label: "analyses",
                count: 1,
                target: ProjectOverviewTarget::Lineage,
                hover: "Focus the lineage graph/table analysis artifacts",
            },
            ProjectOverviewMetric {
                label: "containers",
                count: 2,
                target: ProjectOverviewTarget::Containers,
                hover: "Focus the container section below the lineage view",
            },
            ProjectOverviewMetric {
                label: "arrangements",
                count: 1,
                target: ProjectOverviewTarget::Arrangements,
                hover: "Focus the arrangement section below the lineage view",
            },
        ]
    );
}

#[test]
fn project_overview_focus_targets_adjust_lineage_workspace() {
    let mut app = GENtleApp::default();

    app.focus_project_overview_target(ProjectOverviewTarget::Containers);
    assert!(app.lineage_main_split_fraction < DEFAULT_LINEAGE_MAIN_SPLIT_FRACTION);
    assert!(app.lineage_container_arrangement_split_fraction > 0.5);

    app.focus_project_overview_target(ProjectOverviewTarget::Arrangements);
    assert!(app.lineage_main_split_fraction < DEFAULT_LINEAGE_MAIN_SPLIT_FRACTION);
    assert!(app.lineage_container_arrangement_split_fraction < 0.5);

    app.lineage_graph_view = false;
    app.focus_project_overview_target(ProjectOverviewTarget::Lineage);
    assert!(app.lineage_graph_view);
    assert!(app.lineage_main_split_fraction > DEFAULT_LINEAGE_MAIN_SPLIT_FRACTION);
}

#[test]
fn push_job_event_persists_background_job_history_metadata() {
    let mut app = GENtleApp::default();
    app.next_background_job_id = 9;
    app.push_job_event(
        BackgroundJobKind::PrepareGenome,
        BackgroundJobEventPhase::Started,
        Some(8),
        "prepare started",
    );

    let history = app
        .engine
        .read()
        .unwrap()
        .state()
        .metadata
        .get(BACKGROUND_JOB_HISTORY_METADATA_KEY)
        .cloned()
        .expect("background job history metadata");
    assert_eq!(
        history.get("schema").and_then(|value| value.as_str()),
        Some(BACKGROUND_JOB_HISTORY_SCHEMA)
    );
    assert_eq!(
        history
            .get("next_background_job_id")
            .and_then(|value| value.as_u64()),
        Some(9)
    );
    let events = history
        .get("events")
        .and_then(|value| value.as_array())
        .expect("events array");
    assert_eq!(events.len(), 1);
    let first = &events[0];
    assert_eq!(
        first.get("kind").and_then(|value| value.as_str()),
        Some("prepare_genome")
    );
    assert_eq!(
        first.get("phase").and_then(|value| value.as_str()),
        Some("started")
    );
    assert_eq!(
        first.get("job_id").and_then(|value| value.as_u64()),
        Some(8)
    );
    assert_eq!(
        first.get("summary").and_then(|value| value.as_str()),
        Some("prepare started")
    );
}

#[test]
fn load_background_job_history_from_state_restores_events_and_counter() {
    let mut state = ProjectState::default();
    state.metadata.insert(
        BACKGROUND_JOB_HISTORY_METADATA_KEY.to_string(),
        serde_json::json!({
            "schema": BACKGROUND_JOB_HISTORY_SCHEMA,
            "next_background_job_id": 41,
            "events": [
                {
                    "job_id": 38,
                    "kind": "prepare_genome",
                    "phase": "started",
                    "emitted_at_unix_ms": 1000,
                    "summary": "prepare started"
                },
                {
                    "job_id": 39,
                    "kind": "blast_genome",
                    "phase": "failed",
                    "emitted_at_unix_ms": 1200,
                    "summary": "blast failed"
                }
            ]
        }),
    );

    let mut app = GENtleApp::default();
    app.engine = Arc::new(RwLock::new(GentleEngine::from_state(state)));
    app.load_background_job_history_from_state();

    assert_eq!(app.next_background_job_id, 41);
    assert_eq!(app.job_event_log.len(), 2);
    assert_eq!(app.job_event_log[0].job_id, Some(38));
    assert_eq!(app.job_event_log[0].kind, BackgroundJobKind::PrepareGenome);
    assert_eq!(app.job_event_log[0].phase, BackgroundJobEventPhase::Started);
    assert_eq!(app.job_event_log[1].job_id, Some(39));
    assert_eq!(app.job_event_log[1].kind, BackgroundJobKind::BlastGenome);
    assert_eq!(app.job_event_log[1].phase, BackgroundJobEventPhase::Failed);
    assert_eq!(app.next_retry_snapshot_id, 1);
    assert!(app.retry_argument_snapshots.is_empty());
    assert_eq!(app.next_retry_cleanup_audit_id, 1);
    assert!(app.retry_snapshot_cleanup_audit.is_empty());
}

#[test]
fn load_background_job_history_from_state_restores_retry_snapshots() {
    let mut state = ProjectState::default();
    state.metadata.insert(
        BACKGROUND_JOB_HISTORY_METADATA_KEY.to_string(),
        serde_json::json!({
            "schema": BACKGROUND_JOB_HISTORY_SCHEMA,
            "next_background_job_id": 3,
            "events": [],
            "next_retry_snapshot_id": 12,
            "retry_snapshots": [
                {
                    "snapshot_id": 9,
                    "kind": "prepare_genome",
                    "origin": "background jobs panel",
                    "captured_at_unix_ms": 2000,
                    "arguments": {
                        "schema": "gentle.gui_retry_snapshot.prepare.v1",
                        "genome_id": "hg38"
                    }
                },
                {
                    "snapshot_id": 10,
                    "kind": "blast_genome",
                    "origin": "background jobs panel",
                    "captured_at_unix_ms": 2010,
                    "arguments": {
                        "schema": "gentle.gui_retry_snapshot.blast.v1",
                        "genome_id": "hg38",
                        "source_mode": "manual"
                    }
                }
            ]
        }),
    );

    let mut app = GENtleApp::default();
    app.engine = Arc::new(RwLock::new(GentleEngine::from_state(state)));
    app.load_background_job_history_from_state();

    assert_eq!(app.next_retry_snapshot_id, 12);
    assert_eq!(app.retry_argument_snapshots.len(), 2);
    assert_eq!(app.retry_argument_snapshots[0].snapshot_id, 9);
    assert_eq!(
        app.retry_argument_snapshots[0]
            .arguments
            .get("schema")
            .and_then(|value| value.as_str()),
        Some("gentle.gui_retry_snapshot.prepare.v1")
    );
    assert_eq!(app.retry_argument_snapshots[1].snapshot_id, 10);
    assert_eq!(
        app.retry_argument_snapshots[1].kind,
        BackgroundJobKind::BlastGenome
    );
    assert_eq!(app.next_retry_cleanup_audit_id, 1);
    assert!(app.retry_snapshot_cleanup_audit.is_empty());
}

#[test]
fn capture_retry_snapshot_persists_background_job_history_metadata() {
    let mut app = GENtleApp::default();
    let snapshot_id = app.capture_retry_snapshot(
        BackgroundJobKind::TrackImport,
        "background jobs panel",
        serde_json::json!({
            "schema": "gentle.gui_retry_snapshot.track_import.v1",
            "selected_sequence_id": "seq_1"
        }),
    );
    assert_eq!(snapshot_id, 1);

    let history = app
        .engine
        .read()
        .unwrap()
        .state()
        .metadata
        .get(BACKGROUND_JOB_HISTORY_METADATA_KEY)
        .cloned()
        .expect("background job history metadata");
    assert_eq!(
        history
            .get("next_retry_snapshot_id")
            .and_then(|value| value.as_u64()),
        Some(2)
    );
    let snapshots = history
        .get("retry_snapshots")
        .and_then(|value| value.as_array())
        .expect("retry snapshots array");
    assert_eq!(snapshots.len(), 1);
    assert_eq!(
        snapshots[0]
            .get("snapshot_id")
            .and_then(|value| value.as_u64()),
        Some(1)
    );
    assert_eq!(
        snapshots[0].get("kind").and_then(|value| value.as_str()),
        Some("track_import")
    );
}

#[test]
fn load_background_job_history_from_state_restores_retry_cleanup_audit() {
    let mut state = ProjectState::default();
    state.metadata.insert(
        BACKGROUND_JOB_HISTORY_METADATA_KEY.to_string(),
        serde_json::json!({
            "schema": BACKGROUND_JOB_HISTORY_SCHEMA,
            "next_background_job_id": 3,
            "events": [],
            "next_retry_snapshot_id": 2,
            "retry_snapshots": [],
            "next_retry_cleanup_audit_id": 8,
            "retry_cleanup_audit": [
                {
                    "audit_id": 6,
                    "action": "delete_filtered",
                    "performed_at_unix_ms": 3000,
                    "target_count": 4,
                    "removed_count": 3,
                    "retained_before": 10,
                    "retained_after": 7,
                    "filter_kind": "blast_genome",
                    "filter_text": "origin:panel",
                    "archive_path": null,
                    "summary": "deleted filtered snapshots"
                }
            ]
        }),
    );

    let mut app = GENtleApp::default();
    app.engine = Arc::new(RwLock::new(GentleEngine::from_state(state)));
    app.load_background_job_history_from_state();

    assert_eq!(app.next_retry_cleanup_audit_id, 8);
    assert_eq!(app.retry_snapshot_cleanup_audit.len(), 1);
    assert_eq!(app.retry_snapshot_cleanup_audit[0].audit_id, 6);
    assert_eq!(
        app.retry_snapshot_cleanup_audit[0].action,
        "delete_filtered"
    );
    assert_eq!(
        app.retry_snapshot_cleanup_audit[0].filter_kind,
        "blast_genome".to_string()
    );
}

#[test]
fn record_retry_snapshot_cleanup_audit_persists_background_job_history_metadata() {
    let mut app = GENtleApp::default();
    app.retry_snapshot_kind_filter = RetrySnapshotKindFilter::BlastGenome;
    app.retry_snapshot_text_filter = "origin:panel".to_string();
    let audit_id = app.record_retry_snapshot_cleanup_audit(
        "delete_filtered",
        "deleted filtered snapshots",
        4,
        3,
        10,
        7,
        None,
    );
    assert_eq!(audit_id, 1);

    let history = app
        .engine
        .read()
        .unwrap()
        .state()
        .metadata
        .get(BACKGROUND_JOB_HISTORY_METADATA_KEY)
        .cloned()
        .expect("background job history metadata");
    assert_eq!(
        history
            .get("next_retry_cleanup_audit_id")
            .and_then(|value| value.as_u64()),
        Some(2)
    );
    let audits = history
        .get("retry_cleanup_audit")
        .and_then(|value| value.as_array())
        .expect("retry cleanup audit array");
    assert_eq!(audits.len(), 1);
    assert_eq!(
        audits[0].get("action").and_then(|value| value.as_str()),
        Some("delete_filtered")
    );
    assert_eq!(
        audits[0]
            .get("filter_kind")
            .and_then(|value| value.as_str()),
        Some("blast_genome")
    );
    assert_eq!(
        audits[0]
            .get("removed_count")
            .and_then(|value| value.as_u64()),
        Some(3)
    );
}

#[test]
fn export_retry_cleanup_audit_report_to_path_writes_report_payload_without_self_audit() {
    let temp = tempdir().expect("tempdir");
    let output_path = temp.path().join("retry_cleanup_audit_report.json");
    let mut app = GENtleApp::default();
    app.retry_snapshot_kind_filter = RetrySnapshotKindFilter::BlastGenome;
    app.retry_snapshot_text_filter = "origin:panel".to_string();
    app.record_retry_snapshot_cleanup_audit(
        "delete_filtered",
        "deleted filtered snapshots",
        4,
        3,
        10,
        7,
        None,
    );
    app.record_retry_snapshot_cleanup_audit(
        "archive_delete_filtered",
        "archived and removed filtered snapshots",
        2,
        2,
        7,
        5,
        Some("/tmp/retry_archive.json"),
    );
    let audit_len_before = app.retry_snapshot_cleanup_audit.len();

    let exported = app
        .export_retry_cleanup_audit_report_to_path(output_path.as_path())
        .expect("export retry cleanup audit report");
    assert_eq!(exported, 2);
    assert_eq!(app.retry_snapshot_cleanup_audit.len(), audit_len_before);
    assert_eq!(app.next_retry_cleanup_audit_id, 3);

    let payload = fs::read_to_string(&output_path).expect("read retry cleanup audit report");
    let payload_json =
        serde_json::from_str::<serde_json::Value>(&payload).expect("parse audit report");
    assert_eq!(
        payload_json.get("schema").and_then(|value| value.as_str()),
        Some("gentle.gui_retry_cleanup_audit_report.v1")
    );
    assert_eq!(
        payload_json
            .get("filters")
            .and_then(|value| value.get("action"))
            .and_then(|value| value.as_str()),
        Some("all")
    );
    assert_eq!(
        payload_json
            .get("filters")
            .and_then(|value| value.get("text"))
            .and_then(|value| value.as_str()),
        Some("")
    );
    assert_eq!(
        payload_json
            .get("entry_count")
            .and_then(|value| value.as_u64()),
        Some(2)
    );
    assert_eq!(
        payload_json
            .get("retained_entry_count")
            .and_then(|value| value.as_u64()),
        Some(2)
    );
    assert_eq!(
        payload_json
            .get("action_counts")
            .and_then(|value| value.get("delete_filtered"))
            .and_then(|value| value.as_u64()),
        Some(1)
    );
    let entries = payload_json
        .get("entries")
        .and_then(|value| value.as_array())
        .expect("entries array");
    assert_eq!(entries.len(), 2);
}

#[test]
fn export_retry_cleanup_audit_report_to_path_rejects_empty_audit() {
    let temp = tempdir().expect("tempdir");
    let output_path = temp.path().join("retry_cleanup_audit_report.json");
    let app = GENtleApp::default();
    let err = app
        .export_retry_cleanup_audit_report_to_path(output_path.as_path())
        .expect_err("export should fail when audit is empty");
    assert!(err.contains("No retry cleanup audit entries match current filters"));
    assert!(!output_path.exists());
}

#[test]
fn export_retry_cleanup_audit_report_to_path_writes_filtered_payload() {
    let temp = tempdir().expect("tempdir");
    let output_path = temp.path().join("retry_cleanup_audit_report_filtered.json");
    let mut app = GENtleApp::default();
    app.retry_snapshot_kind_filter = RetrySnapshotKindFilter::BlastGenome;
    app.retry_snapshot_text_filter = "origin:panel args:hg38".to_string();
    app.record_retry_snapshot_cleanup_audit(
        "delete_filtered",
        "Removed blast hg38 panel snapshots",
        3,
        2,
        8,
        6,
        None,
    );
    app.retry_snapshot_kind_filter = RetrySnapshotKindFilter::TrackImport;
    app.retry_snapshot_text_filter = "origin:dialog args:seq_1".to_string();
    app.record_retry_snapshot_cleanup_audit(
        "prune_oldest",
        "Pruned old track-import entries",
        5,
        2,
        7,
        5,
        None,
    );
    app.retry_cleanup_audit_action_filter = RetryCleanupAuditActionFilter::DeleteFiltered;
    app.retry_cleanup_audit_text_filter = "kind:blast summary:hg38".to_string();

    let exported = app
        .export_retry_cleanup_audit_report_to_path(output_path.as_path())
        .expect("export filtered retry cleanup audit report");
    assert_eq!(exported, 1);

    let payload = fs::read_to_string(&output_path).expect("read audit export");
    let payload_json =
        serde_json::from_str::<serde_json::Value>(&payload).expect("parse audit export");
    assert_eq!(
        payload_json
            .get("filters")
            .and_then(|value| value.get("action"))
            .and_then(|value| value.as_str()),
        Some("delete_filtered")
    );
    assert_eq!(
        payload_json
            .get("filters")
            .and_then(|value| value.get("text"))
            .and_then(|value| value.as_str()),
        Some("kind:blast summary:hg38")
    );
    assert_eq!(
        payload_json
            .get("entry_count")
            .and_then(|value| value.as_u64()),
        Some(1)
    );
    assert_eq!(
        payload_json
            .get("retained_entry_count")
            .and_then(|value| value.as_u64()),
        Some(2)
    );
    let entries = payload_json
        .get("entries")
        .and_then(|value| value.as_array())
        .expect("entries array");
    assert_eq!(entries.len(), 1);
    assert_eq!(
        entries[0].get("action").and_then(|value| value.as_str()),
        Some("delete_filtered")
    );
}

#[test]
fn export_retry_cleanup_audit_report_to_path_rejects_empty_filtered_match_set() {
    let temp = tempdir().expect("tempdir");
    let output_path = temp.path().join("retry_cleanup_audit_report_filtered.json");
    let mut app = GENtleApp::default();
    app.record_retry_snapshot_cleanup_audit("delete_filtered", "first", 2, 1, 4, 3, None);
    app.retry_cleanup_audit_action_filter = RetryCleanupAuditActionFilter::PruneOldest;
    app.retry_cleanup_audit_text_filter = "summary:nope".to_string();

    let err = app
        .export_retry_cleanup_audit_report_to_path(output_path.as_path())
        .expect_err("export should fail when no audit entries match filters");
    assert!(err.contains("No retry cleanup audit entries match current filters"));
    assert!(!output_path.exists());
}

#[test]
fn filtered_retry_cleanup_audit_entries_support_action_and_scoped_text_filters() {
    let mut app = GENtleApp::default();
    app.retry_snapshot_kind_filter = RetrySnapshotKindFilter::BlastGenome;
    app.retry_snapshot_text_filter = "origin:panel args:hg38".to_string();
    app.record_retry_snapshot_cleanup_audit(
        "delete_filtered",
        "Removed blast hg38 panel snapshots",
        3,
        2,
        8,
        6,
        None,
    );
    app.retry_snapshot_kind_filter = RetrySnapshotKindFilter::TrackImport;
    app.retry_snapshot_text_filter = "origin:dialog args:seq_1".to_string();
    app.record_retry_snapshot_cleanup_audit(
        "prune_oldest",
        "Pruned old track-import entries",
        5,
        2,
        7,
        5,
        None,
    );

    app.retry_cleanup_audit_action_filter = RetryCleanupAuditActionFilter::DeleteFiltered;
    app.retry_cleanup_audit_text_filter = "kind:blast summary:hg38".to_string();

    let filtered = app.filtered_retry_cleanup_audit_entries();
    assert_eq!(filtered.len(), 1);
    assert_eq!(filtered[0].audit_id, 1);
    assert_eq!(filtered[0].action, "delete_filtered");
}

#[test]
fn record_retry_snapshot_cleanup_audit_honors_audit_retention_limit() {
    let mut app = GENtleApp::default();
    app.retry_cleanup_audit_retain_count = 2;

    app.record_retry_snapshot_cleanup_audit("delete_filtered", "first", 3, 1, 6, 5, None);
    app.record_retry_snapshot_cleanup_audit("prune_oldest", "second", 2, 1, 5, 4, None);
    app.record_retry_snapshot_cleanup_audit("clear_all", "third", 4, 4, 4, 0, None);

    assert_eq!(app.retry_snapshot_cleanup_audit.len(), 2);
    assert_eq!(app.retry_snapshot_cleanup_audit[0].audit_id, 2);
    assert_eq!(app.retry_snapshot_cleanup_audit[1].audit_id, 3);
}

#[test]
fn prune_retry_cleanup_audit_to_limit_removes_oldest_and_persists_metadata() {
    let mut app = GENtleApp::default();
    app.record_retry_snapshot_cleanup_audit("delete_filtered", "first", 3, 1, 6, 5, None);
    app.record_retry_snapshot_cleanup_audit("prune_oldest", "second", 2, 1, 5, 4, None);
    app.record_retry_snapshot_cleanup_audit("clear_all", "third", 4, 4, 4, 0, None);

    let removed = app.prune_retry_cleanup_audit_to_limit(2);
    assert_eq!(removed, 1);
    assert_eq!(app.retry_snapshot_cleanup_audit.len(), 2);
    assert_eq!(app.retry_snapshot_cleanup_audit[0].audit_id, 2);
    assert_eq!(app.retry_snapshot_cleanup_audit[1].audit_id, 3);

    let history = app
        .engine
        .read()
        .unwrap()
        .state()
        .metadata
        .get(BACKGROUND_JOB_HISTORY_METADATA_KEY)
        .cloned()
        .expect("background job history metadata");
    let audits = history
        .get("retry_cleanup_audit")
        .and_then(|value| value.as_array())
        .expect("retry cleanup audit array");
    assert_eq!(audits.len(), 2);
    assert_eq!(
        audits[0].get("audit_id").and_then(|value| value.as_u64()),
        Some(2)
    );
    assert_eq!(
        audits[1].get("audit_id").and_then(|value| value.as_u64()),
        Some(3)
    );
}

#[test]
fn clear_retry_cleanup_audit_removes_all_and_persists_metadata() {
    let mut app = GENtleApp::default();
    app.record_retry_snapshot_cleanup_audit("delete_filtered", "first", 3, 1, 6, 5, None);
    app.record_retry_snapshot_cleanup_audit("prune_oldest", "second", 2, 1, 5, 4, None);

    let removed = app.clear_retry_cleanup_audit();
    assert_eq!(removed, 2);
    assert!(app.retry_snapshot_cleanup_audit.is_empty());

    let history = app
        .engine
        .read()
        .unwrap()
        .state()
        .metadata
        .get(BACKGROUND_JOB_HISTORY_METADATA_KEY)
        .cloned()
        .expect("background job history metadata");
    let audits = history
        .get("retry_cleanup_audit")
        .and_then(|value| value.as_array())
        .expect("retry cleanup audit array");
    assert!(audits.is_empty());
}

#[test]
fn retry_cleanup_audit_clear_confirm_input_matches_expected_phrase_case_insensitive() {
    assert!(GENtleApp::retry_cleanup_audit_clear_confirm_input_matches(
        4,
        "  ClEaR-AuDiT 4  "
    ));
    assert!(!GENtleApp::retry_cleanup_audit_clear_confirm_input_matches(
        4,
        "clear-audit 3"
    ));
}

#[test]
fn filtered_retry_snapshots_support_kind_and_scoped_text_filters() {
    let mut app = GENtleApp::default();
    app.capture_retry_snapshot(
        BackgroundJobKind::PrepareGenome,
        "background jobs panel",
        serde_json::json!({
            "schema": "gentle.gui_retry_snapshot.prepare.v1",
            "genome_id": "hg38"
        }),
    );
    app.capture_retry_snapshot(
        BackgroundJobKind::BlastGenome,
        "background jobs panel",
        serde_json::json!({
            "schema": "gentle.gui_retry_snapshot.blast.v1",
            "genome_id": "hg38",
            "source_mode": "manual"
        }),
    );
    app.capture_retry_snapshot(
        BackgroundJobKind::BlastGenome,
        "dialog retry button",
        serde_json::json!({
            "schema": "gentle.gui_retry_snapshot.blast.v1",
            "genome_id": "mm10",
            "source_mode": "manual"
        }),
    );
    app.retry_snapshot_kind_filter = RetrySnapshotKindFilter::BlastGenome;
    app.retry_snapshot_text_filter = "origin:panel args:hg38".to_string();

    let filtered = app.filtered_retry_snapshots();
    assert_eq!(filtered.len(), 1);
    assert_eq!(filtered[0].snapshot_id, 2);
    assert_eq!(filtered[0].kind, BackgroundJobKind::BlastGenome);
    assert_eq!(filtered[0].origin, "background jobs panel");
}

#[test]
fn export_filtered_retry_snapshots_to_path_writes_filtered_payload() {
    let temp = tempdir().expect("tempdir");
    let output_path = temp.path().join("retry_snapshots.json");
    let mut app = GENtleApp::default();
    app.capture_retry_snapshot(
        BackgroundJobKind::PrepareGenome,
        "background jobs panel",
        serde_json::json!({
            "schema": "gentle.gui_retry_snapshot.prepare.v1",
            "genome_id": "hg38"
        }),
    );
    app.capture_retry_snapshot(
        BackgroundJobKind::BlastGenome,
        "background jobs panel",
        serde_json::json!({
            "schema": "gentle.gui_retry_snapshot.blast.v1",
            "genome_id": "hg38",
            "source_mode": "manual"
        }),
    );
    app.capture_retry_snapshot(
        BackgroundJobKind::BlastGenome,
        "background jobs panel",
        serde_json::json!({
            "schema": "gentle.gui_retry_snapshot.blast.v1",
            "genome_id": "mm10",
            "source_mode": "project_sequence"
        }),
    );
    app.retry_snapshot_kind_filter = RetrySnapshotKindFilter::BlastGenome;
    app.retry_snapshot_text_filter = "args:hg38".to_string();

    let count = app
        .export_filtered_retry_snapshots_to_path(output_path.as_path())
        .expect("retry snapshot export");
    assert_eq!(count, 1);

    let payload = fs::read_to_string(&output_path).expect("read export payload");
    let payload_json =
        serde_json::from_str::<serde_json::Value>(&payload).expect("parse export payload");
    assert_eq!(
        payload_json.get("schema").and_then(|value| value.as_str()),
        Some("gentle.gui_retry_snapshot_export.v1")
    );
    assert_eq!(
        payload_json
            .get("filters")
            .and_then(|value| value.get("kind"))
            .and_then(|value| value.as_str()),
        Some("blast_genome")
    );
    assert_eq!(
        payload_json
            .get("filters")
            .and_then(|value| value.get("text"))
            .and_then(|value| value.as_str()),
        Some("args:hg38")
    );
    assert_eq!(
        payload_json
            .get("snapshot_count")
            .and_then(|value| value.as_u64()),
        Some(1)
    );
    let snapshots = payload_json
        .get("snapshots")
        .and_then(|value| value.as_array())
        .expect("snapshots array");
    assert_eq!(snapshots.len(), 1);
    assert_eq!(
        snapshots[0].get("kind").and_then(|value| value.as_str()),
        Some("blast_genome")
    );
}

#[test]
fn export_filtered_retry_snapshots_to_path_rejects_empty_match_set() {
    let temp = tempdir().expect("tempdir");
    let output_path = temp.path().join("retry_snapshots.json");
    let mut app = GENtleApp::default();
    app.capture_retry_snapshot(
        BackgroundJobKind::PrepareGenome,
        "background jobs panel",
        serde_json::json!({
            "schema": "gentle.gui_retry_snapshot.prepare.v1",
            "genome_id": "hg38"
        }),
    );
    app.retry_snapshot_kind_filter = RetrySnapshotKindFilter::BlastGenome;
    app.retry_snapshot_text_filter = "args:mm10".to_string();

    let err = app
        .export_filtered_retry_snapshots_to_path(output_path.as_path())
        .expect_err("export should fail when no snapshots match filters");
    assert!(err.contains("No retry snapshots match current filters"));
    assert!(!output_path.exists());
}

#[test]
fn capture_retry_snapshot_honors_retention_limit() {
    let mut app = GENtleApp::default();
    app.retry_snapshot_retain_count = 2;
    app.capture_retry_snapshot(
        BackgroundJobKind::PrepareGenome,
        "background jobs panel",
        serde_json::json!({
            "schema": "gentle.gui_retry_snapshot.prepare.v1",
            "genome_id": "hg38"
        }),
    );
    app.capture_retry_snapshot(
        BackgroundJobKind::BlastGenome,
        "background jobs panel",
        serde_json::json!({
            "schema": "gentle.gui_retry_snapshot.blast.v1",
            "genome_id": "hg38"
        }),
    );
    app.capture_retry_snapshot(
        BackgroundJobKind::TrackImport,
        "background jobs panel",
        serde_json::json!({
            "schema": "gentle.gui_retry_snapshot.track_import.v1",
            "selected_sequence_id": "seq_1"
        }),
    );
    assert_eq!(app.retry_argument_snapshots.len(), 2);
    assert_eq!(app.retry_argument_snapshots[0].snapshot_id, 2);
    assert_eq!(app.retry_argument_snapshots[1].snapshot_id, 3);
}

#[test]
fn prune_retry_snapshots_to_limit_removes_oldest_and_persists_metadata() {
    let mut app = GENtleApp::default();
    app.capture_retry_snapshot(
        BackgroundJobKind::PrepareGenome,
        "background jobs panel",
        serde_json::json!({"schema": "gentle.gui_retry_snapshot.prepare.v1"}),
    );
    app.capture_retry_snapshot(
        BackgroundJobKind::BlastGenome,
        "background jobs panel",
        serde_json::json!({"schema": "gentle.gui_retry_snapshot.blast.v1"}),
    );
    app.capture_retry_snapshot(
        BackgroundJobKind::TrackImport,
        "background jobs panel",
        serde_json::json!({"schema": "gentle.gui_retry_snapshot.track_import.v1"}),
    );

    let removed = app.prune_retry_snapshots_to_limit(2);
    assert_eq!(removed, 1);
    assert_eq!(app.retry_argument_snapshots.len(), 2);
    assert_eq!(app.retry_argument_snapshots[0].snapshot_id, 2);
    assert_eq!(app.retry_argument_snapshots[1].snapshot_id, 3);

    let history = app
        .engine
        .read()
        .unwrap()
        .state()
        .metadata
        .get(BACKGROUND_JOB_HISTORY_METADATA_KEY)
        .cloned()
        .expect("background job history metadata");
    let snapshots = history
        .get("retry_snapshots")
        .and_then(|value| value.as_array())
        .expect("retry snapshots array");
    assert_eq!(snapshots.len(), 2);
    assert_eq!(
        snapshots[0]
            .get("snapshot_id")
            .and_then(|value| value.as_u64()),
        Some(2)
    );
    assert_eq!(
        snapshots[1]
            .get("snapshot_id")
            .and_then(|value| value.as_u64()),
        Some(3)
    );
}

#[test]
fn clear_retry_snapshots_removes_all_and_persists_metadata() {
    let mut app = GENtleApp::default();
    app.capture_retry_snapshot(
        BackgroundJobKind::PrepareGenome,
        "background jobs panel",
        serde_json::json!({"schema": "gentle.gui_retry_snapshot.prepare.v1"}),
    );
    app.capture_retry_snapshot(
        BackgroundJobKind::BlastGenome,
        "background jobs panel",
        serde_json::json!({"schema": "gentle.gui_retry_snapshot.blast.v1"}),
    );

    let removed = app.clear_retry_snapshots();
    assert_eq!(removed, 2);
    assert!(app.retry_argument_snapshots.is_empty());

    let history = app
        .engine
        .read()
        .unwrap()
        .state()
        .metadata
        .get(BACKGROUND_JOB_HISTORY_METADATA_KEY)
        .cloned()
        .expect("background job history metadata");
    let snapshots = history
        .get("retry_snapshots")
        .and_then(|value| value.as_array())
        .expect("retry snapshots array");
    assert!(snapshots.is_empty());
}

#[test]
fn remove_filtered_retry_snapshots_deletes_only_matching_entries() {
    let mut app = GENtleApp::default();
    app.capture_retry_snapshot(
        BackgroundJobKind::PrepareGenome,
        "background jobs panel",
        serde_json::json!({
            "schema": "gentle.gui_retry_snapshot.prepare.v1",
            "genome_id": "hg38"
        }),
    );
    app.capture_retry_snapshot(
        BackgroundJobKind::BlastGenome,
        "background jobs panel",
        serde_json::json!({
            "schema": "gentle.gui_retry_snapshot.blast.v1",
            "genome_id": "hg38"
        }),
    );
    app.capture_retry_snapshot(
        BackgroundJobKind::BlastGenome,
        "dialog retry button",
        serde_json::json!({
            "schema": "gentle.gui_retry_snapshot.blast.v1",
            "genome_id": "mm10"
        }),
    );
    app.retry_snapshot_kind_filter = RetrySnapshotKindFilter::BlastGenome;
    app.retry_snapshot_text_filter = "origin:panel args:hg38".to_string();

    let removed = app.remove_filtered_retry_snapshots();
    assert_eq!(removed, 1);
    assert_eq!(app.retry_argument_snapshots.len(), 2);
    assert_eq!(app.retry_argument_snapshots[0].snapshot_id, 1);
    assert_eq!(app.retry_argument_snapshots[1].snapshot_id, 3);
}

#[test]
fn archive_and_delete_filtered_retry_snapshots_to_path_writes_archive_and_removes_matches() {
    let temp = tempdir().expect("tempdir");
    let output_path = temp.path().join("retry_snapshots_archive.json");
    let mut app = GENtleApp::default();
    app.capture_retry_snapshot(
        BackgroundJobKind::PrepareGenome,
        "background jobs panel",
        serde_json::json!({
            "schema": "gentle.gui_retry_snapshot.prepare.v1",
            "genome_id": "hg38"
        }),
    );
    app.capture_retry_snapshot(
        BackgroundJobKind::BlastGenome,
        "background jobs panel",
        serde_json::json!({
            "schema": "gentle.gui_retry_snapshot.blast.v1",
            "genome_id": "hg38"
        }),
    );
    app.capture_retry_snapshot(
        BackgroundJobKind::BlastGenome,
        "dialog retry button",
        serde_json::json!({
            "schema": "gentle.gui_retry_snapshot.blast.v1",
            "genome_id": "mm10"
        }),
    );
    app.retry_snapshot_kind_filter = RetrySnapshotKindFilter::BlastGenome;
    app.retry_snapshot_text_filter = "origin:panel".to_string();

    let removed = app
        .archive_and_delete_filtered_retry_snapshots_to_path(output_path.as_path())
        .expect("archive/delete retry snapshots");
    assert_eq!(removed, 1);
    assert_eq!(app.retry_argument_snapshots.len(), 2);
    assert_eq!(app.retry_argument_snapshots[0].snapshot_id, 1);
    assert_eq!(app.retry_argument_snapshots[1].snapshot_id, 3);

    let payload = fs::read_to_string(&output_path).expect("read archive payload");
    let payload_json =
        serde_json::from_str::<serde_json::Value>(&payload).expect("parse archive payload");
    assert_eq!(
        payload_json.get("schema").and_then(|value| value.as_str()),
        Some("gentle.gui_retry_snapshot_archive.v1")
    );
    assert_eq!(
        payload_json
            .get("archive_mode")
            .and_then(|value| value.as_str()),
        Some("archive_and_delete_filtered")
    );
    assert_eq!(
        payload_json
            .get("archived_count")
            .and_then(|value| value.as_u64()),
        Some(1)
    );
    let archived = payload_json
        .get("snapshots")
        .and_then(|value| value.as_array())
        .expect("archived snapshots array");
    assert_eq!(archived.len(), 1);
    assert_eq!(
        archived[0].get("kind").and_then(|value| value.as_str()),
        Some("blast_genome")
    );
}

#[test]
fn archive_and_delete_filtered_retry_snapshots_to_path_rejects_empty_match_set() {
    let temp = tempdir().expect("tempdir");
    let output_path = temp.path().join("retry_snapshots_archive.json");
    let mut app = GENtleApp::default();
    app.capture_retry_snapshot(
        BackgroundJobKind::PrepareGenome,
        "background jobs panel",
        serde_json::json!({
            "schema": "gentle.gui_retry_snapshot.prepare.v1",
            "genome_id": "hg38"
        }),
    );
    app.retry_snapshot_kind_filter = RetrySnapshotKindFilter::BlastGenome;
    app.retry_snapshot_text_filter = "args:mm10".to_string();

    let err = app
        .archive_and_delete_filtered_retry_snapshots_to_path(output_path.as_path())
        .expect_err("archive/delete should fail when no snapshots match filters");
    assert!(err.contains("No retry snapshots match current filters"));
    assert!(!output_path.exists());
}

#[test]
fn retry_snapshot_dry_run_diff_splits_removed_and_retained_by_filters() {
    let mut app = GENtleApp::default();
    app.capture_retry_snapshot(
        BackgroundJobKind::PrepareGenome,
        "background jobs panel",
        serde_json::json!({"schema": "gentle.gui_retry_snapshot.prepare.v1"}),
    );
    app.capture_retry_snapshot(
        BackgroundJobKind::BlastGenome,
        "background jobs panel",
        serde_json::json!({"schema": "gentle.gui_retry_snapshot.blast.v1", "genome_id": "hg38"}),
    );
    app.capture_retry_snapshot(
        BackgroundJobKind::BlastGenome,
        "dialog retry button",
        serde_json::json!({"schema": "gentle.gui_retry_snapshot.blast.v1", "genome_id": "mm10"}),
    );
    app.retry_snapshot_kind_filter = RetrySnapshotKindFilter::BlastGenome;
    app.retry_snapshot_text_filter = "origin:panel args:hg38".to_string();

    let diff = app.retry_snapshot_dry_run_diff();
    assert_eq!(diff.removed.len(), 1);
    assert_eq!(diff.retained.len(), 2);
    assert_eq!(diff.removed[0].snapshot_id, 2);
    assert_eq!(diff.retained[0].snapshot_id, 1);
    assert_eq!(diff.retained[1].snapshot_id, 3);
}

#[test]
fn retry_snapshot_dry_run_diff_summary_line_reports_counts() {
    let mut app = GENtleApp::default();
    app.capture_retry_snapshot(
        BackgroundJobKind::PrepareGenome,
        "background jobs panel",
        serde_json::json!({"schema": "gentle.gui_retry_snapshot.prepare.v1"}),
    );
    app.capture_retry_snapshot(
        BackgroundJobKind::BlastGenome,
        "background jobs panel",
        serde_json::json!({"schema": "gentle.gui_retry_snapshot.blast.v1"}),
    );
    app.retry_snapshot_kind_filter = RetrySnapshotKindFilter::BlastGenome;
    app.retry_snapshot_text_filter = "origin:panel".to_string();

    let diff = app.retry_snapshot_dry_run_diff();
    assert_eq!(
        diff.summary_line(),
        "Dry-run diff: remove 1 snapshot(s), retain 1"
    );
}

#[test]
fn summarize_retry_snapshot_cleanup_targets_reports_kind_origin_and_id_ranges() {
    let mut app = GENtleApp::default();
    app.capture_retry_snapshot(
        BackgroundJobKind::PrepareGenome,
        "background jobs panel",
        serde_json::json!({"schema": "gentle.gui_retry_snapshot.prepare.v1"}),
    );
    app.capture_retry_snapshot(
        BackgroundJobKind::BlastGenome,
        "background jobs panel",
        serde_json::json!({"schema": "gentle.gui_retry_snapshot.blast.v1"}),
    );
    app.capture_retry_snapshot(
        BackgroundJobKind::TrackImport,
        "dialog retry button",
        serde_json::json!({"schema": "gentle.gui_retry_snapshot.track_import.v1"}),
    );

    let summary =
        GENtleApp::summarize_retry_snapshot_cleanup_targets(&app.retry_argument_snapshots.clone());
    assert!(summary.contains("3 match(es)"));
    assert!(summary.contains("PrepareGenome=1"));
    assert!(summary.contains("BlastGenome=1"));
    assert!(summary.contains("TrackImport=1"));
    assert!(summary.contains("OpenTutorialProject=0"));
    assert!(summary.contains("AgentAssist=0"));
    assert!(summary.contains("background jobs panel=2"));
    assert!(summary.contains("dialog retry button=1"));
    assert!(summary.contains("ids #1..#3"));
}

#[test]
fn format_retry_snapshot_cleanup_status_reports_before_after_counts() {
    let text = GENtleApp::format_retry_snapshot_cleanup_status(
        RetrySnapshotPendingCleanupAction::DeleteFiltered,
        4,
        9,
        5,
    );
    assert_eq!(
        text,
        "Delete filtered confirmed: removed 4 snapshot(s), retained 5 (was 9)"
    );
}

#[test]
fn retry_snapshot_pending_cleanup_action_confirm_phrase_includes_target_count() {
    assert_eq!(
        RetrySnapshotPendingCleanupAction::DeleteFiltered.confirm_phrase(7),
        "delete 7"
    );
    assert_eq!(
        RetrySnapshotPendingCleanupAction::ArchiveAndDeleteFiltered.confirm_phrase(3),
        "archive-delete 3"
    );
}

#[test]
fn retry_snapshot_cleanup_confirm_input_matches_expected_phrase_case_insensitive() {
    assert!(GENtleApp::retry_snapshot_cleanup_confirm_input_matches(
        RetrySnapshotPendingCleanupAction::DeleteFiltered,
        4,
        "  DeLeTe 4  "
    ));
    assert!(!GENtleApp::retry_snapshot_cleanup_confirm_input_matches(
        RetrySnapshotPendingCleanupAction::DeleteFiltered,
        4,
        "delete 3"
    ));
}

#[test]
fn background_jobs_and_history_scroll_area_ids_are_unique_and_stable() {
    let ids = [
        BACKGROUND_JOBS_RECENT_JOB_EVENTS_SCROLL_ID,
        BACKGROUND_JOBS_RETRY_SNAPSHOTS_REMOVED_PREVIEW_SCROLL_ID,
        BACKGROUND_JOBS_RETRY_SNAPSHOTS_RETAINED_PREVIEW_SCROLL_ID,
        BACKGROUND_JOBS_RETRY_SNAPSHOTS_SCROLL_ID,
        BACKGROUND_JOBS_RETRY_CLEANUP_AUDIT_SCROLL_ID,
        OPERATION_HISTORY_SCROLL_ID,
    ];
    let unique: HashSet<_> = ids.into_iter().collect();
    assert_eq!(
        unique.len(),
        ids.len(),
        "duplicate auxiliary scroll id salt"
    );
    assert!(ids.iter().all(|id| id.ends_with("_scroll")));
    assert!(
        ids.iter()
            .all(|id| id.starts_with("background_jobs_") || *id == OPERATION_HISTORY_SCROLL_ID)
    );
}

#[test]
fn lineage_panel_split_clamps_for_small_height_and_round_trips_via_workspace() {
    let expected_graph_height = LINEAGE_MAIN_TOP_PANEL_MIN_HEIGHT;
    let expected_container_height = 120.0;
    let expected_split =
        expected_graph_height / (expected_graph_height + expected_container_height);
    let (graph_h, container_h, split) = GENtleApp::normalized_lineage_panel_heights(80.0, 0.95);
    assert!((graph_h - expected_graph_height).abs() <= 0.0001);
    assert!((container_h - expected_container_height).abs() <= 0.0001);
    assert!((split - expected_split).abs() <= 0.0001);

    let mut app = GENtleApp::default();
    app.lineage_graph_area_height = graph_h;
    app.lineage_container_area_height = container_h;
    app.lineage_main_split_fraction = split;
    app.persist_lineage_graph_workspace_to_state();

    let workspace_value = app
        .engine
        .read()
        .unwrap()
        .state()
        .metadata
        .get(LINEAGE_GRAPH_WORKSPACE_METADATA_KEY)
        .cloned()
        .expect("workspace metadata");

    let mut state = ProjectState::default();
    state.metadata.insert(
        LINEAGE_GRAPH_WORKSPACE_METADATA_KEY.to_string(),
        workspace_value,
    );
    let mut reloaded = GENtleApp::default();
    reloaded.engine = Arc::new(RwLock::new(GentleEngine::from_state(state)));
    reloaded.load_lineage_graph_workspace_from_state();

    assert!((reloaded.lineage_graph_area_height - expected_graph_height).abs() <= 0.0001);
    assert!((reloaded.lineage_container_area_height - expected_container_height).abs() <= 0.0001);
    assert!((reloaded.lineage_main_split_fraction - expected_split).abs() <= 0.0001);
}

#[test]
fn lineage_container_subpanel_split_clamps_for_small_height() {
    let (containers_h, arrangements_h, split) =
        GENtleApp::normalized_lineage_container_subpanel_heights(120.0, 0.95);
    assert!((containers_h - 80.0).abs() <= 0.0001);
    assert!((arrangements_h - 80.0).abs() <= 0.0001);
    assert!((split - (80.0 / 160.0)).abs() <= 0.0001);
}

#[test]
fn lineage_workspace_round_trips_container_arrangement_split_fraction() {
    let mut app = GENtleApp::default();
    app.lineage_container_arrangement_split_fraction = 0.63;
    app.persist_lineage_graph_workspace_to_state();

    let workspace_value = app
        .engine
        .read()
        .unwrap()
        .state()
        .metadata
        .get(LINEAGE_GRAPH_WORKSPACE_METADATA_KEY)
        .cloned()
        .expect("workspace metadata");
    let workspace: PersistedLineageGraphWorkspace =
        serde_json::from_value(workspace_value.clone()).expect("deserialize lineage workspace");
    assert_eq!(workspace.container_arrangement_split_fraction, Some(0.63));

    let mut state = ProjectState::default();
    state.metadata.insert(
        LINEAGE_GRAPH_WORKSPACE_METADATA_KEY.to_string(),
        workspace_value,
    );
    let mut reloaded = GENtleApp::default();
    reloaded.engine = Arc::new(RwLock::new(GentleEngine::from_state(state)));
    reloaded.load_lineage_graph_workspace_from_state();

    assert!((reloaded.lineage_container_arrangement_split_fraction - 0.63).abs() <= 0.0001);
}

#[test]
fn open_sequence_window_enqueues_lazy_window_for_existing_sequence() {
    let mut state = ProjectState::default();
    state.sequences.insert(
        "seq_lazy".to_string(),
        DNAsequence::from_sequence("ACGT").expect("sequence"),
    );
    let engine = Arc::new(RwLock::new(GentleEngine::from_state(state)));
    let mut app = GENtleApp::default();
    app.engine = engine;

    app.open_sequence_window("seq_lazy");

    assert_eq!(app.new_windows.len(), 1);
    assert_eq!(
        app.new_windows[0].sequence_id().as_deref(),
        Some("seq_lazy")
    );
}

#[test]
fn create_sequence_from_text_input_adds_project_sequence_and_window() {
    let mut app = GENtleApp::default();

    let created = app
        .create_sequence_from_text_input(
            "1 acg u\n",
            Some("typed_seq".to_string()),
            Some("Typed sequence".to_string()),
            true,
        )
        .expect("create typed sequence");

    assert_eq!(created, vec!["typed_seq".to_string()]);
    let engine = app.engine.read().unwrap();
    let dna = engine
        .state()
        .sequences
        .get("typed_seq")
        .expect("typed sequence");
    assert_eq!(dna.get_forward_string(), "ACGT");
    assert_eq!(dna.name().as_deref(), Some("Typed sequence"));
    assert!(dna.is_circular());
    assert_eq!(engine.operation_log().len(), 1);
    drop(engine);
    assert_eq!(app.new_windows.len(), 1);
    assert_eq!(
        app.new_windows[0].sequence_id().as_deref(),
        Some("typed_seq")
    );
}

#[test]
fn create_sequence_from_text_input_rejects_invalid_iupac() {
    let mut app = GENtleApp::default();

    let err = app
        .create_sequence_from_text_input("ACGZ", Some("bad_seq".to_string()), None, false)
        .expect_err("invalid base should fail before mutation");

    assert!(
        err.contains("Invalid IUPAC"),
        "unexpected error message: {err}"
    );
    let engine = app.engine.read().unwrap();
    assert!(engine.state().sequences.is_empty());
    assert!(engine.operation_log().is_empty());
    assert!(app.new_windows.is_empty());
}

#[test]
fn register_window_records_one_time_initial_position() {
    let mut app = GENtleApp::default();
    let engine = app.engine.clone();
    let window = Window::new_dna(
        DNAsequence::from_sequence("ACGT").expect("sequence"),
        "seq_a".to_string(),
        engine,
    );

    let viewport_id = app.register_window(window);

    assert_eq!(
        app.pending_window_initial_positions.get(&viewport_id),
        Some(&egui::Pos2 { x: 0.0, y: 0.0 })
    );
    assert!(
        app.pending_window_open_timestamps
            .contains_key(&viewport_id),
        "registered deferred windows should track an open probe"
    );
    assert!(
        app.pending_focus_viewports.contains(&viewport_id),
        "newly registered sequence windows should be raised after first render"
    );
    assert!(
        app.pending_viewport_focus_timestamps
            .contains_key(&viewport_id),
        "foreground focus requests should be tracked for slow-focus diagnostics"
    );
    let consumed = app.pending_window_initial_positions.remove(&viewport_id);
    assert_eq!(consumed, Some(egui::Pos2 { x: 0.0, y: 0.0 }));
    assert!(
        !app.pending_window_initial_positions
            .contains_key(&viewport_id),
        "initial position should be consumed after first deferred viewport render"
    );
}

#[test]
fn deferred_window_position_uses_monotonic_cascade_index() {
    assert_eq!(
        GENtleApp::deferred_window_position(0),
        egui::Pos2 { x: 0.0, y: 0.0 }
    );
    assert_eq!(
        GENtleApp::deferred_window_position(3),
        egui::Pos2 { x: 600.0, y: 600.0 }
    );
}

#[test]
fn deferred_window_initial_commands_emit_one_shot_outer_position() {
    assert_eq!(
        GENtleApp::deferred_window_initial_commands(Some(egui::Pos2 { x: 120.0, y: 240.0 })),
        vec![egui::ViewportCommand::OuterPosition(egui::Pos2 {
            x: 120.0,
            y: 240.0,
        })]
    );
    assert!(
        GENtleApp::deferred_window_initial_commands(None).is_empty(),
        "no initial position should mean no one-shot viewport commands"
    );
}

#[test]
fn macos_prefers_immediate_sequence_viewports() {
    if cfg!(target_os = "macos") {
        assert!(GENtleApp::use_immediate_sequence_viewports());
    } else {
        assert!(!GENtleApp::use_immediate_sequence_viewports());
    }
}

#[test]
fn default_sequence_windows_accept_native_close_requests() {
    let _lock = crate::genomes::genbank_env_lock()
        .lock()
        .unwrap_or_else(|e| e.into_inner());
    let _native_env_guard = EnvVarGuard::set(super::MACOS_NATIVE_CHILD_VIEWPORTS_ENV, "0");
    let _hosted_env_guard = EnvVarGuard::set(super::MACOS_HOSTED_CHILD_VIEWPORTS_ENV, "0");

    assert!(GENtleApp::sequence_window_accepts_native_close_request());
}

#[test]
fn open_configuration_dialog_focuses_existing_window_without_resetting_edits() {
    let mut app = GENtleApp::default();
    app.open_configuration_dialog();
    app.configuration_graphics.feature_details_font_size = 19.25;
    app.configuration_graphics_dirty = true;
    app.configuration_status = "editing".to_string();

    app.open_configuration_dialog();

    assert!(app.show_configuration_dialog);
    assert!(app.configuration_graphics_dirty);
    assert_eq!(app.configuration_graphics.feature_details_font_size, 19.25);
    assert_eq!(app.configuration_status, "editing");
    assert!(
        app.pending_focus_viewports
            .contains(&GENtleApp::configuration_viewport_id())
    );
}

fn collect_rendered_text_from_shape(shape: &egui::epaint::Shape, out: &mut Vec<String>) {
    match shape {
        egui::epaint::Shape::Text(text) => {
            let raw = text.galley.job.text.trim();
            if !raw.is_empty() {
                out.push(raw.to_string());
            }
        }
        egui::epaint::Shape::Vec(shapes) => {
            for shape in shapes {
                collect_rendered_text_from_shape(shape, out);
            }
        }
        _ => {}
    }
}

fn render_routine_assistant_contents_texts(app: &mut GENtleApp) -> Vec<String> {
    let ctx = egui::Context::default();
    ctx.begin_pass(egui::RawInput::default());
    crate::egui_compat::show_central_panel_for_test_context(
        &ctx,
        egui::CentralPanel::default(),
        |ui| {
            app.render_routine_assistant_contents(ui);
        },
    );
    let full_output = ctx.end_pass();
    let mut texts = Vec::new();
    for clipped in full_output.shapes {
        collect_rendered_text_from_shape(&clipped.shape, &mut texts);
    }
    texts
}

#[test]
fn open_routine_assistant_dialog_focuses_existing_window_without_resetting_state() {
    let mut app = GENtleApp::default();
    app.open_routine_assistant_dialog();
    app.routine_assistant_goal = "Assemble TP53 insert with directional control".to_string();
    app.routine_assistant_status = "editing".to_string();

    app.open_routine_assistant_dialog();

    assert!(app.show_routine_assistant_dialog);
    assert_eq!(
        app.routine_assistant_goal,
        "Assemble TP53 insert with directional control"
    );
    assert_eq!(app.routine_assistant_status, "editing");
    assert!(
        app.pending_focus_viewports
            .contains(&GENtleApp::routine_assistant_viewport_id())
    );
}

#[test]
fn open_routine_assistant_dialog_starts_and_persists_decision_trace() {
    let mut app = GENtleApp::default();
    app.routine_assistant_goal = "Assemble reporter".to_string();
    app.routine_assistant_query = "golden gate".to_string();
    app.routine_assistant_candidates = vec![CloningRoutineCatalogRow {
        routine_id: "golden_gate.type_iis_single_insert".to_string(),
        title: "Golden Gate Type IIS Single Insert".to_string(),
        family: "golden_gate".to_string(),
        status: "implemented".to_string(),
        ..Default::default()
    }];

    app.open_routine_assistant_dialog();

    let trace = app
        .routine_assistant_decision_trace
        .as_ref()
        .expect("routine-assistant trace");
    assert_eq!(trace.schema, ROUTINE_DECISION_TRACE_SCHEMA);
    assert_eq!(trace.status, "draft");
    assert_eq!(trace.goal_text, "Assemble reporter");
    assert_eq!(trace.query_text, "golden gate");
    assert!(
        trace
            .candidate_routine_ids
            .iter()
            .any(|id| id == "golden_gate.type_iis_single_insert")
    );
    let store = app
        .engine
        .read()
        .unwrap()
        .state()
        .metadata
        .get(ROUTINE_DECISION_TRACES_METADATA_KEY)
        .cloned()
        .expect("trace store metadata");
    assert_eq!(
        store.get("schema").and_then(|value| value.as_str()),
        Some(ROUTINE_DECISION_TRACE_STORE_SCHEMA)
    );
    let traces = store
        .get("traces")
        .and_then(|value| value.as_array())
        .expect("trace rows");
    assert_eq!(traces.len(), 1);
}

#[test]
fn routine_assistant_trace_persists_planning_context_candidate_scores_and_macro_suggestions() {
    let mut app = GENtleApp::default();
    {
        let mut engine = app.engine.write().expect("engine write");
        engine
            .set_planning_objective(Some(PlanningObjective {
                helper_profile_id: Some("pUC19".to_string()),
                preferred_routine_families: vec!["restriction".to_string()],
                ..PlanningObjective::default()
            }))
            .expect("set planning objective");
        engine
            .apply(Operation::UpsertWorkflowMacroTemplate {
                name: "restriction_setup".to_string(),
                description: Some("Restriction digest and ligation setup".to_string()),
                details_url: Some("https://example.org/macros/restriction".to_string()),
                parameters: vec![],
                input_ports: vec![],
                output_ports: vec![],
                script: "op {\"Note\":{\"message\":\"restriction digest ligation\"}}".to_string(),
            })
            .expect("upsert workflow template");
    }
    app.routine_assistant_goal = "Subclone reporter into pUC19".to_string();
    app.routine_assistant_query = "restriction cloning".to_string();
    app.routine_assistant_selected_routine_id = "restriction.demo_subcloning".to_string();
    app.routine_assistant_candidates = vec![CloningRoutineCatalogRow {
        routine_id: "restriction.demo_subcloning".to_string(),
        title: "Restriction Digest and Ligation".to_string(),
        family: "restriction".to_string(),
        status: "implemented".to_string(),
        estimated_time_hours: Some(6.0),
        estimated_cost: Some(20.0),
        local_fit_score: Some(0.7),
        composite_meta_score: Some(0.9),
        planning_estimate: Some(PlanningEstimate {
            schema: PLANNING_ESTIMATE_SCHEMA.to_string(),
            estimated_time_hours: 6.0,
            estimated_cost: 20.0,
            local_fit_score: 0.7,
            composite_meta_score: 0.9,
            passes_guardrails: true,
            guardrail_failures: vec![],
            explanation: serde_json::json!({
                "routine_family_alignment_bonus": 0.25,
                "routine_family_alignment_sources": ["helper profile"]
            }),
        }),
        ..Default::default()
    }];

    app.open_routine_assistant_dialog();

    let trace = app
        .routine_assistant_decision_trace
        .as_ref()
        .expect("routine assistant trace");
    assert_eq!(
        trace
            .routine_preference_context
            .as_ref()
            .and_then(|row| row.helper_profile_id.as_deref()),
        Some("puc19")
    );
    assert_eq!(trace.candidate_planning_scores.len(), 1);
    assert_eq!(
        trace.candidate_planning_scores[0].routine_family_alignment_bonus,
        Some(0.25)
    );
    assert_eq!(trace.macro_suggestions.len(), 1);
    assert_eq!(
        trace.macro_suggestions[0].template_name,
        "restriction_setup"
    );

    let store = app
        .engine
        .read()
        .unwrap()
        .state()
        .metadata
        .get(ROUTINE_DECISION_TRACES_METADATA_KEY)
        .cloned()
        .expect("trace store metadata");
    let persisted = serde_json::from_value::<RoutineDecisionTraceStore>(store)
        .expect("deserialize routine trace store");
    assert_eq!(persisted.traces.len(), 1);
    assert_eq!(
        persisted.traces[0]
            .routine_preference_context
            .as_ref()
            .and_then(|row| row.helper_profile_id.as_deref()),
        Some("puc19")
    );
    assert_eq!(persisted.traces[0].candidate_planning_scores.len(), 1);
    assert_eq!(persisted.traces[0].macro_suggestions.len(), 1);
}

#[test]
fn open_routine_assistant_dialog_captures_active_sequence_variant_reasoning_context() {
    let mut app = GENtleApp::default();
    let mut dna = DNAsequence::from_sequence("ACGTACGTACGTACGT").expect("sequence");
    dna.features_mut().push(gb_io::seq::Feature {
        kind: "promoter".into(),
        location: gb_io::seq::Location::simple_range(0, 8),
        qualifiers: vec![("label".into(), Some("VKORC1 promoter".to_string()))],
    });
    dna.features_mut().push(gb_io::seq::Feature {
        kind: "TFBS".into(),
        location: gb_io::seq::Location::simple_range(4, 7),
        qualifiers: vec![("label".into(), Some("SP1 motif".to_string()))],
    });
    dna.features_mut().push(gb_io::seq::Feature {
        kind: "variation".into(),
        location: gb_io::seq::Location::simple_range(5, 6),
        qualifiers: vec![
            ("label".into(), Some("rsRoutine".to_string())),
            (
                "gentle_generated".into(),
                Some("dbsnp_variant_marker".to_string()),
            ),
        ],
    });
    dna.update_computed_features();
    app.engine
        .write()
        .unwrap()
        .state_mut()
        .sequences
        .insert("variant_routine_seq".to_string(), dna);
    app.routine_assistant_candidates = vec![CloningRoutineCatalogRow {
        routine_id: "gibson.demo_overlap".to_string(),
        title: "Gibson Demo Overlap".to_string(),
        family: "gibson".to_string(),
        status: "implemented".to_string(),
        ..Default::default()
    }];

    let viewport_id = egui::ViewportId::from_hash_of("routine_assistant_variant_seq");
    let native_key = GENtleApp::native_menu_key_for_viewport(viewport_id);
    app.windows.insert(
        viewport_id,
        Arc::new(RwLock::new(Window::new_dna_lazy(
            "variant_routine_seq".to_string(),
            app.engine.clone(),
        ))),
    );
    app.native_window_key_to_viewport
        .insert(native_key, viewport_id);
    app.active_window_menu_key = Some(native_key);

    app.open_routine_assistant_dialog();

    let context = app
        .routine_assistant_decision_trace
        .as_ref()
        .and_then(|trace| trace.routine_preference_context.as_ref())
        .expect("routine preference context");
    assert_eq!(
        context.construct_reasoning_seq_id.as_deref(),
        Some("variant_routine_seq")
    );
    assert!(
        context
            .variant_effect_tags
            .iter()
            .any(|tag| tag == "promoter_tfbs_regulatory_candidate")
    );
    assert!(
        context
            .variant_suggested_assay_ids
            .iter()
            .any(|assay| assay == "allele_paired_promoter_luciferase_reporter")
    );
    assert!(
        context
            .variant_derived_preferred_routine_families
            .iter()
            .any(|family| family == "gibson")
    );
}

#[test]
fn routine_assistant_compare_stage_renders_sequence_aware_planning_summary() {
    let mut app = GENtleApp::default();
    app.routine_assistant_stage = RoutineAssistantStage::Compare;
    app.routine_assistant_selected_routine_id = "gibson.demo_overlap".to_string();
    app.routine_assistant_candidates = vec![CloningRoutineCatalogRow {
        routine_id: "gibson.demo_overlap".to_string(),
        title: "Gibson Demo Overlap".to_string(),
        family: "gibson".to_string(),
        status: "implemented".to_string(),
        summary: Some("Assembly-compatible demo routine.".to_string()),
        ..Default::default()
    }];
    app.routine_assistant_explain_output = Some(serde_json::json!({
        "planning": {
            "seq_id": "variant_demo_seq",
            "estimate": {
                "composite_meta_score": 0.93,
                "local_fit_score": 0.81,
                "estimated_time_hours": 4.5,
                "estimated_cost": 12.0,
                "explanation": {
                    "routine_family_alignment_bonus": 0.16,
                    "routine_family_alignment_sources": ["variant_derived"]
                }
            }
        },
        "alternatives": []
    }));

    let texts = render_routine_assistant_contents_texts(&mut app);
    assert!(
        texts.iter().any(|text| text.contains(
            "sequence-aware planning: score 0.930 | fit 0.810 | time 4.50 h | cost 12.00"
        )),
        "rendered texts: {texts:?}"
    );
    assert!(
        texts
            .iter()
            .any(|text| text.contains("alignment bonus: +0.16 (variant_derived)")),
        "rendered texts: {texts:?}"
    );
}

#[test]
fn routine_assistant_decision_trace_normalizes_disambiguation_and_preflight_history() {
    let mut app = GENtleApp::default();
    app.ensure_routine_assistant_decision_trace_started();

    let explain_output = serde_json::json!({
        "explanation": {
            "disambiguation_questions": [
                " Question A? ",
                "Question B?",
                "Question A?"
            ]
        }
    });
    let compare_output = serde_json::json!({
        "comparison": {
            "disambiguation_questions": [
                "Question B?",
                "Question C?"
            ]
        }
    });
    let explain_questions =
        GENtleApp::routine_assistant_disambiguation_questions_from_output(&explain_output);
    let compare_questions =
        GENtleApp::routine_assistant_disambiguation_questions_from_output(&compare_output);

    app.update_routine_assistant_decision_trace(|trace| {
        trace.disambiguation_questions_presented = explain_questions;
        GENtleApp::merge_routine_assistant_disambiguation_questions(
            &mut trace.disambiguation_questions_presented,
            compare_questions,
        );
        trace.disambiguation_answers = vec![
            RoutineDecisionTraceDisambiguationAnswer {
                question_id: "question_c".to_string(),
                answer_text: " answer C ".to_string(),
            },
            RoutineDecisionTraceDisambiguationAnswer {
                question_id: "question_a".to_string(),
                answer_text: " answer A ".to_string(),
            },
            RoutineDecisionTraceDisambiguationAnswer {
                question_id: "question_a".to_string(),
                answer_text: " revised answer A ".to_string(),
            },
        ];
        GENtleApp::routine_assistant_commit_preflight_snapshot(
            trace,
            Some(RoutineDecisionTracePreflightSnapshot {
                can_execute: true,
                warnings: vec![" warning_a ".to_string(), "warning_a".to_string()],
                errors: vec![],
                contract_source: Some(" routine_catalog ".to_string()),
            }),
        );
        GENtleApp::routine_assistant_commit_preflight_snapshot(
            trace,
            Some(RoutineDecisionTracePreflightSnapshot {
                can_execute: false,
                warnings: vec![],
                errors: vec![" error_b ".to_string(), "error_b".to_string()],
                contract_source: Some(" routine_catalog ".to_string()),
            }),
        );
    });

    let trace = app
        .routine_assistant_decision_trace
        .as_ref()
        .expect("active trace");
    assert_eq!(trace.disambiguation_questions_presented.len(), 3);
    assert_eq!(
        trace
            .disambiguation_questions_presented
            .iter()
            .map(|row| row.question_id.as_str())
            .collect::<Vec<_>>(),
        vec!["question_a", "question_b", "question_c"]
    );
    assert_eq!(trace.disambiguation_answers.len(), 2);
    assert_eq!(trace.disambiguation_answers[0].question_id, "question_a");
    assert_eq!(
        trace.disambiguation_answers[0].answer_text,
        "revised answer A"
    );
    assert_eq!(trace.disambiguation_answers[1].question_id, "question_c");
    assert_eq!(trace.disambiguation_answers[1].answer_text, "answer C");
    assert_eq!(trace.preflight_history.len(), 2);
    assert_eq!(trace.preflight_history[0].warnings, vec!["warning_a"]);
    assert_eq!(trace.preflight_history[1].errors, vec!["error_b"]);
    let latest_preflight = trace
        .preflight_snapshot
        .as_ref()
        .expect("latest preflight snapshot");
    assert!(!latest_preflight.can_execute);
    assert_eq!(latest_preflight.errors, vec!["error_b"]);
}

#[test]
fn routine_assistant_disambiguation_answers_snapshot_filters_and_orders() {
    let mut app = GENtleApp::default();
    let questions = vec![
        RoutineDecisionTraceDisambiguationQuestion {
            question_id: "question_b".to_string(),
            question_text: "Question B?".to_string(),
        },
        RoutineDecisionTraceDisambiguationQuestion {
            question_id: "question_a".to_string(),
            question_text: "Question A?".to_string(),
        },
    ];
    app.routine_assistant_disambiguation_answers
        .insert("question_b".to_string(), " answer b ".to_string());
    app.routine_assistant_disambiguation_answers
        .insert("question_a".to_string(), " answer a ".to_string());
    app.routine_assistant_disambiguation_answers
        .insert("stale_question".to_string(), "stale".to_string());

    app.sync_routine_assistant_disambiguation_answers_for_questions(&questions, &[]);
    let answers = app.routine_assistant_disambiguation_answers_snapshot(&questions);
    assert_eq!(answers.len(), 2);
    assert_eq!(answers[0].question_id, "question_a");
    assert_eq!(answers[0].answer_text, "answer a");
    assert_eq!(answers[1].question_id, "question_b");
    assert_eq!(answers[1].answer_text, "answer b");
    assert!(
        !app.routine_assistant_disambiguation_answers
            .contains_key("stale_question")
    );
}

#[test]
fn maybe_mark_routine_assistant_trace_aborted_marks_non_terminal_trace() {
    let mut app = GENtleApp::default();
    app.ensure_routine_assistant_decision_trace_started();
    app.update_routine_assistant_decision_trace(|trace| {
        trace.status = "draft".to_string();
    });

    app.maybe_mark_routine_assistant_trace_aborted();

    let trace = app
        .routine_assistant_decision_trace
        .as_ref()
        .expect("active trace");
    assert_eq!(trace.status, "aborted");
    let store = app
        .engine
        .read()
        .unwrap()
        .state()
        .metadata
        .get(ROUTINE_DECISION_TRACES_METADATA_KEY)
        .cloned()
        .expect("trace store metadata");
    let traces = store
        .get("traces")
        .and_then(|value| value.as_array())
        .expect("trace rows");
    assert_eq!(
        traces[0].get("status").and_then(|value| value.as_str()),
        Some("aborted")
    );
}

#[test]
fn maybe_mark_routine_assistant_trace_aborted_preserves_terminal_statuses() {
    for terminal in ["executed", "execution_failed", "aborted", "exported"] {
        let mut app = GENtleApp::default();
        app.ensure_routine_assistant_decision_trace_started();
        app.update_routine_assistant_decision_trace(|trace| {
            trace.status = terminal.to_string();
        });

        app.maybe_mark_routine_assistant_trace_aborted();

        let trace = app
            .routine_assistant_decision_trace
            .as_ref()
            .expect("active trace");
        assert_eq!(trace.status, terminal);
    }
}

#[test]
fn open_planning_dialog_focuses_existing_window_without_resetting_state() {
    let mut app = GENtleApp::default();
    app.open_planning_dialog();
    app.planning_sync_source = "local_lab".to_string();
    app.planning_host_filter = "hsdR".to_string();
    app.planning_host_selected_id = "ecoli_k12_restriction_positive".to_string();
    app.planning_helper_filter = "factor xa".to_string();
    app.planning_helper_selected_id = "pGEX".to_string();
    app.planning_status = "editing".to_string();

    app.open_planning_dialog();

    assert!(app.show_planning_dialog);
    assert_eq!(app.planning_sync_source, "local_lab");
    assert_eq!(app.planning_host_filter, "hsdR");
    assert_eq!(
        app.planning_host_selected_id,
        "ecoli_k12_restriction_positive"
    );
    assert_eq!(app.planning_helper_filter, "factor xa");
    assert_eq!(app.planning_helper_selected_id, "pGEX");
    assert_eq!(app.planning_status, "editing");
    assert!(
        app.pending_focus_viewports
            .contains(&GENtleApp::planning_viewport_id())
    );
}

#[test]
fn open_gibson_dialog_focuses_existing_window_without_resetting_state() {
    let mut app = GENtleApp::default();
    app.open_gibson_dialog();
    app.gibson_destination_seq_id = "destination_vector".to_string();
    app.gibson_status = "editing".to_string();

    app.open_gibson_dialog();

    assert!(app.show_gibson_dialog);
    assert_eq!(app.gibson_destination_seq_id, "destination_vector");
    assert_eq!(app.gibson_status, "editing");
    assert!(
        app.pending_focus_viewports
            .contains(&GENtleApp::gibson_viewport_id())
    );
}

#[test]
fn cancel_gibson_dialog_closes_window_without_resetting_draft() {
    let mut app = GENtleApp::default();
    app.show_gibson_dialog = true;
    app.gibson_destination_seq_id = "destination_vector".to_string();
    app.gibson_insert_seq_id = "insert_x".to_string();
    app.gibson_opening_start_0based = "941".to_string();
    app.gibson_opening_end_0based_exclusive = "941".to_string();

    app.cancel_gibson_dialog();

    assert!(!app.show_gibson_dialog);
    assert_eq!(app.gibson_destination_seq_id, "destination_vector");
    assert_eq!(app.gibson_insert_seq_id, "insert_x");
    assert_eq!(app.gibson_opening_start_0based, "941");
    assert_eq!(app.gibson_opening_end_0based_exclusive, "941");
}

#[test]
fn open_gibson_show_all_unique_cutters_enables_global_scan_mode() {
    let mut app = GENtleApp::default();

    app.open_gibson_show_all_unique_cutters();

    assert!(app.gibson_show_all_unique_cutters);
    assert!(app.gibson_status.contains("unique cleavage sites"));
}

#[test]
fn open_gibson_dialog_prefills_active_sequence() {
    let mut app = GENtleApp::default();
    app.engine.write().unwrap().state_mut().sequences.insert(
        "active_seq".to_string(),
        DNAsequence::from_sequence("ACGTACGTACGT").expect("active sequence"),
    );
    let viewport_id = egui::ViewportId::from_hash_of("gibson_prefill_seq");
    let native_key = GENtleApp::native_menu_key_for_viewport(viewport_id);
    app.windows.insert(
        viewport_id,
        Arc::new(RwLock::new(Window::new_dna_lazy(
            "active_seq".to_string(),
            app.engine.clone(),
        ))),
    );
    app.native_window_key_to_viewport
        .insert(native_key, viewport_id);
    app.active_window_menu_key = Some(native_key);

    app.open_gibson_dialog();

    assert_eq!(app.gibson_destination_seq_id, "active_seq");
}

#[test]
fn gibson_preview_readiness_requires_defined_cut_edges() {
    let mut app = GENtleApp::default();
    app.gibson_destination_seq_id = "destination_vector".to_string();
    app.gibson_insert_seq_id = "insert_x".to_string();
    app.gibson_opening_mode = GibsonUiOpeningMode::DefinedSite;

    let err = app
        .validate_gibson_preview_readiness()
        .expect_err("missing cut edges should block preview");
    assert!(err.contains("choose a cutter suggestion") || err.contains("cut edges"));

    app.gibson_opening_start_0based = "941".to_string();
    app.gibson_opening_end_0based_exclusive = "941".to_string();

    app.validate_gibson_preview_readiness()
        .expect("filled cut edges should allow preview");
}

#[test]
fn build_gibson_plan_from_ui_uses_relaxed_default_max_anneal_hits() {
    let mut app = GENtleApp::default();
    app.engine.write().unwrap().state_mut().sequences.insert(
        "destination_vector".to_string(),
        DNAsequence::from_sequence("ACGTACGTACGTACGTACGT").expect("destination sequence"),
    );
    app.engine.write().unwrap().state_mut().sequences.insert(
        "insert_x".to_string(),
        DNAsequence::from_sequence("ATGCGTACGTTAGCGTACGA").expect("insert sequence"),
    );
    app.gibson_destination_seq_id = "destination_vector".to_string();
    app.gibson_insert_seq_id = "insert_x".to_string();
    app.gibson_opening_mode = GibsonUiOpeningMode::DefinedSite;
    app.gibson_opening_start_0based = "4".to_string();
    app.gibson_opening_end_0based_exclusive = "4".to_string();
    app.gibson_overlap_bp_min = "20".to_string();
    app.gibson_overlap_bp_max = "40".to_string();
    app.gibson_minimum_overlap_tm_celsius = "60.0".to_string();
    app.gibson_priming_tm_min_celsius = "58.0".to_string();
    app.gibson_priming_tm_max_celsius = "68.0".to_string();
    app.gibson_priming_length_min_bp = "18".to_string();
    app.gibson_priming_length_max_bp = "35".to_string();

    let plan = app.build_gibson_plan_from_ui().expect("build Gibson plan");

    assert_eq!(plan.validation_policy.design_targets.max_anneal_hits, 4);
}

#[test]
fn build_gibson_plan_from_ui_supports_multiple_ordered_inserts() {
    let mut app = GENtleApp::default();
    app.engine.write().unwrap().state_mut().sequences.insert(
        "destination_vector".to_string(),
        DNAsequence::from_sequence("ACGTACGTACGTACGTACGT").expect("destination sequence"),
    );
    app.engine.write().unwrap().state_mut().sequences.insert(
        "insert_x".to_string(),
        DNAsequence::from_sequence("ATGCGTACGTTAGCGTACGA").expect("insert sequence"),
    );
    app.engine.write().unwrap().state_mut().sequences.insert(
        "insert_y".to_string(),
        DNAsequence::from_sequence("GCTAGCATCGTACGATCGTA").expect("second insert sequence"),
    );
    app.gibson_destination_seq_id = "destination_vector".to_string();
    app.gibson_insert_seq_id = "insert_x".to_string();
    app.gibson_extra_inserts.push(GibsonUiInsertRow {
        seq_id: "insert_y".to_string(),
        orientation: GibsonUiInsertOrientation::Reverse,
    });
    app.gibson_opening_mode = GibsonUiOpeningMode::DefinedSite;
    app.gibson_opening_start_0based = "4".to_string();
    app.gibson_opening_end_0based_exclusive = "4".to_string();
    app.gibson_overlap_bp_min = "20".to_string();
    app.gibson_overlap_bp_max = "40".to_string();
    app.gibson_minimum_overlap_tm_celsius = "60.0".to_string();
    app.gibson_priming_tm_min_celsius = "58.0".to_string();
    app.gibson_priming_tm_max_celsius = "68.0".to_string();
    app.gibson_priming_length_min_bp = "18".to_string();
    app.gibson_priming_length_max_bp = "35".to_string();

    let plan = app.build_gibson_plan_from_ui().expect("build Gibson plan");

    assert_eq!(plan.fragments.len(), 2);
    assert_eq!(plan.fragments[0].seq_id, "insert_x");
    assert_eq!(plan.fragments[1].seq_id, "insert_y");
    assert_eq!(plan.fragments[1].orientation, "reverse");
    assert_eq!(plan.assembly_order.len(), 4);
    assert_eq!(plan.junctions.len(), 3);
    assert_eq!(plan.junctions[1].id, "junction_1_2");
}

#[test]
fn build_gibson_plan_from_ui_preserves_requested_unique_site() {
    let mut app = GENtleApp::default();
    app.engine.write().unwrap().state_mut().sequences.insert(
        "destination_vector".to_string(),
        DNAsequence::from_sequence("ACGTACGTACGTACGTACGT").expect("destination sequence"),
    );
    app.engine.write().unwrap().state_mut().sequences.insert(
        "insert_x".to_string(),
        DNAsequence::from_sequence("ATGCGTACGTTAGCGTACGA").expect("insert sequence"),
    );
    app.gibson_destination_seq_id = "destination_vector".to_string();
    app.gibson_insert_seq_id = "insert_x".to_string();
    app.gibson_opening_mode = GibsonUiOpeningMode::DefinedSite;
    app.gibson_opening_start_0based = "4".to_string();
    app.gibson_opening_end_0based_exclusive = "4".to_string();
    app.gibson_unique_restriction_site_enzyme_name = "EcoRI".to_string();

    let plan = app.build_gibson_plan_from_ui().expect("build Gibson plan");

    assert_eq!(
        plan.validation_policy
            .desired_unique_restriction_site_enzyme_name
            .as_deref(),
        Some("EcoRI")
    );
}

#[test]
fn load_gibson_plan_into_ui_restores_requested_unique_site() {
    let mut app = GENtleApp::default();
    let mut plan = GibsonAssemblyPlan::default();
    plan.validation_policy
        .desired_unique_restriction_site_enzyme_name = Some("EcoRI".to_string());

    app.load_gibson_plan_into_ui(&plan);

    assert_eq!(app.gibson_unique_restriction_site_enzyme_name, "EcoRI");
}

#[test]
fn gibson_target_review_rows_use_structured_design_adjustments() {
    let app = GENtleApp::default();
    let preview = GibsonAssemblyPreview {
        fragments: vec![GibsonPreviewInsert {
            fragment_id: "insert_1".to_string(),
            ..Default::default()
        }],
        resolved_junctions: vec![
            GibsonResolvedJunctionPreview::default(),
            GibsonResolvedJunctionPreview::default(),
        ],
        errors: vec![
            "Could not derive a Gibson priming segment for the left insert end".to_string(),
        ],
        suggested_design_adjustments: vec![
            GibsonSuggestedDesignAdjustment {
                target: GibsonDesignAdjustmentTarget::PrimingSegmentMaxLengthBp,
                label: "Set max priming length to 40 bp".to_string(),
                summary:
                    "Increase the maximum 3' priming length from 35 bp to 40 bp and rerun preview."
                        .to_string(),
                rationale: "Unused sequence remains beyond the current cap.".to_string(),
                current_value: 35.0,
                suggested_value: 40.0,
            },
            GibsonSuggestedDesignAdjustment {
                target: GibsonDesignAdjustmentTarget::PrimingSegmentTmMinCelsius,
                label: "Set min priming Tm to 52.4 °C".to_string(),
                summary:
                    "Lower the minimum 3' priming Tm from 58.0 °C to 52.4 °C and rerun preview."
                        .to_string(),
                rationale:
                    "The strongest blocked 3' priming segment stays below the current minimum."
                        .to_string(),
                current_value: 58.0,
                suggested_value: 52.4,
            },
        ],
        ..Default::default()
    };

    let rows = app.gibson_target_review_rows(&preview);
    let joined = rows
        .iter()
        .map(|(_, row)| row.as_str())
        .collect::<Vec<_>>()
        .join("\n");

    assert!(joined.contains("set max priming length to 40 bp"));
    assert!(joined.contains("set min priming Tm to 52.4 °C"));
}

#[test]
fn apply_gibson_design_adjustment_updates_target_field() {
    let mut app = GENtleApp::default();
    app.apply_gibson_design_adjustment(&GibsonSuggestedDesignAdjustment {
        target: GibsonDesignAdjustmentTarget::PrimingSegmentMaxLengthBp,
        label: "Set max priming length to 40 bp".to_string(),
        summary: String::new(),
        rationale: String::new(),
        current_value: 35.0,
        suggested_value: 40.0,
    });
    assert_eq!(app.gibson_priming_length_max_bp, "40");

    app.apply_gibson_design_adjustment(&GibsonSuggestedDesignAdjustment {
        target: GibsonDesignAdjustmentTarget::PrimingSegmentTmMinCelsius,
        label: "Set min priming Tm to 52.4 °C".to_string(),
        summary: String::new(),
        rationale: String::new(),
        current_value: 58.0,
        suggested_value: 52.4,
    });
    assert_eq!(app.gibson_priming_tm_min_celsius, "52.4");
}

#[test]
fn humanize_gibson_ui_text_formats_tm_and_celsius() {
    let raw = "minimum overlap Tm 60.0 C";
    assert_eq!(
        GENtleApp::humanize_gibson_ui_text(raw),
        "minimum overlap Tm 60.0 °C"
    );
}

#[test]
fn gibson_opening_detail_text_marks_cut_in_sequence() {
    let mut app = GENtleApp::default();
    app.engine.write().unwrap().state_mut().sequences.insert(
        "destination_vector".to_string(),
        DNAsequence::from_sequence("AAACCCGGGTTT").expect("destination sequence"),
    );
    app.gibson_destination_seq_id = "destination_vector".to_string();
    app.gibson_opening_start_0based = "6".to_string();
    app.gibson_opening_end_0based_exclusive = "6".to_string();

    let detail = app.gibson_opening_detail_text(&[]).expect("opening detail");

    assert!(detail.contains("destination cut: AAACCC|GGGTTT"));
}

#[test]
fn gibson_opened_destination_arm_text_shows_both_arms() {
    let mut app = GENtleApp::default();
    app.engine.write().unwrap().state_mut().sequences.insert(
        "destination_vector".to_string(),
        DNAsequence::from_sequence("AAACCCGGGTTT").expect("destination sequence"),
    );
    app.gibson_destination_seq_id = "destination_vector".to_string();
    app.gibson_opening_start_0based = "6".to_string();
    app.gibson_opening_end_0based_exclusive = "6".to_string();

    let opened = app
        .gibson_opened_destination_arm_text()
        .expect("opened arm detail");

    assert!(opened.contains("left arm  : ...AAACCC|"));
    assert!(opened.contains("right arm : |GGGTTT..."));
}

#[test]
fn gibson_preview_junctions_text_is_copyable_plaintext() {
    let preview = GibsonAssemblyPreview {
        resolved_junctions: vec![GibsonResolvedJunctionPreview {
            junction_id: "junction_left".to_string(),
            left_member_id: "dest_left".to_string(),
            right_member_id: "insert_1".to_string(),
            overlap_bp: 20,
            left_member_bp: 20,
            right_member_bp: 0,
            overlap_tm_celsius: 61.5,
            overlap_sequence: "AAATCGGATCTGATCGAAGG".to_string(),
            overlap_source: "derive_from_destination_left_flank".to_string(),
            distinct_from: vec![],
        }],
        ..Default::default()
    };

    let text = GENtleApp::gibson_preview_junctions_text(&preview);

    assert!(text.contains("junction_left"));
    assert!(text.contains("overlap: 20 bp | AAATCGGATCTGATCGAAGG"));
    assert!(text.contains("Tm: 61.5 °C"));
}

#[test]
fn gibson_preview_primers_text_is_copyable_plaintext() {
    let preview = GibsonAssemblyPreview {
        primer_suggestions: vec![GibsonPrimerSuggestion {
            primer_id: "insert_1_left_insert_primer".to_string(),
            side: "left_insert_primer".to_string(),
            fragment_id: "insert_1".to_string(),
            template_seq_id: "insert_demo".to_string(),
            full_sequence: "AAATCGGATCTGATCGAAGGATGGCTACCGTTAAGCTG".to_string(),
            overlap_5prime: GibsonPrimerSegment {
                sequence: "AAATCGGATCTGATCGAAGG".to_string(),
                length_bp: 20,
                tm_celsius: 61.5,
                gc_fraction: 0.5,
                anneal_hits: 1,
                note: String::new(),
            },
            priming_3prime: GibsonPrimerSegment {
                sequence: "ATGGCTACCGTTAAGCTG".to_string(),
                length_bp: 18,
                tm_celsius: 59.2,
                gc_fraction: 0.5,
                anneal_hits: 1,
                note: String::new(),
            },
        }],
        ..Default::default()
    };

    let text = GENtleApp::gibson_preview_primers_text(&preview);

    assert!(text.contains("insert_1_left_insert_primer"));
    assert!(text.contains("5' overlap: AAATCGGATCTGATCGAAGG (20 bp, 61.5 °C)"));
    assert!(text.contains("3' priming: ATGGCTACCGTTAAGCTG (18 bp, 59.2 °C, hits=1)"));
}

#[test]
fn gibson_destination_suggestion_help_text_marks_mcs_linked_cutter() {
    let suggestion = GibsonDestinationOpeningSuggestion {
        kind: "unique_restriction_site".to_string(),
        label: "SmaI".to_string(),
        summary: "MCS-linked unique cutter".to_string(),
        feature_context: "MCS 934..949 | bla".to_string(),
        start_0based: 941,
        end_0based_exclusive: 941,
        enzyme_name: Some("SmaI".to_string()),
        recognition_sequence: "CCCGGG".to_string(),
        recognition_start_0based: Some(938),
        recognition_end_0based_exclusive: Some(944),
        rebase_cut_summary: "CCCGGG | 3|3".to_string(),
        end_geometry: "blunt".to_string(),
        overhang_bp: None,
        in_mcs_context: true,
    };

    let help = GENtleApp::gibson_destination_suggestion_help_text(&suggestion);
    assert!(help.contains("named in the selected MCS annotation"));
    assert_eq!(
        GENtleApp::gibson_destination_suggestion_display_label(&suggestion),
        "SmaI (MCS)"
    );
    assert_eq!(
        GENtleApp::gibson_destination_suggestion_feature_label(&suggestion),
        "MCS 934..949 | bla"
    );
}

#[test]
fn gibson_destination_suggestion_cut_label_stays_compact() {
    let suggestion = GibsonDestinationOpeningSuggestion {
        kind: "unique_restriction_site".to_string(),
        label: "SmaI".to_string(),
        summary: "Blunt cutpoint".to_string(),
        feature_context: "MCS 934..949 | bla".to_string(),
        start_0based: 941,
        end_0based_exclusive: 941,
        enzyme_name: Some("SmaI".to_string()),
        recognition_sequence: "CCCGGG".to_string(),
        recognition_start_0based: Some(938),
        recognition_end_0based_exclusive: Some(944),
        rebase_cut_summary: "CCCGGG | 3|3".to_string(),
        end_geometry: "blunt".to_string(),
        overhang_bp: None,
        in_mcs_context: true,
    };

    assert_eq!(
        GENtleApp::gibson_destination_suggestion_cut_label(&suggestion),
        "CCCGGG | 3|3"
    );
}

#[test]
fn displayed_gibson_destination_opening_suggestions_hide_non_mcs_rows_by_default() {
    let mcs_row = GibsonDestinationOpeningSuggestion {
        kind: "unique_restriction_site".to_string(),
        label: "SmaI".to_string(),
        summary: "MCS-linked".to_string(),
        feature_context: "MCS 934..949 | bla".to_string(),
        start_0based: 941,
        end_0based_exclusive: 941,
        enzyme_name: Some("SmaI".to_string()),
        recognition_sequence: "CCCGGG".to_string(),
        recognition_start_0based: Some(938),
        recognition_end_0based_exclusive: Some(944),
        rebase_cut_summary: "CCCGGG | 3|3".to_string(),
        end_geometry: "blunt".to_string(),
        overhang_bp: None,
        in_mcs_context: true,
    };
    let other_row = GibsonDestinationOpeningSuggestion {
        in_mcs_context: false,
        label: "PstI".to_string(),
        feature_context: "lacIq".to_string(),
        enzyme_name: Some("PstI".to_string()),
        recognition_sequence: "CTGCAG".to_string(),
        recognition_start_0based: Some(120),
        recognition_end_0based_exclusive: Some(126),
        rebase_cut_summary: "CTGCAG | 5|1".to_string(),
        end_geometry: "3prime_overhang".to_string(),
        overhang_bp: Some(4),
        start_0based: 121,
        end_0based_exclusive: 125,
        summary: "Other unique cutter".to_string(),
        kind: "unique_restriction_site".to_string(),
    };
    let suggestions = vec![mcs_row.clone(), other_row];

    let (displayed, hidden_count) =
        GENtleApp::displayed_gibson_destination_opening_suggestions(&suggestions, false);

    assert_eq!(hidden_count, 1);
    assert_eq!(displayed, vec![mcs_row]);
}

#[test]
fn gibson_apply_blocked_status_surfaces_existing_termini_limitation() {
    let preview = GibsonAssemblyPreview {
        schema: "gentle.gibson_assembly_preview.v1".to_string(),
        plan_id: "gibson_multi_insert".to_string(),
        title: "Gibson preview".to_string(),
        summary: "multi-insert preview".to_string(),
        can_execute: false,
        destination: GibsonPreviewDestination::default(),
        fragments: vec![],
        insert: GibsonPreviewInsert::default(),
        resolved_junctions: vec![],
        primer_suggestions: vec![],
        warnings: vec![],
        errors: vec![GENtleApp::gibson_multi_insert_defined_opening_note().to_string()],
        notes: vec![],
        suggested_design_adjustments: vec![],
        unique_restriction_site: None,
        cartoon: GibsonCartoonPreview::default(),
        routine_handoff: GibsonRoutineHandoffPreview::default(),
    };

    assert_eq!(
        GENtleApp::gibson_apply_blocked_status(&preview),
        format!(
            "Gibson apply blocked: {}",
            GENtleApp::gibson_multi_insert_defined_opening_note()
        )
    );
}

#[test]
fn reopen_gibson_specialist_from_operation_loads_plan_and_preview() {
    let mut app = GENtleApp::default();
    app.engine
        .write()
        .unwrap()
        .state_mut()
        .sequences
        .insert("destination_vector".to_string(), {
            let mut dna =
                DNAsequence::from_sequence("AAACCCGGGTTTAAACCCGGGTTTAAACCCGGGTTTAAACCCGGGTTT")
                    .expect("destination sequence");
            dna.set_circular(true);
            dna
        });
    app.engine.write().unwrap().state_mut().sequences.insert(
        "insert_x_amplicon".to_string(),
        DNAsequence::from_sequence("ATGCGTACGTTAGCGTACGATCGTACGTAGCTAGCTAGCATCGATCGA")
            .expect("insert sequence"),
    );
    let plan_json = r#"{
  "schema": "gentle.gibson_assembly_plan.v1",
  "id": "ui_reopen_test",
  "title": "UI reopen test",
  "summary": "single insert",
  "destination": {
    "seq_id": "destination_vector",
    "topology_before_opening": "circular",
    "opening": {
      "mode": "defined_site",
      "label": "selected window",
      "start_0based": 12,
      "end_0based_exclusive": 18,
      "left_end_id": "dest_left",
      "right_end_id": "dest_right",
      "uniqueness_requirement": "must_be_unambiguous"
    }
  },
  "product": {"topology": "circular", "output_id_hint": "ui_reopen_out"},
  "fragments": [
    {
      "id": "insert_x",
      "seq_id": "insert_x_amplicon",
      "role": "insert",
      "orientation": "forward",
      "left_end_strategy": {"mode": "primer_added_overlap", "target_junction_id": "junction_left"},
      "right_end_strategy": {"mode": "primer_added_overlap", "target_junction_id": "junction_right"}
    }
  ],
  "assembly_order": [
    {"kind": "destination_end", "id": "dest_left"},
    {"kind": "fragment", "id": "insert_x"},
    {"kind": "destination_end", "id": "dest_right"}
  ],
  "junctions": [
    {
      "id": "junction_left",
      "left_member": {"kind": "destination_end", "id": "dest_left"},
      "right_member": {"kind": "fragment", "id": "insert_x"},
      "required_overlap_bp": 20,
      "overlap_partition": {"left_member_bp": 20, "right_member_bp": 0},
      "overlap_source": "derive_from_destination_left_flank",
      "distinct_from": ["junction_right"]
    },
    {
      "id": "junction_right",
      "left_member": {"kind": "fragment", "id": "insert_x"},
      "right_member": {"kind": "destination_end", "id": "dest_right"},
      "required_overlap_bp": 20,
      "overlap_partition": {"left_member_bp": 0, "right_member_bp": 20},
      "overlap_source": "derive_from_destination_right_flank",
      "distinct_from": ["junction_left"]
    }
  ],
  "validation_policy": {
    "require_unambiguous_destination_opening": true,
    "require_distinct_terminal_junctions": true,
    "adjacency_overlap_mismatch": "error",
    "design_targets": {
      "overlap_bp_min": 18,
      "overlap_bp_max": 24,
      "minimum_overlap_tm_celsius": 48.0,
      "priming_segment_tm_min_celsius": 48.0,
      "priming_segment_tm_max_celsius": 70.0,
      "priming_segment_min_length_bp": 18,
      "priming_segment_max_length_bp": 30,
      "max_anneal_hits": 4
    },
    "uniqueness_checks": {
      "destination_context": "warn",
      "participating_fragments": "warn",
      "reference_contexts": []
    }
  }
}"#;
    let result = app
        .engine
        .write()
        .unwrap()
        .apply(Operation::ApplyGibsonAssemblyPlan {
            plan_json: plan_json.to_string(),
        })
        .expect("apply Gibson operation");

    let reopened = app
        .reopen_gibson_specialist_from_operation(&result.op_id)
        .expect("reopen Gibson specialist");

    assert!(reopened);
    assert!(app.show_gibson_dialog);
    assert_eq!(app.gibson_destination_seq_id, "destination_vector");
    assert_eq!(app.gibson_insert_seq_id, "insert_x_amplicon");
    assert_eq!(app.gibson_output_id_hint, "ui_reopen_out");
    assert!(app.gibson_preview_output.is_some());
}

#[test]
fn summarize_operation_apply_gibson_includes_key_context() {
    let summary = GENtleApp::summarize_operation(&Operation::ApplyGibsonAssemblyPlan {
        plan_json: r#"{
  "schema": "gentle.gibson_assembly_plan.v1",
  "destination": {"seq_id": "dest", "opening": {"mode": "defined_site"}},
  "product": {"output_id_hint": "assembled"},
  "fragments": [{"seq_id": "insert"}]
}"#
        .to_string(),
    });
    assert!(summary.contains("Gibson cloning"));
    assert!(summary.contains("destination=dest"));
    assert!(summary.contains("inserts=insert"));
    assert!(summary.contains("output=assembled"));
}

#[test]
fn refresh_lineage_cache_renders_gibson_outputs_as_individual_sequences() {
    let mut app = GENtleApp::default();
    let plan_json = r#"{
  "schema": "gentle.gibson_assembly_plan.v1",
  "id": "lineage_test",
  "title": "Lineage test",
  "summary": "single insert",
  "destination": {
    "seq_id": "destination_vector",
    "topology_before_opening": "circular",
    "opening": {
      "mode": "defined_site",
      "label": "selected window",
      "start_0based": 12,
      "end_0based_exclusive": 18,
      "left_end_id": "dest_left",
      "right_end_id": "dest_right",
      "uniqueness_requirement": "must_be_unambiguous"
    }
  },
  "product": {"topology": "circular", "output_id_hint": "lineage_out"},
  "fragments": [
    {
      "id": "insert_x",
      "seq_id": "insert_x_amplicon",
      "role": "insert",
      "orientation": "forward",
      "left_end_strategy": {"mode": "primer_added_overlap", "target_junction_id": "junction_left"},
      "right_end_strategy": {"mode": "primer_added_overlap", "target_junction_id": "junction_right"}
    }
  ],
  "assembly_order": [
    {"kind": "destination_end", "id": "dest_left"},
    {"kind": "fragment", "id": "insert_x"},
    {"kind": "destination_end", "id": "dest_right"}
  ],
  "junctions": [
    {
      "id": "junction_left",
      "left_member": {"kind": "destination_end", "id": "dest_left"},
      "right_member": {"kind": "fragment", "id": "insert_x"},
      "required_overlap_bp": 20,
      "overlap_partition": {"left_member_bp": 20, "right_member_bp": 0},
      "overlap_source": "derive_from_destination_left_flank",
      "distinct_from": ["junction_right"]
    },
    {
      "id": "junction_right",
      "left_member": {"kind": "fragment", "id": "insert_x"},
      "right_member": {"kind": "destination_end", "id": "dest_right"},
      "required_overlap_bp": 20,
      "overlap_partition": {"left_member_bp": 0, "right_member_bp": 20},
      "overlap_source": "derive_from_destination_right_flank",
      "distinct_from": ["junction_left"]
    }
  ],
  "validation_policy": {
    "require_unambiguous_destination_opening": true,
    "require_distinct_terminal_junctions": true,
    "adjacency_overlap_mismatch": "error",
    "design_targets": {
      "overlap_bp_min": 18,
      "overlap_bp_max": 24,
      "minimum_overlap_tm_celsius": 48.0,
      "priming_segment_tm_min_celsius": 48.0,
      "priming_segment_tm_max_celsius": 70.0,
      "priming_segment_min_length_bp": 18,
      "priming_segment_max_length_bp": 30,
      "max_anneal_hits": 4
    },
    "uniqueness_checks": {
      "destination_context": "warn",
      "participating_fragments": "warn",
      "reference_contexts": []
    }
  }
}"#;
    let op_id = {
        let mut engine = app.engine.write().unwrap();
        let mut destination =
            DNAsequence::from_sequence("AAACCCGGGTTTAAACCCGGGTTTAAACCCGGGTTTAAACCCGGGTTT")
                .expect("destination sequence");
        destination.set_name("destination_vector".to_string());
        destination.set_circular(true);
        let mut insert =
            DNAsequence::from_sequence("ATGCGTACGTTAGCGTACGATCGTACGTAGCTAGCTAGCATCGATCGA")
                .expect("insert sequence");
        insert.set_name("insert_x_amplicon".to_string());
        insert.set_circular(false);
        engine
            .state_mut()
            .sequences
            .insert("destination_vector".to_string(), destination);
        engine
            .state_mut()
            .sequences
            .insert("insert_x_amplicon".to_string(), insert);
        engine.state_mut().display = ProjectState::default().display;
        engine
            .apply(Operation::ApplyGibsonAssemblyPlan {
                plan_json: plan_json.to_string(),
            })
            .expect("apply Gibson operation")
            .op_id
    };

    app.refresh_lineage_cache_if_needed();

    for seq_id in [
        "insert_x_left_insert_primer",
        "insert_x_right_insert_primer",
        "lineage_out",
    ] {
        let row = app
            .lineage_rows
            .iter()
            .find(|row| row.seq_id == seq_id)
            .unwrap_or_else(|| panic!("missing lineage row for {seq_id}"));
        assert_eq!(row.created_by_op, op_id);
        assert_eq!(row.pool_size, 1);
        assert_eq!(row.pool_members, vec![seq_id.to_string()]);
    }
}

#[test]
fn open_reference_genome_retrieve_dialog_focuses_existing_window_without_resetting_state() {
    let mut app = GENtleApp::default();
    app.open_reference_genome_retrieve_dialog();
    app.genome_gene_filter = "TP73".to_string();
    app.genome_output_id = "tp73_roi".to_string();

    app.open_reference_genome_retrieve_dialog();

    assert!(app.show_reference_genome_retrieve_dialog);
    assert_eq!(app.genome_gene_filter, "TP73");
    assert_eq!(app.genome_output_id, "tp73_roi");
    assert!(
        app.pending_focus_viewports
            .contains(&GENtleApp::retrieve_genome_viewport_id())
    );
}

#[test]
fn genome_dialog_scope_preserves_reference_and_helper_paths_independently() {
    let mut app = GENtleApp::default();
    app.open_reference_genome_prepare_dialog();
    app.genome_catalog_path = "custom/reference_catalog.json".to_string();
    app.genome_cache_dir = "custom/reference_cache".to_string();

    app.open_helper_genome_prepare_dialog();
    assert_eq!(app.genome_dialog_scope, GenomeDialogScope::Helper);
    assert_eq!(
        app.genome_catalog_path,
        DEFAULT_HELPER_GENOME_CATALOG_PATH.to_string()
    );
    assert_eq!(
        app.genome_cache_dir,
        DEFAULT_HELPER_GENOME_CACHE_DIR.to_string()
    );

    app.genome_catalog_path = "custom/helper_catalog.json".to_string();
    app.genome_cache_dir = "custom/helper_cache".to_string();

    app.open_reference_genome_prepare_dialog();
    assert_eq!(app.genome_dialog_scope, GenomeDialogScope::Reference);
    assert_eq!(app.genome_catalog_path, "custom/reference_catalog.json");
    assert_eq!(app.genome_cache_dir, "custom/reference_cache");

    app.open_helper_genome_prepare_dialog();
    assert_eq!(app.genome_dialog_scope, GenomeDialogScope::Helper);
    assert_eq!(app.genome_catalog_path, "custom/helper_catalog.json");
    assert_eq!(app.genome_cache_dir, "custom/helper_cache");
}

#[test]
fn queue_prepared_genome_reinstall_uses_active_scope_paths() {
    let mut app = GENtleApp::default();
    app.open_helper_genome_prepare_dialog();
    app.genome_catalog_path = "custom/helper_catalog.json".to_string();
    app.genome_cache_dir = "custom/helper_cache".to_string();

    app.queue_prepared_genome_reinstall(
        "Yeast Helper",
        app.genome_dialog_scope,
        PreparedGenomeReinstallDialogHost::PrepareDialog,
    );

    let request = app
        .pending_prepared_genome_reinstall
        .as_ref()
        .expect("reinstall request");
    assert_eq!(request.genome_id, "Yeast Helper");
    assert_eq!(request.scope, GenomeDialogScope::Helper);
    assert_eq!(request.catalog_path, "custom/helper_catalog.json");
    assert_eq!(request.cache_dir, "custom/helper_cache");
    assert_eq!(
        request.dialog_host,
        PreparedGenomeReinstallDialogHost::PrepareDialog
    );
}

#[test]
fn queue_prepared_genome_removal_uses_none_cache_override_when_input_is_empty() {
    let mut app = GENtleApp::default();
    app.open_helper_genome_prepare_dialog();
    app.genome_catalog_path = "custom/helper_catalog.json".to_string();
    app.genome_cache_dir = "   ".to_string();

    app.queue_prepared_genome_removal(
        "Yeast Helper".to_string(),
        GenomeDialogScope::Helper,
        "/tmp/helper_install".to_string(),
    );

    let request = app
        .pending_prepared_genome_removal
        .as_ref()
        .expect("removal request");
    assert_eq!(request.genome_id, "Yeast Helper");
    assert_eq!(request.scope, GenomeDialogScope::Helper);
    assert_eq!(request.catalog_path, "custom/helper_catalog.json");
    assert!(request.cache_dir.is_none());
    assert_eq!(request.install_dir, "/tmp/helper_install");
}

#[test]
fn helper_catalog_entries_for_scope_apply_filter_and_include_interpretation() {
    let temp = tempdir().expect("tempdir");
    let helper_catalog_path = write_app_test_helper_catalog(temp.path());
    let mut app = GENtleApp::default();
    app.helper_genome_catalog_path = helper_catalog_path;

    let rows = app
        .helper_catalog_entries_for_scope(GenomeDialogScope::Helper, Some("factor xa"))
        .expect("helper rows");

    assert_eq!(rows.len(), 1);
    assert_eq!(rows[0].genome_id, "pGEX_like_vector");
    assert!(rows[0].interpretation.is_some());
    assert_eq!(rows[0].host_system.as_deref(), Some("Escherichia coli"));
}

#[test]
fn helper_vector_cards_and_doctor_are_available_to_gui_browser() {
    let temp = tempdir().expect("tempdir");
    let helper_catalog_path = write_app_test_helper_catalog(temp.path());
    let app = GENtleApp::default();

    let cards = app
        .helper_vector_cards_for_catalog_path(&helper_catalog_path, Some("factor xa"))
        .expect("helper vector cards");
    assert_eq!(cards.len(), 1);
    assert_eq!(cards[0].helper_id, "pGEX_like_vector");
    assert_eq!(
        cards[0].sequence_availability.as_deref(),
        Some("source-backed synthetic test fixture")
    );
    assert!(!cards[0].metadata_only_candidate);

    let issues = app
        .helper_vector_doctor_issues_for_catalog_path(&helper_catalog_path)
        .expect("helper vector doctor");
    assert!(issues.is_empty());
}

#[test]
fn selected_helper_catalog_entry_resolves_aliases_to_interpretation_entry() {
    let temp = tempdir().expect("tempdir");
    let helper_catalog_path = write_app_test_helper_catalog(temp.path());
    let mut app = GENtleApp::default();
    app.open_helper_genome_prepare_dialog();
    app.genome_catalog_path = helper_catalog_path.clone();
    app.helper_genome_catalog_path = helper_catalog_path;
    app.genome_id = "pGEX".to_string();

    let entry = app
        .selected_helper_catalog_entry_for_scope(GenomeDialogScope::Helper)
        .expect("selected helper entry")
        .expect("helper entry");

    assert_eq!(entry.genome_id, "pGEX_like_vector");
    let interpretation = entry.interpretation.expect("helper interpretation");
    assert_eq!(
        interpretation.offered_functions,
        vec![
            "affinity_purification".to_string(),
            "bacterial_expression".to_string(),
            "fusion_tagging".to_string(),
            "insert_cloning".to_string(),
            "protease_cleavage".to_string(),
            "protease_tag_removal".to_string(),
        ]
    );
}

#[test]
fn selected_helper_catalog_entry_uses_active_dialog_catalog_path() {
    let temp = tempdir().expect("tempdir");
    let helper_catalog_path = write_app_test_helper_catalog(temp.path());
    let mut app = GENtleApp::default();
    app.open_helper_genome_prepare_dialog();
    app.helper_genome_catalog_path = "/tmp/old_helper_catalog.json".to_string();
    app.genome_catalog_path = helper_catalog_path;
    app.genome_id = "pGEX".to_string();

    let entry = app
        .selected_helper_catalog_entry_for_scope(GenomeDialogScope::Helper)
        .expect("selected helper entry")
        .expect("helper entry");

    assert_eq!(entry.genome_id, "pGEX_like_vector");
}

#[test]
fn helper_construct_interpretation_lines_include_core_reasoning_fields() {
    let interpretation = HelperConstructInterpretation {
        helper_id: "pGEX_like_vector".to_string(),
        aliases: vec!["pGEX".to_string()],
        helper_kinds: vec!["plasmid_vector".to_string()],
        host_systems: vec!["Escherichia coli".to_string()],
        offered_functions: vec![
            "affinity_purification".to_string(),
            "fusion_tagging".to_string(),
        ],
        constraints: vec!["reading_frame_must_be_preserved".to_string()],
        procurement_channels: vec!["vendor_catalog".to_string()],
        normalized_terms: vec![crate::genomes::HelperConstructNormalizedTerm {
            axis: "component_kind".to_string(),
            value: "fusion_tag".to_string(),
            label: Some("fusion tag".to_string()),
            source: "component:gst_tag".to_string(),
            vocabulary_label: Some("Fusion tag".to_string()),
            vocabulary_description: Some("Frame-sensitive tag".to_string()),
            vocabulary_source: Some("built-in vocabulary".to_string()),
        }],
        routine_hints: vec![crate::genomes::HelperConstructRoutineHint {
            family: "gibson".to_string(),
            rationale: "Overlap assembly keeps the fusion-tag CDS in frame.".to_string(),
            source_terms: vec!["component_kind:fusion_tag".to_string()],
        }],
        local_variant_unpublished: true,
        ..Default::default()
    };

    let lines = GENtleApp::format_helper_construct_interpretation_detail_lines(&interpretation);

    assert!(lines.contains(&"helper id: pGEX_like_vector".to_string()));
    assert!(lines.contains(&"aliases: pGEX".to_string()));
    assert!(lines.contains(&"helper kind: plasmid_vector".to_string()));
    assert!(lines.contains(&"host system: Escherichia coli".to_string()));
    assert!(lines.contains(&"constraints: reading_frame_must_be_preserved".to_string()));
    assert!(lines.contains(&"procurement channels: vendor_catalog".to_string()));
    assert!(lines.contains(&"routine hints: gibson".to_string()));
    assert!(
        lines
            .iter()
            .any(|line| line.contains("normalized terms: component_kind=fusion_tag"))
    );
    assert!(lines.contains(&"local variant: unpublished".to_string()));
}

#[test]
fn prepared_genome_reinstall_confirm_dialog_routes_to_requested_host() {
    let mut app = GENtleApp::default();
    app.queue_prepared_genome_reinstall(
        "GRCh38",
        GenomeDialogScope::Reference,
        PreparedGenomeReinstallDialogHost::PrepareDialog,
    );

    assert!(app.pending_prepared_genome_reinstall_targets_host(
        PreparedGenomeReinstallDialogHost::PrepareDialog
    ));
    assert!(
        !app.pending_prepared_genome_reinstall_targets_host(
            PreparedGenomeReinstallDialogHost::Root
        )
    );
}

#[test]
fn dismiss_pending_prepared_genome_reinstall_only_clears_matching_host() {
    let mut app = GENtleApp::default();
    app.queue_prepared_genome_reinstall(
        "GRCh38",
        GenomeDialogScope::Reference,
        PreparedGenomeReinstallDialogHost::PrepareDialog,
    );

    app.dismiss_pending_prepared_genome_reinstall_for_host(PreparedGenomeReinstallDialogHost::Root);
    assert!(app.pending_prepared_genome_reinstall.is_some());

    app.dismiss_pending_prepared_genome_reinstall_for_host(
        PreparedGenomeReinstallDialogHost::PrepareDialog,
    );
    assert!(app.pending_prepared_genome_reinstall.is_none());
}

#[test]
fn apply_prepared_genome_reinstall_request_updates_scope_and_clears_cached_gene_state() {
    let mut app = GENtleApp::default();
    app.open_reference_genome_prepare_dialog();
    app.genome_catalog_path = "custom/reference_catalog.json".to_string();
    app.genome_cache_dir = "custom/reference_cache".to_string();
    app.sync_active_genome_scope_paths_from_fields();
    app.genome_id = "Reference A".to_string();
    app.genome_genes_error = "stale genes".to_string();
    app.genome_genes_loaded_key = Some("cached".to_string());
    app.genome_selected_gene = Some(3);
    app.genome_gene_filter_page = 7;
    app.genome_biotype_filter
        .insert("protein_coding".to_string(), true);
    app.genome_retrieve_contig_suggestions = vec!["chr17".to_string()];

    app.apply_prepared_genome_reinstall_request(PreparedGenomeReinstallRequest {
        genome_id: "Helper B".to_string(),
        scope: GenomeDialogScope::Helper,
        catalog_path: "custom/helper_catalog.json".to_string(),
        cache_dir: "custom/helper_cache".to_string(),
        dialog_host: PreparedGenomeReinstallDialogHost::Root,
    });

    assert_eq!(app.genome_dialog_scope, GenomeDialogScope::Helper);
    assert_eq!(app.genome_catalog_path, "custom/helper_catalog.json");
    assert_eq!(app.genome_cache_dir, "custom/helper_cache");
    assert_eq!(app.helper_genome_catalog_path, "custom/helper_catalog.json");
    assert_eq!(app.helper_genome_cache_dir, "custom/helper_cache");
    assert_eq!(
        app.reference_genome_catalog_path,
        "custom/reference_catalog.json"
    );
    assert_eq!(app.reference_genome_cache_dir, "custom/reference_cache");
    assert_eq!(app.genome_id, "Helper B");
    assert!(app.genome_genes_error.is_empty());
    assert!(app.genome_genes_loaded_key.is_none());
    assert!(app.genome_selected_gene.is_none());
    assert_eq!(app.genome_gene_filter_page, 0);
    assert!(app.genome_biotype_filter.is_empty());
    assert!(app.genome_retrieve_contig_suggestions.is_empty());
}

#[test]
fn open_cache_cleanup_dialog_prefills_from_active_scope() {
    let mut app = GENtleApp::default();
    app.open_helper_genome_prepare_dialog();
    app.genome_catalog_path = "custom/helper_catalog.json".to_string();
    app.genome_cache_dir = "custom/helper_cache".to_string();
    app.sync_active_genome_scope_paths_from_fields();
    app.reference_genome_cache_dir = "custom/reference_cache".to_string();

    app.open_cache_cleanup_dialog();

    assert!(app.show_cache_cleanup_dialog);
    assert_eq!(app.cache_cleanup_scope, CacheCleanupScope::Helpers);
    assert_eq!(
        app.cache_cleanup_reference_cache_dir,
        "custom/reference_cache"
    );
    assert_eq!(app.cache_cleanup_helper_cache_dir, "custom/helper_cache");
    assert_eq!(
        app.cache_cleanup_mode,
        PreparedCacheCleanupMode::DerivedIndexesOnly
    );
}

fn write_app_test_prepared_cache_install(
    root: &std::path::Path,
    dir_name: &str,
    genome_id: &str,
) -> std::path::PathBuf {
    let install_dir = root.join(dir_name);
    fs::create_dir_all(install_dir.join("blastdb")).expect("create blastdb dir");
    let sequence_path = install_dir.join("sequence.fa");
    let annotation_path = install_dir.join("annotation.gtf");
    let fasta_index_path = install_dir.join("sequence.fa.fai");
    let gene_index_path = install_dir.join("genes.json");
    let transcript_index_path = install_dir.join("transcripts.json");
    let blast_prefix = install_dir.join("blastdb").join("genome");
    fs::write(&sequence_path, ">chr1\nACGTACGT\n").expect("write fasta");
    fs::write(
        &annotation_path,
        "chr1\tsrc\tgene\t1\t8\t.\t+\t.\tgene_id \"GENE1\"; gene_name \"ONE\";\n",
    )
    .expect("write annotation");
    fs::write(&fasta_index_path, "chr1\t8\t6\t8\t9\n").expect("write fasta index");
    fs::write(&gene_index_path, "[]").expect("write gene index");
    fs::write(&transcript_index_path, "[]").expect("write transcript index");
    fs::write(blast_prefix.with_extension("nhr"), "nhr").expect("write blast nhr");
    fs::write(blast_prefix.with_extension("nin"), "nin").expect("write blast nin");
    fs::write(blast_prefix.with_extension("nsq"), "nsq").expect("write blast nsq");
    let manifest = serde_json::json!({
        "genome_id": genome_id,
        "sequence_source": sequence_path.display().to_string(),
        "annotation_source": annotation_path.display().to_string(),
        "sequence_source_type": "local",
        "annotation_source_type": "local",
        "sequence_sha1": "seqsha",
        "annotation_sha1": "annsha",
        "sequence_path": sequence_path.display().to_string(),
        "annotation_path": annotation_path.display().to_string(),
        "fasta_index_path": fasta_index_path.display().to_string(),
        "gene_index_path": gene_index_path.display().to_string(),
        "transcript_index_path": transcript_index_path.display().to_string(),
        "blast_db_prefix": blast_prefix.display().to_string(),
        "blast_index_executable": "makeblastdb",
        "blast_indexed_at_unix_ms": 123,
        "installed_at_unix_ms": 456
    });
    fs::write(
        install_dir.join("manifest.json"),
        serde_json::to_string_pretty(&manifest).expect("serialize manifest"),
    )
    .expect("write manifest");
    install_dir
}

#[test]
fn cache_cleanup_rebuild_candidates_follow_partial_cleanup_results() {
    let mut app = GENtleApp::default();
    app.reference_genome_catalog_path = "custom/reference_catalog.json".to_string();
    app.reference_genome_cache_dir = "custom/reference_cache".to_string();
    app.helper_genome_catalog_path = "custom/helper_catalog.json".to_string();
    app.helper_genome_cache_dir = "custom/helper_cache".to_string();
    app.cache_cleanup_scope = CacheCleanupScope::Both;
    app.cache_cleanup_reference_cache_dir = "custom/reference_cache".to_string();
    app.cache_cleanup_helper_cache_dir = "custom/helper_cache".to_string();
    let report = PreparedCacheCleanupReport {
        schema: "gentle.prepared_cache_cleanup.v1".to_string(),
        mode: PreparedCacheCleanupMode::DerivedIndexesOnly,
        cache_roots: vec![
            "custom/reference_cache".to_string(),
            "custom/helper_cache".to_string(),
        ],
        selected_prepared_ids: vec!["Ref A".to_string(), "Helper B".to_string()],
        selected_prepared_paths: vec![],
        include_orphaned_remnants: false,
        results: vec![
            PreparedCacheCleanupItemReport {
                entry_id: "Ref A".to_string(),
                classification: PreparedCacheEntryKind::PreparedInstall,
                cache_root: "custom/reference_cache".to_string(),
                path: "custom/reference_cache/ref_a".to_string(),
                removed: true,
                removed_artifact_groups: vec![PreparedCacheArtifactGroup::DerivedIndexes],
                removed_bytes: 12,
                removed_file_count: 2,
                skipped_reason: None,
            },
            PreparedCacheCleanupItemReport {
                entry_id: "Helper B".to_string(),
                classification: PreparedCacheEntryKind::PreparedInstall,
                cache_root: "custom/helper_cache".to_string(),
                path: "custom/helper_cache/helper_b".to_string(),
                removed: true,
                removed_artifact_groups: vec![PreparedCacheArtifactGroup::BlastDb],
                removed_bytes: 8,
                removed_file_count: 1,
                skipped_reason: None,
            },
            PreparedCacheCleanupItemReport {
                entry_id: "orphan".to_string(),
                classification: PreparedCacheEntryKind::OrphanedRemnant,
                cache_root: "custom/reference_cache".to_string(),
                path: "custom/reference_cache/orphan".to_string(),
                removed: true,
                removed_artifact_groups: vec![],
                removed_bytes: 4,
                removed_file_count: 1,
                skipped_reason: None,
            },
        ],
        entry_count: 3,
        removed_item_count: 3,
        removed_bytes: 24,
        removed_file_count: 4,
    };

    let candidates = app.cache_cleanup_rebuild_candidates_from_report(&report);

    assert_eq!(candidates.len(), 2);
    assert_eq!(candidates[0].genome_id, "Helper B");
    assert_eq!(candidates[0].scope, GenomeDialogScope::Helper);
    assert_eq!(candidates[0].catalog_path, "custom/helper_catalog.json");
    assert_eq!(candidates[0].cache_dir, "custom/helper_cache");
    assert_eq!(
        candidates[0].dialog_host,
        PreparedGenomeReinstallDialogHost::Root
    );
    assert_eq!(candidates[1].genome_id, "Ref A");
    assert_eq!(candidates[1].scope, GenomeDialogScope::Reference);
    assert_eq!(candidates[1].catalog_path, "custom/reference_catalog.json");
    assert_eq!(candidates[1].cache_dir, "custom/reference_cache");
}

#[test]
fn queue_cache_cleanup_rebuild_candidate_sets_pending_reinstall_request() {
    let mut app = GENtleApp::default();
    app.cache_cleanup_rebuild_candidates = vec![PreparedGenomeReinstallRequest {
        genome_id: "Ref A".to_string(),
        scope: GenomeDialogScope::Reference,
        catalog_path: "custom/reference_catalog.json".to_string(),
        cache_dir: "custom/reference_cache".to_string(),
        dialog_host: PreparedGenomeReinstallDialogHost::Root,
    }];

    app.queue_cache_cleanup_rebuild_candidate(0);

    let request = app
        .pending_prepared_genome_reinstall
        .as_ref()
        .expect("pending rebuild request");
    assert_eq!(request.genome_id, "Ref A");
    assert_eq!(request.scope, GenomeDialogScope::Reference);
    assert_eq!(request.catalog_path, "custom/reference_catalog.json");
    assert_eq!(request.cache_dir, "custom/reference_cache");
    assert_eq!(request.dialog_host, PreparedGenomeReinstallDialogHost::Root);
}

#[test]
fn cache_cleanup_apply_scopes_duplicate_entry_ids_per_root() {
    let temp = tempdir().expect("tempdir");
    let reference_root = temp.path().join("reference_cache");
    let helper_root = temp.path().join("helper_cache");
    let reference_install =
        write_app_test_prepared_cache_install(&reference_root, "shared_ref", "SharedToy");
    let helper_install =
        write_app_test_prepared_cache_install(&helper_root, "shared_helper", "SharedToy");

    let mut app = GENtleApp::default();
    app.cache_cleanup_scope = CacheCleanupScope::Both;
    app.cache_cleanup_reference_cache_dir = reference_root.display().to_string();
    app.cache_cleanup_helper_cache_dir = helper_root.display().to_string();
    app.cache_cleanup_mode = PreparedCacheCleanupMode::DerivedIndexesOnly;
    app.refresh_cache_cleanup_inspection();

    let report = app
        .cache_cleanup_inspection
        .as_ref()
        .expect("inspection report");
    assert_eq!(report.entry_count, 2);
    let helper_entry = report
        .entries
        .iter()
        .find(|entry| {
            std::path::Path::new(&entry.path)
                .file_name()
                .and_then(|value| value.to_str())
                == Some("shared_helper")
        })
        .expect("helper cache row");
    app.cache_cleanup_selected_paths
        .insert(helper_entry.path.clone());

    let (targets, _, _, _) = app.cache_cleanup_preview_targets();
    assert_eq!(targets.len(), 1);
    assert_eq!(targets[0].entry_id, "SharedToy");
    assert_eq!(
        std::path::Path::new(&targets[0].path)
            .file_name()
            .and_then(|value| value.to_str()),
        Some("shared_helper")
    );

    app.apply_cache_cleanup();

    assert!(helper_install.join("sequence.fa").exists());
    assert!(!helper_install.join("sequence.fa.fai").exists());
    assert!(!helper_install.join("genes.json").exists());
    assert!(!helper_install.join("blastdb").join("genome.nsq").exists());

    assert!(reference_install.join("sequence.fa").exists());
    assert!(reference_install.join("sequence.fa.fai").exists());
    assert!(reference_install.join("genes.json").exists());
    assert!(
        reference_install
            .join("blastdb")
            .join("genome.nsq")
            .exists()
    );

    assert!(app.cache_cleanup_status.contains("removed 1 item"));
}

#[test]
fn apply_cache_cleanup_refreshes_gui_state_after_partial_cleanup() {
    let _env_guard = EnvVarGuard::set(
        crate::genomes::MAKEBLASTDB_ENV_BIN,
        "__gentle_makeblastdb_missing_for_test__",
    );
    let td = tempdir().expect("tempdir");
    let (catalog_path, cache_dir) = write_toy_prepare_catalog(td.path());
    let catalog =
        crate::genomes::GenomeCatalog::from_json_file(&catalog_path).expect("load toy catalog");
    catalog
        .prepare_genome_once("ToyGenome")
        .expect("prepare genome");
    let install_dir = std::path::Path::new(&cache_dir).join("toygenome");
    assert!(install_dir.join("sequence.fa.fai").exists());
    assert!(install_dir.join("genes.json").exists());

    let mut app = GENtleApp::default();
    app.genome_dialog_scope = GenomeDialogScope::Reference;
    app.genome_catalog_path = catalog_path.clone();
    app.genome_cache_dir = cache_dir.clone();
    app.reference_genome_catalog_path = catalog_path.clone();
    app.reference_genome_cache_dir = cache_dir.clone();

    app.open_cache_cleanup_dialog();
    assert_eq!(
        app.cache_cleanup_inspection
            .as_ref()
            .expect("initial inspection")
            .entry_count,
        1
    );

    let selected_path = install_dir.display().to_string();
    app.cache_cleanup_selected_paths.insert(selected_path);
    app.cache_cleanup_confirm_pending = true;
    app.cache_cleanup_rebuild_candidates = vec![PreparedGenomeReinstallRequest {
        genome_id: "stale".to_string(),
        scope: GenomeDialogScope::Helper,
        catalog_path: "stale/catalog.json".to_string(),
        cache_dir: "stale/cache".to_string(),
        dialog_host: PreparedGenomeReinstallDialogHost::Root,
    }];

    app.apply_cache_cleanup();

    assert!(!app.cache_cleanup_confirm_pending);
    assert!(app.cache_cleanup_selected_paths.is_empty());
    assert!(
        app.cache_cleanup_status
            .starts_with("Cache cleanup removed 1 item(s), reclaiming ")
    );
    assert!(install_dir.join("sequence.fa").exists());
    assert!(install_dir.join("annotation.gtf").exists());
    assert!(!install_dir.join("sequence.fa.fai").exists());
    assert!(!install_dir.join("genes.json").exists());
    assert_eq!(app.cache_cleanup_rebuild_candidates.len(), 1);
    assert_eq!(
        app.cache_cleanup_rebuild_candidates[0].genome_id,
        "ToyGenome"
    );
    assert_eq!(
        app.cache_cleanup_rebuild_candidates[0].scope,
        GenomeDialogScope::Reference
    );
    assert_eq!(
        app.cache_cleanup_rebuild_candidates[0].catalog_path,
        catalog_path
    );
    assert!(GENtleApp::cache_cleanup_root_matches(
        &app.cache_cleanup_rebuild_candidates[0].cache_dir,
        &cache_dir
    ));
    let refreshed = app
        .cache_cleanup_inspection
        .as_ref()
        .expect("refreshed inspection");
    assert_eq!(refreshed.entry_count, 1);
    let entry = refreshed
        .entries
        .iter()
        .find(|entry| entry.entry_id == "ToyGenome")
        .expect("prepared entry remains after partial cleanup");
    assert!(
        entry
            .artifact_stats
            .iter()
            .all(|stat| stat.group == PreparedCacheArtifactGroup::CachedSources)
    );
}

#[test]
fn cache_cleanup_preview_targets_distinguish_duplicate_ids_by_path() {
    let mut app = GENtleApp::default();
    app.cache_cleanup_mode = PreparedCacheCleanupMode::DerivedIndexesOnly;
    let left_path = "/tmp/cache_a/toygenome".to_string();
    let right_path = "/tmp/cache_b/toygenome".to_string();
    app.cache_cleanup_selected_paths.insert(left_path.clone());
    app.cache_cleanup_inspection = Some(PreparedCacheInspectionReport {
        schema: "gentle.prepared_cache_inspection.v1".to_string(),
        cache_roots: vec!["/tmp/cache_a".to_string(), "/tmp/cache_b".to_string()],
        entries: vec![
            PreparedCacheInspectionEntry {
                entry_id: "ToyGenome".to_string(),
                classification: PreparedCacheEntryKind::PreparedInstall,
                cache_root: "/tmp/cache_a".to_string(),
                path: left_path.clone(),
                artifact_stats: vec![PreparedCacheArtifactStat {
                    group: PreparedCacheArtifactGroup::DerivedIndexes,
                    total_size_bytes: 7,
                    file_count: 1,
                }],
                total_size_bytes: 7,
                file_count: 1,
            },
            PreparedCacheInspectionEntry {
                entry_id: "ToyGenome".to_string(),
                classification: PreparedCacheEntryKind::PreparedInstall,
                cache_root: "/tmp/cache_b".to_string(),
                path: right_path,
                artifact_stats: vec![PreparedCacheArtifactStat {
                    group: PreparedCacheArtifactGroup::DerivedIndexes,
                    total_size_bytes: 11,
                    file_count: 2,
                }],
                total_size_bytes: 11,
                file_count: 2,
            },
        ],
        entry_count: 2,
        total_size_bytes: 18,
        total_file_count: 3,
    });

    let (targets, bytes, root_count, _) = app.cache_cleanup_preview_targets();

    assert_eq!(targets.len(), 1);
    assert_eq!(targets[0].path, left_path);
    assert_eq!(bytes, 7);
    assert_eq!(root_count, 2);
}

#[test]
fn prepare_dialog_primary_action_allows_selecting_installed_genomes_for_reinstall() {
    let all_genomes = vec!["Alpha".to_string(), "Beta".to_string()];
    let preparable_set = HashSet::from(["Alpha".to_string()]);

    assert_eq!(
        GENtleApp::prepare_dialog_primary_action("Alpha", &all_genomes, &preparable_set),
        PrepareGenomeDialogPrimaryAction::Prepare
    );
    assert_eq!(
        GENtleApp::prepare_dialog_primary_action("Beta", &all_genomes, &preparable_set),
        PrepareGenomeDialogPrimaryAction::Reindex
    );
    assert_eq!(
        GENtleApp::prepare_dialog_primary_action("Gamma", &all_genomes, &preparable_set),
        PrepareGenomeDialogPrimaryAction::None
    );
}

#[test]
fn open_window_entries_use_helper_titles_for_helper_scope_dialogs() {
    let mut app = GENtleApp::default();
    app.open_helper_genome_blast_dialog();

    let entries = app.collect_open_window_entries();

    assert!(entries.iter().any(|entry| {
        entry.viewport_id == GENtleApp::blast_genome_viewport_id()
            && entry.title == "BLAST Helper Sequence"
    }));
}

#[test]
fn open_window_entries_include_retrieve_genome_viewport_when_open() {
    let mut app = GENtleApp::default();
    app.show_reference_genome_retrieve_dialog = true;

    let entries = app.collect_open_window_entries();

    assert!(entries.iter().any(|entry| {
        entry.viewport_id == GENtleApp::retrieve_genome_viewport_id()
            && entry.title == "Retrieve Genomic Sequence"
    }));
}

#[test]
fn open_window_entries_include_background_jobs_and_jaspar_when_open() {
    let mut app = GENtleApp::default();
    app.show_jobs_panel = true;
    app.show_jaspar_expert_dialog = true;

    let entries = app.collect_open_window_entries();

    assert!(entries.iter().any(|entry| {
        entry.viewport_id == GENtleApp::background_jobs_viewport_id()
            && entry.title == "Background Jobs"
    }));
    assert!(entries.iter().any(|entry| {
        entry.viewport_id == GENtleApp::jaspar_expert_viewport_id()
            && entry.title == "JASPAR Expert"
    }));
}

#[test]
fn open_configuration_graphics_dialog_focuses_existing_window_without_resetting_edits() {
    let mut app = GENtleApp::default();
    app.open_configuration_dialog();
    app.configuration_graphics.feature_details_font_size = 17.75;
    app.configuration_graphics_dirty = true;
    app.configuration_status = "editing".to_string();

    app.open_configuration_graphics_dialog();

    assert!(app.show_configuration_dialog);
    assert!(matches!(app.configuration_tab, ConfigurationTab::Graphics));
    assert!(app.configuration_graphics_dirty);
    assert_eq!(app.configuration_graphics.feature_details_font_size, 17.75);
    assert_eq!(app.configuration_status, "editing");
    assert!(
        app.pending_focus_viewports
            .contains(&GENtleApp::configuration_viewport_id())
    );
}

#[test]
fn open_help_doc_focuses_help_viewport() {
    let mut app = GENtleApp::default();
    app.open_help_doc(HelpDoc::Gui);
    assert!(app.show_help_dialog);
    assert!(
        app.pending_focus_viewports
            .contains(&GENtleApp::help_viewport_id())
    );
}

#[test]
fn open_help_doc_preserves_loaded_markdown_without_forced_reload() {
    let mut app = GENtleApp::default();
    app.help_gui_markdown = "custom help text".to_string();

    app.open_help_doc(HelpDoc::Gui);

    assert_eq!(app.help_gui_markdown, "custom help text");
}

#[test]
fn rendered_help_markdown_cache_reuses_entry_until_source_changes() {
    let mut app = GENtleApp::default();
    app.help_gui_markdown = "Use `Cmd+K`.".to_string();

    let first = app.rendered_help_markdown_for(HelpDoc::Gui);
    let second = app.rendered_help_markdown_for(HelpDoc::Gui);
    assert_eq!(
        first.as_ref(),
        GENtleApp::help_display_markdown(app.help_gui_markdown.as_str())
    );
    assert!(std::sync::Arc::ptr_eq(&first, &second));

    app.help_gui_markdown = "Use `Cmd+Shift+K`.".to_string();
    let third = app.rendered_help_markdown_for(HelpDoc::Gui);
    assert_eq!(
        third.as_ref(),
        GENtleApp::help_display_markdown(app.help_gui_markdown.as_str())
    );
    assert!(!std::sync::Arc::ptr_eq(&first, &third));
}

#[test]
fn open_help_doc_when_same_tab_open_only_queues_focus() {
    let mut app = GENtleApp::default();
    app.show_help_dialog = true;
    app.help_doc = HelpDoc::Gui;
    app.help_search_query = "persist".to_string();
    app.help_search_matches = vec![HelpSearchMatch {
        line_number: 7,
        snippet: "keep existing matches".to_string(),
    }];

    app.open_help_doc(HelpDoc::Gui);

    assert_eq!(app.help_search_matches.len(), 1);
    assert_eq!(app.help_search_matches[0].line_number, 7);
    assert!(
        app.pending_focus_viewports
            .contains(&GENtleApp::help_viewport_id())
    );
}

#[test]
fn help_content_width_requires_relayout_for_meaningful_resize() {
    assert!(GENtleApp::help_content_width_requires_relayout(0.0, 420.0));
    assert!(GENtleApp::help_content_width_requires_relayout(
        860.0, 620.0
    ));
    assert!(!GENtleApp::help_content_width_requires_relayout(
        620.0, 621.0
    ));
    assert!(!GENtleApp::help_content_width_requires_relayout(
        620.0, 626.0
    ));
    assert!(GENtleApp::help_content_width_requires_relayout(
        620.0, 629.0
    ));
    assert!(!GENtleApp::help_content_width_requires_relayout(
        620.0,
        f32::NAN
    ));
}

#[test]
fn help_topic_combo_width_is_clamped() {
    assert_eq!(GENtleApp::clamp_help_topic_combo_width(200.0), 180.0);
    assert_eq!(GENtleApp::clamp_help_topic_combo_width(1_200.0), 420.0);
    assert_eq!(GENtleApp::clamp_help_topic_combo_width(f32::NAN), 280.0);
}

#[test]
fn configuration_backdrop_path_field_width_is_clamped() {
    assert_eq!(
        GENtleApp::clamp_configuration_backdrop_path_field_width(120.0),
        180.0
    );
    assert_eq!(
        GENtleApp::clamp_configuration_backdrop_path_field_width(1_200.0),
        260.0
    );
    assert_eq!(
        GENtleApp::clamp_configuration_backdrop_path_field_width(f32::NAN),
        240.0
    );
}

#[test]
fn embedded_window_open_state_respects_internal_close_actions() {
    assert!(GENtleApp::reconcile_embedded_window_open_state(true, true));
    assert!(!GENtleApp::reconcile_embedded_window_open_state(
        false, true
    ));
    assert!(!GENtleApp::reconcile_embedded_window_open_state(
        true, false
    ));
    assert!(!GENtleApp::reconcile_embedded_window_open_state(
        false, false
    ));
}

#[test]
fn help_markdown_max_image_width_tracks_available_width() {
    assert_eq!(GENtleApp::help_markdown_max_image_width(200.0), 220);
    assert_eq!(GENtleApp::help_markdown_max_image_width(400.0), 300);
    assert_eq!(GENtleApp::help_markdown_max_image_width(10_000.0), 1600);
}

#[test]
fn rewrite_markdown_inline_code_soft_breaks_preserves_fenced_blocks() {
    let markdown = "inline `File -> Open Tutorial Project...`\n\n```bash\ncargo run --bin gentle_cli -- gibson preview\n```\n";
    let rewritten = GENtleApp::rewrite_markdown_inline_code_soft_breaks(markdown);
    assert!(rewritten.contains("File"));
    assert!(rewritten.contains('\u{200B}'));
    assert!(rewritten.contains("```bash\ncargo run --bin gentle_cli -- gibson preview\n```"));
}

#[test]
fn help_display_markdown_adds_soft_breaks_for_tutorial_style_inline_code() {
    let markdown = "Use `Patterns -> Gibson...` and `docs/examples/workflows/gibson_specialist_testing_baseline.json`.";
    let rendered = GENtleApp::help_display_markdown(markdown);
    assert!(rendered.contains("Patterns"));
    assert!(rendered.contains("docs/"));
    assert!(rendered.contains('\u{200B}'));
}

#[test]
fn help_display_markdown_summarizes_generated_tutorial_front_matter() {
    let markdown = "---\nchapter_id: \"load_branch_reverse_complement_pgex_fasta\"\nsource_example: \"docs/examples/workflows/load_branch_reverse_complement_pgex_fasta.json\"\n---\n\n# Load FASTA\n\nBody line.\n";
    let rendered = GENtleApp::help_display_markdown(markdown);
    let normalized = rendered.replace('\u{200B}', "");

    assert!(!normalized.starts_with("---"));
    assert!(normalized.starts_with("# Load FASTA\n\n_Provenance note:"));
    assert!(normalized.contains("chapter `load_branch_reverse_complement_pgex_fasta`"));
    assert!(normalized.contains(
        "workflow `docs/examples/workflows/load_branch_reverse_complement_pgex_fasta.json`"
    ));
    assert!(normalized.contains("the hands-on walkthrough starts here"));
    assert!(normalized.contains("\n\nBody line."));
}

#[test]
fn open_help_tutorial_doc_switches_to_tutorial_view_and_loads_markdown() {
    let temp = tempdir().expect("tempdir");
    let tutorial_path = temp.path().join("tutorial.md");
    fs::write(&tutorial_path, "# Tutorial Test\n\nBody line.\n").expect("write tutorial");

    let mut app = GENtleApp::default();
    app.help_tutorial_entries = vec![HelpTutorialDocEntry {
        title: "Tutorial Test".to_string(),
        path: tutorial_path.to_string_lossy().to_string(),
        summary: "docs/tutorial/tutorial.md".to_string(),
        audiences: vec![],
        group_label: None,
        group_order: None,
        group_position: None,
        decimal_id: None,
        review_status: None,
        codex_reviewed_at: None,
        human_reviewed_at: None,
        human_reviewer: None,
        review_stale: false,
    }];
    app.help_tutorial_selected = 0;
    app.open_help_tutorial_doc(0);

    assert_eq!(app.help_doc, HelpDoc::Tutorial);
    assert!(app.show_help_dialog);
    assert_eq!(app.help_tutorial_title, "Tutorial Test");
    assert!(app.help_tutorial_markdown.contains("Body line."));
    assert!(
        app.pending_focus_viewports
            .contains(&GENtleApp::help_viewport_id())
    );
}

#[test]
fn tutorial_feedback_context_includes_id_and_source_path() {
    let entries = GENtleApp::discover_help_tutorial_entries();
    let entry = entries
        .iter()
        .find(|entry| entry.summary.contains("simple_pcr_selection_gui"))
        .or_else(|| {
            entries
                .iter()
                .find(|entry| entry.path.contains("simple_pcr_selection_gui"))
        })
        .expect("simple PCR tutorial entry should be discoverable");

    let context = GENtleApp::tutorial_feedback_context_text_for_entry(
        entry,
        Some("line 12: primer"),
        "test-version",
        "test-platform",
    )
    .expect("feedback context");

    assert!(context.contains("Tutorial id: simple_pcr_selection_gui"));
    assert!(context.contains("docs/tutorial/sources/"));
    assert!(context.contains("simple_pcr_selection_gui"));
    assert!(context.contains("Current search section: line 12: primer"));
    assert!(context.contains("GENtle version: test-version"));
}

#[test]
fn tutorial_audience_group_uses_catalog_audiences() {
    let entry = HelpTutorialDocEntry {
        title: "qPCR".to_string(),
        path: "docs/tutorial/04-03_qpcr_exon_junctions_gui.md".to_string(),
        summary: String::new(),
        audiences: vec!["primer_design".to_string()],
        group_label: None,
        group_order: None,
        group_position: None,
        decimal_id: None,
        review_status: None,
        codex_reviewed_at: None,
        human_reviewed_at: None,
        human_reviewer: None,
        review_stale: false,
    };

    assert_eq!(
        GENtleApp::tutorial_audience_group_label(&entry),
        "PCR, qPCR, And Direct Sequence Inspection"
    );
}

#[test]
fn open_help_tutorial_path_adds_dynamic_entry_and_loads_markdown() {
    let temp = tempdir().expect("tempdir");
    let tutorial_path = temp.path().join("tutorial.md");
    fs::write(&tutorial_path, "# Dynamic Tutorial\n\nWalkthrough.\n").expect("write tutorial");

    let mut app = GENtleApp::default();
    app.open_help_tutorial_path(
        tutorial_path.to_string_lossy().as_ref(),
        "Fallback Tutorial",
        "dynamic guide",
    )
    .expect("open direct tutorial path");

    assert_eq!(app.help_doc, HelpDoc::Tutorial);
    assert!(app.show_help_dialog);
    assert_eq!(app.help_tutorial_title, "Dynamic Tutorial");
    assert!(app.help_tutorial_markdown.contains("Walkthrough."));
    assert!(
        app.help_tutorial_entries
            .iter()
            .any(|entry| entry.path == tutorial_path.to_string_lossy())
    );
}

#[test]
fn help_copyable_text_uses_active_document_markdown() {
    let mut app = GENtleApp::default();
    app.help_doc = HelpDoc::Tutorial;
    app.help_tutorial_markdown = "# Tutorial\n\n- `gibson_ui_test_product`\n".to_string();
    assert_eq!(
        app.active_help_copyable_text(),
        "# Tutorial\n\n- `gibson_ui_test_product`\n"
    );

    app.help_doc = HelpDoc::Gui;
    app.help_gui_markdown = "GUI body".to_string();
    assert_eq!(app.active_help_copyable_text(), "GUI body");
}

#[test]
fn command_palette_includes_routine_assistant_entry() {
    let app = GENtleApp::default();
    let entries = app.collect_command_palette_entries();
    assert!(
        entries
            .iter()
            .any(|entry| entry.title == "Routine Assistant")
    );
}

#[test]
fn command_palette_includes_planning_entry() {
    let app = GENtleApp::default();
    let entries = app.collect_command_palette_entries();
    assert!(entries.iter().any(|entry| entry.title == "Planning"));
}

#[test]
fn command_palette_includes_gibson_entry() {
    let app = GENtleApp::default();
    let entries = app.collect_command_palette_entries();
    assert!(entries.iter().any(|entry| entry.title == "Gibson"));
}

#[test]
fn command_palette_includes_lab_assistant_report_export() {
    let app = GENtleApp::default();
    let entries = app.collect_command_palette_entries();
    assert!(
        entries
            .iter()
            .any(|entry| entry.title == "Export Lab Assistant Report"
                && matches!(entry.action, CommandPaletteAction::ExportLabAssistantReport))
    );
}

#[test]
fn command_palette_includes_external_services_entry() {
    let app = GENtleApp::default();
    let entries = app.collect_command_palette_entries();
    assert!(
        entries
            .iter()
            .any(|entry| entry.title == "External Services")
    );
}

#[test]
fn command_palette_includes_new_sequence_entries() {
    let app = GENtleApp::default();
    let entries = app.collect_command_palette_entries();

    assert!(entries.iter().any(|entry| {
        entry.title == "New Sequence" && matches!(entry.action, CommandPaletteAction::NewSequence)
    }));
    assert!(entries.iter().any(|entry| {
        entry.title == "New Sequence from Clipboard"
            && matches!(entry.action, CommandPaletteAction::NewSequenceFromClipboard)
    }));
}

#[test]
fn command_palette_includes_pcr_designer_entry() {
    let app = GENtleApp::default();
    let entries = app.collect_command_palette_entries();
    assert!(entries.iter().any(|entry| entry.title == "PCR Designer"));
}

#[test]
fn command_palette_includes_sequencing_confirmation_entry() {
    let app = GENtleApp::default();
    let entries = app.collect_command_palette_entries();
    assert!(
        entries
            .iter()
            .any(|entry| entry.title == "Sequencing Confirmation")
    );
}

#[test]
fn command_palette_includes_mirna_target_scan_entry() {
    let app = GENtleApp::default();
    let entries = app.collect_command_palette_entries();
    assert!(
        entries
            .iter()
            .any(|entry| entry.title == "microRNA Target Scan")
    );
}

#[test]
fn command_palette_includes_shared_ui_intent_entries() {
    let app = GENtleApp::default();
    let entries = app.collect_command_palette_entries();

    for target in UiIntentTarget::all() {
        assert!(
                entries.iter().any(|entry| {
                    entry.title == target.discoverability_title()
                        && entry.detail == target.discoverability_detail()
                        && entry.keywords.contains(target.discoverability_keywords())
                        && entry.keywords.contains(target.as_str())
                        && matches!(entry.action, CommandPaletteAction::UiIntent(entry_target) if entry_target == *target)
                }),
                "missing command palette entry for {}",
                target.as_str()
            );
    }
}

#[test]
fn command_palette_includes_gui_prominent_glossary_entries() {
    let app = GENtleApp::default();
    let entries = app.collect_command_palette_entries();

    for gui_entry in gui_prominent_glossary_entries() {
        assert!(
            entries
                .iter()
                .any(|entry| entry.title == gui_entry.palette_title),
            "missing command palette entry `{}` for GUI-prominent glossary command `{}`",
            gui_entry.palette_title,
            gui_entry.glossary_path
        );
        if let Some(target) = gui_entry.ui_intent_target {
            assert_eq!(target.discoverability_title(), gui_entry.palette_title);
            assert!(
                gui_entry.menu_path.starts_with(target.menu_path()),
                "GUI-prominent glossary command `{}` should keep its menu path aligned with UI intent target `{}`",
                gui_entry.glossary_path,
                target.as_str()
            );
        }
    }
}

#[test]
fn command_palette_includes_evidence_preparation_direct_action() {
    let app = GENtleApp::default();
    let entries = app.collect_command_palette_entries();
    assert!(entries.iter().any(|entry| {
        entry.title == "Evidence Preparation"
            && matches!(entry.action, CommandPaletteAction::OpenEvidencePreparation)
    }));
}

#[test]
fn evidence_preparation_is_not_a_gui_prominent_glossary_entry() {
    assert!(
        !gui_prominent_glossary_entries()
            .iter()
            .any(|entry| entry.palette_title == "Evidence Preparation")
    );
}

#[test]
fn execute_command_palette_action_opens_routine_assistant_dialog() {
    let mut app = GENtleApp::default();
    app.routine_assistant_candidates
        .push(CloningRoutineCatalogRow::default());

    app.execute_command_palette_action(
        &egui::Context::default(),
        CommandPaletteAction::OpenRoutineAssistant,
    );

    assert!(app.show_routine_assistant_dialog);
}

#[test]
fn execute_command_palette_action_opens_new_sequence_dialog() {
    let mut app = GENtleApp::default();

    app.execute_command_palette_action(
        &egui::Context::default(),
        CommandPaletteAction::NewSequence,
    );

    assert!(app.show_new_sequence_dialog);
    assert!(
        app.pending_focus_viewports
            .contains(&GENtleApp::new_sequence_viewport_id())
    );
}

#[test]
fn execute_command_palette_action_opens_planning_dialog() {
    let mut app = GENtleApp::default();

    app.execute_command_palette_action(
        &egui::Context::default(),
        CommandPaletteAction::OpenPlanning,
    );

    assert!(app.show_planning_dialog);
}

#[test]
fn execute_command_palette_action_opens_external_services_dialog() {
    let mut app = GENtleApp::default();

    app.execute_command_palette_action(
        &egui::Context::default(),
        CommandPaletteAction::OpenExternalServices,
    );

    assert!(app.external_services_ui.show_panel);
    assert!(
        app.external_services_ui
            .provider_catalog_output
            .as_ref()
            .and_then(|catalog| catalog.get("providers"))
            .and_then(serde_json::Value::as_array)
            .map(|providers| providers
                .iter()
                .any(|provider| provider["provider"].as_str() == Some("metabion")))
            .unwrap_or(false)
    );
    assert!(
        app.external_services_ui
            .request_json
            .contains("\"provider\": \"metabion\"")
    );
}

#[test]
fn external_services_request_template_uses_selected_provider_capability() {
    let mut app = GENtleApp::default();
    app.open_external_services_dialog();
    app.external_services_ui.selected_provider = "metabion".to_string();
    app.external_services_ui.selected_service_kind = "dna_fragment".to_string();

    app.reset_external_services_request_from_selection();

    assert!(
        app.external_services_ui
            .request_json
            .contains("\"provider\": \"metabion\"")
    );
    assert!(
        app.external_services_ui
            .request_json
            .contains("\"service_kind\": \"dna_fragment\"")
    );
    assert!(app.external_services_ui.request_json.contains("fragment_1"));
}

#[test]
fn external_services_provider_switch_refreshes_editable_request_template() {
    let mut app = GENtleApp::default();
    app.open_external_services_dialog();
    app.external_services_ui.selected_provider = "metabion".to_string();
    app.external_services_ui.selected_service_kind = "dna_oligo_single_tube".to_string();
    app.reset_external_services_request_from_selection();
    assert!(
        app.external_services_ui
            .request_json
            .contains("\"provider\": \"metabion\"")
    );

    app.external_services_ui.preflight_output = Some(serde_json::json!({"stale": "preflight"}));
    app.external_services_ui.quote_output = Some(serde_json::json!({"stale": "quote"}));
    app.external_services_ui.selected_provider = "geneart".to_string();

    app.reset_external_services_request_from_selection();

    assert!(
        app.external_services_ui
            .request_json
            .contains("\"provider\": \"geneart\"")
    );
    assert!(
        !app.external_services_ui
            .request_json
            .contains("\"provider\": \"metabion\"")
    );
    assert!(
        app.external_services_ui
            .request_json
            .contains("\"service_kind\": \"dna_fragment\"")
    );
    assert!(
        app.external_services_ui
            .quote_output_dir
            .contains("geneart_dna_fragment_handoff")
    );
    assert!(app.external_services_ui.preflight_output.is_none());
    assert!(app.external_services_ui.quote_output.is_none());
}

#[test]
fn external_services_preflight_uses_shared_shell_contract() {
    let mut app = GENtleApp::default();
    app.open_external_services_dialog();
    app.external_services_ui.selected_provider = "metabion".to_string();
    app.external_services_ui.selected_service_kind = "dna_oligo_single_tube".to_string();
    app.reset_external_services_request_from_selection();

    app.run_external_services_preflight();

    let preflight = app
        .external_services_ui
        .preflight_output
        .as_ref()
        .expect("preflight output");
    assert_eq!(
        preflight["schema"].as_str(),
        Some("gentle.external_service_preflight.v1")
    );
    assert_eq!(preflight["provider"].as_str(), Some("metabion"));
    assert_eq!(preflight["eligible"].as_bool(), Some(true));
}

#[test]
fn external_services_project_source_copyable_commands_parse() {
    let mut app = GENtleApp::default();
    app.external_services_ui.project_source_kind = "sequence".to_string();
    app.external_services_ui.project_source_seq_id = "seq_a".to_string();
    app.external_services_ui.project_source_range = "3..18".to_string();
    app.external_services_ui.project_source_as_construct_output = true;
    let sequence_line = app.external_services_project_source_shell_line();
    match parse_shell_line(&sequence_line).expect("parse sequence project-source command") {
        ShellCommand::ServicesRouteProjectSource {
            kind,
            seq_id,
            range,
            source_as,
            ..
        } => {
            assert_eq!(kind, "sequence");
            assert_eq!(seq_id.as_deref(), Some("seq_a"));
            assert_eq!(range.as_deref(), Some("3..18"));
            assert_eq!(source_as.as_deref(), Some("construct-output"));
        }
        other => panic!("unexpected command: {other:?}"),
    }

    app.external_services_ui.project_source_kind = "oligo-form".to_string();
    app.external_services_ui.project_source_form_id = "order_1".to_string();
    assert!(matches!(
        parse_shell_line(&app.external_services_project_source_shell_line())
            .expect("parse oligo-form project-source command"),
        ShellCommand::ServicesRouteProjectSource { kind, form_id, .. }
            if kind == "oligo-form" && form_id.as_deref() == Some("order_1")
    ));

    app.external_services_ui.project_source_kind = "primer-report-rows".to_string();
    app.external_services_ui.project_source_report_id = "report_1".to_string();
    app.external_services_ui.project_source_pair_ranks = "1,2".to_string();
    app.external_services_ui.project_source_form_id = "order_from_report".to_string();
    assert!(matches!(
        parse_shell_line(&app.external_services_project_source_shell_line())
            .expect("parse primer-report project-source command"),
        ShellCommand::ServicesRouteProjectSource { kind, report_id, pair_ranks, form_id, .. }
            if kind == "primer-report-rows"
                && report_id.as_deref() == Some("report_1")
                && pair_ranks == vec![1, 2]
                && form_id.as_deref() == Some("order_from_report")
    ));
}

#[test]
fn external_services_project_source_route_populates_editable_request_json() {
    let mut app = GENtleApp::default();
    {
        let mut engine = app.engine.write().expect("engine lock");
        engine.state_mut().sequences.insert(
            "route_seq".to_string(),
            DNAsequence::from_sequence("ACGTACGTACGTACGTACGT").expect("route sequence"),
        );
    }
    app.open_external_services_dialog();
    app.external_services_ui.project_source_kind = "sequence".to_string();
    app.external_services_ui.project_source_seq_id = "route_seq".to_string();
    app.external_services_ui.project_source_range.clear();
    app.external_services_ui.project_source_as_construct_output = false;

    app.run_external_services_project_source_route();

    let route = app
        .external_services_ui
        .route_project_source_output
        .as_ref()
        .expect("route output");
    assert_eq!(route["status"].as_str(), Some("route_ready"));
    assert_eq!(route["recommended_provider"].as_str(), Some("metabion"));

    app.use_external_services_route_candidate();

    assert!(
        app.external_services_ui
            .request_json
            .contains("\"provider\": \"metabion\"")
    );
    assert!(
        app.external_services_ui
            .request_json
            .contains("\"seq_id\": \"route_seq\"")
    );
    assert!(app.external_services_ui.preflight_output.is_none());
    assert!(app.external_services_ui.quote_output.is_none());
}

#[test]
fn external_services_quote_export_writes_bundle_files() {
    let temp = tempdir().expect("tempdir");
    let mut app = GENtleApp::default();
    app.open_external_services_dialog();
    app.external_services_ui.selected_provider = "metabion".to_string();
    app.external_services_ui.selected_service_kind = "dna_oligo_single_tube".to_string();
    app.reset_external_services_request_from_selection();
    let output_dir = temp.path().join("metabion_gui_export");
    app.external_services_ui.quote_output_dir = output_dir.to_string_lossy().to_string();

    app.export_external_services_quote_bundle();

    assert!(output_dir.join("quote_report.json").is_file());
    assert!(output_dir.join("01_handoff_markdown.md").is_file());
    assert!(
        output_dir
            .join("04_normalized_line_items_csv.csv")
            .is_file()
    );
    let quote = app
        .external_services_ui
        .quote_output
        .as_ref()
        .expect("quote output");
    assert_eq!(quote["quote_status"].as_str(), Some("handoff_ready"));
    assert!(
        quote["service_ready_bundle"]["local_files"]
            .as_array()
            .expect("local files")
            .iter()
            .any(|file| file["artifact_kind"].as_str() == Some("quote_report_json"))
    );
    assert!(
        app.external_services_ui
            .status
            .contains("Exported quote handoff bundle")
    );
}

#[test]
fn execute_command_palette_action_opens_gibson_dialog() {
    let mut app = GENtleApp::default();

    app.execute_command_palette_action(&egui::Context::default(), CommandPaletteAction::OpenGibson);

    assert!(app.show_gibson_dialog);
}

#[test]
fn execute_command_palette_action_opens_mirna_target_scan_dialog() {
    let mut app = GENtleApp::default();

    app.execute_command_palette_action(
        &egui::Context::default(),
        CommandPaletteAction::OpenMirnaTargetScan,
    );

    assert!(app.mirna_panel.show_panel);
}

#[test]
fn execute_command_palette_action_opens_evidence_preparation_dialog() {
    let mut app = GENtleApp::default();

    app.execute_command_palette_action(
        &egui::Context::default(),
        CommandPaletteAction::OpenEvidencePreparation,
    );

    assert!(app.evidence_preparation_panel.show_panel);
}

#[test]
fn execute_command_palette_action_opens_pcr_design_dialog() {
    let mut state = ProjectState::default();
    state.sequences.insert(
        "seq1".to_string(),
        DNAsequence::from_sequence("ACGTACGT").expect("sequence"),
    );
    let mut app = GENtleApp::default();
    app.engine = Arc::new(RwLock::new(GentleEngine::from_state(state)));

    app.execute_command_palette_action(
        &egui::Context::default(),
        CommandPaletteAction::UiIntent(UiIntentTarget::PcrDesign),
    );

    assert!(app.show_pcr_design_dialog);
    assert_eq!(app.pcr_design_seq_id, "seq1");
}

#[test]
fn execute_command_palette_action_opens_sequencing_confirmation_dialog() {
    let mut state = ProjectState::default();
    state.sequences.insert(
        "seq1".to_string(),
        DNAsequence::from_sequence("ACGTACGT").expect("sequence"),
    );
    let mut app = GENtleApp::default();
    app.engine = Arc::new(RwLock::new(GentleEngine::from_state(state)));

    app.execute_command_palette_action(
        &egui::Context::default(),
        CommandPaletteAction::UiIntent(UiIntentTarget::SequencingConfirmation),
    );

    assert!(app.show_sequencing_confirmation_dialog);
    assert_eq!(app.sequencing_confirmation_seq_id, "seq1");
}

#[test]
fn execute_command_palette_ui_intents_open_helper_genome_dialogs() {
    let mut app = GENtleApp::default();

    app.execute_command_palette_action(
        &egui::Context::default(),
        CommandPaletteAction::UiIntent(UiIntentTarget::PrepareHelperGenome),
    );
    assert!(app.show_reference_genome_prepare_dialog);
    assert_eq!(app.genome_dialog_scope, GenomeDialogScope::Helper);

    app.show_reference_genome_prepare_dialog = false;
    app.genome_dialog_scope = GenomeDialogScope::Reference;
    app.execute_command_palette_action(
        &egui::Context::default(),
        CommandPaletteAction::UiIntent(UiIntentTarget::RetrieveHelperSequence),
    );
    assert!(app.show_reference_genome_retrieve_dialog);
    assert_eq!(app.genome_dialog_scope, GenomeDialogScope::Helper);

    app.show_reference_genome_retrieve_dialog = false;
    app.genome_dialog_scope = GenomeDialogScope::Reference;
    app.execute_command_palette_action(
        &egui::Context::default(),
        CommandPaletteAction::UiIntent(UiIntentTarget::BlastHelperSequence),
    );
    assert!(app.show_reference_genome_blast_dialog);
    assert_eq!(app.genome_dialog_scope, GenomeDialogScope::Helper);
}

#[test]
fn open_pcr_design_dialog_focuses_existing_window_without_resetting_context() {
    let mut app = GENtleApp::default();
    app.show_pcr_design_dialog = true;
    app.pcr_design_seq_id = "seq_context".to_string();

    app.open_pcr_design_dialog();

    assert!(app.show_pcr_design_dialog);
    assert_eq!(app.pcr_design_seq_id, "seq_context");
    assert!(
        app.pending_focus_viewports
            .contains(&GENtleApp::pcr_design_viewport_id())
    );
}

#[test]
fn open_pcr_design_dialog_for_seq_id_switches_existing_context() {
    let mut state = ProjectState::default();
    state.sequences.insert(
        "seq_old".to_string(),
        DNAsequence::from_sequence("ACGTACGT").expect("sequence"),
    );
    state.sequences.insert(
        "seq_new".to_string(),
        DNAsequence::from_sequence("TTTTCCCC").expect("sequence"),
    );
    let mut app = GENtleApp::default();
    app.engine = Arc::new(RwLock::new(GentleEngine::from_state(state)));
    app.show_pcr_design_dialog = true;
    app.pcr_design_seq_id = "seq_old".to_string();

    app.open_pcr_design_dialog_for_seq_id("seq_new")
        .expect("switch PCR Designer target");

    assert!(app.show_pcr_design_dialog);
    assert_eq!(app.pcr_design_seq_id, "seq_new");
    assert_eq!(app.new_windows.len(), 1);
    assert!(
        app.pending_focus_viewports
            .contains(&GENtleApp::pcr_design_viewport_id())
    );
}

#[test]
fn open_sequencing_confirmation_dialog_focuses_existing_window_without_resetting_context() {
    let mut app = GENtleApp::default();
    app.show_sequencing_confirmation_dialog = true;
    app.sequencing_confirmation_seq_id = "seq_context".to_string();

    app.open_sequencing_confirmation_dialog();

    assert!(app.show_sequencing_confirmation_dialog);
    assert_eq!(app.sequencing_confirmation_seq_id, "seq_context");
    assert!(
        app.pending_focus_viewports
            .contains(&GENtleApp::sequencing_confirmation_viewport_id())
    );
}

#[test]
fn open_lineage_analysis_artifact_opens_sequencing_confirmation_dialog() {
    let mut state = ProjectState::default();
    state.sequences.insert(
        "seq_confirm".to_string(),
        DNAsequence::from_sequence("ACGTACGT").expect("sequence"),
    );
    let mut app = GENtleApp::default();
    app.engine = Arc::new(RwLock::new(GentleEngine::from_state(state)));

    app.open_lineage_analysis_artifact(
        LineageAnalysisKind::SequencingConfirmation,
        "seq_confirm",
        "report_1",
    );

    assert!(app.show_sequencing_confirmation_dialog);
    assert_eq!(app.sequencing_confirmation_seq_id, "seq_confirm");
    assert_eq!(app.new_windows.len(), 1);
    assert!(
        app.pending_focus_viewports
            .contains(&GENtleApp::sequencing_confirmation_viewport_id())
    );
}

#[test]
fn open_lineage_analysis_artifact_opens_primer_design_dialog() {
    let mut state = ProjectState::default();
    state.sequences.insert(
        "seq_primer".to_string(),
        DNAsequence::from_sequence("ACGTACGT").expect("sequence"),
    );
    let mut app = GENtleApp::default();
    app.engine = Arc::new(RwLock::new(GentleEngine::from_state(state)));

    app.open_lineage_analysis_artifact(
        LineageAnalysisKind::PrimerDesign,
        "seq_primer",
        "primer_report_1",
    );

    assert!(app.show_pcr_design_dialog);
    assert_eq!(app.pcr_design_seq_id, "seq_primer");
    assert_eq!(app.new_windows.len(), 1);
    assert!(
        app.pending_focus_viewports
            .contains(&GENtleApp::pcr_design_viewport_id())
    );
}

#[test]
fn open_lineage_analysis_artifact_opens_qpcr_design_dialog() {
    let mut state = ProjectState::default();
    state.sequences.insert(
        "seq_qpcr".to_string(),
        DNAsequence::from_sequence("ACGTACGT").expect("sequence"),
    );
    let mut app = GENtleApp::default();
    app.engine = Arc::new(RwLock::new(GentleEngine::from_state(state)));

    app.open_lineage_analysis_artifact(
        LineageAnalysisKind::QpcrDesign,
        "seq_qpcr",
        "qpcr_report_1",
    );

    assert!(app.show_pcr_design_dialog);
    assert_eq!(app.pcr_design_seq_id, "seq_qpcr");
    assert_eq!(app.new_windows.len(), 1);
    assert!(
        app.pending_focus_viewports
            .contains(&GENtleApp::pcr_design_viewport_id())
    );
}

#[test]
fn open_lineage_analysis_artifact_opens_restriction_cloning_handoff_in_pcr_designer() {
    let mut state = ProjectState::default();
    state.sequences.insert(
        "seq_handoff".to_string(),
        DNAsequence::from_sequence("ACGTACGT").expect("sequence"),
    );
    let mut app = GENtleApp::default();
    app.engine = Arc::new(RwLock::new(GentleEngine::from_state(state)));

    app.open_lineage_analysis_artifact(
        LineageAnalysisKind::RestrictionCloningPcrHandoff,
        "seq_handoff",
        "restriction_handoff_1",
    );

    assert!(app.show_pcr_design_dialog);
    assert_eq!(app.pcr_design_seq_id, "seq_handoff");
    assert_eq!(app.new_windows.len(), 1);
    assert!(
        app.pending_focus_viewports
            .contains(&GENtleApp::pcr_design_viewport_id())
    );
}

#[test]
fn reopen_pcr_designer_from_operation_opens_specialist_for_template() {
    let mut state = ProjectState::default();
    state.sequences.insert(
        "seq_tpl".to_string(),
        DNAsequence::from_sequence("ACGTACGT").expect("sequence"),
    );
    let mut app = GENtleApp::default();
    app.engine = Arc::new(RwLock::new(GentleEngine::from_state(state)));
    app.lineage_reopenable_pcr_op_seq_ids
        .insert("op_pcr".to_string(), "seq_tpl".to_string());

    let reopened = app
        .reopen_pcr_designer_from_operation("op_pcr")
        .expect("reopen PCR Designer");

    assert!(reopened);
    assert!(app.show_pcr_design_dialog);
    assert_eq!(app.pcr_design_seq_id, "seq_tpl");
    assert_eq!(app.new_windows.len(), 1);
    assert!(
        app.pending_focus_viewports
            .contains(&GENtleApp::pcr_design_viewport_id())
    );
}

#[test]
fn open_lineage_analysis_artifact_opens_rna_read_mapping_workspace() {
    let mut state = ProjectState::default();
    let mut dna = DNAsequence::from_sequence("ACGTACGTACGT").expect("sequence");
    dna.features_mut().push(gb_io::seq::Feature {
        kind: "mRNA".into(),
        location: gb_io::seq::Location::simple_range(0, 12),
        qualifiers: vec![("gene".into(), Some("TP73".to_string()))],
    });
    state.sequences.insert("seq_rna".to_string(), dna);
    let mut app = GENtleApp::default();
    app.engine = Arc::new(RwLock::new(GentleEngine::from_state(state)));
    app.engine
        .write()
        .unwrap()
        .commit_rna_read_report(crate::engine::RnaReadInterpretationReport {
            schema: "gentle.rna_read_report.v1".to_string(),
            report_id: "rna_report_1".to_string(),
            seq_id: "seq_rna".to_string(),
            seed_feature_id: 0,
            generated_at_unix_ms: 1,
            profile: crate::engine::RnaReadInterpretationProfile::NanoporeCdnaV1,
            input_path: "reads.fa".to_string(),
            input_format: crate::engine::RnaReadInputFormat::Fasta,
            scope: crate::engine::SplicingScopePreset::TargetGroupTargetStrand,
            origin_mode: crate::engine::RnaReadOriginMode::SingleGene,
            read_count_total: 4,
            read_count_seed_passed: 2,
            read_count_aligned: 1,
            ..crate::engine::RnaReadInterpretationReport::default()
        })
        .expect("upsert RNA-read report");

    app.open_lineage_analysis_artifact(
        LineageAnalysisKind::RnaReadInterpretation,
        "seq_rna",
        "rna_report_1",
    );

    assert_eq!(app.new_windows.len(), 1);
}

#[test]
fn open_lineage_analysis_artifact_opens_protein_derivation_expert() {
    let mut state = ProjectState::default();
    state.sequences.insert(
        "seq_protein".to_string(),
        DNAsequence::from_sequence("ACGTACGT").expect("sequence"),
    );
    let mut app = GENtleApp::default();
    app.engine = Arc::new(RwLock::new(GentleEngine::from_state(state)));

    app.open_lineage_analysis_artifact(
        LineageAnalysisKind::ProteinDerivation,
        "seq_protein",
        "protein_report_1",
    );

    assert_eq!(app.new_windows.len(), 1);
}

#[test]
fn open_lineage_analysis_artifact_opens_reverse_translation_product_sequence() {
    let mut protein = DNAsequence::from_sequence("MKP").expect("protein");
    protein.set_name("Toy protein");
    protein.set_molecule_type("protein");
    let mut state = ProjectState::default();
    state.sequences.insert("seq_reverse".to_string(), protein);
    let mut app = GENtleApp::default();
    app.engine = Arc::new(RwLock::new(GentleEngine::from_state(state)));
    let report_id = {
        let mut engine = app.engine.write().unwrap();
        engine
            .apply(Operation::ReverseTranslateProteinSequence {
                seq_id: "seq_reverse".to_string(),
                output_id: Some("seq_reverse_coding".to_string()),
                speed_profile: Some(crate::engine::TranslationSpeedProfile::Ecoli),
                speed_mark: Some(crate::engine::TranslationSpeedMark::Slow),
                translation_table: None,
                target_anneal_tm_c: None,
                anneal_window_bp: None,
            })
            .expect("reverse translate")
            .op_id
    };

    app.open_lineage_analysis_artifact(
        LineageAnalysisKind::ReverseTranslation,
        "seq_reverse",
        &report_id,
    );

    assert_eq!(app.new_windows.len(), 1);
}

#[test]
fn open_lineage_analysis_artifact_opens_construct_reasoning_graph() {
    let mut state = ProjectState::default();
    state.sequences.insert(
        "seq_reasoning".to_string(),
        DNAsequence::from_sequence("ACGTACGT").expect("sequence"),
    );
    let mut app = GENtleApp::default();
    app.engine = Arc::new(RwLock::new(GentleEngine::from_state(state)));
    let graph_id = {
        let mut engine = app.engine.write().unwrap();
        engine
            .build_construct_reasoning_graph("seq_reasoning", None, Some("graph_reasoning_1"))
            .expect("build reasoning graph")
            .graph_id
    };

    app.open_lineage_analysis_artifact(
        LineageAnalysisKind::ConstructReasoning,
        "seq_reasoning",
        &graph_id,
    );

    assert_eq!(app.new_windows.len(), 1);
}

#[test]
fn protein_handoff_dialog_builds_reasoning_from_selected_protein_context() {
    let mut protein = DNAsequence::from_sequence("MKP").expect("protein");
    protein.set_name("Toy protein");
    protein.set_molecule_type("protein");
    let mut state = ProjectState::default();
    state.sequences.insert(
        "handoff_target".to_string(),
        DNAsequence::from_sequence("ACGTACGTACGT").expect("dna"),
    );
    state
        .sequences
        .insert("handoff_protein".to_string(), protein);
    let mut app = GENtleApp::default();
    app.engine = Arc::new(RwLock::new(GentleEngine::from_state(state)));
    app.uniprot_map_seq_id = "handoff_target".to_string();
    app.reverse_translate_protein_seq_id = "handoff_protein".to_string();
    app.reverse_translate_speed_profile = Some(TranslationSpeedProfile::Ecoli);
    app.reverse_translate_speed_mark = Some(TranslationSpeedMark::Slow);
    app.reverse_translate_translation_table = "11".to_string();
    app.reverse_translate_target_anneal_tm_c = "58.0".to_string();
    app.reverse_translate_anneal_window_bp = "9".to_string();

    app.build_protein_to_dna_handoff_reasoning_from_dialog();

    let graph = app
        .protein_handoff_graph
        .as_ref()
        .expect("protein handoff graph");
    assert_eq!(graph.seq_id, "handoff_target");
    assert!(!app.protein_handoff_selected_candidate_id.is_empty());
    assert!(
        app.protein_handoff_status
            .contains("Protein-to-DNA handoff reasoning: ok")
    );
    assert!(graph.candidates.iter().any(|candidate| {
        candidate
            .protein_to_dna_handoff
            .as_ref()
            .is_some_and(|detail| {
                detail.strategy == ProteinToDnaHandoffStrategy::ReverseTranslatedSynthetic
            })
    }));
}

#[test]
fn open_lineage_analysis_artifact_opens_protein_handoff_in_protein_evidence() {
    let mut protein = DNAsequence::from_sequence("MKP").expect("protein");
    protein.set_name("Toy protein");
    protein.set_molecule_type("protein");
    let mut state = ProjectState::default();
    state.sequences.insert(
        "handoff_target".to_string(),
        DNAsequence::from_sequence("ACGTACGTACGT").expect("dna"),
    );
    state
        .sequences
        .insert("handoff_protein".to_string(), protein);
    let mut app = GENtleApp::default();
    app.engine = Arc::new(RwLock::new(GentleEngine::from_state(state)));

    let graph_id = {
        let mut engine = app.engine.write().unwrap();
        engine
            .apply(Operation::BuildProteinToDnaHandoffReasoning {
                seq_id: "handoff_target".to_string(),
                protein_seq_id: "handoff_protein".to_string(),
                transcript_filter: None,
                projection_id: None,
                ensembl_entry_id: Some("ENSPTOY1".to_string()),
                feature_query: None,
                ranking_goal: ProteinToDnaHandoffRankingGoal::BalancedProvenance,
                speed_profile: Some(TranslationSpeedProfile::Ecoli),
                speed_mark: Some(TranslationSpeedMark::Slow),
                translation_table: Some(11),
                target_anneal_tm_c: Some(58.0),
                anneal_window_bp: Some(9),
                objective_id: None,
                graph_id: Some("handoff_graph".to_string()),
            })
            .expect("build protein handoff graph")
            .construct_reasoning_graph
            .expect("graph")
            .graph_id
    };

    app.open_lineage_analysis_artifact(
        LineageAnalysisKind::ConstructReasoning,
        "handoff_target",
        &graph_id,
    );

    assert!(app.show_uniprot_dialog);
    assert!(app.new_windows.is_empty());
    assert_eq!(
        app.protein_handoff_graph
            .as_ref()
            .map(|graph| graph.graph_id.as_str()),
        Some(graph_id.as_str())
    );
    assert_eq!(app.reverse_translate_protein_seq_id, "handoff_protein");
    assert!(!app.protein_handoff_selected_candidate_id.is_empty());
    assert!(
        app.protein_handoff_status
            .contains("Opened protein-to-DNA handoff reasoning graph")
    );
}

#[test]
fn open_lineage_analysis_artifact_opens_uniprot_projection_expert() {
    let mut state = ProjectState::default();
    state.sequences.insert(
        "seq_uniprot".to_string(),
        DNAsequence::from_sequence("ACGTACGT").expect("sequence"),
    );
    let mut app = GENtleApp::default();
    app.engine = Arc::new(RwLock::new(GentleEngine::from_state(state)));

    app.open_lineage_analysis_artifact(
        LineageAnalysisKind::UniprotProjection,
        "seq_uniprot",
        "projection_1",
    );

    assert_eq!(app.new_windows.len(), 1);
}

#[test]
fn routine_assistant_detects_circular_gibson_binding() {
    let mut state = ProjectState::default();
    let mut vector = DNAsequence::from_sequence("ACGTACGTACGT").expect("vector");
    vector.set_circular(true);
    state.sequences.insert("vector".to_string(), vector);
    state.sequences.insert(
        "insert".to_string(),
        DNAsequence::from_sequence("TTTTGGGG").expect("insert"),
    );
    let engine = Arc::new(RwLock::new(GentleEngine::from_state(state)));
    let mut app = GENtleApp::default();
    app.engine = engine;

    let routine = make_test_gibson_routine_row();
    app.routine_assistant_bindings
        .insert("left_seq_id".to_string(), "insert".to_string());
    app.routine_assistant_bindings
        .insert("right_seq_id".to_string(), "vector".to_string());

    let circular = app
        .routine_assistant_gibson_circular_binding_for_routine(&routine)
        .expect("circular Gibson binding");
    assert_eq!(circular.port_id, "right_seq_id");
    assert_eq!(circular.seq_id, "vector");
    assert!(circular.circular);
    assert_eq!(circular.length_bp, 12);
}

#[test]
fn routine_assistant_linearize_binding_sequence_rebinds_to_linear_branch() {
    let mut state = ProjectState::default();
    let mut vector = DNAsequence::from_sequence("ACGTACGTACGT").expect("vector");
    vector.set_circular(true);
    state.sequences.insert("vector".to_string(), vector);
    let engine = Arc::new(RwLock::new(GentleEngine::from_state(state)));
    let mut app = GENtleApp::default();
    app.engine = engine.clone();
    app.routine_assistant_bindings
        .insert("right_seq_id".to_string(), "vector".to_string());

    let new_seq_id = app
        .routine_assistant_linearize_binding_sequence("right_seq_id", "vector")
        .expect("linearized sequence");

    assert_eq!(
        app.routine_assistant_bindings.get("right_seq_id"),
        Some(&new_seq_id)
    );
    let guard = engine.read().expect("engine lock");
    let original = guard
        .state()
        .sequences
        .get("vector")
        .expect("original sequence");
    assert!(original.is_circular());
    let linearized = guard
        .state()
        .sequences
        .get(new_seq_id.as_str())
        .expect("linearized sequence");
    assert!(!linearized.is_circular());
}

#[test]
fn routine_assistant_preflight_blocks_gibson_with_circular_input() {
    let mut state = ProjectState::default();
    let mut vector = DNAsequence::from_sequence("ACGTACGTACGT").expect("vector");
    vector.set_circular(true);
    state.sequences.insert("vector".to_string(), vector);
    state.sequences.insert(
        "insert".to_string(),
        DNAsequence::from_sequence("TTTTGGGG").expect("insert"),
    );
    let engine = Arc::new(RwLock::new(GentleEngine::from_state(state)));
    let mut app = GENtleApp::default();
    app.engine = engine;
    let routine = make_test_gibson_routine_row();
    app.routine_assistant_selected_routine_id = routine.routine_id.clone();
    app.routine_assistant_candidates = vec![routine];
    app.routine_assistant_bindings
        .insert("left_seq_id".to_string(), "insert".to_string());
    app.routine_assistant_bindings
        .insert("right_seq_id".to_string(), "vector".to_string());

    app.run_routine_assistant_preflight();

    assert!(matches!(
        app.routine_assistant_stage,
        RoutineAssistantStage::Preflight
    ));
    let output = app
        .routine_assistant_preflight_output
        .as_ref()
        .expect("blocking preflight output");
    assert!(
        !output
            .get("can_execute")
            .and_then(|value| value.as_bool())
            .unwrap_or(true)
    );
    let first_error = output
        .get("preflight")
        .and_then(|value| value.get("errors"))
        .and_then(|value| value.as_array())
        .and_then(|rows| rows.first())
        .and_then(|value| value.as_str())
        .unwrap_or("");
    assert!(first_error.contains("requires linear fragments"));
}

#[test]
fn routine_assistant_execute_blocks_gibson_with_circular_input() {
    let mut state = ProjectState::default();
    let mut vector = DNAsequence::from_sequence("ACGTACGTACGT").expect("vector");
    vector.set_circular(true);
    state.sequences.insert("vector".to_string(), vector);
    state.sequences.insert(
        "insert".to_string(),
        DNAsequence::from_sequence("TTTTGGGG").expect("insert"),
    );
    let engine = Arc::new(RwLock::new(GentleEngine::from_state(state)));
    let mut app = GENtleApp::default();
    app.engine = engine;
    let routine = make_test_gibson_routine_row();
    app.routine_assistant_selected_routine_id = routine.routine_id.clone();
    app.routine_assistant_candidates = vec![routine];
    app.routine_assistant_bindings
        .insert("left_seq_id".to_string(), "insert".to_string());
    app.routine_assistant_bindings
        .insert("right_seq_id".to_string(), "vector".to_string());

    app.run_routine_assistant_execute();

    assert!(matches!(
        app.routine_assistant_stage,
        RoutineAssistantStage::Preflight
    ));
    assert!(app.routine_assistant_execute_output.is_none());
    let output = app
        .routine_assistant_preflight_output
        .as_ref()
        .expect("blocking preflight output");
    assert!(
        !output
            .get("can_execute")
            .and_then(|value| value.as_bool())
            .unwrap_or(true)
    );
}

#[test]
fn native_window_menu_sync_is_blocked_while_open_probe_pending() {
    let mut app = GENtleApp::default();
    let viewport_id = GENtleApp::help_viewport_id();
    assert!(!app.native_window_menu_sync_blocked_by_open_probe());

    app.mark_viewport_open_requested(viewport_id);
    assert!(app.native_window_menu_sync_blocked_by_open_probe());

    app.finalize_viewport_open_probe(viewport_id, "Help");
    assert!(!app.native_window_menu_sync_blocked_by_open_probe());
}

#[test]
fn queue_focus_viewport_records_and_finalizes_focus_probe() {
    let mut app = GENtleApp::default();
    let viewport_id = GENtleApp::help_viewport_id();
    app.queue_focus_viewport(viewport_id);
    assert!(
        app.pending_viewport_focus_timestamps
            .contains_key(&viewport_id)
    );

    app.finalize_viewport_focus_probe(viewport_id);
    assert!(
        !app.pending_viewport_focus_timestamps
            .contains_key(&viewport_id)
    );
}

#[test]
fn embedded_window_foreground_request_is_one_shot_per_render() {
    let mut app = GENtleApp::default();
    for viewport_id in [
        GENtleApp::help_viewport_id(),
        GENtleApp::agent_assistant_viewport_id(),
    ] {
        app.queue_focus_viewport(viewport_id);
        assert!(app.viewport_foreground_requested(viewport_id));

        app.clear_viewport_foreground_request_after_render(viewport_id);

        assert!(
            !app.viewport_foreground_requested(viewport_id),
            "foreground request must clear after one embedded render for {viewport_id:?}"
        );
    }

    let dna = DNAsequence::from_sequence("ACGT").expect("sequence");
    let viewport_id =
        app.register_window(Window::new_dna(dna, "seq1".to_string(), app.engine.clone()));
    app.queue_focus_viewport(viewport_id);
    assert!(app.viewport_foreground_requested(viewport_id));

    app.clear_viewport_foreground_request_after_render(viewport_id);

    assert!(
        !app.viewport_foreground_requested(viewport_id),
        "sequence-window foreground request must clear after one embedded render"
    );
}

#[test]
fn mark_window_open_or_focus_marks_open_and_queues_focus_for_closed_window() {
    let mut app = GENtleApp::default();
    let viewport_id = GENtleApp::planning_viewport_id();

    app.mark_window_open_or_focus(viewport_id, false);

    assert!(
        app.pending_window_open_timestamps
            .contains_key(&viewport_id)
    );
    assert!(app.pending_focus_viewports.contains(&viewport_id));
    assert!(
        app.pending_viewport_focus_timestamps
            .contains_key(&viewport_id)
    );
}

#[test]
fn mark_window_open_or_focus_focuses_existing_window_without_new_open_probe() {
    let mut app = GENtleApp::default();
    let viewport_id = GENtleApp::planning_viewport_id();

    app.mark_window_open_or_focus(viewport_id, true);

    assert!(
        !app.pending_window_open_timestamps
            .contains_key(&viewport_id)
    );
    assert!(app.pending_focus_viewports.contains(&viewport_id));
    assert!(
        app.pending_viewport_focus_timestamps
            .contains_key(&viewport_id)
    );
}

#[test]
fn apply_graphics_settings_to_display_clamps_font_and_opacity_values() {
    let mut source = DisplaySettings::default();
    let mut target = DisplaySettings::default();

    source.show_linear_sequence_panel = true;
    source.feature_details_font_size = 250.0;
    source.sequence_panel_max_text_length_bp = 10_000_000;
    source.linear_external_feature_label_font_size = -50.0;
    source.linear_external_feature_label_background_opacity = 2.0;
    source.reverse_strand_visual_opacity = 10.0;
    source.linear_sequence_helical_phase_offset_bp = 27;
    source.restriction_enzyme_display_mode = RestrictionEnzymeDisplayMode::UniqueOnly;
    source.preferred_restriction_enzymes = vec![
        " EcoRI ".to_string(),
        "ecori".to_string(),
        "BamHI".to_string(),
    ];

    GENtleApp::apply_graphics_settings_to_display(&source, &mut target);

    assert!(target.show_linear_sequence_panel);
    assert_eq!(target.feature_details_font_size, 24.0);
    assert_eq!(target.sequence_panel_max_text_length_bp, 5_000_000);
    assert_eq!(target.linear_external_feature_label_font_size, 8.0);
    assert_eq!(target.linear_external_feature_label_background_opacity, 1.0);
    assert_eq!(target.reverse_strand_visual_opacity, 1.0);
    assert_eq!(target.linear_sequence_helical_phase_offset_bp, 7);
    assert_eq!(
        target.restriction_enzyme_display_mode,
        RestrictionEnzymeDisplayMode::UniqueOnly
    );
    assert_eq!(
        target.preferred_restriction_enzymes,
        vec!["EcoRI".to_string(), "BamHI".to_string()]
    );

    source.feature_details_font_size = f32::NAN;
    source.sequence_panel_max_text_length_bp = 0;
    source.linear_external_feature_label_font_size = f32::NAN;
    source.linear_external_feature_label_background_opacity = f32::NAN;
    source.reverse_strand_visual_opacity = f32::NAN;
    source.linear_sequence_helical_phase_offset_bp = 10;
    GENtleApp::apply_graphics_settings_to_display(&source, &mut target);

    assert_eq!(target.feature_details_font_size, 9.0);
    assert_eq!(target.sequence_panel_max_text_length_bp, 0);
    assert_eq!(target.linear_external_feature_label_font_size, 11.0);
    assert_eq!(target.linear_external_feature_label_background_opacity, 0.9);
    assert_eq!(
        target.reverse_strand_visual_opacity,
        DisplaySettings::default_reverse_strand_visual_opacity()
    );
    assert_eq!(target.linear_sequence_helical_phase_offset_bp, 0);
}

#[test]
fn default_configuration_graphics_prefers_adaptive_linear_letters() {
    let app = GENtleApp::default();
    assert!(!app.configuration_graphics.show_linear_sequence_panel);
    assert_eq!(
        app.configuration_graphics
            .linear_sequence_letter_layout_mode,
        LinearSequenceLetterLayoutMode::AutoAdaptive
    );
    assert!(
        app.configuration_graphics
            .linear_sequence_helical_letters_enabled
    );
}

#[test]
fn apply_configuration_graphics_to_engine_state_keeps_adaptive_linear_defaults() {
    let mut app = GENtleApp::default();
    app.configuration_graphics = DisplaySettings::default();
    app.configuration_graphics_dirty = true;

    app.apply_configuration_graphics_to_engine_state();

    let guard = app.engine.read().expect("engine lock");
    assert_eq!(
        guard.state().display.linear_sequence_letter_layout_mode,
        LinearSequenceLetterLayoutMode::AutoAdaptive
    );
    assert!(
        guard
            .state()
            .display
            .linear_sequence_helical_letters_enabled
    );
    assert!(!app.configuration_graphics_dirty);
}

#[test]
fn upgrade_persisted_configuration_migrates_legacy_linear_letter_defaults() {
    let mut saved = PersistedConfiguration::default();
    saved.schema_version = 1;
    saved
        .graphics_defaults
        .linear_sequence_helical_letters_enabled = false;
    saved.graphics_defaults.linear_sequence_letter_layout_mode =
        LinearSequenceLetterLayoutMode::ContinuousHelical;

    let upgraded = GENtleApp::upgrade_persisted_configuration(&mut saved);

    assert!(upgraded);
    assert_eq!(saved.schema_version, APP_CONFIGURATION_SCHEMA_VERSION);
    assert!(
        saved
            .graphics_defaults
            .linear_sequence_helical_letters_enabled
    );
    assert_eq!(
        saved.graphics_defaults.linear_sequence_letter_layout_mode,
        LinearSequenceLetterLayoutMode::AutoAdaptive
    );
}

#[test]
fn upgrade_persisted_configuration_keeps_default_ui_language_for_legacy_settings() {
    let mut saved = PersistedConfiguration::default();
    saved.schema_version = 2;
    saved.ui_language = UiLanguage::default();

    let upgraded = GENtleApp::upgrade_persisted_configuration(&mut saved);

    assert!(upgraded);
    assert_eq!(saved.schema_version, APP_CONFIGURATION_SCHEMA_VERSION);
    assert_eq!(saved.ui_language, UiLanguage::System);
}

#[test]
fn upgrade_persisted_configuration_is_noop_for_current_schema() {
    let mut saved = PersistedConfiguration::default();
    saved.schema_version = APP_CONFIGURATION_SCHEMA_VERSION;
    saved
        .graphics_defaults
        .linear_sequence_helical_letters_enabled = false;
    saved.graphics_defaults.linear_sequence_letter_layout_mode =
        LinearSequenceLetterLayoutMode::Condensed10Row;

    let upgraded = GENtleApp::upgrade_persisted_configuration(&mut saved);

    assert!(!upgraded);
    assert_eq!(saved.schema_version, APP_CONFIGURATION_SCHEMA_VERSION);
    assert!(
        !saved
            .graphics_defaults
            .linear_sequence_helical_letters_enabled
    );
    assert_eq!(
        saved.graphics_defaults.linear_sequence_letter_layout_mode,
        LinearSequenceLetterLayoutMode::Condensed10Row
    );
}

#[test]
fn load_help_doc_uses_runtime_file_and_rewrites_relative_images() {
    let temp = tempdir().unwrap();
    let docs_dir = temp.path().join("docs");
    let images_dir = docs_dir.join("images");
    fs::create_dir_all(&images_dir).unwrap();

    let image_path = images_dir.join("gui.png");
    fs::write(&image_path, b"fake image").unwrap();

    let markdown_path = docs_dir.join("gui.md");
    fs::write(
        &markdown_path,
        "# Help\n\n![GUI](<images/gui.png> \"GUI Screenshot\")\n",
    )
    .unwrap();

    let loaded = GENtleApp::load_help_doc(markdown_path.to_str().unwrap(), "fallback");
    let abs_image = image_path.canonicalize().unwrap();
    let image_uri = GENtleApp::file_uri_from_path(&abs_image);
    let expected = format!("![GUI](<{}> \"GUI Screenshot\")", image_uri);
    assert!(loaded.contains(&expected), "{loaded}");
}

#[test]
fn rewrite_markdown_relative_image_links_handles_paths_with_parentheses() {
    let temp = tempdir().unwrap();
    let docs_dir = temp.path().join("docs");
    let images_dir = docs_dir.join("images");
    fs::create_dir_all(&images_dir).unwrap();

    let image_path = images_dir.join("gui(1).png");
    fs::write(&image_path, b"fake image").unwrap();

    let markdown = "![Shot](images/gui(1).png)\n";
    let rewritten = GENtleApp::rewrite_markdown_relative_image_links(markdown, &docs_dir);
    let abs_image = image_path.canonicalize().unwrap();
    let expected = format!("![Shot]({})", GENtleApp::file_uri_from_path(&abs_image));
    assert!(rewritten.contains(&expected), "{rewritten}");
}

#[test]
fn rewrite_markdown_relative_image_links_keeps_absolute_and_reference_images() {
    let markdown = "![Web](https://example.com/gui.png)\n![Root](/tmp/gui.png)\n![ByRef][img]\n\n[img]: images/gui.png\n";
    let temp = tempdir().unwrap();
    let rewritten = GENtleApp::rewrite_markdown_relative_image_links(markdown, temp.path());
    assert_eq!(rewritten, markdown);
}

#[test]
fn normalize_recent_project_paths_deduplicates_and_limits() {
    let mut paths = Vec::new();
    for idx in 0..(MAX_RECENT_PROJECTS + 4) {
        paths.push(format!("/tmp/gentle_project_{idx}.gentle.json"));
    }
    paths.insert(2, "/tmp/gentle_project_1.gentle.json".to_string());
    paths.insert(4, "   ".to_string());

    let normalized = GENtleApp::normalize_recent_project_paths(paths);
    assert_eq!(normalized.len(), MAX_RECENT_PROJECTS);
    assert_eq!(
        normalized[0],
        GENtleApp::normalize_project_path("/tmp/gentle_project_0.gentle.json")
    );
    assert_eq!(
        normalized[1],
        GENtleApp::normalize_project_path("/tmp/gentle_project_1.gentle.json")
    );
    let needle = GENtleApp::normalize_project_path("/tmp/gentle_project_1.gentle.json");
    assert_eq!(normalized.iter().filter(|p| *p == &needle).count(), 1);
}

#[test]
fn recent_project_menu_label_shows_name_and_parent() {
    let temp = tempdir().unwrap();
    let project_path = temp.path().join("my_project.gentle.json");
    let label = GENtleApp::recent_project_menu_label(project_path.to_string_lossy().as_ref());
    assert!(label.starts_with("my_project.gentle.json ("));
}

#[test]
fn arrangement_lane_helpers_preserve_order_and_describe_pools() {
    let mut app = GENtleApp::default();
    app.lineage_containers = vec![
        ContainerRow {
            container_id: "container-1".to_string(),
            kind: "Singleton".to_string(),
            declared_contents_exclusive: true,
            member_count: 1,
            representative: "gibson_destination_pgex".to_string(),
            members: vec!["gibson_destination_pgex".to_string()],
        },
        ContainerRow {
            container_id: "container-2".to_string(),
            kind: "Pool".to_string(),
            declared_contents_exclusive: false,
            member_count: 3,
            representative: "pooled_product".to_string(),
            members: vec!["a".to_string(), "b".to_string(), "c".to_string()],
        },
    ];

    let rows = app
        .arrangement_lane_container_rows(&["container-2".to_string(), "container-1".to_string()]);
    assert_eq!(rows.len(), 2);
    assert_eq!(rows[0].container_id, "container-2");
    assert_eq!(rows[1].container_id, "container-1");
    assert_eq!(
        GENtleApp::arrangement_lane_menu_label(rows[0]),
        "container-2 (Pool, 3 seqs)"
    );
    assert_eq!(
        GENtleApp::arrangement_lane_menu_label(rows[1]),
        "container-1 (gibson_destination_pgex)"
    );
}

#[test]
fn container_contents_mode_labels_are_stable() {
    assert_eq!(
        GENtleApp::container_contents_mode_label(true),
        "Declared only"
    );
    assert_eq!(
        GENtleApp::container_contents_mode_label(false),
        "Known subset"
    );
}

#[test]
fn arrangement_gel_preview_coerces_one_sided_choice_to_symmetric_ladder() {
    let mut app = GENtleApp::default();
    app.arrangement_gel_preview.left_ladder_name.clear();
    app.arrangement_gel_preview.right_ladder_name = "NEB 1kb DNA Ladder".to_string();

    app.coerce_arrangement_gel_preview_pair();

    assert_eq!(
        app.arrangement_gel_preview.left_ladder_name,
        "NEB 1kb DNA Ladder"
    );
    assert_eq!(
        app.arrangement_gel_preview_effective_ladders(),
        vec!["NEB 1kb DNA Ladder".to_string()]
    );
}

#[test]
fn open_arrangement_gel_preview_dialog_seeds_ladder_controls_and_preview() {
    let mut state = ProjectState::default();
    state.sequences.insert(
        "seq_a".to_string(),
        DNAsequence::from_sequence(&"ATGC".repeat(90)).unwrap(),
    );
    state.sequences.insert(
        "seq_b".to_string(),
        DNAsequence::from_sequence(&"ATGC".repeat(35)).unwrap(),
    );
    state.container_state.containers.insert(
        "container-1".to_string(),
        Container {
            container_id: "container-1".to_string(),
            kind: ContainerKind::Singleton,
            name: Some("Vector".to_string()),
            members: vec!["seq_a".to_string()],
            declared_contents_exclusive: true,
            created_by_op: None,
            created_at_unix_ms: 0,
        },
    );
    state.container_state.containers.insert(
        "container-2".to_string(),
        Container {
            container_id: "container-2".to_string(),
            kind: ContainerKind::Singleton,
            name: Some("Insert".to_string()),
            members: vec!["seq_b".to_string()],
            declared_contents_exclusive: true,
            created_by_op: None,
            created_at_unix_ms: 0,
        },
    );
    state.container_state.arrangements.insert(
        "arr-1".to_string(),
        Arrangement {
            arrangement_id: "arr-1".to_string(),
            mode: ArrangementMode::Serial,
            name: Some("Demo run".to_string()),
            lane_container_ids: vec!["container-1".to_string(), "container-2".to_string()],
            ladders: vec![
                "NEB 1kb DNA Ladder".to_string(),
                "GeneRuler 100bp DNA Ladder Plus".to_string(),
            ],
            lane_role_labels: vec!["vector".to_string(), "insert_1".to_string()],
            default_rack_id: None,
            created_by_op: None,
            created_at_unix_ms: 0,
        },
    );
    let mut app = GENtleApp::default();
    app.engine = Arc::new(RwLock::new(GentleEngine::from_state(state)));

    app.open_arrangement_gel_preview_dialog("arr-1");

    assert!(app.show_arrangement_gel_preview_dialog);
    assert_eq!(app.arrangement_gel_preview.arrangement_id, "arr-1");
    assert_eq!(app.arrangement_gel_preview.arrangement_title, "Demo run");
    assert_eq!(
        app.arrangement_gel_preview.left_ladder_name,
        "NEB 1kb DNA Ladder"
    );
    assert_eq!(
        app.arrangement_gel_preview.right_ladder_name,
        "GeneRuler 100bp DNA Ladder Plus"
    );
    assert!(!app.arrangement_gel_preview.svg_uri.is_empty());
}

#[test]
fn open_arrangement_labels_preview_dialog_seeds_scope_and_preview() {
    let mut state = ProjectState::default();
    state.sequences.insert(
        "seq_a".to_string(),
        DNAsequence::from_sequence(&"ATGC".repeat(40)).unwrap(),
    );
    state.container_state.containers.insert(
        "container-1".to_string(),
        Container {
            container_id: "container-1".to_string(),
            kind: ContainerKind::Singleton,
            name: Some("Vector".to_string()),
            members: vec!["seq_a".to_string()],
            declared_contents_exclusive: true,
            created_by_op: None,
            created_at_unix_ms: 0,
        },
    );
    state.container_state.arrangements.insert(
        "arr-1".to_string(),
        Arrangement {
            arrangement_id: "arr-1".to_string(),
            mode: ArrangementMode::Serial,
            name: Some("Demo labels".to_string()),
            lane_container_ids: vec!["container-1".to_string()],
            ladders: vec![],
            lane_role_labels: vec!["vector".to_string()],
            default_rack_id: None,
            created_by_op: None,
            created_at_unix_ms: 0,
        },
    );
    let mut app = GENtleApp::default();
    app.engine = Arc::new(RwLock::new(GentleEngine::from_state(state)));

    app.open_arrangement_labels_preview_dialog("arr-1");

    assert!(app.show_rack_labels_preview_dialog);
    assert_eq!(
        app.rack_labels_preview.arrangement_id.as_deref(),
        Some("arr-1")
    );
    assert_eq!(
        app.rack_labels_preview.arrangement_title.as_deref(),
        Some("Demo labels")
    );
    assert!(!app.rack_labels_preview.rack_id.is_empty());
    assert!(!app.rack_labels_preview.svg_uri.is_empty());
}

#[test]
fn open_arrangement_rack_dialog_materializes_default_rack_for_legacy_arrangement() {
    let mut state = ProjectState::default();
    state.sequences.insert(
        "seq_a".to_string(),
        DNAsequence::from_sequence("ATGCATGC").unwrap(),
    );
    state.container_state.containers.insert(
        "container-1".to_string(),
        Container {
            container_id: "container-1".to_string(),
            kind: ContainerKind::Singleton,
            name: Some("Vector".to_string()),
            members: vec!["seq_a".to_string()],
            declared_contents_exclusive: true,
            created_by_op: None,
            created_at_unix_ms: 0,
        },
    );
    state.container_state.arrangements.insert(
        "arr-legacy".to_string(),
        Arrangement {
            arrangement_id: "arr-legacy".to_string(),
            mode: ArrangementMode::Serial,
            name: Some("Legacy".to_string()),
            lane_container_ids: vec!["container-1".to_string()],
            ladders: vec![],
            lane_role_labels: vec!["vector".to_string()],
            default_rack_id: None,
            created_by_op: None,
            created_at_unix_ms: 0,
        },
    );
    let mut app = GENtleApp::default();
    app.engine = Arc::new(RwLock::new(GentleEngine::from_state(state)));

    app.open_arrangement_rack_dialog("arr-legacy");

    assert!(app.show_rack_dialog);
    assert!(!app.rack_view_rack_id.is_empty());
    let engine = app.engine.read().unwrap();
    assert!(
        engine
            .state()
            .container_state
            .racks
            .contains_key(&app.rack_view_rack_id)
    );
}

#[test]
fn open_rack_dialog_prefills_custom_profile_editor_from_snapshot() {
    let mut state = ProjectState::default();
    state.container_state.racks.insert(
        "rack-custom".to_string(),
        Rack {
            rack_id: "rack-custom".to_string(),
            name: "Bench".to_string(),
            profile: RackProfileSnapshot {
                fill_direction: RackFillDirection::ColumnMajor,
                blocked_coordinates: vec!["B2".to_string(), "AA3".to_string()],
                ..RackProfileSnapshot::custom(28, 10)
            },
            placements: vec![],
            created_by_op: None,
            created_at_unix_ms: 0,
        },
    );
    let mut app = GENtleApp::default();
    app.engine = Arc::new(RwLock::new(GentleEngine::from_state(state)));

    app.open_rack_dialog("rack-custom");

    assert!(app.show_rack_dialog);
    assert_eq!(app.rack_view_rack_id, "rack-custom");
    assert_eq!(app.rack_profile_editor_kind, RackProfileKind::Custom);
    assert_eq!(
        app.rack_authoring_template_editor,
        RackAuthoringTemplate::PlateColumns
    );
    assert_eq!(
        app.rack_fill_direction_editor,
        RackFillDirection::ColumnMajor
    );
    assert_eq!(app.rack_custom_profile_rows, "28");
    assert_eq!(app.rack_custom_profile_columns, "10");
    assert_eq!(app.rack_blocked_coordinates_text, "B2, AA3");
}

#[test]
fn open_rack_dialog_infers_edge_avoidance_template_from_perimeter_blocking() {
    let mut state = ProjectState::default();
    state.container_state.racks.insert(
        "rack-edge".to_string(),
        Rack {
            rack_id: "rack-edge".to_string(),
            name: "Plate".to_string(),
            profile: RackProfileSnapshot {
                fill_direction: RackFillDirection::ColumnMajor,
                blocked_coordinates: vec![
                    "A1".to_string(),
                    "B1".to_string(),
                    "C1".to_string(),
                    "D1".to_string(),
                    "A2".to_string(),
                    "D2".to_string(),
                    "A3".to_string(),
                    "D3".to_string(),
                    "A4".to_string(),
                    "B4".to_string(),
                    "C4".to_string(),
                    "D4".to_string(),
                ],
                ..RackProfileSnapshot::custom(4, 4)
            },
            placements: vec![],
            created_by_op: None,
            created_at_unix_ms: 0,
        },
    );
    let mut app = GENtleApp::default();
    app.engine = Arc::new(RwLock::new(GentleEngine::from_state(state)));

    app.open_rack_dialog("rack-edge");

    assert_eq!(
        app.rack_authoring_template_editor,
        RackAuthoringTemplate::PlateEdgeAvoidance
    );
}

#[test]
fn rack_first_coordinate_for_arrangement_prefers_earliest_slot() {
    let entries = vec![
        (
            3,
            "A4".to_string(),
            crate::engine::RackPlacementEntry {
                coordinate: "A4".to_string(),
                occupant: None,
                arrangement_id: "arr-a".to_string(),
                order_index: 1,
                role_label: "insert_1".to_string(),
            },
        ),
        (
            1,
            "A2".to_string(),
            crate::engine::RackPlacementEntry {
                coordinate: "A2".to_string(),
                occupant: None,
                arrangement_id: "arr-a".to_string(),
                order_index: 0,
                role_label: "vector".to_string(),
            },
        ),
        (
            2,
            "A3".to_string(),
            crate::engine::RackPlacementEntry {
                coordinate: "A3".to_string(),
                occupant: None,
                arrangement_id: "arr-b".to_string(),
                order_index: 0,
                role_label: "vector".to_string(),
            },
        ),
    ];
    assert_eq!(
        GENtleApp::rack_first_coordinate_for_arrangement(&entries, "arr-a"),
        Some("A2".to_string())
    );
    assert_eq!(
        GENtleApp::rack_first_coordinate_for_arrangement(&entries, "arr-b"),
        Some("A3".to_string())
    );
    assert_eq!(
        GENtleApp::rack_first_coordinate_for_arrangement(&entries, "missing"),
        None
    );
}

#[test]
fn rack_selected_arrangement_ids_in_order_preserves_rack_order() {
    let entries = vec![
        (
            0,
            "A1".to_string(),
            crate::engine::RackPlacementEntry {
                coordinate: "A1".to_string(),
                occupant: None,
                arrangement_id: "arr-a".to_string(),
                order_index: 0,
                role_label: "lane_1".to_string(),
            },
        ),
        (
            1,
            "A2".to_string(),
            crate::engine::RackPlacementEntry {
                coordinate: "A2".to_string(),
                occupant: None,
                arrangement_id: "arr-b".to_string(),
                order_index: 0,
                role_label: "lane_1".to_string(),
            },
        ),
        (
            2,
            "A3".to_string(),
            crate::engine::RackPlacementEntry {
                coordinate: "A3".to_string(),
                occupant: None,
                arrangement_id: "arr-c".to_string(),
                order_index: 0,
                role_label: "lane_1".to_string(),
            },
        ),
    ];
    let selected = ["arr-c".to_string(), "arr-a".to_string()]
        .into_iter()
        .collect::<std::collections::BTreeSet<_>>();
    assert_eq!(
        GENtleApp::rack_selected_arrangement_ids_in_order(&entries, &selected),
        vec!["arr-a".to_string(), "arr-c".to_string()]
    );
}

#[test]
fn rack_selected_coordinates_in_order_preserves_rack_order() {
    let entries = vec![
        (
            0,
            "A1".to_string(),
            crate::engine::RackPlacementEntry {
                coordinate: "A1".to_string(),
                occupant: None,
                arrangement_id: "arr-a".to_string(),
                order_index: 0,
                role_label: "lane_1".to_string(),
            },
        ),
        (
            1,
            "A2".to_string(),
            crate::engine::RackPlacementEntry {
                coordinate: "A2".to_string(),
                occupant: None,
                arrangement_id: "arr-a".to_string(),
                order_index: 1,
                role_label: "lane_2".to_string(),
            },
        ),
        (
            2,
            "A3".to_string(),
            crate::engine::RackPlacementEntry {
                coordinate: "A3".to_string(),
                occupant: None,
                arrangement_id: "arr-b".to_string(),
                order_index: 0,
                role_label: "lane_1".to_string(),
            },
        ),
    ];
    let selected = ["A3".to_string(), "A1".to_string()]
        .into_iter()
        .collect::<std::collections::BTreeSet<_>>();
    assert_eq!(
        GENtleApp::rack_selected_coordinates_in_order(&entries, &selected),
        vec!["A1".to_string(), "A3".to_string()]
    );
}

#[test]
fn rack_selected_samples_arrangement_id_requires_one_block() {
    let entries = vec![
        (
            0,
            "A1".to_string(),
            crate::engine::RackPlacementEntry {
                coordinate: "A1".to_string(),
                occupant: None,
                arrangement_id: "arr-a".to_string(),
                order_index: 0,
                role_label: "lane_1".to_string(),
            },
        ),
        (
            1,
            "A2".to_string(),
            crate::engine::RackPlacementEntry {
                coordinate: "A2".to_string(),
                occupant: None,
                arrangement_id: "arr-a".to_string(),
                order_index: 1,
                role_label: "lane_2".to_string(),
            },
        ),
        (
            2,
            "A3".to_string(),
            crate::engine::RackPlacementEntry {
                coordinate: "A3".to_string(),
                occupant: None,
                arrangement_id: "arr-b".to_string(),
                order_index: 0,
                role_label: "lane_1".to_string(),
            },
        ),
    ];
    let same_block = ["A1".to_string(), "A2".to_string()]
        .into_iter()
        .collect::<std::collections::BTreeSet<_>>();
    let mixed = ["A1".to_string(), "A3".to_string()]
        .into_iter()
        .collect::<std::collections::BTreeSet<_>>();
    assert_eq!(
        GENtleApp::rack_selected_samples_arrangement_id(&entries, &same_block),
        Some("arr-a".to_string())
    );
    assert_eq!(
        GENtleApp::rack_selected_samples_arrangement_id(&entries, &mixed),
        None
    );
}

#[test]
fn rack_drag_status_text_mentions_drop_target_and_mode() {
    let sample_text = GENtleApp::rack_drag_status_text(
        &RackDragState::Sample {
            from_coordinate: "A2".to_string(),
            arrangement_id: "arr-a".to_string(),
            role_label: "insert_1".to_string(),
        },
        Some("B4"),
    );
    assert!(sample_text.contains("A2"));
    assert!(sample_text.contains("B4"));
    assert!(sample_text.contains("insert_1"));

    let multi_sample_text = GENtleApp::rack_drag_status_text(
        &RackDragState::Samples {
            arrangement_id: "arr-a".to_string(),
            from_coordinates: vec![
                "A1".to_string(),
                "A3".to_string(),
                "A5".to_string(),
                "A7".to_string(),
            ],
        },
        Some("A2"),
    );
    assert!(multi_sample_text.contains("4 samples"));
    assert!(multi_sample_text.contains("A1, A3, A5 (+1 more)"));
    assert!(multi_sample_text.contains("A2"));

    let block_text = GENtleApp::rack_drag_status_text(
        &RackDragState::ArrangementBlock {
            arrangement_id: "arr-b".to_string(),
            from_coordinate: "C1".to_string(),
        },
        Some("D6"),
    );
    assert!(block_text.contains("arr-b"));
    assert!(block_text.contains("C1"));
    assert!(block_text.contains("D6"));
    assert!(block_text.contains("whole block"));

    let multi_block_text = GENtleApp::rack_drag_status_text(
        &RackDragState::ArrangementBlocks {
            arrangement_ids: vec![
                "arr-a".to_string(),
                "arr-b".to_string(),
                "arr-c".to_string(),
                "arr-d".to_string(),
            ],
            from_coordinate: "A1".to_string(),
        },
        Some("B2"),
    );
    assert!(multi_block_text.contains("4 arrangement blocks"));
    assert!(multi_block_text.contains("arr-a, arr-b, arr-c (+1 more)"));
    assert!(multi_block_text.contains("B2"));
}

#[test]
fn rack_help_strip_items_cover_sample_and_block_modes() {
    let items = GENtleApp::rack_help_strip_items();
    assert!(items.iter().any(|(title, _)| *title == "Samples"));
    assert!(items.iter().any(|(title, _)| *title == "Blocks"));
    assert!(items.iter().any(|(_, body)| body.contains("fade")));
}

#[test]
fn rack_help_toggle_label_reflects_collapsed_state() {
    assert_eq!(GENtleApp::rack_help_toggle_label(false), "Hide help");
    assert_eq!(GENtleApp::rack_help_toggle_label(true), "Show help");
}

#[test]
fn rack_help_pin_label_reflects_pinned_state() {
    assert_eq!(GENtleApp::rack_help_pin_label(false), "Pin open");
    assert_eq!(GENtleApp::rack_help_pin_label(true), "Unpin");
}

#[test]
fn persist_rack_workspace_to_state_omits_default_and_stores_collapsed_flag() {
    let mut app = GENtleApp::default();
    app.persist_rack_workspace_to_state();
    assert!(
        app.engine
            .read()
            .unwrap()
            .state()
            .metadata
            .get(RACK_WORKSPACE_METADATA_KEY)
            .is_none()
    );

    app.rack_help_strip_collapsed = true;
    app.persist_rack_workspace_to_state();
    let workspace_value = app
        .engine
        .read()
        .unwrap()
        .state()
        .metadata
        .get(RACK_WORKSPACE_METADATA_KEY)
        .cloned()
        .expect("rack workspace metadata");
    let workspace: PersistedRackWorkspace =
        serde_json::from_value(workspace_value).expect("deserialize rack workspace");
    assert!(workspace.help_strip_collapsed);
}

#[test]
fn load_rack_workspace_from_state_restores_help_strip_collapsed_flag() {
    let mut state = ProjectState::default();
    state.metadata.insert(
        RACK_WORKSPACE_METADATA_KEY.to_string(),
        serde_json::json!({
            "help_strip_collapsed": true
        }),
    );
    let mut app = GENtleApp::default();
    app.engine = Arc::new(RwLock::new(GentleEngine::from_state(state)));
    app.load_rack_workspace_from_state();
    assert!(app.rack_help_strip_collapsed);
}

#[test]
fn rack_help_auto_minimizes_after_success_threshold_unless_pinned() {
    let mut app = GENtleApp::default();
    for _ in 0..RACK_HELP_AUTO_MINIMIZE_MOVE_THRESHOLD {
        app.record_successful_rack_move_and_maybe_autocollapse();
    }
    assert!(app.rack_help_strip_collapsed);
    assert!(app.rack_help_strip_auto_minimized);
    assert_eq!(
        app.rack_help_strip_successful_move_count,
        RACK_HELP_AUTO_MINIMIZE_MOVE_THRESHOLD
    );

    let mut pinned = GENtleApp::default();
    pinned.rack_help_strip_pinned_open = true;
    for _ in 0..RACK_HELP_AUTO_MINIMIZE_MOVE_THRESHOLD {
        pinned.record_successful_rack_move_and_maybe_autocollapse();
    }
    assert!(!pinned.rack_help_strip_collapsed);
    assert!(!pinned.rack_help_strip_auto_minimized);
    assert_eq!(pinned.rack_help_strip_successful_move_count, 0);
}

#[test]
fn rack_help_strip_keycaps_cover_shortcuts() {
    assert_eq!(
        GENtleApp::rack_help_strip_keycaps("Samples"),
        &["Cmd", "Ctrl"]
    );
    assert_eq!(
        GENtleApp::rack_help_strip_keycaps("Blocks"),
        &["Cmd", "Ctrl"]
    );
    assert_eq!(GENtleApp::rack_help_strip_keycaps("Keys"), &["Esc"]);
    assert!(GENtleApp::rack_help_strip_keycaps("Preview").is_empty());
}

#[test]
fn rack_ghost_preview_reorders_sample_within_arrangement_block() {
    let rack = crate::engine::Rack {
        rack_id: "rack-1".to_string(),
        name: "Bench".to_string(),
        profile: crate::engine::RackProfileSnapshot::custom(1, 4),
        placements: vec![
            crate::engine::RackPlacementEntry {
                coordinate: "A1".to_string(),
                occupant: Some(crate::engine::RackOccupant::Container {
                    container_id: "container-vector".to_string(),
                }),
                arrangement_id: "arr-a".to_string(),
                order_index: 0,
                role_label: "vector".to_string(),
            },
            crate::engine::RackPlacementEntry {
                coordinate: "A2".to_string(),
                occupant: Some(crate::engine::RackOccupant::Container {
                    container_id: "container-insert".to_string(),
                }),
                arrangement_id: "arr-a".to_string(),
                order_index: 1,
                role_label: "insert_1".to_string(),
            },
            crate::engine::RackPlacementEntry {
                coordinate: "A3".to_string(),
                occupant: Some(crate::engine::RackOccupant::Container {
                    container_id: "container-product".to_string(),
                }),
                arrangement_id: "arr-a".to_string(),
                order_index: 2,
                role_label: "product".to_string(),
            },
        ],
        created_by_op: None,
        created_at_unix_ms: 0,
    };
    let sorted_entries = GENtleApp::rack_sorted_entries(&rack);
    let preview = GENtleApp::rack_ghost_preview_map(
        &rack,
        &sorted_entries,
        &RackDragState::Sample {
            from_coordinate: "A1".to_string(),
            arrangement_id: "arr-a".to_string(),
            role_label: "vector".to_string(),
        },
        "A3",
    )
    .expect("preview");
    assert_eq!(preview.len(), 3);
    assert_eq!(
        preview
            .get("A1")
            .and_then(|cell| cell.predicted_entry.as_ref())
            .map(|entry| entry.role_label.as_str()),
        Some("insert_1")
    );
    assert_eq!(
        preview
            .get("A2")
            .and_then(|cell| cell.predicted_entry.as_ref())
            .map(|entry| entry.role_label.as_str()),
        Some("product")
    );
    assert_eq!(
        preview
            .get("A3")
            .and_then(|cell| cell.predicted_entry.as_ref())
            .map(|entry| entry.role_label.as_str()),
        Some("vector")
    );
}

#[test]
fn rack_ghost_preview_moves_whole_block_and_shifts_neighbors() {
    let rack = crate::engine::Rack {
        rack_id: "rack-1".to_string(),
        name: "Bench".to_string(),
        profile: crate::engine::RackProfileSnapshot::custom(1, 4),
        placements: vec![
            crate::engine::RackPlacementEntry {
                coordinate: "A1".to_string(),
                occupant: Some(crate::engine::RackOccupant::Container {
                    container_id: "container-a1".to_string(),
                }),
                arrangement_id: "arr-a".to_string(),
                order_index: 0,
                role_label: "vector".to_string(),
            },
            crate::engine::RackPlacementEntry {
                coordinate: "A2".to_string(),
                occupant: Some(crate::engine::RackOccupant::Container {
                    container_id: "container-a2".to_string(),
                }),
                arrangement_id: "arr-a".to_string(),
                order_index: 1,
                role_label: "product".to_string(),
            },
            crate::engine::RackPlacementEntry {
                coordinate: "A3".to_string(),
                occupant: Some(crate::engine::RackOccupant::Container {
                    container_id: "container-b1".to_string(),
                }),
                arrangement_id: "arr-b".to_string(),
                order_index: 0,
                role_label: "vector".to_string(),
            },
        ],
        created_by_op: None,
        created_at_unix_ms: 0,
    };
    let sorted_entries = GENtleApp::rack_sorted_entries(&rack);
    let preview = GENtleApp::rack_ghost_preview_map(
        &rack,
        &sorted_entries,
        &RackDragState::ArrangementBlock {
            arrangement_id: "arr-a".to_string(),
            from_coordinate: "A1".to_string(),
        },
        "A3",
    )
    .expect("preview");
    assert_eq!(
        preview
            .get("A1")
            .and_then(|cell| cell.predicted_entry.as_ref())
            .map(|entry| entry.arrangement_id.as_str()),
        Some("arr-b")
    );
    assert_eq!(
        preview
            .get("A2")
            .and_then(|cell| cell.predicted_entry.as_ref())
            .map(|entry| entry.arrangement_id.as_str()),
        Some("arr-a")
    );
    assert_eq!(
        preview
            .get("A3")
            .and_then(|cell| cell.predicted_entry.as_ref())
            .map(|entry| entry.arrangement_id.as_str()),
        Some("arr-a")
    );
}

#[test]
fn rack_ghost_preview_moves_multiple_blocks_together() {
    let rack = crate::engine::Rack {
        rack_id: "rack-1".to_string(),
        name: "Bench".to_string(),
        profile: crate::engine::RackProfileSnapshot::custom(1, 4),
        placements: vec![
            crate::engine::RackPlacementEntry {
                coordinate: "A1".to_string(),
                occupant: Some(crate::engine::RackOccupant::Container {
                    container_id: "container-a1".to_string(),
                }),
                arrangement_id: "arr-a".to_string(),
                order_index: 0,
                role_label: "lane_1".to_string(),
            },
            crate::engine::RackPlacementEntry {
                coordinate: "A2".to_string(),
                occupant: Some(crate::engine::RackOccupant::Container {
                    container_id: "container-a2".to_string(),
                }),
                arrangement_id: "arr-a".to_string(),
                order_index: 1,
                role_label: "lane_2".to_string(),
            },
            crate::engine::RackPlacementEntry {
                coordinate: "A3".to_string(),
                occupant: Some(crate::engine::RackOccupant::Container {
                    container_id: "container-b1".to_string(),
                }),
                arrangement_id: "arr-b".to_string(),
                order_index: 0,
                role_label: "lane_1".to_string(),
            },
            crate::engine::RackPlacementEntry {
                coordinate: "A4".to_string(),
                occupant: Some(crate::engine::RackOccupant::Container {
                    container_id: "container-c1".to_string(),
                }),
                arrangement_id: "arr-c".to_string(),
                order_index: 0,
                role_label: "lane_1".to_string(),
            },
        ],
        created_by_op: None,
        created_at_unix_ms: 0,
    };
    let sorted_entries = GENtleApp::rack_sorted_entries(&rack);
    let preview = GENtleApp::rack_ghost_preview_map(
        &rack,
        &sorted_entries,
        &RackDragState::ArrangementBlocks {
            arrangement_ids: vec!["arr-c".to_string(), "arr-a".to_string()],
            from_coordinate: "A1".to_string(),
        },
        "A3",
    )
    .expect("preview");
    assert_eq!(preview.len(), 2);
    assert_eq!(
        preview
            .get("A3")
            .and_then(|cell| cell.predicted_entry.as_ref())
            .map(|entry| entry.arrangement_id.as_str()),
        Some("arr-c")
    );
    assert_eq!(
        preview
            .get("A4")
            .and_then(|cell| cell.predicted_entry.as_ref())
            .map(|entry| entry.arrangement_id.as_str()),
        Some("arr-b")
    );
}

#[test]
fn rack_ghost_preview_moves_multiple_samples_together_within_block() {
    let rack = crate::engine::Rack {
        rack_id: "rack-1".to_string(),
        name: "Bench".to_string(),
        profile: crate::engine::RackProfileSnapshot::custom(1, 4),
        placements: vec![
            crate::engine::RackPlacementEntry {
                coordinate: "A1".to_string(),
                occupant: Some(crate::engine::RackOccupant::Container {
                    container_id: "container-a1".to_string(),
                }),
                arrangement_id: "arr-a".to_string(),
                order_index: 0,
                role_label: "lane_1".to_string(),
            },
            crate::engine::RackPlacementEntry {
                coordinate: "A2".to_string(),
                occupant: Some(crate::engine::RackOccupant::Container {
                    container_id: "container-a2".to_string(),
                }),
                arrangement_id: "arr-a".to_string(),
                order_index: 1,
                role_label: "lane_2".to_string(),
            },
            crate::engine::RackPlacementEntry {
                coordinate: "A3".to_string(),
                occupant: Some(crate::engine::RackOccupant::Container {
                    container_id: "container-a3".to_string(),
                }),
                arrangement_id: "arr-a".to_string(),
                order_index: 2,
                role_label: "lane_3".to_string(),
            },
            crate::engine::RackPlacementEntry {
                coordinate: "A4".to_string(),
                occupant: Some(crate::engine::RackOccupant::Container {
                    container_id: "container-a4".to_string(),
                }),
                arrangement_id: "arr-a".to_string(),
                order_index: 3,
                role_label: "lane_4".to_string(),
            },
        ],
        created_by_op: None,
        created_at_unix_ms: 0,
    };
    let sorted_entries = GENtleApp::rack_sorted_entries(&rack);
    let preview = GENtleApp::rack_ghost_preview_map(
        &rack,
        &sorted_entries,
        &RackDragState::Samples {
            arrangement_id: "arr-a".to_string(),
            from_coordinates: vec!["A1".to_string(), "A3".to_string()],
        },
        "A2",
    )
    .expect("preview");
    assert_eq!(preview.len(), 2);
    assert_eq!(
        preview
            .get("A1")
            .and_then(|cell| cell.predicted_entry.as_ref())
            .map(|entry| entry.role_label.as_str()),
        Some("lane_2")
    );
    assert_eq!(
        preview
            .get("A2")
            .and_then(|cell| cell.predicted_entry.as_ref())
            .map(|entry| entry.role_label.as_str()),
        Some("lane_1")
    );
    assert!(!preview.contains_key("A3"));
}

#[test]
fn rack_autoscroll_offset_for_pointer_moves_toward_edges_and_clamps() {
    let viewport = egui::Rect::from_min_max(egui::pos2(100.0, 100.0), egui::pos2(300.0, 260.0));
    let scrolled = GENtleApp::rack_autoscroll_offset_for_pointer(
        egui::vec2(40.0, 50.0),
        egui::pos2(296.0, 256.0),
        viewport,
        120.0,
        140.0,
    );
    assert!(scrolled.x > 40.0);
    assert!(scrolled.y > 50.0);

    let clamped = GENtleApp::rack_autoscroll_offset_for_pointer(
        egui::vec2(2.0, 3.0),
        egui::pos2(90.0, 90.0),
        viewport,
        120.0,
        140.0,
    );
    assert_eq!(clamped.x, 0.0);
    assert_eq!(clamped.y, 0.0);
}

#[test]
fn rack_coordinate_for_slot_supports_multi_letter_rows() {
    let profile = RackProfileSnapshot::custom(28, 3);
    assert_eq!(GENtleApp::rack_coordinate_for_slot(&profile, 26, 0), "AA1");
    assert_eq!(GENtleApp::rack_coordinate_for_slot(&profile, 27, 2), "AB3");
}

#[test]
fn can_close_project_is_disabled_for_empty_untitled_project() {
    let app = GENtleApp::default();
    assert!(!app.can_close_project());
}

#[test]
fn can_close_project_is_enabled_when_project_path_is_set() {
    let mut app = GENtleApp::default();
    app.current_project_path = Some("/tmp/demo.gentle.json".to_string());
    assert!(app.can_close_project());
}

#[test]
fn default_save_file_names_follow_current_project_name() {
    let mut app = GENtleApp::default();
    app.current_project_path = Some("/tmp/demo.gentle.json".to_string());

    assert_eq!(app.default_project_save_file_name(), "demo.gentle.json");
    assert_eq!(app.default_lineage_svg_file_name(), "demo.lineage.svg");
    assert_eq!(
        app.default_lab_assistant_report_file_name(),
        "demo.lab_assistant_report.odt"
    );
}

#[test]
fn can_close_project_is_enabled_for_unsaved_project_with_content() {
    let mut state = ProjectState::default();
    state.sequences.insert(
        "seq1".to_string(),
        DNAsequence::from_sequence("ACGT").expect("sequence"),
    );
    let mut app = GENtleApp::default();
    app.engine = Arc::new(RwLock::new(GentleEngine::from_state(state)));
    app.current_project_path = None;
    assert!(app.can_close_project());
}

#[test]
fn specialist_window_close_hover_text_mentions_window_and_shortcut() {
    assert_eq!(
        GENtleApp::specialist_window_close_hover_text("Prepare Reference Genome"),
        "Close this Prepare Reference Genome window (Cmd/Ctrl+W)"
    );
}

#[test]
fn configure_platform_viewport_mode_sets_expected_embed_flag() {
    let _lock = crate::genomes::genbank_env_lock()
        .lock()
        .unwrap_or_else(|e| e.into_inner());
    let ctx = egui::Context::default();
    ctx.set_embed_viewports(true);
    let _env_guard = EnvVarGuard::set(super::MACOS_NATIVE_CHILD_VIEWPORTS_ENV, "0");
    let _hosted_env_guard = EnvVarGuard::set(super::MACOS_HOSTED_CHILD_VIEWPORTS_ENV, "0");

    GENtleApp::configure_platform_viewport_mode(&ctx);

    assert!(!ctx.embed_viewports());
}

#[test]
fn macos_hosted_child_viewports_env_override_enables_embed_mode() {
    let _lock = crate::genomes::genbank_env_lock()
        .lock()
        .unwrap_or_else(|e| e.into_inner());
    let ctx = egui::Context::default();
    ctx.set_embed_viewports(false);
    let _native_env_guard = EnvVarGuard::set(super::MACOS_NATIVE_CHILD_VIEWPORTS_ENV, "0");
    let _hosted_env_guard = EnvVarGuard::set(super::MACOS_HOSTED_CHILD_VIEWPORTS_ENV, "1");

    GENtleApp::configure_platform_viewport_mode(&ctx);

    assert_eq!(ctx.embed_viewports(), cfg!(target_os = "macos"));
}

#[test]
fn macos_native_child_viewports_env_override_wins_over_hosted_mode() {
    let _lock = crate::genomes::genbank_env_lock()
        .lock()
        .unwrap_or_else(|e| e.into_inner());
    let ctx = egui::Context::default();
    ctx.set_embed_viewports(true);
    let _env_guard = EnvVarGuard::set(super::MACOS_NATIVE_CHILD_VIEWPORTS_ENV, "1");
    let _hosted_env_guard = EnvVarGuard::set(super::MACOS_HOSTED_CHILD_VIEWPORTS_ENV, "1");

    GENtleApp::configure_platform_viewport_mode(&ctx);

    assert!(!ctx.embed_viewports());
}

fn raw_input_with_root_window_state(fullscreen: bool, maximized: bool) -> egui::RawInput {
    let mut raw_input = egui::RawInput::default();
    let root_viewport = raw_input
        .viewports
        .get_mut(&egui::ViewportId::ROOT)
        .expect("root viewport exists in default egui raw input");
    root_viewport.fullscreen = Some(fullscreen);
    root_viewport.maximized = Some(maximized);
    raw_input
}

#[test]
fn root_viewport_fullscreen_or_maximized_tracks_root_viewport_state() {
    let ctx = egui::Context::default();

    ctx.begin_pass(raw_input_with_root_window_state(true, false));
    assert!(GENtleApp::root_viewport_fullscreen_or_maximized(&ctx));
    let _ = ctx.end_pass();

    ctx.begin_pass(raw_input_with_root_window_state(false, true));
    assert!(GENtleApp::root_viewport_fullscreen_or_maximized(&ctx));
    let _ = ctx.end_pass();

    ctx.begin_pass(raw_input_with_root_window_state(false, false));
    assert!(!GENtleApp::root_viewport_fullscreen_or_maximized(&ctx));
    let _ = ctx.end_pass();
}

#[test]
fn pending_native_child_windows_wait_for_regular_macos_root_viewport() {
    let _lock = crate::genomes::genbank_env_lock()
        .lock()
        .unwrap_or_else(|e| e.into_inner());
    let _native_env_guard = EnvVarGuard::set(super::MACOS_NATIVE_CHILD_VIEWPORTS_ENV, "0");
    let _hosted_env_guard = EnvVarGuard::set(super::MACOS_HOSTED_CHILD_VIEWPORTS_ENV, "0");
    let ctx = egui::Context::default();
    let mut app = GENtleApp::default();
    app.new_windows.push(Window::new_dna(
        DNAsequence::from_sequence("ACGT").expect("sequence"),
        "seq1".to_string(),
        app.engine.clone(),
    ));

    ctx.begin_pass(raw_input_with_root_window_state(true, false));
    app.open_pending_sequence_windows(&ctx);
    let full_output = ctx.end_pass();

    if cfg!(target_os = "macos") {
        assert_eq!(app.new_windows.len(), 1);
        assert!(app.windows.is_empty());
        let root_output = full_output
            .viewport_output
            .get(&egui::ViewportId::ROOT)
            .expect("root viewport output");
        assert!(
            root_output
                .commands
                .contains(&egui::ViewportCommand::Fullscreen(false))
        );
        assert!(
            root_output
                .commands
                .contains(&egui::ViewportCommand::Maximized(false))
        );
    } else {
        assert!(app.new_windows.is_empty());
        assert_eq!(app.windows.len(), 1);
    }
}

#[test]
fn hosted_child_viewport_mode_does_not_defer_for_root_window_state() {
    let _lock = crate::genomes::genbank_env_lock()
        .lock()
        .unwrap_or_else(|e| e.into_inner());
    let _native_env_guard = EnvVarGuard::set(super::MACOS_NATIVE_CHILD_VIEWPORTS_ENV, "0");
    let _hosted_env_guard = EnvVarGuard::set(super::MACOS_HOSTED_CHILD_VIEWPORTS_ENV, "1");
    let ctx = egui::Context::default();

    ctx.begin_pass(raw_input_with_root_window_state(true, false));
    assert!(!GENtleApp::should_defer_native_child_viewport_open_for_root_state(&ctx));
    let _ = ctx.end_pass();
}

#[test]
fn native_child_viewport_mode_accepts_native_window_close_request() {
    let _lock = crate::genomes::genbank_env_lock()
        .lock()
        .unwrap_or_else(|e| e.into_inner());
    let _native_env_guard = EnvVarGuard::set(super::MACOS_NATIVE_CHILD_VIEWPORTS_ENV, "0");
    let _hosted_env_guard = EnvVarGuard::set(super::MACOS_HOSTED_CHILD_VIEWPORTS_ENV, "0");

    assert!(GENtleApp::sequence_window_accepts_native_close_request());
}

#[test]
fn hosted_child_viewport_mode_keeps_explicit_close_handling() {
    let _lock = crate::genomes::genbank_env_lock()
        .lock()
        .unwrap_or_else(|e| e.into_inner());
    let _native_env_guard = EnvVarGuard::set(super::MACOS_NATIVE_CHILD_VIEWPORTS_ENV, "0");
    let _hosted_env_guard = EnvVarGuard::set(super::MACOS_HOSTED_CHILD_VIEWPORTS_ENV, "1");

    assert_eq!(
        GENtleApp::sequence_window_accepts_native_close_request(),
        !cfg!(target_os = "macos")
    );
}

#[test]
fn embedded_sequence_viewport_class_uses_outer_shell_directly() {
    assert!(GENtleApp::sequence_viewport_class_has_embedded_shell(
        egui::ViewportClass::EmbeddedWindow
    ));
    assert!(!GENtleApp::sequence_viewport_class_has_embedded_shell(
        egui::ViewportClass::Immediate
    ));
    assert!(!GENtleApp::sequence_viewport_class_has_embedded_shell(
        egui::ViewportClass::Deferred
    ));
}

#[test]
fn render_main_workspace_host_does_not_spawn_hosted_window_area() {
    let ctx = egui::Context::default();
    let mut app = GENtleApp::default();
    let hosted_layer_id = egui::LayerId::new(
        egui::Order::Middle,
        GENtleApp::main_workspace_hosted_window_id(),
    );

    ctx.begin_pass(egui::RawInput::default());
    crate::egui_compat::show_central_panel_for_test_context(
        &ctx,
        egui::CentralPanel::default(),
        |ui| {
            app.render_main_workspace_host(ui, false);
        },
    );
    assert!(!ctx.memory(|mem| mem.areas().is_visible(&hosted_layer_id)));
    let _ = ctx.end_pass();
}

#[test]
fn hosted_main_workspace_window_renders_as_separate_embedded_window() {
    let ctx = egui::Context::default();
    let mut app = GENtleApp::default();
    let hosted_layer_id = egui::LayerId::new(
        egui::Order::Middle,
        GENtleApp::main_workspace_hosted_window_id(),
    );

    ctx.begin_pass(egui::RawInput::default());
    app.render_hosted_main_workspace_window(&ctx, false);
    assert!(ctx.memory(|mem| mem.areas().is_visible(&hosted_layer_id)));
    let _ = ctx.end_pass();
}

#[test]
fn embedded_window_layer_id_for_root_configuration_and_assistants_use_hosted_window_ids() {
    let app = GENtleApp::default();
    assert_eq!(
        app.embedded_window_layer_id_for_viewport(egui::ViewportId::ROOT),
        Some(egui::LayerId::new(
            egui::Order::Middle,
            GENtleApp::main_workspace_hosted_window_id(),
        ))
    );
    assert_eq!(
        app.embedded_window_layer_id_for_viewport(GENtleApp::configuration_viewport_id()),
        Some(egui::LayerId::new(
            egui::Order::Middle,
            GENtleApp::hosted_configuration_window_id(),
        ))
    );
    assert_eq!(
        app.embedded_window_layer_id_for_viewport(GENtleApp::routine_assistant_viewport_id()),
        Some(egui::LayerId::new(
            egui::Order::Middle,
            GENtleApp::hosted_routine_assistant_window_id(),
        ))
    );
    assert_eq!(
        app.embedded_window_layer_id_for_viewport(GENtleApp::agent_assistant_viewport_id()),
        Some(egui::LayerId::new(
            egui::Order::Middle,
            GENtleApp::hosted_agent_assistant_window_id(),
        ))
    );
}

#[test]
fn embedded_window_layer_id_for_root_tools_uses_hosted_window_ids() {
    let app = GENtleApp::default();
    let cases = [
        (
            GENtleApp::prepare_genome_viewport_id(),
            egui::Id::new((
                "hosted_prepare_genome_window",
                GENtleApp::prepare_genome_viewport_id(),
            )),
        ),
        (
            GENtleApp::retrieve_genome_viewport_id(),
            egui::Id::new((
                "hosted_retrieve_genome_window",
                GENtleApp::retrieve_genome_viewport_id(),
            )),
        ),
        (
            GENtleApp::blast_genome_viewport_id(),
            egui::Id::new((
                "hosted_blast_genome_window",
                GENtleApp::blast_genome_viewport_id(),
            )),
        ),
        (
            GENtleApp::bed_track_viewport_id(),
            egui::Id::new((
                "hosted_bed_track_window",
                GENtleApp::bed_track_viewport_id(),
            )),
        ),
        (
            GENtleApp::gibson_viewport_id(),
            egui::Id::new(("hosted_gibson_window", GENtleApp::gibson_viewport_id())),
        ),
        (
            GENtleApp::arrangement_gel_preview_viewport_id(),
            egui::Id::new((
                "hosted_arrangement_gel_preview_window",
                GENtleApp::arrangement_gel_preview_viewport_id(),
            )),
        ),
        (
            GENtleApp::pcr_design_viewport_id(),
            egui::Id::new((
                "hosted_pcr_design_window",
                GENtleApp::pcr_design_viewport_id(),
            )),
        ),
        (
            GENtleApp::sequencing_confirmation_viewport_id(),
            egui::Id::new((
                "hosted_sequencing_confirmation_window",
                GENtleApp::sequencing_confirmation_viewport_id(),
            )),
        ),
        (
            GENtleApp::planning_viewport_id(),
            egui::Id::new(("hosted_planning_window", GENtleApp::planning_viewport_id())),
        ),
        (
            GENtleApp::uniprot_viewport_id(),
            egui::Id::new(("hosted_uniprot_window", GENtleApp::uniprot_viewport_id())),
        ),
    ];

    for (viewport_id, hosted_id) in cases {
        assert_eq!(
            app.embedded_window_layer_id_for_viewport(viewport_id),
            Some(egui::LayerId::new(egui::Order::Middle, hosted_id))
        );
    }
}

#[test]
fn embedded_root_tool_windows_render_as_sibling_hosted_windows() {
    let screen_rect = egui::Rect::from_min_size(egui::Pos2::ZERO, egui::vec2(1600.0, 1000.0));

    let ctx = egui::Context::default();
    ctx.set_embed_viewports(true);
    let mut app = GENtleApp::default();
    app.show_configuration_dialog = true;
    ctx.begin_pass(egui::RawInput {
        screen_rect: Some(screen_rect),
        ..Default::default()
    });
    app.render_configuration_dialog(&ctx);
    assert!(ctx.memory(|mem| {
        mem.areas().is_visible(&egui::LayerId::new(
            egui::Order::Middle,
            GENtleApp::hosted_configuration_window_id(),
        ))
    }));
    assert!(!ctx.memory(|mem| {
        mem.areas().is_visible(&egui::LayerId::new(
            egui::Order::Middle,
            egui::Id::new(GENtleApp::configuration_viewport_id()),
        ))
    }));
    let _ = ctx.end_pass();

    let ctx = egui::Context::default();
    ctx.set_embed_viewports(true);
    let mut app = GENtleApp::default();
    app.show_routine_assistant_dialog = true;
    ctx.begin_pass(egui::RawInput {
        screen_rect: Some(screen_rect),
        ..Default::default()
    });
    app.render_routine_assistant_dialog(&ctx);
    assert!(ctx.memory(|mem| {
        mem.areas().is_visible(&egui::LayerId::new(
            egui::Order::Middle,
            GENtleApp::hosted_routine_assistant_window_id(),
        ))
    }));
    assert!(!ctx.memory(|mem| {
        mem.areas().is_visible(&egui::LayerId::new(
            egui::Order::Middle,
            egui::Id::new(GENtleApp::routine_assistant_viewport_id()),
        ))
    }));
    let _ = ctx.end_pass();

    let ctx = egui::Context::default();
    ctx.set_embed_viewports(true);
    let mut app = GENtleApp::default();
    app.show_command_palette_dialog = true;
    ctx.begin_pass(egui::RawInput {
        screen_rect: Some(screen_rect),
        ..Default::default()
    });
    app.render_command_palette_dialog(&ctx);
    assert!(ctx.memory(|mem| {
        mem.areas().is_visible(&egui::LayerId::new(
            egui::Order::Middle,
            egui::Id::new("Command Palette"),
        ))
    }));
    assert!(!ctx.memory(|mem| {
        mem.areas().is_visible(&egui::LayerId::new(
            egui::Order::Middle,
            egui::Id::new(GENtleApp::command_palette_viewport_id()),
        ))
    }));
    let _ = ctx.end_pass();

    let ctx = egui::Context::default();
    ctx.set_embed_viewports(true);
    let mut app = GENtleApp::default();
    app.history_ui.show_panel = true;
    ctx.begin_pass(egui::RawInput {
        screen_rect: Some(screen_rect),
        ..Default::default()
    });
    app.render_history_panel(&ctx);
    assert!(ctx.memory(|mem| {
        mem.areas().is_visible(&egui::LayerId::new(
            egui::Order::Middle,
            egui::Id::new("Operation History"),
        ))
    }));
    assert!(!ctx.memory(|mem| {
        mem.areas().is_visible(&egui::LayerId::new(
            egui::Order::Middle,
            egui::Id::new(GENtleApp::history_viewport_id()),
        ))
    }));
    let _ = ctx.end_pass();
}

#[test]
fn embedded_window_layer_id_for_sequence_viewport_uses_hosted_sequence_window_id() {
    let dna = DNAsequence::from_sequence("ACGT").expect("sequence");
    let mut app = GENtleApp::default();
    let viewport_id =
        app.register_window(Window::new_dna(dna, "seq1".to_string(), app.engine.clone()));

    assert_eq!(
        app.embedded_window_layer_id_for_viewport(viewport_id),
        Some(egui::LayerId::new(
            egui::Order::Middle,
            egui::Id::new(("hosted_sequence_window_v2", viewport_id)),
        ))
    );
}

#[test]
fn embedded_sequence_viewport_renders_without_legacy_title_bar_window() {
    let ctx = egui::Context::default();
    ctx.set_embed_viewports(true);
    let dna = DNAsequence::from_sequence("ACGT").expect("sequence");
    let mut app = GENtleApp::default();
    let viewport_id =
        app.register_window(Window::new_dna(dna, "seq1".to_string(), app.engine.clone()));
    let window = app
        .windows
        .get(&viewport_id)
        .cloned()
        .expect("registered sequence window");
    let title = window
        .read()
        .map(|guard| guard.name())
        .expect("window name");
    let hosted_layer_id = egui::LayerId::new(
        egui::Order::Middle,
        egui::Id::new(("hosted_sequence_window_v2", viewport_id)),
    );
    let stale_viewport_layer_id =
        egui::LayerId::new(egui::Order::Middle, egui::Id::new(viewport_id));
    let stale_title_layer_id = GENtleApp::stale_hosted_window_title_layer_id(&title);

    ctx.begin_pass(egui::RawInput::default());
    crate::egui_compat::show_legacy_layer_for_tests(&ctx, stale_title_layer_id, |ui| {
        ui.label("legacy title shell");
    });
    assert!(ctx.memory(|mem| mem.areas().is_visible(&stale_title_layer_id)));
    let _ = ctx.end_pass();

    ctx.begin_pass(egui::RawInput::default());
    app.pending_focus_viewports.clear();
    app.pending_viewport_focus_timestamps.clear();

    let initial_position = app.pending_window_initial_positions.remove(&viewport_id);
    app.show_window(&ctx, viewport_id, window, initial_position);
    assert!(ctx.memory(|mem| mem.areas().is_visible(&hosted_layer_id)));
    assert!(!ctx.memory(|mem| mem.areas().is_visible(&stale_viewport_layer_id)));
    let _ = ctx.end_pass();
    assert!(!ctx.memory(|mem| mem.areas().is_visible(&stale_title_layer_id)));
}

#[test]
fn embedded_sequence_viewport_clamps_stale_full_width_shell() {
    let ctx = egui::Context::default();
    ctx.set_embed_viewports(true);
    let screen_rect = egui::Rect::from_min_size(egui::Pos2::ZERO, egui::vec2(1600.0, 1000.0));
    let dna = DNAsequence::from_sequence("ACGTACGT").expect("sequence");
    let mut app = GENtleApp::default();
    let viewport_id =
        app.register_window(Window::new_dna(dna, "seq1".to_string(), app.engine.clone()));
    let window = app
        .windows
        .get(&viewport_id)
        .cloned()
        .expect("registered sequence window");
    const OUTER_CHROME_TOLERANCE_PX: f32 = 96.0;
    let stable_id = egui::Id::new(("hosted_sequence_window_v2", viewport_id));
    let safe = crate::egui_compat::hosted_window_safe_rect_for_rect(screen_rect);
    let max_size = crate::egui_compat::hosted_window_max_inner_size(
        safe,
        egui::vec2(820.0, 520.0),
        egui::vec2(
            super::EMBEDDED_SEQUENCE_WINDOW_DRAG_MARGIN_X_PX,
            super::EMBEDDED_SEQUENCE_WINDOW_DRAG_MARGIN_Y_PX,
        ),
    );
    assert!(
        max_size.x >= safe.width() - 65.0 && max_size.y >= safe.height() - 65.0,
        "embedded sequence windows should reserve only a small chrome/drag allowance: safe={safe:?}, max_size={max_size:?}"
    );

    ctx.begin_pass(egui::RawInput {
        screen_rect: Some(screen_rect),
        ..Default::default()
    });
    let mut stale_open = true;
    crate::egui_compat::show_hosted_window(
        &ctx,
        &crate::egui_compat::HostedWindowSpec::new(
            "seq1",
            stable_id,
            egui::vec2(5000.0, 5000.0),
            egui::vec2(820.0, 520.0),
        ),
        &mut stale_open,
        |ui| {
            ui.set_min_size(egui::vec2(5000.0, 5000.0));
            ui.label("stale full-width sequence shell");
        },
    );
    let stale_rect = ctx
        .memory(|mem| mem.area_rect(stable_id))
        .expect("stale sequence area should be visible");
    assert!(
        stale_rect.width() > max_size.x + OUTER_CHROME_TOLERANCE_PX
            || stale_rect.height() > max_size.y + OUTER_CHROME_TOLERANCE_PX,
        "stale_rect={stale_rect:?}, max_size={max_size:?}"
    );
    let _ = ctx.end_pass();

    ctx.begin_pass(egui::RawInput {
        screen_rect: Some(screen_rect),
        ..Default::default()
    });
    let initial_position = app.pending_window_initial_positions.remove(&viewport_id);
    app.show_window(&ctx, viewport_id, window, initial_position);

    let rect = ctx
        .memory(|mem| mem.area_rect(stable_id))
        .expect("sequence area should be visible");
    assert!(
        rect.width() <= max_size.x + OUTER_CHROME_TOLERANCE_PX,
        "rect={rect:?}, max_size={max_size:?}"
    );
    assert!(
        rect.height() <= max_size.y + OUTER_CHROME_TOLERANCE_PX,
        "rect={rect:?}, max_size={max_size:?}"
    );
    let _ = ctx.end_pass();
}

#[test]
fn embedded_sequence_stale_clamp_keeps_project_window_area() {
    let ctx = egui::Context::default();
    ctx.set_embed_viewports(true);
    let screen_rect = egui::Rect::from_min_size(egui::Pos2::ZERO, egui::vec2(1600.0, 1000.0));
    let dna = DNAsequence::from_sequence("ACGTACGT").expect("sequence");
    let mut app = GENtleApp::default();
    let viewport_id =
        app.register_window(Window::new_dna(dna, "seq1".to_string(), app.engine.clone()));
    let window = app
        .windows
        .get(&viewport_id)
        .cloned()
        .expect("registered sequence window");
    let sequence_stable_id = egui::Id::new(("hosted_sequence_window_v2", viewport_id));
    let project_window_id = GENtleApp::main_workspace_hosted_window_id();
    let project_spec = crate::egui_compat::HostedWindowSpec::new(
        "Project Test",
        project_window_id,
        egui::vec2(980.0, 640.0),
        egui::vec2(900.0, 560.0),
    )
    .initial_pos(Some(egui::pos2(220.0, 140.0)));

    ctx.begin_pass(egui::RawInput {
        screen_rect: Some(screen_rect),
        ..Default::default()
    });
    let mut stale_open = true;
    crate::egui_compat::show_hosted_window(
        &ctx,
        &crate::egui_compat::HostedWindowSpec::new(
            "seq1",
            sequence_stable_id,
            egui::vec2(5000.0, 5000.0),
            egui::vec2(820.0, 520.0),
        ),
        &mut stale_open,
        |ui| {
            ui.set_min_size(egui::vec2(5000.0, 5000.0));
            ui.label("stale full-width sequence shell");
        },
    );
    let _ = ctx.end_pass();

    ctx.begin_pass(egui::RawInput {
        screen_rect: Some(screen_rect),
        ..Default::default()
    });
    let mut project_open = true;
    crate::egui_compat::show_hosted_window(&ctx, &project_spec, &mut project_open, |ui| {
        ui.label("project content");
    });
    let project_rect_before = ctx
        .memory(|mem| mem.area_rect(project_window_id))
        .expect("project area rect before sequence render");
    let initial_position = app.pending_window_initial_positions.remove(&viewport_id);
    app.show_window(&ctx, viewport_id, window, initial_position);

    let (project_rect_after, sequence_rect) = ctx.memory(|mem| {
        (
            mem.area_rect(project_window_id)
                .expect("project area rect after sequence render"),
            mem.area_rect(sequence_stable_id)
                .expect("sequence area rect"),
        )
    });
    assert_eq!(project_rect_after, project_rect_before);
    assert!(
        project_rect_after.min.x > 0.0 && project_rect_after.min.y > 0.0,
        "project_rect_after={project_rect_after:?}"
    );
    assert!(
        sequence_rect.width() >= 820.0 && sequence_rect.height() >= 520.0,
        "sequence_rect={sequence_rect:?}"
    );
    let _ = ctx.end_pass();
}

#[test]
fn embedded_sequence_window_does_not_auto_expand_to_safe_rect_with_overflowing_content() {
    let ctx = egui::Context::default();
    ctx.set_embed_viewports(true);
    let screen_rect = egui::Rect::from_min_size(egui::Pos2::ZERO, egui::vec2(2000.0, 1200.0));
    let safe = crate::egui_compat::hosted_window_safe_rect_for_rect(screen_rect);
    let sequence = "ACGT".repeat(1_250);
    let dna = DNAsequence::from_sequence(&sequence).expect("sequence");
    let mut app = GENtleApp::default();
    let viewport_id =
        app.register_window(Window::new_dna(dna, "seq1".to_string(), app.engine.clone()));
    let sequence_stable_id = egui::Id::new(("hosted_sequence_window_v2", viewport_id));

    for _ in 0..3 {
        let window = app
            .windows
            .get(&viewport_id)
            .cloned()
            .expect("registered sequence window");
        let initial_position = app.pending_window_initial_positions.remove(&viewport_id);
        ctx.begin_pass(egui::RawInput {
            screen_rect: Some(screen_rect),
            ..Default::default()
        });
        app.show_window(&ctx, viewport_id, window, initial_position);
        let _ = ctx.end_pass();
    }

    let rect = ctx
        .memory(|mem| mem.area_rect(sequence_stable_id))
        .expect("sequence area should be visible");
    assert!(
        rect.width() < safe.width() * 0.9,
        "embedded sequence window should keep its default user-resizable width instead of auto-expanding to the safe viewport: rect={rect:?}, safe={safe:?}"
    );
    assert!(
        rect.height() < safe.height() * 0.9,
        "embedded sequence window should keep its default user-resizable height instead of auto-expanding to the safe viewport: rect={rect:?}, safe={safe:?}"
    );
}

#[test]
fn embedded_sequence_window_can_shrink_after_overflowing_content_render() {
    let ctx = egui::Context::default();
    ctx.set_embed_viewports(true);
    let screen_rect = egui::Rect::from_min_size(egui::Pos2::ZERO, egui::vec2(2000.0, 1200.0));
    let sequence = "ACGT".repeat(1_250);
    let dna = DNAsequence::from_sequence(&sequence).expect("sequence");
    let mut app = GENtleApp::default();
    let viewport_id =
        app.register_window(Window::new_dna(dna, "seq1".to_string(), app.engine.clone()));
    let sequence_stable_id = egui::Id::new(("hosted_sequence_window_v2", viewport_id));

    for _ in 0..2 {
        let window = app
            .windows
            .get(&viewport_id)
            .cloned()
            .expect("registered sequence window");
        let initial_position = app.pending_window_initial_positions.remove(&viewport_id);
        ctx.begin_pass(egui::RawInput {
            screen_rect: Some(screen_rect),
            ..Default::default()
        });
        app.show_window(&ctx, viewport_id, window, initial_position);
        let _ = ctx.end_pass();
    }
    let initial_rect = ctx
        .memory(|mem| mem.area_rect(sequence_stable_id))
        .expect("sequence area should be visible");
    let drag_start = initial_rect.right_bottom() - egui::vec2(4.0, 4.0);

    let window = app
        .windows
        .get(&viewport_id)
        .cloned()
        .expect("registered sequence window");
    ctx.begin_pass(egui::RawInput {
        screen_rect: Some(screen_rect),
        events: vec![
            egui::Event::PointerMoved(drag_start),
            egui::Event::PointerButton {
                pos: drag_start,
                button: egui::PointerButton::Primary,
                pressed: true,
                modifiers: egui::Modifiers::default(),
            },
        ],
        ..Default::default()
    });
    app.show_window(&ctx, viewport_id, window, None);
    let _ = ctx.end_pass();

    let drag_mid = drag_start - egui::vec2(80.0, 60.0);
    let window = app
        .windows
        .get(&viewport_id)
        .cloned()
        .expect("registered sequence window");
    ctx.begin_pass(egui::RawInput {
        screen_rect: Some(screen_rect),
        events: vec![egui::Event::PointerMoved(drag_mid)],
        ..Default::default()
    });
    app.show_window(&ctx, viewport_id, window, None);
    let _ = ctx.end_pass();

    let drag_end = drag_start - egui::vec2(180.0, 120.0);
    let window = app
        .windows
        .get(&viewport_id)
        .cloned()
        .expect("registered sequence window");
    ctx.begin_pass(egui::RawInput {
        screen_rect: Some(screen_rect),
        events: vec![egui::Event::PointerMoved(drag_end)],
        ..Default::default()
    });
    app.show_window(&ctx, viewport_id, window, None);
    let shrunken_rect = ctx
        .memory(|mem| mem.area_rect(sequence_stable_id))
        .expect("sequence area should remain visible");
    assert!(
        shrunken_rect.width() <= initial_rect.width() - 70.0,
        "initial_rect={initial_rect:?}, shrunken_rect={shrunken_rect:?}"
    );
    assert!(
        shrunken_rect.height() <= initial_rect.height() - 50.0,
        "initial_rect={initial_rect:?}, shrunken_rect={shrunken_rect:?}"
    );
    let _ = ctx.end_pass();

    let window = app
        .windows
        .get(&viewport_id)
        .cloned()
        .expect("registered sequence window");
    ctx.begin_pass(egui::RawInput {
        screen_rect: Some(screen_rect),
        events: vec![egui::Event::PointerButton {
            pos: drag_end,
            button: egui::PointerButton::Primary,
            pressed: false,
            modifiers: egui::Modifiers::default(),
        }],
        ..Default::default()
    });
    app.show_window(&ctx, viewport_id, window, None);
    let _ = ctx.end_pass();
}

#[test]
fn newly_focused_embedded_sequence_window_renders_above_hosted_project_window() {
    let ctx = egui::Context::default();
    ctx.set_embed_viewports(true);
    let screen_rect = egui::Rect::from_min_size(egui::Pos2::ZERO, egui::vec2(1600.0, 1000.0));
    let dna = DNAsequence::from_sequence("ACGTACGT").expect("sequence");
    let mut app = GENtleApp::default();
    let viewport_id =
        app.register_window(Window::new_dna(dna, "seq1".to_string(), app.engine.clone()));
    let window = app
        .windows
        .get(&viewport_id)
        .cloned()
        .expect("registered sequence window");
    let project_window_id = GENtleApp::main_workspace_hosted_window_id();
    let sequence_layer_id = egui::LayerId::new(
        egui::Order::Foreground,
        egui::Id::new(("hosted_sequence_window_v2", viewport_id)),
    );

    ctx.begin_pass(egui::RawInput {
        screen_rect: Some(screen_rect),
        ..Default::default()
    });
    app.render_hosted_main_workspace_window(&ctx, false);
    let initial_position = app.pending_window_initial_positions.remove(&viewport_id);
    app.show_window(&ctx, viewport_id, window, initial_position);

    let project_layer_id = ctx
        .memory(|mem| {
            mem.areas()
                .visible_layer_ids()
                .into_iter()
                .find(|layer_id| layer_id.id == project_window_id)
        })
        .expect("project layer should be visible");
    assert!(ctx.memory(|mem| mem.areas().is_visible(&project_layer_id)));
    assert!(ctx.memory(|mem| mem.areas().is_visible(&sequence_layer_id)));
    let (project_rect, sequence_rect) = ctx.memory(|mem| {
        (
            mem.area_rect(project_layer_id.id)
                .expect("project area rect"),
            mem.area_rect(sequence_layer_id.id)
                .expect("sequence area rect"),
        )
    });
    let overlap_min = egui::pos2(
        project_rect.min.x.max(sequence_rect.min.x),
        project_rect.min.y.max(sequence_rect.min.y),
    );
    let overlap_max = egui::pos2(
        project_rect.max.x.min(sequence_rect.max.x),
        project_rect.max.y.min(sequence_rect.max.y),
    );
    assert!(
        overlap_min.x < overlap_max.x && overlap_min.y < overlap_max.y,
        "project_rect={project_rect:?}, sequence_rect={sequence_rect:?}"
    );
    let overlap_center = egui::pos2(
        (overlap_min.x + overlap_max.x) / 2.0,
        (overlap_min.y + overlap_max.y) / 2.0,
    );
    assert_eq!(
        ctx.layer_id_at(overlap_center),
        Some(sequence_layer_id),
        "a newly focused sequence window should cover the project window at overlapping points"
    );
    assert!(!app.viewport_foreground_requested(viewport_id));
    let _ = ctx.end_pass();
}

#[test]
fn closing_sequence_window_with_open_rna_mapping_detaches_auxiliary_host() {
    let dna = DNAsequence::from_sequence("ACGT").expect("sequence");
    let mut app = GENtleApp::default();
    let mut window = Window::new_dna(dna, "seq1".to_string(), app.engine.clone());
    window.seed_rna_read_mapping_window_for_tests("seq1", 17, "TP73");
    let viewport_id = app.register_window(window);

    app.windows_to_close.write().unwrap().push(viewport_id);
    app.process_window_close_queue();

    assert!(
        !app.windows.contains_key(&viewport_id),
        "closing the parent sequence window should remove the visible DNA viewport"
    );
    assert!(
        app.detached_auxiliary_window_hosts
            .contains_key(&viewport_id),
        "an open RNA-read Mapping workspace should keep a detached auxiliary host alive"
    );
    let entries = app.collect_open_window_entries();
    assert!(
        entries
            .iter()
            .any(|entry| { entry.title.contains("RNA-read Mapping - TP73 (seq1)") }),
        "detached auxiliary hosts should keep RNA-read Mapping listed in the Windows menu"
    );
    let mapping_viewport_id =
        egui::ViewportId::from_hash_of(("rna_read_mapping_viewport", "seq1", 17usize));
    assert!(
        app.embedded_window_layer_id_for_viewport(mapping_viewport_id)
            .is_some(),
        "detached auxiliary hosts should still resolve embedded layer ids for focus routing"
    );
}

#[test]
fn focusing_rna_mapping_from_windows_menu_queues_child_focus_without_owner_foreground() {
    let ctx = egui::Context::default();
    ctx.set_embed_viewports(true);
    let dna = DNAsequence::from_sequence("ACGT").expect("sequence");
    let mut app = GENtleApp::default();
    let mut window = Window::new_dna(dna, "seq1".to_string(), app.engine.clone());
    window.seed_rna_read_mapping_window_for_tests("seq1", 17, "TP73");
    let owner_viewport_id = app.register_window(window);
    app.pending_focus_viewports.clear();
    let mapping_viewport_id =
        egui::ViewportId::from_hash_of(("rna_read_mapping_viewport", "seq1", 17usize));

    app.focus_window_viewport(&ctx, mapping_viewport_id);

    assert!(
        !app.pending_focus_viewports.contains(&owner_viewport_id),
        "embedded RNA-read Mapping focus should keep the DNA host visible without queuing it above the child workspace"
    );
    let window = app
        .windows
        .get(&owner_viewport_id)
        .expect("registered sequence owner");
    assert!(
        window
            .read()
            .expect("window")
            .rna_read_mapping_focus_requested_for_tests(),
        "the child workspace should render in foreground order on the next frame"
    );
    assert_eq!(
        app.embedded_window_layer_id_for_viewport(mapping_viewport_id),
        Some(egui::LayerId::new(
            egui::Order::Foreground,
            egui::Id::new("rna_read_mapping_window_embedded_seq1_17"),
        ))
    );
}

#[test]
fn focusing_rna_mapping_from_windows_menu_renders_child_above_sequence_window() {
    let ctx = egui::Context::default();
    ctx.set_embed_viewports(true);
    let screen_rect = egui::Rect::from_min_size(egui::Pos2::ZERO, egui::vec2(1600.0, 1000.0));
    let dna = DNAsequence::from_sequence("ACGTACGT").expect("sequence");
    let mut app = GENtleApp::default();
    let mut window = Window::new_dna(dna, "seq1".to_string(), app.engine.clone());
    window.seed_rna_read_mapping_window_for_tests("seq1", 17, "TP73");
    let owner_viewport_id = app.register_window(window);
    app.pending_focus_viewports.clear();
    app.pending_viewport_focus_timestamps.clear();
    let mapping_viewport_id =
        egui::ViewportId::from_hash_of(("rna_read_mapping_viewport", "seq1", 17usize));

    app.focus_window_viewport(&ctx, mapping_viewport_id);

    let window = app
        .windows
        .get(&owner_viewport_id)
        .cloned()
        .expect("registered sequence owner");
    let initial_position = app
        .pending_window_initial_positions
        .remove(&owner_viewport_id);
    let sequence_layer_id = egui::LayerId::new(
        egui::Order::Middle,
        egui::Id::new(("hosted_sequence_window_v2", owner_viewport_id)),
    );
    let foreground_sequence_layer_id = egui::LayerId::new(
        egui::Order::Foreground,
        egui::Id::new(("hosted_sequence_window_v2", owner_viewport_id)),
    );
    let mapping_layer_id = egui::LayerId::new(
        egui::Order::Foreground,
        egui::Id::new("rna_read_mapping_window_embedded_seq1_17"),
    );
    let stale_viewport_layer_id =
        egui::LayerId::new(egui::Order::Middle, egui::Id::new(mapping_viewport_id));

    ctx.begin_pass(egui::RawInput {
        screen_rect: Some(screen_rect),
        ..Default::default()
    });
    app.show_window(&ctx, owner_viewport_id, window, initial_position);

    assert!(ctx.memory(|mem| mem.areas().is_visible(&sequence_layer_id)));
    assert!(!ctx.memory(|mem| mem.areas().is_visible(&foreground_sequence_layer_id)));
    assert!(ctx.memory(|mem| mem.areas().is_visible(&mapping_layer_id)));
    assert!(!ctx.memory(|mem| mem.areas().is_visible(&stale_viewport_layer_id)));
    assert!(
        ctx.memory(|mem| mem.areas().visible_layer_ids().contains(&mapping_layer_id)),
        "focused RNA-read Mapping should render as a foreground hosted layer above its middle-order DNA host"
    );
    let _ = ctx.end_pass();
}

#[test]
fn focusing_splicing_expert_from_windows_menu_queues_child_focus_without_owner_foreground() {
    let ctx = egui::Context::default();
    ctx.set_embed_viewports(true);
    let dna = DNAsequence::from_sequence("ACGT").expect("sequence");
    let mut app = GENtleApp::default();
    let mut window = Window::new_dna(dna, "seq1".to_string(), app.engine.clone());
    window.seed_splicing_expert_window_for_tests("seq1", 17, "TP73");
    let owner_viewport_id = app.register_window(window);
    app.pending_focus_viewports.clear();
    let splicing_viewport_id =
        egui::ViewportId::from_hash_of(("splicing_expert_viewport", "seq1", 17usize));

    app.focus_window_viewport(&ctx, splicing_viewport_id);

    assert!(
        !app.pending_focus_viewports.contains(&owner_viewport_id),
        "embedded Splicing Expert focus should keep the DNA host visible without queuing it above the child workspace"
    );
    let window = app
        .windows
        .get(&owner_viewport_id)
        .expect("registered sequence owner");
    assert!(
        window
            .read()
            .expect("window")
            .splicing_expert_focus_requested_for_tests(),
        "the Splicing Expert workspace should render in foreground order on the next frame"
    );
    assert_eq!(
        app.embedded_window_layer_id_for_viewport(splicing_viewport_id),
        Some(egui::LayerId::new(
            egui::Order::Foreground,
            egui::Id::new("splicing_expert_window_embedded_seq1_17"),
        ))
    );
}

#[test]
fn focusing_promoter_design_from_windows_menu_renders_child_above_sequence_window() {
    let ctx = egui::Context::default();
    ctx.set_embed_viewports(true);
    let screen_rect = egui::Rect::from_min_size(egui::Pos2::ZERO, egui::vec2(1600.0, 1000.0));
    let dna = DNAsequence::from_sequence("ACGTACGT").expect("sequence");
    let mut app = GENtleApp::default();
    let mut window = Window::new_dna(dna, "seq1".to_string(), app.engine.clone());
    window.seed_variant_followup_window_for_tests("seq1", 17, "TP73");
    let owner_viewport_id = app.register_window(window);
    app.pending_focus_viewports.clear();
    app.pending_viewport_focus_timestamps.clear();
    let promoter_viewport_id =
        egui::ViewportId::from_hash_of(("variant_followup_viewport", "seq1", 17usize));

    app.focus_window_viewport(&ctx, promoter_viewport_id);

    assert!(
        !app.pending_focus_viewports.contains(&owner_viewport_id),
        "embedded Promoter design focus should keep the DNA host visible without queuing it above the child workspace"
    );

    let window = app
        .windows
        .get(&owner_viewport_id)
        .cloned()
        .expect("registered sequence owner");
    let initial_position = app
        .pending_window_initial_positions
        .remove(&owner_viewport_id);
    let sequence_layer_id = egui::LayerId::new(
        egui::Order::Middle,
        egui::Id::new(("hosted_sequence_window_v2", owner_viewport_id)),
    );
    let foreground_sequence_layer_id = egui::LayerId::new(
        egui::Order::Foreground,
        egui::Id::new(("hosted_sequence_window_v2", owner_viewport_id)),
    );
    let promoter_layer_id = egui::LayerId::new(
        egui::Order::Foreground,
        egui::Id::new("variant_followup_window_embedded_seq1_17"),
    );
    let stale_viewport_layer_id =
        egui::LayerId::new(egui::Order::Middle, egui::Id::new(promoter_viewport_id));

    ctx.begin_pass(egui::RawInput {
        screen_rect: Some(screen_rect),
        ..Default::default()
    });
    app.show_window(&ctx, owner_viewport_id, window.clone(), initial_position);

    assert!(ctx.memory(|mem| mem.areas().is_visible(&sequence_layer_id)));
    assert!(!ctx.memory(|mem| mem.areas().is_visible(&foreground_sequence_layer_id)));
    assert!(ctx.memory(|mem| mem.areas().is_visible(&promoter_layer_id)));
    assert!(!ctx.memory(|mem| mem.areas().is_visible(&stale_viewport_layer_id)));
    assert!(
        ctx.memory(|mem| mem.areas().visible_layer_ids().contains(&promoter_layer_id)),
        "focused Promoter design should render as a foreground hosted layer above its middle-order DNA host"
    );
    let _ = ctx.end_pass();

    let window = app
        .windows
        .get(&owner_viewport_id)
        .expect("registered sequence owner after render");
    assert!(
        !window
            .read()
            .expect("window")
            .variant_followup_focus_requested_for_tests(),
        "Promoter design foreground requests should be consumed after the first render"
    );
}

#[test]
fn detached_rna_mapping_host_is_removed_when_visible_window_owns_same_workspace() {
    let ctx = egui::Context::default();
    let dna = DNAsequence::from_sequence("ACGT").expect("sequence");
    let mut app = GENtleApp::default();

    let mut detached = Window::new_dna(dna.clone(), "seq1".to_string(), app.engine.clone());
    detached.seed_rna_read_mapping_window_for_tests("seq1", 17, "TP73");
    let detached_viewport_id = app.register_window(detached);
    let detached_host = app
        .windows
        .remove(&detached_viewport_id)
        .expect("registered detached host");
    app.detached_auxiliary_window_hosts
        .insert(detached_viewport_id, detached_host);

    let mut visible = Window::new_dna(dna, "seq1".to_string(), app.engine.clone());
    visible.seed_rna_read_mapping_window_for_tests("seq1", 17, "TP73");
    let visible_viewport_id = app.register_window(visible);

    app.render_detached_auxiliary_window_hosts(&ctx);

    assert!(
        app.windows.contains_key(&visible_viewport_id),
        "visible sequence window should remain registered"
    );
    assert!(
        !app.detached_auxiliary_window_hosts
            .contains_key(&detached_viewport_id),
        "detached host should be discarded once a visible sequence window renders the same RNA-read Mapping workspace"
    );
}

#[test]
fn stale_rna_mapping_title_area_in_root_context_is_reset_when_detected() {
    let ctx = egui::Context::default();
    let dna = DNAsequence::from_sequence("ACGT").expect("sequence");
    let mut app = GENtleApp::default();

    let mut visible = Window::new_dna(dna, "seq1".to_string(), app.engine.clone());
    visible.seed_rna_read_mapping_window_for_tests("seq1", 17, "TP73");
    app.register_window(visible);

    let title = "RNA-read Mapping - TP73 (seq1)";
    let stale_title_layer_id = GENtleApp::stale_hosted_window_title_layer_id(title);

    ctx.begin_pass(egui::RawInput::default());
    crate::egui_compat::show_legacy_layer_for_tests(&ctx, stale_title_layer_id, |ui| {
        ui.label("legacy root-hosted RNA-read Mapping title shell");
    });
    assert!(ctx.memory(|mem| mem.areas().is_visible(&stale_title_layer_id)));

    assert!(app.reset_root_auxiliary_areas_if_legacy_title_layers_visible(&ctx));
    assert!(!ctx.memory(|mem| mem.areas().is_visible(&stale_title_layer_id)));
    let _ = ctx.end_pass();
}

#[test]
fn stale_splicing_expert_title_area_in_root_context_is_reset_when_detected() {
    let ctx = egui::Context::default();
    let dna = DNAsequence::from_sequence("ACGT").expect("sequence");
    let mut app = GENtleApp::default();

    let mut visible = Window::new_dna(dna, "seq1".to_string(), app.engine.clone());
    visible.seed_splicing_expert_window_for_tests("seq1", 17, "TP73");
    app.register_window(visible);

    let title = "Splicing Expert - TP73 (seq1)";
    let stale_title_layer_id = GENtleApp::stale_hosted_window_title_layer_id(title);

    ctx.begin_pass(egui::RawInput::default());
    crate::egui_compat::show_legacy_layer_for_tests(&ctx, stale_title_layer_id, |ui| {
        ui.label("legacy root-hosted Splicing Expert title shell");
    });
    assert!(ctx.memory(|mem| mem.areas().is_visible(&stale_title_layer_id)));

    assert!(app.reset_root_auxiliary_areas_if_legacy_title_layers_visible(&ctx));
    assert!(!ctx.memory(|mem| mem.areas().is_visible(&stale_title_layer_id)));
    let _ = ctx.end_pass();
}

#[test]
fn embedded_help_viewport_renders_without_nested_help_window_area() {
    let ctx = egui::Context::default();
    ctx.set_embed_viewports(true);
    let mut app = GENtleApp::default();
    app.show_help_dialog = true;
    app.mark_viewport_open_requested(GENtleApp::help_viewport_id());
    let hosted_help_layer_id =
        egui::LayerId::new(egui::Order::Middle, GENtleApp::hosted_help_window_id());
    let stale_title_layer_id = GENtleApp::stale_help_title_layer_id("Help - GUI Manual");
    let stale_viewport_layer_id = egui::LayerId::new(
        egui::Order::Middle,
        egui::Id::new(GENtleApp::help_viewport_id()),
    );

    ctx.begin_pass(egui::RawInput::default());
    app.render_help_dialog(&ctx);
    assert!(ctx.memory(|mem| mem.areas().is_visible(&hosted_help_layer_id)));
    assert!(!ctx.memory(|mem| mem.areas().is_visible(&stale_title_layer_id)));
    assert!(!ctx.memory(|mem| mem.areas().is_visible(&stale_viewport_layer_id)));
    let _ = ctx.end_pass();
}

#[test]
fn embedded_help_tutorial_viewport_renders_without_second_title_bar_window() {
    let ctx = egui::Context::default();
    ctx.set_embed_viewports(true);
    let mut app = GENtleApp::default();
    app.show_help_dialog = true;
    app.help_doc = HelpDoc::Tutorial;
    app.help_tutorial_title = "Load pGEX and digest with BamHI/EcoRI".to_string();
    app.help_tutorial_markdown = "# Load pGEX and digest with BamHI/EcoRI\n\nTutorial.".to_string();
    app.mark_viewport_open_requested(GENtleApp::help_viewport_id());
    let hosted_help_layer_id =
        egui::LayerId::new(egui::Order::Middle, GENtleApp::hosted_help_window_id());
    let stale_title_layer_id =
        GENtleApp::stale_help_title_layer_id("Help - Load pGEX and digest with BamHI/EcoRI");
    let stale_viewport_layer_id = egui::LayerId::new(
        egui::Order::Middle,
        egui::Id::new(GENtleApp::help_viewport_id()),
    );

    ctx.begin_pass(egui::RawInput::default());
    app.render_help_dialog(&ctx);
    assert!(ctx.memory(|mem| mem.areas().is_visible(&hosted_help_layer_id)));
    assert!(!ctx.memory(|mem| mem.areas().is_visible(&stale_title_layer_id)));
    assert!(!ctx.memory(|mem| mem.areas().is_visible(&stale_viewport_layer_id)));
    let _ = ctx.end_pass();
}

#[test]
fn focused_embedded_help_window_renders_in_foreground_without_viewport_title_shell() {
    let ctx = egui::Context::default();
    ctx.set_embed_viewports(true);
    let mut app = GENtleApp::default();
    app.show_help_dialog = true;
    app.queue_focus_viewport(GENtleApp::help_viewport_id());
    let foreground_help_layer_id =
        egui::LayerId::new(egui::Order::Foreground, GENtleApp::hosted_help_window_id());
    let stale_viewport_layer_id = egui::LayerId::new(
        egui::Order::Middle,
        egui::Id::new(GENtleApp::help_viewport_id()),
    );

    ctx.begin_pass(egui::RawInput::default());
    app.render_help_dialog(&ctx);
    assert!(ctx.memory(|mem| mem.areas().is_visible(&foreground_help_layer_id)));
    assert!(!ctx.memory(|mem| mem.areas().is_visible(&stale_viewport_layer_id)));
    assert!(!app.viewport_foreground_requested(GENtleApp::help_viewport_id()));
    let _ = ctx.end_pass();
}

#[test]
fn embedded_agent_assistant_renders_as_single_hosted_window_without_viewport_title_shell() {
    let ctx = egui::Context::default();
    ctx.set_embed_viewports(true);
    let mut app = GENtleApp::default();
    app.show_agent_assistant_dialog = true;
    app.queue_focus_viewport(GENtleApp::agent_assistant_viewport_id());
    let hosted_layer_id = egui::LayerId::new(
        egui::Order::Foreground,
        GENtleApp::hosted_agent_assistant_window_id(),
    );
    let stale_viewport_layer_id = egui::LayerId::new(
        egui::Order::Middle,
        egui::Id::new(GENtleApp::agent_assistant_viewport_id()),
    );

    ctx.begin_pass(egui::RawInput::default());
    app.render_agent_assistant_dialog(&ctx);
    assert!(ctx.memory(|mem| mem.areas().is_visible(&hosted_layer_id)));
    assert!(!ctx.memory(|mem| mem.areas().is_visible(&stale_viewport_layer_id)));
    assert!(!app.viewport_foreground_requested(GENtleApp::agent_assistant_viewport_id()));
    let _ = ctx.end_pass();
}

#[test]
fn agent_assistant_content_scrolls_on_small_viewport() {
    let ctx = egui::Context::default();
    let screen_rect = egui::Rect::from_min_size(egui::Pos2::ZERO, egui::vec2(640.0, 260.0));
    let mut app = GENtleApp::default();
    app.agent_system_id = "builtin_echo".to_string();
    app.agent_prompt = (0..10)
        .map(|idx| format!("agent prompt regression line {idx}"))
        .collect::<Vec<_>>()
        .join("\n");

    ctx.begin_pass(egui::RawInput {
        screen_rect: Some(screen_rect),
        ..Default::default()
    });
    let panel_response = crate::egui_compat::show_central_panel_for_test_context(
        &ctx,
        egui::CentralPanel::default(),
        |ui| {
            app.render_agent_assistant_contents_scrollable(
                ui,
                "agent_assistant_scroll_regression_test",
            )
        },
    );
    let scroll_output = panel_response.inner;

    assert!(
        scroll_output.content_size.y > scroll_output.inner_rect.height() + 80.0,
        "agent assistant content should overflow vertically on a small viewport: content={:?}, inner={:?}",
        scroll_output.content_size,
        scroll_output.inner_rect
    );
    assert!(
        scroll_output.inner_rect.height() <= screen_rect.height(),
        "agent assistant scroll area should stay bounded by the viewport: inner={:?}, screen={:?}",
        scroll_output.inner_rect,
        screen_rect
    );
    let _ = ctx.end_pass();
}

#[test]
fn embedded_configuration_graphics_window_width_stays_bounded_across_frames() {
    let ctx = egui::Context::default();
    ctx.set_embed_viewports(true);
    let mut app = GENtleApp::default();
    app.show_configuration_dialog = true;
    app.configuration_tab = ConfigurationTab::Graphics;
    app.mark_viewport_open_requested(GENtleApp::configuration_viewport_id());
    let hosted_window_id = GENtleApp::hosted_configuration_window_id();
    let screen_rect = egui::Rect::from_min_size(egui::Pos2::ZERO, egui::vec2(1200.0, 900.0));
    let mut widths = Vec::new();

    for _ in 0..4 {
        ctx.begin_pass(egui::RawInput {
            screen_rect: Some(screen_rect),
            ..Default::default()
        });
        app.render_configuration_dialog(&ctx);
        let width = ctx.memory(|mem| {
            mem.area_rect(hosted_window_id)
                .map(|rect| rect.width())
                .unwrap_or_default()
        });
        widths.push(width);
        let _ = ctx.end_pass();
    }

    let max_width = widths.iter().copied().fold(0.0, f32::max);
    let min_width = widths.iter().copied().fold(f32::INFINITY, f32::min);
    assert!(widths.iter().all(|width| *width > 0.0), "widths={widths:?}");
    assert!(max_width <= 820.0, "widths={widths:?}");
    assert!(max_width - min_width <= 8.0, "widths={widths:?}");
}

#[test]
fn prepare_dialog_scroll_area_keeps_long_checklist_reachable_on_small_viewport() {
    let td = tempdir().unwrap();
    let (catalog_path, cache_dir) = write_toy_prepare_catalog(td.path());
    let ctx = egui::Context::default();
    let screen_rect = egui::Rect::from_min_size(egui::Pos2::ZERO, egui::vec2(640.0, 220.0));
    let mut app = GENtleApp::default();
    app.genome_catalog_path = catalog_path;
    app.genome_cache_dir = cache_dir;
    app.genome_catalog_genomes = vec!["ToyGenome".to_string()];
    app.genome_id = "ToyGenome".to_string();
    app.reset_prepare_step_state_from_plan(Some(make_prepare_plan(&[
        PrepareGenomeStepId::ResetIndexes,
        PrepareGenomeStepId::Sequence,
        PrepareGenomeStepId::Annotation,
        PrepareGenomeStepId::FastaIndex,
        PrepareGenomeStepId::GeneIndex,
        PrepareGenomeStepId::BlastIndex,
    ])));
    for (idx, step) in app.genome_prepare_steps.iter_mut().enumerate() {
        step.status = if idx < 2 {
            PrepareGenomeUiStepStatus::Completed
        } else if idx == 2 {
            PrepareGenomeUiStepStatus::Running
        } else {
            PrepareGenomeUiStepStatus::Pending
        };
        step.detail = format!(
            "regression detail line {idx}: keep this prepare checklist reachable after resize"
        );
    }
    app.genome_prepare_status = (0..18)
        .map(|idx| format!("diagnostic status line {idx}: index/progress context"))
        .collect::<Vec<_>>()
        .join("\n");

    ctx.begin_pass(egui::RawInput {
        screen_rect: Some(screen_rect),
        ..Default::default()
    });
    let panel_response = crate::egui_compat::show_central_panel_for_test_context(
        &ctx,
        egui::CentralPanel::default(),
        |ui| {
            app.render_reference_genome_prepare_scroll_area(
                ui,
                "prepare_genome_scroll_regression_test",
            )
        },
    );
    let scroll_output = panel_response.inner;

    assert!(
        scroll_output.content_size.y > scroll_output.inner_rect.height() + 80.0,
        "prepare dialog content should overflow vertically on a small viewport: content={:?}, inner={:?}",
        scroll_output.content_size,
        scroll_output.inner_rect
    );
    assert!(
        scroll_output.inner_rect.height() <= screen_rect.height(),
        "scroll area should stay bounded by the viewport: inner={:?}, screen={:?}",
        scroll_output.inner_rect,
        screen_rect
    );
    let _ = ctx.end_pass();
}

#[test]
fn stale_help_window_area_in_root_context_is_reset_when_detected() {
    let ctx = egui::Context::default();
    let title = "Help - Gibson Arrangements Tutorial";
    let stale_help_layer_ids = GENtleApp::legacy_root_help_layer_ids(title);

    ctx.begin_pass(egui::RawInput::default());
    for layer_id in stale_help_layer_ids.iter().copied() {
        crate::egui_compat::show_legacy_layer_for_tests(&ctx, layer_id, |ui| {
            ui.label("legacy nested help shell");
        });
    }
    assert!(
        stale_help_layer_ids
            .iter()
            .any(|layer_id| ctx.memory(|mem| mem.areas().is_visible(layer_id)))
    );

    assert!(GENtleApp::reset_root_help_areas_if_legacy_layers_visible(
        &ctx, title
    ));
    for layer_id in stale_help_layer_ids {
        assert!(!ctx.memory(|mem| mem.areas().is_visible(&layer_id)));
    }
    let _ = ctx.end_pass();
}

#[test]
fn root_help_cleanup_ignores_stable_hosted_help_window_id() {
    let ctx = egui::Context::default();
    let title = "Help - Gibson Arrangements Tutorial";
    let hosted_help_layer_id =
        egui::LayerId::new(egui::Order::Middle, GENtleApp::hosted_help_window_id());

    ctx.begin_pass(egui::RawInput::default());
    egui::Window::new("Help")
        .id(GENtleApp::hosted_help_window_id())
        .show(&ctx, |ui| {
            ui.label("stable hosted help window");
        });
    assert!(ctx.memory(|mem| mem.areas().is_visible(&hosted_help_layer_id)));

    assert!(!GENtleApp::reset_root_help_areas_if_legacy_layers_visible(
        &ctx, title
    ));
    assert!(ctx.memory(|mem| mem.areas().is_visible(&hosted_help_layer_id)));
    let _ = ctx.end_pass();
}

#[test]
fn request_quit_without_unsaved_content_marks_application_for_close() {
    let mut app = GENtleApp::default();
    app.request_project_action(ProjectAction::Quit);
    assert!(app.pending_app_quit);
    assert!(app.pending_project_action.is_none());
}

#[test]
fn request_quit_with_unsaved_content_prompts_before_closing() {
    let mut state = ProjectState::default();
    state.sequences.insert(
        "seq1".to_string(),
        DNAsequence::from_sequence("ACGT").expect("sequence"),
    );
    let mut app = GENtleApp::default();
    app.engine = Arc::new(RwLock::new(GentleEngine::from_state(state)));

    app.request_project_action(ProjectAction::Quit);

    assert!(matches!(
        app.pending_project_action,
        Some(ProjectAction::Quit)
    ));
    assert!(!app.pending_app_quit);
}

#[test]
fn unsaved_changes_dialog_renders_in_foreground_above_hosted_workspace() {
    let ctx = egui::Context::default();
    ctx.set_embed_viewports(true);
    let mut app = GENtleApp::default();
    app.pending_project_action = Some(ProjectAction::Quit);
    let project_layer_id = egui::LayerId::new(
        egui::Order::Middle,
        GENtleApp::main_workspace_hosted_window_id(),
    );
    let unsaved_layer_id = egui::LayerId::new(
        egui::Order::Foreground,
        GENtleApp::unsaved_changes_dialog_id(),
    );

    ctx.begin_pass(egui::RawInput::default());
    app.render_hosted_main_workspace_window(&ctx, true);
    app.render_unsaved_changes_dialog(&ctx);
    assert!(ctx.memory(|mem| mem.areas().is_visible(&project_layer_id)));
    assert!(ctx.memory(|mem| mem.areas().is_visible(&unsaved_layer_id)));
    let _ = ctx.end_pass();
}

#[test]
fn tutorial_project_catalog_loads_chapters() {
    let entries = GENtleApp::load_tutorial_project_entries().expect("tutorial catalog should load");
    assert!(
        !entries.is_empty(),
        "Expected tutorial project entries to be available"
    );
    assert!(
        entries
            .iter()
            .any(|entry| entry.chapter_id == "tp63_anchor_extension_online"),
        "Expected TP63 tutorial chapter to be present"
    );
    let simple_pcr = entries
        .iter()
        .find(|entry| entry.chapter_id == "simple_pcr_selection_gui")
        .expect("simple PCR tutorial chapter should be present");
    assert_eq!(
        simple_pcr.group_label.as_deref(),
        Some("Primers, PCR & qPCR")
    );
    assert_eq!(simple_pcr.decimal_id.as_deref(), Some("04.01"));
}

#[test]
fn tutorial_help_entries_include_agent_interfaces_from_catalog() {
    let entries = GENtleApp::discover_help_tutorial_entries();
    assert!(
        entries.iter().any(|entry| entry.title
            == "GENtle Agent Assistant and Agent Interfaces Tutorial"
            && entry.group_label.as_deref() == Some("Getting Started & Interfaces")
            && entry.decimal_id.as_deref() == Some("01.01")),
        "Expected curated tutorial help entries to include the agent interfaces tutorial"
    );
}

#[test]
fn tutorial_project_guided_walkthroughs_include_agent_interfaces() {
    let entries = GENtleApp::tutorial_project_guided_walkthrough_entries();
    assert!(
        entries.iter().any(|entry| {
            entry.title == "GENtle Agent Assistant and Agent Interfaces Tutorial"
                && entry
                    .path
                    .ends_with("docs/tutorial/01-01_agent_interfaces.md")
        }),
        "Expected Open Tutorial Project guided walkthroughs to surface the agent interfaces tutorial"
    );
}

#[test]
fn tutorial_project_guided_walkthroughs_follow_catalog_reference_entries() {
    let entries = GENtleApp::tutorial_project_guided_walkthrough_entries();
    assert!(
        entries.iter().any(|entry| {
            entry.title == "GENtle Tutorial Landscape Overview"
                && entry.path.ends_with("docs/tutorial/landscape_overview.md")
        }),
        "Expected Open Tutorial Project guided walkthroughs to include tutorial landscape overview from the catalog"
    );
    assert!(
        entries
            .iter()
            .all(|entry| entry.summary.contains("status: manual/reference")),
        "Expected guided walkthrough entries to stay limited to documentation-only catalog entries"
    );
}

#[test]
fn agent_interfaces_tutorial_is_pinned_when_catalog_fallback_is_used() {
    let mut entries = vec![];
    GENtleApp::ensure_agent_interfaces_tutorial_entry(&mut entries);
    GENtleApp::ensure_agent_interfaces_tutorial_entry(&mut entries);

    let matches = entries
        .iter()
        .filter(|entry| {
            entry
                .path
                .ends_with("docs/tutorial/01-01_agent_interfaces.md")
        })
        .collect::<Vec<_>>();
    assert_eq!(matches.len(), 1);
    assert_eq!(
        matches[0].title,
        "GENtle Agent Assistant and Agent Interfaces Tutorial"
    );
    assert!(matches[0].summary.contains("operational_reference"));
}

#[test]
fn open_tutorial_project_chapter_opens_linked_guide() {
    let mut app = GENtleApp::default();

    app.open_tutorial_project_chapter("gibson_specialist_testing_baseline");
    wait_for_tutorial_project_task(&mut app);

    assert!(app.app_status.contains("Opened tutorial project"));
    assert!(app.show_help_dialog);
    assert_eq!(app.help_doc, HelpDoc::Tutorial);
    assert!(
        app.help_tutorial_title
            .contains("Gibson Specialist Testing Tutorial")
    );
    assert!(app.help_tutorial_markdown.contains("Patterns -> Gibson..."));
}

#[test]
fn open_tutorial_project_chapter_arrangements_baseline_opens_arrangement_ready_state() {
    let mut app = GENtleApp::default();

    app.open_tutorial_project_chapter("gibson_arrangements_baseline");
    wait_for_tutorial_project_task(&mut app);

    assert!(app.app_status.contains("Opened tutorial project"));
    assert!(app.show_help_dialog);
    assert_eq!(app.help_doc, HelpDoc::Tutorial);
    assert!(
        app.help_tutorial_title
            .contains("Gibson Arrangements Tutorial")
    );
    assert!(app.help_tutorial_markdown.contains("serial arrangement"));
    let engine = app.engine.read().unwrap();
    assert!(
        engine
            .state()
            .sequences
            .contains_key("gibson_destination_pgex_with_gibson_insert_demo")
    );
    assert_eq!(engine.state().container_state.arrangements.len(), 1);
}

#[test]
fn tutorial_project_progress_updates_status_and_task_title() {
    let (tx, rx2) = mpsc::channel::<TutorialProjectTaskMessage>();
    let mut app = GENtleApp::default();
    app.tutorial_project_task = Some(TutorialProjectTask {
        job_id: 7,
        started: Instant::now(),
        cancel_requested: Arc::new(AtomicBool::new(false)),
        chapter_id: "tp73_old".to_string(),
        chapter_title: "tp73_old".to_string(),
        receiver: rx2,
    });
    tx.send(TutorialProjectTaskMessage::Progress {
        job_id: 7,
        progress: TutorialProjectTaskProgress {
            chapter_id: "tp73_audit".to_string(),
            chapter_title: "TP73 Audit Tutorial".to_string(),
            phase: "execute_workflow".to_string(),
            item: "Waiting for network response".to_string(),
            percent: Some(0.42),
        },
    })
    .expect("send progress");

    app.poll_tutorial_project_task(&egui::Context::default());

    assert_eq!(app.tutorial_project_progress_fraction, Some(0.42));
    assert_eq!(
        app.tutorial_project_progress_label,
        "Waiting for network response"
    );
    assert!(app.tutorial_project_status.contains("TP73 Audit Tutorial"));
    let task = app
        .tutorial_project_task
        .as_ref()
        .expect("tutorial task remains active");
    assert_eq!(task.chapter_id, "tp73_audit");
    assert_eq!(task.chapter_title, "TP73 Audit Tutorial");
}

#[test]
fn tutorial_project_workflow_progress_is_capped_before_opening_phase() {
    assert_eq!(
        GENtleApp::tutorial_project_scale_workflow_percent(Some(0.0)),
        Some(0.15)
    );
    assert_eq!(
        GENtleApp::tutorial_project_scale_workflow_percent(Some(1.0)),
        Some(0.9)
    );
}

#[test]
fn request_tutorial_project_task_cancel_is_idempotent() {
    let (_tx, rx) = mpsc::channel::<TutorialProjectTaskMessage>();
    let cancel_requested = Arc::new(AtomicBool::new(false));
    let mut app = GENtleApp::default();
    app.tutorial_project_task = Some(TutorialProjectTask {
        job_id: 19,
        started: Instant::now(),
        cancel_requested: cancel_requested.clone(),
        chapter_id: "tp73_audit".to_string(),
        chapter_title: "TP73 Audit Tutorial".to_string(),
        receiver: rx,
    });

    app.request_tutorial_project_task_cancel("test");
    assert!(cancel_requested.load(Ordering::Relaxed));
    assert!(
        app.tutorial_project_status
            .contains("Cancellation requested")
    );
    let events_after_first = app.job_event_log.len();

    app.request_tutorial_project_task_cancel("test");
    assert_eq!(app.job_event_log.len(), events_after_first);
    assert!(
        app.tutorial_project_status
            .contains("already requested for tutorial")
    );
}

#[test]
fn tutorial_project_completion_honors_cancel_without_opening_project() {
    let temp = tempdir().expect("tempdir");
    let project_path = temp.path().join("tutorial.project.gentle.json");
    let guide_path = temp.path().join("tutorial.md");
    let mut state = ProjectState::default();
    state.sequences.insert(
        "seq1".to_string(),
        DNAsequence::from_sequence("ACGT").expect("sequence"),
    );
    state
        .save_to_path(project_path.to_string_lossy().as_ref())
        .expect("save project");
    fs::write(&guide_path, "# Tutorial\n\nBody.\n").expect("write guide");

    let (tx, rx) = mpsc::channel::<TutorialProjectTaskMessage>();
    let cancel_requested = Arc::new(AtomicBool::new(true));
    let mut app = GENtleApp::default();
    app.tutorial_project_task = Some(TutorialProjectTask {
        job_id: 31,
        started: Instant::now(),
        cancel_requested: cancel_requested.clone(),
        chapter_id: "tp73_audit".to_string(),
        chapter_title: "TP73 Audit Tutorial".to_string(),
        receiver: rx,
    });
    tx.send(TutorialProjectTaskMessage::Done {
        job_id: 31,
        result: Ok(TutorialProjectOpenOutcome {
            chapter_id: "tp73_audit".to_string(),
            chapter_title: "TP73 Audit Tutorial".to_string(),
            project_path: project_path.to_string_lossy().to_string(),
            guide_path: guide_path.to_string_lossy().to_string(),
            guide_summary: "tutorial guide".to_string(),
        }),
    })
    .expect("send completion");

    app.poll_tutorial_project_task(&egui::Context::default());

    assert!(app.current_project_path.is_none());
    assert!(!app.show_help_dialog);
    assert!(app.tutorial_project_status.contains("cancelled"));
}

#[test]
fn open_new_windows_from_files_imports_selected_files_consecutively() {
    let mut app = GENtleApp::default();
    let temp = tempdir().expect("tempdir");
    let first_path = temp.path().join("first.fa");
    let second_path = temp.path().join("second.gb");
    fs::copy("test_files/pGEX_3X.fa", &first_path).expect("copy FASTA fixture");
    fs::copy("test_files/tp73.ncbi.gb", &second_path).expect("copy GenBank fixture");

    app.open_new_windows_from_files(&[first_path, second_path]);

    let queued_seq_ids = app
        .new_windows
        .iter()
        .map(|window| window.sequence_id().expect("sequence window"))
        .collect::<Vec<_>>();
    assert_eq!(
        queued_seq_ids,
        vec!["first".to_string(), "second".to_string()]
    );
    let engine = app.engine.read().unwrap();
    assert!(engine.state().sequences.contains_key("first"));
    assert!(engine.state().sequences.contains_key("second"));
    assert_eq!(engine.operation_log().len(), 2);
    assert!(
        app.app_status.contains("loaded 2 sequences"),
        "status was: {}",
        app.app_status
    );
}

#[test]
fn open_new_windows_from_files_keeps_successes_when_later_import_fails() {
    let mut app = GENtleApp::default();
    let temp = tempdir().expect("tempdir");
    let first_path = temp.path().join("first.fa");
    let missing_path = temp.path().join("missing.gb");
    fs::copy("test_files/pGEX_3X.fa", &first_path).expect("copy FASTA fixture");

    app.open_new_windows_from_files(&[first_path, missing_path]);

    let queued_seq_ids = app
        .new_windows
        .iter()
        .map(|window| window.sequence_id().expect("sequence window"))
        .collect::<Vec<_>>();
    assert_eq!(queued_seq_ids, vec!["first".to_string()]);
    let engine = app.engine.read().unwrap();
    assert!(engine.state().sequences.contains_key("first"));
    assert_eq!(engine.operation_log().len(), 1);
    assert!(
        app.app_status.contains("1 import failed"),
        "status was: {}",
        app.app_status
    );
}

#[test]
fn open_new_window_from_missing_file_sets_error_status() {
    let mut app = GENtleApp::default();
    let temp = tempdir().expect("tempdir");
    let missing = temp.path().join("missing_sequence_file.gb");

    app.open_new_window_from_file(missing.to_string_lossy().as_ref());

    assert!(
        app.app_status.contains("Open sequence failed"),
        "status was: {}",
        app.app_status
    );
    assert!(app.new_windows.is_empty());
}

#[test]
fn sort_uniprot_entries_recent_prefers_newest_then_entry_id() {
    let rows = vec![
        UniprotEntrySummary {
            entry_id: "P_B".to_string(),
            imported_at_unix_ms: 1_000,
            ..UniprotEntrySummary::default()
        },
        UniprotEntrySummary {
            entry_id: "P_C".to_string(),
            imported_at_unix_ms: 2_000,
            ..UniprotEntrySummary::default()
        },
        UniprotEntrySummary {
            entry_id: "P_A".to_string(),
            imported_at_unix_ms: 2_000,
            ..UniprotEntrySummary::default()
        },
    ];
    let sorted = GENtleApp::sort_uniprot_entries_recent(rows);
    let ordered_ids = sorted
        .into_iter()
        .map(|row| row.entry_id)
        .collect::<Vec<_>>();
    assert_eq!(ordered_ids, vec!["P_A", "P_C", "P_B"]);
}

#[test]
fn uniprot_dialog_import_swiss_prot_persists_entry() {
    let mut app = GENtleApp::default();
    let temp = tempdir().expect("tempdir");
    let swiss_path = temp.path().join("toy_uniprot_dialog.txt");
    let swiss_text = r#"ID   TOY1_HUMAN              Reviewed;         30 AA.
AC   PTEST1;
DE   RecName: Full=Toy DNA-binding protein;
GN   Name=TOY1;
OS   Homo sapiens (Human).
DR   Ensembl; TX1; ENSPTOY1; ENSGTOY1.
FT   DOMAIN          2..8
FT                   /note="toy domain"
SQ   SEQUENCE   30 AA;  3333 MW;  0000000000000000 CRC64;
     MEEPQSDPSV EPPLSQETFSDLWKLLPEN
//
"#;
    fs::write(&swiss_path, swiss_text).expect("write swiss toy");

    app.uniprot_entry_id = "TOY_DIALOG".to_string();
    app.uniprot_swiss_path = swiss_path.display().to_string();
    app.import_uniprot_swiss_prot_from_dialog();

    let entries = app.engine.read().unwrap().list_uniprot_entries();
    assert_eq!(entries.len(), 1);
    assert_eq!(entries[0].entry_id, "TOY_DIALOG");
    assert!(
        app.uniprot_status.contains("UniProt import: ok"),
        "status was: {}",
        app.uniprot_status
    );

    let recent = app.recent_uniprot_entries_for_dialog(5);
    assert_eq!(recent.len(), 1);
    assert_eq!(recent[0].entry_id, "TOY_DIALOG");
}

#[test]
fn uniprot_dialog_select_row_sets_entry_id_without_importing_sequence() {
    let mut app = GENtleApp::default();
    let temp = tempdir().expect("tempdir");
    let swiss_path = temp.path().join("toy_uniprot_dialog_use.txt");
    let swiss_text = r#"ID   TOY2_HUMAN              Reviewed;         12 AA.
AC   PUSE1;
DE   RecName: Full=Toy use protein;
GN   Name=TOY2;
OS   Homo sapiens (Human).
FT   DOMAIN          2..6
FT                   /note="toy segment"
SQ   SEQUENCE   12 AA;  1200 MW;  0000000000000000 CRC64;
     MEEPQSDPSVEP
//
"#;
    fs::write(&swiss_path, swiss_text).expect("write swiss toy");

    app.uniprot_entry_id = "TOY_USE".to_string();
    app.uniprot_swiss_path = swiss_path.display().to_string();
    app.import_uniprot_swiss_prot_from_dialog();

    let recent = app.recent_uniprot_entries_for_dialog(5);
    assert_eq!(recent.len(), 1);
    app.uniprot_entry_id.clear();
    app.uniprot_query.clear();
    let row = recent.first().expect("recent row");
    app.uniprot_entry_id = row.entry_id.clone();
    app.uniprot_query = row.accession.clone();

    let state = app.engine.read().unwrap().state().clone();
    assert!(!state.sequences.contains_key("TOY_USE"));
    assert_eq!(app.uniprot_entry_id, "TOY_USE");
    assert_eq!(app.uniprot_query, "PUSE1");
}

#[test]
fn uniprot_feature_speed_resolution_summary_mentions_source_and_reference() {
    let report = UniprotFeatureCodingDnaQueryReport {
        schema: "gentle.uniprot_feature_coding_dna_query.v1".to_string(),
        projection_id: "P1@seq".to_string(),
        entry_id: "P1".to_string(),
        seq_id: "seq".to_string(),
        feature_query: "DNA-binding".to_string(),
        transcript_filter: Some("TX1".to_string()),
        query_mode: UniprotFeatureCodingDnaQueryMode::Both,
        requested_translation_speed_profile: Some(TranslationSpeedProfile::Mouse),
        resolved_translation_speed_profile: Some(TranslationSpeedProfile::Mouse),
        resolved_translation_speed_profile_source: Some(
            TranslationSpeedProfileSource::SourceOrganismScientificName,
        ),
        resolved_translation_speed_reference_species: Some("Mus musculus domesticus".to_string()),
        match_count: 1,
        matches: vec![],
        warnings: vec![],
    };
    let summary = GENtleApp::format_uniprot_feature_speed_resolution_summary(&report);
    assert!(summary.contains("requested=mouse"), "{summary}");
    assert!(summary.contains("resolved=mouse"), "{summary}");
    assert!(
        summary.contains("source=source_organism_scientific_name"),
        "{summary}"
    );
    assert!(summary.contains("ref=Mus musculus domesticus"), "{summary}");
}

#[test]
fn reverse_translation_speed_resolution_summary_mentions_source_and_reference() {
    let report = ReverseTranslationReport {
        schema: "gentle.reverse_translation_report.v1".to_string(),
        report_id: "rt1".to_string(),
        protein_seq_id: "prot".to_string(),
        coding_seq_id: "prot__coding".to_string(),
        generated_at_unix_ms: 1,
        op_id: Some("op1".to_string()),
        run_id: Some("run1".to_string()),
        requested_output_id: None,
        effective_output_id: "prot__coding".to_string(),
        protein_length_aa: 3,
        coding_length_bp: 9,
        translation_table: 11,
        translation_table_label: "Bacterial, Archaeal and Plant Plastid".to_string(),
        translation_table_source: "feature_qualifier_hint".to_string(),
        translation_context_organism: Some("Escherichia coli".to_string()),
        translation_context_organelle: None,
        requested_speed_profile: Some(TranslationSpeedProfile::Ecoli),
        resolved_speed_profile: Some(TranslationSpeedProfile::Ecoli),
        resolved_speed_profile_source: Some(TranslationSpeedProfileSource::FeatureQualifierHint),
        translation_speed_reference_species: Some("Escherichia coli".to_string()),
        speed_mark: Some(TranslationSpeedMark::Slow),
        target_anneal_tm_c: Some(58.0),
        anneal_window_bp: 9,
        preferred_synonymous_choice_count: 3,
        alternative_synonymous_choice_count: 0,
        fallback_unknown_codon_count: 0,
        gc_fraction: Some(0.5),
        realized_anneal_tm_c: Some(58.0),
        warnings: vec![],
    };
    let summary = GENtleApp::format_reverse_translation_speed_resolution_summary(&report);
    assert!(summary.contains("requested=ecoli"), "{summary}");
    assert!(summary.contains("resolved=ecoli"), "{summary}");
    assert!(
        summary.contains("source=feature_qualifier_hint"),
        "{summary}"
    );
    assert!(summary.contains("ref=Escherichia coli"), "{summary}");
}

#[test]
fn reverse_translate_protein_from_dialog_sets_report_and_status() {
    let mut protein = DNAsequence::from_sequence("MKP").expect("protein");
    protein.set_name("Toy protein");
    protein.set_molecule_type("protein");
    protein.features_mut().push(gb_io::seq::Feature {
        kind: "Protein".into(),
        location: gb_io::seq::Location::simple_range(0, 3),
        qualifiers: vec![
            (
                "translation_speed_profile_hint".into(),
                Some("ecoli".to_string()),
            ),
            ("translation_table".into(), Some("11".to_string())),
        ],
    });
    let mut state = ProjectState::default();
    state.sequences.insert("prot".to_string(), protein);
    let mut app = GENtleApp::default();
    app.engine = Arc::new(RwLock::new(GentleEngine::from_state(state)));
    app.reverse_translate_protein_seq_id = "prot".to_string();
    app.reverse_translate_speed_mark = Some(TranslationSpeedMark::Slow);
    app.reverse_translate_target_anneal_tm_c = "58.0".to_string();
    app.reverse_translate_anneal_window_bp = "9".to_string();

    app.reverse_translate_protein_from_dialog();

    let report = app
        .reverse_translation_report
        .as_ref()
        .expect("reverse translation report");
    assert_eq!(report.protein_seq_id, "prot");
    assert_eq!(report.translation_table, 11);
    assert_eq!(
        report.resolved_speed_profile,
        Some(TranslationSpeedProfile::Ecoli)
    );
    assert_eq!(report.speed_mark, Some(TranslationSpeedMark::Slow));
    assert!(
        app.uniprot_status.contains("translation table: 11"),
        "{}",
        app.uniprot_status
    );
    assert!(
        app.uniprot_status
            .contains("speed profile: requested=auto | resolved=ecoli"),
        "{}",
        app.uniprot_status
    );
    assert!(
        app.engine
            .read()
            .unwrap()
            .state()
            .sequences
            .contains_key(&report.coding_seq_id)
    );
}

#[test]
fn protease_digest_from_dialog_sets_report_and_status() {
    let mut protein = DNAsequence::from_sequence("MAKRPTRKAA").expect("protein");
    protein.set_name("Toy protein");
    protein.set_molecule_type("protein");
    protein.features_mut().push(gb_io::seq::Feature {
        kind: "Protein".into(),
        location: gb_io::seq::Location::simple_range(0, 10),
        qualifiers: vec![("transcript_id".into(), Some("tx1".to_string()))],
    });
    let mut state = ProjectState::default();
    state.sequences.insert("prot".to_string(), protein);
    let mut app = GENtleApp::default();
    app.engine = Arc::new(RwLock::new(GentleEngine::from_state(state)));
    app.reverse_translate_protein_seq_id = "prot".to_string();
    app.protease_digest_names = "Trypsin".to_string();
    app.protease_digest_min_length_aa = "1".to_string();
    app.protease_digest_output_prefix = "prot_trypsin".to_string();
    app.protease_digest_materialize = true;

    app.digest_selected_protein_from_dialog();

    let report = app
        .protease_digest_report
        .as_ref()
        .expect("protease digest report");
    assert_eq!(report.source_seq_id, "prot");
    assert_eq!(report.source_transcript_id.as_deref(), Some("tx1"));
    assert_eq!(report.cleavage_site_count, 3);
    assert_eq!(report.peptide_count, 4);
    assert!(app.uniprot_status.contains("Protease digest: ok"));
    assert!(
        app.engine
            .read()
            .unwrap()
            .state()
            .sequences
            .contains_key(report.created_seq_ids[0].as_str())
    );
}

#[test]
fn resolve_uniprot_projection_id_from_dialog_fields_uses_entry_and_seq_default() {
    let mut app = GENtleApp::default();
    app.uniprot_entry_id = "PTEST1".to_string();
    app.uniprot_map_seq_id = "toy_seq".to_string();
    assert_eq!(
        app.resolve_uniprot_projection_id_from_dialog_fields()
            .as_deref(),
        Some("PTEST1@toy_seq")
    );

    app.uniprot_map_projection_id = "custom_projection".to_string();
    assert_eq!(
        app.resolve_uniprot_projection_id_from_dialog_fields()
            .as_deref(),
        Some("custom_projection")
    );
}

#[test]
fn default_uniprot_projection_svg_file_name_uses_seq_and_projection_id() {
    assert_eq!(
        GENtleApp::default_uniprot_projection_svg_file_name("grch38_tp53", "tp53_uniprot_p04637"),
        "grch38_tp53_tp53_uniprot_p04637.protein_mapping.svg"
    );
}

#[test]
fn default_transcript_protein_svg_file_name_uses_seq_and_optional_transcript() {
    assert_eq!(
        GENtleApp::default_transcript_protein_svg_file_name("grch38_tp53", None),
        "grch38_tp53_derived_proteins.protein_compare.svg"
    );
    assert_eq!(
        GENtleApp::default_transcript_protein_svg_file_name("grch38_tp53", Some("ENST00000269305")),
        "grch38_tp53_ENST00000269305.protein_compare.svg"
    );
}

#[test]
fn default_ensembl_protein_svg_file_name_uses_seq_entry_and_optional_transcript() {
    assert_eq!(
        GENtleApp::default_ensembl_protein_svg_file_name("grch38_tp53", None, "ENSP00000269305"),
        "grch38_tp53_ENSP00000269305.protein_compare.svg"
    );
    assert_eq!(
        GENtleApp::default_ensembl_protein_svg_file_name(
            "grch38_tp53",
            Some("ENST00000269305"),
            "ENSP00000269305"
        ),
        "grch38_tp53_ENST00000269305_ENSP00000269305.protein_compare.svg"
    );
}

#[test]
fn sort_ensembl_protein_entries_recent_prefers_newer_rows_then_entry_id() {
    let rows = vec![
        EnsemblProteinEntrySummary {
            entry_id: "ENSP_B".to_string(),
            protein_id: "ENSP_B".to_string(),
            transcript_id: "ENST_B".to_string(),
            imported_at_unix_ms: 100,
            ..Default::default()
        },
        EnsemblProteinEntrySummary {
            entry_id: "ENSP_A".to_string(),
            protein_id: "ENSP_A".to_string(),
            transcript_id: "ENST_A".to_string(),
            imported_at_unix_ms: 100,
            ..Default::default()
        },
        EnsemblProteinEntrySummary {
            entry_id: "ENSP_NEW".to_string(),
            protein_id: "ENSP_NEW".to_string(),
            transcript_id: "ENST_NEW".to_string(),
            imported_at_unix_ms: 250,
            ..Default::default()
        },
    ];
    let sorted = GENtleApp::sort_ensembl_protein_entries_recent(rows);
    let ordered_ids = sorted
        .iter()
        .map(|row| row.entry_id.as_str())
        .collect::<Vec<_>>();
    assert_eq!(ordered_ids, vec!["ENSP_NEW", "ENSP_A", "ENSP_B"]);
}

#[test]
fn protein_feature_filter_from_dialog_parses_include_and_exclude_keys() {
    let mut app = GENtleApp::default();
    app.protein_feature_key_include = "DOMAIN, DNA_BIND PF02196 DOMAIN".to_string();
    app.protein_feature_key_exclude = "CONFLICT; REGION\nCHAIN".to_string();
    let filter = app.protein_feature_filter_from_dialog();
    assert_eq!(
        filter.include_feature_keys,
        vec!["DOMAIN", "DNA_BIND", "PF02196"]
    );
    assert_eq!(
        filter.exclude_feature_keys,
        vec!["CONFLICT", "REGION", "CHAIN"]
    );
}

#[test]
fn summarize_ensembl_protein_feature_keys_prefers_frequency_then_name() {
    let entry = EnsemblProteinEntry {
        entry_id: "ENSP_TEST".to_string(),
        protein_id: "ENSP_TEST".to_string(),
        transcript_id: "ENST_TEST".to_string(),
        features: vec![
            EnsemblProteinFeature {
                feature_key: "DOMAIN".to_string(),
                ..Default::default()
            },
            EnsemblProteinFeature {
                feature_key: "DOMAIN".to_string(),
                ..Default::default()
            },
            EnsemblProteinFeature {
                feature_key: "REGION".to_string(),
                ..Default::default()
            },
            EnsemblProteinFeature {
                feature_key: String::new(),
                feature_type: "signal_peptide".to_string(),
                ..Default::default()
            },
        ],
        ..Default::default()
    };

    let summary = GENtleApp::summarize_ensembl_protein_feature_keys(&entry, 8);
    assert_eq!(
        summary,
        vec![
            "DOMAIN x2".to_string(),
            "REGION x1".to_string(),
            "signal_peptide x1".to_string()
        ]
    );
}

#[test]
fn selected_ensembl_protein_entry_for_dialog_resolves_aliases() {
    let mut state = ProjectState::default();
    state.metadata.insert(
        "ensembl_protein_entries".to_string(),
        serde_json::json!({
            "schema": "gentle.ensembl_protein_entries.v1",
            "updated_at_unix_ms": 1u128,
            "entries": {
                "ENSP_TEST": {
                    "schema": "gentle.ensembl_protein_entry.v1",
                    "entry_id": "ENSP_TEST",
                    "protein_id": "ENSP_TEST",
                    "transcript_id": "ENST_TEST",
                    "aliases": ["ensp_test_alias"],
                    "sequence": "MSTNPKPQR",
                    "sequence_length": 9,
                    "features": [],
                    "source": "ensembl_rest",
                    "imported_at_unix_ms": 1u128,
                    "transcript_lookup_source_url": "https://example.invalid/transcript",
                    "sequence_source_url": "https://example.invalid/sequence",
                    "feature_source_url": "https://example.invalid/features",
                    "raw_transcript_lookup_json": "{}",
                    "raw_sequence_json": "{}",
                    "raw_feature_json": "[]"
                }
            }
        }),
    );

    let mut app = GENtleApp::default();
    app.engine = Arc::new(RwLock::new(GentleEngine::from_state(state)));
    app.ensembl_protein_entry_id = "ensp_test_alias".to_string();

    let selected = app
        .selected_ensembl_protein_entry_for_dialog()
        .expect("selected entry");
    assert_eq!(selected.entry_id, "ENSP_TEST");
    assert_eq!(selected.transcript_id, "ENST_TEST");
}

#[test]
fn recent_uniprot_projections_for_dialog_filters_selected_entry_and_seq() {
    let temp = tempdir().expect("tempdir");
    let swiss_path = temp.path().join("toy_uniprot_projection_dialog.txt");
    let swiss_text = r#"ID   TOY1_HUMAN              Reviewed;         30 AA.
AC   PTEST1;
DE   RecName: Full=Toy DNA-binding protein;
GN   Name=TOY1;
OS   Homo sapiens (Human).
DR   Ensembl; TX1; ENSPTOY1; ENSGTOY1.
FT   DOMAIN          2..8
FT                   /note="toy domain"
SQ   SEQUENCE   30 AA;  3333 MW;  0000000000000000 CRC64;
     MEEPQSDPSV EPPLSQETFSDLWKLLPEN
//
"#;
    fs::write(&swiss_path, swiss_text).expect("write swiss toy");

    let mut state = ProjectState::default();
    let mut dna = DNAsequence::from_sequence(&"ACGT".repeat(300)).expect("valid DNA");
    dna.features_mut().push(gb_io::seq::Feature {
        kind: "mRNA".into(),
        location: gb_io::seq::Location::simple_range(99, 360),
        qualifiers: vec![
            ("gene".into(), Some("TOY1".to_string())),
            ("transcript_id".into(), Some("TX1".to_string())),
            ("label".into(), Some("TX1".to_string())),
            (
                "cds_ranges_1based".into(),
                Some("100-180,300-360".to_string()),
            ),
        ]
        .into_iter()
        .collect(),
    });
    state.sequences.insert("toy_seq".to_string(), dna);

    let mut app = GENtleApp::default();
    app.engine = Arc::new(RwLock::new(GentleEngine::from_state(state)));
    app.uniprot_swiss_path = swiss_path.display().to_string();
    app.import_uniprot_swiss_prot_from_dialog();
    {
        let mut engine = app.engine.write().expect("engine lock");
        engine
            .apply(Operation::ProjectUniprotToGenome {
                seq_id: "toy_seq".to_string(),
                entry_id: "PTEST1".to_string(),
                projection_id: Some("toy_projection".to_string()),
                transcript_id: None,
            })
            .expect("project uniprot");
    }

    app.uniprot_entry_id = "PTEST1".to_string();
    app.uniprot_map_seq_id = "toy_seq".to_string();
    let rows = app.recent_uniprot_projections_for_dialog(5);
    assert_eq!(rows.len(), 1);
    assert_eq!(rows[0].projection_id, "toy_projection");
    assert_eq!(rows[0].entry_id, "PTEST1");
    assert_eq!(rows[0].seq_id, "toy_seq");
}

#[test]
fn genbank_dialog_fetch_imports_sequence_and_opens_window() {
    let _lock = crate::genomes::genbank_env_lock()
        .lock()
        .unwrap_or_else(|e| e.into_inner());
    let mut app = GENtleApp::default();
    let temp = tempdir().expect("tempdir");
    let mock_dir = temp.path().join("mock");
    fs::create_dir_all(&mock_dir).expect("mock dir");
    fs::copy(
        "test_files/tp73.ncbi.gb",
        mock_dir.join("NC_000001.gbwithparts"),
    )
    .expect("copy genbank fixture");
    let efetch_template = format!("file://{}/{{accession}}.{{rettype}}", mock_dir.display());
    let _env_guard = EnvVarGuard::set("GENTLE_NCBI_EFETCH_URL", &efetch_template);

    app.genbank_accession = "NC_000001".to_string();
    app.genbank_as_id = "tp73_gui_fetch".to_string();
    app.fetch_genbank_accession_from_dialog();

    let state = app.engine.read().unwrap().state().clone();
    assert!(state.sequences.contains_key("tp73_gui_fetch"));
    assert_eq!(
        app.new_windows
            .first()
            .and_then(|window| window.sequence_id()),
        Some("tp73_gui_fetch".to_string())
    );
    assert!(
        app.genbank_status.contains("GenBank fetch: ok"),
        "status was: {}",
        app.genbank_status
    );
}

#[test]
fn open_genbank_dialog_seeds_dbsnp_tutorial_defaults() {
    let mut app = GENtleApp::default();
    app.dbsnp_rs_id.clear();
    app.dbsnp_genome_id.clear();
    app.dbsnp_flank_bp.clear();

    app.open_genbank_dialog();

    assert_eq!(app.dbsnp_rs_id, DEFAULT_DBSNP_TUTORIAL_RS_ID);
    assert_eq!(app.dbsnp_genome_id, "Human GRCh38 Ensembl 113");
    assert_eq!(app.dbsnp_flank_bp, "3000");
}

#[test]
fn poll_dbsnp_fetch_task_updates_runtime_status_from_worker_messages() {
    let mut app = GENtleApp::default();
    let (tx, rx) = mpsc::channel::<DbSnpFetchTaskMessage>();
    app.dbsnp_fetch_task = Some(DbSnpFetchTask {
        started: Instant::now(),
        receiver: rx,
    });
    tx.send(DbSnpFetchTaskMessage::Progress(DbSnpFetchProgress {
        rs_id: "rs123".to_string(),
        genome_id: "ToyGenome".to_string(),
        stage: DbSnpFetchStage::WaitResponse,
        detail: "Waiting for dbSNP response from 'file:///tmp/mock/123.json'".to_string(),
    }))
    .expect("send dbsnp progress");

    app.poll_dbsnp_fetch_task(&egui::Context::default());

    assert!(
        app.dbsnp_status.contains("Waiting for dbSNP response"),
        "status was: {}",
        app.dbsnp_status
    );
    assert!(app.dbsnp_fetch_task.is_some());
}

#[test]
fn dbsnp_dialog_fetch_extracts_region_and_opens_window() {
    let _lock = crate::genomes::genbank_env_lock().lock().unwrap();
    let mut app = GENtleApp::default();
    let temp = tempdir().expect("tempdir");
    let (catalog_path, cache_dir_str) = write_toy_prepare_catalog(temp.path());
    let mock_dir = temp.path().join("mock_dbsnp");
    fs::create_dir_all(&mock_dir).expect("mock dir");
    fs::write(
        mock_dir.join("123.json"),
        r#"{
  "refsnp_id": "123",
  "primary_snapshot_data": {
    "placements_with_allele": [
      {
        "seq_id": "chr1",
        "is_ptlp": true,
        "placement_annot": {
          "seq_id_traits_by_assembly": [
            {
              "assembly_name": "ToyGenome.1",
              "is_top_level": true,
              "is_alt": false,
              "is_patch": false,
              "is_chromosome": true
            }
          ]
        },
        "alleles": [
          {
            "allele": {
              "spdi": {
                "seq_id": "chr1",
                "position": 3,
                "deleted_sequence": "A",
                "inserted_sequence": "G"
              }
            }
          }
        ]
      }
    ],
    "allele_annotations": [
      {
        "assembly_annotation": [
          {
            "genes": [
              {
                "locus": "tagA"
              }
            ]
          }
        ]
      }
    ]
  }
}
"#,
    )
    .expect("write mock dbsnp payload");
    let refsnp_template = format!("file://{}/{{refsnp_id}}.json", mock_dir.display());
    let _env_guard = EnvVarGuard::set("GENTLE_NCBI_DBSNP_REFSNP_URL", &refsnp_template);

    app.reference_genome_catalog_path = catalog_path.clone();
    app.reference_genome_cache_dir = cache_dir_str.clone();
    app.engine
        .write()
        .unwrap()
        .apply(Operation::PrepareGenome {
            genome_id: "ToyGenome".to_string(),
            catalog_path: Some(catalog_path.clone()),
            cache_dir: Some(cache_dir_str.clone()),
            timeout_seconds: None,
        })
        .expect("prepare ToyGenome");

    app.dbsnp_rs_id = "rs123".to_string();
    app.dbsnp_genome_id = "ToyGenome".to_string();
    app.dbsnp_flank_bp = "2".to_string();
    app.dbsnp_output_id = "rs123_gui_fetch".to_string();
    app.fetch_dbsnp_region_from_dialog();
    assert!(
        app.dbsnp_status.contains("Starting dbSNP fetch"),
        "status was: {}",
        app.dbsnp_status
    );
    let wait_started = Instant::now();
    while app.dbsnp_fetch_task.is_some() && wait_started.elapsed() < Duration::from_secs(15) {
        app.poll_dbsnp_fetch_task(&egui::Context::default());
        std::thread::sleep(Duration::from_millis(10));
    }
    if app.dbsnp_fetch_task.is_some() {
        app.poll_dbsnp_fetch_task(&egui::Context::default());
    }
    assert!(
        app.dbsnp_fetch_task.is_none(),
        "dbSNP task did not finish; last status: {}",
        app.dbsnp_status
    );

    let state = app.engine.read().unwrap().state().clone();
    let extracted = state
        .sequences
        .get("rs123_gui_fetch")
        .expect("gui dbsnp extracted sequence");
    let marker = extracted
        .features()
        .iter()
        .find(|feature| {
            feature.kind.to_string().eq_ignore_ascii_case("variation")
                && feature
                    .qualifier_values("label")
                    .any(|value| value == "rs123")
        })
        .expect("gui dbsnp marker");
    assert_eq!(
        marker.location.find_bounds().expect("gui marker bounds"),
        (2, 3)
    );
    assert_eq!(
        app.new_windows
            .first()
            .and_then(|window| window.sequence_id()),
        Some("rs123_gui_fetch".to_string())
    );
    assert!(
        app.dbsnp_status.contains("dbSNP fetch: ok"),
        "status was: {}",
        app.dbsnp_status
    );
    assert!(
        app.dbsnp_status.contains("annotation: requested=full"),
        "status was: {}",
        app.dbsnp_status
    );
}

#[test]
fn request_prepare_cancel_is_idempotent() {
    let mut app = GENtleApp::default();
    let (_tx, rx) = mpsc::channel::<GenomePrepareTaskMessage>();
    app.genome_prepare_task = Some(GenomePrepareTask {
        job_id: 42,
        started: Instant::now(),
        cancel_requested: Arc::new(AtomicBool::new(false)),
        timeout_seconds: None,
        mode: GenomePrepareLaunchMode::Prepare,
        genome_id: "ToyGenome".to_string(),
        scope: GenomeDialogScope::Reference,
        catalog_path: "assets/genomes.json".to_string(),
        cache_dir: "data/genomes".to_string(),
        receiver: rx,
    });

    app.request_prepare_task_cancel("test");
    app.request_prepare_task_cancel("test");

    let cancel_events = app
        .job_event_log
        .iter()
        .filter(|event| {
            event.kind == BackgroundJobKind::PrepareGenome
                && event.phase == BackgroundJobEventPhase::CancelRequested
                && event.job_id == Some(42)
        })
        .count();
    assert_eq!(cancel_events, 1);
    assert!(
        app.genome_prepare_task
            .as_ref()
            .unwrap()
            .cancel_requested
            .load(Ordering::Relaxed)
    );
}

#[test]
fn start_prepare_reference_genome_with_mode_sets_queued_progress_immediately() {
    let mut app = GENtleApp::default();
    app.genome_id = "ToyGenome".to_string();
    app.genome_catalog_path = "/definitely/missing/catalog.json".to_string();

    app.start_prepare_reference_genome_with_mode(GenomePrepareLaunchMode::ReindexCachedFiles);

    let task = app
        .genome_prepare_task
        .as_ref()
        .expect("prepare task should be registered immediately");
    assert_eq!(task.mode, GenomePrepareLaunchMode::ReindexCachedFiles);
    let progress = app
        .genome_prepare_progress
        .as_ref()
        .expect("queued progress should be visible immediately");
    assert_eq!(progress.phase, "queued");
    assert_eq!(progress.genome_id, "ToyGenome");
}

#[test]
fn start_prepare_reference_genome_with_mode_initializes_prepare_steps_immediately() {
    let _env_guard = EnvVarGuard::set(
        crate::genomes::MAKEBLASTDB_ENV_BIN,
        "__gentle_makeblastdb_missing_for_test__",
    );
    let td = tempdir().expect("tempdir");
    let (catalog_path, cache_dir) = write_toy_prepare_catalog(td.path());
    let mut app = GENtleApp::default();
    app.genome_id = "ToyGenome".to_string();
    app.genome_catalog_path = catalog_path;
    app.genome_cache_dir = cache_dir;

    app.start_prepare_reference_genome_with_mode(GenomePrepareLaunchMode::Prepare);

    assert_eq!(app.genome_prepare_steps.len(), 5);
    assert_eq!(
        app.genome_prepare_steps
            .iter()
            .map(|step| step.step_id)
            .collect::<Vec<_>>(),
        vec![
            PrepareGenomeStepId::Sequence,
            PrepareGenomeStepId::Annotation,
            PrepareGenomeStepId::FastaIndex,
            PrepareGenomeStepId::GeneIndex,
            PrepareGenomeStepId::BlastIndex,
        ]
    );
}

#[test]
fn apply_prepare_progress_updates_only_matching_step() {
    let mut app = GENtleApp::default();
    app.reset_prepare_step_state_from_plan(Some(make_prepare_plan(&[
        PrepareGenomeStepId::Sequence,
        PrepareGenomeStepId::Annotation,
        PrepareGenomeStepId::BlastIndex,
    ])));

    app.apply_prepare_progress_to_steps(&PrepareGenomeProgress {
        genome_id: "ToyGenome".to_string(),
        phase: "download_annotation".to_string(),
        item: "annotation.gtf".to_string(),
        bytes_done: 50,
        bytes_total: Some(100),
        percent: Some(50.0),
        step_id: Some(PrepareGenomeStepId::Annotation),
        step_label: Some("Annotation".to_string()),
    });

    assert_eq!(
        app.genome_prepare_steps[0].status,
        PrepareGenomeUiStepStatus::Completed
    );
    assert_eq!(
        app.genome_prepare_steps[1].status,
        PrepareGenomeUiStepStatus::Running
    );
    assert_eq!(app.genome_prepare_steps[1].progress_fraction, Some(0.5));
    assert_eq!(
        app.genome_prepare_steps[2].status,
        PrepareGenomeUiStepStatus::Pending
    );
}

#[test]
fn apply_prepared_genome_inspection_marks_existing_outputs_complete() {
    let mut app = GENtleApp::default();
    app.reset_prepare_step_state_from_plan(Some(make_prepare_plan(&[
        PrepareGenomeStepId::ResetIndexes,
        PrepareGenomeStepId::Sequence,
        PrepareGenomeStepId::Annotation,
        PrepareGenomeStepId::FastaIndex,
        PrepareGenomeStepId::GeneIndex,
        PrepareGenomeStepId::BlastIndex,
    ])));

    let inspection = make_prepared_genome_inspection();
    app.apply_prepared_genome_inspection_to_steps(Some(&inspection));

    assert_eq!(
        app.genome_prepare_steps
            .iter()
            .map(|step| step.status)
            .collect::<Vec<_>>(),
        vec![
            PrepareGenomeUiStepStatus::Pending,
            PrepareGenomeUiStepStatus::Completed,
            PrepareGenomeUiStepStatus::Completed,
            PrepareGenomeUiStepStatus::Completed,
            PrepareGenomeUiStepStatus::Completed,
            PrepareGenomeUiStepStatus::Completed,
        ]
    );
    assert_eq!(
        app.genome_prepare_steps[1].detail,
        "/tmp/toy_install/sequence.fa"
    );
    assert_eq!(
        app.genome_prepare_steps[4].detail,
        "/tmp/toy_install/transcripts.json"
    );
}

#[test]
fn prepare_eta_appears_after_progress_advances() {
    let mut app = GENtleApp::default();
    app.reset_prepare_step_state_from_plan(Some(make_prepare_plan(&[
        PrepareGenomeStepId::Annotation,
    ])));
    let observed_at = Instant::now();

    app.apply_prepare_progress_to_steps_at(
        &PrepareGenomeProgress {
            genome_id: "ToyGenome".to_string(),
            phase: "download_annotation".to_string(),
            item: "annotation.gtf.gz".to_string(),
            bytes_done: 10,
            bytes_total: Some(100),
            percent: Some(10.0),
            step_id: Some(PrepareGenomeStepId::Annotation),
            step_label: Some("Annotation".to_string()),
        },
        observed_at,
    );
    assert_eq!(app.genome_prepare_steps[0].eta_remaining, None);

    app.apply_prepare_progress_to_steps_at(
        &PrepareGenomeProgress {
            genome_id: "ToyGenome".to_string(),
            phase: "download_annotation".to_string(),
            item: "annotation.gtf.gz".to_string(),
            bytes_done: 40,
            bytes_total: Some(100),
            percent: Some(40.0),
            step_id: Some(PrepareGenomeStepId::Annotation),
            step_label: Some("Annotation".to_string()),
        },
        observed_at + Duration::from_secs(6),
    );

    let eta = app.genome_prepare_steps[0]
        .eta_remaining
        .expect("eta should be visible after progress advances");
    assert!(eta.as_secs() > 0, "eta was {eta:?}");
}

#[test]
fn prepare_eta_resets_when_step_changes() {
    let mut app = GENtleApp::default();
    app.reset_prepare_step_state_from_plan(Some(make_prepare_plan(&[
        PrepareGenomeStepId::Sequence,
        PrepareGenomeStepId::Annotation,
    ])));
    let observed_at = Instant::now();

    app.apply_prepare_progress_to_steps_at(
        &PrepareGenomeProgress {
            genome_id: "ToyGenome".to_string(),
            phase: "download_sequence".to_string(),
            item: "sequence.fa.gz".to_string(),
            bytes_done: 25,
            bytes_total: Some(100),
            percent: Some(25.0),
            step_id: Some(PrepareGenomeStepId::Sequence),
            step_label: Some("Sequence".to_string()),
        },
        observed_at,
    );
    app.apply_prepare_progress_to_steps_at(
        &PrepareGenomeProgress {
            genome_id: "ToyGenome".to_string(),
            phase: "download_sequence".to_string(),
            item: "sequence.fa.gz".to_string(),
            bytes_done: 55,
            bytes_total: Some(100),
            percent: Some(55.0),
            step_id: Some(PrepareGenomeStepId::Sequence),
            step_label: Some("Sequence".to_string()),
        },
        observed_at + Duration::from_secs(5),
    );
    assert!(app.genome_prepare_steps[0].eta_remaining.is_some());

    app.apply_prepare_progress_to_steps_at(
        &PrepareGenomeProgress {
            genome_id: "ToyGenome".to_string(),
            phase: "download_annotation".to_string(),
            item: "annotation.gtf.gz".to_string(),
            bytes_done: 5,
            bytes_total: Some(100),
            percent: Some(5.0),
            step_id: Some(PrepareGenomeStepId::Annotation),
            step_label: Some("Annotation".to_string()),
        },
        observed_at + Duration::from_secs(8),
    );

    assert_eq!(app.genome_prepare_steps[0].eta_remaining, None);
    assert_eq!(app.genome_prepare_steps[1].eta_remaining, None);
}

#[test]
fn prepare_eta_hidden_for_completed_and_indeterminate_steps() {
    let mut app = GENtleApp::default();
    app.reset_prepare_step_state_from_plan(Some(make_prepare_plan(&[
        PrepareGenomeStepId::Sequence,
        PrepareGenomeStepId::BlastIndex,
    ])));
    let observed_at = Instant::now();

    app.apply_prepare_progress_to_steps_at(
        &PrepareGenomeProgress {
            genome_id: "ToyGenome".to_string(),
            phase: "download_sequence".to_string(),
            item: "sequence.fa.gz".to_string(),
            bytes_done: 100,
            bytes_total: Some(100),
            percent: Some(100.0),
            step_id: Some(PrepareGenomeStepId::Sequence),
            step_label: Some("Sequence".to_string()),
        },
        observed_at,
    );
    app.apply_prepare_progress_to_steps_at(
        &PrepareGenomeProgress {
            genome_id: "ToyGenome".to_string(),
            phase: "index_blast".to_string(),
            item: "makeblastdb".to_string(),
            bytes_done: 0,
            bytes_total: None,
            percent: None,
            step_id: Some(PrepareGenomeStepId::BlastIndex),
            step_label: Some("BLAST Index".to_string()),
        },
        observed_at + Duration::from_secs(3),
    );

    assert_eq!(app.genome_prepare_steps[0].eta_remaining, None);
    assert_eq!(app.genome_prepare_steps[1].eta_remaining, None);
}

#[test]
fn finalize_prepare_steps_success_marks_all_steps_complete() {
    let mut app = GENtleApp::default();
    app.reset_prepare_step_state_from_plan(Some(make_prepare_plan(&[
        PrepareGenomeStepId::Sequence,
        PrepareGenomeStepId::Annotation,
        PrepareGenomeStepId::BlastIndex,
    ])));

    app.finalize_prepare_steps_success();

    assert!(
        app.genome_prepare_steps
            .iter()
            .all(|step| step.status == PrepareGenomeUiStepStatus::Completed)
    );
}

#[test]
fn finalize_prepare_steps_failure_marks_running_step_and_preserves_checklist() {
    let mut app = GENtleApp::default();
    app.reset_prepare_step_state_from_plan(Some(make_prepare_plan(&[
        PrepareGenomeStepId::Sequence,
        PrepareGenomeStepId::BlastIndex,
    ])));
    app.genome_prepare_steps[0].status = PrepareGenomeUiStepStatus::Completed;
    app.genome_prepare_steps[0].progress_fraction = Some(1.0);
    app.genome_prepare_steps[1].status = PrepareGenomeUiStepStatus::Running;

    app.finalize_prepare_steps_failure(true);

    assert_eq!(app.genome_prepare_steps.len(), 2);
    assert_eq!(
        app.genome_prepare_steps[0].status,
        PrepareGenomeUiStepStatus::Completed
    );
    assert_eq!(
        app.genome_prepare_steps[1].status,
        PrepareGenomeUiStepStatus::Cancelled
    );
}

#[test]
fn prepare_step_summary_reports_completed_count_and_overall_fraction() {
    let mut app = GENtleApp::default();
    app.reset_prepare_step_state_from_plan(Some(make_prepare_plan(&[
        PrepareGenomeStepId::Sequence,
        PrepareGenomeStepId::Annotation,
        PrepareGenomeStepId::BlastIndex,
    ])));
    app.genome_prepare_steps[0].status = PrepareGenomeUiStepStatus::Completed;
    app.genome_prepare_steps[0].progress_fraction = Some(1.0);
    app.genome_prepare_steps[1].status = PrepareGenomeUiStepStatus::Running;
    app.genome_prepare_steps[1].progress_fraction = Some(0.5);

    let (current_step, completed_steps, total_steps, overall, eta) =
        app.prepare_step_summary().expect("summary");

    assert_eq!(current_step, "Annotation");
    assert_eq!(completed_steps, 1);
    assert_eq!(total_steps, 3);
    assert_eq!(eta, None);
    assert!((overall - 0.5).abs() < 0.001, "overall was {overall}");
}

#[test]
fn request_blast_cancel_is_idempotent() {
    let mut app = GENtleApp::default();
    let (_tx, rx) = mpsc::channel::<GenomeBlastTaskMessage>();
    app.genome_blast_task = Some(GenomeBlastTask {
        job_id: 43,
        started: Instant::now(),
        cancel_requested: Arc::new(AtomicBool::new(false)),
        receiver: rx,
    });

    app.request_blast_task_cancel("test");
    app.request_blast_task_cancel("test");

    let cancel_events = app
        .job_event_log
        .iter()
        .filter(|event| {
            event.kind == BackgroundJobKind::BlastGenome
                && event.phase == BackgroundJobEventPhase::CancelRequested
                && event.job_id == Some(43)
        })
        .count();
    assert_eq!(cancel_events, 1);
    assert!(
        app.genome_blast_task
            .as_ref()
            .expect("blast task should still exist")
            .cancel_requested
            .load(Ordering::Relaxed)
    );
}

#[test]
fn poll_prepare_ignores_stale_job_messages() {
    let mut app = GENtleApp::default();
    let (tx, rx) = mpsc::channel::<GenomePrepareTaskMessage>();
    app.genome_prepare_task = Some(GenomePrepareTask {
        job_id: 7,
        started: Instant::now(),
        cancel_requested: Arc::new(AtomicBool::new(false)),
        timeout_seconds: None,
        mode: GenomePrepareLaunchMode::Prepare,
        genome_id: "ToyGenome".to_string(),
        scope: GenomeDialogScope::Reference,
        catalog_path: "assets/genomes.json".to_string(),
        cache_dir: "data/genomes".to_string(),
        receiver: rx,
    });

    tx.send(GenomePrepareTaskMessage::Done {
        job_id: 6,
        result: Err(EngineError {
            code: ErrorCode::Internal,
            message: "stale".to_string(),

            cause_chain: vec![],
        }),
    })
    .unwrap();
    tx.send(GenomePrepareTaskMessage::Done {
        job_id: 7,
        result: Err(EngineError {
            code: ErrorCode::Internal,
            message: "actual".to_string(),

            cause_chain: vec![],
        }),
    })
    .unwrap();

    app.poll_prepare_reference_genome_task(&egui::Context::default());

    assert!(app.genome_prepare_task.is_none());
    assert!(app.job_event_log.iter().any(|event| {
        event.kind == BackgroundJobKind::PrepareGenome
            && event.phase == BackgroundJobEventPhase::IgnoredStale
            && event.job_id == Some(6)
    }));
    assert!(app.job_event_log.iter().any(|event| {
        event.kind == BackgroundJobKind::PrepareGenome
            && event.phase == BackgroundJobEventPhase::Failed
            && event.job_id == Some(7)
    }));
}

#[test]
fn poll_prepare_uses_cancel_request_signal_for_failure_classification() {
    let mut app = GENtleApp::default();
    let (tx, rx) = mpsc::channel::<GenomePrepareTaskMessage>();
    app.genome_prepare_task = Some(GenomePrepareTask {
        job_id: 52,
        started: Instant::now(),
        cancel_requested: Arc::new(AtomicBool::new(true)),
        timeout_seconds: Some(600),
        mode: GenomePrepareLaunchMode::Prepare,
        genome_id: "ToyGenome".to_string(),
        scope: GenomeDialogScope::Reference,
        catalog_path: "assets/genomes.json".to_string(),
        cache_dir: "data/genomes".to_string(),
        receiver: rx,
    });
    tx.send(GenomePrepareTaskMessage::Done {
        job_id: 52,
        result: Err(EngineError {
            code: ErrorCode::Internal,
            message: "worker interrupted".to_string(),

            cause_chain: vec![],
        }),
    })
    .expect("send prepare done");

    app.poll_prepare_reference_genome_task(&egui::Context::default());

    assert!(app.genome_prepare_task.is_none());
    assert!(
        app.genome_prepare_status
            .contains("Prepare genome cancelled after"),
        "status was: {}",
        app.genome_prepare_status
    );
}

#[test]
fn poll_prepare_timebox_overrides_cancel_wording_for_worker_failure() {
    let mut app = GENtleApp::default();
    let (tx, rx) = mpsc::channel::<GenomePrepareTaskMessage>();
    app.genome_prepare_task = Some(GenomePrepareTask {
        job_id: 53,
        started: Instant::now() - Duration::from_secs(3),
        cancel_requested: Arc::new(AtomicBool::new(false)),
        timeout_seconds: Some(1),
        mode: GenomePrepareLaunchMode::Prepare,
        genome_id: "ToyGenome".to_string(),
        scope: GenomeDialogScope::Reference,
        catalog_path: "assets/genomes.json".to_string(),
        cache_dir: "data/genomes".to_string(),
        receiver: rx,
    });
    tx.send(GenomePrepareTaskMessage::Done {
        job_id: 53,
        result: Err(EngineError {
            code: ErrorCode::Io,
            message: "Genome preparation cancelled for 'toy_genome'".to_string(),

            cause_chain: vec![],
        }),
    })
    .expect("send prepare done");

    app.poll_prepare_reference_genome_task(&egui::Context::default());

    assert!(app.genome_prepare_task.is_none());
    assert!(
        app.genome_prepare_status
            .contains("Prepare genome timed out after"),
        "status was: {}",
        app.genome_prepare_status
    );
}

#[test]
fn poll_prepare_captures_reinstall_recovery_for_inconsistent_reindex_failure() {
    let mut app = GENtleApp::default();
    let (tx, rx) = mpsc::channel::<GenomePrepareTaskMessage>();
    app.genome_prepare_task = Some(GenomePrepareTask {
        job_id: 153,
        started: Instant::now(),
        cancel_requested: Arc::new(AtomicBool::new(false)),
        timeout_seconds: None,
        mode: GenomePrepareLaunchMode::ReindexCachedFiles,
        genome_id: "Caenorhabditis elegans WBcel235 Ensembl 115".to_string(),
        scope: GenomeDialogScope::Reference,
        catalog_path: "assets/genomes.json".to_string(),
        cache_dir: "assets/data/genomes".to_string(),
        receiver: rx,
    });
    tx.send(GenomePrepareTaskMessage::Done {
            job_id: 153,
            result: Err(EngineError {
                code: ErrorCode::InvalidInput,
                message: "Prepared genome 'Caenorhabditis elegans WBcel235 Ensembl 115' is inconsistent: annotation gene index '/tmp/genes.json' references contigs missing from prepared sequence '/tmp/sequence.fa'".to_string(),
                cause_chain: vec![],
            }),
        })
        .expect("send prepare failure");

    app.poll_prepare_reference_genome_task(&egui::Context::default());

    let recovery = app
        .genome_prepare_failure_recovery
        .as_ref()
        .expect("failure should suggest reinstall");
    assert_eq!(
        recovery.genome_id,
        "Caenorhabditis elegans WBcel235 Ensembl 115"
    );
    assert_eq!(recovery.cache_dir, "assets/data/genomes");
}

#[test]
fn queue_prepare_failure_reinstall_uses_failed_job_settings() {
    let mut app = GENtleApp::default();
    app.genome_prepare_failure_recovery = Some(PrepareGenomeFailureRecovery {
        genome_id: "ToyGenome".to_string(),
        scope: GenomeDialogScope::Helper,
        catalog_path: "/tmp/helper_catalog.json".to_string(),
        cache_dir: "/tmp/helper_cache".to_string(),
    });

    app.queue_prepare_failure_reinstall(PreparedGenomeReinstallDialogHost::Root);

    let request = app
        .pending_prepared_genome_reinstall
        .as_ref()
        .expect("reinstall request");
    assert_eq!(request.genome_id, "ToyGenome");
    assert_eq!(request.scope, GenomeDialogScope::Helper);
    assert_eq!(request.catalog_path, "/tmp/helper_catalog.json");
    assert_eq!(request.cache_dir, "/tmp/helper_cache");
    assert_eq!(request.dialog_host, PreparedGenomeReinstallDialogHost::Root);
}

#[test]
fn poll_prepare_success_after_cancel_request_reports_completion_prefix() {
    let mut app = GENtleApp::default();
    let (tx, rx) = mpsc::channel::<GenomePrepareTaskMessage>();
    app.genome_prepare_task = Some(GenomePrepareTask {
        job_id: 54,
        started: Instant::now(),
        cancel_requested: Arc::new(AtomicBool::new(true)),
        timeout_seconds: None,
        mode: GenomePrepareLaunchMode::Prepare,
        genome_id: "ToyGenome".to_string(),
        scope: GenomeDialogScope::Reference,
        catalog_path: "assets/genomes.json".to_string(),
        cache_dir: "data/genomes".to_string(),
        receiver: rx,
    });
    tx.send(GenomePrepareTaskMessage::Done {
        job_id: 54,
        result: Ok(OpResult {
            op_id: "background-prepare-genome".to_string(),
            created_seq_ids: vec![],
            changed_seq_ids: vec![],
            warnings: vec![],
            messages: vec!["prepare completed quickly".to_string()],
            protocol_cartoon_preview: None,
            genome_annotation_projection: None,
            sequence_alignment: None,
            protein_derivation_report: None,
            reverse_translation_report: None,
            protease_digest_report: None,
            protein_residue_genomic_coordinates: None,
            exon_skip_selection_plan: None,
            exon_skip_materialization: None,
            cdna_assay_test_report: None,
            cdna_assay_product_materialization: None,
            transcript_qpcr_panel: None,
            primer_specificity_report: None,
            construct_reasoning_graph: None,
            sequencing_confirmation_report: None,
            sequencing_primer_overlay_report: None,
            sequencing_trace_import_report: None,
            sequencing_trace_record: None,
            sequencing_trace_summaries: None,
            cutrun_dataset_list: None,
            cutrun_dataset_status: None,
            cutrun_dataset_projection: None,
            cutrun_read_report: None,
            cutrun_read_report_summaries: None,
            cutrun_read_coverage_export: None,
            cutrun_regulatory_support: None,
            gene_set_resolution: None,
            gene_set_promoter_cohort: None,
            gene_set_cutrun_regulatory_support: None,
            read_acquisition_report: None,
            microarray_projection: None,
            probe_region_evidence_interpretation: None,
            genome_coordinate_projection: None,
            rna_read_gene_support_summary: None,
            rna_read_gene_support_audit: None,
            rna_read_target_quality_export: None,
            rna_read_batch_map_report: None,
            rna_read_isoform_preflight: None,
            tfbs_region_summary: None,
            tfbs_score_tracks: None,
            tfbs_track_similarity: None,
            multi_gene_promoter_tfbs: None,
            repeat_annotation_query: None,
            sequence_repeat_overlaps: None,
            repeat_feature_materialization: None,
            repeat_environment_cohort: None,
            window_cohort_tfbs: None,
            tfbs_hit_scan: None,
            restriction_site_scan: None,
            jaspar_remote_metadata_snapshot: None,
            jaspar_catalog_report: None,
            tf_query_resolution_report: None,
            jaspar_entry_expert_view: None,
            jaspar_registry_benchmark: None,
            jaspar_entry_presentation: None,
            sequence_context_view: None,
            sequence_context_bundle: None,
            alternative_promoter_comparison: None,
            variant_promoter_context: None,
            promoter_evidence_matrix: None,
            isoform_promoter_comparison: None,
            promoter_expression_evidence: None,
            promoter_artifact_manifest: None,
            promoter_reporter_candidates: None,
            reporter_catalog: None,
            reporter_recommendation: None,
            reporter_corpus_export: None,
            reporter_construct_handoff: None,
            uniprot_projection_audit: None,
            uniprot_projection_audit_parity: None,
            lab_assistant_instructions: None,
        }),
    })
    .expect("send prepare done");

    app.poll_prepare_reference_genome_task(&egui::Context::default());

    assert!(app.genome_prepare_task.is_none());
    assert!(
        app.genome_prepare_status
            .contains("Prepare genome finished after cancellation request"),
        "status was: {}",
        app.genome_prepare_status
    );
    assert!(app.job_event_log.iter().any(|event| {
        event.kind == BackgroundJobKind::PrepareGenome
            && event.phase == BackgroundJobEventPhase::Completed
            && event.job_id == Some(54)
            && event
                .summary
                .contains("finished after cancellation request")
    }));
}

#[test]
fn refresh_sequence_windows_for_seq_ids_targets_matching_windows_only() {
    let mut app = make_test_app_with_open_windows(&["seq_a", "seq_b"]);
    let refreshed = app.refresh_sequence_windows_for_seq_ids(&["seq_b".to_string()]);
    assert_eq!(refreshed, 1);
}

#[test]
fn poll_track_import_refreshes_only_changed_sequence_windows() {
    let mut app = make_test_app_with_open_windows(&["seq_a", "seq_b"]);
    let (tx, rx) = mpsc::channel::<GenomeTrackImportTaskMessage>();
    app.genome_track_import_task = Some(GenomeTrackImportTask {
        job_id: 91,
        started: Instant::now(),
        cancel_requested: Arc::new(AtomicBool::new(false)),
        receiver: rx,
    });
    tx.send(GenomeTrackImportTaskMessage::Done {
        job_id: 91,
        result: Ok(OpResult {
            op_id: "op_track_refresh_changed".to_string(),
            created_seq_ids: vec![],
            changed_seq_ids: vec!["seq_b".to_string()],
            warnings: vec![],
            messages: vec!["Imported 5 BED feature(s)".to_string()],
            protocol_cartoon_preview: None,
            genome_annotation_projection: None,
            sequence_alignment: None,
            protein_derivation_report: None,
            reverse_translation_report: None,
            protease_digest_report: None,
            protein_residue_genomic_coordinates: None,
            exon_skip_selection_plan: None,
            exon_skip_materialization: None,
            cdna_assay_test_report: None,
            cdna_assay_product_materialization: None,
            transcript_qpcr_panel: None,
            primer_specificity_report: None,
            construct_reasoning_graph: None,
            sequencing_confirmation_report: None,
            sequencing_primer_overlay_report: None,
            sequencing_trace_import_report: None,
            sequencing_trace_record: None,
            sequencing_trace_summaries: None,
            cutrun_dataset_list: None,
            cutrun_dataset_status: None,
            cutrun_dataset_projection: None,
            cutrun_read_report: None,
            cutrun_read_report_summaries: None,
            cutrun_read_coverage_export: None,
            cutrun_regulatory_support: None,
            gene_set_resolution: None,
            gene_set_promoter_cohort: None,
            gene_set_cutrun_regulatory_support: None,
            read_acquisition_report: None,
            microarray_projection: None,
            probe_region_evidence_interpretation: None,
            genome_coordinate_projection: None,
            rna_read_gene_support_summary: None,
            rna_read_gene_support_audit: None,
            rna_read_target_quality_export: None,
            rna_read_batch_map_report: None,
            rna_read_isoform_preflight: None,
            tfbs_region_summary: None,
            tfbs_score_tracks: None,
            tfbs_track_similarity: None,
            multi_gene_promoter_tfbs: None,
            repeat_annotation_query: None,
            sequence_repeat_overlaps: None,
            repeat_feature_materialization: None,
            repeat_environment_cohort: None,
            window_cohort_tfbs: None,
            tfbs_hit_scan: None,
            restriction_site_scan: None,
            jaspar_remote_metadata_snapshot: None,
            jaspar_catalog_report: None,
            tf_query_resolution_report: None,
            jaspar_entry_expert_view: None,
            jaspar_registry_benchmark: None,
            jaspar_entry_presentation: None,
            sequence_context_view: None,
            sequence_context_bundle: None,
            alternative_promoter_comparison: None,
            variant_promoter_context: None,
            promoter_evidence_matrix: None,
            isoform_promoter_comparison: None,
            promoter_expression_evidence: None,
            promoter_artifact_manifest: None,
            promoter_reporter_candidates: None,
            reporter_catalog: None,
            reporter_recommendation: None,
            reporter_corpus_export: None,
            reporter_construct_handoff: None,
            uniprot_projection_audit: None,
            uniprot_projection_audit_parity: None,
            lab_assistant_instructions: None,
        }),
    })
    .expect("send track import done");

    app.poll_genome_track_import_task(&egui::Context::default());

    assert!(app.genome_track_import_task.is_none());
    assert!(
        app.genome_track_status
            .contains("refreshed 1 open sequence window(s)"),
        "status was: {}",
        app.genome_track_status
    );
}

#[test]
fn poll_track_import_refreshes_all_open_windows_when_changed_ids_missing() {
    let mut app = make_test_app_with_open_windows(&["seq_a", "seq_b"]);
    let (tx, rx) = mpsc::channel::<GenomeTrackImportTaskMessage>();
    app.genome_track_import_task = Some(GenomeTrackImportTask {
        job_id: 92,
        started: Instant::now(),
        cancel_requested: Arc::new(AtomicBool::new(false)),
        receiver: rx,
    });
    tx.send(GenomeTrackImportTaskMessage::Done {
        job_id: 92,
        result: Ok(OpResult {
            op_id: "op_track_refresh_fallback".to_string(),
            created_seq_ids: vec![],
            changed_seq_ids: vec![],
            warnings: vec![],
            messages: vec!["Imported 5 BED feature(s)".to_string()],
            protocol_cartoon_preview: None,
            genome_annotation_projection: None,
            sequence_alignment: None,
            protein_derivation_report: None,
            reverse_translation_report: None,
            protease_digest_report: None,
            protein_residue_genomic_coordinates: None,
            exon_skip_selection_plan: None,
            exon_skip_materialization: None,
            cdna_assay_test_report: None,
            cdna_assay_product_materialization: None,
            transcript_qpcr_panel: None,
            primer_specificity_report: None,
            construct_reasoning_graph: None,
            sequencing_confirmation_report: None,
            sequencing_primer_overlay_report: None,
            sequencing_trace_import_report: None,
            sequencing_trace_record: None,
            sequencing_trace_summaries: None,
            cutrun_dataset_list: None,
            cutrun_dataset_status: None,
            cutrun_dataset_projection: None,
            cutrun_read_report: None,
            cutrun_read_report_summaries: None,
            cutrun_read_coverage_export: None,
            cutrun_regulatory_support: None,
            gene_set_resolution: None,
            gene_set_promoter_cohort: None,
            gene_set_cutrun_regulatory_support: None,
            read_acquisition_report: None,
            microarray_projection: None,
            probe_region_evidence_interpretation: None,
            genome_coordinate_projection: None,
            rna_read_gene_support_summary: None,
            rna_read_gene_support_audit: None,
            rna_read_target_quality_export: None,
            rna_read_batch_map_report: None,
            rna_read_isoform_preflight: None,
            tfbs_region_summary: None,
            tfbs_score_tracks: None,
            tfbs_track_similarity: None,
            multi_gene_promoter_tfbs: None,
            repeat_annotation_query: None,
            sequence_repeat_overlaps: None,
            repeat_feature_materialization: None,
            repeat_environment_cohort: None,
            window_cohort_tfbs: None,
            tfbs_hit_scan: None,
            restriction_site_scan: None,
            jaspar_remote_metadata_snapshot: None,
            jaspar_catalog_report: None,
            tf_query_resolution_report: None,
            jaspar_entry_expert_view: None,
            jaspar_registry_benchmark: None,
            jaspar_entry_presentation: None,
            sequence_context_view: None,
            sequence_context_bundle: None,
            alternative_promoter_comparison: None,
            variant_promoter_context: None,
            promoter_evidence_matrix: None,
            isoform_promoter_comparison: None,
            promoter_expression_evidence: None,
            promoter_artifact_manifest: None,
            promoter_reporter_candidates: None,
            reporter_catalog: None,
            reporter_recommendation: None,
            reporter_corpus_export: None,
            reporter_construct_handoff: None,
            uniprot_projection_audit: None,
            uniprot_projection_audit_parity: None,
            lab_assistant_instructions: None,
        }),
    })
    .expect("send track import done");

    app.poll_genome_track_import_task(&egui::Context::default());

    assert!(app.genome_track_import_task.is_none());
    assert!(
        app.genome_track_status
            .contains("refreshed 2 open sequence window(s)"),
        "status was: {}",
        app.genome_track_status
    );
}

#[test]
fn format_extract_region_status_includes_annotation_fallback_reason() {
    let status = GENtleApp::format_extract_region_status(&OpResult {
        op_id: "op_extract".to_string(),
        created_seq_ids: vec!["grch38_tp73".to_string()],
        changed_seq_ids: vec![],
        warnings: vec![],
        messages: vec!["Extracted region".to_string()],
        protocol_cartoon_preview: None,
        genome_annotation_projection: Some(GenomeAnnotationProjectionTelemetry {
            requested_scope: "full".to_string(),
            effective_scope: "core".to_string(),
            max_features_cap: Some(500),
            candidate_feature_count: 1400,
            attached_feature_count: 480,
            dropped_feature_count: 920,
            genes_attached: 12,
            transcripts_attached: 26,
            exons_attached: 420,
            cds_attached: 22,
            fallback_applied: true,
            fallback_reason: Some(
                "Projected full annotation exceeded cap and fell back to core".to_string(),
            ),
        }),
        sequence_alignment: None,
        protein_derivation_report: None,
        reverse_translation_report: None,
        protease_digest_report: None,
        protein_residue_genomic_coordinates: None,
        exon_skip_selection_plan: None,
        exon_skip_materialization: None,
        cdna_assay_test_report: None,
        cdna_assay_product_materialization: None,
        transcript_qpcr_panel: None,
        primer_specificity_report: None,
        construct_reasoning_graph: None,
        sequencing_confirmation_report: None,
        sequencing_primer_overlay_report: None,
        sequencing_trace_import_report: None,
        sequencing_trace_record: None,
        sequencing_trace_summaries: None,
        cutrun_dataset_list: None,
        cutrun_dataset_status: None,
        cutrun_dataset_projection: None,
        cutrun_read_report: None,
        cutrun_read_report_summaries: None,
        cutrun_read_coverage_export: None,
        cutrun_regulatory_support: None,
        gene_set_resolution: None,
        gene_set_promoter_cohort: None,
        gene_set_cutrun_regulatory_support: None,
        read_acquisition_report: None,
        microarray_projection: None,
        probe_region_evidence_interpretation: None,
        genome_coordinate_projection: None,
        rna_read_gene_support_summary: None,
        rna_read_gene_support_audit: None,
        rna_read_target_quality_export: None,
        rna_read_batch_map_report: None,
        rna_read_isoform_preflight: None,
        tfbs_region_summary: None,
        tfbs_score_tracks: None,
        tfbs_track_similarity: None,
        multi_gene_promoter_tfbs: None,
        repeat_annotation_query: None,
        sequence_repeat_overlaps: None,
        repeat_feature_materialization: None,
        repeat_environment_cohort: None,
        window_cohort_tfbs: None,
        tfbs_hit_scan: None,
        restriction_site_scan: None,
        jaspar_remote_metadata_snapshot: None,
        jaspar_catalog_report: None,
        tf_query_resolution_report: None,
        jaspar_entry_expert_view: None,
        jaspar_registry_benchmark: None,
        jaspar_entry_presentation: None,
        sequence_context_view: None,
        sequence_context_bundle: None,
        alternative_promoter_comparison: None,
        variant_promoter_context: None,
        promoter_evidence_matrix: None,
        isoform_promoter_comparison: None,
        promoter_expression_evidence: None,
        promoter_artifact_manifest: None,
        promoter_reporter_candidates: None,
        reporter_catalog: None,
        reporter_recommendation: None,
        reporter_corpus_export: None,
        reporter_construct_handoff: None,
        uniprot_projection_audit: None,
        uniprot_projection_audit_parity: None,
        lab_assistant_instructions: None,
    });
    assert!(status.contains("annotation: requested=full effective=core"));
    assert!(status.contains("annotation kinds: genes=12 transcripts=26 exons=420 cds=22"));
    assert!(status.contains(
        "annotation fallback reason: Projected full annotation exceeded cap and fell back to core"
    ));
}

#[test]
fn sync_genome_annotation_scope_controls_keeps_checkbox_and_scope_consistent() {
    let mut app = GENtleApp::default();
    app.genome_include_genomic_annotation = false;
    app.genome_annotation_scope = crate::engine::GenomeAnnotationScope::Full;
    app.sync_genome_annotation_scope_controls();
    assert_eq!(
        app.genome_annotation_scope,
        crate::engine::GenomeAnnotationScope::None
    );

    app.genome_include_genomic_annotation = true;
    app.genome_annotation_scope = crate::engine::GenomeAnnotationScope::None;
    app.sync_genome_annotation_scope_controls();
    assert_eq!(
        app.genome_annotation_scope,
        crate::engine::GenomeAnnotationScope::Core
    );
}

#[test]
fn select_gene_record_sets_default_output_id_with_genome_and_interval() {
    let mut app = GENtleApp::default();
    app.genome_id = "Human GRCh38 Ensembl 116".to_string();
    app.genome_genes = vec![crate::genomes::GenomeGeneRecord {
        chromosome: "1".to_string(),
        start_1based: 3652307,
        end_1based: 3736201,
        strand: Some('+'),
        gene_id: Some("ENSG00000078900".to_string()),
        gene_name: Some("TP73".to_string()),
        biotype: Some("protein_coding".to_string()),
    }];

    app.select_gene_record(0);

    assert_eq!(
        app.genome_output_id,
        "human_grch38_ensembl_116_tp73_3652307_3736201"
    );
}

#[test]
fn select_gene_record_keeps_existing_output_id() {
    let mut app = GENtleApp::default();
    app.genome_id = "Human GRCh38 Ensembl 116".to_string();
    app.genome_output_id = "manual_output".to_string();
    app.genome_genes = vec![crate::genomes::GenomeGeneRecord {
        chromosome: "1".to_string(),
        start_1based: 3652307,
        end_1based: 3736201,
        strand: Some('+'),
        gene_id: Some("ENSG00000078900".to_string()),
        gene_name: Some("TP73".to_string()),
        biotype: Some("protein_coding".to_string()),
    }];

    app.select_gene_record(0);

    assert_eq!(app.genome_output_id, "manual_output");
}

#[test]
fn select_gene_record_uses_coding_promoter_output_id_when_mode_is_enabled() {
    let mut app = GENtleApp::default();
    app.genome_id = "Human GRCh38 Ensembl 116".to_string();
    app.genome_gene_extract_mode = GenomeGeneExtractMode::CodingWithPromoter;
    app.genome_gene_promoter_upstream_bp = "2000".to_string();
    app.genome_genes = vec![crate::genomes::GenomeGeneRecord {
        chromosome: "20".to_string(),
        start_1based: 33756336,
        end_1based: 33766457,
        strand: Some('+'),
        gene_id: Some("ENSG00000101412".to_string()),
        gene_name: Some("E2F1".to_string()),
        biotype: Some("protein_coding".to_string()),
    }];

    app.select_gene_record(0);

    assert_eq!(
        app.genome_output_id,
        "human_grch38_ensembl_116_e2f1_coding_promoter_2000bp"
    );
}

#[test]
fn poll_blast_task_updates_runtime_status_from_worker_messages() {
    let mut app = GENtleApp::default();
    let (tx, rx) = mpsc::channel::<GenomeBlastTaskMessage>();
    app.genome_blast_task = Some(GenomeBlastTask {
        job_id: 71,
        started: Instant::now(),
        cancel_requested: Arc::new(AtomicBool::new(false)),
        receiver: rx,
    });
    tx.send(GenomeBlastTaskMessage::Status {
        job_id: 71,
        status: "BLAST query 'q1' running (0.5s)".to_string(),
    })
    .expect("send status");

    app.poll_reference_genome_blast_task(&egui::Context::default());

    assert!(
        app.genome_blast_status.contains("running"),
        "status was: {}",
        app.genome_blast_status
    );
    assert!(app.genome_blast_task.is_some());
}

#[test]
fn blast_options_preset_payloads_are_json_objects() {
    let presets = [
        GenomeBlastOptionsPreset::StrictIdentityCoverage,
        GenomeBlastOptionsPreset::UniqueBestHit,
        GenomeBlastOptionsPreset::HighStringency,
    ];
    for preset in presets {
        let payload = preset
            .as_request_override_json()
            .expect("non-empty preset payload");
        assert!(payload.is_object(), "preset {} payload", preset.label());
    }
}

#[test]
fn build_genome_blast_request_override_merges_preset_and_advanced_json() {
    let mut app = GENtleApp::default();
    app.genome_blast_options_preset = GenomeBlastOptionsPreset::StrictIdentityCoverage;
    app.genome_blast_options_json =
        r#"{"thresholds":{"min_query_coverage_percent":85.0},"max_hits":9}"#.to_string();
    let merged = app
        .build_genome_blast_request_override_json()
        .expect("merge override")
        .expect("override object");
    assert_eq!(merged.get("max_hits").and_then(|v| v.as_u64()), Some(9));
    let thresholds = merged.get("thresholds").expect("thresholds object");
    assert_eq!(
        thresholds
            .get("min_identity_percent")
            .and_then(|v| v.as_f64()),
        Some(97.0)
    );
    assert_eq!(
        thresholds
            .get("min_query_coverage_percent")
            .and_then(|v| v.as_f64()),
        Some(85.0)
    );
}

#[test]
fn build_genome_blast_request_override_rejects_non_object_json() {
    let mut app = GENtleApp::default();
    app.genome_blast_options_json = "[1,2,3]".to_string();
    let err = app
        .build_genome_blast_request_override_json()
        .expect_err("non-object should fail");
    assert!(err.contains("JSON object"), "{err}");
}

#[test]
fn build_genome_blast_request_override_includes_structured_thresholds() {
    let mut app = GENtleApp::default();
    app.genome_blast_threshold_use_min_identity_percent = true;
    app.genome_blast_threshold_min_identity_percent = "98.5".to_string();
    app.genome_blast_threshold_use_min_query_coverage_percent = true;
    app.genome_blast_threshold_min_query_coverage_percent = "82.0".to_string();
    app.genome_blast_threshold_unique_best_hit = true;
    let merged = app
        .build_genome_blast_request_override_json()
        .expect("merge override")
        .expect("override object");
    let thresholds = merged.get("thresholds").expect("thresholds object");
    assert_eq!(
        thresholds
            .get("min_identity_percent")
            .and_then(|v| v.as_f64()),
        Some(98.5)
    );
    assert_eq!(
        thresholds
            .get("min_query_coverage_percent")
            .and_then(|v| v.as_f64()),
        Some(82.0)
    );
    assert_eq!(
        thresholds.get("unique_best_hit").and_then(|v| v.as_bool()),
        Some(true)
    );
}

#[test]
fn build_genome_blast_request_override_rejects_invalid_structured_threshold_value() {
    let mut app = GENtleApp::default();
    app.genome_blast_threshold_use_min_identity_percent = true;
    app.genome_blast_threshold_min_identity_percent = "abc".to_string();
    let err = app
        .build_genome_blast_request_override_json()
        .expect_err("invalid threshold should fail");
    assert!(err.contains("min_identity_percent"), "{err}");
}

fn make_lineage_row(node_id: &str, seq_id: &str) -> LineageRow {
    LineageRow {
        kind: LineageNodeKind::Sequence,
        node_id: node_id.to_string(),
        seq_id: seq_id.to_string(),
        display_name: seq_id.to_string(),
        origin: "Test".to_string(),
        created_by_op: "op_test".to_string(),
        created_at: 1,
        parents: vec![],
        length: 100,
        circular: false,
        pool_size: 1,
        pool_members: vec![seq_id.to_string()],
        arrangement_id: None,
        arrangement_mode: None,
        lane_container_ids: vec![],
        ladders: vec![],
        genome_anchor_summary: None,
        genome_anchor_display: None,
        is_full_genome_sequence: false,
        retrieval_descriptor: None,
        analysis_kind: None,
        analysis_artifact_id: None,
        analysis_reference_seq_id: None,
        analysis_mode: None,
        analysis_status: None,
        analysis_point_count: None,
        analysis_bin_count: None,
        analysis_read_count: None,
        analysis_trace_count: None,
        analysis_target_count: None,
        analysis_variant_count: None,
        macro_instance_id: None,
        macro_routine_id: None,
        macro_template_name: None,
        macro_status: None,
        macro_status_message: None,
        macro_op_ids: vec![],
        macro_inputs: vec![],
        macro_outputs: vec![],
    }
}

#[test]
fn lineage_copy_primary_id_prefers_analysis_artifact_id() {
    let mut row = make_lineage_row("analysis:qpcr:tp73_as3", "tp73_as3_template");
    row.kind = LineageNodeKind::Analysis;
    row.analysis_kind = Some(LineageAnalysisKind::QpcrDesign);
    row.analysis_artifact_id = Some("qpcr_report_tp73_as3".to_string());

    assert_eq!(
        GENtleApp::lineage_copy_payload(&row, LineageCopyPayloadKind::PrimaryId),
        "qpcr_report_tp73_as3"
    );
}

#[test]
fn lineage_copy_display_label_falls_back_to_primary_id() {
    let mut row = make_lineage_row("node:seq1", "seq1");
    row.display_name = "  ".to_string();

    assert_eq!(
        GENtleApp::lineage_copy_payload(&row, LineageCopyPayloadKind::DisplayLabel),
        "seq1"
    );
}

#[test]
fn lineage_copy_row_summary_includes_identifiers_and_analysis_context() {
    let mut row = make_lineage_row("analysis:primer:tp73", "tp73_template");
    row.kind = LineageNodeKind::Analysis;
    row.display_name = "TP73 primer report".to_string();
    row.parents = vec!["node:template".to_string()];
    row.analysis_kind = Some(LineageAnalysisKind::PrimerDesign);
    row.analysis_artifact_id = Some("primer_report_tp73".to_string());
    row.analysis_reference_seq_id = Some("tp73_reference".to_string());
    row.analysis_mode = Some("oligo".to_string());
    row.analysis_target_count = Some(3);

    let summary = GENtleApp::lineage_copy_payload(&row, LineageCopyPayloadKind::RowSummary);

    assert!(
        summary.contains("node_id: analysis:primer:tp73"),
        "{summary}"
    );
    assert!(
        summary.contains("primary_id: primer_report_tp73"),
        "{summary}"
    );
    assert!(
        summary.contains("display_name: TP73 primer report"),
        "{summary}"
    );
    assert!(summary.contains("parents: node:template"), "{summary}");
    assert!(
        summary.contains("analysis_kind: primer_design"),
        "{summary}"
    );
    assert!(
        summary.contains("analysis_reference_seq_id: tp73_reference"),
        "{summary}"
    );
    assert!(summary.contains("analysis_target_count: 3"), "{summary}");
}

fn insert_test_lineage_node(state: &mut ProjectState, node_id: &str, seq_id: &str) {
    state.lineage.nodes.insert(
        node_id.to_string(),
        LineageNode {
            node_id: node_id.to_string(),
            seq_id: seq_id.to_string(),
            origin: SequenceOrigin::ImportedUnknown,
            created_by_op: None,
            created_at_unix_ms: 0,
        },
    );
    state
        .lineage
        .seq_to_node
        .insert(seq_id.to_string(), node_id.to_string());
}

#[test]
fn lineage_analysis_open_payload_infers_missing_metadata_from_node_id() {
    let mut dotplot_row = make_lineage_row("analysis:dotplot:tp73_dp", "seq_dotplot");
    dotplot_row.kind = LineageNodeKind::Analysis;
    dotplot_row.display_name.clear();
    dotplot_row.analysis_kind = None;
    dotplot_row.analysis_artifact_id = None;
    assert_eq!(
        GENtleApp::lineage_analysis_open_payload(&dotplot_row),
        Some((
            LineageAnalysisKind::Dotplot,
            "seq_dotplot".to_string(),
            "tp73_dp".to_string(),
        ))
    );

    let mut flex_row = make_lineage_row("analysis:flex:tp73_fx", "seq_flex");
    flex_row.kind = LineageNodeKind::Analysis;
    flex_row.display_name.clear();
    flex_row.analysis_kind = None;
    flex_row.analysis_artifact_id = None;
    assert_eq!(
        GENtleApp::lineage_analysis_open_payload(&flex_row),
        Some((
            LineageAnalysisKind::FlexibilityTrack,
            "seq_flex".to_string(),
            "tp73_fx".to_string(),
        ))
    );

    let mut rna_row = make_lineage_row("analysis:rna_reads:tp73_rna", "seq_rna");
    rna_row.kind = LineageNodeKind::Analysis;
    rna_row.display_name.clear();
    rna_row.analysis_kind = None;
    rna_row.analysis_artifact_id = None;
    assert_eq!(
        GENtleApp::lineage_analysis_open_payload(&rna_row),
        Some((
            LineageAnalysisKind::RnaReadInterpretation,
            "seq_rna".to_string(),
            "tp73_rna".to_string(),
        ))
    );

    let mut primer_row = make_lineage_row("analysis:primer:tp73_primer", "seq_primer");
    primer_row.kind = LineageNodeKind::Analysis;
    primer_row.display_name.clear();
    primer_row.analysis_kind = None;
    primer_row.analysis_artifact_id = None;
    assert_eq!(
        GENtleApp::lineage_analysis_open_payload(&primer_row),
        Some((
            LineageAnalysisKind::PrimerDesign,
            "seq_primer".to_string(),
            "tp73_primer".to_string(),
        ))
    );

    let mut qpcr_row = make_lineage_row("analysis:qpcr:tp73_qpcr", "seq_qpcr");
    qpcr_row.kind = LineageNodeKind::Analysis;
    qpcr_row.display_name.clear();
    qpcr_row.analysis_kind = None;
    qpcr_row.analysis_artifact_id = None;
    assert_eq!(
        GENtleApp::lineage_analysis_open_payload(&qpcr_row),
        Some((
            LineageAnalysisKind::QpcrDesign,
            "seq_qpcr".to_string(),
            "tp73_qpcr".to_string(),
        ))
    );

    let mut restriction_handoff_row = make_lineage_row(
        "analysis:restriction_cloning_pcr:tp73_clone_handoff",
        "seq_clone",
    );
    restriction_handoff_row.kind = LineageNodeKind::Analysis;
    restriction_handoff_row.display_name.clear();
    restriction_handoff_row.analysis_kind = None;
    restriction_handoff_row.analysis_artifact_id = None;
    assert_eq!(
        GENtleApp::lineage_analysis_open_payload(&restriction_handoff_row),
        Some((
            LineageAnalysisKind::RestrictionCloningPcrHandoff,
            "seq_clone".to_string(),
            "tp73_clone_handoff".to_string(),
        ))
    );

    let mut protein_row = make_lineage_row("analysis:protein_derive:tp73_protein", "seq_protein");
    protein_row.kind = LineageNodeKind::Analysis;
    protein_row.display_name.clear();
    protein_row.analysis_kind = None;
    protein_row.analysis_artifact_id = None;
    assert_eq!(
        GENtleApp::lineage_analysis_open_payload(&protein_row),
        Some((
            LineageAnalysisKind::ProteinDerivation,
            "seq_protein".to_string(),
            "tp73_protein".to_string(),
        ))
    );

    let mut reverse_translation_row =
        make_lineage_row("analysis:reverse_translate:tp73_coding", "seq_reverse");
    reverse_translation_row.kind = LineageNodeKind::Analysis;
    reverse_translation_row.display_name.clear();
    reverse_translation_row.analysis_kind = None;
    reverse_translation_row.analysis_artifact_id = None;
    assert_eq!(
        GENtleApp::lineage_analysis_open_payload(&reverse_translation_row),
        Some((
            LineageAnalysisKind::ReverseTranslation,
            "seq_reverse".to_string(),
            "tp73_coding".to_string(),
        ))
    );

    let mut reasoning_row = make_lineage_row(
        "analysis:construct_reasoning:tp73_reasoning",
        "seq_reasoning",
    );
    reasoning_row.kind = LineageNodeKind::Analysis;
    reasoning_row.display_name.clear();
    reasoning_row.analysis_kind = None;
    reasoning_row.analysis_artifact_id = None;
    assert_eq!(
        GENtleApp::lineage_analysis_open_payload(&reasoning_row),
        Some((
            LineageAnalysisKind::ConstructReasoning,
            "seq_reasoning".to_string(),
            "tp73_reasoning".to_string(),
        ))
    );

    let mut uniprot_row = make_lineage_row("analysis:uniprot:tp53_uniprot_p04637", "seq_uniprot");
    uniprot_row.kind = LineageNodeKind::Analysis;
    uniprot_row.display_name.clear();
    uniprot_row.analysis_kind = None;
    uniprot_row.analysis_artifact_id = None;
    assert_eq!(
        GENtleApp::lineage_analysis_open_payload(&uniprot_row),
        Some((
            LineageAnalysisKind::UniprotProjection,
            "seq_uniprot".to_string(),
            "tp53_uniprot_p04637".to_string(),
        ))
    );

    let mut confirmation_row = make_lineage_row("analysis:seq_confirm:tp73_confirm", "seq_confirm");
    confirmation_row.kind = LineageNodeKind::Analysis;
    confirmation_row.display_name.clear();
    confirmation_row.analysis_kind = None;
    confirmation_row.analysis_artifact_id = None;
    assert_eq!(
        GENtleApp::lineage_analysis_open_payload(&confirmation_row),
        Some((
            LineageAnalysisKind::SequencingConfirmation,
            "seq_confirm".to_string(),
            "tp73_confirm".to_string(),
        ))
    );
}

#[test]
fn sanitize_lineage_node_groups_enforces_disjoint_membership() {
    let valid_node_ids = ["n1", "n2", "n3", "n4"]
        .iter()
        .map(|node_id| (*node_id).to_string())
        .collect();
    let groups = vec![
        PersistedLineageNodeGroup {
            group_id: "grp-a".to_string(),
            label: "A".to_string(),
            representative_node_id: "n1".to_string(),
            member_node_ids: vec!["n2".to_string(), "n5".to_string()],
            collapsed: false,
        },
        PersistedLineageNodeGroup {
            group_id: "grp-b".to_string(),
            label: "B".to_string(),
            representative_node_id: "n2".to_string(),
            member_node_ids: vec!["n3".to_string(), "n4".to_string()],
            collapsed: true,
        },
    ];

    let sanitized = GENtleApp::sanitize_lineage_node_groups(&groups, &valid_node_ids);
    assert_eq!(sanitized.len(), 2);
    assert_eq!(sanitized[0].representative_node_id, "n1");
    assert_eq!(sanitized[0].member_node_ids, vec!["n2".to_string()]);
    assert_eq!(sanitized[1].representative_node_id, "n3");
    assert_eq!(sanitized[1].member_node_ids, vec!["n4".to_string()]);
}

#[test]
fn project_lineage_graph_by_groups_collapses_member_edges_to_representative() {
    let rows = vec![
        make_lineage_row("n1", "seq1"),
        make_lineage_row("n2", "seq2"),
        make_lineage_row("n3", "seq3"),
    ];
    let edges = vec![
        ("n1".to_string(), "n2".to_string(), "op1".to_string()),
        ("n2".to_string(), "n3".to_string(), "op2".to_string()),
    ];
    let groups = vec![PersistedLineageNodeGroup {
        group_id: "grp-1".to_string(),
        label: "Collapsed group".to_string(),
        representative_node_id: "n1".to_string(),
        member_node_ids: vec!["n2".to_string()],
        collapsed: true,
    }];

    let (projected_rows, projected_edges) =
        GENtleApp::project_lineage_graph_by_groups(&rows, &edges, &groups);
    assert_eq!(projected_rows.len(), 2);
    assert!(projected_rows.iter().any(|row| row.node_id == "n1"));
    assert!(projected_rows.iter().any(|row| row.node_id == "n3"));
    assert_eq!(
        projected_edges,
        vec![("n1".to_string(), "n3".to_string(), "op2".to_string())]
    );
}

#[test]
fn project_lineage_graph_operation_hubs_inserts_one_gibson_box() {
    let rows = vec![
        make_lineage_row("n1", "dest"),
        make_lineage_row("n2", "insert"),
        make_lineage_row("n3", "primer_left"),
        make_lineage_row("n4", "primer_right"),
        make_lineage_row("n5", "product"),
    ];
    let edges = vec![
        ("n1".to_string(), "n3".to_string(), "op-gb".to_string()),
        ("n1".to_string(), "n4".to_string(), "op-gb".to_string()),
        ("n1".to_string(), "n5".to_string(), "op-gb".to_string()),
        ("n2".to_string(), "n3".to_string(), "op-gb".to_string()),
        ("n2".to_string(), "n4".to_string(), "op-gb".to_string()),
        ("n2".to_string(), "n5".to_string(), "op-gb".to_string()),
    ];
    let op_labels = HashMap::from([("op-gb".to_string(), "Gibson cloning".to_string())]);
    let hub_ops = HashSet::from(["op-gb".to_string()]);

    let (projected_rows, projected_edges, projected_labels) =
        GENtleApp::project_lineage_graph_operation_hubs(&rows, &edges, &op_labels, &hub_ops);

    let hub = projected_rows
        .iter()
        .find(|row| row.node_id == "operation:op-gb")
        .expect("hub row");
    assert_eq!(hub.origin, "OperationHub");
    assert_eq!(hub.display_name, "Gibson cloning");
    assert!(projected_edges.iter().all(|edge| edge.2 != "op-gb"));
    assert!(projected_edges.contains(&(
        "n1".to_string(),
        "operation:op-gb".to_string(),
        "op-gb::hub_in".to_string(),
    )));
    assert!(projected_edges.contains(&(
        "operation:op-gb".to_string(),
        "n5".to_string(),
        "op-gb::hub_out".to_string(),
    )));
    assert_eq!(projected_edges.len(), 5);
    assert_eq!(
        projected_labels.get("op-gb::hub_in").map(String::as_str),
        Some("Gibson cloning")
    );
}

#[test]
fn render_current_visible_lineage_svg_text_uses_gibson_operation_hub_projection() {
    let mut app = GENtleApp::default();
    app.current_project_path = Some("/tmp/demo.gentle.json".to_string());
    app.lineage_rows = vec![
        make_lineage_row("n1", "dest"),
        make_lineage_row("n2", "insert"),
        make_lineage_row("n3", "primer_left"),
        make_lineage_row("n4", "primer_right"),
        make_lineage_row("n5", "product"),
    ];
    app.lineage_edges = vec![
        ("n1".to_string(), "n3".to_string(), "op-gb".to_string()),
        ("n1".to_string(), "n4".to_string(), "op-gb".to_string()),
        ("n1".to_string(), "n5".to_string(), "op-gb".to_string()),
        ("n2".to_string(), "n3".to_string(), "op-gb".to_string()),
        ("n2".to_string(), "n4".to_string(), "op-gb".to_string()),
        ("n2".to_string(), "n5".to_string(), "op-gb".to_string()),
    ];
    app.lineage_op_label_by_id =
        HashMap::from([("op-gb".to_string(), "Gibson cloning".to_string())]);
    app.lineage_reopenable_gibson_op_ids = HashSet::from(["op-gb".to_string()]);
    app.lineage_cache_valid = true;
    app.lineage_cache_stamp = app.current_lineage_change_stamp();

    let svg = app.render_current_visible_lineage_svg_text();

    assert!(svg.contains("GENtle Lineage (DALG) - demo.gentle.json"));
    assert_eq!(svg.matches("Gibson cloning").count(), 1);
    assert!(!svg.contains("op=op-gb"));
}

#[test]
fn build_lineage_table_entries_indents_group_members() {
    let rows = vec![
        make_lineage_row("n1", "seq1"),
        make_lineage_row("n2", "seq2"),
        make_lineage_row("n3", "seq3"),
    ];
    let expanded_group = vec![PersistedLineageNodeGroup {
        group_id: "grp-1".to_string(),
        label: "Expanded".to_string(),
        representative_node_id: "n1".to_string(),
        member_node_ids: vec!["n2".to_string()],
        collapsed: false,
    }];
    let expanded_entries = GENtleApp::build_lineage_table_entries(&rows, &expanded_group);
    assert_eq!(expanded_entries.len(), 3);
    assert!(expanded_entries[0].is_group_representative);
    assert_eq!(expanded_entries[1].indent_level, 1);
    assert_eq!(expanded_entries[1].row.node_id, "n2");

    let collapsed_group = vec![PersistedLineageNodeGroup {
        collapsed: true,
        ..expanded_group[0].clone()
    }];
    let collapsed_entries = GENtleApp::build_lineage_table_entries(&rows, &collapsed_group);
    assert_eq!(collapsed_entries.len(), 2);
    assert_eq!(collapsed_entries[0].row.node_id, "n1");
    assert_eq!(collapsed_entries[0].hidden_group_member_count, 1);
    assert_eq!(collapsed_entries[1].row.node_id, "n3");
}

#[test]
fn lineage_leaf_node_ids_detects_terminal_nodes() {
    let rows = vec![
        make_lineage_row("n1", "seq1"),
        make_lineage_row("n2", "seq2"),
        make_lineage_row("n3", "seq3"),
    ];
    let edges = vec![
        ("n1".to_string(), "n2".to_string(), "op-1".to_string()),
        ("n2".to_string(), "n3".to_string(), "op-2".to_string()),
    ];
    let leaves = GENtleApp::lineage_leaf_node_ids(&rows, &edges);
    assert!(!leaves.contains("n1"));
    assert!(!leaves.contains("n2"));
    assert!(leaves.contains("n3"));
}

#[test]
fn rename_leaf_lineage_node_updates_sequence_name() {
    let mut app = GENtleApp::default();
    {
        let mut engine = app.engine.write().unwrap();
        let state = engine.state_mut();
        state.sequences.insert(
            "seq_a".to_string(),
            DNAsequence::from_sequence("ATGC").unwrap(),
        );
        state.sequences.insert(
            "seq_b".to_string(),
            DNAsequence::from_sequence("ATGC").unwrap(),
        );
        insert_test_lineage_node(state, "n1", "seq_a");
        insert_test_lineage_node(state, "n2", "seq_b");
        state.lineage.edges.push(LineageEdge {
            from_node_id: "n1".to_string(),
            to_node_id: "n2".to_string(),
            op_id: "op-1".to_string(),
            run_id: "interactive".to_string(),
        });
    }
    let status = app
        .rename_leaf_lineage_node("n2", "Leaf renamed")
        .expect("rename leaf");
    assert!(status.contains("Renamed leaf node 'n2'"));
    let name = app
        .engine
        .read()
        .unwrap()
        .state()
        .sequences
        .get("seq_b")
        .and_then(|dna| dna.name().clone());
    assert_eq!(name.as_deref(), Some("Leaf renamed"));
}

#[test]
fn remove_leaf_lineage_node_rejects_non_leaf_nodes() {
    let mut app = GENtleApp::default();
    {
        let mut engine = app.engine.write().unwrap();
        let state = engine.state_mut();
        state.sequences.insert(
            "seq_a".to_string(),
            DNAsequence::from_sequence("ATGC").unwrap(),
        );
        state.sequences.insert(
            "seq_b".to_string(),
            DNAsequence::from_sequence("ATGC").unwrap(),
        );
        insert_test_lineage_node(state, "n1", "seq_a");
        insert_test_lineage_node(state, "n2", "seq_b");
        state.lineage.edges.push(LineageEdge {
            from_node_id: "n1".to_string(),
            to_node_id: "n2".to_string(),
            op_id: "op-1".to_string(),
            run_id: "interactive".to_string(),
        });
    }
    let err = app
        .remove_leaf_lineage_node("n1")
        .expect_err("non-leaf removal must fail");
    assert!(err.contains("not a leaf"));
    assert!(
        app.engine
            .read()
            .unwrap()
            .state()
            .sequences
            .contains_key("seq_a")
    );
}

#[test]
fn request_remove_leaf_lineage_node_sets_pending_target() {
    let mut app = GENtleApp::default();
    {
        let mut engine = app.engine.write().unwrap();
        let state = engine.state_mut();
        state.sequences.insert(
            "seq_leaf".to_string(),
            DNAsequence::from_sequence("ATGC").unwrap(),
        );
        insert_test_lineage_node(state, "n_leaf", "seq_leaf");
    }
    let status = app
        .request_remove_leaf_lineage_node("n_leaf")
        .expect("request remove leaf");
    assert!(status.contains("confirm in dialog"));
    assert_eq!(app.lineage_node_remove_target_id.as_deref(), Some("n_leaf"));
}

#[test]
fn refresh_lineage_cache_includes_dotplot_and_flexibility_analysis_nodes() {
    let mut app = GENtleApp::default();
    {
        let mut engine = app.engine.write().unwrap();
        let state = engine.state_mut();
        state.sequences.insert(
            "seq_a".to_string(),
            DNAsequence::from_sequence("ATGCGATCGATCGATCGATC").unwrap(),
        );
        state.sequences.insert(
            "seq_b".to_string(),
            DNAsequence::from_sequence("ATGCGATCAATCGATCGATC").unwrap(),
        );
        insert_test_lineage_node(state, "n_seq_a", "seq_a");
        insert_test_lineage_node(state, "n_seq_b", "seq_b");
    }

    let (dotplot_op_id, flex_op_id) = {
        let mut engine = app.engine.write().unwrap();
        let dotplot_result = engine
            .apply(Operation::ComputeDotplot {
                seq_id: "seq_a".to_string(),
                reference_seq_id: Some("seq_b".to_string()),
                span_start_0based: Some(0),
                span_end_0based: Some(20),
                reference_span_start_0based: Some(0),
                reference_span_end_0based: Some(20),
                mode: DotplotMode::PairForward,
                word_size: 4,
                step_bp: 2,
                max_mismatches: 1,
                tile_bp: None,
                store_as: Some("p53_dp".to_string()),
            })
            .expect("compute dotplot");
        let flex_result = engine
            .apply(Operation::ComputeFlexibilityTrack {
                seq_id: "seq_a".to_string(),
                span_start_0based: Some(0),
                span_end_0based: Some(20),
                model: FlexibilityModel::AtRichness,
                bin_bp: 5,
                smoothing_bp: Some(10),
                store_as: Some("p53_fx".to_string()),
            })
            .expect("compute flexibility track");
        (dotplot_result.op_id, flex_result.op_id)
    };

    app.refresh_lineage_cache_if_needed();

    let dotplot_row = app
        .lineage_rows
        .iter()
        .find(|row| row.node_id == "analysis:dotplot:p53_dp")
        .expect("dotplot lineage row");
    assert_eq!(dotplot_row.kind, LineageNodeKind::Analysis);
    assert_eq!(
        dotplot_row.analysis_kind,
        Some(LineageAnalysisKind::Dotplot)
    );
    assert_eq!(dotplot_row.analysis_artifact_id.as_deref(), Some("p53_dp"));
    assert_eq!(
        dotplot_row.analysis_reference_seq_id.as_deref(),
        Some("seq_b")
    );
    assert_eq!(dotplot_row.created_by_op, dotplot_op_id);
    assert!(dotplot_row.analysis_point_count.unwrap_or(0) > 0);

    let flex_row = app
        .lineage_rows
        .iter()
        .find(|row| row.node_id == "analysis:flex:p53_fx")
        .expect("flexibility lineage row");
    assert_eq!(flex_row.kind, LineageNodeKind::Analysis);
    assert_eq!(
        flex_row.analysis_kind,
        Some(LineageAnalysisKind::FlexibilityTrack)
    );
    assert_eq!(flex_row.analysis_artifact_id.as_deref(), Some("p53_fx"));
    assert_eq!(flex_row.created_by_op, flex_op_id);
    assert!(flex_row.analysis_bin_count.unwrap_or(0) > 0);

    assert!(
        app.lineage_edges
            .iter()
            .any(|(from, to, op_id)| from == "n_seq_a"
                && to == "analysis:dotplot:p53_dp"
                && op_id == &dotplot_op_id)
    );
    assert!(
        app.lineage_edges
            .iter()
            .any(|(from, to, op_id)| from == "n_seq_b"
                && to == "analysis:dotplot:p53_dp"
                && op_id == &dotplot_op_id)
    );
    assert!(
        app.lineage_edges
            .iter()
            .any(|(from, to, op_id)| from == "n_seq_a"
                && to == "analysis:flex:p53_fx"
                && op_id == &flex_op_id)
    );
}

#[test]
fn refresh_lineage_cache_includes_primer_and_qpcr_design_analysis_nodes() {
    let mut app = GENtleApp::default();
    {
        let mut engine = app.engine.write().unwrap();
        let state = engine.state_mut();
        state.sequences.insert(
                "tpl".to_string(),
                DNAsequence::from_sequence(
                    "GGGGGGGGGGGGGGGGGGGGCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAAAAATTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT",
                )
                .unwrap(),
            );
        insert_test_lineage_node(state, "n_tpl", "tpl");
    }

    let (primer_op_id, qpcr_op_id) = {
        let mut engine = app.engine.write().unwrap();
        engine.state_mut().parameters.primer_design_backend =
            crate::engine::PrimerDesignBackend::Internal;
        let primer_result = engine
            .apply(Operation::DesignPrimerPairs {
                template: "tpl".to_string(),
                roi_start_0based: 30,
                roi_end_0based: 70,
                forward: PrimerDesignSideConstraint {
                    min_length: 20,
                    max_length: 20,
                    location_0based: Some(5),
                    start_0based: None,
                    end_0based: None,
                    min_tm_c: 40.0,
                    max_tm_c: 90.0,
                    min_gc_fraction: 0.0,
                    max_gc_fraction: 1.0,
                    max_anneal_hits: 10,
                    ..Default::default()
                },
                reverse: PrimerDesignSideConstraint {
                    min_length: 20,
                    max_length: 20,
                    location_0based: Some(60),
                    start_0based: None,
                    end_0based: None,
                    min_tm_c: 40.0,
                    max_tm_c: 90.0,
                    min_gc_fraction: 0.0,
                    max_gc_fraction: 1.0,
                    max_anneal_hits: 10,
                    ..Default::default()
                },
                pair_constraints: PrimerDesignPairConstraint::default(),
                min_amplicon_bp: 40,
                max_amplicon_bp: 130,
                max_tm_delta_c: Some(50.0),
                max_pairs: Some(10),
                report_id: Some("tp73_primer".to_string()),
            })
            .expect("design primer pairs");
        let qpcr_result = engine
            .apply(Operation::DesignQpcrAssays {
                template: "tpl".to_string(),
                roi_start_0based: 30,
                roi_end_0based: 70,
                forward: PrimerDesignSideConstraint {
                    min_length: 20,
                    max_length: 20,
                    location_0based: Some(5),
                    start_0based: None,
                    end_0based: None,
                    min_tm_c: 40.0,
                    max_tm_c: 90.0,
                    min_gc_fraction: 0.0,
                    max_gc_fraction: 1.0,
                    max_anneal_hits: 100,
                    ..Default::default()
                },
                reverse: PrimerDesignSideConstraint {
                    min_length: 20,
                    max_length: 20,
                    location_0based: Some(60),
                    start_0based: None,
                    end_0based: None,
                    min_tm_c: 40.0,
                    max_tm_c: 90.0,
                    min_gc_fraction: 0.0,
                    max_gc_fraction: 1.0,
                    max_anneal_hits: 100,
                    ..Default::default()
                },
                probe: PrimerDesignSideConstraint {
                    min_length: 20,
                    max_length: 20,
                    location_0based: Some(35),
                    start_0based: None,
                    end_0based: None,
                    min_tm_c: 40.0,
                    max_tm_c: 90.0,
                    min_gc_fraction: 0.0,
                    max_gc_fraction: 1.0,
                    max_anneal_hits: 100,
                    ..Default::default()
                },
                pair_constraints: PrimerDesignPairConstraint::default(),
                min_amplicon_bp: 40,
                max_amplicon_bp: 130,
                max_tm_delta_c: Some(50.0),
                max_probe_tm_delta_c: Some(50.0),
                max_assays: Some(10),
                transcript_targeting: None,
                report_id: Some("tp73_qpcr".to_string()),
            })
            .expect("design qpcr assays");
        (primer_result.op_id, qpcr_result.op_id)
    };

    app.refresh_lineage_cache_if_needed();

    let primer_row = app
        .lineage_rows
        .iter()
        .find(|row| row.node_id == "analysis:primer:tp73_primer")
        .expect("primer lineage row");
    assert_eq!(primer_row.kind, LineageNodeKind::Analysis);
    assert_eq!(
        primer_row.analysis_kind,
        Some(LineageAnalysisKind::PrimerDesign)
    );
    assert_eq!(
        primer_row.analysis_artifact_id.as_deref(),
        Some("tp73_primer")
    );
    assert_eq!(primer_row.analysis_mode.as_deref(), Some("internal"));
    assert!(primer_row.analysis_target_count.is_some());
    assert_eq!(primer_row.created_by_op, primer_op_id);

    let qpcr_row = app
        .lineage_rows
        .iter()
        .find(|row| row.node_id == "analysis:qpcr:tp73_qpcr")
        .expect("qpcr lineage row");
    assert_eq!(qpcr_row.kind, LineageNodeKind::Analysis);
    assert_eq!(
        qpcr_row.analysis_kind,
        Some(LineageAnalysisKind::QpcrDesign)
    );
    assert_eq!(qpcr_row.analysis_artifact_id.as_deref(), Some("tp73_qpcr"));
    assert_eq!(qpcr_row.analysis_mode.as_deref(), Some("internal"));
    assert!(qpcr_row.analysis_target_count.is_some());
    assert_eq!(qpcr_row.created_by_op, qpcr_op_id);

    assert!(
        app.lineage_edges
            .iter()
            .any(|(from, to, op_id)| from == "n_tpl"
                && to == "analysis:primer:tp73_primer"
                && op_id == &primer_op_id)
    );
    assert!(
        app.lineage_edges
            .iter()
            .any(|(from, to, op_id)| from == "n_tpl"
                && to == "analysis:qpcr:tp73_qpcr"
                && op_id == &qpcr_op_id)
    );
    assert_eq!(
        app.lineage_reopenable_pcr_op_seq_ids.get(&primer_op_id),
        Some(&"tpl".to_string())
    );
    assert_eq!(
        app.lineage_reopenable_pcr_op_seq_ids.get(&qpcr_op_id),
        Some(&"tpl".to_string())
    );
}

#[test]
fn refresh_lineage_cache_includes_restriction_cloning_pcr_handoff_analysis_nodes() {
    let mut app = GENtleApp::default();
    {
        let mut engine = app.engine.write().unwrap();
        let state = engine.state_mut();
        state.sequences.insert(
                "tpl".to_string(),
                DNAsequence::from_sequence(
                    "ACGTTGCATGTCAGTACGATCGTACGTAGCTAGTCGATCGTACGATCGTAGCTAGCATCGATGCTAGCTAGTACGTAGCATCGATCGTAGCTAGCATGCTAGCTAGTCGATCGATCGTACGATCG",
                )
                .unwrap(),
            );
        let mut vector = DNAsequence::from_sequence("AAAAGAATTCGGGGGAAGCTTTTTT").expect("vector");
        *vector.restriction_enzymes_mut() = crate::enzymes::active_restriction_enzymes();
        let vector_len_i64 = vector.len().try_into().expect("vector len fits i64");
        vector.features_mut().push(gb_io::seq::Feature {
            kind: "misc_feature".into(),
            location: gb_io::seq::Location::simple_range(0, vector_len_i64),
            qualifiers: vec![
                ("label".into(), Some("MCS".to_string())),
                (
                    "mcs_expected_sites".into(),
                    Some("EcoRI,HindIII".to_string()),
                ),
            ],
        });
        vector.update_computed_features();
        state.sequences.insert("vec".to_string(), vector);
        insert_test_lineage_node(state, "n_tpl", "tpl");
    }

    let handoff_op_id = {
        let mut engine = app.engine.write().unwrap();
        engine.state_mut().parameters.primer_design_backend =
            crate::engine::PrimerDesignBackend::Internal;
        engine
            .apply(Operation::DesignPrimerPairs {
                template: "tpl".to_string(),
                roi_start_0based: 40,
                roi_end_0based: 80,
                forward: PrimerDesignSideConstraint {
                    min_length: 20,
                    max_length: 20,
                    location_0based: Some(5),
                    start_0based: None,
                    end_0based: None,
                    min_tm_c: 0.0,
                    max_tm_c: 100.0,
                    min_gc_fraction: 0.0,
                    max_gc_fraction: 1.0,
                    max_anneal_hits: 1000,
                    ..Default::default()
                },
                reverse: PrimerDesignSideConstraint {
                    min_length: 20,
                    max_length: 20,
                    location_0based: Some(90),
                    start_0based: None,
                    end_0based: None,
                    min_tm_c: 0.0,
                    max_tm_c: 100.0,
                    min_gc_fraction: 0.0,
                    max_gc_fraction: 1.0,
                    max_anneal_hits: 1000,
                    ..Default::default()
                },
                pair_constraints: PrimerDesignPairConstraint::default(),
                min_amplicon_bp: 40,
                max_amplicon_bp: 150,
                max_tm_delta_c: Some(100.0),
                max_pairs: Some(10),
                report_id: Some("tp73_primer_handoff".to_string()),
            })
            .expect("design primer pairs");
        let handoff = engine
            .apply(Operation::PrepareRestrictionCloningPcrHandoff {
                template: "tpl".to_string(),
                primer_report_id: "tp73_primer_handoff".to_string(),
                pair_index: 0,
                destination_vector_seq_id: "vec".to_string(),
                mode: RestrictionCloningPcrHandoffMode::DirectedPair,
                forward_enzyme: "EcoRI".to_string(),
                reverse_enzyme: Some("HindIII".to_string()),
                forward_leader_5prime: Some("GC".to_string()),
                reverse_leader_5prime: Some("AT".to_string()),
            })
            .expect("prepare restriction-cloning handoff");
        handoff.op_id
    };

    app.refresh_lineage_cache_if_needed();

    let handoff_row = app
        .lineage_rows
        .iter()
        .find(|row| row.origin == "RestrictionCloningPcrHandoff")
        .expect("restriction-cloning lineage row");
    assert_eq!(handoff_row.kind, LineageNodeKind::Analysis);
    assert_eq!(
        handoff_row.analysis_kind,
        Some(LineageAnalysisKind::RestrictionCloningPcrHandoff)
    );
    assert_eq!(
        handoff_row.analysis_reference_seq_id.as_deref(),
        Some("vec")
    );
    assert_eq!(handoff_row.analysis_mode.as_deref(), Some("directed_pair"));
    assert_eq!(handoff_row.analysis_status.as_deref(), Some("compatible"));
    assert_eq!(handoff_row.created_by_op, handoff_op_id);
    assert!(
        app.lineage_edges
            .iter()
            .any(|(from, to, op_id)| from == "n_tpl"
                && to == &handoff_row.node_id
                && op_id == &handoff_op_id)
    );
}

#[test]
fn refresh_lineage_cache_includes_uniprot_projection_analysis_nodes() {
    let temp = tempdir().expect("tempdir");
    let swiss_path = temp.path().join("toy_uniprot.swiss");
    let swiss_text = r#"
ID   PTEST1_HUMAN             Reviewed;         81 AA.
AC   PTEST1;
DE   RecName: Full=Toy protein 1;
GN   Name=TOY1;
OS   Homo sapiens.
FT   DOMAIN          5..25
FT                   /note="Toy domain"
DR   Ensembl; TX1; ENSPTEST1; ENSTTEST1.
SQ   SEQUENCE   81 AA;  900 MW;  ABC CRC64;
     MAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
//
"#;
    fs::write(&swiss_path, swiss_text).expect("write swiss toy");

    let mut state = ProjectState::default();
    let mut dna = DNAsequence::from_sequence(&"ACGT".repeat(300)).expect("valid DNA");
    dna.features_mut().push(gb_io::seq::Feature {
        kind: "mRNA".into(),
        location: gb_io::seq::Location::simple_range(99, 360),
        qualifiers: vec![
            ("gene".into(), Some("TOY1".to_string())),
            ("transcript_id".into(), Some("TX1".to_string())),
            ("label".into(), Some("TX1".to_string())),
            (
                "cds_ranges_1based".into(),
                Some("100-180,300-360".to_string()),
            ),
        ]
        .into_iter()
        .collect(),
    });
    state.sequences.insert("toy_seq".to_string(), dna);
    insert_test_lineage_node(&mut state, "n_toy", "toy_seq");

    let mut app = GENtleApp::default();
    app.engine = Arc::new(RwLock::new(GentleEngine::from_state(state)));
    app.uniprot_swiss_path = swiss_path.display().to_string();
    app.import_uniprot_swiss_prot_from_dialog();
    let projection_op_id = {
        let mut engine = app.engine.write().expect("engine lock");
        engine
            .apply(Operation::ProjectUniprotToGenome {
                seq_id: "toy_seq".to_string(),
                entry_id: "PTEST1".to_string(),
                projection_id: Some("toy_projection".to_string()),
                transcript_id: None,
            })
            .expect("project uniprot")
            .op_id
    };

    app.refresh_lineage_cache_if_needed();

    let row = app
        .lineage_rows
        .iter()
        .find(|row| row.node_id == "analysis:uniprot:toy_projection")
        .expect("uniprot lineage row");
    assert_eq!(row.kind, LineageNodeKind::Analysis);
    assert_eq!(
        row.analysis_kind,
        Some(LineageAnalysisKind::UniprotProjection)
    );
    assert_eq!(row.analysis_artifact_id.as_deref(), Some("toy_projection"));
    assert_eq!(row.analysis_mode.as_deref(), Some("PTEST1"));
    assert_eq!(row.analysis_status.as_deref(), Some("all_transcripts"));
    assert!(row.analysis_target_count.unwrap_or(0) > 0);
    assert_eq!(row.created_by_op, projection_op_id);

    assert!(
        app.lineage_edges
            .iter()
            .any(|(from, to, op_id)| from == "n_toy"
                && to == "analysis:uniprot:toy_projection"
                && op_id == &projection_op_id)
    );
}

#[test]
fn refresh_lineage_cache_includes_rna_read_analysis_nodes() {
    let mut app = GENtleApp::default();
    {
        let mut engine = app.engine.write().unwrap();
        let state = engine.state_mut();
        let mut dna = DNAsequence::from_sequence("ATGGAATTTACGTACGT").expect("sequence");
        dna.features_mut().push(gb_io::seq::Feature {
            kind: "mRNA".into(),
            location: gb_io::seq::Location::simple_range(0, 18),
            qualifiers: vec![("gene".into(), Some("TP73".to_string()))],
        });
        state.sequences.insert("seq_rna".to_string(), dna);
        insert_test_lineage_node(state, "n_rna", "seq_rna");
        let op_id = engine
            .commit_rna_read_report(crate::engine::RnaReadInterpretationReport {
                schema: "gentle.rna_read_report.v1".to_string(),
                report_id: "tp73_rna".to_string(),
                seq_id: "seq_rna".to_string(),
                seed_feature_id: 0,
                generated_at_unix_ms: 321,
                profile: crate::engine::RnaReadInterpretationProfile::NanoporeCdnaV1,
                input_path: "tp73_reads.fa".to_string(),
                input_format: crate::engine::RnaReadInputFormat::Fasta,
                report_mode: crate::engine::RnaReadReportMode::Full,
                scope: crate::engine::SplicingScopePreset::TargetGroupTargetStrand,
                origin_mode: crate::engine::RnaReadOriginMode::SingleGene,
                target_gene_ids: vec!["TP73".to_string()],
                read_count_total: 24,
                read_count_seed_passed: 7,
                read_count_aligned: 3,
                ..crate::engine::RnaReadInterpretationReport::default()
            })
            .expect("commit RNA-read report")
            .op_id;
        assert!(!op_id.is_empty());
    }

    app.refresh_lineage_cache_if_needed();

    let row = app
        .lineage_rows
        .iter()
        .find(|row| row.node_id == "analysis:rna_reads:tp73_rna")
        .expect("RNA-read lineage row");
    assert_eq!(row.kind, LineageNodeKind::Analysis);
    assert_eq!(
        row.analysis_kind,
        Some(LineageAnalysisKind::RnaReadInterpretation)
    );
    assert_eq!(row.analysis_artifact_id.as_deref(), Some("tp73_rna"));
    assert_eq!(row.analysis_mode.as_deref(), Some("nanopore_cdna_v1"));
    assert_eq!(row.analysis_status.as_deref(), Some("full / single_gene"));
    assert_eq!(row.analysis_point_count, Some(7));
    assert_eq!(row.analysis_read_count, Some(24));
    assert_eq!(row.analysis_variant_count, Some(3));
    assert_eq!(row.analysis_target_count, Some(1));
    assert_ne!(row.created_by_op, "-");

    assert!(
        app.lineage_edges
            .iter()
            .any(|(from, to, op_id)| from == "n_rna"
                && to == "analysis:rna_reads:tp73_rna"
                && op_id == &row.created_by_op)
    );
}

#[test]
fn refresh_lineage_cache_includes_protein_derivation_analysis_nodes() {
    let mut app = GENtleApp::default();
    {
        let mut engine = app.engine.write().unwrap();
        let state = engine.state_mut();
        let mut dna = DNAsequence::from_sequence("ATGAAACCCTAA").expect("sequence");
        dna.features_mut().push(gb_io::seq::Feature {
            kind: "source".into(),
            location: gb_io::seq::Location::simple_range(0, 12),
            qualifiers: vec![("organism".into(), Some("Escherichia coli".to_string()))],
        });
        dna.features_mut().push(gb_io::seq::Feature {
            kind: "mRNA".into(),
            location: gb_io::seq::Location::simple_range(0, 12),
            qualifiers: vec![
                ("gene".into(), Some("toyA".to_string())),
                ("transcript_id".into(), Some("TX_TOY".to_string())),
                ("label".into(), Some("TX_TOY".to_string())),
            ],
        });
        dna.features_mut().push(gb_io::seq::Feature {
            kind: "CDS".into(),
            location: gb_io::seq::Location::simple_range(0, 12),
            qualifiers: vec![
                ("transcript_id".into(), Some("TX_TOY".to_string())),
                ("product".into(), Some("Toy enzyme".to_string())),
                ("protein_id".into(), Some("PROT_TOY".to_string())),
                ("transl_table".into(), Some("11".to_string())),
            ],
        });
        state.sequences.insert("protein_source".to_string(), dna);
        insert_test_lineage_node(state, "n_protein_source", "protein_source");
    }

    let protein_op_id = {
        let mut engine = app.engine.write().unwrap();
        engine
            .apply(Operation::DeriveProteinSequences {
                seq_id: "protein_source".to_string(),
                feature_ids: vec![1],
                feature_query: None,
                scope: None,
                output_prefix: Some("tp73_protein".to_string()),
                report_id: None,
            })
            .expect("derive proteins")
            .op_id
    };

    app.refresh_lineage_cache_if_needed();

    let row = app
        .lineage_rows
        .iter()
        .find(|row| row.node_id == format!("analysis:protein_derive:{protein_op_id}"))
        .expect("protein-derivation lineage row");
    assert_eq!(row.kind, LineageNodeKind::Analysis);
    assert_eq!(
        row.analysis_kind,
        Some(LineageAnalysisKind::ProteinDerivation)
    );
    assert_eq!(
        row.analysis_artifact_id.as_deref(),
        Some(protein_op_id.as_str())
    );
    assert_eq!(row.analysis_mode.as_deref(), Some("annotated_cds"));
    assert_eq!(row.analysis_target_count, Some(1));
    assert_eq!(row.created_by_op, protein_op_id);

    assert!(
        app.lineage_edges
            .iter()
            .any(|(from, to, op_id)| from == "n_protein_source"
                && to == &format!("analysis:protein_derive:{protein_op_id}")
                && op_id == &protein_op_id)
    );
}

#[test]
fn refresh_lineage_cache_includes_reverse_translation_analysis_nodes() {
    let mut app = GENtleApp::default();
    {
        let mut engine = app.engine.write().unwrap();
        let state = engine.state_mut();
        let mut protein = DNAsequence::from_sequence("MKP").expect("protein");
        protein.set_name("Toy protein");
        protein.set_molecule_type("protein");
        state.sequences.insert("seq_reverse".to_string(), protein);
        insert_test_lineage_node(state, "n_reverse", "seq_reverse");
    }

    let reverse_op_id = {
        let mut engine = app.engine.write().unwrap();
        engine
            .apply(Operation::ReverseTranslateProteinSequence {
                seq_id: "seq_reverse".to_string(),
                output_id: Some("seq_reverse_coding".to_string()),
                speed_profile: Some(crate::engine::TranslationSpeedProfile::Ecoli),
                speed_mark: Some(crate::engine::TranslationSpeedMark::Slow),
                translation_table: Some(11),
                target_anneal_tm_c: Some(58.0),
                anneal_window_bp: Some(9),
            })
            .expect("reverse translate")
            .op_id
    };

    app.refresh_lineage_cache_if_needed();

    let row = app
        .lineage_rows
        .iter()
        .find(|row| row.node_id == format!("analysis:reverse_translate:{reverse_op_id}"))
        .expect("reverse-translation lineage row");
    assert_eq!(row.kind, LineageNodeKind::Analysis);
    assert_eq!(
        row.analysis_kind,
        Some(LineageAnalysisKind::ReverseTranslation)
    );
    assert_eq!(
        row.analysis_artifact_id.as_deref(),
        Some(reverse_op_id.as_str())
    );
    assert_eq!(
        row.analysis_reference_seq_id.as_deref(),
        Some("seq_reverse_coding")
    );
    assert_eq!(row.analysis_mode.as_deref(), Some("ecoli:slow"));
    assert!(
        row.analysis_status
            .as_deref()
            .is_some_and(|value| value.contains("gc="))
    );
    assert_eq!(row.analysis_point_count, Some(11));
    assert_eq!(row.analysis_target_count, Some(9));
    assert_eq!(row.analysis_variant_count, Some(3));
    assert_eq!(row.created_by_op, reverse_op_id);

    assert!(
        app.lineage_edges
            .iter()
            .any(|(from, to, op_id)| from == "n_reverse"
                && to == &format!("analysis:reverse_translate:{reverse_op_id}")
                && op_id == &reverse_op_id)
    );
}

#[test]
fn refresh_lineage_cache_includes_construct_reasoning_analysis_nodes() {
    let mut app = GENtleApp::default();
    {
        let mut engine = app.engine.write().unwrap();
        let state = engine.state_mut();
        let mut dna = DNAsequence::from_sequence("ATGGAATTTACGTACGT").expect("sequence");
        dna.features_mut().push(gb_io::seq::Feature {
            kind: "variation".into(),
            location: gb_io::seq::Location::simple_range(3, 4),
            qualifiers: vec![
                ("label".into(), Some("rsGui".to_string())),
                (
                    "gentle_generated".into(),
                    Some("genome_vcf_track".to_string()),
                ),
                ("vcf_ref".into(), Some("G".to_string())),
                ("vcf_alt".into(), Some("A".to_string())),
            ],
        });
        state.sequences.insert("seq_reasoning".to_string(), dna);
        insert_test_lineage_node(state, "n_reasoning", "seq_reasoning");
    }

    let graph_id = {
        let mut engine = app.engine.write().unwrap();
        engine
            .build_construct_reasoning_graph("seq_reasoning", None, Some("reasoning_demo"))
            .expect("build reasoning graph")
            .graph_id
    };

    app.refresh_lineage_cache_if_needed();

    let row = app
        .lineage_rows
        .iter()
        .find(|row| row.node_id == format!("analysis:construct_reasoning:{graph_id}"))
        .expect("construct-reasoning lineage row");
    assert_eq!(row.kind, LineageNodeKind::Analysis);
    assert_eq!(
        row.analysis_kind,
        Some(LineageAnalysisKind::ConstructReasoning)
    );
    assert_eq!(row.analysis_artifact_id.as_deref(), Some(graph_id.as_str()));
    assert_eq!(
        row.analysis_mode.as_deref(),
        Some("construct_objective_seq_reasoning")
    );
    assert!(row.analysis_status.is_some());
    assert!(row.analysis_point_count.unwrap_or(0) > 0);
    assert!(row.analysis_variant_count.unwrap_or(0) > 0);
    assert_eq!(row.created_by_op, "-");

    assert!(
        app.lineage_edges
            .iter()
            .any(|(from, to, op_id)| from == "n_reasoning"
                && to == &format!("analysis:construct_reasoning:{graph_id}")
                && op_id == &format!("analysis:construct_reasoning:{graph_id}"))
    );
}

#[test]
fn refresh_lineage_cache_marks_protein_handoff_construct_reasoning_nodes_distinctly() {
    let mut app = GENtleApp::default();
    {
        let mut engine = app.engine.write().unwrap();
        let state = engine.state_mut();
        state.sequences.insert(
            "handoff_target".to_string(),
            DNAsequence::from_sequence("ACGTACGTACGT").expect("dna"),
        );
        let mut protein = DNAsequence::from_sequence("MKP").expect("protein");
        protein.set_name("Toy protein");
        protein.set_molecule_type("protein");
        state
            .sequences
            .insert("handoff_protein".to_string(), protein);
        insert_test_lineage_node(state, "n_handoff", "handoff_target");
    }

    let graph_id = {
        let mut engine = app.engine.write().unwrap();
        engine
            .apply(Operation::BuildProteinToDnaHandoffReasoning {
                seq_id: "handoff_target".to_string(),
                protein_seq_id: "handoff_protein".to_string(),
                transcript_filter: None,
                projection_id: None,
                ensembl_entry_id: Some("ENSPTOY1".to_string()),
                feature_query: None,
                ranking_goal: ProteinToDnaHandoffRankingGoal::BalancedProvenance,
                speed_profile: Some(TranslationSpeedProfile::Ecoli),
                speed_mark: Some(TranslationSpeedMark::Slow),
                translation_table: Some(11),
                target_anneal_tm_c: Some(58.0),
                anneal_window_bp: Some(9),
                objective_id: None,
                graph_id: Some("handoff_lineage_demo".to_string()),
            })
            .expect("build protein handoff graph")
            .construct_reasoning_graph
            .expect("graph")
            .graph_id
    };

    app.refresh_lineage_cache_if_needed();

    let row = app
        .lineage_rows
        .iter()
        .find(|row| row.node_id == format!("analysis:construct_reasoning:{graph_id}"))
        .expect("protein handoff lineage row");
    assert_eq!(row.kind, LineageNodeKind::Analysis);
    assert_eq!(
        row.analysis_kind,
        Some(LineageAnalysisKind::ConstructReasoning)
    );
    assert_eq!(row.origin, "ConstructReasoning / ProteinToDnaHandoff");
    assert_eq!(row.analysis_artifact_id.as_deref(), Some(graph_id.as_str()));
}

#[test]
fn refresh_lineage_cache_includes_sequencing_confirmation_analysis_nodes() {
    let mut app = GENtleApp::default();
    {
        let mut engine = app.engine.write().unwrap();
        let state = engine.state_mut();
        state.sequences.insert(
            "expected".to_string(),
            DNAsequence::from_sequence("ACGGACGT").unwrap(),
        );
        state.sequences.insert(
            "baseline".to_string(),
            DNAsequence::from_sequence("ACGTACGT").unwrap(),
        );
        state.sequences.insert(
            "read_expected".to_string(),
            DNAsequence::from_sequence("ACGGACGT").unwrap(),
        );
        insert_test_lineage_node(state, "n_expected", "expected");
        insert_test_lineage_node(state, "n_baseline", "baseline");
    }

    let confirmation_op_id = {
        let mut engine = app.engine.write().unwrap();
        engine
            .apply(Operation::ConfirmConstructReads {
                expected_seq_id: "expected".to_string(),
                baseline_seq_id: Some("baseline".to_string()),
                read_seq_ids: vec!["read_expected".to_string()],
                trace_ids: vec![],
                targets: vec![],
                alignment_mode: PairwiseAlignmentMode::Local,
                match_score: 2,
                mismatch_score: -3,
                gap_open: -5,
                gap_extend: -1,
                min_identity_fraction: 0.8,
                min_target_coverage_fraction: 0.8,
                allow_reverse_complement: true,
                report_id: Some("tp73_confirm".to_string()),
            })
            .expect("confirm construct")
            .op_id
    };

    app.refresh_lineage_cache_if_needed();

    let row = app
        .lineage_rows
        .iter()
        .find(|row| row.node_id == "analysis:seq_confirm:tp73_confirm")
        .expect("sequencing-confirmation lineage row");
    assert_eq!(row.kind, LineageNodeKind::Analysis);
    assert_eq!(
        row.analysis_kind,
        Some(LineageAnalysisKind::SequencingConfirmation)
    );
    assert_eq!(row.analysis_artifact_id.as_deref(), Some("tp73_confirm"));
    assert_eq!(row.analysis_reference_seq_id.as_deref(), Some("baseline"));
    assert_eq!(row.analysis_status.as_deref(), Some("confirmed"));
    assert_eq!(row.analysis_read_count, Some(1));
    assert_eq!(row.analysis_trace_count, Some(0));
    assert_eq!(row.created_by_op, confirmation_op_id);
    assert!(row.analysis_target_count.unwrap_or(0) > 0);
    assert!(row.analysis_variant_count.unwrap_or(0) > 0);

    assert!(
        app.lineage_edges
            .iter()
            .any(|(from, to, op_id)| from == "n_expected"
                && to == "analysis:seq_confirm:tp73_confirm"
                && op_id == &confirmation_op_id)
    );
    assert!(
        app.lineage_edges
            .iter()
            .any(|(from, to, op_id)| from == "n_baseline"
                && to == "analysis:seq_confirm:tp73_confirm"
                && op_id == &confirmation_op_id)
    );
}

#[test]
fn refresh_lineage_cache_includes_svg_export_analysis_nodes() {
    let mut app = GENtleApp::default();
    {
        let mut engine = app.engine.write().unwrap();
        let state = engine.state_mut();
        state.sequences.insert(
            "seq_a".to_string(),
            DNAsequence::from_sequence("ATGCGATCGATCGATCGATC").unwrap(),
        );
        state.sequences.insert(
            "seq_b".to_string(),
            DNAsequence::from_sequence("ATGCGATCAATCGATCGATC").unwrap(),
        );
        insert_test_lineage_node(state, "n_seq_a", "seq_a");
        insert_test_lineage_node(state, "n_seq_b", "seq_b");
    }
    let temp_sequence_svg = tempfile::NamedTempFile::new().expect("temp sequence svg");
    let sequence_svg_path = temp_sequence_svg.path().with_extension("seq.linear.svg");
    let sequence_svg_path_text = sequence_svg_path.display().to_string();
    let temp_dotplot_svg = tempfile::NamedTempFile::new().expect("temp dotplot svg");
    let dotplot_svg_path = temp_dotplot_svg.path().with_extension("seq.dotplot.svg");
    let dotplot_svg_path_text = dotplot_svg_path.display().to_string();

    let (render_sequence_op_id, render_dotplot_op_id) = {
        let mut engine = app.engine.write().unwrap();
        let dotplot_result = engine
            .apply(Operation::ComputeDotplot {
                seq_id: "seq_a".to_string(),
                reference_seq_id: Some("seq_b".to_string()),
                span_start_0based: Some(0),
                span_end_0based: Some(20),
                reference_span_start_0based: Some(0),
                reference_span_end_0based: Some(20),
                mode: DotplotMode::PairForward,
                word_size: 4,
                step_bp: 2,
                max_mismatches: 1,
                tile_bp: None,
                store_as: Some("svg_dp".to_string()),
            })
            .expect("compute dotplot");
        let render_sequence_result = engine
            .apply(Operation::RenderSequenceSvg {
                seq_id: "seq_a".to_string(),
                mode: RenderSvgMode::Linear,
                path: sequence_svg_path_text.clone(),
            })
            .expect("render sequence svg");
        let render_dotplot_result = engine
            .apply(Operation::RenderDotplotSvg {
                seq_id: "seq_a".to_string(),
                dotplot_id: "svg_dp".to_string(),
                path: dotplot_svg_path_text.clone(),
                flex_track_id: None,
                display_density_threshold: None,
                display_intensity_gain: None,
                overlay_x_axis_mode: Default::default(),
                overlay_anchor_exon: None,
            })
            .expect("render dotplot svg");
        assert!(dotplot_result.messages.iter().any(|m| m.contains("svg_dp")));
        (render_sequence_result.op_id, render_dotplot_result.op_id)
    };

    app.refresh_lineage_cache_if_needed();

    let sequence_export_node = format!("analysis:export:{render_sequence_op_id}");
    let sequence_export_row = app
        .lineage_rows
        .iter()
        .find(|row| row.node_id == sequence_export_node)
        .expect("sequence SVG export row");
    assert_eq!(sequence_export_row.kind, LineageNodeKind::Analysis);
    assert_eq!(sequence_export_row.origin, "SequenceSvgExport");
    assert_eq!(sequence_export_row.seq_id, "seq_a");
    assert_eq!(sequence_export_row.analysis_mode.as_deref(), Some("linear"));
    assert!(
        app.lineage_edges
            .iter()
            .any(|(from, to, op_id)| from == "n_seq_a"
                && to == &sequence_export_node
                && op_id == &render_sequence_op_id)
    );

    let dotplot_export_node = format!("analysis:export:{render_dotplot_op_id}");
    let dotplot_export_row = app
        .lineage_rows
        .iter()
        .find(|row| row.node_id == dotplot_export_node)
        .expect("dotplot SVG export row");
    assert_eq!(dotplot_export_row.kind, LineageNodeKind::Analysis);
    assert_eq!(dotplot_export_row.origin, "DotplotSvgExport");
    assert_eq!(
        dotplot_export_row.analysis_kind,
        Some(LineageAnalysisKind::Dotplot)
    );
    assert_eq!(
        dotplot_export_row.analysis_artifact_id.as_deref(),
        Some("svg_dp")
    );
    assert!(
        app.lineage_edges
            .iter()
            .any(|(from, to, op_id)| from == "n_seq_a"
                && to == &dotplot_export_node
                && op_id == &render_dotplot_op_id)
    );
    assert!(
        app.lineage_edges
            .iter()
            .any(|(from, to, op_id)| from == "n_seq_b"
                && to == &dotplot_export_node
                && op_id == &render_dotplot_op_id)
    );
}

#[test]
fn lineage_operation_symbol_prefers_common_operation_families() {
    assert_eq!(GENtleApp::lineage_operation_symbol("Digest: input=a"), "D");
    assert_eq!(
        GENtleApp::lineage_operation_symbol("Ligation: inputs=a,b"),
        "L"
    );
    assert_eq!(
        GENtleApp::lineage_operation_symbol("PCR advanced: template=x"),
        "P"
    );
    assert_eq!(
        GENtleApp::lineage_operation_symbol("Reverse complement: input=x"),
        "RC"
    );
    assert_eq!(
        GENtleApp::lineage_operation_symbol("Molecular weight filter: inputs=x"),
        "F"
    );
    assert_eq!(
        GENtleApp::lineage_operation_symbol("Compute dotplot: seq_id=a"),
        "DP"
    );
    assert_eq!(
        GENtleApp::lineage_operation_symbol("Compute flexibility track: seq_id=a"),
        "FX"
    );
    assert_eq!(
        GENtleApp::lineage_operation_symbol("Sequencing confirmation: expected=a"),
        "SC"
    );
    assert_eq!(GENtleApp::lineage_operation_symbol("unknown op"), "U");
}

#[test]
fn summarize_operation_compute_dotplot_includes_mode_and_spans() {
    let op = Operation::ComputeDotplot {
        seq_id: "seq_a".to_string(),
        reference_seq_id: Some("seq_b".to_string()),
        span_start_0based: Some(10),
        span_end_0based: Some(120),
        reference_span_start_0based: Some(20),
        reference_span_end_0based: Some(180),
        mode: DotplotMode::PairForward,
        word_size: 11,
        step_bp: 40,
        max_mismatches: 1,
        tile_bp: Some(500),
        store_as: Some("family_dp".to_string()),
    };
    let summary = GENtleApp::summarize_operation(&op);
    assert!(summary.contains("Compute dotplot"));
    assert!(summary.contains("mode=pair_forward"));
    assert!(summary.contains("query_span=10..120"));
    assert!(summary.contains("reference_span=20..180"));
    assert!(summary.contains("store_as=family_dp"));
}

#[test]
fn summarize_operation_compute_flexibility_track_includes_model_and_binning() {
    let op = Operation::ComputeFlexibilityTrack {
        seq_id: "seq_a".to_string(),
        span_start_0based: Some(0),
        span_end_0based: Some(250),
        model: FlexibilityModel::AtSkew,
        bin_bp: 25,
        smoothing_bp: Some(75),
        store_as: Some("family_fx".to_string()),
    };
    let summary = GENtleApp::summarize_operation(&op);
    assert!(summary.contains("Compute flexibility track"));
    assert!(summary.contains("model=at_skew"));
    assert!(summary.contains("span=0..250"));
    assert!(summary.contains("bin_bp=25"));
    assert!(summary.contains("smoothing_bp=75"));
    assert!(summary.contains("store_as=family_fx"));
}

#[test]
fn summarize_operation_blast_import_includes_invocation_preview() {
    let op = Operation::ImportBlastHitsTrack {
        seq_id: "query".to_string(),
        hits: vec![BlastHitFeatureInput {
            subject_id: "chr1".to_string(),
            query_start_1based: 1,
            query_end_1based: 8,
            subject_start_1based: 100,
            subject_end_1based: 107,
            identity_percent: 99.0,
            bit_score: 42.0,
            evalue: 1e-8,
            query_coverage_percent: Some(100.0),
        }],
        track_name: Some("blast_hits".to_string()),
        clear_existing: Some(true),
        blast_provenance: Some(BlastInvocationProvenance {
            genome_id: "grch38".to_string(),
            query_label: "query".to_string(),
            query_length: 8,
            max_hits: 20,
            task: "blastn-short".to_string(),
            blastn_executable: "blastn".to_string(),
            blast_db_prefix: "/tmp/db".to_string(),
            command: vec![
                "-db".to_string(),
                "/tmp/db".to_string(),
                "-query".to_string(),
                "/tmp/query.fa".to_string(),
            ],
            command_line: "blastn -db /tmp/db -query /tmp/query.fa".to_string(),
            catalog_path: None,
            cache_dir: None,
            options_override_json: None,
            effective_options_json: None,
        }),
    };
    let summary = GENtleApp::summarize_operation(&op);
    assert!(summary.contains("invocation="));
    assert!(summary.contains("blastn -db /tmp/db"));
}

#[test]
fn lineage_operation_glyph_collapses_mixed_edge_families_to_count_symbol() {
    let labels = vec![
        "Digest: input=a".to_string(),
        "Ligation: inputs=a,b".to_string(),
    ];
    let (symbol, color) = GENtleApp::lineage_operation_glyph(&labels);
    assert_eq!(symbol, "2");
    assert_eq!(color, egui::Color32::from_rgb(102, 102, 102));
}

#[test]
fn lineage_edge_groups_collapses_parallel_edges_by_endpoints() {
    let edges = vec![
        ("n1".to_string(), "n2".to_string(), "opA".to_string()),
        ("n1".to_string(), "n2".to_string(), "opB".to_string()),
        ("n2".to_string(), "n3".to_string(), "opC".to_string()),
        ("n1".to_string(), "n2".to_string(), "opA".to_string()),
    ];
    let groups = GENtleApp::lineage_edge_groups(&edges);
    assert_eq!(groups.len(), 2);
    assert_eq!(groups[0].0, "n1");
    assert_eq!(groups[0].1, "n2");
    assert_eq!(groups[0].2, vec!["opA".to_string(), "opB".to_string()]);
    assert_eq!(groups[1].0, "n2");
    assert_eq!(groups[1].1, "n3");
    assert_eq!(groups[1].2, vec!["opC".to_string()]);
}

#[test]
fn lineage_collapsed_group_hidden_op_badges_summarizes_member_hidden_ops() {
    let edges = vec![
        ("n2".to_string(), "n3".to_string(), "op_digest".to_string()),
        ("n2".to_string(), "n4".to_string(), "op_digest".to_string()),
        (
            "n3".to_string(),
            "n2".to_string(),
            "op_ligation".to_string(),
        ),
        ("n2".to_string(), "n1".to_string(), "op_reverse".to_string()),
        ("n1".to_string(), "n3".to_string(), "op_visible".to_string()),
    ];
    let groups = vec![
        PersistedLineageNodeGroup {
            group_id: "grp-1".to_string(),
            label: "Collapsed".to_string(),
            representative_node_id: "n1".to_string(),
            member_node_ids: vec!["n2".to_string()],
            collapsed: true,
        },
        PersistedLineageNodeGroup {
            group_id: "grp-2".to_string(),
            label: "Expanded".to_string(),
            representative_node_id: "n3".to_string(),
            member_node_ids: vec!["n4".to_string()],
            collapsed: false,
        },
    ];
    let op_label_by_id: HashMap<String, String> = [
        ("op_digest".to_string(), "Digest: input=a".to_string()),
        (
            "op_ligation".to_string(),
            "Ligation: inputs=a,b".to_string(),
        ),
        (
            "op_reverse".to_string(),
            "Reverse complement: input=x".to_string(),
        ),
        ("op_visible".to_string(), "PCR: template=z".to_string()),
    ]
    .into_iter()
    .collect();

    let badges =
        GENtleApp::lineage_collapsed_group_hidden_op_badges(&edges, &groups, &op_label_by_id);
    assert_eq!(badges.len(), 1);
    let badge = badges.get("n1").expect("collapsed group badge");
    assert_eq!(badge.total_ops, 3);
    let family_counts: HashMap<String, usize> = badge.families.iter().cloned().collect();
    assert_eq!(family_counts.get("digest"), Some(&1));
    assert_eq!(family_counts.get("ligation"), Some(&1));
    assert_eq!(family_counts.get("reverse_complement"), Some(&1));
    assert_eq!(family_counts.get("pcr"), None);
    let summary = GENtleApp::lineage_hidden_op_families_summary(badge, 2);
    assert!(summary.contains("D"));
    assert!(summary.contains("L"));
    assert!(summary.contains("+1"));
}

#[test]
fn start_lineage_group_draft_from_marked_populates_form() {
    let mut app = GENtleApp::default();
    app.lineage_group_marked_nodes.insert("n2".to_string());
    app.lineage_group_marked_nodes.insert("n3".to_string());
    let valid_node_ids: HashSet<String> = ["n1", "n2", "n3"]
        .iter()
        .map(|node_id| (*node_id).to_string())
        .collect();

    let count = app
        .start_lineage_group_draft_from_marked("n1", &valid_node_ids)
        .expect("draft from marked");

    assert_eq!(count, 2);
    assert_eq!(app.lineage_group_form_representative, "n1");
    assert_eq!(app.lineage_group_form_members, "n2, n3");
}

#[test]
fn apply_lineage_group_form_clears_marked_nodes_on_success() {
    let mut app = GENtleApp::default();
    app.lineage_group_form_label = "My Group".to_string();
    app.lineage_group_form_representative = "n1".to_string();
    app.lineage_group_form_members = "n2, n3".to_string();
    app.lineage_group_marked_nodes.insert("n2".to_string());
    app.lineage_group_marked_nodes.insert("n3".to_string());
    let valid_node_ids: HashSet<String> = ["n1", "n2", "n3"]
        .iter()
        .map(|node_id| (*node_id).to_string())
        .collect();

    let status = app
        .apply_lineage_group_form(&valid_node_ids)
        .expect("apply group form");

    assert!(status.contains("Saved group"));
    assert_eq!(app.lineage_node_groups.len(), 1);
    assert!(app.lineage_group_marked_nodes.is_empty());
}

#[test]
fn prepare_dialog_ensembl_update_action_depends_on_catalog_template_metadata() {
    let td = tempdir().unwrap();
    let root = td.path();
    let with_template = root.join("with_template.json");
    let without_template = root.join("without_template.json");
    fs::write(
            &with_template,
            r#"{
  "Mouse GRCm39 Ensembl 116": {
    "sequence_remote": "https://ftp.ensembl.org/pub/release-116/vertebrates/fasta/mus_musculus/dna/Mus_musculus.GRCm39.dna_sm.toplevel.fa.gz",
    "annotations_remote": "https://ftp.ensembl.org/pub/release-116/vertebrates/gtf/mus_musculus/Mus_musculus.GRCm39.116.gtf.gz",
    "ensembl_template": {
      "provider": "ensembl",
      "collection": "vertebrates",
      "species_dir": "mus_musculus",
      "file_stem": "Mus_musculus.GRCm39",
      "release": 116
    }
  }
}"#,
        )
        .unwrap();
    fs::write(
        &without_template,
        r#"{
  "LocalProject": {
    "sequence_local": "../test_files/fixtures/genomes/AB011549.2.fa",
    "annotations_local": "../test_files/fixtures/genomes/AB011549.2.gb"
  }
}"#,
    )
    .unwrap();

    let mut app = GENtleApp::default();
    app.genome_catalog_path = with_template.display().to_string();
    assert!(app.selected_genome_catalog_has_ensembl_templates());

    app.genome_catalog_path = without_template.display().to_string();
    assert!(!app.selected_genome_catalog_has_ensembl_templates());
}

#[test]
fn clear_prepare_dialog_ephemeral_state_drops_pending_ensembl_auxiliary_state() {
    let mut app = GENtleApp::default();
    app.pending_ensembl_catalog_update = Some(PendingEnsemblCatalogUpdateDialog {
        scope: GenomeDialogScope::Reference,
        catalog_path: "assets/genomes.json".to_string(),
        preview: EnsemblCatalogUpdatePreview::default(),
        output_catalog_path: String::new(),
    });
    app.pending_ensembl_installable_genomes = Some(PendingEnsemblInstallableGenomeDialog {
        scope: GenomeDialogScope::Reference,
        collection_filter: "all".to_string(),
        filter: "human".to_string(),
        report: EnsemblInstallableGenomeCatalog::default(),
    });
    app.pending_ensembl_quick_install = Some(PendingEnsemblQuickInstallDialog {
        scope: GenomeDialogScope::Reference,
        collection: "vertebrates".to_string(),
        species_dir: "mus_musculus".to_string(),
        genome_id: "Mouse GRCm39 Ensembl 116".to_string(),
        output_catalog_path: "data/genomes_overlay.json".to_string(),
        preview: EnsemblQuickInstallPreview {
            collection: "vertebrates".to_string(),
            species_dir: "mus_musculus".to_string(),
            display_name: "Mouse".to_string(),
            file_stem: "Mus_musculus.GRCm39".to_string(),
            release: 116,
            genome_id: "Mouse GRCm39 Ensembl 116".to_string(),
            catalog_origin_label: "assets/genomes.json".to_string(),
            output_catalog_path: "data/genomes_overlay.json".to_string(),
            catalog_write_mode: "overlay_entry".to_string(),
            catalog_entry_action: "add_overlay_entry".to_string(),
            sequence_remote: "https://example.invalid/mouse.fa.gz".to_string(),
            annotations_remote: "https://example.invalid/mouse.gtf.gz".to_string(),
            warnings: vec![],
        },
    });
    app.genome_prepare_steps = vec![PrepareGenomeUiStepState {
        step_id: PrepareGenomeStepId::Sequence,
        label: "Sequence".to_string(),
        operation_summary: "Download sequence".to_string(),
        determinate_hint: true,
        status: PrepareGenomeUiStepStatus::Running,
        progress_fraction: Some(0.5),
        detail: "demo".to_string(),
        raw_phase: Some("download_sequence".to_string()),
        bytes_done: Some(5),
        bytes_total: Some(10),
        eta_remaining: Some(Duration::from_secs(5)),
    }];

    app.clear_prepare_dialog_ephemeral_state();

    assert!(app.pending_ensembl_catalog_update.is_none());
    assert!(app.pending_ensembl_installable_genomes.is_none());
    assert!(app.pending_ensembl_quick_install.is_none());
    assert!(app.genome_prepare_steps.is_empty());
}

#[test]
fn prepared_reference_catalog_entry_removal_requires_writable_catalog() {
    let td = tempdir().unwrap();
    let root = td.path();
    let catalog_path = root.join("catalog.json");
    fs::write(
        &catalog_path,
        r#"{
  "LocalProject": {
    "sequence_local": "../test_files/fixtures/genomes/AB011549.2.fa",
    "annotations_local": "../test_files/fixtures/genomes/AB011549.2.gb"
  }
}"#,
    )
    .unwrap();

    let mut app = GENtleApp::default();
    app.genome_catalog_path = catalog_path.display().to_string();
    assert!(app.selected_genome_catalog_is_writable());

    let mut perms = fs::metadata(&catalog_path).unwrap().permissions();
    perms.set_readonly(true);
    fs::set_permissions(&catalog_path, perms).unwrap();
    assert!(!app.selected_genome_catalog_is_writable());
}
