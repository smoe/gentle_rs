//! Gene-set source authoring, resolution inspection, and collection operations.
//!
//! This GUI is deliberately thin: it constructs the shared `ResolveGeneSet`
//! operation, inspects persisted portable reports, requires explicit
//! primer-report bindings, and executes collection specificity through the
//! shared shell/engine path on detached engine snapshots.

use super::*;
use crate::{
    engine::{
        CollectionMemberOutcome, CollectionOperationReport, CollectionSubjectRef,
        GeneSetProvenanceRow, GeneSetRequest, GeneSetResolutionReport,
        GeneSetResolutionReviewStatus, PrimerDesignReportSummary,
        PrimerSpecificityCollectionMemberBinding, PrimerSpecificityPolicy,
    },
    engine_shell::{ShellCommand, execute_shell_command},
};

fn gene_set_review_status_label(status: GeneSetResolutionReviewStatus) -> &'static str {
    match status {
        GeneSetResolutionReviewStatus::Unreviewed => "unreviewed",
        GeneSetResolutionReviewStatus::Reviewed => "reviewed",
        GeneSetResolutionReviewStatus::Included => "included",
        GeneSetResolutionReviewStatus::Draft => "draft",
        GeneSetResolutionReviewStatus::Deprecated => "deprecated",
    }
}

#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
enum GeneSetInspectorMode {
    Resolve,
    #[default]
    Inspect,
}

#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
enum GeneSetResolveSourceKind {
    #[default]
    CatalogGroup,
    ExplicitMembers,
    ExternalMapping,
    GenomicNeighbors,
    Random,
}

impl GeneSetResolveSourceKind {
    fn label(self) -> &'static str {
        match self {
            Self::CatalogGroup => "Catalog group",
            Self::ExplicitMembers => "Explicit members",
            Self::ExternalMapping => "External mapping",
            Self::GenomicNeighbors => "Genomic neighbors",
            Self::Random => "Random sampling",
        }
    }
}

#[derive(Clone, Debug)]
struct GeneSetResolveFormState {
    source_kind: GeneSetResolveSourceKind,
    group_id: String,
    members: String,
    external_namespace: String,
    external_id: String,
    anchor_gene: String,
    flank_genes: String,
    flank_bp: String,
    exclude_anchor: bool,
    random_size: String,
    random_seed: String,
    exclude_members: String,
    genome_id: String,
    gene_group_catalog_path: String,
    genome_catalog_path: String,
    cache_dir: String,
    allow_draft: bool,
    allow_deprecated: bool,
    output_path: String,
}

impl Default for GeneSetResolveFormState {
    fn default() -> Self {
        Self {
            source_kind: GeneSetResolveSourceKind::CatalogGroup,
            group_id: String::new(),
            members: String::new(),
            external_namespace: "GO".to_string(),
            external_id: String::new(),
            anchor_gene: String::new(),
            flank_genes: "5".to_string(),
            flank_bp: String::new(),
            exclude_anchor: false,
            random_size: "10".to_string(),
            random_seed: "0".to_string(),
            exclude_members: String::new(),
            genome_id: String::new(),
            gene_group_catalog_path: String::new(),
            genome_catalog_path: String::new(),
            cache_dir: String::new(),
            allow_draft: false,
            allow_deprecated: false,
            output_path: String::new(),
        }
    }
}

fn optional_trimmed(value: &str) -> Option<String> {
    let value = value.trim();
    (!value.is_empty()).then(|| value.to_string())
}

fn parse_gene_set_member_list(
    value: &str,
    label: &str,
    required: bool,
) -> Result<Vec<String>, String> {
    let value = value.trim();
    if value.is_empty() {
        return if required {
            Err(format!("{label} must not be empty"))
        } else {
            Ok(vec![])
        };
    }
    let mut members = Vec::new();
    for token in value.split(',') {
        let token = token.trim();
        if token.is_empty() {
            return Err(format!(
                "{label} contains an empty entry; remove double or trailing commas"
            ));
        }
        members.push(token.to_string());
    }
    Ok(members)
}

fn parse_optional_positive_usize(value: &str, label: &str) -> Result<Option<usize>, String> {
    let value = value.trim();
    if value.is_empty() {
        return Ok(None);
    }
    let parsed = value
        .parse::<usize>()
        .map_err(|_| format!("{label} must be a positive integer"))?;
    if parsed == 0 {
        return Err(format!("{label} must be at least 1"));
    }
    Ok(Some(parsed))
}

fn gene_set_resolve_form_to_operation(form: &GeneSetResolveFormState) -> Result<Operation, String> {
    let source = match form.source_kind {
        GeneSetResolveSourceKind::CatalogGroup => GeneSetRequest::CatalogGroup {
            query: optional_trimmed(&form.group_id)
                .ok_or_else(|| "Group ID must not be empty".to_string())?,
        },
        GeneSetResolveSourceKind::ExplicitMembers => GeneSetRequest::ExplicitMembers {
            members: parse_gene_set_member_list(&form.members, "Members", true)?,
        },
        GeneSetResolveSourceKind::ExternalMapping => {
            let namespace = optional_trimmed(&form.external_namespace)
                .ok_or_else(|| "External namespace must not be empty".to_string())?
                .to_ascii_uppercase();
            let raw_id = optional_trimmed(&form.external_id)
                .ok_or_else(|| "External ID must not be empty".to_string())?;
            let id = if let Some((prefix, suffix)) = raw_id.split_once(':') {
                if !prefix.trim().eq_ignore_ascii_case(&namespace) {
                    return Err(format!(
                        "External ID namespace '{}' does not match selected namespace '{}'",
                        prefix.trim(),
                        namespace
                    ));
                }
                let suffix = suffix.trim();
                if suffix.is_empty() {
                    return Err("External ID suffix must not be empty".to_string());
                }
                raw_id.clone()
            } else {
                format!("{namespace}:{raw_id}")
            };
            GeneSetRequest::ExternalMapping { namespace, id }
        }
        GeneSetResolveSourceKind::GenomicNeighbors => {
            let anchor = optional_trimmed(&form.anchor_gene)
                .ok_or_else(|| "Anchor gene must not be empty".to_string())?;
            let flank_gene_count = parse_optional_positive_usize(&form.flank_genes, "Flank genes")?;
            let flank_bp = parse_optional_positive_usize(&form.flank_bp, "Flank bp")?;
            if flank_gene_count.is_none() && flank_bp.is_none() {
                return Err("Enter flank genes, flank bp, or both".to_string());
            }
            GeneSetRequest::GenomicNeighbors {
                anchor,
                flank_gene_count,
                flank_bp,
                exclude_anchor: form.exclude_anchor,
            }
        }
        GeneSetResolveSourceKind::Random => {
            let count = parse_optional_positive_usize(&form.random_size, "Sample size")?
                .ok_or_else(|| "Sample size must not be empty".to_string())?;
            let random_seed = match optional_trimmed(&form.random_seed) {
                Some(value) => value
                    .parse::<u64>()
                    .map_err(|_| "Random seed must be a non-negative integer".to_string())?,
                None => 0,
            };
            GeneSetRequest::Random {
                count,
                random_seed,
                exclude_members: parse_gene_set_member_list(
                    &form.exclude_members,
                    "Excluded members",
                    false,
                )?,
            }
        }
    };
    Ok(Operation::ResolveGeneSet {
        source,
        genome_id: optional_trimmed(&form.genome_id),
        gene_group_catalog_path: optional_trimmed(&form.gene_group_catalog_path),
        genome_catalog_path: optional_trimmed(&form.genome_catalog_path),
        cache_dir: optional_trimmed(&form.cache_dir),
        allow_draft: form.allow_draft,
        allow_deprecated: form.allow_deprecated,
        path: optional_trimmed(&form.output_path),
    })
}

fn gene_set_provenance_summary(rows: &[GeneSetProvenanceRow]) -> String {
    if rows.is_empty() {
        return "-".to_string();
    }
    rows.iter()
        .map(|row| {
            let mut summary = format!("{}:{}", row.source_kind, row.source_id);
            if let Some(label) = row
                .source_label
                .as_deref()
                .filter(|value| !value.is_empty())
            {
                summary.push_str(&format!(" ({label})"));
            }
            summary
        })
        .collect::<Vec<_>>()
        .join("; ")
}

#[derive(Clone, Debug)]
struct GeneSetResolutionChoice {
    report_id: String,
    report: GeneSetResolutionReport,
}

#[derive(Clone, Debug, Default)]
struct PrimerSpecificityChildPresentation {
    report_id: String,
    primary_seq_id: Option<String>,
    status: String,
    specificity_pass: bool,
    amplicon_count: usize,
    failing_unintended_amplicon_count: usize,
    summary: String,
}

#[derive(Debug)]
struct GeneSetSpecificityTaskResult {
    report: CollectionOperationReport,
    child_reports: BTreeMap<String, PrimerSpecificityChildPresentation>,
}

struct GeneSetSpecificityTask {
    job_id: u64,
    resolution_id: String,
    started: Instant,
    cancel_requested: Arc<AtomicBool>,
    runtime_frame: RuntimeStatusGuard,
    receiver: mpsc::Receiver<Result<GeneSetSpecificityTaskResult, EngineError>>,
}

struct GeneSetResolveTask {
    started: Instant,
    runtime_frame: RuntimeStatusGuard,
    receiver: mpsc::Receiver<Result<GeneSetResolutionReport, EngineError>>,
}

pub(super) struct GeneSetInspectorUiState {
    pub(super) show_panel: bool,
    mode: GeneSetInspectorMode,
    catalog_engine_identity: usize,
    catalog_execution_revision: Option<u64>,
    resolutions: Vec<GeneSetResolutionChoice>,
    primer_reports: Vec<PrimerDesignReportSummary>,
    selected_resolution_id: String,
    member_bindings: BTreeMap<String, String>,
    pair_rank_1based: String,
    target_genome_id: String,
    catalog_path: String,
    cache_dir: String,
    status: String,
    resolve_form: GeneSetResolveFormState,
    resolve_status: String,
    resolve_task: Option<GeneSetResolveTask>,
    task: Option<GeneSetSpecificityTask>,
    last_report: Option<CollectionOperationReport>,
    child_reports: BTreeMap<String, PrimerSpecificityChildPresentation>,
}

impl Default for GeneSetInspectorUiState {
    fn default() -> Self {
        Self {
            show_panel: false,
            mode: GeneSetInspectorMode::Inspect,
            catalog_engine_identity: 0,
            catalog_execution_revision: None,
            resolutions: vec![],
            primer_reports: vec![],
            selected_resolution_id: String::new(),
            member_bindings: BTreeMap::new(),
            pair_rank_1based: "1".to_string(),
            target_genome_id: String::new(),
            catalog_path: String::new(),
            cache_dir: String::new(),
            status: String::new(),
            resolve_form: GeneSetResolveFormState::default(),
            resolve_status: String::new(),
            resolve_task: None,
            task: None,
            last_report: None,
            child_reports: BTreeMap::new(),
        }
    }
}

impl GENtleApp {
    pub(super) fn open_gene_set_inspector_dialog(&mut self) {
        let was_open = self.gene_set_inspector.show_panel;
        self.gene_set_inspector.show_panel = true;
        if self
            .gene_set_inspector
            .resolve_form
            .genome_id
            .trim()
            .is_empty()
        {
            self.gene_set_inspector.resolve_form.genome_id = self.genome_id.clone();
        }
        if self
            .gene_set_inspector
            .resolve_form
            .genome_catalog_path
            .trim()
            .is_empty()
        {
            self.gene_set_inspector.resolve_form.genome_catalog_path =
                self.genome_catalog_path.clone();
        }
        if self
            .gene_set_inspector
            .resolve_form
            .cache_dir
            .trim()
            .is_empty()
        {
            self.gene_set_inspector.resolve_form.cache_dir = self.genome_cache_dir.clone();
        }
        self.refresh_gene_set_inspector_catalog(true);
        if self.gene_set_inspector.resolutions.is_empty() {
            self.gene_set_inspector.mode = GeneSetInspectorMode::Resolve;
        }
        self.mark_window_open_or_focus(Self::gene_set_inspector_viewport_id(), was_open);
    }

    fn refresh_gene_set_inspector_catalog(&mut self, force: bool) {
        let engine_identity = Arc::as_ptr(&self.engine) as usize;
        let Ok(engine) = self.engine.read() else {
            self.gene_set_inspector.status =
                "Could not read project state: engine lock poisoned".to_string();
            return;
        };
        let execution_revision = engine.execution_revision();
        if !force
            && self.gene_set_inspector.catalog_engine_identity == engine_identity
            && self.gene_set_inspector.catalog_execution_revision == Some(execution_revision)
        {
            return;
        }

        let mut resolutions =
            GentleEngine::gene_set_resolution_artifacts_from_state(engine.state())
                .into_iter()
                .map(|report| GeneSetResolutionChoice {
                    report_id: GentleEngine::gene_set_resolution_artifact_id(&report),
                    report,
                })
                .collect::<Vec<_>>();
        resolutions.sort_by(|left, right| left.report_id.cmp(&right.report_id));
        let primer_reports = engine.list_primer_design_reports();
        drop(engine);

        self.gene_set_inspector.catalog_engine_identity = engine_identity;
        self.gene_set_inspector.catalog_execution_revision = Some(execution_revision);
        self.gene_set_inspector.resolutions = resolutions;
        self.gene_set_inspector.primer_reports = primer_reports;
        if !self
            .gene_set_inspector
            .resolutions
            .iter()
            .any(|choice| choice.report_id == self.gene_set_inspector.selected_resolution_id)
        {
            self.gene_set_inspector.selected_resolution_id = self
                .gene_set_inspector
                .resolutions
                .first()
                .map(|choice| choice.report_id.clone())
                .unwrap_or_default();
            self.gene_set_inspector.last_report = None;
            self.gene_set_inspector.child_reports.clear();
        }
        self.reconcile_gene_set_inspector_bindings();
    }

    fn selected_gene_set_resolution(&self) -> Option<&GeneSetResolutionChoice> {
        self.gene_set_inspector
            .resolutions
            .iter()
            .find(|choice| choice.report_id == self.gene_set_inspector.selected_resolution_id)
    }

    fn start_gene_set_resolution_task(&mut self) {
        if self.gene_set_inspector.resolve_task.is_some() {
            self.gene_set_inspector.resolve_status =
                "Gene-set resolution is already running".to_string();
            return;
        }
        let operation =
            match gene_set_resolve_form_to_operation(&self.gene_set_inspector.resolve_form) {
                Ok(operation) => operation,
                Err(error) => {
                    self.gene_set_inspector.resolve_status = error;
                    return;
                }
            };
        let source_label = match &operation {
            Operation::ResolveGeneSet { source, .. } => source.source_kind_label(),
            _ => unreachable!("gene-set form constructs ResolveGeneSet"),
        };
        let engine = self.engine.clone();
        let (tx, receiver) = mpsc::channel::<Result<GeneSetResolutionReport, EngineError>>();
        let runtime_frame = runtime_status_registry().push_with_detail(
            RuntimeStatusFrameKind::BackgroundJob,
            "Resolve gene set",
            Some(format!("source={source_label}")),
        );
        runtime_frame.update_phase("running");
        std::thread::spawn(move || {
            let result =
                crate::background_engine::execute_on_engine_snapshot(&engine, move |snapshot| {
                    snapshot
                        .apply(operation)?
                        .gene_set_resolution
                        .ok_or_else(|| {
                            EngineError::new(
                                ErrorCode::Internal,
                                "ResolveGeneSet completed without a gene-set resolution report",
                            )
                        })
                });
            let _ = tx.send(result);
        });
        self.gene_set_inspector.resolve_status =
            format!("Resolving {source_label} source in a detached project snapshot...");
        self.gene_set_inspector.resolve_task = Some(GeneSetResolveTask {
            started: Instant::now(),
            runtime_frame,
            receiver,
        });
    }

    pub(super) fn poll_gene_set_resolution_task(&mut self, ctx: &egui::Context) {
        let Some(task) = self.gene_set_inspector.resolve_task.as_ref() else {
            return;
        };
        ctx.request_repaint_after(std::time::Duration::from_millis(100));
        let outcome = match task.receiver.try_recv() {
            Ok(result) => Some(result),
            Err(mpsc::TryRecvError::Empty) => None,
            Err(mpsc::TryRecvError::Disconnected) => Some(Err(EngineError::new(
                ErrorCode::Io,
                "Gene-set resolution background worker disconnected",
            ))),
        };
        let Some(outcome) = outcome else {
            return;
        };
        let Some(task) = self.gene_set_inspector.resolve_task.take() else {
            return;
        };
        let elapsed = task.started.elapsed().as_secs_f64();
        match outcome {
            Ok(report) => {
                task.runtime_frame.update_phase("completed");
                let report_id = GentleEngine::gene_set_resolution_artifact_id(&report);
                let resolved_count = report.resolved_member_count;
                let unresolved_count = report.unresolved_member_count;
                self.refresh_gene_set_inspector_catalog(true);
                self.gene_set_inspector.selected_resolution_id = report_id.clone();
                self.gene_set_inspector.last_report = None;
                self.gene_set_inspector.child_reports.clear();
                self.reconcile_gene_set_inspector_bindings();
                self.gene_set_inspector.mode = GeneSetInspectorMode::Inspect;
                self.gene_set_inspector.resolve_status = format!(
                    "Resolved and persisted {resolved_count} member(s) in {elapsed:.1}s; {unresolved_count} unresolved."
                );
                self.gene_set_inspector.status = format!(
                    "Selected newly resolved gene set '{report_id}' for inspection and collection operations."
                );
            }
            Err(error) => {
                let stale = error.message.contains("became stale");
                if stale {
                    task.runtime_frame.cancel(error.message.clone());
                    self.gene_set_inspector.resolve_status = format!(
                        "Gene-set resolution became stale after {elapsed:.1}s; project state was not changed."
                    );
                } else {
                    task.runtime_frame.fail(error.message.clone());
                    self.gene_set_inspector.resolve_status = format!(
                        "Gene-set resolution failed after {elapsed:.1}s: {}",
                        error.message
                    );
                }
            }
        }
        ctx.request_repaint();
    }

    pub(super) fn has_active_gene_set_resolution_task(&self) -> bool {
        self.gene_set_inspector.resolve_task.is_some()
    }

    fn reconcile_gene_set_inspector_bindings(&mut self) {
        let member_ids = self
            .selected_gene_set_resolution()
            .map(|choice| {
                choice
                    .report
                    .resolved_members
                    .iter()
                    .map(|member| member.dedup_key.clone())
                    .collect::<BTreeSet<_>>()
            })
            .unwrap_or_default();
        let report_ids = self
            .gene_set_inspector
            .primer_reports
            .iter()
            .map(|report| report.report_id.as_str())
            .collect::<BTreeSet<_>>();
        self.gene_set_inspector
            .member_bindings
            .retain(|member_id, report_id| {
                member_ids.contains(member_id)
                    && (report_id.is_empty() || report_ids.contains(report_id.as_str()))
            });
        for member_id in member_ids {
            self.gene_set_inspector
                .member_bindings
                .entry(member_id)
                .or_default();
        }
    }

    fn gene_set_specificity_shell_command(&self) -> Result<ShellCommand, String> {
        let choice = self
            .selected_gene_set_resolution()
            .ok_or_else(|| "Select a persisted gene-set resolution first".to_string())?;
        if choice.report.resolved_members.is_empty() {
            return Err("The selected gene set has no resolved members".to_string());
        }
        let pair_rank = self
            .gene_set_inspector
            .pair_rank_1based
            .trim()
            .parse::<usize>()
            .map_err(|_| "Pair rank must be a positive integer".to_string())?;
        if pair_rank == 0 {
            return Err("Pair rank must be at least 1".to_string());
        }
        let target_genome_id = self.gene_set_inspector.target_genome_id.trim();
        if target_genome_id.is_empty() {
            return Err("Enter a prepared target genome or transcriptome id".to_string());
        }
        let reports_by_id = self
            .gene_set_inspector
            .primer_reports
            .iter()
            .map(|report| (report.report_id.as_str(), report))
            .collect::<BTreeMap<_, _>>();
        let mut bindings = Vec::with_capacity(choice.report.resolved_members.len());
        let mut missing_members = Vec::new();
        for member in &choice.report.resolved_members {
            let report_id = self
                .gene_set_inspector
                .member_bindings
                .get(&member.dedup_key)
                .map(String::as_str)
                .unwrap_or_default()
                .trim();
            if report_id.is_empty() {
                missing_members.push(member.symbol.clone());
                continue;
            }
            let report = reports_by_id.get(report_id).ok_or_else(|| {
                format!(
                    "Primer-design report '{}' selected for '{}' is no longer present",
                    report_id, member.symbol
                )
            })?;
            if pair_rank > report.pair_count {
                return Err(format!(
                    "Primer-design report '{}' for '{}' has {} pair(s), so rank {} is unavailable",
                    report_id, member.symbol, report.pair_count, pair_rank
                ));
            }
            bindings.push(PrimerSpecificityCollectionMemberBinding {
                stable_member_id: member.dedup_key.clone(),
                primer_report_id: report_id.to_string(),
            });
        }
        if !missing_members.is_empty() {
            return Err(format!(
                "Bind a primer-design report for every member; missing: {}",
                missing_members.join(", ")
            ));
        }
        let mut policy = PrimerSpecificityPolicy::default();
        policy.specificity_target_genome_id = Some(target_genome_id.to_string());
        let optional_path = |value: &str| {
            let trimmed = value.trim();
            (!trimmed.is_empty()).then(|| trimmed.to_string())
        };
        Ok(ShellCommand::CollectionsRunPrimerSpecificity {
            collection_subject: CollectionSubjectRef::GeneSetResolution {
                report_id: choice.report_id.clone(),
            },
            member_bindings: bindings,
            pair_rank: Some(pair_rank),
            pair_index: None,
            target_genome_id: target_genome_id.to_string(),
            policy,
            catalog_path: optional_path(&self.gene_set_inspector.catalog_path),
            cache_dir: optional_path(&self.gene_set_inspector.cache_dir),
            path: None,
        })
    }

    fn start_gene_set_specificity_task(&mut self) {
        if self.gene_set_inspector.task.is_some() {
            self.gene_set_inspector.status =
                "Gene-set primer specificity is already running".to_string();
            return;
        }
        let command = match self.gene_set_specificity_shell_command() {
            Ok(command) => command,
            Err(error) => {
                self.gene_set_inspector.status = error;
                return;
            }
        };
        let resolution_id = self.gene_set_inspector.selected_resolution_id.clone();
        let member_count = self
            .selected_gene_set_resolution()
            .map(|choice| choice.report.resolved_members.len())
            .unwrap_or_default();
        let job_id = self.alloc_background_job_id();
        let cancel_requested = Arc::new(AtomicBool::new(false));
        let worker_cancel = cancel_requested.clone();
        let engine = self.engine.clone();
        let (tx, receiver) = mpsc::channel::<Result<GeneSetSpecificityTaskResult, EngineError>>();
        let runtime_frame = Self::push_runtime_background_job_frame(
            BackgroundJobKind::PrimerSpecificityCollection,
            job_id,
            format!("resolution='{resolution_id}', members={member_count}"),
        );
        runtime_frame.update_phase("running");
        self.push_job_event(
            BackgroundJobKind::PrimerSpecificityCollection,
            BackgroundJobEventPhase::Started,
            Some(job_id),
            format!(
                "Collection primer specificity started: resolution='{}', members={}",
                resolution_id, member_count
            ),
        );
        std::thread::spawn(move || {
            let result = crate::background_engine::execute_on_engine_snapshot(
                &engine,
                move |snapshot| {
                    if worker_cancel.load(Ordering::Relaxed) {
                        return Err(EngineError::invalid_input(
                            "Collection primer specificity was discarded before execution",
                        ));
                    }
                    let shell_result = execute_shell_command(snapshot, &command)
                        .map_err(|message| EngineError::new(ErrorCode::InvalidInput, message))?;
                    if worker_cancel.load(Ordering::Relaxed) {
                        return Err(EngineError::invalid_input(
                            "Collection primer specificity finished after discard was requested; detached result was not committed",
                        ));
                    }
                    let report = serde_json::from_value::<CollectionOperationReport>(
                        shell_result
                            .output
                            .get("report")
                            .cloned()
                            .unwrap_or(serde_json::Value::Null),
                    )
                    .map_err(|error| {
                        EngineError::new(
                            ErrorCode::Internal,
                            format!(
                                "Collection specificity command returned an invalid report: {error}"
                            ),
                        )
                    })?;
                    let mut child_reports = BTreeMap::new();
                    for report_id in report
                        .per_member_status
                        .iter()
                        .flat_map(|row| row.produced_report_ids.iter())
                    {
                        if let Ok(child) = snapshot.get_primer_specificity_report(report_id) {
                            child_reports.insert(
                                report_id.clone(),
                                PrimerSpecificityChildPresentation {
                                    report_id: report_id.clone(),
                                    primary_seq_id: child.primary_seq_id,
                                    status: child.summary.status,
                                    specificity_pass: child.summary.specificity_pass,
                                    amplicon_count: child.summary.amplicon_count,
                                    failing_unintended_amplicon_count: child
                                        .summary
                                        .failing_unintended_amplicon_count,
                                    summary: child.summary.summary,
                                },
                            );
                        }
                    }
                    Ok(GeneSetSpecificityTaskResult {
                        report,
                        child_reports,
                    })
                },
            );
            let _ = tx.send(result);
        });
        self.gene_set_inspector.status = format!(
            "Running primer specificity for {member_count} gene-set member(s) in background..."
        );
        self.gene_set_inspector.task = Some(GeneSetSpecificityTask {
            job_id,
            resolution_id,
            started: Instant::now(),
            cancel_requested,
            runtime_frame,
            receiver,
        });
    }

    fn request_gene_set_specificity_discard(&mut self, origin: &str) {
        let Some(task) = self.gene_set_inspector.task.as_mut() else {
            self.gene_set_inspector.status =
                "No running gene-set specificity task to discard".to_string();
            return;
        };
        if task.cancel_requested.swap(true, Ordering::Relaxed) {
            self.gene_set_inspector.status =
                "Discard was already requested; waiting for the detached run to stop".to_string();
            return;
        }
        let job_id = task.job_id;
        task.runtime_frame.update_phase("discard_requested");
        self.gene_set_inspector.status =
            "Discard requested. An active BLAST process may finish before GENtle drops the detached result."
                .to_string();
        self.push_job_event(
            BackgroundJobKind::PrimerSpecificityCollection,
            BackgroundJobEventPhase::CancelRequested,
            Some(job_id),
            format!("Detached result discard requested from {origin}"),
        );
    }

    pub(super) fn poll_gene_set_specificity_task(&mut self, ctx: &egui::Context) {
        let Some(task) = self.gene_set_inspector.task.as_ref() else {
            return;
        };
        ctx.request_repaint_after(std::time::Duration::from_millis(100));
        let outcome = match task.receiver.try_recv() {
            Ok(result) => Some(result),
            Err(mpsc::TryRecvError::Empty) => None,
            Err(mpsc::TryRecvError::Disconnected) => Some(Err(EngineError::new(
                ErrorCode::Io,
                "Gene-set specificity background worker disconnected",
            ))),
        };
        let Some(outcome) = outcome else {
            return;
        };
        let Some(task) = self.gene_set_inspector.task.take() else {
            return;
        };
        let elapsed = task.started.elapsed().as_secs_f64();
        match outcome {
            Ok(result) => {
                let succeeded = result
                    .report
                    .per_member_status
                    .iter()
                    .filter(|row| row.outcome == CollectionMemberOutcome::Succeeded)
                    .count();
                let failed = result
                    .report
                    .per_member_status
                    .iter()
                    .filter(|row| row.outcome == CollectionMemberOutcome::Failed)
                    .count();
                self.gene_set_inspector.status = format!(
                    "Collection specificity completed in {elapsed:.1}s: {succeeded} executed, {failed} failed. Biological pass/fail is shown per child report."
                );
                self.gene_set_inspector.child_reports = result.child_reports;
                self.gene_set_inspector.last_report = Some(result.report);
                self.push_job_event(
                    BackgroundJobKind::PrimerSpecificityCollection,
                    BackgroundJobEventPhase::Completed,
                    Some(task.job_id),
                    format!(
                        "Collection primer specificity completed in {elapsed:.1}s for '{}': {succeeded} executed, {failed} failed",
                        task.resolution_id
                    ),
                );
            }
            Err(error) => {
                let stale = error.message.contains("became stale");
                let discarded = error.message.contains("discard");
                self.gene_set_inspector.status = if stale {
                    format!(
                        "Collection specificity produced a stale result after {elapsed:.1}s; project state was not changed."
                    )
                } else if discarded {
                    format!(
                        "Collection specificity result discarded after {elapsed:.1}s; project state was not changed."
                    )
                } else {
                    format!(
                        "Collection specificity failed after {elapsed:.1}s: {}",
                        error.message
                    )
                };
                self.push_job_event(
                    BackgroundJobKind::PrimerSpecificityCollection,
                    if stale || discarded {
                        BackgroundJobEventPhase::IgnoredStale
                    } else {
                        BackgroundJobEventPhase::Failed
                    },
                    Some(task.job_id),
                    self.gene_set_inspector.status.clone(),
                );
            }
        }
        ctx.request_repaint();
    }

    pub(super) fn has_active_gene_set_specificity_task(&self) -> bool {
        self.gene_set_inspector.task.is_some()
    }

    pub(super) fn render_gene_set_specificity_background_job(&mut self, ui: &mut Ui) {
        ui.separator();
        ui.strong("Gene-set Primer Specificity");
        let mut discard_clicked = false;
        if let Some(task) = self.gene_set_inspector.task.as_ref() {
            ui.horizontal_wrapped(|ui| {
                ui.add(egui::Spinner::new());
                ui.label(format!(
                    "Running #{} for '{}' ({:.1}s)",
                    task.job_id,
                    task.resolution_id,
                    task.started.elapsed().as_secs_f32()
                ));
                if task.cancel_requested.load(Ordering::Relaxed) {
                    ui.small("Discard requested...");
                } else if ui
                    .button("Discard result")
                    .on_hover_text(
                        "Stop waiting for the detached collection result. A running BLAST process may finish before GENtle can drop it.",
                    )
                    .clicked()
                {
                    discard_clicked = true;
                }
            });
        } else {
            ui.horizontal_wrapped(|ui| {
                ui.small("Idle");
                if ui
                    .button("Open inspector")
                    .on_hover_text(
                        "Open the binding-aware gene-set inspector to configure a collection specificity run",
                    )
                    .clicked()
                {
                    self.open_gene_set_inspector_dialog();
                }
            });
        }
        if discard_clicked {
            self.request_gene_set_specificity_discard("background jobs panel");
        }
        if !self.gene_set_inspector.status.trim().is_empty() {
            ui.small(self.gene_set_inspector.status.clone());
        }
    }

    fn render_gene_set_resolution_form(&mut self, ui: &mut Ui) {
        ui.label(
            "Resolve a catalog group, explicit list, local external mapping, genomic neighborhood, or deterministic random sample through the shared engine operation.",
        );
        ui.small(
            "External mappings use configured local gene-group catalogs; this form does not query an online ontology service.",
        );
        ui.separator();

        let mut choose_output_path = false;
        {
            let form = &mut self.gene_set_inspector.resolve_form;
            ui.horizontal_wrapped(|ui| {
                ui.label("Source kind");
                egui::ComboBox::from_id_salt("gene_set_resolve_source_kind")
                    .selected_text(form.source_kind.label())
                    .width(220.0)
                    .show_ui(ui, |ui| {
                        for source_kind in [
                            GeneSetResolveSourceKind::CatalogGroup,
                            GeneSetResolveSourceKind::ExplicitMembers,
                            GeneSetResolveSourceKind::ExternalMapping,
                            GeneSetResolveSourceKind::GenomicNeighbors,
                            GeneSetResolveSourceKind::Random,
                        ] {
                            ui.selectable_value(
                                &mut form.source_kind,
                                source_kind,
                                source_kind.label(),
                            );
                        }
                    });
            });

            ui.add_space(4.0);
            match form.source_kind {
                GeneSetResolveSourceKind::CatalogGroup => {
                    ui.horizontal_wrapped(|ui| {
                        ui.label("Group ID");
                        ui.add(
                            egui::TextEdit::singleline(&mut form.group_id)
                                .desired_width(460.0)
                                .hint_text("Catalog group ID or alias"),
                        )
                        .on_hover_text("Resolve one group from the configured gene-group catalog");
                    });
                }
                GeneSetResolveSourceKind::ExplicitMembers => {
                    ui.label("Members");
                    ui.add(
                        egui::TextEdit::multiline(&mut form.members)
                            .desired_rows(3)
                            .desired_width(f32::INFINITY)
                            .hint_text("TP53, TP73, TP63"),
                    )
                    .on_hover_text(
                        "Comma-separated gene symbols or stable IDs; empty entries are rejected",
                    );
                }
                GeneSetResolveSourceKind::ExternalMapping => {
                    ui.horizontal_wrapped(|ui| {
                        ui.label("Namespace");
                        egui::ComboBox::from_id_salt("gene_set_external_namespace_quick")
                            .selected_text(if form.external_namespace.trim().is_empty() {
                                "Choose common namespace"
                            } else {
                                form.external_namespace.as_str()
                            })
                            .width(170.0)
                            .show_ui(ui, |ui| {
                                for namespace in ["GO", "REACTOME", "KEGG"] {
                                    ui.selectable_value(
                                        &mut form.external_namespace,
                                        namespace.to_string(),
                                        namespace,
                                    );
                                }
                            });
                        ui.add(
                            egui::TextEdit::singleline(&mut form.external_namespace)
                                .desired_width(140.0)
                                .hint_text("custom namespace"),
                        )
                        .on_hover_text(
                            "Namespace prefix used by the configured local gene-group mapping",
                        );
                        ui.label("ID");
                        ui.add(
                            egui::TextEdit::singleline(&mut form.external_id)
                                .desired_width(260.0)
                                .hint_text("GO:0008152 or 0008152"),
                        )
                        .on_hover_text(
                            "External identifier; an included prefix must match the namespace",
                        );
                    });
                }
                GeneSetResolveSourceKind::GenomicNeighbors => {
                    ui.horizontal_wrapped(|ui| {
                        ui.label("Anchor gene");
                        ui.add(
                            egui::TextEdit::singleline(&mut form.anchor_gene)
                                .desired_width(240.0)
                                .hint_text("TP53"),
                        );
                        ui.label("Flank genes");
                        ui.add(
                            egui::TextEdit::singleline(&mut form.flank_genes)
                                .desired_width(70.0),
                        )
                        .on_hover_text("Positive count on each side; leave empty to use bp only");
                        ui.label("Flank bp");
                        ui.add(
                            egui::TextEdit::singleline(&mut form.flank_bp).desired_width(100.0),
                        )
                        .on_hover_text(
                            "Positive genomic distance on each side; leave empty to use gene count only",
                        );
                        ui.checkbox(&mut form.exclude_anchor, "Exclude anchor")
                            .on_hover_text("Omit the anchor gene from the resolved set");
                    });
                }
                GeneSetResolveSourceKind::Random => {
                    ui.horizontal_wrapped(|ui| {
                        ui.label("Sample size");
                        ui.add(
                            egui::TextEdit::singleline(&mut form.random_size)
                                .desired_width(80.0),
                        );
                        ui.label("Random seed");
                        ui.add(
                            egui::TextEdit::singleline(&mut form.random_seed)
                                .desired_width(120.0),
                        )
                        .on_hover_text(
                            "Non-negative deterministic seed; the same universe and seed reproduce the sample",
                        );
                    });
                    ui.horizontal_wrapped(|ui| {
                        ui.label("Exclude members");
                        ui.add(
                            egui::TextEdit::singleline(&mut form.exclude_members)
                                .desired_width(620.0)
                                .hint_text("optional comma-separated symbols or IDs"),
                        );
                    });
                }
            }

            ui.separator();
            ui.strong("Resolution context");
            ui.horizontal_wrapped(|ui| {
                ui.label("Genome ID");
                ui.add(
                    egui::TextEdit::singleline(&mut form.genome_id)
                        .desired_width(360.0)
                        .hint_text("optional prepared genome ID"),
                )
                .on_hover_text(
                    "Genome used for stable gene resolution, neighborhoods, and random sampling",
                );
                ui.checkbox(&mut form.allow_draft, "Allow draft")
                    .on_hover_text("Allow draft gene-group catalog entries");
                ui.checkbox(&mut form.allow_deprecated, "Allow deprecated")
                    .on_hover_text("Allow deprecated gene-group catalog entries");
            });
            egui::CollapsingHeader::new("Advanced local paths and JSON output")
                .default_open(false)
                .show(ui, |ui| {
                    egui::Grid::new("gene_set_resolve_advanced_paths")
                        .num_columns(3)
                        .spacing([8.0, 5.0])
                        .show(ui, |ui| {
                            ui.label("Gene-group catalog");
                            ui.add(
                                egui::TextEdit::singleline(&mut form.gene_group_catalog_path)
                                    .desired_width(560.0)
                                    .hint_text("empty = configured/default catalog"),
                            );
                            ui.end_row();
                            ui.label("Genome catalog");
                            ui.add(
                                egui::TextEdit::singleline(&mut form.genome_catalog_path)
                                    .desired_width(560.0)
                                    .hint_text("empty = configured/default catalog"),
                            );
                            ui.end_row();
                            ui.label("Genome cache");
                            ui.add(
                                egui::TextEdit::singleline(&mut form.cache_dir)
                                    .desired_width(560.0)
                                    .hint_text("empty = configured/default cache"),
                            );
                            ui.end_row();
                            ui.label("Output JSON");
                            ui.add(
                                egui::TextEdit::singleline(&mut form.output_path)
                                    .desired_width(480.0)
                                    .hint_text("optional report path"),
                            );
                            if ui
                                .button("Choose...")
                                .on_hover_text(
                                    "Choose an optional path for the engine-written resolution JSON",
                                )
                                .clicked()
                            {
                                choose_output_path = true;
                            }
                            ui.end_row();
                        });
                });
        }
        if choose_output_path
            && let Some(path) = rfd::FileDialog::new()
                .add_filter("JSON report", &["json"])
                .set_file_name("gene_set_resolution.json")
                .save_file()
        {
            self.gene_set_inspector.resolve_form.output_path = path.display().to_string();
        }

        let readiness = gene_set_resolve_form_to_operation(&self.gene_set_inspector.resolve_form);
        let running = self.gene_set_inspector.resolve_task.is_some();
        let mut resolve_clicked = false;
        ui.horizontal_wrapped(|ui| {
            if ui
                .add_enabled(
                    !running && readiness.is_ok(),
                    egui::Button::new("Resolve gene set"),
                )
                .on_hover_text(
                    readiness
                        .as_ref()
                        .map(|_| {
                            "Run the shared ResolveGeneSet engine operation in a detached project snapshot"
                        })
                        .unwrap_or_else(|error| error.as_str()),
                )
                .clicked()
            {
                resolve_clicked = true;
            }
            if let Some(task) = self.gene_set_inspector.resolve_task.as_ref() {
                ui.add(egui::Spinner::new());
                ui.label(format!(
                    "Resolving ({:.1}s)",
                    task.started.elapsed().as_secs_f32()
                ));
            }
        });
        if let Err(error) = &readiness {
            ui.small(format!("Not ready: {error}"));
        }
        if !self.gene_set_inspector.resolve_status.trim().is_empty() {
            ui.small(self.gene_set_inspector.resolve_status.clone());
        }
        if resolve_clicked {
            self.start_gene_set_resolution_task();
        }
    }

    fn render_gene_set_inspector_contents(&mut self, ui: &mut Ui) {
        let close_hover = Self::specialist_window_close_hover_text("Gene Set Inspector");
        if self.render_specialist_window_nav_with_close(ui, Some(("Close", close_hover.as_str()))) {
            self.gene_set_inspector.show_panel = false;
            return;
        }
        ui.heading("Gene Set Inspector");
        ui.horizontal_wrapped(|ui| {
            ui.selectable_value(
                &mut self.gene_set_inspector.mode,
                GeneSetInspectorMode::Resolve,
                "Resolve new set",
            )
            .on_hover_text("Author and resolve a new gene set through the shared engine operation");
            ui.selectable_value(
                &mut self.gene_set_inspector.mode,
                GeneSetInspectorMode::Inspect,
                "Inspect & run",
            )
            .on_hover_text("Inspect persisted resolutions and bind them to collection operations");
        });
        ui.separator();
        if self.gene_set_inspector.mode == GeneSetInspectorMode::Resolve {
            self.render_gene_set_resolution_form(ui);
            return;
        }
        ui.label(
            "Bind each logical gene-set member to one exact persisted primer-design report, then run the shared collection Map operation.",
        );
        ui.small(
            "Execution success means GENtle produced a specificity report. It does not mean the primer pair passed biological specificity.",
        );
        ui.separator();

        let resolution_choices = self.gene_set_inspector.resolutions.clone();
        let mut selection_changed = false;
        ui.horizontal_wrapped(|ui| {
            ui.label("Resolved gene set");
            egui::ComboBox::from_id_salt("gene_set_inspector_resolution")
                .selected_text(
                    resolution_choices
                        .iter()
                        .find(|choice| {
                            choice.report_id == self.gene_set_inspector.selected_resolution_id
                        })
                        .map(|choice| {
                            format!(
                                "{} ({} members)",
                                choice.report_id, choice.report.resolved_member_count
                            )
                        })
                        .unwrap_or_else(|| "No persisted gene-set resolution".to_string()),
                )
                .width(430.0)
                .show_ui(ui, |ui| {
                    for choice in &resolution_choices {
                        let label = format!(
                            "{} | resolved={} unresolved={} | genome={}",
                            choice.report_id,
                            choice.report.resolved_member_count,
                            choice.report.unresolved_member_count,
                            choice.report.genome_id.as_deref().unwrap_or("-")
                        );
                        selection_changed |= ui
                            .selectable_value(
                                &mut self.gene_set_inspector.selected_resolution_id,
                                choice.report_id.clone(),
                                label,
                            )
                            .changed();
                    }
                });
            if ui
                .button("Reload")
                .on_hover_text(
                    "Reload persisted gene-set resolutions and primer-design reports from project state",
                )
                .clicked()
            {
                self.refresh_gene_set_inspector_catalog(true);
            }
        });
        if selection_changed {
            self.gene_set_inspector.last_report = None;
            self.gene_set_inspector.child_reports.clear();
            self.reconcile_gene_set_inspector_bindings();
        }

        let selected = self.selected_gene_set_resolution().cloned();
        let Some(choice) = selected else {
            ui.separator();
            ui.label("No persisted gene-set resolutions are available.");
            ui.small(
                "Create one here or resolve and persist a set with `gene-sets resolve ...` in another shared adapter.",
            );
            if ui
                .button("Resolve a new set")
                .on_hover_text("Switch to the gene-set source authoring form")
                .clicked()
            {
                self.gene_set_inspector.mode = GeneSetInspectorMode::Resolve;
            }
            if !self.gene_set_inspector.status.trim().is_empty() {
                ui.separator();
                ui.small(self.gene_set_inspector.status.clone());
            }
            return;
        };

        ui.horizontal_wrapped(|ui| {
            ui.small(format!(
                "source={} | review={} | namespace={} | organism={} | resolved={} | unresolved={}",
                choice.report.request.source_kind_label(),
                gene_set_review_status_label(choice.report.review_status),
                choice.report.symbol_namespace.as_deref().unwrap_or("-"),
                choice.report.organism.as_deref().unwrap_or("-"),
                choice.report.resolved_member_count,
                choice.report.unresolved_member_count
            ));
        });
        if !choice.report.warnings.is_empty() {
            egui::CollapsingHeader::new(format!(
                "Resolution warnings ({})",
                choice.report.warnings.len()
            ))
            .show(ui, |ui| {
                for warning in &choice.report.warnings {
                    ui.small(format!("• {warning}"));
                }
            });
        }
        if !self.gene_set_inspector.resolve_status.trim().is_empty() {
            ui.small(self.gene_set_inspector.resolve_status.clone());
        }

        let resolution_json = serde_json::to_string_pretty(&choice.report).ok();
        let mut export_resolution = false;
        ui.horizontal_wrapped(|ui| {
            if ui
                .add_enabled(
                    resolution_json.is_some(),
                    egui::Button::new("Copy resolution JSON"),
                )
                .on_hover_text("Copy the selected portable GeneSetResolutionReport as JSON")
                .clicked()
                && let Some(json) = resolution_json.as_ref()
            {
                ui.ctx().copy_text(json.clone());
            }
            if ui
                .add_enabled(
                    resolution_json.is_some(),
                    egui::Button::new("Export resolution JSON..."),
                )
                .on_hover_text("Write the selected portable resolution report to a JSON file")
                .clicked()
            {
                export_resolution = true;
            }
        });
        if export_resolution
            && let Some(path) = rfd::FileDialog::new()
                .add_filter("JSON report", &["json"])
                .set_file_name("gene_set_resolution.json")
                .save_file()
            && let Some(json) = resolution_json.as_ref()
        {
            match std::fs::write(&path, format!("{json}\n")) {
                Ok(()) => {
                    self.gene_set_inspector.status =
                        format!("Exported gene-set resolution to '{}'", path.display());
                }
                Err(error) => {
                    self.gene_set_inspector.status = format!(
                        "Could not export gene-set resolution to '{}': {error}",
                        path.display()
                    );
                }
            }
        }

        egui::CollapsingHeader::new(format!(
            "Resolved members ({})",
            choice.report.resolved_members.len()
        ))
        .default_open(true)
        .show(ui, |ui| {
            egui::ScrollArea::both()
                .id_salt("gene_set_resolution_members")
                .max_height(220.0)
                .show(ui, |ui| {
                    egui::Grid::new("gene_set_resolution_member_grid")
                        .striped(true)
                        .spacing([14.0, 5.0])
                        .show(ui, |ui| {
                            ui.strong("Stable member ID");
                            ui.strong("Symbol");
                            ui.strong("Gene ID");
                            ui.strong("Member annotation");
                            ui.strong("Source provenance");
                            ui.end_row();
                            for member in &choice.report.resolved_members {
                                ui.monospace(&member.dedup_key);
                                ui.label(&member.symbol);
                                ui.monospace(member.gene_id.as_deref().unwrap_or("-"));
                                let annotation = match (
                                    member.member_status.as_deref(),
                                    member.confidence.as_deref(),
                                ) {
                                    (Some(status), Some(confidence)) => {
                                        format!("{status}; confidence={confidence}")
                                    }
                                    (Some(status), None) => status.to_string(),
                                    (None, Some(confidence)) => {
                                        format!("confidence={confidence}")
                                    }
                                    (None, None) => "-".to_string(),
                                };
                                ui.label(annotation).on_hover_text(
                                    "Member-owned status/confidence metadata; report review state is shown above",
                                );
                                ui.label(gene_set_provenance_summary(&member.provenance));
                                ui.end_row();
                            }
                        });
                });
        });
        if !choice.report.unresolved_members.is_empty() {
            egui::CollapsingHeader::new(format!(
                "Unresolved members ({})",
                choice.report.unresolved_members.len()
            ))
            .default_open(true)
            .show(ui, |ui| {
                egui::Grid::new("gene_set_resolution_unresolved_grid")
                    .striped(true)
                    .spacing([14.0, 5.0])
                    .show(ui, |ui| {
                        ui.strong("Query");
                        ui.strong("Reason");
                        ui.strong("Source");
                        ui.end_row();
                        for member in &choice.report.unresolved_members {
                            ui.monospace(&member.query);
                            ui.label(&member.reason);
                            ui.label(match member.source_id.as_deref() {
                                Some(source_id) => {
                                    format!("{}:{source_id}", member.source_kind)
                                }
                                None => member.source_kind.clone(),
                            });
                            ui.end_row();
                        }
                    });
            });
        }
        if !choice.report.provenance.is_empty() {
            egui::CollapsingHeader::new(format!(
                "Resolution provenance ({})",
                choice.report.provenance.len()
            ))
            .show(ui, |ui| {
                for row in &choice.report.provenance {
                    let mut details = gene_set_provenance_summary(std::slice::from_ref(row));
                    if let Some(path) = row.source_path.as_deref() {
                        details.push_str(&format!(" | {path}"));
                    }
                    if let Some(note) = row.note.as_deref() {
                        details.push_str(&format!(" | {note}"));
                    }
                    ui.small(details);
                }
            });
        }

        ui.separator();
        ui.strong("Primer-report bindings");
        ui.small(
            "Bindings are explicit and auditable. GENtle never guesses an assay from a gene symbol.",
        );
        let reports = self.gene_set_inspector.primer_reports.clone();
        ui.horizontal_wrapped(|ui| {
            ui.label(format!("Available primer reports: {}", reports.len()));
            if ui
                .button("Clear bindings")
                .on_hover_text("Clear every member-to-primer-report selection")
                .clicked()
            {
                for report_id in self.gene_set_inspector.member_bindings.values_mut() {
                    report_id.clear();
                }
            }
        });
        ui.horizontal(|ui| {
            ui.strong("Member");
            ui.add_space(150.0);
            ui.strong("Primer-design report");
        });
        egui::ScrollArea::vertical()
            .id_salt("gene_set_inspector_member_bindings")
            .max_height(300.0)
            .show_rows(
                ui,
                34.0,
                choice.report.resolved_members.len(),
                |ui, row_range| {
                    for row_index in row_range {
                        let member = &choice.report.resolved_members[row_index];
                        ui.horizontal(|ui| {
                            ui.vertical(|ui| {
                                ui.label(&member.symbol);
                                ui.small(
                                    member
                                        .gene_id
                                        .as_deref()
                                        .unwrap_or(member.dedup_key.as_str()),
                                );
                            });
                            ui.add_space(12.0);
                            let binding = self
                                .gene_set_inspector
                                .member_bindings
                                .entry(member.dedup_key.clone())
                                .or_default();
                            egui::ComboBox::from_id_salt((
                                "gene_set_member_primer_report",
                                &member.dedup_key,
                            ))
                            .selected_text(if binding.is_empty() {
                                "Select exact primer report..."
                            } else {
                                binding.as_str()
                            })
                            .width(520.0)
                            .show_ui(ui, |ui| {
                                ui.selectable_value(binding, String::new(), "Not selected");
                                for report in &reports {
                                    ui.selectable_value(
                                        binding,
                                        report.report_id.clone(),
                                        format!(
                                            "{} | template={} | pairs={} | {}",
                                            report.report_id,
                                            report.template,
                                            report.pair_count,
                                            report.backend_used
                                        ),
                                    );
                                }
                            });
                        });
                    }
                },
            );

        ui.separator();
        ui.strong("Specificity run");
        ui.horizontal_wrapped(|ui| {
            ui.label("Pair rank");
            ui.add(
                egui::TextEdit::singleline(&mut self.gene_set_inspector.pair_rank_1based)
                    .desired_width(70.0),
            )
            .on_hover_text("One-based accepted pair rank, applied to every bound report");
            ui.label("Prepared target genome/transcriptome");
            ui.add(
                egui::TextEdit::singleline(&mut self.gene_set_inspector.target_genome_id)
                    .desired_width(280.0)
                    .hint_text("Genome catalog ID"),
            )
            .on_hover_text(
                "Prepared local BLAST target id. GENtle uses the same specificity policy and database checks as the single-report command.",
            );
        });
        egui::CollapsingHeader::new("Advanced local paths")
            .default_open(false)
            .show(ui, |ui| {
                ui.horizontal_wrapped(|ui| {
                    ui.label("Genome catalog");
                    ui.add(
                        egui::TextEdit::singleline(&mut self.gene_set_inspector.catalog_path)
                            .desired_width(420.0)
                            .hint_text("empty = configured/default catalog"),
                    );
                });
                ui.horizontal_wrapped(|ui| {
                    ui.label("Genome cache");
                    ui.add(
                        egui::TextEdit::singleline(&mut self.gene_set_inspector.cache_dir)
                            .desired_width(420.0)
                            .hint_text("empty = configured/default cache"),
                    );
                });
            });

        let readiness = self.gene_set_specificity_shell_command();
        let running = self.gene_set_inspector.task.is_some();
        let mut run_clicked = false;
        let mut discard_clicked = false;
        ui.horizontal_wrapped(|ui| {
            if ui
                .add_enabled(
                    !running && readiness.is_ok(),
                    egui::Button::new("Run specificity on gene set"),
                )
                .on_hover_text(
                    readiness
                        .as_ref()
                        .map(|_| {
                            "Run the shared collections primer-specificity command in a background engine snapshot"
                        })
                        .unwrap_or_else(|error| error.as_str()),
                )
                .clicked()
            {
                run_clicked = true;
            }
            if let Some(task) = self.gene_set_inspector.task.as_ref() {
                ui.add(egui::Spinner::new());
                ui.label(format!(
                    "Running #{} ({:.1}s)",
                    task.job_id,
                    task.started.elapsed().as_secs_f32()
                ));
                if task.cancel_requested.load(Ordering::Relaxed) {
                    ui.small("Discard requested...");
                } else if ui
                    .button("Discard result")
                    .on_hover_text(
                        "Request that GENtle drop the detached result. An active BLAST subprocess may still need to finish.",
                    )
                    .clicked()
                {
                    discard_clicked = true;
                }
            }
        });
        if let Err(readiness_error) = &readiness {
            ui.small(format!("Not ready: {readiness_error}"));
        }
        if !self.gene_set_inspector.status.trim().is_empty() {
            ui.small(self.gene_set_inspector.status.clone());
        }

        let mut open_child: Option<(String, String)> = None;
        if let Some(report) = self.gene_set_inspector.last_report.clone() {
            ui.separator();
            let succeeded = report
                .per_member_status
                .iter()
                .filter(|row| row.outcome == CollectionMemberOutcome::Succeeded)
                .count();
            let failed = report
                .per_member_status
                .iter()
                .filter(|row| row.outcome == CollectionMemberOutcome::Failed)
                .count();
            ui.strong(format!(
                "Latest collection result: {} executed, {} failed",
                succeeded, failed
            ));
            ui.small(format!(
                "report_id={} | lift={:?} | fingerprint={}",
                report.report_id,
                report.lifting_mode,
                report.collection_membership_fingerprint_sha256
            ));
            if ui
                .button("Copy collection report JSON")
                .on_hover_text("Copy the portable collection-operation report")
                .clicked()
                && let Ok(json) = serde_json::to_string_pretty(&report)
            {
                ui.ctx().copy_text(json);
            }
            egui::ScrollArea::vertical()
                .id_salt("gene_set_inspector_result_rows")
                .max_height(280.0)
                .show_rows(
                    ui,
                    52.0,
                    report.per_member_status.len(),
                    |ui, row_range| {
                        for row_index in row_range {
                            let row = &report.per_member_status[row_index];
                            ui.horizontal_wrapped(|ui| {
                                ui.strong(
                                    row.member
                                        .gene_symbol
                                        .as_deref()
                                        .unwrap_or(row.member.stable_member_id.as_str()),
                                );
                                ui.label(format!("execution={:?}", row.outcome));
                                if let Some(error) = &row.error {
                                    ui.colored_label(
                                        ui.visuals().error_fg_color,
                                        error.message.clone(),
                                    );
                                }
                                for report_id in &row.produced_report_ids {
                                    if let Some(child) =
                                        self.gene_set_inspector.child_reports.get(report_id)
                                    {
                                        ui.label(format!(
                                            "biology={} | pass={} | products={} | failing off-targets={}",
                                            child.status,
                                            child.specificity_pass,
                                            child.amplicon_count,
                                            child.failing_unintended_amplicon_count
                                        ))
                                        .on_hover_text(child.summary.clone());
                                        if let Some(seq_id) = child.primary_seq_id.as_ref()
                                            && ui
                                                .button("Open report")
                                                .on_hover_text(
                                                    "Open this persisted child specificity report in PCR Designer",
                                                )
                                                .clicked()
                                        {
                                            open_child =
                                                Some((seq_id.clone(), child.report_id.clone()));
                                        }
                                    } else {
                                        ui.monospace(report_id);
                                    }
                                }
                            });
                        }
                    },
                );
            if !report.aggregate_warnings.is_empty() {
                egui::CollapsingHeader::new(format!(
                    "Aggregate warnings ({})",
                    report.aggregate_warnings.len()
                ))
                .show(ui, |ui| {
                    for warning in &report.aggregate_warnings {
                        ui.small(format!("• {warning}"));
                    }
                });
            }
        }

        if run_clicked {
            self.start_gene_set_specificity_task();
        }
        if discard_clicked {
            self.request_gene_set_specificity_discard("gene-set inspector");
        }
        if let Some((seq_id, report_id)) = open_child {
            self.open_sequence_window_for_primer_specificity_report(&seq_id, &report_id);
        }
    }

    pub(super) fn render_gene_set_inspector_dialog(&mut self, ctx: &egui::Context) {
        if !self.gene_set_inspector.show_panel {
            return;
        }
        self.refresh_gene_set_inspector_catalog(false);
        let mut open = self.gene_set_inspector.show_panel;
        let viewport_id = Self::gene_set_inspector_viewport_id();
        let spec = self.hosted_window_spec_for_viewport(
            "Gene Set Inspector",
            Self::hosted_gene_set_inspector_window_id(),
            viewport_id,
            Vec2::new(980.0, 760.0),
            Vec2::new(680.0, 460.0),
        );
        if ctx.embed_viewports() {
            crate::egui_compat::show_hosted_window(ctx, &spec, &mut open, |ui| {
                self.render_gene_set_inspector_contents(ui);
            });
            self.clear_viewport_foreground_request_after_render(viewport_id);
            self.finalize_viewport_open_probe(viewport_id, "Gene Set Inspector");
            self.gene_set_inspector.show_panel = open && self.gene_set_inspector.show_panel;
            return;
        }

        let builder = crate::egui_compat::viewport_builder_for_hosted_window(&spec);
        ctx.show_viewport_immediate(viewport_id, builder, |ctx, class| {
            self.note_viewport_focus_if_active(ctx, viewport_id);
            if class == egui::ViewportClass::EmbeddedWindow {
                crate::egui_compat::show_hosted_window(&mut *ctx, &spec, &mut open, |ui| {
                    self.render_gene_set_inspector_contents(ui);
                });
            } else {
                crate::egui_compat::show_central_panel(
                    &mut *ctx,
                    egui::CentralPanel::default(),
                    |ui| self.render_gene_set_inspector_contents(ui),
                );
                if Self::viewport_close_requested_or_shortcut(ctx) {
                    open = false;
                }
            }
        });
        self.gene_set_inspector.show_panel = open && self.gene_set_inspector.show_panel;
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{engine::GeneSetResolvedMember, engine_shell::parse_shell_line};

    fn resolve_operation_from_shell(command: ShellCommand) -> Operation {
        match command {
            ShellCommand::GeneSetsResolve {
                source,
                genome_id,
                gene_group_catalog_path,
                genome_catalog_path,
                cache_dir,
                allow_draft,
                allow_deprecated,
                output,
            } => Operation::ResolveGeneSet {
                source,
                genome_id,
                gene_group_catalog_path,
                genome_catalog_path,
                cache_dir,
                allow_draft,
                allow_deprecated,
                path: output,
            },
            other => panic!("unexpected command: {other:?}"),
        }
    }

    fn assert_resolve_form_matches_shell(form: GeneSetResolveFormState, shell_line: &str) {
        let gui_operation = gene_set_resolve_form_to_operation(&form).expect("GUI operation");
        let shell_operation = resolve_operation_from_shell(
            parse_shell_line(shell_line).expect("shared shell parser command"),
        );
        assert_eq!(
            serde_json::to_value(gui_operation).expect("serialize GUI operation"),
            serde_json::to_value(shell_operation).expect("serialize shell operation")
        );
    }

    #[test]
    fn gene_set_resolve_parity_catalog_group() {
        assert_resolve_form_matches_shell(
            GeneSetResolveFormState {
                group_id: "my_group".to_string(),
                ..GeneSetResolveFormState::default()
            },
            "gene-sets resolve my_group",
        );
    }

    #[test]
    fn gene_set_resolve_parity_explicit_members() {
        assert_resolve_form_matches_shell(
            GeneSetResolveFormState {
                source_kind: GeneSetResolveSourceKind::ExplicitMembers,
                members: "A, B, C".to_string(),
                ..GeneSetResolveFormState::default()
            },
            "gene-sets resolve --members A,B,C",
        );
    }

    #[test]
    fn gene_set_resolve_parity_external_mapping() {
        assert_resolve_form_matches_shell(
            GeneSetResolveFormState {
                source_kind: GeneSetResolveSourceKind::ExternalMapping,
                external_namespace: "go".to_string(),
                external_id: "GO:0008152".to_string(),
                ..GeneSetResolveFormState::default()
            },
            "gene-sets resolve --go GO:0008152",
        );
    }

    #[test]
    fn gene_set_resolve_parity_external_mapping_preserves_prefixed_id() {
        assert_resolve_form_matches_shell(
            GeneSetResolveFormState {
                source_kind: GeneSetResolveSourceKind::ExternalMapping,
                external_namespace: "reactome".to_string(),
                external_id: "Reactome:R-HSA-123".to_string(),
                ..GeneSetResolveFormState::default()
            },
            "gene-sets resolve --external-mapping Reactome:R-HSA-123",
        );
    }

    #[test]
    fn gene_set_resolve_parity_genomic_neighbors() {
        assert_resolve_form_matches_shell(
            GeneSetResolveFormState {
                source_kind: GeneSetResolveSourceKind::GenomicNeighbors,
                anchor_gene: "TP53".to_string(),
                flank_genes: "5".to_string(),
                flank_bp: String::new(),
                ..GeneSetResolveFormState::default()
            },
            "gene-sets resolve --neighbors TP53 --flank-genes 5",
        );
    }

    #[test]
    fn gene_set_resolve_parity_random() {
        assert_resolve_form_matches_shell(
            GeneSetResolveFormState {
                source_kind: GeneSetResolveSourceKind::Random,
                random_size: "10".to_string(),
                random_seed: "42".to_string(),
                ..GeneSetResolveFormState::default()
            },
            "gene-sets resolve --random-size 10 --seed 42",
        );
    }

    #[test]
    fn gene_set_resolve_parity_random_uses_shell_default_seed() {
        assert_resolve_form_matches_shell(
            GeneSetResolveFormState {
                source_kind: GeneSetResolveSourceKind::Random,
                random_size: "10".to_string(),
                random_seed: String::new(),
                ..GeneSetResolveFormState::default()
            },
            "gene-sets resolve --random-size 10",
        );
    }

    #[test]
    fn gene_set_resolve_form_rejects_incomplete_or_invalid_sources() {
        assert!(
            gene_set_resolve_form_to_operation(&GeneSetResolveFormState::default())
                .expect_err("empty group")
                .contains("Group ID")
        );
        assert!(
            gene_set_resolve_form_to_operation(&GeneSetResolveFormState {
                source_kind: GeneSetResolveSourceKind::ExplicitMembers,
                members: "A,,B".to_string(),
                ..GeneSetResolveFormState::default()
            })
            .expect_err("empty member token")
            .contains("empty entry")
        );
        assert!(
            gene_set_resolve_form_to_operation(&GeneSetResolveFormState {
                source_kind: GeneSetResolveSourceKind::ExternalMapping,
                external_namespace: "GO".to_string(),
                external_id: "KEGG:123".to_string(),
                ..GeneSetResolveFormState::default()
            })
            .expect_err("namespace mismatch")
            .contains("does not match")
        );
        assert!(
            gene_set_resolve_form_to_operation(&GeneSetResolveFormState {
                source_kind: GeneSetResolveSourceKind::GenomicNeighbors,
                anchor_gene: "TP53".to_string(),
                flank_genes: String::new(),
                flank_bp: String::new(),
                ..GeneSetResolveFormState::default()
            })
            .expect_err("missing flank")
            .contains("flank genes")
        );
        assert!(
            gene_set_resolve_form_to_operation(&GeneSetResolveFormState {
                source_kind: GeneSetResolveSourceKind::Random,
                random_size: "0".to_string(),
                ..GeneSetResolveFormState::default()
            })
            .expect_err("zero sample")
            .contains("at least 1")
        );
    }

    #[test]
    fn gene_set_resolve_task_persists_and_selects_exact_report() {
        let mut app = GENtleApp::default();
        app.gene_set_inspector.resolve_form = GeneSetResolveFormState {
            source_kind: GeneSetResolveSourceKind::ExplicitMembers,
            members: "TP53".to_string(),
            ..GeneSetResolveFormState::default()
        };
        app.start_gene_set_resolution_task();
        let ctx = egui::Context::default();
        for _ in 0..200 {
            app.poll_gene_set_resolution_task(&ctx);
            if app.gene_set_inspector.resolve_task.is_none() {
                break;
            }
            std::thread::sleep(std::time::Duration::from_millis(5));
        }

        assert!(app.gene_set_inspector.resolve_task.is_none());
        assert_eq!(app.gene_set_inspector.mode, GeneSetInspectorMode::Inspect);
        let selected = app
            .selected_gene_set_resolution()
            .expect("persisted selected resolution");
        assert_eq!(
            selected.report.request,
            GeneSetRequest::ExplicitMembers {
                members: vec!["TP53".to_string()]
            }
        );
        assert!(app.gene_set_inspector.resolve_status.contains("persisted"));
    }

    fn seeded_app() -> GENtleApp {
        let mut app = GENtleApp::default();
        app.gene_set_inspector.resolutions = vec![GeneSetResolutionChoice {
            report_id: "resolution:genes".to_string(),
            report: GeneSetResolutionReport {
                review_status: GeneSetResolutionReviewStatus::Reviewed,
                resolved_member_count: 2,
                resolved_members: vec![
                    GeneSetResolvedMember {
                        dedup_key: "GENE1".to_string(),
                        symbol: "GENE1".to_string(),
                        ..GeneSetResolvedMember::default()
                    },
                    GeneSetResolvedMember {
                        dedup_key: "GENE2".to_string(),
                        symbol: "GENE2".to_string(),
                        ..GeneSetResolvedMember::default()
                    },
                ],
                ..GeneSetResolutionReport::default()
            },
        }];
        app.gene_set_inspector.selected_resolution_id = "resolution:genes".to_string();
        app.gene_set_inspector.primer_reports = vec![
            PrimerDesignReportSummary {
                report_id: "report_1".to_string(),
                template: "seq_1".to_string(),
                pair_count: 2,
                ..PrimerDesignReportSummary::default()
            },
            PrimerDesignReportSummary {
                report_id: "report_2".to_string(),
                template: "seq_2".to_string(),
                pair_count: 1,
                ..PrimerDesignReportSummary::default()
            },
        ];
        app.gene_set_inspector
            .member_bindings
            .insert("GENE1".to_string(), "report_1".to_string());
        app.gene_set_inspector
            .member_bindings
            .insert("GENE2".to_string(), "report_2".to_string());
        app.gene_set_inspector.target_genome_id = "GRCh38".to_string();
        app
    }

    fn command_projection(command: &ShellCommand) -> serde_json::Value {
        match command {
            ShellCommand::CollectionsRunPrimerSpecificity {
                collection_subject,
                member_bindings,
                pair_rank,
                pair_index,
                target_genome_id,
                policy,
                catalog_path,
                cache_dir,
                path,
                ..
            } => serde_json::json!({
                "collection_subject": collection_subject,
                "member_bindings": member_bindings,
                "pair_rank": pair_rank,
                "pair_index": pair_index,
                "target_genome_id": target_genome_id,
                "policy": policy,
                "catalog_path": catalog_path,
                "cache_dir": cache_dir,
                "path": path,
            }),
            other => panic!("unexpected command: {other:?}"),
        }
    }

    #[test]
    fn gene_set_inspector_command_matches_shared_shell_parser() {
        let app = seeded_app();
        let gui_command = app
            .gene_set_specificity_shell_command()
            .expect("GUI command");
        let parsed = parse_shell_line(
            "collections run primer-specificity resolution:genes \
             --member-report GENE1=report_1 \
             --member-report GENE2=report_2 \
             --pair-rank 1 --target-genome GRCh38",
        )
        .expect("shared shell command");

        assert_eq!(
            command_projection(&gui_command),
            command_projection(&parsed)
        );
    }

    #[test]
    fn gene_set_inspector_requires_complete_explicit_bindings_and_valid_rank() {
        let mut app = seeded_app();
        app.gene_set_inspector.member_bindings.remove("GENE2");
        assert!(
            app.gene_set_specificity_shell_command()
                .expect_err("missing binding")
                .contains("missing: GENE2")
        );

        app.gene_set_inspector
            .member_bindings
            .insert("GENE2".to_string(), "report_2".to_string());
        app.gene_set_inspector.pair_rank_1based = "2".to_string();
        assert!(
            app.gene_set_specificity_shell_command()
                .expect_err("rank unavailable for report_2")
                .contains("rank 2 is unavailable")
        );
    }

    #[test]
    fn gene_set_specificity_task_completion_preserves_execution_and_biology_states() {
        let mut app = GENtleApp::default();
        let (tx, receiver) = mpsc::channel();
        let runtime_frame = GENtleApp::push_runtime_background_job_frame(
            BackgroundJobKind::PrimerSpecificityCollection,
            9,
            "test",
        );
        app.gene_set_inspector.task = Some(GeneSetSpecificityTask {
            job_id: 9,
            resolution_id: "resolution:genes".to_string(),
            started: Instant::now(),
            cancel_requested: Arc::new(AtomicBool::new(false)),
            runtime_frame,
            receiver,
        });
        let child = PrimerSpecificityChildPresentation {
            report_id: "specificity:1".to_string(),
            status: "fail".to_string(),
            specificity_pass: false,
            summary: "Synthetic biological failure".to_string(),
            ..PrimerSpecificityChildPresentation::default()
        };
        tx.send(Ok(GeneSetSpecificityTaskResult {
            report: CollectionOperationReport {
                report_id: "collection:1".to_string(),
                per_member_status: vec![crate::engine::CollectionMemberStatusRow {
                    outcome: CollectionMemberOutcome::Succeeded,
                    produced_report_ids: vec!["specificity:1".to_string()],
                    ..crate::engine::CollectionMemberStatusRow::default()
                }],
                ..CollectionOperationReport::default()
            },
            child_reports: BTreeMap::from([("specificity:1".to_string(), child)]),
        }))
        .expect("send result");

        app.poll_gene_set_specificity_task(&egui::Context::default());

        assert!(app.gene_set_inspector.task.is_none());
        assert!(
            app.gene_set_inspector
                .status
                .contains("1 executed, 0 failed")
        );
        assert_eq!(
            app.gene_set_inspector.child_reports["specificity:1"].status,
            "fail"
        );
        assert!(
            !app.gene_set_inspector.child_reports["specificity:1"].specificity_pass,
            "execution success must not be conflated with biological pass"
        );
    }
}
