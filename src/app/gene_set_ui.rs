//! Gene-set source authoring, resolution inspection, and collection operations.
//!
//! This GUI is deliberately thin: it constructs the shared `ResolveGeneSet`
//! operation, inspects persisted portable reports, requires explicit
//! primer-report bindings, and executes typed collection adapters through the
//! shared shell/engine path on detached engine snapshots.

use super::collection_operations_ui::{
    CollectionLauncherAdapter, CollectionLauncherReadiness, CollectionLauncherRow,
    collection_launcher_rows,
};
use super::*;
use crate::{
    engine::{
        CollectionContextRequirement, CollectionDigestReport, CollectionLiftRejectionReason,
        CollectionMemberOutcome, CollectionMemberRef, CollectionOperationReport,
        CollectionRestrictionSiteScanReport, CollectionSubjectKind, CollectionSubjectRef,
        CollectionTfbsHitScanReport, ContainerKind, DEFAULT_PROMOTER_WINDOW_DOWNSTREAM_BP,
        DEFAULT_PROMOTER_WINDOW_UPSTREAM_BP, DigestCollectionMemberBinding,
        GeneSetCohortRelationship, GeneSetPoolCreationReport, GeneSetPoolMemberBinding,
        GeneSetPromoterCohortReport, GeneSetProvenanceRow, GeneSetRequest, GeneSetResolutionReport,
        GeneSetResolutionReviewStatus, PrimerDesignReportSummary,
        PrimerSpecificityCollectionMemberBinding, PrimerSpecificityPolicy,
        RestrictionSiteScanCollectionMemberBinding, TfThresholdOverride,
        TfbsHitScanCollectionMemberBinding, homogeneous_collection_biological_context,
        validate_collection_context_target_genome,
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

fn gene_set_relationship_label(relationship: GeneSetCohortRelationship) -> &'static str {
    match relationship {
        GeneSetCohortRelationship::Unspecified => "Unspecified",
        GeneSetCohortRelationship::Manual => "Manual",
        GeneSetCohortRelationship::CoRegulated => "Co-regulated",
        GeneSetCohortRelationship::AntiCoRegulated => "Anti-co-regulated",
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

fn parse_optional_usize(value: &str, label: &str) -> Result<Option<usize>, String> {
    let value = value.trim();
    if value.is_empty() {
        return Ok(None);
    }
    let parsed = value
        .parse::<usize>()
        .map_err(|_| format!("{label} must be a non-negative integer"))?;
    Ok(Some(parsed))
}

fn parse_optional_f64(value: &str, label: &str) -> Result<Option<f64>, String> {
    let value = value.trim();
    if value.is_empty() {
        return Ok(None);
    }
    value
        .parse::<f64>()
        .map(Some)
        .map_err(|_| format!("{label} must be a number"))
}

fn parse_tfbs_threshold_overrides(
    min_bits: &str,
    min_quantiles: &str,
) -> Result<Vec<TfThresholdOverride>, String> {
    let mut overrides = BTreeMap::<String, TfThresholdOverride>::new();
    for (raw, label, quantile) in [
        (min_bits, "Per-TF minimum bits", false),
        (min_quantiles, "Per-TF minimum quantiles", true),
    ] {
        for entry in raw
            .split([',', '\n'])
            .map(str::trim)
            .filter(|entry| !entry.is_empty())
        {
            let (tf, value) = entry
                .split_once('=')
                .ok_or_else(|| format!("{label} entries require TF=VALUE"))?;
            let tf = tf.trim();
            let value = value
                .trim()
                .parse::<f64>()
                .map_err(|_| format!("{label} value for '{tf}' must be a number"))?;
            if tf.is_empty() {
                return Err(format!("{label} entries require a non-empty TF name"));
            }
            let row = overrides
                .entry(tf.to_string())
                .or_insert_with(|| TfThresholdOverride {
                    tf: tf.to_string(),
                    min_llr_bits: None,
                    min_llr_quantile: None,
                });
            if quantile {
                row.min_llr_quantile = Some(value);
            } else {
                row.min_llr_bits = Some(value);
            }
        }
    }
    Ok(overrides.into_values().collect())
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
            let flank_gene_count = parse_optional_usize(&form.flank_genes, "Flank genes")?;
            let flank_bp = parse_optional_usize(&form.flank_bp, "Flank bp")?;
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
            let count = parse_optional_usize(&form.random_size, "Sample size")?
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

#[derive(Clone, Debug)]
struct GeneSetSpecificityFormState {
    member_bindings: BTreeMap<String, String>,
    pair_rank_1based: String,
    target_genome_id: String,
    catalog_path: String,
    cache_dir: String,
}

impl Default for GeneSetSpecificityFormState {
    fn default() -> Self {
        Self {
            member_bindings: BTreeMap::new(),
            pair_rank_1based: "1".to_string(),
            target_genome_id: String::new(),
            catalog_path: String::new(),
            cache_dir: String::new(),
        }
    }
}

#[derive(Clone, Debug)]
struct GeneSetPromoterFormState {
    genome_id: String,
    relationship: GeneSetCohortRelationship,
    upstream_bp: String,
    downstream_bp: String,
    gene_group_catalog_path: String,
    genome_catalog_path: String,
    cache_dir: String,
    output_path: String,
}

#[derive(Clone, Debug)]
struct GeneSetRestrictionScanFormState {
    member_bindings: BTreeMap<String, String>,
    enzymes: String,
    max_sites_per_enzyme: String,
    include_cut_geometry: bool,
}

#[derive(Clone, Debug, Default)]
struct GeneSetTfbsScanFormState {
    member_bindings: BTreeMap<String, String>,
    motifs: String,
    min_llr_bits: String,
    min_llr_quantile: String,
    per_tf_min_llr_bits: String,
    per_tf_min_llr_quantile: String,
    max_hits_per_member: String,
}

#[derive(Clone, Debug, Default)]
struct GeneSetDigestFormState {
    member_bindings: BTreeMap<String, String>,
    enzymes: String,
    output_prefix: String,
    apply_change: bool,
    expected_plan_fingerprint_sha256: String,
}

#[derive(Clone, Debug, Default)]
struct GeneSetPoolFormState {
    member_bindings: BTreeMap<String, String>,
    output_prefix: String,
    container_name: String,
    apply_change: bool,
    expected_plan_fingerprint_sha256: String,
}

impl Default for GeneSetRestrictionScanFormState {
    fn default() -> Self {
        Self {
            member_bindings: BTreeMap::new(),
            enzymes: String::new(),
            max_sites_per_enzyme: String::new(),
            include_cut_geometry: true,
        }
    }
}

impl Default for GeneSetPromoterFormState {
    fn default() -> Self {
        Self {
            genome_id: String::new(),
            relationship: GeneSetCohortRelationship::Unspecified,
            upstream_bp: DEFAULT_PROMOTER_WINDOW_UPSTREAM_BP.to_string(),
            downstream_bp: DEFAULT_PROMOTER_WINDOW_DOWNSTREAM_BP.to_string(),
            gene_group_catalog_path: String::new(),
            genome_catalog_path: String::new(),
            cache_dir: String::new(),
            output_path: String::new(),
        }
    }
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

#[derive(Clone, Debug)]
enum CollectionOperationPayload {
    PrimerSpecificity {
        child_reports: BTreeMap<String, PrimerSpecificityChildPresentation>,
    },
    RestrictionScan(Box<CollectionRestrictionSiteScanReport>),
    TfbsScan(Box<CollectionTfbsHitScanReport>),
    Digest(Box<CollectionDigestReport>),
    CreatePool(Box<GeneSetPoolCreationReport>),
    PromoterCohort(Box<GeneSetPromoterCohortReport>),
}

#[derive(Clone, Debug)]
struct CollectionOperationTaskResult {
    adapter: CollectionLauncherAdapter,
    report: CollectionOperationReport,
    payload: CollectionOperationPayload,
}

struct GeneSetCollectionOperationTask {
    job_id: u64,
    adapter: CollectionLauncherAdapter,
    resolution_id: String,
    started: Instant,
    cancel_requested: Arc<AtomicBool>,
    runtime_frame: RuntimeStatusGuard,
    receiver: mpsc::Receiver<Result<CollectionOperationTaskResult, EngineError>>,
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
    selected_operation: CollectionLauncherAdapter,
    specificity_form: GeneSetSpecificityFormState,
    restriction_scan_form: GeneSetRestrictionScanFormState,
    tfbs_scan_form: GeneSetTfbsScanFormState,
    digest_form: GeneSetDigestFormState,
    pool_form: GeneSetPoolFormState,
    promoter_form: GeneSetPromoterFormState,
    status: String,
    resolve_form: GeneSetResolveFormState,
    resolve_status: String,
    resolve_task: Option<GeneSetResolveTask>,
    task: Option<GeneSetCollectionOperationTask>,
    last_result: Option<CollectionOperationTaskResult>,
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
            selected_operation: CollectionLauncherAdapter::default(),
            specificity_form: GeneSetSpecificityFormState::default(),
            restriction_scan_form: GeneSetRestrictionScanFormState::default(),
            tfbs_scan_form: GeneSetTfbsScanFormState::default(),
            digest_form: GeneSetDigestFormState::default(),
            pool_form: GeneSetPoolFormState::default(),
            promoter_form: GeneSetPromoterFormState::default(),
            status: String::new(),
            resolve_form: GeneSetResolveFormState::default(),
            resolve_status: String::new(),
            resolve_task: None,
            task: None,
            last_result: None,
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
        self.seed_gene_set_promoter_form_from_selection(false);
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
        let selection_changed = !self
            .gene_set_inspector
            .resolutions
            .iter()
            .any(|choice| choice.report_id == self.gene_set_inspector.selected_resolution_id);
        if selection_changed {
            self.gene_set_inspector.selected_resolution_id = self
                .gene_set_inspector
                .resolutions
                .first()
                .map(|choice| choice.report_id.clone())
                .unwrap_or_default();
            self.gene_set_inspector.last_result = None;
        }
        self.reconcile_gene_set_inspector_bindings();
        if selection_changed {
            self.seed_gene_set_promoter_form_from_selection(true);
        }
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
                self.gene_set_inspector.last_result = None;
                self.reconcile_gene_set_inspector_bindings();
                self.seed_gene_set_promoter_form_from_selection(true);
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

    fn gene_set_pool_source_container_choices(&self) -> Vec<(String, String)> {
        let mut choices = self
            .engine
            .read()
            .ok()
            .map(|engine| {
                engine
                    .state()
                    .container_state
                    .containers
                    .values()
                    .filter(|container| {
                        matches!(&container.kind, ContainerKind::Singleton)
                            && container.declared_contents_exclusive
                            && container.members.len() == 1
                            && engine.state().sequences.contains_key(&container.members[0])
                    })
                    .map(|container| {
                        let label = format!(
                            "{} | {}{}",
                            container.container_id,
                            container.members[0],
                            container
                                .name
                                .as_deref()
                                .filter(|name| !name.trim().is_empty())
                                .map(|name| format!(" | {name}"))
                                .unwrap_or_default()
                        );
                        (container.container_id.clone(), label)
                    })
                    .collect::<Vec<_>>()
            })
            .unwrap_or_default();
        choices.sort_by(|left, right| left.0.cmp(&right.0));
        choices
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
            .specificity_form
            .member_bindings
            .retain(|member_id, report_id| {
                member_ids.contains(member_id)
                    && (report_id.is_empty() || report_ids.contains(report_id.as_str()))
            });
        for member_id in member_ids {
            self.gene_set_inspector
                .specificity_form
                .member_bindings
                .entry(member_id)
                .or_default();
        }
        let loaded_seq_ids = self
            .engine
            .read()
            .ok()
            .map(|engine| {
                engine
                    .state()
                    .sequences
                    .keys()
                    .cloned()
                    .collect::<BTreeSet<_>>()
            })
            .unwrap_or_default();
        let selected_member_ids = self
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
        self.gene_set_inspector
            .restriction_scan_form
            .member_bindings
            .retain(|member_id, seq_id| {
                selected_member_ids.contains(member_id)
                    && (seq_id.is_empty() || loaded_seq_ids.contains(seq_id))
            });
        for member_id in &selected_member_ids {
            self.gene_set_inspector
                .restriction_scan_form
                .member_bindings
                .entry(member_id.clone())
                .or_default();
        }
        let pool_source_container_ids = self
            .gene_set_pool_source_container_choices()
            .into_iter()
            .map(|(container_id, _)| container_id)
            .collect::<BTreeSet<_>>();
        self.gene_set_inspector
            .pool_form
            .member_bindings
            .retain(|member_id, container_id| {
                selected_member_ids.contains(member_id)
                    && (container_id.is_empty() || pool_source_container_ids.contains(container_id))
            });
        for member_id in &selected_member_ids {
            self.gene_set_inspector
                .pool_form
                .member_bindings
                .entry(member_id.clone())
                .or_default();
        }
        self.gene_set_inspector
            .digest_form
            .member_bindings
            .retain(|member_id, seq_id| {
                selected_member_ids.contains(member_id)
                    && (seq_id.is_empty() || loaded_seq_ids.contains(seq_id))
            });
        for member_id in &selected_member_ids {
            self.gene_set_inspector
                .digest_form
                .member_bindings
                .entry(member_id.clone())
                .or_default();
        }
        self.gene_set_inspector
            .tfbs_scan_form
            .member_bindings
            .retain(|member_id, seq_id| {
                selected_member_ids.contains(member_id)
                    && (seq_id.is_empty() || loaded_seq_ids.contains(seq_id))
            });
        for member_id in selected_member_ids {
            self.gene_set_inspector
                .tfbs_scan_form
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
            .specificity_form
            .pair_rank_1based
            .trim()
            .parse::<usize>()
            .map_err(|_| "Pair rank must be a positive integer".to_string())?;
        if pair_rank == 0 {
            return Err("Pair rank must be at least 1".to_string());
        }
        let target_genome_id = self
            .gene_set_inspector
            .specificity_form
            .target_genome_id
            .trim();
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
                .specificity_form
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
            catalog_path: optional_path(&self.gene_set_inspector.specificity_form.catalog_path),
            cache_dir: optional_path(&self.gene_set_inspector.specificity_form.cache_dir),
            path: None,
        })
    }

    fn gene_set_restriction_scan_shell_command(&self) -> Result<ShellCommand, String> {
        let choice = self
            .selected_gene_set_resolution()
            .ok_or_else(|| "Select a persisted gene-set resolution first".to_string())?;
        if choice.report.resolved_members.is_empty() {
            return Err("The selected gene set has no resolved members".to_string());
        }
        let mut member_bindings = Vec::with_capacity(choice.report.resolved_members.len());
        let mut missing_members = Vec::new();
        for member in &choice.report.resolved_members {
            let seq_id = self
                .gene_set_inspector
                .restriction_scan_form
                .member_bindings
                .get(&member.dedup_key)
                .map(String::as_str)
                .unwrap_or_default()
                .trim();
            if seq_id.is_empty() {
                missing_members.push(member.symbol.clone());
            } else {
                member_bindings.push(RestrictionSiteScanCollectionMemberBinding {
                    stable_member_id: member.dedup_key.clone(),
                    seq_id: seq_id.to_string(),
                });
            }
        }
        if !missing_members.is_empty() {
            return Err(format!(
                "Bind one loaded DNA sequence for every member; missing: {}",
                missing_members.join(", ")
            ));
        }
        let enzymes = self
            .gene_set_inspector
            .restriction_scan_form
            .enzymes
            .split([',', '\n'])
            .map(str::trim)
            .filter(|value| !value.is_empty())
            .map(str::to_string)
            .collect::<Vec<_>>();
        let max_sites_per_enzyme = match self
            .gene_set_inspector
            .restriction_scan_form
            .max_sites_per_enzyme
            .trim()
        {
            "" => None,
            value => Some(value.parse::<usize>().map_err(|_| {
                "Maximum sites per enzyme must be a non-negative integer".to_string()
            })?),
        };
        Ok(ShellCommand::CollectionsRunRestrictionScan {
            collection_subject: CollectionSubjectRef::GeneSetResolution {
                report_id: choice.report_id.clone(),
            },
            member_bindings,
            enzymes,
            max_sites_per_enzyme,
            include_cut_geometry: self
                .gene_set_inspector
                .restriction_scan_form
                .include_cut_geometry,
            path: None,
        })
    }

    fn gene_set_tfbs_scan_shell_command(&self) -> Result<ShellCommand, String> {
        let choice = self
            .selected_gene_set_resolution()
            .ok_or_else(|| "Select a persisted gene-set resolution first".to_string())?;
        if choice.report.resolved_members.is_empty() {
            return Err("The selected gene set has no resolved members".to_string());
        }
        let mut member_bindings = Vec::with_capacity(choice.report.resolved_members.len());
        let mut missing_members = Vec::new();
        for member in &choice.report.resolved_members {
            let seq_id = self
                .gene_set_inspector
                .tfbs_scan_form
                .member_bindings
                .get(&member.dedup_key)
                .map(String::as_str)
                .unwrap_or_default()
                .trim();
            if seq_id.is_empty() {
                missing_members.push(member.symbol.clone());
            } else {
                member_bindings.push(TfbsHitScanCollectionMemberBinding {
                    stable_member_id: member.dedup_key.clone(),
                    seq_id: seq_id.to_string(),
                });
            }
        }
        if !missing_members.is_empty() {
            return Err(format!(
                "Bind one loaded DNA sequence for every member; missing: {}",
                missing_members.join(", ")
            ));
        }
        let form = &self.gene_set_inspector.tfbs_scan_form;
        let motifs = form
            .motifs
            .split([',', '\n'])
            .map(str::trim)
            .filter(|value| !value.is_empty())
            .map(str::to_string)
            .collect::<Vec<_>>();
        if motifs.is_empty() {
            return Err("Enter at least one TF or motif query".to_string());
        }
        let max_hits_per_member = parse_optional_usize(
            &form.max_hits_per_member,
            "Maximum retained hits per member",
        )?;
        if max_hits_per_member == Some(0) {
            return Err("Maximum retained hits per member must be at least 1".to_string());
        }
        Ok(ShellCommand::CollectionsRunTfbsScan {
            collection_subject: CollectionSubjectRef::GeneSetResolution {
                report_id: choice.report_id.clone(),
            },
            member_bindings,
            motifs,
            min_llr_bits: parse_optional_f64(&form.min_llr_bits, "Minimum LLR bits")?,
            min_llr_quantile: parse_optional_f64(&form.min_llr_quantile, "Minimum LLR quantile")?,
            per_tf_thresholds: parse_tfbs_threshold_overrides(
                &form.per_tf_min_llr_bits,
                &form.per_tf_min_llr_quantile,
            )?,
            max_hits_per_member,
            path: None,
        })
    }

    fn gene_set_digest_shell_command(&self) -> Result<ShellCommand, String> {
        let choice = self
            .selected_gene_set_resolution()
            .ok_or_else(|| "Select a persisted gene-set resolution first".to_string())?;
        if choice.report.resolved_members.is_empty() {
            return Err("The selected gene set has no resolved members".to_string());
        }
        let form = &self.gene_set_inspector.digest_form;
        let mut member_bindings = Vec::with_capacity(choice.report.resolved_members.len());
        let mut missing_members = Vec::new();
        for member in &choice.report.resolved_members {
            let seq_id = form
                .member_bindings
                .get(&member.dedup_key)
                .map(String::as_str)
                .unwrap_or_default()
                .trim();
            if seq_id.is_empty() {
                missing_members.push(member.symbol.clone());
            } else {
                member_bindings.push(DigestCollectionMemberBinding {
                    stable_member_id: member.dedup_key.clone(),
                    seq_id: seq_id.to_string(),
                });
            }
        }
        if !missing_members.is_empty() {
            return Err(format!(
                "Bind one loaded DNA sequence for every member; missing: {}",
                missing_members.join(", ")
            ));
        }
        let enzymes = form
            .enzymes
            .split([',', '\n'])
            .map(str::trim)
            .filter(|value| !value.is_empty())
            .map(str::to_string)
            .collect::<Vec<_>>();
        if enzymes.is_empty() {
            return Err("Enter at least one restriction enzyme".to_string());
        }
        let expected_plan_fingerprint_sha256 = if form.apply_change {
            Some(
                optional_trimmed(&form.expected_plan_fingerprint_sha256)
                    .ok_or_else(|| "Preview this digest before applying it".to_string())?,
            )
        } else {
            None
        };
        Ok(ShellCommand::CollectionsRunDigest {
            collection_subject: CollectionSubjectRef::GeneSetResolution {
                report_id: choice.report_id.clone(),
            },
            member_bindings,
            enzymes,
            output_prefix: optional_trimmed(&form.output_prefix),
            dry_run: !form.apply_change,
            expected_plan_fingerprint_sha256,
            path: None,
        })
    }

    fn gene_set_create_pool_shell_command(&self) -> Result<ShellCommand, String> {
        let choice = self
            .selected_gene_set_resolution()
            .ok_or_else(|| "Select a persisted gene-set resolution first".to_string())?;
        if choice.report.resolved_members.is_empty() {
            return Err("The selected gene set has no resolved members".to_string());
        }
        let unresolved_count = choice
            .report
            .unresolved_member_count
            .max(choice.report.unresolved_members.len());
        if unresolved_count > 0 {
            return Err(format!(
                "Resolve all members before creating a physical pool; {unresolved_count} remain unresolved"
            ));
        }
        let form = &self.gene_set_inspector.pool_form;
        let mut member_bindings = Vec::with_capacity(choice.report.resolved_members.len());
        let mut missing_members = Vec::new();
        let mut bound_source_containers = BTreeSet::new();
        for member in &choice.report.resolved_members {
            let source_container_id = form
                .member_bindings
                .get(&member.dedup_key)
                .map(String::as_str)
                .unwrap_or_default()
                .trim();
            if source_container_id.is_empty() {
                missing_members.push(member.symbol.clone());
            } else if !bound_source_containers.insert(source_container_id.to_string()) {
                return Err(format!(
                    "Source container '{source_container_id}' is selected more than once; bind one distinct physical tube per member"
                ));
            } else {
                member_bindings.push(GeneSetPoolMemberBinding {
                    stable_member_id: member.dedup_key.clone(),
                    source_container_id: source_container_id.to_string(),
                });
            }
        }
        if !missing_members.is_empty() {
            return Err(format!(
                "Bind one exclusive singleton source container for every member; missing: {}",
                missing_members.join(", ")
            ));
        }
        let expected_plan_fingerprint_sha256 = if form.apply_change {
            Some(
                optional_trimmed(&form.expected_plan_fingerprint_sha256)
                    .ok_or_else(|| "Preview this physical pool before creating it".to_string())?,
            )
        } else {
            None
        };
        Ok(ShellCommand::GeneSetsCreatePool {
            resolution_id: choice.report_id.clone(),
            member_bindings,
            output_prefix: optional_trimmed(&form.output_prefix),
            container_name: optional_trimmed(&form.container_name),
            dry_run: !form.apply_change,
            expected_plan_fingerprint_sha256,
            path: None,
        })
    }

    fn seed_gene_set_promoter_form_from_selection(&mut self, selection_changed: bool) {
        let selected_genome_id = self
            .selected_gene_set_resolution()
            .and_then(|choice| choice.report.genome_id.clone())
            .filter(|value| !value.trim().is_empty())
            .unwrap_or_else(|| self.genome_id.clone());
        let form = &mut self.gene_set_inspector.promoter_form;
        if selection_changed || form.genome_id.trim().is_empty() {
            form.genome_id = selected_genome_id;
        }
        if form.genome_catalog_path.trim().is_empty() {
            form.genome_catalog_path = self.genome_catalog_path.clone();
        }
        if form.cache_dir.trim().is_empty() {
            form.cache_dir = self.genome_cache_dir.clone();
        }
    }

    fn gene_set_promoter_shell_command(&self) -> Result<ShellCommand, String> {
        let choice = self
            .selected_gene_set_resolution()
            .ok_or_else(|| "Select a persisted gene-set resolution first".to_string())?;
        if choice.report.resolved_members.is_empty() {
            return Err("The selected gene set has no resolved members".to_string());
        }
        let form = &self.gene_set_inspector.promoter_form;
        let genome_id = optional_trimmed(&form.genome_id)
            .ok_or_else(|| "Enter the annotation genome for promoter derivation".to_string())?;
        let parse_span = |value: &str, label: &str| {
            value
                .trim()
                .parse::<usize>()
                .map_err(|_| format!("{label} must be a non-negative integer"))
        };
        Ok(ShellCommand::GeneSetsPromoterCohort {
            genome_id,
            source: None,
            resolution: Some(Box::new(choice.report.clone())),
            relationship: form.relationship,
            gene_group_catalog_path: optional_trimmed(&form.gene_group_catalog_path),
            genome_catalog_path: optional_trimmed(&form.genome_catalog_path),
            cache_dir: optional_trimmed(&form.cache_dir),
            upstream_bp: parse_span(&form.upstream_bp, "Upstream bp")?,
            downstream_bp: parse_span(&form.downstream_bp, "Downstream bp")?,
            allow_draft: false,
            allow_deprecated: false,
            output: optional_trimmed(&form.output_path),
        })
    }

    fn gene_set_collection_context_readiness(
        report: &GeneSetResolutionReport,
        requirement: CollectionContextRequirement,
        target_genome_id: &str,
    ) -> CollectionLauncherReadiness {
        match requirement {
            CollectionContextRequirement::ContextAgnostic => CollectionLauncherReadiness::Ready,
            CollectionContextRequirement::Homogeneous => {
                let members = report
                    .resolved_members
                    .iter()
                    .map(|member| CollectionMemberRef {
                        stable_member_id: member.dedup_key.clone(),
                        context_id: member.context_id.clone(),
                        ..CollectionMemberRef::default()
                    })
                    .collect::<Vec<_>>();
                match homogeneous_collection_biological_context(
                    &report.biological_contexts,
                    &members,
                )
                .and_then(|context| {
                    validate_collection_context_target_genome(&context, target_genome_id)
                }) {
                    Ok(()) => CollectionLauncherReadiness::Ready,
                    Err(error) => CollectionLauncherReadiness::PolicyRejected {
                        reason: error.reason,
                        detail: error.detail,
                    },
                }
            }
            CollectionContextRequirement::NotReviewed
            | CollectionContextRequirement::Partitionable
            | CollectionContextRequirement::ExplicitCrossContext => {
                CollectionLauncherReadiness::AdapterUnavailable {
                    detail: format!(
                        "The GUI launcher does not implement {:?} collection-context readiness",
                        requirement
                    ),
                }
            }
        }
    }

    fn gene_set_collection_operation_readiness(
        &self,
        row: &CollectionLauncherRow,
        report: &GeneSetResolutionReport,
    ) -> CollectionLauncherReadiness {
        let baseline = row.baseline_readiness();
        if !baseline.is_ready() {
            return baseline;
        }
        let Some(adapter) = row.adapter else {
            return baseline;
        };
        if adapter == CollectionLauncherAdapter::CreatePool {
            let unresolved_count = report
                .unresolved_member_count
                .max(report.unresolved_members.len());
            if unresolved_count > 0 {
                return CollectionLauncherReadiness::NeedsInput {
                    detail: format!(
                        "Resolve all gene-set members before creating a physical pool; {unresolved_count} remain unresolved"
                    ),
                };
            }
        }
        if matches!(
            adapter,
            CollectionLauncherAdapter::PrimerSpecificity
                | CollectionLauncherAdapter::RestrictionScan
                | CollectionLauncherAdapter::TfbsScan
                | CollectionLauncherAdapter::Digest
                | CollectionLauncherAdapter::CreatePool
        ) {
            let bindings = match adapter {
                CollectionLauncherAdapter::PrimerSpecificity => {
                    &self.gene_set_inspector.specificity_form.member_bindings
                }
                CollectionLauncherAdapter::RestrictionScan => {
                    &self
                        .gene_set_inspector
                        .restriction_scan_form
                        .member_bindings
                }
                CollectionLauncherAdapter::TfbsScan => {
                    &self.gene_set_inspector.tfbs_scan_form.member_bindings
                }
                CollectionLauncherAdapter::Digest => {
                    &self.gene_set_inspector.digest_form.member_bindings
                }
                CollectionLauncherAdapter::CreatePool => {
                    &self.gene_set_inspector.pool_form.member_bindings
                }
                CollectionLauncherAdapter::PromoterCohort => unreachable!(),
            };
            let missing = report
                .resolved_members
                .iter()
                .filter(|member| {
                    bindings
                        .get(&member.dedup_key)
                        .is_none_or(|value| value.trim().is_empty())
                })
                .map(|member| member.symbol.clone())
                .collect::<Vec<_>>();
            if !missing.is_empty() {
                return CollectionLauncherReadiness::NeedsBindings {
                    detail: format!(
                        "Bind one exact {} for: {}",
                        match adapter {
                            CollectionLauncherAdapter::PrimerSpecificity => "primer-design report",
                            CollectionLauncherAdapter::RestrictionScan => "loaded DNA sequence",
                            CollectionLauncherAdapter::TfbsScan => "loaded DNA sequence",
                            CollectionLauncherAdapter::Digest => "loaded DNA sequence",
                            CollectionLauncherAdapter::CreatePool => {
                                "exclusive singleton source container"
                            }
                            CollectionLauncherAdapter::PromoterCohort => unreachable!(),
                        },
                        missing.join(", ")
                    ),
                };
            }
        }
        let command = match adapter {
            CollectionLauncherAdapter::PrimerSpecificity => {
                self.gene_set_specificity_shell_command()
            }
            CollectionLauncherAdapter::RestrictionScan => {
                self.gene_set_restriction_scan_shell_command()
            }
            CollectionLauncherAdapter::TfbsScan => self.gene_set_tfbs_scan_shell_command(),
            CollectionLauncherAdapter::Digest => self.gene_set_digest_shell_command(),
            CollectionLauncherAdapter::CreatePool => self.gene_set_create_pool_shell_command(),
            CollectionLauncherAdapter::PromoterCohort => self.gene_set_promoter_shell_command(),
        };
        let command = match command {
            Ok(command) => command,
            Err(detail) => return CollectionLauncherReadiness::NeedsInput { detail },
        };
        let target_genome_id = match &command {
            ShellCommand::CollectionsRunPrimerSpecificity {
                target_genome_id, ..
            }
            | ShellCommand::GeneSetsPromoterCohort {
                genome_id: target_genome_id,
                ..
            } => target_genome_id,
            ShellCommand::CollectionsRunRestrictionScan { .. }
            | ShellCommand::CollectionsRunTfbsScan { .. }
            | ShellCommand::CollectionsRunDigest { .. }
            | ShellCommand::GeneSetsCreatePool { .. } => "",
            _ => unreachable!("collection launcher adapters construct known shell commands"),
        };
        let normalized_report = if adapter == CollectionLauncherAdapter::PromoterCohort {
            let mut candidate = report.clone();
            if let Err(error) = candidate.ensure_default_biological_context() {
                let reason = match error {
                    gentle_protocol::BiologicalContextResolutionError::ConflictingField {
                        ..
                    }
                    | gentle_protocol::BiologicalContextResolutionError::DuplicateContextId {
                        ..
                    } => CollectionLiftRejectionReason::MixedBiologicalContext,
                    gentle_protocol::BiologicalContextResolutionError::EmptyContextId
                    | gentle_protocol::BiologicalContextResolutionError::UnknownContextId {
                        ..
                    } => CollectionLiftRejectionReason::MissingBiologicalContext,
                };
                return CollectionLauncherReadiness::PolicyRejected {
                    reason,
                    detail: error.to_string(),
                };
            }
            Some(candidate)
        } else {
            None
        };
        Self::gene_set_collection_context_readiness(
            normalized_report.as_ref().unwrap_or(report),
            row.policy.context_requirement,
            target_genome_id,
        )
    }

    fn selected_gene_set_collection_shell_command(&self) -> Result<ShellCommand, String> {
        match self.gene_set_inspector.selected_operation {
            CollectionLauncherAdapter::PrimerSpecificity => {
                self.gene_set_specificity_shell_command()
            }
            CollectionLauncherAdapter::RestrictionScan => {
                self.gene_set_restriction_scan_shell_command()
            }
            CollectionLauncherAdapter::TfbsScan => self.gene_set_tfbs_scan_shell_command(),
            CollectionLauncherAdapter::Digest => self.gene_set_digest_shell_command(),
            CollectionLauncherAdapter::CreatePool => self.gene_set_create_pool_shell_command(),
            CollectionLauncherAdapter::PromoterCohort => self.gene_set_promoter_shell_command(),
        }
    }

    fn collection_operation_task_result_from_output(
        snapshot: &GentleEngine,
        adapter: CollectionLauncherAdapter,
        expected_subject: &CollectionSubjectRef,
        output: serde_json::Value,
    ) -> Result<CollectionOperationTaskResult, EngineError> {
        let invalid_output = |message: String| EngineError::new(ErrorCode::Internal, message);
        match adapter {
            CollectionLauncherAdapter::PrimerSpecificity => {
                let report = serde_json::from_value::<CollectionOperationReport>(
                    output
                        .get("report")
                        .cloned()
                        .unwrap_or(serde_json::Value::Null),
                )
                .map_err(|error| {
                    invalid_output(format!(
                        "Collection specificity command returned an invalid report: {error}"
                    ))
                })?;
                if &report.collection_subject != expected_subject {
                    return Err(invalid_output(format!(
                        "Collection specificity returned subject {:?}, expected {:?}",
                        report.collection_subject, expected_subject
                    )));
                }
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
                Ok(CollectionOperationTaskResult {
                    adapter,
                    report,
                    payload: CollectionOperationPayload::PrimerSpecificity { child_reports },
                })
            }
            CollectionLauncherAdapter::RestrictionScan => {
                let restriction_scan = serde_json::from_value::<
                    CollectionRestrictionSiteScanReport,
                >(
                    output
                        .get("report")
                        .cloned()
                        .unwrap_or(serde_json::Value::Null),
                )
                .map_err(|error| {
                    invalid_output(format!(
                        "Collection restriction-scan command returned an invalid report: {error}"
                    ))
                })?;
                if &restriction_scan.collection_operation.collection_subject != expected_subject {
                    return Err(invalid_output(format!(
                        "Collection restriction scan returned subject {:?}, expected {:?}",
                        restriction_scan.collection_operation.collection_subject, expected_subject
                    )));
                }
                Ok(CollectionOperationTaskResult {
                    adapter,
                    report: restriction_scan.collection_operation.clone(),
                    payload: CollectionOperationPayload::RestrictionScan(Box::new(
                        restriction_scan,
                    )),
                })
            }
            CollectionLauncherAdapter::TfbsScan => {
                let tfbs_scan = serde_json::from_value::<CollectionTfbsHitScanReport>(
                    output
                        .get("report")
                        .cloned()
                        .unwrap_or(serde_json::Value::Null),
                )
                .map_err(|error| {
                    invalid_output(format!(
                        "Collection TFBS-scan command returned an invalid report: {error}"
                    ))
                })?;
                if &tfbs_scan.collection_operation.collection_subject != expected_subject {
                    return Err(invalid_output(format!(
                        "Collection TFBS scan returned subject {:?}, expected {:?}",
                        tfbs_scan.collection_operation.collection_subject, expected_subject
                    )));
                }
                Ok(CollectionOperationTaskResult {
                    adapter,
                    report: tfbs_scan.collection_operation.clone(),
                    payload: CollectionOperationPayload::TfbsScan(Box::new(tfbs_scan)),
                })
            }
            CollectionLauncherAdapter::Digest => {
                let digest = serde_json::from_value::<CollectionDigestReport>(
                    output
                        .get("report")
                        .cloned()
                        .unwrap_or(serde_json::Value::Null),
                )
                .map_err(|error| {
                    invalid_output(format!(
                        "Collection digest command returned an invalid report: {error}"
                    ))
                })?;
                if &digest.collection_operation.collection_subject != expected_subject {
                    return Err(invalid_output(format!(
                        "Collection digest returned subject {:?}, expected {:?}",
                        digest.collection_operation.collection_subject, expected_subject
                    )));
                }
                Ok(CollectionOperationTaskResult {
                    adapter,
                    report: digest.collection_operation.clone(),
                    payload: CollectionOperationPayload::Digest(Box::new(digest)),
                })
            }
            CollectionLauncherAdapter::CreatePool => {
                let pool = serde_json::from_value::<GeneSetPoolCreationReport>(
                    output
                        .get("report")
                        .cloned()
                        .unwrap_or(serde_json::Value::Null),
                )
                .map_err(|error| {
                    invalid_output(format!(
                        "Gene-set pool command returned an invalid report: {error}"
                    ))
                })?;
                if &pool.collection_operation.collection_subject != expected_subject {
                    return Err(invalid_output(format!(
                        "Gene-set pool returned subject {:?}, expected {:?}",
                        pool.collection_operation.collection_subject, expected_subject
                    )));
                }
                Ok(CollectionOperationTaskResult {
                    adapter,
                    report: pool.collection_operation.clone(),
                    payload: CollectionOperationPayload::CreatePool(Box::new(pool)),
                })
            }
            CollectionLauncherAdapter::PromoterCohort => {
                let promoter = serde_json::from_value::<GeneSetPromoterCohortReport>(output)
                    .map_err(|error| {
                        invalid_output(format!(
                            "Promoter-cohort command returned an invalid report: {error}"
                        ))
                    })?;
                let report = promoter
                    .collection_operation
                    .as_deref()
                    .cloned()
                    .ok_or_else(|| {
                        invalid_output(
                            "Promoter-cohort report omitted its collection-operation wrapper"
                                .to_string(),
                        )
                    })?;
                if &report.collection_subject != expected_subject {
                    return Err(invalid_output(format!(
                        "Promoter cohort returned subject {:?}, expected {:?}",
                        report.collection_subject, expected_subject
                    )));
                }
                Ok(CollectionOperationTaskResult {
                    adapter,
                    report,
                    payload: CollectionOperationPayload::PromoterCohort(Box::new(promoter)),
                })
            }
        }
    }

    fn start_gene_set_collection_operation_task(&mut self) {
        if self.gene_set_inspector.task.is_some() {
            self.gene_set_inspector.status =
                "A gene-set collection operation is already running".to_string();
            return;
        }
        let command = match self.selected_gene_set_collection_shell_command() {
            Ok(command) => command,
            Err(error) => {
                self.gene_set_inspector.status = error;
                return;
            }
        };
        let adapter = self.gene_set_inspector.selected_operation;
        let resolution_id = self.gene_set_inspector.selected_resolution_id.clone();
        let expected_subject = CollectionSubjectRef::GeneSetResolution {
            report_id: resolution_id.clone(),
        };
        let member_count = self
            .selected_gene_set_resolution()
            .map(|choice| choice.report.resolved_members.len())
            .unwrap_or_default();
        let job_id = self.alloc_background_job_id();
        let cancel_requested = Arc::new(AtomicBool::new(false));
        let worker_cancel = cancel_requested.clone();
        let engine = self.engine.clone();
        let (tx, receiver) = mpsc::channel::<Result<CollectionOperationTaskResult, EngineError>>();
        let runtime_frame = Self::push_runtime_background_job_frame(
            BackgroundJobKind::CollectionOperation,
            job_id,
            format!(
                "operation={}, resolution='{resolution_id}', members={member_count}",
                adapter.capability_name()
            ),
        );
        runtime_frame.update_phase("running");
        self.push_job_event(
            BackgroundJobKind::CollectionOperation,
            BackgroundJobEventPhase::Started,
            Some(job_id),
            format!(
                "Collection operation {} started: resolution='{}', members={}",
                adapter.capability_name(),
                resolution_id,
                member_count
            ),
        );
        std::thread::spawn(move || {
            let result = crate::background_engine::execute_on_engine_snapshot(
                &engine,
                move |snapshot| {
                    if worker_cancel.load(Ordering::Relaxed) {
                        return Err(EngineError::invalid_input(
                            "Collection operation was discarded before execution",
                        ));
                    }
                    let shell_result = execute_shell_command(snapshot, &command)
                        .map_err(|message| EngineError::new(ErrorCode::InvalidInput, message))?;
                    if worker_cancel.load(Ordering::Relaxed) {
                        return Err(EngineError::invalid_input(
                            "Collection operation finished after discard was requested; detached result was not committed",
                        ));
                    }
                    let task_result = Self::collection_operation_task_result_from_output(
                        snapshot,
                        adapter,
                        &expected_subject,
                        shell_result.output,
                    )?;
                    if worker_cancel.load(Ordering::Relaxed) {
                        return Err(EngineError::invalid_input(
                            "Collection operation finished after discard was requested; detached result was not committed",
                        ));
                    }
                    Ok(task_result)
                },
            );
            let _ = tx.send(result);
        });
        self.gene_set_inspector.status = format!(
            "Running {} for {member_count} gene-set member(s) in background...",
            adapter.label()
        );
        self.gene_set_inspector.task = Some(GeneSetCollectionOperationTask {
            job_id,
            adapter,
            resolution_id,
            started: Instant::now(),
            cancel_requested,
            runtime_frame,
            receiver,
        });
    }

    fn request_gene_set_collection_operation_discard(&mut self, origin: &str) {
        let Some(task) = self.gene_set_inspector.task.as_mut() else {
            self.gene_set_inspector.status =
                "No running gene-set collection operation to discard".to_string();
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
            "Discard requested. Active external work may finish before GENtle drops the detached result; output files already written by the shared command are not rolled back."
                .to_string();
        self.push_job_event(
            BackgroundJobKind::CollectionOperation,
            BackgroundJobEventPhase::CancelRequested,
            Some(job_id),
            format!("Detached result discard requested from {origin}"),
        );
    }

    pub(super) fn poll_gene_set_collection_operation_task(&mut self, ctx: &egui::Context) {
        let Some(task) = self.gene_set_inspector.task.as_ref() else {
            return;
        };
        ctx.request_repaint_after(std::time::Duration::from_millis(100));
        let outcome = match task.receiver.try_recv() {
            Ok(result) => Some(result),
            Err(mpsc::TryRecvError::Empty) => None,
            Err(mpsc::TryRecvError::Disconnected) => Some(Err(EngineError::new(
                ErrorCode::Io,
                "Gene-set collection-operation background worker disconnected",
            ))),
        };
        let Some(outcome) = outcome else {
            return;
        };
        let Some(task) = self.gene_set_inspector.task.take() else {
            return;
        };
        let elapsed = task.started.elapsed().as_secs_f64();
        if task.cancel_requested.load(Ordering::Relaxed) {
            task.runtime_frame
                .cancel("Collection operation result discarded".to_string());
            self.gene_set_inspector.status = format!(
                "Collection operation result discarded after {elapsed:.1}s; no result was displayed. Output files already written by the shared command are not rolled back."
            );
            self.push_job_event(
                BackgroundJobKind::CollectionOperation,
                BackgroundJobEventPhase::IgnoredStale,
                Some(task.job_id),
                self.gene_set_inspector.status.clone(),
            );
            ctx.request_repaint();
            return;
        }
        match outcome {
            Ok(result) => {
                task.runtime_frame.update_phase("completed");
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
                self.gene_set_inspector.status = match &result.payload {
                    CollectionOperationPayload::PrimerSpecificity { .. } => format!(
                        "Collection specificity completed in {elapsed:.1}s: {succeeded} executed, {failed} failed. Biological pass/fail is shown per child report."
                    ),
                    CollectionOperationPayload::RestrictionScan(report) => format!(
                        "Restriction scan completed in {elapsed:.1}s: {succeeded} executed, {failed} failed, {} total matched site(s).",
                        report.total_matched_site_count
                    ),
                    CollectionOperationPayload::TfbsScan(report) => format!(
                        "TFBS hit scan completed in {elapsed:.1}s: {succeeded} executed, {failed} failed, {} retained hit(s), aggregate counts {}.",
                        report.total_retained_hit_count,
                        if report.aggregate_counts_complete {
                            "complete"
                        } else {
                            "incomplete"
                        }
                    ),
                    CollectionOperationPayload::Digest(report) => format!(
                        "Collection digest {} in {elapsed:.1}s: {succeeded} executed, {failed} failed, {} fragment(s) {}.",
                        if report.collection_operation.applied {
                            "applied"
                        } else {
                            "previewed"
                        },
                        report.total_planned_fragment_count,
                        if report.collection_operation.applied {
                            "created"
                        } else {
                            "planned"
                        }
                    ),
                    CollectionOperationPayload::CreatePool(report) => format!(
                        "Physical gene-set pool {} in {elapsed:.1}s: {} member aliquot(s) {}, source containers retained.",
                        if report.collection_operation.applied {
                            "created"
                        } else {
                            "previewed"
                        },
                        report.planned_member_count,
                        if report.collection_operation.applied {
                            "materialized"
                        } else {
                            "planned"
                        }
                    ),
                    CollectionOperationPayload::PromoterCohort(report) => format!(
                        "Promoter cohort completed in {elapsed:.1}s: {} window(s), {} unresolved. Derivation stores the cohort and may normalize missing source-context metadata.",
                        report.returned_window_count,
                        report.unresolved_members.len()
                    ),
                };
                if let CollectionOperationPayload::Digest(report) = &result.payload {
                    self.gene_set_inspector
                        .digest_form
                        .expected_plan_fingerprint_sha256 = report.plan_fingerprint_sha256.clone();
                    self.gene_set_inspector.digest_form.apply_change = false;
                }
                if let CollectionOperationPayload::CreatePool(report) = &result.payload {
                    self.gene_set_inspector
                        .pool_form
                        .expected_plan_fingerprint_sha256 = report.plan_fingerprint_sha256.clone();
                    self.gene_set_inspector.pool_form.apply_change = false;
                }
                self.refresh_gene_set_inspector_catalog(true);
                let selection_changed =
                    self.gene_set_inspector.selected_resolution_id != task.resolution_id;
                self.gene_set_inspector.selected_resolution_id = task.resolution_id.clone();
                self.reconcile_gene_set_inspector_bindings();
                if selection_changed {
                    self.seed_gene_set_promoter_form_from_selection(true);
                }
                self.gene_set_inspector.last_result = Some(result);
                self.push_job_event(
                    BackgroundJobKind::CollectionOperation,
                    BackgroundJobEventPhase::Completed,
                    Some(task.job_id),
                    format!(
                        "Collection operation {} completed in {elapsed:.1}s for '{}': {succeeded} executed, {failed} failed",
                        task.adapter.capability_name(), task.resolution_id
                    ),
                );
            }
            Err(error) => {
                let stale = error.message.contains("became stale");
                let discarded = error.message.contains("discard");
                self.gene_set_inspector.status = if stale {
                    format!(
                        "Collection operation produced a stale result after {elapsed:.1}s; project state was not changed."
                    )
                } else if discarded {
                    format!(
                        "Collection operation result discarded after {elapsed:.1}s; project state was not changed."
                    )
                } else {
                    format!(
                        "Collection operation {} failed after {elapsed:.1}s: {}",
                        task.adapter.label(),
                        error.message
                    )
                };
                if stale || discarded {
                    task.runtime_frame.cancel(error.message.clone());
                } else {
                    task.runtime_frame.fail(error.message.clone());
                }
                self.push_job_event(
                    BackgroundJobKind::CollectionOperation,
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

    pub(super) fn has_active_gene_set_collection_operation_task(&self) -> bool {
        self.gene_set_inspector.task.is_some()
    }

    pub(super) fn render_gene_set_collection_operation_background_job(&mut self, ui: &mut Ui) {
        ui.separator();
        ui.strong("Gene-set Collection Operation");
        let mut discard_clicked = false;
        if let Some(task) = self.gene_set_inspector.task.as_ref() {
            ui.horizontal_wrapped(|ui| {
                ui.add(egui::Spinner::new());
                ui.label(format!(
                    "Running #{} {} for '{}' ({:.1}s)",
                    task.job_id,
                    task.adapter.label(),
                    task.resolution_id,
                    task.started.elapsed().as_secs_f32()
                ));
                if task.cancel_requested.load(Ordering::Relaxed) {
                    ui.small("Discard requested...");
                } else if ui
                    .button("Discard result")
                    .on_hover_text(
                        "Drop the detached collection result. Active external work may finish, and output files already written by the shared command are not rolled back.",
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
                        "Open the gene-set inspector to configure a collection operation",
                    )
                    .clicked()
                {
                    self.open_gene_set_inspector_dialog();
                }
            });
        }
        if discard_clicked {
            self.request_gene_set_collection_operation_discard("background jobs panel");
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

    fn render_gene_set_specificity_form(&mut self, ui: &mut Ui, choice: &GeneSetResolutionChoice) {
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
                for report_id in self
                    .gene_set_inspector
                    .specificity_form
                    .member_bindings
                    .values_mut()
                {
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
                                .specificity_form
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
        ui.strong("Specificity parameters");
        let form = &mut self.gene_set_inspector.specificity_form;
        ui.horizontal_wrapped(|ui| {
            ui.label("Pair rank");
            ui.add(
                egui::TextEdit::singleline(&mut form.pair_rank_1based).desired_width(70.0),
            )
            .on_hover_text("One-based accepted pair rank, applied to every bound report");
            ui.label("Prepared target genome/transcriptome");
            ui.add(
                egui::TextEdit::singleline(&mut form.target_genome_id)
                    .desired_width(280.0)
                    .hint_text("Genome catalog ID"),
            )
            .on_hover_text(
                "Prepared local BLAST target id. GENtle uses the same specificity policy and database checks as the single-report command.",
            );
        });
        egui::CollapsingHeader::new("Specificity paths")
            .default_open(false)
            .show(ui, |ui| {
                ui.horizontal_wrapped(|ui| {
                    ui.label("Genome catalog");
                    ui.add(
                        egui::TextEdit::singleline(&mut form.catalog_path)
                            .desired_width(420.0)
                            .hint_text("empty = configured/default catalog"),
                    );
                });
                ui.horizontal_wrapped(|ui| {
                    ui.label("Genome cache");
                    ui.add(
                        egui::TextEdit::singleline(&mut form.cache_dir)
                            .desired_width(420.0)
                            .hint_text("empty = configured/default cache"),
                    );
                });
            });
    }

    fn render_gene_set_restriction_scan_form(
        &mut self,
        ui: &mut Ui,
        choice: &GeneSetResolutionChoice,
    ) {
        let sequence_ids = self
            .engine
            .read()
            .ok()
            .map(|engine| engine.state().sequences.keys().cloned().collect::<Vec<_>>())
            .unwrap_or_default();
        ui.strong("Sequence bindings");
        ui.small("Each resolved gene is bound explicitly to one loaded DNA sequence.");
        ui.horizontal_wrapped(|ui| {
            ui.label(format!("Loaded DNA sequences: {}", sequence_ids.len()));
            if ui
                .button("Clear bindings")
                .on_hover_text("Clear every member-to-sequence selection")
                .clicked()
            {
                for seq_id in self
                    .gene_set_inspector
                    .restriction_scan_form
                    .member_bindings
                    .values_mut()
                {
                    seq_id.clear();
                }
            }
        });
        egui::ScrollArea::vertical()
            .id_salt("gene_set_restriction_scan_bindings")
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
                                .restriction_scan_form
                                .member_bindings
                                .entry(member.dedup_key.clone())
                                .or_default();
                            egui::ComboBox::from_id_salt((
                                "gene_set_member_restriction_sequence",
                                &member.dedup_key,
                            ))
                            .selected_text(if binding.is_empty() {
                                "Select loaded sequence..."
                            } else {
                                binding.as_str()
                            })
                            .width(420.0)
                            .show_ui(ui, |ui| {
                                ui.selectable_value(binding, String::new(), "Not selected");
                                for seq_id in &sequence_ids {
                                    ui.selectable_value(binding, seq_id.clone(), seq_id);
                                }
                            });
                        });
                    }
                },
            );

        ui.separator();
        ui.strong("Restriction scan parameters");
        let form = &mut self.gene_set_inspector.restriction_scan_form;
        ui.horizontal_wrapped(|ui| {
            ui.label("Enzymes");
            ui.add(
                egui::TextEdit::singleline(&mut form.enzymes)
                    .desired_width(300.0)
                    .hint_text("Configured preferred enzymes"),
            )
            .on_hover_text("Comma-separated names. Leave empty to use configured preferred enzymes, then the built-in default set.");
            ui.label("Maximum sites / enzyme");
            ui.add(
                egui::TextEdit::singleline(&mut form.max_sites_per_enzyme)
                    .desired_width(100.0)
                    .hint_text("No limit"),
            );
            ui.checkbox(&mut form.include_cut_geometry, "Cut geometry");
        });
    }

    fn render_gene_set_tfbs_scan_form(&mut self, ui: &mut Ui, choice: &GeneSetResolutionChoice) {
        let sequence_ids = self
            .engine
            .read()
            .ok()
            .map(|engine| engine.state().sequences.keys().cloned().collect::<Vec<_>>())
            .unwrap_or_default();
        ui.strong("Sequence bindings");
        ui.small("Each resolved gene is bound explicitly to one loaded DNA sequence.");
        ui.horizontal_wrapped(|ui| {
            ui.label(format!("Loaded DNA sequences: {}", sequence_ids.len()));
            if ui
                .button("Clear bindings")
                .on_hover_text("Clear every member-to-sequence selection")
                .clicked()
            {
                for seq_id in self
                    .gene_set_inspector
                    .tfbs_scan_form
                    .member_bindings
                    .values_mut()
                {
                    seq_id.clear();
                }
            }
        });
        egui::ScrollArea::vertical()
            .id_salt("gene_set_tfbs_scan_bindings")
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
                                .tfbs_scan_form
                                .member_bindings
                                .entry(member.dedup_key.clone())
                                .or_default();
                            egui::ComboBox::from_id_salt((
                                "gene_set_member_tfbs_sequence",
                                &member.dedup_key,
                            ))
                            .selected_text(if binding.is_empty() {
                                "Select loaded sequence..."
                            } else {
                                binding.as_str()
                            })
                            .width(420.0)
                            .show_ui(ui, |ui| {
                                ui.selectable_value(binding, String::new(), "Not selected");
                                for seq_id in &sequence_ids {
                                    ui.selectable_value(binding, seq_id.clone(), seq_id);
                                }
                            });
                        });
                    }
                },
            );

        ui.separator();
        ui.strong("TFBS hit-scan parameters");
        let form = &mut self.gene_set_inspector.tfbs_scan_form;
        ui.horizontal_wrapped(|ui| {
            ui.label("Motifs");
            ui.add(
                egui::TextEdit::singleline(&mut form.motifs)
                    .desired_width(360.0)
                    .hint_text("SP1, CTCF, MA0139.1"),
            )
            .on_hover_text(
                "Comma-separated motif IDs, TF names, groups, family prefixes, or IUPAC motifs",
            );
            ui.label("Max hits / member");
            ui.add(
                egui::TextEdit::singleline(&mut form.max_hits_per_member)
                    .desired_width(100.0)
                    .hint_text("Unlimited"),
            )
            .on_hover_text(
                "A cap can stop motif iteration; leave empty for complete aggregate counts",
            );
        });
        ui.horizontal_wrapped(|ui| {
            ui.label("Minimum LLR bits");
            ui.add(
                egui::TextEdit::singleline(&mut form.min_llr_bits)
                    .desired_width(100.0)
                    .hint_text("No minimum"),
            );
            ui.label("Minimum LLR quantile");
            ui.add(
                egui::TextEdit::singleline(&mut form.min_llr_quantile)
                    .desired_width(100.0)
                    .hint_text("0.0"),
            );
        });
        egui::CollapsingHeader::new("Per-TF thresholds")
            .default_open(false)
            .show(ui, |ui| {
                ui.horizontal_wrapped(|ui| {
                    ui.label("Minimum bits");
                    ui.add(
                        egui::TextEdit::singleline(&mut form.per_tf_min_llr_bits)
                            .desired_width(360.0)
                            .hint_text("SP1=4.0, CTCF=5.0"),
                    );
                });
                ui.horizontal_wrapped(|ui| {
                    ui.label("Minimum quantiles");
                    ui.add(
                        egui::TextEdit::singleline(&mut form.per_tf_min_llr_quantile)
                            .desired_width(360.0)
                            .hint_text("SP1=0.95, CTCF=0.99"),
                    );
                });
            });
    }

    fn render_gene_set_digest_form(&mut self, ui: &mut Ui, choice: &GeneSetResolutionChoice) {
        let sequence_ids = self
            .engine
            .read()
            .ok()
            .map(|engine| engine.state().sequences.keys().cloned().collect::<Vec<_>>())
            .unwrap_or_default();
        let mut plan_inputs_changed = false;
        ui.strong("Sequence bindings");
        ui.small("Each resolved gene is bound explicitly to one loaded DNA sequence.");
        ui.horizontal_wrapped(|ui| {
            ui.label(format!("Loaded DNA sequences: {}", sequence_ids.len()));
            if ui
                .button("Clear bindings")
                .on_hover_text("Clear every member-to-sequence selection")
                .clicked()
            {
                for seq_id in self
                    .gene_set_inspector
                    .digest_form
                    .member_bindings
                    .values_mut()
                {
                    seq_id.clear();
                }
                plan_inputs_changed = true;
            }
        });
        egui::ScrollArea::vertical()
            .id_salt("gene_set_digest_bindings")
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
                                .digest_form
                                .member_bindings
                                .entry(member.dedup_key.clone())
                                .or_default();
                            egui::ComboBox::from_id_salt((
                                "gene_set_member_digest_sequence",
                                &member.dedup_key,
                            ))
                            .selected_text(if binding.is_empty() {
                                "Select loaded sequence..."
                            } else {
                                binding.as_str()
                            })
                            .width(420.0)
                            .show_ui(ui, |ui| {
                                plan_inputs_changed |= ui
                                    .selectable_value(binding, String::new(), "Not selected")
                                    .changed();
                                for seq_id in &sequence_ids {
                                    plan_inputs_changed |= ui
                                        .selectable_value(binding, seq_id.clone(), seq_id)
                                        .changed();
                                }
                            });
                        });
                    }
                },
            );

        ui.separator();
        ui.strong("Restriction digest parameters");
        let form = &mut self.gene_set_inspector.digest_form;
        ui.horizontal_wrapped(|ui| {
            ui.label("Enzymes");
            plan_inputs_changed |= ui
                .add(
                    egui::TextEdit::singleline(&mut form.enzymes)
                        .desired_width(300.0)
                        .hint_text("EcoRI, BamHI"),
                )
                .changed();
            ui.label("Output prefix");
            plan_inputs_changed |= ui
                .add(
                    egui::TextEdit::singleline(&mut form.output_prefix)
                        .desired_width(220.0)
                        .hint_text("Generated per member"),
                )
                .changed();
        });
        if plan_inputs_changed {
            form.expected_plan_fingerprint_sha256.clear();
            form.apply_change = false;
        }
        ui.horizontal_wrapped(|ui| {
            let has_preview = !form.expected_plan_fingerprint_sha256.trim().is_empty();
            ui.add_enabled_ui(has_preview, |ui| {
                ui.checkbox(&mut form.apply_change, "Apply previewed digest")
                    .on_hover_text("Materialize the exact fingerprint-locked fragment plan");
            });
            if has_preview {
                let fingerprint = form.expected_plan_fingerprint_sha256.as_str();
                let compact = fingerprint.chars().take(24).collect::<String>();
                ui.monospace(format!("plan {compact}..."))
                    .on_hover_text(fingerprint);
            } else {
                ui.small("No locked preview");
            }
        });
    }

    fn render_gene_set_create_pool_form(&mut self, ui: &mut Ui, choice: &GeneSetResolutionChoice) {
        let container_choices = self.gene_set_pool_source_container_choices();
        let mut plan_inputs_changed = false;
        ui.strong("Physical source tubes");
        ui.small(
            "Bind each resolved member to one exclusive singleton container. Pooling creates derived aliquots and retains every source container.",
        );
        ui.horizontal_wrapped(|ui| {
            ui.label(format!(
                "Eligible singleton containers: {}",
                container_choices.len()
            ));
            if ui
                .button("Clear bindings")
                .on_hover_text("Clear every member-to-source-container selection")
                .clicked()
            {
                for container_id in self
                    .gene_set_inspector
                    .pool_form
                    .member_bindings
                    .values_mut()
                {
                    container_id.clear();
                }
                plan_inputs_changed = true;
            }
        });
        egui::ScrollArea::vertical()
            .id_salt("gene_set_pool_bindings")
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
                                .pool_form
                                .member_bindings
                                .entry(member.dedup_key.clone())
                                .or_default();
                            let selected_label = container_choices
                                .iter()
                                .find(|(container_id, _)| container_id == binding)
                                .map(|(_, label)| label.as_str())
                                .unwrap_or("Select source tube...");
                            egui::ComboBox::from_id_salt((
                                "gene_set_member_pool_container",
                                &member.dedup_key,
                            ))
                            .selected_text(selected_label)
                            .width(460.0)
                            .show_ui(ui, |ui| {
                                plan_inputs_changed |= ui
                                    .selectable_value(binding, String::new(), "Not selected")
                                    .changed();
                                for (container_id, label) in &container_choices {
                                    plan_inputs_changed |= ui
                                        .selectable_value(binding, container_id.clone(), label)
                                        .changed();
                                }
                            });
                        });
                    }
                },
            );

        ui.separator();
        ui.strong("Pool output");
        let form = &mut self.gene_set_inspector.pool_form;
        ui.horizontal_wrapped(|ui| {
            ui.label("Aliquot ID prefix");
            plan_inputs_changed |= ui
                .add(
                    egui::TextEdit::singleline(&mut form.output_prefix)
                        .desired_width(220.0)
                        .hint_text("gene_set_pool"),
                )
                .changed();
            ui.label("Container name");
            plan_inputs_changed |= ui
                .add(
                    egui::TextEdit::singleline(&mut form.container_name)
                        .desired_width(260.0)
                        .hint_text("Generated from resolution ID"),
                )
                .changed();
        });
        if plan_inputs_changed {
            form.expected_plan_fingerprint_sha256.clear();
            form.apply_change = false;
        }
        ui.horizontal_wrapped(|ui| {
            let has_preview = !form.expected_plan_fingerprint_sha256.trim().is_empty();
            ui.add_enabled_ui(has_preview, |ui| {
                ui.checkbox(&mut form.apply_change, "Create previewed physical pool")
                    .on_hover_text("Apply the exact fingerprint-locked aliquot and container plan");
            });
            if has_preview {
                let fingerprint = form.expected_plan_fingerprint_sha256.as_str();
                let compact = fingerprint.chars().take(24).collect::<String>();
                ui.monospace(format!("plan {compact}..."))
                    .on_hover_text(fingerprint);
            } else {
                ui.small("No locked preview");
            }
        });
    }

    fn render_gene_set_promoter_form(&mut self, ui: &mut Ui) {
        ui.strong("Promoter-cohort parameters");
        ui.small(
            "Derive strand-aware promoter windows through the shared promoter/TSS resolver. Relationship expectations are non-blocking evidence framing, not proof of co-regulation.",
        );
        let form = &mut self.gene_set_inspector.promoter_form;
        ui.horizontal_wrapped(|ui| {
            ui.label("Annotation genome");
            ui.add(
                egui::TextEdit::singleline(&mut form.genome_id)
                    .desired_width(260.0)
                    .hint_text("Genome catalog ID"),
            );
            ui.label("Relationship");
            egui::ComboBox::from_id_salt("gene_set_promoter_relationship")
                .selected_text(gene_set_relationship_label(form.relationship))
                .show_ui(ui, |ui| {
                    for relationship in [
                        GeneSetCohortRelationship::Unspecified,
                        GeneSetCohortRelationship::Manual,
                        GeneSetCohortRelationship::CoRegulated,
                        GeneSetCohortRelationship::AntiCoRegulated,
                    ] {
                        ui.selectable_value(
                            &mut form.relationship,
                            relationship,
                            gene_set_relationship_label(relationship),
                        );
                    }
                });
        });
        ui.horizontal_wrapped(|ui| {
            ui.label("Upstream bp");
            ui.add(egui::TextEdit::singleline(&mut form.upstream_bp).desired_width(90.0));
            ui.label("Downstream bp");
            ui.add(egui::TextEdit::singleline(&mut form.downstream_bp).desired_width(90.0));
        });
        egui::CollapsingHeader::new("Promoter catalogs and output")
            .default_open(false)
            .show(ui, |ui| {
                for (label, value, hint) in [
                    (
                        "Gene-group catalog",
                        &mut form.gene_group_catalog_path,
                        "empty = configured/default catalog",
                    ),
                    (
                        "Genome catalog",
                        &mut form.genome_catalog_path,
                        "empty = configured/default catalog",
                    ),
                    (
                        "Genome cache",
                        &mut form.cache_dir,
                        "empty = configured/default cache",
                    ),
                    (
                        "Output JSON",
                        &mut form.output_path,
                        "optional portable cohort path",
                    ),
                ] {
                    ui.horizontal_wrapped(|ui| {
                        ui.label(label);
                        ui.add(
                            egui::TextEdit::singleline(value)
                                .desired_width(440.0)
                                .hint_text(hint),
                        );
                    });
                }
            });
    }

    fn render_gene_set_collection_result(&mut self, ui: &mut Ui) {
        let Some(result) = self.gene_set_inspector.last_result.clone() else {
            return;
        };
        let report = &result.report;
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
            "Latest result: {} | {} executed, {} failed",
            result.adapter.label(),
            succeeded,
            failed
        ));
        ui.small(format!(
            "report_id={} | lift={:?} | fingerprint={}",
            report.report_id, report.lifting_mode, report.collection_membership_fingerprint_sha256
        ));
        if ui
            .button("Copy collection report JSON")
            .on_hover_text("Copy the portable collection-operation report")
            .clicked()
            && let Ok(json) = serde_json::to_string_pretty(report)
        {
            ui.ctx().copy_text(json);
        }
        let child_reports = match &result.payload {
            CollectionOperationPayload::PrimerSpecificity { child_reports } => Some(child_reports),
            CollectionOperationPayload::RestrictionScan(_) => None,
            CollectionOperationPayload::TfbsScan(_) => None,
            CollectionOperationPayload::Digest(_) => None,
            CollectionOperationPayload::CreatePool(_) => None,
            CollectionOperationPayload::PromoterCohort(_) => None,
        };
        let mut open_child: Option<(String, String)> = None;
        let mut open_digest_sequence: Option<String> = None;
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
                            if let Some(parent) = row.member.parent_member_id.as_deref() {
                                ui.small(format!("derived from {parent}"));
                            }
                            if let Some(error) = &row.error {
                                ui.colored_label(
                                    ui.visuals().error_fg_color,
                                    error.message.clone(),
                                );
                            }
                            for report_id in &row.produced_report_ids {
                                if let Some(child) = child_reports
                                    .and_then(|children| children.get(report_id))
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
                    ui.small(format!("- {warning}"));
                }
            });
        }

        if let CollectionOperationPayload::RestrictionScan(restriction_scan) = &result.payload {
            ui.separator();
            ui.strong(format!(
                "Restriction sites: {} across {} successful sequence(s)",
                restriction_scan.total_matched_site_count,
                restriction_scan.member_reports.len()
            ));
            ui.small(format!(
                "enzymes={} | cut geometry={}",
                restriction_scan.effective_enzymes.join(", "),
                restriction_scan.include_cut_geometry
            ));
            if ui
                .button("Copy restriction scan JSON")
                .on_hover_text("Copy the complete portable collection restriction-scan report")
                .clicked()
                && let Ok(json) = serde_json::to_string_pretty(restriction_scan.as_ref())
            {
                ui.ctx().copy_text(json);
            }
        }

        if let CollectionOperationPayload::TfbsScan(tfbs_scan) = &result.payload {
            ui.separator();
            ui.strong(format!(
                "TFBS hits retained: {} across {} successful sequence(s)",
                tfbs_scan.total_retained_hit_count,
                tfbs_scan.member_reports.len()
            ));
            let completeness = if tfbs_scan.aggregate_counts_complete {
                "complete"
            } else {
                "incomplete"
            };
            ui.small(format!(
                "motifs={} | aggregate counts={} | truncated members={}",
                tfbs_scan.effective_motif_ids.join(", "),
                completeness,
                tfbs_scan.truncated_member_count
            ));
            if !tfbs_scan.retained_hit_counts_by_tf_id.is_empty() {
                ui.small(format!(
                    "retained by motif: {}",
                    tfbs_scan
                        .retained_hit_counts_by_tf_id
                        .iter()
                        .map(|(tf_id, count)| format!("{tf_id}={count}"))
                        .collect::<Vec<_>>()
                        .join(", ")
                ));
            }
            if ui
                .button("Copy TFBS scan JSON")
                .on_hover_text("Copy the complete portable collection TFBS hit-scan report")
                .clicked()
                && let Ok(json) = serde_json::to_string_pretty(tfbs_scan.as_ref())
            {
                ui.ctx().copy_text(json);
            }
            egui::ScrollArea::vertical()
                .id_salt("gene_set_tfbs_member_results")
                .max_height(220.0)
                .show(ui, |ui| {
                    egui::Grid::new("gene_set_tfbs_member_result_grid")
                        .striped(true)
                        .spacing([12.0, 5.0])
                        .show(ui, |ui| {
                            ui.strong("Member");
                            ui.strong("Sequence");
                            ui.strong("Retained hits");
                            ui.strong("Motifs scanned");
                            ui.strong("Complete");
                            ui.end_row();
                            for member_report in tfbs_scan.member_reports.iter().take(100) {
                                ui.label(&member_report.stable_member_id);
                                ui.monospace(&member_report.seq_id);
                                ui.label(member_report.report.matched_hit_count.to_string());
                                ui.label(member_report.report.motifs_scanned.len().to_string());
                                ui.label(
                                    (!member_report.report.truncated_at_max_hits
                                        && tfbs_scan.effective_motif_ids.iter().all(|tf_id| {
                                            member_report.report.motifs_scanned.contains(tf_id)
                                        }))
                                    .to_string(),
                                );
                                ui.end_row();
                            }
                        });
                });
            for member_report in tfbs_scan.member_reports.iter().take(20) {
                egui::CollapsingHeader::new(format!(
                    "{} hits ({})",
                    member_report.stable_member_id, member_report.report.matched_hit_count
                ))
                .default_open(false)
                .show(ui, |ui| {
                    egui::ScrollArea::both()
                        .id_salt(("gene_set_tfbs_hits", &member_report.stable_member_id))
                        .max_height(220.0)
                        .show(ui, |ui| {
                            egui::Grid::new((
                                "gene_set_tfbs_hit_grid",
                                &member_report.stable_member_id,
                            ))
                            .striped(true)
                            .spacing([12.0, 5.0])
                            .show(ui, |ui| {
                                ui.strong("TF");
                                ui.strong("Source range");
                                ui.strong("Strand");
                                ui.strong("LLR bits");
                                ui.strong("Quantile");
                                ui.end_row();
                                for hit in member_report.report.rows.iter().take(100) {
                                    ui.monospace(&hit.tf_id);
                                    ui.monospace(format!(
                                        "{}..{}",
                                        hit.source_match_start_0based,
                                        hit.source_match_end_0based_exclusive
                                    ));
                                    ui.label(if hit.forward_strand { "+" } else { "-" });
                                    ui.label(format!("{:.3}", hit.llr_bits));
                                    ui.label(format!("{:.3}", hit.llr_quantile));
                                    ui.end_row();
                                }
                            });
                        });
                });
            }
        }

        if let CollectionOperationPayload::Digest(digest) = &result.payload {
            ui.separator();
            ui.strong(format!(
                "Digest fragments: {} planned, {} created",
                digest.total_planned_fragment_count, digest.total_created_fragment_count
            ));
            ui.small(format!(
                "enzymes={} | mode={} | aggregate counts={}",
                digest.effective_enzymes.join(", "),
                if digest.collection_operation.applied {
                    "applied"
                } else {
                    "preview"
                },
                if digest.aggregate_counts_complete {
                    "complete"
                } else {
                    "incomplete"
                }
            ));
            ui.monospace(format!("plan={}", digest.plan_fingerprint_sha256));
            if ui
                .button("Copy digest JSON")
                .on_hover_text("Copy the complete portable collection digest report")
                .clicked()
                && let Ok(json) = serde_json::to_string_pretty(digest.as_ref())
            {
                ui.ctx().copy_text(json);
            }
            egui::ScrollArea::vertical()
                .id_salt("gene_set_digest_member_results")
                .max_height(240.0)
                .show(ui, |ui| {
                    egui::Grid::new("gene_set_digest_member_result_grid")
                        .striped(true)
                        .spacing([12.0, 5.0])
                        .show(ui, |ui| {
                            ui.strong("Member");
                            ui.strong("Source");
                            ui.strong("Fragments");
                            ui.strong("Sequence IDs");
                            ui.strong("Open");
                            ui.end_row();
                            for member in digest.member_reports.iter().take(100) {
                                ui.label(&member.stable_member_id);
                                ui.monospace(&member.source_seq_id);
                                ui.label(member.fragments.len().to_string());
                                let all_ids = member
                                    .fragments
                                    .iter()
                                    .map(|fragment| fragment.seq_id.as_str())
                                    .collect::<Vec<_>>();
                                let visible_ids = all_ids
                                    .iter()
                                    .take(3)
                                    .copied()
                                    .collect::<Vec<_>>()
                                    .join(", ");
                                ui.monospace(if all_ids.len() > 3 {
                                    format!("{visible_ids}, ...")
                                } else {
                                    visible_ids
                                })
                                .on_hover_text(all_ids.join("\n"));
                                ui.add_enabled_ui(
                                    member
                                        .fragments
                                        .iter()
                                        .any(|fragment| fragment.materialized),
                                    |ui| {
                                        ui.menu_button("Open", |ui| {
                                            for fragment in member
                                                .fragments
                                                .iter()
                                                .filter(|fragment| fragment.materialized)
                                            {
                                                if ui.button(&fragment.seq_id).clicked() {
                                                    open_digest_sequence =
                                                        Some(fragment.seq_id.clone());
                                                    ui.close();
                                                }
                                            }
                                        });
                                    },
                                );
                                ui.end_row();
                            }
                        });
                });
        }

        if let CollectionOperationPayload::CreatePool(pool) = &result.payload {
            ui.separator();
            ui.strong(format!(
                "Physical pool: {} member aliquot(s) {}",
                pool.planned_member_count,
                if pool.collection_operation.applied {
                    "materialized"
                } else {
                    "planned"
                }
            ));
            ui.small(format!(
                "container={} | name={} | source containers retained={}",
                pool.created_container_id
                    .as_deref()
                    .unwrap_or(pool.planned_container_id.as_str()),
                pool.container_name,
                pool.source_containers_retained
            ));
            ui.monospace(format!("plan={}", pool.plan_fingerprint_sha256));
            if ui
                .button("Copy pool creation JSON")
                .on_hover_text("Copy the complete portable gene-set pool creation report")
                .clicked()
                && let Ok(json) = serde_json::to_string_pretty(pool.as_ref())
            {
                ui.ctx().copy_text(json);
            }
            egui::ScrollArea::vertical()
                .id_salt("gene_set_pool_member_results")
                .max_height(240.0)
                .show(ui, |ui| {
                    egui::Grid::new("gene_set_pool_member_result_grid")
                        .striped(true)
                        .spacing([12.0, 5.0])
                        .show(ui, |ui| {
                            ui.strong("Member");
                            ui.strong("Source tube");
                            ui.strong("Source sequence");
                            ui.strong("Pool aliquot");
                            ui.strong("Open");
                            ui.end_row();
                            for member in pool.member_reports.iter().take(100) {
                                ui.label(&member.stable_member_id);
                                ui.monospace(&member.source_container_id);
                                ui.monospace(&member.source_seq_id);
                                ui.monospace(&member.output_seq_id);
                                if ui
                                    .add_enabled(member.materialized, egui::Button::new("Open"))
                                    .on_disabled_hover_text(
                                        "Previewed aliquots are not project sequences until apply",
                                    )
                                    .clicked()
                                {
                                    open_digest_sequence = Some(member.output_seq_id.clone());
                                }
                                ui.end_row();
                            }
                        });
                });
        }

        if let CollectionOperationPayload::PromoterCohort(promoter) = &result.payload {
            ui.separator();
            ui.strong(format!(
                "Promoter windows: {} returned from {} requested member(s)",
                promoter.returned_window_count, promoter.requested_member_count
            ));
            ui.small(format!(
                "genome={} | upstream={} bp | downstream={} bp | relationship={}",
                promoter.genome_id,
                promoter.upstream_bp,
                promoter.downstream_bp,
                gene_set_relationship_label(promoter.relationship)
            ));
            if ui
                .button("Copy promoter cohort JSON")
                .on_hover_text("Copy the complete portable promoter-cohort report")
                .clicked()
                && let Ok(json) = serde_json::to_string_pretty(promoter.as_ref())
            {
                ui.ctx().copy_text(json);
            }
            egui::ScrollArea::both()
                .id_salt("gene_set_promoter_window_rows")
                .max_height(240.0)
                .show(ui, |ui| {
                    egui::Grid::new("gene_set_promoter_window_grid")
                        .striped(true)
                        .spacing([12.0, 5.0])
                        .show(ui, |ui| {
                            ui.strong("Gene");
                            ui.strong("Transcript");
                            ui.strong("TSS");
                            ui.strong("Promoter window");
                            ui.strong("Strand");
                            ui.end_row();
                            for window in &promoter.windows {
                                ui.label(&window.display_label);
                                ui.monospace(&window.transcript_id);
                                ui.monospace(window.tss_1based.to_string());
                                ui.monospace(format!(
                                    "{}:{}-{}",
                                    window.chromosome,
                                    window.promoter_start_1based,
                                    window.promoter_end_1based
                                ));
                                ui.label(&window.strand);
                                ui.end_row();
                            }
                        });
                });
            if !promoter.unresolved_members.is_empty() {
                egui::CollapsingHeader::new(format!(
                    "Unresolved promoter members ({})",
                    promoter.unresolved_members.len()
                ))
                .default_open(true)
                .show(ui, |ui| {
                    for member in &promoter.unresolved_members {
                        ui.small(format!("{}: {}", member.query, member.reason));
                    }
                });
            }
            if !promoter.relationship_flags.is_empty() {
                egui::CollapsingHeader::new(format!(
                    "Relationship evidence flags ({})",
                    promoter.relationship_flags.len()
                ))
                .default_open(true)
                .show(ui, |ui| {
                    for flag in &promoter.relationship_flags {
                        ui.label(format!("{} | {}", flag.flag_kind, flag.evidence_kind))
                            .on_hover_text(flag.detail.clone());
                    }
                    ui.small(
                        "These flags are non-blocking association evidence, not regulatory proof.",
                    );
                });
            }
            if !promoter.warnings.is_empty() {
                egui::CollapsingHeader::new(format!(
                    "Promoter warnings ({})",
                    promoter.warnings.len()
                ))
                .show(ui, |ui| {
                    for warning in &promoter.warnings {
                        ui.small(format!("- {warning}"));
                    }
                });
            }
        }
        if let Some((seq_id, report_id)) = open_child {
            self.open_sequence_window_for_primer_specificity_report(&seq_id, &report_id);
        }
        if let Some(seq_id) = open_digest_sequence {
            self.open_sequence_window(&seq_id);
        }
    }

    fn render_gene_set_collection_launcher(
        &mut self,
        ui: &mut Ui,
        choice: &GeneSetResolutionChoice,
    ) {
        ui.separator();
        ui.strong("Collection operations");
        ui.small(
            "The operation catalog comes from shared engine lifting policies. Logical gene sets are never silently converted into physical pools or gel lanes.",
        );
        let rows = collection_launcher_rows(CollectionSubjectKind::GeneSetResolution);
        let presented = rows
            .iter()
            .cloned()
            .map(|row| {
                let readiness = self.gene_set_collection_operation_readiness(&row, &choice.report);
                (row, readiness)
            })
            .collect::<Vec<_>>();
        egui::ScrollArea::horizontal()
            .id_salt("gene_set_collection_launcher_scroll")
            .show(ui, |ui| {
                egui::Grid::new("gene_set_collection_launcher_grid")
                    .striped(true)
                    .spacing([14.0, 7.0])
                    .show(ui, |ui| {
                ui.strong("Operation");
                ui.strong("Mode");
                ui.strong("Readiness");
                ui.strong("Result payload");
                ui.end_row();
                for (row, readiness) in &presented {
                    let title = row.title.as_str();
                    if let Some(adapter) = row.adapter
                        && !matches!(
                            readiness,
                            CollectionLauncherReadiness::AdapterUnavailable { .. }
                                | CollectionLauncherReadiness::RequiresMaterialization { .. }
                                | CollectionLauncherReadiness::RequiresPhysicalPool { .. }
                                | CollectionLauncherReadiness::PolicyRejected { .. }
                        )
                    {
                        ui.selectable_value(
                            &mut self.gene_set_inspector.selected_operation,
                            adapter,
                            title,
                        )
                        .on_hover_text(row.description.clone());
                    } else {
                        ui.add_enabled(false, egui::Button::new(title))
                            .on_disabled_hover_text(row.description.clone());
                    }
                    ui.label(row.lifting_mode_label());
                    let readiness_color = if readiness.is_ready() {
                        ui.visuals().strong_text_color()
                    } else if matches!(
                        readiness,
                        CollectionLauncherReadiness::NeedsInput { .. }
                            | CollectionLauncherReadiness::NeedsBindings { .. }
                    ) {
                        ui.visuals().warn_fg_color
                    } else {
                        ui.visuals().error_fg_color
                    };
                    let mut readiness_text = readiness.label().to_string();
                    if let Some(reason) = readiness.rejection_reason() {
                        readiness_text.push_str(&format!(" ({})", reason.as_str()));
                    }
                    let response = ui.colored_label(readiness_color, readiness_text);
                    if let Some(detail) = readiness.detail() {
                        response.on_hover_text(detail);
                    }
                    ui.monospace(row.result_payload_kind());
                    ui.end_row();
                    }
                });
            });

        ui.separator();
        match self.gene_set_inspector.selected_operation {
            CollectionLauncherAdapter::PrimerSpecificity => {
                self.render_gene_set_specificity_form(ui, choice)
            }
            CollectionLauncherAdapter::RestrictionScan => {
                self.render_gene_set_restriction_scan_form(ui, choice)
            }
            CollectionLauncherAdapter::TfbsScan => self.render_gene_set_tfbs_scan_form(ui, choice),
            CollectionLauncherAdapter::Digest => self.render_gene_set_digest_form(ui, choice),
            CollectionLauncherAdapter::CreatePool => {
                self.render_gene_set_create_pool_form(ui, choice)
            }
            CollectionLauncherAdapter::PromoterCohort => self.render_gene_set_promoter_form(ui),
        }
        let readiness = presented
            .iter()
            .find(|(row, _)| row.adapter == Some(self.gene_set_inspector.selected_operation))
            .map(|(_, readiness)| readiness.clone())
            .unwrap_or_else(|| CollectionLauncherReadiness::AdapterUnavailable {
                detail: "The selected operation is absent from the shared policy registry"
                    .to_string(),
            });
        let running = self.gene_set_inspector.task.is_some();
        let mut run_clicked = false;
        let mut discard_clicked = false;
        ui.horizontal_wrapped(|ui| {
            if ui
                .add_enabled(
                    !running && readiness.is_ready(),
                    egui::Button::new(format!(
                        "Run {}",
                        self.gene_set_inspector.selected_operation.label()
                    )),
                )
                .on_hover_text(
                    readiness
                        .detail()
                        .unwrap_or("Run the shared shell/engine operation in a detached project snapshot"),
                )
                .clicked()
            {
                run_clicked = true;
            }
            if let Some(task) = self.gene_set_inspector.task.as_ref() {
                ui.add(egui::Spinner::new());
                ui.label(format!(
                    "Running #{} {} ({:.1}s)",
                    task.job_id,
                    task.adapter.label(),
                    task.started.elapsed().as_secs_f32()
                ));
                if task.cancel_requested.load(Ordering::Relaxed) {
                    ui.small("Discard requested...");
                } else if ui
                    .button("Discard result")
                    .on_hover_text(
                        "Request that GENtle drop the detached result. Active external work may finish, and output files already written by the shared command are not rolled back.",
                    )
                    .clicked()
                {
                    discard_clicked = true;
                }
            }
        });
        if !readiness.is_ready() {
            ui.small(format!(
                "{}: {}",
                readiness.label(),
                readiness.detail().unwrap_or("not ready")
            ));
        }
        if !self.gene_set_inspector.status.trim().is_empty() {
            ui.small(self.gene_set_inspector.status.clone());
        }
        self.render_gene_set_collection_result(ui);

        if run_clicked {
            self.start_gene_set_collection_operation_task();
        }
        if discard_clicked {
            self.request_gene_set_collection_operation_discard("gene-set inspector");
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
            "Inspect one persisted logical gene set, then choose a shared collection operation with its declared lifting policy.",
        );
        ui.small(
            "Execution success means GENtle produced the requested report. Biological specificity, relationship evidence, and unresolved members remain separate result states.",
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
                self.seed_gene_set_promoter_form_from_selection(false);
            }
        });
        if selection_changed {
            self.gene_set_inspector.last_result = None;
            self.gene_set_inspector
                .digest_form
                .expected_plan_fingerprint_sha256
                .clear();
            self.gene_set_inspector.digest_form.apply_change = false;
            self.gene_set_inspector
                .pool_form
                .expected_plan_fingerprint_sha256
                .clear();
            self.gene_set_inspector.pool_form.apply_change = false;
            self.reconcile_gene_set_inspector_bindings();
            self.seed_gene_set_promoter_form_from_selection(true);
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

        self.render_gene_set_collection_launcher(ui, &choice);
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
    use crate::{
        engine::{BiologicalContext, BiologicalContextRegistry, GeneSetResolvedMember},
        engine_shell::parse_shell_line,
    };

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
    fn gene_set_resolve_parity_genomic_neighbors_allows_zero_flanks() {
        assert_resolve_form_matches_shell(
            GeneSetResolveFormState {
                source_kind: GeneSetResolveSourceKind::GenomicNeighbors,
                anchor_gene: "TP53".to_string(),
                flank_genes: "0".to_string(),
                flank_bp: "0".to_string(),
                ..GeneSetResolveFormState::default()
            },
            "gene-sets resolve --neighbors TP53 --flank-genes 0 --flank-bp 0",
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
    fn gene_set_resolve_parity_random_allows_zero_size() {
        assert_resolve_form_matches_shell(
            GeneSetResolveFormState {
                source_kind: GeneSetResolveSourceKind::Random,
                random_size: "0".to_string(),
                ..GeneSetResolveFormState::default()
            },
            "gene-sets resolve --random-size 0",
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
                random_size: "-1".to_string(),
                ..GeneSetResolveFormState::default()
            })
            .expect_err("negative sample")
            .contains("non-negative integer")
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
        let mut resolution = GeneSetResolutionReport {
            op_id: Some("resolution:genes".to_string()),
            review_status: GeneSetResolutionReviewStatus::Reviewed,
            genome_id: Some("GRCh38".to_string()),
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
        };
        resolution
            .ensure_default_biological_context()
            .expect("seed biological context");
        {
            let mut engine = app.engine.write().expect("engine");
            engine
                .upsert_gene_set_resolution_artifact(resolution)
                .expect("persist resolution");
            engine.state_mut().metadata.insert(
                crate::engine::PRIMER_DESIGN_REPORTS_METADATA_KEY.to_string(),
                serde_json::json!({
                    "schema": "gentle.primer_design_reports.v1",
                    "reports": {
                        "report_1": {
                            "report_id": "report_1",
                            "template": "seq_1",
                            "pair_count": 2,
                            "backend": {"used": "internal"}
                        },
                        "report_2": {
                            "report_id": "report_2",
                            "template": "seq_2",
                            "pair_count": 1,
                            "backend": {"used": "internal"}
                        }
                    }
                }),
            );
        }
        app.refresh_gene_set_inspector_catalog(true);
        app.gene_set_inspector.selected_resolution_id = "resolution:genes".to_string();
        app.gene_set_inspector
            .specificity_form
            .member_bindings
            .insert("GENE1".to_string(), "report_1".to_string());
        app.gene_set_inspector
            .specificity_form
            .member_bindings
            .insert("GENE2".to_string(), "report_2".to_string());
        app.gene_set_inspector.specificity_form.target_genome_id = "GRCh38".to_string();
        app.gene_set_inspector.promoter_form.genome_id = "GRCh38".to_string();
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

    fn promoter_command_projection(command: &ShellCommand) -> serde_json::Value {
        match command {
            ShellCommand::GeneSetsPromoterCohort {
                genome_id,
                source,
                resolution,
                relationship,
                gene_group_catalog_path,
                genome_catalog_path,
                cache_dir,
                upstream_bp,
                downstream_bp,
                allow_draft,
                allow_deprecated,
                output,
            } => serde_json::json!({
                "genome_id": genome_id,
                "source": source,
                "resolution": resolution,
                "relationship": relationship,
                "gene_group_catalog_path": gene_group_catalog_path,
                "genome_catalog_path": genome_catalog_path,
                "cache_dir": cache_dir,
                "upstream_bp": upstream_bp,
                "downstream_bp": downstream_bp,
                "allow_draft": allow_draft,
                "allow_deprecated": allow_deprecated,
                "output": output,
            }),
            other => panic!("unexpected command: {other:?}"),
        }
    }

    #[test]
    fn gene_set_promoter_genome_follows_selection_changes_only() {
        let mut app = seeded_app();
        let mut mouse_report = app
            .selected_gene_set_resolution()
            .expect("selected resolution")
            .report
            .clone();
        mouse_report.op_id = Some("resolution:mouse".to_string());
        mouse_report.genome_id = Some("GRCm39".to_string());
        mouse_report.biological_contexts = BiologicalContextRegistry::default();
        for member in &mut mouse_report.resolved_members {
            member.context_id = None;
        }
        mouse_report
            .ensure_default_biological_context()
            .expect("mouse biological context");
        app.gene_set_inspector
            .resolutions
            .push(GeneSetResolutionChoice {
                report_id: "resolution:mouse".to_string(),
                report: mouse_report,
            });

        app.gene_set_inspector.promoter_form.genome_id = "manual-genome".to_string();
        app.seed_gene_set_promoter_form_from_selection(false);
        assert_eq!(
            app.gene_set_inspector.promoter_form.genome_id, "manual-genome",
            "opening or reloading the same selection preserves an explicit form value"
        );

        app.gene_set_inspector.selected_resolution_id = "resolution:mouse".to_string();
        app.seed_gene_set_promoter_form_from_selection(true);
        assert_eq!(app.gene_set_inspector.promoter_form.genome_id, "GRCm39");
    }

    fn restriction_scan_command_projection(command: &ShellCommand) -> serde_json::Value {
        match command {
            ShellCommand::CollectionsRunRestrictionScan {
                collection_subject,
                member_bindings,
                enzymes,
                max_sites_per_enzyme,
                include_cut_geometry,
                path,
            } => serde_json::json!({
                "collection_subject": collection_subject,
                "member_bindings": member_bindings,
                "enzymes": enzymes,
                "max_sites_per_enzyme": max_sites_per_enzyme,
                "include_cut_geometry": include_cut_geometry,
                "path": path,
            }),
            other => panic!("unexpected command: {other:?}"),
        }
    }

    fn tfbs_scan_command_projection(command: &ShellCommand) -> serde_json::Value {
        match command {
            ShellCommand::CollectionsRunTfbsScan {
                collection_subject,
                member_bindings,
                motifs,
                min_llr_bits,
                min_llr_quantile,
                per_tf_thresholds,
                max_hits_per_member,
                path,
            } => serde_json::json!({
                "collection_subject": collection_subject,
                "member_bindings": member_bindings,
                "motifs": motifs,
                "min_llr_bits": min_llr_bits,
                "min_llr_quantile": min_llr_quantile,
                "per_tf_thresholds": per_tf_thresholds,
                "max_hits_per_member": max_hits_per_member,
                "path": path,
            }),
            other => panic!("unexpected command: {other:?}"),
        }
    }

    fn digest_command_projection(command: &ShellCommand) -> serde_json::Value {
        match command {
            ShellCommand::CollectionsRunDigest {
                collection_subject,
                member_bindings,
                enzymes,
                output_prefix,
                dry_run,
                expected_plan_fingerprint_sha256,
                path,
            } => serde_json::json!({
                "collection_subject": collection_subject,
                "member_bindings": member_bindings,
                "enzymes": enzymes,
                "output_prefix": output_prefix,
                "dry_run": dry_run,
                "expected_plan_fingerprint_sha256": expected_plan_fingerprint_sha256,
                "path": path,
            }),
            other => panic!("unexpected command: {other:?}"),
        }
    }

    fn pool_command_projection(command: &ShellCommand) -> serde_json::Value {
        match command {
            ShellCommand::GeneSetsCreatePool {
                resolution_id,
                member_bindings,
                output_prefix,
                container_name,
                dry_run,
                expected_plan_fingerprint_sha256,
                path,
            } => serde_json::json!({
                "resolution_id": resolution_id,
                "member_bindings": member_bindings,
                "output_prefix": output_prefix,
                "container_name": container_name,
                "dry_run": dry_run,
                "expected_plan_fingerprint_sha256": expected_plan_fingerprint_sha256,
                "path": path,
            }),
            other => panic!("unexpected command: {other:?}"),
        }
    }

    #[test]
    fn collection_operation_launcher_specificity_command_matches_shared_shell_parser_once() {
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
        let ShellCommand::CollectionsRunPrimerSpecificity {
            member_bindings, ..
        } = gui_command
        else {
            unreachable!("specificity adapter returns one collection command")
        };
        assert_eq!(member_bindings.len(), 2);
    }

    #[test]
    fn collection_operation_launcher_restriction_scan_matches_shared_shell_parser_once() {
        let mut app = seeded_app();
        {
            let mut engine = app.engine.write().expect("engine");
            engine.state_mut().sequences.insert(
                "seq_1".to_string(),
                DNAsequence::from_sequence("AAGAATTCCCGG").expect("sequence 1"),
            );
            engine.state_mut().sequences.insert(
                "seq_2".to_string(),
                DNAsequence::from_sequence("TTTGGATCCAAA").expect("sequence 2"),
            );
        }
        app.gene_set_inspector
            .restriction_scan_form
            .member_bindings
            .insert("GENE1".to_string(), "seq_1".to_string());
        app.gene_set_inspector
            .restriction_scan_form
            .member_bindings
            .insert("GENE2".to_string(), "seq_2".to_string());
        app.gene_set_inspector.restriction_scan_form.enzymes = "EcoRI, BamHI".to_string();
        app.gene_set_inspector
            .restriction_scan_form
            .max_sites_per_enzyme = "3".to_string();
        app.gene_set_inspector
            .restriction_scan_form
            .include_cut_geometry = false;

        let gui_command = app
            .gene_set_restriction_scan_shell_command()
            .expect("GUI command");
        let parsed = parse_shell_line(
            "collections run restriction-scan resolution:genes \
             --member-sequence GENE1=seq_1 --member-sequence GENE2=seq_2 \
             --enzyme EcoRI --enzyme BamHI --max-sites-per-enzyme 3 --no-cut-geometry",
        )
        .expect("shared shell command");
        assert_eq!(
            restriction_scan_command_projection(&gui_command),
            restriction_scan_command_projection(&parsed)
        );

        app.gene_set_inspector
            .restriction_scan_form
            .member_bindings
            .remove("GENE2");
        assert!(
            app.gene_set_restriction_scan_shell_command()
                .expect_err("missing binding")
                .contains("missing: GENE2")
        );
    }

    #[test]
    fn collection_operation_launcher_tfbs_scan_matches_shared_shell_parser_once() {
        let mut app = seeded_app();
        app.gene_set_inspector
            .tfbs_scan_form
            .member_bindings
            .insert("GENE1".to_string(), "seq_1".to_string());
        app.gene_set_inspector
            .tfbs_scan_form
            .member_bindings
            .insert("GENE2".to_string(), "seq_2".to_string());
        app.gene_set_inspector.tfbs_scan_form.motifs = "SP1, CTCF".to_string();
        app.gene_set_inspector.tfbs_scan_form.min_llr_bits = "3.5".to_string();
        app.gene_set_inspector.tfbs_scan_form.min_llr_quantile = "0.9".to_string();
        app.gene_set_inspector.tfbs_scan_form.per_tf_min_llr_bits = "SP1=4.5".to_string();
        app.gene_set_inspector
            .tfbs_scan_form
            .per_tf_min_llr_quantile = "CTCF=0.95".to_string();
        app.gene_set_inspector.tfbs_scan_form.max_hits_per_member = "40".to_string();

        let gui_command = app
            .gene_set_tfbs_scan_shell_command()
            .expect("GUI TFBS command");
        let parsed = parse_shell_line(
            "collections run tfbs-scan resolution:genes \
             --member-sequence GENE1=seq_1 --member-sequence GENE2=seq_2 \
             --motif SP1 --motif CTCF --min-llr-bits 3.5 \
             --min-llr-quantile 0.9 --per-tf-min-llr-bits SP1=4.5 \
             --per-tf-min-llr-quantile CTCF=0.95 --max-hits 40",
        )
        .expect("shared shell TFBS command");
        assert_eq!(
            tfbs_scan_command_projection(&gui_command),
            tfbs_scan_command_projection(&parsed)
        );

        app.gene_set_inspector
            .tfbs_scan_form
            .member_bindings
            .remove("GENE2");
        assert!(
            app.gene_set_tfbs_scan_shell_command()
                .expect_err("missing binding")
                .contains("missing: GENE2")
        );
        app.gene_set_inspector
            .tfbs_scan_form
            .member_bindings
            .insert("GENE2".to_string(), "seq_2".to_string());
        app.gene_set_inspector.tfbs_scan_form.motifs.clear();
        assert!(
            app.gene_set_tfbs_scan_shell_command()
                .expect_err("missing motifs")
                .contains("at least one")
        );
    }

    #[test]
    fn collection_operation_launcher_digest_matches_shared_shell_preview_and_apply() {
        let mut app = seeded_app();
        app.gene_set_inspector
            .digest_form
            .member_bindings
            .insert("GENE1".to_string(), "seq_1".to_string());
        app.gene_set_inspector
            .digest_form
            .member_bindings
            .insert("GENE2".to_string(), "seq_2".to_string());
        app.gene_set_inspector.digest_form.enzymes = "EcoRI, BamHI".to_string();
        app.gene_set_inspector.digest_form.output_prefix = "batch".to_string();

        let gui_preview = app
            .gene_set_digest_shell_command()
            .expect("GUI preview command");
        let parsed_preview = parse_shell_line(
            "collections run digest resolution:genes \
             --member-sequence GENE1=seq_1 --member-sequence GENE2=seq_2 \
             --enzyme EcoRI --enzyme BamHI --output-prefix batch",
        )
        .expect("shared shell preview command");
        assert_eq!(
            digest_command_projection(&gui_preview),
            digest_command_projection(&parsed_preview)
        );

        app.gene_set_inspector.digest_form.apply_change = true;
        assert!(
            app.gene_set_digest_shell_command()
                .expect_err("apply requires a preview lock")
                .contains("Preview this digest")
        );
        app.gene_set_inspector
            .digest_form
            .expected_plan_fingerprint_sha256 = "sha256:locked".to_string();
        let gui_apply = app
            .gene_set_digest_shell_command()
            .expect("GUI apply command");
        let parsed_apply = parse_shell_line(
            "collections run digest resolution:genes \
             --member-sequence GENE1=seq_1 --member-sequence GENE2=seq_2 \
             --enzyme EcoRI --enzyme BamHI --output-prefix batch --apply \
             --expected-plan-fingerprint-sha256 sha256:locked",
        )
        .expect("shared shell apply command");
        assert_eq!(
            digest_command_projection(&gui_apply),
            digest_command_projection(&parsed_apply)
        );

        app.gene_set_inspector
            .digest_form
            .member_bindings
            .remove("GENE2");
        assert!(
            app.gene_set_digest_shell_command()
                .expect_err("missing binding")
                .contains("missing: GENE2")
        );
    }

    #[test]
    fn collection_operation_launcher_create_pool_matches_shell_and_requires_complete_bindings() {
        let mut app = seeded_app();
        let (container_1, container_2) = {
            let mut engine = app.engine.write().expect("engine");
            engine
                .apply(Operation::CreateSequenceFromText {
                    sequence_text: "ACGTACGT".to_string(),
                    output_id: Some("pool_source_1".to_string()),
                    name: Some("Source 1".to_string()),
                    circular: false,
                })
                .expect("create source 1");
            engine
                .apply(Operation::CreateSequenceFromText {
                    sequence_text: "TTGGAACC".to_string(),
                    output_id: Some("pool_source_2".to_string()),
                    name: Some("Source 2".to_string()),
                    circular: false,
                })
                .expect("create source 2");
            (
                engine.state().container_state.seq_to_latest_container["pool_source_1"].clone(),
                engine.state().container_state.seq_to_latest_container["pool_source_2"].clone(),
            )
        };
        app.reconcile_gene_set_inspector_bindings();
        app.gene_set_inspector
            .pool_form
            .member_bindings
            .insert("GENE1".to_string(), container_1.clone());
        app.gene_set_inspector
            .pool_form
            .member_bindings
            .insert("GENE2".to_string(), container_2.clone());
        app.gene_set_inspector.pool_form.output_prefix = "cohort".to_string();
        app.gene_set_inspector.pool_form.container_name = "Cohort stock".to_string();

        let gui_preview = app
            .gene_set_create_pool_shell_command()
            .expect("GUI pool preview command");
        let parsed_preview = parse_shell_line(&format!(
            "gene-sets create-pool resolution:genes --member-container GENE1={container_1} \
             --member-container GENE2={container_2} --output-prefix cohort \
             --container-name 'Cohort stock'"
        ))
        .expect("shared shell pool preview command");
        assert_eq!(
            pool_command_projection(&gui_preview),
            pool_command_projection(&parsed_preview)
        );

        let row = collection_launcher_rows(CollectionSubjectKind::GeneSetResolution)
            .into_iter()
            .find(|row| row.adapter == Some(CollectionLauncherAdapter::CreatePool))
            .expect("CreateGeneSetPool launcher row");
        let report = app
            .selected_gene_set_resolution()
            .expect("selected resolution")
            .report
            .clone();
        assert!(matches!(
            app.gene_set_collection_operation_readiness(&row, &report),
            CollectionLauncherReadiness::Ready
        ));

        app.gene_set_inspector.pool_form.apply_change = true;
        assert!(
            app.gene_set_create_pool_shell_command()
                .expect_err("apply requires preview fingerprint")
                .contains("Preview this physical pool")
        );
        app.gene_set_inspector
            .pool_form
            .expected_plan_fingerprint_sha256 = "sha256:locked".to_string();
        let gui_apply = app
            .gene_set_create_pool_shell_command()
            .expect("GUI pool apply command");
        let parsed_apply = parse_shell_line(&format!(
            "gene-sets create-pool resolution:genes --member-container GENE1={container_1} \
             --member-container GENE2={container_2} --output-prefix cohort \
             --container-name 'Cohort stock' --apply \
             --expected-plan-fingerprint-sha256 sha256:locked"
        ))
        .expect("shared shell pool apply command");
        assert_eq!(
            pool_command_projection(&gui_apply),
            pool_command_projection(&parsed_apply)
        );

        app.gene_set_inspector
            .pool_form
            .member_bindings
            .remove("GENE2");
        assert!(
            app.gene_set_create_pool_shell_command()
                .expect_err("missing pool binding")
                .contains("missing: GENE2")
        );

        let mut incomplete = report;
        incomplete.unresolved_member_count = 1;
        assert!(matches!(
            app.gene_set_collection_operation_readiness(&row, &incomplete),
            CollectionLauncherReadiness::NeedsInput { .. }
        ));
    }

    #[test]
    fn collection_operation_launcher_requires_complete_explicit_bindings_and_valid_rank() {
        let mut app = seeded_app();
        app.gene_set_inspector
            .specificity_form
            .member_bindings
            .remove("GENE2");
        assert!(
            app.gene_set_specificity_shell_command()
                .expect_err("missing binding")
                .contains("missing: GENE2")
        );

        app.gene_set_inspector
            .specificity_form
            .member_bindings
            .insert("GENE2".to_string(), "report_2".to_string());
        app.gene_set_inspector.specificity_form.pair_rank_1based = "2".to_string();
        assert!(
            app.gene_set_specificity_shell_command()
                .expect_err("rank unavailable for report_2")
                .contains("rank 2 is unavailable")
        );
    }

    #[test]
    fn collection_operation_launcher_promoter_command_matches_shared_shell_parser() {
        let mut app = seeded_app();
        app.gene_set_inspector.promoter_form = GeneSetPromoterFormState {
            genome_id: "GRCh38".to_string(),
            relationship: GeneSetCohortRelationship::CoRegulated,
            upstream_bp: "750".to_string(),
            downstream_bp: "125".to_string(),
            gene_group_catalog_path: "groups.json".to_string(),
            genome_catalog_path: "genomes.json".to_string(),
            cache_dir: "cache".to_string(),
            output_path: "cohort.json".to_string(),
        };
        let gui_command = app.gene_set_promoter_shell_command().expect("GUI command");
        let resolution_json = serde_json::to_string(
            &app.selected_gene_set_resolution()
                .expect("selected resolution")
                .report,
        )
        .expect("resolution JSON");
        let parsed = parse_shell_line(&format!(
            "gene-sets promoter-cohort GRCh38 --resolution '{resolution_json}' \
             --relationship co-regulated --upstream-bp 750 --downstream-bp 125 \
             --catalog groups.json --genome-catalog genomes.json --cache-dir cache \
             --output cohort.json"
        ))
        .expect("shared shell command");

        assert_eq!(
            promoter_command_projection(&gui_command),
            promoter_command_projection(&parsed)
        );
        let ShellCommand::GeneSetsPromoterCohort {
            resolution: Some(resolution),
            ..
        } = &gui_command
        else {
            panic!("promoter adapter must inline one complete resolution")
        };
        assert_eq!(resolution.resolved_members.len(), 2);
        assert_eq!(
            app.gene_set_inspector.specificity_form.target_genome_id, "GRCh38",
            "promoter and specificity parameters remain separate form state"
        );
    }

    #[test]
    fn collection_operation_launcher_missing_bindings_is_not_materialization() {
        let mut app = seeded_app();
        app.gene_set_inspector
            .specificity_form
            .member_bindings
            .remove("GENE2");
        let report = app
            .selected_gene_set_resolution()
            .expect("selected resolution")
            .report
            .clone();
        let row = collection_launcher_rows(CollectionSubjectKind::GeneSetResolution)
            .into_iter()
            .find(|row| row.adapter == Some(CollectionLauncherAdapter::PrimerSpecificity))
            .expect("specificity row");
        assert!(matches!(
            app.gene_set_collection_operation_readiness(&row, &report),
            CollectionLauncherReadiness::NeedsBindings { .. }
        ));
    }

    #[test]
    fn collection_operation_launcher_mixed_context_readiness_matches_engine_rejection() {
        let context = |context_id: &str, genome_id: &str| BiologicalContext {
            context_id: context_id.to_string(),
            genome_id: Some(genome_id.to_string()),
            ..BiologicalContext::default()
        };
        let mut mixed = GeneSetResolutionReport {
            op_id: Some("resolution:mixed".to_string()),
            resolved_member_count: 2,
            resolved_members: vec![
                GeneSetResolvedMember {
                    dedup_key: "GENE1".to_string(),
                    symbol: "GENE1".to_string(),
                    context_id: Some("human".to_string()),
                    ..GeneSetResolvedMember::default()
                },
                GeneSetResolvedMember {
                    dedup_key: "GENE2".to_string(),
                    symbol: "GENE2".to_string(),
                    context_id: Some("mouse".to_string()),
                    ..GeneSetResolvedMember::default()
                },
            ],
            biological_contexts: BiologicalContextRegistry {
                contexts: vec![context("human", "GRCh38"), context("mouse", "GRCm39")],
                default_context_id: Some("human".to_string()),
            },
            ..GeneSetResolutionReport::default()
        };
        let mut app = seeded_app();
        app.gene_set_inspector.resolutions = vec![GeneSetResolutionChoice {
            report_id: "resolution:mixed".to_string(),
            report: mixed.clone(),
        }];
        app.gene_set_inspector.selected_resolution_id = "resolution:mixed".to_string();
        let row = collection_launcher_rows(CollectionSubjectKind::GeneSetResolution)
            .into_iter()
            .find(|row| row.adapter == Some(CollectionLauncherAdapter::PrimerSpecificity))
            .expect("specificity row");
        let readiness = app.gene_set_collection_operation_readiness(&row, &mixed);
        assert!(matches!(
            readiness,
            CollectionLauncherReadiness::PolicyRejected {
                reason: CollectionLiftRejectionReason::MixedBiologicalContext,
                ..
            }
        ));

        app.gene_set_inspector
            .restriction_scan_form
            .member_bindings
            .insert("GENE1".to_string(), "seq_1".to_string());
        app.gene_set_inspector
            .restriction_scan_form
            .member_bindings
            .insert("GENE2".to_string(), "seq_2".to_string());
        let restriction = collection_launcher_rows(CollectionSubjectKind::GeneSetResolution)
            .into_iter()
            .find(|row| row.adapter == Some(CollectionLauncherAdapter::RestrictionScan))
            .expect("restriction-scan row");
        assert_eq!(
            app.gene_set_collection_operation_readiness(&restriction, &mixed),
            CollectionLauncherReadiness::Ready,
            "restriction scanning has no genome-scoped parameter and is context agnostic"
        );

        app.gene_set_inspector
            .tfbs_scan_form
            .member_bindings
            .insert("GENE1".to_string(), "seq_1".to_string());
        app.gene_set_inspector
            .tfbs_scan_form
            .member_bindings
            .insert("GENE2".to_string(), "seq_2".to_string());
        app.gene_set_inspector.tfbs_scan_form.motifs = "SP1".to_string();
        let tfbs = collection_launcher_rows(CollectionSubjectKind::GeneSetResolution)
            .into_iter()
            .find(|row| row.adapter == Some(CollectionLauncherAdapter::TfbsScan))
            .expect("TFBS-scan row");
        assert_eq!(
            app.gene_set_collection_operation_readiness(&tfbs, &mixed),
            CollectionLauncherReadiness::Ready,
            "TFBS hit scanning has no genome-scoped parameter and is context agnostic"
        );

        app.gene_set_inspector
            .digest_form
            .member_bindings
            .insert("GENE1".to_string(), "seq_1".to_string());
        app.gene_set_inspector
            .digest_form
            .member_bindings
            .insert("GENE2".to_string(), "seq_2".to_string());
        app.gene_set_inspector.digest_form.enzymes = "EcoRI".to_string();
        let digest = collection_launcher_rows(CollectionSubjectKind::GeneSetResolution)
            .into_iter()
            .find(|row| row.adapter == Some(CollectionLauncherAdapter::Digest))
            .expect("digest row");
        assert_eq!(
            app.gene_set_collection_operation_readiness(&digest, &mixed),
            CollectionLauncherReadiness::Ready,
            "restriction digestion is context agnostic"
        );

        let mut engine = GentleEngine::new();
        engine
            .upsert_gene_set_resolution_artifact(std::mem::take(&mut mixed))
            .expect("persist mixed resolution");
        let error = engine
            .apply(Operation::AssessPrimerPairSpecificityCollection {
                collection_subject: CollectionSubjectRef::GeneSetResolution {
                    report_id: "resolution:mixed".to_string(),
                },
                member_bindings: vec![],
                pair_rank: Some(1),
                pair_index: None,
                target_genome_id: "GRCh38".to_string(),
                policy: PrimerSpecificityPolicy::default(),
                catalog_path: None,
                cache_dir: None,
                path: None,
            })
            .expect_err("mixed context must fail before report bindings or BLAST");
        assert!(error.message.contains("mixed_biological_context"));
    }

    #[test]
    fn collection_operation_launcher_promoter_readiness_matches_legacy_context_normalization() {
        let app = seeded_app();
        let mut legacy = app
            .selected_gene_set_resolution()
            .expect("selected resolution")
            .report
            .clone();
        legacy.biological_contexts = BiologicalContextRegistry::default();
        for member in &mut legacy.resolved_members {
            member.context_id = None;
        }
        let rows = collection_launcher_rows(CollectionSubjectKind::GeneSetResolution);
        let promoter = rows
            .iter()
            .find(|row| row.adapter == Some(CollectionLauncherAdapter::PromoterCohort))
            .expect("promoter row");
        assert_eq!(
            app.gene_set_collection_operation_readiness(promoter, &legacy),
            CollectionLauncherReadiness::Ready,
            "promoter engine promotes legacy report-level genome context"
        );

        let specificity = rows
            .iter()
            .find(|row| row.adapter == Some(CollectionLauncherAdapter::PrimerSpecificity))
            .expect("specificity row");
        assert!(matches!(
            app.gene_set_collection_operation_readiness(specificity, &legacy),
            CollectionLauncherReadiness::PolicyRejected {
                reason: CollectionLiftRejectionReason::MissingBiologicalContext,
                ..
            }
        ));
    }

    #[test]
    fn collection_operation_launcher_reports_unimplemented_context_as_adapter_gap() {
        let app = seeded_app();
        let report = &app
            .selected_gene_set_resolution()
            .expect("selected resolution")
            .report;
        let readiness = GENtleApp::gene_set_collection_context_readiness(
            report,
            CollectionContextRequirement::Partitionable,
            "GRCh38",
        );
        assert!(matches!(
            readiness,
            CollectionLauncherReadiness::AdapterUnavailable { detail }
                if detail.contains("Partitionable")
        ));
    }

    #[test]
    fn collection_operation_launcher_promoter_subject_identity_and_success_refresh_are_preserved() {
        let mut app = seeded_app();
        app.gene_set_inspector.selected_operation = CollectionLauncherAdapter::PromoterCohort;
        app.gene_set_inspector.promoter_form.upstream_bp = "777".to_string();
        app.gene_set_inspector.promoter_form.downstream_bp = "123".to_string();
        let selected_report = app
            .selected_gene_set_resolution()
            .expect("selected resolution")
            .report
            .clone();
        let expected_subject = CollectionSubjectRef::GeneSetResolution {
            report_id: GentleEngine::gene_set_resolution_artifact_id(&selected_report),
        };
        let mut promoter_report = GeneSetPromoterCohortReport {
            genome_id: "GRCh38".to_string(),
            upstream_bp: 777,
            downstream_bp: 123,
            gene_set_resolution: selected_report,
            requested_member_count: 2,
            ..GeneSetPromoterCohortReport::default()
        };
        let collection_report = GentleEngine::build_gene_set_promoter_collection_operation_report(
            &promoter_report,
            "promoter_cohort:test",
            "promoter:test",
            "run:test",
        )
        .expect("build production collection wrapper");
        assert_eq!(collection_report.collection_subject, expected_subject);
        promoter_report.collection_operation = Some(Box::new(collection_report));
        let parsed = {
            let engine = app.engine.read().expect("engine");
            GENtleApp::collection_operation_task_result_from_output(
                &engine,
                CollectionLauncherAdapter::PromoterCohort,
                &expected_subject,
                serde_json::to_value(&promoter_report).expect("promoter output"),
            )
            .expect("normalize promoter output")
        };

        let mut mismatched = promoter_report.clone();
        mismatched
            .collection_operation
            .as_mut()
            .expect("wrapper")
            .collection_subject = CollectionSubjectRef::GeneSetResolution {
            report_id: "resolution:other".to_string(),
        };
        let mismatch_error = {
            let engine = app.engine.read().expect("engine");
            GENtleApp::collection_operation_task_result_from_output(
                &engine,
                CollectionLauncherAdapter::PromoterCohort,
                &expected_subject,
                serde_json::to_value(mismatched).expect("mismatched output"),
            )
            .expect_err("mismatched source identity must fail")
        };
        assert!(mismatch_error.message.contains("expected"));

        let (tx, receiver) = mpsc::channel();
        let runtime_frame = GENtleApp::push_runtime_background_job_frame(
            BackgroundJobKind::CollectionOperation,
            10,
            "test promoter",
        );
        app.gene_set_inspector.task = Some(GeneSetCollectionOperationTask {
            job_id: 10,
            adapter: CollectionLauncherAdapter::PromoterCohort,
            resolution_id: "resolution:genes".to_string(),
            started: Instant::now(),
            cancel_requested: Arc::new(AtomicBool::new(false)),
            runtime_frame,
            receiver,
        });
        tx.send(Ok(parsed)).expect("send promoter result");
        app.poll_gene_set_collection_operation_task(&egui::Context::default());

        assert_eq!(
            app.gene_set_inspector.selected_resolution_id,
            "resolution:genes"
        );
        assert_eq!(
            app.gene_set_inspector.specificity_form.member_bindings["GENE1"],
            "report_1"
        );
        assert_eq!(app.gene_set_inspector.promoter_form.upstream_bp, "777");
        assert_eq!(app.gene_set_inspector.promoter_form.downstream_bp, "123");
        assert!(matches!(
            app.gene_set_inspector.last_result.as_ref(),
            Some(CollectionOperationTaskResult {
                payload: CollectionOperationPayload::PromoterCohort(_),
                ..
            })
        ));
    }

    #[test]
    fn collection_operation_launcher_discard_drops_an_already_queued_success() {
        let mut app = seeded_app();
        app.gene_set_inspector.selected_operation = CollectionLauncherAdapter::PromoterCohort;
        app.gene_set_inspector.promoter_form.upstream_bp = "901".to_string();
        app.gene_set_inspector.promoter_form.downstream_bp = "101".to_string();
        let (tx, receiver) = mpsc::channel();
        let cancel_requested = Arc::new(AtomicBool::new(false));
        app.gene_set_inspector.task = Some(GeneSetCollectionOperationTask {
            job_id: 11,
            adapter: CollectionLauncherAdapter::PromoterCohort,
            resolution_id: "resolution:genes".to_string(),
            started: Instant::now(),
            cancel_requested: cancel_requested.clone(),
            runtime_frame: GENtleApp::push_runtime_background_job_frame(
                BackgroundJobKind::CollectionOperation,
                11,
                "test discard",
            ),
            receiver,
        });
        tx.send(Ok(CollectionOperationTaskResult {
            adapter: CollectionLauncherAdapter::PrimerSpecificity,
            report: CollectionOperationReport {
                report_id: "collection:discarded".to_string(),
                collection_subject: CollectionSubjectRef::GeneSetResolution {
                    report_id: "resolution:genes".to_string(),
                },
                ..CollectionOperationReport::default()
            },
            payload: CollectionOperationPayload::PrimerSpecificity {
                child_reports: BTreeMap::new(),
            },
        }))
        .expect("queue successful result");
        app.request_gene_set_collection_operation_discard("test");
        app.poll_gene_set_collection_operation_task(&egui::Context::default());

        assert!(cancel_requested.load(Ordering::Relaxed));
        assert!(app.gene_set_inspector.task.is_none());
        assert!(app.gene_set_inspector.last_result.is_none());
        assert!(
            app.gene_set_inspector
                .status
                .contains("no result was displayed")
        );
        assert_eq!(
            app.gene_set_inspector.selected_operation,
            CollectionLauncherAdapter::PromoterCohort
        );
        assert_eq!(app.gene_set_inspector.promoter_form.upstream_bp, "901");
        assert_eq!(app.gene_set_inspector.promoter_form.downstream_bp, "101");
    }

    #[test]
    fn collection_operation_launcher_failure_preserves_selection_and_parameters() {
        let mut app = seeded_app();
        app.gene_set_inspector.selected_operation = CollectionLauncherAdapter::PromoterCohort;
        app.gene_set_inspector.promoter_form.upstream_bp = "811".to_string();
        app.gene_set_inspector.promoter_form.downstream_bp = "211".to_string();
        let (tx, receiver) = mpsc::channel();
        app.gene_set_inspector.task = Some(GeneSetCollectionOperationTask {
            job_id: 12,
            adapter: CollectionLauncherAdapter::PromoterCohort,
            resolution_id: "resolution:genes".to_string(),
            started: Instant::now(),
            cancel_requested: Arc::new(AtomicBool::new(false)),
            runtime_frame: GENtleApp::push_runtime_background_job_frame(
                BackgroundJobKind::CollectionOperation,
                12,
                "test failure",
            ),
            receiver,
        });
        tx.send(Err(EngineError::invalid_input("synthetic failure")))
            .expect("send failed result");
        app.poll_gene_set_collection_operation_task(&egui::Context::default());

        assert!(app.gene_set_inspector.task.is_none());
        assert!(app.gene_set_inspector.last_result.is_none());
        assert!(app.gene_set_inspector.status.contains("synthetic failure"));
        assert_eq!(
            app.gene_set_inspector.selected_operation,
            CollectionLauncherAdapter::PromoterCohort
        );
        assert_eq!(app.gene_set_inspector.promoter_form.upstream_bp, "811");
        assert_eq!(app.gene_set_inspector.promoter_form.downstream_bp, "211");
    }

    #[test]
    fn collection_operation_launcher_specificity_completion_preserves_biology_states() {
        let mut app = seeded_app();
        let (tx, receiver) = mpsc::channel();
        let runtime_frame = GENtleApp::push_runtime_background_job_frame(
            BackgroundJobKind::CollectionOperation,
            9,
            "test",
        );
        app.gene_set_inspector.task = Some(GeneSetCollectionOperationTask {
            job_id: 9,
            adapter: CollectionLauncherAdapter::PrimerSpecificity,
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
        tx.send(Ok(CollectionOperationTaskResult {
            adapter: CollectionLauncherAdapter::PrimerSpecificity,
            report: CollectionOperationReport {
                report_id: "collection:1".to_string(),
                collection_subject: CollectionSubjectRef::GeneSetResolution {
                    report_id: "resolution:genes".to_string(),
                },
                per_member_status: vec![crate::engine::CollectionMemberStatusRow {
                    outcome: CollectionMemberOutcome::Succeeded,
                    produced_report_ids: vec!["specificity:1".to_string()],
                    ..crate::engine::CollectionMemberStatusRow::default()
                }],
                ..CollectionOperationReport::default()
            },
            payload: CollectionOperationPayload::PrimerSpecificity {
                child_reports: BTreeMap::from([("specificity:1".to_string(), child)]),
            },
        }))
        .expect("send result");

        app.poll_gene_set_collection_operation_task(&egui::Context::default());

        assert!(app.gene_set_inspector.task.is_none());
        assert!(
            app.gene_set_inspector
                .status
                .contains("1 executed, 0 failed")
        );
        let Some(CollectionOperationTaskResult {
            payload: CollectionOperationPayload::PrimerSpecificity { child_reports },
            ..
        }) = app.gene_set_inspector.last_result.as_ref()
        else {
            panic!("specificity result")
        };
        assert_eq!(child_reports["specificity:1"].status, "fail");
        assert!(
            !child_reports["specificity:1"].specificity_pass,
            "execution success must not be conflated with biological pass"
        );
    }
}
