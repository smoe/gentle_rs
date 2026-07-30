//! Binding-aware gene-set collection inspector.
//!
//! This GUI is deliberately thin: it selects one persisted gene-set
//! resolution, requires an explicit primer-report binding for every logical
//! member, and executes the shared `collections run primer-specificity`
//! command on a detached engine snapshot.

use super::*;
use crate::{
    engine::{
        CollectionMemberOutcome, CollectionOperationReport, CollectionSubjectRef,
        GeneSetResolutionReport, GeneSetResolutionReviewStatus, PrimerDesignReportSummary,
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

pub(super) struct GeneSetInspectorUiState {
    pub(super) show_panel: bool,
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
    task: Option<GeneSetSpecificityTask>,
    last_report: Option<CollectionOperationReport>,
    child_reports: BTreeMap<String, PrimerSpecificityChildPresentation>,
}

impl Default for GeneSetInspectorUiState {
    fn default() -> Self {
        Self {
            show_panel: false,
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
        self.refresh_gene_set_inspector_catalog(true);
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

    fn render_gene_set_inspector_contents(&mut self, ui: &mut Ui) {
        let close_hover = Self::specialist_window_close_hover_text("Gene Set Inspector");
        if self.render_specialist_window_nav_with_close(ui, Some(("Close", close_hover.as_str()))) {
            self.gene_set_inspector.show_panel = false;
            return;
        }
        ui.heading("Gene Set Inspector");
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
                "Resolve and persist a set with `gene-sets resolve ...` in GUI Shell, CLI, MCP, or another shared adapter, then reload this inspector.",
            );
            if !self.gene_set_inspector.status.trim().is_empty() {
                ui.separator();
                ui.small(self.gene_set_inspector.status.clone());
            }
            return;
        };

        ui.horizontal_wrapped(|ui| {
            ui.small(format!(
                "review={} | namespace={} | organism={} | resolved={} | unresolved={}",
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
