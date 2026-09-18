//! Thin, session-scoped study review over shared planning, execution and dossier routes.
//! Approvals are never persisted with a window or inferred from a saved plan.

use super::*;
use crate::command_execution::CommandExecutionService;
use crate::digest_utils::sha256_prefixed_bytes;
use crate::engine::{GeneIsoformAssayStudyPlanReport, GeneIsoformAssayStudyPlanRequest};
use crate::engine_shell::{ShellProgressCallback, ShellRunResult};
use std::path::PathBuf;

#[derive(Clone, Copy, Debug, PartialEq)]
enum StudyJobKind {
    ReadRequest,
    ReadPlan,
    PanelSummaries,
    Normalize,
    Plan,
    Execute,
    Publish,
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::engine::{ProjectState, TranscriptAssayPanelReport};
    use serde_json::json;

    fn area() -> MainAreaDna {
        let mut engine = GentleEngine::new();
        engine
            .apply(Operation::LoadFile {
                path:
                    "test_files/fixtures/transcript_assay_panel/patz1/patz1_assay_minus_strand.gb"
                        .into(),
                as_id: Some("patz1_transcript_assay_demo".into()),
            })
            .unwrap();
        engine
            .apply(Operation::SetParameter {
                name: "primer_design_backend".into(),
                value: json!("internal"),
            })
            .unwrap();
        let dna = engine.state().sequences["patz1_transcript_assay_demo"].clone();
        MainAreaDna::new(
            dna,
            Some("patz1_transcript_assay_demo".into()),
            Some(Arc::new(RwLock::new(engine))),
        )
    }

    fn drain(area: &mut MainAreaDna) {
        let ctx = egui::Context::default();
        let deadline = Instant::now() + Duration::from_secs(20);
        while area.gene_assay_study_ui.job.is_some() {
            assert!(Instant::now() < deadline, "study worker timed out");
            area.poll_gene_assay_study_task(&ctx);
            std::thread::sleep(Duration::from_millis(5));
        }
    }

    fn planned_area(directory: &std::path::Path) -> MainAreaDna {
        let mut area = area();
        area.gene_assay_study_ui.request_json =
            include_str!("../../docs/tutorial/inputs/transcript_assay_followup_study.json").into();
        area.normalize_study_draft();
        drain(&mut area);
        assert!(
            area.gene_assay_study_ui.normalized.is_some(),
            "{}",
            area.gene_assay_study_ui.status
        );
        area.gene_assay_study_ui.output_directory = directory.to_string_lossy().into();
        // Merely having a normalized request does not authorize planning.
        area.plan_reviewed_study();
        assert!(area.gene_assay_study_ui.job.is_none());
        area.gene_assay_study_ui.planning_reviewed = true;
        area.plan_reviewed_study();
        drain(&mut area);
        assert!(
            area.gene_assay_study_ui.plan.is_some(),
            "{}",
            area.gene_assay_study_ui.status
        );
        area
    }

    #[test]
    fn study_normalize_plan_matches_shared_engine_and_preserves_two_reviews() {
        let temp = tempfile::tempdir().unwrap();
        let mut area = planned_area(&temp.path().join("study"));
        let plan = area.gene_assay_study_ui.plan.clone().unwrap();
        let mut engine = area
            .engine
            .as_ref()
            .unwrap()
            .read()
            .unwrap()
            .clone_without_history();
        let direct = engine
            .apply(Operation::PlanGeneIsoformAssayStudy {
                request: plan.normalized_request.clone(),
                path: None,
                workflow_path: None,
            })
            .unwrap()
            .gene_isoform_assay_study_plan
            .unwrap();
        assert_eq!(plan.as_ref(), direct.as_ref());
        assert!(!area.gene_assay_study_ui.execution_reviewed);
        area.execute_reviewed_study();
        assert!(area.gene_assay_study_ui.job.is_none());
        assert!(
            area.engine
                .as_ref()
                .unwrap()
                .read()
                .unwrap()
                .list_transcript_assay_panel_reports()
                .is_empty()
        );
        assert!(
            plan.evidence_summary
                .missing_evidence
                .iter()
                .any(|s| s.contains("threshold"))
        );
        assert!(temp.path().join("study/workflow.json").is_file());
        area.gene_assay_study_ui.request_json.push(' ');
        assert!(
            !area
                .gene_assay_study_ui
                .normalized_is_current(area.study_baseline())
        );
    }

    #[test]
    fn study_tampered_workflow_fails_without_attaching_design() {
        let temp = tempfile::tempdir().unwrap();
        let dir = temp.path().join("study");
        let mut area = planned_area(&dir);
        let workflow = std::fs::read_to_string(dir.join("workflow.json")).unwrap();
        std::fs::write(dir.join("workflow.json"), format!("{workflow}\n")).unwrap();
        area.gene_assay_study_ui.execution_reviewed = true;
        area.execute_reviewed_study();
        drain(&mut area);
        assert!(
            area.gene_assay_study_ui.status.contains("digest mismatch"),
            "{}",
            area.gene_assay_study_ui.status
        );
        assert!(
            area.engine
                .as_ref()
                .unwrap()
                .read()
                .unwrap()
                .list_transcript_assay_panel_reports()
                .is_empty()
        );
    }

    #[test]
    fn study_plan_does_not_overwrite_previous_output() {
        let temp = tempfile::tempdir().unwrap();
        let dir = temp.path().join("study");
        let mut area = planned_area(&dir);
        let before = std::fs::read(dir.join("plan.json")).unwrap();
        area.plan_reviewed_study();
        drain(&mut area);
        assert!(
            area.gene_assay_study_ui
                .status
                .contains("never overwritten")
        );
        assert_eq!(before, std::fs::read(dir.join("plan.json")).unwrap());
        assert!(!area.gene_assay_study_ui.executable_plan);
    }

    #[test]
    fn study_execution_refuses_existing_panel_identity() {
        let temp = tempfile::tempdir().unwrap();
        let mut area = planned_area(&temp.path().join("study"));
        let plan = area.gene_assay_study_ui.plan.clone().unwrap();
        let id = plan.planned_operations[0].operation["DesignTranscriptAssayPanel"]["report_id"]
            .as_str()
            .unwrap();
        let workflow: serde_json::Value = serde_json::from_str(include_str!(
            "../../docs/examples/workflows/patz1_endpoint_sybr_transcript_assay_panel_offline.json"
        ))
        .unwrap();
        let mut operation = workflow["workflow"]["ops"][2].clone();
        operation["DesignTranscriptAssayPanel"]["report_id"] = json!(id);
        operation["DesignTranscriptAssayPanel"]["path"] = serde_json::Value::Null;
        area.engine
            .as_ref()
            .unwrap()
            .write()
            .unwrap()
            .apply(serde_json::from_value(operation).unwrap())
            .unwrap();
        // A report added since planning invalidates the first approval already.
        area.gene_assay_study_ui.execution_reviewed = true;
        area.execute_reviewed_study();
        assert!(area.gene_assay_study_ui.job.is_none());
        // Even after re-normalizing/replanning, an existing ID is never replaced.
        area.normalize_study_draft();
        drain(&mut area);
        area.gene_assay_study_ui.output_directory =
            temp.path().join("iteration2").to_string_lossy().into();
        area.gene_assay_study_ui.planning_reviewed = true;
        area.plan_reviewed_study();
        drain(&mut area);
        area.gene_assay_study_ui.execution_reviewed = true;
        area.execute_reviewed_study();
        drain(&mut area);
        assert!(
            area.gene_assay_study_ui.status.contains("already exists"),
            "{}",
            area.gene_assay_study_ui.status
        );
        assert_eq!(
            area.engine
                .as_ref()
                .unwrap()
                .read()
                .unwrap()
                .list_transcript_assay_panel_reports()
                .len(),
            1
        );
    }

    #[test]
    fn study_dossier_matches_shared_export_and_remains_pending() {
        let temp = tempfile::tempdir().unwrap();
        let mut area = planned_area(&temp.path().join("study"));
        let plan_path = temp.path().join("study/plan.json");
        let request_path = temp.path().join("publication.json");
        std::fs::write(&request_path, serde_json::to_vec(&json!({
            "schema": "gentle.gene_isoform_assay_publication_request.v1",
            "report_id": "study_gui_test", "title": "Synthetic planning only",
            "genes": [{"gene_symbol": "PATZ1", "status": "pending", "status_reason": "Not executed; no specificity or order approval",
                "study_plan": {"path": "study/plan.json", "expected_sha256": sha256_prefixed_bytes(&std::fs::read(&plan_path).unwrap())}}]
        })).unwrap()).unwrap();
        area.gene_assay_study_ui.publication_request = request_path.to_string_lossy().into();
        let output = temp.path().join("gui-dossier");
        area.gene_assay_study_ui.publication_directory = output.to_string_lossy().into();
        area.publish_study();
        drain(&mut area);
        assert!(
            area.gene_assay_study_ui
                .status
                .contains("Publish completed"),
            "{}",
            area.gene_assay_study_ui.status
        );
        let canonical = std::fs::read(output.join("canonical-report.json")).unwrap();
        let reference = temp.path().join("shared-dossier");
        run_command(
            &mut area.engine.as_ref().unwrap().write().unwrap(),
            &ShellCommand::PrimersPublishGeneIsoformAssayStudy {
                request_path: request_path.to_string_lossy().into(),
                output_directory: reference.to_string_lossy().into(),
                profile: None,
                block_ids: vec![],
                generate_pdf: false,
            },
            Arc::new(Mutex::new(Box::new(|_| true))),
        )
        .unwrap();
        let actual: serde_json::Value = serde_json::from_slice(&canonical).unwrap();
        let expected: serde_json::Value = serde_json::from_slice(
            &std::fs::read(reference.join("canonical-report.json")).unwrap(),
        )
        .unwrap();
        assert_eq!(actual["genes"], expected["genes"]);
        assert_eq!(actual["complete"], false);
        assert_eq!(actual["pending_gene_count"], 1);
        assert_eq!(actual["genes"][0]["handoffs"], json!([]));
        area.publish_study();
        drain(&mut area);
        assert!(
            area.gene_assay_study_ui
                .status
                .contains("new dossier output directory")
        );
        assert_eq!(
            std::fs::read(output.join("canonical-report.json")).unwrap(),
            canonical
        );
    }

    #[test]
    fn saved_plan_rejects_foreign_locus_and_non_plan_payloads() {
        let plan = GeneIsoformAssayStudyPlanReport {
            schema: "gentle.gene_isoform_assay_study_plan.v1".into(),
            plan_id: "p".into(),
            seq_id: "foreign".into(),
            ..Default::default()
        };
        assert!(
            parse_saved_plan(&serde_json::to_vec(&plan).unwrap(), "local")
                .unwrap_err()
                .contains("active sequence")
        );
        assert!(parse_saved_plan(b"{}", "").is_err());
    }

    #[test]
    fn study_result_does_not_attach_to_replaced_project_or_edited_draft() {
        let mut area = area();
        let (entered, wait) = mpsc::channel();
        let (release, resume) = mpsc::channel();
        area.start_study_work(StudyJobKind::Normalize, "synthetic".into(), move |_, _| {
            entered.send(()).unwrap();
            resume.recv().unwrap();
            value_result(json!(GeneIsoformAssayStudyPlanRequest::default()))
        });
        wait.recv_timeout(Duration::from_secs(3)).unwrap();
        area.gene_assay_study_ui.request_json = "edited while running".into();
        *area.engine.as_ref().unwrap().write().unwrap() =
            GentleEngine::from_state(ProjectState::default());
        release.send(()).unwrap();
        drain(&mut area);
        assert!(area.gene_assay_study_ui.normalized.is_none());
        assert!(area.gene_assay_study_ui.status.contains("not attached"));
    }

    #[test]
    fn study_cancelled_work_never_grants_review() {
        let mut area = area();
        let (entered, wait) = mpsc::channel();
        let (release, resume) = mpsc::channel();
        area.start_study_work(StudyJobKind::Normalize, "synthetic".into(), move |_, _| {
            entered.send(()).unwrap();
            resume.recv().unwrap();
            value_result(json!(GeneIsoformAssayStudyPlanRequest::default()))
        });
        wait.recv_timeout(Duration::from_secs(3)).unwrap();
        assert!(
            area.gene_assay_study_ui
                .service
                .cancel(area.gene_assay_study_ui.job.as_ref().unwrap().id)
        );
        release.send(()).unwrap();
        drain(&mut area);
        assert!(area.gene_assay_study_ui.normalized.is_none());
        assert!(area.gene_assay_study_ui.status.contains("cancelled"));
    }

    #[test]
    fn study_result_rejects_edits_after_worker_commit_before_ui_poll() {
        let mut area = area();
        area.start_study_work(StudyJobKind::Execute, "synthetic".into(), |_, _| {
            value_result(json!({}))
        });
        let job = area.gene_assay_study_ui.job.as_ref().unwrap().id;
        let deadline = Instant::now() + Duration::from_secs(10);
        loop {
            let receipt = area.gene_assay_study_ui.service.status(job).unwrap();
            if receipt.result_structural_revision.is_some() {
                break;
            }
            assert!(
                Instant::now() < deadline,
                "worker did not commit: {receipt:?}"
            );
            std::thread::sleep(Duration::from_millis(5));
        }
        area.engine
            .as_ref()
            .unwrap()
            .write()
            .unwrap()
            .apply(Operation::CreateSequenceFromText {
                sequence_text: "ATGC".into(),
                output_id: Some("later_edit".into()),
                name: None,
                circular: false,
            })
            .unwrap();
        drain(&mut area);
        assert!(area.gene_assay_study_ui.status.contains("not attached"));
        assert!(!area.gene_assay_study_ui.execution_reviewed);
    }

    #[test]
    fn selected_pair_keeps_minus_locus_cdna_footprints_and_missing_checks() {
        let report: TranscriptAssayPanelReport = serde_json::from_str(include_str!("../../docs/tutorial/generated/artifacts/patz1_transcript_assay_panels_cli/artifacts/patz1_sybr_juc_panel.report.json")).unwrap();
        assert_eq!(report.strand, "-");
        let mut selected = report.selected_assays.last().unwrap().assay_id.clone();
        let before = serde_json::to_value(&report).unwrap();
        let ctx = egui::Context::default();
        ctx.begin_pass(egui::RawInput {
            screen_rect: Some(egui::Rect::from_min_size(
                egui::Pos2::ZERO,
                egui::vec2(1600.0, 2000.0),
            )),
            ..Default::default()
        });
        crate::egui_compat::show_central_panel_for_test_context(
            &ctx,
            egui::CentralPanel::default(),
            |ui| MainAreaDna::render_selected_transcript_pair(ui, &report, &mut selected),
        );
        let output = crate::egui_compat::end_test_pass(&ctx);
        fn text(shape: &egui::epaint::Shape) -> String {
            match shape {
                egui::epaint::Shape::Text(t) => t.galley.job.text.clone(),
                egui::epaint::Shape::Vec(shapes) => {
                    shapes.iter().map(text).collect::<Vec<_>>().join("\n")
                }
                _ => String::new(),
            }
        }
        let labels = output
            .shapes
            .iter()
            .map(|s| text(&s.shape))
            .collect::<Vec<_>>()
            .join("\n");
        let pair = &report.selected_assays.last().unwrap().primer_pair;
        assert!(
            labels.contains(&format!(
                "Forward: {} | [{}, {})",
                pair.forward.sequence, pair.forward.start_0based, pair.forward.end_0based_exclusive
            )),
            "{labels}"
        );
        assert!(labels.contains(&format!(
            "Reverse: {} | [{}, {})",
            pair.reverse.sequence, pair.reverse.start_0based, pair.reverse.end_0based_exclusive
        )));
        for transcript in &report.transcript_rows {
            assert!(labels.contains(&transcript.transcript_id));
        }
        assert_eq!(selected, report.selected_assays.last().unwrap().assay_id);
        assert_eq!(before, serde_json::to_value(&report).unwrap());
        assert!(report.specificity_acceptance.is_none());
    }
}

#[derive(Clone, Debug)]
struct StudyJob {
    id: u64,
    kind: StudyJobKind,
    seq_id: String,
    baseline: (u64, u64),
    draft: String,
}

#[derive(Clone, Default)]
pub(super) struct GeneAssayStudyUi {
    service: CommandExecutionService,
    job: Option<StudyJob>,
    request_json: String,
    normalized: Option<GeneIsoformAssayStudyPlanRequest>,
    normalized_draft: String,
    baseline: Option<(u64, u64)>,
    planning_reviewed: bool,
    execution_reviewed: bool,
    plan: Option<Arc<GeneIsoformAssayStudyPlanReport>>,
    executable_plan: bool,
    output_directory: String,
    planned_directory: Option<PathBuf>,
    publication_request: String,
    publication_directory: String,
    publication_pdf: bool,
    pub(super) selected_pair: String,
    status: String,
    last_receipt: Option<crate::command_execution::CommandReceipt>,
    panel_summary_attempt: Option<(u64, u64)>,
    panel_summaries: Option<(
        (u64, u64),
        Vec<crate::engine::TranscriptAssayPanelReportSummary>,
    )>,
}

impl std::fmt::Debug for GeneAssayStudyUi {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("GeneAssayStudyUi")
            .field("job", &self.job)
            .finish_non_exhaustive()
    }
}

impl GeneAssayStudyUi {
    fn invalidate_review(&mut self) {
        self.normalized = None;
        self.planning_reviewed = false;
        self.execution_reviewed = false;
        self.executable_plan = false;
    }

    fn normalized_is_current(&self, baseline: Option<(u64, u64)>) -> bool {
        self.normalized.is_some()
            && self.normalized_draft == self.request_json
            && self.baseline.is_some()
            && self.baseline == baseline
    }
}

fn parse_saved_plan(bytes: &[u8], seq_id: &str) -> Result<GeneIsoformAssayStudyPlanReport, String> {
    let plan: GeneIsoformAssayStudyPlanReport =
        serde_json::from_slice(bytes).map_err(|e| format!("Not a study-plan JSON: {e}"))?;
    if plan.schema != "gentle.gene_isoform_assay_study_plan.v1" || plan.plan_id.is_empty() {
        return Err(
            "Select a typed gene isoform assay study plan, not a PDF or publication receipt".into(),
        );
    }
    if plan.seq_id != seq_id {
        return Err(format!(
            "This plan belongs to '{}', not the active sequence '{seq_id}'. Open its source sequence first.",
            plan.seq_id
        ));
    }
    Ok(plan)
}

fn run_command(
    engine: &mut GentleEngine,
    command: &ShellCommand,
    progress: ShellProgressCallback,
) -> Result<ShellRunResult, String> {
    execute_shell_command_with_options(
        engine,
        command,
        &ShellExecutionOptions {
            allow_agent_commands: false,
            progress_callback: Some(progress),
            ..Default::default()
        },
    )
}

fn value_result(value: serde_json::Value) -> Result<ShellRunResult, String> {
    Ok(ShellRunResult {
        state_changed: false,
        output: value,
    })
}

fn study_button(
    ui: &mut egui::Ui,
    enabled: bool,
    label: &str,
    _id: &'static str,
    _seq_id: &str,
) -> bool {
    let response = ui.add_enabled(enabled, egui::Button::new(label));
    #[cfg(feature = "gui-test-support")]
    crate::gui_test_support::register_response(
        &response,
        _id,
        "window.pcr_design",
        Some(&crate::gui_test_support::pseudonymous_subject_scope(&[
            _seq_id,
        ])),
        crate::gui_test_support::GuiTestWidgetKind::Button,
        false,
    );
    response.clicked()
}

fn show_json(ui: &mut egui::Ui, title: &str, value: &impl serde::Serialize) {
    egui::CollapsingHeader::new(title).show(ui, |ui| {
        if let Ok(text) = serde_json::to_string_pretty(value) {
            ui.monospace(text);
        }
    });
}

impl MainAreaDna {
    fn study_baseline(&self) -> Option<(u64, u64)> {
        let engine = self.engine.as_ref()?.try_read().ok()?;
        Some((engine.instance_id(), engine.structural_revision()))
    }

    fn start_study_work<F>(&mut self, kind: StudyJobKind, input_binding: String, work: F)
    where
        F: FnOnce(&mut GentleEngine, ShellProgressCallback) -> Result<ShellRunResult, String>
            + Send
            + 'static,
    {
        let Some(engine) = self.engine.clone() else {
            return;
        };
        let Some(baseline) = self.study_baseline() else {
            self.gene_assay_study_ui.status = "Project busy; retry shortly".into();
            return;
        };
        if self.gene_assay_study_ui.job.is_some() {
            return;
        }
        let draft = self.gene_assay_study_ui.request_json.clone();
        let digest =
            sha256_prefixed_bytes(format!("gene-study:{kind:?}:{input_binding}").as_bytes());
        match self
            .gene_assay_study_ui
            .service
            .submit_work(engine, digest, work)
        {
            Ok(id) => {
                self.gene_assay_study_ui.job = Some(StudyJob {
                    id,
                    kind,
                    seq_id: self.seq_id.clone().unwrap_or_default(),
                    baseline,
                    draft,
                });
                self.gene_assay_study_ui.status =
                    format!("{kind:?} admitted; waiting for completion");
            }
            Err(error) => self.gene_assay_study_ui.status = error,
        }
    }

    pub(super) fn poll_gene_assay_study_task(&mut self, ctx: &egui::Context) {
        let Some(job) = self.gene_assay_study_ui.job.clone() else {
            return;
        };
        let Some(engine) = self.engine.clone() else {
            return;
        };
        let Ok(live) = engine.try_read() else {
            ctx.request_repaint_after(Duration::from_millis(160));
            return;
        };
        // Observe the receipt and project revision under the same read lock;
        // a worker commit between the two observations must not look stale.
        let current = (live.instance_id(), live.structural_revision());
        let Some(result) = self.gene_assay_study_ui.service.take_result(job.id) else {
            ctx.request_repaint_after(Duration::from_millis(160));
            return;
        };
        self.gene_assay_study_ui.job = None;
        self.gene_assay_study_ui.last_receipt = self.gene_assay_study_ui.service.status(job.id);
        let current = Some(current);
        let committed = self
            .gene_assay_study_ui
            .last_receipt
            .as_ref()
            .and_then(|r| Some((r.result_instance?, r.result_structural_revision?)));
        let expected = committed.unwrap_or(job.baseline);
        if self.seq_id.as_deref() != Some(job.seq_id.as_str())
            || current.map(|b| b.0) != Some(job.baseline.0)
            || current != Some(expected)
            || self.gene_assay_study_ui.request_json != job.draft
        {
            self.gene_assay_study_ui.invalidate_review();
            self.gene_assay_study_ui.status = "Study result not attached: project, sequence or request changed. Files already written are not rolled back.".into();
            return;
        }
        let result = match result {
            Ok(result) => result,
            Err(error) => {
                self.gene_assay_study_ui.execution_reviewed = false;
                self.gene_assay_study_ui.status = error;
                return;
            }
        };
        let accepted = (|| -> Result<(), String> {
            match job.kind {
                StudyJobKind::ReadRequest => {
                    let request: GeneIsoformAssayStudyPlanRequest =
                        serde_json::from_value(result.output).map_err(|e| e.to_string())?;
                    self.gene_assay_study_ui.request_json =
                        serde_json::to_string_pretty(&request).map_err(|e| e.to_string())?;
                    self.gene_assay_study_ui.invalidate_review();
                }
                StudyJobKind::Normalize => {
                    let request =
                        serde_json::from_value(result.output).map_err(|e| e.to_string())?;
                    self.gene_assay_study_ui.normalized = Some(request);
                    self.gene_assay_study_ui.normalized_draft = job.draft;
                    self.gene_assay_study_ui.baseline = current;
                    self.gene_assay_study_ui.planning_reviewed = false;
                }
                StudyJobKind::ReadPlan | StudyJobKind::Plan => {
                    let value = if job.kind == StudyJobKind::Plan {
                        result
                            .output
                            .get("report")
                            .cloned()
                            .ok_or("Planner returned no report")?
                    } else {
                        result.output
                    };
                    let plan = parse_saved_plan(
                        &serde_json::to_vec(&value).map_err(|e| e.to_string())?,
                        &job.seq_id,
                    )?;
                    self.gene_assay_study_ui.plan = Some(Arc::new(plan));
                    self.gene_assay_study_ui.execution_reviewed = false;
                    self.gene_assay_study_ui.executable_plan = job.kind == StudyJobKind::Plan;
                    self.gene_assay_study_ui.baseline = current;
                }
                StudyJobKind::Execute => {
                    self.gene_assay_study_ui.invalidate_review();
                    self.cached_transcript_assay_panel_report = None;
                    self.cached_experimental_assay_handoff = None;
                }
                StudyJobKind::PanelSummaries => {
                    let summaries =
                        serde_json::from_value(result.output).map_err(|e| e.to_string())?;
                    self.gene_assay_study_ui.panel_summaries = current.map(|key| (key, summaries));
                }
                StudyJobKind::Publish => {}
            }
            Ok(())
        })();
        self.gene_assay_study_ui.status = match accepted {
            Ok(()) => format!(
                "{:?} completed. This is not experimental validation or order approval.",
                job.kind
            ),
            Err(error) => {
                self.gene_assay_study_ui.invalidate_review();
                error
            }
        };
    }

    fn load_study_file(&mut self, kind: StudyJobKind) {
        let Some(path) = rfd::FileDialog::new()
            .add_filter("JSON", &["json"])
            .pick_file()
        else {
            return;
        };
        let seq_id = self.seq_id.clone().unwrap_or_default();
        self.gene_assay_study_ui.invalidate_review();
        self.start_study_work(kind, path.to_string_lossy().into(), move |_, _| {
            let bytes = std::fs::read(&path).map_err(|e| e.to_string())?;
            let value = if kind == StudyJobKind::ReadPlan {
                serde_json::to_value(parse_saved_plan(&bytes, &seq_id)?)
                    .map_err(|e| e.to_string())?
            } else {
                let request: GeneIsoformAssayStudyPlanRequest =
                    serde_json::from_slice(&bytes).map_err(|e| e.to_string())?;
                if request.schema != "gentle.gene_isoform_assay_study_plan_request.v1" {
                    return Err("Select a study request, not a saved report".into());
                }
                serde_json::to_value(request).map_err(|e| e.to_string())?
            };
            value_result(value)
        });
    }

    fn normalize_study_draft(&mut self) {
        self.gene_assay_study_ui.invalidate_review();
        let draft = self.gene_assay_study_ui.request_json.clone();
        self.start_study_work(
            StudyJobKind::Normalize,
            draft.clone(),
            move |engine, progress| {
                run_command(
                    engine,
                    &ShellCommand::PrimersPlanGeneIsoformAssayStudy {
                        request_json: draft,
                        normalize_only: true,
                        normalized_request_path: None,
                        path: None,
                        workflow_path: None,
                    },
                    progress,
                )
            },
        );
    }

    fn plan_reviewed_study(&mut self) {
        let state = &self.gene_assay_study_ui;
        if !state.planning_reviewed || !state.normalized_is_current(self.study_baseline()) {
            return;
        }
        let request = state.normalized.clone().unwrap();
        let directory = PathBuf::from(state.output_directory.trim());
        if !directory.is_absolute() {
            self.gene_assay_study_ui.status =
                "Choose a new absolute output directory, not an existing study directory".into();
            return;
        }
        self.gene_assay_study_ui.planned_directory = Some(directory.clone());
        self.gene_assay_study_ui.executable_plan = false;
        let binding = format!(
            "{}|{}",
            serde_json::to_string(&request).unwrap(),
            directory.display()
        );
        self.start_study_work(StudyJobKind::Plan, binding, move |engine, progress| {
            if engine
                .normalize_gene_isoform_assay_study_request(request.clone())
                .map_err(|e| e.to_string())?
                != request
            {
                return Err(
                    "Normalized input changed after review; normalize and review again".into(),
                );
            }
            std::fs::create_dir(&directory).map_err(|e| {
                format!("Use a new directory; previous studies are never overwritten: {e}")
            })?;
            run_command(
                engine,
                &ShellCommand::PrimersPlanGeneIsoformAssayStudy {
                    request_json: serde_json::to_string(&request).map_err(|e| e.to_string())?,
                    normalize_only: false,
                    normalized_request_path: Some(
                        directory.join("request.json").to_string_lossy().into(),
                    ),
                    path: Some(directory.join("plan.json").to_string_lossy().into()),
                    workflow_path: Some(directory.join("workflow.json").to_string_lossy().into()),
                },
                progress,
            )
        });
    }

    fn execute_reviewed_study(&mut self) {
        let state = &self.gene_assay_study_ui;
        if !state.execution_reviewed
            || !state.executable_plan
            || !state.normalized_is_current(self.study_baseline())
        {
            return;
        }
        let Some(plan) = state.plan.clone() else {
            return;
        };
        let Some(directory) = state.planned_directory.clone() else {
            return;
        };
        let plan_path = directory.join("plan.json");
        let expected = match serde_json::to_vec_pretty(plan.as_ref()) {
            Ok(bytes) => sha256_prefixed_bytes(&bytes),
            Err(error) => {
                self.gene_assay_study_ui.status = error.to_string();
                return;
            }
        };
        self.gene_assay_study_ui.execution_reviewed = false;
        self.gene_assay_study_ui.executable_plan = false;
        let binding = format!(
            "{}|{}|{}",
            expected,
            plan.approved_workflow_sha256,
            directory.display()
        );
        self.start_study_work(StudyJobKind::Execute, binding, move |engine, progress| {
            let bytes = std::fs::read(&plan_path).map_err(|e| e.to_string())?;
            if sha256_prefixed_bytes(&bytes) != expected {
                return Err("Plan bytes changed after review; no design executed".into());
            }
            let normalized = engine.normalize_gene_isoform_assay_study_request(plan.normalized_request.clone()).map_err(|e| e.to_string())?;
            if normalized != plan.normalized_request {
                return Err("Evidence or effective request changed after review; no design executed".into());
            }
            for row in &plan.planned_operations {
                let op = &row.operation;
                let spec = op.get("DesignTranscriptAssayPanel")
                    .or_else(|| op.pointer("/DesignTranscriptAssayPanelWithFallback/strict_operation/DesignTranscriptAssayPanel"));
                if let Some(id) = spec.and_then(|v| v.get("report_id")).and_then(|v| v.as_str())
                    && engine.get_transcript_assay_panel_report(id).is_ok() {
                    return Err(format!("Panel '{id}' already exists; choose a new plan_id rather than overwrite results"));
                }
                if let Some(id) = op.pointer("/DesignTranscriptAssayPanelWithFallback/fallback_submission/fallback_report_id").and_then(|v| v.as_str())
                    && engine.get_transcript_assay_panel_report(id).is_ok() {
                    return Err(format!("Fallback panel '{id}' already exists; use a new identity"));
                }
            }
            run_command(engine, &ShellCommand::PrimersExecuteGeneIsoformAssayStudyWorkflow {
                plan_json: String::from_utf8(bytes).map_err(|e| e.to_string())?,
                workflow_json: format!("@{}", directory.join("workflow.json").display()),
            }, progress)
        });
    }

    pub(super) fn render_gene_assay_study(&mut self, ui: &mut egui::Ui) {
        self.render_pcr_designer_mode_selector(ui);
        let current = self.study_baseline();
        // Store deserialization can be large; do it once per project revision
        // on the worker, never on each egui repaint.
        if self.gene_assay_study_ui.job.is_none()
            && self.gene_assay_study_ui.plan.is_some()
            && current.is_some()
            && self.gene_assay_study_ui.panel_summary_attempt != current
            && self
                .gene_assay_study_ui
                .panel_summaries
                .as_ref()
                .map(|(key, _)| *key)
                != current
        {
            self.gene_assay_study_ui.panel_summary_attempt = current;
            self.gene_assay_study_ui.panel_summaries = None;
            self.cached_transcript_assay_panel_report = None;
            self.start_study_work(
                StudyJobKind::PanelSummaries,
                format!("{current:?}"),
                |engine, _| {
                    value_result(
                        serde_json::to_value(engine.list_transcript_assay_panel_reports())
                            .map_err(|e| e.to_string())?,
                    )
                },
            );
        }
        let busy = self.gene_assay_study_ui.job.is_some();
        if self.gene_assay_study_ui.baseline.is_some()
            && self.gene_assay_study_ui.baseline != current
            && !busy
        {
            self.gene_assay_study_ui.invalidate_review();
        }
        ui.heading("Gene-informed primer-pair study");
        ui.label("Only declared evidence is assessed. Annotation, observations, predictions and missing checks remain distinct.");
        let seq_id = self.seq_id.clone().unwrap_or_default();
        ui.horizontal_wrapped(|ui| {
            if study_button(
                ui,
                !busy,
                "Open study plan...",
                "assay_study.open_plan",
                &seq_id,
            ) {
                self.load_study_file(StudyJobKind::ReadPlan);
            }
            if study_button(
                ui,
                !busy,
                "Open planning request...",
                "assay_study.open_request",
                &seq_id,
            ) {
                self.load_study_file(StudyJobKind::ReadRequest);
            }
        });
        if let Some(job) = &self.gene_assay_study_ui.job {
            if let Some(receipt) = self.gene_assay_study_ui.service.status(job.id) {
                ui.label(format!(
                    "{:?}: {} | {:?}/{:?} steps computed (not yet accepted)",
                    job.kind, receipt.phase, receipt.completed_steps, receipt.total_steps
                ));
                if ui.button("Cancel study task").clicked() {
                    self.gene_assay_study_ui.service.cancel(job.id);
                }
            }
        }
        ui.label(&self.gene_assay_study_ui.status);
        if let Some(receipt) = &self.gene_assay_study_ui.last_receipt {
            show_json(
                ui,
                "Last task receipt (execution, not scientific acceptance)",
                receipt,
            );
        }
        egui::ScrollArea::vertical().id_salt("gene_assay_study_scroll").show(ui, |ui| {
            if let Some(plan) = self.gene_assay_study_ui.plan.clone() {
                Self::render_study_plan(ui, &plan);
                if ui.button("Inspect transcript architecture in Splicing Expert").clicked() {
                    self.open_splicing_expert_for_feature(plan.source_feature_id, "gene assay study");
                }
                if !self.gene_assay_study_ui.executable_plan {
                    ui.small("Saved/historical plan inspection is read-only. Load its normalized request for a new reviewed iteration; inspection does not authorize execution.");
                }
                self.render_study_panels(ui, &plan);
            }
            egui::CollapsingHeader::new("Design and review a new iteration").default_open(self.gene_assay_study_ui.plan.is_none()).show(ui, |ui| {
                ui.small("Load a typed request with an explicit evidence path. Paths resolve as in the shared Shell, not relative to this window. Use a new plan_id for each iteration.");
                let mut changed = false;
                ui.add_enabled_ui(!busy, |ui| {
                    if let Ok(mut request) = serde_json::from_str::<GeneIsoformAssayStudyPlanRequest>(&self.gene_assay_study_ui.request_json) {
                        ui.horizontal(|ui| {
                            ui.label("Short-product maximum (bp; hard requirement)");
                            if ui.add(egui::DragValue::new(&mut request.policy.short_max_amplicon_bp).range(1..=100_000)).changed() {
                                self.gene_assay_study_ui.request_json = serde_json::to_string_pretty(&request).unwrap();
                                changed = true;
                            }
                        });
                    }
                    changed |= ui.add(egui::TextEdit::multiline(&mut self.gene_assay_study_ui.request_json).desired_rows(8).desired_width(f32::INFINITY).code_editor()).changed();
                });
                if changed { self.gene_assay_study_ui.invalidate_review(); }
                if study_button(ui, !busy && !self.gene_assay_study_ui.request_json.trim().is_empty(), "Normalize request", "assay_study.normalize", &seq_id) {
                    self.normalize_study_draft();
                }
                if let Some(normalized) = &self.gene_assay_study_ui.normalized {
                    show_json(ui, "Effective normalized request (review before planning)", normalized);
                }
                let ready = !busy && self.gene_assay_study_ui.normalized_is_current(current);
                ui.add_enabled_ui(ready, |ui| {
                    ui.checkbox(&mut self.gene_assay_study_ui.planning_reviewed, "I reviewed this effective request and its bound evidence (planning only)");
                    ui.label("New absolute output directory (parent must exist)");
                    ui.text_edit_singleline(&mut self.gene_assay_study_ui.output_directory);
                });
                if study_button(ui, ready && self.gene_assay_study_ui.planning_reviewed, "Plan reviewed study", "assay_study.plan", &seq_id) { self.plan_reviewed_study(); }
                ui.add_enabled_ui(ready && self.gene_assay_study_ui.executable_plan, |ui| {
                    ui.checkbox(&mut self.gene_assay_study_ui.execution_reviewed, "I separately reviewed the exact ordered operations and workflow digest (execute design)");
                });
                if study_button(ui, ready && self.gene_assay_study_ui.executable_plan && self.gene_assay_study_ui.execution_reviewed, "Execute reviewed workflow", "assay_study.execute", &seq_id) { self.execute_reviewed_study(); }
                ui.small("Cancellation or stale results discard engine changes, not files already written. Prior output directories and panel identities are never intentionally reused.");
            });
            self.render_study_publication(ui, busy);
        });
    }

    fn render_study_plan(ui: &mut egui::Ui, plan: &GeneIsoformAssayStudyPlanReport) {
        ui.heading(&plan.label);
        ui.label(format!(
            "{} | {} | source {} / feature n-{} | annotation {}",
            plan.gene_symbol,
            plan.plan_id,
            plan.seq_id,
            plan.source_feature_id + 1,
            plan.annotation_release.as_deref().unwrap_or("not supplied")
        ));
        ui.label(format!(
            "{} transcript records / {} exact cDNA groups; scope {:?}",
            plan.evidence_summary.transcript_count,
            plan.evidence_summary.exact_cdna_equivalence_group_count,
            plan.coverage_resolution.universe.kind
        ));
        ui.label(format!(
            "Automatic recommendation: {}. Selected: {}.",
            plan.recommended_profile.as_str(),
            plan.selected_profile.as_str()
        ));
        if let Some(value) = &plan.profile_override {
            ui.label(format!("Explicit override: {}", value.reason));
        }
        ui.small("Reference/assembly identity is not inferred from the gene symbol. Inspect source evidence and the panel's stored genome anchor below.");
        show_json(
            ui,
            "Declared evidence inventory and input digests",
            &plan.resolved_evidence_inputs,
        );
        show_json(
            ui,
            "Evidence summary (assessed inputs only)",
            &plan.evidence_summary,
        );
        show_json(
            ui,
            "Transcript coverage scope and exclusions",
            &plan.coverage_resolution,
        );
        for missing in &plan.evidence_summary.missing_evidence {
            ui.label(format!("Missing: {missing}"));
        }
        for factor in &plan.decision_factors {
            ui.label(format!(
                "{} [{}]: {}",
                factor.rule_id, factor.triggered, factor.summary
            ));
        }
        for warning in &plan.warnings {
            ui.colored_label(egui::Color32::DARK_RED, warning);
        }
        ui.monospace(format!(
            "Request: {}\nOperations: {}\nWorkflow bytes: {}",
            plan.request_sha256, plan.operation_batch_sha256, plan.approved_workflow_sha256
        ));
        show_json(
            ui,
            "Exact ordered operations (second review)",
            &plan.planned_operations,
        );
    }

    fn render_study_panels(&mut self, ui: &mut egui::Ui, plan: &GeneIsoformAssayStudyPlanReport) {
        if ui
            .add_enabled(
                self.gene_assay_study_ui.job.is_none(),
                egui::Button::new("Refresh panel list"),
            )
            .clicked()
        {
            self.gene_assay_study_ui.panel_summary_attempt = None;
            self.gene_assay_study_ui.panel_summaries = None;
        }
        let Some((_, reports)) = self.gene_assay_study_ui.panel_summaries.clone() else {
            ui.label(if self.gene_assay_study_ui.job.is_some() {
                "Loading persisted panel summaries in the background..."
            } else {
                "Panel summaries unavailable; refresh to retry. This does not mean that no panels exist."
            });
            return;
        };
        let matching = reports
            .iter()
            .filter(|r| r.source_seq_id == plan.seq_id)
            .collect::<Vec<_>>();
        ui.separator();
        ui.label(
            "Persisted panels on this source sequence (not automatically members of this plan)",
        );
        if matching.is_empty() {
            ui.label("No matching panels in this project. A study plan alone contains no designed primer pairs.");
        }
        for report in matching {
            if ui
                .button(format!(
                    "Inspect {} ({} pairs, feature n-{})",
                    report.report_id,
                    report.selected_assay_count,
                    report.source_feature_id + 1
                ))
                .clicked()
            {
                self.show_transcript_assay_panel_report(&report.report_id);
                self.pcr_designer_mode = PcrDesignerMode::GeneAssayStudy;
            }
        }
        if let Some(report) = self.cached_transcript_assay_panel_report.clone()
            && report.source_seq_id == plan.seq_id
        {
            if report.source_feature_id != plan.source_feature_id {
                ui.label("Different source feature: this is comparison context, not proof that this panel belongs to the study. Canonical publication revalidates the exact plan/handoff binding.");
            }
            Self::render_selected_transcript_pair(
                ui,
                &report,
                &mut self.gene_assay_study_ui.selected_pair,
            );
            if ui
                .button("Open panel design, matrix and readiness checks")
                .clicked()
            {
                self.pcr_designer_mode = PcrDesignerMode::TranscriptPanels;
            }
        }
    }

    fn render_study_publication(&mut self, ui: &mut egui::Ui, busy: bool) {
        egui::CollapsingHeader::new("Canonical dossier export").show(ui, |ui| {
            ui.small("Use the same gentle.gene_isoform_assay_publication_request.v1 as CLI/OpenClaw. It names exact plan/handoff/order hashes. No GUI-authored scientific narrative and no ordering.");
            ui.add_enabled_ui(!busy, |ui| {
                ui.label("Publication request JSON (absolute path)");
                ui.text_edit_singleline(&mut self.gene_assay_study_ui.publication_request);
                ui.label("New absolute output directory");
                ui.text_edit_singleline(&mut self.gene_assay_study_ui.publication_directory);
                ui.checkbox(&mut self.gene_assay_study_ui.publication_pdf, "Also PDF (requires installed browser)");
            });
            if study_button(ui, !busy, "Export declared dossier", "assay_study.publish", self.seq_id.as_deref().unwrap_or("unnamed")) {
                self.publish_study();
            }
        });
    }

    fn publish_study(&mut self) {
        let request_path = self
            .gene_assay_study_ui
            .publication_request
            .trim()
            .to_string();
        let output_directory = self
            .gene_assay_study_ui
            .publication_directory
            .trim()
            .to_string();
        let generate_pdf = self.gene_assay_study_ui.publication_pdf;
        let binding = format!("{request_path}|{output_directory}|pdf={generate_pdf}");
        self.start_study_work(StudyJobKind::Publish, binding, move |engine, progress| {
            if !PathBuf::from(&request_path).is_absolute()
                || !PathBuf::from(&output_directory).is_absolute()
            {
                return Err("Supply absolute request/output paths".into());
            }
            std::fs::create_dir(&output_directory)
                .map_err(|e| format!("Use a new dossier output directory: {e}"))?;
            run_command(
                engine,
                &ShellCommand::PrimersPublishGeneIsoformAssayStudy {
                    request_path,
                    output_directory,
                    profile: None,
                    block_ids: vec![],
                    generate_pdf,
                },
                progress,
            )
        });
    }

    pub(super) fn render_selected_transcript_pair(
        ui: &mut egui::Ui,
        report: &TranscriptAssayPanelReport,
        selected: &mut String,
    ) {
        ui.heading("Selected primer pair");
        if !report
            .selected_assays
            .iter()
            .any(|a| a.assay_id == *selected)
        {
            *selected = report
                .selected_assays
                .first()
                .map(|a| a.assay_id.clone())
                .unwrap_or_default();
        }
        egui::ComboBox::from_id_salt("study_selected_pair")
            .selected_text(selected.as_str())
            .show_ui(ui, |ui| {
                for assay in &report.selected_assays {
                    ui.selectable_value(
                        selected,
                        assay.assay_id.clone(),
                        format!(
                            "A{} {}",
                            assay.rank, assay.primer_pair_summary.display_label
                        ),
                    );
                }
            });
        let Some(assay) = report
            .selected_assays
            .iter()
            .find(|a| a.assay_id == *selected)
        else {
            ui.label("No selected pair");
            return;
        };
        let summary = &assay.primer_pair_summary;
        ui.label(format!(
            "Design transcript: {} | locus strand {}",
            assay.design_transcript_id, report.strand
        ));
        ui.small("Primer footprints below are mature-cDNA 0-based half-open coordinates, not genomic coordinates. Both oligo sequences are written 5' to 3'. Tm is not annealing temperature.");
        for (role, primer) in [
            ("Forward", &assay.primer_pair.forward),
            ("Reverse", &assay.primer_pair.reverse),
        ] {
            ui.monospace(format!(
                "{role}: {} | [{}, {}) | Tm {:.1} C",
                primer.sequence, primer.start_0based, primer.end_0based_exclusive, primer.tm_c
            ));
        }
        ui.label(format!(
            "Designed product: {} bp; cross-transcript products are separate predictions",
            assay.primer_pair.amplicon_length_bp
        ));
        if report.specificity_acceptance.is_none() {
            ui.label("Specificity not assessed: predicted products in this panel are not a whole-reference pass or order approval.");
        } else {
            ui.label("A stored specificity assessment is available below; its scope and outcome must be reviewed separately from candidate selection.");
        }
        if let Some(group) = report
            .equivalence_groups
            .iter()
            .find(|g| g.equivalence_group_id == assay.design_equivalence_group_id)
            && group.cdna_length_bp > 0
        {
            let (rect, _) = ui.allocate_exact_size(
                egui::vec2(ui.available_width().max(120.0), 68.0),
                egui::Sense::hover(),
            );
            let x = |bp: usize| {
                rect.left()
                    + rect.width() * bp.min(group.cdna_length_bp) as f32
                        / group.cdna_length_bp as f32
            };
            let painter = ui.painter();
            painter.line_segment(
                [
                    egui::pos2(rect.left(), rect.center().y),
                    egui::pos2(rect.right(), rect.center().y),
                ],
                egui::Stroke::new(1.0, egui::Color32::GRAY),
            );
            for (primer, direction, color) in [
                (
                    &assay.primer_pair.forward,
                    1.0,
                    egui::Color32::from_rgb(0, 114, 178),
                ),
                (
                    &assay.primer_pair.reverse,
                    -1.0,
                    egui::Color32::from_rgb(213, 94, 0),
                ),
            ] {
                let left = x(primer.start_0based);
                let right = x(primer.end_0based_exclusive);
                let (from, to) = if direction > 0.0 {
                    (left, right)
                } else {
                    (right, left)
                };
                painter.arrow(
                    egui::pos2(from, rect.center().y),
                    egui::vec2(to - from, 0.0),
                    egui::Stroke::new(3.0, color),
                );
            }
            painter.text(
                rect.left_top(),
                egui::Align2::LEFT_TOP,
                "0 | mature cDNA 5'",
                egui::FontId::monospace(11.0),
                ui.visuals().text_color(),
            );
            painter.text(
                rect.right_top(),
                egui::Align2::RIGHT_TOP,
                format!("{} | 3'", group.cdna_length_bp),
                egui::FontId::monospace(11.0),
                ui.visuals().text_color(),
            );
            ui.small("Blue: forward; orange: reverse. Source-locus strand does not reverse this cDNA axis.");
        }
        for reason in &summary.selection_reasons {
            ui.label(&reason.message);
        }
        for transcript in &report.transcript_rows {
            let cell = report.detection_matrix.iter().find(|c| {
                c.assay_id == assay.assay_id
                    && c.transcript_feature_id == transcript.transcript_feature_id
            });
            let result = cell
                .map(|c| match c.status {
                    TranscriptAssayDetectionStatus::NoProduct => "no predicted product".into(),
                    _ => Self::transcript_assay_detection_label(c.status, &c.amplicon_lengths_bp),
                })
                .unwrap_or_else(|| "not assessed (matrix cell missing)".into());
            ui.label(format!(
                "{}: {} | class {}",
                transcript.transcript_id, result, transcript.equivalence_group_id
            ));
        }
        show_json(
            ui,
            "Pair selection evidence (observations, projections, missing thresholds)",
            &summary.selection_evidence,
        );
        show_json(
            ui,
            "Reference identity captured at design",
            &report.source_genome_anchor,
        );
        show_json(
            ui,
            "Specificity assessment (null means not assessed)",
            &report.specificity_acceptance,
        );
        ui.small("A predicted product does not establish isoform abundance or experimental specificity. Candidate sequences are not an order approval.");
    }
}
