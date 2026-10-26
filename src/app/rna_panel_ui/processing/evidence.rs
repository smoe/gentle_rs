//! Thin projections of shared report/row inspection and ambiguity-preserving summary export.

use super::*;
use crate::engine::{
    ExportRnaReadPanelSummaryRequest, RNA_PANEL_SUMMARY_REQUEST_SCHEMA,
    RnaReadPanelArtifactVerification, RnaReadPanelInspection, RnaReadPanelResultBinding,
    RnaReadPanelScreenReport, RnaReadPanelSummary,
};

#[derive(Default, serde::Deserialize)]
struct Content {
    rna_read_panel: Option<RnaReadPanel>,
    rna_read_panel_screen: Option<RnaReadPanelScreenReport>,
    rna_read_panel_summary: Option<RnaReadPanelSummary>,
    rna_read_panel_inspection: Option<RnaReadPanelInspection>,
    rna_read_panel_verification: Option<RnaReadPanelArtifactVerification>,
}

pub(super) struct InspectionRequest {
    path: String,
    binding: Option<RnaReadPanelResultBinding>,
    expected_content: Option<String>,
    file_id: Option<String>,
    rows: usize,
    verify: bool,
}

impl InspectionRequest {
    pub(super) fn bound(binding: RnaReadPanelResultBinding) -> Self {
        Self {
            path: binding.report_path.clone(),
            expected_content: Some(binding.report_content_sha256.clone()),
            binding: Some(binding),
            file_id: None,
            rows: 0,
            verify: false,
        }
    }
}

pub(super) struct Loaded {
    path: String,
    binding: Option<RnaReadPanelResultBinding>,
    content: Content,
    content_sha256: String,
    json: String,
}

impl Loaded {
    pub(super) fn from_result(
        request: InspectionRequest,
        mut output: serde_json::Value,
    ) -> Result<Self, String> {
        let content: Content =
            serde_json::from_value(output["result"].take()).map_err(|e| e.to_string())?;
        if usize::from(content.rna_read_panel.is_some())
            + usize::from(content.rna_read_panel_screen.is_some())
            + usize::from(content.rna_read_panel_summary.is_some())
            != 1
        {
            return Err("Shared inspection did not return exactly one scientific record".into());
        }
        let (hash, json) = if let Some(panel) = &content.rna_read_panel {
            (
                panel.content_sha256.clone(),
                serde_json::to_string_pretty(panel),
            )
        } else if let Some(report) = &content.rna_read_panel_screen {
            (
                report.content_sha256.clone(),
                serde_json::to_string_pretty(report),
            )
        } else {
            let report = content
                .rna_read_panel_summary
                .as_ref()
                .ok_or("Missing summary")?;
            (
                report.content_sha256.clone(),
                serde_json::to_string_pretty(report),
            )
        };
        if request
            .expected_content
            .as_ref()
            .is_some_and(|h| h != &hash)
        {
            return Err("Report changed since the inspected snapshot; prior view retained, refresh explicitly".into());
        }
        Ok(Self {
            path: request.path,
            binding: request.binding,
            content,
            content_sha256: hash,
            json: json.map_err(|e| e.to_string())?,
        })
    }

    fn inspect_again(
        &self,
        file_id: Option<String>,
        rows: usize,
        verify: bool,
    ) -> InspectionRequest {
        InspectionRequest {
            path: self.path.clone(),
            binding: self.binding.clone(),
            expected_content: Some(self.content_sha256.clone()),
            file_id,
            rows,
            verify,
        }
    }
}

pub(super) struct Evidence {
    pub(super) loading: bool,
    pub(super) loaded: Option<Loaded>,
    pub(super) error: Option<String>,
    pub(super) export_detail: String,
    path: String,
    binding_json: String,
    rows: usize,
    file_id: String,
    published: Vec<RnaReadPanelResultBinding>,
    request: ExportRnaReadPanelSummaryRequest,
    reviewed_export: Option<(ExportRnaReadPanelSummaryRequest, String)>,
    export_job: Option<u64>,
}

impl Default for Evidence {
    fn default() -> Self {
        Self {
            loading: false,
            loaded: None,
            error: None,
            export_detail: String::new(),
            path: String::new(),
            binding_json: String::new(),
            rows: 20,
            file_id: String::new(),
            published: vec![],
            request: ExportRnaReadPanelSummaryRequest {
                schema: RNA_PANEL_SUMMARY_REQUEST_SCHEMA.into(),
                report_paths: vec![],
                output_dir: String::new(),
            },
            reviewed_export: None,
            export_job: None,
        }
    }
}

impl Evidence {
    pub(super) fn published(&mut self, binding: RnaReadPanelResultBinding) {
        if !self.published.iter().any(|b| b == &binding) {
            if self.published.len() == 64 {
                self.published.remove(0);
            }
            self.published.push(binding);
        }
    }
    fn export_request(&self) -> Result<ExportRnaReadPanelSummaryRequest, String> {
        let (request, _) = self
            .reviewed_export
            .as_ref()
            .ok_or("Review the export request first")?;
        if serde_json::to_value(request).ok() != serde_json::to_value(&self.request).ok() {
            return Err("Export draft changed; review again".into());
        }
        Ok(request.clone())
    }
}

/// One virtualized, unframed table; no silent cap or reclassification of stored rows.
fn table(ui: &mut egui::Ui, id: &str, rows: usize, mut row: impl FnMut(&mut egui::Ui, usize)) {
    egui::ScrollArea::both()
        .id_salt(id)
        .max_height(260.0)
        .show_rows(
            ui,
            ui.text_style_height(&egui::TextStyle::Body),
            rows,
            |ui, range| {
                ui.style_mut().wrap_mode = Some(egui::TextWrapMode::Extend);
                for i in range {
                    row(ui, i);
                }
            },
        );
}

impl GENtleApp {
    pub(super) fn start_rna_inspection(
        &mut self,
        request: InspectionRequest,
    ) -> Result<(), String> {
        if self.rna_panel_ui.processing.evidence.loading {
            return Err("An inspection is already running".into());
        }
        self.submit_rna_operation(
            Operation::InspectRnaReadPanel {
                path: request.path.clone(),
                expected_result: request.binding.clone().map(Box::new),
                verify_artifacts: request.verify,
                file_id: request.file_id.clone(),
                rows: request.rows,
            },
            Task::Inspect(request),
        )?;
        self.rna_panel_ui.processing.evidence.loading = true;
        self.rna_panel_ui.page = Page::Inspect;
        Ok(())
    }

    pub(super) fn render_rna_inspector(&mut self, ui: &mut egui::Ui) {
        let mut action = None;
        let state = &mut self.rna_panel_ui.processing.evidence;
        if let Some(error) = &state.error {
            ui.colored_label(egui::Color32::RED, error);
        }
        ui.add_enabled_ui(!state.loading, |ui| {
            ui.horizontal_wrapped(|ui| {
                input(
                    ui,
                    "Saved panel / report",
                    &mut state.path,
                    "rna_panel.inspect.path",
                );
                if button(ui, "Choose report...", "rna_panel.inspect.choose")
                    && let Some(path) = rfd::FileDialog::new().pick_file()
                {
                    state.path = path.display().to_string();
                    state.binding_json.clear();
                }
                if button(ui, "Inspect saved file", "rna_panel.inspect.open") {
                    let binding = if state.binding_json.trim().is_empty() {
                        Ok(None)
                    } else {
                        serde_json::from_str::<RnaReadPanelResultBinding>(&state.binding_json)
                            .map(Some)
                            .map_err(|e| e.to_string())
                    };
                    match binding {
                        Ok(binding) => {
                            action = Some(InspectionRequest {
                                path: state.path.clone(),
                                expected_content: binding
                                    .as_ref()
                                    .map(|b| b.report_content_sha256.clone()),
                                binding,
                                file_id: None,
                                rows: 0,
                                verify: false,
                            })
                        }
                        Err(error) => state.error = Some(error),
                    }
                }
            });
            ui.collapsing("Retained result binding (optional JSON)", |ui| {
                let response =
                    ui.add(egui::TextEdit::multiline(&mut state.binding_json).desired_rows(4));
                capture_input(&response, "rna_panel.inspect.binding".into());
            });
            ui.collapsing("This session's published results", |ui| {
                for (i, binding) in state.published.iter().enumerate() {
                    if button(
                        ui,
                        &format!("{}: {}", binding.outcome, binding.report_path),
                        &format!("rna_panel.inspect.published.{i}"),
                    ) {
                        action = Some(InspectionRequest::bound(binding.clone()));
                    }
                }
            });
        });
        if state.loading {
            ui.spinner();
        }
        let Some(loaded) = &state.loaded else {
            if let Some(request) = action
                && let Err(error) = self.start_rna_inspection(request)
            {
                self.rna_panel_ui.processing.evidence.error = Some(error);
            }
            return;
        };
        ui.separator();
        ui.label(format!("Inspected: {}", loaded.path));
        ui.monospace(&loaded.content_sha256);
        ui.label(if loaded.binding.is_some() {
            "Exact published request/result binding checked"
        } else {
            "Saved metadata checked; no submitted-service binding supplied"
        });
        ui.horizontal_wrapped(|ui| {
            if button(ui, "Copy report JSON", "rna_panel.inspect.copy_json") {
                ui.ctx().copy_text(loaded.json.clone());
            }
            if let Some(binding) = &loaded.binding
                && button(ui, "Copy result binding", "rna_panel.inspect.copy_binding")
            {
                if let Ok(json) = serde_json::to_string_pretty(binding) {
                    ui.ctx().copy_text(json);
                }
            }
            if loaded.content.rna_read_panel_screen.is_some()
                && button(ui, "Add screen to summary", "rna_panel.inspect.add_summary")
            {
                if !state.request.report_paths.contains(&loaded.path) {
                    state.request.report_paths.push(loaded.path.clone());
                    state.reviewed_export = None;
                }
            }
        });
        egui::ScrollArea::vertical().id_salt("rna_inspection_content").show(ui, |ui| {
            if let Some(panel) = &loaded.content.rna_read_panel {
                ui.heading(format!("{} genes / {} loci / {} selected of {} transcripts", panel.gene_count, panel.locus_count, panel.selected_transcript_count, panel.annotated_transcript_count));
                if !panel.overlapping_source_pairs.is_empty() { ui.label(format!("Overlapping source pairs: {:?}", panel.overlapping_source_pairs)); }
                if let Some(shared) = &panel.shared_transcript_sequences {
                    ui.label("Shared mature sequences: review source provenance; sequence equality does not establish one gene.");
                    table(ui, "rna_shared_transcript_sequences", shared.len(), |ui, i| { let group = &shared[i]; ui.label(format!("{} | loci {:?} | members {:?}", group.sequence_sha256, group.locus_ids, group.member_ids)); });
                } else { ui.label("Shared mature-sequence diagnostic: not recorded in this historical panel."); }
                table(ui, "rna_panel_sources", panel.loci.len(), |ui, i| { let l = &panel.loci[i]; ui.label(format!("{} | {} | {} | {} | {} | {} | {}", l.source.locus_id, l.source.gene_symbol, l.source.gene_id, l.source.assembly, l.source.annotation_release, l.source.source, l.source_record_sha256)); });
                table(ui, "rna_panel_members", panel.templates.len(), |ui, i| { let t = &panel.templates[i]; ui.label(format!("{} | {} | {} | {} bases | junctions {:?} | source exons {:?}", t.member_id, t.transcript_id, t.strand, t.sequence.len(), t.junction_offsets, t.source_exons)); });
            }
            if let Some(report) = &loaded.content.rna_read_panel_screen {
                let outcome = ui.heading(format!("Scientific outcome: {}", report.status));
                capture_status(&outcome, "rna_panel.inspect.outcome", Some(&report.content_sha256), &report.status);
                ui.label(format!("Panel: {}", report.panel_sha256));
                if let Some(profile) = &report.index_work_profile { ui.label(format!("Index work (counts, not timings): {:?}", profile)); }
                else { ui.label("Index work profile: not recorded (historical report)."); }
                table(ui, "rna_file_work_profiles", report.files_completed.len(), |ui, i| { let file = &report.files_completed[i]; if let Some(profile) = &file.work_profile { ui.label(format!("{} work: {:?}", file.file_id, profile)); } else { ui.label(format!("{} work: not recorded", file.file_id)); } });
                ui.label(format!("Submitted request: {}", report.submitted_request_sha256.as_deref().unwrap_or("not recorded (legacy)")));
                if let Some(binding) = &loaded.binding { ui.label(format!("Effective request: {}", binding.effective_request_sha256)); }
                ui.label(format!("{} verified file units / {} requested", report.files_completed.len(), report.request.files.len()));
                if let Some(error) = &report.control_error {
                    ui.colored_label(egui::Color32::RED, format!("{}: {} | {} | {:?} | {:?}", error.code, error.message, error.marker_path, error.phase, error.file_id));
                }
                table(ui, "rna_verified_files", report.files_completed.len(), |ui, i| { let file = &report.files_completed[i]; ui.label(format!("{} | {} / {} / {} | {} reads | {} compressed bytes | MD5 {} | SHA {} | {:?}", file.file_id, file.input.sample_id, file.input.barcode, file.input.priming, file.reads, file.compressed_bytes, file.observed_md5, file.observed_sha256, file.classes)); });
                if let Some(missing) = &report.inputs_not_completed {
                    table(ui, "rna_missing_files", missing.len(), |ui, i| { let f = &missing[i]; ui.label(format!("Unassessed: {} / {} / {} | {}", f.sample_id, f.barcode, f.priming, f.path)); });
                } else if !report.files_not_completed.is_empty() { ui.label(format!("Unfinished legacy file IDs: {:?}; sample association unavailable", report.files_not_completed)); }
                ui.add_enabled_ui(!state.loading && !report.files_completed.is_empty(), |ui| {
                    let picker = egui::ComboBox::from_id_salt("rna_inspect_file").selected_text(if state.file_id.is_empty() { "Choose verified file" } else { &state.file_id }).show_ui(ui, |ui| {
                        for file in &report.files_completed { capture(&ui.selectable_value(&mut state.file_id, file.file_id.clone(), &file.file_id), &format!("rna_panel.inspect.file.{}", file.file_id)); }
                    });
                    capture(&picker.response, "rna_panel.inspect.file");
                    let response = ui.add(egui::DragValue::new(&mut state.rows).range(1..=1000));
                    capture_input(&response, "rna_panel.inspect.rows".into());
                    if button(ui, "Inspect nomination prefix", "rna_panel.inspect.prefix") && !state.file_id.is_empty() { action = Some(loaded.inspect_again(Some(state.file_id.clone()), state.rows, false)); }
                });
                if let Some(prefix) = &loaded.content.rna_read_panel_inspection {
                    ui.strong(format!("{} rows, artifact_content_verified={} (bounded prefix, not a whole-table verification)", prefix.rows.len(), prefix.artifact_content_verified));
                    // Each row may contain several candidates; keep every class/member/junction association.
                    egui::ScrollArea::both().id_salt("rna_prefix").max_height(300.0).show(ui, |ui| {
                        for row in &prefix.rows {
                            ui.label(format!("{} / {} | {} | record {} | offset {} | {}", row.sample_id, row.barcode, row.read_id, row.record_index, row.decompressed_record_offset, row.class));
                            for candidate in &row.candidates { ui.label(format!("{} | {} / {} | {} | {} chain bases | read {:?} | junctions {:?}", candidate.member_id, candidate.locus_id, candidate.transcript_id, candidate.read_orientation, candidate.chain_bases, candidate.read_interval_0based, candidate.supported_junction_offsets)); }
                        }
                    });
                }
                for nonclaim in &report.nonclaims { ui.label(nonclaim); }
            }
            if let Some(summary) = &loaded.content.rna_read_panel_summary {
                let outcome = ui.heading(format!("Summary coverage: {}", summary.completeness));
                capture_status(&outcome, "rna_panel.inspect.outcome", Some(&summary.content_sha256), &summary.completeness);
                ui.label(format!("Panel: {}", summary.panel_sha256));
                ui.label(format!("Original inputs reverified: {}; nomination sidecars reverified: {}", summary.original_inputs_reverified, summary.nomination_sidecars_reverified));
                ui.label(format!("Class/member counts complete: {:?}", summary.class_member_counts_complete));
                table(ui, "rna_summary_runs", summary.source_runs.len(), |ui, i| { let r = &summary.source_runs[i]; ui.label(format!("{} | {} | {} | verified {} | missing {:?} | {:?}", r.report_path, r.report_sha256, r.status, r.completed_files, r.files_not_completed, r.control_error)); });
                if let Some(coverage) = &summary.sample_coverage {
                    table(ui, "rna_summary_coverage", coverage.len(), |ui, i| { let r = &coverage[i]; ui.label(format!("{} / {} / {} | {} | files {}/{} | reads {}", r.sample_id, r.barcode, r.priming, r.status, r.verified_files, r.requested_files, r.reads.map(|n| n.to_string()).unwrap_or_else(|| "unassessed".into()))); });
                } else { ui.label("Legacy summary: requested sample/input coverage unavailable"); }
                table(ui, "rna_summary_samples", summary.samples.len(), |ui, i| { let r = &summary.samples[i]; ui.label(format!("{} | {} verified reads | {} single-candidate | co-nominated {:?} | historical ambiguous {:?} | {:?}", r.sample_id, r.reads, r.single_candidate_reads, r.co_nominated_reads, r.ambiguous_reads, r.classes)); });
                table(ui, "rna_summary_transcripts", summary.transcripts.len(), |ui, i| { let r = &summary.transcripts[i]; ui.label(format!("{} | {} / {} | {} | {} | single {} | co-nominated {:?} | historical {:?}", r.sample_id, r.gene_symbol, r.gene_id, r.member_id, r.transcript_id, r.single_candidate_reads, r.co_nominated_reads, r.ambiguous_compatible_reads)); });
                if let Some(sets) = &summary.co_nomination_sets {
                    table(ui, "rna_summary_sets", sets.len(), |ui, i| { let r = &sets[i]; ui.label(format!("{} | {} | {} reads | {:?} | sequence sharing {:?}", r.sample_id, r.class, r.reads, r.member_ids, r.sequence_sharing)); });
                }
                if let Some(sets) = &summary.ambiguous_sets {
                    ui.label("legacy_class_unresolved: class association cannot be recovered");
                    table(ui, "rna_summary_legacy_sets", sets.len(), |ui, i| { let r = &sets[i]; ui.label(format!("{} | {} reads | {:?}", r.sample_id, r.reads, r.member_ids)); });
                }
                for artifact in &summary.artifacts { ui.label(format!("{} | {} bytes | {}", artifact.path, artifact.bytes, artifact.sha256)); }
                if let Some(verification) = &loaded.content.rna_read_panel_verification {
                    let outcome = ui.strong(format!("Tables verified: {} bytes | {}", verification.verified_bytes, verification.scope));
                    capture_status(&outcome, "rna_panel.inspect.verification", Some(&verification.summary_sha256), "verified_tables_only");
                }
                if !state.loading && button(ui, "Verify summary tables", "rna_panel.inspect.verify") { action = Some(loaded.inspect_again(None, 0, true)); }
                for nonclaim in &summary.nonclaims { ui.label(nonclaim); }
            }
        });
        if let Some(request) = action
            && let Err(error) = self.start_rna_inspection(request)
        {
            self.rna_panel_ui.processing.evidence.error = Some(error);
        }
    }

    pub(super) fn render_rna_export(&mut self, ui: &mut egui::Ui) {
        let mut export = false;
        let mut open = None;
        let state = &mut self.rna_panel_ui.processing.evidence;
        if let Some(error) = &state.error {
            ui.colored_label(egui::Color32::RED, error);
        }
        ui.label(&state.export_detail);
        egui::ScrollArea::vertical()
            .id_salt("rna_summary_export")
            .show(ui, |ui| {
                ui.horizontal_wrapped(|ui| {
                    if button(ui, "Add screen reports...", "rna_panel.export.add")
                        && let Some(paths) = rfd::FileDialog::new()
                            .add_filter("Screen report", &["json"])
                            .pick_files()
                    {
                        for path in paths {
                            let path = path.display().to_string();
                            if !state.request.report_paths.contains(&path) {
                                state.request.report_paths.push(path);
                            }
                        }
                        state.reviewed_export = None;
                    }
                    if button(ui, "Choose summary parent...", "rna_panel.export.parent")
                        && let Some(path) = rfd::FileDialog::new().pick_folder()
                    {
                        state.request.output_dir = path
                            .join(format!("rna-summary-{}", Self::now_unix_ms()))
                            .display()
                            .to_string();
                        state.reviewed_export = None;
                    }
                });
                let mut remove = None;
                for (i, path) in state.request.report_paths.iter().enumerate() {
                    ui.horizontal_wrapped(|ui| {
                        ui.label(path);
                        if button(ui, "Remove", &format!("rna_panel.export.remove.{i}")) {
                            remove = Some(i);
                        }
                    });
                }
                if let Some(i) = remove {
                    state.request.report_paths.remove(i);
                    state.reviewed_export = None;
                }
                if input(
                    ui,
                    "New summary directory",
                    &mut state.request.output_dir,
                    "rna_panel.export.output",
                ) {
                    state.reviewed_export = None;
                }
                if button(ui, "Review summary request", "rna_panel.export.review") {
                    if state.request.report_paths.is_empty() || state.request.output_dir.is_empty()
                    {
                        state.error = Some(
                            "Select explicit screen reports and a fresh output directory".into(),
                        );
                    } else if let Ok(json) = serde_json::to_string_pretty(&state.request) {
                        state.reviewed_export = Some((state.request.clone(), json));
                        state.error = None;
                    }
                }
                if let Some((_, json)) = &state.reviewed_export {
                    ui.monospace(json);
                    if button(ui, "Copy summary request", "rna_panel.export.copy") {
                        ui.ctx().copy_text(json.clone());
                    }
                    export = button(
                        ui,
                        "Export reviewed summary (writes files)",
                        "rna_panel.export.run",
                    );
                }
                if let Some(id) = state.export_job
                    && let Some(receipt) = self.agent_command_service.status(id)
                {
                    ui.label(format!(
                        "Export {id}: {:?} | {}",
                        receipt.state, receipt.phase
                    ));
                    if let Some(binding) = receipt.artifact_result
                        && button(ui, "Open bound summary", "rna_panel.export.inspect")
                    {
                        open = Some(binding);
                    }
                }
            });
        if export {
            if let Err(error) = self.start_reviewed_rna_export() {
                self.rna_panel_ui.processing.evidence.error = Some(error);
            }
        }
        if let Some(binding) = open
            && let Err(error) = self.start_rna_inspection(InspectionRequest::bound(binding))
        {
            self.rna_panel_ui.processing.evidence.error = Some(error);
        }
    }

    fn start_reviewed_rna_export(&mut self) -> Result<(), String> {
        let request = self.rna_panel_ui.processing.evidence.export_request()?;
        let id = self.submit_rna_operation(
            Operation::ExportRnaReadPanelSummary {
                request: Box::new(request),
            },
            Task::Export,
        )?;
        let state = &mut self.rna_panel_ui.processing.evidence;
        state.export_job = Some(id);
        state.reviewed_export = None;
        state.export_detail =
            "Export admitted; source compatibility is checked by the shared engine".into();
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn drain(app: &mut GENtleApp) {
        let deadline = std::time::Instant::now() + std::time::Duration::from_secs(15);
        while app.rna_panel_pending() && std::time::Instant::now() < deadline {
            app.poll_rna_panel_action();
            std::thread::sleep(std::time::Duration::from_millis(1));
        }
        assert!(
            !app.rna_panel_pending(),
            "native command-service test timed out"
        );
        assert!(
            app.rna_panel_ui.processing.evidence.error.is_none(),
            "{:?}",
            app.rna_panel_ui.processing.evidence.error
        );
    }

    #[test]
    fn rna_native_screen_monitor_bound_reopen_and_summary_work_without_the_live_project_lock() {
        // Complete app-control path over the public fixed-LCG rehearsal, not native acceptance.
        let td = tempfile::tempdir().unwrap();
        let root = td.path().join("native synthetic");
        GentleEngine::default()
            .prepare_rna_read_panel_demo(&root.display().to_string())
            .unwrap();
        let read = |file: &str| std::fs::read(root.join("artifacts").join(file)).unwrap();
        let engine = GentleEngine::from_state(
            serde_json::from_slice::<crate::engine::ProjectState>(&read("starter.json")).unwrap(),
        );
        let build =
            serde_json::from_slice::<BuildRnaReadPanelRequest>(&read("build.json")).unwrap();
        let panel = engine.build_rna_read_panel(&build).unwrap();
        assert_eq!(
            engine
                .rna_read_panel_gene_set_sources(&build.gene_set.unwrap().resolution_id)
                .unwrap()
                .1
                .len(),
            5
        );
        let mut request: ScreenRnaReadPanelRequest =
            serde_json::from_slice(&read("screen.json")).unwrap();
        request.expected_panel_sha256 = Some(panel.content_sha256);
        let mut app = GENtleApp::default();
        app.rna_panel_ui.processing.draft = request;
        app.rna_panel_ui.processing.review().unwrap();
        let shared = app.engine.clone();
        let mut live = shared.write().unwrap();
        live.auxiliary_metadata_mut()
            .insert("unrelated_live_edit".into(), serde_json::json!(true));
        let journal = live.journal_len();
        app.start_rna_screen().unwrap();
        app.close_rna_panel_dialog();
        drain(&mut app);
        let loaded = app
            .rna_panel_ui
            .processing
            .evidence
            .loaded
            .as_ref()
            .unwrap();
        let report = loaded.content.rna_read_panel_screen.as_ref().unwrap();
        assert_eq!(report.status, "completed");
        assert_eq!(report.files_completed.len(), 2);
        assert!(loaded.binding.is_some());
        let state = &mut app.rna_panel_ui.processing.evidence;
        state.request = serde_json::from_slice(&read("summary.json")).unwrap();
        state.reviewed_export = Some((
            state.request.clone(),
            serde_json::to_string(&state.request).unwrap(),
        ));
        app.start_reviewed_rna_export().unwrap();
        drain(&mut app);
        let loaded = app
            .rna_panel_ui
            .processing
            .evidence
            .loaded
            .as_ref()
            .unwrap();
        let summary = loaded.content.rna_read_panel_summary.as_ref().unwrap();
        assert_eq!(summary.completeness, "complete");
        assert_eq!(summary.samples.len(), 2);
        assert!(summary.samples.iter().all(|s| s.reads == 4));
        assert!(!summary.co_nomination_sets.as_ref().unwrap().is_empty());
        assert!(loaded.binding.is_some());
        assert_eq!(live.state().metadata["unrelated_live_edit"], true);
        assert_eq!(
            live.journal_len(),
            journal,
            "file-only native runs never merge their private journal"
        );
    }

    fn inspect(engine: &mut GentleEngine, request: &InspectionRequest) -> serde_json::Value {
        let (command, _) = op_command(Operation::InspectRnaReadPanel {
            path: request.path.clone(),
            expected_result: request.binding.clone().map(Box::new),
            verify_artifacts: request.verify,
            file_id: request.file_id.clone(),
            rows: request.rows,
        })
        .unwrap();
        crate::engine_shell::execute_shell_command(engine, &command)
            .unwrap()
            .output
    }

    #[test]
    fn rna_native_inspection_and_export_project_the_shared_reports_without_resolving_sets() {
        // Existing fixed-LCG demo: tiny synthetic loci/reads, never human or private input.
        let td = tempfile::tempdir().unwrap();
        let root = td.path().join("synthetic demo");
        GentleEngine::default()
            .prepare_rna_read_panel_demo(&root.display().to_string())
            .unwrap();
        let read = |name: &str| std::fs::read(root.join("artifacts").join(name)).unwrap();
        let state =
            serde_json::from_slice::<crate::engine::ProjectState>(&read("starter.json")).unwrap();
        let mut engine = GentleEngine::from_state(state);
        let build =
            serde_json::from_slice::<BuildRnaReadPanelRequest>(&read("build.json")).unwrap();
        let panel = engine.build_rna_read_panel(&build).unwrap();
        let mut request =
            serde_json::from_slice::<ScreenRnaReadPanelRequest>(&read("screen.json")).unwrap();
        request.expected_panel_sha256 = Some(panel.content_sha256);
        let report = engine
            .screen_rna_read_panel(request.clone(), &mut |_| true)
            .unwrap();
        let binding =
            GentleEngine::rna_read_panel_screen_result_binding(&request, &report).unwrap();
        let mut inspection = InspectionRequest::bound(binding);
        let loaded = Loaded::from_result(
            InspectionRequest::bound(inspection.binding.clone().unwrap()),
            inspect(&mut engine, &inspection),
        )
        .unwrap();
        assert_eq!(
            serde_json::from_str::<serde_json::Value>(&loaded.json).unwrap(),
            serde_json::to_value(&report).unwrap()
        );
        inspection.file_id = Some(report.files_completed[0].file_id.clone());
        inspection.rows = 4;
        let loaded = Loaded::from_result(
            loaded.inspect_again(inspection.file_id.clone(), 4, false),
            inspect(&mut engine, &inspection),
        )
        .unwrap();
        let prefix = loaded.content.rna_read_panel_inspection.as_ref().unwrap();
        assert_eq!(prefix.rows.len(), 4);
        assert!(!prefix.artifact_content_verified);
        assert!(
            prefix.rows[0].candidates.len() >= 2,
            "shared exon does not earn an isoform winner"
        );
        let summary_request =
            serde_json::from_slice::<ExportRnaReadPanelSummaryRequest>(&read("summary.json"))
                .unwrap();
        let summary = engine
            .export_rna_read_panel_summary(&summary_request)
            .unwrap();
        let binding =
            GentleEngine::rna_read_panel_summary_result_binding(&summary_request, &summary)
                .unwrap();
        let inspection = InspectionRequest::bound(binding.clone());
        let loaded = Loaded::from_result(
            InspectionRequest::bound(binding),
            inspect(&mut engine, &inspection),
        )
        .unwrap();
        assert_eq!(
            serde_json::from_str::<serde_json::Value>(&loaded.json).unwrap(),
            serde_json::to_value(&summary).unwrap()
        );
        assert_eq!(
            loaded
                .content
                .rna_read_panel_summary
                .as_ref()
                .unwrap()
                .co_nomination_sets
                .as_ref()
                .unwrap()
                .len(),
            summary.co_nomination_sets.as_ref().unwrap().len()
        );
        let request = loaded.inspect_again(None, 0, true);
        let verified = Loaded::from_result(
            loaded.inspect_again(None, 0, true),
            inspect(&mut engine, &request),
        )
        .unwrap();
        assert!(verified.content.rna_read_panel_verification.is_some());
        let mut stale = loaded.inspect_again(None, 0, false);
        stale.expected_content = Some("sha256:replaced".into());
        assert!(Loaded::from_result(stale, inspect(&mut engine, &inspection)).is_err());
    }

    #[test]
    fn rna_native_export_review_is_invalidated_by_path_changes() {
        let mut state = Evidence::default();
        state
            .request
            .report_paths
            .push("synthetic/report.json".into());
        state.request.output_dir = "new-summary".into();
        assert!(state.export_request().is_err());
        state.reviewed_export = Some((
            state.request.clone(),
            serde_json::to_string(&state.request).unwrap(),
        ));
        assert!(state.export_request().is_ok());
        state.request.report_paths.push("other/report.json".into());
        assert!(state.export_request().is_err());
        state.request.report_paths.pop();
        state.request.output_dir.push_str("-changed");
        assert!(state.export_request().is_err());
    }
}
