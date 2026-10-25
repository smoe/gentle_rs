//! Native request editor and nonmodal monitor over the app's shared command service.

use super::*;
use crate::command_execution::CommandReceipt;
use crate::engine::{
    RNA_PANEL_SCREEN_REQUEST_SCHEMA, RnaReadPanelInput, RnaReadPanelPolicy,
    RnaReadPanelScreenDraft, ScreenRnaReadPanelRequest,
};
use crate::engine_shell::{ShellExecutionOptions, quote_shell_arg};
use crate::runtime_status::RuntimeStatusFrameState;

#[path = "processing/evidence.rs"]
mod evidence;

struct Reviewed {
    request: ScreenRnaReadPanelRequest,
    json: String,
}

enum Task {
    Panel(String),
    Draft { load: bool },
    Screen,
    Cancel(u64),
    Inspect(evidence::InspectionRequest),
    Export,
    Demo(String),
}

struct PendingTask {
    id: u64,
    task: Task,
}

struct Run {
    id: u64,
    request: ScreenRnaReadPanelRequest,
    receipt: Option<CommandReceipt>,
    detail: String,
    cancel_pending: bool,
    cancel_requested: bool,
}

pub(super) struct Processing {
    draft: ScreenRnaReadPanelRequest,
    review: Option<Reviewed>,
    pending: Vec<PendingTask>,
    runs: Vec<Run>,
    sample: String,
    barcode: String,
    priming: String,
    parent: String,
    serial: u64,
    status: String,
    error: Option<String>,
    evidence: evidence::Evidence,
    demo: Option<String>,
}

impl Default for Processing {
    fn default() -> Self {
        Self {
            draft: ScreenRnaReadPanelRequest {
                schema: RNA_PANEL_SCREEN_REQUEST_SCHEMA.into(),
                panel_path: String::new(),
                expected_panel_sha256: None,
                files: vec![],
                md5_manifest_path: None,
                output_dir: String::new(),
                policy: RnaReadPanelPolicy::default(),
                cancel_path: None,
            },
            review: None,
            pending: vec![],
            runs: vec![],
            sample: String::new(),
            barcode: String::new(),
            priming: "custom_only".into(),
            parent: String::new(),
            serial: 0,
            status: String::new(),
            error: None,
            evidence: evidence::Evidence::default(),
            demo: None,
        }
    }
}

impl Processing {
    pub(super) fn pending(&self) -> bool {
        !self.pending.is_empty()
    }

    fn loading(&self) -> bool {
        self.pending
            .iter()
            .any(|p| matches!(p.task, Task::Panel(_) | Task::Draft { .. }))
    }

    fn fresh_paths(&mut self) {
        self.serial += 1;
        let stem = format!("rna-screen-{}-{}", GENtleApp::now_unix_ms(), self.serial);
        let parent = std::path::Path::new(&self.parent);
        self.draft.output_dir = parent.join(&stem).display().to_string();
        self.draft.cancel_path = Some(
            parent
                .join(format!("{stem}.cancel.json"))
                .display()
                .to_string(),
        );
        self.review = None;
    }

    fn review(&mut self) -> Result<(), String> {
        if self.draft.expected_panel_sha256.is_none() {
            return Err("Inspect the reference panel before reviewing this run".into());
        }
        if self.draft.files.is_empty() || self.draft.files.iter().any(|f| !f.declared_complete) {
            return Err("Select files and explicitly confirm that each copy has finished".into());
        }
        if self.draft.output_dir.is_empty()
            || self.draft.cancel_path.as_deref().is_none_or(str::is_empty)
        {
            return Err("Choose fresh run and control paths".into());
        }
        let json = serde_json::to_string_pretty(&self.draft).map_err(|e| e.to_string())?;
        self.review = Some(Reviewed {
            request: self.draft.clone(),
            json,
        });
        Ok(())
    }

    fn reviewed_request(&self) -> Result<ScreenRnaReadPanelRequest, String> {
        let reviewed = self.review.as_ref().ok_or("Review the request first")?;
        if serde_json::to_value(&reviewed.request).ok() != serde_json::to_value(&self.draft).ok() {
            return Err("Draft changed; review the new request before running".into());
        }
        Ok(reviewed.request.clone())
    }
}

fn terminal(receipt: &CommandReceipt) -> bool {
    matches!(
        receipt.state,
        RuntimeStatusFrameState::Completed
            | RuntimeStatusFrameState::Failed
            | RuntimeStatusFrameState::Cancelled
    )
}

fn button(ui: &mut egui::Ui, text: &str, id: &str) -> bool {
    let response = ui.button(text);
    capture(&response, id);
    response.clicked()
}

fn input(ui: &mut egui::Ui, label: &str, value: &mut String, id: &str) -> bool {
    ui.label(label);
    let response = ui.add(egui::TextEdit::singleline(value).desired_width(300.0));
    capture_input(&response, id.into());
    response.changed()
}

fn optional_input(ui: &mut egui::Ui, label: &str, value: &mut Option<String>, id: &str) -> bool {
    let mut text = value.clone().unwrap_or_default();
    let changed = input(ui, label, &mut text, id);
    if changed {
        *value = (!text.is_empty()).then_some(text);
    }
    changed
}

fn op_command(op: Operation) -> Result<(ShellCommand, String), String> {
    let payload = serde_json::to_string(&op).map_err(|e| e.to_string())?;
    let text = format!("op {}", quote_shell_arg(&payload));
    Ok((ShellCommand::Op { payload }, text))
}

impl GENtleApp {
    fn submit_rna_processing(
        &mut self,
        command: ShellCommand,
        text: String,
        task: Task,
    ) -> Result<u64, String> {
        let id = self.agent_command_service.submit_parsed(
            self.engine.clone(),
            command,
            text,
            ShellExecutionOptions::default(),
        )?;
        self.rna_panel_ui
            .processing
            .pending
            .push(PendingTask { id, task });
        Ok(id)
    }

    fn submit_rna_operation(&mut self, op: Operation, task: Task) -> Result<u64, String> {
        let (command, text) = op_command(op)?;
        self.submit_rna_processing(command, text, task)
    }

    fn start_rna_screen(&mut self) -> Result<(), String> {
        let request = self.rna_panel_ui.processing.reviewed_request()?;
        let id = self.submit_rna_operation(
            Operation::ScreenRnaReadPanel {
                request: Box::new(request.clone()),
            },
            Task::Screen,
        )?;
        let state = &mut self.rna_panel_ui.processing;
        // Bound session history, never delete published files or an active worker.
        if state.runs.len() >= 32 {
            if let Some(index) = state
                .runs
                .iter()
                .position(|r| r.receipt.as_ref().is_some_and(terminal))
            {
                state.runs.remove(index);
            }
        }
        state.runs.push(Run {
            id,
            request,
            receipt: None,
            detail: String::new(),
            cancel_pending: false,
            cancel_requested: false,
        });
        state.review = None;
        self.rna_panel_ui.page = Page::Runs;
        Ok(())
    }

    pub(super) fn poll_rna_processing_jobs(&mut self) {
        let mut published = vec![];
        for run in &mut self.rna_panel_ui.processing.runs {
            if let Some(receipt) = self.agent_command_service.status(run.id) {
                run.receipt = Some(receipt);
            }
        }
        let mut index = 0;
        while index < self.rna_panel_ui.processing.pending.len() {
            let id = self.rna_panel_ui.processing.pending[index].id;
            let Some(result) = self.agent_command_service.take_result(id) else {
                index += 1;
                continue;
            };
            let task = self.rna_panel_ui.processing.pending.remove(index).task;
            let receipt = self.agent_command_service.status(id);
            let state = &mut self.rna_panel_ui.processing;
            match task {
                Task::Screen => {
                    if let Some(run) = state.runs.iter_mut().find(|r| r.id == id) {
                        run.receipt = receipt;
                        if let Some(binding) =
                            run.receipt.as_ref().and_then(|r| r.artifact_result.clone())
                        {
                            published.push(binding);
                        }
                        run.detail = match result {
                            Ok(_) => run
                                .receipt
                                .as_ref()
                                .and_then(|r| r.error.clone())
                                .unwrap_or_else(|| {
                                    "Published result retained; no project attachment".into()
                                }),
                            Err(error) => error,
                        };
                    }
                }
                Task::Cancel(run_id) => {
                    if let Some(run) = state.runs.iter_mut().find(|r| r.id == run_id) {
                        run.cancel_pending = false;
                        match result {
                            Ok(_) => {
                                run.cancel_requested = true;
                                run.detail = if run.receipt.as_ref().is_some_and(terminal) {
                                    "Cancellation marker written after termination; published outcome unchanged"
                                } else { "Cancellation marker written; termination is not yet confirmed" }.into();
                            }
                            Err(error) => {
                                run.detail = format!(
                                    "Cancellation request failed: {error}; worker outcome remains unchanged"
                                );
                            }
                        }
                    }
                }
                Task::Panel(path) => {
                    let parsed = result.and_then(|r| {
                        serde_json::from_value::<RnaReadPanel>(
                            r.output["result"]["rna_read_panel"].clone(),
                        )
                        .map_err(|e| e.to_string())
                    });
                    match parsed {
                        Ok(panel) => {
                            state.draft.panel_path = path;
                            state.draft.expected_panel_sha256 = Some(panel.content_sha256.clone());
                            state.review = None;
                            state.status = format!(
                                "Reviewed panel: {} genes, {} selected transcripts; {}",
                                panel.gene_count,
                                panel.selected_transcript_count,
                                panel.content_sha256
                            );
                        }
                        Err(error) => {
                            state.draft.expected_panel_sha256 = None;
                            state.error = Some(error);
                        }
                    }
                }
                Task::Draft { load } => match result.and_then(|r| {
                    serde_json::from_value::<RnaReadPanelScreenDraft>(r.output)
                        .map_err(|e| e.to_string())
                }) {
                    Ok(draft) => {
                        if load {
                            state.draft = draft.request;
                            // A serialized digest is a binding, not evidence that this session inspected the panel.
                            state.draft.expected_panel_sha256 = None;
                            state.review = None;
                        }
                        state.status = draft
                            .artifact
                            .map(|a| format!("Saved request {} ({})", a.path, a.sha256))
                            .unwrap_or_else(|| {
                                "Draft loaded; inspect the panel and review every file binding"
                                    .into()
                            });
                    }
                    Err(error) => state.error = Some(error),
                },
                Task::Inspect(request) => {
                    state.evidence.loading = false;
                    match result.and_then(|run| evidence::Loaded::from_result(request, run.output))
                    {
                        Ok(loaded) => {
                            state.evidence.loaded = Some(loaded);
                            state.evidence.error = None;
                        }
                        Err(error) => state.evidence.error = Some(error),
                    }
                }
                Task::Export => {
                    if let Some(binding) = receipt.as_ref().and_then(|r| r.artifact_result.clone())
                    {
                        published.push(binding);
                    }
                    state.evidence.export_detail = match result {
                        Ok(_) => receipt
                            .as_ref()
                            .and_then(|r| r.error.clone())
                            .unwrap_or_else(|| "Published summary retained".into()),
                        Err(error) => error,
                    };
                }
                Task::Demo(path) => match result {
                    Ok(_) => {
                        state.demo = Some(path.clone());
                        state.status = format!("Synthetic rehearsal saved: {path}");
                    }
                    Err(error) => state.error = Some(error),
                },
            }
        }
        for binding in &published {
            self.rna_panel_ui
                .processing
                .evidence
                .published(binding.clone());
        }
        if let Some(binding) = published.pop() {
            if !self.rna_panel_ui.processing.evidence.loading {
                if let Err(error) =
                    self.start_rna_inspection(evidence::InspectionRequest::bound(binding))
                {
                    self.rna_panel_ui.processing.evidence.error = Some(error);
                }
            }
        }
    }

    pub(super) fn render_rna_processing(&mut self, ui: &mut egui::Ui) {
        if self.rna_panel_ui.page == Page::Inspect {
            self.render_rna_inspector(ui);
            return;
        }
        if self.rna_panel_ui.page == Page::Export {
            self.render_rna_export(ui);
            return;
        }
        if self.rna_panel_ui.page == Page::Runs {
            self.render_rna_runs(ui);
            return;
        }
        let mut operation = None;
        let mut load = None;
        let mut save = None;
        let mut run = false;
        let prepared = self
            .rna_panel_ui
            .base_path
            .as_ref()
            .zip(self.rna_panel_ui.panel.as_ref())
            .map(|(path, panel)| (path.clone(), panel.content_sha256.clone()));
        let state = &mut self.rna_panel_ui.processing;
        if let Some(error) = &state.error {
            ui.colored_label(egui::Color32::RED, error);
        }
        ui.label(&state.status);
        let loading = state.loading();
        egui::ScrollArea::vertical()
            .id_salt("rna_screen_form")
            .show(ui, |ui| {
                ui.add_enabled_ui(!loading, |ui| {
                    let mut changed = false;
                    ui.horizontal_wrapped(|ui| {
                        if button(ui, "Load request...", "rna_panel.screen.load") {
                            load = rfd::FileDialog::new()
                                .add_filter("Screen request", &["json"])
                                .pick_file();
                        }
                        if button(ui, "Save request...", "rna_panel.screen.save") {
                            save = rfd::FileDialog::new()
                                .add_filter("Screen request", &["json"])
                                .set_file_name("screen.json")
                                .save_file();
                        }
                        if button(ui, "Use prepared panel", "rna_panel.screen.use_panel") {
                            if let Some((path, hash)) = &prepared {
                                state.draft.panel_path = path.clone();
                                state.draft.expected_panel_sha256 = Some(hash.clone());
                                changed = true;
                            } else {
                                state.error = Some("Open or build a panel in Prepare first".into());
                            }
                        }
                    });
                    ui.horizontal_wrapped(|ui| {
                        if input(
                            ui,
                            "Reference panel",
                            &mut state.draft.panel_path,
                            "rna_panel.screen.panel",
                        ) {
                            state.draft.expected_panel_sha256 = None;
                            changed = true;
                        }
                        if button(ui, "Choose...", "rna_panel.screen.choose_panel")
                            && let Some(path) = rfd::FileDialog::new().pick_file()
                        {
                            state.draft.panel_path = path.display().to_string();
                            state.draft.expected_panel_sha256 = None;
                            changed = true;
                        }
                        if button(ui, "Inspect panel", "rna_panel.screen.inspect_panel") {
                            operation = Some((
                                Operation::InspectRnaReadPanel {
                                    path: state.draft.panel_path.clone(),
                                    expected_result: None,
                                    verify_artifacts: false,
                                    file_id: None,
                                    rows: 0,
                                },
                                Task::Panel(state.draft.panel_path.clone()),
                            ));
                        }
                    });
                    let pin = ui.monospace(
                        state
                            .draft
                            .expected_panel_sha256
                            .as_deref()
                            .unwrap_or("Panel not inspected"),
                    );
                    capture_status(
                        &pin,
                        "rna_panel.screen.panel_binding",
                        state.draft.expected_panel_sha256.as_deref(),
                        if state.draft.expected_panel_sha256.is_some() {
                            "reviewed"
                        } else {
                            "not_evaluated"
                        },
                    );
                    ui.separator();
                    ui.horizontal_wrapped(|ui| {
                        input(
                            ui,
                            "New files: sample",
                            &mut state.sample,
                            "rna_panel.screen.sample",
                        );
                        input(
                            ui,
                            "Barcode",
                            &mut state.barcode,
                            "rna_panel.screen.barcode",
                        );
                        input(
                            ui,
                            "Priming",
                            &mut state.priming,
                            "rna_panel.screen.priming",
                        );
                    });
                    ui.horizontal_wrapped(|ui| {
                        for (sample, barcode) in
                            [("26R266-1", "barcode13"), ("26R266-2", "barcode14")]
                        {
                            if button(
                                ui,
                                &format!("{sample} / {barcode}"),
                                &format!("rna_panel.screen.preset.{barcode}"),
                            ) {
                                state.sample = sample.into();
                                state.barcode = barcode.into();
                                state.priming = "custom_only".into();
                            }
                        }
                        if button(ui, "Add FASTQ files...", "rna_panel.screen.add_files")
                            && let Some(paths) = rfd::FileDialog::new()
                                .add_filter("FASTQ / gzip", &["fastq", "fq", "gz"])
                                .pick_files()
                        {
                            for path in paths {
                                state.draft.files.push(RnaReadPanelInput {
                                    sample_id: state.sample.clone(),
                                    barcode: state.barcode.clone(),
                                    priming: state.priming.clone(),
                                    path: path.display().to_string(),
                                    declared_complete: false,
                                    expected_md5: None,
                                    expected_sha256: None,
                                    checksum_manifest_key: None,
                                });
                            }
                            changed = true;
                        }
                    });
                    let mut remove = None;
                    for (i, file) in state.draft.files.iter_mut().enumerate() {
                        ui.push_id(i, |ui| {
                            ui.separator();
                            changed |= input(
                                ui,
                                "FASTQ",
                                &mut file.path,
                                &format!("rna_panel.screen.file.{i}.path"),
                            );
                            ui.horizontal_wrapped(|ui| {
                                changed |= input(
                                    ui,
                                    "Sample",
                                    &mut file.sample_id,
                                    &format!("rna_panel.screen.file.{i}.sample"),
                                );
                                changed |= input(
                                    ui,
                                    "Barcode",
                                    &mut file.barcode,
                                    &format!("rna_panel.screen.file.{i}.barcode"),
                                );
                                changed |= input(
                                    ui,
                                    "Priming",
                                    &mut file.priming,
                                    &format!("rna_panel.screen.file.{i}.priming"),
                                );
                            });
                            let complete =
                                ui.checkbox(&mut file.declared_complete, "File copy finished");
                            capture_selection(
                                &complete,
                                &format!("rna_panel.screen.file.{i}.complete"),
                                file.declared_complete,
                            );
                            changed |= complete.changed();
                            ui.collapsing("Checksum binding", |ui| {
                                changed |= optional_input(
                                    ui,
                                    "MD5",
                                    &mut file.expected_md5,
                                    &format!("rna_panel.screen.file.{i}.md5"),
                                );
                                changed |= optional_input(
                                    ui,
                                    "SHA-256",
                                    &mut file.expected_sha256,
                                    &format!("rna_panel.screen.file.{i}.sha256"),
                                );
                                changed |= optional_input(
                                    ui,
                                    "Exact manifest key",
                                    &mut file.checksum_manifest_key,
                                    &format!("rna_panel.screen.file.{i}.key"),
                                );
                            });
                            if button(
                                ui,
                                "Remove file",
                                &format!("rna_panel.screen.file.{i}.remove"),
                            ) {
                                remove = Some(i);
                            }
                        });
                    }
                    if let Some(i) = remove {
                        state.draft.files.remove(i);
                        changed = true;
                    }
                    ui.horizontal_wrapped(|ui| {
                        changed |= optional_input(
                            ui,
                            "MD5 manifest",
                            &mut state.draft.md5_manifest_path,
                            "rna_panel.screen.manifest",
                        );
                        if button(ui, "Choose manifest...", "rna_panel.screen.choose_manifest")
                            && let Some(path) = rfd::FileDialog::new().pick_file()
                        {
                            state.draft.md5_manifest_path = Some(path.display().to_string());
                            changed = true;
                        }
                    });
                    ui.collapsing("Screen policy", |ui| {
                        macro_rules! count {
                            ($field:ident, $label:literal) => {
                                ui.horizontal_wrapped(|ui| {
                                    ui.label($label);
                                    let response = ui
                                        .add(egui::DragValue::new(&mut state.draft.policy.$field));
                                    capture_input(
                                        &response,
                                        concat!("rna_panel.screen.policy.", stringify!($field))
                                            .into(),
                                    );
                                    changed |= response.changed();
                                });
                            };
                        }
                        count!(kmer_len, "Seed length");
                        count!(min_chain_bases, "Minimum chain bases");
                        count!(min_junction_flank_bp, "Junction flank bases");
                        count!(max_index_postings, "Posting budget");
                        count!(max_anchors_per_read, "Anchors per read budget");
                        count!(max_equivalence_classes, "Nomination-set budget");
                        changed |= ui
                            .checkbox(
                                &mut state.draft.policy.write_locus_fastas,
                                "Write per-locus FASTA (reads may occur in multiple loci)",
                            )
                            .changed();
                    });
                    ui.horizontal_wrapped(|ui| {
                        if button(
                            ui,
                            "Choose output parent...",
                            "rna_panel.screen.output_parent",
                        ) && let Some(path) = rfd::FileDialog::new().pick_folder()
                        {
                            state.parent = path.display().to_string();
                            state.fresh_paths();
                            changed = true;
                        }
                        if !state.parent.is_empty()
                            && button(ui, "Fresh run paths", "rna_panel.screen.fresh")
                        {
                            state.fresh_paths();
                            changed = true;
                        }
                    });
                    changed |= input(
                        ui,
                        "New output directory",
                        &mut state.draft.output_dir,
                        "rna_panel.screen.output",
                    );
                    changed |= optional_input(
                        ui,
                        "New cancellation marker",
                        &mut state.draft.cancel_path,
                        "rna_panel.screen.cancel_path",
                    );
                    if changed {
                        state.review = None;
                    }
                    if button(ui, "Review request", "rna_panel.screen.review") {
                        state.error = state.review().err();
                    }
                    if let Some(review) = &state.review {
                        ui.collapsing("Reviewed request JSON", |ui| {
                            ui.monospace(&review.json);
                        });
                        if button(ui, "Copy JSON", "rna_panel.screen.copy_json") {
                            ui.ctx().copy_text(review.json.clone());
                        }
                        if button(ui, "Copy CLI command", "rna_panel.screen.copy_command") {
                            ui.ctx().copy_text(format!(
                                "rna-reads panel-screen {}",
                                quote_shell_arg(&review.json)
                            ));
                        }
                        run = button(
                            ui,
                            "Run reviewed screen (writes files)",
                            "rna_panel.screen.run",
                        );
                    }
                });
            });
        let result = if run {
            self.start_rna_screen().map(|_| ())
        } else if let Some((op, task)) = operation {
            self.submit_rna_operation(op, task).map(|_| ())
        } else if let Some(path) = load {
            let text = format!(
                "rna-reads panel-screen --request {}",
                quote_shell_arg(&format!("@{}", path.display()))
            );
            self.agent_command_service
                .submit(self.engine.clone(), text, ShellExecutionOptions::default())
                .map(|id| {
                    self.rna_panel_ui.processing.pending.push(PendingTask {
                        id,
                        task: Task::Draft { load: true },
                    });
                })
        } else if let Some(path) = save {
            let request = self.rna_panel_ui.processing.draft.clone();
            let output_path = path.display().to_string();
            let text = format!(
                "rna-reads panel-screen --request {} --save {}",
                quote_shell_arg(&serde_json::to_string(&request).unwrap_or_default()),
                quote_shell_arg(&output_path)
            );
            self.submit_rna_processing(
                ShellCommand::RnaReadPanelScreenDraft {
                    request: Box::new(request),
                    output_path: Some(output_path),
                },
                text,
                Task::Draft { load: false },
            )
            .map(|_| ())
        } else {
            Ok(())
        };
        if let Err(error) = result {
            self.rna_panel_ui.processing.error = Some(error);
        }
    }

    pub(super) fn render_rna_rehearsal_controls(&mut self, ui: &mut egui::Ui) {
        let mut prepare = None;
        ui.horizontal_wrapped(|ui| {
            if button(
                ui,
                "Prepare synthetic rehearsal... (writes files)",
                "rna_panel.demo.prepare",
            ) && let Some(parent) = rfd::FileDialog::new().pick_folder()
            {
                prepare = Some(
                    parent
                        .join(format!("rna-panel-demo-{}", Self::now_unix_ms()))
                        .display()
                        .to_string(),
                );
            }
            if let Some(root) = &self.rna_panel_ui.processing.demo {
                ui.label(root);
                if button(ui, "Copy starter path", "rna_panel.demo.starter") {
                    ui.ctx().copy_text(
                        std::path::Path::new(root)
                            .join("artifacts/starter.json")
                            .display()
                            .to_string(),
                    );
                }
                if button(ui, "Copy screen-request path", "rna_panel.demo.request") {
                    ui.ctx().copy_text(
                        std::path::Path::new(root)
                            .join("artifacts/screen.json")
                            .display()
                            .to_string(),
                    );
                }
            }
        });
        if let Some(path) = prepare
            && let Err(error) = self.submit_rna_operation(
                Operation::PrepareRnaReadPanelDemo {
                    output_dir: path.clone(),
                },
                Task::Demo(path),
            )
        {
            self.rna_panel_ui.processing.error = Some(error);
        }
        if let Some(error) = &self.rna_panel_ui.processing.error {
            ui.colored_label(egui::Color32::RED, error);
        }
    }

    fn render_rna_runs(&mut self, ui: &mut egui::Ui) {
        let mut cancel = None;
        let mut inspect = None;
        if let Some(error) = &self.rna_panel_ui.processing.error {
            ui.colored_label(egui::Color32::RED, error);
        }
        egui::ScrollArea::vertical().id_salt("rna_runs").show(ui, |ui| {
            for run in self.rna_panel_ui.processing.runs.iter().rev() {
                ui.push_id(run.id, |ui| {
                    ui.separator();
                    ui.heading(format!("Run {}", run.id));
                    ui.label(&run.request.output_dir);
                    if let Some(receipt) = &run.receipt {
                        ui.label(format!("{:?}: {}", receipt.state, receipt.phase));
                        if let Some(progress) = &receipt.rna_progress {
                            ui.label(format!("File {}/{} | {} | {} | {}", progress.file_ordinal, progress.files_total, progress.sample_id, progress.file_id, progress.phase));
                            ui.label(format!("{} reads in current file; compressed bytes {}/{} (decoder read-ahead, not verified completion)", progress.processed_reads, progress.compressed_bytes_read, progress.compressed_bytes_total));
                        }
                        if let Some(binding) = &receipt.artifact_result {
                            let outcome = ui.strong(format!("Published outcome: {}", binding.outcome));
                            capture_status(&outcome, &format!("rna_panel.run.{}.outcome", run.id), Some(&binding.report_content_sha256), &binding.outcome);
                            ui.label(&binding.report_path);
                            if button(ui, "Copy report binding", &format!("rna_panel.run.{}.binding", run.id)) {
                                if let Ok(json) = serde_json::to_string_pretty(binding) { ui.ctx().copy_text(json); }
                            }
                            if button(ui, "Open bound result", &format!("rna_panel.run.{}.inspect", run.id)) { inspect = Some(binding.clone()); }
                        }
                        if !terminal(receipt) && !run.cancel_pending && !run.cancel_requested
                            && button(ui, "Request cancellation", &format!("rna_panel.run.{}.cancel", run.id)) { cancel = Some(run.id); }
                    }
                    ui.label(&run.detail);
                });
            }
            ui.separator();
            ui.heading("Other shared-service activity (session-local)");
            for receipt in self.agent_command_service.receipts() {
                if receipt.rna_progress.is_some() || receipt.artifact_result.is_some() {
                    ui.label(format!("{}: {:?} / {}", receipt.job_id, receipt.state, receipt.phase));
                    if let Some(binding) = receipt.artifact_result {
                        if button(ui, "Open bound result", &format!("rna_panel.shared.{}.inspect", receipt.job_id)) { inspect = Some(binding); }
                    }
                }
            }
        });
        if let Some(binding) = inspect {
            self.rna_panel_ui
                .processing
                .evidence
                .published(binding.clone());
            if let Err(error) =
                self.start_rna_inspection(evidence::InspectionRequest::bound(binding))
            {
                self.rna_panel_ui.processing.evidence.error = Some(error);
            }
        }
        if let Some(id) = cancel {
            let request = self
                .rna_panel_ui
                .processing
                .runs
                .iter()
                .find(|r| r.id == id)
                .map(|r| r.request.clone());
            if let Some(request) = request {
                match self.submit_rna_operation(
                    Operation::CancelRnaReadPanel {
                        request: Box::new(request),
                    },
                    Task::Cancel(id),
                ) {
                    Ok(_) => {
                        if let Some(run) = self
                            .rna_panel_ui
                            .processing
                            .runs
                            .iter_mut()
                            .find(|r| r.id == id)
                        {
                            run.cancel_pending = true;
                        }
                    }
                    Err(error) => self.rna_panel_ui.processing.error = Some(error),
                }
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn rna_screen_review_is_explicit_immutable_and_edit_sensitive() {
        // Synthetic paths only: review must never open files or claim admission/checksum success.
        let mut state = Processing::default();
        assert!(state.review().is_err());
        state.draft.expected_panel_sha256 = Some("sha256:reviewed_synthetic_panel".into());
        state.draft.output_dir = "unopened-output".into();
        state.draft.cancel_path = Some("unopened-control".into());
        state.draft.files.push(RnaReadPanelInput {
            sample_id: "a".into(),
            barcode: "b".into(),
            priming: "custom_only".into(),
            path: "never-open.fastq.gz".into(),
            declared_complete: false,
            expected_md5: None,
            expected_sha256: None,
            checksum_manifest_key: None,
        });
        assert!(state.review().is_err());
        state.draft.files[0].declared_complete = true;
        state.review().unwrap();
        let reviewed = state.reviewed_request().unwrap();
        state.draft.files[0].sample_id = "changed".into();
        assert!(state.reviewed_request().is_err());
        assert_eq!(reviewed.files[0].sample_id, "a");
        state.draft.files[0].sample_id = "a".into();
        state.draft.policy.kmer_len += 2;
        assert!(state.reviewed_request().is_err());
    }

    #[test]
    fn rna_screen_window_close_does_not_cancel_or_discard_pending_work() {
        let mut app = GENtleApp::default();
        app.rna_panel_ui.open = true;
        app.rna_panel_ui.processing.pending.push(PendingTask {
            id: 42,
            task: Task::Screen,
        });
        assert!(app.rna_panel_pending());
        assert!(app.close_rna_panel_dialog());
        assert!(app.rna_panel_pending());
        assert_eq!(app.rna_panel_ui.processing.pending[0].id, 42);
    }

    #[cfg(feature = "gui-test-support")]
    #[test]
    fn rna_screen_native_semantic_targets_are_actual_controls() {
        let mut app = GENtleApp::default();
        app.rna_panel_ui.page = Page::Screen;
        let ctx = egui::Context::default();
        ctx.begin_pass(egui::RawInput {
            screen_rect: Some(egui::Rect::from_min_size(
                egui::Pos2::ZERO,
                Vec2::new(1200.0, 1600.0),
            )),
            ..Default::default()
        });
        crate::gui_test_support::begin_frame(&ctx);
        crate::egui_compat::show_central_panel_for_test_context(
            &ctx,
            egui::CentralPanel::default(),
            |ui| {
                app.render_rna_panel_contents(ui);
            },
        );
        let snapshot = crate::gui_test_support::snapshot(&ctx);
        for id in [
            "rna_panel.tab.screen",
            "rna_panel.tab.runs",
            "rna_panel.screen.load",
            "rna_panel.screen.inspect_panel",
            "rna_panel.screen.add_files",
            "rna_panel.screen.review",
        ] {
            assert!(
                snapshot.items.iter().any(|item| item.semantic_id == id),
                "{id}"
            );
        }
        crate::egui_compat::discard_test_pass_output(&ctx);
    }
}
