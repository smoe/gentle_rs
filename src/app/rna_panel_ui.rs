//! Native RNA-panel preparation, reviewed screening, monitoring and evidence inspection.
//!
//! Drafts are session-local. Preparation uses history-free project snapshots;
//! file-only processing uses the shared command service's independent policy.
//! Saved panels and earlier read-screen reports are never edited in place.

use super::*;
use crate::engine::{
    BuildRnaReadPanelRequest, RNA_PANEL_BUILD_SCHEMA, RnaReadPanel, RnaReadPanelGeneSetBinding,
    RnaReadPanelGeneSetChoice, RnaReadPanelSourceChoice,
};
use std::sync::mpsc::{self, Receiver, TryRecvError};

#[path = "rna_panel_ui/processing.rs"]
mod processing;

#[derive(Default, PartialEq, Eq, Clone, Copy)]
enum Page {
    #[default]
    Prepare,
    Screen,
    Runs,
    Inspect,
    Export,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
struct Binding {
    instance: u64,
    revision: u64,
}

impl Binding {
    fn new(engine: &GentleEngine) -> Self {
        Self {
            instance: engine.instance_id(),
            revision: engine.structural_revision(),
        }
    }
    fn matches(self, engine: &GentleEngine) -> bool {
        self == Self::new(engine)
    }
}

enum Action {
    Refresh,
    BindGeneSet(String),
    Open(String),
    Build(BuildRnaReadPanelRequest),
}
enum Outcome {
    Choices(
        Vec<RnaReadPanelSourceChoice>,
        Vec<RnaReadPanelGeneSetChoice>,
    ),
    Bound(RnaReadPanelGeneSetBinding, Vec<RnaReadPanelSourceChoice>),
    Opened(String, RnaReadPanel),
    Built(String, RnaReadPanel),
}
struct Pending {
    binding: Binding,
    commits: bool,
    receiver: Receiver<Result<Outcome, String>>,
}

#[derive(Default)]
pub(super) struct RnaPanelUiState {
    pub(super) open: bool,
    page: Page,
    processing: processing::Processing,
    choices: Vec<RnaReadPanelSourceChoice>,
    gene_sets: Vec<RnaReadPanelGeneSetChoice>,
    gene_set: Option<RnaReadPanelGeneSetBinding>,
    file_only: bool,
    rows: Vec<RnaReadPanelSourceChoice>,
    filter: String,
    source_label: String,
    assembly: String,
    release: String,
    base_path: Option<String>,
    base_count: usize,
    output_path: String,
    pending: Option<Pending>,
    panel: Option<RnaReadPanel>,
    status: String,
    error: Option<String>,
}

impl RnaPanelUiState {
    fn add(&mut self, choice: RnaReadPanelSourceChoice) -> Result<(), String> {
        if self.rows.iter().any(|r| {
            r.source.locus_id == choice.source.locus_id
                || (r.source.seq_id == choice.source.seq_id
                    && r.source.source_feature_id == choice.source.source_feature_id)
        }) {
            return Err("This source is already in the panel draft".into());
        }
        if self.rows.len() >= 128 {
            return Err("Panels admit at most 128 loci".into());
        }
        self.rows.push(choice);
        self.panel = None;
        Ok(())
    }

    fn opened(&mut self, path: String, panel: RnaReadPanel) {
        self.page = Page::Prepare;
        self.gene_set = panel
            .gene_set_resolution
            .as_deref()
            .and_then(|report| GentleEngine::rna_read_panel_gene_set_binding(report).ok());
        self.rows = panel
            .loci
            .iter()
            .map(|locus| RnaReadPanelSourceChoice {
                source: locus.source.clone(),
                transcript_ids: locus.annotated_transcript_ids.clone(),
                diagnostic: None,
            })
            .collect();
        self.base_count = self.rows.len();
        self.base_path = Some(path);
        self.output_path.clear();
        self.panel = Some(panel);
    }

    fn request(&self) -> Result<BuildRnaReadPanelRequest, String> {
        if self.output_path.trim().is_empty() {
            return Err("Choose a new panel output file".into());
        }
        if self.rows.is_empty() {
            return Err("Add at least one annotated locus".into());
        }
        if self.gene_set.is_none() {
            return Err(
                "Select an existing gene set; create or resolve it in Gene Set Inspector first"
                    .into(),
            );
        }
        Ok(BuildRnaReadPanelRequest {
            schema: RNA_PANEL_BUILD_SCHEMA.into(),
            sources: self.rows[self.base_count..]
                .iter()
                .map(|row| row.source.clone())
                .collect(),
            output_path: self.output_path.trim().into(),
            base_panel_path: self.base_path.clone(),
            gene_set: self.gene_set.clone(),
            register_in_project: Some(!self.file_only),
        })
    }

    fn recompose(&mut self) {
        self.base_path = None;
        self.base_count = 0;
        self.panel = None;
    }
}

fn execute(snapshot: &mut GentleEngine, action: Action) -> Result<Outcome, String> {
    match action {
        Action::Refresh | Action::BindGeneSet(_) => {
            let gene_set_id = match action {
                Action::BindGeneSet(id) => Some(id),
                _ => None,
            };
            let result = crate::engine_shell::execute_shell_command(
                snapshot,
                &ShellCommand::RnaReadPanelSources { gene_set_id },
            )?;
            let choices = serde_json::from_value(result.output["sources"].clone())
                .map_err(|e| e.to_string())?;
            if !result.output["selected_gene_set"].is_null() {
                let binding = serde_json::from_value(result.output["selected_gene_set"].clone())
                    .map_err(|e| e.to_string())?;
                Ok(Outcome::Bound(binding, choices))
            } else {
                let sets = serde_json::from_value(result.output["gene_sets"].clone())
                    .map_err(|e| e.to_string())?;
                Ok(Outcome::Choices(choices, sets))
            }
        }
        Action::Open(path) => {
            let (panel, _, _, _) = snapshot
                .inspect_rna_read_panel(&path, None, 0)
                .map_err(|e| e.to_string())?;
            Ok(Outcome::Opened(
                path,
                panel.ok_or("Choose a reference panel, not a screen or summary report")?,
            ))
        }
        Action::Build(request) => {
            let path = request.output_path.clone();
            let command = ShellCommand::Op {
                payload: serde_json::to_string(&Operation::BuildRnaReadPanel {
                    request: Box::new(request),
                })
                .map_err(|e| e.to_string())?,
            };
            let result = crate::engine_shell::execute_shell_command(snapshot, &command)?;
            let panel = serde_json::from_value(
                result
                    .output
                    .pointer("/result/rna_read_panel")
                    .ok_or("Shared build returned no panel")?
                    .clone(),
            )
            .map_err(|e| e.to_string())?;
            Ok(Outcome::Built(path, panel))
        }
    }
}

fn capture(_response: &egui::Response, _id: &str) {
    #[cfg(feature = "gui-test-support")]
    crate::gui_test_support::register_response_with_outcome(
        _response,
        _id,
        "window.rna_read_panel",
        None,
        if _id.ends_with(".filter") || _id.ends_with(".output") {
            crate::gui_test_support::GuiTestWidgetKind::TextInput
        } else if _id.starts_with("rna_panel.tab.") {
            crate::gui_test_support::GuiTestWidgetKind::Tab
        } else if _id.ends_with(".complete") {
            crate::gui_test_support::GuiTestWidgetKind::Checkbox
        } else {
            crate::gui_test_support::GuiTestWidgetKind::Button
        },
        false,
        None,
    );
}

fn capture_input(_response: &egui::Response, _id: String) {
    #[cfg(feature = "gui-test-support")]
    crate::gui_test_support::register_response_with_outcome(
        _response,
        _id,
        "window.rna_read_panel",
        None,
        crate::gui_test_support::GuiTestWidgetKind::TextInput,
        false,
        None,
    );
}

fn capture_status(_response: &egui::Response, _id: &str, _subject: Option<&str>, _outcome: &str) {
    #[cfg(feature = "gui-test-support")]
    crate::gui_test_support::register_response_with_outcome(
        _response,
        _id,
        "window.rna_read_panel",
        _subject,
        crate::gui_test_support::GuiTestWidgetKind::Status,
        false,
        Some(_outcome),
    );
}

fn capture_selection(_response: &egui::Response, _id: &str, _selected: bool) {
    #[cfg(feature = "gui-test-support")]
    crate::gui_test_support::register_response_with_outcome(
        _response,
        _id,
        "window.rna_read_panel",
        None,
        if _id.starts_with("rna_panel.tab.") {
            crate::gui_test_support::GuiTestWidgetKind::Tab
        } else {
            crate::gui_test_support::GuiTestWidgetKind::Checkbox
        },
        _selected,
        None,
    );
}

fn source_editor(ui: &mut egui::Ui, row: &mut RnaReadPanelSourceChoice) {
    let source = &mut row.source;
    ui.horizontal_wrapped(|ui| {
        ui.label(format!(
            "{} / feature {}",
            source.seq_id, source.source_feature_id
        ));
        ui.label(format!(
            "{} annotated transcripts",
            row.transcript_ids.len()
        ));
    });
    egui::Grid::new("identity").num_columns(2).show(ui, |ui| {
        for (label, value) in [
            ("Locus ID", &mut source.locus_id),
            ("Gene symbol", &mut source.gene_symbol),
            ("Stable gene ID", &mut source.gene_id),
            ("Source", &mut source.source),
            ("Assembly", &mut source.assembly),
            ("Annotation release", &mut source.annotation_release),
        ] {
            ui.label(label);
            let response = ui.add(egui::TextEdit::singleline(value).desired_width(240.0));
            capture_input(
                &response,
                format!(
                    "rna_panel.{}.{}.{}",
                    source.seq_id,
                    source.source_feature_id,
                    label.replace(' ', "_")
                ),
            );
            ui.end_row();
        }
    });
    let mut all = source.selected_transcript_ids.is_empty();
    if ui.checkbox(&mut all, "All annotated transcripts").changed() {
        source.selected_transcript_ids = if all {
            vec![]
        } else {
            row.transcript_ids.clone()
        };
    }
    if !all {
        ui.colored_label(
            egui::Color32::YELLOW,
            "Restricted competitor set; single-candidate support is relative to this selection.",
        );
        for id in &row.transcript_ids {
            let mut selected = source.selected_transcript_ids.contains(id);
            if ui
                .add_enabled(
                    !selected || source.selected_transcript_ids.len() > 1,
                    egui::Checkbox::new(&mut selected, id),
                )
                .on_hover_text("At least one transcript is required")
                .changed()
            {
                if selected {
                    source.selected_transcript_ids.push(id.clone());
                } else if source.selected_transcript_ids.len() > 1 {
                    source.selected_transcript_ids.retain(|s| s != id);
                }
                source.selected_transcript_ids.sort();
            }
        }
    }
    if let Some(error) = &row.diagnostic {
        ui.colored_label(egui::Color32::RED, error);
    }
}

impl GENtleApp {
    pub(super) fn append_rna_panel_lineage_rows(
        state: &crate::engine::ProjectState,
        rows: &mut Vec<LineageRow>,
        edges: &mut Vec<(String, String, String)>,
        labels: &mut HashMap<String, String>,
    ) {
        let entries = match GentleEngine::rna_read_panel_project_entries(state) {
            Ok(entries) => entries,
            Err(error) => {
                rows.push(rna_panel_lineage_row(
                    "analysis:rna_panel_error".into(),
                    "RNA panel store unavailable".into(),
                    error.to_string(),
                    None,
                    0,
                ));
                return;
            }
        };
        let present = entries
            .iter()
            .map(|e| e.panel.content_sha256.clone())
            .collect::<HashSet<_>>();
        let mut gene_sets = HashSet::new();
        for entry in entries {
            let panel = &entry.panel;
            let node = format!("analysis:rna_panel:{}", panel.content_sha256);
            let mut parents = vec![];
            if let Some(report) = &panel.gene_set_resolution {
                let id = GentleEngine::gene_set_resolution_artifact_id(report);
                let source = format!("gene_set:{id}");
                if gene_sets.insert(source.clone()) {
                    rows.push(rna_panel_lineage_row(
                        source.clone(),
                        "Target gene set".into(),
                        format!(
                            "{} resolved / {} unresolved: {}",
                            report.resolved_member_count,
                            report.unresolved_member_count,
                            report
                                .resolved_members
                                .iter()
                                .map(|m| m.symbol.as_str())
                                .collect::<Vec<_>>()
                                .join(", ")
                        ),
                        Some(id),
                        report.resolved_member_count,
                    ));
                }
                parents.push(source);
            }
            if let Some(parent) = &panel.parent_panel_sha256 {
                if present.contains(parent) {
                    parents.push(format!("analysis:rna_panel:{parent}"));
                }
            }
            for locus in &panel.loci {
                if let Some(parent) = state.lineage.seq_to_node.get(&locus.source.seq_id) {
                    parents.push(parent.clone());
                }
            }
            parents.sort();
            parents.dedup();
            let mut row = rna_panel_lineage_row(
                node.clone(),
                format!(
                    "RNA panel: {} genes / {} transcripts",
                    panel.gene_count, panel.selected_transcript_count
                ),
                "Reference definition, not read support".into(),
                Some(panel.content_sha256.clone()),
                panel.selected_transcript_count,
            );
            row.parents = parents.clone();
            row.created_by_op = entry.op_id.clone();
            row.created_at = entry.created_at_unix_ms;
            labels.insert(entry.op_id.clone(), "Build RNA reference panel".into());
            for parent in parents {
                edges.push((parent, node.clone(), entry.op_id.clone()));
            }
            rows.push(row);
        }
    }

    pub(super) fn open_registered_rna_panel(&mut self, id: &str) {
        let entry = self.engine.try_read().ok().and_then(|engine| {
            GentleEngine::rna_read_panel_project_entries(engine.state())
                .ok()?
                .into_iter()
                .find(|e| e.panel.content_sha256 == id)
        });
        if let Some(entry) = entry {
            self.rna_panel_ui.opened(entry.output_path, entry.panel);
            self.open_rna_panel_dialog();
        } else {
            self.app_status = "Saved RNA panel unavailable; inspect project evidence".into();
        }
    }

    pub(super) fn try_open_rna_panel_graph_node(&mut self, node: &str) -> bool {
        if let Some(id) = node.strip_prefix("analysis:rna_panel:") {
            self.open_registered_rna_panel(id);
            true
        } else if node.starts_with("gene_set:") {
            self.open_gene_set_inspector_dialog();
            true
        } else {
            node == "analysis:rna_panel_error"
        }
    }

    pub(super) fn rna_panel_viewport_id() -> ViewportId {
        ViewportId::from_hash_of("GENtle RNA Read Panel Viewport")
    }
    pub(super) fn hosted_rna_panel_window_id() -> egui::Id {
        egui::Id::new("hosted_rna_panel")
    }
    pub(super) fn rna_panel_pending(&self) -> bool {
        self.rna_panel_ui.pending.is_some() || self.rna_panel_ui.processing.pending()
    }
    pub(super) fn open_rna_panel_dialog(&mut self) {
        let was_open = self.rna_panel_ui.open;
        self.rna_panel_ui.open = true;
        self.mark_window_open_or_focus(Self::rna_panel_viewport_id(), was_open);
    }
    pub(super) fn close_rna_panel_dialog(&mut self) -> bool {
        let was_open = self.rna_panel_ui.open;
        self.rna_panel_ui.open = false;
        // Closing a window does not cancel an already approved artifact write.
        was_open
    }

    fn start_rna_panel_action(&mut self, ctx: &egui::Context, action: Action) {
        if self.rna_panel_ui.pending.is_some() {
            return;
        }
        let binding = match self.engine.try_read() {
            Ok(engine) => Binding::new(&engine),
            Err(_) => {
                self.rna_panel_ui.error = Some("Project is busy; retry when available".into());
                return;
            }
        };
        self.rna_panel_ui.error = None;
        self.rna_panel_ui.status = "Preparing panel".into();
        let shared = self.engine.clone();
        let commits = matches!(&action, Action::Build(request) if request.register_in_project.unwrap_or(false));
        let output_path = match &action {
            Action::Build(request) => Some(request.output_path.clone()),
            _ => None,
        };
        let repaint = ctx.clone();
        let (sender, receiver) = mpsc::channel();
        match std::thread::Builder::new()
            .name("rna-panel-editor".into())
            .stack_size(16 * 1024 * 1024)
            .spawn(move || {
                let result = if commits {
                    crate::background_engine::execute_on_engine_snapshot(&shared, |snapshot| {
                        if !binding.matches(snapshot) { return Err(crate::engine::EngineError::invalid_input("Project changed before panel build")); }
                        execute(snapshot, action).map_err(crate::engine::EngineError::invalid_input)
                    }).map_err(|error| format!("{error}; project attachment did not complete. A panel file may already exist at {}; inspect before retrying.", output_path.as_deref().unwrap_or("the requested path")))
                } else { (|| {
                    let mut snapshot = {
                        let engine = shared.read().map_err(|_| "Engine lock poisoned")?;
                        if !binding.matches(&engine) {
                            return Err("Project changed; refresh the panel sources".into());
                        }
                        engine.clone_without_history()
                    };
                    execute(&mut snapshot, action)
                })() };
                let _ = sender.send(result);
                repaint.request_repaint();
                repaint.request_repaint_of(ViewportId::ROOT);
            }) {
            Ok(_) => self.rna_panel_ui.pending = Some(Pending { binding, commits, receiver }),
            Err(error) => self.rna_panel_ui.error = Some(error.to_string()),
        }
    }

    pub(super) fn poll_rna_panel_action(&mut self) {
        self.poll_rna_processing_jobs();
        let Some(pending) = &self.rna_panel_ui.pending else {
            return;
        };
        let Ok(engine) = self.engine.try_read() else {
            return;
        };
        let current = pending.binding.matches(&engine)
            || (pending.commits && pending.binding.instance == engine.instance_id());
        drop(engine);
        let outcome = match pending.receiver.try_recv() {
            Ok(result) => result,
            Err(TryRecvError::Empty) => return,
            Err(TryRecvError::Disconnected) => Err("Panel worker stopped".into()),
        };
        self.rna_panel_ui.pending = None;
        let state = &mut self.rna_panel_ui;
        match outcome {
            Ok(Outcome::Built(path, panel)) => {
                state.status = format!(
                    "Saved {path}; {} genes, {} selected / {} annotated transcripts",
                    panel.gene_count,
                    panel.selected_transcript_count,
                    panel.annotated_transcript_count
                );
                if current {
                    state.opened(path, panel);
                } else {
                    state.error = Some(
                        "Saved from the earlier project snapshot; refresh before extending it"
                            .into(),
                    );
                }
            }
            Ok(_) if !current => {
                state.error = Some("Project changed; refresh or reopen the panel".into())
            }
            Ok(Outcome::Choices(choices, sets)) => {
                state.choices = choices;
                state.gene_sets = sets;
                state.status = format!("{} loaded source groups", state.choices.len());
            }
            Ok(Outcome::Bound(binding, sources)) => {
                state.gene_set = Some(binding);
                state.rows.truncate(state.base_count);
                for source in sources {
                    if !state
                        .rows
                        .iter()
                        .any(|row| row.source.gene_id == source.source.gene_id)
                    {
                        if let Err(error) = state.add(source) {
                            state.error = Some(error);
                        }
                    }
                }
                state.status =
                    "Reference bindings prepared from exact gene-set member IDs; review provenance"
                        .into();
            }
            Ok(Outcome::Opened(path, panel)) => {
                state.opened(path, panel);
                state.status = "Saved sources retained; add loci or explicitly recompose".into();
            }
            Err(error) => state.error = Some(error),
        }
    }

    fn render_rna_panel_contents(&mut self, ui: &mut egui::Ui) {
        #[cfg(feature = "gui-test-support")]
        crate::gui_test_support::register_rect(
            ui.ctx().clone(),
            "window.rna_read_panel",
            "window.rna_read_panel",
            None,
            crate::gui_test_support::GuiTestWidgetKind::Window,
            ui.max_rect(),
            true,
            true,
            false,
            None,
        );
        let mut action = None;
        let mut retrieve = false;
        let mut inspect_sets = false;
        let busy = self.rna_panel_ui.pending.is_some();
        ui.heading("RNA Read Panel");
        ui.horizontal_wrapped(|ui| {
            for (page, label, id) in [
                (Page::Prepare, "Prepare", "rna_panel.tab.prepare"),
                (Page::Screen, "Screen", "rna_panel.tab.screen"),
                (Page::Runs, "Runs", "rna_panel.tab.runs"),
                (Page::Inspect, "Inspect", "rna_panel.tab.inspect"),
                (Page::Export, "Export", "rna_panel.tab.export"),
            ] {
                let response = ui.selectable_value(&mut self.rna_panel_ui.page, page, label);
                capture_selection(&response, id, self.rna_panel_ui.page == page);
            }
        });
        self.render_rna_rehearsal_controls(ui);
        if self.rna_panel_ui.page != Page::Prepare {
            self.render_rna_processing(ui);
            return;
        }
        ui.horizontal_wrapped(|ui| {
            let refresh = ui.add_enabled(!busy, egui::Button::new("Refresh loaded loci"));
            capture(&refresh, "rna_panel.refresh");
            if refresh.clicked() {
                action = Some(Action::Refresh);
            }
            let open = ui.add_enabled(!busy, egui::Button::new("Open panel..."));
            capture(&open, "rna_panel.open");
            if open.clicked()
                && let Some(path) = rfd::FileDialog::new()
                    .add_filter("RNA reference panel", &["json"])
                    .pick_file()
            {
                action = Some(Action::Open(path.display().to_string()));
            }
            if ui.button("Retrieve gene...").clicked() {
                retrieve = true;
            }
            if ui.button("Gene sets...").clicked() {
                inspect_sets = true;
            }
            if busy {
                ui.spinner();
            }
        });
        let state = &mut self.rna_panel_ui;
        ui.label(&state.status);
        if let Some(error) = &state.error {
            ui.colored_label(egui::Color32::RED, error);
        }
        egui::ScrollArea::vertical()
            .id_salt("rna_panel_scroll")
            .show(ui, |ui| {
                ui.add_enabled_ui(!busy, |ui| {
                    ui.heading("Target gene set");
                    egui::ComboBox::from_id_salt("rna_panel_gene_set")
                        .selected_text(
                            state
                                .gene_set
                                .as_ref()
                                .map(|b| b.resolution_id.as_str())
                                .unwrap_or("Choose stored gene set"),
                        )
                        .show_ui(ui, |ui| {
                            for choice in &state.gene_sets {
                                ui.selectable_value(
                                    &mut state.gene_set,
                                    Some(choice.binding.clone()),
                                    format!(
                                        "{} / {} resolved, {} unresolved",
                                        choice.binding.resolution_id,
                                        choice.resolution.resolved_member_count,
                                        choice.resolution.unresolved_member_count
                                    ),
                                );
                            }
                        });
                    if let Some(binding) = &state.gene_set {
                        if ui.button("Bind loaded gene-set members").clicked() {
                            action = Some(Action::BindGeneSet(binding.resolution_id.clone()));
                        }
                        if let Some(choice) = state.gene_sets.iter().find(|c| &c.binding == binding)
                        {
                            for member in &choice.resolution.resolved_members {
                                ui.label(format!(
                                    "{} / {}",
                                    member.symbol,
                                    member.gene_id.as_deref().unwrap_or("stable ID unassessed")
                                ));
                            }
                            for member in &choice.resolution.unresolved_members {
                                ui.colored_label(
                                    egui::Color32::RED,
                                    format!("{}: {}", member.query, member.reason),
                                );
                            }
                        }
                    }
                    ui.heading("Reference and transcript bindings");
                    ui.horizontal_wrapped(|ui| {
                        ui.label("Find loaded locus");
                        let input = ui.text_edit_singleline(&mut state.filter);
                        capture(&input, "rna_panel.filter");
                    });
                    let filter = state.filter.to_lowercase();
                    let mut add = None;
                    let picker = egui::ComboBox::from_id_salt("rna_panel_source_picker")
                        .selected_text("Add loaded locus")
                        .width(320.0)
                        .show_ui(ui, |ui| {
                            for choice in &state.choices {
                                let label = format!(
                                    "{} | {} | feature {} | {} transcripts",
                                    choice.source.gene_symbol,
                                    choice.source.seq_id,
                                    choice.source.source_feature_id,
                                    choice.transcript_ids.len()
                                );
                                if !label.to_lowercase().contains(&filter) {
                                    continue;
                                }
                                let button = ui
                                    .add_enabled(
                                        choice.diagnostic.is_none(),
                                        egui::Button::new(&label),
                                    )
                                    .on_hover_text(
                                        choice
                                            .diagnostic
                                            .as_deref()
                                            .unwrap_or(&choice.source.gene_id),
                                    );
                                capture(
                                    &button,
                                    &format!(
                                        "rna_panel.add.{}.{}",
                                        choice.source.seq_id, choice.source.source_feature_id
                                    ),
                                );
                                if button.clicked() {
                                    add = Some(choice.clone());
                                }
                            }
                        });
                    capture(&picker.response, "rna_panel.add_menu");
                    if let Some(choice) = add {
                        if let Err(error) = state.add(choice) {
                            state.error = Some(error);
                        }
                    }
                    for choice in state.choices.iter().filter(|c| c.diagnostic.is_some()) {
                        ui.colored_label(
                            egui::Color32::RED,
                            format!(
                                "{} / feature {}: {}",
                                choice.source.seq_id,
                                choice.source.source_feature_id,
                                choice.diagnostic.as_deref().unwrap_or_default()
                            ),
                        );
                    }
                    ui.collapsing("Declared provenance defaults", |ui| {
                        for (label, value) in [
                            ("Source", &mut state.source_label),
                            ("Assembly", &mut state.assembly),
                            ("Annotation release", &mut state.release),
                        ] {
                            ui.horizontal_wrapped(|ui| {
                                ui.label(label);
                                ui.text_edit_singleline(value);
                            });
                        }
                        if ui.button("Fill missing draft provenance").clicked() {
                            for row in &mut state.rows[state.base_count..] {
                                for (field, value) in [
                                    (&mut row.source.source, &state.source_label),
                                    (&mut row.source.assembly, &state.assembly),
                                    (&mut row.source.annotation_release, &state.release),
                                ] {
                                    if field.is_empty() {
                                        *field = value.trim().into();
                                    }
                                }
                            }
                        }
                    });
                    if let Some(path) = &state.base_path {
                        ui.label(format!("Base panel: {path}"));
                        if ui.button("Recompose all sources").clicked() {
                            state.recompose();
                        }
                    }
                    let mut remove = None;
                    for (i, row) in state.rows.iter_mut().enumerate() {
                        let label = format!("{} / {}", row.source.gene_symbol, row.source.seq_id);
                        ui.push_id(i, |ui| {
                            ui.separator();
                            ui.horizontal_wrapped(|ui| {
                                ui.strong(&label);
                                if i < state.base_count {
                                    ui.label("Retained from saved panel");
                                } else if ui
                                    .small_button("-")
                                    .on_hover_text("Remove locus from draft")
                                    .clicked()
                                {
                                    remove = Some(i);
                                }
                            });
                            ui.collapsing("Identity and transcript selection", |ui| {
                                ui.add_enabled_ui(i >= state.base_count, |ui| {
                                    source_editor(ui, row)
                                });
                            });
                        });
                    }
                    if let Some(i) = remove {
                        state.rows.remove(i);
                        state.panel = None;
                    }
                    ui.separator();
                    ui.label(format!("{} loci in draft", state.rows.len()));
                    ui.horizontal_wrapped(|ui| {
                        ui.label("New panel file");
                        let path = ui.text_edit_singleline(&mut state.output_path);
                        capture(&path, "rna_panel.output");
                        if ui.button("Choose...").clicked()
                            && let Some(path) = rfd::FileDialog::new()
                                .add_filter("RNA reference panel", &["json"])
                                .save_file()
                        {
                            state.output_path = path.display().to_string();
                        }
                    });
                    let mut keep_in_graph = !state.file_only;
                    if ui
                        .checkbox(&mut keep_in_graph, "Keep panel in project graph")
                        .changed()
                    {
                        state.file_only = !keep_in_graph;
                    }
                    let build = ui
                        .add_enabled(!state.rows.is_empty(), egui::Button::new("Build new panel"));
                    capture(&build, "rna_panel.build");
                    if build.clicked() {
                        match state.request() {
                            Ok(request) => action = Some(Action::Build(request)),
                            Err(error) => state.error = Some(error),
                        }
                    }
                    if ui.button("Copy build request").clicked() {
                        match state.request().and_then(|request| {
                            serde_json::to_string_pretty(&request).map_err(|e| e.to_string())
                        }) {
                            Ok(json) => ui.ctx().copy_text(json),
                            Err(error) => state.error = Some(error),
                        }
                    }
                    if let Some(panel) = &state.panel {
                        ui.monospace(&panel.content_sha256);
                        if !panel.overlapping_source_pairs.is_empty() {
                            ui.colored_label(
                                egui::Color32::RED,
                                format!(
                                    "{} overlapping source pairs; screening is blocked",
                                    panel.overlapping_source_pairs.len()
                                ),
                            );
                        }
                    }
                });
            });
        if retrieve {
            self.open_reference_genome_retrieve_dialog();
        }
        if inspect_sets {
            self.open_gene_set_inspector_dialog();
        }
        if let Some(action) = action {
            self.start_rna_panel_action(ui.ctx(), action);
        }
    }

    pub(super) fn render_rna_panel_dialog(&mut self, ctx: &egui::Context) {
        if !self.rna_panel_ui.open {
            return;
        }
        let mut open = true;
        let viewport_id = Self::rna_panel_viewport_id();
        let spec = self.hosted_window_spec_for_viewport(
            "RNA Read Panel",
            Self::hosted_rna_panel_window_id(),
            viewport_id,
            Vec2::new(860.0, 720.0),
            Vec2::new(420.0, 320.0),
        );
        if ctx.embed_viewports() {
            crate::egui_compat::show_hosted_window(ctx, &spec, &mut open, |ui| {
                self.render_rna_panel_contents(ui)
            });
            self.clear_viewport_foreground_request_after_render(viewport_id);
            self.finalize_viewport_open_probe(viewport_id, "RNA Read Panel");
        } else {
            ctx.show_viewport_immediate(
                viewport_id,
                crate::egui_compat::viewport_builder_for_hosted_window(&spec),
                |ctx, class| {
                    self.note_viewport_focus_if_active(ctx, viewport_id);
                    if class == egui::ViewportClass::EmbeddedWindow {
                        crate::egui_compat::show_hosted_window(&mut *ctx, &spec, &mut open, |ui| {
                            self.render_rna_panel_contents(ui)
                        });
                    } else {
                        crate::egui_compat::show_central_panel(
                            &mut *ctx,
                            egui::CentralPanel::default(),
                            |ui| self.render_rna_panel_contents(ui),
                        );
                        if Self::viewport_close_requested_or_shortcut(ctx) {
                            open = false;
                        }
                    }
                },
            );
        }
        if !open {
            self.close_rna_panel_dialog();
        }
    }
}

fn rna_panel_lineage_row(
    node_id: String,
    display_name: String,
    status: String,
    artifact: Option<String>,
    count: usize,
) -> LineageRow {
    LineageRow {
        kind: LineageNodeKind::Analysis,
        seq_id: node_id.clone(),
        node_id,
        display_name,
        origin: "RnaReadPanelDefinition".into(),
        created_by_op: "-".into(),
        created_at: 0,
        parents: vec![],
        length: 0,
        circular: false,
        pool_size: 0,
        pool_members: vec![],
        arrangement_id: None,
        arrangement_mode: None,
        lane_container_ids: vec![],
        ladders: vec![],
        genome_anchor_summary: None,
        genome_anchor_display: None,
        is_full_genome_sequence: false,
        retrieval_descriptor: None,
        analysis_kind: None,
        analysis_artifact_id: artifact,
        analysis_reference_seq_id: None,
        analysis_mode: Some("rna_panel_definition".into()),
        analysis_status: Some(status),
        analysis_point_count: None,
        analysis_bin_count: None,
        analysis_read_count: None,
        analysis_trace_count: None,
        analysis_target_count: Some(count),
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

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{
        dna_sequence::DNAsequence,
        engine::{GeneSetRequest, GeneSetResolutionReport, GeneSetResolvedMember, ProjectState},
    };
    use std::path::Path;

    fn synthetic_engine() -> GentleEngine {
        // Inline hand-crafted two-exon toy records; not human genes or experimental evidence.
        let mut state = ProjectState::default();
        for number in 0..3 {
            let mut dna = DNAsequence::from_sequence(&"ACGT".repeat(40)).unwrap();
            for transcript in 0..2 {
                dna.features_mut().push(gb_io::seq::Feature {
                    kind: "mRNA".into(),
                    location: gb_io::seq::Location::Join(vec![
                        gb_io::seq::Location::simple_range(0, 40),
                        gb_io::seq::Location::simple_range(
                            80 + 10 * transcript,
                            120 + 10 * transcript,
                        ),
                    ]),
                    qualifiers: vec![
                        ("gene".into(), Some(format!("TOY_{number}"))),
                        ("gene_id".into(), Some(format!("SYN_{number}"))),
                        ("transcript_id".into(), Some(format!("tx_{transcript}"))),
                    ],
                });
            }
            state.sequences.insert(format!("toy_{number}"), dna);
        }
        let mut engine = GentleEngine::from_state(state);
        for count in 1..=3 {
            let report = GeneSetResolutionReport {
                schema: "gentle.gene_set_resolution.v1".into(),
                generated_at_unix_ms: count as u128,
                op_id: Some(format!("synthetic_set_{count}")),
                request: GeneSetRequest::ExplicitMembers {
                    members: (0..count).map(|n| format!("TOY_{n}")).collect(),
                },
                requested_member_count: count,
                resolved_member_count: count,
                resolved_members: (0..count)
                    .map(|n| GeneSetResolvedMember {
                        dedup_key: format!("SYN_{n}"),
                        symbol: format!("TOY_{n}"),
                        gene_id: Some(format!("SYN_{n}")),
                        ..Default::default()
                    })
                    .collect(),
                ..Default::default()
            };
            engine.upsert_gene_set_resolution_artifact(report).unwrap();
        }
        engine
    }

    fn set_binding(engine: &GentleEngine, count: usize) -> RnaReadPanelGeneSetBinding {
        engine
            .rna_read_panel_gene_sets()
            .unwrap()
            .into_iter()
            .find(|c| c.resolution.resolved_member_count == count)
            .unwrap()
            .binding
    }

    fn reviewed(mut choice: RnaReadPanelSourceChoice) -> RnaReadPanelSourceChoice {
        choice.source.source = "hand-crafted synthetic".into();
        choice.source.assembly = "synthetic".into();
        choice.source.annotation_release = "v1".into();
        choice
    }

    #[test]
    fn rna_panel_editor_discovery_and_extension_match_shared_shell() {
        let td = tempfile::tempdir().unwrap();
        let mut engine = synthetic_engine();
        let before = serde_json::to_value(engine.state()).unwrap();
        let choices = match execute(&mut engine, Action::Refresh).unwrap() {
            Outcome::Choices(choices, sets) => {
                assert_eq!(sets.len(), 3);
                choices
            }
            _ => panic!("expected source choices"),
        };
        assert_eq!(choices.len(), 3);
        let set_id = set_binding(&engine, 1).resolution_id;
        match execute(&mut engine, Action::BindGeneSet(set_id)).unwrap() {
            Outcome::Bound(binding, sources) => {
                assert_eq!(binding, set_binding(&engine, 1));
                assert_eq!(sources.len(), 1);
                assert_eq!(sources[0].source.gene_id, "SYN_0");
            }
            _ => panic!("expected exact gene-set source binding"),
        }
        let mut draft = RnaPanelUiState::default();
        draft.gene_set = Some(set_binding(&engine, 1));
        draft.file_only = true;
        draft.add(reviewed(choices[0].clone())).unwrap();
        assert!(draft.add(reviewed(choices[0].clone())).is_err());
        draft.output_path = td.path().join("first panel.json").display().to_string();
        let (path, base) =
            match execute(&mut engine, Action::Build(draft.request().unwrap())).unwrap() {
                Outcome::Built(path, panel) => (path, panel),
                _ => panic!("expected panel"),
            };
        let old_bytes = std::fs::read(&path).unwrap();
        draft.opened(path.clone(), base.clone());
        draft.gene_set = Some(set_binding(&engine, 2));
        draft.add(reviewed(choices[1].clone())).unwrap();
        draft.output_path = td.path().join("extended panel.json").display().to_string();
        let request = draft.request().unwrap();
        assert_eq!(request.base_panel_path.as_deref(), Some(path.as_str()));
        assert_eq!(
            request.sources.len(),
            1,
            "retained sources are supplied by the bound base"
        );
        let panel = match execute(&mut engine, Action::Build(request)).unwrap() {
            Outcome::Built(_, panel) => panel,
            _ => panic!("expected extension"),
        };
        assert_eq!((panel.gene_count, panel.selected_transcript_count), (2, 4));
        assert_ne!(panel.content_sha256, base.content_sha256);
        assert_eq!(std::fs::read(&path).unwrap(), old_bytes);
        assert_eq!(serde_json::to_value(engine.state()).unwrap(), before);
        draft.recompose();
        assert!(draft.request().unwrap().base_panel_path.is_none());
        assert_eq!(draft.request().unwrap().sources.len(), 2);
    }

    #[test]
    fn rna_panel_editor_ui_intent_opens_focuses_and_closes_without_writing() {
        let mut app = GENtleApp::default();
        let before = serde_json::to_value(&*app.engine.read().unwrap()).unwrap();
        for action in ["open", "focus", "close"] {
            let command =
                crate::engine_shell::parse_shell_line(&format!("ui {action} rna-read-panel"))
                    .unwrap();
            assert!(app.try_apply_shell_ui_intent(&command).is_some());
            assert_eq!(app.rna_panel_ui.open, action != "close");
        }
        assert_eq!(
            serde_json::to_value(&*app.engine.read().unwrap()).unwrap(),
            before
        );
        assert!(!app.close_rna_panel_dialog());
    }

    #[test]
    fn rna_panel_gene_set_registration_graph_save_reload_and_undo_agree() {
        let td = tempfile::tempdir().unwrap();
        let mut engine = synthetic_engine();
        let binding = set_binding(&engine, 1);
        let (bound, sources) = engine
            .rna_read_panel_gene_set_sources(&binding.resolution_id)
            .unwrap();
        assert_eq!(bound, binding);
        let mut draft = RnaPanelUiState::default();
        draft.gene_set = Some(binding.clone());
        draft.add(reviewed(sources[0].clone())).unwrap();
        draft.output_path = td
            .path()
            .join("registered panel.json")
            .display()
            .to_string();
        let base = match execute(&mut engine, Action::Build(draft.request().unwrap())).unwrap() {
            Outcome::Built(_, panel) => panel,
            _ => panic!("expected panel"),
        };
        assert_eq!(
            GentleEngine::rna_read_panel_project_entries(engine.state())
                .unwrap()
                .len(),
            1
        );
        draft.opened(draft.output_path.clone(), base.clone());
        draft.gene_set = Some(set_binding(&engine, 2));
        let (_, sources) = engine
            .rna_read_panel_gene_set_sources(&draft.gene_set.as_ref().unwrap().resolution_id)
            .unwrap();
        draft.add(reviewed(sources[1].clone())).unwrap();
        draft.output_path = td
            .path()
            .join("registered extension.json")
            .display()
            .to_string();
        let extended = match execute(&mut engine, Action::Build(draft.request().unwrap())).unwrap()
        {
            Outcome::Built(_, panel) => panel,
            _ => panic!("expected panel"),
        };
        let reloaded = GentleEngine::from_state(
            serde_json::from_value(serde_json::to_value(engine.state()).unwrap()).unwrap(),
        );
        assert_eq!(
            GentleEngine::rna_read_panel_project_entries(reloaded.state())
                .unwrap()
                .len(),
            2
        );
        let mut rows = vec![];
        let mut edges = vec![];
        GENtleApp::append_rna_panel_lineage_rows(
            reloaded.state(),
            &mut rows,
            &mut edges,
            &mut HashMap::new(),
        );
        let base_node = format!("analysis:rna_panel:{}", base.content_sha256);
        let child_node = format!("analysis:rna_panel:{}", extended.content_sha256);
        let set_node = format!("gene_set:{}", binding.resolution_id);
        assert_eq!(rows.iter().filter(|r| r.node_id == base_node).count(), 1);
        assert!(
            edges
                .iter()
                .any(|(from, to, _)| from == &set_node && to == &base_node)
        );
        assert!(
            edges
                .iter()
                .any(|(from, to, _)| from == &base_node && to == &child_node)
        );
        let (svg_nodes, svg_edges) = crate::lineage_export::build_lineage_svg_graph(
            reloaded.state(),
            reloaded.operation_log(),
        );
        for row in &rows {
            assert!(svg_nodes.iter().any(|n| n.node_id == row.node_id));
        }
        for (from, to, _) in &edges {
            assert!(
                svg_edges
                    .iter()
                    .any(|e| &e.from_node_id == from && &e.to_node_id == to)
            );
        }
        assert!(
            svg_nodes
                .iter()
                .any(|n| n.node_id == child_node && n.subtitle.contains("reference definition"))
        );
        let mut app = GENtleApp::default();
        *app.engine.write().unwrap() = reloaded;
        assert!(app.try_open_rna_panel_graph_node(&child_node));
        assert_eq!(
            app.rna_panel_ui.panel.as_ref().unwrap().content_sha256,
            extended.content_sha256
        );
        assert!(!app.try_open_rna_panel_graph_node("ordinary_sequence"));
        engine.undo_last_operation().unwrap();
        assert_eq!(
            GentleEngine::rna_read_panel_project_entries(engine.state())
                .unwrap()
                .len(),
            1
        );
        assert!(
            Path::new(&draft.output_path).exists(),
            "undo changes project state, not external files"
        );
        engine.redo_last_operation().unwrap();
        assert_eq!(
            GentleEngine::rna_read_panel_project_entries(engine.state())
                .unwrap()
                .len(),
            2
        );
    }

    #[test]
    fn rna_panel_gene_set_binding_rejects_stale_unresolved_and_relabelled_sources() {
        let td = tempfile::tempdir().unwrap();
        let mut engine = synthetic_engine();
        let binding = set_binding(&engine, 1);
        let mut draft = RnaPanelUiState::default();
        draft.gene_set = Some(binding.clone());
        let choices = engine.rna_read_panel_source_choices().unwrap();
        // A label is not a stable-ID relation: source 1 must not impersonate set member 0.
        let mut wrong = reviewed(choices[1].clone());
        wrong.source.gene_id = "SYN_0".into();
        draft.add(wrong).unwrap();
        draft.output_path = td.path().join("wrong.json").display().to_string();
        assert!(
            execute(&mut engine, Action::Build(draft.request().unwrap()))
                .err()
                .unwrap()
                .contains("panel_gene_set_source_id_unbound")
        );
        assert!(!Path::new(&draft.output_path).exists());
        draft.rows = vec![reviewed(choices[0].clone())];
        draft.gene_set.as_mut().unwrap().resolution_sha256 = "sha256:stale".into();
        assert!(
            execute(&mut engine, Action::Build(draft.request().unwrap()))
                .err()
                .unwrap()
                .contains("panel_gene_set_changed")
        );
        let mut report = engine
            .get_gene_set_resolution_artifact(&binding.resolution_id)
            .unwrap();
        report.unresolved_member_count = 1;
        report
            .unresolved_members
            .push(crate::engine::GeneSetUnresolvedMember {
                query: "UNKNOWN".into(),
                reason: "synthetic unresolved target".into(),
                ..Default::default()
            });
        engine.upsert_gene_set_resolution_artifact(report).unwrap();
        assert!(
            engine
                .rna_read_panel_gene_set_sources(&binding.resolution_id)
                .unwrap_err()
                .message
                .contains("unresolved")
        );
        draft.gene_set = Some(set_binding(&engine, 1));
        assert!(execute(&mut engine, Action::Build(draft.request().unwrap())).is_err());
        assert!(!Path::new(&draft.output_path).exists());
    }

    #[test]
    fn rna_panel_gene_set_mapping_refuses_ambiguous_sources_and_unknown_reference() {
        let td = tempfile::tempdir().unwrap();
        let mut engine = synthetic_engine();
        let binding = set_binding(&engine, 1);
        let duplicate = engine.state().sequences["toy_0"].clone();
        engine
            .state_mut()
            .sequences
            .insert("duplicate_view".into(), duplicate);
        assert!(
            engine
                .rna_read_panel_gene_set_sources(&binding.resolution_id)
                .unwrap_err()
                .message
                .contains("ambiguous_or_missing_source")
        );
        engine.state_mut().sequences.remove("duplicate_view");
        let mut report = engine
            .get_gene_set_resolution_artifact(&binding.resolution_id)
            .unwrap();
        report.genome_id = Some("synthetic_reference_without_a_bound_anchor".into());
        report.biological_contexts = Default::default();
        engine.upsert_gene_set_resolution_artifact(report).unwrap();
        let mut draft = RnaPanelUiState::default();
        draft.gene_set = Some(set_binding(&engine, 1));
        let source = engine.rna_read_panel_source_choices().unwrap().remove(0);
        draft.add(reviewed(source)).unwrap();
        draft.output_path = td
            .path()
            .join("unknown reference.json")
            .display()
            .to_string();
        assert!(
            execute(&mut engine, Action::Build(draft.request().unwrap()))
                .err()
                .unwrap()
                .contains("reference_mismatch_or_unassessed")
        );
        assert!(!Path::new(&draft.output_path).exists());
    }

    #[test]
    fn rna_panel_editor_drops_late_inventory_after_project_replacement() {
        let mut app = GENtleApp::default();
        let binding = Binding::new(&app.engine.read().unwrap());
        let (sender, receiver) = mpsc::channel();
        sender
            .send(Ok(Outcome::Choices(
                synthetic_engine().rna_read_panel_source_choices().unwrap(),
                vec![],
            )))
            .unwrap();
        app.rna_panel_ui.pending = Some(Pending {
            binding,
            commits: false,
            receiver,
        });
        *app.engine.write().unwrap() = GentleEngine::default();
        app.poll_rna_panel_action();
        assert!(app.rna_panel_ui.choices.is_empty());
        assert!(
            app.rna_panel_ui
                .error
                .as_deref()
                .unwrap()
                .contains("Project changed")
        );
    }
}
