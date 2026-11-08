//! Thin reviewed requests and virtualized shared-renderer lanes; no file I/O or biological projection.

use super::*;
use crate::engine_shell::quote_shell_arg;
use crate::observed_transcript_structures::read_ends::pileup::{
    self, Channels, Filter, SortKey, Status, View,
};
use crate::rna_read_ends::{Command, Format};

pub(super) struct UiState {
    draft: String,
    reviewed: Option<Command>,
    approved: bool,
    report_json: String,
    export_dir: String,
    view: View,
    view_dirty: bool,
    status_filter: usize,
    active: Option<u64>,
    active_sets_report: bool,
    active_inspection_binding: Option<serde_json::Value>,
    status: String,
    error: Option<String>,
    page: Option<pileup::Page>,
    scene: Option<gentle_render::read_end_pileup::Scene>,
    lane_heights: Vec<u32>,
    max_lane_height: f32,
    model_height: f32,
    request_text: String,
}
impl Default for UiState {
    fn default() -> Self {
        Self {
            draft: String::new(),
            reviewed: None,
            approved: false,
            report_json: String::new(),
            export_dir: String::new(),
            view: View {
                schema: pileup::VIEW_SCHEMA.into(),
                sort: SortKey::UnitKey,
                descending: false,
                filters: vec![],
                offset: 0,
                count: 100,
                max_sublanes: 2000,
                width: 1400,
                scale: "genomic".into(),
                channels: Channels {
                    utr5: true,
                    utr3: true,
                    read_length: true,
                    clips: true,
                    composition: true,
                    tail: true,
                },
            },
            view_dirty: false,
            status_filter: 0,
            active: None,
            active_sets_report: false,
            active_inspection_binding: None,
            status: String::new(),
            error: None,
            page: None,
            scene: None,
            lane_heights: vec![],
            max_lane_height: 54.0,
            model_height: 54.0,
            request_text: String::new(),
        }
    }
}

fn parse_draft(text: &str) -> Result<Command, String> {
    if text.len() > 1024 * 1024 {
        return Err("Use a bounded command draft (at most 1 MiB)".into());
    }
    serde_json::from_str(text).map_err(|e| e.to_string())
}

fn matches_draft(state: &UiState, view_sha256: &str) -> bool {
    pileup::hash(&state.view).ok().as_deref() == Some(view_sha256)
        && serde_json::from_str::<serde_json::Value>(&state.report_json).ok()
            == state.active_inspection_binding
}

impl GENtleApp {
    fn submit_read_end(&mut self, command: Command) -> Result<(), String> {
        let sets_report = matches!(&command, Command::Build { .. });
        let inspection = if let Command::Inspect { report, view } = &command {
            Some((
                serde_json::to_value(report).map_err(|e| e.to_string())?,
                view.clone(),
            ))
        } else {
            None
        };
        let payload = serde_json::to_string(&Operation::RnaReadEnds {
            command: Box::new(command),
        })
        .map_err(|e| e.to_string())?;
        let text = format!("op {}", quote_shell_arg(&payload));
        let id = self.agent_command_service.submit_parsed(
            self.engine.clone(),
            ShellCommand::Op { payload },
            text.clone(),
            ShellExecutionOptions::default(),
        )?;
        let state = &mut self.rna_panel_ui.read_ends;
        state.active = Some(id);
        state.active_sets_report = sets_report;
        state.active_inspection_binding = inspection.as_ref().map(|(binding, _)| binding.clone());
        if let Some((binding, view)) = inspection {
            state.report_json = serde_json::to_string(&binding).map_err(|e| e.to_string())?;
            state.view = view;
            state.view_dirty = true;
        }
        state.request_text = text;
        state.approved = false;
        state.status = "Admitted; awaiting an actual terminal receipt".into();
        state.error = None;
        Ok(())
    }

    pub(super) fn poll_read_end_jobs(&mut self) {
        let Some(id) = self.rna_panel_ui.read_ends.active else {
            return;
        };
        let receipt = self.agent_command_service.status(id);
        let state = &mut self.rna_panel_ui.read_ends;
        if let Some(receipt) = &receipt {
            state.status = format!(
                "Job {}: {:?}, {}; cancellation requested={}, publication={}",
                id, receipt.state, receipt.phase, receipt.cancel_requested, receipt.publication
            );
        }
        let Some(result) = self.agent_command_service.take_result(id) else {
            return;
        };
        state.active = None;
        match result {
            Err(error) => state.error = Some(error),
            Ok(run) => {
                let value = &run.output["result"]["rna_read_ends"];
                if let Some(page) = value.get("page") {
                    let loaded = (|| {
                        let page: pileup::Page =
                            serde_json::from_value(page.clone()).map_err(|e| e.to_string())?;
                        let scene: gentle_render::read_end_pileup::Scene =
                            serde_json::from_value(value["scene"].clone())
                                .map_err(|e| e.to_string())?;
                        let heights = gentle_render::read_end_pileup::lane_heights(&scene)?;
                        let model_height = gentle_render::read_end_pileup::model_height(&scene)?;
                        Ok::<_, String>((page, scene, heights, model_height))
                    })();
                    match loaded {
                        Ok((page, scene, heights, model_height)) => {
                            state.view_dirty = !matches_draft(state, &page.view_sha256);
                            state.max_lane_height =
                                heights.iter().copied().max().unwrap_or(54) as f32;
                            state.lane_heights = heights;
                            state.model_height = model_height as f32;
                            state.page = Some(page);
                            state.scene = Some(scene);
                        }
                        Err(error) => state.error = Some(error),
                    }
                }
                if state.active_sets_report
                    && let Some(report) = value.get("report")
                {
                    state.report_json = serde_json::to_string(report).unwrap_or_default();
                    state.page = None;
                    state.scene = None;
                    state.view_dirty = true;
                }
                if state.active_sets_report
                    && let Some(completion) = receipt.and_then(|r| r.read_end_result)
                {
                    match completion.file_binding().and_then(|b| {
                        serde_json::to_string(&b)
                            .map_err(|e| crate::engine::EngineError::invalid_input(e.to_string()))
                    }) {
                        Ok(binding) => {
                            state.report_json = binding;
                            state.page = None;
                            state.scene = None;
                            state.view_dirty = true;
                        }
                        Err(error) => state.error = Some(error.to_string()),
                    }
                }
            }
        }
        state.active_sets_report = false;
        state.active_inspection_binding = None;
    }

    pub(super) fn render_read_end_pileup(&mut self, ui: &mut egui::Ui) {
        let mut action = None;
        let state = &mut self.rna_panel_ui.read_ends;
        ui.heading("Saved read ends and transcript pileup");
        ui.label("Explicit source / gene / transcript / ordered samples. No representative isoform; declared SAM units are not molecules.");
        ui.label("Paste the typed synthetic or independently reviewed command draft. Reviewing a draft performs no file I/O.");
        let settings_height = (ui.available_height() * 0.45).clamp(80.0, 320.0);
        egui::CollapsingHeader::new("Requests, explicit selections, view and export settings").id_salt("read_end_settings").default_open(false).show(ui,|ui| {
        egui::ScrollArea::vertical().id_salt("read_end_settings_scroll").max_height(settings_height).show(ui,|ui| {
        let changed = ui
            .add(
                egui::TextEdit::multiline(&mut state.draft)
                    .desired_rows(4)
                    .code_editor(),
            )
            .changed();
        if changed {
            state.reviewed = None;
            state.approved = false;
        }
        if ui.button("Review exact command draft").clicked() {
            match parse_draft(&state.draft) {
                Ok(command) => {
                    state.reviewed = Some(command);
                    state.approved = false;
                    state.error = None;
                }
                Err(error) => state.error = Some(error),
            }
        }
        if let Some(command) = &mut state.reviewed {
            let mut selection_changed=false;
            if let Command::Build { request } = command {
                ui.horizontal_wrapped(|ui| {
                    ui.label("Source");
                    selection_changed |= ui.text_edit_singleline(&mut request.source_id).changed();
                    ui.label("Gene");
                    selection_changed |= ui.text_edit_singleline(&mut request.gene_id).changed();
                    ui.label("Selected transcript");
                    selection_changed |= ui
                        .text_edit_singleline(&mut request.transcript_id)
                        .changed();
                });
                ui.label("Ordered sample bindings (edit the draft to add/remove/reorder or change file bindings):");
                for sample in &request.samples {
                    ui.label(format!(
                        "{}: {} | model {} | annotation {}",
                        sample.sample_id,
                        sample.analysis.report_path,
                        sample.model.sha256,
                        sample.normalized_annotation.sha256
                    ));
                }
            }
            if selection_changed {state.approved=false;state.draft=serde_json::to_string_pretty(command).unwrap_or_default();}
            if command.requires_confirmation() {
                ui.checkbox(&mut state.approved,"I approve this exact file/tool/write request; no biological acceptance is implied");
            }
            if ui
                .add_enabled(
                    state.active.is_none() && (!command.requires_confirmation() || state.approved),
                    egui::Button::new("Run reviewed command"),
                )
                .clicked()
            {
                action = Some(command.clone());
            }
        }
        ui.separator();
        ui.label("Completed read-end report binding (path, raw SHA-256, bytes):");
        if ui.text_edit_multiline(&mut state.report_json).changed() {
            state.approved = false;
            state.view_dirty=true;
        }
        let mut view_changed = false;
        ui.horizontal_wrapped(|ui| {
            ui.label("Sort");
            egui::ComboBox::from_id_salt("read_end_sort")
                .selected_text(format!("{:?}", state.view.sort))
                .show_ui(ui, |ui| {
                    for key in [
                        SortKey::Sample,
                        SortKey::End5,
                        SortKey::End3,
                        SortKey::Status5,
                        SortKey::Status3,
                        SortKey::Value5,
                        SortKey::Value3,
                        SortKey::Compatibility,
                        SortKey::ReadLength,
                        SortKey::Clip5,
                        SortKey::Clip3,
                        SortKey::UnitKey,
                    ] {
                        view_changed |= ui
                            .selectable_value(&mut state.view.sort, key, format!("{key:?}"))
                            .changed();
                    }
                });
            view_changed |= ui
                .checkbox(&mut state.view.descending, "Descending")
                .changed();
            ui.label("First lane (0-based)");
            view_changed |= ui
                .add(egui::DragValue::new(&mut state.view.offset).range(0..=100_000))
                .changed();
            ui.label("Page units");
            view_changed |= ui
                .add(egui::DragValue::new(&mut state.view.count).range(1..=2000))
                .changed();
            ui.label("Maximum hypothesis sublanes");
            view_changed |= ui
                .add(egui::DragValue::new(&mut state.view.max_sublanes).range(1..=20_000))
                .changed();
            ui.label("Export width (px)");
            view_changed |= ui
                .add(egui::DragValue::new(&mut state.view.width).range(900..=8000))
                .changed();
        });
        ui.horizontal_wrapped(|ui| {
            let c = &mut state.view.channels;
            for (label, flag) in [
                ("5' UTR", &mut c.utr5),
                ("3' UTR", &mut c.utr3),
                ("Read length", &mut c.read_length),
                ("Clips", &mut c.clips),
                ("A/T", &mut c.composition),
                ("Signal estimate", &mut c.tail),
            ] {
                view_changed |= ui.checkbox(flag, label).changed();
            }
        });
        let old_filter = state.status_filter;
        egui::ComboBox::from_id_salt("read_end_status_filter")
            .selected_text(
                ["All statuses", "3' observed", "3' censored", "3' missing"][state.status_filter],
            )
            .show_ui(ui, |ui| {
                for (index, label) in ["All statuses", "3' observed", "3' censored", "3' missing"]
                    .into_iter()
                    .enumerate()
                {
                    ui.selectable_value(&mut state.status_filter, index, label);
                }
            });
        if old_filter != state.status_filter {
            state.view.filters = match state.status_filter {
                1 => vec![Filter::Status3 {
                    value: Status::Observed,
                }],
                2 => vec![Filter::Status3 {
                    value: Status::RightCensored,
                }],
                3 => vec![Filter::Status3 {
                    value: Status::Missing,
                }],
                _ => vec![],
            };
            view_changed = true;
        }
        ui.label(format!("Draft filters: {:?}", state.view.filters));
        if ui.button("Clear all filters").clicked() {
            state.view.filters.clear();
            state.status_filter = 0;
            view_changed = true;
        }
        if view_changed {
            state.view_dirty = true;
            state.approved = false;
        }
        if state.view_dirty {
            ui.label("View draft changed; the retained page below is the previous applied view. Apply to refresh.");
        }
        if ui
            .add_enabled(
                state.active.is_none() && !state.report_json.is_empty(),
                egui::Button::new("Apply view / inspect saved evidence"),
            )
            .clicked()
        {
            match serde_json::from_str(&state.report_json) {
                Ok(report) => {
                    action = Some(Command::Inspect {
                        report,
                        view: state.view.clone(),
                    })
                }
                Err(error) => state.error = Some(error.to_string()),
            }
        }
        ui.horizontal_wrapped(|ui| {
            ui.label("Fresh export directory");
            if ui.text_edit_singleline(&mut state.export_dir).changed() {
                state.approved = false;
            }
            if ui
                .add_enabled(
                    !state.report_json.is_empty() && !state.export_dir.is_empty(),
                    egui::Button::new("Prepare JSON / TSV / SVG / PNG / PDF export for review"),
                )
                .clicked()
            {
                match serde_json::from_str(&state.report_json) {
                    Ok(report) => {
                        let command = Command::Export {
                            report,
                            view: state.view.clone(),
                            formats: vec![
                                Format::Json,
                                Format::Tsv,
                                Format::Svg,
                                Format::Png,
                                Format::Pdf,
                            ],
                            output_dir: state.export_dir.clone(),
                            max_output_bytes: 128 * 1024 * 1024,
                            cancel_path: None,
                        };
                        state.draft = serde_json::to_string_pretty(&command).unwrap_or_default();
                        state.reviewed = Some(command);
                        state.approved = false;
                    }
                    Err(error) => state.error = Some(error.to_string()),
                }
            }
        });
        });
        });
        if let Some(id) = state.active {
            if ui.button("Request cancellation").clicked() {
                self.agent_command_service.cancel(id);
            }
            ui.label(
                "A request alone does not confirm termination. Await the terminal job receipt.",
            );
            ui.ctx()
                .request_repaint_after(std::time::Duration::from_millis(100));
        }
        if state.view_dirty && state.page.is_some() {
            ui.label("Retained page belongs to its displayed report/view hashes; apply the current draft to refresh.");
        }
        ui.label(&state.status).on_hover_text(&state.request_text);
        if let Some(error) = &state.error {
            ui.colored_label(ui.visuals().error_fg_color, error);
        }
        if let (Some(page), Some(scene)) = (&state.page, &state.scene) {
            ui.strong(format!(
                "All samples: gene-assigned {} | selected-transcript-compatible {}",
                page.total_gene_assigned, page.total_selected_transcript_compatible
            ));
            egui::ScrollArea::vertical().id_salt("read_end_sample_denominators").max_height(72.0).show(ui,|ui| {
            for sample in &page.coverage.samples {
                ui.strong(format!("{}: gene-assigned {} | selected-transcript-compatible {} | declared {} | assessed {}",sample.sample_id,sample.gene_assigned_units,sample.selected_transcript_compatible_units,sample.declared_units,sample.assessed_units));
            }
            });
            ui.label(format!(
                "Shown {}..{} / {} matching / {} gene-assigned; view {}",
                page.shown_start,
                page.shown_end,
                page.total_matching,
                page.total_gene_assigned,
                page.view_sha256
            ));
            if let Some(scope) = scene.header.first() {
                ui.label(scope);
            }
            egui::ScrollArea::horizontal()
                .id_salt("read_end_shared_horizontal_frame")
                .show(ui, |ui| {
                    ui.set_min_width(scene.width as f32);
                    if let Ok(svg) = gentle_render::read_end_pileup::render_model(scene) {
                        let id = format!(
                            "bytes://read-end-header-{}-{}.svg",
                            page.report_sha256, page.view_sha256
                        );
                        ui.add(
                            egui::Image::from_bytes(id, svg.into_bytes()).fit_to_exact_size(
                                egui::vec2(scene.width as f32, state.model_height),
                            ),
                        );
                    }
                    egui::ScrollArea::vertical()
                        .id_salt("read_end_virtual_lanes")
                        .max_height(ui.available_height().max(1.0))
                        .show_rows(ui, state.max_lane_height, scene.rows.len(), |ui, range| {
                            ui.set_min_width(scene.width as f32);
                            for index in range {
                                if let Ok(svg) =
                                    gentle_render::read_end_pileup::render_lane(scene, index)
                                {
                                    let id = format!(
                                        "bytes://read-end-{}-{}-{index}.svg",
                                        page.report_sha256, page.view_sha256
                                    );
                                    ui.vertical(|ui| {
                                        ui.set_min_height(state.max_lane_height);
                                        ui.add(
                                            egui::Image::from_bytes(id, svg.into_bytes())
                                                .fit_to_exact_size(egui::vec2(
                                                    scene.width as f32,
                                                    state.lane_heights[index] as f32,
                                                )),
                                        )
                                        .on_hover_text(&scene.rows[index].details);
                                    });
                                }
                            }
                        });
                });
        }
        if let Some(command) = action {
            if let Err(error) = self.submit_read_end(command) {
                self.rna_panel_ui.read_ends.error = Some(error);
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn read_end_gui_starts_without_transcript_or_sample_selection_and_does_not_read_files() {
        let state = UiState::default();
        assert!(state.draft.is_empty() && state.reviewed.is_none() && state.report_json.is_empty());
        assert!(!state.approved);
        assert!(parse_draft("@not-a-file-read").is_err());
        assert!(parse_draft(&" ".repeat(1024 * 1024 + 1)).is_err());
    }
    #[test]
    fn read_end_gui_view_keeps_separate_channels_and_explicit_paging() {
        let state = UiState::default();
        assert_eq!(state.view.count, 100);
        assert_eq!(state.view.scale, "genomic");
        assert!(state.view.filters.is_empty());
        assert!(state.view.channels.utr5 && state.view.channels.utr3 && state.view.channels.clips);
    }
    #[test]
    fn read_end_gui_late_inspection_cannot_relabel_an_edited_binding_or_view() {
        let mut state = UiState::default();
        let binding = serde_json::json!({"path":"synthetic report","bytes":1,"sha256":"synthetic bound hash"});
        state.report_json = serde_json::to_string(&binding).unwrap();
        state.active_inspection_binding = Some(binding);
        let hash = pileup::hash(&state.view).unwrap();
        assert!(matches_draft(&state, &hash));
        state.view.channels.utr3 = false;
        assert!(!matches_draft(&state, &hash));
        state.view.channels.utr3 = true;
        state.report_json = "{}".into();
        assert!(!matches_draft(&state, &hash));
    }
}
