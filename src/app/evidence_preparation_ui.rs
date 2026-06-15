//! Guided TP73 evidence-viewer preparation assistant.

use super::*;

const TP73_EVIDENCE_WORKFLOW_JSON: &str =
    include_str!("../../docs/examples/workflows/tp73_genome_evidence_viewer_release_proof.json");

#[derive(Clone)]
pub(super) struct EvidencePreparationPanelState {
    pub(super) show_panel: bool,
    defaults: EvidencePreparationDefaults,
    seq_id: String,
    sequence_path: String,
    rmsk_index_path: String,
    repeat_report_path: String,
    repeat_max_features: String,
    array_manifest_path: String,
    array_contrasts: String,
    array_level: String,
    array_max_features: String,
    bed_path: String,
    bed_track_name: String,
    tfbs_motifs: String,
    tfbs_min_llr_quantile: String,
    tfbs_max_hits: String,
    sequence_svg_path: String,
    feature_expert_svg_path: String,
    tfbs_score_tracks_path: String,
    tfbs_score_tracks_svg_path: String,
    probe_region_output_dir: String,
    status: String,
    error: Option<String>,
}

impl Default for EvidencePreparationPanelState {
    fn default() -> Self {
        let defaults = tp73_evidence_preparation_defaults()
            .unwrap_or_else(|error| EvidencePreparationDefaults::fallback(&error));
        Self::from_defaults(defaults)
    }
}

impl EvidencePreparationPanelState {
    fn from_defaults(defaults: EvidencePreparationDefaults) -> Self {
        Self {
            show_panel: false,
            seq_id: defaults.seq_id.clone(),
            sequence_path: defaults.sequence_path.clone(),
            rmsk_index_path: defaults.rmsk_index_path.clone(),
            repeat_report_path: defaults.repeat_report_path.clone(),
            repeat_max_features: defaults.repeat_max_features.to_string(),
            array_manifest_path: defaults.array_manifest_path.clone(),
            array_contrasts: defaults.array_contrasts.join(","),
            array_level: defaults.array_level.clone(),
            array_max_features: defaults.array_max_features.to_string(),
            bed_path: defaults.bed_path.clone(),
            bed_track_name: defaults.bed_track_name.clone(),
            tfbs_motifs: defaults.tfbs_motifs.join(","),
            tfbs_min_llr_quantile: defaults.tfbs_min_llr_quantile.to_string(),
            tfbs_max_hits: defaults.tfbs_max_hits.to_string(),
            sequence_svg_path: defaults.sequence_svg_path.clone(),
            feature_expert_svg_path: defaults.feature_expert_svg_path.clone(),
            tfbs_score_tracks_path: defaults.tfbs_score_tracks_path.clone(),
            tfbs_score_tracks_svg_path: defaults.tfbs_score_tracks_svg_path.clone(),
            probe_region_output_dir: "analysis/probe_regions".to_string(),
            status: "Ready. Defaults are loaded from the TP73 evidence-viewer proof workflow."
                .to_string(),
            error: None,
            defaults,
        }
    }

    fn reload_defaults(&mut self) {
        *self = Self::default();
        self.show_panel = true;
    }
}

#[derive(Clone, Debug, PartialEq)]
pub(super) struct EvidencePreparationDefaults {
    pub(super) seq_id: String,
    pub(super) sequence_path: String,
    pub(super) rmsk_index_path: String,
    pub(super) repeat_report_path: String,
    pub(super) repeat_max_features: usize,
    pub(super) array_manifest_path: String,
    pub(super) array_contrasts: Vec<String>,
    pub(super) array_level: String,
    pub(super) array_max_features: usize,
    pub(super) bed_path: String,
    pub(super) bed_track_name: String,
    pub(super) tfbs_motifs: Vec<String>,
    pub(super) tfbs_min_llr_quantile: f64,
    pub(super) tfbs_max_hits: usize,
    pub(super) sequence_svg_path: String,
    pub(super) feature_expert_feature_id: usize,
    pub(super) feature_expert_svg_path: String,
    pub(super) tfbs_score_tracks_path: String,
    pub(super) tfbs_score_tracks_svg_path: String,
    parse_warning: Option<String>,
}

impl EvidencePreparationDefaults {
    fn fallback(error: &str) -> Self {
        Self {
            seq_id: "tp73_evidence_viewer".to_string(),
            sequence_path: "test_files/tp73.ncbi.gb".to_string(),
            rmsk_index_path:
                "test_files/fixtures/evidence_viewer/tp73_evidence_viewer.rmsk.hg38.interval-index.json"
                    .to_string(),
            repeat_report_path:
                "artifacts/tp73_evidence_viewer/tp73_evidence_viewer.repeat_materialization.json"
                    .to_string(),
            repeat_max_features: 10,
            array_manifest_path:
                "test_files/fixtures/evidence_viewer/clariomd.tp73_evidence_viewer.manifest.json"
                    .to_string(),
            array_contrasts: vec![
                "AdTAp73alpha-AdGFP".to_string(),
                "AdTAp73beta-AdGFP".to_string(),
            ],
            array_level: "probeset".to_string(),
            array_max_features: 20,
            bed_path: "test_files/fixtures/evidence_viewer/tp73_cutrun_demo.bed".to_string(),
            bed_track_name: "TP73 CUT&RUN proof BED".to_string(),
            tfbs_motifs: vec!["TP73".to_string(), "TP53".to_string(), "TP63".to_string()],
            tfbs_min_llr_quantile: 0.9995,
            tfbs_max_hits: 200,
            sequence_svg_path:
                "artifacts/tp73_evidence_viewer/tp73_evidence_viewer.linear.svg".to_string(),
            feature_expert_feature_id: 4,
            feature_expert_svg_path:
                "artifacts/tp73_evidence_viewer/tp73_evidence_viewer.splicing.expert.svg"
                    .to_string(),
            tfbs_score_tracks_path:
                "artifacts/tp73_evidence_viewer/tp73_evidence_viewer.tfbs_score_tracks.json"
                    .to_string(),
            tfbs_score_tracks_svg_path:
                "artifacts/tp73_evidence_viewer/tp73_evidence_viewer.tfbs_score_tracks.svg"
                    .to_string(),
            parse_warning: Some(format!("Could not parse proof workflow defaults: {error}")),
        }
    }
}

pub(super) fn tp73_evidence_preparation_defaults()
-> std::result::Result<EvidencePreparationDefaults, String> {
    let value: serde_json::Value = serde_json::from_str(TP73_EVIDENCE_WORKFLOW_JSON)
        .map_err(|error| format!("workflow JSON parse failed: {error}"))?;
    let ops = value
        .pointer("/workflow/ops")
        .and_then(serde_json::Value::as_array)
        .ok_or_else(|| "workflow.ops is missing or is not an array".to_string())?;

    let load = workflow_op_object(ops, "LoadFile")?;
    let repeat = workflow_op_object(ops, "MaterializeRepeatFeatures")?;
    let array = workflow_op_object(ops, "ProjectMicroarrayTrack")?;
    let bed = workflow_op_object(ops, "ImportGenomeBedTrack")?;
    let tfbs = workflow_op_object(ops, "AnnotateTfbs")?;
    let sequence_svg = workflow_op_object(ops, "RenderSequenceSvg")?;
    let feature_svg = workflow_op_object(ops, "RenderFeatureExpertSvg")?;
    let tfbs_score = workflow_op_object(ops, "SummarizeTfbsScoreTracks")?;
    let tfbs_svg = workflow_op_object(ops, "RenderTfbsScoreTracksSvg")?;

    Ok(EvidencePreparationDefaults {
        seq_id: workflow_string(load, "as_id")?,
        sequence_path: workflow_string(load, "path")?,
        rmsk_index_path: workflow_string(repeat, "rmsk_index_path")?,
        repeat_report_path: workflow_string(repeat, "path")?,
        repeat_max_features: workflow_usize(repeat, "max_features")?,
        array_manifest_path: workflow_string(array, "manifest_path")?,
        array_contrasts: workflow_string_array(array, "contrasts")?,
        array_level: workflow_string(array, "level")?,
        array_max_features: workflow_usize(array, "max_features")?,
        bed_path: workflow_string(bed, "path")?,
        bed_track_name: workflow_string(bed, "track_name")?,
        tfbs_motifs: workflow_string_array(tfbs, "motifs")?,
        tfbs_min_llr_quantile: workflow_f64(tfbs, "min_llr_quantile")?,
        tfbs_max_hits: workflow_usize(tfbs, "max_hits")?,
        sequence_svg_path: workflow_string(sequence_svg, "path")?,
        feature_expert_feature_id: feature_expert_feature_id(feature_svg)?,
        feature_expert_svg_path: workflow_string(feature_svg, "path")?,
        tfbs_score_tracks_path: workflow_string(tfbs_score, "path")?,
        tfbs_score_tracks_svg_path: workflow_string(tfbs_svg, "path")?,
        parse_warning: None,
    })
}

impl GENtleApp {
    pub(super) fn open_evidence_preparation_dialog(&mut self) {
        let was_open = self.evidence_preparation_panel.show_panel;
        self.evidence_preparation_panel.show_panel = true;
        if let Some((seq_id, _)) = self.active_dna_window_context()
            && self.evidence_preparation_panel.seq_id.trim().is_empty()
        {
            self.evidence_preparation_panel.seq_id = seq_id;
        }
        self.mark_window_open_or_focus(Self::evidence_preparation_viewport_id(), was_open);
    }

    pub(super) fn render_evidence_preparation_dialog(&mut self, ctx: &egui::Context) {
        if !self.evidence_preparation_panel.show_panel {
            return;
        }
        let mut open = self.evidence_preparation_panel.show_panel;
        let viewport_id = Self::evidence_preparation_viewport_id();
        let spec = self.hosted_window_spec_for_viewport(
            "Evidence Preparation",
            Self::hosted_evidence_preparation_window_id(),
            viewport_id,
            Vec2::new(1120.0, 760.0),
            Vec2::new(760.0, 520.0),
        );
        crate::egui_compat::show_hosted_window(ctx, &spec, &mut open, |ui| {
            self.render_evidence_preparation_contents(ui);
        });
        self.clear_viewport_foreground_request_after_render(viewport_id);
        self.finalize_viewport_open_probe(viewport_id, "Evidence Preparation");
        self.evidence_preparation_panel.show_panel =
            open && self.evidence_preparation_panel.show_panel;
    }

    fn render_evidence_preparation_contents(&mut self, ui: &mut Ui) {
        if self.render_specialist_window_nav_with_close(
            ui,
            Some(("Close", "Close Evidence Preparation")),
        ) {
            self.evidence_preparation_panel.show_panel = false;
        }
        ui.heading("TP73 Evidence Preparation");
        ui.small("Prepare the local proof material used by the TP73 genome-anchored evidence viewer. Runnable buttons call shared GENtle operations; external R/CEL steps stay copy-command handoffs.");
        if let Some(warning) = self
            .evidence_preparation_panel
            .defaults
            .parse_warning
            .as_deref()
        {
            ui.colored_label(egui::Color32::from_rgb(150, 96, 0), warning);
        }
        if let Some(error) = &self.evidence_preparation_panel.error {
            ui.colored_label(egui::Color32::RED, error);
        }
        if !self.evidence_preparation_panel.status.is_empty() {
            ui.small(self.evidence_preparation_panel.status.as_str());
        }
        ui.separator();
        self.render_evidence_preparation_inputs(ui);
        ui.separator();
        egui::ScrollArea::vertical()
            .id_salt("evidence_preparation_steps_scroll")
            .show(ui, |ui| {
                self.render_evidence_preparation_runnable_cards(ui);
                ui.separator();
                self.render_evidence_preparation_manual_cards(ui);
            });
    }

    fn render_evidence_preparation_inputs(&mut self, ui: &mut Ui) {
        ui.horizontal(|ui| {
            ui.label("sequence id");
            ui.text_edit_singleline(&mut self.evidence_preparation_panel.seq_id);
            if ui.button("Use active sequence").clicked()
                && let Some((seq_id, _)) = self.active_dna_window_context()
            {
                self.evidence_preparation_panel.seq_id = seq_id;
            }
            if ui.button("Reload proof defaults").clicked() {
                self.evidence_preparation_panel.reload_defaults();
            }
        });
        ui.horizontal(|ui| {
            ui.label("TP73 GenBank");
            ui.text_edit_singleline(&mut self.evidence_preparation_panel.sequence_path);
        });
        ui.horizontal(|ui| {
            ui.label("output SVG");
            ui.text_edit_singleline(&mut self.evidence_preparation_panel.sequence_svg_path);
        });
    }

    fn render_evidence_preparation_runnable_cards(&mut self, ui: &mut Ui) {
        self.evidence_card(ui, "1. Anchored sequence", |app, ui| {
            ui.small("Load the committed GRCh38.p14 TP73 GenBank locus with the proof sequence id.");
            let command = format!(
                "cargo run --bin gentle_cli -- op '{{\"LoadFile\":{{\"path\":\"{}\",\"as_id\":\"{}\"}}}}'",
                app.evidence_preparation_panel.sequence_path,
                app.evidence_preparation_panel.seq_id
            );
            app.render_copyable_command(ui, &command);
            ui.horizontal(|ui| {
                if ui.button("Load TP73 proof sequence").clicked() {
                    app.run_evidence_operation(
                        "Load TP73 proof sequence",
                        app.evidence_load_sequence_operation(),
                    );
                }
                if ui.button("Open DNA viewer").clicked() {
                    let seq_id = app.evidence_preparation_panel.seq_id.trim().to_string();
                    app.open_sequence_window(&seq_id);
                }
            });
        });

        self.evidence_card(ui, "2. Repeats", |app, ui| {
            ui.label("rmsk interval index");
            ui.text_edit_singleline(&mut app.evidence_preparation_panel.rmsk_index_path);
            ui.horizontal(|ui| {
                ui.label("max features");
                ui.add(
                    egui::TextEdit::singleline(
                        &mut app.evidence_preparation_panel.repeat_max_features,
                    )
                    .desired_width(72.0),
                );
                ui.label("report path");
                ui.text_edit_singleline(&mut app.evidence_preparation_panel.repeat_report_path);
            });
            let command = app.evidence_repeat_command();
            app.render_copyable_command(ui, &command);
            if ui.button("Materialize repeats").clicked() {
                match app.evidence_repeat_operation() {
                    Ok(op) => app.run_evidence_operation("Materialize repeats", op),
                    Err(error) => app.evidence_preparation_panel.error = Some(error),
                }
            }
        });

        self.evidence_card(ui, "3. Clariom D array track", |app, ui| {
            ui.label("manifest");
            ui.text_edit_singleline(&mut app.evidence_preparation_panel.array_manifest_path);
            ui.horizontal(|ui| {
                ui.label("contrasts");
                ui.text_edit_singleline(&mut app.evidence_preparation_panel.array_contrasts);
            });
            ui.horizontal(|ui| {
                ui.label("level");
                ui.text_edit_singleline(&mut app.evidence_preparation_panel.array_level);
                ui.label("max features");
                ui.add(
                    egui::TextEdit::singleline(
                        &mut app.evidence_preparation_panel.array_max_features,
                    )
                    .desired_width(72.0),
                );
            });
            app.render_copyable_command(ui, &app.evidence_array_inspect_command());
            app.render_copyable_command(ui, &app.evidence_array_project_command());
            ui.horizontal(|ui| {
                if ui.button("Inspect manifest").clicked() {
                    app.run_evidence_shell(
                        "Inspect Clariom D manifest",
                        app.evidence_array_inspect_shell_command(),
                    );
                }
                if ui.button("Project array track").clicked() {
                    match app.evidence_array_project_operation() {
                        Ok(op) => app.run_evidence_operation("Project Clariom D array track", op),
                        Err(error) => app.evidence_preparation_panel.error = Some(error),
                    }
                }
            });
        });

        self.evidence_card(ui, "4. CUT&RUN-style BED", |app, ui| {
            ui.label("BED path");
            ui.text_edit_singleline(&mut app.evidence_preparation_panel.bed_path);
            ui.label("track name");
            ui.text_edit_singleline(&mut app.evidence_preparation_panel.bed_track_name);
            app.render_copyable_command(ui, &app.evidence_bed_command());
            if ui.button("Import BED track").clicked() {
                app.run_evidence_operation(
                    "Import CUT&RUN-style BED",
                    app.evidence_bed_operation(),
                );
            }
        });

        self.evidence_card(ui, "5. TFBS and visibility", |app, ui| {
            ui.horizontal(|ui| {
                ui.label("motifs");
                ui.text_edit_singleline(&mut app.evidence_preparation_panel.tfbs_motifs);
            });
            ui.horizontal(|ui| {
                ui.label("min LLR quantile");
                ui.add(
                    egui::TextEdit::singleline(
                        &mut app.evidence_preparation_panel.tfbs_min_llr_quantile,
                    )
                    .desired_width(88.0),
                );
                ui.label("max hits");
                ui.add(
                    egui::TextEdit::singleline(&mut app.evidence_preparation_panel.tfbs_max_hits)
                        .desired_width(72.0),
                );
            });
            app.render_copyable_command(ui, &app.evidence_tfbs_command());
            ui.horizontal(|ui| {
                if ui.button("Annotate TFBS").clicked() {
                    match app.evidence_tfbs_operation() {
                        Ok(op) => app.run_evidence_operation("Annotate TFBS", op),
                        Err(error) => app.evidence_preparation_panel.error = Some(error),
                    }
                }
                if ui.button("Show evidence layers").clicked() {
                    app.run_evidence_visibility_operations();
                }
            });
        });

        self.evidence_card(ui, "6. Proof exports", |app, ui| {
            ui.horizontal(|ui| {
                ui.label("splicing SVG");
                ui.text_edit_singleline(
                    &mut app.evidence_preparation_panel.feature_expert_svg_path,
                );
            });
            ui.horizontal(|ui| {
                ui.label("TFBS JSON");
                ui.text_edit_singleline(&mut app.evidence_preparation_panel.tfbs_score_tracks_path);
            });
            ui.horizontal(|ui| {
                ui.label("TFBS SVG");
                ui.text_edit_singleline(
                    &mut app.evidence_preparation_panel.tfbs_score_tracks_svg_path,
                );
            });
            app.render_copyable_command(ui, &app.evidence_render_sequence_command());
            if ui.button("Render proof artifacts").clicked() {
                match app.evidence_render_operations() {
                    Ok(ops) => app.run_evidence_operations("Render proof artifacts", ops),
                    Err(error) => app.evidence_preparation_panel.error = Some(error),
                }
            }
        });
    }

    fn render_evidence_preparation_manual_cards(&mut self, ui: &mut Ui) {
        self.evidence_card(ui, "External CEL/probeset preparation", |app, ui| {
            ui.small("GENtle does not auto-download login-walled vendor assets or run CEL summarization here. Use these handoff commands, then inspect the output before projecting selected rows.");
            app.render_copyable_command(ui, &app.evidence_probe_regions_dry_run_command());
            app.render_copyable_command(ui, &app.evidence_probe_regions_r_command());
            ui.horizontal(|ui| {
                ui.label("output dir");
                ui.text_edit_singleline(&mut app.evidence_preparation_panel.probe_region_output_dir);
            });
            app.render_copyable_command(ui, &app.evidence_probe_regions_inspect_command());
            ui.small("Thermo Fisher Clariom D na36 hg38 support ZIPs are manually staged under data/resources/affymetrix/clariom_d_human_na36_hg38/. External Services remains the broader vendor handoff surface; this card stays limited to array-analysis preparation commands.");
            if ui.button("Run probe-region preflight").clicked() {
                app.run_evidence_shell(
                    "Probe-region preflight",
                    app.evidence_probe_regions_shell_command(),
                );
            }
            if ui.button("Inspect helper output").clicked() {
                app.run_evidence_shell(
                    "Inspect probe-region helper output",
                    app.evidence_probe_regions_inspect_shell_command(),
                );
            }
        });
    }

    fn evidence_card(
        &mut self,
        ui: &mut Ui,
        title: &str,
        add_contents: impl FnOnce(&mut GENtleApp, &mut Ui),
    ) {
        egui::Frame::group(ui.style()).show(ui, |ui| {
            ui.strong(title);
            add_contents(self, ui);
        });
        ui.add_space(6.0);
    }

    fn render_copyable_command(&mut self, ui: &mut Ui, command: &str) {
        let mut display = command.to_string();
        ui.horizontal(|ui| {
            ui.add(
                egui::TextEdit::singleline(&mut display)
                    .desired_width(f32::INFINITY)
                    .code_editor(),
            );
            if ui.button("Copy").clicked() {
                ui.ctx().copy_text(command.to_string());
                self.evidence_preparation_panel.status = "Copied command.".to_string();
            }
        });
    }

    fn run_evidence_operation(&mut self, label: &str, op: Operation) {
        self.run_evidence_operations(label, vec![op]);
    }

    fn run_evidence_operations(&mut self, label: &str, ops: Vec<Operation>) {
        let mut created = vec![];
        let mut changed = vec![];
        let mut warnings = vec![];
        let mut messages = vec![];
        let result = self
            .engine
            .write()
            .map_err(|_| "Could not lock GENtle engine for evidence preparation.".to_string())
            .and_then(|mut engine| {
                for op in ops {
                    let result = engine.apply(op).map_err(|error| error.message)?;
                    created.extend(result.created_seq_ids);
                    changed.extend(result.changed_seq_ids);
                    warnings.extend(result.warnings);
                    messages.extend(result.messages);
                }
                Ok(())
            });
        match result {
            Ok(()) => {
                self.evidence_preparation_panel.error = None;
                let mut refresh_ids = created.clone();
                refresh_ids.extend(changed.clone());
                let refreshed = self.refresh_sequence_windows_for_seq_ids(&refresh_ids);
                for seq_id in &created {
                    self.open_sequence_window(seq_id);
                }
                self.evidence_preparation_panel.status = format!(
                    "{}\nrefreshed windows: {}",
                    Self::format_op_result_status(label, &created, &warnings, &messages),
                    refreshed
                );
                self.app_status = format!("{label}: complete");
            }
            Err(error) => {
                self.evidence_preparation_panel.error = Some(format!("{label} failed: {error}"));
                self.evidence_preparation_panel.status.clear();
            }
        }
    }

    fn run_evidence_shell(&mut self, label: &str, command: ShellCommand) {
        let result = self
            .engine
            .write()
            .map_err(|_| "Could not lock GENtle engine for evidence shell command.".to_string())
            .and_then(|mut engine| {
                execute_shell_command_with_options(
                    &mut engine,
                    &command,
                    &ShellExecutionOptions::default(),
                )
            });
        match result {
            Ok(run) => {
                self.evidence_preparation_panel.error = None;
                self.evidence_preparation_panel.status = format!(
                    "{label}: complete (state_changed={})\n{}",
                    run.state_changed,
                    serde_json::to_string_pretty(&run.output)
                        .unwrap_or_else(|_| "<output unavailable>".to_string())
                );
            }
            Err(error) => {
                self.evidence_preparation_panel.error = Some(format!("{label} failed: {error}"));
                self.evidence_preparation_panel.status.clear();
            }
        }
    }

    pub(super) fn evidence_load_sequence_operation(&self) -> Operation {
        Operation::LoadFile {
            path: self
                .evidence_preparation_panel
                .sequence_path
                .trim()
                .to_string(),
            as_id: Some(self.evidence_preparation_panel.seq_id.trim().to_string()),
        }
    }

    pub(super) fn evidence_repeat_operation(&self) -> std::result::Result<Operation, String> {
        Ok(Operation::MaterializeRepeatFeatures {
            seq_id: self.evidence_seq_id(),
            rmsk_index_path: self
                .evidence_preparation_panel
                .rmsk_index_path
                .trim()
                .to_string(),
            max_features: Some(self.parse_evidence_usize("repeat max features")?),
            clear_existing: Some(true),
            path: evidence_nonempty_text(&self.evidence_preparation_panel.repeat_report_path),
        })
    }

    pub(super) fn evidence_array_project_operation(
        &self,
    ) -> std::result::Result<Operation, String> {
        Ok(Operation::ProjectMicroarrayTrack {
            seq_id: self.evidence_seq_id(),
            manifest_path: self
                .evidence_preparation_panel
                .array_manifest_path
                .trim()
                .to_string(),
            contrasts: split_evidence_csv(&self.evidence_preparation_panel.array_contrasts),
            level: evidence_nonempty_text(&self.evidence_preparation_panel.array_level),
            min_abs_logfc: Some(0.0),
            max_adj_p: Some(1.0),
            max_features: Some(self.parse_evidence_usize("array max features")?),
            clear_existing: Some(true),
        })
    }

    pub(super) fn evidence_bed_operation(&self) -> Operation {
        Operation::ImportGenomeBedTrack {
            seq_id: self.evidence_seq_id(),
            path: self.evidence_preparation_panel.bed_path.trim().to_string(),
            track_name: evidence_nonempty_text(&self.evidence_preparation_panel.bed_track_name),
            min_score: None,
            max_score: None,
            clear_existing: Some(true),
        }
    }

    pub(super) fn evidence_tfbs_operation(&self) -> std::result::Result<Operation, String> {
        Ok(Operation::AnnotateTfbs {
            seq_id: self.evidence_seq_id(),
            motifs: split_evidence_csv(&self.evidence_preparation_panel.tfbs_motifs),
            min_llr_bits: None,
            min_llr_quantile: Some(self.parse_evidence_f64("TFBS min LLR quantile")?),
            per_tf_thresholds: vec![],
            clear_existing: Some(true),
            max_hits: Some(self.parse_evidence_usize("TFBS max hits")?),
        })
    }

    pub(super) fn evidence_render_operations(&self) -> std::result::Result<Vec<Operation>, String> {
        let target = SequenceScanTarget::SeqId {
            seq_id: self.evidence_seq_id(),
            span_start_0based: Some(0),
            span_end_0based_exclusive: Some(1200),
        };
        Ok(vec![
            Operation::SetLinearViewport {
                start_bp: 0,
                span_bp: 1200,
            },
            Operation::RenderSequenceSvg {
                seq_id: self.evidence_seq_id(),
                mode: RenderSvgMode::Linear,
                path: self
                    .evidence_preparation_panel
                    .sequence_svg_path
                    .trim()
                    .to_string(),
            },
            Operation::RenderFeatureExpertSvg {
                seq_id: self.evidence_seq_id(),
                target: FeatureExpertTarget::SplicingFeature {
                    feature_id: self
                        .evidence_preparation_panel
                        .defaults
                        .feature_expert_feature_id,
                    scope: SplicingScopePreset::AllOverlappingAnyStrand,
                },
                path: self
                    .evidence_preparation_panel
                    .feature_expert_svg_path
                    .trim()
                    .to_string(),
            },
            Operation::SummarizeTfbsScoreTracks {
                target: target.clone(),
                motifs: split_evidence_csv(&self.evidence_preparation_panel.tfbs_motifs),
                score_kind: TfbsScoreTrackValueKind::LlrBackgroundTailLog10,
                clip_negative: false,
                path: evidence_nonempty_text(
                    &self.evidence_preparation_panel.tfbs_score_tracks_path,
                ),
            },
            Operation::RenderTfbsScoreTracksSvg {
                target,
                motifs: split_evidence_csv(&self.evidence_preparation_panel.tfbs_motifs),
                score_kind: TfbsScoreTrackValueKind::LlrBackgroundTailLog10,
                clip_negative: false,
                path: self
                    .evidence_preparation_panel
                    .tfbs_score_tracks_svg_path
                    .trim()
                    .to_string(),
            },
        ])
    }

    fn run_evidence_visibility_operations(&mut self) {
        self.run_evidence_operations(
            "Show evidence layers",
            vec![
                Operation::SetDisplayVisibility {
                    target: DisplayTarget::RepeatFeatures,
                    visible: true,
                },
                Operation::SetDisplayVisibility {
                    target: DisplayTarget::ArrayFeatures,
                    visible: true,
                },
                Operation::SetDisplayVisibility {
                    target: DisplayTarget::Tfbs,
                    visible: true,
                },
            ],
        );
    }

    pub(super) fn evidence_array_inspect_shell_command(&self) -> ShellCommand {
        ShellCommand::ArraysInspectMicroarrayTrack {
            manifest_path: self
                .evidence_preparation_panel
                .array_manifest_path
                .trim()
                .to_string(),
        }
    }

    pub(super) fn evidence_probe_regions_shell_command(&self) -> ShellCommand {
        ShellCommand::ArraysProbeRegions {
            cel_paths: vec![],
            dataset: Some("E-MTAB-14704".to_string()),
            metadata_path: None,
            genes: vec!["TP73".to_string()],
            loci: vec![],
            transcript_cluster_ids: vec![],
            probeset_ids: vec![],
            platform: Some("Clariom_D_Human".to_string()),
            annotation_library_path: None,
            condition_column: None,
            sample_column: None,
            block_column: None,
            paired_by_replicate_suffix: false,
            plot: false,
            normalization: "rma".to_string(),
            output_dir: Some(
                self.evidence_preparation_panel
                    .probe_region_output_dir
                    .clone(),
            ),
            cache_dir: None,
            dry_run: true,
        }
    }

    pub(super) fn evidence_probe_regions_inspect_shell_command(&self) -> ShellCommand {
        ShellCommand::ArraysInspectProbeRegionOutput {
            output_dir: self
                .evidence_preparation_panel
                .probe_region_output_dir
                .trim()
                .to_string(),
        }
    }

    fn evidence_seq_id(&self) -> String {
        self.evidence_preparation_panel.seq_id.trim().to_string()
    }

    fn parse_evidence_usize(&self, field: &str) -> std::result::Result<usize, String> {
        let raw = match field {
            "repeat max features" => &self.evidence_preparation_panel.repeat_max_features,
            "array max features" => &self.evidence_preparation_panel.array_max_features,
            "TFBS max hits" => &self.evidence_preparation_panel.tfbs_max_hits,
            _ => return Err(format!("Unknown numeric field '{field}'")),
        };
        raw.trim()
            .parse::<usize>()
            .map_err(|error| format!("Invalid {field} '{}': {error}", raw.trim()))
    }

    fn parse_evidence_f64(&self, field: &str) -> std::result::Result<f64, String> {
        let raw = match field {
            "TFBS min LLR quantile" => &self.evidence_preparation_panel.tfbs_min_llr_quantile,
            _ => return Err(format!("Unknown numeric field '{field}'")),
        };
        raw.trim()
            .parse::<f64>()
            .map_err(|error| format!("Invalid {field} '{}': {error}", raw.trim()))
    }

    pub(super) fn evidence_array_inspect_command(&self) -> String {
        format!(
            "arrays inspect-microarray-track {}",
            self.evidence_preparation_panel.array_manifest_path.trim()
        )
    }

    pub(super) fn evidence_array_project_command(&self) -> String {
        format!(
            "arrays project-microarray-track {} {} --contrasts {} --level {} --max-features {} --clear-existing",
            self.evidence_seq_id(),
            self.evidence_preparation_panel.array_manifest_path.trim(),
            self.evidence_preparation_panel.array_contrasts.trim(),
            self.evidence_preparation_panel.array_level.trim(),
            self.evidence_preparation_panel.array_max_features.trim()
        )
    }

    pub(super) fn evidence_probe_regions_dry_run_command(&self) -> String {
        "arrays probe-regions --dataset E-MTAB-14704 --gene TP73 --platform Clariom_D_Human --output analysis/probe_regions --dry-run".to_string()
    }

    pub(super) fn evidence_probe_regions_r_command(&self) -> String {
        "Rscript scripts/probe_regions_oligo.R --dataset E-MTAB-14704 --gene TP73 --platform Clariom_D_Human --output analysis/probe_regions".to_string()
    }

    pub(super) fn evidence_probe_regions_inspect_command(&self) -> String {
        format!(
            "arrays inspect-probe-region-output {}",
            self.evidence_preparation_panel
                .probe_region_output_dir
                .trim()
        )
    }

    fn evidence_repeat_command(&self) -> String {
        format!(
            "op '{{\"MaterializeRepeatFeatures\":{{\"seq_id\":\"{}\",\"rmsk_index_path\":\"{}\",\"max_features\":{},\"clear_existing\":true,\"path\":\"{}\"}}}}'",
            self.evidence_seq_id(),
            self.evidence_preparation_panel.rmsk_index_path.trim(),
            self.evidence_preparation_panel.repeat_max_features.trim(),
            self.evidence_preparation_panel.repeat_report_path.trim()
        )
    }

    fn evidence_bed_command(&self) -> String {
        format!(
            "tracks import-bed {} {} --name \"{}\" --clear-existing",
            self.evidence_seq_id(),
            self.evidence_preparation_panel.bed_path.trim(),
            self.evidence_preparation_panel.bed_track_name.trim()
        )
    }

    fn evidence_tfbs_command(&self) -> String {
        format!(
            "op '{{\"AnnotateTfbs\":{{\"seq_id\":\"{}\",\"motifs\":{},\"min_llr_quantile\":{},\"clear_existing\":true,\"max_hits\":{}}}}}'",
            self.evidence_seq_id(),
            serde_json::to_string(&split_evidence_csv(
                &self.evidence_preparation_panel.tfbs_motifs
            ))
            .unwrap_or_else(|_| "[]".to_string()),
            self.evidence_preparation_panel.tfbs_min_llr_quantile.trim(),
            self.evidence_preparation_panel.tfbs_max_hits.trim()
        )
    }

    fn evidence_render_sequence_command(&self) -> String {
        format!(
            "render-svg {} linear {}",
            self.evidence_seq_id(),
            self.evidence_preparation_panel.sequence_svg_path.trim()
        )
    }
}

fn workflow_op_object<'a>(
    ops: &'a [serde_json::Value],
    op_name: &str,
) -> std::result::Result<&'a serde_json::Map<String, serde_json::Value>, String> {
    ops.iter()
        .find_map(|op| op.get(op_name).and_then(serde_json::Value::as_object))
        .ok_or_else(|| format!("workflow op '{op_name}' not found"))
}

fn workflow_string(
    object: &serde_json::Map<String, serde_json::Value>,
    key: &str,
) -> std::result::Result<String, String> {
    object
        .get(key)
        .and_then(serde_json::Value::as_str)
        .map(ToString::to_string)
        .ok_or_else(|| format!("workflow field '{key}' is missing or not a string"))
}

fn workflow_string_array(
    object: &serde_json::Map<String, serde_json::Value>,
    key: &str,
) -> std::result::Result<Vec<String>, String> {
    object
        .get(key)
        .and_then(serde_json::Value::as_array)
        .ok_or_else(|| format!("workflow field '{key}' is missing or not an array"))?
        .iter()
        .map(|value| {
            value
                .as_str()
                .map(ToString::to_string)
                .ok_or_else(|| format!("workflow field '{key}' contains a non-string value"))
        })
        .collect()
}

fn workflow_usize(
    object: &serde_json::Map<String, serde_json::Value>,
    key: &str,
) -> std::result::Result<usize, String> {
    let value = object
        .get(key)
        .and_then(serde_json::Value::as_u64)
        .ok_or_else(|| format!("workflow field '{key}' is missing or not an unsigned integer"))?;
    usize::try_from(value).map_err(|_| format!("workflow field '{key}' is too large"))
}

fn workflow_f64(
    object: &serde_json::Map<String, serde_json::Value>,
    key: &str,
) -> std::result::Result<f64, String> {
    object
        .get(key)
        .and_then(serde_json::Value::as_f64)
        .ok_or_else(|| format!("workflow field '{key}' is missing or not a number"))
}

fn feature_expert_feature_id(
    object: &serde_json::Map<String, serde_json::Value>,
) -> std::result::Result<usize, String> {
    let value = object
        .get("target")
        .and_then(|target| target.get("splicing_feature"))
        .and_then(|target| target.get("feature_id"))
        .and_then(serde_json::Value::as_u64)
        .ok_or_else(|| {
            "RenderFeatureExpertSvg target.splicing_feature.feature_id missing".to_string()
        })?;
    usize::try_from(value).map_err(|_| "splicing feature id is too large".to_string())
}

fn split_evidence_csv(raw: &str) -> Vec<String> {
    raw.split(',')
        .map(str::trim)
        .filter(|value| !value.is_empty())
        .map(ToString::to_string)
        .collect()
}

fn evidence_nonempty_text(value: &str) -> Option<String> {
    let trimmed = value.trim();
    (!trimmed.is_empty()).then(|| trimmed.to_string())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn evidence_preparation_defaults_parse_proof_workflow_json() {
        let defaults = tp73_evidence_preparation_defaults().expect("parse proof workflow");
        assert_eq!(defaults.seq_id, "tp73_evidence_viewer");
        assert_eq!(defaults.sequence_path, "test_files/tp73.ncbi.gb");
        assert_eq!(
            defaults.array_contrasts,
            vec!["AdTAp73alpha-AdGFP", "AdTAp73beta-AdGFP"]
        );
        assert_eq!(defaults.feature_expert_feature_id, 4);
        assert_eq!(defaults.tfbs_motifs, vec!["TP73", "TP53", "TP63"]);
    }

    #[test]
    fn evidence_preparation_copyable_commands_parse_shared_shell() {
        let app = GENtleApp::default();
        for command in [
            app.evidence_array_inspect_command(),
            app.evidence_array_project_command(),
            app.evidence_probe_regions_dry_run_command(),
            app.evidence_probe_regions_inspect_command(),
        ] {
            parse_shell_line(&command).unwrap_or_else(|error| {
                panic!("copyable command should parse: {command}\n{error}")
            });
        }
    }

    #[test]
    fn evidence_preparation_array_projection_matches_direct_operation_route() {
        let mut app = GENtleApp::default();
        let load = app.evidence_load_sequence_operation();
        app.run_evidence_operation("load", load);
        let project = app
            .evidence_array_project_operation()
            .expect("assistant project operation");
        app.run_evidence_operation("project", project.clone());

        let assistant_count = app
            .engine
            .read()
            .expect("engine")
            .state()
            .sequences
            .get("tp73_evidence_viewer")
            .expect("assistant sequence")
            .features()
            .iter()
            .filter(|feature| {
                feature
                    .qualifier_values("gentle_track_source")
                    .any(|value| value == "Array")
            })
            .count();

        let mut direct = GentleEngine::default();
        let defaults = tp73_evidence_preparation_defaults().expect("defaults");
        direct
            .apply(Operation::LoadFile {
                path: defaults.sequence_path,
                as_id: Some(defaults.seq_id.clone()),
            })
            .expect("direct load");
        direct.apply(project).expect("direct project");
        let direct_count = direct
            .state()
            .sequences
            .get("tp73_evidence_viewer")
            .expect("direct sequence")
            .features()
            .iter()
            .filter(|feature| {
                feature
                    .qualifier_values("gentle_track_source")
                    .any(|value| value == "Array")
            })
            .count();

        assert_eq!(assistant_count, direct_count);
        assert_eq!(assistant_count, 4);
    }
}
