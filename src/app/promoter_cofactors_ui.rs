//! Thin, background-loaded collaborator-package browser. No GUI-side SQL.

use super::*;
use gentle_protocol::promoter_cofactors::*;
use serde_json::Value;

pub(super) struct CofactorBrowser {
    pub open: bool,
    request: PromoterCofactorRequest,
    motif_filter: String,
    gene_filter: String,
    species_filter: String,
    duckdb_executable: String,
    chromosome: String,
    start: String,
    end: String,
    q_filter: bool,
    q_value: f64,
    report: Option<PromoterCofactorReport>,
    report_json: String,
    task: Option<mpsc::Receiver<Result<OpResult, EngineError>>>,
    status: String,
    submitted_form: Option<Value>,
    projection_seq_id: String,
    saved_feature_region: Option<gentle_protocol::GenomicRegionOfInterest>,
    saved_feature_regions: Vec<gentle_protocol::GenomicRegionOfInterest>,
    feature_seq_id: String,
    feature_preview: Option<gentle_protocol::GenomicRegionFeatureReport>,
    feature_confirmed: bool,
}

impl Default for CofactorBrowser {
    fn default() -> Self {
        Self {
            open: false,
            request: Default::default(),
            motif_filter: String::new(),
            gene_filter: String::new(),
            species_filter: String::new(),
            duckdb_executable: String::new(),
            chromosome: String::new(),
            start: String::new(),
            end: String::new(),
            q_filter: false,
            q_value: 0.05,
            report: None,
            report_json: String::new(),
            task: None,
            status: String::new(),
            submitted_form: None,
            projection_seq_id: String::new(),
            saved_feature_region: None,
            saved_feature_regions: vec![],
            feature_seq_id: String::new(),
            feature_preview: None,
            feature_confirmed: false,
        }
    }
}

impl CofactorBrowser {
    fn form_snapshot(&self) -> Value {
        serde_json::json!({"request":self.request,"motif":self.motif_filter,
            "gene":self.gene_filter,"species":self.species_filter,"duckdb":self.duckdb_executable,
            "chromosome":self.chromosome,"start":self.start,"end":self.end,
            "q_filter":self.q_filter,"q_value":self.q_value})
    }

    fn report_is_stale(&self) -> bool {
        self.report.is_some() && self.submitted_form.as_ref() != Some(&self.form_snapshot())
    }

    fn region_request(
        &self,
        target: CofactorRegionTarget,
    ) -> Option<gentle_protocol::GenomicRegionCaptureRequest> {
        Some(gentle_protocol::GenomicRegionCaptureRequest {
            set_id: "promoter_cofactors".into(),
            set_label: Some("Selected promoter-cofactor evidence".into()),
            label: Some(format!("{target:?}")),
            source: gentle_protocol::GenomicRegionCaptureSource::PromoterCofactor {
                report: Box::new(self.report.as_ref()?.clone()),
                target,
                seq_id: (!self.projection_seq_id.trim().is_empty())
                    .then(|| self.projection_seq_id.trim().into()),
            },
            ..Default::default()
        })
    }
}

fn region_buttons(
    ui: &mut egui::Ui,
    target: CofactorRegionTarget,
    action: &mut Option<(CofactorRegionTarget, bool)>,
) {
    if ui.button("Save evidence region").on_hover_text("Save the displayed interval and its original report in this project; no DNA or coverage import").clicked() {
        *action = Some((target.clone(), false));
    }
    if ui.button("Copy region request").clicked() {
        *action = Some((target, true));
    }
}

fn number(value: Option<f64>) -> String {
    value
        .map(|v| format!("{v:.4}"))
        .unwrap_or_else(|| "unavailable".into())
}

impl GENtleApp {
    fn reload_cofactor_feature_regions(&mut self) -> Result<(), String> {
        let (value, _) = self.execute_shared_shell_command_json(&ShellCommand::GenomicRegions {
            operation: Operation::ListGenomicRegions {
                request: gentle_protocol::GenomicRegionListRequest {
                    set_id: Some("promoter_cofactors".into()),
                },
            },
        })?;
        let report: gentle_protocol::GenomicRegionOperationReport =
            serde_json::from_value(value)
                .map_err(|e| format!("Invalid saved-region response: {e}"))?;
        let s = &mut self.cofactor_browser;
        s.saved_feature_regions = report
            .set
            .map(|set| set.regions)
            .unwrap_or_default()
            .into_iter()
            .filter(|region| {
                region.evidence.iter().any(|e| {
                    e.source_kind == "promoter_cofactor_query"
                        && e.feature_or_window_id
                            .as_deref()
                            .is_some_and(|id| id.starts_with("hit:"))
                })
            })
            .collect();
        s.saved_feature_regions
            .sort_by(|a, b| a.region_id.cmp(&b.region_id));
        let selected_id = s
            .saved_feature_region
            .as_ref()
            .map(|r| r.region_id.as_str());
        s.saved_feature_region = s
            .saved_feature_regions
            .iter()
            .find(|r| Some(r.region_id.as_str()) == selected_id)
            .or_else(|| s.saved_feature_regions.first())
            .cloned();
        s.feature_preview = None;
        s.feature_confirmed = false;
        s.status = format!(
            "Loaded {} saved motif region(s); preview again before attachment.",
            s.saved_feature_regions.len()
        );
        Ok(())
    }

    fn run_cofactor_feature_action(&mut self, apply: bool) -> Result<(), String> {
        let s = &self.cofactor_browser;
        if apply && !s.feature_confirmed {
            return Err("Confirm the previewed DNA annotation before attaching".into());
        }
        let region = s
            .saved_feature_region
            .as_ref()
            .ok_or("Save a retained motif hit first")?;
        let request = gentle_protocol::GenomicRegionFeatureRequest {
            set_id: "promoter_cofactors".into(),
            region_id: region.region_id.clone(),
            seq_id: s.feature_seq_id.clone(),
            expected_region_content_sha256: region.content_sha256.clone(),
            catalog_path: None,
            cache_dir: None,
            expected_approval_sha256: if apply {
                Some(
                    s.feature_preview
                        .as_ref()
                        .ok_or("Preview the DNA annotation first")?
                        .approval_sha256
                        .clone(),
                )
            } else {
                None
            },
        };
        self.cofactor_browser.feature_preview = None;
        self.cofactor_browser.feature_confirmed = false;
        let operation = if apply {
            Operation::MaterializeGenomicRegionFeature { request }
        } else {
            Operation::PreviewGenomicRegionFeature { request }
        };
        let (value, changed) =
            self.execute_shared_shell_command_json(&ShellCommand::GenomicRegions { operation })?;
        let report: gentle_protocol::GenomicRegionOperationReport =
            serde_json::from_value(value).map_err(|e| format!("Invalid feature response: {e}"))?;
        if changed {
            self.refresh_sequence_windows_from_engine_state();
            self.cofactor_browser.status = "Attached motif evidence to DNA. Undo is available; this does not establish occupancy or regulation.".into();
        } else {
            self.cofactor_browser.feature_preview = report.feature_materialization;
            self.cofactor_browser.status =
                "Reference verified; review the proposed annotation before attaching.".into();
        }
        Ok(())
    }

    fn start_cofactor_query(&mut self, query: CofactorQuery, anchor_id: Option<u64>) {
        let state = &mut self.cofactor_browser;
        if state.task.is_some() {
            return;
        }
        let mut request = state.request.clone();
        request.duckdb_executable = (!state.duckdb_executable.trim().is_empty())
            .then(|| state.duckdb_executable.trim().into());
        request.query = query;
        request.motif =
            (!state.motif_filter.trim().is_empty()).then(|| state.motif_filter.trim().into());
        request.anchor_id = anchor_id;
        if !matches!(query, CofactorQuery::Rankings | CofactorQuery::AnchorDetail) {
            request.motif = None;
            request.distance_band = None;
        }
        if matches!(query, CofactorQuery::Inspect | CofactorQuery::Candidates) {
            request.presence_threshold = 0.0;
        }
        request.region = None;
        request.gene_id = None;
        request.source_species = None;
        request.max_q_value = None;
        if query == CofactorQuery::Rankings {
            request.source_species = (!state.species_filter.trim().is_empty())
                .then(|| state.species_filter.trim().into());
            request.max_q_value = state.q_filter.then_some(state.q_value);
        } else if matches!(query, CofactorQuery::Anchors | CofactorQuery::Promoters) {
            request.gene_id =
                (!state.gene_filter.trim().is_empty()).then(|| state.gene_filter.trim().into());
            if !state.chromosome.trim().is_empty() {
                let (Ok(start), Ok(end)) = (state.start.parse(), state.end.parse()) else {
                    state.status = "Enter BED start/end coordinates".into();
                    return;
                };
                request.region = Some(CofactorRegion {
                    chromosome: state.chromosome.trim().into(),
                    start_0based: start,
                    end_0based_exclusive: end,
                });
            }
        }
        if let Err(e) = crate::promoter_cofactors::validate_request(&request) {
            state.status = e;
            return;
        }
        let (tx, rx) = mpsc::channel();
        let engine = self.engine.clone();
        std::thread::spawn(move || {
            let result =
                crate::background_engine::execute_on_engine_snapshot(&engine, move |snapshot| {
                    snapshot.apply(Operation::QueryPromoterCofactors { request })
                });
            let _ = tx.send(result);
        });
        state.task = Some(rx);
        state.submitted_form = Some(state.form_snapshot());
        state.status = "Verifying package and querying bounded evidence...".into();
        state.report = None;
        state.report_json.clear();
    }

    pub(super) fn render_promoter_cofactor_browser(&mut self, ctx: &egui::Context) {
        if let Some(receiver) = &self.cofactor_browser.task {
            match receiver.try_recv() {
                Ok(result) => {
                    self.cofactor_browser.task = None;
                    match result {
                        Ok(result) => {
                            if let Some(report) = result.promoter_cofactors {
                                self.cofactor_browser.status = format!(
                                    "{:?}: {}",
                                    report.availability,
                                    report.diagnostic.as_deref().unwrap_or(&report.report_id)
                                );
                                self.cofactor_browser.report_json =
                                    serde_json::to_string_pretty(&report).unwrap_or_default();
                                self.cofactor_browser.report = Some(report);
                            }
                        }
                        Err(e) => self.cofactor_browser.status = e.to_string(),
                    }
                }
                Err(mpsc::TryRecvError::Disconnected) => {
                    self.cofactor_browser.task = None;
                    self.cofactor_browser.status =
                        "Package query worker disconnected; no result published".into();
                }
                Err(mpsc::TryRecvError::Empty) => {
                    ctx.request_repaint_after(Duration::from_millis(100))
                }
            }
        }
        if !self.cofactor_browser.open {
            return;
        }
        let mut open = true;
        let viewport_id = ViewportId::from_hash_of("promoter_cofactor_browser");
        let spec = self.hosted_window_spec_for_viewport(
            "Promoter Cofactors",
            egui::Id::new("promoter_cofactors"),
            viewport_id,
            egui::vec2(1250.0, 850.0),
            egui::vec2(700.0, 450.0),
        );
        let mut action = None;
        let mut navigate = None;
        let mut region_action = None;
        let mut feature_action = None;
        let mut reload_features = false;
        let mut sequence_ids = self
            .engine
            .read()
            .ok()
            .map(|engine| engine.state().sequences.keys().cloned().collect::<Vec<_>>())
            .unwrap_or_default();
        sequence_ids.sort();
        crate::egui_compat::show_hosted_window(ctx, &spec, &mut open, |ui| {
            let s = &mut self.cofactor_browser;
            ui.horizontal_wrapped(|ui| {
                reload_features = ui.button("Reload saved motifs").clicked();
                if !s.saved_feature_regions.is_empty() {
                    egui::ComboBox::from_id_salt("cofactor_saved_feature")
                        .selected_text(
                            s.saved_feature_region
                                .as_ref()
                                .map(|r| r.label.as_deref().unwrap_or(&r.region_id))
                                .unwrap_or("Saved motif"),
                        )
                        .show_ui(ui, |ui| {
                            for region in &s.saved_feature_regions {
                                if ui
                                    .selectable_label(
                                        s.saved_feature_region
                                            .as_ref()
                                            .is_some_and(|r| r.region_id == region.region_id),
                                        region.label.as_deref().unwrap_or(&region.region_id),
                                    )
                                    .clicked()
                                {
                                    s.saved_feature_region = Some(region.clone());
                                    s.feature_preview = None;
                                    s.feature_confirmed = false;
                                }
                            }
                        });
                }
            });
            if let Some(region) = &s.saved_feature_region {
                ui.collapsing("Saved motif: DNA annotation", |ui| {
                    if let Some(label) = &region.label {
                        ui.strong(label);
                    }
                    ui.label(format!(
                        "{} | {}:[{}, {})",
                        region.region_id,
                        region.interval.reference.contig_name,
                        region.interval.start_0based,
                        region.interval.end_0based_exclusive
                    ));
                    for evidence in &region.evidence {
                        ui.label(&evidence.evidence_statement);
                    }
                    let before = s.feature_seq_id.clone();
                    egui::ComboBox::from_id_salt("cofactor_feature_sequence")
                        .selected_text(if s.feature_seq_id.is_empty() {
                            "Target DNA"
                        } else {
                            &s.feature_seq_id
                        })
                        .show_ui(ui, |ui| {
                            for id in &sequence_ids {
                                ui.selectable_value(&mut s.feature_seq_id, id.clone(), id);
                            }
                        });
                    if before != s.feature_seq_id {
                        s.feature_preview = None;
                        s.feature_confirmed = false;
                    }
                    if ui
                        .add_enabled(
                            !s.feature_seq_id.is_empty(),
                            egui::Button::new("Preview DNA annotation"),
                        )
                        .clicked()
                    {
                        feature_action = Some(false);
                    }
                    if let Some(preview) = &s.feature_preview {
                        ui.label(format!(
                            "{} | local [{}, {}) | {:?}",
                            preview.catalog_entry_id,
                            preview.projection.local_start_0based,
                            preview.projection.local_end_0based_exclusive,
                            preview.projection.local_strand
                        ));
                        ui.small(
                            "Motif association evidence only; not occupancy or causal regulation.",
                        );
                        ui.checkbox(
                            &mut s.feature_confirmed,
                            "I approve this annotation on the selected DNA",
                        );
                        if ui
                            .add_enabled(
                                s.feature_confirmed,
                                egui::Button::new("Attach DNA annotation"),
                            )
                            .clicked()
                        {
                            feature_action = Some(true);
                        }
                    }
                });
            }
            ui.add_enabled_ui(s.task.is_none(), |ui| {
                ui.horizontal_wrapped(|ui| {
                    ui.label("Package");
                    let path_changed = ui
                        .add(
                            egui::TextEdit::singleline(&mut s.request.package_path)
                                .desired_width(450.0),
                        )
                        .changed();
                    if ui.button("Browse...").clicked()
                        && let Some(path) = rfd::FileDialog::new().pick_folder()
                    {
                        s.request.package_path = path.display().to_string();
                        s.report = None;
                        s.report_json.clear();
                    }
                    ui.label("Assembly");
                    let assembly_changed = ui
                        .add(
                            egui::TextEdit::singleline(&mut s.request.assembly).desired_width(90.0),
                        )
                        .changed();
                    if path_changed || assembly_changed {
                        s.report = None;
                        s.report_json.clear();
                        s.status.clear();
                    }
                    if ui.button("Inspect package").clicked() {
                        action = Some((CofactorQuery::Inspect, None));
                    }
                });
                ui.horizontal_wrapped(|ui| {
                    ui.label("Motif / TF name");
                    ui.add(egui::TextEdit::singleline(&mut s.motif_filter).desired_width(125.0));
                    egui::ComboBox::from_id_salt("cofactor_band")
                        .selected_text(s.request.distance_band.as_deref().unwrap_or("All bands"))
                        .show_ui(ui, |ui| {
                            ui.selectable_value(&mut s.request.distance_band, None, "All bands");
                            for band in [
                                "overlap",
                                "adjacent_0_5",
                                "gap_6_20",
                                "gap_21_50",
                                "gap_51_100",
                                "gap_101_150",
                            ] {
                                ui.selectable_value(
                                    &mut s.request.distance_band,
                                    Some(band.into()),
                                    band,
                                );
                            }
                        });
                    egui::ComboBox::from_id_salt("cofactor_rank")
                        .selected_text(format!("{:?}", s.request.ranking))
                        .show_ui(ui, |ui| {
                            for (rank, label) in [
                                (CofactorRanking::TaEnriched, "TA enriched"),
                                (CofactorRanking::DnEnriched, "DN enriched"),
                                (CofactorRanking::TaDepleted, "TA depleted"),
                                (CofactorRanking::DnDepleted, "DN depleted"),
                                (
                                    CofactorRanking::IsoformDifference,
                                    "Strongest isoform difference",
                                ),
                            ] {
                                ui.selectable_value(&mut s.request.ranking, rank, label);
                            }
                        });
                    ui.label("Source species");
                    ui.add(egui::TextEdit::singleline(&mut s.species_filter).desired_width(130.0));
                    ui.checkbox(&mut s.q_filter, "BH q <=");
                    ui.add(
                        egui::DragValue::new(&mut s.q_value)
                            .range(0.0..=1.0)
                            .speed(0.01),
                    );
                    if ui.button("Rank motifs").clicked() {
                        action = Some((CofactorQuery::Rankings, None));
                    }
                });
                ui.horizontal_wrapped(|ui| {
                    ui.label("Promoter gene ID");
                    ui.add(egui::TextEdit::singleline(&mut s.gene_filter).desired_width(155.0));
                    ui.label("Chr");
                    ui.add(egui::TextEdit::singleline(&mut s.chromosome).desired_width(45.0));
                    ui.label("BED [start,end)");
                    ui.add(egui::TextEdit::singleline(&mut s.start).desired_width(95.0));
                    ui.add(egui::TextEdit::singleline(&mut s.end).desired_width(95.0));
                    if ui.button("Find anchors").clicked() {
                        action = Some((CofactorQuery::Anchors, None));
                    }
                    if ui.button("Find promoters").clicked() {
                        action = Some((CofactorQuery::Promoters, None));
                    }
                });
                ui.collapsing("Query limits", |ui| {
                    ui.horizontal_wrapped(|ui| {
                        ui.label("DuckDB executable");
                        ui.add(
                            egui::TextEdit::singleline(&mut s.duckdb_executable)
                                .desired_width(350.0)
                                .hint_text("GENTLE_DUCKDB_BIN or duckdb on PATH"),
                        );
                    });
                    ui.horizontal_wrapped(|ui| {
                        ui.label("Combined row limit");
                        ui.add(egui::DragValue::new(&mut s.request.max_rows).range(1..=2000));
                        ui.label("Timeout (s)");
                        ui.add(egui::DragValue::new(&mut s.request.timeout_seconds).range(1..=120));
                        ui.label("Presence score >=");
                        ui.add(
                            egui::DragValue::new(&mut s.request.presence_threshold)
                                .range(-1.0..=1000.0)
                                .speed(0.1),
                        );
                    });
                });
            });
            ui.horizontal_wrapped(|ui| {
                if s.task.is_some() {
                    ui.spinner();
                }
                ui.label(&s.status);
                if ui
                    .add_enabled(
                        !s.report_json.is_empty(),
                        egui::Button::new("Copy report JSON"),
                    )
                    .clicked()
                {
                    ui.ctx().copy_text(s.report_json.clone());
                }
                if ui
                    .add_enabled(s.report.is_some(), egui::Button::new("Copy exact request"))
                    .clicked()
                    && let Some(report) = &s.report
                    && let Ok(json) = serde_json::to_string_pretty(&report.request)
                {
                    ui.ctx().copy_text(json);
                }
            });
            if s.report_is_stale() {
                ui.colored_label(ui.visuals().warn_fg_color,
                    "Form changed: results below belong to the previous request. Rerun the query to update them. Copy/save actions use the displayed report, not the edited form.");
            }
            ui.horizontal_wrapped(|ui| {
                ui.label("Optional loaded sequence ID for saved-region projection");
                ui.text_edit_singleline(&mut s.projection_seq_id);
            });
            ui.separator();
            egui::ScrollArea::both()
                .id_salt("cofactor_report")
                .auto_shrink([false, false])
                .show(ui, |ui| {
                let Some(report) = &s.report else { return; };
                ui.label(format!("Displayed result: {:?}; motif {}; band {}; presence threshold >= {}",
                    report.request.query, report.request.motif.as_deref().unwrap_or("not filtered"),
                    report.request.distance_band.as_deref().unwrap_or("all"), report.request.presence_threshold));
                ui.label("Presence threshold changes only the selected-site presence flag. Cohort statistics and counts at -1 / 0 remain original. Unavailable is not zero.");
                if let Some(c) = &report.coverage {
                    ui.label(format!(
                        "{} | {} | source floor {} | positive >= {} | positional detail: {} motifs",
                        c.assembly, c.retention, c.source_score_floor, c.positive_threshold,
                        c.detailed_motif_ids.len()
                    ));
                    ui.collapsing("Coverage, provenance and interpretation", |ui| {
                        ui.label(format!("Package: {}", report.request.package_path));
                        ui.label(format!("Chromosomes: {}", c.chromosomes.join(", ")));
                        ui.monospace(format!(
                            "Manifest SHA-256: {}",
                            report.package_manifest_sha256.as_deref().unwrap_or("unavailable")
                        ));
                        ui.label(format!("H3K4me3 model effects: {}", c.h3k4me3_model_effects));
                        for (key, value) in &c.score_configuration {
                            ui.label(format!("{key}: {value}"));
                        }
                        for text in &report.non_claims {
                            ui.label(text);
                        }
                    });
                    if matches!(report.request.query, CofactorQuery::Inspect | CofactorQuery::Candidates) {
                        egui::Grid::new("cofactor_candidates").striped(true).show(ui, |ui| {
                            for heading in ["Requested factor", "Provenance group", "Matrices", "Availability"] {
                                ui.strong(heading);
                            }
                            ui.end_row();
                            for row in &c.requested_candidates {
                                ui.label(&row.gene);
                                ui.label(&row.group);
                                ui.label(row.motif_ids.join(", "));
                                ui.label(&row.status);
                                ui.end_row();
                            }
                        });
                    }
                }
                if !report.rankings.is_empty() {
                    ui.label("Cohort: included TP73 anchors overlapping extended regulatory promoters; statistics are not recomputed for a selected gene or region.");
                    if report.more_rankings_available {
                        ui.label(format!("First {} ranked rows; additional matches are not displayed.", report.rankings.len()));
                    }
                    egui::Grid::new("cofactor_rankings").striped(true).show(ui, |ui| {
                        for heading in [
                            "Motif / detail", "Band", "TA OR [95% CI]", "DN OR [95% CI]",
                            "TA/DN OR [95% CI]", "BH q: TA / DN / difference",
                            "Frequency >=0 (positive/eligible)", "Source species", "Support / status",
                        ] {
                            ui.strong(heading);
                        }
                        ui.end_row();
                        for row in &report.rankings {
                            let detail = report.coverage.as_ref()
                                .is_some_and(|c| c.detailed_motif_ids.contains(&row.motif_id));
                            let label = format!("{} {} ({})", row.motif_name, row.motif_id,
                                if detail { "positions" } else { "overview only" });
                            if ui.selectable_label(s.motif_filter == row.motif_id, label).clicked() {
                                s.motif_filter = row.motif_id.clone();
                            }
                            ui.label(&row.distance_band);
                            for (v, lo, hi) in [
                                (row.ta_adjusted_odds_ratio, row.ta_confidence_interval_95_lower, row.ta_confidence_interval_95_upper),
                                (row.dn_adjusted_odds_ratio, row.dn_confidence_interval_95_lower, row.dn_confidence_interval_95_upper),
                                (row.ta_vs_dn_odds_ratio_ratio, row.confidence_interval_95_lower, row.confidence_interval_95_upper),
                            ] {
                                ui.label(format!("{} [{}, {}]", number(v), number(lo), number(hi)));
                            }
                            ui.label(format!("{} / {} / {}", number(row.ta_q_value_bh_tax_group),
                                number(row.dn_q_value_bh_tax_group), number(row.q_value_bh_tax_group)));
                            ui.label(format!("{}% ({}/{})", number(row.positive_anchor_fraction.map(|v| v * 100.0)),
                                row.anchors_positive, row.anchors_total));
                            ui.label(&row.source_species);
                            ui.label(format!("{}; TA={} DN={} difference={}", row.class_support_flag,
                                row.ta_evaluation_status, row.dn_evaluation_status, row.evaluation_status));
                            ui.end_row();
                        }
                    });
                }
                for anchor in &report.anchors {
                    ui.horizontal_wrapped(|ui| {
                        ui.strong(format!("Anchor {}: {}:[{}, {}) | TP73 score {:.5}",
                            anchor.anchor_id, anchor.chrom, anchor.anchor_start, anchor.anchor_end, anchor.anchor_score));
                        if ui.add_enabled(s.task.is_none() && !s.motif_filter.is_empty(),
                            egui::Button::new("Cofactor detail")).clicked()
                        {
                            action = Some((CofactorQuery::AnchorDetail, Some(anchor.anchor_id)));
                        }
                        if ui.button("Open region...").clicked() {
                            navigate = Some((anchor.chrom.clone(), anchor.anchor_start, anchor.anchor_end, report.coverage.clone()));
                        }
                        region_buttons(ui, CofactorRegionTarget::Anchor { anchor_id: anchor.anchor_id }, &mut region_action);
                    });
                    egui::Grid::new(("cofactor_support", anchor.anchor_id)).striped(true).show(ui, |ui| {
                        for heading in ["Sample / series", "TP73 support / max depth", "Control support / max depth"] {
                            ui.strong(heading);
                        }
                        ui.end_row();
                        for series in ["saos2_TA", "saos2_DN", "skmel29_2_TA", "skmel29_2_DN"] {
                            ui.label(series);
                            for kind in ["tp73", "negative_control"] {
                                let support = anchor.source_fields.get(&format!("supported_{kind}_{series}"))
                                    .map(Value::to_string).unwrap_or_else(|| "unavailable".into());
                                let depth = anchor.source_fields.get(&format!("depth_{kind}_{series}"))
                                    .and_then(Value::as_f64);
                                ui.label(format!("{support} / {}", number(depth)));
                            }
                            ui.end_row();
                        }
                    });
                }
                if !report.details.is_empty() {
                    egui::Grid::new("cofactor_detail").striped(true).show(ui, |ui| {
                        for h in ["Band", "Full BED interval", "Score / strand", "Gap / genomic side",
                            "Counts >=-1 / >=0", "Presence at displayed threshold", "Selected evidence"]
                        {
                            ui.strong(h);
                        }
                        ui.end_row();
                        for d in &report.details {
                            ui.label(&d.distance_band);
                            ui.label(match (d.hit_start, d.hit_end) {
                                (Some(a), Some(b)) => format!("[{a}, {b})"),
                                _ => "no retained locus".into(),
                            });
                            ui.label(format!("{} / {}", number(d.best_score), d.best_strand.as_deref().unwrap_or("unavailable")));
                            ui.label(format!("{} / {}", d.interval_distance_bp.map(|v| v.to_string()).unwrap_or_else(|| "unavailable".into()),
                                d.genomic_side.as_deref().unwrap_or("unavailable")));
                            ui.label(format!("{} / {}", d.n_source_loci, d.n_score_zero_loci));
                            ui.label(d.present_at_requested_threshold.to_string());
                            if d.hit_start.is_some() {
                                ui.horizontal_wrapped(|ui| {
                                    region_buttons(ui, CofactorRegionTarget::Hit { anchor_id: d.anchor_id,
                                        motif_id: d.motif_id.clone(), distance_band: d.distance_band.clone() }, &mut region_action);
                                    if ui.button("Open hit region...").clicked()
                                        && let Some(anchor) = report.anchors.iter().find(|a| a.anchor_id == d.anchor_id)
                                        && let (Some(start), Some(end)) = (d.hit_start, d.hit_end)
                                    {
                                        navigate = Some((anchor.chrom.clone(), start, end, report.coverage.clone()));
                                    }
                                });
                            } else { ui.label("No retained locus to save"); }
                            ui.end_row();
                        }
                    });
                }
                for p in &report.promoters {
                    ui.horizontal_wrapped(|ui| {
                        ui.strong(format!("Promoter {} | {}:[{}, {})", p.regulatory_feature_id, p.chrom, p.extended_start, p.extended_end));
                        if ui.button("Open region...").clicked() {
                            navigate = Some((p.chrom.clone(), p.extended_start, p.extended_end, report.coverage.clone()));
                        }
                        region_buttons(ui, CofactorRegionTarget::Promoter { regulatory_feature_id: p.regulatory_feature_id.clone() }, &mut region_action);
                    });
                    for g in &p.gene_links {
                        ui.label(format!("{} | {} | {}", g.gene_id, g.link_source, g.annotation_release));
                    }
                    let anchors = report.memberships.iter()
                        .filter(|m| m.regulatory_feature_id == p.regulatory_feature_id)
                        .map(|m| m.anchor_id.to_string()).collect::<Vec<_>>();
                    ui.label(format!("Member anchors: {}", anchors.join(", ")));
                }
            });
        });
        self.cofactor_browser.open = open;
        if reload_features && let Err(error) = self.reload_cofactor_feature_regions() {
            self.cofactor_browser.status = error;
        }
        if let Some((kind, anchor)) = action {
            self.start_cofactor_query(kind, anchor);
        }
        if let Some((target, copy)) = region_action
            && let Some(request) = self.cofactor_browser.region_request(target)
        {
            if copy {
                if let Ok(json) = serde_json::to_string(&request) {
                    ctx.copy_text(json);
                }
            } else {
                let command = ShellCommand::GenomicRegions {
                    operation: Operation::CaptureGenomicRegion { request },
                };
                self.cofactor_browser.status = match self
                    .execute_shared_shell_command_json(&command)
                {
                    Ok((value, _)) => {
                        self.cofactor_browser.feature_preview = None;
                        self.cofactor_browser.feature_confirmed = false;
                        self.cofactor_browser.saved_feature_region = serde_json::from_value::<
                            gentle_protocol::GenomicRegionOperationReport,
                        >(
                            value
                        )
                        .ok()
                        .and_then(|report| report.region)
                        .filter(|region| {
                            region.evidence.iter().any(|e| {
                                e.feature_or_window_id
                                    .as_deref()
                                    .is_some_and(|id| id.starts_with("hit:"))
                            })
                        });
                        "Saved in region set promoter_cofactors. No DNA annotations or coverage were imported. Stranded hits can be previewed in Saved motif: DNA annotation.".into()
                    }
                    Err(error) => error,
                };
            }
        }
        if let Some(apply) = feature_action
            && let Err(error) = self.run_cofactor_feature_action(apply)
        {
            self.cofactor_browser.status = error;
        }
        if let Some((chrom, start, end, Some(coverage))) = navigate {
            self.open_reference_genome_retrieve_dialog();
            self.genome_id = coverage
                .score_configuration
                .get("genome_id")
                .and_then(Value::as_str)
                .unwrap_or_default()
                .into();
            self.genome_chromosome = chrom;
            self.genome_start_1based = (start + 1).to_string();
            self.genome_end_1based = end.to_string();
            self.genome_retrieve_status = format!(
                "Promoter-cofactor region on {}: BED [start,end) converted to 1-based inclusive. Select the matching prepared reference before explicit extraction; no genome was downloaded.",
                coverage.assembly
            );
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn feature_gui_frame(
        app: &mut GENtleApp,
        ctx: &egui::Context,
        events: Vec<egui::Event>,
    ) -> Vec<(String, egui::Pos2)> {
        fn collect(shape: &egui::epaint::Shape, text: &mut Vec<(String, egui::Pos2)>) {
            match shape {
                egui::epaint::Shape::Text(t) => {
                    text.push((t.galley.text().into(), t.pos + t.galley.size() * 0.5))
                }
                egui::epaint::Shape::Vec(shapes) => {
                    for shape in shapes {
                        collect(shape, text);
                    }
                }
                _ => {}
            }
        }
        ctx.begin_pass(egui::RawInput {
            screen_rect: Some(egui::Rect::from_min_size(
                egui::Pos2::ZERO,
                egui::vec2(1400.0, 1000.0),
            )),
            events,
            ..Default::default()
        });
        app.render_promoter_cofactor_browser(ctx);
        let mut text = vec![];
        for shape in crate::egui_compat::end_test_pass(ctx).shapes {
            collect(&shape.shape, &mut text);
        }
        text
    }

    fn click_feature_gui_text(app: &mut GENtleApp, ctx: &egui::Context, label: &str) {
        let mut text = vec![];
        for _ in 0..3 {
            text = feature_gui_frame(app, ctx, vec![]);
        }
        let pos = text
            .iter()
            .find(|(t, _)| t == label)
            .unwrap_or_else(|| panic!("Missing {label}: {text:?}"))
            .1;
        for pressed in [true, false] {
            feature_gui_frame(
                app,
                ctx,
                vec![
                    egui::Event::PointerMoved(pos),
                    egui::Event::PointerButton {
                        pos,
                        button: egui::PointerButton::Primary,
                        pressed,
                        modifiers: egui::Modifiers::NONE,
                    },
                ],
            );
        }
    }

    #[test]
    fn promoter_cofactors_rendered_gui_clicks_preview_confirm_and_attach() {
        let (_dir, engine, request) = crate::promoter_cofactors::tests::feature_handoff_engine("+");
        let region = engine.genomic_region_store_snapshot().unwrap().sets[0].regions[0].clone();
        let mut app = GENtleApp::default();
        *app.engine.write().unwrap() = engine;
        app.cofactor_browser.open = true;
        app.cofactor_browser.saved_feature_region = Some(region);
        app.cofactor_browser.feature_seq_id = request.seq_id;
        let ctx = egui::Context::default();
        ctx.style_mut_of(egui::Theme::Dark, |style| style.animation_time = 0.0);
        ctx.style_mut_of(egui::Theme::Light, |style| style.animation_time = 0.0);
        app.cofactor_browser.saved_feature_region = None;
        click_feature_gui_text(&mut app, &ctx, "Reload saved motifs");
        assert!(app.cofactor_browser.saved_feature_region.is_some());
        click_feature_gui_text(&mut app, &ctx, "Saved motif: DNA annotation");
        click_feature_gui_text(&mut app, &ctx, "Preview DNA annotation");
        assert!(
            app.cofactor_browser.feature_preview.is_some(),
            "{}",
            app.cofactor_browser.status
        );
        click_feature_gui_text(&mut app, &ctx, "Attach DNA annotation");
        assert!(
            app.engine.read().unwrap().state().sequences["demo"]
                .features()
                .is_empty()
        );
        click_feature_gui_text(
            &mut app,
            &ctx,
            "I approve this annotation on the selected DNA",
        );
        assert!(app.cofactor_browser.feature_confirmed);
        click_feature_gui_text(&mut app, &ctx, "Attach DNA annotation");
        assert_eq!(
            app.engine.read().unwrap().state().sequences["demo"]
                .features()
                .len(),
            1
        );
    }

    #[test]
    fn promoter_cofactors_gui_feature_preview_and_confirm_use_shared_engine() {
        let (_dir, engine, request) = crate::promoter_cofactors::tests::feature_handoff_engine("-");
        let region = engine.genomic_region_store_snapshot().unwrap().sets[0].regions[0].clone();
        let mut app = GENtleApp::default();
        *app.engine.write().unwrap() = engine;
        app.cofactor_browser.saved_feature_region = Some(region);
        app.cofactor_browser.feature_seq_id = request.seq_id;
        assert!(app.run_cofactor_feature_action(true).is_err());
        app.run_cofactor_feature_action(false).unwrap();
        let preview = app.cofactor_browser.feature_preview.as_ref().unwrap();
        assert_eq!(preview.projection.local_start_0based, 52);
        assert_eq!(
            preview.projection.local_strand,
            gentle_protocol::GenomicRegionStrand::Minus
        );
        assert!(
            app.engine.read().unwrap().state().sequences["demo"]
                .features()
                .is_empty()
        );
        app.cofactor_browser.feature_confirmed = true;
        app.run_cofactor_feature_action(true).unwrap();
        assert_eq!(
            app.engine.read().unwrap().state().sequences["demo"]
                .features()
                .len(),
            1
        );
        assert!(app.cofactor_browser.feature_preview.is_none());
        assert!(!app.cofactor_browser.feature_confirmed);
    }

    #[test]
    fn promoter_cofactors_gui_feature_target_change_rejects_stale_preview() {
        let (_dir, engine, request) = crate::promoter_cofactors::tests::feature_handoff_engine("+");
        let region = engine.genomic_region_store_snapshot().unwrap().sets[0].regions[0].clone();
        let mut app = GENtleApp::default();
        *app.engine.write().unwrap() = engine;
        app.cofactor_browser.saved_feature_region = Some(region);
        app.cofactor_browser.feature_seq_id = request.seq_id;
        app.run_cofactor_feature_action(false).unwrap();
        app.cofactor_browser.feature_seq_id = "missing".into();
        app.cofactor_browser.feature_confirmed = true;
        assert!(app.run_cofactor_feature_action(true).is_err());
        assert!(app.cofactor_browser.feature_preview.is_none());
        assert!(
            app.engine.read().unwrap().state().sequences["demo"]
                .features()
                .is_empty()
        );
    }

    #[test]
    fn promoter_cofactors_form_edits_mark_results_stale_but_keep_copy_and_capture_bound() {
        let mut browser = CofactorBrowser::default();
        browser.submitted_form = Some(browser.form_snapshot());
        let report = crate::promoter_cofactors::tests::handoff_report();
        let original = serde_json::to_value(&report.request).unwrap();
        browser.report = Some(report);
        assert!(!browser.report_is_stale());
        browser.request.presence_threshold = 10.0;
        assert!(browser.report_is_stale());
        assert_eq!(
            serde_json::to_value(&browser.report.as_ref().unwrap().request).unwrap(),
            original
        );
        let request = browser
            .region_request(CofactorRegionTarget::Anchor { anchor_id: 1 })
            .unwrap();
        let gentle_protocol::GenomicRegionCaptureSource::PromoterCofactor { report, .. } =
            request.source
        else {
            panic!("capture source")
        };
        assert_eq!(report.request.presence_threshold, 0.0);
        browser.request.presence_threshold = 0.0;
        assert!(!browser.report_is_stale());
        browser.motif_filter = "another motif".into();
        assert!(browser.report_is_stale());
    }

    #[test]
    fn promoter_cofactors_browser_opens_without_genome_or_runtime() {
        let mut app = GENtleApp::default();
        let before = serde_json::to_value(app.engine.read().unwrap().state()).unwrap();
        let ctx = egui::Context::default();
        assert!(app.collect_command_palette_entries().iter().any(|entry| {
            entry.title == "Promoter Cofactors"
                && matches!(entry.action, CommandPaletteAction::OpenPromoterCofactors)
        }));
        app.execute_command_palette_action(&ctx, CommandPaletteAction::OpenPromoterCofactors);
        assert!(app.cofactor_browser.open);
        ctx.begin_pass(egui::RawInput::default());
        app.render_promoter_cofactor_browser(&ctx);
        let mut output = ctx.end_pass();
        output.textures_delta.clear();
        assert!(!output.shapes.is_empty());
        assert!(app.cofactor_browser.task.is_none());
        assert_eq!(
            before,
            serde_json::to_value(app.engine.read().unwrap().state()).unwrap()
        );
        assert_eq!(app.cofactor_browser.q_value, 0.05);
        app.cofactor_browser.submitted_form = Some(app.cofactor_browser.form_snapshot());
        app.cofactor_browser.report = Some(crate::promoter_cofactors::tests::handoff_report());
        app.cofactor_browser.request.presence_threshold = 5.0;
        assert!(app.cofactor_browser.report_is_stale());
        ctx.begin_pass(egui::RawInput::default());
        app.render_promoter_cofactor_browser(&ctx);
        let mut output = ctx.end_pass();
        output.textures_delta.clear();
        assert!(!output.shapes.is_empty());
        assert!(app.cofactor_browser.task.is_none());
        assert_eq!(
            before,
            serde_json::to_value(app.engine.read().unwrap().state()).unwrap()
        );
    }
}
