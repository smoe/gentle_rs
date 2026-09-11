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
        }
    }
}

fn number(value: Option<f64>) -> String {
    value
        .map(|v| format!("{v:.4}"))
        .unwrap_or_else(|| "unavailable".into())
}

impl GENtleApp {
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
        crate::egui_compat::show_hosted_window(ctx, &spec, &mut open, |ui| {
            let s = &mut self.cofactor_browser;
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
            ui.separator();
            egui::ScrollArea::both()
                .id_salt("cofactor_report")
                .auto_shrink([false, false])
                .show(ui, |ui| {
                let Some(report) = &s.report else { return; };
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
                            "Counts >=-1 / >=0", "Presence at requested threshold"]
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
        if let Some((kind, anchor)) = action {
            self.start_cofactor_query(kind, anchor);
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
    }
}
