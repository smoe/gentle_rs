//! Thin DNA-viewer workspace for the shared three-source TATA evidence report.

use super::*;
use gentle_protocol::tata_boxes::*;

#[derive(Debug, Clone)]
struct TataTask {
    receiver: Arc<Mutex<Receiver<Result<OpResult, EngineError>>>>,
    started: Instant,
    baseline_revision: u64,
    mutating: bool,
}

#[derive(Debug, Clone, Default)]
pub(super) struct TataBoxWorkspace {
    pub open: bool,
    request: TataBoxScreenRequest,
    epd_json: String,
    additional_tss_json: String,
    task: Option<TataTask>,
    report: Option<Arc<TataBoxScreenReport>>,
    report_revision: Option<u64>,
    selected: BTreeSet<String>,
    status: String,
}

impl MainAreaDna {
    pub(crate) fn open_tata_boxes(&mut self) {
        self.tata_ui.open = true;
        self.tata_ui.request.seq_id = self.seq_id.clone().unwrap_or_default();
    }

    #[cfg(test)]
    pub(crate) fn tata_workspace_open(&self) -> bool {
        self.tata_ui.open
    }

    fn tata_request(&self) -> Result<TataBoxScreenRequest, String> {
        let mut request = self.tata_ui.request.clone();
        request.seq_id = self.seq_id.clone().ok_or("No sequence")?;
        request.epd = if self.tata_ui.epd_json.trim().is_empty() {
            None
        } else {
            Some(
                serde_json::from_str(self.tata_ui.epd_json.trim())
                    .map_err(|e| format!("EPD: {e}"))?,
            )
        };
        request.additional_tss = if self.tata_ui.additional_tss_json.trim().is_empty() {
            vec![]
        } else {
            serde_json::from_str(self.tata_ui.additional_tss_json.trim())
                .map_err(|e| format!("TSS: {e}"))?
        };
        Ok(request)
    }

    fn tata_report_matches_request(&self, report: &TataBoxScreenReport) -> bool {
        self.tata_request().is_ok_and(|mut request| {
            if request.end_0based_exclusive.is_none() {
                request.end_0based_exclusive = self.dna.try_read().ok().map(|dna| dna.len());
            }
            request == report.request
        })
    }

    fn start_tata_operation(&mut self, operation: Operation, mutating: bool) {
        if self.tata_ui.task.is_some() {
            return;
        }
        let Some(engine) = self.engine.clone() else {
            self.tata_ui.status = "No engine".into();
            return;
        };
        let revision = match engine.try_read() {
            Ok(guard) => guard.structural_revision(),
            Err(_) => {
                self.tata_ui.status = Self::tr("tata.busy");
                return;
            }
        };
        let (tx, rx) = mpsc::channel();
        self.tata_ui.task = Some(TataTask {
            receiver: Arc::new(Mutex::new(rx)),
            started: Instant::now(),
            baseline_revision: revision,
            mutating,
        });
        self.tata_ui.status = Self::tr("tata.running");
        std::thread::spawn(move || {
            let result = if mutating {
                crate::background_engine::execute_on_engine_snapshot(&engine, |snapshot| {
                    snapshot.apply(operation)
                })
            } else {
                crate::background_engine::execute_read_only_operation_on_engine_snapshot(
                    &engine, operation,
                )
            };
            let _ = tx.send(result);
        });
    }

    pub(super) fn poll_tata_task(&mut self, ctx: &egui::Context) {
        let Some(task) = self.tata_ui.task.clone() else {
            return;
        };
        let outcome = match task.receiver.lock() {
            Ok(rx) => match rx.try_recv() {
                Ok(result) => Some(result),
                Err(TryRecvError::Empty) => None,
                Err(TryRecvError::Disconnected) => Some(Err(EngineError::new(
                    ErrorCode::Internal,
                    "TATA worker disconnected",
                ))),
            },
            Err(_) => Some(Err(EngineError::new(
                ErrorCode::Internal,
                "TATA worker channel unavailable",
            ))),
        };
        let Some(outcome) = outcome else {
            ctx.request_repaint_after(Duration::from_millis(100));
            return;
        };
        self.tata_ui.task = None;
        match outcome {
            Ok(result) if task.mutating => {
                self.handle_operation_success(result, task.started);
                self.tata_ui.report = None;
                self.tata_ui.selected.clear();
                self.tata_ui.status = Self::tr("tata.added");
            }
            Ok(mut result) => {
                let revision = self
                    .engine
                    .as_ref()
                    .and_then(|e| e.try_read().ok().map(|g| g.structural_revision()));
                if revision != Some(task.baseline_revision) {
                    self.tata_ui.report = None;
                    self.tata_ui.status = Self::tr("tata.stale");
                    return;
                }
                self.tata_ui.report = result.tata_box_screen.take().map(|r| Arc::new(*r));
                self.tata_ui.report_revision = revision;
                self.tata_ui.selected.clear();
                self.tata_ui.status = if self.tata_ui.report.is_some() {
                    Self::tr("tata.complete")
                } else {
                    "TATA operation returned no report".into()
                };
            }
            Err(error) => {
                self.tata_ui.report = None;
                self.tata_ui.status = error.to_string();
            }
        }
        ctx.request_repaint();
    }

    pub(super) fn render_tata_workspace(&mut self, ctx: &egui::Context) {
        if !self.tata_ui.open {
            return;
        }
        let mut open = true;
        let spec = crate::egui_compat::HostedWindowSpec::new(
            format!(
                "{}: {}",
                Self::tr("tata.title"),
                self.seq_id.as_deref().unwrap_or("")
            ),
            egui::Id::new(("tata_boxes", self.panel_scope_key())),
            Vec2::new(1050.0, 650.0),
            Vec2::new(620.0, 360.0),
        );
        crate::egui_compat::show_hosted_window(ctx, &spec, &mut open, |ui| {
            egui::ScrollArea::vertical()
                .id_salt("tata_body")
                .show(ui, |ui| {
                    let idle = self.tata_ui.task.is_none();
                    ui.add_enabled_ui(idle, |ui| {
                        ui.horizontal_wrapped(|ui| {
                            ui.checkbox(
                                &mut self.tata_ui.request.include_annotations,
                                Self::tr("tata.annotations"),
                            );
                            ui.checkbox(
                                &mut self.tata_ui.request.predict,
                                Self::tr("tata.predictions"),
                            );
                            ui.checkbox(
                                &mut self.tata_ui.request.scan_without_tss,
                                Self::tr("tata.without_tss"),
                            );
                        });
                        ui.horizontal_wrapped(|ui| {
                            ui.label(Self::tr("tata.matrix"));
                            ui.add(
                                egui::TextEdit::singleline(&mut self.tata_ui.request.motif_id)
                                    .desired_width(100.0),
                            );
                            ui.label(Self::tr("tata.score"));
                            ui.add(
                                egui::DragValue::new(&mut self.tata_ui.request.minimum_llr_bits)
                                    .speed(0.1),
                            );
                            ui.label(Self::tr("tata.offset"));
                            ui.add(egui::DragValue::new(
                                &mut self.tata_ui.request.minimum_tss_offset,
                            ));
                            ui.label("..");
                            ui.add(egui::DragValue::new(
                                &mut self.tata_ui.request.maximum_tss_offset,
                            ));
                        });
                        ui.horizontal_wrapped(|ui| {
                            if ui.button(Self::tr("tata.whole")).clicked() {
                                self.tata_ui.request.start_0based = 0;
                                self.tata_ui.request.end_0based_exclusive = None;
                            }
                            if ui.button(Self::tr("tata.selection")).clicked()
                                && let Some((start, end)) = self.current_selection_range_0based()
                            {
                                self.tata_ui.request.start_0based = start;
                                self.tata_ui.request.end_0based_exclusive = Some(end);
                            }
                            let end = self
                                .tata_ui
                                .request
                                .end_0based_exclusive
                                .or_else(|| self.dna.try_read().ok().map(|d| d.len()))
                                .unwrap_or(0);
                            ui.label(format!(
                                "{}..{} (1-based)",
                                self.tata_ui.request.start_0based + 1,
                                end
                            ));
                        });
                        egui::CollapsingHeader::new(Self::tr("tata.epd")).show(ui, |ui| {
                            ui.hyperlink_to(
                                "EPDnew",
                                "https://epd.expasy.org/epd/EPDnew_select.php",
                            );
                            ui.add(
                                egui::TextEdit::multiline(&mut self.tata_ui.epd_json)
                                    .code_editor()
                                    .desired_rows(4)
                                    .desired_width(f32::INFINITY),
                            );
                            if ui.button(Self::tr("tata.epd_template")).clicked() {
                                let source = TataBoxEpdSource {
                                    bed_path: "data/resources/epd/Hs_EPDnew_006_hg38.bed".into(),
                                    motifs_path: "data/resources/epd/promoter_motifs.txt".into(),
                                    assembly: "GRCh38".into(),
                                    taxon_id: 9606,
                                    release: "006".into(),
                                    source_url: "https://epd.expasy.org/ftp/epdnew/H_sapiens/006/"
                                        .into(),
                                    required: true,
                                    expected_bed_sha256: None,
                                    expected_motifs_sha256: None,
                                };
                                if let Ok(json) = serde_json::to_string_pretty(&source) {
                                    self.tata_ui.epd_json = json;
                                }
                            }
                        });
                        egui::CollapsingHeader::new(Self::tr("tata.tss")).show(ui, |ui| {
                            ui.add(
                                egui::TextEdit::multiline(&mut self.tata_ui.additional_tss_json)
                                    .code_editor()
                                    .desired_rows(3)
                                    .desired_width(f32::INFINITY),
                            );
                        });
                        ui.horizontal(|ui| {
                            if ui.button(Self::tr("tata.inspect")).clicked() {
                                match self.tata_request() {
                                    Ok(request) => self.start_tata_operation(
                                        Operation::ScreenTataBoxes {
                                            request,
                                            path: None,
                                        },
                                        false,
                                    ),
                                    Err(error) => self.tata_ui.status = error,
                                }
                            }
                            if ui.button(Self::tr("tata.copy_request")).clicked() {
                                match self.tata_request().and_then(|r| {
                                    serde_json::to_string(&r).map_err(|e| e.to_string())
                                }) {
                                    Ok(json) => ui.ctx().copy_text(format!(
                                        "promoters tata-screen {}",
                                        crate::engine_shell::shell_quote(&json)
                                    )),
                                    Err(error) => self.tata_ui.status = error,
                                }
                            }
                        });
                    });
                    if let Some(task) = &self.tata_ui.task {
                        ui.horizontal(|ui| {
                            ui.spinner();
                            ui.label(format!(
                                "{} ({:.1}s)",
                                Self::tr("tata.running"),
                                task.started.elapsed().as_secs_f32()
                            ));
                        });
                    }
                    ui.label(&self.tata_ui.status);
                    let Some(report) = self.tata_ui.report.clone() else {
                        return;
                    };
                    let revision = self
                        .engine
                        .as_ref()
                        .and_then(|e| e.try_read().ok().map(|g| g.structural_revision()));
                    let same_request = self.tata_report_matches_request(&report);
                    let current = same_request
                        && revision.is_some()
                        && revision == self.tata_ui.report_revision;
                    if !current {
                        ui.colored_label(egui::Color32::RED, Self::tr("tata.stale"));
                    }
                    ui.small(&report.non_claim);
                    ui.label(format!(
                        "EPD: {} | {} {}",
                        report.epd_status,
                        report.scored_windows,
                        Self::tr("tata.windows")
                    ));
                    for warning in &report.warnings {
                        ui.colored_label(egui::Color32::DARK_RED, warning);
                    }
                    ui.horizontal(|ui| {
                        if ui.button(Self::tr("tata.copy_report")).clicked()
                            && let Ok(json) = serde_json::to_string_pretty(&*report)
                        {
                            ui.ctx().copy_text(json);
                        }
                        if ui
                            .add_enabled(
                                idle && current && !self.tata_ui.selected.is_empty(),
                                egui::Button::new(Self::tr("tata.add")),
                            )
                            .clicked()
                        {
                            self.start_tata_operation(
                                Operation::MaterializeTataBoxFeatures {
                                    request: TataBoxMaterializeRequest {
                                        screen: report.request.clone(),
                                        expected_report_sha256: report.content_sha256.clone(),
                                        row_ids: self.tata_ui.selected.iter().cloned().collect(),
                                    },
                                },
                                true,
                            );
                        }
                    });
                    egui::ScrollArea::both()
                        .id_salt("tata_rows")
                        .max_height(400.0)
                        .auto_shrink([false, false])
                        .show_rows(ui, 28.0, report.rows.len(), |ui, range| {
                            for row in &report.rows[range] {
                                ui.horizontal(|ui| {
                                    let mut selected = self.tata_ui.selected.contains(&row.row_id);
                                    if ui
                                        .add_enabled(
                                            row.evidence_kind
                                                != TataBoxEvidenceKind::SourceAnnotation,
                                            egui::Checkbox::new(&mut selected, ""),
                                        )
                                        .changed()
                                    {
                                        if selected {
                                            self.tata_ui.selected.insert(row.row_id.clone());
                                        } else {
                                            self.tata_ui.selected.remove(&row.row_id);
                                        }
                                    }
                                    let kind = match row.evidence_kind {
                                        TataBoxEvidenceKind::SourceAnnotation => {
                                            Self::tr("tata.annotations")
                                        }
                                        TataBoxEvidenceKind::EpdClassification => format!(
                                            "EPD: {}",
                                            match row.tata_positive {
                                                Some(true) => "+",
                                                Some(false) => "-",
                                                None => "?",
                                            }
                                        ),
                                        TataBoxEvidenceKind::MotifPrediction => {
                                            Self::tr("tata.predictions")
                                        }
                                    };
                                    ui.label(kind);
                                    if ui
                                        .add_enabled(
                                            current,
                                            egui::Button::new(format!(
                                                "{}..{} {}",
                                                row.start_0based + 1,
                                                row.end_0based_exclusive,
                                                if row.reverse { "-" } else { "+" }
                                            )),
                                        )
                                        .clicked()
                                    {
                                        let _ = self.inspect_sequence_span_0based(
                                            row.start_0based,
                                            row.end_0based_exclusive,
                                            &row.label,
                                        );
                                    }
                                    ui.label(&row.label).on_hover_text(
                                        serde_json::to_string_pretty(row).unwrap_or_default(),
                                    );
                                    if let Some(score) = row.llr_bits {
                                        ui.label(format!("{score:.2} bits"));
                                    }
                                    for tss in &row.tss_associations {
                                        ui.label(format!(
                                            "{} {:+} bp",
                                            tss.tss_id, tss.signed_distance_bp
                                        ));
                                    }
                                });
                            }
                        });
                });
        });
        self.tata_ui.open = open;
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::engine::ProjectState;

    #[test]
    fn tata_gui_worker_is_read_only_and_uses_typed_request() {
        let _lock = crate::tf_motifs::test_registry_lock()
            .lock()
            .unwrap_or_else(|poisoned| poisoned.into_inner());
        crate::tf_motifs::reload_builtin_for_test();
        // Hand-crafted motif-bearing DNA; no experimental sequence.
        let dna = DNAsequence::from_sequence("CCCCTATAAAACCCC").unwrap();
        let mut state = ProjectState::default();
        state.sequences.insert("toy".into(), dna.clone());
        let engine = GentleEngine::from_state(state);
        let shared = Arc::new(RwLock::new(engine));
        let mut area = MainAreaDna::new(dna, Some("toy".into()), Some(shared.clone()));
        let before = serde_json::to_value(shared.read().unwrap().state()).unwrap();
        area.open_tata_boxes();
        area.tata_ui.request.scan_without_tss = true;
        let request = area.tata_request().unwrap();
        let expected = shared.read().unwrap().screen_tata_boxes(&request).unwrap();
        area.start_tata_operation(
            Operation::ScreenTataBoxes {
                request,
                path: None,
            },
            false,
        );
        let ctx = egui::Context::default();
        let deadline = Instant::now() + Duration::from_secs(10);
        while area.tata_ui.task.is_some() && Instant::now() < deadline {
            area.poll_tata_task(&ctx);
            std::thread::sleep(Duration::from_millis(10));
        }
        assert!(area.tata_ui.task.is_none());
        assert_eq!(area.tata_ui.report.as_deref(), Some(&expected));
        assert_eq!(
            serde_json::to_value(shared.read().unwrap().state()).unwrap(),
            before
        );
        for _ in 0..2 {
            let mut output = ctx.run_ui(egui::RawInput::default(), |ui| {
                area.render_tata_workspace(ui.ctx())
            });
            output.textures_delta.clear();
        }
        assert!(area.tata_ui.open);
        assert!(area.tata_report_matches_request(&expected));
        area.tata_ui.request.minimum_llr_bits += 1.0;
        assert!(!area.tata_report_matches_request(&expected));
        area.tata_ui.request.minimum_llr_bits -= 1.0;
        let mut partial = expected.clone();
        partial.request.end_0based_exclusive = Some(10);
        assert!(!area.tata_report_matches_request(&partial));
        area.tata_ui.request.end_0based_exclusive = Some(10);
        assert!(area.tata_report_matches_request(&partial));
        assert_eq!(
            serde_json::to_value(shared.read().unwrap().state()).unwrap(),
            before
        );
    }

    #[test]
    fn tata_gui_refuses_late_results_after_sequence_change() {
        let _lock = crate::tf_motifs::test_registry_lock()
            .lock()
            .unwrap_or_else(|poisoned| poisoned.into_inner());
        crate::tf_motifs::reload_builtin_for_test();
        let dna = DNAsequence::from_sequence("CCCCTATAAAACCCC").unwrap();
        let mut state = ProjectState::default();
        state.sequences.insert("toy".into(), dna.clone());
        let mut engine = GentleEngine::from_state(state);
        let result = engine
            .apply(Operation::ScreenTataBoxes {
                request: TataBoxScreenRequest {
                    seq_id: "toy".into(),
                    scan_without_tss: true,
                    ..Default::default()
                },
                path: None,
            })
            .unwrap();
        let revision = engine.structural_revision();
        engine
            .state_mut()
            .sequences
            .insert("other".into(), dna.clone());
        let shared = Arc::new(RwLock::new(engine));
        let mut area = MainAreaDna::new(dna, Some("toy".into()), Some(shared));
        let (tx, rx) = mpsc::channel();
        tx.send(Ok(result)).unwrap();
        area.tata_ui.task = Some(TataTask {
            receiver: Arc::new(Mutex::new(rx)),
            started: Instant::now(),
            baseline_revision: revision,
            mutating: false,
        });
        area.poll_tata_task(&egui::Context::default());
        assert!(area.tata_ui.report.is_none());
        assert_eq!(area.tata_ui.status, MainAreaDna::tr("tata.stale"));
    }
}
