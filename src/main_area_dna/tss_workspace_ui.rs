//! DNA-viewer preview/approval form for the engine-owned TSS collection workflow.

use super::*;
use gentle_protocol::tss_workspace::*;

#[derive(Debug, Clone)]
struct TssTask {
    receiver: Arc<Mutex<Receiver<Result<OpResult, EngineError>>>>,
    revision: u64,
    mutating: bool,
    started: Instant,
}

#[derive(Debug, Clone)]
pub(super) struct TssWorkspace {
    open: bool,
    gene: String,
    collection: String,
    upstream: usize,
    downstream: usize,
    selected: BTreeSet<String>,
    report: Option<Arc<TssInventoryReport>>,
    task: Option<TssTask>,
    materialized_collection: Option<String>,
    status: String,
}

impl Default for TssWorkspace {
    fn default() -> Self {
        Self {
            open: false,
            gene: String::new(),
            collection: "tss_windows".into(),
            upstream: 500,
            downstream: 200,
            selected: BTreeSet::new(),
            report: None,
            task: None,
            materialized_collection: None,
            status: String::new(),
        }
    }
}

impl MainAreaDna {
    fn tss_inventory_request(&self) -> TssInventoryRequest {
        TssInventoryRequest {
            seq_id: self.seq_id.clone().unwrap_or_default(),
            gene_query: self.tss_inventory_ui.gene.trim().into(),
            collection_id: self.tss_inventory_ui.collection.trim().into(),
            upstream_bp: self.tss_inventory_ui.upstream,
            downstream_bp: self.tss_inventory_ui.downstream,
        }
    }

    fn start_tss_operation(&mut self, operation: Operation, mutating: bool) {
        if self.tss_inventory_ui.task.is_some() {
            return;
        }
        let Some(engine) = self.engine.clone() else {
            self.tss_inventory_ui.status = "No project engine".into();
            return;
        };
        let revision = match engine.try_read() {
            Ok(e) => e.structural_revision(),
            Err(_) => {
                self.tss_inventory_ui.status = "Project busy; retry shortly".into();
                return;
            }
        };
        let (tx, rx) = mpsc::channel();
        self.tss_inventory_ui.task = Some(TssTask {
            receiver: Arc::new(Mutex::new(rx)),
            revision,
            mutating,
            started: Instant::now(),
        });
        self.tss_inventory_ui.status = "TSS operation running in the background".into();
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

    pub(super) fn poll_tss_task(&mut self, ctx: &egui::Context) {
        let Some(task) = self.tss_inventory_ui.task.clone() else {
            return;
        };
        let outcome = match task.receiver.lock() {
            Ok(rx) => match rx.try_recv() {
                Ok(result) => Some(result),
                Err(TryRecvError::Empty) => None,
                Err(TryRecvError::Disconnected) => {
                    Some(Err(EngineError::internal("TSS worker disconnected")))
                }
            },
            Err(_) => Some(Err(EngineError::internal("TSS worker channel unavailable"))),
        };
        let Some(outcome) = outcome else {
            ctx.request_repaint_after(Duration::from_millis(100));
            return;
        };
        self.tss_inventory_ui.task = None;
        match outcome {
            Ok(result) if task.mutating => {
                self.tss_inventory_ui.materialized_collection = result
                    .tss_collection
                    .as_ref()
                    .map(|r| r.collection_id.clone());
                self.handle_operation_success(result, task.started);
                self.tss_inventory_ui.status = "Collection is stored. Open its TSS windows below; no TFBS or occupancy analysis has been run.".into();
            }
            Ok(mut result) => {
                let revision = self
                    .engine
                    .as_ref()
                    .and_then(|e| e.try_read().ok().map(|g| g.structural_revision()));
                if revision != Some(task.revision) {
                    self.tss_inventory_ui.report = None;
                    self.tss_inventory_ui.status =
                        "Project changed during preview; inspect again".into();
                    return;
                }
                self.tss_inventory_ui.report = result.tss_inventory.take().map(|r| Arc::new(*r));
                self.tss_inventory_ui.selected.clear();
                self.tss_inventory_ui.status =
                    "Preview ready. Select the starts to materialize and approve below.".into();
            }
            Err(error) => self.tss_inventory_ui.status = error.to_string(),
        }
        ctx.request_repaint();
    }

    pub(super) fn render_tss_workspace(&mut self, ctx: &egui::Context) {
        if !self.tss_inventory_ui.open {
            return;
        }
        let mut open = true;
        let spec = crate::egui_compat::HostedWindowSpec::new(
            format!(
                "Transcript starts: {}",
                self.seq_id.as_deref().unwrap_or("")
            ),
            egui::Id::new(("tss_inventory", self.panel_scope_key())),
            Vec2::new(1000.0, 600.0),
            Vec2::new(620.0, 360.0),
        );
        crate::egui_compat::show_hosted_window(ctx, &spec, &mut open, |ui| {
            ui.label("Inspect exact annotated transcript starts on this loaded locus. Shared starts become one window; distinct starts and strands stay separate.");
            ui.small("Requires an anchored sequence with mRNA/transcript annotations. If flanks are missing, extend the locus first. These are annotation-derived candidates, not measured initiation sites.");
            let idle = self.tss_inventory_ui.task.is_none();
            ui.add_enabled_ui(idle, |ui| {
                ui.horizontal_wrapped(|ui| {
                    ui.label("Gene");
                    ui.text_edit_singleline(&mut self.tss_inventory_ui.gene);
                    ui.label("Collection ID");
                    ui.text_edit_singleline(&mut self.tss_inventory_ui.collection);
                });
                ui.horizontal_wrapped(|ui| {
                    ui.label("Upstream bp");
                    ui.add(
                        egui::DragValue::new(&mut self.tss_inventory_ui.upstream)
                            .range(0..=1_000_000),
                    );
                    ui.label("Downstream bp");
                    ui.add(
                        egui::DragValue::new(&mut self.tss_inventory_ui.downstream)
                            .range(0..=1_000_000),
                    );
                    if ui.button("Inspect starts (no changes)").clicked() {
                        self.start_tss_operation(
                            Operation::InspectTssInventory {
                                request: self.tss_inventory_request(),
                            },
                            false,
                        );
                    }
                });
            });
            if !idle {
                ui.spinner();
                ctx.request_repaint_after(Duration::from_millis(100));
            }
            ui.label(&self.tss_inventory_ui.status);
            if let Some(id) = self.tss_inventory_ui.materialized_collection.clone() {
                if ui
                    .add_enabled(
                        idle,
                        egui::Button::new(format!("Open TSS collection '{id}' (up to 32 windows)")),
                    )
                    .clicked()
                {
                    let queued = self
                        .engine
                        .as_ref()
                        .is_some_and(|e| crate::app::tss_collection_ui::request_open(e, &id));
                    self.tss_inventory_ui.status = if queued {
                        "Window request queued; validation and loading continue in the background"
                    } else {
                        "Another window request is pending; retry shortly"
                    }
                    .into();
                }
            }
            let Some(report) = self.tss_inventory_ui.report.clone() else {
                return;
            };
            for warning in &report.warnings {
                ui.small(warning);
            }
            let current = self.tss_inventory_request() == report.request;
            if !current {
                ui.colored_label(
                    egui::Color32::DARK_RED,
                    "Parameters changed: inspect again before approval",
                );
            }
            ui.add_enabled_ui(idle && current, |ui| {
                ui.horizontal(|ui| {
                    if ui.button("Select available").clicked() {
                        self.tss_inventory_ui.selected = report
                            .rows
                            .iter()
                            .filter(|r| r.availability == TssWindowAvailability::Available)
                            .map(|r| r.tss_id.clone())
                            .collect();
                    }
                    if ui.button("Clear selection").clicked() {
                        self.tss_inventory_ui.selected.clear();
                    }
                    if ui
                        .add_enabled(
                            !self.tss_inventory_ui.selected.is_empty(),
                            egui::Button::new("Approve and create selected windows"),
                        )
                        .clicked()
                    {
                        self.start_tss_operation(
                            Operation::MaterializeTssWindows {
                                request: TssMaterializeRequest {
                                    inventory: report.request.clone(),
                                    expected_approval_sha256: report.approval_sha256.clone(),
                                    selected_tss_ids: self
                                        .tss_inventory_ui
                                        .selected
                                        .iter()
                                        .cloned()
                                        .collect(),
                                },
                            },
                            true,
                        );
                    }
                });
                egui::ScrollArea::vertical()
                    .id_salt("tss_inventory_rows")
                    .show(ui, |ui| {
                        for row in &report.rows {
                            let mut selected = self.tss_inventory_ui.selected.contains(&row.tss_id);
                            let label = format!(
                                "{}:{} ({:?}) | {} transcript records | {}",
                                row.genomic_tss.reference.contig_name,
                                row.genomic_tss.start_0based + 1,
                                row.genomic_tss.strand,
                                row.transcript_feature_ids.len(),
                                row.explanation
                            );
                            if ui
                                .add_enabled(
                                    row.availability == TssWindowAvailability::Available,
                                    egui::Checkbox::new(&mut selected, label),
                                )
                                .on_hover_text(row.transcript_ids.join(", "))
                                .changed()
                            {
                                if selected {
                                    self.tss_inventory_ui.selected.insert(row.tss_id.clone());
                                } else {
                                    self.tss_inventory_ui.selected.remove(&row.tss_id);
                                }
                            }
                        }
                    });
            });
        });
        self.tss_inventory_ui.open = open;
    }

    pub(crate) fn open_tss_inventory(&mut self) {
        self.tss_inventory_ui.open = true;
    }

    #[cfg(test)]
    pub(crate) fn tss_inventory_workspace_open(&self) -> bool {
        self.tss_inventory_ui.open
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn tss_workspace_preview_is_background_read_only_and_form_bound() {
        let engine = crate::engine::synthetic_tss_engine(false);
        let request = crate::engine::synthetic_tss_approval(&engine).inventory;
        let dna = engine.state().sequences["locus"].clone();
        let expected = engine.inspect_tss_inventory(&request).unwrap();
        let shared = Arc::new(RwLock::new(engine));
        let mut area = MainAreaDna::new(dna, Some("locus".into()), Some(shared.clone()));
        // Viewer construction prepares its reasoning overlay before our read-only action.
        let before = serde_json::to_value(shared.read().unwrap().state()).unwrap();
        area.tss_inventory_ui.gene = request.gene_query.clone();
        area.tss_inventory_ui.collection = request.collection_id.clone();
        area.tss_inventory_ui.upstream = request.upstream_bp;
        area.tss_inventory_ui.downstream = request.downstream_bp;
        area.open_tss_inventory();
        area.start_tss_operation(Operation::InspectTssInventory { request }, false);
        let ctx = egui::Context::default();
        let deadline = Instant::now() + Duration::from_secs(10);
        while area.tss_inventory_ui.task.is_some() && Instant::now() < deadline {
            area.poll_tss_task(&ctx);
            std::thread::sleep(Duration::from_millis(5));
        }
        assert_eq!(area.tss_inventory_ui.report.as_deref(), Some(&expected));
        assert_eq!(
            serde_json::to_value(shared.read().unwrap().state()).unwrap(),
            before
        );
        assert_eq!(area.tss_inventory_request(), expected.request);
        area.tss_inventory_ui.upstream += 1;
        assert_ne!(area.tss_inventory_request(), expected.request);
        for _ in 0..2 {
            let mut output = ctx.run_ui(egui::RawInput::default(), |ui| {
                area.render_tss_workspace(ui.ctx())
            });
            output.textures_delta.clear();
        }
    }
}
