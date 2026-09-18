//! DNA-viewer preview/approval form for the engine-owned TSS collection workflow.

use super::*;
use crate::tutorial_gui_semantics::*;
use gentle_protocol::tss_workspace::*;

fn tss_control(
    response: egui::Response,
    id: &'static str,
    seq_id: &str,
    outcome: Option<&str>,
) -> egui::Response {
    #[cfg(feature = "gui-test-support")]
    {
        use crate::gui_test_support::{
            GuiTestWidgetKind, pseudonymous_subject_scope, register_response_with_outcome,
        };
        let spec = tutorial_gui_control(id).expect("TSS control must be catalogued");
        let kind = if spec.text_policy.is_some() {
            GuiTestWidgetKind::TextInput
        } else if spec.authority == TutorialGuiControlAuthority::Observe {
            GuiTestWidgetKind::Status
        } else {
            GuiTestWidgetKind::Button
        };
        register_response_with_outcome(
            &response,
            id,
            WINDOW_TSS_WORKSPACE,
            Some(&pseudonymous_subject_scope(&[seq_id])),
            kind,
            false,
            outcome,
        );
    }
    #[cfg(not(feature = "gui-test-support"))]
    let _ = (id, seq_id, outcome);
    response
}

#[derive(Debug, Clone)]
struct TssTask {
    receiver: Arc<Mutex<Receiver<Result<OpResult, EngineError>>>>,
    revision: u64,
    owner: u64,
    mutating: bool,
    started: Instant,
    forgotten_collection: Option<String>,
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
    inspected_collection: Option<Arc<TssCollectionReport>>,
    collection_list: Option<Arc<TssCollectionListReport>>,
    inspected_at: Option<(u64, u64)>,
    operation_failed: bool,
    forget_confirmation: Option<String>,
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
            inspected_collection: None,
            collection_list: None,
            inspected_at: None,
            operation_failed: false,
            forget_confirmation: None,
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
        let (owner, revision) = match engine.try_read() {
            Ok(e) => (e.instance_id(), e.structural_revision()),
            Err(_) => {
                self.tss_inventory_ui.status = "Project busy; retry shortly".into();
                return;
            }
        };
        let (tx, rx) = mpsc::channel();
        self.tss_inventory_ui.task = Some(TssTask {
            receiver: Arc::new(Mutex::new(rx)),
            revision,
            owner,
            mutating,
            started: Instant::now(),
            forgotten_collection: match &operation {
                Operation::ForgetTssCollection { collection_id } => Some(collection_id.clone()),
                _ => None,
            },
        });
        self.tss_inventory_ui.inspected_collection = None;
        self.tss_inventory_ui.inspected_at = None;
        self.tss_inventory_ui.operation_failed = false;
        self.tss_inventory_ui.forget_confirmation = None;
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
                self.tss_inventory_ui.collection_list = None;
                self.tss_inventory_ui.materialized_collection = result
                    .tss_collection
                    .as_ref()
                    .map(|r| r.collection_id.clone());
                self.handle_operation_success(result, task.started);
                self.tss_inventory_ui.status = if let Some(id) = task.forgotten_collection {
                    format!(
                        "Forgot registry entry '{id}'. Sequences, windows and lineage are retained. Use a new collection ID to avoid overwriting retained sequences."
                    )
                } else {
                    "Collection is stored. Open its TSS windows below; no TFBS or occupancy analysis has been run.".into()
                };
            }
            Ok(mut result) => {
                let identity = self.engine.as_ref().and_then(|e| {
                    e.try_read()
                        .ok()
                        .map(|g| (g.instance_id(), g.structural_revision()))
                });
                if identity != Some((task.owner, task.revision)) {
                    self.tss_inventory_ui.operation_failed = true;
                    self.tss_inventory_ui.report = None;
                    self.tss_inventory_ui.status =
                        "Project changed during preview; inspect again".into();
                    return;
                }
                if let Some(list) = result.tss_collection_list.take() {
                    self.tss_inventory_ui.collection_list = Some(Arc::new(*list));
                    self.tss_inventory_ui.status = "Registry refreshed; member validation not checked. Select an entry, then inspect it explicitly.".into();
                    ctx.request_repaint();
                    return;
                }
                if let Some(collection) = result.tss_collection.take() {
                    self.tss_inventory_ui.inspected_at = identity;
                    self.tss_inventory_ui.materialized_collection =
                        Some(collection.collection_id.clone());
                    self.tss_inventory_ui.inspected_collection = Some(Arc::new(*collection));
                    self.tss_inventory_ui.status = "Collection and all member snapshots validated. This does not establish TFBS, occupancy or primer specificity.".into();
                    ctx.request_repaint();
                    return;
                }
                self.tss_inventory_ui.report = result.tss_inventory.take().map(|r| Arc::new(*r));
                self.tss_inventory_ui.selected.clear();
                self.tss_inventory_ui.status =
                    "Preview ready. Select the starts to materialize and approve below.".into();
            }
            Err(error) => {
                self.tss_inventory_ui.operation_failed = true;
                self.tss_inventory_ui.materialized_collection = None;
                self.tss_inventory_ui.status = error.to_string();
            }
        }
        ctx.request_repaint();
    }

    fn request_tss_forget(&mut self) {
        let id = self.tss_inventory_ui.collection.trim();
        if self.tss_inventory_ui.task.is_none() && !id.is_empty() {
            self.tss_inventory_ui.forget_confirmation = Some(id.into());
        }
    }

    fn select_tss_collection(&mut self, id: String) {
        self.tss_inventory_ui.operation_failed = false;
        self.tss_inventory_ui.collection = id;
        self.tss_inventory_ui.inspected_collection = None;
        self.tss_inventory_ui.inspected_at = None;
        self.tss_inventory_ui.materialized_collection = None;
        self.tss_inventory_ui.forget_confirmation = None;
        self.tss_inventory_ui.status =
            "Selected registry entry; member validation not checked".into();
    }

    fn invalidate_tss_inspection(&mut self) {
        let identity = self.engine.as_ref().and_then(|e| {
            e.try_read()
                .ok()
                .map(|g| (g.instance_id(), g.structural_revision()))
        });
        if self.tss_inventory_ui.inspected_at.is_some()
            && self.tss_inventory_ui.inspected_at != identity
        {
            self.tss_inventory_ui.inspected_collection = None;
            self.tss_inventory_ui.inspected_at = None;
            self.tss_inventory_ui.materialized_collection = None;
            self.tss_inventory_ui.status =
                "Project changed since validation; inspect the collection again".into();
        }
    }

    fn confirm_tss_forget(&mut self) {
        if let Some(id) = self.tss_inventory_ui.forget_confirmation.take()
            && id == self.tss_inventory_ui.collection.trim()
            && self.tss_inventory_ui.task.is_none()
        {
            self.start_tss_operation(Operation::ForgetTssCollection { collection_id: id }, true);
        }
    }

    pub(super) fn render_tss_workspace(&mut self, ctx: &egui::Context) {
        if !self.tss_inventory_ui.open {
            return;
        }
        self.invalidate_tss_inspection();
        let scope_seq = self.seq_id.clone().unwrap_or_default();
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
            #[cfg(feature = "gui-test-support")]
            crate::gui_test_support::register_rect(
                ui.ctx().clone(),
                WINDOW_TSS_WORKSPACE,
                crate::tutorial_gui_semantics::WINDOW_DNA_VIEWER,
                Some(&crate::gui_test_support::pseudonymous_subject_scope(&[
                    &scope_seq,
                ])),
                crate::gui_test_support::GuiTestWidgetKind::Status,
                ui.max_rect(),
                true,
                true,
                true,
                Some("ready"),
            );
            egui::ScrollArea::vertical().id_salt("tss_workspace_scroll").show(ui, |ui| {
            ui.label("Inspect exact annotated transcript starts on this loaded locus. Shared starts become one window; distinct starts and strands stay separate.");
            ui.small("Requires an anchored sequence with mRNA/transcript annotations. If flanks are missing, extend the locus first. These are annotation-derived candidates, not measured initiation sites.");
            let idle = self.tss_inventory_ui.task.is_none();
            ui.add_enabled_ui(idle, |ui| {
                ui.horizontal_wrapped(|ui| {
                    ui.label("Gene");
                    tss_control(ui.text_edit_singleline(&mut self.tss_inventory_ui.gene), TSS_GENE, &scope_seq, None);
                    ui.label("Collection ID");
                    if tss_control(ui.text_edit_singleline(&mut self.tss_inventory_ui.collection), TSS_COLLECTION, &scope_seq, None).changed() {
                        self.select_tss_collection(self.tss_inventory_ui.collection.clone());
                    }
                });
                ui.horizontal_wrapped(|ui| {
                    if tss_control(ui.button("Refresh collections"), TSS_REFRESH, &scope_seq, None).clicked() {
                        self.start_tss_operation(Operation::ListTssCollections {}, false);
                    }
                    let has_id = !self.tss_inventory_ui.collection.trim().is_empty();
                    if tss_control(ui.add_enabled(has_id, egui::Button::new("Inspect stored collection")), TSS_INSPECT, &scope_seq, None).clicked() {
                        self.start_tss_operation(Operation::GetTssCollection {
                            collection_id: self.tss_inventory_ui.collection.trim().into(),
                        }, false);
                    }
                    if tss_control(ui.add_enabled(has_id, egui::Button::new("Forget registry entry...")), TSS_FORGET, &scope_seq, None).clicked() {
                        self.request_tss_forget();
                    }
                });
                if let Some(list) = self.tss_inventory_ui.collection_list.clone() {
                    if list.collections.is_empty() { ui.label("No stored TSS collections"); }
                    egui::ScrollArea::vertical().id_salt("tss_collection_browser").max_height(110.0).show(ui, |ui| {
                        for entry in &list.collections {
                            let label = format!("{} | gene {} | locus {} | {} windows | {:?} / not checked",
                                entry.collection_id,
                                entry.gene_query.as_deref().unwrap_or("unavailable"),
                                entry.source_seq_id.as_deref().unwrap_or("unavailable"),
                                entry.window_count.map(|n| n.to_string()).unwrap_or_else(|| "unavailable".into()),
                                entry.record_status);
                            let selected = self.tss_inventory_ui.collection == entry.collection_id;
                            let row = ui.selectable_label(selected, label);
                            #[cfg(feature = "gui-test-support")]
                            crate::gui_test_support::register_response(&row, TSS_COLLECTION_ROW, WINDOW_TSS_WORKSPACE,
                                Some(&crate::gui_test_support::pseudonymous_subject_scope(&[&scope_seq, &entry.collection_id])), crate::gui_test_support::GuiTestWidgetKind::Row, selected);
                            if row.clicked() { self.select_tss_collection(entry.collection_id.clone()); }
                            if let Some(note) = &entry.diagnostic { row.on_hover_text(note); }
                        }
                    });
                }
                if let Some(id) = self.tss_inventory_ui.forget_confirmation.clone() {
                    ui.label(format!("Forget '{id}'? Only registry metadata is removed; sequences and lineage remain. Re-derivation normally needs a new collection ID."));
                    ui.horizontal(|ui| {
                        if tss_control(ui.button("Confirm forget registry entry"), TSS_CONFIRM_FORGET, &scope_seq, None).clicked() { self.confirm_tss_forget(); }
                        if ui.button("Cancel").clicked() { self.tss_inventory_ui.forget_confirmation = None; }
                    });
                }
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
                    if tss_control(ui.button("Inspect starts (no changes)"), TSS_PREVIEW, &scope_seq, None).clicked() {
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
            let outcome = if !idle { "running" }
                else if self.tss_inventory_ui.operation_failed { "failed" }
                else if self.tss_inventory_ui.inspected_collection.is_some() { "validated" }
                else if self.tss_inventory_ui.collection_list.is_some() { "not_checked" }
                else { "idle" };
            tss_control(ui.label(&self.tss_inventory_ui.status), TSS_STATUS, &scope_seq, Some(outcome));
            if let Some(collection) = self.tss_inventory_ui.inspected_collection.clone() {
                ui.label(format!(
                    "{}: {} validated windows; gene {}",
                    collection.collection_id,
                    collection.members.len(),
                    collection.inventory.request.gene_query
                ));
                ui.monospace(&collection.collection_membership_fingerprint_sha256);
                if ui.button("Copy collection JSON").clicked()
                    && let Ok(json) = serde_json::to_string_pretty(collection.as_ref())
                {
                    ui.ctx().copy_text(json);
                }
                egui::ScrollArea::vertical()
                    .id_salt("tss_collection_members")
                    .max_height(160.0)
                    .show(ui, |ui| {
                        for member in &collection.members {
                            ui.label(format!(
                                "{}: {}:{} ({:?}) | {}",
                                member.tss.output_seq_id,
                                member.tss.genomic_tss.reference.contig_name,
                                member.tss.genomic_tss.start_0based + 1,
                                member.tss.genomic_tss.strand,
                                member.tss.transcript_ids.join(", ")
                            ));
                        }
                    });
            }
            if let Some(id) = self.tss_inventory_ui.materialized_collection.clone() {
                if tss_control(ui
                    .add_enabled(
                        idle,
                        egui::Button::new(format!("Open TSS collection '{id}' (up to 32 windows)")),
                    ), TSS_OPEN_WINDOWS, &scope_seq, None)
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
                    if tss_control(ui.button("Select available"), TSS_SELECT_AVAILABLE, &scope_seq, None).clicked() {
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
                    if tss_control(ui
                        .add_enabled(
                            !self.tss_inventory_ui.selected.is_empty(),
                            egui::Button::new("Approve and create selected windows"),
                        ), TSS_MATERIALIZE, &scope_seq, None)
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
                        for excluded in &report.excluded_transcripts {
                            ui.add_enabled(
                                false,
                                egui::Label::new(format!(
                                    "{} (feature {}): {}",
                                    excluded.transcript_id,
                                    excluded.feature_id,
                                    excluded.explanation
                                )),
                            );
                        }
                        for unassigned in &report.unassigned_transcripts {
                            ui.add_enabled(false, egui::Label::new(format!("Unassigned in locus: {} (feature {}); gene linkage unavailable", unassigned.transcript_id, unassigned.feature_id)));
                        }
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

    fn finish_task(area: &mut MainAreaDna, ctx: &egui::Context) {
        let deadline = Instant::now() + Duration::from_secs(10);
        while area.tss_inventory_ui.task.is_some() && Instant::now() < deadline {
            area.poll_tss_task(ctx);
            std::thread::sleep(Duration::from_millis(5));
        }
        assert!(area.tss_inventory_ui.task.is_none());
    }

    #[test]
    fn tss_browser_lists_without_validation_and_clears_old_inspection() {
        let mut engine = crate::engine::synthetic_tss_engine(false);
        let request = crate::engine::synthetic_tss_approval(&engine);
        engine
            .apply(Operation::MaterializeTssWindows { request })
            .unwrap();
        let dna = engine.state().sequences["locus"].clone();
        let shared = Arc::new(RwLock::new(engine));
        let mut area = MainAreaDna::new(dna, Some("locus".into()), Some(shared.clone()));
        let ctx = egui::Context::default();
        let before = serde_json::to_value(shared.read().unwrap().state()).unwrap();
        area.start_tss_operation(Operation::ListTssCollections {}, false);
        finish_task(&mut area, &ctx);
        assert_eq!(
            area.tss_inventory_ui
                .collection_list
                .as_ref()
                .unwrap()
                .collections
                .len(),
            1
        );
        assert!(area.tss_inventory_ui.inspected_collection.is_none());
        area.select_tss_collection("toy_tss".into());
        assert!(area.tss_inventory_ui.task.is_none());
        assert!(area.tss_inventory_ui.materialized_collection.is_none());
        area.start_tss_operation(
            Operation::GetTssCollection {
                collection_id: "toy_tss".into(),
            },
            false,
        );
        finish_task(&mut area, &ctx);
        assert!(area.tss_inventory_ui.inspected_collection.is_some());
        assert_eq!(
            serde_json::to_value(shared.read().unwrap().state()).unwrap(),
            before
        );
        // A replacement project can have the same revision: owner identity also matters.
        let state = shared.read().unwrap().state().clone();
        *shared.write().unwrap() = GentleEngine::from_state(state);
        area.invalidate_tss_inspection();
        assert!(area.tss_inventory_ui.inspected_collection.is_none());
        assert!(
            area.tss_inventory_ui
                .status
                .contains("inspect the collection again")
        );
    }

    #[test]
    fn tss_workspace_inspects_stale_failure_and_explicitly_forgets_without_deleting() {
        let mut engine = crate::engine::synthetic_tss_engine(false);
        let request = crate::engine::synthetic_tss_approval(&engine);
        let collection = engine
            .apply(Operation::MaterializeTssWindows { request })
            .unwrap()
            .tss_collection
            .unwrap();
        let member_id = collection.members[0].tss.output_seq_id.clone();
        let dna = engine.state().sequences["locus"].clone();
        let shared = Arc::new(RwLock::new(engine));
        let mut area = MainAreaDna::new(dna, Some("locus".into()), Some(shared.clone()));
        let ctx = egui::Context::default();
        area.open_tss_inventory();
        area.tss_inventory_ui.collection = "toy_tss".into();
        let before = serde_json::to_value(shared.read().unwrap().state()).unwrap();
        area.start_tss_operation(
            Operation::GetTssCollection {
                collection_id: "toy_tss".into(),
            },
            false,
        );
        finish_task(&mut area, &ctx);
        assert_eq!(
            area.tss_inventory_ui
                .inspected_collection
                .as_ref()
                .unwrap()
                .members
                .len(),
            2
        );
        assert_eq!(
            serde_json::to_value(shared.read().unwrap().state()).unwrap(),
            before
        );
        shared
            .write()
            .unwrap()
            .state_mut()
            .sequences
            .get_mut(&member_id)
            .unwrap()
            .features_mut()
            .clear();
        area.start_tss_operation(
            Operation::GetTssCollection {
                collection_id: "toy_tss".into(),
            },
            false,
        );
        finish_task(&mut area, &ctx);
        assert!(area.tss_inventory_ui.inspected_collection.is_none());
        assert!(area.tss_inventory_ui.status.contains("edited"));
        let before = serde_json::to_value(shared.read().unwrap().state()).unwrap();
        area.request_tss_forget();
        assert!(area.tss_inventory_ui.task.is_none());
        assert_eq!(
            serde_json::to_value(shared.read().unwrap().state()).unwrap(),
            before
        );
        area.tss_inventory_ui.collection = "different_id".into();
        area.confirm_tss_forget();
        assert!(area.tss_inventory_ui.task.is_none());
        assert_eq!(
            serde_json::to_value(shared.read().unwrap().state()).unwrap(),
            before
        );
        area.tss_inventory_ui.collection = "toy_tss".into();
        area.request_tss_forget();
        for _ in 0..2 {
            let mut output = ctx.run_ui(egui::RawInput::default(), |ui| {
                area.render_tss_workspace(ui.ctx())
            });
            output.textures_delta.clear();
        }
        area.confirm_tss_forget();
        finish_task(&mut area, &ctx);
        assert!(
            area.tss_inventory_ui
                .status
                .contains("Forgot registry entry")
        );
        let mut expected = before;
        expected["metadata"]["tss_collections_v1"]
            .as_object_mut()
            .unwrap()
            .remove("toy_tss");
        assert_eq!(
            serde_json::to_value(shared.read().unwrap().state()).unwrap(),
            expected
        );
        shared.write().unwrap().undo_last_operation().unwrap();
        assert!(
            shared.read().unwrap().state().metadata["tss_collections_v1"]
                .get("toy_tss")
                .is_some()
        );
    }

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
