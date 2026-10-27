//! Read-only presentation of the shared protein-expression handoff report.
//!
//! Request fields are session-local. Workers discard their history-free engine
//! snapshot; neither this window nor its copy-only actions apply a service job.

use super::*;
use crate::engine::{PlanningProfileScope, ProteinExpressionHandoffReport};
use crate::engine_shell::execute_shell_command;
use std::sync::mpsc::{self, Receiver, TryRecvError};

#[cfg(feature = "gui-test-support")]
const CAPTURE_WINDOW: &str = "window.protein_expression_handoff";

#[cfg(feature = "gui-test-support")]
fn capture_control(
    response: &egui::Response,
    id: &'static str,
    kind: crate::gui_test_support::GuiTestWidgetKind,
    outcome: Option<&str>,
) {
    crate::gui_test_support::register_response_with_outcome(
        response,
        id,
        CAPTURE_WINDOW,
        None,
        kind,
        false,
        outcome,
    );
}

#[derive(Clone, Debug, Default, PartialEq, Eq)]
struct Inputs {
    seq_id: String,
    requirements: String,
    provider: String,
}

impl Inputs {
    fn command(&self) -> Result<ShellCommand, String> {
        let requirements = self.requirements.trim();
        // The GUI accepts pasted JSON only. File loading remains an explicit
        // shared-shell action rather than an implicit path in an editable field.
        if !requirements.is_empty() {
            serde_json::from_str::<serde_json::Value>(requirements)
                .map_err(|err| format!("Reviewed requirements must be inline JSON: {err}"))?;
        }
        let optional = |value: &str| {
            let value = value.trim();
            (!value.is_empty()).then(|| value.to_string())
        };
        Ok(ShellCommand::PlanningProteinExpressionHandoff {
            seq_id: optional(&self.seq_id),
            objective_json: None,
            requirements_json: optional(requirements),
            provider: optional(&self.provider),
            profile_scope: PlanningProfileScope::Effective,
            output_format: "text".into(),
        })
    }
}

#[derive(Clone, Debug, PartialEq, Eq)]
struct Binding {
    instance: u64,
    revision: u64,
    inputs: Inputs,
}

impl Binding {
    fn new(engine: &GentleEngine, inputs: Inputs) -> Self {
        Self {
            instance: engine.instance_id(),
            revision: engine.mutation_revision(),
            inputs,
        }
    }

    fn matches(&self, engine: &GentleEngine, inputs: &Inputs) -> bool {
        self.instance == engine.instance_id()
            && self.revision == engine.mutation_revision()
            && &self.inputs == inputs
    }
}

struct Pending {
    binding: Binding,
    receiver: Receiver<Result<ProteinExpressionHandoffReport, String>>,
}

#[derive(Default)]
pub(super) struct ProteinExpressionUiState {
    pub(super) open: bool,
    inputs: Inputs,
    pending: Option<Pending>,
    report: Option<ProteinExpressionHandoffReport>,
    binding: Option<Binding>,
    error: Option<String>,
}

fn inspect_snapshot(
    shared: &Arc<RwLock<GentleEngine>>,
    binding: &Binding,
) -> Result<ProteinExpressionHandoffReport, String> {
    let command = binding.inputs.command()?;
    let mut snapshot = {
        let engine = shared.read().map_err(|_| "Engine lock poisoned")?;
        if !binding.matches(&engine, &binding.inputs) {
            return Err("Project changed before inspection; inspect again".into());
        }
        engine.clone_without_history()
    };
    let result = execute_shell_command(&mut snapshot, &command)?;
    if result.state_changed {
        return Err("Read-only protein-expression inspection unexpectedly changed state".into());
    }
    let report: ProteinExpressionHandoffReport = serde_json::from_value(result.output)
        .map_err(|err| format!("Invalid shared protein-expression report: {err}"))?;
    if report.schema != "gentle.protein_expression_handoff.v2" {
        return Err(format!(
            "Unsupported protein-expression report: {}",
            report.schema
        ));
    }
    Ok(report)
}

fn count(value: Option<usize>) -> String {
    value
        .map(|value| value.to_string())
        .unwrap_or_else(|| "unavailable".into())
}

fn details(ui: &mut egui::Ui, label: &str, value: &impl serde::Serialize) {
    ui.collapsing(label, |ui| {
        if let Ok(text) = serde_json::to_string_pretty(value) {
            ui.add(
                egui::Label::new(egui::RichText::new(text).monospace())
                    .selectable(true)
                    .wrap(),
            );
        }
    });
}

impl GENtleApp {
    pub(super) fn protein_expression_inspection_pending(&self) -> bool {
        self.protein_expression_ui.pending.is_some()
    }

    pub(super) fn protein_expression_viewport_id() -> ViewportId {
        ViewportId::from_hash_of("GENtle Protein Expression Handoff Viewport")
    }

    pub(super) fn hosted_protein_expression_window_id() -> egui::Id {
        egui::Id::new((
            "hosted_protein_expression_window",
            Self::protein_expression_viewport_id(),
        ))
    }

    pub(super) fn open_protein_expression_dialog(&mut self) {
        let was_open = self.protein_expression_ui.open;
        self.protein_expression_ui.open = true;
        self.mark_window_open_or_focus(Self::protein_expression_viewport_id(), was_open);
    }

    pub(super) fn close_protein_expression_dialog(&mut self) -> bool {
        let was_open = self.protein_expression_ui.open;
        self.protein_expression_ui.open = false;
        // Closing abandons the receiver, not a live project transaction.
        self.protein_expression_ui.pending = None;
        was_open
    }

    fn start_protein_expression_report(&mut self, ctx: &egui::Context) {
        let state = &mut self.protein_expression_ui;
        if state.pending.is_some() {
            return;
        }
        state.report = None;
        state.binding = None;
        state.error = None;
        let binding = match self.engine.try_read() {
            Ok(engine) => Binding::new(&engine, state.inputs.clone()),
            Err(_) => {
                state.error = Some("Project is busy; inspect again when available".into());
                return;
            }
        };
        let shared = self.engine.clone();
        let worker_binding = binding.clone();
        let repaint = ctx.clone();
        let (sender, receiver) = mpsc::channel();
        match std::thread::Builder::new()
            .name("protein-expression-inspector".into())
            .stack_size(16 * 1024 * 1024)
            .spawn(move || {
                let _ = sender.send(inspect_snapshot(&shared, &worker_binding));
                repaint.request_repaint();
                repaint.request_repaint_of(ViewportId::ROOT);
            }) {
            Ok(_) => state.pending = Some(Pending { binding, receiver }),
            Err(err) => state.error = Some(format!("Could not start inspection: {err}")),
        }
    }

    pub(super) fn poll_protein_expression_report(&mut self) {
        let state = &mut self.protein_expression_ui;
        let Ok(engine) = self.engine.try_read() else {
            return;
        };
        let stale = state
            .binding
            .as_ref()
            .is_some_and(|binding| !binding.matches(&engine, &state.inputs))
            || state
                .pending
                .as_ref()
                .is_some_and(|pending| !pending.binding.matches(&engine, &state.inputs));
        if stale {
            state.report = None;
            state.binding = None;
            state.pending = None;
            state.error =
                Some("Project or request changed; inspect again for a current report".into());
            return;
        }
        let Some(pending) = &state.pending else {
            return;
        };
        match pending.receiver.try_recv() {
            Ok(result) => {
                state.binding = Some(pending.binding.clone());
                state.pending = None;
                match result {
                    Ok(report) => state.report = Some(report),
                    Err(err) => state.error = Some(err),
                }
            }
            Err(TryRecvError::Empty) => {}
            Err(TryRecvError::Disconnected) => {
                state.pending = None;
                state.error = Some("Inspection worker stopped; inspect again".into());
            }
        }
    }

    fn render_protein_expression_contents(&mut self, ui: &mut egui::Ui) {
        help_topics::button(ui, HelpTopic::ProteinExpression);
        self.poll_protein_expression_report();
        #[cfg(feature = "gui-test-support")]
        crate::gui_test_support::register_rect(
            ui.ctx().clone(),
            CAPTURE_WINDOW,
            CAPTURE_WINDOW,
            None,
            crate::gui_test_support::GuiTestWidgetKind::Window,
            ui.max_rect(),
            true,
            true,
            false,
            None,
        );
        if self.engine.try_read().is_err() {
            ui.label("Project is busy; waiting to validate the current report");
            ui.ctx()
                .request_repaint_after(std::time::Duration::from_millis(100));
            return;
        }
        egui::ScrollArea::vertical()
            .id_salt("protein_expression_inspector")
            .show(ui, |ui| {
                let state = &mut self.protein_expression_ui;
                let before = state.inputs.clone();
                ui.heading("Protein Expression Handoff");
                ui.horizontal_wrapped(|ui| {
                    ui.label("Product ID");
                    let _product = ui.add(
                        egui::TextEdit::singleline(&mut state.inputs.seq_id).desired_width(220.0),
                    );
                    #[cfg(feature = "gui-test-support")]
                    capture_control(
                        &_product,
                        "protein_expression.product",
                        crate::gui_test_support::GuiTestWidgetKind::TextInput,
                        None,
                    );
                    let sequences = self
                        .engine
                        .try_read()
                        .ok()
                        .map(|engine| {
                            let mut rows =
                                engine.state().sequences.keys().cloned().collect::<Vec<_>>();
                            rows.sort();
                            rows
                        })
                        .unwrap_or_default();
                    egui::ComboBox::from_id_salt("protein_expression_product")
                        .selected_text("Loaded products")
                        .show_ui(ui, |ui| {
                            ui.selectable_value(
                                &mut state.inputs.seq_id,
                                String::new(),
                                "No product",
                            );
                            for id in sequences {
                                ui.selectable_value(&mut state.inputs.seq_id, id.clone(), &id);
                            }
                        });
                });
                ui.horizontal_wrapped(|ui| {
                    ui.label("Provider (optional)");
                    let _provider = ui.add(
                        egui::TextEdit::singleline(&mut state.inputs.provider).desired_width(220.0),
                    );
                    #[cfg(feature = "gui-test-support")]
                    capture_control(
                        &_provider,
                        "protein_expression.provider",
                        crate::gui_test_support::GuiTestWidgetKind::TextInput,
                        None,
                    );
                });
                ui.collapsing("Reviewed requirements JSON (optional)", |ui| {
                    let _requirements = ui.add(
                        egui::TextEdit::multiline(&mut state.inputs.requirements)
                            .code_editor()
                            .desired_rows(8)
                            .desired_width(f32::INFINITY),
                    );
                    #[cfg(feature = "gui-test-support")]
                    capture_control(
                        &_requirements,
                        "protein_expression.requirements",
                        crate::gui_test_support::GuiTestWidgetKind::TextInput,
                        None,
                    );
                });
                if before != state.inputs {
                    state.report = None;
                    state.binding = None;
                    state.pending = None;
                    state.error = None;
                }
                let inspect = ui
                    .horizontal_wrapped(|ui| {
                        let inspect =
                            ui.add_enabled(state.pending.is_none(), egui::Button::new("Inspect"));
                        #[cfg(feature = "gui-test-support")]
                        capture_control(
                            &inspect,
                            "protein_expression.inspect",
                            crate::gui_test_support::GuiTestWidgetKind::Button,
                            None,
                        );
                        if state.pending.is_some() {
                            ui.spinner();
                            ui.label("Inspecting product context");
                        }
                        inspect.clicked()
                    })
                    .inner;
                if let Some(error) = &state.error {
                    ui.colored_label(egui::Color32::RED, error);
                }
                if let Some(report) = &state.report {
                    ui.separator();
                    let _status = ui.heading(&report.status);
                    #[cfg(feature = "gui-test-support")]
                    capture_control(
                        &_status,
                        "protein_expression.status",
                        crate::gui_test_support::GuiTestWidgetKind::Status,
                        Some(&report.status),
                    );
                    ui.label(&report.product_readiness.review_gate);
                    if let Some(name) = &report.product_definition.sequence_name {
                        ui.label(name);
                    }
                    ui.label(format!(
                        "Translation possible: {}",
                        report.product_readiness.translation_possible
                    ));
                    for blocker in &report.product_readiness.blockers {
                        ui.label(blocker);
                    }
                    details(
                        ui,
                        "Readiness and product definition",
                        &(&report.product_readiness, &report.product_definition),
                    );
                    egui::CollapsingHeader::new("Unanswered requirements")
                        .default_open(true)
                        .show(ui, |ui| {
                            if report.missing_questions.is_empty() {
                                ui.label("None reported");
                            }
                            for question in &report.missing_questions {
                                ui.strong(format!("{}: {}", question.question_id, question.prompt));
                                ui.label(&question.reason);
                            }
                        });
                    details(ui, "Reviewed requirements", &report.requirements);
                    ui.collapsing("Sequence and annotation evidence", |ui| {
                        if let Some(context) = &report.sequence_context {
                            ui.label(format!(
                                "Product ID: {}",
                                context.seq_id.as_deref().unwrap_or("unavailable")
                            ));
                            ui.label(format!(
                                "DNA length (bp): {}",
                                count(context.nucleotide_length)
                            ));
                            ui.label(format!(
                                "Protein length (aa): {}",
                                count(context.protein_length_aa)
                            ));
                            ui.label(format!("Annotations: {}", context.feature_count));
                            ui.label(format!("Ambiguous bases: {}", context.ambiguous_base_count));
                        } else {
                            ui.label("Sequence context unavailable");
                        }
                        details(
                            ui,
                            "Exact sequence/annotation fields",
                            &report.sequence_context,
                        );
                    });
                    ui.collapsing("CDS and translation evidence", |ui| {
                        let cds = &report.cds_assessment;
                        ui.label(format!("Context source: {}", cds.context_source));
                        ui.label(format!("CDS length (nt): {}", count(cds.nucleotide_length)));
                        ui.label(format!(
                            "Protein length (aa): {}",
                            count(cds.protein_length_aa)
                        ));
                        ui.label(format!(
                            "Start codon: {}",
                            cds.start_codon.as_deref().unwrap_or("unavailable")
                        ));
                        ui.label(format!(
                            "Stop codon: {}",
                            cds.stop_codon.as_deref().unwrap_or("unavailable")
                        ));
                        for warning in &cds.warnings {
                            ui.label(warning);
                        }
                        details(ui, "Exact CDS/translation fields", cds);
                    });
                    ui.collapsing("Tag and signal-peptide evidence", |ui| {
                        let tags = &report.tag_assessment;
                        ui.label(format!("Annotated tags: {}", tags.annotated_tag_count));
                        ui.label(format!(
                            "Annotated signal peptides: {}",
                            tags.annotated_signal_peptide_count
                        ));
                        ui.label(format!("Tag policy known: {}", tags.tag_policy_known));
                        for line in tags.missing_inputs.iter().chain(&tags.warnings) {
                            ui.label(line);
                        }
                        details(ui, "Exact tag/signal fields", tags);
                    });
                    ui.collapsing("Chassis review candidates", |ui| {
                        for row in &report.host_chassis_candidates {
                            ui.strong(format!("{}. {} [{}]", row.rank, row.label, row.status));
                            for line in row.rationale.iter().chain(&row.warnings) {
                                ui.label(line);
                            }
                        }
                    });
                    ui.collapsing("Vector-route review candidates", |ui| {
                        for row in &report.vector_route_candidates {
                            ui.strong(format!("{}. {} [{}]", row.rank, row.label, row.status));
                            for line in row.rationale.iter().chain(&row.missing_inputs) {
                                ui.label(line);
                            }
                        }
                    });
                    let binding_status = serde_json::to_value(report.provider_binding.status)
                        .ok()
                        .and_then(|value| value.as_str().map(str::to_string))
                        .unwrap_or_default();
                    let _provider_binding =
                        ui.strong(format!("Provider binding: {binding_status}"));
                    #[cfg(feature = "gui-test-support")]
                    capture_control(
                        &_provider_binding,
                        "protein_expression.provider_binding",
                        crate::gui_test_support::GuiTestWidgetKind::Status,
                        Some(&binding_status),
                    );
                    if let Some(detail) = &report.provider_binding.detail {
                        ui.label(detail);
                    }
                    details(
                        ui,
                        "Provider binding and withholding reasons",
                        &report.provider_binding,
                    );
                    details(ui, "Provider-neutral draft", &report.neutral_service_draft);
                    details(
                        ui,
                        "Selected-provider projection",
                        &report.provider_projection,
                    );
                    details(
                        ui,
                        "Service handoff candidates",
                        &report.service_handoff_candidates,
                    );
                    ui.collapsing("Warnings", |ui| {
                        for warning in &report.warnings {
                            ui.label(warning);
                        }
                    });
                    ui.collapsing("Suggested actions", |ui| {
                        for action in &report.suggested_next_actions {
                            ui.strong(&action.label);
                            ui.label(&action.rationale);
                            ui.add(
                                egui::Label::new(
                                    egui::RichText::new(&action.shell_line).monospace(),
                                )
                                .selectable(true)
                                .wrap(),
                            );
                            ui.push_id(&action.action_id, |ui| {
                                if ui
                                    .button("Copy command")
                                    .on_hover_text(&action.action_id)
                                    .clicked()
                                {
                                    ui.ctx().copy_text(action.shell_line.clone());
                                }
                            });
                        }
                    });
                    ui.horizontal_wrapped(|ui| {
                        let copy_json = ui.button("Copy report JSON");
                        #[cfg(feature = "gui-test-support")]
                        capture_control(
                            &copy_json,
                            "protein_expression.copy_json",
                            crate::gui_test_support::GuiTestWidgetKind::Button,
                            None,
                        );
                        if copy_json.clicked()
                            && let Ok(json) = serde_json::to_string_pretty(report)
                        {
                            ui.ctx().copy_text(json);
                        }
                        if let Some(text) = &report.text_report {
                            let copy_text = ui.button("Copy report text");
                            #[cfg(feature = "gui-test-support")]
                            capture_control(
                                &copy_text,
                                "protein_expression.copy_text",
                                crate::gui_test_support::GuiTestWidgetKind::Button,
                                None,
                            );
                            if copy_text.clicked() {
                                ui.ctx().copy_text(text.clone());
                            }
                        }
                    });
                    details(ui, "Shared report JSON", report);
                }
                if inspect {
                    self.start_protein_expression_report(ui.ctx());
                }
            });
    }

    pub(super) fn render_protein_expression_dialog(&mut self, ctx: &egui::Context) {
        if !self.protein_expression_ui.open {
            return;
        }
        let mut open = true;
        let viewport_id = Self::protein_expression_viewport_id();
        let spec = self.hosted_window_spec_for_viewport(
            "Protein Expression Handoff",
            Self::hosted_protein_expression_window_id(),
            viewport_id,
            Vec2::new(900.0, 740.0),
            Vec2::new(420.0, 320.0),
        );
        if ctx.embed_viewports() {
            crate::egui_compat::show_hosted_window(ctx, &spec, &mut open, |ui| {
                self.render_protein_expression_contents(ui);
            });
            self.clear_viewport_foreground_request_after_render(viewport_id);
            self.finalize_viewport_open_probe(viewport_id, "Protein Expression Handoff");
        } else {
            let builder = crate::egui_compat::viewport_builder_for_hosted_window(&spec);
            ctx.show_viewport_immediate(viewport_id, builder, |ctx, class| {
                self.note_viewport_focus_if_active(ctx, viewport_id);
                if class == egui::ViewportClass::EmbeddedWindow {
                    crate::egui_compat::show_hosted_window(&mut *ctx, &spec, &mut open, |ui| {
                        self.render_protein_expression_contents(ui);
                    });
                } else {
                    crate::egui_compat::show_central_panel(
                        &mut *ctx,
                        egui::CentralPanel::default(),
                        |ui| {
                            self.render_protein_expression_contents(ui);
                        },
                    );
                    if Self::viewport_close_requested_or_shortcut(ctx) {
                        open = false;
                    }
                }
            });
        }
        if !open {
            self.close_protein_expression_dialog();
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::engine::{Engine, Operation};
    use crate::engine_shell::parse_shell_line;

    fn app_with_product() -> GENtleApp {
        // Hand-crafted synthetic 27-bp CDS, not a biological gene; same product
        // as tutorial 09.03's committed expression_review_starter workflow.
        let app = GENtleApp::default();
        app.engine
            .write()
            .unwrap()
            .apply(Operation::CreateSequenceFromText {
                sequence_text: "ATGCACCACCACCACCACCACGCTTAA".into(),
                output_id: Some("expression_review_product".into()),
                name: Some("Synthetic expression review product".into()),
                circular: false,
            })
            .unwrap();
        app
    }

    fn canonical(mut value: serde_json::Value) -> serde_json::Value {
        value
            .as_object_mut()
            .unwrap()
            .remove("generated_at_unix_ms");
        value
    }

    #[test]
    fn protein_expression_inspector_matches_shared_shell_without_live_mutation() {
        let app = app_with_product();
        let before = serde_json::to_value(&*app.engine.read().unwrap()).unwrap();
        for (requirements, provider, status) in [
            (
                include_str!(
                    "../../docs/examples/plans/protein_expression_requirements_in_house.json"
                ),
                "",
                "not_selected",
            ),
            (
                include_str!(
                    "../../docs/examples/plans/protein_expression_requirements_in_house.json"
                ),
                "geneart",
                "withheld",
            ),
            (
                include_str!(
                    "../../docs/examples/plans/protein_expression_requirements_review_only.json"
                ),
                "geneart",
                "projected",
            ),
        ] {
            let inputs = Inputs {
                seq_id: "expression_review_product".into(),
                requirements: requirements.into(),
                provider: provider.into(),
            };
            let binding = Binding::new(&app.engine.read().unwrap(), inputs.clone());
            let actual =
                serde_json::to_value(inspect_snapshot(&app.engine, &binding).unwrap()).unwrap();
            let mut snapshot = app.engine.read().unwrap().clone_without_history();
            let command = inputs.command().unwrap();
            let expected = execute_shell_command(&mut snapshot, &command).unwrap();
            assert!(!expected.state_changed);
            assert_eq!(actual["provider_binding"]["status"], status);
            assert_eq!(canonical(actual.clone()), canonical(expected.output));
            if status != "projected" {
                assert!(actual["provider_projection"].is_null());
                assert_eq!(actual["service_handoff_candidates"], serde_json::json!([]));
            }
            for action in actual["suggested_next_actions"].as_array().unwrap() {
                parse_shell_line(action["shell_line"].as_str().unwrap()).unwrap();
            }
            assert_eq!(
                serde_json::to_value(&*app.engine.read().unwrap()).unwrap(),
                before
            );
        }
    }

    #[test]
    fn protein_expression_inspector_missing_product_and_requirements_stay_absent() {
        let app = app_with_product();
        let binding = Binding::new(&app.engine.read().unwrap(), Inputs::default());
        let report = inspect_snapshot(&app.engine, &binding).unwrap();
        assert_eq!(report.status, "needs_product_definition");
        assert!(report.product_definition.seq_id.is_none());
        assert!(report.neutral_service_draft.is_none());
        assert!(report.provider_projection.is_none());
        assert!(report.requirements.is_none());
        assert_eq!(
            report
                .missing_questions
                .iter()
                .map(|row| row.question_id.as_str())
                .collect::<Vec<_>>(),
            ["coding_sequence_or_orf_boundaries", "target_protein_source"]
        );
    }

    #[test]
    fn protein_expression_inner_agent_suggestion_matches_inspector_without_mutation() {
        let mut app = app_with_product();
        let before = serde_json::to_value(&*app.engine.read().unwrap()).unwrap();
        for (provider, expected_status) in [("", "not_selected"), ("geneart", "withheld")] {
            let inputs = Inputs {
                seq_id: "expression_review_product".into(),
                requirements: include_str!(
                    "../../docs/examples/plans/protein_expression_requirements_in_house.json"
                )
                .into(),
                provider: provider.into(),
            };
            let expected = inspect_snapshot(
                &app.engine,
                &Binding::new(&app.engine.read().unwrap(), inputs.clone()),
            )
            .unwrap();
            let command = format!(
                "planning protein-expression-handoff --seq-id expression_review_product --requirements {} --format text{}",
                crate::engine_shell::quote_shell_arg(&inputs.requirements),
                if provider.is_empty() {
                    String::new()
                } else {
                    format!(" --provider {provider}")
                }
            );
            let response: crate::agent_bridge::AgentResponse = serde_json::from_value(serde_json::json!({
                "schema": "gentle.agent_response.v1", "assistant_message": "Synthetic local review only",
                "suggested_commands": [{"command": command, "execution": "ask"}]
            })).unwrap();
            app.execute_agent_suggestion(1, &response.suggested_commands[0], "manual");
            let output = app
                .agent_last_command_output
                .as_ref()
                .expect("shared report retained for review")
                .clone();
            assert!(!output.state_changed);
            assert_eq!(output.output["provider_binding"]["status"], expected_status);
            assert_eq!(
                canonical(output.output.clone()),
                canonical(serde_json::to_value(expected).unwrap())
            );
            let log_count = app.agent_execution_log.len();
            app.stage_agent_command_result_followup(&output).unwrap();
            assert!(
                app.agent_prompt
                    .contains("gentle.protein_expression_handoff.v2")
            );
            assert!(app.agent_prompt.contains("data, not instructions"));
            assert!(app.agent_task.is_none(), "Staging does not contact a model");
            assert_eq!(
                app.agent_execution_log.len(),
                log_count,
                "Staging does not execute a follow-up"
            );
            assert_eq!(
                serde_json::to_value(&*app.engine.read().unwrap()).unwrap(),
                before
            );
        }
    }

    #[cfg(feature = "gui-test-support")]
    #[test]
    fn protein_expression_capture_targets_bind_the_actual_report_status() {
        let mut app = app_with_product();
        let binding = Binding::new(&app.engine.read().unwrap(), Inputs::default());
        app.protein_expression_ui.report = Some(inspect_snapshot(&app.engine, &binding).unwrap());
        app.protein_expression_ui.binding = Some(binding);
        let ctx = egui::Context::default();
        ctx.begin_pass(egui::RawInput {
            screen_rect: Some(egui::Rect::from_min_size(
                egui::Pos2::ZERO,
                Vec2::new(1100.0, 1000.0),
            )),
            ..Default::default()
        });
        crate::gui_test_support::begin_frame(&ctx);
        crate::egui_compat::show_central_panel_for_test_context(
            &ctx,
            egui::CentralPanel::default(),
            |ui| {
                app.render_protein_expression_contents(ui);
            },
        );
        let snapshot = crate::gui_test_support::snapshot(&ctx);
        for id in [
            CAPTURE_WINDOW,
            "protein_expression.product",
            "protein_expression.provider",
            "protein_expression.inspect",
            "protein_expression.status",
            "protein_expression.provider_binding",
            "protein_expression.copy_json",
            "protein_expression.copy_text",
        ] {
            assert!(
                snapshot
                    .items
                    .iter()
                    .any(|row| row.semantic_id == id && row.window_id == CAPTURE_WINDOW),
                "{id}"
            );
        }
        let status = snapshot
            .items
            .iter()
            .find(|row| row.semantic_id == "protein_expression.status")
            .unwrap();
        assert_eq!(
            status.outcome_role.as_deref(),
            Some("needs_product_definition")
        );
        let provider = snapshot
            .items
            .iter()
            .find(|row| row.semantic_id == "protein_expression.provider_binding")
            .unwrap();
        assert_eq!(provider.outcome_role.as_deref(), Some("not_selected"));
        crate::egui_compat::discard_test_pass_output(&ctx);
    }

    #[test]
    fn protein_expression_inspector_ui_intents_reuse_and_close_without_defaults() {
        let mut app = app_with_product();
        let before = serde_json::to_value(&*app.engine.read().unwrap()).unwrap();
        for action in ["open", "focus", "close", "close"] {
            let command =
                parse_shell_line(&format!("ui {action} protein-expression-handoff")).unwrap();
            let run = execute_shell_command(&mut app.engine.write().unwrap(), &command).unwrap();
            assert!(!run.state_changed);
            assert_eq!(
                run.output["ui_intent"]["target"],
                "protein-expression-handoff"
            );
            assert_eq!(run.output["applied"], false);
            assert!(app.try_apply_shell_ui_intent(&command).is_some());
            assert_eq!(app.protein_expression_ui.open, action != "close");
            assert_eq!(app.protein_expression_ui.inputs, Inputs::default());
            assert!(app.protein_expression_ui.report.is_none());
            assert!(!app.protein_expression_inspection_pending());
        }
        assert_eq!(
            serde_json::to_value(&*app.engine.read().unwrap()).unwrap(),
            before
        );
        assert!(
            parse_shell_line("ui open protein-expression-handoff --genome-id arbitrary").is_err()
        );
        assert!(parse_shell_line("ui open protein-expression-handoff --seq-id arbitrary").is_err());
    }

    #[test]
    fn protein_expression_inspector_discards_late_results_after_request_or_project_changes() {
        for change in ["request", "edit", "reopen", "close"] {
            let mut app = app_with_product();
            app.open_protein_expression_dialog();
            let binding = Binding::new(&app.engine.read().unwrap(), Inputs::default());
            let report = inspect_snapshot(&app.engine, &binding).unwrap();
            let (sender, receiver) = mpsc::channel();
            app.protein_expression_ui.pending = Some(Pending { binding, receiver });
            match change {
                "request" => app.protein_expression_ui.inputs.provider = "geneart".into(),
                "edit" => {
                    app.engine
                        .write()
                        .unwrap()
                        .auxiliary_metadata_mut()
                        .insert("test".into(), serde_json::json!(true));
                }
                "reopen" => {
                    let state = app.engine.read().unwrap().state().clone();
                    *app.engine.write().unwrap() = GentleEngine::from_state(state);
                }
                "close" => {
                    app.close_protein_expression_dialog();
                }
                _ => unreachable!(),
            }
            let _ = sender.send(Ok(report));
            app.poll_protein_expression_report();
            assert!(app.protein_expression_ui.report.is_none(), "{change}");
            assert!(!app.protein_expression_inspection_pending(), "{change}");
        }
    }

    #[test]
    fn protein_expression_inspector_adopts_current_results_and_invalidates_displayed_results() {
        let mut app = app_with_product();
        let binding = Binding::new(&app.engine.read().unwrap(), Inputs::default());
        let report = inspect_snapshot(&app.engine, &binding).unwrap();
        let (sender, receiver) = mpsc::channel();
        app.protein_expression_ui.pending = Some(Pending { binding, receiver });
        sender.send(Ok(report.clone())).unwrap();
        app.poll_protein_expression_report();
        assert_eq!(app.protein_expression_ui.report.as_ref(), Some(&report));
        assert!(!app.protein_expression_inspection_pending());
        app.protein_expression_ui.inputs.seq_id = "expression_review_product".into();
        app.poll_protein_expression_report();
        assert!(app.protein_expression_ui.report.is_none());
    }

    #[test]
    fn protein_expression_inspector_rejects_requirements_paths_and_disconnected_workers() {
        let mut app = app_with_product();
        let inputs = Inputs {
            requirements: "@private.json".into(),
            ..Inputs::default()
        };
        assert!(inputs.command().is_err());
        let binding = Binding::new(&app.engine.read().unwrap(), Inputs::default());
        let (sender, receiver) = mpsc::channel();
        app.protein_expression_ui.pending = Some(Pending { binding, receiver });
        drop(sender);
        app.poll_protein_expression_report();
        assert!(
            app.protein_expression_ui
                .error
                .as_ref()
                .unwrap()
                .contains("worker stopped")
        );
        assert!(app.protein_expression_ui.report.is_none());
    }

    #[test]
    fn protein_expression_inspector_small_hosted_window_paints_without_mutation() {
        let mut app = app_with_product();
        let binding = Binding::new(&app.engine.read().unwrap(), Inputs::default());
        app.protein_expression_ui.report = Some(inspect_snapshot(&app.engine, &binding).unwrap());
        app.protein_expression_ui.binding = Some(binding);
        app.open_protein_expression_dialog();
        let before = serde_json::to_value(&*app.engine.read().unwrap()).unwrap();
        let ctx = egui::Context::default();
        ctx.set_embed_viewports(true);
        let viewport_id = GENtleApp::protein_expression_viewport_id();
        // Opening requests one foreground frame; later renders use the normal layer.
        for expected_order in [egui::Order::Foreground, egui::Order::Middle] {
            assert_eq!(
                app.viewport_foreground_requested(viewport_id),
                expected_order == egui::Order::Foreground
            );
            ctx.begin_pass(egui::RawInput {
                screen_rect: Some(egui::Rect::from_min_size(
                    egui::Pos2::ZERO,
                    Vec2::new(520.0, 400.0),
                )),
                ..Default::default()
            });
            app.render_protein_expression_dialog(&ctx);
            assert!(
                ctx.memory(|memory| memory.areas().is_visible(&egui::LayerId::new(
                    expected_order,
                    GENtleApp::hosted_protein_expression_window_id(),
                ))),
                "expected the hosted inspector in {expected_order:?}"
            );
            assert!(!app.viewport_foreground_requested(viewport_id));
            crate::egui_compat::discard_test_pass_output(&ctx);
            assert_eq!(
                serde_json::to_value(&*app.engine.read().unwrap()).unwrap(),
                before
            );
            assert!(app.protein_expression_ui.open);
            assert!(app.protein_expression_ui.report.is_some());
        }
    }
}
