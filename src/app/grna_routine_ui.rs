//! Explicit gRNA routine setup through the existing shared macro preflight/executor.

use super::*;
use crate::engine_shell::routine_bindings::{
    grna_binding_readiness, routine_parameter_name, validate_routine_template_parameters,
};

#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub(super) enum GrnaRoutine {
    Anchor,
    Priority,
    Oligos,
}

impl GrnaRoutine {
    pub(super) const ALL: [Self; 3] = [Self::Anchor, Self::Priority, Self::Oligos];
    pub(super) fn id(self) -> &'static str {
        match self {
            Self::Anchor => "crispr.guides.anchor_window_scan",
            Self::Priority => "crispr.guides.candidate_priority_scan",
            Self::Oligos => "crispr.guides.practical_filter_and_oligos",
        }
    }
    pub(super) fn template(self) -> &'static str {
        match self {
            Self::Anchor => "grna_anchor_window_scan",
            Self::Priority => "grna_candidate_priority_scan",
            Self::Oligos => "grna_practical_filter_and_oligos",
        }
    }
    fn title(self) -> &'static str {
        match self {
            Self::Anchor => "Use gRNA Anchor Window Scan...",
            Self::Priority => "Use gRNA Candidate Priority Scan...",
            Self::Oligos => "Use gRNA Practical Filter and Oligos...",
        }
    }
    pub(super) fn validate(self, row: &CloningRoutineCatalogRow) -> Result<(), String> {
        if row.routine_id != self.id()
            || row.template_name != self.template()
            || row.status != "implemented"
        {
            return Err("Catalog routine identity or implementation status has changed".into());
        }
        let ports: &[(&str, &str)] = match self {
            Self::Anchor => &[
                ("seq_id", "sequence"),
                ("anchor_a", "sequence_anchor"),
                ("anchor_b", "sequence_anchor"),
            ],
            Self::Priority => &[("seq_id", "sequence")],
            Self::Oligos => &[("guide_set_id", "guide_set")],
        };
        for (id, kind) in ports {
            if !row.input_ports.iter().any(|port| {
                port["port_id"] == *id && port["kind"] == *kind && port["required"] == true
            }) {
                return Err(format!(
                    "Routine '{}' has no required {kind} port '{id}'",
                    self.id()
                ));
            }
        }
        let params = ports
            .iter()
            .map(|(port, _)| routine_parameter_name(self.template(), port))
            .collect::<Vec<_>>();
        validate_routine_template_parameters(
            row.template_path
                .as_deref()
                .ok_or("Missing template file")?,
            self.template(),
            &params,
        )
    }
}

pub(super) struct GrnaPreflightToken {
    engine: std::sync::Weak<RwLock<GentleEngine>>,
    revision: u64,
    template: String,
    bindings: HashMap<String, String>,
}

impl GENtleApp {
    pub(super) fn grna_action_readiness(
        &self,
        routine: GrnaRoutine,
        context: &LaunchContext,
    ) -> ActionReadiness {
        let input = context.readiness(CommandPaletteAction::UseGrnaRoutine(routine));
        if !input.is_ready() {
            return input;
        }
        self.pattern_catalog
            .grna_routine(routine)
            .map(|_| ActionReadiness::Ready)
            .unwrap_or_else(|reason| reason)
    }

    pub(super) fn palette_action_readiness(&self, action: CommandPaletteAction) -> ActionReadiness {
        let context = LaunchContext::capture(self, self.command_palette_dispatching);
        match action {
            CommandPaletteAction::UseGrnaRoutine(routine) => {
                self.grna_action_readiness(routine, &context)
            }
            _ => context.readiness(action),
        }
    }

    pub(super) fn append_grna_palette_entries(
        &self,
        entries: &mut Vec<CommandPaletteEntry>,
        context: &LaunchContext,
    ) {
        for routine in GrnaRoutine::ALL {
            entries.push(CommandPaletteEntry {
                title: routine.title().into(),
                detail: "Open bound Routine Assistant setup; no automatic execution".into(),
                keywords: "CRISPR guide RNA anchors oligos routine preflight".into(),
                action: CommandPaletteAction::UseGrnaRoutine(routine),
                readiness: self.grna_action_readiness(routine, context),
            });
        }
    }

    pub(super) fn render_grna_use_menu(&mut self, ui: &mut Ui, context: &LaunchContext) {
        self.pattern_catalog.ensure_started(ui.ctx());
        let states = GrnaRoutine::ALL.map(|routine| self.grna_action_readiness(routine, context));
        let parent = ActionReadiness::any_child(states.clone());
        ui.add_enabled_ui(parent.is_ready(), |ui| {
            ui.menu_button("Use gRNA routines", |ui| {
                let scans = ActionReadiness::any_child(states[..2].iter().cloned());
                ui.add_enabled_ui(scans.is_ready(), |ui| {
                    ui.menu_button("Candidate scans with selected DNA", |ui| {
                        for (routine, state) in GrnaRoutine::ALL[..2].iter().zip(&states) {
                            if state.button(ui, routine.title()).clicked() {
                                self.open_grna_routine(*routine);
                                ui.close();
                            }
                        }
                    });
                })
                .response
                .on_disabled_hover_text(scans.detail().unwrap_or_default());
                if states[2].button(ui, GrnaRoutine::Oligos.title()).clicked() {
                    self.open_grna_routine(GrnaRoutine::Oligos);
                    ui.close();
                }
            });
        })
        .response
        .on_disabled_hover_text(parent.detail().unwrap_or_default());
        if !parent.is_ready() {
            ui.small(
                "Use routines unavailable (hover for details). Template import remains available.",
            )
            .on_hover_text(parent.detail().unwrap_or_default());
        }
    }

    pub(super) fn open_grna_routine(&mut self, routine: GrnaRoutine) {
        let readiness =
            self.palette_action_readiness(CommandPaletteAction::UseGrnaRoutine(routine));
        if !readiness.is_ready() {
            self.app_status = readiness.detail().unwrap_or_default().into();
            return;
        }
        let Ok(row) = self.pattern_catalog.grna_routine(routine) else {
            return;
        };
        // Revalidate the file only at invocation, never in a menu frame.
        if let Err(error) = routine.validate(&row) {
            self.app_status =
                format!("Routine unavailable: {error}. Refresh the template catalog.");
            return;
        }
        let mut bindings = BTreeMap::new();
        if routine != GrnaRoutine::Oligos {
            let Ok((seq_id, _)) = self.selected_sequence_context(true) else {
                return;
            };
            bindings.insert("seq_id".into(), seq_id);
        }
        self.routine_assistant_candidates = vec![row];
        self.routine_assistant_selected_routine_id = routine.id().into();
        self.routine_assistant_bindings = bindings;
        self.routine_assistant_preflight_output = None;
        self.routine_assistant_execute_output = None;
        self.grna_preflight_token = None;
        self.routine_assistant_explain_output = None;
        self.routine_assistant_compare_output = None;
        self.routine_assistant_decision_trace = None;
        self.routine_assistant_stage = RoutineAssistantStage::Parameters;
        self.sync_routine_assistant_bindings_for_selected();
        self.open_routine_assistant_dialog();
        self.routine_assistant_status = "Choose bindings, then run preflight. Candidate scans are generic preselection, not PAM-aware guide design or specificity confirmation.".into();
        self.app_status = format!("Opened {} setup", routine.title().trim_end_matches("..."));
    }

    pub(super) fn grna_form_readiness(&self) -> ActionReadiness {
        let Some(row) = self.routine_assistant_selected_routine() else {
            return ActionReadiness::Ready;
        };
        if !GrnaRoutine::ALL
            .iter()
            .any(|r| r.template() == row.template_name)
        {
            return ActionReadiness::Ready;
        }
        let Ok(engine) = self.engine.try_read() else {
            return ActionReadiness::Checking {
                detail: "Project is busy".into(),
            };
        };
        match grna_binding_readiness(
            &engine,
            &row.template_name,
            &self.routine_assistant_bindings_compact(),
        ) {
            Ok(()) => ActionReadiness::Ready,
            Err(detail) => ActionReadiness::NeedsBindings { detail },
        }
    }

    pub(super) fn capture_grna_preflight(&mut self) {
        self.grna_preflight_token = self.routine_assistant_selected_routine().and_then(|row| {
            if !GrnaRoutine::ALL
                .iter()
                .any(|r| r.template() == row.template_name)
            {
                return None;
            }
            Some(GrnaPreflightToken {
                engine: Arc::downgrade(&self.engine),
                revision: self.engine.try_read().ok()?.mutation_revision(),
                template: row.template_name,
                bindings: self.routine_assistant_bindings_compact(),
            })
        });
    }

    pub(super) fn grna_preflight_current(&self) -> bool {
        let Some(row) = self.routine_assistant_selected_routine() else {
            return true;
        };
        if !GrnaRoutine::ALL
            .iter()
            .any(|r| r.template() == row.template_name)
        {
            return true;
        }
        self.routine_assistant_preflight_output
            .as_ref()
            .and_then(|v| v.get("can_execute"))
            .and_then(|v| v.as_bool())
            == Some(true)
            && self.grna_form_readiness().is_ready()
            && self.grna_preflight_token.as_ref().is_some_and(|token| {
                token.engine.ptr_eq(&Arc::downgrade(&self.engine))
                    && token.template == row.template_name
                    && token.bindings == self.routine_assistant_bindings_compact()
                    && self
                        .engine
                        .try_read()
                        .is_ok_and(|e| e.mutation_revision() == token.revision)
            })
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn app_with_dna() -> GENtleApp {
        let mut state = ProjectState::default();
        state.sequences.insert(
            "dna".into(),
            DNAsequence::from_sequence(&"ACGT".repeat(12)).unwrap(),
        );
        let mut app = GENtleApp {
            engine: Arc::new(RwLock::new(GentleEngine::from_state(state))),
            pattern_catalog: pattern_catalog_ui::load_test_catalog(),
            ..Default::default()
        };
        app.lineage_graph_selected_node_id = app
            .engine
            .read()
            .unwrap()
            .state()
            .lineage
            .seq_to_node
            .get("dna")
            .cloned();
        app
    }

    #[test]
    fn catalog_identity_validation_preserves_ids_and_rejects_wrong_ports() {
        let cache = pattern_catalog_ui::load_test_catalog();
        for routine in GrnaRoutine::ALL {
            let mut row = cache.grna_routine(routine).unwrap();
            routine.validate(&row).unwrap();
            assert_eq!(row.routine_id, routine.id());
            row.input_ports.clear();
            assert!(routine.validate(&row).is_err());
        }
    }

    #[test]
    fn anchor_setup_requires_bindings_then_shared_preflight_and_rejects_stale_run() {
        let mut app = app_with_dna();
        app.open_grna_routine(GrnaRoutine::Anchor);
        assert!(app.show_routine_assistant_dialog);
        assert_eq!(app.routine_assistant_bindings["seq_id"], "dna");
        assert!(!app.grna_form_readiness().is_ready());
        let before = app.engine.read().unwrap().mutation_revision();
        app.run_routine_assistant_preflight();
        assert_eq!(before, app.engine.read().unwrap().mutation_revision());
        app.routine_assistant_bindings
            .insert("anchor_a".into(), "0".into());
        app.routine_assistant_bindings
            .insert("anchor_b".into(), "48".into());
        assert!(app.grna_form_readiness().is_ready());
        let bindings = app.routine_assistant_bindings_compact();
        assert_eq!(bindings["anchor_a_pos"], "0");
        assert!(!bindings.contains_key("anchor_a"));
        app.run_routine_assistant_preflight();
        assert!(
            app.routine_assistant_can_execute(),
            "{} {:?}",
            app.routine_assistant_status,
            app.routine_assistant_preflight_output
        );
        app.engine
            .write()
            .unwrap()
            .state_mut()
            .metadata
            .insert("synthetic_edit".into(), serde_json::json!(true));
        assert!(!app.routine_assistant_can_execute());
        let before = app.engine.read().unwrap().mutation_revision();
        app.run_routine_assistant_execute();
        assert_eq!(before, app.engine.read().unwrap().mutation_revision());
        app.run_routine_assistant_preflight();
        assert!(
            app.routine_assistant_can_execute(),
            "{}",
            app.routine_assistant_status
        );
        app.run_routine_assistant_execute();
        assert!(
            app.routine_assistant_execute_output.is_some(),
            "{}",
            app.routine_assistant_status
        );
    }

    #[test]
    fn priority_scan_launch_uses_selected_dna_without_anchor_requirement() {
        let mut app = app_with_dna();
        app.open_grna_routine(GrnaRoutine::Priority);
        assert!(app.grna_form_readiness().is_ready());
        app.run_routine_assistant_preflight();
        assert!(
            app.routine_assistant_can_execute(),
            "{}",
            app.routine_assistant_status
        );
    }

    #[test]
    fn oligo_route_requires_guide_set_but_not_dna_and_never_picks_first_set() {
        let mut app = app_with_dna();
        assert!(
            !app.grna_action_readiness(GrnaRoutine::Oligos, &LaunchContext::capture(&app, false))
                .is_ready()
        );
        let mut engine = GentleEngine::from_state(ProjectState::default());
        engine
            .apply(Operation::UpsertGuideSet {
                guide_set_id: "guides".into(),
                guides: vec![crate::engine::protocol::GuideCandidate {
                    guide_id: "g1".into(),
                    seq_id: "external".into(),
                    start_0based: 0,
                    end_0based_exclusive: 20,
                    strand: "+".into(),
                    protospacer: "GACCTGTTGACGATGTTCCA".into(),
                    pam: "AGG".into(),
                    nuclease: "SpCas9".into(),
                    cut_offset_from_protospacer_start: 17,
                    rank: None,
                }],
            })
            .unwrap();
        app.engine = Arc::new(RwLock::new(engine));
        app.lineage_graph_selected_node_id = None;
        assert!(
            app.grna_action_readiness(GrnaRoutine::Oligos, &LaunchContext::capture(&app, false))
                .is_ready()
        );
        app.open_grna_routine(GrnaRoutine::Oligos);
        assert_eq!(app.routine_assistant_bindings["guide_set_id"], "");
        assert!(!app.grna_form_readiness().is_ready());
        app.routine_assistant_bindings
            .insert("guide_set_id".into(), "guides".into());
        assert!(app.grna_form_readiness().is_ready());
        app.run_routine_assistant_preflight();
        assert!(
            app.routine_assistant_can_execute(),
            "{} {:?}",
            app.routine_assistant_status,
            app.routine_assistant_preflight_output
        );
    }
}
