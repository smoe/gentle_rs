//! Shared launcher presentation; availability is not permission or scientific preflight.

use super::*;
use gentle_protocol::CollectionLiftRejectionReason;

impl CommandPaletteEntry {
    pub(super) fn ready_action(&self) -> Option<CommandPaletteAction> {
        self.readiness.is_ready().then_some(self.action)
    }

    pub(super) fn render_row(&self, ui: &mut Ui, selected: bool, label: String) -> egui::Response {
        let response = ui
            .add_enabled(
                self.readiness.is_ready(),
                egui::Button::new(label).selected(selected),
            )
            .on_disabled_hover_text(self.readiness.detail().unwrap_or_default());
        if selected && let Some(detail) = self.readiness.detail() {
            ui.small(format!("{}: {detail}", self.readiness.label()));
        }
        response
    }
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub(super) enum ActionReadiness {
    Ready,
    Checking {
        detail: String,
    },
    NeedsInput {
        detail: String,
    },
    NeedsBindings {
        detail: String,
    },
    RequiresMaterialization {
        reason: CollectionLiftRejectionReason,
        detail: String,
    },
    RequiresPhysicalPool {
        reason: CollectionLiftRejectionReason,
        detail: String,
    },
    PolicyRejected {
        reason: CollectionLiftRejectionReason,
        detail: String,
    },
    AdapterUnavailable {
        detail: String,
    },
}

impl ActionReadiness {
    pub(super) fn is_ready(&self) -> bool {
        matches!(self, Self::Ready)
    }

    pub(super) fn label(&self) -> &'static str {
        match self {
            Self::Ready => "Ready",
            Self::Checking { .. } => "Checking",
            Self::NeedsInput { .. } => "Needs input",
            Self::NeedsBindings { .. } => "Needs bindings",
            Self::RequiresMaterialization { .. } => "Requires materialization",
            Self::RequiresPhysicalPool { .. } => "Requires physical pool",
            Self::PolicyRejected { .. } => "Unsupported",
            Self::AdapterUnavailable { .. } => "GUI adapter unavailable",
        }
    }

    pub(super) fn detail(&self) -> Option<&str> {
        match self {
            Self::Ready => None,
            Self::Checking { detail }
            | Self::NeedsInput { detail }
            | Self::NeedsBindings { detail }
            | Self::RequiresMaterialization { detail, .. }
            | Self::RequiresPhysicalPool { detail, .. }
            | Self::PolicyRejected { detail, .. }
            | Self::AdapterUnavailable { detail } => Some(detail),
        }
    }

    pub(super) fn rejection_reason(&self) -> Option<CollectionLiftRejectionReason> {
        match self {
            Self::RequiresMaterialization { reason, .. }
            | Self::RequiresPhysicalPool { reason, .. }
            | Self::PolicyRejected { reason, .. } => Some(*reason),
            _ => None,
        }
    }

    pub(super) fn any_child(children: impl IntoIterator<Item = Self>) -> Self {
        let children = children.into_iter().collect::<Vec<_>>();
        if children.iter().any(Self::is_ready) {
            return Self::Ready;
        }
        let detail = children
            .iter()
            .filter_map(Self::detail)
            .collect::<Vec<_>>()
            .join("; ");
        Self::NeedsInput {
            detail: if detail.is_empty() {
                "No available actions".into()
            } else {
                detail
            },
        }
    }

    pub(super) fn button(&self, ui: &mut Ui, label: impl Into<egui::WidgetText>) -> egui::Response {
        ui.add_enabled(self.is_ready(), egui::Button::new(label))
            .on_disabled_hover_text(self.detail().unwrap_or("No available action"))
    }
}

// Small per-surface snapshot: no project graph, template parsing or host probes.
pub(super) struct LaunchContext {
    pub(super) dna: ActionReadiness,
    sequence: ActionReadiness,
    pub(super) guides: ActionReadiness,
    empty: bool,
    pcr_open: bool,
    confirmation_open: bool,
}

impl LaunchContext {
    pub(super) fn capture(app: &GENtleApp, palette: bool) -> Self {
        let selected = app.selected_sequence_context_with_origin(false, palette);
        let engine = app.engine.try_read();
        let sequence = match &selected {
            Ok(_) => ActionReadiness::Ready,
            Err(detail) => ActionReadiness::NeedsInput {
                detail: detail.clone(),
            },
        };
        let (dna, guides, empty) = match engine {
            Ok(engine) => {
                let dna = match &selected {
                    Ok((id, _)) if engine.sequence_kind(id) == Some("dna") => {
                        ActionReadiness::Ready
                    }
                    Ok((id, _)) => ActionReadiness::NeedsInput {
                        detail: format!(
                            "'{id}' is not DNA; select DNA in the project graph or its viewer"
                        ),
                    },
                    Err(_) => sequence.clone(),
                };
                let guides = if engine.guide_set_input_ids().is_empty() {
                    ActionReadiness::NeedsInput { detail: "A stored guide set is required, not a DNA sequence or generic candidate set. Import one with guides put, then choose its ID in the form.".into() }
                } else {
                    ActionReadiness::Ready
                };
                (dna, guides, engine.state().sequences.is_empty())
            }
            Err(_) => {
                let checking = ActionReadiness::Checking {
                    detail: "Project is busy; try again".into(),
                };
                (checking.clone(), checking, false)
            }
        };
        Self {
            dna,
            sequence,
            guides,
            empty,
            pcr_open: app.show_pcr_design_dialog && !app.pcr_design_seq_id.is_empty(),
            confirmation_open: app.show_sequencing_confirmation_dialog,
        }
    }

    pub(super) fn readiness(&self, action: CommandPaletteAction) -> ActionReadiness {
        match action {
            CommandPaletteAction::OpenCrypticSplicingScreen
            | CommandPaletteAction::OpenTataBoxes
            | CommandPaletteAction::OpenTssInventory => self.dna.clone(),
            CommandPaletteAction::OpenGenomicRegionConservation if !self.empty => self.dna.clone(),
            CommandPaletteAction::UiIntent(UiIntentTarget::FeatureLocationEditor) => {
                self.sequence.clone()
            }
            CommandPaletteAction::UiIntent(UiIntentTarget::SavedGenomicRegions) => self.dna.clone(),
            CommandPaletteAction::UiIntent(UiIntentTarget::PcrDesign) if !self.pcr_open => {
                self.dna.clone()
            }
            CommandPaletteAction::UiIntent(UiIntentTarget::SequencingConfirmation)
                if !self.confirmation_open =>
            {
                self.dna.clone()
            }
            CommandPaletteAction::UseGrnaRoutine(grna_routine_ui::GrnaRoutine::Oligos) => {
                self.guides.clone()
            }
            CommandPaletteAction::UseGrnaRoutine(_) => self.dna.clone(),
            // Non-migrated actions retain their existing checks; this is not a global registry.
            _ => ActionReadiness::Ready,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn select(app: &mut GENtleApp, id: &str) {
        app.lineage_graph_selected_node_id = app
            .engine
            .read()
            .unwrap()
            .state()
            .lineage
            .seq_to_node
            .get(id)
            .cloned();
    }

    #[test]
    fn project_transitions_recompute_cheap_readiness_without_mutation() {
        let mut state = ProjectState::default();
        let mut protein = DNAsequence::from_sequence("ACGTACGT").unwrap();
        protein.set_molecule_type("protein");
        state.sequences.insert("a_protein".into(), protein);
        state.sequences.insert(
            "dna".into(),
            DNAsequence::from_sequence("ACGTACGT").unwrap(),
        );
        let mut app = GENtleApp {
            engine: Arc::new(RwLock::new(GentleEngine::from_state(state))),
            ..Default::default()
        };
        let action = CommandPaletteAction::OpenTataBoxes;
        assert!(
            !LaunchContext::capture(&app, false)
                .readiness(action)
                .is_ready()
        );
        select(&mut app, "a_protein");
        assert!(
            !LaunchContext::capture(&app, false)
                .readiness(action)
                .is_ready()
        );
        select(&mut app, "dna");
        let before = app.engine.read().unwrap().mutation_revision();
        for _ in 0..20 {
            assert!(
                LaunchContext::capture(&app, false)
                    .readiness(action)
                    .is_ready()
            );
        }
        assert_eq!(before, app.engine.read().unwrap().mutation_revision());
        let engine = app.engine.clone();
        let locked = engine.write().unwrap();
        assert!(matches!(
            LaunchContext::capture(&app, false).dna,
            ActionReadiness::Checking { .. }
        ));
        drop(locked);
        app.engine
            .write()
            .unwrap()
            .state_mut()
            .sequences
            .remove("dna");
        assert!(!app.execute_command_palette_action(&egui::Context::default(), action));
    }

    #[test]
    fn recovery_children_keep_parent_enabled() {
        let blocked = ActionReadiness::NeedsInput {
            detail: "Select DNA".into(),
        };
        assert!(!ActionReadiness::any_child([]).is_ready());
        assert!(!ActionReadiness::any_child([blocked.clone()]).is_ready());
        assert!(ActionReadiness::any_child([blocked, ActionReadiness::Ready]).is_ready());
    }

    #[test]
    fn disabled_palette_row_rejects_mouse_and_enter() {
        let entry = CommandPaletteEntry {
            title: "Blocked".into(),
            detail: "test".into(),
            keywords: String::new(),
            action: CommandPaletteAction::OpenTataBoxes,
            readiness: ActionReadiness::NeedsInput {
                detail: "Select DNA".into(),
            },
        };
        assert!(entry.ready_action().is_none());
        let ctx = egui::Context::default();
        let mut pos = egui::Pos2::ZERO;
        ctx.run_ui(Default::default(), |ui| {
            pos = entry
                .render_row(ui, true, entry.title.clone())
                .rect
                .center();
        })
        .drop_without_applying_deltas();
        for pressed in [true, false] {
            let input = egui::RawInput {
                events: vec![
                    egui::Event::PointerMoved(pos),
                    egui::Event::PointerButton {
                        pos,
                        button: egui::PointerButton::Primary,
                        pressed,
                        modifiers: egui::Modifiers::NONE,
                    },
                ],
                ..Default::default()
            };
            ctx.run_ui(input, |ui| {
                let response = entry.render_row(ui, true, entry.title.clone());
                assert!(!response.enabled());
                assert!(!response.clicked());
            })
            .drop_without_applying_deltas();
        }
    }

    #[test]
    fn palette_enter_in_both_window_paths_keeps_blocked_action_open() {
        // Egui's immediate-renderer callback is thread-local; do not leak it to other tests.
        std::thread::spawn(|| {
            for embedded in [true, false] {
                let ctx = egui::Context::default();
                ctx.set_embed_viewports(embedded);
                let calls = std::rc::Rc::new(std::cell::Cell::new(0));
                let counted = calls.clone();
                egui::Context::set_immediate_viewport_renderer(move |ctx, mut viewport| {
                    counted.set(counted.get() + 1);
                    let mut input = egui::RawInput {
                        viewport_id: viewport.ids.this,
                        events: vec![egui::Event::Key {
                            key: Key::Enter,
                            physical_key: None,
                            pressed: true,
                            repeat: false,
                            modifiers: egui::Modifiers::NONE,
                        }],
                        ..Default::default()
                    };
                    input.viewports.insert(
                        viewport.ids.this,
                        egui::ViewportInfo {
                            parent: Some(viewport.ids.parent),
                            ..Default::default()
                        },
                    );
                    ctx.run_ui(input, |ui| (viewport.viewport_ui_cb)(ui))
                        .drop_without_applying_deltas();
                });
                let mut app = GENtleApp::default();
                app.open_command_palette_dialog();
                app.command_palette_query = "TATA".into();
                let before = app.engine.read().unwrap().mutation_revision();
                let input = egui::RawInput {
                    events: vec![egui::Event::Key {
                        key: Key::Enter,
                        physical_key: None,
                        pressed: true,
                        repeat: false,
                        modifiers: egui::Modifiers::NONE,
                    }],
                    ..Default::default()
                };
                ctx.run_ui(input, |_| app.render_command_palette_dialog(&ctx))
                    .drop_without_applying_deltas();
                assert!(app.show_command_palette_dialog);
                assert_eq!(before, app.engine.read().unwrap().mutation_revision());
                assert!(
                    !app.execute_command_palette_action(&ctx, CommandPaletteAction::OpenTataBoxes)
                );
                assert_eq!(calls.get() > 0, !embedded);
            }
        })
        .join()
        .unwrap();
    }
}
