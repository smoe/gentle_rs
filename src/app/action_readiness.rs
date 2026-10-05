//! Shared launcher presentation; availability is not permission or scientific preflight.

use super::*;
use gentle_protocol::CollectionLiftRejectionReason;

impl CommandPaletteEntry {
    pub(super) fn ready_action(&self) -> Option<CommandPaletteAction> {
        self.readiness.is_ready().then_some(self.action)
    }

    pub(super) fn render_row(&self, ui: &mut Ui, selected: bool) -> egui::Response {
        ui.add_enabled(
            self.readiness.is_ready(),
            egui::Button::new(&self.title).selected(selected),
        )
        .on_hover_text(&self.detail)
        .on_disabled_hover_text(self.readiness.detail().unwrap_or_default())
    }

    fn render_detail(&self, ui: &mut Ui) {
        ui.add(egui::Label::new(egui::RichText::new(&self.title).strong()).wrap());
        ui.add(egui::Label::new(&self.detail).wrap());
        if let Some(detail) = self.readiness.detail() {
            ui.add(egui::Label::new(format!("{}: {detail}", self.readiness.label())).wrap());
        }
        ui.small(if self.readiness.is_ready() {
            "Enter to choose this action · Esc to close"
        } else {
            "Choose another action, or follow the guidance above · Esc to close"
        });
    }
}

impl GENtleApp {
    // Both hosted and native palettes use this presentation and keyboard policy.
    pub(super) fn render_command_palette_contents(
        &mut self,
        ui: &mut Ui,
        entries: &[CommandPaletteEntry],
    ) -> Option<CommandPaletteAction> {
        crate::agent_help::render_agent_help_button(
            ui,
            "Command Palette",
            "window.command_palette",
        );
        ui.separator();
        ui.label("Search actions, settings, and help topics");
        let input_id = ui.make_persistent_id("gentle_command_palette_search");
        let search_response = ui.add(
            egui::TextEdit::singleline(&mut self.command_palette_query)
                .id(input_id)
                .desired_width(f32::INFINITY)
                .hint_text("Type action name (Cmd/Ctrl+K)"),
        );
        if self.command_palette_focus_query {
            search_response.request_focus();
            self.command_palette_focus_query = false;
        }
        let previous_selection = self.command_palette_selected;
        let (down, up, enter, pointer_moved) = ui.input(|i| {
            (
                i.key_pressed(Key::ArrowDown),
                i.key_pressed(Key::ArrowUp),
                i.key_pressed(Key::Enter),
                i.events
                    .iter()
                    .any(|event| matches!(event, egui::Event::PointerMoved(_))),
            )
        });
        if !entries.is_empty() {
            self.command_palette_selected = self.command_palette_selected.min(entries.len() - 1);
            if down {
                self.command_palette_selected = (self.command_palette_selected + 1) % entries.len();
            }
            if up {
                self.command_palette_selected =
                    (self.command_palette_selected + entries.len() - 1) % entries.len();
            }
        }

        ui.separator();
        if entries.is_empty() {
            self.command_palette_selected = 0;
            ui.small("No matching commands");
            return None;
        }
        let available = ui.available_height();
        let detail_height = (available * 0.4).clamp(96.0, 150.0).min(available);
        let results_height =
            (available - detail_height - ui.spacing().item_spacing.y * 2.0).max(0.0);
        let mut chosen = None;
        egui::ScrollArea::vertical()
            .id_salt("command_palette_results_scroll")
            .max_height(results_height)
            .auto_shrink([false, false])
            .show(ui, |ui| {
                scroll_input_policy::apply_scrollarea_keyboard_navigation(
                    ui,
                    scroll_input_policy::DEFAULT_SCROLLAREA_KEYBOARD_STEP,
                );
                for (idx, entry) in entries.iter().enumerate() {
                    let selected = self.command_palette_selected == idx;
                    let response = entry.render_row(ui, selected);
                    if selected && (down || up) {
                        response.scroll_to_me(Some(egui::Align::Center));
                    }
                    // Disabled rows can explain themselves without becoming executable.
                    // A stationary pointer must not undo a keyboard selection.
                    if response.contains_pointer() && pointer_moved && !down && !up {
                        self.command_palette_selected = idx;
                        self.hover_status_name = format!("Command palette: {}", entry.title);
                    }
                    if response.clicked() {
                        chosen = entry.ready_action();
                    }
                }
            });
        ui.separator();
        let mut details = egui::ScrollArea::vertical()
            .id_salt("command_palette_detail_scroll")
            .max_height(detail_height);
        if self.command_palette_selected != previous_selection {
            details = details.vertical_scroll_offset(0.0);
        }
        details.show(ui, |ui| {
            entries[self.command_palette_selected].render_detail(ui)
        });
        chosen.or_else(|| {
            enter
                .then(|| entries[self.command_palette_selected].ready_action())
                .flatten()
        })
    }
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub(super) enum ActionReadiness {
    Ready,
    Blocked {
        detail: String,
    },
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
            Self::Blocked { .. } => "Unavailable",
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
            Self::Blocked { detail }
            | Self::Checking { detail }
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

    fn history_action(available: usize, background_jobs: bool, action: &str) -> Self {
        if background_jobs {
            Self::Blocked {
                detail: format!(
                    "{action} is unavailable while background jobs are active. Wait for them to finish, or inspect/cancel them in Background Jobs."
                ),
            }
        } else if available == 0 {
            Self::Blocked {
                detail: if action == "Undo" {
                    "Nothing to undo in this session. Make a project change first; saved history is not restored when reopening a project.".into()
                } else {
                    "Nothing to redo in this session. Undo a change first; making a new change clears redo.".into()
                },
            }
        } else {
            Self::Ready
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
    undo: ActionReadiness,
    redo: ActionReadiness,
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
        let (dna, guides, empty, undo, redo) = match engine {
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
                let background_jobs = app.has_active_background_jobs();
                (
                    dna,
                    guides,
                    engine.state().sequences.is_empty(),
                    ActionReadiness::history_action(
                        engine.undo_available(),
                        background_jobs,
                        "Undo",
                    ),
                    ActionReadiness::history_action(
                        engine.redo_available(),
                        background_jobs,
                        "Redo",
                    ),
                )
            }
            Err(_) => {
                let checking = ActionReadiness::Checking {
                    detail: "Project is busy; try again".into(),
                };
                (
                    checking.clone(),
                    checking.clone(),
                    false,
                    checking.clone(),
                    checking,
                )
            }
        };
        Self {
            dna,
            sequence,
            guides,
            empty,
            pcr_open: app.show_pcr_design_dialog && !app.pcr_design_seq_id.is_empty(),
            confirmation_open: app.show_sequencing_confirmation_dialog,
            undo,
            redo,
        }
    }

    pub(super) fn readiness(&self, action: CommandPaletteAction) -> ActionReadiness {
        match action {
            CommandPaletteAction::Undo => self.undo.clone(),
            CommandPaletteAction::Redo => self.redo.clone(),
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

    // Synthetic in-memory project/history and egui inputs, recreated by these tests.
    // No file fixtures, provider calls or screenshot capture are involved.
    #[test]
    fn history_actions_follow_checkpoints_and_reject_stale_palette_entries() {
        let ctx = egui::Context::default();
        let mut app = GENtleApp::default();
        for action in [CommandPaletteAction::Undo, CommandPaletteAction::Redo] {
            let before = app.engine.read().unwrap().mutation_revision();
            assert!(!app.palette_action_readiness(action).is_ready());
            assert!(!app.execute_command_palette_action(&ctx, action));
            assert!(app.app_status.contains("Nothing to"));
            assert_eq!(before, app.engine.read().unwrap().mutation_revision());
        }
        app.engine
            .write()
            .unwrap()
            .apply(Operation::SetDisplayVisibility {
                target: crate::engine::DisplayTarget::Features,
                visible: false,
            })
            .unwrap();
        let entries = app.collect_command_palette_entries();
        assert!(
            entries
                .iter()
                .find(|entry| matches!(entry.action, CommandPaletteAction::Undo))
                .unwrap()
                .ready_action()
                .is_some()
        );
        assert!(
            entries
                .iter()
                .find(|entry| matches!(entry.action, CommandPaletteAction::Redo))
                .unwrap()
                .ready_action()
                .is_none()
        );
        assert!(app.execute_command_palette_action(&ctx, CommandPaletteAction::Undo));
        assert!(app.engine.read().unwrap().state().display.show_features);
        assert!(
            !app.palette_action_readiness(CommandPaletteAction::Undo)
                .is_ready()
        );
        assert!(
            app.palette_action_readiness(CommandPaletteAction::Redo)
                .is_ready()
        );
        assert!(app.execute_command_palette_action(&ctx, CommandPaletteAction::Redo));
        assert!(!app.engine.read().unwrap().state().display.show_features);

        assert!(app.execute_command_palette_action(&ctx, CommandPaletteAction::Undo));
        let stale_redo = app
            .collect_command_palette_entries()
            .into_iter()
            .find(|entry| matches!(entry.action, CommandPaletteAction::Redo))
            .unwrap();
        assert!(stale_redo.ready_action().is_some());
        app.engine
            .write()
            .unwrap()
            .apply(Operation::SetDisplayVisibility {
                target: crate::engine::DisplayTarget::Features,
                visible: false,
            })
            .unwrap();
        let before = app.engine.read().unwrap().mutation_revision();
        assert!(!app.execute_command_palette_action(&ctx, stale_redo.action));
        assert!(app.app_status.contains("new change clears redo"));
        assert_eq!(before, app.engine.read().unwrap().mutation_revision());
        assert!(app.new_windows.is_empty());
    }

    #[test]
    fn history_actions_explain_background_jobs_and_never_wait_for_busy_project() {
        let ctx = egui::Context::default();
        let mut app = GENtleApp::default();
        for visible in [false, true] {
            app.engine
                .write()
                .unwrap()
                .apply(Operation::SetDisplayVisibility {
                    target: crate::engine::DisplayTarget::Features,
                    visible,
                })
                .unwrap();
        }
        app.engine.write().unwrap().undo_last_operation().unwrap();
        assert!(
            app.palette_action_readiness(CommandPaletteAction::Undo)
                .is_ready()
        );
        assert!(
            app.palette_action_readiness(CommandPaletteAction::Redo)
                .is_ready()
        );
        let (_tx, receiver) = mpsc::channel();
        app.agent_task = Some(AgentAskTask {
            job_id: 1,
            prompt: "synthetic background job".into(),
            attachment_summaries: vec![],
            _attachment_files: vec![],
            started: Instant::now(),
            runtime_frame: crate::runtime_status::runtime_status_registry().push_with_detail(
                crate::runtime_status::RuntimeStatusFrameKind::BackgroundJob,
                "synthetic history readiness test",
                None,
            ),
            receiver,
        });
        let before = app.engine.read().unwrap().mutation_revision();
        for action in [CommandPaletteAction::Undo, CommandPaletteAction::Redo] {
            let menu = LaunchContext::capture(&app, false).readiness(action);
            let palette = app.palette_action_readiness(action);
            assert_eq!(menu, palette);
            assert_eq!(palette.label(), "Unavailable");
            assert!(palette.detail().unwrap().contains("Background Jobs"));
            assert!(!app.execute_command_palette_action(&ctx, action));
        }
        assert_eq!(before, app.engine.read().unwrap().mutation_revision());
        assert!(
            app.palette_action_readiness(CommandPaletteAction::OpenSequence)
                .is_ready()
        );
        app.agent_task = None;
        let engine = app.engine.clone();
        let locked = engine.write().unwrap();
        for action in [CommandPaletteAction::Undo, CommandPaletteAction::Redo] {
            assert!(matches!(
                app.palette_action_readiness(action),
                ActionReadiness::Checking { .. }
            ));
            assert!(!app.execute_command_palette_action(&ctx, action));
        }
        drop(locked);
        assert!(
            app.palette_action_readiness(CommandPaletteAction::Undo)
                .is_ready()
        );
        assert!(
            app.palette_action_readiness(CommandPaletteAction::Redo)
                .is_ready()
        );
    }

    fn palette_entry(
        title: &str,
        action: CommandPaletteAction,
        readiness: ActionReadiness,
    ) -> CommandPaletteEntry {
        CommandPaletteEntry {
            title: title.into(),
            detail: format!("Open {title} setup; no analysis runs automatically."),
            keywords: String::new(),
            action,
            readiness,
        }
    }

    fn small_palette_input(events: Vec<egui::Event>) -> egui::RawInput {
        egui::RawInput {
            screen_rect: Some(egui::Rect::from_min_size(
                egui::Pos2::ZERO,
                Vec2::new(500.0, 320.0),
            )),
            events,
            ..Default::default()
        }
    }

    fn key_event(key: Key) -> egui::Event {
        egui::Event::Key {
            key,
            physical_key: None,
            pressed: true,
            repeat: false,
            modifiers: egui::Modifiers::NONE,
        }
    }

    fn text_shapes(shape: &egui::epaint::Shape, out: &mut Vec<(String, egui::Rect)>) {
        match shape {
            egui::epaint::Shape::Text(text) => out.push((
                text.galley.job.text.clone(),
                egui::Rect::from_min_size(text.pos, text.galley.size()),
            )),
            egui::epaint::Shape::Vec(shapes) => {
                for shape in shapes {
                    text_shapes(shape, out);
                }
            }
            _ => {}
        }
    }

    fn render_small_palette(
        app: &mut GENtleApp,
        ctx: &egui::Context,
        entries: &[CommandPaletteEntry],
        events: Vec<egui::Event>,
    ) -> (Option<CommandPaletteAction>, Vec<(String, egui::Rect)>) {
        let mut chosen = None;
        let output = ctx.run_ui(small_palette_input(events), |ui| {
            chosen = app.render_command_palette_contents(ui, entries);
        });
        let mut texts = vec![];
        for clipped in &output.shapes {
            let mut shapes = vec![];
            text_shapes(&clipped.shape, &mut shapes);
            texts.extend(
                shapes
                    .into_iter()
                    .filter(|(_, rect)| clipped.clip_rect.intersects(*rect)),
            );
        }
        output.drop_without_applying_deltas();
        (chosen, texts)
    }

    #[test]
    fn palette_detail_shows_purpose_and_recovery_at_minimum_window_size() {
        let ctx = egui::Context::default();
        let mut app = GENtleApp::default();
        let entry = palette_entry(
            "TATA-box Evidence",
            CommandPaletteAction::OpenTataBoxes,
            ActionReadiness::NeedsInput {
                detail:
                    "Select DNA in the project graph, or focus its sequence window, then retry."
                        .into(),
            },
        );
        // Warm egui's layout; this is a headless frame, not a native presentation verdict.
        render_small_palette(&mut app, &ctx, std::slice::from_ref(&entry), vec![]);
        let (chosen, texts) = render_small_palette(
            &mut app,
            &ctx,
            std::slice::from_ref(&entry),
            vec![key_event(Key::Enter)],
        );
        assert!(chosen.is_none());
        let purpose = texts
            .iter()
            .find(|(text, _)| text == &entry.detail)
            .unwrap();
        let recovery = texts
            .iter()
            .find(|(text, _)| text.starts_with("Needs input: Select DNA"))
            .unwrap();
        assert!(purpose.1.bottom() <= 320.0, "{:?}", purpose.1);
        assert!(recovery.1.bottom() <= 320.0, "{:?}", recovery.1);
        assert!(
            texts
                .iter()
                .any(|(text, _)| text.contains("Choose another action"))
        );
    }

    #[test]
    fn disabled_hover_shows_details_and_stationary_pointer_does_not_override_keyboard() {
        let ctx = egui::Context::default();
        let mut app = GENtleApp::default();
        let entries = [
            palette_entry(
                "Open Sequence",
                CommandPaletteAction::OpenSequence,
                ActionReadiness::Ready,
            ),
            palette_entry(
                "TATA-box Evidence",
                CommandPaletteAction::OpenTataBoxes,
                ActionReadiness::NeedsInput {
                    detail: "Select DNA first".into(),
                },
            ),
            palette_entry(
                "Configuration",
                CommandPaletteAction::OpenConfiguration,
                ActionReadiness::Ready,
            ),
        ];
        render_small_palette(&mut app, &ctx, &entries, vec![]);
        let (_, texts) = render_small_palette(&mut app, &ctx, &entries, vec![]);
        let pos = texts
            .iter()
            .find(|(text, _)| text == "TATA-box Evidence")
            .unwrap()
            .1
            .center();
        let (chosen, texts) = render_small_palette(
            &mut app,
            &ctx,
            &entries,
            vec![egui::Event::PointerMoved(pos)],
        );
        assert!(chosen.is_none());
        assert_eq!(app.command_palette_selected, 1);
        assert!(
            texts
                .iter()
                .any(|(text, _)| text == "Needs input: Select DNA first")
        );
        let (chosen, texts) = render_small_palette(
            &mut app,
            &ctx,
            &entries,
            vec![key_event(Key::ArrowDown), key_event(Key::Enter)],
        );
        assert_eq!(app.command_palette_selected, 2);
        assert!(matches!(
            chosen,
            Some(CommandPaletteAction::OpenConfiguration)
        ));
        assert!(texts.iter().any(|(text, _)| text == &entries[2].detail));
    }

    #[test]
    fn palette_details_clear_when_no_entries_match() {
        let ctx = egui::Context::default();
        let mut app = GENtleApp::default();
        app.command_palette_selected = usize::MAX;
        let (chosen, texts) =
            render_small_palette(&mut app, &ctx, &[], vec![key_event(Key::Enter)]);
        assert!(chosen.is_none());
        assert!(texts.iter().any(|(text, _)| text == "No matching commands"));
        assert!(
            !texts
                .iter()
                .any(|(text, _)| text.contains("Enter to choose"))
        );
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
            pos = entry.render_row(ui, true).rect.center();
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
                let response = entry.render_row(ui, true);
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
