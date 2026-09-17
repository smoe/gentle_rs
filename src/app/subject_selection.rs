//! Explicit launch subjects for sequence tools; project inventory is not selection.

use super::*;

pub(super) struct PaletteSubject {
    context: Result<(String, Option<(usize, usize)>), String>,
    viewport: Option<ViewportId>,
    engine: std::sync::Weak<RwLock<GentleEngine>>,
}

impl GENtleApp {
    pub(super) fn capture_palette_subject(&self) -> PaletteSubject {
        PaletteSubject {
            context: self.selected_sequence_context(false),
            viewport: self
                .active_window_menu_key
                .and_then(|key| self.native_window_key_to_viewport.get(&key).copied())
                .filter(|id| self.windows.contains_key(id)),
            engine: Arc::downgrade(&self.engine),
        }
    }

    pub(super) fn selected_sequence_context(
        &self,
        require_dna: bool,
    ) -> Result<(String, Option<(usize, usize)>), String> {
        self.selected_sequence_context_with_origin(require_dna, self.command_palette_dispatching)
    }

    pub(super) fn selected_sequence_context_with_origin(
        &self,
        require_dna: bool,
        palette: bool,
    ) -> Result<(String, Option<(usize, usize)>), String> {
        let active = if palette && self.command_palette_subject.is_some() {
            None
        } else if let Some(window) = self
            .active_window_menu_key
            .and_then(|key| self.native_window_key_to_viewport.get(&key))
            .and_then(|id| self.windows.get(id))
        {
            let guard = window
                .try_read()
                .map_err(|_| "Sequence viewer is busy; try again".to_string())?;
            guard
                .sequence_id()
                .map(|id| (id, guard.selection_range_0based()))
        } else {
            None
        };
        let engine = self
            .engine
            .try_read()
            .map_err(|_| "Project is busy; try again".to_string())?;
        let context = if palette && let Some(subject) = &self.command_palette_subject {
            if !subject.engine.ptr_eq(&Arc::downgrade(&self.engine))
                || subject
                    .viewport
                    .is_some_and(|id| !self.windows.contains_key(&id))
            {
                return Err("The command palette's originating project or sequence window has closed; reopen it from the intended subject".into());
            }
            subject.context.clone()?
        } else if let Some(context) = active {
            context
        } else {
            let root = Self::native_menu_key_for_viewport(ViewportId::ROOT);
            if self.active_window_menu_key.is_some_and(|key| key != root) {
                return Err("Focus the intended sequence window, or select a sequence in the main project graph".into());
            }
            if self.lineage_group_marked_nodes.len() > 1 {
                return Err(
                    "Choose one sequence for this action; multiple graph nodes are marked".into(),
                );
            }
            let seq_id = self.lineage_graph_selected_node_id.as_ref()
                .and_then(|id| engine.state().lineage.nodes.get(id))
                .map(|node| node.seq_id.clone())
                .ok_or_else(|| "Select a sequence in the main project graph, or focus its sequence window, then retry".to_string())?;
            (seq_id, None)
        };
        let kind = engine.sequence_kind(&context.0).ok_or_else(|| {
            format!(
                "Selected sequence '{}' is no longer in the project",
                context.0
            )
        })?;
        if require_dna && kind != "dna" {
            return Err(format!(
                "Selected sequence '{}' is {kind}; select DNA for this action",
                context.0
            ));
        }
        Ok(context)
    }

    pub(super) fn sequence_subject_for_action(
        &mut self,
        action: &str,
        require_dna: bool,
    ) -> Option<String> {
        match self.selected_sequence_context(require_dna) {
            Ok((seq_id, _)) => Some(seq_id),
            Err(reason) => {
                self.app_status = format!("Cannot open {action}: {reason}");
                None
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    // Hand-crafted in-memory fixtures; no reference downloads or assay execution.
    fn mixed_project() -> GENtleApp {
        let mut state = ProjectState::default();
        for (id, molecule) in [
            ("a_protein", "protein"),
            ("b_rna", "RNA"),
            ("z_dna", "dsDNA"),
        ] {
            let mut seq = DNAsequence::from_sequence("ACGTACGT").unwrap();
            seq.set_molecule_type(molecule);
            state.sequences.insert(id.into(), seq);
        }
        GENtleApp {
            engine: Arc::new(RwLock::new(GentleEngine::from_state(state))),
            ..Default::default()
        }
    }

    fn focus(app: &mut GENtleApp, viewport: ViewportId) {
        let key = GENtleApp::native_menu_key_for_viewport(viewport);
        // Model the native window registry without publishing global test focus.
        app.native_window_key_to_viewport.insert(key, viewport);
        app.active_window_menu_key = Some(key);
    }

    fn select(app: &mut GENtleApp, seq_id: &str) {
        app.lineage_graph_selected_node_id = app
            .engine
            .read()
            .unwrap()
            .state()
            .lineage
            .seq_to_node
            .get(seq_id)
            .cloned();
        focus(app, ViewportId::ROOT);
    }

    #[test]
    fn inventory_is_not_selection_and_kind_is_checked() {
        let mut app = mixed_project();
        assert!(app.selected_sequence_context(true).is_err());
        select(&mut app, "a_protein");
        assert!(
            app.selected_sequence_context(true)
                .unwrap_err()
                .contains("protein")
        );
        assert_eq!(app.selected_sequence_context(false).unwrap().0, "a_protein");
        select(&mut app, "b_rna");
        assert!(
            app.selected_sequence_context(true)
                .unwrap_err()
                .contains("rna")
        );
        select(&mut app, "z_dna");
        assert_eq!(app.selected_sequence_context(true).unwrap().0, "z_dna");
        app.engine
            .write()
            .unwrap()
            .state_mut()
            .sequences
            .remove("z_dna");
        assert!(
            app.selected_sequence_context(true)
                .unwrap_err()
                .contains("no longer")
        );
    }

    #[test]
    fn dna_launchers_do_not_fall_back_to_first_sequence() {
        let actions: [fn(&mut GENtleApp); 7] = [
            GENtleApp::open_cryptic_splicing_screen,
            GENtleApp::open_tata_box_workspace,
            GENtleApp::open_tss_inventory_workspace,
            GENtleApp::open_saved_genomic_regions,
            GENtleApp::open_genomic_region_conservation,
            GENtleApp::open_pcr_design_dialog,
            GENtleApp::open_sequencing_confirmation_dialog,
        ];
        for action in actions {
            let mut app = mixed_project();
            action(&mut app);
            assert!(app.new_windows.is_empty());
            assert!(!app.show_pcr_design_dialog && !app.show_sequencing_confirmation_dialog);
            assert!(
                app.app_status.contains("Select a sequence"),
                "{}",
                app.app_status
            );
            select(&mut app, "a_protein");
            action(&mut app);
            assert!(app.app_status.contains("protein"), "{}", app.app_status);
            assert!(app.new_windows.is_empty());
            app.engine
                .write()
                .unwrap()
                .state_mut()
                .sequences
                .retain(|id, _| id == "a_protein");
            action(&mut app);
            assert!(
                app.new_windows.is_empty(),
                "protein-only project must not launch a DNA viewer"
            );
        }
    }

    #[test]
    fn explicit_pcr_binding_checks_kind_and_opens_selected_dna() {
        let mut app = mixed_project();
        assert!(
            app.open_pcr_design_dialog_for_seq_id("a_protein")
                .unwrap_err()
                .contains("requires DNA")
        );
        assert!(app.open_pcr_design_dialog_for_seq_id("b_rna").is_err());
        assert!(app.new_windows.is_empty());
        select(&mut app, "z_dna");
        app.open_pcr_design_dialog();
        assert!(app.show_pcr_design_dialog);
        assert_eq!(app.pcr_design_seq_id, "z_dna");
        assert_eq!(app.new_windows.len(), 1);
    }

    #[test]
    fn gibson_prefills_only_explicit_destination_not_arbitrary_insert() {
        let mut app = mixed_project();
        app.prefill_gibson_from_active_context();
        assert!(app.gibson_destination_seq_id.is_empty());
        select(&mut app, "z_dna");
        app.prefill_gibson_from_active_context();
        assert_eq!(app.gibson_destination_seq_id, "z_dna");
        assert!(app.gibson_insert_seq_id.is_empty());
    }

    #[test]
    fn root_closed_viewer_and_multiple_marked_subjects_do_not_borrow_context() {
        let mut app = mixed_project();
        let viewport =
            app.register_window(Window::new_dna_lazy("z_dna".into(), app.engine.clone()));
        focus(&mut app, viewport);
        assert_eq!(app.selected_sequence_context(true).unwrap().0, "z_dna");
        focus(&mut app, ViewportId::ROOT);
        assert!(app.selected_sequence_context(true).is_err());
        select(&mut app, "z_dna");
        app.lineage_group_marked_nodes
            .extend(["one".into(), "two".into()]);
        assert!(
            app.selected_sequence_context(true)
                .unwrap_err()
                .contains("multiple")
        );
        app.lineage_group_marked_nodes.clear();
        focus(&mut app, viewport);
        app.windows.remove(&viewport);
        assert!(app.selected_sequence_context(true).is_err());
    }

    #[test]
    fn specialist_paint_does_not_silently_fill_empty_targets() {
        let mut app = mixed_project();
        let ctx = egui::Context::default();
        ctx.run_ui(egui::RawInput::default(), |ui| {
            app.render_pcr_design_contents(ui, &ctx);
            app.render_sequencing_confirmation_contents(ui, &ctx);
        })
        .drop_without_applying_deltas();
        assert!(app.pcr_design_seq_id.is_empty());
        assert!(app.sequencing_confirmation_seq_id.is_empty());
        assert!(app.new_windows.is_empty());
    }

    #[test]
    fn palette_keeps_explicit_origin_across_embedded_and_separate_focus() {
        for palette_focus in [ViewportId::ROOT, GENtleApp::command_palette_viewport_id()] {
            let mut app = mixed_project();
            let viewport =
                app.register_window(Window::new_dna_lazy("z_dna".into(), app.engine.clone()));
            focus(&mut app, viewport);
            app.open_command_palette_dialog();
            select(&mut app, "a_protein");
            focus(&mut app, palette_focus);
            app.execute_command_palette_action(
                &egui::Context::default(),
                CommandPaletteAction::OpenTataBoxes,
            );
            assert!(
                app.windows[&viewport]
                    .read()
                    .unwrap()
                    .tata_workspace_open_or_pending()
            );
            assert!(app.app_status.contains("z_dna"), "{}", app.app_status);
            assert!(!app.command_palette_dispatching);
        }
    }

    #[test]
    fn palette_rejects_closed_origin_and_project_replacement() {
        for replace_project in [false, true] {
            let mut app = mixed_project();
            let viewport =
                app.register_window(Window::new_dna_lazy("z_dna".into(), app.engine.clone()));
            focus(&mut app, viewport);
            app.open_command_palette_dialog();
            if replace_project {
                app.engine = mixed_project().engine;
            } else {
                app.windows.remove(&viewport);
            }
            focus(&mut app, GENtleApp::command_palette_viewport_id());
            app.execute_command_palette_action(
                &egui::Context::default(),
                CommandPaletteAction::OpenTataBoxes,
            );
            assert!(app.app_status.contains("has closed"), "{}", app.app_status);
            assert!(app.new_windows.is_empty());
        }
    }
}
