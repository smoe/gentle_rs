//! Root workspace orchestration helpers.
//!
//! This module starts the gradual extraction of `render_root_ui` by moving
//! keyboard shortcut dispatch out of the monolithic update body without
//! changing render order or project/window ownership.

use super::*;

impl GENtleApp {
    pub(super) fn handle_root_keyboard_shortcuts(&mut self, ctx: &egui::Context) {
        let open_project = KeyboardShortcut::new(Modifiers::COMMAND, Key::O);
        let new_project = KeyboardShortcut::new(Modifiers::COMMAND, Key::N);
        let open_sequence = KeyboardShortcut::new(Modifiers::COMMAND | Modifiers::SHIFT, Key::O);
        let open_retrieve_genome =
            KeyboardShortcut::new(Modifiers::COMMAND | Modifiers::SHIFT, Key::G);
        let open_prepare_genome =
            KeyboardShortcut::new(Modifiers::COMMAND | Modifiers::SHIFT, Key::P);
        let open_blast_genome =
            KeyboardShortcut::new(Modifiers::COMMAND | Modifiers::SHIFT, Key::L);
        let open_import_bed_track =
            KeyboardShortcut::new(Modifiers::COMMAND | Modifiers::SHIFT, Key::B);
        let open_agent_assistant =
            KeyboardShortcut::new(Modifiers::COMMAND | Modifiers::SHIFT, Key::A);
        let save_project = KeyboardShortcut::new(Modifiers::COMMAND, Key::S);
        let close_project = KeyboardShortcut::new(Modifiers::COMMAND | Modifiers::SHIFT, Key::W);
        let quit_application = KeyboardShortcut::new(Modifiers::COMMAND, Key::Q);
        let open_configuration = KeyboardShortcut::new(Modifiers::COMMAND, Key::Comma);
        let focus_main_window = KeyboardShortcut::new(Modifiers::COMMAND, Key::Backtick);
        let open_command_palette = KeyboardShortcut::new(Modifiers::COMMAND, Key::K);
        let open_help_f1 = KeyboardShortcut::new(Modifiers::NONE, Key::F1);
        let open_help_ctrl_f1 = KeyboardShortcut::new(Modifiers::CTRL, Key::F1);
        let open_help_cmd_shift_slash =
            KeyboardShortcut::new(Modifiers::COMMAND | Modifiers::SHIFT, Key::Slash);
        let undo_shortcut = KeyboardShortcut::new(Modifiers::COMMAND, Key::Z);
        let redo_shortcut_shift =
            KeyboardShortcut::new(Modifiers::COMMAND | Modifiers::SHIFT, Key::Z);
        let redo_shortcut_y = KeyboardShortcut::new(Modifiers::COMMAND, Key::Y);

        if ctx.input_mut(|i| i.consume_shortcut(&new_project)) {
            self.request_project_action(ProjectAction::New);
        }
        if ctx.input_mut(|i| i.consume_shortcut(&open_project)) {
            self.request_project_action(ProjectAction::Open);
        }
        if ctx.input_mut(|i| i.consume_shortcut(&open_sequence)) {
            self.prompt_open_sequence();
        }
        if ctx.input_mut(|i| i.consume_shortcut(&open_retrieve_genome)) {
            self.open_reference_genome_retrieve_dialog();
        }
        if ctx.input_mut(|i| i.consume_shortcut(&open_prepare_genome)) {
            self.open_reference_genome_prepare_dialog();
        }
        if ctx.input_mut(|i| i.consume_shortcut(&open_blast_genome)) {
            self.open_reference_genome_blast_dialog();
        }
        if ctx.input_mut(|i| i.consume_shortcut(&open_import_bed_track)) {
            self.open_genome_bed_track_dialog();
        }
        if ctx.input_mut(|i| i.consume_shortcut(&open_agent_assistant)) {
            self.open_agent_assistant_dialog();
        }
        if ctx.input_mut(|i| i.consume_shortcut(&save_project)) {
            let _ = self.save_current_project();
        }
        if ctx.input_mut(|i| i.consume_shortcut(&close_project)) {
            self.request_project_action(ProjectAction::Close);
        }
        if ctx.input_mut(|i| i.consume_shortcut(&quit_application))
            || Self::consume_command_or_ctrl_shortcut(ctx, Key::Q)
        {
            self.request_project_action(ProjectAction::Quit);
        }
        if ctx.input_mut(|i| i.consume_shortcut(&open_configuration)) {
            self.open_configuration_dialog();
        }
        if ctx.input_mut(|i| i.consume_shortcut(&focus_main_window)) {
            self.queue_focus_viewport(ViewportId::ROOT);
        }
        if ctx.input_mut(|i| i.consume_shortcut(&open_command_palette)) {
            self.open_command_palette_dialog();
        }
        if ctx.input_mut(|i| i.consume_shortcut(&open_help_f1))
            || ctx.input_mut(|i| i.consume_shortcut(&open_help_ctrl_f1))
            || ctx.input_mut(|i| i.consume_shortcut(&open_help_cmd_shift_slash))
        {
            self.open_help_doc(HelpDoc::Gui);
        }
        if ctx.input_mut(|i| i.consume_shortcut(&undo_shortcut)) {
            self.undo_last_operation();
        }
        if ctx.input_mut(|i| i.consume_shortcut(&redo_shortcut_shift))
            || ctx.input_mut(|i| i.consume_shortcut(&redo_shortcut_y))
        {
            self.redo_last_operation();
        }
    }
}
