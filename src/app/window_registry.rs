//! Open-window registry and deferred child-window placement helpers.
//!
//! This submodule keeps the application-wide window listing/focus plumbing in
//! one place so the main app file can focus more on menu/dialog orchestration.

use super::*;
use crate::engine_shell::UiIntentTarget;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct GuiProminentGlossaryEntry {
    pub glossary_path: &'static str,
    pub menu_path: &'static str,
    pub palette_title: &'static str,
    pub ui_intent_target: Option<UiIntentTarget>,
}

const GUI_PROMINENT_GLOSSARY_ENTRIES: &[GuiProminentGlossaryEntry] = &[
    GuiProminentGlossaryEntry {
        glossary_path: "load-project",
        menu_path: "File > Open Project...",
        palette_title: "Open Project",
        ui_intent_target: None,
    },
    GuiProminentGlossaryEntry {
        glossary_path: "save-project",
        menu_path: "File > Save Project...",
        palette_title: "Save Project",
        ui_intent_target: None,
    },
    GuiProminentGlossaryEntry {
        glossary_path: "history undo",
        menu_path: "Edit > Undo",
        palette_title: "Undo Last Operation",
        ui_intent_target: None,
    },
    GuiProminentGlossaryEntry {
        glossary_path: "history redo",
        menu_path: "Edit > Redo",
        palette_title: "Redo Last Operation",
        ui_intent_target: None,
    },
    GuiProminentGlossaryEntry {
        glossary_path: "gibson preview",
        menu_path: "Patterns > Gibson...",
        palette_title: "Gibson",
        ui_intent_target: None,
    },
    GuiProminentGlossaryEntry {
        glossary_path: "gibson apply",
        menu_path: "Patterns > Gibson...",
        palette_title: "Gibson",
        ui_intent_target: None,
    },
    GuiProminentGlossaryEntry {
        glossary_path: "genomes prepare",
        menu_path: "Genome > Prepare Reference Genome...",
        palette_title: "Prepare Reference Genome",
        ui_intent_target: Some(UiIntentTarget::PrepareReferenceGenome),
    },
    GuiProminentGlossaryEntry {
        glossary_path: "genomes extract-region",
        menu_path: "Genome > Retrieve Genomic Sequence...",
        palette_title: "Retrieve Genomic Sequence",
        ui_intent_target: Some(UiIntentTarget::RetrieveGenomeSequence),
    },
    GuiProminentGlossaryEntry {
        glossary_path: "genomes extract-gene",
        menu_path: "Genome > Retrieve Genomic Sequence...",
        palette_title: "Retrieve Genomic Sequence",
        ui_intent_target: Some(UiIntentTarget::RetrieveGenomeSequence),
    },
    GuiProminentGlossaryEntry {
        glossary_path: "genomes blast-track",
        menu_path: "Genome > BLAST Genome Sequence...",
        palette_title: "BLAST Genome Sequence",
        ui_intent_target: Some(UiIntentTarget::BlastGenomeSequence),
    },
    GuiProminentGlossaryEntry {
        glossary_path: "helpers prepare",
        menu_path: "Genome > Prepare Helper Genome...",
        palette_title: "Prepare Helper Genome",
        ui_intent_target: Some(UiIntentTarget::PrepareHelperGenome),
    },
    GuiProminentGlossaryEntry {
        glossary_path: "helpers extract-region",
        menu_path: "Genome > Retrieve Helper Sequence...",
        palette_title: "Retrieve Helper Sequence",
        ui_intent_target: Some(UiIntentTarget::RetrieveHelperSequence),
    },
    GuiProminentGlossaryEntry {
        glossary_path: "helpers extract-gene",
        menu_path: "Genome > Retrieve Helper Sequence...",
        palette_title: "Retrieve Helper Sequence",
        ui_intent_target: Some(UiIntentTarget::RetrieveHelperSequence),
    },
    GuiProminentGlossaryEntry {
        glossary_path: "helpers blast-track",
        menu_path: "Genome > BLAST Helper Sequence...",
        palette_title: "BLAST Helper Sequence",
        ui_intent_target: Some(UiIntentTarget::BlastHelperSequence),
    },
    GuiProminentGlossaryEntry {
        glossary_path: "tracks import-bed",
        menu_path: "Genome > Import Genome Track...",
        palette_title: "Import Genome Track",
        ui_intent_target: Some(UiIntentTarget::ImportGenomeTrack),
    },
    GuiProminentGlossaryEntry {
        glossary_path: "tracks import-bigwig",
        menu_path: "Genome > Import Genome Track...",
        palette_title: "Import Genome Track",
        ui_intent_target: Some(UiIntentTarget::ImportGenomeTrack),
    },
    GuiProminentGlossaryEntry {
        glossary_path: "tracks import-vcf",
        menu_path: "Genome > Import Genome Track...",
        palette_title: "Import Genome Track",
        ui_intent_target: Some(UiIntentTarget::ImportGenomeTrack),
    },
    GuiProminentGlossaryEntry {
        glossary_path: "primers design",
        menu_path: "Patterns > PCR Designer...",
        palette_title: "PCR Designer",
        ui_intent_target: Some(UiIntentTarget::PcrDesign),
    },
    GuiProminentGlossaryEntry {
        glossary_path: "primers design-qpcr",
        menu_path: "Patterns > PCR Designer...",
        palette_title: "PCR Designer",
        ui_intent_target: Some(UiIntentTarget::PcrDesign),
    },
    GuiProminentGlossaryEntry {
        glossary_path: "seq-confirm run",
        menu_path: "Patterns > Sequencing Confirmation...",
        palette_title: "Sequencing Confirmation",
        ui_intent_target: Some(UiIntentTarget::SequencingConfirmation),
    },
    GuiProminentGlossaryEntry {
        glossary_path: "mirna scan-target",
        menu_path: "Patterns > microRNA Target Scan...",
        palette_title: "microRNA Target Scan",
        ui_intent_target: None,
    },
];

pub fn gui_prominent_glossary_entries() -> &'static [GuiProminentGlossaryEntry] {
    GUI_PROMINENT_GLOSSARY_ENTRIES
}

#[derive(Clone)]
pub(super) struct OpenWindowEntry {
    pub(super) native_menu_key: u64,
    pub(super) viewport_id: ViewportId,
    pub(super) title: String,
    pub(super) detail: String,
}

impl GENtleApp {
    fn embedded_window_layer_from_window_id(window_id: egui::Id) -> egui::LayerId {
        egui::LayerId::new(egui::Order::Middle, window_id)
    }

    fn embedded_window_layer_from_title(title: impl Into<String>) -> egui::LayerId {
        Self::embedded_window_layer_from_window_id(egui::Id::new(title.into()))
    }

    pub(crate) fn embedded_window_layer_id_for_viewport(
        &self,
        viewport_id: ViewportId,
    ) -> Option<egui::LayerId> {
        if viewport_id == ViewportId::ROOT {
            return Some(Self::embedded_window_layer_from_window_id(
                Self::main_workspace_hosted_window_id(),
            ));
        }
        if viewport_id == Self::configuration_viewport_id() {
            return Some(Self::embedded_window_layer_from_window_id(
                Self::hosted_configuration_window_id(),
            ));
        }
        if viewport_id == Self::help_viewport_id() {
            return Some(Self::embedded_window_layer_from_window_id(
                Self::hosted_help_window_id(),
            ));
        }
        if viewport_id == Self::command_palette_viewport_id() {
            return Some(Self::embedded_window_layer_from_title("Command Palette"));
        }
        if viewport_id == Self::history_viewport_id() {
            return Some(Self::embedded_window_layer_from_title("Operation History"));
        }
        if viewport_id == Self::prepare_genome_viewport_id() {
            return Some(Self::embedded_window_layer_from_window_id(egui::Id::new((
                "hosted_prepare_genome_window",
                viewport_id,
            ))));
        }
        if viewport_id == Self::retrieve_genome_viewport_id() {
            return Some(Self::embedded_window_layer_from_window_id(egui::Id::new((
                "hosted_retrieve_genome_window",
                viewport_id,
            ))));
        }
        if viewport_id == Self::blast_genome_viewport_id() {
            return Some(Self::embedded_window_layer_from_window_id(egui::Id::new((
                "hosted_blast_genome_window",
                viewport_id,
            ))));
        }
        if viewport_id == Self::bed_track_viewport_id() {
            return Some(Self::embedded_window_layer_from_window_id(egui::Id::new((
                "hosted_bed_track_window",
                viewport_id,
            ))));
        }
        if viewport_id == Self::gibson_viewport_id() {
            return Some(Self::embedded_window_layer_from_window_id(egui::Id::new((
                "hosted_gibson_window",
                viewport_id,
            ))));
        }
        if viewport_id == Self::arrangement_gel_preview_viewport_id() {
            return Some(Self::embedded_window_layer_from_window_id(egui::Id::new((
                "hosted_arrangement_gel_preview_window",
                viewport_id,
            ))));
        }
        if viewport_id == Self::pcr_design_viewport_id() {
            return Some(Self::embedded_window_layer_from_window_id(egui::Id::new((
                "hosted_pcr_design_window",
                viewport_id,
            ))));
        }
        if viewport_id == Self::sequencing_confirmation_viewport_id() {
            return Some(Self::embedded_window_layer_from_window_id(egui::Id::new((
                "hosted_sequencing_confirmation_window",
                viewport_id,
            ))));
        }
        if viewport_id == Self::planning_viewport_id() {
            return Some(Self::embedded_window_layer_from_window_id(egui::Id::new((
                "hosted_planning_window",
                viewport_id,
            ))));
        }
        if viewport_id == Self::routine_assistant_viewport_id() {
            return Some(Self::embedded_window_layer_from_window_id(
                Self::hosted_routine_assistant_window_id(),
            ));
        }
        if viewport_id == Self::agent_assistant_viewport_id() {
            return Some(Self::embedded_window_layer_from_window_id(
                Self::hosted_agent_assistant_window_id(),
            ));
        }
        if viewport_id == Self::uniprot_viewport_id() {
            return Some(Self::embedded_window_layer_from_window_id(egui::Id::new((
                "hosted_uniprot_window",
                viewport_id,
            ))));
        }
        if viewport_id == Self::genbank_viewport_id() {
            return Some(Self::embedded_window_layer_from_title(
                "GenBank / dbSNP Fetch",
            ));
        }
        if self.windows.contains_key(&viewport_id) {
            return Some(Self::embedded_window_layer_from_window_id(egui::Id::new((
                "hosted_sequence_window_v2",
                viewport_id,
            ))));
        }
        self.windows
            .values()
            .find_map(|window| {
                window
                    .read()
                    .ok()
                    .and_then(|guard| guard.embedded_auxiliary_window_layer_id(viewport_id))
            })
            .or_else(|| {
                self.detached_auxiliary_window_hosts
                    .values()
                    .find_map(|window| {
                        window
                            .read()
                            .ok()
                            .and_then(|guard| guard.embedded_auxiliary_window_layer_id(viewport_id))
                    })
            })
    }

    pub(super) fn collect_open_window_entries(&self) -> Vec<OpenWindowEntry> {
        let mut entries = vec![OpenWindowEntry {
            native_menu_key: Self::native_menu_key_for_viewport(ViewportId::ROOT),
            viewport_id: ViewportId::ROOT,
            title: format!("Main Window — {}", self.current_project_name()),
            detail: "Project workspace".to_string(),
        }];

        if self.show_configuration_dialog {
            entries.push(OpenWindowEntry {
                native_menu_key: Self::native_menu_key_for_viewport(
                    Self::configuration_viewport_id(),
                ),
                viewport_id: Self::configuration_viewport_id(),
                title: "Configuration".to_string(),
                detail: "External tools and graphics defaults".to_string(),
            });
        }
        if self.show_help_dialog {
            entries.push(OpenWindowEntry {
                native_menu_key: Self::native_menu_key_for_viewport(Self::help_viewport_id()),
                viewport_id: Self::help_viewport_id(),
                title: format!("Help — {}", self.active_help_title()),
                detail: "GUI/CLI manual".to_string(),
            });
        }
        if self.show_command_palette_dialog {
            entries.push(OpenWindowEntry {
                native_menu_key: Self::native_menu_key_for_viewport(
                    Self::command_palette_viewport_id(),
                ),
                viewport_id: Self::command_palette_viewport_id(),
                title: "Command Palette".to_string(),
                detail: "Action launcher".to_string(),
            });
        }
        if self.external_services_ui.show_panel {
            entries.push(OpenWindowEntry {
                native_menu_key: Self::native_menu_key_for_viewport(
                    Self::external_services_viewport_id(),
                ),
                viewport_id: Self::external_services_viewport_id(),
                title: "External Services".to_string(),
                detail: "Provider catalog, preflight, and quote handoff workspace".to_string(),
            });
        }
        if self.evidence_preparation_panel.show_panel {
            entries.push(OpenWindowEntry {
                native_menu_key: Self::native_menu_key_for_viewport(
                    Self::evidence_preparation_viewport_id(),
                ),
                viewport_id: Self::evidence_preparation_viewport_id(),
                title: "Evidence Preparation".to_string(),
                detail: "TP73 evidence-viewer proof material preparation".to_string(),
            });
        }
        if self.history_ui.show_panel {
            entries.push(OpenWindowEntry {
                native_menu_key: Self::native_menu_key_for_viewport(Self::history_viewport_id()),
                viewport_id: Self::history_viewport_id(),
                title: "Operation History".to_string(),
                detail: "Undo/redo operation log".to_string(),
            });
        }
        if self.show_jobs_panel {
            entries.push(OpenWindowEntry {
                native_menu_key: Self::native_menu_key_for_viewport(
                    Self::background_jobs_viewport_id(),
                ),
                viewport_id: Self::background_jobs_viewport_id(),
                title: "Background Jobs".to_string(),
                detail: "Progress, cancellation, and retry snapshots".to_string(),
            });
        }
        if self.show_reference_genome_prepare_dialog {
            let scope = self.genome_dialog_scope;
            entries.push(OpenWindowEntry {
                native_menu_key: Self::native_menu_key_for_viewport(
                    Self::prepare_genome_viewport_id(),
                ),
                viewport_id: Self::prepare_genome_viewport_id(),
                title: scope.prepare_title().to_string(),
                detail: format!("{} preparation", scope.description()),
            });
        }
        if self.show_reference_genome_retrieve_dialog {
            let scope = self.genome_dialog_scope;
            entries.push(OpenWindowEntry {
                native_menu_key: Self::native_menu_key_for_viewport(
                    Self::retrieve_genome_viewport_id(),
                ),
                viewport_id: Self::retrieve_genome_viewport_id(),
                title: scope.retrieve_title().to_string(),
                detail: format!(
                    "Gene/region extraction from prepared {} indices",
                    scope.description()
                ),
            });
        }
        if self.show_reference_genome_blast_dialog {
            let scope = self.genome_dialog_scope;
            entries.push(OpenWindowEntry {
                native_menu_key: Self::native_menu_key_for_viewport(
                    Self::blast_genome_viewport_id(),
                ),
                viewport_id: Self::blast_genome_viewport_id(),
                title: scope.blast_title().to_string(),
                detail: format!("BLAST against prepared {} indices", scope.description()),
            });
        }
        if self.show_genome_bed_track_dialog {
            entries.push(OpenWindowEntry {
                native_menu_key: Self::native_menu_key_for_viewport(Self::bed_track_viewport_id()),
                viewport_id: Self::bed_track_viewport_id(),
                title: "Import Genome Tracks".to_string(),
                detail: "Import BED/BigWig/VCF track overlays".to_string(),
            });
        }
        if self.show_gibson_dialog {
            entries.push(OpenWindowEntry {
                native_menu_key: Self::native_menu_key_for_viewport(Self::gibson_viewport_id()),
                viewport_id: Self::gibson_viewport_id(),
                title: "Gibson".to_string(),
                detail: "Destination-first Gibson planning with primer and cartoon preview"
                    .to_string(),
            });
        }
        if self.show_arrangement_gel_preview_dialog {
            entries.push(OpenWindowEntry {
                native_menu_key: Self::native_menu_key_for_viewport(
                    Self::arrangement_gel_preview_viewport_id(),
                ),
                viewport_id: Self::arrangement_gel_preview_viewport_id(),
                title: self.arrangement_gel_preview_title(),
                detail: "Serial gel preview with live ladder selection".to_string(),
            });
        }
        if self.show_pcr_design_dialog {
            entries.push(OpenWindowEntry {
                native_menu_key: Self::native_menu_key_for_viewport(Self::pcr_design_viewport_id()),
                viewport_id: Self::pcr_design_viewport_id(),
                title: if self.pcr_design_seq_id.trim().is_empty() {
                    "PCR Designer".to_string()
                } else {
                    format!("PCR Designer — {}", self.pcr_design_seq_id.trim())
                },
                detail: "Paint-first pair-PCR specialist".to_string(),
            });
        }
        if self.show_sequencing_confirmation_dialog {
            entries.push(OpenWindowEntry {
                native_menu_key: Self::native_menu_key_for_viewport(
                    Self::sequencing_confirmation_viewport_id(),
                ),
                viewport_id: Self::sequencing_confirmation_viewport_id(),
                title: if self.sequencing_confirmation_seq_id.trim().is_empty() {
                    "Sequencing Confirmation".to_string()
                } else {
                    format!(
                        "Sequencing Confirmation — {}",
                        self.sequencing_confirmation_seq_id.trim()
                    )
                },
                detail: "Called-read construct confirmation specialist".to_string(),
            });
        }
        if self.show_planning_dialog {
            entries.push(OpenWindowEntry {
                native_menu_key: Self::native_menu_key_for_viewport(Self::planning_viewport_id()),
                viewport_id: Self::planning_viewport_id(),
                title: "Planning".to_string(),
                detail: "Planning profiles/objectives and sync suggestions".to_string(),
            });
        }
        if self.show_routine_assistant_dialog {
            entries.push(OpenWindowEntry {
                native_menu_key: Self::native_menu_key_for_viewport(
                    Self::routine_assistant_viewport_id(),
                ),
                viewport_id: Self::routine_assistant_viewport_id(),
                title: "Routine Assistant".to_string(),
                detail: "Staged routine selection, preflight, and execution".to_string(),
            });
        }
        if self.show_agent_assistant_dialog {
            entries.push(OpenWindowEntry {
                native_menu_key: Self::native_menu_key_for_viewport(
                    Self::agent_assistant_viewport_id(),
                ),
                viewport_id: Self::agent_assistant_viewport_id(),
                title: "Agent Assistant".to_string(),
                detail: "Agent chat and per-reply command execution".to_string(),
            });
        }
        if self.show_jaspar_expert_dialog {
            entries.push(OpenWindowEntry {
                native_menu_key: Self::native_menu_key_for_viewport(
                    Self::jaspar_expert_viewport_id(),
                ),
                viewport_id: Self::jaspar_expert_viewport_id(),
                title: "JASPAR Expert".to_string(),
                detail: "Motif registry, sequence-logo, and scoring diagnostics".to_string(),
            });
        }
        if self.show_new_sequence_dialog {
            entries.push(OpenWindowEntry {
                native_menu_key: Self::native_menu_key_for_viewport(
                    Self::new_sequence_viewport_id(),
                ),
                viewport_id: Self::new_sequence_viewport_id(),
                title: "New Sequence".to_string(),
                detail: "Typed/clipboard sequence creation".to_string(),
            });
        }
        if self.show_uniprot_dialog {
            entries.push(OpenWindowEntry {
                native_menu_key: Self::native_menu_key_for_viewport(Self::uniprot_viewport_id()),
                viewport_id: Self::uniprot_viewport_id(),
                title: "Protein Evidence".to_string(),
                detail: "UniProt/Ensembl protein evidence fetch, import, and compare tools"
                    .to_string(),
            });
        }
        if self.show_genbank_dialog {
            entries.push(OpenWindowEntry {
                native_menu_key: Self::native_menu_key_for_viewport(Self::genbank_viewport_id()),
                viewport_id: Self::genbank_viewport_id(),
                title: "GenBank / dbSNP Fetch".to_string(),
                detail: "GenBank accession and dbSNP retrieval tool".to_string(),
            });
        }

        let mut sequence_windows = self
            .windows
            .iter()
            .map(|(viewport_id, window)| {
                let title = window
                    .read()
                    .map(|w| w.name())
                    .unwrap_or_else(|_| "Sequence window".to_string());
                OpenWindowEntry {
                    native_menu_key: Self::native_menu_key_for_viewport(*viewport_id),
                    viewport_id: *viewport_id,
                    title,
                    detail: "Sequence map window".to_string(),
                }
            })
            .collect::<Vec<_>>();
        sequence_windows.sort_by(|left, right| left.title.cmp(&right.title));
        entries.extend(sequence_windows);

        let mut auxiliary_windows = Vec::new();
        for window in self.windows.values() {
            let Ok(window_guard) = window.read() else {
                continue;
            };
            for (viewport_id, title, detail) in window_guard.collect_open_auxiliary_window_entries()
            {
                auxiliary_windows.push(OpenWindowEntry {
                    native_menu_key: Self::native_menu_key_for_viewport(viewport_id),
                    viewport_id,
                    title,
                    detail,
                });
            }
        }
        for window in self.detached_auxiliary_window_hosts.values() {
            let Ok(window_guard) = window.read() else {
                continue;
            };
            for (viewport_id, title, detail) in window_guard.collect_open_auxiliary_window_entries()
            {
                auxiliary_windows.push(OpenWindowEntry {
                    native_menu_key: Self::native_menu_key_for_viewport(viewport_id),
                    viewport_id,
                    title,
                    detail,
                });
            }
        }
        auxiliary_windows.sort_by(|left, right| {
            left.title
                .cmp(&right.title)
                .then(left.detail.cmp(&right.detail))
        });
        entries.extend(auxiliary_windows);
        entries
    }

    fn collect_open_auxiliary_window_titles(&self) -> Vec<String> {
        let mut titles = Vec::new();
        let mut seen = HashSet::new();
        for window_map in [&self.windows, &self.detached_auxiliary_window_hosts] {
            for window in window_map.values() {
                let Ok(window_guard) = window.read() else {
                    continue;
                };
                for (_, title, _) in window_guard.collect_open_auxiliary_window_entries() {
                    if seen.insert(title.clone()) {
                        titles.push(title);
                    }
                }
            }
        }
        titles
    }

    pub(super) fn reset_root_auxiliary_areas_if_legacy_title_layers_visible(
        &self,
        ctx: &egui::Context,
    ) -> bool {
        let has_stale_auxiliary_title_layer = self
            .collect_open_auxiliary_window_titles()
            .iter()
            .any(|title| {
                let layer_id = Self::stale_hosted_window_title_layer_id(title);
                ctx.memory(|mem| mem.areas().is_visible(&layer_id))
            });
        if has_stale_auxiliary_title_layer {
            ctx.memory_mut(|mem| mem.reset_areas());
            ctx.request_repaint();
            true
        } else {
            false
        }
    }

    pub(super) fn process_window_close_queue(&mut self) {
        let Ok(mut to_close) = self.windows_to_close.write() else {
            eprintln!("W GENtleApp: close-queue lock poisoned");
            return;
        };
        let queued_ids: Vec<ViewportId> = to_close.drain(..).collect();
        drop(to_close);

        for viewport_id in queued_ids {
            self.pending_focus_viewports.retain(|id| *id != viewport_id);
            self.pending_window_initial_positions.remove(&viewport_id);
            self.pending_window_open_timestamps.remove(&viewport_id);
            self.pending_viewport_focus_timestamps.remove(&viewport_id);

            let Some(window) = self.windows.remove(&viewport_id) else {
                self.detached_auxiliary_window_hosts.remove(&viewport_id);
                continue;
            };
            let keep_auxiliary_host_alive = window
                .read()
                .ok()
                .map(|guard| guard.has_open_auxiliary_windows())
                .unwrap_or(false);
            if keep_auxiliary_host_alive {
                self.detached_auxiliary_window_hosts
                    .insert(viewport_id, window);
            }
        }
    }

    pub(super) fn render_detached_auxiliary_window_hosts(&mut self, ctx: &egui::Context) {
        let visible_auxiliary_viewports = self
            .windows
            .values()
            .flat_map(|window| {
                window
                    .read()
                    .ok()
                    .map(|guard| {
                        guard
                            .collect_open_auxiliary_window_entries()
                            .into_iter()
                            .map(|(viewport_id, _, _)| viewport_id)
                            .collect::<Vec<_>>()
                    })
                    .unwrap_or_default()
            })
            .collect::<HashSet<_>>();
        let detached_hosts = self
            .detached_auxiliary_window_hosts
            .iter()
            .map(|(viewport_id, window)| (*viewport_id, window.clone()))
            .collect::<Vec<_>>();
        let mut stale_hosts = Vec::new();
        for (viewport_id, window) in detached_hosts {
            let detached_auxiliary_viewports = window
                .read()
                .ok()
                .map(|guard| {
                    guard
                        .collect_open_auxiliary_window_entries()
                        .into_iter()
                        .map(|(aux_viewport_id, _, _)| aux_viewport_id)
                        .collect::<Vec<_>>()
                })
                .unwrap_or_default();
            if detached_auxiliary_viewports
                .iter()
                .any(|aux_viewport_id| visible_auxiliary_viewports.contains(aux_viewport_id))
            {
                stale_hosts.push(viewport_id);
                continue;
            }
            let update_result = catch_unwind(AssertUnwindSafe(|| {
                if let Ok(mut guard) = window.write() {
                    guard.update_auxiliary_windows_only(ctx);
                } else {
                    eprintln!("W GENtleApp: detached auxiliary host lock poisoned");
                }
            }));
            if update_result.is_err() {
                eprintln!(
                    "E GENtleApp: recovered from panic while updating detached auxiliary host"
                );
            }
            let keep_host = window
                .read()
                .ok()
                .map(|guard| guard.has_open_auxiliary_windows())
                .unwrap_or(false);
            if !keep_host {
                stale_hosts.push(viewport_id);
            }
        }
        for viewport_id in stale_hosts {
            self.detached_auxiliary_window_hosts.remove(&viewport_id);
            self.pending_focus_viewports.retain(|id| *id != viewport_id);
        }
    }

    fn window_owns_auxiliary_viewport(
        window: &Arc<RwLock<Window>>,
        viewport_id: ViewportId,
    ) -> bool {
        window
            .read()
            .ok()
            .map(|guard| {
                guard
                    .collect_open_auxiliary_window_entries()
                    .iter()
                    .any(|(aux_viewport_id, _, _)| *aux_viewport_id == viewport_id)
            })
            .unwrap_or(false)
    }

    fn request_auxiliary_window_focus(&mut self, viewport_id: ViewportId) -> Option<ViewportId> {
        let visible_owner = self
            .windows
            .iter()
            .find(|(_, window)| Self::window_owns_auxiliary_viewport(window, viewport_id))
            .map(|(owner_viewport_id, window)| (*owner_viewport_id, window.clone()));
        if let Some((owner_viewport_id, window)) = visible_owner {
            if let Ok(mut guard) = window.write() {
                guard.request_focus_auxiliary_window(viewport_id);
            }
            return Some(owner_viewport_id);
        }

        let detached_owner = self
            .detached_auxiliary_window_hosts
            .values()
            .find(|window| Self::window_owns_auxiliary_viewport(window, viewport_id))
            .cloned();
        if let Some(window) = detached_owner
            && let Ok(mut guard) = window.write()
        {
            guard.request_focus_auxiliary_window(viewport_id);
        }
        None
    }

    pub(super) fn focus_window_viewport(&mut self, ctx: &egui::Context, viewport_id: ViewportId) {
        let auxiliary_owner_viewport = self.request_auxiliary_window_focus(viewport_id);
        let embedded_auxiliary_focus = ctx.embed_viewports() && auxiliary_owner_viewport.is_some();
        if let Some(owner_viewport_id) = auxiliary_owner_viewport {
            if ctx.embed_viewports() {
                if let Some(layer_id) =
                    self.embedded_window_layer_id_for_viewport(owner_viewport_id)
                {
                    ctx.move_to_top(layer_id);
                }
                ctx.send_viewport_cmd_to(ViewportId::ROOT, egui::ViewportCommand::Visible(true));
                ctx.send_viewport_cmd_to(ViewportId::ROOT, egui::ViewportCommand::Focus);
            } else {
                self.queue_focus_viewport(owner_viewport_id);
                ctx.send_viewport_cmd_to(owner_viewport_id, egui::ViewportCommand::Visible(true));
                ctx.send_viewport_cmd_to(owner_viewport_id, egui::ViewportCommand::Focus);
            }
        }

        if viewport_id == Self::configuration_viewport_id() {
            self.show_configuration_dialog = true;
        } else if viewport_id == Self::help_viewport_id() {
            self.show_help_dialog = true;
        } else if viewport_id == Self::command_palette_viewport_id() {
            self.show_command_palette_dialog = true;
            self.command_palette_focus_query = true;
        } else if viewport_id == Self::external_services_viewport_id() {
            self.external_services_ui.show_panel = true;
        } else if viewport_id == Self::history_viewport_id() {
            self.history_ui.show_panel = true;
        } else if viewport_id == Self::prepare_genome_viewport_id() {
            self.show_reference_genome_prepare_dialog = true;
        } else if viewport_id == Self::retrieve_genome_viewport_id() {
            self.show_reference_genome_retrieve_dialog = true;
        } else if viewport_id == Self::blast_genome_viewport_id() {
            self.show_reference_genome_blast_dialog = true;
        } else if viewport_id == Self::bed_track_viewport_id() {
            self.show_genome_bed_track_dialog = true;
        } else if viewport_id == Self::gibson_viewport_id() {
            self.show_gibson_dialog = true;
        } else if viewport_id == Self::pcr_design_viewport_id() {
            self.show_pcr_design_dialog = true;
        } else if viewport_id == Self::sequencing_confirmation_viewport_id() {
            self.show_sequencing_confirmation_dialog = true;
        } else if viewport_id == Self::planning_viewport_id() {
            self.show_planning_dialog = true;
        } else if viewport_id == Self::routine_assistant_viewport_id() {
            self.show_routine_assistant_dialog = true;
        } else if viewport_id == Self::agent_assistant_viewport_id() {
            self.show_agent_assistant_dialog = true;
        } else if viewport_id == Self::uniprot_viewport_id() {
            self.show_uniprot_dialog = true;
        } else if viewport_id == Self::new_sequence_viewport_id() {
            self.show_new_sequence_dialog = true;
        } else if viewport_id == Self::genbank_viewport_id() {
            self.show_genbank_dialog = true;
        }

        let focused_embedded_host = if ctx.embed_viewports()
            && let Some(layer_id) = self.embedded_window_layer_id_for_viewport(viewport_id)
        {
            ctx.move_to_top(layer_id);
            ctx.send_viewport_cmd_to(ViewportId::ROOT, egui::ViewportCommand::Visible(true));
            ctx.send_viewport_cmd_to(ViewportId::ROOT, egui::ViewportCommand::Focus);
            true
        } else {
            false
        };
        if !embedded_auxiliary_focus && !focused_embedded_host {
            ctx.send_viewport_cmd_to(viewport_id, egui::ViewportCommand::Visible(true));
            ctx.send_viewport_cmd_to(viewport_id, egui::ViewportCommand::Focus);
        }
        ctx.request_repaint();
        self.set_active_window_viewport(viewport_id);
    }

    pub(super) fn deferred_window_position(index: usize) -> Pos2 {
        Pos2 {
            x: index as f32 * 200.0,
            y: index as f32 * 200.0,
        }
    }

    pub(super) fn register_window(&mut self, mut window: Window) -> ViewportId {
        let id = format!("Viewport {}", self.viewport_id_counter);
        let id = ViewportId::from_hash_of(id);
        let position = Self::deferred_window_position(self.viewport_id_counter);
        self.viewport_id_counter += 1;
        self.mark_viewport_open_requested(id);
        self.pending_window_initial_positions.insert(id, position);
        window.set_window_scope_id(format!("{id:?}"));
        self.windows.insert(id, Arc::new(RwLock::new(window)));
        self.queue_focus_viewport(id);
        id
    }
}
