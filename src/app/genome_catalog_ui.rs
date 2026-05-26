//! Reference-genome catalog, cache, and track dialog helpers.
//!
//! This module is a move-only extraction from `app.rs`: it keeps the
//! reference-genome catalog dialogs, helper panels, cache cleanup UI, and
//! genome-track rendering helpers close to `GENtleApp` while reducing the
//! top-level app monolith.

use super::*;

impl GENtleApp {
    pub(super) fn render_host_profile_record_panel(
        ui: &mut Ui,
        profile: &HostProfileRecord,
        heading: &str,
    ) {
        ui.separator();
        ui.heading(heading);
        ui.monospace(format!("profile id: {}", profile.profile_id));
        ui.small(format!("species: {}", profile.species));
        ui.small(format!("strain: {}", profile.strain));
        if !profile.aliases.is_empty() {
            ui.small(format!("aliases: {}", profile.aliases.join(", ")));
        }
        if !profile.genotype_tags.is_empty() {
            ui.small(format!(
                "genotype tags: {}",
                profile.genotype_tags.join(", ")
            ));
        }
        if !profile.phenotype_tags.is_empty() {
            ui.small(format!(
                "phenotype tags: {}",
                profile.phenotype_tags.join(", ")
            ));
        }
        if !profile.notes.is_empty() {
            ui.collapsing("Notes", |ui| {
                for note in &profile.notes {
                    ui.small(note);
                }
            });
        }
        if !profile.source_notes.is_empty() {
            ui.collapsing("Source Notes", |ui| {
                for note in &profile.source_notes {
                    ui.small(note);
                }
            });
        }
    }

    fn format_helper_vector_card_component_summary(
        component: &crate::genomes::HelperConstructComponent,
    ) -> String {
        let title = component
            .label
            .as_deref()
            .filter(|value| !value.trim().is_empty())
            .unwrap_or(component.id.as_str());
        let mut parts = vec![format!("{title} [{}]", component.kind)];
        if title != component.id {
            parts.push(format!("id={}", component.id));
        }
        if !component.tags.is_empty() {
            parts.push(format!("tags={}", component.tags.join(",")));
        }
        if !component.attributes.is_empty() {
            let attrs = component
                .attributes
                .iter()
                .map(|(key, value)| format!("{key}={value}"))
                .collect::<Vec<_>>()
                .join(", ");
            parts.push(format!("attrs={attrs}"));
        }
        if let Some(description) = component
            .description
            .as_deref()
            .map(str::trim)
            .filter(|value| !value.is_empty())
        {
            parts.push(description.to_string());
        }
        parts.join(" | ")
    }

    pub(super) fn render_helper_vector_card_panel(ui: &mut Ui, card: &HelperVectorCard) {
        ui.separator();
        ui.heading("Vector Card");
        ui.monospace(format!("helper id: {}", card.helper_id));
        if let Some(helper_kind) = card
            .helper_kind
            .as_deref()
            .map(str::trim)
            .filter(|value| !value.is_empty())
        {
            ui.small(format!("helper kind: {helper_kind}"));
        }
        if let Some(host_system) = card
            .host_system
            .as_deref()
            .map(str::trim)
            .filter(|value| !value.is_empty())
        {
            ui.small(format!("host system: {host_system}"));
        }
        ui.small(format!(
            "metadata-only candidate: {}",
            card.metadata_only_candidate
        ));
        if let Some(usable) = card.usable_as_empty_backbone {
            ui.small(format!("usable as empty backbone: {usable}"));
        }
        if let Some(sequence_availability) = card
            .sequence_availability
            .as_deref()
            .map(str::trim)
            .filter(|value| !value.is_empty())
        {
            ui.small(format!("sequence availability: {sequence_availability}"));
        }
        if let Some(redistribution_status) = card
            .redistribution_status
            .as_deref()
            .map(str::trim)
            .filter(|value| !value.is_empty())
        {
            ui.small(format!("redistribution status: {redistribution_status}"));
        }
        if let Some(note) = card
            .biological_safety_note
            .as_deref()
            .map(str::trim)
            .filter(|value| !value.is_empty())
        {
            ui.small(format!("biological-safety note: {note}"));
        }
        if !card.affordances.is_empty() {
            ui.small(format!("affordances: {}", card.affordances.join(", ")));
        }
        if !card.constraints.is_empty() {
            ui.small(format!("constraints: {}", card.constraints.join(", ")));
        }
        if !card.components.is_empty() {
            ui.collapsing("Vector Card Components", |ui| {
                for component in &card.components {
                    ui.monospace(Self::format_helper_vector_card_component_summary(component));
                }
            });
        }
        if !card.relationships.is_empty() {
            ui.collapsing("Vector Card Relationships", |ui| {
                for relationship in &card.relationships {
                    ui.monospace(Self::format_helper_relationship_summary(relationship));
                }
            });
        }
    }

    pub(super) fn render_helper_catalog_entry_panel(
        ui: &mut Ui,
        entry: &GenomeCatalogListEntry,
        heading: &str,
    ) {
        ui.separator();
        ui.heading(heading);
        if let Some(summary) = entry
            .summary
            .as_deref()
            .map(str::trim)
            .filter(|value| !value.is_empty())
        {
            ui.small(summary);
        }
        if let Some(description) = entry
            .description
            .as_deref()
            .map(str::trim)
            .filter(|value| !value.is_empty())
        {
            let summary_trimmed = entry.summary.as_deref().map(str::trim).unwrap_or_default();
            if description != summary_trimmed {
                ui.small(description);
            }
        }
        if let Some(procurement) = entry.procurement.as_ref() {
            if let Some(vendor) = procurement
                .vendor_name
                .as_deref()
                .map(str::trim)
                .filter(|value| !value.is_empty())
            {
                ui.small(format!("vendor: {vendor}"));
            }
            if let Some(catalog_number) = procurement
                .catalog_number
                .as_deref()
                .map(str::trim)
                .filter(|value| !value.is_empty())
            {
                ui.small(format!("catalog number: {catalog_number}"));
            }
            if let Some(order_url) = procurement
                .order_url
                .as_deref()
                .map(str::trim)
                .filter(|value| !value.is_empty())
            {
                ui.small(format!("order URL: {order_url}"));
            }
            if let Some(reference_url) = procurement
                .reference_url
                .as_deref()
                .map(str::trim)
                .filter(|value| !value.is_empty())
            {
                ui.small(format!("reference URL: {reference_url}"));
            }
            if let Some(notes) = procurement
                .notes
                .as_deref()
                .map(str::trim)
                .filter(|value| !value.is_empty())
            {
                ui.small(format!("procurement notes: {notes}"));
            }
        }
        match entry.interpretation.as_ref() {
            Some(interpretation) => {
                for line in
                    Self::format_helper_construct_interpretation_detail_lines(interpretation)
                {
                    ui.small(line);
                }
                if !interpretation.components.is_empty() {
                    ui.collapsing("Components", |ui| {
                        for component in &interpretation.components {
                            ui.monospace(Self::format_helper_component_summary(component));
                        }
                    });
                }
                if !interpretation.relationships.is_empty() {
                    ui.collapsing("Relationships", |ui| {
                        for relationship in &interpretation.relationships {
                            ui.monospace(Self::format_helper_relationship_summary(relationship));
                        }
                    });
                }
            }
            None => {
                ui.small("No structured interpretation is available for this helper yet.");
            }
        }
    }

    pub(super) fn format_genome_source_plan_summary(plan: &GenomeSourcePlan) -> String {
        let mut parts = vec![format!(
            "sources: sequence={} | annotation={}",
            plan.sequence_source_type, plan.annotation_source_type
        )];
        if let Some(length_bp) = plan.nucleotide_length_bp {
            parts.push(format!("length={length_bp} bp"));
        }
        if let Some(mass_da) = plan.molecular_mass_da {
            let source = plan
                .molecular_mass_source
                .as_deref()
                .unwrap_or("unknown_source");
            parts.push(format!("mass={mass_da:.3e} Da ({source})"));
        }
        parts.join(" | ")
    }

    pub(super) fn choose_genome_from_catalog(
        ui: &mut Ui,
        genome_id: &mut String,
        names: &[String],
    ) -> bool {
        let mut changed = false;
        let selected_text = if names.iter().any(|name| name == genome_id) {
            genome_id.clone()
        } else if names.is_empty() {
            "(no genomes available)".to_string()
        } else {
            "(choose genome)".to_string()
        };
        egui::ComboBox::from_label("genome")
            .selected_text(selected_text)
            .show_ui(ui, |ui| {
                for name in names {
                    if ui.selectable_value(genome_id, name.clone(), name).changed() {
                        changed = true;
                    }
                }
            });
        changed
    }

    pub(super) fn collect_prepared_genome_inspections(
        &self,
    ) -> Result<(Vec<PreparedGenomeInspection>, Vec<String>), String> {
        let catalog_path = self.genome_catalog_path_resolved();
        let catalog = GenomeCatalog::from_json_file(&catalog_path)?;
        let cache_dir = self.genome_cache_dir_opt();
        let mut inspections: Vec<PreparedGenomeInspection> = vec![];
        let mut errors: Vec<String> = vec![];
        for genome_id in catalog.list_genomes() {
            match catalog.inspect_prepared_genome(&genome_id, cache_dir.as_deref()) {
                Ok(Some(inspection)) => inspections.push(inspection),
                Ok(None) => {}
                Err(e) => errors.push(format!("{genome_id}: {e}")),
            }
        }
        inspections.sort_by(|a, b| a.genome_id.cmp(&b.genome_id));
        Ok((inspections, errors))
    }

    pub(super) fn prepared_genome_chromosome_records(
        &self,
        genome_id: &str,
    ) -> Result<Vec<GenomeChromosomeRecord>, String> {
        let catalog_path = self.genome_catalog_path_resolved();
        let catalog = GenomeCatalog::from_json_file(&catalog_path)?;
        let cache_dir = self.genome_cache_dir_opt();
        catalog.list_chromosome_lengths(genome_id, cache_dir.as_deref())
    }

    pub(super) fn render_chromosome_length_lines(
        ui: &mut Ui,
        chromosomes: &[GenomeChromosomeRecord],
    ) {
        if chromosomes.is_empty() {
            ui.small("No chromosomes/contigs found in FASTA index.");
            return;
        }
        let longest_bp = chromosomes
            .iter()
            .map(|record| record.length_bp)
            .max()
            .unwrap_or(1)
            .max(1);
        let total_bp: u128 = chromosomes.iter().fold(0u128, |acc, record| {
            acc.saturating_add(record.length_bp as u128)
        });
        ui.small(format!(
            "contigs: {} | longest: {} ({} bp) | total span: {} bp",
            chromosomes.len(),
            chromosomes[0].chromosome,
            chromosomes[0].length_bp,
            total_bp
        ));
        egui::ScrollArea::vertical()
            .max_height(280.0)
            .show(ui, |ui| {
                scroll_input_policy::apply_scrollarea_keyboard_navigation(
                    ui,
                    scroll_input_policy::DEFAULT_SCROLLAREA_KEYBOARD_STEP,
                );
                for record in chromosomes {
                    ui.horizontal(|ui| {
                        ui.add_sized(
                            [220.0, 0.0],
                            egui::Label::new(
                                egui::RichText::new(record.chromosome.clone()).monospace(),
                            ),
                        );
                        ui.add_sized(
                            [120.0, 0.0],
                            egui::Label::new(format!("{} bp", record.length_bp)),
                        );
                        let width = ui.available_width().max(80.0);
                        let (rect, response) =
                            ui.allocate_exact_size(egui::vec2(width, 14.0), egui::Sense::hover());
                        let ratio = (record.length_bp as f32 / longest_bp as f32).clamp(0.0, 1.0);
                        let line_width = (rect.width() * ratio.max(0.005)).max(1.0);
                        let y = rect.center().y;
                        ui.painter().line_segment(
                            [
                                Pos2::new(rect.left(), y),
                                Pos2::new(rect.left() + line_width, y),
                            ],
                            egui::Stroke::new(2.0, egui::Color32::from_rgb(80, 130, 175)),
                        );
                        response.on_hover_text(format!(
                            "{} bp ({:.2}% of longest)",
                            record.length_bp,
                            ratio * 100.0
                        ));
                    });
                }
            });
    }

    pub(super) fn format_bytes_compact(bytes: u64) -> String {
        const UNITS: [&str; 5] = ["B", "KB", "MB", "GB", "TB"];
        let mut value = bytes as f64;
        let mut unit = 0usize;
        while value >= 1024.0 && unit + 1 < UNITS.len() {
            value /= 1024.0;
            unit += 1;
        }
        if unit == 0 {
            format!("{bytes} {}", UNITS[unit])
        } else {
            format!("{value:.2} {}", UNITS[unit])
        }
    }

    pub(super) fn format_bp_compact(bp: u64) -> String {
        const UNITS: [&str; 4] = ["bp", "kbp", "Mbp", "Gbp"];
        let mut value = bp as f64;
        let mut unit = 0usize;
        while value >= 1000.0 && unit + 1 < UNITS.len() {
            value /= 1000.0;
            unit += 1;
        }
        if unit == 0 {
            format!("{bp} {}", UNITS[unit])
        } else {
            format!("{value:.2} {}", UNITS[unit])
        }
    }

    pub(super) fn format_short_sha1(value: &Option<String>) -> String {
        value
            .as_ref()
            .map(|v| v.trim())
            .filter(|v| !v.is_empty())
            .map(|v| {
                if v.len() > 12 {
                    v[..12].to_string()
                } else {
                    v.to_string()
                }
            })
            .unwrap_or_else(|| "-".to_string())
    }

    pub(super) fn format_prepared_cache_summary(inspection: &PreparedGenomeInspection) -> String {
        if inspection.cached_contig_count == 0 {
            return "-".to_string();
        }
        let longest = inspection.cached_longest_contig.as_deref().unwrap_or("?");
        let longest_bp = inspection.cached_longest_contig_bp.unwrap_or(0);
        format!(
            "{} ctgs | {} {}",
            inspection.cached_contig_count,
            longest,
            Self::format_bp_compact(longest_bp)
        )
    }

    pub(super) fn format_prepared_cache_hover_text(
        inspection: &PreparedGenomeInspection,
    ) -> String {
        if inspection.cached_contig_count == 0 {
            return "No cached contig summary available from the FASTA index.".to_string();
        }
        let longest = inspection.cached_longest_contig.as_deref().unwrap_or("?");
        let longest_bp = inspection.cached_longest_contig_bp.unwrap_or(0);
        let preview = if inspection.cached_contig_preview.is_empty() {
            "-".to_string()
        } else {
            inspection.cached_contig_preview.join(", ")
        };
        format!(
            "contigs: {}\ntotal span: {}\nlongest: {} ({})\npreview: {}",
            inspection.cached_contig_count,
            Self::format_bp_compact(inspection.cached_total_span_bp),
            longest,
            Self::format_bp_compact(longest_bp),
            preview
        )
    }

    pub(super) fn prepare_step_summary(
        &self,
    ) -> Option<(String, usize, usize, f32, Option<Duration>)> {
        let total_steps = self.genome_prepare_steps.len();
        if total_steps == 0 {
            return None;
        }
        let completed_steps = self
            .genome_prepare_steps
            .iter()
            .filter(|step| step.status == PrepareGenomeUiStepStatus::Completed)
            .count();
        let running_step = self
            .genome_prepare_steps
            .iter()
            .find(|step| step.status == PrepareGenomeUiStepStatus::Running);
        let current_label = running_step
            .map(|step| step.label.clone())
            .or_else(|| {
                self.genome_prepare_progress
                    .as_ref()
                    .and_then(|progress| progress.step_label.clone())
            })
            .or_else(|| {
                self.genome_prepare_steps
                    .iter()
                    .rev()
                    .find(|step| !matches!(step.status, PrepareGenomeUiStepStatus::Pending))
                    .map(|step| step.label.clone())
            })
            .unwrap_or_else(|| "Queued".to_string());
        let running_fraction = running_step
            .and_then(|step| step.progress_fraction)
            .unwrap_or(0.0);
        let overall = ((completed_steps as f32) + running_fraction) / total_steps as f32;
        let eta_remaining = running_step.and_then(|step| step.eta_remaining);
        Some((
            current_label,
            completed_steps,
            total_steps,
            overall.clamp(0.0, 1.0),
            eta_remaining,
        ))
    }

    pub(super) fn prepare_step_state_status_label(
        status: PrepareGenomeUiStepStatus,
    ) -> &'static str {
        match status {
            PrepareGenomeUiStepStatus::Pending => "pending",
            PrepareGenomeUiStepStatus::Running => "running",
            PrepareGenomeUiStepStatus::Completed => "completed",
            PrepareGenomeUiStepStatus::Failed => "failed",
            PrepareGenomeUiStepStatus::Cancelled => "cancelled",
        }
    }

    pub(super) fn prepare_step_state_icon(status: PrepareGenomeUiStepStatus) -> &'static str {
        match status {
            PrepareGenomeUiStepStatus::Pending => "○",
            PrepareGenomeUiStepStatus::Running => "◌",
            PrepareGenomeUiStepStatus::Completed => "✓",
            PrepareGenomeUiStepStatus::Failed => "!",
            PrepareGenomeUiStepStatus::Cancelled => "-",
        }
    }

    pub(super) fn prepare_step_state_color(status: PrepareGenomeUiStepStatus) -> egui::Color32 {
        match status {
            PrepareGenomeUiStepStatus::Pending => egui::Color32::from_gray(130),
            PrepareGenomeUiStepStatus::Running => egui::Color32::from_rgb(50, 110, 180),
            PrepareGenomeUiStepStatus::Completed => egui::Color32::from_rgb(50, 135, 70),
            PrepareGenomeUiStepStatus::Failed => egui::Color32::from_rgb(190, 70, 70),
            PrepareGenomeUiStepStatus::Cancelled => egui::Color32::from_rgb(170, 120, 50),
        }
    }

    pub(super) fn render_prepare_step_checklist(&self, ui: &mut Ui) {
        if self.genome_prepare_steps.is_empty() {
            return;
        }
        ui.separator();
        ui.strong("Planned steps");
        if let Some((current_step, completed_steps, total_steps, overall, eta_remaining)) =
            self.prepare_step_summary()
        {
            let running_step = self
                .genome_prepare_steps
                .iter()
                .find(|step| step.status == PrepareGenomeUiStepStatus::Running);
            if let Some(_step) = running_step {
                let mut summary = format!(
                    "Current step: {} ({}/{})",
                    current_step, completed_steps, total_steps
                );
                if let Some(eta) = eta_remaining {
                    summary.push_str(&format!(" • ETA {}", Self::format_duration_compact(eta)));
                }
                ui.small(summary);
            } else if completed_steps == total_steps {
                ui.colored_label(
                    egui::Color32::from_rgb(50, 135, 70),
                    "All planned steps completed.",
                );
            } else {
                ui.small(format!("Completed: {completed_steps}/{total_steps}"));
            }
            ui.add(
                egui::ProgressBar::new(overall)
                    .show_percentage()
                    .text(format!("{completed_steps}/{total_steps} steps")),
            );
        }
        for step in &self.genome_prepare_steps {
            ui.add_space(4.0);
            let color = Self::prepare_step_state_color(step.status);
            let indeterminate_active = step.status == PrepareGenomeUiStepStatus::Running
                && (!step.determinate_hint || step.progress_fraction.is_none());
            ui.horizontal(|ui| {
                if indeterminate_active {
                    ui.add(egui::Spinner::new());
                } else {
                    ui.colored_label(color, Self::prepare_step_state_icon(step.status));
                }
                ui.label(egui::RichText::new(&step.label).strong());
                ui.label(
                    egui::RichText::new(Self::prepare_step_state_status_label(step.status))
                        .small()
                        .color(color),
                );
            });
            let detail = if !step.detail.trim().is_empty() {
                step.detail.as_str()
            } else {
                step.operation_summary.as_str()
            };
            ui.horizontal(|ui| {
                ui.add_space(18.0);
                ui.label(
                    egui::RichText::new(detail)
                        .small()
                        .color(egui::Color32::from_gray(120)),
                );
            });
            if step.status == PrepareGenomeUiStepStatus::Running
                && let Some(bytes_total) = step.bytes_total.filter(|total| *total > 0)
            {
                let bytes_done = step.bytes_done.unwrap_or(0).min(bytes_total);
                ui.horizontal(|ui| {
                    ui.add_space(18.0);
                    ui.label(
                        egui::RichText::new(Self::format_prepare_byte_progress(
                            bytes_done,
                            bytes_total,
                            step.eta_remaining,
                        ))
                        .small()
                        .color(egui::Color32::from_gray(120)),
                    );
                });
            }
        }
    }

    pub(super) fn render_reference_genome_prepare_contents(&mut self, ui: &mut Ui) -> bool {
        self.refresh_genome_catalog_list();
        let scope = self.genome_dialog_scope;
        let mut close_requested = false;
        let close_hover = Self::specialist_window_close_hover_text(scope.prepare_title());
        if self.render_specialist_window_nav_with_close(ui, Some(("Close", close_hover.as_str()))) {
            close_requested = true;
        }
        let mut prepare_context_changed = false;
        ui.label(format!(
            "Download and index a {} once.",
            scope.description()
        ));
        ui.horizontal(|ui| {
            ui.label("catalog");
            let response = ui.text_edit_singleline(&mut self.genome_catalog_path);
            if response.changed() {
                prepare_context_changed = true;
            }
            if ui
                .button("Browse...")
                .on_hover_text("Browse filesystem and fill this path")
                .clicked()
                && let Some(path) = rfd::FileDialog::new()
                    .add_filter("JSON", &["json"])
                    .pick_file()
            {
                self.genome_catalog_path = path.display().to_string();
                prepare_context_changed = true;
            }
        });
        ui.horizontal(|ui| {
            ui.label("cache_dir");
            let response = ui.text_edit_singleline(&mut self.genome_cache_dir);
            if response.changed() {
                prepare_context_changed = true;
            }
            if ui
                .button("Browse...")
                .on_hover_text("Browse filesystem and fill this path")
                .clicked()
                && let Some(path) = rfd::FileDialog::new().pick_folder()
            {
                self.genome_cache_dir = path.display().to_string();
                prepare_context_changed = true;
            }
        });
        ui.horizontal(|ui| {
            ui.label("timeout_sec");
            ui.add(
                egui::TextEdit::singleline(&mut self.genome_prepare_timeout_secs)
                    .desired_width(90.0),
            );
            ui.small("optional; empty or 0 means no timebox");
        });
        if !self.genome_catalog_error.is_empty() {
            ui.colored_label(
                egui::Color32::from_rgb(190, 70, 70),
                format!("Catalog error: {}", self.genome_catalog_error),
            );
        }
        let all_genomes = self.genome_catalog_genomes.clone();
        let preparable_genomes = match self.unprepared_genomes_for_prepare_dialog() {
            Ok(names) => names,
            Err(e) => {
                ui.colored_label(
                    egui::Color32::from_rgb(190, 70, 70),
                    format!("Prepared-state check error: {e}"),
                );
                vec![]
            }
        };
        let preparable_set: HashSet<String> = preparable_genomes.iter().cloned().collect();
        let mut selection_changed = false;
        let selected_text = if self.genome_id.trim().is_empty() {
            "Select genome".to_string()
        } else {
            self.genome_id.clone()
        };
        egui::ComboBox::from_label("genome")
            .selected_text(selected_text)
            .show_ui(ui, |ui| {
                for genome_name in &all_genomes {
                    let is_preparable = preparable_set.contains(genome_name);
                    let item_label = if is_preparable {
                        genome_name.clone()
                    } else {
                        format!("{genome_name} (already prepared)")
                    };
                    if ui
                        .selectable_label(self.genome_id == *genome_name, item_label)
                        .clicked()
                    {
                        self.genome_id = genome_name.clone();
                        selection_changed = true;
                        ui.close();
                    }
                }
            });
        if selection_changed {
            self.invalidate_genome_genes();
            prepare_context_changed = true;
        }
        if prepare_context_changed && self.genome_prepare_task.is_none() {
            if !selection_changed {
                self.invalidate_genome_genes();
            }
            self.reset_prepare_dialog_preview_state();
        }
        match self.selected_genome_source_plan() {
            Ok(Some(plan)) => {
                ui.small(Self::format_genome_source_plan_summary(&plan));
            }
            Ok(None) => {}
            Err(e) => {
                ui.colored_label(
                    egui::Color32::from_rgb(190, 70, 70),
                    format!("Source-plan error: {e}"),
                );
            }
        }
        if matches!(scope, GenomeDialogScope::Helper) {
            match self.selected_helper_catalog_entry_for_scope(scope) {
                Ok(Some(entry)) => {
                    Self::render_helper_catalog_entry_panel(
                        ui,
                        &entry,
                        "Helper Construct Interpretation",
                    );
                }
                Ok(None) => {}
                Err(e) => {
                    ui.colored_label(
                        egui::Color32::from_rgb(190, 70, 70),
                        format!("Helper interpretation error: {e}"),
                    );
                }
            }
        }
        if preparable_genomes.is_empty() {
            ui.label("All genomes in this catalog are already prepared.");
        }
        let running = self.genome_prepare_task.is_some();
        let primary_action =
            Self::prepare_dialog_primary_action(&self.genome_id, &all_genomes, &preparable_set);
        if !running {
            self.ensure_prepare_step_preview_for_current_selection(primary_action);
        }
        if matches!(primary_action, PrepareGenomeDialogPrimaryAction::Reindex) {
            ui.label(
                "Selected genome is already prepared. Use the main action below to reindex it from cached local files.",
            );
            ui.small(
                "Reindex keeps the downloaded local files by default. Removing cached files and re-downloading now requires explicit confirmation.",
            );
        } else if !self.genome_id.trim().is_empty()
            && matches!(primary_action, PrepareGenomeDialogPrimaryAction::None)
        {
            ui.label("Selected genome is not present in this catalog.");
        }
        if let Some(task) = &self.genome_prepare_task {
            ui.horizontal(|ui| {
                ui.add(egui::Spinner::new());
                let mut status = format!(
                    "{} task running ({:.1}s)",
                    task.mode.progress_label(),
                    task.started.elapsed().as_secs_f32()
                );
                if let Some(timeout) = task.timeout_seconds {
                    status.push_str(&format!(", timeout={}s", timeout));
                }
                if task.cancel_requested.load(Ordering::Relaxed) {
                    status.push_str(", cancellation requested");
                }
                ui.label(status);
            });
        }
        ui.horizontal(|ui| {
            let (primary_label, primary_hover) = match primary_action {
                PrepareGenomeDialogPrimaryAction::Prepare => (
                    "Prepare Genome",
                    format!("Download and index the selected {}", scope.description()),
                ),
                PrepareGenomeDialogPrimaryAction::Reindex => (
                    "Reindex Selected...",
                    format!(
                        "Reuse cached local files and rebuild indexes for the selected {}. A confirmation dialog will be shown.",
                        scope.description()
                    ),
                ),
                PrepareGenomeDialogPrimaryAction::None => (
                    "Prepare Genome",
                    format!("Select a {} from this catalog first", scope.description()),
                ),
            };
            if ui
                .add_enabled(
                    !running
                        && !matches!(primary_action, PrepareGenomeDialogPrimaryAction::None),
                    egui::Button::new(primary_label),
                )
                .on_hover_text(primary_hover)
                .clicked()
            {
                match primary_action {
                    PrepareGenomeDialogPrimaryAction::Prepare => {
                        self.start_prepare_reference_genome();
                    }
                    PrepareGenomeDialogPrimaryAction::Reindex => {
                        self.queue_prepared_genome_reinstall(
                            self.genome_id.clone(),
                            scope,
                            PreparedGenomeReinstallDialogHost::PrepareDialog,
                        );
                    }
                    PrepareGenomeDialogPrimaryAction::None => {}
                }
            }
            if self.genome_prepare_task.is_some()
                && ui
                    .button("Cancel Prepare")
                    .on_hover_text("Request cancellation of the running prepare task.")
                    .clicked()
                {
                    self.request_prepare_task_cancel("prepare dialog");
                }
        });
        ui.horizontal(|ui| {
            let can_update_ensembl = !running && self.selected_genome_catalog_has_ensembl_templates();
            let hover = if can_update_ensembl {
                "Preview and apply pinned Ensembl URL/spec refreshes for catalog entries that carry Ensembl template metadata."
            } else {
                "Current catalog has no Ensembl template metadata entries to update."
            };
            if ui
                .add_enabled(can_update_ensembl, egui::Button::new("Update Ensembl Specs..."))
                .on_hover_text(hover)
                .clicked()
            {
                self.queue_ensembl_catalog_update_preview();
            }
            if ui
                .add_enabled(
                    !running,
                    egui::Button::new("Browse Ensembl Candidates..."),
                )
                .on_hover_text(
                    "Inspect species currently exposed by Ensembl/Ensembl Metazoa where both FASTA and GTF listings are present. This does not edit the catalog.",
                )
                .clicked()
            {
                self.queue_ensembl_installable_genome_discovery(Some("all"), None);
            }
        });
        if !self.genome_prepare_steps.is_empty() {
            self.render_prepare_step_checklist(ui);
        } else if let Some(progress) = &self.genome_prepare_progress {
            let fraction = Self::prepare_progress_fraction(progress).unwrap_or(0.0);
            let mut bar = egui::ProgressBar::new(fraction)
                .text(format!("{}: {}", progress.phase, progress.item));
            if progress.percent.is_some() || progress.bytes_total.is_some() {
                bar = bar.show_percentage();
            }
            ui.add(bar);
            let bytes_total = progress
                .bytes_total
                .map(|b| b.to_string())
                .unwrap_or_else(|| "?".to_string());
            if progress.bytes_done > 0 || progress.bytes_total.is_some() {
                ui.small(format!("bytes: {} / {}", progress.bytes_done, bytes_total));
            }
        }
        if self.genome_prepare_task.is_none() && self.genome_prepare_failure_recovery.is_some() {
            ui.separator();
            ui.colored_label(
                egui::Color32::from_rgb(170, 120, 50),
                "Cached prepared files are inconsistent; reinstall from sources is recommended.",
            );
            if ui
                .button("Reinstall From Sources...")
                .on_hover_text(
                    "Discard the stale prepared install for this genome and rebuild it from the configured sources.",
                )
                .clicked()
            {
                self.queue_prepare_failure_reinstall(
                    PreparedGenomeReinstallDialogHost::PrepareDialog,
                );
            }
        }
        if !self.genome_prepare_status.is_empty() {
            ui.separator();
            ui.monospace(&self.genome_prepare_status);
        }
        close_requested
    }

    pub(super) fn render_reference_genome_prepare_scroll_area(
        &mut self,
        ui: &mut Ui,
        id_salt: &'static str,
    ) -> egui::containers::scroll_area::ScrollAreaOutput<bool> {
        egui::ScrollArea::vertical()
            .id_salt(id_salt)
            .auto_shrink([false, false])
            .show(ui, |ui| self.render_reference_genome_prepare_contents(ui))
    }

    pub(super) fn render_reference_genome_prepare_dialog(&mut self, ctx: &egui::Context) {
        if !self.show_reference_genome_prepare_dialog {
            return;
        }
        let title = self.genome_dialog_scope.prepare_title();
        let viewport_id = Self::prepare_genome_viewport_id();
        let spec = self.hosted_window_spec_for_viewport(
            title,
            egui::Id::new(("hosted_prepare_genome_window", viewport_id)),
            viewport_id,
            Vec2::new(760.0, 560.0),
            Vec2::new(520.0, 360.0),
        );

        if ctx.embed_viewports() {
            let mut open = self.show_reference_genome_prepare_dialog;
            let mut close_requested = false;
            crate::egui_compat::show_hosted_window(ctx, &spec, &mut open, |ui| {
                close_requested = self
                    .render_reference_genome_prepare_scroll_area(
                        ui,
                        "prepare_genome_embedded_scroll",
                    )
                    .inner;
            });
            self.clear_viewport_foreground_request_after_render(viewport_id);
            if close_requested {
                open = false;
                self.dismiss_pending_prepared_genome_reinstall_for_host(
                    PreparedGenomeReinstallDialogHost::PrepareDialog,
                );
            }
            self.render_prepared_genome_reinstall_confirm_dialog(
                ctx,
                PreparedGenomeReinstallDialogHost::PrepareDialog,
            );
            self.show_reference_genome_prepare_dialog = open;
            if !self.show_reference_genome_prepare_dialog {
                self.dismiss_pending_prepared_genome_reinstall_for_host(
                    PreparedGenomeReinstallDialogHost::PrepareDialog,
                );
                self.clear_prepare_dialog_ephemeral_state();
            }
            return;
        }

        let builder = crate::egui_compat::viewport_builder_for_hosted_window(&spec);
        ctx.show_viewport_immediate(viewport_id, builder, |ctx, class| {
            self.note_viewport_focus_if_active(ctx, viewport_id);
            if class == egui::ViewportClass::EmbeddedWindow {
                let mut open = self.show_reference_genome_prepare_dialog;
                let mut close_requested = false;
                crate::egui_compat::show_hosted_window(ctx, &spec, &mut open, |ui| {
                    close_requested = self
                        .render_reference_genome_prepare_scroll_area(
                            ui,
                            "prepare_genome_embedded_scroll",
                        )
                        .inner;
                });
                if close_requested {
                    open = false;
                    self.dismiss_pending_prepared_genome_reinstall_for_host(
                        PreparedGenomeReinstallDialogHost::PrepareDialog,
                    );
                }
                self.render_prepared_genome_reinstall_confirm_dialog(
                    ctx,
                    PreparedGenomeReinstallDialogHost::PrepareDialog,
                );
                self.show_reference_genome_prepare_dialog = open;
                if !self.show_reference_genome_prepare_dialog {
                    self.dismiss_pending_prepared_genome_reinstall_for_host(
                        PreparedGenomeReinstallDialogHost::PrepareDialog,
                    );
                    self.clear_prepare_dialog_ephemeral_state();
                }
                return;
            }

            let mut close_requested = false;
            crate::egui_compat::show_central_panel(ctx, egui::CentralPanel::default(), |ui| {
                close_requested = self
                    .render_reference_genome_prepare_scroll_area(
                        ui,
                        "prepare_genome_viewport_scroll",
                    )
                    .inner;
            });
            self.render_prepared_genome_reinstall_confirm_dialog(
                ctx,
                PreparedGenomeReinstallDialogHost::PrepareDialog,
            );

            if close_requested || Self::viewport_close_requested_or_shortcut(ctx) {
                self.dismiss_pending_prepared_genome_reinstall_for_host(
                    PreparedGenomeReinstallDialogHost::PrepareDialog,
                );
                self.show_reference_genome_prepare_dialog = false;
                self.clear_prepare_dialog_ephemeral_state();
            }
        });
    }

    pub(super) fn render_prepared_genome_reinstall_confirm_dialog(
        &mut self,
        ctx: &egui::Context,
        dialog_host: PreparedGenomeReinstallDialogHost,
    ) {
        if !self.pending_prepared_genome_reinstall_targets_host(dialog_host) {
            return;
        }
        let Some(request) = self.pending_prepared_genome_reinstall.clone() else {
            return;
        };
        let resolved_catalog = if request.catalog_path.trim().is_empty() {
            request.scope.default_catalog_path().to_string()
        } else {
            request.catalog_path.trim().to_string()
        };
        let resolved_cache = if request.cache_dir.trim().is_empty() {
            request.scope.default_cache_dir().to_string()
        } else {
            request.cache_dir.trim().to_string()
        };
        let mut open = true;
        let spec = crate::egui_compat::ModalWindowSpec::new(
            "Reindex Prepared Genome",
            egui::Id::new((
                "reindex_prepared_genome_modal",
                matches!(request.scope, GenomeDialogScope::Helper),
                request.genome_id.clone(),
            )),
        );
        crate::egui_compat::show_modal_window(ctx, &spec, &mut open, |ui| {
            ui.label(format!(
                "Reindex prepared {} '{}'?",
                request.scope.description(),
                request.genome_id
            ));
            ui.small(
                    "Reindex keeps the cached local sequence/annotation files and rebuilds local indexes in the selected cache.",
                );
            ui.small(
                    "Removing cached files and re-downloading from sources is available below, but only as an explicit action.",
                );
            ui.separator();
            ui.label(format!("catalog: {resolved_catalog}"));
            ui.label(format!("cache_dir: {resolved_cache}"));
            if self.genome_prepare_task.is_some() {
                ui.separator();
                ui.colored_label(
                        egui::Color32::from_rgb(190, 70, 70),
                        "Another genome prepare task is already running. Finish or cancel it before starting a reindex or refresh.",
                    );
            }
            ui.horizontal(|ui| {
                    if ui
                        .add_enabled(
                            self.genome_prepare_task.is_none(),
                            egui::Button::new("Reindex Using Cached Files"),
                        )
                        .on_hover_text(
                            "Keep the cached local files and rebuild indexes using the shown catalog/cache settings. This may take some time.",
                        )
                        .clicked()
                    {
                        self.confirm_prepared_genome_reinstall(
                            GenomePrepareLaunchMode::ReindexCachedFiles,
                        );
                    }
                    if ui
                        .add_enabled(
                            self.genome_prepare_task.is_none(),
                            egui::Button::new("Remove Cached Files + Re-download"),
                        )
                        .on_hover_text(
                            "Explicitly delete cached local files, download sources again, and rebuild all indexes.",
                        )
                        .clicked()
                    {
                        self.confirm_prepared_genome_reinstall(
                            GenomePrepareLaunchMode::RefreshFromSources,
                        );
                    }
                    if ui
                        .button("Cancel")
                        .on_hover_text("Keep the current installation unchanged")
                        .clicked()
                    {
                        self.pending_prepared_genome_reinstall = None;
                    }
                });
        });
        if !open {
            self.pending_prepared_genome_reinstall = None;
        }
    }

    pub(super) fn render_pending_ensembl_catalog_update_dialog(&mut self, ctx: &egui::Context) {
        let Some(mut dialog) = self.pending_ensembl_catalog_update.clone() else {
            return;
        };
        let mut open = true;
        let spec = crate::egui_compat::ModalWindowSpec::new(
            "Update Ensembl Specs",
            egui::Id::new((
                "update_ensembl_specs_modal",
                matches!(dialog.scope, GenomeDialogScope::Helper),
                dialog.catalog_path.clone(),
            )),
        )
        .default_size(Vec2::new(760.0, 520.0))
        .min_size(Vec2::new(560.0, 320.0))
        .resizable(true);
        crate::egui_compat::show_modal_window(ctx, &spec, &mut open, |ui| {
            ui.label(format!(
                "Refresh pinned Ensembl catalog specs for the {} catalog?",
                dialog.scope.description()
            ));
            ui.small(
                    "Older pinned releases stay in the catalog. Newer release rows are added or refreshed without downloading/preparing genomes.",
                );
            ui.separator();
            ui.label(format!("catalog: {}", dialog.catalog_path));
            if dialog.preview.writable_catalog {
                ui.small("This catalog is writable and can be updated in place.");
            } else {
                ui.colored_label(
                        egui::Color32::from_rgb(170, 120, 50),
                        "This catalog is not writable in place. Choose a save path for an updated copy.",
                    );
            }
            ui.horizontal(|ui| {
                ui.label("output");
                ui.text_edit_singleline(&mut dialog.output_catalog_path);
                if ui
                    .button("Browse...")
                    .on_hover_text("Pick a destination for an updated catalog copy")
                    .clicked()
                    && let Some(path) = rfd::FileDialog::new()
                        .add_filter("JSON", &["json"])
                        .save_file()
                {
                    dialog.output_catalog_path = path.display().to_string();
                }
            });
            if !dialog.preview.collection_latest_releases.is_empty() {
                let mut collections: Vec<_> =
                    dialog.preview.collection_latest_releases.iter().collect();
                collections.sort_by(|left, right| left.0.cmp(right.0));
                for (collection, release) in collections {
                    ui.small(format!("latest {collection} release seen: {release}"));
                }
            }
            ui.separator();
            ui.label(format!(
                "Updates: {} | unchanged: {} | skipped: {}",
                dialog.preview.updates.len(),
                dialog.preview.unchanged_entries.len(),
                dialog.preview.skipped_entries.len()
            ));
            egui::ScrollArea::vertical()
                .id_salt("ensembl_catalog_update_updates_scroll")
                .max_height(220.0)
                .show(ui, |ui| {
                    scroll_input_policy::apply_scrollarea_keyboard_navigation(
                        ui,
                        scroll_input_policy::DEFAULT_SCROLLAREA_KEYBOARD_STEP,
                    );
                    for update in &dialog.preview.updates {
                        ui.group(|ui| {
                            ui.label(format!(
                                "{} -> {} ({}, {} -> {})",
                                update.source_genome_id,
                                update.target_genome_id,
                                update.action,
                                update.old_release,
                                update.new_release
                            ));
                            ui.small(format!("sequence: {}", update.new_sequence_remote));
                            ui.small(format!("annotation: {}", update.new_annotations_remote));
                        });
                    }
                });
            if !dialog.preview.skipped_entries.is_empty() {
                ui.separator();
                ui.label("Skipped:");
                for row in &dialog.preview.skipped_entries {
                    ui.small(row);
                }
            }
            if !dialog.preview.warnings.is_empty() {
                ui.separator();
                ui.label("Warnings:");
                for warning in &dialog.preview.warnings {
                    ui.small(warning);
                }
            }
            ui.separator();
            let needs_output_copy =
                !dialog.preview.writable_catalog && dialog.output_catalog_path.trim().is_empty();
            ui.horizontal(|ui| {
                if ui
                    .add_enabled(
                        !needs_output_copy,
                        egui::Button::new(if dialog.preview.writable_catalog {
                            "Update Catalog"
                        } else {
                            "Save Updated Copy"
                        }),
                    )
                    .on_hover_text(
                        "Write the refreshed Ensembl entries to the selected catalog path.",
                    )
                    .clicked()
                {
                    self.pending_ensembl_catalog_update = Some(dialog.clone());
                    self.apply_pending_ensembl_catalog_update();
                }
                if ui
                    .button("Cancel")
                    .on_hover_text("Close this preview without changing the catalog")
                    .clicked()
                {
                    self.pending_ensembl_catalog_update = None;
                }
            });
            if self.pending_ensembl_catalog_update.is_some() {
                self.pending_ensembl_catalog_update = Some(dialog);
            }
        });
        if !open {
            self.pending_ensembl_catalog_update = None;
        }
    }

    pub(super) fn render_pending_ensembl_installable_genomes_dialog(
        &mut self,
        ctx: &egui::Context,
    ) {
        let Some(mut dialog) = self.pending_ensembl_installable_genomes.clone() else {
            return;
        };
        let mut reload_requested = false;
        let mut close_requested = false;
        let mut open = true;
        let spec = crate::egui_compat::ModalWindowSpec::new(
            "Browse Ensembl Candidates",
            egui::Id::new((
                "browse_ensembl_candidates_modal",
                matches!(dialog.scope, GenomeDialogScope::Helper),
            )),
        )
        .default_size(Vec2::new(820.0, 520.0))
        .min_size(Vec2::new(560.0, 320.0))
        .resizable(true);
        crate::egui_compat::show_modal_window(ctx, &spec, &mut open, |ui| {
            ui.label(format!(
                "Inspect Ensembl-discovered candidates for the {} workflow.",
                dialog.scope.description()
            ));
            ui.small(
                    "This is a read-only discovery view. Candidates appear when both current FASTA and GTF species-directory listings are present; no catalog rows are created here.",
                );
            ui.small(format!(
                "Availability basis: {}",
                dialog.report.availability_basis
            ));
            ui.separator();
            ui.horizontal(|ui| {
                ui.label("collection");
                egui::ComboBox::from_id_salt("ensembl_installable_collection_filter")
                    .selected_text(dialog.collection_filter.as_str())
                    .show_ui(ui, |ui| {
                        ui.selectable_value(
                            &mut dialog.collection_filter,
                            "all".to_string(),
                            "all",
                        );
                        ui.selectable_value(
                            &mut dialog.collection_filter,
                            "vertebrates".to_string(),
                            "vertebrates",
                        );
                        ui.selectable_value(
                            &mut dialog.collection_filter,
                            "metazoa".to_string(),
                            "metazoa",
                        );
                    });
                ui.label("filter");
                ui.add(egui::TextEdit::singleline(&mut dialog.filter).desired_width(f32::INFINITY));
                if ui
                    .button("Clear")
                    .on_hover_text("Clear the Ensembl candidate text filter")
                    .clicked()
                {
                    dialog.filter.clear();
                }
            });
            ui.horizontal(|ui| {
                if ui
                    .button("Reload")
                    .on_hover_text(
                        "Refetch current Ensembl listings with the selected collection/filter",
                    )
                    .clicked()
                {
                    reload_requested = true;
                }
                if ui
                    .button("Close")
                    .on_hover_text("Close the Ensembl candidate browser")
                    .clicked()
                {
                    close_requested = true;
                }
            });
            if !dialog.report.collection_latest_releases.is_empty() {
                let mut collections: Vec<_> =
                    dialog.report.collection_latest_releases.iter().collect();
                collections.sort_by(|left, right| left.0.cmp(right.0));
                for (collection, release) in collections {
                    ui.small(format!("latest {collection} release seen: {release}"));
                }
            }
            ui.separator();
            ui.label(format!(
                "Candidates: {} (loaded for collection='{}')",
                dialog.report.candidates.len(),
                dialog.report.collection_filter
            ));
            egui::ScrollArea::vertical()
                    .id_salt("ensembl_candidates_scroll")
                    .max_height(320.0)
                    .show(ui, |ui| {
                        scroll_input_policy::apply_scrollarea_keyboard_navigation(
                            ui,
                            scroll_input_policy::DEFAULT_SCROLLAREA_KEYBOARD_STEP,
                        );
                        if dialog.report.candidates.is_empty() {
                            ui.small("No Ensembl candidates matched the current query.");
                        } else {
                            for candidate in &dialog.report.candidates {
                                ui.group(|ui| {
                                    ui.horizontal(|ui| {
                                        ui.label(format!(
                                            "{} ({})",
                                            candidate.display_name, candidate.collection
                                        ));
                                        if ui
                                            .add_enabled(
                                                self.genome_prepare_task.is_none(),
                                                egui::Button::new("Quick Install..."),
                                            )
                                            .on_hover_text(
                                                "Resolve concrete Ensembl URLs for this species, write a local catalog entry, and then start the normal prepare workflow.",
                                            )
                                            .clicked()
                                        {
                                            self.pending_ensembl_installable_genomes =
                                                Some(dialog.clone());
                                            self.queue_ensembl_quick_install_preview(
                                                dialog.scope,
                                                &candidate.collection,
                                                &candidate.species_dir,
                                            );
                                        }
                                    });
                                    ui.small(format!("species_dir: {}", candidate.species_dir));
                                    ui.small(format!("latest release seen: {}", candidate.latest_release));
                                    ui.small(format!(
                                        "FASTA listing: {}",
                                        candidate.current_fasta_listing_url
                                    ));
                                    ui.small(format!(
                                        "GTF listing: {}",
                                        candidate.current_gtf_listing_url
                                    ));
                                });
                            }
                        }
                    });
            if !dialog.report.warnings.is_empty() {
                ui.separator();
                ui.label("Warnings:");
                for warning in &dialog.report.warnings {
                    ui.small(warning);
                }
            }
        });
        if !open || close_requested {
            self.pending_ensembl_installable_genomes = None;
            return;
        }
        if reload_requested {
            let collection_filter = dialog.collection_filter.clone();
            let filter = dialog.filter.clone();
            self.queue_ensembl_installable_genome_discovery(
                Some(collection_filter.as_str()),
                Some(filter.as_str()),
            );
            return;
        }
        self.pending_ensembl_installable_genomes = Some(dialog);
    }

    pub(super) fn render_pending_ensembl_quick_install_dialog(&mut self, ctx: &egui::Context) {
        let Some(mut dialog) = self.pending_ensembl_quick_install.clone() else {
            return;
        };
        let mut open = true;
        let spec = crate::egui_compat::ModalWindowSpec::new(
            "Quick Install Ensembl Candidate",
            egui::Id::new((
                "quick_install_ensembl_candidate_modal",
                dialog.preview.collection.clone(),
                dialog.preview.species_dir.clone(),
            )),
        )
        .default_size(Vec2::new(760.0, 520.0))
        .min_size(Vec2::new(560.0, 320.0))
        .resizable(true);
        crate::egui_compat::show_modal_window(ctx, &spec, &mut open, |ui| {
            ui.label(format!(
                "Install '{}' from the current Ensembl {} listing?",
                dialog.preview.display_name, dialog.preview.collection
            ));
            ui.small(
                    "This writes a real catalog entry first and then starts the normal prepare workflow. Downloads and indexing may take some time.",
                );
            ui.separator();
            ui.label(format!("species_dir: {}", dialog.preview.species_dir));
            ui.label(format!("resolved release: {}", dialog.preview.release));
            ui.label(format!("file_stem: {}", dialog.preview.file_stem));
            ui.horizontal(|ui| {
                ui.label("genome_id");
                ui.text_edit_singleline(&mut dialog.genome_id);
            });
            ui.horizontal(|ui| {
                ui.label("output");
                ui.text_edit_singleline(&mut dialog.output_catalog_path);
                if ui
                    .button("Browse...")
                    .on_hover_text("Choose a catalog file or overlay directory for the new entry")
                    .clicked()
                    && let Some(path) = rfd::FileDialog::new().save_file()
                {
                    dialog.output_catalog_path = path.display().to_string();
                }
            });
            ui.small(format!(
                "catalog origin: {}",
                dialog.preview.catalog_origin_label
            ));
            ui.small(format!(
                "catalog write mode: {} ({})",
                dialog.preview.catalog_write_mode, dialog.preview.catalog_entry_action
            ));
            ui.separator();
            ui.small(format!("sequence: {}", dialog.preview.sequence_remote));
            ui.small(format!("annotation: {}", dialog.preview.annotations_remote));
            if !dialog.preview.warnings.is_empty() {
                ui.separator();
                ui.label("Warnings:");
                for warning in &dialog.preview.warnings {
                    ui.small(warning);
                }
            }
            if self.genome_prepare_task.is_some() {
                ui.separator();
                ui.colored_label(
                        egui::Color32::from_rgb(190, 70, 70),
                        "Another genome prepare task is already running. Finish or cancel it before starting a quick install.",
                    );
            }
            ui.separator();
            ui.horizontal(|ui| {
                    if ui
                        .add_enabled(
                            self.genome_prepare_task.is_none(),
                            egui::Button::new("Install + Prepare"),
                        )
                        .on_hover_text(
                            "Write the catalog entry using the shown output path and launch the prepare workflow.",
                        )
                        .clicked()
                    {
                        self.pending_ensembl_quick_install = Some(dialog.clone());
                        self.apply_pending_ensembl_quick_install();
                    }
                    if ui
                        .button("Cancel")
                        .on_hover_text("Close without writing a catalog entry")
                        .clicked()
                    {
                        self.pending_ensembl_quick_install = None;
                    }
                });
            if self.pending_ensembl_quick_install.is_some() {
                self.pending_ensembl_quick_install = Some(dialog);
            }
        });
        if !open {
            self.pending_ensembl_quick_install = None;
        }
    }

    pub(super) fn render_pending_prepared_genome_removal_dialog(&mut self, ctx: &egui::Context) {
        let Some(request) = self.pending_prepared_genome_removal.clone() else {
            return;
        };
        let mut open = true;
        let spec = crate::egui_compat::ModalWindowSpec::new(
            "Remove Prepared Genome",
            egui::Id::new((
                "remove_prepared_genome_modal",
                matches!(request.scope, GenomeDialogScope::Helper),
                request.genome_id.clone(),
            )),
        );
        crate::egui_compat::show_modal_window(ctx, &spec, &mut open, |ui| {
            ui.label(format!(
                "Remove prepared {} '{}' from the cache?",
                request.scope.description(),
                request.genome_id
            ));
            ui.small("This deletes the prepared install only. The catalog entry stays available.");
            ui.separator();
            ui.label(format!("catalog: {}", request.catalog_path));
            let cache_display = request
                .cache_dir
                .as_deref()
                .filter(|value| !value.trim().is_empty())
                .unwrap_or("(catalog/default)");
            ui.label(format!("cache_dir: {cache_display}"));
            ui.label(format!("install_dir: {}", request.install_dir));
            ui.horizontal(|ui| {
                    if ui
                        .button("Remove Prepared")
                        .on_hover_text("Delete the prepared install directory and all cached/indexed files for this genome")
                        .clicked()
                    {
                        self.confirm_prepared_genome_removal();
                    }
                    if ui
                        .button("Cancel")
                        .on_hover_text("Keep the prepared install unchanged")
                        .clicked()
                    {
                        self.pending_prepared_genome_removal = None;
                    }
                });
        });
        if !open {
            self.pending_prepared_genome_removal = None;
        }
    }

    pub(super) fn render_pending_catalog_entry_removal_dialog(&mut self, ctx: &egui::Context) {
        let Some(request) = self.pending_genome_catalog_entry_removal.clone() else {
            return;
        };
        let mut open = true;
        let spec = crate::egui_compat::ModalWindowSpec::new(
            "Remove Catalog Entry",
            egui::Id::new((
                "remove_catalog_entry_modal",
                matches!(request.scope, GenomeDialogScope::Helper),
                request.genome_id.clone(),
            )),
        );
        crate::egui_compat::show_modal_window(ctx, &spec, &mut open, |ui| {
            ui.label(format!(
                "Remove {} catalog entry '{}'?",
                request.scope.description(),
                request.genome_id
            ));
            ui.small(
                    "This edits only the catalog JSON. Prepared cache files, if any, are not deleted unless you also use 'Remove Prepared...'.",
                );
            ui.separator();
            ui.label(format!("catalog: {}", request.catalog_path));
            ui.horizontal(|ui| {
                if ui
                    .button("Remove Catalog Entry")
                    .on_hover_text("Delete this genome row from the active writable catalog")
                    .clicked()
                {
                    self.confirm_catalog_entry_removal();
                }
                if ui
                    .button("Cancel")
                    .on_hover_text("Keep the catalog entry unchanged")
                    .clicked()
                {
                    self.pending_genome_catalog_entry_removal = None;
                }
            });
        });
        if !open {
            self.pending_genome_catalog_entry_removal = None;
        }
    }

    pub(super) fn render_reference_genome_retrieve_contents(&mut self, ui: &mut Ui) -> bool {
        let mut close_requested = false;
        ui.push_id("retrieve_genome_dialog_contents", |ui| {
            let scope = self.genome_dialog_scope;
            let close_hover = Self::specialist_window_close_hover_text(scope.retrieve_title());
            if self
                .render_specialist_window_nav_with_close(ui, Some(("Close", close_hover.as_str())))
            {
                close_requested = true;
            }
            ui.label(format!(
                "Retrieve region sequences from a prepared {}.",
                scope.description()
            ));
            ui.horizontal(|ui| {
                ui.label("catalog");
                ui.text_edit_singleline(&mut self.genome_catalog_path);
                if ui
                    .button("Browse...")
                    .on_hover_text("Browse filesystem and fill this path")
                    .clicked()
                    && let Some(path) = rfd::FileDialog::new()
                        .add_filter("JSON", &["json"])
                        .pick_file()
                    {
                        self.genome_catalog_path = path.display().to_string();
                    }
            });
            ui.horizontal(|ui| {
                ui.label("cache_dir");
                ui.text_edit_singleline(&mut self.genome_cache_dir);
                if ui
                    .button("Browse...")
                    .on_hover_text("Browse filesystem and fill this path")
                    .clicked()
                    && let Some(path) = rfd::FileDialog::new().pick_folder() {
                        self.genome_cache_dir = path.display().to_string();
                    }
            });
            if !self.genome_catalog_error.is_empty() {
                ui.colored_label(
                    egui::Color32::from_rgb(190, 70, 70),
                    format!("Catalog error: {}", self.genome_catalog_error),
                );
            }
            let prepared_genomes = match self.prepared_genomes_for_retrieve_dialog() {
                Ok(names) => names,
                Err(e) => {
                    ui.colored_label(
                        egui::Color32::from_rgb(190, 70, 70),
                        format!("Prepared-state check error: {e}"),
                    );
                    vec![]
                }
            };
            if !prepared_genomes.is_empty() && !prepared_genomes.contains(&self.genome_id) {
                self.genome_id = prepared_genomes[0].clone();
                self.invalidate_genome_genes();
            }
            let selection_changed =
                Self::choose_genome_from_catalog(ui, &mut self.genome_id, &prepared_genomes);
            if selection_changed {
                self.invalidate_genome_genes();
            }
            match self.selected_genome_source_plan() {
                Ok(Some(plan)) => {
                    ui.small(Self::format_genome_source_plan_summary(&plan));
                }
                Ok(None) => {}
                Err(e) => {
                    ui.colored_label(
                        egui::Color32::from_rgb(190, 70, 70),
                        format!("Source-plan error: {e}"),
                    );
                }
            }
            if matches!(scope, GenomeDialogScope::Helper) {
                match self.selected_helper_catalog_entry_for_scope(scope) {
                    Ok(Some(entry)) => {
                        Self::render_helper_catalog_entry_panel(
                            ui,
                            &entry,
                            "Helper Construct Interpretation",
                        );
                    }
                    Ok(None) => {}
                    Err(e) => {
                        ui.colored_label(
                            egui::Color32::from_rgb(190, 70, 70),
                            format!("Helper interpretation error: {e}"),
                        );
                    }
                }
            }
            if prepared_genomes.is_empty() {
                ui.label(scope.empty_prepared_hint());
            } else {
                self.ensure_genome_genes_loaded();
                if !self.genome_genes_error.is_empty() {
                    ui.colored_label(
                        egui::Color32::from_rgb(190, 70, 70),
                        format!("Gene index error: {}", self.genome_genes_error),
                    );
                }
                ui.separator();
                ui.horizontal(|ui| {
                    ui.label("Gene filter");
                    let response = ui.text_edit_singleline(&mut self.genome_gene_filter);
                    if response.changed() {
                        self.genome_gene_filter_page = 0;
                    }
                    if ui
                        .button("Clear")
                        .on_hover_text("Clear the gene filter text")
                        .clicked()
                    {
                        self.genome_gene_filter.clear();
                        self.genome_gene_filter_page = 0;
                    }
                });
                ui.small("Supports case-insensitive regex (for example: ^TP53$).");
                ui.horizontal(|ui| {
                    ui.label("Top matches");
                    if ui
                        .add(
                            egui::DragValue::new(&mut self.genome_gene_filter_limit)
                                .speed(100.0)
                                .range(100..=100_000),
                        )
                        .changed()
                    {
                        self.genome_gene_filter_page = 0;
                    }
                });
                if !self.genome_biotype_filter.is_empty() {
                    ui.separator();
                    ui.label("Biotype filter");
                    ui.horizontal(|ui| {
                        if ui
                            .button("All")
                            .on_hover_text("Enable all biotypes")
                            .clicked()
                        {
                            for enabled in self.genome_biotype_filter.values_mut() {
                                *enabled = true;
                            }
                            self.genome_gene_filter_page = 0;
                        }
                        if ui
                            .button("None")
                            .on_hover_text("Disable all biotypes")
                            .clicked()
                        {
                            for enabled in self.genome_biotype_filter.values_mut() {
                                *enabled = false;
                            }
                            self.genome_gene_filter_page = 0;
                        }
                    });
                    let mut biotypes: Vec<String> =
                        self.genome_biotype_filter.keys().cloned().collect();
                    biotypes.sort();
                    egui::ScrollArea::vertical()
                        .id_salt("retrieve_genome_biotype_scroll")
                        .max_height(120.0)
                        .show(ui, |ui| {
                            scroll_input_policy::apply_scrollarea_keyboard_navigation(
                                ui,
                                scroll_input_policy::DEFAULT_SCROLLAREA_KEYBOARD_STEP,
                            );
                            for biotype in biotypes {
                                if let Some(enabled) = self.genome_biotype_filter.get_mut(&biotype)
                                    && ui.checkbox(enabled, &biotype).changed() {
                                        self.genome_gene_filter_page = 0;
                                    }
                            }
                        });
                }
                let (filtered_indices, overflow, regex_error) = self.filtered_gene_candidate_indices();
                if let Some(error) = regex_error {
                    ui.colored_label(egui::Color32::from_rgb(190, 70, 70), error);
                }
                const GENE_PAGE_SIZE: usize = 100;
                let page_count =
                    filtered_indices.len().div_ceil(GENE_PAGE_SIZE).max(1);
                if self.genome_gene_filter_page >= page_count {
                    self.genome_gene_filter_page = page_count.saturating_sub(1);
                }
                let page_start = self.genome_gene_filter_page * GENE_PAGE_SIZE;
                let page_end = (page_start + GENE_PAGE_SIZE).min(filtered_indices.len());
                let shown_note = if overflow {
                    format!(" (limited to top {} matches)", self.genome_gene_filter_limit)
                } else {
                    String::new()
                };
                ui.label(format!(
                    "Genes: {} candidate{} / {} total",
                    filtered_indices.len(),
                    shown_note,
                    self.genome_genes.len()
                ));
                ui.horizontal(|ui| {
                    if ui
                        .add_enabled(
                            self.genome_gene_filter_page > 0,
                            egui::Button::new("Prev"),
                        )
                        .on_hover_text("Show previous page of gene matches")
                        .clicked()
                    {
                        self.genome_gene_filter_page -= 1;
                    }
                    ui.label(format!(
                        "Page {}/{}",
                        self.genome_gene_filter_page + 1,
                        page_count
                    ));
                    if ui
                        .add_enabled(
                            self.genome_gene_filter_page + 1 < page_count,
                            egui::Button::new("Next"),
                        )
                        .on_hover_text("Show next page of gene matches")
                        .clicked()
                    {
                        self.genome_gene_filter_page += 1;
                    }
                });
                egui::ScrollArea::vertical()
                    .id_salt("retrieve_genome_gene_list_scroll")
                    .max_height(220.0)
                    .show(ui, |ui| {
                        scroll_input_policy::apply_scrollarea_keyboard_navigation(
                            ui,
                            scroll_input_policy::DEFAULT_SCROLLAREA_KEYBOARD_STEP,
                        );
                        for idx in filtered_indices
                            .iter()
                            .copied()
                            .skip(page_start)
                            .take(page_end.saturating_sub(page_start))
                        {
                            let label = Self::gene_record_label(&self.genome_genes[idx]);
                            let selected = self.genome_selected_gene == Some(idx);
                            if ui.selectable_label(selected, label).clicked() {
                                self.select_gene_record(idx);
                            }
                        }
                    });
                ui.separator();
                ui.horizontal(|ui| {
                    ui.label("chr");
                    if ui
                        .text_edit_singleline(&mut self.genome_chromosome)
                        .changed()
                    {
                        self.genome_retrieve_contig_suggestions.clear();
                    }
                    ui.label("start_1based");
                    let start_changed = ui
                        .add(
                            egui::TextEdit::singleline(&mut self.genome_start_1based)
                                .id_salt("retrieve_genome_start_1based_field")
                                .desired_width(120.0)
                                .char_limit(10),
                        )
                        .changed();
                    if start_changed {
                        Self::normalize_coordinate_field(&mut self.genome_start_1based);
                    }
                    ui.label("end_1based");
                    let end_changed = ui
                        .add(
                            egui::TextEdit::singleline(&mut self.genome_end_1based)
                                .id_salt("retrieve_genome_end_1based_field")
                                .desired_width(120.0)
                                .char_limit(10),
                        )
                        .changed();
                    if end_changed {
                        Self::normalize_coordinate_field(&mut self.genome_end_1based);
                    }
                });
                ui.horizontal(|ui| {
                    ui.label("selected gene extract");
                    let previous_mode = self.genome_gene_extract_mode;
                    egui::ComboBox::from_id_salt("retrieve_genome_gene_extract_mode")
                        .selected_text(match self.genome_gene_extract_mode {
                            GenomeGeneExtractMode::Gene => "gene span",
                            GenomeGeneExtractMode::CodingWithPromoter => "CDS + promoter",
                        })
                        .show_ui(ui, |ui| {
                            ui.selectable_value(
                                &mut self.genome_gene_extract_mode,
                                GenomeGeneExtractMode::Gene,
                                "gene span",
                            );
                            ui.selectable_value(
                                &mut self.genome_gene_extract_mode,
                                GenomeGeneExtractMode::CodingWithPromoter,
                                "CDS + promoter",
                            );
                        })
                        .response
                        .on_hover_text(
                            "Choose whether Extract Selected Gene uses the full gene span or the CDS span plus an additional 5' promoter flank",
                        );
                    ui.label("promoter bp before CDS");
                    let promoter_changed = ui
                        .add_enabled(
                            matches!(
                                self.genome_gene_extract_mode,
                                GenomeGeneExtractMode::CodingWithPromoter
                            ),
                            egui::TextEdit::singleline(
                                &mut self.genome_gene_promoter_upstream_bp,
                            )
                            .desired_width(120.0)
                            .hint_text("0 = CDS only"),
                        )
                        .on_hover_text(
                            "Additional bases to include on the gene's 5' side before the first coding base (strand-aware)",
                        )
                        .changed();
                    if promoter_changed {
                        self.genome_gene_promoter_upstream_bp
                            .retain(|ch| ch.is_ascii_digit());
                    }
                    if self.genome_gene_extract_mode != previous_mode || promoter_changed {
                        self.refresh_selected_gene_output_id_if_autofilled();
                    }
                });
                ui.small(
                    "Extract Selected Gene uses the mode above. Extract Region always uses the explicit chr/start/end interval.",
                );
                ui.horizontal(|ui| {
                    ui.label("output_id");
                    if ui.text_edit_singleline(&mut self.genome_output_id).changed() {
                        self.genome_output_id_autofilled = false;
                    }
                    if ui
                        .checkbox(
                            &mut self.genome_include_genomic_annotation,
                            "include genomic annotation",
                        )
                        .on_hover_text(
                            "Attach overlapping genomic annotation when extracting selected genes or explicit regions",
                        )
                        .changed()
                    {
                        self.sync_genome_annotation_scope_controls();
                    }
                    if ui
                        .add_enabled(
                            self.genome_selected_gene.is_some(),
                            egui::Button::new("Extract Selected Gene"),
                        )
                        .on_hover_text(
                            "Extract the currently selected gene from the prepared reference",
                        )
                        .clicked()
                    {
                        self.extract_reference_genome_gene();
                    }
                    if ui
                        .button("Extract Region")
                        .on_hover_text("Extract the explicit chromosome start/end interval")
                        .clicked()
                    {
                        self.extract_reference_genome_region();
                    }
                });
                ui.horizontal(|ui| {
                    ui.label("annotation scope");
                    let previous_scope = self.genome_annotation_scope;
                    let scope_combo =
                        egui::ComboBox::from_id_salt("retrieve_genome_annotation_scope")
                            .selected_text(self.genome_annotation_scope.as_str())
                            .show_ui(ui, |ui| {
                                ui.selectable_value(
                                    &mut self.genome_annotation_scope,
                                    GenomeAnnotationScope::None,
                                    "none",
                                );
                                ui.selectable_value(
                                    &mut self.genome_annotation_scope,
                                    GenomeAnnotationScope::Core,
                                    "core",
                                );
                                ui.selectable_value(
                                    &mut self.genome_annotation_scope,
                                    GenomeAnnotationScope::Full,
                                    "full",
                                );
                            });
                    scope_combo.response.on_hover_text(
                        "Projection policy for annotation transfer: none, core gene/transcript context, or full",
                    );
                    if self.genome_annotation_scope != previous_scope {
                        self.genome_include_genomic_annotation =
                            !matches!(self.genome_annotation_scope, GenomeAnnotationScope::None);
                    }
                    ui.label("max annotation features");
                    let cap_changed = ui
                        .add(
                            egui::TextEdit::singleline(&mut self.genome_max_annotation_features)
                                .desired_width(120.0)
                                .hint_text("empty = unlimited"),
                        )
                        .on_hover_text("Optional safety cap. 0 means unlimited explicit cap.")
                        .changed();
                    if cap_changed {
                        self.genome_max_annotation_features
                            .retain(|ch| ch.is_ascii_digit());
                    }
                });
            }
            if !self.genome_retrieve_status.is_empty() {
                ui.separator();
                ui.monospace(&self.genome_retrieve_status);
                if let Some(suggested) = self.genome_retrieve_contig_suggestions.first().cloned() {
                    ui.horizontal_wrapped(|ui| {
                        ui.monospace(format!("Suggested contig: '{suggested}'"));
                        if ui
                            .button("Apply")
                            .on_hover_text(
                                "Apply the best matching prepared contig/chromosome to the chr field",
                            )
                            .clicked()
                        {
                            self.genome_chromosome = suggested.clone();
                            self.genome_retrieve_status = format!(
                                "Applied suggested contig '{}' to chr field. Re-run extraction.",
                                suggested
                            );
                            self.genome_retrieve_contig_suggestions.clear();
                        }
                        if ui
                            .button("Copy")
                            .on_hover_text("Copy the suggested contig/chromosome to the clipboard")
                            .clicked()
                        {
                            ui.ctx().copy_text(suggested.clone());
                        }
                        let alternatives = self
                            .genome_retrieve_contig_suggestions
                            .iter()
                            .skip(1)
                            .take(4)
                            .cloned()
                            .collect::<Vec<_>>();
                        if !alternatives.is_empty() {
                            ui.small(format!("alternatives: {}", alternatives.join(", ")));
                        }
                    });
                }
            }
        });
        close_requested
    }

    pub(super) fn render_reference_genome_retrieve_dialog(&mut self, ctx: &egui::Context) {
        if !self.show_reference_genome_retrieve_dialog {
            return;
        }
        self.refresh_genome_catalog_list();
        let title = self.genome_dialog_scope.retrieve_title();
        let viewport_id = Self::retrieve_genome_viewport_id();
        let spec = self.hosted_window_spec_for_viewport(
            title,
            egui::Id::new(("hosted_retrieve_genome_window", viewport_id)),
            viewport_id,
            Vec2::new(980.0, 760.0),
            Vec2::new(860.0, 520.0),
        );

        if ctx.embed_viewports() {
            let mut open = self.show_reference_genome_retrieve_dialog;
            let mut close_requested = false;
            crate::egui_compat::show_hosted_window(ctx, &spec, &mut open, |ui| {
                close_requested = self.render_reference_genome_retrieve_contents(ui);
            });
            self.clear_viewport_foreground_request_after_render(viewport_id);
            if close_requested {
                open = false;
            }
            self.show_reference_genome_retrieve_dialog = open;
            return;
        }

        let builder = crate::egui_compat::viewport_builder_for_hosted_window(&spec);
        ctx.show_viewport_immediate(viewport_id, builder, |ctx, class| {
            self.note_viewport_focus_if_active(ctx, viewport_id);
            if class == egui::ViewportClass::EmbeddedWindow {
                let mut open = self.show_reference_genome_retrieve_dialog;
                let mut close_requested = false;
                crate::egui_compat::show_hosted_window(ctx, &spec, &mut open, |ui| {
                    close_requested = self.render_reference_genome_retrieve_contents(ui);
                });
                if close_requested {
                    open = false;
                }
                self.show_reference_genome_retrieve_dialog = open;
                return;
            }

            let mut close_requested = false;
            crate::egui_compat::show_central_panel(ctx, egui::CentralPanel::default(), |ui| {
                close_requested = self.render_reference_genome_retrieve_contents(ui);
            });

            if close_requested || Self::viewport_close_requested_or_shortcut(ctx) {
                self.show_reference_genome_retrieve_dialog = false;
            }
        });
    }

    pub(super) fn format_blast_command_line(executable: &str, args: &[String]) -> String {
        let exe = executable.trim();
        if args.is_empty() {
            return exe.to_string();
        }
        if exe.is_empty() {
            return args.join(" ");
        }
        format!("{exe} {}", args.join(" "))
    }

    pub(super) fn render_blast_query_result(ui: &mut Ui, result: &GenomeBlastQueryResult) {
        let report = &result.report;
        ui.label(format!(
            "Query: {} ({} bp) | Genome: {} | Hits: {}",
            result.query_label, result.query_length, report.genome_id, report.hit_count
        ));
        ui.label(format!(
            "Task: {} | max_hits: {}",
            report.task, report.max_hits
        ));
        ui.monospace(format!("blastn: {}", report.blastn_executable));
        ui.monospace(format!("db: {}", report.blast_db_prefix));
        let invocation =
            Self::format_blast_command_line(&report.blastn_executable, &report.command);
        if !invocation.trim().is_empty() {
            ui.monospace(format!("invocation: {invocation}"));
        }
        if let Some(request) = report.options_override_json.as_ref()
            && let Ok(compact) = serde_json::to_string(request)
        {
            ui.monospace(format!("request_options: {compact}"));
        }
        if let Some(effective) = report.effective_options_json.as_ref()
            && let Ok(compact) = serde_json::to_string(effective)
        {
            ui.monospace(format!("effective_options: {compact}"));
        }
        if !report.warnings.is_empty() {
            ui.separator();
            ui.label("Warnings");
            for warning in &report.warnings {
                ui.monospace(format!("- {}", warning));
            }
        }
        if !report.stderr.trim().is_empty() {
            ui.separator();
            ui.label("BLAST stderr");
            ui.monospace(report.stderr.trim().to_string());
        }
        ui.separator();
        ui.label("Hits");
        let max_rows = 200usize;
        let shown_hits = report.hits.len().min(max_rows);
        if shown_hits < report.hits.len() {
            ui.small(format!(
                "Showing first {} of {} hits",
                shown_hits,
                report.hits.len()
            ));
        }
        egui::ScrollArea::both().max_height(280.0).show(ui, |ui| {
            scroll_input_policy::apply_scrollarea_keyboard_navigation(
                ui,
                scroll_input_policy::DEFAULT_SCROLLAREA_KEYBOARD_STEP,
            );
            egui::Grid::new(format!("blast_hits_{}", result.query_label))
                .striped(true)
                .show(ui, |ui| {
                    ui.strong("subject");
                    ui.strong("%identity");
                    ui.strong("aln_len");
                    ui.strong("q_range");
                    ui.strong("s_range");
                    ui.strong("evalue");
                    ui.strong("bit_score");
                    ui.strong("qcov%");
                    ui.end_row();
                    for hit in report.hits.iter().take(max_rows) {
                        ui.monospace(hit.subject_id.clone());
                        ui.label(format!("{:.2}", hit.identity_percent));
                        ui.label(hit.alignment_length.to_string());
                        ui.label(format!("{}-{}", hit.query_start, hit.query_end));
                        ui.label(format!("{}-{}", hit.subject_start, hit.subject_end));
                        ui.label(format!("{:.3e}", hit.evalue));
                        ui.label(format!("{:.2}", hit.bit_score));
                        ui.label(
                            hit.query_coverage_percent
                                .map(|v| format!("{v:.1}"))
                                .unwrap_or_else(|| "-".to_string()),
                        );
                        ui.end_row();
                    }
                });
        });
    }

    pub(super) fn render_reference_genome_blast_contents(&mut self, ui: &mut Ui) -> bool {
        self.refresh_genome_catalog_list();
        let scope = self.genome_dialog_scope;
        let mut close_requested = false;
        let close_hover = Self::specialist_window_close_hover_text(scope.blast_title());
        if self.render_specialist_window_nav_with_close(ui, Some(("Close", close_hover.as_str()))) {
            close_requested = true;
        }
        let sequence_ids = self.project_sequence_ids_for_blast();
        if !sequence_ids.is_empty() && !sequence_ids.contains(&self.genome_blast_query_seq_id) {
            self.genome_blast_query_seq_id = sequence_ids[0].clone();
        }
        let pool_options = self.blast_pool_options();
        if !pool_options
            .iter()
            .any(|pool| pool.container_id == self.genome_blast_query_pool_id)
        {
            self.genome_blast_query_pool_id = pool_options
                .first()
                .map(|pool| pool.container_id.clone())
                .unwrap_or_default();
        }

        ui.label(format!(
            "Run BLAST against a prepared {} index.",
            scope.description()
        ));
        ui.separator();
        ui.strong("1. Target");
        ui.horizontal(|ui| {
            ui.label("catalog");
            ui.text_edit_singleline(&mut self.genome_catalog_path);
            if ui
                .button("Browse...")
                .on_hover_text("Browse filesystem and fill this path")
                .clicked()
                && let Some(path) = rfd::FileDialog::new()
                    .add_filter("JSON", &["json"])
                    .pick_file()
            {
                self.genome_catalog_path = path.display().to_string();
            }
        });
        ui.horizontal(|ui| {
            ui.label("cache_dir");
            ui.text_edit_singleline(&mut self.genome_cache_dir);
            if ui
                .button("Browse...")
                .on_hover_text("Browse filesystem and fill this path")
                .clicked()
                && let Some(path) = rfd::FileDialog::new().pick_folder()
            {
                self.genome_cache_dir = path.display().to_string();
            }
        });
        if !self.genome_catalog_error.is_empty() {
            ui.colored_label(
                egui::Color32::from_rgb(190, 70, 70),
                format!("Catalog error: {}", self.genome_catalog_error),
            );
        }

        let prepared_genomes = match self.prepared_genomes_for_retrieve_dialog() {
            Ok(names) => names,
            Err(e) => {
                ui.colored_label(
                    egui::Color32::from_rgb(190, 70, 70),
                    format!("Prepared-state check error: {e}"),
                );
                vec![]
            }
        };
        if !prepared_genomes.is_empty() && !prepared_genomes.contains(&self.genome_id) {
            self.genome_id = prepared_genomes[0].clone();
            self.invalidate_genome_genes();
        }
        let selection_changed =
            Self::choose_genome_from_catalog(ui, &mut self.genome_id, &prepared_genomes);
        if selection_changed {
            self.invalidate_genome_genes();
        }
        match self.selected_genome_source_plan() {
            Ok(Some(plan)) => {
                ui.small(Self::format_genome_source_plan_summary(&plan));
            }
            Ok(None) => {}
            Err(e) => {
                ui.colored_label(
                    egui::Color32::from_rgb(190, 70, 70),
                    format!("Source-plan error: {e}"),
                );
            }
        }
        if matches!(scope, GenomeDialogScope::Helper) {
            match self.selected_helper_catalog_entry_for_scope(scope) {
                Ok(Some(entry)) => {
                    Self::render_helper_catalog_entry_panel(
                        ui,
                        &entry,
                        "Helper Construct Interpretation",
                    );
                }
                Ok(None) => {}
                Err(e) => {
                    ui.colored_label(
                        egui::Color32::from_rgb(190, 70, 70),
                        format!("Helper interpretation error: {e}"),
                    );
                }
            }
        }
        if prepared_genomes.is_empty() {
            ui.label(scope.empty_prepared_hint());
        }

        ui.separator();
        ui.strong("2. Input");
        ui.label("Query source");
        ui.horizontal(|ui| {
            ui.selectable_value(
                &mut self.genome_blast_source_mode,
                GenomeBlastSourceMode::Manual,
                "Manual sequence",
            );
            ui.selectable_value(
                &mut self.genome_blast_source_mode,
                GenomeBlastSourceMode::ProjectSequence,
                "Project sequence",
            );
            ui.selectable_value(
                &mut self.genome_blast_source_mode,
                GenomeBlastSourceMode::ProjectPool,
                "Project pool",
            );
        });
        match self.genome_blast_source_mode {
            GenomeBlastSourceMode::Manual => {
                ui.label("Query sequence (IUPAC DNA letters)");
                ui.add(
                    egui::TextEdit::multiline(&mut self.genome_blast_query_manual)
                        .desired_rows(4)
                        .hint_text("Paste query sequence here"),
                );
            }
            GenomeBlastSourceMode::ProjectSequence => {
                if sequence_ids.is_empty() {
                    ui.label("No sequences are loaded in this project.");
                } else {
                    egui::ComboBox::from_label("sequence")
                        .selected_text(if self.genome_blast_query_seq_id.trim().is_empty() {
                            "(choose sequence)"
                        } else {
                            self.genome_blast_query_seq_id.as_str()
                        })
                        .show_ui(ui, |ui| {
                            for seq_id in &sequence_ids {
                                ui.selectable_value(
                                    &mut self.genome_blast_query_seq_id,
                                    seq_id.clone(),
                                    seq_id,
                                );
                            }
                        });
                }
            }
            GenomeBlastSourceMode::ProjectPool => {
                if pool_options.is_empty() {
                    ui.label("No sequence pools/containers with >1 sequence are available.");
                } else {
                    egui::ComboBox::from_label("pool/container")
                        .selected_text(if self.genome_blast_query_pool_id.trim().is_empty() {
                            "(choose pool)"
                        } else {
                            self.genome_blast_query_pool_id.as_str()
                        })
                        .show_ui(ui, |ui| {
                            for pool in &pool_options {
                                ui.selectable_value(
                                    &mut self.genome_blast_query_pool_id,
                                    pool.container_id.clone(),
                                    pool.label.clone(),
                                );
                            }
                        });
                    if let Some(pool) = pool_options
                        .iter()
                        .find(|pool| pool.container_id == self.genome_blast_query_pool_id)
                    {
                        ui.small(format!(
                            "Selected pool will run {} query sequence(s).",
                            pool.members.len()
                        ));
                    }
                }
            }
        }

        ui.separator();
        ui.strong("3. Options");
        ui.horizontal(|ui| {
            ui.label("max_hits");
            ui.add(
                egui::DragValue::new(&mut self.genome_blast_max_hits)
                    .range(1..=1000)
                    .speed(1.0),
            );
            ui.label("task");
            egui::ComboBox::from_id_salt("genome_blast_task_combo")
                .selected_text(self.genome_blast_task_name.clone())
                .show_ui(ui, |ui| {
                    ui.selectable_value(
                        &mut self.genome_blast_task_name,
                        "blastn-short".to_string(),
                        "blastn-short",
                    );
                    ui.selectable_value(
                        &mut self.genome_blast_task_name,
                        "blastn".to_string(),
                        "blastn",
                    );
                });
        });
        ui.horizontal(|ui| {
            ui.label("preset");
            egui::ComboBox::from_id_salt("genome_blast_options_preset_combo")
                .selected_text(self.genome_blast_options_preset.label())
                .show_ui(ui, |ui| {
                    ui.selectable_value(
                        &mut self.genome_blast_options_preset,
                        GenomeBlastOptionsPreset::None,
                        GenomeBlastOptionsPreset::None.label(),
                    );
                    ui.selectable_value(
                        &mut self.genome_blast_options_preset,
                        GenomeBlastOptionsPreset::StrictIdentityCoverage,
                        GenomeBlastOptionsPreset::StrictIdentityCoverage.label(),
                    );
                    ui.selectable_value(
                        &mut self.genome_blast_options_preset,
                        GenomeBlastOptionsPreset::UniqueBestHit,
                        GenomeBlastOptionsPreset::UniqueBestHit.label(),
                    );
                    ui.selectable_value(
                        &mut self.genome_blast_options_preset,
                        GenomeBlastOptionsPreset::HighStringency,
                        GenomeBlastOptionsPreset::HighStringency.label(),
                    );
                });
            if ui
                .button("Load JSON...")
                .on_hover_text("Load advanced BLAST options JSON from file")
                .clicked()
                && let Some(path) = rfd::FileDialog::new()
                    .add_filter("JSON", &["json"])
                    .pick_file()
            {
                match fs::read_to_string(&path) {
                    Ok(text) => {
                        self.genome_blast_options_json = text;
                    }
                    Err(e) => {
                        self.genome_blast_status = format!(
                            "Could not read BLAST options file '{}': {}",
                            path.display(),
                            e
                        );
                    }
                }
            }
        });
        ui.group(|ui| {
            ui.label("Structured threshold controls");
            egui::Grid::new("genome_blast_structured_threshold_grid")
                .num_columns(3)
                .spacing([8.0, 4.0])
                .show(ui, |ui| {
                    ui.checkbox(
                        &mut self.genome_blast_threshold_use_max_evalue,
                        "max_evalue",
                    );
                    ui.add_enabled(
                        self.genome_blast_threshold_use_max_evalue,
                        egui::TextEdit::singleline(&mut self.genome_blast_threshold_max_evalue)
                            .desired_width(140.0)
                            .hint_text("1e-6"),
                    );
                    ui.small("finite >= 0");
                    ui.end_row();

                    ui.checkbox(
                        &mut self.genome_blast_threshold_use_min_identity_percent,
                        "min_identity_percent",
                    );
                    ui.add_enabled(
                        self.genome_blast_threshold_use_min_identity_percent,
                        egui::TextEdit::singleline(
                            &mut self.genome_blast_threshold_min_identity_percent,
                        )
                        .desired_width(140.0)
                        .hint_text("97.0"),
                    );
                    ui.small("0..100");
                    ui.end_row();

                    ui.checkbox(
                        &mut self.genome_blast_threshold_use_min_query_coverage_percent,
                        "min_query_coverage_percent",
                    );
                    ui.add_enabled(
                        self.genome_blast_threshold_use_min_query_coverage_percent,
                        egui::TextEdit::singleline(
                            &mut self.genome_blast_threshold_min_query_coverage_percent,
                        )
                        .desired_width(140.0)
                        .hint_text("80.0"),
                    );
                    ui.small("0..100");
                    ui.end_row();

                    ui.checkbox(
                        &mut self.genome_blast_threshold_use_min_alignment_length_bp,
                        "min_alignment_length_bp",
                    );
                    ui.add_enabled(
                        self.genome_blast_threshold_use_min_alignment_length_bp,
                        egui::TextEdit::singleline(
                            &mut self.genome_blast_threshold_min_alignment_length_bp,
                        )
                        .desired_width(140.0)
                        .hint_text("20"),
                    );
                    ui.small(">= 1");
                    ui.end_row();

                    ui.checkbox(
                        &mut self.genome_blast_threshold_use_min_bit_score,
                        "min_bit_score",
                    );
                    ui.add_enabled(
                        self.genome_blast_threshold_use_min_bit_score,
                        egui::TextEdit::singleline(&mut self.genome_blast_threshold_min_bit_score)
                            .desired_width(140.0)
                            .hint_text("40.0"),
                    );
                    ui.small("finite >= 0");
                    ui.end_row();

                    ui.checkbox(
                        &mut self.genome_blast_threshold_unique_best_hit,
                        "unique_best_hit",
                    );
                    ui.small("-");
                    ui.small("require exactly one hit after filtering");
                    ui.end_row();
                });
            ui.horizontal(|ui| {
                if ui
                    .button("Reset thresholds")
                    .on_hover_text("Clear all structured threshold controls")
                    .clicked()
                {
                    self.genome_blast_threshold_use_max_evalue = false;
                    self.genome_blast_threshold_max_evalue.clear();
                    self.genome_blast_threshold_use_min_identity_percent = false;
                    self.genome_blast_threshold_min_identity_percent.clear();
                    self.genome_blast_threshold_use_min_query_coverage_percent = false;
                    self.genome_blast_threshold_min_query_coverage_percent
                        .clear();
                    self.genome_blast_threshold_use_min_alignment_length_bp = false;
                    self.genome_blast_threshold_min_alignment_length_bp.clear();
                    self.genome_blast_threshold_use_min_bit_score = false;
                    self.genome_blast_threshold_min_bit_score.clear();
                    self.genome_blast_threshold_unique_best_hit = false;
                }
            });
        });
        ui.label("advanced options JSON (object; merged after preset + structured controls)");
        ui.add(
            egui::TextEdit::multiline(&mut self.genome_blast_options_json)
                .desired_rows(4)
                .hint_text("{\"thresholds\":{\"min_identity_percent\":97.0}}"),
        );
        match self.build_genome_blast_request_override_json() {
            Ok(request_override_json) => {
                if let Some(v) = request_override_json.as_ref() {
                    match self.resolve_genome_blast_options_preview(Some(v)) {
                        Ok(resolved) => {
                            ui.small(format!(
                                "Effective options preview: task='{}', max_hits={}, thresholds={}",
                                resolved.task,
                                resolved.max_hits,
                                Self::format_blast_thresholds_summary(&resolved.thresholds)
                            ));
                        }
                        Err(e) => {
                            ui.colored_label(
                                egui::Color32::from_rgb(190, 70, 70),
                                format!("Options validation error: {e}"),
                            );
                        }
                    }
                } else {
                    match self.resolve_genome_blast_options_preview(None) {
                        Ok(resolved) => {
                            ui.small(format!(
                                "Effective options preview: task='{}', max_hits={}, thresholds={}",
                                resolved.task,
                                resolved.max_hits,
                                Self::format_blast_thresholds_summary(&resolved.thresholds)
                            ));
                        }
                        Err(e) => {
                            ui.colored_label(
                                egui::Color32::from_rgb(190, 70, 70),
                                format!("Options validation error: {e}"),
                            );
                        }
                    }
                }
            }
            Err(e) => {
                ui.colored_label(
                    egui::Color32::from_rgb(190, 70, 70),
                    format!("Options JSON error: {e}"),
                );
            }
        }

        ui.separator();
        ui.strong("4. Execution");
        let running = self.genome_blast_task.is_some();
        ui.horizontal(|ui| {
            if ui
                .add_enabled(
                    !running && !prepared_genomes.is_empty(),
                    egui::Button::new("Run BLAST"),
                )
                .on_hover_text("Run BLAST query and collect hits for the selected prepared genome")
                .clicked()
            {
                self.start_reference_genome_blast();
            }
            if ui
                .add_enabled(running, egui::Button::new("Cancel BLAST"))
                .on_hover_text("Request cancellation of the running BLAST task")
                .clicked()
            {
                self.request_blast_task_cancel("blast dialog");
            }
            if ui
                .button("Clear Results")
                .on_hover_text("Clear all BLAST query results from this dialog")
                .clicked()
            {
                self.genome_blast_results.clear();
                self.genome_blast_selected_result = 0;
            }
        });

        if let Some(fraction) = self.genome_blast_progress_fraction {
            ui.add(
                egui::ProgressBar::new(fraction.clamp(0.0, 1.0))
                    .show_percentage()
                    .text(self.genome_blast_progress_label.clone()),
            );
        }
        if !self.genome_blast_status.is_empty() {
            ui.separator();
            ui.monospace(self.genome_blast_status.clone());
        }
        if !self.genome_blast_results.is_empty() {
            ui.separator();
            ui.strong("5. Results");
            if self.genome_blast_selected_result >= self.genome_blast_results.len() {
                self.genome_blast_selected_result = 0;
            }
            if self.genome_blast_results.len() > 1 {
                egui::ComboBox::from_label("Query result")
                    .selected_text(
                        self.genome_blast_results[self.genome_blast_selected_result]
                            .query_label
                            .clone(),
                    )
                    .show_ui(ui, |ui| {
                        for (idx, result) in self.genome_blast_results.iter().enumerate() {
                            let label = format!(
                                "{} (hits={})",
                                result.query_label, result.report.hit_count
                            );
                            ui.selectable_value(&mut self.genome_blast_selected_result, idx, label);
                        }
                    });
            }
            if let Some(result) = self
                .genome_blast_results
                .get(self.genome_blast_selected_result)
                .cloned()
            {
                let target_seq_id = self.target_seq_id_for_blast_result(&result);
                ui.separator();
                ui.label("Import BLAST hits as features");
                ui.horizontal(|ui| {
                    ui.label("track_name");
                    ui.text_edit_singleline(&mut self.genome_blast_import_track_name);
                    ui.checkbox(
                        &mut self.genome_blast_import_clear_existing,
                        "Clear existing BLAST-hit features first",
                    );
                });
                let can_import = target_seq_id.is_some() && !result.report.hits.is_empty();
                let button = ui
                    .add_enabled(can_import, egui::Button::new("Import Hits To Sequence"))
                    .on_hover_text(
                        "Import current BLAST hit intervals as features on the target sequence",
                    );
                if button.clicked() {
                    self.import_selected_blast_hits_as_track();
                }
                if let Some(seq_id) = target_seq_id {
                    ui.small(format!("Import target sequence: {seq_id}"));
                } else {
                    ui.small(
                        "Hit import requires query label to match an existing project sequence ID.",
                    );
                }
                Self::render_blast_query_result(ui, &result);
            }
        }
        close_requested
    }

    pub(super) fn render_reference_genome_blast_dialog(&mut self, ctx: &egui::Context) {
        if !self.show_reference_genome_blast_dialog {
            return;
        }
        let title = self.genome_dialog_scope.blast_title();
        let viewport_id = Self::blast_genome_viewport_id();
        let spec = self.hosted_window_spec_for_viewport(
            title,
            egui::Id::new(("hosted_blast_genome_window", viewport_id)),
            viewport_id,
            Vec2::new(980.0, 700.0),
            Vec2::new(640.0, 420.0),
        );

        if ctx.embed_viewports() {
            let mut open = self.show_reference_genome_blast_dialog;
            let mut close_requested = false;
            crate::egui_compat::show_hosted_window(ctx, &spec, &mut open, |ui| {
                close_requested = self.render_reference_genome_blast_contents(ui);
            });
            self.clear_viewport_foreground_request_after_render(viewport_id);
            if close_requested {
                open = false;
            }
            self.show_reference_genome_blast_dialog = open;
            return;
        }

        let builder = crate::egui_compat::viewport_builder_for_hosted_window(&spec);
        ctx.show_viewport_immediate(viewport_id, builder, |ctx, class| {
            self.note_viewport_focus_if_active(ctx, viewport_id);
            if class == egui::ViewportClass::EmbeddedWindow {
                let mut open = self.show_reference_genome_blast_dialog;
                let mut close_requested = false;
                crate::egui_compat::show_hosted_window(ctx, &spec, &mut open, |ui| {
                    close_requested = self.render_reference_genome_blast_contents(ui);
                });
                if close_requested {
                    open = false;
                }
                self.show_reference_genome_blast_dialog = open;
                return;
            }

            let mut close_requested = false;
            crate::egui_compat::show_central_panel(ctx, egui::CentralPanel::default(), |ui| {
                close_requested = self.render_reference_genome_blast_contents(ui);
            });

            if close_requested || Self::viewport_close_requested_or_shortcut(ctx) {
                self.show_reference_genome_blast_dialog = false;
            }
        });
    }

    pub(super) fn render_reference_genome_inspector_dialog(&mut self, ctx: &egui::Context) {
        if !self.show_reference_genome_inspector_dialog {
            return;
        }
        self.refresh_genome_catalog_list();
        let scope = self.genome_dialog_scope;
        let title = if matches!(scope, GenomeDialogScope::Helper) {
            "Prepared Helper Genomes"
        } else {
            "Prepared Genome References"
        };
        let mut open = self.show_reference_genome_inspector_dialog;
        let mut close_requested = false;
        let spec = crate::egui_compat::HostedWindowSpec::new(
            title,
            egui::Id::new((
                "hosted_reference_genome_inspector_window",
                matches!(scope, GenomeDialogScope::Helper),
            )),
            Vec2::new(900.0, 620.0),
            Vec2::new(620.0, 360.0),
        );
        crate::egui_compat::show_hosted_window(ctx, &spec, &mut open, |ui| {
            let close_hover = Self::specialist_window_close_hover_text(title);
            if self
                .render_specialist_window_nav_with_close(ui, Some(("Close", close_hover.as_str())))
            {
                close_requested = true;
            }
            ui.label(format!(
                "Inspect prepared {} installations and integrity metadata.",
                scope.description()
            ));
            ui.horizontal(|ui| {
                ui.label("catalog");
                ui.text_edit_singleline(&mut self.genome_catalog_path);
                if ui
                    .button("Browse...")
                    .on_hover_text("Browse filesystem and fill this path")
                    .clicked()
                    && let Some(path) = rfd::FileDialog::new()
                        .add_filter("JSON", &["json"])
                        .pick_file()
                {
                    self.genome_catalog_path = path.display().to_string();
                }
            });
            ui.horizontal(|ui| {
                ui.label("cache_dir");
                ui.text_edit_singleline(&mut self.genome_cache_dir);
                if ui
                    .button("Browse...")
                    .on_hover_text("Browse filesystem and fill this path")
                    .clicked()
                    && let Some(path) = rfd::FileDialog::new().pick_folder()
                {
                    self.genome_cache_dir = path.display().to_string();
                }
            });
            if !self.genome_catalog_error.is_empty() {
                ui.colored_label(
                    egui::Color32::from_rgb(190, 70, 70),
                    format!("Catalog error: {}", self.genome_catalog_error),
                );
            }
            match self.collect_prepared_genome_inspections() {
                Ok((inspections, errors)) => {
                    let total_size: u64 = inspections.iter().map(|r| r.total_size_bytes).sum();
                    let catalog_writable = self.selected_genome_catalog_is_writable();
                    ui.label(format!(
                        "Prepared references: {} | total size: {}",
                        inspections.len(),
                        Self::format_bytes_compact(total_size)
                    ));
                    if ui
                        .small_button("Clear Caches...")
                        .on_hover_text("Open the conservative prepared-cache cleanup specialist")
                        .clicked()
                    {
                        self.open_cache_cleanup_dialog();
                    }
                    if !catalog_writable {
                        ui.small(
                                "Catalog-entry removal is hidden here because the active catalog file is not writable.",
                            );
                    }
                    if inspections.is_empty() {
                        ui.label("No prepared references found for this catalog/cache.");
                    } else {
                        egui::ScrollArea::vertical()
                                .id_salt("prepared_genome_inspector_scroll")
                                .max_height(320.0)
                                .show(ui, |ui| {
                                    scroll_input_policy::apply_scrollarea_keyboard_navigation(
                                        ui,
                                        scroll_input_policy::DEFAULT_SCROLLAREA_KEYBOARD_STEP,
                                    );
                                    egui::Grid::new("prepared_genome_inspector_grid")
                                        .striped(true)
                                        .num_columns(9)
                                        .show(ui, |ui| {
                                            ui.strong("Genome");
                                            ui.strong("Size");
                                            ui.strong("Ready");
                                            ui.strong("Cache");
                                            ui.strong("Sources");
                                            ui.strong("SHA1 seq/ann");
                                            ui.strong("Installed");
                                            ui.strong("Path");
                                            ui.strong("");
                                            ui.end_row();
                                            for inspection in &inspections {
                                                ui.label(&inspection.genome_id);
                                                ui.label(Self::format_bytes_compact(
                                                    inspection.total_size_bytes,
                                                ));
                                                ui.label(format!(
                                                    "seq:{} ann:{} fai:{} gene:{}",
                                                    if inspection.sequence_present {
                                                        "y"
                                                    } else {
                                                        "n"
                                                    },
                                                    if inspection.annotation_present {
                                                        "y"
                                                    } else {
                                                        "n"
                                                    },
                                                    if inspection.fasta_index_ready {
                                                        "y"
                                                    } else {
                                                        "n"
                                                    },
                                                    if inspection.gene_index_ready {
                                                        "y"
                                                    } else {
                                                        "n"
                                                    }
                                                ));
                                                ui.label(
                                                    egui::RichText::new(
                                                        Self::format_prepared_cache_summary(
                                                            inspection,
                                                        ),
                                                    )
                                                    .small(),
                                                )
                                                .on_hover_text(
                                                    Self::format_prepared_cache_hover_text(
                                                        inspection,
                                                    ),
                                                );
                                                ui.label(format!(
                                                    "{}/{}",
                                                    inspection.sequence_source_type,
                                                    inspection.annotation_source_type
                                                ));
                                                ui.label(format!(
                                                    "{}/{}",
                                                    Self::format_short_sha1(
                                                        &inspection.sequence_sha1
                                                    ),
                                                    Self::format_short_sha1(
                                                        &inspection.annotation_sha1
                                                    )
                                                ));
                                                ui.label(
                                                    inspection.installed_at_unix_ms.to_string(),
                                                );
                                                let path_label = inspection.install_dir.clone();
                                                ui.label(
                                                    egui::RichText::new(path_label)
                                                        .monospace()
                                                        .small(),
                                                );
                                                ui.horizontal(|ui| {
                                                    if ui
                                                        .small_button("Retrieve")
                                                        .on_hover_text(
                                                            "Open retrieval dialog preselected for this prepared genome",
                                                        )
                                                        .clicked()
                                                    {
                                                        self.genome_id = inspection.genome_id.clone();
                                                        self.invalidate_genome_genes();
                                                        self.open_retrieve_genome_dialog_for_scope(
                                                            self.genome_dialog_scope,
                                                        );
                                                    }
                                                    if ui
                                                        .add_enabled(
                                                            self.genome_prepare_task.is_none(),
                                                            egui::Button::new("Reindex..."),
                                                        )
                                                        .on_hover_text(
                                                            "Reuse cached local files and rebuild this prepared genome's indexes using the current catalog/cache settings. A confirmation dialog will be shown because it may take some time.",
                                                        )
                                                        .clicked()
                                                    {
                                                        self.queue_prepared_genome_reinstall(
                                                            inspection.genome_id.clone(),
                                                            scope,
                                                            PreparedGenomeReinstallDialogHost::Root,
                                                        );
                                                    }
                                                    if ui
                                                        .small_button("Remove Prepared...")
                                                        .on_hover_text(
                                                            "Delete this prepared install from the cache without removing the catalog entry.",
                                                        )
                                                        .clicked()
                                                    {
                                                        self.queue_prepared_genome_removal(
                                                            inspection.genome_id.clone(),
                                                            scope,
                                                            inspection.install_dir.clone(),
                                                        );
                                                    }
                                                    if catalog_writable
                                                        && ui
                                                            .small_button("Remove Catalog Entry...")
                                                            .on_hover_text(
                                                                "Delete this row from the active writable catalog. Prepared cache files remain until 'Remove Prepared...' is run explicitly.",
                                                            )
                                                            .clicked()
                                                    {
                                                        self.queue_catalog_entry_removal(
                                                            inspection.genome_id.clone(),
                                                            scope,
                                                        );
                                                    }
                                                });
                                                ui.end_row();
                                            }
                                        });
                                });

                        ui.separator();
                        ui.label("Chromosome inspector");
                        let prepared_ids: Vec<String> =
                            inspections.iter().map(|i| i.genome_id.clone()).collect();
                        if !prepared_ids.is_empty()
                            && !prepared_ids.iter().any(|id| id == &self.genome_id)
                        {
                            self.genome_id = prepared_ids[0].clone();
                        }
                        ui.horizontal(|ui| {
                            ui.label("genome");
                            egui::ComboBox::from_id_salt("prepared_genome_chromosome_inspect")
                                .selected_text(if self.genome_id.trim().is_empty() {
                                    "<select genome>"
                                } else {
                                    self.genome_id.as_str()
                                })
                                .show_ui(ui, |ui| {
                                    for genome_id in &prepared_ids {
                                        ui.selectable_value(
                                            &mut self.genome_id,
                                            genome_id.clone(),
                                            genome_id,
                                        );
                                    }
                                });
                            if ui
                                .small_button("Use In Retrieve")
                                .on_hover_text("Open retrieval dialog preselected for this genome")
                                .clicked()
                            {
                                self.invalidate_genome_genes();
                                self.open_retrieve_genome_dialog_for_scope(
                                    self.genome_dialog_scope,
                                );
                            }
                        });
                        if !self.genome_id.trim().is_empty() {
                            match self.prepared_genome_chromosome_records(&self.genome_id) {
                                Ok(chromosomes) => {
                                    Self::render_chromosome_length_lines(ui, &chromosomes);
                                }
                                Err(e) => {
                                    ui.colored_label(
                                        egui::Color32::from_rgb(190, 70, 70),
                                        format!("Chromosome inspector error: {e}"),
                                    );
                                }
                            }
                        }
                        if matches!(scope, GenomeDialogScope::Helper) {
                            match self.selected_helper_catalog_entry_for_scope(scope) {
                                Ok(Some(entry)) => {
                                    Self::render_helper_catalog_entry_panel(
                                        ui,
                                        &entry,
                                        "Selected Helper Construct",
                                    );
                                }
                                Ok(None) => {}
                                Err(e) => {
                                    ui.colored_label(
                                        egui::Color32::from_rgb(190, 70, 70),
                                        format!("Helper interpretation error: {e}"),
                                    );
                                }
                            }
                        }
                    }
                    if !errors.is_empty() {
                        ui.separator();
                        ui.colored_label(
                            egui::Color32::from_rgb(190, 70, 70),
                            "Inspection errors:",
                        );
                        for err in errors {
                            ui.colored_label(egui::Color32::from_rgb(190, 70, 70), err);
                        }
                    }
                }
                Err(e) => {
                    ui.colored_label(
                        egui::Color32::from_rgb(190, 70, 70),
                        format!("Inspector error: {e}"),
                    );
                }
            }
        });
        if close_requested {
            open = false;
        }
        self.show_reference_genome_inspector_dialog = open;
    }

    pub(super) fn render_genome_bed_track_contents(&mut self, ui: &mut Ui) -> bool {
        let mut close_requested = false;
        let close_hover = Self::specialist_window_close_hover_text("Import Genome Tracks");
        if self.render_specialist_window_nav_with_close(ui, Some(("Close", close_hover.as_str()))) {
            close_requested = true;
        }
        let anchor_summaries = self.anchored_sequence_anchor_summaries_for_tracks();
        let anchored_seq_ids = anchor_summaries
            .iter()
            .map(|summary| summary.seq_id.clone())
            .collect::<Vec<_>>();
        let import_running = self.genome_track_import_task.is_some();
        if !anchored_seq_ids.is_empty() && !anchored_seq_ids.contains(&self.genome_track_seq_id) {
            self.genome_track_seq_id = anchored_seq_ids[0].clone();
        }

        ui.label("Map genome signal tracks onto genome-anchored sequences.");
        ui.small("Supports .bed/.bed.gz and .bw/.bigWig input files.");
        ui.small(format!(
            "Anchored sequences available: {}",
            anchored_seq_ids.len()
        ));

        if anchored_seq_ids.is_empty() {
            ui.colored_label(
                egui::Color32::from_rgb(190, 70, 70),
                "No genome-anchored sequence is available. Extract a genome region or gene first.",
            );
        } else {
            ui.horizontal(|ui| {
                ui.label("sequence");
                egui::ComboBox::from_id_salt("genome_track_seq_id_combo")
                    .selected_text(if self.genome_track_seq_id.trim().is_empty() {
                        "<select sequence>"
                    } else {
                        self.genome_track_seq_id.as_str()
                    })
                    .show_ui(ui, |ui| {
                        for seq_id in &anchored_seq_ids {
                            ui.selectable_value(
                                &mut self.genome_track_seq_id,
                                seq_id.clone(),
                                seq_id,
                            );
                        }
                    });
                let selected_resp = self.track_hover_status(
                    ui.add_enabled(
                        !anchored_seq_ids.is_empty() && !import_running,
                        egui::Button::new("Import To Selected"),
                    )
                    .on_hover_text(
                        "Import this BED/BigWig/VCF signal file onto only the currently selected anchored sequence.",
                    ),
                    "Genome Tracks > Import Selected",
                );
                if selected_resp.clicked() {
                    self.import_genome_bed_track_for_selected_sequence();
                }
            });
            if let Some(anchor) = self.describe_sequence_genome_anchor(&self.genome_track_seq_id) {
                ui.small(format!("selected anchor: {anchor}"));
            }
            ui.small(
                "Tip: Keep this window open. Change 'sequence' and click 'Import To Selected' again to reuse the same track file for another anchored sequence.",
            );
        }

        let detected_source = GenomeTrackSource::from_path(&self.genome_track_path);
        ui.horizontal(|ui| {
            ui.label("source");
            egui::ComboBox::from_id_salt("genome_track_source_selection")
                .selected_text(self.genome_track_source_selection.label())
                .show_ui(ui, |ui| {
                    ui.selectable_value(
                        &mut self.genome_track_source_selection,
                        GenomeTrackSourceSelection::Auto,
                        GenomeTrackSourceSelection::Auto.label(),
                    );
                    ui.selectable_value(
                        &mut self.genome_track_source_selection,
                        GenomeTrackSourceSelection::Bed,
                        GenomeTrackSourceSelection::Bed.label(),
                    );
                    ui.selectable_value(
                        &mut self.genome_track_source_selection,
                        GenomeTrackSourceSelection::BigWig,
                        GenomeTrackSourceSelection::BigWig.label(),
                    );
                    ui.selectable_value(
                        &mut self.genome_track_source_selection,
                        GenomeTrackSourceSelection::Vcf,
                        GenomeTrackSourceSelection::Vcf.label(),
                    );
                });
            ui.small(format!("Detected extension: {}", detected_source.label()));
        });

        ui.horizontal(|ui| {
            ui.label("track_path");
            ui.text_edit_singleline(&mut self.genome_track_path);
            let browse_track_resp = self.track_hover_status(
                ui.button("Browse...")
                    .on_hover_text("Browse filesystem and fill this path"),
                "Genome Tracks > Browse Track Path",
            );
            if browse_track_resp.clicked()
                && let Some(path) = rfd::FileDialog::new()
                    .add_filter("Signal tracks", &["bed", "gz", "bw", "bigwig", "vcf"])
                    .pick_file()
            {
                self.genome_track_path = path.display().to_string();
                if self.genome_track_name.trim().is_empty() {
                    self.genome_track_name = Self::track_name_default_from_path(&path);
                }
            }
        });
        let resolved_source = self
            .genome_track_source_selection
            .resolve(&self.genome_track_path);
        ui.small(format!(
            "Resolved source format for import: {}",
            resolved_source.label()
        ));
        ui.horizontal(|ui| {
            ui.label("track_name");
            ui.text_edit_singleline(&mut self.genome_track_name);
            ui.small("(optional)");
        });
        ui.horizontal(|ui| {
            ui.label("min_score");
            ui.add(
                egui::TextEdit::singleline(&mut self.genome_track_min_score).desired_width(90.0),
            );
            ui.label("max_score");
            ui.add(
                egui::TextEdit::singleline(&mut self.genome_track_max_score).desired_width(90.0),
            );
        });
        ui.checkbox(
            &mut self.genome_track_clear_existing,
            "Clear existing imported track features first",
        );
        if let Some(selected_anchor) = anchor_summaries
            .iter()
            .find(|summary| summary.seq_id == self.genome_track_seq_id)
        {
            let matching_anchor_count = anchor_summaries
                .iter()
                .filter(|candidate| Self::anchors_share_mapping_group(candidate, selected_anchor))
                .count();
            let projected_targets = anchored_seq_ids.len();
            let mapping_status = if matching_anchor_count == projected_targets {
                "all anchored sequences match selected genome/chromosome"
            } else {
                "mixed anchors detected; import still applies to all anchored sequences"
            };
            let track_path = self.genome_track_path.trim().to_string();
            let path_exists = !track_path.is_empty() && Path::new(&track_path).is_file();
            let projected_track_name = {
                let trimmed = self.genome_track_name.trim();
                if trimmed.is_empty() {
                    resolved_source.label().to_string()
                } else {
                    trimmed.to_string()
                }
            };

            ui.group(|ui| {
                ui.strong("Import Preflight");
                ui.small(format!(
                    "anchor: {}:{}:{}..{} (strand {})",
                    selected_anchor.genome_id,
                    selected_anchor.chromosome,
                    selected_anchor.start_1based,
                    selected_anchor.end_1based,
                    selected_anchor
                        .strand
                        .map(|s| s.to_string())
                        .unwrap_or_else(|| "?".to_string())
                ));
                ui.small(format!(
                    "matching status: {} / {} anchored sequences share this mapping group ({})",
                    matching_anchor_count, projected_targets, mapping_status
                ));
                ui.small(format!(
                    "projected tracks: {} target sequence(s) with track '{}'",
                    projected_targets, projected_track_name
                ));
                if path_exists {
                    ui.small(format!("track file: {}", track_path));
                } else {
                    ui.colored_label(
                        egui::Color32::from_rgb(190, 70, 70),
                        "track file: missing or unreadable path",
                    );
                }
                ui.horizontal(|ui| {
                    let toggle_resp = ui
                        .checkbox(
                            &mut self.genome_track_preflight_track_subscription,
                            "Track this file for auto-sync",
                        )
                        .on_hover_text(
                            "When enabled, this file is added to tracked subscriptions and auto-applied to newly anchored sequences.",
                        );
                    if toggle_resp.hovered() {
                        self.hover_status_name =
                            "Genome Tracks > Track Subscription Toggle".to_string();
                    }
                    let response = ui
                        .add_enabled(
                            !anchored_seq_ids.is_empty() && !import_running && path_exists,
                            egui::Button::new("Apply To All Anchored Now"),
                        )
                        .on_hover_text(
                            "One-click import to all anchored sequences using the current preflight settings.",
                        );
                    let response = self.track_hover_status(
                        response,
                        "Genome Tracks > Apply To All Anchored Now",
                    );
                    if response.clicked() {
                        self.import_genome_bed_track_for_all_anchored_sequences(
                            self.genome_track_preflight_track_subscription,
                        );
                    }
                });
            });
        }
        if self.genome_track_import_task.is_some() {
            let elapsed_secs = self
                .genome_track_import_task
                .as_ref()
                .map(|task| task.started.elapsed().as_secs_f32())
                .unwrap_or(0.0);
            let progress_snapshot = self.genome_track_import_progress.clone();
            ui.horizontal(|ui| {
                ui.add(egui::Spinner::new());
                if let Some(progress) = &progress_snapshot {
                    ui.label(format!(
                        "Running import: {} '{}' parsed={} imported={} skipped={} ({:.1}s)",
                        progress.source,
                        progress.path,
                        progress.parsed_records,
                        progress.imported_features,
                        progress.skipped_records,
                        elapsed_secs
                    ));
                } else {
                    ui.label(format!("Running import task... ({:.1}s)", elapsed_secs));
                }
                let cancel_import_resp = self.track_hover_status(
                    ui.button("Cancel Import").on_hover_text(
                        "Request cancellation. Imported features up to the cancellation point are kept.",
                    ),
                    "Genome Tracks > Cancel Import",
                );
                if cancel_import_resp.clicked() {
                    self.request_track_import_task_cancel("genome tracks dialog");
                }
            });
        }
        ui.horizontal(|ui| {
            let all_once_resp = self.track_hover_status(
                ui.add_enabled(
                    !anchored_seq_ids.is_empty() && !import_running,
                    egui::Button::new("Import To All Anchored (One-Time)"),
                )
                .on_hover_text(
                    "Import once onto every currently anchored sequence. No subscription is saved for future extracts.",
                ),
                "Genome Tracks > Import All Anchored One-Time",
            );
            if all_once_resp.clicked()
            {
                self.import_genome_bed_track_for_all_anchored_sequences(false);
            }
            let all_track_resp = self.track_hover_status(
                ui.add_enabled(
                    !anchored_seq_ids.is_empty() && !import_running,
                    egui::Button::new("Import To All Anchored + Track"),
                )
                .on_hover_text(
                    "Import onto all current anchored sequences and save this file as a tracked subscription for automatic future auto-sync.",
                ),
                "Genome Tracks > Import All Anchored And Track",
            );
            if all_track_resp.clicked()
            {
                self.import_genome_bed_track_for_all_anchored_sequences(true);
            }
        });

        ui.separator();
        ui.label("Tracked genome signal files for auto-sync to new anchored sequences");
        ui.horizontal(|ui| {
            ui.label("Filter")
                .on_hover_text(Self::tracked_file_filter_help_text());
            ui.add(
                egui::TextEdit::singleline(&mut self.genome_track_subscription_filter)
                    .hint_text("path, track, source")
                    .desired_width(280.0),
            )
            .on_hover_text(Self::tracked_file_filter_help_text());
            if ui
                .add_enabled(
                    !self.genome_track_subscription_filter.trim().is_empty(),
                    egui::Button::new("Clear"),
                )
                .clicked()
            {
                self.genome_track_subscription_filter.clear();
            }
        });
        ui.horizontal_wrapped(|ui| {
            ui.small("Presets:");
            for preset in [
                "source:bed",
                "source:bigwig",
                "source:vcf",
                "path:.bed",
                "track:chip",
            ] {
                if ui.small_button(preset).clicked() {
                    Self::append_filter_term(&mut self.genome_track_subscription_filter, preset);
                }
            }
        });

        let mut apply_now_index: Option<usize> = None;
        let mut remove_index: Option<usize> = None;
        if self.genome_bed_track_subscriptions.is_empty() {
            ui.small("No tracked files yet.");
        } else {
            let subscription_rows = self.genome_bed_track_subscriptions.clone();
            let filter_text = self.genome_track_subscription_filter.trim().to_string();
            let filtered_rows = subscription_rows
                .iter()
                .enumerate()
                .filter(|(_, subscription)| {
                    Self::subscription_matches_filter(subscription, &filter_text)
                })
                .collect::<Vec<_>>();
            ui.small(format!(
                "Showing {} of {} tracked files",
                filtered_rows.len(),
                subscription_rows.len()
            ));
            if filtered_rows.is_empty() {
                ui.small("No tracked files match the current filter.");
            }
            egui::Grid::new("genome_bed_track_subscriptions_grid")
                .striped(true)
                .show(ui, |ui| {
                    ui.strong("Type");
                    ui.strong("Path");
                    ui.strong("Track");
                    ui.strong("Score Range");
                    ui.strong("Clear Existing");
                    ui.strong("Actions");
                    ui.end_row();
                    for (index, subscription) in filtered_rows {
                        ui.label(subscription.source.label());
                        ui.monospace(subscription.path.as_str());
                        ui.label(
                            subscription
                                .track_name
                                .clone()
                                .unwrap_or_else(|| "-".to_string()),
                        );
                        ui.label(format!(
                            "{}..{}",
                            subscription
                                .min_score
                                .map(|v| format!("{v:.3}"))
                                .unwrap_or_else(|| "-".to_string()),
                            subscription
                                .max_score
                                .map(|v| format!("{v:.3}"))
                                .unwrap_or_else(|| "-".to_string())
                        ));
                        ui.label(if subscription.clear_existing {
                            "yes"
                        } else {
                            "no"
                        });
                        ui.horizontal(|ui| {
                            let apply_resp = self.track_hover_status(
                                ui.add_enabled(!import_running, egui::Button::new("Apply now"))
                                    .on_hover_text(
                                        "Re-apply this tracked file to all currently anchored sequences now.",
                                    ),
                                "Genome Tracks > Apply Tracked File",
                            );
                            if apply_resp.clicked()
                            {
                                apply_now_index = Some(index);
                            }
                            let remove_resp = self.track_hover_status(
                                ui.add_enabled(!import_running, egui::Button::new("Remove"))
                                    .on_hover_text(
                                        "Remove this tracked subscription (does not remove already imported features).",
                                    ),
                                "Genome Tracks > Remove Tracked File",
                            );
                            if remove_resp.clicked()
                            {
                                remove_index = Some(index);
                            }
                        });
                        ui.end_row();
                    }
                });
        }

        if let Some(index) = apply_now_index {
            self.apply_tracked_bed_subscription_to_all_anchored(index);
        }
        if let Some(index) = remove_index {
            let result = self
                .engine
                .write()
                .unwrap()
                .remove_genome_track_subscription(index);
            match result {
                Ok(removed) => {
                    self.load_bed_track_subscriptions_from_state();
                    self.genome_track_status = format!(
                        "Removed tracked {} '{}'",
                        removed.source.label(),
                        removed.path
                    );
                }
                Err(e) => {
                    self.genome_track_status =
                        format!("Could not remove tracked file: {}", e.message);
                }
            }
        }
        if !self.genome_bed_track_subscriptions.is_empty() {
            let clear_resp = self.track_hover_status(
                ui.add_enabled(!import_running, egui::Button::new("Clear Tracked Files"))
                    .on_hover_text(
                        "Remove all tracked subscriptions (does not remove already imported features).",
                    ),
                "Genome Tracks > Clear Tracked Files",
            );
            if clear_resp.clicked() {
                self.engine
                    .write()
                    .unwrap()
                    .clear_genome_track_subscriptions();
                self.load_bed_track_subscriptions_from_state();
                self.genome_track_status = "Cleared all tracked files".to_string();
            }
        }

        if !self.genome_track_autosync_status.is_empty() {
            ui.separator();
            ui.monospace(&self.genome_track_autosync_status);
        }
        if !self.genome_track_status.is_empty() {
            ui.separator();
            ui.monospace(&self.genome_track_status);
        }
        close_requested
    }

    pub(super) fn render_genome_bed_track_dialog(&mut self, ctx: &egui::Context) {
        if !self.show_genome_bed_track_dialog {
            return;
        }
        let viewport_id = Self::bed_track_viewport_id();
        let spec = self.hosted_window_spec_for_viewport(
            "Import Genome Tracks",
            egui::Id::new(("hosted_bed_track_window", viewport_id)),
            viewport_id,
            Vec2::new(980.0, 620.0),
            Vec2::new(620.0, 320.0),
        );
        if ctx.embed_viewports() {
            let mut open = self.show_genome_bed_track_dialog;
            let mut close_requested = false;
            crate::egui_compat::show_hosted_window(ctx, &spec, &mut open, |ui| {
                close_requested = self.render_genome_bed_track_contents(ui)
            });
            self.clear_viewport_foreground_request_after_render(viewport_id);
            if close_requested {
                open = false;
            }
            self.show_genome_bed_track_dialog = open;
            return;
        }
        let builder = crate::egui_compat::viewport_builder_for_hosted_window(&spec);
        ctx.show_viewport_immediate(viewport_id, builder, |ctx, class| {
            self.note_viewport_focus_if_active(ctx, viewport_id);
            if class == egui::ViewportClass::EmbeddedWindow {
                let mut open = self.show_genome_bed_track_dialog;
                let mut close_requested = false;
                crate::egui_compat::show_hosted_window(ctx, &spec, &mut open, |ui| {
                    close_requested = self.render_genome_bed_track_contents(ui)
                });
                if close_requested {
                    open = false;
                }
                self.show_genome_bed_track_dialog = open;
                return;
            }

            let mut close_requested = false;
            crate::egui_compat::show_central_panel(ctx, egui::CentralPanel::default(), |ui| {
                close_requested = self.render_genome_bed_track_contents(ui);
            });

            if close_requested || Self::viewport_close_requested_or_shortcut(ctx) {
                self.show_genome_bed_track_dialog = false;
            }
        });
    }

    pub(super) fn render_cache_cleanup_dialog(&mut self, ctx: &egui::Context) {
        if !self.show_cache_cleanup_dialog {
            return;
        }
        let mut open = self.show_cache_cleanup_dialog;
        let mut close_requested = false;
        let spec = crate::egui_compat::ModalWindowSpec::new(
            "Clear Caches",
            egui::Id::new("clear_caches_modal"),
        )
        .default_size(Vec2::new(980.0, 720.0))
        .min_size(Vec2::new(700.0, 420.0))
        .resizable(true);
        crate::egui_compat::show_modal_window(ctx, &spec, &mut open, |ui| {
            let close_hover = Self::specialist_window_close_hover_text("Genome > Clear Caches");
            if self
                .render_specialist_window_nav_with_close(ui, Some(("Close", close_hover.as_str())))
            {
                close_requested = true;
                return;
            }
            ui.label(
                    "Inspect and conservatively clean prepared genome/helper caches. Project state files, catalog JSON, MCP/runtime files, and build artifacts are out of scope.",
                );
            ui.separator();
            ui.horizontal(|ui| {
                ui.label("scope");
                for scope in [
                    CacheCleanupScope::References,
                    CacheCleanupScope::Helpers,
                    CacheCleanupScope::Both,
                ] {
                    let clicked = ui
                        .radio_value(&mut self.cache_cleanup_scope, scope, scope.label())
                        .clicked();
                    if clicked {
                        self.cache_cleanup_confirm_pending = false;
                        self.cache_cleanup_rebuild_candidates.clear();
                    }
                }
            });
            ui.horizontal(|ui| {
                ui.label("reference cache");
                let edited = ui
                    .text_edit_singleline(&mut self.cache_cleanup_reference_cache_dir)
                    .changed();
                if ui
                    .button("Browse...")
                    .on_hover_text("Browse filesystem and fill this path")
                    .clicked()
                    && let Some(path) = rfd::FileDialog::new().pick_folder()
                {
                    self.cache_cleanup_reference_cache_dir = path.display().to_string();
                }
                if edited {
                    self.cache_cleanup_confirm_pending = false;
                    self.cache_cleanup_rebuild_candidates.clear();
                }
            });
            ui.horizontal(|ui| {
                ui.label("helper cache");
                let edited = ui
                    .text_edit_singleline(&mut self.cache_cleanup_helper_cache_dir)
                    .changed();
                if ui
                    .button("Browse...")
                    .on_hover_text("Browse filesystem and fill this path")
                    .clicked()
                    && let Some(path) = rfd::FileDialog::new().pick_folder()
                {
                    self.cache_cleanup_helper_cache_dir = path.display().to_string();
                }
                if edited {
                    self.cache_cleanup_confirm_pending = false;
                    self.cache_cleanup_rebuild_candidates.clear();
                }
            });
            ui.horizontal(|ui| {
                    if ui
                        .button("Inspect")
                        .on_hover_text("Inspect prepared installs and orphaned remnants in the selected cache roots")
                        .clicked()
                    {
                        self.refresh_cache_cleanup_inspection();
                    }
                    ui.small(format!(
                        "Selected roots: {}",
                        self.cache_cleanup_selected_roots().join(", ")
                    ));
                });
            ui.separator();
            ui.label("Cleanup mode");
            let modes = [
                (
                    PreparedCacheCleanupMode::DerivedIndexesOnly,
                    "Remove rebuildable indexes",
                ),
                (
                    PreparedCacheCleanupMode::BlastDbOnly,
                    "Remove BLAST databases only",
                ),
                (
                    PreparedCacheCleanupMode::SelectedPreparedInstalls,
                    "Remove selected prepared installs",
                ),
                (
                    PreparedCacheCleanupMode::AllPreparedInCache,
                    "Remove all prepared installs in cache",
                ),
            ];
            for (mode, label) in modes {
                if ui
                    .radio_value(&mut self.cache_cleanup_mode, mode, label)
                    .clicked()
                {
                    self.cache_cleanup_confirm_pending = false;
                    self.cache_cleanup_rebuild_candidates.clear();
                }
            }
            let include_orphans_enabled = self.cache_cleanup_mode.allows_orphaned_remnants();
            let include_changed = ui
                .add_enabled(
                    include_orphans_enabled,
                    egui::Checkbox::new(
                        &mut self.cache_cleanup_include_orphans,
                        "Include orphaned remnants",
                    ),
                )
                .changed();
            if include_changed {
                self.cache_cleanup_confirm_pending = false;
                self.cache_cleanup_rebuild_candidates.clear();
            }
            if !include_orphans_enabled {
                self.cache_cleanup_include_orphans = false;
            }
            ui.separator();
            ui.label("Discovered items");
            if let Some(report) = &self.cache_cleanup_inspection {
                if report.entries.is_empty() {
                    ui.small("Nothing to clean in selected cache roots.");
                } else {
                    egui::ScrollArea::vertical()
                            .id_salt("cache_cleanup_discovered_items_scroll")
                            .max_height(260.0)
                            .show(ui, |ui| {
                                egui::Grid::new("cache_cleanup_grid")
                                    .striped(true)
                                    .num_columns(7)
                                    .show(ui, |ui| {
                                        ui.strong("");
                                        ui.strong("Id");
                                        ui.strong("Kind");
                                        ui.strong("Cache root");
                                        ui.strong("Bytes");
                                        ui.strong("Artifacts");
                                        ui.strong("Path");
                                        ui.end_row();
                                        for entry in &report.entries {
                                            let can_select = self.cache_cleanup_mode
                                                != PreparedCacheCleanupMode::AllPreparedInCache;
                                            let mut selected = self
                                                .cache_cleanup_selected_paths
                                                .contains(&entry.path);
                                            if ui
                                                .add_enabled(
                                                    can_select,
                                                    egui::Checkbox::new(&mut selected, ""),
                                                )
                                                .changed()
                                            {
                                                self.cache_cleanup_confirm_pending = false;
                                                self.cache_cleanup_rebuild_candidates.clear();
                                                if selected {
                                                    self.cache_cleanup_selected_paths
                                                        .insert(entry.path.clone());
                                                } else {
                                                    self.cache_cleanup_selected_paths
                                                        .remove(&entry.path);
                                                }
                                            }
                                            ui.label(&entry.entry_id);
                                            ui.label(match entry.classification {
                                                crate::genomes::PreparedCacheEntryKind::PreparedInstall => {
                                                    "prepared install"
                                                }
                                                crate::genomes::PreparedCacheEntryKind::OrphanedRemnant => {
                                                    "orphaned remnant"
                                                }
                                            });
                                            ui.label(
                                                Path::new(&entry.cache_root)
                                                    .file_name()
                                                    .and_then(|value| value.to_str())
                                                    .unwrap_or(entry.cache_root.as_str()),
                                            );
                                            ui.label(Self::format_bytes_compact(
                                                entry.total_size_bytes,
                                            ));
                                            ui.label(Self::format_cache_cleanup_groups(entry));
                                            ui.label(
                                                egui::RichText::new(entry.path.clone())
                                                    .small()
                                                    .monospace(),
                                            );
                                            ui.end_row();
                                        }
                                    });
                            });
                }
                let (targets, preview_bytes, root_count, explanation) =
                    self.cache_cleanup_preview_targets();
                let target_count = targets.len();
                let has_targets = target_count > 0;
                let roots_text = self.cache_cleanup_selected_roots().join(", ");
                ui.separator();
                ui.label("Preview");
                ui.small(format!(
                    "Targets: {} item(s) across {} selected cache root(s) | reclaimable: {}",
                    target_count,
                    root_count,
                    Self::format_bytes_compact(preview_bytes)
                ));
                ui.small(explanation);
                ui.small(format!("Roots: {roots_text}"));
                ui.horizontal(|ui| {
                    if ui
                        .add_enabled(has_targets, egui::Button::new("Review Cleanup"))
                        .on_hover_text("Arm the confirm step for the current cleanup preview")
                        .clicked()
                    {
                        self.cache_cleanup_confirm_pending = true;
                    }
                    if ui
                        .add_enabled(
                            self.cache_cleanup_confirm_pending && has_targets,
                            egui::Button::new("Confirm Cleanup"),
                        )
                        .on_hover_text("Delete the currently previewed cache artifacts")
                        .clicked()
                    {
                        self.apply_cache_cleanup();
                    }
                    if ui.button("Cancel").clicked() {
                        close_requested = true;
                    }
                });
                if !self.cache_cleanup_rebuild_candidates.is_empty() {
                    ui.separator();
                    ui.label("Rebuild from cached files");
                    ui.small(
                            "These prepared installs kept their cached FASTA/annotation sources. Reindex them through the existing cached-files rebuild flow if you want the indexes back now.",
                        );
                    for index in 0..self.cache_cleanup_rebuild_candidates.len() {
                        let request = self.cache_cleanup_rebuild_candidates[index].clone();
                        ui.horizontal(|ui| {
                                if ui
                                    .add_enabled(
                                        self.genome_prepare_task.is_none(),
                                        egui::Button::new(format!(
                                            "Rebuild '{}'",
                                            request.genome_id
                                        )),
                                    )
                                    .on_hover_text(
                                        "Open the standard cached-files reindex confirmation for this prepared install.",
                                    )
                                    .clicked()
                                {
                                    self.queue_cache_cleanup_rebuild_candidate(index);
                                }
                                ui.small(format!(
                                    "{} | cache {}",
                                    request.scope.description(),
                                    request.cache_dir
                                ));
                            });
                    }
                    if self.genome_prepare_task.is_some() {
                        ui.small(
                                "Finish or cancel the running genome prepare task before launching another rebuild.",
                            );
                    }
                }
            } else {
                ui.small("Inspect selected cache roots to preview removable items.");
            }
            if !self.cache_cleanup_status.trim().is_empty() {
                ui.separator();
                ui.monospace(self.cache_cleanup_status.clone());
            }
        });
        self.show_cache_cleanup_dialog = open && !close_requested;
    }
}
