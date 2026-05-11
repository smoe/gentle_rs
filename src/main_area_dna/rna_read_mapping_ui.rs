//! RNA-read mapping workspace UI support for `MainAreaDna`.
//!
//! This submodule owns the dedicated RNA-read mapping workspace controls,
//! evidence panels, task/status helpers, and report preview rendering. The
//! shared RNA-read state/cache records remain in `rna_read_support`.

use super::*;

impl MainAreaDna {
    pub(super) fn render_rna_read_mapping_workspace_controls(
        &mut self,
        ui: &mut egui::Ui,
        view: &SplicingExpertView,
    ) {
        let previous_op_status = self.op_status.clone();
        let mut persist_ui_state = false;
        self.render_rna_read_mapping_status(ui);
        let workspace_report_id = self.current_rna_read_mapping_workspace_report_id();
        let workspace_raw_report = workspace_report_id
            .as_deref()
            .and_then(|report_id| self.get_saved_rna_read_report_by_id(report_id));
        let workspace_report_view_mismatch = workspace_raw_report
            .as_ref()
            .filter(|report| !Self::rna_read_report_matches_splicing_view(report.as_ref(), view))
            .map(|report| Self::rna_read_report_view_mismatch_message(report.as_ref(), view));
        let workspace_saved_report = workspace_raw_report
            .as_ref()
            .filter(|report| Self::rna_read_report_matches_splicing_view(report.as_ref(), view))
            .cloned();
        let report_id_usable_for_view = workspace_report_view_mismatch.is_none();
        if let Some(message) = workspace_report_view_mismatch.as_ref() {
            ui.small(egui::RichText::new(message).color(egui::Color32::from_rgb(180, 83, 9)));
        }
        ui.add_space(4.0);
        ui.horizontal_wrapped(|ui| {
            ui.label(
                egui::RichText::new(Self::rna_read_mapping_parameter_section_title())
                    .strong()
                    .color(egui::Color32::from_rgb(51, 65, 85)),
            );
            ui.small(
                egui::RichText::new(format!(
                    "Current profile: {}",
                    self.rna_reads_ui.profile.as_str()
                ))
                .color(egui::Color32::from_rgb(100, 116, 139)),
            );
        });
        ui.horizontal_wrapped(|ui| {
            ui.label(
                egui::RichText::new("Mapping guide [?]")
                    .size(9.0)
                    .color(egui::Color32::from_rgb(71, 85, 105)),
            )
            .on_hover_text(Self::splicing_nanopore_cdna_panel_help_text());
            ui.label(
                egui::RichText::new(
                    "Phase 1 keeps a retained top-hit report; phase 2 aligns that saved report and refreshes exon/junction support.",
                )
                .size(9.0)
                .color(egui::Color32::from_rgb(100, 116, 139)),
            );
        });
        ui.label(
            egui::RichText::new(
                "Phase-1 path: FASTA input (.fa/.fasta, optional .gz); .sra requires external conversion.",
            )
            .size(9.0)
            .color(egui::Color32::from_rgb(100, 116, 139)),
        )
        .on_hover_text(Self::splicing_nanopore_cdna_panel_help_text());
        let controls_enabled = self.rna_read_task.is_none();
        ui.add_enabled_ui(controls_enabled, |ui| {
            let mut refresh_auto_report_id = self.rna_reads_ui.report_id_auto_sync
                && self.rna_reads_ui.report_id.trim().is_empty()
                && !self.rna_reads_ui.input_path.trim().is_empty();
            ui.horizontal(|ui| {
                ui.label("Input FASTA").on_hover_text(
                    "Path to phase-1 input reads in FASTA format (.fa/.fasta, optional .gz). Reads are streamed sequentially from this file; .sra must be converted externally first.",
                );
                if ui
                    .text_edit_singleline(&mut self.rna_reads_ui.input_path)
                    .on_hover_text(
                        "InterpretRnaReads streams reads from this file during phase 1. The retained report stores scored rows, not the original reads file itself, so keep the path if you plan to rerun with different thresholds.",
                    )
                    .changed()
                {
                    persist_ui_state = true;
                    refresh_auto_report_id = true;
                }
                if ui
                    .button("Browse...")
                    .on_hover_text(
                        "Open a file chooser for FASTA/FASTA.gz input. In Auto mode, Report ID is refreshed from the current locus + input + scope/origin settings.",
                    )
                    .clicked()
                {
                    if let Some(path) = rfd::FileDialog::new()
                        .add_filter("FASTA", &["fa", "fasta", "gz"])
                        .pick_file()
                    {
                        self.rna_reads_ui.input_path = path.display().to_string();
                        persist_ui_state = true;
                        refresh_auto_report_id = true;
                    }
                }
            });
            ui.horizontal(|ui| {
                ui.label("Report ID").on_hover_text(
                    "Identifier used to store and retrieve the retained top-hit report produced by phase 1. The same ID is reused by phase-2 alignment, inspection, and TSV/SVG export actions.",
                );
                let report_id_changed = ui
                    .text_edit_singleline(&mut self.rna_reads_ui.report_id)
                    .on_hover_text(
                        "Auto mode keeps this ID synchronized with the current locus, input file, scope, and origin settings. Typing a custom value switches to manual mode until you re-enable Auto or press Refresh ID.",
                    )
                    .changed();
                if report_id_changed {
                    if self.rna_reads_ui.report_id.trim().is_empty() {
                        self.rna_reads_ui.report_id_auto_sync = true;
                        refresh_auto_report_id = true;
                    } else {
                        self.rna_reads_ui.report_id_auto_sync = false;
                    }
                    persist_ui_state = true;
                }
                let auto_toggled = ui
                    .checkbox(&mut self.rna_reads_ui.report_id_auto_sync, "Auto")
                    .on_hover_text(
                        "Keep Report ID synchronized with the current locus, input, and run settings instead of treating it as a manual override.",
                    )
                    .changed();
                if auto_toggled {
                    persist_ui_state = true;
                    if self.rna_reads_ui.report_id_auto_sync {
                        refresh_auto_report_id = true;
                    }
                }
                if ui
                    .button("Refresh ID")
                    .on_hover_text(
                        "Regenerate the suggested Report ID from the current locus, input, scope, and origin settings without changing any other mapping controls.",
                    )
                    .clicked()
                {
                    self.rna_reads_ui.report_id =
                        Self::default_rna_read_report_id(view, &self.rna_reads_ui);
                    persist_ui_state = true;
                }
            });
            ui.horizontal(|ui| {
                ui.label("Region / gene scope").on_hover_text(
                    "Controls which genes/transcript templates contribute exon-body and junction seed hashes. Any-strand scopes admit target-gene-strand and antisense/opposite-strand templates; target-gene-strand scopes keep the run on the selected target gene/group's annotated strand.",
                );
                let mut scope_changed = false;
                egui::ComboBox::from_id_salt(format!(
                    "rna_read_scope_{}_{}",
                    view.seq_id, view.target_feature_id
                ))
                .selected_text(Self::rna_read_scope_selection_label(
                    view,
                    self.rna_reads_ui.scope,
                ))
                .show_ui(ui, |ui| {
                    scope_changed |= ui
                        .selectable_value(
                            &mut self.rna_reads_ui.scope,
                            SplicingScopePreset::AllOverlappingAnyStrand,
                            "All overlapping incl. antisense (any strand)",
                        )
                        .on_hover_text(
                            "Broadest mode: include all overlapping transcripts on any strand, including antisense/opposite-strand genes relative to the selected target gene/group.",
                        )
                        .changed();
                    scope_changed |= ui
                        .selectable_value(
                            &mut self.rna_reads_ui.scope,
                            SplicingScopePreset::TargetGroupAnyStrand,
                            "Target gene/group (any strand)",
                        )
                        .on_hover_text(
                            "Restrict to the current target gene/group, but allow any annotated strand inside that group.",
                        )
                        .changed();
                    scope_changed |= ui
                        .selectable_value(
                            &mut self.rna_reads_ui.scope,
                            SplicingScopePreset::AllOverlappingTargetStrand,
                            "All overlapping (target-gene strand)",
                        )
                        .on_hover_text(
                            "Include all overlapping transcripts only on the selected target gene/group's annotated strand; antisense/opposite-strand genes are excluded.",
                        )
                        .changed();
                    scope_changed |= ui
                        .selectable_value(
                            &mut self.rna_reads_ui.scope,
                            SplicingScopePreset::TargetGroupTargetStrand,
                            "Target gene/group (target-gene strand)",
                        )
                        .on_hover_text(
                            "Most specific mode: the current target gene/group on the selected target gene/group's annotated strand only.",
                        )
                        .changed();
                });
                persist_ui_state |= scope_changed;
                refresh_auto_report_id |= scope_changed;
                persist_ui_state |= ui
                    .checkbox(&mut self.rna_reads_ui.show_advanced, "Show advanced")
                    .on_hover_text(
                        "Show deterministic seed-gate, origin-expansion, checkpoint/resume, and phase-2 alignment controls shared with InterpretRnaReads/AlignRnaReadReport across GUI, CLI, JS, and Lua.",
                    )
                    .changed();
            });
            if self.rna_reads_ui.show_advanced {
                ui.horizontal_wrapped(|ui| {
                    ui.label("Gene expansion mode").on_hover_text(
                        "Controls how transcript templates are gathered before hashing. This is still one run, not multiple invocations: `single_gene` uses only the current target gene/group; `multi_gene_sparse` adds transcript templates from the explicit Target genes list. This remains local annotation-driven, not a genome-wide search.",
                    );
                    let mut origin_changed = false;
                    egui::ComboBox::from_id_salt(format!(
                        "rna_read_origin_mode_{}_{}",
                        view.seq_id, view.target_feature_id
                    ))
                    .selected_text(Self::rna_read_origin_mode_selection_label(
                        view,
                        self.rna_reads_ui.origin_mode,
                    ))
                    .show_ui(ui, |ui| {
                        origin_changed |= ui
                            .selectable_value(
                                &mut self.rna_reads_ui.origin_mode,
                                RnaReadOriginMode::SingleGene,
                                "single_gene (baseline)",
                        )
                        .on_hover_text(
                            "Current deterministic baseline: use current seed feature/scope only.",
                        )
                        .changed();
                        origin_changed |= ui
                            .selectable_value(
                                &mut self.rna_reads_ui.origin_mode,
                                RnaReadOriginMode::MultiGeneSparse,
                                "multi_gene_sparse",
                        )
                        .on_hover_text(
                            "Expand transcript templates from target_gene_ids (local annotation, deterministic). ROI seed-capture remains a planned follow-up.",
                        )
                        .changed();
                    });
                    persist_ui_state |= origin_changed;
                    refresh_auto_report_id |= origin_changed;
                    ui.label("Target genes").on_hover_text(
                        "Optional gene IDs for multi-gene sparse mode. Comma/space/semicolon separated.",
                    );
                    let target_genes_changed = ui
                        .add(
                            egui::TextEdit::singleline(&mut self.rna_reads_ui.target_gene_ids)
                                .desired_width(280.0)
                                .hint_text(format!("e.g. {}, TP53", view.group_label)),
                        )
                        .on_hover_text(
                            "Example: TP73, TP53. Applied when origin mode is multi_gene_sparse and persisted in the report payload.",
                        )
                        .changed();
                    persist_ui_state |= target_genes_changed;
                    refresh_auto_report_id |= target_genes_changed;
                    persist_ui_state |= ui
                        .checkbox(
                            &mut self.rna_reads_ui.roi_seed_capture_enabled,
                            "ROI seed capture (planned)",
                        )
                        .on_hover_text(
                            "Planned annotation-independent ROI seed-capture layer. Persisted now; engine emits a deterministic warning until implemented.",
                        )
                        .changed();
                });
                ui.small(
                    egui::RichText::new(
                        "Sparse-origin settings are persisted in report metadata. multi_gene_sparse actively expands local transcript templates; ROI seed capture remains a planned follow-up.",
                    )
                    .color(egui::Color32::from_rgb(100, 116, 139)),
                );
                ui.add_space(4.0);
                ui.horizontal_wrapped(|ui| {
                    ui.label("Report mode").on_hover_text(
                        "Controls how much of the retained top-hit set is persisted under Report ID. This affects later inspection/export size, not the live seed scoring decisions themselves.",
                    );
                    egui::ComboBox::from_id_salt(format!(
                        "rna_read_report_mode_{}_{}",
                        view.seq_id, view.target_feature_id
                    ))
                    .selected_text(self.rna_reads_ui.report_mode.as_str())
                    .show_ui(ui, |ui| {
                        persist_ui_state |= ui
                            .selectable_value(
                                &mut self.rna_reads_ui.report_mode,
                                RnaReadReportMode::Full,
                                "full",
                            )
                            .on_hover_text(
                                "Persist retained hits exactly as ranked by the seed-stage retention policy.",
                            )
                            .changed();
                        persist_ui_state |= ui
                            .selectable_value(
                                &mut self.rna_reads_ui.report_mode,
                                RnaReadReportMode::SeedPassedOnly,
                                "seed_passed_only",
                            )
                            .on_hover_text(
                                "Persist a smaller retained report that still keeps rows useful for later review: composite seed-pass hits plus retained rows at or above raw min_hit.",
                            )
                            .changed();
                    });
                    ui.label("Checkpoint path").on_hover_text(
                        "Optional JSON checkpoint file for deterministic pause/resume.",
                    );
                    persist_ui_state |= ui
                        .add(
                            egui::TextEdit::singleline(&mut self.rna_reads_ui.checkpoint_path)
                                .desired_width(260.0),
                        )
                        .on_hover_text("If set, periodic checkpoint snapshots are written here.")
                        .changed();
                    ui.label("every reads").on_hover_text(
                        "Checkpoint write cadence in processed reads (must be > 0).",
                    );
                    persist_ui_state |= ui
                        .add(
                            egui::TextEdit::singleline(
                                &mut self.rna_reads_ui.checkpoint_every_reads,
                            )
                            .desired_width(68.0),
                        )
                        .on_hover_text(
                            "Example: 10000 writes a checkpoint every 10k processed reads.",
                        )
                        .changed();
                    persist_ui_state |= ui
                        .checkbox(&mut self.rna_reads_ui.resume_from_checkpoint, "Resume")
                        .on_hover_text(
                            "Resume from checkpoint_path. Requires checkpoint_path to be set.",
                        )
                        .changed();
                });
                ui.small(
                    egui::RichText::new(
                        "Checkpoint+resume settings map directly to InterpretRnaReads runtime options and are shared with CLI/JS/Lua.",
                    )
                    .color(egui::Color32::from_rgb(100, 116, 139)),
                );
            }
            let (gene_scope_summary, gene_scope_color) = Self::rna_read_gene_scope_summary(
                view,
                self.rna_reads_ui.scope,
                self.rna_reads_ui.origin_mode,
                &self.rna_reads_ui.target_gene_ids,
            );
            ui.small(egui::RichText::new(gene_scope_summary).color(gene_scope_color));
            ui.horizontal_wrapped(|ui| {
                persist_ui_state |= ui
                    .checkbox(
                        &mut self.rna_reads_ui.cdna_poly_t_flip_enabled,
                        "Input is cDNA (normalize T-rich 5' head)",
                    )
                    .on_hover_text(
                        "If enabled, reads with a strong T-rich 5' head are reverse-complement normalized before scoring so cDNA reads are compared in transcript orientation. Disable this for direct RNA or when input orientation is already known to be correct.",
                    )
                    .changed();
                if ui
                    .button("Apply demo specificity preset")
                    .on_hover_text(
                        "Apply stricter TP73-focused defaults: target-group/target-strand scope plus tighter chain, gap, and transition thresholds for focused pilot filtering.",
                    )
                    .clicked()
                {
                    self.apply_rna_reads_demo_specificity_preset();
                    persist_ui_state = true;
                    refresh_auto_report_id = true;
                }
                ui.small(
                    egui::RichText::new(
                        "Unchecked = direct RNA mode (no automatic reverse-complement normalization).",
                    )
                    .color(egui::Color32::from_rgb(100, 116, 139)),
                );
            });
            let (scope_title, scope_details, scope_strand_note) =
                Self::splicing_scope_description(self.rna_reads_ui.scope);
            ui.small(
                egui::RichText::new(format!("Scope detail: {scope_title}"))
                    .color(egui::Color32::from_rgb(71, 85, 105)),
            );
            ui.small(
                egui::RichText::new(scope_details)
                    .color(egui::Color32::from_rgb(100, 116, 139)),
            );
            ui.small(
                egui::RichText::new(scope_strand_note)
                    .color(egui::Color32::from_rgb(100, 116, 139)),
            );
            self.render_rna_read_scope_eligibility_sketch(
                ui,
                view,
                self.rna_reads_ui.scope,
                self.rna_reads_ui.origin_mode,
            );
            ui.small(
                egui::RichText::new(
                    "Seed-index note: annotated exon bodies and exon-exon junction transitions are indexed for transcripts admitted by the selected scope.",
                )
                .color(egui::Color32::from_rgb(100, 116, 139)),
            );
            ui.label(
                egui::RichText::new(format!(
                    "Profile: {} | Input: fasta | Origin mode: {} | Target genes: {} | Read mode: {} | Seed gate: raw>=min hit AND weighted>=min weighted AND unique>=min(min unique, tested kmers) AND chain>=min chain AND median transcript gap<=max gap AND confirmed transitions>=min transitions AND transition fraction>=min transition frac | use 'Run alignment phase' for retained-hit mapping",
                    self.rna_reads_ui.profile.as_str(),
                    self.rna_reads_ui.origin_mode.as_str(),
                    Self::parse_rna_target_gene_ids(&self.rna_reads_ui.target_gene_ids).len(),
                    if self.rna_reads_ui.cdna_poly_t_flip_enabled {
                        "cDNA (T-rich 5' head normalization enabled; minor interruptions tolerated)"
                    } else {
                        "direct RNA (no poly-T flip)"
                    }
                ))
                .size(9.0)
                .color(egui::Color32::from_rgb(100, 116, 139)),
            );
            ui.horizontal_wrapped(|ui| {
                ui.small(
                    egui::RichText::new(self.rna_read_hash_parameter_summary())
                        .color(egui::Color32::from_rgb(71, 85, 105)),
                );
                ui.separator();
                ui.small(
                    egui::RichText::new(self.rna_read_dotplot_parameter_summary())
                        .color(egui::Color32::from_rgb(71, 85, 105)),
                );
            });
            ui.horizontal_wrapped(|ui| {
                if ui
                    .small_button("Dense 9-mer similarity preset")
                    .on_hover_text(
                        "Set RNA-read hashing to exact dense 9-mers (k=9/stride=1) and RNA-read dotplots to exact dense 9-mers too (word=9/step=1/max mismatches=0).",
                    )
                    .clicked()
                {
                    self.apply_rna_read_dense_similarity_preset();
                    persist_ui_state = true;
                }
                if ui
                    .small_button("Reset dotplot defaults")
                    .on_hover_text(
                        "Restore the shared dotplot defaults used for new dotplot workspaces.",
                    )
                    .clicked()
                {
                    self.reset_rna_read_dotplot_parameters_to_defaults();
                    persist_ui_state = true;
                }
                if !self.rna_reads_ui.show_advanced
                    && ui
                        .small_button("Show tuning knobs")
                        .on_hover_text(
                            "Reveal the editable hashing and RNA-read dotplot parameters below.",
                        )
                        .clicked()
                {
                    self.rna_reads_ui.show_advanced = true;
                    persist_ui_state = true;
                }
                if ui
                    .small_button("Open Dotplot workspace")
                    .on_hover_text(
                        "Open the full dotplot workspace. RNA-read dotplot exports reuse the same word/step/mismatch/tile settings shown here.",
                    )
                    .clicked()
                {
                    self.open_dotplot_window();
                }
            });
            if self.rna_reads_ui.show_advanced {
                ui.horizontal_wrapped(|ui| {
                    ui.label("k-mer")
                        .on_hover_text("Seed length used for phase-1 hash matching.");
                    persist_ui_state |= ui
                        .add(
                            egui::TextEdit::singleline(&mut self.rna_reads_ui.kmer_len)
                                .desired_width(46.0),
                        )
                        .on_hover_text("Seed hash length in bases.")
                        .changed();
                    ui.label("hash stride").on_hover_text(
                        "Seed-start spacing in bases along each read. 1 hashes every possible start; higher values make the initial hash screen sparser and faster.",
                    );
                    persist_ui_state |= ui
                        .add(
                            egui::TextEdit::singleline(&mut self.rna_reads_ui.seed_stride_bp)
                                .desired_width(52.0),
                        )
                        .on_hover_text(
                            "Phase-1 hash-density knob. Default 1 = one seed start per base.",
                        )
                        .changed();
                    ui.label("min hit").on_hover_text(
                        "Minimum raw matched/tested seed fraction required to pass.",
                    );
                    persist_ui_state |= ui
                        .add(
                            egui::TextEdit::singleline(
                                &mut self.rna_reads_ui.min_seed_hit_fraction,
                            )
                            .desired_width(56.0),
                        )
                        .on_hover_text("Raw seed-hit threshold in [0,1], for example 0.30.")
                        .changed();
                    ui.label("min weighted").on_hover_text(
                        "Minimum occurrence-weighted seed-hit fraction required to pass.",
                    );
                    persist_ui_state |= ui
                        .add(
                            egui::TextEdit::singleline(
                                &mut self.rna_reads_ui.min_weighted_seed_hit_fraction,
                            )
                            .desired_width(56.0),
                        )
                        .on_hover_text(
                            "Weighted threshold in [0,1]; downweights highly repeated seed bits.",
                        )
                        .changed();
                    ui.label("min unique")
                        .on_hover_text("Minimum number of distinct matched seed hashes.");
                    persist_ui_state |= ui
                        .add(
                            egui::TextEdit::singleline(
                                &mut self.rna_reads_ui.min_unique_matched_kmers,
                            )
                            .desired_width(56.0),
                        )
                        .on_hover_text(
                            "Prevents low-complexity reads from passing on repetitive seeds alone.",
                        )
                        .changed();
                    ui.label("max median gap").on_hover_text(
                        "Maximum allowed median distance between matched seed positions in the inferred transcript chain.",
                    );
                    persist_ui_state |= ui
                        .add(
                            egui::TextEdit::singleline(
                                &mut self.rna_reads_ui.max_median_transcript_gap,
                            )
                            .desired_width(64.0),
                        )
                        .on_hover_text(
                            "Lower values require denser, more contiguous seed placement.",
                        )
                        .changed();
                    ui.label("min chain").on_hover_text(
                        "Minimum coherent-chain support fraction for matched seed observations.",
                    );
                    persist_ui_state |= ui
                        .add(
                            egui::TextEdit::singleline(
                                &mut self.rna_reads_ui.min_chain_consistency_fraction,
                            )
                            .desired_width(62.0),
                        )
                        .on_hover_text(
                            "Fraction in [0,1]; higher values reject dispersed local matches.",
                        )
                        .changed();
                    ui.label("min transitions").on_hover_text(
                        "Minimum number of confirmed exon-exon transitions in inferred exon path.",
                    );
                    persist_ui_state |= ui
                        .add(
                            egui::TextEdit::singleline(
                                &mut self.rna_reads_ui.min_confirmed_exon_transitions,
                            )
                            .desired_width(52.0),
                        )
                        .on_hover_text(
                            "Require at least this many confirmed junction transitions for pass.",
                        )
                        .changed();
                    ui.label("min transition frac").on_hover_text(
                        "Minimum confirmed/expected transition fraction in inferred exon path.",
                    );
                    persist_ui_state |= ui
                        .add(
                            egui::TextEdit::singleline(
                                &mut self.rna_reads_ui.min_transition_support_fraction,
                            )
                            .desired_width(62.0),
                        )
                        .on_hover_text(
                            "Fraction in [0,1] controlling transition-consistency strictness.",
                        )
                        .changed();
                    if self.rna_reads_ui.cdna_poly_t_flip_enabled {
                        ui.label("poly-T head min T-bp").on_hover_text(
                            "Minimum T support used by the tolerant 5' poly-T-head detector for cDNA normalization.",
                        );
                        persist_ui_state |= ui
                            .add(
                                egui::TextEdit::singleline(
                                    &mut self.rna_reads_ui.poly_t_prefix_min_bp,
                                )
                                .desired_width(56.0),
                            )
                            .on_hover_text(
                                "Higher values require stronger T-rich heads before reverse-complement normalization.",
                            )
                            .changed();
                    }
                });
                ui.horizontal_wrapped(|ui| {
                    ui.label("dotplot word").on_hover_text(
                        "Word size used when exporting RNA-read sequence dotplots from the read-effects/detail panel.",
                    );
                    persist_ui_state |= ui
                        .add(
                            egui::TextEdit::singleline(&mut self.dotplot_ui.word_size)
                                .desired_width(46.0),
                        )
                        .on_hover_text(
                            "Smaller values are more sensitive; larger values are stricter.",
                        )
                        .changed();
                    ui.label("dotplot step").on_hover_text(
                        "Sampling stride used by RNA-read dotplot export. Smaller values draw denser dotplots.",
                    );
                    persist_ui_state |= ui
                        .add(
                            egui::TextEdit::singleline(&mut self.dotplot_ui.step_bp)
                                .desired_width(46.0),
                        )
                        .on_hover_text(
                            "Default 1 = every possible start; larger values make the dotplot sparser and faster.",
                        )
                        .changed();
                    ui.label("dotplot mismatches (0=exact)").on_hover_text(
                        "Allowed mismatches per exported dotplot word. Keep this at 0 for exact words only; raise it only when you deliberately want a more tolerant, slower search.",
                    );
                    persist_ui_state |= ui
                        .add(
                            egui::TextEdit::singleline(&mut self.dotplot_ui.max_mismatches)
                                .desired_width(52.0),
                        )
                        .on_hover_text(
                            "0 = exact words only. Increasing this allows inexact word matches and usually slows the dotplot search.",
                        )
                        .changed();
                    ui.label("dotplot tile").on_hover_text(
                        "Optional tiling chunk size for exported RNA-read dotplots. Leave empty to avoid tiling.",
                    );
                    persist_ui_state |= ui
                        .add(
                            egui::TextEdit::singleline(&mut self.dotplot_ui.tile_bp)
                                .desired_width(56.0),
                        )
                        .on_hover_text(
                            "Usually left empty for RNA-read dotplots unless the compared spans become very large.",
                        )
                        .changed();
                });
                ui.small(
                    egui::RichText::new(
                        "These dotplot knobs are the exact settings used by `Export dotplot...` and `Export dotplots for selected reads...` in the read-effects panel.",
                    )
                    .color(egui::Color32::from_rgb(100, 116, 139)),
                );
                ui.horizontal_wrapped(|ui| {
                    ui.label("align band").on_hover_text(
                        "Phase-2 alignment band width used by `Run Alignment Phase` (and optional shell override).",
                    );
                    persist_ui_state |= ui
                        .add(
                            egui::TextEdit::singleline(
                                &mut self.rna_reads_ui.align_band_width_bp,
                            )
                            .desired_width(54.0),
                        )
                        .on_hover_text(
                            "Banded aligner width used when aligning retained reads in phase 2.",
                        )
                        .changed();
                    ui.label("min identity").on_hover_text(
                        "Minimum identity fraction for phase-2 retained-read alignment.",
                    );
                    persist_ui_state |= ui
                        .add(
                            egui::TextEdit::singleline(
                                &mut self.rna_reads_ui.align_min_identity_fraction,
                            )
                            .desired_width(56.0),
                        )
                        .on_hover_text("Alignments below this identity are discarded.")
                        .changed();
                    ui.label("max secondary").on_hover_text(
                        "Maximum number of secondary mappings kept per read in phase 2.",
                    );
                    persist_ui_state |= ui
                        .add(
                            egui::TextEdit::singleline(
                                &mut self.rna_reads_ui.align_max_secondary_mappings,
                            )
                            .desired_width(42.0),
                        )
                        .on_hover_text("Set to 0 to keep only the best mapping.")
                        .changed();
                    ui.label("align selection").on_hover_text(
                        "Which retained-hit subset from the saved report is re-aligned in phase 2. The default `seed_passed` setting is the narrower/faster rerun mode; `all retained` is broader; `already_aligned` is mainly for rerunning phase 2 on rows that already received a mapping in an earlier pass.",
                    );
                    egui::ComboBox::from_id_salt(format!(
                        "rna_read_align_selection_{}_{}",
                        view.seq_id, view.target_feature_id
                    ))
                    .selected_text(Self::rna_read_align_selection_ui_label(
                        self.rna_reads_ui.align_phase_selection,
                    ))
                    .show_ui(ui, |ui| {
                        persist_ui_state |= ui
                            .selectable_value(
                                &mut self.rna_reads_ui.align_phase_selection,
                                RnaReadHitSelection::SeedPassed,
                                "seed_passed",
                            )
                            .on_hover_text(
                                "Default: align only retained reads that currently pass the composite seed gate. If none do, phase 2 falls back to retained rows at or above raw min_hit.",
                            )
                            .changed();
                        persist_ui_state |= ui
                            .selectable_value(
                                &mut self.rna_reads_ui.align_phase_selection,
                                RnaReadHitSelection::All,
                                "all retained",
                            )
                            .on_hover_text(
                                "Align every retained row in the saved report, including rescued high-score rows that failed the composite seed gate.",
                            )
                            .changed();
                        persist_ui_state |= ui
                            .selectable_value(
                                &mut self.rna_reads_ui.align_phase_selection,
                                RnaReadHitSelection::Aligned,
                                "already_aligned",
                            )
                            .on_hover_text(
                                "Re-align only rows that already have a stored phase-2 mapping from an earlier alignment pass. This is mostly useful when you want to rerun phase 2 with different band/identity settings without broadening the working set.",
                            )
                            .changed();
                    });
                });
                ui.small(
                    egui::RichText::new(
                        "Phase-2 algorithm: reference-guided pairwise alignment against each admitted transcript template. We try banded semiglobal and local alignment first, then fall back to deterministic dense semiglobal/local alignment if the banded pass yields no hit; the best mapping is ranked by score, then identity and query coverage.",
                    )
                    .color(egui::Color32::from_rgb(100, 116, 139)),
                );
                if let Some(report) = workspace_saved_report.as_ref() {
                    ui.small(
                        egui::RichText::new(
                            Self::format_rna_read_alignment_selection_summary(
                                report.as_ref(),
                                self.rna_reads_ui.align_phase_selection,
                            ),
                        )
                        .color(egui::Color32::from_rgb(100, 116, 139)),
                    );
                }
                ui.horizontal_wrapped(|ui| {
                    if ui
                        .small_button("Prep concatemer review")
                        .on_hover_text(
                            "Set a pragmatic phase-2 review preset for concatemer/fragment audits: full report mode, all retained rows, and up to 5 secondary mappings per read.",
                        )
                        .clicked()
                    {
                        self.apply_rna_read_concatemer_review_preset();
                        persist_ui_state = true;
                    }
                    ui.small(
                        egui::RichText::new(
                            "For `nanopore_cdna_v1`, phase 1 intentionally stores `max_secondary_mappings=0`; use phase 2 with retained secondaries before treating concatemer audits as decisive.",
                        )
                        .color(egui::Color32::from_rgb(100, 116, 139)),
                    );
                });
            }
            if refresh_auto_report_id {
                persist_ui_state |= self.refresh_auto_rna_read_report_id(view);
            }
        });
        if !controls_enabled {
            ui.small(
                egui::RichText::new("Input/advanced controls are locked while a run is active.")
                    .color(egui::Color32::from_rgb(100, 116, 139)),
            );
        }
        let mut highlight_selection_update: Option<Option<usize>> = None;
        if report_id_usable_for_view {
            self.sync_rna_read_evidence_selection_to_mapping_report();
        }
        let progress_snapshot = self.current_rna_read_mapping_progress_for_view(view);
        let align_selection_label =
            Self::rna_read_align_selection_ui_label(self.rna_reads_ui.align_phase_selection)
                .to_string();
        let active_alignment_selection_summary = if self
            .rna_read_task
            .as_ref()
            .is_some_and(|task| task.operation_label == "Nanopore alignment phase")
        {
            workspace_saved_report.as_ref().map(|report| {
                Self::format_rna_read_alignment_selection_summary(
                    report.as_ref(),
                    self.rna_reads_ui.align_phase_selection,
                )
            })
        } else {
            None
        };
        if let Some(task) = &self.rna_read_task {
            let compression = Self::rna_reads_input_compression_label(&task.input_path);
            let task_label = task.operation_label.as_str();
            ui.horizontal(|ui| {
                        ui.add(egui::Spinner::new());
                        if let Some(progress) = progress_snapshot.as_ref() {
                            let elapsed_s = task.started.elapsed().as_secs_f64().max(0.001);
                            let reads_per_sec = progress.reads_processed as f64 / elapsed_s;
                            let bp_per_sec = progress.read_bases_processed as f64 / elapsed_s;
                            let seed_pass_pct = if progress.reads_processed == 0 {
                                0.0
                            } else {
                                (progress.seed_passed as f64 / progress.reads_processed as f64)
                                    * 100.0
                            };
                            let percent = if progress.reads_total == 0 {
                                0.0
                            } else {
                                (progress.reads_processed as f64 / progress.reads_total as f64)
                                    * 100.0
                            };
                            let row_eta_suffix =
                                Self::format_rna_read_progress_eta(
                                    progress.reads_processed,
                                    progress.reads_total,
                                    elapsed_s,
                                )
                                .unwrap_or_else(|| "ETA: n/a".to_string());
                            if progress.reads_total == 0 {
                                let byte_progress = if progress.input_bytes_total == 0 {
                                    "unknown".to_string()
                                } else {
                                    let fraction = (progress.input_bytes_processed as f64
                                        / progress.input_bytes_total as f64)
                                        .clamp(0.0, 1.0);
                                    format!(
                                        "{}/{} ({:.1}%)",
                                        Self::format_bytes_compact(progress.input_bytes_processed),
                                        Self::format_bytes_compact(progress.input_bytes_total),
                                        fraction * 100.0
                                    )
                                };
                                let eta_suffix = self
                                    .rna_stream_eta_text
                                    .clone()
                                    .unwrap_or_else(|| "ETA: n/a".to_string());
                                ui.label(format!(
                                    "{task_label} running: streaming FASTA ({compression}), sequences-read={}, bytes={}, {}, {:.1} reads/s, {}, len mean/med/p95={:.1}/{}/{}, seed-passed={} ({:.2}%), aligned={}, elapsed {:.1}s",
                                    progress.reads_processed,
                                    byte_progress,
                                    eta_suffix,
                                    reads_per_sec,
                                    Self::format_bp_rate_compact(bp_per_sec),
                                    progress.mean_read_length_bp,
                                    progress.median_read_length_bp,
                                    progress.p95_read_length_bp,
                                    progress.seed_passed,
                                    seed_pass_pct,
                                    progress.aligned,
                                    task.started.elapsed().as_secs_f32()
                                ));
                            } else {
                                let selection_suffix = if task.operation_label
                                    == "Nanopore alignment phase"
                                {
                                    format!(
                                        ", selection={}",
                                        align_selection_label.as_str()
                                    )
                                } else {
                                    String::new()
                                };
                                ui.label(format!(
                                    "{task_label} running: reads {}/{} ({:.1}%){}{}, input={compression}, {}, {:.1} reads/s, {}, len mean/med/p95={:.1}/{}/{}, seed-passed={} ({:.2}%), aligned={}, elapsed {:.1}s",
                                    progress.reads_processed,
                                    progress.reads_total,
                                    percent,
                                    if task.operation_label == "Nanopore alignment phase" {
                                        " retained rows"
                                    } else {
                                        ""
                                    },
                                    selection_suffix,
                                    row_eta_suffix,
                                    reads_per_sec,
                                    Self::format_bp_rate_compact(bp_per_sec),
                                    progress.mean_read_length_bp,
                                    progress.median_read_length_bp,
                                    progress.p95_read_length_bp,
                                    progress.seed_passed,
                                    seed_pass_pct,
                                    progress.aligned,
                                    task.started.elapsed().as_secs_f32()
                                ));
                            }
                        } else {
                            ui.label(format!(
                                "{task_label} running: input='{}' ({compression}), {:.1}s elapsed",
                                task.input_path,
                                task.started.elapsed().as_secs_f32()
                            ));
                        }
                    });
            ui.horizontal(|ui| {
                let cancel_requested = task.cancel_requested.load(AtomicOrdering::Relaxed);
                if cancel_requested {
                    ui.small(
                        egui::RichText::new("Cancel requested... waiting for worker to stop.")
                            .color(egui::Color32::from_rgb(220, 38, 38)),
                    );
                } else if ui
                    .button(format!("Cancel {}", task.operation_label))
                    .clicked()
                {
                    task.cancel_requested.store(true, AtomicOrdering::Relaxed);
                    self.op_status = format!("Cancel requested for {}", task.operation_label);
                }
            });
            if let Some(progress) = progress_snapshot.as_ref() {
                let percent = if progress.reads_total == 0 {
                    0.0
                } else {
                    (progress.reads_processed as f32 / progress.reads_total as f32).clamp(0.0, 1.0)
                };
                if progress.reads_total == 0 {
                    if progress.input_bytes_total > 0 {
                        let byte_fraction = (progress.input_bytes_processed as f32
                            / progress.input_bytes_total as f32)
                            .clamp(0.0, 1.0);
                        ui.add(
                                    egui::ProgressBar::new(byte_fraction)
                                        .show_percentage()
                                        .text(format!(
                                            "Streaming input FASTA ({compression}) bytes {}/{} | records-read={} (update stride: {})",
                                            Self::format_bytes_compact(progress.input_bytes_processed),
                                            Self::format_bytes_compact(progress.input_bytes_total),
                                            progress.reads_processed,
                                            progress.update_every_reads,
                                        )),
                                );
                    } else {
                        ui.add(
                                    egui::ProgressBar::new(0.0).animate(true).text(format!(
                                        "Streaming input FASTA ({compression}) records-read={} (update stride: {})",
                                        progress.reads_processed,
                                        progress.update_every_reads,
                                    )),
                                );
                    }
                } else {
                    if task.operation_label == "Nanopore alignment phase" {
                        ui.add(
                                    egui::ProgressBar::new(percent)
                                        .show_percentage()
                                        .text(format!(
                                            "Retained rows processed: {}/{} (selection={} | update stride: {})",
                                            progress.reads_processed,
                                            progress.reads_total,
                                            align_selection_label.as_str(),
                                            progress.update_every_reads,
                                        )),
                                );
                    } else {
                        ui.add(
                            egui::ProgressBar::new(percent)
                                .show_percentage()
                                .text(format!(
                                    "Reads processed: {}/{} (update stride: {})",
                                    progress.reads_processed,
                                    progress.reads_total,
                                    progress.update_every_reads,
                                )),
                        );
                    }
                }
                if task.operation_label == "Nanopore alignment phase"
                    && let Some(summary) = active_alignment_selection_summary.as_ref()
                {
                    ui.small(
                        egui::RichText::new(summary).color(egui::Color32::from_rgb(100, 116, 139)),
                    );
                }
                let matched_ratio = if progress.tested_kmers == 0 {
                    0.0
                } else {
                    progress.matched_kmers as f64 / progress.tested_kmers as f64
                };
                ui.small(format!(
                    "Cumulative seed confirmations: {}/{} ({:.3})",
                    progress.matched_kmers, progress.tested_kmers, matched_ratio
                ));
                ui.small(format!(
                    "THROUGHPUT: reads={} bases={} | mean len={:.1} bp median={} bp p95={} bp",
                    progress.reads_processed,
                    progress.read_bases_processed,
                    progress.mean_read_length_bp,
                    progress.median_read_length_bp,
                    progress.p95_read_length_bp
                ));
                if !progress.origin_class_counts.is_empty() {
                    let class_parts = progress
                        .origin_class_counts
                        .iter()
                        .map(|(class, count)| format!("{class}={count}"))
                        .collect::<Vec<_>>();
                    ui.small(format!("Origin classes: {}", class_parts.join(" | ")));
                }
                let elapsed_ms = task.started.elapsed().as_secs_f64() * 1000.0;
                let seed_ms = progress.seed_compute_ms.max(0.0);
                let align_ms = progress.align_compute_ms.max(0.0);
                let io_ms = progress.io_read_ms.max(0.0);
                let parse_ms = progress.fasta_parse_ms.max(0.0);
                let normalize_ms = progress.normalize_compute_ms.max(0.0);
                let inference_ms = progress.inference_compute_ms.max(0.0);
                let emit_ms = progress.progress_emit_ms.max(0.0);
                let overhead_ms = (elapsed_ms
                    - seed_ms
                    - align_ms
                    - io_ms
                    - parse_ms
                    - normalize_ms
                    - inference_ms
                    - emit_ms)
                    .max(0.0);
                ui.small(format!(
                            "COMPUTE: seed={:.2}s align={:.2}s io={:.2}s parse={:.2}s norm={:.2}s infer={:.2}s emit={:.2}s other={:.2}s",
                            seed_ms / 1000.0,
                            align_ms / 1000.0,
                            io_ms / 1000.0,
                            parse_ms / 1000.0,
                            normalize_ms / 1000.0,
                            inference_ms / 1000.0,
                            emit_ms / 1000.0,
                            overhead_ms / 1000.0
                        ));
                ui.small(self.rna_alignment_debug_line(progress));
                ui.horizontal(|ui| {
                    ui.small("Overlay guides:");
                    ui.checkbox(&mut self.rna_seed_overlay_show_exons, "Exons");
                    ui.checkbox(&mut self.rna_seed_overlay_show_introns, "Introns");
                    ui.checkbox(&mut self.rna_seed_overlay_exonic_coords, "Exonic coords");
                });
                self.render_rna_read_seed_histogram(
                    ui,
                    progress,
                    &self.rna_seed_catalog_preview,
                    &self.rna_seed_template_audit_preview,
                    view,
                    self.rna_seed_overlay_show_exons,
                    self.rna_seed_overlay_show_introns,
                    self.rna_seed_overlay_exonic_coords,
                );
                ui.horizontal(|ui| {
                            ui.small("Score density scale:");
                            persist_ui_state |= ui
                                .selectable_value(
                                    &mut self.rna_read_evidence_ui.score_density_use_log_scale,
                                    false,
                                    "Linear",
                                )
                                .on_hover_text(
                                    "Draw score-density bar heights proportionally to raw bin counts.",
                                )
                                .changed();
                            persist_ui_state |= ui
                                .selectable_value(
                                    &mut self.rna_read_evidence_ui.score_density_use_log_scale,
                                    true,
                                    "Log",
                                )
                                .on_hover_text(
                                    "Draw score-density bar heights using log(1+count) to emphasize sparse bins.",
                                )
                                .changed();
                            ui.separator();
                            ui.small("Population:");
                            persist_ui_state |= ui
                                .selectable_value(
                                    &mut self.rna_read_evidence_ui.score_density_variant,
                                    RnaReadScoreDensityVariant::AllScored,
                                    "All scored",
                                )
                                .on_hover_text(
                                    "Show the phase-1 score distribution across all scored reads, before the full composite seed gate is applied.",
                                )
                                .changed();
                            persist_ui_state |= ui
                                .selectable_value(
                                    &mut self.rna_read_evidence_ui.score_density_variant,
                                    RnaReadScoreDensityVariant::CompositeSeedGate,
                                    "Composite gate",
                                )
                                .on_hover_text(
                                    "Show only reads that passed the full composite seed gate under the recorded phase-1 thresholds for this run.",
                                )
                                .changed();
                            persist_ui_state |= ui
                                .selectable_value(
                                    &mut self.rna_read_evidence_ui.score_density_variant,
                                    RnaReadScoreDensityVariant::RetainedReplayCurrentControls,
                                    "Retained replay",
                                )
                                .on_hover_text(
                                    "Replay the full composite seed gate instantly over retained saved-report rows using the current controls. Reads never retained by the original run are not revisited.",
                                )
                                .changed();
                        });
                self.render_rna_read_score_density_plot(
                    ui,
                    progress,
                    self.rna_read_evidence_ui.score_density_use_log_scale,
                );
                self.render_rna_read_length_distributions_panel(
                    ui,
                    view,
                    workspace_saved_report.as_deref(),
                );
                self.render_rna_read_statistics_tabs(ui, view, progress, true);
                highlight_selection_update =
                    self.render_rna_read_top_hits_preview(ui, view, progress, true);
            }
        } else if let Some(progress) = progress_snapshot.as_ref() {
            let reads_denominator = progress.reads_total.max(progress.reads_processed).max(1);
            let seed_pass_pct = (progress.seed_passed as f64 / reads_denominator as f64) * 100.0;
            ui.small(format!(
                        "Last run: reads {}/{} | seed-passed={} ({:.2}%) | aligned={} | matched/tested={}/{}",
                        progress.reads_processed,
                        progress.reads_total,
                        progress.seed_passed,
                        seed_pass_pct,
                        progress.aligned,
                        progress.matched_kmers,
                        progress.tested_kmers
                    ));
            ui.small(format!(
                "THROUGHPUT (cumulative): bases={} | mean len={:.1} bp median={} bp p95={} bp",
                progress.read_bases_processed,
                progress.mean_read_length_bp,
                progress.median_read_length_bp,
                progress.p95_read_length_bp
            ));
            ui.small(format!(
                        "COMPUTE (cumulative): seed={:.2}s align={:.2}s io={:.2}s parse={:.2}s norm={:.2}s infer={:.2}s emit={:.2}s",
                        progress.seed_compute_ms.max(0.0) / 1000.0,
                        progress.align_compute_ms.max(0.0) / 1000.0,
                        progress.io_read_ms.max(0.0) / 1000.0,
                        progress.fasta_parse_ms.max(0.0) / 1000.0,
                        progress.normalize_compute_ms.max(0.0) / 1000.0,
                        progress.inference_compute_ms.max(0.0) / 1000.0,
                        progress.progress_emit_ms.max(0.0) / 1000.0,
                    ));
            ui.small(self.rna_alignment_debug_line(progress));
            ui.horizontal(|ui| {
                ui.small("Overlay guides:");
                ui.checkbox(&mut self.rna_seed_overlay_show_exons, "Exons");
                ui.checkbox(&mut self.rna_seed_overlay_show_introns, "Introns");
                ui.checkbox(&mut self.rna_seed_overlay_exonic_coords, "Exonic coords");
            });
            self.render_rna_read_seed_histogram(
                ui,
                progress,
                &self.rna_seed_catalog_preview,
                &self.rna_seed_template_audit_preview,
                view,
                self.rna_seed_overlay_show_exons,
                self.rna_seed_overlay_show_introns,
                self.rna_seed_overlay_exonic_coords,
            );
            ui.horizontal(|ui| {
                        ui.small("Score density scale:");
                        persist_ui_state |= ui
                            .selectable_value(
                                &mut self.rna_read_evidence_ui.score_density_use_log_scale,
                                false,
                                "Linear",
                            )
                            .on_hover_text(
                                "Draw score-density bar heights proportionally to raw bin counts.",
                            )
                            .changed();
                        persist_ui_state |= ui
                            .selectable_value(
                                &mut self.rna_read_evidence_ui.score_density_use_log_scale,
                                true,
                                "Log",
                            )
                            .on_hover_text(
                                "Draw score-density bar heights using log(1+count) to emphasize sparse bins.",
                            )
                            .changed();
                        ui.separator();
                        ui.small("Population:");
                        persist_ui_state |= ui
                            .selectable_value(
                                &mut self.rna_read_evidence_ui.score_density_variant,
                                RnaReadScoreDensityVariant::AllScored,
                                "All scored",
                            )
                            .on_hover_text(
                                "Show the phase-1 score distribution across all scored reads, before the full composite seed gate is applied.",
                            )
                            .changed();
                        persist_ui_state |= ui
                            .selectable_value(
                                &mut self.rna_read_evidence_ui.score_density_variant,
                                RnaReadScoreDensityVariant::CompositeSeedGate,
                                "Composite gate",
                            )
                            .on_hover_text(
                                "Show only reads that passed the full composite seed gate under the recorded phase-1 thresholds for this run.",
                            )
                            .changed();
                        persist_ui_state |= ui
                            .selectable_value(
                                &mut self.rna_read_evidence_ui.score_density_variant,
                                RnaReadScoreDensityVariant::RetainedReplayCurrentControls,
                                "Retained replay",
                            )
                            .on_hover_text(
                                "Replay the full composite seed gate instantly over retained saved-report rows using the current controls. Reads never retained by the original run are not revisited.",
                            )
                            .changed();
                    });
            self.render_rna_read_score_density_plot(
                ui,
                progress,
                self.rna_read_evidence_ui.score_density_use_log_scale,
            );
            self.render_rna_read_length_distributions_panel(
                ui,
                view,
                workspace_saved_report.as_deref(),
            );
            self.render_rna_read_statistics_tabs(ui, view, progress, true);
            highlight_selection_update =
                self.render_rna_read_top_hits_preview(ui, view, progress, true);
        }
        if let Some(next_selection) = highlight_selection_update {
            self.rna_seed_highlight_record_index = next_selection;
        }
        if self.rna_read_task.is_none() {
            if let Some(report) = workspace_saved_report.as_deref() {
                let report_seed_pass_pct = if report.read_count_total == 0 {
                    0.0
                } else {
                    (report.read_count_seed_passed as f64 / report.read_count_total as f64) * 100.0
                };
                ui.small(format!(
                            "Report '{}': mode={} targets={} roi_capture={} | retained_hits={} / total_reads={} | seed-passed={} ({:.2}%) | aligned={} | msa-eligible(retained)={}",
                            report.report_id,
                            report.origin_mode.as_str(),
                            report.target_gene_ids.len(),
                            report.roi_seed_capture_enabled,
                            report.hits.len(),
                            report.read_count_total,
                            report.read_count_seed_passed,
                            report_seed_pass_pct,
                            report.read_count_aligned,
                            report.retained_count_msa_eligible
                        ));
                if !report.target_gene_ids.is_empty() {
                    let preview = report
                        .target_gene_ids
                        .iter()
                        .take(8)
                        .cloned()
                        .collect::<Vec<_>>()
                        .join(", ");
                    let suffix = if report.target_gene_ids.len() > 8 {
                        format!(" (+{} more)", report.target_gene_ids.len() - 8)
                    } else {
                        String::new()
                    };
                    ui.small(format!("Target genes (report): {}{}", preview, suffix));
                }
                egui::CollapsingHeader::new("Saved report details")
                            .default_open(false)
                            .show(ui, |ui| {
                                if !report.warnings.is_empty() {
                                    for warning in report.warnings.iter().take(2) {
                                        ui.small(
                                            egui::RichText::new(format!("warning: {warning}"))
                                                .color(egui::Color32::from_rgb(180, 83, 9)),
                                        );
                                    }
                                }
                                if !report.origin_class_counts.is_empty() {
                                    let class_parts = report
                                        .origin_class_counts
                                        .iter()
                                        .map(|(class, count)| format!("{class}={count}"))
                                        .collect::<Vec<_>>();
                                    ui.small(format!("Origin classes: {}", class_parts.join(" | ")));
                                }
                                ui.collapsing("Retained read preview (top 20)", |ui| {
                                    egui::ScrollArea::vertical()
                                        .max_height(140.0)
                                        .show(ui, |ui| {
                                            for hit in report.hits.iter().take(20) {
                                                let seq_preview = hit
                                                    .sequence
                                                    .chars()
                                                    .take(48)
                                                    .collect::<String>();
                                                ui.monospace(format!(
                                                    "#{} {} score={:.3} wscore={:.4} gap-med={} gap-n={} chain={:.2}/{} class={} oconf={:.2} sconf={:.2} matched/tested={}/{} pass={} msa={} seq={}",
                                                    hit.record_index + 1,
                                                    hit.header_id,
                                                    hit.seed_hit_fraction,
                                                    hit.weighted_seed_hit_fraction,
                                                    if hit.seed_transcript_gap_count == 0 {
                                                        "na".to_string()
                                                    } else {
                                                        format!("{:.2}", hit.seed_median_transcript_gap)
                                                    },
                                                    hit.seed_transcript_gap_count,
                                                    hit.seed_chain_support_fraction,
                                                    hit.seed_chain_support_kmers,
                                                    hit.origin_class.as_str(),
                                                    hit.origin_confidence,
                                                    hit.strand_confidence,
                                                    hit.matched_kmers,
                                                    hit.tested_kmers,
                                                    hit.passed_seed_filter,
                                                    hit.msa_eligible,
                                                    seq_preview,
                                                ));
                                            }
                                        });
                                });
                            });
            }
        }
        let report_action_hover = workspace_report_view_mismatch.as_deref().unwrap_or(
            "Start asynchronous phase-1 interpretation with the current settings. Reads are optionally cDNA-normalized, scored against the admitted transcript/junction seed index, and written into the current Report ID for later inspection, alignment, and export.",
        );
        if ui
            .add_enabled(
                self.rna_read_task.is_none() && report_id_usable_for_view,
                egui::Button::new(Self::rna_read_mapping_run_button_label()),
            )
            .on_hover_text(report_action_hover)
            .clicked()
        {
            self.run_splicing_rna_read_interpretation(view);
        }
        let align_action_hover = workspace_report_view_mismatch.as_deref().unwrap_or(
            "Reopen the retained report stored under Report ID, run phase-2 pairwise alignment on the selected retained rows, and refresh mapping summaries plus exon-transition/isoform support tables. This does not reread the FASTA input unless you rerun phase 1.",
        );
        if ui
            .add_enabled(
                self.rna_read_task.is_none() && report_id_usable_for_view,
                egui::Button::new("Run alignment phase (retained report)"),
            )
            .on_hover_text(align_action_hover)
            .clicked()
        {
            self.run_splicing_rna_read_alignment_phase(view);
        }
        ui.horizontal_wrapped(|ui| {
                    if ui
                        .add_enabled(
                            self.rna_read_task.is_none() && report_id_usable_for_view,
                            egui::Button::new("Prepare Workflow Op"),
                        )
                        .on_hover_text(workspace_report_view_mismatch.as_deref().unwrap_or(
                            "Build the same InterpretRnaReads payload from this panel and stage it in Engine Ops -> Workflow runner.",
                        ))
                        .clicked()
                    {
                        self.prepare_splicing_rna_read_interpretation_workflow(view);
                    }
                    if ui
                        .add_enabled(
                            self.rna_read_task.is_none() && report_id_usable_for_view,
                            egui::Button::new("Copy Workflow JSON"),
                        )
                        .on_hover_text(workspace_report_view_mismatch.as_deref().unwrap_or(
                            "Copy a full workflow JSON object (run_id + ops) for gentle_cli workflow/shell use.",
                        ))
                        .clicked()
                    {
                        self.copy_splicing_rna_read_interpretation_workflow_json(view, ui.ctx());
                    }
                });
        if ui
                    .add_enabled(
                        self.rna_read_task.is_none(),
                        egui::Button::new("Export Seed Hash Catalog (TSV)..."),
                    )
                    .on_hover_text(
                        "Export indexed seed hashes (sequence, positions, transcript context) for auditing.",
                    )
                    .clicked()
                {
                    self.export_splicing_seed_hash_catalog(view);
                }
        if ui
            .add_enabled(
                self.rna_read_task.is_none() && report_id_usable_for_view,
                egui::Button::new("Export Retained Top Reads (FASTA)..."),
            )
            .on_hover_text(workspace_report_view_mismatch.as_deref().unwrap_or(
                "Export retained top-ranked reads with seed diagnostics in FASTA headers.",
            ))
            .clicked()
        {
            self.export_retained_rna_hits_fasta();
        }
        if ui
            .add_enabled(
                self.rna_read_task.is_none() && report_id_usable_for_view,
                egui::Button::new("Export Exon Paths (TSV)..."),
            )
            .on_hover_text(workspace_report_view_mismatch.as_deref().unwrap_or(
                "Export per-read inferred exon paths and seed metrics for downstream review.",
            ))
            .clicked()
        {
            self.export_rna_read_exon_paths_tsv();
        }
        if ui
            .add_enabled(
                self.rna_read_task.is_none() && report_id_usable_for_view,
                egui::Button::new("Export Exon Abundance (TSV)..."),
            )
            .on_hover_text(workspace_report_view_mismatch.as_deref().unwrap_or(
                "Export exon/transition support abundance aggregated across retained reads.",
            ))
            .clicked()
        {
            self.export_rna_read_exon_abundance_tsv();
        }
        if ui
                    .add_enabled(
                        self.rna_read_task.is_none() && report_id_usable_for_view,
                        egui::Button::new("Export Score Density (SVG)..."),
                    )
                    .on_hover_text(workspace_report_view_mismatch.as_deref().unwrap_or(
                        "Export the seed-hit score-density chart as SVG using the current Linear/Log scale toggle.",
                    ))
                    .clicked()
                {
                    self.export_rna_read_score_density_svg();
                }
        if ui
                    .button("Export RNA sample sheet (all reports for current sequence)...")
                    .on_hover_text(
                        "Export one TSV row per RNA-read report with summary metrics for cohort annotation.",
                    )
                    .clicked()
                {
                    let default_name = format!("{}_rna_read_sample_sheet.tsv", view.seq_id);
                    if let Some(path) = rfd::FileDialog::new()
                        .set_file_name(&default_name)
                        .add_filter("TSV", &["tsv", "txt"])
                        .save_file()
                    {
                        self.apply_operation_with_feedback(Operation::ExportRnaReadSampleSheet {
                            path: path.display().to_string(),
                            seq_id: Some(view.seq_id.clone()),
                            report_ids: vec![],
                            gene_ids: vec![],
                            complete_rule: RnaReadGeneSupportCompleteRule::Near,
                            append: false,
                        });
                    }
                }
        if persist_ui_state {
            self.save_engine_ops_state();
        }
        self.capture_rna_read_mapping_workspace_status(previous_op_status, true);
        ui.separator();
    }

    pub(super) fn render_splicing_rna_read_evidence_section(
        &mut self,
        ui: &mut egui::Ui,
        view: &SplicingExpertView,
    ) {
        let mut persist_ui_state = false;
        let summaries = self.ensure_selected_rna_read_evidence_report_for_view(view);
        ui.horizontal_wrapped(|ui| {
            ui.label(
                egui::RichText::new("RNA-read evidence")
                    .strong()
                    .color(egui::Color32::from_rgb(30, 41, 59)),
            );
            ui.separator();
            ui.label(
                egui::RichText::new(
                    "Select a saved report for this splicing group, or open the dedicated RNA-read Mapping workspace to run or rerun interpretation.",
                )
                .size(9.0)
                .color(egui::Color32::from_rgb(100, 116, 139)),
            );
        });
        ui.horizontal_wrapped(|ui| {
            let cta = ui
                .button("Open RNA-read Mapping Workspace...")
                .on_hover_text(
                    "Open the dedicated RNA-read Mapping workspace seeded from this splicing locus.",
                );
            if cta.clicked() {
                self.open_rna_read_mapping_workspace_for_view(view);
                self.focus_rna_read_mapping_workspace_view(ui.ctx(), view);
            }
            if self.active_rna_read_task_matches_splicing_view(view) {
                ui.small(
                    egui::RichText::new(
                        "A mapping task is currently active for this locus; the dedicated workspace shows run controls and cancellation.",
                    )
                    .color(egui::Color32::from_rgb(100, 116, 139)),
                );
            }
        });
        if summaries.is_empty() {
            ui.small(
                egui::RichText::new(
                    "No RNA-read report for this splicing group yet. Open RNA-read Mapping Workspace... to seed the locus and run InterpretRnaReads / AlignRnaReadReport there.",
                )
                .color(egui::Color32::from_rgb(180, 83, 9)),
            );
            ui.separator();
            return;
        }

        ui.horizontal_wrapped(|ui| {
            ui.label("Report").on_hover_text(
                "Saved RNA-read reports filtered to the current sequence and splicing seed feature.",
            );
            egui::ComboBox::from_id_salt((
                "splicing_rna_read_report_selector",
                view.seq_id.as_str(),
                view.target_feature_id,
            ))
            .selected_text(
                summaries
                    .iter()
                    .find(|row| {
                        row.report_id
                            .eq_ignore_ascii_case(self.rna_read_evidence_ui.selected_report_id.trim())
                    })
                    .map(Self::format_rna_read_report_summary_picker_label)
                    .unwrap_or_else(|| "<none>".to_string()),
            )
            .show_ui(ui, |ui| {
                let mut next_report_id = self.rna_read_evidence_ui.selected_report_id.clone();
                let mut report_selection_changed = false;
                for row in summaries.iter() {
                    report_selection_changed |= ui
                        .selectable_value(
                            &mut next_report_id,
                            row.report_id.clone(),
                            Self::format_rna_read_report_summary_picker_label(row),
                        )
                        .on_hover_text(
                            "Switch the saved RNA-read report used by this splicing evidence panel.",
                        )
                        .changed();
                }
                if report_selection_changed {
                    persist_ui_state |= self.set_selected_rna_read_evidence_report_id(next_report_id);
                }
            });
        });

        let selected_report = self.current_saved_rna_read_report();
        if let Some(report) = selected_report.as_ref() {
            let report_seed_pass_pct = if report.read_count_total == 0 {
                0.0
            } else {
                (report.read_count_seed_passed as f64 / report.read_count_total as f64) * 100.0
            };
            ui.small(format!(
                "Report '{}': mode={} targets={} roi_capture={} | retained_hits={} / total_reads={} | seed-passed={} ({:.2}%) | aligned={} | msa-eligible(retained)={}",
                report.report_id,
                report.origin_mode.as_str(),
                report.target_gene_ids.len(),
                report.roi_seed_capture_enabled,
                report.hits.len(),
                report.read_count_total,
                report.read_count_seed_passed,
                report_seed_pass_pct,
                report.read_count_aligned,
                report.retained_count_msa_eligible
            ));
            if !report.target_gene_ids.is_empty() {
                let preview = report
                    .target_gene_ids
                    .iter()
                    .take(8)
                    .cloned()
                    .collect::<Vec<_>>()
                    .join(", ");
                let suffix = if report.target_gene_ids.len() > 8 {
                    format!(" (+{} more)", report.target_gene_ids.len() - 8)
                } else {
                    String::new()
                };
                ui.small(format!("Target genes (report): {}{}", preview, suffix));
            }
        }

        if let Some(report) = selected_report.as_deref() {
            persist_ui_state |= self.render_rna_read_concatemer_review_section(ui, report);
        }

        let progress =
            self.current_rna_read_evidence_progress_for_view(view, selected_report.as_deref());
        if let Some(progress) = progress.as_deref() {
            if progress.bins.is_empty() {
                ui.small(
                    egui::RichText::new(
                        "This report-view uses saved report state. Score density, thresholded support, mapped support, and read-effects inspection are available; the live seed-position histogram is only available while a matching run is still in memory.",
                    )
                    .color(egui::Color32::from_rgb(100, 116, 139)),
                );
            }
            ui.horizontal(|ui| {
                ui.small("Overlay guides:");
                ui.checkbox(&mut self.rna_seed_overlay_show_exons, "Exons");
                ui.checkbox(&mut self.rna_seed_overlay_show_introns, "Introns");
                ui.checkbox(&mut self.rna_seed_overlay_exonic_coords, "Exonic coords");
            });
            self.render_rna_read_seed_histogram(
                ui,
                progress,
                &self.rna_seed_catalog_preview,
                &self.rna_seed_template_audit_preview,
                view,
                self.rna_seed_overlay_show_exons,
                self.rna_seed_overlay_show_introns,
                self.rna_seed_overlay_exonic_coords,
            );
            ui.horizontal(|ui| {
                ui.small("Score density scale:");
                persist_ui_state |= ui
                    .selectable_value(
                        &mut self.rna_read_evidence_ui.score_density_use_log_scale,
                        false,
                        "Linear",
                    )
                    .changed();
                persist_ui_state |= ui
                    .selectable_value(
                        &mut self.rna_read_evidence_ui.score_density_use_log_scale,
                        true,
                        "Log",
                    )
                    .changed();
                ui.separator();
                ui.small("Population:");
                persist_ui_state |= ui
                    .selectable_value(
                        &mut self.rna_read_evidence_ui.score_density_variant,
                        RnaReadScoreDensityVariant::AllScored,
                        "All scored",
                    )
                    .changed();
                persist_ui_state |= ui
                    .selectable_value(
                        &mut self.rna_read_evidence_ui.score_density_variant,
                        RnaReadScoreDensityVariant::CompositeSeedGate,
                        "Composite gate",
                    )
                    .changed();
                persist_ui_state |= ui
                    .selectable_value(
                        &mut self.rna_read_evidence_ui.score_density_variant,
                        RnaReadScoreDensityVariant::RetainedReplayCurrentControls,
                        "Retained replay",
                    )
                    .changed();
            });
            self.render_rna_read_score_density_plot(
                ui,
                progress,
                self.rna_read_evidence_ui.score_density_use_log_scale,
            );
            self.render_rna_read_length_distributions_panel(ui, view, selected_report.as_deref());
            self.render_rna_read_statistics_tabs(ui, view, progress, false);
            if let Some(next_selection) =
                self.render_rna_read_top_hits_preview(ui, view, progress, false)
            {
                self.rna_seed_highlight_record_index = next_selection;
            }
        }

        if persist_ui_state {
            self.save_engine_ops_state();
        }
        ui.separator();
    }

    pub(super) fn render_rna_read_concatemer_review_section(
        &mut self,
        ui: &mut egui::Ui,
        report: &RnaReadInterpretationReport,
    ) -> bool {
        let mut persist_ui_state = false;
        egui::CollapsingHeader::new("Concatemer review / partner census")
            .default_open(false)
            .show(ui, |ui| {
                ui.horizontal_wrapped(|ui| {
                    ui.label(
                        egui::RichText::new(
                            "Run the shared concatemer suspicion audit here, then review recurring partner genes/transcripts without leaving the RNA-read GUI.",
                        )
                        .size(9.0)
                        .color(egui::Color32::from_rgb(100, 116, 139)),
                    );
                });
                if report.align_config.max_secondary_mappings == 0 {
                    ui.small(
                        egui::RichText::new(
                            "This report was aligned with max_secondary_mappings=0, so concatemer evidence is weaker than it could be. Re-align with retained secondaries before treating negative results as decisive.",
                        )
                        .color(egui::Color32::from_rgb(180, 83, 9)),
                    );
                }

                let selected_row_count = self.selected_rna_record_indices().len();
                ui.horizontal_wrapped(|ui| {
                    ui.label("Selection");
                    persist_ui_state |= ui
                        .selectable_value(
                            &mut self.rna_read_concatemer_ui.selection,
                            RnaReadHitSelection::All,
                            "all retained",
                        )
                        .changed();
                    persist_ui_state |= ui
                        .selectable_value(
                            &mut self.rna_read_concatemer_ui.selection,
                            RnaReadHitSelection::Aligned,
                            "already aligned",
                        )
                        .changed();
                    persist_ui_state |= ui
                        .selectable_value(
                            &mut self.rna_read_concatemer_ui.selection,
                            RnaReadHitSelection::SeedPassed,
                            "seed passed",
                        )
                        .changed();
                    ui.separator();
                    ui.label("Subset");
                    egui::ComboBox::from_id_salt((
                        "rna_read_concatemer_subset_mode",
                        report.report_id.as_str(),
                    ))
                    .selected_text(Self::rna_read_concatemer_subset_mode_label(
                        self.rna_read_concatemer_ui.subset_mode,
                    ))
                    .show_ui(ui, |ui| {
                        persist_ui_state |= ui
                            .selectable_value(
                                &mut self.rna_read_concatemer_ui.subset_mode,
                                RnaReadConcatemerSubsetMode::AllMatchingRows,
                                Self::rna_read_concatemer_subset_mode_label(
                                    RnaReadConcatemerSubsetMode::AllMatchingRows,
                                ),
                            )
                            .changed();
                        persist_ui_state |= ui
                            .selectable_value(
                                &mut self.rna_read_concatemer_ui.subset_mode,
                                RnaReadConcatemerSubsetMode::SelectedRowsOnly,
                                Self::rna_read_concatemer_subset_mode_label(
                                    RnaReadConcatemerSubsetMode::SelectedRowsOnly,
                                ),
                            )
                            .changed();
                    });
                    ui.small(format!("selected rows={selected_row_count}"));
                    ui.separator();
                    ui.label("Limit");
                    persist_ui_state |= ui
                        .add(
                            egui::TextEdit::singleline(&mut self.rna_read_concatemer_ui.limit)
                                .desired_width(50.0),
                        )
                        .changed();
                    persist_ui_state |= ui
                        .checkbox(&mut self.rna_read_concatemer_ui.show_advanced, "Advanced")
                        .changed();
                });

                ui.horizontal_wrapped(|ui| {
                    ui.label("Adapter/barcode FASTA");
                    persist_ui_state |= ui
                        .add(
                            egui::TextEdit::singleline(
                                &mut self.rna_read_concatemer_ui.adapter_fasta_path,
                            )
                            .desired_width(420.0),
                        )
                        .changed();
                });
                ui.horizontal_wrapped(|ui| {
                    ui.label("Transcript FASTA(s)");
                    ui.small(
                        egui::RichText::new("one path per line")
                            .color(egui::Color32::from_rgb(100, 116, 139)),
                    );
                });
                persist_ui_state |= ui
                    .add(
                        egui::TextEdit::multiline(
                            &mut self.rna_read_concatemer_ui.transcript_fasta_paths_text,
                        )
                        .desired_rows(2)
                        .desired_width(f32::INFINITY),
                    )
                    .changed();
                ui.horizontal_wrapped(|ui| {
                    ui.label("Prepared transcript index(es)");
                    ui.small(
                        egui::RichText::new("one path per line")
                            .color(egui::Color32::from_rgb(100, 116, 139)),
                    );
                });
                persist_ui_state |= ui
                    .add(
                        egui::TextEdit::multiline(
                            &mut self.rna_read_concatemer_ui.transcript_index_paths_text,
                        )
                        .desired_rows(2)
                        .desired_width(f32::INFINITY),
                    )
                    .changed();

                if self.rna_read_concatemer_ui.show_advanced {
                    ui.horizontal_wrapped(|ui| {
                        ui.label("internal poly bp");
                        persist_ui_state |= ui
                            .add(
                                egui::TextEdit::singleline(
                                    &mut self
                                        .rna_read_concatemer_ui
                                        .internal_homopolymer_min_bp,
                                )
                                .desired_width(44.0),
                            )
                            .changed();
                        ui.label("end margin");
                        persist_ui_state |= ui
                            .add(
                                egui::TextEdit::singleline(
                                    &mut self.rna_read_concatemer_ui.end_margin_bp,
                                )
                                .desired_width(44.0),
                            )
                            .changed();
                        ui.label("max primary qcov");
                        persist_ui_state |= ui
                            .add(
                                egui::TextEdit::singleline(
                                    &mut self
                                        .rna_read_concatemer_ui
                                        .max_primary_query_coverage_fraction,
                                )
                                .desired_width(52.0),
                            )
                            .changed();
                        ui.label("min secondary id");
                        persist_ui_state |= ui
                            .add(
                                egui::TextEdit::singleline(
                                    &mut self
                                        .rna_read_concatemer_ui
                                        .min_secondary_identity_fraction,
                                )
                                .desired_width(52.0),
                            )
                            .changed();
                        ui.label("max secondary overlap");
                        persist_ui_state |= ui
                            .add(
                                egui::TextEdit::singleline(
                                    &mut self
                                        .rna_read_concatemer_ui
                                        .max_secondary_query_overlap_fraction,
                                )
                                .desired_width(52.0),
                            )
                            .changed();
                    });
                    ui.horizontal_wrapped(|ui| {
                        ui.label("adapter match bp");
                        persist_ui_state |= ui
                            .add(
                                egui::TextEdit::singleline(
                                    &mut self.rna_read_concatemer_ui.adapter_min_match_bp,
                                )
                                .desired_width(44.0),
                            )
                            .changed();
                        ui.label("fragment min bp");
                        persist_ui_state |= ui
                            .add(
                                egui::TextEdit::singleline(
                                    &mut self.rna_read_concatemer_ui.fragment_min_bp,
                                )
                                .desired_width(44.0),
                            )
                            .changed();
                        ui.label("fragment max parts");
                        persist_ui_state |= ui
                            .add(
                                egui::TextEdit::singleline(
                                    &mut self.rna_read_concatemer_ui.fragment_max_parts,
                                )
                                .desired_width(44.0),
                            )
                            .changed();
                        ui.label("fragment min id");
                        persist_ui_state |= ui
                            .add(
                                egui::TextEdit::singleline(
                                    &mut self
                                        .rna_read_concatemer_ui
                                        .fragment_min_identity_fraction,
                                )
                                .desired_width(52.0),
                            )
                            .changed();
                        ui.label("fragment min qcov");
                        persist_ui_state |= ui
                            .add(
                                egui::TextEdit::singleline(
                                    &mut self
                                        .rna_read_concatemer_ui
                                        .fragment_min_query_coverage_fraction,
                                )
                                .desired_width(52.0),
                            )
                            .changed();
                    });
                }

                let request = self.current_rna_read_concatemer_request();
                let current_cache_key = request.as_ref().ok().map(
                    |(selection, limit, selected_record_indices, settings)| {
                        Self::rna_read_concatemer_inspection_cache_key(
                            &report.report_id,
                            *selection,
                            *limit,
                            selected_record_indices,
                            settings,
                        )
                    },
                );
                let mut inspection_result = current_cache_key.as_ref().and_then(|cache_key| {
                    self.cached_rna_read_concatemer_inspection
                        .as_ref()
                        .filter(|cached| cached.cache_key == *cache_key)
                        .map(|cached| {
                            cached
                                .result
                                .as_ref()
                                .map(Arc::clone)
                                .map_err(Clone::clone)
                        })
                });

                ui.horizontal_wrapped(|ui| {
                    let run_clicked = ui
                        .button("Inspect concatemers")
                        .on_hover_text(
                            "Run the shared fragment/concatemer suspicion audit for the selected saved report and current filters.",
                        )
                        .clicked();
                    if run_clicked {
                        inspection_result = Some(match &request {
                            Ok((selection, limit, selected_record_indices, settings)) => self
                                .saved_rna_read_concatemer_inspection_for_report_id(
                                    &report.report_id,
                                    *selection,
                                    *limit,
                                    selected_record_indices.clone(),
                                    settings.clone(),
                                ),
                            Err(err) => Err(err.clone()),
                        });
                    }
                    if let Some(result) = inspection_result.as_ref()
                        && let Ok(inspection) = result
                    {
                        if ui
                            .button("Focus suspicious rows")
                            .on_hover_text(
                                "Switch the mapped cDNA read-effects view to the suspicious rows returned by this audit.",
                            )
                            .clicked()
                        {
                            self.focus_rna_read_alignment_effect_record_indices(
                                inspection
                                    .rows
                                    .iter()
                                    .map(|row| row.record_index)
                                    .collect::<Vec<_>>(),
                                "concatemer review",
                            );
                        }
                        let strong_indices = inspection
                            .rows
                            .iter()
                            .filter(|row| row.suspicion_level.as_str() == "strong")
                            .map(|row| row.record_index)
                            .collect::<Vec<_>>();
                        if ui
                            .add_enabled(
                                !strong_indices.is_empty(),
                                egui::Button::new("Focus strong rows"),
                            )
                            .on_hover_text(
                                "Focus only the strongest concatemer-suspicion rows in mapped cDNA -> Read effects.",
                            )
                            .clicked()
                        {
                            self.focus_rna_read_alignment_effect_record_indices(
                                strong_indices,
                                "strong concatemer suspicion",
                            );
                        }
                    }
                });

                if let Err(err) = &request {
                    ui.small(
                        egui::RichText::new(format!("Concatemer review settings: {err}"))
                            .color(egui::Color32::from_rgb(180, 83, 9)),
                    );
                }

                let Some(result) = inspection_result else {
                    ui.small(
                        egui::RichText::new(
                            "Click 'Inspect concatemers' to populate the suspicion ranking and recurring partner census for this report.",
                        )
                        .color(egui::Color32::from_rgb(100, 116, 139)),
                    );
                    ui.separator();
                    return;
                };

                let inspection = match result {
                    Ok(inspection) => inspection,
                    Err(err) => {
                        ui.small(
                            egui::RichText::new(format!(
                                "Concatemer review unavailable: {err}"
                            ))
                            .color(egui::Color32::from_rgb(180, 83, 9)),
                        );
                        ui.separator();
                        return;
                    }
                };

                ui.small(format!(
                    "Subset matches={} inspected={} suspicious={} strong={} | low-qcov={} internal polyA={} internal polyT={} adapter-matches={} disjoint-secondaries={} phase1-partial={} multi-gene={} | selection={} limit={}",
                    inspection.subset_match_count,
                    inspection.inspected_count,
                    inspection.suspicious_count,
                    inspection.strong_count,
                    inspection.low_query_coverage_count,
                    inspection.internal_poly_a_count,
                    inspection.internal_poly_t_count,
                    inspection.internal_adapter_match_count,
                    inspection.disjoint_secondary_mapping_count,
                    inspection.phase1_partial_origin_count,
                    inspection.multi_gene_fragment_count,
                    inspection.selection.as_str(),
                    inspection.limit,
                ));
                ui.small(format!(
                    "Partners: genes={} transcripts={} | settings: adapter={} transcript FASTA={} transcript index={}",
                    inspection.partner_gene_summaries.len(),
                    inspection.partner_transcript_summaries.len(),
                    inspection
                        .settings
                        .adapter_fasta_path
                        .as_deref()
                        .unwrap_or("<none>"),
                    inspection.settings.transcript_fasta_paths.len(),
                    inspection.settings.transcript_index_paths.len(),
                ));
                for warning in inspection.warnings.iter().take(3) {
                    ui.small(
                        egui::RichText::new(warning)
                            .color(egui::Color32::from_rgb(180, 83, 9)),
                    );
                }

                if !inspection.partner_gene_summaries.is_empty() {
                    ui.label(egui::RichText::new("Recurring partner genes").strong());
                    egui::Grid::new(format!(
                        "rna_concatemer_partner_gene_grid_{}",
                        report.report_id
                    ))
                    .striped(true)
                    .show(ui, |ui| {
                        ui.small("Rank");
                        ui.small("Gene");
                        ui.small("Suspicious reads");
                        ui.small("Fragments");
                        ui.end_row();
                        for row in inspection.partner_gene_summaries.iter().take(12) {
                            ui.small(row.rank.to_string());
                            ui.small(&row.gene_id);
                            ui.small(row.suspicious_read_count.to_string());
                            ui.small(row.fragment_count.to_string());
                            ui.end_row();
                        }
                    });
                }

                if !inspection.partner_transcript_summaries.is_empty() {
                    ui.label(egui::RichText::new("Recurring partner transcripts").strong());
                    egui::Grid::new(format!(
                        "rna_concatemer_partner_transcript_grid_{}",
                        report.report_id
                    ))
                    .striped(true)
                    .show(ui, |ui| {
                        ui.small("Rank");
                        ui.small("Transcript");
                        ui.small("Gene");
                        ui.small("Suspicious reads");
                        ui.small("Fragments");
                        ui.end_row();
                        for row in inspection.partner_transcript_summaries.iter().take(12) {
                            ui.small(row.rank.to_string());
                            ui.small(format!("{} ({})", row.transcript_id, row.transcript_label));
                            ui.small(&row.gene_id);
                            ui.small(row.suspicious_read_count.to_string());
                            ui.small(row.fragment_count.to_string());
                            ui.end_row();
                        }
                    });
                }

                ui.label(egui::RichText::new("Suspicious rows").strong());
                egui::ScrollArea::vertical()
                    .id_salt(format!("rna_concatemer_rows_{}", report.report_id))
                    .max_height(220.0)
                    .show(ui, |ui| {
                        egui::Grid::new(format!("rna_concatemer_grid_{}", report.report_id))
                            .striped(true)
                            .show(ui, |ui| {
                                ui.small("Focus");
                                ui.small("Rank");
                                ui.small("Read");
                                ui.small("Level");
                                ui.small("Best");
                                ui.small("Partners");
                                ui.small("Signals");
                                ui.end_row();
                                for row in inspection.rows.iter().take(64) {
                                    if ui.small_button("Show").clicked() {
                                        self.focus_rna_read_alignment_effect_record_indices(
                                            vec![row.record_index],
                                            "concatemer review row",
                                        );
                                    }
                                    ui.small(format!("#{}", row.rank));
                                    ui.small(format!(
                                        "{} (r{})",
                                        row.header_id,
                                        row.record_index + 1
                                    ));
                                    let level_color = match row.suspicion_level.as_str() {
                                        "strong" => egui::Color32::from_rgb(185, 28, 28),
                                        "moderate" => egui::Color32::from_rgb(180, 83, 9),
                                        "weak" => egui::Color32::from_rgb(100, 116, 139),
                                        _ => egui::Color32::from_rgb(100, 116, 139),
                                    };
                                    ui.small(
                                        egui::RichText::new(format!(
                                            "{} ({})",
                                            row.suspicion_level.as_str(),
                                            row.suspicion_score
                                        ))
                                        .color(level_color),
                                    );
                                    ui.small(
                                        row.best_gene_id
                                            .as_deref()
                                            .or(row.best_transcript_id.as_deref())
                                            .unwrap_or("<none>"),
                                    );
                                    let partner_preview = if row.partner_gene_ids.is_empty() {
                                        row.partner_transcript_ids
                                            .iter()
                                            .take(3)
                                            .cloned()
                                            .collect::<Vec<_>>()
                                            .join(", ")
                                    } else {
                                        row.partner_gene_ids
                                            .iter()
                                            .take(3)
                                            .cloned()
                                            .collect::<Vec<_>>()
                                            .join(", ")
                                    };
                                    let partner_suffix = if row.partner_gene_count > 3
                                        || row.partner_transcript_count > 3
                                    {
                                        " …"
                                    } else {
                                        ""
                                    };
                                    ui.small(format!(
                                        "{}{}",
                                        if partner_preview.is_empty() {
                                            "<none>"
                                        } else {
                                            &partner_preview
                                        },
                                        partner_suffix
                                    ));
                                    ui.small(Self::summarize_status_items(
                                        &row.suspicion_signals,
                                        3,
                                        " | ",
                                    ));
                                    ui.end_row();
                                }
                            });
                    });
                ui.separator();
            });
        persist_ui_state
    }

    pub(super) fn render_rna_read_mapping_window(&mut self, ctx: &egui::Context) {
        if !self.show_rna_read_mapping_window {
            return;
        }
        let Some(view) = self.rna_read_mapping_window_view.clone() else {
            self.show_rna_read_mapping_window = false;
            return;
        };
        let pending_initial_render = self.rna_read_mapping_window_pending_initial_render;
        self.log_rna_read_mapping_status(&view, "render begin", pending_initial_render);
        let title = Self::rna_read_mapping_window_title(&view);
        let viewport_id = Self::rna_read_mapping_viewport_id(&view.seq_id, view.target_feature_id);
        let default_size = Self::rna_read_mapping_window_default_size();
        let min_size = Self::rna_read_mapping_window_min_size();
        let content_min_size = Self::rna_read_mapping_window_content_min_size();
        let repaint_delay = Self::async_task_repaint_delay(1, false);
        if ctx.embed_viewports() {
            self.render_rna_read_mapping_embedded_window_shell(
                ctx,
                &title,
                &view,
                default_size,
                content_min_size,
                pending_initial_render,
                self.rna_read_mapping_window_focus_requested,
            );
            if self.active_rna_read_task_matches_splicing_view(&view) {
                ctx.request_repaint_after(repaint_delay);
            }
            if self.rna_read_mapping_window_focus_requested {
                self.focus_rna_read_mapping_workspace_view(ctx, &view);
                self.rna_read_mapping_window_focus_requested = false;
            }
            return;
        }
        let builder = egui::ViewportBuilder::default()
            .with_title(title.clone())
            .with_inner_size([default_size.x, default_size.y])
            .with_min_inner_size([min_size.x, min_size.y]);
        ctx.show_viewport_immediate(viewport_id, builder, |ctx, class| {
            if class == egui::ViewportClass::EmbeddedWindow {
                self.render_rna_read_mapping_embedded_window_shell(
                    ctx,
                    &title,
                    &view,
                    default_size,
                    content_min_size,
                    pending_initial_render,
                    self.rna_read_mapping_window_focus_requested,
                );
                if self.active_rna_read_task_matches_splicing_view(&view) {
                    ctx.request_repaint_after(repaint_delay);
                }
                return;
            }

            crate::egui_compat::show_central_panel(ctx, egui::CentralPanel::default(), |ui| {
                let backdrop_settings = current_window_backdrop_settings();
                paint_window_backdrop(ui, WindowBackdropKind::Splicing, &backdrop_settings);
                egui::ScrollArea::both()
                    .id_salt(format!(
                        "rna_read_mapping_scroll_viewport_{}_{}",
                        view.seq_id, view.target_feature_id
                    ))
                    .auto_shrink([false, false])
                    .show(ui, |ui| {
                        scroll_input_policy::apply_scrollarea_keyboard_navigation(
                            ui,
                            scroll_input_policy::DEFAULT_SCROLLAREA_KEYBOARD_STEP,
                        );
                        ui.set_min_size(content_min_size);
                        self.render_rna_read_mapping_window_body(
                            ctx,
                            ui,
                            &view,
                            pending_initial_render,
                        );
                    });
            });
            if self.active_rna_read_task_matches_splicing_view(&view) {
                ctx.request_repaint_after(repaint_delay);
            }

            if crate::app::GENtleApp::viewport_close_requested_or_shortcut(ctx) {
                self.show_rna_read_mapping_window = false;
            }
        });
        if self.rna_read_mapping_window_focus_requested {
            self.focus_rna_read_mapping_workspace_view(ctx, &view);
            self.rna_read_mapping_window_focus_requested = false;
        }
    }

    pub(super) fn rna_read_mapping_workspace_matches_task(
        &self,
        seq_id: &str,
        seed_feature_id: usize,
    ) -> bool {
        self.show_rna_read_mapping_window
            && self
                .rna_read_mapping_window_view
                .as_ref()
                .map(|view| view.seq_id == seq_id && view.target_feature_id == seed_feature_id)
                .unwrap_or(false)
    }

    pub(super) fn capture_rna_read_mapping_workspace_status(
        &mut self,
        previous_op_status: String,
        route_to_workspace: bool,
    ) {
        if route_to_workspace && self.op_status != previous_op_status {
            self.rna_read_mapping_status = self.op_status.clone();
            self.op_status = previous_op_status;
        }
    }

    pub(super) fn render_rna_read_mapping_status(&self, ui: &mut egui::Ui) {
        if self.rna_read_mapping_status.trim().is_empty() {
            return;
        }
        ui.group(|ui| {
            ui.horizontal_wrapped(|ui| {
                ui.label(
                    egui::RichText::new("Workspace status")
                        .strong()
                        .color(egui::Color32::from_rgb(71, 85, 105)),
                );
            });
            ui.add(
                egui::Label::new(
                    egui::RichText::new(&self.rna_read_mapping_status)
                        .monospace()
                        .size(self.feature_details_font_size()),
                )
                .wrap(),
            );
        });
    }

    pub(super) fn parse_rna_target_gene_ids(raw: &str) -> Vec<String> {
        let mut out = Vec::<String>::new();
        let mut seen = HashSet::<String>::new();
        for token in raw
            .split(|ch: char| ch == ',' || ch == ';' || ch.is_ascii_whitespace())
            .map(str::trim)
            .filter(|token| !token.is_empty())
        {
            let key = token.to_ascii_lowercase();
            if seen.insert(key) {
                out.push(token.to_string());
            }
        }
        out
    }

    pub(super) fn sanitize_workflow_run_id_component(raw: &str) -> String {
        let mut out = String::with_capacity(raw.len());
        for ch in raw.chars() {
            if ch.is_ascii_alphanumeric() || matches!(ch, '_' | '-') {
                out.push(ch.to_ascii_lowercase());
            } else if !out.ends_with('_') {
                out.push('_');
            }
        }
        out.trim_matches('_').to_string()
    }

    pub(super) fn rna_read_report_group_token(view: &SplicingExpertView) -> String {
        let group = Self::sanitize_workflow_run_id_component(&view.group_label);
        if group.is_empty() {
            let seq = Self::sanitize_workflow_run_id_component(&view.seq_id);
            if seq.is_empty() {
                "seq".to_string()
            } else {
                seq
            }
        } else {
            group
        }
    }

    pub(super) fn rna_read_report_source_token(input_path: &str) -> String {
        Path::new(input_path)
            .file_name()
            .and_then(|name| name.to_str())
            .map(|name| {
                let mut stem = name.to_string();
                for suffix in [
                    ".fasta.gz",
                    ".fa.gz",
                    ".fasta",
                    ".fa",
                    ".fastq.gz",
                    ".fq.gz",
                    ".fastq",
                    ".fq",
                    ".gz",
                ] {
                    if let Some(stripped) = stem.strip_suffix(suffix) {
                        stem = stripped.to_string();
                        break;
                    }
                }
                Self::sanitize_workflow_run_id_component(&stem)
            })
            .filter(|token| !token.is_empty())
            .unwrap_or_else(|| "reads".to_string())
    }

    pub(super) fn rna_read_profile_report_token(
        profile: RnaReadInterpretationProfile,
    ) -> &'static str {
        match profile {
            RnaReadInterpretationProfile::NanoporeCdnaV1 => "ncdna",
            RnaReadInterpretationProfile::ShortReadV1 => "short",
            RnaReadInterpretationProfile::TransposonV1 => "transposon",
        }
    }

    pub(super) fn rna_read_scope_report_token(scope: SplicingScopePreset) -> &'static str {
        match scope {
            SplicingScopePreset::AllOverlappingAnyStrand => "aoas",
            SplicingScopePreset::TargetGroupAnyStrand => "tgas",
            SplicingScopePreset::AllOverlappingTargetStrand => "aots",
            SplicingScopePreset::TargetGroupTargetStrand => "tgts",
        }
    }

    pub(super) fn rna_read_origin_report_token(origin_mode: RnaReadOriginMode) -> &'static str {
        match origin_mode {
            RnaReadOriginMode::SingleGene => "sg",
            RnaReadOriginMode::MultiGeneSparse => "mgs",
        }
    }

    pub(super) fn rna_read_target_gene_report_suffix(raw: &str) -> Option<String> {
        let genes = Self::parse_rna_target_gene_ids(raw);
        if genes.is_empty() {
            return None;
        }
        if genes.len() <= 2 {
            let compact = genes
                .iter()
                .map(|gene| Self::sanitize_workflow_run_id_component(gene))
                .filter(|token| !token.is_empty())
                .collect::<Vec<_>>();
            if !compact.is_empty() {
                return Some(format!("tg_{}", compact.join("_")));
            }
        }
        Some(format!("tg{}", genes.len()))
    }

    pub(super) fn default_rna_read_report_id(
        view: &SplicingExpertView,
        ui_state: &RnaReadInterpretOpsUiState,
    ) -> String {
        let mut parts = vec![
            "rna".to_string(),
            Self::rna_read_report_group_token(view),
            Self::rna_read_report_source_token(&ui_state.input_path),
            Self::rna_read_profile_report_token(ui_state.profile).to_string(),
            Self::rna_read_scope_report_token(ui_state.scope).to_string(),
            Self::rna_read_origin_report_token(ui_state.origin_mode).to_string(),
        ];
        if matches!(ui_state.origin_mode, RnaReadOriginMode::MultiGeneSparse)
            && let Some(target_suffix) =
                Self::rna_read_target_gene_report_suffix(&ui_state.target_gene_ids)
        {
            parts.push(target_suffix);
        }
        parts.join("_")
    }

    pub(super) fn refresh_auto_rna_read_report_id(&mut self, view: &SplicingExpertView) -> bool {
        if !self.rna_reads_ui.report_id_auto_sync {
            return false;
        }
        let suggested = Self::default_rna_read_report_id(view, &self.rna_reads_ui);
        if self.rna_reads_ui.report_id != suggested {
            self.rna_reads_ui.report_id = suggested;
            return true;
        }
        false
    }

    pub(super) fn format_rna_read_report_summary_picker_label(
        row: &RnaReadInterpretationReportSummary,
    ) -> String {
        let target_suffix = if row.target_gene_count == 0 {
            "tg=0".to_string()
        } else {
            format!("tg={}", row.target_gene_count)
        };
        format!(
            "{} | {} | {} | reads={} aligned={} | {} | {}",
            row.report_id,
            Self::rna_read_profile_report_token(row.profile),
            row.origin_mode.as_str(),
            row.read_count_total,
            row.read_count_aligned,
            Self::rna_read_scope_report_token(row.scope),
            target_suffix
        )
    }

    pub(super) fn default_rna_read_workflow_run_id(seq_id: &str, report_id: &str) -> String {
        let seq = Self::sanitize_workflow_run_id_component(seq_id);
        let report = Self::sanitize_workflow_run_id_component(report_id);
        let seq = if seq.is_empty() {
            "seq".to_string()
        } else {
            seq
        };
        let report = if report.is_empty() {
            "cdna".to_string()
        } else {
            report
        };
        format!("workflow_rna_reads_{}_{}", seq, report)
    }

    pub(super) fn apply_rna_reads_demo_specificity_preset(&mut self) {
        self.rna_reads_ui.scope = SplicingScopePreset::TargetGroupTargetStrand;
        self.rna_reads_ui.cdna_poly_t_flip_enabled = true;
        self.rna_reads_ui.kmer_len = "10".to_string();
        self.rna_reads_ui.seed_stride_bp = "1".to_string();
        self.rna_reads_ui.min_seed_hit_fraction = "0.30".to_string();
        self.rna_reads_ui.min_weighted_seed_hit_fraction = "0.07".to_string();
        self.rna_reads_ui.min_unique_matched_kmers = "16".to_string();
        self.rna_reads_ui.min_chain_consistency_fraction = "0.60".to_string();
        self.rna_reads_ui.max_median_transcript_gap = "3.0".to_string();
        self.rna_reads_ui.min_confirmed_exon_transitions = "1".to_string();
        self.rna_reads_ui.min_transition_support_fraction = "0.10".to_string();
        self.rna_reads_ui.show_advanced = true;
        self.op_status = "Applied RNA-read demo specificity preset for seed filtering".to_string();
    }

    pub(super) fn apply_rna_read_dense_similarity_preset(&mut self) {
        self.rna_reads_ui.kmer_len = "9".to_string();
        self.rna_reads_ui.seed_stride_bp = "1".to_string();
        self.dotplot_ui.word_size = "9".to_string();
        self.dotplot_ui.step_bp = "1".to_string();
        self.dotplot_ui.max_mismatches = "0".to_string();
        self.dotplot_ui.tile_bp.clear();
        self.rna_reads_ui.show_advanced = true;
        self.op_status =
            "Applied dense exact 9-mer preset for RNA-read hashing and dotplots".to_string();
    }

    pub(super) fn apply_rna_read_concatemer_review_preset(&mut self) {
        self.rna_reads_ui.align_max_secondary_mappings = "5".to_string();
        self.rna_reads_ui.align_phase_selection = RnaReadHitSelection::All;
        self.rna_reads_ui.report_mode = RnaReadReportMode::Full;
        self.rna_reads_ui.show_advanced = true;
        self.op_status = "Prepared RNA-read alignment settings for concatemer review (all retained, max secondary=5, full report mode).".to_string();
    }

    pub(super) fn rna_read_concatemer_subset_mode_label(
        mode: RnaReadConcatemerSubsetMode,
    ) -> &'static str {
        match mode {
            RnaReadConcatemerSubsetMode::AllMatchingRows => "all matching rows",
            RnaReadConcatemerSubsetMode::SelectedRowsOnly => "selected rows only",
        }
    }

    pub(super) fn parse_multiline_nonempty_paths(text: &str) -> Vec<String> {
        text.lines()
            .map(str::trim)
            .filter(|line| !line.is_empty())
            .map(ToString::to_string)
            .collect()
    }

    pub(super) fn current_rna_read_concatemer_request(
        &self,
    ) -> Result<
        (
            RnaReadHitSelection,
            usize,
            Vec<usize>,
            RnaReadConcatemerInspectionSettings,
        ),
        String,
    > {
        let limit = Self::parse_positive_usize_text(
            &self.rna_read_concatemer_ui.limit,
            "concatemer limit",
        )?;
        let selected_record_indices = match self.rna_read_concatemer_ui.subset_mode {
            RnaReadConcatemerSubsetMode::AllMatchingRows => vec![],
            RnaReadConcatemerSubsetMode::SelectedRowsOnly => {
                let selected = self.selected_rna_record_indices();
                if selected.is_empty() {
                    return Err(
                        "Select one or more mapped cDNA rows first before using 'selected rows only'"
                            .to_string(),
                    );
                }
                selected
            }
        };
        Ok((
            self.rna_read_concatemer_ui.selection,
            limit,
            selected_record_indices,
            RnaReadConcatemerInspectionSettings {
                internal_homopolymer_min_bp: Self::parse_positive_usize_text(
                    &self.rna_read_concatemer_ui.internal_homopolymer_min_bp,
                    "internal_homopolymer_min_bp",
                )?,
                end_margin_bp: Self::parse_positive_usize_text(
                    &self.rna_read_concatemer_ui.end_margin_bp,
                    "end_margin_bp",
                )?,
                max_primary_query_coverage_fraction: Self::parse_optional_f64_text(
                    &self
                        .rna_read_concatemer_ui
                        .max_primary_query_coverage_fraction,
                    "max_primary_query_coverage_fraction",
                )?
                .unwrap_or_else(|| {
                    RnaReadConcatemerInspectionSettings::default()
                        .max_primary_query_coverage_fraction
                }),
                min_secondary_identity_fraction: Self::parse_optional_f64_text(
                    &self.rna_read_concatemer_ui.min_secondary_identity_fraction,
                    "min_secondary_identity_fraction",
                )?
                .unwrap_or_else(|| {
                    RnaReadConcatemerInspectionSettings::default().min_secondary_identity_fraction
                }),
                max_secondary_query_overlap_fraction: Self::parse_optional_f64_text(
                    &self
                        .rna_read_concatemer_ui
                        .max_secondary_query_overlap_fraction,
                    "max_secondary_query_overlap_fraction",
                )?
                .unwrap_or_else(|| {
                    RnaReadConcatemerInspectionSettings::default()
                        .max_secondary_query_overlap_fraction
                }),
                adapter_fasta_path: {
                    let trimmed = self.rna_read_concatemer_ui.adapter_fasta_path.trim();
                    (!trimmed.is_empty()).then_some(trimmed.to_string())
                },
                adapter_min_match_bp: Self::parse_positive_usize_text(
                    &self.rna_read_concatemer_ui.adapter_min_match_bp,
                    "adapter_min_match_bp",
                )?,
                fragment_min_bp: Self::parse_positive_usize_text(
                    &self.rna_read_concatemer_ui.fragment_min_bp,
                    "fragment_min_bp",
                )?,
                fragment_max_parts: self
                    .rna_read_concatemer_ui
                    .fragment_max_parts
                    .trim()
                    .parse::<usize>()
                    .map_err(|_| "fragment_max_parts must be a non-negative integer".to_string())?,
                fragment_min_identity_fraction: Self::parse_optional_f64_text(
                    &self.rna_read_concatemer_ui.fragment_min_identity_fraction,
                    "fragment_min_identity_fraction",
                )?
                .unwrap_or_else(|| {
                    RnaReadConcatemerInspectionSettings::default().fragment_min_identity_fraction
                }),
                fragment_min_query_coverage_fraction: Self::parse_optional_f64_text(
                    &self
                        .rna_read_concatemer_ui
                        .fragment_min_query_coverage_fraction,
                    "fragment_min_query_coverage_fraction",
                )?
                .unwrap_or_else(|| {
                    RnaReadConcatemerInspectionSettings::default()
                        .fragment_min_query_coverage_fraction
                }),
                transcript_fasta_paths: Self::parse_multiline_nonempty_paths(
                    &self.rna_read_concatemer_ui.transcript_fasta_paths_text,
                ),
                transcript_index_paths: Self::parse_multiline_nonempty_paths(
                    &self.rna_read_concatemer_ui.transcript_index_paths_text,
                ),
            },
        ))
    }

    pub(super) fn reset_rna_read_dotplot_parameters_to_defaults(&mut self) {
        let defaults = DotplotOpsUiState::default();
        self.dotplot_ui.word_size = defaults.word_size;
        self.dotplot_ui.step_bp = defaults.step_bp;
        self.dotplot_ui.max_mismatches = defaults.max_mismatches;
        self.dotplot_ui.tile_bp = defaults.tile_bp;
        self.op_status =
            "Reset RNA-read dotplot parameters to the shared dotplot defaults".to_string();
    }

    pub(super) fn describe_ordered_window_overlap(window_len: usize, step: usize) -> String {
        let safe_step = step.max(1);
        let overlap_bp = window_len.saturating_sub(safe_step);
        let windows_per_base = (window_len.saturating_add(safe_step).saturating_sub(1)) / safe_step;
        if safe_step == 1 {
            format!(
                "adjacent windows overlap by {overlap_bp} bp; each interior base participates in {windows_per_base} consecutive ordered windows (dense sliding)"
            )
        } else if safe_step < window_len {
            format!(
                "adjacent windows overlap by {overlap_bp} bp; each interior base participates in up to {windows_per_base} consecutive ordered windows (subsampled sliding)"
            )
        } else if safe_step == window_len {
            "adjacent windows do not overlap; each interior base participates in at most 1 ordered window (edge-touching sampling)".to_string()
        } else {
            format!(
                "adjacent windows do not overlap and can leave up to {} bp unsampled between starts; each interior base participates in at most 1 ordered window (sparse sampling)",
                safe_step - window_len
            )
        }
    }

    pub(super) fn rna_read_hash_parameter_summary(&self) -> String {
        match (
            self.rna_reads_ui.kmer_len.trim().parse::<usize>(),
            self.rna_reads_ui.seed_stride_bp.trim().parse::<usize>(),
        ) {
            (Ok(kmer_len), Ok(seed_stride_bp)) if kmer_len > 0 && seed_stride_bp > 0 => format!(
                "Hashing now: k={} stride={} | {}",
                kmer_len,
                seed_stride_bp,
                Self::describe_ordered_window_overlap(kmer_len, seed_stride_bp),
            ),
            _ => format!(
                "Hashing now: k={} stride={} | enter valid integers to compute overlap/order density",
                self.rna_reads_ui.kmer_len.trim(),
                self.rna_reads_ui.seed_stride_bp.trim(),
            ),
        }
    }

    pub(super) fn rna_read_dotplot_parameter_summary(&self) -> String {
        let tile = self.dotplot_ui.tile_bp.trim();
        match (
            self.dotplot_ui.word_size.trim().parse::<usize>(),
            self.dotplot_ui.step_bp.trim().parse::<usize>(),
        ) {
            (Ok(word_size), Ok(step_bp)) if word_size > 0 && step_bp > 0 => format!(
                "RNA-read dotplots now: word={} step={} mismatches={} ({}) tile={} | {}",
                word_size,
                step_bp,
                self.dotplot_ui.max_mismatches.trim(),
                if self.dotplot_ui.max_mismatches.trim() == "0" {
                    "exact words only"
                } else {
                    "inexact words allowed"
                },
                if tile.is_empty() { "off" } else { tile },
                Self::describe_ordered_window_overlap(word_size, step_bp),
            ),
            _ => format!(
                "RNA-read dotplots now: word={} step={} mismatches={} tile={} | enter valid integers to compute overlap/order density",
                self.dotplot_ui.word_size.trim(),
                self.dotplot_ui.step_bp.trim(),
                self.dotplot_ui.max_mismatches.trim(),
                if tile.is_empty() { "off" } else { tile }
            ),
        }
    }

    pub(super) fn splicing_scope_description(
        scope: SplicingScopePreset,
    ) -> (&'static str, &'static str, &'static str) {
        match scope {
            SplicingScopePreset::AllOverlappingAnyStrand => (
                "All overlapping transcripts on any strand are indexed.",
                "Includes every annotated transcript lane overlapping the ROI, including antisense/opposite-strand genes relative to the selected target gene/group.",
                "Strand note: scoring uses the union of indexed templates from the target-gene strand and antisense/opposite strand; if a seed exists in both orientations, it can contribute in this mode.",
            ),
            SplicingScopePreset::TargetGroupAnyStrand => (
                "Only the target transcript group is indexed (any strand allowed).",
                "Restricts indexing to the selected target group but keeps every annotated strand direction if the group contains them.",
                "Strand note: scoring still uses an any-strand union within the selected group.",
            ),
            SplicingScopePreset::AllOverlappingTargetStrand => (
                "All overlapping transcripts on the target-gene strand are indexed.",
                "Includes all overlapping transcripts but only those matching the selected target gene/group's annotated strand.",
                "Strand note: antisense/opposite-strand seeds are excluded from the score in this mode.",
            ),
            SplicingScopePreset::TargetGroupTargetStrand => (
                "Only the target group on the target-gene strand is indexed.",
                "Most restrictive mode: selected group and selected target gene/group strand only.",
                "Strand note: antisense/opposite-strand seeds cannot contribute to the score in this mode.",
            ),
        }
    }

    pub(super) fn rna_reads_input_compression_label(path: &str) -> &'static str {
        if path.trim().to_ascii_lowercase().ends_with(".gz") {
            "gzip"
        } else {
            "plain"
        }
    }

    pub(super) fn refresh_rna_seed_catalog_preview(
        &mut self,
        seq_id: &str,
        seed_feature_id: usize,
        scope: SplicingScopePreset,
        kmer_len: usize,
    ) {
        let Some(engine) = self.engine.clone() else {
            self.rna_seed_catalog_preview.clear();
            self.rna_seed_template_audit_preview.clear();
            return;
        };
        let mut seed_filter = RnaReadSeedFilterConfig::default();
        seed_filter.kmer_len = kmer_len;
        let preview = match engine.read() {
            Ok(guard) => guard
                .collect_rna_seed_hash_catalog(seq_id, seed_feature_id, scope, &seed_filter)
                .and_then(|rows| {
                    guard
                        .collect_rna_seed_hash_template_audit(
                            seq_id,
                            seed_feature_id,
                            scope,
                            &seed_filter,
                        )
                        .map(|templates| (rows, templates))
                }),
            Err(_) => Err(EngineError {
                code: ErrorCode::Internal,
                message: "Engine lock poisoned while preparing RNA seed-hash preview".to_string(),
            }),
        };
        match preview {
            Ok((rows, templates)) => {
                self.rna_seed_catalog_preview = rows;
                self.rna_seed_template_audit_preview = templates;
            }
            Err(err) => {
                self.rna_seed_catalog_preview.clear();
                self.rna_seed_template_audit_preview.clear();
                self.op_status = format!("RNA seed-hash preview unavailable: {}", err.message);
            }
        }
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

    pub(super) fn format_bp_rate_compact(bp_per_sec: f64) -> String {
        let value = bp_per_sec.max(0.0);
        if value >= 1_000_000_000.0 {
            format!("{:.2} Gbp/s", value / 1_000_000_000.0)
        } else if value >= 1_000_000.0 {
            format!("{:.2} Mbp/s", value / 1_000_000.0)
        } else if value >= 1_000.0 {
            format!("{:.2} Kbp/s", value / 1_000.0)
        } else {
            format!("{value:.0} bp/s")
        }
    }

    pub(super) fn format_count_compact_km(count: u64) -> String {
        if count >= 1_000_000 {
            let millions = count as f64 / 1_000_000.0;
            let mut label = if millions < 10.0 {
                format!("{millions:.1}M")
            } else {
                format!("{millions:.0}M")
            };
            if label.len() > 4 {
                label = format!("{:.0}M", millions.min(999.0));
            }
            label
        } else if count >= 1_000 {
            let thousands = count as f64 / 1_000.0;
            let mut label = if thousands < 10.0 {
                format!("{thousands:.1}k")
            } else {
                format!("{thousands:.0}k")
            };
            if label.len() > 4 {
                label = format!("{:.0}k", thousands.min(999.0));
            }
            label
        } else {
            count.to_string()
        }
    }

    pub(super) fn format_duration_compact(seconds: f64) -> String {
        let mut total = seconds.max(0.0).round() as u64;
        let hours = total / 3600;
        total %= 3600;
        let minutes = total / 60;
        let secs = total % 60;
        if hours > 0 {
            format!("{hours}h {minutes:02}m {secs:02}s")
        } else if minutes > 0 {
            format!("{minutes}m {secs:02}s")
        } else {
            format!("{secs}s")
        }
    }

    pub(super) fn format_rna_read_progress_eta(
        reads_processed: usize,
        reads_total: usize,
        elapsed_s: f64,
    ) -> Option<String> {
        if reads_total == 0 || reads_processed == 0 {
            return None;
        }
        let elapsed_s = elapsed_s.max(0.001);
        let fraction_done = (reads_processed as f64 / reads_total as f64).clamp(0.0, 1.0);
        if fraction_done <= 0.0 {
            return None;
        }
        let estimated_total_s = elapsed_s / fraction_done;
        Some(format!(
            "ETA: {}",
            Self::format_duration_compact((estimated_total_s - elapsed_s).max(0.0))
        ))
    }

    pub(super) fn rna_read_align_selection_ui_label(
        selection: RnaReadHitSelection,
    ) -> &'static str {
        match selection {
            RnaReadHitSelection::SeedPassed => "seed_passed",
            RnaReadHitSelection::All => "all retained",
            RnaReadHitSelection::Aligned => "already_aligned",
        }
    }

    pub(super) fn rna_seed_base_to_bits(base: u8) -> Option<u32> {
        match base.to_ascii_uppercase() {
            b'A' => Some(0),
            b'C' => Some(1),
            b'G' => Some(2),
            b'T' | b'U' => Some(3),
            _ => None,
        }
    }

    pub(super) fn encode_rna_seed_bits(window: &[u8]) -> Option<u32> {
        let mut bits = 0u32;
        for base in window {
            let value = Self::rna_seed_base_to_bits(*base)?;
            bits = (bits << 2) | value;
        }
        Some(bits)
    }

    pub(super) fn collect_read_seed_bit_counts(
        sequence: &str,
        kmer_len: usize,
        seed_stride_bp: usize,
    ) -> HashMap<u32, usize> {
        let bytes = sequence.as_bytes();
        if kmer_len == 0 || bytes.len() < kmer_len {
            return HashMap::new();
        }
        let stride = seed_stride_bp.max(1);
        let mut out = HashMap::<u32, usize>::new();
        for start in (0..=bytes.len() - kmer_len).step_by(stride) {
            if let Some(bits) = Self::encode_rna_seed_bits(&bytes[start..start + kmer_len]) {
                *out.entry(bits).or_insert(0) += 1;
            }
        }
        out
    }

    pub(super) fn merged_splicing_exon_ranges(view: &SplicingExpertView) -> Vec<(usize, usize)> {
        let mut ranges = view
            .unique_exons
            .iter()
            .map(|row| {
                (
                    row.start_1based.min(row.end_1based),
                    row.start_1based.max(row.end_1based),
                )
            })
            .collect::<Vec<_>>();
        ranges.sort_by_key(|row| (row.0, row.1));
        let mut merged = Vec::<(usize, usize)>::new();
        for (start, end) in ranges {
            if let Some(last) = merged.last_mut() {
                if start <= last.1.saturating_add(1) {
                    last.1 = last.1.max(end);
                    continue;
                }
            }
            merged.push((start, end));
        }
        merged
    }

    pub(super) fn genomic_to_exonic_pos_1based(
        merged_exons: &[(usize, usize)],
        genomic_pos_1based: usize,
    ) -> Option<usize> {
        let mut offset = 0usize;
        for (start, end) in merged_exons {
            if genomic_pos_1based >= *start && genomic_pos_1based <= *end {
                return Some(
                    offset
                        .saturating_add(genomic_pos_1based.saturating_sub(*start))
                        .saturating_add(1),
                );
            }
            offset = offset.saturating_add(end.saturating_sub(*start).saturating_add(1));
        }
        None
    }

    pub(super) fn project_genomic_interval_to_exonic(
        merged_exons: &[(usize, usize)],
        interval_start_1based: usize,
        interval_end_1based: usize,
    ) -> Vec<(usize, usize)> {
        let start = interval_start_1based.min(interval_end_1based);
        let end = interval_start_1based.max(interval_end_1based);
        let mut out = Vec::<(usize, usize)>::new();
        let mut offset = 0usize;
        for (exon_start, exon_end) in merged_exons {
            let overlap_start = start.max(*exon_start);
            let overlap_end = end.min(*exon_end);
            if overlap_start <= overlap_end {
                let exo_start = offset
                    .saturating_add(overlap_start.saturating_sub(*exon_start))
                    .saturating_add(1);
                let exo_end = offset
                    .saturating_add(overlap_end.saturating_sub(*exon_start))
                    .saturating_add(1);
                out.push((exo_start, exo_end));
            }
            offset = offset.saturating_add(exon_end.saturating_sub(*exon_start).saturating_add(1));
        }
        out
    }

    pub(super) fn selected_rna_top_hit_preview<'a>(
        &self,
        progress: &'a RnaReadInterpretProgress,
    ) -> Option<&'a RnaReadTopHitPreview> {
        let selected_index = self.rna_seed_highlight_record_index?;
        progress
            .top_hits_preview
            .iter()
            .find(|row| row.record_index == selected_index)
    }

    pub(super) fn summarize_status_items(
        items: &[String],
        max_visible: usize,
        separator: &str,
    ) -> String {
        if items.is_empty() {
            return "-".to_string();
        }
        let max_visible = max_visible.max(1);
        let visible = items.iter().take(max_visible).cloned().collect::<Vec<_>>();
        if items.len() > max_visible {
            format!(
                "{} (+{} more)",
                visible.join(separator),
                items.len() - max_visible
            )
        } else {
            visible.join(separator)
        }
    }

    pub(super) fn set_transcript_derivation_status(
        &mut self,
        source_label: &str,
        result: &OpResult,
    ) {
        let created_count = result.created_seq_ids.len();
        let mut status_lines = if created_count == 1 {
            let derived_seq_id = result
                .created_seq_ids
                .first()
                .cloned()
                .unwrap_or_else(|| "<derived transcript>".to_string());
            vec![format!(
                "Derived transcript '{}' from {} and switched this window to the new cDNA sequence.",
                derived_seq_id, source_label
            )]
        } else {
            vec![
                format!(
                    "Derived {created_count} transcript sequence(s) from {}. Engine Ops pool/export targets now point to this set.",
                    source_label
                ),
                format!(
                    "created: {}",
                    Self::summarize_status_items(&result.created_seq_ids, 3, ", ")
                ),
            ]
        };
        if !result.warnings.is_empty() {
            status_lines.push(format!(
                "warnings: {}",
                Self::summarize_status_items(&result.warnings, 1, " | ")
            ));
        }
        self.op_status = status_lines.join("\n");
        self.op_error_popup = None;
    }

    pub(super) fn derive_transcript_sequences_with_compact_feedback(
        &mut self,
        seq_id: String,
        feature_ids: Vec<usize>,
        scope: Option<SplicingScopePreset>,
        source_label: &str,
    ) -> Option<OpResult> {
        let result =
            self.apply_operation_with_feedback_and_result(Operation::DeriveTranscriptSequences {
                seq_id: seq_id.clone(),
                feature_ids,
                scope,
                output_prefix: Some(format!("{seq_id}__mrna")),
            })?;
        self.set_transcript_derivation_status(source_label, &result);
        Some(result)
    }

    pub(super) fn derive_transcript_sequence_for_dotplot_reference(
        &mut self,
        seq_id: String,
        transcript_feature_id: usize,
        source_label: &str,
    ) -> Option<String> {
        let Some(engine) = self.engine.clone() else {
            self.op_status = "No engine attached".to_string();
            return None;
        };
        let result = {
            let Ok(mut guard) = engine.write() else {
                self.op_status =
                    "Engine lock poisoned while deriving transcript reference".to_string();
                return None;
            };
            match guard.apply(Operation::DeriveTranscriptSequences {
                seq_id: seq_id.clone(),
                feature_ids: vec![transcript_feature_id],
                scope: None,
                output_prefix: Some(format!("{seq_id}__mrna")),
            }) {
                Ok(result) => result,
                Err(error) => {
                    self.op_status = error.message.clone();
                    self.op_error_popup = Some(error.message);
                    return None;
                }
            }
        };
        if !result.created_seq_ids.is_empty() {
            self.last_created_seq_ids = result.created_seq_ids.clone();
            self.export_pool_inputs_text = self.last_created_seq_ids.join(", ");
        }
        let Some(derived_seq_id) = result.created_seq_ids.first().cloned() else {
            self.op_status = format!(
                "Transcript derivation from {} completed without a created reference sequence",
                source_label
            );
            return None;
        };
        let mut status_lines = vec![format!(
            "Derived transcript '{}' from {} for dotplot reference use. The active query window was left unchanged.",
            derived_seq_id, source_label
        )];
        if !result.warnings.is_empty() {
            status_lines.push(format!(
                "warnings: {}",
                Self::summarize_status_items(&result.warnings, 1, " | ")
            ));
        }
        self.op_status = status_lines.join("\n");
        self.op_error_popup = None;
        Some(derived_seq_id)
    }

    pub(super) fn ensure_materialized_rna_read_dotplot_query_sequence(
        &mut self,
        hit: &RnaReadInterpretationHit,
    ) -> Option<String> {
        let report_id = self.rna_reads_ui.report_id.trim().to_string();
        if report_id.is_empty() {
            self.op_status =
                "Load/save a Report ID before opening RNA-read dotplot workspace".to_string();
            return None;
        }
        let Some(engine) = self.engine.clone() else {
            self.op_status = "No engine attached".to_string();
            return None;
        };
        let report_token = Self::sanitize_workflow_run_id_component(&report_id);
        let header_token = Self::sanitize_export_name_component(&hit.header_id, "read");
        let preferred_seq_id = format!(
            "{}_dotplot_r{}_{}",
            report_token,
            hit.record_index + 1,
            header_token
        );
        if engine
            .read()
            .ok()
            .and_then(|guard| {
                guard
                    .state()
                    .sequences
                    .contains_key(&preferred_seq_id)
                    .then_some(())
            })
            .is_some()
        {
            return Some(preferred_seq_id);
        }
        let mut guard = match engine.write() {
            Ok(guard) => guard,
            Err(_) => {
                self.op_status =
                    "Engine lock poisoned while materializing RNA-read dotplot query".to_string();
                return None;
            }
        };
        let result = match guard.apply(Operation::MaterializeRnaReadHitSequences {
            report_id: report_id.clone(),
            selection: RnaReadHitSelection::All,
            selected_record_indices: vec![hit.record_index],
            output_prefix: Some(format!("{report_token}_dotplot")),
        }) {
            Ok(result) => result,
            Err(error) => {
                self.op_status = error.message.clone();
                self.op_error_popup = Some(error.message);
                return None;
            }
        };
        let seq_id = result
            .created_seq_ids
            .first()
            .cloned()
            .unwrap_or(preferred_seq_id);
        self.last_created_seq_ids = result.created_seq_ids.clone();
        self.export_pool_inputs_text = self.last_created_seq_ids.join(", ");
        Some(seq_id)
    }

    pub(super) fn open_rna_read_dotplot_workspace(
        &mut self,
        view: &SplicingExpertView,
        hit: &RnaReadInterpretationHit,
    ) {
        let Some(query_seq_id) = self.ensure_materialized_rna_read_dotplot_query_sequence(hit)
        else {
            return;
        };
        let reference_span_start_0based = view.region_start_1based.saturating_sub(1);
        let reference_span_end_0based = view.region_end_1based;
        if reference_span_end_0based <= reference_span_start_0based {
            self.op_status = "Splicing ROI is invalid for RNA-read dotplot workspace".to_string();
            return;
        }
        let query_len = self
            .engine
            .as_ref()
            .and_then(|engine| {
                engine.read().ok().and_then(|guard| {
                    guard
                        .state()
                        .sequences
                        .get(&query_seq_id)
                        .map(|dna| dna.len())
                })
            })
            .unwrap_or(hit.sequence.len());
        if query_len == 0 {
            self.op_status = format!(
                "RNA-read dotplot query '{}' is empty; workspace not opened",
                query_seq_id
            );
            return;
        }
        let requested_word_size =
            Self::parse_positive_usize_text(&self.dotplot_ui.word_size, "dotplot word_size")
                .unwrap_or(7);
        let requested_step_bp =
            Self::parse_positive_usize_text(&self.dotplot_ui.step_bp, "dotplot step_bp")
                .unwrap_or(1);
        let reference_span_bp = reference_span_end_0based
            .saturating_sub(reference_span_start_0based)
            .max(1);
        let adjusted_step_bp = Self::recommend_pair_dotplot_step(
            query_len,
            reference_span_bp,
            requested_word_size,
            requested_step_bp,
        );
        self.dotplot_query_override_seq_id = query_seq_id.clone();
        self.dotplot_query_override_source_label = format!(
            "RNA-read #{} {} from report {}",
            hit.record_index + 1,
            hit.header_id,
            self.rna_reads_ui.report_id.trim()
        );
        self.dotplot_ui.mode = Self::rna_read_dotplot_mode_from_hit(hit);
        self.dotplot_ui.reference_seq_id = view.seq_id.clone();
        self.dotplot_ui.reference_span_start_0based = reference_span_start_0based.to_string();
        self.dotplot_ui.reference_span_end_0based = reference_span_end_0based.to_string();
        self.dotplot_ui.step_bp = adjusted_step_bp.to_string();
        self.dotplot_ui.dotplot_id = format!(
            "{}_vs_{}_r{}",
            Self::normalize_operation_id_token(&query_seq_id),
            Self::normalize_operation_id_token(&view.seq_id),
            hit.record_index + 1
        );
        self.dotplot_locked_crosshair_bp = None;
        self.dotplot_hover_crosshair_bp = None;
        self.invalidate_dotplot_cache();
        self.compute_primary_dotplot();
        self.open_dotplot_window();
        self.op_status = if adjusted_step_bp != requested_step_bp {
            format!(
                "Opened interactive dotplot workspace for RNA-read #{} with auto-adjusted step {} (requested {}).",
                hit.record_index + 1,
                adjusted_step_bp,
                requested_step_bp
            )
        } else {
            format!(
                "Opened interactive dotplot workspace for RNA-read #{} ({})",
                hit.record_index + 1,
                hit.header_id
            )
        };
    }

    pub(super) fn rna_read_dotplot_mode_from_hit(hit: &RnaReadInterpretationHit) -> DotplotMode {
        let reverse_complement_against_genome = hit
            .best_mapping
            .as_ref()
            .map(|mapping| {
                let transcript_reverse = mapping.strand.trim() == "-";
                transcript_reverse ^ mapping.query_reverse_complemented
            })
            .unwrap_or_else(|| hit.strand_diagnostics.selected_strand.trim() == "-");
        if reverse_complement_against_genome {
            DotplotMode::PairReverseComplement
        } else {
            DotplotMode::PairForward
        }
    }

    pub(super) fn current_rna_read_alignment_display(
        &mut self,
        record_index: usize,
    ) -> Result<Arc<RnaReadAlignmentDisplay>, String> {
        let Some(report) = self.current_saved_rna_read_report() else {
            return Err(
                "Load/save a Report ID before inspecting phase-2 pairwise alignment".to_string(),
            );
        };
        self.saved_rna_read_alignment_display_for_report_id(&report.report_id, record_index)
    }

    pub(super) fn dotplot_parameter_tag(
        word_size: usize,
        step_bp: usize,
        max_mismatches: usize,
        tile_bp: Option<usize>,
    ) -> String {
        let tile = tile_bp
            .map(|value| value.to_string())
            .unwrap_or_else(|| "auto".to_string());
        format!("w{word_size}_s{step_bp}_mm{max_mismatches}_tile{tile}")
    }

    pub(super) fn resolve_rna_read_sequence_dotplot_word_and_step(
        query_len: usize,
        reference_span_bp: usize,
        requested_word_size: usize,
        requested_step_bp: usize,
    ) -> (usize, usize) {
        let word_size = requested_word_size
            .min(query_len.min(reference_span_bp))
            .max(1);
        let step_bp = Self::recommend_pair_dotplot_step(
            query_len,
            reference_span_bp,
            word_size,
            requested_step_bp,
        );
        (word_size, step_bp)
    }

    pub(super) fn default_rna_read_sequence_dotplot_svg_file_name(
        report_id: &str,
        hit: &RnaReadInterpretationHit,
        word_size: usize,
        step_bp: usize,
        max_mismatches: usize,
        tile_bp: Option<usize>,
    ) -> String {
        let report_token = Self::sanitize_export_name_component(report_id, "report");
        let header_token = Self::sanitize_export_name_component(&hit.header_id, "read");
        let params = Self::dotplot_parameter_tag(word_size, step_bp, max_mismatches, tile_bp);
        format!(
            "rna_read_dotplot_{}_r{}_{}_{}.svg",
            report_token,
            hit.record_index + 1,
            header_token,
            params
        )
    }

    pub(super) fn export_rna_read_sequence_dotplot_svgs(
        &mut self,
        view: &SplicingExpertView,
        hits: Vec<&RnaReadInterpretationHit>,
        output_paths: Vec<PathBuf>,
    ) {
        if hits.is_empty() {
            self.op_status = "No RNA-read hits selected for dotplot export".to_string();
            return;
        }
        if hits.len() != output_paths.len() {
            self.op_status =
                "RNA-read dotplot export path count did not match selected hit count".to_string();
            return;
        }
        let Some(engine) = self.engine.clone() else {
            self.op_status = "No engine attached".to_string();
            return;
        };
        let report_id = self.rna_reads_ui.report_id.trim().to_string();
        if report_id.is_empty() {
            self.op_status = "RNA-read dotplot export requires a saved Report ID".to_string();
            return;
        }
        let reference_seq_id = view.seq_id.clone();
        let reference_span_start_0based = view.region_start_1based.saturating_sub(1);
        let reference_span_end_0based = view.region_end_1based;
        let reference_span_bp = reference_span_end_0based
            .saturating_sub(reference_span_start_0based)
            .max(1);
        let requested_word_size =
            Self::parse_positive_usize_text(&self.dotplot_ui.word_size, "dotplot word_size")
                .unwrap_or(7);
        let requested_step_bp =
            Self::parse_positive_usize_text(&self.dotplot_ui.step_bp, "dotplot step_bp")
                .unwrap_or(1);
        let max_mismatches = Self::parse_optional_usize_text(
            &self.dotplot_ui.max_mismatches,
            "dotplot max_mismatches",
        )
        .ok()
        .flatten()
        .unwrap_or(0);
        let tile_bp = Self::parse_optional_usize_text(&self.dotplot_ui.tile_bp, "dotplot tile_bp")
            .ok()
            .flatten()
            .filter(|value| *value > 0);

        let explicit_indices = hits.iter().map(|hit| hit.record_index).collect::<Vec<_>>();
        let output_prefix = Some(format!(
            "{}_dotplot_reads",
            Self::sanitize_workflow_run_id_component(&report_id)
        ));
        let mut guard = match engine.write() {
            Ok(guard) => guard,
            Err(_) => {
                self.op_status =
                    "Engine lock poisoned while exporting RNA-read dotplots".to_string();
                return;
            }
        };
        let materialize_result = match guard.apply(Operation::MaterializeRnaReadHitSequences {
            report_id: report_id.clone(),
            selection: RnaReadHitSelection::All,
            selected_record_indices: explicit_indices,
            output_prefix,
        }) {
            Ok(result) => result,
            Err(error) => {
                self.op_status = error.message.clone();
                self.op_error_popup = Some(error.message);
                return;
            }
        };
        if materialize_result.created_seq_ids.len() != hits.len() {
            self.op_status = format!(
                "RNA-read dotplot export expected {} created sequence(s) but got {}",
                hits.len(),
                materialize_result.created_seq_ids.len()
            );
            return;
        }
        self.last_created_seq_ids = materialize_result.created_seq_ids.clone();
        self.export_pool_inputs_text = self.last_created_seq_ids.join(", ");

        for ((hit, created_seq_id), output_path) in hits
            .iter()
            .zip(materialize_result.created_seq_ids.iter())
            .zip(output_paths.iter())
        {
            let query_len = hit.sequence.len();
            if query_len == 0 {
                continue;
            }
            let (word_size, step_bp) = Self::resolve_rna_read_sequence_dotplot_word_and_step(
                query_len,
                reference_span_bp,
                requested_word_size,
                requested_step_bp,
            );
            let dotplot_id = format!(
                "{}_vs_{}_r{}",
                Self::normalize_operation_id_token(created_seq_id),
                Self::normalize_operation_id_token(&reference_seq_id),
                hit.record_index + 1
            );
            if let Err(error) = guard.apply(Operation::ComputeDotplot {
                seq_id: created_seq_id.clone(),
                reference_seq_id: Some(reference_seq_id.clone()),
                span_start_0based: Some(0),
                span_end_0based: Some(query_len),
                reference_span_start_0based: Some(reference_span_start_0based),
                reference_span_end_0based: Some(reference_span_end_0based),
                mode: Self::rna_read_dotplot_mode_from_hit(hit),
                word_size,
                step_bp,
                max_mismatches,
                tile_bp,
                store_as: Some(dotplot_id.clone()),
            }) {
                self.op_status = error.message.clone();
                self.op_error_popup = Some(error.message);
                return;
            }
            if let Err(error) = guard.apply(Operation::RenderDotplotSvg {
                seq_id: created_seq_id.clone(),
                dotplot_id,
                path: output_path.display().to_string(),
                flex_track_id: None,
                display_density_threshold: Some(self.dotplot_ui.display_density_threshold),
                display_intensity_gain: Some(self.dotplot_ui.display_intensity_gain),
                overlay_x_axis_mode: DotplotOverlayXAxisMode::PercentLength,
                overlay_anchor_exon: None,
            }) {
                self.op_status = error.message.clone();
                self.op_error_popup = Some(error.message);
                return;
            }
        }
        let folder = output_paths
            .first()
            .and_then(|path| path.parent())
            .map(|path| path.display().to_string())
            .unwrap_or_else(|| "<unknown>".to_string());
        self.op_status = format!(
            "Exported {} RNA-read sequence dotplot SVG(s) to '{}'",
            output_paths.len(),
            folder
        );
    }

    pub(super) fn export_highlighted_rna_read_sequence_dotplot_svg(
        &mut self,
        view: &SplicingExpertView,
        progress: &RnaReadInterpretProgress,
    ) {
        let highlighted_record_index = self.rna_seed_highlight_record_index.or_else(|| {
            self.selected_rna_top_hit_preview(progress)
                .map(|row| row.record_index)
        });
        let Some(highlighted_record_index) = highlighted_record_index else {
            self.op_status = "Highlight one RNA-read row before exporting a dotplot".to_string();
            return;
        };
        let Some(report) = self.current_saved_rna_read_report() else {
            self.op_status = "Load/save a Report ID before exporting RNA-read dotplots".to_string();
            return;
        };
        let Some(hit) = report
            .hits
            .iter()
            .find(|hit| hit.record_index == highlighted_record_index)
        else {
            self.op_status = format!(
                "Highlighted read #{} is not present in the saved report",
                highlighted_record_index + 1
            );
            return;
        };
        let reference_span_bp = view
            .region_end_1based
            .saturating_sub(view.region_start_1based.saturating_sub(1))
            .max(1);
        let requested_word_size =
            Self::parse_positive_usize_text(&self.dotplot_ui.word_size, "dotplot word_size")
                .unwrap_or(7);
        let requested_step_bp =
            Self::parse_positive_usize_text(&self.dotplot_ui.step_bp, "dotplot step_bp")
                .unwrap_or(1);
        let max_mismatches = Self::parse_optional_usize_text(
            &self.dotplot_ui.max_mismatches,
            "dotplot max_mismatches",
        )
        .ok()
        .flatten()
        .unwrap_or(0);
        let tile_bp = Self::parse_optional_usize_text(&self.dotplot_ui.tile_bp, "dotplot tile_bp")
            .ok()
            .flatten()
            .filter(|value| *value > 0);
        let (word_size, step_bp) = Self::resolve_rna_read_sequence_dotplot_word_and_step(
            hit.sequence.len(),
            reference_span_bp,
            requested_word_size,
            requested_step_bp,
        );
        let default_name = Self::default_rna_read_sequence_dotplot_svg_file_name(
            &report.report_id,
            hit,
            word_size,
            step_bp,
            max_mismatches,
            tile_bp,
        );
        let Some(path) = rfd::FileDialog::new()
            .set_file_name(&default_name)
            .add_filter("SVG", &["svg"])
            .save_file()
        else {
            return;
        };
        self.export_rna_read_sequence_dotplot_svgs(view, vec![hit], vec![path]);
    }

    pub(super) fn export_selected_rna_read_sequence_dotplot_svgs(
        &mut self,
        view: &SplicingExpertView,
    ) {
        let Some(report) = self.current_saved_rna_read_report() else {
            self.op_status = "Load/save a Report ID before exporting RNA-read dotplots".to_string();
            return;
        };
        let hits = Self::selected_rna_report_hits(&report, &self.rna_seed_selected_record_indices);
        if hits.is_empty() {
            self.op_status =
                "Select at least one RNA-read row before exporting dotplots".to_string();
            return;
        }
        let Some(folder) = rfd::FileDialog::new().pick_folder() else {
            return;
        };
        let reference_span_bp = view
            .region_end_1based
            .saturating_sub(view.region_start_1based.saturating_sub(1))
            .max(1);
        let requested_word_size =
            Self::parse_positive_usize_text(&self.dotplot_ui.word_size, "dotplot word_size")
                .unwrap_or(7);
        let requested_step_bp =
            Self::parse_positive_usize_text(&self.dotplot_ui.step_bp, "dotplot step_bp")
                .unwrap_or(1);
        let max_mismatches = Self::parse_optional_usize_text(
            &self.dotplot_ui.max_mismatches,
            "dotplot max_mismatches",
        )
        .ok()
        .flatten()
        .unwrap_or(0);
        let tile_bp = Self::parse_optional_usize_text(&self.dotplot_ui.tile_bp, "dotplot tile_bp")
            .ok()
            .flatten()
            .filter(|value| *value > 0);
        let output_paths = hits
            .iter()
            .map(|hit| {
                let (word_size, step_bp) = Self::resolve_rna_read_sequence_dotplot_word_and_step(
                    hit.sequence.len(),
                    reference_span_bp,
                    requested_word_size,
                    requested_step_bp,
                );
                folder.join(Self::default_rna_read_sequence_dotplot_svg_file_name(
                    &report.report_id,
                    hit,
                    word_size,
                    step_bp,
                    max_mismatches,
                    tile_bp,
                ))
            })
            .collect::<Vec<_>>();
        self.export_rna_read_sequence_dotplot_svgs(view, hits, output_paths);
    }

    pub(super) fn rna_top_hit_alignment_summary(row: &RnaReadTopHitPreview) -> String {
        if !row.aligned {
            return "none".to_string();
        }
        let tx = if row.best_alignment_transcript_id.is_empty() {
            "unknown".to_string()
        } else if row.best_alignment_transcript_label.is_empty()
            || row.best_alignment_transcript_label == row.best_alignment_transcript_id
        {
            row.best_alignment_transcript_id.clone()
        } else {
            format!(
                "{}/{}",
                row.best_alignment_transcript_id, row.best_alignment_transcript_label
            )
        };
        format!(
            "{} tx={} strand={} t={}-{} id={:.3} cov={:.3} score={} sec={}",
            if row.best_alignment_mode.is_empty() {
                "aligned"
            } else {
                row.best_alignment_mode.as_str()
            },
            tx,
            if row.best_alignment_strand.is_empty() {
                "na"
            } else {
                row.best_alignment_strand.as_str()
            },
            row.best_alignment_target_start_1based,
            row.best_alignment_target_end_1based,
            row.best_alignment_identity_fraction,
            row.best_alignment_query_coverage_fraction,
            row.best_alignment_score,
            row.secondary_mapping_count
        )
    }

    pub(super) fn format_rna_top_hit_preview_fasta_entry(row: &RnaReadTopHitPreview) -> String {
        let header_id = row.header_id.trim().replace(['\n', '\r', '\t'], " ");
        let gap_median = if row.seed_transcript_gap_count == 0 {
            "na".to_string()
        } else {
            format!("{:.2}", row.seed_median_transcript_gap)
        };
        let alignment_summary = Self::rna_top_hit_alignment_summary(row).replace(' ', "_");
        format!(
            ">{header_id} record_index={} score={:.3} wscore={:.4} wsupport={:.2} gap_med={} gap_n={} chain={:.2}/{} tx={} class={} origin_conf={:.3} strand_conf={:.3} strand={} opp={} ambig={} matched/tested={}/{} pass={} rc={} msa_eligible={} msa_reason={} align={} len={}\n{}",
            row.record_index + 1,
            row.seed_hit_fraction,
            row.weighted_seed_hit_fraction,
            row.weighted_matched_kmers,
            gap_median,
            row.seed_transcript_gap_count,
            row.seed_chain_support_fraction,
            row.seed_chain_support_kmers,
            if row.seed_chain_transcript_id.is_empty() {
                "none"
            } else {
                row.seed_chain_transcript_id.as_str()
            },
            row.origin_class.as_str(),
            row.origin_confidence,
            row.strand_confidence,
            if row.selected_strand.is_empty() {
                "na"
            } else {
                row.selected_strand.as_str()
            },
            row.competing_opposite_strand,
            row.ambiguous_strand_tie,
            row.matched_kmers,
            row.tested_kmers,
            row.passed_seed_filter,
            row.reverse_complement_applied,
            row.msa_eligible,
            row.msa_eligibility_reason.trim(),
            alignment_summary,
            row.read_length_bp,
            row.sequence,
        )
    }

    pub(super) fn copy_rna_top_hit_previews_as_fasta(
        &mut self,
        ui: &egui::Ui,
        rows: &[RnaReadTopHitPreview],
        source_label: &str,
    ) {
        if rows.is_empty() {
            self.op_status = format!("No reads available to copy from {source_label}");
            return;
        }
        let fasta = rows
            .iter()
            .map(Self::format_rna_top_hit_preview_fasta_entry)
            .collect::<Vec<_>>()
            .join("\n");
        ui.ctx().copy_text(fasta);
        self.op_status = format!(
            "Copied {} read sequence(s) as FASTA from {source_label}",
            rows.len()
        );
    }

    pub(super) fn run_splicing_rna_read_interpretation(&mut self, view: &SplicingExpertView) {
        if self.rna_read_task.is_some() {
            self.op_status = "RNA-read task is already running".to_string();
            return;
        }
        let op = match self.build_splicing_rna_read_interpret_operation(view) {
            Ok(op) => op,
            Err(message) => {
                self.op_status = message.clone();
                self.op_error_popup = Some(message);
                return;
            }
        };
        if let Operation::InterpretRnaReads {
            seq_id,
            scope,
            seed_filter,
            ..
        } = &op
        {
            self.save_engine_ops_state();
            self.refresh_rna_seed_catalog_preview(
                seq_id,
                view.target_feature_id,
                *scope,
                seed_filter.kmer_len,
            );
        }
        self.start_rna_read_interpretation(op);
    }

    pub(super) fn run_splicing_rna_read_alignment_phase(&mut self, view: &SplicingExpertView) {
        self.run_splicing_rna_read_alignment_phase_with_records(view, None);
    }

    pub(super) fn run_splicing_rna_read_alignment_phase_for_selected(
        &mut self,
        view: &SplicingExpertView,
        record_indices: Vec<usize>,
    ) {
        self.run_splicing_rna_read_alignment_phase_with_records(view, Some(record_indices));
    }

    pub(super) fn run_splicing_rna_read_alignment_phase_with_records(
        &mut self,
        view: &SplicingExpertView,
        record_indices: Option<Vec<usize>>,
    ) {
        if self.rna_read_task.is_some() {
            self.op_status = "RNA-read alignment phase is already running".to_string();
            return;
        }
        if let Some(indices) = record_indices.as_ref() {
            if indices.is_empty() {
                self.op_status =
                    "No top-hit rows selected; select at least one row first".to_string();
                return;
            }
        }
        let op = match self.build_splicing_rna_read_align_operation(view, record_indices) {
            Ok(op) => op,
            Err(message) => {
                self.op_status = message.clone();
                self.op_error_popup = Some(message);
                return;
            }
        };
        self.save_engine_ops_state();
        self.start_rna_read_interpretation(op);
    }

    pub(super) fn prepare_splicing_rna_read_interpretation_workflow(
        &mut self,
        view: &SplicingExpertView,
    ) {
        let workflow = match self.build_splicing_rna_read_workflow(view) {
            Ok(workflow) => workflow,
            Err(message) => {
                self.op_status = message.clone();
                self.op_error_popup = Some(message);
                return;
            }
        };
        self.workflow_run_id = workflow.run_id.clone();
        match serde_json::to_string_pretty(&workflow.ops) {
            Ok(ops_json) => {
                self.workflow_ops_json = ops_json;
                self.save_engine_ops_state();
                self.op_status = format!(
                    "Prepared RNA-read workflow '{}' in Engine Ops workflow runner (1 operation).",
                    self.workflow_run_id
                );
            }
            Err(e) => {
                let message = format!("Could not serialize RNA-read workflow ops JSON: {e}");
                self.op_status = message.clone();
                self.op_error_popup = Some(message);
            }
        }
    }

    pub(super) fn copy_splicing_rna_read_interpretation_workflow_json(
        &mut self,
        view: &SplicingExpertView,
        ctx: &egui::Context,
    ) {
        let workflow = match self.build_splicing_rna_read_workflow(view) {
            Ok(workflow) => workflow,
            Err(message) => {
                self.op_status = message.clone();
                self.op_error_popup = Some(message);
                return;
            }
        };
        match serde_json::to_string_pretty(&workflow) {
            Ok(text) => {
                ctx.copy_text(text);
                self.op_status = format!(
                    "Copied RNA-read workflow JSON for run_id='{}' (1 operation).",
                    workflow.run_id
                );
            }
            Err(e) => {
                let message = format!("Could not serialize RNA-read workflow JSON: {e}");
                self.op_status = message.clone();
                self.op_error_popup = Some(message);
            }
        }
    }

    pub(super) fn build_splicing_rna_read_workflow(
        &mut self,
        view: &SplicingExpertView,
    ) -> Result<Workflow, String> {
        let op = self.build_splicing_rna_read_interpret_operation(view)?;
        let (seq_id, report_id) = match &op {
            Operation::InterpretRnaReads {
                seq_id, report_id, ..
            } => (seq_id.as_str(), report_id.as_deref().unwrap_or("cdna")),
            _ => ("seq", "cdna"),
        };
        let run_id = if self.workflow_run_id.trim().is_empty() {
            Self::default_rna_read_workflow_run_id(seq_id, report_id)
        } else {
            self.workflow_run_id.trim().to_string()
        };
        Ok(Workflow {
            run_id,
            ops: vec![op],
        })
    }

    pub(super) fn build_splicing_rna_read_interpret_operation(
        &mut self,
        view: &SplicingExpertView,
    ) -> Result<Operation, String> {
        let seq_id = self
            .seq_id
            .clone()
            .ok_or_else(|| "No active sequence selected".to_string())?;
        let input_path = self.rna_reads_ui.input_path.trim().to_string();
        if input_path.is_empty() {
            return Err("RNA-read interpretation requires an input FASTA path".to_string());
        }
        let seed_filter = self.parse_rna_seed_filter_from_ui()?;
        let parsed_align_config = self.parse_rna_align_config_from_ui()?;
        let max_secondary_mappings = if matches!(
            self.rna_reads_ui.profile,
            RnaReadInterpretationProfile::NanoporeCdnaV1
        ) {
            0
        } else {
            parsed_align_config.max_secondary_mappings
        };
        let report_id = {
            let raw = self.rna_reads_ui.report_id.trim();
            if self.rna_reads_ui.report_id_auto_sync || raw.is_empty() {
                let default_id = Self::default_rna_read_report_id(view, &self.rna_reads_ui);
                self.rna_reads_ui.report_id = default_id.clone();
                Some(default_id)
            } else {
                Some(raw.to_string())
            }
        };
        if let Some(report_id) = report_id.as_deref() {
            self.validate_rna_read_report_id_for_splicing_view(report_id, view)?;
        }
        let checkpoint_every_reads = Self::parse_positive_usize_text(
            &self.rna_reads_ui.checkpoint_every_reads,
            "checkpoint_every_reads",
        )?;
        let checkpoint_path = {
            let trimmed = self.rna_reads_ui.checkpoint_path.trim();
            if trimmed.is_empty() {
                None
            } else {
                Some(trimmed.to_string())
            }
        };
        if self.rna_reads_ui.resume_from_checkpoint && checkpoint_path.is_none() {
            return Err(
                "resume_from_checkpoint=true requires a non-empty checkpoint_path".to_string(),
            );
        }
        let target_gene_ids = Self::parse_rna_target_gene_ids(&self.rna_reads_ui.target_gene_ids);
        Ok(Operation::InterpretRnaReads {
            seq_id,
            seed_feature_id: view.target_feature_id,
            profile: self.rna_reads_ui.profile,
            input_path,
            input_format: self.rna_reads_ui.input_format,
            scope: self.rna_reads_ui.scope,
            origin_mode: self.rna_reads_ui.origin_mode,
            target_gene_ids,
            roi_seed_capture_enabled: self.rna_reads_ui.roi_seed_capture_enabled,
            seed_filter,
            align_config: RnaReadAlignConfig {
                band_width_bp: parsed_align_config.band_width_bp,
                min_identity_fraction: parsed_align_config.min_identity_fraction,
                max_secondary_mappings,
            },
            report_id,
            report_mode: self.rna_reads_ui.report_mode,
            checkpoint_path,
            checkpoint_every_reads,
            resume_from_checkpoint: self.rna_reads_ui.resume_from_checkpoint,
        })
    }

    pub(super) fn parse_rna_seed_filter_from_ui(&self) -> Result<RnaReadSeedFilterConfig, String> {
        let kmer_len = Self::parse_positive_usize_text(&self.rna_reads_ui.kmer_len, "k-mer")?;
        if !(1..=16).contains(&kmer_len) {
            return Err("Invalid k-mer: expected within 1..=16".to_string());
        }
        let seed_stride_bp =
            Self::parse_positive_usize_text(&self.rna_reads_ui.seed_stride_bp, "seed_stride_bp")?;
        let min_seed_hit_fraction = match Self::parse_required_f64_text(
            &self.rna_reads_ui.min_seed_hit_fraction,
            "min_seed_hit_fraction",
        ) {
            Ok(value) if (0.0..=1.0).contains(&value) => value,
            Ok(_) => {
                return Err("Invalid min_seed_hit_fraction: expected within 0.0..=1.0".to_string());
            }
            Err(message) => return Err(message),
        };
        let min_weighted_seed_hit_fraction = match Self::parse_required_f64_text(
            &self.rna_reads_ui.min_weighted_seed_hit_fraction,
            "min_weighted_seed_hit_fraction",
        ) {
            Ok(value) if (0.0..=1.0).contains(&value) => value,
            Ok(_) => {
                return Err(
                    "Invalid min_weighted_seed_hit_fraction: expected within 0.0..=1.0".to_string(),
                );
            }
            Err(message) => return Err(message),
        };
        let min_unique_matched_kmers = Self::parse_required_usize_text(
            &self.rna_reads_ui.min_unique_matched_kmers,
            "min_unique_matched_kmers",
        )?;
        let max_median_transcript_gap = match Self::parse_required_f64_text(
            &self.rna_reads_ui.max_median_transcript_gap,
            "max_median_transcript_gap",
        ) {
            Ok(value) if value.is_finite() && value >= 1.0 => value,
            Ok(_) => {
                return Err(
                    "Invalid max_median_transcript_gap: expected finite value >= 1.0".to_string(),
                );
            }
            Err(message) => return Err(message),
        };
        let min_chain_consistency_fraction = match Self::parse_required_f64_text(
            &self.rna_reads_ui.min_chain_consistency_fraction,
            "min_chain_consistency_fraction",
        ) {
            Ok(value) if (0.0..=1.0).contains(&value) => value,
            Ok(_) => {
                return Err(
                    "Invalid min_chain_consistency_fraction: expected within 0.0..=1.0".to_string(),
                );
            }
            Err(message) => return Err(message),
        };
        let min_confirmed_exon_transitions = Self::parse_required_usize_text(
            &self.rna_reads_ui.min_confirmed_exon_transitions,
            "min_confirmed_exon_transitions",
        )?;
        let min_transition_support_fraction = match Self::parse_required_f64_text(
            &self.rna_reads_ui.min_transition_support_fraction,
            "min_transition_support_fraction",
        ) {
            Ok(value) if (0.0..=1.0).contains(&value) => value,
            Ok(_) => {
                return Err(
                    "Invalid min_transition_support_fraction: expected within 0.0..=1.0"
                        .to_string(),
                );
            }
            Err(message) => return Err(message),
        };
        let poly_t_prefix_min_bp = Self::parse_positive_usize_text(
            &self.rna_reads_ui.poly_t_prefix_min_bp,
            "poly_t_prefix_min_bp",
        )?;
        Ok(RnaReadSeedFilterConfig {
            kmer_len,
            seed_stride_bp,
            min_seed_hit_fraction,
            min_weighted_seed_hit_fraction,
            min_unique_matched_kmers,
            max_median_transcript_gap,
            min_chain_consistency_fraction,
            min_confirmed_exon_transitions,
            min_transition_support_fraction,
            cdna_poly_t_flip_enabled: self.rna_reads_ui.cdna_poly_t_flip_enabled,
            poly_t_prefix_min_bp,
        })
    }

    pub(super) fn parse_rna_align_config_from_ui(&self) -> Result<RnaReadAlignConfig, String> {
        let band_width_bp = Self::parse_positive_usize_text(
            &self.rna_reads_ui.align_band_width_bp,
            "align_band_width_bp",
        )?;
        let min_identity_fraction = match Self::parse_required_f64_text(
            &self.rna_reads_ui.align_min_identity_fraction,
            "align_min_identity_fraction",
        ) {
            Ok(value) if (0.0..=1.0).contains(&value) => value,
            Ok(_) => {
                return Err(
                    "Invalid align_min_identity_fraction: expected within 0.0..=1.0".to_string(),
                );
            }
            Err(message) => return Err(message),
        };
        let max_secondary_mappings = Self::parse_required_usize_text(
            &self.rna_reads_ui.align_max_secondary_mappings,
            "align_max_secondary_mappings",
        )?;
        Ok(RnaReadAlignConfig {
            band_width_bp,
            min_identity_fraction,
            max_secondary_mappings,
        })
    }

    pub(super) fn rna_alignment_debug_line(&self, progress: &RnaReadInterpretProgress) -> String {
        let requested_k = self
            .rna_reads_ui
            .kmer_len
            .trim()
            .parse::<usize>()
            .ok()
            .filter(|value| *value > 0)
            .unwrap_or(RnaReadSeedFilterConfig::default().kmer_len);
        let requested_w = self
            .rna_reads_ui
            .align_band_width_bp
            .trim()
            .parse::<usize>()
            .ok()
            .filter(|value| *value > 0)
            .unwrap_or(RnaReadAlignConfig::default().band_width_bp);
        let read_len = progress.median_read_length_bp.max(1);
        let effective_k = if read_len <= 2 {
            1
        } else {
            requested_k.clamp(3, read_len)
        };
        let active_note = if progress.align_compute_ms > 0.0 {
            "active"
        } else {
            "idle (seed-only phase or no aligned hits yet)"
        };
        format!(
            "ALIGN DEBUG: backend=pairwise::banded (dense fallback on empty band hit) | mode=semiglobal>local | k(req/eff~median-read)={requested_k}/{effective_k} | w={requested_w} | status={active_note}"
        )
    }

    pub(super) fn build_splicing_rna_read_align_operation(
        &mut self,
        view: &SplicingExpertView,
        selected_record_indices: Option<Vec<usize>>,
    ) -> Result<Operation, String> {
        let report_id = self.rna_reads_ui.report_id.trim();
        if report_id.is_empty() {
            return Err("Set a Report ID first before running alignment phase".to_string());
        }
        let report_id = report_id.to_string();
        self.validate_rna_read_report_id_for_splicing_view(&report_id, view)?;
        let align_config = self.parse_rna_align_config_from_ui()?;
        let selected_record_indices = selected_record_indices
            .unwrap_or_default()
            .into_iter()
            .collect::<BTreeSet<_>>()
            .into_iter()
            .collect::<Vec<_>>();
        Ok(Operation::AlignRnaReadReport {
            report_id,
            selection: self.rna_reads_ui.align_phase_selection,
            align_config_override: Some(align_config),
            selected_record_indices,
        })
    }

    pub(super) fn export_splicing_seed_hash_catalog(&mut self, view: &SplicingExpertView) {
        if self.rna_read_task.is_some() {
            self.op_status =
                "Seed-hash catalog export is disabled while RNA-read interpretation is running"
                    .to_string();
            return;
        }
        let Some(seq_id) = self.seq_id.clone() else {
            self.op_status = "No active sequence selected".to_string();
            return;
        };
        let Some(engine) = self.engine.clone() else {
            self.op_status = "No engine attached".to_string();
            return;
        };
        let kmer_len = match Self::parse_positive_usize_text(&self.rna_reads_ui.kmer_len, "k-mer") {
            Ok(value) if (1..=16).contains(&value) => value,
            Ok(_) => {
                let message = "Invalid k-mer: expected within 1..=16".to_string();
                self.op_status = message.clone();
                self.op_error_popup = Some(message);
                return;
            }
            Err(message) => {
                self.op_status = message.clone();
                self.op_error_popup = Some(message);
                return;
            }
        };
        let default_name = format!("{seq_id}_seed_hash_catalog_k{kmer_len}.tsv");
        let Some(path) = rfd::FileDialog::new()
            .set_file_name(&default_name)
            .add_filter("TSV", &["tsv", "txt"])
            .save_file()
        else {
            return;
        };
        let path_text = path.display().to_string();
        let mut seed_filter = RnaReadSeedFilterConfig::default();
        seed_filter.kmer_len = kmer_len;
        let export_result = match engine.read() {
            Ok(guard) => guard.export_rna_seed_hash_catalog(
                &seq_id,
                view.target_feature_id,
                self.rna_reads_ui.scope,
                &seed_filter,
                &path_text,
            ),
            Err(_) => Err(EngineError {
                code: ErrorCode::Internal,
                message: "Engine lock poisoned while exporting seed-hash catalog".to_string(),
            }),
        };
        match export_result {
            Ok((rows, unique_hashes)) => {
                self.op_status = format!(
                    "Exported RNA seed-hash catalog: rows={rows}, unique_hashes={unique_hashes}, path='{path_text}'"
                );
            }
            Err(err) => {
                self.op_status = err.message.clone();
                self.op_error_popup = Some(err.message);
            }
        }
    }

    pub(super) fn export_retained_rna_hits_fasta(&mut self) {
        let report_id = self.rna_reads_ui.report_id.trim().to_string();
        if report_id.is_empty() {
            let message = "Set a Report ID first to export retained reads".to_string();
            self.op_status = message.clone();
            self.op_error_popup = Some(message);
            return;
        }
        let default_name = format!("{report_id}_retained_reads.fa");
        let Some(path) = rfd::FileDialog::new()
            .set_file_name(&default_name)
            .add_filter("FASTA", &["fa", "fasta"])
            .save_file()
        else {
            return;
        };
        self.apply_operation_with_feedback(Operation::ExportRnaReadHitsFasta {
            report_id,
            path: path.display().to_string(),
            selection: RnaReadHitSelection::All,
            selected_record_indices: vec![],
            subset_spec: None,
        });
    }

    pub(super) fn export_rna_read_exon_paths_tsv(&mut self) {
        let report_id = self.rna_reads_ui.report_id.trim().to_string();
        if report_id.is_empty() {
            let message = "Set a Report ID first to export exon paths".to_string();
            self.op_status = message.clone();
            self.op_error_popup = Some(message);
            return;
        }
        let default_name = format!("{report_id}_exon_paths.tsv");
        let Some(path) = rfd::FileDialog::new()
            .set_file_name(&default_name)
            .add_filter("TSV", &["tsv", "txt"])
            .save_file()
        else {
            return;
        };
        self.apply_operation_with_feedback(Operation::ExportRnaReadExonPathsTsv {
            report_id,
            path: path.display().to_string(),
            selection: RnaReadHitSelection::All,
            selected_record_indices: vec![],
            subset_spec: None,
        });
    }

    pub(super) fn export_rna_read_exon_abundance_tsv(&mut self) {
        let report_id = self.rna_reads_ui.report_id.trim().to_string();
        if report_id.is_empty() {
            let message = "Set a Report ID first to export exon abundance".to_string();
            self.op_status = message.clone();
            self.op_error_popup = Some(message);
            return;
        }
        let default_name = format!("{report_id}_exon_abundance.tsv");
        let Some(path) = rfd::FileDialog::new()
            .set_file_name(&default_name)
            .add_filter("TSV", &["tsv", "txt"])
            .save_file()
        else {
            return;
        };
        self.apply_operation_with_feedback(Operation::ExportRnaReadExonAbundanceTsv {
            report_id,
            path: path.display().to_string(),
            selection: RnaReadHitSelection::All,
            selected_record_indices: vec![],
            subset_spec: None,
        });
    }

    pub(super) fn export_rna_read_score_density_svg(&mut self) {
        let report_id = self.rna_reads_ui.report_id.trim().to_string();
        if report_id.is_empty() {
            let message = "Set a Report ID first to export score-density SVG".to_string();
            self.op_status = message.clone();
            self.op_error_popup = Some(message);
            return;
        }
        let variant_suffix = match self.rna_read_evidence_ui.score_density_variant {
            RnaReadScoreDensityVariant::AllScored => "all_scored",
            RnaReadScoreDensityVariant::CompositeSeedGate => "composite_gate",
            RnaReadScoreDensityVariant::RetainedReplayCurrentControls => "retained_replay",
        };
        let default_name = format!("{report_id}_score_density_{variant_suffix}.svg");
        let Some(path) = rfd::FileDialog::new()
            .set_file_name(&default_name)
            .add_filter("SVG", &["svg"])
            .save_file()
        else {
            return;
        };
        let scale = if self.rna_read_evidence_ui.score_density_use_log_scale {
            RnaReadScoreDensityScale::Log
        } else {
            RnaReadScoreDensityScale::Linear
        };
        let variant = self.rna_read_evidence_ui.score_density_variant;
        let seed_filter_override =
            match self.current_rna_read_score_density_seed_filter_override(variant) {
                Ok(override_config) => override_config,
                Err(message) => {
                    self.op_status = message.clone();
                    self.op_error_popup = Some(message);
                    return;
                }
            };
        self.apply_operation_with_feedback(Operation::ExportRnaReadScoreDensitySvg {
            report_id,
            path: path.display().to_string(),
            scale,
            variant,
            seed_filter_override,
        });
    }

    pub(super) fn export_rna_read_target_quality(
        &mut self,
        report_id: String,
        target_gene_ids: Vec<String>,
        complete_rule: RnaReadGeneSupportCompleteRule,
    ) {
        if report_id.trim().is_empty() {
            let message = "Set a Report ID first to export target-quality comparison".to_string();
            self.op_status = message.clone();
            self.op_error_popup = Some(message);
            return;
        }
        if target_gene_ids.is_empty() {
            let message =
                "No target gene/group is available yet for target-quality export".to_string();
            self.op_status = message.clone();
            self.op_error_popup = Some(message);
            return;
        }
        let gene_token = target_gene_ids
            .join("_")
            .chars()
            .map(|ch| {
                if ch.is_ascii_alphanumeric() || matches!(ch, '_' | '-' | '.') {
                    ch
                } else {
                    '_'
                }
            })
            .collect::<String>();
        let default_name = format!("{report_id}_{gene_token}_target_quality.svg");
        let Some(path) = rfd::FileDialog::new()
            .set_file_name(&default_name)
            .add_filter("SVG", &["svg"])
            .add_filter("JSON", &["json"])
            .save_file()
        else {
            return;
        };
        self.apply_operation_with_feedback(Operation::ExportRnaReadTargetQuality {
            report_id,
            path: path.display().to_string(),
            gene_ids: target_gene_ids,
            complete_rule,
        });
    }

    pub(super) fn selected_rna_record_indices(&self) -> Vec<usize> {
        self.rna_seed_selected_record_indices
            .iter()
            .copied()
            .collect::<Vec<_>>()
    }

    pub(super) fn rna_read_record_indices_all_selected<I>(&self, record_indices: I) -> bool
    where
        I: IntoIterator<Item = usize>,
    {
        let mut saw_any = false;
        for record_index in record_indices {
            saw_any = true;
            if !self
                .rna_seed_selected_record_indices
                .contains(&record_index)
            {
                return false;
            }
        }
        saw_any
    }

    pub(super) fn set_rna_read_record_indices_selected<I>(
        &mut self,
        record_indices: I,
        selected: bool,
    ) -> usize
    where
        I: IntoIterator<Item = usize>,
    {
        let mut affected = 0usize;
        for record_index in record_indices {
            affected += 1;
            if selected {
                self.rna_seed_selected_record_indices.insert(record_index);
            } else {
                self.rna_seed_selected_record_indices.remove(&record_index);
            }
        }
        affected
    }

    pub(super) fn export_selected_rna_read_subset(&mut self, kind: RnaReadSelectedExportKind) {
        let selected_record_indices = self.selected_rna_record_indices();
        self.export_rna_read_subset_with_record_indices(
            kind,
            selected_record_indices,
            None,
            kind.empty_selection_message(),
        );
    }

    pub(super) fn export_rna_read_subset_with_record_indices(
        &mut self,
        kind: RnaReadSelectedExportKind,
        selected_record_indices: Vec<usize>,
        subset_spec: Option<String>,
        empty_selection_message: &str,
    ) {
        let Some(report_id) = self.selected_rna_read_evidence_report_id() else {
            let message = kind.missing_report_message().to_string();
            self.op_status = message.clone();
            self.op_error_popup = Some(message);
            return;
        };
        if selected_record_indices.is_empty() {
            let message = empty_selection_message.to_string();
            self.op_status = message.clone();
            self.op_error_popup = Some(message);
            return;
        }
        let default_name = kind.default_file_name(&report_id);
        let Some(path) = kind.configure_dialog(&default_name).save_file() else {
            return;
        };
        self.apply_operation_with_feedback(kind.build_operation(
            report_id,
            path.display().to_string(),
            selected_record_indices,
            subset_spec,
        ));
    }

    pub(super) fn start_rna_read_interpretation(&mut self, op: Operation) {
        let Some(engine) = self.engine.clone() else {
            self.op_status = "No engine attached".to_string();
            return;
        };
        if self.rna_read_task.is_some() {
            self.op_status = "RNA-read task is already running".to_string();
            return;
        }
        let input_path = match &op {
            Operation::InterpretRnaReads { input_path, .. } => input_path.clone(),
            Operation::AlignRnaReadReport { report_id, .. } => format!("report:{report_id}"),
            _ => "<unknown>".to_string(),
        };
        let (task_seq_id, task_seed_feature_id, task_report_id_hint) = match &op {
            Operation::InterpretRnaReads {
                seq_id,
                seed_feature_id,
                report_id,
                ..
            } => (seq_id.clone(), *seed_feature_id, report_id.clone()),
            Operation::AlignRnaReadReport { report_id, .. } => (
                self.rna_read_mapping_window_view
                    .as_ref()
                    .map(|view| view.seq_id.clone())
                    .unwrap_or_default(),
                self.rna_read_mapping_window_view
                    .as_ref()
                    .map(|view| view.target_feature_id)
                    .unwrap_or_default(),
                Some(report_id.clone()),
            ),
            _ => (String::new(), 0, None),
        };
        let operation_label = match &op {
            Operation::InterpretRnaReads { .. } => "Nanopore interpretation".to_string(),
            Operation::AlignRnaReadReport { .. } => "Nanopore alignment phase".to_string(),
            _ => "RNA-read task".to_string(),
        };
        let started = Instant::now();
        let (tx, rx) = mpsc::channel::<RnaReadTaskMessage>();
        self.rna_read_progress = None;
        self.clear_rna_read_report_scoped_selection_state();
        self.rna_stream_eta_text = None;
        self.rna_stream_eta_reads_processed = 0;
        if let Some(report_id) = task_report_id_hint.as_deref()
            && !report_id.trim().is_empty()
        {
            self.set_selected_rna_read_evidence_report_id(report_id.to_string());
        }
        self.op_status = match &op {
            Operation::InterpretRnaReads {
                seed_filter,
                origin_mode,
                target_gene_ids,
                roi_seed_capture_enabled,
                report_mode,
                checkpoint_path,
                checkpoint_every_reads,
                resume_from_checkpoint,
                ..
            } => format!(
                "{} started ({} mode, origin={}, targets={}, roi_capture={}, report_mode={}, checkpoint='{}', every={}, resume={}): '{}'",
                operation_label,
                if seed_filter.cdna_poly_t_flip_enabled {
                    "cDNA"
                } else {
                    "direct RNA"
                },
                origin_mode.as_str(),
                target_gene_ids.len(),
                roi_seed_capture_enabled,
                report_mode.as_str(),
                checkpoint_path
                    .as_deref()
                    .filter(|path| !path.trim().is_empty())
                    .unwrap_or("none"),
                checkpoint_every_reads,
                resume_from_checkpoint,
                input_path
            ),
            Operation::AlignRnaReadReport {
                report_id,
                selection,
                align_config_override,
                selected_record_indices,
            } => {
                let align_note = align_config_override
                    .as_ref()
                    .map(|cfg| {
                        format!(
                            "band={} min_identity={:.2} max_secondary={}",
                            cfg.band_width_bp,
                            cfg.min_identity_fraction,
                            cfg.max_secondary_mappings
                        )
                    })
                    .unwrap_or_else(|| "report-default".to_string());
                format!(
                    "{} started for report '{}' (selection={}{}align={}): '{}'",
                    operation_label,
                    report_id,
                    selection.as_str(),
                    if selected_record_indices.is_empty() {
                        ", ".to_string()
                    } else {
                        format!(
                            ", selected_record_indices={}, ",
                            selected_record_indices.len()
                        )
                    },
                    align_note,
                    input_path
                )
            }
            _ => format!("{} started: '{}'", operation_label, input_path),
        };
        self.rna_read_task = Some(RnaReadTask {
            started,
            seq_id: task_seq_id,
            seed_feature_id: task_seed_feature_id,
            report_id_hint: task_report_id_hint,
            input_path: input_path.clone(),
            operation_label: operation_label.clone(),
            cancel_requested: Arc::new(AtomicBool::new(false)),
            receiver: Arc::new(Mutex::new(rx)),
        });
        let cancel_requested = self
            .rna_read_task
            .as_ref()
            .map(|task| task.cancel_requested.clone())
            .expect("rna-read task just set");
        match op {
            Operation::InterpretRnaReads {
                seq_id,
                seed_feature_id,
                profile,
                input_path,
                input_format,
                scope,
                origin_mode,
                target_gene_ids,
                roi_seed_capture_enabled,
                seed_filter,
                align_config,
                report_id,
                report_mode,
                checkpoint_path,
                checkpoint_every_reads,
                resume_from_checkpoint,
            } => {
                std::thread::spawn(move || {
                    let tx_progress = tx.clone();
                    let cancel_for_progress = Arc::clone(&cancel_requested);
                    let cancel_for_compute = Arc::clone(&cancel_requested);
                    let outcome = match engine.read() {
                        Ok(guard) => {
                            let mut should_continue =
                                move || !cancel_for_compute.load(AtomicOrdering::Relaxed);
                            guard
                                .compute_rna_read_report_with_runtime_options_and_progress_and_cancel(
                                    &seq_id,
                                    seed_feature_id,
                                    profile,
                                    &input_path,
                                    input_format,
                                    scope,
                                    origin_mode,
                                    &target_gene_ids,
                                    roi_seed_capture_enabled,
                                    &seed_filter,
                                    &align_config,
                                    report_id.as_deref(),
                                    report_mode,
                                    checkpoint_path.as_deref(),
                                    checkpoint_every_reads,
                                    resume_from_checkpoint,
                                    &mut move |progress| {
                                        if cancel_for_progress
                                            .load(AtomicOrdering::Relaxed)
                                        {
                                            return false;
                                        }
                                        let OperationProgress::RnaReadInterpret(p) = progress
                                        else {
                                            return true;
                                        };
                                        let _ =
                                            tx_progress.send(RnaReadTaskMessage::Progress(p));
                                        !cancel_for_progress.load(AtomicOrdering::Relaxed)
                                    },
                                    &mut should_continue,
                                )
                                .map(RnaReadTaskOutcome::Interpret)
                        }
                        Err(_) => Err(EngineError {
                            code: ErrorCode::Internal,
                            message: "Engine lock poisoned while preparing RNA-read interpretation"
                                .to_string(),
                        }),
                    };
                    let _ = tx.send(RnaReadTaskMessage::Done(outcome));
                });
            }
            Operation::AlignRnaReadReport {
                report_id,
                selection,
                align_config_override,
                selected_record_indices,
            } => {
                std::thread::spawn(move || {
                    let tx_progress = tx.clone();
                    let cancel_for_progress = Arc::clone(&cancel_requested);
                    let cancel_for_compute = Arc::clone(&cancel_requested);
                    let outcome = match engine.read() {
                        Ok(guard) => {
                            let mut should_continue =
                                move || !cancel_for_compute.load(AtomicOrdering::Relaxed);
                            guard
                                .align_rna_read_report_with_progress_and_cancel(
                                    &report_id,
                                    selection,
                                    align_config_override.clone(),
                                    &selected_record_indices,
                                    &mut move |progress| {
                                        if cancel_for_progress.load(AtomicOrdering::Relaxed) {
                                            return false;
                                        }
                                        let OperationProgress::RnaReadInterpret(p) = progress
                                        else {
                                            return true;
                                        };
                                        let _ = tx_progress.send(RnaReadTaskMessage::Progress(p));
                                        !cancel_for_progress.load(AtomicOrdering::Relaxed)
                                    },
                                    &mut should_continue,
                                )
                                .map(|report| RnaReadTaskOutcome::Align {
                                    report,
                                    selection,
                                    align_config_override,
                                    selected_record_indices,
                                })
                        }
                        Err(_) => Err(EngineError {
                            code: ErrorCode::Internal,
                            message: "Engine lock poisoned while preparing RNA-read alignment"
                                .to_string(),
                        }),
                    };
                    let _ = tx.send(RnaReadTaskMessage::Done(outcome));
                });
            }
            other => {
                self.op_status = format!("Unsupported RNA-read worker operation: {other:?}");
                self.rna_read_task = None;
            }
        }
    }

    pub(super) fn async_task_repaint_delay(
        processed_progress_messages: usize,
        hit_progress_cap: bool,
    ) -> Duration {
        if hit_progress_cap {
            Duration::from_millis(20)
        } else if processed_progress_messages > 0 {
            Duration::from_millis(120)
        } else {
            Duration::from_millis(700)
        }
    }

    pub(super) fn commit_completed_rna_read_task_outcome(
        &mut self,
        outcome: RnaReadTaskOutcome,
    ) -> Result<OpResult, EngineError> {
        let result = {
            let Some(engine) = self.engine.as_ref() else {
                return Err(EngineError {
                    code: ErrorCode::Internal,
                    message: "No engine attached while finalizing RNA-read task".to_string(),
                });
            };
            let mut guard = engine.write().map_err(|_| EngineError {
                code: ErrorCode::Internal,
                message: "Engine lock poisoned while finalizing RNA-read task".to_string(),
            })?;
            match outcome {
                RnaReadTaskOutcome::Interpret(report) => guard.commit_rna_read_report(report),
                RnaReadTaskOutcome::Align {
                    report,
                    selection,
                    align_config_override,
                    selected_record_indices,
                } => guard.commit_aligned_rna_read_report(
                    report,
                    selection,
                    align_config_override,
                    selected_record_indices,
                ),
            }
        };
        if result.is_ok() {
            self.invalidate_rna_read_report_display_cache();
        }
        result
    }

    pub(super) fn poll_rna_read_task(&mut self, ctx: &egui::Context) {
        if self.rna_read_task.is_none() {
            return;
        }
        let mut done: Option<Result<RnaReadTaskOutcome, EngineError>> = None;
        let mut processed_progress_messages = 0usize;
        let mut hit_progress_cap = false;
        let mut task_still_running = false;
        if let Some(task) = &self.rna_read_task {
            match task.receiver.lock() {
                Ok(rx) => {
                    const MAX_PROGRESS_MESSAGES_PER_TICK: usize = 512;
                    let mut processed_progress = 0usize;
                    loop {
                        match rx.try_recv() {
                            Ok(RnaReadTaskMessage::Progress(progress)) => {
                                if let Some(selected_index) = self.rna_seed_highlight_record_index {
                                    let still_visible = progress
                                        .top_hits_preview
                                        .iter()
                                        .any(|row| row.record_index == selected_index);
                                    if !still_visible {
                                        self.rna_seed_highlight_record_index = None;
                                    }
                                }
                                if progress.reads_total == 0
                                    && progress.input_bytes_total > 0
                                    && progress.reads_processed
                                        != self.rna_stream_eta_reads_processed
                                {
                                    let elapsed_s = task.started.elapsed().as_secs_f64().max(0.001);
                                    let fraction = (progress.input_bytes_processed as f64
                                        / progress.input_bytes_total as f64)
                                        .clamp(0.0, 1.0);
                                    self.rna_stream_eta_text =
                                        if progress.input_bytes_processed > 0 && fraction > 0.0 {
                                            let estimated_total_s = elapsed_s / fraction;
                                            Some(format!(
                                                "ETA: {}",
                                                Self::format_duration_compact(
                                                    (estimated_total_s - elapsed_s).max(0.0)
                                                )
                                            ))
                                        } else {
                                            None
                                        };
                                    self.rna_stream_eta_reads_processed = progress.reads_processed;
                                }
                                self.rna_read_progress = Some(progress);
                                processed_progress = processed_progress.saturating_add(1);
                                if processed_progress >= MAX_PROGRESS_MESSAGES_PER_TICK {
                                    hit_progress_cap = true;
                                    break;
                                }
                            }
                            Ok(RnaReadTaskMessage::Done(res)) => {
                                done = Some(res);
                                break;
                            }
                            Err(TryRecvError::Disconnected) => {
                                done = Some(Err(EngineError {
                                    code: ErrorCode::Internal,
                                    message: "RNA-read worker disconnected unexpectedly"
                                        .to_string(),
                                }));
                                break;
                            }
                            Err(TryRecvError::Empty) => break,
                        }
                    }
                    processed_progress_messages = processed_progress;
                    task_still_running = done.is_none();
                }
                Err(_) => {
                    done = Some(Err(EngineError {
                        code: ErrorCode::Internal,
                        message: "RNA-read progress channel lock poisoned".to_string(),
                    }));
                }
            }
        }

        if done.is_none() && task_still_running && !self.show_rna_read_mapping_window {
            ctx.request_repaint_after(Self::async_task_repaint_delay(
                processed_progress_messages,
                hit_progress_cap,
            ));
        }

        if let Some(done) = done {
            let (started, cancel_requested, operation_label, route_status_to_workspace) = self
                .rna_read_task
                .as_ref()
                .map(|task| {
                    (
                        task.started,
                        task.cancel_requested.load(AtomicOrdering::Relaxed),
                        task.operation_label.clone(),
                        self.rna_read_mapping_workspace_matches_task(
                            &task.seq_id,
                            task.seed_feature_id,
                        ),
                    )
                })
                .unwrap_or_else(|| (Instant::now(), false, "RNA-read task".to_string(), false));
            let previous_op_status = self.op_status.clone();
            self.rna_read_task = None;
            match done {
                Ok(outcome) => match self.commit_completed_rna_read_task_outcome(outcome) {
                    Ok(result) => self.handle_operation_success(result, started),
                    Err(err) => self.handle_operation_error(err, started),
                },
                Err(err) => {
                    let cancelled =
                        cancel_requested && err.message.to_ascii_lowercase().contains("cancel");
                    if cancelled {
                        self.op_status = format!(
                            "{} cancelled after {:.1}s",
                            operation_label,
                            started.elapsed().as_secs_f32()
                        );
                    } else {
                        self.handle_operation_error(err, started);
                    }
                }
            }
            self.capture_rna_read_mapping_workspace_status(
                previous_op_status,
                route_status_to_workspace,
            );
        }
    }

    pub(super) fn render_rna_read_seed_histogram(
        &self,
        ui: &mut egui::Ui,
        progress: &RnaReadInterpretProgress,
        seed_catalog: &[RnaSeedHashCatalogEntry],
        template_audit: &[RnaSeedHashTemplateAuditEntry],
        view: &SplicingExpertView,
        show_exons: bool,
        show_introns: bool,
        exonic_coords: bool,
    ) {
        if progress.bins.is_empty() {
            ui.small("No seed-confirmation histogram bins available yet.");
            return;
        }
        let desired = Vec2::new(ui.available_width().max(280.0), 150.0);
        let (rect, response) = ui.allocate_exact_size(desired, egui::Sense::hover());
        let painter = ui.painter_at(rect);
        let border = egui::Stroke::new(1.0, egui::Color32::from_gray(120));
        painter.rect_stroke(rect, 4.0, border, egui::StrokeKind::Inside);

        let mid_y = rect.center().y;
        painter.line_segment(
            [
                egui::pos2(rect.left(), mid_y),
                egui::pos2(rect.right(), mid_y),
            ],
            egui::Stroke::new(1.0, egui::Color32::from_gray(120)),
        );

        let max_count = progress
            .bins
            .iter()
            .map(|bin| bin.confirmed_plus.max(bin.confirmed_minus))
            .max()
            .unwrap_or(0)
            .max(1);
        let max_count_sqrt = (max_count as f32).sqrt().max(1.0);
        let half_h = (rect.height() * 0.44).max(1.0);
        let genomic_left_label = progress
            .bins
            .first()
            .map(|bin| bin.start_1based)
            .unwrap_or(1);
        let genomic_right_label = progress
            .bins
            .last()
            .map(|bin| bin.end_1based)
            .unwrap_or(genomic_left_label);
        let merged_exons = Self::merged_splicing_exon_ranges(view);
        let exonic_total_len = merged_exons.iter().fold(0usize, |acc, (start, end)| {
            acc.saturating_add(end.saturating_sub(*start).saturating_add(1))
        });
        let use_exonic_coords = exonic_coords && exonic_total_len > 0;
        let left_label = if use_exonic_coords {
            1
        } else {
            genomic_left_label
        };
        let right_label = if use_exonic_coords {
            exonic_total_len.max(1)
        } else {
            genomic_right_label.max(genomic_left_label)
        };
        let left_f = left_label as f32;
        let right_f = right_label as f32;
        let span = (right_f - left_f).max(1.0);
        let draw_strand_bar =
            |x0: f32, x1: f32, count: u64, strand_minus: bool, painter: &egui::Painter| {
                if count == 0 || x1 <= x0 {
                    return;
                }
                let h = (((count as f32).sqrt() / max_count_sqrt) * half_h).max(1.0);
                if strand_minus {
                    let y0 = mid_y + 0.6;
                    let mut y1 = mid_y + h;
                    if y1 <= y0 {
                        y1 = y0 + 1.0;
                    }
                    painter.rect_filled(
                        egui::Rect::from_min_max(egui::pos2(x0, y0), egui::pos2(x1, y1)),
                        0.0,
                        egui::Color32::from_rgb(249, 115, 22),
                    );
                    if (x1 - x0) >= 22.0 && h >= 12.0 {
                        painter.text(
                            egui::pos2((x0 + x1) * 0.5, y1 - 1.0),
                            egui::Align2::CENTER_BOTTOM,
                            Self::format_count_compact_km(count),
                            egui::FontId::monospace(8.0),
                            egui::Color32::from_rgb(12, 12, 12),
                        );
                    }
                } else {
                    let y0 = mid_y - h;
                    let mut y1 = mid_y - 0.6;
                    if y1 <= y0 {
                        y1 = y0 + 1.0;
                    }
                    painter.rect_filled(
                        egui::Rect::from_min_max(egui::pos2(x0, y0), egui::pos2(x1, y1)),
                        0.0,
                        egui::Color32::from_rgb(37, 99, 235),
                    );
                    if (x1 - x0) >= 22.0 && h >= 12.0 {
                        painter.text(
                            egui::pos2((x0 + x1) * 0.5, y0 + 1.0),
                            egui::Align2::CENTER_TOP,
                            Self::format_count_compact_km(count),
                            egui::FontId::monospace(8.0),
                            egui::Color32::from_rgb(245, 245, 245),
                        );
                    }
                }
            };
        if use_exonic_coords {
            let merged_exons_for_bars = Self::merged_splicing_exon_ranges(view);
            let exonic_total = merged_exons_for_bars
                .iter()
                .fold(0usize, |acc, (start, end)| {
                    acc.saturating_add(end.saturating_sub(*start).saturating_add(1))
                });
            let left_f_exo = 1.0f32;
            let right_f_exo = exonic_total.max(1) as f32;
            let span_exo = (right_f_exo - left_f_exo).max(1.0);
            for bin in &progress.bins {
                let spans = Self::project_genomic_interval_to_exonic(
                    &merged_exons_for_bars,
                    bin.start_1based,
                    bin.end_1based,
                );
                for (span_start, span_end) in spans {
                    let x0 = rect.left()
                        + ((span_start as f32 - left_f_exo) / span_exo).clamp(0.0, 1.0)
                            * rect.width()
                        + 0.8;
                    let x1 = rect.left()
                        + ((span_end as f32 - left_f_exo) / span_exo).clamp(0.0, 1.0)
                            * rect.width()
                        - 0.8;
                    draw_strand_bar(x0, x1, bin.confirmed_plus, false, &painter);
                    draw_strand_bar(x0, x1, bin.confirmed_minus, true, &painter);
                }
            }
        } else {
            let bin_w = rect.width() / progress.bins.len().max(1) as f32;
            for (idx, bin) in progress.bins.iter().enumerate() {
                let x0 = rect.left() + idx as f32 * bin_w + 0.8;
                let x1 = if idx + 1 == progress.bins.len() {
                    rect.right() - 0.8
                } else {
                    rect.left() + (idx + 1) as f32 * bin_w - 0.8
                };
                draw_strand_bar(x0, x1, bin.confirmed_plus, false, &painter);
                draw_strand_bar(x0, x1, bin.confirmed_minus, true, &painter);
            }
        }
        let sorted_exons = merged_exons.clone();
        if show_introns {
            if use_exonic_coords {
                let mut running = 0usize;
                for (idx, (start, end)) in sorted_exons.iter().enumerate() {
                    running = running.saturating_add(end.saturating_sub(*start).saturating_add(1));
                    if idx + 1 >= sorted_exons.len() {
                        continue;
                    }
                    let t = ((running as f32 - left_f) / span).clamp(0.0, 1.0);
                    let x = rect.left() + t * rect.width();
                    painter.line_segment(
                        [
                            egui::pos2(x, rect.top() + 2.0),
                            egui::pos2(x, rect.top() + 11.0),
                        ],
                        egui::Stroke::new(1.0, egui::Color32::from_gray(140)),
                    );
                }
            } else {
                for pair in sorted_exons.windows(2) {
                    let left = pair[0];
                    let right = pair[1];
                    let intron_start = left.1.saturating_add(1);
                    if right.0 <= intron_start {
                        continue;
                    }
                    let intron_end = right.0.saturating_sub(1);
                    if intron_end < genomic_left_label || intron_start > genomic_right_label {
                        continue;
                    }
                    let start_t = ((intron_start as f32 - left_f) / span).clamp(0.0, 1.0);
                    let end_t = ((intron_end as f32 - left_f) / span).clamp(0.0, 1.0);
                    let x0 = rect.left() + start_t * rect.width();
                    let x1 = rect.left() + end_t * rect.width();
                    painter.line_segment(
                        [
                            egui::pos2(x0, rect.top() + 7.0),
                            egui::pos2(x1, rect.top() + 7.0),
                        ],
                        egui::Stroke::new(2.0, egui::Color32::from_gray(140)),
                    );
                }
            }
        }
        if show_exons {
            if use_exonic_coords {
                let exon_rect = egui::Rect::from_min_max(
                    egui::pos2(rect.left(), rect.top() + 3.0),
                    egui::pos2(rect.right(), rect.top() + 10.0),
                );
                painter.rect_filled(exon_rect, 2.0, egui::Color32::from_rgb(34, 197, 94));
            } else {
                for exon in &sorted_exons {
                    if exon.1 < genomic_left_label || exon.0 > genomic_right_label {
                        continue;
                    }
                    let start_t = ((exon.0 as f32 - left_f) / span).clamp(0.0, 1.0);
                    let end_t = ((exon.1 as f32 - left_f) / span).clamp(0.0, 1.0);
                    let x0 = rect.left() + start_t * rect.width();
                    let mut x1 = rect.left() + end_t * rect.width();
                    if x1 <= x0 {
                        x1 = x0 + 1.0;
                    }
                    let exon_rect = egui::Rect::from_min_max(
                        egui::pos2(x0, rect.top() + 3.0),
                        egui::pos2(x1, rect.top() + 10.0),
                    );
                    painter.rect_filled(exon_rect, 2.0, egui::Color32::from_rgb(34, 197, 94));
                }
            }
        }
        painter.text(
            egui::pos2(rect.left() + 4.0, rect.bottom() - 14.0),
            egui::Align2::LEFT_BOTTOM,
            format!("{left_label}"),
            egui::TextStyle::Small.resolve(ui.style()),
            egui::Color32::from_gray(120),
        );
        painter.text(
            egui::pos2(rect.right() - 4.0, rect.bottom() - 14.0),
            egui::Align2::RIGHT_BOTTOM,
            format!("{right_label}"),
            egui::TextStyle::Small.resolve(ui.style()),
            egui::Color32::from_gray(120),
        );
        painter.text(
            egui::pos2(rect.right() - 4.0, rect.top() + 2.0),
            egui::Align2::RIGHT_TOP,
            if use_exonic_coords {
                "coords: exonic-only"
            } else {
                "coords: genomic"
            },
            egui::FontId::monospace(10.0),
            egui::Color32::from_gray(110),
        );
        if !seed_catalog.is_empty() {
            let selected_top_hit = self.selected_rna_top_hit_preview(progress);
            let template_audit_by_feature_id = template_audit
                .iter()
                .map(|row| (row.transcript_feature_id, row))
                .collect::<HashMap<_, _>>();
            let selected_kmer_len = seed_catalog
                .first()
                .map(|row| row.kmer_sequence.len())
                .filter(|len| *len > 0)
                .unwrap_or_else(|| {
                    self.rna_reads_ui
                        .kmer_len
                        .trim()
                        .parse::<usize>()
                        .ok()
                        .filter(|len| *len > 0)
                        .unwrap_or(9)
                });
            let selected_seed_stride = self
                .rna_reads_ui
                .seed_stride_bp
                .trim()
                .parse::<usize>()
                .ok()
                .filter(|value| *value > 0)
                .unwrap_or_else(|| RnaReadSeedFilterConfig::default().seed_stride_bp);
            let recompute_started = Instant::now();
            let selected_seed_counts = selected_top_hit
                .map(|row| {
                    Self::collect_read_seed_bit_counts(
                        &row.sequence,
                        selected_kmer_len,
                        selected_seed_stride,
                    )
                })
                .unwrap_or_default();
            let mut selected_seed_budget = selected_seed_counts.clone();
            let recompute_elapsed_ms = recompute_started.elapsed().as_secs_f64() * 1000.0;
            let mut seed_occurrence_counts = HashMap::<u32, usize>::new();
            for entry in seed_catalog {
                *seed_occurrence_counts.entry(entry.seed_bits).or_insert(0) += 1;
            }
            let dot_offset = 6.0;
            let dot_color = egui::Color32::from_rgb(220, 38, 38);
            let selected_color = egui::Color32::from_rgb(22, 163, 74);
            let mut hovered: Option<(f32, &RnaSeedHashCatalogEntry)> = None;
            let pointer = response.hover_pos();
            let mut plus_unique = HashSet::<usize>::new();
            let mut minus_unique = HashSet::<usize>::new();
            let mut rendered_pixel_buckets = HashSet::<(i32, bool)>::new();
            let mut selected_supported_positions = 0usize;
            let mut selected_supported_pixel_buckets = HashSet::<(i32, bool)>::new();
            for entry in seed_catalog {
                let maybe_display_pos = if use_exonic_coords {
                    Self::genomic_to_exonic_pos_1based(&merged_exons, entry.genomic_pos_1based)
                } else if entry.genomic_pos_1based >= genomic_left_label
                    && entry.genomic_pos_1based <= genomic_right_label
                {
                    Some(entry.genomic_pos_1based)
                } else {
                    None
                };
                let Some(display_pos_1based) = maybe_display_pos else {
                    continue;
                };
                let t = ((display_pos_1based as f32 - left_f) / span).clamp(0.0, 1.0);
                let x = rect.left() + t * rect.width();
                let strand_minus = entry.strand.trim() == "-";
                let y = if strand_minus {
                    mid_y + dot_offset
                } else {
                    mid_y - dot_offset
                };
                let pixel_x = ((x - rect.left()).round() as i32).clamp(0, rect.width() as i32);
                rendered_pixel_buckets.insert((pixel_x, strand_minus));
                if strand_minus {
                    minus_unique.insert(entry.genomic_pos_1based);
                } else {
                    plus_unique.insert(entry.genomic_pos_1based);
                }
                painter.circle_filled(egui::pos2(x, y), 1.6, dot_color);
                let mut highlight_selected = false;
                if let Some(remaining_budget) = selected_seed_budget.get_mut(&entry.seed_bits) {
                    if *remaining_budget > 0 {
                        *remaining_budget = remaining_budget.saturating_sub(1);
                        highlight_selected = true;
                    }
                }
                if highlight_selected {
                    selected_supported_positions = selected_supported_positions.saturating_add(1);
                    selected_supported_pixel_buckets.insert((pixel_x, strand_minus));
                    let spike_tip = if strand_minus { y + 5.0 } else { y - 5.0 };
                    painter.line_segment(
                        [egui::pos2(x, y), egui::pos2(x, spike_tip)],
                        egui::Stroke::new(1.2, selected_color),
                    );
                    painter.circle_filled(egui::pos2(x, y), 2.5, selected_color);
                }
                if let Some(pointer_pos) = pointer {
                    let dx = pointer_pos.x - x;
                    let dy = pointer_pos.y - y;
                    let dist_sq = dx * dx + dy * dy;
                    if dist_sq <= 81.0 {
                        match hovered {
                            Some((best_dist_sq, _)) if dist_sq >= best_dist_sq => {}
                            _ => hovered = Some((dist_sq, entry)),
                        }
                    }
                }
            }
            if let Some((_, entry)) = hovered {
                response.clone().on_hover_ui_at_pointer(|ui| {
                    ui.monospace(format!(
                        "seed={} ({:08X}) seq={}",
                        entry.seed_bits, entry.seed_bits, entry.kmer_sequence
                    ));
                    ui.monospace(format!(
                        "genomic={} strand={} template_offset={}",
                        entry.genomic_pos_1based, entry.strand, entry.template_offset_0based
                    ));
                    ui.label(format!(
                        "transcript n-{} {}",
                        entry.transcript_feature_id, entry.transcript_id
                    ));
                    if let Some(template) =
                        template_audit_by_feature_id.get(&entry.transcript_feature_id)
                    {
                        ui.monospace(format!(
                            "hash window={}..{} of {} bp template",
                            entry.template_offset_0based,
                            entry
                                .template_offset_0based
                                .saturating_add(entry.kmer_sequence.len()),
                            template.template_length_bp
                        ));
                        ui.monospace(format!(
                            "genomic path={}→{} rc_from_genome={}",
                            template.template_first_genomic_pos_1based,
                            template.template_last_genomic_pos_1based,
                            template.reverse_complemented_from_genome
                        ));
                    }
                });
            }
            ui.small(format!(
                "Seed hashes indexed: {} rows | unique genomic starts: {} (+{} / -{}) | rendered dot buckets at current width: {} | coord-mode={}",
                seed_catalog.len(),
                plus_unique.len() + minus_unique.len(),
                plus_unique.len(),
                minus_unique.len(),
                rendered_pixel_buckets.len(),
                if use_exonic_coords { "exonic" } else { "genomic" },
            ));
            let repeated_seed_bits = seed_occurrence_counts
                .values()
                .filter(|count| **count > 1)
                .count();
            let max_seed_occurrence = seed_occurrence_counts.values().copied().max().unwrap_or(0);
            ui.small(format!(
                "Seed bit diversity: {} unique bits | repeated bits={} | max occurrences per bit={}",
                seed_occurrence_counts.len(),
                repeated_seed_bits,
                max_seed_occurrence
            ));
            if !template_audit.is_empty() {
                ui.small(format!(
                    "Template audit: {} transcript-oriented template sequence(s) loaded; '-' strand templates are reverse-complemented exon-concatenated transcript sequences. Hover a seed dot for its template window and inspect the template-audit section below for the exact full sequence.",
                    template_audit.len()
                ));
            }
            if let Some(selected) = selected_top_hit {
                let alignment_summary = Self::rna_top_hit_alignment_summary(selected);
                ui.small(format!(
                    "Selected top read #{} supports {} hash positions ({} visible pixel buckets); hash recompute {:.2} ms | score={:.3} wscore={:.4} gap-med={} gap-n={} chain={:.2}/{} tx={} class={} oconf={:.2} sconf={:.2} align={}",
                    selected.record_index + 1,
                    selected_supported_positions,
                    selected_supported_pixel_buckets.len(),
                    recompute_elapsed_ms,
                    selected.seed_hit_fraction,
                    selected.weighted_seed_hit_fraction,
                    if selected.seed_transcript_gap_count == 0 {
                        "na".to_string()
                    } else {
                        format!("{:.2}", selected.seed_median_transcript_gap)
                    },
                    selected.seed_transcript_gap_count,
                    selected.seed_chain_support_fraction,
                    selected.seed_chain_support_kmers,
                    if selected.seed_chain_transcript_id.is_empty() {
                        "none"
                    } else {
                        selected.seed_chain_transcript_id.as_str()
                    },
                    selected.origin_class.as_str(),
                    selected.origin_confidence,
                    selected.strand_confidence,
                    alignment_summary,
                ));
            } else {
                ui.small("Select one of the best-performing reads to light supported hash positions in green.");
            }
        }
        ui.small(
            "Seed-confirmation frequency by genomic position (blue/orange bars, sqrt-scaled). Red dots mark indexed seed-hash starts (+ above baseline, - below). Optional exon/intron guides are shown at the top (green/gray). Bar labels use compact counts (k/M).",
        );
        ui.collapsing(
            format!(
                "Template sequences used for hashing ({} indexed transcript templates)",
                template_audit.len()
            ),
            |ui| {
                if template_audit.is_empty() {
                    ui.small("No transcript-template audit rows are loaded.");
                    return;
                }
                ui.small(
                    "These are the exact transcript-oriented exon-concatenated sequences used to derive the seed hashes. Minus-strand templates are reverse-complemented from the genomic exon path before hashing.",
                );
                egui::ScrollArea::vertical()
                    .max_height(200.0)
                    .show(ui, |ui| {
                        for entry in template_audit.iter().take(24) {
                            ui.monospace(format!(
                                "n-{}\t{}\tstrand={}\tlen={}\tgenomic_path={}→{}\trc_from_genome={}",
                                entry.transcript_feature_id,
                                if entry.transcript_label.is_empty() {
                                    entry.transcript_id.as_str()
                                } else {
                                    entry.transcript_label.as_str()
                                },
                                entry.strand,
                                entry.template_length_bp,
                                entry.template_first_genomic_pos_1based,
                                entry.template_last_genomic_pos_1based,
                                entry.reverse_complemented_from_genome,
                            ));
                            ui.add(
                                egui::Label::new(
                                    egui::RichText::new(&entry.template_sequence).monospace(),
                                )
                                .wrap(),
                            );
                            ui.add_space(6.0);
                        }
                        if template_audit.len() > 24 {
                            ui.small(format!(
                                "Showing first 24 of {} template sequences.",
                                template_audit.len()
                            ));
                        }
                    });
            },
        );
        ui.collapsing(
            format!(
                "Seed hash preview (first 120 of {} seed-passed reads)",
                progress.seed_passed
            ),
            |ui| {
            if seed_catalog.is_empty() {
                ui.small("No seed-hash catalog loaded for preview.");
                return;
            }
            ui.small(format!(
                "Loaded seed hashes: {} (showing first 120; seed-passed reads so far: {}; export TSV for full list).",
                seed_catalog.len(),
                progress.seed_passed
            ));
            egui::ScrollArea::vertical()
                .max_height(120.0)
                .show(ui, |ui| {
                    for entry in seed_catalog.iter().take(120) {
                        ui.monospace(format!(
                            "{}\t{}\t{}\t{}\t{}",
                            entry.genomic_pos_1based,
                            entry.strand,
                            entry.kmer_sequence,
                            entry.seed_bits,
                            entry.transcript_id,
                        ));
                    }
                });
        });
        ui.small(format!(
            "Composite seed-pass summary: {} reads passed the full current seed gate.",
            progress.seed_passed
        ));
    }

    pub(super) fn render_rna_read_score_density_plot(
        &mut self,
        ui: &mut egui::Ui,
        progress: &RnaReadInterpretProgress,
        use_log_scale: bool,
    ) {
        let score_density_variant = self.rna_read_evidence_ui.score_density_variant;
        let saved_report = self.current_saved_rna_read_report();
        let score_density_seed_filter_override = match self
            .current_rna_read_score_density_seed_filter_override(score_density_variant)
        {
            Ok(override_config) => override_config,
            Err(message) => {
                ui.small(format!(
                        "Current histogram population is unavailable until the seed controls are valid: {message}"
                    ));
                return;
            }
        };
        let score_density_bins = if let Some(report) = saved_report.as_ref() {
            GentleEngine::score_density_bins_for_report_with_override(
                report,
                score_density_variant,
                score_density_seed_filter_override.as_ref(),
            )
            .0
        } else {
            Self::rna_read_score_density_bins_for_progress(progress, score_density_variant)
        };
        if score_density_bins.is_empty() {
            ui.small(match score_density_variant {
                RnaReadScoreDensityVariant::RetainedReplayCurrentControls => {
                    "No retained-replay bins available yet. This mode needs a saved report and only replays over retained rows."
                }
                _ => "No score-density bins available yet.",
            });
            return;
        }
        let scale_label = if use_log_scale {
            "log-scale counts"
        } else {
            "linear counts"
        };
        ui.small(format!(
            "Seed-hit score density across {} ({scale_label})",
            Self::rna_read_score_density_variant_label(score_density_variant)
        ));
        ui.small(Self::rna_read_score_density_variant_note(
            score_density_variant,
        ));
        let desired = Vec2::new(ui.available_width().max(280.0), 120.0);
        let (rect, response) = ui.allocate_exact_size(desired, egui::Sense::click());
        let painter = ui.painter_at(rect);
        painter.rect_stroke(
            rect,
            4.0,
            egui::Stroke::new(1.0, egui::Color32::from_gray(120)),
            egui::StrokeKind::Inside,
        );
        let max_count = score_density_bins.iter().copied().max().unwrap_or(0);
        if max_count == 0 {
            painter.text(
                rect.center(),
                egui::Align2::CENTER_CENTER,
                "Awaiting scored reads",
                egui::FontId::monospace(11.0),
                egui::Color32::from_gray(120),
            );
            return;
        }
        let axis_bottom = rect.bottom() - 16.0;
        let plot_top = rect.top() + 4.0;
        let plot_h = (axis_bottom - plot_top).max(1.0);
        let max_scaled_count = if use_log_scale {
            (max_count as f32 + 1.0).ln().max(1.0)
        } else {
            max_count as f32
        };
        let bin_count = score_density_bins.len().max(1);
        let min_seed_hit_fraction = saved_report
            .as_ref()
            .map(|report| {
                score_density_seed_filter_override
                    .as_ref()
                    .map(|seed_filter| seed_filter.min_seed_hit_fraction)
                    .unwrap_or(report.seed_filter.min_seed_hit_fraction)
            })
            .or_else(|| {
                Self::parse_required_f64_text(
                    &self.rna_reads_ui.min_seed_hit_fraction,
                    "min_seed_hit_fraction",
                )
                .ok()
            })
            .unwrap_or(0.30)
            .clamp(0.0, 1.0);
        let selected_bin_index = self
            .rna_read_alignment_effect_score_bin_index
            .map(|idx| idx.min(bin_count.saturating_sub(1)));
        let bin_w = rect.width() / bin_count as f32;
        for (idx, count) in score_density_bins.iter().enumerate() {
            if *count == 0 {
                continue;
            }
            let x0 = rect.left() + idx as f32 * bin_w + 0.8;
            let x1 = if idx + 1 == score_density_bins.len() {
                rect.right() - 0.8
            } else {
                rect.left() + (idx + 1) as f32 * bin_w - 0.8
            };
            if x1 <= x0 {
                continue;
            }
            let scaled_count = if use_log_scale {
                (*count as f32 + 1.0).ln()
            } else {
                *count as f32
            };
            let h = (scaled_count / max_scaled_count) * plot_h;
            let bar = egui::Rect::from_min_max(
                egui::pos2(x0, axis_bottom - h),
                egui::pos2(x1, axis_bottom - 0.6),
            );
            let is_selected = selected_bin_index == Some(idx);
            let fill = if is_selected {
                egui::Color32::from_rgb(5, 150, 105)
            } else {
                egui::Color32::from_rgb(16, 185, 129)
            };
            painter.rect_filled(bar, 0.0, fill);
            if is_selected {
                painter.rect_stroke(
                    bar.expand(1.2),
                    0.0,
                    egui::Stroke::new(2.0, egui::Color32::from_rgb(8, 47, 73)),
                    egui::StrokeKind::Inside,
                );
            }
            if (x1 - x0) >= 22.0 {
                let label = Self::format_count_compact_km(*count);
                if h >= 12.0 {
                    painter.text(
                        egui::pos2((x0 + x1) * 0.5, (axis_bottom - h + 1.0).max(plot_top + 1.0)),
                        egui::Align2::CENTER_TOP,
                        label,
                        egui::FontId::monospace(8.0),
                        egui::Color32::from_rgb(7, 24, 17),
                    );
                } else {
                    painter.text(
                        egui::pos2((x0 + x1) * 0.5, (axis_bottom - h - 1.0).max(plot_top + 1.0)),
                        egui::Align2::CENTER_BOTTOM,
                        label,
                        egui::FontId::monospace(8.0),
                        egui::Color32::from_rgb(18, 18, 18),
                    );
                }
            }
        }
        painter.line_segment(
            [
                egui::pos2(rect.left(), axis_bottom),
                egui::pos2(rect.right(), axis_bottom),
            ],
            egui::Stroke::new(1.0, egui::Color32::from_gray(120)),
        );
        let threshold_x = rect.left() + rect.width() * min_seed_hit_fraction as f32;
        painter.line_segment(
            [
                egui::pos2(threshold_x, plot_top),
                egui::pos2(threshold_x, axis_bottom),
            ],
            egui::Stroke::new(1.0, egui::Color32::from_rgb(220, 38, 38)),
        );
        painter.text(
            egui::pos2(rect.left() + 4.0, rect.bottom() - 13.0),
            egui::Align2::LEFT_TOP,
            "0.0",
            egui::FontId::monospace(10.0),
            egui::Color32::from_gray(130),
        );
        painter.text(
            egui::pos2(rect.right() - 4.0, rect.bottom() - 13.0),
            egui::Align2::RIGHT_TOP,
            "1.0",
            egui::FontId::monospace(10.0),
            egui::Color32::from_gray(130),
        );
        painter.text(
            egui::pos2(threshold_x + 2.0, plot_top + 1.0),
            egui::Align2::LEFT_TOP,
            format!("{min_seed_hit_fraction:.2}"),
            egui::FontId::monospace(9.0),
            egui::Color32::from_rgb(220, 38, 38),
        );
        let mut pointed_bin: Option<usize> = None;
        if let Some(pointer) = response.hover_pos() {
            if pointer.x >= rect.left()
                && pointer.x <= rect.right()
                && pointer.y >= plot_top
                && pointer.y <= axis_bottom
            {
                let mut idx =
                    ((pointer.x - rect.left()) / rect.width() * bin_count as f32).floor() as usize;
                if idx >= bin_count {
                    idx = bin_count.saturating_sub(1);
                }
                pointed_bin = Some(idx);
                if let Some(count) = score_density_bins.get(idx) {
                    let left = idx as f64 / bin_count as f64;
                    let right = (idx + 1) as f64 / bin_count as f64;
                    let retained_count = saved_report
                        .as_ref()
                        .map(|report| {
                            Self::select_rna_read_report_score_bin_record_indices(
                                report,
                                idx,
                                bin_count,
                                score_density_variant,
                                score_density_seed_filter_override.as_ref(),
                            )
                            .len()
                        })
                        .unwrap_or(0);
                    response.clone().on_hover_ui_at_pointer(|ui| {
                        ui.monospace(format!("bin {idx}: [{left:.3}, {right:.3})"));
                        ui.monospace(format!(
                            "{} in histogram: {count}",
                            Self::rna_read_score_density_variant_label(score_density_variant)
                        ));
                        if saved_report.is_some() {
                            ui.monospace(format!(
                                "retained saved-report rows in this bin: {retained_count}"
                            ));
                        }
                        ui.small(
                            "Click to focus the retained saved-report rows from this score bin.",
                        );
                    });
                }
            }
        }
        if response.clicked()
            && let Some(idx) = pointed_bin
        {
            let tested_count = score_density_bins.get(idx).copied().unwrap_or(0);
            if let Some(report) = saved_report {
                let selected_record_indices = Self::select_rna_read_report_score_bin_record_indices(
                    &report,
                    idx,
                    bin_count,
                    score_density_variant,
                    score_density_seed_filter_override.as_ref(),
                );
                self.rna_read_alignment_effect_score_bin_index = Some(idx);
                self.rna_read_alignment_effect_filter = RnaReadAlignmentEffectFilter::SelectedOnly;
                self.rna_read_alignment_effect_sort_key = RnaReadAlignmentEffectSortKey::Score;
                self.rna_read_alignment_effect_search.clear();
                self.rna_seed_selected_record_indices =
                    selected_record_indices.iter().copied().collect();
                self.rna_seed_highlight_record_index = selected_record_indices.first().copied();
                if report.read_count_aligned > 0 {
                    self.rna_read_statistics_tab = RnaReadEvidenceSourceTab::MappedCdna;
                    self.rna_read_mapped_cdna_subview = RnaReadMappedCdnaSubview::ReadEffects;
                } else {
                    self.rna_read_statistics_tab = RnaReadEvidenceSourceTab::ThresholdedCdna;
                }
                self.op_status = if selected_record_indices.is_empty() {
                    format!(
                        "Score bin {} contains {} {} in the histogram, but 0 retained saved-report rows. The histogram and saved report use different populations here; the report currently stores only the retained top {} row(s).",
                        Self::format_rna_read_score_bin_spec(idx, bin_count),
                        tested_count,
                        Self::rna_read_score_density_variant_label(score_density_variant),
                        report.hits.len()
                    )
                } else if report.read_count_aligned > 0 {
                    format!(
                        "Focused mapped read effects on {} retained saved-report row(s) from {} score bin {}",
                        selected_record_indices.len(),
                        Self::rna_read_score_density_variant_label(score_density_variant),
                        Self::format_rna_read_score_bin_spec(idx, bin_count)
                    )
                } else {
                    format!(
                        "Selected {} retained saved-report row(s) from {} score bin {}. Phase 2 has not been run yet, so mapped read effects stay empty until alignment.",
                        selected_record_indices.len(),
                        Self::rna_read_score_density_variant_label(score_density_variant),
                        Self::format_rna_read_score_bin_spec(idx, bin_count)
                    )
                };
            } else {
                self.op_status =
                    "Save or load a Report ID before focusing a score-density bin".to_string();
            }
        }
        ui.small(format!(
            "Score bins: {} | max bin count: {}",
            score_density_bins.len(),
            max_count
        ));
        if let Some(idx) = selected_bin_index {
            ui.small(format!(
                "Selected score bin: {}",
                Self::format_rna_read_score_bin_spec(idx, bin_count)
            ));
        } else {
            ui.small("Selected score bin: <none>");
        }
        match score_density_variant {
            RnaReadScoreDensityVariant::AllScored => ui.small(format!(
                "Red line = raw min_hit threshold ({min_seed_hit_fraction:.2}) only. Crossing it is necessary but not sufficient for the full composite seed gate; weighted support, unique matched k-mers, chain consistency, transcript-gap, and transition requirements are checked separately."
            )),
            RnaReadScoreDensityVariant::CompositeSeedGate => ui.small(format!(
                "Bars already respect the full composite seed gate recorded for this run. The red line still marks the raw min_hit component ({min_seed_hit_fraction:.2}) inside that stricter gate."
            )),
            RnaReadScoreDensityVariant::RetainedReplayCurrentControls => ui.small(format!(
                "Bars replay the full composite seed gate over retained saved-report rows only, using the current controls. The red line marks the current raw min_hit component ({min_seed_hit_fraction:.2}); reads never retained by the original run are not revisited."
            )),
        };
    }

    pub(super) fn render_rna_read_statistics_tabs(
        &mut self,
        ui: &mut egui::Ui,
        view: &SplicingExpertView,
        progress: &RnaReadInterpretProgress,
        allow_mapping_actions: bool,
    ) {
        let has_mapped_support = progress.aligned > 0
            || !progress.mapped_exon_support_frequencies.is_empty()
            || !progress.mapped_junction_support_frequencies.is_empty()
            || !progress.mapped_isoform_support_rows.is_empty();
        ui.horizontal_wrapped(|ui| {
            ui.small("Analysis source:");
            ui.selectable_value(
                &mut self.rna_read_statistics_tab,
                RnaReadEvidenceSourceTab::ReportedTranscript,
                "Reported transcript",
            )
            .on_hover_text(
                "Show the annotated transcript model for the selected locus: reported exons, reported junctions, and the transcript catalogue admitted by the current splicing view.",
            );
            ui.selectable_value(
                &mut self.rna_read_statistics_tab,
                RnaReadEvidenceSourceTab::ThresholdedCdna,
                "Thresholded cDNA",
            )
            .on_hover_text(
                "Show phase-1 retained/seed-passed cDNA support. These rows are useful for early triage but remain inference-driven rather than alignment-confirmed.",
            );
            ui.selectable_value(
                &mut self.rna_read_statistics_tab,
                RnaReadEvidenceSourceTab::MappedCdna,
                "Mapped cDNA",
            )
            .on_hover_text(
                "Show phase-2 support derived from retained-read alignments (best mappings, exon/junction overlap, mapped isoform ranking).",
            );
            if !has_mapped_support {
                ui.small(
                    egui::RichText::new(
                        "Mapped cDNA stays empty until phase 2 produces aligned reads.",
                    )
                    .color(egui::Color32::from_rgb(100, 116, 139)),
                );
            }
        });
        match self.rna_read_statistics_tab {
            RnaReadEvidenceSourceTab::ReportedTranscript => {
                self.render_reported_transcript_support_tables(ui, view);
            }
            RnaReadEvidenceSourceTab::ThresholdedCdna => {
                self.render_thresholded_cdna_support_tables(ui, view, progress);
            }
            RnaReadEvidenceSourceTab::MappedCdna => {
                self.render_rna_read_mapped_cdna_panels(ui, view, progress, allow_mapping_actions);
            }
        }
    }

    pub(super) fn default_rna_read_effect_table_height(ui: &egui::Ui) -> f32 {
        let row_height = ui.text_style_height(&egui::TextStyle::Monospace).max(14.0) + 4.0;
        (row_height * 15.0 + 36.0).clamp(320.0, 560.0)
    }

    pub(super) fn default_rna_read_support_table_height(ui: &egui::Ui) -> f32 {
        let row_height = ui.text_style_height(&egui::TextStyle::Monospace).max(14.0) + 4.0;
        (row_height * 8.0 + 36.0).clamp(180.0, 260.0)
    }

    pub(super) fn target_rna_read_preview_body_rows(row_count: usize) -> usize {
        row_count.clamp(12, 16)
    }

    pub(super) fn default_rna_read_preview_table_height_for_rows(
        ui: &egui::Ui,
        row_count: usize,
    ) -> f32 {
        let row_height = ui.text_style_height(&egui::TextStyle::Monospace).max(14.0) + 4.0;
        let visible_rows = Self::target_rna_read_preview_body_rows(row_count) as f32 + 1.0;
        (row_height * visible_rows + 36.0).clamp(260.0, 520.0)
    }

    pub(super) fn virtual_rna_read_table_row_height(ui: &egui::Ui) -> f32 {
        ui.spacing().interact_size.y.max(18.0)
    }

    pub(super) fn format_rna_read_target_coverage_summary(
        aligned_target_bp: usize,
        target_length_bp: usize,
        target_coverage_fraction: f64,
    ) -> String {
        format!(
            "{} / {} bp ({:.1}%)",
            aligned_target_bp,
            target_length_bp,
            target_coverage_fraction * 100.0
        )
    }

    pub(super) fn rna_read_full_length_tooltip() -> &'static str {
        "Full-length classes: exact = 100% template coverage; near = >=95% template coverage; strict_end = near plus both template ends within 15 bp and identity >= active alignment threshold."
    }

    pub(super) fn rna_read_full_length_class_color(label: &str) -> egui::Color32 {
        match label {
            "exact" => egui::Color32::from_rgb(4, 120, 87),
            "strict_end" => egui::Color32::from_rgb(22, 163, 74),
            "near" => egui::Color32::from_rgb(59, 130, 246),
            _ => egui::Color32::from_rgb(100, 116, 139),
        }
    }

    pub(super) fn rna_read_target_gene_ids_for_view(
        report: &RnaReadInterpretationReport,
        view: &SplicingExpertView,
    ) -> Vec<String> {
        let mut gene_ids = if report.target_gene_ids.is_empty() {
            let trimmed = view.group_label.trim();
            if trimmed.is_empty() {
                vec![]
            } else {
                vec![trimmed.to_string()]
            }
        } else {
            report.target_gene_ids.clone()
        };
        gene_ids.sort_by(|left, right| left.to_ascii_lowercase().cmp(&right.to_ascii_lowercase()));
        gene_ids.dedup_by(|left, right| left.eq_ignore_ascii_case(right));
        gene_ids
    }

    pub(super) fn format_length_distribution_stats_line(
        summary: &RnaReadLengthDistributionSummary,
    ) -> String {
        format!(
            "mean={:.1} bp | q0={} bp | q25={} bp | q50={} bp | q75={} bp | q100={} bp",
            summary.mean_length_bp,
            summary.min_length_bp,
            summary.q25_length_bp,
            summary.median_length_bp,
            summary.q75_length_bp,
            summary.max_length_bp,
        )
    }

    pub(super) fn render_rna_read_target_share_by_length_series(
        &self,
        ui: &mut egui::Ui,
        label: &str,
        total_length_counts: &[u64],
        target_length_counts: &[u64],
        bar_color: egui::Color32,
    ) {
        let bins = GentleEngine::auto_bin_read_length_counts(total_length_counts, 24);
        ui.horizontal_wrapped(|ui| {
            ui.small(
                egui::RichText::new(label)
                    .strong()
                    .color(egui::Color32::from_rgb(51, 65, 85)),
            );
            ui.separator();
            ui.small("y = target-positive share within each read-length bin");
        });
        if bins.is_empty() {
            ui.small("No reads in this subset yet.");
            ui.add_space(4.0);
            return;
        }
        let desired = Vec2::new(ui.available_width().max(300.0), 52.0);
        let (rect, response) = ui.allocate_exact_size(desired, egui::Sense::hover());
        let painter = ui.painter_at(rect);
        painter.rect_stroke(
            rect,
            2.0,
            egui::Stroke::new(1.0, egui::Color32::from_gray(130)),
            egui::StrokeKind::Inside,
        );
        let axis_bottom = rect.bottom() - 1.5;
        let plot_top = rect.top() + 2.0;
        let plot_h = (axis_bottom - plot_top).max(1.0);
        let bin_w = rect.width() / bins.len().max(1) as f32;
        for (idx, (start, end, total_count)) in bins.iter().enumerate() {
            if *total_count == 0 {
                continue;
            }
            let target_count = target_length_counts
                .iter()
                .enumerate()
                .skip(*start)
                .take(end.saturating_sub(*start).saturating_add(1))
                .map(|(_, count)| *count)
                .sum::<u64>();
            let share = target_count as f32 / (*total_count).max(1) as f32;
            let x0 = rect.left() + idx as f32 * bin_w + 0.6;
            let x1 = if idx + 1 == bins.len() {
                rect.right() - 0.6
            } else {
                rect.left() + (idx + 1) as f32 * bin_w - 0.6
            };
            if x1 <= x0 {
                continue;
            }
            let bar_rect = egui::Rect::from_min_max(
                egui::pos2(x0, axis_bottom - share * plot_h),
                egui::pos2(x1, axis_bottom),
            );
            painter.rect_filled(bar_rect, 0.0, bar_color);
        }
        if let Some(pointer) = response.hover_pos()
            && pointer.x >= rect.left()
            && pointer.x <= rect.right()
            && pointer.y >= plot_top
            && pointer.y <= axis_bottom
        {
            let mut idx =
                ((pointer.x - rect.left()) / rect.width() * bins.len() as f32).floor() as usize;
            idx = idx.min(bins.len().saturating_sub(1));
            if let Some((start, end, total_count)) = bins.get(idx) {
                let target_count = target_length_counts
                    .iter()
                    .enumerate()
                    .skip(*start)
                    .take(end.saturating_sub(*start).saturating_add(1))
                    .map(|(_, count)| *count)
                    .sum::<u64>();
                let share = if *total_count == 0 {
                    0.0
                } else {
                    target_count as f64 / *total_count as f64
                };
                response.on_hover_ui_at_pointer(|ui| {
                    if start == end {
                        ui.monospace(format!("{start} bp"));
                    } else {
                        ui.monospace(format!("{start}..{end} bp"));
                    }
                    ui.monospace(format!(
                        "target-positive: {} / {} ({:.1}%)",
                        target_count,
                        total_count,
                        share * 100.0
                    ));
                });
            }
        }
        ui.add_space(4.0);
    }

    pub(super) fn render_rna_read_fraction_distribution_series(
        &self,
        ui: &mut egui::Ui,
        label: &str,
        bin_counts: &[u64],
        bar_color: egui::Color32,
    ) {
        let total = bin_counts.iter().copied().sum::<u64>();
        ui.horizontal_wrapped(|ui| {
            ui.small(
                egui::RichText::new(label)
                    .strong()
                    .color(egui::Color32::from_rgb(51, 65, 85)),
            );
            ui.separator();
            ui.small(format!("samples={}", Self::format_count_compact_km(total)));
            ui.separator();
            ui.small("x = fraction of the read covered by target-assigned fragment");
        });
        if total == 0 || bin_counts.is_empty() {
            ui.small("No target-positive fragments yet.");
            ui.add_space(4.0);
            return;
        }
        let desired = Vec2::new(ui.available_width().max(300.0), 52.0);
        let (rect, response) = ui.allocate_exact_size(desired, egui::Sense::hover());
        let painter = ui.painter_at(rect);
        painter.rect_stroke(
            rect,
            2.0,
            egui::Stroke::new(1.0, egui::Color32::from_gray(130)),
            egui::StrokeKind::Inside,
        );
        let max_count = bin_counts.iter().copied().max().unwrap_or(1) as f32;
        let axis_bottom = rect.bottom() - 1.5;
        let plot_top = rect.top() + 2.0;
        let plot_h = (axis_bottom - plot_top).max(1.0);
        let bin_w = rect.width() / bin_counts.len().max(1) as f32;
        for (idx, count) in bin_counts.iter().enumerate() {
            if *count == 0 {
                continue;
            }
            let x0 = rect.left() + idx as f32 * bin_w + 0.6;
            let x1 = if idx + 1 == bin_counts.len() {
                rect.right() - 0.6
            } else {
                rect.left() + (idx + 1) as f32 * bin_w - 0.6
            };
            if x1 <= x0 {
                continue;
            }
            let h = (*count as f32 / max_count) * plot_h;
            let bar_rect = egui::Rect::from_min_max(
                egui::pos2(x0, axis_bottom - h),
                egui::pos2(x1, axis_bottom),
            );
            painter.rect_filled(bar_rect, 0.0, bar_color);
        }
        if let Some(pointer) = response.hover_pos()
            && pointer.x >= rect.left()
            && pointer.x <= rect.right()
            && pointer.y >= plot_top
            && pointer.y <= axis_bottom
        {
            let mut idx = ((pointer.x - rect.left()) / rect.width() * bin_counts.len() as f32)
                .floor() as usize;
            idx = idx.min(bin_counts.len().saturating_sub(1));
            let left = idx as f64 / bin_counts.len().saturating_sub(1).max(1) as f64;
            let right = ((idx + 1).min(bin_counts.len().saturating_sub(1))) as f64
                / bin_counts.len().saturating_sub(1).max(1) as f64;
            response.on_hover_ui_at_pointer(|ui| {
                ui.monospace(format!("{:.0}%..{:.0}%", left * 100.0, right * 100.0));
                ui.monospace(format!("fragments: {}", bin_counts[idx]));
            });
        }
        ui.add_space(4.0);
    }

    pub(super) fn render_rna_read_length_distributions_panel(
        &mut self,
        ui: &mut egui::Ui,
        view: &SplicingExpertView,
        report: Option<&RnaReadInterpretationReport>,
    ) {
        egui::CollapsingHeader::new("Read length distributions")
            .default_open(true)
            .show(ui, |ui| {
                let Some(report) = report else {
                    ui.small(
                        "Load or create a saved RNA-read report to view read-length distributions.",
                    );
                    return;
                };
                ui.small(
                    "Engine stores exact 1-bp read-length counts; this panel renders adaptive bins for quick visual inspection.",
                );
                self.render_rna_read_length_distribution_series(
                    ui,
                    "All encountered reads",
                    &report.read_length_counts_all,
                    egui::Color32::from_rgb(71, 85, 105),
                );
                self.render_rna_read_length_distribution_series(
                    ui,
                    "Seed-passed subset",
                    &report.read_length_counts_seed_passed,
                    egui::Color32::from_rgb(2, 132, 199),
                );
                self.render_rna_read_length_distribution_series(
                    ui,
                    "Aligned subset",
                    &report.read_length_counts_aligned,
                    egui::Color32::from_rgb(22, 163, 74),
                );
                self.render_rna_read_length_distribution_series(
                    ui,
                    "Full-length exact",
                    &report.read_length_counts_full_length_exact,
                    egui::Color32::from_rgb(4, 120, 87),
                );
                self.render_rna_read_length_distribution_series(
                    ui,
                    "Full-length near",
                    &report.read_length_counts_full_length_near,
                    egui::Color32::from_rgb(14, 165, 233),
                );
                self.render_rna_read_length_distribution_series(
                    ui,
                    "Full-length strict_end",
                    &report.read_length_counts_full_length_strict,
                    egui::Color32::from_rgb(21, 128, 61),
                );

                let target_gene_ids = Self::rna_read_target_gene_ids_for_view(report, view);
                if target_gene_ids.is_empty() {
                    ui.small(
                        egui::RichText::new(
                            "No target gene/group label is available for the target-fragment quality summary yet.",
                        )
                        .color(egui::Color32::from_rgb(180, 83, 9)),
                    );
                    return;
                }
                let target_summary = match self.saved_rna_read_gene_support_summary_for_report_id(
                    &report.report_id,
                    target_gene_ids.clone(),
                    vec![],
                    RnaReadGeneSupportCompleteRule::Near,
                ) {
                    Ok(summary) => summary,
                    Err(err) => {
                        ui.small(
                            egui::RichText::new(format!(
                                "Target-fragment quality summary unavailable: {err}"
                            ))
                            .color(egui::Color32::from_rgb(180, 83, 9)),
                        );
                        return;
                    }
                };
                let target_label = if target_summary.matched_gene_ids.is_empty() {
                    target_gene_ids.join(", ")
                } else {
                    target_summary.matched_gene_ids.join(", ")
                };
                ui.separator();
                ui.small(
                    egui::RichText::new(format!(
                        "Target-fragment quality for {target_label} (accepted target-support reads under complete_rule=near).",
                    ))
                    .color(egui::Color32::from_rgb(30, 41, 59)),
                );
                ui.small(format!(
                    "Target fragment/read coverage: n={} mean={:.1}% median={:.1}% p95={:.1}%",
                    Self::format_count_compact_km(
                        target_summary.accepted_target_query_coverage.sample_count as u64
                    ),
                    target_summary.accepted_target_query_coverage.mean_fraction * 100.0,
                    target_summary.accepted_target_query_coverage.median_fraction * 100.0,
                    target_summary.accepted_target_query_coverage.p95_fraction * 100.0,
                ));
                if ui
                    .small_button("Export target quality...")
                    .on_hover_text(
                        "Export this target-fragment quality summary as SVG or JSON through the shared engine comparison route.",
                    )
                    .clicked()
                {
                    self.export_rna_read_target_quality(
                        report.report_id.clone(),
                        target_gene_ids.clone(),
                        RnaReadGeneSupportCompleteRule::Near,
                    );
                }
                self.render_rna_read_target_share_by_length_series(
                    ui,
                    &format!("{target_label} share by read length"),
                    &report.read_length_counts_all,
                    &target_summary.accepted_target_read_lengths.length_counts,
                    egui::Color32::from_rgb(5, 150, 105),
                );
                self.render_rna_read_length_distribution_series(
                    ui,
                    &format!("{target_label} positive reads"),
                    &target_summary.accepted_target_read_lengths.length_counts,
                    egui::Color32::from_rgb(5, 150, 105),
                );
                self.render_rna_read_length_distribution_series(
                    ui,
                    &format!("{target_label} fragment lengths"),
                    &target_summary.accepted_target_fragment_lengths.length_counts,
                    egui::Color32::from_rgb(13, 148, 136),
                );
                self.render_rna_read_fraction_distribution_series(
                    ui,
                    &format!("{target_label} fragment/read coverage"),
                    &target_summary.accepted_target_query_coverage.bin_counts,
                    egui::Color32::from_rgb(2, 132, 199),
                );
            });
    }

    pub(super) fn render_rna_read_length_distribution_series(
        &self,
        ui: &mut egui::Ui,
        label: &str,
        length_counts: &[u64],
        bar_color: egui::Color32,
    ) {
        let summary = GentleEngine::summarize_read_length_distribution(length_counts);
        let bins = GentleEngine::auto_bin_read_length_counts(length_counts, 24);
        let compact = GentleEngine::format_read_length_distribution_compact(length_counts, 24, 8);
        ui.horizontal_wrapped(|ui| {
            ui.small(
                egui::RichText::new(label)
                    .strong()
                    .color(egui::Color32::from_rgb(51, 65, 85)),
            );
            ui.separator();
            ui.small(format!(
                "reads={}",
                Self::format_count_compact_km(summary.sample_count as u64)
            ));
            if let (Some(first), Some(last)) = (bins.first(), bins.last()) {
                ui.separator();
                ui.small(format!("range={}..{} bp", first.0, last.1));
            }
            ui.separator();
            ui.small(format!("bins={}", bins.len()));
        });
        if bins.is_empty() {
            ui.small("No reads in this subset yet.");
            ui.add_space(2.0);
            return;
        }
        ui.small(Self::format_length_distribution_stats_line(&summary));

        let desired = Vec2::new(ui.available_width().max(300.0), 52.0);
        let (rect, response) = ui.allocate_exact_size(desired, egui::Sense::hover());
        let painter = ui.painter_at(rect);
        painter.rect_stroke(
            rect,
            2.0,
            egui::Stroke::new(1.0, egui::Color32::from_gray(130)),
            egui::StrokeKind::Inside,
        );
        let max_count = bins.iter().map(|(_, _, count)| *count).max().unwrap_or(1) as f32;
        let axis_bottom = rect.bottom() - 1.5;
        let plot_top = rect.top() + 2.0;
        let plot_h = (axis_bottom - plot_top).max(1.0);
        let bin_w = rect.width() / bins.len().max(1) as f32;
        for (idx, (_, _, count)) in bins.iter().enumerate() {
            if *count == 0 {
                continue;
            }
            let x0 = rect.left() + idx as f32 * bin_w + 0.6;
            let x1 = if idx + 1 == bins.len() {
                rect.right() - 0.6
            } else {
                rect.left() + (idx + 1) as f32 * bin_w - 0.6
            };
            if x1 <= x0 {
                continue;
            }
            let h = (*count as f32 / max_count) * plot_h;
            let bar_rect = egui::Rect::from_min_max(
                egui::pos2(x0, axis_bottom - h),
                egui::pos2(x1, axis_bottom),
            );
            painter.rect_filled(bar_rect, 0.0, bar_color);
        }
        if let Some(pointer) = response.hover_pos()
            && pointer.x >= rect.left()
            && pointer.x <= rect.right()
            && pointer.y >= plot_top
            && pointer.y <= axis_bottom
        {
            let mut idx =
                ((pointer.x - rect.left()) / rect.width() * bins.len() as f32).floor() as usize;
            idx = idx.min(bins.len().saturating_sub(1));
            if let Some((start, end, count)) = bins.get(idx) {
                response.on_hover_ui_at_pointer(|ui| {
                    if start == end {
                        ui.monospace(format!("{start} bp"));
                    } else {
                        ui.monospace(format!("{start}..{end} bp"));
                    }
                    ui.monospace(format!("reads: {count}"));
                });
            }
        }
        ui.small(format!("Auto bins: {compact}"));
        ui.add_space(4.0);
    }

    pub(super) fn rna_read_fragment_alignment_note(
        target_coverage_fraction: f64,
        aligned_target_bp: usize,
        target_length_bp: usize,
    ) -> Option<String> {
        if target_length_bp == 0 || target_coverage_fraction >= 0.20 {
            return None;
        }
        Some(format!(
            "Fragment-only confirmation: the aligned transcript span covers only {} of the transcript template. This can still be biologically useful when the fragment crosses an exon boundary, but it is not a near-full-length transcript confirmation.",
            Self::format_rna_read_target_coverage_summary(
                aligned_target_bp,
                target_length_bp,
                target_coverage_fraction,
            )
        ))
    }

    pub(super) fn render_rna_read_mapped_cdna_panels(
        &mut self,
        ui: &mut egui::Ui,
        view: &SplicingExpertView,
        progress: &RnaReadInterpretProgress,
        allow_mapping_actions: bool,
    ) {
        ui.horizontal_wrapped(|ui| {
            ui.small("Mapped cDNA view:");
            ui.selectable_value(
                &mut self.rna_read_mapped_cdna_subview,
                RnaReadMappedCdnaSubview::ReadEffects,
                "Read effects",
            )
            .on_hover_text(
                "Read-first aligned cDNA inspection: one row per aligned retained read, showing phase-1 interpretation, phase-2 alignment outcome, and mapped exon/junction contributions.",
            );
            ui.selectable_value(
                &mut self.rna_read_mapped_cdna_subview,
                RnaReadMappedCdnaSubview::AggregateSupport,
                "Aggregate support",
            )
            .on_hover_text(
                "Aggregate exon/junction/isoform support derived from phase-2 best mappings.",
            );
        });

        let saved_report = self.current_saved_rna_read_report();
        match self.rna_read_mapped_cdna_subview {
            RnaReadMappedCdnaSubview::ReadEffects => {
                self.render_rna_read_mapped_read_effects(
                    ui,
                    view,
                    progress,
                    saved_report.as_deref(),
                    allow_mapping_actions,
                );
            }
            RnaReadMappedCdnaSubview::AggregateSupport => {
                self.render_rna_read_mapped_aggregate_support_tables(
                    ui,
                    progress,
                    saved_report.as_deref(),
                );
            }
        }
    }

    pub(super) fn render_rna_read_mapped_read_effects(
        &mut self,
        ui: &mut egui::Ui,
        view: &SplicingExpertView,
        progress: &RnaReadInterpretProgress,
        report: Option<&RnaReadInterpretationReport>,
        allow_mapping_actions: bool,
    ) {
        ui.small(
            "Read effects compare phase-1 thresholded interpretation against phase-2 pairwise transcript alignment for each aligned retained read.",
        );
        ui.small(
            "Workflow: 1) click a score-density bar above to define a formal score_bin subset, 2) refine with filter/search, 3) select a row to inspect the pairwise alignment, then open the external dotplot only when you actually want it.",
        );
        let Some(report) = report else {
            ui.small(
                egui::RichText::new(
                    "No saved report is loaded yet, so aligned read effects cannot be inspected.",
                )
                .color(egui::Color32::from_rgb(180, 83, 9)),
            );
            return;
        };
        let inspection = match self
            .current_saved_rna_read_alignment_inspection(report.hits.len(), None)
        {
            Ok(inspection) => inspection,
            Err(message) => {
                ui.small(egui::RichText::new(message).color(egui::Color32::from_rgb(180, 83, 9)));
                return;
            }
        };
        let summary = Self::summarize_rna_read_alignment_effects(report, &inspection);
        let hits_by_record_index = report
            .hits
            .iter()
            .map(|hit| (hit.record_index, hit))
            .collect::<HashMap<_, _>>();
        ui.horizontal_wrapped(|ui| {
            ui.small(format!("Aligned rows: {}", summary.aligned_rows));
            ui.separator();
            ui.small(format!("Confirmed: {}", summary.confirmed_assignments));
            ui.separator();
            ui.small(format!("Reassigned: {}", summary.reassigned_transcripts));
            ui.separator();
            ui.small(format!(
                "Aligned with no phase-1 tx: {}",
                summary.aligned_without_phase1_assignment
            ));
            ui.separator();
            ui.small(format!(
                "Seed-passed but unaligned: {}",
                summary.seed_passed_but_unaligned
            ));
        });

        let mut focus_disagreements = false;
        let mut focus_max_score_outliers = false;
        let mut focus_rightmost_bin = false;
        let mut reset_read_effect_view = false;
        ui.horizontal_wrapped(|ui| {
            ui.small("Subset shortcuts:");
            if ui
                .small_button("Disagreements")
                .on_hover_text(
                    "Focus non-confirmed aligned reads: transcript reassignments plus rows aligned without a phase-1 transcript assignment.",
                )
                .clicked()
            {
                focus_disagreements = true;
            }
            if ui
                .small_button("Max-score ties")
                .on_hover_text(
                    "Select only saved-report reads tied at the maximal phase-1 seed score, then show that subset in the read-effects table.",
                )
                .clicked()
            {
                focus_max_score_outliers = true;
            }
            if ui
                .small_button("Rightmost score bin")
                .on_hover_text(
                    "Select the broader outlier subset from the rightmost non-empty score-density bin, then focus the table on that subset.",
                )
                .clicked()
            {
                focus_rightmost_bin = true;
            }
            if ui
                .small_button("Reset view")
                .on_hover_text(
                    "Return to the full aligned-read view with rank ordering and no search filter.",
                )
                .clicked()
            {
                reset_read_effect_view = true;
            }
        });

        if reset_read_effect_view {
            self.rna_read_alignment_effect_filter = RnaReadAlignmentEffectFilter::AllAligned;
            self.rna_read_alignment_effect_sort_key = RnaReadAlignmentEffectSortKey::Rank;
            self.rna_read_alignment_effect_search.clear();
            self.clear_rna_read_alignment_score_bin_focus();
        }
        if focus_disagreements {
            self.rna_read_alignment_effect_filter = RnaReadAlignmentEffectFilter::DisagreementOnly;
            self.rna_read_alignment_effect_sort_key = RnaReadAlignmentEffectSortKey::Score;
            self.rna_read_alignment_effect_search.clear();
        }
        if focus_max_score_outliers || focus_rightmost_bin {
            let score_density_seed_filter_override = self
                .current_rna_read_score_density_seed_filter_override(
                    self.rna_read_evidence_ui.score_density_variant,
                )
                .ok()
                .flatten();
            let selected_record_indices = if focus_max_score_outliers {
                Self::select_rna_read_report_max_score_record_indices(report)
            } else {
                Self::select_rna_read_report_rightmost_score_bin_record_indices(
                    report,
                    self.rna_read_evidence_ui.score_density_variant,
                    score_density_seed_filter_override.as_ref(),
                )
            };
            if focus_rightmost_bin {
                let bin_count = report.score_density_bins.len().max(40);
                let score_bin_index = selected_record_indices.first().and_then(|record_index| {
                    report
                        .hits
                        .iter()
                        .find(|hit| hit.record_index == *record_index)
                        .map(|hit| {
                            Self::rna_read_score_density_bin_index(hit.seed_hit_fraction, bin_count)
                        })
                });
                self.rna_read_alignment_effect_score_bin_index = score_bin_index;
            } else {
                self.clear_rna_read_alignment_score_bin_focus();
            }
            self.rna_seed_selected_record_indices = selected_record_indices
                .iter()
                .copied()
                .collect::<BTreeSet<_>>();
            self.rna_read_alignment_effect_filter = RnaReadAlignmentEffectFilter::SelectedOnly;
            self.rna_read_alignment_effect_sort_key = RnaReadAlignmentEffectSortKey::Score;
            self.rna_read_alignment_effect_search.clear();
            self.rna_seed_highlight_record_index = selected_record_indices.first().copied();
            self.op_status = if focus_max_score_outliers {
                format!(
                    "Focused read effects on {} maximal-score outlier read(s)",
                    self.rna_seed_selected_record_indices.len()
                )
            } else {
                format!(
                    "Focused read effects on {} rightmost-bin outlier read(s)",
                    self.rna_seed_selected_record_indices.len()
                )
            };
        }

        ui.horizontal_wrapped(|ui| {
            ui.small("Filter:");
            egui::ComboBox::from_id_salt(format!(
                "rna_alignment_effect_filter_{}",
                report.report_id
            ))
            .selected_text(Self::rna_read_alignment_effect_filter_label(
                self.rna_read_alignment_effect_filter,
            ))
            .show_ui(ui, |ui| {
                for filter in [
                    RnaReadAlignmentEffectFilter::AllAligned,
                    RnaReadAlignmentEffectFilter::ConfirmedOnly,
                    RnaReadAlignmentEffectFilter::DisagreementOnly,
                    RnaReadAlignmentEffectFilter::ReassignedOnly,
                    RnaReadAlignmentEffectFilter::NoPhase1Only,
                    RnaReadAlignmentEffectFilter::SelectedOnly,
                ] {
                    ui.selectable_value(
                        &mut self.rna_read_alignment_effect_filter,
                        filter,
                        Self::rna_read_alignment_effect_filter_label(filter),
                    );
                }
            });
            ui.separator();
            ui.small("Sort:");
            egui::ComboBox::from_id_salt(format!("rna_alignment_effect_sort_{}", report.report_id))
                .selected_text(Self::rna_read_alignment_effect_sort_key_label(
                    self.rna_read_alignment_effect_sort_key,
                ))
                .show_ui(ui, |ui| {
                    for sort_key in [
                        RnaReadAlignmentEffectSortKey::Rank,
                        RnaReadAlignmentEffectSortKey::Identity,
                        RnaReadAlignmentEffectSortKey::Coverage,
                        RnaReadAlignmentEffectSortKey::Score,
                    ] {
                        ui.selectable_value(
                            &mut self.rna_read_alignment_effect_sort_key,
                            sort_key,
                            Self::rna_read_alignment_effect_sort_key_label(sort_key),
                        );
                    }
                });
            ui.separator();
            ui.small("Search:");
            ui.add(
                egui::TextEdit::singleline(&mut self.rna_read_alignment_effect_search)
                    .desired_width(220.0)
                    .hint_text("read, tx, effect, #index"),
            );
            if !self.rna_read_alignment_effect_search.trim().is_empty()
                && ui.small_button("Clear").clicked()
            {
                self.rna_read_alignment_effect_search.clear();
            }
        });
        let current_subset_spec = self.current_rna_read_alignment_subset_spec();
        let filtered_inspection = match self.current_saved_rna_read_alignment_inspection(
            report.hits.len(),
            Some(current_subset_spec),
        ) {
            Ok(inspection) => inspection,
            Err(message) => {
                ui.small(egui::RichText::new(message).color(egui::Color32::from_rgb(180, 83, 9)));
                return;
            }
        };
        if focus_disagreements {
            self.rna_seed_highlight_record_index =
                filtered_inspection.rows.first().map(|row| row.record_index);
            self.op_status = format!(
                "Focused read effects on {} non-confirmed aligned row(s)",
                filtered_inspection.subset_match_count
            );
        }
        if focus_max_score_outliers || focus_rightmost_bin {
            self.rna_seed_highlight_record_index =
                filtered_inspection.rows.first().map(|row| row.record_index);
        }
        let displayed_rows = filtered_inspection.rows.iter().collect::<Vec<_>>();
        let displayed_record_indices =
            Self::collect_rna_read_alignment_effect_record_indices(&displayed_rows);
        let visible_count = filtered_inspection.subset_match_count;
        let score_bin_count = report.score_density_bins.len().max(40);
        let filtered_subset_spec = Self::rna_read_alignment_effect_subset_spec(
            self.rna_read_alignment_effect_filter,
            &self.rna_read_alignment_effect_search,
            self.rna_read_alignment_effect_sort_key,
            self.rna_read_evidence_ui.score_density_variant,
            self.rna_read_alignment_effect_score_bin_index,
            score_bin_count,
        );
        ui.small(format!(
            "Showing {} of {} aligned rows",
            filtered_inspection.subset_match_count, inspection.aligned_count
        ));
        ui.small(format!("Current filtered subset: {filtered_subset_spec}"));
        if let Some(score_bin_index) = self.rna_read_alignment_effect_score_bin_index {
            let selected_from_bin = Self::select_rna_read_report_score_bin_record_indices(
                report,
                score_bin_index,
                score_bin_count,
                self.rna_read_evidence_ui.score_density_variant,
                filtered_inspection
                    .subset_spec
                    .score_density_seed_filter_override
                    .as_ref(),
            );
            ui.small(format!(
                "Score-bin provenance: {} saved-report read(s) in {} bin {}, {} aligned row(s) currently match the subset.",
                selected_from_bin.len(),
                Self::rna_read_score_density_variant_label(
                    self.rna_read_evidence_ui.score_density_variant
                ),
                Self::format_rna_read_score_bin_spec(score_bin_index, score_bin_count),
                filtered_inspection.subset_match_count
            ));
        }

        ui.horizontal_wrapped(|ui| {
            ui.menu_button("Selection tools", |ui| {
                if ui
                    .button("Select max-score ties")
                    .on_hover_text(
                        "Select all saved-report rows tied for the highest retained phase-1 seed score.",
                    )
                    .clicked()
                {
                    self.clear_rna_read_alignment_score_bin_focus();
                    self.rna_seed_selected_record_indices =
                        Self::select_rna_read_report_max_score_record_indices(report)
                            .into_iter()
                            .collect::<BTreeSet<_>>();
                    self.op_status = format!(
                        "Selected {} saved-report read(s) tied at maximal seed score",
                        self.rna_seed_selected_record_indices.len()
                    );
                    ui.close();
                }
                if ui
                    .button("Select rightmost score bin")
                    .on_hover_text(
                        "Select saved-report rows from the highest non-empty score-density bin under the current density mode.",
                    )
                    .clicked()
                {
                    let bin_count = report.score_density_bins.len().max(40);
                    let selected_record_indices =
                        Self::select_rna_read_report_rightmost_score_bin_record_indices(
                            report,
                            self.rna_read_evidence_ui.score_density_variant,
                            filtered_inspection
                                .subset_spec
                                .score_density_seed_filter_override
                                .as_ref(),
                        );
                    self.rna_read_alignment_effect_score_bin_index =
                        selected_record_indices.first().and_then(|record_index| {
                            report
                                .hits
                                .iter()
                                .find(|hit| hit.record_index == *record_index)
                                .map(|hit| {
                                    Self::rna_read_score_density_bin_index(
                                        hit.seed_hit_fraction,
                                        bin_count,
                                    )
                                })
                        });
                    self.rna_seed_selected_record_indices =
                        selected_record_indices.into_iter().collect::<BTreeSet<_>>();
                    self.op_status = format!(
                        "Selected {} saved-report read(s) from the rightmost non-empty score bin",
                        self.rna_seed_selected_record_indices.len()
                    );
                    ui.close();
                }
                if ui
                    .button("Select filtered rows")
                    .on_hover_text(
                        "Select every aligned row currently visible after the active filter and search text.",
                    )
                    .clicked()
                {
                    self.clear_rna_read_alignment_score_bin_focus();
                    self.rna_seed_selected_record_indices =
                        displayed_record_indices.iter().copied().collect::<BTreeSet<_>>();
                    self.op_status = format!(
                        "Selected {} aligned read(s) from the current filtered subset ({filtered_subset_spec})",
                        self.rna_seed_selected_record_indices.len()
                    );
                    ui.close();
                }
                if ui
                    .button("Select aligned rows")
                    .on_hover_text(
                        "Select all phase-2 aligned rows in this report, ignoring the current table filter.",
                    )
                    .clicked()
                {
                    self.clear_rna_read_alignment_score_bin_focus();
                    self.rna_seed_selected_record_indices = inspection
                        .rows
                        .iter()
                        .map(|row| row.record_index)
                        .collect::<BTreeSet<_>>();
                    ui.close();
                }
                if ui
                    .button("Clear selected")
                    .on_hover_text("Clear the current saved-report row selection.")
                    .clicked()
                {
                    self.clear_rna_read_alignment_score_bin_focus();
                    self.rna_seed_selected_record_indices.clear();
                    ui.close();
                }
            });
            let selected_count = self.rna_seed_selected_record_indices.len();
            let selected_report_hits_available = selected_count > 0
                && report.hits.iter().any(|hit| {
                    self.rna_seed_selected_record_indices
                        .contains(&hit.record_index)
                });
            let highlighted_report_hit_available =
                self.selected_highlighted_rna_report_hit(report).is_some();
            if allow_mapping_actions
                && ui
                    .add_enabled(
                        self.rna_read_task.is_none() && selected_report_hits_available,
                        egui::Button::new(format!(
                            "Evaluate Selected (phase-2) [{selected_count}]"
                        )),
                    )
                    .on_hover_text(
                        "Align only checkbox-selected saved-report rows by record_index and refresh mapped read-effect summaries.",
                    )
                    .clicked()
            {
                let selected_indices = self
                    .rna_seed_selected_record_indices
                    .iter()
                    .copied()
                    .collect::<Vec<_>>();
                self.run_splicing_rna_read_alignment_phase_for_selected(view, selected_indices);
            }
            if ui
                .add_enabled(
                    selected_report_hits_available,
                    egui::Button::new(format!("Copy selected FASTA ({selected_count})")),
                )
                .on_hover_text(
                    "Copy FASTA for selected rows that are present in the active saved report.",
                )
                .clicked()
            {
                let hits =
                    Self::selected_rna_report_hits(report, &self.rna_seed_selected_record_indices);
                self.copy_rna_report_hits_as_fasta(ui, &hits, "selected report reads");
            }
            if ui
                .add_enabled(
                    highlighted_report_hit_available,
                    egui::Button::new("Copy highlighted FASTA"),
                )
                .on_hover_text(
                    "Copy FASTA for the highlighted row when it is present in the active saved report.",
                )
                .clicked()
            {
                let hits = self
                    .selected_highlighted_rna_report_hit(report)
                    .into_iter()
                    .collect::<Vec<_>>();
                self.copy_rna_report_hits_as_fasta(ui, &hits, "highlighted aligned read");
            }
            ui.menu_button(format!("Export selected ({selected_count})..."), |ui| {
                if selected_count == 0 {
                    ui.small("Select one or more aligned/audited rows first.");
                    ui.separator();
                }
                for export_kind in [
                    RnaReadSelectedExportKind::Fasta,
                    RnaReadSelectedExportKind::AlignmentsTsv,
                    RnaReadSelectedExportKind::ExonPathsTsv,
                    RnaReadSelectedExportKind::ExonAbundanceTsv,
                ] {
                    if ui
                        .add_enabled(
                            selected_report_hits_available,
                            egui::Button::new(export_kind.menu_label()),
                        )
                        .on_hover_text(export_kind.hover_text())
                        .clicked()
                    {
                        self.export_selected_rna_read_subset(export_kind);
                        ui.close();
                    }
                }
            });
            ui.menu_button(format!("Export filtered ({visible_count})..."), |ui| {
                if visible_count == 0 {
                    ui.small("No aligned rows match the current filter/search.");
                    ui.separator();
                }
                for export_kind in [
                    RnaReadSelectedExportKind::Fasta,
                    RnaReadSelectedExportKind::AlignmentsTsv,
                    RnaReadSelectedExportKind::ExonPathsTsv,
                    RnaReadSelectedExportKind::ExonAbundanceTsv,
                ] {
                    if ui
                        .add_enabled(
                            visible_count > 0,
                            egui::Button::new(export_kind.menu_label()),
                        )
                        .on_hover_text(export_kind.hover_text())
                        .clicked()
                    {
                        self.export_rna_read_subset_with_record_indices(
                            export_kind,
                            displayed_record_indices.clone(),
                            Some(filtered_subset_spec.clone()),
                            "No aligned rows match the current filter/search.",
                        );
                        ui.close();
                    }
                }
            });
            if ui
                .add_enabled(
                    selected_report_hits_available,
                    egui::Button::new("Materialize selected"),
                )
                .on_hover_text(
                    "Create project sequence entries from selected rows in the active saved report.",
                )
                .clicked()
            {
                self.materialize_selected_rna_read_report_hits();
            }
            if ui
                .add_enabled(
                    highlighted_report_hit_available,
                    egui::Button::new("Materialize highlighted"),
                )
                .on_hover_text(
                    "Create a project sequence entry from the highlighted row in the active saved report.",
                )
                .clicked()
            {
                self.materialize_highlighted_rna_read_report_hit();
            }
            ui.menu_button("Dotplots", |ui| {
                if ui
                    .button("Export dotplot for highlighted read...")
                    .on_hover_text(
                        "Write an SVG dotplot for the currently highlighted read against its mapped transcript template.",
                    )
                    .clicked()
                {
                    self.export_highlighted_rna_read_sequence_dotplot_svg(view, progress);
                    ui.close();
                }
                if ui
                    .button("Export dotplots for selected reads...")
                    .on_hover_text(
                        "Write SVG dotplots for all selected reads that can be resolved in the active report.",
                    )
                    .clicked()
                {
                    self.export_selected_rna_read_sequence_dotplot_svgs(view);
                    ui.close();
                }
            });
        });

        if inspection.rows.is_empty() {
            ui.small(
                egui::RichText::new(
                    "No aligned rows are available yet for read-level phase-1 vs phase-2 inspection.",
                )
                .color(egui::Color32::from_rgb(180, 83, 9)),
            );
            return;
        }

        let mut select_filtered_rows =
            self.rna_read_record_indices_all_selected(displayed_record_indices.iter().copied());
        let row_height = Self::virtual_rna_read_table_row_height(ui);
        let table_height = Self::default_rna_read_effect_table_height(ui);
        egui_extras::TableBuilder::new(ui)
            .id_salt(format!("rna_alignment_effects_table_{}", report.report_id))
            .striped(true)
            .max_scroll_height(table_height)
            .min_scrolled_height(table_height)
            .auto_shrink([false, false])
            .column(egui_extras::Column::exact(28.0))
            .column(egui_extras::Column::exact(44.0))
            .column(egui_extras::Column::exact(220.0))
            .column(egui_extras::Column::exact(48.0))
            .column(egui_extras::Column::exact(120.0))
            .column(egui_extras::Column::exact(96.0))
            .column(egui_extras::Column::exact(150.0))
            .column(egui_extras::Column::exact(44.0))
            .column(egui_extras::Column::exact(52.0))
            .column(egui_extras::Column::exact(52.0))
            .column(egui_extras::Column::exact(52.0))
            .column(egui_extras::Column::exact(74.0))
            .column(egui_extras::Column::exact(56.0))
            .column(egui_extras::Column::exact(46.0))
            .column(egui_extras::Column::exact(40.0))
            .header(row_height, |mut header| {
                header.col(|ui| {
                    if ui
                        .checkbox(&mut select_filtered_rows, "")
                        .on_hover_text(
                            "Select or clear all aligned rows currently shown by the active filter/search.",
                        )
                        .changed()
                    {
                        let affected = self.set_rna_read_record_indices_selected(
                            displayed_record_indices.iter().copied(),
                            select_filtered_rows,
                        );
                        self.op_status = if select_filtered_rows {
                            format!(
                                "Selected {affected} aligned read(s) from the current filtered subset"
                            )
                        } else {
                            format!(
                                "Cleared {affected} aligned read(s) from the current filtered subset"
                            )
                        };
                    }
                });
                header.col(|ui| {
                    ui.small("Rank");
                });
                header.col(|ui| {
                    ui.small("Read");
                });
                header.col(|ui| {
                    ui.small("Len");
                });
                header.col(|ui| {
                    ui.small("Phase 1");
                });
                header.col(|ui| {
                    ui.small("Effect");
                });
                header.col(|ui| {
                    ui.small("Phase 2");
                });
                header.col(|ui| {
                    ui.small("Str");
                });
                header.col(|ui| {
                    ui.small("Id%");
                });
                header.col(|ui| {
                    ui.small("Cov%");
                });
                header.col(|ui| {
                    ui.small("Tx%");
                });
                header.col(|ui| {
                    ui.label("FL")
                        .on_hover_text(Self::rna_read_full_length_tooltip());
                });
                header.col(|ui| {
                    ui.small("Score");
                });
                header.col(|ui| {
                    ui.small("Exons");
                });
                header.col(|ui| {
                    ui.small("Jx");
                });
            })
            .body(|body| {
                body.rows(row_height, displayed_rows.len(), |mut table_row| {
                    let row = displayed_rows[table_row.index()];
                    table_row.col(|ui| {
                        let mut include_for_copy = self
                            .rna_seed_selected_record_indices
                            .contains(&row.record_index);
                        if ui.checkbox(&mut include_for_copy, "").changed() {
                            if include_for_copy {
                                self.rna_seed_selected_record_indices.insert(row.record_index);
                            } else {
                                self.rna_seed_selected_record_indices.remove(&row.record_index);
                            }
                        }
                    });
                    table_row.col(|ui| {
                        ui.monospace(row.rank.to_string());
                    });
                    table_row.col(|ui| {
                        let selected =
                            self.rna_seed_highlight_record_index == Some(row.record_index);
                        let response = ui
                            .selectable_label(
                                selected,
                                format!("#{} {}", row.record_index + 1, row.header_id),
                            )
                            .on_hover_text(
                                "Select this aligned read to inspect its phase-1 and phase-2 details below.",
                            );
                        if response.clicked() {
                            self.rna_seed_highlight_record_index = Some(row.record_index);
                        }
                    });
                    table_row.col(|ui| {
                        let read_len = hits_by_record_index
                            .get(&row.record_index)
                            .map(|hit| hit.read_length_bp)
                            .unwrap_or(0);
                        ui.monospace(read_len.to_string());
                    });
                    table_row.col(|ui| {
                        ui.monospace(if row.phase1_primary_transcript_id.trim().is_empty() {
                            "none".to_string()
                        } else {
                            row.phase1_primary_transcript_id.clone()
                        });
                    });
                    table_row.col(|ui| {
                        ui.monospace(Self::rna_read_alignment_effect_label(
                            row.alignment_effect,
                        ));
                    });
                    table_row.col(|ui| {
                        ui.monospace(Self::compact_rna_read_transcript_label(
                            &row.transcript_id,
                            &row.transcript_label,
                        ));
                    });
                    table_row.col(|ui| {
                        ui.monospace(if row.strand.trim().is_empty() {
                            "na".to_string()
                        } else {
                            row.strand.clone()
                        });
                    });
                    table_row.col(|ui| {
                        ui.monospace(format!("{:.1}", row.identity_fraction * 100.0));
                    });
                    table_row.col(|ui| {
                        ui.monospace(format!("{:.1}", row.query_coverage_fraction * 100.0));
                    });
                    table_row.col(|ui| {
                        ui.monospace(format!("{:.1}", row.target_coverage_fraction * 100.0));
                    });
                    table_row.col(|ui| {
                        let full_length_class = GentleEngine::rna_read_full_length_class_label(
                            row.full_length_exact,
                            row.full_length_near,
                            row.full_length_strict,
                        );
                        ui.colored_label(
                            Self::rna_read_full_length_class_color(full_length_class),
                            full_length_class,
                        )
                        .on_hover_text(Self::rna_read_full_length_tooltip());
                    });
                    table_row.col(|ui| {
                        ui.monospace(row.score.to_string());
                    });
                    table_row.col(|ui| {
                        ui.monospace(row.mapped_exon_support.len().to_string());
                    });
                    table_row.col(|ui| {
                        ui.monospace(row.mapped_junction_support.len().to_string());
                    });
                });
            });

        if displayed_rows.is_empty() {
            ui.small(
                egui::RichText::new(
                    "No aligned rows match the current filter/search. Clear the controls or select a broader subset.",
                )
                .color(egui::Color32::from_rgb(180, 83, 9)),
            );
            return;
        }

        ui.separator();
        let Some(selected_record_index) = self.rna_seed_highlight_record_index else {
            ui.small(
                "Select one aligned read row above to inspect its alignment effect in detail.",
            );
            return;
        };
        let Some(selected_row) = displayed_rows
            .iter()
            .copied()
            .find(|row| row.record_index == selected_record_index)
        else {
            ui.small(
                "The currently highlighted read is outside the current filtered subset. Select a row from the filtered table or widen the filter/search.",
            );
            return;
        };
        let Some(selected_hit) = hits_by_record_index.get(&selected_record_index).copied() else {
            ui.small("The highlighted aligned read could not be found in the saved report.");
            return;
        };
        self.render_rna_read_alignment_effect_detail(
            ui,
            &report.report_id,
            view,
            progress,
            report.align_config.min_identity_fraction,
            selected_row,
            selected_hit,
        );
    }

    pub(super) fn render_rna_read_alignment_effect_detail(
        &mut self,
        ui: &mut egui::Ui,
        report_id: &str,
        view: &SplicingExpertView,
        progress: &RnaReadInterpretProgress,
        min_identity_fraction: f64,
        row: &RnaReadAlignmentInspectionRow,
        hit: &RnaReadInterpretationHit,
    ) {
        ui.group(|ui| {
            ui.label(
                egui::RichText::new(format!(
                    "Selected aligned read: #{} {}",
                    row.record_index + 1,
                    row.header_id
                ))
                .strong(),
            );
            ui.horizontal_wrapped(|ui| {
                ui.small(format!(
                    "Effect: {}",
                    Self::rna_read_alignment_effect_label(row.alignment_effect)
                ));
                ui.separator();
                ui.small(format!("Length: {} bp", hit.read_length_bp));
                ui.separator();
                ui.small(format!("Seed pass: {}", row.passed_seed_filter));
                ui.separator();
                ui.small(format!("MSA eligible: {}", row.msa_eligible));
                ui.separator();
                ui.small(format!("Origin: {}", row.origin_class.as_str()));
            });
            ui.horizontal_wrapped(|ui| {
                if ui.button("Copy highlighted FASTA").clicked() {
                    self.copy_rna_report_hits_as_fasta(
                        ui,
                        &[hit],
                        "highlighted aligned read detail",
                    );
                }
                if ui.button("Materialize highlighted").clicked() {
                    self.materialize_highlighted_rna_read_report_hit();
                }
                let alignment_cache_key =
                    Self::rna_read_alignment_detail_cache_key(report_id, row.record_index);
                let alignment_visible = self
                    .rna_read_alignment_detail_visible_key
                    .as_deref()
                    == Some(alignment_cache_key.as_str());
                if ui
                    .button(if alignment_visible {
                        "Hide alignment"
                    } else {
                        "Show alignment"
                    })
                    .on_hover_text(
                        "Inspect the exact phase-2 read-vs-transcript-template pairwise alignment that justified this mapping.",
                    )
                    .clicked()
                {
                    if alignment_visible {
                        self.rna_read_alignment_detail_visible_key = None;
                    } else {
                        self.rna_read_alignment_detail_visible_key =
                            Some(alignment_cache_key.clone());
                    }
                }
                if ui.button("Open interactive dotplot").clicked() {
                    self.open_rna_read_dotplot_workspace(view, hit);
                }
                if ui.button("Export dotplot...").clicked() {
                    self.export_highlighted_rna_read_sequence_dotplot_svg(view, progress);
                }
            });
            ui.separator();
            ui.label(egui::RichText::new("Phase 1").strong());
            ui.small(format!(
                "Primary transcript guess: {}",
                if row.phase1_primary_transcript_id.trim().is_empty() {
                    "none".to_string()
                } else {
                    row.phase1_primary_transcript_id.clone()
                }
            ));
            ui.small(format!(
                "Seed-chain transcript: {}",
                if row.seed_chain_transcript_id.trim().is_empty() {
                    "none".to_string()
                } else {
                    row.seed_chain_transcript_id.clone()
                }
            ));
            ui.small(format!(
                "Exon-path transcript/path: {} | {}",
                if row.exon_path_transcript_id.trim().is_empty() {
                    "none".to_string()
                } else {
                    row.exon_path_transcript_id.clone()
                },
                if row.exon_path.trim().is_empty() {
                    "none".to_string()
                } else {
                    row.exon_path.clone()
                }
            ));
            ui.small(format!(
                "Confirmed exon transitions: {}/{} | selected strand={} | rc_applied={}",
                row.exon_transitions_confirmed,
                row.exon_transitions_total,
                if row.selected_strand.trim().is_empty() {
                    "na".to_string()
                } else {
                    row.selected_strand.clone()
                },
                row.reverse_complement_applied
            ));
            ui.separator();
            ui.label(egui::RichText::new("Phase 2").strong());
            ui.small(format!(
                "Aligned transcript: {}",
                Self::compact_rna_read_transcript_label(&row.transcript_id, &row.transcript_label)
            ));
            ui.small(format!(
                "Mode={} strand={} query={} target={}..{} id={:.1}% cov={:.1}% tx_cov={:.1}% score={} secondary={}",
                row.alignment_mode.as_str(),
                if row.strand.trim().is_empty() {
                    "na".to_string()
                } else {
                    row.strand.clone()
                },
                hit.best_mapping
                    .as_ref()
                    .map(|mapping| {
                        if mapping.query_reverse_complemented {
                            "reverse-complemented"
                        } else {
                            "as stored"
                        }
                    })
                    .unwrap_or("as stored"),
                row.target_start_1based,
                row.target_end_1based,
                row.identity_fraction * 100.0,
                row.query_coverage_fraction * 100.0,
                row.target_coverage_fraction * 100.0,
                row.score,
                row.secondary_mapping_count
            ));
            let aligned_target_bp = row
                .target_end_1based
                .saturating_sub(row.target_start_1based)
                .saturating_add(1);
            ui.small(format!(
                "Transcript span: {}",
                Self::format_rna_read_target_coverage_summary(
                    aligned_target_bp,
                    row.target_length_bp,
                    row.target_coverage_fraction
                )
            ));
            let full_length_class = GentleEngine::rna_read_full_length_class_label(
                row.full_length_exact,
                row.full_length_near,
                row.full_length_strict,
            );
            ui.colored_label(
                Self::rna_read_full_length_class_color(full_length_class),
                format!("Full-length class: {full_length_class}"),
            )
            .on_hover_text(Self::rna_read_full_length_tooltip());
            ui.small(format!(
                "exact (100% template coverage): {} | near (>=95%): {} | strict_end (near + both ends <=15 bp + id >= {:.1}%): {}",
                row.full_length_exact,
                row.full_length_near,
                min_identity_fraction * 100.0,
                row.full_length_strict
            ));
            if let Some(note) = Self::rna_read_fragment_alignment_note(
                row.target_coverage_fraction,
                aligned_target_bp,
                row.target_length_bp,
            ) {
                ui.small(
                    egui::RichText::new(note).color(egui::Color32::from_rgb(180, 83, 9)),
                );
            }
            ui.small(format!(
                "Phase-1 normalization={} | Phase-2 query orientation={}",
                if row.reverse_complement_applied {
                    "reverse-complemented before scoring"
                } else {
                    "kept in input orientation"
                },
                hit.best_mapping
                    .as_ref()
                    .map(|mapping| {
                        if mapping.query_reverse_complemented {
                            "reverse-complemented to fit the transcript template"
                        } else {
                            "stored query already fits the transcript template"
                        }
                    })
                    .unwrap_or("unknown"),
            ));
            egui::CollapsingHeader::new("Phase-2 pairwise alignment")
                .default_open(true)
                .show(ui, |ui| match self.current_rna_read_alignment_display(row.record_index) {
                    Ok(display) => {
                        ui.small(format!(
                            "Exact rust-bio pairwise alignment: mode={} | query={} | aligned columns={} | matches={} mismatches={} insertions={} deletions={}",
                            display.alignment_mode.as_str(),
                            if display.query_reverse_complemented {
                                "reverse-complemented"
                            } else {
                                "as stored"
                            },
                            display.aligned_columns,
                            display.matches,
                            display.mismatches,
                            display.insertions,
                            display.deletions,
                        ));
                        ui.small(format!(
                            "Query span {}..{} of {} bp | template offsets {}..{} of {} bp ({:.1}%) | genomic {}..{}",
                            display.query_start_0based,
                            display.query_end_0based_exclusive,
                            hit.read_length_bp,
                            display.target_start_offset_0based,
                            display.target_end_offset_0based_exclusive,
                            display.target_length_bp,
                            display.target_coverage_fraction * 100.0,
                            display.target_start_1based,
                            display.target_end_1based,
                        ));
                        ui.small(
                            "Legend: `|` exact match, `.` mismatch, blank = insertion/deletion.",
                        );
                        egui::ScrollArea::both()
                            .id_salt(format!(
                                "rna_alignment_effect_pairwise_{}",
                                row.record_index
                            ))
                            .max_height(200.0)
                            .auto_shrink([false, false])
                            .show(ui, |ui| {
                                for (query_chunk, midline_chunk, target_chunk) in
                                    Self::rna_read_alignment_display_chunks(&display, 96)
                                {
                                    let query_label = if display.query_reverse_complemented {
                                        "Q(rc)"
                                    } else {
                                        "Q"
                                    };
                                    ui.monospace(format!("{query_label:<5} {query_chunk}"));
                                    ui.monospace(format!("{:<5} {}", "", midline_chunk));
                                    ui.monospace(format!("{:<5} {}", "T", target_chunk));
                                    ui.add_space(4.0);
                                }
                            });
                    }
                    Err(err) => {
                        ui.small(
                            egui::RichText::new(err).color(egui::Color32::from_rgb(180, 83, 9)),
                        );
                    }
                });
            ui.separator();
            ui.label(egui::RichText::new("Mapped support contribution").strong());
            ui.small(format!(
                "Exons: {}",
                Self::format_rna_read_mapped_exon_support_compact(row)
            ));
            ui.small(format!(
                "Junctions: {}",
                Self::format_rna_read_mapped_junction_support_compact(row)
            ));
            ui.separator();
            let alignment_cache_key = Self::rna_read_alignment_detail_cache_key(report_id, row.record_index);
            if self
                .rna_read_alignment_detail_visible_key
                .as_deref()
                == Some(alignment_cache_key.as_str())
            {
                ui.label(egui::RichText::new("Phase-2 alignment detail").strong());
                match self.saved_rna_read_alignment_detail_for_report_id(report_id, row.record_index)
                {
                    Ok(detail) => {
                        ui.small(format!(
                            "Backend={} mode={} cigar={} aligned={} columns | query {}..{} / {} bp | template offsets {}..{} / {} bp ({:.1}%) | genomic {}..{}",
                            detail.backend.as_str(),
                            detail.alignment_mode.as_str(),
                            detail.cigar,
                            detail.aligned_columns,
                            detail.aligned_query_start_0based,
                            detail.aligned_query_end_0based_exclusive,
                            detail.query_length_bp,
                            detail.aligned_target_start_offset_0based,
                            detail.aligned_target_end_offset_0based_exclusive,
                            detail.target_length_bp,
                            detail.target_coverage_fraction * 100.0,
                            detail.target_start_1based,
                            detail.target_end_1based
                        ));
                        ui.small(format!(
                            "Matches={} mismatches={} insertions={} deletions={} | id={:.1}% cov={:.1}% score={}",
                            detail.matches,
                            detail.mismatches,
                            detail.insertions,
                            detail.deletions,
                            detail.identity_fraction * 100.0,
                            detail.query_coverage_fraction * 100.0,
                            detail.score
                        ));
                        ui.small(
                            "Only the aligned query span contributes to phase-2 confirmation; clipped read tails remain outside the local/semiglobal alignment.",
                        );
                        egui::ScrollArea::both()
                            .id_salt(format!("rna_alignment_effect_detail_{}", row.record_index))
                            .max_height(240.0)
                            .auto_shrink([false, false])
                            .show(ui, |ui| {
                                Self::render_rna_read_pairwise_alignment_detail(ui, &detail);
                            });
                    }
                    Err(message) => {
                        ui.small(
                            egui::RichText::new(message)
                                .color(egui::Color32::from_rgb(180, 83, 9)),
                        );
                    }
                }
            } else {
                ui.small(
                    "Alignment detail is available on demand. Use `Show alignment` to inspect the exact phase-2 read-vs-transcript-template comparison, and `Open interactive dotplot` only when you want the external visual view.",
                );
            }
            ui.separator();
            ui.label(egui::RichText::new("Sequence").strong());
            egui::ScrollArea::both()
                .id_salt(format!(
                    "rna_alignment_effect_sequence_{}",
                    row.record_index
                ))
                .max_height(96.0)
                .auto_shrink([false, false])
                .show(ui, |ui| {
                    ui.monospace(hit.sequence.as_str());
                });
        });
    }

    pub(super) fn render_rna_read_pairwise_alignment_detail(
        ui: &mut egui::Ui,
        detail: &RnaReadPairwiseAlignmentDetail,
    ) {
        const ALIGNMENT_BLOCK_WIDTH: usize = 96;
        let query_bytes = detail.aligned_query.as_bytes();
        let relation_bytes = detail.aligned_relation.as_bytes();
        let target_bytes = detail.aligned_target.as_bytes();
        if query_bytes.is_empty() || target_bytes.is_empty() {
            ui.small("No aligned columns are available for this detail view.");
            return;
        }
        let mut offset = 0usize;
        while offset < query_bytes.len() {
            let end = (offset + ALIGNMENT_BLOCK_WIDTH).min(query_bytes.len());
            let query_chunk = String::from_utf8_lossy(&query_bytes[offset..end]);
            let relation_chunk = String::from_utf8_lossy(&relation_bytes[offset..end]);
            let target_chunk = String::from_utf8_lossy(&target_bytes[offset..end]);
            ui.monospace(format!("Q {offset:>5}  {query_chunk}"));
            ui.monospace(format!("           {relation_chunk}"));
            ui.monospace(format!("T {offset:>5}  {target_chunk}"));
            ui.add_space(4.0);
            offset = end;
        }
    }

    pub(super) fn rna_read_alignment_display_chunks(
        display: &RnaReadAlignmentDisplay,
        chunk_width: usize,
    ) -> Vec<(String, String, String)> {
        let width = chunk_width.max(1);
        let query_chunks = display
            .aligned_query
            .as_bytes()
            .chunks(width)
            .map(|chunk| String::from_utf8_lossy(chunk).to_string())
            .collect::<Vec<_>>();
        let midline_chunks = display
            .aligned_midline
            .as_bytes()
            .chunks(width)
            .map(|chunk| String::from_utf8_lossy(chunk).to_string())
            .collect::<Vec<_>>();
        let target_chunks = display
            .aligned_target
            .as_bytes()
            .chunks(width)
            .map(|chunk| String::from_utf8_lossy(chunk).to_string())
            .collect::<Vec<_>>();
        let chunk_count = query_chunks
            .len()
            .max(midline_chunks.len())
            .max(target_chunks.len());
        let mut rows = Vec::<(String, String, String)>::with_capacity(chunk_count);
        for idx in 0..chunk_count {
            rows.push((
                query_chunks.get(idx).cloned().unwrap_or_default(),
                midline_chunks.get(idx).cloned().unwrap_or_default(),
                target_chunks.get(idx).cloned().unwrap_or_default(),
            ));
        }
        rows
    }

    pub(super) fn collect_thresholded_cdna_exon_support_rows(
        view: &SplicingExpertView,
        isoform_support_rows: &[RnaReadIsoformSupportRow],
        seed_passed_denominator: usize,
    ) -> Vec<RnaReadExonSupportFrequency> {
        let exon_index = view
            .unique_exons
            .iter()
            .enumerate()
            .map(|(idx, exon)| ((exon.start_1based, exon.end_1based), idx))
            .collect::<HashMap<_, _>>();
        let mut support_counts = vec![0usize; view.unique_exons.len()];
        for row in isoform_support_rows {
            if row.reads_seed_passed == 0 {
                continue;
            }
            let Some(transcript) = view
                .transcripts
                .iter()
                .find(|lane| lane.transcript_feature_id == row.transcript_feature_id)
            else {
                continue;
            };
            for exon in &transcript.exons {
                let Some(exon_idx) = exon_index
                    .get(&(exon.start_1based, exon.end_1based))
                    .copied()
                else {
                    continue;
                };
                support_counts[exon_idx] =
                    support_counts[exon_idx].saturating_add(row.reads_seed_passed);
            }
        }
        view.unique_exons
            .iter()
            .enumerate()
            .map(|(idx, exon)| {
                let support_read_count = support_counts[idx];
                let support_fraction = if seed_passed_denominator == 0 {
                    0.0
                } else {
                    support_read_count as f64 / seed_passed_denominator as f64
                };
                RnaReadExonSupportFrequency {
                    start_1based: exon.start_1based,
                    end_1based: exon.end_1based,
                    support_read_count,
                    support_fraction,
                }
            })
            .collect::<Vec<_>>()
    }

    pub(super) fn render_reported_transcript_support_tables(
        &self,
        ui: &mut egui::Ui,
        view: &SplicingExpertView,
    ) {
        ui.small(format!(
            "Reported transcript annotation: transcripts={} | unique exons={} | junctions={} | strand={}",
            view.transcript_count,
            view.unique_exons.len(),
            view.junctions.len(),
            view.strand
        ));

        ui.collapsing("Reported exon support", |ui| {
            if view.unique_exons.is_empty() {
                ui.small("No reported exon rows are available in the current annotation scope.");
                return;
            }
            egui::ScrollArea::vertical()
                .max_height(Self::default_rna_read_support_table_height(ui))
                .min_scrolled_height(Self::default_rna_read_support_table_height(ui))
                .show(ui, |ui| {
                    egui::Grid::new(format!("reported_exon_support_grid_{}", view.seq_id))
                        .striped(true)
                        .num_columns(5)
                        .show(ui, |ui| {
                            ui.small("Exon span");
                            ui.small("Transcripts");
                            ui.small("Transcript %");
                            ui.small("Constitutive");
                            ui.small("Fraction");
                            ui.end_row();
                            for exon in &view.unique_exons {
                                let pct = if view.transcript_count == 0 {
                                    0.0
                                } else {
                                    (exon.support_transcript_count as f64
                                        / view.transcript_count as f64)
                                        * 100.0
                                };
                                ui.monospace(format!("{}..{}", exon.start_1based, exon.end_1based));
                                ui.monospace(exon.support_transcript_count.to_string());
                                ui.monospace(format!("{pct:.2}%"));
                                ui.monospace(if exon.constitutive { "yes" } else { "no" });
                                ui.monospace(format!(
                                    "{:.4}",
                                    if view.transcript_count == 0 {
                                        0.0
                                    } else {
                                        exon.support_transcript_count as f64
                                            / view.transcript_count as f64
                                    }
                                ));
                                ui.end_row();
                            }
                        });
                });
        });

        ui.collapsing("Reported junction support", |ui| {
            if view.junctions.is_empty() {
                ui.small("No reported exon-exon junction rows are available in the current annotation scope.");
                return;
            }
            egui::ScrollArea::vertical()
                .max_height(Self::default_rna_read_support_table_height(ui))
                .min_scrolled_height(Self::default_rna_read_support_table_height(ui))
                .show(ui, |ui| {
                    egui::Grid::new(format!("reported_junction_support_grid_{}", view.seq_id))
                        .striped(true)
                        .num_columns(5)
                        .show(ui, |ui| {
                            ui.small("Donor");
                            ui.small("Acceptor");
                            ui.small("Transcripts");
                            ui.small("Transcript %");
                            ui.small("Fraction");
                            ui.end_row();
                            for junction in &view.junctions {
                                let pct = if view.transcript_count == 0 {
                                    0.0
                                } else {
                                    (junction.support_transcript_count as f64
                                        / view.transcript_count as f64)
                                        * 100.0
                                };
                                ui.monospace(junction.donor_1based.to_string());
                                ui.monospace(junction.acceptor_1based.to_string());
                                ui.monospace(junction.support_transcript_count.to_string());
                                ui.monospace(format!("{pct:.2}%"));
                                ui.monospace(format!(
                                    "{:.4}",
                                    if view.transcript_count == 0 {
                                        0.0
                                    } else {
                                        junction.support_transcript_count as f64
                                            / view.transcript_count as f64
                                    }
                                ));
                                ui.end_row();
                            }
                        });
                });
        });

        ui.collapsing("Reported isoform catalogue", |ui| {
            if view.transcripts.is_empty() {
                ui.small("No reported transcripts are available in the current annotation scope.");
                return;
            }
            egui::ScrollArea::vertical()
                .max_height(Self::default_rna_read_support_table_height(ui))
                .min_scrolled_height(Self::default_rna_read_support_table_height(ui))
                .show(ui, |ui| {
                    egui::Grid::new(format!("reported_isoform_catalogue_grid_{}", view.seq_id))
                        .striped(true)
                        .num_columns(5)
                        .show(ui, |ui| {
                            ui.small("Transcript");
                            ui.small("Strand");
                            ui.small("Exons");
                            ui.small("Expected jx");
                            ui.small("Target");
                            ui.end_row();
                            for transcript in &view.transcripts {
                                let row_color = if transcript.has_target_feature {
                                    egui::Color32::from_rgb(30, 64, 175)
                                } else {
                                    egui::Color32::from_gray(80)
                                };
                                ui.colored_label(
                                    row_color,
                                    format!("{} ({})", transcript.label, transcript.transcript_id),
                                );
                                ui.monospace(transcript.strand.as_str());
                                ui.monospace(transcript.exons.len().to_string());
                                ui.monospace(transcript.exons.len().saturating_sub(1).to_string());
                                ui.monospace(if transcript.has_target_feature {
                                    "yes"
                                } else {
                                    "no"
                                });
                                ui.end_row();
                            }
                        });
                });
            ui.small(
                "These rows are the reported transcript model from the annotation currently loaded into the splicing view.",
            );
        });
    }

    pub(super) fn render_thresholded_cdna_support_tables(
        &self,
        ui: &mut egui::Ui,
        view: &SplicingExpertView,
        progress: &RnaReadInterpretProgress,
    ) {
        let thresholded_exon_rows = Self::collect_thresholded_cdna_exon_support_rows(
            view,
            &progress.isoform_support_rows,
            progress.seed_passed,
        );
        ui.small(
            "Thresholded cDNA uses phase-1 seed-passed reads only. Exon rows are inferred from assigned transcript paths, and junction rows come from confirmed seed-supported transitions.",
        );
        ui.collapsing("Thresholded cDNA exon support", |ui| {
            if thresholded_exon_rows.is_empty() {
                ui.small("No thresholded exon-support rows are available yet.");
                return;
            }
            egui::ScrollArea::vertical()
                .max_height(Self::default_rna_read_support_table_height(ui))
                .min_scrolled_height(Self::default_rna_read_support_table_height(ui))
                .show(ui, |ui| {
                    egui::Grid::new(format!(
                        "thresholded_cdna_exon_support_grid_{}",
                        progress.seq_id
                    ))
                    .striped(true)
                    .num_columns(4)
                    .show(ui, |ui| {
                        ui.small("Exon span");
                        ui.small("Reads");
                        ui.small("Seed-pass %");
                        ui.small("Fraction");
                        ui.end_row();
                        for row in &thresholded_exon_rows {
                            ui.monospace(format!("{}..{}", row.start_1based, row.end_1based));
                            ui.monospace(row.support_read_count.to_string());
                            ui.monospace(format!("{:.2}%", row.support_fraction * 100.0));
                            ui.monospace(format!("{:.4}", row.support_fraction));
                            ui.end_row();
                        }
                    });
                });
        });
        self.render_rna_read_transition_support_table(ui, progress);
        self.render_rna_read_isoform_support_table(ui, progress);
    }

    pub(super) fn render_rna_read_transition_support_table(
        &self,
        ui: &mut egui::Ui,
        progress: &RnaReadInterpretProgress,
    ) {
        let seed_passed_denominator = progress.seed_passed;
        let reads_with_support_pct = if seed_passed_denominator == 0 {
            0.0
        } else {
            (progress.reads_with_transition_support as f64 / seed_passed_denominator as f64) * 100.0
        };
        ui.small(format!(
            "Thresholded cDNA junction support: indexed junction-crossing bits={} | reads with confirmed exon-exon transitions (seed-passed)={}/{} ({:.2}%) | confirmed transitions total={}",
            progress.junction_crossing_seed_bits_indexed,
            progress.reads_with_transition_support,
            seed_passed_denominator,
            reads_with_support_pct,
            progress.transition_confirmations,
        ));
        ui.collapsing("Thresholded cDNA junction support", |ui| {
            if progress.transition_support_rows.is_empty() {
                ui.small("No exon-exon transition catalog rows in current scope.");
                return;
            }
            egui::ScrollArea::vertical()
                .max_height(Self::default_rna_read_support_table_height(ui))
                .min_scrolled_height(Self::default_rna_read_support_table_height(ui))
                .show(ui, |ui| {
                    egui::Grid::new(format!("rna_transition_support_grid_{}", progress.seq_id))
                        .striped(true)
                        .num_columns(5)
                        .show(ui, |ui| {
                            ui.small("Transition");
                            ui.small("From (bp)");
                            ui.small("To (bp)");
                            ui.small("Reads");
                            ui.small("Read %");
                            ui.end_row();
                            for row in &progress.transition_support_rows {
                                let pct = if seed_passed_denominator == 0 {
                                    0.0
                                } else {
                                    (row.support_read_count as f64 / seed_passed_denominator as f64)
                                        * 100.0
                                };
                                ui.monospace(format!(
                                    "E{} -> E{}",
                                    row.from_exon_ordinal, row.to_exon_ordinal
                                ));
                                ui.monospace(format!(
                                    "{}..{}",
                                    row.from_start_1based, row.from_end_1based
                                ));
                                ui.monospace(format!(
                                    "{}..{}",
                                    row.to_start_1based, row.to_end_1based
                                ));
                                ui.monospace(row.support_read_count.to_string());
                                ui.monospace(format!("{pct:.2}%"));
                                ui.end_row();
                            }
                        });
                });
        });
    }

    pub(super) fn render_rna_read_isoform_support_table(
        &self,
        ui: &mut egui::Ui,
        progress: &RnaReadInterpretProgress,
    ) {
        if progress.isoform_support_rows.is_empty() {
            ui.small("No isoform support rows available yet.");
            return;
        }
        let auto_pick = progress
            .isoform_support_rows
            .iter()
            .find(|row| row.reads_seed_passed > 0)
            .unwrap_or(&progress.isoform_support_rows[0]);
        ui.small(format!(
            "Top thresholded cDNA isoform: {} ({}) strand={} | seed-passed={}/{} | transition-coverage={:.1}%",
            auto_pick.transcript_label,
            auto_pick.transcript_id,
            auto_pick.strand,
            auto_pick.reads_seed_passed,
            auto_pick.reads_assigned,
            auto_pick.transition_rows_supported_fraction * 100.0
        ));
        ui.collapsing("Thresholded cDNA isoform ranking", |ui| {
            egui::ScrollArea::vertical()
                .max_height(Self::default_rna_read_support_table_height(ui))
                .min_scrolled_height(Self::default_rna_read_support_table_height(ui))
                .show(ui, |ui| {
                    egui::Grid::new(format!("rna_isoform_support_grid_{}", progress.seq_id))
                        .striped(true)
                        .num_columns(14)
                        .show(ui, |ui| {
                            ui.small("Transcript");
                            ui.small("Strand");
                            ui.small("Exons");
                            ui.small("Expected jx");
                            ui.small("Assigned");
                            ui.small("Seed-pass");
                            ui.small("Jx supported");
                            ui.small("Jx cov%");
                            ui.small("Mean jx frac");
                            ui.small("Mean gap");
                            ui.small("Best score");
                            ui.small("Chain=same");
                            ui.small("Opposite");
                            ui.small("Ambig");
                            ui.end_row();
                            for row in &progress.isoform_support_rows {
                                let is_auto_pick = row.transcript_id == auto_pick.transcript_id;
                                let row_color = if is_auto_pick {
                                    egui::Color32::from_rgb(30, 64, 175)
                                } else {
                                    egui::Color32::from_gray(80)
                                };
                                ui.colored_label(
                                    row_color,
                                    format!("{} ({})", row.transcript_label, row.transcript_id),
                                );
                                ui.monospace(row.strand.as_str());
                                ui.monospace(row.exon_count.to_string());
                                ui.monospace(row.expected_transition_count.to_string());
                                ui.monospace(row.reads_assigned.to_string());
                                ui.monospace(row.reads_seed_passed.to_string());
                                ui.monospace(row.transition_rows_supported.to_string());
                                ui.monospace(format!(
                                    "{:.1}",
                                    row.transition_rows_supported_fraction * 100.0
                                ));
                                ui.monospace(format!(
                                    "{:.2}",
                                    row.mean_confirmed_transition_fraction
                                ));
                                if row.mean_seed_median_gap < 0.0 {
                                    ui.monospace("na");
                                } else {
                                    ui.monospace(format!("{:.2}", row.mean_seed_median_gap));
                                }
                                ui.monospace(format!(
                                    "{:.3}/{:.3}",
                                    row.best_seed_hit_fraction, row.best_weighted_seed_hit_fraction
                                ));
                                ui.monospace(row.reads_chain_same_strand.to_string());
                                ui.monospace(
                                    row.reads_with_opposite_strand_competition.to_string(),
                                );
                                ui.monospace(row.reads_ambiguous_strand_ties.to_string());
                                ui.end_row();
                            }
                        });
                });
            ui.small(
                "Rows aggregate phase-1 thresholded cDNA evidence over one joint run across all transcripts admitted by scope (including reverse strand when selected).",
            );
        });
    }

    pub(super) fn render_rna_read_mapped_aggregate_support_tables(
        &mut self,
        ui: &mut egui::Ui,
        progress: &RnaReadInterpretProgress,
        report: Option<&RnaReadInterpretationReport>,
    ) {
        ui.small(
            "Mapped cDNA uses phase-2 best mappings only; these rows are the evidence to use for exon/junction/isoform interpretation after alignment.",
        );
        if progress.aligned == 0 {
            ui.small(
                egui::RichText::new(
                    "No mapping-derived exon or isoform support is available yet for this run/report.",
                )
                .color(egui::Color32::from_rgb(180, 83, 9)),
            );
            return;
        }

        ui.small(format!(
            "Aligned reads contributing to mapped support: {}",
            progress.aligned
        ));
        let inspection = report.and_then(|report| {
            self.current_saved_rna_read_alignment_inspection(report.hits.len(), None)
                .ok()
        });
        let exon_contributors = inspection
            .as_ref()
            .map(|inspection| Self::collect_rna_read_mapped_exon_contributors(inspection.as_ref()))
            .unwrap_or_default();
        let junction_contributors = inspection
            .as_ref()
            .map(|inspection| {
                Self::collect_rna_read_mapped_junction_contributors(inspection.as_ref())
            })
            .unwrap_or_default();
        let isoform_contributors = inspection
            .as_ref()
            .map(|inspection| {
                Self::collect_rna_read_mapped_isoform_contributors(inspection.as_ref())
            })
            .unwrap_or_default();
        if inspection.is_some() {
            ui.small(
                egui::RichText::new(
                    "Use Audit to jump from an aggregate row back to its exact aligned reads, or Export... to write that contributor subset directly.",
                )
                .color(egui::Color32::from_rgb(100, 116, 139)),
            );
        } else {
            ui.small(
                egui::RichText::new(
                    "Load a saved report to audit each aggregate row back to its contributing aligned reads.",
                )
                .color(egui::Color32::from_rgb(100, 116, 139)),
            );
        }

        ui.collapsing("Mapped cDNA exon support", |ui| {
            if progress.mapped_exon_support_frequencies.is_empty() {
                ui.small("No mapped exon-overlap rows are available.");
                return;
            }
            egui::ScrollArea::vertical()
                .max_height(Self::default_rna_read_support_table_height(ui))
                .min_scrolled_height(Self::default_rna_read_support_table_height(ui))
                .show(ui, |ui| {
                    egui::Grid::new(format!("rna_mapped_exon_support_grid_{}", progress.seq_id))
                        .striped(true)
                        .num_columns(6)
                        .show(ui, |ui| {
                            ui.small("Exon span");
                            ui.small("Reads");
                            ui.small("Aligned %");
                            ui.small("Fraction");
                            ui.small("Audit");
                            ui.small("Export");
                            ui.end_row();
                            for row in &progress.mapped_exon_support_frequencies {
                                let contributors = exon_contributors
                                    .get(&(row.start_1based, row.end_1based))
                                    .cloned()
                                    .unwrap_or_default();
                                ui.monospace(format!("{}..{}", row.start_1based, row.end_1based));
                                ui.monospace(row.support_read_count.to_string());
                                ui.monospace(format!("{:.2}%", row.support_fraction * 100.0));
                                ui.monospace(format!("{:.4}", row.support_fraction));
                                let mut response = ui.add_enabled(
                                    !contributors.is_empty(),
                                    egui::Button::new(format!("Audit ({})", contributors.len())),
                                );
                                if let Some(inspection) = inspection.as_ref() {
                                    let hover = Self::format_rna_read_contributor_hover_text(
                                        inspection,
                                        &contributors,
                                    );
                                    if !hover.is_empty() {
                                        response = response.on_hover_text(hover);
                                    }
                                }
                                if response.clicked() {
                                    self.focus_rna_read_alignment_effect_record_indices(
                                        contributors.clone(),
                                        &format!(
                                            "mapped exon {}..{}",
                                            row.start_1based, row.end_1based
                                        ),
                                    );
                                }
                                ui.add_enabled_ui(!contributors.is_empty(), |ui| {
                                    ui.menu_button("Export...", |ui| {
                                        for export_kind in [
                                            RnaReadSelectedExportKind::Fasta,
                                            RnaReadSelectedExportKind::AlignmentsTsv,
                                            RnaReadSelectedExportKind::ExonPathsTsv,
                                            RnaReadSelectedExportKind::ExonAbundanceTsv,
                                        ] {
                                            if ui.button(export_kind.menu_label()).clicked() {
                                                self.export_rna_read_subset_with_record_indices(
                                                    export_kind,
                                                    contributors.clone(),
                                                    None,
                                                    "No aligned contributor rows are available to export for this mapped exon.",
                                                );
                                                ui.close();
                                            }
                                        }
                                    });
                                });
                                ui.end_row();
                            }
                        });
                });
        });

        ui.collapsing("Mapped cDNA junction support", |ui| {
            if progress.mapped_junction_support_frequencies.is_empty() {
                ui.small("No mapped junction-overlap rows are available.");
                return;
            }
            egui::ScrollArea::vertical()
                .max_height(Self::default_rna_read_support_table_height(ui))
                .min_scrolled_height(Self::default_rna_read_support_table_height(ui))
                .show(ui, |ui| {
                    egui::Grid::new(format!(
                        "rna_mapped_junction_support_grid_{}",
                        progress.seq_id
                    ))
                    .striped(true)
                    .num_columns(7)
                    .show(ui, |ui| {
                        ui.small("Donor");
                        ui.small("Acceptor");
                        ui.small("Reads");
                        ui.small("Aligned %");
                        ui.small("Fraction");
                        ui.small("Audit");
                        ui.small("Export");
                        ui.end_row();
                        for row in &progress.mapped_junction_support_frequencies {
                            let contributors = junction_contributors
                                .get(&(row.donor_1based, row.acceptor_1based))
                                .cloned()
                                .unwrap_or_default();
                            ui.monospace(row.donor_1based.to_string());
                            ui.monospace(row.acceptor_1based.to_string());
                            ui.monospace(row.support_read_count.to_string());
                            ui.monospace(format!("{:.2}%", row.support_fraction * 100.0));
                            ui.monospace(format!("{:.4}", row.support_fraction));
                            let mut response = ui.add_enabled(
                                !contributors.is_empty(),
                                egui::Button::new(format!("Audit ({})", contributors.len())),
                            );
                            if let Some(inspection) = inspection.as_ref() {
                                let hover = Self::format_rna_read_contributor_hover_text(
                                    inspection,
                                    &contributors,
                                );
                                if !hover.is_empty() {
                                    response = response.on_hover_text(hover);
                                }
                            }
                            if response.clicked() {
                                self.focus_rna_read_alignment_effect_record_indices(
                                    contributors.clone(),
                                    &format!(
                                        "mapped junction {}->{}",
                                        row.donor_1based, row.acceptor_1based
                                    ),
                                );
                            }
                            ui.add_enabled_ui(!contributors.is_empty(), |ui| {
                                ui.menu_button("Export...", |ui| {
                                    for export_kind in [
                                        RnaReadSelectedExportKind::Fasta,
                                        RnaReadSelectedExportKind::AlignmentsTsv,
                                        RnaReadSelectedExportKind::ExonPathsTsv,
                                        RnaReadSelectedExportKind::ExonAbundanceTsv,
                                    ] {
                                        if ui.button(export_kind.menu_label()).clicked() {
                                            self.export_rna_read_subset_with_record_indices(
                                                export_kind,
                                                contributors.clone(),
                                                None,
                                                "No aligned contributor rows are available to export for this mapped junction.",
                                            );
                                            ui.close();
                                        }
                                    }
                                });
                            });
                            ui.end_row();
                        }
                    });
                });
        });

        if progress.mapped_isoform_support_rows.is_empty() {
            ui.small("No mapped isoform rows are available yet.");
            return;
        }
        let auto_pick = &progress.mapped_isoform_support_rows[0];
        ui.small(format!(
            "Top mapped cDNA isoform: {} ({}) strand={} | aligned={} | msa-eligible={} | mean identity={:.1}% | mean query coverage={:.1}%",
            auto_pick.transcript_label,
            auto_pick.transcript_id,
            auto_pick.strand,
            auto_pick.aligned_read_count,
            auto_pick.msa_eligible_read_count,
            auto_pick.mean_identity_fraction * 100.0,
            auto_pick.mean_query_coverage_fraction * 100.0
        ));
        ui.collapsing("Mapped cDNA isoform ranking", |ui| {
            egui::ScrollArea::vertical()
                .max_height(Self::default_rna_read_support_table_height(ui))
                .min_scrolled_height(Self::default_rna_read_support_table_height(ui))
                .show(ui, |ui| {
                    egui::Grid::new(format!(
                        "rna_mapped_isoform_support_grid_{}",
                        progress.seq_id
                    ))
                    .striped(true)
                    .num_columns(10)
                    .show(ui, |ui| {
                        ui.small("Transcript");
                        ui.small("Strand");
                        ui.small("Aligned");
                        ui.small("MSA");
                        ui.small("Mean id%");
                        ui.small("Mean cov%");
                        ui.small("Best score");
                        ui.small("Secondary");
                        ui.small("Audit");
                        ui.small("Export");
                        ui.end_row();
                        for row in &progress.mapped_isoform_support_rows {
                            let contributors = isoform_contributors
                                .get(&row.transcript_id)
                                .cloned()
                                .unwrap_or_default();
                            let is_auto_pick = row.transcript_id == auto_pick.transcript_id;
                            let row_color = if is_auto_pick {
                                egui::Color32::from_rgb(22, 101, 52)
                            } else {
                                egui::Color32::from_gray(80)
                            };
                            ui.colored_label(
                                row_color,
                                format!("{} ({})", row.transcript_label, row.transcript_id),
                            );
                            ui.monospace(row.strand.as_str());
                            ui.monospace(row.aligned_read_count.to_string());
                            ui.monospace(row.msa_eligible_read_count.to_string());
                            ui.monospace(format!("{:.1}", row.mean_identity_fraction * 100.0));
                            ui.monospace(format!(
                                "{:.1}",
                                row.mean_query_coverage_fraction * 100.0
                            ));
                            ui.monospace(row.best_alignment_score.to_string());
                            ui.monospace(row.secondary_mapping_total.to_string());
                            let mut response = ui.add_enabled(
                                !contributors.is_empty(),
                                egui::Button::new(format!("Audit ({})", contributors.len())),
                            );
                            if let Some(inspection) = inspection.as_ref() {
                                let hover = Self::format_rna_read_contributor_hover_text(
                                    inspection,
                                    &contributors,
                                );
                                if !hover.is_empty() {
                                    response = response.on_hover_text(hover);
                                }
                            }
                            if response.clicked() {
                                self.focus_rna_read_alignment_effect_record_indices(
                                    contributors.clone(),
                                    &format!("mapped isoform {}", row.transcript_id),
                                );
                            }
                            ui.add_enabled_ui(!contributors.is_empty(), |ui| {
                                ui.menu_button("Export...", |ui| {
                                    for export_kind in [
                                        RnaReadSelectedExportKind::Fasta,
                                        RnaReadSelectedExportKind::AlignmentsTsv,
                                        RnaReadSelectedExportKind::ExonPathsTsv,
                                        RnaReadSelectedExportKind::ExonAbundanceTsv,
                                    ] {
                                        if ui.button(export_kind.menu_label()).clicked() {
                                            self.export_rna_read_subset_with_record_indices(
                                                export_kind,
                                                contributors.clone(),
                                                None,
                                                "No aligned contributor rows are available to export for this mapped isoform.",
                                            );
                                            ui.close();
                                        }
                                    }
                                });
                            });
                            ui.end_row();
                        }
                    });
                });
            ui.small(
                "Rows aggregate best-mapping evidence only; they are intentionally separate from the thresholded cDNA ranking above.",
            );
        });
    }

    pub(super) fn render_rna_read_top_hits_preview(
        &mut self,
        ui: &mut egui::Ui,
        view: &SplicingExpertView,
        progress: &RnaReadInterpretProgress,
        allow_mapping_actions: bool,
    ) -> Option<Option<usize>> {
        let saved_report = self.current_saved_rna_read_report();
        let score_density_seed_filter_override = self
            .current_rna_read_score_density_seed_filter_override(
                self.rna_read_evidence_ui.score_density_variant,
            )
            .ok()
            .flatten();
        let (mut preview_rows, using_saved_report_score_bin, preview_header, preview_note) = if self
            .rna_read_task
            .is_none()
        {
            if let (Some(report), Some(score_bin_index)) = (
                saved_report.as_ref(),
                self.rna_read_alignment_effect_score_bin_index,
            ) {
                let bin_count = report
                    .score_density_bins
                    .len()
                    .max(progress.score_density_bins.len())
                    .max(40);
                let rows = Self::collect_rna_read_top_hit_previews_for_score_bin(
                    report,
                    score_bin_index,
                    bin_count,
                    self.rna_read_evidence_ui.score_density_variant,
                    score_density_seed_filter_override.as_ref(),
                );
                let score_density_bins = GentleEngine::score_density_bins_for_report_with_override(
                    report,
                    self.rna_read_evidence_ui.score_density_variant,
                    score_density_seed_filter_override.as_ref(),
                )
                .0;
                let tested_count = score_density_bins
                    .get(score_bin_index.min(score_density_bins.len().saturating_sub(1)))
                    .copied()
                    .unwrap_or(0);
                if !rows.is_empty() {
                    let row_count_label = Self::format_count_compact_km(rows.len() as u64);
                    (
                        rows,
                        true,
                        format!(
                            "Score-bin preview ({} retained saved-report rows)",
                            row_count_label
                        ),
                        format!(
                            "This table shows every retained saved-report row in {} score bin {} and is sorted by phase-1 Score. Ret.rank remains the saved-report retention rank. Histogram count for this bin: {}.",
                            Self::rna_read_score_density_variant_label(
                                self.rna_read_evidence_ui.score_density_variant
                            ),
                            Self::format_rna_read_score_bin_spec(score_bin_index, bin_count),
                            tested_count
                        ),
                    )
                } else {
                    (
                        Vec::new(),
                        true,
                        format!(
                            "Score-bin preview (0 retained saved-report rows for {})",
                            Self::format_rna_read_score_bin_spec(score_bin_index, bin_count)
                        ),
                        format!(
                            "The {} histogram bin {} contains {} read(s), but none survived into the retained saved-report subset. The saved report currently stores only the retained top {} row(s).",
                            Self::rna_read_score_density_variant_label(
                                self.rna_read_evidence_ui.score_density_variant
                            ),
                            Self::format_rna_read_score_bin_spec(score_bin_index, bin_count),
                            tested_count,
                            report.hits.len()
                        ),
                    )
                }
            } else {
                (
                        progress.top_hits_preview.clone(),
                        false,
                        format!(
                            "Live preview (top {} retained rows from run)",
                            progress.top_hits_preview.len()
                        ),
                        "This is the capped retained-hit preview from the running/saved report. Rows are sorted here by phase-1 Score; Ret.rank remains the saved-report retention rank. Use mapped `Read effects` above for the full aligned-read inspection surface. Id%=phase-2 pairwise alignment identity; Cov%=query coverage.".to_string(),
                    )
            }
        } else {
            (
                    progress.top_hits_preview.clone(),
                    false,
                    format!(
                        "Live preview (top {} retained rows from run)",
                        progress.top_hits_preview.len()
                    ),
                    "This is the capped retained-hit preview while the run is active. Rows are sorted here by phase-1 Score; Ret.rank remains the live retention rank. Id%=phase-2 pairwise alignment identity; Cov%=query coverage.".to_string(),
                )
        };
        Self::sort_rna_read_top_hit_previews_by_phase1_score(&mut preview_rows);
        let visible_record_indices = preview_rows
            .iter()
            .map(|row| row.record_index)
            .collect::<BTreeSet<_>>();
        let rank_by_record_index = saved_report
            .as_ref()
            .map(|report| {
                report
                    .hits
                    .iter()
                    .enumerate()
                    .map(|(idx, hit)| (hit.record_index, idx + 1))
                    .collect::<HashMap<_, _>>()
            })
            .unwrap_or_default();
        let mut next_selection: Option<Option<usize>> = None;
        egui::CollapsingHeader::new(preview_header)
        .default_open(self.rna_read_task.is_some() || using_saved_report_score_bin)
        .show(ui, |ui| {
                ui.horizontal(|ui| {
                    ui.small(&preview_note);
                    if ui
                        .button("Clear highlight")
                        .on_hover_text("Clear the currently highlighted top-hit/read-effect row.")
                        .clicked()
                    {
                        next_selection = Some(None);
                    }
                    ui.menu_button("Selection tools", |ui| {
                        if let Some(report) = saved_report.as_ref() {
                            if ui
                                .button("Select max-score ties")
                                .on_hover_text(
                                    "Select all saved-report rows tied for the highest retained phase-1 seed score.",
                                )
                                .clicked()
                            {
                                self.rna_seed_selected_record_indices =
                                    Self::select_rna_read_report_max_score_record_indices(report)
                                        .into_iter()
                                        .collect::<BTreeSet<_>>();
                                self.op_status = format!(
                                    "Selected {} saved-report read(s) tied at maximal seed score",
                                    self.rna_seed_selected_record_indices.len()
                                );
                                ui.close();
                            }
                            if ui
                                .button("Select rightmost score bin")
                                .on_hover_text(
                                    "Select saved-report rows from the highest non-empty score-density bin under the current density mode.",
                                )
                                .clicked()
                            {
                                self.rna_seed_selected_record_indices =
                                    Self::select_rna_read_report_rightmost_score_bin_record_indices(
                                        &report,
                                        self.rna_read_evidence_ui.score_density_variant,
                                        score_density_seed_filter_override.as_ref(),
                                    )
                                        .into_iter()
                                        .collect::<BTreeSet<_>>();
                                self.op_status = format!(
                                    "Selected {} saved-report read(s) from the rightmost non-empty score bin",
                                    self.rna_seed_selected_record_indices.len()
                                );
                                ui.close();
                            }
                        } else {
                            ui.small("Save/load a report before using saved-report selection helpers.");
                        }
                        if ui
                            .button("Select listed rows")
                            .on_hover_text(
                                "Select every row currently listed in the top-hit preview table.",
                            )
                            .clicked()
                        {
                            self.rna_seed_selected_record_indices = visible_record_indices.clone();
                            ui.close();
                        }
                        if ui
                            .button("Clear selected")
                            .on_hover_text("Clear the current top-hit/read-effect row selection.")
                            .clicked()
                        {
                            self.rna_seed_selected_record_indices.clear();
                            ui.close();
                        }
                    });
                    if allow_mapping_actions
                        && ui
                            .add_enabled(
                                self.rna_read_task.is_none(),
                                egui::Button::new("Evaluate Top Hits (phase-2)"),
                            )
                            .on_hover_text(format!(
                                "Runs AlignRnaReadReport for report '{}' using current align selection '{}', then refreshes top-hit alignment columns and similarity metrics.",
                                self.rna_reads_ui.report_id.trim(),
                                self.rna_reads_ui.align_phase_selection.as_str()
                            ))
                            .clicked()
                    {
                        self.run_splicing_rna_read_alignment_phase(view);
                    }
                    let selected_count = self.rna_seed_selected_record_indices.len();
                    let selected_top_hits_available = selected_count > 0
                        && saved_report.as_ref().map_or_else(
                            || {
                                preview_rows.iter().any(|row| {
                                    self.rna_seed_selected_record_indices
                                        .contains(&row.record_index)
                                })
                            },
                            |report| {
                                report.hits.iter().any(|hit| {
                                    self.rna_seed_selected_record_indices
                                        .contains(&hit.record_index)
                                })
                            },
                        );
                    let highlighted_top_hit_available = saved_report.as_ref().map_or_else(
                        || {
                            preview_rows.iter().any(|row| {
                                self.rna_seed_highlight_record_index == Some(row.record_index)
                            })
                        },
                        |report| self.selected_highlighted_rna_report_hit(report).is_some(),
                    );
                    if allow_mapping_actions
                        && ui
                            .add_enabled(
                                self.rna_read_task.is_none() && selected_top_hits_available,
                                egui::Button::new(format!(
                                    "Evaluate Selected (phase-2) [{selected_count}]"
                                )),
                            )
                            .on_hover_text(
                                "Align only selected top-hit rows by record_index and refresh inline mapping metrics.",
                            )
                            .clicked()
                    {
                        let selected_indices = self
                            .rna_seed_selected_record_indices
                            .iter()
                            .copied()
                            .collect::<Vec<_>>();
                        self.run_splicing_rna_read_alignment_phase_for_selected(
                            view,
                            selected_indices,
                        );
                    }
                    if ui
                        .add_enabled(
                            selected_top_hits_available,
                            egui::Button::new(format!("Copy selected FASTA ({selected_count})")),
                        )
                        .on_hover_text(
                            "Copy FASTA for selected rows that are present in the current preview/report.",
                        )
                        .clicked()
                    {
                        if let Some(report) = saved_report.as_ref() {
                            let hits = Self::selected_rna_report_hits(
                                report,
                                &self.rna_seed_selected_record_indices,
                            );
                            self.copy_rna_report_hits_as_fasta(ui, &hits, "selected report reads");
                        } else {
                            let rows = preview_rows
                                .iter()
                                .filter(|row| {
                                    self.rna_seed_selected_record_indices
                                        .contains(&row.record_index)
                                })
                                .cloned()
                                .collect::<Vec<_>>();
                            self.copy_rna_top_hit_previews_as_fasta(
                                ui,
                                &rows,
                                "selected top reads",
                            );
                        }
                    }
                    if ui
                        .add_enabled(
                            highlighted_top_hit_available,
                            egui::Button::new("Copy highlighted FASTA"),
                        )
                        .on_hover_text(
                            "Copy FASTA for the highlighted row when it is present in the current preview/report.",
                        )
                        .clicked()
                    {
                        if let Some(report) = saved_report.as_ref() {
                            let hits = self
                                .selected_highlighted_rna_report_hit(report)
                                .into_iter()
                                .collect::<Vec<_>>();
                            self.copy_rna_report_hits_as_fasta(
                                ui,
                                &hits,
                                "highlighted report read",
                            );
                        } else {
                            let rows = preview_rows
                                .iter()
                                .find(|row| {
                                    self.rna_seed_highlight_record_index == Some(row.record_index)
                                })
                                .cloned()
                                .into_iter()
                                .collect::<Vec<_>>();
                            self.copy_rna_top_hit_previews_as_fasta(
                                ui,
                                &rows,
                                "highlighted top read",
                            );
                        }
                    }
                    ui.menu_button("Dotplots", |ui| {
                        if ui
                            .button("Export dotplot for highlighted read...")
                            .on_hover_text(
                                "Write an SVG dotplot for the currently highlighted read against its mapped transcript template.",
                            )
                            .clicked()
                        {
                            self.export_highlighted_rna_read_sequence_dotplot_svg(view, progress);
                            ui.close();
                        }
                        if ui
                            .button("Export dotplots for selected reads...")
                            .on_hover_text(
                                "Write SVG dotplots for all selected reads that can be resolved in the current preview/report.",
                            )
                            .clicked()
                        {
                            self.export_selected_rna_read_sequence_dotplot_svgs(view);
                            ui.close();
                        }
                    });
                    ui.small("Use ↑ / ↓ to iterate listed rows.");
                });
                if preview_rows.is_empty() {
                    ui.small(
                        egui::RichText::new(
                            "No retained saved-report rows match the currently selected histogram bin. Clear the score-bin focus or choose a bin with retained rows.",
                        )
                        .color(egui::Color32::from_rgb(180, 83, 9)),
                    );
                    return;
                }
                let hidden_selected_count = self
                    .rna_seed_selected_record_indices
                    .iter()
                    .filter(|record_index| !visible_record_indices.contains(record_index))
                    .count();
                if hidden_selected_count > 0 {
                    ui.small(format!(
                        "{} selected read(s) are outside the currently listed preview rows but remain active for phase-2 alignment, FASTA copy, and dotplot export.",
                        hidden_selected_count,
                    ));
                }
                ui.small(
                    "Rows are sorted by phase-1 Score. Ret.rank is retention rank. Id%=phase-2 alignment identity. Cov%=query coverage.",
                );
                let mut select_listed_rows = self
                    .rna_read_record_indices_all_selected(visible_record_indices.iter().copied());
                let row_height = Self::virtual_rna_read_table_row_height(ui);
                let table_height =
                    Self::default_rna_read_preview_table_height_for_rows(ui, preview_rows.len());
                egui_extras::TableBuilder::new(ui)
                    .id_salt(format!("rna_top_hits_table_{}", progress.seq_id))
                    .striped(true)
                    .max_scroll_height(table_height)
                    .min_scrolled_height(table_height)
                    .auto_shrink([false, false])
                    .column(egui_extras::Column::exact(28.0))
                    .column(egui_extras::Column::exact(44.0))
                    .column(egui_extras::Column::exact(56.0))
                    .column(egui_extras::Column::exact(220.0))
                    .column(egui_extras::Column::exact(58.0))
                    .column(egui_extras::Column::exact(130.0))
                    .column(egui_extras::Column::exact(52.0))
                    .column(egui_extras::Column::exact(52.0))
                    .column(egui_extras::Column::exact(46.0))
                    .column(egui_extras::Column::exact(150.0))
                    .column(egui_extras::Column::exact(48.0))
                    .header(row_height, |mut header| {
                        header.col(|ui| {
                            if ui
                                .checkbox(&mut select_listed_rows, "")
                                .on_hover_text(
                                    "Select or clear all currently listed top-hit rows without affecting hidden selections.",
                                )
                                .changed()
                            {
                                let affected = self.set_rna_read_record_indices_selected(
                                    visible_record_indices.iter().copied(),
                                    select_listed_rows,
                                );
                                self.op_status = if select_listed_rows {
                                    format!("Selected {affected} listed top-hit read(s)")
                                } else {
                                    format!("Cleared {affected} listed top-hit read(s)")
                                };
                            }
                        });
                        header.col(|ui| {
                            ui.small("Run");
                        });
                        header.col(|ui| {
                            ui.small("Ret.rank");
                        });
                        header.col(|ui| {
                            ui.small("Read");
                        });
                        header.col(|ui| {
                            ui.small("Score");
                        });
                        header.col(|ui| {
                            ui.small("Phase 1 tx");
                        });
                        header.col(|ui| {
                            ui.small("Id%");
                        });
                        header.col(|ui| {
                            ui.small("Cov%");
                        });
                        header.col(|ui| {
                            ui.small("Pass");
                        });
                        header.col(|ui| {
                            ui.small("Phase 2");
                        });
                        header.col(|ui| {
                            ui.small("Len");
                        });
                    })
                    .body(|body| {
                        body.rows(row_height, preview_rows.len(), |mut table_row| {
                            let display_idx = table_row.index();
                            let row = &preview_rows[display_idx];
                            let selected =
                                self.rna_seed_highlight_record_index == Some(row.record_index);
                            let rank = rank_by_record_index
                                .get(&row.record_index)
                                .copied()
                                .unwrap_or(display_idx + 1);
                            let phase1_tx = if row.seed_chain_transcript_id.is_empty() {
                                "none".to_string()
                            } else {
                                row.seed_chain_transcript_id.clone()
                            };
                            let phase2_tx = if row.aligned {
                                Self::compact_rna_read_transcript_label(
                                    &row.best_alignment_transcript_id,
                                    &row.best_alignment_transcript_label,
                                )
                            } else {
                                "none".to_string()
                            };
                            let alignment_summary = Self::rna_top_hit_alignment_summary(row);
                            let gap_median = if row.seed_transcript_gap_count == 0 {
                                "na".to_string()
                            } else {
                                format!("{:.2}", row.seed_median_transcript_gap)
                            };
                            let hover_text = format!(
                                "rank={rank} score={:.3} wscore={:.4} wsupport={:.2} gap-med={} gap-n={} chain={:.2}/{} phase1_tx={} class={} oconf={:.2} sconf={:.2} strand={} opp={} ambig={} matched/tested={}/{} pass={} rc={} msa={} align={} len={} seq={}",
                                row.seed_hit_fraction,
                                row.weighted_seed_hit_fraction,
                                row.weighted_matched_kmers,
                                gap_median,
                                row.seed_transcript_gap_count,
                                row.seed_chain_support_fraction,
                                row.seed_chain_support_kmers,
                                phase1_tx,
                                row.origin_class.as_str(),
                                row.origin_confidence,
                                row.strand_confidence,
                                if row.selected_strand.is_empty() {
                                    "na"
                                } else {
                                    row.selected_strand.as_str()
                                },
                                row.competing_opposite_strand,
                                row.ambiguous_strand_tie,
                                row.matched_kmers,
                                row.tested_kmers,
                                row.passed_seed_filter,
                                row.reverse_complement_applied,
                                row.msa_eligible,
                                alignment_summary,
                                row.read_length_bp,
                                row.sequence_preview,
                            );
                            table_row.col(|ui| {
                                let mut include_for_copy = self
                                    .rna_seed_selected_record_indices
                                    .contains(&row.record_index);
                                if ui.checkbox(&mut include_for_copy, "").changed() {
                                    if include_for_copy {
                                        self.rna_seed_selected_record_indices
                                            .insert(row.record_index);
                                    } else {
                                        self.rna_seed_selected_record_indices
                                            .remove(&row.record_index);
                                    }
                                }
                            });
                            table_row.col(|ui| {
                                if ui
                                    .add_enabled(
                                        self.rna_read_task.is_none(),
                                        egui::Button::new("Run"),
                                    )
                                    .on_hover_text(
                                        "Align only this saved-report row by record_index and refresh phase-2 columns.",
                                    )
                                    .clicked()
                                {
                                    self.rna_seed_selected_record_indices
                                        .insert(row.record_index);
                                    next_selection = Some(Some(row.record_index));
                                    self.run_splicing_rna_read_alignment_phase_for_selected(
                                        view,
                                        vec![row.record_index],
                                    );
                                }
                            });
                            table_row.col(|ui| {
                                ui.monospace(rank.to_string());
                            });
                            table_row.col(|ui| {
                                let response = ui
                                    .selectable_label(
                                        selected,
                                        format!("#{} {}", row.record_index + 1, row.header_id),
                                    )
                                    .on_hover_text(hover_text.clone());
                                if response.clicked() {
                                    next_selection = Some(Some(row.record_index));
                                }
                                response.context_menu(|ui| {
                                    if ui.button("Evaluate this read (phase-2)").clicked() {
                                        self.rna_seed_selected_record_indices
                                            .insert(row.record_index);
                                        next_selection = Some(Some(row.record_index));
                                        self.run_splicing_rna_read_alignment_phase_for_selected(
                                            view,
                                            vec![row.record_index],
                                        );
                                        ui.close();
                                    }
                                    if ui.button("Copy FASTA (this read)").clicked() {
                                        self.copy_rna_top_hit_previews_as_fasta(
                                            ui,
                                            std::slice::from_ref(row),
                                            "context-menu top read",
                                        );
                                        ui.close();
                                    }
                                    if ui.button("Copy FASTA (selected reads)").clicked() {
                                        if let Some(report) = saved_report.as_ref() {
                                            let hits = Self::selected_rna_report_hits(
                                                report,
                                                &self.rna_seed_selected_record_indices,
                                            );
                                            self.copy_rna_report_hits_as_fasta(
                                                ui,
                                                &hits,
                                                "context-menu selected report reads",
                                            );
                                        } else {
                                            let rows = preview_rows
                                                .iter()
                                                .filter(|candidate| {
                                                    self.rna_seed_selected_record_indices
                                                        .contains(&candidate.record_index)
                                                })
                                                .cloned()
                                                .collect::<Vec<_>>();
                                            self.copy_rna_top_hit_previews_as_fasta(
                                                ui,
                                                &rows,
                                                "context-menu selected top reads",
                                            );
                                        }
                                        ui.close();
                                    }
                                    if ui.button("Toggle selected").clicked() {
                                        if self
                                            .rna_seed_selected_record_indices
                                            .contains(&row.record_index)
                                        {
                                            self.rna_seed_selected_record_indices
                                                .remove(&row.record_index);
                                        } else {
                                            self.rna_seed_selected_record_indices
                                                .insert(row.record_index);
                                        }
                                        ui.close();
                                    }
                                });
                            });
                            table_row.col(|ui| {
                                ui.monospace(format!("{:.3}", row.seed_hit_fraction));
                            });
                            table_row.col(|ui| {
                                ui.monospace(phase1_tx);
                            });
                            table_row.col(|ui| {
                                ui.monospace(if row.aligned {
                                    format!("{:.1}", row.best_alignment_identity_fraction * 100.0)
                                } else {
                                    "na".to_string()
                                });
                            });
                            table_row.col(|ui| {
                                ui.monospace(if row.aligned {
                                    format!(
                                        "{:.1}",
                                        row.best_alignment_query_coverage_fraction * 100.0
                                    )
                                } else {
                                    "na".to_string()
                                });
                            });
                            table_row.col(|ui| {
                                ui.monospace(if row.passed_seed_filter {
                                    "yes".to_string()
                                } else {
                                    "no".to_string()
                                });
                            });
                            table_row.col(|ui| {
                                ui.monospace(phase2_tx);
                            });
                            table_row.col(|ui| {
                                ui.monospace(row.read_length_bp.to_string());
                            });
                        });
                    });
                let arrow_down = ui.input(|i| i.key_pressed(egui::Key::ArrowDown));
                let arrow_up = ui.input(|i| i.key_pressed(egui::Key::ArrowUp));
                if (arrow_down || arrow_up) && !preview_rows.is_empty() {
                    let current_pos = self
                        .rna_seed_highlight_record_index
                        .and_then(|record_index| {
                            preview_rows
                                .iter()
                                .position(|row| row.record_index == record_index)
                        });
                    let len = preview_rows.len();
                    let next_pos = if arrow_down {
                        current_pos.map_or(0, |idx| (idx + 1) % len)
                    } else {
                        current_pos.map_or(len.saturating_sub(1), |idx| {
                            if idx == 0 { len - 1 } else { idx - 1 }
                        })
                    };
                    next_selection = Some(Some(preview_rows[next_pos].record_index));
                }
                let copy_shortcut_pressed = ui.input(|i| i.key_pressed(egui::Key::C) && i.modifiers.command);
                if copy_shortcut_pressed && !ui.ctx().egui_wants_keyboard_input() {
                    if let Some(report) = saved_report.as_ref() {
                        let hits = Self::selected_rna_report_hits(
                            report,
                            &self.rna_seed_selected_record_indices,
                        );
                        if !hits.is_empty() {
                            self.copy_rna_report_hits_as_fasta(
                                ui,
                                &hits,
                                "Ctrl/Cmd+C selected report reads",
                            );
                        } else if let Some(highlighted) =
                            self.selected_highlighted_rna_report_hit(&report)
                        {
                            self.copy_rna_report_hits_as_fasta(
                                ui,
                                std::slice::from_ref(&highlighted),
                                "Ctrl/Cmd+C highlighted report read",
                            );
                        } else {
                            self.op_status =
                                "Ctrl/Cmd+C: select or highlight a top read first".to_string();
                        }
                    } else if let Some(highlighted) = preview_rows
                        .iter()
                        .find(|row| self.rna_seed_highlight_record_index == Some(row.record_index))
                        .cloned()
                    {
                        self.copy_rna_top_hit_previews_as_fasta(
                            ui,
                            std::slice::from_ref(&highlighted),
                            "Ctrl/Cmd+C highlighted top read",
                        );
                    } else {
                        self.op_status =
                            "Ctrl/Cmd+C: select or highlight a top read first".to_string();
                    }
                }
            });
        next_selection
    }
}
