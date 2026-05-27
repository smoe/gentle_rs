//! Main lineage workspace rendering helpers.
//!
//! This module is a move-only extraction from `app.rs`: it keeps the
//! project overview strip, lineage graph/table workspace, and lineage modal
//! helpers close to `GENtleApp` while reducing the top-level app monolith.

use super::*;

fn lineage_graph_canvas_width(base_width: f32, graph_zoom: f32, viewport_width: f32) -> f32 {
    (base_width * graph_zoom + 280.0 * graph_zoom).max(viewport_width.max(1.0))
}

impl GENtleApp {
    pub(super) fn render_lineage_node_rename_window(
        &mut self,
        ctx: &egui::Context,
        persist_workspace_after_frame: &mut bool,
    ) {
        let Some(target_node_id) = self.lineage_node_rename_target_id.clone() else {
            return;
        };
        let mut open = true;
        let mut apply_clicked = false;
        let spec = crate::egui_compat::ModalWindowSpec::new(
            "Rename Leaf Node",
            egui::Id::new(("rename_leaf_node_modal", target_node_id.clone())),
        )
        .default_size(Vec2::new(420.0, 120.0));
        crate::egui_compat::show_modal_window(ctx, &spec, &mut open, |ui| {
            ui.monospace(format!("node: {}", target_node_id));
            ui.horizontal(|ui| {
                ui.label("new_name");
                ui.text_edit_singleline(&mut self.lineage_node_rename_text)
                    .on_hover_text("Display name shown in lineage graph/table");
            });
            ui.horizontal(|ui| {
                if ui
                    .button("Apply")
                    .on_hover_text("Apply leaf-node rename")
                    .clicked()
                {
                    apply_clicked = true;
                }
                if ui
                    .button("Cancel")
                    .on_hover_text("Close rename window without changes")
                    .clicked()
                {
                    self.lineage_node_rename_target_id = None;
                    self.lineage_node_rename_text.clear();
                }
            });
        });
        if !open {
            self.lineage_node_rename_target_id = None;
            self.lineage_node_rename_text.clear();
            return;
        }
        if apply_clicked {
            let rename_text = self.lineage_node_rename_text.clone();
            self.lineage_group_status =
                match self.rename_leaf_lineage_node(&target_node_id, &rename_text) {
                    Ok(status) => {
                        *persist_workspace_after_frame = true;
                        self.lineage_node_rename_target_id = None;
                        self.lineage_node_rename_text.clear();
                        status
                    }
                    Err(err) => err,
                };
        }
    }

    pub(super) fn render_lineage_node_remove_confirm_window(
        &mut self,
        ctx: &egui::Context,
        persist_workspace_after_frame: &mut bool,
    ) {
        let Some(target_node_id) = self.lineage_node_remove_target_id.clone() else {
            return;
        };
        let seq_id_hint = {
            let engine = self.engine.read().unwrap();
            let state = engine.state();
            state
                .lineage
                .nodes
                .get(&target_node_id)
                .map(|node| node.seq_id.clone())
                .unwrap_or_else(|| "-".to_string())
        };
        let mut open = true;
        let mut remove_clicked = false;
        let spec = crate::egui_compat::ModalWindowSpec::new(
            "Remove Leaf Node",
            egui::Id::new(("remove_leaf_node_modal", target_node_id.clone())),
        )
        .default_size(Vec2::new(460.0, 130.0));
        crate::egui_compat::show_modal_window(ctx, &spec, &mut open, |ui| {
            ui.colored_label(
                egui::Color32::from_rgb(150, 20, 20),
                "This action removes the sequence/node from the current project.",
            );
            ui.monospace(format!("node: {}", target_node_id));
            ui.monospace(format!("sequence: {}", seq_id_hint));
            ui.horizontal(|ui| {
                if ui
                    .button("Remove")
                    .on_hover_text("Confirm removal of this leaf node")
                    .clicked()
                {
                    remove_clicked = true;
                }
                if ui
                    .button("Cancel")
                    .on_hover_text("Cancel removal")
                    .clicked()
                {
                    self.lineage_node_remove_target_id = None;
                }
            });
        });
        if !open {
            self.lineage_node_remove_target_id = None;
            return;
        }
        if remove_clicked {
            self.lineage_group_status = match self.remove_leaf_lineage_node(&target_node_id) {
                Ok(status) => {
                    *persist_workspace_after_frame = true;
                    self.lineage_node_remove_target_id = None;
                    status
                }
                Err(err) => err,
            };
        }
    }

    pub(super) fn normalized_lineage_panel_heights(
        total_available_height: f32,
        preferred_split_fraction: f32,
    ) -> (f32, f32, f32) {
        let min_graph_height = LINEAGE_MAIN_TOP_PANEL_MIN_HEIGHT;
        let min_container_height = 120.0;
        let min_total_height = min_graph_height + min_container_height;
        let total_height = if total_available_height.is_finite() {
            total_available_height.clamp(min_total_height, 4000.0)
        } else {
            min_total_height
        };
        let clamped_split = preferred_split_fraction.clamp(0.2, 0.9);
        let graph_height = (total_height * clamped_split)
            .clamp(min_graph_height, total_height - min_container_height);
        let container_height = (total_height - graph_height)
            .clamp(min_container_height, total_height - min_graph_height);
        let realized_split = (graph_height / total_height.max(1.0)).clamp(0.2, 0.9);
        (graph_height, container_height, realized_split)
    }

    pub(super) fn normalized_lineage_container_subpanel_heights(
        total_available_height: f32,
        preferred_split_fraction: f32,
    ) -> (f32, f32, f32) {
        let min_containers_height = 80.0;
        let min_arrangements_height = 80.0;
        let min_total_height = min_containers_height + min_arrangements_height;
        let total_height = if total_available_height.is_finite() {
            total_available_height.clamp(min_total_height, 2400.0)
        } else {
            min_total_height
        };
        let clamped_split = preferred_split_fraction.clamp(0.2, 0.8);
        let containers_height = (total_height * clamped_split).clamp(
            min_containers_height,
            total_height - min_arrangements_height,
        );
        let arrangements_height = (total_height - containers_height).clamp(
            min_arrangements_height,
            total_height - min_containers_height,
        );
        let realized_split = (containers_height / total_height.max(1.0)).clamp(0.2, 0.8);
        (containers_height, arrangements_height, realized_split)
    }

    pub(super) fn infer_lineage_analysis_kind_from_row(
        row: &LineageRow,
    ) -> Option<LineageAnalysisKind> {
        if let Some(kind) = row.analysis_kind {
            return Some(kind);
        }
        if row.node_id.starts_with("analysis:dotplot:")
            || row.origin.eq_ignore_ascii_case("dotplot")
        {
            return Some(LineageAnalysisKind::Dotplot);
        }
        if row.node_id.starts_with("analysis:flex:")
            || row.node_id.starts_with("analysis:flexibility:")
            || row.origin.eq_ignore_ascii_case("flexibilitytrack")
        {
            return Some(LineageAnalysisKind::FlexibilityTrack);
        }
        if row.node_id.starts_with("analysis:rna_reads:")
            || row.node_id.starts_with("analysis:rna_read:")
            || row.origin.eq_ignore_ascii_case("rnareadinterpretation")
        {
            return Some(LineageAnalysisKind::RnaReadInterpretation);
        }
        if row.node_id.starts_with("analysis:primer:")
            || row.origin.eq_ignore_ascii_case("primerdesign")
        {
            return Some(LineageAnalysisKind::PrimerDesign);
        }
        if row.node_id.starts_with("analysis:qpcr:")
            || row.origin.eq_ignore_ascii_case("qpcrdesign")
        {
            return Some(LineageAnalysisKind::QpcrDesign);
        }
        if row.node_id.starts_with("analysis:restriction_cloning_pcr:")
            || row
                .origin
                .eq_ignore_ascii_case("restrictioncloningpcrhandoff")
        {
            return Some(LineageAnalysisKind::RestrictionCloningPcrHandoff);
        }
        if row.node_id.starts_with("analysis:protein_derive:")
            || row.node_id.starts_with("analysis:protein_derivation:")
            || row.origin.eq_ignore_ascii_case("proteinderivation")
        {
            return Some(LineageAnalysisKind::ProteinDerivation);
        }
        if row.node_id.starts_with("analysis:reverse_translate:")
            || row.node_id.starts_with("analysis:reverse_translation:")
            || row.origin.eq_ignore_ascii_case("reversetranslation")
        {
            return Some(LineageAnalysisKind::ReverseTranslation);
        }
        if row.node_id.starts_with("analysis:construct_reasoning:")
            || row.node_id.starts_with("analysis:reasoning:")
            || row.origin.eq_ignore_ascii_case("constructreasoning")
        {
            return Some(LineageAnalysisKind::ConstructReasoning);
        }
        if row.node_id.starts_with("analysis:uniprot:")
            || row.node_id.starts_with("analysis:uniprot_projection:")
            || row.origin.eq_ignore_ascii_case("uniprotprojection")
        {
            return Some(LineageAnalysisKind::UniprotProjection);
        }
        if row.node_id.starts_with("analysis:seq_confirm:")
            || row.node_id.starts_with("analysis:sequencing_confirmation:")
            || row.origin.eq_ignore_ascii_case("sequencingconfirmation")
        {
            return Some(LineageAnalysisKind::SequencingConfirmation);
        }
        None
    }

    pub(super) fn infer_lineage_analysis_artifact_id_from_row(row: &LineageRow) -> Option<String> {
        if let Some(artifact_id) = row
            .analysis_artifact_id
            .as_deref()
            .map(str::trim)
            .filter(|id| !id.is_empty())
        {
            return Some(artifact_id.to_string());
        }
        if let Some(rest) = row.node_id.strip_prefix("analysis:dotplot:") {
            let id = rest.trim();
            if !id.is_empty() {
                return Some(id.to_string());
            }
        }
        if let Some(rest) = row.node_id.strip_prefix("analysis:flex:") {
            let id = rest.trim();
            if !id.is_empty() {
                return Some(id.to_string());
            }
        }
        if let Some(rest) = row.node_id.strip_prefix("analysis:flexibility:") {
            let id = rest.trim();
            if !id.is_empty() {
                return Some(id.to_string());
            }
        }
        if let Some(rest) = row.node_id.strip_prefix("analysis:rna_reads:") {
            let id = rest.trim();
            if !id.is_empty() {
                return Some(id.to_string());
            }
        }
        if let Some(rest) = row.node_id.strip_prefix("analysis:rna_read:") {
            let id = rest.trim();
            if !id.is_empty() {
                return Some(id.to_string());
            }
        }
        if let Some(rest) = row.node_id.strip_prefix("analysis:primer:") {
            let id = rest.trim();
            if !id.is_empty() {
                return Some(id.to_string());
            }
        }
        if let Some(rest) = row.node_id.strip_prefix("analysis:qpcr:") {
            let id = rest.trim();
            if !id.is_empty() {
                return Some(id.to_string());
            }
        }
        if let Some(rest) = row
            .node_id
            .strip_prefix("analysis:restriction_cloning_pcr:")
        {
            let id = rest.trim();
            if !id.is_empty() {
                return Some(id.to_string());
            }
        }
        if let Some(rest) = row.node_id.strip_prefix("analysis:protein_derive:") {
            let id = rest.trim();
            if !id.is_empty() {
                return Some(id.to_string());
            }
        }
        if let Some(rest) = row.node_id.strip_prefix("analysis:protein_derivation:") {
            let id = rest.trim();
            if !id.is_empty() {
                return Some(id.to_string());
            }
        }
        if let Some(rest) = row.node_id.strip_prefix("analysis:reverse_translate:") {
            let id = rest.trim();
            if !id.is_empty() {
                return Some(id.to_string());
            }
        }
        if let Some(rest) = row.node_id.strip_prefix("analysis:reverse_translation:") {
            let id = rest.trim();
            if !id.is_empty() {
                return Some(id.to_string());
            }
        }
        if let Some(rest) = row.node_id.strip_prefix("analysis:construct_reasoning:") {
            let id = rest.trim();
            if !id.is_empty() {
                return Some(id.to_string());
            }
        }
        if let Some(rest) = row.node_id.strip_prefix("analysis:reasoning:") {
            let id = rest.trim();
            if !id.is_empty() {
                return Some(id.to_string());
            }
        }
        if let Some(rest) = row.node_id.strip_prefix("analysis:uniprot:") {
            let id = rest.trim();
            if !id.is_empty() {
                return Some(id.to_string());
            }
        }
        if let Some(rest) = row.node_id.strip_prefix("analysis:uniprot_projection:") {
            let id = rest.trim();
            if !id.is_empty() {
                return Some(id.to_string());
            }
        }
        if let Some(rest) = row.node_id.strip_prefix("analysis:seq_confirm:") {
            let id = rest.trim();
            if !id.is_empty() {
                return Some(id.to_string());
            }
        }
        if let Some(rest) = row
            .node_id
            .strip_prefix("analysis:sequencing_confirmation:")
        {
            let id = rest.trim();
            if !id.is_empty() {
                return Some(id.to_string());
            }
        }
        let display_name = row.display_name.trim();
        if !display_name.is_empty() {
            return Some(display_name.to_string());
        }
        None
    }

    pub(super) fn lineage_analysis_open_payload(
        row: &LineageRow,
    ) -> Option<(LineageAnalysisKind, String, String)> {
        let kind = Self::infer_lineage_analysis_kind_from_row(row)?;
        let artifact_id = Self::infer_lineage_analysis_artifact_id_from_row(row)?;
        Some((kind, row.seq_id.clone(), artifact_id))
    }

    pub(super) fn project_overview_metrics(
        rows: &[LineageRow],
        container_count: usize,
        arrangement_count: usize,
    ) -> [ProjectOverviewMetric; 5] {
        let sequence_count = rows
            .iter()
            .filter(|row| row.kind == LineageNodeKind::Sequence && row.pool_size <= 1)
            .count();
        let pool_count = rows
            .iter()
            .filter(|row| row.kind == LineageNodeKind::Sequence && row.pool_size > 1)
            .count();
        let analysis_count = rows
            .iter()
            .filter(|row| row.kind == LineageNodeKind::Analysis)
            .count();
        [
            ProjectOverviewMetric {
                label: "sequences",
                count: sequence_count,
                target: ProjectOverviewTarget::Lineage,
                hover: "Focus the lineage graph/table sequence nodes",
            },
            ProjectOverviewMetric {
                label: "pools",
                count: pool_count,
                target: ProjectOverviewTarget::Lineage,
                hover: "Focus the lineage graph/table pool nodes",
            },
            ProjectOverviewMetric {
                label: "analyses",
                count: analysis_count,
                target: ProjectOverviewTarget::Lineage,
                hover: "Focus the lineage graph/table analysis artifacts",
            },
            ProjectOverviewMetric {
                label: "containers",
                count: container_count,
                target: ProjectOverviewTarget::Containers,
                hover: "Focus the container section below the lineage view",
            },
            ProjectOverviewMetric {
                label: "arrangements",
                count: arrangement_count,
                target: ProjectOverviewTarget::Arrangements,
                hover: "Focus the arrangement section below the lineage view",
            },
        ]
    }

    pub(super) fn focus_project_overview_target(&mut self, target: ProjectOverviewTarget) {
        match target {
            ProjectOverviewTarget::Lineage => {
                self.lineage_graph_view = true;
                self.lineage_main_split_fraction = 0.78;
            }
            ProjectOverviewTarget::Containers => {
                self.lineage_main_split_fraction = 0.36;
                self.lineage_container_arrangement_split_fraction = 0.72;
            }
            ProjectOverviewTarget::Arrangements => {
                self.lineage_main_split_fraction = 0.36;
                self.lineage_container_arrangement_split_fraction = 0.28;
            }
        }
        self.lineage_graph_area_height = 420.0;
        self.lineage_container_area_height = 420.0;
    }

    pub(super) fn render_project_overview_strip(&mut self, ui: &mut Ui) -> bool {
        let metrics = Self::project_overview_metrics(
            &self.lineage_rows,
            self.lineage_containers.len(),
            self.lineage_arrangements.len(),
        );
        let mut changed = false;
        ui.horizontal_wrapped(|ui| {
            ui.label(egui::RichText::new("Overview").strong());
            for metric in metrics {
                let response = ui
                    .button(format!("{} {}", metric.count, metric.label))
                    .on_hover_text(metric.hover);
                if response.clicked() {
                    self.focus_project_overview_target(metric.target);
                    self.app_status = format!("Focused project overview: {}", metric.label);
                    changed = true;
                }
            }
        });
        changed
    }

    pub(super) fn render_main_lineage(&mut self, ui: &mut Ui) {
        self.refresh_lineage_cache_if_needed();

        ui.heading("Lineage Graph");
        ui.label(format!("Project: {}", self.current_project_name()));
        ui.label("Project-level sequence lineage (branch and merge aware)");
        if self.render_project_overview_strip(ui) {
            self.persist_lineage_graph_workspace_to_state();
        }
        ui.horizontal_wrapped(|ui| {
            let table = self.track_hover_status(
                ui.selectable_label(!self.lineage_graph_view, "Table")
                    .on_hover_text("Show lineage as table"),
                "Lineage > View Table",
            );
            if table.clicked() {
                self.lineage_graph_view = false;
                self.persist_lineage_graph_workspace_to_state();
            }
            let graph = self.track_hover_status(
                ui.selectable_label(self.lineage_graph_view, "Graph")
                    .on_hover_text("Show lineage as node graph"),
                "Lineage > View Graph",
            );
            if graph.clicked() {
                self.lineage_graph_view = true;
                self.persist_lineage_graph_workspace_to_state();
            }
        });
        ui.separator();

        let mut open_seq: Option<String> = None;
        let mut open_pool: Option<(String, Vec<String>)> = None;
        let mut open_lineage_analysis: Option<(LineageAnalysisKind, String, String)> = None;
        let mut open_lineage_retrieval: Option<LineageRetrievalDescriptor> = None;
        let mut open_lane_containers: Option<Vec<String>> = None;
        let mut export_container_gel: Option<(String, Vec<String>)> = None;
        let mut set_container_exclusivity: Option<(String, bool)> = None;
        let mut export_arrangement_gel: Option<(String, String)> = None;
        let mut open_arrangement_gel_preview: Option<String> = None;
        let mut export_lineage_dotplot_svg: Option<(String, String, Option<String>)> = None;
        let mut request_export_lineage_svg = false;
        let mut reopen_gibson_from_operation: Option<String> = None;
        let mut reopen_pcr_from_operation: Option<String> = None;
        let mut select_candidate_from: Option<String> = None;
        let mut graph_zoom = self.lineage_graph_zoom.clamp(0.35, 4.0);
        let mut graph_area_height = self
            .lineage_graph_area_height
            .clamp(LINEAGE_MAIN_TOP_PANEL_MIN_HEIGHT, 2400.0);
        let mut container_area_height = self.lineage_container_area_height.clamp(120.0, 1600.0);
        let mut container_arrangement_split_fraction = self
            .lineage_container_arrangement_split_fraction
            .clamp(0.2, 0.8);
        let mut main_split_fraction = self.lineage_main_split_fraction.clamp(0.2, 0.9);
        let mut graph_scroll_offset = Vec2::new(
            if self.lineage_graph_scroll_offset.x.is_finite() {
                self.lineage_graph_scroll_offset.x.max(0.0)
            } else {
                0.0
            },
            if self.lineage_graph_scroll_offset.y.is_finite() {
                self.lineage_graph_scroll_offset.y.max(0.0)
            } else {
                0.0
            },
        );
        let mut graph_compact_labels = self.lineage_graph_compact_labels;
        let mut persist_workspace_after_frame = false;
        let mut graph_rows = self.lineage_rows.clone();
        let mut graph_edges = self.lineage_edges.clone();
        let mut graph_op_label_by_id = self.lineage_op_label_by_id.clone();
        let seq_node_by_seq_id: HashMap<String, String> = self
            .lineage_rows
            .iter()
            .filter(|row| row.kind == LineageNodeKind::Sequence)
            .map(|row| (row.seq_id.clone(), row.node_id.clone()))
            .collect();
        let container_members_by_id: HashMap<String, Vec<String>> = self
            .lineage_containers
            .iter()
            .map(|row| (row.container_id.clone(), row.members.clone()))
            .collect();
        for arrangement in &self.lineage_arrangements {
            let arrangement_node_id = format!("arr:{}", arrangement.arrangement_id);
            let arrangement_edge_op_id = format!(
                "{}::arrangement:{}",
                arrangement.created_by_op, arrangement.arrangement_id
            );
            let mut source_node_ids: Vec<String> = vec![];
            let mut seen_sources: HashSet<String> = HashSet::new();
            for container_id in &arrangement.lane_container_ids {
                if let Some(members) = container_members_by_id.get(container_id) {
                    for seq_id in members {
                        if let Some(source_node_id) = seq_node_by_seq_id.get(seq_id).cloned()
                            && seen_sources.insert(source_node_id.clone())
                        {
                            source_node_ids.push(source_node_id.clone());
                            graph_edges.push((
                                source_node_id,
                                arrangement_node_id.clone(),
                                arrangement_edge_op_id.clone(),
                            ));
                        }
                    }
                }
            }
            graph_op_label_by_id
                .entry(arrangement_edge_op_id)
                .or_insert_with(|| "Arrange serial lanes".to_string());
            graph_rows.push(LineageRow {
                kind: LineageNodeKind::Arrangement,
                node_id: arrangement_node_id,
                seq_id: arrangement.arrangement_id.clone(),
                display_name: if arrangement.name.trim().is_empty() {
                    arrangement.arrangement_id.clone()
                } else {
                    arrangement.name.clone()
                },
                origin: "Arrangement".to_string(),
                created_by_op: arrangement.created_by_op.clone(),
                created_at: arrangement.created_at,
                parents: source_node_ids,
                length: 0,
                circular: false,
                pool_size: 1,
                pool_members: vec![],
                arrangement_id: Some(arrangement.arrangement_id.clone()),
                arrangement_mode: Some(arrangement.mode.clone()),
                lane_container_ids: arrangement.lane_container_ids.clone(),
                ladders: arrangement.ladders.clone(),
                genome_anchor_summary: None,
                genome_anchor_display: None,
                is_full_genome_sequence: false,
                retrieval_descriptor: None,
                analysis_kind: None,
                analysis_artifact_id: None,
                analysis_reference_seq_id: None,
                analysis_mode: None,
                analysis_status: None,
                analysis_point_count: None,
                analysis_bin_count: None,
                analysis_read_count: None,
                analysis_trace_count: None,
                analysis_target_count: None,
                analysis_variant_count: None,
                macro_instance_id: None,
                macro_routine_id: None,
                macro_template_name: None,
                macro_status: None,
                macro_status_message: None,
                macro_op_ids: vec![],
                macro_inputs: vec![],
                macro_outputs: vec![],
            });
        }
        graph_rows.sort_by(|a, b| {
            a.created_at
                .cmp(&b.created_at)
                .then(a.node_id.cmp(&b.node_id))
        });

        let valid_lineage_node_ids: HashSet<String> =
            graph_rows.iter().map(|row| row.node_id.clone()).collect();
        self.lineage_group_marked_nodes
            .retain(|node_id| valid_lineage_node_ids.contains(node_id));
        let mut sanitized_groups =
            Self::sanitize_lineage_node_groups(&self.lineage_node_groups, &valid_lineage_node_ids);
        if sanitized_groups != self.lineage_node_groups {
            self.lineage_node_groups = sanitized_groups.clone();
            if self
                .lineage_group_form_editing_id
                .as_ref()
                .is_some_and(|group_id| {
                    !self
                        .lineage_node_groups
                        .iter()
                        .any(|group| &group.group_id == group_id)
                })
            {
                self.lineage_group_form_editing_id = None;
            }
            persist_workspace_after_frame = true;
        }

        let mut splitter_dragging = false;
        egui::ScrollArea::both()
            .id_salt("lineage_main_scroll")
            .auto_shrink([false, false])
            .scroll_source(egui::containers::scroll_area::ScrollSource {
                drag: false,
                ..Default::default()
            })
            .show(ui, |ui| {
                scroll_input_policy::apply_scrollarea_keyboard_navigation(
                    ui,
                    scroll_input_policy::DEFAULT_SCROLLAREA_KEYBOARD_STEP,
                );
                ui.collapsing("Advanced: Node Groups", |ui| {
                    ui.small(
                        "Group lineage nodes into disjoint sets. Members are shown indented in table view and can collapse to the representative node in graph view.",
                    );
                    if self.lineage_group_marked_nodes.is_empty() {
                        ui.small("Marked nodes: none");
                    } else {
                        let mut marked = self
                            .lineage_group_marked_nodes
                            .iter()
                            .cloned()
                            .collect::<Vec<_>>();
                        marked.sort();
                        ui.horizontal_wrapped(|ui| {
                            ui.small(format!("Marked nodes ({})", marked.len()));
                            ui.monospace(marked.join(", "));
                            if ui
                                .small_button("Clear marked")
                                .on_hover_text("Clear marked candidate nodes for grouping")
                                .clicked()
                            {
                                self.lineage_group_marked_nodes.clear();
                                self.lineage_group_status = "Cleared marked nodes".to_string();
                            }
                        });
                    }
                    if !self.lineage_group_status.trim().is_empty() {
                        ui.small(self.lineage_group_status.clone());
                    }
                    if self.lineage_node_groups.is_empty() {
                        ui.small("No node groups yet.");
                    } else {
                        egui::Grid::new("lineage_node_groups_grid")
                            .striped(true)
                            .min_col_width(48.0)
                            .show(ui, |ui| {
                                ui.strong("Collapse");
                                ui.strong("Group");
                                ui.strong("Representative");
                                ui.strong("Members");
                                ui.strong("Actions");
                                ui.end_row();
                                let mut edit_group_id: Option<String> = None;
                                let mut delete_group_id: Option<String> = None;
                                for group in &mut self.lineage_node_groups {
                                    if ui.checkbox(&mut group.collapsed, "").changed() {
                                        persist_workspace_after_frame = true;
                                    }
                                    ui.label(format!("{} ({})", group.label, group.group_id));
                                    ui.monospace(group.representative_node_id.clone());
                                    ui.monospace(group.member_node_ids.join(", "));
                                    ui.horizontal(|ui| {
                                        if ui
                                            .small_button("Edit")
                                            .on_hover_text("Edit this node group")
                                            .clicked()
                                        {
                                            edit_group_id = Some(group.group_id.clone());
                                        }
                                        if ui
                                            .small_button("Delete")
                                            .on_hover_text("Delete this node group")
                                            .clicked()
                                        {
                                            delete_group_id = Some(group.group_id.clone());
                                        }
                                    });
                                    ui.end_row();
                                }
                                if let Some(group_id) = edit_group_id
                                    && let Some(group) = self
                                        .lineage_node_groups
                                        .iter()
                                        .find(|group| group.group_id == group_id)
                                    {
                                        self.lineage_group_form_editing_id =
                                            Some(group.group_id.clone());
                                        self.lineage_group_form_label = group.label.clone();
                                        self.lineage_group_form_representative =
                                            group.representative_node_id.clone();
                                        self.lineage_group_form_members =
                                            group.member_node_ids.join(", ");
                                        self.lineage_group_status =
                                            format!("Editing group '{}'", group.label);
                                    }
                                if let Some(group_id) = delete_group_id {
                                    let before = self.lineage_node_groups.len();
                                    self.lineage_node_groups
                                        .retain(|group| group.group_id != group_id);
                                    if self.lineage_node_groups.len() != before {
                                        if self
                                            .lineage_group_form_editing_id
                                            .as_ref()
                                            .is_some_and(|editing_id| editing_id == &group_id)
                                        {
                                            self.lineage_group_form_editing_id = None;
                                        }
                                        self.lineage_group_status =
                                            format!("Deleted group '{group_id}'");
                                        persist_workspace_after_frame = true;
                                    }
                                }
                            });
                    }
                    ui.separator();
                    let editing = self.lineage_group_form_editing_id.is_some();
                    ui.label(if editing {
                        "Edit node group"
                    } else {
                        "Create node group"
                    });
                    ui.horizontal(|ui| {
                        ui.label("Label");
                        ui.text_edit_singleline(&mut self.lineage_group_form_label);
                    });
                    ui.horizontal(|ui| {
                        ui.label("Representative");
                        ui.text_edit_singleline(&mut self.lineage_group_form_representative);
                        if ui
                            .small_button("Use selected")
                            .on_hover_text("Use currently selected graph node as representative")
                            .clicked()
                            && let Some(selected_node_id) = &self.lineage_graph_selected_node_id {
                                self.lineage_group_form_representative = selected_node_id.clone();
                            }
                    });
                    ui.horizontal(|ui| {
                        ui.label("Members");
                        ui.text_edit_singleline(&mut self.lineage_group_form_members);
                    });
                    ui.small(
                        "Members: comma/space/newline separated node ids (representative excluded).",
                    );
                    ui.horizontal(|ui| {
                        let save_label = if editing { "Apply Group" } else { "Add Group" };
                        if ui.button(save_label).clicked() {
                            match self.apply_lineage_group_form(&valid_lineage_node_ids) {
                                Ok(status) => {
                                    self.lineage_group_status = status;
                                    persist_workspace_after_frame = true;
                                }
                                Err(err) => {
                                    self.lineage_group_status = err;
                                }
                            }
                        }
                        if ui.button("Clear Form").clicked() {
                            self.lineage_group_form_editing_id = None;
                            self.lineage_group_form_label.clear();
                            self.lineage_group_form_members.clear();
                            self.lineage_group_form_representative.clear();
                            self.lineage_group_status.clear();
                            self.lineage_group_marked_nodes.clear();
                        }
                    });
                });
                ui.separator();

                sanitized_groups = Self::sanitize_lineage_node_groups(
                    &self.lineage_node_groups,
                    &valid_lineage_node_ids,
                );
                if sanitized_groups != self.lineage_node_groups {
                    self.lineage_node_groups = sanitized_groups.clone();
                    persist_workspace_after_frame = true;
                }
                let lineage_groups = self.lineage_node_groups.clone();
                let (group_by_id, node_to_group_id) = Self::lineage_node_group_maps(&lineage_groups);
                let (projected_graph_rows, projected_graph_edges) =
                    Self::project_lineage_graph_by_groups(&graph_rows, &graph_edges, &lineage_groups);
                let collapsed_group_badges = Self::lineage_collapsed_group_hidden_op_badges(
                    &graph_edges,
                    &lineage_groups,
                    &graph_op_label_by_id,
                );
                let table_entries = Self::build_lineage_table_entries(&graph_rows, &lineage_groups);
                let (
                    graph_view_rows,
                    graph_view_edges,
                    graph_view_op_label_by_id,
                ) = Self::project_lineage_graph_operation_hubs(
                    &projected_graph_rows,
                    &projected_graph_edges,
                    &graph_op_label_by_id,
                    &self.lineage_reopenable_gibson_op_ids,
                );
                let leaf_node_ids = Self::lineage_leaf_node_ids(&graph_rows, &graph_edges);
                let visible_node_ids: HashSet<String> = graph_view_rows
                    .iter()
                    .map(|row| row.node_id.clone())
                    .collect();
                if let Some(selected_node_id) = self.lineage_graph_selected_node_id.clone()
                    && !visible_node_ids.contains(&selected_node_id) {
                        let remapped = node_to_group_id
                            .get(&selected_node_id)
                            .and_then(|group_id| group_by_id.get(group_id))
                            .map(|group| group.representative_node_id.clone());
                        self.lineage_graph_selected_node_id = remapped;
                    }
                // Keep graph + containers filling the lineage viewport proportionally.
                let visible_height = ui.clip_rect().height();
                let panel_budget_height = if visible_height.is_finite() && visible_height > 0.0 {
                    ui.available_height().min(visible_height)
                } else {
                    ui.available_height()
                };
                let (target_graph_height, target_container_height, normalized_split) =
                    Self::normalized_lineage_panel_heights(panel_budget_height, main_split_fraction);
                graph_area_height = target_graph_height;
                container_area_height = target_container_height;
                main_split_fraction = normalized_split;

                if self.lineage_graph_view {
            if self.lineage_graph_offsets_synced_stamp != self.lineage_cache_stamp {
                let offset_count_before = self.lineage_graph_node_offsets.len();
                let active_node_ids = graph_rows
                    .iter()
                    .map(|row| row.node_id.clone())
                    .collect::<std::collections::HashSet<_>>();
                self.lineage_graph_node_offsets
                    .retain(|node_id, _| active_node_ids.contains(node_id));
                if self
                    .lineage_graph_drag_origin
                    .as_ref()
                    .is_some_and(|(node_id, _)| !active_node_ids.contains(node_id))
                {
                    self.lineage_graph_drag_origin = None;
                }
                if self
                    .lineage_graph_pan_origin
                    .as_ref()
                    .is_some_and(|_| active_node_ids.is_empty())
                {
                    self.lineage_graph_pan_origin = None;
                }
                if self
                    .lineage_graph_selected_node_id
                    .as_ref()
                    .is_some_and(|node_id| !active_node_ids.contains(node_id))
                {
                    self.lineage_graph_selected_node_id = None;
                }
                if self.lineage_graph_node_offsets.len() != offset_count_before {
                    persist_workspace_after_frame = true;
                }
                self.lineage_graph_offsets_synced_stamp = self.lineage_cache_stamp;
            }
            let mut request_fit_zoom = false;
            let mut request_fit_origin = false;
            let science_palette = SciencePalette::default();
            ui.horizontal_wrapped(|ui| {
                ui.label("Legend:");
                ui.colored_label(science_palette.sequence_node, "● single sequence");
                ui.colored_label(science_palette.pool_node, "◆ pool");
                ui.colored_label(science_palette.macro_node, "▭ macro run");
                ui.colored_label(science_palette.arrangement_node, "▭ arrangement");
                ui.colored_label(science_palette.retrieval_node, "▣ retrieval pattern");
                ui.colored_label(science_palette.node_group, "◻ node group");
            });
            ui.horizontal_wrapped(|ui| {
                let zoom_out_resp = self.track_hover_status(
                    ui.button("−").on_hover_text("Zoom out"),
                    "Lineage Graph > Zoom Out",
                );
                if zoom_out_resp.clicked() {
                    graph_zoom = (graph_zoom / 1.15).clamp(0.35, 4.0);
                }
                let zoom_in_resp = self.track_hover_status(
                    ui.button("+").on_hover_text("Zoom in"),
                    "Lineage Graph > Zoom In",
                );
                if zoom_in_resp.clicked() {
                    graph_zoom = (graph_zoom * 1.15).clamp(0.35, 4.0);
                }
                let zoom_reset_resp = self.track_hover_status(
                    ui.button("Reset").on_hover_text("Reset zoom"),
                    "Lineage Graph > Zoom Reset",
                );
                if zoom_reset_resp.clicked() {
                    graph_zoom = 1.0;
                }
                let fit_resp = self.track_hover_status(
                    ui.button("Fit")
                        .on_hover_text("Fit full graph content into current graph area"),
                    "Lineage Graph > Fit",
                );
                if fit_resp.clicked() {
                    request_fit_zoom = true;
                    request_fit_origin = true;
                }
                let reset_layout_resp = self.track_hover_status(
                    ui.button("Reset Layout")
                        .on_hover_text("Reset manually moved node positions"),
                    "Lineage Graph > Reset Layout",
                );
                if reset_layout_resp.clicked() {
                    self.lineage_graph_node_offsets.clear();
                    self.lineage_graph_drag_origin = None;
                    self.lineage_graph_pan_origin = None;
                    persist_workspace_after_frame = true;
                }
                if ui
                    .checkbox(&mut graph_compact_labels, "Compact labels")
                    .on_hover_text(
                        "Reduce label density in crowded graphs by simplifying operation and node labels",
                    )
                    .changed()
                {
                    persist_workspace_after_frame = true;
                }
                ui.add(
                    egui::Slider::new(&mut graph_zoom, 0.35..=4.0)
                        .logarithmic(true)
                        .text("Zoom"),
                );
                ui.label(format!("{:.0}%", graph_zoom * 100.0));
            });
            ui.separator();
            let rows = &graph_view_rows;
            let lineage_edges = &graph_view_edges;
            let op_label_by_id = &graph_view_op_label_by_id;
            let graph_panel_width = ui.available_width().max(1.0);
            let graph_canvas_viewport_width = (graph_panel_width - 16.0).max(1.0);
            let graph_panel_height = graph_area_height.max(LINEAGE_MAIN_TOP_PANEL_MIN_HEIGHT);
            ui.allocate_ui_with_layout(
                egui::vec2(graph_panel_width, graph_panel_height),
                egui::Layout::top_down(egui::Align::Min),
                |ui| {
                    ui.set_max_width(graph_panel_width);
                    let dark_mode = ui.visuals().dark_mode;
                    theme::canvas_frame(dark_mode).show(ui, |ui| {
                        ui.set_max_width(graph_canvas_viewport_width);
                        ui.small(
                            "Shift+scroll zooms (Cmd/Ctrl+scroll still works). Option+drag pans; Space+drag pans on background. Drag nodes to reposition.",
                        );
                        let graph_scroll_output = egui::ScrollArea::both()
                        .id_salt("lineage_graph_scroll")
                        .auto_shrink([false, false])
                        .scroll_offset(graph_scroll_offset)
                        .scroll_source(egui::containers::scroll_area::ScrollSource {
                            drag: false,
                            ..Default::default()
                        })
                        .max_height(ui.available_height())
                        .show(ui, |ui| {
                            let (layout_by_node, layer_count, max_nodes_in_layer) =
                                Self::compute_lineage_dag_layout(rows, lineage_edges);
                            let base_width = (layer_count.max(1) as f32) * 220.0 + 220.0;
                            let base_height = (max_nodes_in_layer.max(1) as f32) * 110.0 + 300.0;
                            if request_fit_zoom {
                                let available_height = ui.available_height();
                                let fit_x =
                                    ((graph_canvas_viewport_width - 24.0).max(120.0)
                                        / base_width.max(1.0))
                                    .max(0.01);
                                let fit_y =
                                    ((available_height - 24.0).max(120.0)
                                        / base_height.max(1.0))
                                        .max(0.01);
                                graph_zoom = fit_x.min(fit_y).clamp(0.35, 4.0);
                                request_fit_zoom = false;
                            }
                            if request_fit_origin {
                                graph_scroll_offset = Vec2::ZERO;
                                self.lineage_graph_pan_origin = None;
                                request_fit_origin = false;
                            }
                            let width = lineage_graph_canvas_width(
                                base_width,
                                graph_zoom,
                                graph_canvas_viewport_width,
                            );
                            let height = base_height * graph_zoom;
                            let (resp, painter) = ui.allocate_painter(
                                Vec2::new(width, height),
                                egui::Sense::click_and_drag(),
                            );
                            let mut graph_wheel_intent = WheelIntent::None;
                            if resp.hovered() {
                                let (wheel_delta, modifiers) =
                                    ui.input(|i| (i.smooth_scroll_delta, i.modifiers));
                                graph_wheel_intent = scroll_input_policy::lineage_graph_wheel_intent(
                                    wheel_delta,
                                    modifiers,
                                );
                                if let WheelIntent::Zoom { direction, amount } = graph_wheel_intent
                                {
                                    let signed_amount = if direction == ZoomDirection::In {
                                        amount
                                    } else {
                                        -amount
                                    };
                                    let factor = (1.0 + signed_amount * 0.0015).clamp(0.8, 1.25);
                                    graph_zoom = (graph_zoom * factor).clamp(0.35, 4.0);
                                }
                                if !ui.ctx().egui_wants_keyboard_input() {
                                    let keyboard_pan_delta = ui.input(|i| {
                                        scroll_input_policy::canvas_keyboard_pan_delta(
                                            i,
                                            scroll_input_policy::DEFAULT_CANVAS_PAN_STEP,
                                        )
                                    });
                                    if keyboard_pan_delta != Vec2::ZERO {
                                        graph_scroll_offset.x =
                                            (graph_scroll_offset.x + keyboard_pan_delta.x).max(0.0);
                                        graph_scroll_offset.y =
                                            (graph_scroll_offset.y + keyboard_pan_delta.y).max(0.0);
                                    }
                                    let zoom_steps =
                                        ui.input(scroll_input_policy::canvas_keyboard_zoom_steps);
                                    if zoom_steps != 0 {
                                        for _ in 0..zoom_steps.unsigned_abs() {
                                            graph_zoom = if zoom_steps > 0 {
                                                (graph_zoom * 1.12).clamp(0.35, 4.0)
                                            } else {
                                                (graph_zoom / 1.12).clamp(0.35, 4.0)
                                            };
                                        }
                                    }
                                }
                            }
                            let dense_graph = rows.len() >= 22 || lineage_edges.len() >= 30;
                            let simplify_labels = graph_compact_labels && dense_graph;
                            let edge_label_stride = if simplify_labels {
                                if lineage_edges.len() > 180 {
                                    4
                                } else if lineage_edges.len() > 120 {
                                    3
                                } else if lineage_edges.len() > 60 {
                                    2
                                } else {
                                    1
                                }
                            } else {
                                1
                            };
                            let rect = resp.rect;
                            let edge_stroke_width = (1.0 * graph_zoom).clamp(1.0, 2.5);
                            let node_id_font_size = (10.0 * graph_zoom).clamp(8.0, 15.0);
                            let op_font_size = if simplify_labels {
                                (9.0 * graph_zoom).clamp(8.0, 14.0)
                            } else {
                                (10.0 * graph_zoom).clamp(9.0, 16.0)
                            };
                            let name_font_size = if simplify_labels {
                                (11.0 * graph_zoom).clamp(9.0, 15.0)
                            } else {
                                (12.0 * graph_zoom).clamp(10.0, 18.0)
                            };
                            let details_font_size = (10.0 * graph_zoom).clamp(9.0, 15.0);
                            let node_radius = 16.0 * graph_zoom;
                            let mut pos_by_node: HashMap<String, Pos2> = HashMap::new();
                            for (fallback_rank, row) in rows.iter().enumerate() {
                                let (layer, rank) = layout_by_node
                                    .get(&row.node_id)
                                    .copied()
                                    .unwrap_or((0, fallback_rank));
                                let manual_offset = self
                                    .lineage_graph_node_offsets
                                    .get(&row.node_id)
                                    .copied()
                                    .unwrap_or(Vec2::ZERO);
                                let x = rect.left()
                                    + (120.0 + layer as f32 * 220.0) * graph_zoom
                                    + manual_offset.x;
                                let y = rect.top()
                                    + (120.0 + rank as f32 * 110.0) * graph_zoom
                                    + manual_offset.y;
                                pos_by_node.insert(row.node_id.clone(), Pos2::new(x, y));
                            }
                            for group in &lineage_groups {
                                let mut group_nodes: Vec<String> =
                                    Vec::with_capacity(group.member_node_ids.len() + 1);
                                group_nodes.push(group.representative_node_id.clone());
                                group_nodes.extend(group.member_node_ids.clone());
                                let mut bounds: Option<egui::Rect> = None;
                                let mut visible_nodes = 0usize;
                                for node_id in group_nodes {
                                    let Some(pos) = pos_by_node.get(&node_id).copied() else {
                                        continue;
                                    };
                                    visible_nodes += 1;
                                    let glyph_rect = egui::Rect::from_center_size(
                                        pos,
                                        Vec2::new(52.0 * graph_zoom, 42.0 * graph_zoom),
                                    );
                                    bounds = Some(if let Some(existing) = bounds {
                                        existing.union(glyph_rect)
                                    } else {
                                        glyph_rect
                                    });
                                }
                                let Some(mut bounds) = bounds else {
                                    continue;
                                };
                                bounds = bounds.expand(24.0 * graph_zoom);
                                let group_color = Self::lineage_group_color(&group.group_id);
                                painter.rect_filled(
                                    bounds,
                                    10.0 * graph_zoom,
                                    egui::Color32::from_rgba_premultiplied(
                                        group_color.r(),
                                        group_color.g(),
                                        group_color.b(),
                                        18,
                                    ),
                                );
                                painter.rect_stroke(
                                    bounds,
                                    10.0 * graph_zoom,
                                    egui::Stroke::new((1.4 * graph_zoom).clamp(1.0, 2.2), group_color),
                                    egui::StrokeKind::Inside,
                                );
                                let mut label = group.label.clone();
                                if group.collapsed {
                                    label.push_str(&format!(" (collapsed, {} hidden)", group.member_node_ids.len()));
                                } else {
                                    label.push_str(&format!(" ({visible_nodes} node(s))"));
                                }
                                painter.text(
                                    Pos2::new(bounds.left() + 8.0 * graph_zoom, bounds.top() - 2.0 * graph_zoom),
                                    egui::Align2::LEFT_BOTTOM,
                                    label,
                                    egui::FontId::proportional((10.0 * graph_zoom).clamp(9.0, 14.0)),
                                    group_color,
                                );
                            }
                            let pointer = resp.interact_pointer_pos();
                            let is_node_hit = |row: &LineageRow, pointer: Pos2| -> bool {
                                let Some(pos) = pos_by_node.get(&row.node_id).copied() else {
                                    return false;
                                };
                                let glyph_hit = match row.kind {
                                    LineageNodeKind::Arrangement => egui::Rect::from_center_size(
                                        pos,
                                        Vec2::new(56.0 * graph_zoom, 28.0 * graph_zoom),
                                    )
                                    .contains(pointer),
                                    LineageNodeKind::Macro => egui::Rect::from_center_size(
                                        pos,
                                        Vec2::new(64.0 * graph_zoom, 30.0 * graph_zoom),
                                    )
                                    .contains(pointer),
                                    LineageNodeKind::Analysis => egui::Rect::from_center_size(
                                        pos,
                                        if Self::is_lineage_operation_hub(row) {
                                            Vec2::new(74.0 * graph_zoom, 30.0 * graph_zoom)
                                        } else {
                                            Vec2::new(62.0 * graph_zoom, 28.0 * graph_zoom)
                                        },
                                    )
                                    .contains(pointer),
                                    LineageNodeKind::Sequence if row.pool_size > 1 => {
                                        egui::Rect::from_center_size(
                                            pos,
                                            Vec2::new(34.0 * graph_zoom, 28.0 * graph_zoom),
                                        )
                                        .contains(pointer)
                                    }
                                    LineageNodeKind::Sequence => pointer.distance(pos) <= 19.0 * graph_zoom,
                                };
                                if glyph_hit {
                                    return true;
                                }
                                let display_name = if simplify_labels {
                                    Self::compact_lineage_node_label(&row.display_name, 26)
                                } else {
                                    row.display_name.clone()
                                };
                                let name_galley = painter.layout_no_wrap(
                                    display_name,
                                    egui::FontId::proportional(name_font_size),
                                    egui::Color32::BLACK,
                                );
                                let name_anchor = pos + Vec2::new(22.0 * graph_zoom, -4.0 * graph_zoom);
                                let name_rect = egui::Rect::from_min_size(
                                    Pos2::new(name_anchor.x, name_anchor.y - name_galley.size().y),
                                    name_galley.size(),
                                )
                                .expand(3.0 * graph_zoom);
                                if name_rect.contains(pointer) {
                                    return true;
                                }
                                if !simplify_labels {
                                    let detail_text = match row.kind {
                                        LineageNodeKind::Arrangement => {
                                            let ladders = if row.ladders.is_empty() {
                                                "auto".to_string()
                                            } else {
                                                row.ladders.join(" + ")
                                            };
                                            format!(
                                                "{} lane(s) | mode={} | ladders={}",
                                                row.lane_container_ids.len(),
                                                row.arrangement_mode
                                                    .as_deref()
                                                    .unwrap_or("Serial")
                                                    .to_lowercase(),
                                                ladders
                                            )
                                        }
                                        LineageNodeKind::Macro => format!(
                                            "ops={} | template={} | status={}",
                                            row.macro_op_ids.len(),
                                            row.macro_template_name
                                                .as_deref()
                                                .unwrap_or("-"),
                                            row.macro_status.as_deref().unwrap_or("ok")
                                        ),
                                        LineageNodeKind::Analysis => {
                                            if Self::is_lineage_operation_hub(row) {
                                                format!("operation={}", row.created_by_op)
                                            } else {
                                                let artifact_id = row
                                                    .analysis_artifact_id
                                                    .as_deref()
                                                    .unwrap_or(&row.display_name);
                                                let mode = row.analysis_mode.as_deref().unwrap_or("-");
                                                match row.analysis_kind {
                                                    Some(LineageAnalysisKind::Dotplot) => format!(
                                                        "{} | mode={} | points={}",
                                                        artifact_id,
                                                        mode,
                                                        row.analysis_point_count.unwrap_or(0)
                                                    ),
                                                    Some(LineageAnalysisKind::FlexibilityTrack) => {
                                                        format!(
                                                            "{} | model={} | bins={}",
                                                            artifact_id,
                                                            mode,
                                                            row.analysis_bin_count.unwrap_or(0)
                                                        )
                                                    }
                                                    Some(LineageAnalysisKind::RnaReadInterpretation) => {
                                                        format!(
                                                            "{} | profile={} | status={} | reads={} | seed_passed={} | aligned={}",
                                                            artifact_id,
                                                            mode,
                                                            row.analysis_status
                                                                .as_deref()
                                                                .unwrap_or("-"),
                                                            row.analysis_read_count.unwrap_or(0),
                                                            row.analysis_point_count.unwrap_or(0),
                                                            row.analysis_variant_count.unwrap_or(0)
                                                        )
                                                    }
                                                    Some(LineageAnalysisKind::PrimerDesign) => {
                                                        format!(
                                                            "{} | backend={} | pairs={}",
                                                            artifact_id,
                                                            mode,
                                                            row.analysis_target_count.unwrap_or(0)
                                                        )
                                                    }
                                                    Some(LineageAnalysisKind::QpcrDesign) => {
                                                        format!(
                                                            "{} | backend={} | assays={}",
                                                            artifact_id,
                                                            mode,
                                                            row.analysis_target_count.unwrap_or(0)
                                                        )
                                                    }
                                                    Some(LineageAnalysisKind::RestrictionCloningPcrHandoff) => {
                                                        format!(
                                                            "{} | mode={} | status={} | vector={}",
                                                            artifact_id,
                                                            mode,
                                                            row.analysis_status
                                                                .as_deref()
                                                                .unwrap_or("-"),
                                                            row.analysis_reference_seq_id
                                                                .as_deref()
                                                                .unwrap_or("-")
                                                        )
                                                    }
                                                    Some(LineageAnalysisKind::ProteinDerivation) => {
                                                        format!(
                                                            "{} | mode={} | proteins={}",
                                                            artifact_id,
                                                            mode,
                                                            row.analysis_target_count.unwrap_or(0)
                                                        )
                                                    }
                                                    Some(LineageAnalysisKind::ReverseTranslation) => {
                                                        format!(
                                                            "{} | profile={} | table={} | bp={} | aa={}",
                                                            artifact_id,
                                                            mode,
                                                            row.analysis_point_count.unwrap_or(0),
                                                            row.analysis_target_count.unwrap_or(0),
                                                            row.analysis_variant_count.unwrap_or(0)
                                                        )
                                                    }
                                                    Some(LineageAnalysisKind::ConstructReasoning) => {
                                                        format!(
                                                            "{} | objective={} | evidence={} | decisions={} | candidates={}",
                                                            artifact_id,
                                                            mode,
                                                            row.analysis_point_count.unwrap_or(0),
                                                            row.analysis_variant_count.unwrap_or(0),
                                                            row.analysis_target_count.unwrap_or(0)
                                                        )
                                                    }
                                                    Some(LineageAnalysisKind::UniprotProjection) => {
                                                        format!(
                                                            "{} | entry={} | transcripts={}",
                                                            artifact_id,
                                                            mode,
                                                            row.analysis_target_count.unwrap_or(0)
                                                        )
                                                    }
                                                    Some(LineageAnalysisKind::SequencingConfirmation) => {
                                                        format!(
                                                            "{} | status={} | reads={} | traces={} | targets={} | variants={}",
                                                            artifact_id,
                                                            row.analysis_status
                                                                .as_deref()
                                                                .unwrap_or("-"),
                                                            row.analysis_read_count.unwrap_or(0),
                                                            row.analysis_trace_count.unwrap_or(0),
                                                            row.analysis_target_count.unwrap_or(0),
                                                            row.analysis_variant_count.unwrap_or(0)
                                                        )
                                                    }
                                                    None => artifact_id.to_string(),
                                                }
                                            }
                                        }
                                        LineageNodeKind::Sequence if row.pool_size > 1 => {
                                            format!(
                                                "{} ({} bp) | pool={}",
                                                row.seq_id, row.length, row.pool_size
                                            )
                                        }
                                        LineageNodeKind::Sequence => {
                                            format!("{} ({} bp)", row.seq_id, row.length)
                                        }
                                    };
                                    let details_galley = painter.layout_no_wrap(
                                        detail_text,
                                        egui::FontId::proportional(details_font_size),
                                        egui::Color32::BLACK,
                                    );
                                    let details_anchor =
                                        pos + Vec2::new(22.0 * graph_zoom, 10.0 * graph_zoom);
                                    let details_rect =
                                        egui::Rect::from_min_size(details_anchor, details_galley.size())
                                            .expand(3.0 * graph_zoom);
                                    if details_rect.contains(pointer) {
                                        return true;
                                    }
                                }
                                false
                            };
                            let hovered_node_id = pointer.and_then(|pointer| {
                                rows.iter()
                                    .filter_map(|row| {
                                        if !is_node_hit(row, pointer) {
                                            return None;
                                        }
                                        pos_by_node
                                            .get(&row.node_id)
                                            .map(|pos| (row.node_id.clone(), pointer.distance_sq(*pos)))
                                    })
                                    .min_by(|a, b| {
                                        a.1.partial_cmp(&b.1)
                                            .unwrap_or(std::cmp::Ordering::Equal)
                                    })
                                    .map(|(node_id, _)| node_id)
                            });
                            let mut used_label_rects: Vec<egui::Rect> = Vec::new();
                            let mut op_label_galleys: HashMap<String, Arc<egui::Galley>> =
                                HashMap::new();
                            let mut retrieval_badge_rects: HashMap<String, egui::Rect> =
                                HashMap::new();
                            let edge_groups = Self::lineage_edge_groups(lineage_edges);
                            for (edge_idx, (from_node, to_node, op_ids)) in
                                edge_groups.iter().enumerate()
                            {
                                let Some(from) = pos_by_node.get(from_node) else {
                                    continue;
                                };
                                let Some(to) = pos_by_node.get(to_node) else {
                                    continue;
                                };
                                let mid = Pos2::new((from.x + to.x) * 0.5, (from.y + to.y) * 0.5);
                                let op_labels = Self::lineage_edge_operation_labels(op_ids, op_label_by_id);
                                let op_label = Self::lineage_edge_display_label(&op_labels);
                                let edge = *to - *from;
                                let edge_len = edge.length();
                                if edge_len < 0.1 {
                                    continue;
                                }
                                let edge_dir = edge / edge_len;
                                if Self::is_lineage_operation_hub_connector(op_ids) {
                                    painter.line_segment(
                                        [*from, *to],
                                        egui::Stroke::new(edge_stroke_width, egui::Color32::GRAY),
                                    );
                                    continue;
                                }
                                let op_node_size = (14.0 * graph_zoom).clamp(10.0, 20.0);
                                let op_rect =
                                    egui::Rect::from_center_size(mid, Vec2::splat(op_node_size));
                                let is_reopenable_gibson_op = op_ids.len() == 1
                                    && self
                                        .lineage_reopenable_gibson_op_ids
                                        .contains(&op_ids[0]);
                                let is_reopenable_pcr_op = op_ids.len() == 1
                                    && self
                                        .lineage_reopenable_pcr_op_seq_ids
                                        .contains_key(&op_ids[0]);
                                let mut op_response = ui.interact(
                                    op_rect,
                                    ui.id().with(("lineage_op", edge_idx)),
                                    egui::Sense::click(),
                                );
                                let op_accessible_name =
                                    Self::lineage_edge_accessible_name(&op_labels);
                                op_response = if is_reopenable_gibson_op {
                                    op_response.on_hover_text(format!(
                                        "{}\nClick to reopen the Gibson specialist for this cloning operation.",
                                        op_accessible_name
                                    ))
                                } else if is_reopenable_pcr_op {
                                    op_response.on_hover_text(format!(
                                        "{}\nClick to open the PCR Designer for this PCR-related operation.",
                                        op_accessible_name
                                    ))
                                } else {
                                    op_response.on_hover_text(op_accessible_name)
                                };
                                if is_reopenable_gibson_op && op_response.clicked() {
                                    reopen_gibson_from_operation = Some(op_ids[0].clone());
                                } else if is_reopenable_pcr_op && op_response.clicked() {
                                    reopen_pcr_from_operation = Some(op_ids[0].clone());
                                }
                                let op_gap = op_node_size * 0.5 + 2.0 * graph_zoom;
                                if edge_len > op_gap * 2.0 + 2.0 * graph_zoom {
                                    painter.line_segment(
                                        [*from, mid - edge_dir * op_gap],
                                        egui::Stroke::new(edge_stroke_width, egui::Color32::GRAY),
                                    );
                                    painter.line_segment(
                                        [mid + edge_dir * op_gap, *to],
                                        egui::Stroke::new(edge_stroke_width, egui::Color32::GRAY),
                                    );
                                } else {
                                    painter.line_segment(
                                        [*from, *to],
                                        egui::Stroke::new(edge_stroke_width, egui::Color32::GRAY),
                                    );
                                }
                                let (op_symbol, op_fill_color) = Self::lineage_operation_glyph(&op_labels);
                                painter.rect_filled(op_rect, 2.0 * graph_zoom, op_fill_color);
                                painter.rect_stroke(
                                    op_rect,
                                    2.0 * graph_zoom,
                                        egui::Stroke::new(
                                        if (is_reopenable_gibson_op || is_reopenable_pcr_op)
                                            && op_response.hovered()
                                        {
                                            (1.8 * graph_zoom).clamp(1.2, 2.6)
                                        } else {
                                            (1.1 * graph_zoom).clamp(1.0, 2.0)
                                        },
                                        if (is_reopenable_gibson_op || is_reopenable_pcr_op)
                                            && op_response.hovered()
                                        {
                                            egui::Color32::from_rgb(250, 220, 80)
                                        } else {
                                            egui::Color32::from_rgb(55, 55, 55)
                                        },
                                    ),
                                    egui::StrokeKind::Inside,
                                );
                                let op_symbol_font_size = if op_symbol.len() > 1 {
                                    (8.4 * graph_zoom).clamp(7.0, 10.5)
                                } else {
                                    (9.4 * graph_zoom).clamp(7.8, 11.8)
                                };
                                painter.text(
                                    mid,
                                    egui::Align2::CENTER_CENTER,
                                    op_symbol,
                                    egui::FontId::monospace(op_symbol_font_size),
                                    egui::Color32::WHITE,
                                );
                                used_label_rects.push(op_rect.expand(2.0 * graph_zoom));
                                if edge_label_stride > 1 && edge_idx % edge_label_stride != 0 {
                                    continue;
                                }
                                let display = if simplify_labels {
                                    Self::compact_lineage_op_label(&op_label)
                                } else {
                                    Self::lineage_edge_accessible_name(&op_labels)
                                };
                                let galley = op_label_galleys
                                    .entry(display.clone())
                                    .or_insert_with(|| {
                                        painter.layout_no_wrap(
                                            display.clone(),
                                            egui::FontId::proportional(op_font_size),
                                            egui::Color32::BLACK,
                                        )
                                    })
                                    .clone();
                                let perp = Vec2::new(-edge_dir.y, edge_dir.x);
                                let mut placed = None;
                                for idx in 0..14 {
                                    let sign = if idx % 2 == 0 { 1.0 } else { -1.0 };
                                    let step = (idx / 2) as f32;
                                    let candidate_center =
                                        mid + perp * (12.0 + step * 12.0) * graph_zoom * sign;
                                    let candidate_rect = egui::Rect::from_center_size(
                                        candidate_center,
                                        Vec2::new(
                                            galley.size().x + 10.0 * graph_zoom,
                                            galley.size().y + 6.0 * graph_zoom,
                                        ),
                                    );
                                    if !used_label_rects
                                        .iter()
                                        .any(|r| r.intersects(candidate_rect))
                                    {
                                        placed = Some((candidate_center, candidate_rect));
                                        break;
                                    }
                                }
                                let (label_center, bg_rect) = placed.unwrap_or_else(|| {
                                    let fallback_center = mid + perp * 12.0 * graph_zoom;
                                    let rect = egui::Rect::from_center_size(
                                        fallback_center,
                                        Vec2::new(
                                            galley.size().x + 10.0 * graph_zoom,
                                            galley.size().y + 6.0 * graph_zoom,
                                        ),
                                    );
                                    (fallback_center, rect)
                                });
                                used_label_rects.push(bg_rect);
                                painter.rect_filled(
                                    bg_rect,
                                    3.0 * graph_zoom,
                                    egui::Color32::from_rgba_premultiplied(245, 245, 245, 235),
                                );
                                painter.text(
                                    label_center,
                                    egui::Align2::CENTER_CENTER,
                                    display,
                                    egui::FontId::proportional(op_font_size),
                                    egui::Color32::BLACK,
                                );
                            }
                            for row in rows {
                                let Some(pos) = pos_by_node.get(&row.node_id).cloned() else {
                                    continue;
                                };
                                let is_selected = self
                                    .lineage_graph_selected_node_id
                                    .as_ref()
                                    .is_some_and(|node_id| node_id == &row.node_id);
                                let is_hovered = hovered_node_id
                                    .as_ref()
                                    .is_some_and(|node_id| node_id == &row.node_id);
                                let highlight_stroke = if is_selected {
                                    egui::Stroke::new(
                                        (2.4 * graph_zoom).clamp(1.6, 4.0),
                                        egui::Color32::from_rgb(250, 220, 80),
                                    )
                                } else if is_hovered {
                                    egui::Stroke::new(
                                        (2.0 * graph_zoom).clamp(1.4, 3.2),
                                        egui::Color32::from_rgb(230, 230, 150),
                                    )
                                } else {
                                    egui::Stroke::new(
                                        edge_stroke_width,
                                        egui::Color32::from_rgb(235, 196, 150),
                                    )
                                };
                                match row.kind {
                                    LineageNodeKind::Arrangement => {
                                        let rect = egui::Rect::from_center_size(
                                            pos,
                                            Vec2::new(56.0 * graph_zoom, 28.0 * graph_zoom),
                                        );
                                        painter.rect_filled(
                                            rect,
                                            5.0 * graph_zoom,
                                            if is_selected {
                                                egui::Color32::from_rgb(98, 140, 112)
                                            } else {
                                                egui::Color32::from_rgb(108, 154, 122)
                                            },
                                        );
                                        painter.rect_stroke(
                                            rect,
                                            5.0 * graph_zoom,
                                            highlight_stroke,
                                            egui::StrokeKind::Inside,
                                        );
                                    }
                                    LineageNodeKind::Macro => {
                                        let rect = egui::Rect::from_center_size(
                                            pos,
                                            Vec2::new(64.0 * graph_zoom, 30.0 * graph_zoom),
                                        );
                                        painter.rect_filled(
                                            rect,
                                            4.0 * graph_zoom,
                                            if is_selected {
                                                egui::Color32::from_rgb(118, 118, 128)
                                            } else {
                                                egui::Color32::from_rgb(98, 98, 108)
                                            },
                                        );
                                        painter.rect_stroke(
                                            rect,
                                            4.0 * graph_zoom,
                                            highlight_stroke,
                                            egui::StrokeKind::Inside,
                                        );
                                    }
                                    LineageNodeKind::Analysis => {
                                        let rect = egui::Rect::from_center_size(
                                            pos,
                                            if Self::is_lineage_operation_hub(row) {
                                                Vec2::new(74.0 * graph_zoom, 30.0 * graph_zoom)
                                            } else {
                                                Vec2::new(62.0 * graph_zoom, 28.0 * graph_zoom)
                                            },
                                        );
                                        painter.rect_filled(
                                            rect,
                                            4.0 * graph_zoom,
                                            if Self::is_lineage_operation_hub(row) {
                                                if is_selected {
                                                    egui::Color32::from_rgb(76, 162, 142)
                                                } else {
                                                    egui::Color32::from_rgb(62, 142, 124)
                                                }
                                            } else if is_selected {
                                                egui::Color32::from_rgb(156, 112, 184)
                                            } else {
                                                egui::Color32::from_rgb(140, 98, 172)
                                            },
                                        );
                                        painter.rect_stroke(
                                            rect,
                                            4.0 * graph_zoom,
                                            highlight_stroke,
                                            egui::StrokeKind::Inside,
                                        );
                                    }
                                    LineageNodeKind::Sequence if row.pool_size > 1 => {
                                        let points = vec![
                                            pos + Vec2::new(0.0, -16.0 * graph_zoom),
                                            pos + Vec2::new(16.0 * graph_zoom, 0.0),
                                            pos + Vec2::new(0.0, 16.0 * graph_zoom),
                                            pos + Vec2::new(-16.0 * graph_zoom, 0.0),
                                        ];
                                        painter.add(egui::Shape::convex_polygon(
                                            points,
                                            if is_selected {
                                                egui::Color32::from_rgb(205, 140, 80)
                                            } else {
                                                egui::Color32::from_rgb(180, 120, 70)
                                            },
                                            highlight_stroke,
                                        ));
                                    }
                                    LineageNodeKind::Sequence => {
                                        painter.circle_filled(
                                            pos,
                                            node_radius,
                                            if is_selected {
                                                egui::Color32::from_rgb(70, 125, 215)
                                            } else {
                                                egui::Color32::from_rgb(90, 140, 210)
                                            },
                                        );
                                        painter.circle_stroke(pos, node_radius, highlight_stroke);
                                    }
                                }
                                let node_id_label = match row.kind {
                                    LineageNodeKind::Arrangement => row
                                        .arrangement_id
                                        .as_ref()
                                        .map(|id| Self::compact_lineage_node_label(id, 12))
                                        .unwrap_or_else(|| {
                                            Self::compact_lineage_node_label(&row.node_id, 12)
                                        }),
                                    LineageNodeKind::Macro => row
                                        .macro_instance_id
                                        .as_ref()
                                        .map(|id| Self::compact_lineage_node_label(id, 12))
                                        .unwrap_or_else(|| {
                                            Self::compact_lineage_node_label(&row.node_id, 12)
                                        }),
                                    LineageNodeKind::Analysis => row
                                        .analysis_artifact_id
                                        .as_ref()
                                        .map(|id| Self::compact_lineage_node_label(id, 12))
                                        .unwrap_or_else(|| {
                                            if Self::is_lineage_operation_hub(row) {
                                                "GB".to_string()
                                            } else {
                                                Self::compact_lineage_node_label(&row.node_id, 12)
                                            }
                                        }),
                                    LineageNodeKind::Sequence => {
                                        Self::compact_lineage_node_label(&row.node_id, 10)
                                    }
                                };
                                painter.text(
                                    pos,
                                    egui::Align2::CENTER_CENTER,
                                    node_id_label,
                                    egui::FontId::monospace(node_id_font_size),
                                    egui::Color32::WHITE,
                                );
                                if let Some(descriptor) = &row.retrieval_descriptor {
                                    let badge_center = pos + Vec2::new(16.0 * graph_zoom, -18.0 * graph_zoom);
                                    let badge_rect = egui::Rect::from_center_size(
                                        badge_center,
                                        Vec2::new(
                                            (28.0 * graph_zoom).clamp(18.0, 36.0),
                                            (14.0 * graph_zoom).clamp(10.0, 20.0),
                                        ),
                                    );
                                    painter.rect_filled(
                                        badge_rect,
                                        3.0 * graph_zoom,
                                        egui::Color32::from_rgb(188, 146, 48),
                                    );
                                    painter.rect_stroke(
                                        badge_rect,
                                        3.0 * graph_zoom,
                                        egui::Stroke::new(
                                            (1.0 * graph_zoom).clamp(0.8, 1.8),
                                            egui::Color32::from_rgb(80, 62, 20),
                                        ),
                                        egui::StrokeKind::Inside,
                                    );
                                    painter.text(
                                        badge_rect.center(),
                                        egui::Align2::CENTER_CENTER,
                                        Self::lineage_retrieval_pattern_label(descriptor),
                                        egui::FontId::monospace((8.0 * graph_zoom).clamp(6.8, 10.5)),
                                        egui::Color32::BLACK,
                                    );
                                    retrieval_badge_rects.insert(row.node_id.clone(), badge_rect);
                                }
                                if let Some(badge) = collapsed_group_badges.get(&row.node_id) {
                                    let chip_font =
                                        egui::FontId::monospace((8.2 * graph_zoom).clamp(6.8, 10.8));
                                    let chip_spacing = 3.0 * graph_zoom;
                                    let mut chips: Vec<(String, egui::Color32)> =
                                        Vec::with_capacity(5);
                                    chips.push((
                                        badge.total_ops.min(99).to_string(),
                                        egui::Color32::from_rgb(52, 52, 52),
                                    ));
                                    for (family, count) in badge.families.iter().take(3) {
                                        let symbol = Self::lineage_operation_symbol_for_family(family);
                                        let label = if *count > 1 {
                                            format!("{}{}", symbol, (*count).min(99))
                                        } else {
                                            symbol.to_string()
                                        };
                                        chips.push((label, Self::lineage_operation_color_for_family(family)));
                                    }
                                    if badge.families.len() > 3 {
                                        chips.push((
                                            format!("+{}", badge.families.len() - 3),
                                            egui::Color32::from_rgb(95, 95, 95),
                                        ));
                                    }

                                    let mut chip_sizes: Vec<Vec2> = Vec::with_capacity(chips.len());
                                    let mut total_width = 0.0f32;
                                    for (label, _) in &chips {
                                        let galley = painter.layout_no_wrap(
                                            label.clone(),
                                            chip_font.clone(),
                                            egui::Color32::WHITE,
                                        );
                                        let size = Vec2::new(
                                            (galley.size().x + 8.0 * graph_zoom)
                                                .max(12.0 * graph_zoom),
                                            (galley.size().y + 4.0 * graph_zoom)
                                                .max(10.0 * graph_zoom),
                                        );
                                        total_width += size.x;
                                        chip_sizes.push(size);
                                    }
                                    total_width += chip_spacing * chips.len().saturating_sub(1) as f32;
                                    let mut chip_x = pos.x - total_width * 0.5;
                                    let chip_y = pos.y - 24.0 * graph_zoom;
                                    for ((label, fill), size) in chips.into_iter().zip(chip_sizes) {
                                        let chip_rect = egui::Rect::from_min_size(
                                            Pos2::new(chip_x, chip_y),
                                            size,
                                        );
                                        painter.rect_filled(chip_rect, 2.8 * graph_zoom, fill);
                                        painter.rect_stroke(
                                            chip_rect,
                                            2.8 * graph_zoom,
                                            egui::Stroke::new(
                                                (0.9 * graph_zoom).clamp(0.7, 1.6),
                                                egui::Color32::from_rgb(40, 40, 40),
                                            ),
                                            egui::StrokeKind::Inside,
                                        );
                                        painter.text(
                                            chip_rect.center(),
                                            egui::Align2::CENTER_CENTER,
                                            label,
                                            chip_font.clone(),
                                            egui::Color32::WHITE,
                                        );
                                        chip_x += size.x + chip_spacing;
                                    }
                                }
                                let display_name = if simplify_labels {
                                    Self::compact_lineage_node_label(&row.display_name, 26)
                                } else {
                                    row.display_name.clone()
                                };
                                painter.text(
                                    pos + Vec2::new(22.0 * graph_zoom, -4.0 * graph_zoom),
                                    egui::Align2::LEFT_BOTTOM,
                                    display_name,
                                    egui::FontId::proportional(name_font_size),
                                    egui::Color32::BLACK,
                                );
                                if !simplify_labels {
                                    let detail_text = match row.kind {
                                        LineageNodeKind::Arrangement => {
                                            let ladders = if row.ladders.is_empty() {
                                                "auto".to_string()
                                            } else {
                                                row.ladders.join(" + ")
                                            };
                                            format!(
                                                "{} lane(s) | mode={} | ladders={}",
                                                row.lane_container_ids.len(),
                                                row.arrangement_mode
                                                    .as_deref()
                                                    .unwrap_or("Serial")
                                                    .to_lowercase(),
                                                ladders
                                            )
                                        }
                                        LineageNodeKind::Macro => format!(
                                            "ops={} | template={} | status={}",
                                            row.macro_op_ids.len(),
                                            row.macro_template_name
                                                .as_deref()
                                                .unwrap_or("-"),
                                            row.macro_status.as_deref().unwrap_or("ok")
                                        ),
                                        LineageNodeKind::Analysis => {
                                            if Self::is_lineage_operation_hub(row) {
                                                format!("op={}", row.created_by_op)
                                            } else {
                                                let artifact_id = row
                                                    .analysis_artifact_id
                                                    .as_deref()
                                                    .unwrap_or(&row.display_name);
                                                let mode = row.analysis_mode.as_deref().unwrap_or("-");
                                                match row.analysis_kind {
                                                    Some(LineageAnalysisKind::Dotplot) => format!(
                                                        "dotplot={} | mode={} | points={}",
                                                        artifact_id,
                                                        mode,
                                                        row.analysis_point_count.unwrap_or(0)
                                                    ),
                                                    Some(LineageAnalysisKind::FlexibilityTrack) => {
                                                        format!(
                                                            "track={} | model={} | bins={}",
                                                            artifact_id,
                                                            mode,
                                                            row.analysis_bin_count.unwrap_or(0)
                                                        )
                                                    }
                                                    Some(LineageAnalysisKind::RnaReadInterpretation) => {
                                                        format!(
                                                            "rna_report={} | profile={} | status={} | reads={} | seed_passed={} | aligned={}",
                                                            artifact_id,
                                                            mode,
                                                            row.analysis_status
                                                                .as_deref()
                                                                .unwrap_or("-"),
                                                            row.analysis_read_count.unwrap_or(0),
                                                            row.analysis_point_count.unwrap_or(0),
                                                            row.analysis_variant_count.unwrap_or(0)
                                                        )
                                                    }
                                                    Some(LineageAnalysisKind::PrimerDesign) => {
                                                        format!(
                                                            "primer_report={} | backend={} | pairs={}",
                                                            artifact_id,
                                                            mode,
                                                            row.analysis_target_count.unwrap_or(0)
                                                        )
                                                    }
                                                    Some(LineageAnalysisKind::QpcrDesign) => {
                                                        format!(
                                                            "qpcr_report={} | backend={} | assays={}",
                                                            artifact_id,
                                                            mode,
                                                            row.analysis_target_count.unwrap_or(0)
                                                        )
                                                    }
                                                    Some(LineageAnalysisKind::RestrictionCloningPcrHandoff) => {
                                                        format!(
                                                            "restriction_handoff={} | mode={} | status={} | vector={}",
                                                            artifact_id,
                                                            mode,
                                                            row.analysis_status
                                                                .as_deref()
                                                                .unwrap_or("-"),
                                                            row.analysis_reference_seq_id
                                                                .as_deref()
                                                                .unwrap_or("-")
                                                        )
                                                    }
                                                    Some(LineageAnalysisKind::ProteinDerivation) => {
                                                        format!(
                                                            "protein_derivation={} | mode={} | proteins={}",
                                                            artifact_id,
                                                            mode,
                                                            row.analysis_target_count.unwrap_or(0)
                                                        )
                                                    }
                                                    Some(LineageAnalysisKind::ReverseTranslation) => {
                                                        format!(
                                                            "reverse_translation={} | profile={} | table={} | coding_bp={} | protein_aa={}",
                                                            artifact_id,
                                                            mode,
                                                            row.analysis_point_count.unwrap_or(0),
                                                            row.analysis_target_count.unwrap_or(0),
                                                            row.analysis_variant_count.unwrap_or(0)
                                                        )
                                                    }
                                                    Some(LineageAnalysisKind::ConstructReasoning) => {
                                                        format!(
                                                            "construct_reasoning={} | objective={} | evidence={} | decisions={} | candidates={}",
                                                            artifact_id,
                                                            mode,
                                                            row.analysis_point_count.unwrap_or(0),
                                                            row.analysis_variant_count.unwrap_or(0),
                                                            row.analysis_target_count.unwrap_or(0)
                                                        )
                                                    }
                                                    Some(LineageAnalysisKind::UniprotProjection) => {
                                                        format!(
                                                            "uniprot_projection={} | entry={} | transcripts={}",
                                                            artifact_id,
                                                            mode,
                                                            row.analysis_target_count.unwrap_or(0)
                                                        )
                                                    }
                                                    Some(LineageAnalysisKind::SequencingConfirmation) => {
                                                        format!(
                                                            "confirmation={} | status={} | reads={} | traces={} | targets={} | variants={}",
                                                            artifact_id,
                                                            row.analysis_status
                                                                .as_deref()
                                                                .unwrap_or("-"),
                                                            row.analysis_read_count.unwrap_or(0),
                                                            row.analysis_trace_count.unwrap_or(0),
                                                            row.analysis_target_count.unwrap_or(0),
                                                            row.analysis_variant_count.unwrap_or(0)
                                                        )
                                                    }
                                                    None => format!("analysis={artifact_id}"),
                                                }
                                            }
                                        }
                                        LineageNodeKind::Sequence if row.pool_size > 1 => {
                                            format!(
                                                "{} ({} bp) | pool={}",
                                                row.seq_id, row.length, row.pool_size
                                            )
                                        }
                                        LineageNodeKind::Sequence => {
                                            format!("{} ({} bp)", row.seq_id, row.length)
                                        }
                                    };
                                    painter.text(
                                        pos + Vec2::new(22.0 * graph_zoom, 10.0 * graph_zoom),
                                        egui::Align2::LEFT_TOP,
                                        detail_text,
                                        egui::FontId::proportional(details_font_size),
                                        egui::Color32::BLACK,
                                    );
                                }
                            }

                            let hit_row = hovered_node_id
                                .as_ref()
                                .and_then(|node_id| rows.iter().find(|row| &row.node_id == node_id));
                            let hovered_retrieval_node_id = pointer.and_then(|pointer| {
                                retrieval_badge_rects
                                    .iter()
                                    .find_map(|(node_id, rect)| rect.contains(pointer).then_some(node_id.clone()))
                            });
                            let hovered_retrieval_row = hovered_retrieval_node_id
                                .as_ref()
                                .and_then(|node_id| rows.iter().find(|row| &row.node_id == node_id));
                            let hover_row = hit_row.or(hovered_retrieval_row);
                            if let Some(row) = hover_row {
                                let mut hover_pool_range: Option<(usize, usize)> = None;
                                let mut hover_ladder_hint: Option<String> = None;
                                if row.kind == LineageNodeKind::Sequence && row.pool_size > 1 {
                                    let member_lengths: Vec<crate::pool_gel::GelSampleMember> = {
                                        let engine = self.engine.read().unwrap();
                                        row.pool_members
                                            .iter()
                                            .filter_map(|seq_id| {
                                                engine
                                                    .state()
                                                    .sequences
                                                    .get(seq_id)
                                                    .map(|dna| crate::pool_gel::GelSampleMember {
                                                        seq_id: seq_id.clone(),
                                                        bp: dna.len(),
                                                        topology_form: crate::engine::GentleEngine::infer_gel_topology_form_from_dna(dna),
                                                    })
                                            })
                                            .collect()
                                    };
                                    if !member_lengths.is_empty() {
                                        let min_bp = member_lengths
                                            .iter()
                                            .map(|member| member.bp)
                                            .min()
                                            .unwrap_or(0);
                                        let max_bp = member_lengths
                                            .iter()
                                            .map(|member| member.bp)
                                            .max()
                                            .unwrap_or(min_bp);
                                        hover_pool_range = Some((min_bp, max_bp));
                                        if let Ok(layout) =
                                            crate::pool_gel::build_pool_gel_layout(
                                                &member_lengths,
                                                &[],
                                                None,
                                            )
                                            && !layout.selected_ladders.is_empty() {
                                                hover_ladder_hint = Some(layout.selected_ladders.join(" + "));
                                            }
                                    }
                                }
                                let tooltip_member_preview = if row.kind == LineageNodeKind::Sequence
                                    && row.pool_size > 1
                                {
                                    row.pool_members
                                        .iter()
                                        .take(6)
                                        .cloned()
                                        .collect::<Vec<_>>()
                                        .join(", ")
                                } else {
                                    String::new()
                                };
                                resp.clone().on_hover_ui_at_pointer(|ui| {
                                    ui.strong(&row.display_name);
                                    match row.kind {
                                        LineageNodeKind::Arrangement => {
                                            let ladders = if row.ladders.is_empty() {
                                                "auto".to_string()
                                            } else {
                                                row.ladders.join(" + ")
                                            };
                                            ui.monospace(format!(
                                                "{} | mode={} | lanes={} | ladders={}",
                                                row.arrangement_id
                                                    .as_deref()
                                                    .unwrap_or(&row.seq_id),
                                                row.arrangement_mode
                                                    .as_deref()
                                                    .unwrap_or("Serial")
                                                    .to_lowercase(),
                                                row.lane_container_ids.len(),
                                                ladders
                                            ));
                                        }
                                        LineageNodeKind::Macro => {
                                            ui.monospace(format!(
                                                "{} | template={} | status={} | ops={}",
                                                row.macro_instance_id
                                                    .as_deref()
                                                    .unwrap_or(&row.seq_id),
                                                row.macro_template_name
                                                    .as_deref()
                                                    .unwrap_or("-"),
                                                row.macro_status.as_deref().unwrap_or("ok"),
                                                row.macro_op_ids.len()
                                            ));
                                            if let Some(routine_id) = &row.macro_routine_id {
                                                ui.small(format!("routine_id={routine_id}"));
                                            }
                                            if let Some(status_message) = &row.macro_status_message {
                                                ui.small(format!("status_message={status_message}"));
                                            }
                                        }
                                        LineageNodeKind::Analysis => {
                                            if Self::is_lineage_operation_hub(row) {
                                                ui.monospace(format!(
                                                    "{} | op={}",
                                                    row.display_name, row.created_by_op
                                                ));
                                                ui.small(
                                                    "Click to reopen the Gibson specialist for this cloning operation.",
                                                );
                                            } else {
                                                let artifact_id = row
                                                    .analysis_artifact_id
                                                    .as_deref()
                                                    .unwrap_or(&row.display_name);
                                                let mode = row.analysis_mode.as_deref().unwrap_or("-");
                                                let analysis_kind = row
                                                    .analysis_kind
                                                    .map(LineageAnalysisKind::as_str)
                                                    .unwrap_or("analysis");
                                                ui.monospace(format!(
                                                    "{} | kind={} | mode={}",
                                                    artifact_id, analysis_kind, mode
                                                ));
                                                match row.analysis_kind {
                                                    Some(LineageAnalysisKind::Dotplot) => {
                                                        ui.small(format!(
                                                            "query={} | reference={} | points={}",
                                                            row.seq_id,
                                                            row.analysis_reference_seq_id
                                                                .as_deref()
                                                                .unwrap_or(&row.seq_id),
                                                            row.analysis_point_count.unwrap_or(0)
                                                        ));
                                                    }
                                                    Some(LineageAnalysisKind::FlexibilityTrack) => {
                                                        ui.small(format!(
                                                            "seq={} | bins={}",
                                                            row.seq_id,
                                                            row.analysis_bin_count.unwrap_or(0)
                                                        ));
                                                    }
                                                    Some(LineageAnalysisKind::RnaReadInterpretation) => {
                                                        ui.small(format!(
                                                            "seq={} | profile={} | status={} | reads={} | seed_passed={} | aligned={} | target_genes={}",
                                                            row.seq_id,
                                                            row.analysis_mode
                                                                .as_deref()
                                                                .unwrap_or("-"),
                                                            row.analysis_status
                                                                .as_deref()
                                                                .unwrap_or("-"),
                                                            row.analysis_read_count.unwrap_or(0),
                                                            row.analysis_point_count.unwrap_or(0),
                                                            row.analysis_variant_count.unwrap_or(0),
                                                            row.analysis_target_count.unwrap_or(0)
                                                        ));
                                                    }
                                                    Some(LineageAnalysisKind::PrimerDesign) => {
                                                        ui.small(format!(
                                                            "template={} | backend={} | pairs={}",
                                                            row.seq_id,
                                                            row.analysis_mode
                                                                .as_deref()
                                                                .unwrap_or("-"),
                                                            row.analysis_target_count.unwrap_or(0)
                                                        ));
                                                    }
                                                    Some(LineageAnalysisKind::QpcrDesign) => {
                                                        ui.small(format!(
                                                            "template={} | backend={} | assays={}",
                                                            row.seq_id,
                                                            row.analysis_mode
                                                                .as_deref()
                                                                .unwrap_or("-"),
                                                            row.analysis_target_count.unwrap_or(0)
                                                        ));
                                                    }
                                                    Some(LineageAnalysisKind::RestrictionCloningPcrHandoff) => {
                                                        ui.small(format!(
                                                            "template={} | vector={} | mode={} | status={}",
                                                            row.seq_id,
                                                            row.analysis_reference_seq_id
                                                                .as_deref()
                                                                .unwrap_or("-"),
                                                            row.analysis_mode
                                                                .as_deref()
                                                                .unwrap_or("-"),
                                                            row.analysis_status
                                                                .as_deref()
                                                                .unwrap_or("-")
                                                        ));
                                                    }
                                                    Some(LineageAnalysisKind::ProteinDerivation) => {
                                                        ui.small(format!(
                                                            "source={} | mode={} | proteins={}",
                                                            row.seq_id,
                                                            row.analysis_mode
                                                                .as_deref()
                                                                .unwrap_or("-"),
                                                            row.analysis_target_count.unwrap_or(0)
                                                        ));
                                                    }
                                                    Some(LineageAnalysisKind::ReverseTranslation) => {
                                                        ui.small(format!(
                                                            "protein={} | product={} | profile={} | table={} | coding_bp={} | protein_aa={}",
                                                            row.seq_id,
                                                            row.analysis_reference_seq_id
                                                                .as_deref()
                                                                .unwrap_or("-"),
                                                            row.analysis_mode
                                                                .as_deref()
                                                                .unwrap_or("-"),
                                                            row.analysis_point_count.unwrap_or(0),
                                                            row.analysis_target_count.unwrap_or(0),
                                                            row.analysis_variant_count.unwrap_or(0)
                                                        ));
                                                    }
                                                    Some(LineageAnalysisKind::ConstructReasoning) => {
                                                        ui.small(format!(
                                                            "seq={} | objective={} | goal={} | evidence={} | decisions={} | candidates={}",
                                                            row.seq_id,
                                                            row.analysis_mode
                                                                .as_deref()
                                                                .unwrap_or("-"),
                                                            row.analysis_status
                                                                .as_deref()
                                                                .unwrap_or("-"),
                                                            row.analysis_point_count.unwrap_or(0),
                                                            row.analysis_variant_count.unwrap_or(0),
                                                            row.analysis_target_count.unwrap_or(0)
                                                        ));
                                                    }
                                                    Some(LineageAnalysisKind::UniprotProjection) => {
                                                        ui.small(format!(
                                                            "seq={} | entry={} | transcript_filter={} | transcripts={}",
                                                            row.seq_id,
                                                            row.analysis_mode
                                                                .as_deref()
                                                                .unwrap_or("-"),
                                                            row.analysis_status
                                                                .as_deref()
                                                                .unwrap_or("all_transcripts"),
                                                            row.analysis_target_count.unwrap_or(0)
                                                        ));
                                                    }
                                                    Some(LineageAnalysisKind::SequencingConfirmation) => {
                                                        ui.small(format!(
                                                            "expected={} | baseline={} | status={} | reads={} | traces={} | targets={} | variants={}",
                                                            row.seq_id,
                                                            row.analysis_reference_seq_id
                                                                .as_deref()
                                                                .unwrap_or("-"),
                                                            row.analysis_status
                                                                .as_deref()
                                                                .unwrap_or("-"),
                                                            row.analysis_read_count.unwrap_or(0),
                                                            row.analysis_trace_count.unwrap_or(0),
                                                            row.analysis_target_count.unwrap_or(0),
                                                            row.analysis_variant_count.unwrap_or(0)
                                                        ));
                                                    }
                                                    None => {}
                                                }
                                            }
                                        }
                                        LineageNodeKind::Sequence => {
                                            ui.monospace(format!(
                                                "{} | {} bp | {}",
                                                row.seq_id,
                                                row.length,
                                                if row.circular { "circular" } else { "linear" }
                                            ));
                                        }
                                    }
                                    ui.small(format!("node={} | origin={}", row.node_id, row.origin));
                                    ui.small(format!("parents={} | op={}", row.parents.len(), row.created_by_op));
                                    if let Some(badge) = collapsed_group_badges.get(&row.node_id) {
                                        ui.separator();
                                        ui.small(format!(
                                            "collapsed hidden operations={}",
                                            badge.total_ops
                                        ));
                                        ui.small(format!(
                                            "hidden families={}",
                                            Self::lineage_hidden_op_families_summary(badge, 5)
                                        ));
                                    }
                                    if let Some(anchor) = &row.genome_anchor_summary {
                                        let strand = anchor.strand.unwrap_or('+');
                                        let mut anchor_text = format!(
                                            "anchor={} {}:{}-{} (strand {})",
                                            anchor.genome_id,
                                            anchor.chromosome,
                                            anchor.start_1based,
                                            anchor.end_1based,
                                            strand
                                        );
                                        if row.is_full_genome_sequence {
                                            anchor_text
                                                .push_str(" | full chromosome/genome sequence");
                                        }
                                        ui.small(anchor_text);
                                    }
                                    if let Some(descriptor) = row.retrieval_descriptor.as_ref() {
                                        ui.separator();
                                        ui.small(Self::lineage_retrieval_pattern_tooltip(
                                            descriptor,
                                            row.genome_anchor_summary.as_ref(),
                                        ));
                                    }
                                    match row.kind {
                                        LineageNodeKind::Arrangement => {
                                            if !row.lane_container_ids.is_empty() {
                                                ui.separator();
                                                let preview = row
                                                    .lane_container_ids
                                                    .iter()
                                                    .take(6)
                                                    .cloned()
                                                    .collect::<Vec<_>>()
                                                    .join(", ");
                                                ui.small(format!("lane containers: {preview}"));
                                                if row.lane_container_ids.len() > 6 {
                                                    ui.small(format!(
                                                        "... and {} more",
                                                        row.lane_container_ids.len() - 6
                                                    ));
                                                }
                                            }
                                        }
                                        LineageNodeKind::Macro => {
                                            if !row.macro_inputs.is_empty() || !row.macro_outputs.is_empty() {
                                                ui.separator();
                                                if !row.macro_inputs.is_empty() {
                                                    let preview = row
                                                        .macro_inputs
                                                        .iter()
                                                        .map(|binding| {
                                                            format!(
                                                                "{}={}",
                                                                binding.port_id,
                                                                binding.values.join(",")
                                                            )
                                                        })
                                                        .take(4)
                                                        .collect::<Vec<_>>()
                                                        .join(" | ");
                                                    if !preview.is_empty() {
                                                        ui.small(format!("inputs: {preview}"));
                                                    }
                                                }
                                                if !row.macro_outputs.is_empty() {
                                                    let preview = row
                                                        .macro_outputs
                                                        .iter()
                                                        .map(|binding| {
                                                            format!(
                                                                "{}={}",
                                                                binding.port_id,
                                                                binding.values.join(",")
                                                            )
                                                        })
                                                        .take(4)
                                                        .collect::<Vec<_>>()
                                                        .join(" | ");
                                                    if !preview.is_empty() {
                                                        ui.small(format!("outputs: {preview}"));
                                                    }
                                                }
                                            }
                                            if !row.macro_op_ids.is_empty() {
                                                let preview = row
                                                    .macro_op_ids
                                                    .iter()
                                                    .take(6)
                                                    .cloned()
                                                    .collect::<Vec<_>>()
                                                    .join(", ");
                                                ui.small(format!("op_ids: {preview}"));
                                            }
                                        }
                                        LineageNodeKind::Analysis => {
                                            if let Some(artifact_id) =
                                                row.analysis_artifact_id.as_deref()
                                            {
                                                ui.separator();
                                                ui.small(format!("artifact_id={artifact_id}"));
                                            }
                                            if let Some(reference_seq_id) =
                                                row.analysis_reference_seq_id.as_deref()
                                            {
                                                ui.small(format!("reference_seq_id={reference_seq_id}"));
                                            }
                                        }
                                        LineageNodeKind::Sequence if row.pool_size > 1 => {
                                            ui.separator();
                                            ui.small(format!("pool members={}", row.pool_size));
                                            if let Some((min_bp, max_bp)) = hover_pool_range {
                                                ui.small(format!(
                                                    "pool range={}..{} bp",
                                                    min_bp, max_bp
                                                ));
                                            }
                                            if let Some(ladders) = &hover_ladder_hint {
                                                ui.small(format!("suggested ladders={}", ladders));
                                            }
                                            if !tooltip_member_preview.is_empty() {
                                                ui.small(format!(
                                                    "members: {}",
                                                    tooltip_member_preview
                                                ));
                                            }
                                        }
                                        LineageNodeKind::Sequence => {}
                                    }
                                });
                            }
                            let context_node_id = hover_row.map(|row| row.node_id.clone());
                            let context_row = hover_row.cloned();
                            resp.clone().context_menu(|ui| {
                                if let Some(node_id) = context_node_id.as_ref() {
                                    self.render_lineage_group_context_menu(
                                        ui,
                                        node_id,
                                        context_row.as_ref(),
                                        &leaf_node_ids,
                                        &valid_lineage_node_ids,
                                        &mut persist_workspace_after_frame,
                                    );
                                } else {
                                    ui.label("No node under cursor");
                                }
                                ui.separator();
                                if ui
                                    .button("Save Graph as SVG...")
                                    .on_hover_text(
                                        "Export the currently visible lineage graph view as SVG",
                                    )
                                    .clicked()
                                {
                                    request_export_lineage_svg = true;
                                    ui.close();
                                }
                            });
                            let (space_pan_requested, option_pan_requested, modifiers) = ui
                                .input(|i| {
                                    (
                                        i.key_down(Key::Space),
                                        scroll_input_policy::option_pan_modifier_active(i.modifiers),
                                        i.modifiers,
                                    )
                                });
                            if resp.drag_started() {
                                if option_pan_requested || (space_pan_requested && hover_row.is_none())
                                {
                                    self.lineage_graph_pan_origin = Some(graph_scroll_offset);
                                    self.lineage_graph_drag_origin = None;
                                } else if let Some(row) = hover_row {
                                    let start_offset = self
                                        .lineage_graph_node_offsets
                                        .get(&row.node_id)
                                        .copied()
                                        .unwrap_or(Vec2::ZERO);
                                    self.lineage_graph_drag_origin =
                                        Some((row.node_id.clone(), start_offset));
                                    self.lineage_graph_pan_origin = None;
                                    self.lineage_graph_selected_node_id = Some(row.node_id.clone());
                                } else {
                                    self.lineage_graph_drag_origin = None;
                                    self.lineage_graph_pan_origin = None;
                                }
                            }
                            if let Some(start_scroll_offset) = self.lineage_graph_pan_origin {
                                if resp.dragged() {
                                    let next_offset = start_scroll_offset - resp.drag_delta();
                                    graph_scroll_offset = Vec2::new(
                                        next_offset.x.max(0.0),
                                        next_offset.y.max(0.0),
                                    );
                                }
                                if resp.drag_stopped() {
                                    self.lineage_graph_pan_origin = None;
                                    persist_workspace_after_frame = true;
                                }
                            }
                            if self.lineage_graph_pan_origin.is_none()
                                && let Some((node_id, start_offset)) =
                                    self.lineage_graph_drag_origin.clone()
                                {
                                    if resp.dragged() {
                                        self.lineage_graph_node_offsets
                                            .insert(node_id.clone(), start_offset + resp.drag_delta());
                                    }
                                    if resp.drag_stopped() {
                                        self.lineage_graph_drag_origin = None;
                                        persist_workspace_after_frame = true;
                                    }
                                }
                            if let Some(cursor) = scroll_input_policy::canvas_hover_cursor(
                                modifiers,
                                true,
                                self.lineage_graph_pan_origin.is_some(),
                                self.lineage_graph_drag_origin.is_some(),
                                graph_wheel_intent,
                            ) {
                                ui.output_mut(|o| o.cursor_icon = cursor);
                            } else if space_pan_requested && hover_row.is_none() {
                                ui.output_mut(|o| o.cursor_icon = egui::CursorIcon::Grab);
                            }
                            if self.lineage_graph_pan_origin.is_none()
                                && self.lineage_graph_drag_origin.is_none()
                            {
                                let mut retrieval_click_consumed = false;
                                if resp.clicked()
                                    && let Some(row) = hovered_retrieval_row
                                        && let Some(descriptor) = row.retrieval_descriptor.clone() {
                                            self.lineage_graph_selected_node_id =
                                                Some(row.node_id.clone());
                                            open_lineage_retrieval = Some(descriptor);
                                            retrieval_click_consumed = true;
                                        }
                                if !retrieval_click_consumed {
                                    if resp.clicked() {
                                        self.lineage_graph_selected_node_id =
                                            hover_row.map(|row| row.node_id.clone());
                                        if let Some(row) = hover_row
                                            && Self::is_lineage_operation_hub(row)
                                                && self
                                                    .lineage_reopenable_gibson_op_ids
                                                    .contains(&row.created_by_op)
                                            {
                                                reopen_gibson_from_operation =
                                                    Some(row.created_by_op.clone());
                                            }
                                    }
                                    if let Some(row) = hover_row
                                        && resp.double_clicked() {
                                            match row.kind {
                                                LineageNodeKind::Arrangement => {}
                                                LineageNodeKind::Macro => {}
                                                LineageNodeKind::Analysis => {
                                                    if Self::is_lineage_operation_hub(row) {
                                                        reopen_gibson_from_operation =
                                                            Some(row.created_by_op.clone());
                                                    } else {
                                                        if let Some((kind, seq_id, artifact_id)) =
                                                            Self::lineage_analysis_open_payload(row)
                                                        {
                                                            open_lineage_analysis =
                                                                Some((kind, seq_id, artifact_id));
                                                        } else {
                                                            open_seq = Some(row.seq_id.clone());
                                                        }
                                                    }
                                                }
                                                LineageNodeKind::Sequence if row.pool_size > 1 => {
                                                    open_pool = Some((
                                                        row.seq_id.clone(),
                                                        row.pool_members.clone(),
                                                    ));
                                                }
                                                LineageNodeKind::Sequence => {
                                                    open_seq = Some(row.seq_id.clone());
                                                }
                                            }
                                        }
                                }
                            }
                            if self.lineage_graph_drag_origin.is_some()
                                && !ui.input(|i| i.pointer.primary_down())
                            {
                                self.lineage_graph_drag_origin = None;
                                persist_workspace_after_frame = true;
                            }
                            if self.lineage_graph_pan_origin.is_some()
                                && !ui.input(|i| i.pointer.primary_down())
                            {
                                self.lineage_graph_pan_origin = None;
                                persist_workspace_after_frame = true;
                            }
                        });
                    let max_scroll_x =
                        (graph_scroll_output.content_size.x - graph_scroll_output.inner_rect.width())
                            .max(0.0);
                    let max_scroll_y = (graph_scroll_output.content_size.y
                        - graph_scroll_output.inner_rect.height())
                    .max(0.0);
                    if self.lineage_graph_pan_origin.is_some() {
                        graph_scroll_offset.x = graph_scroll_offset.x.clamp(0.0, max_scroll_x);
                        graph_scroll_offset.y = graph_scroll_offset.y.clamp(0.0, max_scroll_y);
                    } else {
                        let mut measured_offset = graph_scroll_output.state.offset;
                        measured_offset.x = measured_offset.x.clamp(0.0, max_scroll_x);
                        measured_offset.y = measured_offset.y.clamp(0.0, max_scroll_y);
                        graph_scroll_offset = measured_offset;
                    }
                    });
                },
            );
        } else {
            let table_panel_width = ui.available_width().max(1.0);
            let table_panel_height = graph_area_height.max(LINEAGE_MAIN_TOP_PANEL_MIN_HEIGHT);
            let mut toggle_group_collapse_id: Option<String> = None;
            ui.allocate_ui_with_layout(
                egui::vec2(table_panel_width, table_panel_height),
                egui::Layout::top_down(egui::Align::Min),
                |ui| {
                    ui.set_min_width(table_panel_width);
                    ui.set_max_width(table_panel_width);
                    egui::ScrollArea::both()
                        .id_salt("lineage_table_scroll")
                        .auto_shrink([false, false])
                        .max_height(ui.available_height())
                        .show(ui, |ui| {
                    scroll_input_policy::apply_scrollarea_keyboard_navigation(
                        ui,
                        scroll_input_policy::DEFAULT_SCROLLAREA_KEYBOARD_STEP,
                    );
                    egui::Grid::new("lineage_table_grid")
                        .striped(true)
                        .min_col_width(72.0)
                        .spacing(egui::vec2(10.0, 6.0))
                        .show(ui, |ui| {
                            ui.strong("Node");
                            ui.strong("Sequence");
                            ui.strong("Parents");
                            ui.strong("Origin");
                            ui.strong("Op");
                            ui.strong("Length");
                            ui.strong("Topology");
                            ui.strong("Genome anchor");
                            ui.strong("Pattern");
                            ui.strong("Action");
                            ui.end_row();
                            for entry in &table_entries {
                                let row = &entry.row;
                                let node_display = Self::compact_lineage_node_label(&row.node_id, 10);
                                let mut node_response = ui
                                    .allocate_ui_with_layout(
                                        egui::vec2(160.0, ui.spacing().interact_size.y),
                                        egui::Layout::left_to_right(egui::Align::Center),
                                        |ui| {
                                            ui.horizontal(|ui| {
                                                ui.add_space(14.0 * entry.indent_level as f32);
                                                if entry.is_group_representative {
                                                    if let Some(group_id) = entry.group_id.as_ref()
                                                    {
                                                        let collapsed =
                                                            entry.hidden_group_member_count > 0;
                                                        if ui
                                                            .small_button(if collapsed {
                                                                "▸"
                                                            } else {
                                                                "▾"
                                                            })
                                                            .on_hover_text(
                                                                "Collapse/expand this node group",
                                                            )
                                                            .clicked()
                                                        {
                                                            toggle_group_collapse_id =
                                                                Some(group_id.clone());
                                                        }
                                                    }
                                                } else if entry.indent_level > 0 {
                                                    ui.label("↳");
                                                }
                                                ui.add(
                                                    egui::Label::new(
                                                        egui::RichText::new(node_display.clone())
                                                            .monospace(),
                                                    )
                                                    .truncate(),
                                                );
                                            });
                                        },
                                    )
                                    .response;
                                if node_display != row.node_id {
                                    node_response = node_response.on_hover_text(row.node_id.clone());
                                }
                                node_response.context_menu(|ui| {
                                    self.render_lineage_group_context_menu(
                                        ui,
                                        &row.node_id,
                                        Some(row),
                                        &leaf_node_ids,
                                        &valid_lineage_node_ids,
                                        &mut persist_workspace_after_frame,
                                    );
                                });
                                if node_response.clicked() {
                                    self.lineage_graph_selected_node_id = Some(row.node_id.clone());
                                }
                                match row.kind {
                                    LineageNodeKind::Arrangement => {
                                        ui.monospace(
                                            row.arrangement_id
                                                .as_deref()
                                                .unwrap_or(&row.seq_id)
                                                .to_string(),
                                        );
                                    }
                                    LineageNodeKind::Macro => {
                                        ui.monospace(
                                            row.macro_instance_id
                                                .as_deref()
                                                .unwrap_or(&row.seq_id)
                                                .to_string(),
                                        );
                                    }
                                    LineageNodeKind::Analysis => {
                                        ui.vertical(|ui| {
                                            ui.monospace(
                                                row.analysis_artifact_id
                                                    .as_deref()
                                                    .unwrap_or(&row.display_name),
                                            );
                                            if ui
                                                .small_button(row.seq_id.clone())
                                                .on_hover_text(
                                                    "Open the query/source sequence in a dedicated window",
                                                )
                                                .clicked()
                                            {
                                                open_seq = Some(row.seq_id.clone());
                                            }
                                        });
                                    }
                                    LineageNodeKind::Sequence => {
                                        if ui
                                            .button(&row.seq_id)
                                            .on_hover_text("Open this sequence in a dedicated window")
                                            .clicked()
                                        {
                                            open_seq = Some(row.seq_id.clone());
                                        }
                                    }
                                }
                                ui.label(if row.parents.is_empty() {
                                    "-".to_string()
                                } else {
                                    row.parents.join(" + ")
                                });
                                if entry.is_group_representative {
                                    if let Some(group_label) = entry.group_label.as_ref() {
                                        if entry.hidden_group_member_count > 0 {
                                            ui.label(format!(
                                                "{} [{} | {} hidden]",
                                                row.origin, group_label, entry.hidden_group_member_count
                                            ));
                                        } else {
                                            ui.label(format!("{} [{}]", row.origin, group_label));
                                        }
                                    } else {
                                        ui.label(&row.origin);
                                    }
                                } else {
                                    ui.label(&row.origin);
                                }
                                let op_display = Self::compact_lineage_node_label(&row.created_by_op, 11);
                                let is_reopenable_gibson_op =
                                    self.lineage_reopenable_gibson_op_ids.contains(&row.created_by_op);
                                let is_reopenable_pcr_op = self
                                    .lineage_reopenable_pcr_op_seq_ids
                                    .contains_key(&row.created_by_op);
                                if is_reopenable_gibson_op || is_reopenable_pcr_op {
                                    let response = ui.add_sized(
                                        egui::vec2(92.0, ui.spacing().interact_size.y),
                                        egui::Button::new(
                                            egui::RichText::new(op_display.clone()).monospace(),
                                        )
                                        .small(),
                                    );
                                    let response = if is_reopenable_gibson_op {
                                        response.on_hover_text(format!(
                                            "{}\nClick to reopen the Gibson specialist for this cloning operation.",
                                            row.created_by_op
                                        ))
                                    } else {
                                        response.on_hover_text(format!(
                                            "{}\nClick to open the PCR Designer for this PCR-related operation.",
                                            row.created_by_op
                                        ))
                                    };
                                    if response.clicked() {
                                        if is_reopenable_gibson_op {
                                            reopen_gibson_from_operation =
                                                Some(row.created_by_op.clone());
                                        } else {
                                            reopen_pcr_from_operation =
                                                Some(row.created_by_op.clone());
                                        }
                                    }
                                } else {
                                    let op_response = ui
                                        .allocate_ui_with_layout(
                                            egui::vec2(92.0, ui.spacing().interact_size.y),
                                            egui::Layout::left_to_right(egui::Align::Center),
                                            |ui| {
                                                ui.add(
                                                    egui::Label::new(
                                                        egui::RichText::new(op_display.clone())
                                                            .monospace(),
                                                    )
                                                    .truncate(),
                                                )
                                            },
                                        )
                                        .inner;
                                    if op_display != row.created_by_op {
                                        op_response.on_hover_text(row.created_by_op.clone());
                                    }
                                }
                                if row.kind == LineageNodeKind::Arrangement
                                    || row.kind == LineageNodeKind::Macro
                                    || row.kind == LineageNodeKind::Analysis
                                {
                                    ui.label("-");
                                    ui.label("-");
                                } else {
                                    ui.monospace(format!("{} bp", row.length));
                                    ui.label(if row.circular { "circular" } else { "linear" });
                                }
                                match row.kind {
                                    LineageNodeKind::Arrangement => {
                                        ui.label("-");
                                    }
                                    LineageNodeKind::Macro => {
                                        if ui
                                            .button("Inspect")
                                            .on_hover_text(
                                                "Show persistent details for this macro node below",
                                            )
                                            .clicked()
                                        {
                                            self.lineage_graph_selected_node_id =
                                                Some(row.node_id.clone());
                                        }
                                    }
                                    LineageNodeKind::Analysis => {
                                        ui.monospace(
                                            row.analysis_reference_seq_id
                                                .as_ref()
                                                .map(|reference| {
                                                    format!("query={} ref={reference}", row.seq_id)
                                                })
                                                .unwrap_or_else(|| {
                                                    format!("seq={}", row.seq_id)
                                                }),
                                        );
                                    }
                                    LineageNodeKind::Sequence => {
                                        let mut anchor_text = row
                                            .genome_anchor_display
                                            .clone()
                                            .unwrap_or_else(|| "-".to_string());
                                        if row.is_full_genome_sequence {
                                            anchor_text.push_str(" [full chromosome/genome]");
                                        }
                                        let response = ui.monospace(anchor_text);
                                        if let Some(anchor) = &row.genome_anchor_summary {
                                            let strand = anchor.strand.unwrap_or('+');
                                            let mut hover = format!(
                                                "Genome anchor: {}:{}-{} (strand {}, genome '{}')",
                                                anchor.chromosome,
                                                anchor.start_1based,
                                                anchor.end_1based,
                                                strand,
                                                anchor.genome_id
                                            );
                                            if row.is_full_genome_sequence {
                                                hover.push_str(
                                                    "\nCovers full chromosome/genome length in prepared reference.",
                                                );
                                            }
                                            response.on_hover_text(hover);
                                        }
                                    }
                                }
                                if let Some(descriptor) = row.retrieval_descriptor.as_ref() {
                                    let button = ui
                                        .small_button(Self::lineage_retrieval_pattern_label(
                                            descriptor,
                                        ))
                                        .on_hover_text(Self::lineage_retrieval_pattern_tooltip(
                                            descriptor,
                                            row.genome_anchor_summary.as_ref(),
                                        ));
                                    if button.clicked() {
                                        open_lineage_retrieval = Some(descriptor.clone());
                                        self.lineage_graph_selected_node_id =
                                            Some(row.node_id.clone());
                                    }
                                } else {
                                    ui.label("-");
                                }
                                match row.kind {
                                    LineageNodeKind::Arrangement => {
                                        if !row.lane_container_ids.is_empty() {
                                            self.render_open_lanes_menu(
                                                ui,
                                                &row.lane_container_ids,
                                                &mut open_lane_containers,
                                            );
                                        } else {
                                            ui.label("-");
                                        }
                                    }
                                    LineageNodeKind::Macro => {
                                        ui.label("-");
                                    }
                                    LineageNodeKind::Analysis => {
                                        let analysis_payload =
                                            Self::lineage_analysis_open_payload(row);
                                        let (label, hover) = match analysis_payload
                                            .as_ref()
                                            .map(|(kind, _, _)| *kind)
                                        {
                                            Some(LineageAnalysisKind::Dotplot) => (
                                                "Open Dotplot",
                                                "Open the dotplot analysis view for this sequence",
                                            ),
                                            Some(LineageAnalysisKind::FlexibilityTrack) => (
                                                "Open Track",
                                                "Open the flexibility-track analysis view for this sequence",
                                            ),
                                            Some(LineageAnalysisKind::RnaReadInterpretation) => (
                                                "Open RNA-read Mapping",
                                                "Open the RNA-read Mapping workspace on this persisted interpretation report",
                                            ),
                                            Some(LineageAnalysisKind::PrimerDesign) => (
                                                "Open Primer Report",
                                                "Open the PCR Designer on this persisted primer-design report",
                                            ),
                                            Some(LineageAnalysisKind::QpcrDesign) => (
                                                "Open qPCR Report",
                                                "Open the PCR Designer on this persisted qPCR-design report",
                                            ),
                                            Some(LineageAnalysisKind::RestrictionCloningPcrHandoff) => (
                                                "Open Cloning Handoff",
                                                "Open the PCR Designer on this persisted restriction-site cloning handoff",
                                            ),
                                            Some(LineageAnalysisKind::ProteinDerivation) => (
                                                "Open Derived Protein Expert",
                                                "Open the transcript-native Protein Expert on this persisted protein-derivation artifact",
                                            ),
                                            Some(LineageAnalysisKind::ReverseTranslation) => (
                                                "Open Coding Sequence",
                                                "Open the created coding-DNA product for this persisted reverse-translation artifact",
                                            ),
                                            Some(LineageAnalysisKind::ConstructReasoning) => (
                                                "Open Construct Reasoning",
                                                "Open the sequence window on this persisted construct-reasoning graph",
                                            ),
                                            Some(LineageAnalysisKind::UniprotProjection) => (
                                                "Open UniProt Projection",
                                                "Open the UniProt protein expert on this persisted genome projection",
                                            ),
                                            Some(LineageAnalysisKind::SequencingConfirmation) => (
                                                "Open Confirmation",
                                                "Open the sequencing-confirmation specialist on this persisted report",
                                            ),
                                            None => (
                                                "Open Seq",
                                                "Open the query/source sequence in a dedicated window",
                                            ),
                                        };
                                        ui.horizontal(|ui| {
                                            if ui
                                                .button(label)
                                                .on_hover_text(hover)
                                                .clicked()
                                            {
                                                if let Some((kind, seq_id, artifact_id)) =
                                                    analysis_payload.clone()
                                                {
                                                    open_lineage_analysis =
                                                        Some((kind, seq_id, artifact_id));
                                                } else {
                                                    open_seq = Some(row.seq_id.clone());
                                                }
                                            }
                                            if matches!(
                                                analysis_payload.as_ref().map(|(kind, _, _)| *kind),
                                                Some(LineageAnalysisKind::Dotplot)
                                            ) {
                                                let dotplot_id = row
                                                    .analysis_artifact_id
                                                    .clone()
                                                    .or_else(|| {
                                                        Self::infer_lineage_analysis_artifact_id_from_row(
                                                            row,
                                                        )
                                                    })
                                                    .unwrap_or_default();
                                                if !dotplot_id.trim().is_empty()
                                                    && ui
                                                        .button("Dotplot SVG")
                                                        .on_hover_text(
                                                            "Export this dotplot analysis as an SVG artifact",
                                                        )
                                                        .clicked()
                                                {
                                                    export_lineage_dotplot_svg = Some((
                                                        row.seq_id.clone(),
                                                        dotplot_id,
                                                        row.analysis_reference_seq_id.clone(),
                                                    ));
                                                }
                                            }
                                        });
                                    }
                                    LineageNodeKind::Sequence => {
                                        ui.allocate_ui_with_layout(
                                            egui::vec2(280.0, ui.spacing().interact_size.y),
                                            egui::Layout::left_to_right(egui::Align::Center),
                                            |ui| {
                                                if ui
                                                    .button("Select")
                                                    .on_hover_text(
                                                        "Run candidate selection operation using this sequence as input",
                                                    )
                                                    .clicked()
                                                {
                                                    select_candidate_from = Some(row.seq_id.clone());
                                                }
                                                if ui
                                                    .button("Copy")
                                                    .on_hover_text(
                                                        "Copy this sequence as plain DNA/RNA text",
                                                    )
                                                    .clicked()
                                                {
                                                    self.copy_project_sequence_to_clipboard(
                                                        ui.ctx(),
                                                        &row.seq_id,
                                                        false,
                                                    );
                                                }
                                                if ui
                                                    .button("FASTA")
                                                    .on_hover_text(
                                                        "Copy this sequence as FASTA text",
                                                    )
                                                    .clicked()
                                                {
                                                    self.copy_project_sequence_to_clipboard(
                                                        ui.ctx(),
                                                        &row.seq_id,
                                                        true,
                                                    );
                                                }
                                                if ui
                                                    .button("Remove")
                                                    .on_hover_text(
                                                        "Sequence-removal from lineage table is not yet implemented in engine operations.",
                                                    )
                                                    .clicked()
                                                {
                                                    self.app_status = format!(
                                                        "Remove for '{}' is not yet available",
                                                        row.seq_id
                                                    );
                                                }
                                            },
                                        );
                                    }
                                }
                                ui.end_row();
                            }
                        });
                });
                },
            );
            if let Some(group_id) = toggle_group_collapse_id
                && let Some(group) = self
                    .lineage_node_groups
                    .iter_mut()
                    .find(|group| group.group_id == group_id)
                {
                    group.collapsed = !group.collapsed;
                    persist_workspace_after_frame = true;
                }
        }
        if let Some(selected_node_id) = self.lineage_graph_selected_node_id.as_ref()
            && let Some(selected_row) = graph_rows
                .iter()
                .find(|row| row.node_id == *selected_node_id)
                .cloned()
            {
                if selected_row.kind == LineageNodeKind::Macro {
                    ui.separator();
                    ui.heading("Selected Macro Node");
                    ui.small(format!(
                        "{} ({})",
                        selected_row
                            .macro_instance_id
                            .as_deref()
                            .unwrap_or(&selected_row.node_id),
                        selected_row.display_name
                    ));
                    ui.small(format!(
                        "status={} | template={} | routine={} | ops={}",
                        selected_row.macro_status.as_deref().unwrap_or("ok"),
                        selected_row
                            .macro_template_name
                            .as_deref()
                            .unwrap_or("-"),
                        selected_row
                            .macro_routine_id
                            .as_deref()
                            .unwrap_or("-"),
                        selected_row.macro_op_ids.len()
                    ));
                    if let Some(status_message) = &selected_row.macro_status_message {
                        ui.colored_label(
                            egui::Color32::from_rgb(150, 20, 20),
                            format!("status message: {status_message}"),
                        );
                    }
                    let render_bindings = |ui: &mut Ui,
                                           title: &str,
                                           bindings: &[LineageMacroPortBinding]| {
                        ui.label(title);
                        if bindings.is_empty() {
                            ui.small("- none -");
                            return;
                        }
                        egui::Grid::new(format!(
                            "lineage_macro_bindings_{}_{}",
                            selected_row.node_id, title
                        ))
                        .striped(true)
                        .min_col_width(64.0)
                        .show(ui, |ui| {
                            ui.strong("Port");
                            ui.strong("Kind");
                            ui.strong("Values");
                            ui.end_row();
                            for binding in bindings {
                                ui.monospace(&binding.port_id);
                                ui.label(binding.kind.clone());
                                let values = if binding.values.is_empty() {
                                    "-".to_string()
                                } else {
                                    binding.values.join(", ")
                                };
                                ui.monospace(values);
                                ui.end_row();
                            }
                        });
                    };
                    render_bindings(ui, "Inputs", &selected_row.macro_inputs);
                    render_bindings(ui, "Outputs", &selected_row.macro_outputs);

                    ui.label("Emitted operations");
                    if selected_row.macro_op_ids.is_empty() {
                        ui.small("- none -");
                    } else {
                        egui::Grid::new(format!("lineage_macro_ops_{}", selected_row.node_id))
                            .striped(true)
                            .min_col_width(64.0)
                            .show(ui, |ui| {
                                ui.strong("Op ID");
                                ui.strong("Summary");
                                ui.end_row();
                                for op_id in &selected_row.macro_op_ids {
                                    ui.monospace(op_id);
                                    ui.label(
                                        graph_op_label_by_id
                                            .get(op_id)
                                            .cloned()
                                            .unwrap_or_else(|| "-".to_string()),
                                    );
                                    ui.end_row();
                                }
                            });
                    }

                    let mut output_sequences = selected_row
                        .macro_outputs
                        .iter()
                        .filter(|binding| binding.kind.eq_ignore_ascii_case("sequence"))
                        .flat_map(|binding| binding.values.iter().cloned())
                        .collect::<Vec<_>>();
                    output_sequences.sort();
                    output_sequences.dedup();
                    if !output_sequences.is_empty() {
                        ui.label("Open output sequence");
                        ui.horizontal_wrapped(|ui| {
                            for seq_id in output_sequences {
                                if ui
                                    .small_button(seq_id.clone())
                                    .on_hover_text(
                                        "Open this sequence output in a dedicated sequence window",
                                    )
                                    .clicked()
                                {
                                    open_seq = Some(seq_id);
                                }
                            }
                        });
                    }
                } else if selected_row.kind == LineageNodeKind::Analysis {
                    ui.separator();
                    ui.heading("Selected Analysis Node");
                    let kind = Self::infer_lineage_analysis_kind_from_row(&selected_row)
                        .map(LineageAnalysisKind::as_str)
                        .unwrap_or("analysis");
                    ui.small(format!(
                        "{} ({})",
                        Self::infer_lineage_analysis_artifact_id_from_row(&selected_row)
                            .unwrap_or_else(|| selected_row.node_id.clone()),
                        kind
                    ));
                    ui.small(format!("source sequence={}", selected_row.seq_id));
                    if let Some(reference_seq_id) = selected_row.analysis_reference_seq_id.as_deref()
                    {
                        let reference_label = match Self::infer_lineage_analysis_kind_from_row(
                            &selected_row,
                        ) {
                            Some(LineageAnalysisKind::SequencingConfirmation) => "baseline sequence",
                            Some(LineageAnalysisKind::ReverseTranslation) => "product sequence",
                            _ => "reference sequence",
                        };
                        ui.small(format!("{reference_label}={reference_seq_id}"));
                    }
                    if let Some(mode) = selected_row.analysis_mode.as_deref() {
                        let mode_label = match Self::infer_lineage_analysis_kind_from_row(
                            &selected_row,
                        ) {
                            Some(LineageAnalysisKind::RnaReadInterpretation) => "profile",
                            Some(LineageAnalysisKind::PrimerDesign)
                            | Some(LineageAnalysisKind::QpcrDesign) => "backend",
                            Some(LineageAnalysisKind::RestrictionCloningPcrHandoff) => {
                                "handoff_mode"
                            }
                            Some(LineageAnalysisKind::ProteinDerivation) => "derivation_mode",
                            Some(LineageAnalysisKind::ReverseTranslation) => "speed_profile",
                            Some(LineageAnalysisKind::ConstructReasoning) => "objective_id",
                            Some(LineageAnalysisKind::UniprotProjection) => "entry_id",
                            _ => "mode/model",
                        };
                        ui.small(format!("{mode_label}={mode}"));
                    }
                    if let Some(status) = selected_row.analysis_status.as_deref() {
                        let status_label = match Self::infer_lineage_analysis_kind_from_row(
                            &selected_row,
                        ) {
                            Some(LineageAnalysisKind::RnaReadInterpretation) => {
                                "report_mode/origin_mode"
                            }
                            Some(LineageAnalysisKind::ReverseTranslation) => "diagnostics",
                            Some(LineageAnalysisKind::ConstructReasoning) => "goal",
                            Some(LineageAnalysisKind::UniprotProjection) => "transcript_filter",
                            Some(LineageAnalysisKind::RestrictionCloningPcrHandoff) => {
                                "compatibility_status"
                            }
                            _ => "status",
                        };
                        ui.small(format!("{status_label}={status}"));
                    }
                    if let Some(points) = selected_row.analysis_point_count {
                        let point_label = match Self::infer_lineage_analysis_kind_from_row(
                            &selected_row,
                        ) {
                            Some(LineageAnalysisKind::RnaReadInterpretation) => {
                                "seed_passed_count"
                            }
                            Some(LineageAnalysisKind::ReverseTranslation) => "translation_table",
                            Some(LineageAnalysisKind::ConstructReasoning) => "evidence_count",
                            _ => "point_count",
                        };
                        ui.small(format!("{point_label}={points}"));
                    }
                    if let Some(bin_count) = selected_row.analysis_bin_count {
                        ui.small(format!("bin_count={bin_count}"));
                    }
                    if let Some(read_count) = selected_row.analysis_read_count {
                        ui.small(format!("read_count={read_count}"));
                    }
                    if let Some(trace_count) = selected_row.analysis_trace_count {
                        ui.small(format!("trace_count={trace_count}"));
                    }
                    if let Some(target_count) = selected_row.analysis_target_count {
                        let count_label = match Self::infer_lineage_analysis_kind_from_row(
                            &selected_row,
                        ) {
                            Some(LineageAnalysisKind::RnaReadInterpretation) => {
                                "target_gene_count"
                            }
                            Some(LineageAnalysisKind::PrimerDesign) => "pair_count",
                            Some(LineageAnalysisKind::QpcrDesign) => "assay_count",
                            Some(LineageAnalysisKind::RestrictionCloningPcrHandoff) => {
                                "handoff_count"
                            }
                            Some(LineageAnalysisKind::ProteinDerivation) => "protein_count",
                            Some(LineageAnalysisKind::ReverseTranslation) => "coding_bp",
                            Some(LineageAnalysisKind::ConstructReasoning) => "candidate_count",
                            Some(LineageAnalysisKind::UniprotProjection) => "transcript_count",
                            _ => "target_count",
                        };
                        ui.small(format!("{count_label}={target_count}"));
                    }
                    if let Some(variant_count) = selected_row.analysis_variant_count {
                        let variant_label = match Self::infer_lineage_analysis_kind_from_row(
                            &selected_row,
                        ) {
                            Some(LineageAnalysisKind::RnaReadInterpretation) => {
                                "aligned_read_count"
                            }
                            Some(LineageAnalysisKind::ReverseTranslation) => "protein_aa",
                            Some(LineageAnalysisKind::ConstructReasoning) => "decision_count",
                            _ => "variant_count",
                        };
                        ui.small(format!("{variant_label}={variant_count}"));
                    }
                    ui.horizontal(|ui| {
                        let analysis_payload =
                            Self::lineage_analysis_open_payload(&selected_row);
                        let inferred_dotplot_id = selected_row
                            .analysis_artifact_id
                            .clone()
                            .or_else(|| {
                                Self::infer_lineage_analysis_artifact_id_from_row(&selected_row)
                            });
                        let open_analysis_button_label = match analysis_payload
                            .as_ref()
                            .map(|(kind, _, _)| *kind)
                        {
                            Some(LineageAnalysisKind::Dotplot) => "Open Dotplot View",
                            Some(LineageAnalysisKind::FlexibilityTrack) => "Open Track View",
                            Some(LineageAnalysisKind::RnaReadInterpretation) => {
                                "Open RNA-read Mapping"
                            }
                            Some(LineageAnalysisKind::PrimerDesign) => "Open Primer Report",
                            Some(LineageAnalysisKind::QpcrDesign) => "Open qPCR Report",
                            Some(LineageAnalysisKind::RestrictionCloningPcrHandoff) => {
                                "Open Cloning Handoff"
                            }
                            Some(LineageAnalysisKind::ProteinDerivation) => {
                                "Open Derived Protein Expert"
                            }
                            Some(LineageAnalysisKind::ReverseTranslation) => {
                                "Open Coding Sequence"
                            }
                            Some(LineageAnalysisKind::ConstructReasoning) => {
                                "Open Construct Reasoning"
                            }
                            Some(LineageAnalysisKind::UniprotProjection) => {
                                "Open UniProt Projection"
                            }
                            Some(LineageAnalysisKind::SequencingConfirmation) => {
                                "Open Confirmation View"
                            }
                            None => "Open Analysis View",
                        };
                        if ui
                            .button(open_analysis_button_label)
                            .on_hover_text("Open this analysis artifact view for the sequence")
                            .clicked()
                        {
                            if let Some((kind, seq_id, artifact_id)) = analysis_payload.clone() {
                                open_lineage_analysis = Some((kind, seq_id, artifact_id));
                            } else {
                                open_seq = Some(selected_row.seq_id.clone());
                            }
                        }
                        if ui
                            .button("Open Source Sequence")
                            .on_hover_text("Open the analysis source sequence")
                            .clicked()
                        {
                            open_seq = Some(selected_row.seq_id.clone());
                        }
                        if let Some(reference_seq_id) =
                            selected_row.analysis_reference_seq_id.as_deref()
                        {
                            let reference_button_label =
                                match Self::infer_lineage_analysis_kind_from_row(&selected_row) {
                                    Some(LineageAnalysisKind::SequencingConfirmation) => {
                                        "Open Baseline Sequence"
                                    }
                                    _ => "Open Reference Sequence",
                                };
                            let reference_hover_text =
                                match Self::infer_lineage_analysis_kind_from_row(&selected_row) {
                                    Some(LineageAnalysisKind::SequencingConfirmation) => {
                                        "Open the baseline/reference sequence used to classify intended edits and reversions"
                                    }
                                    _ => "Open the analysis reference sequence",
                                };
                            if ui
                                .button(reference_button_label)
                                .on_hover_text(reference_hover_text)
                                .clicked()
                            {
                                open_seq = Some(reference_seq_id.to_string());
                            }
                        }
                        if matches!(
                            analysis_payload.as_ref().map(|(kind, _, _)| *kind),
                            Some(LineageAnalysisKind::Dotplot)
                        ) {
                            let dotplot_id = inferred_dotplot_id
                                .as_deref()
                                .map(str::trim)
                                .filter(|value| !value.is_empty())
                                .map(|value| value.to_string());
                            if let Some(dotplot_id) = dotplot_id
                                && ui
                                    .button("Dotplot SVG")
                                    .on_hover_text(
                                        "Export this dotplot analysis as an SVG artifact",
                                    )
                                    .clicked()
                                {
                                    export_lineage_dotplot_svg = Some((
                                        selected_row.seq_id.clone(),
                                        dotplot_id,
                                        selected_row.analysis_reference_seq_id.clone(),
                                    ));
                                }
                        }
                    });
                }
            }
        ui.separator();
        let splitter_width = ui.available_width().max(1.0);
        let splitter_response = ui
            .add_sized(
                egui::vec2(splitter_width, 14.0),
                egui::Label::new(egui::RichText::new("⋮⋮").monospace())
                    .sense(egui::Sense::click_and_drag()),
            )
            .on_hover_text("Drag to resize graph/table and container sections");
        let splitter_rect = splitter_response.rect;
        let splitter_active = splitter_response.hovered() || splitter_response.dragged();
        let splitter_color = if splitter_active {
            egui::Color32::from_rgb(120, 120, 120)
        } else {
            egui::Color32::from_rgb(150, 150, 150)
        };
        ui.painter().line_segment(
            [
                egui::pos2(splitter_rect.left(), splitter_rect.center().y),
                egui::pos2(splitter_rect.right(), splitter_rect.center().y),
            ],
            egui::Stroke::new(1.2, splitter_color),
        );
        ui.painter().text(
            splitter_rect.center(),
            egui::Align2::CENTER_CENTER,
            "⋮⋮",
            egui::FontId::monospace(10.0),
            splitter_color,
        );
        if splitter_active {
            ui.output_mut(|o| o.cursor_icon = egui::CursorIcon::ResizeVertical);
        }
        if splitter_response.drag_started() {
            self.lineage_main_split_drag_origin =
                Some(self.lineage_main_split_fraction.clamp(0.2, 0.9));
            self.lineage_main_split_drag_start_y =
                ui.input(|i| i.pointer.interact_pos().map(|pos| pos.y));
        }
        let pointer_down = ui.input(|i| i.pointer.primary_down());
        let pointer_y = ui.input(|i| i.pointer.interact_pos().map(|pos| pos.y));
        splitter_dragging =
            splitter_response.dragged() || (pointer_down && self.lineage_main_split_drag_origin.is_some());
        let splitter_drag_finished = splitter_response.drag_stopped();
        let splitter_drag_active = splitter_dragging || splitter_drag_finished;
        if splitter_drag_active
            && let (Some(drag_origin_split), Some(start_y), Some(current_y)) = (
                self.lineage_main_split_drag_origin,
                self.lineage_main_split_drag_start_y,
                pointer_y,
            ) {
                let drag_delta_y = current_y - start_y;
                let total_height = (graph_area_height + container_area_height).max(400.0);
                let drag_base_graph_height =
                    (total_height * drag_origin_split.clamp(0.2, 0.9))
                        .clamp(LINEAGE_MAIN_TOP_PANEL_MIN_HEIGHT, total_height - 80.0);
                let next_graph_height =
                    (drag_base_graph_height + drag_delta_y)
                        .clamp(LINEAGE_MAIN_TOP_PANEL_MIN_HEIGHT, total_height - 80.0);
                if (next_graph_height - graph_area_height).abs() > f32::EPSILON {
                    graph_area_height = next_graph_height;
                    container_area_height =
                        (total_height - graph_area_height)
                            .clamp(80.0, total_height - LINEAGE_MAIN_TOP_PANEL_MIN_HEIGHT);
                    main_split_fraction = (graph_area_height / total_height).clamp(0.2, 0.9);
                    self.lineage_main_split_fraction = main_split_fraction;
                }
            }
        if splitter_drag_finished || (!pointer_down && self.lineage_main_split_drag_origin.is_some()) {
            self.lineage_main_split_drag_origin = None;
            self.lineage_main_split_drag_start_y = None;
            persist_workspace_after_frame = true;
        }
        ui.small("Drag split bar to resize graph/table and container sections");
        let container_panel_width = ui.available_width().max(1.0);
        let subpanel_total_height = container_area_height.max(220.0);
        let (mut containers_panel_height, mut arrangements_panel_height, normalized_subsplit) =
            Self::normalized_lineage_container_subpanel_heights(
                subpanel_total_height,
                container_arrangement_split_fraction,
            );
        container_arrangement_split_fraction = normalized_subsplit;

        ui.allocate_ui_with_layout(
            egui::vec2(container_panel_width, containers_panel_height),
            egui::Layout::top_down(egui::Align::Min),
            |ui| {
                ui.set_min_width(container_panel_width);
                ui.set_max_width(container_panel_width);
                ui.heading("Containers");
                ui.label("Container-level view of candidate sequence sets");
                egui::ScrollArea::both()
                    .id_salt("lineage_container_grid_scroll")
                    .auto_shrink([false, false])
                    .max_height(ui.available_height().max(40.0))
                    .show(ui, |ui| {
                        scroll_input_policy::apply_scrollarea_keyboard_navigation(
                            ui,
                            scroll_input_policy::DEFAULT_SCROLLAREA_KEYBOARD_STEP,
                        );
                        egui::Grid::new("container_grid")
                            .striped(true)
                            .show(ui, |ui| {
                                ui.strong("Container");
                                ui.strong("Kind");
                                ui.strong("Contents");
                                ui.strong("Members");
                                ui.strong("Representative");
                                ui.strong("Action");
                                ui.end_row();
                                for c in &self.lineage_containers {
                                    ui.monospace(&c.container_id);
                                    ui.label(&c.kind);
                                    ui.vertical(|ui| {
                                        ui.label(Self::container_contents_mode_label(
                                            c.declared_contents_exclusive,
                                        ))
                                        .on_hover_text(
                                            Self::container_contents_mode_hover(
                                                c.declared_contents_exclusive,
                                            ),
                                        );
                                        let toggle_label = if c.declared_contents_exclusive {
                                            "Mark subset"
                                        } else {
                                            "Mark exhaustive"
                                        };
                                        if ui
                                            .small_button(toggle_label)
                                            .on_hover_text(if c.declared_contents_exclusive {
                                                "Switch this container to a non-exclusive measured-subset interpretation."
                                            } else {
                                                "Mark this container as declared exhaustive again."
                                            })
                                            .clicked()
                                        {
                                            set_container_exclusivity = Some((
                                                c.container_id.clone(),
                                                !c.declared_contents_exclusive,
                                            ));
                                        }
                                    });
                                    ui.monospace(format!("{}", c.member_count));
                                    ui.monospace(&c.representative);
                                    ui.horizontal(|ui| {
                                        if c.member_count > 1 {
                                            if ui
                                                .button("Open Pool")
                                                .on_hover_text("Open this container as a pool view")
                                                .clicked()
                                            {
                                                open_pool = Some((
                                                    c.representative.clone(),
                                                    c.members.clone(),
                                                ));
                                            }
                                        } else if !c.representative.is_empty() {
                                            if ui
                                                .button("Open Seq")
                                                .on_hover_text(
                                                    "Open this representative sequence",
                                                )
                                                .clicked()
                                            {
                                                open_seq = Some(c.representative.clone());
                                            }
                                        } else {
                                            ui.label("-");
                                        }
                                        if c.member_count > 0
                                            && ui
                                                .button("Gel SVG")
                                                .on_hover_text(
                                                    "Export a serial gel lane for this container",
                                                )
                                                .clicked()
                                        {
                                            export_container_gel = Some((
                                                c.container_id.clone(),
                                                vec![c.container_id.clone()],
                                            ));
                                        }
                                    });
                                    ui.end_row();
                                }
                            });
                    });
            },
        );

        let sub_split_response = ui
            .add_sized(
                egui::vec2(container_panel_width, 14.0),
                egui::Label::new(egui::RichText::new("⋮⋮").monospace())
                    .sense(egui::Sense::click_and_drag()),
            )
            .on_hover_text("Drag to resize containers and arrangements");
        let sub_split_rect = sub_split_response.rect;
        let sub_split_active = sub_split_response.hovered() || sub_split_response.dragged();
        let sub_split_color = if sub_split_active {
            egui::Color32::from_rgb(120, 120, 120)
        } else {
            egui::Color32::from_rgb(150, 150, 150)
        };
        ui.painter().line_segment(
            [
                egui::pos2(sub_split_rect.left(), sub_split_rect.center().y),
                egui::pos2(sub_split_rect.right(), sub_split_rect.center().y),
            ],
            egui::Stroke::new(1.2, sub_split_color),
        );
        ui.painter().text(
            sub_split_rect.center(),
            egui::Align2::CENTER_CENTER,
            "⋮⋮",
            egui::FontId::monospace(10.0),
            sub_split_color,
        );
        if sub_split_active {
            ui.output_mut(|o| o.cursor_icon = egui::CursorIcon::ResizeVertical);
        }
        if sub_split_response.drag_started() {
            self.lineage_container_arrangement_split_drag_origin =
                Some(container_arrangement_split_fraction.clamp(0.2, 0.8));
            self.lineage_container_arrangement_split_drag_start_y =
                ui.input(|i| i.pointer.interact_pos().map(|pos| pos.y));
        }
        let sub_pointer_down = ui.input(|i| i.pointer.primary_down());
        let sub_pointer_y = ui.input(|i| i.pointer.interact_pos().map(|pos| pos.y));
        let sub_split_dragging = sub_split_response.dragged()
            || (sub_pointer_down && self.lineage_container_arrangement_split_drag_origin.is_some());
        let sub_split_finished = sub_split_response.drag_stopped();
        let sub_split_active_drag = sub_split_dragging || sub_split_finished;
        if sub_split_active_drag
            && let (Some(drag_origin_split), Some(start_y), Some(current_y)) = (
                self.lineage_container_arrangement_split_drag_origin,
                self.lineage_container_arrangement_split_drag_start_y,
                sub_pointer_y,
            ) {
                let drag_delta_y = current_y - start_y;
                let total_height = (containers_panel_height + arrangements_panel_height).max(160.0);
                let drag_base_containers_height =
                    (total_height * drag_origin_split.clamp(0.2, 0.8)).clamp(80.0, total_height - 80.0);
                let next_containers_height =
                    (drag_base_containers_height + drag_delta_y).clamp(80.0, total_height - 80.0);
                if (next_containers_height - containers_panel_height).abs() > f32::EPSILON {
                    containers_panel_height = next_containers_height;
                    arrangements_panel_height =
                        (total_height - containers_panel_height).clamp(80.0, total_height - 80.0);
                    container_arrangement_split_fraction =
                        (containers_panel_height / total_height).clamp(0.2, 0.8);
                }
            }
        if sub_split_finished
            || (!sub_pointer_down
                && self.lineage_container_arrangement_split_drag_origin.is_some())
        {
            self.lineage_container_arrangement_split_drag_origin = None;
            self.lineage_container_arrangement_split_drag_start_y = None;
            persist_workspace_after_frame = true;
        }
        ui.small("Drag split bar to resize containers and arrangements");

        ui.allocate_ui_with_layout(
            egui::vec2(container_panel_width, arrangements_panel_height),
            egui::Layout::top_down(egui::Align::Min),
            |ui| {
                ui.set_min_width(container_panel_width);
                ui.set_max_width(container_panel_width);
                ui.heading("Arrangements");
                ui.label("Serial lane setups from one or more containers");
                ui.horizontal(|ui| {
                    ui.label("Label preset");
                    egui::ComboBox::from_id_salt("arrangement_label_sheet_preset_combo")
                        .selected_text(Self::rack_label_sheet_preset_label(
                            self.rack_label_sheet_preset,
                        ))
                        .show_ui(ui, |ui| {
                            for preset in [
                                RackLabelSheetPreset::CompactCards,
                                RackLabelSheetPreset::PrintA4,
                                RackLabelSheetPreset::WideCards,
                            ] {
                                ui.selectable_value(
                                    &mut self.rack_label_sheet_preset,
                                    preset,
                                    Self::rack_label_sheet_preset_label(preset),
                                );
                            }
                        });
                    ui.small(
                        "Used by arrangement-scoped and rack-wide label SVG export.",
                    );
                });
                ui.horizontal(|ui| {
                    ui.label("Carrier preset");
                    egui::ComboBox::from_id_salt("arrangement_carrier_label_preset_combo")
                        .selected_text(Self::rack_carrier_label_preset_label(
                            self.rack_carrier_label_preset,
                        ))
                        .show_ui(ui, |ui| {
                            for preset in [
                                RackCarrierLabelPreset::FrontStripAndCards,
                                RackCarrierLabelPreset::FrontStripOnly,
                                RackCarrierLabelPreset::ModuleCardsOnly,
                            ] {
                                ui.selectable_value(
                                    &mut self.rack_carrier_label_preset,
                                    preset,
                                    Self::rack_carrier_label_preset_label(preset),
                                );
                            }
                        });
                    ui.small(
                        "Used by arrangement-scoped carrier/front-strip SVG export.",
                    );
                });
                egui::ScrollArea::both()
                    .id_salt("lineage_arrangement_grid_scroll")
                    .auto_shrink([false, false])
                    .max_height(ui.available_height().max(40.0))
                    .show(ui, |ui| {
                        scroll_input_policy::apply_scrollarea_keyboard_navigation(
                            ui,
                            scroll_input_policy::DEFAULT_SCROLLAREA_KEYBOARD_STEP,
                        );
                        if self.lineage_arrangements.is_empty() {
                            ui.label("No arrangements recorded");
                        } else {
                            let arrangement_rows = self.lineage_arrangements.clone();
                            egui::Grid::new("arrangement_grid")
                                .striped(true)
                                .show(ui, |ui| {
                                    ui.strong("Arrangement");
                                    ui.strong("Mode");
                                    ui.strong("Name");
                                    ui.strong("Lanes");
                                    ui.strong("Lane containers");
                                    ui.strong("Ladders");
                                    ui.strong("Rack");
                                    ui.strong("Action");
                                    ui.end_row();
                                    for arrangement in &arrangement_rows {
                                        ui.monospace(&arrangement.arrangement_id);
                                        ui.label(&arrangement.mode);
                                        ui.label(if arrangement.name.trim().is_empty() {
                                            "-".to_string()
                                        } else {
                                            arrangement.name.clone()
                                        });
                                        ui.monospace(arrangement.lane_count.to_string());
                                        ui.label(if arrangement.lane_container_ids.is_empty() {
                                            "-".to_string()
                                        } else {
                                            arrangement.lane_container_ids.join(", ")
                                        });
                                        ui.label(if arrangement.ladders.is_empty() {
                                            "auto".to_string()
                                        } else {
                                            arrangement.ladders.join(", ")
                                        });
                                        ui.label(
                                            arrangement
                                                .default_rack_id
                                                .clone()
                                                .unwrap_or_else(|| "draft on demand".to_string()),
                                        );
                                        ui.horizontal(|ui| {
                                            if ui
                                                .button("Preview Gel")
                                                .on_hover_text(
                                                    "Open an in-app gel preview where left/right ladder choices update the visual result together",
                                                )
                                                .clicked()
                                            {
                                                open_arrangement_gel_preview =
                                                    Some(arrangement.arrangement_id.clone());
                                            }
                                            if ui
                                                .button("Export Gel")
                                                .on_hover_text(
                                                    "Export one serial gel using this arrangement",
                                                )
                                                .clicked()
                                            {
                                                let stem = if arrangement.name.trim().is_empty() {
                                                    arrangement.arrangement_id.clone()
                                                } else {
                                                    arrangement.name.clone()
                                                };
                                                export_arrangement_gel = Some((
                                                    stem,
                                                    arrangement.arrangement_id.clone(),
                                                ));
                                            }
                                            if !arrangement.lane_container_ids.is_empty() {
                                                self.render_open_lanes_menu(
                                                    ui,
                                                    &arrangement.lane_container_ids,
                                                    &mut open_lane_containers,
                                                );
                                            }
                                            if ui
                                                .button("Open Rack")
                                                .on_hover_text(
                                                    "Open the linked physical rack draft for this arrangement, creating the default one when needed",
                                                )
                                                .clicked()
                                            {
                                                self.open_arrangement_rack_dialog(
                                                    &arrangement.arrangement_id,
                                                );
                                            }
                                            if ui
                                                .button("Preview Labels")
                                                .on_hover_text(
                                                    "Open an in-app preview of the deterministic label sheet for this arrangement from its linked rack draft before exporting it",
                                                )
                                                .clicked()
                                            {
                                                self.open_arrangement_labels_preview_dialog(
                                                    &arrangement.arrangement_id,
                                                );
                                            }
                                            if ui
                                                .button("Carrier SVG")
                                                .on_hover_text(
                                                    "Export carrier-matched front-strip and module-label SVGs for this arrangement from its linked rack draft using the current physical template",
                                                )
                                                .clicked()
                                            {
                                                self.prompt_export_arrangement_carrier_labels_svg(
                                                    &arrangement.arrangement_id,
                                                );
                                            }
                                            if !self.lineage_racks.is_empty()
                                                && ui
                                                    .button("Place on Existing Rack...")
                                                    .on_hover_text(
                                                        "Append this arrangement as one contiguous block onto another saved rack",
                                                    )
                                                    .clicked()
                                            {
                                                self.open_place_arrangement_on_rack_dialog(
                                                    &arrangement.arrangement_id,
                                                );
                                            }
                                        });
                                        ui.end_row();
                                    }
                                });
                        }
                    });
            },
        );
            });

        if let Some((stem, container_ids)) = export_container_gel.take() {
            self.prompt_export_serial_gel_svg(&stem, Some(container_ids), None, None, None);
        }
        if let Some((container_id, exclusive)) = set_container_exclusivity.take() {
            let update_result = {
                self.engine.write().unwrap().apply(
                    Operation::SetContainerDeclaredContentsExclusive {
                        container_id: container_id.clone(),
                        exclusive,
                    },
                )
            };
            match update_result {
                Ok(op_result) => {
                    self.app_status = op_result.messages.first().cloned().unwrap_or_else(|| {
                        format!(
                            "Updated container '{}' contents mode to {}",
                            container_id,
                            if exclusive {
                                "declared-only"
                            } else {
                                "known-subset"
                            }
                        )
                    });
                    self.lineage_cache_valid = false;
                    self.refresh_lineage_cache_if_needed();
                }
                Err(err) => {
                    self.app_status = format!(
                        "Could not update container '{}' contents mode: {}",
                        container_id, err.message
                    );
                }
            }
        }
        if let Some((stem, arrangement_id)) = export_arrangement_gel.take() {
            self.prompt_export_serial_gel_svg(&stem, None, Some(arrangement_id), None, None);
        }
        if let Some(arrangement_id) = open_arrangement_gel_preview.take() {
            self.open_arrangement_gel_preview_dialog(&arrangement_id);
        }
        if let Some((seq_id, dotplot_id, reference_seq_id)) = export_lineage_dotplot_svg.take() {
            self.prompt_export_dotplot_svg_from_lineage(
                &seq_id,
                &dotplot_id,
                reference_seq_id.as_deref(),
            );
        }
        if request_export_lineage_svg {
            self.prompt_export_lineage_svg();
        }
        if let Some(container_ids) = open_lane_containers.take() {
            for container_id in &container_ids {
                if let Some(container_row) = self
                    .lineage_containers
                    .iter()
                    .find(|row| row.container_id == *container_id)
                    .cloned()
                {
                    if container_row.member_count > 1 {
                        self.open_pool_window_compact(
                            &container_row.representative,
                            container_row.members.clone(),
                        );
                    } else if !container_row.representative.is_empty() {
                        self.open_sequence_window_compact(&container_row.representative);
                    }
                }
            }
        }

        if (self.lineage_graph_zoom - graph_zoom).abs() > 0.0001 {
            self.lineage_graph_zoom = graph_zoom;
            persist_workspace_after_frame = true;
        }
        if (self.lineage_graph_area_height - graph_area_height).abs() > 0.5 {
            self.lineage_graph_area_height = graph_area_height;
            persist_workspace_after_frame = true;
        }
        if (self.lineage_container_area_height - container_area_height).abs() > 0.5 {
            self.lineage_container_area_height = container_area_height;
            persist_workspace_after_frame = true;
        }
        if (self.lineage_container_arrangement_split_fraction
            - container_arrangement_split_fraction)
            .abs()
            > 0.0001
        {
            self.lineage_container_arrangement_split_fraction =
                container_arrangement_split_fraction;
            persist_workspace_after_frame = true;
        }
        let combined_panel_height =
            (graph_area_height.max(0.0) + container_area_height.max(0.0)).max(1.0);
        let computed_main_split_fraction =
            (graph_area_height.max(0.0) / combined_panel_height).clamp(0.2, 0.9);
        main_split_fraction = computed_main_split_fraction;
        if !splitter_dragging
            && (self.lineage_main_split_fraction - main_split_fraction).abs() > 0.0001
        {
            self.lineage_main_split_fraction = main_split_fraction;
            persist_workspace_after_frame = true;
        }
        if self.lineage_graph_compact_labels != graph_compact_labels {
            self.lineage_graph_compact_labels = graph_compact_labels;
            persist_workspace_after_frame = true;
        }
        if (self.lineage_graph_scroll_offset - graph_scroll_offset).length_sq() > 0.25 {
            self.lineage_graph_scroll_offset = graph_scroll_offset;
            persist_workspace_after_frame = true;
        }
        self.render_lineage_node_rename_window(ui.ctx(), &mut persist_workspace_after_frame);
        self.render_lineage_node_remove_confirm_window(
            ui.ctx(),
            &mut persist_workspace_after_frame,
        );
        if persist_workspace_after_frame {
            self.persist_lineage_graph_workspace_to_state();
        }

        if let Some(input) = select_candidate_from {
            let criterion = format!("gui_lineage_select:{input}");
            let result = self
                .engine
                .write()
                .unwrap()
                .apply(Operation::SelectCandidate {
                    input: input.clone(),
                    criterion,
                    output_id: None,
                });
            if let Ok(op_result) = result
                && let Some(seq_id) = op_result.created_seq_ids.first()
            {
                open_seq = Some(seq_id.clone());
            }
        }

        if let Some(descriptor) = open_lineage_retrieval {
            self.open_lineage_retrieval_descriptor(&descriptor);
        }
        if let Some((kind, seq_id, artifact_id)) = open_lineage_analysis {
            self.open_lineage_analysis_artifact(kind, &seq_id, &artifact_id);
        }
        if let Some(op_id) = reopen_gibson_from_operation {
            match self.reopen_gibson_specialist_from_operation(&op_id) {
                Ok(true) => {}
                Ok(false) => {
                    self.app_status = format!(
                        "Operation '{}' does not carry a reopenable Gibson plan payload",
                        op_id
                    );
                }
                Err(err) => {
                    self.app_status = format!(
                        "Could not reopen Gibson specialist from operation '{}': {}",
                        op_id, err
                    );
                }
            }
        }
        if let Some(op_id) = reopen_pcr_from_operation {
            match self.reopen_pcr_designer_from_operation(&op_id) {
                Ok(true) => {}
                Ok(false) => {
                    self.app_status = format!(
                        "Operation '{}' does not carry a reopenable PCR template context",
                        op_id
                    );
                }
                Err(err) => {
                    self.app_status = format!(
                        "Could not open PCR Designer from operation '{}': {}",
                        op_id, err
                    );
                }
            }
        }
        if let Some(seq_id) = open_seq {
            self.open_sequence_window(&seq_id);
        }
        if let Some((representative, pool_members)) = open_pool {
            self.open_pool_window(&representative, pool_members);
        }
    }
}

#[cfg(test)]
mod tests {
    use super::lineage_graph_canvas_width;

    #[test]
    fn lineage_graph_canvas_width_uses_parent_viewport_snapshot() {
        let base_width = 440.0;
        let zoom = 1.0;
        let parent_viewport_width = 720.0;
        let inflated_child_available_width = 1280.0;

        assert_eq!(
            lineage_graph_canvas_width(base_width, zoom, parent_viewport_width),
            parent_viewport_width
        );
        assert_ne!(
            lineage_graph_canvas_width(base_width, zoom, parent_viewport_width),
            inflated_child_available_width
        );
    }
}
